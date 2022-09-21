//
// 
// This code is generated from CPPINTS, a C++ program to generate the 
// analytical integrals based on Gaussian form primitive functions. 
// Copyright (C) 2015 The State University of New York at Buffalo 
// This software uses the MIT license as below: 
// 
// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, 
// and/or sell copies of the Software, and to permit persons to whom the Software 
// is furnished to do so, subject to the following conditions: 
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software. 
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  
// DEALINGS IN THE SOFTWARE. 
// 
// 

#include "constants.h"
#include <cstddef>
#include <math.h>

typedef int             Int;
typedef size_t          UInt;
#define THRESHOLD_MATH  0.00000000000001
#ifdef WITH_SINGLE_PRECISION
typedef float           Double;
#else
typedef double          Double;
#endif

//
//  here below is a list of variables used in the program
//
//  alpha is the bra1's exponent
//  beta  is the bra2's exponent
//  gamma is the ket1's exponent
//  delta is the ket2's exponent
//  A is the nuclear center for bra1
//  B is the nuclear center for bra2
//  C is the nuclear center for ket1
//  D is the nuclear center for ket2
//  P is the new center after bra1 combined with bra2
//  Q is the new center after ket1 combined with ket2
//  W is the new center after P combined with Q
//  pMax is maximum value of corresponding density matrix block(or value pair), used for ERI
//  omega is the exponential factor used for operator in form of erf(omega*r12)/r12
//  also omega could be the exponential factor used for operator in form of e^(-omega*r12^2)
//
//  variables:
//
//  zeta      = alpha + beta
//  eta       = gamma + delta
//  oned2z    = 1/(2*zeta)
//  oned2e    = 1/(2*eta)
//  onedz     = 1/zeta
//  onede     = 1/eta
//  kappa     = zeta + eta
//  onedk     = 1/kappa
//  oned2zeta = 1/(2*(alpha+beta+gamma))
//  xi        = alpha*beta*onedz
//  twoxi     = 2*alpha*beta*onedz
//  rho       = zeta*eta*onedk
//  rhod2zsq  = rho/(2*zeta*zeta)
//  rhod2esq  = rho/(2*eta*eta)
//  odorho    = omega/(rho+omega)
//  rhodorho  = rho/(rho+omega)
//  orhod2z2  = (rho/(2*zeta*zeta))*(omega/(rho+omega))
//  orhod2e2  = (rho/(2*eta*eta))*(omega/(rho+omega))
//  od2k      = (1/(2*kappa))*(omega/(rho+omega))
//  adz       = alpha*onedz
//  bdz       = beta*onedz
//  gde       = gamma*onede
//  gde       = delta*onede
//
//  input parameters based on primitive functions pair:
//
//  bra side shell pair is index as i
//  inp2  is the number of primitive pairs
//  iexp  is the array of 1/(alpha+beta)
//  icoe  is the array of ic_bra1*ic_bra2
//  ifac  is the array of pre-factor on bra side 
//  for (SS|SS)^{m} etc. type of integrals
//  ket side shell pair is index as j
//  jnp2  is the number of primitive pairs
//  jexp  is the array of 1/(gamma+delta)
//  jcoe  is the array of jc_ket1*jc_ket2
//  jfac  is the array of pre-factor on ket side 
//  for (SS|SS)^{m} etc. type of integrals
//

void hgp_os_eri_g_s_g_s(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
    Double ic2   = icoe[ip2];
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt jp2=0; jp2<jnp2; jp2++) {
      Double onede = jexp[jp2];
      Double jc2   = jcoe[jp2];
      Double fket  = jfac[jp2];
      Double pref      = fbra*fket;
      Double prefactor = ic2*jc2*pref;

      // 
      // here below the code is performing significance test for integrals on
      // primitive integrals. Here we use the overlap integrals to roughly 
      // estimate the order of the result integrals
      // the threshold value should be for primitive function quartet, we compare
      // the value against machine precision for significance test
      // 
      Double I_ERI_S_S_S_S_vrr_IntegralTest = pref;
      if (fabs(ic2*jc2)>1.0E0) {
        I_ERI_S_S_S_S_vrr_IntegralTest = prefactor;
      }

      // test the integrals with the pMax, which is the maximum value
      // of the corresponding density matrix block(or it may be maximum
      // value pair of the corresponding density matrix block)
      if(fabs(I_ERI_S_S_S_S_vrr_IntegralTest*pMax)<THRESHOLD_MATH) continue;


      UInt offsetQ  = 3*jp2;
      Double QX    = Q[offsetQ  ];
      Double QY    = Q[offsetQ+1];
      Double QZ    = Q[offsetQ+2];
      Double rho   = 1.0E0/(onedz+onede);
      Double sqrho = sqrt(rho);
      Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);
      Double u     = rho*PQ2;
      if (withErfR12) u = PQ2/(1.0E0/(omega*omega)+1.0E0/rho);
      Double squ   = sqrt(u);
      Double QCX   = QX - C[0];
      Double QCY   = QY - C[1];
      Double QCZ   = QZ - C[2];
      Double WX    = rho*(PX*onede + QX*onedz);
      Double WY    = rho*(PY*onede + QY*onedz);
      Double WZ    = rho*(PZ*onede + QZ*onedz);
      Double oned2k= 0.5E0*rho*onede*onedz;
      Double WPX   = WX - PX;
      Double WPY   = WY - PY;
      Double WPZ   = WZ - PZ;
      Double rhod2zsq = rho*oned2z*onedz;
      Double WQX   = WX - QX;
      Double WQY   = WY - QY;
      Double WQZ   = WZ - QZ;
      Double oned2e= 0.5E0*onede;
      Double rhod2esq= rho*oned2e*onede;


      //
      //
      // now here for maxM>0 to compute the infamous incomplete Gamma function f_{m}(u)
      // the implementation is divided in two situations:
      // 1  if u <=1.8; use power series to get f_{Mmax}(u), then use down recursive
      //    relation to get the rest of incomplete Gamma functions;
      // 2  for u >1.8 and M <= 10 we calculate erf(u), then use up recursive
      //    relation to calculate the rest of results
      // 3  for u> 1.8 and M >  10 we calculate f_{Mmax}(u) then use down 
      //    recursive relation to get rest of incomplete Gamma functions 
      // The above procedure is tested for u between 0 to 40 with step length 1.0E-6
      // (or 1.0E-5 for float double data), for up recursive relation it shows the error
      // within 1.0E-12 (for M_limit = 12 or error within 1.0E-6 for float type of data
      // For the polynomial expansion and down recursive procedure the error is within 
      // 1.0E-14. All of the testing details please refer to the fmt_test folder
      // 
      // There's one thing need to note for up recursive process. We found that the up
      // recursive procedure is only stable for maxM<=10 and u>1.8 with double
      // precision data, single precision data will lose accuracy quickly so the result
      // for single precision calculation is not doable. Therefore if the "WITH_SINGLE_PRECISION"
      // is defined, then for erf function calculation as well as up recursive
      // process we will use the double type of data
      // 
      //

      Double I_ERI_S_S_S_S_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M1_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M2_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M3_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M4_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M5_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M6_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M7_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M8_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M5_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M6_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M7_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M8_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER51;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER49*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = ONEOVER17*I_ERI_S_S_S_S_M8_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M8_vrr  = f*I_ERI_S_S_S_S_M8_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M7_vrr  = ONEOVER15*(u2*I_ERI_S_S_S_S_M8_vrr+f);
        I_ERI_S_S_S_S_M6_vrr  = ONEOVER13*(u2*I_ERI_S_S_S_S_M7_vrr+f);
        I_ERI_S_S_S_S_M5_vrr  = ONEOVER11*(u2*I_ERI_S_S_S_S_M6_vrr+f);
        I_ERI_S_S_S_S_M4_vrr  = ONEOVER9*(u2*I_ERI_S_S_S_S_M5_vrr+f);
        I_ERI_S_S_S_S_M3_vrr  = ONEOVER7*(u2*I_ERI_S_S_S_S_M4_vrr+f);
        I_ERI_S_S_S_S_M2_vrr  = ONEOVER5*(u2*I_ERI_S_S_S_S_M3_vrr+f);
        I_ERI_S_S_S_S_M1_vrr  = ONEOVER3*(u2*I_ERI_S_S_S_S_M2_vrr+f);
        I_ERI_S_S_S_S_vrr  = ONEOVER1*(u2*I_ERI_S_S_S_S_M1_vrr+f);

      }else{
#ifdef WITH_SINGLE_PRECISION

        // recompute the variable in terms of double accuracy
        double u_d     = u;
        double rho_d   = rho;
        double fac_d   = prefactor;
        double sqrho_d = sqrt(rho_d);
        double squ_d   = sqrt(u_d);

        // use erf function to get (SS|SS)^{0}
        if (fabs(u_d)<THRESHOLD_MATH) {
          I_ERI_S_S_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ERI_S_S_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ERI_S_S_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ERI_S_S_S_S_vrr_d-f);
        I_ERI_S_S_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ERI_S_S_S_S_M1_vrr_d-f);
        I_ERI_S_S_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ERI_S_S_S_S_M2_vrr_d-f);
        I_ERI_S_S_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ERI_S_S_S_S_M3_vrr_d-f);
        I_ERI_S_S_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ERI_S_S_S_S_M4_vrr_d-f);
        I_ERI_S_S_S_S_M6_vrr_d = oneO2u*(11.0E0*I_ERI_S_S_S_S_M5_vrr_d-f);
        I_ERI_S_S_S_S_M7_vrr_d = oneO2u*(13.0E0*I_ERI_S_S_S_S_M6_vrr_d-f);
        I_ERI_S_S_S_S_M8_vrr_d = oneO2u*(15.0E0*I_ERI_S_S_S_S_M7_vrr_d-f);

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);
        I_ERI_S_S_S_S_M5_vrr = static_cast<Double>(I_ERI_S_S_S_S_M5_vrr_d);
        I_ERI_S_S_S_S_M6_vrr = static_cast<Double>(I_ERI_S_S_S_S_M6_vrr_d);
        I_ERI_S_S_S_S_M7_vrr = static_cast<Double>(I_ERI_S_S_S_S_M7_vrr_d);
        I_ERI_S_S_S_S_M8_vrr = static_cast<Double>(I_ERI_S_S_S_S_M8_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ERI_S_S_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ERI_S_S_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M1_vrr = oneO2u*(1.0E0*I_ERI_S_S_S_S_vrr-f);
        I_ERI_S_S_S_S_M2_vrr = oneO2u*(3.0E0*I_ERI_S_S_S_S_M1_vrr-f);
        I_ERI_S_S_S_S_M3_vrr = oneO2u*(5.0E0*I_ERI_S_S_S_S_M2_vrr-f);
        I_ERI_S_S_S_S_M4_vrr = oneO2u*(7.0E0*I_ERI_S_S_S_S_M3_vrr-f);
        I_ERI_S_S_S_S_M5_vrr = oneO2u*(9.0E0*I_ERI_S_S_S_S_M4_vrr-f);
        I_ERI_S_S_S_S_M6_vrr = oneO2u*(11.0E0*I_ERI_S_S_S_S_M5_vrr-f);
        I_ERI_S_S_S_S_M7_vrr = oneO2u*(13.0E0*I_ERI_S_S_S_S_M6_vrr-f);
        I_ERI_S_S_S_S_M8_vrr = oneO2u*(15.0E0*I_ERI_S_S_S_S_M7_vrr-f);

#endif

      }


      // now scale the bottom integral if oper in erf(r12)/r12 form
      if (withErfR12) {
        Double erfPref0   = 1.0E0+rho/(omega*omega);
        Double erfPref1   = 1.0E0/erfPref0;
        Double erfp       = sqrt(erfPref1);
        Double erfp2      = erfp*erfp;
        Double erfPref_1  = erfp;
        I_ERI_S_S_S_S_vrr = I_ERI_S_S_S_S_vrr*erfPref_1;
        Double erfPref_3 = erfPref_1*erfp2;
        Double erfPref_5 = erfPref_3*erfp2;
        Double erfPref_7 = erfPref_5*erfp2;
        Double erfPref_9 = erfPref_7*erfp2;
        Double erfPref_11 = erfPref_9*erfp2;
        Double erfPref_13 = erfPref_11*erfp2;
        Double erfPref_15 = erfPref_13*erfp2;
        Double erfPref_17 = erfPref_15*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
        I_ERI_S_S_S_S_M8_vrr = I_ERI_S_S_S_S_M8_vrr*erfPref_17;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M7
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_S_S_Px_S_M7_vrr = QCX*I_ERI_S_S_S_S_M7_vrr+WQX*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Py_S_M7_vrr = QCY*I_ERI_S_S_S_S_M7_vrr+WQY*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Pz_S_M7_vrr = QCZ*I_ERI_S_S_S_S_M7_vrr+WQZ*I_ERI_S_S_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_Px_S_M6_vrr = QCX*I_ERI_S_S_S_S_M6_vrr+WQX*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Py_S_M6_vrr = QCY*I_ERI_S_S_S_S_M6_vrr+WQY*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Pz_S_M6_vrr = QCZ*I_ERI_S_S_S_S_M6_vrr+WQZ*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M6_vrr = QCX*I_ERI_S_S_Px_S_M6_vrr+WQX*I_ERI_S_S_Px_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_D2y_S_M6_vrr = QCY*I_ERI_S_S_Py_S_M6_vrr+WQY*I_ERI_S_S_Py_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_D2z_S_M6_vrr = QCZ*I_ERI_S_S_Pz_S_M6_vrr+WQZ*I_ERI_S_S_Pz_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_Px_S_M5_vrr = QCX*I_ERI_S_S_S_S_M5_vrr+WQX*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Py_S_M5_vrr = QCY*I_ERI_S_S_S_S_M5_vrr+WQY*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Pz_S_M5_vrr = QCZ*I_ERI_S_S_S_S_M5_vrr+WQZ*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M5_vrr = QCX*I_ERI_S_S_Px_S_M5_vrr+WQX*I_ERI_S_S_Px_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Dxy_S_M5_vrr = QCY*I_ERI_S_S_Px_S_M5_vrr+WQY*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_D2y_S_M5_vrr = QCY*I_ERI_S_S_Py_S_M5_vrr+WQY*I_ERI_S_S_Py_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_D2z_S_M5_vrr = QCZ*I_ERI_S_S_Pz_S_M5_vrr+WQZ*I_ERI_S_S_Pz_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M5_vrr = QCX*I_ERI_S_S_D2x_S_M5_vrr+WQX*I_ERI_S_S_D2x_S_M6_vrr+2*oned2e*I_ERI_S_S_Px_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_F2xy_S_M5_vrr = QCY*I_ERI_S_S_D2x_S_M5_vrr+WQY*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_S_S_F2xz_S_M5_vrr = QCZ*I_ERI_S_S_D2x_S_M5_vrr+WQZ*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_S_S_Fx2y_S_M5_vrr = QCX*I_ERI_S_S_D2y_S_M5_vrr+WQX*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_S_S_Fx2z_S_M5_vrr = QCX*I_ERI_S_S_D2z_S_M5_vrr+WQX*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_S_S_F3y_S_M5_vrr = QCY*I_ERI_S_S_D2y_S_M5_vrr+WQY*I_ERI_S_S_D2y_S_M6_vrr+2*oned2e*I_ERI_S_S_Py_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_S_S_F2yz_S_M5_vrr = QCZ*I_ERI_S_S_D2y_S_M5_vrr+WQZ*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_S_S_F3z_S_M5_vrr = QCZ*I_ERI_S_S_D2z_S_M5_vrr+WQZ*I_ERI_S_S_D2z_S_M6_vrr+2*oned2e*I_ERI_S_S_Pz_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_Px_S_M4_vrr = QCX*I_ERI_S_S_S_S_M4_vrr+WQX*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Py_S_M4_vrr = QCY*I_ERI_S_S_S_S_M4_vrr+WQY*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Pz_S_M4_vrr = QCZ*I_ERI_S_S_S_S_M4_vrr+WQZ*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M4_vrr = QCX*I_ERI_S_S_Px_S_M4_vrr+WQX*I_ERI_S_S_Px_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dxy_S_M4_vrr = QCY*I_ERI_S_S_Px_S_M4_vrr+WQY*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_Dxz_S_M4_vrr = QCZ*I_ERI_S_S_Px_S_M4_vrr+WQZ*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_D2y_S_M4_vrr = QCY*I_ERI_S_S_Py_S_M4_vrr+WQY*I_ERI_S_S_Py_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dyz_S_M4_vrr = QCZ*I_ERI_S_S_Py_S_M4_vrr+WQZ*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_D2z_S_M4_vrr = QCZ*I_ERI_S_S_Pz_S_M4_vrr+WQZ*I_ERI_S_S_Pz_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M4_vrr = QCX*I_ERI_S_S_D2x_S_M4_vrr+WQX*I_ERI_S_S_D2x_S_M5_vrr+2*oned2e*I_ERI_S_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_F2xy_S_M4_vrr = QCY*I_ERI_S_S_D2x_S_M4_vrr+WQY*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_F2xz_S_M4_vrr = QCZ*I_ERI_S_S_D2x_S_M4_vrr+WQZ*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_Fx2y_S_M4_vrr = QCX*I_ERI_S_S_D2y_S_M4_vrr+WQX*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fxyz_S_M4_vrr = QCZ*I_ERI_S_S_Dxy_S_M4_vrr+WQZ*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_S_S_Fx2z_S_M4_vrr = QCX*I_ERI_S_S_D2z_S_M4_vrr+WQX*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3y_S_M4_vrr = QCY*I_ERI_S_S_D2y_S_M4_vrr+WQY*I_ERI_S_S_D2y_S_M5_vrr+2*oned2e*I_ERI_S_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_F2yz_S_M4_vrr = QCZ*I_ERI_S_S_D2y_S_M4_vrr+WQZ*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fy2z_S_M4_vrr = QCY*I_ERI_S_S_D2z_S_M4_vrr+WQY*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3z_S_M4_vrr = QCZ*I_ERI_S_S_D2z_S_M4_vrr+WQZ*I_ERI_S_S_D2z_S_M5_vrr+2*oned2e*I_ERI_S_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M4_vrr = QCX*I_ERI_S_S_F3x_S_M4_vrr+WQX*I_ERI_S_S_F3x_S_M5_vrr+3*oned2e*I_ERI_S_S_D2x_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_G3xy_S_M4_vrr = QCY*I_ERI_S_S_F3x_S_M4_vrr+WQY*I_ERI_S_S_F3x_S_M5_vrr;
      Double I_ERI_S_S_G3xz_S_M4_vrr = QCZ*I_ERI_S_S_F3x_S_M4_vrr+WQZ*I_ERI_S_S_F3x_S_M5_vrr;
      Double I_ERI_S_S_G2x2y_S_M4_vrr = QCY*I_ERI_S_S_F2xy_S_M4_vrr+WQY*I_ERI_S_S_F2xy_S_M5_vrr+oned2e*I_ERI_S_S_D2x_S_M4_vrr-rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_G2xyz_S_M4_vrr = QCZ*I_ERI_S_S_F2xy_S_M4_vrr+WQZ*I_ERI_S_S_F2xy_S_M5_vrr;
      Double I_ERI_S_S_G2x2z_S_M4_vrr = QCZ*I_ERI_S_S_F2xz_S_M4_vrr+WQZ*I_ERI_S_S_F2xz_S_M5_vrr+oned2e*I_ERI_S_S_D2x_S_M4_vrr-rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_Gx3y_S_M4_vrr = QCX*I_ERI_S_S_F3y_S_M4_vrr+WQX*I_ERI_S_S_F3y_S_M5_vrr;
      Double I_ERI_S_S_Gx2yz_S_M4_vrr = QCZ*I_ERI_S_S_Fx2y_S_M4_vrr+WQZ*I_ERI_S_S_Fx2y_S_M5_vrr;
      Double I_ERI_S_S_Gxy2z_S_M4_vrr = QCY*I_ERI_S_S_Fx2z_S_M4_vrr+WQY*I_ERI_S_S_Fx2z_S_M5_vrr;
      Double I_ERI_S_S_Gx3z_S_M4_vrr = QCX*I_ERI_S_S_F3z_S_M4_vrr+WQX*I_ERI_S_S_F3z_S_M5_vrr;
      Double I_ERI_S_S_G4y_S_M4_vrr = QCY*I_ERI_S_S_F3y_S_M4_vrr+WQY*I_ERI_S_S_F3y_S_M5_vrr+3*oned2e*I_ERI_S_S_D2y_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_G3yz_S_M4_vrr = QCZ*I_ERI_S_S_F3y_S_M4_vrr+WQZ*I_ERI_S_S_F3y_S_M5_vrr;
      Double I_ERI_S_S_G2y2z_S_M4_vrr = QCZ*I_ERI_S_S_F2yz_S_M4_vrr+WQZ*I_ERI_S_S_F2yz_S_M5_vrr+oned2e*I_ERI_S_S_D2y_S_M4_vrr-rhod2esq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Gy3z_S_M4_vrr = QCY*I_ERI_S_S_F3z_S_M4_vrr+WQY*I_ERI_S_S_F3z_S_M5_vrr;
      Double I_ERI_S_S_G4z_S_M4_vrr = QCZ*I_ERI_S_S_F3z_S_M4_vrr+WQZ*I_ERI_S_S_F3z_S_M5_vrr+3*oned2e*I_ERI_S_S_D2z_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_S_S_Px_S_M3_vrr = QCX*I_ERI_S_S_S_S_M3_vrr+WQX*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Py_S_M3_vrr = QCY*I_ERI_S_S_S_S_M3_vrr+WQY*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Pz_S_M3_vrr = QCZ*I_ERI_S_S_S_S_M3_vrr+WQZ*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M3_vrr = PAX*I_ERI_S_S_Px_S_M3_vrr+WPX*I_ERI_S_S_Px_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Py_S_Px_S_M3_vrr = PAY*I_ERI_S_S_Px_S_M3_vrr+WPY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Px_S_M3_vrr = PAZ*I_ERI_S_S_Px_S_M3_vrr+WPZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_Py_S_M3_vrr = PAX*I_ERI_S_S_Py_S_M3_vrr+WPX*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Py_S_M3_vrr = PAY*I_ERI_S_S_Py_S_M3_vrr+WPY*I_ERI_S_S_Py_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Pz_S_Py_S_M3_vrr = PAZ*I_ERI_S_S_Py_S_M3_vrr+WPZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_Pz_S_M3_vrr = PAX*I_ERI_S_S_Pz_S_M3_vrr+WPX*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Pz_S_M3_vrr = PAY*I_ERI_S_S_Pz_S_M3_vrr+WPY*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Pz_S_M3_vrr = PAZ*I_ERI_S_S_Pz_S_M3_vrr+WPZ*I_ERI_S_S_Pz_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M3_vrr = QCX*I_ERI_S_S_Px_S_M3_vrr+WQX*I_ERI_S_S_Px_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dxy_S_M3_vrr = QCY*I_ERI_S_S_Px_S_M3_vrr+WQY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_Dxz_S_M3_vrr = QCZ*I_ERI_S_S_Px_S_M3_vrr+WQZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_D2y_S_M3_vrr = QCY*I_ERI_S_S_Py_S_M3_vrr+WQY*I_ERI_S_S_Py_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dyz_S_M3_vrr = QCZ*I_ERI_S_S_Py_S_M3_vrr+WQZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_D2z_S_M3_vrr = QCZ*I_ERI_S_S_Pz_S_M3_vrr+WQZ*I_ERI_S_S_Pz_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M3_vrr = PAX*I_ERI_S_S_D2x_S_M3_vrr+WPX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Py_S_D2x_S_M3_vrr = PAY*I_ERI_S_S_D2x_S_M3_vrr+WPY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_D2x_S_M3_vrr = PAZ*I_ERI_S_S_D2x_S_M3_vrr+WPZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Dxy_S_M3_vrr = PAX*I_ERI_S_S_Dxy_S_M3_vrr+WPX*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Dxy_S_M3_vrr = PAY*I_ERI_S_S_Dxy_S_M3_vrr+WPY*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Dxy_S_M3_vrr = PAZ*I_ERI_S_S_Dxy_S_M3_vrr+WPZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Dxz_S_M3_vrr = PAX*I_ERI_S_S_Dxz_S_M3_vrr+WPX*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Dxz_S_M3_vrr = PAY*I_ERI_S_S_Dxz_S_M3_vrr+WPY*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Dxz_S_M3_vrr = PAZ*I_ERI_S_S_Dxz_S_M3_vrr+WPZ*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_D2y_S_M3_vrr = PAX*I_ERI_S_S_D2y_S_M3_vrr+WPX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_D2y_S_M3_vrr = PAY*I_ERI_S_S_D2y_S_M3_vrr+WPY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Pz_S_D2y_S_M3_vrr = PAZ*I_ERI_S_S_D2y_S_M3_vrr+WPZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Dyz_S_M3_vrr = PAX*I_ERI_S_S_Dyz_S_M3_vrr+WPX*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Dyz_S_M3_vrr = PAY*I_ERI_S_S_Dyz_S_M3_vrr+WPY*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Dyz_S_M3_vrr = PAZ*I_ERI_S_S_Dyz_S_M3_vrr+WPZ*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_D2z_S_M3_vrr = PAX*I_ERI_S_S_D2z_S_M3_vrr+WPX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_D2z_S_M3_vrr = PAY*I_ERI_S_S_D2z_S_M3_vrr+WPY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_D2z_S_M3_vrr = PAZ*I_ERI_S_S_D2z_S_M3_vrr+WPZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M3_vrr = QCX*I_ERI_S_S_D2x_S_M3_vrr+WQX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2e*I_ERI_S_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_F2xy_S_M3_vrr = QCY*I_ERI_S_S_D2x_S_M3_vrr+WQY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_F2xz_S_M3_vrr = QCZ*I_ERI_S_S_D2x_S_M3_vrr+WQZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_Fx2y_S_M3_vrr = QCX*I_ERI_S_S_D2y_S_M3_vrr+WQX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fxyz_S_M3_vrr = QCZ*I_ERI_S_S_Dxy_S_M3_vrr+WQZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_S_S_Fx2z_S_M3_vrr = QCX*I_ERI_S_S_D2z_S_M3_vrr+WQX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3y_S_M3_vrr = QCY*I_ERI_S_S_D2y_S_M3_vrr+WQY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2e*I_ERI_S_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_F2yz_S_M3_vrr = QCZ*I_ERI_S_S_D2y_S_M3_vrr+WQZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fy2z_S_M3_vrr = QCY*I_ERI_S_S_D2z_S_M3_vrr+WQY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3z_S_M3_vrr = QCZ*I_ERI_S_S_D2z_S_M3_vrr+WQZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2e*I_ERI_S_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M3_vrr = PAX*I_ERI_S_S_F3x_S_M3_vrr+WPX*I_ERI_S_S_F3x_S_M4_vrr+3*oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F3x_S_M3_vrr = PAY*I_ERI_S_S_F3x_S_M3_vrr+WPY*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Pz_S_F3x_S_M3_vrr = PAZ*I_ERI_S_S_F3x_S_M3_vrr+WPZ*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Px_S_F2xy_S_M3_vrr = PAX*I_ERI_S_S_F2xy_S_M3_vrr+WPX*I_ERI_S_S_F2xy_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Py_S_F2xy_S_M3_vrr = PAY*I_ERI_S_S_F2xy_S_M3_vrr+WPY*I_ERI_S_S_F2xy_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_F2xy_S_M3_vrr = PAZ*I_ERI_S_S_F2xy_S_M3_vrr+WPZ*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Px_S_F2xz_S_M3_vrr = PAX*I_ERI_S_S_F2xz_S_M3_vrr+WPX*I_ERI_S_S_F2xz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Py_S_F2xz_S_M3_vrr = PAY*I_ERI_S_S_F2xz_S_M3_vrr+WPY*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Pz_S_F2xz_S_M3_vrr = PAZ*I_ERI_S_S_F2xz_S_M3_vrr+WPZ*I_ERI_S_S_F2xz_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Fx2y_S_M3_vrr = PAX*I_ERI_S_S_Fx2y_S_M3_vrr+WPX*I_ERI_S_S_Fx2y_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_Fx2y_S_M3_vrr = PAY*I_ERI_S_S_Fx2y_S_M3_vrr+WPY*I_ERI_S_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M3_vrr = PAZ*I_ERI_S_S_Fx2y_S_M3_vrr+WPZ*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Px_S_Fxyz_S_M3_vrr = PAX*I_ERI_S_S_Fxyz_S_M3_vrr+WPX*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Fxyz_S_M3_vrr = PAY*I_ERI_S_S_Fxyz_S_M3_vrr+WPY*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M3_vrr = PAZ*I_ERI_S_S_Fxyz_S_M3_vrr+WPZ*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Fx2z_S_M3_vrr = PAX*I_ERI_S_S_Fx2z_S_M3_vrr+WPX*I_ERI_S_S_Fx2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_Fx2z_S_M3_vrr = PAY*I_ERI_S_S_Fx2z_S_M3_vrr+WPY*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M3_vrr = PAZ*I_ERI_S_S_Fx2z_S_M3_vrr+WPZ*I_ERI_S_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Px_S_F3y_S_M3_vrr = PAX*I_ERI_S_S_F3y_S_M3_vrr+WPX*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Py_S_F3y_S_M3_vrr = PAY*I_ERI_S_S_F3y_S_M3_vrr+WPY*I_ERI_S_S_F3y_S_M4_vrr+3*oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_F3y_S_M3_vrr = PAZ*I_ERI_S_S_F3y_S_M3_vrr+WPZ*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Px_S_F2yz_S_M3_vrr = PAX*I_ERI_S_S_F2yz_S_M3_vrr+WPX*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Py_S_F2yz_S_M3_vrr = PAY*I_ERI_S_S_F2yz_S_M3_vrr+WPY*I_ERI_S_S_F2yz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Pz_S_F2yz_S_M3_vrr = PAZ*I_ERI_S_S_F2yz_S_M3_vrr+WPZ*I_ERI_S_S_F2yz_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Fy2z_S_M3_vrr = PAX*I_ERI_S_S_Fy2z_S_M3_vrr+WPX*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Py_S_Fy2z_S_M3_vrr = PAY*I_ERI_S_S_Fy2z_S_M3_vrr+WPY*I_ERI_S_S_Fy2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M3_vrr = PAZ*I_ERI_S_S_Fy2z_S_M3_vrr+WPZ*I_ERI_S_S_Fy2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Px_S_F3z_S_M3_vrr = PAX*I_ERI_S_S_F3z_S_M3_vrr+WPX*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Py_S_F3z_S_M3_vrr = PAY*I_ERI_S_S_F3z_S_M3_vrr+WPY*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Pz_S_F3z_S_M3_vrr = PAZ*I_ERI_S_S_F3z_S_M3_vrr+WPZ*I_ERI_S_S_F3z_S_M4_vrr+3*oned2k*I_ERI_S_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M3_vrr = QCX*I_ERI_S_S_F3x_S_M3_vrr+WQX*I_ERI_S_S_F3x_S_M4_vrr+3*oned2e*I_ERI_S_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_G3xy_S_M3_vrr = QCY*I_ERI_S_S_F3x_S_M3_vrr+WQY*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_S_S_G3xz_S_M3_vrr = QCZ*I_ERI_S_S_F3x_S_M3_vrr+WQZ*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_S_S_G2x2y_S_M3_vrr = QCY*I_ERI_S_S_F2xy_S_M3_vrr+WQY*I_ERI_S_S_F2xy_S_M4_vrr+oned2e*I_ERI_S_S_D2x_S_M3_vrr-rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_G2xyz_S_M3_vrr = QCZ*I_ERI_S_S_F2xy_S_M3_vrr+WQZ*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_S_S_G2x2z_S_M3_vrr = QCZ*I_ERI_S_S_F2xz_S_M3_vrr+WQZ*I_ERI_S_S_F2xz_S_M4_vrr+oned2e*I_ERI_S_S_D2x_S_M3_vrr-rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_Gx3y_S_M3_vrr = QCX*I_ERI_S_S_F3y_S_M3_vrr+WQX*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_S_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_S_S_Fx2y_S_M3_vrr+WQZ*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_S_S_Gxy2z_S_M3_vrr = QCY*I_ERI_S_S_Fx2z_S_M3_vrr+WQY*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_S_S_Gx3z_S_M3_vrr = QCX*I_ERI_S_S_F3z_S_M3_vrr+WQX*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_S_S_G4y_S_M3_vrr = QCY*I_ERI_S_S_F3y_S_M3_vrr+WQY*I_ERI_S_S_F3y_S_M4_vrr+3*oned2e*I_ERI_S_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_G3yz_S_M3_vrr = QCZ*I_ERI_S_S_F3y_S_M3_vrr+WQZ*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_S_S_G2y2z_S_M3_vrr = QCZ*I_ERI_S_S_F2yz_S_M3_vrr+WQZ*I_ERI_S_S_F2yz_S_M4_vrr+oned2e*I_ERI_S_S_D2y_S_M3_vrr-rhod2esq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Gy3z_S_M3_vrr = QCY*I_ERI_S_S_F3z_S_M3_vrr+WQY*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_S_S_G4z_S_M3_vrr = QCZ*I_ERI_S_S_F3z_S_M3_vrr+WQZ*I_ERI_S_S_F3z_S_M4_vrr+3*oned2e*I_ERI_S_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M3_vrr = PAX*I_ERI_S_S_G4x_S_M3_vrr+WPX*I_ERI_S_S_G4x_S_M4_vrr+4*oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Py_S_G4x_S_M3_vrr = PAY*I_ERI_S_S_G4x_S_M3_vrr+WPY*I_ERI_S_S_G4x_S_M4_vrr;
      Double I_ERI_Pz_S_G4x_S_M3_vrr = PAZ*I_ERI_S_S_G4x_S_M3_vrr+WPZ*I_ERI_S_S_G4x_S_M4_vrr;
      Double I_ERI_Px_S_G3xy_S_M3_vrr = PAX*I_ERI_S_S_G3xy_S_M3_vrr+WPX*I_ERI_S_S_G3xy_S_M4_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Py_S_G3xy_S_M3_vrr = PAY*I_ERI_S_S_G3xy_S_M3_vrr+WPY*I_ERI_S_S_G3xy_S_M4_vrr+oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Pz_S_G3xy_S_M3_vrr = PAZ*I_ERI_S_S_G3xy_S_M3_vrr+WPZ*I_ERI_S_S_G3xy_S_M4_vrr;
      Double I_ERI_Px_S_G3xz_S_M3_vrr = PAX*I_ERI_S_S_G3xz_S_M3_vrr+WPX*I_ERI_S_S_G3xz_S_M4_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Py_S_G3xz_S_M3_vrr = PAY*I_ERI_S_S_G3xz_S_M3_vrr+WPY*I_ERI_S_S_G3xz_S_M4_vrr;
      Double I_ERI_Pz_S_G3xz_S_M3_vrr = PAZ*I_ERI_S_S_G3xz_S_M3_vrr+WPZ*I_ERI_S_S_G3xz_S_M4_vrr+oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Px_S_G2x2y_S_M3_vrr = PAX*I_ERI_S_S_G2x2y_S_M3_vrr+WPX*I_ERI_S_S_G2x2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Py_S_G2x2y_S_M3_vrr = PAY*I_ERI_S_S_G2x2y_S_M3_vrr+WPY*I_ERI_S_S_G2x2y_S_M4_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M3_vrr = PAZ*I_ERI_S_S_G2x2y_S_M3_vrr+WPZ*I_ERI_S_S_G2x2y_S_M4_vrr;
      Double I_ERI_Px_S_G2xyz_S_M3_vrr = PAX*I_ERI_S_S_G2xyz_S_M3_vrr+WPX*I_ERI_S_S_G2xyz_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Py_S_G2xyz_S_M3_vrr = PAY*I_ERI_S_S_G2xyz_S_M3_vrr+WPY*I_ERI_S_S_G2xyz_S_M4_vrr+oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M3_vrr = PAZ*I_ERI_S_S_G2xyz_S_M3_vrr+WPZ*I_ERI_S_S_G2xyz_S_M4_vrr+oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Px_S_G2x2z_S_M3_vrr = PAX*I_ERI_S_S_G2x2z_S_M3_vrr+WPX*I_ERI_S_S_G2x2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Py_S_G2x2z_S_M3_vrr = PAY*I_ERI_S_S_G2x2z_S_M3_vrr+WPY*I_ERI_S_S_G2x2z_S_M4_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M3_vrr = PAZ*I_ERI_S_S_G2x2z_S_M3_vrr+WPZ*I_ERI_S_S_G2x2z_S_M4_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Px_S_Gx3y_S_M3_vrr = PAX*I_ERI_S_S_Gx3y_S_M3_vrr+WPX*I_ERI_S_S_Gx3y_S_M4_vrr+oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Py_S_Gx3y_S_M3_vrr = PAY*I_ERI_S_S_Gx3y_S_M3_vrr+WPY*I_ERI_S_S_Gx3y_S_M4_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M3_vrr = PAZ*I_ERI_S_S_Gx3y_S_M3_vrr+WPZ*I_ERI_S_S_Gx3y_S_M4_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M3_vrr = PAX*I_ERI_S_S_Gx2yz_S_M3_vrr+WPX*I_ERI_S_S_Gx2yz_S_M4_vrr+oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M3_vrr = PAY*I_ERI_S_S_Gx2yz_S_M3_vrr+WPY*I_ERI_S_S_Gx2yz_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M3_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M3_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M4_vrr+oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M3_vrr = PAX*I_ERI_S_S_Gxy2z_S_M3_vrr+WPX*I_ERI_S_S_Gxy2z_S_M4_vrr+oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M3_vrr = PAY*I_ERI_S_S_Gxy2z_S_M3_vrr+WPY*I_ERI_S_S_Gxy2z_S_M4_vrr+oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M3_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M3_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Px_S_Gx3z_S_M3_vrr = PAX*I_ERI_S_S_Gx3z_S_M3_vrr+WPX*I_ERI_S_S_Gx3z_S_M4_vrr+oned2k*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Py_S_Gx3z_S_M3_vrr = PAY*I_ERI_S_S_Gx3z_S_M3_vrr+WPY*I_ERI_S_S_Gx3z_S_M4_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M3_vrr = PAZ*I_ERI_S_S_Gx3z_S_M3_vrr+WPZ*I_ERI_S_S_Gx3z_S_M4_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Px_S_G4y_S_M3_vrr = PAX*I_ERI_S_S_G4y_S_M3_vrr+WPX*I_ERI_S_S_G4y_S_M4_vrr;
      Double I_ERI_Py_S_G4y_S_M3_vrr = PAY*I_ERI_S_S_G4y_S_M3_vrr+WPY*I_ERI_S_S_G4y_S_M4_vrr+4*oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Pz_S_G4y_S_M3_vrr = PAZ*I_ERI_S_S_G4y_S_M3_vrr+WPZ*I_ERI_S_S_G4y_S_M4_vrr;
      Double I_ERI_Px_S_G3yz_S_M3_vrr = PAX*I_ERI_S_S_G3yz_S_M3_vrr+WPX*I_ERI_S_S_G3yz_S_M4_vrr;
      Double I_ERI_Py_S_G3yz_S_M3_vrr = PAY*I_ERI_S_S_G3yz_S_M3_vrr+WPY*I_ERI_S_S_G3yz_S_M4_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Pz_S_G3yz_S_M3_vrr = PAZ*I_ERI_S_S_G3yz_S_M3_vrr+WPZ*I_ERI_S_S_G3yz_S_M4_vrr+oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Px_S_G2y2z_S_M3_vrr = PAX*I_ERI_S_S_G2y2z_S_M3_vrr+WPX*I_ERI_S_S_G2y2z_S_M4_vrr;
      Double I_ERI_Py_S_G2y2z_S_M3_vrr = PAY*I_ERI_S_S_G2y2z_S_M3_vrr+WPY*I_ERI_S_S_G2y2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M3_vrr = PAZ*I_ERI_S_S_G2y2z_S_M3_vrr+WPZ*I_ERI_S_S_G2y2z_S_M4_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Px_S_Gy3z_S_M3_vrr = PAX*I_ERI_S_S_Gy3z_S_M3_vrr+WPX*I_ERI_S_S_Gy3z_S_M4_vrr;
      Double I_ERI_Py_S_Gy3z_S_M3_vrr = PAY*I_ERI_S_S_Gy3z_S_M3_vrr+WPY*I_ERI_S_S_Gy3z_S_M4_vrr+oned2k*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M3_vrr = PAZ*I_ERI_S_S_Gy3z_S_M3_vrr+WPZ*I_ERI_S_S_Gy3z_S_M4_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Px_S_G4z_S_M3_vrr = PAX*I_ERI_S_S_G4z_S_M3_vrr+WPX*I_ERI_S_S_G4z_S_M4_vrr;
      Double I_ERI_Py_S_G4z_S_M3_vrr = PAY*I_ERI_S_S_G4z_S_M3_vrr+WPY*I_ERI_S_S_G4z_S_M4_vrr;
      Double I_ERI_Pz_S_G4z_S_M3_vrr = PAZ*I_ERI_S_S_G4z_S_M3_vrr+WPZ*I_ERI_S_S_G4z_S_M4_vrr+4*oned2k*I_ERI_S_S_F3z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_S_Px_S_M2_vrr = QCX*I_ERI_S_S_S_S_M2_vrr+WQX*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Py_S_M2_vrr = QCY*I_ERI_S_S_S_S_M2_vrr+WQY*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Pz_S_M2_vrr = QCZ*I_ERI_S_S_S_S_M2_vrr+WQZ*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M2_vrr = QCX*I_ERI_S_S_Px_S_M2_vrr+WQX*I_ERI_S_S_Px_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dxy_S_M2_vrr = QCY*I_ERI_S_S_Px_S_M2_vrr+WQY*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_Dxz_S_M2_vrr = QCZ*I_ERI_S_S_Px_S_M2_vrr+WQZ*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_D2y_S_M2_vrr = QCY*I_ERI_S_S_Py_S_M2_vrr+WQY*I_ERI_S_S_Py_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dyz_S_M2_vrr = QCZ*I_ERI_S_S_Py_S_M2_vrr+WQZ*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_D2z_S_M2_vrr = QCZ*I_ERI_S_S_Pz_S_M2_vrr+WQZ*I_ERI_S_S_Pz_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M2_vrr = PAX*I_ERI_S_S_D2x_S_M2_vrr+WPX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Py_S_D2x_S_M2_vrr = PAY*I_ERI_S_S_D2x_S_M2_vrr+WPY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_D2x_S_M2_vrr = PAZ*I_ERI_S_S_D2x_S_M2_vrr+WPZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Dxy_S_M2_vrr = PAX*I_ERI_S_S_Dxy_S_M2_vrr+WPX*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Py_S_Dxy_S_M2_vrr = PAY*I_ERI_S_S_Dxy_S_M2_vrr+WPY*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Pz_S_Dxy_S_M2_vrr = PAZ*I_ERI_S_S_Dxy_S_M2_vrr+WPZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Dxz_S_M2_vrr = PAX*I_ERI_S_S_Dxz_S_M2_vrr+WPX*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Py_S_Dxz_S_M2_vrr = PAY*I_ERI_S_S_Dxz_S_M2_vrr+WPY*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Dxz_S_M2_vrr = PAZ*I_ERI_S_S_Dxz_S_M2_vrr+WPZ*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Px_S_D2y_S_M2_vrr = PAX*I_ERI_S_S_D2y_S_M2_vrr+WPX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_D2y_S_M2_vrr = PAY*I_ERI_S_S_D2y_S_M2_vrr+WPY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Pz_S_D2y_S_M2_vrr = PAZ*I_ERI_S_S_D2y_S_M2_vrr+WPZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Dyz_S_M2_vrr = PAX*I_ERI_S_S_Dyz_S_M2_vrr+WPX*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Dyz_S_M2_vrr = PAY*I_ERI_S_S_Dyz_S_M2_vrr+WPY*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Pz_S_Dyz_S_M2_vrr = PAZ*I_ERI_S_S_Dyz_S_M2_vrr+WPZ*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Px_S_D2z_S_M2_vrr = PAX*I_ERI_S_S_D2z_S_M2_vrr+WPX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_D2z_S_M2_vrr = PAY*I_ERI_S_S_D2z_S_M2_vrr+WPY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_D2z_S_M2_vrr = PAZ*I_ERI_S_S_D2z_S_M2_vrr+WPZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M2_vrr = QCX*I_ERI_S_S_D2x_S_M2_vrr+WQX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2e*I_ERI_S_S_Px_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_F2xy_S_M2_vrr = QCY*I_ERI_S_S_D2x_S_M2_vrr+WQY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_F2xz_S_M2_vrr = QCZ*I_ERI_S_S_D2x_S_M2_vrr+WQZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_Fx2y_S_M2_vrr = QCX*I_ERI_S_S_D2y_S_M2_vrr+WQX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fxyz_S_M2_vrr = QCZ*I_ERI_S_S_Dxy_S_M2_vrr+WQZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_S_S_Fx2z_S_M2_vrr = QCX*I_ERI_S_S_D2z_S_M2_vrr+WQX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3y_S_M2_vrr = QCY*I_ERI_S_S_D2y_S_M2_vrr+WQY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2e*I_ERI_S_S_Py_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_F2yz_S_M2_vrr = QCZ*I_ERI_S_S_D2y_S_M2_vrr+WQZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fy2z_S_M2_vrr = QCY*I_ERI_S_S_D2z_S_M2_vrr+WQY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3z_S_M2_vrr = QCZ*I_ERI_S_S_D2z_S_M2_vrr+WQZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2e*I_ERI_S_S_Pz_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 18 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M2_vrr = PAX*I_ERI_Px_S_D2x_S_M2_vrr+WPX*I_ERI_Px_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Px_S_Px_S_M3_vrr;
      Double I_ERI_D2y_S_D2x_S_M2_vrr = PAY*I_ERI_Py_S_D2x_S_M2_vrr+WPY*I_ERI_Py_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_D2x_S_M2_vrr = PAZ*I_ERI_Pz_S_D2x_S_M2_vrr+WPZ*I_ERI_Pz_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Dxy_S_M2_vrr = PAX*I_ERI_Px_S_Dxy_S_M2_vrr+WPX*I_ERI_Px_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Dxy_S_M2_vrr = PAY*I_ERI_Py_S_Dxy_S_M2_vrr+WPY*I_ERI_Py_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Dxy_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Dxz_S_M2_vrr = PAX*I_ERI_Px_S_Dxz_S_M2_vrr+WPX*I_ERI_Px_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_D2y_S_Dxz_S_M2_vrr = PAY*I_ERI_Py_S_Dxz_S_M2_vrr+WPY*I_ERI_Py_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Dxz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_D2y_S_M2_vrr = PAX*I_ERI_Px_S_D2y_S_M2_vrr+WPX*I_ERI_Px_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_D2y_S_M2_vrr = PAY*I_ERI_Py_S_D2y_S_M2_vrr+WPY*I_ERI_Py_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Py_S_M3_vrr;
      Double I_ERI_D2z_S_D2y_S_M2_vrr = PAZ*I_ERI_Pz_S_D2y_S_M2_vrr+WPZ*I_ERI_Pz_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Dyz_S_M2_vrr = PAX*I_ERI_Px_S_Dyz_S_M2_vrr+WPX*I_ERI_Px_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Dyz_S_M2_vrr = PAY*I_ERI_Py_S_Dyz_S_M2_vrr+WPY*I_ERI_Py_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Py_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_Dyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_D2z_S_M2_vrr = PAX*I_ERI_Px_S_D2z_S_M2_vrr+WPX*I_ERI_Px_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_D2z_S_M2_vrr = PAY*I_ERI_Py_S_D2z_S_M2_vrr+WPY*I_ERI_Py_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_D2z_S_M2_vrr = PAZ*I_ERI_Pz_S_D2z_S_M2_vrr+WPZ*I_ERI_Pz_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M2_vrr = PAX*I_ERI_S_S_F3x_S_M2_vrr+WPX*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Py_S_F3x_S_M2_vrr = PAY*I_ERI_S_S_F3x_S_M2_vrr+WPY*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_F3x_S_M2_vrr = PAZ*I_ERI_S_S_F3x_S_M2_vrr+WPZ*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_F2xy_S_M2_vrr = PAX*I_ERI_S_S_F2xy_S_M2_vrr+WPX*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Py_S_F2xy_S_M2_vrr = PAY*I_ERI_S_S_F2xy_S_M2_vrr+WPY*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_F2xy_S_M2_vrr = PAZ*I_ERI_S_S_F2xy_S_M2_vrr+WPZ*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_F2xz_S_M2_vrr = PAX*I_ERI_S_S_F2xz_S_M2_vrr+WPX*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Py_S_F2xz_S_M2_vrr = PAY*I_ERI_S_S_F2xz_S_M2_vrr+WPY*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Pz_S_F2xz_S_M2_vrr = PAZ*I_ERI_S_S_F2xz_S_M2_vrr+WPZ*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Fx2y_S_M2_vrr = PAX*I_ERI_S_S_Fx2y_S_M2_vrr+WPX*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_Fx2y_S_M2_vrr = PAY*I_ERI_S_S_Fx2y_S_M2_vrr+WPY*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_S_S_Fx2y_S_M2_vrr+WPZ*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Fxyz_S_M2_vrr = PAX*I_ERI_S_S_Fxyz_S_M2_vrr+WPX*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Fxyz_S_M2_vrr = PAY*I_ERI_S_S_Fxyz_S_M2_vrr+WPY*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_S_S_Fxyz_S_M2_vrr+WPZ*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Fx2z_S_M2_vrr = PAX*I_ERI_S_S_Fx2z_S_M2_vrr+WPX*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_Fx2z_S_M2_vrr = PAY*I_ERI_S_S_Fx2z_S_M2_vrr+WPY*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_S_S_Fx2z_S_M2_vrr+WPZ*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Px_S_F3y_S_M2_vrr = PAX*I_ERI_S_S_F3y_S_M2_vrr+WPX*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_F3y_S_M2_vrr = PAY*I_ERI_S_S_F3y_S_M2_vrr+WPY*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Pz_S_F3y_S_M2_vrr = PAZ*I_ERI_S_S_F3y_S_M2_vrr+WPZ*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_F2yz_S_M2_vrr = PAX*I_ERI_S_S_F2yz_S_M2_vrr+WPX*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Py_S_F2yz_S_M2_vrr = PAY*I_ERI_S_S_F2yz_S_M2_vrr+WPY*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Pz_S_F2yz_S_M2_vrr = PAZ*I_ERI_S_S_F2yz_S_M2_vrr+WPZ*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Fy2z_S_M2_vrr = PAX*I_ERI_S_S_Fy2z_S_M2_vrr+WPX*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Py_S_Fy2z_S_M2_vrr = PAY*I_ERI_S_S_Fy2z_S_M2_vrr+WPY*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_S_S_Fy2z_S_M2_vrr+WPZ*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Px_S_F3z_S_M2_vrr = PAX*I_ERI_S_S_F3z_S_M2_vrr+WPX*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_F3z_S_M2_vrr = PAY*I_ERI_S_S_F3z_S_M2_vrr+WPY*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_F3z_S_M2_vrr = PAZ*I_ERI_S_S_F3z_S_M2_vrr+WPZ*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_S_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M2_vrr = QCX*I_ERI_S_S_F3x_S_M2_vrr+WQX*I_ERI_S_S_F3x_S_M3_vrr+3*oned2e*I_ERI_S_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_G3xy_S_M2_vrr = QCY*I_ERI_S_S_F3x_S_M2_vrr+WQY*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_S_S_G3xz_S_M2_vrr = QCZ*I_ERI_S_S_F3x_S_M2_vrr+WQZ*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_S_S_G2x2y_S_M2_vrr = QCY*I_ERI_S_S_F2xy_S_M2_vrr+WQY*I_ERI_S_S_F2xy_S_M3_vrr+oned2e*I_ERI_S_S_D2x_S_M2_vrr-rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_G2xyz_S_M2_vrr = QCZ*I_ERI_S_S_F2xy_S_M2_vrr+WQZ*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_S_S_G2x2z_S_M2_vrr = QCZ*I_ERI_S_S_F2xz_S_M2_vrr+WQZ*I_ERI_S_S_F2xz_S_M3_vrr+oned2e*I_ERI_S_S_D2x_S_M2_vrr-rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_Gx3y_S_M2_vrr = QCX*I_ERI_S_S_F3y_S_M2_vrr+WQX*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_S_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_S_S_Fx2y_S_M2_vrr+WQZ*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_S_S_Gxy2z_S_M2_vrr = QCY*I_ERI_S_S_Fx2z_S_M2_vrr+WQY*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_S_S_Gx3z_S_M2_vrr = QCX*I_ERI_S_S_F3z_S_M2_vrr+WQX*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_S_S_G4y_S_M2_vrr = QCY*I_ERI_S_S_F3y_S_M2_vrr+WQY*I_ERI_S_S_F3y_S_M3_vrr+3*oned2e*I_ERI_S_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_G3yz_S_M2_vrr = QCZ*I_ERI_S_S_F3y_S_M2_vrr+WQZ*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_S_S_G2y2z_S_M2_vrr = QCZ*I_ERI_S_S_F2yz_S_M2_vrr+WQZ*I_ERI_S_S_F2yz_S_M3_vrr+oned2e*I_ERI_S_S_D2y_S_M2_vrr-rhod2esq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Gy3z_S_M2_vrr = QCY*I_ERI_S_S_F3z_S_M2_vrr+WQY*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_S_S_G4z_S_M2_vrr = QCZ*I_ERI_S_S_F3z_S_M2_vrr+WQZ*I_ERI_S_S_F3z_S_M3_vrr+3*oned2e*I_ERI_S_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M2_vrr = PAX*I_ERI_Px_S_F3x_S_M2_vrr+WPX*I_ERI_Px_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_F3x_S_M2_vrr = PAY*I_ERI_Py_S_F3x_S_M2_vrr+WPY*I_ERI_Py_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_F3x_S_M2_vrr = PAZ*I_ERI_Pz_S_F3x_S_M2_vrr+WPZ*I_ERI_Pz_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_F2xy_S_M2_vrr = PAX*I_ERI_Px_S_F2xy_S_M2_vrr+WPX*I_ERI_Px_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M3_vrr;
      Double I_ERI_D2y_S_F2xy_S_M2_vrr = PAY*I_ERI_Py_S_F2xy_S_M2_vrr+WPY*I_ERI_Py_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_F2xy_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_F2xz_S_M2_vrr = PAX*I_ERI_Px_S_F2xz_S_M2_vrr+WPX*I_ERI_Px_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M3_vrr;
      Double I_ERI_D2y_S_F2xz_S_M2_vrr = PAY*I_ERI_Py_S_F2xz_S_M2_vrr+WPY*I_ERI_Py_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_D2z_S_F2xz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M2_vrr = PAX*I_ERI_Px_S_Fx2y_S_M2_vrr+WPX*I_ERI_Px_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M2_vrr = PAY*I_ERI_Py_S_Fx2y_S_M2_vrr+WPY*I_ERI_Py_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M2_vrr = PAX*I_ERI_Px_S_Fxyz_S_M2_vrr+WPX*I_ERI_Px_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Px_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M2_vrr = PAY*I_ERI_Py_S_Fxyz_S_M2_vrr+WPY*I_ERI_Py_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Py_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M2_vrr = PAX*I_ERI_Px_S_Fx2z_S_M2_vrr+WPX*I_ERI_Px_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M2_vrr = PAY*I_ERI_Py_S_Fx2z_S_M2_vrr+WPY*I_ERI_Py_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M3_vrr;
      Double I_ERI_D2x_S_F3y_S_M2_vrr = PAX*I_ERI_Px_S_F3y_S_M2_vrr+WPX*I_ERI_Px_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_F3y_S_M2_vrr = PAY*I_ERI_Py_S_F3y_S_M2_vrr+WPY*I_ERI_Py_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M3_vrr;
      Double I_ERI_D2z_S_F3y_S_M2_vrr = PAZ*I_ERI_Pz_S_F3y_S_M2_vrr+WPZ*I_ERI_Pz_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_F2yz_S_M2_vrr = PAX*I_ERI_Px_S_F2yz_S_M2_vrr+WPX*I_ERI_Px_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_D2y_S_F2yz_S_M2_vrr = PAY*I_ERI_Py_S_F2yz_S_M2_vrr+WPY*I_ERI_Py_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M3_vrr;
      Double I_ERI_D2z_S_F2yz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M2_vrr = PAX*I_ERI_Px_S_Fy2z_S_M2_vrr+WPX*I_ERI_Px_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M2_vrr = PAY*I_ERI_Py_S_Fy2z_S_M2_vrr+WPY*I_ERI_Py_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_Py_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M3_vrr;
      Double I_ERI_D2x_S_F3z_S_M2_vrr = PAX*I_ERI_Px_S_F3z_S_M2_vrr+WPX*I_ERI_Px_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_F3z_S_M2_vrr = PAY*I_ERI_Py_S_F3z_S_M2_vrr+WPY*I_ERI_Py_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_F3z_S_M2_vrr = PAZ*I_ERI_Pz_S_F3z_S_M2_vrr+WPZ*I_ERI_Pz_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M2_vrr = PAX*I_ERI_S_S_G4x_S_M2_vrr+WPX*I_ERI_S_S_G4x_S_M3_vrr+4*oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Py_S_G4x_S_M2_vrr = PAY*I_ERI_S_S_G4x_S_M2_vrr+WPY*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_Pz_S_G4x_S_M2_vrr = PAZ*I_ERI_S_S_G4x_S_M2_vrr+WPZ*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_Px_S_G3xy_S_M2_vrr = PAX*I_ERI_S_S_G3xy_S_M2_vrr+WPX*I_ERI_S_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Py_S_G3xy_S_M2_vrr = PAY*I_ERI_S_S_G3xy_S_M2_vrr+WPY*I_ERI_S_S_G3xy_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_G3xy_S_M2_vrr = PAZ*I_ERI_S_S_G3xy_S_M2_vrr+WPZ*I_ERI_S_S_G3xy_S_M3_vrr;
      Double I_ERI_Px_S_G3xz_S_M2_vrr = PAX*I_ERI_S_S_G3xz_S_M2_vrr+WPX*I_ERI_S_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Py_S_G3xz_S_M2_vrr = PAY*I_ERI_S_S_G3xz_S_M2_vrr+WPY*I_ERI_S_S_G3xz_S_M3_vrr;
      Double I_ERI_Pz_S_G3xz_S_M2_vrr = PAZ*I_ERI_S_S_G3xz_S_M2_vrr+WPZ*I_ERI_S_S_G3xz_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_G2x2y_S_M2_vrr = PAX*I_ERI_S_S_G2x2y_S_M2_vrr+WPX*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Py_S_G2x2y_S_M2_vrr = PAY*I_ERI_S_S_G2x2y_S_M2_vrr+WPY*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M2_vrr = PAZ*I_ERI_S_S_G2x2y_S_M2_vrr+WPZ*I_ERI_S_S_G2x2y_S_M3_vrr;
      Double I_ERI_Px_S_G2xyz_S_M2_vrr = PAX*I_ERI_S_S_G2xyz_S_M2_vrr+WPX*I_ERI_S_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Py_S_G2xyz_S_M2_vrr = PAY*I_ERI_S_S_G2xyz_S_M2_vrr+WPY*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M2_vrr = PAZ*I_ERI_S_S_G2xyz_S_M2_vrr+WPZ*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_G2x2z_S_M2_vrr = PAX*I_ERI_S_S_G2x2z_S_M2_vrr+WPX*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Py_S_G2x2z_S_M2_vrr = PAY*I_ERI_S_S_G2x2z_S_M2_vrr+WPY*I_ERI_S_S_G2x2z_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M2_vrr = PAZ*I_ERI_S_S_G2x2z_S_M2_vrr+WPZ*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Px_S_Gx3y_S_M2_vrr = PAX*I_ERI_S_S_Gx3y_S_M2_vrr+WPX*I_ERI_S_S_Gx3y_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_Gx3y_S_M2_vrr = PAY*I_ERI_S_S_Gx3y_S_M2_vrr+WPY*I_ERI_S_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M2_vrr = PAZ*I_ERI_S_S_Gx3y_S_M2_vrr+WPZ*I_ERI_S_S_Gx3y_S_M3_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M2_vrr = PAX*I_ERI_S_S_Gx2yz_S_M2_vrr+WPX*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M2_vrr = PAY*I_ERI_S_S_Gx2yz_S_M2_vrr+WPY*I_ERI_S_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M2_vrr = PAX*I_ERI_S_S_Gxy2z_S_M2_vrr+WPX*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M2_vrr = PAY*I_ERI_S_S_Gxy2z_S_M2_vrr+WPY*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Px_S_Gx3z_S_M2_vrr = PAX*I_ERI_S_S_Gx3z_S_M2_vrr+WPX*I_ERI_S_S_Gx3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_Gx3z_S_M2_vrr = PAY*I_ERI_S_S_Gx3z_S_M2_vrr+WPY*I_ERI_S_S_Gx3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M2_vrr = PAZ*I_ERI_S_S_Gx3z_S_M2_vrr+WPZ*I_ERI_S_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Px_S_G4y_S_M2_vrr = PAX*I_ERI_S_S_G4y_S_M2_vrr+WPX*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_Py_S_G4y_S_M2_vrr = PAY*I_ERI_S_S_G4y_S_M2_vrr+WPY*I_ERI_S_S_G4y_S_M3_vrr+4*oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Pz_S_G4y_S_M2_vrr = PAZ*I_ERI_S_S_G4y_S_M2_vrr+WPZ*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_Px_S_G3yz_S_M2_vrr = PAX*I_ERI_S_S_G3yz_S_M2_vrr+WPX*I_ERI_S_S_G3yz_S_M3_vrr;
      Double I_ERI_Py_S_G3yz_S_M2_vrr = PAY*I_ERI_S_S_G3yz_S_M2_vrr+WPY*I_ERI_S_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Pz_S_G3yz_S_M2_vrr = PAZ*I_ERI_S_S_G3yz_S_M2_vrr+WPZ*I_ERI_S_S_G3yz_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_G2y2z_S_M2_vrr = PAX*I_ERI_S_S_G2y2z_S_M2_vrr+WPX*I_ERI_S_S_G2y2z_S_M3_vrr;
      Double I_ERI_Py_S_G2y2z_S_M2_vrr = PAY*I_ERI_S_S_G2y2z_S_M2_vrr+WPY*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M2_vrr = PAZ*I_ERI_S_S_G2y2z_S_M2_vrr+WPZ*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Px_S_Gy3z_S_M2_vrr = PAX*I_ERI_S_S_Gy3z_S_M2_vrr+WPX*I_ERI_S_S_Gy3z_S_M3_vrr;
      Double I_ERI_Py_S_Gy3z_S_M2_vrr = PAY*I_ERI_S_S_Gy3z_S_M2_vrr+WPY*I_ERI_S_S_Gy3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M2_vrr = PAZ*I_ERI_S_S_Gy3z_S_M2_vrr+WPZ*I_ERI_S_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Px_S_G4z_S_M2_vrr = PAX*I_ERI_S_S_G4z_S_M2_vrr+WPX*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_Py_S_G4z_S_M2_vrr = PAY*I_ERI_S_S_G4z_S_M2_vrr+WPY*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_Pz_S_G4z_S_M2_vrr = PAZ*I_ERI_S_S_G4z_S_M2_vrr+WPZ*I_ERI_S_S_G4z_S_M3_vrr+4*oned2k*I_ERI_S_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M2_vrr = PAX*I_ERI_Px_S_G4x_S_M2_vrr+WPX*I_ERI_Px_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_G4x_S_M2_vrr = PAY*I_ERI_Py_S_G4x_S_M2_vrr+WPY*I_ERI_Py_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_D2z_S_G4x_S_M2_vrr = PAZ*I_ERI_Pz_S_G4x_S_M2_vrr+WPZ*I_ERI_Pz_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_D2x_S_G3xy_S_M2_vrr = PAX*I_ERI_Px_S_G3xy_S_M2_vrr+WPX*I_ERI_Px_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M3_vrr;
      Double I_ERI_D2y_S_G3xy_S_M2_vrr = PAY*I_ERI_Py_S_G3xy_S_M2_vrr+WPY*I_ERI_Py_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr+oned2k*I_ERI_Py_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_G3xy_S_M2_vrr = PAZ*I_ERI_Pz_S_G3xy_S_M2_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr;
      Double I_ERI_D2x_S_G3xz_S_M2_vrr = PAX*I_ERI_Px_S_G3xz_S_M2_vrr+WPX*I_ERI_Px_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M3_vrr;
      Double I_ERI_D2y_S_G3xz_S_M2_vrr = PAY*I_ERI_Py_S_G3xz_S_M2_vrr+WPY*I_ERI_Py_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr;
      Double I_ERI_D2z_S_G3xz_S_M2_vrr = PAZ*I_ERI_Pz_S_G3xz_S_M2_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr+oned2k*I_ERI_Pz_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M2_vrr = PAX*I_ERI_Px_S_G2x2y_S_M2_vrr+WPX*I_ERI_Px_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M2_vrr = PAY*I_ERI_Py_S_G2x2y_S_M2_vrr+WPY*I_ERI_Py_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M2_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_M2_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M2_vrr = PAX*I_ERI_Px_S_G2xyz_S_M2_vrr+WPX*I_ERI_Px_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M2_vrr = PAY*I_ERI_Py_S_G2xyz_S_M2_vrr+WPY*I_ERI_Py_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_Py_S_F2xz_S_M3_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M2_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_M2_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M2_vrr = PAX*I_ERI_Px_S_G2x2z_S_M2_vrr+WPX*I_ERI_Px_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M2_vrr = PAY*I_ERI_Py_S_G2x2z_S_M2_vrr+WPY*I_ERI_Py_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M2_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_M2_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M2_vrr = PAX*I_ERI_Px_S_Gx3y_S_M2_vrr+WPX*I_ERI_Px_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr+oned2k*I_ERI_Px_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M2_vrr = PAY*I_ERI_Py_S_Gx3y_S_M2_vrr+WPY*I_ERI_Py_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M2_vrr = PAX*I_ERI_Px_S_Gx2yz_S_M2_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_Px_S_F2yz_S_M3_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M2_vrr = PAY*I_ERI_Py_S_Gx2yz_S_M2_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M2_vrr = PAX*I_ERI_Px_S_Gxy2z_S_M2_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M2_vrr = PAY*I_ERI_Py_S_Gxy2z_S_M2_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M2_vrr = PAX*I_ERI_Px_S_Gx3z_S_M2_vrr+WPX*I_ERI_Px_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr+oned2k*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M2_vrr = PAY*I_ERI_Py_S_Gx3z_S_M2_vrr+WPY*I_ERI_Py_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2x_S_G4y_S_M2_vrr = PAX*I_ERI_Px_S_G4y_S_M2_vrr+WPX*I_ERI_Px_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_D2y_S_G4y_S_M2_vrr = PAY*I_ERI_Py_S_G4y_S_M2_vrr+WPY*I_ERI_Py_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M3_vrr;
      Double I_ERI_D2z_S_G4y_S_M2_vrr = PAZ*I_ERI_Pz_S_G4y_S_M2_vrr+WPZ*I_ERI_Pz_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_D2x_S_G3yz_S_M2_vrr = PAX*I_ERI_Px_S_G3yz_S_M2_vrr+WPX*I_ERI_Px_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr;
      Double I_ERI_D2y_S_G3yz_S_M2_vrr = PAY*I_ERI_Py_S_G3yz_S_M2_vrr+WPY*I_ERI_Py_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M3_vrr;
      Double I_ERI_D2z_S_G3yz_S_M2_vrr = PAZ*I_ERI_Pz_S_G3yz_S_M2_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr+oned2k*I_ERI_Pz_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M2_vrr = PAX*I_ERI_Px_S_G2y2z_S_M2_vrr+WPX*I_ERI_Px_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M2_vrr = PAY*I_ERI_Py_S_G2y2z_S_M2_vrr+WPY*I_ERI_Py_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M2_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_M2_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M3_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M2_vrr = PAX*I_ERI_Px_S_Gy3z_S_M2_vrr+WPX*I_ERI_Px_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M2_vrr = PAY*I_ERI_Py_S_Gy3z_S_M2_vrr+WPY*I_ERI_Py_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr+oned2k*I_ERI_Py_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2x_S_G4z_S_M2_vrr = PAX*I_ERI_Px_S_G4z_S_M2_vrr+WPX*I_ERI_Px_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_D2y_S_G4z_S_M2_vrr = PAY*I_ERI_Py_S_G4z_S_M2_vrr+WPY*I_ERI_Py_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_D2z_S_G4z_S_M2_vrr = PAZ*I_ERI_Pz_S_G4z_S_M2_vrr+WPZ*I_ERI_Pz_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_Px_S_M1_vrr = QCX*I_ERI_S_S_S_S_M1_vrr+WQX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Py_S_M1_vrr = QCY*I_ERI_S_S_S_S_M1_vrr+WQY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Pz_S_M1_vrr = QCZ*I_ERI_S_S_S_S_M1_vrr+WQZ*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M1_vrr = QCX*I_ERI_S_S_Px_S_M1_vrr+WQX*I_ERI_S_S_Px_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dxy_S_M1_vrr = QCY*I_ERI_S_S_Px_S_M1_vrr+WQY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_D2y_S_M1_vrr = QCY*I_ERI_S_S_Py_S_M1_vrr+WQY*I_ERI_S_S_Py_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_D2z_S_M1_vrr = QCZ*I_ERI_S_S_Pz_S_M1_vrr+WQZ*I_ERI_S_S_Pz_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M1_vrr = QCX*I_ERI_S_S_D2x_S_M1_vrr+WQX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2e*I_ERI_S_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_F2xy_S_M1_vrr = QCY*I_ERI_S_S_D2x_S_M1_vrr+WQY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_F2xz_S_M1_vrr = QCZ*I_ERI_S_S_D2x_S_M1_vrr+WQZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_Fx2y_S_M1_vrr = QCX*I_ERI_S_S_D2y_S_M1_vrr+WQX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fxyz_S_M1_vrr = QCZ*I_ERI_S_S_Dxy_S_M1_vrr+WQZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_S_S_Fx2z_S_M1_vrr = QCX*I_ERI_S_S_D2z_S_M1_vrr+WQX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3y_S_M1_vrr = QCY*I_ERI_S_S_D2y_S_M1_vrr+WQY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2e*I_ERI_S_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_S_F2yz_S_M1_vrr = QCZ*I_ERI_S_S_D2y_S_M1_vrr+WQZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fy2z_S_M1_vrr = QCY*I_ERI_S_S_D2z_S_M1_vrr+WQY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3z_S_M1_vrr = QCZ*I_ERI_S_S_D2z_S_M1_vrr+WQZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2e*I_ERI_S_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M1_vrr = PAX*I_ERI_S_S_F3x_S_M1_vrr+WPX*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Py_S_F3x_S_M1_vrr = PAY*I_ERI_S_S_F3x_S_M1_vrr+WPY*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Pz_S_F3x_S_M1_vrr = PAZ*I_ERI_S_S_F3x_S_M1_vrr+WPZ*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Px_S_F2xy_S_M1_vrr = PAX*I_ERI_S_S_F2xy_S_M1_vrr+WPX*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Py_S_F2xy_S_M1_vrr = PAY*I_ERI_S_S_F2xy_S_M1_vrr+WPY*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Pz_S_F2xy_S_M1_vrr = PAZ*I_ERI_S_S_F2xy_S_M1_vrr+WPZ*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Px_S_F2xz_S_M1_vrr = PAX*I_ERI_S_S_F2xz_S_M1_vrr+WPX*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Py_S_F2xz_S_M1_vrr = PAY*I_ERI_S_S_F2xz_S_M1_vrr+WPY*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Pz_S_F2xz_S_M1_vrr = PAZ*I_ERI_S_S_F2xz_S_M1_vrr+WPZ*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Px_S_Fx2y_S_M1_vrr = PAX*I_ERI_S_S_Fx2y_S_M1_vrr+WPX*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Py_S_Fx2y_S_M1_vrr = PAY*I_ERI_S_S_Fx2y_S_M1_vrr+WPY*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_S_S_Fx2y_S_M1_vrr+WPZ*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Px_S_Fxyz_S_M1_vrr = PAX*I_ERI_S_S_Fxyz_S_M1_vrr+WPX*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Py_S_Fxyz_S_M1_vrr = PAY*I_ERI_S_S_Fxyz_S_M1_vrr+WPY*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_S_S_Fxyz_S_M1_vrr+WPZ*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Px_S_Fx2z_S_M1_vrr = PAX*I_ERI_S_S_Fx2z_S_M1_vrr+WPX*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Py_S_Fx2z_S_M1_vrr = PAY*I_ERI_S_S_Fx2z_S_M1_vrr+WPY*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_S_S_Fx2z_S_M1_vrr+WPZ*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Px_S_F3y_S_M1_vrr = PAX*I_ERI_S_S_F3y_S_M1_vrr+WPX*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Py_S_F3y_S_M1_vrr = PAY*I_ERI_S_S_F3y_S_M1_vrr+WPY*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Pz_S_F3y_S_M1_vrr = PAZ*I_ERI_S_S_F3y_S_M1_vrr+WPZ*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Px_S_F2yz_S_M1_vrr = PAX*I_ERI_S_S_F2yz_S_M1_vrr+WPX*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Py_S_F2yz_S_M1_vrr = PAY*I_ERI_S_S_F2yz_S_M1_vrr+WPY*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Pz_S_F2yz_S_M1_vrr = PAZ*I_ERI_S_S_F2yz_S_M1_vrr+WPZ*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Px_S_Fy2z_S_M1_vrr = PAX*I_ERI_S_S_Fy2z_S_M1_vrr+WPX*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Py_S_Fy2z_S_M1_vrr = PAY*I_ERI_S_S_Fy2z_S_M1_vrr+WPY*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_S_S_Fy2z_S_M1_vrr+WPZ*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Px_S_F3z_S_M1_vrr = PAX*I_ERI_S_S_F3z_S_M1_vrr+WPX*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Py_S_F3z_S_M1_vrr = PAY*I_ERI_S_S_F3z_S_M1_vrr+WPY*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Pz_S_F3z_S_M1_vrr = PAZ*I_ERI_S_S_F3z_S_M1_vrr+WPZ*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_S_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M1_vrr = QCX*I_ERI_S_S_F3x_S_M1_vrr+WQX*I_ERI_S_S_F3x_S_M2_vrr+3*oned2e*I_ERI_S_S_D2x_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_G3xy_S_M1_vrr = QCY*I_ERI_S_S_F3x_S_M1_vrr+WQY*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_S_S_G3xz_S_M1_vrr = QCZ*I_ERI_S_S_F3x_S_M1_vrr+WQZ*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_S_S_G2x2y_S_M1_vrr = QCY*I_ERI_S_S_F2xy_S_M1_vrr+WQY*I_ERI_S_S_F2xy_S_M2_vrr+oned2e*I_ERI_S_S_D2x_S_M1_vrr-rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_G2xyz_S_M1_vrr = QCZ*I_ERI_S_S_F2xy_S_M1_vrr+WQZ*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_S_S_G2x2z_S_M1_vrr = QCZ*I_ERI_S_S_F2xz_S_M1_vrr+WQZ*I_ERI_S_S_F2xz_S_M2_vrr+oned2e*I_ERI_S_S_D2x_S_M1_vrr-rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_Gx3y_S_M1_vrr = QCX*I_ERI_S_S_F3y_S_M1_vrr+WQX*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_S_S_Gx2yz_S_M1_vrr = QCZ*I_ERI_S_S_Fx2y_S_M1_vrr+WQZ*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_S_S_Gxy2z_S_M1_vrr = QCY*I_ERI_S_S_Fx2z_S_M1_vrr+WQY*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_S_S_Gx3z_S_M1_vrr = QCX*I_ERI_S_S_F3z_S_M1_vrr+WQX*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_S_S_G4y_S_M1_vrr = QCY*I_ERI_S_S_F3y_S_M1_vrr+WQY*I_ERI_S_S_F3y_S_M2_vrr+3*oned2e*I_ERI_S_S_D2y_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_G3yz_S_M1_vrr = QCZ*I_ERI_S_S_F3y_S_M1_vrr+WQZ*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_S_S_G2y2z_S_M1_vrr = QCZ*I_ERI_S_S_F2yz_S_M1_vrr+WQZ*I_ERI_S_S_F2yz_S_M2_vrr+oned2e*I_ERI_S_S_D2y_S_M1_vrr-rhod2esq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Gy3z_S_M1_vrr = QCY*I_ERI_S_S_F3z_S_M1_vrr+WQY*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_S_S_G4z_S_M1_vrr = QCZ*I_ERI_S_S_F3z_S_M1_vrr+WQZ*I_ERI_S_S_F3z_S_M2_vrr+3*oned2e*I_ERI_S_S_D2z_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M1_vrr = PAX*I_ERI_Px_S_F3x_S_M1_vrr+WPX*I_ERI_Px_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F3x_S_M1_vrr = PAY*I_ERI_Py_S_F3x_S_M1_vrr+WPY*I_ERI_Py_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_F3x_S_M1_vrr = PAZ*I_ERI_Pz_S_F3x_S_M1_vrr+WPZ*I_ERI_Pz_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_F2xy_S_M1_vrr = PAX*I_ERI_Px_S_F2xy_S_M1_vrr+WPX*I_ERI_Px_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_D2y_S_F2xy_S_M1_vrr = PAY*I_ERI_Py_S_F2xy_S_M1_vrr+WPY*I_ERI_Py_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_F2xz_S_M1_vrr = PAX*I_ERI_Px_S_F2xz_S_M1_vrr+WPX*I_ERI_Px_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_F2xz_S_M1_vrr = PAY*I_ERI_Py_S_F2xz_S_M1_vrr+WPY*I_ERI_Py_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M1_vrr = PAX*I_ERI_Px_S_Fx2y_S_M1_vrr+WPX*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_Py_S_Fx2y_S_M1_vrr+WPY*I_ERI_Py_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M1_vrr = PAX*I_ERI_Px_S_Fxyz_S_M1_vrr+WPX*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_Py_S_Fxyz_S_M1_vrr+WPY*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M1_vrr = PAX*I_ERI_Px_S_Fx2z_S_M1_vrr+WPX*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_Py_S_Fx2z_S_M1_vrr+WPY*I_ERI_Py_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M2_vrr;
      Double I_ERI_D2x_S_F3y_S_M1_vrr = PAX*I_ERI_Px_S_F3y_S_M1_vrr+WPX*I_ERI_Px_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_F3y_S_M1_vrr = PAY*I_ERI_Py_S_F3y_S_M1_vrr+WPY*I_ERI_Py_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F3y_S_M1_vrr = PAZ*I_ERI_Pz_S_F3y_S_M1_vrr+WPZ*I_ERI_Pz_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_F2yz_S_M1_vrr = PAX*I_ERI_Px_S_F2yz_S_M1_vrr+WPX*I_ERI_Px_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_D2y_S_F2yz_S_M1_vrr = PAY*I_ERI_Py_S_F2yz_S_M1_vrr+WPY*I_ERI_Py_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M2_vrr;
      Double I_ERI_D2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M1_vrr = PAX*I_ERI_Px_S_Fy2z_S_M1_vrr+WPX*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_Py_S_Fy2z_S_M1_vrr+WPY*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M2_vrr;
      Double I_ERI_D2x_S_F3z_S_M1_vrr = PAX*I_ERI_Px_S_F3z_S_M1_vrr+WPX*I_ERI_Px_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_F3z_S_M1_vrr = PAY*I_ERI_Py_S_F3z_S_M1_vrr+WPY*I_ERI_Py_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_F3z_S_M1_vrr = PAZ*I_ERI_Pz_S_F3z_S_M1_vrr+WPZ*I_ERI_Pz_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M1_vrr = PAX*I_ERI_S_S_G4x_S_M1_vrr+WPX*I_ERI_S_S_G4x_S_M2_vrr+4*oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Py_S_G4x_S_M1_vrr = PAY*I_ERI_S_S_G4x_S_M1_vrr+WPY*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_Pz_S_G4x_S_M1_vrr = PAZ*I_ERI_S_S_G4x_S_M1_vrr+WPZ*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_Px_S_G3xy_S_M1_vrr = PAX*I_ERI_S_S_G3xy_S_M1_vrr+WPX*I_ERI_S_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Py_S_G3xy_S_M1_vrr = PAY*I_ERI_S_S_G3xy_S_M1_vrr+WPY*I_ERI_S_S_G3xy_S_M2_vrr+oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Pz_S_G3xy_S_M1_vrr = PAZ*I_ERI_S_S_G3xy_S_M1_vrr+WPZ*I_ERI_S_S_G3xy_S_M2_vrr;
      Double I_ERI_Px_S_G3xz_S_M1_vrr = PAX*I_ERI_S_S_G3xz_S_M1_vrr+WPX*I_ERI_S_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Py_S_G3xz_S_M1_vrr = PAY*I_ERI_S_S_G3xz_S_M1_vrr+WPY*I_ERI_S_S_G3xz_S_M2_vrr;
      Double I_ERI_Pz_S_G3xz_S_M1_vrr = PAZ*I_ERI_S_S_G3xz_S_M1_vrr+WPZ*I_ERI_S_S_G3xz_S_M2_vrr+oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Px_S_G2x2y_S_M1_vrr = PAX*I_ERI_S_S_G2x2y_S_M1_vrr+WPX*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Py_S_G2x2y_S_M1_vrr = PAY*I_ERI_S_S_G2x2y_S_M1_vrr+WPY*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_S_S_G2x2y_S_M1_vrr+WPZ*I_ERI_S_S_G2x2y_S_M2_vrr;
      Double I_ERI_Px_S_G2xyz_S_M1_vrr = PAX*I_ERI_S_S_G2xyz_S_M1_vrr+WPX*I_ERI_S_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Py_S_G2xyz_S_M1_vrr = PAY*I_ERI_S_S_G2xyz_S_M1_vrr+WPY*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_S_S_G2xyz_S_M1_vrr+WPZ*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Px_S_G2x2z_S_M1_vrr = PAX*I_ERI_S_S_G2x2z_S_M1_vrr+WPX*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Py_S_G2x2z_S_M1_vrr = PAY*I_ERI_S_S_G2x2z_S_M1_vrr+WPY*I_ERI_S_S_G2x2z_S_M2_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_S_S_G2x2z_S_M1_vrr+WPZ*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Px_S_Gx3y_S_M1_vrr = PAX*I_ERI_S_S_Gx3y_S_M1_vrr+WPX*I_ERI_S_S_Gx3y_S_M2_vrr+oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Py_S_Gx3y_S_M1_vrr = PAY*I_ERI_S_S_Gx3y_S_M1_vrr+WPY*I_ERI_S_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_S_S_Gx3y_S_M1_vrr+WPZ*I_ERI_S_S_Gx3y_S_M2_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M1_vrr = PAX*I_ERI_S_S_Gx2yz_S_M1_vrr+WPX*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M1_vrr = PAY*I_ERI_S_S_Gx2yz_S_M1_vrr+WPY*I_ERI_S_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M1_vrr = PAX*I_ERI_S_S_Gxy2z_S_M1_vrr+WPX*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M1_vrr = PAY*I_ERI_S_S_Gxy2z_S_M1_vrr+WPY*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Px_S_Gx3z_S_M1_vrr = PAX*I_ERI_S_S_Gx3z_S_M1_vrr+WPX*I_ERI_S_S_Gx3z_S_M2_vrr+oned2k*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Py_S_Gx3z_S_M1_vrr = PAY*I_ERI_S_S_Gx3z_S_M1_vrr+WPY*I_ERI_S_S_Gx3z_S_M2_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_S_S_Gx3z_S_M1_vrr+WPZ*I_ERI_S_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Px_S_G4y_S_M1_vrr = PAX*I_ERI_S_S_G4y_S_M1_vrr+WPX*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_Py_S_G4y_S_M1_vrr = PAY*I_ERI_S_S_G4y_S_M1_vrr+WPY*I_ERI_S_S_G4y_S_M2_vrr+4*oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Pz_S_G4y_S_M1_vrr = PAZ*I_ERI_S_S_G4y_S_M1_vrr+WPZ*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_Px_S_G3yz_S_M1_vrr = PAX*I_ERI_S_S_G3yz_S_M1_vrr+WPX*I_ERI_S_S_G3yz_S_M2_vrr;
      Double I_ERI_Py_S_G3yz_S_M1_vrr = PAY*I_ERI_S_S_G3yz_S_M1_vrr+WPY*I_ERI_S_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Pz_S_G3yz_S_M1_vrr = PAZ*I_ERI_S_S_G3yz_S_M1_vrr+WPZ*I_ERI_S_S_G3yz_S_M2_vrr+oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Px_S_G2y2z_S_M1_vrr = PAX*I_ERI_S_S_G2y2z_S_M1_vrr+WPX*I_ERI_S_S_G2y2z_S_M2_vrr;
      Double I_ERI_Py_S_G2y2z_S_M1_vrr = PAY*I_ERI_S_S_G2y2z_S_M1_vrr+WPY*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_S_S_G2y2z_S_M1_vrr+WPZ*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Px_S_Gy3z_S_M1_vrr = PAX*I_ERI_S_S_Gy3z_S_M1_vrr+WPX*I_ERI_S_S_Gy3z_S_M2_vrr;
      Double I_ERI_Py_S_Gy3z_S_M1_vrr = PAY*I_ERI_S_S_Gy3z_S_M1_vrr+WPY*I_ERI_S_S_Gy3z_S_M2_vrr+oned2k*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_S_S_Gy3z_S_M1_vrr+WPZ*I_ERI_S_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Px_S_G4z_S_M1_vrr = PAX*I_ERI_S_S_G4z_S_M1_vrr+WPX*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_Py_S_G4z_S_M1_vrr = PAY*I_ERI_S_S_G4z_S_M1_vrr+WPY*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_Pz_S_G4z_S_M1_vrr = PAZ*I_ERI_S_S_G4z_S_M1_vrr+WPZ*I_ERI_S_S_G4z_S_M2_vrr+4*oned2k*I_ERI_S_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 20 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M1_vrr = PAX*I_ERI_D2x_S_F3x_S_M1_vrr+WPX*I_ERI_D2x_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F3x_S_M1_vrr = PAY*I_ERI_D2x_S_F3x_S_M1_vrr+WPY*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2x_S_F3x_S_M1_vrr+WPZ*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3x_S_M1_vrr = PAX*I_ERI_D2y_S_F3x_S_M1_vrr+WPX*I_ERI_D2y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3x_S_M1_vrr = PAX*I_ERI_D2z_S_F3x_S_M1_vrr+WPX*I_ERI_D2z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F3x_S_M1_vrr = PAY*I_ERI_D2y_S_F3x_S_M1_vrr+WPY*I_ERI_D2y_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_F2yz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2y_S_F3x_S_M1_vrr+WPZ*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F3z_S_F3x_S_M1_vrr = PAZ*I_ERI_D2z_S_F3x_S_M1_vrr+WPZ*I_ERI_D2z_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_F3x_S_F2xy_S_M1_vrr = PAX*I_ERI_D2x_S_F2xy_S_M1_vrr+WPX*I_ERI_D2x_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M1_vrr = PAY*I_ERI_D2x_S_F2xy_S_M1_vrr+WPY*I_ERI_D2x_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_M1_vrr = PAX*I_ERI_D2y_S_F2xy_S_M1_vrr+WPX*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_M1_vrr = PAX*I_ERI_D2z_S_F2xy_S_M1_vrr+WPX*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3y_S_F2xy_S_M1_vrr = PAY*I_ERI_D2y_S_F2xy_S_M1_vrr+WPY*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F3z_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_F3x_S_F2xz_S_M1_vrr = PAX*I_ERI_D2x_S_F2xz_S_M1_vrr+WPX*I_ERI_D2x_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M1_vrr = PAY*I_ERI_D2x_S_F2xz_S_M1_vrr+WPY*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_M1_vrr = PAX*I_ERI_D2y_S_F2xz_S_M1_vrr+WPX*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_M1_vrr = PAX*I_ERI_D2z_S_F2xz_S_M1_vrr+WPX*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3y_S_F2xz_S_M1_vrr = PAY*I_ERI_D2y_S_F2xz_S_M1_vrr+WPY*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3y_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3x_S_F3y_S_M1_vrr = PAX*I_ERI_D2x_S_F3y_S_M1_vrr+WPX*I_ERI_D2x_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_F2xy_S_F3y_S_M1_vrr = PAY*I_ERI_D2x_S_F3y_S_M1_vrr+WPY*I_ERI_D2x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2x_S_F3y_S_M1_vrr+WPZ*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3y_S_M1_vrr = PAX*I_ERI_D2y_S_F3y_S_M1_vrr+WPX*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3y_S_M1_vrr = PAX*I_ERI_D2z_S_F3y_S_M1_vrr+WPX*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3y_S_F3y_S_M1_vrr = PAY*I_ERI_D2y_S_F3y_S_M1_vrr+WPY*I_ERI_D2y_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2y_S_F3y_S_M1_vrr+WPZ*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F3z_S_F3y_S_M1_vrr = PAZ*I_ERI_D2z_S_F3y_S_M1_vrr+WPZ*I_ERI_D2z_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_F3x_S_F2yz_S_M1_vrr = PAX*I_ERI_D2x_S_F2yz_S_M1_vrr+WPX*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M1_vrr = PAY*I_ERI_D2x_S_F2yz_S_M1_vrr+WPY*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_M1_vrr = PAX*I_ERI_D2y_S_F2yz_S_M1_vrr+WPX*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_M1_vrr = PAX*I_ERI_D2z_S_F2yz_S_M1_vrr+WPX*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3y_S_F2yz_S_M1_vrr = PAY*I_ERI_D2y_S_F2yz_S_M1_vrr+WPY*I_ERI_D2y_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F3z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3x_S_F3z_S_M1_vrr = PAX*I_ERI_D2x_S_F3z_S_M1_vrr+WPX*I_ERI_D2x_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_F2xy_S_F3z_S_M1_vrr = PAY*I_ERI_D2x_S_F3z_S_M1_vrr+WPY*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2x_S_F3z_S_M1_vrr+WPZ*I_ERI_D2x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3z_S_M1_vrr = PAX*I_ERI_D2y_S_F3z_S_M1_vrr+WPX*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3z_S_M1_vrr = PAX*I_ERI_D2z_S_F3z_S_M1_vrr+WPX*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3y_S_F3z_S_M1_vrr = PAY*I_ERI_D2y_S_F3z_S_M1_vrr+WPY*I_ERI_D2y_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_F2yz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2y_S_F3z_S_M1_vrr+WPZ*I_ERI_D2y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_F3z_S_F3z_S_M1_vrr = PAZ*I_ERI_D2z_S_F3z_S_M1_vrr+WPZ*I_ERI_D2z_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M1_vrr = PAX*I_ERI_Px_S_G4x_S_M1_vrr+WPX*I_ERI_Px_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_D2y_S_G4x_S_M1_vrr = PAY*I_ERI_Py_S_G4x_S_M1_vrr+WPY*I_ERI_Py_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2z_S_G4x_S_M1_vrr = PAZ*I_ERI_Pz_S_G4x_S_M1_vrr+WPZ*I_ERI_Pz_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2x_S_G3xy_S_M1_vrr = PAX*I_ERI_Px_S_G3xy_S_M1_vrr+WPX*I_ERI_Px_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M2_vrr;
      Double I_ERI_D2y_S_G3xy_S_M1_vrr = PAY*I_ERI_Py_S_G3xy_S_M1_vrr+WPY*I_ERI_Py_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+oned2k*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xy_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr;
      Double I_ERI_D2x_S_G3xz_S_M1_vrr = PAX*I_ERI_Px_S_G3xz_S_M1_vrr+WPX*I_ERI_Px_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_D2y_S_G3xz_S_M1_vrr = PAY*I_ERI_Py_S_G3xz_S_M1_vrr+WPY*I_ERI_Py_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr;
      Double I_ERI_D2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M1_vrr = PAX*I_ERI_Px_S_G2x2y_S_M1_vrr+WPX*I_ERI_Px_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M1_vrr = PAY*I_ERI_Py_S_G2x2y_S_M1_vrr+WPY*I_ERI_Py_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M1_vrr = PAX*I_ERI_Px_S_G2xyz_S_M1_vrr+WPX*I_ERI_Px_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M1_vrr = PAY*I_ERI_Py_S_G2xyz_S_M1_vrr+WPY*I_ERI_Py_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M1_vrr = PAX*I_ERI_Px_S_G2x2z_S_M1_vrr+WPX*I_ERI_Px_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M1_vrr = PAY*I_ERI_Py_S_G2x2z_S_M1_vrr+WPY*I_ERI_Py_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M1_vrr = PAX*I_ERI_Px_S_Gx3y_S_M1_vrr+WPX*I_ERI_Px_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+oned2k*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M1_vrr = PAY*I_ERI_Py_S_Gx3y_S_M1_vrr+WPY*I_ERI_Py_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_Px_S_Gx2yz_S_M1_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_Py_S_Gx2yz_S_M1_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_Px_S_Gxy2z_S_M1_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_Py_S_Gxy2z_S_M1_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M1_vrr = PAX*I_ERI_Px_S_Gx3z_S_M1_vrr+WPX*I_ERI_Px_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+oned2k*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M1_vrr = PAY*I_ERI_Py_S_Gx3z_S_M1_vrr+WPY*I_ERI_Py_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4y_S_M1_vrr = PAX*I_ERI_Px_S_G4y_S_M1_vrr+WPX*I_ERI_Px_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_D2y_S_G4y_S_M1_vrr = PAY*I_ERI_Py_S_G4y_S_M1_vrr+WPY*I_ERI_Py_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M2_vrr;
      Double I_ERI_D2z_S_G4y_S_M1_vrr = PAZ*I_ERI_Pz_S_G4y_S_M1_vrr+WPZ*I_ERI_Pz_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_D2x_S_G3yz_S_M1_vrr = PAX*I_ERI_Px_S_G3yz_S_M1_vrr+WPX*I_ERI_Px_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr;
      Double I_ERI_D2y_S_G3yz_S_M1_vrr = PAY*I_ERI_Py_S_G3yz_S_M1_vrr+WPY*I_ERI_Py_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M2_vrr;
      Double I_ERI_D2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3yz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M1_vrr = PAX*I_ERI_Px_S_G2y2z_S_M1_vrr+WPX*I_ERI_Px_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M1_vrr = PAY*I_ERI_Py_S_G2y2z_S_M1_vrr+WPY*I_ERI_Py_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M2_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M1_vrr = PAX*I_ERI_Px_S_Gy3z_S_M1_vrr+WPX*I_ERI_Px_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M1_vrr = PAY*I_ERI_Py_S_Gy3z_S_M1_vrr+WPY*I_ERI_Py_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+oned2k*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4z_S_M1_vrr = PAX*I_ERI_Px_S_G4z_S_M1_vrr+WPX*I_ERI_Px_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_D2y_S_G4z_S_M1_vrr = PAY*I_ERI_Py_S_G4z_S_M1_vrr+WPY*I_ERI_Py_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_D2z_S_G4z_S_M1_vrr = PAZ*I_ERI_Pz_S_G4z_S_M1_vrr+WPZ*I_ERI_Pz_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_M1_vrr = PAX*I_ERI_D2x_S_G4x_S_M1_vrr+WPX*I_ERI_D2x_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xy_S_G4x_S_M1_vrr = PAY*I_ERI_D2x_S_G4x_S_M1_vrr+WPY*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_F2xz_S_G4x_S_M1_vrr = PAZ*I_ERI_D2x_S_G4x_S_M1_vrr+WPZ*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4x_S_M1_vrr = PAX*I_ERI_D2y_S_G4x_S_M1_vrr+WPX*I_ERI_D2y_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4x_S_M1_vrr = PAX*I_ERI_D2z_S_G4x_S_M1_vrr+WPX*I_ERI_D2z_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_F3y_S_G4x_S_M1_vrr = PAY*I_ERI_D2y_S_G4x_S_M1_vrr+WPY*I_ERI_D2y_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4x_S_M2_vrr;
      Double I_ERI_F2yz_S_G4x_S_M1_vrr = PAZ*I_ERI_D2y_S_G4x_S_M1_vrr+WPZ*I_ERI_D2y_S_G4x_S_M2_vrr;
      Double I_ERI_F3z_S_G4x_S_M1_vrr = PAZ*I_ERI_D2z_S_G4x_S_M1_vrr+WPZ*I_ERI_D2z_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4x_S_M2_vrr;
      Double I_ERI_F3x_S_G3xy_S_M1_vrr = PAX*I_ERI_D2x_S_G3xy_S_M1_vrr+WPX*I_ERI_D2x_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_F2xy_S_G3xy_S_M1_vrr = PAY*I_ERI_D2x_S_G3xy_S_M1_vrr+WPY*I_ERI_D2x_S_G3xy_S_M2_vrr+oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xz_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2x_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2x_S_G3xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_M1_vrr = PAX*I_ERI_D2y_S_G3xy_S_M1_vrr+WPX*I_ERI_D2y_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_M1_vrr = PAX*I_ERI_D2z_S_G3xy_S_M1_vrr+WPX*I_ERI_D2z_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M2_vrr;
      Double I_ERI_F3y_S_G3xy_S_M1_vrr = PAY*I_ERI_D2y_S_G3xy_S_M1_vrr+WPY*I_ERI_D2y_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3xy_S_M2_vrr+oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F2yz_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2y_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2y_S_G3xy_S_M2_vrr;
      Double I_ERI_F3z_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2z_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2z_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xy_S_M2_vrr;
      Double I_ERI_F3x_S_G3xz_S_M1_vrr = PAX*I_ERI_D2x_S_G3xz_S_M1_vrr+WPX*I_ERI_D2x_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xy_S_G3xz_S_M1_vrr = PAY*I_ERI_D2x_S_G3xz_S_M1_vrr+WPY*I_ERI_D2x_S_G3xz_S_M2_vrr;
      Double I_ERI_F2xz_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2x_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2x_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_M1_vrr = PAX*I_ERI_D2y_S_G3xz_S_M1_vrr+WPX*I_ERI_D2y_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_M1_vrr = PAX*I_ERI_D2z_S_G3xz_S_M1_vrr+WPX*I_ERI_D2z_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M2_vrr;
      Double I_ERI_F3y_S_G3xz_S_M1_vrr = PAY*I_ERI_D2y_S_G3xz_S_M1_vrr+WPY*I_ERI_D2y_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3xz_S_M2_vrr;
      Double I_ERI_F2yz_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2y_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2y_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F3z_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2z_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2z_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_F3x_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_M1_vrr = PAY*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPY*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2x_S_G2x2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2z_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3y_S_G2x2y_S_M1_vrr = PAY*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPY*I_ERI_D2y_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2y_S_G2x2y_S_M2_vrr;
      Double I_ERI_F3z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2z_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2z_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2y_S_M2_vrr;
      Double I_ERI_F3x_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2x_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_M1_vrr = PAY*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPY*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2y_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2z_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2z_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3y_S_G2xyz_S_M1_vrr = PAY*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPY*I_ERI_D2y_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F3z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2z_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2z_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M2_vrr;
      Double I_ERI_F3x_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_M1_vrr = PAY*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPY*I_ERI_D2x_S_G2x2z_S_M2_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2z_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3y_S_G2x2z_S_M1_vrr = PAY*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPY*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2z_S_M2_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_F3z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2z_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2z_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_F2xz_S_M2_vrr;
      Double I_ERI_F3x_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2x_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_M1_vrr = PAY*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPY*I_ERI_D2x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2y_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2z_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2z_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3y_S_Gx3y_S_M1_vrr = PAY*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPY*I_ERI_D2y_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx3y_S_M2_vrr;
      Double I_ERI_F3z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3y_S_M2_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2x_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_M1_vrr = PAY*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPY*I_ERI_D2x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2z_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPY*I_ERI_D2y_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2x_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_M1_vrr = PAY*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPY*I_ERI_D2x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2z_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPY*I_ERI_D2y_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3x_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2x_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_M1_vrr = PAY*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPY*I_ERI_D2x_S_Gx3z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2y_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2z_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2z_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3y_S_Gx3z_S_M1_vrr = PAY*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPY*I_ERI_D2y_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3x_S_G4y_S_M1_vrr = PAX*I_ERI_D2x_S_G4y_S_M1_vrr+WPX*I_ERI_D2x_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4y_S_M2_vrr;
      Double I_ERI_F2xy_S_G4y_S_M1_vrr = PAY*I_ERI_D2x_S_G4y_S_M1_vrr+WPY*I_ERI_D2x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_F2xz_S_G4y_S_M1_vrr = PAZ*I_ERI_D2x_S_G4y_S_M1_vrr+WPZ*I_ERI_D2x_S_G4y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4y_S_M1_vrr = PAX*I_ERI_D2y_S_G4y_S_M1_vrr+WPX*I_ERI_D2y_S_G4y_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4y_S_M1_vrr = PAX*I_ERI_D2z_S_G4y_S_M1_vrr+WPX*I_ERI_D2z_S_G4y_S_M2_vrr;
      Double I_ERI_F3y_S_G4y_S_M1_vrr = PAY*I_ERI_D2y_S_G4y_S_M1_vrr+WPY*I_ERI_D2y_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4y_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F2yz_S_G4y_S_M1_vrr = PAZ*I_ERI_D2y_S_G4y_S_M1_vrr+WPZ*I_ERI_D2y_S_G4y_S_M2_vrr;
      Double I_ERI_F3z_S_G4y_S_M1_vrr = PAZ*I_ERI_D2z_S_G4y_S_M1_vrr+WPZ*I_ERI_D2z_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4y_S_M2_vrr;
      Double I_ERI_F3x_S_G3yz_S_M1_vrr = PAX*I_ERI_D2x_S_G3yz_S_M1_vrr+WPX*I_ERI_D2x_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3yz_S_M2_vrr;
      Double I_ERI_F2xy_S_G3yz_S_M1_vrr = PAY*I_ERI_D2x_S_G3yz_S_M1_vrr+WPY*I_ERI_D2x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xz_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2x_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2x_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_M1_vrr = PAX*I_ERI_D2y_S_G3yz_S_M1_vrr+WPX*I_ERI_D2y_S_G3yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_M1_vrr = PAX*I_ERI_D2z_S_G3yz_S_M1_vrr+WPX*I_ERI_D2z_S_G3yz_S_M2_vrr;
      Double I_ERI_F3y_S_G3yz_S_M1_vrr = PAY*I_ERI_D2y_S_G3yz_S_M1_vrr+WPY*I_ERI_D2y_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_F2yz_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2y_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2y_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F3z_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2z_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2z_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3x_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2y2z_S_M2_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_M1_vrr = PAY*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPY*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2y_S_G2y2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2z_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2z_S_G2y2z_S_M2_vrr;
      Double I_ERI_F3y_S_G2y2z_S_M1_vrr = PAY*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPY*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_F3z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2z_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2z_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3x_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2x_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gy3z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_M1_vrr = PAY*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPY*I_ERI_D2x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2y_S_Gy3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2z_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2z_S_Gy3z_S_M2_vrr;
      Double I_ERI_F3y_S_Gy3z_S_M1_vrr = PAY*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPY*I_ERI_D2y_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gy3z_S_M2_vrr+oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3x_S_G4z_S_M1_vrr = PAX*I_ERI_D2x_S_G4z_S_M1_vrr+WPX*I_ERI_D2x_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4z_S_M2_vrr;
      Double I_ERI_F2xy_S_G4z_S_M1_vrr = PAY*I_ERI_D2x_S_G4z_S_M1_vrr+WPY*I_ERI_D2x_S_G4z_S_M2_vrr;
      Double I_ERI_F2xz_S_G4z_S_M1_vrr = PAZ*I_ERI_D2x_S_G4z_S_M1_vrr+WPZ*I_ERI_D2x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4z_S_M1_vrr = PAX*I_ERI_D2y_S_G4z_S_M1_vrr+WPX*I_ERI_D2y_S_G4z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4z_S_M1_vrr = PAX*I_ERI_D2z_S_G4z_S_M1_vrr+WPX*I_ERI_D2z_S_G4z_S_M2_vrr;
      Double I_ERI_F3y_S_G4z_S_M1_vrr = PAY*I_ERI_D2y_S_G4z_S_M1_vrr+WPY*I_ERI_D2y_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4z_S_M2_vrr;
      Double I_ERI_F2yz_S_G4z_S_M1_vrr = PAZ*I_ERI_D2y_S_G4z_S_M1_vrr+WPZ*I_ERI_D2y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_F3z_S_G4z_S_M1_vrr = PAZ*I_ERI_D2z_S_G4z_S_M1_vrr+WPZ*I_ERI_D2z_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2z_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_Px_S_vrr = QCX*I_ERI_S_S_S_S_vrr+WQX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Py_S_vrr = QCY*I_ERI_S_S_S_S_vrr+WQY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Pz_S_vrr = QCZ*I_ERI_S_S_S_S_vrr+WQZ*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_S_F3x_S_vrr = QCX*I_ERI_S_S_D2x_S_vrr+WQX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2e*I_ERI_S_S_Px_S_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_F2xy_S_vrr = QCY*I_ERI_S_S_D2x_S_vrr+WQY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_F2xz_S_vrr = QCZ*I_ERI_S_S_D2x_S_vrr+WQZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Fx2y_S_vrr = QCX*I_ERI_S_S_D2y_S_vrr+WQX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fx2z_S_vrr = QCX*I_ERI_S_S_D2z_S_vrr+WQX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3y_S_vrr = QCY*I_ERI_S_S_D2y_S_vrr+WQY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2e*I_ERI_S_S_Py_S_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_F2yz_S_vrr = QCZ*I_ERI_S_S_D2y_S_vrr+WQZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_F3z_S_vrr = QCZ*I_ERI_S_S_D2z_S_vrr+WQZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2e*I_ERI_S_S_Pz_S_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_S_S_G4x_S_vrr = QCX*I_ERI_S_S_F3x_S_vrr+WQX*I_ERI_S_S_F3x_S_M1_vrr+3*oned2e*I_ERI_S_S_D2x_S_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_G3xy_S_vrr = QCY*I_ERI_S_S_F3x_S_vrr+WQY*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_S_S_G3xz_S_vrr = QCZ*I_ERI_S_S_F3x_S_vrr+WQZ*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_S_S_G2x2y_S_vrr = QCY*I_ERI_S_S_F2xy_S_vrr+WQY*I_ERI_S_S_F2xy_S_M1_vrr+oned2e*I_ERI_S_S_D2x_S_vrr-rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_G2xyz_S_vrr = QCZ*I_ERI_S_S_F2xy_S_vrr+WQZ*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_S_S_G2x2z_S_vrr = QCZ*I_ERI_S_S_F2xz_S_vrr+WQZ*I_ERI_S_S_F2xz_S_M1_vrr+oned2e*I_ERI_S_S_D2x_S_vrr-rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Gx3y_S_vrr = QCX*I_ERI_S_S_F3y_S_vrr+WQX*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_S_S_Gx2yz_S_vrr = QCZ*I_ERI_S_S_Fx2y_S_vrr+WQZ*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_S_S_Gxy2z_S_vrr = QCY*I_ERI_S_S_Fx2z_S_vrr+WQY*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_S_S_Gx3z_S_vrr = QCX*I_ERI_S_S_F3z_S_vrr+WQX*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_S_S_G4y_S_vrr = QCY*I_ERI_S_S_F3y_S_vrr+WQY*I_ERI_S_S_F3y_S_M1_vrr+3*oned2e*I_ERI_S_S_D2y_S_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_G3yz_S_vrr = QCZ*I_ERI_S_S_F3y_S_vrr+WQZ*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_S_S_G2y2z_S_vrr = QCZ*I_ERI_S_S_F2yz_S_vrr+WQZ*I_ERI_S_S_F2yz_S_M1_vrr+oned2e*I_ERI_S_S_D2y_S_vrr-rhod2esq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Gy3z_S_vrr = QCY*I_ERI_S_S_F3z_S_vrr+WQY*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_S_S_G4z_S_vrr = QCZ*I_ERI_S_S_F3z_S_vrr+WQZ*I_ERI_S_S_F3z_S_M1_vrr+3*oned2e*I_ERI_S_S_D2z_S_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_vrr = PAX*I_ERI_S_S_G4x_S_vrr+WPX*I_ERI_S_S_G4x_S_M1_vrr+4*oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Py_S_G4x_S_vrr = PAY*I_ERI_S_S_G4x_S_vrr+WPY*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_Pz_S_G4x_S_vrr = PAZ*I_ERI_S_S_G4x_S_vrr+WPZ*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_Px_S_G3xy_S_vrr = PAX*I_ERI_S_S_G3xy_S_vrr+WPX*I_ERI_S_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Py_S_G3xy_S_vrr = PAY*I_ERI_S_S_G3xy_S_vrr+WPY*I_ERI_S_S_G3xy_S_M1_vrr+oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Pz_S_G3xy_S_vrr = PAZ*I_ERI_S_S_G3xy_S_vrr+WPZ*I_ERI_S_S_G3xy_S_M1_vrr;
      Double I_ERI_Px_S_G3xz_S_vrr = PAX*I_ERI_S_S_G3xz_S_vrr+WPX*I_ERI_S_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Py_S_G3xz_S_vrr = PAY*I_ERI_S_S_G3xz_S_vrr+WPY*I_ERI_S_S_G3xz_S_M1_vrr;
      Double I_ERI_Pz_S_G3xz_S_vrr = PAZ*I_ERI_S_S_G3xz_S_vrr+WPZ*I_ERI_S_S_G3xz_S_M1_vrr+oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Px_S_G2x2y_S_vrr = PAX*I_ERI_S_S_G2x2y_S_vrr+WPX*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Py_S_G2x2y_S_vrr = PAY*I_ERI_S_S_G2x2y_S_vrr+WPY*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Pz_S_G2x2y_S_vrr = PAZ*I_ERI_S_S_G2x2y_S_vrr+WPZ*I_ERI_S_S_G2x2y_S_M1_vrr;
      Double I_ERI_Px_S_G2xyz_S_vrr = PAX*I_ERI_S_S_G2xyz_S_vrr+WPX*I_ERI_S_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Py_S_G2xyz_S_vrr = PAY*I_ERI_S_S_G2xyz_S_vrr+WPY*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Pz_S_G2xyz_S_vrr = PAZ*I_ERI_S_S_G2xyz_S_vrr+WPZ*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Px_S_G2x2z_S_vrr = PAX*I_ERI_S_S_G2x2z_S_vrr+WPX*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Py_S_G2x2z_S_vrr = PAY*I_ERI_S_S_G2x2z_S_vrr+WPY*I_ERI_S_S_G2x2z_S_M1_vrr;
      Double I_ERI_Pz_S_G2x2z_S_vrr = PAZ*I_ERI_S_S_G2x2z_S_vrr+WPZ*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Px_S_Gx3y_S_vrr = PAX*I_ERI_S_S_Gx3y_S_vrr+WPX*I_ERI_S_S_Gx3y_S_M1_vrr+oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Py_S_Gx3y_S_vrr = PAY*I_ERI_S_S_Gx3y_S_vrr+WPY*I_ERI_S_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Pz_S_Gx3y_S_vrr = PAZ*I_ERI_S_S_Gx3y_S_vrr+WPZ*I_ERI_S_S_Gx3y_S_M1_vrr;
      Double I_ERI_Px_S_Gx2yz_S_vrr = PAX*I_ERI_S_S_Gx2yz_S_vrr+WPX*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Py_S_Gx2yz_S_vrr = PAY*I_ERI_S_S_Gx2yz_S_vrr+WPY*I_ERI_S_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_vrr = PAZ*I_ERI_S_S_Gx2yz_S_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Px_S_Gxy2z_S_vrr = PAX*I_ERI_S_S_Gxy2z_S_vrr+WPX*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Py_S_Gxy2z_S_vrr = PAY*I_ERI_S_S_Gxy2z_S_vrr+WPY*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_vrr = PAZ*I_ERI_S_S_Gxy2z_S_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Px_S_Gx3z_S_vrr = PAX*I_ERI_S_S_Gx3z_S_vrr+WPX*I_ERI_S_S_Gx3z_S_M1_vrr+oned2k*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Py_S_Gx3z_S_vrr = PAY*I_ERI_S_S_Gx3z_S_vrr+WPY*I_ERI_S_S_Gx3z_S_M1_vrr;
      Double I_ERI_Pz_S_Gx3z_S_vrr = PAZ*I_ERI_S_S_Gx3z_S_vrr+WPZ*I_ERI_S_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Px_S_G4y_S_vrr = PAX*I_ERI_S_S_G4y_S_vrr+WPX*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_Py_S_G4y_S_vrr = PAY*I_ERI_S_S_G4y_S_vrr+WPY*I_ERI_S_S_G4y_S_M1_vrr+4*oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Pz_S_G4y_S_vrr = PAZ*I_ERI_S_S_G4y_S_vrr+WPZ*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_Px_S_G3yz_S_vrr = PAX*I_ERI_S_S_G3yz_S_vrr+WPX*I_ERI_S_S_G3yz_S_M1_vrr;
      Double I_ERI_Py_S_G3yz_S_vrr = PAY*I_ERI_S_S_G3yz_S_vrr+WPY*I_ERI_S_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Pz_S_G3yz_S_vrr = PAZ*I_ERI_S_S_G3yz_S_vrr+WPZ*I_ERI_S_S_G3yz_S_M1_vrr+oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Px_S_G2y2z_S_vrr = PAX*I_ERI_S_S_G2y2z_S_vrr+WPX*I_ERI_S_S_G2y2z_S_M1_vrr;
      Double I_ERI_Py_S_G2y2z_S_vrr = PAY*I_ERI_S_S_G2y2z_S_vrr+WPY*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Pz_S_G2y2z_S_vrr = PAZ*I_ERI_S_S_G2y2z_S_vrr+WPZ*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Px_S_Gy3z_S_vrr = PAX*I_ERI_S_S_Gy3z_S_vrr+WPX*I_ERI_S_S_Gy3z_S_M1_vrr;
      Double I_ERI_Py_S_Gy3z_S_vrr = PAY*I_ERI_S_S_Gy3z_S_vrr+WPY*I_ERI_S_S_Gy3z_S_M1_vrr+oned2k*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Pz_S_Gy3z_S_vrr = PAZ*I_ERI_S_S_Gy3z_S_vrr+WPZ*I_ERI_S_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Px_S_G4z_S_vrr = PAX*I_ERI_S_S_G4z_S_vrr+WPX*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Py_S_G4z_S_vrr = PAY*I_ERI_S_S_G4z_S_vrr+WPY*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Pz_S_G4z_S_vrr = PAZ*I_ERI_S_S_G4z_S_vrr+WPZ*I_ERI_S_S_G4z_S_M1_vrr+4*oned2k*I_ERI_S_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_vrr = PAX*I_ERI_Px_S_G4x_S_vrr+WPX*I_ERI_Px_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G4x_S_vrr = PAY*I_ERI_Py_S_G4x_S_vrr+WPY*I_ERI_Py_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_D2z_S_G4x_S_vrr = PAZ*I_ERI_Pz_S_G4x_S_vrr+WPZ*I_ERI_Pz_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_D2x_S_G3xy_S_vrr = PAX*I_ERI_Px_S_G3xy_S_vrr+WPX*I_ERI_Px_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_G3xy_S_vrr = PAY*I_ERI_Py_S_G3xy_S_vrr+WPY*I_ERI_Py_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_G3xy_S_vrr = PAZ*I_ERI_Pz_S_G3xy_S_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr;
      Double I_ERI_D2x_S_G3xz_S_vrr = PAX*I_ERI_Px_S_G3xz_S_vrr+WPX*I_ERI_Px_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_D2y_S_G3xz_S_vrr = PAY*I_ERI_Py_S_G3xz_S_vrr+WPY*I_ERI_Py_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr;
      Double I_ERI_D2z_S_G3xz_S_vrr = PAZ*I_ERI_Pz_S_G3xz_S_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2y_S_vrr = PAX*I_ERI_Px_S_G2x2y_S_vrr+WPX*I_ERI_Px_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2y_S_vrr = PAY*I_ERI_Py_S_G2x2y_S_vrr+WPY*I_ERI_Py_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2y_S_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr;
      Double I_ERI_D2x_S_G2xyz_S_vrr = PAX*I_ERI_Px_S_G2xyz_S_vrr+WPX*I_ERI_Px_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2y_S_G2xyz_S_vrr = PAY*I_ERI_Py_S_G2xyz_S_vrr+WPY*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_D2z_S_G2xyz_S_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2z_S_vrr = PAX*I_ERI_Px_S_G2x2z_S_vrr+WPX*I_ERI_Px_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2z_S_vrr = PAY*I_ERI_Py_S_G2x2z_S_vrr+WPY*I_ERI_Py_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2z_S_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3y_S_vrr = PAX*I_ERI_Px_S_Gx3y_S_vrr+WPX*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3y_S_vrr = PAY*I_ERI_Py_S_Gx3y_S_vrr+WPY*I_ERI_Py_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3y_S_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_vrr = PAX*I_ERI_Px_S_Gx2yz_S_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_vrr = PAY*I_ERI_Py_S_Gx2yz_S_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_vrr = PAX*I_ERI_Px_S_Gxy2z_S_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_vrr = PAY*I_ERI_Py_S_Gxy2z_S_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3z_S_vrr = PAX*I_ERI_Px_S_Gx3z_S_vrr+WPX*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3z_S_vrr = PAY*I_ERI_Py_S_Gx3z_S_vrr+WPY*I_ERI_Py_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3z_S_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4y_S_vrr = PAX*I_ERI_Px_S_G4y_S_vrr+WPX*I_ERI_Px_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_D2y_S_G4y_S_vrr = PAY*I_ERI_Py_S_G4y_S_vrr+WPY*I_ERI_Py_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_G4y_S_vrr = PAZ*I_ERI_Pz_S_G4y_S_vrr+WPZ*I_ERI_Pz_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_D2x_S_G3yz_S_vrr = PAX*I_ERI_Px_S_G3yz_S_vrr+WPX*I_ERI_Px_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr;
      Double I_ERI_D2y_S_G3yz_S_vrr = PAY*I_ERI_Py_S_G3yz_S_vrr+WPY*I_ERI_Py_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M1_vrr;
      Double I_ERI_D2z_S_G3yz_S_vrr = PAZ*I_ERI_Pz_S_G3yz_S_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_G2y2z_S_vrr = PAX*I_ERI_Px_S_G2y2z_S_vrr+WPX*I_ERI_Px_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr;
      Double I_ERI_D2y_S_G2y2z_S_vrr = PAY*I_ERI_Py_S_G2y2z_S_vrr+WPY*I_ERI_Py_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2z_S_G2y2z_S_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M1_vrr;
      Double I_ERI_D2x_S_Gy3z_S_vrr = PAX*I_ERI_Px_S_Gy3z_S_vrr+WPX*I_ERI_Px_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gy3z_S_vrr = PAY*I_ERI_Py_S_Gy3z_S_vrr+WPY*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gy3z_S_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4z_S_vrr = PAX*I_ERI_Px_S_G4z_S_vrr+WPX*I_ERI_Px_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_D2y_S_G4z_S_vrr = PAY*I_ERI_Py_S_G4z_S_vrr+WPY*I_ERI_Py_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_D2z_S_G4z_S_vrr = PAZ*I_ERI_Pz_S_G4z_S_vrr+WPZ*I_ERI_Pz_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_G_S
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_vrr = PAX*I_ERI_D2x_S_G4x_S_vrr+WPX*I_ERI_D2x_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xy_S_G4x_S_vrr = PAY*I_ERI_D2x_S_G4x_S_vrr+WPY*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_F2xz_S_G4x_S_vrr = PAZ*I_ERI_D2x_S_G4x_S_vrr+WPZ*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4x_S_vrr = PAX*I_ERI_D2y_S_G4x_S_vrr+WPX*I_ERI_D2y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4x_S_vrr = PAX*I_ERI_D2z_S_G4x_S_vrr+WPX*I_ERI_D2z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G4x_S_vrr = PAY*I_ERI_D2y_S_G4x_S_vrr+WPY*I_ERI_D2y_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4x_S_M1_vrr;
      Double I_ERI_F2yz_S_G4x_S_vrr = PAZ*I_ERI_D2y_S_G4x_S_vrr+WPZ*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_F3z_S_G4x_S_vrr = PAZ*I_ERI_D2z_S_G4x_S_vrr+WPZ*I_ERI_D2z_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4x_S_M1_vrr;
      Double I_ERI_F3x_S_G3xy_S_vrr = PAX*I_ERI_D2x_S_G3xy_S_vrr+WPX*I_ERI_D2x_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xy_S_vrr = PAY*I_ERI_D2x_S_G3xy_S_vrr+WPY*I_ERI_D2x_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xy_S_vrr = PAZ*I_ERI_D2x_S_G3xy_S_vrr+WPZ*I_ERI_D2x_S_G3xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_vrr = PAX*I_ERI_D2y_S_G3xy_S_vrr+WPX*I_ERI_D2y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_vrr = PAX*I_ERI_D2z_S_G3xy_S_vrr+WPX*I_ERI_D2z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3y_S_G3xy_S_vrr = PAY*I_ERI_D2y_S_G3xy_S_vrr+WPY*I_ERI_D2y_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xy_S_vrr = PAZ*I_ERI_D2y_S_G3xy_S_vrr+WPZ*I_ERI_D2y_S_G3xy_S_M1_vrr;
      Double I_ERI_F3z_S_G3xy_S_vrr = PAZ*I_ERI_D2z_S_G3xy_S_vrr+WPZ*I_ERI_D2z_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xy_S_M1_vrr;
      Double I_ERI_F3x_S_G3xz_S_vrr = PAX*I_ERI_D2x_S_G3xz_S_vrr+WPX*I_ERI_D2x_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xz_S_vrr = PAY*I_ERI_D2x_S_G3xz_S_vrr+WPY*I_ERI_D2x_S_G3xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xz_S_vrr = PAZ*I_ERI_D2x_S_G3xz_S_vrr+WPZ*I_ERI_D2x_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_vrr = PAX*I_ERI_D2y_S_G3xz_S_vrr+WPX*I_ERI_D2y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_vrr = PAX*I_ERI_D2z_S_G3xz_S_vrr+WPX*I_ERI_D2z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3y_S_G3xz_S_vrr = PAY*I_ERI_D2y_S_G3xz_S_vrr+WPY*I_ERI_D2y_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xz_S_vrr = PAZ*I_ERI_D2y_S_G3xz_S_vrr+WPZ*I_ERI_D2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_G3xz_S_vrr = PAZ*I_ERI_D2z_S_G3xz_S_vrr+WPZ*I_ERI_D2z_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2y_S_vrr = PAX*I_ERI_D2x_S_G2x2y_S_vrr+WPX*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_vrr = PAY*I_ERI_D2x_S_G2x2y_S_vrr+WPY*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_vrr = PAZ*I_ERI_D2x_S_G2x2y_S_vrr+WPZ*I_ERI_D2x_S_G2x2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_vrr = PAX*I_ERI_D2y_S_G2x2y_S_vrr+WPX*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_vrr = PAX*I_ERI_D2z_S_G2x2y_S_vrr+WPX*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2y_S_vrr = PAY*I_ERI_D2y_S_G2x2y_S_vrr+WPY*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_vrr = PAZ*I_ERI_D2y_S_G2x2y_S_vrr+WPZ*I_ERI_D2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2y_S_vrr = PAZ*I_ERI_D2z_S_G2x2y_S_vrr+WPZ*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2y_S_M1_vrr;
      Double I_ERI_F3x_S_G2xyz_S_vrr = PAX*I_ERI_D2x_S_G2xyz_S_vrr+WPX*I_ERI_D2x_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_vrr = PAY*I_ERI_D2x_S_G2xyz_S_vrr+WPY*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_vrr = PAZ*I_ERI_D2x_S_G2xyz_S_vrr+WPZ*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_vrr = PAX*I_ERI_D2y_S_G2xyz_S_vrr+WPX*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_vrr = PAX*I_ERI_D2z_S_G2xyz_S_vrr+WPX*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3y_S_G2xyz_S_vrr = PAY*I_ERI_D2y_S_G2xyz_S_vrr+WPY*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_vrr = PAZ*I_ERI_D2y_S_G2xyz_S_vrr+WPZ*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F3z_S_G2xyz_S_vrr = PAZ*I_ERI_D2z_S_G2xyz_S_vrr+WPZ*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2z_S_vrr = PAX*I_ERI_D2x_S_G2x2z_S_vrr+WPX*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_vrr = PAY*I_ERI_D2x_S_G2x2z_S_vrr+WPY*I_ERI_D2x_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_vrr = PAZ*I_ERI_D2x_S_G2x2z_S_vrr+WPZ*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_vrr = PAX*I_ERI_D2y_S_G2x2z_S_vrr+WPX*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_vrr = PAX*I_ERI_D2z_S_G2x2z_S_vrr+WPX*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2z_S_vrr = PAY*I_ERI_D2y_S_G2x2z_S_vrr+WPY*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_vrr = PAZ*I_ERI_D2y_S_G2x2z_S_vrr+WPZ*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2z_S_vrr = PAZ*I_ERI_D2z_S_G2x2z_S_vrr+WPZ*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3y_S_vrr = PAX*I_ERI_D2x_S_Gx3y_S_vrr+WPX*I_ERI_D2x_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_vrr = PAY*I_ERI_D2x_S_Gx3y_S_vrr+WPY*I_ERI_D2x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_vrr = PAZ*I_ERI_D2x_S_Gx3y_S_vrr+WPZ*I_ERI_D2x_S_Gx3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_vrr = PAX*I_ERI_D2y_S_Gx3y_S_vrr+WPX*I_ERI_D2y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_vrr = PAX*I_ERI_D2z_S_Gx3y_S_vrr+WPX*I_ERI_D2z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3y_S_vrr = PAY*I_ERI_D2y_S_Gx3y_S_vrr+WPY*I_ERI_D2y_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_vrr = PAZ*I_ERI_D2y_S_Gx3y_S_vrr+WPZ*I_ERI_D2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3y_S_vrr = PAZ*I_ERI_D2z_S_Gx3y_S_vrr+WPZ*I_ERI_D2z_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3y_S_M1_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_vrr = PAX*I_ERI_D2x_S_Gx2yz_S_vrr+WPX*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_vrr = PAY*I_ERI_D2x_S_Gx2yz_S_vrr+WPY*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2x_S_Gx2yz_S_vrr+WPZ*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_vrr = PAX*I_ERI_D2y_S_Gx2yz_S_vrr+WPX*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_vrr = PAX*I_ERI_D2z_S_Gx2yz_S_vrr+WPX*I_ERI_D2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_vrr = PAY*I_ERI_D2y_S_Gx2yz_S_vrr+WPY*I_ERI_D2y_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2y_S_Gx2yz_S_vrr+WPZ*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_vrr = PAZ*I_ERI_D2z_S_Gx2yz_S_vrr+WPZ*I_ERI_D2z_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_vrr = PAX*I_ERI_D2x_S_Gxy2z_S_vrr+WPX*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_vrr = PAY*I_ERI_D2x_S_Gxy2z_S_vrr+WPY*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2x_S_Gxy2z_S_vrr+WPZ*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_vrr = PAX*I_ERI_D2y_S_Gxy2z_S_vrr+WPX*I_ERI_D2y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_vrr = PAX*I_ERI_D2z_S_Gxy2z_S_vrr+WPX*I_ERI_D2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_vrr = PAY*I_ERI_D2y_S_Gxy2z_S_vrr+WPY*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2y_S_Gxy2z_S_vrr+WPZ*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_vrr = PAZ*I_ERI_D2z_S_Gxy2z_S_vrr+WPZ*I_ERI_D2z_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3z_S_vrr = PAX*I_ERI_D2x_S_Gx3z_S_vrr+WPX*I_ERI_D2x_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_vrr = PAY*I_ERI_D2x_S_Gx3z_S_vrr+WPY*I_ERI_D2x_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_vrr = PAZ*I_ERI_D2x_S_Gx3z_S_vrr+WPZ*I_ERI_D2x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_vrr = PAX*I_ERI_D2y_S_Gx3z_S_vrr+WPX*I_ERI_D2y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_vrr = PAX*I_ERI_D2z_S_Gx3z_S_vrr+WPX*I_ERI_D2z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3z_S_vrr = PAY*I_ERI_D2y_S_Gx3z_S_vrr+WPY*I_ERI_D2y_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_vrr = PAZ*I_ERI_D2y_S_Gx3z_S_vrr+WPZ*I_ERI_D2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3z_S_vrr = PAZ*I_ERI_D2z_S_Gx3z_S_vrr+WPZ*I_ERI_D2z_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4y_S_vrr = PAX*I_ERI_D2x_S_G4y_S_vrr+WPX*I_ERI_D2x_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4y_S_M1_vrr;
      Double I_ERI_F2xy_S_G4y_S_vrr = PAY*I_ERI_D2x_S_G4y_S_vrr+WPY*I_ERI_D2x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_G4y_S_vrr = PAZ*I_ERI_D2x_S_G4y_S_vrr+WPZ*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4y_S_vrr = PAX*I_ERI_D2y_S_G4y_S_vrr+WPX*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4y_S_vrr = PAX*I_ERI_D2z_S_G4y_S_vrr+WPX*I_ERI_D2z_S_G4y_S_M1_vrr;
      Double I_ERI_F3y_S_G4y_S_vrr = PAY*I_ERI_D2y_S_G4y_S_vrr+WPY*I_ERI_D2y_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_G4y_S_vrr = PAZ*I_ERI_D2y_S_G4y_S_vrr+WPZ*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_F3z_S_G4y_S_vrr = PAZ*I_ERI_D2z_S_G4y_S_vrr+WPZ*I_ERI_D2z_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4y_S_M1_vrr;
      Double I_ERI_F3x_S_G3yz_S_vrr = PAX*I_ERI_D2x_S_G3yz_S_vrr+WPX*I_ERI_D2x_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3yz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3yz_S_vrr = PAY*I_ERI_D2x_S_G3yz_S_vrr+WPY*I_ERI_D2x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3yz_S_vrr = PAZ*I_ERI_D2x_S_G3yz_S_vrr+WPZ*I_ERI_D2x_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_vrr = PAX*I_ERI_D2y_S_G3yz_S_vrr+WPX*I_ERI_D2y_S_G3yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_vrr = PAX*I_ERI_D2z_S_G3yz_S_vrr+WPX*I_ERI_D2z_S_G3yz_S_M1_vrr;
      Double I_ERI_F3y_S_G3yz_S_vrr = PAY*I_ERI_D2y_S_G3yz_S_vrr+WPY*I_ERI_D2y_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3yz_S_vrr = PAZ*I_ERI_D2y_S_G3yz_S_vrr+WPZ*I_ERI_D2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_G3yz_S_vrr = PAZ*I_ERI_D2z_S_G3yz_S_vrr+WPZ*I_ERI_D2z_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_G2y2z_S_vrr = PAX*I_ERI_D2x_S_G2y2z_S_vrr+WPX*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2y2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_vrr = PAY*I_ERI_D2x_S_G2y2z_S_vrr+WPY*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_vrr = PAZ*I_ERI_D2x_S_G2y2z_S_vrr+WPZ*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_vrr = PAX*I_ERI_D2y_S_G2y2z_S_vrr+WPX*I_ERI_D2y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_vrr = PAX*I_ERI_D2z_S_G2y2z_S_vrr+WPX*I_ERI_D2z_S_G2y2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2y2z_S_vrr = PAY*I_ERI_D2y_S_G2y2z_S_vrr+WPY*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_vrr = PAZ*I_ERI_D2y_S_G2y2z_S_vrr+WPZ*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_F3z_S_G2y2z_S_vrr = PAZ*I_ERI_D2z_S_G2y2z_S_vrr+WPZ*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3x_S_Gy3z_S_vrr = PAX*I_ERI_D2x_S_Gy3z_S_vrr+WPX*I_ERI_D2x_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gy3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_vrr = PAY*I_ERI_D2x_S_Gy3z_S_vrr+WPY*I_ERI_D2x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_vrr = PAZ*I_ERI_D2x_S_Gy3z_S_vrr+WPZ*I_ERI_D2x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_vrr = PAX*I_ERI_D2y_S_Gy3z_S_vrr+WPX*I_ERI_D2y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_vrr = PAX*I_ERI_D2z_S_Gy3z_S_vrr+WPX*I_ERI_D2z_S_Gy3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gy3z_S_vrr = PAY*I_ERI_D2y_S_Gy3z_S_vrr+WPY*I_ERI_D2y_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_vrr = PAZ*I_ERI_D2y_S_Gy3z_S_vrr+WPZ*I_ERI_D2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gy3z_S_vrr = PAZ*I_ERI_D2z_S_Gy3z_S_vrr+WPZ*I_ERI_D2z_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4z_S_vrr = PAX*I_ERI_D2x_S_G4z_S_vrr+WPX*I_ERI_D2x_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4z_S_M1_vrr;
      Double I_ERI_F2xy_S_G4z_S_vrr = PAY*I_ERI_D2x_S_G4z_S_vrr+WPY*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_F2xz_S_G4z_S_vrr = PAZ*I_ERI_D2x_S_G4z_S_vrr+WPZ*I_ERI_D2x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4z_S_vrr = PAX*I_ERI_D2y_S_G4z_S_vrr+WPX*I_ERI_D2y_S_G4z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4z_S_vrr = PAX*I_ERI_D2z_S_G4z_S_vrr+WPX*I_ERI_D2z_S_G4z_S_M1_vrr;
      Double I_ERI_F3y_S_G4z_S_vrr = PAY*I_ERI_D2y_S_G4z_S_vrr+WPY*I_ERI_D2y_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4z_S_M1_vrr;
      Double I_ERI_F2yz_S_G4z_S_vrr = PAZ*I_ERI_D2y_S_G4z_S_vrr+WPZ*I_ERI_D2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_G4z_S_vrr = PAZ*I_ERI_D2z_S_G4z_S_vrr+WPZ*I_ERI_D2z_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_G_S
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_G_S
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_G4x_S_vrr = PAX*I_ERI_F3x_S_G4x_S_vrr+WPX*I_ERI_F3x_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xy_S_G4x_S_vrr = PAY*I_ERI_F3x_S_G4x_S_vrr+WPY*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_G3xz_S_G4x_S_vrr = PAZ*I_ERI_F3x_S_G4x_S_vrr+WPZ*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4x_S_vrr = PAY*I_ERI_F2xy_S_G4x_S_vrr+WPY*I_ERI_F2xy_S_G4x_S_M1_vrr+oned2z*I_ERI_D2x_S_G4x_S_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4x_S_vrr = PAZ*I_ERI_F2xy_S_G4x_S_vrr+WPZ*I_ERI_F2xy_S_G4x_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4x_S_vrr = PAZ*I_ERI_F2xz_S_G4x_S_vrr+WPZ*I_ERI_F2xz_S_G4x_S_M1_vrr+oned2z*I_ERI_D2x_S_G4x_S_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4x_S_vrr = PAX*I_ERI_F3y_S_G4x_S_vrr+WPX*I_ERI_F3y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4x_S_vrr = PAZ*I_ERI_Fx2y_S_G4x_S_vrr+WPZ*I_ERI_Fx2y_S_G4x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4x_S_vrr = PAY*I_ERI_Fx2z_S_G4x_S_vrr+WPY*I_ERI_Fx2z_S_G4x_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4x_S_vrr = PAX*I_ERI_F3z_S_G4x_S_vrr+WPX*I_ERI_F3z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4y_S_G4x_S_vrr = PAY*I_ERI_F3y_S_G4x_S_vrr+WPY*I_ERI_F3y_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_G3yz_S_G4x_S_vrr = PAZ*I_ERI_F3y_S_G4x_S_vrr+WPZ*I_ERI_F3y_S_G4x_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4x_S_vrr = PAZ*I_ERI_F2yz_S_G4x_S_vrr+WPZ*I_ERI_F2yz_S_G4x_S_M1_vrr+oned2z*I_ERI_D2y_S_G4x_S_vrr-rhod2zsq*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4x_S_vrr = PAY*I_ERI_F3z_S_G4x_S_vrr+WPY*I_ERI_F3z_S_G4x_S_M1_vrr;
      Double I_ERI_G4z_S_G4x_S_vrr = PAZ*I_ERI_F3z_S_G4x_S_vrr+WPZ*I_ERI_F3z_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4x_S_M1_vrr;
      Double I_ERI_G4x_S_G3xy_S_vrr = PAX*I_ERI_F3x_S_G3xy_S_vrr+WPX*I_ERI_F3x_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G3xy_S_G3xy_S_vrr = PAY*I_ERI_F3x_S_G3xy_S_vrr+WPY*I_ERI_F3x_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xz_S_G3xy_S_vrr = PAZ*I_ERI_F3x_S_G3xy_S_vrr+WPZ*I_ERI_F3x_S_G3xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3xy_S_vrr = PAY*I_ERI_F2xy_S_G3xy_S_vrr+WPY*I_ERI_F2xy_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3xy_S_vrr = PAZ*I_ERI_F2xy_S_G3xy_S_vrr+WPZ*I_ERI_F2xy_S_G3xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3xy_S_vrr = PAZ*I_ERI_F2xz_S_G3xy_S_vrr+WPZ*I_ERI_F2xz_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3xy_S_vrr = PAX*I_ERI_F3y_S_G3xy_S_vrr+WPX*I_ERI_F3y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3xy_S_vrr = PAZ*I_ERI_Fx2y_S_G3xy_S_vrr+WPZ*I_ERI_Fx2y_S_G3xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3xy_S_vrr = PAY*I_ERI_Fx2z_S_G3xy_S_vrr+WPY*I_ERI_Fx2z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3xy_S_vrr = PAX*I_ERI_F3z_S_G3xy_S_vrr+WPX*I_ERI_F3z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4y_S_G3xy_S_vrr = PAY*I_ERI_F3y_S_G3xy_S_vrr+WPY*I_ERI_F3y_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G3yz_S_G3xy_S_vrr = PAZ*I_ERI_F3y_S_G3xy_S_vrr+WPZ*I_ERI_F3y_S_G3xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3xy_S_vrr = PAZ*I_ERI_F2yz_S_G3xy_S_vrr+WPZ*I_ERI_F2yz_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2y_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2y_S_G3xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3xy_S_vrr = PAY*I_ERI_F3z_S_G3xy_S_vrr+WPY*I_ERI_F3z_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4z_S_G3xy_S_vrr = PAZ*I_ERI_F3z_S_G3xy_S_vrr+WPZ*I_ERI_F3z_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xy_S_M1_vrr;
      Double I_ERI_G4x_S_G3xz_S_vrr = PAX*I_ERI_F3x_S_G3xz_S_vrr+WPX*I_ERI_F3x_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xy_S_G3xz_S_vrr = PAY*I_ERI_F3x_S_G3xz_S_vrr+WPY*I_ERI_F3x_S_G3xz_S_M1_vrr;
      Double I_ERI_G3xz_S_G3xz_S_vrr = PAZ*I_ERI_F3x_S_G3xz_S_vrr+WPZ*I_ERI_F3x_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3xz_S_vrr = PAY*I_ERI_F2xy_S_G3xz_S_vrr+WPY*I_ERI_F2xy_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3xz_S_vrr = PAZ*I_ERI_F2xy_S_G3xz_S_vrr+WPZ*I_ERI_F2xy_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3xz_S_vrr = PAZ*I_ERI_F2xz_S_G3xz_S_vrr+WPZ*I_ERI_F2xz_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3xz_S_vrr = PAX*I_ERI_F3y_S_G3xz_S_vrr+WPX*I_ERI_F3y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3xz_S_vrr = PAZ*I_ERI_Fx2y_S_G3xz_S_vrr+WPZ*I_ERI_Fx2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3xz_S_vrr = PAY*I_ERI_Fx2z_S_G3xz_S_vrr+WPY*I_ERI_Fx2z_S_G3xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3xz_S_vrr = PAX*I_ERI_F3z_S_G3xz_S_vrr+WPX*I_ERI_F3z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4y_S_G3xz_S_vrr = PAY*I_ERI_F3y_S_G3xz_S_vrr+WPY*I_ERI_F3y_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xz_S_M1_vrr;
      Double I_ERI_G3yz_S_G3xz_S_vrr = PAZ*I_ERI_F3y_S_G3xz_S_vrr+WPZ*I_ERI_F3y_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3xz_S_vrr = PAZ*I_ERI_F2yz_S_G3xz_S_vrr+WPZ*I_ERI_F2yz_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2y_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F3x_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3xz_S_vrr = PAY*I_ERI_F3z_S_G3xz_S_vrr+WPY*I_ERI_F3z_S_G3xz_S_M1_vrr;
      Double I_ERI_G4z_S_G3xz_S_vrr = PAZ*I_ERI_F3z_S_G3xz_S_vrr+WPZ*I_ERI_F3z_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4x_S_G2x2y_S_vrr = PAX*I_ERI_F3x_S_G2x2y_S_vrr+WPX*I_ERI_F3x_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3xy_S_G2x2y_S_vrr = PAY*I_ERI_F3x_S_G2x2y_S_vrr+WPY*I_ERI_F3x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G3xz_S_G2x2y_S_vrr = PAZ*I_ERI_F3x_S_G2x2y_S_vrr+WPZ*I_ERI_F3x_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2x2y_S_vrr = PAY*I_ERI_F2xy_S_G2x2y_S_vrr+WPY*I_ERI_F2xy_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2x2y_S_vrr = PAZ*I_ERI_F2xy_S_G2x2y_S_vrr+WPZ*I_ERI_F2xy_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2x2y_S_vrr = PAZ*I_ERI_F2xz_S_G2x2y_S_vrr+WPZ*I_ERI_F2xz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2x2y_S_vrr = PAX*I_ERI_F3y_S_G2x2y_S_vrr+WPX*I_ERI_F3y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2x2y_S_vrr = PAZ*I_ERI_Fx2y_S_G2x2y_S_vrr+WPZ*I_ERI_Fx2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2x2y_S_vrr = PAY*I_ERI_Fx2z_S_G2x2y_S_vrr+WPY*I_ERI_Fx2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2x2y_S_vrr = PAX*I_ERI_F3z_S_G2x2y_S_vrr+WPX*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4y_S_G2x2y_S_vrr = PAY*I_ERI_F3y_S_G2x2y_S_vrr+WPY*I_ERI_F3y_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G3yz_S_G2x2y_S_vrr = PAZ*I_ERI_F3y_S_G2x2y_S_vrr+WPZ*I_ERI_F3y_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2x2y_S_vrr = PAZ*I_ERI_F2yz_S_G2x2y_S_vrr+WPZ*I_ERI_F2yz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2x2y_S_vrr = PAY*I_ERI_F3z_S_G2x2y_S_vrr+WPY*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4z_S_G2x2y_S_vrr = PAZ*I_ERI_F3z_S_G2x2y_S_vrr+WPZ*I_ERI_F3z_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2y_S_M1_vrr;
      Double I_ERI_G4x_S_G2xyz_S_vrr = PAX*I_ERI_F3x_S_G2xyz_S_vrr+WPX*I_ERI_F3x_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3xy_S_G2xyz_S_vrr = PAY*I_ERI_F3x_S_G2xyz_S_vrr+WPY*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xz_S_G2xyz_S_vrr = PAZ*I_ERI_F3x_S_G2xyz_S_vrr+WPZ*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2xyz_S_vrr = PAY*I_ERI_F2xy_S_G2xyz_S_vrr+WPY*I_ERI_F2xy_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F2xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2xyz_S_vrr = PAZ*I_ERI_F2xy_S_G2xyz_S_vrr+WPZ*I_ERI_F2xy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2xyz_S_vrr = PAZ*I_ERI_F2xz_S_G2xyz_S_vrr+WPZ*I_ERI_F2xz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2xyz_S_vrr = PAX*I_ERI_F3y_S_G2xyz_S_vrr+WPX*I_ERI_F3y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2xyz_S_vrr = PAZ*I_ERI_Fx2y_S_G2xyz_S_vrr+WPZ*I_ERI_Fx2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2xyz_S_vrr = PAY*I_ERI_Fx2z_S_G2xyz_S_vrr+WPY*I_ERI_Fx2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2xyz_S_vrr = PAX*I_ERI_F3z_S_G2xyz_S_vrr+WPX*I_ERI_F3z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4y_S_G2xyz_S_vrr = PAY*I_ERI_F3y_S_G2xyz_S_vrr+WPY*I_ERI_F3y_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_G3yz_S_G2xyz_S_vrr = PAZ*I_ERI_F3y_S_G2xyz_S_vrr+WPZ*I_ERI_F3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2xyz_S_vrr = PAZ*I_ERI_F2yz_S_G2xyz_S_vrr+WPZ*I_ERI_F2yz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F2xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2xyz_S_vrr = PAY*I_ERI_F3z_S_G2xyz_S_vrr+WPY*I_ERI_F3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4z_S_G2xyz_S_vrr = PAZ*I_ERI_F3z_S_G2xyz_S_vrr+WPZ*I_ERI_F3z_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4x_S_G2x2z_S_vrr = PAX*I_ERI_F3x_S_G2x2z_S_vrr+WPX*I_ERI_F3x_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xy_S_G2x2z_S_vrr = PAY*I_ERI_F3x_S_G2x2z_S_vrr+WPY*I_ERI_F3x_S_G2x2z_S_M1_vrr;
      Double I_ERI_G3xz_S_G2x2z_S_vrr = PAZ*I_ERI_F3x_S_G2x2z_S_vrr+WPZ*I_ERI_F3x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2x2z_S_vrr = PAY*I_ERI_F2xy_S_G2x2z_S_vrr+WPY*I_ERI_F2xy_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2x2z_S_vrr = PAZ*I_ERI_F2xy_S_G2x2z_S_vrr+WPZ*I_ERI_F2xy_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2xz_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2x2z_S_vrr = PAZ*I_ERI_F2xz_S_G2x2z_S_vrr+WPZ*I_ERI_F2xz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2x2z_S_vrr = PAX*I_ERI_F3y_S_G2x2z_S_vrr+WPX*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2x2z_S_vrr = PAZ*I_ERI_Fx2y_S_G2x2z_S_vrr+WPZ*I_ERI_Fx2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2x2z_S_vrr = PAY*I_ERI_Fx2z_S_G2x2z_S_vrr+WPY*I_ERI_Fx2z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2x2z_S_vrr = PAX*I_ERI_F3z_S_G2x2z_S_vrr+WPX*I_ERI_F3z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4y_S_G2x2z_S_vrr = PAY*I_ERI_F3y_S_G2x2z_S_vrr+WPY*I_ERI_F3y_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2z_S_M1_vrr;
      Double I_ERI_G3yz_S_G2x2z_S_vrr = PAZ*I_ERI_F3y_S_G2x2z_S_vrr+WPZ*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2x2z_S_vrr = PAZ*I_ERI_F2yz_S_G2x2z_S_vrr+WPZ*I_ERI_F2yz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_F2xz_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2x2z_S_vrr = PAY*I_ERI_F3z_S_G2x2z_S_vrr+WPY*I_ERI_F3z_S_G2x2z_S_M1_vrr;
      Double I_ERI_G4z_S_G2x2z_S_vrr = PAZ*I_ERI_F3z_S_G2x2z_S_vrr+WPZ*I_ERI_F3z_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4x_S_Gx3y_S_vrr = PAX*I_ERI_F3x_S_Gx3y_S_vrr+WPX*I_ERI_F3x_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx3y_S_vrr = PAY*I_ERI_F3x_S_Gx3y_S_vrr+WPY*I_ERI_F3x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx3y_S_vrr = PAZ*I_ERI_F3x_S_Gx3y_S_vrr+WPZ*I_ERI_F3x_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx3y_S_vrr = PAY*I_ERI_F2xy_S_Gx3y_S_vrr+WPY*I_ERI_F2xy_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx3y_S_vrr = PAZ*I_ERI_F2xy_S_Gx3y_S_vrr+WPZ*I_ERI_F2xy_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx3y_S_vrr = PAZ*I_ERI_F2xz_S_Gx3y_S_vrr+WPZ*I_ERI_F2xz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx3y_S_vrr = PAX*I_ERI_F3y_S_Gx3y_S_vrr+WPX*I_ERI_F3y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx3y_S_vrr = PAZ*I_ERI_Fx2y_S_Gx3y_S_vrr+WPZ*I_ERI_Fx2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx3y_S_vrr = PAY*I_ERI_Fx2z_S_Gx3y_S_vrr+WPY*I_ERI_Fx2z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx3y_S_vrr = PAX*I_ERI_F3z_S_Gx3y_S_vrr+WPX*I_ERI_F3z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4y_S_Gx3y_S_vrr = PAY*I_ERI_F3y_S_Gx3y_S_vrr+WPY*I_ERI_F3y_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx3y_S_vrr = PAZ*I_ERI_F3y_S_Gx3y_S_vrr+WPZ*I_ERI_F3y_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx3y_S_vrr = PAZ*I_ERI_F2yz_S_Gx3y_S_vrr+WPZ*I_ERI_F2yz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx3y_S_vrr = PAY*I_ERI_F3z_S_Gx3y_S_vrr+WPY*I_ERI_F3z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4z_S_Gx3y_S_vrr = PAZ*I_ERI_F3z_S_Gx3y_S_vrr+WPZ*I_ERI_F3z_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3y_S_M1_vrr;
      Double I_ERI_G4x_S_Gx2yz_S_vrr = PAX*I_ERI_F3x_S_Gx2yz_S_vrr+WPX*I_ERI_F3x_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx2yz_S_vrr = PAY*I_ERI_F3x_S_Gx2yz_S_vrr+WPY*I_ERI_F3x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx2yz_S_vrr = PAZ*I_ERI_F3x_S_Gx2yz_S_vrr+WPZ*I_ERI_F3x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx2yz_S_vrr = PAY*I_ERI_F2xy_S_Gx2yz_S_vrr+WPY*I_ERI_F2xy_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx2yz_S_vrr = PAZ*I_ERI_F2xy_S_Gx2yz_S_vrr+WPZ*I_ERI_F2xy_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx2yz_S_vrr = PAZ*I_ERI_F2xz_S_Gx2yz_S_vrr+WPZ*I_ERI_F2xz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx2yz_S_vrr = PAX*I_ERI_F3y_S_Gx2yz_S_vrr+WPX*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_Fx2y_S_Gx2yz_S_vrr+WPZ*I_ERI_Fx2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx2yz_S_vrr = PAY*I_ERI_Fx2z_S_Gx2yz_S_vrr+WPY*I_ERI_Fx2z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx2yz_S_vrr = PAX*I_ERI_F3z_S_Gx2yz_S_vrr+WPX*I_ERI_F3z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4y_S_Gx2yz_S_vrr = PAY*I_ERI_F3y_S_Gx2yz_S_vrr+WPY*I_ERI_F3y_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx2yz_S_vrr = PAZ*I_ERI_F3y_S_Gx2yz_S_vrr+WPZ*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx2yz_S_vrr = PAZ*I_ERI_F2yz_S_Gx2yz_S_vrr+WPZ*I_ERI_F2yz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx2yz_S_vrr = PAY*I_ERI_F3z_S_Gx2yz_S_vrr+WPY*I_ERI_F3z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4z_S_Gx2yz_S_vrr = PAZ*I_ERI_F3z_S_Gx2yz_S_vrr+WPZ*I_ERI_F3z_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4x_S_Gxy2z_S_vrr = PAX*I_ERI_F3x_S_Gxy2z_S_vrr+WPX*I_ERI_F3x_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gxy2z_S_vrr = PAY*I_ERI_F3x_S_Gxy2z_S_vrr+WPY*I_ERI_F3x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gxy2z_S_vrr = PAZ*I_ERI_F3x_S_Gxy2z_S_vrr+WPZ*I_ERI_F3x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gxy2z_S_vrr = PAY*I_ERI_F2xy_S_Gxy2z_S_vrr+WPY*I_ERI_F2xy_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F2xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gxy2z_S_vrr = PAZ*I_ERI_F2xy_S_Gxy2z_S_vrr+WPZ*I_ERI_F2xy_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gxy2z_S_vrr = PAZ*I_ERI_F2xz_S_Gxy2z_S_vrr+WPZ*I_ERI_F2xz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gxy2z_S_vrr = PAX*I_ERI_F3y_S_Gxy2z_S_vrr+WPX*I_ERI_F3y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_Fx2y_S_Gxy2z_S_vrr+WPZ*I_ERI_Fx2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gxy2z_S_vrr = PAY*I_ERI_Fx2z_S_Gxy2z_S_vrr+WPY*I_ERI_Fx2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gxy2z_S_vrr = PAX*I_ERI_F3z_S_Gxy2z_S_vrr+WPX*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4y_S_Gxy2z_S_vrr = PAY*I_ERI_F3y_S_Gxy2z_S_vrr+WPY*I_ERI_F3y_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gxy2z_S_vrr = PAZ*I_ERI_F3y_S_Gxy2z_S_vrr+WPZ*I_ERI_F3y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gxy2z_S_vrr = PAZ*I_ERI_F2yz_S_Gxy2z_S_vrr+WPZ*I_ERI_F2yz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gxy2z_S_vrr = PAY*I_ERI_F3z_S_Gxy2z_S_vrr+WPY*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4z_S_Gxy2z_S_vrr = PAZ*I_ERI_F3z_S_Gxy2z_S_vrr+WPZ*I_ERI_F3z_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4x_S_Gx3z_S_vrr = PAX*I_ERI_F3x_S_Gx3z_S_vrr+WPX*I_ERI_F3x_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx3z_S_vrr = PAY*I_ERI_F3x_S_Gx3z_S_vrr+WPY*I_ERI_F3x_S_Gx3z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx3z_S_vrr = PAZ*I_ERI_F3x_S_Gx3z_S_vrr+WPZ*I_ERI_F3x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx3z_S_vrr = PAY*I_ERI_F2xy_S_Gx3z_S_vrr+WPY*I_ERI_F2xy_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx3z_S_vrr = PAZ*I_ERI_F2xy_S_Gx3z_S_vrr+WPZ*I_ERI_F2xy_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx3z_S_vrr = PAZ*I_ERI_F2xz_S_Gx3z_S_vrr+WPZ*I_ERI_F2xz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx3z_S_vrr = PAX*I_ERI_F3y_S_Gx3z_S_vrr+WPX*I_ERI_F3y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx3z_S_vrr = PAZ*I_ERI_Fx2y_S_Gx3z_S_vrr+WPZ*I_ERI_Fx2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx3z_S_vrr = PAY*I_ERI_Fx2z_S_Gx3z_S_vrr+WPY*I_ERI_Fx2z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx3z_S_vrr = PAX*I_ERI_F3z_S_Gx3z_S_vrr+WPX*I_ERI_F3z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4y_S_Gx3z_S_vrr = PAY*I_ERI_F3y_S_Gx3z_S_vrr+WPY*I_ERI_F3y_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx3z_S_vrr = PAZ*I_ERI_F3y_S_Gx3z_S_vrr+WPZ*I_ERI_F3y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx3z_S_vrr = PAZ*I_ERI_F2yz_S_Gx3z_S_vrr+WPZ*I_ERI_F2yz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx3z_S_vrr = PAY*I_ERI_F3z_S_Gx3z_S_vrr+WPY*I_ERI_F3z_S_Gx3z_S_M1_vrr;
      Double I_ERI_G4z_S_Gx3z_S_vrr = PAZ*I_ERI_F3z_S_Gx3z_S_vrr+WPZ*I_ERI_F3z_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4x_S_G4y_S_vrr = PAX*I_ERI_F3x_S_G4y_S_vrr+WPX*I_ERI_F3x_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_G3xy_S_G4y_S_vrr = PAY*I_ERI_F3x_S_G4y_S_vrr+WPY*I_ERI_F3x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xz_S_G4y_S_vrr = PAZ*I_ERI_F3x_S_G4y_S_vrr+WPZ*I_ERI_F3x_S_G4y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4y_S_vrr = PAY*I_ERI_F2xy_S_G4y_S_vrr+WPY*I_ERI_F2xy_S_G4y_S_M1_vrr+oned2z*I_ERI_D2x_S_G4y_S_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4y_S_vrr = PAZ*I_ERI_F2xy_S_G4y_S_vrr+WPZ*I_ERI_F2xy_S_G4y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4y_S_vrr = PAZ*I_ERI_F2xz_S_G4y_S_vrr+WPZ*I_ERI_F2xz_S_G4y_S_M1_vrr+oned2z*I_ERI_D2x_S_G4y_S_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4y_S_vrr = PAX*I_ERI_F3y_S_G4y_S_vrr+WPX*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4y_S_vrr = PAZ*I_ERI_Fx2y_S_G4y_S_vrr+WPZ*I_ERI_Fx2y_S_G4y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4y_S_vrr = PAY*I_ERI_Fx2z_S_G4y_S_vrr+WPY*I_ERI_Fx2z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Fx2z_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4y_S_vrr = PAX*I_ERI_F3z_S_G4y_S_vrr+WPX*I_ERI_F3z_S_G4y_S_M1_vrr;
      Double I_ERI_G4y_S_G4y_S_vrr = PAY*I_ERI_F3y_S_G4y_S_vrr+WPY*I_ERI_F3y_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G3yz_S_G4y_S_vrr = PAZ*I_ERI_F3y_S_G4y_S_vrr+WPZ*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4y_S_vrr = PAZ*I_ERI_F2yz_S_G4y_S_vrr+WPZ*I_ERI_F2yz_S_G4y_S_M1_vrr+oned2z*I_ERI_D2y_S_G4y_S_vrr-rhod2zsq*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4y_S_vrr = PAY*I_ERI_F3z_S_G4y_S_vrr+WPY*I_ERI_F3z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4z_S_G4y_S_vrr = PAZ*I_ERI_F3z_S_G4y_S_vrr+WPZ*I_ERI_F3z_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4y_S_M1_vrr;
      Double I_ERI_G4x_S_G3yz_S_vrr = PAX*I_ERI_F3x_S_G3yz_S_vrr+WPX*I_ERI_F3x_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr;
      Double I_ERI_G3xy_S_G3yz_S_vrr = PAY*I_ERI_F3x_S_G3yz_S_vrr+WPY*I_ERI_F3x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xz_S_G3yz_S_vrr = PAZ*I_ERI_F3x_S_G3yz_S_vrr+WPZ*I_ERI_F3x_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3yz_S_vrr = PAY*I_ERI_F2xy_S_G3yz_S_vrr+WPY*I_ERI_F2xy_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_F2yz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3yz_S_vrr = PAZ*I_ERI_F2xy_S_G3yz_S_vrr+WPZ*I_ERI_F2xy_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3yz_S_vrr = PAZ*I_ERI_F2xz_S_G3yz_S_vrr+WPZ*I_ERI_F2xz_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3yz_S_vrr = PAX*I_ERI_F3y_S_G3yz_S_vrr+WPX*I_ERI_F3y_S_G3yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3yz_S_vrr = PAZ*I_ERI_Fx2y_S_G3yz_S_vrr+WPZ*I_ERI_Fx2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3yz_S_vrr = PAY*I_ERI_Fx2z_S_G3yz_S_vrr+WPY*I_ERI_Fx2z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3yz_S_vrr = PAX*I_ERI_F3z_S_G3yz_S_vrr+WPX*I_ERI_F3z_S_G3yz_S_M1_vrr;
      Double I_ERI_G4y_S_G3yz_S_vrr = PAY*I_ERI_F3y_S_G3yz_S_vrr+WPY*I_ERI_F3y_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_G3yz_S_G3yz_S_vrr = PAZ*I_ERI_F3y_S_G3yz_S_vrr+WPZ*I_ERI_F3y_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3yz_S_vrr = PAZ*I_ERI_F2yz_S_G3yz_S_vrr+WPZ*I_ERI_F2yz_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2y_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3yz_S_vrr = PAY*I_ERI_F3z_S_G3yz_S_vrr+WPY*I_ERI_F3z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4z_S_G3yz_S_vrr = PAZ*I_ERI_F3z_S_G3yz_S_vrr+WPZ*I_ERI_F3z_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4x_S_G2y2z_S_vrr = PAX*I_ERI_F3x_S_G2y2z_S_vrr+WPX*I_ERI_F3x_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr;
      Double I_ERI_G3xy_S_G2y2z_S_vrr = PAY*I_ERI_F3x_S_G2y2z_S_vrr+WPY*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xz_S_G2y2z_S_vrr = PAZ*I_ERI_F3x_S_G2y2z_S_vrr+WPZ*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2y2z_S_vrr = PAY*I_ERI_F2xy_S_G2y2z_S_vrr+WPY*I_ERI_F2xy_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2y2z_S_vrr = PAZ*I_ERI_F2xy_S_G2y2z_S_vrr+WPZ*I_ERI_F2xy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2yz_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2y2z_S_vrr = PAZ*I_ERI_F2xz_S_G2y2z_S_vrr+WPZ*I_ERI_F2xz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2y2z_S_vrr = PAX*I_ERI_F3y_S_G2y2z_S_vrr+WPX*I_ERI_F3y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2y2z_S_vrr = PAZ*I_ERI_Fx2y_S_G2y2z_S_vrr+WPZ*I_ERI_Fx2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2y2z_S_vrr = PAY*I_ERI_Fx2z_S_G2y2z_S_vrr+WPY*I_ERI_Fx2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2y2z_S_vrr = PAX*I_ERI_F3z_S_G2y2z_S_vrr+WPX*I_ERI_F3z_S_G2y2z_S_M1_vrr;
      Double I_ERI_G4y_S_G2y2z_S_vrr = PAY*I_ERI_F3y_S_G2y2z_S_vrr+WPY*I_ERI_F3y_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3yz_S_G2y2z_S_vrr = PAZ*I_ERI_F3y_S_G2y2z_S_vrr+WPZ*I_ERI_F3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2y2z_S_vrr = PAZ*I_ERI_F2yz_S_G2y2z_S_vrr+WPZ*I_ERI_F2yz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_F2yz_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2y2z_S_vrr = PAY*I_ERI_F3z_S_G2y2z_S_vrr+WPY*I_ERI_F3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4z_S_G2y2z_S_vrr = PAZ*I_ERI_F3z_S_G2y2z_S_vrr+WPZ*I_ERI_F3z_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4x_S_Gy3z_S_vrr = PAX*I_ERI_F3x_S_Gy3z_S_vrr+WPX*I_ERI_F3x_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gy3z_S_vrr = PAY*I_ERI_F3x_S_Gy3z_S_vrr+WPY*I_ERI_F3x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gy3z_S_vrr = PAZ*I_ERI_F3x_S_Gy3z_S_vrr+WPZ*I_ERI_F3x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gy3z_S_vrr = PAY*I_ERI_F2xy_S_Gy3z_S_vrr+WPY*I_ERI_F2xy_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gy3z_S_vrr = PAZ*I_ERI_F2xy_S_Gy3z_S_vrr+WPZ*I_ERI_F2xy_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gy3z_S_vrr = PAZ*I_ERI_F2xz_S_Gy3z_S_vrr+WPZ*I_ERI_F2xz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gy3z_S_vrr = PAX*I_ERI_F3y_S_Gy3z_S_vrr+WPX*I_ERI_F3y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gy3z_S_vrr = PAZ*I_ERI_Fx2y_S_Gy3z_S_vrr+WPZ*I_ERI_Fx2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gy3z_S_vrr = PAY*I_ERI_Fx2z_S_Gy3z_S_vrr+WPY*I_ERI_Fx2z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gy3z_S_vrr = PAX*I_ERI_F3z_S_Gy3z_S_vrr+WPX*I_ERI_F3z_S_Gy3z_S_M1_vrr;
      Double I_ERI_G4y_S_Gy3z_S_vrr = PAY*I_ERI_F3y_S_Gy3z_S_vrr+WPY*I_ERI_F3y_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gy3z_S_vrr = PAZ*I_ERI_F3y_S_Gy3z_S_vrr+WPZ*I_ERI_F3y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gy3z_S_vrr = PAZ*I_ERI_F2yz_S_Gy3z_S_vrr+WPZ*I_ERI_F2yz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gy3z_S_vrr = PAY*I_ERI_F3z_S_Gy3z_S_vrr+WPY*I_ERI_F3z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4z_S_Gy3z_S_vrr = PAZ*I_ERI_F3z_S_Gy3z_S_vrr+WPZ*I_ERI_F3z_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4x_S_G4z_S_vrr = PAX*I_ERI_F3x_S_G4z_S_vrr+WPX*I_ERI_F3x_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_G3xy_S_G4z_S_vrr = PAY*I_ERI_F3x_S_G4z_S_vrr+WPY*I_ERI_F3x_S_G4z_S_M1_vrr;
      Double I_ERI_G3xz_S_G4z_S_vrr = PAZ*I_ERI_F3x_S_G4z_S_vrr+WPZ*I_ERI_F3x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4z_S_vrr = PAY*I_ERI_F2xy_S_G4z_S_vrr+WPY*I_ERI_F2xy_S_G4z_S_M1_vrr+oned2z*I_ERI_D2x_S_G4z_S_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4z_S_vrr = PAZ*I_ERI_F2xy_S_G4z_S_vrr+WPZ*I_ERI_F2xy_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2xy_S_F3z_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4z_S_vrr = PAZ*I_ERI_F2xz_S_G4z_S_vrr+WPZ*I_ERI_F2xz_S_G4z_S_M1_vrr+oned2z*I_ERI_D2x_S_G4z_S_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2xz_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4z_S_vrr = PAX*I_ERI_F3y_S_G4z_S_vrr+WPX*I_ERI_F3y_S_G4z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4z_S_vrr = PAZ*I_ERI_Fx2y_S_G4z_S_vrr+WPZ*I_ERI_Fx2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Fx2y_S_F3z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4z_S_vrr = PAY*I_ERI_Fx2z_S_G4z_S_vrr+WPY*I_ERI_Fx2z_S_G4z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4z_S_vrr = PAX*I_ERI_F3z_S_G4z_S_vrr+WPX*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_G4y_S_G4z_S_vrr = PAY*I_ERI_F3y_S_G4z_S_vrr+WPY*I_ERI_F3y_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4z_S_M1_vrr;
      Double I_ERI_G3yz_S_G4z_S_vrr = PAZ*I_ERI_F3y_S_G4z_S_vrr+WPZ*I_ERI_F3y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4z_S_vrr = PAZ*I_ERI_F2yz_S_G4z_S_vrr+WPZ*I_ERI_F2yz_S_G4z_S_M1_vrr+oned2z*I_ERI_D2y_S_G4z_S_vrr-rhod2zsq*I_ERI_D2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2yz_S_F3z_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4z_S_vrr = PAY*I_ERI_F3z_S_G4z_S_vrr+WPY*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_G4z_S_G4z_S_vrr = PAZ*I_ERI_F3z_S_G4z_S_vrr+WPZ*I_ERI_F3z_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      abcd[0] += I_ERI_G4x_S_G4x_S_vrr;
      abcd[1] += I_ERI_G3xy_S_G4x_S_vrr;
      abcd[2] += I_ERI_G3xz_S_G4x_S_vrr;
      abcd[3] += I_ERI_G2x2y_S_G4x_S_vrr;
      abcd[4] += I_ERI_G2xyz_S_G4x_S_vrr;
      abcd[5] += I_ERI_G2x2z_S_G4x_S_vrr;
      abcd[6] += I_ERI_Gx3y_S_G4x_S_vrr;
      abcd[7] += I_ERI_Gx2yz_S_G4x_S_vrr;
      abcd[8] += I_ERI_Gxy2z_S_G4x_S_vrr;
      abcd[9] += I_ERI_Gx3z_S_G4x_S_vrr;
      abcd[10] += I_ERI_G4y_S_G4x_S_vrr;
      abcd[11] += I_ERI_G3yz_S_G4x_S_vrr;
      abcd[12] += I_ERI_G2y2z_S_G4x_S_vrr;
      abcd[13] += I_ERI_Gy3z_S_G4x_S_vrr;
      abcd[14] += I_ERI_G4z_S_G4x_S_vrr;
      abcd[15] += I_ERI_G4x_S_G3xy_S_vrr;
      abcd[16] += I_ERI_G3xy_S_G3xy_S_vrr;
      abcd[17] += I_ERI_G3xz_S_G3xy_S_vrr;
      abcd[18] += I_ERI_G2x2y_S_G3xy_S_vrr;
      abcd[19] += I_ERI_G2xyz_S_G3xy_S_vrr;
      abcd[20] += I_ERI_G2x2z_S_G3xy_S_vrr;
      abcd[21] += I_ERI_Gx3y_S_G3xy_S_vrr;
      abcd[22] += I_ERI_Gx2yz_S_G3xy_S_vrr;
      abcd[23] += I_ERI_Gxy2z_S_G3xy_S_vrr;
      abcd[24] += I_ERI_Gx3z_S_G3xy_S_vrr;
      abcd[25] += I_ERI_G4y_S_G3xy_S_vrr;
      abcd[26] += I_ERI_G3yz_S_G3xy_S_vrr;
      abcd[27] += I_ERI_G2y2z_S_G3xy_S_vrr;
      abcd[28] += I_ERI_Gy3z_S_G3xy_S_vrr;
      abcd[29] += I_ERI_G4z_S_G3xy_S_vrr;
      abcd[30] += I_ERI_G4x_S_G3xz_S_vrr;
      abcd[31] += I_ERI_G3xy_S_G3xz_S_vrr;
      abcd[32] += I_ERI_G3xz_S_G3xz_S_vrr;
      abcd[33] += I_ERI_G2x2y_S_G3xz_S_vrr;
      abcd[34] += I_ERI_G2xyz_S_G3xz_S_vrr;
      abcd[35] += I_ERI_G2x2z_S_G3xz_S_vrr;
      abcd[36] += I_ERI_Gx3y_S_G3xz_S_vrr;
      abcd[37] += I_ERI_Gx2yz_S_G3xz_S_vrr;
      abcd[38] += I_ERI_Gxy2z_S_G3xz_S_vrr;
      abcd[39] += I_ERI_Gx3z_S_G3xz_S_vrr;
      abcd[40] += I_ERI_G4y_S_G3xz_S_vrr;
      abcd[41] += I_ERI_G3yz_S_G3xz_S_vrr;
      abcd[42] += I_ERI_G2y2z_S_G3xz_S_vrr;
      abcd[43] += I_ERI_Gy3z_S_G3xz_S_vrr;
      abcd[44] += I_ERI_G4z_S_G3xz_S_vrr;
      abcd[45] += I_ERI_G4x_S_G2x2y_S_vrr;
      abcd[46] += I_ERI_G3xy_S_G2x2y_S_vrr;
      abcd[47] += I_ERI_G3xz_S_G2x2y_S_vrr;
      abcd[48] += I_ERI_G2x2y_S_G2x2y_S_vrr;
      abcd[49] += I_ERI_G2xyz_S_G2x2y_S_vrr;
      abcd[50] += I_ERI_G2x2z_S_G2x2y_S_vrr;
      abcd[51] += I_ERI_Gx3y_S_G2x2y_S_vrr;
      abcd[52] += I_ERI_Gx2yz_S_G2x2y_S_vrr;
      abcd[53] += I_ERI_Gxy2z_S_G2x2y_S_vrr;
      abcd[54] += I_ERI_Gx3z_S_G2x2y_S_vrr;
      abcd[55] += I_ERI_G4y_S_G2x2y_S_vrr;
      abcd[56] += I_ERI_G3yz_S_G2x2y_S_vrr;
      abcd[57] += I_ERI_G2y2z_S_G2x2y_S_vrr;
      abcd[58] += I_ERI_Gy3z_S_G2x2y_S_vrr;
      abcd[59] += I_ERI_G4z_S_G2x2y_S_vrr;
      abcd[60] += I_ERI_G4x_S_G2xyz_S_vrr;
      abcd[61] += I_ERI_G3xy_S_G2xyz_S_vrr;
      abcd[62] += I_ERI_G3xz_S_G2xyz_S_vrr;
      abcd[63] += I_ERI_G2x2y_S_G2xyz_S_vrr;
      abcd[64] += I_ERI_G2xyz_S_G2xyz_S_vrr;
      abcd[65] += I_ERI_G2x2z_S_G2xyz_S_vrr;
      abcd[66] += I_ERI_Gx3y_S_G2xyz_S_vrr;
      abcd[67] += I_ERI_Gx2yz_S_G2xyz_S_vrr;
      abcd[68] += I_ERI_Gxy2z_S_G2xyz_S_vrr;
      abcd[69] += I_ERI_Gx3z_S_G2xyz_S_vrr;
      abcd[70] += I_ERI_G4y_S_G2xyz_S_vrr;
      abcd[71] += I_ERI_G3yz_S_G2xyz_S_vrr;
      abcd[72] += I_ERI_G2y2z_S_G2xyz_S_vrr;
      abcd[73] += I_ERI_Gy3z_S_G2xyz_S_vrr;
      abcd[74] += I_ERI_G4z_S_G2xyz_S_vrr;
      abcd[75] += I_ERI_G4x_S_G2x2z_S_vrr;
      abcd[76] += I_ERI_G3xy_S_G2x2z_S_vrr;
      abcd[77] += I_ERI_G3xz_S_G2x2z_S_vrr;
      abcd[78] += I_ERI_G2x2y_S_G2x2z_S_vrr;
      abcd[79] += I_ERI_G2xyz_S_G2x2z_S_vrr;
      abcd[80] += I_ERI_G2x2z_S_G2x2z_S_vrr;
      abcd[81] += I_ERI_Gx3y_S_G2x2z_S_vrr;
      abcd[82] += I_ERI_Gx2yz_S_G2x2z_S_vrr;
      abcd[83] += I_ERI_Gxy2z_S_G2x2z_S_vrr;
      abcd[84] += I_ERI_Gx3z_S_G2x2z_S_vrr;
      abcd[85] += I_ERI_G4y_S_G2x2z_S_vrr;
      abcd[86] += I_ERI_G3yz_S_G2x2z_S_vrr;
      abcd[87] += I_ERI_G2y2z_S_G2x2z_S_vrr;
      abcd[88] += I_ERI_Gy3z_S_G2x2z_S_vrr;
      abcd[89] += I_ERI_G4z_S_G2x2z_S_vrr;
      abcd[90] += I_ERI_G4x_S_Gx3y_S_vrr;
      abcd[91] += I_ERI_G3xy_S_Gx3y_S_vrr;
      abcd[92] += I_ERI_G3xz_S_Gx3y_S_vrr;
      abcd[93] += I_ERI_G2x2y_S_Gx3y_S_vrr;
      abcd[94] += I_ERI_G2xyz_S_Gx3y_S_vrr;
      abcd[95] += I_ERI_G2x2z_S_Gx3y_S_vrr;
      abcd[96] += I_ERI_Gx3y_S_Gx3y_S_vrr;
      abcd[97] += I_ERI_Gx2yz_S_Gx3y_S_vrr;
      abcd[98] += I_ERI_Gxy2z_S_Gx3y_S_vrr;
      abcd[99] += I_ERI_Gx3z_S_Gx3y_S_vrr;
      abcd[100] += I_ERI_G4y_S_Gx3y_S_vrr;
      abcd[101] += I_ERI_G3yz_S_Gx3y_S_vrr;
      abcd[102] += I_ERI_G2y2z_S_Gx3y_S_vrr;
      abcd[103] += I_ERI_Gy3z_S_Gx3y_S_vrr;
      abcd[104] += I_ERI_G4z_S_Gx3y_S_vrr;
      abcd[105] += I_ERI_G4x_S_Gx2yz_S_vrr;
      abcd[106] += I_ERI_G3xy_S_Gx2yz_S_vrr;
      abcd[107] += I_ERI_G3xz_S_Gx2yz_S_vrr;
      abcd[108] += I_ERI_G2x2y_S_Gx2yz_S_vrr;
      abcd[109] += I_ERI_G2xyz_S_Gx2yz_S_vrr;
      abcd[110] += I_ERI_G2x2z_S_Gx2yz_S_vrr;
      abcd[111] += I_ERI_Gx3y_S_Gx2yz_S_vrr;
      abcd[112] += I_ERI_Gx2yz_S_Gx2yz_S_vrr;
      abcd[113] += I_ERI_Gxy2z_S_Gx2yz_S_vrr;
      abcd[114] += I_ERI_Gx3z_S_Gx2yz_S_vrr;
      abcd[115] += I_ERI_G4y_S_Gx2yz_S_vrr;
      abcd[116] += I_ERI_G3yz_S_Gx2yz_S_vrr;
      abcd[117] += I_ERI_G2y2z_S_Gx2yz_S_vrr;
      abcd[118] += I_ERI_Gy3z_S_Gx2yz_S_vrr;
      abcd[119] += I_ERI_G4z_S_Gx2yz_S_vrr;
      abcd[120] += I_ERI_G4x_S_Gxy2z_S_vrr;
      abcd[121] += I_ERI_G3xy_S_Gxy2z_S_vrr;
      abcd[122] += I_ERI_G3xz_S_Gxy2z_S_vrr;
      abcd[123] += I_ERI_G2x2y_S_Gxy2z_S_vrr;
      abcd[124] += I_ERI_G2xyz_S_Gxy2z_S_vrr;
      abcd[125] += I_ERI_G2x2z_S_Gxy2z_S_vrr;
      abcd[126] += I_ERI_Gx3y_S_Gxy2z_S_vrr;
      abcd[127] += I_ERI_Gx2yz_S_Gxy2z_S_vrr;
      abcd[128] += I_ERI_Gxy2z_S_Gxy2z_S_vrr;
      abcd[129] += I_ERI_Gx3z_S_Gxy2z_S_vrr;
      abcd[130] += I_ERI_G4y_S_Gxy2z_S_vrr;
      abcd[131] += I_ERI_G3yz_S_Gxy2z_S_vrr;
      abcd[132] += I_ERI_G2y2z_S_Gxy2z_S_vrr;
      abcd[133] += I_ERI_Gy3z_S_Gxy2z_S_vrr;
      abcd[134] += I_ERI_G4z_S_Gxy2z_S_vrr;
      abcd[135] += I_ERI_G4x_S_Gx3z_S_vrr;
      abcd[136] += I_ERI_G3xy_S_Gx3z_S_vrr;
      abcd[137] += I_ERI_G3xz_S_Gx3z_S_vrr;
      abcd[138] += I_ERI_G2x2y_S_Gx3z_S_vrr;
      abcd[139] += I_ERI_G2xyz_S_Gx3z_S_vrr;
      abcd[140] += I_ERI_G2x2z_S_Gx3z_S_vrr;
      abcd[141] += I_ERI_Gx3y_S_Gx3z_S_vrr;
      abcd[142] += I_ERI_Gx2yz_S_Gx3z_S_vrr;
      abcd[143] += I_ERI_Gxy2z_S_Gx3z_S_vrr;
      abcd[144] += I_ERI_Gx3z_S_Gx3z_S_vrr;
      abcd[145] += I_ERI_G4y_S_Gx3z_S_vrr;
      abcd[146] += I_ERI_G3yz_S_Gx3z_S_vrr;
      abcd[147] += I_ERI_G2y2z_S_Gx3z_S_vrr;
      abcd[148] += I_ERI_Gy3z_S_Gx3z_S_vrr;
      abcd[149] += I_ERI_G4z_S_Gx3z_S_vrr;
      abcd[150] += I_ERI_G4x_S_G4y_S_vrr;
      abcd[151] += I_ERI_G3xy_S_G4y_S_vrr;
      abcd[152] += I_ERI_G3xz_S_G4y_S_vrr;
      abcd[153] += I_ERI_G2x2y_S_G4y_S_vrr;
      abcd[154] += I_ERI_G2xyz_S_G4y_S_vrr;
      abcd[155] += I_ERI_G2x2z_S_G4y_S_vrr;
      abcd[156] += I_ERI_Gx3y_S_G4y_S_vrr;
      abcd[157] += I_ERI_Gx2yz_S_G4y_S_vrr;
      abcd[158] += I_ERI_Gxy2z_S_G4y_S_vrr;
      abcd[159] += I_ERI_Gx3z_S_G4y_S_vrr;
      abcd[160] += I_ERI_G4y_S_G4y_S_vrr;
      abcd[161] += I_ERI_G3yz_S_G4y_S_vrr;
      abcd[162] += I_ERI_G2y2z_S_G4y_S_vrr;
      abcd[163] += I_ERI_Gy3z_S_G4y_S_vrr;
      abcd[164] += I_ERI_G4z_S_G4y_S_vrr;
      abcd[165] += I_ERI_G4x_S_G3yz_S_vrr;
      abcd[166] += I_ERI_G3xy_S_G3yz_S_vrr;
      abcd[167] += I_ERI_G3xz_S_G3yz_S_vrr;
      abcd[168] += I_ERI_G2x2y_S_G3yz_S_vrr;
      abcd[169] += I_ERI_G2xyz_S_G3yz_S_vrr;
      abcd[170] += I_ERI_G2x2z_S_G3yz_S_vrr;
      abcd[171] += I_ERI_Gx3y_S_G3yz_S_vrr;
      abcd[172] += I_ERI_Gx2yz_S_G3yz_S_vrr;
      abcd[173] += I_ERI_Gxy2z_S_G3yz_S_vrr;
      abcd[174] += I_ERI_Gx3z_S_G3yz_S_vrr;
      abcd[175] += I_ERI_G4y_S_G3yz_S_vrr;
      abcd[176] += I_ERI_G3yz_S_G3yz_S_vrr;
      abcd[177] += I_ERI_G2y2z_S_G3yz_S_vrr;
      abcd[178] += I_ERI_Gy3z_S_G3yz_S_vrr;
      abcd[179] += I_ERI_G4z_S_G3yz_S_vrr;
      abcd[180] += I_ERI_G4x_S_G2y2z_S_vrr;
      abcd[181] += I_ERI_G3xy_S_G2y2z_S_vrr;
      abcd[182] += I_ERI_G3xz_S_G2y2z_S_vrr;
      abcd[183] += I_ERI_G2x2y_S_G2y2z_S_vrr;
      abcd[184] += I_ERI_G2xyz_S_G2y2z_S_vrr;
      abcd[185] += I_ERI_G2x2z_S_G2y2z_S_vrr;
      abcd[186] += I_ERI_Gx3y_S_G2y2z_S_vrr;
      abcd[187] += I_ERI_Gx2yz_S_G2y2z_S_vrr;
      abcd[188] += I_ERI_Gxy2z_S_G2y2z_S_vrr;
      abcd[189] += I_ERI_Gx3z_S_G2y2z_S_vrr;
      abcd[190] += I_ERI_G4y_S_G2y2z_S_vrr;
      abcd[191] += I_ERI_G3yz_S_G2y2z_S_vrr;
      abcd[192] += I_ERI_G2y2z_S_G2y2z_S_vrr;
      abcd[193] += I_ERI_Gy3z_S_G2y2z_S_vrr;
      abcd[194] += I_ERI_G4z_S_G2y2z_S_vrr;
      abcd[195] += I_ERI_G4x_S_Gy3z_S_vrr;
      abcd[196] += I_ERI_G3xy_S_Gy3z_S_vrr;
      abcd[197] += I_ERI_G3xz_S_Gy3z_S_vrr;
      abcd[198] += I_ERI_G2x2y_S_Gy3z_S_vrr;
      abcd[199] += I_ERI_G2xyz_S_Gy3z_S_vrr;
      abcd[200] += I_ERI_G2x2z_S_Gy3z_S_vrr;
      abcd[201] += I_ERI_Gx3y_S_Gy3z_S_vrr;
      abcd[202] += I_ERI_Gx2yz_S_Gy3z_S_vrr;
      abcd[203] += I_ERI_Gxy2z_S_Gy3z_S_vrr;
      abcd[204] += I_ERI_Gx3z_S_Gy3z_S_vrr;
      abcd[205] += I_ERI_G4y_S_Gy3z_S_vrr;
      abcd[206] += I_ERI_G3yz_S_Gy3z_S_vrr;
      abcd[207] += I_ERI_G2y2z_S_Gy3z_S_vrr;
      abcd[208] += I_ERI_Gy3z_S_Gy3z_S_vrr;
      abcd[209] += I_ERI_G4z_S_Gy3z_S_vrr;
      abcd[210] += I_ERI_G4x_S_G4z_S_vrr;
      abcd[211] += I_ERI_G3xy_S_G4z_S_vrr;
      abcd[212] += I_ERI_G3xz_S_G4z_S_vrr;
      abcd[213] += I_ERI_G2x2y_S_G4z_S_vrr;
      abcd[214] += I_ERI_G2xyz_S_G4z_S_vrr;
      abcd[215] += I_ERI_G2x2z_S_G4z_S_vrr;
      abcd[216] += I_ERI_Gx3y_S_G4z_S_vrr;
      abcd[217] += I_ERI_Gx2yz_S_G4z_S_vrr;
      abcd[218] += I_ERI_Gxy2z_S_G4z_S_vrr;
      abcd[219] += I_ERI_Gx3z_S_G4z_S_vrr;
      abcd[220] += I_ERI_G4y_S_G4z_S_vrr;
      abcd[221] += I_ERI_G3yz_S_G4z_S_vrr;
      abcd[222] += I_ERI_G2y2z_S_G4z_S_vrr;
      abcd[223] += I_ERI_Gy3z_S_G4z_S_vrr;
      abcd[224] += I_ERI_G4z_S_G4z_S_vrr;
    }
  }
}
