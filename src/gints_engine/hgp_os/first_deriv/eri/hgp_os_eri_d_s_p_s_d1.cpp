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

//
// print out the information regarding of derivatives 
// here below we count on all of RHS integrals(including the repeat ones)
// this is used to simulate the FLOPS counting
// BRA1 as redundant position, total RHS integrals evaluated as: 12330
// BRA2 as redundant position, total RHS integrals evaluated as: 12249
// KET1 as redundant position, total RHS integrals evaluated as: 12339
// KET2 as redundant position, total RHS integrals evaluated as: 12249
// the redundant position is: KET2
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####
// BRA2
// X
// Y
// Z
// ####
// KET1
// X
// Y
// Z
// ####

void hgp_os_eri_d_s_p_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_Px_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_a = 0.0E0;
  Double I_ERI_Px_S_Px_S = 0.0E0;
  Double I_ERI_Py_S_Px_S = 0.0E0;
  Double I_ERI_Pz_S_Px_S = 0.0E0;
  Double I_ERI_Px_S_Py_S = 0.0E0;
  Double I_ERI_Py_S_Py_S = 0.0E0;
  Double I_ERI_Pz_S_Py_S = 0.0E0;
  Double I_ERI_Px_S_Pz_S = 0.0E0;
  Double I_ERI_Py_S_Pz_S = 0.0E0;
  Double I_ERI_Pz_S_Pz_S = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_c = 0.0E0;
  Double I_ERI_D2x_S_S_S = 0.0E0;
  Double I_ERI_Dxy_S_S_S = 0.0E0;
  Double I_ERI_Dxz_S_S_S = 0.0E0;
  Double I_ERI_D2y_S_S_S = 0.0E0;
  Double I_ERI_Dyz_S_S_S = 0.0E0;
  Double I_ERI_D2z_S_S_S = 0.0E0;
  Double I_ERI_F3x_S_Px_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_b = 0.0E0;

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
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
      Double eta   = 1.0E0/onede;
      Double ediff = jexpdiff[jp2];
      Double gamma = 0.5E0*(eta+ediff);
      Double delta = 0.5E0*(eta-ediff);
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
      isSignificant = true;


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

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER43;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER17*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER15*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER13*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER11*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = ONEOVER9*I_ERI_S_S_S_S_M4_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M4_vrr  = f*I_ERI_S_S_S_S_M4_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);

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
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_Px_S_S_S_M3_vrr = PAX*I_ERI_S_S_S_S_M3_vrr+WPX*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Py_S_S_S_M3_vrr = PAY*I_ERI_S_S_S_S_M3_vrr+WPY*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Pz_S_S_S_M3_vrr = PAZ*I_ERI_S_S_S_S_M3_vrr+WPZ*I_ERI_S_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_Px_S_S_S_M2_vrr = PAX*I_ERI_S_S_S_S_M2_vrr+WPX*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Py_S_S_S_M2_vrr = PAY*I_ERI_S_S_S_S_M2_vrr+WPY*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Pz_S_S_S_M2_vrr = PAZ*I_ERI_S_S_S_S_M2_vrr+WPZ*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M2_vrr = PAX*I_ERI_Px_S_S_S_M2_vrr+WPX*I_ERI_Px_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_S_S_M2_vrr = PAY*I_ERI_Px_S_S_S_M2_vrr+WPY*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_S_S_M2_vrr = PAZ*I_ERI_Px_S_S_S_M2_vrr+WPZ*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_S_S_M2_vrr = PAY*I_ERI_Py_S_S_S_M2_vrr+WPY*I_ERI_Py_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_S_S_M2_vrr = PAZ*I_ERI_Py_S_S_S_M2_vrr+WPZ*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_S_S_M2_vrr = PAZ*I_ERI_Pz_S_S_S_M2_vrr+WPZ*I_ERI_Pz_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_Px_S_S_S_M1_vrr = PAX*I_ERI_S_S_S_S_M1_vrr+WPX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Py_S_S_S_M1_vrr = PAY*I_ERI_S_S_S_S_M1_vrr+WPY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_S_S_M1_vrr = PAZ*I_ERI_S_S_S_S_M1_vrr+WPZ*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M1_vrr = PAX*I_ERI_S_S_Px_S_M1_vrr+WPX*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Py_S_Px_S_M1_vrr = PAY*I_ERI_S_S_Px_S_M1_vrr+WPY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Px_S_M1_vrr = PAZ*I_ERI_S_S_Px_S_M1_vrr+WPZ*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_Py_S_M1_vrr = PAX*I_ERI_S_S_Py_S_M1_vrr+WPX*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Py_S_M1_vrr = PAY*I_ERI_S_S_Py_S_M1_vrr+WPY*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_Py_S_M1_vrr = PAZ*I_ERI_S_S_Py_S_M1_vrr+WPZ*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_Pz_S_M1_vrr = PAX*I_ERI_S_S_Pz_S_M1_vrr+WPX*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Py_S_Pz_S_M1_vrr = PAY*I_ERI_S_S_Pz_S_M1_vrr+WPY*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Pz_S_Pz_S_M1_vrr = PAZ*I_ERI_S_S_Pz_S_M1_vrr+WPZ*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M1_vrr = PAX*I_ERI_Px_S_S_S_M1_vrr+WPX*I_ERI_Px_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_S_S_M1_vrr = PAY*I_ERI_Px_S_S_S_M1_vrr+WPY*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_S_S_M1_vrr = PAZ*I_ERI_Px_S_S_S_M1_vrr+WPZ*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_S_S_M1_vrr = PAY*I_ERI_Py_S_S_S_M1_vrr+WPY*I_ERI_Py_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_S_S_M1_vrr = PAZ*I_ERI_Py_S_S_S_M1_vrr+WPZ*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_S_S_M1_vrr = PAZ*I_ERI_Pz_S_S_S_M1_vrr+WPZ*I_ERI_Pz_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M1_vrr = QCX*I_ERI_D2x_S_S_S_M1_vrr+WQX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Px_S_M1_vrr = QCX*I_ERI_Dxy_S_S_S_M1_vrr+WQX*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Px_S_M1_vrr = QCX*I_ERI_Dxz_S_S_S_M1_vrr+WQX*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Px_S_M1_vrr = QCX*I_ERI_D2y_S_S_S_M1_vrr+WQX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Px_S_M1_vrr = QCX*I_ERI_Dyz_S_S_S_M1_vrr+WQX*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Px_S_M1_vrr = QCX*I_ERI_D2z_S_S_S_M1_vrr+WQX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Py_S_M1_vrr = QCY*I_ERI_D2x_S_S_S_M1_vrr+WQY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Py_S_M1_vrr = QCY*I_ERI_Dxy_S_S_S_M1_vrr+WQY*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Py_S_M1_vrr = QCY*I_ERI_Dxz_S_S_S_M1_vrr+WQY*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Py_S_M1_vrr = QCY*I_ERI_D2y_S_S_S_M1_vrr+WQY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Py_S_M1_vrr = QCY*I_ERI_Dyz_S_S_S_M1_vrr+WQY*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Py_S_M1_vrr = QCY*I_ERI_D2z_S_S_S_M1_vrr+WQY*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Pz_S_M1_vrr = QCZ*I_ERI_D2x_S_S_S_M1_vrr+WQZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Pz_S_M1_vrr = QCZ*I_ERI_Dxy_S_S_S_M1_vrr+WQZ*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Pz_S_M1_vrr = QCZ*I_ERI_Dxz_S_S_S_M1_vrr+WQZ*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Pz_S_M1_vrr = QCZ*I_ERI_D2y_S_S_S_M1_vrr+WQZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Pz_S_M1_vrr = QCZ*I_ERI_Dyz_S_S_S_M1_vrr+WQZ*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Pz_S_M1_vrr = QCZ*I_ERI_D2z_S_S_S_M1_vrr+WQZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Pz_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_Px_S_S_S_vrr = PAX*I_ERI_S_S_S_S_vrr+WPX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Py_S_S_S_vrr = PAY*I_ERI_S_S_S_S_vrr+WPY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_S_S_vrr = PAZ*I_ERI_S_S_S_S_vrr+WPZ*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_Px_S_Px_S_vrr = PAX*I_ERI_S_S_Px_S_vrr+WPX*I_ERI_S_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Py_S_Px_S_vrr = PAY*I_ERI_S_S_Px_S_vrr+WPY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Px_S_vrr = PAZ*I_ERI_S_S_Px_S_vrr+WPZ*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_Py_S_vrr = PAX*I_ERI_S_S_Py_S_vrr+WPX*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_Py_S_vrr = PAY*I_ERI_S_S_Py_S_vrr+WPY*I_ERI_S_S_Py_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_Py_S_vrr = PAZ*I_ERI_S_S_Py_S_vrr+WPZ*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_Pz_S_vrr = PAX*I_ERI_S_S_Pz_S_vrr+WPX*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Py_S_Pz_S_vrr = PAY*I_ERI_S_S_Pz_S_vrr+WPY*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Pz_S_Pz_S_vrr = PAZ*I_ERI_S_S_Pz_S_vrr+WPZ*I_ERI_S_S_Pz_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_S_vrr = PAY*I_ERI_Px_S_S_S_vrr+WPY*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_S_S_vrr = PAZ*I_ERI_Px_S_S_S_vrr+WPZ*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_S_S_vrr = PAZ*I_ERI_Py_S_S_S_vrr+WPZ*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_S_vrr = PAZ*I_ERI_Pz_S_S_S_vrr+WPZ*I_ERI_Pz_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_vrr = QCX*I_ERI_D2x_S_S_S_vrr+WQX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Px_S_vrr = QCX*I_ERI_Dxy_S_S_S_vrr+WQX*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Px_S_vrr = QCX*I_ERI_Dxz_S_S_S_vrr+WQX*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Px_S_vrr = QCX*I_ERI_D2y_S_S_S_vrr+WQX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Px_S_vrr = QCX*I_ERI_Dyz_S_S_S_vrr+WQX*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Px_S_vrr = QCX*I_ERI_D2z_S_S_S_vrr+WQX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Py_S_vrr = QCY*I_ERI_D2x_S_S_S_vrr+WQY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Py_S_vrr = QCY*I_ERI_Dxy_S_S_S_vrr+WQY*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Py_S_vrr = QCY*I_ERI_Dxz_S_S_S_vrr+WQY*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Py_S_vrr = QCY*I_ERI_D2y_S_S_S_vrr+WQY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Py_S_vrr = QCY*I_ERI_Dyz_S_S_S_vrr+WQY*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Py_S_vrr = QCY*I_ERI_D2z_S_S_S_vrr+WQY*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Pz_S_vrr = QCZ*I_ERI_D2x_S_S_S_vrr+WQZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Pz_S_vrr = QCZ*I_ERI_Dxy_S_S_S_vrr+WQZ*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Pz_S_vrr = QCZ*I_ERI_Dxz_S_S_S_vrr+WQZ*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Pz_S_vrr = QCZ*I_ERI_D2y_S_S_S_vrr+WQZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Pz_S_vrr = QCZ*I_ERI_Dyz_S_S_S_vrr+WQZ*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Pz_S_vrr = QCZ*I_ERI_D2z_S_S_S_vrr+WQZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_vrr = QCX*I_ERI_D2x_S_Px_S_vrr+WQX*I_ERI_D2x_S_Px_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_D2x_S_vrr = QCX*I_ERI_Dxy_S_Px_S_vrr+WQX*I_ERI_Dxy_S_Px_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_D2x_S_vrr = QCX*I_ERI_Dxz_S_Px_S_vrr+WQX*I_ERI_Dxz_S_Px_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_D2x_S_vrr = QCX*I_ERI_D2y_S_Px_S_vrr+WQX*I_ERI_D2y_S_Px_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_D2x_S_vrr = QCX*I_ERI_Dyz_S_Px_S_vrr+WQX*I_ERI_Dyz_S_Px_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_D2x_S_vrr = QCX*I_ERI_D2z_S_Px_S_vrr+WQX*I_ERI_D2z_S_Px_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dxy_S_vrr = QCY*I_ERI_D2x_S_Px_S_vrr+WQY*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxy_S_vrr = QCY*I_ERI_Dxy_S_Px_S_vrr+WQY*I_ERI_Dxy_S_Px_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxy_S_vrr = QCY*I_ERI_Dxz_S_Px_S_vrr+WQY*I_ERI_Dxz_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxy_S_vrr = QCY*I_ERI_D2y_S_Px_S_vrr+WQY*I_ERI_D2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxy_S_vrr = QCY*I_ERI_Dyz_S_Px_S_vrr+WQY*I_ERI_Dyz_S_Px_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxy_S_vrr = QCY*I_ERI_D2z_S_Px_S_vrr+WQY*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_Dxz_S_vrr = QCZ*I_ERI_D2x_S_Px_S_vrr+WQZ*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxz_S_vrr = QCZ*I_ERI_Dxy_S_Px_S_vrr+WQZ*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxz_S_vrr = QCZ*I_ERI_Dxz_S_Px_S_vrr+WQZ*I_ERI_Dxz_S_Px_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxz_S_vrr = QCZ*I_ERI_D2y_S_Px_S_vrr+WQZ*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxz_S_vrr = QCZ*I_ERI_Dyz_S_Px_S_vrr+WQZ*I_ERI_Dyz_S_Px_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxz_S_vrr = QCZ*I_ERI_D2z_S_Px_S_vrr+WQZ*I_ERI_D2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_D2y_S_vrr = QCY*I_ERI_D2x_S_Py_S_vrr+WQY*I_ERI_D2x_S_Py_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_D2y_S_vrr = QCY*I_ERI_Dxy_S_Py_S_vrr+WQY*I_ERI_Dxy_S_Py_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_D2y_S_vrr = QCY*I_ERI_Dxz_S_Py_S_vrr+WQY*I_ERI_Dxz_S_Py_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = QCY*I_ERI_D2y_S_Py_S_vrr+WQY*I_ERI_D2y_S_Py_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_D2y_S_vrr = QCY*I_ERI_Dyz_S_Py_S_vrr+WQY*I_ERI_Dyz_S_Py_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = QCY*I_ERI_D2z_S_Py_S_vrr+WQY*I_ERI_D2z_S_Py_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dyz_S_vrr = QCZ*I_ERI_D2x_S_Py_S_vrr+WQZ*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_Dyz_S_vrr = QCZ*I_ERI_Dxy_S_Py_S_vrr+WQZ*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_Dyz_S_vrr = QCZ*I_ERI_Dxz_S_Py_S_vrr+WQZ*I_ERI_Dxz_S_Py_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Dyz_S_vrr = QCZ*I_ERI_D2y_S_Py_S_vrr+WQZ*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_Dyz_S_vrr = QCZ*I_ERI_Dyz_S_Py_S_vrr+WQZ*I_ERI_Dyz_S_Py_S_M1_vrr+oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Dyz_S_vrr = QCZ*I_ERI_D2z_S_Py_S_vrr+WQZ*I_ERI_D2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = QCZ*I_ERI_D2x_S_Pz_S_vrr+WQZ*I_ERI_D2x_S_Pz_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_D2z_S_vrr = QCZ*I_ERI_Dxy_S_Pz_S_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_D2z_S_vrr = QCZ*I_ERI_Dxz_S_Pz_S_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = QCZ*I_ERI_D2y_S_Pz_S_vrr+WQZ*I_ERI_D2y_S_Pz_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_D2z_S_vrr = QCZ*I_ERI_Dyz_S_Pz_S_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = QCZ*I_ERI_D2z_S_Pz_S_vrr+WQZ*I_ERI_D2z_S_Pz_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_vrr = PAX*I_ERI_D2x_S_Px_S_vrr+WPX*I_ERI_D2x_S_Px_S_M1_vrr+2*oned2z*I_ERI_Px_S_Px_S_vrr-2*rhod2zsq*I_ERI_Px_S_Px_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_Px_S_vrr = PAY*I_ERI_D2x_S_Px_S_vrr+WPY*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Px_S_vrr = PAZ*I_ERI_D2x_S_Px_S_vrr+WPZ*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Px_S_vrr = PAX*I_ERI_D2y_S_Px_S_vrr+WPX*I_ERI_D2y_S_Px_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_Px_S_vrr = PAZ*I_ERI_Dxy_S_Px_S_vrr+WPZ*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Px_S_vrr = PAX*I_ERI_D2z_S_Px_S_vrr+WPX*I_ERI_D2z_S_Px_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_Px_S_vrr = PAY*I_ERI_D2y_S_Px_S_vrr+WPY*I_ERI_D2y_S_Px_S_M1_vrr+2*oned2z*I_ERI_Py_S_Px_S_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Px_S_vrr = PAZ*I_ERI_D2y_S_Px_S_vrr+WPZ*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Px_S_vrr = PAY*I_ERI_D2z_S_Px_S_vrr+WPY*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Px_S_vrr = PAZ*I_ERI_D2z_S_Px_S_vrr+WPZ*I_ERI_D2z_S_Px_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Px_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_Py_S_vrr = PAX*I_ERI_D2x_S_Py_S_vrr+WPX*I_ERI_D2x_S_Py_S_M1_vrr+2*oned2z*I_ERI_Px_S_Py_S_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_Py_S_vrr = PAY*I_ERI_D2x_S_Py_S_vrr+WPY*I_ERI_D2x_S_Py_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_Py_S_vrr = PAZ*I_ERI_D2x_S_Py_S_vrr+WPZ*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_Py_S_vrr = PAX*I_ERI_D2y_S_Py_S_vrr+WPX*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_Py_S_vrr = PAZ*I_ERI_Dxy_S_Py_S_vrr+WPZ*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_Py_S_vrr = PAX*I_ERI_D2z_S_Py_S_vrr+WPX*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_Py_S_vrr = PAY*I_ERI_D2y_S_Py_S_vrr+WPY*I_ERI_D2y_S_Py_S_M1_vrr+2*oned2z*I_ERI_Py_S_Py_S_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_Py_S_vrr = PAZ*I_ERI_D2y_S_Py_S_vrr+WPZ*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_Py_S_vrr = PAY*I_ERI_D2z_S_Py_S_vrr+WPY*I_ERI_D2z_S_Py_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_Py_S_vrr = PAZ*I_ERI_D2z_S_Py_S_vrr+WPZ*I_ERI_D2z_S_Py_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Py_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_Pz_S_vrr = PAX*I_ERI_D2x_S_Pz_S_vrr+WPX*I_ERI_D2x_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_F2xy_S_Pz_S_vrr = PAY*I_ERI_D2x_S_Pz_S_vrr+WPY*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xz_S_Pz_S_vrr = PAZ*I_ERI_D2x_S_Pz_S_vrr+WPZ*I_ERI_D2x_S_Pz_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_Pz_S_vrr = PAX*I_ERI_D2y_S_Pz_S_vrr+WPX*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Pz_S_vrr = PAZ*I_ERI_Dxy_S_Pz_S_vrr+WPZ*I_ERI_Dxy_S_Pz_S_M1_vrr+oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_Pz_S_vrr = PAX*I_ERI_D2z_S_Pz_S_vrr+WPX*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3y_S_Pz_S_vrr = PAY*I_ERI_D2y_S_Pz_S_vrr+WPY*I_ERI_D2y_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_F2yz_S_Pz_S_vrr = PAZ*I_ERI_D2y_S_Pz_S_vrr+WPZ*I_ERI_D2y_S_Pz_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_Pz_S_vrr = PAY*I_ERI_D2z_S_Pz_S_vrr+WPY*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3z_S_Pz_S_vrr = PAZ*I_ERI_D2z_S_Pz_S_vrr+WPZ*I_ERI_D2z_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_a_coefs = alpha;
      I_ERI_F3x_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_Px_S += I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S += I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S += I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S += I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S += I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S += I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S += I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S += I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S += I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_c_coefs = gamma;
      I_ERI_D2x_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_S_S += I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S += I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S += I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S += I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S += I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S += I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_b_coefs = beta;
      I_ERI_F3x_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_b_coefs = beta;
      I_ERI_D2x_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Pz_S_vrr;
    }
  }

  /************************************************************
   * let's see the significance test result. if VRR result is
   * insignificant, there's no need to do following codes
   ************************************************************/
  if (! isSignificant) return;

  /************************************************************
   * declare the HRR1 result shell quartets in array form
   ************************************************************/

  /************************************************************
   * initilize the HRR steps : build the AB/CD variables
   ************************************************************/
  Double ABX = A[0] - B[0];
  Double ABY = A[1] - B[1];
  Double ABZ = A[2] - B[2];

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_b = I_ERI_F3x_S_Px_S_b+ABX*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Px_Px_S_b = I_ERI_F2xy_S_Px_S_b+ABX*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Px_Px_S_b = I_ERI_F2xz_S_Px_S_b+ABX*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Px_Px_S_b = I_ERI_Fx2y_S_Px_S_b+ABX*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Px_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABX*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Px_Px_S_b = I_ERI_Fx2z_S_Px_S_b+ABX*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Py_Px_S_b = I_ERI_F2xy_S_Px_S_b+ABY*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Py_Px_S_b = I_ERI_Fx2y_S_Px_S_b+ABY*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Py_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABY*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Py_Px_S_b = I_ERI_F3y_S_Px_S_b+ABY*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Py_Px_S_b = I_ERI_F2yz_S_Px_S_b+ABY*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Py_Px_S_b = I_ERI_Fy2z_S_Px_S_b+ABY*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Pz_Px_S_b = I_ERI_F2xz_S_Px_S_b+ABZ*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Pz_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABZ*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Pz_Px_S_b = I_ERI_Fx2z_S_Px_S_b+ABZ*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Pz_Px_S_b = I_ERI_F2yz_S_Px_S_b+ABZ*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Pz_Px_S_b = I_ERI_Fy2z_S_Px_S_b+ABZ*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Pz_Px_S_b = I_ERI_F3z_S_Px_S_b+ABZ*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Px_Py_S_b = I_ERI_F3x_S_Py_S_b+ABX*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Px_Py_S_b = I_ERI_F2xy_S_Py_S_b+ABX*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Px_Py_S_b = I_ERI_F2xz_S_Py_S_b+ABX*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Px_Py_S_b = I_ERI_Fx2y_S_Py_S_b+ABX*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Px_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABX*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Px_Py_S_b = I_ERI_Fx2z_S_Py_S_b+ABX*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Py_Py_S_b = I_ERI_F2xy_S_Py_S_b+ABY*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Py_Py_S_b = I_ERI_Fx2y_S_Py_S_b+ABY*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Py_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABY*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Py_Py_S_b = I_ERI_F3y_S_Py_S_b+ABY*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Py_Py_S_b = I_ERI_F2yz_S_Py_S_b+ABY*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Py_Py_S_b = I_ERI_Fy2z_S_Py_S_b+ABY*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Pz_Py_S_b = I_ERI_F2xz_S_Py_S_b+ABZ*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Pz_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABZ*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Pz_Py_S_b = I_ERI_Fx2z_S_Py_S_b+ABZ*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Pz_Py_S_b = I_ERI_F2yz_S_Py_S_b+ABZ*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Pz_Py_S_b = I_ERI_Fy2z_S_Py_S_b+ABZ*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Pz_Py_S_b = I_ERI_F3z_S_Py_S_b+ABZ*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Px_Pz_S_b = I_ERI_F3x_S_Pz_S_b+ABX*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Px_Pz_S_b = I_ERI_F2xy_S_Pz_S_b+ABX*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Px_Pz_S_b = I_ERI_F2xz_S_Pz_S_b+ABX*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Px_Pz_S_b = I_ERI_Fx2y_S_Pz_S_b+ABX*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Px_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABX*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Px_Pz_S_b = I_ERI_Fx2z_S_Pz_S_b+ABX*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_Py_Pz_S_b = I_ERI_F2xy_S_Pz_S_b+ABY*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Py_Pz_S_b = I_ERI_Fx2y_S_Pz_S_b+ABY*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Py_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABY*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Py_Pz_S_b = I_ERI_F3y_S_Pz_S_b+ABY*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Py_Pz_S_b = I_ERI_F2yz_S_Pz_S_b+ABY*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Py_Pz_S_b = I_ERI_Fy2z_S_Pz_S_b+ABY*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_Pz_Pz_S_b = I_ERI_F2xz_S_Pz_S_b+ABZ*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Pz_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABZ*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Pz_Pz_S_b = I_ERI_Fx2z_S_Pz_S_b+ABZ*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Pz_Pz_S_b = I_ERI_F2yz_S_Pz_S_b+ABZ*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Pz_Pz_S_b = I_ERI_Fy2z_S_Pz_S_b+ABZ*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Pz_Pz_S_b = I_ERI_F3z_S_Pz_S_b+ABZ*I_ERI_D2z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_F3x_S_Px_S_a-2*I_ERI_Px_S_Px_S;
  abcd[1] = 2.0E0*I_ERI_F2xy_S_Px_S_a-1*I_ERI_Py_S_Px_S;
  abcd[2] = 2.0E0*I_ERI_F2xz_S_Px_S_a-1*I_ERI_Pz_S_Px_S;
  abcd[3] = 2.0E0*I_ERI_Fx2y_S_Px_S_a;
  abcd[4] = 2.0E0*I_ERI_Fxyz_S_Px_S_a;
  abcd[5] = 2.0E0*I_ERI_Fx2z_S_Px_S_a;
  abcd[6] = 2.0E0*I_ERI_F3x_S_Py_S_a-2*I_ERI_Px_S_Py_S;
  abcd[7] = 2.0E0*I_ERI_F2xy_S_Py_S_a-1*I_ERI_Py_S_Py_S;
  abcd[8] = 2.0E0*I_ERI_F2xz_S_Py_S_a-1*I_ERI_Pz_S_Py_S;
  abcd[9] = 2.0E0*I_ERI_Fx2y_S_Py_S_a;
  abcd[10] = 2.0E0*I_ERI_Fxyz_S_Py_S_a;
  abcd[11] = 2.0E0*I_ERI_Fx2z_S_Py_S_a;
  abcd[12] = 2.0E0*I_ERI_F3x_S_Pz_S_a-2*I_ERI_Px_S_Pz_S;
  abcd[13] = 2.0E0*I_ERI_F2xy_S_Pz_S_a-1*I_ERI_Py_S_Pz_S;
  abcd[14] = 2.0E0*I_ERI_F2xz_S_Pz_S_a-1*I_ERI_Pz_S_Pz_S;
  abcd[15] = 2.0E0*I_ERI_Fx2y_S_Pz_S_a;
  abcd[16] = 2.0E0*I_ERI_Fxyz_S_Pz_S_a;
  abcd[17] = 2.0E0*I_ERI_Fx2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[18] = 2.0E0*I_ERI_F2xy_S_Px_S_a;
  abcd[19] = 2.0E0*I_ERI_Fx2y_S_Px_S_a-1*I_ERI_Px_S_Px_S;
  abcd[20] = 2.0E0*I_ERI_Fxyz_S_Px_S_a;
  abcd[21] = 2.0E0*I_ERI_F3y_S_Px_S_a-2*I_ERI_Py_S_Px_S;
  abcd[22] = 2.0E0*I_ERI_F2yz_S_Px_S_a-1*I_ERI_Pz_S_Px_S;
  abcd[23] = 2.0E0*I_ERI_Fy2z_S_Px_S_a;
  abcd[24] = 2.0E0*I_ERI_F2xy_S_Py_S_a;
  abcd[25] = 2.0E0*I_ERI_Fx2y_S_Py_S_a-1*I_ERI_Px_S_Py_S;
  abcd[26] = 2.0E0*I_ERI_Fxyz_S_Py_S_a;
  abcd[27] = 2.0E0*I_ERI_F3y_S_Py_S_a-2*I_ERI_Py_S_Py_S;
  abcd[28] = 2.0E0*I_ERI_F2yz_S_Py_S_a-1*I_ERI_Pz_S_Py_S;
  abcd[29] = 2.0E0*I_ERI_Fy2z_S_Py_S_a;
  abcd[30] = 2.0E0*I_ERI_F2xy_S_Pz_S_a;
  abcd[31] = 2.0E0*I_ERI_Fx2y_S_Pz_S_a-1*I_ERI_Px_S_Pz_S;
  abcd[32] = 2.0E0*I_ERI_Fxyz_S_Pz_S_a;
  abcd[33] = 2.0E0*I_ERI_F3y_S_Pz_S_a-2*I_ERI_Py_S_Pz_S;
  abcd[34] = 2.0E0*I_ERI_F2yz_S_Pz_S_a-1*I_ERI_Pz_S_Pz_S;
  abcd[35] = 2.0E0*I_ERI_Fy2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[36] = 2.0E0*I_ERI_F2xz_S_Px_S_a;
  abcd[37] = 2.0E0*I_ERI_Fxyz_S_Px_S_a;
  abcd[38] = 2.0E0*I_ERI_Fx2z_S_Px_S_a-1*I_ERI_Px_S_Px_S;
  abcd[39] = 2.0E0*I_ERI_F2yz_S_Px_S_a;
  abcd[40] = 2.0E0*I_ERI_Fy2z_S_Px_S_a-1*I_ERI_Py_S_Px_S;
  abcd[41] = 2.0E0*I_ERI_F3z_S_Px_S_a-2*I_ERI_Pz_S_Px_S;
  abcd[42] = 2.0E0*I_ERI_F2xz_S_Py_S_a;
  abcd[43] = 2.0E0*I_ERI_Fxyz_S_Py_S_a;
  abcd[44] = 2.0E0*I_ERI_Fx2z_S_Py_S_a-1*I_ERI_Px_S_Py_S;
  abcd[45] = 2.0E0*I_ERI_F2yz_S_Py_S_a;
  abcd[46] = 2.0E0*I_ERI_Fy2z_S_Py_S_a-1*I_ERI_Py_S_Py_S;
  abcd[47] = 2.0E0*I_ERI_F3z_S_Py_S_a-2*I_ERI_Pz_S_Py_S;
  abcd[48] = 2.0E0*I_ERI_F2xz_S_Pz_S_a;
  abcd[49] = 2.0E0*I_ERI_Fxyz_S_Pz_S_a;
  abcd[50] = 2.0E0*I_ERI_Fx2z_S_Pz_S_a-1*I_ERI_Px_S_Pz_S;
  abcd[51] = 2.0E0*I_ERI_F2yz_S_Pz_S_a;
  abcd[52] = 2.0E0*I_ERI_Fy2z_S_Pz_S_a-1*I_ERI_Py_S_Pz_S;
  abcd[53] = 2.0E0*I_ERI_F3z_S_Pz_S_a-2*I_ERI_Pz_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[54] = 2.0E0*I_ERI_D2x_Px_Px_S_b;
  abcd[55] = 2.0E0*I_ERI_Dxy_Px_Px_S_b;
  abcd[56] = 2.0E0*I_ERI_Dxz_Px_Px_S_b;
  abcd[57] = 2.0E0*I_ERI_D2y_Px_Px_S_b;
  abcd[58] = 2.0E0*I_ERI_Dyz_Px_Px_S_b;
  abcd[59] = 2.0E0*I_ERI_D2z_Px_Px_S_b;
  abcd[60] = 2.0E0*I_ERI_D2x_Px_Py_S_b;
  abcd[61] = 2.0E0*I_ERI_Dxy_Px_Py_S_b;
  abcd[62] = 2.0E0*I_ERI_Dxz_Px_Py_S_b;
  abcd[63] = 2.0E0*I_ERI_D2y_Px_Py_S_b;
  abcd[64] = 2.0E0*I_ERI_Dyz_Px_Py_S_b;
  abcd[65] = 2.0E0*I_ERI_D2z_Px_Py_S_b;
  abcd[66] = 2.0E0*I_ERI_D2x_Px_Pz_S_b;
  abcd[67] = 2.0E0*I_ERI_Dxy_Px_Pz_S_b;
  abcd[68] = 2.0E0*I_ERI_Dxz_Px_Pz_S_b;
  abcd[69] = 2.0E0*I_ERI_D2y_Px_Pz_S_b;
  abcd[70] = 2.0E0*I_ERI_Dyz_Px_Pz_S_b;
  abcd[71] = 2.0E0*I_ERI_D2z_Px_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[72] = 2.0E0*I_ERI_D2x_Py_Px_S_b;
  abcd[73] = 2.0E0*I_ERI_Dxy_Py_Px_S_b;
  abcd[74] = 2.0E0*I_ERI_Dxz_Py_Px_S_b;
  abcd[75] = 2.0E0*I_ERI_D2y_Py_Px_S_b;
  abcd[76] = 2.0E0*I_ERI_Dyz_Py_Px_S_b;
  abcd[77] = 2.0E0*I_ERI_D2z_Py_Px_S_b;
  abcd[78] = 2.0E0*I_ERI_D2x_Py_Py_S_b;
  abcd[79] = 2.0E0*I_ERI_Dxy_Py_Py_S_b;
  abcd[80] = 2.0E0*I_ERI_Dxz_Py_Py_S_b;
  abcd[81] = 2.0E0*I_ERI_D2y_Py_Py_S_b;
  abcd[82] = 2.0E0*I_ERI_Dyz_Py_Py_S_b;
  abcd[83] = 2.0E0*I_ERI_D2z_Py_Py_S_b;
  abcd[84] = 2.0E0*I_ERI_D2x_Py_Pz_S_b;
  abcd[85] = 2.0E0*I_ERI_Dxy_Py_Pz_S_b;
  abcd[86] = 2.0E0*I_ERI_Dxz_Py_Pz_S_b;
  abcd[87] = 2.0E0*I_ERI_D2y_Py_Pz_S_b;
  abcd[88] = 2.0E0*I_ERI_Dyz_Py_Pz_S_b;
  abcd[89] = 2.0E0*I_ERI_D2z_Py_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[90] = 2.0E0*I_ERI_D2x_Pz_Px_S_b;
  abcd[91] = 2.0E0*I_ERI_Dxy_Pz_Px_S_b;
  abcd[92] = 2.0E0*I_ERI_Dxz_Pz_Px_S_b;
  abcd[93] = 2.0E0*I_ERI_D2y_Pz_Px_S_b;
  abcd[94] = 2.0E0*I_ERI_Dyz_Pz_Px_S_b;
  abcd[95] = 2.0E0*I_ERI_D2z_Pz_Px_S_b;
  abcd[96] = 2.0E0*I_ERI_D2x_Pz_Py_S_b;
  abcd[97] = 2.0E0*I_ERI_Dxy_Pz_Py_S_b;
  abcd[98] = 2.0E0*I_ERI_Dxz_Pz_Py_S_b;
  abcd[99] = 2.0E0*I_ERI_D2y_Pz_Py_S_b;
  abcd[100] = 2.0E0*I_ERI_Dyz_Pz_Py_S_b;
  abcd[101] = 2.0E0*I_ERI_D2z_Pz_Py_S_b;
  abcd[102] = 2.0E0*I_ERI_D2x_Pz_Pz_S_b;
  abcd[103] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[104] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[105] = 2.0E0*I_ERI_D2y_Pz_Pz_S_b;
  abcd[106] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[107] = 2.0E0*I_ERI_D2z_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[108] = 2.0E0*I_ERI_D2x_S_D2x_S_c-1*I_ERI_D2x_S_S_S;
  abcd[109] = 2.0E0*I_ERI_Dxy_S_D2x_S_c-1*I_ERI_Dxy_S_S_S;
  abcd[110] = 2.0E0*I_ERI_Dxz_S_D2x_S_c-1*I_ERI_Dxz_S_S_S;
  abcd[111] = 2.0E0*I_ERI_D2y_S_D2x_S_c-1*I_ERI_D2y_S_S_S;
  abcd[112] = 2.0E0*I_ERI_Dyz_S_D2x_S_c-1*I_ERI_Dyz_S_S_S;
  abcd[113] = 2.0E0*I_ERI_D2z_S_D2x_S_c-1*I_ERI_D2z_S_S_S;
  abcd[114] = 2.0E0*I_ERI_D2x_S_Dxy_S_c;
  abcd[115] = 2.0E0*I_ERI_Dxy_S_Dxy_S_c;
  abcd[116] = 2.0E0*I_ERI_Dxz_S_Dxy_S_c;
  abcd[117] = 2.0E0*I_ERI_D2y_S_Dxy_S_c;
  abcd[118] = 2.0E0*I_ERI_Dyz_S_Dxy_S_c;
  abcd[119] = 2.0E0*I_ERI_D2z_S_Dxy_S_c;
  abcd[120] = 2.0E0*I_ERI_D2x_S_Dxz_S_c;
  abcd[121] = 2.0E0*I_ERI_Dxy_S_Dxz_S_c;
  abcd[122] = 2.0E0*I_ERI_Dxz_S_Dxz_S_c;
  abcd[123] = 2.0E0*I_ERI_D2y_S_Dxz_S_c;
  abcd[124] = 2.0E0*I_ERI_Dyz_S_Dxz_S_c;
  abcd[125] = 2.0E0*I_ERI_D2z_S_Dxz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[126] = 2.0E0*I_ERI_D2x_S_Dxy_S_c;
  abcd[127] = 2.0E0*I_ERI_Dxy_S_Dxy_S_c;
  abcd[128] = 2.0E0*I_ERI_Dxz_S_Dxy_S_c;
  abcd[129] = 2.0E0*I_ERI_D2y_S_Dxy_S_c;
  abcd[130] = 2.0E0*I_ERI_Dyz_S_Dxy_S_c;
  abcd[131] = 2.0E0*I_ERI_D2z_S_Dxy_S_c;
  abcd[132] = 2.0E0*I_ERI_D2x_S_D2y_S_c-1*I_ERI_D2x_S_S_S;
  abcd[133] = 2.0E0*I_ERI_Dxy_S_D2y_S_c-1*I_ERI_Dxy_S_S_S;
  abcd[134] = 2.0E0*I_ERI_Dxz_S_D2y_S_c-1*I_ERI_Dxz_S_S_S;
  abcd[135] = 2.0E0*I_ERI_D2y_S_D2y_S_c-1*I_ERI_D2y_S_S_S;
  abcd[136] = 2.0E0*I_ERI_Dyz_S_D2y_S_c-1*I_ERI_Dyz_S_S_S;
  abcd[137] = 2.0E0*I_ERI_D2z_S_D2y_S_c-1*I_ERI_D2z_S_S_S;
  abcd[138] = 2.0E0*I_ERI_D2x_S_Dyz_S_c;
  abcd[139] = 2.0E0*I_ERI_Dxy_S_Dyz_S_c;
  abcd[140] = 2.0E0*I_ERI_Dxz_S_Dyz_S_c;
  abcd[141] = 2.0E0*I_ERI_D2y_S_Dyz_S_c;
  abcd[142] = 2.0E0*I_ERI_Dyz_S_Dyz_S_c;
  abcd[143] = 2.0E0*I_ERI_D2z_S_Dyz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[144] = 2.0E0*I_ERI_D2x_S_Dxz_S_c;
  abcd[145] = 2.0E0*I_ERI_Dxy_S_Dxz_S_c;
  abcd[146] = 2.0E0*I_ERI_Dxz_S_Dxz_S_c;
  abcd[147] = 2.0E0*I_ERI_D2y_S_Dxz_S_c;
  abcd[148] = 2.0E0*I_ERI_Dyz_S_Dxz_S_c;
  abcd[149] = 2.0E0*I_ERI_D2z_S_Dxz_S_c;
  abcd[150] = 2.0E0*I_ERI_D2x_S_Dyz_S_c;
  abcd[151] = 2.0E0*I_ERI_Dxy_S_Dyz_S_c;
  abcd[152] = 2.0E0*I_ERI_Dxz_S_Dyz_S_c;
  abcd[153] = 2.0E0*I_ERI_D2y_S_Dyz_S_c;
  abcd[154] = 2.0E0*I_ERI_Dyz_S_Dyz_S_c;
  abcd[155] = 2.0E0*I_ERI_D2z_S_Dyz_S_c;
  abcd[156] = 2.0E0*I_ERI_D2x_S_D2z_S_c-1*I_ERI_D2x_S_S_S;
  abcd[157] = 2.0E0*I_ERI_Dxy_S_D2z_S_c-1*I_ERI_Dxy_S_S_S;
  abcd[158] = 2.0E0*I_ERI_Dxz_S_D2z_S_c-1*I_ERI_Dxz_S_S_S;
  abcd[159] = 2.0E0*I_ERI_D2y_S_D2z_S_c-1*I_ERI_D2y_S_S_S;
  abcd[160] = 2.0E0*I_ERI_Dyz_S_D2z_S_c-1*I_ERI_Dyz_S_S_S;
  abcd[161] = 2.0E0*I_ERI_D2z_S_D2z_S_c-1*I_ERI_D2z_S_S_S;
}
