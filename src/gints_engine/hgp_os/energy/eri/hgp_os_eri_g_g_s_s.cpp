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

void hgp_os_eri_g_g_s_s(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_L8x_S_S_S = 0.0E0;
  Double I_ERI_L7xy_S_S_S = 0.0E0;
  Double I_ERI_L7xz_S_S_S = 0.0E0;
  Double I_ERI_L6x2y_S_S_S = 0.0E0;
  Double I_ERI_L6xyz_S_S_S = 0.0E0;
  Double I_ERI_L6x2z_S_S_S = 0.0E0;
  Double I_ERI_L5x3y_S_S_S = 0.0E0;
  Double I_ERI_L5x2yz_S_S_S = 0.0E0;
  Double I_ERI_L5xy2z_S_S_S = 0.0E0;
  Double I_ERI_L5x3z_S_S_S = 0.0E0;
  Double I_ERI_L4x4y_S_S_S = 0.0E0;
  Double I_ERI_L4x3yz_S_S_S = 0.0E0;
  Double I_ERI_L4x2y2z_S_S_S = 0.0E0;
  Double I_ERI_L4xy3z_S_S_S = 0.0E0;
  Double I_ERI_L4x4z_S_S_S = 0.0E0;
  Double I_ERI_L3x5y_S_S_S = 0.0E0;
  Double I_ERI_L3x4yz_S_S_S = 0.0E0;
  Double I_ERI_L3x3y2z_S_S_S = 0.0E0;
  Double I_ERI_L3x2y3z_S_S_S = 0.0E0;
  Double I_ERI_L3xy4z_S_S_S = 0.0E0;
  Double I_ERI_L3x5z_S_S_S = 0.0E0;
  Double I_ERI_L2x6y_S_S_S = 0.0E0;
  Double I_ERI_L2x5yz_S_S_S = 0.0E0;
  Double I_ERI_L2x4y2z_S_S_S = 0.0E0;
  Double I_ERI_L2x3y3z_S_S_S = 0.0E0;
  Double I_ERI_L2x2y4z_S_S_S = 0.0E0;
  Double I_ERI_L2xy5z_S_S_S = 0.0E0;
  Double I_ERI_L2x6z_S_S_S = 0.0E0;
  Double I_ERI_Lx7y_S_S_S = 0.0E0;
  Double I_ERI_Lx6yz_S_S_S = 0.0E0;
  Double I_ERI_Lx5y2z_S_S_S = 0.0E0;
  Double I_ERI_Lx4y3z_S_S_S = 0.0E0;
  Double I_ERI_Lx3y4z_S_S_S = 0.0E0;
  Double I_ERI_Lx2y5z_S_S_S = 0.0E0;
  Double I_ERI_Lxy6z_S_S_S = 0.0E0;
  Double I_ERI_Lx7z_S_S_S = 0.0E0;
  Double I_ERI_L8y_S_S_S = 0.0E0;
  Double I_ERI_L7yz_S_S_S = 0.0E0;
  Double I_ERI_L6y2z_S_S_S = 0.0E0;
  Double I_ERI_L5y3z_S_S_S = 0.0E0;
  Double I_ERI_L4y4z_S_S_S = 0.0E0;
  Double I_ERI_L3y5z_S_S_S = 0.0E0;
  Double I_ERI_L2y6z_S_S_S = 0.0E0;
  Double I_ERI_Ly7z_S_S_S = 0.0E0;
  Double I_ERI_L8z_S_S_S = 0.0E0;
  Double I_ERI_K7x_S_S_S = 0.0E0;
  Double I_ERI_K6xy_S_S_S = 0.0E0;
  Double I_ERI_K6xz_S_S_S = 0.0E0;
  Double I_ERI_K5x2y_S_S_S = 0.0E0;
  Double I_ERI_K5xyz_S_S_S = 0.0E0;
  Double I_ERI_K5x2z_S_S_S = 0.0E0;
  Double I_ERI_K4x3y_S_S_S = 0.0E0;
  Double I_ERI_K4x2yz_S_S_S = 0.0E0;
  Double I_ERI_K4xy2z_S_S_S = 0.0E0;
  Double I_ERI_K4x3z_S_S_S = 0.0E0;
  Double I_ERI_K3x4y_S_S_S = 0.0E0;
  Double I_ERI_K3x3yz_S_S_S = 0.0E0;
  Double I_ERI_K3x2y2z_S_S_S = 0.0E0;
  Double I_ERI_K3xy3z_S_S_S = 0.0E0;
  Double I_ERI_K3x4z_S_S_S = 0.0E0;
  Double I_ERI_K2x5y_S_S_S = 0.0E0;
  Double I_ERI_K2x4yz_S_S_S = 0.0E0;
  Double I_ERI_K2x3y2z_S_S_S = 0.0E0;
  Double I_ERI_K2x2y3z_S_S_S = 0.0E0;
  Double I_ERI_K2xy4z_S_S_S = 0.0E0;
  Double I_ERI_K2x5z_S_S_S = 0.0E0;
  Double I_ERI_Kx6y_S_S_S = 0.0E0;
  Double I_ERI_Kx5yz_S_S_S = 0.0E0;
  Double I_ERI_Kx4y2z_S_S_S = 0.0E0;
  Double I_ERI_Kx3y3z_S_S_S = 0.0E0;
  Double I_ERI_Kx2y4z_S_S_S = 0.0E0;
  Double I_ERI_Kxy5z_S_S_S = 0.0E0;
  Double I_ERI_Kx6z_S_S_S = 0.0E0;
  Double I_ERI_K7y_S_S_S = 0.0E0;
  Double I_ERI_K6yz_S_S_S = 0.0E0;
  Double I_ERI_K5y2z_S_S_S = 0.0E0;
  Double I_ERI_K4y3z_S_S_S = 0.0E0;
  Double I_ERI_K3y4z_S_S_S = 0.0E0;
  Double I_ERI_K2y5z_S_S_S = 0.0E0;
  Double I_ERI_Ky6z_S_S_S = 0.0E0;
  Double I_ERI_K7z_S_S_S = 0.0E0;
  Double I_ERI_I6x_S_S_S = 0.0E0;
  Double I_ERI_I5xy_S_S_S = 0.0E0;
  Double I_ERI_I5xz_S_S_S = 0.0E0;
  Double I_ERI_I4x2y_S_S_S = 0.0E0;
  Double I_ERI_I4xyz_S_S_S = 0.0E0;
  Double I_ERI_I4x2z_S_S_S = 0.0E0;
  Double I_ERI_I3x3y_S_S_S = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S = 0.0E0;
  Double I_ERI_I3x3z_S_S_S = 0.0E0;
  Double I_ERI_I2x4y_S_S_S = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S = 0.0E0;
  Double I_ERI_I2x4z_S_S_S = 0.0E0;
  Double I_ERI_Ix5y_S_S_S = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S = 0.0E0;
  Double I_ERI_Ix5z_S_S_S = 0.0E0;
  Double I_ERI_I6y_S_S_S = 0.0E0;
  Double I_ERI_I5yz_S_S_S = 0.0E0;
  Double I_ERI_I4y2z_S_S_S = 0.0E0;
  Double I_ERI_I3y3z_S_S_S = 0.0E0;
  Double I_ERI_I2y4z_S_S_S = 0.0E0;
  Double I_ERI_Iy5z_S_S_S = 0.0E0;
  Double I_ERI_I6z_S_S_S = 0.0E0;
  Double I_ERI_H5x_S_S_S = 0.0E0;
  Double I_ERI_H4xy_S_S_S = 0.0E0;
  Double I_ERI_H4xz_S_S_S = 0.0E0;
  Double I_ERI_H3x2y_S_S_S = 0.0E0;
  Double I_ERI_H3xyz_S_S_S = 0.0E0;
  Double I_ERI_H3x2z_S_S_S = 0.0E0;
  Double I_ERI_H2x3y_S_S_S = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S = 0.0E0;
  Double I_ERI_H2x3z_S_S_S = 0.0E0;
  Double I_ERI_Hx4y_S_S_S = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S = 0.0E0;
  Double I_ERI_Hx4z_S_S_S = 0.0E0;
  Double I_ERI_H5y_S_S_S = 0.0E0;
  Double I_ERI_H4yz_S_S_S = 0.0E0;
  Double I_ERI_H3y2z_S_S_S = 0.0E0;
  Double I_ERI_H2y3z_S_S_S = 0.0E0;
  Double I_ERI_Hy4z_S_S_S = 0.0E0;
  Double I_ERI_H5z_S_S_S = 0.0E0;
  Double I_ERI_G4x_S_S_S = 0.0E0;
  Double I_ERI_G3xy_S_S_S = 0.0E0;
  Double I_ERI_G3xz_S_S_S = 0.0E0;
  Double I_ERI_G2x2y_S_S_S = 0.0E0;
  Double I_ERI_G2xyz_S_S_S = 0.0E0;
  Double I_ERI_G2x2z_S_S_S = 0.0E0;
  Double I_ERI_Gx3y_S_S_S = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S = 0.0E0;
  Double I_ERI_Gx3z_S_S_S = 0.0E0;
  Double I_ERI_G4y_S_S_S = 0.0E0;
  Double I_ERI_G3yz_S_S_S = 0.0E0;
  Double I_ERI_G2y2z_S_S_S = 0.0E0;
  Double I_ERI_Gy3z_S_S_S = 0.0E0;
  Double I_ERI_G4z_S_S_S = 0.0E0;

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

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
      Double WX    = rho*(PX*onede + QX*onedz);
      Double WY    = rho*(PY*onede + QY*onedz);
      Double WZ    = rho*(PZ*onede + QZ*onedz);
      Double oned2k= 0.5E0*rho*onede*onedz;
      Double WPX   = WX - PX;
      Double WPY   = WY - PY;
      Double WPZ   = WZ - PZ;
      Double rhod2zsq = rho*oned2z*onedz;


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
       * shell quartet name: SQ_ERI_P_S_S_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_Px_S_S_S_M7_vrr = PAX*I_ERI_S_S_S_S_M7_vrr+WPX*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_Py_S_S_S_M7_vrr = PAY*I_ERI_S_S_S_S_M7_vrr+WPY*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_Pz_S_S_S_M7_vrr = PAZ*I_ERI_S_S_S_S_M7_vrr+WPZ*I_ERI_S_S_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_Px_S_S_S_M6_vrr = PAX*I_ERI_S_S_S_S_M6_vrr+WPX*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Py_S_S_S_M6_vrr = PAY*I_ERI_S_S_S_S_M6_vrr+WPY*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Pz_S_S_S_M6_vrr = PAZ*I_ERI_S_S_S_S_M6_vrr+WPZ*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M6_vrr = PAX*I_ERI_Px_S_S_S_M6_vrr+WPX*I_ERI_Px_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_D2y_S_S_S_M6_vrr = PAY*I_ERI_Py_S_S_S_M6_vrr+WPY*I_ERI_Py_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_D2z_S_S_S_M6_vrr = PAZ*I_ERI_Pz_S_S_S_M6_vrr+WPZ*I_ERI_Pz_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_Px_S_S_S_M5_vrr = PAX*I_ERI_S_S_S_S_M5_vrr+WPX*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Py_S_S_S_M5_vrr = PAY*I_ERI_S_S_S_S_M5_vrr+WPY*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Pz_S_S_S_M5_vrr = PAZ*I_ERI_S_S_S_S_M5_vrr+WPZ*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M5_vrr = PAX*I_ERI_Px_S_S_S_M5_vrr+WPX*I_ERI_Px_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_D2y_S_S_S_M5_vrr = PAY*I_ERI_Py_S_S_S_M5_vrr+WPY*I_ERI_Py_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_D2z_S_S_S_M5_vrr = PAZ*I_ERI_Pz_S_S_S_M5_vrr+WPZ*I_ERI_Pz_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M5_vrr = PAX*I_ERI_D2x_S_S_S_M5_vrr+WPX*I_ERI_D2x_S_S_S_M6_vrr+2*oned2z*I_ERI_Px_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M6_vrr;
      Double I_ERI_F3y_S_S_S_M5_vrr = PAY*I_ERI_D2y_S_S_S_M5_vrr+WPY*I_ERI_D2y_S_S_S_M6_vrr+2*oned2z*I_ERI_Py_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M6_vrr;
      Double I_ERI_F3z_S_S_S_M5_vrr = PAZ*I_ERI_D2z_S_S_S_M5_vrr+WPZ*I_ERI_D2z_S_S_S_M6_vrr+2*oned2z*I_ERI_Pz_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_Px_S_S_S_M4_vrr = PAX*I_ERI_S_S_S_S_M4_vrr+WPX*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Py_S_S_S_M4_vrr = PAY*I_ERI_S_S_S_S_M4_vrr+WPY*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Pz_S_S_S_M4_vrr = PAZ*I_ERI_S_S_S_S_M4_vrr+WPZ*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M4_vrr = PAX*I_ERI_Px_S_S_S_M4_vrr+WPX*I_ERI_Px_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_D2y_S_S_S_M4_vrr = PAY*I_ERI_Py_S_S_S_M4_vrr+WPY*I_ERI_Py_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_D2z_S_S_S_M4_vrr = PAZ*I_ERI_Pz_S_S_S_M4_vrr+WPZ*I_ERI_Pz_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M4_vrr = PAX*I_ERI_F3x_S_S_S_M4_vrr+WPX*I_ERI_F3x_S_S_S_M5_vrr+3*oned2z*I_ERI_D2x_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_G3xy_S_S_S_M4_vrr = PAY*I_ERI_F3x_S_S_S_M4_vrr+WPY*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_G3xz_S_S_S_M4_vrr = PAZ*I_ERI_F3x_S_S_S_M4_vrr+WPZ*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_G4y_S_S_S_M4_vrr = PAY*I_ERI_F3y_S_S_S_M4_vrr+WPY*I_ERI_F3y_S_S_S_M5_vrr+3*oned2z*I_ERI_D2y_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_G3yz_S_S_S_M4_vrr = PAZ*I_ERI_F3y_S_S_S_M4_vrr+WPZ*I_ERI_F3y_S_S_S_M5_vrr;
      Double I_ERI_G4z_S_S_S_M4_vrr = PAZ*I_ERI_F3z_S_S_S_M4_vrr+WPZ*I_ERI_F3z_S_S_S_M5_vrr+3*oned2z*I_ERI_D2z_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M3_vrr = PAX*I_ERI_Px_S_S_S_M3_vrr+WPX*I_ERI_Px_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_S_S_M3_vrr = PAY*I_ERI_Py_S_S_S_M3_vrr+WPY*I_ERI_Py_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_S_S_M3_vrr = PAZ*I_ERI_Pz_S_S_S_M3_vrr+WPZ*I_ERI_Pz_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M3_vrr = PAX*I_ERI_D2x_S_S_S_M3_vrr+WPX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2z*I_ERI_Px_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_S_S_M3_vrr = PAY*I_ERI_D2x_S_S_S_M3_vrr+WPY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_S_S_M3_vrr = PAY*I_ERI_D2y_S_S_S_M3_vrr+WPY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2z*I_ERI_Py_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_S_S_M3_vrr = PAZ*I_ERI_D2z_S_S_S_M3_vrr+WPZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2z*I_ERI_Pz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M3_vrr = PAX*I_ERI_G4x_S_S_S_M3_vrr+WPX*I_ERI_G4x_S_S_S_M4_vrr+4*oned2z*I_ERI_F3x_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_H4xy_S_S_S_M3_vrr = PAY*I_ERI_G4x_S_S_S_M3_vrr+WPY*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_H4xz_S_S_S_M3_vrr = PAZ*I_ERI_G4x_S_S_S_M3_vrr+WPZ*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_H3x2y_S_S_S_M3_vrr = PAY*I_ERI_G3xy_S_S_S_M3_vrr+WPY*I_ERI_G3xy_S_S_S_M4_vrr+oned2z*I_ERI_F3x_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_H3x2z_S_S_S_M3_vrr = PAZ*I_ERI_G3xz_S_S_S_M3_vrr+WPZ*I_ERI_G3xz_S_S_S_M4_vrr+oned2z*I_ERI_F3x_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_Hx4y_S_S_S_M3_vrr = PAX*I_ERI_G4y_S_S_S_M3_vrr+WPX*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_Hx4z_S_S_S_M3_vrr = PAX*I_ERI_G4z_S_S_S_M3_vrr+WPX*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_H5y_S_S_S_M3_vrr = PAY*I_ERI_G4y_S_S_S_M3_vrr+WPY*I_ERI_G4y_S_S_S_M4_vrr+4*oned2z*I_ERI_F3y_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_H4yz_S_S_S_M3_vrr = PAZ*I_ERI_G4y_S_S_S_M3_vrr+WPZ*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_H3y2z_S_S_S_M3_vrr = PAZ*I_ERI_G3yz_S_S_S_M3_vrr+WPZ*I_ERI_G3yz_S_S_S_M4_vrr+oned2z*I_ERI_F3y_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Hy4z_S_S_S_M3_vrr = PAY*I_ERI_G4z_S_S_S_M3_vrr+WPY*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_H5z_S_S_S_M3_vrr = PAZ*I_ERI_G4z_S_S_S_M3_vrr+WPZ*I_ERI_G4z_S_S_S_M4_vrr+4*oned2z*I_ERI_F3z_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M4_vrr;

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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M2_vrr = PAX*I_ERI_Px_S_S_S_M2_vrr+WPX*I_ERI_Px_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_S_S_M2_vrr = PAY*I_ERI_Py_S_S_S_M2_vrr+WPY*I_ERI_Py_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_S_S_M2_vrr = PAZ*I_ERI_Pz_S_S_S_M2_vrr+WPZ*I_ERI_Pz_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M2_vrr = PAX*I_ERI_D2x_S_S_S_M2_vrr+WPX*I_ERI_D2x_S_S_S_M3_vrr+2*oned2z*I_ERI_Px_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_S_S_M2_vrr = PAY*I_ERI_D2x_S_S_S_M2_vrr+WPY*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_S_S_M2_vrr = PAZ*I_ERI_D2x_S_S_S_M2_vrr+WPZ*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_S_S_M2_vrr = PAY*I_ERI_D2y_S_S_S_M2_vrr+WPY*I_ERI_D2y_S_S_S_M3_vrr+2*oned2z*I_ERI_Py_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_S_S_M2_vrr = PAZ*I_ERI_D2y_S_S_S_M2_vrr+WPZ*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_S_S_M2_vrr = PAZ*I_ERI_D2z_S_S_S_M2_vrr+WPZ*I_ERI_D2z_S_S_S_M3_vrr+2*oned2z*I_ERI_Pz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M2_vrr = PAX*I_ERI_F3x_S_S_S_M2_vrr+WPX*I_ERI_F3x_S_S_S_M3_vrr+3*oned2z*I_ERI_D2x_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_S_S_M2_vrr = PAY*I_ERI_F3x_S_S_S_M2_vrr+WPY*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_S_S_M2_vrr = PAZ*I_ERI_F3x_S_S_S_M2_vrr+WPZ*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_S_S_M2_vrr = PAY*I_ERI_F2xy_S_S_S_M2_vrr+WPY*I_ERI_F2xy_S_S_S_M3_vrr+oned2z*I_ERI_D2x_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_S_S_M2_vrr = PAX*I_ERI_F3y_S_S_S_M2_vrr+WPX*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_S_S_M2_vrr = PAX*I_ERI_F3z_S_S_S_M2_vrr+WPX*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_S_S_M2_vrr = PAY*I_ERI_F3y_S_S_S_M2_vrr+WPY*I_ERI_F3y_S_S_S_M3_vrr+3*oned2z*I_ERI_D2y_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_S_S_M2_vrr = PAZ*I_ERI_F3y_S_S_S_M2_vrr+WPZ*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_S_S_M2_vrr = PAY*I_ERI_F3z_S_S_S_M2_vrr+WPY*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_S_S_M2_vrr = PAZ*I_ERI_F3z_S_S_S_M2_vrr+WPZ*I_ERI_F3z_S_S_S_M3_vrr+3*oned2z*I_ERI_D2z_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M2_vrr = PAX*I_ERI_G4x_S_S_S_M2_vrr+WPX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2z*I_ERI_F3x_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H4xy_S_S_S_M2_vrr = PAY*I_ERI_G4x_S_S_S_M2_vrr+WPY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H4xz_S_S_S_M2_vrr = PAZ*I_ERI_G4x_S_S_S_M2_vrr+WPZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H3x2y_S_S_S_M2_vrr = PAY*I_ERI_G3xy_S_S_S_M2_vrr+WPY*I_ERI_G3xy_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H3x2z_S_S_S_M2_vrr = PAZ*I_ERI_G3xz_S_S_S_M2_vrr+WPZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H2x3y_S_S_S_M2_vrr = PAX*I_ERI_Gx3y_S_S_S_M2_vrr+WPX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2x3z_S_S_S_M2_vrr = PAX*I_ERI_Gx3z_S_S_S_M2_vrr+WPX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4y_S_S_S_M2_vrr = PAX*I_ERI_G4y_S_S_S_M2_vrr+WPX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_Hx4z_S_S_S_M2_vrr = PAX*I_ERI_G4z_S_S_S_M2_vrr+WPX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5y_S_S_S_M2_vrr = PAY*I_ERI_G4y_S_S_S_M2_vrr+WPY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2z*I_ERI_F3y_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H4yz_S_S_S_M2_vrr = PAZ*I_ERI_G4y_S_S_S_M2_vrr+WPZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_H3y2z_S_S_S_M2_vrr = PAZ*I_ERI_G3yz_S_S_S_M2_vrr+WPZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Hy4z_S_S_S_M2_vrr = PAY*I_ERI_G4z_S_S_S_M2_vrr+WPY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5z_S_S_S_M2_vrr = PAZ*I_ERI_G4z_S_S_S_M2_vrr+WPZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2z*I_ERI_F3z_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 10 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M2_vrr = PAX*I_ERI_H5x_S_S_S_M2_vrr+WPX*I_ERI_H5x_S_S_S_M3_vrr+5*oned2z*I_ERI_G4x_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I5xy_S_S_S_M2_vrr = PAY*I_ERI_H5x_S_S_S_M2_vrr+WPY*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_I5xz_S_S_S_M2_vrr = PAZ*I_ERI_H5x_S_S_S_M2_vrr+WPZ*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_I4x2y_S_S_S_M2_vrr = PAY*I_ERI_H4xy_S_S_S_M2_vrr+WPY*I_ERI_H4xy_S_S_S_M3_vrr+oned2z*I_ERI_G4x_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I4x2z_S_S_S_M2_vrr = PAZ*I_ERI_H4xz_S_S_S_M2_vrr+WPZ*I_ERI_H4xz_S_S_S_M3_vrr+oned2z*I_ERI_G4x_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I3x3y_S_S_S_M2_vrr = PAY*I_ERI_H3x2y_S_S_S_M2_vrr+WPY*I_ERI_H3x2y_S_S_S_M3_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_I3x3z_S_S_S_M2_vrr = PAZ*I_ERI_H3x2z_S_S_S_M2_vrr+WPZ*I_ERI_H3x2z_S_S_S_M3_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M3_vrr;
      Double I_ERI_I2x4y_S_S_S_M2_vrr = PAX*I_ERI_Hx4y_S_S_S_M2_vrr+WPX*I_ERI_Hx4y_S_S_S_M3_vrr+oned2z*I_ERI_G4y_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I2x4z_S_S_S_M2_vrr = PAX*I_ERI_Hx4z_S_S_S_M2_vrr+WPX*I_ERI_Hx4z_S_S_S_M3_vrr+oned2z*I_ERI_G4z_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_Ix5y_S_S_S_M2_vrr = PAX*I_ERI_H5y_S_S_S_M2_vrr+WPX*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_Ix5z_S_S_S_M2_vrr = PAX*I_ERI_H5z_S_S_S_M2_vrr+WPX*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_I6y_S_S_S_M2_vrr = PAY*I_ERI_H5y_S_S_S_M2_vrr+WPY*I_ERI_H5y_S_S_S_M3_vrr+5*oned2z*I_ERI_G4y_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I5yz_S_S_S_M2_vrr = PAZ*I_ERI_H5y_S_S_S_M2_vrr+WPZ*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_I4y2z_S_S_S_M2_vrr = PAZ*I_ERI_H4yz_S_S_S_M2_vrr+WPZ*I_ERI_H4yz_S_S_S_M3_vrr+oned2z*I_ERI_G4y_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I3y3z_S_S_S_M2_vrr = PAZ*I_ERI_H3y2z_S_S_S_M2_vrr+WPZ*I_ERI_H3y2z_S_S_S_M3_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M3_vrr;
      Double I_ERI_I2y4z_S_S_S_M2_vrr = PAY*I_ERI_Hy4z_S_S_S_M2_vrr+WPY*I_ERI_Hy4z_S_S_S_M3_vrr+oned2z*I_ERI_G4z_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_Iy5z_S_S_S_M2_vrr = PAY*I_ERI_H5z_S_S_S_M2_vrr+WPY*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_I6z_S_S_S_M2_vrr = PAZ*I_ERI_H5z_S_S_S_M2_vrr+WPZ*I_ERI_H5z_S_S_S_M3_vrr+5*oned2z*I_ERI_G4z_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M1_vrr = PAX*I_ERI_Px_S_S_S_M1_vrr+WPX*I_ERI_Px_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_S_S_M1_vrr = PAY*I_ERI_Py_S_S_S_M1_vrr+WPY*I_ERI_Py_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_S_S_M1_vrr = PAZ*I_ERI_Pz_S_S_S_M1_vrr+WPZ*I_ERI_Pz_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M1_vrr = PAX*I_ERI_D2x_S_S_S_M1_vrr+WPX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2z*I_ERI_Px_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_S_S_M1_vrr = PAY*I_ERI_D2x_S_S_S_M1_vrr+WPY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_S_S_M1_vrr = PAZ*I_ERI_D2x_S_S_S_M1_vrr+WPZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_S_S_M1_vrr = PAX*I_ERI_D2y_S_S_S_M1_vrr+WPX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_S_S_M1_vrr = PAX*I_ERI_D2z_S_S_S_M1_vrr+WPX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_S_S_M1_vrr = PAY*I_ERI_D2y_S_S_S_M1_vrr+WPY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2z*I_ERI_Py_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_S_S_M1_vrr = PAZ*I_ERI_D2y_S_S_S_M1_vrr+WPZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_S_S_M1_vrr = PAZ*I_ERI_D2z_S_S_S_M1_vrr+WPZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2z*I_ERI_Pz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M1_vrr = PAX*I_ERI_F3x_S_S_S_M1_vrr+WPX*I_ERI_F3x_S_S_S_M2_vrr+3*oned2z*I_ERI_D2x_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_S_S_M1_vrr = PAY*I_ERI_F3x_S_S_S_M1_vrr+WPY*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_S_S_M1_vrr = PAZ*I_ERI_F3x_S_S_S_M1_vrr+WPZ*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_S_S_M1_vrr = PAY*I_ERI_F2xy_S_S_S_M1_vrr+WPY*I_ERI_F2xy_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_S_S_M1_vrr = PAZ*I_ERI_F2xz_S_S_S_M1_vrr+WPZ*I_ERI_F2xz_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_S_S_M1_vrr = PAX*I_ERI_F3y_S_S_S_M1_vrr+WPX*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_S_S_M1_vrr = PAX*I_ERI_F3z_S_S_S_M1_vrr+WPX*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_S_S_M1_vrr = PAY*I_ERI_F3y_S_S_S_M1_vrr+WPY*I_ERI_F3y_S_S_S_M2_vrr+3*oned2z*I_ERI_D2y_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_S_S_M1_vrr = PAZ*I_ERI_F3y_S_S_S_M1_vrr+WPZ*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_S_S_M1_vrr = PAZ*I_ERI_F2yz_S_S_S_M1_vrr+WPZ*I_ERI_F2yz_S_S_S_M2_vrr+oned2z*I_ERI_D2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_S_S_M1_vrr = PAY*I_ERI_F3z_S_S_S_M1_vrr+WPY*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_S_S_M1_vrr = PAZ*I_ERI_F3z_S_S_S_M1_vrr+WPZ*I_ERI_F3z_S_S_S_M2_vrr+3*oned2z*I_ERI_D2z_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M1_vrr = PAX*I_ERI_G4x_S_S_S_M1_vrr+WPX*I_ERI_G4x_S_S_S_M2_vrr+4*oned2z*I_ERI_F3x_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_S_S_M1_vrr = PAY*I_ERI_G4x_S_S_S_M1_vrr+WPY*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_S_S_M1_vrr = PAZ*I_ERI_G4x_S_S_S_M1_vrr+WPZ*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_S_S_M1_vrr = PAY*I_ERI_G3xy_S_S_S_M1_vrr+WPY*I_ERI_G3xy_S_S_S_M2_vrr+oned2z*I_ERI_F3x_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_S_S_M1_vrr = PAZ*I_ERI_G3xz_S_S_S_M1_vrr+WPZ*I_ERI_G3xz_S_S_S_M2_vrr+oned2z*I_ERI_F3x_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_S_S_M1_vrr = PAX*I_ERI_Gx3y_S_S_S_M1_vrr+WPX*I_ERI_Gx3y_S_S_S_M2_vrr+oned2z*I_ERI_F3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_S_S_M1_vrr = PAZ*I_ERI_G2x2y_S_S_S_M1_vrr+WPZ*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_S_S_M1_vrr = PAX*I_ERI_Gx3z_S_S_S_M1_vrr+WPX*I_ERI_Gx3z_S_S_S_M2_vrr+oned2z*I_ERI_F3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_S_S_M1_vrr = PAX*I_ERI_G4y_S_S_S_M1_vrr+WPX*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_S_S_M1_vrr = PAX*I_ERI_G4z_S_S_S_M1_vrr+WPX*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_S_S_M1_vrr = PAY*I_ERI_G4y_S_S_S_M1_vrr+WPY*I_ERI_G4y_S_S_S_M2_vrr+4*oned2z*I_ERI_F3y_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_S_S_M1_vrr = PAZ*I_ERI_G4y_S_S_S_M1_vrr+WPZ*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_S_S_M1_vrr = PAZ*I_ERI_G3yz_S_S_S_M1_vrr+WPZ*I_ERI_G3yz_S_S_S_M2_vrr+oned2z*I_ERI_F3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_S_S_M1_vrr = PAY*I_ERI_Gy3z_S_S_S_M1_vrr+WPY*I_ERI_Gy3z_S_S_S_M2_vrr+oned2z*I_ERI_F3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_S_S_M1_vrr = PAY*I_ERI_G4z_S_S_S_M1_vrr+WPY*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_S_S_M1_vrr = PAZ*I_ERI_G4z_S_S_S_M1_vrr+WPZ*I_ERI_G4z_S_S_S_M2_vrr+4*oned2z*I_ERI_F3z_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M1_vrr = PAX*I_ERI_H5x_S_S_S_M1_vrr+WPX*I_ERI_H5x_S_S_S_M2_vrr+5*oned2z*I_ERI_G4x_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I5xy_S_S_S_M1_vrr = PAY*I_ERI_H5x_S_S_S_M1_vrr+WPY*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I5xz_S_S_S_M1_vrr = PAZ*I_ERI_H5x_S_S_S_M1_vrr+WPZ*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I4x2y_S_S_S_M1_vrr = PAY*I_ERI_H4xy_S_S_S_M1_vrr+WPY*I_ERI_H4xy_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I4x2z_S_S_S_M1_vrr = PAZ*I_ERI_H4xz_S_S_S_M1_vrr+WPZ*I_ERI_H4xz_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I3x3y_S_S_S_M1_vrr = PAY*I_ERI_H3x2y_S_S_S_M1_vrr+WPY*I_ERI_H3x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_I3x2yz_S_S_S_M1_vrr = PAZ*I_ERI_H3x2y_S_S_S_M1_vrr+WPZ*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_I3x3z_S_S_S_M1_vrr = PAZ*I_ERI_H3x2z_S_S_S_M1_vrr+WPZ*I_ERI_H3x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_I2x4y_S_S_S_M1_vrr = PAX*I_ERI_Hx4y_S_S_S_M1_vrr+WPX*I_ERI_Hx4y_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I2x3yz_S_S_S_M1_vrr = PAZ*I_ERI_H2x3y_S_S_S_M1_vrr+WPZ*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_I2xy3z_S_S_S_M1_vrr = PAY*I_ERI_H2x3z_S_S_S_M1_vrr+WPY*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_I2x4z_S_S_S_M1_vrr = PAX*I_ERI_Hx4z_S_S_S_M1_vrr+WPX*I_ERI_Hx4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Ix5y_S_S_S_M1_vrr = PAX*I_ERI_H5y_S_S_S_M1_vrr+WPX*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_Ix5z_S_S_S_M1_vrr = PAX*I_ERI_H5z_S_S_S_M1_vrr+WPX*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6y_S_S_S_M1_vrr = PAY*I_ERI_H5y_S_S_S_M1_vrr+WPY*I_ERI_H5y_S_S_S_M2_vrr+5*oned2z*I_ERI_G4y_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I5yz_S_S_S_M1_vrr = PAZ*I_ERI_H5y_S_S_S_M1_vrr+WPZ*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_I4y2z_S_S_S_M1_vrr = PAZ*I_ERI_H4yz_S_S_S_M1_vrr+WPZ*I_ERI_H4yz_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I3y3z_S_S_S_M1_vrr = PAZ*I_ERI_H3y2z_S_S_S_M1_vrr+WPZ*I_ERI_H3y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_I2y4z_S_S_S_M1_vrr = PAY*I_ERI_Hy4z_S_S_S_M1_vrr+WPY*I_ERI_Hy4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Iy5z_S_S_S_M1_vrr = PAY*I_ERI_H5z_S_S_S_M1_vrr+WPY*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6z_S_S_S_M1_vrr = PAZ*I_ERI_H5z_S_S_S_M1_vrr+WPZ*I_ERI_H5z_S_S_S_M2_vrr+5*oned2z*I_ERI_G4z_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       ************************************************************/
      Double I_ERI_K7x_S_S_S_M1_vrr = PAX*I_ERI_I6x_S_S_S_M1_vrr+WPX*I_ERI_I6x_S_S_S_M2_vrr+6*oned2z*I_ERI_H5x_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K6xy_S_S_S_M1_vrr = PAY*I_ERI_I6x_S_S_S_M1_vrr+WPY*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_K6xz_S_S_S_M1_vrr = PAZ*I_ERI_I6x_S_S_S_M1_vrr+WPZ*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_K5x2y_S_S_S_M1_vrr = PAY*I_ERI_I5xy_S_S_S_M1_vrr+WPY*I_ERI_I5xy_S_S_S_M2_vrr+oned2z*I_ERI_H5x_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K5x2z_S_S_S_M1_vrr = PAZ*I_ERI_I5xz_S_S_S_M1_vrr+WPZ*I_ERI_I5xz_S_S_S_M2_vrr+oned2z*I_ERI_H5x_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K4x3y_S_S_S_M1_vrr = PAY*I_ERI_I4x2y_S_S_S_M1_vrr+WPY*I_ERI_I4x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_H4xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M2_vrr;
      Double I_ERI_K4x2yz_S_S_S_M1_vrr = PAZ*I_ERI_I4x2y_S_S_S_M1_vrr+WPZ*I_ERI_I4x2y_S_S_S_M2_vrr;
      Double I_ERI_K4x3z_S_S_S_M1_vrr = PAZ*I_ERI_I4x2z_S_S_S_M1_vrr+WPZ*I_ERI_I4x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_H4xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M2_vrr;
      Double I_ERI_K3x4y_S_S_S_M1_vrr = PAX*I_ERI_I2x4y_S_S_S_M1_vrr+WPX*I_ERI_I2x4y_S_S_S_M2_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M2_vrr;
      Double I_ERI_K3x3yz_S_S_S_M1_vrr = PAZ*I_ERI_I3x3y_S_S_S_M1_vrr+WPZ*I_ERI_I3x3y_S_S_S_M2_vrr;
      Double I_ERI_K3xy3z_S_S_S_M1_vrr = PAY*I_ERI_I3x3z_S_S_S_M1_vrr+WPY*I_ERI_I3x3z_S_S_S_M2_vrr;
      Double I_ERI_K3x4z_S_S_S_M1_vrr = PAX*I_ERI_I2x4z_S_S_S_M1_vrr+WPX*I_ERI_I2x4z_S_S_S_M2_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M2_vrr;
      Double I_ERI_K2x5y_S_S_S_M1_vrr = PAX*I_ERI_Ix5y_S_S_S_M1_vrr+WPX*I_ERI_Ix5y_S_S_S_M2_vrr+oned2z*I_ERI_H5y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K2x4yz_S_S_S_M1_vrr = PAZ*I_ERI_I2x4y_S_S_S_M1_vrr+WPZ*I_ERI_I2x4y_S_S_S_M2_vrr;
      Double I_ERI_K2xy4z_S_S_S_M1_vrr = PAY*I_ERI_I2x4z_S_S_S_M1_vrr+WPY*I_ERI_I2x4z_S_S_S_M2_vrr;
      Double I_ERI_K2x5z_S_S_S_M1_vrr = PAX*I_ERI_Ix5z_S_S_S_M1_vrr+WPX*I_ERI_Ix5z_S_S_S_M2_vrr+oned2z*I_ERI_H5z_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_Kx6y_S_S_S_M1_vrr = PAX*I_ERI_I6y_S_S_S_M1_vrr+WPX*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_Kx3y3z_S_S_S_M1_vrr = PAX*I_ERI_I3y3z_S_S_S_M1_vrr+WPX*I_ERI_I3y3z_S_S_S_M2_vrr;
      Double I_ERI_Kx6z_S_S_S_M1_vrr = PAX*I_ERI_I6z_S_S_S_M1_vrr+WPX*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_K7y_S_S_S_M1_vrr = PAY*I_ERI_I6y_S_S_S_M1_vrr+WPY*I_ERI_I6y_S_S_S_M2_vrr+6*oned2z*I_ERI_H5y_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K6yz_S_S_S_M1_vrr = PAZ*I_ERI_I6y_S_S_S_M1_vrr+WPZ*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_K5y2z_S_S_S_M1_vrr = PAZ*I_ERI_I5yz_S_S_S_M1_vrr+WPZ*I_ERI_I5yz_S_S_S_M2_vrr+oned2z*I_ERI_H5y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K4y3z_S_S_S_M1_vrr = PAZ*I_ERI_I4y2z_S_S_S_M1_vrr+WPZ*I_ERI_I4y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_H4yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M2_vrr;
      Double I_ERI_K3y4z_S_S_S_M1_vrr = PAY*I_ERI_I2y4z_S_S_S_M1_vrr+WPY*I_ERI_I2y4z_S_S_S_M2_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M2_vrr;
      Double I_ERI_K2y5z_S_S_S_M1_vrr = PAY*I_ERI_Iy5z_S_S_S_M1_vrr+WPY*I_ERI_Iy5z_S_S_S_M2_vrr+oned2z*I_ERI_H5z_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_Ky6z_S_S_S_M1_vrr = PAY*I_ERI_I6z_S_S_S_M1_vrr+WPY*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_K7z_S_S_S_M1_vrr = PAZ*I_ERI_I6z_S_S_S_M1_vrr+WPZ*I_ERI_I6z_S_S_S_M2_vrr+6*oned2z*I_ERI_H5z_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_S_vrr = PAZ*I_ERI_Pz_S_S_S_vrr+WPZ*I_ERI_Pz_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_S_S_vrr = PAX*I_ERI_D2x_S_S_S_vrr+WPX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2z*I_ERI_Px_S_S_S_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_S_vrr = PAY*I_ERI_D2x_S_S_S_vrr+WPY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_S_vrr = PAZ*I_ERI_D2x_S_S_S_vrr+WPZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_S_vrr = PAX*I_ERI_D2y_S_S_S_vrr+WPX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_S_vrr = PAX*I_ERI_D2z_S_S_S_vrr+WPX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_S_vrr = PAY*I_ERI_D2y_S_S_S_vrr+WPY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2z*I_ERI_Py_S_S_S_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_S_vrr = PAZ*I_ERI_D2y_S_S_S_vrr+WPZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_S_vrr = PAZ*I_ERI_D2z_S_S_S_vrr+WPZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2z*I_ERI_Pz_S_S_S_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_S_S_vrr = PAX*I_ERI_F3x_S_S_S_vrr+WPX*I_ERI_F3x_S_S_S_M1_vrr+3*oned2z*I_ERI_D2x_S_S_S_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_S_S_vrr = PAY*I_ERI_F3x_S_S_S_vrr+WPY*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_S_S_vrr = PAZ*I_ERI_F3x_S_S_S_vrr+WPZ*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_S_S_vrr = PAY*I_ERI_F2xy_S_S_S_vrr+WPY*I_ERI_F2xy_S_S_S_M1_vrr+oned2z*I_ERI_D2x_S_S_S_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_S_S_vrr = PAZ*I_ERI_F2xy_S_S_S_vrr+WPZ*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_S_S_vrr = PAZ*I_ERI_F2xz_S_S_S_vrr+WPZ*I_ERI_F2xz_S_S_S_M1_vrr+oned2z*I_ERI_D2x_S_S_S_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_S_S_vrr = PAX*I_ERI_F3y_S_S_S_vrr+WPX*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_S_S_vrr = PAZ*I_ERI_Fx2y_S_S_S_vrr+WPZ*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_S_S_vrr = PAY*I_ERI_Fx2z_S_S_S_vrr+WPY*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_S_S_vrr = PAX*I_ERI_F3z_S_S_S_vrr+WPX*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_S_S_vrr = PAY*I_ERI_F3y_S_S_S_vrr+WPY*I_ERI_F3y_S_S_S_M1_vrr+3*oned2z*I_ERI_D2y_S_S_S_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_S_S_vrr = PAZ*I_ERI_F3y_S_S_S_vrr+WPZ*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_S_S_vrr = PAZ*I_ERI_F2yz_S_S_S_vrr+WPZ*I_ERI_F2yz_S_S_S_M1_vrr+oned2z*I_ERI_D2y_S_S_S_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_S_S_vrr = PAY*I_ERI_F3z_S_S_S_vrr+WPY*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_S_S_vrr = PAZ*I_ERI_F3z_S_S_S_vrr+WPZ*I_ERI_F3z_S_S_S_M1_vrr+3*oned2z*I_ERI_D2z_S_S_S_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_S_S_vrr = PAX*I_ERI_G4x_S_S_S_vrr+WPX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2z*I_ERI_F3x_S_S_S_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_S_S_vrr = PAY*I_ERI_G4x_S_S_S_vrr+WPY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_S_S_vrr = PAZ*I_ERI_G4x_S_S_S_vrr+WPZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_S_S_vrr = PAY*I_ERI_G3xy_S_S_S_vrr+WPY*I_ERI_G3xy_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_S_S_vrr = PAZ*I_ERI_G3xy_S_S_S_vrr+WPZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_S_S_vrr = PAZ*I_ERI_G3xz_S_S_S_vrr+WPZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_S_S_vrr = PAX*I_ERI_Gx3y_S_S_S_vrr+WPX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_S_S_vrr = PAZ*I_ERI_G2x2y_S_S_S_vrr+WPZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_S_S_vrr = PAY*I_ERI_G2x2z_S_S_S_vrr+WPY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_S_S_vrr = PAX*I_ERI_Gx3z_S_S_S_vrr+WPX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_S_S_vrr = PAX*I_ERI_G4y_S_S_S_vrr+WPX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_S_S_vrr = PAZ*I_ERI_Gx3y_S_S_S_vrr+WPZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_S_S_vrr = PAX*I_ERI_G2y2z_S_S_S_vrr+WPX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_S_S_vrr = PAY*I_ERI_Gx3z_S_S_S_vrr+WPY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_S_S_vrr = PAX*I_ERI_G4z_S_S_S_vrr+WPX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_S_S_vrr = PAY*I_ERI_G4y_S_S_S_vrr+WPY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2z*I_ERI_F3y_S_S_S_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_S_S_vrr = PAZ*I_ERI_G4y_S_S_S_vrr+WPZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_S_S_vrr = PAZ*I_ERI_G3yz_S_S_S_vrr+WPZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_S_S_vrr = PAY*I_ERI_Gy3z_S_S_S_vrr+WPY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_S_S_vrr = PAY*I_ERI_G4z_S_S_S_vrr+WPY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_S_S_vrr = PAZ*I_ERI_G4z_S_S_S_vrr+WPZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2z*I_ERI_F3z_S_S_S_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_I6x_S_S_S_vrr = PAX*I_ERI_H5x_S_S_S_vrr+WPX*I_ERI_H5x_S_S_S_M1_vrr+5*oned2z*I_ERI_G4x_S_S_S_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_S_S_vrr = PAY*I_ERI_H5x_S_S_S_vrr+WPY*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_S_S_vrr = PAZ*I_ERI_H5x_S_S_S_vrr+WPZ*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_S_S_vrr = PAY*I_ERI_H4xy_S_S_S_vrr+WPY*I_ERI_H4xy_S_S_S_M1_vrr+oned2z*I_ERI_G4x_S_S_S_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_S_S_vrr = PAZ*I_ERI_H4xy_S_S_S_vrr+WPZ*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_S_S_vrr = PAZ*I_ERI_H4xz_S_S_S_vrr+WPZ*I_ERI_H4xz_S_S_S_M1_vrr+oned2z*I_ERI_G4x_S_S_S_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_S_S_vrr = PAY*I_ERI_H3x2y_S_S_S_vrr+WPY*I_ERI_H3x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_S_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_S_S_vrr = PAZ*I_ERI_H3x2y_S_S_S_vrr+WPZ*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_S_S_vrr = PAY*I_ERI_H3x2z_S_S_S_vrr+WPY*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_S_S_vrr = PAZ*I_ERI_H3x2z_S_S_S_vrr+WPZ*I_ERI_H3x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_S_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_S_S_vrr = PAX*I_ERI_Hx4y_S_S_S_vrr+WPX*I_ERI_Hx4y_S_S_S_M1_vrr+oned2z*I_ERI_G4y_S_S_S_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_S_S_vrr = PAZ*I_ERI_H2x3y_S_S_S_vrr+WPZ*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_S_S_vrr = PAZ*I_ERI_H2x2yz_S_S_S_vrr+WPZ*I_ERI_H2x2yz_S_S_S_M1_vrr+oned2z*I_ERI_G2x2y_S_S_S_vrr-rhod2zsq*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_S_S_vrr = PAY*I_ERI_H2x3z_S_S_S_vrr+WPY*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_S_S_vrr = PAX*I_ERI_Hx4z_S_S_S_vrr+WPX*I_ERI_Hx4z_S_S_S_M1_vrr+oned2z*I_ERI_G4z_S_S_S_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_S_S_vrr = PAX*I_ERI_H5y_S_S_S_vrr+WPX*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_S_S_vrr = PAZ*I_ERI_Hx4y_S_S_S_vrr+WPZ*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_S_S_vrr = PAX*I_ERI_H3y2z_S_S_S_vrr+WPX*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_S_S_vrr = PAX*I_ERI_H2y3z_S_S_S_vrr+WPX*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_S_S_vrr = PAY*I_ERI_Hx4z_S_S_S_vrr+WPY*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_S_S_vrr = PAX*I_ERI_H5z_S_S_S_vrr+WPX*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_S_S_vrr = PAY*I_ERI_H5y_S_S_S_vrr+WPY*I_ERI_H5y_S_S_S_M1_vrr+5*oned2z*I_ERI_G4y_S_S_S_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_S_S_vrr = PAZ*I_ERI_H5y_S_S_S_vrr+WPZ*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_S_S_vrr = PAZ*I_ERI_H4yz_S_S_S_vrr+WPZ*I_ERI_H4yz_S_S_S_M1_vrr+oned2z*I_ERI_G4y_S_S_S_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_S_S_vrr = PAZ*I_ERI_H3y2z_S_S_S_vrr+WPZ*I_ERI_H3y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_S_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_S_S_vrr = PAY*I_ERI_Hy4z_S_S_S_vrr+WPY*I_ERI_Hy4z_S_S_S_M1_vrr+oned2z*I_ERI_G4z_S_S_S_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_S_S_vrr = PAY*I_ERI_H5z_S_S_S_vrr+WPY*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_S_S_vrr = PAZ*I_ERI_H5z_S_S_S_vrr+WPZ*I_ERI_H5z_S_S_S_M1_vrr+5*oned2z*I_ERI_G4z_S_S_S_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       ************************************************************/
      Double I_ERI_K7x_S_S_S_vrr = PAX*I_ERI_I6x_S_S_S_vrr+WPX*I_ERI_I6x_S_S_S_M1_vrr+6*oned2z*I_ERI_H5x_S_S_S_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_S_S_vrr = PAY*I_ERI_I6x_S_S_S_vrr+WPY*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_S_S_vrr = PAZ*I_ERI_I6x_S_S_S_vrr+WPZ*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_S_S_vrr = PAY*I_ERI_I5xy_S_S_S_vrr+WPY*I_ERI_I5xy_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_S_S_vrr = PAZ*I_ERI_I5xy_S_S_S_vrr+WPZ*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_S_S_vrr = PAZ*I_ERI_I5xz_S_S_S_vrr+WPZ*I_ERI_I5xz_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_S_S_vrr = PAY*I_ERI_I4x2y_S_S_S_vrr+WPY*I_ERI_I4x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xy_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_S_S_vrr = PAZ*I_ERI_I4x2y_S_S_S_vrr+WPZ*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_S_S_vrr = PAY*I_ERI_I4x2z_S_S_S_vrr+WPY*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_S_S_vrr = PAZ*I_ERI_I4x2z_S_S_S_vrr+WPZ*I_ERI_I4x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_S_S_vrr = PAX*I_ERI_I2x4y_S_S_S_vrr+WPX*I_ERI_I2x4y_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_S_S_vrr = PAZ*I_ERI_I3x3y_S_S_S_vrr+WPZ*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_S_S_vrr = PAZ*I_ERI_I3x2yz_S_S_S_vrr+WPZ*I_ERI_I3x2yz_S_S_S_M1_vrr+oned2z*I_ERI_H3x2y_S_S_S_vrr-rhod2zsq*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_S_S_vrr = PAY*I_ERI_I3x3z_S_S_S_vrr+WPY*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_S_S_vrr = PAX*I_ERI_I2x4z_S_S_S_vrr+WPX*I_ERI_I2x4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_S_S_vrr = PAX*I_ERI_Ix5y_S_S_S_vrr+WPX*I_ERI_Ix5y_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_S_S_vrr = PAZ*I_ERI_I2x4y_S_S_S_vrr+WPZ*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_S_S_vrr = PAZ*I_ERI_I2x3yz_S_S_S_vrr+WPZ*I_ERI_I2x3yz_S_S_S_M1_vrr+oned2z*I_ERI_H2x3y_S_S_S_vrr-rhod2zsq*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_S_S_vrr = PAY*I_ERI_I2xy3z_S_S_S_vrr+WPY*I_ERI_I2xy3z_S_S_S_M1_vrr+oned2z*I_ERI_H2x3z_S_S_S_vrr-rhod2zsq*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_S_S_vrr = PAY*I_ERI_I2x4z_S_S_S_vrr+WPY*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_S_S_vrr = PAX*I_ERI_Ix5z_S_S_S_vrr+WPX*I_ERI_Ix5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_S_S_vrr = PAX*I_ERI_I6y_S_S_S_vrr+WPX*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_S_S_vrr = PAZ*I_ERI_Ix5y_S_S_S_vrr+WPZ*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_S_S_vrr = PAX*I_ERI_I4y2z_S_S_S_vrr+WPX*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_S_S_vrr = PAX*I_ERI_I3y3z_S_S_S_vrr+WPX*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_S_S_vrr = PAX*I_ERI_I2y4z_S_S_S_vrr+WPX*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_S_S_vrr = PAY*I_ERI_Ix5z_S_S_S_vrr+WPY*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_S_S_vrr = PAX*I_ERI_I6z_S_S_S_vrr+WPX*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_S_S_vrr = PAY*I_ERI_I6y_S_S_S_vrr+WPY*I_ERI_I6y_S_S_S_M1_vrr+6*oned2z*I_ERI_H5y_S_S_S_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_S_S_vrr = PAZ*I_ERI_I6y_S_S_S_vrr+WPZ*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_S_S_vrr = PAZ*I_ERI_I5yz_S_S_S_vrr+WPZ*I_ERI_I5yz_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_S_S_vrr = PAZ*I_ERI_I4y2z_S_S_S_vrr+WPZ*I_ERI_I4y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4yz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_S_S_vrr = PAY*I_ERI_I2y4z_S_S_S_vrr+WPY*I_ERI_I2y4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_S_S_vrr = PAY*I_ERI_Iy5z_S_S_S_vrr+WPY*I_ERI_Iy5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_S_S_vrr = PAY*I_ERI_I6z_S_S_S_vrr+WPY*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_S_S_vrr = PAZ*I_ERI_I6z_S_S_S_vrr+WPZ*I_ERI_I6z_S_S_S_M1_vrr+6*oned2z*I_ERI_H5z_S_S_S_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_K_S_S_S
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       ************************************************************/
      Double I_ERI_L8x_S_S_S_vrr = PAX*I_ERI_K7x_S_S_S_vrr+WPX*I_ERI_K7x_S_S_S_M1_vrr+7*oned2z*I_ERI_I6x_S_S_S_vrr-7*rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L7xy_S_S_S_vrr = PAY*I_ERI_K7x_S_S_S_vrr+WPY*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L7xz_S_S_S_vrr = PAZ*I_ERI_K7x_S_S_S_vrr+WPZ*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L6x2y_S_S_S_vrr = PAY*I_ERI_K6xy_S_S_S_vrr+WPY*I_ERI_K6xy_S_S_S_M1_vrr+oned2z*I_ERI_I6x_S_S_S_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L6xyz_S_S_S_vrr = PAZ*I_ERI_K6xy_S_S_S_vrr+WPZ*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_L6x2z_S_S_S_vrr = PAZ*I_ERI_K6xz_S_S_S_vrr+WPZ*I_ERI_K6xz_S_S_S_M1_vrr+oned2z*I_ERI_I6x_S_S_S_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L5x3y_S_S_S_vrr = PAY*I_ERI_K5x2y_S_S_S_vrr+WPY*I_ERI_K5x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_I5xy_S_S_S_vrr-2*rhod2zsq*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_L5x2yz_S_S_S_vrr = PAZ*I_ERI_K5x2y_S_S_S_vrr+WPZ*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_L5xy2z_S_S_S_vrr = PAY*I_ERI_K5x2z_S_S_S_vrr+WPY*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3z_S_S_S_vrr = PAZ*I_ERI_K5x2z_S_S_S_vrr+WPZ*I_ERI_K5x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_I5xz_S_S_S_vrr-2*rhod2zsq*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_L4x4y_S_S_S_vrr = PAY*I_ERI_K4x3y_S_S_S_vrr+WPY*I_ERI_K4x3y_S_S_S_M1_vrr+3*oned2z*I_ERI_I4x2y_S_S_S_vrr-3*rhod2zsq*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_L4x3yz_S_S_S_vrr = PAZ*I_ERI_K4x3y_S_S_S_vrr+WPZ*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_L4x2y2z_S_S_S_vrr = PAZ*I_ERI_K4x2yz_S_S_S_vrr+WPZ*I_ERI_K4x2yz_S_S_S_M1_vrr+oned2z*I_ERI_I4x2y_S_S_S_vrr-rhod2zsq*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_L4xy3z_S_S_S_vrr = PAY*I_ERI_K4x3z_S_S_S_vrr+WPY*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4z_S_S_S_vrr = PAZ*I_ERI_K4x3z_S_S_S_vrr+WPZ*I_ERI_K4x3z_S_S_S_M1_vrr+3*oned2z*I_ERI_I4x2z_S_S_S_vrr-3*rhod2zsq*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_L3x5y_S_S_S_vrr = PAX*I_ERI_K2x5y_S_S_S_vrr+WPX*I_ERI_K2x5y_S_S_S_M1_vrr+2*oned2z*I_ERI_Ix5y_S_S_S_vrr-2*rhod2zsq*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_L3x4yz_S_S_S_vrr = PAZ*I_ERI_K3x4y_S_S_S_vrr+WPZ*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_L3x3y2z_S_S_S_vrr = PAZ*I_ERI_K3x3yz_S_S_S_vrr+WPZ*I_ERI_K3x3yz_S_S_S_M1_vrr+oned2z*I_ERI_I3x3y_S_S_S_vrr-rhod2zsq*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_L3x2y3z_S_S_S_vrr = PAY*I_ERI_K3xy3z_S_S_S_vrr+WPY*I_ERI_K3xy3z_S_S_S_M1_vrr+oned2z*I_ERI_I3x3z_S_S_S_vrr-rhod2zsq*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_L3xy4z_S_S_S_vrr = PAY*I_ERI_K3x4z_S_S_S_vrr+WPY*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5z_S_S_S_vrr = PAX*I_ERI_K2x5z_S_S_S_vrr+WPX*I_ERI_K2x5z_S_S_S_M1_vrr+2*oned2z*I_ERI_Ix5z_S_S_S_vrr-2*rhod2zsq*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6y_S_S_S_vrr = PAX*I_ERI_Kx6y_S_S_S_vrr+WPX*I_ERI_Kx6y_S_S_S_M1_vrr+oned2z*I_ERI_I6y_S_S_S_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L2x5yz_S_S_S_vrr = PAZ*I_ERI_K2x5y_S_S_S_vrr+WPZ*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_L2x4y2z_S_S_S_vrr = PAZ*I_ERI_K2x4yz_S_S_S_vrr+WPZ*I_ERI_K2x4yz_S_S_S_M1_vrr+oned2z*I_ERI_I2x4y_S_S_S_vrr-rhod2zsq*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_L2x3y3z_S_S_S_vrr = PAX*I_ERI_Kx3y3z_S_S_S_vrr+WPX*I_ERI_Kx3y3z_S_S_S_M1_vrr+oned2z*I_ERI_I3y3z_S_S_S_vrr-rhod2zsq*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_L2x2y4z_S_S_S_vrr = PAY*I_ERI_K2xy4z_S_S_S_vrr+WPY*I_ERI_K2xy4z_S_S_S_M1_vrr+oned2z*I_ERI_I2x4z_S_S_S_vrr-rhod2zsq*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_L2xy5z_S_S_S_vrr = PAY*I_ERI_K2x5z_S_S_S_vrr+WPY*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6z_S_S_S_vrr = PAX*I_ERI_Kx6z_S_S_S_vrr+WPX*I_ERI_Kx6z_S_S_S_M1_vrr+oned2z*I_ERI_I6z_S_S_S_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7y_S_S_S_vrr = PAX*I_ERI_K7y_S_S_S_vrr+WPX*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_Lx6yz_S_S_S_vrr = PAZ*I_ERI_Kx6y_S_S_S_vrr+WPZ*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_Lx5y2z_S_S_S_vrr = PAX*I_ERI_K5y2z_S_S_S_vrr+WPX*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_Lx4y3z_S_S_S_vrr = PAX*I_ERI_K4y3z_S_S_S_vrr+WPX*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_Lx3y4z_S_S_S_vrr = PAX*I_ERI_K3y4z_S_S_S_vrr+WPX*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_Lx2y5z_S_S_S_vrr = PAX*I_ERI_K2y5z_S_S_S_vrr+WPX*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_Lxy6z_S_S_S_vrr = PAY*I_ERI_Kx6z_S_S_S_vrr+WPY*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7z_S_S_S_vrr = PAX*I_ERI_K7z_S_S_S_vrr+WPX*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8y_S_S_S_vrr = PAY*I_ERI_K7y_S_S_S_vrr+WPY*I_ERI_K7y_S_S_S_M1_vrr+7*oned2z*I_ERI_I6y_S_S_S_vrr-7*rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L7yz_S_S_S_vrr = PAZ*I_ERI_K7y_S_S_S_vrr+WPZ*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_L6y2z_S_S_S_vrr = PAZ*I_ERI_K6yz_S_S_S_vrr+WPZ*I_ERI_K6yz_S_S_S_M1_vrr+oned2z*I_ERI_I6y_S_S_S_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L5y3z_S_S_S_vrr = PAZ*I_ERI_K5y2z_S_S_S_vrr+WPZ*I_ERI_K5y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_I5yz_S_S_S_vrr-2*rhod2zsq*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_L4y4z_S_S_S_vrr = PAZ*I_ERI_K4y3z_S_S_S_vrr+WPZ*I_ERI_K4y3z_S_S_S_M1_vrr+3*oned2z*I_ERI_I4y2z_S_S_S_vrr-3*rhod2zsq*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_L3y5z_S_S_S_vrr = PAY*I_ERI_K2y5z_S_S_S_vrr+WPY*I_ERI_K2y5z_S_S_S_M1_vrr+2*oned2z*I_ERI_Iy5z_S_S_S_vrr-2*rhod2zsq*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_L2y6z_S_S_S_vrr = PAY*I_ERI_Ky6z_S_S_S_vrr+WPY*I_ERI_Ky6z_S_S_S_M1_vrr+oned2z*I_ERI_I6z_S_S_S_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_Ly7z_S_S_S_vrr = PAY*I_ERI_K7z_S_S_S_vrr+WPY*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8z_S_S_S_vrr = PAZ*I_ERI_K7z_S_S_S_vrr+WPZ*I_ERI_K7z_S_S_S_M1_vrr+7*oned2z*I_ERI_I6z_S_S_S_vrr-7*rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_L8x_S_S_S += I_ERI_L8x_S_S_S_vrr;
      I_ERI_L7xy_S_S_S += I_ERI_L7xy_S_S_S_vrr;
      I_ERI_L7xz_S_S_S += I_ERI_L7xz_S_S_S_vrr;
      I_ERI_L6x2y_S_S_S += I_ERI_L6x2y_S_S_S_vrr;
      I_ERI_L6xyz_S_S_S += I_ERI_L6xyz_S_S_S_vrr;
      I_ERI_L6x2z_S_S_S += I_ERI_L6x2z_S_S_S_vrr;
      I_ERI_L5x3y_S_S_S += I_ERI_L5x3y_S_S_S_vrr;
      I_ERI_L5x2yz_S_S_S += I_ERI_L5x2yz_S_S_S_vrr;
      I_ERI_L5xy2z_S_S_S += I_ERI_L5xy2z_S_S_S_vrr;
      I_ERI_L5x3z_S_S_S += I_ERI_L5x3z_S_S_S_vrr;
      I_ERI_L4x4y_S_S_S += I_ERI_L4x4y_S_S_S_vrr;
      I_ERI_L4x3yz_S_S_S += I_ERI_L4x3yz_S_S_S_vrr;
      I_ERI_L4x2y2z_S_S_S += I_ERI_L4x2y2z_S_S_S_vrr;
      I_ERI_L4xy3z_S_S_S += I_ERI_L4xy3z_S_S_S_vrr;
      I_ERI_L4x4z_S_S_S += I_ERI_L4x4z_S_S_S_vrr;
      I_ERI_L3x5y_S_S_S += I_ERI_L3x5y_S_S_S_vrr;
      I_ERI_L3x4yz_S_S_S += I_ERI_L3x4yz_S_S_S_vrr;
      I_ERI_L3x3y2z_S_S_S += I_ERI_L3x3y2z_S_S_S_vrr;
      I_ERI_L3x2y3z_S_S_S += I_ERI_L3x2y3z_S_S_S_vrr;
      I_ERI_L3xy4z_S_S_S += I_ERI_L3xy4z_S_S_S_vrr;
      I_ERI_L3x5z_S_S_S += I_ERI_L3x5z_S_S_S_vrr;
      I_ERI_L2x6y_S_S_S += I_ERI_L2x6y_S_S_S_vrr;
      I_ERI_L2x5yz_S_S_S += I_ERI_L2x5yz_S_S_S_vrr;
      I_ERI_L2x4y2z_S_S_S += I_ERI_L2x4y2z_S_S_S_vrr;
      I_ERI_L2x3y3z_S_S_S += I_ERI_L2x3y3z_S_S_S_vrr;
      I_ERI_L2x2y4z_S_S_S += I_ERI_L2x2y4z_S_S_S_vrr;
      I_ERI_L2xy5z_S_S_S += I_ERI_L2xy5z_S_S_S_vrr;
      I_ERI_L2x6z_S_S_S += I_ERI_L2x6z_S_S_S_vrr;
      I_ERI_Lx7y_S_S_S += I_ERI_Lx7y_S_S_S_vrr;
      I_ERI_Lx6yz_S_S_S += I_ERI_Lx6yz_S_S_S_vrr;
      I_ERI_Lx5y2z_S_S_S += I_ERI_Lx5y2z_S_S_S_vrr;
      I_ERI_Lx4y3z_S_S_S += I_ERI_Lx4y3z_S_S_S_vrr;
      I_ERI_Lx3y4z_S_S_S += I_ERI_Lx3y4z_S_S_S_vrr;
      I_ERI_Lx2y5z_S_S_S += I_ERI_Lx2y5z_S_S_S_vrr;
      I_ERI_Lxy6z_S_S_S += I_ERI_Lxy6z_S_S_S_vrr;
      I_ERI_Lx7z_S_S_S += I_ERI_Lx7z_S_S_S_vrr;
      I_ERI_L8y_S_S_S += I_ERI_L8y_S_S_S_vrr;
      I_ERI_L7yz_S_S_S += I_ERI_L7yz_S_S_S_vrr;
      I_ERI_L6y2z_S_S_S += I_ERI_L6y2z_S_S_S_vrr;
      I_ERI_L5y3z_S_S_S += I_ERI_L5y3z_S_S_S_vrr;
      I_ERI_L4y4z_S_S_S += I_ERI_L4y4z_S_S_S_vrr;
      I_ERI_L3y5z_S_S_S += I_ERI_L3y5z_S_S_S_vrr;
      I_ERI_L2y6z_S_S_S += I_ERI_L2y6z_S_S_S_vrr;
      I_ERI_Ly7z_S_S_S += I_ERI_Ly7z_S_S_S_vrr;
      I_ERI_L8z_S_S_S += I_ERI_L8z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_K7x_S_S_S += I_ERI_K7x_S_S_S_vrr;
      I_ERI_K6xy_S_S_S += I_ERI_K6xy_S_S_S_vrr;
      I_ERI_K6xz_S_S_S += I_ERI_K6xz_S_S_S_vrr;
      I_ERI_K5x2y_S_S_S += I_ERI_K5x2y_S_S_S_vrr;
      I_ERI_K5xyz_S_S_S += I_ERI_K5xyz_S_S_S_vrr;
      I_ERI_K5x2z_S_S_S += I_ERI_K5x2z_S_S_S_vrr;
      I_ERI_K4x3y_S_S_S += I_ERI_K4x3y_S_S_S_vrr;
      I_ERI_K4x2yz_S_S_S += I_ERI_K4x2yz_S_S_S_vrr;
      I_ERI_K4xy2z_S_S_S += I_ERI_K4xy2z_S_S_S_vrr;
      I_ERI_K4x3z_S_S_S += I_ERI_K4x3z_S_S_S_vrr;
      I_ERI_K3x4y_S_S_S += I_ERI_K3x4y_S_S_S_vrr;
      I_ERI_K3x3yz_S_S_S += I_ERI_K3x3yz_S_S_S_vrr;
      I_ERI_K3x2y2z_S_S_S += I_ERI_K3x2y2z_S_S_S_vrr;
      I_ERI_K3xy3z_S_S_S += I_ERI_K3xy3z_S_S_S_vrr;
      I_ERI_K3x4z_S_S_S += I_ERI_K3x4z_S_S_S_vrr;
      I_ERI_K2x5y_S_S_S += I_ERI_K2x5y_S_S_S_vrr;
      I_ERI_K2x4yz_S_S_S += I_ERI_K2x4yz_S_S_S_vrr;
      I_ERI_K2x3y2z_S_S_S += I_ERI_K2x3y2z_S_S_S_vrr;
      I_ERI_K2x2y3z_S_S_S += I_ERI_K2x2y3z_S_S_S_vrr;
      I_ERI_K2xy4z_S_S_S += I_ERI_K2xy4z_S_S_S_vrr;
      I_ERI_K2x5z_S_S_S += I_ERI_K2x5z_S_S_S_vrr;
      I_ERI_Kx6y_S_S_S += I_ERI_Kx6y_S_S_S_vrr;
      I_ERI_Kx5yz_S_S_S += I_ERI_Kx5yz_S_S_S_vrr;
      I_ERI_Kx4y2z_S_S_S += I_ERI_Kx4y2z_S_S_S_vrr;
      I_ERI_Kx3y3z_S_S_S += I_ERI_Kx3y3z_S_S_S_vrr;
      I_ERI_Kx2y4z_S_S_S += I_ERI_Kx2y4z_S_S_S_vrr;
      I_ERI_Kxy5z_S_S_S += I_ERI_Kxy5z_S_S_S_vrr;
      I_ERI_Kx6z_S_S_S += I_ERI_Kx6z_S_S_S_vrr;
      I_ERI_K7y_S_S_S += I_ERI_K7y_S_S_S_vrr;
      I_ERI_K6yz_S_S_S += I_ERI_K6yz_S_S_S_vrr;
      I_ERI_K5y2z_S_S_S += I_ERI_K5y2z_S_S_S_vrr;
      I_ERI_K4y3z_S_S_S += I_ERI_K4y3z_S_S_S_vrr;
      I_ERI_K3y4z_S_S_S += I_ERI_K3y4z_S_S_S_vrr;
      I_ERI_K2y5z_S_S_S += I_ERI_K2y5z_S_S_S_vrr;
      I_ERI_Ky6z_S_S_S += I_ERI_Ky6z_S_S_S_vrr;
      I_ERI_K7z_S_S_S += I_ERI_K7z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_I6x_S_S_S += I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S += I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S += I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S += I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S += I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S += I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S += I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S += I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S += I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S += I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S += I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S += I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S += I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S += I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S += I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S += I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S += I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S += I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S += I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S += I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S += I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S += I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S += I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S += I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S += I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S += I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S += I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S += I_ERI_I6z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_H5x_S_S_S += I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S += I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S += I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S += I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S += I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S += I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S += I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S += I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S += I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S += I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S += I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S += I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S += I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S += I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S += I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S += I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S += I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S += I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S += I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S += I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S += I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_S_S += I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S += I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S += I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S += I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S += I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S += I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S += I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S += I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S += I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S += I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S += I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S += I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S += I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S += I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S += I_ERI_G4z_S_S_S_vrr;
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
   * shell quartet name: SQ_ERI_G_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S
   * RHS shell quartet name: SQ_ERI_G_S_S_S
   ************************************************************/
  Double I_ERI_G4x_Px_S_S = I_ERI_H5x_S_S_S+ABX*I_ERI_G4x_S_S_S;
  Double I_ERI_G3xy_Px_S_S = I_ERI_H4xy_S_S_S+ABX*I_ERI_G3xy_S_S_S;
  Double I_ERI_G3xz_Px_S_S = I_ERI_H4xz_S_S_S+ABX*I_ERI_G3xz_S_S_S;
  Double I_ERI_G2x2y_Px_S_S = I_ERI_H3x2y_S_S_S+ABX*I_ERI_G2x2y_S_S_S;
  Double I_ERI_G2xyz_Px_S_S = I_ERI_H3xyz_S_S_S+ABX*I_ERI_G2xyz_S_S_S;
  Double I_ERI_G2x2z_Px_S_S = I_ERI_H3x2z_S_S_S+ABX*I_ERI_G2x2z_S_S_S;
  Double I_ERI_Gx3y_Px_S_S = I_ERI_H2x3y_S_S_S+ABX*I_ERI_Gx3y_S_S_S;
  Double I_ERI_Gx2yz_Px_S_S = I_ERI_H2x2yz_S_S_S+ABX*I_ERI_Gx2yz_S_S_S;
  Double I_ERI_Gxy2z_Px_S_S = I_ERI_H2xy2z_S_S_S+ABX*I_ERI_Gxy2z_S_S_S;
  Double I_ERI_Gx3z_Px_S_S = I_ERI_H2x3z_S_S_S+ABX*I_ERI_Gx3z_S_S_S;
  Double I_ERI_G4y_Px_S_S = I_ERI_Hx4y_S_S_S+ABX*I_ERI_G4y_S_S_S;
  Double I_ERI_G3yz_Px_S_S = I_ERI_Hx3yz_S_S_S+ABX*I_ERI_G3yz_S_S_S;
  Double I_ERI_G2y2z_Px_S_S = I_ERI_Hx2y2z_S_S_S+ABX*I_ERI_G2y2z_S_S_S;
  Double I_ERI_Gy3z_Px_S_S = I_ERI_Hxy3z_S_S_S+ABX*I_ERI_Gy3z_S_S_S;
  Double I_ERI_G4z_Px_S_S = I_ERI_Hx4z_S_S_S+ABX*I_ERI_G4z_S_S_S;
  Double I_ERI_G4x_Py_S_S = I_ERI_H4xy_S_S_S+ABY*I_ERI_G4x_S_S_S;
  Double I_ERI_G3xy_Py_S_S = I_ERI_H3x2y_S_S_S+ABY*I_ERI_G3xy_S_S_S;
  Double I_ERI_G3xz_Py_S_S = I_ERI_H3xyz_S_S_S+ABY*I_ERI_G3xz_S_S_S;
  Double I_ERI_G2x2y_Py_S_S = I_ERI_H2x3y_S_S_S+ABY*I_ERI_G2x2y_S_S_S;
  Double I_ERI_G2xyz_Py_S_S = I_ERI_H2x2yz_S_S_S+ABY*I_ERI_G2xyz_S_S_S;
  Double I_ERI_G2x2z_Py_S_S = I_ERI_H2xy2z_S_S_S+ABY*I_ERI_G2x2z_S_S_S;
  Double I_ERI_Gx3y_Py_S_S = I_ERI_Hx4y_S_S_S+ABY*I_ERI_Gx3y_S_S_S;
  Double I_ERI_Gx2yz_Py_S_S = I_ERI_Hx3yz_S_S_S+ABY*I_ERI_Gx2yz_S_S_S;
  Double I_ERI_Gxy2z_Py_S_S = I_ERI_Hx2y2z_S_S_S+ABY*I_ERI_Gxy2z_S_S_S;
  Double I_ERI_Gx3z_Py_S_S = I_ERI_Hxy3z_S_S_S+ABY*I_ERI_Gx3z_S_S_S;
  Double I_ERI_G4y_Py_S_S = I_ERI_H5y_S_S_S+ABY*I_ERI_G4y_S_S_S;
  Double I_ERI_G3yz_Py_S_S = I_ERI_H4yz_S_S_S+ABY*I_ERI_G3yz_S_S_S;
  Double I_ERI_G2y2z_Py_S_S = I_ERI_H3y2z_S_S_S+ABY*I_ERI_G2y2z_S_S_S;
  Double I_ERI_Gy3z_Py_S_S = I_ERI_H2y3z_S_S_S+ABY*I_ERI_Gy3z_S_S_S;
  Double I_ERI_G4z_Py_S_S = I_ERI_Hy4z_S_S_S+ABY*I_ERI_G4z_S_S_S;
  Double I_ERI_G4x_Pz_S_S = I_ERI_H4xz_S_S_S+ABZ*I_ERI_G4x_S_S_S;
  Double I_ERI_G3xy_Pz_S_S = I_ERI_H3xyz_S_S_S+ABZ*I_ERI_G3xy_S_S_S;
  Double I_ERI_G3xz_Pz_S_S = I_ERI_H3x2z_S_S_S+ABZ*I_ERI_G3xz_S_S_S;
  Double I_ERI_G2x2y_Pz_S_S = I_ERI_H2x2yz_S_S_S+ABZ*I_ERI_G2x2y_S_S_S;
  Double I_ERI_G2xyz_Pz_S_S = I_ERI_H2xy2z_S_S_S+ABZ*I_ERI_G2xyz_S_S_S;
  Double I_ERI_G2x2z_Pz_S_S = I_ERI_H2x3z_S_S_S+ABZ*I_ERI_G2x2z_S_S_S;
  Double I_ERI_Gx3y_Pz_S_S = I_ERI_Hx3yz_S_S_S+ABZ*I_ERI_Gx3y_S_S_S;
  Double I_ERI_Gx2yz_Pz_S_S = I_ERI_Hx2y2z_S_S_S+ABZ*I_ERI_Gx2yz_S_S_S;
  Double I_ERI_Gxy2z_Pz_S_S = I_ERI_Hxy3z_S_S_S+ABZ*I_ERI_Gxy2z_S_S_S;
  Double I_ERI_Gx3z_Pz_S_S = I_ERI_Hx4z_S_S_S+ABZ*I_ERI_Gx3z_S_S_S;
  Double I_ERI_G4y_Pz_S_S = I_ERI_H4yz_S_S_S+ABZ*I_ERI_G4y_S_S_S;
  Double I_ERI_G3yz_Pz_S_S = I_ERI_H3y2z_S_S_S+ABZ*I_ERI_G3yz_S_S_S;
  Double I_ERI_G2y2z_Pz_S_S = I_ERI_H2y3z_S_S_S+ABZ*I_ERI_G2y2z_S_S_S;
  Double I_ERI_Gy3z_Pz_S_S = I_ERI_Hy4z_S_S_S+ABZ*I_ERI_Gy3z_S_S_S;
  Double I_ERI_G4z_Pz_S_S = I_ERI_H5z_S_S_S+ABZ*I_ERI_G4z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S
   * RHS shell quartet name: SQ_ERI_H_S_S_S
   ************************************************************/
  Double I_ERI_H5x_Px_S_S = I_ERI_I6x_S_S_S+ABX*I_ERI_H5x_S_S_S;
  Double I_ERI_H4xy_Px_S_S = I_ERI_I5xy_S_S_S+ABX*I_ERI_H4xy_S_S_S;
  Double I_ERI_H4xz_Px_S_S = I_ERI_I5xz_S_S_S+ABX*I_ERI_H4xz_S_S_S;
  Double I_ERI_H3x2y_Px_S_S = I_ERI_I4x2y_S_S_S+ABX*I_ERI_H3x2y_S_S_S;
  Double I_ERI_H3xyz_Px_S_S = I_ERI_I4xyz_S_S_S+ABX*I_ERI_H3xyz_S_S_S;
  Double I_ERI_H3x2z_Px_S_S = I_ERI_I4x2z_S_S_S+ABX*I_ERI_H3x2z_S_S_S;
  Double I_ERI_H2x3y_Px_S_S = I_ERI_I3x3y_S_S_S+ABX*I_ERI_H2x3y_S_S_S;
  Double I_ERI_H2x2yz_Px_S_S = I_ERI_I3x2yz_S_S_S+ABX*I_ERI_H2x2yz_S_S_S;
  Double I_ERI_H2xy2z_Px_S_S = I_ERI_I3xy2z_S_S_S+ABX*I_ERI_H2xy2z_S_S_S;
  Double I_ERI_H2x3z_Px_S_S = I_ERI_I3x3z_S_S_S+ABX*I_ERI_H2x3z_S_S_S;
  Double I_ERI_Hx4y_Px_S_S = I_ERI_I2x4y_S_S_S+ABX*I_ERI_Hx4y_S_S_S;
  Double I_ERI_Hx3yz_Px_S_S = I_ERI_I2x3yz_S_S_S+ABX*I_ERI_Hx3yz_S_S_S;
  Double I_ERI_Hx2y2z_Px_S_S = I_ERI_I2x2y2z_S_S_S+ABX*I_ERI_Hx2y2z_S_S_S;
  Double I_ERI_Hxy3z_Px_S_S = I_ERI_I2xy3z_S_S_S+ABX*I_ERI_Hxy3z_S_S_S;
  Double I_ERI_Hx4z_Px_S_S = I_ERI_I2x4z_S_S_S+ABX*I_ERI_Hx4z_S_S_S;
  Double I_ERI_H5y_Px_S_S = I_ERI_Ix5y_S_S_S+ABX*I_ERI_H5y_S_S_S;
  Double I_ERI_H4yz_Px_S_S = I_ERI_Ix4yz_S_S_S+ABX*I_ERI_H4yz_S_S_S;
  Double I_ERI_H3y2z_Px_S_S = I_ERI_Ix3y2z_S_S_S+ABX*I_ERI_H3y2z_S_S_S;
  Double I_ERI_H2y3z_Px_S_S = I_ERI_Ix2y3z_S_S_S+ABX*I_ERI_H2y3z_S_S_S;
  Double I_ERI_Hy4z_Px_S_S = I_ERI_Ixy4z_S_S_S+ABX*I_ERI_Hy4z_S_S_S;
  Double I_ERI_H5z_Px_S_S = I_ERI_Ix5z_S_S_S+ABX*I_ERI_H5z_S_S_S;
  Double I_ERI_H5x_Py_S_S = I_ERI_I5xy_S_S_S+ABY*I_ERI_H5x_S_S_S;
  Double I_ERI_H4xy_Py_S_S = I_ERI_I4x2y_S_S_S+ABY*I_ERI_H4xy_S_S_S;
  Double I_ERI_H4xz_Py_S_S = I_ERI_I4xyz_S_S_S+ABY*I_ERI_H4xz_S_S_S;
  Double I_ERI_H3x2y_Py_S_S = I_ERI_I3x3y_S_S_S+ABY*I_ERI_H3x2y_S_S_S;
  Double I_ERI_H3xyz_Py_S_S = I_ERI_I3x2yz_S_S_S+ABY*I_ERI_H3xyz_S_S_S;
  Double I_ERI_H3x2z_Py_S_S = I_ERI_I3xy2z_S_S_S+ABY*I_ERI_H3x2z_S_S_S;
  Double I_ERI_H2x3y_Py_S_S = I_ERI_I2x4y_S_S_S+ABY*I_ERI_H2x3y_S_S_S;
  Double I_ERI_H2x2yz_Py_S_S = I_ERI_I2x3yz_S_S_S+ABY*I_ERI_H2x2yz_S_S_S;
  Double I_ERI_H2xy2z_Py_S_S = I_ERI_I2x2y2z_S_S_S+ABY*I_ERI_H2xy2z_S_S_S;
  Double I_ERI_H2x3z_Py_S_S = I_ERI_I2xy3z_S_S_S+ABY*I_ERI_H2x3z_S_S_S;
  Double I_ERI_Hx4y_Py_S_S = I_ERI_Ix5y_S_S_S+ABY*I_ERI_Hx4y_S_S_S;
  Double I_ERI_Hx3yz_Py_S_S = I_ERI_Ix4yz_S_S_S+ABY*I_ERI_Hx3yz_S_S_S;
  Double I_ERI_Hx2y2z_Py_S_S = I_ERI_Ix3y2z_S_S_S+ABY*I_ERI_Hx2y2z_S_S_S;
  Double I_ERI_Hxy3z_Py_S_S = I_ERI_Ix2y3z_S_S_S+ABY*I_ERI_Hxy3z_S_S_S;
  Double I_ERI_Hx4z_Py_S_S = I_ERI_Ixy4z_S_S_S+ABY*I_ERI_Hx4z_S_S_S;
  Double I_ERI_H5y_Py_S_S = I_ERI_I6y_S_S_S+ABY*I_ERI_H5y_S_S_S;
  Double I_ERI_H4yz_Py_S_S = I_ERI_I5yz_S_S_S+ABY*I_ERI_H4yz_S_S_S;
  Double I_ERI_H3y2z_Py_S_S = I_ERI_I4y2z_S_S_S+ABY*I_ERI_H3y2z_S_S_S;
  Double I_ERI_H2y3z_Py_S_S = I_ERI_I3y3z_S_S_S+ABY*I_ERI_H2y3z_S_S_S;
  Double I_ERI_Hy4z_Py_S_S = I_ERI_I2y4z_S_S_S+ABY*I_ERI_Hy4z_S_S_S;
  Double I_ERI_H5z_Py_S_S = I_ERI_Iy5z_S_S_S+ABY*I_ERI_H5z_S_S_S;
  Double I_ERI_H5x_Pz_S_S = I_ERI_I5xz_S_S_S+ABZ*I_ERI_H5x_S_S_S;
  Double I_ERI_H4xy_Pz_S_S = I_ERI_I4xyz_S_S_S+ABZ*I_ERI_H4xy_S_S_S;
  Double I_ERI_H4xz_Pz_S_S = I_ERI_I4x2z_S_S_S+ABZ*I_ERI_H4xz_S_S_S;
  Double I_ERI_H3x2y_Pz_S_S = I_ERI_I3x2yz_S_S_S+ABZ*I_ERI_H3x2y_S_S_S;
  Double I_ERI_H3xyz_Pz_S_S = I_ERI_I3xy2z_S_S_S+ABZ*I_ERI_H3xyz_S_S_S;
  Double I_ERI_H3x2z_Pz_S_S = I_ERI_I3x3z_S_S_S+ABZ*I_ERI_H3x2z_S_S_S;
  Double I_ERI_H2x3y_Pz_S_S = I_ERI_I2x3yz_S_S_S+ABZ*I_ERI_H2x3y_S_S_S;
  Double I_ERI_H2x2yz_Pz_S_S = I_ERI_I2x2y2z_S_S_S+ABZ*I_ERI_H2x2yz_S_S_S;
  Double I_ERI_H2xy2z_Pz_S_S = I_ERI_I2xy3z_S_S_S+ABZ*I_ERI_H2xy2z_S_S_S;
  Double I_ERI_H2x3z_Pz_S_S = I_ERI_I2x4z_S_S_S+ABZ*I_ERI_H2x3z_S_S_S;
  Double I_ERI_Hx4y_Pz_S_S = I_ERI_Ix4yz_S_S_S+ABZ*I_ERI_Hx4y_S_S_S;
  Double I_ERI_Hx3yz_Pz_S_S = I_ERI_Ix3y2z_S_S_S+ABZ*I_ERI_Hx3yz_S_S_S;
  Double I_ERI_Hx2y2z_Pz_S_S = I_ERI_Ix2y3z_S_S_S+ABZ*I_ERI_Hx2y2z_S_S_S;
  Double I_ERI_Hxy3z_Pz_S_S = I_ERI_Ixy4z_S_S_S+ABZ*I_ERI_Hxy3z_S_S_S;
  Double I_ERI_Hx4z_Pz_S_S = I_ERI_Ix5z_S_S_S+ABZ*I_ERI_Hx4z_S_S_S;
  Double I_ERI_H5y_Pz_S_S = I_ERI_I5yz_S_S_S+ABZ*I_ERI_H5y_S_S_S;
  Double I_ERI_H4yz_Pz_S_S = I_ERI_I4y2z_S_S_S+ABZ*I_ERI_H4yz_S_S_S;
  Double I_ERI_H3y2z_Pz_S_S = I_ERI_I3y3z_S_S_S+ABZ*I_ERI_H3y2z_S_S_S;
  Double I_ERI_H2y3z_Pz_S_S = I_ERI_I2y4z_S_S_S+ABZ*I_ERI_H2y3z_S_S_S;
  Double I_ERI_Hy4z_Pz_S_S = I_ERI_Iy5z_S_S_S+ABZ*I_ERI_Hy4z_S_S_S;
  Double I_ERI_H5z_Pz_S_S = I_ERI_I6z_S_S_S+ABZ*I_ERI_H5z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_S_S
   * RHS shell quartet name: SQ_ERI_G_P_S_S
   ************************************************************/
  Double I_ERI_G4x_D2x_S_S = I_ERI_H5x_Px_S_S+ABX*I_ERI_G4x_Px_S_S;
  Double I_ERI_G3xy_D2x_S_S = I_ERI_H4xy_Px_S_S+ABX*I_ERI_G3xy_Px_S_S;
  Double I_ERI_G3xz_D2x_S_S = I_ERI_H4xz_Px_S_S+ABX*I_ERI_G3xz_Px_S_S;
  Double I_ERI_G2x2y_D2x_S_S = I_ERI_H3x2y_Px_S_S+ABX*I_ERI_G2x2y_Px_S_S;
  Double I_ERI_G2xyz_D2x_S_S = I_ERI_H3xyz_Px_S_S+ABX*I_ERI_G2xyz_Px_S_S;
  Double I_ERI_G2x2z_D2x_S_S = I_ERI_H3x2z_Px_S_S+ABX*I_ERI_G2x2z_Px_S_S;
  Double I_ERI_Gx3y_D2x_S_S = I_ERI_H2x3y_Px_S_S+ABX*I_ERI_Gx3y_Px_S_S;
  Double I_ERI_Gx2yz_D2x_S_S = I_ERI_H2x2yz_Px_S_S+ABX*I_ERI_Gx2yz_Px_S_S;
  Double I_ERI_Gxy2z_D2x_S_S = I_ERI_H2xy2z_Px_S_S+ABX*I_ERI_Gxy2z_Px_S_S;
  Double I_ERI_Gx3z_D2x_S_S = I_ERI_H2x3z_Px_S_S+ABX*I_ERI_Gx3z_Px_S_S;
  Double I_ERI_G4y_D2x_S_S = I_ERI_Hx4y_Px_S_S+ABX*I_ERI_G4y_Px_S_S;
  Double I_ERI_G3yz_D2x_S_S = I_ERI_Hx3yz_Px_S_S+ABX*I_ERI_G3yz_Px_S_S;
  Double I_ERI_G2y2z_D2x_S_S = I_ERI_Hx2y2z_Px_S_S+ABX*I_ERI_G2y2z_Px_S_S;
  Double I_ERI_Gy3z_D2x_S_S = I_ERI_Hxy3z_Px_S_S+ABX*I_ERI_Gy3z_Px_S_S;
  Double I_ERI_G4z_D2x_S_S = I_ERI_Hx4z_Px_S_S+ABX*I_ERI_G4z_Px_S_S;
  Double I_ERI_G4x_D2y_S_S = I_ERI_H4xy_Py_S_S+ABY*I_ERI_G4x_Py_S_S;
  Double I_ERI_G3xy_D2y_S_S = I_ERI_H3x2y_Py_S_S+ABY*I_ERI_G3xy_Py_S_S;
  Double I_ERI_G3xz_D2y_S_S = I_ERI_H3xyz_Py_S_S+ABY*I_ERI_G3xz_Py_S_S;
  Double I_ERI_G2x2y_D2y_S_S = I_ERI_H2x3y_Py_S_S+ABY*I_ERI_G2x2y_Py_S_S;
  Double I_ERI_G2xyz_D2y_S_S = I_ERI_H2x2yz_Py_S_S+ABY*I_ERI_G2xyz_Py_S_S;
  Double I_ERI_G2x2z_D2y_S_S = I_ERI_H2xy2z_Py_S_S+ABY*I_ERI_G2x2z_Py_S_S;
  Double I_ERI_Gx3y_D2y_S_S = I_ERI_Hx4y_Py_S_S+ABY*I_ERI_Gx3y_Py_S_S;
  Double I_ERI_Gx2yz_D2y_S_S = I_ERI_Hx3yz_Py_S_S+ABY*I_ERI_Gx2yz_Py_S_S;
  Double I_ERI_Gxy2z_D2y_S_S = I_ERI_Hx2y2z_Py_S_S+ABY*I_ERI_Gxy2z_Py_S_S;
  Double I_ERI_Gx3z_D2y_S_S = I_ERI_Hxy3z_Py_S_S+ABY*I_ERI_Gx3z_Py_S_S;
  Double I_ERI_G4y_D2y_S_S = I_ERI_H5y_Py_S_S+ABY*I_ERI_G4y_Py_S_S;
  Double I_ERI_G3yz_D2y_S_S = I_ERI_H4yz_Py_S_S+ABY*I_ERI_G3yz_Py_S_S;
  Double I_ERI_G2y2z_D2y_S_S = I_ERI_H3y2z_Py_S_S+ABY*I_ERI_G2y2z_Py_S_S;
  Double I_ERI_Gy3z_D2y_S_S = I_ERI_H2y3z_Py_S_S+ABY*I_ERI_Gy3z_Py_S_S;
  Double I_ERI_G4z_D2y_S_S = I_ERI_Hy4z_Py_S_S+ABY*I_ERI_G4z_Py_S_S;
  Double I_ERI_G4x_D2z_S_S = I_ERI_H4xz_Pz_S_S+ABZ*I_ERI_G4x_Pz_S_S;
  Double I_ERI_G3xy_D2z_S_S = I_ERI_H3xyz_Pz_S_S+ABZ*I_ERI_G3xy_Pz_S_S;
  Double I_ERI_G3xz_D2z_S_S = I_ERI_H3x2z_Pz_S_S+ABZ*I_ERI_G3xz_Pz_S_S;
  Double I_ERI_G2x2y_D2z_S_S = I_ERI_H2x2yz_Pz_S_S+ABZ*I_ERI_G2x2y_Pz_S_S;
  Double I_ERI_G2xyz_D2z_S_S = I_ERI_H2xy2z_Pz_S_S+ABZ*I_ERI_G2xyz_Pz_S_S;
  Double I_ERI_G2x2z_D2z_S_S = I_ERI_H2x3z_Pz_S_S+ABZ*I_ERI_G2x2z_Pz_S_S;
  Double I_ERI_Gx3y_D2z_S_S = I_ERI_Hx3yz_Pz_S_S+ABZ*I_ERI_Gx3y_Pz_S_S;
  Double I_ERI_Gx2yz_D2z_S_S = I_ERI_Hx2y2z_Pz_S_S+ABZ*I_ERI_Gx2yz_Pz_S_S;
  Double I_ERI_Gxy2z_D2z_S_S = I_ERI_Hxy3z_Pz_S_S+ABZ*I_ERI_Gxy2z_Pz_S_S;
  Double I_ERI_Gx3z_D2z_S_S = I_ERI_Hx4z_Pz_S_S+ABZ*I_ERI_Gx3z_Pz_S_S;
  Double I_ERI_G4y_D2z_S_S = I_ERI_H4yz_Pz_S_S+ABZ*I_ERI_G4y_Pz_S_S;
  Double I_ERI_G3yz_D2z_S_S = I_ERI_H3y2z_Pz_S_S+ABZ*I_ERI_G3yz_Pz_S_S;
  Double I_ERI_G2y2z_D2z_S_S = I_ERI_H2y3z_Pz_S_S+ABZ*I_ERI_G2y2z_Pz_S_S;
  Double I_ERI_Gy3z_D2z_S_S = I_ERI_Hy4z_Pz_S_S+ABZ*I_ERI_Gy3z_Pz_S_S;
  Double I_ERI_G4z_D2z_S_S = I_ERI_H5z_Pz_S_S+ABZ*I_ERI_G4z_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_I_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_S_S_S
   * RHS shell quartet name: SQ_ERI_I_S_S_S
   ************************************************************/
  Double I_ERI_I6x_Px_S_S = I_ERI_K7x_S_S_S+ABX*I_ERI_I6x_S_S_S;
  Double I_ERI_I5xy_Px_S_S = I_ERI_K6xy_S_S_S+ABX*I_ERI_I5xy_S_S_S;
  Double I_ERI_I5xz_Px_S_S = I_ERI_K6xz_S_S_S+ABX*I_ERI_I5xz_S_S_S;
  Double I_ERI_I4x2y_Px_S_S = I_ERI_K5x2y_S_S_S+ABX*I_ERI_I4x2y_S_S_S;
  Double I_ERI_I4xyz_Px_S_S = I_ERI_K5xyz_S_S_S+ABX*I_ERI_I4xyz_S_S_S;
  Double I_ERI_I4x2z_Px_S_S = I_ERI_K5x2z_S_S_S+ABX*I_ERI_I4x2z_S_S_S;
  Double I_ERI_I3x3y_Px_S_S = I_ERI_K4x3y_S_S_S+ABX*I_ERI_I3x3y_S_S_S;
  Double I_ERI_I3x2yz_Px_S_S = I_ERI_K4x2yz_S_S_S+ABX*I_ERI_I3x2yz_S_S_S;
  Double I_ERI_I3xy2z_Px_S_S = I_ERI_K4xy2z_S_S_S+ABX*I_ERI_I3xy2z_S_S_S;
  Double I_ERI_I3x3z_Px_S_S = I_ERI_K4x3z_S_S_S+ABX*I_ERI_I3x3z_S_S_S;
  Double I_ERI_I2x4y_Px_S_S = I_ERI_K3x4y_S_S_S+ABX*I_ERI_I2x4y_S_S_S;
  Double I_ERI_I2x3yz_Px_S_S = I_ERI_K3x3yz_S_S_S+ABX*I_ERI_I2x3yz_S_S_S;
  Double I_ERI_I2x2y2z_Px_S_S = I_ERI_K3x2y2z_S_S_S+ABX*I_ERI_I2x2y2z_S_S_S;
  Double I_ERI_I2xy3z_Px_S_S = I_ERI_K3xy3z_S_S_S+ABX*I_ERI_I2xy3z_S_S_S;
  Double I_ERI_I2x4z_Px_S_S = I_ERI_K3x4z_S_S_S+ABX*I_ERI_I2x4z_S_S_S;
  Double I_ERI_Ix5y_Px_S_S = I_ERI_K2x5y_S_S_S+ABX*I_ERI_Ix5y_S_S_S;
  Double I_ERI_Ix4yz_Px_S_S = I_ERI_K2x4yz_S_S_S+ABX*I_ERI_Ix4yz_S_S_S;
  Double I_ERI_Ix3y2z_Px_S_S = I_ERI_K2x3y2z_S_S_S+ABX*I_ERI_Ix3y2z_S_S_S;
  Double I_ERI_Ix2y3z_Px_S_S = I_ERI_K2x2y3z_S_S_S+ABX*I_ERI_Ix2y3z_S_S_S;
  Double I_ERI_Ixy4z_Px_S_S = I_ERI_K2xy4z_S_S_S+ABX*I_ERI_Ixy4z_S_S_S;
  Double I_ERI_Ix5z_Px_S_S = I_ERI_K2x5z_S_S_S+ABX*I_ERI_Ix5z_S_S_S;
  Double I_ERI_I6y_Px_S_S = I_ERI_Kx6y_S_S_S+ABX*I_ERI_I6y_S_S_S;
  Double I_ERI_I5yz_Px_S_S = I_ERI_Kx5yz_S_S_S+ABX*I_ERI_I5yz_S_S_S;
  Double I_ERI_I4y2z_Px_S_S = I_ERI_Kx4y2z_S_S_S+ABX*I_ERI_I4y2z_S_S_S;
  Double I_ERI_I3y3z_Px_S_S = I_ERI_Kx3y3z_S_S_S+ABX*I_ERI_I3y3z_S_S_S;
  Double I_ERI_I2y4z_Px_S_S = I_ERI_Kx2y4z_S_S_S+ABX*I_ERI_I2y4z_S_S_S;
  Double I_ERI_Iy5z_Px_S_S = I_ERI_Kxy5z_S_S_S+ABX*I_ERI_Iy5z_S_S_S;
  Double I_ERI_I6z_Px_S_S = I_ERI_Kx6z_S_S_S+ABX*I_ERI_I6z_S_S_S;
  Double I_ERI_I5xy_Py_S_S = I_ERI_K5x2y_S_S_S+ABY*I_ERI_I5xy_S_S_S;
  Double I_ERI_I5xz_Py_S_S = I_ERI_K5xyz_S_S_S+ABY*I_ERI_I5xz_S_S_S;
  Double I_ERI_I4x2y_Py_S_S = I_ERI_K4x3y_S_S_S+ABY*I_ERI_I4x2y_S_S_S;
  Double I_ERI_I4xyz_Py_S_S = I_ERI_K4x2yz_S_S_S+ABY*I_ERI_I4xyz_S_S_S;
  Double I_ERI_I4x2z_Py_S_S = I_ERI_K4xy2z_S_S_S+ABY*I_ERI_I4x2z_S_S_S;
  Double I_ERI_I3x3y_Py_S_S = I_ERI_K3x4y_S_S_S+ABY*I_ERI_I3x3y_S_S_S;
  Double I_ERI_I3x2yz_Py_S_S = I_ERI_K3x3yz_S_S_S+ABY*I_ERI_I3x2yz_S_S_S;
  Double I_ERI_I3xy2z_Py_S_S = I_ERI_K3x2y2z_S_S_S+ABY*I_ERI_I3xy2z_S_S_S;
  Double I_ERI_I3x3z_Py_S_S = I_ERI_K3xy3z_S_S_S+ABY*I_ERI_I3x3z_S_S_S;
  Double I_ERI_I2x4y_Py_S_S = I_ERI_K2x5y_S_S_S+ABY*I_ERI_I2x4y_S_S_S;
  Double I_ERI_I2x3yz_Py_S_S = I_ERI_K2x4yz_S_S_S+ABY*I_ERI_I2x3yz_S_S_S;
  Double I_ERI_I2x2y2z_Py_S_S = I_ERI_K2x3y2z_S_S_S+ABY*I_ERI_I2x2y2z_S_S_S;
  Double I_ERI_I2xy3z_Py_S_S = I_ERI_K2x2y3z_S_S_S+ABY*I_ERI_I2xy3z_S_S_S;
  Double I_ERI_I2x4z_Py_S_S = I_ERI_K2xy4z_S_S_S+ABY*I_ERI_I2x4z_S_S_S;
  Double I_ERI_Ix5y_Py_S_S = I_ERI_Kx6y_S_S_S+ABY*I_ERI_Ix5y_S_S_S;
  Double I_ERI_Ix4yz_Py_S_S = I_ERI_Kx5yz_S_S_S+ABY*I_ERI_Ix4yz_S_S_S;
  Double I_ERI_Ix3y2z_Py_S_S = I_ERI_Kx4y2z_S_S_S+ABY*I_ERI_Ix3y2z_S_S_S;
  Double I_ERI_Ix2y3z_Py_S_S = I_ERI_Kx3y3z_S_S_S+ABY*I_ERI_Ix2y3z_S_S_S;
  Double I_ERI_Ixy4z_Py_S_S = I_ERI_Kx2y4z_S_S_S+ABY*I_ERI_Ixy4z_S_S_S;
  Double I_ERI_Ix5z_Py_S_S = I_ERI_Kxy5z_S_S_S+ABY*I_ERI_Ix5z_S_S_S;
  Double I_ERI_I6y_Py_S_S = I_ERI_K7y_S_S_S+ABY*I_ERI_I6y_S_S_S;
  Double I_ERI_I5yz_Py_S_S = I_ERI_K6yz_S_S_S+ABY*I_ERI_I5yz_S_S_S;
  Double I_ERI_I4y2z_Py_S_S = I_ERI_K5y2z_S_S_S+ABY*I_ERI_I4y2z_S_S_S;
  Double I_ERI_I3y3z_Py_S_S = I_ERI_K4y3z_S_S_S+ABY*I_ERI_I3y3z_S_S_S;
  Double I_ERI_I2y4z_Py_S_S = I_ERI_K3y4z_S_S_S+ABY*I_ERI_I2y4z_S_S_S;
  Double I_ERI_Iy5z_Py_S_S = I_ERI_K2y5z_S_S_S+ABY*I_ERI_Iy5z_S_S_S;
  Double I_ERI_I6z_Py_S_S = I_ERI_Ky6z_S_S_S+ABY*I_ERI_I6z_S_S_S;
  Double I_ERI_I5xy_Pz_S_S = I_ERI_K5xyz_S_S_S+ABZ*I_ERI_I5xy_S_S_S;
  Double I_ERI_I5xz_Pz_S_S = I_ERI_K5x2z_S_S_S+ABZ*I_ERI_I5xz_S_S_S;
  Double I_ERI_I4x2y_Pz_S_S = I_ERI_K4x2yz_S_S_S+ABZ*I_ERI_I4x2y_S_S_S;
  Double I_ERI_I4xyz_Pz_S_S = I_ERI_K4xy2z_S_S_S+ABZ*I_ERI_I4xyz_S_S_S;
  Double I_ERI_I4x2z_Pz_S_S = I_ERI_K4x3z_S_S_S+ABZ*I_ERI_I4x2z_S_S_S;
  Double I_ERI_I3x3y_Pz_S_S = I_ERI_K3x3yz_S_S_S+ABZ*I_ERI_I3x3y_S_S_S;
  Double I_ERI_I3x2yz_Pz_S_S = I_ERI_K3x2y2z_S_S_S+ABZ*I_ERI_I3x2yz_S_S_S;
  Double I_ERI_I3xy2z_Pz_S_S = I_ERI_K3xy3z_S_S_S+ABZ*I_ERI_I3xy2z_S_S_S;
  Double I_ERI_I3x3z_Pz_S_S = I_ERI_K3x4z_S_S_S+ABZ*I_ERI_I3x3z_S_S_S;
  Double I_ERI_I2x4y_Pz_S_S = I_ERI_K2x4yz_S_S_S+ABZ*I_ERI_I2x4y_S_S_S;
  Double I_ERI_I2x3yz_Pz_S_S = I_ERI_K2x3y2z_S_S_S+ABZ*I_ERI_I2x3yz_S_S_S;
  Double I_ERI_I2x2y2z_Pz_S_S = I_ERI_K2x2y3z_S_S_S+ABZ*I_ERI_I2x2y2z_S_S_S;
  Double I_ERI_I2xy3z_Pz_S_S = I_ERI_K2xy4z_S_S_S+ABZ*I_ERI_I2xy3z_S_S_S;
  Double I_ERI_I2x4z_Pz_S_S = I_ERI_K2x5z_S_S_S+ABZ*I_ERI_I2x4z_S_S_S;
  Double I_ERI_Ix5y_Pz_S_S = I_ERI_Kx5yz_S_S_S+ABZ*I_ERI_Ix5y_S_S_S;
  Double I_ERI_Ix4yz_Pz_S_S = I_ERI_Kx4y2z_S_S_S+ABZ*I_ERI_Ix4yz_S_S_S;
  Double I_ERI_Ix3y2z_Pz_S_S = I_ERI_Kx3y3z_S_S_S+ABZ*I_ERI_Ix3y2z_S_S_S;
  Double I_ERI_Ix2y3z_Pz_S_S = I_ERI_Kx2y4z_S_S_S+ABZ*I_ERI_Ix2y3z_S_S_S;
  Double I_ERI_Ixy4z_Pz_S_S = I_ERI_Kxy5z_S_S_S+ABZ*I_ERI_Ixy4z_S_S_S;
  Double I_ERI_Ix5z_Pz_S_S = I_ERI_Kx6z_S_S_S+ABZ*I_ERI_Ix5z_S_S_S;
  Double I_ERI_I5yz_Pz_S_S = I_ERI_K5y2z_S_S_S+ABZ*I_ERI_I5yz_S_S_S;
  Double I_ERI_I4y2z_Pz_S_S = I_ERI_K4y3z_S_S_S+ABZ*I_ERI_I4y2z_S_S_S;
  Double I_ERI_I3y3z_Pz_S_S = I_ERI_K3y4z_S_S_S+ABZ*I_ERI_I3y3z_S_S_S;
  Double I_ERI_I2y4z_Pz_S_S = I_ERI_K2y5z_S_S_S+ABZ*I_ERI_I2y4z_S_S_S;
  Double I_ERI_Iy5z_Pz_S_S = I_ERI_Ky6z_S_S_S+ABZ*I_ERI_Iy5z_S_S_S;
  Double I_ERI_I6z_Pz_S_S = I_ERI_K7z_S_S_S+ABZ*I_ERI_I6z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_D_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_P_S_S
   * RHS shell quartet name: SQ_ERI_H_P_S_S
   ************************************************************/
  Double I_ERI_H5x_D2x_S_S = I_ERI_I6x_Px_S_S+ABX*I_ERI_H5x_Px_S_S;
  Double I_ERI_H4xy_D2x_S_S = I_ERI_I5xy_Px_S_S+ABX*I_ERI_H4xy_Px_S_S;
  Double I_ERI_H4xz_D2x_S_S = I_ERI_I5xz_Px_S_S+ABX*I_ERI_H4xz_Px_S_S;
  Double I_ERI_H3x2y_D2x_S_S = I_ERI_I4x2y_Px_S_S+ABX*I_ERI_H3x2y_Px_S_S;
  Double I_ERI_H3xyz_D2x_S_S = I_ERI_I4xyz_Px_S_S+ABX*I_ERI_H3xyz_Px_S_S;
  Double I_ERI_H3x2z_D2x_S_S = I_ERI_I4x2z_Px_S_S+ABX*I_ERI_H3x2z_Px_S_S;
  Double I_ERI_H2x3y_D2x_S_S = I_ERI_I3x3y_Px_S_S+ABX*I_ERI_H2x3y_Px_S_S;
  Double I_ERI_H2x2yz_D2x_S_S = I_ERI_I3x2yz_Px_S_S+ABX*I_ERI_H2x2yz_Px_S_S;
  Double I_ERI_H2xy2z_D2x_S_S = I_ERI_I3xy2z_Px_S_S+ABX*I_ERI_H2xy2z_Px_S_S;
  Double I_ERI_H2x3z_D2x_S_S = I_ERI_I3x3z_Px_S_S+ABX*I_ERI_H2x3z_Px_S_S;
  Double I_ERI_Hx4y_D2x_S_S = I_ERI_I2x4y_Px_S_S+ABX*I_ERI_Hx4y_Px_S_S;
  Double I_ERI_Hx3yz_D2x_S_S = I_ERI_I2x3yz_Px_S_S+ABX*I_ERI_Hx3yz_Px_S_S;
  Double I_ERI_Hx2y2z_D2x_S_S = I_ERI_I2x2y2z_Px_S_S+ABX*I_ERI_Hx2y2z_Px_S_S;
  Double I_ERI_Hxy3z_D2x_S_S = I_ERI_I2xy3z_Px_S_S+ABX*I_ERI_Hxy3z_Px_S_S;
  Double I_ERI_Hx4z_D2x_S_S = I_ERI_I2x4z_Px_S_S+ABX*I_ERI_Hx4z_Px_S_S;
  Double I_ERI_H5y_D2x_S_S = I_ERI_Ix5y_Px_S_S+ABX*I_ERI_H5y_Px_S_S;
  Double I_ERI_H4yz_D2x_S_S = I_ERI_Ix4yz_Px_S_S+ABX*I_ERI_H4yz_Px_S_S;
  Double I_ERI_H3y2z_D2x_S_S = I_ERI_Ix3y2z_Px_S_S+ABX*I_ERI_H3y2z_Px_S_S;
  Double I_ERI_H2y3z_D2x_S_S = I_ERI_Ix2y3z_Px_S_S+ABX*I_ERI_H2y3z_Px_S_S;
  Double I_ERI_Hy4z_D2x_S_S = I_ERI_Ixy4z_Px_S_S+ABX*I_ERI_Hy4z_Px_S_S;
  Double I_ERI_H5z_D2x_S_S = I_ERI_Ix5z_Px_S_S+ABX*I_ERI_H5z_Px_S_S;
  Double I_ERI_H5x_D2y_S_S = I_ERI_I5xy_Py_S_S+ABY*I_ERI_H5x_Py_S_S;
  Double I_ERI_H4xy_D2y_S_S = I_ERI_I4x2y_Py_S_S+ABY*I_ERI_H4xy_Py_S_S;
  Double I_ERI_H4xz_D2y_S_S = I_ERI_I4xyz_Py_S_S+ABY*I_ERI_H4xz_Py_S_S;
  Double I_ERI_H3x2y_D2y_S_S = I_ERI_I3x3y_Py_S_S+ABY*I_ERI_H3x2y_Py_S_S;
  Double I_ERI_H3xyz_D2y_S_S = I_ERI_I3x2yz_Py_S_S+ABY*I_ERI_H3xyz_Py_S_S;
  Double I_ERI_H3x2z_D2y_S_S = I_ERI_I3xy2z_Py_S_S+ABY*I_ERI_H3x2z_Py_S_S;
  Double I_ERI_H2x3y_D2y_S_S = I_ERI_I2x4y_Py_S_S+ABY*I_ERI_H2x3y_Py_S_S;
  Double I_ERI_H2x2yz_D2y_S_S = I_ERI_I2x3yz_Py_S_S+ABY*I_ERI_H2x2yz_Py_S_S;
  Double I_ERI_H2xy2z_D2y_S_S = I_ERI_I2x2y2z_Py_S_S+ABY*I_ERI_H2xy2z_Py_S_S;
  Double I_ERI_H2x3z_D2y_S_S = I_ERI_I2xy3z_Py_S_S+ABY*I_ERI_H2x3z_Py_S_S;
  Double I_ERI_Hx4y_D2y_S_S = I_ERI_Ix5y_Py_S_S+ABY*I_ERI_Hx4y_Py_S_S;
  Double I_ERI_Hx3yz_D2y_S_S = I_ERI_Ix4yz_Py_S_S+ABY*I_ERI_Hx3yz_Py_S_S;
  Double I_ERI_Hx2y2z_D2y_S_S = I_ERI_Ix3y2z_Py_S_S+ABY*I_ERI_Hx2y2z_Py_S_S;
  Double I_ERI_Hxy3z_D2y_S_S = I_ERI_Ix2y3z_Py_S_S+ABY*I_ERI_Hxy3z_Py_S_S;
  Double I_ERI_Hx4z_D2y_S_S = I_ERI_Ixy4z_Py_S_S+ABY*I_ERI_Hx4z_Py_S_S;
  Double I_ERI_H5y_D2y_S_S = I_ERI_I6y_Py_S_S+ABY*I_ERI_H5y_Py_S_S;
  Double I_ERI_H4yz_D2y_S_S = I_ERI_I5yz_Py_S_S+ABY*I_ERI_H4yz_Py_S_S;
  Double I_ERI_H3y2z_D2y_S_S = I_ERI_I4y2z_Py_S_S+ABY*I_ERI_H3y2z_Py_S_S;
  Double I_ERI_H2y3z_D2y_S_S = I_ERI_I3y3z_Py_S_S+ABY*I_ERI_H2y3z_Py_S_S;
  Double I_ERI_Hy4z_D2y_S_S = I_ERI_I2y4z_Py_S_S+ABY*I_ERI_Hy4z_Py_S_S;
  Double I_ERI_H5z_D2y_S_S = I_ERI_Iy5z_Py_S_S+ABY*I_ERI_H5z_Py_S_S;
  Double I_ERI_H5x_D2z_S_S = I_ERI_I5xz_Pz_S_S+ABZ*I_ERI_H5x_Pz_S_S;
  Double I_ERI_H4xy_D2z_S_S = I_ERI_I4xyz_Pz_S_S+ABZ*I_ERI_H4xy_Pz_S_S;
  Double I_ERI_H4xz_D2z_S_S = I_ERI_I4x2z_Pz_S_S+ABZ*I_ERI_H4xz_Pz_S_S;
  Double I_ERI_H3x2y_D2z_S_S = I_ERI_I3x2yz_Pz_S_S+ABZ*I_ERI_H3x2y_Pz_S_S;
  Double I_ERI_H3xyz_D2z_S_S = I_ERI_I3xy2z_Pz_S_S+ABZ*I_ERI_H3xyz_Pz_S_S;
  Double I_ERI_H3x2z_D2z_S_S = I_ERI_I3x3z_Pz_S_S+ABZ*I_ERI_H3x2z_Pz_S_S;
  Double I_ERI_H2x3y_D2z_S_S = I_ERI_I2x3yz_Pz_S_S+ABZ*I_ERI_H2x3y_Pz_S_S;
  Double I_ERI_H2x2yz_D2z_S_S = I_ERI_I2x2y2z_Pz_S_S+ABZ*I_ERI_H2x2yz_Pz_S_S;
  Double I_ERI_H2xy2z_D2z_S_S = I_ERI_I2xy3z_Pz_S_S+ABZ*I_ERI_H2xy2z_Pz_S_S;
  Double I_ERI_H2x3z_D2z_S_S = I_ERI_I2x4z_Pz_S_S+ABZ*I_ERI_H2x3z_Pz_S_S;
  Double I_ERI_Hx4y_D2z_S_S = I_ERI_Ix4yz_Pz_S_S+ABZ*I_ERI_Hx4y_Pz_S_S;
  Double I_ERI_Hx3yz_D2z_S_S = I_ERI_Ix3y2z_Pz_S_S+ABZ*I_ERI_Hx3yz_Pz_S_S;
  Double I_ERI_Hx2y2z_D2z_S_S = I_ERI_Ix2y3z_Pz_S_S+ABZ*I_ERI_Hx2y2z_Pz_S_S;
  Double I_ERI_Hxy3z_D2z_S_S = I_ERI_Ixy4z_Pz_S_S+ABZ*I_ERI_Hxy3z_Pz_S_S;
  Double I_ERI_Hx4z_D2z_S_S = I_ERI_Ix5z_Pz_S_S+ABZ*I_ERI_Hx4z_Pz_S_S;
  Double I_ERI_H5y_D2z_S_S = I_ERI_I5yz_Pz_S_S+ABZ*I_ERI_H5y_Pz_S_S;
  Double I_ERI_H4yz_D2z_S_S = I_ERI_I4y2z_Pz_S_S+ABZ*I_ERI_H4yz_Pz_S_S;
  Double I_ERI_H3y2z_D2z_S_S = I_ERI_I3y3z_Pz_S_S+ABZ*I_ERI_H3y2z_Pz_S_S;
  Double I_ERI_H2y3z_D2z_S_S = I_ERI_I2y4z_Pz_S_S+ABZ*I_ERI_H2y3z_Pz_S_S;
  Double I_ERI_Hy4z_D2z_S_S = I_ERI_Iy5z_Pz_S_S+ABZ*I_ERI_Hy4z_Pz_S_S;
  Double I_ERI_H5z_D2z_S_S = I_ERI_I6z_Pz_S_S+ABZ*I_ERI_H5z_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_F_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S
   * RHS shell quartet name: SQ_ERI_G_D_S_S
   ************************************************************/
  Double I_ERI_G4x_F3x_S_S = I_ERI_H5x_D2x_S_S+ABX*I_ERI_G4x_D2x_S_S;
  Double I_ERI_G3xy_F3x_S_S = I_ERI_H4xy_D2x_S_S+ABX*I_ERI_G3xy_D2x_S_S;
  Double I_ERI_G3xz_F3x_S_S = I_ERI_H4xz_D2x_S_S+ABX*I_ERI_G3xz_D2x_S_S;
  Double I_ERI_G2x2y_F3x_S_S = I_ERI_H3x2y_D2x_S_S+ABX*I_ERI_G2x2y_D2x_S_S;
  Double I_ERI_G2xyz_F3x_S_S = I_ERI_H3xyz_D2x_S_S+ABX*I_ERI_G2xyz_D2x_S_S;
  Double I_ERI_G2x2z_F3x_S_S = I_ERI_H3x2z_D2x_S_S+ABX*I_ERI_G2x2z_D2x_S_S;
  Double I_ERI_Gx3y_F3x_S_S = I_ERI_H2x3y_D2x_S_S+ABX*I_ERI_Gx3y_D2x_S_S;
  Double I_ERI_Gx2yz_F3x_S_S = I_ERI_H2x2yz_D2x_S_S+ABX*I_ERI_Gx2yz_D2x_S_S;
  Double I_ERI_Gxy2z_F3x_S_S = I_ERI_H2xy2z_D2x_S_S+ABX*I_ERI_Gxy2z_D2x_S_S;
  Double I_ERI_Gx3z_F3x_S_S = I_ERI_H2x3z_D2x_S_S+ABX*I_ERI_Gx3z_D2x_S_S;
  Double I_ERI_G4y_F3x_S_S = I_ERI_Hx4y_D2x_S_S+ABX*I_ERI_G4y_D2x_S_S;
  Double I_ERI_G3yz_F3x_S_S = I_ERI_Hx3yz_D2x_S_S+ABX*I_ERI_G3yz_D2x_S_S;
  Double I_ERI_G2y2z_F3x_S_S = I_ERI_Hx2y2z_D2x_S_S+ABX*I_ERI_G2y2z_D2x_S_S;
  Double I_ERI_Gy3z_F3x_S_S = I_ERI_Hxy3z_D2x_S_S+ABX*I_ERI_Gy3z_D2x_S_S;
  Double I_ERI_G4z_F3x_S_S = I_ERI_Hx4z_D2x_S_S+ABX*I_ERI_G4z_D2x_S_S;
  Double I_ERI_G4x_F2xy_S_S = I_ERI_H4xy_D2x_S_S+ABY*I_ERI_G4x_D2x_S_S;
  Double I_ERI_G3xy_F2xy_S_S = I_ERI_H3x2y_D2x_S_S+ABY*I_ERI_G3xy_D2x_S_S;
  Double I_ERI_G3xz_F2xy_S_S = I_ERI_H3xyz_D2x_S_S+ABY*I_ERI_G3xz_D2x_S_S;
  Double I_ERI_G2x2y_F2xy_S_S = I_ERI_H2x3y_D2x_S_S+ABY*I_ERI_G2x2y_D2x_S_S;
  Double I_ERI_G2xyz_F2xy_S_S = I_ERI_H2x2yz_D2x_S_S+ABY*I_ERI_G2xyz_D2x_S_S;
  Double I_ERI_G2x2z_F2xy_S_S = I_ERI_H2xy2z_D2x_S_S+ABY*I_ERI_G2x2z_D2x_S_S;
  Double I_ERI_Gx3y_F2xy_S_S = I_ERI_Hx4y_D2x_S_S+ABY*I_ERI_Gx3y_D2x_S_S;
  Double I_ERI_Gx2yz_F2xy_S_S = I_ERI_Hx3yz_D2x_S_S+ABY*I_ERI_Gx2yz_D2x_S_S;
  Double I_ERI_Gxy2z_F2xy_S_S = I_ERI_Hx2y2z_D2x_S_S+ABY*I_ERI_Gxy2z_D2x_S_S;
  Double I_ERI_Gx3z_F2xy_S_S = I_ERI_Hxy3z_D2x_S_S+ABY*I_ERI_Gx3z_D2x_S_S;
  Double I_ERI_G4y_F2xy_S_S = I_ERI_H5y_D2x_S_S+ABY*I_ERI_G4y_D2x_S_S;
  Double I_ERI_G3yz_F2xy_S_S = I_ERI_H4yz_D2x_S_S+ABY*I_ERI_G3yz_D2x_S_S;
  Double I_ERI_G2y2z_F2xy_S_S = I_ERI_H3y2z_D2x_S_S+ABY*I_ERI_G2y2z_D2x_S_S;
  Double I_ERI_Gy3z_F2xy_S_S = I_ERI_H2y3z_D2x_S_S+ABY*I_ERI_Gy3z_D2x_S_S;
  Double I_ERI_G4z_F2xy_S_S = I_ERI_Hy4z_D2x_S_S+ABY*I_ERI_G4z_D2x_S_S;
  Double I_ERI_G4x_F2xz_S_S = I_ERI_H4xz_D2x_S_S+ABZ*I_ERI_G4x_D2x_S_S;
  Double I_ERI_G3xy_F2xz_S_S = I_ERI_H3xyz_D2x_S_S+ABZ*I_ERI_G3xy_D2x_S_S;
  Double I_ERI_G3xz_F2xz_S_S = I_ERI_H3x2z_D2x_S_S+ABZ*I_ERI_G3xz_D2x_S_S;
  Double I_ERI_G2x2y_F2xz_S_S = I_ERI_H2x2yz_D2x_S_S+ABZ*I_ERI_G2x2y_D2x_S_S;
  Double I_ERI_G2xyz_F2xz_S_S = I_ERI_H2xy2z_D2x_S_S+ABZ*I_ERI_G2xyz_D2x_S_S;
  Double I_ERI_G2x2z_F2xz_S_S = I_ERI_H2x3z_D2x_S_S+ABZ*I_ERI_G2x2z_D2x_S_S;
  Double I_ERI_Gx3y_F2xz_S_S = I_ERI_Hx3yz_D2x_S_S+ABZ*I_ERI_Gx3y_D2x_S_S;
  Double I_ERI_Gx2yz_F2xz_S_S = I_ERI_Hx2y2z_D2x_S_S+ABZ*I_ERI_Gx2yz_D2x_S_S;
  Double I_ERI_Gxy2z_F2xz_S_S = I_ERI_Hxy3z_D2x_S_S+ABZ*I_ERI_Gxy2z_D2x_S_S;
  Double I_ERI_Gx3z_F2xz_S_S = I_ERI_Hx4z_D2x_S_S+ABZ*I_ERI_Gx3z_D2x_S_S;
  Double I_ERI_G4y_F2xz_S_S = I_ERI_H4yz_D2x_S_S+ABZ*I_ERI_G4y_D2x_S_S;
  Double I_ERI_G3yz_F2xz_S_S = I_ERI_H3y2z_D2x_S_S+ABZ*I_ERI_G3yz_D2x_S_S;
  Double I_ERI_G2y2z_F2xz_S_S = I_ERI_H2y3z_D2x_S_S+ABZ*I_ERI_G2y2z_D2x_S_S;
  Double I_ERI_Gy3z_F2xz_S_S = I_ERI_Hy4z_D2x_S_S+ABZ*I_ERI_Gy3z_D2x_S_S;
  Double I_ERI_G4z_F2xz_S_S = I_ERI_H5z_D2x_S_S+ABZ*I_ERI_G4z_D2x_S_S;
  Double I_ERI_G4x_Fx2y_S_S = I_ERI_H5x_D2y_S_S+ABX*I_ERI_G4x_D2y_S_S;
  Double I_ERI_G3xy_Fx2y_S_S = I_ERI_H4xy_D2y_S_S+ABX*I_ERI_G3xy_D2y_S_S;
  Double I_ERI_G3xz_Fx2y_S_S = I_ERI_H4xz_D2y_S_S+ABX*I_ERI_G3xz_D2y_S_S;
  Double I_ERI_G2x2y_Fx2y_S_S = I_ERI_H3x2y_D2y_S_S+ABX*I_ERI_G2x2y_D2y_S_S;
  Double I_ERI_G2xyz_Fx2y_S_S = I_ERI_H3xyz_D2y_S_S+ABX*I_ERI_G2xyz_D2y_S_S;
  Double I_ERI_G2x2z_Fx2y_S_S = I_ERI_H3x2z_D2y_S_S+ABX*I_ERI_G2x2z_D2y_S_S;
  Double I_ERI_Gx3y_Fx2y_S_S = I_ERI_H2x3y_D2y_S_S+ABX*I_ERI_Gx3y_D2y_S_S;
  Double I_ERI_Gx2yz_Fx2y_S_S = I_ERI_H2x2yz_D2y_S_S+ABX*I_ERI_Gx2yz_D2y_S_S;
  Double I_ERI_Gxy2z_Fx2y_S_S = I_ERI_H2xy2z_D2y_S_S+ABX*I_ERI_Gxy2z_D2y_S_S;
  Double I_ERI_Gx3z_Fx2y_S_S = I_ERI_H2x3z_D2y_S_S+ABX*I_ERI_Gx3z_D2y_S_S;
  Double I_ERI_G4y_Fx2y_S_S = I_ERI_Hx4y_D2y_S_S+ABX*I_ERI_G4y_D2y_S_S;
  Double I_ERI_G3yz_Fx2y_S_S = I_ERI_Hx3yz_D2y_S_S+ABX*I_ERI_G3yz_D2y_S_S;
  Double I_ERI_G2y2z_Fx2y_S_S = I_ERI_Hx2y2z_D2y_S_S+ABX*I_ERI_G2y2z_D2y_S_S;
  Double I_ERI_Gy3z_Fx2y_S_S = I_ERI_Hxy3z_D2y_S_S+ABX*I_ERI_Gy3z_D2y_S_S;
  Double I_ERI_G4z_Fx2y_S_S = I_ERI_Hx4z_D2y_S_S+ABX*I_ERI_G4z_D2y_S_S;
  Double I_ERI_G4x_Fx2z_S_S = I_ERI_H5x_D2z_S_S+ABX*I_ERI_G4x_D2z_S_S;
  Double I_ERI_G3xy_Fx2z_S_S = I_ERI_H4xy_D2z_S_S+ABX*I_ERI_G3xy_D2z_S_S;
  Double I_ERI_G3xz_Fx2z_S_S = I_ERI_H4xz_D2z_S_S+ABX*I_ERI_G3xz_D2z_S_S;
  Double I_ERI_G2x2y_Fx2z_S_S = I_ERI_H3x2y_D2z_S_S+ABX*I_ERI_G2x2y_D2z_S_S;
  Double I_ERI_G2xyz_Fx2z_S_S = I_ERI_H3xyz_D2z_S_S+ABX*I_ERI_G2xyz_D2z_S_S;
  Double I_ERI_G2x2z_Fx2z_S_S = I_ERI_H3x2z_D2z_S_S+ABX*I_ERI_G2x2z_D2z_S_S;
  Double I_ERI_Gx3y_Fx2z_S_S = I_ERI_H2x3y_D2z_S_S+ABX*I_ERI_Gx3y_D2z_S_S;
  Double I_ERI_Gx2yz_Fx2z_S_S = I_ERI_H2x2yz_D2z_S_S+ABX*I_ERI_Gx2yz_D2z_S_S;
  Double I_ERI_Gxy2z_Fx2z_S_S = I_ERI_H2xy2z_D2z_S_S+ABX*I_ERI_Gxy2z_D2z_S_S;
  Double I_ERI_Gx3z_Fx2z_S_S = I_ERI_H2x3z_D2z_S_S+ABX*I_ERI_Gx3z_D2z_S_S;
  Double I_ERI_G4y_Fx2z_S_S = I_ERI_Hx4y_D2z_S_S+ABX*I_ERI_G4y_D2z_S_S;
  Double I_ERI_G3yz_Fx2z_S_S = I_ERI_Hx3yz_D2z_S_S+ABX*I_ERI_G3yz_D2z_S_S;
  Double I_ERI_G2y2z_Fx2z_S_S = I_ERI_Hx2y2z_D2z_S_S+ABX*I_ERI_G2y2z_D2z_S_S;
  Double I_ERI_Gy3z_Fx2z_S_S = I_ERI_Hxy3z_D2z_S_S+ABX*I_ERI_Gy3z_D2z_S_S;
  Double I_ERI_G4z_Fx2z_S_S = I_ERI_Hx4z_D2z_S_S+ABX*I_ERI_G4z_D2z_S_S;
  Double I_ERI_G4x_F3y_S_S = I_ERI_H4xy_D2y_S_S+ABY*I_ERI_G4x_D2y_S_S;
  Double I_ERI_G3xy_F3y_S_S = I_ERI_H3x2y_D2y_S_S+ABY*I_ERI_G3xy_D2y_S_S;
  Double I_ERI_G3xz_F3y_S_S = I_ERI_H3xyz_D2y_S_S+ABY*I_ERI_G3xz_D2y_S_S;
  Double I_ERI_G2x2y_F3y_S_S = I_ERI_H2x3y_D2y_S_S+ABY*I_ERI_G2x2y_D2y_S_S;
  Double I_ERI_G2xyz_F3y_S_S = I_ERI_H2x2yz_D2y_S_S+ABY*I_ERI_G2xyz_D2y_S_S;
  Double I_ERI_G2x2z_F3y_S_S = I_ERI_H2xy2z_D2y_S_S+ABY*I_ERI_G2x2z_D2y_S_S;
  Double I_ERI_Gx3y_F3y_S_S = I_ERI_Hx4y_D2y_S_S+ABY*I_ERI_Gx3y_D2y_S_S;
  Double I_ERI_Gx2yz_F3y_S_S = I_ERI_Hx3yz_D2y_S_S+ABY*I_ERI_Gx2yz_D2y_S_S;
  Double I_ERI_Gxy2z_F3y_S_S = I_ERI_Hx2y2z_D2y_S_S+ABY*I_ERI_Gxy2z_D2y_S_S;
  Double I_ERI_Gx3z_F3y_S_S = I_ERI_Hxy3z_D2y_S_S+ABY*I_ERI_Gx3z_D2y_S_S;
  Double I_ERI_G4y_F3y_S_S = I_ERI_H5y_D2y_S_S+ABY*I_ERI_G4y_D2y_S_S;
  Double I_ERI_G3yz_F3y_S_S = I_ERI_H4yz_D2y_S_S+ABY*I_ERI_G3yz_D2y_S_S;
  Double I_ERI_G2y2z_F3y_S_S = I_ERI_H3y2z_D2y_S_S+ABY*I_ERI_G2y2z_D2y_S_S;
  Double I_ERI_Gy3z_F3y_S_S = I_ERI_H2y3z_D2y_S_S+ABY*I_ERI_Gy3z_D2y_S_S;
  Double I_ERI_G4z_F3y_S_S = I_ERI_Hy4z_D2y_S_S+ABY*I_ERI_G4z_D2y_S_S;
  Double I_ERI_G4x_F2yz_S_S = I_ERI_H4xz_D2y_S_S+ABZ*I_ERI_G4x_D2y_S_S;
  Double I_ERI_G3xy_F2yz_S_S = I_ERI_H3xyz_D2y_S_S+ABZ*I_ERI_G3xy_D2y_S_S;
  Double I_ERI_G3xz_F2yz_S_S = I_ERI_H3x2z_D2y_S_S+ABZ*I_ERI_G3xz_D2y_S_S;
  Double I_ERI_G2x2y_F2yz_S_S = I_ERI_H2x2yz_D2y_S_S+ABZ*I_ERI_G2x2y_D2y_S_S;
  Double I_ERI_G2xyz_F2yz_S_S = I_ERI_H2xy2z_D2y_S_S+ABZ*I_ERI_G2xyz_D2y_S_S;
  Double I_ERI_G2x2z_F2yz_S_S = I_ERI_H2x3z_D2y_S_S+ABZ*I_ERI_G2x2z_D2y_S_S;
  Double I_ERI_Gx3y_F2yz_S_S = I_ERI_Hx3yz_D2y_S_S+ABZ*I_ERI_Gx3y_D2y_S_S;
  Double I_ERI_Gx2yz_F2yz_S_S = I_ERI_Hx2y2z_D2y_S_S+ABZ*I_ERI_Gx2yz_D2y_S_S;
  Double I_ERI_Gxy2z_F2yz_S_S = I_ERI_Hxy3z_D2y_S_S+ABZ*I_ERI_Gxy2z_D2y_S_S;
  Double I_ERI_Gx3z_F2yz_S_S = I_ERI_Hx4z_D2y_S_S+ABZ*I_ERI_Gx3z_D2y_S_S;
  Double I_ERI_G4y_F2yz_S_S = I_ERI_H4yz_D2y_S_S+ABZ*I_ERI_G4y_D2y_S_S;
  Double I_ERI_G3yz_F2yz_S_S = I_ERI_H3y2z_D2y_S_S+ABZ*I_ERI_G3yz_D2y_S_S;
  Double I_ERI_G2y2z_F2yz_S_S = I_ERI_H2y3z_D2y_S_S+ABZ*I_ERI_G2y2z_D2y_S_S;
  Double I_ERI_Gy3z_F2yz_S_S = I_ERI_Hy4z_D2y_S_S+ABZ*I_ERI_Gy3z_D2y_S_S;
  Double I_ERI_G4z_F2yz_S_S = I_ERI_H5z_D2y_S_S+ABZ*I_ERI_G4z_D2y_S_S;
  Double I_ERI_G4x_F3z_S_S = I_ERI_H4xz_D2z_S_S+ABZ*I_ERI_G4x_D2z_S_S;
  Double I_ERI_G3xy_F3z_S_S = I_ERI_H3xyz_D2z_S_S+ABZ*I_ERI_G3xy_D2z_S_S;
  Double I_ERI_G3xz_F3z_S_S = I_ERI_H3x2z_D2z_S_S+ABZ*I_ERI_G3xz_D2z_S_S;
  Double I_ERI_G2x2y_F3z_S_S = I_ERI_H2x2yz_D2z_S_S+ABZ*I_ERI_G2x2y_D2z_S_S;
  Double I_ERI_G2xyz_F3z_S_S = I_ERI_H2xy2z_D2z_S_S+ABZ*I_ERI_G2xyz_D2z_S_S;
  Double I_ERI_G2x2z_F3z_S_S = I_ERI_H2x3z_D2z_S_S+ABZ*I_ERI_G2x2z_D2z_S_S;
  Double I_ERI_Gx3y_F3z_S_S = I_ERI_Hx3yz_D2z_S_S+ABZ*I_ERI_Gx3y_D2z_S_S;
  Double I_ERI_Gx2yz_F3z_S_S = I_ERI_Hx2y2z_D2z_S_S+ABZ*I_ERI_Gx2yz_D2z_S_S;
  Double I_ERI_Gxy2z_F3z_S_S = I_ERI_Hxy3z_D2z_S_S+ABZ*I_ERI_Gxy2z_D2z_S_S;
  Double I_ERI_Gx3z_F3z_S_S = I_ERI_Hx4z_D2z_S_S+ABZ*I_ERI_Gx3z_D2z_S_S;
  Double I_ERI_G4y_F3z_S_S = I_ERI_H4yz_D2z_S_S+ABZ*I_ERI_G4y_D2z_S_S;
  Double I_ERI_G3yz_F3z_S_S = I_ERI_H3y2z_D2z_S_S+ABZ*I_ERI_G3yz_D2z_S_S;
  Double I_ERI_G2y2z_F3z_S_S = I_ERI_H2y3z_D2z_S_S+ABZ*I_ERI_G2y2z_D2z_S_S;
  Double I_ERI_Gy3z_F3z_S_S = I_ERI_Hy4z_D2z_S_S+ABZ*I_ERI_Gy3z_D2z_S_S;
  Double I_ERI_G4z_F3z_S_S = I_ERI_H5z_D2z_S_S+ABZ*I_ERI_G4z_D2z_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_K_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 27 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_L_S_S_S
   * RHS shell quartet name: SQ_ERI_K_S_S_S
   ************************************************************/
  Double I_ERI_K7x_Px_S_S = I_ERI_L8x_S_S_S+ABX*I_ERI_K7x_S_S_S;
  Double I_ERI_K6xy_Px_S_S = I_ERI_L7xy_S_S_S+ABX*I_ERI_K6xy_S_S_S;
  Double I_ERI_K6xz_Px_S_S = I_ERI_L7xz_S_S_S+ABX*I_ERI_K6xz_S_S_S;
  Double I_ERI_K5x2y_Px_S_S = I_ERI_L6x2y_S_S_S+ABX*I_ERI_K5x2y_S_S_S;
  Double I_ERI_K5xyz_Px_S_S = I_ERI_L6xyz_S_S_S+ABX*I_ERI_K5xyz_S_S_S;
  Double I_ERI_K5x2z_Px_S_S = I_ERI_L6x2z_S_S_S+ABX*I_ERI_K5x2z_S_S_S;
  Double I_ERI_K4x3y_Px_S_S = I_ERI_L5x3y_S_S_S+ABX*I_ERI_K4x3y_S_S_S;
  Double I_ERI_K4x2yz_Px_S_S = I_ERI_L5x2yz_S_S_S+ABX*I_ERI_K4x2yz_S_S_S;
  Double I_ERI_K4xy2z_Px_S_S = I_ERI_L5xy2z_S_S_S+ABX*I_ERI_K4xy2z_S_S_S;
  Double I_ERI_K4x3z_Px_S_S = I_ERI_L5x3z_S_S_S+ABX*I_ERI_K4x3z_S_S_S;
  Double I_ERI_K3x4y_Px_S_S = I_ERI_L4x4y_S_S_S+ABX*I_ERI_K3x4y_S_S_S;
  Double I_ERI_K3x3yz_Px_S_S = I_ERI_L4x3yz_S_S_S+ABX*I_ERI_K3x3yz_S_S_S;
  Double I_ERI_K3x2y2z_Px_S_S = I_ERI_L4x2y2z_S_S_S+ABX*I_ERI_K3x2y2z_S_S_S;
  Double I_ERI_K3xy3z_Px_S_S = I_ERI_L4xy3z_S_S_S+ABX*I_ERI_K3xy3z_S_S_S;
  Double I_ERI_K3x4z_Px_S_S = I_ERI_L4x4z_S_S_S+ABX*I_ERI_K3x4z_S_S_S;
  Double I_ERI_K2x5y_Px_S_S = I_ERI_L3x5y_S_S_S+ABX*I_ERI_K2x5y_S_S_S;
  Double I_ERI_K2x4yz_Px_S_S = I_ERI_L3x4yz_S_S_S+ABX*I_ERI_K2x4yz_S_S_S;
  Double I_ERI_K2x3y2z_Px_S_S = I_ERI_L3x3y2z_S_S_S+ABX*I_ERI_K2x3y2z_S_S_S;
  Double I_ERI_K2x2y3z_Px_S_S = I_ERI_L3x2y3z_S_S_S+ABX*I_ERI_K2x2y3z_S_S_S;
  Double I_ERI_K2xy4z_Px_S_S = I_ERI_L3xy4z_S_S_S+ABX*I_ERI_K2xy4z_S_S_S;
  Double I_ERI_K2x5z_Px_S_S = I_ERI_L3x5z_S_S_S+ABX*I_ERI_K2x5z_S_S_S;
  Double I_ERI_Kx6y_Px_S_S = I_ERI_L2x6y_S_S_S+ABX*I_ERI_Kx6y_S_S_S;
  Double I_ERI_Kx5yz_Px_S_S = I_ERI_L2x5yz_S_S_S+ABX*I_ERI_Kx5yz_S_S_S;
  Double I_ERI_Kx4y2z_Px_S_S = I_ERI_L2x4y2z_S_S_S+ABX*I_ERI_Kx4y2z_S_S_S;
  Double I_ERI_Kx3y3z_Px_S_S = I_ERI_L2x3y3z_S_S_S+ABX*I_ERI_Kx3y3z_S_S_S;
  Double I_ERI_Kx2y4z_Px_S_S = I_ERI_L2x2y4z_S_S_S+ABX*I_ERI_Kx2y4z_S_S_S;
  Double I_ERI_Kxy5z_Px_S_S = I_ERI_L2xy5z_S_S_S+ABX*I_ERI_Kxy5z_S_S_S;
  Double I_ERI_Kx6z_Px_S_S = I_ERI_L2x6z_S_S_S+ABX*I_ERI_Kx6z_S_S_S;
  Double I_ERI_K5x2y_Py_S_S = I_ERI_L5x3y_S_S_S+ABY*I_ERI_K5x2y_S_S_S;
  Double I_ERI_K5xyz_Py_S_S = I_ERI_L5x2yz_S_S_S+ABY*I_ERI_K5xyz_S_S_S;
  Double I_ERI_K4x3y_Py_S_S = I_ERI_L4x4y_S_S_S+ABY*I_ERI_K4x3y_S_S_S;
  Double I_ERI_K4x2yz_Py_S_S = I_ERI_L4x3yz_S_S_S+ABY*I_ERI_K4x2yz_S_S_S;
  Double I_ERI_K4xy2z_Py_S_S = I_ERI_L4x2y2z_S_S_S+ABY*I_ERI_K4xy2z_S_S_S;
  Double I_ERI_K3x4y_Py_S_S = I_ERI_L3x5y_S_S_S+ABY*I_ERI_K3x4y_S_S_S;
  Double I_ERI_K3x3yz_Py_S_S = I_ERI_L3x4yz_S_S_S+ABY*I_ERI_K3x3yz_S_S_S;
  Double I_ERI_K3x2y2z_Py_S_S = I_ERI_L3x3y2z_S_S_S+ABY*I_ERI_K3x2y2z_S_S_S;
  Double I_ERI_K3xy3z_Py_S_S = I_ERI_L3x2y3z_S_S_S+ABY*I_ERI_K3xy3z_S_S_S;
  Double I_ERI_K2x5y_Py_S_S = I_ERI_L2x6y_S_S_S+ABY*I_ERI_K2x5y_S_S_S;
  Double I_ERI_K2x4yz_Py_S_S = I_ERI_L2x5yz_S_S_S+ABY*I_ERI_K2x4yz_S_S_S;
  Double I_ERI_K2x3y2z_Py_S_S = I_ERI_L2x4y2z_S_S_S+ABY*I_ERI_K2x3y2z_S_S_S;
  Double I_ERI_K2x2y3z_Py_S_S = I_ERI_L2x3y3z_S_S_S+ABY*I_ERI_K2x2y3z_S_S_S;
  Double I_ERI_K2xy4z_Py_S_S = I_ERI_L2x2y4z_S_S_S+ABY*I_ERI_K2xy4z_S_S_S;
  Double I_ERI_Kx6y_Py_S_S = I_ERI_Lx7y_S_S_S+ABY*I_ERI_Kx6y_S_S_S;
  Double I_ERI_Kx5yz_Py_S_S = I_ERI_Lx6yz_S_S_S+ABY*I_ERI_Kx5yz_S_S_S;
  Double I_ERI_Kx4y2z_Py_S_S = I_ERI_Lx5y2z_S_S_S+ABY*I_ERI_Kx4y2z_S_S_S;
  Double I_ERI_Kx3y3z_Py_S_S = I_ERI_Lx4y3z_S_S_S+ABY*I_ERI_Kx3y3z_S_S_S;
  Double I_ERI_Kx2y4z_Py_S_S = I_ERI_Lx3y4z_S_S_S+ABY*I_ERI_Kx2y4z_S_S_S;
  Double I_ERI_Kxy5z_Py_S_S = I_ERI_Lx2y5z_S_S_S+ABY*I_ERI_Kxy5z_S_S_S;
  Double I_ERI_K7y_Py_S_S = I_ERI_L8y_S_S_S+ABY*I_ERI_K7y_S_S_S;
  Double I_ERI_K6yz_Py_S_S = I_ERI_L7yz_S_S_S+ABY*I_ERI_K6yz_S_S_S;
  Double I_ERI_K5y2z_Py_S_S = I_ERI_L6y2z_S_S_S+ABY*I_ERI_K5y2z_S_S_S;
  Double I_ERI_K4y3z_Py_S_S = I_ERI_L5y3z_S_S_S+ABY*I_ERI_K4y3z_S_S_S;
  Double I_ERI_K3y4z_Py_S_S = I_ERI_L4y4z_S_S_S+ABY*I_ERI_K3y4z_S_S_S;
  Double I_ERI_K2y5z_Py_S_S = I_ERI_L3y5z_S_S_S+ABY*I_ERI_K2y5z_S_S_S;
  Double I_ERI_Ky6z_Py_S_S = I_ERI_L2y6z_S_S_S+ABY*I_ERI_Ky6z_S_S_S;
  Double I_ERI_K5xyz_Pz_S_S = I_ERI_L5xy2z_S_S_S+ABZ*I_ERI_K5xyz_S_S_S;
  Double I_ERI_K5x2z_Pz_S_S = I_ERI_L5x3z_S_S_S+ABZ*I_ERI_K5x2z_S_S_S;
  Double I_ERI_K4x2yz_Pz_S_S = I_ERI_L4x2y2z_S_S_S+ABZ*I_ERI_K4x2yz_S_S_S;
  Double I_ERI_K4xy2z_Pz_S_S = I_ERI_L4xy3z_S_S_S+ABZ*I_ERI_K4xy2z_S_S_S;
  Double I_ERI_K4x3z_Pz_S_S = I_ERI_L4x4z_S_S_S+ABZ*I_ERI_K4x3z_S_S_S;
  Double I_ERI_K3x3yz_Pz_S_S = I_ERI_L3x3y2z_S_S_S+ABZ*I_ERI_K3x3yz_S_S_S;
  Double I_ERI_K3x2y2z_Pz_S_S = I_ERI_L3x2y3z_S_S_S+ABZ*I_ERI_K3x2y2z_S_S_S;
  Double I_ERI_K3xy3z_Pz_S_S = I_ERI_L3xy4z_S_S_S+ABZ*I_ERI_K3xy3z_S_S_S;
  Double I_ERI_K3x4z_Pz_S_S = I_ERI_L3x5z_S_S_S+ABZ*I_ERI_K3x4z_S_S_S;
  Double I_ERI_K2x4yz_Pz_S_S = I_ERI_L2x4y2z_S_S_S+ABZ*I_ERI_K2x4yz_S_S_S;
  Double I_ERI_K2x3y2z_Pz_S_S = I_ERI_L2x3y3z_S_S_S+ABZ*I_ERI_K2x3y2z_S_S_S;
  Double I_ERI_K2x2y3z_Pz_S_S = I_ERI_L2x2y4z_S_S_S+ABZ*I_ERI_K2x2y3z_S_S_S;
  Double I_ERI_K2xy4z_Pz_S_S = I_ERI_L2xy5z_S_S_S+ABZ*I_ERI_K2xy4z_S_S_S;
  Double I_ERI_K2x5z_Pz_S_S = I_ERI_L2x6z_S_S_S+ABZ*I_ERI_K2x5z_S_S_S;
  Double I_ERI_Kx5yz_Pz_S_S = I_ERI_Lx5y2z_S_S_S+ABZ*I_ERI_Kx5yz_S_S_S;
  Double I_ERI_Kx4y2z_Pz_S_S = I_ERI_Lx4y3z_S_S_S+ABZ*I_ERI_Kx4y2z_S_S_S;
  Double I_ERI_Kx3y3z_Pz_S_S = I_ERI_Lx3y4z_S_S_S+ABZ*I_ERI_Kx3y3z_S_S_S;
  Double I_ERI_Kx2y4z_Pz_S_S = I_ERI_Lx2y5z_S_S_S+ABZ*I_ERI_Kx2y4z_S_S_S;
  Double I_ERI_Kxy5z_Pz_S_S = I_ERI_Lxy6z_S_S_S+ABZ*I_ERI_Kxy5z_S_S_S;
  Double I_ERI_Kx6z_Pz_S_S = I_ERI_Lx7z_S_S_S+ABZ*I_ERI_Kx6z_S_S_S;
  Double I_ERI_K5y2z_Pz_S_S = I_ERI_L5y3z_S_S_S+ABZ*I_ERI_K5y2z_S_S_S;
  Double I_ERI_K4y3z_Pz_S_S = I_ERI_L4y4z_S_S_S+ABZ*I_ERI_K4y3z_S_S_S;
  Double I_ERI_K3y4z_Pz_S_S = I_ERI_L3y5z_S_S_S+ABZ*I_ERI_K3y4z_S_S_S;
  Double I_ERI_K2y5z_Pz_S_S = I_ERI_L2y6z_S_S_S+ABZ*I_ERI_K2y5z_S_S_S;
  Double I_ERI_Ky6z_Pz_S_S = I_ERI_Ly7z_S_S_S+ABZ*I_ERI_Ky6z_S_S_S;
  Double I_ERI_K7z_Pz_S_S = I_ERI_L8z_S_S_S+ABZ*I_ERI_K7z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_I_D_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 87 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_P_S_S
   * RHS shell quartet name: SQ_ERI_I_P_S_S
   ************************************************************/
  Double I_ERI_I6x_D2x_S_S = I_ERI_K7x_Px_S_S+ABX*I_ERI_I6x_Px_S_S;
  Double I_ERI_I5xy_D2x_S_S = I_ERI_K6xy_Px_S_S+ABX*I_ERI_I5xy_Px_S_S;
  Double I_ERI_I5xz_D2x_S_S = I_ERI_K6xz_Px_S_S+ABX*I_ERI_I5xz_Px_S_S;
  Double I_ERI_I4x2y_D2x_S_S = I_ERI_K5x2y_Px_S_S+ABX*I_ERI_I4x2y_Px_S_S;
  Double I_ERI_I4xyz_D2x_S_S = I_ERI_K5xyz_Px_S_S+ABX*I_ERI_I4xyz_Px_S_S;
  Double I_ERI_I4x2z_D2x_S_S = I_ERI_K5x2z_Px_S_S+ABX*I_ERI_I4x2z_Px_S_S;
  Double I_ERI_I3x3y_D2x_S_S = I_ERI_K4x3y_Px_S_S+ABX*I_ERI_I3x3y_Px_S_S;
  Double I_ERI_I3x2yz_D2x_S_S = I_ERI_K4x2yz_Px_S_S+ABX*I_ERI_I3x2yz_Px_S_S;
  Double I_ERI_I3xy2z_D2x_S_S = I_ERI_K4xy2z_Px_S_S+ABX*I_ERI_I3xy2z_Px_S_S;
  Double I_ERI_I3x3z_D2x_S_S = I_ERI_K4x3z_Px_S_S+ABX*I_ERI_I3x3z_Px_S_S;
  Double I_ERI_I2x4y_D2x_S_S = I_ERI_K3x4y_Px_S_S+ABX*I_ERI_I2x4y_Px_S_S;
  Double I_ERI_I2x3yz_D2x_S_S = I_ERI_K3x3yz_Px_S_S+ABX*I_ERI_I2x3yz_Px_S_S;
  Double I_ERI_I2x2y2z_D2x_S_S = I_ERI_K3x2y2z_Px_S_S+ABX*I_ERI_I2x2y2z_Px_S_S;
  Double I_ERI_I2xy3z_D2x_S_S = I_ERI_K3xy3z_Px_S_S+ABX*I_ERI_I2xy3z_Px_S_S;
  Double I_ERI_I2x4z_D2x_S_S = I_ERI_K3x4z_Px_S_S+ABX*I_ERI_I2x4z_Px_S_S;
  Double I_ERI_Ix5y_D2x_S_S = I_ERI_K2x5y_Px_S_S+ABX*I_ERI_Ix5y_Px_S_S;
  Double I_ERI_Ix4yz_D2x_S_S = I_ERI_K2x4yz_Px_S_S+ABX*I_ERI_Ix4yz_Px_S_S;
  Double I_ERI_Ix3y2z_D2x_S_S = I_ERI_K2x3y2z_Px_S_S+ABX*I_ERI_Ix3y2z_Px_S_S;
  Double I_ERI_Ix2y3z_D2x_S_S = I_ERI_K2x2y3z_Px_S_S+ABX*I_ERI_Ix2y3z_Px_S_S;
  Double I_ERI_Ixy4z_D2x_S_S = I_ERI_K2xy4z_Px_S_S+ABX*I_ERI_Ixy4z_Px_S_S;
  Double I_ERI_Ix5z_D2x_S_S = I_ERI_K2x5z_Px_S_S+ABX*I_ERI_Ix5z_Px_S_S;
  Double I_ERI_I6y_D2x_S_S = I_ERI_Kx6y_Px_S_S+ABX*I_ERI_I6y_Px_S_S;
  Double I_ERI_I5yz_D2x_S_S = I_ERI_Kx5yz_Px_S_S+ABX*I_ERI_I5yz_Px_S_S;
  Double I_ERI_I4y2z_D2x_S_S = I_ERI_Kx4y2z_Px_S_S+ABX*I_ERI_I4y2z_Px_S_S;
  Double I_ERI_I3y3z_D2x_S_S = I_ERI_Kx3y3z_Px_S_S+ABX*I_ERI_I3y3z_Px_S_S;
  Double I_ERI_I2y4z_D2x_S_S = I_ERI_Kx2y4z_Px_S_S+ABX*I_ERI_I2y4z_Px_S_S;
  Double I_ERI_Iy5z_D2x_S_S = I_ERI_Kxy5z_Px_S_S+ABX*I_ERI_Iy5z_Px_S_S;
  Double I_ERI_I6z_D2x_S_S = I_ERI_Kx6z_Px_S_S+ABX*I_ERI_I6z_Px_S_S;
  Double I_ERI_I5xy_D2y_S_S = I_ERI_K5x2y_Py_S_S+ABY*I_ERI_I5xy_Py_S_S;
  Double I_ERI_I5xz_D2y_S_S = I_ERI_K5xyz_Py_S_S+ABY*I_ERI_I5xz_Py_S_S;
  Double I_ERI_I4x2y_D2y_S_S = I_ERI_K4x3y_Py_S_S+ABY*I_ERI_I4x2y_Py_S_S;
  Double I_ERI_I4xyz_D2y_S_S = I_ERI_K4x2yz_Py_S_S+ABY*I_ERI_I4xyz_Py_S_S;
  Double I_ERI_I4x2z_D2y_S_S = I_ERI_K4xy2z_Py_S_S+ABY*I_ERI_I4x2z_Py_S_S;
  Double I_ERI_I3x3y_D2y_S_S = I_ERI_K3x4y_Py_S_S+ABY*I_ERI_I3x3y_Py_S_S;
  Double I_ERI_I3x2yz_D2y_S_S = I_ERI_K3x3yz_Py_S_S+ABY*I_ERI_I3x2yz_Py_S_S;
  Double I_ERI_I3xy2z_D2y_S_S = I_ERI_K3x2y2z_Py_S_S+ABY*I_ERI_I3xy2z_Py_S_S;
  Double I_ERI_I3x3z_D2y_S_S = I_ERI_K3xy3z_Py_S_S+ABY*I_ERI_I3x3z_Py_S_S;
  Double I_ERI_I2x4y_D2y_S_S = I_ERI_K2x5y_Py_S_S+ABY*I_ERI_I2x4y_Py_S_S;
  Double I_ERI_I2x3yz_D2y_S_S = I_ERI_K2x4yz_Py_S_S+ABY*I_ERI_I2x3yz_Py_S_S;
  Double I_ERI_I2x2y2z_D2y_S_S = I_ERI_K2x3y2z_Py_S_S+ABY*I_ERI_I2x2y2z_Py_S_S;
  Double I_ERI_I2xy3z_D2y_S_S = I_ERI_K2x2y3z_Py_S_S+ABY*I_ERI_I2xy3z_Py_S_S;
  Double I_ERI_I2x4z_D2y_S_S = I_ERI_K2xy4z_Py_S_S+ABY*I_ERI_I2x4z_Py_S_S;
  Double I_ERI_Ix5y_D2y_S_S = I_ERI_Kx6y_Py_S_S+ABY*I_ERI_Ix5y_Py_S_S;
  Double I_ERI_Ix4yz_D2y_S_S = I_ERI_Kx5yz_Py_S_S+ABY*I_ERI_Ix4yz_Py_S_S;
  Double I_ERI_Ix3y2z_D2y_S_S = I_ERI_Kx4y2z_Py_S_S+ABY*I_ERI_Ix3y2z_Py_S_S;
  Double I_ERI_Ix2y3z_D2y_S_S = I_ERI_Kx3y3z_Py_S_S+ABY*I_ERI_Ix2y3z_Py_S_S;
  Double I_ERI_Ixy4z_D2y_S_S = I_ERI_Kx2y4z_Py_S_S+ABY*I_ERI_Ixy4z_Py_S_S;
  Double I_ERI_Ix5z_D2y_S_S = I_ERI_Kxy5z_Py_S_S+ABY*I_ERI_Ix5z_Py_S_S;
  Double I_ERI_I6y_D2y_S_S = I_ERI_K7y_Py_S_S+ABY*I_ERI_I6y_Py_S_S;
  Double I_ERI_I5yz_D2y_S_S = I_ERI_K6yz_Py_S_S+ABY*I_ERI_I5yz_Py_S_S;
  Double I_ERI_I4y2z_D2y_S_S = I_ERI_K5y2z_Py_S_S+ABY*I_ERI_I4y2z_Py_S_S;
  Double I_ERI_I3y3z_D2y_S_S = I_ERI_K4y3z_Py_S_S+ABY*I_ERI_I3y3z_Py_S_S;
  Double I_ERI_I2y4z_D2y_S_S = I_ERI_K3y4z_Py_S_S+ABY*I_ERI_I2y4z_Py_S_S;
  Double I_ERI_Iy5z_D2y_S_S = I_ERI_K2y5z_Py_S_S+ABY*I_ERI_Iy5z_Py_S_S;
  Double I_ERI_I6z_D2y_S_S = I_ERI_Ky6z_Py_S_S+ABY*I_ERI_I6z_Py_S_S;
  Double I_ERI_I5xy_D2z_S_S = I_ERI_K5xyz_Pz_S_S+ABZ*I_ERI_I5xy_Pz_S_S;
  Double I_ERI_I5xz_D2z_S_S = I_ERI_K5x2z_Pz_S_S+ABZ*I_ERI_I5xz_Pz_S_S;
  Double I_ERI_I4x2y_D2z_S_S = I_ERI_K4x2yz_Pz_S_S+ABZ*I_ERI_I4x2y_Pz_S_S;
  Double I_ERI_I4xyz_D2z_S_S = I_ERI_K4xy2z_Pz_S_S+ABZ*I_ERI_I4xyz_Pz_S_S;
  Double I_ERI_I4x2z_D2z_S_S = I_ERI_K4x3z_Pz_S_S+ABZ*I_ERI_I4x2z_Pz_S_S;
  Double I_ERI_I3x3y_D2z_S_S = I_ERI_K3x3yz_Pz_S_S+ABZ*I_ERI_I3x3y_Pz_S_S;
  Double I_ERI_I3x2yz_D2z_S_S = I_ERI_K3x2y2z_Pz_S_S+ABZ*I_ERI_I3x2yz_Pz_S_S;
  Double I_ERI_I3xy2z_D2z_S_S = I_ERI_K3xy3z_Pz_S_S+ABZ*I_ERI_I3xy2z_Pz_S_S;
  Double I_ERI_I3x3z_D2z_S_S = I_ERI_K3x4z_Pz_S_S+ABZ*I_ERI_I3x3z_Pz_S_S;
  Double I_ERI_I2x4y_D2z_S_S = I_ERI_K2x4yz_Pz_S_S+ABZ*I_ERI_I2x4y_Pz_S_S;
  Double I_ERI_I2x3yz_D2z_S_S = I_ERI_K2x3y2z_Pz_S_S+ABZ*I_ERI_I2x3yz_Pz_S_S;
  Double I_ERI_I2x2y2z_D2z_S_S = I_ERI_K2x2y3z_Pz_S_S+ABZ*I_ERI_I2x2y2z_Pz_S_S;
  Double I_ERI_I2xy3z_D2z_S_S = I_ERI_K2xy4z_Pz_S_S+ABZ*I_ERI_I2xy3z_Pz_S_S;
  Double I_ERI_I2x4z_D2z_S_S = I_ERI_K2x5z_Pz_S_S+ABZ*I_ERI_I2x4z_Pz_S_S;
  Double I_ERI_Ix5y_D2z_S_S = I_ERI_Kx5yz_Pz_S_S+ABZ*I_ERI_Ix5y_Pz_S_S;
  Double I_ERI_Ix4yz_D2z_S_S = I_ERI_Kx4y2z_Pz_S_S+ABZ*I_ERI_Ix4yz_Pz_S_S;
  Double I_ERI_Ix3y2z_D2z_S_S = I_ERI_Kx3y3z_Pz_S_S+ABZ*I_ERI_Ix3y2z_Pz_S_S;
  Double I_ERI_Ix2y3z_D2z_S_S = I_ERI_Kx2y4z_Pz_S_S+ABZ*I_ERI_Ix2y3z_Pz_S_S;
  Double I_ERI_Ixy4z_D2z_S_S = I_ERI_Kxy5z_Pz_S_S+ABZ*I_ERI_Ixy4z_Pz_S_S;
  Double I_ERI_Ix5z_D2z_S_S = I_ERI_Kx6z_Pz_S_S+ABZ*I_ERI_Ix5z_Pz_S_S;
  Double I_ERI_I5yz_D2z_S_S = I_ERI_K5y2z_Pz_S_S+ABZ*I_ERI_I5yz_Pz_S_S;
  Double I_ERI_I4y2z_D2z_S_S = I_ERI_K4y3z_Pz_S_S+ABZ*I_ERI_I4y2z_Pz_S_S;
  Double I_ERI_I3y3z_D2z_S_S = I_ERI_K3y4z_Pz_S_S+ABZ*I_ERI_I3y3z_Pz_S_S;
  Double I_ERI_I2y4z_D2z_S_S = I_ERI_K2y5z_Pz_S_S+ABZ*I_ERI_I2y4z_Pz_S_S;
  Double I_ERI_Iy5z_D2z_S_S = I_ERI_Ky6z_Pz_S_S+ABZ*I_ERI_Iy5z_Pz_S_S;
  Double I_ERI_I6z_D2z_S_S = I_ERI_K7z_Pz_S_S+ABZ*I_ERI_I6z_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_F_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 67 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_D_S_S
   * RHS shell quartet name: SQ_ERI_H_D_S_S
   ************************************************************/
  Double I_ERI_H5x_F3x_S_S = I_ERI_I6x_D2x_S_S+ABX*I_ERI_H5x_D2x_S_S;
  Double I_ERI_H4xy_F3x_S_S = I_ERI_I5xy_D2x_S_S+ABX*I_ERI_H4xy_D2x_S_S;
  Double I_ERI_H4xz_F3x_S_S = I_ERI_I5xz_D2x_S_S+ABX*I_ERI_H4xz_D2x_S_S;
  Double I_ERI_H3x2y_F3x_S_S = I_ERI_I4x2y_D2x_S_S+ABX*I_ERI_H3x2y_D2x_S_S;
  Double I_ERI_H3xyz_F3x_S_S = I_ERI_I4xyz_D2x_S_S+ABX*I_ERI_H3xyz_D2x_S_S;
  Double I_ERI_H3x2z_F3x_S_S = I_ERI_I4x2z_D2x_S_S+ABX*I_ERI_H3x2z_D2x_S_S;
  Double I_ERI_H2x3y_F3x_S_S = I_ERI_I3x3y_D2x_S_S+ABX*I_ERI_H2x3y_D2x_S_S;
  Double I_ERI_H2x2yz_F3x_S_S = I_ERI_I3x2yz_D2x_S_S+ABX*I_ERI_H2x2yz_D2x_S_S;
  Double I_ERI_H2xy2z_F3x_S_S = I_ERI_I3xy2z_D2x_S_S+ABX*I_ERI_H2xy2z_D2x_S_S;
  Double I_ERI_H2x3z_F3x_S_S = I_ERI_I3x3z_D2x_S_S+ABX*I_ERI_H2x3z_D2x_S_S;
  Double I_ERI_Hx4y_F3x_S_S = I_ERI_I2x4y_D2x_S_S+ABX*I_ERI_Hx4y_D2x_S_S;
  Double I_ERI_Hx3yz_F3x_S_S = I_ERI_I2x3yz_D2x_S_S+ABX*I_ERI_Hx3yz_D2x_S_S;
  Double I_ERI_Hx2y2z_F3x_S_S = I_ERI_I2x2y2z_D2x_S_S+ABX*I_ERI_Hx2y2z_D2x_S_S;
  Double I_ERI_Hxy3z_F3x_S_S = I_ERI_I2xy3z_D2x_S_S+ABX*I_ERI_Hxy3z_D2x_S_S;
  Double I_ERI_Hx4z_F3x_S_S = I_ERI_I2x4z_D2x_S_S+ABX*I_ERI_Hx4z_D2x_S_S;
  Double I_ERI_H5y_F3x_S_S = I_ERI_Ix5y_D2x_S_S+ABX*I_ERI_H5y_D2x_S_S;
  Double I_ERI_H4yz_F3x_S_S = I_ERI_Ix4yz_D2x_S_S+ABX*I_ERI_H4yz_D2x_S_S;
  Double I_ERI_H3y2z_F3x_S_S = I_ERI_Ix3y2z_D2x_S_S+ABX*I_ERI_H3y2z_D2x_S_S;
  Double I_ERI_H2y3z_F3x_S_S = I_ERI_Ix2y3z_D2x_S_S+ABX*I_ERI_H2y3z_D2x_S_S;
  Double I_ERI_Hy4z_F3x_S_S = I_ERI_Ixy4z_D2x_S_S+ABX*I_ERI_Hy4z_D2x_S_S;
  Double I_ERI_H5z_F3x_S_S = I_ERI_Ix5z_D2x_S_S+ABX*I_ERI_H5z_D2x_S_S;
  Double I_ERI_H4xy_F2xy_S_S = I_ERI_I4x2y_D2x_S_S+ABY*I_ERI_H4xy_D2x_S_S;
  Double I_ERI_H4xz_F2xy_S_S = I_ERI_I4xyz_D2x_S_S+ABY*I_ERI_H4xz_D2x_S_S;
  Double I_ERI_H3x2y_F2xy_S_S = I_ERI_I3x3y_D2x_S_S+ABY*I_ERI_H3x2y_D2x_S_S;
  Double I_ERI_H3xyz_F2xy_S_S = I_ERI_I3x2yz_D2x_S_S+ABY*I_ERI_H3xyz_D2x_S_S;
  Double I_ERI_H3x2z_F2xy_S_S = I_ERI_I3xy2z_D2x_S_S+ABY*I_ERI_H3x2z_D2x_S_S;
  Double I_ERI_H2x3y_F2xy_S_S = I_ERI_I2x4y_D2x_S_S+ABY*I_ERI_H2x3y_D2x_S_S;
  Double I_ERI_H2x2yz_F2xy_S_S = I_ERI_I2x3yz_D2x_S_S+ABY*I_ERI_H2x2yz_D2x_S_S;
  Double I_ERI_H2xy2z_F2xy_S_S = I_ERI_I2x2y2z_D2x_S_S+ABY*I_ERI_H2xy2z_D2x_S_S;
  Double I_ERI_H2x3z_F2xy_S_S = I_ERI_I2xy3z_D2x_S_S+ABY*I_ERI_H2x3z_D2x_S_S;
  Double I_ERI_Hx4y_F2xy_S_S = I_ERI_Ix5y_D2x_S_S+ABY*I_ERI_Hx4y_D2x_S_S;
  Double I_ERI_Hx3yz_F2xy_S_S = I_ERI_Ix4yz_D2x_S_S+ABY*I_ERI_Hx3yz_D2x_S_S;
  Double I_ERI_Hx2y2z_F2xy_S_S = I_ERI_Ix3y2z_D2x_S_S+ABY*I_ERI_Hx2y2z_D2x_S_S;
  Double I_ERI_Hxy3z_F2xy_S_S = I_ERI_Ix2y3z_D2x_S_S+ABY*I_ERI_Hxy3z_D2x_S_S;
  Double I_ERI_Hx4z_F2xy_S_S = I_ERI_Ixy4z_D2x_S_S+ABY*I_ERI_Hx4z_D2x_S_S;
  Double I_ERI_H5y_F2xy_S_S = I_ERI_I6y_D2x_S_S+ABY*I_ERI_H5y_D2x_S_S;
  Double I_ERI_H4yz_F2xy_S_S = I_ERI_I5yz_D2x_S_S+ABY*I_ERI_H4yz_D2x_S_S;
  Double I_ERI_H3y2z_F2xy_S_S = I_ERI_I4y2z_D2x_S_S+ABY*I_ERI_H3y2z_D2x_S_S;
  Double I_ERI_H2y3z_F2xy_S_S = I_ERI_I3y3z_D2x_S_S+ABY*I_ERI_H2y3z_D2x_S_S;
  Double I_ERI_Hy4z_F2xy_S_S = I_ERI_I2y4z_D2x_S_S+ABY*I_ERI_Hy4z_D2x_S_S;
  Double I_ERI_H5z_F2xy_S_S = I_ERI_Iy5z_D2x_S_S+ABY*I_ERI_H5z_D2x_S_S;
  Double I_ERI_H4xz_F2xz_S_S = I_ERI_I4x2z_D2x_S_S+ABZ*I_ERI_H4xz_D2x_S_S;
  Double I_ERI_H3xyz_F2xz_S_S = I_ERI_I3xy2z_D2x_S_S+ABZ*I_ERI_H3xyz_D2x_S_S;
  Double I_ERI_H3x2z_F2xz_S_S = I_ERI_I3x3z_D2x_S_S+ABZ*I_ERI_H3x2z_D2x_S_S;
  Double I_ERI_H2x2yz_F2xz_S_S = I_ERI_I2x2y2z_D2x_S_S+ABZ*I_ERI_H2x2yz_D2x_S_S;
  Double I_ERI_H2xy2z_F2xz_S_S = I_ERI_I2xy3z_D2x_S_S+ABZ*I_ERI_H2xy2z_D2x_S_S;
  Double I_ERI_H2x3z_F2xz_S_S = I_ERI_I2x4z_D2x_S_S+ABZ*I_ERI_H2x3z_D2x_S_S;
  Double I_ERI_Hx3yz_F2xz_S_S = I_ERI_Ix3y2z_D2x_S_S+ABZ*I_ERI_Hx3yz_D2x_S_S;
  Double I_ERI_Hx2y2z_F2xz_S_S = I_ERI_Ix2y3z_D2x_S_S+ABZ*I_ERI_Hx2y2z_D2x_S_S;
  Double I_ERI_Hxy3z_F2xz_S_S = I_ERI_Ixy4z_D2x_S_S+ABZ*I_ERI_Hxy3z_D2x_S_S;
  Double I_ERI_Hx4z_F2xz_S_S = I_ERI_Ix5z_D2x_S_S+ABZ*I_ERI_Hx4z_D2x_S_S;
  Double I_ERI_H4yz_F2xz_S_S = I_ERI_I4y2z_D2x_S_S+ABZ*I_ERI_H4yz_D2x_S_S;
  Double I_ERI_H3y2z_F2xz_S_S = I_ERI_I3y3z_D2x_S_S+ABZ*I_ERI_H3y2z_D2x_S_S;
  Double I_ERI_H2y3z_F2xz_S_S = I_ERI_I2y4z_D2x_S_S+ABZ*I_ERI_H2y3z_D2x_S_S;
  Double I_ERI_Hy4z_F2xz_S_S = I_ERI_Iy5z_D2x_S_S+ABZ*I_ERI_Hy4z_D2x_S_S;
  Double I_ERI_H5z_F2xz_S_S = I_ERI_I6z_D2x_S_S+ABZ*I_ERI_H5z_D2x_S_S;
  Double I_ERI_H4xz_Fx2y_S_S = I_ERI_I5xz_D2y_S_S+ABX*I_ERI_H4xz_D2y_S_S;
  Double I_ERI_H3xyz_Fx2y_S_S = I_ERI_I4xyz_D2y_S_S+ABX*I_ERI_H3xyz_D2y_S_S;
  Double I_ERI_H3x2z_Fx2y_S_S = I_ERI_I4x2z_D2y_S_S+ABX*I_ERI_H3x2z_D2y_S_S;
  Double I_ERI_H2x2yz_Fx2y_S_S = I_ERI_I3x2yz_D2y_S_S+ABX*I_ERI_H2x2yz_D2y_S_S;
  Double I_ERI_H2xy2z_Fx2y_S_S = I_ERI_I3xy2z_D2y_S_S+ABX*I_ERI_H2xy2z_D2y_S_S;
  Double I_ERI_H2x3z_Fx2y_S_S = I_ERI_I3x3z_D2y_S_S+ABX*I_ERI_H2x3z_D2y_S_S;
  Double I_ERI_Hx3yz_Fx2y_S_S = I_ERI_I2x3yz_D2y_S_S+ABX*I_ERI_Hx3yz_D2y_S_S;
  Double I_ERI_Hx2y2z_Fx2y_S_S = I_ERI_I2x2y2z_D2y_S_S+ABX*I_ERI_Hx2y2z_D2y_S_S;
  Double I_ERI_Hxy3z_Fx2y_S_S = I_ERI_I2xy3z_D2y_S_S+ABX*I_ERI_Hxy3z_D2y_S_S;
  Double I_ERI_Hx4z_Fx2y_S_S = I_ERI_I2x4z_D2y_S_S+ABX*I_ERI_Hx4z_D2y_S_S;
  Double I_ERI_H4yz_Fx2y_S_S = I_ERI_Ix4yz_D2y_S_S+ABX*I_ERI_H4yz_D2y_S_S;
  Double I_ERI_H3y2z_Fx2y_S_S = I_ERI_Ix3y2z_D2y_S_S+ABX*I_ERI_H3y2z_D2y_S_S;
  Double I_ERI_H2y3z_Fx2y_S_S = I_ERI_Ix2y3z_D2y_S_S+ABX*I_ERI_H2y3z_D2y_S_S;
  Double I_ERI_Hy4z_Fx2y_S_S = I_ERI_Ixy4z_D2y_S_S+ABX*I_ERI_Hy4z_D2y_S_S;
  Double I_ERI_H5z_Fx2y_S_S = I_ERI_Ix5z_D2y_S_S+ABX*I_ERI_H5z_D2y_S_S;
  Double I_ERI_H4xy_Fx2z_S_S = I_ERI_I5xy_D2z_S_S+ABX*I_ERI_H4xy_D2z_S_S;
  Double I_ERI_H3x2y_Fx2z_S_S = I_ERI_I4x2y_D2z_S_S+ABX*I_ERI_H3x2y_D2z_S_S;
  Double I_ERI_H3xyz_Fx2z_S_S = I_ERI_I4xyz_D2z_S_S+ABX*I_ERI_H3xyz_D2z_S_S;
  Double I_ERI_H2x3y_Fx2z_S_S = I_ERI_I3x3y_D2z_S_S+ABX*I_ERI_H2x3y_D2z_S_S;
  Double I_ERI_H2x2yz_Fx2z_S_S = I_ERI_I3x2yz_D2z_S_S+ABX*I_ERI_H2x2yz_D2z_S_S;
  Double I_ERI_H2xy2z_Fx2z_S_S = I_ERI_I3xy2z_D2z_S_S+ABX*I_ERI_H2xy2z_D2z_S_S;
  Double I_ERI_Hx4y_Fx2z_S_S = I_ERI_I2x4y_D2z_S_S+ABX*I_ERI_Hx4y_D2z_S_S;
  Double I_ERI_Hx3yz_Fx2z_S_S = I_ERI_I2x3yz_D2z_S_S+ABX*I_ERI_Hx3yz_D2z_S_S;
  Double I_ERI_Hx2y2z_Fx2z_S_S = I_ERI_I2x2y2z_D2z_S_S+ABX*I_ERI_Hx2y2z_D2z_S_S;
  Double I_ERI_Hxy3z_Fx2z_S_S = I_ERI_I2xy3z_D2z_S_S+ABX*I_ERI_Hxy3z_D2z_S_S;
  Double I_ERI_H5y_Fx2z_S_S = I_ERI_Ix5y_D2z_S_S+ABX*I_ERI_H5y_D2z_S_S;
  Double I_ERI_H4yz_Fx2z_S_S = I_ERI_Ix4yz_D2z_S_S+ABX*I_ERI_H4yz_D2z_S_S;
  Double I_ERI_H3y2z_Fx2z_S_S = I_ERI_Ix3y2z_D2z_S_S+ABX*I_ERI_H3y2z_D2z_S_S;
  Double I_ERI_H2y3z_Fx2z_S_S = I_ERI_Ix2y3z_D2z_S_S+ABX*I_ERI_H2y3z_D2z_S_S;
  Double I_ERI_Hy4z_Fx2z_S_S = I_ERI_Ixy4z_D2z_S_S+ABX*I_ERI_Hy4z_D2z_S_S;
  Double I_ERI_H5x_F3y_S_S = I_ERI_I5xy_D2y_S_S+ABY*I_ERI_H5x_D2y_S_S;
  Double I_ERI_H4xy_F3y_S_S = I_ERI_I4x2y_D2y_S_S+ABY*I_ERI_H4xy_D2y_S_S;
  Double I_ERI_H4xz_F3y_S_S = I_ERI_I4xyz_D2y_S_S+ABY*I_ERI_H4xz_D2y_S_S;
  Double I_ERI_H3x2y_F3y_S_S = I_ERI_I3x3y_D2y_S_S+ABY*I_ERI_H3x2y_D2y_S_S;
  Double I_ERI_H3xyz_F3y_S_S = I_ERI_I3x2yz_D2y_S_S+ABY*I_ERI_H3xyz_D2y_S_S;
  Double I_ERI_H3x2z_F3y_S_S = I_ERI_I3xy2z_D2y_S_S+ABY*I_ERI_H3x2z_D2y_S_S;
  Double I_ERI_H2x3y_F3y_S_S = I_ERI_I2x4y_D2y_S_S+ABY*I_ERI_H2x3y_D2y_S_S;
  Double I_ERI_H2x2yz_F3y_S_S = I_ERI_I2x3yz_D2y_S_S+ABY*I_ERI_H2x2yz_D2y_S_S;
  Double I_ERI_H2xy2z_F3y_S_S = I_ERI_I2x2y2z_D2y_S_S+ABY*I_ERI_H2xy2z_D2y_S_S;
  Double I_ERI_H2x3z_F3y_S_S = I_ERI_I2xy3z_D2y_S_S+ABY*I_ERI_H2x3z_D2y_S_S;
  Double I_ERI_Hx4y_F3y_S_S = I_ERI_Ix5y_D2y_S_S+ABY*I_ERI_Hx4y_D2y_S_S;
  Double I_ERI_Hx3yz_F3y_S_S = I_ERI_Ix4yz_D2y_S_S+ABY*I_ERI_Hx3yz_D2y_S_S;
  Double I_ERI_Hx2y2z_F3y_S_S = I_ERI_Ix3y2z_D2y_S_S+ABY*I_ERI_Hx2y2z_D2y_S_S;
  Double I_ERI_Hxy3z_F3y_S_S = I_ERI_Ix2y3z_D2y_S_S+ABY*I_ERI_Hxy3z_D2y_S_S;
  Double I_ERI_Hx4z_F3y_S_S = I_ERI_Ixy4z_D2y_S_S+ABY*I_ERI_Hx4z_D2y_S_S;
  Double I_ERI_H5y_F3y_S_S = I_ERI_I6y_D2y_S_S+ABY*I_ERI_H5y_D2y_S_S;
  Double I_ERI_H4yz_F3y_S_S = I_ERI_I5yz_D2y_S_S+ABY*I_ERI_H4yz_D2y_S_S;
  Double I_ERI_H3y2z_F3y_S_S = I_ERI_I4y2z_D2y_S_S+ABY*I_ERI_H3y2z_D2y_S_S;
  Double I_ERI_H2y3z_F3y_S_S = I_ERI_I3y3z_D2y_S_S+ABY*I_ERI_H2y3z_D2y_S_S;
  Double I_ERI_Hy4z_F3y_S_S = I_ERI_I2y4z_D2y_S_S+ABY*I_ERI_Hy4z_D2y_S_S;
  Double I_ERI_H5z_F3y_S_S = I_ERI_Iy5z_D2y_S_S+ABY*I_ERI_H5z_D2y_S_S;
  Double I_ERI_H4xz_F2yz_S_S = I_ERI_I4x2z_D2y_S_S+ABZ*I_ERI_H4xz_D2y_S_S;
  Double I_ERI_H3xyz_F2yz_S_S = I_ERI_I3xy2z_D2y_S_S+ABZ*I_ERI_H3xyz_D2y_S_S;
  Double I_ERI_H3x2z_F2yz_S_S = I_ERI_I3x3z_D2y_S_S+ABZ*I_ERI_H3x2z_D2y_S_S;
  Double I_ERI_H2x2yz_F2yz_S_S = I_ERI_I2x2y2z_D2y_S_S+ABZ*I_ERI_H2x2yz_D2y_S_S;
  Double I_ERI_H2xy2z_F2yz_S_S = I_ERI_I2xy3z_D2y_S_S+ABZ*I_ERI_H2xy2z_D2y_S_S;
  Double I_ERI_H2x3z_F2yz_S_S = I_ERI_I2x4z_D2y_S_S+ABZ*I_ERI_H2x3z_D2y_S_S;
  Double I_ERI_Hx3yz_F2yz_S_S = I_ERI_Ix3y2z_D2y_S_S+ABZ*I_ERI_Hx3yz_D2y_S_S;
  Double I_ERI_Hx2y2z_F2yz_S_S = I_ERI_Ix2y3z_D2y_S_S+ABZ*I_ERI_Hx2y2z_D2y_S_S;
  Double I_ERI_Hxy3z_F2yz_S_S = I_ERI_Ixy4z_D2y_S_S+ABZ*I_ERI_Hxy3z_D2y_S_S;
  Double I_ERI_Hx4z_F2yz_S_S = I_ERI_Ix5z_D2y_S_S+ABZ*I_ERI_Hx4z_D2y_S_S;
  Double I_ERI_H4yz_F2yz_S_S = I_ERI_I4y2z_D2y_S_S+ABZ*I_ERI_H4yz_D2y_S_S;
  Double I_ERI_H3y2z_F2yz_S_S = I_ERI_I3y3z_D2y_S_S+ABZ*I_ERI_H3y2z_D2y_S_S;
  Double I_ERI_H2y3z_F2yz_S_S = I_ERI_I2y4z_D2y_S_S+ABZ*I_ERI_H2y3z_D2y_S_S;
  Double I_ERI_Hy4z_F2yz_S_S = I_ERI_Iy5z_D2y_S_S+ABZ*I_ERI_Hy4z_D2y_S_S;
  Double I_ERI_H5z_F2yz_S_S = I_ERI_I6z_D2y_S_S+ABZ*I_ERI_H5z_D2y_S_S;
  Double I_ERI_H5x_F3z_S_S = I_ERI_I5xz_D2z_S_S+ABZ*I_ERI_H5x_D2z_S_S;
  Double I_ERI_H4xy_F3z_S_S = I_ERI_I4xyz_D2z_S_S+ABZ*I_ERI_H4xy_D2z_S_S;
  Double I_ERI_H4xz_F3z_S_S = I_ERI_I4x2z_D2z_S_S+ABZ*I_ERI_H4xz_D2z_S_S;
  Double I_ERI_H3x2y_F3z_S_S = I_ERI_I3x2yz_D2z_S_S+ABZ*I_ERI_H3x2y_D2z_S_S;
  Double I_ERI_H3xyz_F3z_S_S = I_ERI_I3xy2z_D2z_S_S+ABZ*I_ERI_H3xyz_D2z_S_S;
  Double I_ERI_H3x2z_F3z_S_S = I_ERI_I3x3z_D2z_S_S+ABZ*I_ERI_H3x2z_D2z_S_S;
  Double I_ERI_H2x3y_F3z_S_S = I_ERI_I2x3yz_D2z_S_S+ABZ*I_ERI_H2x3y_D2z_S_S;
  Double I_ERI_H2x2yz_F3z_S_S = I_ERI_I2x2y2z_D2z_S_S+ABZ*I_ERI_H2x2yz_D2z_S_S;
  Double I_ERI_H2xy2z_F3z_S_S = I_ERI_I2xy3z_D2z_S_S+ABZ*I_ERI_H2xy2z_D2z_S_S;
  Double I_ERI_H2x3z_F3z_S_S = I_ERI_I2x4z_D2z_S_S+ABZ*I_ERI_H2x3z_D2z_S_S;
  Double I_ERI_Hx4y_F3z_S_S = I_ERI_Ix4yz_D2z_S_S+ABZ*I_ERI_Hx4y_D2z_S_S;
  Double I_ERI_Hx3yz_F3z_S_S = I_ERI_Ix3y2z_D2z_S_S+ABZ*I_ERI_Hx3yz_D2z_S_S;
  Double I_ERI_Hx2y2z_F3z_S_S = I_ERI_Ix2y3z_D2z_S_S+ABZ*I_ERI_Hx2y2z_D2z_S_S;
  Double I_ERI_Hxy3z_F3z_S_S = I_ERI_Ixy4z_D2z_S_S+ABZ*I_ERI_Hxy3z_D2z_S_S;
  Double I_ERI_Hx4z_F3z_S_S = I_ERI_Ix5z_D2z_S_S+ABZ*I_ERI_Hx4z_D2z_S_S;
  Double I_ERI_H5y_F3z_S_S = I_ERI_I5yz_D2z_S_S+ABZ*I_ERI_H5y_D2z_S_S;
  Double I_ERI_H4yz_F3z_S_S = I_ERI_I4y2z_D2z_S_S+ABZ*I_ERI_H4yz_D2z_S_S;
  Double I_ERI_H3y2z_F3z_S_S = I_ERI_I3y3z_D2z_S_S+ABZ*I_ERI_H3y2z_D2z_S_S;
  Double I_ERI_H2y3z_F3z_S_S = I_ERI_I2y4z_D2z_S_S+ABZ*I_ERI_H2y3z_D2z_S_S;
  Double I_ERI_Hy4z_F3z_S_S = I_ERI_Iy5z_D2z_S_S+ABZ*I_ERI_Hy4z_D2z_S_S;
  Double I_ERI_H5z_F3z_S_S = I_ERI_I6z_D2z_S_S+ABZ*I_ERI_H5z_D2z_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_G_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_F_S_S
   * RHS shell quartet name: SQ_ERI_G_F_S_S
   ************************************************************/
  abcd[0] = I_ERI_H5x_F3x_S_S+ABX*I_ERI_G4x_F3x_S_S;
  abcd[1] = I_ERI_H4xy_F3x_S_S+ABX*I_ERI_G3xy_F3x_S_S;
  abcd[2] = I_ERI_H4xz_F3x_S_S+ABX*I_ERI_G3xz_F3x_S_S;
  abcd[3] = I_ERI_H3x2y_F3x_S_S+ABX*I_ERI_G2x2y_F3x_S_S;
  abcd[4] = I_ERI_H3xyz_F3x_S_S+ABX*I_ERI_G2xyz_F3x_S_S;
  abcd[5] = I_ERI_H3x2z_F3x_S_S+ABX*I_ERI_G2x2z_F3x_S_S;
  abcd[6] = I_ERI_H2x3y_F3x_S_S+ABX*I_ERI_Gx3y_F3x_S_S;
  abcd[7] = I_ERI_H2x2yz_F3x_S_S+ABX*I_ERI_Gx2yz_F3x_S_S;
  abcd[8] = I_ERI_H2xy2z_F3x_S_S+ABX*I_ERI_Gxy2z_F3x_S_S;
  abcd[9] = I_ERI_H2x3z_F3x_S_S+ABX*I_ERI_Gx3z_F3x_S_S;
  abcd[10] = I_ERI_Hx4y_F3x_S_S+ABX*I_ERI_G4y_F3x_S_S;
  abcd[11] = I_ERI_Hx3yz_F3x_S_S+ABX*I_ERI_G3yz_F3x_S_S;
  abcd[12] = I_ERI_Hx2y2z_F3x_S_S+ABX*I_ERI_G2y2z_F3x_S_S;
  abcd[13] = I_ERI_Hxy3z_F3x_S_S+ABX*I_ERI_Gy3z_F3x_S_S;
  abcd[14] = I_ERI_Hx4z_F3x_S_S+ABX*I_ERI_G4z_F3x_S_S;
  abcd[15] = I_ERI_H4xy_F3x_S_S+ABY*I_ERI_G4x_F3x_S_S;
  abcd[16] = I_ERI_H3x2y_F3x_S_S+ABY*I_ERI_G3xy_F3x_S_S;
  abcd[17] = I_ERI_H3xyz_F3x_S_S+ABY*I_ERI_G3xz_F3x_S_S;
  abcd[18] = I_ERI_H2x3y_F3x_S_S+ABY*I_ERI_G2x2y_F3x_S_S;
  abcd[19] = I_ERI_H2x2yz_F3x_S_S+ABY*I_ERI_G2xyz_F3x_S_S;
  abcd[20] = I_ERI_H2xy2z_F3x_S_S+ABY*I_ERI_G2x2z_F3x_S_S;
  abcd[21] = I_ERI_Hx4y_F3x_S_S+ABY*I_ERI_Gx3y_F3x_S_S;
  abcd[22] = I_ERI_Hx3yz_F3x_S_S+ABY*I_ERI_Gx2yz_F3x_S_S;
  abcd[23] = I_ERI_Hx2y2z_F3x_S_S+ABY*I_ERI_Gxy2z_F3x_S_S;
  abcd[24] = I_ERI_Hxy3z_F3x_S_S+ABY*I_ERI_Gx3z_F3x_S_S;
  abcd[25] = I_ERI_H5y_F3x_S_S+ABY*I_ERI_G4y_F3x_S_S;
  abcd[26] = I_ERI_H4yz_F3x_S_S+ABY*I_ERI_G3yz_F3x_S_S;
  abcd[27] = I_ERI_H3y2z_F3x_S_S+ABY*I_ERI_G2y2z_F3x_S_S;
  abcd[28] = I_ERI_H2y3z_F3x_S_S+ABY*I_ERI_Gy3z_F3x_S_S;
  abcd[29] = I_ERI_Hy4z_F3x_S_S+ABY*I_ERI_G4z_F3x_S_S;
  abcd[30] = I_ERI_H4xz_F3x_S_S+ABZ*I_ERI_G4x_F3x_S_S;
  abcd[31] = I_ERI_H3xyz_F3x_S_S+ABZ*I_ERI_G3xy_F3x_S_S;
  abcd[32] = I_ERI_H3x2z_F3x_S_S+ABZ*I_ERI_G3xz_F3x_S_S;
  abcd[33] = I_ERI_H2x2yz_F3x_S_S+ABZ*I_ERI_G2x2y_F3x_S_S;
  abcd[34] = I_ERI_H2xy2z_F3x_S_S+ABZ*I_ERI_G2xyz_F3x_S_S;
  abcd[35] = I_ERI_H2x3z_F3x_S_S+ABZ*I_ERI_G2x2z_F3x_S_S;
  abcd[36] = I_ERI_Hx3yz_F3x_S_S+ABZ*I_ERI_Gx3y_F3x_S_S;
  abcd[37] = I_ERI_Hx2y2z_F3x_S_S+ABZ*I_ERI_Gx2yz_F3x_S_S;
  abcd[38] = I_ERI_Hxy3z_F3x_S_S+ABZ*I_ERI_Gxy2z_F3x_S_S;
  abcd[39] = I_ERI_Hx4z_F3x_S_S+ABZ*I_ERI_Gx3z_F3x_S_S;
  abcd[40] = I_ERI_H4yz_F3x_S_S+ABZ*I_ERI_G4y_F3x_S_S;
  abcd[41] = I_ERI_H3y2z_F3x_S_S+ABZ*I_ERI_G3yz_F3x_S_S;
  abcd[42] = I_ERI_H2y3z_F3x_S_S+ABZ*I_ERI_G2y2z_F3x_S_S;
  abcd[43] = I_ERI_Hy4z_F3x_S_S+ABZ*I_ERI_Gy3z_F3x_S_S;
  abcd[44] = I_ERI_H5z_F3x_S_S+ABZ*I_ERI_G4z_F3x_S_S;
  abcd[45] = I_ERI_H4xy_F2xy_S_S+ABY*I_ERI_G4x_F2xy_S_S;
  abcd[46] = I_ERI_H3x2y_F2xy_S_S+ABY*I_ERI_G3xy_F2xy_S_S;
  abcd[47] = I_ERI_H3xyz_F2xy_S_S+ABY*I_ERI_G3xz_F2xy_S_S;
  abcd[48] = I_ERI_H2x3y_F2xy_S_S+ABY*I_ERI_G2x2y_F2xy_S_S;
  abcd[49] = I_ERI_H2x2yz_F2xy_S_S+ABY*I_ERI_G2xyz_F2xy_S_S;
  abcd[50] = I_ERI_H2xy2z_F2xy_S_S+ABY*I_ERI_G2x2z_F2xy_S_S;
  abcd[51] = I_ERI_Hx4y_F2xy_S_S+ABY*I_ERI_Gx3y_F2xy_S_S;
  abcd[52] = I_ERI_Hx3yz_F2xy_S_S+ABY*I_ERI_Gx2yz_F2xy_S_S;
  abcd[53] = I_ERI_Hx2y2z_F2xy_S_S+ABY*I_ERI_Gxy2z_F2xy_S_S;
  abcd[54] = I_ERI_Hxy3z_F2xy_S_S+ABY*I_ERI_Gx3z_F2xy_S_S;
  abcd[55] = I_ERI_H5y_F2xy_S_S+ABY*I_ERI_G4y_F2xy_S_S;
  abcd[56] = I_ERI_H4yz_F2xy_S_S+ABY*I_ERI_G3yz_F2xy_S_S;
  abcd[57] = I_ERI_H3y2z_F2xy_S_S+ABY*I_ERI_G2y2z_F2xy_S_S;
  abcd[58] = I_ERI_H2y3z_F2xy_S_S+ABY*I_ERI_Gy3z_F2xy_S_S;
  abcd[59] = I_ERI_Hy4z_F2xy_S_S+ABY*I_ERI_G4z_F2xy_S_S;
  abcd[60] = I_ERI_H4xz_F2xy_S_S+ABZ*I_ERI_G4x_F2xy_S_S;
  abcd[61] = I_ERI_H3xyz_F2xy_S_S+ABZ*I_ERI_G3xy_F2xy_S_S;
  abcd[62] = I_ERI_H3x2z_F2xy_S_S+ABZ*I_ERI_G3xz_F2xy_S_S;
  abcd[63] = I_ERI_H2x2yz_F2xy_S_S+ABZ*I_ERI_G2x2y_F2xy_S_S;
  abcd[64] = I_ERI_H2xy2z_F2xy_S_S+ABZ*I_ERI_G2xyz_F2xy_S_S;
  abcd[65] = I_ERI_H2x3z_F2xy_S_S+ABZ*I_ERI_G2x2z_F2xy_S_S;
  abcd[66] = I_ERI_Hx3yz_F2xy_S_S+ABZ*I_ERI_Gx3y_F2xy_S_S;
  abcd[67] = I_ERI_Hx2y2z_F2xy_S_S+ABZ*I_ERI_Gx2yz_F2xy_S_S;
  abcd[68] = I_ERI_Hxy3z_F2xy_S_S+ABZ*I_ERI_Gxy2z_F2xy_S_S;
  abcd[69] = I_ERI_Hx4z_F2xy_S_S+ABZ*I_ERI_Gx3z_F2xy_S_S;
  abcd[70] = I_ERI_H4yz_F2xy_S_S+ABZ*I_ERI_G4y_F2xy_S_S;
  abcd[71] = I_ERI_H3y2z_F2xy_S_S+ABZ*I_ERI_G3yz_F2xy_S_S;
  abcd[72] = I_ERI_H2y3z_F2xy_S_S+ABZ*I_ERI_G2y2z_F2xy_S_S;
  abcd[73] = I_ERI_Hy4z_F2xy_S_S+ABZ*I_ERI_Gy3z_F2xy_S_S;
  abcd[74] = I_ERI_H5z_F2xy_S_S+ABZ*I_ERI_G4z_F2xy_S_S;
  abcd[75] = I_ERI_H4xz_F2xz_S_S+ABZ*I_ERI_G4x_F2xz_S_S;
  abcd[76] = I_ERI_H3xyz_F2xz_S_S+ABZ*I_ERI_G3xy_F2xz_S_S;
  abcd[77] = I_ERI_H3x2z_F2xz_S_S+ABZ*I_ERI_G3xz_F2xz_S_S;
  abcd[78] = I_ERI_H2x2yz_F2xz_S_S+ABZ*I_ERI_G2x2y_F2xz_S_S;
  abcd[79] = I_ERI_H2xy2z_F2xz_S_S+ABZ*I_ERI_G2xyz_F2xz_S_S;
  abcd[80] = I_ERI_H2x3z_F2xz_S_S+ABZ*I_ERI_G2x2z_F2xz_S_S;
  abcd[81] = I_ERI_Hx3yz_F2xz_S_S+ABZ*I_ERI_Gx3y_F2xz_S_S;
  abcd[82] = I_ERI_Hx2y2z_F2xz_S_S+ABZ*I_ERI_Gx2yz_F2xz_S_S;
  abcd[83] = I_ERI_Hxy3z_F2xz_S_S+ABZ*I_ERI_Gxy2z_F2xz_S_S;
  abcd[84] = I_ERI_Hx4z_F2xz_S_S+ABZ*I_ERI_Gx3z_F2xz_S_S;
  abcd[85] = I_ERI_H4yz_F2xz_S_S+ABZ*I_ERI_G4y_F2xz_S_S;
  abcd[86] = I_ERI_H3y2z_F2xz_S_S+ABZ*I_ERI_G3yz_F2xz_S_S;
  abcd[87] = I_ERI_H2y3z_F2xz_S_S+ABZ*I_ERI_G2y2z_F2xz_S_S;
  abcd[88] = I_ERI_Hy4z_F2xz_S_S+ABZ*I_ERI_Gy3z_F2xz_S_S;
  abcd[89] = I_ERI_H5z_F2xz_S_S+ABZ*I_ERI_G4z_F2xz_S_S;
  abcd[90] = I_ERI_H5x_F3y_S_S+ABX*I_ERI_G4x_F3y_S_S;
  abcd[91] = I_ERI_H4xy_F3y_S_S+ABX*I_ERI_G3xy_F3y_S_S;
  abcd[92] = I_ERI_H4xz_F3y_S_S+ABX*I_ERI_G3xz_F3y_S_S;
  abcd[93] = I_ERI_H3x2y_F3y_S_S+ABX*I_ERI_G2x2y_F3y_S_S;
  abcd[94] = I_ERI_H3xyz_F3y_S_S+ABX*I_ERI_G2xyz_F3y_S_S;
  abcd[95] = I_ERI_H3x2z_F3y_S_S+ABX*I_ERI_G2x2z_F3y_S_S;
  abcd[96] = I_ERI_H2x3y_F3y_S_S+ABX*I_ERI_Gx3y_F3y_S_S;
  abcd[97] = I_ERI_H2x2yz_F3y_S_S+ABX*I_ERI_Gx2yz_F3y_S_S;
  abcd[98] = I_ERI_H2xy2z_F3y_S_S+ABX*I_ERI_Gxy2z_F3y_S_S;
  abcd[99] = I_ERI_H2x3z_F3y_S_S+ABX*I_ERI_Gx3z_F3y_S_S;
  abcd[100] = I_ERI_Hx4y_F3y_S_S+ABX*I_ERI_G4y_F3y_S_S;
  abcd[101] = I_ERI_Hx3yz_F3y_S_S+ABX*I_ERI_G3yz_F3y_S_S;
  abcd[102] = I_ERI_Hx2y2z_F3y_S_S+ABX*I_ERI_G2y2z_F3y_S_S;
  abcd[103] = I_ERI_Hxy3z_F3y_S_S+ABX*I_ERI_Gy3z_F3y_S_S;
  abcd[104] = I_ERI_Hx4z_F3y_S_S+ABX*I_ERI_G4z_F3y_S_S;
  abcd[105] = I_ERI_H4xz_Fx2y_S_S+ABZ*I_ERI_G4x_Fx2y_S_S;
  abcd[106] = I_ERI_H3xyz_Fx2y_S_S+ABZ*I_ERI_G3xy_Fx2y_S_S;
  abcd[107] = I_ERI_H3x2z_Fx2y_S_S+ABZ*I_ERI_G3xz_Fx2y_S_S;
  abcd[108] = I_ERI_H2x2yz_Fx2y_S_S+ABZ*I_ERI_G2x2y_Fx2y_S_S;
  abcd[109] = I_ERI_H2xy2z_Fx2y_S_S+ABZ*I_ERI_G2xyz_Fx2y_S_S;
  abcd[110] = I_ERI_H2x3z_Fx2y_S_S+ABZ*I_ERI_G2x2z_Fx2y_S_S;
  abcd[111] = I_ERI_Hx3yz_Fx2y_S_S+ABZ*I_ERI_Gx3y_Fx2y_S_S;
  abcd[112] = I_ERI_Hx2y2z_Fx2y_S_S+ABZ*I_ERI_Gx2yz_Fx2y_S_S;
  abcd[113] = I_ERI_Hxy3z_Fx2y_S_S+ABZ*I_ERI_Gxy2z_Fx2y_S_S;
  abcd[114] = I_ERI_Hx4z_Fx2y_S_S+ABZ*I_ERI_Gx3z_Fx2y_S_S;
  abcd[115] = I_ERI_H4yz_Fx2y_S_S+ABZ*I_ERI_G4y_Fx2y_S_S;
  abcd[116] = I_ERI_H3y2z_Fx2y_S_S+ABZ*I_ERI_G3yz_Fx2y_S_S;
  abcd[117] = I_ERI_H2y3z_Fx2y_S_S+ABZ*I_ERI_G2y2z_Fx2y_S_S;
  abcd[118] = I_ERI_Hy4z_Fx2y_S_S+ABZ*I_ERI_Gy3z_Fx2y_S_S;
  abcd[119] = I_ERI_H5z_Fx2y_S_S+ABZ*I_ERI_G4z_Fx2y_S_S;
  abcd[120] = I_ERI_H4xy_Fx2z_S_S+ABY*I_ERI_G4x_Fx2z_S_S;
  abcd[121] = I_ERI_H3x2y_Fx2z_S_S+ABY*I_ERI_G3xy_Fx2z_S_S;
  abcd[122] = I_ERI_H3xyz_Fx2z_S_S+ABY*I_ERI_G3xz_Fx2z_S_S;
  abcd[123] = I_ERI_H2x3y_Fx2z_S_S+ABY*I_ERI_G2x2y_Fx2z_S_S;
  abcd[124] = I_ERI_H2x2yz_Fx2z_S_S+ABY*I_ERI_G2xyz_Fx2z_S_S;
  abcd[125] = I_ERI_H2xy2z_Fx2z_S_S+ABY*I_ERI_G2x2z_Fx2z_S_S;
  abcd[126] = I_ERI_Hx4y_Fx2z_S_S+ABY*I_ERI_Gx3y_Fx2z_S_S;
  abcd[127] = I_ERI_Hx3yz_Fx2z_S_S+ABY*I_ERI_Gx2yz_Fx2z_S_S;
  abcd[128] = I_ERI_Hx2y2z_Fx2z_S_S+ABY*I_ERI_Gxy2z_Fx2z_S_S;
  abcd[129] = I_ERI_Hxy3z_Fx2z_S_S+ABY*I_ERI_Gx3z_Fx2z_S_S;
  abcd[130] = I_ERI_H5y_Fx2z_S_S+ABY*I_ERI_G4y_Fx2z_S_S;
  abcd[131] = I_ERI_H4yz_Fx2z_S_S+ABY*I_ERI_G3yz_Fx2z_S_S;
  abcd[132] = I_ERI_H3y2z_Fx2z_S_S+ABY*I_ERI_G2y2z_Fx2z_S_S;
  abcd[133] = I_ERI_H2y3z_Fx2z_S_S+ABY*I_ERI_Gy3z_Fx2z_S_S;
  abcd[134] = I_ERI_Hy4z_Fx2z_S_S+ABY*I_ERI_G4z_Fx2z_S_S;
  abcd[135] = I_ERI_H5x_F3z_S_S+ABX*I_ERI_G4x_F3z_S_S;
  abcd[136] = I_ERI_H4xy_F3z_S_S+ABX*I_ERI_G3xy_F3z_S_S;
  abcd[137] = I_ERI_H4xz_F3z_S_S+ABX*I_ERI_G3xz_F3z_S_S;
  abcd[138] = I_ERI_H3x2y_F3z_S_S+ABX*I_ERI_G2x2y_F3z_S_S;
  abcd[139] = I_ERI_H3xyz_F3z_S_S+ABX*I_ERI_G2xyz_F3z_S_S;
  abcd[140] = I_ERI_H3x2z_F3z_S_S+ABX*I_ERI_G2x2z_F3z_S_S;
  abcd[141] = I_ERI_H2x3y_F3z_S_S+ABX*I_ERI_Gx3y_F3z_S_S;
  abcd[142] = I_ERI_H2x2yz_F3z_S_S+ABX*I_ERI_Gx2yz_F3z_S_S;
  abcd[143] = I_ERI_H2xy2z_F3z_S_S+ABX*I_ERI_Gxy2z_F3z_S_S;
  abcd[144] = I_ERI_H2x3z_F3z_S_S+ABX*I_ERI_Gx3z_F3z_S_S;
  abcd[145] = I_ERI_Hx4y_F3z_S_S+ABX*I_ERI_G4y_F3z_S_S;
  abcd[146] = I_ERI_Hx3yz_F3z_S_S+ABX*I_ERI_G3yz_F3z_S_S;
  abcd[147] = I_ERI_Hx2y2z_F3z_S_S+ABX*I_ERI_G2y2z_F3z_S_S;
  abcd[148] = I_ERI_Hxy3z_F3z_S_S+ABX*I_ERI_Gy3z_F3z_S_S;
  abcd[149] = I_ERI_Hx4z_F3z_S_S+ABX*I_ERI_G4z_F3z_S_S;
  abcd[150] = I_ERI_H4xy_F3y_S_S+ABY*I_ERI_G4x_F3y_S_S;
  abcd[151] = I_ERI_H3x2y_F3y_S_S+ABY*I_ERI_G3xy_F3y_S_S;
  abcd[152] = I_ERI_H3xyz_F3y_S_S+ABY*I_ERI_G3xz_F3y_S_S;
  abcd[153] = I_ERI_H2x3y_F3y_S_S+ABY*I_ERI_G2x2y_F3y_S_S;
  abcd[154] = I_ERI_H2x2yz_F3y_S_S+ABY*I_ERI_G2xyz_F3y_S_S;
  abcd[155] = I_ERI_H2xy2z_F3y_S_S+ABY*I_ERI_G2x2z_F3y_S_S;
  abcd[156] = I_ERI_Hx4y_F3y_S_S+ABY*I_ERI_Gx3y_F3y_S_S;
  abcd[157] = I_ERI_Hx3yz_F3y_S_S+ABY*I_ERI_Gx2yz_F3y_S_S;
  abcd[158] = I_ERI_Hx2y2z_F3y_S_S+ABY*I_ERI_Gxy2z_F3y_S_S;
  abcd[159] = I_ERI_Hxy3z_F3y_S_S+ABY*I_ERI_Gx3z_F3y_S_S;
  abcd[160] = I_ERI_H5y_F3y_S_S+ABY*I_ERI_G4y_F3y_S_S;
  abcd[161] = I_ERI_H4yz_F3y_S_S+ABY*I_ERI_G3yz_F3y_S_S;
  abcd[162] = I_ERI_H3y2z_F3y_S_S+ABY*I_ERI_G2y2z_F3y_S_S;
  abcd[163] = I_ERI_H2y3z_F3y_S_S+ABY*I_ERI_Gy3z_F3y_S_S;
  abcd[164] = I_ERI_Hy4z_F3y_S_S+ABY*I_ERI_G4z_F3y_S_S;
  abcd[165] = I_ERI_H4xz_F3y_S_S+ABZ*I_ERI_G4x_F3y_S_S;
  abcd[166] = I_ERI_H3xyz_F3y_S_S+ABZ*I_ERI_G3xy_F3y_S_S;
  abcd[167] = I_ERI_H3x2z_F3y_S_S+ABZ*I_ERI_G3xz_F3y_S_S;
  abcd[168] = I_ERI_H2x2yz_F3y_S_S+ABZ*I_ERI_G2x2y_F3y_S_S;
  abcd[169] = I_ERI_H2xy2z_F3y_S_S+ABZ*I_ERI_G2xyz_F3y_S_S;
  abcd[170] = I_ERI_H2x3z_F3y_S_S+ABZ*I_ERI_G2x2z_F3y_S_S;
  abcd[171] = I_ERI_Hx3yz_F3y_S_S+ABZ*I_ERI_Gx3y_F3y_S_S;
  abcd[172] = I_ERI_Hx2y2z_F3y_S_S+ABZ*I_ERI_Gx2yz_F3y_S_S;
  abcd[173] = I_ERI_Hxy3z_F3y_S_S+ABZ*I_ERI_Gxy2z_F3y_S_S;
  abcd[174] = I_ERI_Hx4z_F3y_S_S+ABZ*I_ERI_Gx3z_F3y_S_S;
  abcd[175] = I_ERI_H4yz_F3y_S_S+ABZ*I_ERI_G4y_F3y_S_S;
  abcd[176] = I_ERI_H3y2z_F3y_S_S+ABZ*I_ERI_G3yz_F3y_S_S;
  abcd[177] = I_ERI_H2y3z_F3y_S_S+ABZ*I_ERI_G2y2z_F3y_S_S;
  abcd[178] = I_ERI_Hy4z_F3y_S_S+ABZ*I_ERI_Gy3z_F3y_S_S;
  abcd[179] = I_ERI_H5z_F3y_S_S+ABZ*I_ERI_G4z_F3y_S_S;
  abcd[180] = I_ERI_H4xz_F2yz_S_S+ABZ*I_ERI_G4x_F2yz_S_S;
  abcd[181] = I_ERI_H3xyz_F2yz_S_S+ABZ*I_ERI_G3xy_F2yz_S_S;
  abcd[182] = I_ERI_H3x2z_F2yz_S_S+ABZ*I_ERI_G3xz_F2yz_S_S;
  abcd[183] = I_ERI_H2x2yz_F2yz_S_S+ABZ*I_ERI_G2x2y_F2yz_S_S;
  abcd[184] = I_ERI_H2xy2z_F2yz_S_S+ABZ*I_ERI_G2xyz_F2yz_S_S;
  abcd[185] = I_ERI_H2x3z_F2yz_S_S+ABZ*I_ERI_G2x2z_F2yz_S_S;
  abcd[186] = I_ERI_Hx3yz_F2yz_S_S+ABZ*I_ERI_Gx3y_F2yz_S_S;
  abcd[187] = I_ERI_Hx2y2z_F2yz_S_S+ABZ*I_ERI_Gx2yz_F2yz_S_S;
  abcd[188] = I_ERI_Hxy3z_F2yz_S_S+ABZ*I_ERI_Gxy2z_F2yz_S_S;
  abcd[189] = I_ERI_Hx4z_F2yz_S_S+ABZ*I_ERI_Gx3z_F2yz_S_S;
  abcd[190] = I_ERI_H4yz_F2yz_S_S+ABZ*I_ERI_G4y_F2yz_S_S;
  abcd[191] = I_ERI_H3y2z_F2yz_S_S+ABZ*I_ERI_G3yz_F2yz_S_S;
  abcd[192] = I_ERI_H2y3z_F2yz_S_S+ABZ*I_ERI_G2y2z_F2yz_S_S;
  abcd[193] = I_ERI_Hy4z_F2yz_S_S+ABZ*I_ERI_Gy3z_F2yz_S_S;
  abcd[194] = I_ERI_H5z_F2yz_S_S+ABZ*I_ERI_G4z_F2yz_S_S;
  abcd[195] = I_ERI_H4xy_F3z_S_S+ABY*I_ERI_G4x_F3z_S_S;
  abcd[196] = I_ERI_H3x2y_F3z_S_S+ABY*I_ERI_G3xy_F3z_S_S;
  abcd[197] = I_ERI_H3xyz_F3z_S_S+ABY*I_ERI_G3xz_F3z_S_S;
  abcd[198] = I_ERI_H2x3y_F3z_S_S+ABY*I_ERI_G2x2y_F3z_S_S;
  abcd[199] = I_ERI_H2x2yz_F3z_S_S+ABY*I_ERI_G2xyz_F3z_S_S;
  abcd[200] = I_ERI_H2xy2z_F3z_S_S+ABY*I_ERI_G2x2z_F3z_S_S;
  abcd[201] = I_ERI_Hx4y_F3z_S_S+ABY*I_ERI_Gx3y_F3z_S_S;
  abcd[202] = I_ERI_Hx3yz_F3z_S_S+ABY*I_ERI_Gx2yz_F3z_S_S;
  abcd[203] = I_ERI_Hx2y2z_F3z_S_S+ABY*I_ERI_Gxy2z_F3z_S_S;
  abcd[204] = I_ERI_Hxy3z_F3z_S_S+ABY*I_ERI_Gx3z_F3z_S_S;
  abcd[205] = I_ERI_H5y_F3z_S_S+ABY*I_ERI_G4y_F3z_S_S;
  abcd[206] = I_ERI_H4yz_F3z_S_S+ABY*I_ERI_G3yz_F3z_S_S;
  abcd[207] = I_ERI_H3y2z_F3z_S_S+ABY*I_ERI_G2y2z_F3z_S_S;
  abcd[208] = I_ERI_H2y3z_F3z_S_S+ABY*I_ERI_Gy3z_F3z_S_S;
  abcd[209] = I_ERI_Hy4z_F3z_S_S+ABY*I_ERI_G4z_F3z_S_S;
  abcd[210] = I_ERI_H4xz_F3z_S_S+ABZ*I_ERI_G4x_F3z_S_S;
  abcd[211] = I_ERI_H3xyz_F3z_S_S+ABZ*I_ERI_G3xy_F3z_S_S;
  abcd[212] = I_ERI_H3x2z_F3z_S_S+ABZ*I_ERI_G3xz_F3z_S_S;
  abcd[213] = I_ERI_H2x2yz_F3z_S_S+ABZ*I_ERI_G2x2y_F3z_S_S;
  abcd[214] = I_ERI_H2xy2z_F3z_S_S+ABZ*I_ERI_G2xyz_F3z_S_S;
  abcd[215] = I_ERI_H2x3z_F3z_S_S+ABZ*I_ERI_G2x2z_F3z_S_S;
  abcd[216] = I_ERI_Hx3yz_F3z_S_S+ABZ*I_ERI_Gx3y_F3z_S_S;
  abcd[217] = I_ERI_Hx2y2z_F3z_S_S+ABZ*I_ERI_Gx2yz_F3z_S_S;
  abcd[218] = I_ERI_Hxy3z_F3z_S_S+ABZ*I_ERI_Gxy2z_F3z_S_S;
  abcd[219] = I_ERI_Hx4z_F3z_S_S+ABZ*I_ERI_Gx3z_F3z_S_S;
  abcd[220] = I_ERI_H4yz_F3z_S_S+ABZ*I_ERI_G4y_F3z_S_S;
  abcd[221] = I_ERI_H3y2z_F3z_S_S+ABZ*I_ERI_G3yz_F3z_S_S;
  abcd[222] = I_ERI_H2y3z_F3z_S_S+ABZ*I_ERI_G2y2z_F3z_S_S;
  abcd[223] = I_ERI_Hy4z_F3z_S_S+ABZ*I_ERI_Gy3z_F3z_S_S;
  abcd[224] = I_ERI_H5z_F3z_S_S+ABZ*I_ERI_G4z_F3z_S_S;
}
