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

void hgp_os_nai_h_f(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_L8x_S = 0.0E0;
  Double I_NAI_L7xy_S = 0.0E0;
  Double I_NAI_L7xz_S = 0.0E0;
  Double I_NAI_L6x2y_S = 0.0E0;
  Double I_NAI_L6xyz_S = 0.0E0;
  Double I_NAI_L6x2z_S = 0.0E0;
  Double I_NAI_L5x3y_S = 0.0E0;
  Double I_NAI_L5x2yz_S = 0.0E0;
  Double I_NAI_L5xy2z_S = 0.0E0;
  Double I_NAI_L5x3z_S = 0.0E0;
  Double I_NAI_L4x4y_S = 0.0E0;
  Double I_NAI_L4x3yz_S = 0.0E0;
  Double I_NAI_L4x2y2z_S = 0.0E0;
  Double I_NAI_L4xy3z_S = 0.0E0;
  Double I_NAI_L4x4z_S = 0.0E0;
  Double I_NAI_L3x5y_S = 0.0E0;
  Double I_NAI_L3x4yz_S = 0.0E0;
  Double I_NAI_L3x3y2z_S = 0.0E0;
  Double I_NAI_L3x2y3z_S = 0.0E0;
  Double I_NAI_L3xy4z_S = 0.0E0;
  Double I_NAI_L3x5z_S = 0.0E0;
  Double I_NAI_L2x6y_S = 0.0E0;
  Double I_NAI_L2x5yz_S = 0.0E0;
  Double I_NAI_L2x4y2z_S = 0.0E0;
  Double I_NAI_L2x3y3z_S = 0.0E0;
  Double I_NAI_L2x2y4z_S = 0.0E0;
  Double I_NAI_L2xy5z_S = 0.0E0;
  Double I_NAI_L2x6z_S = 0.0E0;
  Double I_NAI_Lx7y_S = 0.0E0;
  Double I_NAI_Lx6yz_S = 0.0E0;
  Double I_NAI_Lx5y2z_S = 0.0E0;
  Double I_NAI_Lx4y3z_S = 0.0E0;
  Double I_NAI_Lx3y4z_S = 0.0E0;
  Double I_NAI_Lx2y5z_S = 0.0E0;
  Double I_NAI_Lxy6z_S = 0.0E0;
  Double I_NAI_Lx7z_S = 0.0E0;
  Double I_NAI_L8y_S = 0.0E0;
  Double I_NAI_L7yz_S = 0.0E0;
  Double I_NAI_L6y2z_S = 0.0E0;
  Double I_NAI_L5y3z_S = 0.0E0;
  Double I_NAI_L4y4z_S = 0.0E0;
  Double I_NAI_L3y5z_S = 0.0E0;
  Double I_NAI_L2y6z_S = 0.0E0;
  Double I_NAI_Ly7z_S = 0.0E0;
  Double I_NAI_L8z_S = 0.0E0;
  Double I_NAI_K7x_S = 0.0E0;
  Double I_NAI_K6xy_S = 0.0E0;
  Double I_NAI_K6xz_S = 0.0E0;
  Double I_NAI_K5x2y_S = 0.0E0;
  Double I_NAI_K5xyz_S = 0.0E0;
  Double I_NAI_K5x2z_S = 0.0E0;
  Double I_NAI_K4x3y_S = 0.0E0;
  Double I_NAI_K4x2yz_S = 0.0E0;
  Double I_NAI_K4xy2z_S = 0.0E0;
  Double I_NAI_K4x3z_S = 0.0E0;
  Double I_NAI_K3x4y_S = 0.0E0;
  Double I_NAI_K3x3yz_S = 0.0E0;
  Double I_NAI_K3x2y2z_S = 0.0E0;
  Double I_NAI_K3xy3z_S = 0.0E0;
  Double I_NAI_K3x4z_S = 0.0E0;
  Double I_NAI_K2x5y_S = 0.0E0;
  Double I_NAI_K2x4yz_S = 0.0E0;
  Double I_NAI_K2x3y2z_S = 0.0E0;
  Double I_NAI_K2x2y3z_S = 0.0E0;
  Double I_NAI_K2xy4z_S = 0.0E0;
  Double I_NAI_K2x5z_S = 0.0E0;
  Double I_NAI_Kx6y_S = 0.0E0;
  Double I_NAI_Kx5yz_S = 0.0E0;
  Double I_NAI_Kx4y2z_S = 0.0E0;
  Double I_NAI_Kx3y3z_S = 0.0E0;
  Double I_NAI_Kx2y4z_S = 0.0E0;
  Double I_NAI_Kxy5z_S = 0.0E0;
  Double I_NAI_Kx6z_S = 0.0E0;
  Double I_NAI_K7y_S = 0.0E0;
  Double I_NAI_K6yz_S = 0.0E0;
  Double I_NAI_K5y2z_S = 0.0E0;
  Double I_NAI_K4y3z_S = 0.0E0;
  Double I_NAI_K3y4z_S = 0.0E0;
  Double I_NAI_K2y5z_S = 0.0E0;
  Double I_NAI_Ky6z_S = 0.0E0;
  Double I_NAI_K7z_S = 0.0E0;
  Double I_NAI_I6x_S = 0.0E0;
  Double I_NAI_I5xy_S = 0.0E0;
  Double I_NAI_I5xz_S = 0.0E0;
  Double I_NAI_I4x2y_S = 0.0E0;
  Double I_NAI_I4xyz_S = 0.0E0;
  Double I_NAI_I4x2z_S = 0.0E0;
  Double I_NAI_I3x3y_S = 0.0E0;
  Double I_NAI_I3x2yz_S = 0.0E0;
  Double I_NAI_I3xy2z_S = 0.0E0;
  Double I_NAI_I3x3z_S = 0.0E0;
  Double I_NAI_I2x4y_S = 0.0E0;
  Double I_NAI_I2x3yz_S = 0.0E0;
  Double I_NAI_I2x2y2z_S = 0.0E0;
  Double I_NAI_I2xy3z_S = 0.0E0;
  Double I_NAI_I2x4z_S = 0.0E0;
  Double I_NAI_Ix5y_S = 0.0E0;
  Double I_NAI_Ix4yz_S = 0.0E0;
  Double I_NAI_Ix3y2z_S = 0.0E0;
  Double I_NAI_Ix2y3z_S = 0.0E0;
  Double I_NAI_Ixy4z_S = 0.0E0;
  Double I_NAI_Ix5z_S = 0.0E0;
  Double I_NAI_I6y_S = 0.0E0;
  Double I_NAI_I5yz_S = 0.0E0;
  Double I_NAI_I4y2z_S = 0.0E0;
  Double I_NAI_I3y3z_S = 0.0E0;
  Double I_NAI_I2y4z_S = 0.0E0;
  Double I_NAI_Iy5z_S = 0.0E0;
  Double I_NAI_I6z_S = 0.0E0;
  Double I_NAI_H5x_S = 0.0E0;
  Double I_NAI_H4xy_S = 0.0E0;
  Double I_NAI_H4xz_S = 0.0E0;
  Double I_NAI_H3x2y_S = 0.0E0;
  Double I_NAI_H3xyz_S = 0.0E0;
  Double I_NAI_H3x2z_S = 0.0E0;
  Double I_NAI_H2x3y_S = 0.0E0;
  Double I_NAI_H2x2yz_S = 0.0E0;
  Double I_NAI_H2xy2z_S = 0.0E0;
  Double I_NAI_H2x3z_S = 0.0E0;
  Double I_NAI_Hx4y_S = 0.0E0;
  Double I_NAI_Hx3yz_S = 0.0E0;
  Double I_NAI_Hx2y2z_S = 0.0E0;
  Double I_NAI_Hxy3z_S = 0.0E0;
  Double I_NAI_Hx4z_S = 0.0E0;
  Double I_NAI_H5y_S = 0.0E0;
  Double I_NAI_H4yz_S = 0.0E0;
  Double I_NAI_H3y2z_S = 0.0E0;
  Double I_NAI_H2y3z_S = 0.0E0;
  Double I_NAI_Hy4z_S = 0.0E0;
  Double I_NAI_H5z_S = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double rho   = 1.0E0/onedz;
    Double sqrho = sqrt(rho);
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
      Double PNX   = PX - N[iAtom*3  ];
      Double PNY   = PY - N[iAtom*3+1];
      Double PNZ   = PZ - N[iAtom*3+2];
      Double PN2   = PNX*PNX+PNY*PNY+PNZ*PNZ;
      Double charge= Z[iAtom];
      Double u     = rho*PN2;
      Double squ   = sqrt(u);
      Double prefactor = -ic2*charge*fbra;

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

      Double I_NAI_S_S_vrr  = 0.0E0;
      Double I_NAI_S_S_M1_vrr  = 0.0E0;
      Double I_NAI_S_S_M2_vrr  = 0.0E0;
      Double I_NAI_S_S_M3_vrr  = 0.0E0;
      Double I_NAI_S_S_M4_vrr  = 0.0E0;
      Double I_NAI_S_S_M5_vrr  = 0.0E0;
      Double I_NAI_S_S_M6_vrr  = 0.0E0;
      Double I_NAI_S_S_M7_vrr  = 0.0E0;
      Double I_NAI_S_S_M8_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
      double I_NAI_S_S_M5_vrr_d  = 0.0E0;
      double I_NAI_S_S_M6_vrr_d  = 0.0E0;
      double I_NAI_S_S_M7_vrr_d  = 0.0E0;
      double I_NAI_S_S_M8_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER51;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER49*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER47*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M8_vrr;
        I_NAI_S_S_M8_vrr = ONEOVER17*I_NAI_S_S_M8_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M8_vrr  = f*I_NAI_S_S_M8_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_NAI_S_S_M7_vrr  = ONEOVER15*(u2*I_NAI_S_S_M8_vrr+f);
        I_NAI_S_S_M6_vrr  = ONEOVER13*(u2*I_NAI_S_S_M7_vrr+f);
        I_NAI_S_S_M5_vrr  = ONEOVER11*(u2*I_NAI_S_S_M6_vrr+f);
        I_NAI_S_S_M4_vrr  = ONEOVER9*(u2*I_NAI_S_S_M5_vrr+f);
        I_NAI_S_S_M3_vrr  = ONEOVER7*(u2*I_NAI_S_S_M4_vrr+f);
        I_NAI_S_S_M2_vrr  = ONEOVER5*(u2*I_NAI_S_S_M3_vrr+f);
        I_NAI_S_S_M1_vrr  = ONEOVER3*(u2*I_NAI_S_S_M2_vrr+f);
        I_NAI_S_S_vrr  = ONEOVER1*(u2*I_NAI_S_S_M1_vrr+f);

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
          I_NAI_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_NAI_S_S_M1_vrr_d = oneO2u*(1.0E0*I_NAI_S_S_vrr_d-f);
        I_NAI_S_S_M2_vrr_d = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr_d-f);
        I_NAI_S_S_M3_vrr_d = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr_d-f);
        I_NAI_S_S_M4_vrr_d = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr_d-f);
        I_NAI_S_S_M5_vrr_d = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr_d-f);
        I_NAI_S_S_M6_vrr_d = oneO2u*(11.0E0*I_NAI_S_S_M5_vrr_d-f);
        I_NAI_S_S_M7_vrr_d = oneO2u*(13.0E0*I_NAI_S_S_M6_vrr_d-f);
        I_NAI_S_S_M8_vrr_d = oneO2u*(15.0E0*I_NAI_S_S_M7_vrr_d-f);

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);
        I_NAI_S_S_M5_vrr = static_cast<Double>(I_NAI_S_S_M5_vrr_d);
        I_NAI_S_S_M6_vrr = static_cast<Double>(I_NAI_S_S_M6_vrr_d);
        I_NAI_S_S_M7_vrr = static_cast<Double>(I_NAI_S_S_M7_vrr_d);
        I_NAI_S_S_M8_vrr = static_cast<Double>(I_NAI_S_S_M8_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_NAI_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M1_vrr = oneO2u*(1.0E0*I_NAI_S_S_vrr-f);
        I_NAI_S_S_M2_vrr = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr-f);
        I_NAI_S_S_M3_vrr = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr-f);
        I_NAI_S_S_M4_vrr = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr-f);
        I_NAI_S_S_M5_vrr = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr-f);
        I_NAI_S_S_M6_vrr = oneO2u*(11.0E0*I_NAI_S_S_M5_vrr-f);
        I_NAI_S_S_M7_vrr = oneO2u*(13.0E0*I_NAI_S_S_M6_vrr-f);
        I_NAI_S_S_M8_vrr = oneO2u*(15.0E0*I_NAI_S_S_M7_vrr-f);

#endif

      }


      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M7
       * RHS shell quartet name: SQ_NAI_S_S_M8
       ************************************************************/
      Double I_NAI_Px_S_M7_vrr = PAX*I_NAI_S_S_M7_vrr-PNX*I_NAI_S_S_M8_vrr;
      Double I_NAI_Py_S_M7_vrr = PAY*I_NAI_S_S_M7_vrr-PNY*I_NAI_S_S_M8_vrr;
      Double I_NAI_Pz_S_M7_vrr = PAZ*I_NAI_S_S_M7_vrr-PNZ*I_NAI_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M6
       * RHS shell quartet name: SQ_NAI_S_S_M7
       ************************************************************/
      Double I_NAI_Px_S_M6_vrr = PAX*I_NAI_S_S_M6_vrr-PNX*I_NAI_S_S_M7_vrr;
      Double I_NAI_Py_S_M6_vrr = PAY*I_NAI_S_S_M6_vrr-PNY*I_NAI_S_S_M7_vrr;
      Double I_NAI_Pz_S_M6_vrr = PAZ*I_NAI_S_S_M6_vrr-PNZ*I_NAI_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M6
       * RHS shell quartet name: SQ_NAI_P_S_M7
       * RHS shell quartet name: SQ_NAI_S_S_M6
       * RHS shell quartet name: SQ_NAI_S_S_M7
       ************************************************************/
      Double I_NAI_D2x_S_M6_vrr = PAX*I_NAI_Px_S_M6_vrr-PNX*I_NAI_Px_S_M7_vrr+oned2z*I_NAI_S_S_M6_vrr-oned2z*I_NAI_S_S_M7_vrr;
      Double I_NAI_D2y_S_M6_vrr = PAY*I_NAI_Py_S_M6_vrr-PNY*I_NAI_Py_S_M7_vrr+oned2z*I_NAI_S_S_M6_vrr-oned2z*I_NAI_S_S_M7_vrr;
      Double I_NAI_D2z_S_M6_vrr = PAZ*I_NAI_Pz_S_M6_vrr-PNZ*I_NAI_Pz_S_M7_vrr+oned2z*I_NAI_S_S_M6_vrr-oned2z*I_NAI_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M6
       ************************************************************/
      Double I_NAI_Px_S_M5_vrr = PAX*I_NAI_S_S_M5_vrr-PNX*I_NAI_S_S_M6_vrr;
      Double I_NAI_Py_S_M5_vrr = PAY*I_NAI_S_S_M5_vrr-PNY*I_NAI_S_S_M6_vrr;
      Double I_NAI_Pz_S_M5_vrr = PAZ*I_NAI_S_S_M5_vrr-PNZ*I_NAI_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M5
       * RHS shell quartet name: SQ_NAI_P_S_M6
       * RHS shell quartet name: SQ_NAI_S_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M6
       ************************************************************/
      Double I_NAI_D2x_S_M5_vrr = PAX*I_NAI_Px_S_M5_vrr-PNX*I_NAI_Px_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;
      Double I_NAI_D2y_S_M5_vrr = PAY*I_NAI_Py_S_M5_vrr-PNY*I_NAI_Py_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;
      Double I_NAI_D2z_S_M5_vrr = PAZ*I_NAI_Pz_S_M5_vrr-PNZ*I_NAI_Pz_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M5
       * RHS shell quartet name: SQ_NAI_D_S_M6
       * RHS shell quartet name: SQ_NAI_P_S_M5
       * RHS shell quartet name: SQ_NAI_P_S_M6
       ************************************************************/
      Double I_NAI_F3x_S_M5_vrr = PAX*I_NAI_D2x_S_M5_vrr-PNX*I_NAI_D2x_S_M6_vrr+2*oned2z*I_NAI_Px_S_M5_vrr-2*oned2z*I_NAI_Px_S_M6_vrr;
      Double I_NAI_F3y_S_M5_vrr = PAY*I_NAI_D2y_S_M5_vrr-PNY*I_NAI_D2y_S_M6_vrr+2*oned2z*I_NAI_Py_S_M5_vrr-2*oned2z*I_NAI_Py_S_M6_vrr;
      Double I_NAI_F3z_S_M5_vrr = PAZ*I_NAI_D2z_S_M5_vrr-PNZ*I_NAI_D2z_S_M6_vrr+2*oned2z*I_NAI_Pz_S_M5_vrr-2*oned2z*I_NAI_Pz_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M5
       ************************************************************/
      Double I_NAI_Px_S_M4_vrr = PAX*I_NAI_S_S_M4_vrr-PNX*I_NAI_S_S_M5_vrr;
      Double I_NAI_Py_S_M4_vrr = PAY*I_NAI_S_S_M4_vrr-PNY*I_NAI_S_S_M5_vrr;
      Double I_NAI_Pz_S_M4_vrr = PAZ*I_NAI_S_S_M4_vrr-PNZ*I_NAI_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M5
       ************************************************************/
      Double I_NAI_D2x_S_M4_vrr = PAX*I_NAI_Px_S_M4_vrr-PNX*I_NAI_Px_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;
      Double I_NAI_D2y_S_M4_vrr = PAY*I_NAI_Py_S_M4_vrr-PNY*I_NAI_Py_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;
      Double I_NAI_D2z_S_M4_vrr = PAZ*I_NAI_Pz_S_M4_vrr-PNZ*I_NAI_Pz_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M4
       * RHS shell quartet name: SQ_NAI_D_S_M5
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M5
       ************************************************************/
      Double I_NAI_F3x_S_M4_vrr = PAX*I_NAI_D2x_S_M4_vrr-PNX*I_NAI_D2x_S_M5_vrr+2*oned2z*I_NAI_Px_S_M4_vrr-2*oned2z*I_NAI_Px_S_M5_vrr;
      Double I_NAI_F3y_S_M4_vrr = PAY*I_NAI_D2y_S_M4_vrr-PNY*I_NAI_D2y_S_M5_vrr+2*oned2z*I_NAI_Py_S_M4_vrr-2*oned2z*I_NAI_Py_S_M5_vrr;
      Double I_NAI_F3z_S_M4_vrr = PAZ*I_NAI_D2z_S_M4_vrr-PNZ*I_NAI_D2z_S_M5_vrr+2*oned2z*I_NAI_Pz_S_M4_vrr-2*oned2z*I_NAI_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M4
       * RHS shell quartet name: SQ_NAI_F_S_M5
       * RHS shell quartet name: SQ_NAI_D_S_M4
       * RHS shell quartet name: SQ_NAI_D_S_M5
       ************************************************************/
      Double I_NAI_G4x_S_M4_vrr = PAX*I_NAI_F3x_S_M4_vrr-PNX*I_NAI_F3x_S_M5_vrr+3*oned2z*I_NAI_D2x_S_M4_vrr-3*oned2z*I_NAI_D2x_S_M5_vrr;
      Double I_NAI_G3xy_S_M4_vrr = PAY*I_NAI_F3x_S_M4_vrr-PNY*I_NAI_F3x_S_M5_vrr;
      Double I_NAI_G3xz_S_M4_vrr = PAZ*I_NAI_F3x_S_M4_vrr-PNZ*I_NAI_F3x_S_M5_vrr;
      Double I_NAI_G4y_S_M4_vrr = PAY*I_NAI_F3y_S_M4_vrr-PNY*I_NAI_F3y_S_M5_vrr+3*oned2z*I_NAI_D2y_S_M4_vrr-3*oned2z*I_NAI_D2y_S_M5_vrr;
      Double I_NAI_G3yz_S_M4_vrr = PAZ*I_NAI_F3y_S_M4_vrr-PNZ*I_NAI_F3y_S_M5_vrr;
      Double I_NAI_G4z_S_M4_vrr = PAZ*I_NAI_F3z_S_M4_vrr-PNZ*I_NAI_F3z_S_M5_vrr+3*oned2z*I_NAI_D2z_S_M4_vrr-3*oned2z*I_NAI_D2z_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_Px_S_M3_vrr = PAX*I_NAI_S_S_M3_vrr-PNX*I_NAI_S_S_M4_vrr;
      Double I_NAI_Py_S_M3_vrr = PAY*I_NAI_S_S_M3_vrr-PNY*I_NAI_S_S_M4_vrr;
      Double I_NAI_Pz_S_M3_vrr = PAZ*I_NAI_S_S_M3_vrr-PNZ*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_D2x_S_M3_vrr = PAX*I_NAI_Px_S_M3_vrr-PNX*I_NAI_Px_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2y_S_M3_vrr = PAY*I_NAI_Py_S_M3_vrr-PNY*I_NAI_Py_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2z_S_M3_vrr = PAZ*I_NAI_Pz_S_M3_vrr-PNZ*I_NAI_Pz_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M4
       ************************************************************/
      Double I_NAI_F3x_S_M3_vrr = PAX*I_NAI_D2x_S_M3_vrr-PNX*I_NAI_D2x_S_M4_vrr+2*oned2z*I_NAI_Px_S_M3_vrr-2*oned2z*I_NAI_Px_S_M4_vrr;
      Double I_NAI_F2xy_S_M3_vrr = PAY*I_NAI_D2x_S_M3_vrr-PNY*I_NAI_D2x_S_M4_vrr;
      Double I_NAI_F3y_S_M3_vrr = PAY*I_NAI_D2y_S_M3_vrr-PNY*I_NAI_D2y_S_M4_vrr+2*oned2z*I_NAI_Py_S_M3_vrr-2*oned2z*I_NAI_Py_S_M4_vrr;
      Double I_NAI_F3z_S_M3_vrr = PAZ*I_NAI_D2z_S_M3_vrr-PNZ*I_NAI_D2z_S_M4_vrr+2*oned2z*I_NAI_Pz_S_M3_vrr-2*oned2z*I_NAI_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M3
       * RHS shell quartet name: SQ_NAI_F_S_M4
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M4
       ************************************************************/
      Double I_NAI_G4x_S_M3_vrr = PAX*I_NAI_F3x_S_M3_vrr-PNX*I_NAI_F3x_S_M4_vrr+3*oned2z*I_NAI_D2x_S_M3_vrr-3*oned2z*I_NAI_D2x_S_M4_vrr;
      Double I_NAI_G3xy_S_M3_vrr = PAY*I_NAI_F3x_S_M3_vrr-PNY*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_G3xz_S_M3_vrr = PAZ*I_NAI_F3x_S_M3_vrr-PNZ*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_Gx3y_S_M3_vrr = PAX*I_NAI_F3y_S_M3_vrr-PNX*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_Gx3z_S_M3_vrr = PAX*I_NAI_F3z_S_M3_vrr-PNX*I_NAI_F3z_S_M4_vrr;
      Double I_NAI_G4y_S_M3_vrr = PAY*I_NAI_F3y_S_M3_vrr-PNY*I_NAI_F3y_S_M4_vrr+3*oned2z*I_NAI_D2y_S_M3_vrr-3*oned2z*I_NAI_D2y_S_M4_vrr;
      Double I_NAI_G3yz_S_M3_vrr = PAZ*I_NAI_F3y_S_M3_vrr-PNZ*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_G4z_S_M3_vrr = PAZ*I_NAI_F3z_S_M3_vrr-PNZ*I_NAI_F3z_S_M4_vrr+3*oned2z*I_NAI_D2z_S_M3_vrr-3*oned2z*I_NAI_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M3
       * RHS shell quartet name: SQ_NAI_G_S_M4
       * RHS shell quartet name: SQ_NAI_F_S_M3
       * RHS shell quartet name: SQ_NAI_F_S_M4
       ************************************************************/
      Double I_NAI_H5x_S_M3_vrr = PAX*I_NAI_G4x_S_M3_vrr-PNX*I_NAI_G4x_S_M4_vrr+4*oned2z*I_NAI_F3x_S_M3_vrr-4*oned2z*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_H4xy_S_M3_vrr = PAY*I_NAI_G4x_S_M3_vrr-PNY*I_NAI_G4x_S_M4_vrr;
      Double I_NAI_H4xz_S_M3_vrr = PAZ*I_NAI_G4x_S_M3_vrr-PNZ*I_NAI_G4x_S_M4_vrr;
      Double I_NAI_H3x2y_S_M3_vrr = PAY*I_NAI_G3xy_S_M3_vrr-PNY*I_NAI_G3xy_S_M4_vrr+oned2z*I_NAI_F3x_S_M3_vrr-oned2z*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_H3x2z_S_M3_vrr = PAZ*I_NAI_G3xz_S_M3_vrr-PNZ*I_NAI_G3xz_S_M4_vrr+oned2z*I_NAI_F3x_S_M3_vrr-oned2z*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_Hx4y_S_M3_vrr = PAX*I_NAI_G4y_S_M3_vrr-PNX*I_NAI_G4y_S_M4_vrr;
      Double I_NAI_Hx4z_S_M3_vrr = PAX*I_NAI_G4z_S_M3_vrr-PNX*I_NAI_G4z_S_M4_vrr;
      Double I_NAI_H5y_S_M3_vrr = PAY*I_NAI_G4y_S_M3_vrr-PNY*I_NAI_G4y_S_M4_vrr+4*oned2z*I_NAI_F3y_S_M3_vrr-4*oned2z*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_H4yz_S_M3_vrr = PAZ*I_NAI_G4y_S_M3_vrr-PNZ*I_NAI_G4y_S_M4_vrr;
      Double I_NAI_H3y2z_S_M3_vrr = PAZ*I_NAI_G3yz_S_M3_vrr-PNZ*I_NAI_G3yz_S_M4_vrr+oned2z*I_NAI_F3y_S_M3_vrr-oned2z*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_Hy4z_S_M3_vrr = PAY*I_NAI_G4z_S_M3_vrr-PNY*I_NAI_G4z_S_M4_vrr;
      Double I_NAI_H5z_S_M3_vrr = PAZ*I_NAI_G4z_S_M3_vrr-PNZ*I_NAI_G4z_S_M4_vrr+4*oned2z*I_NAI_F3z_S_M3_vrr-4*oned2z*I_NAI_F3z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_Px_S_M2_vrr = PAX*I_NAI_S_S_M2_vrr-PNX*I_NAI_S_S_M3_vrr;
      Double I_NAI_Py_S_M2_vrr = PAY*I_NAI_S_S_M2_vrr-PNY*I_NAI_S_S_M3_vrr;
      Double I_NAI_Pz_S_M2_vrr = PAZ*I_NAI_S_S_M2_vrr-PNZ*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_D2x_S_M2_vrr = PAX*I_NAI_Px_S_M2_vrr-PNX*I_NAI_Px_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2y_S_M2_vrr = PAY*I_NAI_Py_S_M2_vrr-PNY*I_NAI_Py_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2z_S_M2_vrr = PAZ*I_NAI_Pz_S_M2_vrr-PNZ*I_NAI_Pz_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       ************************************************************/
      Double I_NAI_F3x_S_M2_vrr = PAX*I_NAI_D2x_S_M2_vrr-PNX*I_NAI_D2x_S_M3_vrr+2*oned2z*I_NAI_Px_S_M2_vrr-2*oned2z*I_NAI_Px_S_M3_vrr;
      Double I_NAI_F2xy_S_M2_vrr = PAY*I_NAI_D2x_S_M2_vrr-PNY*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F2xz_S_M2_vrr = PAZ*I_NAI_D2x_S_M2_vrr-PNZ*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F3y_S_M2_vrr = PAY*I_NAI_D2y_S_M2_vrr-PNY*I_NAI_D2y_S_M3_vrr+2*oned2z*I_NAI_Py_S_M2_vrr-2*oned2z*I_NAI_Py_S_M3_vrr;
      Double I_NAI_F2yz_S_M2_vrr = PAZ*I_NAI_D2y_S_M2_vrr-PNZ*I_NAI_D2y_S_M3_vrr;
      Double I_NAI_F3z_S_M2_vrr = PAZ*I_NAI_D2z_S_M2_vrr-PNZ*I_NAI_D2z_S_M3_vrr+2*oned2z*I_NAI_Pz_S_M2_vrr-2*oned2z*I_NAI_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M3
       ************************************************************/
      Double I_NAI_G4x_S_M2_vrr = PAX*I_NAI_F3x_S_M2_vrr-PNX*I_NAI_F3x_S_M3_vrr+3*oned2z*I_NAI_D2x_S_M2_vrr-3*oned2z*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_G3xy_S_M2_vrr = PAY*I_NAI_F3x_S_M2_vrr-PNY*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_G3xz_S_M2_vrr = PAZ*I_NAI_F3x_S_M2_vrr-PNZ*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_G2x2y_S_M2_vrr = PAY*I_NAI_F2xy_S_M2_vrr-PNY*I_NAI_F2xy_S_M3_vrr+oned2z*I_NAI_D2x_S_M2_vrr-oned2z*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_Gx3y_S_M2_vrr = PAX*I_NAI_F3y_S_M2_vrr-PNX*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Gx3z_S_M2_vrr = PAX*I_NAI_F3z_S_M2_vrr-PNX*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_G4y_S_M2_vrr = PAY*I_NAI_F3y_S_M2_vrr-PNY*I_NAI_F3y_S_M3_vrr+3*oned2z*I_NAI_D2y_S_M2_vrr-3*oned2z*I_NAI_D2y_S_M3_vrr;
      Double I_NAI_G3yz_S_M2_vrr = PAZ*I_NAI_F3y_S_M2_vrr-PNZ*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Gy3z_S_M2_vrr = PAY*I_NAI_F3z_S_M2_vrr-PNY*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_G4z_S_M2_vrr = PAZ*I_NAI_F3z_S_M2_vrr-PNZ*I_NAI_F3z_S_M3_vrr+3*oned2z*I_NAI_D2z_S_M2_vrr-3*oned2z*I_NAI_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M2
       * RHS shell quartet name: SQ_NAI_G_S_M3
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M3
       ************************************************************/
      Double I_NAI_H5x_S_M2_vrr = PAX*I_NAI_G4x_S_M2_vrr-PNX*I_NAI_G4x_S_M3_vrr+4*oned2z*I_NAI_F3x_S_M2_vrr-4*oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H4xy_S_M2_vrr = PAY*I_NAI_G4x_S_M2_vrr-PNY*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_H4xz_S_M2_vrr = PAZ*I_NAI_G4x_S_M2_vrr-PNZ*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_H3x2y_S_M2_vrr = PAY*I_NAI_G3xy_S_M2_vrr-PNY*I_NAI_G3xy_S_M3_vrr+oned2z*I_NAI_F3x_S_M2_vrr-oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H3x2z_S_M2_vrr = PAZ*I_NAI_G3xz_S_M2_vrr-PNZ*I_NAI_G3xz_S_M3_vrr+oned2z*I_NAI_F3x_S_M2_vrr-oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H2x3y_S_M2_vrr = PAX*I_NAI_Gx3y_S_M2_vrr-PNX*I_NAI_Gx3y_S_M3_vrr+oned2z*I_NAI_F3y_S_M2_vrr-oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_H2x3z_S_M2_vrr = PAX*I_NAI_Gx3z_S_M2_vrr-PNX*I_NAI_Gx3z_S_M3_vrr+oned2z*I_NAI_F3z_S_M2_vrr-oned2z*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_Hx4y_S_M2_vrr = PAX*I_NAI_G4y_S_M2_vrr-PNX*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_Hx4z_S_M2_vrr = PAX*I_NAI_G4z_S_M2_vrr-PNX*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_H5y_S_M2_vrr = PAY*I_NAI_G4y_S_M2_vrr-PNY*I_NAI_G4y_S_M3_vrr+4*oned2z*I_NAI_F3y_S_M2_vrr-4*oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_H4yz_S_M2_vrr = PAZ*I_NAI_G4y_S_M2_vrr-PNZ*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_H3y2z_S_M2_vrr = PAZ*I_NAI_G3yz_S_M2_vrr-PNZ*I_NAI_G3yz_S_M3_vrr+oned2z*I_NAI_F3y_S_M2_vrr-oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Hy4z_S_M2_vrr = PAY*I_NAI_G4z_S_M2_vrr-PNY*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_H5z_S_M2_vrr = PAZ*I_NAI_G4z_S_M2_vrr-PNZ*I_NAI_G4z_S_M3_vrr+4*oned2z*I_NAI_F3z_S_M2_vrr-4*oned2z*I_NAI_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 10 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S_M2
       * RHS shell quartet name: SQ_NAI_H_S_M3
       * RHS shell quartet name: SQ_NAI_G_S_M2
       * RHS shell quartet name: SQ_NAI_G_S_M3
       ************************************************************/
      Double I_NAI_I6x_S_M2_vrr = PAX*I_NAI_H5x_S_M2_vrr-PNX*I_NAI_H5x_S_M3_vrr+5*oned2z*I_NAI_G4x_S_M2_vrr-5*oned2z*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_I5xy_S_M2_vrr = PAY*I_NAI_H5x_S_M2_vrr-PNY*I_NAI_H5x_S_M3_vrr;
      Double I_NAI_I5xz_S_M2_vrr = PAZ*I_NAI_H5x_S_M2_vrr-PNZ*I_NAI_H5x_S_M3_vrr;
      Double I_NAI_I4x2y_S_M2_vrr = PAY*I_NAI_H4xy_S_M2_vrr-PNY*I_NAI_H4xy_S_M3_vrr+oned2z*I_NAI_G4x_S_M2_vrr-oned2z*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_I4x2z_S_M2_vrr = PAZ*I_NAI_H4xz_S_M2_vrr-PNZ*I_NAI_H4xz_S_M3_vrr+oned2z*I_NAI_G4x_S_M2_vrr-oned2z*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_I3x3y_S_M2_vrr = PAY*I_NAI_H3x2y_S_M2_vrr-PNY*I_NAI_H3x2y_S_M3_vrr+2*oned2z*I_NAI_G3xy_S_M2_vrr-2*oned2z*I_NAI_G3xy_S_M3_vrr;
      Double I_NAI_I3x3z_S_M2_vrr = PAZ*I_NAI_H3x2z_S_M2_vrr-PNZ*I_NAI_H3x2z_S_M3_vrr+2*oned2z*I_NAI_G3xz_S_M2_vrr-2*oned2z*I_NAI_G3xz_S_M3_vrr;
      Double I_NAI_I2x4y_S_M2_vrr = PAX*I_NAI_Hx4y_S_M2_vrr-PNX*I_NAI_Hx4y_S_M3_vrr+oned2z*I_NAI_G4y_S_M2_vrr-oned2z*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_I2x4z_S_M2_vrr = PAX*I_NAI_Hx4z_S_M2_vrr-PNX*I_NAI_Hx4z_S_M3_vrr+oned2z*I_NAI_G4z_S_M2_vrr-oned2z*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_Ix5y_S_M2_vrr = PAX*I_NAI_H5y_S_M2_vrr-PNX*I_NAI_H5y_S_M3_vrr;
      Double I_NAI_Ix5z_S_M2_vrr = PAX*I_NAI_H5z_S_M2_vrr-PNX*I_NAI_H5z_S_M3_vrr;
      Double I_NAI_I6y_S_M2_vrr = PAY*I_NAI_H5y_S_M2_vrr-PNY*I_NAI_H5y_S_M3_vrr+5*oned2z*I_NAI_G4y_S_M2_vrr-5*oned2z*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_I5yz_S_M2_vrr = PAZ*I_NAI_H5y_S_M2_vrr-PNZ*I_NAI_H5y_S_M3_vrr;
      Double I_NAI_I4y2z_S_M2_vrr = PAZ*I_NAI_H4yz_S_M2_vrr-PNZ*I_NAI_H4yz_S_M3_vrr+oned2z*I_NAI_G4y_S_M2_vrr-oned2z*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_I3y3z_S_M2_vrr = PAZ*I_NAI_H3y2z_S_M2_vrr-PNZ*I_NAI_H3y2z_S_M3_vrr+2*oned2z*I_NAI_G3yz_S_M2_vrr-2*oned2z*I_NAI_G3yz_S_M3_vrr;
      Double I_NAI_I2y4z_S_M2_vrr = PAY*I_NAI_Hy4z_S_M2_vrr-PNY*I_NAI_Hy4z_S_M3_vrr+oned2z*I_NAI_G4z_S_M2_vrr-oned2z*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_Iy5z_S_M2_vrr = PAY*I_NAI_H5z_S_M2_vrr-PNY*I_NAI_H5z_S_M3_vrr;
      Double I_NAI_I6z_S_M2_vrr = PAZ*I_NAI_H5z_S_M2_vrr-PNZ*I_NAI_H5z_S_M3_vrr+5*oned2z*I_NAI_G4z_S_M2_vrr-5*oned2z*I_NAI_G4z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_Px_S_M1_vrr = PAX*I_NAI_S_S_M1_vrr-PNX*I_NAI_S_S_M2_vrr;
      Double I_NAI_Py_S_M1_vrr = PAY*I_NAI_S_S_M1_vrr-PNY*I_NAI_S_S_M2_vrr;
      Double I_NAI_Pz_S_M1_vrr = PAZ*I_NAI_S_S_M1_vrr-PNZ*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_D2x_S_M1_vrr = PAX*I_NAI_Px_S_M1_vrr-PNX*I_NAI_Px_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_D2y_S_M1_vrr = PAY*I_NAI_Py_S_M1_vrr-PNY*I_NAI_Py_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_D2z_S_M1_vrr = PAZ*I_NAI_Pz_S_M1_vrr-PNZ*I_NAI_Pz_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       ************************************************************/
      Double I_NAI_F3x_S_M1_vrr = PAX*I_NAI_D2x_S_M1_vrr-PNX*I_NAI_D2x_S_M2_vrr+2*oned2z*I_NAI_Px_S_M1_vrr-2*oned2z*I_NAI_Px_S_M2_vrr;
      Double I_NAI_F2xy_S_M1_vrr = PAY*I_NAI_D2x_S_M1_vrr-PNY*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_F2xz_S_M1_vrr = PAZ*I_NAI_D2x_S_M1_vrr-PNZ*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_F3y_S_M1_vrr = PAY*I_NAI_D2y_S_M1_vrr-PNY*I_NAI_D2y_S_M2_vrr+2*oned2z*I_NAI_Py_S_M1_vrr-2*oned2z*I_NAI_Py_S_M2_vrr;
      Double I_NAI_F2yz_S_M1_vrr = PAZ*I_NAI_D2y_S_M1_vrr-PNZ*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_F3z_S_M1_vrr = PAZ*I_NAI_D2z_S_M1_vrr-PNZ*I_NAI_D2z_S_M2_vrr+2*oned2z*I_NAI_Pz_S_M1_vrr-2*oned2z*I_NAI_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       ************************************************************/
      Double I_NAI_G4x_S_M1_vrr = PAX*I_NAI_F3x_S_M1_vrr-PNX*I_NAI_F3x_S_M2_vrr+3*oned2z*I_NAI_D2x_S_M1_vrr-3*oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G3xy_S_M1_vrr = PAY*I_NAI_F3x_S_M1_vrr-PNY*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G3xz_S_M1_vrr = PAZ*I_NAI_F3x_S_M1_vrr-PNZ*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G2x2y_S_M1_vrr = PAY*I_NAI_F2xy_S_M1_vrr-PNY*I_NAI_F2xy_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G2x2z_S_M1_vrr = PAZ*I_NAI_F2xz_S_M1_vrr-PNZ*I_NAI_F2xz_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Gx3y_S_M1_vrr = PAX*I_NAI_F3y_S_M1_vrr-PNX*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_Gx3z_S_M1_vrr = PAX*I_NAI_F3z_S_M1_vrr-PNX*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4y_S_M1_vrr = PAY*I_NAI_F3y_S_M1_vrr-PNY*I_NAI_F3y_S_M2_vrr+3*oned2z*I_NAI_D2y_S_M1_vrr-3*oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_G3yz_S_M1_vrr = PAZ*I_NAI_F3y_S_M1_vrr-PNZ*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_G2y2z_S_M1_vrr = PAZ*I_NAI_F2yz_S_M1_vrr-PNZ*I_NAI_F2yz_S_M2_vrr+oned2z*I_NAI_D2y_S_M1_vrr-oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Gy3z_S_M1_vrr = PAY*I_NAI_F3z_S_M1_vrr-PNY*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4z_S_M1_vrr = PAZ*I_NAI_F3z_S_M1_vrr-PNZ*I_NAI_F3z_S_M2_vrr+3*oned2z*I_NAI_D2z_S_M1_vrr-3*oned2z*I_NAI_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_G_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_F_S_M2
       ************************************************************/
      Double I_NAI_H5x_S_M1_vrr = PAX*I_NAI_G4x_S_M1_vrr-PNX*I_NAI_G4x_S_M2_vrr+4*oned2z*I_NAI_F3x_S_M1_vrr-4*oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H4xy_S_M1_vrr = PAY*I_NAI_G4x_S_M1_vrr-PNY*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_H4xz_S_M1_vrr = PAZ*I_NAI_G4x_S_M1_vrr-PNZ*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_H3x2y_S_M1_vrr = PAY*I_NAI_G3xy_S_M1_vrr-PNY*I_NAI_G3xy_S_M2_vrr+oned2z*I_NAI_F3x_S_M1_vrr-oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H3x2z_S_M1_vrr = PAZ*I_NAI_G3xz_S_M1_vrr-PNZ*I_NAI_G3xz_S_M2_vrr+oned2z*I_NAI_F3x_S_M1_vrr-oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H2x3y_S_M1_vrr = PAX*I_NAI_Gx3y_S_M1_vrr-PNX*I_NAI_Gx3y_S_M2_vrr+oned2z*I_NAI_F3y_S_M1_vrr-oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H2x2yz_S_M1_vrr = PAZ*I_NAI_G2x2y_S_M1_vrr-PNZ*I_NAI_G2x2y_S_M2_vrr;
      Double I_NAI_H2x3z_S_M1_vrr = PAX*I_NAI_Gx3z_S_M1_vrr-PNX*I_NAI_Gx3z_S_M2_vrr+oned2z*I_NAI_F3z_S_M1_vrr-oned2z*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_Hx4y_S_M1_vrr = PAX*I_NAI_G4y_S_M1_vrr-PNX*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_Hx4z_S_M1_vrr = PAX*I_NAI_G4z_S_M1_vrr-PNX*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_H5y_S_M1_vrr = PAY*I_NAI_G4y_S_M1_vrr-PNY*I_NAI_G4y_S_M2_vrr+4*oned2z*I_NAI_F3y_S_M1_vrr-4*oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H4yz_S_M1_vrr = PAZ*I_NAI_G4y_S_M1_vrr-PNZ*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_H3y2z_S_M1_vrr = PAZ*I_NAI_G3yz_S_M1_vrr-PNZ*I_NAI_G3yz_S_M2_vrr+oned2z*I_NAI_F3y_S_M1_vrr-oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H2y3z_S_M1_vrr = PAY*I_NAI_Gy3z_S_M1_vrr-PNY*I_NAI_Gy3z_S_M2_vrr+oned2z*I_NAI_F3z_S_M1_vrr-oned2z*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_Hy4z_S_M1_vrr = PAY*I_NAI_G4z_S_M1_vrr-PNY*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_H5z_S_M1_vrr = PAZ*I_NAI_G4z_S_M1_vrr-PNZ*I_NAI_G4z_S_M2_vrr+4*oned2z*I_NAI_F3z_S_M1_vrr-4*oned2z*I_NAI_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S_M1
       * RHS shell quartet name: SQ_NAI_H_S_M2
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_G_S_M2
       ************************************************************/
      Double I_NAI_I6x_S_M1_vrr = PAX*I_NAI_H5x_S_M1_vrr-PNX*I_NAI_H5x_S_M2_vrr+5*oned2z*I_NAI_G4x_S_M1_vrr-5*oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I5xy_S_M1_vrr = PAY*I_NAI_H5x_S_M1_vrr-PNY*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_I5xz_S_M1_vrr = PAZ*I_NAI_H5x_S_M1_vrr-PNZ*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_I4x2y_S_M1_vrr = PAY*I_NAI_H4xy_S_M1_vrr-PNY*I_NAI_H4xy_S_M2_vrr+oned2z*I_NAI_G4x_S_M1_vrr-oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I4x2z_S_M1_vrr = PAZ*I_NAI_H4xz_S_M1_vrr-PNZ*I_NAI_H4xz_S_M2_vrr+oned2z*I_NAI_G4x_S_M1_vrr-oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I3x3y_S_M1_vrr = PAY*I_NAI_H3x2y_S_M1_vrr-PNY*I_NAI_H3x2y_S_M2_vrr+2*oned2z*I_NAI_G3xy_S_M1_vrr-2*oned2z*I_NAI_G3xy_S_M2_vrr;
      Double I_NAI_I3x2yz_S_M1_vrr = PAZ*I_NAI_H3x2y_S_M1_vrr-PNZ*I_NAI_H3x2y_S_M2_vrr;
      Double I_NAI_I3x3z_S_M1_vrr = PAZ*I_NAI_H3x2z_S_M1_vrr-PNZ*I_NAI_H3x2z_S_M2_vrr+2*oned2z*I_NAI_G3xz_S_M1_vrr-2*oned2z*I_NAI_G3xz_S_M2_vrr;
      Double I_NAI_I2x4y_S_M1_vrr = PAX*I_NAI_Hx4y_S_M1_vrr-PNX*I_NAI_Hx4y_S_M2_vrr+oned2z*I_NAI_G4y_S_M1_vrr-oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I2x3yz_S_M1_vrr = PAZ*I_NAI_H2x3y_S_M1_vrr-PNZ*I_NAI_H2x3y_S_M2_vrr;
      Double I_NAI_I2xy3z_S_M1_vrr = PAY*I_NAI_H2x3z_S_M1_vrr-PNY*I_NAI_H2x3z_S_M2_vrr;
      Double I_NAI_I2x4z_S_M1_vrr = PAX*I_NAI_Hx4z_S_M1_vrr-PNX*I_NAI_Hx4z_S_M2_vrr+oned2z*I_NAI_G4z_S_M1_vrr-oned2z*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_Ix5y_S_M1_vrr = PAX*I_NAI_H5y_S_M1_vrr-PNX*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_Ix5z_S_M1_vrr = PAX*I_NAI_H5z_S_M1_vrr-PNX*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_I6y_S_M1_vrr = PAY*I_NAI_H5y_S_M1_vrr-PNY*I_NAI_H5y_S_M2_vrr+5*oned2z*I_NAI_G4y_S_M1_vrr-5*oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I5yz_S_M1_vrr = PAZ*I_NAI_H5y_S_M1_vrr-PNZ*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_I4y2z_S_M1_vrr = PAZ*I_NAI_H4yz_S_M1_vrr-PNZ*I_NAI_H4yz_S_M2_vrr+oned2z*I_NAI_G4y_S_M1_vrr-oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I3y3z_S_M1_vrr = PAZ*I_NAI_H3y2z_S_M1_vrr-PNZ*I_NAI_H3y2z_S_M2_vrr+2*oned2z*I_NAI_G3yz_S_M1_vrr-2*oned2z*I_NAI_G3yz_S_M2_vrr;
      Double I_NAI_I2y4z_S_M1_vrr = PAY*I_NAI_Hy4z_S_M1_vrr-PNY*I_NAI_Hy4z_S_M2_vrr+oned2z*I_NAI_G4z_S_M1_vrr-oned2z*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_Iy5z_S_M1_vrr = PAY*I_NAI_H5z_S_M1_vrr-PNY*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_I6z_S_M1_vrr = PAZ*I_NAI_H5z_S_M1_vrr-PNZ*I_NAI_H5z_S_M2_vrr+5*oned2z*I_NAI_G4z_S_M1_vrr-5*oned2z*I_NAI_G4z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_I_S_M1
       * RHS shell quartet name: SQ_NAI_I_S_M2
       * RHS shell quartet name: SQ_NAI_H_S_M1
       * RHS shell quartet name: SQ_NAI_H_S_M2
       ************************************************************/
      Double I_NAI_K7x_S_M1_vrr = PAX*I_NAI_I6x_S_M1_vrr-PNX*I_NAI_I6x_S_M2_vrr+6*oned2z*I_NAI_H5x_S_M1_vrr-6*oned2z*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_K6xy_S_M1_vrr = PAY*I_NAI_I6x_S_M1_vrr-PNY*I_NAI_I6x_S_M2_vrr;
      Double I_NAI_K6xz_S_M1_vrr = PAZ*I_NAI_I6x_S_M1_vrr-PNZ*I_NAI_I6x_S_M2_vrr;
      Double I_NAI_K5x2y_S_M1_vrr = PAY*I_NAI_I5xy_S_M1_vrr-PNY*I_NAI_I5xy_S_M2_vrr+oned2z*I_NAI_H5x_S_M1_vrr-oned2z*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_K5x2z_S_M1_vrr = PAZ*I_NAI_I5xz_S_M1_vrr-PNZ*I_NAI_I5xz_S_M2_vrr+oned2z*I_NAI_H5x_S_M1_vrr-oned2z*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_K4x3y_S_M1_vrr = PAY*I_NAI_I4x2y_S_M1_vrr-PNY*I_NAI_I4x2y_S_M2_vrr+2*oned2z*I_NAI_H4xy_S_M1_vrr-2*oned2z*I_NAI_H4xy_S_M2_vrr;
      Double I_NAI_K4x2yz_S_M1_vrr = PAZ*I_NAI_I4x2y_S_M1_vrr-PNZ*I_NAI_I4x2y_S_M2_vrr;
      Double I_NAI_K4x3z_S_M1_vrr = PAZ*I_NAI_I4x2z_S_M1_vrr-PNZ*I_NAI_I4x2z_S_M2_vrr+2*oned2z*I_NAI_H4xz_S_M1_vrr-2*oned2z*I_NAI_H4xz_S_M2_vrr;
      Double I_NAI_K3x4y_S_M1_vrr = PAX*I_NAI_I2x4y_S_M1_vrr-PNX*I_NAI_I2x4y_S_M2_vrr+2*oned2z*I_NAI_Hx4y_S_M1_vrr-2*oned2z*I_NAI_Hx4y_S_M2_vrr;
      Double I_NAI_K3x3yz_S_M1_vrr = PAZ*I_NAI_I3x3y_S_M1_vrr-PNZ*I_NAI_I3x3y_S_M2_vrr;
      Double I_NAI_K3xy3z_S_M1_vrr = PAY*I_NAI_I3x3z_S_M1_vrr-PNY*I_NAI_I3x3z_S_M2_vrr;
      Double I_NAI_K3x4z_S_M1_vrr = PAX*I_NAI_I2x4z_S_M1_vrr-PNX*I_NAI_I2x4z_S_M2_vrr+2*oned2z*I_NAI_Hx4z_S_M1_vrr-2*oned2z*I_NAI_Hx4z_S_M2_vrr;
      Double I_NAI_K2x5y_S_M1_vrr = PAX*I_NAI_Ix5y_S_M1_vrr-PNX*I_NAI_Ix5y_S_M2_vrr+oned2z*I_NAI_H5y_S_M1_vrr-oned2z*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_K2x4yz_S_M1_vrr = PAZ*I_NAI_I2x4y_S_M1_vrr-PNZ*I_NAI_I2x4y_S_M2_vrr;
      Double I_NAI_K2xy4z_S_M1_vrr = PAY*I_NAI_I2x4z_S_M1_vrr-PNY*I_NAI_I2x4z_S_M2_vrr;
      Double I_NAI_K2x5z_S_M1_vrr = PAX*I_NAI_Ix5z_S_M1_vrr-PNX*I_NAI_Ix5z_S_M2_vrr+oned2z*I_NAI_H5z_S_M1_vrr-oned2z*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_Kx6y_S_M1_vrr = PAX*I_NAI_I6y_S_M1_vrr-PNX*I_NAI_I6y_S_M2_vrr;
      Double I_NAI_Kx3y3z_S_M1_vrr = PAX*I_NAI_I3y3z_S_M1_vrr-PNX*I_NAI_I3y3z_S_M2_vrr;
      Double I_NAI_Kx6z_S_M1_vrr = PAX*I_NAI_I6z_S_M1_vrr-PNX*I_NAI_I6z_S_M2_vrr;
      Double I_NAI_K7y_S_M1_vrr = PAY*I_NAI_I6y_S_M1_vrr-PNY*I_NAI_I6y_S_M2_vrr+6*oned2z*I_NAI_H5y_S_M1_vrr-6*oned2z*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_K6yz_S_M1_vrr = PAZ*I_NAI_I6y_S_M1_vrr-PNZ*I_NAI_I6y_S_M2_vrr;
      Double I_NAI_K5y2z_S_M1_vrr = PAZ*I_NAI_I5yz_S_M1_vrr-PNZ*I_NAI_I5yz_S_M2_vrr+oned2z*I_NAI_H5y_S_M1_vrr-oned2z*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_K4y3z_S_M1_vrr = PAZ*I_NAI_I4y2z_S_M1_vrr-PNZ*I_NAI_I4y2z_S_M2_vrr+2*oned2z*I_NAI_H4yz_S_M1_vrr-2*oned2z*I_NAI_H4yz_S_M2_vrr;
      Double I_NAI_K3y4z_S_M1_vrr = PAY*I_NAI_I2y4z_S_M1_vrr-PNY*I_NAI_I2y4z_S_M2_vrr+2*oned2z*I_NAI_Hy4z_S_M1_vrr-2*oned2z*I_NAI_Hy4z_S_M2_vrr;
      Double I_NAI_K2y5z_S_M1_vrr = PAY*I_NAI_Iy5z_S_M1_vrr-PNY*I_NAI_Iy5z_S_M2_vrr+oned2z*I_NAI_H5z_S_M1_vrr-oned2z*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_Ky6z_S_M1_vrr = PAY*I_NAI_I6z_S_M1_vrr-PNY*I_NAI_I6z_S_M2_vrr;
      Double I_NAI_K7z_S_M1_vrr = PAZ*I_NAI_I6z_S_M1_vrr-PNZ*I_NAI_I6z_S_M2_vrr+6*oned2z*I_NAI_H5z_S_M1_vrr-6*oned2z*I_NAI_H5z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_Px_S_vrr = PAX*I_NAI_S_S_vrr-PNX*I_NAI_S_S_M1_vrr;
      Double I_NAI_Py_S_vrr = PAY*I_NAI_S_S_vrr-PNY*I_NAI_S_S_M1_vrr;
      Double I_NAI_Pz_S_vrr = PAZ*I_NAI_S_S_vrr-PNZ*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_D2z_S_vrr = PAZ*I_NAI_Pz_S_vrr-PNZ*I_NAI_Pz_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       ************************************************************/
      Double I_NAI_F3x_S_vrr = PAX*I_NAI_D2x_S_vrr-PNX*I_NAI_D2x_S_M1_vrr+2*oned2z*I_NAI_Px_S_vrr-2*oned2z*I_NAI_Px_S_M1_vrr;
      Double I_NAI_F2xy_S_vrr = PAY*I_NAI_D2x_S_vrr-PNY*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F2xz_S_vrr = PAZ*I_NAI_D2x_S_vrr-PNZ*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F3y_S_vrr = PAY*I_NAI_D2y_S_vrr-PNY*I_NAI_D2y_S_M1_vrr+2*oned2z*I_NAI_Py_S_vrr-2*oned2z*I_NAI_Py_S_M1_vrr;
      Double I_NAI_F2yz_S_vrr = PAZ*I_NAI_D2y_S_vrr-PNZ*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_F3z_S_vrr = PAZ*I_NAI_D2z_S_vrr-PNZ*I_NAI_D2z_S_M1_vrr+2*oned2z*I_NAI_Pz_S_vrr-2*oned2z*I_NAI_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       ************************************************************/
      Double I_NAI_G4x_S_vrr = PAX*I_NAI_F3x_S_vrr-PNX*I_NAI_F3x_S_M1_vrr+3*oned2z*I_NAI_D2x_S_vrr-3*oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G3xy_S_vrr = PAY*I_NAI_F3x_S_vrr-PNY*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G3xz_S_vrr = PAZ*I_NAI_F3x_S_vrr-PNZ*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G2x2y_S_vrr = PAY*I_NAI_F2xy_S_vrr-PNY*I_NAI_F2xy_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G2x2z_S_vrr = PAZ*I_NAI_F2xz_S_vrr-PNZ*I_NAI_F2xz_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Gx3y_S_vrr = PAX*I_NAI_F3y_S_vrr-PNX*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_Gx3z_S_vrr = PAX*I_NAI_F3z_S_vrr-PNX*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4y_S_vrr = PAY*I_NAI_F3y_S_vrr-PNY*I_NAI_F3y_S_M1_vrr+3*oned2z*I_NAI_D2y_S_vrr-3*oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_G3yz_S_vrr = PAZ*I_NAI_F3y_S_vrr-PNZ*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_G2y2z_S_vrr = PAZ*I_NAI_F2yz_S_vrr-PNZ*I_NAI_F2yz_S_M1_vrr+oned2z*I_NAI_D2y_S_vrr-oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Gy3z_S_vrr = PAY*I_NAI_F3z_S_vrr-PNY*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4z_S_vrr = PAZ*I_NAI_F3z_S_vrr-PNZ*I_NAI_F3z_S_M1_vrr+3*oned2z*I_NAI_D2z_S_vrr-3*oned2z*I_NAI_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       ************************************************************/
      Double I_NAI_H5x_S_vrr = PAX*I_NAI_G4x_S_vrr-PNX*I_NAI_G4x_S_M1_vrr+4*oned2z*I_NAI_F3x_S_vrr-4*oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H4xy_S_vrr = PAY*I_NAI_G4x_S_vrr-PNY*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H4xz_S_vrr = PAZ*I_NAI_G4x_S_vrr-PNZ*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H3x2y_S_vrr = PAY*I_NAI_G3xy_S_vrr-PNY*I_NAI_G3xy_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H3xyz_S_vrr = PAZ*I_NAI_G3xy_S_vrr-PNZ*I_NAI_G3xy_S_M1_vrr;
      Double I_NAI_H3x2z_S_vrr = PAZ*I_NAI_G3xz_S_vrr-PNZ*I_NAI_G3xz_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H2x3y_S_vrr = PAX*I_NAI_Gx3y_S_vrr-PNX*I_NAI_Gx3y_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2x2yz_S_vrr = PAZ*I_NAI_G2x2y_S_vrr-PNZ*I_NAI_G2x2y_S_M1_vrr;
      Double I_NAI_H2xy2z_S_vrr = PAY*I_NAI_G2x2z_S_vrr-PNY*I_NAI_G2x2z_S_M1_vrr;
      Double I_NAI_H2x3z_S_vrr = PAX*I_NAI_Gx3z_S_vrr-PNX*I_NAI_Gx3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hx4y_S_vrr = PAX*I_NAI_G4y_S_vrr-PNX*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_Hx3yz_S_vrr = PAZ*I_NAI_Gx3y_S_vrr-PNZ*I_NAI_Gx3y_S_M1_vrr;
      Double I_NAI_Hx2y2z_S_vrr = PAX*I_NAI_G2y2z_S_vrr-PNX*I_NAI_G2y2z_S_M1_vrr;
      Double I_NAI_Hxy3z_S_vrr = PAY*I_NAI_Gx3z_S_vrr-PNY*I_NAI_Gx3z_S_M1_vrr;
      Double I_NAI_Hx4z_S_vrr = PAX*I_NAI_G4z_S_vrr-PNX*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5y_S_vrr = PAY*I_NAI_G4y_S_vrr-PNY*I_NAI_G4y_S_M1_vrr+4*oned2z*I_NAI_F3y_S_vrr-4*oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H4yz_S_vrr = PAZ*I_NAI_G4y_S_vrr-PNZ*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_H3y2z_S_vrr = PAZ*I_NAI_G3yz_S_vrr-PNZ*I_NAI_G3yz_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2y3z_S_vrr = PAY*I_NAI_Gy3z_S_vrr-PNY*I_NAI_Gy3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hy4z_S_vrr = PAY*I_NAI_G4z_S_vrr-PNY*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5z_S_vrr = PAZ*I_NAI_G4z_S_vrr-PNZ*I_NAI_G4z_S_M1_vrr+4*oned2z*I_NAI_F3z_S_vrr-4*oned2z*I_NAI_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S
       * RHS shell quartet name: SQ_NAI_H_S_M1
       * RHS shell quartet name: SQ_NAI_G_S
       * RHS shell quartet name: SQ_NAI_G_S_M1
       ************************************************************/
      Double I_NAI_I6x_S_vrr = PAX*I_NAI_H5x_S_vrr-PNX*I_NAI_H5x_S_M1_vrr+5*oned2z*I_NAI_G4x_S_vrr-5*oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I5xy_S_vrr = PAY*I_NAI_H5x_S_vrr-PNY*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_I5xz_S_vrr = PAZ*I_NAI_H5x_S_vrr-PNZ*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_I4x2y_S_vrr = PAY*I_NAI_H4xy_S_vrr-PNY*I_NAI_H4xy_S_M1_vrr+oned2z*I_NAI_G4x_S_vrr-oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I4xyz_S_vrr = PAZ*I_NAI_H4xy_S_vrr-PNZ*I_NAI_H4xy_S_M1_vrr;
      Double I_NAI_I4x2z_S_vrr = PAZ*I_NAI_H4xz_S_vrr-PNZ*I_NAI_H4xz_S_M1_vrr+oned2z*I_NAI_G4x_S_vrr-oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I3x3y_S_vrr = PAY*I_NAI_H3x2y_S_vrr-PNY*I_NAI_H3x2y_S_M1_vrr+2*oned2z*I_NAI_G3xy_S_vrr-2*oned2z*I_NAI_G3xy_S_M1_vrr;
      Double I_NAI_I3x2yz_S_vrr = PAZ*I_NAI_H3x2y_S_vrr-PNZ*I_NAI_H3x2y_S_M1_vrr;
      Double I_NAI_I3xy2z_S_vrr = PAY*I_NAI_H3x2z_S_vrr-PNY*I_NAI_H3x2z_S_M1_vrr;
      Double I_NAI_I3x3z_S_vrr = PAZ*I_NAI_H3x2z_S_vrr-PNZ*I_NAI_H3x2z_S_M1_vrr+2*oned2z*I_NAI_G3xz_S_vrr-2*oned2z*I_NAI_G3xz_S_M1_vrr;
      Double I_NAI_I2x4y_S_vrr = PAX*I_NAI_Hx4y_S_vrr-PNX*I_NAI_Hx4y_S_M1_vrr+oned2z*I_NAI_G4y_S_vrr-oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I2x3yz_S_vrr = PAZ*I_NAI_H2x3y_S_vrr-PNZ*I_NAI_H2x3y_S_M1_vrr;
      Double I_NAI_I2x2y2z_S_vrr = PAZ*I_NAI_H2x2yz_S_vrr-PNZ*I_NAI_H2x2yz_S_M1_vrr+oned2z*I_NAI_G2x2y_S_vrr-oned2z*I_NAI_G2x2y_S_M1_vrr;
      Double I_NAI_I2xy3z_S_vrr = PAY*I_NAI_H2x3z_S_vrr-PNY*I_NAI_H2x3z_S_M1_vrr;
      Double I_NAI_I2x4z_S_vrr = PAX*I_NAI_Hx4z_S_vrr-PNX*I_NAI_Hx4z_S_M1_vrr+oned2z*I_NAI_G4z_S_vrr-oned2z*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_Ix5y_S_vrr = PAX*I_NAI_H5y_S_vrr-PNX*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_Ix4yz_S_vrr = PAZ*I_NAI_Hx4y_S_vrr-PNZ*I_NAI_Hx4y_S_M1_vrr;
      Double I_NAI_Ix3y2z_S_vrr = PAX*I_NAI_H3y2z_S_vrr-PNX*I_NAI_H3y2z_S_M1_vrr;
      Double I_NAI_Ix2y3z_S_vrr = PAX*I_NAI_H2y3z_S_vrr-PNX*I_NAI_H2y3z_S_M1_vrr;
      Double I_NAI_Ixy4z_S_vrr = PAY*I_NAI_Hx4z_S_vrr-PNY*I_NAI_Hx4z_S_M1_vrr;
      Double I_NAI_Ix5z_S_vrr = PAX*I_NAI_H5z_S_vrr-PNX*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_I6y_S_vrr = PAY*I_NAI_H5y_S_vrr-PNY*I_NAI_H5y_S_M1_vrr+5*oned2z*I_NAI_G4y_S_vrr-5*oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I5yz_S_vrr = PAZ*I_NAI_H5y_S_vrr-PNZ*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_I4y2z_S_vrr = PAZ*I_NAI_H4yz_S_vrr-PNZ*I_NAI_H4yz_S_M1_vrr+oned2z*I_NAI_G4y_S_vrr-oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I3y3z_S_vrr = PAZ*I_NAI_H3y2z_S_vrr-PNZ*I_NAI_H3y2z_S_M1_vrr+2*oned2z*I_NAI_G3yz_S_vrr-2*oned2z*I_NAI_G3yz_S_M1_vrr;
      Double I_NAI_I2y4z_S_vrr = PAY*I_NAI_Hy4z_S_vrr-PNY*I_NAI_Hy4z_S_M1_vrr+oned2z*I_NAI_G4z_S_vrr-oned2z*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_Iy5z_S_vrr = PAY*I_NAI_H5z_S_vrr-PNY*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_I6z_S_vrr = PAZ*I_NAI_H5z_S_vrr-PNZ*I_NAI_H5z_S_M1_vrr+5*oned2z*I_NAI_G4z_S_vrr-5*oned2z*I_NAI_G4z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_I_S
       * RHS shell quartet name: SQ_NAI_I_S_M1
       * RHS shell quartet name: SQ_NAI_H_S
       * RHS shell quartet name: SQ_NAI_H_S_M1
       ************************************************************/
      Double I_NAI_K7x_S_vrr = PAX*I_NAI_I6x_S_vrr-PNX*I_NAI_I6x_S_M1_vrr+6*oned2z*I_NAI_H5x_S_vrr-6*oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K6xy_S_vrr = PAY*I_NAI_I6x_S_vrr-PNY*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_K6xz_S_vrr = PAZ*I_NAI_I6x_S_vrr-PNZ*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_K5x2y_S_vrr = PAY*I_NAI_I5xy_S_vrr-PNY*I_NAI_I5xy_S_M1_vrr+oned2z*I_NAI_H5x_S_vrr-oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K5xyz_S_vrr = PAZ*I_NAI_I5xy_S_vrr-PNZ*I_NAI_I5xy_S_M1_vrr;
      Double I_NAI_K5x2z_S_vrr = PAZ*I_NAI_I5xz_S_vrr-PNZ*I_NAI_I5xz_S_M1_vrr+oned2z*I_NAI_H5x_S_vrr-oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K4x3y_S_vrr = PAY*I_NAI_I4x2y_S_vrr-PNY*I_NAI_I4x2y_S_M1_vrr+2*oned2z*I_NAI_H4xy_S_vrr-2*oned2z*I_NAI_H4xy_S_M1_vrr;
      Double I_NAI_K4x2yz_S_vrr = PAZ*I_NAI_I4x2y_S_vrr-PNZ*I_NAI_I4x2y_S_M1_vrr;
      Double I_NAI_K4xy2z_S_vrr = PAY*I_NAI_I4x2z_S_vrr-PNY*I_NAI_I4x2z_S_M1_vrr;
      Double I_NAI_K4x3z_S_vrr = PAZ*I_NAI_I4x2z_S_vrr-PNZ*I_NAI_I4x2z_S_M1_vrr+2*oned2z*I_NAI_H4xz_S_vrr-2*oned2z*I_NAI_H4xz_S_M1_vrr;
      Double I_NAI_K3x4y_S_vrr = PAX*I_NAI_I2x4y_S_vrr-PNX*I_NAI_I2x4y_S_M1_vrr+2*oned2z*I_NAI_Hx4y_S_vrr-2*oned2z*I_NAI_Hx4y_S_M1_vrr;
      Double I_NAI_K3x3yz_S_vrr = PAZ*I_NAI_I3x3y_S_vrr-PNZ*I_NAI_I3x3y_S_M1_vrr;
      Double I_NAI_K3x2y2z_S_vrr = PAZ*I_NAI_I3x2yz_S_vrr-PNZ*I_NAI_I3x2yz_S_M1_vrr+oned2z*I_NAI_H3x2y_S_vrr-oned2z*I_NAI_H3x2y_S_M1_vrr;
      Double I_NAI_K3xy3z_S_vrr = PAY*I_NAI_I3x3z_S_vrr-PNY*I_NAI_I3x3z_S_M1_vrr;
      Double I_NAI_K3x4z_S_vrr = PAX*I_NAI_I2x4z_S_vrr-PNX*I_NAI_I2x4z_S_M1_vrr+2*oned2z*I_NAI_Hx4z_S_vrr-2*oned2z*I_NAI_Hx4z_S_M1_vrr;
      Double I_NAI_K2x5y_S_vrr = PAX*I_NAI_Ix5y_S_vrr-PNX*I_NAI_Ix5y_S_M1_vrr+oned2z*I_NAI_H5y_S_vrr-oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K2x4yz_S_vrr = PAZ*I_NAI_I2x4y_S_vrr-PNZ*I_NAI_I2x4y_S_M1_vrr;
      Double I_NAI_K2x3y2z_S_vrr = PAZ*I_NAI_I2x3yz_S_vrr-PNZ*I_NAI_I2x3yz_S_M1_vrr+oned2z*I_NAI_H2x3y_S_vrr-oned2z*I_NAI_H2x3y_S_M1_vrr;
      Double I_NAI_K2x2y3z_S_vrr = PAY*I_NAI_I2xy3z_S_vrr-PNY*I_NAI_I2xy3z_S_M1_vrr+oned2z*I_NAI_H2x3z_S_vrr-oned2z*I_NAI_H2x3z_S_M1_vrr;
      Double I_NAI_K2xy4z_S_vrr = PAY*I_NAI_I2x4z_S_vrr-PNY*I_NAI_I2x4z_S_M1_vrr;
      Double I_NAI_K2x5z_S_vrr = PAX*I_NAI_Ix5z_S_vrr-PNX*I_NAI_Ix5z_S_M1_vrr+oned2z*I_NAI_H5z_S_vrr-oned2z*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_Kx6y_S_vrr = PAX*I_NAI_I6y_S_vrr-PNX*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_Kx5yz_S_vrr = PAZ*I_NAI_Ix5y_S_vrr-PNZ*I_NAI_Ix5y_S_M1_vrr;
      Double I_NAI_Kx4y2z_S_vrr = PAX*I_NAI_I4y2z_S_vrr-PNX*I_NAI_I4y2z_S_M1_vrr;
      Double I_NAI_Kx3y3z_S_vrr = PAX*I_NAI_I3y3z_S_vrr-PNX*I_NAI_I3y3z_S_M1_vrr;
      Double I_NAI_Kx2y4z_S_vrr = PAX*I_NAI_I2y4z_S_vrr-PNX*I_NAI_I2y4z_S_M1_vrr;
      Double I_NAI_Kxy5z_S_vrr = PAY*I_NAI_Ix5z_S_vrr-PNY*I_NAI_Ix5z_S_M1_vrr;
      Double I_NAI_Kx6z_S_vrr = PAX*I_NAI_I6z_S_vrr-PNX*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_K7y_S_vrr = PAY*I_NAI_I6y_S_vrr-PNY*I_NAI_I6y_S_M1_vrr+6*oned2z*I_NAI_H5y_S_vrr-6*oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K6yz_S_vrr = PAZ*I_NAI_I6y_S_vrr-PNZ*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_K5y2z_S_vrr = PAZ*I_NAI_I5yz_S_vrr-PNZ*I_NAI_I5yz_S_M1_vrr+oned2z*I_NAI_H5y_S_vrr-oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K4y3z_S_vrr = PAZ*I_NAI_I4y2z_S_vrr-PNZ*I_NAI_I4y2z_S_M1_vrr+2*oned2z*I_NAI_H4yz_S_vrr-2*oned2z*I_NAI_H4yz_S_M1_vrr;
      Double I_NAI_K3y4z_S_vrr = PAY*I_NAI_I2y4z_S_vrr-PNY*I_NAI_I2y4z_S_M1_vrr+2*oned2z*I_NAI_Hy4z_S_vrr-2*oned2z*I_NAI_Hy4z_S_M1_vrr;
      Double I_NAI_K2y5z_S_vrr = PAY*I_NAI_Iy5z_S_vrr-PNY*I_NAI_Iy5z_S_M1_vrr+oned2z*I_NAI_H5z_S_vrr-oned2z*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_Ky6z_S_vrr = PAY*I_NAI_I6z_S_vrr-PNY*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_K7z_S_vrr = PAZ*I_NAI_I6z_S_vrr-PNZ*I_NAI_I6z_S_M1_vrr+6*oned2z*I_NAI_H5z_S_vrr-6*oned2z*I_NAI_H5z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_K_S
       * RHS shell quartet name: SQ_NAI_K_S_M1
       * RHS shell quartet name: SQ_NAI_I_S
       * RHS shell quartet name: SQ_NAI_I_S_M1
       ************************************************************/
      Double I_NAI_L8x_S_vrr = PAX*I_NAI_K7x_S_vrr-PNX*I_NAI_K7x_S_M1_vrr+7*oned2z*I_NAI_I6x_S_vrr-7*oned2z*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_L7xy_S_vrr = PAY*I_NAI_K7x_S_vrr-PNY*I_NAI_K7x_S_M1_vrr;
      Double I_NAI_L7xz_S_vrr = PAZ*I_NAI_K7x_S_vrr-PNZ*I_NAI_K7x_S_M1_vrr;
      Double I_NAI_L6x2y_S_vrr = PAY*I_NAI_K6xy_S_vrr-PNY*I_NAI_K6xy_S_M1_vrr+oned2z*I_NAI_I6x_S_vrr-oned2z*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_L6xyz_S_vrr = PAZ*I_NAI_K6xy_S_vrr-PNZ*I_NAI_K6xy_S_M1_vrr;
      Double I_NAI_L6x2z_S_vrr = PAZ*I_NAI_K6xz_S_vrr-PNZ*I_NAI_K6xz_S_M1_vrr+oned2z*I_NAI_I6x_S_vrr-oned2z*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_L5x3y_S_vrr = PAY*I_NAI_K5x2y_S_vrr-PNY*I_NAI_K5x2y_S_M1_vrr+2*oned2z*I_NAI_I5xy_S_vrr-2*oned2z*I_NAI_I5xy_S_M1_vrr;
      Double I_NAI_L5x2yz_S_vrr = PAZ*I_NAI_K5x2y_S_vrr-PNZ*I_NAI_K5x2y_S_M1_vrr;
      Double I_NAI_L5xy2z_S_vrr = PAY*I_NAI_K5x2z_S_vrr-PNY*I_NAI_K5x2z_S_M1_vrr;
      Double I_NAI_L5x3z_S_vrr = PAZ*I_NAI_K5x2z_S_vrr-PNZ*I_NAI_K5x2z_S_M1_vrr+2*oned2z*I_NAI_I5xz_S_vrr-2*oned2z*I_NAI_I5xz_S_M1_vrr;
      Double I_NAI_L4x4y_S_vrr = PAY*I_NAI_K4x3y_S_vrr-PNY*I_NAI_K4x3y_S_M1_vrr+3*oned2z*I_NAI_I4x2y_S_vrr-3*oned2z*I_NAI_I4x2y_S_M1_vrr;
      Double I_NAI_L4x3yz_S_vrr = PAZ*I_NAI_K4x3y_S_vrr-PNZ*I_NAI_K4x3y_S_M1_vrr;
      Double I_NAI_L4x2y2z_S_vrr = PAZ*I_NAI_K4x2yz_S_vrr-PNZ*I_NAI_K4x2yz_S_M1_vrr+oned2z*I_NAI_I4x2y_S_vrr-oned2z*I_NAI_I4x2y_S_M1_vrr;
      Double I_NAI_L4xy3z_S_vrr = PAY*I_NAI_K4x3z_S_vrr-PNY*I_NAI_K4x3z_S_M1_vrr;
      Double I_NAI_L4x4z_S_vrr = PAZ*I_NAI_K4x3z_S_vrr-PNZ*I_NAI_K4x3z_S_M1_vrr+3*oned2z*I_NAI_I4x2z_S_vrr-3*oned2z*I_NAI_I4x2z_S_M1_vrr;
      Double I_NAI_L3x5y_S_vrr = PAX*I_NAI_K2x5y_S_vrr-PNX*I_NAI_K2x5y_S_M1_vrr+2*oned2z*I_NAI_Ix5y_S_vrr-2*oned2z*I_NAI_Ix5y_S_M1_vrr;
      Double I_NAI_L3x4yz_S_vrr = PAZ*I_NAI_K3x4y_S_vrr-PNZ*I_NAI_K3x4y_S_M1_vrr;
      Double I_NAI_L3x3y2z_S_vrr = PAZ*I_NAI_K3x3yz_S_vrr-PNZ*I_NAI_K3x3yz_S_M1_vrr+oned2z*I_NAI_I3x3y_S_vrr-oned2z*I_NAI_I3x3y_S_M1_vrr;
      Double I_NAI_L3x2y3z_S_vrr = PAY*I_NAI_K3xy3z_S_vrr-PNY*I_NAI_K3xy3z_S_M1_vrr+oned2z*I_NAI_I3x3z_S_vrr-oned2z*I_NAI_I3x3z_S_M1_vrr;
      Double I_NAI_L3xy4z_S_vrr = PAY*I_NAI_K3x4z_S_vrr-PNY*I_NAI_K3x4z_S_M1_vrr;
      Double I_NAI_L3x5z_S_vrr = PAX*I_NAI_K2x5z_S_vrr-PNX*I_NAI_K2x5z_S_M1_vrr+2*oned2z*I_NAI_Ix5z_S_vrr-2*oned2z*I_NAI_Ix5z_S_M1_vrr;
      Double I_NAI_L2x6y_S_vrr = PAX*I_NAI_Kx6y_S_vrr-PNX*I_NAI_Kx6y_S_M1_vrr+oned2z*I_NAI_I6y_S_vrr-oned2z*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_L2x5yz_S_vrr = PAZ*I_NAI_K2x5y_S_vrr-PNZ*I_NAI_K2x5y_S_M1_vrr;
      Double I_NAI_L2x4y2z_S_vrr = PAZ*I_NAI_K2x4yz_S_vrr-PNZ*I_NAI_K2x4yz_S_M1_vrr+oned2z*I_NAI_I2x4y_S_vrr-oned2z*I_NAI_I2x4y_S_M1_vrr;
      Double I_NAI_L2x3y3z_S_vrr = PAX*I_NAI_Kx3y3z_S_vrr-PNX*I_NAI_Kx3y3z_S_M1_vrr+oned2z*I_NAI_I3y3z_S_vrr-oned2z*I_NAI_I3y3z_S_M1_vrr;
      Double I_NAI_L2x2y4z_S_vrr = PAY*I_NAI_K2xy4z_S_vrr-PNY*I_NAI_K2xy4z_S_M1_vrr+oned2z*I_NAI_I2x4z_S_vrr-oned2z*I_NAI_I2x4z_S_M1_vrr;
      Double I_NAI_L2xy5z_S_vrr = PAY*I_NAI_K2x5z_S_vrr-PNY*I_NAI_K2x5z_S_M1_vrr;
      Double I_NAI_L2x6z_S_vrr = PAX*I_NAI_Kx6z_S_vrr-PNX*I_NAI_Kx6z_S_M1_vrr+oned2z*I_NAI_I6z_S_vrr-oned2z*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_Lx7y_S_vrr = PAX*I_NAI_K7y_S_vrr-PNX*I_NAI_K7y_S_M1_vrr;
      Double I_NAI_Lx6yz_S_vrr = PAZ*I_NAI_Kx6y_S_vrr-PNZ*I_NAI_Kx6y_S_M1_vrr;
      Double I_NAI_Lx5y2z_S_vrr = PAX*I_NAI_K5y2z_S_vrr-PNX*I_NAI_K5y2z_S_M1_vrr;
      Double I_NAI_Lx4y3z_S_vrr = PAX*I_NAI_K4y3z_S_vrr-PNX*I_NAI_K4y3z_S_M1_vrr;
      Double I_NAI_Lx3y4z_S_vrr = PAX*I_NAI_K3y4z_S_vrr-PNX*I_NAI_K3y4z_S_M1_vrr;
      Double I_NAI_Lx2y5z_S_vrr = PAX*I_NAI_K2y5z_S_vrr-PNX*I_NAI_K2y5z_S_M1_vrr;
      Double I_NAI_Lxy6z_S_vrr = PAY*I_NAI_Kx6z_S_vrr-PNY*I_NAI_Kx6z_S_M1_vrr;
      Double I_NAI_Lx7z_S_vrr = PAX*I_NAI_K7z_S_vrr-PNX*I_NAI_K7z_S_M1_vrr;
      Double I_NAI_L8y_S_vrr = PAY*I_NAI_K7y_S_vrr-PNY*I_NAI_K7y_S_M1_vrr+7*oned2z*I_NAI_I6y_S_vrr-7*oned2z*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_L7yz_S_vrr = PAZ*I_NAI_K7y_S_vrr-PNZ*I_NAI_K7y_S_M1_vrr;
      Double I_NAI_L6y2z_S_vrr = PAZ*I_NAI_K6yz_S_vrr-PNZ*I_NAI_K6yz_S_M1_vrr+oned2z*I_NAI_I6y_S_vrr-oned2z*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_L5y3z_S_vrr = PAZ*I_NAI_K5y2z_S_vrr-PNZ*I_NAI_K5y2z_S_M1_vrr+2*oned2z*I_NAI_I5yz_S_vrr-2*oned2z*I_NAI_I5yz_S_M1_vrr;
      Double I_NAI_L4y4z_S_vrr = PAZ*I_NAI_K4y3z_S_vrr-PNZ*I_NAI_K4y3z_S_M1_vrr+3*oned2z*I_NAI_I4y2z_S_vrr-3*oned2z*I_NAI_I4y2z_S_M1_vrr;
      Double I_NAI_L3y5z_S_vrr = PAY*I_NAI_K2y5z_S_vrr-PNY*I_NAI_K2y5z_S_M1_vrr+2*oned2z*I_NAI_Iy5z_S_vrr-2*oned2z*I_NAI_Iy5z_S_M1_vrr;
      Double I_NAI_L2y6z_S_vrr = PAY*I_NAI_Ky6z_S_vrr-PNY*I_NAI_Ky6z_S_M1_vrr+oned2z*I_NAI_I6z_S_vrr-oned2z*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_Ly7z_S_vrr = PAY*I_NAI_K7z_S_vrr-PNY*I_NAI_K7z_S_M1_vrr;
      Double I_NAI_L8z_S_vrr = PAZ*I_NAI_K7z_S_vrr-PNZ*I_NAI_K7z_S_M1_vrr+7*oned2z*I_NAI_I6z_S_vrr-7*oned2z*I_NAI_I6z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_L8x_S += I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S += I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S += I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S += I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S += I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S += I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S += I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S += I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S += I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S += I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S += I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S += I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S += I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S += I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S += I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S += I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S += I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S += I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S += I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S += I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S += I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S += I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S += I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S += I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S += I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S += I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S += I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S += I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S += I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S += I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S += I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S += I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S += I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S += I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S += I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S += I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S += I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S += I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S += I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S += I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S += I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S += I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S += I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S += I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S += I_NAI_L8z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_K7x_S += I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S += I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S += I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S += I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S += I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S += I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S += I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S += I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S += I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S += I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S += I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S += I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S += I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S += I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S += I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S += I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S += I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S += I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S += I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S += I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S += I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S += I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S += I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S += I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S += I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S += I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S += I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S += I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S += I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S += I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S += I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S += I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S += I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S += I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S += I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S += I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_I6x_S += I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S += I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S += I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S += I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S += I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S += I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S += I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S += I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S += I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S += I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S += I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S += I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S += I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S += I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S += I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S += I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S += I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S += I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S += I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S += I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S += I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S += I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S += I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S += I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S += I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S += I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S += I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S += I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_H5x_S += I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S += I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S += I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S += I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S += I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S += I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S += I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S += I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S += I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S += I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S += I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S += I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S += I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S += I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S += I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S += I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S += I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S += I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S += I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S += I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S += I_NAI_H5z_S_vrr;
    }
  }

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
   * shell quartet name: SQ_NAI_H_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S
   * RHS shell quartet name: SQ_NAI_H_S
   ************************************************************/
  Double I_NAI_H5x_Px = I_NAI_I6x_S+ABX*I_NAI_H5x_S;
  Double I_NAI_H4xy_Px = I_NAI_I5xy_S+ABX*I_NAI_H4xy_S;
  Double I_NAI_H4xz_Px = I_NAI_I5xz_S+ABX*I_NAI_H4xz_S;
  Double I_NAI_H3x2y_Px = I_NAI_I4x2y_S+ABX*I_NAI_H3x2y_S;
  Double I_NAI_H3xyz_Px = I_NAI_I4xyz_S+ABX*I_NAI_H3xyz_S;
  Double I_NAI_H3x2z_Px = I_NAI_I4x2z_S+ABX*I_NAI_H3x2z_S;
  Double I_NAI_H2x3y_Px = I_NAI_I3x3y_S+ABX*I_NAI_H2x3y_S;
  Double I_NAI_H2x2yz_Px = I_NAI_I3x2yz_S+ABX*I_NAI_H2x2yz_S;
  Double I_NAI_H2xy2z_Px = I_NAI_I3xy2z_S+ABX*I_NAI_H2xy2z_S;
  Double I_NAI_H2x3z_Px = I_NAI_I3x3z_S+ABX*I_NAI_H2x3z_S;
  Double I_NAI_Hx4y_Px = I_NAI_I2x4y_S+ABX*I_NAI_Hx4y_S;
  Double I_NAI_Hx3yz_Px = I_NAI_I2x3yz_S+ABX*I_NAI_Hx3yz_S;
  Double I_NAI_Hx2y2z_Px = I_NAI_I2x2y2z_S+ABX*I_NAI_Hx2y2z_S;
  Double I_NAI_Hxy3z_Px = I_NAI_I2xy3z_S+ABX*I_NAI_Hxy3z_S;
  Double I_NAI_Hx4z_Px = I_NAI_I2x4z_S+ABX*I_NAI_Hx4z_S;
  Double I_NAI_H5y_Px = I_NAI_Ix5y_S+ABX*I_NAI_H5y_S;
  Double I_NAI_H4yz_Px = I_NAI_Ix4yz_S+ABX*I_NAI_H4yz_S;
  Double I_NAI_H3y2z_Px = I_NAI_Ix3y2z_S+ABX*I_NAI_H3y2z_S;
  Double I_NAI_H2y3z_Px = I_NAI_Ix2y3z_S+ABX*I_NAI_H2y3z_S;
  Double I_NAI_Hy4z_Px = I_NAI_Ixy4z_S+ABX*I_NAI_Hy4z_S;
  Double I_NAI_H5z_Px = I_NAI_Ix5z_S+ABX*I_NAI_H5z_S;
  Double I_NAI_H5x_Py = I_NAI_I5xy_S+ABY*I_NAI_H5x_S;
  Double I_NAI_H4xy_Py = I_NAI_I4x2y_S+ABY*I_NAI_H4xy_S;
  Double I_NAI_H4xz_Py = I_NAI_I4xyz_S+ABY*I_NAI_H4xz_S;
  Double I_NAI_H3x2y_Py = I_NAI_I3x3y_S+ABY*I_NAI_H3x2y_S;
  Double I_NAI_H3xyz_Py = I_NAI_I3x2yz_S+ABY*I_NAI_H3xyz_S;
  Double I_NAI_H3x2z_Py = I_NAI_I3xy2z_S+ABY*I_NAI_H3x2z_S;
  Double I_NAI_H2x3y_Py = I_NAI_I2x4y_S+ABY*I_NAI_H2x3y_S;
  Double I_NAI_H2x2yz_Py = I_NAI_I2x3yz_S+ABY*I_NAI_H2x2yz_S;
  Double I_NAI_H2xy2z_Py = I_NAI_I2x2y2z_S+ABY*I_NAI_H2xy2z_S;
  Double I_NAI_H2x3z_Py = I_NAI_I2xy3z_S+ABY*I_NAI_H2x3z_S;
  Double I_NAI_Hx4y_Py = I_NAI_Ix5y_S+ABY*I_NAI_Hx4y_S;
  Double I_NAI_Hx3yz_Py = I_NAI_Ix4yz_S+ABY*I_NAI_Hx3yz_S;
  Double I_NAI_Hx2y2z_Py = I_NAI_Ix3y2z_S+ABY*I_NAI_Hx2y2z_S;
  Double I_NAI_Hxy3z_Py = I_NAI_Ix2y3z_S+ABY*I_NAI_Hxy3z_S;
  Double I_NAI_Hx4z_Py = I_NAI_Ixy4z_S+ABY*I_NAI_Hx4z_S;
  Double I_NAI_H5y_Py = I_NAI_I6y_S+ABY*I_NAI_H5y_S;
  Double I_NAI_H4yz_Py = I_NAI_I5yz_S+ABY*I_NAI_H4yz_S;
  Double I_NAI_H3y2z_Py = I_NAI_I4y2z_S+ABY*I_NAI_H3y2z_S;
  Double I_NAI_H2y3z_Py = I_NAI_I3y3z_S+ABY*I_NAI_H2y3z_S;
  Double I_NAI_Hy4z_Py = I_NAI_I2y4z_S+ABY*I_NAI_Hy4z_S;
  Double I_NAI_H5z_Py = I_NAI_Iy5z_S+ABY*I_NAI_H5z_S;
  Double I_NAI_H5x_Pz = I_NAI_I5xz_S+ABZ*I_NAI_H5x_S;
  Double I_NAI_H4xy_Pz = I_NAI_I4xyz_S+ABZ*I_NAI_H4xy_S;
  Double I_NAI_H4xz_Pz = I_NAI_I4x2z_S+ABZ*I_NAI_H4xz_S;
  Double I_NAI_H3x2y_Pz = I_NAI_I3x2yz_S+ABZ*I_NAI_H3x2y_S;
  Double I_NAI_H3xyz_Pz = I_NAI_I3xy2z_S+ABZ*I_NAI_H3xyz_S;
  Double I_NAI_H3x2z_Pz = I_NAI_I3x3z_S+ABZ*I_NAI_H3x2z_S;
  Double I_NAI_H2x3y_Pz = I_NAI_I2x3yz_S+ABZ*I_NAI_H2x3y_S;
  Double I_NAI_H2x2yz_Pz = I_NAI_I2x2y2z_S+ABZ*I_NAI_H2x2yz_S;
  Double I_NAI_H2xy2z_Pz = I_NAI_I2xy3z_S+ABZ*I_NAI_H2xy2z_S;
  Double I_NAI_H2x3z_Pz = I_NAI_I2x4z_S+ABZ*I_NAI_H2x3z_S;
  Double I_NAI_Hx4y_Pz = I_NAI_Ix4yz_S+ABZ*I_NAI_Hx4y_S;
  Double I_NAI_Hx3yz_Pz = I_NAI_Ix3y2z_S+ABZ*I_NAI_Hx3yz_S;
  Double I_NAI_Hx2y2z_Pz = I_NAI_Ix2y3z_S+ABZ*I_NAI_Hx2y2z_S;
  Double I_NAI_Hxy3z_Pz = I_NAI_Ixy4z_S+ABZ*I_NAI_Hxy3z_S;
  Double I_NAI_Hx4z_Pz = I_NAI_Ix5z_S+ABZ*I_NAI_Hx4z_S;
  Double I_NAI_H5y_Pz = I_NAI_I5yz_S+ABZ*I_NAI_H5y_S;
  Double I_NAI_H4yz_Pz = I_NAI_I4y2z_S+ABZ*I_NAI_H4yz_S;
  Double I_NAI_H3y2z_Pz = I_NAI_I3y3z_S+ABZ*I_NAI_H3y2z_S;
  Double I_NAI_H2y3z_Pz = I_NAI_I2y4z_S+ABZ*I_NAI_H2y3z_S;
  Double I_NAI_Hy4z_Pz = I_NAI_Iy5z_S+ABZ*I_NAI_Hy4z_S;
  Double I_NAI_H5z_Pz = I_NAI_I6z_S+ABZ*I_NAI_H5z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S
   * RHS shell quartet name: SQ_NAI_I_S
   ************************************************************/
  Double I_NAI_I6x_Px = I_NAI_K7x_S+ABX*I_NAI_I6x_S;
  Double I_NAI_I5xy_Px = I_NAI_K6xy_S+ABX*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Px = I_NAI_K6xz_S+ABX*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Px = I_NAI_K5x2y_S+ABX*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Px = I_NAI_K5xyz_S+ABX*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Px = I_NAI_K5x2z_S+ABX*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Px = I_NAI_K4x3y_S+ABX*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Px = I_NAI_K4x2yz_S+ABX*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Px = I_NAI_K4xy2z_S+ABX*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Px = I_NAI_K4x3z_S+ABX*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Px = I_NAI_K3x4y_S+ABX*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Px = I_NAI_K3x3yz_S+ABX*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Px = I_NAI_K3x2y2z_S+ABX*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Px = I_NAI_K3xy3z_S+ABX*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Px = I_NAI_K3x4z_S+ABX*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Px = I_NAI_K2x5y_S+ABX*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Px = I_NAI_K2x4yz_S+ABX*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Px = I_NAI_K2x3y2z_S+ABX*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Px = I_NAI_K2x2y3z_S+ABX*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Px = I_NAI_K2xy4z_S+ABX*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Px = I_NAI_K2x5z_S+ABX*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Px = I_NAI_Kx6y_S+ABX*I_NAI_I6y_S;
  Double I_NAI_I5yz_Px = I_NAI_Kx5yz_S+ABX*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Px = I_NAI_Kx4y2z_S+ABX*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Px = I_NAI_Kx3y3z_S+ABX*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Px = I_NAI_Kx2y4z_S+ABX*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Px = I_NAI_Kxy5z_S+ABX*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Px = I_NAI_Kx6z_S+ABX*I_NAI_I6z_S;
  Double I_NAI_I6x_Py = I_NAI_K6xy_S+ABY*I_NAI_I6x_S;
  Double I_NAI_I5xy_Py = I_NAI_K5x2y_S+ABY*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Py = I_NAI_K5xyz_S+ABY*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Py = I_NAI_K4x3y_S+ABY*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Py = I_NAI_K4x2yz_S+ABY*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Py = I_NAI_K4xy2z_S+ABY*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Py = I_NAI_K3x4y_S+ABY*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Py = I_NAI_K3x3yz_S+ABY*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Py = I_NAI_K3x2y2z_S+ABY*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Py = I_NAI_K3xy3z_S+ABY*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Py = I_NAI_K2x5y_S+ABY*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Py = I_NAI_K2x4yz_S+ABY*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Py = I_NAI_K2x3y2z_S+ABY*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Py = I_NAI_K2x2y3z_S+ABY*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Py = I_NAI_K2xy4z_S+ABY*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Py = I_NAI_Kx6y_S+ABY*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Py = I_NAI_Kx5yz_S+ABY*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Py = I_NAI_Kx4y2z_S+ABY*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Py = I_NAI_Kx3y3z_S+ABY*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Py = I_NAI_Kx2y4z_S+ABY*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Py = I_NAI_Kxy5z_S+ABY*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Py = I_NAI_K7y_S+ABY*I_NAI_I6y_S;
  Double I_NAI_I5yz_Py = I_NAI_K6yz_S+ABY*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Py = I_NAI_K5y2z_S+ABY*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Py = I_NAI_K4y3z_S+ABY*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Py = I_NAI_K3y4z_S+ABY*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Py = I_NAI_K2y5z_S+ABY*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Py = I_NAI_Ky6z_S+ABY*I_NAI_I6z_S;
  Double I_NAI_I6x_Pz = I_NAI_K6xz_S+ABZ*I_NAI_I6x_S;
  Double I_NAI_I5xy_Pz = I_NAI_K5xyz_S+ABZ*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Pz = I_NAI_K5x2z_S+ABZ*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Pz = I_NAI_K4x2yz_S+ABZ*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Pz = I_NAI_K4xy2z_S+ABZ*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Pz = I_NAI_K4x3z_S+ABZ*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Pz = I_NAI_K3x3yz_S+ABZ*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Pz = I_NAI_K3x2y2z_S+ABZ*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Pz = I_NAI_K3xy3z_S+ABZ*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Pz = I_NAI_K3x4z_S+ABZ*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Pz = I_NAI_K2x4yz_S+ABZ*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Pz = I_NAI_K2x3y2z_S+ABZ*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Pz = I_NAI_K2x2y3z_S+ABZ*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Pz = I_NAI_K2xy4z_S+ABZ*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Pz = I_NAI_K2x5z_S+ABZ*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Pz = I_NAI_Kx5yz_S+ABZ*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Pz = I_NAI_Kx4y2z_S+ABZ*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Pz = I_NAI_Kx3y3z_S+ABZ*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Pz = I_NAI_Kx2y4z_S+ABZ*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Pz = I_NAI_Kxy5z_S+ABZ*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Pz = I_NAI_Kx6z_S+ABZ*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Pz = I_NAI_K6yz_S+ABZ*I_NAI_I6y_S;
  Double I_NAI_I5yz_Pz = I_NAI_K5y2z_S+ABZ*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Pz = I_NAI_K4y3z_S+ABZ*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Pz = I_NAI_K3y4z_S+ABZ*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Pz = I_NAI_K2y5z_S+ABZ*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Pz = I_NAI_Ky6z_S+ABZ*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Pz = I_NAI_K7z_S+ABZ*I_NAI_I6z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 42 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P
   * RHS shell quartet name: SQ_NAI_H_P
   ************************************************************/
  Double I_NAI_H5x_D2x = I_NAI_I6x_Px+ABX*I_NAI_H5x_Px;
  Double I_NAI_H4xy_D2x = I_NAI_I5xy_Px+ABX*I_NAI_H4xy_Px;
  Double I_NAI_H4xz_D2x = I_NAI_I5xz_Px+ABX*I_NAI_H4xz_Px;
  Double I_NAI_H3x2y_D2x = I_NAI_I4x2y_Px+ABX*I_NAI_H3x2y_Px;
  Double I_NAI_H3xyz_D2x = I_NAI_I4xyz_Px+ABX*I_NAI_H3xyz_Px;
  Double I_NAI_H3x2z_D2x = I_NAI_I4x2z_Px+ABX*I_NAI_H3x2z_Px;
  Double I_NAI_H2x3y_D2x = I_NAI_I3x3y_Px+ABX*I_NAI_H2x3y_Px;
  Double I_NAI_H2x2yz_D2x = I_NAI_I3x2yz_Px+ABX*I_NAI_H2x2yz_Px;
  Double I_NAI_H2xy2z_D2x = I_NAI_I3xy2z_Px+ABX*I_NAI_H2xy2z_Px;
  Double I_NAI_H2x3z_D2x = I_NAI_I3x3z_Px+ABX*I_NAI_H2x3z_Px;
  Double I_NAI_Hx4y_D2x = I_NAI_I2x4y_Px+ABX*I_NAI_Hx4y_Px;
  Double I_NAI_Hx3yz_D2x = I_NAI_I2x3yz_Px+ABX*I_NAI_Hx3yz_Px;
  Double I_NAI_Hx2y2z_D2x = I_NAI_I2x2y2z_Px+ABX*I_NAI_Hx2y2z_Px;
  Double I_NAI_Hxy3z_D2x = I_NAI_I2xy3z_Px+ABX*I_NAI_Hxy3z_Px;
  Double I_NAI_Hx4z_D2x = I_NAI_I2x4z_Px+ABX*I_NAI_Hx4z_Px;
  Double I_NAI_H5y_D2x = I_NAI_Ix5y_Px+ABX*I_NAI_H5y_Px;
  Double I_NAI_H4yz_D2x = I_NAI_Ix4yz_Px+ABX*I_NAI_H4yz_Px;
  Double I_NAI_H3y2z_D2x = I_NAI_Ix3y2z_Px+ABX*I_NAI_H3y2z_Px;
  Double I_NAI_H2y3z_D2x = I_NAI_Ix2y3z_Px+ABX*I_NAI_H2y3z_Px;
  Double I_NAI_Hy4z_D2x = I_NAI_Ixy4z_Px+ABX*I_NAI_Hy4z_Px;
  Double I_NAI_H5z_D2x = I_NAI_Ix5z_Px+ABX*I_NAI_H5z_Px;
  Double I_NAI_H5x_Dxy = I_NAI_I5xy_Px+ABY*I_NAI_H5x_Px;
  Double I_NAI_H4xy_Dxy = I_NAI_I4x2y_Px+ABY*I_NAI_H4xy_Px;
  Double I_NAI_H4xz_Dxy = I_NAI_I4xyz_Px+ABY*I_NAI_H4xz_Px;
  Double I_NAI_H3x2y_Dxy = I_NAI_I3x3y_Px+ABY*I_NAI_H3x2y_Px;
  Double I_NAI_H3xyz_Dxy = I_NAI_I3x2yz_Px+ABY*I_NAI_H3xyz_Px;
  Double I_NAI_H3x2z_Dxy = I_NAI_I3xy2z_Px+ABY*I_NAI_H3x2z_Px;
  Double I_NAI_H2x3y_Dxy = I_NAI_I2x4y_Px+ABY*I_NAI_H2x3y_Px;
  Double I_NAI_H2x2yz_Dxy = I_NAI_I2x3yz_Px+ABY*I_NAI_H2x2yz_Px;
  Double I_NAI_H2xy2z_Dxy = I_NAI_I2x2y2z_Px+ABY*I_NAI_H2xy2z_Px;
  Double I_NAI_H2x3z_Dxy = I_NAI_I2xy3z_Px+ABY*I_NAI_H2x3z_Px;
  Double I_NAI_Hx4y_Dxy = I_NAI_Ix5y_Px+ABY*I_NAI_Hx4y_Px;
  Double I_NAI_Hx3yz_Dxy = I_NAI_Ix4yz_Px+ABY*I_NAI_Hx3yz_Px;
  Double I_NAI_Hx2y2z_Dxy = I_NAI_Ix3y2z_Px+ABY*I_NAI_Hx2y2z_Px;
  Double I_NAI_Hxy3z_Dxy = I_NAI_Ix2y3z_Px+ABY*I_NAI_Hxy3z_Px;
  Double I_NAI_Hx4z_Dxy = I_NAI_Ixy4z_Px+ABY*I_NAI_Hx4z_Px;
  Double I_NAI_H5y_Dxy = I_NAI_I6y_Px+ABY*I_NAI_H5y_Px;
  Double I_NAI_H4yz_Dxy = I_NAI_I5yz_Px+ABY*I_NAI_H4yz_Px;
  Double I_NAI_H3y2z_Dxy = I_NAI_I4y2z_Px+ABY*I_NAI_H3y2z_Px;
  Double I_NAI_H2y3z_Dxy = I_NAI_I3y3z_Px+ABY*I_NAI_H2y3z_Px;
  Double I_NAI_Hy4z_Dxy = I_NAI_I2y4z_Px+ABY*I_NAI_Hy4z_Px;
  Double I_NAI_H5z_Dxy = I_NAI_Iy5z_Px+ABY*I_NAI_H5z_Px;
  Double I_NAI_H5x_D2y = I_NAI_I5xy_Py+ABY*I_NAI_H5x_Py;
  Double I_NAI_H4xy_D2y = I_NAI_I4x2y_Py+ABY*I_NAI_H4xy_Py;
  Double I_NAI_H4xz_D2y = I_NAI_I4xyz_Py+ABY*I_NAI_H4xz_Py;
  Double I_NAI_H3x2y_D2y = I_NAI_I3x3y_Py+ABY*I_NAI_H3x2y_Py;
  Double I_NAI_H3xyz_D2y = I_NAI_I3x2yz_Py+ABY*I_NAI_H3xyz_Py;
  Double I_NAI_H3x2z_D2y = I_NAI_I3xy2z_Py+ABY*I_NAI_H3x2z_Py;
  Double I_NAI_H2x3y_D2y = I_NAI_I2x4y_Py+ABY*I_NAI_H2x3y_Py;
  Double I_NAI_H2x2yz_D2y = I_NAI_I2x3yz_Py+ABY*I_NAI_H2x2yz_Py;
  Double I_NAI_H2xy2z_D2y = I_NAI_I2x2y2z_Py+ABY*I_NAI_H2xy2z_Py;
  Double I_NAI_H2x3z_D2y = I_NAI_I2xy3z_Py+ABY*I_NAI_H2x3z_Py;
  Double I_NAI_Hx4y_D2y = I_NAI_Ix5y_Py+ABY*I_NAI_Hx4y_Py;
  Double I_NAI_Hx3yz_D2y = I_NAI_Ix4yz_Py+ABY*I_NAI_Hx3yz_Py;
  Double I_NAI_Hx2y2z_D2y = I_NAI_Ix3y2z_Py+ABY*I_NAI_Hx2y2z_Py;
  Double I_NAI_Hxy3z_D2y = I_NAI_Ix2y3z_Py+ABY*I_NAI_Hxy3z_Py;
  Double I_NAI_Hx4z_D2y = I_NAI_Ixy4z_Py+ABY*I_NAI_Hx4z_Py;
  Double I_NAI_H5y_D2y = I_NAI_I6y_Py+ABY*I_NAI_H5y_Py;
  Double I_NAI_H4yz_D2y = I_NAI_I5yz_Py+ABY*I_NAI_H4yz_Py;
  Double I_NAI_H3y2z_D2y = I_NAI_I4y2z_Py+ABY*I_NAI_H3y2z_Py;
  Double I_NAI_H2y3z_D2y = I_NAI_I3y3z_Py+ABY*I_NAI_H2y3z_Py;
  Double I_NAI_Hy4z_D2y = I_NAI_I2y4z_Py+ABY*I_NAI_Hy4z_Py;
  Double I_NAI_H5z_D2y = I_NAI_Iy5z_Py+ABY*I_NAI_H5z_Py;
  Double I_NAI_H5x_D2z = I_NAI_I5xz_Pz+ABZ*I_NAI_H5x_Pz;
  Double I_NAI_H4xy_D2z = I_NAI_I4xyz_Pz+ABZ*I_NAI_H4xy_Pz;
  Double I_NAI_H4xz_D2z = I_NAI_I4x2z_Pz+ABZ*I_NAI_H4xz_Pz;
  Double I_NAI_H3x2y_D2z = I_NAI_I3x2yz_Pz+ABZ*I_NAI_H3x2y_Pz;
  Double I_NAI_H3xyz_D2z = I_NAI_I3xy2z_Pz+ABZ*I_NAI_H3xyz_Pz;
  Double I_NAI_H3x2z_D2z = I_NAI_I3x3z_Pz+ABZ*I_NAI_H3x2z_Pz;
  Double I_NAI_H2x3y_D2z = I_NAI_I2x3yz_Pz+ABZ*I_NAI_H2x3y_Pz;
  Double I_NAI_H2x2yz_D2z = I_NAI_I2x2y2z_Pz+ABZ*I_NAI_H2x2yz_Pz;
  Double I_NAI_H2xy2z_D2z = I_NAI_I2xy3z_Pz+ABZ*I_NAI_H2xy2z_Pz;
  Double I_NAI_H2x3z_D2z = I_NAI_I2x4z_Pz+ABZ*I_NAI_H2x3z_Pz;
  Double I_NAI_Hx4y_D2z = I_NAI_Ix4yz_Pz+ABZ*I_NAI_Hx4y_Pz;
  Double I_NAI_Hx3yz_D2z = I_NAI_Ix3y2z_Pz+ABZ*I_NAI_Hx3yz_Pz;
  Double I_NAI_Hx2y2z_D2z = I_NAI_Ix2y3z_Pz+ABZ*I_NAI_Hx2y2z_Pz;
  Double I_NAI_Hxy3z_D2z = I_NAI_Ixy4z_Pz+ABZ*I_NAI_Hxy3z_Pz;
  Double I_NAI_Hx4z_D2z = I_NAI_Ix5z_Pz+ABZ*I_NAI_Hx4z_Pz;
  Double I_NAI_H5y_D2z = I_NAI_I5yz_Pz+ABZ*I_NAI_H5y_Pz;
  Double I_NAI_H4yz_D2z = I_NAI_I4y2z_Pz+ABZ*I_NAI_H4yz_Pz;
  Double I_NAI_H3y2z_D2z = I_NAI_I3y3z_Pz+ABZ*I_NAI_H3y2z_Pz;
  Double I_NAI_H2y3z_D2z = I_NAI_I2y4z_Pz+ABZ*I_NAI_H2y3z_Pz;
  Double I_NAI_Hy4z_D2z = I_NAI_Iy5z_Pz+ABZ*I_NAI_Hy4z_Pz;
  Double I_NAI_H5z_D2z = I_NAI_I6z_Pz+ABZ*I_NAI_H5z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_K_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S
   * RHS shell quartet name: SQ_NAI_K_S
   ************************************************************/
  Double I_NAI_K7x_Px = I_NAI_L8x_S+ABX*I_NAI_K7x_S;
  Double I_NAI_K6xy_Px = I_NAI_L7xy_S+ABX*I_NAI_K6xy_S;
  Double I_NAI_K6xz_Px = I_NAI_L7xz_S+ABX*I_NAI_K6xz_S;
  Double I_NAI_K5x2y_Px = I_NAI_L6x2y_S+ABX*I_NAI_K5x2y_S;
  Double I_NAI_K5xyz_Px = I_NAI_L6xyz_S+ABX*I_NAI_K5xyz_S;
  Double I_NAI_K5x2z_Px = I_NAI_L6x2z_S+ABX*I_NAI_K5x2z_S;
  Double I_NAI_K4x3y_Px = I_NAI_L5x3y_S+ABX*I_NAI_K4x3y_S;
  Double I_NAI_K4x2yz_Px = I_NAI_L5x2yz_S+ABX*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Px = I_NAI_L5xy2z_S+ABX*I_NAI_K4xy2z_S;
  Double I_NAI_K4x3z_Px = I_NAI_L5x3z_S+ABX*I_NAI_K4x3z_S;
  Double I_NAI_K3x4y_Px = I_NAI_L4x4y_S+ABX*I_NAI_K3x4y_S;
  Double I_NAI_K3x3yz_Px = I_NAI_L4x3yz_S+ABX*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Px = I_NAI_L4x2y2z_S+ABX*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Px = I_NAI_L4xy3z_S+ABX*I_NAI_K3xy3z_S;
  Double I_NAI_K3x4z_Px = I_NAI_L4x4z_S+ABX*I_NAI_K3x4z_S;
  Double I_NAI_K2x5y_Px = I_NAI_L3x5y_S+ABX*I_NAI_K2x5y_S;
  Double I_NAI_K2x4yz_Px = I_NAI_L3x4yz_S+ABX*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Px = I_NAI_L3x3y2z_S+ABX*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Px = I_NAI_L3x2y3z_S+ABX*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Px = I_NAI_L3xy4z_S+ABX*I_NAI_K2xy4z_S;
  Double I_NAI_K2x5z_Px = I_NAI_L3x5z_S+ABX*I_NAI_K2x5z_S;
  Double I_NAI_Kx6y_Px = I_NAI_L2x6y_S+ABX*I_NAI_Kx6y_S;
  Double I_NAI_Kx5yz_Px = I_NAI_L2x5yz_S+ABX*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Px = I_NAI_L2x4y2z_S+ABX*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Px = I_NAI_L2x3y3z_S+ABX*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Px = I_NAI_L2x2y4z_S+ABX*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Px = I_NAI_L2xy5z_S+ABX*I_NAI_Kxy5z_S;
  Double I_NAI_Kx6z_Px = I_NAI_L2x6z_S+ABX*I_NAI_Kx6z_S;
  Double I_NAI_K6yz_Px = I_NAI_Lx6yz_S+ABX*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Px = I_NAI_Lx5y2z_S+ABX*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Px = I_NAI_Lx4y3z_S+ABX*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Px = I_NAI_Lx3y4z_S+ABX*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Px = I_NAI_Lx2y5z_S+ABX*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Px = I_NAI_Lxy6z_S+ABX*I_NAI_Ky6z_S;
  Double I_NAI_K6xy_Py = I_NAI_L6x2y_S+ABY*I_NAI_K6xy_S;
  Double I_NAI_K5x2y_Py = I_NAI_L5x3y_S+ABY*I_NAI_K5x2y_S;
  Double I_NAI_K5xyz_Py = I_NAI_L5x2yz_S+ABY*I_NAI_K5xyz_S;
  Double I_NAI_K4x3y_Py = I_NAI_L4x4y_S+ABY*I_NAI_K4x3y_S;
  Double I_NAI_K4x2yz_Py = I_NAI_L4x3yz_S+ABY*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Py = I_NAI_L4x2y2z_S+ABY*I_NAI_K4xy2z_S;
  Double I_NAI_K3x4y_Py = I_NAI_L3x5y_S+ABY*I_NAI_K3x4y_S;
  Double I_NAI_K3x3yz_Py = I_NAI_L3x4yz_S+ABY*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Py = I_NAI_L3x3y2z_S+ABY*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Py = I_NAI_L3x2y3z_S+ABY*I_NAI_K3xy3z_S;
  Double I_NAI_K2x5y_Py = I_NAI_L2x6y_S+ABY*I_NAI_K2x5y_S;
  Double I_NAI_K2x4yz_Py = I_NAI_L2x5yz_S+ABY*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Py = I_NAI_L2x4y2z_S+ABY*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Py = I_NAI_L2x3y3z_S+ABY*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Py = I_NAI_L2x2y4z_S+ABY*I_NAI_K2xy4z_S;
  Double I_NAI_Kx6y_Py = I_NAI_Lx7y_S+ABY*I_NAI_Kx6y_S;
  Double I_NAI_Kx5yz_Py = I_NAI_Lx6yz_S+ABY*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Py = I_NAI_Lx5y2z_S+ABY*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Py = I_NAI_Lx4y3z_S+ABY*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Py = I_NAI_Lx3y4z_S+ABY*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Py = I_NAI_Lx2y5z_S+ABY*I_NAI_Kxy5z_S;
  Double I_NAI_K7y_Py = I_NAI_L8y_S+ABY*I_NAI_K7y_S;
  Double I_NAI_K6yz_Py = I_NAI_L7yz_S+ABY*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Py = I_NAI_L6y2z_S+ABY*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Py = I_NAI_L5y3z_S+ABY*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Py = I_NAI_L4y4z_S+ABY*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Py = I_NAI_L3y5z_S+ABY*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Py = I_NAI_L2y6z_S+ABY*I_NAI_Ky6z_S;
  Double I_NAI_K6xz_Pz = I_NAI_L6x2z_S+ABZ*I_NAI_K6xz_S;
  Double I_NAI_K5xyz_Pz = I_NAI_L5xy2z_S+ABZ*I_NAI_K5xyz_S;
  Double I_NAI_K5x2z_Pz = I_NAI_L5x3z_S+ABZ*I_NAI_K5x2z_S;
  Double I_NAI_K4x2yz_Pz = I_NAI_L4x2y2z_S+ABZ*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Pz = I_NAI_L4xy3z_S+ABZ*I_NAI_K4xy2z_S;
  Double I_NAI_K4x3z_Pz = I_NAI_L4x4z_S+ABZ*I_NAI_K4x3z_S;
  Double I_NAI_K3x3yz_Pz = I_NAI_L3x3y2z_S+ABZ*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Pz = I_NAI_L3x2y3z_S+ABZ*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Pz = I_NAI_L3xy4z_S+ABZ*I_NAI_K3xy3z_S;
  Double I_NAI_K3x4z_Pz = I_NAI_L3x5z_S+ABZ*I_NAI_K3x4z_S;
  Double I_NAI_K2x4yz_Pz = I_NAI_L2x4y2z_S+ABZ*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Pz = I_NAI_L2x3y3z_S+ABZ*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Pz = I_NAI_L2x2y4z_S+ABZ*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Pz = I_NAI_L2xy5z_S+ABZ*I_NAI_K2xy4z_S;
  Double I_NAI_K2x5z_Pz = I_NAI_L2x6z_S+ABZ*I_NAI_K2x5z_S;
  Double I_NAI_Kx5yz_Pz = I_NAI_Lx5y2z_S+ABZ*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Pz = I_NAI_Lx4y3z_S+ABZ*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Pz = I_NAI_Lx3y4z_S+ABZ*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Pz = I_NAI_Lx2y5z_S+ABZ*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Pz = I_NAI_Lxy6z_S+ABZ*I_NAI_Kxy5z_S;
  Double I_NAI_Kx6z_Pz = I_NAI_Lx7z_S+ABZ*I_NAI_Kx6z_S;
  Double I_NAI_K6yz_Pz = I_NAI_L6y2z_S+ABZ*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Pz = I_NAI_L5y3z_S+ABZ*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Pz = I_NAI_L4y4z_S+ABZ*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Pz = I_NAI_L3y5z_S+ABZ*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Pz = I_NAI_L2y6z_S+ABZ*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Pz = I_NAI_Ly7z_S+ABZ*I_NAI_Ky6z_S;
  Double I_NAI_K7z_Pz = I_NAI_L8z_S+ABZ*I_NAI_K7z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P
   * RHS shell quartet name: SQ_NAI_I_P
   ************************************************************/
  Double I_NAI_I6x_D2x = I_NAI_K7x_Px+ABX*I_NAI_I6x_Px;
  Double I_NAI_I5xy_D2x = I_NAI_K6xy_Px+ABX*I_NAI_I5xy_Px;
  Double I_NAI_I5xz_D2x = I_NAI_K6xz_Px+ABX*I_NAI_I5xz_Px;
  Double I_NAI_I4x2y_D2x = I_NAI_K5x2y_Px+ABX*I_NAI_I4x2y_Px;
  Double I_NAI_I4xyz_D2x = I_NAI_K5xyz_Px+ABX*I_NAI_I4xyz_Px;
  Double I_NAI_I4x2z_D2x = I_NAI_K5x2z_Px+ABX*I_NAI_I4x2z_Px;
  Double I_NAI_I3x3y_D2x = I_NAI_K4x3y_Px+ABX*I_NAI_I3x3y_Px;
  Double I_NAI_I3x2yz_D2x = I_NAI_K4x2yz_Px+ABX*I_NAI_I3x2yz_Px;
  Double I_NAI_I3xy2z_D2x = I_NAI_K4xy2z_Px+ABX*I_NAI_I3xy2z_Px;
  Double I_NAI_I3x3z_D2x = I_NAI_K4x3z_Px+ABX*I_NAI_I3x3z_Px;
  Double I_NAI_I2x4y_D2x = I_NAI_K3x4y_Px+ABX*I_NAI_I2x4y_Px;
  Double I_NAI_I2x3yz_D2x = I_NAI_K3x3yz_Px+ABX*I_NAI_I2x3yz_Px;
  Double I_NAI_I2x2y2z_D2x = I_NAI_K3x2y2z_Px+ABX*I_NAI_I2x2y2z_Px;
  Double I_NAI_I2xy3z_D2x = I_NAI_K3xy3z_Px+ABX*I_NAI_I2xy3z_Px;
  Double I_NAI_I2x4z_D2x = I_NAI_K3x4z_Px+ABX*I_NAI_I2x4z_Px;
  Double I_NAI_Ix5y_D2x = I_NAI_K2x5y_Px+ABX*I_NAI_Ix5y_Px;
  Double I_NAI_Ix4yz_D2x = I_NAI_K2x4yz_Px+ABX*I_NAI_Ix4yz_Px;
  Double I_NAI_Ix3y2z_D2x = I_NAI_K2x3y2z_Px+ABX*I_NAI_Ix3y2z_Px;
  Double I_NAI_Ix2y3z_D2x = I_NAI_K2x2y3z_Px+ABX*I_NAI_Ix2y3z_Px;
  Double I_NAI_Ixy4z_D2x = I_NAI_K2xy4z_Px+ABX*I_NAI_Ixy4z_Px;
  Double I_NAI_Ix5z_D2x = I_NAI_K2x5z_Px+ABX*I_NAI_Ix5z_Px;
  Double I_NAI_I6y_D2x = I_NAI_Kx6y_Px+ABX*I_NAI_I6y_Px;
  Double I_NAI_I5yz_D2x = I_NAI_Kx5yz_Px+ABX*I_NAI_I5yz_Px;
  Double I_NAI_I4y2z_D2x = I_NAI_Kx4y2z_Px+ABX*I_NAI_I4y2z_Px;
  Double I_NAI_I3y3z_D2x = I_NAI_Kx3y3z_Px+ABX*I_NAI_I3y3z_Px;
  Double I_NAI_I2y4z_D2x = I_NAI_Kx2y4z_Px+ABX*I_NAI_I2y4z_Px;
  Double I_NAI_Iy5z_D2x = I_NAI_Kxy5z_Px+ABX*I_NAI_Iy5z_Px;
  Double I_NAI_I6z_D2x = I_NAI_Kx6z_Px+ABX*I_NAI_I6z_Px;
  Double I_NAI_I5xz_Dxy = I_NAI_K5xyz_Px+ABY*I_NAI_I5xz_Px;
  Double I_NAI_I4xyz_Dxy = I_NAI_K4x2yz_Px+ABY*I_NAI_I4xyz_Px;
  Double I_NAI_I4x2z_Dxy = I_NAI_K4xy2z_Px+ABY*I_NAI_I4x2z_Px;
  Double I_NAI_I3x2yz_Dxy = I_NAI_K3x3yz_Px+ABY*I_NAI_I3x2yz_Px;
  Double I_NAI_I3xy2z_Dxy = I_NAI_K3x2y2z_Px+ABY*I_NAI_I3xy2z_Px;
  Double I_NAI_I3x3z_Dxy = I_NAI_K3xy3z_Px+ABY*I_NAI_I3x3z_Px;
  Double I_NAI_I2x3yz_Dxy = I_NAI_K2x4yz_Px+ABY*I_NAI_I2x3yz_Px;
  Double I_NAI_I2x2y2z_Dxy = I_NAI_K2x3y2z_Px+ABY*I_NAI_I2x2y2z_Px;
  Double I_NAI_I2xy3z_Dxy = I_NAI_K2x2y3z_Px+ABY*I_NAI_I2xy3z_Px;
  Double I_NAI_I2x4z_Dxy = I_NAI_K2xy4z_Px+ABY*I_NAI_I2x4z_Px;
  Double I_NAI_Ix4yz_Dxy = I_NAI_Kx5yz_Px+ABY*I_NAI_Ix4yz_Px;
  Double I_NAI_Ix3y2z_Dxy = I_NAI_Kx4y2z_Px+ABY*I_NAI_Ix3y2z_Px;
  Double I_NAI_Ix2y3z_Dxy = I_NAI_Kx3y3z_Px+ABY*I_NAI_Ix2y3z_Px;
  Double I_NAI_Ixy4z_Dxy = I_NAI_Kx2y4z_Px+ABY*I_NAI_Ixy4z_Px;
  Double I_NAI_Ix5z_Dxy = I_NAI_Kxy5z_Px+ABY*I_NAI_Ix5z_Px;
  Double I_NAI_I5yz_Dxy = I_NAI_K6yz_Px+ABY*I_NAI_I5yz_Px;
  Double I_NAI_I4y2z_Dxy = I_NAI_K5y2z_Px+ABY*I_NAI_I4y2z_Px;
  Double I_NAI_I3y3z_Dxy = I_NAI_K4y3z_Px+ABY*I_NAI_I3y3z_Px;
  Double I_NAI_I2y4z_Dxy = I_NAI_K3y4z_Px+ABY*I_NAI_I2y4z_Px;
  Double I_NAI_Iy5z_Dxy = I_NAI_K2y5z_Px+ABY*I_NAI_Iy5z_Px;
  Double I_NAI_I6z_Dxy = I_NAI_Ky6z_Px+ABY*I_NAI_I6z_Px;
  Double I_NAI_I6x_D2y = I_NAI_K6xy_Py+ABY*I_NAI_I6x_Py;
  Double I_NAI_I5xy_D2y = I_NAI_K5x2y_Py+ABY*I_NAI_I5xy_Py;
  Double I_NAI_I5xz_D2y = I_NAI_K5xyz_Py+ABY*I_NAI_I5xz_Py;
  Double I_NAI_I4x2y_D2y = I_NAI_K4x3y_Py+ABY*I_NAI_I4x2y_Py;
  Double I_NAI_I4xyz_D2y = I_NAI_K4x2yz_Py+ABY*I_NAI_I4xyz_Py;
  Double I_NAI_I4x2z_D2y = I_NAI_K4xy2z_Py+ABY*I_NAI_I4x2z_Py;
  Double I_NAI_I3x3y_D2y = I_NAI_K3x4y_Py+ABY*I_NAI_I3x3y_Py;
  Double I_NAI_I3x2yz_D2y = I_NAI_K3x3yz_Py+ABY*I_NAI_I3x2yz_Py;
  Double I_NAI_I3xy2z_D2y = I_NAI_K3x2y2z_Py+ABY*I_NAI_I3xy2z_Py;
  Double I_NAI_I3x3z_D2y = I_NAI_K3xy3z_Py+ABY*I_NAI_I3x3z_Py;
  Double I_NAI_I2x4y_D2y = I_NAI_K2x5y_Py+ABY*I_NAI_I2x4y_Py;
  Double I_NAI_I2x3yz_D2y = I_NAI_K2x4yz_Py+ABY*I_NAI_I2x3yz_Py;
  Double I_NAI_I2x2y2z_D2y = I_NAI_K2x3y2z_Py+ABY*I_NAI_I2x2y2z_Py;
  Double I_NAI_I2xy3z_D2y = I_NAI_K2x2y3z_Py+ABY*I_NAI_I2xy3z_Py;
  Double I_NAI_I2x4z_D2y = I_NAI_K2xy4z_Py+ABY*I_NAI_I2x4z_Py;
  Double I_NAI_Ix5y_D2y = I_NAI_Kx6y_Py+ABY*I_NAI_Ix5y_Py;
  Double I_NAI_Ix4yz_D2y = I_NAI_Kx5yz_Py+ABY*I_NAI_Ix4yz_Py;
  Double I_NAI_Ix3y2z_D2y = I_NAI_Kx4y2z_Py+ABY*I_NAI_Ix3y2z_Py;
  Double I_NAI_Ix2y3z_D2y = I_NAI_Kx3y3z_Py+ABY*I_NAI_Ix2y3z_Py;
  Double I_NAI_Ixy4z_D2y = I_NAI_Kx2y4z_Py+ABY*I_NAI_Ixy4z_Py;
  Double I_NAI_Ix5z_D2y = I_NAI_Kxy5z_Py+ABY*I_NAI_Ix5z_Py;
  Double I_NAI_I6y_D2y = I_NAI_K7y_Py+ABY*I_NAI_I6y_Py;
  Double I_NAI_I5yz_D2y = I_NAI_K6yz_Py+ABY*I_NAI_I5yz_Py;
  Double I_NAI_I4y2z_D2y = I_NAI_K5y2z_Py+ABY*I_NAI_I4y2z_Py;
  Double I_NAI_I3y3z_D2y = I_NAI_K4y3z_Py+ABY*I_NAI_I3y3z_Py;
  Double I_NAI_I2y4z_D2y = I_NAI_K3y4z_Py+ABY*I_NAI_I2y4z_Py;
  Double I_NAI_Iy5z_D2y = I_NAI_K2y5z_Py+ABY*I_NAI_Iy5z_Py;
  Double I_NAI_I6z_D2y = I_NAI_Ky6z_Py+ABY*I_NAI_I6z_Py;
  Double I_NAI_I6x_D2z = I_NAI_K6xz_Pz+ABZ*I_NAI_I6x_Pz;
  Double I_NAI_I5xy_D2z = I_NAI_K5xyz_Pz+ABZ*I_NAI_I5xy_Pz;
  Double I_NAI_I5xz_D2z = I_NAI_K5x2z_Pz+ABZ*I_NAI_I5xz_Pz;
  Double I_NAI_I4x2y_D2z = I_NAI_K4x2yz_Pz+ABZ*I_NAI_I4x2y_Pz;
  Double I_NAI_I4xyz_D2z = I_NAI_K4xy2z_Pz+ABZ*I_NAI_I4xyz_Pz;
  Double I_NAI_I4x2z_D2z = I_NAI_K4x3z_Pz+ABZ*I_NAI_I4x2z_Pz;
  Double I_NAI_I3x3y_D2z = I_NAI_K3x3yz_Pz+ABZ*I_NAI_I3x3y_Pz;
  Double I_NAI_I3x2yz_D2z = I_NAI_K3x2y2z_Pz+ABZ*I_NAI_I3x2yz_Pz;
  Double I_NAI_I3xy2z_D2z = I_NAI_K3xy3z_Pz+ABZ*I_NAI_I3xy2z_Pz;
  Double I_NAI_I3x3z_D2z = I_NAI_K3x4z_Pz+ABZ*I_NAI_I3x3z_Pz;
  Double I_NAI_I2x4y_D2z = I_NAI_K2x4yz_Pz+ABZ*I_NAI_I2x4y_Pz;
  Double I_NAI_I2x3yz_D2z = I_NAI_K2x3y2z_Pz+ABZ*I_NAI_I2x3yz_Pz;
  Double I_NAI_I2x2y2z_D2z = I_NAI_K2x2y3z_Pz+ABZ*I_NAI_I2x2y2z_Pz;
  Double I_NAI_I2xy3z_D2z = I_NAI_K2xy4z_Pz+ABZ*I_NAI_I2xy3z_Pz;
  Double I_NAI_I2x4z_D2z = I_NAI_K2x5z_Pz+ABZ*I_NAI_I2x4z_Pz;
  Double I_NAI_Ix5y_D2z = I_NAI_Kx5yz_Pz+ABZ*I_NAI_Ix5y_Pz;
  Double I_NAI_Ix4yz_D2z = I_NAI_Kx4y2z_Pz+ABZ*I_NAI_Ix4yz_Pz;
  Double I_NAI_Ix3y2z_D2z = I_NAI_Kx3y3z_Pz+ABZ*I_NAI_Ix3y2z_Pz;
  Double I_NAI_Ix2y3z_D2z = I_NAI_Kx2y4z_Pz+ABZ*I_NAI_Ix2y3z_Pz;
  Double I_NAI_Ixy4z_D2z = I_NAI_Kxy5z_Pz+ABZ*I_NAI_Ixy4z_Pz;
  Double I_NAI_Ix5z_D2z = I_NAI_Kx6z_Pz+ABZ*I_NAI_Ix5z_Pz;
  Double I_NAI_I6y_D2z = I_NAI_K6yz_Pz+ABZ*I_NAI_I6y_Pz;
  Double I_NAI_I5yz_D2z = I_NAI_K5y2z_Pz+ABZ*I_NAI_I5yz_Pz;
  Double I_NAI_I4y2z_D2z = I_NAI_K4y3z_Pz+ABZ*I_NAI_I4y2z_Pz;
  Double I_NAI_I3y3z_D2z = I_NAI_K3y4z_Pz+ABZ*I_NAI_I3y3z_Pz;
  Double I_NAI_I2y4z_D2z = I_NAI_K2y5z_Pz+ABZ*I_NAI_I2y4z_Pz;
  Double I_NAI_Iy5z_D2z = I_NAI_Ky6z_Pz+ABZ*I_NAI_Iy5z_Pz;
  Double I_NAI_I6z_D2z = I_NAI_K7z_Pz+ABZ*I_NAI_I6z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D
   * RHS shell quartet name: SQ_NAI_H_D
   ************************************************************/
  abcd[0] = I_NAI_I6x_D2x+ABX*I_NAI_H5x_D2x;
  abcd[1] = I_NAI_I5xy_D2x+ABX*I_NAI_H4xy_D2x;
  abcd[2] = I_NAI_I5xz_D2x+ABX*I_NAI_H4xz_D2x;
  abcd[3] = I_NAI_I4x2y_D2x+ABX*I_NAI_H3x2y_D2x;
  abcd[4] = I_NAI_I4xyz_D2x+ABX*I_NAI_H3xyz_D2x;
  abcd[5] = I_NAI_I4x2z_D2x+ABX*I_NAI_H3x2z_D2x;
  abcd[6] = I_NAI_I3x3y_D2x+ABX*I_NAI_H2x3y_D2x;
  abcd[7] = I_NAI_I3x2yz_D2x+ABX*I_NAI_H2x2yz_D2x;
  abcd[8] = I_NAI_I3xy2z_D2x+ABX*I_NAI_H2xy2z_D2x;
  abcd[9] = I_NAI_I3x3z_D2x+ABX*I_NAI_H2x3z_D2x;
  abcd[10] = I_NAI_I2x4y_D2x+ABX*I_NAI_Hx4y_D2x;
  abcd[11] = I_NAI_I2x3yz_D2x+ABX*I_NAI_Hx3yz_D2x;
  abcd[12] = I_NAI_I2x2y2z_D2x+ABX*I_NAI_Hx2y2z_D2x;
  abcd[13] = I_NAI_I2xy3z_D2x+ABX*I_NAI_Hxy3z_D2x;
  abcd[14] = I_NAI_I2x4z_D2x+ABX*I_NAI_Hx4z_D2x;
  abcd[15] = I_NAI_Ix5y_D2x+ABX*I_NAI_H5y_D2x;
  abcd[16] = I_NAI_Ix4yz_D2x+ABX*I_NAI_H4yz_D2x;
  abcd[17] = I_NAI_Ix3y2z_D2x+ABX*I_NAI_H3y2z_D2x;
  abcd[18] = I_NAI_Ix2y3z_D2x+ABX*I_NAI_H2y3z_D2x;
  abcd[19] = I_NAI_Ixy4z_D2x+ABX*I_NAI_Hy4z_D2x;
  abcd[20] = I_NAI_Ix5z_D2x+ABX*I_NAI_H5z_D2x;
  abcd[21] = I_NAI_I5xy_D2x+ABY*I_NAI_H5x_D2x;
  abcd[22] = I_NAI_I4x2y_D2x+ABY*I_NAI_H4xy_D2x;
  abcd[23] = I_NAI_I4xyz_D2x+ABY*I_NAI_H4xz_D2x;
  abcd[24] = I_NAI_I3x3y_D2x+ABY*I_NAI_H3x2y_D2x;
  abcd[25] = I_NAI_I3x2yz_D2x+ABY*I_NAI_H3xyz_D2x;
  abcd[26] = I_NAI_I3xy2z_D2x+ABY*I_NAI_H3x2z_D2x;
  abcd[27] = I_NAI_I2x4y_D2x+ABY*I_NAI_H2x3y_D2x;
  abcd[28] = I_NAI_I2x3yz_D2x+ABY*I_NAI_H2x2yz_D2x;
  abcd[29] = I_NAI_I2x2y2z_D2x+ABY*I_NAI_H2xy2z_D2x;
  abcd[30] = I_NAI_I2xy3z_D2x+ABY*I_NAI_H2x3z_D2x;
  abcd[31] = I_NAI_Ix5y_D2x+ABY*I_NAI_Hx4y_D2x;
  abcd[32] = I_NAI_Ix4yz_D2x+ABY*I_NAI_Hx3yz_D2x;
  abcd[33] = I_NAI_Ix3y2z_D2x+ABY*I_NAI_Hx2y2z_D2x;
  abcd[34] = I_NAI_Ix2y3z_D2x+ABY*I_NAI_Hxy3z_D2x;
  abcd[35] = I_NAI_Ixy4z_D2x+ABY*I_NAI_Hx4z_D2x;
  abcd[36] = I_NAI_I6y_D2x+ABY*I_NAI_H5y_D2x;
  abcd[37] = I_NAI_I5yz_D2x+ABY*I_NAI_H4yz_D2x;
  abcd[38] = I_NAI_I4y2z_D2x+ABY*I_NAI_H3y2z_D2x;
  abcd[39] = I_NAI_I3y3z_D2x+ABY*I_NAI_H2y3z_D2x;
  abcd[40] = I_NAI_I2y4z_D2x+ABY*I_NAI_Hy4z_D2x;
  abcd[41] = I_NAI_Iy5z_D2x+ABY*I_NAI_H5z_D2x;
  abcd[42] = I_NAI_I5xz_D2x+ABZ*I_NAI_H5x_D2x;
  abcd[43] = I_NAI_I4xyz_D2x+ABZ*I_NAI_H4xy_D2x;
  abcd[44] = I_NAI_I4x2z_D2x+ABZ*I_NAI_H4xz_D2x;
  abcd[45] = I_NAI_I3x2yz_D2x+ABZ*I_NAI_H3x2y_D2x;
  abcd[46] = I_NAI_I3xy2z_D2x+ABZ*I_NAI_H3xyz_D2x;
  abcd[47] = I_NAI_I3x3z_D2x+ABZ*I_NAI_H3x2z_D2x;
  abcd[48] = I_NAI_I2x3yz_D2x+ABZ*I_NAI_H2x3y_D2x;
  abcd[49] = I_NAI_I2x2y2z_D2x+ABZ*I_NAI_H2x2yz_D2x;
  abcd[50] = I_NAI_I2xy3z_D2x+ABZ*I_NAI_H2xy2z_D2x;
  abcd[51] = I_NAI_I2x4z_D2x+ABZ*I_NAI_H2x3z_D2x;
  abcd[52] = I_NAI_Ix4yz_D2x+ABZ*I_NAI_Hx4y_D2x;
  abcd[53] = I_NAI_Ix3y2z_D2x+ABZ*I_NAI_Hx3yz_D2x;
  abcd[54] = I_NAI_Ix2y3z_D2x+ABZ*I_NAI_Hx2y2z_D2x;
  abcd[55] = I_NAI_Ixy4z_D2x+ABZ*I_NAI_Hxy3z_D2x;
  abcd[56] = I_NAI_Ix5z_D2x+ABZ*I_NAI_Hx4z_D2x;
  abcd[57] = I_NAI_I5yz_D2x+ABZ*I_NAI_H5y_D2x;
  abcd[58] = I_NAI_I4y2z_D2x+ABZ*I_NAI_H4yz_D2x;
  abcd[59] = I_NAI_I3y3z_D2x+ABZ*I_NAI_H3y2z_D2x;
  abcd[60] = I_NAI_I2y4z_D2x+ABZ*I_NAI_H2y3z_D2x;
  abcd[61] = I_NAI_Iy5z_D2x+ABZ*I_NAI_Hy4z_D2x;
  abcd[62] = I_NAI_I6z_D2x+ABZ*I_NAI_H5z_D2x;
  abcd[63] = I_NAI_I6x_D2y+ABX*I_NAI_H5x_D2y;
  abcd[64] = I_NAI_I5xy_D2y+ABX*I_NAI_H4xy_D2y;
  abcd[65] = I_NAI_I5xz_D2y+ABX*I_NAI_H4xz_D2y;
  abcd[66] = I_NAI_I4x2y_D2y+ABX*I_NAI_H3x2y_D2y;
  abcd[67] = I_NAI_I4xyz_D2y+ABX*I_NAI_H3xyz_D2y;
  abcd[68] = I_NAI_I4x2z_D2y+ABX*I_NAI_H3x2z_D2y;
  abcd[69] = I_NAI_I3x3y_D2y+ABX*I_NAI_H2x3y_D2y;
  abcd[70] = I_NAI_I3x2yz_D2y+ABX*I_NAI_H2x2yz_D2y;
  abcd[71] = I_NAI_I3xy2z_D2y+ABX*I_NAI_H2xy2z_D2y;
  abcd[72] = I_NAI_I3x3z_D2y+ABX*I_NAI_H2x3z_D2y;
  abcd[73] = I_NAI_I2x4y_D2y+ABX*I_NAI_Hx4y_D2y;
  abcd[74] = I_NAI_I2x3yz_D2y+ABX*I_NAI_Hx3yz_D2y;
  abcd[75] = I_NAI_I2x2y2z_D2y+ABX*I_NAI_Hx2y2z_D2y;
  abcd[76] = I_NAI_I2xy3z_D2y+ABX*I_NAI_Hxy3z_D2y;
  abcd[77] = I_NAI_I2x4z_D2y+ABX*I_NAI_Hx4z_D2y;
  abcd[78] = I_NAI_Ix5y_D2y+ABX*I_NAI_H5y_D2y;
  abcd[79] = I_NAI_Ix4yz_D2y+ABX*I_NAI_H4yz_D2y;
  abcd[80] = I_NAI_Ix3y2z_D2y+ABX*I_NAI_H3y2z_D2y;
  abcd[81] = I_NAI_Ix2y3z_D2y+ABX*I_NAI_H2y3z_D2y;
  abcd[82] = I_NAI_Ixy4z_D2y+ABX*I_NAI_Hy4z_D2y;
  abcd[83] = I_NAI_Ix5z_D2y+ABX*I_NAI_H5z_D2y;
  abcd[84] = I_NAI_I5xz_Dxy+ABZ*I_NAI_H5x_Dxy;
  abcd[85] = I_NAI_I4xyz_Dxy+ABZ*I_NAI_H4xy_Dxy;
  abcd[86] = I_NAI_I4x2z_Dxy+ABZ*I_NAI_H4xz_Dxy;
  abcd[87] = I_NAI_I3x2yz_Dxy+ABZ*I_NAI_H3x2y_Dxy;
  abcd[88] = I_NAI_I3xy2z_Dxy+ABZ*I_NAI_H3xyz_Dxy;
  abcd[89] = I_NAI_I3x3z_Dxy+ABZ*I_NAI_H3x2z_Dxy;
  abcd[90] = I_NAI_I2x3yz_Dxy+ABZ*I_NAI_H2x3y_Dxy;
  abcd[91] = I_NAI_I2x2y2z_Dxy+ABZ*I_NAI_H2x2yz_Dxy;
  abcd[92] = I_NAI_I2xy3z_Dxy+ABZ*I_NAI_H2xy2z_Dxy;
  abcd[93] = I_NAI_I2x4z_Dxy+ABZ*I_NAI_H2x3z_Dxy;
  abcd[94] = I_NAI_Ix4yz_Dxy+ABZ*I_NAI_Hx4y_Dxy;
  abcd[95] = I_NAI_Ix3y2z_Dxy+ABZ*I_NAI_Hx3yz_Dxy;
  abcd[96] = I_NAI_Ix2y3z_Dxy+ABZ*I_NAI_Hx2y2z_Dxy;
  abcd[97] = I_NAI_Ixy4z_Dxy+ABZ*I_NAI_Hxy3z_Dxy;
  abcd[98] = I_NAI_Ix5z_Dxy+ABZ*I_NAI_Hx4z_Dxy;
  abcd[99] = I_NAI_I5yz_Dxy+ABZ*I_NAI_H5y_Dxy;
  abcd[100] = I_NAI_I4y2z_Dxy+ABZ*I_NAI_H4yz_Dxy;
  abcd[101] = I_NAI_I3y3z_Dxy+ABZ*I_NAI_H3y2z_Dxy;
  abcd[102] = I_NAI_I2y4z_Dxy+ABZ*I_NAI_H2y3z_Dxy;
  abcd[103] = I_NAI_Iy5z_Dxy+ABZ*I_NAI_Hy4z_Dxy;
  abcd[104] = I_NAI_I6z_Dxy+ABZ*I_NAI_H5z_Dxy;
  abcd[105] = I_NAI_I6x_D2z+ABX*I_NAI_H5x_D2z;
  abcd[106] = I_NAI_I5xy_D2z+ABX*I_NAI_H4xy_D2z;
  abcd[107] = I_NAI_I5xz_D2z+ABX*I_NAI_H4xz_D2z;
  abcd[108] = I_NAI_I4x2y_D2z+ABX*I_NAI_H3x2y_D2z;
  abcd[109] = I_NAI_I4xyz_D2z+ABX*I_NAI_H3xyz_D2z;
  abcd[110] = I_NAI_I4x2z_D2z+ABX*I_NAI_H3x2z_D2z;
  abcd[111] = I_NAI_I3x3y_D2z+ABX*I_NAI_H2x3y_D2z;
  abcd[112] = I_NAI_I3x2yz_D2z+ABX*I_NAI_H2x2yz_D2z;
  abcd[113] = I_NAI_I3xy2z_D2z+ABX*I_NAI_H2xy2z_D2z;
  abcd[114] = I_NAI_I3x3z_D2z+ABX*I_NAI_H2x3z_D2z;
  abcd[115] = I_NAI_I2x4y_D2z+ABX*I_NAI_Hx4y_D2z;
  abcd[116] = I_NAI_I2x3yz_D2z+ABX*I_NAI_Hx3yz_D2z;
  abcd[117] = I_NAI_I2x2y2z_D2z+ABX*I_NAI_Hx2y2z_D2z;
  abcd[118] = I_NAI_I2xy3z_D2z+ABX*I_NAI_Hxy3z_D2z;
  abcd[119] = I_NAI_I2x4z_D2z+ABX*I_NAI_Hx4z_D2z;
  abcd[120] = I_NAI_Ix5y_D2z+ABX*I_NAI_H5y_D2z;
  abcd[121] = I_NAI_Ix4yz_D2z+ABX*I_NAI_H4yz_D2z;
  abcd[122] = I_NAI_Ix3y2z_D2z+ABX*I_NAI_H3y2z_D2z;
  abcd[123] = I_NAI_Ix2y3z_D2z+ABX*I_NAI_H2y3z_D2z;
  abcd[124] = I_NAI_Ixy4z_D2z+ABX*I_NAI_Hy4z_D2z;
  abcd[125] = I_NAI_Ix5z_D2z+ABX*I_NAI_H5z_D2z;
  abcd[126] = I_NAI_I5xy_D2y+ABY*I_NAI_H5x_D2y;
  abcd[127] = I_NAI_I4x2y_D2y+ABY*I_NAI_H4xy_D2y;
  abcd[128] = I_NAI_I4xyz_D2y+ABY*I_NAI_H4xz_D2y;
  abcd[129] = I_NAI_I3x3y_D2y+ABY*I_NAI_H3x2y_D2y;
  abcd[130] = I_NAI_I3x2yz_D2y+ABY*I_NAI_H3xyz_D2y;
  abcd[131] = I_NAI_I3xy2z_D2y+ABY*I_NAI_H3x2z_D2y;
  abcd[132] = I_NAI_I2x4y_D2y+ABY*I_NAI_H2x3y_D2y;
  abcd[133] = I_NAI_I2x3yz_D2y+ABY*I_NAI_H2x2yz_D2y;
  abcd[134] = I_NAI_I2x2y2z_D2y+ABY*I_NAI_H2xy2z_D2y;
  abcd[135] = I_NAI_I2xy3z_D2y+ABY*I_NAI_H2x3z_D2y;
  abcd[136] = I_NAI_Ix5y_D2y+ABY*I_NAI_Hx4y_D2y;
  abcd[137] = I_NAI_Ix4yz_D2y+ABY*I_NAI_Hx3yz_D2y;
  abcd[138] = I_NAI_Ix3y2z_D2y+ABY*I_NAI_Hx2y2z_D2y;
  abcd[139] = I_NAI_Ix2y3z_D2y+ABY*I_NAI_Hxy3z_D2y;
  abcd[140] = I_NAI_Ixy4z_D2y+ABY*I_NAI_Hx4z_D2y;
  abcd[141] = I_NAI_I6y_D2y+ABY*I_NAI_H5y_D2y;
  abcd[142] = I_NAI_I5yz_D2y+ABY*I_NAI_H4yz_D2y;
  abcd[143] = I_NAI_I4y2z_D2y+ABY*I_NAI_H3y2z_D2y;
  abcd[144] = I_NAI_I3y3z_D2y+ABY*I_NAI_H2y3z_D2y;
  abcd[145] = I_NAI_I2y4z_D2y+ABY*I_NAI_Hy4z_D2y;
  abcd[146] = I_NAI_Iy5z_D2y+ABY*I_NAI_H5z_D2y;
  abcd[147] = I_NAI_I5xz_D2y+ABZ*I_NAI_H5x_D2y;
  abcd[148] = I_NAI_I4xyz_D2y+ABZ*I_NAI_H4xy_D2y;
  abcd[149] = I_NAI_I4x2z_D2y+ABZ*I_NAI_H4xz_D2y;
  abcd[150] = I_NAI_I3x2yz_D2y+ABZ*I_NAI_H3x2y_D2y;
  abcd[151] = I_NAI_I3xy2z_D2y+ABZ*I_NAI_H3xyz_D2y;
  abcd[152] = I_NAI_I3x3z_D2y+ABZ*I_NAI_H3x2z_D2y;
  abcd[153] = I_NAI_I2x3yz_D2y+ABZ*I_NAI_H2x3y_D2y;
  abcd[154] = I_NAI_I2x2y2z_D2y+ABZ*I_NAI_H2x2yz_D2y;
  abcd[155] = I_NAI_I2xy3z_D2y+ABZ*I_NAI_H2xy2z_D2y;
  abcd[156] = I_NAI_I2x4z_D2y+ABZ*I_NAI_H2x3z_D2y;
  abcd[157] = I_NAI_Ix4yz_D2y+ABZ*I_NAI_Hx4y_D2y;
  abcd[158] = I_NAI_Ix3y2z_D2y+ABZ*I_NAI_Hx3yz_D2y;
  abcd[159] = I_NAI_Ix2y3z_D2y+ABZ*I_NAI_Hx2y2z_D2y;
  abcd[160] = I_NAI_Ixy4z_D2y+ABZ*I_NAI_Hxy3z_D2y;
  abcd[161] = I_NAI_Ix5z_D2y+ABZ*I_NAI_Hx4z_D2y;
  abcd[162] = I_NAI_I5yz_D2y+ABZ*I_NAI_H5y_D2y;
  abcd[163] = I_NAI_I4y2z_D2y+ABZ*I_NAI_H4yz_D2y;
  abcd[164] = I_NAI_I3y3z_D2y+ABZ*I_NAI_H3y2z_D2y;
  abcd[165] = I_NAI_I2y4z_D2y+ABZ*I_NAI_H2y3z_D2y;
  abcd[166] = I_NAI_Iy5z_D2y+ABZ*I_NAI_Hy4z_D2y;
  abcd[167] = I_NAI_I6z_D2y+ABZ*I_NAI_H5z_D2y;
  abcd[168] = I_NAI_I5xy_D2z+ABY*I_NAI_H5x_D2z;
  abcd[169] = I_NAI_I4x2y_D2z+ABY*I_NAI_H4xy_D2z;
  abcd[170] = I_NAI_I4xyz_D2z+ABY*I_NAI_H4xz_D2z;
  abcd[171] = I_NAI_I3x3y_D2z+ABY*I_NAI_H3x2y_D2z;
  abcd[172] = I_NAI_I3x2yz_D2z+ABY*I_NAI_H3xyz_D2z;
  abcd[173] = I_NAI_I3xy2z_D2z+ABY*I_NAI_H3x2z_D2z;
  abcd[174] = I_NAI_I2x4y_D2z+ABY*I_NAI_H2x3y_D2z;
  abcd[175] = I_NAI_I2x3yz_D2z+ABY*I_NAI_H2x2yz_D2z;
  abcd[176] = I_NAI_I2x2y2z_D2z+ABY*I_NAI_H2xy2z_D2z;
  abcd[177] = I_NAI_I2xy3z_D2z+ABY*I_NAI_H2x3z_D2z;
  abcd[178] = I_NAI_Ix5y_D2z+ABY*I_NAI_Hx4y_D2z;
  abcd[179] = I_NAI_Ix4yz_D2z+ABY*I_NAI_Hx3yz_D2z;
  abcd[180] = I_NAI_Ix3y2z_D2z+ABY*I_NAI_Hx2y2z_D2z;
  abcd[181] = I_NAI_Ix2y3z_D2z+ABY*I_NAI_Hxy3z_D2z;
  abcd[182] = I_NAI_Ixy4z_D2z+ABY*I_NAI_Hx4z_D2z;
  abcd[183] = I_NAI_I6y_D2z+ABY*I_NAI_H5y_D2z;
  abcd[184] = I_NAI_I5yz_D2z+ABY*I_NAI_H4yz_D2z;
  abcd[185] = I_NAI_I4y2z_D2z+ABY*I_NAI_H3y2z_D2z;
  abcd[186] = I_NAI_I3y3z_D2z+ABY*I_NAI_H2y3z_D2z;
  abcd[187] = I_NAI_I2y4z_D2z+ABY*I_NAI_Hy4z_D2z;
  abcd[188] = I_NAI_Iy5z_D2z+ABY*I_NAI_H5z_D2z;
  abcd[189] = I_NAI_I5xz_D2z+ABZ*I_NAI_H5x_D2z;
  abcd[190] = I_NAI_I4xyz_D2z+ABZ*I_NAI_H4xy_D2z;
  abcd[191] = I_NAI_I4x2z_D2z+ABZ*I_NAI_H4xz_D2z;
  abcd[192] = I_NAI_I3x2yz_D2z+ABZ*I_NAI_H3x2y_D2z;
  abcd[193] = I_NAI_I3xy2z_D2z+ABZ*I_NAI_H3xyz_D2z;
  abcd[194] = I_NAI_I3x3z_D2z+ABZ*I_NAI_H3x2z_D2z;
  abcd[195] = I_NAI_I2x3yz_D2z+ABZ*I_NAI_H2x3y_D2z;
  abcd[196] = I_NAI_I2x2y2z_D2z+ABZ*I_NAI_H2x2yz_D2z;
  abcd[197] = I_NAI_I2xy3z_D2z+ABZ*I_NAI_H2xy2z_D2z;
  abcd[198] = I_NAI_I2x4z_D2z+ABZ*I_NAI_H2x3z_D2z;
  abcd[199] = I_NAI_Ix4yz_D2z+ABZ*I_NAI_Hx4y_D2z;
  abcd[200] = I_NAI_Ix3y2z_D2z+ABZ*I_NAI_Hx3yz_D2z;
  abcd[201] = I_NAI_Ix2y3z_D2z+ABZ*I_NAI_Hx2y2z_D2z;
  abcd[202] = I_NAI_Ixy4z_D2z+ABZ*I_NAI_Hxy3z_D2z;
  abcd[203] = I_NAI_Ix5z_D2z+ABZ*I_NAI_Hx4z_D2z;
  abcd[204] = I_NAI_I5yz_D2z+ABZ*I_NAI_H5y_D2z;
  abcd[205] = I_NAI_I4y2z_D2z+ABZ*I_NAI_H4yz_D2z;
  abcd[206] = I_NAI_I3y3z_D2z+ABZ*I_NAI_H3y2z_D2z;
  abcd[207] = I_NAI_I2y4z_D2z+ABZ*I_NAI_H2y3z_D2z;
  abcd[208] = I_NAI_Iy5z_D2z+ABZ*I_NAI_Hy4z_D2z;
  abcd[209] = I_NAI_I6z_D2z+ABZ*I_NAI_H5z_D2z;
}
