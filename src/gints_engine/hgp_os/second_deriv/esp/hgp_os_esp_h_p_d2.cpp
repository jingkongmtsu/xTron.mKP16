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
// BRA1 as redundant position, total RHS integrals evaluated as: 5430
// BRA2 as redundant position, total RHS integrals evaluated as: 4896
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_esp_h_p_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_L8x_S_aa = 0.0E0;
    Double I_ESP_L7xy_S_aa = 0.0E0;
    Double I_ESP_L7xz_S_aa = 0.0E0;
    Double I_ESP_L6x2y_S_aa = 0.0E0;
    Double I_ESP_L6xyz_S_aa = 0.0E0;
    Double I_ESP_L6x2z_S_aa = 0.0E0;
    Double I_ESP_L5x3y_S_aa = 0.0E0;
    Double I_ESP_L5x2yz_S_aa = 0.0E0;
    Double I_ESP_L5xy2z_S_aa = 0.0E0;
    Double I_ESP_L5x3z_S_aa = 0.0E0;
    Double I_ESP_L4x4y_S_aa = 0.0E0;
    Double I_ESP_L4x3yz_S_aa = 0.0E0;
    Double I_ESP_L4x2y2z_S_aa = 0.0E0;
    Double I_ESP_L4xy3z_S_aa = 0.0E0;
    Double I_ESP_L4x4z_S_aa = 0.0E0;
    Double I_ESP_L3x5y_S_aa = 0.0E0;
    Double I_ESP_L3x4yz_S_aa = 0.0E0;
    Double I_ESP_L3x3y2z_S_aa = 0.0E0;
    Double I_ESP_L3x2y3z_S_aa = 0.0E0;
    Double I_ESP_L3xy4z_S_aa = 0.0E0;
    Double I_ESP_L3x5z_S_aa = 0.0E0;
    Double I_ESP_L2x6y_S_aa = 0.0E0;
    Double I_ESP_L2x5yz_S_aa = 0.0E0;
    Double I_ESP_L2x4y2z_S_aa = 0.0E0;
    Double I_ESP_L2x3y3z_S_aa = 0.0E0;
    Double I_ESP_L2x2y4z_S_aa = 0.0E0;
    Double I_ESP_L2xy5z_S_aa = 0.0E0;
    Double I_ESP_L2x6z_S_aa = 0.0E0;
    Double I_ESP_Lx7y_S_aa = 0.0E0;
    Double I_ESP_Lx6yz_S_aa = 0.0E0;
    Double I_ESP_Lx5y2z_S_aa = 0.0E0;
    Double I_ESP_Lx4y3z_S_aa = 0.0E0;
    Double I_ESP_Lx3y4z_S_aa = 0.0E0;
    Double I_ESP_Lx2y5z_S_aa = 0.0E0;
    Double I_ESP_Lxy6z_S_aa = 0.0E0;
    Double I_ESP_Lx7z_S_aa = 0.0E0;
    Double I_ESP_L8y_S_aa = 0.0E0;
    Double I_ESP_L7yz_S_aa = 0.0E0;
    Double I_ESP_L6y2z_S_aa = 0.0E0;
    Double I_ESP_L5y3z_S_aa = 0.0E0;
    Double I_ESP_L4y4z_S_aa = 0.0E0;
    Double I_ESP_L3y5z_S_aa = 0.0E0;
    Double I_ESP_L2y6z_S_aa = 0.0E0;
    Double I_ESP_Ly7z_S_aa = 0.0E0;
    Double I_ESP_L8z_S_aa = 0.0E0;
    Double I_ESP_K7x_S_aa = 0.0E0;
    Double I_ESP_K6xy_S_aa = 0.0E0;
    Double I_ESP_K6xz_S_aa = 0.0E0;
    Double I_ESP_K5x2y_S_aa = 0.0E0;
    Double I_ESP_K5xyz_S_aa = 0.0E0;
    Double I_ESP_K5x2z_S_aa = 0.0E0;
    Double I_ESP_K4x3y_S_aa = 0.0E0;
    Double I_ESP_K4x2yz_S_aa = 0.0E0;
    Double I_ESP_K4xy2z_S_aa = 0.0E0;
    Double I_ESP_K4x3z_S_aa = 0.0E0;
    Double I_ESP_K3x4y_S_aa = 0.0E0;
    Double I_ESP_K3x3yz_S_aa = 0.0E0;
    Double I_ESP_K3x2y2z_S_aa = 0.0E0;
    Double I_ESP_K3xy3z_S_aa = 0.0E0;
    Double I_ESP_K3x4z_S_aa = 0.0E0;
    Double I_ESP_K2x5y_S_aa = 0.0E0;
    Double I_ESP_K2x4yz_S_aa = 0.0E0;
    Double I_ESP_K2x3y2z_S_aa = 0.0E0;
    Double I_ESP_K2x2y3z_S_aa = 0.0E0;
    Double I_ESP_K2xy4z_S_aa = 0.0E0;
    Double I_ESP_K2x5z_S_aa = 0.0E0;
    Double I_ESP_Kx6y_S_aa = 0.0E0;
    Double I_ESP_Kx5yz_S_aa = 0.0E0;
    Double I_ESP_Kx4y2z_S_aa = 0.0E0;
    Double I_ESP_Kx3y3z_S_aa = 0.0E0;
    Double I_ESP_Kx2y4z_S_aa = 0.0E0;
    Double I_ESP_Kxy5z_S_aa = 0.0E0;
    Double I_ESP_Kx6z_S_aa = 0.0E0;
    Double I_ESP_K7y_S_aa = 0.0E0;
    Double I_ESP_K6yz_S_aa = 0.0E0;
    Double I_ESP_K5y2z_S_aa = 0.0E0;
    Double I_ESP_K4y3z_S_aa = 0.0E0;
    Double I_ESP_K3y4z_S_aa = 0.0E0;
    Double I_ESP_K2y5z_S_aa = 0.0E0;
    Double I_ESP_Ky6z_S_aa = 0.0E0;
    Double I_ESP_K7z_S_aa = 0.0E0;
    Double I_ESP_I6x_S_a = 0.0E0;
    Double I_ESP_I5xy_S_a = 0.0E0;
    Double I_ESP_I5xz_S_a = 0.0E0;
    Double I_ESP_I4x2y_S_a = 0.0E0;
    Double I_ESP_I4xyz_S_a = 0.0E0;
    Double I_ESP_I4x2z_S_a = 0.0E0;
    Double I_ESP_I3x3y_S_a = 0.0E0;
    Double I_ESP_I3x2yz_S_a = 0.0E0;
    Double I_ESP_I3xy2z_S_a = 0.0E0;
    Double I_ESP_I3x3z_S_a = 0.0E0;
    Double I_ESP_I2x4y_S_a = 0.0E0;
    Double I_ESP_I2x3yz_S_a = 0.0E0;
    Double I_ESP_I2x2y2z_S_a = 0.0E0;
    Double I_ESP_I2xy3z_S_a = 0.0E0;
    Double I_ESP_I2x4z_S_a = 0.0E0;
    Double I_ESP_Ix5y_S_a = 0.0E0;
    Double I_ESP_Ix4yz_S_a = 0.0E0;
    Double I_ESP_Ix3y2z_S_a = 0.0E0;
    Double I_ESP_Ix2y3z_S_a = 0.0E0;
    Double I_ESP_Ixy4z_S_a = 0.0E0;
    Double I_ESP_Ix5z_S_a = 0.0E0;
    Double I_ESP_I6y_S_a = 0.0E0;
    Double I_ESP_I5yz_S_a = 0.0E0;
    Double I_ESP_I4y2z_S_a = 0.0E0;
    Double I_ESP_I3y3z_S_a = 0.0E0;
    Double I_ESP_I2y4z_S_a = 0.0E0;
    Double I_ESP_Iy5z_S_a = 0.0E0;
    Double I_ESP_I6z_S_a = 0.0E0;
    Double I_ESP_H5x_S_a = 0.0E0;
    Double I_ESP_H4xy_S_a = 0.0E0;
    Double I_ESP_H4xz_S_a = 0.0E0;
    Double I_ESP_H3x2y_S_a = 0.0E0;
    Double I_ESP_H3xyz_S_a = 0.0E0;
    Double I_ESP_H3x2z_S_a = 0.0E0;
    Double I_ESP_H2x3y_S_a = 0.0E0;
    Double I_ESP_H2x2yz_S_a = 0.0E0;
    Double I_ESP_H2xy2z_S_a = 0.0E0;
    Double I_ESP_H2x3z_S_a = 0.0E0;
    Double I_ESP_Hx4y_S_a = 0.0E0;
    Double I_ESP_Hx3yz_S_a = 0.0E0;
    Double I_ESP_Hx2y2z_S_a = 0.0E0;
    Double I_ESP_Hxy3z_S_a = 0.0E0;
    Double I_ESP_Hx4z_S_a = 0.0E0;
    Double I_ESP_H5y_S_a = 0.0E0;
    Double I_ESP_H4yz_S_a = 0.0E0;
    Double I_ESP_H3y2z_S_a = 0.0E0;
    Double I_ESP_H2y3z_S_a = 0.0E0;
    Double I_ESP_Hy4z_S_a = 0.0E0;
    Double I_ESP_H5z_S_a = 0.0E0;
    Double I_ESP_G4x_S = 0.0E0;
    Double I_ESP_G3xy_S = 0.0E0;
    Double I_ESP_G3xz_S = 0.0E0;
    Double I_ESP_G2x2y_S = 0.0E0;
    Double I_ESP_G2xyz_S = 0.0E0;
    Double I_ESP_G2x2z_S = 0.0E0;
    Double I_ESP_Gx3y_S = 0.0E0;
    Double I_ESP_Gx2yz_S = 0.0E0;
    Double I_ESP_Gxy2z_S = 0.0E0;
    Double I_ESP_Gx3z_S = 0.0E0;
    Double I_ESP_G4y_S = 0.0E0;
    Double I_ESP_G3yz_S = 0.0E0;
    Double I_ESP_G2y2z_S = 0.0E0;
    Double I_ESP_Gy3z_S = 0.0E0;
    Double I_ESP_G4z_S = 0.0E0;
    Double I_ESP_F3x_S = 0.0E0;
    Double I_ESP_F2xy_S = 0.0E0;
    Double I_ESP_F2xz_S = 0.0E0;
    Double I_ESP_Fx2y_S = 0.0E0;
    Double I_ESP_Fxyz_S = 0.0E0;
    Double I_ESP_Fx2z_S = 0.0E0;
    Double I_ESP_F3y_S = 0.0E0;
    Double I_ESP_F2yz_S = 0.0E0;
    Double I_ESP_Fy2z_S = 0.0E0;
    Double I_ESP_F3z_S = 0.0E0;

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
      Double onedz = iexp[ip2];
      Double rho   = 1.0E0/onedz;
      Double zeta  = rho;
      Double zdiff = iexpdiff[ip2];
      Double alpha = 0.5E0*(zeta+zdiff);
      Double beta  = 0.5E0*(zeta-zdiff);
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
      Double PRX   = PX - R[iGrid*3  ];
      Double PRY   = PY - R[iGrid*3+1];
      Double PRZ   = PZ - R[iGrid*3+2];
      Double PR2   = PRX*PRX+PRY*PRY+PRZ*PRZ;
      Double u     = rho*PR2;
      Double squ   = sqrt(u);
      Double prefactor = ic2*fbra;

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

      Double I_ESP_S_S_vrr  = 0.0E0;
      Double I_ESP_S_S_M1_vrr  = 0.0E0;
      Double I_ESP_S_S_M2_vrr  = 0.0E0;
      Double I_ESP_S_S_M3_vrr  = 0.0E0;
      Double I_ESP_S_S_M4_vrr  = 0.0E0;
      Double I_ESP_S_S_M5_vrr  = 0.0E0;
      Double I_ESP_S_S_M6_vrr  = 0.0E0;
      Double I_ESP_S_S_M7_vrr  = 0.0E0;
      Double I_ESP_S_S_M8_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ESP_S_S_vrr_d  = 0.0E0;
      double I_ESP_S_S_M1_vrr_d  = 0.0E0;
      double I_ESP_S_S_M2_vrr_d  = 0.0E0;
      double I_ESP_S_S_M3_vrr_d  = 0.0E0;
      double I_ESP_S_S_M4_vrr_d  = 0.0E0;
      double I_ESP_S_S_M5_vrr_d  = 0.0E0;
      double I_ESP_S_S_M6_vrr_d  = 0.0E0;
      double I_ESP_S_S_M7_vrr_d  = 0.0E0;
      double I_ESP_S_S_M8_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER51;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER49*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = 1.0E0+u2*ONEOVER19*I_ESP_S_S_M8_vrr;
        I_ESP_S_S_M8_vrr = ONEOVER17*I_ESP_S_S_M8_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M8_vrr  = f*I_ESP_S_S_M8_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M7_vrr  = ONEOVER15*(u2*I_ESP_S_S_M8_vrr+f);
        I_ESP_S_S_M6_vrr  = ONEOVER13*(u2*I_ESP_S_S_M7_vrr+f);
        I_ESP_S_S_M5_vrr  = ONEOVER11*(u2*I_ESP_S_S_M6_vrr+f);
        I_ESP_S_S_M4_vrr  = ONEOVER9*(u2*I_ESP_S_S_M5_vrr+f);
        I_ESP_S_S_M3_vrr  = ONEOVER7*(u2*I_ESP_S_S_M4_vrr+f);
        I_ESP_S_S_M2_vrr  = ONEOVER5*(u2*I_ESP_S_S_M3_vrr+f);
        I_ESP_S_S_M1_vrr  = ONEOVER3*(u2*I_ESP_S_S_M2_vrr+f);
        I_ESP_S_S_vrr  = ONEOVER1*(u2*I_ESP_S_S_M1_vrr+f);

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
          I_ESP_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ESP_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ESP_S_S_vrr_d-f);
        I_ESP_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr_d-f);
        I_ESP_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr_d-f);
        I_ESP_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr_d-f);
        I_ESP_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr_d-f);
        I_ESP_S_S_M6_vrr_d = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr_d-f);
        I_ESP_S_S_M7_vrr_d = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr_d-f);
        I_ESP_S_S_M8_vrr_d = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr_d-f);

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);
        I_ESP_S_S_M6_vrr = static_cast<Double>(I_ESP_S_S_M6_vrr_d);
        I_ESP_S_S_M7_vrr = static_cast<Double>(I_ESP_S_S_M7_vrr_d);
        I_ESP_S_S_M8_vrr = static_cast<Double>(I_ESP_S_S_M8_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ESP_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M1_vrr = oneO2u*(1.0E0*I_ESP_S_S_vrr-f);
        I_ESP_S_S_M2_vrr = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr-f);
        I_ESP_S_S_M3_vrr = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr-f);
        I_ESP_S_S_M4_vrr = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr-f);
        I_ESP_S_S_M5_vrr = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr-f);
        I_ESP_S_S_M6_vrr = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr-f);
        I_ESP_S_S_M7_vrr = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr-f);
        I_ESP_S_S_M8_vrr = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M8
         ************************************************************/
        Double I_ESP_Px_S_M7_vrr = PAX*I_ESP_S_S_M7_vrr-PRX*I_ESP_S_S_M8_vrr;
        Double I_ESP_Py_S_M7_vrr = PAY*I_ESP_S_S_M7_vrr-PRY*I_ESP_S_S_M8_vrr;
        Double I_ESP_Pz_S_M7_vrr = PAZ*I_ESP_S_S_M7_vrr-PRZ*I_ESP_S_S_M8_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_Px_S_M6_vrr = PAX*I_ESP_S_S_M6_vrr-PRX*I_ESP_S_S_M7_vrr;
        Double I_ESP_Py_S_M6_vrr = PAY*I_ESP_S_S_M6_vrr-PRY*I_ESP_S_S_M7_vrr;
        Double I_ESP_Pz_S_M6_vrr = PAZ*I_ESP_S_S_M6_vrr-PRZ*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_D2x_S_M6_vrr = PAX*I_ESP_Px_S_M6_vrr-PRX*I_ESP_Px_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2y_S_M6_vrr = PAY*I_ESP_Py_S_M6_vrr-PRY*I_ESP_Py_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2z_S_M6_vrr = PAZ*I_ESP_Pz_S_M6_vrr-PRZ*I_ESP_Pz_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_Px_S_M5_vrr = PAX*I_ESP_S_S_M5_vrr-PRX*I_ESP_S_S_M6_vrr;
        Double I_ESP_Py_S_M5_vrr = PAY*I_ESP_S_S_M5_vrr-PRY*I_ESP_S_S_M6_vrr;
        Double I_ESP_Pz_S_M5_vrr = PAZ*I_ESP_S_S_M5_vrr-PRZ*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_D2x_S_M5_vrr = PAX*I_ESP_Px_S_M5_vrr-PRX*I_ESP_Px_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2y_S_M5_vrr = PAY*I_ESP_Py_S_M5_vrr-PRY*I_ESP_Py_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2z_S_M5_vrr = PAZ*I_ESP_Pz_S_M5_vrr-PRZ*I_ESP_Pz_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         ************************************************************/
        Double I_ESP_F3x_S_M5_vrr = PAX*I_ESP_D2x_S_M5_vrr-PRX*I_ESP_D2x_S_M6_vrr+2*oned2z*I_ESP_Px_S_M5_vrr-2*oned2z*I_ESP_Px_S_M6_vrr;
        Double I_ESP_F3y_S_M5_vrr = PAY*I_ESP_D2y_S_M5_vrr-PRY*I_ESP_D2y_S_M6_vrr+2*oned2z*I_ESP_Py_S_M5_vrr-2*oned2z*I_ESP_Py_S_M6_vrr;
        Double I_ESP_F3z_S_M5_vrr = PAZ*I_ESP_D2z_S_M5_vrr-PRZ*I_ESP_D2z_S_M6_vrr+2*oned2z*I_ESP_Pz_S_M5_vrr-2*oned2z*I_ESP_Pz_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_Px_S_M4_vrr = PAX*I_ESP_S_S_M4_vrr-PRX*I_ESP_S_S_M5_vrr;
        Double I_ESP_Py_S_M4_vrr = PAY*I_ESP_S_S_M4_vrr-PRY*I_ESP_S_S_M5_vrr;
        Double I_ESP_Pz_S_M4_vrr = PAZ*I_ESP_S_S_M4_vrr-PRZ*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_D2x_S_M4_vrr = PAX*I_ESP_Px_S_M4_vrr-PRX*I_ESP_Px_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2y_S_M4_vrr = PAY*I_ESP_Py_S_M4_vrr-PRY*I_ESP_Py_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2z_S_M4_vrr = PAZ*I_ESP_Pz_S_M4_vrr-PRZ*I_ESP_Pz_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         ************************************************************/
        Double I_ESP_F3x_S_M4_vrr = PAX*I_ESP_D2x_S_M4_vrr-PRX*I_ESP_D2x_S_M5_vrr+2*oned2z*I_ESP_Px_S_M4_vrr-2*oned2z*I_ESP_Px_S_M5_vrr;
        Double I_ESP_F3y_S_M4_vrr = PAY*I_ESP_D2y_S_M4_vrr-PRY*I_ESP_D2y_S_M5_vrr+2*oned2z*I_ESP_Py_S_M4_vrr-2*oned2z*I_ESP_Py_S_M5_vrr;
        Double I_ESP_F3z_S_M4_vrr = PAZ*I_ESP_D2z_S_M4_vrr-PRZ*I_ESP_D2z_S_M5_vrr+2*oned2z*I_ESP_Pz_S_M4_vrr-2*oned2z*I_ESP_Pz_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         ************************************************************/
        Double I_ESP_G4x_S_M4_vrr = PAX*I_ESP_F3x_S_M4_vrr-PRX*I_ESP_F3x_S_M5_vrr+3*oned2z*I_ESP_D2x_S_M4_vrr-3*oned2z*I_ESP_D2x_S_M5_vrr;
        Double I_ESP_G3xy_S_M4_vrr = PAY*I_ESP_F3x_S_M4_vrr-PRY*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G3xz_S_M4_vrr = PAZ*I_ESP_F3x_S_M4_vrr-PRZ*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G4y_S_M4_vrr = PAY*I_ESP_F3y_S_M4_vrr-PRY*I_ESP_F3y_S_M5_vrr+3*oned2z*I_ESP_D2y_S_M4_vrr-3*oned2z*I_ESP_D2y_S_M5_vrr;
        Double I_ESP_G3yz_S_M4_vrr = PAZ*I_ESP_F3y_S_M4_vrr-PRZ*I_ESP_F3y_S_M5_vrr;
        Double I_ESP_G4z_S_M4_vrr = PAZ*I_ESP_F3z_S_M4_vrr-PRZ*I_ESP_F3z_S_M5_vrr+3*oned2z*I_ESP_D2z_S_M4_vrr-3*oned2z*I_ESP_D2z_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_Px_S_M3_vrr = PAX*I_ESP_S_S_M3_vrr-PRX*I_ESP_S_S_M4_vrr;
        Double I_ESP_Py_S_M3_vrr = PAY*I_ESP_S_S_M3_vrr-PRY*I_ESP_S_S_M4_vrr;
        Double I_ESP_Pz_S_M3_vrr = PAZ*I_ESP_S_S_M3_vrr-PRZ*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_D2x_S_M3_vrr = PAX*I_ESP_Px_S_M3_vrr-PRX*I_ESP_Px_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2y_S_M3_vrr = PAY*I_ESP_Py_S_M3_vrr-PRY*I_ESP_Py_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2z_S_M3_vrr = PAZ*I_ESP_Pz_S_M3_vrr-PRZ*I_ESP_Pz_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 6 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         ************************************************************/
        Double I_ESP_F3x_S_M3_vrr = PAX*I_ESP_D2x_S_M3_vrr-PRX*I_ESP_D2x_S_M4_vrr+2*oned2z*I_ESP_Px_S_M3_vrr-2*oned2z*I_ESP_Px_S_M4_vrr;
        Double I_ESP_F2xy_S_M3_vrr = PAY*I_ESP_D2x_S_M3_vrr-PRY*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_F3y_S_M3_vrr = PAY*I_ESP_D2y_S_M3_vrr-PRY*I_ESP_D2y_S_M4_vrr+2*oned2z*I_ESP_Py_S_M3_vrr-2*oned2z*I_ESP_Py_S_M4_vrr;
        Double I_ESP_F3z_S_M3_vrr = PAZ*I_ESP_D2z_S_M3_vrr-PRZ*I_ESP_D2z_S_M4_vrr+2*oned2z*I_ESP_Pz_S_M3_vrr-2*oned2z*I_ESP_Pz_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         ************************************************************/
        Double I_ESP_G4x_S_M3_vrr = PAX*I_ESP_F3x_S_M3_vrr-PRX*I_ESP_F3x_S_M4_vrr+3*oned2z*I_ESP_D2x_S_M3_vrr-3*oned2z*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_G3xy_S_M3_vrr = PAY*I_ESP_F3x_S_M3_vrr-PRY*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_G3xz_S_M3_vrr = PAZ*I_ESP_F3x_S_M3_vrr-PRZ*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Gx3y_S_M3_vrr = PAX*I_ESP_F3y_S_M3_vrr-PRX*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Gx3z_S_M3_vrr = PAX*I_ESP_F3z_S_M3_vrr-PRX*I_ESP_F3z_S_M4_vrr;
        Double I_ESP_G4y_S_M3_vrr = PAY*I_ESP_F3y_S_M3_vrr-PRY*I_ESP_F3y_S_M4_vrr+3*oned2z*I_ESP_D2y_S_M3_vrr-3*oned2z*I_ESP_D2y_S_M4_vrr;
        Double I_ESP_G3yz_S_M3_vrr = PAZ*I_ESP_F3y_S_M3_vrr-PRZ*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_G4z_S_M3_vrr = PAZ*I_ESP_F3z_S_M3_vrr-PRZ*I_ESP_F3z_S_M4_vrr+3*oned2z*I_ESP_D2z_S_M3_vrr-3*oned2z*I_ESP_D2z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         ************************************************************/
        Double I_ESP_H5x_S_M3_vrr = PAX*I_ESP_G4x_S_M3_vrr-PRX*I_ESP_G4x_S_M4_vrr+4*oned2z*I_ESP_F3x_S_M3_vrr-4*oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H4xy_S_M3_vrr = PAY*I_ESP_G4x_S_M3_vrr-PRY*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H4xz_S_M3_vrr = PAZ*I_ESP_G4x_S_M3_vrr-PRZ*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H3x2y_S_M3_vrr = PAY*I_ESP_G3xy_S_M3_vrr-PRY*I_ESP_G3xy_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H3x2z_S_M3_vrr = PAZ*I_ESP_G3xz_S_M3_vrr-PRZ*I_ESP_G3xz_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Hx4y_S_M3_vrr = PAX*I_ESP_G4y_S_M3_vrr-PRX*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_Hx4z_S_M3_vrr = PAX*I_ESP_G4z_S_M3_vrr-PRX*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5y_S_M3_vrr = PAY*I_ESP_G4y_S_M3_vrr-PRY*I_ESP_G4y_S_M4_vrr+4*oned2z*I_ESP_F3y_S_M3_vrr-4*oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_H4yz_S_M3_vrr = PAZ*I_ESP_G4y_S_M3_vrr-PRZ*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_H3y2z_S_M3_vrr = PAZ*I_ESP_G3yz_S_M3_vrr-PRZ*I_ESP_G3yz_S_M4_vrr+oned2z*I_ESP_F3y_S_M3_vrr-oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Hy4z_S_M3_vrr = PAY*I_ESP_G4z_S_M3_vrr-PRY*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5z_S_M3_vrr = PAZ*I_ESP_G4z_S_M3_vrr-PRZ*I_ESP_G4z_S_M4_vrr+4*oned2z*I_ESP_F3z_S_M3_vrr-4*oned2z*I_ESP_F3z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_Px_S_M2_vrr = PAX*I_ESP_S_S_M2_vrr-PRX*I_ESP_S_S_M3_vrr;
        Double I_ESP_Py_S_M2_vrr = PAY*I_ESP_S_S_M2_vrr-PRY*I_ESP_S_S_M3_vrr;
        Double I_ESP_Pz_S_M2_vrr = PAZ*I_ESP_S_S_M2_vrr-PRZ*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_D2x_S_M2_vrr = PAX*I_ESP_Px_S_M2_vrr-PRX*I_ESP_Px_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2y_S_M2_vrr = PAY*I_ESP_Py_S_M2_vrr-PRY*I_ESP_Py_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2z_S_M2_vrr = PAZ*I_ESP_Pz_S_M2_vrr-PRZ*I_ESP_Pz_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         ************************************************************/
        Double I_ESP_F3x_S_M2_vrr = PAX*I_ESP_D2x_S_M2_vrr-PRX*I_ESP_D2x_S_M3_vrr+2*oned2z*I_ESP_Px_S_M2_vrr-2*oned2z*I_ESP_Px_S_M3_vrr;
        Double I_ESP_F2xy_S_M2_vrr = PAY*I_ESP_D2x_S_M2_vrr-PRY*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F2xz_S_M2_vrr = PAZ*I_ESP_D2x_S_M2_vrr-PRZ*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F3y_S_M2_vrr = PAY*I_ESP_D2y_S_M2_vrr-PRY*I_ESP_D2y_S_M3_vrr+2*oned2z*I_ESP_Py_S_M2_vrr-2*oned2z*I_ESP_Py_S_M3_vrr;
        Double I_ESP_F2yz_S_M2_vrr = PAZ*I_ESP_D2y_S_M2_vrr-PRZ*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_F3z_S_M2_vrr = PAZ*I_ESP_D2z_S_M2_vrr-PRZ*I_ESP_D2z_S_M3_vrr+2*oned2z*I_ESP_Pz_S_M2_vrr-2*oned2z*I_ESP_Pz_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         ************************************************************/
        Double I_ESP_G4x_S_M2_vrr = PAX*I_ESP_F3x_S_M2_vrr-PRX*I_ESP_F3x_S_M3_vrr+3*oned2z*I_ESP_D2x_S_M2_vrr-3*oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_G3xy_S_M2_vrr = PAY*I_ESP_F3x_S_M2_vrr-PRY*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G3xz_S_M2_vrr = PAZ*I_ESP_F3x_S_M2_vrr-PRZ*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G2x2y_S_M2_vrr = PAY*I_ESP_F2xy_S_M2_vrr-PRY*I_ESP_F2xy_S_M3_vrr+oned2z*I_ESP_D2x_S_M2_vrr-oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_Gx3y_S_M2_vrr = PAX*I_ESP_F3y_S_M2_vrr-PRX*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gx3z_S_M2_vrr = PAX*I_ESP_F3z_S_M2_vrr-PRX*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4y_S_M2_vrr = PAY*I_ESP_F3y_S_M2_vrr-PRY*I_ESP_F3y_S_M3_vrr+3*oned2z*I_ESP_D2y_S_M2_vrr-3*oned2z*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_G3yz_S_M2_vrr = PAZ*I_ESP_F3y_S_M2_vrr-PRZ*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gy3z_S_M2_vrr = PAY*I_ESP_F3z_S_M2_vrr-PRY*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4z_S_M2_vrr = PAZ*I_ESP_F3z_S_M2_vrr-PRZ*I_ESP_F3z_S_M3_vrr+3*oned2z*I_ESP_D2z_S_M2_vrr-3*oned2z*I_ESP_D2z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         ************************************************************/
        Double I_ESP_H5x_S_M2_vrr = PAX*I_ESP_G4x_S_M2_vrr-PRX*I_ESP_G4x_S_M3_vrr+4*oned2z*I_ESP_F3x_S_M2_vrr-4*oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H4xy_S_M2_vrr = PAY*I_ESP_G4x_S_M2_vrr-PRY*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H4xz_S_M2_vrr = PAZ*I_ESP_G4x_S_M2_vrr-PRZ*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H3x2y_S_M2_vrr = PAY*I_ESP_G3xy_S_M2_vrr-PRY*I_ESP_G3xy_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H3x2z_S_M2_vrr = PAZ*I_ESP_G3xz_S_M2_vrr-PRZ*I_ESP_G3xz_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H2x3y_S_M2_vrr = PAX*I_ESP_Gx3y_S_M2_vrr-PRX*I_ESP_Gx3y_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H2x3z_S_M2_vrr = PAX*I_ESP_Gx3z_S_M2_vrr-PRX*I_ESP_Gx3z_S_M3_vrr+oned2z*I_ESP_F3z_S_M2_vrr-oned2z*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_Hx4y_S_M2_vrr = PAX*I_ESP_G4y_S_M2_vrr-PRX*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_Hx4z_S_M2_vrr = PAX*I_ESP_G4z_S_M2_vrr-PRX*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5y_S_M2_vrr = PAY*I_ESP_G4y_S_M2_vrr-PRY*I_ESP_G4y_S_M3_vrr+4*oned2z*I_ESP_F3y_S_M2_vrr-4*oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H4yz_S_M2_vrr = PAZ*I_ESP_G4y_S_M2_vrr-PRZ*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_H3y2z_S_M2_vrr = PAZ*I_ESP_G3yz_S_M2_vrr-PRZ*I_ESP_G3yz_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Hy4z_S_M2_vrr = PAY*I_ESP_G4z_S_M2_vrr-PRY*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5z_S_M2_vrr = PAZ*I_ESP_G4z_S_M2_vrr-PRZ*I_ESP_G4z_S_M3_vrr+4*oned2z*I_ESP_F3z_S_M2_vrr-4*oned2z*I_ESP_F3z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 10 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         ************************************************************/
        Double I_ESP_I6x_S_M2_vrr = PAX*I_ESP_H5x_S_M2_vrr-PRX*I_ESP_H5x_S_M3_vrr+5*oned2z*I_ESP_G4x_S_M2_vrr-5*oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I5xy_S_M2_vrr = PAY*I_ESP_H5x_S_M2_vrr-PRY*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I5xz_S_M2_vrr = PAZ*I_ESP_H5x_S_M2_vrr-PRZ*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I4x2y_S_M2_vrr = PAY*I_ESP_H4xy_S_M2_vrr-PRY*I_ESP_H4xy_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I4x2z_S_M2_vrr = PAZ*I_ESP_H4xz_S_M2_vrr-PRZ*I_ESP_H4xz_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I3x3y_S_M2_vrr = PAY*I_ESP_H3x2y_S_M2_vrr-PRY*I_ESP_H3x2y_S_M3_vrr+2*oned2z*I_ESP_G3xy_S_M2_vrr-2*oned2z*I_ESP_G3xy_S_M3_vrr;
        Double I_ESP_I3x3z_S_M2_vrr = PAZ*I_ESP_H3x2z_S_M2_vrr-PRZ*I_ESP_H3x2z_S_M3_vrr+2*oned2z*I_ESP_G3xz_S_M2_vrr-2*oned2z*I_ESP_G3xz_S_M3_vrr;
        Double I_ESP_I2x4y_S_M2_vrr = PAX*I_ESP_Hx4y_S_M2_vrr-PRX*I_ESP_Hx4y_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I2x4z_S_M2_vrr = PAX*I_ESP_Hx4z_S_M2_vrr-PRX*I_ESP_Hx4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Ix5y_S_M2_vrr = PAX*I_ESP_H5y_S_M2_vrr-PRX*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_Ix5z_S_M2_vrr = PAX*I_ESP_H5z_S_M2_vrr-PRX*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6y_S_M2_vrr = PAY*I_ESP_H5y_S_M2_vrr-PRY*I_ESP_H5y_S_M3_vrr+5*oned2z*I_ESP_G4y_S_M2_vrr-5*oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I5yz_S_M2_vrr = PAZ*I_ESP_H5y_S_M2_vrr-PRZ*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_I4y2z_S_M2_vrr = PAZ*I_ESP_H4yz_S_M2_vrr-PRZ*I_ESP_H4yz_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I3y3z_S_M2_vrr = PAZ*I_ESP_H3y2z_S_M2_vrr-PRZ*I_ESP_H3y2z_S_M3_vrr+2*oned2z*I_ESP_G3yz_S_M2_vrr-2*oned2z*I_ESP_G3yz_S_M3_vrr;
        Double I_ESP_I2y4z_S_M2_vrr = PAY*I_ESP_Hy4z_S_M2_vrr-PRY*I_ESP_Hy4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Iy5z_S_M2_vrr = PAY*I_ESP_H5z_S_M2_vrr-PRY*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6z_S_M2_vrr = PAZ*I_ESP_H5z_S_M2_vrr-PRZ*I_ESP_H5z_S_M3_vrr+5*oned2z*I_ESP_G4z_S_M2_vrr-5*oned2z*I_ESP_G4z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_Px_S_M1_vrr = PAX*I_ESP_S_S_M1_vrr-PRX*I_ESP_S_S_M2_vrr;
        Double I_ESP_Py_S_M1_vrr = PAY*I_ESP_S_S_M1_vrr-PRY*I_ESP_S_S_M2_vrr;
        Double I_ESP_Pz_S_M1_vrr = PAZ*I_ESP_S_S_M1_vrr-PRZ*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_Dxy_S_M1_vrr = PAY*I_ESP_Px_S_M1_vrr-PRY*I_ESP_Px_S_M2_vrr;
        Double I_ESP_D2y_S_M1_vrr = PAY*I_ESP_Py_S_M1_vrr-PRY*I_ESP_Py_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2z_S_M1_vrr = PAZ*I_ESP_Pz_S_M1_vrr-PRZ*I_ESP_Pz_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Fx2y_S_M1_vrr = PAX*I_ESP_D2y_S_M1_vrr-PRX*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Fx2z_S_M1_vrr = PAX*I_ESP_D2z_S_M1_vrr-PRX*I_ESP_D2z_S_M2_vrr;
        Double I_ESP_F3y_S_M1_vrr = PAY*I_ESP_D2y_S_M1_vrr-PRY*I_ESP_D2y_S_M2_vrr+2*oned2z*I_ESP_Py_S_M1_vrr-2*oned2z*I_ESP_Py_S_M2_vrr;
        Double I_ESP_F2yz_S_M1_vrr = PAZ*I_ESP_D2y_S_M1_vrr-PRZ*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_F3z_S_M1_vrr = PAZ*I_ESP_D2z_S_M1_vrr-PRZ*I_ESP_D2z_S_M2_vrr+2*oned2z*I_ESP_Pz_S_M1_vrr-2*oned2z*I_ESP_Pz_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         ************************************************************/
        Double I_ESP_G4x_S_M1_vrr = PAX*I_ESP_F3x_S_M1_vrr-PRX*I_ESP_F3x_S_M2_vrr+3*oned2z*I_ESP_D2x_S_M1_vrr-3*oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G3xy_S_M1_vrr = PAY*I_ESP_F3x_S_M1_vrr-PRY*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G3xz_S_M1_vrr = PAZ*I_ESP_F3x_S_M1_vrr-PRZ*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G2x2y_S_M1_vrr = PAY*I_ESP_F2xy_S_M1_vrr-PRY*I_ESP_F2xy_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G2x2z_S_M1_vrr = PAZ*I_ESP_F2xz_S_M1_vrr-PRZ*I_ESP_F2xz_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Gx3y_S_M1_vrr = PAX*I_ESP_F3y_S_M1_vrr-PRX*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_Gx3z_S_M1_vrr = PAX*I_ESP_F3z_S_M1_vrr-PRX*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4y_S_M1_vrr = PAY*I_ESP_F3y_S_M1_vrr-PRY*I_ESP_F3y_S_M2_vrr+3*oned2z*I_ESP_D2y_S_M1_vrr-3*oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_G3yz_S_M1_vrr = PAZ*I_ESP_F3y_S_M1_vrr-PRZ*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_G2y2z_S_M1_vrr = PAZ*I_ESP_F2yz_S_M1_vrr-PRZ*I_ESP_F2yz_S_M2_vrr+oned2z*I_ESP_D2y_S_M1_vrr-oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Gy3z_S_M1_vrr = PAY*I_ESP_F3z_S_M1_vrr-PRY*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4z_S_M1_vrr = PAZ*I_ESP_F3z_S_M1_vrr-PRZ*I_ESP_F3z_S_M2_vrr+3*oned2z*I_ESP_D2z_S_M1_vrr-3*oned2z*I_ESP_D2z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         ************************************************************/
        Double I_ESP_H5x_S_M1_vrr = PAX*I_ESP_G4x_S_M1_vrr-PRX*I_ESP_G4x_S_M2_vrr+4*oned2z*I_ESP_F3x_S_M1_vrr-4*oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H4xy_S_M1_vrr = PAY*I_ESP_G4x_S_M1_vrr-PRY*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H4xz_S_M1_vrr = PAZ*I_ESP_G4x_S_M1_vrr-PRZ*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H3x2y_S_M1_vrr = PAY*I_ESP_G3xy_S_M1_vrr-PRY*I_ESP_G3xy_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H3x2z_S_M1_vrr = PAZ*I_ESP_G3xz_S_M1_vrr-PRZ*I_ESP_G3xz_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H2x3y_S_M1_vrr = PAX*I_ESP_Gx3y_S_M1_vrr-PRX*I_ESP_Gx3y_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2x2yz_S_M1_vrr = PAZ*I_ESP_G2x2y_S_M1_vrr-PRZ*I_ESP_G2x2y_S_M2_vrr;
        Double I_ESP_H2x3z_S_M1_vrr = PAX*I_ESP_Gx3z_S_M1_vrr-PRX*I_ESP_Gx3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hx4y_S_M1_vrr = PAX*I_ESP_G4y_S_M1_vrr-PRX*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_Hx4z_S_M1_vrr = PAX*I_ESP_G4z_S_M1_vrr-PRX*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5y_S_M1_vrr = PAY*I_ESP_G4y_S_M1_vrr-PRY*I_ESP_G4y_S_M2_vrr+4*oned2z*I_ESP_F3y_S_M1_vrr-4*oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H4yz_S_M1_vrr = PAZ*I_ESP_G4y_S_M1_vrr-PRZ*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_H3y2z_S_M1_vrr = PAZ*I_ESP_G3yz_S_M1_vrr-PRZ*I_ESP_G3yz_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2y3z_S_M1_vrr = PAY*I_ESP_Gy3z_S_M1_vrr-PRY*I_ESP_Gy3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hy4z_S_M1_vrr = PAY*I_ESP_G4z_S_M1_vrr-PRY*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5z_S_M1_vrr = PAZ*I_ESP_G4z_S_M1_vrr-PRZ*I_ESP_G4z_S_M2_vrr+4*oned2z*I_ESP_F3z_S_M1_vrr-4*oned2z*I_ESP_F3z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         ************************************************************/
        Double I_ESP_I6x_S_M1_vrr = PAX*I_ESP_H5x_S_M1_vrr-PRX*I_ESP_H5x_S_M2_vrr+5*oned2z*I_ESP_G4x_S_M1_vrr-5*oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I5xy_S_M1_vrr = PAY*I_ESP_H5x_S_M1_vrr-PRY*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I5xz_S_M1_vrr = PAZ*I_ESP_H5x_S_M1_vrr-PRZ*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I4x2y_S_M1_vrr = PAY*I_ESP_H4xy_S_M1_vrr-PRY*I_ESP_H4xy_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I4x2z_S_M1_vrr = PAZ*I_ESP_H4xz_S_M1_vrr-PRZ*I_ESP_H4xz_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I3x3y_S_M1_vrr = PAY*I_ESP_H3x2y_S_M1_vrr-PRY*I_ESP_H3x2y_S_M2_vrr+2*oned2z*I_ESP_G3xy_S_M1_vrr-2*oned2z*I_ESP_G3xy_S_M2_vrr;
        Double I_ESP_I3x2yz_S_M1_vrr = PAZ*I_ESP_H3x2y_S_M1_vrr-PRZ*I_ESP_H3x2y_S_M2_vrr;
        Double I_ESP_I3x3z_S_M1_vrr = PAZ*I_ESP_H3x2z_S_M1_vrr-PRZ*I_ESP_H3x2z_S_M2_vrr+2*oned2z*I_ESP_G3xz_S_M1_vrr-2*oned2z*I_ESP_G3xz_S_M2_vrr;
        Double I_ESP_I2x4y_S_M1_vrr = PAX*I_ESP_Hx4y_S_M1_vrr-PRX*I_ESP_Hx4y_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I2x3yz_S_M1_vrr = PAZ*I_ESP_H2x3y_S_M1_vrr-PRZ*I_ESP_H2x3y_S_M2_vrr;
        Double I_ESP_I2xy3z_S_M1_vrr = PAY*I_ESP_H2x3z_S_M1_vrr-PRY*I_ESP_H2x3z_S_M2_vrr;
        Double I_ESP_I2x4z_S_M1_vrr = PAX*I_ESP_Hx4z_S_M1_vrr-PRX*I_ESP_Hx4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Ix5y_S_M1_vrr = PAX*I_ESP_H5y_S_M1_vrr-PRX*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_Ix5z_S_M1_vrr = PAX*I_ESP_H5z_S_M1_vrr-PRX*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6y_S_M1_vrr = PAY*I_ESP_H5y_S_M1_vrr-PRY*I_ESP_H5y_S_M2_vrr+5*oned2z*I_ESP_G4y_S_M1_vrr-5*oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I5yz_S_M1_vrr = PAZ*I_ESP_H5y_S_M1_vrr-PRZ*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_I4y2z_S_M1_vrr = PAZ*I_ESP_H4yz_S_M1_vrr-PRZ*I_ESP_H4yz_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I3y3z_S_M1_vrr = PAZ*I_ESP_H3y2z_S_M1_vrr-PRZ*I_ESP_H3y2z_S_M2_vrr+2*oned2z*I_ESP_G3yz_S_M1_vrr-2*oned2z*I_ESP_G3yz_S_M2_vrr;
        Double I_ESP_I2y4z_S_M1_vrr = PAY*I_ESP_Hy4z_S_M1_vrr-PRY*I_ESP_Hy4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Iy5z_S_M1_vrr = PAY*I_ESP_H5z_S_M1_vrr-PRY*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6z_S_M1_vrr = PAZ*I_ESP_H5z_S_M1_vrr-PRZ*I_ESP_H5z_S_M2_vrr+5*oned2z*I_ESP_G4z_S_M1_vrr-5*oned2z*I_ESP_G4z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         ************************************************************/
        Double I_ESP_K7x_S_M1_vrr = PAX*I_ESP_I6x_S_M1_vrr-PRX*I_ESP_I6x_S_M2_vrr+6*oned2z*I_ESP_H5x_S_M1_vrr-6*oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K6xy_S_M1_vrr = PAY*I_ESP_I6x_S_M1_vrr-PRY*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K6xz_S_M1_vrr = PAZ*I_ESP_I6x_S_M1_vrr-PRZ*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K5x2y_S_M1_vrr = PAY*I_ESP_I5xy_S_M1_vrr-PRY*I_ESP_I5xy_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K5x2z_S_M1_vrr = PAZ*I_ESP_I5xz_S_M1_vrr-PRZ*I_ESP_I5xz_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K4x3y_S_M1_vrr = PAY*I_ESP_I4x2y_S_M1_vrr-PRY*I_ESP_I4x2y_S_M2_vrr+2*oned2z*I_ESP_H4xy_S_M1_vrr-2*oned2z*I_ESP_H4xy_S_M2_vrr;
        Double I_ESP_K4x2yz_S_M1_vrr = PAZ*I_ESP_I4x2y_S_M1_vrr-PRZ*I_ESP_I4x2y_S_M2_vrr;
        Double I_ESP_K4x3z_S_M1_vrr = PAZ*I_ESP_I4x2z_S_M1_vrr-PRZ*I_ESP_I4x2z_S_M2_vrr+2*oned2z*I_ESP_H4xz_S_M1_vrr-2*oned2z*I_ESP_H4xz_S_M2_vrr;
        Double I_ESP_K3x4y_S_M1_vrr = PAX*I_ESP_I2x4y_S_M1_vrr-PRX*I_ESP_I2x4y_S_M2_vrr+2*oned2z*I_ESP_Hx4y_S_M1_vrr-2*oned2z*I_ESP_Hx4y_S_M2_vrr;
        Double I_ESP_K3x3yz_S_M1_vrr = PAZ*I_ESP_I3x3y_S_M1_vrr-PRZ*I_ESP_I3x3y_S_M2_vrr;
        Double I_ESP_K3xy3z_S_M1_vrr = PAY*I_ESP_I3x3z_S_M1_vrr-PRY*I_ESP_I3x3z_S_M2_vrr;
        Double I_ESP_K3x4z_S_M1_vrr = PAX*I_ESP_I2x4z_S_M1_vrr-PRX*I_ESP_I2x4z_S_M2_vrr+2*oned2z*I_ESP_Hx4z_S_M1_vrr-2*oned2z*I_ESP_Hx4z_S_M2_vrr;
        Double I_ESP_K2x5y_S_M1_vrr = PAX*I_ESP_Ix5y_S_M1_vrr-PRX*I_ESP_Ix5y_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K2x4yz_S_M1_vrr = PAZ*I_ESP_I2x4y_S_M1_vrr-PRZ*I_ESP_I2x4y_S_M2_vrr;
        Double I_ESP_K2xy4z_S_M1_vrr = PAY*I_ESP_I2x4z_S_M1_vrr-PRY*I_ESP_I2x4z_S_M2_vrr;
        Double I_ESP_K2x5z_S_M1_vrr = PAX*I_ESP_Ix5z_S_M1_vrr-PRX*I_ESP_Ix5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Kx6y_S_M1_vrr = PAX*I_ESP_I6y_S_M1_vrr-PRX*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_Kx3y3z_S_M1_vrr = PAX*I_ESP_I3y3z_S_M1_vrr-PRX*I_ESP_I3y3z_S_M2_vrr;
        Double I_ESP_Kx6z_S_M1_vrr = PAX*I_ESP_I6z_S_M1_vrr-PRX*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7y_S_M1_vrr = PAY*I_ESP_I6y_S_M1_vrr-PRY*I_ESP_I6y_S_M2_vrr+6*oned2z*I_ESP_H5y_S_M1_vrr-6*oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K6yz_S_M1_vrr = PAZ*I_ESP_I6y_S_M1_vrr-PRZ*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_K5y2z_S_M1_vrr = PAZ*I_ESP_I5yz_S_M1_vrr-PRZ*I_ESP_I5yz_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K4y3z_S_M1_vrr = PAZ*I_ESP_I4y2z_S_M1_vrr-PRZ*I_ESP_I4y2z_S_M2_vrr+2*oned2z*I_ESP_H4yz_S_M1_vrr-2*oned2z*I_ESP_H4yz_S_M2_vrr;
        Double I_ESP_K3y4z_S_M1_vrr = PAY*I_ESP_I2y4z_S_M1_vrr-PRY*I_ESP_I2y4z_S_M2_vrr+2*oned2z*I_ESP_Hy4z_S_M1_vrr-2*oned2z*I_ESP_Hy4z_S_M2_vrr;
        Double I_ESP_K2y5z_S_M1_vrr = PAY*I_ESP_Iy5z_S_M1_vrr-PRY*I_ESP_Iy5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Ky6z_S_M1_vrr = PAY*I_ESP_I6z_S_M1_vrr-PRY*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7z_S_M1_vrr = PAZ*I_ESP_I6z_S_M1_vrr-PRZ*I_ESP_I6z_S_M2_vrr+6*oned2z*I_ESP_H5z_S_M1_vrr-6*oned2z*I_ESP_H5z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_Px_S_vrr = PAX*I_ESP_S_S_vrr-PRX*I_ESP_S_S_M1_vrr;
        Double I_ESP_Py_S_vrr = PAY*I_ESP_S_S_vrr-PRY*I_ESP_S_S_M1_vrr;
        Double I_ESP_Pz_S_vrr = PAZ*I_ESP_S_S_vrr-PRZ*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fxyz_S_vrr = PAZ*I_ESP_Dxy_S_vrr-PRZ*I_ESP_Dxy_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fy2z_S_vrr = PAY*I_ESP_D2z_S_vrr-PRY*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2xyz_S_vrr = PAZ*I_ESP_F2xy_S_vrr-PRZ*I_ESP_F2xy_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_Gx2yz_S_vrr = PAZ*I_ESP_Fx2y_S_vrr-PRZ*I_ESP_Fx2y_S_M1_vrr;
        Double I_ESP_Gxy2z_S_vrr = PAY*I_ESP_Fx2z_S_vrr-PRY*I_ESP_Fx2z_S_M1_vrr;
        Double I_ESP_Gx3z_S_vrr = PAX*I_ESP_F3z_S_vrr-PRX*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4y_S_vrr = PAY*I_ESP_F3y_S_vrr-PRY*I_ESP_F3y_S_M1_vrr+3*oned2z*I_ESP_D2y_S_vrr-3*oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_G3yz_S_vrr = PAZ*I_ESP_F3y_S_vrr-PRZ*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_G2y2z_S_vrr = PAZ*I_ESP_F2yz_S_vrr-PRZ*I_ESP_F2yz_S_M1_vrr+oned2z*I_ESP_D2y_S_vrr-oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Gy3z_S_vrr = PAY*I_ESP_F3z_S_vrr-PRY*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4z_S_vrr = PAZ*I_ESP_F3z_S_vrr-PRZ*I_ESP_F3z_S_M1_vrr+3*oned2z*I_ESP_D2z_S_vrr-3*oned2z*I_ESP_D2z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         ************************************************************/
        Double I_ESP_H5x_S_vrr = PAX*I_ESP_G4x_S_vrr-PRX*I_ESP_G4x_S_M1_vrr+4*oned2z*I_ESP_F3x_S_vrr-4*oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H4xy_S_vrr = PAY*I_ESP_G4x_S_vrr-PRY*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H4xz_S_vrr = PAZ*I_ESP_G4x_S_vrr-PRZ*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H3x2y_S_vrr = PAY*I_ESP_G3xy_S_vrr-PRY*I_ESP_G3xy_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H3xyz_S_vrr = PAZ*I_ESP_G3xy_S_vrr-PRZ*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_H3x2z_S_vrr = PAZ*I_ESP_G3xz_S_vrr-PRZ*I_ESP_G3xz_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H2x3y_S_vrr = PAX*I_ESP_Gx3y_S_vrr-PRX*I_ESP_Gx3y_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2x2yz_S_vrr = PAZ*I_ESP_G2x2y_S_vrr-PRZ*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_H2xy2z_S_vrr = PAY*I_ESP_G2x2z_S_vrr-PRY*I_ESP_G2x2z_S_M1_vrr;
        Double I_ESP_H2x3z_S_vrr = PAX*I_ESP_Gx3z_S_vrr-PRX*I_ESP_Gx3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hx4y_S_vrr = PAX*I_ESP_G4y_S_vrr-PRX*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_Hx3yz_S_vrr = PAZ*I_ESP_Gx3y_S_vrr-PRZ*I_ESP_Gx3y_S_M1_vrr;
        Double I_ESP_Hx2y2z_S_vrr = PAX*I_ESP_G2y2z_S_vrr-PRX*I_ESP_G2y2z_S_M1_vrr;
        Double I_ESP_Hxy3z_S_vrr = PAY*I_ESP_Gx3z_S_vrr-PRY*I_ESP_Gx3z_S_M1_vrr;
        Double I_ESP_Hx4z_S_vrr = PAX*I_ESP_G4z_S_vrr-PRX*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5y_S_vrr = PAY*I_ESP_G4y_S_vrr-PRY*I_ESP_G4y_S_M1_vrr+4*oned2z*I_ESP_F3y_S_vrr-4*oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H4yz_S_vrr = PAZ*I_ESP_G4y_S_vrr-PRZ*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_H3y2z_S_vrr = PAZ*I_ESP_G3yz_S_vrr-PRZ*I_ESP_G3yz_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2y3z_S_vrr = PAY*I_ESP_Gy3z_S_vrr-PRY*I_ESP_Gy3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hy4z_S_vrr = PAY*I_ESP_G4z_S_vrr-PRY*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5z_S_vrr = PAZ*I_ESP_G4z_S_vrr-PRZ*I_ESP_G4z_S_M1_vrr+4*oned2z*I_ESP_F3z_S_vrr-4*oned2z*I_ESP_F3z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         ************************************************************/
        Double I_ESP_I6x_S_vrr = PAX*I_ESP_H5x_S_vrr-PRX*I_ESP_H5x_S_M1_vrr+5*oned2z*I_ESP_G4x_S_vrr-5*oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I5xy_S_vrr = PAY*I_ESP_H5x_S_vrr-PRY*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I5xz_S_vrr = PAZ*I_ESP_H5x_S_vrr-PRZ*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I4x2y_S_vrr = PAY*I_ESP_H4xy_S_vrr-PRY*I_ESP_H4xy_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I4xyz_S_vrr = PAZ*I_ESP_H4xy_S_vrr-PRZ*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_I4x2z_S_vrr = PAZ*I_ESP_H4xz_S_vrr-PRZ*I_ESP_H4xz_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I3x3y_S_vrr = PAY*I_ESP_H3x2y_S_vrr-PRY*I_ESP_H3x2y_S_M1_vrr+2*oned2z*I_ESP_G3xy_S_vrr-2*oned2z*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_I3x2yz_S_vrr = PAZ*I_ESP_H3x2y_S_vrr-PRZ*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_I3xy2z_S_vrr = PAY*I_ESP_H3x2z_S_vrr-PRY*I_ESP_H3x2z_S_M1_vrr;
        Double I_ESP_I3x3z_S_vrr = PAZ*I_ESP_H3x2z_S_vrr-PRZ*I_ESP_H3x2z_S_M1_vrr+2*oned2z*I_ESP_G3xz_S_vrr-2*oned2z*I_ESP_G3xz_S_M1_vrr;
        Double I_ESP_I2x4y_S_vrr = PAX*I_ESP_Hx4y_S_vrr-PRX*I_ESP_Hx4y_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I2x3yz_S_vrr = PAZ*I_ESP_H2x3y_S_vrr-PRZ*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_I2x2y2z_S_vrr = PAZ*I_ESP_H2x2yz_S_vrr-PRZ*I_ESP_H2x2yz_S_M1_vrr+oned2z*I_ESP_G2x2y_S_vrr-oned2z*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_I2xy3z_S_vrr = PAY*I_ESP_H2x3z_S_vrr-PRY*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_I2x4z_S_vrr = PAX*I_ESP_Hx4z_S_vrr-PRX*I_ESP_Hx4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Ix5y_S_vrr = PAX*I_ESP_H5y_S_vrr-PRX*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_Ix4yz_S_vrr = PAZ*I_ESP_Hx4y_S_vrr-PRZ*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_Ix3y2z_S_vrr = PAX*I_ESP_H3y2z_S_vrr-PRX*I_ESP_H3y2z_S_M1_vrr;
        Double I_ESP_Ix2y3z_S_vrr = PAX*I_ESP_H2y3z_S_vrr-PRX*I_ESP_H2y3z_S_M1_vrr;
        Double I_ESP_Ixy4z_S_vrr = PAY*I_ESP_Hx4z_S_vrr-PRY*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_Ix5z_S_vrr = PAX*I_ESP_H5z_S_vrr-PRX*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6y_S_vrr = PAY*I_ESP_H5y_S_vrr-PRY*I_ESP_H5y_S_M1_vrr+5*oned2z*I_ESP_G4y_S_vrr-5*oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I5yz_S_vrr = PAZ*I_ESP_H5y_S_vrr-PRZ*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_I4y2z_S_vrr = PAZ*I_ESP_H4yz_S_vrr-PRZ*I_ESP_H4yz_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I3y3z_S_vrr = PAZ*I_ESP_H3y2z_S_vrr-PRZ*I_ESP_H3y2z_S_M1_vrr+2*oned2z*I_ESP_G3yz_S_vrr-2*oned2z*I_ESP_G3yz_S_M1_vrr;
        Double I_ESP_I2y4z_S_vrr = PAY*I_ESP_Hy4z_S_vrr-PRY*I_ESP_Hy4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Iy5z_S_vrr = PAY*I_ESP_H5z_S_vrr-PRY*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6z_S_vrr = PAZ*I_ESP_H5z_S_vrr-PRZ*I_ESP_H5z_S_M1_vrr+5*oned2z*I_ESP_G4z_S_vrr-5*oned2z*I_ESP_G4z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         ************************************************************/
        Double I_ESP_K7x_S_vrr = PAX*I_ESP_I6x_S_vrr-PRX*I_ESP_I6x_S_M1_vrr+6*oned2z*I_ESP_H5x_S_vrr-6*oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K6xy_S_vrr = PAY*I_ESP_I6x_S_vrr-PRY*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K6xz_S_vrr = PAZ*I_ESP_I6x_S_vrr-PRZ*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K5x2y_S_vrr = PAY*I_ESP_I5xy_S_vrr-PRY*I_ESP_I5xy_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K5xyz_S_vrr = PAZ*I_ESP_I5xy_S_vrr-PRZ*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_K5x2z_S_vrr = PAZ*I_ESP_I5xz_S_vrr-PRZ*I_ESP_I5xz_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K4x3y_S_vrr = PAY*I_ESP_I4x2y_S_vrr-PRY*I_ESP_I4x2y_S_M1_vrr+2*oned2z*I_ESP_H4xy_S_vrr-2*oned2z*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_K4x2yz_S_vrr = PAZ*I_ESP_I4x2y_S_vrr-PRZ*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_K4xy2z_S_vrr = PAY*I_ESP_I4x2z_S_vrr-PRY*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_K4x3z_S_vrr = PAZ*I_ESP_I4x2z_S_vrr-PRZ*I_ESP_I4x2z_S_M1_vrr+2*oned2z*I_ESP_H4xz_S_vrr-2*oned2z*I_ESP_H4xz_S_M1_vrr;
        Double I_ESP_K3x4y_S_vrr = PAX*I_ESP_I2x4y_S_vrr-PRX*I_ESP_I2x4y_S_M1_vrr+2*oned2z*I_ESP_Hx4y_S_vrr-2*oned2z*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_K3x3yz_S_vrr = PAZ*I_ESP_I3x3y_S_vrr-PRZ*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_K3x2y2z_S_vrr = PAZ*I_ESP_I3x2yz_S_vrr-PRZ*I_ESP_I3x2yz_S_M1_vrr+oned2z*I_ESP_H3x2y_S_vrr-oned2z*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_K3xy3z_S_vrr = PAY*I_ESP_I3x3z_S_vrr-PRY*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_K3x4z_S_vrr = PAX*I_ESP_I2x4z_S_vrr-PRX*I_ESP_I2x4z_S_M1_vrr+2*oned2z*I_ESP_Hx4z_S_vrr-2*oned2z*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_K2x5y_S_vrr = PAX*I_ESP_Ix5y_S_vrr-PRX*I_ESP_Ix5y_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K2x4yz_S_vrr = PAZ*I_ESP_I2x4y_S_vrr-PRZ*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_K2x3y2z_S_vrr = PAZ*I_ESP_I2x3yz_S_vrr-PRZ*I_ESP_I2x3yz_S_M1_vrr+oned2z*I_ESP_H2x3y_S_vrr-oned2z*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_K2x2y3z_S_vrr = PAY*I_ESP_I2xy3z_S_vrr-PRY*I_ESP_I2xy3z_S_M1_vrr+oned2z*I_ESP_H2x3z_S_vrr-oned2z*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_K2xy4z_S_vrr = PAY*I_ESP_I2x4z_S_vrr-PRY*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_K2x5z_S_vrr = PAX*I_ESP_Ix5z_S_vrr-PRX*I_ESP_Ix5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Kx6y_S_vrr = PAX*I_ESP_I6y_S_vrr-PRX*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_Kx5yz_S_vrr = PAZ*I_ESP_Ix5y_S_vrr-PRZ*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_Kx4y2z_S_vrr = PAX*I_ESP_I4y2z_S_vrr-PRX*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_Kx3y3z_S_vrr = PAX*I_ESP_I3y3z_S_vrr-PRX*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_Kx2y4z_S_vrr = PAX*I_ESP_I2y4z_S_vrr-PRX*I_ESP_I2y4z_S_M1_vrr;
        Double I_ESP_Kxy5z_S_vrr = PAY*I_ESP_Ix5z_S_vrr-PRY*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_Kx6z_S_vrr = PAX*I_ESP_I6z_S_vrr-PRX*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7y_S_vrr = PAY*I_ESP_I6y_S_vrr-PRY*I_ESP_I6y_S_M1_vrr+6*oned2z*I_ESP_H5y_S_vrr-6*oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K6yz_S_vrr = PAZ*I_ESP_I6y_S_vrr-PRZ*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_K5y2z_S_vrr = PAZ*I_ESP_I5yz_S_vrr-PRZ*I_ESP_I5yz_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K4y3z_S_vrr = PAZ*I_ESP_I4y2z_S_vrr-PRZ*I_ESP_I4y2z_S_M1_vrr+2*oned2z*I_ESP_H4yz_S_vrr-2*oned2z*I_ESP_H4yz_S_M1_vrr;
        Double I_ESP_K3y4z_S_vrr = PAY*I_ESP_I2y4z_S_vrr-PRY*I_ESP_I2y4z_S_M1_vrr+2*oned2z*I_ESP_Hy4z_S_vrr-2*oned2z*I_ESP_Hy4z_S_M1_vrr;
        Double I_ESP_K2y5z_S_vrr = PAY*I_ESP_Iy5z_S_vrr-PRY*I_ESP_Iy5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Ky6z_S_vrr = PAY*I_ESP_I6z_S_vrr-PRY*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7z_S_vrr = PAZ*I_ESP_I6z_S_vrr-PRZ*I_ESP_I6z_S_M1_vrr+6*oned2z*I_ESP_H5z_S_vrr-6*oned2z*I_ESP_H5z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         ************************************************************/
        Double I_ESP_L8x_S_vrr = PAX*I_ESP_K7x_S_vrr-PRX*I_ESP_K7x_S_M1_vrr+7*oned2z*I_ESP_I6x_S_vrr-7*oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L7xy_S_vrr = PAY*I_ESP_K7x_S_vrr-PRY*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L7xz_S_vrr = PAZ*I_ESP_K7x_S_vrr-PRZ*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L6x2y_S_vrr = PAY*I_ESP_K6xy_S_vrr-PRY*I_ESP_K6xy_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L6xyz_S_vrr = PAZ*I_ESP_K6xy_S_vrr-PRZ*I_ESP_K6xy_S_M1_vrr;
        Double I_ESP_L6x2z_S_vrr = PAZ*I_ESP_K6xz_S_vrr-PRZ*I_ESP_K6xz_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L5x3y_S_vrr = PAY*I_ESP_K5x2y_S_vrr-PRY*I_ESP_K5x2y_S_M1_vrr+2*oned2z*I_ESP_I5xy_S_vrr-2*oned2z*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_L5x2yz_S_vrr = PAZ*I_ESP_K5x2y_S_vrr-PRZ*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_L5xy2z_S_vrr = PAY*I_ESP_K5x2z_S_vrr-PRY*I_ESP_K5x2z_S_M1_vrr;
        Double I_ESP_L5x3z_S_vrr = PAZ*I_ESP_K5x2z_S_vrr-PRZ*I_ESP_K5x2z_S_M1_vrr+2*oned2z*I_ESP_I5xz_S_vrr-2*oned2z*I_ESP_I5xz_S_M1_vrr;
        Double I_ESP_L4x4y_S_vrr = PAY*I_ESP_K4x3y_S_vrr-PRY*I_ESP_K4x3y_S_M1_vrr+3*oned2z*I_ESP_I4x2y_S_vrr-3*oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4x3yz_S_vrr = PAZ*I_ESP_K4x3y_S_vrr-PRZ*I_ESP_K4x3y_S_M1_vrr;
        Double I_ESP_L4x2y2z_S_vrr = PAZ*I_ESP_K4x2yz_S_vrr-PRZ*I_ESP_K4x2yz_S_M1_vrr+oned2z*I_ESP_I4x2y_S_vrr-oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4xy3z_S_vrr = PAY*I_ESP_K4x3z_S_vrr-PRY*I_ESP_K4x3z_S_M1_vrr;
        Double I_ESP_L4x4z_S_vrr = PAZ*I_ESP_K4x3z_S_vrr-PRZ*I_ESP_K4x3z_S_M1_vrr+3*oned2z*I_ESP_I4x2z_S_vrr-3*oned2z*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_L3x5y_S_vrr = PAX*I_ESP_K2x5y_S_vrr-PRX*I_ESP_K2x5y_S_M1_vrr+2*oned2z*I_ESP_Ix5y_S_vrr-2*oned2z*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_L3x4yz_S_vrr = PAZ*I_ESP_K3x4y_S_vrr-PRZ*I_ESP_K3x4y_S_M1_vrr;
        Double I_ESP_L3x3y2z_S_vrr = PAZ*I_ESP_K3x3yz_S_vrr-PRZ*I_ESP_K3x3yz_S_M1_vrr+oned2z*I_ESP_I3x3y_S_vrr-oned2z*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_L3x2y3z_S_vrr = PAY*I_ESP_K3xy3z_S_vrr-PRY*I_ESP_K3xy3z_S_M1_vrr+oned2z*I_ESP_I3x3z_S_vrr-oned2z*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_L3xy4z_S_vrr = PAY*I_ESP_K3x4z_S_vrr-PRY*I_ESP_K3x4z_S_M1_vrr;
        Double I_ESP_L3x5z_S_vrr = PAX*I_ESP_K2x5z_S_vrr-PRX*I_ESP_K2x5z_S_M1_vrr+2*oned2z*I_ESP_Ix5z_S_vrr-2*oned2z*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_L2x6y_S_vrr = PAX*I_ESP_Kx6y_S_vrr-PRX*I_ESP_Kx6y_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L2x5yz_S_vrr = PAZ*I_ESP_K2x5y_S_vrr-PRZ*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_L2x4y2z_S_vrr = PAZ*I_ESP_K2x4yz_S_vrr-PRZ*I_ESP_K2x4yz_S_M1_vrr+oned2z*I_ESP_I2x4y_S_vrr-oned2z*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_L2x3y3z_S_vrr = PAX*I_ESP_Kx3y3z_S_vrr-PRX*I_ESP_Kx3y3z_S_M1_vrr+oned2z*I_ESP_I3y3z_S_vrr-oned2z*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_L2x2y4z_S_vrr = PAY*I_ESP_K2xy4z_S_vrr-PRY*I_ESP_K2xy4z_S_M1_vrr+oned2z*I_ESP_I2x4z_S_vrr-oned2z*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_L2xy5z_S_vrr = PAY*I_ESP_K2x5z_S_vrr-PRY*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_L2x6z_S_vrr = PAX*I_ESP_Kx6z_S_vrr-PRX*I_ESP_Kx6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Lx7y_S_vrr = PAX*I_ESP_K7y_S_vrr-PRX*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_Lx6yz_S_vrr = PAZ*I_ESP_Kx6y_S_vrr-PRZ*I_ESP_Kx6y_S_M1_vrr;
        Double I_ESP_Lx5y2z_S_vrr = PAX*I_ESP_K5y2z_S_vrr-PRX*I_ESP_K5y2z_S_M1_vrr;
        Double I_ESP_Lx4y3z_S_vrr = PAX*I_ESP_K4y3z_S_vrr-PRX*I_ESP_K4y3z_S_M1_vrr;
        Double I_ESP_Lx3y4z_S_vrr = PAX*I_ESP_K3y4z_S_vrr-PRX*I_ESP_K3y4z_S_M1_vrr;
        Double I_ESP_Lx2y5z_S_vrr = PAX*I_ESP_K2y5z_S_vrr-PRX*I_ESP_K2y5z_S_M1_vrr;
        Double I_ESP_Lxy6z_S_vrr = PAY*I_ESP_Kx6z_S_vrr-PRY*I_ESP_Kx6z_S_M1_vrr;
        Double I_ESP_Lx7z_S_vrr = PAX*I_ESP_K7z_S_vrr-PRX*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8y_S_vrr = PAY*I_ESP_K7y_S_vrr-PRY*I_ESP_K7y_S_M1_vrr+7*oned2z*I_ESP_I6y_S_vrr-7*oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L7yz_S_vrr = PAZ*I_ESP_K7y_S_vrr-PRZ*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_L6y2z_S_vrr = PAZ*I_ESP_K6yz_S_vrr-PRZ*I_ESP_K6yz_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L5y3z_S_vrr = PAZ*I_ESP_K5y2z_S_vrr-PRZ*I_ESP_K5y2z_S_M1_vrr+2*oned2z*I_ESP_I5yz_S_vrr-2*oned2z*I_ESP_I5yz_S_M1_vrr;
        Double I_ESP_L4y4z_S_vrr = PAZ*I_ESP_K4y3z_S_vrr-PRZ*I_ESP_K4y3z_S_M1_vrr+3*oned2z*I_ESP_I4y2z_S_vrr-3*oned2z*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_L3y5z_S_vrr = PAY*I_ESP_K2y5z_S_vrr-PRY*I_ESP_K2y5z_S_M1_vrr+2*oned2z*I_ESP_Iy5z_S_vrr-2*oned2z*I_ESP_Iy5z_S_M1_vrr;
        Double I_ESP_L2y6z_S_vrr = PAY*I_ESP_Ky6z_S_vrr-PRY*I_ESP_Ky6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Ly7z_S_vrr = PAY*I_ESP_K7z_S_vrr-PRY*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8z_S_vrr = PAZ*I_ESP_K7z_S_vrr-PRZ*I_ESP_K7z_S_M1_vrr+7*oned2z*I_ESP_I6z_S_vrr-7*oned2z*I_ESP_I6z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_L_S_aa_coefs = alpha*alpha;
        I_ESP_L8x_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_K_S_aa_coefs = alpha*alpha;
        I_ESP_K7x_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_a_coefs = alpha;
        I_ESP_I6x_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_a_coefs = alpha;
        I_ESP_H5x_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_G4x_S += I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S += I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S += I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S += I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S += I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S += I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S += I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S += I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S += I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S += I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S += I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S += I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S += I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S += I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S += I_ESP_G4z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_F3x_S += I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S += I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S += I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S += I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S += I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S += I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S += I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S += I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S += I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S += I_ESP_F3z_S_vrr;
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
     * shell quartet name: SQ_ESP_F_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_S
     * RHS shell quartet name: SQ_ESP_F_S
     ************************************************************/
    Double I_ESP_F3x_Px = I_ESP_G4x_S+ABX*I_ESP_F3x_S;
    Double I_ESP_F2xy_Px = I_ESP_G3xy_S+ABX*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Px = I_ESP_G3xz_S+ABX*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Px = I_ESP_G2x2y_S+ABX*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Px = I_ESP_G2xyz_S+ABX*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Px = I_ESP_G2x2z_S+ABX*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Px = I_ESP_Gx3y_S+ABX*I_ESP_F3y_S;
    Double I_ESP_F2yz_Px = I_ESP_Gx2yz_S+ABX*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Px = I_ESP_Gxy2z_S+ABX*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Px = I_ESP_Gx3z_S+ABX*I_ESP_F3z_S;
    Double I_ESP_F3x_Py = I_ESP_G3xy_S+ABY*I_ESP_F3x_S;
    Double I_ESP_F2xy_Py = I_ESP_G2x2y_S+ABY*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Py = I_ESP_G2xyz_S+ABY*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Py = I_ESP_Gx3y_S+ABY*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Py = I_ESP_Gx2yz_S+ABY*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Py = I_ESP_Gxy2z_S+ABY*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Py = I_ESP_G4y_S+ABY*I_ESP_F3y_S;
    Double I_ESP_F2yz_Py = I_ESP_G3yz_S+ABY*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Py = I_ESP_G2y2z_S+ABY*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Py = I_ESP_Gy3z_S+ABY*I_ESP_F3z_S;
    Double I_ESP_F3x_Pz = I_ESP_G3xz_S+ABZ*I_ESP_F3x_S;
    Double I_ESP_F2xy_Pz = I_ESP_G2xyz_S+ABZ*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Pz = I_ESP_G2x2z_S+ABZ*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Pz = I_ESP_Gx2yz_S+ABZ*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Pz = I_ESP_Gxy2z_S+ABZ*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Pz = I_ESP_Gx3z_S+ABZ*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Pz = I_ESP_G3yz_S+ABZ*I_ESP_F3y_S;
    Double I_ESP_F2yz_Pz = I_ESP_G2y2z_S+ABZ*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Pz = I_ESP_Gy3z_S+ABZ*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Pz = I_ESP_G4z_S+ABZ*I_ESP_F3z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_a
     * RHS shell quartet name: SQ_ESP_H_S_a
     ************************************************************/
    Double I_ESP_H5x_Px_a = I_ESP_I6x_S_a+ABX*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Px_a = I_ESP_I5xy_S_a+ABX*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Px_a = I_ESP_I5xz_S_a+ABX*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Px_a = I_ESP_I4x2y_S_a+ABX*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Px_a = I_ESP_I4xyz_S_a+ABX*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Px_a = I_ESP_I4x2z_S_a+ABX*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Px_a = I_ESP_I3x3y_S_a+ABX*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Px_a = I_ESP_I3x2yz_S_a+ABX*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Px_a = I_ESP_I3xy2z_S_a+ABX*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Px_a = I_ESP_I3x3z_S_a+ABX*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Px_a = I_ESP_I2x4y_S_a+ABX*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Px_a = I_ESP_I2x3yz_S_a+ABX*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Px_a = I_ESP_I2x2y2z_S_a+ABX*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Px_a = I_ESP_I2xy3z_S_a+ABX*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Px_a = I_ESP_I2x4z_S_a+ABX*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Px_a = I_ESP_Ix5y_S_a+ABX*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Px_a = I_ESP_Ix4yz_S_a+ABX*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Px_a = I_ESP_Ix3y2z_S_a+ABX*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Px_a = I_ESP_Ix2y3z_S_a+ABX*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Px_a = I_ESP_Ixy4z_S_a+ABX*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Px_a = I_ESP_Ix5z_S_a+ABX*I_ESP_H5z_S_a;
    Double I_ESP_H5x_Py_a = I_ESP_I5xy_S_a+ABY*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Py_a = I_ESP_I4x2y_S_a+ABY*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Py_a = I_ESP_I4xyz_S_a+ABY*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Py_a = I_ESP_I3x3y_S_a+ABY*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Py_a = I_ESP_I3x2yz_S_a+ABY*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Py_a = I_ESP_I3xy2z_S_a+ABY*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Py_a = I_ESP_I2x4y_S_a+ABY*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Py_a = I_ESP_I2x3yz_S_a+ABY*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Py_a = I_ESP_I2x2y2z_S_a+ABY*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Py_a = I_ESP_I2xy3z_S_a+ABY*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Py_a = I_ESP_Ix5y_S_a+ABY*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Py_a = I_ESP_Ix4yz_S_a+ABY*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Py_a = I_ESP_Ix3y2z_S_a+ABY*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Py_a = I_ESP_Ix2y3z_S_a+ABY*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Py_a = I_ESP_Ixy4z_S_a+ABY*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Py_a = I_ESP_I6y_S_a+ABY*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Py_a = I_ESP_I5yz_S_a+ABY*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Py_a = I_ESP_I4y2z_S_a+ABY*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Py_a = I_ESP_I3y3z_S_a+ABY*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Py_a = I_ESP_I2y4z_S_a+ABY*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Py_a = I_ESP_Iy5z_S_a+ABY*I_ESP_H5z_S_a;
    Double I_ESP_H5x_Pz_a = I_ESP_I5xz_S_a+ABZ*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Pz_a = I_ESP_I4xyz_S_a+ABZ*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Pz_a = I_ESP_I4x2z_S_a+ABZ*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Pz_a = I_ESP_I3x2yz_S_a+ABZ*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Pz_a = I_ESP_I3xy2z_S_a+ABZ*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Pz_a = I_ESP_I3x3z_S_a+ABZ*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Pz_a = I_ESP_I2x3yz_S_a+ABZ*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Pz_a = I_ESP_I2x2y2z_S_a+ABZ*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Pz_a = I_ESP_I2xy3z_S_a+ABZ*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Pz_a = I_ESP_I2x4z_S_a+ABZ*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Pz_a = I_ESP_Ix4yz_S_a+ABZ*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Pz_a = I_ESP_Ix3y2z_S_a+ABZ*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Pz_a = I_ESP_Ix2y3z_S_a+ABZ*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Pz_a = I_ESP_Ixy4z_S_a+ABZ*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Pz_a = I_ESP_Ix5z_S_a+ABZ*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Pz_a = I_ESP_I5yz_S_a+ABZ*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Pz_a = I_ESP_I4y2z_S_a+ABZ*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Pz_a = I_ESP_I3y3z_S_a+ABZ*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Pz_a = I_ESP_I2y4z_S_a+ABZ*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Pz_a = I_ESP_Iy5z_S_a+ABZ*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Pz_a = I_ESP_I6z_S_a+ABZ*I_ESP_H5z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S_aa
     * RHS shell quartet name: SQ_ESP_K_S_aa
     ************************************************************/
    Double I_ESP_K7x_Px_aa = I_ESP_L8x_S_aa+ABX*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Px_aa = I_ESP_L7xy_S_aa+ABX*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Px_aa = I_ESP_L7xz_S_aa+ABX*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Px_aa = I_ESP_L6x2y_S_aa+ABX*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Px_aa = I_ESP_L6xyz_S_aa+ABX*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Px_aa = I_ESP_L6x2z_S_aa+ABX*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Px_aa = I_ESP_L5x3y_S_aa+ABX*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Px_aa = I_ESP_L5x2yz_S_aa+ABX*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Px_aa = I_ESP_L5xy2z_S_aa+ABX*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Px_aa = I_ESP_L5x3z_S_aa+ABX*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Px_aa = I_ESP_L4x4y_S_aa+ABX*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Px_aa = I_ESP_L4x3yz_S_aa+ABX*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Px_aa = I_ESP_L4x2y2z_S_aa+ABX*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Px_aa = I_ESP_L4xy3z_S_aa+ABX*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Px_aa = I_ESP_L4x4z_S_aa+ABX*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Px_aa = I_ESP_L3x5y_S_aa+ABX*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Px_aa = I_ESP_L3x4yz_S_aa+ABX*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Px_aa = I_ESP_L3x3y2z_S_aa+ABX*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Px_aa = I_ESP_L3x2y3z_S_aa+ABX*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Px_aa = I_ESP_L3xy4z_S_aa+ABX*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Px_aa = I_ESP_L3x5z_S_aa+ABX*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Px_aa = I_ESP_L2x6y_S_aa+ABX*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Px_aa = I_ESP_L2x5yz_S_aa+ABX*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Px_aa = I_ESP_L2x4y2z_S_aa+ABX*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Px_aa = I_ESP_L2x3y3z_S_aa+ABX*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Px_aa = I_ESP_L2x2y4z_S_aa+ABX*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Px_aa = I_ESP_L2xy5z_S_aa+ABX*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Px_aa = I_ESP_L2x6z_S_aa+ABX*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Px_aa = I_ESP_Lx7y_S_aa+ABX*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Px_aa = I_ESP_Lx6yz_S_aa+ABX*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Px_aa = I_ESP_Lx5y2z_S_aa+ABX*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Px_aa = I_ESP_Lx4y3z_S_aa+ABX*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Px_aa = I_ESP_Lx3y4z_S_aa+ABX*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Px_aa = I_ESP_Lx2y5z_S_aa+ABX*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Px_aa = I_ESP_Lxy6z_S_aa+ABX*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Px_aa = I_ESP_Lx7z_S_aa+ABX*I_ESP_K7z_S_aa;
    Double I_ESP_K7x_Py_aa = I_ESP_L7xy_S_aa+ABY*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Py_aa = I_ESP_L6x2y_S_aa+ABY*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Py_aa = I_ESP_L6xyz_S_aa+ABY*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Py_aa = I_ESP_L5x3y_S_aa+ABY*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Py_aa = I_ESP_L5x2yz_S_aa+ABY*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Py_aa = I_ESP_L5xy2z_S_aa+ABY*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Py_aa = I_ESP_L4x4y_S_aa+ABY*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Py_aa = I_ESP_L4x3yz_S_aa+ABY*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Py_aa = I_ESP_L4x2y2z_S_aa+ABY*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Py_aa = I_ESP_L4xy3z_S_aa+ABY*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Py_aa = I_ESP_L3x5y_S_aa+ABY*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Py_aa = I_ESP_L3x4yz_S_aa+ABY*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Py_aa = I_ESP_L3x3y2z_S_aa+ABY*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Py_aa = I_ESP_L3x2y3z_S_aa+ABY*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Py_aa = I_ESP_L3xy4z_S_aa+ABY*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Py_aa = I_ESP_L2x6y_S_aa+ABY*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Py_aa = I_ESP_L2x5yz_S_aa+ABY*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Py_aa = I_ESP_L2x4y2z_S_aa+ABY*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Py_aa = I_ESP_L2x3y3z_S_aa+ABY*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Py_aa = I_ESP_L2x2y4z_S_aa+ABY*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Py_aa = I_ESP_L2xy5z_S_aa+ABY*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Py_aa = I_ESP_Lx7y_S_aa+ABY*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Py_aa = I_ESP_Lx6yz_S_aa+ABY*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Py_aa = I_ESP_Lx5y2z_S_aa+ABY*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Py_aa = I_ESP_Lx4y3z_S_aa+ABY*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Py_aa = I_ESP_Lx3y4z_S_aa+ABY*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Py_aa = I_ESP_Lx2y5z_S_aa+ABY*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Py_aa = I_ESP_Lxy6z_S_aa+ABY*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Py_aa = I_ESP_L8y_S_aa+ABY*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Py_aa = I_ESP_L7yz_S_aa+ABY*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Py_aa = I_ESP_L6y2z_S_aa+ABY*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Py_aa = I_ESP_L5y3z_S_aa+ABY*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Py_aa = I_ESP_L4y4z_S_aa+ABY*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Py_aa = I_ESP_L3y5z_S_aa+ABY*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Py_aa = I_ESP_L2y6z_S_aa+ABY*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Py_aa = I_ESP_Ly7z_S_aa+ABY*I_ESP_K7z_S_aa;
    Double I_ESP_K7x_Pz_aa = I_ESP_L7xz_S_aa+ABZ*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Pz_aa = I_ESP_L6xyz_S_aa+ABZ*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Pz_aa = I_ESP_L6x2z_S_aa+ABZ*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Pz_aa = I_ESP_L5x2yz_S_aa+ABZ*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Pz_aa = I_ESP_L5xy2z_S_aa+ABZ*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Pz_aa = I_ESP_L5x3z_S_aa+ABZ*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Pz_aa = I_ESP_L4x3yz_S_aa+ABZ*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Pz_aa = I_ESP_L4x2y2z_S_aa+ABZ*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Pz_aa = I_ESP_L4xy3z_S_aa+ABZ*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Pz_aa = I_ESP_L4x4z_S_aa+ABZ*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Pz_aa = I_ESP_L3x4yz_S_aa+ABZ*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Pz_aa = I_ESP_L3x3y2z_S_aa+ABZ*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Pz_aa = I_ESP_L3x2y3z_S_aa+ABZ*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Pz_aa = I_ESP_L3xy4z_S_aa+ABZ*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Pz_aa = I_ESP_L3x5z_S_aa+ABZ*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Pz_aa = I_ESP_L2x5yz_S_aa+ABZ*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Pz_aa = I_ESP_L2x4y2z_S_aa+ABZ*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Pz_aa = I_ESP_L2x3y3z_S_aa+ABZ*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Pz_aa = I_ESP_L2x2y4z_S_aa+ABZ*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Pz_aa = I_ESP_L2xy5z_S_aa+ABZ*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Pz_aa = I_ESP_L2x6z_S_aa+ABZ*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Pz_aa = I_ESP_Lx6yz_S_aa+ABZ*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Pz_aa = I_ESP_Lx5y2z_S_aa+ABZ*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Pz_aa = I_ESP_Lx4y3z_S_aa+ABZ*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Pz_aa = I_ESP_Lx3y4z_S_aa+ABZ*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Pz_aa = I_ESP_Lx2y5z_S_aa+ABZ*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Pz_aa = I_ESP_Lxy6z_S_aa+ABZ*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Pz_aa = I_ESP_Lx7z_S_aa+ABZ*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Pz_aa = I_ESP_L7yz_S_aa+ABZ*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Pz_aa = I_ESP_L6y2z_S_aa+ABZ*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Pz_aa = I_ESP_L5y3z_S_aa+ABZ*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Pz_aa = I_ESP_L4y4z_S_aa+ABZ*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Pz_aa = I_ESP_L3y5z_S_aa+ABZ*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Pz_aa = I_ESP_L2y6z_S_aa+ABZ*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Pz_aa = I_ESP_Ly7z_S_aa+ABZ*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Pz_aa = I_ESP_L8z_S_aa+ABZ*I_ESP_K7z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+0] = 4.0E0*I_ESP_K7x_Px_aa-2.0E0*5*I_ESP_H5x_Px_a-2.0E0*6*I_ESP_H5x_Px_a+5*4*I_ESP_F3x_Px;
    abcd[iGrid*378+1] = 4.0E0*I_ESP_K6xy_Px_aa-2.0E0*4*I_ESP_H4xy_Px_a-2.0E0*5*I_ESP_H4xy_Px_a+4*3*I_ESP_F2xy_Px;
    abcd[iGrid*378+2] = 4.0E0*I_ESP_K6xz_Px_aa-2.0E0*4*I_ESP_H4xz_Px_a-2.0E0*5*I_ESP_H4xz_Px_a+4*3*I_ESP_F2xz_Px;
    abcd[iGrid*378+3] = 4.0E0*I_ESP_K5x2y_Px_aa-2.0E0*3*I_ESP_H3x2y_Px_a-2.0E0*4*I_ESP_H3x2y_Px_a+3*2*I_ESP_Fx2y_Px;
    abcd[iGrid*378+4] = 4.0E0*I_ESP_K5xyz_Px_aa-2.0E0*3*I_ESP_H3xyz_Px_a-2.0E0*4*I_ESP_H3xyz_Px_a+3*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+5] = 4.0E0*I_ESP_K5x2z_Px_aa-2.0E0*3*I_ESP_H3x2z_Px_a-2.0E0*4*I_ESP_H3x2z_Px_a+3*2*I_ESP_Fx2z_Px;
    abcd[iGrid*378+6] = 4.0E0*I_ESP_K4x3y_Px_aa-2.0E0*2*I_ESP_H2x3y_Px_a-2.0E0*3*I_ESP_H2x3y_Px_a+2*1*I_ESP_F3y_Px;
    abcd[iGrid*378+7] = 4.0E0*I_ESP_K4x2yz_Px_aa-2.0E0*2*I_ESP_H2x2yz_Px_a-2.0E0*3*I_ESP_H2x2yz_Px_a+2*1*I_ESP_F2yz_Px;
    abcd[iGrid*378+8] = 4.0E0*I_ESP_K4xy2z_Px_aa-2.0E0*2*I_ESP_H2xy2z_Px_a-2.0E0*3*I_ESP_H2xy2z_Px_a+2*1*I_ESP_Fy2z_Px;
    abcd[iGrid*378+9] = 4.0E0*I_ESP_K4x3z_Px_aa-2.0E0*2*I_ESP_H2x3z_Px_a-2.0E0*3*I_ESP_H2x3z_Px_a+2*1*I_ESP_F3z_Px;
    abcd[iGrid*378+10] = 4.0E0*I_ESP_K3x4y_Px_aa-2.0E0*1*I_ESP_Hx4y_Px_a-2.0E0*2*I_ESP_Hx4y_Px_a;
    abcd[iGrid*378+11] = 4.0E0*I_ESP_K3x3yz_Px_aa-2.0E0*1*I_ESP_Hx3yz_Px_a-2.0E0*2*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+12] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_Hx2y2z_Px_a-2.0E0*2*I_ESP_Hx2y2z_Px_a;
    abcd[iGrid*378+13] = 4.0E0*I_ESP_K3xy3z_Px_aa-2.0E0*1*I_ESP_Hxy3z_Px_a-2.0E0*2*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+14] = 4.0E0*I_ESP_K3x4z_Px_aa-2.0E0*1*I_ESP_Hx4z_Px_a-2.0E0*2*I_ESP_Hx4z_Px_a;
    abcd[iGrid*378+15] = 4.0E0*I_ESP_K2x5y_Px_aa-2.0E0*1*I_ESP_H5y_Px_a;
    abcd[iGrid*378+16] = 4.0E0*I_ESP_K2x4yz_Px_aa-2.0E0*1*I_ESP_H4yz_Px_a;
    abcd[iGrid*378+17] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*1*I_ESP_H3y2z_Px_a;
    abcd[iGrid*378+18] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*1*I_ESP_H2y3z_Px_a;
    abcd[iGrid*378+19] = 4.0E0*I_ESP_K2xy4z_Px_aa-2.0E0*1*I_ESP_Hy4z_Px_a;
    abcd[iGrid*378+20] = 4.0E0*I_ESP_K2x5z_Px_aa-2.0E0*1*I_ESP_H5z_Px_a;
    abcd[iGrid*378+21] = 4.0E0*I_ESP_K7x_Py_aa-2.0E0*5*I_ESP_H5x_Py_a-2.0E0*6*I_ESP_H5x_Py_a+5*4*I_ESP_F3x_Py;
    abcd[iGrid*378+22] = 4.0E0*I_ESP_K6xy_Py_aa-2.0E0*4*I_ESP_H4xy_Py_a-2.0E0*5*I_ESP_H4xy_Py_a+4*3*I_ESP_F2xy_Py;
    abcd[iGrid*378+23] = 4.0E0*I_ESP_K6xz_Py_aa-2.0E0*4*I_ESP_H4xz_Py_a-2.0E0*5*I_ESP_H4xz_Py_a+4*3*I_ESP_F2xz_Py;
    abcd[iGrid*378+24] = 4.0E0*I_ESP_K5x2y_Py_aa-2.0E0*3*I_ESP_H3x2y_Py_a-2.0E0*4*I_ESP_H3x2y_Py_a+3*2*I_ESP_Fx2y_Py;
    abcd[iGrid*378+25] = 4.0E0*I_ESP_K5xyz_Py_aa-2.0E0*3*I_ESP_H3xyz_Py_a-2.0E0*4*I_ESP_H3xyz_Py_a+3*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+26] = 4.0E0*I_ESP_K5x2z_Py_aa-2.0E0*3*I_ESP_H3x2z_Py_a-2.0E0*4*I_ESP_H3x2z_Py_a+3*2*I_ESP_Fx2z_Py;
    abcd[iGrid*378+27] = 4.0E0*I_ESP_K4x3y_Py_aa-2.0E0*2*I_ESP_H2x3y_Py_a-2.0E0*3*I_ESP_H2x3y_Py_a+2*1*I_ESP_F3y_Py;
    abcd[iGrid*378+28] = 4.0E0*I_ESP_K4x2yz_Py_aa-2.0E0*2*I_ESP_H2x2yz_Py_a-2.0E0*3*I_ESP_H2x2yz_Py_a+2*1*I_ESP_F2yz_Py;
    abcd[iGrid*378+29] = 4.0E0*I_ESP_K4xy2z_Py_aa-2.0E0*2*I_ESP_H2xy2z_Py_a-2.0E0*3*I_ESP_H2xy2z_Py_a+2*1*I_ESP_Fy2z_Py;
    abcd[iGrid*378+30] = 4.0E0*I_ESP_K4x3z_Py_aa-2.0E0*2*I_ESP_H2x3z_Py_a-2.0E0*3*I_ESP_H2x3z_Py_a+2*1*I_ESP_F3z_Py;
    abcd[iGrid*378+31] = 4.0E0*I_ESP_K3x4y_Py_aa-2.0E0*1*I_ESP_Hx4y_Py_a-2.0E0*2*I_ESP_Hx4y_Py_a;
    abcd[iGrid*378+32] = 4.0E0*I_ESP_K3x3yz_Py_aa-2.0E0*1*I_ESP_Hx3yz_Py_a-2.0E0*2*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+33] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_Hx2y2z_Py_a-2.0E0*2*I_ESP_Hx2y2z_Py_a;
    abcd[iGrid*378+34] = 4.0E0*I_ESP_K3xy3z_Py_aa-2.0E0*1*I_ESP_Hxy3z_Py_a-2.0E0*2*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+35] = 4.0E0*I_ESP_K3x4z_Py_aa-2.0E0*1*I_ESP_Hx4z_Py_a-2.0E0*2*I_ESP_Hx4z_Py_a;
    abcd[iGrid*378+36] = 4.0E0*I_ESP_K2x5y_Py_aa-2.0E0*1*I_ESP_H5y_Py_a;
    abcd[iGrid*378+37] = 4.0E0*I_ESP_K2x4yz_Py_aa-2.0E0*1*I_ESP_H4yz_Py_a;
    abcd[iGrid*378+38] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*1*I_ESP_H3y2z_Py_a;
    abcd[iGrid*378+39] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*1*I_ESP_H2y3z_Py_a;
    abcd[iGrid*378+40] = 4.0E0*I_ESP_K2xy4z_Py_aa-2.0E0*1*I_ESP_Hy4z_Py_a;
    abcd[iGrid*378+41] = 4.0E0*I_ESP_K2x5z_Py_aa-2.0E0*1*I_ESP_H5z_Py_a;
    abcd[iGrid*378+42] = 4.0E0*I_ESP_K7x_Pz_aa-2.0E0*5*I_ESP_H5x_Pz_a-2.0E0*6*I_ESP_H5x_Pz_a+5*4*I_ESP_F3x_Pz;
    abcd[iGrid*378+43] = 4.0E0*I_ESP_K6xy_Pz_aa-2.0E0*4*I_ESP_H4xy_Pz_a-2.0E0*5*I_ESP_H4xy_Pz_a+4*3*I_ESP_F2xy_Pz;
    abcd[iGrid*378+44] = 4.0E0*I_ESP_K6xz_Pz_aa-2.0E0*4*I_ESP_H4xz_Pz_a-2.0E0*5*I_ESP_H4xz_Pz_a+4*3*I_ESP_F2xz_Pz;
    abcd[iGrid*378+45] = 4.0E0*I_ESP_K5x2y_Pz_aa-2.0E0*3*I_ESP_H3x2y_Pz_a-2.0E0*4*I_ESP_H3x2y_Pz_a+3*2*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+46] = 4.0E0*I_ESP_K5xyz_Pz_aa-2.0E0*3*I_ESP_H3xyz_Pz_a-2.0E0*4*I_ESP_H3xyz_Pz_a+3*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+47] = 4.0E0*I_ESP_K5x2z_Pz_aa-2.0E0*3*I_ESP_H3x2z_Pz_a-2.0E0*4*I_ESP_H3x2z_Pz_a+3*2*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+48] = 4.0E0*I_ESP_K4x3y_Pz_aa-2.0E0*2*I_ESP_H2x3y_Pz_a-2.0E0*3*I_ESP_H2x3y_Pz_a+2*1*I_ESP_F3y_Pz;
    abcd[iGrid*378+49] = 4.0E0*I_ESP_K4x2yz_Pz_aa-2.0E0*2*I_ESP_H2x2yz_Pz_a-2.0E0*3*I_ESP_H2x2yz_Pz_a+2*1*I_ESP_F2yz_Pz;
    abcd[iGrid*378+50] = 4.0E0*I_ESP_K4xy2z_Pz_aa-2.0E0*2*I_ESP_H2xy2z_Pz_a-2.0E0*3*I_ESP_H2xy2z_Pz_a+2*1*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+51] = 4.0E0*I_ESP_K4x3z_Pz_aa-2.0E0*2*I_ESP_H2x3z_Pz_a-2.0E0*3*I_ESP_H2x3z_Pz_a+2*1*I_ESP_F3z_Pz;
    abcd[iGrid*378+52] = 4.0E0*I_ESP_K3x4y_Pz_aa-2.0E0*1*I_ESP_Hx4y_Pz_a-2.0E0*2*I_ESP_Hx4y_Pz_a;
    abcd[iGrid*378+53] = 4.0E0*I_ESP_K3x3yz_Pz_aa-2.0E0*1*I_ESP_Hx3yz_Pz_a-2.0E0*2*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+54] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_Hx2y2z_Pz_a-2.0E0*2*I_ESP_Hx2y2z_Pz_a;
    abcd[iGrid*378+55] = 4.0E0*I_ESP_K3xy3z_Pz_aa-2.0E0*1*I_ESP_Hxy3z_Pz_a-2.0E0*2*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+56] = 4.0E0*I_ESP_K3x4z_Pz_aa-2.0E0*1*I_ESP_Hx4z_Pz_a-2.0E0*2*I_ESP_Hx4z_Pz_a;
    abcd[iGrid*378+57] = 4.0E0*I_ESP_K2x5y_Pz_aa-2.0E0*1*I_ESP_H5y_Pz_a;
    abcd[iGrid*378+58] = 4.0E0*I_ESP_K2x4yz_Pz_aa-2.0E0*1*I_ESP_H4yz_Pz_a;
    abcd[iGrid*378+59] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*1*I_ESP_H3y2z_Pz_a;
    abcd[iGrid*378+60] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*1*I_ESP_H2y3z_Pz_a;
    abcd[iGrid*378+61] = 4.0E0*I_ESP_K2xy4z_Pz_aa-2.0E0*1*I_ESP_Hy4z_Pz_a;
    abcd[iGrid*378+62] = 4.0E0*I_ESP_K2x5z_Pz_aa-2.0E0*1*I_ESP_H5z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+63] = 4.0E0*I_ESP_K6xy_Px_aa-2.0E0*5*I_ESP_H4xy_Px_a;
    abcd[iGrid*378+64] = 4.0E0*I_ESP_K5x2y_Px_aa-2.0E0*1*I_ESP_H5x_Px_a-2.0E0*4*I_ESP_H3x2y_Px_a+4*1*I_ESP_F3x_Px;
    abcd[iGrid*378+65] = 4.0E0*I_ESP_K5xyz_Px_aa-2.0E0*4*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+66] = 4.0E0*I_ESP_K4x3y_Px_aa-2.0E0*2*I_ESP_H4xy_Px_a-2.0E0*3*I_ESP_H2x3y_Px_a+3*2*I_ESP_F2xy_Px;
    abcd[iGrid*378+67] = 4.0E0*I_ESP_K4x2yz_Px_aa-2.0E0*1*I_ESP_H4xz_Px_a-2.0E0*3*I_ESP_H2x2yz_Px_a+3*1*I_ESP_F2xz_Px;
    abcd[iGrid*378+68] = 4.0E0*I_ESP_K4xy2z_Px_aa-2.0E0*3*I_ESP_H2xy2z_Px_a;
    abcd[iGrid*378+69] = 4.0E0*I_ESP_K3x4y_Px_aa-2.0E0*3*I_ESP_H3x2y_Px_a-2.0E0*2*I_ESP_Hx4y_Px_a+2*3*I_ESP_Fx2y_Px;
    abcd[iGrid*378+70] = 4.0E0*I_ESP_K3x3yz_Px_aa-2.0E0*2*I_ESP_H3xyz_Px_a-2.0E0*2*I_ESP_Hx3yz_Px_a+2*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+71] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_H3x2z_Px_a-2.0E0*2*I_ESP_Hx2y2z_Px_a+2*1*I_ESP_Fx2z_Px;
    abcd[iGrid*378+72] = 4.0E0*I_ESP_K3xy3z_Px_aa-2.0E0*2*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+73] = 4.0E0*I_ESP_K2x5y_Px_aa-2.0E0*4*I_ESP_H2x3y_Px_a-2.0E0*1*I_ESP_H5y_Px_a+4*I_ESP_F3y_Px;
    abcd[iGrid*378+74] = 4.0E0*I_ESP_K2x4yz_Px_aa-2.0E0*3*I_ESP_H2x2yz_Px_a-2.0E0*1*I_ESP_H4yz_Px_a+3*I_ESP_F2yz_Px;
    abcd[iGrid*378+75] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*2*I_ESP_H2xy2z_Px_a-2.0E0*1*I_ESP_H3y2z_Px_a+2*I_ESP_Fy2z_Px;
    abcd[iGrid*378+76] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*1*I_ESP_H2x3z_Px_a-2.0E0*1*I_ESP_H2y3z_Px_a+1*I_ESP_F3z_Px;
    abcd[iGrid*378+77] = 4.0E0*I_ESP_K2xy4z_Px_aa-2.0E0*1*I_ESP_Hy4z_Px_a;
    abcd[iGrid*378+78] = 4.0E0*I_ESP_Kx6y_Px_aa-2.0E0*5*I_ESP_Hx4y_Px_a;
    abcd[iGrid*378+79] = 4.0E0*I_ESP_Kx5yz_Px_aa-2.0E0*4*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+80] = 4.0E0*I_ESP_Kx4y2z_Px_aa-2.0E0*3*I_ESP_Hx2y2z_Px_a;
    abcd[iGrid*378+81] = 4.0E0*I_ESP_Kx3y3z_Px_aa-2.0E0*2*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+82] = 4.0E0*I_ESP_Kx2y4z_Px_aa-2.0E0*1*I_ESP_Hx4z_Px_a;
    abcd[iGrid*378+83] = 4.0E0*I_ESP_Kxy5z_Px_aa;
    abcd[iGrid*378+84] = 4.0E0*I_ESP_K6xy_Py_aa-2.0E0*5*I_ESP_H4xy_Py_a;
    abcd[iGrid*378+85] = 4.0E0*I_ESP_K5x2y_Py_aa-2.0E0*1*I_ESP_H5x_Py_a-2.0E0*4*I_ESP_H3x2y_Py_a+4*1*I_ESP_F3x_Py;
    abcd[iGrid*378+86] = 4.0E0*I_ESP_K5xyz_Py_aa-2.0E0*4*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+87] = 4.0E0*I_ESP_K4x3y_Py_aa-2.0E0*2*I_ESP_H4xy_Py_a-2.0E0*3*I_ESP_H2x3y_Py_a+3*2*I_ESP_F2xy_Py;
    abcd[iGrid*378+88] = 4.0E0*I_ESP_K4x2yz_Py_aa-2.0E0*1*I_ESP_H4xz_Py_a-2.0E0*3*I_ESP_H2x2yz_Py_a+3*1*I_ESP_F2xz_Py;
    abcd[iGrid*378+89] = 4.0E0*I_ESP_K4xy2z_Py_aa-2.0E0*3*I_ESP_H2xy2z_Py_a;
    abcd[iGrid*378+90] = 4.0E0*I_ESP_K3x4y_Py_aa-2.0E0*3*I_ESP_H3x2y_Py_a-2.0E0*2*I_ESP_Hx4y_Py_a+2*3*I_ESP_Fx2y_Py;
    abcd[iGrid*378+91] = 4.0E0*I_ESP_K3x3yz_Py_aa-2.0E0*2*I_ESP_H3xyz_Py_a-2.0E0*2*I_ESP_Hx3yz_Py_a+2*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+92] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_H3x2z_Py_a-2.0E0*2*I_ESP_Hx2y2z_Py_a+2*1*I_ESP_Fx2z_Py;
    abcd[iGrid*378+93] = 4.0E0*I_ESP_K3xy3z_Py_aa-2.0E0*2*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+94] = 4.0E0*I_ESP_K2x5y_Py_aa-2.0E0*4*I_ESP_H2x3y_Py_a-2.0E0*1*I_ESP_H5y_Py_a+4*I_ESP_F3y_Py;
    abcd[iGrid*378+95] = 4.0E0*I_ESP_K2x4yz_Py_aa-2.0E0*3*I_ESP_H2x2yz_Py_a-2.0E0*1*I_ESP_H4yz_Py_a+3*I_ESP_F2yz_Py;
    abcd[iGrid*378+96] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*2*I_ESP_H2xy2z_Py_a-2.0E0*1*I_ESP_H3y2z_Py_a+2*I_ESP_Fy2z_Py;
    abcd[iGrid*378+97] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*1*I_ESP_H2x3z_Py_a-2.0E0*1*I_ESP_H2y3z_Py_a+1*I_ESP_F3z_Py;
    abcd[iGrid*378+98] = 4.0E0*I_ESP_K2xy4z_Py_aa-2.0E0*1*I_ESP_Hy4z_Py_a;
    abcd[iGrid*378+99] = 4.0E0*I_ESP_Kx6y_Py_aa-2.0E0*5*I_ESP_Hx4y_Py_a;
    abcd[iGrid*378+100] = 4.0E0*I_ESP_Kx5yz_Py_aa-2.0E0*4*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+101] = 4.0E0*I_ESP_Kx4y2z_Py_aa-2.0E0*3*I_ESP_Hx2y2z_Py_a;
    abcd[iGrid*378+102] = 4.0E0*I_ESP_Kx3y3z_Py_aa-2.0E0*2*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+103] = 4.0E0*I_ESP_Kx2y4z_Py_aa-2.0E0*1*I_ESP_Hx4z_Py_a;
    abcd[iGrid*378+104] = 4.0E0*I_ESP_Kxy5z_Py_aa;
    abcd[iGrid*378+105] = 4.0E0*I_ESP_K6xy_Pz_aa-2.0E0*5*I_ESP_H4xy_Pz_a;
    abcd[iGrid*378+106] = 4.0E0*I_ESP_K5x2y_Pz_aa-2.0E0*1*I_ESP_H5x_Pz_a-2.0E0*4*I_ESP_H3x2y_Pz_a+4*1*I_ESP_F3x_Pz;
    abcd[iGrid*378+107] = 4.0E0*I_ESP_K5xyz_Pz_aa-2.0E0*4*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+108] = 4.0E0*I_ESP_K4x3y_Pz_aa-2.0E0*2*I_ESP_H4xy_Pz_a-2.0E0*3*I_ESP_H2x3y_Pz_a+3*2*I_ESP_F2xy_Pz;
    abcd[iGrid*378+109] = 4.0E0*I_ESP_K4x2yz_Pz_aa-2.0E0*1*I_ESP_H4xz_Pz_a-2.0E0*3*I_ESP_H2x2yz_Pz_a+3*1*I_ESP_F2xz_Pz;
    abcd[iGrid*378+110] = 4.0E0*I_ESP_K4xy2z_Pz_aa-2.0E0*3*I_ESP_H2xy2z_Pz_a;
    abcd[iGrid*378+111] = 4.0E0*I_ESP_K3x4y_Pz_aa-2.0E0*3*I_ESP_H3x2y_Pz_a-2.0E0*2*I_ESP_Hx4y_Pz_a+2*3*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+112] = 4.0E0*I_ESP_K3x3yz_Pz_aa-2.0E0*2*I_ESP_H3xyz_Pz_a-2.0E0*2*I_ESP_Hx3yz_Pz_a+2*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+113] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_H3x2z_Pz_a-2.0E0*2*I_ESP_Hx2y2z_Pz_a+2*1*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+114] = 4.0E0*I_ESP_K3xy3z_Pz_aa-2.0E0*2*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+115] = 4.0E0*I_ESP_K2x5y_Pz_aa-2.0E0*4*I_ESP_H2x3y_Pz_a-2.0E0*1*I_ESP_H5y_Pz_a+4*I_ESP_F3y_Pz;
    abcd[iGrid*378+116] = 4.0E0*I_ESP_K2x4yz_Pz_aa-2.0E0*3*I_ESP_H2x2yz_Pz_a-2.0E0*1*I_ESP_H4yz_Pz_a+3*I_ESP_F2yz_Pz;
    abcd[iGrid*378+117] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*2*I_ESP_H2xy2z_Pz_a-2.0E0*1*I_ESP_H3y2z_Pz_a+2*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+118] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*1*I_ESP_H2x3z_Pz_a-2.0E0*1*I_ESP_H2y3z_Pz_a+1*I_ESP_F3z_Pz;
    abcd[iGrid*378+119] = 4.0E0*I_ESP_K2xy4z_Pz_aa-2.0E0*1*I_ESP_Hy4z_Pz_a;
    abcd[iGrid*378+120] = 4.0E0*I_ESP_Kx6y_Pz_aa-2.0E0*5*I_ESP_Hx4y_Pz_a;
    abcd[iGrid*378+121] = 4.0E0*I_ESP_Kx5yz_Pz_aa-2.0E0*4*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+122] = 4.0E0*I_ESP_Kx4y2z_Pz_aa-2.0E0*3*I_ESP_Hx2y2z_Pz_a;
    abcd[iGrid*378+123] = 4.0E0*I_ESP_Kx3y3z_Pz_aa-2.0E0*2*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+124] = 4.0E0*I_ESP_Kx2y4z_Pz_aa-2.0E0*1*I_ESP_Hx4z_Pz_a;
    abcd[iGrid*378+125] = 4.0E0*I_ESP_Kxy5z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+126] = 4.0E0*I_ESP_K6xz_Px_aa-2.0E0*5*I_ESP_H4xz_Px_a;
    abcd[iGrid*378+127] = 4.0E0*I_ESP_K5xyz_Px_aa-2.0E0*4*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+128] = 4.0E0*I_ESP_K5x2z_Px_aa-2.0E0*1*I_ESP_H5x_Px_a-2.0E0*4*I_ESP_H3x2z_Px_a+4*1*I_ESP_F3x_Px;
    abcd[iGrid*378+129] = 4.0E0*I_ESP_K4x2yz_Px_aa-2.0E0*3*I_ESP_H2x2yz_Px_a;
    abcd[iGrid*378+130] = 4.0E0*I_ESP_K4xy2z_Px_aa-2.0E0*1*I_ESP_H4xy_Px_a-2.0E0*3*I_ESP_H2xy2z_Px_a+3*1*I_ESP_F2xy_Px;
    abcd[iGrid*378+131] = 4.0E0*I_ESP_K4x3z_Px_aa-2.0E0*2*I_ESP_H4xz_Px_a-2.0E0*3*I_ESP_H2x3z_Px_a+3*2*I_ESP_F2xz_Px;
    abcd[iGrid*378+132] = 4.0E0*I_ESP_K3x3yz_Px_aa-2.0E0*2*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+133] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_H3x2y_Px_a-2.0E0*2*I_ESP_Hx2y2z_Px_a+2*1*I_ESP_Fx2y_Px;
    abcd[iGrid*378+134] = 4.0E0*I_ESP_K3xy3z_Px_aa-2.0E0*2*I_ESP_H3xyz_Px_a-2.0E0*2*I_ESP_Hxy3z_Px_a+2*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+135] = 4.0E0*I_ESP_K3x4z_Px_aa-2.0E0*3*I_ESP_H3x2z_Px_a-2.0E0*2*I_ESP_Hx4z_Px_a+2*3*I_ESP_Fx2z_Px;
    abcd[iGrid*378+136] = 4.0E0*I_ESP_K2x4yz_Px_aa-2.0E0*1*I_ESP_H4yz_Px_a;
    abcd[iGrid*378+137] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*1*I_ESP_H2x3y_Px_a-2.0E0*1*I_ESP_H3y2z_Px_a+1*I_ESP_F3y_Px;
    abcd[iGrid*378+138] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*2*I_ESP_H2x2yz_Px_a-2.0E0*1*I_ESP_H2y3z_Px_a+2*I_ESP_F2yz_Px;
    abcd[iGrid*378+139] = 4.0E0*I_ESP_K2xy4z_Px_aa-2.0E0*3*I_ESP_H2xy2z_Px_a-2.0E0*1*I_ESP_Hy4z_Px_a+3*I_ESP_Fy2z_Px;
    abcd[iGrid*378+140] = 4.0E0*I_ESP_K2x5z_Px_aa-2.0E0*4*I_ESP_H2x3z_Px_a-2.0E0*1*I_ESP_H5z_Px_a+4*I_ESP_F3z_Px;
    abcd[iGrid*378+141] = 4.0E0*I_ESP_Kx5yz_Px_aa;
    abcd[iGrid*378+142] = 4.0E0*I_ESP_Kx4y2z_Px_aa-2.0E0*1*I_ESP_Hx4y_Px_a;
    abcd[iGrid*378+143] = 4.0E0*I_ESP_Kx3y3z_Px_aa-2.0E0*2*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+144] = 4.0E0*I_ESP_Kx2y4z_Px_aa-2.0E0*3*I_ESP_Hx2y2z_Px_a;
    abcd[iGrid*378+145] = 4.0E0*I_ESP_Kxy5z_Px_aa-2.0E0*4*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+146] = 4.0E0*I_ESP_Kx6z_Px_aa-2.0E0*5*I_ESP_Hx4z_Px_a;
    abcd[iGrid*378+147] = 4.0E0*I_ESP_K6xz_Py_aa-2.0E0*5*I_ESP_H4xz_Py_a;
    abcd[iGrid*378+148] = 4.0E0*I_ESP_K5xyz_Py_aa-2.0E0*4*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+149] = 4.0E0*I_ESP_K5x2z_Py_aa-2.0E0*1*I_ESP_H5x_Py_a-2.0E0*4*I_ESP_H3x2z_Py_a+4*1*I_ESP_F3x_Py;
    abcd[iGrid*378+150] = 4.0E0*I_ESP_K4x2yz_Py_aa-2.0E0*3*I_ESP_H2x2yz_Py_a;
    abcd[iGrid*378+151] = 4.0E0*I_ESP_K4xy2z_Py_aa-2.0E0*1*I_ESP_H4xy_Py_a-2.0E0*3*I_ESP_H2xy2z_Py_a+3*1*I_ESP_F2xy_Py;
    abcd[iGrid*378+152] = 4.0E0*I_ESP_K4x3z_Py_aa-2.0E0*2*I_ESP_H4xz_Py_a-2.0E0*3*I_ESP_H2x3z_Py_a+3*2*I_ESP_F2xz_Py;
    abcd[iGrid*378+153] = 4.0E0*I_ESP_K3x3yz_Py_aa-2.0E0*2*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+154] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_H3x2y_Py_a-2.0E0*2*I_ESP_Hx2y2z_Py_a+2*1*I_ESP_Fx2y_Py;
    abcd[iGrid*378+155] = 4.0E0*I_ESP_K3xy3z_Py_aa-2.0E0*2*I_ESP_H3xyz_Py_a-2.0E0*2*I_ESP_Hxy3z_Py_a+2*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+156] = 4.0E0*I_ESP_K3x4z_Py_aa-2.0E0*3*I_ESP_H3x2z_Py_a-2.0E0*2*I_ESP_Hx4z_Py_a+2*3*I_ESP_Fx2z_Py;
    abcd[iGrid*378+157] = 4.0E0*I_ESP_K2x4yz_Py_aa-2.0E0*1*I_ESP_H4yz_Py_a;
    abcd[iGrid*378+158] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*1*I_ESP_H2x3y_Py_a-2.0E0*1*I_ESP_H3y2z_Py_a+1*I_ESP_F3y_Py;
    abcd[iGrid*378+159] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*2*I_ESP_H2x2yz_Py_a-2.0E0*1*I_ESP_H2y3z_Py_a+2*I_ESP_F2yz_Py;
    abcd[iGrid*378+160] = 4.0E0*I_ESP_K2xy4z_Py_aa-2.0E0*3*I_ESP_H2xy2z_Py_a-2.0E0*1*I_ESP_Hy4z_Py_a+3*I_ESP_Fy2z_Py;
    abcd[iGrid*378+161] = 4.0E0*I_ESP_K2x5z_Py_aa-2.0E0*4*I_ESP_H2x3z_Py_a-2.0E0*1*I_ESP_H5z_Py_a+4*I_ESP_F3z_Py;
    abcd[iGrid*378+162] = 4.0E0*I_ESP_Kx5yz_Py_aa;
    abcd[iGrid*378+163] = 4.0E0*I_ESP_Kx4y2z_Py_aa-2.0E0*1*I_ESP_Hx4y_Py_a;
    abcd[iGrid*378+164] = 4.0E0*I_ESP_Kx3y3z_Py_aa-2.0E0*2*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+165] = 4.0E0*I_ESP_Kx2y4z_Py_aa-2.0E0*3*I_ESP_Hx2y2z_Py_a;
    abcd[iGrid*378+166] = 4.0E0*I_ESP_Kxy5z_Py_aa-2.0E0*4*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+167] = 4.0E0*I_ESP_Kx6z_Py_aa-2.0E0*5*I_ESP_Hx4z_Py_a;
    abcd[iGrid*378+168] = 4.0E0*I_ESP_K6xz_Pz_aa-2.0E0*5*I_ESP_H4xz_Pz_a;
    abcd[iGrid*378+169] = 4.0E0*I_ESP_K5xyz_Pz_aa-2.0E0*4*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+170] = 4.0E0*I_ESP_K5x2z_Pz_aa-2.0E0*1*I_ESP_H5x_Pz_a-2.0E0*4*I_ESP_H3x2z_Pz_a+4*1*I_ESP_F3x_Pz;
    abcd[iGrid*378+171] = 4.0E0*I_ESP_K4x2yz_Pz_aa-2.0E0*3*I_ESP_H2x2yz_Pz_a;
    abcd[iGrid*378+172] = 4.0E0*I_ESP_K4xy2z_Pz_aa-2.0E0*1*I_ESP_H4xy_Pz_a-2.0E0*3*I_ESP_H2xy2z_Pz_a+3*1*I_ESP_F2xy_Pz;
    abcd[iGrid*378+173] = 4.0E0*I_ESP_K4x3z_Pz_aa-2.0E0*2*I_ESP_H4xz_Pz_a-2.0E0*3*I_ESP_H2x3z_Pz_a+3*2*I_ESP_F2xz_Pz;
    abcd[iGrid*378+174] = 4.0E0*I_ESP_K3x3yz_Pz_aa-2.0E0*2*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+175] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_H3x2y_Pz_a-2.0E0*2*I_ESP_Hx2y2z_Pz_a+2*1*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+176] = 4.0E0*I_ESP_K3xy3z_Pz_aa-2.0E0*2*I_ESP_H3xyz_Pz_a-2.0E0*2*I_ESP_Hxy3z_Pz_a+2*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+177] = 4.0E0*I_ESP_K3x4z_Pz_aa-2.0E0*3*I_ESP_H3x2z_Pz_a-2.0E0*2*I_ESP_Hx4z_Pz_a+2*3*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+178] = 4.0E0*I_ESP_K2x4yz_Pz_aa-2.0E0*1*I_ESP_H4yz_Pz_a;
    abcd[iGrid*378+179] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*1*I_ESP_H2x3y_Pz_a-2.0E0*1*I_ESP_H3y2z_Pz_a+1*I_ESP_F3y_Pz;
    abcd[iGrid*378+180] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*2*I_ESP_H2x2yz_Pz_a-2.0E0*1*I_ESP_H2y3z_Pz_a+2*I_ESP_F2yz_Pz;
    abcd[iGrid*378+181] = 4.0E0*I_ESP_K2xy4z_Pz_aa-2.0E0*3*I_ESP_H2xy2z_Pz_a-2.0E0*1*I_ESP_Hy4z_Pz_a+3*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+182] = 4.0E0*I_ESP_K2x5z_Pz_aa-2.0E0*4*I_ESP_H2x3z_Pz_a-2.0E0*1*I_ESP_H5z_Pz_a+4*I_ESP_F3z_Pz;
    abcd[iGrid*378+183] = 4.0E0*I_ESP_Kx5yz_Pz_aa;
    abcd[iGrid*378+184] = 4.0E0*I_ESP_Kx4y2z_Pz_aa-2.0E0*1*I_ESP_Hx4y_Pz_a;
    abcd[iGrid*378+185] = 4.0E0*I_ESP_Kx3y3z_Pz_aa-2.0E0*2*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+186] = 4.0E0*I_ESP_Kx2y4z_Pz_aa-2.0E0*3*I_ESP_Hx2y2z_Pz_a;
    abcd[iGrid*378+187] = 4.0E0*I_ESP_Kxy5z_Pz_aa-2.0E0*4*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+188] = 4.0E0*I_ESP_Kx6z_Pz_aa-2.0E0*5*I_ESP_Hx4z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+189] = 4.0E0*I_ESP_K5x2y_Px_aa-2.0E0*1*I_ESP_H5x_Px_a;
    abcd[iGrid*378+190] = 4.0E0*I_ESP_K4x3y_Px_aa-2.0E0*1*I_ESP_H4xy_Px_a-2.0E0*2*I_ESP_H4xy_Px_a;
    abcd[iGrid*378+191] = 4.0E0*I_ESP_K4x2yz_Px_aa-2.0E0*1*I_ESP_H4xz_Px_a;
    abcd[iGrid*378+192] = 4.0E0*I_ESP_K3x4y_Px_aa-2.0E0*2*I_ESP_H3x2y_Px_a-2.0E0*3*I_ESP_H3x2y_Px_a+2*1*I_ESP_F3x_Px;
    abcd[iGrid*378+193] = 4.0E0*I_ESP_K3x3yz_Px_aa-2.0E0*1*I_ESP_H3xyz_Px_a-2.0E0*2*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+194] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_H3x2z_Px_a;
    abcd[iGrid*378+195] = 4.0E0*I_ESP_K2x5y_Px_aa-2.0E0*3*I_ESP_H2x3y_Px_a-2.0E0*4*I_ESP_H2x3y_Px_a+3*2*I_ESP_F2xy_Px;
    abcd[iGrid*378+196] = 4.0E0*I_ESP_K2x4yz_Px_aa-2.0E0*2*I_ESP_H2x2yz_Px_a-2.0E0*3*I_ESP_H2x2yz_Px_a+2*1*I_ESP_F2xz_Px;
    abcd[iGrid*378+197] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*1*I_ESP_H2xy2z_Px_a-2.0E0*2*I_ESP_H2xy2z_Px_a;
    abcd[iGrid*378+198] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*1*I_ESP_H2x3z_Px_a;
    abcd[iGrid*378+199] = 4.0E0*I_ESP_Kx6y_Px_aa-2.0E0*4*I_ESP_Hx4y_Px_a-2.0E0*5*I_ESP_Hx4y_Px_a+4*3*I_ESP_Fx2y_Px;
    abcd[iGrid*378+200] = 4.0E0*I_ESP_Kx5yz_Px_aa-2.0E0*3*I_ESP_Hx3yz_Px_a-2.0E0*4*I_ESP_Hx3yz_Px_a+3*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+201] = 4.0E0*I_ESP_Kx4y2z_Px_aa-2.0E0*2*I_ESP_Hx2y2z_Px_a-2.0E0*3*I_ESP_Hx2y2z_Px_a+2*1*I_ESP_Fx2z_Px;
    abcd[iGrid*378+202] = 4.0E0*I_ESP_Kx3y3z_Px_aa-2.0E0*1*I_ESP_Hxy3z_Px_a-2.0E0*2*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+203] = 4.0E0*I_ESP_Kx2y4z_Px_aa-2.0E0*1*I_ESP_Hx4z_Px_a;
    abcd[iGrid*378+204] = 4.0E0*I_ESP_K7y_Px_aa-2.0E0*5*I_ESP_H5y_Px_a-2.0E0*6*I_ESP_H5y_Px_a+5*4*I_ESP_F3y_Px;
    abcd[iGrid*378+205] = 4.0E0*I_ESP_K6yz_Px_aa-2.0E0*4*I_ESP_H4yz_Px_a-2.0E0*5*I_ESP_H4yz_Px_a+4*3*I_ESP_F2yz_Px;
    abcd[iGrid*378+206] = 4.0E0*I_ESP_K5y2z_Px_aa-2.0E0*3*I_ESP_H3y2z_Px_a-2.0E0*4*I_ESP_H3y2z_Px_a+3*2*I_ESP_Fy2z_Px;
    abcd[iGrid*378+207] = 4.0E0*I_ESP_K4y3z_Px_aa-2.0E0*2*I_ESP_H2y3z_Px_a-2.0E0*3*I_ESP_H2y3z_Px_a+2*1*I_ESP_F3z_Px;
    abcd[iGrid*378+208] = 4.0E0*I_ESP_K3y4z_Px_aa-2.0E0*1*I_ESP_Hy4z_Px_a-2.0E0*2*I_ESP_Hy4z_Px_a;
    abcd[iGrid*378+209] = 4.0E0*I_ESP_K2y5z_Px_aa-2.0E0*1*I_ESP_H5z_Px_a;
    abcd[iGrid*378+210] = 4.0E0*I_ESP_K5x2y_Py_aa-2.0E0*1*I_ESP_H5x_Py_a;
    abcd[iGrid*378+211] = 4.0E0*I_ESP_K4x3y_Py_aa-2.0E0*1*I_ESP_H4xy_Py_a-2.0E0*2*I_ESP_H4xy_Py_a;
    abcd[iGrid*378+212] = 4.0E0*I_ESP_K4x2yz_Py_aa-2.0E0*1*I_ESP_H4xz_Py_a;
    abcd[iGrid*378+213] = 4.0E0*I_ESP_K3x4y_Py_aa-2.0E0*2*I_ESP_H3x2y_Py_a-2.0E0*3*I_ESP_H3x2y_Py_a+2*1*I_ESP_F3x_Py;
    abcd[iGrid*378+214] = 4.0E0*I_ESP_K3x3yz_Py_aa-2.0E0*1*I_ESP_H3xyz_Py_a-2.0E0*2*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+215] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_H3x2z_Py_a;
    abcd[iGrid*378+216] = 4.0E0*I_ESP_K2x5y_Py_aa-2.0E0*3*I_ESP_H2x3y_Py_a-2.0E0*4*I_ESP_H2x3y_Py_a+3*2*I_ESP_F2xy_Py;
    abcd[iGrid*378+217] = 4.0E0*I_ESP_K2x4yz_Py_aa-2.0E0*2*I_ESP_H2x2yz_Py_a-2.0E0*3*I_ESP_H2x2yz_Py_a+2*1*I_ESP_F2xz_Py;
    abcd[iGrid*378+218] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*1*I_ESP_H2xy2z_Py_a-2.0E0*2*I_ESP_H2xy2z_Py_a;
    abcd[iGrid*378+219] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*1*I_ESP_H2x3z_Py_a;
    abcd[iGrid*378+220] = 4.0E0*I_ESP_Kx6y_Py_aa-2.0E0*4*I_ESP_Hx4y_Py_a-2.0E0*5*I_ESP_Hx4y_Py_a+4*3*I_ESP_Fx2y_Py;
    abcd[iGrid*378+221] = 4.0E0*I_ESP_Kx5yz_Py_aa-2.0E0*3*I_ESP_Hx3yz_Py_a-2.0E0*4*I_ESP_Hx3yz_Py_a+3*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+222] = 4.0E0*I_ESP_Kx4y2z_Py_aa-2.0E0*2*I_ESP_Hx2y2z_Py_a-2.0E0*3*I_ESP_Hx2y2z_Py_a+2*1*I_ESP_Fx2z_Py;
    abcd[iGrid*378+223] = 4.0E0*I_ESP_Kx3y3z_Py_aa-2.0E0*1*I_ESP_Hxy3z_Py_a-2.0E0*2*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+224] = 4.0E0*I_ESP_Kx2y4z_Py_aa-2.0E0*1*I_ESP_Hx4z_Py_a;
    abcd[iGrid*378+225] = 4.0E0*I_ESP_K7y_Py_aa-2.0E0*5*I_ESP_H5y_Py_a-2.0E0*6*I_ESP_H5y_Py_a+5*4*I_ESP_F3y_Py;
    abcd[iGrid*378+226] = 4.0E0*I_ESP_K6yz_Py_aa-2.0E0*4*I_ESP_H4yz_Py_a-2.0E0*5*I_ESP_H4yz_Py_a+4*3*I_ESP_F2yz_Py;
    abcd[iGrid*378+227] = 4.0E0*I_ESP_K5y2z_Py_aa-2.0E0*3*I_ESP_H3y2z_Py_a-2.0E0*4*I_ESP_H3y2z_Py_a+3*2*I_ESP_Fy2z_Py;
    abcd[iGrid*378+228] = 4.0E0*I_ESP_K4y3z_Py_aa-2.0E0*2*I_ESP_H2y3z_Py_a-2.0E0*3*I_ESP_H2y3z_Py_a+2*1*I_ESP_F3z_Py;
    abcd[iGrid*378+229] = 4.0E0*I_ESP_K3y4z_Py_aa-2.0E0*1*I_ESP_Hy4z_Py_a-2.0E0*2*I_ESP_Hy4z_Py_a;
    abcd[iGrid*378+230] = 4.0E0*I_ESP_K2y5z_Py_aa-2.0E0*1*I_ESP_H5z_Py_a;
    abcd[iGrid*378+231] = 4.0E0*I_ESP_K5x2y_Pz_aa-2.0E0*1*I_ESP_H5x_Pz_a;
    abcd[iGrid*378+232] = 4.0E0*I_ESP_K4x3y_Pz_aa-2.0E0*1*I_ESP_H4xy_Pz_a-2.0E0*2*I_ESP_H4xy_Pz_a;
    abcd[iGrid*378+233] = 4.0E0*I_ESP_K4x2yz_Pz_aa-2.0E0*1*I_ESP_H4xz_Pz_a;
    abcd[iGrid*378+234] = 4.0E0*I_ESP_K3x4y_Pz_aa-2.0E0*2*I_ESP_H3x2y_Pz_a-2.0E0*3*I_ESP_H3x2y_Pz_a+2*1*I_ESP_F3x_Pz;
    abcd[iGrid*378+235] = 4.0E0*I_ESP_K3x3yz_Pz_aa-2.0E0*1*I_ESP_H3xyz_Pz_a-2.0E0*2*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+236] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_H3x2z_Pz_a;
    abcd[iGrid*378+237] = 4.0E0*I_ESP_K2x5y_Pz_aa-2.0E0*3*I_ESP_H2x3y_Pz_a-2.0E0*4*I_ESP_H2x3y_Pz_a+3*2*I_ESP_F2xy_Pz;
    abcd[iGrid*378+238] = 4.0E0*I_ESP_K2x4yz_Pz_aa-2.0E0*2*I_ESP_H2x2yz_Pz_a-2.0E0*3*I_ESP_H2x2yz_Pz_a+2*1*I_ESP_F2xz_Pz;
    abcd[iGrid*378+239] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*1*I_ESP_H2xy2z_Pz_a-2.0E0*2*I_ESP_H2xy2z_Pz_a;
    abcd[iGrid*378+240] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*1*I_ESP_H2x3z_Pz_a;
    abcd[iGrid*378+241] = 4.0E0*I_ESP_Kx6y_Pz_aa-2.0E0*4*I_ESP_Hx4y_Pz_a-2.0E0*5*I_ESP_Hx4y_Pz_a+4*3*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+242] = 4.0E0*I_ESP_Kx5yz_Pz_aa-2.0E0*3*I_ESP_Hx3yz_Pz_a-2.0E0*4*I_ESP_Hx3yz_Pz_a+3*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+243] = 4.0E0*I_ESP_Kx4y2z_Pz_aa-2.0E0*2*I_ESP_Hx2y2z_Pz_a-2.0E0*3*I_ESP_Hx2y2z_Pz_a+2*1*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+244] = 4.0E0*I_ESP_Kx3y3z_Pz_aa-2.0E0*1*I_ESP_Hxy3z_Pz_a-2.0E0*2*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+245] = 4.0E0*I_ESP_Kx2y4z_Pz_aa-2.0E0*1*I_ESP_Hx4z_Pz_a;
    abcd[iGrid*378+246] = 4.0E0*I_ESP_K7y_Pz_aa-2.0E0*5*I_ESP_H5y_Pz_a-2.0E0*6*I_ESP_H5y_Pz_a+5*4*I_ESP_F3y_Pz;
    abcd[iGrid*378+247] = 4.0E0*I_ESP_K6yz_Pz_aa-2.0E0*4*I_ESP_H4yz_Pz_a-2.0E0*5*I_ESP_H4yz_Pz_a+4*3*I_ESP_F2yz_Pz;
    abcd[iGrid*378+248] = 4.0E0*I_ESP_K5y2z_Pz_aa-2.0E0*3*I_ESP_H3y2z_Pz_a-2.0E0*4*I_ESP_H3y2z_Pz_a+3*2*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+249] = 4.0E0*I_ESP_K4y3z_Pz_aa-2.0E0*2*I_ESP_H2y3z_Pz_a-2.0E0*3*I_ESP_H2y3z_Pz_a+2*1*I_ESP_F3z_Pz;
    abcd[iGrid*378+250] = 4.0E0*I_ESP_K3y4z_Pz_aa-2.0E0*1*I_ESP_Hy4z_Pz_a-2.0E0*2*I_ESP_Hy4z_Pz_a;
    abcd[iGrid*378+251] = 4.0E0*I_ESP_K2y5z_Pz_aa-2.0E0*1*I_ESP_H5z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+252] = 4.0E0*I_ESP_K5xyz_Px_aa;
    abcd[iGrid*378+253] = 4.0E0*I_ESP_K4x2yz_Px_aa-2.0E0*1*I_ESP_H4xz_Px_a;
    abcd[iGrid*378+254] = 4.0E0*I_ESP_K4xy2z_Px_aa-2.0E0*1*I_ESP_H4xy_Px_a;
    abcd[iGrid*378+255] = 4.0E0*I_ESP_K3x3yz_Px_aa-2.0E0*2*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+256] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_H3x2y_Px_a-2.0E0*1*I_ESP_H3x2z_Px_a+1*I_ESP_F3x_Px;
    abcd[iGrid*378+257] = 4.0E0*I_ESP_K3xy3z_Px_aa-2.0E0*2*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+258] = 4.0E0*I_ESP_K2x4yz_Px_aa-2.0E0*3*I_ESP_H2x2yz_Px_a;
    abcd[iGrid*378+259] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*1*I_ESP_H2x3y_Px_a-2.0E0*2*I_ESP_H2xy2z_Px_a+2*1*I_ESP_F2xy_Px;
    abcd[iGrid*378+260] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*2*I_ESP_H2x2yz_Px_a-2.0E0*1*I_ESP_H2x3z_Px_a+2*I_ESP_F2xz_Px;
    abcd[iGrid*378+261] = 4.0E0*I_ESP_K2xy4z_Px_aa-2.0E0*3*I_ESP_H2xy2z_Px_a;
    abcd[iGrid*378+262] = 4.0E0*I_ESP_Kx5yz_Px_aa-2.0E0*4*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+263] = 4.0E0*I_ESP_Kx4y2z_Px_aa-2.0E0*1*I_ESP_Hx4y_Px_a-2.0E0*3*I_ESP_Hx2y2z_Px_a+3*1*I_ESP_Fx2y_Px;
    abcd[iGrid*378+264] = 4.0E0*I_ESP_Kx3y3z_Px_aa-2.0E0*2*I_ESP_Hx3yz_Px_a-2.0E0*2*I_ESP_Hxy3z_Px_a+2*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+265] = 4.0E0*I_ESP_Kx2y4z_Px_aa-2.0E0*3*I_ESP_Hx2y2z_Px_a-2.0E0*1*I_ESP_Hx4z_Px_a+3*I_ESP_Fx2z_Px;
    abcd[iGrid*378+266] = 4.0E0*I_ESP_Kxy5z_Px_aa-2.0E0*4*I_ESP_Hxy3z_Px_a;
    abcd[iGrid*378+267] = 4.0E0*I_ESP_K6yz_Px_aa-2.0E0*5*I_ESP_H4yz_Px_a;
    abcd[iGrid*378+268] = 4.0E0*I_ESP_K5y2z_Px_aa-2.0E0*1*I_ESP_H5y_Px_a-2.0E0*4*I_ESP_H3y2z_Px_a+4*1*I_ESP_F3y_Px;
    abcd[iGrid*378+269] = 4.0E0*I_ESP_K4y3z_Px_aa-2.0E0*2*I_ESP_H4yz_Px_a-2.0E0*3*I_ESP_H2y3z_Px_a+3*2*I_ESP_F2yz_Px;
    abcd[iGrid*378+270] = 4.0E0*I_ESP_K3y4z_Px_aa-2.0E0*3*I_ESP_H3y2z_Px_a-2.0E0*2*I_ESP_Hy4z_Px_a+2*3*I_ESP_Fy2z_Px;
    abcd[iGrid*378+271] = 4.0E0*I_ESP_K2y5z_Px_aa-2.0E0*4*I_ESP_H2y3z_Px_a-2.0E0*1*I_ESP_H5z_Px_a+4*I_ESP_F3z_Px;
    abcd[iGrid*378+272] = 4.0E0*I_ESP_Ky6z_Px_aa-2.0E0*5*I_ESP_Hy4z_Px_a;
    abcd[iGrid*378+273] = 4.0E0*I_ESP_K5xyz_Py_aa;
    abcd[iGrid*378+274] = 4.0E0*I_ESP_K4x2yz_Py_aa-2.0E0*1*I_ESP_H4xz_Py_a;
    abcd[iGrid*378+275] = 4.0E0*I_ESP_K4xy2z_Py_aa-2.0E0*1*I_ESP_H4xy_Py_a;
    abcd[iGrid*378+276] = 4.0E0*I_ESP_K3x3yz_Py_aa-2.0E0*2*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+277] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_H3x2y_Py_a-2.0E0*1*I_ESP_H3x2z_Py_a+1*I_ESP_F3x_Py;
    abcd[iGrid*378+278] = 4.0E0*I_ESP_K3xy3z_Py_aa-2.0E0*2*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+279] = 4.0E0*I_ESP_K2x4yz_Py_aa-2.0E0*3*I_ESP_H2x2yz_Py_a;
    abcd[iGrid*378+280] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*1*I_ESP_H2x3y_Py_a-2.0E0*2*I_ESP_H2xy2z_Py_a+2*1*I_ESP_F2xy_Py;
    abcd[iGrid*378+281] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*2*I_ESP_H2x2yz_Py_a-2.0E0*1*I_ESP_H2x3z_Py_a+2*I_ESP_F2xz_Py;
    abcd[iGrid*378+282] = 4.0E0*I_ESP_K2xy4z_Py_aa-2.0E0*3*I_ESP_H2xy2z_Py_a;
    abcd[iGrid*378+283] = 4.0E0*I_ESP_Kx5yz_Py_aa-2.0E0*4*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+284] = 4.0E0*I_ESP_Kx4y2z_Py_aa-2.0E0*1*I_ESP_Hx4y_Py_a-2.0E0*3*I_ESP_Hx2y2z_Py_a+3*1*I_ESP_Fx2y_Py;
    abcd[iGrid*378+285] = 4.0E0*I_ESP_Kx3y3z_Py_aa-2.0E0*2*I_ESP_Hx3yz_Py_a-2.0E0*2*I_ESP_Hxy3z_Py_a+2*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+286] = 4.0E0*I_ESP_Kx2y4z_Py_aa-2.0E0*3*I_ESP_Hx2y2z_Py_a-2.0E0*1*I_ESP_Hx4z_Py_a+3*I_ESP_Fx2z_Py;
    abcd[iGrid*378+287] = 4.0E0*I_ESP_Kxy5z_Py_aa-2.0E0*4*I_ESP_Hxy3z_Py_a;
    abcd[iGrid*378+288] = 4.0E0*I_ESP_K6yz_Py_aa-2.0E0*5*I_ESP_H4yz_Py_a;
    abcd[iGrid*378+289] = 4.0E0*I_ESP_K5y2z_Py_aa-2.0E0*1*I_ESP_H5y_Py_a-2.0E0*4*I_ESP_H3y2z_Py_a+4*1*I_ESP_F3y_Py;
    abcd[iGrid*378+290] = 4.0E0*I_ESP_K4y3z_Py_aa-2.0E0*2*I_ESP_H4yz_Py_a-2.0E0*3*I_ESP_H2y3z_Py_a+3*2*I_ESP_F2yz_Py;
    abcd[iGrid*378+291] = 4.0E0*I_ESP_K3y4z_Py_aa-2.0E0*3*I_ESP_H3y2z_Py_a-2.0E0*2*I_ESP_Hy4z_Py_a+2*3*I_ESP_Fy2z_Py;
    abcd[iGrid*378+292] = 4.0E0*I_ESP_K2y5z_Py_aa-2.0E0*4*I_ESP_H2y3z_Py_a-2.0E0*1*I_ESP_H5z_Py_a+4*I_ESP_F3z_Py;
    abcd[iGrid*378+293] = 4.0E0*I_ESP_Ky6z_Py_aa-2.0E0*5*I_ESP_Hy4z_Py_a;
    abcd[iGrid*378+294] = 4.0E0*I_ESP_K5xyz_Pz_aa;
    abcd[iGrid*378+295] = 4.0E0*I_ESP_K4x2yz_Pz_aa-2.0E0*1*I_ESP_H4xz_Pz_a;
    abcd[iGrid*378+296] = 4.0E0*I_ESP_K4xy2z_Pz_aa-2.0E0*1*I_ESP_H4xy_Pz_a;
    abcd[iGrid*378+297] = 4.0E0*I_ESP_K3x3yz_Pz_aa-2.0E0*2*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+298] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_H3x2y_Pz_a-2.0E0*1*I_ESP_H3x2z_Pz_a+1*I_ESP_F3x_Pz;
    abcd[iGrid*378+299] = 4.0E0*I_ESP_K3xy3z_Pz_aa-2.0E0*2*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+300] = 4.0E0*I_ESP_K2x4yz_Pz_aa-2.0E0*3*I_ESP_H2x2yz_Pz_a;
    abcd[iGrid*378+301] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*1*I_ESP_H2x3y_Pz_a-2.0E0*2*I_ESP_H2xy2z_Pz_a+2*1*I_ESP_F2xy_Pz;
    abcd[iGrid*378+302] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*2*I_ESP_H2x2yz_Pz_a-2.0E0*1*I_ESP_H2x3z_Pz_a+2*I_ESP_F2xz_Pz;
    abcd[iGrid*378+303] = 4.0E0*I_ESP_K2xy4z_Pz_aa-2.0E0*3*I_ESP_H2xy2z_Pz_a;
    abcd[iGrid*378+304] = 4.0E0*I_ESP_Kx5yz_Pz_aa-2.0E0*4*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+305] = 4.0E0*I_ESP_Kx4y2z_Pz_aa-2.0E0*1*I_ESP_Hx4y_Pz_a-2.0E0*3*I_ESP_Hx2y2z_Pz_a+3*1*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+306] = 4.0E0*I_ESP_Kx3y3z_Pz_aa-2.0E0*2*I_ESP_Hx3yz_Pz_a-2.0E0*2*I_ESP_Hxy3z_Pz_a+2*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+307] = 4.0E0*I_ESP_Kx2y4z_Pz_aa-2.0E0*3*I_ESP_Hx2y2z_Pz_a-2.0E0*1*I_ESP_Hx4z_Pz_a+3*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+308] = 4.0E0*I_ESP_Kxy5z_Pz_aa-2.0E0*4*I_ESP_Hxy3z_Pz_a;
    abcd[iGrid*378+309] = 4.0E0*I_ESP_K6yz_Pz_aa-2.0E0*5*I_ESP_H4yz_Pz_a;
    abcd[iGrid*378+310] = 4.0E0*I_ESP_K5y2z_Pz_aa-2.0E0*1*I_ESP_H5y_Pz_a-2.0E0*4*I_ESP_H3y2z_Pz_a+4*1*I_ESP_F3y_Pz;
    abcd[iGrid*378+311] = 4.0E0*I_ESP_K4y3z_Pz_aa-2.0E0*2*I_ESP_H4yz_Pz_a-2.0E0*3*I_ESP_H2y3z_Pz_a+3*2*I_ESP_F2yz_Pz;
    abcd[iGrid*378+312] = 4.0E0*I_ESP_K3y4z_Pz_aa-2.0E0*3*I_ESP_H3y2z_Pz_a-2.0E0*2*I_ESP_Hy4z_Pz_a+2*3*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+313] = 4.0E0*I_ESP_K2y5z_Pz_aa-2.0E0*4*I_ESP_H2y3z_Pz_a-2.0E0*1*I_ESP_H5z_Pz_a+4*I_ESP_F3z_Pz;
    abcd[iGrid*378+314] = 4.0E0*I_ESP_Ky6z_Pz_aa-2.0E0*5*I_ESP_Hy4z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    abcd[iGrid*378+315] = 4.0E0*I_ESP_K5x2z_Px_aa-2.0E0*1*I_ESP_H5x_Px_a;
    abcd[iGrid*378+316] = 4.0E0*I_ESP_K4xy2z_Px_aa-2.0E0*1*I_ESP_H4xy_Px_a;
    abcd[iGrid*378+317] = 4.0E0*I_ESP_K4x3z_Px_aa-2.0E0*1*I_ESP_H4xz_Px_a-2.0E0*2*I_ESP_H4xz_Px_a;
    abcd[iGrid*378+318] = 4.0E0*I_ESP_K3x2y2z_Px_aa-2.0E0*1*I_ESP_H3x2y_Px_a;
    abcd[iGrid*378+319] = 4.0E0*I_ESP_K3xy3z_Px_aa-2.0E0*1*I_ESP_H3xyz_Px_a-2.0E0*2*I_ESP_H3xyz_Px_a;
    abcd[iGrid*378+320] = 4.0E0*I_ESP_K3x4z_Px_aa-2.0E0*2*I_ESP_H3x2z_Px_a-2.0E0*3*I_ESP_H3x2z_Px_a+2*1*I_ESP_F3x_Px;
    abcd[iGrid*378+321] = 4.0E0*I_ESP_K2x3y2z_Px_aa-2.0E0*1*I_ESP_H2x3y_Px_a;
    abcd[iGrid*378+322] = 4.0E0*I_ESP_K2x2y3z_Px_aa-2.0E0*1*I_ESP_H2x2yz_Px_a-2.0E0*2*I_ESP_H2x2yz_Px_a;
    abcd[iGrid*378+323] = 4.0E0*I_ESP_K2xy4z_Px_aa-2.0E0*2*I_ESP_H2xy2z_Px_a-2.0E0*3*I_ESP_H2xy2z_Px_a+2*1*I_ESP_F2xy_Px;
    abcd[iGrid*378+324] = 4.0E0*I_ESP_K2x5z_Px_aa-2.0E0*3*I_ESP_H2x3z_Px_a-2.0E0*4*I_ESP_H2x3z_Px_a+3*2*I_ESP_F2xz_Px;
    abcd[iGrid*378+325] = 4.0E0*I_ESP_Kx4y2z_Px_aa-2.0E0*1*I_ESP_Hx4y_Px_a;
    abcd[iGrid*378+326] = 4.0E0*I_ESP_Kx3y3z_Px_aa-2.0E0*1*I_ESP_Hx3yz_Px_a-2.0E0*2*I_ESP_Hx3yz_Px_a;
    abcd[iGrid*378+327] = 4.0E0*I_ESP_Kx2y4z_Px_aa-2.0E0*2*I_ESP_Hx2y2z_Px_a-2.0E0*3*I_ESP_Hx2y2z_Px_a+2*1*I_ESP_Fx2y_Px;
    abcd[iGrid*378+328] = 4.0E0*I_ESP_Kxy5z_Px_aa-2.0E0*3*I_ESP_Hxy3z_Px_a-2.0E0*4*I_ESP_Hxy3z_Px_a+3*2*I_ESP_Fxyz_Px;
    abcd[iGrid*378+329] = 4.0E0*I_ESP_Kx6z_Px_aa-2.0E0*4*I_ESP_Hx4z_Px_a-2.0E0*5*I_ESP_Hx4z_Px_a+4*3*I_ESP_Fx2z_Px;
    abcd[iGrid*378+330] = 4.0E0*I_ESP_K5y2z_Px_aa-2.0E0*1*I_ESP_H5y_Px_a;
    abcd[iGrid*378+331] = 4.0E0*I_ESP_K4y3z_Px_aa-2.0E0*1*I_ESP_H4yz_Px_a-2.0E0*2*I_ESP_H4yz_Px_a;
    abcd[iGrid*378+332] = 4.0E0*I_ESP_K3y4z_Px_aa-2.0E0*2*I_ESP_H3y2z_Px_a-2.0E0*3*I_ESP_H3y2z_Px_a+2*1*I_ESP_F3y_Px;
    abcd[iGrid*378+333] = 4.0E0*I_ESP_K2y5z_Px_aa-2.0E0*3*I_ESP_H2y3z_Px_a-2.0E0*4*I_ESP_H2y3z_Px_a+3*2*I_ESP_F2yz_Px;
    abcd[iGrid*378+334] = 4.0E0*I_ESP_Ky6z_Px_aa-2.0E0*4*I_ESP_Hy4z_Px_a-2.0E0*5*I_ESP_Hy4z_Px_a+4*3*I_ESP_Fy2z_Px;
    abcd[iGrid*378+335] = 4.0E0*I_ESP_K7z_Px_aa-2.0E0*5*I_ESP_H5z_Px_a-2.0E0*6*I_ESP_H5z_Px_a+5*4*I_ESP_F3z_Px;
    abcd[iGrid*378+336] = 4.0E0*I_ESP_K5x2z_Py_aa-2.0E0*1*I_ESP_H5x_Py_a;
    abcd[iGrid*378+337] = 4.0E0*I_ESP_K4xy2z_Py_aa-2.0E0*1*I_ESP_H4xy_Py_a;
    abcd[iGrid*378+338] = 4.0E0*I_ESP_K4x3z_Py_aa-2.0E0*1*I_ESP_H4xz_Py_a-2.0E0*2*I_ESP_H4xz_Py_a;
    abcd[iGrid*378+339] = 4.0E0*I_ESP_K3x2y2z_Py_aa-2.0E0*1*I_ESP_H3x2y_Py_a;
    abcd[iGrid*378+340] = 4.0E0*I_ESP_K3xy3z_Py_aa-2.0E0*1*I_ESP_H3xyz_Py_a-2.0E0*2*I_ESP_H3xyz_Py_a;
    abcd[iGrid*378+341] = 4.0E0*I_ESP_K3x4z_Py_aa-2.0E0*2*I_ESP_H3x2z_Py_a-2.0E0*3*I_ESP_H3x2z_Py_a+2*1*I_ESP_F3x_Py;
    abcd[iGrid*378+342] = 4.0E0*I_ESP_K2x3y2z_Py_aa-2.0E0*1*I_ESP_H2x3y_Py_a;
    abcd[iGrid*378+343] = 4.0E0*I_ESP_K2x2y3z_Py_aa-2.0E0*1*I_ESP_H2x2yz_Py_a-2.0E0*2*I_ESP_H2x2yz_Py_a;
    abcd[iGrid*378+344] = 4.0E0*I_ESP_K2xy4z_Py_aa-2.0E0*2*I_ESP_H2xy2z_Py_a-2.0E0*3*I_ESP_H2xy2z_Py_a+2*1*I_ESP_F2xy_Py;
    abcd[iGrid*378+345] = 4.0E0*I_ESP_K2x5z_Py_aa-2.0E0*3*I_ESP_H2x3z_Py_a-2.0E0*4*I_ESP_H2x3z_Py_a+3*2*I_ESP_F2xz_Py;
    abcd[iGrid*378+346] = 4.0E0*I_ESP_Kx4y2z_Py_aa-2.0E0*1*I_ESP_Hx4y_Py_a;
    abcd[iGrid*378+347] = 4.0E0*I_ESP_Kx3y3z_Py_aa-2.0E0*1*I_ESP_Hx3yz_Py_a-2.0E0*2*I_ESP_Hx3yz_Py_a;
    abcd[iGrid*378+348] = 4.0E0*I_ESP_Kx2y4z_Py_aa-2.0E0*2*I_ESP_Hx2y2z_Py_a-2.0E0*3*I_ESP_Hx2y2z_Py_a+2*1*I_ESP_Fx2y_Py;
    abcd[iGrid*378+349] = 4.0E0*I_ESP_Kxy5z_Py_aa-2.0E0*3*I_ESP_Hxy3z_Py_a-2.0E0*4*I_ESP_Hxy3z_Py_a+3*2*I_ESP_Fxyz_Py;
    abcd[iGrid*378+350] = 4.0E0*I_ESP_Kx6z_Py_aa-2.0E0*4*I_ESP_Hx4z_Py_a-2.0E0*5*I_ESP_Hx4z_Py_a+4*3*I_ESP_Fx2z_Py;
    abcd[iGrid*378+351] = 4.0E0*I_ESP_K5y2z_Py_aa-2.0E0*1*I_ESP_H5y_Py_a;
    abcd[iGrid*378+352] = 4.0E0*I_ESP_K4y3z_Py_aa-2.0E0*1*I_ESP_H4yz_Py_a-2.0E0*2*I_ESP_H4yz_Py_a;
    abcd[iGrid*378+353] = 4.0E0*I_ESP_K3y4z_Py_aa-2.0E0*2*I_ESP_H3y2z_Py_a-2.0E0*3*I_ESP_H3y2z_Py_a+2*1*I_ESP_F3y_Py;
    abcd[iGrid*378+354] = 4.0E0*I_ESP_K2y5z_Py_aa-2.0E0*3*I_ESP_H2y3z_Py_a-2.0E0*4*I_ESP_H2y3z_Py_a+3*2*I_ESP_F2yz_Py;
    abcd[iGrid*378+355] = 4.0E0*I_ESP_Ky6z_Py_aa-2.0E0*4*I_ESP_Hy4z_Py_a-2.0E0*5*I_ESP_Hy4z_Py_a+4*3*I_ESP_Fy2z_Py;
    abcd[iGrid*378+356] = 4.0E0*I_ESP_K7z_Py_aa-2.0E0*5*I_ESP_H5z_Py_a-2.0E0*6*I_ESP_H5z_Py_a+5*4*I_ESP_F3z_Py;
    abcd[iGrid*378+357] = 4.0E0*I_ESP_K5x2z_Pz_aa-2.0E0*1*I_ESP_H5x_Pz_a;
    abcd[iGrid*378+358] = 4.0E0*I_ESP_K4xy2z_Pz_aa-2.0E0*1*I_ESP_H4xy_Pz_a;
    abcd[iGrid*378+359] = 4.0E0*I_ESP_K4x3z_Pz_aa-2.0E0*1*I_ESP_H4xz_Pz_a-2.0E0*2*I_ESP_H4xz_Pz_a;
    abcd[iGrid*378+360] = 4.0E0*I_ESP_K3x2y2z_Pz_aa-2.0E0*1*I_ESP_H3x2y_Pz_a;
    abcd[iGrid*378+361] = 4.0E0*I_ESP_K3xy3z_Pz_aa-2.0E0*1*I_ESP_H3xyz_Pz_a-2.0E0*2*I_ESP_H3xyz_Pz_a;
    abcd[iGrid*378+362] = 4.0E0*I_ESP_K3x4z_Pz_aa-2.0E0*2*I_ESP_H3x2z_Pz_a-2.0E0*3*I_ESP_H3x2z_Pz_a+2*1*I_ESP_F3x_Pz;
    abcd[iGrid*378+363] = 4.0E0*I_ESP_K2x3y2z_Pz_aa-2.0E0*1*I_ESP_H2x3y_Pz_a;
    abcd[iGrid*378+364] = 4.0E0*I_ESP_K2x2y3z_Pz_aa-2.0E0*1*I_ESP_H2x2yz_Pz_a-2.0E0*2*I_ESP_H2x2yz_Pz_a;
    abcd[iGrid*378+365] = 4.0E0*I_ESP_K2xy4z_Pz_aa-2.0E0*2*I_ESP_H2xy2z_Pz_a-2.0E0*3*I_ESP_H2xy2z_Pz_a+2*1*I_ESP_F2xy_Pz;
    abcd[iGrid*378+366] = 4.0E0*I_ESP_K2x5z_Pz_aa-2.0E0*3*I_ESP_H2x3z_Pz_a-2.0E0*4*I_ESP_H2x3z_Pz_a+3*2*I_ESP_F2xz_Pz;
    abcd[iGrid*378+367] = 4.0E0*I_ESP_Kx4y2z_Pz_aa-2.0E0*1*I_ESP_Hx4y_Pz_a;
    abcd[iGrid*378+368] = 4.0E0*I_ESP_Kx3y3z_Pz_aa-2.0E0*1*I_ESP_Hx3yz_Pz_a-2.0E0*2*I_ESP_Hx3yz_Pz_a;
    abcd[iGrid*378+369] = 4.0E0*I_ESP_Kx2y4z_Pz_aa-2.0E0*2*I_ESP_Hx2y2z_Pz_a-2.0E0*3*I_ESP_Hx2y2z_Pz_a+2*1*I_ESP_Fx2y_Pz;
    abcd[iGrid*378+370] = 4.0E0*I_ESP_Kxy5z_Pz_aa-2.0E0*3*I_ESP_Hxy3z_Pz_a-2.0E0*4*I_ESP_Hxy3z_Pz_a+3*2*I_ESP_Fxyz_Pz;
    abcd[iGrid*378+371] = 4.0E0*I_ESP_Kx6z_Pz_aa-2.0E0*4*I_ESP_Hx4z_Pz_a-2.0E0*5*I_ESP_Hx4z_Pz_a+4*3*I_ESP_Fx2z_Pz;
    abcd[iGrid*378+372] = 4.0E0*I_ESP_K5y2z_Pz_aa-2.0E0*1*I_ESP_H5y_Pz_a;
    abcd[iGrid*378+373] = 4.0E0*I_ESP_K4y3z_Pz_aa-2.0E0*1*I_ESP_H4yz_Pz_a-2.0E0*2*I_ESP_H4yz_Pz_a;
    abcd[iGrid*378+374] = 4.0E0*I_ESP_K3y4z_Pz_aa-2.0E0*2*I_ESP_H3y2z_Pz_a-2.0E0*3*I_ESP_H3y2z_Pz_a+2*1*I_ESP_F3y_Pz;
    abcd[iGrid*378+375] = 4.0E0*I_ESP_K2y5z_Pz_aa-2.0E0*3*I_ESP_H2y3z_Pz_a-2.0E0*4*I_ESP_H2y3z_Pz_a+3*2*I_ESP_F2yz_Pz;
    abcd[iGrid*378+376] = 4.0E0*I_ESP_Ky6z_Pz_aa-2.0E0*4*I_ESP_Hy4z_Pz_a-2.0E0*5*I_ESP_Hy4z_Pz_a+4*3*I_ESP_Fy2z_Pz;
    abcd[iGrid*378+377] = 4.0E0*I_ESP_K7z_Pz_aa-2.0E0*5*I_ESP_H5z_Pz_a-2.0E0*6*I_ESP_H5z_Pz_a+5*4*I_ESP_F3z_Pz;
  }
}
