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
// BRA1 as redundant position, total RHS integrals evaluated as: 4295
// BRA2 as redundant position, total RHS integrals evaluated as: 3594
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

void hgp_os_esp_g_sp_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_I6x_S_C4_aa = 0.0E0;
    Double I_ESP_I5xy_S_C4_aa = 0.0E0;
    Double I_ESP_I5xz_S_C4_aa = 0.0E0;
    Double I_ESP_I4x2y_S_C4_aa = 0.0E0;
    Double I_ESP_I4xyz_S_C4_aa = 0.0E0;
    Double I_ESP_I4x2z_S_C4_aa = 0.0E0;
    Double I_ESP_I3x3y_S_C4_aa = 0.0E0;
    Double I_ESP_I3x2yz_S_C4_aa = 0.0E0;
    Double I_ESP_I3xy2z_S_C4_aa = 0.0E0;
    Double I_ESP_I3x3z_S_C4_aa = 0.0E0;
    Double I_ESP_I2x4y_S_C4_aa = 0.0E0;
    Double I_ESP_I2x3yz_S_C4_aa = 0.0E0;
    Double I_ESP_I2x2y2z_S_C4_aa = 0.0E0;
    Double I_ESP_I2xy3z_S_C4_aa = 0.0E0;
    Double I_ESP_I2x4z_S_C4_aa = 0.0E0;
    Double I_ESP_Ix5y_S_C4_aa = 0.0E0;
    Double I_ESP_Ix4yz_S_C4_aa = 0.0E0;
    Double I_ESP_Ix3y2z_S_C4_aa = 0.0E0;
    Double I_ESP_Ix2y3z_S_C4_aa = 0.0E0;
    Double I_ESP_Ixy4z_S_C4_aa = 0.0E0;
    Double I_ESP_Ix5z_S_C4_aa = 0.0E0;
    Double I_ESP_I6y_S_C4_aa = 0.0E0;
    Double I_ESP_I5yz_S_C4_aa = 0.0E0;
    Double I_ESP_I4y2z_S_C4_aa = 0.0E0;
    Double I_ESP_I3y3z_S_C4_aa = 0.0E0;
    Double I_ESP_I2y4z_S_C4_aa = 0.0E0;
    Double I_ESP_Iy5z_S_C4_aa = 0.0E0;
    Double I_ESP_I6z_S_C4_aa = 0.0E0;
    Double I_ESP_G4x_S_C4_a = 0.0E0;
    Double I_ESP_G3xy_S_C4_a = 0.0E0;
    Double I_ESP_G3xz_S_C4_a = 0.0E0;
    Double I_ESP_G2x2y_S_C4_a = 0.0E0;
    Double I_ESP_G2xyz_S_C4_a = 0.0E0;
    Double I_ESP_G2x2z_S_C4_a = 0.0E0;
    Double I_ESP_Gx3y_S_C4_a = 0.0E0;
    Double I_ESP_Gx2yz_S_C4_a = 0.0E0;
    Double I_ESP_Gxy2z_S_C4_a = 0.0E0;
    Double I_ESP_Gx3z_S_C4_a = 0.0E0;
    Double I_ESP_G4y_S_C4_a = 0.0E0;
    Double I_ESP_G3yz_S_C4_a = 0.0E0;
    Double I_ESP_G2y2z_S_C4_a = 0.0E0;
    Double I_ESP_Gy3z_S_C4_a = 0.0E0;
    Double I_ESP_G4z_S_C4_a = 0.0E0;
    Double I_ESP_D2x_S_C4 = 0.0E0;
    Double I_ESP_Dxy_S_C4 = 0.0E0;
    Double I_ESP_Dxz_S_C4 = 0.0E0;
    Double I_ESP_D2y_S_C4 = 0.0E0;
    Double I_ESP_Dyz_S_C4 = 0.0E0;
    Double I_ESP_D2z_S_C4 = 0.0E0;
    Double I_ESP_K7x_S_C1004_aa = 0.0E0;
    Double I_ESP_K6xy_S_C1004_aa = 0.0E0;
    Double I_ESP_K6xz_S_C1004_aa = 0.0E0;
    Double I_ESP_K5x2y_S_C1004_aa = 0.0E0;
    Double I_ESP_K5xyz_S_C1004_aa = 0.0E0;
    Double I_ESP_K5x2z_S_C1004_aa = 0.0E0;
    Double I_ESP_K4x3y_S_C1004_aa = 0.0E0;
    Double I_ESP_K4x2yz_S_C1004_aa = 0.0E0;
    Double I_ESP_K4xy2z_S_C1004_aa = 0.0E0;
    Double I_ESP_K4x3z_S_C1004_aa = 0.0E0;
    Double I_ESP_K3x4y_S_C1004_aa = 0.0E0;
    Double I_ESP_K3x3yz_S_C1004_aa = 0.0E0;
    Double I_ESP_K3x2y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_K3xy3z_S_C1004_aa = 0.0E0;
    Double I_ESP_K3x4z_S_C1004_aa = 0.0E0;
    Double I_ESP_K2x5y_S_C1004_aa = 0.0E0;
    Double I_ESP_K2x4yz_S_C1004_aa = 0.0E0;
    Double I_ESP_K2x3y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_K2x2y3z_S_C1004_aa = 0.0E0;
    Double I_ESP_K2xy4z_S_C1004_aa = 0.0E0;
    Double I_ESP_K2x5z_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx6y_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx5yz_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx4y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx3y3z_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx2y4z_S_C1004_aa = 0.0E0;
    Double I_ESP_Kxy5z_S_C1004_aa = 0.0E0;
    Double I_ESP_Kx6z_S_C1004_aa = 0.0E0;
    Double I_ESP_K7y_S_C1004_aa = 0.0E0;
    Double I_ESP_K6yz_S_C1004_aa = 0.0E0;
    Double I_ESP_K5y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_K4y3z_S_C1004_aa = 0.0E0;
    Double I_ESP_K3y4z_S_C1004_aa = 0.0E0;
    Double I_ESP_K2y5z_S_C1004_aa = 0.0E0;
    Double I_ESP_Ky6z_S_C1004_aa = 0.0E0;
    Double I_ESP_K7z_S_C1004_aa = 0.0E0;
    Double I_ESP_I6x_S_C1004_aa = 0.0E0;
    Double I_ESP_I5xy_S_C1004_aa = 0.0E0;
    Double I_ESP_I5xz_S_C1004_aa = 0.0E0;
    Double I_ESP_I4x2y_S_C1004_aa = 0.0E0;
    Double I_ESP_I4xyz_S_C1004_aa = 0.0E0;
    Double I_ESP_I4x2z_S_C1004_aa = 0.0E0;
    Double I_ESP_I3x3y_S_C1004_aa = 0.0E0;
    Double I_ESP_I3x2yz_S_C1004_aa = 0.0E0;
    Double I_ESP_I3xy2z_S_C1004_aa = 0.0E0;
    Double I_ESP_I3x3z_S_C1004_aa = 0.0E0;
    Double I_ESP_I2x4y_S_C1004_aa = 0.0E0;
    Double I_ESP_I2x3yz_S_C1004_aa = 0.0E0;
    Double I_ESP_I2x2y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_I2xy3z_S_C1004_aa = 0.0E0;
    Double I_ESP_I2x4z_S_C1004_aa = 0.0E0;
    Double I_ESP_Ix5y_S_C1004_aa = 0.0E0;
    Double I_ESP_Ix4yz_S_C1004_aa = 0.0E0;
    Double I_ESP_Ix3y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_Ix2y3z_S_C1004_aa = 0.0E0;
    Double I_ESP_Ixy4z_S_C1004_aa = 0.0E0;
    Double I_ESP_Ix5z_S_C1004_aa = 0.0E0;
    Double I_ESP_I6y_S_C1004_aa = 0.0E0;
    Double I_ESP_I5yz_S_C1004_aa = 0.0E0;
    Double I_ESP_I4y2z_S_C1004_aa = 0.0E0;
    Double I_ESP_I3y3z_S_C1004_aa = 0.0E0;
    Double I_ESP_I2y4z_S_C1004_aa = 0.0E0;
    Double I_ESP_Iy5z_S_C1004_aa = 0.0E0;
    Double I_ESP_I6z_S_C1004_aa = 0.0E0;
    Double I_ESP_H5x_S_C1004_a = 0.0E0;
    Double I_ESP_H4xy_S_C1004_a = 0.0E0;
    Double I_ESP_H4xz_S_C1004_a = 0.0E0;
    Double I_ESP_H3x2y_S_C1004_a = 0.0E0;
    Double I_ESP_H3xyz_S_C1004_a = 0.0E0;
    Double I_ESP_H3x2z_S_C1004_a = 0.0E0;
    Double I_ESP_H2x3y_S_C1004_a = 0.0E0;
    Double I_ESP_H2x2yz_S_C1004_a = 0.0E0;
    Double I_ESP_H2xy2z_S_C1004_a = 0.0E0;
    Double I_ESP_H2x3z_S_C1004_a = 0.0E0;
    Double I_ESP_Hx4y_S_C1004_a = 0.0E0;
    Double I_ESP_Hx3yz_S_C1004_a = 0.0E0;
    Double I_ESP_Hx2y2z_S_C1004_a = 0.0E0;
    Double I_ESP_Hxy3z_S_C1004_a = 0.0E0;
    Double I_ESP_Hx4z_S_C1004_a = 0.0E0;
    Double I_ESP_H5y_S_C1004_a = 0.0E0;
    Double I_ESP_H4yz_S_C1004_a = 0.0E0;
    Double I_ESP_H3y2z_S_C1004_a = 0.0E0;
    Double I_ESP_H2y3z_S_C1004_a = 0.0E0;
    Double I_ESP_Hy4z_S_C1004_a = 0.0E0;
    Double I_ESP_H5z_S_C1004_a = 0.0E0;
    Double I_ESP_G4x_S_C1004_a = 0.0E0;
    Double I_ESP_G3xy_S_C1004_a = 0.0E0;
    Double I_ESP_G3xz_S_C1004_a = 0.0E0;
    Double I_ESP_G2x2y_S_C1004_a = 0.0E0;
    Double I_ESP_G2xyz_S_C1004_a = 0.0E0;
    Double I_ESP_G2x2z_S_C1004_a = 0.0E0;
    Double I_ESP_Gx3y_S_C1004_a = 0.0E0;
    Double I_ESP_Gx2yz_S_C1004_a = 0.0E0;
    Double I_ESP_Gxy2z_S_C1004_a = 0.0E0;
    Double I_ESP_Gx3z_S_C1004_a = 0.0E0;
    Double I_ESP_G4y_S_C1004_a = 0.0E0;
    Double I_ESP_G3yz_S_C1004_a = 0.0E0;
    Double I_ESP_G2y2z_S_C1004_a = 0.0E0;
    Double I_ESP_Gy3z_S_C1004_a = 0.0E0;
    Double I_ESP_G4z_S_C1004_a = 0.0E0;
    Double I_ESP_F3x_S_C1004 = 0.0E0;
    Double I_ESP_F2xy_S_C1004 = 0.0E0;
    Double I_ESP_F2xz_S_C1004 = 0.0E0;
    Double I_ESP_Fx2y_S_C1004 = 0.0E0;
    Double I_ESP_Fxyz_S_C1004 = 0.0E0;
    Double I_ESP_Fx2z_S_C1004 = 0.0E0;
    Double I_ESP_F3y_S_C1004 = 0.0E0;
    Double I_ESP_F2yz_S_C1004 = 0.0E0;
    Double I_ESP_Fy2z_S_C1004 = 0.0E0;
    Double I_ESP_F3z_S_C1004 = 0.0E0;
    Double I_ESP_D2x_S_C1004 = 0.0E0;
    Double I_ESP_Dxy_S_C1004 = 0.0E0;
    Double I_ESP_Dxz_S_C1004 = 0.0E0;
    Double I_ESP_D2y_S_C1004 = 0.0E0;
    Double I_ESP_Dyz_S_C1004 = 0.0E0;
    Double I_ESP_D2z_S_C1004 = 0.0E0;

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
      Double ic2_1 = icoe[ip2+1*inp2];
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
      Double prefactor = fbra;

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
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER49;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER19*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = 1.0E0+u2*ONEOVER17*I_ESP_S_S_M7_vrr;
        I_ESP_S_S_M7_vrr = ONEOVER15*I_ESP_S_S_M7_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M7_vrr  = f*I_ESP_S_S_M7_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);
        I_ESP_S_S_M6_vrr = static_cast<Double>(I_ESP_S_S_M6_vrr_d);
        I_ESP_S_S_M7_vrr = static_cast<Double>(I_ESP_S_S_M7_vrr_d);

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

#endif

      }


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
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_Dxz_S_vrr = PAZ*I_ESP_Px_S_vrr-PRZ*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dyz_S_vrr = PAZ*I_ESP_Py_S_vrr-PRZ*I_ESP_Py_S_M1_vrr;
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
         * shell quartet name: SQ_ESP_I_S_C4_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_C4_aa_coefs = ic2*alpha*alpha;
        I_ESP_I6x_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_C4_aa += SQ_ESP_I_S_C4_aa_coefs*I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_C4_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_G_S_C4_a_coefs = ic2*alpha;
        I_ESP_G4x_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S_C4_a += SQ_ESP_G_S_C4_a_coefs*I_ESP_G4z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_C4
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_D_S_C4_coefs = ic2;
        I_ESP_D2x_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_D2x_S_vrr;
        I_ESP_Dxy_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_Dxy_S_vrr;
        I_ESP_Dxz_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_Dxz_S_vrr;
        I_ESP_D2y_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_D2y_S_vrr;
        I_ESP_Dyz_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_Dyz_S_vrr;
        I_ESP_D2z_S_C4 += SQ_ESP_D_S_C4_coefs*I_ESP_D2z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_C1004_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_K_S_C1004_aa_coefs = ic2_1*alpha*alpha;
        I_ESP_K7x_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S_C1004_aa += SQ_ESP_K_S_C1004_aa_coefs*I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_C1004_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_C1004_aa_coefs = ic2_1*alpha*alpha;
        I_ESP_I6x_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_C1004_aa += SQ_ESP_I_S_C1004_aa_coefs*I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_C1004_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_C1004_a_coefs = ic2_1*alpha;
        I_ESP_H5x_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_C1004_a += SQ_ESP_H_S_C1004_a_coefs*I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_C1004_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_G_S_C1004_a_coefs = ic2_1*alpha;
        I_ESP_G4x_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S_C1004_a += SQ_ESP_G_S_C1004_a_coefs*I_ESP_G4z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_C1004
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_F_S_C1004_coefs = ic2_1;
        I_ESP_F3x_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S_C1004 += SQ_ESP_F_S_C1004_coefs*I_ESP_F3z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_C1004
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_D_S_C1004_coefs = ic2_1;
        I_ESP_D2x_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_D2x_S_vrr;
        I_ESP_Dxy_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_Dxy_S_vrr;
        I_ESP_Dxz_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_Dxz_S_vrr;
        I_ESP_D2y_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_D2y_S_vrr;
        I_ESP_Dyz_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_Dyz_S_vrr;
        I_ESP_D2z_S_C1004 += SQ_ESP_D_S_C1004_coefs*I_ESP_D2z_S_vrr;
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
     * shell quartet name: SQ_ESP_D_P_C1004
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_S_C1004
     * RHS shell quartet name: SQ_ESP_D_S_C1004
     ************************************************************/
    Double I_ESP_D2x_Px_C1004 = I_ESP_F3x_S_C1004+ABX*I_ESP_D2x_S_C1004;
    Double I_ESP_Dxy_Px_C1004 = I_ESP_F2xy_S_C1004+ABX*I_ESP_Dxy_S_C1004;
    Double I_ESP_Dxz_Px_C1004 = I_ESP_F2xz_S_C1004+ABX*I_ESP_Dxz_S_C1004;
    Double I_ESP_D2y_Px_C1004 = I_ESP_Fx2y_S_C1004+ABX*I_ESP_D2y_S_C1004;
    Double I_ESP_Dyz_Px_C1004 = I_ESP_Fxyz_S_C1004+ABX*I_ESP_Dyz_S_C1004;
    Double I_ESP_D2z_Px_C1004 = I_ESP_Fx2z_S_C1004+ABX*I_ESP_D2z_S_C1004;
    Double I_ESP_D2x_Py_C1004 = I_ESP_F2xy_S_C1004+ABY*I_ESP_D2x_S_C1004;
    Double I_ESP_Dxy_Py_C1004 = I_ESP_Fx2y_S_C1004+ABY*I_ESP_Dxy_S_C1004;
    Double I_ESP_Dxz_Py_C1004 = I_ESP_Fxyz_S_C1004+ABY*I_ESP_Dxz_S_C1004;
    Double I_ESP_D2y_Py_C1004 = I_ESP_F3y_S_C1004+ABY*I_ESP_D2y_S_C1004;
    Double I_ESP_Dyz_Py_C1004 = I_ESP_F2yz_S_C1004+ABY*I_ESP_Dyz_S_C1004;
    Double I_ESP_D2z_Py_C1004 = I_ESP_Fy2z_S_C1004+ABY*I_ESP_D2z_S_C1004;
    Double I_ESP_D2x_Pz_C1004 = I_ESP_F2xz_S_C1004+ABZ*I_ESP_D2x_S_C1004;
    Double I_ESP_Dxy_Pz_C1004 = I_ESP_Fxyz_S_C1004+ABZ*I_ESP_Dxy_S_C1004;
    Double I_ESP_Dxz_Pz_C1004 = I_ESP_Fx2z_S_C1004+ABZ*I_ESP_Dxz_S_C1004;
    Double I_ESP_D2y_Pz_C1004 = I_ESP_F2yz_S_C1004+ABZ*I_ESP_D2y_S_C1004;
    Double I_ESP_Dyz_Pz_C1004 = I_ESP_Fy2z_S_C1004+ABZ*I_ESP_Dyz_S_C1004;
    Double I_ESP_D2z_Pz_C1004 = I_ESP_F3z_S_C1004+ABZ*I_ESP_D2z_S_C1004;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_C1004_a
     * RHS shell quartet name: SQ_ESP_G_S_C1004_a
     ************************************************************/
    Double I_ESP_G4x_Px_C1004_a = I_ESP_H5x_S_C1004_a+ABX*I_ESP_G4x_S_C1004_a;
    Double I_ESP_G3xy_Px_C1004_a = I_ESP_H4xy_S_C1004_a+ABX*I_ESP_G3xy_S_C1004_a;
    Double I_ESP_G3xz_Px_C1004_a = I_ESP_H4xz_S_C1004_a+ABX*I_ESP_G3xz_S_C1004_a;
    Double I_ESP_G2x2y_Px_C1004_a = I_ESP_H3x2y_S_C1004_a+ABX*I_ESP_G2x2y_S_C1004_a;
    Double I_ESP_G2xyz_Px_C1004_a = I_ESP_H3xyz_S_C1004_a+ABX*I_ESP_G2xyz_S_C1004_a;
    Double I_ESP_G2x2z_Px_C1004_a = I_ESP_H3x2z_S_C1004_a+ABX*I_ESP_G2x2z_S_C1004_a;
    Double I_ESP_Gx3y_Px_C1004_a = I_ESP_H2x3y_S_C1004_a+ABX*I_ESP_Gx3y_S_C1004_a;
    Double I_ESP_Gx2yz_Px_C1004_a = I_ESP_H2x2yz_S_C1004_a+ABX*I_ESP_Gx2yz_S_C1004_a;
    Double I_ESP_Gxy2z_Px_C1004_a = I_ESP_H2xy2z_S_C1004_a+ABX*I_ESP_Gxy2z_S_C1004_a;
    Double I_ESP_Gx3z_Px_C1004_a = I_ESP_H2x3z_S_C1004_a+ABX*I_ESP_Gx3z_S_C1004_a;
    Double I_ESP_G4y_Px_C1004_a = I_ESP_Hx4y_S_C1004_a+ABX*I_ESP_G4y_S_C1004_a;
    Double I_ESP_G3yz_Px_C1004_a = I_ESP_Hx3yz_S_C1004_a+ABX*I_ESP_G3yz_S_C1004_a;
    Double I_ESP_G2y2z_Px_C1004_a = I_ESP_Hx2y2z_S_C1004_a+ABX*I_ESP_G2y2z_S_C1004_a;
    Double I_ESP_Gy3z_Px_C1004_a = I_ESP_Hxy3z_S_C1004_a+ABX*I_ESP_Gy3z_S_C1004_a;
    Double I_ESP_G4z_Px_C1004_a = I_ESP_Hx4z_S_C1004_a+ABX*I_ESP_G4z_S_C1004_a;
    Double I_ESP_G4x_Py_C1004_a = I_ESP_H4xy_S_C1004_a+ABY*I_ESP_G4x_S_C1004_a;
    Double I_ESP_G3xy_Py_C1004_a = I_ESP_H3x2y_S_C1004_a+ABY*I_ESP_G3xy_S_C1004_a;
    Double I_ESP_G3xz_Py_C1004_a = I_ESP_H3xyz_S_C1004_a+ABY*I_ESP_G3xz_S_C1004_a;
    Double I_ESP_G2x2y_Py_C1004_a = I_ESP_H2x3y_S_C1004_a+ABY*I_ESP_G2x2y_S_C1004_a;
    Double I_ESP_G2xyz_Py_C1004_a = I_ESP_H2x2yz_S_C1004_a+ABY*I_ESP_G2xyz_S_C1004_a;
    Double I_ESP_G2x2z_Py_C1004_a = I_ESP_H2xy2z_S_C1004_a+ABY*I_ESP_G2x2z_S_C1004_a;
    Double I_ESP_Gx3y_Py_C1004_a = I_ESP_Hx4y_S_C1004_a+ABY*I_ESP_Gx3y_S_C1004_a;
    Double I_ESP_Gx2yz_Py_C1004_a = I_ESP_Hx3yz_S_C1004_a+ABY*I_ESP_Gx2yz_S_C1004_a;
    Double I_ESP_Gxy2z_Py_C1004_a = I_ESP_Hx2y2z_S_C1004_a+ABY*I_ESP_Gxy2z_S_C1004_a;
    Double I_ESP_Gx3z_Py_C1004_a = I_ESP_Hxy3z_S_C1004_a+ABY*I_ESP_Gx3z_S_C1004_a;
    Double I_ESP_G4y_Py_C1004_a = I_ESP_H5y_S_C1004_a+ABY*I_ESP_G4y_S_C1004_a;
    Double I_ESP_G3yz_Py_C1004_a = I_ESP_H4yz_S_C1004_a+ABY*I_ESP_G3yz_S_C1004_a;
    Double I_ESP_G2y2z_Py_C1004_a = I_ESP_H3y2z_S_C1004_a+ABY*I_ESP_G2y2z_S_C1004_a;
    Double I_ESP_Gy3z_Py_C1004_a = I_ESP_H2y3z_S_C1004_a+ABY*I_ESP_Gy3z_S_C1004_a;
    Double I_ESP_G4z_Py_C1004_a = I_ESP_Hy4z_S_C1004_a+ABY*I_ESP_G4z_S_C1004_a;
    Double I_ESP_G4x_Pz_C1004_a = I_ESP_H4xz_S_C1004_a+ABZ*I_ESP_G4x_S_C1004_a;
    Double I_ESP_G3xy_Pz_C1004_a = I_ESP_H3xyz_S_C1004_a+ABZ*I_ESP_G3xy_S_C1004_a;
    Double I_ESP_G3xz_Pz_C1004_a = I_ESP_H3x2z_S_C1004_a+ABZ*I_ESP_G3xz_S_C1004_a;
    Double I_ESP_G2x2y_Pz_C1004_a = I_ESP_H2x2yz_S_C1004_a+ABZ*I_ESP_G2x2y_S_C1004_a;
    Double I_ESP_G2xyz_Pz_C1004_a = I_ESP_H2xy2z_S_C1004_a+ABZ*I_ESP_G2xyz_S_C1004_a;
    Double I_ESP_G2x2z_Pz_C1004_a = I_ESP_H2x3z_S_C1004_a+ABZ*I_ESP_G2x2z_S_C1004_a;
    Double I_ESP_Gx3y_Pz_C1004_a = I_ESP_Hx3yz_S_C1004_a+ABZ*I_ESP_Gx3y_S_C1004_a;
    Double I_ESP_Gx2yz_Pz_C1004_a = I_ESP_Hx2y2z_S_C1004_a+ABZ*I_ESP_Gx2yz_S_C1004_a;
    Double I_ESP_Gxy2z_Pz_C1004_a = I_ESP_Hxy3z_S_C1004_a+ABZ*I_ESP_Gxy2z_S_C1004_a;
    Double I_ESP_Gx3z_Pz_C1004_a = I_ESP_Hx4z_S_C1004_a+ABZ*I_ESP_Gx3z_S_C1004_a;
    Double I_ESP_G4y_Pz_C1004_a = I_ESP_H4yz_S_C1004_a+ABZ*I_ESP_G4y_S_C1004_a;
    Double I_ESP_G3yz_Pz_C1004_a = I_ESP_H3y2z_S_C1004_a+ABZ*I_ESP_G3yz_S_C1004_a;
    Double I_ESP_G2y2z_Pz_C1004_a = I_ESP_H2y3z_S_C1004_a+ABZ*I_ESP_G2y2z_S_C1004_a;
    Double I_ESP_Gy3z_Pz_C1004_a = I_ESP_Hy4z_S_C1004_a+ABZ*I_ESP_Gy3z_S_C1004_a;
    Double I_ESP_G4z_Pz_C1004_a = I_ESP_H5z_S_C1004_a+ABZ*I_ESP_G4z_S_C1004_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P_C1004_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S_C1004_aa
     * RHS shell quartet name: SQ_ESP_I_S_C1004_aa
     ************************************************************/
    Double I_ESP_I6x_Px_C1004_aa = I_ESP_K7x_S_C1004_aa+ABX*I_ESP_I6x_S_C1004_aa;
    Double I_ESP_I5xy_Px_C1004_aa = I_ESP_K6xy_S_C1004_aa+ABX*I_ESP_I5xy_S_C1004_aa;
    Double I_ESP_I5xz_Px_C1004_aa = I_ESP_K6xz_S_C1004_aa+ABX*I_ESP_I5xz_S_C1004_aa;
    Double I_ESP_I4x2y_Px_C1004_aa = I_ESP_K5x2y_S_C1004_aa+ABX*I_ESP_I4x2y_S_C1004_aa;
    Double I_ESP_I4xyz_Px_C1004_aa = I_ESP_K5xyz_S_C1004_aa+ABX*I_ESP_I4xyz_S_C1004_aa;
    Double I_ESP_I4x2z_Px_C1004_aa = I_ESP_K5x2z_S_C1004_aa+ABX*I_ESP_I4x2z_S_C1004_aa;
    Double I_ESP_I3x3y_Px_C1004_aa = I_ESP_K4x3y_S_C1004_aa+ABX*I_ESP_I3x3y_S_C1004_aa;
    Double I_ESP_I3x2yz_Px_C1004_aa = I_ESP_K4x2yz_S_C1004_aa+ABX*I_ESP_I3x2yz_S_C1004_aa;
    Double I_ESP_I3xy2z_Px_C1004_aa = I_ESP_K4xy2z_S_C1004_aa+ABX*I_ESP_I3xy2z_S_C1004_aa;
    Double I_ESP_I3x3z_Px_C1004_aa = I_ESP_K4x3z_S_C1004_aa+ABX*I_ESP_I3x3z_S_C1004_aa;
    Double I_ESP_I2x4y_Px_C1004_aa = I_ESP_K3x4y_S_C1004_aa+ABX*I_ESP_I2x4y_S_C1004_aa;
    Double I_ESP_I2x3yz_Px_C1004_aa = I_ESP_K3x3yz_S_C1004_aa+ABX*I_ESP_I2x3yz_S_C1004_aa;
    Double I_ESP_I2x2y2z_Px_C1004_aa = I_ESP_K3x2y2z_S_C1004_aa+ABX*I_ESP_I2x2y2z_S_C1004_aa;
    Double I_ESP_I2xy3z_Px_C1004_aa = I_ESP_K3xy3z_S_C1004_aa+ABX*I_ESP_I2xy3z_S_C1004_aa;
    Double I_ESP_I2x4z_Px_C1004_aa = I_ESP_K3x4z_S_C1004_aa+ABX*I_ESP_I2x4z_S_C1004_aa;
    Double I_ESP_Ix5y_Px_C1004_aa = I_ESP_K2x5y_S_C1004_aa+ABX*I_ESP_Ix5y_S_C1004_aa;
    Double I_ESP_Ix4yz_Px_C1004_aa = I_ESP_K2x4yz_S_C1004_aa+ABX*I_ESP_Ix4yz_S_C1004_aa;
    Double I_ESP_Ix3y2z_Px_C1004_aa = I_ESP_K2x3y2z_S_C1004_aa+ABX*I_ESP_Ix3y2z_S_C1004_aa;
    Double I_ESP_Ix2y3z_Px_C1004_aa = I_ESP_K2x2y3z_S_C1004_aa+ABX*I_ESP_Ix2y3z_S_C1004_aa;
    Double I_ESP_Ixy4z_Px_C1004_aa = I_ESP_K2xy4z_S_C1004_aa+ABX*I_ESP_Ixy4z_S_C1004_aa;
    Double I_ESP_Ix5z_Px_C1004_aa = I_ESP_K2x5z_S_C1004_aa+ABX*I_ESP_Ix5z_S_C1004_aa;
    Double I_ESP_I6y_Px_C1004_aa = I_ESP_Kx6y_S_C1004_aa+ABX*I_ESP_I6y_S_C1004_aa;
    Double I_ESP_I5yz_Px_C1004_aa = I_ESP_Kx5yz_S_C1004_aa+ABX*I_ESP_I5yz_S_C1004_aa;
    Double I_ESP_I4y2z_Px_C1004_aa = I_ESP_Kx4y2z_S_C1004_aa+ABX*I_ESP_I4y2z_S_C1004_aa;
    Double I_ESP_I3y3z_Px_C1004_aa = I_ESP_Kx3y3z_S_C1004_aa+ABX*I_ESP_I3y3z_S_C1004_aa;
    Double I_ESP_I2y4z_Px_C1004_aa = I_ESP_Kx2y4z_S_C1004_aa+ABX*I_ESP_I2y4z_S_C1004_aa;
    Double I_ESP_Iy5z_Px_C1004_aa = I_ESP_Kxy5z_S_C1004_aa+ABX*I_ESP_Iy5z_S_C1004_aa;
    Double I_ESP_I6z_Px_C1004_aa = I_ESP_Kx6z_S_C1004_aa+ABX*I_ESP_I6z_S_C1004_aa;
    Double I_ESP_I6x_Py_C1004_aa = I_ESP_K6xy_S_C1004_aa+ABY*I_ESP_I6x_S_C1004_aa;
    Double I_ESP_I5xy_Py_C1004_aa = I_ESP_K5x2y_S_C1004_aa+ABY*I_ESP_I5xy_S_C1004_aa;
    Double I_ESP_I5xz_Py_C1004_aa = I_ESP_K5xyz_S_C1004_aa+ABY*I_ESP_I5xz_S_C1004_aa;
    Double I_ESP_I4x2y_Py_C1004_aa = I_ESP_K4x3y_S_C1004_aa+ABY*I_ESP_I4x2y_S_C1004_aa;
    Double I_ESP_I4xyz_Py_C1004_aa = I_ESP_K4x2yz_S_C1004_aa+ABY*I_ESP_I4xyz_S_C1004_aa;
    Double I_ESP_I4x2z_Py_C1004_aa = I_ESP_K4xy2z_S_C1004_aa+ABY*I_ESP_I4x2z_S_C1004_aa;
    Double I_ESP_I3x3y_Py_C1004_aa = I_ESP_K3x4y_S_C1004_aa+ABY*I_ESP_I3x3y_S_C1004_aa;
    Double I_ESP_I3x2yz_Py_C1004_aa = I_ESP_K3x3yz_S_C1004_aa+ABY*I_ESP_I3x2yz_S_C1004_aa;
    Double I_ESP_I3xy2z_Py_C1004_aa = I_ESP_K3x2y2z_S_C1004_aa+ABY*I_ESP_I3xy2z_S_C1004_aa;
    Double I_ESP_I3x3z_Py_C1004_aa = I_ESP_K3xy3z_S_C1004_aa+ABY*I_ESP_I3x3z_S_C1004_aa;
    Double I_ESP_I2x4y_Py_C1004_aa = I_ESP_K2x5y_S_C1004_aa+ABY*I_ESP_I2x4y_S_C1004_aa;
    Double I_ESP_I2x3yz_Py_C1004_aa = I_ESP_K2x4yz_S_C1004_aa+ABY*I_ESP_I2x3yz_S_C1004_aa;
    Double I_ESP_I2x2y2z_Py_C1004_aa = I_ESP_K2x3y2z_S_C1004_aa+ABY*I_ESP_I2x2y2z_S_C1004_aa;
    Double I_ESP_I2xy3z_Py_C1004_aa = I_ESP_K2x2y3z_S_C1004_aa+ABY*I_ESP_I2xy3z_S_C1004_aa;
    Double I_ESP_I2x4z_Py_C1004_aa = I_ESP_K2xy4z_S_C1004_aa+ABY*I_ESP_I2x4z_S_C1004_aa;
    Double I_ESP_Ix5y_Py_C1004_aa = I_ESP_Kx6y_S_C1004_aa+ABY*I_ESP_Ix5y_S_C1004_aa;
    Double I_ESP_Ix4yz_Py_C1004_aa = I_ESP_Kx5yz_S_C1004_aa+ABY*I_ESP_Ix4yz_S_C1004_aa;
    Double I_ESP_Ix3y2z_Py_C1004_aa = I_ESP_Kx4y2z_S_C1004_aa+ABY*I_ESP_Ix3y2z_S_C1004_aa;
    Double I_ESP_Ix2y3z_Py_C1004_aa = I_ESP_Kx3y3z_S_C1004_aa+ABY*I_ESP_Ix2y3z_S_C1004_aa;
    Double I_ESP_Ixy4z_Py_C1004_aa = I_ESP_Kx2y4z_S_C1004_aa+ABY*I_ESP_Ixy4z_S_C1004_aa;
    Double I_ESP_Ix5z_Py_C1004_aa = I_ESP_Kxy5z_S_C1004_aa+ABY*I_ESP_Ix5z_S_C1004_aa;
    Double I_ESP_I6y_Py_C1004_aa = I_ESP_K7y_S_C1004_aa+ABY*I_ESP_I6y_S_C1004_aa;
    Double I_ESP_I5yz_Py_C1004_aa = I_ESP_K6yz_S_C1004_aa+ABY*I_ESP_I5yz_S_C1004_aa;
    Double I_ESP_I4y2z_Py_C1004_aa = I_ESP_K5y2z_S_C1004_aa+ABY*I_ESP_I4y2z_S_C1004_aa;
    Double I_ESP_I3y3z_Py_C1004_aa = I_ESP_K4y3z_S_C1004_aa+ABY*I_ESP_I3y3z_S_C1004_aa;
    Double I_ESP_I2y4z_Py_C1004_aa = I_ESP_K3y4z_S_C1004_aa+ABY*I_ESP_I2y4z_S_C1004_aa;
    Double I_ESP_Iy5z_Py_C1004_aa = I_ESP_K2y5z_S_C1004_aa+ABY*I_ESP_Iy5z_S_C1004_aa;
    Double I_ESP_I6z_Py_C1004_aa = I_ESP_Ky6z_S_C1004_aa+ABY*I_ESP_I6z_S_C1004_aa;
    Double I_ESP_I6x_Pz_C1004_aa = I_ESP_K6xz_S_C1004_aa+ABZ*I_ESP_I6x_S_C1004_aa;
    Double I_ESP_I5xy_Pz_C1004_aa = I_ESP_K5xyz_S_C1004_aa+ABZ*I_ESP_I5xy_S_C1004_aa;
    Double I_ESP_I5xz_Pz_C1004_aa = I_ESP_K5x2z_S_C1004_aa+ABZ*I_ESP_I5xz_S_C1004_aa;
    Double I_ESP_I4x2y_Pz_C1004_aa = I_ESP_K4x2yz_S_C1004_aa+ABZ*I_ESP_I4x2y_S_C1004_aa;
    Double I_ESP_I4xyz_Pz_C1004_aa = I_ESP_K4xy2z_S_C1004_aa+ABZ*I_ESP_I4xyz_S_C1004_aa;
    Double I_ESP_I4x2z_Pz_C1004_aa = I_ESP_K4x3z_S_C1004_aa+ABZ*I_ESP_I4x2z_S_C1004_aa;
    Double I_ESP_I3x3y_Pz_C1004_aa = I_ESP_K3x3yz_S_C1004_aa+ABZ*I_ESP_I3x3y_S_C1004_aa;
    Double I_ESP_I3x2yz_Pz_C1004_aa = I_ESP_K3x2y2z_S_C1004_aa+ABZ*I_ESP_I3x2yz_S_C1004_aa;
    Double I_ESP_I3xy2z_Pz_C1004_aa = I_ESP_K3xy3z_S_C1004_aa+ABZ*I_ESP_I3xy2z_S_C1004_aa;
    Double I_ESP_I3x3z_Pz_C1004_aa = I_ESP_K3x4z_S_C1004_aa+ABZ*I_ESP_I3x3z_S_C1004_aa;
    Double I_ESP_I2x4y_Pz_C1004_aa = I_ESP_K2x4yz_S_C1004_aa+ABZ*I_ESP_I2x4y_S_C1004_aa;
    Double I_ESP_I2x3yz_Pz_C1004_aa = I_ESP_K2x3y2z_S_C1004_aa+ABZ*I_ESP_I2x3yz_S_C1004_aa;
    Double I_ESP_I2x2y2z_Pz_C1004_aa = I_ESP_K2x2y3z_S_C1004_aa+ABZ*I_ESP_I2x2y2z_S_C1004_aa;
    Double I_ESP_I2xy3z_Pz_C1004_aa = I_ESP_K2xy4z_S_C1004_aa+ABZ*I_ESP_I2xy3z_S_C1004_aa;
    Double I_ESP_I2x4z_Pz_C1004_aa = I_ESP_K2x5z_S_C1004_aa+ABZ*I_ESP_I2x4z_S_C1004_aa;
    Double I_ESP_Ix5y_Pz_C1004_aa = I_ESP_Kx5yz_S_C1004_aa+ABZ*I_ESP_Ix5y_S_C1004_aa;
    Double I_ESP_Ix4yz_Pz_C1004_aa = I_ESP_Kx4y2z_S_C1004_aa+ABZ*I_ESP_Ix4yz_S_C1004_aa;
    Double I_ESP_Ix3y2z_Pz_C1004_aa = I_ESP_Kx3y3z_S_C1004_aa+ABZ*I_ESP_Ix3y2z_S_C1004_aa;
    Double I_ESP_Ix2y3z_Pz_C1004_aa = I_ESP_Kx2y4z_S_C1004_aa+ABZ*I_ESP_Ix2y3z_S_C1004_aa;
    Double I_ESP_Ixy4z_Pz_C1004_aa = I_ESP_Kxy5z_S_C1004_aa+ABZ*I_ESP_Ixy4z_S_C1004_aa;
    Double I_ESP_Ix5z_Pz_C1004_aa = I_ESP_Kx6z_S_C1004_aa+ABZ*I_ESP_Ix5z_S_C1004_aa;
    Double I_ESP_I6y_Pz_C1004_aa = I_ESP_K6yz_S_C1004_aa+ABZ*I_ESP_I6y_S_C1004_aa;
    Double I_ESP_I5yz_Pz_C1004_aa = I_ESP_K5y2z_S_C1004_aa+ABZ*I_ESP_I5yz_S_C1004_aa;
    Double I_ESP_I4y2z_Pz_C1004_aa = I_ESP_K4y3z_S_C1004_aa+ABZ*I_ESP_I4y2z_S_C1004_aa;
    Double I_ESP_I3y3z_Pz_C1004_aa = I_ESP_K3y4z_S_C1004_aa+ABZ*I_ESP_I3y3z_S_C1004_aa;
    Double I_ESP_I2y4z_Pz_C1004_aa = I_ESP_K2y5z_S_C1004_aa+ABZ*I_ESP_I2y4z_S_C1004_aa;
    Double I_ESP_Iy5z_Pz_C1004_aa = I_ESP_Ky6z_S_C1004_aa+ABZ*I_ESP_Iy5z_S_C1004_aa;
    Double I_ESP_I6z_Pz_C1004_aa = I_ESP_K7z_S_C1004_aa+ABZ*I_ESP_I6z_S_C1004_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+0] = 4.0E0*I_ESP_I6x_S_C4_aa-2.0E0*4*I_ESP_G4x_S_C4_a-2.0E0*5*I_ESP_G4x_S_C4_a+4*3*I_ESP_D2x_S_C4;
    abcd[iGrid*360+1] = 4.0E0*I_ESP_I5xy_S_C4_aa-2.0E0*3*I_ESP_G3xy_S_C4_a-2.0E0*4*I_ESP_G3xy_S_C4_a+3*2*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+2] = 4.0E0*I_ESP_I5xz_S_C4_aa-2.0E0*3*I_ESP_G3xz_S_C4_a-2.0E0*4*I_ESP_G3xz_S_C4_a+3*2*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+3] = 4.0E0*I_ESP_I4x2y_S_C4_aa-2.0E0*2*I_ESP_G2x2y_S_C4_a-2.0E0*3*I_ESP_G2x2y_S_C4_a+2*1*I_ESP_D2y_S_C4;
    abcd[iGrid*360+4] = 4.0E0*I_ESP_I4xyz_S_C4_aa-2.0E0*2*I_ESP_G2xyz_S_C4_a-2.0E0*3*I_ESP_G2xyz_S_C4_a+2*1*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+5] = 4.0E0*I_ESP_I4x2z_S_C4_aa-2.0E0*2*I_ESP_G2x2z_S_C4_a-2.0E0*3*I_ESP_G2x2z_S_C4_a+2*1*I_ESP_D2z_S_C4;
    abcd[iGrid*360+6] = 4.0E0*I_ESP_I3x3y_S_C4_aa-2.0E0*1*I_ESP_Gx3y_S_C4_a-2.0E0*2*I_ESP_Gx3y_S_C4_a;
    abcd[iGrid*360+7] = 4.0E0*I_ESP_I3x2yz_S_C4_aa-2.0E0*1*I_ESP_Gx2yz_S_C4_a-2.0E0*2*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+8] = 4.0E0*I_ESP_I3xy2z_S_C4_aa-2.0E0*1*I_ESP_Gxy2z_S_C4_a-2.0E0*2*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+9] = 4.0E0*I_ESP_I3x3z_S_C4_aa-2.0E0*1*I_ESP_Gx3z_S_C4_a-2.0E0*2*I_ESP_Gx3z_S_C4_a;
    abcd[iGrid*360+10] = 4.0E0*I_ESP_I2x4y_S_C4_aa-2.0E0*1*I_ESP_G4y_S_C4_a;
    abcd[iGrid*360+11] = 4.0E0*I_ESP_I2x3yz_S_C4_aa-2.0E0*1*I_ESP_G3yz_S_C4_a;
    abcd[iGrid*360+12] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2y2z_S_C4_a;
    abcd[iGrid*360+13] = 4.0E0*I_ESP_I2xy3z_S_C4_aa-2.0E0*1*I_ESP_Gy3z_S_C4_a;
    abcd[iGrid*360+14] = 4.0E0*I_ESP_I2x4z_S_C4_aa-2.0E0*1*I_ESP_G4z_S_C4_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+15] = 4.0E0*I_ESP_I6x_Px_C1004_aa-2.0E0*4*I_ESP_G4x_Px_C1004_a-2.0E0*5*I_ESP_G4x_Px_C1004_a+4*3*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+16] = 4.0E0*I_ESP_I5xy_Px_C1004_aa-2.0E0*3*I_ESP_G3xy_Px_C1004_a-2.0E0*4*I_ESP_G3xy_Px_C1004_a+3*2*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+17] = 4.0E0*I_ESP_I5xz_Px_C1004_aa-2.0E0*3*I_ESP_G3xz_Px_C1004_a-2.0E0*4*I_ESP_G3xz_Px_C1004_a+3*2*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+18] = 4.0E0*I_ESP_I4x2y_Px_C1004_aa-2.0E0*2*I_ESP_G2x2y_Px_C1004_a-2.0E0*3*I_ESP_G2x2y_Px_C1004_a+2*1*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+19] = 4.0E0*I_ESP_I4xyz_Px_C1004_aa-2.0E0*2*I_ESP_G2xyz_Px_C1004_a-2.0E0*3*I_ESP_G2xyz_Px_C1004_a+2*1*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+20] = 4.0E0*I_ESP_I4x2z_Px_C1004_aa-2.0E0*2*I_ESP_G2x2z_Px_C1004_a-2.0E0*3*I_ESP_G2x2z_Px_C1004_a+2*1*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+21] = 4.0E0*I_ESP_I3x3y_Px_C1004_aa-2.0E0*1*I_ESP_Gx3y_Px_C1004_a-2.0E0*2*I_ESP_Gx3y_Px_C1004_a;
    abcd[iGrid*360+22] = 4.0E0*I_ESP_I3x2yz_Px_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Px_C1004_a-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+23] = 4.0E0*I_ESP_I3xy2z_Px_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Px_C1004_a-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+24] = 4.0E0*I_ESP_I3x3z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3z_Px_C1004_a-2.0E0*2*I_ESP_Gx3z_Px_C1004_a;
    abcd[iGrid*360+25] = 4.0E0*I_ESP_I2x4y_Px_C1004_aa-2.0E0*1*I_ESP_G4y_Px_C1004_a;
    abcd[iGrid*360+26] = 4.0E0*I_ESP_I2x3yz_Px_C1004_aa-2.0E0*1*I_ESP_G3yz_Px_C1004_a;
    abcd[iGrid*360+27] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2y2z_Px_C1004_a;
    abcd[iGrid*360+28] = 4.0E0*I_ESP_I2xy3z_Px_C1004_aa-2.0E0*1*I_ESP_Gy3z_Px_C1004_a;
    abcd[iGrid*360+29] = 4.0E0*I_ESP_I2x4z_Px_C1004_aa-2.0E0*1*I_ESP_G4z_Px_C1004_a;
    abcd[iGrid*360+30] = 4.0E0*I_ESP_I6x_Py_C1004_aa-2.0E0*4*I_ESP_G4x_Py_C1004_a-2.0E0*5*I_ESP_G4x_Py_C1004_a+4*3*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+31] = 4.0E0*I_ESP_I5xy_Py_C1004_aa-2.0E0*3*I_ESP_G3xy_Py_C1004_a-2.0E0*4*I_ESP_G3xy_Py_C1004_a+3*2*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+32] = 4.0E0*I_ESP_I5xz_Py_C1004_aa-2.0E0*3*I_ESP_G3xz_Py_C1004_a-2.0E0*4*I_ESP_G3xz_Py_C1004_a+3*2*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+33] = 4.0E0*I_ESP_I4x2y_Py_C1004_aa-2.0E0*2*I_ESP_G2x2y_Py_C1004_a-2.0E0*3*I_ESP_G2x2y_Py_C1004_a+2*1*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+34] = 4.0E0*I_ESP_I4xyz_Py_C1004_aa-2.0E0*2*I_ESP_G2xyz_Py_C1004_a-2.0E0*3*I_ESP_G2xyz_Py_C1004_a+2*1*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+35] = 4.0E0*I_ESP_I4x2z_Py_C1004_aa-2.0E0*2*I_ESP_G2x2z_Py_C1004_a-2.0E0*3*I_ESP_G2x2z_Py_C1004_a+2*1*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+36] = 4.0E0*I_ESP_I3x3y_Py_C1004_aa-2.0E0*1*I_ESP_Gx3y_Py_C1004_a-2.0E0*2*I_ESP_Gx3y_Py_C1004_a;
    abcd[iGrid*360+37] = 4.0E0*I_ESP_I3x2yz_Py_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Py_C1004_a-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+38] = 4.0E0*I_ESP_I3xy2z_Py_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Py_C1004_a-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+39] = 4.0E0*I_ESP_I3x3z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3z_Py_C1004_a-2.0E0*2*I_ESP_Gx3z_Py_C1004_a;
    abcd[iGrid*360+40] = 4.0E0*I_ESP_I2x4y_Py_C1004_aa-2.0E0*1*I_ESP_G4y_Py_C1004_a;
    abcd[iGrid*360+41] = 4.0E0*I_ESP_I2x3yz_Py_C1004_aa-2.0E0*1*I_ESP_G3yz_Py_C1004_a;
    abcd[iGrid*360+42] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2y2z_Py_C1004_a;
    abcd[iGrid*360+43] = 4.0E0*I_ESP_I2xy3z_Py_C1004_aa-2.0E0*1*I_ESP_Gy3z_Py_C1004_a;
    abcd[iGrid*360+44] = 4.0E0*I_ESP_I2x4z_Py_C1004_aa-2.0E0*1*I_ESP_G4z_Py_C1004_a;
    abcd[iGrid*360+45] = 4.0E0*I_ESP_I6x_Pz_C1004_aa-2.0E0*4*I_ESP_G4x_Pz_C1004_a-2.0E0*5*I_ESP_G4x_Pz_C1004_a+4*3*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+46] = 4.0E0*I_ESP_I5xy_Pz_C1004_aa-2.0E0*3*I_ESP_G3xy_Pz_C1004_a-2.0E0*4*I_ESP_G3xy_Pz_C1004_a+3*2*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+47] = 4.0E0*I_ESP_I5xz_Pz_C1004_aa-2.0E0*3*I_ESP_G3xz_Pz_C1004_a-2.0E0*4*I_ESP_G3xz_Pz_C1004_a+3*2*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+48] = 4.0E0*I_ESP_I4x2y_Pz_C1004_aa-2.0E0*2*I_ESP_G2x2y_Pz_C1004_a-2.0E0*3*I_ESP_G2x2y_Pz_C1004_a+2*1*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+49] = 4.0E0*I_ESP_I4xyz_Pz_C1004_aa-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a-2.0E0*3*I_ESP_G2xyz_Pz_C1004_a+2*1*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+50] = 4.0E0*I_ESP_I4x2z_Pz_C1004_aa-2.0E0*2*I_ESP_G2x2z_Pz_C1004_a-2.0E0*3*I_ESP_G2x2z_Pz_C1004_a+2*1*I_ESP_D2z_Pz_C1004;
    abcd[iGrid*360+51] = 4.0E0*I_ESP_I3x3y_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3y_Pz_C1004_a-2.0E0*2*I_ESP_Gx3y_Pz_C1004_a;
    abcd[iGrid*360+52] = 4.0E0*I_ESP_I3x2yz_Pz_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Pz_C1004_a-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+53] = 4.0E0*I_ESP_I3xy2z_Pz_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Pz_C1004_a-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+54] = 4.0E0*I_ESP_I3x3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3z_Pz_C1004_a-2.0E0*2*I_ESP_Gx3z_Pz_C1004_a;
    abcd[iGrid*360+55] = 4.0E0*I_ESP_I2x4y_Pz_C1004_aa-2.0E0*1*I_ESP_G4y_Pz_C1004_a;
    abcd[iGrid*360+56] = 4.0E0*I_ESP_I2x3yz_Pz_C1004_aa-2.0E0*1*I_ESP_G3yz_Pz_C1004_a;
    abcd[iGrid*360+57] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2y2z_Pz_C1004_a;
    abcd[iGrid*360+58] = 4.0E0*I_ESP_I2xy3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gy3z_Pz_C1004_a;
    abcd[iGrid*360+59] = 4.0E0*I_ESP_I2x4z_Pz_C1004_aa-2.0E0*1*I_ESP_G4z_Pz_C1004_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+60] = 4.0E0*I_ESP_I5xy_S_C4_aa-2.0E0*4*I_ESP_G3xy_S_C4_a;
    abcd[iGrid*360+61] = 4.0E0*I_ESP_I4x2y_S_C4_aa-2.0E0*1*I_ESP_G4x_S_C4_a-2.0E0*3*I_ESP_G2x2y_S_C4_a+3*1*I_ESP_D2x_S_C4;
    abcd[iGrid*360+62] = 4.0E0*I_ESP_I4xyz_S_C4_aa-2.0E0*3*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+63] = 4.0E0*I_ESP_I3x3y_S_C4_aa-2.0E0*2*I_ESP_G3xy_S_C4_a-2.0E0*2*I_ESP_Gx3y_S_C4_a+2*2*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+64] = 4.0E0*I_ESP_I3x2yz_S_C4_aa-2.0E0*1*I_ESP_G3xz_S_C4_a-2.0E0*2*I_ESP_Gx2yz_S_C4_a+2*1*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+65] = 4.0E0*I_ESP_I3xy2z_S_C4_aa-2.0E0*2*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+66] = 4.0E0*I_ESP_I2x4y_S_C4_aa-2.0E0*3*I_ESP_G2x2y_S_C4_a-2.0E0*1*I_ESP_G4y_S_C4_a+3*I_ESP_D2y_S_C4;
    abcd[iGrid*360+67] = 4.0E0*I_ESP_I2x3yz_S_C4_aa-2.0E0*2*I_ESP_G2xyz_S_C4_a-2.0E0*1*I_ESP_G3yz_S_C4_a+2*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+68] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2x2z_S_C4_a-2.0E0*1*I_ESP_G2y2z_S_C4_a+1*I_ESP_D2z_S_C4;
    abcd[iGrid*360+69] = 4.0E0*I_ESP_I2xy3z_S_C4_aa-2.0E0*1*I_ESP_Gy3z_S_C4_a;
    abcd[iGrid*360+70] = 4.0E0*I_ESP_Ix5y_S_C4_aa-2.0E0*4*I_ESP_Gx3y_S_C4_a;
    abcd[iGrid*360+71] = 4.0E0*I_ESP_Ix4yz_S_C4_aa-2.0E0*3*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+72] = 4.0E0*I_ESP_Ix3y2z_S_C4_aa-2.0E0*2*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+73] = 4.0E0*I_ESP_Ix2y3z_S_C4_aa-2.0E0*1*I_ESP_Gx3z_S_C4_a;
    abcd[iGrid*360+74] = 4.0E0*I_ESP_Ixy4z_S_C4_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+75] = 4.0E0*I_ESP_I5xy_Px_C1004_aa-2.0E0*4*I_ESP_G3xy_Px_C1004_a;
    abcd[iGrid*360+76] = 4.0E0*I_ESP_I4x2y_Px_C1004_aa-2.0E0*1*I_ESP_G4x_Px_C1004_a-2.0E0*3*I_ESP_G2x2y_Px_C1004_a+3*1*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+77] = 4.0E0*I_ESP_I4xyz_Px_C1004_aa-2.0E0*3*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+78] = 4.0E0*I_ESP_I3x3y_Px_C1004_aa-2.0E0*2*I_ESP_G3xy_Px_C1004_a-2.0E0*2*I_ESP_Gx3y_Px_C1004_a+2*2*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+79] = 4.0E0*I_ESP_I3x2yz_Px_C1004_aa-2.0E0*1*I_ESP_G3xz_Px_C1004_a-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a+2*1*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+80] = 4.0E0*I_ESP_I3xy2z_Px_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+81] = 4.0E0*I_ESP_I2x4y_Px_C1004_aa-2.0E0*3*I_ESP_G2x2y_Px_C1004_a-2.0E0*1*I_ESP_G4y_Px_C1004_a+3*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+82] = 4.0E0*I_ESP_I2x3yz_Px_C1004_aa-2.0E0*2*I_ESP_G2xyz_Px_C1004_a-2.0E0*1*I_ESP_G3yz_Px_C1004_a+2*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+83] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2x2z_Px_C1004_a-2.0E0*1*I_ESP_G2y2z_Px_C1004_a+1*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+84] = 4.0E0*I_ESP_I2xy3z_Px_C1004_aa-2.0E0*1*I_ESP_Gy3z_Px_C1004_a;
    abcd[iGrid*360+85] = 4.0E0*I_ESP_Ix5y_Px_C1004_aa-2.0E0*4*I_ESP_Gx3y_Px_C1004_a;
    abcd[iGrid*360+86] = 4.0E0*I_ESP_Ix4yz_Px_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+87] = 4.0E0*I_ESP_Ix3y2z_Px_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+88] = 4.0E0*I_ESP_Ix2y3z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3z_Px_C1004_a;
    abcd[iGrid*360+89] = 4.0E0*I_ESP_Ixy4z_Px_C1004_aa;
    abcd[iGrid*360+90] = 4.0E0*I_ESP_I5xy_Py_C1004_aa-2.0E0*4*I_ESP_G3xy_Py_C1004_a;
    abcd[iGrid*360+91] = 4.0E0*I_ESP_I4x2y_Py_C1004_aa-2.0E0*1*I_ESP_G4x_Py_C1004_a-2.0E0*3*I_ESP_G2x2y_Py_C1004_a+3*1*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+92] = 4.0E0*I_ESP_I4xyz_Py_C1004_aa-2.0E0*3*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+93] = 4.0E0*I_ESP_I3x3y_Py_C1004_aa-2.0E0*2*I_ESP_G3xy_Py_C1004_a-2.0E0*2*I_ESP_Gx3y_Py_C1004_a+2*2*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+94] = 4.0E0*I_ESP_I3x2yz_Py_C1004_aa-2.0E0*1*I_ESP_G3xz_Py_C1004_a-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a+2*1*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+95] = 4.0E0*I_ESP_I3xy2z_Py_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+96] = 4.0E0*I_ESP_I2x4y_Py_C1004_aa-2.0E0*3*I_ESP_G2x2y_Py_C1004_a-2.0E0*1*I_ESP_G4y_Py_C1004_a+3*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+97] = 4.0E0*I_ESP_I2x3yz_Py_C1004_aa-2.0E0*2*I_ESP_G2xyz_Py_C1004_a-2.0E0*1*I_ESP_G3yz_Py_C1004_a+2*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+98] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2x2z_Py_C1004_a-2.0E0*1*I_ESP_G2y2z_Py_C1004_a+1*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+99] = 4.0E0*I_ESP_I2xy3z_Py_C1004_aa-2.0E0*1*I_ESP_Gy3z_Py_C1004_a;
    abcd[iGrid*360+100] = 4.0E0*I_ESP_Ix5y_Py_C1004_aa-2.0E0*4*I_ESP_Gx3y_Py_C1004_a;
    abcd[iGrid*360+101] = 4.0E0*I_ESP_Ix4yz_Py_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+102] = 4.0E0*I_ESP_Ix3y2z_Py_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+103] = 4.0E0*I_ESP_Ix2y3z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3z_Py_C1004_a;
    abcd[iGrid*360+104] = 4.0E0*I_ESP_Ixy4z_Py_C1004_aa;
    abcd[iGrid*360+105] = 4.0E0*I_ESP_I5xy_Pz_C1004_aa-2.0E0*4*I_ESP_G3xy_Pz_C1004_a;
    abcd[iGrid*360+106] = 4.0E0*I_ESP_I4x2y_Pz_C1004_aa-2.0E0*1*I_ESP_G4x_Pz_C1004_a-2.0E0*3*I_ESP_G2x2y_Pz_C1004_a+3*1*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+107] = 4.0E0*I_ESP_I4xyz_Pz_C1004_aa-2.0E0*3*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+108] = 4.0E0*I_ESP_I3x3y_Pz_C1004_aa-2.0E0*2*I_ESP_G3xy_Pz_C1004_a-2.0E0*2*I_ESP_Gx3y_Pz_C1004_a+2*2*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+109] = 4.0E0*I_ESP_I3x2yz_Pz_C1004_aa-2.0E0*1*I_ESP_G3xz_Pz_C1004_a-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a+2*1*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+110] = 4.0E0*I_ESP_I3xy2z_Pz_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+111] = 4.0E0*I_ESP_I2x4y_Pz_C1004_aa-2.0E0*3*I_ESP_G2x2y_Pz_C1004_a-2.0E0*1*I_ESP_G4y_Pz_C1004_a+3*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+112] = 4.0E0*I_ESP_I2x3yz_Pz_C1004_aa-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a-2.0E0*1*I_ESP_G3yz_Pz_C1004_a+2*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+113] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2x2z_Pz_C1004_a-2.0E0*1*I_ESP_G2y2z_Pz_C1004_a+1*I_ESP_D2z_Pz_C1004;
    abcd[iGrid*360+114] = 4.0E0*I_ESP_I2xy3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gy3z_Pz_C1004_a;
    abcd[iGrid*360+115] = 4.0E0*I_ESP_Ix5y_Pz_C1004_aa-2.0E0*4*I_ESP_Gx3y_Pz_C1004_a;
    abcd[iGrid*360+116] = 4.0E0*I_ESP_Ix4yz_Pz_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+117] = 4.0E0*I_ESP_Ix3y2z_Pz_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+118] = 4.0E0*I_ESP_Ix2y3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3z_Pz_C1004_a;
    abcd[iGrid*360+119] = 4.0E0*I_ESP_Ixy4z_Pz_C1004_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+120] = 4.0E0*I_ESP_I5xz_S_C4_aa-2.0E0*4*I_ESP_G3xz_S_C4_a;
    abcd[iGrid*360+121] = 4.0E0*I_ESP_I4xyz_S_C4_aa-2.0E0*3*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+122] = 4.0E0*I_ESP_I4x2z_S_C4_aa-2.0E0*1*I_ESP_G4x_S_C4_a-2.0E0*3*I_ESP_G2x2z_S_C4_a+3*1*I_ESP_D2x_S_C4;
    abcd[iGrid*360+123] = 4.0E0*I_ESP_I3x2yz_S_C4_aa-2.0E0*2*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+124] = 4.0E0*I_ESP_I3xy2z_S_C4_aa-2.0E0*1*I_ESP_G3xy_S_C4_a-2.0E0*2*I_ESP_Gxy2z_S_C4_a+2*1*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+125] = 4.0E0*I_ESP_I3x3z_S_C4_aa-2.0E0*2*I_ESP_G3xz_S_C4_a-2.0E0*2*I_ESP_Gx3z_S_C4_a+2*2*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+126] = 4.0E0*I_ESP_I2x3yz_S_C4_aa-2.0E0*1*I_ESP_G3yz_S_C4_a;
    abcd[iGrid*360+127] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2x2y_S_C4_a-2.0E0*1*I_ESP_G2y2z_S_C4_a+1*I_ESP_D2y_S_C4;
    abcd[iGrid*360+128] = 4.0E0*I_ESP_I2xy3z_S_C4_aa-2.0E0*2*I_ESP_G2xyz_S_C4_a-2.0E0*1*I_ESP_Gy3z_S_C4_a+2*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+129] = 4.0E0*I_ESP_I2x4z_S_C4_aa-2.0E0*3*I_ESP_G2x2z_S_C4_a-2.0E0*1*I_ESP_G4z_S_C4_a+3*I_ESP_D2z_S_C4;
    abcd[iGrid*360+130] = 4.0E0*I_ESP_Ix4yz_S_C4_aa;
    abcd[iGrid*360+131] = 4.0E0*I_ESP_Ix3y2z_S_C4_aa-2.0E0*1*I_ESP_Gx3y_S_C4_a;
    abcd[iGrid*360+132] = 4.0E0*I_ESP_Ix2y3z_S_C4_aa-2.0E0*2*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+133] = 4.0E0*I_ESP_Ixy4z_S_C4_aa-2.0E0*3*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+134] = 4.0E0*I_ESP_Ix5z_S_C4_aa-2.0E0*4*I_ESP_Gx3z_S_C4_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+135] = 4.0E0*I_ESP_I5xz_Px_C1004_aa-2.0E0*4*I_ESP_G3xz_Px_C1004_a;
    abcd[iGrid*360+136] = 4.0E0*I_ESP_I4xyz_Px_C1004_aa-2.0E0*3*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+137] = 4.0E0*I_ESP_I4x2z_Px_C1004_aa-2.0E0*1*I_ESP_G4x_Px_C1004_a-2.0E0*3*I_ESP_G2x2z_Px_C1004_a+3*1*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+138] = 4.0E0*I_ESP_I3x2yz_Px_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+139] = 4.0E0*I_ESP_I3xy2z_Px_C1004_aa-2.0E0*1*I_ESP_G3xy_Px_C1004_a-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a+2*1*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+140] = 4.0E0*I_ESP_I3x3z_Px_C1004_aa-2.0E0*2*I_ESP_G3xz_Px_C1004_a-2.0E0*2*I_ESP_Gx3z_Px_C1004_a+2*2*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+141] = 4.0E0*I_ESP_I2x3yz_Px_C1004_aa-2.0E0*1*I_ESP_G3yz_Px_C1004_a;
    abcd[iGrid*360+142] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2x2y_Px_C1004_a-2.0E0*1*I_ESP_G2y2z_Px_C1004_a+1*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+143] = 4.0E0*I_ESP_I2xy3z_Px_C1004_aa-2.0E0*2*I_ESP_G2xyz_Px_C1004_a-2.0E0*1*I_ESP_Gy3z_Px_C1004_a+2*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+144] = 4.0E0*I_ESP_I2x4z_Px_C1004_aa-2.0E0*3*I_ESP_G2x2z_Px_C1004_a-2.0E0*1*I_ESP_G4z_Px_C1004_a+3*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+145] = 4.0E0*I_ESP_Ix4yz_Px_C1004_aa;
    abcd[iGrid*360+146] = 4.0E0*I_ESP_Ix3y2z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3y_Px_C1004_a;
    abcd[iGrid*360+147] = 4.0E0*I_ESP_Ix2y3z_Px_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+148] = 4.0E0*I_ESP_Ixy4z_Px_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+149] = 4.0E0*I_ESP_Ix5z_Px_C1004_aa-2.0E0*4*I_ESP_Gx3z_Px_C1004_a;
    abcd[iGrid*360+150] = 4.0E0*I_ESP_I5xz_Py_C1004_aa-2.0E0*4*I_ESP_G3xz_Py_C1004_a;
    abcd[iGrid*360+151] = 4.0E0*I_ESP_I4xyz_Py_C1004_aa-2.0E0*3*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+152] = 4.0E0*I_ESP_I4x2z_Py_C1004_aa-2.0E0*1*I_ESP_G4x_Py_C1004_a-2.0E0*3*I_ESP_G2x2z_Py_C1004_a+3*1*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+153] = 4.0E0*I_ESP_I3x2yz_Py_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+154] = 4.0E0*I_ESP_I3xy2z_Py_C1004_aa-2.0E0*1*I_ESP_G3xy_Py_C1004_a-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a+2*1*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+155] = 4.0E0*I_ESP_I3x3z_Py_C1004_aa-2.0E0*2*I_ESP_G3xz_Py_C1004_a-2.0E0*2*I_ESP_Gx3z_Py_C1004_a+2*2*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+156] = 4.0E0*I_ESP_I2x3yz_Py_C1004_aa-2.0E0*1*I_ESP_G3yz_Py_C1004_a;
    abcd[iGrid*360+157] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2x2y_Py_C1004_a-2.0E0*1*I_ESP_G2y2z_Py_C1004_a+1*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+158] = 4.0E0*I_ESP_I2xy3z_Py_C1004_aa-2.0E0*2*I_ESP_G2xyz_Py_C1004_a-2.0E0*1*I_ESP_Gy3z_Py_C1004_a+2*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+159] = 4.0E0*I_ESP_I2x4z_Py_C1004_aa-2.0E0*3*I_ESP_G2x2z_Py_C1004_a-2.0E0*1*I_ESP_G4z_Py_C1004_a+3*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+160] = 4.0E0*I_ESP_Ix4yz_Py_C1004_aa;
    abcd[iGrid*360+161] = 4.0E0*I_ESP_Ix3y2z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3y_Py_C1004_a;
    abcd[iGrid*360+162] = 4.0E0*I_ESP_Ix2y3z_Py_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+163] = 4.0E0*I_ESP_Ixy4z_Py_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+164] = 4.0E0*I_ESP_Ix5z_Py_C1004_aa-2.0E0*4*I_ESP_Gx3z_Py_C1004_a;
    abcd[iGrid*360+165] = 4.0E0*I_ESP_I5xz_Pz_C1004_aa-2.0E0*4*I_ESP_G3xz_Pz_C1004_a;
    abcd[iGrid*360+166] = 4.0E0*I_ESP_I4xyz_Pz_C1004_aa-2.0E0*3*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+167] = 4.0E0*I_ESP_I4x2z_Pz_C1004_aa-2.0E0*1*I_ESP_G4x_Pz_C1004_a-2.0E0*3*I_ESP_G2x2z_Pz_C1004_a+3*1*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+168] = 4.0E0*I_ESP_I3x2yz_Pz_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+169] = 4.0E0*I_ESP_I3xy2z_Pz_C1004_aa-2.0E0*1*I_ESP_G3xy_Pz_C1004_a-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a+2*1*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+170] = 4.0E0*I_ESP_I3x3z_Pz_C1004_aa-2.0E0*2*I_ESP_G3xz_Pz_C1004_a-2.0E0*2*I_ESP_Gx3z_Pz_C1004_a+2*2*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+171] = 4.0E0*I_ESP_I2x3yz_Pz_C1004_aa-2.0E0*1*I_ESP_G3yz_Pz_C1004_a;
    abcd[iGrid*360+172] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2x2y_Pz_C1004_a-2.0E0*1*I_ESP_G2y2z_Pz_C1004_a+1*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+173] = 4.0E0*I_ESP_I2xy3z_Pz_C1004_aa-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a-2.0E0*1*I_ESP_Gy3z_Pz_C1004_a+2*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+174] = 4.0E0*I_ESP_I2x4z_Pz_C1004_aa-2.0E0*3*I_ESP_G2x2z_Pz_C1004_a-2.0E0*1*I_ESP_G4z_Pz_C1004_a+3*I_ESP_D2z_Pz_C1004;
    abcd[iGrid*360+175] = 4.0E0*I_ESP_Ix4yz_Pz_C1004_aa;
    abcd[iGrid*360+176] = 4.0E0*I_ESP_Ix3y2z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3y_Pz_C1004_a;
    abcd[iGrid*360+177] = 4.0E0*I_ESP_Ix2y3z_Pz_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+178] = 4.0E0*I_ESP_Ixy4z_Pz_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+179] = 4.0E0*I_ESP_Ix5z_Pz_C1004_aa-2.0E0*4*I_ESP_Gx3z_Pz_C1004_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+180] = 4.0E0*I_ESP_I4x2y_S_C4_aa-2.0E0*1*I_ESP_G4x_S_C4_a;
    abcd[iGrid*360+181] = 4.0E0*I_ESP_I3x3y_S_C4_aa-2.0E0*1*I_ESP_G3xy_S_C4_a-2.0E0*2*I_ESP_G3xy_S_C4_a;
    abcd[iGrid*360+182] = 4.0E0*I_ESP_I3x2yz_S_C4_aa-2.0E0*1*I_ESP_G3xz_S_C4_a;
    abcd[iGrid*360+183] = 4.0E0*I_ESP_I2x4y_S_C4_aa-2.0E0*2*I_ESP_G2x2y_S_C4_a-2.0E0*3*I_ESP_G2x2y_S_C4_a+2*1*I_ESP_D2x_S_C4;
    abcd[iGrid*360+184] = 4.0E0*I_ESP_I2x3yz_S_C4_aa-2.0E0*1*I_ESP_G2xyz_S_C4_a-2.0E0*2*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+185] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2x2z_S_C4_a;
    abcd[iGrid*360+186] = 4.0E0*I_ESP_Ix5y_S_C4_aa-2.0E0*3*I_ESP_Gx3y_S_C4_a-2.0E0*4*I_ESP_Gx3y_S_C4_a+3*2*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+187] = 4.0E0*I_ESP_Ix4yz_S_C4_aa-2.0E0*2*I_ESP_Gx2yz_S_C4_a-2.0E0*3*I_ESP_Gx2yz_S_C4_a+2*1*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+188] = 4.0E0*I_ESP_Ix3y2z_S_C4_aa-2.0E0*1*I_ESP_Gxy2z_S_C4_a-2.0E0*2*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+189] = 4.0E0*I_ESP_Ix2y3z_S_C4_aa-2.0E0*1*I_ESP_Gx3z_S_C4_a;
    abcd[iGrid*360+190] = 4.0E0*I_ESP_I6y_S_C4_aa-2.0E0*4*I_ESP_G4y_S_C4_a-2.0E0*5*I_ESP_G4y_S_C4_a+4*3*I_ESP_D2y_S_C4;
    abcd[iGrid*360+191] = 4.0E0*I_ESP_I5yz_S_C4_aa-2.0E0*3*I_ESP_G3yz_S_C4_a-2.0E0*4*I_ESP_G3yz_S_C4_a+3*2*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+192] = 4.0E0*I_ESP_I4y2z_S_C4_aa-2.0E0*2*I_ESP_G2y2z_S_C4_a-2.0E0*3*I_ESP_G2y2z_S_C4_a+2*1*I_ESP_D2z_S_C4;
    abcd[iGrid*360+193] = 4.0E0*I_ESP_I3y3z_S_C4_aa-2.0E0*1*I_ESP_Gy3z_S_C4_a-2.0E0*2*I_ESP_Gy3z_S_C4_a;
    abcd[iGrid*360+194] = 4.0E0*I_ESP_I2y4z_S_C4_aa-2.0E0*1*I_ESP_G4z_S_C4_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+195] = 4.0E0*I_ESP_I4x2y_Px_C1004_aa-2.0E0*1*I_ESP_G4x_Px_C1004_a;
    abcd[iGrid*360+196] = 4.0E0*I_ESP_I3x3y_Px_C1004_aa-2.0E0*1*I_ESP_G3xy_Px_C1004_a-2.0E0*2*I_ESP_G3xy_Px_C1004_a;
    abcd[iGrid*360+197] = 4.0E0*I_ESP_I3x2yz_Px_C1004_aa-2.0E0*1*I_ESP_G3xz_Px_C1004_a;
    abcd[iGrid*360+198] = 4.0E0*I_ESP_I2x4y_Px_C1004_aa-2.0E0*2*I_ESP_G2x2y_Px_C1004_a-2.0E0*3*I_ESP_G2x2y_Px_C1004_a+2*1*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+199] = 4.0E0*I_ESP_I2x3yz_Px_C1004_aa-2.0E0*1*I_ESP_G2xyz_Px_C1004_a-2.0E0*2*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+200] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2x2z_Px_C1004_a;
    abcd[iGrid*360+201] = 4.0E0*I_ESP_Ix5y_Px_C1004_aa-2.0E0*3*I_ESP_Gx3y_Px_C1004_a-2.0E0*4*I_ESP_Gx3y_Px_C1004_a+3*2*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+202] = 4.0E0*I_ESP_Ix4yz_Px_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a-2.0E0*3*I_ESP_Gx2yz_Px_C1004_a+2*1*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+203] = 4.0E0*I_ESP_Ix3y2z_Px_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Px_C1004_a-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+204] = 4.0E0*I_ESP_Ix2y3z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3z_Px_C1004_a;
    abcd[iGrid*360+205] = 4.0E0*I_ESP_I6y_Px_C1004_aa-2.0E0*4*I_ESP_G4y_Px_C1004_a-2.0E0*5*I_ESP_G4y_Px_C1004_a+4*3*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+206] = 4.0E0*I_ESP_I5yz_Px_C1004_aa-2.0E0*3*I_ESP_G3yz_Px_C1004_a-2.0E0*4*I_ESP_G3yz_Px_C1004_a+3*2*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+207] = 4.0E0*I_ESP_I4y2z_Px_C1004_aa-2.0E0*2*I_ESP_G2y2z_Px_C1004_a-2.0E0*3*I_ESP_G2y2z_Px_C1004_a+2*1*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+208] = 4.0E0*I_ESP_I3y3z_Px_C1004_aa-2.0E0*1*I_ESP_Gy3z_Px_C1004_a-2.0E0*2*I_ESP_Gy3z_Px_C1004_a;
    abcd[iGrid*360+209] = 4.0E0*I_ESP_I2y4z_Px_C1004_aa-2.0E0*1*I_ESP_G4z_Px_C1004_a;
    abcd[iGrid*360+210] = 4.0E0*I_ESP_I4x2y_Py_C1004_aa-2.0E0*1*I_ESP_G4x_Py_C1004_a;
    abcd[iGrid*360+211] = 4.0E0*I_ESP_I3x3y_Py_C1004_aa-2.0E0*1*I_ESP_G3xy_Py_C1004_a-2.0E0*2*I_ESP_G3xy_Py_C1004_a;
    abcd[iGrid*360+212] = 4.0E0*I_ESP_I3x2yz_Py_C1004_aa-2.0E0*1*I_ESP_G3xz_Py_C1004_a;
    abcd[iGrid*360+213] = 4.0E0*I_ESP_I2x4y_Py_C1004_aa-2.0E0*2*I_ESP_G2x2y_Py_C1004_a-2.0E0*3*I_ESP_G2x2y_Py_C1004_a+2*1*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+214] = 4.0E0*I_ESP_I2x3yz_Py_C1004_aa-2.0E0*1*I_ESP_G2xyz_Py_C1004_a-2.0E0*2*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+215] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2x2z_Py_C1004_a;
    abcd[iGrid*360+216] = 4.0E0*I_ESP_Ix5y_Py_C1004_aa-2.0E0*3*I_ESP_Gx3y_Py_C1004_a-2.0E0*4*I_ESP_Gx3y_Py_C1004_a+3*2*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+217] = 4.0E0*I_ESP_Ix4yz_Py_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a-2.0E0*3*I_ESP_Gx2yz_Py_C1004_a+2*1*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+218] = 4.0E0*I_ESP_Ix3y2z_Py_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Py_C1004_a-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+219] = 4.0E0*I_ESP_Ix2y3z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3z_Py_C1004_a;
    abcd[iGrid*360+220] = 4.0E0*I_ESP_I6y_Py_C1004_aa-2.0E0*4*I_ESP_G4y_Py_C1004_a-2.0E0*5*I_ESP_G4y_Py_C1004_a+4*3*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+221] = 4.0E0*I_ESP_I5yz_Py_C1004_aa-2.0E0*3*I_ESP_G3yz_Py_C1004_a-2.0E0*4*I_ESP_G3yz_Py_C1004_a+3*2*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+222] = 4.0E0*I_ESP_I4y2z_Py_C1004_aa-2.0E0*2*I_ESP_G2y2z_Py_C1004_a-2.0E0*3*I_ESP_G2y2z_Py_C1004_a+2*1*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+223] = 4.0E0*I_ESP_I3y3z_Py_C1004_aa-2.0E0*1*I_ESP_Gy3z_Py_C1004_a-2.0E0*2*I_ESP_Gy3z_Py_C1004_a;
    abcd[iGrid*360+224] = 4.0E0*I_ESP_I2y4z_Py_C1004_aa-2.0E0*1*I_ESP_G4z_Py_C1004_a;
    abcd[iGrid*360+225] = 4.0E0*I_ESP_I4x2y_Pz_C1004_aa-2.0E0*1*I_ESP_G4x_Pz_C1004_a;
    abcd[iGrid*360+226] = 4.0E0*I_ESP_I3x3y_Pz_C1004_aa-2.0E0*1*I_ESP_G3xy_Pz_C1004_a-2.0E0*2*I_ESP_G3xy_Pz_C1004_a;
    abcd[iGrid*360+227] = 4.0E0*I_ESP_I3x2yz_Pz_C1004_aa-2.0E0*1*I_ESP_G3xz_Pz_C1004_a;
    abcd[iGrid*360+228] = 4.0E0*I_ESP_I2x4y_Pz_C1004_aa-2.0E0*2*I_ESP_G2x2y_Pz_C1004_a-2.0E0*3*I_ESP_G2x2y_Pz_C1004_a+2*1*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+229] = 4.0E0*I_ESP_I2x3yz_Pz_C1004_aa-2.0E0*1*I_ESP_G2xyz_Pz_C1004_a-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+230] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2x2z_Pz_C1004_a;
    abcd[iGrid*360+231] = 4.0E0*I_ESP_Ix5y_Pz_C1004_aa-2.0E0*3*I_ESP_Gx3y_Pz_C1004_a-2.0E0*4*I_ESP_Gx3y_Pz_C1004_a+3*2*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+232] = 4.0E0*I_ESP_Ix4yz_Pz_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a-2.0E0*3*I_ESP_Gx2yz_Pz_C1004_a+2*1*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+233] = 4.0E0*I_ESP_Ix3y2z_Pz_C1004_aa-2.0E0*1*I_ESP_Gxy2z_Pz_C1004_a-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+234] = 4.0E0*I_ESP_Ix2y3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3z_Pz_C1004_a;
    abcd[iGrid*360+235] = 4.0E0*I_ESP_I6y_Pz_C1004_aa-2.0E0*4*I_ESP_G4y_Pz_C1004_a-2.0E0*5*I_ESP_G4y_Pz_C1004_a+4*3*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+236] = 4.0E0*I_ESP_I5yz_Pz_C1004_aa-2.0E0*3*I_ESP_G3yz_Pz_C1004_a-2.0E0*4*I_ESP_G3yz_Pz_C1004_a+3*2*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+237] = 4.0E0*I_ESP_I4y2z_Pz_C1004_aa-2.0E0*2*I_ESP_G2y2z_Pz_C1004_a-2.0E0*3*I_ESP_G2y2z_Pz_C1004_a+2*1*I_ESP_D2z_Pz_C1004;
    abcd[iGrid*360+238] = 4.0E0*I_ESP_I3y3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gy3z_Pz_C1004_a-2.0E0*2*I_ESP_Gy3z_Pz_C1004_a;
    abcd[iGrid*360+239] = 4.0E0*I_ESP_I2y4z_Pz_C1004_aa-2.0E0*1*I_ESP_G4z_Pz_C1004_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+240] = 4.0E0*I_ESP_I4xyz_S_C4_aa;
    abcd[iGrid*360+241] = 4.0E0*I_ESP_I3x2yz_S_C4_aa-2.0E0*1*I_ESP_G3xz_S_C4_a;
    abcd[iGrid*360+242] = 4.0E0*I_ESP_I3xy2z_S_C4_aa-2.0E0*1*I_ESP_G3xy_S_C4_a;
    abcd[iGrid*360+243] = 4.0E0*I_ESP_I2x3yz_S_C4_aa-2.0E0*2*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+244] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2x2y_S_C4_a-2.0E0*1*I_ESP_G2x2z_S_C4_a+1*I_ESP_D2x_S_C4;
    abcd[iGrid*360+245] = 4.0E0*I_ESP_I2xy3z_S_C4_aa-2.0E0*2*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+246] = 4.0E0*I_ESP_Ix4yz_S_C4_aa-2.0E0*3*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+247] = 4.0E0*I_ESP_Ix3y2z_S_C4_aa-2.0E0*1*I_ESP_Gx3y_S_C4_a-2.0E0*2*I_ESP_Gxy2z_S_C4_a+2*1*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+248] = 4.0E0*I_ESP_Ix2y3z_S_C4_aa-2.0E0*2*I_ESP_Gx2yz_S_C4_a-2.0E0*1*I_ESP_Gx3z_S_C4_a+2*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+249] = 4.0E0*I_ESP_Ixy4z_S_C4_aa-2.0E0*3*I_ESP_Gxy2z_S_C4_a;
    abcd[iGrid*360+250] = 4.0E0*I_ESP_I5yz_S_C4_aa-2.0E0*4*I_ESP_G3yz_S_C4_a;
    abcd[iGrid*360+251] = 4.0E0*I_ESP_I4y2z_S_C4_aa-2.0E0*1*I_ESP_G4y_S_C4_a-2.0E0*3*I_ESP_G2y2z_S_C4_a+3*1*I_ESP_D2y_S_C4;
    abcd[iGrid*360+252] = 4.0E0*I_ESP_I3y3z_S_C4_aa-2.0E0*2*I_ESP_G3yz_S_C4_a-2.0E0*2*I_ESP_Gy3z_S_C4_a+2*2*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+253] = 4.0E0*I_ESP_I2y4z_S_C4_aa-2.0E0*3*I_ESP_G2y2z_S_C4_a-2.0E0*1*I_ESP_G4z_S_C4_a+3*I_ESP_D2z_S_C4;
    abcd[iGrid*360+254] = 4.0E0*I_ESP_Iy5z_S_C4_aa-2.0E0*4*I_ESP_Gy3z_S_C4_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+255] = 4.0E0*I_ESP_I4xyz_Px_C1004_aa;
    abcd[iGrid*360+256] = 4.0E0*I_ESP_I3x2yz_Px_C1004_aa-2.0E0*1*I_ESP_G3xz_Px_C1004_a;
    abcd[iGrid*360+257] = 4.0E0*I_ESP_I3xy2z_Px_C1004_aa-2.0E0*1*I_ESP_G3xy_Px_C1004_a;
    abcd[iGrid*360+258] = 4.0E0*I_ESP_I2x3yz_Px_C1004_aa-2.0E0*2*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+259] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2x2y_Px_C1004_a-2.0E0*1*I_ESP_G2x2z_Px_C1004_a+1*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+260] = 4.0E0*I_ESP_I2xy3z_Px_C1004_aa-2.0E0*2*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+261] = 4.0E0*I_ESP_Ix4yz_Px_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+262] = 4.0E0*I_ESP_Ix3y2z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3y_Px_C1004_a-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a+2*1*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+263] = 4.0E0*I_ESP_Ix2y3z_Px_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a-2.0E0*1*I_ESP_Gx3z_Px_C1004_a+2*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+264] = 4.0E0*I_ESP_Ixy4z_Px_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Px_C1004_a;
    abcd[iGrid*360+265] = 4.0E0*I_ESP_I5yz_Px_C1004_aa-2.0E0*4*I_ESP_G3yz_Px_C1004_a;
    abcd[iGrid*360+266] = 4.0E0*I_ESP_I4y2z_Px_C1004_aa-2.0E0*1*I_ESP_G4y_Px_C1004_a-2.0E0*3*I_ESP_G2y2z_Px_C1004_a+3*1*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+267] = 4.0E0*I_ESP_I3y3z_Px_C1004_aa-2.0E0*2*I_ESP_G3yz_Px_C1004_a-2.0E0*2*I_ESP_Gy3z_Px_C1004_a+2*2*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+268] = 4.0E0*I_ESP_I2y4z_Px_C1004_aa-2.0E0*3*I_ESP_G2y2z_Px_C1004_a-2.0E0*1*I_ESP_G4z_Px_C1004_a+3*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+269] = 4.0E0*I_ESP_Iy5z_Px_C1004_aa-2.0E0*4*I_ESP_Gy3z_Px_C1004_a;
    abcd[iGrid*360+270] = 4.0E0*I_ESP_I4xyz_Py_C1004_aa;
    abcd[iGrid*360+271] = 4.0E0*I_ESP_I3x2yz_Py_C1004_aa-2.0E0*1*I_ESP_G3xz_Py_C1004_a;
    abcd[iGrid*360+272] = 4.0E0*I_ESP_I3xy2z_Py_C1004_aa-2.0E0*1*I_ESP_G3xy_Py_C1004_a;
    abcd[iGrid*360+273] = 4.0E0*I_ESP_I2x3yz_Py_C1004_aa-2.0E0*2*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+274] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2x2y_Py_C1004_a-2.0E0*1*I_ESP_G2x2z_Py_C1004_a+1*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+275] = 4.0E0*I_ESP_I2xy3z_Py_C1004_aa-2.0E0*2*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+276] = 4.0E0*I_ESP_Ix4yz_Py_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+277] = 4.0E0*I_ESP_Ix3y2z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3y_Py_C1004_a-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a+2*1*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+278] = 4.0E0*I_ESP_Ix2y3z_Py_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a-2.0E0*1*I_ESP_Gx3z_Py_C1004_a+2*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+279] = 4.0E0*I_ESP_Ixy4z_Py_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Py_C1004_a;
    abcd[iGrid*360+280] = 4.0E0*I_ESP_I5yz_Py_C1004_aa-2.0E0*4*I_ESP_G3yz_Py_C1004_a;
    abcd[iGrid*360+281] = 4.0E0*I_ESP_I4y2z_Py_C1004_aa-2.0E0*1*I_ESP_G4y_Py_C1004_a-2.0E0*3*I_ESP_G2y2z_Py_C1004_a+3*1*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+282] = 4.0E0*I_ESP_I3y3z_Py_C1004_aa-2.0E0*2*I_ESP_G3yz_Py_C1004_a-2.0E0*2*I_ESP_Gy3z_Py_C1004_a+2*2*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+283] = 4.0E0*I_ESP_I2y4z_Py_C1004_aa-2.0E0*3*I_ESP_G2y2z_Py_C1004_a-2.0E0*1*I_ESP_G4z_Py_C1004_a+3*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+284] = 4.0E0*I_ESP_Iy5z_Py_C1004_aa-2.0E0*4*I_ESP_Gy3z_Py_C1004_a;
    abcd[iGrid*360+285] = 4.0E0*I_ESP_I4xyz_Pz_C1004_aa;
    abcd[iGrid*360+286] = 4.0E0*I_ESP_I3x2yz_Pz_C1004_aa-2.0E0*1*I_ESP_G3xz_Pz_C1004_a;
    abcd[iGrid*360+287] = 4.0E0*I_ESP_I3xy2z_Pz_C1004_aa-2.0E0*1*I_ESP_G3xy_Pz_C1004_a;
    abcd[iGrid*360+288] = 4.0E0*I_ESP_I2x3yz_Pz_C1004_aa-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+289] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2x2y_Pz_C1004_a-2.0E0*1*I_ESP_G2x2z_Pz_C1004_a+1*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+290] = 4.0E0*I_ESP_I2xy3z_Pz_C1004_aa-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+291] = 4.0E0*I_ESP_Ix4yz_Pz_C1004_aa-2.0E0*3*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+292] = 4.0E0*I_ESP_Ix3y2z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3y_Pz_C1004_a-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a+2*1*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+293] = 4.0E0*I_ESP_Ix2y3z_Pz_C1004_aa-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a-2.0E0*1*I_ESP_Gx3z_Pz_C1004_a+2*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+294] = 4.0E0*I_ESP_Ixy4z_Pz_C1004_aa-2.0E0*3*I_ESP_Gxy2z_Pz_C1004_a;
    abcd[iGrid*360+295] = 4.0E0*I_ESP_I5yz_Pz_C1004_aa-2.0E0*4*I_ESP_G3yz_Pz_C1004_a;
    abcd[iGrid*360+296] = 4.0E0*I_ESP_I4y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G4y_Pz_C1004_a-2.0E0*3*I_ESP_G2y2z_Pz_C1004_a+3*1*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+297] = 4.0E0*I_ESP_I3y3z_Pz_C1004_aa-2.0E0*2*I_ESP_G3yz_Pz_C1004_a-2.0E0*2*I_ESP_Gy3z_Pz_C1004_a+2*2*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+298] = 4.0E0*I_ESP_I2y4z_Pz_C1004_aa-2.0E0*3*I_ESP_G2y2z_Pz_C1004_a-2.0E0*1*I_ESP_G4z_Pz_C1004_a+3*I_ESP_D2z_Pz_C1004;
    abcd[iGrid*360+299] = 4.0E0*I_ESP_Iy5z_Pz_C1004_aa-2.0E0*4*I_ESP_Gy3z_Pz_C1004_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_S_C4_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_C4_aa
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_G_S_C4_a
     * RHS shell quartet name: SQ_ESP_D_S_C4
     ************************************************************/
    abcd[iGrid*360+300] = 4.0E0*I_ESP_I4x2z_S_C4_aa-2.0E0*1*I_ESP_G4x_S_C4_a;
    abcd[iGrid*360+301] = 4.0E0*I_ESP_I3xy2z_S_C4_aa-2.0E0*1*I_ESP_G3xy_S_C4_a;
    abcd[iGrid*360+302] = 4.0E0*I_ESP_I3x3z_S_C4_aa-2.0E0*1*I_ESP_G3xz_S_C4_a-2.0E0*2*I_ESP_G3xz_S_C4_a;
    abcd[iGrid*360+303] = 4.0E0*I_ESP_I2x2y2z_S_C4_aa-2.0E0*1*I_ESP_G2x2y_S_C4_a;
    abcd[iGrid*360+304] = 4.0E0*I_ESP_I2xy3z_S_C4_aa-2.0E0*1*I_ESP_G2xyz_S_C4_a-2.0E0*2*I_ESP_G2xyz_S_C4_a;
    abcd[iGrid*360+305] = 4.0E0*I_ESP_I2x4z_S_C4_aa-2.0E0*2*I_ESP_G2x2z_S_C4_a-2.0E0*3*I_ESP_G2x2z_S_C4_a+2*1*I_ESP_D2x_S_C4;
    abcd[iGrid*360+306] = 4.0E0*I_ESP_Ix3y2z_S_C4_aa-2.0E0*1*I_ESP_Gx3y_S_C4_a;
    abcd[iGrid*360+307] = 4.0E0*I_ESP_Ix2y3z_S_C4_aa-2.0E0*1*I_ESP_Gx2yz_S_C4_a-2.0E0*2*I_ESP_Gx2yz_S_C4_a;
    abcd[iGrid*360+308] = 4.0E0*I_ESP_Ixy4z_S_C4_aa-2.0E0*2*I_ESP_Gxy2z_S_C4_a-2.0E0*3*I_ESP_Gxy2z_S_C4_a+2*1*I_ESP_Dxy_S_C4;
    abcd[iGrid*360+309] = 4.0E0*I_ESP_Ix5z_S_C4_aa-2.0E0*3*I_ESP_Gx3z_S_C4_a-2.0E0*4*I_ESP_Gx3z_S_C4_a+3*2*I_ESP_Dxz_S_C4;
    abcd[iGrid*360+310] = 4.0E0*I_ESP_I4y2z_S_C4_aa-2.0E0*1*I_ESP_G4y_S_C4_a;
    abcd[iGrid*360+311] = 4.0E0*I_ESP_I3y3z_S_C4_aa-2.0E0*1*I_ESP_G3yz_S_C4_a-2.0E0*2*I_ESP_G3yz_S_C4_a;
    abcd[iGrid*360+312] = 4.0E0*I_ESP_I2y4z_S_C4_aa-2.0E0*2*I_ESP_G2y2z_S_C4_a-2.0E0*3*I_ESP_G2y2z_S_C4_a+2*1*I_ESP_D2y_S_C4;
    abcd[iGrid*360+313] = 4.0E0*I_ESP_Iy5z_S_C4_aa-2.0E0*3*I_ESP_Gy3z_S_C4_a-2.0E0*4*I_ESP_Gy3z_S_C4_a+3*2*I_ESP_Dyz_S_C4;
    abcd[iGrid*360+314] = 4.0E0*I_ESP_I6z_S_C4_aa-2.0E0*4*I_ESP_G4z_S_C4_a-2.0E0*5*I_ESP_G4z_S_C4_a+4*3*I_ESP_D2z_S_C4;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_C1004_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_C1004_aa
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_G_P_C1004_a
     * RHS shell quartet name: SQ_ESP_D_P_C1004
     ************************************************************/
    abcd[iGrid*360+315] = 4.0E0*I_ESP_I4x2z_Px_C1004_aa-2.0E0*1*I_ESP_G4x_Px_C1004_a;
    abcd[iGrid*360+316] = 4.0E0*I_ESP_I3xy2z_Px_C1004_aa-2.0E0*1*I_ESP_G3xy_Px_C1004_a;
    abcd[iGrid*360+317] = 4.0E0*I_ESP_I3x3z_Px_C1004_aa-2.0E0*1*I_ESP_G3xz_Px_C1004_a-2.0E0*2*I_ESP_G3xz_Px_C1004_a;
    abcd[iGrid*360+318] = 4.0E0*I_ESP_I2x2y2z_Px_C1004_aa-2.0E0*1*I_ESP_G2x2y_Px_C1004_a;
    abcd[iGrid*360+319] = 4.0E0*I_ESP_I2xy3z_Px_C1004_aa-2.0E0*1*I_ESP_G2xyz_Px_C1004_a-2.0E0*2*I_ESP_G2xyz_Px_C1004_a;
    abcd[iGrid*360+320] = 4.0E0*I_ESP_I2x4z_Px_C1004_aa-2.0E0*2*I_ESP_G2x2z_Px_C1004_a-2.0E0*3*I_ESP_G2x2z_Px_C1004_a+2*1*I_ESP_D2x_Px_C1004;
    abcd[iGrid*360+321] = 4.0E0*I_ESP_Ix3y2z_Px_C1004_aa-2.0E0*1*I_ESP_Gx3y_Px_C1004_a;
    abcd[iGrid*360+322] = 4.0E0*I_ESP_Ix2y3z_Px_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Px_C1004_a-2.0E0*2*I_ESP_Gx2yz_Px_C1004_a;
    abcd[iGrid*360+323] = 4.0E0*I_ESP_Ixy4z_Px_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Px_C1004_a-2.0E0*3*I_ESP_Gxy2z_Px_C1004_a+2*1*I_ESP_Dxy_Px_C1004;
    abcd[iGrid*360+324] = 4.0E0*I_ESP_Ix5z_Px_C1004_aa-2.0E0*3*I_ESP_Gx3z_Px_C1004_a-2.0E0*4*I_ESP_Gx3z_Px_C1004_a+3*2*I_ESP_Dxz_Px_C1004;
    abcd[iGrid*360+325] = 4.0E0*I_ESP_I4y2z_Px_C1004_aa-2.0E0*1*I_ESP_G4y_Px_C1004_a;
    abcd[iGrid*360+326] = 4.0E0*I_ESP_I3y3z_Px_C1004_aa-2.0E0*1*I_ESP_G3yz_Px_C1004_a-2.0E0*2*I_ESP_G3yz_Px_C1004_a;
    abcd[iGrid*360+327] = 4.0E0*I_ESP_I2y4z_Px_C1004_aa-2.0E0*2*I_ESP_G2y2z_Px_C1004_a-2.0E0*3*I_ESP_G2y2z_Px_C1004_a+2*1*I_ESP_D2y_Px_C1004;
    abcd[iGrid*360+328] = 4.0E0*I_ESP_Iy5z_Px_C1004_aa-2.0E0*3*I_ESP_Gy3z_Px_C1004_a-2.0E0*4*I_ESP_Gy3z_Px_C1004_a+3*2*I_ESP_Dyz_Px_C1004;
    abcd[iGrid*360+329] = 4.0E0*I_ESP_I6z_Px_C1004_aa-2.0E0*4*I_ESP_G4z_Px_C1004_a-2.0E0*5*I_ESP_G4z_Px_C1004_a+4*3*I_ESP_D2z_Px_C1004;
    abcd[iGrid*360+330] = 4.0E0*I_ESP_I4x2z_Py_C1004_aa-2.0E0*1*I_ESP_G4x_Py_C1004_a;
    abcd[iGrid*360+331] = 4.0E0*I_ESP_I3xy2z_Py_C1004_aa-2.0E0*1*I_ESP_G3xy_Py_C1004_a;
    abcd[iGrid*360+332] = 4.0E0*I_ESP_I3x3z_Py_C1004_aa-2.0E0*1*I_ESP_G3xz_Py_C1004_a-2.0E0*2*I_ESP_G3xz_Py_C1004_a;
    abcd[iGrid*360+333] = 4.0E0*I_ESP_I2x2y2z_Py_C1004_aa-2.0E0*1*I_ESP_G2x2y_Py_C1004_a;
    abcd[iGrid*360+334] = 4.0E0*I_ESP_I2xy3z_Py_C1004_aa-2.0E0*1*I_ESP_G2xyz_Py_C1004_a-2.0E0*2*I_ESP_G2xyz_Py_C1004_a;
    abcd[iGrid*360+335] = 4.0E0*I_ESP_I2x4z_Py_C1004_aa-2.0E0*2*I_ESP_G2x2z_Py_C1004_a-2.0E0*3*I_ESP_G2x2z_Py_C1004_a+2*1*I_ESP_D2x_Py_C1004;
    abcd[iGrid*360+336] = 4.0E0*I_ESP_Ix3y2z_Py_C1004_aa-2.0E0*1*I_ESP_Gx3y_Py_C1004_a;
    abcd[iGrid*360+337] = 4.0E0*I_ESP_Ix2y3z_Py_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Py_C1004_a-2.0E0*2*I_ESP_Gx2yz_Py_C1004_a;
    abcd[iGrid*360+338] = 4.0E0*I_ESP_Ixy4z_Py_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Py_C1004_a-2.0E0*3*I_ESP_Gxy2z_Py_C1004_a+2*1*I_ESP_Dxy_Py_C1004;
    abcd[iGrid*360+339] = 4.0E0*I_ESP_Ix5z_Py_C1004_aa-2.0E0*3*I_ESP_Gx3z_Py_C1004_a-2.0E0*4*I_ESP_Gx3z_Py_C1004_a+3*2*I_ESP_Dxz_Py_C1004;
    abcd[iGrid*360+340] = 4.0E0*I_ESP_I4y2z_Py_C1004_aa-2.0E0*1*I_ESP_G4y_Py_C1004_a;
    abcd[iGrid*360+341] = 4.0E0*I_ESP_I3y3z_Py_C1004_aa-2.0E0*1*I_ESP_G3yz_Py_C1004_a-2.0E0*2*I_ESP_G3yz_Py_C1004_a;
    abcd[iGrid*360+342] = 4.0E0*I_ESP_I2y4z_Py_C1004_aa-2.0E0*2*I_ESP_G2y2z_Py_C1004_a-2.0E0*3*I_ESP_G2y2z_Py_C1004_a+2*1*I_ESP_D2y_Py_C1004;
    abcd[iGrid*360+343] = 4.0E0*I_ESP_Iy5z_Py_C1004_aa-2.0E0*3*I_ESP_Gy3z_Py_C1004_a-2.0E0*4*I_ESP_Gy3z_Py_C1004_a+3*2*I_ESP_Dyz_Py_C1004;
    abcd[iGrid*360+344] = 4.0E0*I_ESP_I6z_Py_C1004_aa-2.0E0*4*I_ESP_G4z_Py_C1004_a-2.0E0*5*I_ESP_G4z_Py_C1004_a+4*3*I_ESP_D2z_Py_C1004;
    abcd[iGrid*360+345] = 4.0E0*I_ESP_I4x2z_Pz_C1004_aa-2.0E0*1*I_ESP_G4x_Pz_C1004_a;
    abcd[iGrid*360+346] = 4.0E0*I_ESP_I3xy2z_Pz_C1004_aa-2.0E0*1*I_ESP_G3xy_Pz_C1004_a;
    abcd[iGrid*360+347] = 4.0E0*I_ESP_I3x3z_Pz_C1004_aa-2.0E0*1*I_ESP_G3xz_Pz_C1004_a-2.0E0*2*I_ESP_G3xz_Pz_C1004_a;
    abcd[iGrid*360+348] = 4.0E0*I_ESP_I2x2y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G2x2y_Pz_C1004_a;
    abcd[iGrid*360+349] = 4.0E0*I_ESP_I2xy3z_Pz_C1004_aa-2.0E0*1*I_ESP_G2xyz_Pz_C1004_a-2.0E0*2*I_ESP_G2xyz_Pz_C1004_a;
    abcd[iGrid*360+350] = 4.0E0*I_ESP_I2x4z_Pz_C1004_aa-2.0E0*2*I_ESP_G2x2z_Pz_C1004_a-2.0E0*3*I_ESP_G2x2z_Pz_C1004_a+2*1*I_ESP_D2x_Pz_C1004;
    abcd[iGrid*360+351] = 4.0E0*I_ESP_Ix3y2z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx3y_Pz_C1004_a;
    abcd[iGrid*360+352] = 4.0E0*I_ESP_Ix2y3z_Pz_C1004_aa-2.0E0*1*I_ESP_Gx2yz_Pz_C1004_a-2.0E0*2*I_ESP_Gx2yz_Pz_C1004_a;
    abcd[iGrid*360+353] = 4.0E0*I_ESP_Ixy4z_Pz_C1004_aa-2.0E0*2*I_ESP_Gxy2z_Pz_C1004_a-2.0E0*3*I_ESP_Gxy2z_Pz_C1004_a+2*1*I_ESP_Dxy_Pz_C1004;
    abcd[iGrid*360+354] = 4.0E0*I_ESP_Ix5z_Pz_C1004_aa-2.0E0*3*I_ESP_Gx3z_Pz_C1004_a-2.0E0*4*I_ESP_Gx3z_Pz_C1004_a+3*2*I_ESP_Dxz_Pz_C1004;
    abcd[iGrid*360+355] = 4.0E0*I_ESP_I4y2z_Pz_C1004_aa-2.0E0*1*I_ESP_G4y_Pz_C1004_a;
    abcd[iGrid*360+356] = 4.0E0*I_ESP_I3y3z_Pz_C1004_aa-2.0E0*1*I_ESP_G3yz_Pz_C1004_a-2.0E0*2*I_ESP_G3yz_Pz_C1004_a;
    abcd[iGrid*360+357] = 4.0E0*I_ESP_I2y4z_Pz_C1004_aa-2.0E0*2*I_ESP_G2y2z_Pz_C1004_a-2.0E0*3*I_ESP_G2y2z_Pz_C1004_a+2*1*I_ESP_D2y_Pz_C1004;
    abcd[iGrid*360+358] = 4.0E0*I_ESP_Iy5z_Pz_C1004_aa-2.0E0*3*I_ESP_Gy3z_Pz_C1004_a-2.0E0*4*I_ESP_Gy3z_Pz_C1004_a+3*2*I_ESP_Dyz_Pz_C1004;
    abcd[iGrid*360+359] = 4.0E0*I_ESP_I6z_Pz_C1004_aa-2.0E0*4*I_ESP_G4z_Pz_C1004_a-2.0E0*5*I_ESP_G4z_Pz_C1004_a+4*3*I_ESP_D2z_Pz_C1004;
  }
}
