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
// BRA1 as redundant position, total RHS integrals evaluated as: 4687
// BRA2 as redundant position, total RHS integrals evaluated as: 3926
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

void hgp_os_esp_g_d_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
    Double I_ESP_I6x_S_aa = 0.0E0;
    Double I_ESP_I5xy_S_aa = 0.0E0;
    Double I_ESP_I5xz_S_aa = 0.0E0;
    Double I_ESP_I4x2y_S_aa = 0.0E0;
    Double I_ESP_I4xyz_S_aa = 0.0E0;
    Double I_ESP_I4x2z_S_aa = 0.0E0;
    Double I_ESP_I3x3y_S_aa = 0.0E0;
    Double I_ESP_I3x2yz_S_aa = 0.0E0;
    Double I_ESP_I3xy2z_S_aa = 0.0E0;
    Double I_ESP_I3x3z_S_aa = 0.0E0;
    Double I_ESP_I2x4y_S_aa = 0.0E0;
    Double I_ESP_I2x3yz_S_aa = 0.0E0;
    Double I_ESP_I2x2y2z_S_aa = 0.0E0;
    Double I_ESP_I2xy3z_S_aa = 0.0E0;
    Double I_ESP_I2x4z_S_aa = 0.0E0;
    Double I_ESP_Ix5y_S_aa = 0.0E0;
    Double I_ESP_Ix4yz_S_aa = 0.0E0;
    Double I_ESP_Ix3y2z_S_aa = 0.0E0;
    Double I_ESP_Ix2y3z_S_aa = 0.0E0;
    Double I_ESP_Ixy4z_S_aa = 0.0E0;
    Double I_ESP_Ix5z_S_aa = 0.0E0;
    Double I_ESP_I6y_S_aa = 0.0E0;
    Double I_ESP_I5yz_S_aa = 0.0E0;
    Double I_ESP_I4y2z_S_aa = 0.0E0;
    Double I_ESP_I3y3z_S_aa = 0.0E0;
    Double I_ESP_I2y4z_S_aa = 0.0E0;
    Double I_ESP_Iy5z_S_aa = 0.0E0;
    Double I_ESP_I6z_S_aa = 0.0E0;
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
    Double I_ESP_G4x_S_a = 0.0E0;
    Double I_ESP_G3xy_S_a = 0.0E0;
    Double I_ESP_G3xz_S_a = 0.0E0;
    Double I_ESP_G2x2y_S_a = 0.0E0;
    Double I_ESP_G2xyz_S_a = 0.0E0;
    Double I_ESP_G2x2z_S_a = 0.0E0;
    Double I_ESP_Gx3y_S_a = 0.0E0;
    Double I_ESP_Gx2yz_S_a = 0.0E0;
    Double I_ESP_Gxy2z_S_a = 0.0E0;
    Double I_ESP_Gx3z_S_a = 0.0E0;
    Double I_ESP_G4y_S_a = 0.0E0;
    Double I_ESP_G3yz_S_a = 0.0E0;
    Double I_ESP_G2y2z_S_a = 0.0E0;
    Double I_ESP_Gy3z_S_a = 0.0E0;
    Double I_ESP_G4z_S_a = 0.0E0;
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
    Double I_ESP_D2x_S = 0.0E0;
    Double I_ESP_Dxy_S = 0.0E0;
    Double I_ESP_Dxz_S = 0.0E0;
    Double I_ESP_D2y_S = 0.0E0;
    Double I_ESP_Dyz_S = 0.0E0;
    Double I_ESP_D2z_S = 0.0E0;

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
         * shell quartet name: SQ_ESP_I_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_aa_coefs = alpha*alpha;
        I_ESP_I6x_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6z_S_vrr;

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
         * shell quartet name: SQ_ESP_G_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_G_S_a_coefs = alpha;
        I_ESP_G4x_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4z_S_vrr;

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

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_D2x_S += I_ESP_D2x_S_vrr;
        I_ESP_Dxy_S += I_ESP_Dxy_S_vrr;
        I_ESP_Dxz_S += I_ESP_Dxz_S_vrr;
        I_ESP_D2y_S += I_ESP_D2y_S_vrr;
        I_ESP_Dyz_S += I_ESP_Dyz_S_vrr;
        I_ESP_D2z_S += I_ESP_D2z_S_vrr;
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
     * shell quartet name: SQ_ESP_D_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_S
     * RHS shell quartet name: SQ_ESP_D_S
     ************************************************************/
    Double I_ESP_D2x_Px = I_ESP_F3x_S+ABX*I_ESP_D2x_S;
    Double I_ESP_Dxy_Px = I_ESP_F2xy_S+ABX*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Px = I_ESP_F2xz_S+ABX*I_ESP_Dxz_S;
    Double I_ESP_D2y_Px = I_ESP_Fx2y_S+ABX*I_ESP_D2y_S;
    Double I_ESP_Dyz_Px = I_ESP_Fxyz_S+ABX*I_ESP_Dyz_S;
    Double I_ESP_D2z_Px = I_ESP_Fx2z_S+ABX*I_ESP_D2z_S;
    Double I_ESP_D2x_Py = I_ESP_F2xy_S+ABY*I_ESP_D2x_S;
    Double I_ESP_Dxy_Py = I_ESP_Fx2y_S+ABY*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Py = I_ESP_Fxyz_S+ABY*I_ESP_Dxz_S;
    Double I_ESP_D2y_Py = I_ESP_F3y_S+ABY*I_ESP_D2y_S;
    Double I_ESP_Dyz_Py = I_ESP_F2yz_S+ABY*I_ESP_Dyz_S;
    Double I_ESP_D2z_Py = I_ESP_Fy2z_S+ABY*I_ESP_D2z_S;
    Double I_ESP_D2x_Pz = I_ESP_F2xz_S+ABZ*I_ESP_D2x_S;
    Double I_ESP_Dxy_Pz = I_ESP_Fxyz_S+ABZ*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Pz = I_ESP_Fx2z_S+ABZ*I_ESP_Dxz_S;
    Double I_ESP_D2y_Pz = I_ESP_F2yz_S+ABZ*I_ESP_D2y_S;
    Double I_ESP_Dyz_Pz = I_ESP_Fy2z_S+ABZ*I_ESP_Dyz_S;
    Double I_ESP_D2z_Pz = I_ESP_F3z_S+ABZ*I_ESP_D2z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 5 integrals are omitted 
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
    Double I_ESP_F2xy_Py = I_ESP_G2x2y_S+ABY*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Py = I_ESP_G2xyz_S+ABY*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Py = I_ESP_Gx3y_S+ABY*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Py = I_ESP_Gx2yz_S+ABY*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Py = I_ESP_Gxy2z_S+ABY*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Py = I_ESP_G4y_S+ABY*I_ESP_F3y_S;
    Double I_ESP_F2yz_Py = I_ESP_G3yz_S+ABY*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Py = I_ESP_G2y2z_S+ABY*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Py = I_ESP_Gy3z_S+ABY*I_ESP_F3z_S;
    Double I_ESP_F2xz_Pz = I_ESP_G2x2z_S+ABZ*I_ESP_F2xz_S;
    Double I_ESP_Fxyz_Pz = I_ESP_Gxy2z_S+ABZ*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Pz = I_ESP_Gx3z_S+ABZ*I_ESP_Fx2z_S;
    Double I_ESP_F2yz_Pz = I_ESP_G2y2z_S+ABZ*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Pz = I_ESP_Gy3z_S+ABZ*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Pz = I_ESP_G4z_S+ABZ*I_ESP_F3z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_D_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_P
     * RHS shell quartet name: SQ_ESP_D_P
     ************************************************************/
    Double I_ESP_D2x_D2x = I_ESP_F3x_Px+ABX*I_ESP_D2x_Px;
    Double I_ESP_Dxy_D2x = I_ESP_F2xy_Px+ABX*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_D2x = I_ESP_F2xz_Px+ABX*I_ESP_Dxz_Px;
    Double I_ESP_D2y_D2x = I_ESP_Fx2y_Px+ABX*I_ESP_D2y_Px;
    Double I_ESP_Dyz_D2x = I_ESP_Fxyz_Px+ABX*I_ESP_Dyz_Px;
    Double I_ESP_D2z_D2x = I_ESP_Fx2z_Px+ABX*I_ESP_D2z_Px;
    Double I_ESP_D2x_Dxy = I_ESP_F2xy_Px+ABY*I_ESP_D2x_Px;
    Double I_ESP_Dxy_Dxy = I_ESP_Fx2y_Px+ABY*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_Dxy = I_ESP_Fxyz_Px+ABY*I_ESP_Dxz_Px;
    Double I_ESP_D2y_Dxy = I_ESP_F3y_Px+ABY*I_ESP_D2y_Px;
    Double I_ESP_Dyz_Dxy = I_ESP_F2yz_Px+ABY*I_ESP_Dyz_Px;
    Double I_ESP_D2z_Dxy = I_ESP_Fy2z_Px+ABY*I_ESP_D2z_Px;
    Double I_ESP_D2x_Dxz = I_ESP_F2xz_Px+ABZ*I_ESP_D2x_Px;
    Double I_ESP_Dxy_Dxz = I_ESP_Fxyz_Px+ABZ*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_Dxz = I_ESP_Fx2z_Px+ABZ*I_ESP_Dxz_Px;
    Double I_ESP_D2y_Dxz = I_ESP_F2yz_Px+ABZ*I_ESP_D2y_Px;
    Double I_ESP_Dyz_Dxz = I_ESP_Fy2z_Px+ABZ*I_ESP_Dyz_Px;
    Double I_ESP_D2z_Dxz = I_ESP_F3z_Px+ABZ*I_ESP_D2z_Px;
    Double I_ESP_D2x_D2y = I_ESP_F2xy_Py+ABY*I_ESP_D2x_Py;
    Double I_ESP_Dxy_D2y = I_ESP_Fx2y_Py+ABY*I_ESP_Dxy_Py;
    Double I_ESP_Dxz_D2y = I_ESP_Fxyz_Py+ABY*I_ESP_Dxz_Py;
    Double I_ESP_D2y_D2y = I_ESP_F3y_Py+ABY*I_ESP_D2y_Py;
    Double I_ESP_Dyz_D2y = I_ESP_F2yz_Py+ABY*I_ESP_Dyz_Py;
    Double I_ESP_D2z_D2y = I_ESP_Fy2z_Py+ABY*I_ESP_D2z_Py;
    Double I_ESP_D2x_Dyz = I_ESP_F2xz_Py+ABZ*I_ESP_D2x_Py;
    Double I_ESP_Dxy_Dyz = I_ESP_Fxyz_Py+ABZ*I_ESP_Dxy_Py;
    Double I_ESP_Dxz_Dyz = I_ESP_Fx2z_Py+ABZ*I_ESP_Dxz_Py;
    Double I_ESP_D2y_Dyz = I_ESP_F2yz_Py+ABZ*I_ESP_D2y_Py;
    Double I_ESP_Dyz_Dyz = I_ESP_Fy2z_Py+ABZ*I_ESP_Dyz_Py;
    Double I_ESP_D2z_Dyz = I_ESP_F3z_Py+ABZ*I_ESP_D2z_Py;
    Double I_ESP_D2x_D2z = I_ESP_F2xz_Pz+ABZ*I_ESP_D2x_Pz;
    Double I_ESP_Dxy_D2z = I_ESP_Fxyz_Pz+ABZ*I_ESP_Dxy_Pz;
    Double I_ESP_Dxz_D2z = I_ESP_Fx2z_Pz+ABZ*I_ESP_Dxz_Pz;
    Double I_ESP_D2y_D2z = I_ESP_F2yz_Pz+ABZ*I_ESP_D2y_Pz;
    Double I_ESP_Dyz_D2z = I_ESP_Fy2z_Pz+ABZ*I_ESP_Dyz_Pz;
    Double I_ESP_D2z_D2z = I_ESP_F3z_Pz+ABZ*I_ESP_D2z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_a
     * RHS shell quartet name: SQ_ESP_G_S_a
     ************************************************************/
    Double I_ESP_G4x_Px_a = I_ESP_H5x_S_a+ABX*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Px_a = I_ESP_H4xy_S_a+ABX*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Px_a = I_ESP_H4xz_S_a+ABX*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Px_a = I_ESP_H3x2y_S_a+ABX*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Px_a = I_ESP_H3xyz_S_a+ABX*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Px_a = I_ESP_H3x2z_S_a+ABX*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Px_a = I_ESP_H2x3y_S_a+ABX*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Px_a = I_ESP_H2x2yz_S_a+ABX*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Px_a = I_ESP_H2xy2z_S_a+ABX*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Px_a = I_ESP_H2x3z_S_a+ABX*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Px_a = I_ESP_Hx4y_S_a+ABX*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Px_a = I_ESP_Hx3yz_S_a+ABX*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Px_a = I_ESP_Hx2y2z_S_a+ABX*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Px_a = I_ESP_Hxy3z_S_a+ABX*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Px_a = I_ESP_Hx4z_S_a+ABX*I_ESP_G4z_S_a;
    Double I_ESP_G4x_Py_a = I_ESP_H4xy_S_a+ABY*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Py_a = I_ESP_H3x2y_S_a+ABY*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Py_a = I_ESP_H3xyz_S_a+ABY*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Py_a = I_ESP_H2x3y_S_a+ABY*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Py_a = I_ESP_H2x2yz_S_a+ABY*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Py_a = I_ESP_H2xy2z_S_a+ABY*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Py_a = I_ESP_Hx4y_S_a+ABY*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Py_a = I_ESP_Hx3yz_S_a+ABY*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Py_a = I_ESP_Hx2y2z_S_a+ABY*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Py_a = I_ESP_Hxy3z_S_a+ABY*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Py_a = I_ESP_H5y_S_a+ABY*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Py_a = I_ESP_H4yz_S_a+ABY*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Py_a = I_ESP_H3y2z_S_a+ABY*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Py_a = I_ESP_H2y3z_S_a+ABY*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Py_a = I_ESP_Hy4z_S_a+ABY*I_ESP_G4z_S_a;
    Double I_ESP_G4x_Pz_a = I_ESP_H4xz_S_a+ABZ*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Pz_a = I_ESP_H3xyz_S_a+ABZ*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Pz_a = I_ESP_H3x2z_S_a+ABZ*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Pz_a = I_ESP_H2x2yz_S_a+ABZ*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Pz_a = I_ESP_H2xy2z_S_a+ABZ*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Pz_a = I_ESP_H2x3z_S_a+ABZ*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Pz_a = I_ESP_Hx3yz_S_a+ABZ*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Pz_a = I_ESP_Hx2y2z_S_a+ABZ*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Pz_a = I_ESP_Hxy3z_S_a+ABZ*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Pz_a = I_ESP_Hx4z_S_a+ABZ*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Pz_a = I_ESP_H4yz_S_a+ABZ*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Pz_a = I_ESP_H3y2z_S_a+ABZ*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Pz_a = I_ESP_H2y3z_S_a+ABZ*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Pz_a = I_ESP_Hy4z_S_a+ABZ*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Pz_a = I_ESP_H5z_S_a+ABZ*I_ESP_G4z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 7 integrals are omitted 
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
    Double I_ESP_H4xz_Pz_a = I_ESP_I4x2z_S_a+ABZ*I_ESP_H4xz_S_a;
    Double I_ESP_H3xyz_Pz_a = I_ESP_I3xy2z_S_a+ABZ*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Pz_a = I_ESP_I3x3z_S_a+ABZ*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x2yz_Pz_a = I_ESP_I2x2y2z_S_a+ABZ*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Pz_a = I_ESP_I2xy3z_S_a+ABZ*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Pz_a = I_ESP_I2x4z_S_a+ABZ*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx3yz_Pz_a = I_ESP_Ix3y2z_S_a+ABZ*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Pz_a = I_ESP_Ix2y3z_S_a+ABZ*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Pz_a = I_ESP_Ixy4z_S_a+ABZ*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Pz_a = I_ESP_Ix5z_S_a+ABZ*I_ESP_Hx4z_S_a;
    Double I_ESP_H4yz_Pz_a = I_ESP_I4y2z_S_a+ABZ*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Pz_a = I_ESP_I3y3z_S_a+ABZ*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Pz_a = I_ESP_I2y4z_S_a+ABZ*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Pz_a = I_ESP_Iy5z_S_a+ABZ*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Pz_a = I_ESP_I6z_S_a+ABZ*I_ESP_H5z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_G_P_a
     ************************************************************/
    Double I_ESP_G4x_D2x_a = I_ESP_H5x_Px_a+ABX*I_ESP_G4x_Px_a;
    Double I_ESP_G3xy_D2x_a = I_ESP_H4xy_Px_a+ABX*I_ESP_G3xy_Px_a;
    Double I_ESP_G3xz_D2x_a = I_ESP_H4xz_Px_a+ABX*I_ESP_G3xz_Px_a;
    Double I_ESP_G2x2y_D2x_a = I_ESP_H3x2y_Px_a+ABX*I_ESP_G2x2y_Px_a;
    Double I_ESP_G2xyz_D2x_a = I_ESP_H3xyz_Px_a+ABX*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_D2x_a = I_ESP_H3x2z_Px_a+ABX*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx3y_D2x_a = I_ESP_H2x3y_Px_a+ABX*I_ESP_Gx3y_Px_a;
    Double I_ESP_Gx2yz_D2x_a = I_ESP_H2x2yz_Px_a+ABX*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_D2x_a = I_ESP_H2xy2z_Px_a+ABX*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_D2x_a = I_ESP_H2x3z_Px_a+ABX*I_ESP_Gx3z_Px_a;
    Double I_ESP_G4y_D2x_a = I_ESP_Hx4y_Px_a+ABX*I_ESP_G4y_Px_a;
    Double I_ESP_G3yz_D2x_a = I_ESP_Hx3yz_Px_a+ABX*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_D2x_a = I_ESP_Hx2y2z_Px_a+ABX*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_D2x_a = I_ESP_Hxy3z_Px_a+ABX*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_D2x_a = I_ESP_Hx4z_Px_a+ABX*I_ESP_G4z_Px_a;
    Double I_ESP_G4x_Dxy_a = I_ESP_H4xy_Px_a+ABY*I_ESP_G4x_Px_a;
    Double I_ESP_G3xy_Dxy_a = I_ESP_H3x2y_Px_a+ABY*I_ESP_G3xy_Px_a;
    Double I_ESP_G3xz_Dxy_a = I_ESP_H3xyz_Px_a+ABY*I_ESP_G3xz_Px_a;
    Double I_ESP_G2x2y_Dxy_a = I_ESP_H2x3y_Px_a+ABY*I_ESP_G2x2y_Px_a;
    Double I_ESP_G2xyz_Dxy_a = I_ESP_H2x2yz_Px_a+ABY*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_Dxy_a = I_ESP_H2xy2z_Px_a+ABY*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx3y_Dxy_a = I_ESP_Hx4y_Px_a+ABY*I_ESP_Gx3y_Px_a;
    Double I_ESP_Gx2yz_Dxy_a = I_ESP_Hx3yz_Px_a+ABY*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_Dxy_a = I_ESP_Hx2y2z_Px_a+ABY*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_Dxy_a = I_ESP_Hxy3z_Px_a+ABY*I_ESP_Gx3z_Px_a;
    Double I_ESP_G4y_Dxy_a = I_ESP_H5y_Px_a+ABY*I_ESP_G4y_Px_a;
    Double I_ESP_G3yz_Dxy_a = I_ESP_H4yz_Px_a+ABY*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_Dxy_a = I_ESP_H3y2z_Px_a+ABY*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_Dxy_a = I_ESP_H2y3z_Px_a+ABY*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_Dxy_a = I_ESP_Hy4z_Px_a+ABY*I_ESP_G4z_Px_a;
    Double I_ESP_G4x_Dxz_a = I_ESP_H4xz_Px_a+ABZ*I_ESP_G4x_Px_a;
    Double I_ESP_G3xy_Dxz_a = I_ESP_H3xyz_Px_a+ABZ*I_ESP_G3xy_Px_a;
    Double I_ESP_G3xz_Dxz_a = I_ESP_H3x2z_Px_a+ABZ*I_ESP_G3xz_Px_a;
    Double I_ESP_G2x2y_Dxz_a = I_ESP_H2x2yz_Px_a+ABZ*I_ESP_G2x2y_Px_a;
    Double I_ESP_G2xyz_Dxz_a = I_ESP_H2xy2z_Px_a+ABZ*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_Dxz_a = I_ESP_H2x3z_Px_a+ABZ*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx3y_Dxz_a = I_ESP_Hx3yz_Px_a+ABZ*I_ESP_Gx3y_Px_a;
    Double I_ESP_Gx2yz_Dxz_a = I_ESP_Hx2y2z_Px_a+ABZ*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_Dxz_a = I_ESP_Hxy3z_Px_a+ABZ*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_Dxz_a = I_ESP_Hx4z_Px_a+ABZ*I_ESP_Gx3z_Px_a;
    Double I_ESP_G4y_Dxz_a = I_ESP_H4yz_Px_a+ABZ*I_ESP_G4y_Px_a;
    Double I_ESP_G3yz_Dxz_a = I_ESP_H3y2z_Px_a+ABZ*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_Dxz_a = I_ESP_H2y3z_Px_a+ABZ*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_Dxz_a = I_ESP_Hy4z_Px_a+ABZ*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_Dxz_a = I_ESP_H5z_Px_a+ABZ*I_ESP_G4z_Px_a;
    Double I_ESP_G4x_D2y_a = I_ESP_H4xy_Py_a+ABY*I_ESP_G4x_Py_a;
    Double I_ESP_G3xy_D2y_a = I_ESP_H3x2y_Py_a+ABY*I_ESP_G3xy_Py_a;
    Double I_ESP_G3xz_D2y_a = I_ESP_H3xyz_Py_a+ABY*I_ESP_G3xz_Py_a;
    Double I_ESP_G2x2y_D2y_a = I_ESP_H2x3y_Py_a+ABY*I_ESP_G2x2y_Py_a;
    Double I_ESP_G2xyz_D2y_a = I_ESP_H2x2yz_Py_a+ABY*I_ESP_G2xyz_Py_a;
    Double I_ESP_G2x2z_D2y_a = I_ESP_H2xy2z_Py_a+ABY*I_ESP_G2x2z_Py_a;
    Double I_ESP_Gx3y_D2y_a = I_ESP_Hx4y_Py_a+ABY*I_ESP_Gx3y_Py_a;
    Double I_ESP_Gx2yz_D2y_a = I_ESP_Hx3yz_Py_a+ABY*I_ESP_Gx2yz_Py_a;
    Double I_ESP_Gxy2z_D2y_a = I_ESP_Hx2y2z_Py_a+ABY*I_ESP_Gxy2z_Py_a;
    Double I_ESP_Gx3z_D2y_a = I_ESP_Hxy3z_Py_a+ABY*I_ESP_Gx3z_Py_a;
    Double I_ESP_G4y_D2y_a = I_ESP_H5y_Py_a+ABY*I_ESP_G4y_Py_a;
    Double I_ESP_G3yz_D2y_a = I_ESP_H4yz_Py_a+ABY*I_ESP_G3yz_Py_a;
    Double I_ESP_G2y2z_D2y_a = I_ESP_H3y2z_Py_a+ABY*I_ESP_G2y2z_Py_a;
    Double I_ESP_Gy3z_D2y_a = I_ESP_H2y3z_Py_a+ABY*I_ESP_Gy3z_Py_a;
    Double I_ESP_G4z_D2y_a = I_ESP_Hy4z_Py_a+ABY*I_ESP_G4z_Py_a;
    Double I_ESP_G4x_Dyz_a = I_ESP_H4xz_Py_a+ABZ*I_ESP_G4x_Py_a;
    Double I_ESP_G3xy_Dyz_a = I_ESP_H3xyz_Py_a+ABZ*I_ESP_G3xy_Py_a;
    Double I_ESP_G3xz_Dyz_a = I_ESP_H3x2z_Py_a+ABZ*I_ESP_G3xz_Py_a;
    Double I_ESP_G2x2y_Dyz_a = I_ESP_H2x2yz_Py_a+ABZ*I_ESP_G2x2y_Py_a;
    Double I_ESP_G2xyz_Dyz_a = I_ESP_H2xy2z_Py_a+ABZ*I_ESP_G2xyz_Py_a;
    Double I_ESP_G2x2z_Dyz_a = I_ESP_H2x3z_Py_a+ABZ*I_ESP_G2x2z_Py_a;
    Double I_ESP_Gx3y_Dyz_a = I_ESP_Hx3yz_Py_a+ABZ*I_ESP_Gx3y_Py_a;
    Double I_ESP_Gx2yz_Dyz_a = I_ESP_Hx2y2z_Py_a+ABZ*I_ESP_Gx2yz_Py_a;
    Double I_ESP_Gxy2z_Dyz_a = I_ESP_Hxy3z_Py_a+ABZ*I_ESP_Gxy2z_Py_a;
    Double I_ESP_Gx3z_Dyz_a = I_ESP_Hx4z_Py_a+ABZ*I_ESP_Gx3z_Py_a;
    Double I_ESP_G4y_Dyz_a = I_ESP_H4yz_Py_a+ABZ*I_ESP_G4y_Py_a;
    Double I_ESP_G3yz_Dyz_a = I_ESP_H3y2z_Py_a+ABZ*I_ESP_G3yz_Py_a;
    Double I_ESP_G2y2z_Dyz_a = I_ESP_H2y3z_Py_a+ABZ*I_ESP_G2y2z_Py_a;
    Double I_ESP_Gy3z_Dyz_a = I_ESP_Hy4z_Py_a+ABZ*I_ESP_Gy3z_Py_a;
    Double I_ESP_G4z_Dyz_a = I_ESP_H5z_Py_a+ABZ*I_ESP_G4z_Py_a;
    Double I_ESP_G4x_D2z_a = I_ESP_H4xz_Pz_a+ABZ*I_ESP_G4x_Pz_a;
    Double I_ESP_G3xy_D2z_a = I_ESP_H3xyz_Pz_a+ABZ*I_ESP_G3xy_Pz_a;
    Double I_ESP_G3xz_D2z_a = I_ESP_H3x2z_Pz_a+ABZ*I_ESP_G3xz_Pz_a;
    Double I_ESP_G2x2y_D2z_a = I_ESP_H2x2yz_Pz_a+ABZ*I_ESP_G2x2y_Pz_a;
    Double I_ESP_G2xyz_D2z_a = I_ESP_H2xy2z_Pz_a+ABZ*I_ESP_G2xyz_Pz_a;
    Double I_ESP_G2x2z_D2z_a = I_ESP_H2x3z_Pz_a+ABZ*I_ESP_G2x2z_Pz_a;
    Double I_ESP_Gx3y_D2z_a = I_ESP_Hx3yz_Pz_a+ABZ*I_ESP_Gx3y_Pz_a;
    Double I_ESP_Gx2yz_D2z_a = I_ESP_Hx2y2z_Pz_a+ABZ*I_ESP_Gx2yz_Pz_a;
    Double I_ESP_Gxy2z_D2z_a = I_ESP_Hxy3z_Pz_a+ABZ*I_ESP_Gxy2z_Pz_a;
    Double I_ESP_Gx3z_D2z_a = I_ESP_Hx4z_Pz_a+ABZ*I_ESP_Gx3z_Pz_a;
    Double I_ESP_G4y_D2z_a = I_ESP_H4yz_Pz_a+ABZ*I_ESP_G4y_Pz_a;
    Double I_ESP_G3yz_D2z_a = I_ESP_H3y2z_Pz_a+ABZ*I_ESP_G3yz_Pz_a;
    Double I_ESP_G2y2z_D2z_a = I_ESP_H2y3z_Pz_a+ABZ*I_ESP_G2y2z_Pz_a;
    Double I_ESP_Gy3z_D2z_a = I_ESP_Hy4z_Pz_a+ABZ*I_ESP_Gy3z_Pz_a;
    Double I_ESP_G4z_D2z_a = I_ESP_H5z_Pz_a+ABZ*I_ESP_G4z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S_aa
     * RHS shell quartet name: SQ_ESP_I_S_aa
     ************************************************************/
    Double I_ESP_I6x_Px_aa = I_ESP_K7x_S_aa+ABX*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Px_aa = I_ESP_K6xy_S_aa+ABX*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Px_aa = I_ESP_K6xz_S_aa+ABX*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Px_aa = I_ESP_K5x2y_S_aa+ABX*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Px_aa = I_ESP_K5xyz_S_aa+ABX*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Px_aa = I_ESP_K5x2z_S_aa+ABX*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Px_aa = I_ESP_K4x3y_S_aa+ABX*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Px_aa = I_ESP_K4x2yz_S_aa+ABX*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Px_aa = I_ESP_K4xy2z_S_aa+ABX*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Px_aa = I_ESP_K4x3z_S_aa+ABX*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Px_aa = I_ESP_K3x4y_S_aa+ABX*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Px_aa = I_ESP_K3x3yz_S_aa+ABX*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Px_aa = I_ESP_K3x2y2z_S_aa+ABX*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Px_aa = I_ESP_K3xy3z_S_aa+ABX*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Px_aa = I_ESP_K3x4z_S_aa+ABX*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Px_aa = I_ESP_K2x5y_S_aa+ABX*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Px_aa = I_ESP_K2x4yz_S_aa+ABX*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Px_aa = I_ESP_K2x3y2z_S_aa+ABX*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Px_aa = I_ESP_K2x2y3z_S_aa+ABX*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Px_aa = I_ESP_K2xy4z_S_aa+ABX*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Px_aa = I_ESP_K2x5z_S_aa+ABX*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Px_aa = I_ESP_Kx6y_S_aa+ABX*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Px_aa = I_ESP_Kx5yz_S_aa+ABX*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Px_aa = I_ESP_Kx4y2z_S_aa+ABX*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Px_aa = I_ESP_Kx3y3z_S_aa+ABX*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Px_aa = I_ESP_Kx2y4z_S_aa+ABX*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Px_aa = I_ESP_Kxy5z_S_aa+ABX*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Px_aa = I_ESP_Kx6z_S_aa+ABX*I_ESP_I6z_S_aa;
    Double I_ESP_I6x_Py_aa = I_ESP_K6xy_S_aa+ABY*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Py_aa = I_ESP_K5x2y_S_aa+ABY*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Py_aa = I_ESP_K5xyz_S_aa+ABY*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Py_aa = I_ESP_K4x3y_S_aa+ABY*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Py_aa = I_ESP_K4x2yz_S_aa+ABY*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Py_aa = I_ESP_K4xy2z_S_aa+ABY*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Py_aa = I_ESP_K3x4y_S_aa+ABY*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Py_aa = I_ESP_K3x3yz_S_aa+ABY*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Py_aa = I_ESP_K3x2y2z_S_aa+ABY*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Py_aa = I_ESP_K3xy3z_S_aa+ABY*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Py_aa = I_ESP_K2x5y_S_aa+ABY*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Py_aa = I_ESP_K2x4yz_S_aa+ABY*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Py_aa = I_ESP_K2x3y2z_S_aa+ABY*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Py_aa = I_ESP_K2x2y3z_S_aa+ABY*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Py_aa = I_ESP_K2xy4z_S_aa+ABY*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Py_aa = I_ESP_Kx6y_S_aa+ABY*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Py_aa = I_ESP_Kx5yz_S_aa+ABY*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Py_aa = I_ESP_Kx4y2z_S_aa+ABY*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Py_aa = I_ESP_Kx3y3z_S_aa+ABY*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Py_aa = I_ESP_Kx2y4z_S_aa+ABY*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Py_aa = I_ESP_Kxy5z_S_aa+ABY*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Py_aa = I_ESP_K7y_S_aa+ABY*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Py_aa = I_ESP_K6yz_S_aa+ABY*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Py_aa = I_ESP_K5y2z_S_aa+ABY*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Py_aa = I_ESP_K4y3z_S_aa+ABY*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Py_aa = I_ESP_K3y4z_S_aa+ABY*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Py_aa = I_ESP_K2y5z_S_aa+ABY*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Py_aa = I_ESP_Ky6z_S_aa+ABY*I_ESP_I6z_S_aa;
    Double I_ESP_I6x_Pz_aa = I_ESP_K6xz_S_aa+ABZ*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Pz_aa = I_ESP_K5xyz_S_aa+ABZ*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Pz_aa = I_ESP_K5x2z_S_aa+ABZ*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Pz_aa = I_ESP_K4x2yz_S_aa+ABZ*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Pz_aa = I_ESP_K4xy2z_S_aa+ABZ*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Pz_aa = I_ESP_K4x3z_S_aa+ABZ*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Pz_aa = I_ESP_K3x3yz_S_aa+ABZ*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Pz_aa = I_ESP_K3x2y2z_S_aa+ABZ*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Pz_aa = I_ESP_K3xy3z_S_aa+ABZ*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Pz_aa = I_ESP_K3x4z_S_aa+ABZ*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Pz_aa = I_ESP_K2x4yz_S_aa+ABZ*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Pz_aa = I_ESP_K2x3y2z_S_aa+ABZ*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Pz_aa = I_ESP_K2x2y3z_S_aa+ABZ*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Pz_aa = I_ESP_K2xy4z_S_aa+ABZ*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Pz_aa = I_ESP_K2x5z_S_aa+ABZ*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Pz_aa = I_ESP_Kx5yz_S_aa+ABZ*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Pz_aa = I_ESP_Kx4y2z_S_aa+ABZ*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Pz_aa = I_ESP_Kx3y3z_S_aa+ABZ*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Pz_aa = I_ESP_Kx2y4z_S_aa+ABZ*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Pz_aa = I_ESP_Kxy5z_S_aa+ABZ*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Pz_aa = I_ESP_Kx6z_S_aa+ABZ*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Pz_aa = I_ESP_K6yz_S_aa+ABZ*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Pz_aa = I_ESP_K5y2z_S_aa+ABZ*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Pz_aa = I_ESP_K4y3z_S_aa+ABZ*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Pz_aa = I_ESP_K3y4z_S_aa+ABZ*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Pz_aa = I_ESP_K2y5z_S_aa+ABZ*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Pz_aa = I_ESP_Ky6z_S_aa+ABZ*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Pz_aa = I_ESP_K7z_S_aa+ABZ*I_ESP_I6z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 9 integrals are omitted 
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
    Double I_ESP_K6xz_Pz_aa = I_ESP_L6x2z_S_aa+ABZ*I_ESP_K6xz_S_aa;
    Double I_ESP_K5xyz_Pz_aa = I_ESP_L5xy2z_S_aa+ABZ*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Pz_aa = I_ESP_L5x3z_S_aa+ABZ*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x2yz_Pz_aa = I_ESP_L4x2y2z_S_aa+ABZ*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Pz_aa = I_ESP_L4xy3z_S_aa+ABZ*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Pz_aa = I_ESP_L4x4z_S_aa+ABZ*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x3yz_Pz_aa = I_ESP_L3x3y2z_S_aa+ABZ*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Pz_aa = I_ESP_L3x2y3z_S_aa+ABZ*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Pz_aa = I_ESP_L3xy4z_S_aa+ABZ*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Pz_aa = I_ESP_L3x5z_S_aa+ABZ*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x4yz_Pz_aa = I_ESP_L2x4y2z_S_aa+ABZ*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Pz_aa = I_ESP_L2x3y3z_S_aa+ABZ*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Pz_aa = I_ESP_L2x2y4z_S_aa+ABZ*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Pz_aa = I_ESP_L2xy5z_S_aa+ABZ*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Pz_aa = I_ESP_L2x6z_S_aa+ABZ*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx5yz_Pz_aa = I_ESP_Lx5y2z_S_aa+ABZ*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Pz_aa = I_ESP_Lx4y3z_S_aa+ABZ*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Pz_aa = I_ESP_Lx3y4z_S_aa+ABZ*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Pz_aa = I_ESP_Lx2y5z_S_aa+ABZ*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Pz_aa = I_ESP_Lxy6z_S_aa+ABZ*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Pz_aa = I_ESP_Lx7z_S_aa+ABZ*I_ESP_Kx6z_S_aa;
    Double I_ESP_K6yz_Pz_aa = I_ESP_L6y2z_S_aa+ABZ*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Pz_aa = I_ESP_L5y3z_S_aa+ABZ*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Pz_aa = I_ESP_L4y4z_S_aa+ABZ*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Pz_aa = I_ESP_L3y5z_S_aa+ABZ*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Pz_aa = I_ESP_L2y6z_S_aa+ABZ*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Pz_aa = I_ESP_Ly7z_S_aa+ABZ*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Pz_aa = I_ESP_L8z_S_aa+ABZ*I_ESP_K7z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_I_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_I_P_aa
     ************************************************************/
    Double I_ESP_I6x_D2x_aa = I_ESP_K7x_Px_aa+ABX*I_ESP_I6x_Px_aa;
    Double I_ESP_I5xy_D2x_aa = I_ESP_K6xy_Px_aa+ABX*I_ESP_I5xy_Px_aa;
    Double I_ESP_I5xz_D2x_aa = I_ESP_K6xz_Px_aa+ABX*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4x2y_D2x_aa = I_ESP_K5x2y_Px_aa+ABX*I_ESP_I4x2y_Px_aa;
    Double I_ESP_I4xyz_D2x_aa = I_ESP_K5xyz_Px_aa+ABX*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_D2x_aa = I_ESP_K5x2z_Px_aa+ABX*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x3y_D2x_aa = I_ESP_K4x3y_Px_aa+ABX*I_ESP_I3x3y_Px_aa;
    Double I_ESP_I3x2yz_D2x_aa = I_ESP_K4x2yz_Px_aa+ABX*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_D2x_aa = I_ESP_K4xy2z_Px_aa+ABX*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_D2x_aa = I_ESP_K4x3z_Px_aa+ABX*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x4y_D2x_aa = I_ESP_K3x4y_Px_aa+ABX*I_ESP_I2x4y_Px_aa;
    Double I_ESP_I2x3yz_D2x_aa = I_ESP_K3x3yz_Px_aa+ABX*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_D2x_aa = I_ESP_K3x2y2z_Px_aa+ABX*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_D2x_aa = I_ESP_K3xy3z_Px_aa+ABX*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_D2x_aa = I_ESP_K3x4z_Px_aa+ABX*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix5y_D2x_aa = I_ESP_K2x5y_Px_aa+ABX*I_ESP_Ix5y_Px_aa;
    Double I_ESP_Ix4yz_D2x_aa = I_ESP_K2x4yz_Px_aa+ABX*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_D2x_aa = I_ESP_K2x3y2z_Px_aa+ABX*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_D2x_aa = I_ESP_K2x2y3z_Px_aa+ABX*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_D2x_aa = I_ESP_K2xy4z_Px_aa+ABX*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_D2x_aa = I_ESP_K2x5z_Px_aa+ABX*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I6y_D2x_aa = I_ESP_Kx6y_Px_aa+ABX*I_ESP_I6y_Px_aa;
    Double I_ESP_I5yz_D2x_aa = I_ESP_Kx5yz_Px_aa+ABX*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_D2x_aa = I_ESP_Kx4y2z_Px_aa+ABX*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_D2x_aa = I_ESP_Kx3y3z_Px_aa+ABX*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_D2x_aa = I_ESP_Kx2y4z_Px_aa+ABX*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_D2x_aa = I_ESP_Kxy5z_Px_aa+ABX*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_D2x_aa = I_ESP_Kx6z_Px_aa+ABX*I_ESP_I6z_Px_aa;
    Double I_ESP_I6x_Dxy_aa = I_ESP_K6xy_Px_aa+ABY*I_ESP_I6x_Px_aa;
    Double I_ESP_I5xy_Dxy_aa = I_ESP_K5x2y_Px_aa+ABY*I_ESP_I5xy_Px_aa;
    Double I_ESP_I5xz_Dxy_aa = I_ESP_K5xyz_Px_aa+ABY*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4x2y_Dxy_aa = I_ESP_K4x3y_Px_aa+ABY*I_ESP_I4x2y_Px_aa;
    Double I_ESP_I4xyz_Dxy_aa = I_ESP_K4x2yz_Px_aa+ABY*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_Dxy_aa = I_ESP_K4xy2z_Px_aa+ABY*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x3y_Dxy_aa = I_ESP_K3x4y_Px_aa+ABY*I_ESP_I3x3y_Px_aa;
    Double I_ESP_I3x2yz_Dxy_aa = I_ESP_K3x3yz_Px_aa+ABY*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_Dxy_aa = I_ESP_K3x2y2z_Px_aa+ABY*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_Dxy_aa = I_ESP_K3xy3z_Px_aa+ABY*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x4y_Dxy_aa = I_ESP_K2x5y_Px_aa+ABY*I_ESP_I2x4y_Px_aa;
    Double I_ESP_I2x3yz_Dxy_aa = I_ESP_K2x4yz_Px_aa+ABY*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_Dxy_aa = I_ESP_K2x3y2z_Px_aa+ABY*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_Dxy_aa = I_ESP_K2x2y3z_Px_aa+ABY*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_Dxy_aa = I_ESP_K2xy4z_Px_aa+ABY*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix5y_Dxy_aa = I_ESP_Kx6y_Px_aa+ABY*I_ESP_Ix5y_Px_aa;
    Double I_ESP_Ix4yz_Dxy_aa = I_ESP_Kx5yz_Px_aa+ABY*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_Dxy_aa = I_ESP_Kx4y2z_Px_aa+ABY*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_Dxy_aa = I_ESP_Kx3y3z_Px_aa+ABY*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_Dxy_aa = I_ESP_Kx2y4z_Px_aa+ABY*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_Dxy_aa = I_ESP_Kxy5z_Px_aa+ABY*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I6y_Dxy_aa = I_ESP_K7y_Px_aa+ABY*I_ESP_I6y_Px_aa;
    Double I_ESP_I5yz_Dxy_aa = I_ESP_K6yz_Px_aa+ABY*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_Dxy_aa = I_ESP_K5y2z_Px_aa+ABY*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_Dxy_aa = I_ESP_K4y3z_Px_aa+ABY*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_Dxy_aa = I_ESP_K3y4z_Px_aa+ABY*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_Dxy_aa = I_ESP_K2y5z_Px_aa+ABY*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_Dxy_aa = I_ESP_Ky6z_Px_aa+ABY*I_ESP_I6z_Px_aa;
    Double I_ESP_I6x_Dxz_aa = I_ESP_K6xz_Px_aa+ABZ*I_ESP_I6x_Px_aa;
    Double I_ESP_I5xy_Dxz_aa = I_ESP_K5xyz_Px_aa+ABZ*I_ESP_I5xy_Px_aa;
    Double I_ESP_I5xz_Dxz_aa = I_ESP_K5x2z_Px_aa+ABZ*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4x2y_Dxz_aa = I_ESP_K4x2yz_Px_aa+ABZ*I_ESP_I4x2y_Px_aa;
    Double I_ESP_I4xyz_Dxz_aa = I_ESP_K4xy2z_Px_aa+ABZ*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_Dxz_aa = I_ESP_K4x3z_Px_aa+ABZ*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x3y_Dxz_aa = I_ESP_K3x3yz_Px_aa+ABZ*I_ESP_I3x3y_Px_aa;
    Double I_ESP_I3x2yz_Dxz_aa = I_ESP_K3x2y2z_Px_aa+ABZ*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_Dxz_aa = I_ESP_K3xy3z_Px_aa+ABZ*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_Dxz_aa = I_ESP_K3x4z_Px_aa+ABZ*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x4y_Dxz_aa = I_ESP_K2x4yz_Px_aa+ABZ*I_ESP_I2x4y_Px_aa;
    Double I_ESP_I2x3yz_Dxz_aa = I_ESP_K2x3y2z_Px_aa+ABZ*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_Dxz_aa = I_ESP_K2x2y3z_Px_aa+ABZ*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_Dxz_aa = I_ESP_K2xy4z_Px_aa+ABZ*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_Dxz_aa = I_ESP_K2x5z_Px_aa+ABZ*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix5y_Dxz_aa = I_ESP_Kx5yz_Px_aa+ABZ*I_ESP_Ix5y_Px_aa;
    Double I_ESP_Ix4yz_Dxz_aa = I_ESP_Kx4y2z_Px_aa+ABZ*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_Dxz_aa = I_ESP_Kx3y3z_Px_aa+ABZ*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_Dxz_aa = I_ESP_Kx2y4z_Px_aa+ABZ*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_Dxz_aa = I_ESP_Kxy5z_Px_aa+ABZ*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_Dxz_aa = I_ESP_Kx6z_Px_aa+ABZ*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I6y_Dxz_aa = I_ESP_K6yz_Px_aa+ABZ*I_ESP_I6y_Px_aa;
    Double I_ESP_I5yz_Dxz_aa = I_ESP_K5y2z_Px_aa+ABZ*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_Dxz_aa = I_ESP_K4y3z_Px_aa+ABZ*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_Dxz_aa = I_ESP_K3y4z_Px_aa+ABZ*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_Dxz_aa = I_ESP_K2y5z_Px_aa+ABZ*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_Dxz_aa = I_ESP_Ky6z_Px_aa+ABZ*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_Dxz_aa = I_ESP_K7z_Px_aa+ABZ*I_ESP_I6z_Px_aa;
    Double I_ESP_I6x_D2y_aa = I_ESP_K6xy_Py_aa+ABY*I_ESP_I6x_Py_aa;
    Double I_ESP_I5xy_D2y_aa = I_ESP_K5x2y_Py_aa+ABY*I_ESP_I5xy_Py_aa;
    Double I_ESP_I5xz_D2y_aa = I_ESP_K5xyz_Py_aa+ABY*I_ESP_I5xz_Py_aa;
    Double I_ESP_I4x2y_D2y_aa = I_ESP_K4x3y_Py_aa+ABY*I_ESP_I4x2y_Py_aa;
    Double I_ESP_I4xyz_D2y_aa = I_ESP_K4x2yz_Py_aa+ABY*I_ESP_I4xyz_Py_aa;
    Double I_ESP_I4x2z_D2y_aa = I_ESP_K4xy2z_Py_aa+ABY*I_ESP_I4x2z_Py_aa;
    Double I_ESP_I3x3y_D2y_aa = I_ESP_K3x4y_Py_aa+ABY*I_ESP_I3x3y_Py_aa;
    Double I_ESP_I3x2yz_D2y_aa = I_ESP_K3x3yz_Py_aa+ABY*I_ESP_I3x2yz_Py_aa;
    Double I_ESP_I3xy2z_D2y_aa = I_ESP_K3x2y2z_Py_aa+ABY*I_ESP_I3xy2z_Py_aa;
    Double I_ESP_I3x3z_D2y_aa = I_ESP_K3xy3z_Py_aa+ABY*I_ESP_I3x3z_Py_aa;
    Double I_ESP_I2x4y_D2y_aa = I_ESP_K2x5y_Py_aa+ABY*I_ESP_I2x4y_Py_aa;
    Double I_ESP_I2x3yz_D2y_aa = I_ESP_K2x4yz_Py_aa+ABY*I_ESP_I2x3yz_Py_aa;
    Double I_ESP_I2x2y2z_D2y_aa = I_ESP_K2x3y2z_Py_aa+ABY*I_ESP_I2x2y2z_Py_aa;
    Double I_ESP_I2xy3z_D2y_aa = I_ESP_K2x2y3z_Py_aa+ABY*I_ESP_I2xy3z_Py_aa;
    Double I_ESP_I2x4z_D2y_aa = I_ESP_K2xy4z_Py_aa+ABY*I_ESP_I2x4z_Py_aa;
    Double I_ESP_Ix5y_D2y_aa = I_ESP_Kx6y_Py_aa+ABY*I_ESP_Ix5y_Py_aa;
    Double I_ESP_Ix4yz_D2y_aa = I_ESP_Kx5yz_Py_aa+ABY*I_ESP_Ix4yz_Py_aa;
    Double I_ESP_Ix3y2z_D2y_aa = I_ESP_Kx4y2z_Py_aa+ABY*I_ESP_Ix3y2z_Py_aa;
    Double I_ESP_Ix2y3z_D2y_aa = I_ESP_Kx3y3z_Py_aa+ABY*I_ESP_Ix2y3z_Py_aa;
    Double I_ESP_Ixy4z_D2y_aa = I_ESP_Kx2y4z_Py_aa+ABY*I_ESP_Ixy4z_Py_aa;
    Double I_ESP_Ix5z_D2y_aa = I_ESP_Kxy5z_Py_aa+ABY*I_ESP_Ix5z_Py_aa;
    Double I_ESP_I6y_D2y_aa = I_ESP_K7y_Py_aa+ABY*I_ESP_I6y_Py_aa;
    Double I_ESP_I5yz_D2y_aa = I_ESP_K6yz_Py_aa+ABY*I_ESP_I5yz_Py_aa;
    Double I_ESP_I4y2z_D2y_aa = I_ESP_K5y2z_Py_aa+ABY*I_ESP_I4y2z_Py_aa;
    Double I_ESP_I3y3z_D2y_aa = I_ESP_K4y3z_Py_aa+ABY*I_ESP_I3y3z_Py_aa;
    Double I_ESP_I2y4z_D2y_aa = I_ESP_K3y4z_Py_aa+ABY*I_ESP_I2y4z_Py_aa;
    Double I_ESP_Iy5z_D2y_aa = I_ESP_K2y5z_Py_aa+ABY*I_ESP_Iy5z_Py_aa;
    Double I_ESP_I6z_D2y_aa = I_ESP_Ky6z_Py_aa+ABY*I_ESP_I6z_Py_aa;
    Double I_ESP_I6x_Dyz_aa = I_ESP_K6xz_Py_aa+ABZ*I_ESP_I6x_Py_aa;
    Double I_ESP_I5xy_Dyz_aa = I_ESP_K5xyz_Py_aa+ABZ*I_ESP_I5xy_Py_aa;
    Double I_ESP_I5xz_Dyz_aa = I_ESP_K5x2z_Py_aa+ABZ*I_ESP_I5xz_Py_aa;
    Double I_ESP_I4x2y_Dyz_aa = I_ESP_K4x2yz_Py_aa+ABZ*I_ESP_I4x2y_Py_aa;
    Double I_ESP_I4xyz_Dyz_aa = I_ESP_K4xy2z_Py_aa+ABZ*I_ESP_I4xyz_Py_aa;
    Double I_ESP_I4x2z_Dyz_aa = I_ESP_K4x3z_Py_aa+ABZ*I_ESP_I4x2z_Py_aa;
    Double I_ESP_I3x3y_Dyz_aa = I_ESP_K3x3yz_Py_aa+ABZ*I_ESP_I3x3y_Py_aa;
    Double I_ESP_I3x2yz_Dyz_aa = I_ESP_K3x2y2z_Py_aa+ABZ*I_ESP_I3x2yz_Py_aa;
    Double I_ESP_I3xy2z_Dyz_aa = I_ESP_K3xy3z_Py_aa+ABZ*I_ESP_I3xy2z_Py_aa;
    Double I_ESP_I3x3z_Dyz_aa = I_ESP_K3x4z_Py_aa+ABZ*I_ESP_I3x3z_Py_aa;
    Double I_ESP_I2x4y_Dyz_aa = I_ESP_K2x4yz_Py_aa+ABZ*I_ESP_I2x4y_Py_aa;
    Double I_ESP_I2x3yz_Dyz_aa = I_ESP_K2x3y2z_Py_aa+ABZ*I_ESP_I2x3yz_Py_aa;
    Double I_ESP_I2x2y2z_Dyz_aa = I_ESP_K2x2y3z_Py_aa+ABZ*I_ESP_I2x2y2z_Py_aa;
    Double I_ESP_I2xy3z_Dyz_aa = I_ESP_K2xy4z_Py_aa+ABZ*I_ESP_I2xy3z_Py_aa;
    Double I_ESP_I2x4z_Dyz_aa = I_ESP_K2x5z_Py_aa+ABZ*I_ESP_I2x4z_Py_aa;
    Double I_ESP_Ix5y_Dyz_aa = I_ESP_Kx5yz_Py_aa+ABZ*I_ESP_Ix5y_Py_aa;
    Double I_ESP_Ix4yz_Dyz_aa = I_ESP_Kx4y2z_Py_aa+ABZ*I_ESP_Ix4yz_Py_aa;
    Double I_ESP_Ix3y2z_Dyz_aa = I_ESP_Kx3y3z_Py_aa+ABZ*I_ESP_Ix3y2z_Py_aa;
    Double I_ESP_Ix2y3z_Dyz_aa = I_ESP_Kx2y4z_Py_aa+ABZ*I_ESP_Ix2y3z_Py_aa;
    Double I_ESP_Ixy4z_Dyz_aa = I_ESP_Kxy5z_Py_aa+ABZ*I_ESP_Ixy4z_Py_aa;
    Double I_ESP_Ix5z_Dyz_aa = I_ESP_Kx6z_Py_aa+ABZ*I_ESP_Ix5z_Py_aa;
    Double I_ESP_I6y_Dyz_aa = I_ESP_K6yz_Py_aa+ABZ*I_ESP_I6y_Py_aa;
    Double I_ESP_I5yz_Dyz_aa = I_ESP_K5y2z_Py_aa+ABZ*I_ESP_I5yz_Py_aa;
    Double I_ESP_I4y2z_Dyz_aa = I_ESP_K4y3z_Py_aa+ABZ*I_ESP_I4y2z_Py_aa;
    Double I_ESP_I3y3z_Dyz_aa = I_ESP_K3y4z_Py_aa+ABZ*I_ESP_I3y3z_Py_aa;
    Double I_ESP_I2y4z_Dyz_aa = I_ESP_K2y5z_Py_aa+ABZ*I_ESP_I2y4z_Py_aa;
    Double I_ESP_Iy5z_Dyz_aa = I_ESP_Ky6z_Py_aa+ABZ*I_ESP_Iy5z_Py_aa;
    Double I_ESP_I6z_Dyz_aa = I_ESP_K7z_Py_aa+ABZ*I_ESP_I6z_Py_aa;
    Double I_ESP_I6x_D2z_aa = I_ESP_K6xz_Pz_aa+ABZ*I_ESP_I6x_Pz_aa;
    Double I_ESP_I5xy_D2z_aa = I_ESP_K5xyz_Pz_aa+ABZ*I_ESP_I5xy_Pz_aa;
    Double I_ESP_I5xz_D2z_aa = I_ESP_K5x2z_Pz_aa+ABZ*I_ESP_I5xz_Pz_aa;
    Double I_ESP_I4x2y_D2z_aa = I_ESP_K4x2yz_Pz_aa+ABZ*I_ESP_I4x2y_Pz_aa;
    Double I_ESP_I4xyz_D2z_aa = I_ESP_K4xy2z_Pz_aa+ABZ*I_ESP_I4xyz_Pz_aa;
    Double I_ESP_I4x2z_D2z_aa = I_ESP_K4x3z_Pz_aa+ABZ*I_ESP_I4x2z_Pz_aa;
    Double I_ESP_I3x3y_D2z_aa = I_ESP_K3x3yz_Pz_aa+ABZ*I_ESP_I3x3y_Pz_aa;
    Double I_ESP_I3x2yz_D2z_aa = I_ESP_K3x2y2z_Pz_aa+ABZ*I_ESP_I3x2yz_Pz_aa;
    Double I_ESP_I3xy2z_D2z_aa = I_ESP_K3xy3z_Pz_aa+ABZ*I_ESP_I3xy2z_Pz_aa;
    Double I_ESP_I3x3z_D2z_aa = I_ESP_K3x4z_Pz_aa+ABZ*I_ESP_I3x3z_Pz_aa;
    Double I_ESP_I2x4y_D2z_aa = I_ESP_K2x4yz_Pz_aa+ABZ*I_ESP_I2x4y_Pz_aa;
    Double I_ESP_I2x3yz_D2z_aa = I_ESP_K2x3y2z_Pz_aa+ABZ*I_ESP_I2x3yz_Pz_aa;
    Double I_ESP_I2x2y2z_D2z_aa = I_ESP_K2x2y3z_Pz_aa+ABZ*I_ESP_I2x2y2z_Pz_aa;
    Double I_ESP_I2xy3z_D2z_aa = I_ESP_K2xy4z_Pz_aa+ABZ*I_ESP_I2xy3z_Pz_aa;
    Double I_ESP_I2x4z_D2z_aa = I_ESP_K2x5z_Pz_aa+ABZ*I_ESP_I2x4z_Pz_aa;
    Double I_ESP_Ix5y_D2z_aa = I_ESP_Kx5yz_Pz_aa+ABZ*I_ESP_Ix5y_Pz_aa;
    Double I_ESP_Ix4yz_D2z_aa = I_ESP_Kx4y2z_Pz_aa+ABZ*I_ESP_Ix4yz_Pz_aa;
    Double I_ESP_Ix3y2z_D2z_aa = I_ESP_Kx3y3z_Pz_aa+ABZ*I_ESP_Ix3y2z_Pz_aa;
    Double I_ESP_Ix2y3z_D2z_aa = I_ESP_Kx2y4z_Pz_aa+ABZ*I_ESP_Ix2y3z_Pz_aa;
    Double I_ESP_Ixy4z_D2z_aa = I_ESP_Kxy5z_Pz_aa+ABZ*I_ESP_Ixy4z_Pz_aa;
    Double I_ESP_Ix5z_D2z_aa = I_ESP_Kx6z_Pz_aa+ABZ*I_ESP_Ix5z_Pz_aa;
    Double I_ESP_I6y_D2z_aa = I_ESP_K6yz_Pz_aa+ABZ*I_ESP_I6y_Pz_aa;
    Double I_ESP_I5yz_D2z_aa = I_ESP_K5y2z_Pz_aa+ABZ*I_ESP_I5yz_Pz_aa;
    Double I_ESP_I4y2z_D2z_aa = I_ESP_K4y3z_Pz_aa+ABZ*I_ESP_I4y2z_Pz_aa;
    Double I_ESP_I3y3z_D2z_aa = I_ESP_K3y4z_Pz_aa+ABZ*I_ESP_I3y3z_Pz_aa;
    Double I_ESP_I2y4z_D2z_aa = I_ESP_K2y5z_Pz_aa+ABZ*I_ESP_I2y4z_Pz_aa;
    Double I_ESP_Iy5z_D2z_aa = I_ESP_Ky6z_Pz_aa+ABZ*I_ESP_Iy5z_Pz_aa;
    Double I_ESP_I6z_D2z_aa = I_ESP_K7z_Pz_aa+ABZ*I_ESP_I6z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+0] = 4.0E0*I_ESP_I6x_D2x_aa-2.0E0*4*I_ESP_G4x_D2x_a-2.0E0*5*I_ESP_G4x_D2x_a+4*3*I_ESP_D2x_D2x;
    abcd[iGrid*540+1] = 4.0E0*I_ESP_I5xy_D2x_aa-2.0E0*3*I_ESP_G3xy_D2x_a-2.0E0*4*I_ESP_G3xy_D2x_a+3*2*I_ESP_Dxy_D2x;
    abcd[iGrid*540+2] = 4.0E0*I_ESP_I5xz_D2x_aa-2.0E0*3*I_ESP_G3xz_D2x_a-2.0E0*4*I_ESP_G3xz_D2x_a+3*2*I_ESP_Dxz_D2x;
    abcd[iGrid*540+3] = 4.0E0*I_ESP_I4x2y_D2x_aa-2.0E0*2*I_ESP_G2x2y_D2x_a-2.0E0*3*I_ESP_G2x2y_D2x_a+2*1*I_ESP_D2y_D2x;
    abcd[iGrid*540+4] = 4.0E0*I_ESP_I4xyz_D2x_aa-2.0E0*2*I_ESP_G2xyz_D2x_a-2.0E0*3*I_ESP_G2xyz_D2x_a+2*1*I_ESP_Dyz_D2x;
    abcd[iGrid*540+5] = 4.0E0*I_ESP_I4x2z_D2x_aa-2.0E0*2*I_ESP_G2x2z_D2x_a-2.0E0*3*I_ESP_G2x2z_D2x_a+2*1*I_ESP_D2z_D2x;
    abcd[iGrid*540+6] = 4.0E0*I_ESP_I3x3y_D2x_aa-2.0E0*1*I_ESP_Gx3y_D2x_a-2.0E0*2*I_ESP_Gx3y_D2x_a;
    abcd[iGrid*540+7] = 4.0E0*I_ESP_I3x2yz_D2x_aa-2.0E0*1*I_ESP_Gx2yz_D2x_a-2.0E0*2*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+8] = 4.0E0*I_ESP_I3xy2z_D2x_aa-2.0E0*1*I_ESP_Gxy2z_D2x_a-2.0E0*2*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+9] = 4.0E0*I_ESP_I3x3z_D2x_aa-2.0E0*1*I_ESP_Gx3z_D2x_a-2.0E0*2*I_ESP_Gx3z_D2x_a;
    abcd[iGrid*540+10] = 4.0E0*I_ESP_I2x4y_D2x_aa-2.0E0*1*I_ESP_G4y_D2x_a;
    abcd[iGrid*540+11] = 4.0E0*I_ESP_I2x3yz_D2x_aa-2.0E0*1*I_ESP_G3yz_D2x_a;
    abcd[iGrid*540+12] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2y2z_D2x_a;
    abcd[iGrid*540+13] = 4.0E0*I_ESP_I2xy3z_D2x_aa-2.0E0*1*I_ESP_Gy3z_D2x_a;
    abcd[iGrid*540+14] = 4.0E0*I_ESP_I2x4z_D2x_aa-2.0E0*1*I_ESP_G4z_D2x_a;
    abcd[iGrid*540+15] = 4.0E0*I_ESP_I6x_Dxy_aa-2.0E0*4*I_ESP_G4x_Dxy_a-2.0E0*5*I_ESP_G4x_Dxy_a+4*3*I_ESP_D2x_Dxy;
    abcd[iGrid*540+16] = 4.0E0*I_ESP_I5xy_Dxy_aa-2.0E0*3*I_ESP_G3xy_Dxy_a-2.0E0*4*I_ESP_G3xy_Dxy_a+3*2*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+17] = 4.0E0*I_ESP_I5xz_Dxy_aa-2.0E0*3*I_ESP_G3xz_Dxy_a-2.0E0*4*I_ESP_G3xz_Dxy_a+3*2*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+18] = 4.0E0*I_ESP_I4x2y_Dxy_aa-2.0E0*2*I_ESP_G2x2y_Dxy_a-2.0E0*3*I_ESP_G2x2y_Dxy_a+2*1*I_ESP_D2y_Dxy;
    abcd[iGrid*540+19] = 4.0E0*I_ESP_I4xyz_Dxy_aa-2.0E0*2*I_ESP_G2xyz_Dxy_a-2.0E0*3*I_ESP_G2xyz_Dxy_a+2*1*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+20] = 4.0E0*I_ESP_I4x2z_Dxy_aa-2.0E0*2*I_ESP_G2x2z_Dxy_a-2.0E0*3*I_ESP_G2x2z_Dxy_a+2*1*I_ESP_D2z_Dxy;
    abcd[iGrid*540+21] = 4.0E0*I_ESP_I3x3y_Dxy_aa-2.0E0*1*I_ESP_Gx3y_Dxy_a-2.0E0*2*I_ESP_Gx3y_Dxy_a;
    abcd[iGrid*540+22] = 4.0E0*I_ESP_I3x2yz_Dxy_aa-2.0E0*1*I_ESP_Gx2yz_Dxy_a-2.0E0*2*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+23] = 4.0E0*I_ESP_I3xy2z_Dxy_aa-2.0E0*1*I_ESP_Gxy2z_Dxy_a-2.0E0*2*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+24] = 4.0E0*I_ESP_I3x3z_Dxy_aa-2.0E0*1*I_ESP_Gx3z_Dxy_a-2.0E0*2*I_ESP_Gx3z_Dxy_a;
    abcd[iGrid*540+25] = 4.0E0*I_ESP_I2x4y_Dxy_aa-2.0E0*1*I_ESP_G4y_Dxy_a;
    abcd[iGrid*540+26] = 4.0E0*I_ESP_I2x3yz_Dxy_aa-2.0E0*1*I_ESP_G3yz_Dxy_a;
    abcd[iGrid*540+27] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2y2z_Dxy_a;
    abcd[iGrid*540+28] = 4.0E0*I_ESP_I2xy3z_Dxy_aa-2.0E0*1*I_ESP_Gy3z_Dxy_a;
    abcd[iGrid*540+29] = 4.0E0*I_ESP_I2x4z_Dxy_aa-2.0E0*1*I_ESP_G4z_Dxy_a;
    abcd[iGrid*540+30] = 4.0E0*I_ESP_I6x_Dxz_aa-2.0E0*4*I_ESP_G4x_Dxz_a-2.0E0*5*I_ESP_G4x_Dxz_a+4*3*I_ESP_D2x_Dxz;
    abcd[iGrid*540+31] = 4.0E0*I_ESP_I5xy_Dxz_aa-2.0E0*3*I_ESP_G3xy_Dxz_a-2.0E0*4*I_ESP_G3xy_Dxz_a+3*2*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+32] = 4.0E0*I_ESP_I5xz_Dxz_aa-2.0E0*3*I_ESP_G3xz_Dxz_a-2.0E0*4*I_ESP_G3xz_Dxz_a+3*2*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+33] = 4.0E0*I_ESP_I4x2y_Dxz_aa-2.0E0*2*I_ESP_G2x2y_Dxz_a-2.0E0*3*I_ESP_G2x2y_Dxz_a+2*1*I_ESP_D2y_Dxz;
    abcd[iGrid*540+34] = 4.0E0*I_ESP_I4xyz_Dxz_aa-2.0E0*2*I_ESP_G2xyz_Dxz_a-2.0E0*3*I_ESP_G2xyz_Dxz_a+2*1*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+35] = 4.0E0*I_ESP_I4x2z_Dxz_aa-2.0E0*2*I_ESP_G2x2z_Dxz_a-2.0E0*3*I_ESP_G2x2z_Dxz_a+2*1*I_ESP_D2z_Dxz;
    abcd[iGrid*540+36] = 4.0E0*I_ESP_I3x3y_Dxz_aa-2.0E0*1*I_ESP_Gx3y_Dxz_a-2.0E0*2*I_ESP_Gx3y_Dxz_a;
    abcd[iGrid*540+37] = 4.0E0*I_ESP_I3x2yz_Dxz_aa-2.0E0*1*I_ESP_Gx2yz_Dxz_a-2.0E0*2*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+38] = 4.0E0*I_ESP_I3xy2z_Dxz_aa-2.0E0*1*I_ESP_Gxy2z_Dxz_a-2.0E0*2*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+39] = 4.0E0*I_ESP_I3x3z_Dxz_aa-2.0E0*1*I_ESP_Gx3z_Dxz_a-2.0E0*2*I_ESP_Gx3z_Dxz_a;
    abcd[iGrid*540+40] = 4.0E0*I_ESP_I2x4y_Dxz_aa-2.0E0*1*I_ESP_G4y_Dxz_a;
    abcd[iGrid*540+41] = 4.0E0*I_ESP_I2x3yz_Dxz_aa-2.0E0*1*I_ESP_G3yz_Dxz_a;
    abcd[iGrid*540+42] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2y2z_Dxz_a;
    abcd[iGrid*540+43] = 4.0E0*I_ESP_I2xy3z_Dxz_aa-2.0E0*1*I_ESP_Gy3z_Dxz_a;
    abcd[iGrid*540+44] = 4.0E0*I_ESP_I2x4z_Dxz_aa-2.0E0*1*I_ESP_G4z_Dxz_a;
    abcd[iGrid*540+45] = 4.0E0*I_ESP_I6x_D2y_aa-2.0E0*4*I_ESP_G4x_D2y_a-2.0E0*5*I_ESP_G4x_D2y_a+4*3*I_ESP_D2x_D2y;
    abcd[iGrid*540+46] = 4.0E0*I_ESP_I5xy_D2y_aa-2.0E0*3*I_ESP_G3xy_D2y_a-2.0E0*4*I_ESP_G3xy_D2y_a+3*2*I_ESP_Dxy_D2y;
    abcd[iGrid*540+47] = 4.0E0*I_ESP_I5xz_D2y_aa-2.0E0*3*I_ESP_G3xz_D2y_a-2.0E0*4*I_ESP_G3xz_D2y_a+3*2*I_ESP_Dxz_D2y;
    abcd[iGrid*540+48] = 4.0E0*I_ESP_I4x2y_D2y_aa-2.0E0*2*I_ESP_G2x2y_D2y_a-2.0E0*3*I_ESP_G2x2y_D2y_a+2*1*I_ESP_D2y_D2y;
    abcd[iGrid*540+49] = 4.0E0*I_ESP_I4xyz_D2y_aa-2.0E0*2*I_ESP_G2xyz_D2y_a-2.0E0*3*I_ESP_G2xyz_D2y_a+2*1*I_ESP_Dyz_D2y;
    abcd[iGrid*540+50] = 4.0E0*I_ESP_I4x2z_D2y_aa-2.0E0*2*I_ESP_G2x2z_D2y_a-2.0E0*3*I_ESP_G2x2z_D2y_a+2*1*I_ESP_D2z_D2y;
    abcd[iGrid*540+51] = 4.0E0*I_ESP_I3x3y_D2y_aa-2.0E0*1*I_ESP_Gx3y_D2y_a-2.0E0*2*I_ESP_Gx3y_D2y_a;
    abcd[iGrid*540+52] = 4.0E0*I_ESP_I3x2yz_D2y_aa-2.0E0*1*I_ESP_Gx2yz_D2y_a-2.0E0*2*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+53] = 4.0E0*I_ESP_I3xy2z_D2y_aa-2.0E0*1*I_ESP_Gxy2z_D2y_a-2.0E0*2*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+54] = 4.0E0*I_ESP_I3x3z_D2y_aa-2.0E0*1*I_ESP_Gx3z_D2y_a-2.0E0*2*I_ESP_Gx3z_D2y_a;
    abcd[iGrid*540+55] = 4.0E0*I_ESP_I2x4y_D2y_aa-2.0E0*1*I_ESP_G4y_D2y_a;
    abcd[iGrid*540+56] = 4.0E0*I_ESP_I2x3yz_D2y_aa-2.0E0*1*I_ESP_G3yz_D2y_a;
    abcd[iGrid*540+57] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2y2z_D2y_a;
    abcd[iGrid*540+58] = 4.0E0*I_ESP_I2xy3z_D2y_aa-2.0E0*1*I_ESP_Gy3z_D2y_a;
    abcd[iGrid*540+59] = 4.0E0*I_ESP_I2x4z_D2y_aa-2.0E0*1*I_ESP_G4z_D2y_a;
    abcd[iGrid*540+60] = 4.0E0*I_ESP_I6x_Dyz_aa-2.0E0*4*I_ESP_G4x_Dyz_a-2.0E0*5*I_ESP_G4x_Dyz_a+4*3*I_ESP_D2x_Dyz;
    abcd[iGrid*540+61] = 4.0E0*I_ESP_I5xy_Dyz_aa-2.0E0*3*I_ESP_G3xy_Dyz_a-2.0E0*4*I_ESP_G3xy_Dyz_a+3*2*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+62] = 4.0E0*I_ESP_I5xz_Dyz_aa-2.0E0*3*I_ESP_G3xz_Dyz_a-2.0E0*4*I_ESP_G3xz_Dyz_a+3*2*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+63] = 4.0E0*I_ESP_I4x2y_Dyz_aa-2.0E0*2*I_ESP_G2x2y_Dyz_a-2.0E0*3*I_ESP_G2x2y_Dyz_a+2*1*I_ESP_D2y_Dyz;
    abcd[iGrid*540+64] = 4.0E0*I_ESP_I4xyz_Dyz_aa-2.0E0*2*I_ESP_G2xyz_Dyz_a-2.0E0*3*I_ESP_G2xyz_Dyz_a+2*1*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+65] = 4.0E0*I_ESP_I4x2z_Dyz_aa-2.0E0*2*I_ESP_G2x2z_Dyz_a-2.0E0*3*I_ESP_G2x2z_Dyz_a+2*1*I_ESP_D2z_Dyz;
    abcd[iGrid*540+66] = 4.0E0*I_ESP_I3x3y_Dyz_aa-2.0E0*1*I_ESP_Gx3y_Dyz_a-2.0E0*2*I_ESP_Gx3y_Dyz_a;
    abcd[iGrid*540+67] = 4.0E0*I_ESP_I3x2yz_Dyz_aa-2.0E0*1*I_ESP_Gx2yz_Dyz_a-2.0E0*2*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+68] = 4.0E0*I_ESP_I3xy2z_Dyz_aa-2.0E0*1*I_ESP_Gxy2z_Dyz_a-2.0E0*2*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+69] = 4.0E0*I_ESP_I3x3z_Dyz_aa-2.0E0*1*I_ESP_Gx3z_Dyz_a-2.0E0*2*I_ESP_Gx3z_Dyz_a;
    abcd[iGrid*540+70] = 4.0E0*I_ESP_I2x4y_Dyz_aa-2.0E0*1*I_ESP_G4y_Dyz_a;
    abcd[iGrid*540+71] = 4.0E0*I_ESP_I2x3yz_Dyz_aa-2.0E0*1*I_ESP_G3yz_Dyz_a;
    abcd[iGrid*540+72] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2y2z_Dyz_a;
    abcd[iGrid*540+73] = 4.0E0*I_ESP_I2xy3z_Dyz_aa-2.0E0*1*I_ESP_Gy3z_Dyz_a;
    abcd[iGrid*540+74] = 4.0E0*I_ESP_I2x4z_Dyz_aa-2.0E0*1*I_ESP_G4z_Dyz_a;
    abcd[iGrid*540+75] = 4.0E0*I_ESP_I6x_D2z_aa-2.0E0*4*I_ESP_G4x_D2z_a-2.0E0*5*I_ESP_G4x_D2z_a+4*3*I_ESP_D2x_D2z;
    abcd[iGrid*540+76] = 4.0E0*I_ESP_I5xy_D2z_aa-2.0E0*3*I_ESP_G3xy_D2z_a-2.0E0*4*I_ESP_G3xy_D2z_a+3*2*I_ESP_Dxy_D2z;
    abcd[iGrid*540+77] = 4.0E0*I_ESP_I5xz_D2z_aa-2.0E0*3*I_ESP_G3xz_D2z_a-2.0E0*4*I_ESP_G3xz_D2z_a+3*2*I_ESP_Dxz_D2z;
    abcd[iGrid*540+78] = 4.0E0*I_ESP_I4x2y_D2z_aa-2.0E0*2*I_ESP_G2x2y_D2z_a-2.0E0*3*I_ESP_G2x2y_D2z_a+2*1*I_ESP_D2y_D2z;
    abcd[iGrid*540+79] = 4.0E0*I_ESP_I4xyz_D2z_aa-2.0E0*2*I_ESP_G2xyz_D2z_a-2.0E0*3*I_ESP_G2xyz_D2z_a+2*1*I_ESP_Dyz_D2z;
    abcd[iGrid*540+80] = 4.0E0*I_ESP_I4x2z_D2z_aa-2.0E0*2*I_ESP_G2x2z_D2z_a-2.0E0*3*I_ESP_G2x2z_D2z_a+2*1*I_ESP_D2z_D2z;
    abcd[iGrid*540+81] = 4.0E0*I_ESP_I3x3y_D2z_aa-2.0E0*1*I_ESP_Gx3y_D2z_a-2.0E0*2*I_ESP_Gx3y_D2z_a;
    abcd[iGrid*540+82] = 4.0E0*I_ESP_I3x2yz_D2z_aa-2.0E0*1*I_ESP_Gx2yz_D2z_a-2.0E0*2*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+83] = 4.0E0*I_ESP_I3xy2z_D2z_aa-2.0E0*1*I_ESP_Gxy2z_D2z_a-2.0E0*2*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+84] = 4.0E0*I_ESP_I3x3z_D2z_aa-2.0E0*1*I_ESP_Gx3z_D2z_a-2.0E0*2*I_ESP_Gx3z_D2z_a;
    abcd[iGrid*540+85] = 4.0E0*I_ESP_I2x4y_D2z_aa-2.0E0*1*I_ESP_G4y_D2z_a;
    abcd[iGrid*540+86] = 4.0E0*I_ESP_I2x3yz_D2z_aa-2.0E0*1*I_ESP_G3yz_D2z_a;
    abcd[iGrid*540+87] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2y2z_D2z_a;
    abcd[iGrid*540+88] = 4.0E0*I_ESP_I2xy3z_D2z_aa-2.0E0*1*I_ESP_Gy3z_D2z_a;
    abcd[iGrid*540+89] = 4.0E0*I_ESP_I2x4z_D2z_aa-2.0E0*1*I_ESP_G4z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+90] = 4.0E0*I_ESP_I5xy_D2x_aa-2.0E0*4*I_ESP_G3xy_D2x_a;
    abcd[iGrid*540+91] = 4.0E0*I_ESP_I4x2y_D2x_aa-2.0E0*1*I_ESP_G4x_D2x_a-2.0E0*3*I_ESP_G2x2y_D2x_a+3*1*I_ESP_D2x_D2x;
    abcd[iGrid*540+92] = 4.0E0*I_ESP_I4xyz_D2x_aa-2.0E0*3*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+93] = 4.0E0*I_ESP_I3x3y_D2x_aa-2.0E0*2*I_ESP_G3xy_D2x_a-2.0E0*2*I_ESP_Gx3y_D2x_a+2*2*I_ESP_Dxy_D2x;
    abcd[iGrid*540+94] = 4.0E0*I_ESP_I3x2yz_D2x_aa-2.0E0*1*I_ESP_G3xz_D2x_a-2.0E0*2*I_ESP_Gx2yz_D2x_a+2*1*I_ESP_Dxz_D2x;
    abcd[iGrid*540+95] = 4.0E0*I_ESP_I3xy2z_D2x_aa-2.0E0*2*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+96] = 4.0E0*I_ESP_I2x4y_D2x_aa-2.0E0*3*I_ESP_G2x2y_D2x_a-2.0E0*1*I_ESP_G4y_D2x_a+3*I_ESP_D2y_D2x;
    abcd[iGrid*540+97] = 4.0E0*I_ESP_I2x3yz_D2x_aa-2.0E0*2*I_ESP_G2xyz_D2x_a-2.0E0*1*I_ESP_G3yz_D2x_a+2*I_ESP_Dyz_D2x;
    abcd[iGrid*540+98] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2x2z_D2x_a-2.0E0*1*I_ESP_G2y2z_D2x_a+1*I_ESP_D2z_D2x;
    abcd[iGrid*540+99] = 4.0E0*I_ESP_I2xy3z_D2x_aa-2.0E0*1*I_ESP_Gy3z_D2x_a;
    abcd[iGrid*540+100] = 4.0E0*I_ESP_Ix5y_D2x_aa-2.0E0*4*I_ESP_Gx3y_D2x_a;
    abcd[iGrid*540+101] = 4.0E0*I_ESP_Ix4yz_D2x_aa-2.0E0*3*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+102] = 4.0E0*I_ESP_Ix3y2z_D2x_aa-2.0E0*2*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+103] = 4.0E0*I_ESP_Ix2y3z_D2x_aa-2.0E0*1*I_ESP_Gx3z_D2x_a;
    abcd[iGrid*540+104] = 4.0E0*I_ESP_Ixy4z_D2x_aa;
    abcd[iGrid*540+105] = 4.0E0*I_ESP_I5xy_Dxy_aa-2.0E0*4*I_ESP_G3xy_Dxy_a;
    abcd[iGrid*540+106] = 4.0E0*I_ESP_I4x2y_Dxy_aa-2.0E0*1*I_ESP_G4x_Dxy_a-2.0E0*3*I_ESP_G2x2y_Dxy_a+3*1*I_ESP_D2x_Dxy;
    abcd[iGrid*540+107] = 4.0E0*I_ESP_I4xyz_Dxy_aa-2.0E0*3*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+108] = 4.0E0*I_ESP_I3x3y_Dxy_aa-2.0E0*2*I_ESP_G3xy_Dxy_a-2.0E0*2*I_ESP_Gx3y_Dxy_a+2*2*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+109] = 4.0E0*I_ESP_I3x2yz_Dxy_aa-2.0E0*1*I_ESP_G3xz_Dxy_a-2.0E0*2*I_ESP_Gx2yz_Dxy_a+2*1*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+110] = 4.0E0*I_ESP_I3xy2z_Dxy_aa-2.0E0*2*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+111] = 4.0E0*I_ESP_I2x4y_Dxy_aa-2.0E0*3*I_ESP_G2x2y_Dxy_a-2.0E0*1*I_ESP_G4y_Dxy_a+3*I_ESP_D2y_Dxy;
    abcd[iGrid*540+112] = 4.0E0*I_ESP_I2x3yz_Dxy_aa-2.0E0*2*I_ESP_G2xyz_Dxy_a-2.0E0*1*I_ESP_G3yz_Dxy_a+2*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+113] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2x2z_Dxy_a-2.0E0*1*I_ESP_G2y2z_Dxy_a+1*I_ESP_D2z_Dxy;
    abcd[iGrid*540+114] = 4.0E0*I_ESP_I2xy3z_Dxy_aa-2.0E0*1*I_ESP_Gy3z_Dxy_a;
    abcd[iGrid*540+115] = 4.0E0*I_ESP_Ix5y_Dxy_aa-2.0E0*4*I_ESP_Gx3y_Dxy_a;
    abcd[iGrid*540+116] = 4.0E0*I_ESP_Ix4yz_Dxy_aa-2.0E0*3*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+117] = 4.0E0*I_ESP_Ix3y2z_Dxy_aa-2.0E0*2*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+118] = 4.0E0*I_ESP_Ix2y3z_Dxy_aa-2.0E0*1*I_ESP_Gx3z_Dxy_a;
    abcd[iGrid*540+119] = 4.0E0*I_ESP_Ixy4z_Dxy_aa;
    abcd[iGrid*540+120] = 4.0E0*I_ESP_I5xy_Dxz_aa-2.0E0*4*I_ESP_G3xy_Dxz_a;
    abcd[iGrid*540+121] = 4.0E0*I_ESP_I4x2y_Dxz_aa-2.0E0*1*I_ESP_G4x_Dxz_a-2.0E0*3*I_ESP_G2x2y_Dxz_a+3*1*I_ESP_D2x_Dxz;
    abcd[iGrid*540+122] = 4.0E0*I_ESP_I4xyz_Dxz_aa-2.0E0*3*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+123] = 4.0E0*I_ESP_I3x3y_Dxz_aa-2.0E0*2*I_ESP_G3xy_Dxz_a-2.0E0*2*I_ESP_Gx3y_Dxz_a+2*2*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+124] = 4.0E0*I_ESP_I3x2yz_Dxz_aa-2.0E0*1*I_ESP_G3xz_Dxz_a-2.0E0*2*I_ESP_Gx2yz_Dxz_a+2*1*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+125] = 4.0E0*I_ESP_I3xy2z_Dxz_aa-2.0E0*2*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+126] = 4.0E0*I_ESP_I2x4y_Dxz_aa-2.0E0*3*I_ESP_G2x2y_Dxz_a-2.0E0*1*I_ESP_G4y_Dxz_a+3*I_ESP_D2y_Dxz;
    abcd[iGrid*540+127] = 4.0E0*I_ESP_I2x3yz_Dxz_aa-2.0E0*2*I_ESP_G2xyz_Dxz_a-2.0E0*1*I_ESP_G3yz_Dxz_a+2*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+128] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2x2z_Dxz_a-2.0E0*1*I_ESP_G2y2z_Dxz_a+1*I_ESP_D2z_Dxz;
    abcd[iGrid*540+129] = 4.0E0*I_ESP_I2xy3z_Dxz_aa-2.0E0*1*I_ESP_Gy3z_Dxz_a;
    abcd[iGrid*540+130] = 4.0E0*I_ESP_Ix5y_Dxz_aa-2.0E0*4*I_ESP_Gx3y_Dxz_a;
    abcd[iGrid*540+131] = 4.0E0*I_ESP_Ix4yz_Dxz_aa-2.0E0*3*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+132] = 4.0E0*I_ESP_Ix3y2z_Dxz_aa-2.0E0*2*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+133] = 4.0E0*I_ESP_Ix2y3z_Dxz_aa-2.0E0*1*I_ESP_Gx3z_Dxz_a;
    abcd[iGrid*540+134] = 4.0E0*I_ESP_Ixy4z_Dxz_aa;
    abcd[iGrid*540+135] = 4.0E0*I_ESP_I5xy_D2y_aa-2.0E0*4*I_ESP_G3xy_D2y_a;
    abcd[iGrid*540+136] = 4.0E0*I_ESP_I4x2y_D2y_aa-2.0E0*1*I_ESP_G4x_D2y_a-2.0E0*3*I_ESP_G2x2y_D2y_a+3*1*I_ESP_D2x_D2y;
    abcd[iGrid*540+137] = 4.0E0*I_ESP_I4xyz_D2y_aa-2.0E0*3*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+138] = 4.0E0*I_ESP_I3x3y_D2y_aa-2.0E0*2*I_ESP_G3xy_D2y_a-2.0E0*2*I_ESP_Gx3y_D2y_a+2*2*I_ESP_Dxy_D2y;
    abcd[iGrid*540+139] = 4.0E0*I_ESP_I3x2yz_D2y_aa-2.0E0*1*I_ESP_G3xz_D2y_a-2.0E0*2*I_ESP_Gx2yz_D2y_a+2*1*I_ESP_Dxz_D2y;
    abcd[iGrid*540+140] = 4.0E0*I_ESP_I3xy2z_D2y_aa-2.0E0*2*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+141] = 4.0E0*I_ESP_I2x4y_D2y_aa-2.0E0*3*I_ESP_G2x2y_D2y_a-2.0E0*1*I_ESP_G4y_D2y_a+3*I_ESP_D2y_D2y;
    abcd[iGrid*540+142] = 4.0E0*I_ESP_I2x3yz_D2y_aa-2.0E0*2*I_ESP_G2xyz_D2y_a-2.0E0*1*I_ESP_G3yz_D2y_a+2*I_ESP_Dyz_D2y;
    abcd[iGrid*540+143] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2x2z_D2y_a-2.0E0*1*I_ESP_G2y2z_D2y_a+1*I_ESP_D2z_D2y;
    abcd[iGrid*540+144] = 4.0E0*I_ESP_I2xy3z_D2y_aa-2.0E0*1*I_ESP_Gy3z_D2y_a;
    abcd[iGrid*540+145] = 4.0E0*I_ESP_Ix5y_D2y_aa-2.0E0*4*I_ESP_Gx3y_D2y_a;
    abcd[iGrid*540+146] = 4.0E0*I_ESP_Ix4yz_D2y_aa-2.0E0*3*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+147] = 4.0E0*I_ESP_Ix3y2z_D2y_aa-2.0E0*2*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+148] = 4.0E0*I_ESP_Ix2y3z_D2y_aa-2.0E0*1*I_ESP_Gx3z_D2y_a;
    abcd[iGrid*540+149] = 4.0E0*I_ESP_Ixy4z_D2y_aa;
    abcd[iGrid*540+150] = 4.0E0*I_ESP_I5xy_Dyz_aa-2.0E0*4*I_ESP_G3xy_Dyz_a;
    abcd[iGrid*540+151] = 4.0E0*I_ESP_I4x2y_Dyz_aa-2.0E0*1*I_ESP_G4x_Dyz_a-2.0E0*3*I_ESP_G2x2y_Dyz_a+3*1*I_ESP_D2x_Dyz;
    abcd[iGrid*540+152] = 4.0E0*I_ESP_I4xyz_Dyz_aa-2.0E0*3*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+153] = 4.0E0*I_ESP_I3x3y_Dyz_aa-2.0E0*2*I_ESP_G3xy_Dyz_a-2.0E0*2*I_ESP_Gx3y_Dyz_a+2*2*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+154] = 4.0E0*I_ESP_I3x2yz_Dyz_aa-2.0E0*1*I_ESP_G3xz_Dyz_a-2.0E0*2*I_ESP_Gx2yz_Dyz_a+2*1*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+155] = 4.0E0*I_ESP_I3xy2z_Dyz_aa-2.0E0*2*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+156] = 4.0E0*I_ESP_I2x4y_Dyz_aa-2.0E0*3*I_ESP_G2x2y_Dyz_a-2.0E0*1*I_ESP_G4y_Dyz_a+3*I_ESP_D2y_Dyz;
    abcd[iGrid*540+157] = 4.0E0*I_ESP_I2x3yz_Dyz_aa-2.0E0*2*I_ESP_G2xyz_Dyz_a-2.0E0*1*I_ESP_G3yz_Dyz_a+2*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+158] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2x2z_Dyz_a-2.0E0*1*I_ESP_G2y2z_Dyz_a+1*I_ESP_D2z_Dyz;
    abcd[iGrid*540+159] = 4.0E0*I_ESP_I2xy3z_Dyz_aa-2.0E0*1*I_ESP_Gy3z_Dyz_a;
    abcd[iGrid*540+160] = 4.0E0*I_ESP_Ix5y_Dyz_aa-2.0E0*4*I_ESP_Gx3y_Dyz_a;
    abcd[iGrid*540+161] = 4.0E0*I_ESP_Ix4yz_Dyz_aa-2.0E0*3*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+162] = 4.0E0*I_ESP_Ix3y2z_Dyz_aa-2.0E0*2*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+163] = 4.0E0*I_ESP_Ix2y3z_Dyz_aa-2.0E0*1*I_ESP_Gx3z_Dyz_a;
    abcd[iGrid*540+164] = 4.0E0*I_ESP_Ixy4z_Dyz_aa;
    abcd[iGrid*540+165] = 4.0E0*I_ESP_I5xy_D2z_aa-2.0E0*4*I_ESP_G3xy_D2z_a;
    abcd[iGrid*540+166] = 4.0E0*I_ESP_I4x2y_D2z_aa-2.0E0*1*I_ESP_G4x_D2z_a-2.0E0*3*I_ESP_G2x2y_D2z_a+3*1*I_ESP_D2x_D2z;
    abcd[iGrid*540+167] = 4.0E0*I_ESP_I4xyz_D2z_aa-2.0E0*3*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+168] = 4.0E0*I_ESP_I3x3y_D2z_aa-2.0E0*2*I_ESP_G3xy_D2z_a-2.0E0*2*I_ESP_Gx3y_D2z_a+2*2*I_ESP_Dxy_D2z;
    abcd[iGrid*540+169] = 4.0E0*I_ESP_I3x2yz_D2z_aa-2.0E0*1*I_ESP_G3xz_D2z_a-2.0E0*2*I_ESP_Gx2yz_D2z_a+2*1*I_ESP_Dxz_D2z;
    abcd[iGrid*540+170] = 4.0E0*I_ESP_I3xy2z_D2z_aa-2.0E0*2*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+171] = 4.0E0*I_ESP_I2x4y_D2z_aa-2.0E0*3*I_ESP_G2x2y_D2z_a-2.0E0*1*I_ESP_G4y_D2z_a+3*I_ESP_D2y_D2z;
    abcd[iGrid*540+172] = 4.0E0*I_ESP_I2x3yz_D2z_aa-2.0E0*2*I_ESP_G2xyz_D2z_a-2.0E0*1*I_ESP_G3yz_D2z_a+2*I_ESP_Dyz_D2z;
    abcd[iGrid*540+173] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2x2z_D2z_a-2.0E0*1*I_ESP_G2y2z_D2z_a+1*I_ESP_D2z_D2z;
    abcd[iGrid*540+174] = 4.0E0*I_ESP_I2xy3z_D2z_aa-2.0E0*1*I_ESP_Gy3z_D2z_a;
    abcd[iGrid*540+175] = 4.0E0*I_ESP_Ix5y_D2z_aa-2.0E0*4*I_ESP_Gx3y_D2z_a;
    abcd[iGrid*540+176] = 4.0E0*I_ESP_Ix4yz_D2z_aa-2.0E0*3*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+177] = 4.0E0*I_ESP_Ix3y2z_D2z_aa-2.0E0*2*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+178] = 4.0E0*I_ESP_Ix2y3z_D2z_aa-2.0E0*1*I_ESP_Gx3z_D2z_a;
    abcd[iGrid*540+179] = 4.0E0*I_ESP_Ixy4z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+180] = 4.0E0*I_ESP_I5xz_D2x_aa-2.0E0*4*I_ESP_G3xz_D2x_a;
    abcd[iGrid*540+181] = 4.0E0*I_ESP_I4xyz_D2x_aa-2.0E0*3*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+182] = 4.0E0*I_ESP_I4x2z_D2x_aa-2.0E0*1*I_ESP_G4x_D2x_a-2.0E0*3*I_ESP_G2x2z_D2x_a+3*1*I_ESP_D2x_D2x;
    abcd[iGrid*540+183] = 4.0E0*I_ESP_I3x2yz_D2x_aa-2.0E0*2*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+184] = 4.0E0*I_ESP_I3xy2z_D2x_aa-2.0E0*1*I_ESP_G3xy_D2x_a-2.0E0*2*I_ESP_Gxy2z_D2x_a+2*1*I_ESP_Dxy_D2x;
    abcd[iGrid*540+185] = 4.0E0*I_ESP_I3x3z_D2x_aa-2.0E0*2*I_ESP_G3xz_D2x_a-2.0E0*2*I_ESP_Gx3z_D2x_a+2*2*I_ESP_Dxz_D2x;
    abcd[iGrid*540+186] = 4.0E0*I_ESP_I2x3yz_D2x_aa-2.0E0*1*I_ESP_G3yz_D2x_a;
    abcd[iGrid*540+187] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2x2y_D2x_a-2.0E0*1*I_ESP_G2y2z_D2x_a+1*I_ESP_D2y_D2x;
    abcd[iGrid*540+188] = 4.0E0*I_ESP_I2xy3z_D2x_aa-2.0E0*2*I_ESP_G2xyz_D2x_a-2.0E0*1*I_ESP_Gy3z_D2x_a+2*I_ESP_Dyz_D2x;
    abcd[iGrid*540+189] = 4.0E0*I_ESP_I2x4z_D2x_aa-2.0E0*3*I_ESP_G2x2z_D2x_a-2.0E0*1*I_ESP_G4z_D2x_a+3*I_ESP_D2z_D2x;
    abcd[iGrid*540+190] = 4.0E0*I_ESP_Ix4yz_D2x_aa;
    abcd[iGrid*540+191] = 4.0E0*I_ESP_Ix3y2z_D2x_aa-2.0E0*1*I_ESP_Gx3y_D2x_a;
    abcd[iGrid*540+192] = 4.0E0*I_ESP_Ix2y3z_D2x_aa-2.0E0*2*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+193] = 4.0E0*I_ESP_Ixy4z_D2x_aa-2.0E0*3*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+194] = 4.0E0*I_ESP_Ix5z_D2x_aa-2.0E0*4*I_ESP_Gx3z_D2x_a;
    abcd[iGrid*540+195] = 4.0E0*I_ESP_I5xz_Dxy_aa-2.0E0*4*I_ESP_G3xz_Dxy_a;
    abcd[iGrid*540+196] = 4.0E0*I_ESP_I4xyz_Dxy_aa-2.0E0*3*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+197] = 4.0E0*I_ESP_I4x2z_Dxy_aa-2.0E0*1*I_ESP_G4x_Dxy_a-2.0E0*3*I_ESP_G2x2z_Dxy_a+3*1*I_ESP_D2x_Dxy;
    abcd[iGrid*540+198] = 4.0E0*I_ESP_I3x2yz_Dxy_aa-2.0E0*2*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+199] = 4.0E0*I_ESP_I3xy2z_Dxy_aa-2.0E0*1*I_ESP_G3xy_Dxy_a-2.0E0*2*I_ESP_Gxy2z_Dxy_a+2*1*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+200] = 4.0E0*I_ESP_I3x3z_Dxy_aa-2.0E0*2*I_ESP_G3xz_Dxy_a-2.0E0*2*I_ESP_Gx3z_Dxy_a+2*2*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+201] = 4.0E0*I_ESP_I2x3yz_Dxy_aa-2.0E0*1*I_ESP_G3yz_Dxy_a;
    abcd[iGrid*540+202] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2x2y_Dxy_a-2.0E0*1*I_ESP_G2y2z_Dxy_a+1*I_ESP_D2y_Dxy;
    abcd[iGrid*540+203] = 4.0E0*I_ESP_I2xy3z_Dxy_aa-2.0E0*2*I_ESP_G2xyz_Dxy_a-2.0E0*1*I_ESP_Gy3z_Dxy_a+2*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+204] = 4.0E0*I_ESP_I2x4z_Dxy_aa-2.0E0*3*I_ESP_G2x2z_Dxy_a-2.0E0*1*I_ESP_G4z_Dxy_a+3*I_ESP_D2z_Dxy;
    abcd[iGrid*540+205] = 4.0E0*I_ESP_Ix4yz_Dxy_aa;
    abcd[iGrid*540+206] = 4.0E0*I_ESP_Ix3y2z_Dxy_aa-2.0E0*1*I_ESP_Gx3y_Dxy_a;
    abcd[iGrid*540+207] = 4.0E0*I_ESP_Ix2y3z_Dxy_aa-2.0E0*2*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+208] = 4.0E0*I_ESP_Ixy4z_Dxy_aa-2.0E0*3*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+209] = 4.0E0*I_ESP_Ix5z_Dxy_aa-2.0E0*4*I_ESP_Gx3z_Dxy_a;
    abcd[iGrid*540+210] = 4.0E0*I_ESP_I5xz_Dxz_aa-2.0E0*4*I_ESP_G3xz_Dxz_a;
    abcd[iGrid*540+211] = 4.0E0*I_ESP_I4xyz_Dxz_aa-2.0E0*3*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+212] = 4.0E0*I_ESP_I4x2z_Dxz_aa-2.0E0*1*I_ESP_G4x_Dxz_a-2.0E0*3*I_ESP_G2x2z_Dxz_a+3*1*I_ESP_D2x_Dxz;
    abcd[iGrid*540+213] = 4.0E0*I_ESP_I3x2yz_Dxz_aa-2.0E0*2*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+214] = 4.0E0*I_ESP_I3xy2z_Dxz_aa-2.0E0*1*I_ESP_G3xy_Dxz_a-2.0E0*2*I_ESP_Gxy2z_Dxz_a+2*1*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+215] = 4.0E0*I_ESP_I3x3z_Dxz_aa-2.0E0*2*I_ESP_G3xz_Dxz_a-2.0E0*2*I_ESP_Gx3z_Dxz_a+2*2*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+216] = 4.0E0*I_ESP_I2x3yz_Dxz_aa-2.0E0*1*I_ESP_G3yz_Dxz_a;
    abcd[iGrid*540+217] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2x2y_Dxz_a-2.0E0*1*I_ESP_G2y2z_Dxz_a+1*I_ESP_D2y_Dxz;
    abcd[iGrid*540+218] = 4.0E0*I_ESP_I2xy3z_Dxz_aa-2.0E0*2*I_ESP_G2xyz_Dxz_a-2.0E0*1*I_ESP_Gy3z_Dxz_a+2*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+219] = 4.0E0*I_ESP_I2x4z_Dxz_aa-2.0E0*3*I_ESP_G2x2z_Dxz_a-2.0E0*1*I_ESP_G4z_Dxz_a+3*I_ESP_D2z_Dxz;
    abcd[iGrid*540+220] = 4.0E0*I_ESP_Ix4yz_Dxz_aa;
    abcd[iGrid*540+221] = 4.0E0*I_ESP_Ix3y2z_Dxz_aa-2.0E0*1*I_ESP_Gx3y_Dxz_a;
    abcd[iGrid*540+222] = 4.0E0*I_ESP_Ix2y3z_Dxz_aa-2.0E0*2*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+223] = 4.0E0*I_ESP_Ixy4z_Dxz_aa-2.0E0*3*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+224] = 4.0E0*I_ESP_Ix5z_Dxz_aa-2.0E0*4*I_ESP_Gx3z_Dxz_a;
    abcd[iGrid*540+225] = 4.0E0*I_ESP_I5xz_D2y_aa-2.0E0*4*I_ESP_G3xz_D2y_a;
    abcd[iGrid*540+226] = 4.0E0*I_ESP_I4xyz_D2y_aa-2.0E0*3*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+227] = 4.0E0*I_ESP_I4x2z_D2y_aa-2.0E0*1*I_ESP_G4x_D2y_a-2.0E0*3*I_ESP_G2x2z_D2y_a+3*1*I_ESP_D2x_D2y;
    abcd[iGrid*540+228] = 4.0E0*I_ESP_I3x2yz_D2y_aa-2.0E0*2*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+229] = 4.0E0*I_ESP_I3xy2z_D2y_aa-2.0E0*1*I_ESP_G3xy_D2y_a-2.0E0*2*I_ESP_Gxy2z_D2y_a+2*1*I_ESP_Dxy_D2y;
    abcd[iGrid*540+230] = 4.0E0*I_ESP_I3x3z_D2y_aa-2.0E0*2*I_ESP_G3xz_D2y_a-2.0E0*2*I_ESP_Gx3z_D2y_a+2*2*I_ESP_Dxz_D2y;
    abcd[iGrid*540+231] = 4.0E0*I_ESP_I2x3yz_D2y_aa-2.0E0*1*I_ESP_G3yz_D2y_a;
    abcd[iGrid*540+232] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2x2y_D2y_a-2.0E0*1*I_ESP_G2y2z_D2y_a+1*I_ESP_D2y_D2y;
    abcd[iGrid*540+233] = 4.0E0*I_ESP_I2xy3z_D2y_aa-2.0E0*2*I_ESP_G2xyz_D2y_a-2.0E0*1*I_ESP_Gy3z_D2y_a+2*I_ESP_Dyz_D2y;
    abcd[iGrid*540+234] = 4.0E0*I_ESP_I2x4z_D2y_aa-2.0E0*3*I_ESP_G2x2z_D2y_a-2.0E0*1*I_ESP_G4z_D2y_a+3*I_ESP_D2z_D2y;
    abcd[iGrid*540+235] = 4.0E0*I_ESP_Ix4yz_D2y_aa;
    abcd[iGrid*540+236] = 4.0E0*I_ESP_Ix3y2z_D2y_aa-2.0E0*1*I_ESP_Gx3y_D2y_a;
    abcd[iGrid*540+237] = 4.0E0*I_ESP_Ix2y3z_D2y_aa-2.0E0*2*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+238] = 4.0E0*I_ESP_Ixy4z_D2y_aa-2.0E0*3*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+239] = 4.0E0*I_ESP_Ix5z_D2y_aa-2.0E0*4*I_ESP_Gx3z_D2y_a;
    abcd[iGrid*540+240] = 4.0E0*I_ESP_I5xz_Dyz_aa-2.0E0*4*I_ESP_G3xz_Dyz_a;
    abcd[iGrid*540+241] = 4.0E0*I_ESP_I4xyz_Dyz_aa-2.0E0*3*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+242] = 4.0E0*I_ESP_I4x2z_Dyz_aa-2.0E0*1*I_ESP_G4x_Dyz_a-2.0E0*3*I_ESP_G2x2z_Dyz_a+3*1*I_ESP_D2x_Dyz;
    abcd[iGrid*540+243] = 4.0E0*I_ESP_I3x2yz_Dyz_aa-2.0E0*2*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+244] = 4.0E0*I_ESP_I3xy2z_Dyz_aa-2.0E0*1*I_ESP_G3xy_Dyz_a-2.0E0*2*I_ESP_Gxy2z_Dyz_a+2*1*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+245] = 4.0E0*I_ESP_I3x3z_Dyz_aa-2.0E0*2*I_ESP_G3xz_Dyz_a-2.0E0*2*I_ESP_Gx3z_Dyz_a+2*2*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+246] = 4.0E0*I_ESP_I2x3yz_Dyz_aa-2.0E0*1*I_ESP_G3yz_Dyz_a;
    abcd[iGrid*540+247] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2x2y_Dyz_a-2.0E0*1*I_ESP_G2y2z_Dyz_a+1*I_ESP_D2y_Dyz;
    abcd[iGrid*540+248] = 4.0E0*I_ESP_I2xy3z_Dyz_aa-2.0E0*2*I_ESP_G2xyz_Dyz_a-2.0E0*1*I_ESP_Gy3z_Dyz_a+2*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+249] = 4.0E0*I_ESP_I2x4z_Dyz_aa-2.0E0*3*I_ESP_G2x2z_Dyz_a-2.0E0*1*I_ESP_G4z_Dyz_a+3*I_ESP_D2z_Dyz;
    abcd[iGrid*540+250] = 4.0E0*I_ESP_Ix4yz_Dyz_aa;
    abcd[iGrid*540+251] = 4.0E0*I_ESP_Ix3y2z_Dyz_aa-2.0E0*1*I_ESP_Gx3y_Dyz_a;
    abcd[iGrid*540+252] = 4.0E0*I_ESP_Ix2y3z_Dyz_aa-2.0E0*2*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+253] = 4.0E0*I_ESP_Ixy4z_Dyz_aa-2.0E0*3*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+254] = 4.0E0*I_ESP_Ix5z_Dyz_aa-2.0E0*4*I_ESP_Gx3z_Dyz_a;
    abcd[iGrid*540+255] = 4.0E0*I_ESP_I5xz_D2z_aa-2.0E0*4*I_ESP_G3xz_D2z_a;
    abcd[iGrid*540+256] = 4.0E0*I_ESP_I4xyz_D2z_aa-2.0E0*3*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+257] = 4.0E0*I_ESP_I4x2z_D2z_aa-2.0E0*1*I_ESP_G4x_D2z_a-2.0E0*3*I_ESP_G2x2z_D2z_a+3*1*I_ESP_D2x_D2z;
    abcd[iGrid*540+258] = 4.0E0*I_ESP_I3x2yz_D2z_aa-2.0E0*2*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+259] = 4.0E0*I_ESP_I3xy2z_D2z_aa-2.0E0*1*I_ESP_G3xy_D2z_a-2.0E0*2*I_ESP_Gxy2z_D2z_a+2*1*I_ESP_Dxy_D2z;
    abcd[iGrid*540+260] = 4.0E0*I_ESP_I3x3z_D2z_aa-2.0E0*2*I_ESP_G3xz_D2z_a-2.0E0*2*I_ESP_Gx3z_D2z_a+2*2*I_ESP_Dxz_D2z;
    abcd[iGrid*540+261] = 4.0E0*I_ESP_I2x3yz_D2z_aa-2.0E0*1*I_ESP_G3yz_D2z_a;
    abcd[iGrid*540+262] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2x2y_D2z_a-2.0E0*1*I_ESP_G2y2z_D2z_a+1*I_ESP_D2y_D2z;
    abcd[iGrid*540+263] = 4.0E0*I_ESP_I2xy3z_D2z_aa-2.0E0*2*I_ESP_G2xyz_D2z_a-2.0E0*1*I_ESP_Gy3z_D2z_a+2*I_ESP_Dyz_D2z;
    abcd[iGrid*540+264] = 4.0E0*I_ESP_I2x4z_D2z_aa-2.0E0*3*I_ESP_G2x2z_D2z_a-2.0E0*1*I_ESP_G4z_D2z_a+3*I_ESP_D2z_D2z;
    abcd[iGrid*540+265] = 4.0E0*I_ESP_Ix4yz_D2z_aa;
    abcd[iGrid*540+266] = 4.0E0*I_ESP_Ix3y2z_D2z_aa-2.0E0*1*I_ESP_Gx3y_D2z_a;
    abcd[iGrid*540+267] = 4.0E0*I_ESP_Ix2y3z_D2z_aa-2.0E0*2*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+268] = 4.0E0*I_ESP_Ixy4z_D2z_aa-2.0E0*3*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+269] = 4.0E0*I_ESP_Ix5z_D2z_aa-2.0E0*4*I_ESP_Gx3z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+270] = 4.0E0*I_ESP_I4x2y_D2x_aa-2.0E0*1*I_ESP_G4x_D2x_a;
    abcd[iGrid*540+271] = 4.0E0*I_ESP_I3x3y_D2x_aa-2.0E0*1*I_ESP_G3xy_D2x_a-2.0E0*2*I_ESP_G3xy_D2x_a;
    abcd[iGrid*540+272] = 4.0E0*I_ESP_I3x2yz_D2x_aa-2.0E0*1*I_ESP_G3xz_D2x_a;
    abcd[iGrid*540+273] = 4.0E0*I_ESP_I2x4y_D2x_aa-2.0E0*2*I_ESP_G2x2y_D2x_a-2.0E0*3*I_ESP_G2x2y_D2x_a+2*1*I_ESP_D2x_D2x;
    abcd[iGrid*540+274] = 4.0E0*I_ESP_I2x3yz_D2x_aa-2.0E0*1*I_ESP_G2xyz_D2x_a-2.0E0*2*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+275] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2x2z_D2x_a;
    abcd[iGrid*540+276] = 4.0E0*I_ESP_Ix5y_D2x_aa-2.0E0*3*I_ESP_Gx3y_D2x_a-2.0E0*4*I_ESP_Gx3y_D2x_a+3*2*I_ESP_Dxy_D2x;
    abcd[iGrid*540+277] = 4.0E0*I_ESP_Ix4yz_D2x_aa-2.0E0*2*I_ESP_Gx2yz_D2x_a-2.0E0*3*I_ESP_Gx2yz_D2x_a+2*1*I_ESP_Dxz_D2x;
    abcd[iGrid*540+278] = 4.0E0*I_ESP_Ix3y2z_D2x_aa-2.0E0*1*I_ESP_Gxy2z_D2x_a-2.0E0*2*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+279] = 4.0E0*I_ESP_Ix2y3z_D2x_aa-2.0E0*1*I_ESP_Gx3z_D2x_a;
    abcd[iGrid*540+280] = 4.0E0*I_ESP_I6y_D2x_aa-2.0E0*4*I_ESP_G4y_D2x_a-2.0E0*5*I_ESP_G4y_D2x_a+4*3*I_ESP_D2y_D2x;
    abcd[iGrid*540+281] = 4.0E0*I_ESP_I5yz_D2x_aa-2.0E0*3*I_ESP_G3yz_D2x_a-2.0E0*4*I_ESP_G3yz_D2x_a+3*2*I_ESP_Dyz_D2x;
    abcd[iGrid*540+282] = 4.0E0*I_ESP_I4y2z_D2x_aa-2.0E0*2*I_ESP_G2y2z_D2x_a-2.0E0*3*I_ESP_G2y2z_D2x_a+2*1*I_ESP_D2z_D2x;
    abcd[iGrid*540+283] = 4.0E0*I_ESP_I3y3z_D2x_aa-2.0E0*1*I_ESP_Gy3z_D2x_a-2.0E0*2*I_ESP_Gy3z_D2x_a;
    abcd[iGrid*540+284] = 4.0E0*I_ESP_I2y4z_D2x_aa-2.0E0*1*I_ESP_G4z_D2x_a;
    abcd[iGrid*540+285] = 4.0E0*I_ESP_I4x2y_Dxy_aa-2.0E0*1*I_ESP_G4x_Dxy_a;
    abcd[iGrid*540+286] = 4.0E0*I_ESP_I3x3y_Dxy_aa-2.0E0*1*I_ESP_G3xy_Dxy_a-2.0E0*2*I_ESP_G3xy_Dxy_a;
    abcd[iGrid*540+287] = 4.0E0*I_ESP_I3x2yz_Dxy_aa-2.0E0*1*I_ESP_G3xz_Dxy_a;
    abcd[iGrid*540+288] = 4.0E0*I_ESP_I2x4y_Dxy_aa-2.0E0*2*I_ESP_G2x2y_Dxy_a-2.0E0*3*I_ESP_G2x2y_Dxy_a+2*1*I_ESP_D2x_Dxy;
    abcd[iGrid*540+289] = 4.0E0*I_ESP_I2x3yz_Dxy_aa-2.0E0*1*I_ESP_G2xyz_Dxy_a-2.0E0*2*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+290] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2x2z_Dxy_a;
    abcd[iGrid*540+291] = 4.0E0*I_ESP_Ix5y_Dxy_aa-2.0E0*3*I_ESP_Gx3y_Dxy_a-2.0E0*4*I_ESP_Gx3y_Dxy_a+3*2*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+292] = 4.0E0*I_ESP_Ix4yz_Dxy_aa-2.0E0*2*I_ESP_Gx2yz_Dxy_a-2.0E0*3*I_ESP_Gx2yz_Dxy_a+2*1*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+293] = 4.0E0*I_ESP_Ix3y2z_Dxy_aa-2.0E0*1*I_ESP_Gxy2z_Dxy_a-2.0E0*2*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+294] = 4.0E0*I_ESP_Ix2y3z_Dxy_aa-2.0E0*1*I_ESP_Gx3z_Dxy_a;
    abcd[iGrid*540+295] = 4.0E0*I_ESP_I6y_Dxy_aa-2.0E0*4*I_ESP_G4y_Dxy_a-2.0E0*5*I_ESP_G4y_Dxy_a+4*3*I_ESP_D2y_Dxy;
    abcd[iGrid*540+296] = 4.0E0*I_ESP_I5yz_Dxy_aa-2.0E0*3*I_ESP_G3yz_Dxy_a-2.0E0*4*I_ESP_G3yz_Dxy_a+3*2*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+297] = 4.0E0*I_ESP_I4y2z_Dxy_aa-2.0E0*2*I_ESP_G2y2z_Dxy_a-2.0E0*3*I_ESP_G2y2z_Dxy_a+2*1*I_ESP_D2z_Dxy;
    abcd[iGrid*540+298] = 4.0E0*I_ESP_I3y3z_Dxy_aa-2.0E0*1*I_ESP_Gy3z_Dxy_a-2.0E0*2*I_ESP_Gy3z_Dxy_a;
    abcd[iGrid*540+299] = 4.0E0*I_ESP_I2y4z_Dxy_aa-2.0E0*1*I_ESP_G4z_Dxy_a;
    abcd[iGrid*540+300] = 4.0E0*I_ESP_I4x2y_Dxz_aa-2.0E0*1*I_ESP_G4x_Dxz_a;
    abcd[iGrid*540+301] = 4.0E0*I_ESP_I3x3y_Dxz_aa-2.0E0*1*I_ESP_G3xy_Dxz_a-2.0E0*2*I_ESP_G3xy_Dxz_a;
    abcd[iGrid*540+302] = 4.0E0*I_ESP_I3x2yz_Dxz_aa-2.0E0*1*I_ESP_G3xz_Dxz_a;
    abcd[iGrid*540+303] = 4.0E0*I_ESP_I2x4y_Dxz_aa-2.0E0*2*I_ESP_G2x2y_Dxz_a-2.0E0*3*I_ESP_G2x2y_Dxz_a+2*1*I_ESP_D2x_Dxz;
    abcd[iGrid*540+304] = 4.0E0*I_ESP_I2x3yz_Dxz_aa-2.0E0*1*I_ESP_G2xyz_Dxz_a-2.0E0*2*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+305] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2x2z_Dxz_a;
    abcd[iGrid*540+306] = 4.0E0*I_ESP_Ix5y_Dxz_aa-2.0E0*3*I_ESP_Gx3y_Dxz_a-2.0E0*4*I_ESP_Gx3y_Dxz_a+3*2*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+307] = 4.0E0*I_ESP_Ix4yz_Dxz_aa-2.0E0*2*I_ESP_Gx2yz_Dxz_a-2.0E0*3*I_ESP_Gx2yz_Dxz_a+2*1*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+308] = 4.0E0*I_ESP_Ix3y2z_Dxz_aa-2.0E0*1*I_ESP_Gxy2z_Dxz_a-2.0E0*2*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+309] = 4.0E0*I_ESP_Ix2y3z_Dxz_aa-2.0E0*1*I_ESP_Gx3z_Dxz_a;
    abcd[iGrid*540+310] = 4.0E0*I_ESP_I6y_Dxz_aa-2.0E0*4*I_ESP_G4y_Dxz_a-2.0E0*5*I_ESP_G4y_Dxz_a+4*3*I_ESP_D2y_Dxz;
    abcd[iGrid*540+311] = 4.0E0*I_ESP_I5yz_Dxz_aa-2.0E0*3*I_ESP_G3yz_Dxz_a-2.0E0*4*I_ESP_G3yz_Dxz_a+3*2*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+312] = 4.0E0*I_ESP_I4y2z_Dxz_aa-2.0E0*2*I_ESP_G2y2z_Dxz_a-2.0E0*3*I_ESP_G2y2z_Dxz_a+2*1*I_ESP_D2z_Dxz;
    abcd[iGrid*540+313] = 4.0E0*I_ESP_I3y3z_Dxz_aa-2.0E0*1*I_ESP_Gy3z_Dxz_a-2.0E0*2*I_ESP_Gy3z_Dxz_a;
    abcd[iGrid*540+314] = 4.0E0*I_ESP_I2y4z_Dxz_aa-2.0E0*1*I_ESP_G4z_Dxz_a;
    abcd[iGrid*540+315] = 4.0E0*I_ESP_I4x2y_D2y_aa-2.0E0*1*I_ESP_G4x_D2y_a;
    abcd[iGrid*540+316] = 4.0E0*I_ESP_I3x3y_D2y_aa-2.0E0*1*I_ESP_G3xy_D2y_a-2.0E0*2*I_ESP_G3xy_D2y_a;
    abcd[iGrid*540+317] = 4.0E0*I_ESP_I3x2yz_D2y_aa-2.0E0*1*I_ESP_G3xz_D2y_a;
    abcd[iGrid*540+318] = 4.0E0*I_ESP_I2x4y_D2y_aa-2.0E0*2*I_ESP_G2x2y_D2y_a-2.0E0*3*I_ESP_G2x2y_D2y_a+2*1*I_ESP_D2x_D2y;
    abcd[iGrid*540+319] = 4.0E0*I_ESP_I2x3yz_D2y_aa-2.0E0*1*I_ESP_G2xyz_D2y_a-2.0E0*2*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+320] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2x2z_D2y_a;
    abcd[iGrid*540+321] = 4.0E0*I_ESP_Ix5y_D2y_aa-2.0E0*3*I_ESP_Gx3y_D2y_a-2.0E0*4*I_ESP_Gx3y_D2y_a+3*2*I_ESP_Dxy_D2y;
    abcd[iGrid*540+322] = 4.0E0*I_ESP_Ix4yz_D2y_aa-2.0E0*2*I_ESP_Gx2yz_D2y_a-2.0E0*3*I_ESP_Gx2yz_D2y_a+2*1*I_ESP_Dxz_D2y;
    abcd[iGrid*540+323] = 4.0E0*I_ESP_Ix3y2z_D2y_aa-2.0E0*1*I_ESP_Gxy2z_D2y_a-2.0E0*2*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+324] = 4.0E0*I_ESP_Ix2y3z_D2y_aa-2.0E0*1*I_ESP_Gx3z_D2y_a;
    abcd[iGrid*540+325] = 4.0E0*I_ESP_I6y_D2y_aa-2.0E0*4*I_ESP_G4y_D2y_a-2.0E0*5*I_ESP_G4y_D2y_a+4*3*I_ESP_D2y_D2y;
    abcd[iGrid*540+326] = 4.0E0*I_ESP_I5yz_D2y_aa-2.0E0*3*I_ESP_G3yz_D2y_a-2.0E0*4*I_ESP_G3yz_D2y_a+3*2*I_ESP_Dyz_D2y;
    abcd[iGrid*540+327] = 4.0E0*I_ESP_I4y2z_D2y_aa-2.0E0*2*I_ESP_G2y2z_D2y_a-2.0E0*3*I_ESP_G2y2z_D2y_a+2*1*I_ESP_D2z_D2y;
    abcd[iGrid*540+328] = 4.0E0*I_ESP_I3y3z_D2y_aa-2.0E0*1*I_ESP_Gy3z_D2y_a-2.0E0*2*I_ESP_Gy3z_D2y_a;
    abcd[iGrid*540+329] = 4.0E0*I_ESP_I2y4z_D2y_aa-2.0E0*1*I_ESP_G4z_D2y_a;
    abcd[iGrid*540+330] = 4.0E0*I_ESP_I4x2y_Dyz_aa-2.0E0*1*I_ESP_G4x_Dyz_a;
    abcd[iGrid*540+331] = 4.0E0*I_ESP_I3x3y_Dyz_aa-2.0E0*1*I_ESP_G3xy_Dyz_a-2.0E0*2*I_ESP_G3xy_Dyz_a;
    abcd[iGrid*540+332] = 4.0E0*I_ESP_I3x2yz_Dyz_aa-2.0E0*1*I_ESP_G3xz_Dyz_a;
    abcd[iGrid*540+333] = 4.0E0*I_ESP_I2x4y_Dyz_aa-2.0E0*2*I_ESP_G2x2y_Dyz_a-2.0E0*3*I_ESP_G2x2y_Dyz_a+2*1*I_ESP_D2x_Dyz;
    abcd[iGrid*540+334] = 4.0E0*I_ESP_I2x3yz_Dyz_aa-2.0E0*1*I_ESP_G2xyz_Dyz_a-2.0E0*2*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+335] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2x2z_Dyz_a;
    abcd[iGrid*540+336] = 4.0E0*I_ESP_Ix5y_Dyz_aa-2.0E0*3*I_ESP_Gx3y_Dyz_a-2.0E0*4*I_ESP_Gx3y_Dyz_a+3*2*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+337] = 4.0E0*I_ESP_Ix4yz_Dyz_aa-2.0E0*2*I_ESP_Gx2yz_Dyz_a-2.0E0*3*I_ESP_Gx2yz_Dyz_a+2*1*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+338] = 4.0E0*I_ESP_Ix3y2z_Dyz_aa-2.0E0*1*I_ESP_Gxy2z_Dyz_a-2.0E0*2*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+339] = 4.0E0*I_ESP_Ix2y3z_Dyz_aa-2.0E0*1*I_ESP_Gx3z_Dyz_a;
    abcd[iGrid*540+340] = 4.0E0*I_ESP_I6y_Dyz_aa-2.0E0*4*I_ESP_G4y_Dyz_a-2.0E0*5*I_ESP_G4y_Dyz_a+4*3*I_ESP_D2y_Dyz;
    abcd[iGrid*540+341] = 4.0E0*I_ESP_I5yz_Dyz_aa-2.0E0*3*I_ESP_G3yz_Dyz_a-2.0E0*4*I_ESP_G3yz_Dyz_a+3*2*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+342] = 4.0E0*I_ESP_I4y2z_Dyz_aa-2.0E0*2*I_ESP_G2y2z_Dyz_a-2.0E0*3*I_ESP_G2y2z_Dyz_a+2*1*I_ESP_D2z_Dyz;
    abcd[iGrid*540+343] = 4.0E0*I_ESP_I3y3z_Dyz_aa-2.0E0*1*I_ESP_Gy3z_Dyz_a-2.0E0*2*I_ESP_Gy3z_Dyz_a;
    abcd[iGrid*540+344] = 4.0E0*I_ESP_I2y4z_Dyz_aa-2.0E0*1*I_ESP_G4z_Dyz_a;
    abcd[iGrid*540+345] = 4.0E0*I_ESP_I4x2y_D2z_aa-2.0E0*1*I_ESP_G4x_D2z_a;
    abcd[iGrid*540+346] = 4.0E0*I_ESP_I3x3y_D2z_aa-2.0E0*1*I_ESP_G3xy_D2z_a-2.0E0*2*I_ESP_G3xy_D2z_a;
    abcd[iGrid*540+347] = 4.0E0*I_ESP_I3x2yz_D2z_aa-2.0E0*1*I_ESP_G3xz_D2z_a;
    abcd[iGrid*540+348] = 4.0E0*I_ESP_I2x4y_D2z_aa-2.0E0*2*I_ESP_G2x2y_D2z_a-2.0E0*3*I_ESP_G2x2y_D2z_a+2*1*I_ESP_D2x_D2z;
    abcd[iGrid*540+349] = 4.0E0*I_ESP_I2x3yz_D2z_aa-2.0E0*1*I_ESP_G2xyz_D2z_a-2.0E0*2*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+350] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2x2z_D2z_a;
    abcd[iGrid*540+351] = 4.0E0*I_ESP_Ix5y_D2z_aa-2.0E0*3*I_ESP_Gx3y_D2z_a-2.0E0*4*I_ESP_Gx3y_D2z_a+3*2*I_ESP_Dxy_D2z;
    abcd[iGrid*540+352] = 4.0E0*I_ESP_Ix4yz_D2z_aa-2.0E0*2*I_ESP_Gx2yz_D2z_a-2.0E0*3*I_ESP_Gx2yz_D2z_a+2*1*I_ESP_Dxz_D2z;
    abcd[iGrid*540+353] = 4.0E0*I_ESP_Ix3y2z_D2z_aa-2.0E0*1*I_ESP_Gxy2z_D2z_a-2.0E0*2*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+354] = 4.0E0*I_ESP_Ix2y3z_D2z_aa-2.0E0*1*I_ESP_Gx3z_D2z_a;
    abcd[iGrid*540+355] = 4.0E0*I_ESP_I6y_D2z_aa-2.0E0*4*I_ESP_G4y_D2z_a-2.0E0*5*I_ESP_G4y_D2z_a+4*3*I_ESP_D2y_D2z;
    abcd[iGrid*540+356] = 4.0E0*I_ESP_I5yz_D2z_aa-2.0E0*3*I_ESP_G3yz_D2z_a-2.0E0*4*I_ESP_G3yz_D2z_a+3*2*I_ESP_Dyz_D2z;
    abcd[iGrid*540+357] = 4.0E0*I_ESP_I4y2z_D2z_aa-2.0E0*2*I_ESP_G2y2z_D2z_a-2.0E0*3*I_ESP_G2y2z_D2z_a+2*1*I_ESP_D2z_D2z;
    abcd[iGrid*540+358] = 4.0E0*I_ESP_I3y3z_D2z_aa-2.0E0*1*I_ESP_Gy3z_D2z_a-2.0E0*2*I_ESP_Gy3z_D2z_a;
    abcd[iGrid*540+359] = 4.0E0*I_ESP_I2y4z_D2z_aa-2.0E0*1*I_ESP_G4z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+360] = 4.0E0*I_ESP_I4xyz_D2x_aa;
    abcd[iGrid*540+361] = 4.0E0*I_ESP_I3x2yz_D2x_aa-2.0E0*1*I_ESP_G3xz_D2x_a;
    abcd[iGrid*540+362] = 4.0E0*I_ESP_I3xy2z_D2x_aa-2.0E0*1*I_ESP_G3xy_D2x_a;
    abcd[iGrid*540+363] = 4.0E0*I_ESP_I2x3yz_D2x_aa-2.0E0*2*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+364] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2x2y_D2x_a-2.0E0*1*I_ESP_G2x2z_D2x_a+1*I_ESP_D2x_D2x;
    abcd[iGrid*540+365] = 4.0E0*I_ESP_I2xy3z_D2x_aa-2.0E0*2*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+366] = 4.0E0*I_ESP_Ix4yz_D2x_aa-2.0E0*3*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+367] = 4.0E0*I_ESP_Ix3y2z_D2x_aa-2.0E0*1*I_ESP_Gx3y_D2x_a-2.0E0*2*I_ESP_Gxy2z_D2x_a+2*1*I_ESP_Dxy_D2x;
    abcd[iGrid*540+368] = 4.0E0*I_ESP_Ix2y3z_D2x_aa-2.0E0*2*I_ESP_Gx2yz_D2x_a-2.0E0*1*I_ESP_Gx3z_D2x_a+2*I_ESP_Dxz_D2x;
    abcd[iGrid*540+369] = 4.0E0*I_ESP_Ixy4z_D2x_aa-2.0E0*3*I_ESP_Gxy2z_D2x_a;
    abcd[iGrid*540+370] = 4.0E0*I_ESP_I5yz_D2x_aa-2.0E0*4*I_ESP_G3yz_D2x_a;
    abcd[iGrid*540+371] = 4.0E0*I_ESP_I4y2z_D2x_aa-2.0E0*1*I_ESP_G4y_D2x_a-2.0E0*3*I_ESP_G2y2z_D2x_a+3*1*I_ESP_D2y_D2x;
    abcd[iGrid*540+372] = 4.0E0*I_ESP_I3y3z_D2x_aa-2.0E0*2*I_ESP_G3yz_D2x_a-2.0E0*2*I_ESP_Gy3z_D2x_a+2*2*I_ESP_Dyz_D2x;
    abcd[iGrid*540+373] = 4.0E0*I_ESP_I2y4z_D2x_aa-2.0E0*3*I_ESP_G2y2z_D2x_a-2.0E0*1*I_ESP_G4z_D2x_a+3*I_ESP_D2z_D2x;
    abcd[iGrid*540+374] = 4.0E0*I_ESP_Iy5z_D2x_aa-2.0E0*4*I_ESP_Gy3z_D2x_a;
    abcd[iGrid*540+375] = 4.0E0*I_ESP_I4xyz_Dxy_aa;
    abcd[iGrid*540+376] = 4.0E0*I_ESP_I3x2yz_Dxy_aa-2.0E0*1*I_ESP_G3xz_Dxy_a;
    abcd[iGrid*540+377] = 4.0E0*I_ESP_I3xy2z_Dxy_aa-2.0E0*1*I_ESP_G3xy_Dxy_a;
    abcd[iGrid*540+378] = 4.0E0*I_ESP_I2x3yz_Dxy_aa-2.0E0*2*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+379] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2x2y_Dxy_a-2.0E0*1*I_ESP_G2x2z_Dxy_a+1*I_ESP_D2x_Dxy;
    abcd[iGrid*540+380] = 4.0E0*I_ESP_I2xy3z_Dxy_aa-2.0E0*2*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+381] = 4.0E0*I_ESP_Ix4yz_Dxy_aa-2.0E0*3*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+382] = 4.0E0*I_ESP_Ix3y2z_Dxy_aa-2.0E0*1*I_ESP_Gx3y_Dxy_a-2.0E0*2*I_ESP_Gxy2z_Dxy_a+2*1*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+383] = 4.0E0*I_ESP_Ix2y3z_Dxy_aa-2.0E0*2*I_ESP_Gx2yz_Dxy_a-2.0E0*1*I_ESP_Gx3z_Dxy_a+2*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+384] = 4.0E0*I_ESP_Ixy4z_Dxy_aa-2.0E0*3*I_ESP_Gxy2z_Dxy_a;
    abcd[iGrid*540+385] = 4.0E0*I_ESP_I5yz_Dxy_aa-2.0E0*4*I_ESP_G3yz_Dxy_a;
    abcd[iGrid*540+386] = 4.0E0*I_ESP_I4y2z_Dxy_aa-2.0E0*1*I_ESP_G4y_Dxy_a-2.0E0*3*I_ESP_G2y2z_Dxy_a+3*1*I_ESP_D2y_Dxy;
    abcd[iGrid*540+387] = 4.0E0*I_ESP_I3y3z_Dxy_aa-2.0E0*2*I_ESP_G3yz_Dxy_a-2.0E0*2*I_ESP_Gy3z_Dxy_a+2*2*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+388] = 4.0E0*I_ESP_I2y4z_Dxy_aa-2.0E0*3*I_ESP_G2y2z_Dxy_a-2.0E0*1*I_ESP_G4z_Dxy_a+3*I_ESP_D2z_Dxy;
    abcd[iGrid*540+389] = 4.0E0*I_ESP_Iy5z_Dxy_aa-2.0E0*4*I_ESP_Gy3z_Dxy_a;
    abcd[iGrid*540+390] = 4.0E0*I_ESP_I4xyz_Dxz_aa;
    abcd[iGrid*540+391] = 4.0E0*I_ESP_I3x2yz_Dxz_aa-2.0E0*1*I_ESP_G3xz_Dxz_a;
    abcd[iGrid*540+392] = 4.0E0*I_ESP_I3xy2z_Dxz_aa-2.0E0*1*I_ESP_G3xy_Dxz_a;
    abcd[iGrid*540+393] = 4.0E0*I_ESP_I2x3yz_Dxz_aa-2.0E0*2*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+394] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2x2y_Dxz_a-2.0E0*1*I_ESP_G2x2z_Dxz_a+1*I_ESP_D2x_Dxz;
    abcd[iGrid*540+395] = 4.0E0*I_ESP_I2xy3z_Dxz_aa-2.0E0*2*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+396] = 4.0E0*I_ESP_Ix4yz_Dxz_aa-2.0E0*3*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+397] = 4.0E0*I_ESP_Ix3y2z_Dxz_aa-2.0E0*1*I_ESP_Gx3y_Dxz_a-2.0E0*2*I_ESP_Gxy2z_Dxz_a+2*1*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+398] = 4.0E0*I_ESP_Ix2y3z_Dxz_aa-2.0E0*2*I_ESP_Gx2yz_Dxz_a-2.0E0*1*I_ESP_Gx3z_Dxz_a+2*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+399] = 4.0E0*I_ESP_Ixy4z_Dxz_aa-2.0E0*3*I_ESP_Gxy2z_Dxz_a;
    abcd[iGrid*540+400] = 4.0E0*I_ESP_I5yz_Dxz_aa-2.0E0*4*I_ESP_G3yz_Dxz_a;
    abcd[iGrid*540+401] = 4.0E0*I_ESP_I4y2z_Dxz_aa-2.0E0*1*I_ESP_G4y_Dxz_a-2.0E0*3*I_ESP_G2y2z_Dxz_a+3*1*I_ESP_D2y_Dxz;
    abcd[iGrid*540+402] = 4.0E0*I_ESP_I3y3z_Dxz_aa-2.0E0*2*I_ESP_G3yz_Dxz_a-2.0E0*2*I_ESP_Gy3z_Dxz_a+2*2*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+403] = 4.0E0*I_ESP_I2y4z_Dxz_aa-2.0E0*3*I_ESP_G2y2z_Dxz_a-2.0E0*1*I_ESP_G4z_Dxz_a+3*I_ESP_D2z_Dxz;
    abcd[iGrid*540+404] = 4.0E0*I_ESP_Iy5z_Dxz_aa-2.0E0*4*I_ESP_Gy3z_Dxz_a;
    abcd[iGrid*540+405] = 4.0E0*I_ESP_I4xyz_D2y_aa;
    abcd[iGrid*540+406] = 4.0E0*I_ESP_I3x2yz_D2y_aa-2.0E0*1*I_ESP_G3xz_D2y_a;
    abcd[iGrid*540+407] = 4.0E0*I_ESP_I3xy2z_D2y_aa-2.0E0*1*I_ESP_G3xy_D2y_a;
    abcd[iGrid*540+408] = 4.0E0*I_ESP_I2x3yz_D2y_aa-2.0E0*2*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+409] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2x2y_D2y_a-2.0E0*1*I_ESP_G2x2z_D2y_a+1*I_ESP_D2x_D2y;
    abcd[iGrid*540+410] = 4.0E0*I_ESP_I2xy3z_D2y_aa-2.0E0*2*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+411] = 4.0E0*I_ESP_Ix4yz_D2y_aa-2.0E0*3*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+412] = 4.0E0*I_ESP_Ix3y2z_D2y_aa-2.0E0*1*I_ESP_Gx3y_D2y_a-2.0E0*2*I_ESP_Gxy2z_D2y_a+2*1*I_ESP_Dxy_D2y;
    abcd[iGrid*540+413] = 4.0E0*I_ESP_Ix2y3z_D2y_aa-2.0E0*2*I_ESP_Gx2yz_D2y_a-2.0E0*1*I_ESP_Gx3z_D2y_a+2*I_ESP_Dxz_D2y;
    abcd[iGrid*540+414] = 4.0E0*I_ESP_Ixy4z_D2y_aa-2.0E0*3*I_ESP_Gxy2z_D2y_a;
    abcd[iGrid*540+415] = 4.0E0*I_ESP_I5yz_D2y_aa-2.0E0*4*I_ESP_G3yz_D2y_a;
    abcd[iGrid*540+416] = 4.0E0*I_ESP_I4y2z_D2y_aa-2.0E0*1*I_ESP_G4y_D2y_a-2.0E0*3*I_ESP_G2y2z_D2y_a+3*1*I_ESP_D2y_D2y;
    abcd[iGrid*540+417] = 4.0E0*I_ESP_I3y3z_D2y_aa-2.0E0*2*I_ESP_G3yz_D2y_a-2.0E0*2*I_ESP_Gy3z_D2y_a+2*2*I_ESP_Dyz_D2y;
    abcd[iGrid*540+418] = 4.0E0*I_ESP_I2y4z_D2y_aa-2.0E0*3*I_ESP_G2y2z_D2y_a-2.0E0*1*I_ESP_G4z_D2y_a+3*I_ESP_D2z_D2y;
    abcd[iGrid*540+419] = 4.0E0*I_ESP_Iy5z_D2y_aa-2.0E0*4*I_ESP_Gy3z_D2y_a;
    abcd[iGrid*540+420] = 4.0E0*I_ESP_I4xyz_Dyz_aa;
    abcd[iGrid*540+421] = 4.0E0*I_ESP_I3x2yz_Dyz_aa-2.0E0*1*I_ESP_G3xz_Dyz_a;
    abcd[iGrid*540+422] = 4.0E0*I_ESP_I3xy2z_Dyz_aa-2.0E0*1*I_ESP_G3xy_Dyz_a;
    abcd[iGrid*540+423] = 4.0E0*I_ESP_I2x3yz_Dyz_aa-2.0E0*2*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+424] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2x2y_Dyz_a-2.0E0*1*I_ESP_G2x2z_Dyz_a+1*I_ESP_D2x_Dyz;
    abcd[iGrid*540+425] = 4.0E0*I_ESP_I2xy3z_Dyz_aa-2.0E0*2*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+426] = 4.0E0*I_ESP_Ix4yz_Dyz_aa-2.0E0*3*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+427] = 4.0E0*I_ESP_Ix3y2z_Dyz_aa-2.0E0*1*I_ESP_Gx3y_Dyz_a-2.0E0*2*I_ESP_Gxy2z_Dyz_a+2*1*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+428] = 4.0E0*I_ESP_Ix2y3z_Dyz_aa-2.0E0*2*I_ESP_Gx2yz_Dyz_a-2.0E0*1*I_ESP_Gx3z_Dyz_a+2*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+429] = 4.0E0*I_ESP_Ixy4z_Dyz_aa-2.0E0*3*I_ESP_Gxy2z_Dyz_a;
    abcd[iGrid*540+430] = 4.0E0*I_ESP_I5yz_Dyz_aa-2.0E0*4*I_ESP_G3yz_Dyz_a;
    abcd[iGrid*540+431] = 4.0E0*I_ESP_I4y2z_Dyz_aa-2.0E0*1*I_ESP_G4y_Dyz_a-2.0E0*3*I_ESP_G2y2z_Dyz_a+3*1*I_ESP_D2y_Dyz;
    abcd[iGrid*540+432] = 4.0E0*I_ESP_I3y3z_Dyz_aa-2.0E0*2*I_ESP_G3yz_Dyz_a-2.0E0*2*I_ESP_Gy3z_Dyz_a+2*2*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+433] = 4.0E0*I_ESP_I2y4z_Dyz_aa-2.0E0*3*I_ESP_G2y2z_Dyz_a-2.0E0*1*I_ESP_G4z_Dyz_a+3*I_ESP_D2z_Dyz;
    abcd[iGrid*540+434] = 4.0E0*I_ESP_Iy5z_Dyz_aa-2.0E0*4*I_ESP_Gy3z_Dyz_a;
    abcd[iGrid*540+435] = 4.0E0*I_ESP_I4xyz_D2z_aa;
    abcd[iGrid*540+436] = 4.0E0*I_ESP_I3x2yz_D2z_aa-2.0E0*1*I_ESP_G3xz_D2z_a;
    abcd[iGrid*540+437] = 4.0E0*I_ESP_I3xy2z_D2z_aa-2.0E0*1*I_ESP_G3xy_D2z_a;
    abcd[iGrid*540+438] = 4.0E0*I_ESP_I2x3yz_D2z_aa-2.0E0*2*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+439] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2x2y_D2z_a-2.0E0*1*I_ESP_G2x2z_D2z_a+1*I_ESP_D2x_D2z;
    abcd[iGrid*540+440] = 4.0E0*I_ESP_I2xy3z_D2z_aa-2.0E0*2*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+441] = 4.0E0*I_ESP_Ix4yz_D2z_aa-2.0E0*3*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+442] = 4.0E0*I_ESP_Ix3y2z_D2z_aa-2.0E0*1*I_ESP_Gx3y_D2z_a-2.0E0*2*I_ESP_Gxy2z_D2z_a+2*1*I_ESP_Dxy_D2z;
    abcd[iGrid*540+443] = 4.0E0*I_ESP_Ix2y3z_D2z_aa-2.0E0*2*I_ESP_Gx2yz_D2z_a-2.0E0*1*I_ESP_Gx3z_D2z_a+2*I_ESP_Dxz_D2z;
    abcd[iGrid*540+444] = 4.0E0*I_ESP_Ixy4z_D2z_aa-2.0E0*3*I_ESP_Gxy2z_D2z_a;
    abcd[iGrid*540+445] = 4.0E0*I_ESP_I5yz_D2z_aa-2.0E0*4*I_ESP_G3yz_D2z_a;
    abcd[iGrid*540+446] = 4.0E0*I_ESP_I4y2z_D2z_aa-2.0E0*1*I_ESP_G4y_D2z_a-2.0E0*3*I_ESP_G2y2z_D2z_a+3*1*I_ESP_D2y_D2z;
    abcd[iGrid*540+447] = 4.0E0*I_ESP_I3y3z_D2z_aa-2.0E0*2*I_ESP_G3yz_D2z_a-2.0E0*2*I_ESP_Gy3z_D2z_a+2*2*I_ESP_Dyz_D2z;
    abcd[iGrid*540+448] = 4.0E0*I_ESP_I2y4z_D2z_aa-2.0E0*3*I_ESP_G2y2z_D2z_a-2.0E0*1*I_ESP_G4z_D2z_a+3*I_ESP_D2z_D2z;
    abcd[iGrid*540+449] = 4.0E0*I_ESP_Iy5z_D2z_aa-2.0E0*4*I_ESP_Gy3z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    abcd[iGrid*540+450] = 4.0E0*I_ESP_I4x2z_D2x_aa-2.0E0*1*I_ESP_G4x_D2x_a;
    abcd[iGrid*540+451] = 4.0E0*I_ESP_I3xy2z_D2x_aa-2.0E0*1*I_ESP_G3xy_D2x_a;
    abcd[iGrid*540+452] = 4.0E0*I_ESP_I3x3z_D2x_aa-2.0E0*1*I_ESP_G3xz_D2x_a-2.0E0*2*I_ESP_G3xz_D2x_a;
    abcd[iGrid*540+453] = 4.0E0*I_ESP_I2x2y2z_D2x_aa-2.0E0*1*I_ESP_G2x2y_D2x_a;
    abcd[iGrid*540+454] = 4.0E0*I_ESP_I2xy3z_D2x_aa-2.0E0*1*I_ESP_G2xyz_D2x_a-2.0E0*2*I_ESP_G2xyz_D2x_a;
    abcd[iGrid*540+455] = 4.0E0*I_ESP_I2x4z_D2x_aa-2.0E0*2*I_ESP_G2x2z_D2x_a-2.0E0*3*I_ESP_G2x2z_D2x_a+2*1*I_ESP_D2x_D2x;
    abcd[iGrid*540+456] = 4.0E0*I_ESP_Ix3y2z_D2x_aa-2.0E0*1*I_ESP_Gx3y_D2x_a;
    abcd[iGrid*540+457] = 4.0E0*I_ESP_Ix2y3z_D2x_aa-2.0E0*1*I_ESP_Gx2yz_D2x_a-2.0E0*2*I_ESP_Gx2yz_D2x_a;
    abcd[iGrid*540+458] = 4.0E0*I_ESP_Ixy4z_D2x_aa-2.0E0*2*I_ESP_Gxy2z_D2x_a-2.0E0*3*I_ESP_Gxy2z_D2x_a+2*1*I_ESP_Dxy_D2x;
    abcd[iGrid*540+459] = 4.0E0*I_ESP_Ix5z_D2x_aa-2.0E0*3*I_ESP_Gx3z_D2x_a-2.0E0*4*I_ESP_Gx3z_D2x_a+3*2*I_ESP_Dxz_D2x;
    abcd[iGrid*540+460] = 4.0E0*I_ESP_I4y2z_D2x_aa-2.0E0*1*I_ESP_G4y_D2x_a;
    abcd[iGrid*540+461] = 4.0E0*I_ESP_I3y3z_D2x_aa-2.0E0*1*I_ESP_G3yz_D2x_a-2.0E0*2*I_ESP_G3yz_D2x_a;
    abcd[iGrid*540+462] = 4.0E0*I_ESP_I2y4z_D2x_aa-2.0E0*2*I_ESP_G2y2z_D2x_a-2.0E0*3*I_ESP_G2y2z_D2x_a+2*1*I_ESP_D2y_D2x;
    abcd[iGrid*540+463] = 4.0E0*I_ESP_Iy5z_D2x_aa-2.0E0*3*I_ESP_Gy3z_D2x_a-2.0E0*4*I_ESP_Gy3z_D2x_a+3*2*I_ESP_Dyz_D2x;
    abcd[iGrid*540+464] = 4.0E0*I_ESP_I6z_D2x_aa-2.0E0*4*I_ESP_G4z_D2x_a-2.0E0*5*I_ESP_G4z_D2x_a+4*3*I_ESP_D2z_D2x;
    abcd[iGrid*540+465] = 4.0E0*I_ESP_I4x2z_Dxy_aa-2.0E0*1*I_ESP_G4x_Dxy_a;
    abcd[iGrid*540+466] = 4.0E0*I_ESP_I3xy2z_Dxy_aa-2.0E0*1*I_ESP_G3xy_Dxy_a;
    abcd[iGrid*540+467] = 4.0E0*I_ESP_I3x3z_Dxy_aa-2.0E0*1*I_ESP_G3xz_Dxy_a-2.0E0*2*I_ESP_G3xz_Dxy_a;
    abcd[iGrid*540+468] = 4.0E0*I_ESP_I2x2y2z_Dxy_aa-2.0E0*1*I_ESP_G2x2y_Dxy_a;
    abcd[iGrid*540+469] = 4.0E0*I_ESP_I2xy3z_Dxy_aa-2.0E0*1*I_ESP_G2xyz_Dxy_a-2.0E0*2*I_ESP_G2xyz_Dxy_a;
    abcd[iGrid*540+470] = 4.0E0*I_ESP_I2x4z_Dxy_aa-2.0E0*2*I_ESP_G2x2z_Dxy_a-2.0E0*3*I_ESP_G2x2z_Dxy_a+2*1*I_ESP_D2x_Dxy;
    abcd[iGrid*540+471] = 4.0E0*I_ESP_Ix3y2z_Dxy_aa-2.0E0*1*I_ESP_Gx3y_Dxy_a;
    abcd[iGrid*540+472] = 4.0E0*I_ESP_Ix2y3z_Dxy_aa-2.0E0*1*I_ESP_Gx2yz_Dxy_a-2.0E0*2*I_ESP_Gx2yz_Dxy_a;
    abcd[iGrid*540+473] = 4.0E0*I_ESP_Ixy4z_Dxy_aa-2.0E0*2*I_ESP_Gxy2z_Dxy_a-2.0E0*3*I_ESP_Gxy2z_Dxy_a+2*1*I_ESP_Dxy_Dxy;
    abcd[iGrid*540+474] = 4.0E0*I_ESP_Ix5z_Dxy_aa-2.0E0*3*I_ESP_Gx3z_Dxy_a-2.0E0*4*I_ESP_Gx3z_Dxy_a+3*2*I_ESP_Dxz_Dxy;
    abcd[iGrid*540+475] = 4.0E0*I_ESP_I4y2z_Dxy_aa-2.0E0*1*I_ESP_G4y_Dxy_a;
    abcd[iGrid*540+476] = 4.0E0*I_ESP_I3y3z_Dxy_aa-2.0E0*1*I_ESP_G3yz_Dxy_a-2.0E0*2*I_ESP_G3yz_Dxy_a;
    abcd[iGrid*540+477] = 4.0E0*I_ESP_I2y4z_Dxy_aa-2.0E0*2*I_ESP_G2y2z_Dxy_a-2.0E0*3*I_ESP_G2y2z_Dxy_a+2*1*I_ESP_D2y_Dxy;
    abcd[iGrid*540+478] = 4.0E0*I_ESP_Iy5z_Dxy_aa-2.0E0*3*I_ESP_Gy3z_Dxy_a-2.0E0*4*I_ESP_Gy3z_Dxy_a+3*2*I_ESP_Dyz_Dxy;
    abcd[iGrid*540+479] = 4.0E0*I_ESP_I6z_Dxy_aa-2.0E0*4*I_ESP_G4z_Dxy_a-2.0E0*5*I_ESP_G4z_Dxy_a+4*3*I_ESP_D2z_Dxy;
    abcd[iGrid*540+480] = 4.0E0*I_ESP_I4x2z_Dxz_aa-2.0E0*1*I_ESP_G4x_Dxz_a;
    abcd[iGrid*540+481] = 4.0E0*I_ESP_I3xy2z_Dxz_aa-2.0E0*1*I_ESP_G3xy_Dxz_a;
    abcd[iGrid*540+482] = 4.0E0*I_ESP_I3x3z_Dxz_aa-2.0E0*1*I_ESP_G3xz_Dxz_a-2.0E0*2*I_ESP_G3xz_Dxz_a;
    abcd[iGrid*540+483] = 4.0E0*I_ESP_I2x2y2z_Dxz_aa-2.0E0*1*I_ESP_G2x2y_Dxz_a;
    abcd[iGrid*540+484] = 4.0E0*I_ESP_I2xy3z_Dxz_aa-2.0E0*1*I_ESP_G2xyz_Dxz_a-2.0E0*2*I_ESP_G2xyz_Dxz_a;
    abcd[iGrid*540+485] = 4.0E0*I_ESP_I2x4z_Dxz_aa-2.0E0*2*I_ESP_G2x2z_Dxz_a-2.0E0*3*I_ESP_G2x2z_Dxz_a+2*1*I_ESP_D2x_Dxz;
    abcd[iGrid*540+486] = 4.0E0*I_ESP_Ix3y2z_Dxz_aa-2.0E0*1*I_ESP_Gx3y_Dxz_a;
    abcd[iGrid*540+487] = 4.0E0*I_ESP_Ix2y3z_Dxz_aa-2.0E0*1*I_ESP_Gx2yz_Dxz_a-2.0E0*2*I_ESP_Gx2yz_Dxz_a;
    abcd[iGrid*540+488] = 4.0E0*I_ESP_Ixy4z_Dxz_aa-2.0E0*2*I_ESP_Gxy2z_Dxz_a-2.0E0*3*I_ESP_Gxy2z_Dxz_a+2*1*I_ESP_Dxy_Dxz;
    abcd[iGrid*540+489] = 4.0E0*I_ESP_Ix5z_Dxz_aa-2.0E0*3*I_ESP_Gx3z_Dxz_a-2.0E0*4*I_ESP_Gx3z_Dxz_a+3*2*I_ESP_Dxz_Dxz;
    abcd[iGrid*540+490] = 4.0E0*I_ESP_I4y2z_Dxz_aa-2.0E0*1*I_ESP_G4y_Dxz_a;
    abcd[iGrid*540+491] = 4.0E0*I_ESP_I3y3z_Dxz_aa-2.0E0*1*I_ESP_G3yz_Dxz_a-2.0E0*2*I_ESP_G3yz_Dxz_a;
    abcd[iGrid*540+492] = 4.0E0*I_ESP_I2y4z_Dxz_aa-2.0E0*2*I_ESP_G2y2z_Dxz_a-2.0E0*3*I_ESP_G2y2z_Dxz_a+2*1*I_ESP_D2y_Dxz;
    abcd[iGrid*540+493] = 4.0E0*I_ESP_Iy5z_Dxz_aa-2.0E0*3*I_ESP_Gy3z_Dxz_a-2.0E0*4*I_ESP_Gy3z_Dxz_a+3*2*I_ESP_Dyz_Dxz;
    abcd[iGrid*540+494] = 4.0E0*I_ESP_I6z_Dxz_aa-2.0E0*4*I_ESP_G4z_Dxz_a-2.0E0*5*I_ESP_G4z_Dxz_a+4*3*I_ESP_D2z_Dxz;
    abcd[iGrid*540+495] = 4.0E0*I_ESP_I4x2z_D2y_aa-2.0E0*1*I_ESP_G4x_D2y_a;
    abcd[iGrid*540+496] = 4.0E0*I_ESP_I3xy2z_D2y_aa-2.0E0*1*I_ESP_G3xy_D2y_a;
    abcd[iGrid*540+497] = 4.0E0*I_ESP_I3x3z_D2y_aa-2.0E0*1*I_ESP_G3xz_D2y_a-2.0E0*2*I_ESP_G3xz_D2y_a;
    abcd[iGrid*540+498] = 4.0E0*I_ESP_I2x2y2z_D2y_aa-2.0E0*1*I_ESP_G2x2y_D2y_a;
    abcd[iGrid*540+499] = 4.0E0*I_ESP_I2xy3z_D2y_aa-2.0E0*1*I_ESP_G2xyz_D2y_a-2.0E0*2*I_ESP_G2xyz_D2y_a;
    abcd[iGrid*540+500] = 4.0E0*I_ESP_I2x4z_D2y_aa-2.0E0*2*I_ESP_G2x2z_D2y_a-2.0E0*3*I_ESP_G2x2z_D2y_a+2*1*I_ESP_D2x_D2y;
    abcd[iGrid*540+501] = 4.0E0*I_ESP_Ix3y2z_D2y_aa-2.0E0*1*I_ESP_Gx3y_D2y_a;
    abcd[iGrid*540+502] = 4.0E0*I_ESP_Ix2y3z_D2y_aa-2.0E0*1*I_ESP_Gx2yz_D2y_a-2.0E0*2*I_ESP_Gx2yz_D2y_a;
    abcd[iGrid*540+503] = 4.0E0*I_ESP_Ixy4z_D2y_aa-2.0E0*2*I_ESP_Gxy2z_D2y_a-2.0E0*3*I_ESP_Gxy2z_D2y_a+2*1*I_ESP_Dxy_D2y;
    abcd[iGrid*540+504] = 4.0E0*I_ESP_Ix5z_D2y_aa-2.0E0*3*I_ESP_Gx3z_D2y_a-2.0E0*4*I_ESP_Gx3z_D2y_a+3*2*I_ESP_Dxz_D2y;
    abcd[iGrid*540+505] = 4.0E0*I_ESP_I4y2z_D2y_aa-2.0E0*1*I_ESP_G4y_D2y_a;
    abcd[iGrid*540+506] = 4.0E0*I_ESP_I3y3z_D2y_aa-2.0E0*1*I_ESP_G3yz_D2y_a-2.0E0*2*I_ESP_G3yz_D2y_a;
    abcd[iGrid*540+507] = 4.0E0*I_ESP_I2y4z_D2y_aa-2.0E0*2*I_ESP_G2y2z_D2y_a-2.0E0*3*I_ESP_G2y2z_D2y_a+2*1*I_ESP_D2y_D2y;
    abcd[iGrid*540+508] = 4.0E0*I_ESP_Iy5z_D2y_aa-2.0E0*3*I_ESP_Gy3z_D2y_a-2.0E0*4*I_ESP_Gy3z_D2y_a+3*2*I_ESP_Dyz_D2y;
    abcd[iGrid*540+509] = 4.0E0*I_ESP_I6z_D2y_aa-2.0E0*4*I_ESP_G4z_D2y_a-2.0E0*5*I_ESP_G4z_D2y_a+4*3*I_ESP_D2z_D2y;
    abcd[iGrid*540+510] = 4.0E0*I_ESP_I4x2z_Dyz_aa-2.0E0*1*I_ESP_G4x_Dyz_a;
    abcd[iGrid*540+511] = 4.0E0*I_ESP_I3xy2z_Dyz_aa-2.0E0*1*I_ESP_G3xy_Dyz_a;
    abcd[iGrid*540+512] = 4.0E0*I_ESP_I3x3z_Dyz_aa-2.0E0*1*I_ESP_G3xz_Dyz_a-2.0E0*2*I_ESP_G3xz_Dyz_a;
    abcd[iGrid*540+513] = 4.0E0*I_ESP_I2x2y2z_Dyz_aa-2.0E0*1*I_ESP_G2x2y_Dyz_a;
    abcd[iGrid*540+514] = 4.0E0*I_ESP_I2xy3z_Dyz_aa-2.0E0*1*I_ESP_G2xyz_Dyz_a-2.0E0*2*I_ESP_G2xyz_Dyz_a;
    abcd[iGrid*540+515] = 4.0E0*I_ESP_I2x4z_Dyz_aa-2.0E0*2*I_ESP_G2x2z_Dyz_a-2.0E0*3*I_ESP_G2x2z_Dyz_a+2*1*I_ESP_D2x_Dyz;
    abcd[iGrid*540+516] = 4.0E0*I_ESP_Ix3y2z_Dyz_aa-2.0E0*1*I_ESP_Gx3y_Dyz_a;
    abcd[iGrid*540+517] = 4.0E0*I_ESP_Ix2y3z_Dyz_aa-2.0E0*1*I_ESP_Gx2yz_Dyz_a-2.0E0*2*I_ESP_Gx2yz_Dyz_a;
    abcd[iGrid*540+518] = 4.0E0*I_ESP_Ixy4z_Dyz_aa-2.0E0*2*I_ESP_Gxy2z_Dyz_a-2.0E0*3*I_ESP_Gxy2z_Dyz_a+2*1*I_ESP_Dxy_Dyz;
    abcd[iGrid*540+519] = 4.0E0*I_ESP_Ix5z_Dyz_aa-2.0E0*3*I_ESP_Gx3z_Dyz_a-2.0E0*4*I_ESP_Gx3z_Dyz_a+3*2*I_ESP_Dxz_Dyz;
    abcd[iGrid*540+520] = 4.0E0*I_ESP_I4y2z_Dyz_aa-2.0E0*1*I_ESP_G4y_Dyz_a;
    abcd[iGrid*540+521] = 4.0E0*I_ESP_I3y3z_Dyz_aa-2.0E0*1*I_ESP_G3yz_Dyz_a-2.0E0*2*I_ESP_G3yz_Dyz_a;
    abcd[iGrid*540+522] = 4.0E0*I_ESP_I2y4z_Dyz_aa-2.0E0*2*I_ESP_G2y2z_Dyz_a-2.0E0*3*I_ESP_G2y2z_Dyz_a+2*1*I_ESP_D2y_Dyz;
    abcd[iGrid*540+523] = 4.0E0*I_ESP_Iy5z_Dyz_aa-2.0E0*3*I_ESP_Gy3z_Dyz_a-2.0E0*4*I_ESP_Gy3z_Dyz_a+3*2*I_ESP_Dyz_Dyz;
    abcd[iGrid*540+524] = 4.0E0*I_ESP_I6z_Dyz_aa-2.0E0*4*I_ESP_G4z_Dyz_a-2.0E0*5*I_ESP_G4z_Dyz_a+4*3*I_ESP_D2z_Dyz;
    abcd[iGrid*540+525] = 4.0E0*I_ESP_I4x2z_D2z_aa-2.0E0*1*I_ESP_G4x_D2z_a;
    abcd[iGrid*540+526] = 4.0E0*I_ESP_I3xy2z_D2z_aa-2.0E0*1*I_ESP_G3xy_D2z_a;
    abcd[iGrid*540+527] = 4.0E0*I_ESP_I3x3z_D2z_aa-2.0E0*1*I_ESP_G3xz_D2z_a-2.0E0*2*I_ESP_G3xz_D2z_a;
    abcd[iGrid*540+528] = 4.0E0*I_ESP_I2x2y2z_D2z_aa-2.0E0*1*I_ESP_G2x2y_D2z_a;
    abcd[iGrid*540+529] = 4.0E0*I_ESP_I2xy3z_D2z_aa-2.0E0*1*I_ESP_G2xyz_D2z_a-2.0E0*2*I_ESP_G2xyz_D2z_a;
    abcd[iGrid*540+530] = 4.0E0*I_ESP_I2x4z_D2z_aa-2.0E0*2*I_ESP_G2x2z_D2z_a-2.0E0*3*I_ESP_G2x2z_D2z_a+2*1*I_ESP_D2x_D2z;
    abcd[iGrid*540+531] = 4.0E0*I_ESP_Ix3y2z_D2z_aa-2.0E0*1*I_ESP_Gx3y_D2z_a;
    abcd[iGrid*540+532] = 4.0E0*I_ESP_Ix2y3z_D2z_aa-2.0E0*1*I_ESP_Gx2yz_D2z_a-2.0E0*2*I_ESP_Gx2yz_D2z_a;
    abcd[iGrid*540+533] = 4.0E0*I_ESP_Ixy4z_D2z_aa-2.0E0*2*I_ESP_Gxy2z_D2z_a-2.0E0*3*I_ESP_Gxy2z_D2z_a+2*1*I_ESP_Dxy_D2z;
    abcd[iGrid*540+534] = 4.0E0*I_ESP_Ix5z_D2z_aa-2.0E0*3*I_ESP_Gx3z_D2z_a-2.0E0*4*I_ESP_Gx3z_D2z_a+3*2*I_ESP_Dxz_D2z;
    abcd[iGrid*540+535] = 4.0E0*I_ESP_I4y2z_D2z_aa-2.0E0*1*I_ESP_G4y_D2z_a;
    abcd[iGrid*540+536] = 4.0E0*I_ESP_I3y3z_D2z_aa-2.0E0*1*I_ESP_G3yz_D2z_a-2.0E0*2*I_ESP_G3yz_D2z_a;
    abcd[iGrid*540+537] = 4.0E0*I_ESP_I2y4z_D2z_aa-2.0E0*2*I_ESP_G2y2z_D2z_a-2.0E0*3*I_ESP_G2y2z_D2z_a+2*1*I_ESP_D2y_D2z;
    abcd[iGrid*540+538] = 4.0E0*I_ESP_Iy5z_D2z_aa-2.0E0*3*I_ESP_Gy3z_D2z_a-2.0E0*4*I_ESP_Gy3z_D2z_a+3*2*I_ESP_Dyz_D2z;
    abcd[iGrid*540+539] = 4.0E0*I_ESP_I6z_D2z_aa-2.0E0*4*I_ESP_G4z_D2z_a-2.0E0*5*I_ESP_G4z_D2z_a+4*3*I_ESP_D2z_D2z;
  }
}
