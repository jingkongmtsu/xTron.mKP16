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
// BRA1 as redundant position, total RHS integrals evaluated as: 4102
// BRA2 as redundant position, total RHS integrals evaluated as: 3774
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####

void hgp_os_esp_g_f_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_L8x_S_a = 0.0E0;
    Double I_ESP_L7xy_S_a = 0.0E0;
    Double I_ESP_L7xz_S_a = 0.0E0;
    Double I_ESP_L6x2y_S_a = 0.0E0;
    Double I_ESP_L6xyz_S_a = 0.0E0;
    Double I_ESP_L6x2z_S_a = 0.0E0;
    Double I_ESP_L5x3y_S_a = 0.0E0;
    Double I_ESP_L5x2yz_S_a = 0.0E0;
    Double I_ESP_L5xy2z_S_a = 0.0E0;
    Double I_ESP_L5x3z_S_a = 0.0E0;
    Double I_ESP_L4x4y_S_a = 0.0E0;
    Double I_ESP_L4x3yz_S_a = 0.0E0;
    Double I_ESP_L4x2y2z_S_a = 0.0E0;
    Double I_ESP_L4xy3z_S_a = 0.0E0;
    Double I_ESP_L4x4z_S_a = 0.0E0;
    Double I_ESP_L3x5y_S_a = 0.0E0;
    Double I_ESP_L3x4yz_S_a = 0.0E0;
    Double I_ESP_L3x3y2z_S_a = 0.0E0;
    Double I_ESP_L3x2y3z_S_a = 0.0E0;
    Double I_ESP_L3xy4z_S_a = 0.0E0;
    Double I_ESP_L3x5z_S_a = 0.0E0;
    Double I_ESP_L2x6y_S_a = 0.0E0;
    Double I_ESP_L2x5yz_S_a = 0.0E0;
    Double I_ESP_L2x4y2z_S_a = 0.0E0;
    Double I_ESP_L2x3y3z_S_a = 0.0E0;
    Double I_ESP_L2x2y4z_S_a = 0.0E0;
    Double I_ESP_L2xy5z_S_a = 0.0E0;
    Double I_ESP_L2x6z_S_a = 0.0E0;
    Double I_ESP_Lx7y_S_a = 0.0E0;
    Double I_ESP_Lx6yz_S_a = 0.0E0;
    Double I_ESP_Lx5y2z_S_a = 0.0E0;
    Double I_ESP_Lx4y3z_S_a = 0.0E0;
    Double I_ESP_Lx3y4z_S_a = 0.0E0;
    Double I_ESP_Lx2y5z_S_a = 0.0E0;
    Double I_ESP_Lxy6z_S_a = 0.0E0;
    Double I_ESP_Lx7z_S_a = 0.0E0;
    Double I_ESP_L8y_S_a = 0.0E0;
    Double I_ESP_L7yz_S_a = 0.0E0;
    Double I_ESP_L6y2z_S_a = 0.0E0;
    Double I_ESP_L5y3z_S_a = 0.0E0;
    Double I_ESP_L4y4z_S_a = 0.0E0;
    Double I_ESP_L3y5z_S_a = 0.0E0;
    Double I_ESP_L2y6z_S_a = 0.0E0;
    Double I_ESP_Ly7z_S_a = 0.0E0;
    Double I_ESP_L8z_S_a = 0.0E0;
    Double I_ESP_K7x_S_a = 0.0E0;
    Double I_ESP_K6xy_S_a = 0.0E0;
    Double I_ESP_K6xz_S_a = 0.0E0;
    Double I_ESP_K5x2y_S_a = 0.0E0;
    Double I_ESP_K5xyz_S_a = 0.0E0;
    Double I_ESP_K5x2z_S_a = 0.0E0;
    Double I_ESP_K4x3y_S_a = 0.0E0;
    Double I_ESP_K4x2yz_S_a = 0.0E0;
    Double I_ESP_K4xy2z_S_a = 0.0E0;
    Double I_ESP_K4x3z_S_a = 0.0E0;
    Double I_ESP_K3x4y_S_a = 0.0E0;
    Double I_ESP_K3x3yz_S_a = 0.0E0;
    Double I_ESP_K3x2y2z_S_a = 0.0E0;
    Double I_ESP_K3xy3z_S_a = 0.0E0;
    Double I_ESP_K3x4z_S_a = 0.0E0;
    Double I_ESP_K2x5y_S_a = 0.0E0;
    Double I_ESP_K2x4yz_S_a = 0.0E0;
    Double I_ESP_K2x3y2z_S_a = 0.0E0;
    Double I_ESP_K2x2y3z_S_a = 0.0E0;
    Double I_ESP_K2xy4z_S_a = 0.0E0;
    Double I_ESP_K2x5z_S_a = 0.0E0;
    Double I_ESP_Kx6y_S_a = 0.0E0;
    Double I_ESP_Kx5yz_S_a = 0.0E0;
    Double I_ESP_Kx4y2z_S_a = 0.0E0;
    Double I_ESP_Kx3y3z_S_a = 0.0E0;
    Double I_ESP_Kx2y4z_S_a = 0.0E0;
    Double I_ESP_Kxy5z_S_a = 0.0E0;
    Double I_ESP_Kx6z_S_a = 0.0E0;
    Double I_ESP_K7y_S_a = 0.0E0;
    Double I_ESP_K6yz_S_a = 0.0E0;
    Double I_ESP_K5y2z_S_a = 0.0E0;
    Double I_ESP_K4y3z_S_a = 0.0E0;
    Double I_ESP_K3y4z_S_a = 0.0E0;
    Double I_ESP_K2y5z_S_a = 0.0E0;
    Double I_ESP_Ky6z_S_a = 0.0E0;
    Double I_ESP_K7z_S_a = 0.0E0;
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
    Double I_ESP_I6x_S = 0.0E0;
    Double I_ESP_I5xy_S = 0.0E0;
    Double I_ESP_I5xz_S = 0.0E0;
    Double I_ESP_I4x2y_S = 0.0E0;
    Double I_ESP_I4xyz_S = 0.0E0;
    Double I_ESP_I4x2z_S = 0.0E0;
    Double I_ESP_I3x3y_S = 0.0E0;
    Double I_ESP_I3x2yz_S = 0.0E0;
    Double I_ESP_I3xy2z_S = 0.0E0;
    Double I_ESP_I3x3z_S = 0.0E0;
    Double I_ESP_I2x4y_S = 0.0E0;
    Double I_ESP_I2x3yz_S = 0.0E0;
    Double I_ESP_I2x2y2z_S = 0.0E0;
    Double I_ESP_I2xy3z_S = 0.0E0;
    Double I_ESP_I2x4z_S = 0.0E0;
    Double I_ESP_Ix5y_S = 0.0E0;
    Double I_ESP_Ix4yz_S = 0.0E0;
    Double I_ESP_Ix3y2z_S = 0.0E0;
    Double I_ESP_Ix2y3z_S = 0.0E0;
    Double I_ESP_Ixy4z_S = 0.0E0;
    Double I_ESP_Ix5z_S = 0.0E0;
    Double I_ESP_I6y_S = 0.0E0;
    Double I_ESP_I5yz_S = 0.0E0;
    Double I_ESP_I4y2z_S = 0.0E0;
    Double I_ESP_I3y3z_S = 0.0E0;
    Double I_ESP_I2y4z_S = 0.0E0;
    Double I_ESP_Iy5z_S = 0.0E0;
    Double I_ESP_I6z_S = 0.0E0;
    Double I_ESP_H5x_S = 0.0E0;
    Double I_ESP_H4xy_S = 0.0E0;
    Double I_ESP_H4xz_S = 0.0E0;
    Double I_ESP_H3x2y_S = 0.0E0;
    Double I_ESP_H3xyz_S = 0.0E0;
    Double I_ESP_H3x2z_S = 0.0E0;
    Double I_ESP_H2x3y_S = 0.0E0;
    Double I_ESP_H2x2yz_S = 0.0E0;
    Double I_ESP_H2xy2z_S = 0.0E0;
    Double I_ESP_H2x3z_S = 0.0E0;
    Double I_ESP_Hx4y_S = 0.0E0;
    Double I_ESP_Hx3yz_S = 0.0E0;
    Double I_ESP_Hx2y2z_S = 0.0E0;
    Double I_ESP_Hxy3z_S = 0.0E0;
    Double I_ESP_Hx4z_S = 0.0E0;
    Double I_ESP_H5y_S = 0.0E0;
    Double I_ESP_H4yz_S = 0.0E0;
    Double I_ESP_H3y2z_S = 0.0E0;
    Double I_ESP_H2y3z_S = 0.0E0;
    Double I_ESP_Hy4z_S = 0.0E0;
    Double I_ESP_H5z_S = 0.0E0;
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
         * shell quartet name: SQ_ESP_L_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_L_S_a_coefs = alpha;
        I_ESP_L8x_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_K_S_a_coefs = alpha;
        I_ESP_K7x_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7z_S_vrr;

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
         * shell quartet name: SQ_ESP_I_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_I6x_S += I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S += I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S += I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S += I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S += I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S += I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S += I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S += I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S += I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S += I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S += I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S += I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S += I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S += I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S += I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S += I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S += I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S += I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S += I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S += I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S += I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S += I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S += I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S += I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S += I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S += I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S += I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S += I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_H5x_S += I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S += I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S += I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S += I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S += I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S += I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S += I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S += I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S += I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S += I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S += I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S += I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S += I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S += I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S += I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S += I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S += I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S += I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S += I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S += I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S += I_ESP_H5z_S_vrr;

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
     * shell quartet name: SQ_ESP_G_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S
     * RHS shell quartet name: SQ_ESP_G_S
     ************************************************************/
    Double I_ESP_G4x_Px = I_ESP_H5x_S+ABX*I_ESP_G4x_S;
    Double I_ESP_G3xy_Px = I_ESP_H4xy_S+ABX*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Px = I_ESP_H4xz_S+ABX*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Px = I_ESP_H3x2y_S+ABX*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Px = I_ESP_H3xyz_S+ABX*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Px = I_ESP_H3x2z_S+ABX*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Px = I_ESP_H2x3y_S+ABX*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Px = I_ESP_H2x2yz_S+ABX*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Px = I_ESP_H2xy2z_S+ABX*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Px = I_ESP_H2x3z_S+ABX*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Px = I_ESP_Hx4y_S+ABX*I_ESP_G4y_S;
    Double I_ESP_G3yz_Px = I_ESP_Hx3yz_S+ABX*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Px = I_ESP_Hx2y2z_S+ABX*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Px = I_ESP_Hxy3z_S+ABX*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Px = I_ESP_Hx4z_S+ABX*I_ESP_G4z_S;
    Double I_ESP_G4x_Py = I_ESP_H4xy_S+ABY*I_ESP_G4x_S;
    Double I_ESP_G3xy_Py = I_ESP_H3x2y_S+ABY*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Py = I_ESP_H3xyz_S+ABY*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Py = I_ESP_H2x3y_S+ABY*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Py = I_ESP_H2x2yz_S+ABY*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Py = I_ESP_H2xy2z_S+ABY*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Py = I_ESP_Hx4y_S+ABY*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Py = I_ESP_Hx3yz_S+ABY*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Py = I_ESP_Hx2y2z_S+ABY*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Py = I_ESP_Hxy3z_S+ABY*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Py = I_ESP_H5y_S+ABY*I_ESP_G4y_S;
    Double I_ESP_G3yz_Py = I_ESP_H4yz_S+ABY*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Py = I_ESP_H3y2z_S+ABY*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Py = I_ESP_H2y3z_S+ABY*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Py = I_ESP_Hy4z_S+ABY*I_ESP_G4z_S;
    Double I_ESP_G4x_Pz = I_ESP_H4xz_S+ABZ*I_ESP_G4x_S;
    Double I_ESP_G3xy_Pz = I_ESP_H3xyz_S+ABZ*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Pz = I_ESP_H3x2z_S+ABZ*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Pz = I_ESP_H2x2yz_S+ABZ*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Pz = I_ESP_H2xy2z_S+ABZ*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Pz = I_ESP_H2x3z_S+ABZ*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Pz = I_ESP_Hx3yz_S+ABZ*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Pz = I_ESP_Hx2y2z_S+ABZ*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Pz = I_ESP_Hxy3z_S+ABZ*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Pz = I_ESP_Hx4z_S+ABZ*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Pz = I_ESP_H4yz_S+ABZ*I_ESP_G4y_S;
    Double I_ESP_G3yz_Pz = I_ESP_H3y2z_S+ABZ*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Pz = I_ESP_H2y3z_S+ABZ*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Pz = I_ESP_Hy4z_S+ABZ*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Pz = I_ESP_H5z_S+ABZ*I_ESP_G4z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_F_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_P
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    Double I_ESP_F3x_D2x = I_ESP_G4x_Px+ABX*I_ESP_F3x_Px;
    Double I_ESP_F2xy_D2x = I_ESP_G3xy_Px+ABX*I_ESP_F2xy_Px;
    Double I_ESP_F2xz_D2x = I_ESP_G3xz_Px+ABX*I_ESP_F2xz_Px;
    Double I_ESP_Fx2y_D2x = I_ESP_G2x2y_Px+ABX*I_ESP_Fx2y_Px;
    Double I_ESP_Fxyz_D2x = I_ESP_G2xyz_Px+ABX*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_D2x = I_ESP_G2x2z_Px+ABX*I_ESP_Fx2z_Px;
    Double I_ESP_F3y_D2x = I_ESP_Gx3y_Px+ABX*I_ESP_F3y_Px;
    Double I_ESP_F2yz_D2x = I_ESP_Gx2yz_Px+ABX*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_D2x = I_ESP_Gxy2z_Px+ABX*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_D2x = I_ESP_Gx3z_Px+ABX*I_ESP_F3z_Px;
    Double I_ESP_F3x_Dxy = I_ESP_G3xy_Px+ABY*I_ESP_F3x_Px;
    Double I_ESP_F2xy_Dxy = I_ESP_G2x2y_Px+ABY*I_ESP_F2xy_Px;
    Double I_ESP_F2xz_Dxy = I_ESP_G2xyz_Px+ABY*I_ESP_F2xz_Px;
    Double I_ESP_Fx2y_Dxy = I_ESP_Gx3y_Px+ABY*I_ESP_Fx2y_Px;
    Double I_ESP_Fxyz_Dxy = I_ESP_Gx2yz_Px+ABY*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_Dxy = I_ESP_Gxy2z_Px+ABY*I_ESP_Fx2z_Px;
    Double I_ESP_F3y_Dxy = I_ESP_G4y_Px+ABY*I_ESP_F3y_Px;
    Double I_ESP_F2yz_Dxy = I_ESP_G3yz_Px+ABY*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_Dxy = I_ESP_G2y2z_Px+ABY*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_Dxy = I_ESP_Gy3z_Px+ABY*I_ESP_F3z_Px;
    Double I_ESP_F3x_D2y = I_ESP_G3xy_Py+ABY*I_ESP_F3x_Py;
    Double I_ESP_F2xy_D2y = I_ESP_G2x2y_Py+ABY*I_ESP_F2xy_Py;
    Double I_ESP_F2xz_D2y = I_ESP_G2xyz_Py+ABY*I_ESP_F2xz_Py;
    Double I_ESP_Fx2y_D2y = I_ESP_Gx3y_Py+ABY*I_ESP_Fx2y_Py;
    Double I_ESP_Fxyz_D2y = I_ESP_Gx2yz_Py+ABY*I_ESP_Fxyz_Py;
    Double I_ESP_Fx2z_D2y = I_ESP_Gxy2z_Py+ABY*I_ESP_Fx2z_Py;
    Double I_ESP_F3y_D2y = I_ESP_G4y_Py+ABY*I_ESP_F3y_Py;
    Double I_ESP_F2yz_D2y = I_ESP_G3yz_Py+ABY*I_ESP_F2yz_Py;
    Double I_ESP_Fy2z_D2y = I_ESP_G2y2z_Py+ABY*I_ESP_Fy2z_Py;
    Double I_ESP_F3z_D2y = I_ESP_Gy3z_Py+ABY*I_ESP_F3z_Py;
    Double I_ESP_F3x_D2z = I_ESP_G3xz_Pz+ABZ*I_ESP_F3x_Pz;
    Double I_ESP_F2xy_D2z = I_ESP_G2xyz_Pz+ABZ*I_ESP_F2xy_Pz;
    Double I_ESP_F2xz_D2z = I_ESP_G2x2z_Pz+ABZ*I_ESP_F2xz_Pz;
    Double I_ESP_Fx2y_D2z = I_ESP_Gx2yz_Pz+ABZ*I_ESP_Fx2y_Pz;
    Double I_ESP_Fxyz_D2z = I_ESP_Gxy2z_Pz+ABZ*I_ESP_Fxyz_Pz;
    Double I_ESP_Fx2z_D2z = I_ESP_Gx3z_Pz+ABZ*I_ESP_Fx2z_Pz;
    Double I_ESP_F3y_D2z = I_ESP_G3yz_Pz+ABZ*I_ESP_F3y_Pz;
    Double I_ESP_F2yz_D2z = I_ESP_G2y2z_Pz+ABZ*I_ESP_F2yz_Pz;
    Double I_ESP_Fy2z_D2z = I_ESP_Gy3z_Pz+ABZ*I_ESP_Fy2z_Pz;
    Double I_ESP_F3z_D2z = I_ESP_G4z_Pz+ABZ*I_ESP_F3z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 14 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S
     * RHS shell quartet name: SQ_ESP_H_S
     ************************************************************/
    Double I_ESP_H5x_Px = I_ESP_I6x_S+ABX*I_ESP_H5x_S;
    Double I_ESP_H4xy_Px = I_ESP_I5xy_S+ABX*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Px = I_ESP_I5xz_S+ABX*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Px = I_ESP_I4x2y_S+ABX*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Px = I_ESP_I4xyz_S+ABX*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Px = I_ESP_I4x2z_S+ABX*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Px = I_ESP_I3x3y_S+ABX*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Px = I_ESP_I3x2yz_S+ABX*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Px = I_ESP_I3xy2z_S+ABX*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Px = I_ESP_I3x3z_S+ABX*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Px = I_ESP_I2x4y_S+ABX*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Px = I_ESP_I2x3yz_S+ABX*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Px = I_ESP_I2x2y2z_S+ABX*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Px = I_ESP_I2xy3z_S+ABX*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Px = I_ESP_I2x4z_S+ABX*I_ESP_Hx4z_S;
    Double I_ESP_H4yz_Px = I_ESP_Ix4yz_S+ABX*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Px = I_ESP_Ix3y2z_S+ABX*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Px = I_ESP_Ix2y3z_S+ABX*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Px = I_ESP_Ixy4z_S+ABX*I_ESP_Hy4z_S;
    Double I_ESP_H4xy_Py = I_ESP_I4x2y_S+ABY*I_ESP_H4xy_S;
    Double I_ESP_H3x2y_Py = I_ESP_I3x3y_S+ABY*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Py = I_ESP_I3x2yz_S+ABY*I_ESP_H3xyz_S;
    Double I_ESP_H2x3y_Py = I_ESP_I2x4y_S+ABY*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Py = I_ESP_I2x3yz_S+ABY*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Py = I_ESP_I2x2y2z_S+ABY*I_ESP_H2xy2z_S;
    Double I_ESP_Hx4y_Py = I_ESP_Ix5y_S+ABY*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Py = I_ESP_Ix4yz_S+ABY*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Py = I_ESP_Ix3y2z_S+ABY*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Py = I_ESP_Ix2y3z_S+ABY*I_ESP_Hxy3z_S;
    Double I_ESP_H5y_Py = I_ESP_I6y_S+ABY*I_ESP_H5y_S;
    Double I_ESP_H4yz_Py = I_ESP_I5yz_S+ABY*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Py = I_ESP_I4y2z_S+ABY*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Py = I_ESP_I3y3z_S+ABY*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Py = I_ESP_I2y4z_S+ABY*I_ESP_Hy4z_S;
    Double I_ESP_H4xz_Pz = I_ESP_I4x2z_S+ABZ*I_ESP_H4xz_S;
    Double I_ESP_H3xyz_Pz = I_ESP_I3xy2z_S+ABZ*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Pz = I_ESP_I3x3z_S+ABZ*I_ESP_H3x2z_S;
    Double I_ESP_H2x2yz_Pz = I_ESP_I2x2y2z_S+ABZ*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Pz = I_ESP_I2xy3z_S+ABZ*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Pz = I_ESP_I2x4z_S+ABZ*I_ESP_H2x3z_S;
    Double I_ESP_Hx3yz_Pz = I_ESP_Ix3y2z_S+ABZ*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Pz = I_ESP_Ix2y3z_S+ABZ*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Pz = I_ESP_Ixy4z_S+ABZ*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Pz = I_ESP_Ix5z_S+ABZ*I_ESP_Hx4z_S;
    Double I_ESP_H4yz_Pz = I_ESP_I4y2z_S+ABZ*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Pz = I_ESP_I3y3z_S+ABZ*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Pz = I_ESP_I2y4z_S+ABZ*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Pz = I_ESP_Iy5z_S+ABZ*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Pz = I_ESP_I6z_S+ABZ*I_ESP_H5z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 35 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P
     * RHS shell quartet name: SQ_ESP_G_P
     ************************************************************/
    Double I_ESP_G4x_D2x = I_ESP_H5x_Px+ABX*I_ESP_G4x_Px;
    Double I_ESP_G3xy_D2x = I_ESP_H4xy_Px+ABX*I_ESP_G3xy_Px;
    Double I_ESP_G3xz_D2x = I_ESP_H4xz_Px+ABX*I_ESP_G3xz_Px;
    Double I_ESP_G2x2y_D2x = I_ESP_H3x2y_Px+ABX*I_ESP_G2x2y_Px;
    Double I_ESP_G2xyz_D2x = I_ESP_H3xyz_Px+ABX*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_D2x = I_ESP_H3x2z_Px+ABX*I_ESP_G2x2z_Px;
    Double I_ESP_Gx3y_D2x = I_ESP_H2x3y_Px+ABX*I_ESP_Gx3y_Px;
    Double I_ESP_Gx2yz_D2x = I_ESP_H2x2yz_Px+ABX*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_D2x = I_ESP_H2xy2z_Px+ABX*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_D2x = I_ESP_H2x3z_Px+ABX*I_ESP_Gx3z_Px;
    Double I_ESP_G4y_D2x = I_ESP_Hx4y_Px+ABX*I_ESP_G4y_Px;
    Double I_ESP_G3yz_D2x = I_ESP_Hx3yz_Px+ABX*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_D2x = I_ESP_Hx2y2z_Px+ABX*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_D2x = I_ESP_Hxy3z_Px+ABX*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_D2x = I_ESP_Hx4z_Px+ABX*I_ESP_G4z_Px;
    Double I_ESP_G3xz_Dxy = I_ESP_H3xyz_Px+ABY*I_ESP_G3xz_Px;
    Double I_ESP_G2xyz_Dxy = I_ESP_H2x2yz_Px+ABY*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_Dxy = I_ESP_H2xy2z_Px+ABY*I_ESP_G2x2z_Px;
    Double I_ESP_Gx2yz_Dxy = I_ESP_Hx3yz_Px+ABY*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_Dxy = I_ESP_Hx2y2z_Px+ABY*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_Dxy = I_ESP_Hxy3z_Px+ABY*I_ESP_Gx3z_Px;
    Double I_ESP_G3yz_Dxy = I_ESP_H4yz_Px+ABY*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_Dxy = I_ESP_H3y2z_Px+ABY*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_Dxy = I_ESP_H2y3z_Px+ABY*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_Dxy = I_ESP_Hy4z_Px+ABY*I_ESP_G4z_Px;
    Double I_ESP_G4x_D2y = I_ESP_H4xy_Py+ABY*I_ESP_G4x_Py;
    Double I_ESP_G3xy_D2y = I_ESP_H3x2y_Py+ABY*I_ESP_G3xy_Py;
    Double I_ESP_G3xz_D2y = I_ESP_H3xyz_Py+ABY*I_ESP_G3xz_Py;
    Double I_ESP_G2x2y_D2y = I_ESP_H2x3y_Py+ABY*I_ESP_G2x2y_Py;
    Double I_ESP_G2xyz_D2y = I_ESP_H2x2yz_Py+ABY*I_ESP_G2xyz_Py;
    Double I_ESP_G2x2z_D2y = I_ESP_H2xy2z_Py+ABY*I_ESP_G2x2z_Py;
    Double I_ESP_Gx3y_D2y = I_ESP_Hx4y_Py+ABY*I_ESP_Gx3y_Py;
    Double I_ESP_Gx2yz_D2y = I_ESP_Hx3yz_Py+ABY*I_ESP_Gx2yz_Py;
    Double I_ESP_Gxy2z_D2y = I_ESP_Hx2y2z_Py+ABY*I_ESP_Gxy2z_Py;
    Double I_ESP_Gx3z_D2y = I_ESP_Hxy3z_Py+ABY*I_ESP_Gx3z_Py;
    Double I_ESP_G4y_D2y = I_ESP_H5y_Py+ABY*I_ESP_G4y_Py;
    Double I_ESP_G3yz_D2y = I_ESP_H4yz_Py+ABY*I_ESP_G3yz_Py;
    Double I_ESP_G2y2z_D2y = I_ESP_H3y2z_Py+ABY*I_ESP_G2y2z_Py;
    Double I_ESP_Gy3z_D2y = I_ESP_H2y3z_Py+ABY*I_ESP_Gy3z_Py;
    Double I_ESP_G4z_D2y = I_ESP_Hy4z_Py+ABY*I_ESP_G4z_Py;
    Double I_ESP_G4x_D2z = I_ESP_H4xz_Pz+ABZ*I_ESP_G4x_Pz;
    Double I_ESP_G3xy_D2z = I_ESP_H3xyz_Pz+ABZ*I_ESP_G3xy_Pz;
    Double I_ESP_G3xz_D2z = I_ESP_H3x2z_Pz+ABZ*I_ESP_G3xz_Pz;
    Double I_ESP_G2x2y_D2z = I_ESP_H2x2yz_Pz+ABZ*I_ESP_G2x2y_Pz;
    Double I_ESP_G2xyz_D2z = I_ESP_H2xy2z_Pz+ABZ*I_ESP_G2xyz_Pz;
    Double I_ESP_G2x2z_D2z = I_ESP_H2x3z_Pz+ABZ*I_ESP_G2x2z_Pz;
    Double I_ESP_Gx3y_D2z = I_ESP_Hx3yz_Pz+ABZ*I_ESP_Gx3y_Pz;
    Double I_ESP_Gx2yz_D2z = I_ESP_Hx2y2z_Pz+ABZ*I_ESP_Gx2yz_Pz;
    Double I_ESP_Gxy2z_D2z = I_ESP_Hxy3z_Pz+ABZ*I_ESP_Gxy2z_Pz;
    Double I_ESP_Gx3z_D2z = I_ESP_Hx4z_Pz+ABZ*I_ESP_Gx3z_Pz;
    Double I_ESP_G4y_D2z = I_ESP_H4yz_Pz+ABZ*I_ESP_G4y_Pz;
    Double I_ESP_G3yz_D2z = I_ESP_H3y2z_Pz+ABZ*I_ESP_G3yz_Pz;
    Double I_ESP_G2y2z_D2z = I_ESP_H2y3z_Pz+ABZ*I_ESP_G2y2z_Pz;
    Double I_ESP_Gy3z_D2z = I_ESP_Hy4z_Pz+ABZ*I_ESP_Gy3z_Pz;
    Double I_ESP_G4z_D2z = I_ESP_H5z_Pz+ABZ*I_ESP_G4z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_D
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    Double I_ESP_F3x_F3x = I_ESP_G4x_D2x+ABX*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F3x = I_ESP_G3xy_D2x+ABX*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F3x = I_ESP_G3xz_D2x+ABX*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F3x = I_ESP_G2x2y_D2x+ABX*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F3x = I_ESP_G2xyz_D2x+ABX*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F3x = I_ESP_G2x2z_D2x+ABX*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F3x = I_ESP_Gx3y_D2x+ABX*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F3x = I_ESP_Gx2yz_D2x+ABX*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F3x = I_ESP_Gxy2z_D2x+ABX*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F3x = I_ESP_Gx3z_D2x+ABX*I_ESP_F3z_D2x;
    Double I_ESP_F3x_F2xy = I_ESP_G3xy_D2x+ABY*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F2xy = I_ESP_G2x2y_D2x+ABY*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F2xy = I_ESP_G2xyz_D2x+ABY*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F2xy = I_ESP_Gx3y_D2x+ABY*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F2xy = I_ESP_Gx2yz_D2x+ABY*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xy = I_ESP_Gxy2z_D2x+ABY*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F2xy = I_ESP_G4y_D2x+ABY*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F2xy = I_ESP_G3yz_D2x+ABY*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xy = I_ESP_G2y2z_D2x+ABY*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xy = I_ESP_Gy3z_D2x+ABY*I_ESP_F3z_D2x;
    Double I_ESP_F3x_F2xz = I_ESP_G3xz_D2x+ABZ*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F2xz = I_ESP_G2xyz_D2x+ABZ*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F2xz = I_ESP_G2x2z_D2x+ABZ*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F2xz = I_ESP_Gx2yz_D2x+ABZ*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F2xz = I_ESP_Gxy2z_D2x+ABZ*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xz = I_ESP_Gx3z_D2x+ABZ*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F2xz = I_ESP_G3yz_D2x+ABZ*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F2xz = I_ESP_G2y2z_D2x+ABZ*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xz = I_ESP_Gy3z_D2x+ABZ*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xz = I_ESP_G4z_D2x+ABZ*I_ESP_F3z_D2x;
    Double I_ESP_F3x_Fx2y = I_ESP_G4x_D2y+ABX*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_Fx2y = I_ESP_G3xy_D2y+ABX*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_Fx2y = I_ESP_G3xz_D2y+ABX*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_Fx2y = I_ESP_G2x2y_D2y+ABX*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_Fx2y = I_ESP_G2xyz_D2y+ABX*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_Fx2y = I_ESP_G2x2z_D2y+ABX*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_Fx2y = I_ESP_Gx3y_D2y+ABX*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_Fx2y = I_ESP_Gx2yz_D2y+ABX*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_Fx2y = I_ESP_Gxy2z_D2y+ABX*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_Fx2y = I_ESP_Gx3z_D2y+ABX*I_ESP_F3z_D2y;
    Double I_ESP_F3x_Fxyz = I_ESP_G3xz_Dxy+ABZ*I_ESP_F3x_Dxy;
    Double I_ESP_F2xy_Fxyz = I_ESP_G2xyz_Dxy+ABZ*I_ESP_F2xy_Dxy;
    Double I_ESP_F2xz_Fxyz = I_ESP_G2x2z_Dxy+ABZ*I_ESP_F2xz_Dxy;
    Double I_ESP_Fx2y_Fxyz = I_ESP_Gx2yz_Dxy+ABZ*I_ESP_Fx2y_Dxy;
    Double I_ESP_Fxyz_Fxyz = I_ESP_Gxy2z_Dxy+ABZ*I_ESP_Fxyz_Dxy;
    Double I_ESP_Fx2z_Fxyz = I_ESP_Gx3z_Dxy+ABZ*I_ESP_Fx2z_Dxy;
    Double I_ESP_F3y_Fxyz = I_ESP_G3yz_Dxy+ABZ*I_ESP_F3y_Dxy;
    Double I_ESP_F2yz_Fxyz = I_ESP_G2y2z_Dxy+ABZ*I_ESP_F2yz_Dxy;
    Double I_ESP_Fy2z_Fxyz = I_ESP_Gy3z_Dxy+ABZ*I_ESP_Fy2z_Dxy;
    Double I_ESP_F3z_Fxyz = I_ESP_G4z_Dxy+ABZ*I_ESP_F3z_Dxy;
    Double I_ESP_F3x_Fx2z = I_ESP_G4x_D2z+ABX*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_Fx2z = I_ESP_G3xy_D2z+ABX*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_Fx2z = I_ESP_G3xz_D2z+ABX*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_Fx2z = I_ESP_G2x2y_D2z+ABX*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_Fx2z = I_ESP_G2xyz_D2z+ABX*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_Fx2z = I_ESP_G2x2z_D2z+ABX*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_Fx2z = I_ESP_Gx3y_D2z+ABX*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_Fx2z = I_ESP_Gx2yz_D2z+ABX*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_Fx2z = I_ESP_Gxy2z_D2z+ABX*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_Fx2z = I_ESP_Gx3z_D2z+ABX*I_ESP_F3z_D2z;
    Double I_ESP_F3x_F3y = I_ESP_G3xy_D2y+ABY*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_F3y = I_ESP_G2x2y_D2y+ABY*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_F3y = I_ESP_G2xyz_D2y+ABY*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_F3y = I_ESP_Gx3y_D2y+ABY*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_F3y = I_ESP_Gx2yz_D2y+ABY*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_F3y = I_ESP_Gxy2z_D2y+ABY*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_F3y = I_ESP_G4y_D2y+ABY*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_F3y = I_ESP_G3yz_D2y+ABY*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_F3y = I_ESP_G2y2z_D2y+ABY*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_F3y = I_ESP_Gy3z_D2y+ABY*I_ESP_F3z_D2y;
    Double I_ESP_F3x_F2yz = I_ESP_G3xz_D2y+ABZ*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_F2yz = I_ESP_G2xyz_D2y+ABZ*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_F2yz = I_ESP_G2x2z_D2y+ABZ*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_F2yz = I_ESP_Gx2yz_D2y+ABZ*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_F2yz = I_ESP_Gxy2z_D2y+ABZ*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_F2yz = I_ESP_Gx3z_D2y+ABZ*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_F2yz = I_ESP_G3yz_D2y+ABZ*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_F2yz = I_ESP_G2y2z_D2y+ABZ*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_F2yz = I_ESP_Gy3z_D2y+ABZ*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_F2yz = I_ESP_G4z_D2y+ABZ*I_ESP_F3z_D2y;
    Double I_ESP_F3x_Fy2z = I_ESP_G3xy_D2z+ABY*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_Fy2z = I_ESP_G2x2y_D2z+ABY*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_Fy2z = I_ESP_G2xyz_D2z+ABY*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_Fy2z = I_ESP_Gx3y_D2z+ABY*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_Fy2z = I_ESP_Gx2yz_D2z+ABY*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_Fy2z = I_ESP_Gxy2z_D2z+ABY*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_Fy2z = I_ESP_G4y_D2z+ABY*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_Fy2z = I_ESP_G3yz_D2z+ABY*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_Fy2z = I_ESP_G2y2z_D2z+ABY*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_Fy2z = I_ESP_Gy3z_D2z+ABY*I_ESP_F3z_D2z;
    Double I_ESP_F3x_F3z = I_ESP_G3xz_D2z+ABZ*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_F3z = I_ESP_G2xyz_D2z+ABZ*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_F3z = I_ESP_G2x2z_D2z+ABZ*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_F3z = I_ESP_Gx2yz_D2z+ABZ*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_F3z = I_ESP_Gxy2z_D2z+ABZ*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_F3z = I_ESP_Gx3z_D2z+ABZ*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_F3z = I_ESP_G3yz_D2z+ABZ*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_F3z = I_ESP_G2y2z_D2z+ABZ*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_F3z = I_ESP_Gy3z_D2z+ABZ*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_F3z = I_ESP_G4z_D2z+ABZ*I_ESP_F3z_D2z;

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
     * shell quartet name: SQ_ESP_I_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S_a
     * RHS shell quartet name: SQ_ESP_I_S_a
     ************************************************************/
    Double I_ESP_I6x_Px_a = I_ESP_K7x_S_a+ABX*I_ESP_I6x_S_a;
    Double I_ESP_I5xy_Px_a = I_ESP_K6xy_S_a+ABX*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Px_a = I_ESP_K6xz_S_a+ABX*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Px_a = I_ESP_K5x2y_S_a+ABX*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Px_a = I_ESP_K5xyz_S_a+ABX*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Px_a = I_ESP_K5x2z_S_a+ABX*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Px_a = I_ESP_K4x3y_S_a+ABX*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Px_a = I_ESP_K4x2yz_S_a+ABX*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Px_a = I_ESP_K4xy2z_S_a+ABX*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Px_a = I_ESP_K4x3z_S_a+ABX*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Px_a = I_ESP_K3x4y_S_a+ABX*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Px_a = I_ESP_K3x3yz_S_a+ABX*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Px_a = I_ESP_K3x2y2z_S_a+ABX*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Px_a = I_ESP_K3xy3z_S_a+ABX*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Px_a = I_ESP_K3x4z_S_a+ABX*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Px_a = I_ESP_K2x5y_S_a+ABX*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Px_a = I_ESP_K2x4yz_S_a+ABX*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Px_a = I_ESP_K2x3y2z_S_a+ABX*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Px_a = I_ESP_K2x2y3z_S_a+ABX*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Px_a = I_ESP_K2xy4z_S_a+ABX*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Px_a = I_ESP_K2x5z_S_a+ABX*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Px_a = I_ESP_Kx6y_S_a+ABX*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Px_a = I_ESP_Kx5yz_S_a+ABX*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Px_a = I_ESP_Kx4y2z_S_a+ABX*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Px_a = I_ESP_Kx3y3z_S_a+ABX*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Px_a = I_ESP_Kx2y4z_S_a+ABX*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Px_a = I_ESP_Kxy5z_S_a+ABX*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Px_a = I_ESP_Kx6z_S_a+ABX*I_ESP_I6z_S_a;
    Double I_ESP_I6x_Py_a = I_ESP_K6xy_S_a+ABY*I_ESP_I6x_S_a;
    Double I_ESP_I5xy_Py_a = I_ESP_K5x2y_S_a+ABY*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Py_a = I_ESP_K5xyz_S_a+ABY*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Py_a = I_ESP_K4x3y_S_a+ABY*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Py_a = I_ESP_K4x2yz_S_a+ABY*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Py_a = I_ESP_K4xy2z_S_a+ABY*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Py_a = I_ESP_K3x4y_S_a+ABY*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Py_a = I_ESP_K3x3yz_S_a+ABY*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Py_a = I_ESP_K3x2y2z_S_a+ABY*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Py_a = I_ESP_K3xy3z_S_a+ABY*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Py_a = I_ESP_K2x5y_S_a+ABY*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Py_a = I_ESP_K2x4yz_S_a+ABY*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Py_a = I_ESP_K2x3y2z_S_a+ABY*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Py_a = I_ESP_K2x2y3z_S_a+ABY*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Py_a = I_ESP_K2xy4z_S_a+ABY*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Py_a = I_ESP_Kx6y_S_a+ABY*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Py_a = I_ESP_Kx5yz_S_a+ABY*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Py_a = I_ESP_Kx4y2z_S_a+ABY*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Py_a = I_ESP_Kx3y3z_S_a+ABY*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Py_a = I_ESP_Kx2y4z_S_a+ABY*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Py_a = I_ESP_Kxy5z_S_a+ABY*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Py_a = I_ESP_K7y_S_a+ABY*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Py_a = I_ESP_K6yz_S_a+ABY*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Py_a = I_ESP_K5y2z_S_a+ABY*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Py_a = I_ESP_K4y3z_S_a+ABY*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Py_a = I_ESP_K3y4z_S_a+ABY*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Py_a = I_ESP_K2y5z_S_a+ABY*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Py_a = I_ESP_Ky6z_S_a+ABY*I_ESP_I6z_S_a;
    Double I_ESP_I6x_Pz_a = I_ESP_K6xz_S_a+ABZ*I_ESP_I6x_S_a;
    Double I_ESP_I5xy_Pz_a = I_ESP_K5xyz_S_a+ABZ*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Pz_a = I_ESP_K5x2z_S_a+ABZ*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Pz_a = I_ESP_K4x2yz_S_a+ABZ*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Pz_a = I_ESP_K4xy2z_S_a+ABZ*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Pz_a = I_ESP_K4x3z_S_a+ABZ*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Pz_a = I_ESP_K3x3yz_S_a+ABZ*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Pz_a = I_ESP_K3x2y2z_S_a+ABZ*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Pz_a = I_ESP_K3xy3z_S_a+ABZ*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Pz_a = I_ESP_K3x4z_S_a+ABZ*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Pz_a = I_ESP_K2x4yz_S_a+ABZ*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Pz_a = I_ESP_K2x3y2z_S_a+ABZ*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Pz_a = I_ESP_K2x2y3z_S_a+ABZ*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Pz_a = I_ESP_K2xy4z_S_a+ABZ*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Pz_a = I_ESP_K2x5z_S_a+ABZ*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Pz_a = I_ESP_Kx5yz_S_a+ABZ*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Pz_a = I_ESP_Kx4y2z_S_a+ABZ*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Pz_a = I_ESP_Kx3y3z_S_a+ABZ*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Pz_a = I_ESP_Kx2y4z_S_a+ABZ*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Pz_a = I_ESP_Kxy5z_S_a+ABZ*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Pz_a = I_ESP_Kx6z_S_a+ABZ*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Pz_a = I_ESP_K6yz_S_a+ABZ*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Pz_a = I_ESP_K5y2z_S_a+ABZ*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Pz_a = I_ESP_K4y3z_S_a+ABZ*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Pz_a = I_ESP_K3y4z_S_a+ABZ*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Pz_a = I_ESP_K2y5z_S_a+ABZ*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Pz_a = I_ESP_Ky6z_S_a+ABZ*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Pz_a = I_ESP_K7z_S_a+ABZ*I_ESP_I6z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     ************************************************************/
    Double I_ESP_H5x_D2x_a = I_ESP_I6x_Px_a+ABX*I_ESP_H5x_Px_a;
    Double I_ESP_H4xy_D2x_a = I_ESP_I5xy_Px_a+ABX*I_ESP_H4xy_Px_a;
    Double I_ESP_H4xz_D2x_a = I_ESP_I5xz_Px_a+ABX*I_ESP_H4xz_Px_a;
    Double I_ESP_H3x2y_D2x_a = I_ESP_I4x2y_Px_a+ABX*I_ESP_H3x2y_Px_a;
    Double I_ESP_H3xyz_D2x_a = I_ESP_I4xyz_Px_a+ABX*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_D2x_a = I_ESP_I4x2z_Px_a+ABX*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x3y_D2x_a = I_ESP_I3x3y_Px_a+ABX*I_ESP_H2x3y_Px_a;
    Double I_ESP_H2x2yz_D2x_a = I_ESP_I3x2yz_Px_a+ABX*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_D2x_a = I_ESP_I3xy2z_Px_a+ABX*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_D2x_a = I_ESP_I3x3z_Px_a+ABX*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx4y_D2x_a = I_ESP_I2x4y_Px_a+ABX*I_ESP_Hx4y_Px_a;
    Double I_ESP_Hx3yz_D2x_a = I_ESP_I2x3yz_Px_a+ABX*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_D2x_a = I_ESP_I2x2y2z_Px_a+ABX*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_D2x_a = I_ESP_I2xy3z_Px_a+ABX*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_D2x_a = I_ESP_I2x4z_Px_a+ABX*I_ESP_Hx4z_Px_a;
    Double I_ESP_H5y_D2x_a = I_ESP_Ix5y_Px_a+ABX*I_ESP_H5y_Px_a;
    Double I_ESP_H4yz_D2x_a = I_ESP_Ix4yz_Px_a+ABX*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_D2x_a = I_ESP_Ix3y2z_Px_a+ABX*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_D2x_a = I_ESP_Ix2y3z_Px_a+ABX*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_D2x_a = I_ESP_Ixy4z_Px_a+ABX*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_D2x_a = I_ESP_Ix5z_Px_a+ABX*I_ESP_H5z_Px_a;
    Double I_ESP_H5x_Dxy_a = I_ESP_I5xy_Px_a+ABY*I_ESP_H5x_Px_a;
    Double I_ESP_H4xy_Dxy_a = I_ESP_I4x2y_Px_a+ABY*I_ESP_H4xy_Px_a;
    Double I_ESP_H4xz_Dxy_a = I_ESP_I4xyz_Px_a+ABY*I_ESP_H4xz_Px_a;
    Double I_ESP_H3x2y_Dxy_a = I_ESP_I3x3y_Px_a+ABY*I_ESP_H3x2y_Px_a;
    Double I_ESP_H3xyz_Dxy_a = I_ESP_I3x2yz_Px_a+ABY*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_Dxy_a = I_ESP_I3xy2z_Px_a+ABY*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x3y_Dxy_a = I_ESP_I2x4y_Px_a+ABY*I_ESP_H2x3y_Px_a;
    Double I_ESP_H2x2yz_Dxy_a = I_ESP_I2x3yz_Px_a+ABY*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_Dxy_a = I_ESP_I2x2y2z_Px_a+ABY*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_Dxy_a = I_ESP_I2xy3z_Px_a+ABY*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx4y_Dxy_a = I_ESP_Ix5y_Px_a+ABY*I_ESP_Hx4y_Px_a;
    Double I_ESP_Hx3yz_Dxy_a = I_ESP_Ix4yz_Px_a+ABY*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_Dxy_a = I_ESP_Ix3y2z_Px_a+ABY*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_Dxy_a = I_ESP_Ix2y3z_Px_a+ABY*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_Dxy_a = I_ESP_Ixy4z_Px_a+ABY*I_ESP_Hx4z_Px_a;
    Double I_ESP_H5y_Dxy_a = I_ESP_I6y_Px_a+ABY*I_ESP_H5y_Px_a;
    Double I_ESP_H4yz_Dxy_a = I_ESP_I5yz_Px_a+ABY*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_Dxy_a = I_ESP_I4y2z_Px_a+ABY*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_Dxy_a = I_ESP_I3y3z_Px_a+ABY*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_Dxy_a = I_ESP_I2y4z_Px_a+ABY*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_Dxy_a = I_ESP_Iy5z_Px_a+ABY*I_ESP_H5z_Px_a;
    Double I_ESP_H5x_D2y_a = I_ESP_I5xy_Py_a+ABY*I_ESP_H5x_Py_a;
    Double I_ESP_H4xy_D2y_a = I_ESP_I4x2y_Py_a+ABY*I_ESP_H4xy_Py_a;
    Double I_ESP_H4xz_D2y_a = I_ESP_I4xyz_Py_a+ABY*I_ESP_H4xz_Py_a;
    Double I_ESP_H3x2y_D2y_a = I_ESP_I3x3y_Py_a+ABY*I_ESP_H3x2y_Py_a;
    Double I_ESP_H3xyz_D2y_a = I_ESP_I3x2yz_Py_a+ABY*I_ESP_H3xyz_Py_a;
    Double I_ESP_H3x2z_D2y_a = I_ESP_I3xy2z_Py_a+ABY*I_ESP_H3x2z_Py_a;
    Double I_ESP_H2x3y_D2y_a = I_ESP_I2x4y_Py_a+ABY*I_ESP_H2x3y_Py_a;
    Double I_ESP_H2x2yz_D2y_a = I_ESP_I2x3yz_Py_a+ABY*I_ESP_H2x2yz_Py_a;
    Double I_ESP_H2xy2z_D2y_a = I_ESP_I2x2y2z_Py_a+ABY*I_ESP_H2xy2z_Py_a;
    Double I_ESP_H2x3z_D2y_a = I_ESP_I2xy3z_Py_a+ABY*I_ESP_H2x3z_Py_a;
    Double I_ESP_Hx4y_D2y_a = I_ESP_Ix5y_Py_a+ABY*I_ESP_Hx4y_Py_a;
    Double I_ESP_Hx3yz_D2y_a = I_ESP_Ix4yz_Py_a+ABY*I_ESP_Hx3yz_Py_a;
    Double I_ESP_Hx2y2z_D2y_a = I_ESP_Ix3y2z_Py_a+ABY*I_ESP_Hx2y2z_Py_a;
    Double I_ESP_Hxy3z_D2y_a = I_ESP_Ix2y3z_Py_a+ABY*I_ESP_Hxy3z_Py_a;
    Double I_ESP_Hx4z_D2y_a = I_ESP_Ixy4z_Py_a+ABY*I_ESP_Hx4z_Py_a;
    Double I_ESP_H5y_D2y_a = I_ESP_I6y_Py_a+ABY*I_ESP_H5y_Py_a;
    Double I_ESP_H4yz_D2y_a = I_ESP_I5yz_Py_a+ABY*I_ESP_H4yz_Py_a;
    Double I_ESP_H3y2z_D2y_a = I_ESP_I4y2z_Py_a+ABY*I_ESP_H3y2z_Py_a;
    Double I_ESP_H2y3z_D2y_a = I_ESP_I3y3z_Py_a+ABY*I_ESP_H2y3z_Py_a;
    Double I_ESP_Hy4z_D2y_a = I_ESP_I2y4z_Py_a+ABY*I_ESP_Hy4z_Py_a;
    Double I_ESP_H5z_D2y_a = I_ESP_Iy5z_Py_a+ABY*I_ESP_H5z_Py_a;
    Double I_ESP_H5x_D2z_a = I_ESP_I5xz_Pz_a+ABZ*I_ESP_H5x_Pz_a;
    Double I_ESP_H4xy_D2z_a = I_ESP_I4xyz_Pz_a+ABZ*I_ESP_H4xy_Pz_a;
    Double I_ESP_H4xz_D2z_a = I_ESP_I4x2z_Pz_a+ABZ*I_ESP_H4xz_Pz_a;
    Double I_ESP_H3x2y_D2z_a = I_ESP_I3x2yz_Pz_a+ABZ*I_ESP_H3x2y_Pz_a;
    Double I_ESP_H3xyz_D2z_a = I_ESP_I3xy2z_Pz_a+ABZ*I_ESP_H3xyz_Pz_a;
    Double I_ESP_H3x2z_D2z_a = I_ESP_I3x3z_Pz_a+ABZ*I_ESP_H3x2z_Pz_a;
    Double I_ESP_H2x3y_D2z_a = I_ESP_I2x3yz_Pz_a+ABZ*I_ESP_H2x3y_Pz_a;
    Double I_ESP_H2x2yz_D2z_a = I_ESP_I2x2y2z_Pz_a+ABZ*I_ESP_H2x2yz_Pz_a;
    Double I_ESP_H2xy2z_D2z_a = I_ESP_I2xy3z_Pz_a+ABZ*I_ESP_H2xy2z_Pz_a;
    Double I_ESP_H2x3z_D2z_a = I_ESP_I2x4z_Pz_a+ABZ*I_ESP_H2x3z_Pz_a;
    Double I_ESP_Hx4y_D2z_a = I_ESP_Ix4yz_Pz_a+ABZ*I_ESP_Hx4y_Pz_a;
    Double I_ESP_Hx3yz_D2z_a = I_ESP_Ix3y2z_Pz_a+ABZ*I_ESP_Hx3yz_Pz_a;
    Double I_ESP_Hx2y2z_D2z_a = I_ESP_Ix2y3z_Pz_a+ABZ*I_ESP_Hx2y2z_Pz_a;
    Double I_ESP_Hxy3z_D2z_a = I_ESP_Ixy4z_Pz_a+ABZ*I_ESP_Hxy3z_Pz_a;
    Double I_ESP_Hx4z_D2z_a = I_ESP_Ix5z_Pz_a+ABZ*I_ESP_Hx4z_Pz_a;
    Double I_ESP_H5y_D2z_a = I_ESP_I5yz_Pz_a+ABZ*I_ESP_H5y_Pz_a;
    Double I_ESP_H4yz_D2z_a = I_ESP_I4y2z_Pz_a+ABZ*I_ESP_H4yz_Pz_a;
    Double I_ESP_H3y2z_D2z_a = I_ESP_I3y3z_Pz_a+ABZ*I_ESP_H3y2z_Pz_a;
    Double I_ESP_H2y3z_D2z_a = I_ESP_I2y4z_Pz_a+ABZ*I_ESP_H2y3z_Pz_a;
    Double I_ESP_Hy4z_D2z_a = I_ESP_Iy5z_Pz_a+ABZ*I_ESP_Hy4z_Pz_a;
    Double I_ESP_H5z_D2z_a = I_ESP_I6z_Pz_a+ABZ*I_ESP_H5z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S_a
     * RHS shell quartet name: SQ_ESP_K_S_a
     ************************************************************/
    Double I_ESP_K7x_Px_a = I_ESP_L8x_S_a+ABX*I_ESP_K7x_S_a;
    Double I_ESP_K6xy_Px_a = I_ESP_L7xy_S_a+ABX*I_ESP_K6xy_S_a;
    Double I_ESP_K6xz_Px_a = I_ESP_L7xz_S_a+ABX*I_ESP_K6xz_S_a;
    Double I_ESP_K5x2y_Px_a = I_ESP_L6x2y_S_a+ABX*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Px_a = I_ESP_L6xyz_S_a+ABX*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Px_a = I_ESP_L6x2z_S_a+ABX*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x3y_Px_a = I_ESP_L5x3y_S_a+ABX*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Px_a = I_ESP_L5x2yz_S_a+ABX*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Px_a = I_ESP_L5xy2z_S_a+ABX*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Px_a = I_ESP_L5x3z_S_a+ABX*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x4y_Px_a = I_ESP_L4x4y_S_a+ABX*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Px_a = I_ESP_L4x3yz_S_a+ABX*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Px_a = I_ESP_L4x2y2z_S_a+ABX*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Px_a = I_ESP_L4xy3z_S_a+ABX*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Px_a = I_ESP_L4x4z_S_a+ABX*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x5y_Px_a = I_ESP_L3x5y_S_a+ABX*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Px_a = I_ESP_L3x4yz_S_a+ABX*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Px_a = I_ESP_L3x3y2z_S_a+ABX*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Px_a = I_ESP_L3x2y3z_S_a+ABX*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Px_a = I_ESP_L3xy4z_S_a+ABX*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Px_a = I_ESP_L3x5z_S_a+ABX*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx6y_Px_a = I_ESP_L2x6y_S_a+ABX*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Px_a = I_ESP_L2x5yz_S_a+ABX*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Px_a = I_ESP_L2x4y2z_S_a+ABX*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Px_a = I_ESP_L2x3y3z_S_a+ABX*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Px_a = I_ESP_L2x2y4z_S_a+ABX*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Px_a = I_ESP_L2xy5z_S_a+ABX*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Px_a = I_ESP_L2x6z_S_a+ABX*I_ESP_Kx6z_S_a;
    Double I_ESP_K6yz_Px_a = I_ESP_Lx6yz_S_a+ABX*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Px_a = I_ESP_Lx5y2z_S_a+ABX*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Px_a = I_ESP_Lx4y3z_S_a+ABX*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Px_a = I_ESP_Lx3y4z_S_a+ABX*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Px_a = I_ESP_Lx2y5z_S_a+ABX*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Px_a = I_ESP_Lxy6z_S_a+ABX*I_ESP_Ky6z_S_a;
    Double I_ESP_K6xy_Py_a = I_ESP_L6x2y_S_a+ABY*I_ESP_K6xy_S_a;
    Double I_ESP_K5x2y_Py_a = I_ESP_L5x3y_S_a+ABY*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Py_a = I_ESP_L5x2yz_S_a+ABY*I_ESP_K5xyz_S_a;
    Double I_ESP_K4x3y_Py_a = I_ESP_L4x4y_S_a+ABY*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Py_a = I_ESP_L4x3yz_S_a+ABY*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Py_a = I_ESP_L4x2y2z_S_a+ABY*I_ESP_K4xy2z_S_a;
    Double I_ESP_K3x4y_Py_a = I_ESP_L3x5y_S_a+ABY*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Py_a = I_ESP_L3x4yz_S_a+ABY*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Py_a = I_ESP_L3x3y2z_S_a+ABY*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Py_a = I_ESP_L3x2y3z_S_a+ABY*I_ESP_K3xy3z_S_a;
    Double I_ESP_K2x5y_Py_a = I_ESP_L2x6y_S_a+ABY*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Py_a = I_ESP_L2x5yz_S_a+ABY*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Py_a = I_ESP_L2x4y2z_S_a+ABY*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Py_a = I_ESP_L2x3y3z_S_a+ABY*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Py_a = I_ESP_L2x2y4z_S_a+ABY*I_ESP_K2xy4z_S_a;
    Double I_ESP_Kx6y_Py_a = I_ESP_Lx7y_S_a+ABY*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Py_a = I_ESP_Lx6yz_S_a+ABY*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Py_a = I_ESP_Lx5y2z_S_a+ABY*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Py_a = I_ESP_Lx4y3z_S_a+ABY*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Py_a = I_ESP_Lx3y4z_S_a+ABY*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Py_a = I_ESP_Lx2y5z_S_a+ABY*I_ESP_Kxy5z_S_a;
    Double I_ESP_K7y_Py_a = I_ESP_L8y_S_a+ABY*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Py_a = I_ESP_L7yz_S_a+ABY*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Py_a = I_ESP_L6y2z_S_a+ABY*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Py_a = I_ESP_L5y3z_S_a+ABY*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Py_a = I_ESP_L4y4z_S_a+ABY*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Py_a = I_ESP_L3y5z_S_a+ABY*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Py_a = I_ESP_L2y6z_S_a+ABY*I_ESP_Ky6z_S_a;
    Double I_ESP_K6xz_Pz_a = I_ESP_L6x2z_S_a+ABZ*I_ESP_K6xz_S_a;
    Double I_ESP_K5xyz_Pz_a = I_ESP_L5xy2z_S_a+ABZ*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Pz_a = I_ESP_L5x3z_S_a+ABZ*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x2yz_Pz_a = I_ESP_L4x2y2z_S_a+ABZ*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Pz_a = I_ESP_L4xy3z_S_a+ABZ*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Pz_a = I_ESP_L4x4z_S_a+ABZ*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x3yz_Pz_a = I_ESP_L3x3y2z_S_a+ABZ*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Pz_a = I_ESP_L3x2y3z_S_a+ABZ*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Pz_a = I_ESP_L3xy4z_S_a+ABZ*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Pz_a = I_ESP_L3x5z_S_a+ABZ*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x4yz_Pz_a = I_ESP_L2x4y2z_S_a+ABZ*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Pz_a = I_ESP_L2x3y3z_S_a+ABZ*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Pz_a = I_ESP_L2x2y4z_S_a+ABZ*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Pz_a = I_ESP_L2xy5z_S_a+ABZ*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Pz_a = I_ESP_L2x6z_S_a+ABZ*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx5yz_Pz_a = I_ESP_Lx5y2z_S_a+ABZ*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Pz_a = I_ESP_Lx4y3z_S_a+ABZ*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Pz_a = I_ESP_Lx3y4z_S_a+ABZ*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Pz_a = I_ESP_Lx2y5z_S_a+ABZ*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Pz_a = I_ESP_Lxy6z_S_a+ABZ*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Pz_a = I_ESP_Lx7z_S_a+ABZ*I_ESP_Kx6z_S_a;
    Double I_ESP_K6yz_Pz_a = I_ESP_L6y2z_S_a+ABZ*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Pz_a = I_ESP_L5y3z_S_a+ABZ*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Pz_a = I_ESP_L4y4z_S_a+ABZ*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Pz_a = I_ESP_L3y5z_S_a+ABZ*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Pz_a = I_ESP_L2y6z_S_a+ABZ*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Pz_a = I_ESP_Ly7z_S_a+ABZ*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Pz_a = I_ESP_L8z_S_a+ABZ*I_ESP_K7z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 63 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_a
     * RHS shell quartet name: SQ_ESP_I_P_a
     ************************************************************/
    Double I_ESP_I6x_D2x_a = I_ESP_K7x_Px_a+ABX*I_ESP_I6x_Px_a;
    Double I_ESP_I5xy_D2x_a = I_ESP_K6xy_Px_a+ABX*I_ESP_I5xy_Px_a;
    Double I_ESP_I5xz_D2x_a = I_ESP_K6xz_Px_a+ABX*I_ESP_I5xz_Px_a;
    Double I_ESP_I4x2y_D2x_a = I_ESP_K5x2y_Px_a+ABX*I_ESP_I4x2y_Px_a;
    Double I_ESP_I4xyz_D2x_a = I_ESP_K5xyz_Px_a+ABX*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_D2x_a = I_ESP_K5x2z_Px_a+ABX*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x3y_D2x_a = I_ESP_K4x3y_Px_a+ABX*I_ESP_I3x3y_Px_a;
    Double I_ESP_I3x2yz_D2x_a = I_ESP_K4x2yz_Px_a+ABX*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_D2x_a = I_ESP_K4xy2z_Px_a+ABX*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_D2x_a = I_ESP_K4x3z_Px_a+ABX*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x4y_D2x_a = I_ESP_K3x4y_Px_a+ABX*I_ESP_I2x4y_Px_a;
    Double I_ESP_I2x3yz_D2x_a = I_ESP_K3x3yz_Px_a+ABX*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_D2x_a = I_ESP_K3x2y2z_Px_a+ABX*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_D2x_a = I_ESP_K3xy3z_Px_a+ABX*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_D2x_a = I_ESP_K3x4z_Px_a+ABX*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix5y_D2x_a = I_ESP_K2x5y_Px_a+ABX*I_ESP_Ix5y_Px_a;
    Double I_ESP_Ix4yz_D2x_a = I_ESP_K2x4yz_Px_a+ABX*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_D2x_a = I_ESP_K2x3y2z_Px_a+ABX*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_D2x_a = I_ESP_K2x2y3z_Px_a+ABX*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_D2x_a = I_ESP_K2xy4z_Px_a+ABX*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_D2x_a = I_ESP_K2x5z_Px_a+ABX*I_ESP_Ix5z_Px_a;
    Double I_ESP_I6y_D2x_a = I_ESP_Kx6y_Px_a+ABX*I_ESP_I6y_Px_a;
    Double I_ESP_I5yz_D2x_a = I_ESP_Kx5yz_Px_a+ABX*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_D2x_a = I_ESP_Kx4y2z_Px_a+ABX*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_D2x_a = I_ESP_Kx3y3z_Px_a+ABX*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_D2x_a = I_ESP_Kx2y4z_Px_a+ABX*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_D2x_a = I_ESP_Kxy5z_Px_a+ABX*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_D2x_a = I_ESP_Kx6z_Px_a+ABX*I_ESP_I6z_Px_a;
    Double I_ESP_I5xz_Dxy_a = I_ESP_K5xyz_Px_a+ABY*I_ESP_I5xz_Px_a;
    Double I_ESP_I4xyz_Dxy_a = I_ESP_K4x2yz_Px_a+ABY*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_Dxy_a = I_ESP_K4xy2z_Px_a+ABY*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x2yz_Dxy_a = I_ESP_K3x3yz_Px_a+ABY*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_Dxy_a = I_ESP_K3x2y2z_Px_a+ABY*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_Dxy_a = I_ESP_K3xy3z_Px_a+ABY*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x3yz_Dxy_a = I_ESP_K2x4yz_Px_a+ABY*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_Dxy_a = I_ESP_K2x3y2z_Px_a+ABY*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_Dxy_a = I_ESP_K2x2y3z_Px_a+ABY*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_Dxy_a = I_ESP_K2xy4z_Px_a+ABY*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix4yz_Dxy_a = I_ESP_Kx5yz_Px_a+ABY*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_Dxy_a = I_ESP_Kx4y2z_Px_a+ABY*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_Dxy_a = I_ESP_Kx3y3z_Px_a+ABY*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_Dxy_a = I_ESP_Kx2y4z_Px_a+ABY*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_Dxy_a = I_ESP_Kxy5z_Px_a+ABY*I_ESP_Ix5z_Px_a;
    Double I_ESP_I5yz_Dxy_a = I_ESP_K6yz_Px_a+ABY*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_Dxy_a = I_ESP_K5y2z_Px_a+ABY*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_Dxy_a = I_ESP_K4y3z_Px_a+ABY*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_Dxy_a = I_ESP_K3y4z_Px_a+ABY*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_Dxy_a = I_ESP_K2y5z_Px_a+ABY*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_Dxy_a = I_ESP_Ky6z_Px_a+ABY*I_ESP_I6z_Px_a;
    Double I_ESP_I6x_D2y_a = I_ESP_K6xy_Py_a+ABY*I_ESP_I6x_Py_a;
    Double I_ESP_I5xy_D2y_a = I_ESP_K5x2y_Py_a+ABY*I_ESP_I5xy_Py_a;
    Double I_ESP_I5xz_D2y_a = I_ESP_K5xyz_Py_a+ABY*I_ESP_I5xz_Py_a;
    Double I_ESP_I4x2y_D2y_a = I_ESP_K4x3y_Py_a+ABY*I_ESP_I4x2y_Py_a;
    Double I_ESP_I4xyz_D2y_a = I_ESP_K4x2yz_Py_a+ABY*I_ESP_I4xyz_Py_a;
    Double I_ESP_I4x2z_D2y_a = I_ESP_K4xy2z_Py_a+ABY*I_ESP_I4x2z_Py_a;
    Double I_ESP_I3x3y_D2y_a = I_ESP_K3x4y_Py_a+ABY*I_ESP_I3x3y_Py_a;
    Double I_ESP_I3x2yz_D2y_a = I_ESP_K3x3yz_Py_a+ABY*I_ESP_I3x2yz_Py_a;
    Double I_ESP_I3xy2z_D2y_a = I_ESP_K3x2y2z_Py_a+ABY*I_ESP_I3xy2z_Py_a;
    Double I_ESP_I3x3z_D2y_a = I_ESP_K3xy3z_Py_a+ABY*I_ESP_I3x3z_Py_a;
    Double I_ESP_I2x4y_D2y_a = I_ESP_K2x5y_Py_a+ABY*I_ESP_I2x4y_Py_a;
    Double I_ESP_I2x3yz_D2y_a = I_ESP_K2x4yz_Py_a+ABY*I_ESP_I2x3yz_Py_a;
    Double I_ESP_I2x2y2z_D2y_a = I_ESP_K2x3y2z_Py_a+ABY*I_ESP_I2x2y2z_Py_a;
    Double I_ESP_I2xy3z_D2y_a = I_ESP_K2x2y3z_Py_a+ABY*I_ESP_I2xy3z_Py_a;
    Double I_ESP_I2x4z_D2y_a = I_ESP_K2xy4z_Py_a+ABY*I_ESP_I2x4z_Py_a;
    Double I_ESP_Ix5y_D2y_a = I_ESP_Kx6y_Py_a+ABY*I_ESP_Ix5y_Py_a;
    Double I_ESP_Ix4yz_D2y_a = I_ESP_Kx5yz_Py_a+ABY*I_ESP_Ix4yz_Py_a;
    Double I_ESP_Ix3y2z_D2y_a = I_ESP_Kx4y2z_Py_a+ABY*I_ESP_Ix3y2z_Py_a;
    Double I_ESP_Ix2y3z_D2y_a = I_ESP_Kx3y3z_Py_a+ABY*I_ESP_Ix2y3z_Py_a;
    Double I_ESP_Ixy4z_D2y_a = I_ESP_Kx2y4z_Py_a+ABY*I_ESP_Ixy4z_Py_a;
    Double I_ESP_Ix5z_D2y_a = I_ESP_Kxy5z_Py_a+ABY*I_ESP_Ix5z_Py_a;
    Double I_ESP_I6y_D2y_a = I_ESP_K7y_Py_a+ABY*I_ESP_I6y_Py_a;
    Double I_ESP_I5yz_D2y_a = I_ESP_K6yz_Py_a+ABY*I_ESP_I5yz_Py_a;
    Double I_ESP_I4y2z_D2y_a = I_ESP_K5y2z_Py_a+ABY*I_ESP_I4y2z_Py_a;
    Double I_ESP_I3y3z_D2y_a = I_ESP_K4y3z_Py_a+ABY*I_ESP_I3y3z_Py_a;
    Double I_ESP_I2y4z_D2y_a = I_ESP_K3y4z_Py_a+ABY*I_ESP_I2y4z_Py_a;
    Double I_ESP_Iy5z_D2y_a = I_ESP_K2y5z_Py_a+ABY*I_ESP_Iy5z_Py_a;
    Double I_ESP_I6z_D2y_a = I_ESP_Ky6z_Py_a+ABY*I_ESP_I6z_Py_a;
    Double I_ESP_I6x_D2z_a = I_ESP_K6xz_Pz_a+ABZ*I_ESP_I6x_Pz_a;
    Double I_ESP_I5xy_D2z_a = I_ESP_K5xyz_Pz_a+ABZ*I_ESP_I5xy_Pz_a;
    Double I_ESP_I5xz_D2z_a = I_ESP_K5x2z_Pz_a+ABZ*I_ESP_I5xz_Pz_a;
    Double I_ESP_I4x2y_D2z_a = I_ESP_K4x2yz_Pz_a+ABZ*I_ESP_I4x2y_Pz_a;
    Double I_ESP_I4xyz_D2z_a = I_ESP_K4xy2z_Pz_a+ABZ*I_ESP_I4xyz_Pz_a;
    Double I_ESP_I4x2z_D2z_a = I_ESP_K4x3z_Pz_a+ABZ*I_ESP_I4x2z_Pz_a;
    Double I_ESP_I3x3y_D2z_a = I_ESP_K3x3yz_Pz_a+ABZ*I_ESP_I3x3y_Pz_a;
    Double I_ESP_I3x2yz_D2z_a = I_ESP_K3x2y2z_Pz_a+ABZ*I_ESP_I3x2yz_Pz_a;
    Double I_ESP_I3xy2z_D2z_a = I_ESP_K3xy3z_Pz_a+ABZ*I_ESP_I3xy2z_Pz_a;
    Double I_ESP_I3x3z_D2z_a = I_ESP_K3x4z_Pz_a+ABZ*I_ESP_I3x3z_Pz_a;
    Double I_ESP_I2x4y_D2z_a = I_ESP_K2x4yz_Pz_a+ABZ*I_ESP_I2x4y_Pz_a;
    Double I_ESP_I2x3yz_D2z_a = I_ESP_K2x3y2z_Pz_a+ABZ*I_ESP_I2x3yz_Pz_a;
    Double I_ESP_I2x2y2z_D2z_a = I_ESP_K2x2y3z_Pz_a+ABZ*I_ESP_I2x2y2z_Pz_a;
    Double I_ESP_I2xy3z_D2z_a = I_ESP_K2xy4z_Pz_a+ABZ*I_ESP_I2xy3z_Pz_a;
    Double I_ESP_I2x4z_D2z_a = I_ESP_K2x5z_Pz_a+ABZ*I_ESP_I2x4z_Pz_a;
    Double I_ESP_Ix5y_D2z_a = I_ESP_Kx5yz_Pz_a+ABZ*I_ESP_Ix5y_Pz_a;
    Double I_ESP_Ix4yz_D2z_a = I_ESP_Kx4y2z_Pz_a+ABZ*I_ESP_Ix4yz_Pz_a;
    Double I_ESP_Ix3y2z_D2z_a = I_ESP_Kx3y3z_Pz_a+ABZ*I_ESP_Ix3y2z_Pz_a;
    Double I_ESP_Ix2y3z_D2z_a = I_ESP_Kx2y4z_Pz_a+ABZ*I_ESP_Ix2y3z_Pz_a;
    Double I_ESP_Ixy4z_D2z_a = I_ESP_Kxy5z_Pz_a+ABZ*I_ESP_Ixy4z_Pz_a;
    Double I_ESP_Ix5z_D2z_a = I_ESP_Kx6z_Pz_a+ABZ*I_ESP_Ix5z_Pz_a;
    Double I_ESP_I6y_D2z_a = I_ESP_K6yz_Pz_a+ABZ*I_ESP_I6y_Pz_a;
    Double I_ESP_I5yz_D2z_a = I_ESP_K5y2z_Pz_a+ABZ*I_ESP_I5yz_Pz_a;
    Double I_ESP_I4y2z_D2z_a = I_ESP_K4y3z_Pz_a+ABZ*I_ESP_I4y2z_Pz_a;
    Double I_ESP_I3y3z_D2z_a = I_ESP_K3y4z_Pz_a+ABZ*I_ESP_I3y3z_Pz_a;
    Double I_ESP_I2y4z_D2z_a = I_ESP_K2y5z_Pz_a+ABZ*I_ESP_I2y4z_Pz_a;
    Double I_ESP_Iy5z_D2z_a = I_ESP_Ky6z_Pz_a+ABZ*I_ESP_Iy5z_Pz_a;
    Double I_ESP_I6z_D2z_a = I_ESP_K7z_Pz_a+ABZ*I_ESP_I6z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     ************************************************************/
    Double I_ESP_H5x_F3x_a = I_ESP_I6x_D2x_a+ABX*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F3x_a = I_ESP_I5xy_D2x_a+ABX*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F3x_a = I_ESP_I5xz_D2x_a+ABX*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F3x_a = I_ESP_I4x2y_D2x_a+ABX*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F3x_a = I_ESP_I4xyz_D2x_a+ABX*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F3x_a = I_ESP_I4x2z_D2x_a+ABX*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F3x_a = I_ESP_I3x3y_D2x_a+ABX*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F3x_a = I_ESP_I3x2yz_D2x_a+ABX*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F3x_a = I_ESP_I3xy2z_D2x_a+ABX*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F3x_a = I_ESP_I3x3z_D2x_a+ABX*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F3x_a = I_ESP_I2x4y_D2x_a+ABX*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F3x_a = I_ESP_I2x3yz_D2x_a+ABX*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F3x_a = I_ESP_I2x2y2z_D2x_a+ABX*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F3x_a = I_ESP_I2xy3z_D2x_a+ABX*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F3x_a = I_ESP_I2x4z_D2x_a+ABX*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F3x_a = I_ESP_Ix5y_D2x_a+ABX*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F3x_a = I_ESP_Ix4yz_D2x_a+ABX*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F3x_a = I_ESP_Ix3y2z_D2x_a+ABX*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F3x_a = I_ESP_Ix2y3z_D2x_a+ABX*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F3x_a = I_ESP_Ixy4z_D2x_a+ABX*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F3x_a = I_ESP_Ix5z_D2x_a+ABX*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_F2xy_a = I_ESP_I5xy_D2x_a+ABY*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F2xy_a = I_ESP_I4x2y_D2x_a+ABY*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F2xy_a = I_ESP_I4xyz_D2x_a+ABY*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F2xy_a = I_ESP_I3x3y_D2x_a+ABY*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F2xy_a = I_ESP_I3x2yz_D2x_a+ABY*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F2xy_a = I_ESP_I3xy2z_D2x_a+ABY*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F2xy_a = I_ESP_I2x4y_D2x_a+ABY*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F2xy_a = I_ESP_I2x3yz_D2x_a+ABY*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F2xy_a = I_ESP_I2x2y2z_D2x_a+ABY*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F2xy_a = I_ESP_I2xy3z_D2x_a+ABY*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F2xy_a = I_ESP_Ix5y_D2x_a+ABY*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F2xy_a = I_ESP_Ix4yz_D2x_a+ABY*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F2xy_a = I_ESP_Ix3y2z_D2x_a+ABY*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F2xy_a = I_ESP_Ix2y3z_D2x_a+ABY*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F2xy_a = I_ESP_Ixy4z_D2x_a+ABY*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F2xy_a = I_ESP_I6y_D2x_a+ABY*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F2xy_a = I_ESP_I5yz_D2x_a+ABY*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F2xy_a = I_ESP_I4y2z_D2x_a+ABY*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F2xy_a = I_ESP_I3y3z_D2x_a+ABY*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F2xy_a = I_ESP_I2y4z_D2x_a+ABY*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F2xy_a = I_ESP_Iy5z_D2x_a+ABY*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_F2xz_a = I_ESP_I5xz_D2x_a+ABZ*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F2xz_a = I_ESP_I4xyz_D2x_a+ABZ*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F2xz_a = I_ESP_I4x2z_D2x_a+ABZ*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F2xz_a = I_ESP_I3x2yz_D2x_a+ABZ*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F2xz_a = I_ESP_I3xy2z_D2x_a+ABZ*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F2xz_a = I_ESP_I3x3z_D2x_a+ABZ*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F2xz_a = I_ESP_I2x3yz_D2x_a+ABZ*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F2xz_a = I_ESP_I2x2y2z_D2x_a+ABZ*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F2xz_a = I_ESP_I2xy3z_D2x_a+ABZ*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F2xz_a = I_ESP_I2x4z_D2x_a+ABZ*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F2xz_a = I_ESP_Ix4yz_D2x_a+ABZ*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F2xz_a = I_ESP_Ix3y2z_D2x_a+ABZ*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F2xz_a = I_ESP_Ix2y3z_D2x_a+ABZ*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F2xz_a = I_ESP_Ixy4z_D2x_a+ABZ*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F2xz_a = I_ESP_Ix5z_D2x_a+ABZ*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F2xz_a = I_ESP_I5yz_D2x_a+ABZ*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F2xz_a = I_ESP_I4y2z_D2x_a+ABZ*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F2xz_a = I_ESP_I3y3z_D2x_a+ABZ*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F2xz_a = I_ESP_I2y4z_D2x_a+ABZ*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F2xz_a = I_ESP_Iy5z_D2x_a+ABZ*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F2xz_a = I_ESP_I6z_D2x_a+ABZ*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_Fx2y_a = I_ESP_I6x_D2y_a+ABX*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_Fx2y_a = I_ESP_I5xy_D2y_a+ABX*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_Fx2y_a = I_ESP_I5xz_D2y_a+ABX*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_Fx2y_a = I_ESP_I4x2y_D2y_a+ABX*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_Fx2y_a = I_ESP_I4xyz_D2y_a+ABX*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_Fx2y_a = I_ESP_I4x2z_D2y_a+ABX*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_Fx2y_a = I_ESP_I3x3y_D2y_a+ABX*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_Fx2y_a = I_ESP_I3x2yz_D2y_a+ABX*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_Fx2y_a = I_ESP_I3xy2z_D2y_a+ABX*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_Fx2y_a = I_ESP_I3x3z_D2y_a+ABX*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_Fx2y_a = I_ESP_I2x4y_D2y_a+ABX*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_Fx2y_a = I_ESP_I2x3yz_D2y_a+ABX*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_Fx2y_a = I_ESP_I2x2y2z_D2y_a+ABX*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_Fx2y_a = I_ESP_I2xy3z_D2y_a+ABX*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_Fx2y_a = I_ESP_I2x4z_D2y_a+ABX*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_Fx2y_a = I_ESP_Ix5y_D2y_a+ABX*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_Fx2y_a = I_ESP_Ix4yz_D2y_a+ABX*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_Fx2y_a = I_ESP_Ix3y2z_D2y_a+ABX*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_Fx2y_a = I_ESP_Ix2y3z_D2y_a+ABX*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_Fx2y_a = I_ESP_Ixy4z_D2y_a+ABX*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_Fx2y_a = I_ESP_Ix5z_D2y_a+ABX*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_Fxyz_a = I_ESP_I5xz_Dxy_a+ABZ*I_ESP_H5x_Dxy_a;
    Double I_ESP_H4xy_Fxyz_a = I_ESP_I4xyz_Dxy_a+ABZ*I_ESP_H4xy_Dxy_a;
    Double I_ESP_H4xz_Fxyz_a = I_ESP_I4x2z_Dxy_a+ABZ*I_ESP_H4xz_Dxy_a;
    Double I_ESP_H3x2y_Fxyz_a = I_ESP_I3x2yz_Dxy_a+ABZ*I_ESP_H3x2y_Dxy_a;
    Double I_ESP_H3xyz_Fxyz_a = I_ESP_I3xy2z_Dxy_a+ABZ*I_ESP_H3xyz_Dxy_a;
    Double I_ESP_H3x2z_Fxyz_a = I_ESP_I3x3z_Dxy_a+ABZ*I_ESP_H3x2z_Dxy_a;
    Double I_ESP_H2x3y_Fxyz_a = I_ESP_I2x3yz_Dxy_a+ABZ*I_ESP_H2x3y_Dxy_a;
    Double I_ESP_H2x2yz_Fxyz_a = I_ESP_I2x2y2z_Dxy_a+ABZ*I_ESP_H2x2yz_Dxy_a;
    Double I_ESP_H2xy2z_Fxyz_a = I_ESP_I2xy3z_Dxy_a+ABZ*I_ESP_H2xy2z_Dxy_a;
    Double I_ESP_H2x3z_Fxyz_a = I_ESP_I2x4z_Dxy_a+ABZ*I_ESP_H2x3z_Dxy_a;
    Double I_ESP_Hx4y_Fxyz_a = I_ESP_Ix4yz_Dxy_a+ABZ*I_ESP_Hx4y_Dxy_a;
    Double I_ESP_Hx3yz_Fxyz_a = I_ESP_Ix3y2z_Dxy_a+ABZ*I_ESP_Hx3yz_Dxy_a;
    Double I_ESP_Hx2y2z_Fxyz_a = I_ESP_Ix2y3z_Dxy_a+ABZ*I_ESP_Hx2y2z_Dxy_a;
    Double I_ESP_Hxy3z_Fxyz_a = I_ESP_Ixy4z_Dxy_a+ABZ*I_ESP_Hxy3z_Dxy_a;
    Double I_ESP_Hx4z_Fxyz_a = I_ESP_Ix5z_Dxy_a+ABZ*I_ESP_Hx4z_Dxy_a;
    Double I_ESP_H5y_Fxyz_a = I_ESP_I5yz_Dxy_a+ABZ*I_ESP_H5y_Dxy_a;
    Double I_ESP_H4yz_Fxyz_a = I_ESP_I4y2z_Dxy_a+ABZ*I_ESP_H4yz_Dxy_a;
    Double I_ESP_H3y2z_Fxyz_a = I_ESP_I3y3z_Dxy_a+ABZ*I_ESP_H3y2z_Dxy_a;
    Double I_ESP_H2y3z_Fxyz_a = I_ESP_I2y4z_Dxy_a+ABZ*I_ESP_H2y3z_Dxy_a;
    Double I_ESP_Hy4z_Fxyz_a = I_ESP_Iy5z_Dxy_a+ABZ*I_ESP_Hy4z_Dxy_a;
    Double I_ESP_H5z_Fxyz_a = I_ESP_I6z_Dxy_a+ABZ*I_ESP_H5z_Dxy_a;
    Double I_ESP_H5x_Fx2z_a = I_ESP_I6x_D2z_a+ABX*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_Fx2z_a = I_ESP_I5xy_D2z_a+ABX*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_Fx2z_a = I_ESP_I5xz_D2z_a+ABX*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_Fx2z_a = I_ESP_I4x2y_D2z_a+ABX*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_Fx2z_a = I_ESP_I4xyz_D2z_a+ABX*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_Fx2z_a = I_ESP_I4x2z_D2z_a+ABX*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_Fx2z_a = I_ESP_I3x3y_D2z_a+ABX*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_Fx2z_a = I_ESP_I3x2yz_D2z_a+ABX*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_Fx2z_a = I_ESP_I3xy2z_D2z_a+ABX*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_Fx2z_a = I_ESP_I3x3z_D2z_a+ABX*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_Fx2z_a = I_ESP_I2x4y_D2z_a+ABX*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_Fx2z_a = I_ESP_I2x3yz_D2z_a+ABX*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_Fx2z_a = I_ESP_I2x2y2z_D2z_a+ABX*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_Fx2z_a = I_ESP_I2xy3z_D2z_a+ABX*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_Fx2z_a = I_ESP_I2x4z_D2z_a+ABX*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_Fx2z_a = I_ESP_Ix5y_D2z_a+ABX*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_Fx2z_a = I_ESP_Ix4yz_D2z_a+ABX*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_Fx2z_a = I_ESP_Ix3y2z_D2z_a+ABX*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_Fx2z_a = I_ESP_Ix2y3z_D2z_a+ABX*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_Fx2z_a = I_ESP_Ixy4z_D2z_a+ABX*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_Fx2z_a = I_ESP_Ix5z_D2z_a+ABX*I_ESP_H5z_D2z_a;
    Double I_ESP_H5x_F3y_a = I_ESP_I5xy_D2y_a+ABY*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_F3y_a = I_ESP_I4x2y_D2y_a+ABY*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_F3y_a = I_ESP_I4xyz_D2y_a+ABY*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_F3y_a = I_ESP_I3x3y_D2y_a+ABY*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_F3y_a = I_ESP_I3x2yz_D2y_a+ABY*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_F3y_a = I_ESP_I3xy2z_D2y_a+ABY*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_F3y_a = I_ESP_I2x4y_D2y_a+ABY*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_F3y_a = I_ESP_I2x3yz_D2y_a+ABY*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_F3y_a = I_ESP_I2x2y2z_D2y_a+ABY*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_F3y_a = I_ESP_I2xy3z_D2y_a+ABY*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_F3y_a = I_ESP_Ix5y_D2y_a+ABY*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_F3y_a = I_ESP_Ix4yz_D2y_a+ABY*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_F3y_a = I_ESP_Ix3y2z_D2y_a+ABY*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_F3y_a = I_ESP_Ix2y3z_D2y_a+ABY*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_F3y_a = I_ESP_Ixy4z_D2y_a+ABY*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_F3y_a = I_ESP_I6y_D2y_a+ABY*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_F3y_a = I_ESP_I5yz_D2y_a+ABY*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_F3y_a = I_ESP_I4y2z_D2y_a+ABY*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_F3y_a = I_ESP_I3y3z_D2y_a+ABY*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_F3y_a = I_ESP_I2y4z_D2y_a+ABY*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_F3y_a = I_ESP_Iy5z_D2y_a+ABY*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_F2yz_a = I_ESP_I5xz_D2y_a+ABZ*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_F2yz_a = I_ESP_I4xyz_D2y_a+ABZ*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_F2yz_a = I_ESP_I4x2z_D2y_a+ABZ*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_F2yz_a = I_ESP_I3x2yz_D2y_a+ABZ*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_F2yz_a = I_ESP_I3xy2z_D2y_a+ABZ*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_F2yz_a = I_ESP_I3x3z_D2y_a+ABZ*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_F2yz_a = I_ESP_I2x3yz_D2y_a+ABZ*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_F2yz_a = I_ESP_I2x2y2z_D2y_a+ABZ*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_F2yz_a = I_ESP_I2xy3z_D2y_a+ABZ*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_F2yz_a = I_ESP_I2x4z_D2y_a+ABZ*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_F2yz_a = I_ESP_Ix4yz_D2y_a+ABZ*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_F2yz_a = I_ESP_Ix3y2z_D2y_a+ABZ*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_F2yz_a = I_ESP_Ix2y3z_D2y_a+ABZ*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_F2yz_a = I_ESP_Ixy4z_D2y_a+ABZ*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_F2yz_a = I_ESP_Ix5z_D2y_a+ABZ*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_F2yz_a = I_ESP_I5yz_D2y_a+ABZ*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_F2yz_a = I_ESP_I4y2z_D2y_a+ABZ*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_F2yz_a = I_ESP_I3y3z_D2y_a+ABZ*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_F2yz_a = I_ESP_I2y4z_D2y_a+ABZ*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_F2yz_a = I_ESP_Iy5z_D2y_a+ABZ*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_F2yz_a = I_ESP_I6z_D2y_a+ABZ*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_Fy2z_a = I_ESP_I5xy_D2z_a+ABY*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_Fy2z_a = I_ESP_I4x2y_D2z_a+ABY*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_Fy2z_a = I_ESP_I4xyz_D2z_a+ABY*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_Fy2z_a = I_ESP_I3x3y_D2z_a+ABY*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_Fy2z_a = I_ESP_I3x2yz_D2z_a+ABY*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_Fy2z_a = I_ESP_I3xy2z_D2z_a+ABY*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_Fy2z_a = I_ESP_I2x4y_D2z_a+ABY*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_Fy2z_a = I_ESP_I2x3yz_D2z_a+ABY*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_Fy2z_a = I_ESP_I2x2y2z_D2z_a+ABY*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_Fy2z_a = I_ESP_I2xy3z_D2z_a+ABY*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_Fy2z_a = I_ESP_Ix5y_D2z_a+ABY*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_Fy2z_a = I_ESP_Ix4yz_D2z_a+ABY*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_Fy2z_a = I_ESP_Ix3y2z_D2z_a+ABY*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_Fy2z_a = I_ESP_Ix2y3z_D2z_a+ABY*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_Fy2z_a = I_ESP_Ixy4z_D2z_a+ABY*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_Fy2z_a = I_ESP_I6y_D2z_a+ABY*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_Fy2z_a = I_ESP_I5yz_D2z_a+ABY*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_Fy2z_a = I_ESP_I4y2z_D2z_a+ABY*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_Fy2z_a = I_ESP_I3y3z_D2z_a+ABY*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_Fy2z_a = I_ESP_I2y4z_D2z_a+ABY*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_Fy2z_a = I_ESP_Iy5z_D2z_a+ABY*I_ESP_H5z_D2z_a;
    Double I_ESP_H5x_F3z_a = I_ESP_I5xz_D2z_a+ABZ*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_F3z_a = I_ESP_I4xyz_D2z_a+ABZ*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_F3z_a = I_ESP_I4x2z_D2z_a+ABZ*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_F3z_a = I_ESP_I3x2yz_D2z_a+ABZ*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_F3z_a = I_ESP_I3xy2z_D2z_a+ABZ*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_F3z_a = I_ESP_I3x3z_D2z_a+ABZ*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_F3z_a = I_ESP_I2x3yz_D2z_a+ABZ*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_F3z_a = I_ESP_I2x2y2z_D2z_a+ABZ*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_F3z_a = I_ESP_I2xy3z_D2z_a+ABZ*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_F3z_a = I_ESP_I2x4z_D2z_a+ABZ*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_F3z_a = I_ESP_Ix4yz_D2z_a+ABZ*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_F3z_a = I_ESP_Ix3y2z_D2z_a+ABZ*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_F3z_a = I_ESP_Ix2y3z_D2z_a+ABZ*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_F3z_a = I_ESP_Ixy4z_D2z_a+ABZ*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_F3z_a = I_ESP_Ix5z_D2z_a+ABZ*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_F3z_a = I_ESP_I5yz_D2z_a+ABZ*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_F3z_a = I_ESP_I4y2z_D2z_a+ABZ*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_F3z_a = I_ESP_I3y3z_D2z_a+ABZ*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_F3z_a = I_ESP_I2y4z_D2z_a+ABZ*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_F3z_a = I_ESP_Iy5z_D2z_a+ABZ*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_F3z_a = I_ESP_I6z_D2z_a+ABZ*I_ESP_H5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*450+0] = 2.0E0*I_ESP_H5x_F3x_a-4*I_ESP_F3x_F3x;
    abcd[iGrid*450+1] = 2.0E0*I_ESP_H4xy_F3x_a-3*I_ESP_F2xy_F3x;
    abcd[iGrid*450+2] = 2.0E0*I_ESP_H4xz_F3x_a-3*I_ESP_F2xz_F3x;
    abcd[iGrid*450+3] = 2.0E0*I_ESP_H3x2y_F3x_a-2*I_ESP_Fx2y_F3x;
    abcd[iGrid*450+4] = 2.0E0*I_ESP_H3xyz_F3x_a-2*I_ESP_Fxyz_F3x;
    abcd[iGrid*450+5] = 2.0E0*I_ESP_H3x2z_F3x_a-2*I_ESP_Fx2z_F3x;
    abcd[iGrid*450+6] = 2.0E0*I_ESP_H2x3y_F3x_a-1*I_ESP_F3y_F3x;
    abcd[iGrid*450+7] = 2.0E0*I_ESP_H2x2yz_F3x_a-1*I_ESP_F2yz_F3x;
    abcd[iGrid*450+8] = 2.0E0*I_ESP_H2xy2z_F3x_a-1*I_ESP_Fy2z_F3x;
    abcd[iGrid*450+9] = 2.0E0*I_ESP_H2x3z_F3x_a-1*I_ESP_F3z_F3x;
    abcd[iGrid*450+10] = 2.0E0*I_ESP_Hx4y_F3x_a;
    abcd[iGrid*450+11] = 2.0E0*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*450+12] = 2.0E0*I_ESP_Hx2y2z_F3x_a;
    abcd[iGrid*450+13] = 2.0E0*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*450+14] = 2.0E0*I_ESP_Hx4z_F3x_a;
    abcd[iGrid*450+15] = 2.0E0*I_ESP_H5x_F2xy_a-4*I_ESP_F3x_F2xy;
    abcd[iGrid*450+16] = 2.0E0*I_ESP_H4xy_F2xy_a-3*I_ESP_F2xy_F2xy;
    abcd[iGrid*450+17] = 2.0E0*I_ESP_H4xz_F2xy_a-3*I_ESP_F2xz_F2xy;
    abcd[iGrid*450+18] = 2.0E0*I_ESP_H3x2y_F2xy_a-2*I_ESP_Fx2y_F2xy;
    abcd[iGrid*450+19] = 2.0E0*I_ESP_H3xyz_F2xy_a-2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*450+20] = 2.0E0*I_ESP_H3x2z_F2xy_a-2*I_ESP_Fx2z_F2xy;
    abcd[iGrid*450+21] = 2.0E0*I_ESP_H2x3y_F2xy_a-1*I_ESP_F3y_F2xy;
    abcd[iGrid*450+22] = 2.0E0*I_ESP_H2x2yz_F2xy_a-1*I_ESP_F2yz_F2xy;
    abcd[iGrid*450+23] = 2.0E0*I_ESP_H2xy2z_F2xy_a-1*I_ESP_Fy2z_F2xy;
    abcd[iGrid*450+24] = 2.0E0*I_ESP_H2x3z_F2xy_a-1*I_ESP_F3z_F2xy;
    abcd[iGrid*450+25] = 2.0E0*I_ESP_Hx4y_F2xy_a;
    abcd[iGrid*450+26] = 2.0E0*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*450+27] = 2.0E0*I_ESP_Hx2y2z_F2xy_a;
    abcd[iGrid*450+28] = 2.0E0*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*450+29] = 2.0E0*I_ESP_Hx4z_F2xy_a;
    abcd[iGrid*450+30] = 2.0E0*I_ESP_H5x_F2xz_a-4*I_ESP_F3x_F2xz;
    abcd[iGrid*450+31] = 2.0E0*I_ESP_H4xy_F2xz_a-3*I_ESP_F2xy_F2xz;
    abcd[iGrid*450+32] = 2.0E0*I_ESP_H4xz_F2xz_a-3*I_ESP_F2xz_F2xz;
    abcd[iGrid*450+33] = 2.0E0*I_ESP_H3x2y_F2xz_a-2*I_ESP_Fx2y_F2xz;
    abcd[iGrid*450+34] = 2.0E0*I_ESP_H3xyz_F2xz_a-2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*450+35] = 2.0E0*I_ESP_H3x2z_F2xz_a-2*I_ESP_Fx2z_F2xz;
    abcd[iGrid*450+36] = 2.0E0*I_ESP_H2x3y_F2xz_a-1*I_ESP_F3y_F2xz;
    abcd[iGrid*450+37] = 2.0E0*I_ESP_H2x2yz_F2xz_a-1*I_ESP_F2yz_F2xz;
    abcd[iGrid*450+38] = 2.0E0*I_ESP_H2xy2z_F2xz_a-1*I_ESP_Fy2z_F2xz;
    abcd[iGrid*450+39] = 2.0E0*I_ESP_H2x3z_F2xz_a-1*I_ESP_F3z_F2xz;
    abcd[iGrid*450+40] = 2.0E0*I_ESP_Hx4y_F2xz_a;
    abcd[iGrid*450+41] = 2.0E0*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*450+42] = 2.0E0*I_ESP_Hx2y2z_F2xz_a;
    abcd[iGrid*450+43] = 2.0E0*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*450+44] = 2.0E0*I_ESP_Hx4z_F2xz_a;
    abcd[iGrid*450+45] = 2.0E0*I_ESP_H5x_Fx2y_a-4*I_ESP_F3x_Fx2y;
    abcd[iGrid*450+46] = 2.0E0*I_ESP_H4xy_Fx2y_a-3*I_ESP_F2xy_Fx2y;
    abcd[iGrid*450+47] = 2.0E0*I_ESP_H4xz_Fx2y_a-3*I_ESP_F2xz_Fx2y;
    abcd[iGrid*450+48] = 2.0E0*I_ESP_H3x2y_Fx2y_a-2*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*450+49] = 2.0E0*I_ESP_H3xyz_Fx2y_a-2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*450+50] = 2.0E0*I_ESP_H3x2z_Fx2y_a-2*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*450+51] = 2.0E0*I_ESP_H2x3y_Fx2y_a-1*I_ESP_F3y_Fx2y;
    abcd[iGrid*450+52] = 2.0E0*I_ESP_H2x2yz_Fx2y_a-1*I_ESP_F2yz_Fx2y;
    abcd[iGrid*450+53] = 2.0E0*I_ESP_H2xy2z_Fx2y_a-1*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*450+54] = 2.0E0*I_ESP_H2x3z_Fx2y_a-1*I_ESP_F3z_Fx2y;
    abcd[iGrid*450+55] = 2.0E0*I_ESP_Hx4y_Fx2y_a;
    abcd[iGrid*450+56] = 2.0E0*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*450+57] = 2.0E0*I_ESP_Hx2y2z_Fx2y_a;
    abcd[iGrid*450+58] = 2.0E0*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*450+59] = 2.0E0*I_ESP_Hx4z_Fx2y_a;
    abcd[iGrid*450+60] = 2.0E0*I_ESP_H5x_Fxyz_a-4*I_ESP_F3x_Fxyz;
    abcd[iGrid*450+61] = 2.0E0*I_ESP_H4xy_Fxyz_a-3*I_ESP_F2xy_Fxyz;
    abcd[iGrid*450+62] = 2.0E0*I_ESP_H4xz_Fxyz_a-3*I_ESP_F2xz_Fxyz;
    abcd[iGrid*450+63] = 2.0E0*I_ESP_H3x2y_Fxyz_a-2*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*450+64] = 2.0E0*I_ESP_H3xyz_Fxyz_a-2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*450+65] = 2.0E0*I_ESP_H3x2z_Fxyz_a-2*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*450+66] = 2.0E0*I_ESP_H2x3y_Fxyz_a-1*I_ESP_F3y_Fxyz;
    abcd[iGrid*450+67] = 2.0E0*I_ESP_H2x2yz_Fxyz_a-1*I_ESP_F2yz_Fxyz;
    abcd[iGrid*450+68] = 2.0E0*I_ESP_H2xy2z_Fxyz_a-1*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*450+69] = 2.0E0*I_ESP_H2x3z_Fxyz_a-1*I_ESP_F3z_Fxyz;
    abcd[iGrid*450+70] = 2.0E0*I_ESP_Hx4y_Fxyz_a;
    abcd[iGrid*450+71] = 2.0E0*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*450+72] = 2.0E0*I_ESP_Hx2y2z_Fxyz_a;
    abcd[iGrid*450+73] = 2.0E0*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*450+74] = 2.0E0*I_ESP_Hx4z_Fxyz_a;
    abcd[iGrid*450+75] = 2.0E0*I_ESP_H5x_Fx2z_a-4*I_ESP_F3x_Fx2z;
    abcd[iGrid*450+76] = 2.0E0*I_ESP_H4xy_Fx2z_a-3*I_ESP_F2xy_Fx2z;
    abcd[iGrid*450+77] = 2.0E0*I_ESP_H4xz_Fx2z_a-3*I_ESP_F2xz_Fx2z;
    abcd[iGrid*450+78] = 2.0E0*I_ESP_H3x2y_Fx2z_a-2*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*450+79] = 2.0E0*I_ESP_H3xyz_Fx2z_a-2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*450+80] = 2.0E0*I_ESP_H3x2z_Fx2z_a-2*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*450+81] = 2.0E0*I_ESP_H2x3y_Fx2z_a-1*I_ESP_F3y_Fx2z;
    abcd[iGrid*450+82] = 2.0E0*I_ESP_H2x2yz_Fx2z_a-1*I_ESP_F2yz_Fx2z;
    abcd[iGrid*450+83] = 2.0E0*I_ESP_H2xy2z_Fx2z_a-1*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*450+84] = 2.0E0*I_ESP_H2x3z_Fx2z_a-1*I_ESP_F3z_Fx2z;
    abcd[iGrid*450+85] = 2.0E0*I_ESP_Hx4y_Fx2z_a;
    abcd[iGrid*450+86] = 2.0E0*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*450+87] = 2.0E0*I_ESP_Hx2y2z_Fx2z_a;
    abcd[iGrid*450+88] = 2.0E0*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*450+89] = 2.0E0*I_ESP_Hx4z_Fx2z_a;
    abcd[iGrid*450+90] = 2.0E0*I_ESP_H5x_F3y_a-4*I_ESP_F3x_F3y;
    abcd[iGrid*450+91] = 2.0E0*I_ESP_H4xy_F3y_a-3*I_ESP_F2xy_F3y;
    abcd[iGrid*450+92] = 2.0E0*I_ESP_H4xz_F3y_a-3*I_ESP_F2xz_F3y;
    abcd[iGrid*450+93] = 2.0E0*I_ESP_H3x2y_F3y_a-2*I_ESP_Fx2y_F3y;
    abcd[iGrid*450+94] = 2.0E0*I_ESP_H3xyz_F3y_a-2*I_ESP_Fxyz_F3y;
    abcd[iGrid*450+95] = 2.0E0*I_ESP_H3x2z_F3y_a-2*I_ESP_Fx2z_F3y;
    abcd[iGrid*450+96] = 2.0E0*I_ESP_H2x3y_F3y_a-1*I_ESP_F3y_F3y;
    abcd[iGrid*450+97] = 2.0E0*I_ESP_H2x2yz_F3y_a-1*I_ESP_F2yz_F3y;
    abcd[iGrid*450+98] = 2.0E0*I_ESP_H2xy2z_F3y_a-1*I_ESP_Fy2z_F3y;
    abcd[iGrid*450+99] = 2.0E0*I_ESP_H2x3z_F3y_a-1*I_ESP_F3z_F3y;
    abcd[iGrid*450+100] = 2.0E0*I_ESP_Hx4y_F3y_a;
    abcd[iGrid*450+101] = 2.0E0*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*450+102] = 2.0E0*I_ESP_Hx2y2z_F3y_a;
    abcd[iGrid*450+103] = 2.0E0*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*450+104] = 2.0E0*I_ESP_Hx4z_F3y_a;
    abcd[iGrid*450+105] = 2.0E0*I_ESP_H5x_F2yz_a-4*I_ESP_F3x_F2yz;
    abcd[iGrid*450+106] = 2.0E0*I_ESP_H4xy_F2yz_a-3*I_ESP_F2xy_F2yz;
    abcd[iGrid*450+107] = 2.0E0*I_ESP_H4xz_F2yz_a-3*I_ESP_F2xz_F2yz;
    abcd[iGrid*450+108] = 2.0E0*I_ESP_H3x2y_F2yz_a-2*I_ESP_Fx2y_F2yz;
    abcd[iGrid*450+109] = 2.0E0*I_ESP_H3xyz_F2yz_a-2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*450+110] = 2.0E0*I_ESP_H3x2z_F2yz_a-2*I_ESP_Fx2z_F2yz;
    abcd[iGrid*450+111] = 2.0E0*I_ESP_H2x3y_F2yz_a-1*I_ESP_F3y_F2yz;
    abcd[iGrid*450+112] = 2.0E0*I_ESP_H2x2yz_F2yz_a-1*I_ESP_F2yz_F2yz;
    abcd[iGrid*450+113] = 2.0E0*I_ESP_H2xy2z_F2yz_a-1*I_ESP_Fy2z_F2yz;
    abcd[iGrid*450+114] = 2.0E0*I_ESP_H2x3z_F2yz_a-1*I_ESP_F3z_F2yz;
    abcd[iGrid*450+115] = 2.0E0*I_ESP_Hx4y_F2yz_a;
    abcd[iGrid*450+116] = 2.0E0*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*450+117] = 2.0E0*I_ESP_Hx2y2z_F2yz_a;
    abcd[iGrid*450+118] = 2.0E0*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*450+119] = 2.0E0*I_ESP_Hx4z_F2yz_a;
    abcd[iGrid*450+120] = 2.0E0*I_ESP_H5x_Fy2z_a-4*I_ESP_F3x_Fy2z;
    abcd[iGrid*450+121] = 2.0E0*I_ESP_H4xy_Fy2z_a-3*I_ESP_F2xy_Fy2z;
    abcd[iGrid*450+122] = 2.0E0*I_ESP_H4xz_Fy2z_a-3*I_ESP_F2xz_Fy2z;
    abcd[iGrid*450+123] = 2.0E0*I_ESP_H3x2y_Fy2z_a-2*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*450+124] = 2.0E0*I_ESP_H3xyz_Fy2z_a-2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*450+125] = 2.0E0*I_ESP_H3x2z_Fy2z_a-2*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*450+126] = 2.0E0*I_ESP_H2x3y_Fy2z_a-1*I_ESP_F3y_Fy2z;
    abcd[iGrid*450+127] = 2.0E0*I_ESP_H2x2yz_Fy2z_a-1*I_ESP_F2yz_Fy2z;
    abcd[iGrid*450+128] = 2.0E0*I_ESP_H2xy2z_Fy2z_a-1*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*450+129] = 2.0E0*I_ESP_H2x3z_Fy2z_a-1*I_ESP_F3z_Fy2z;
    abcd[iGrid*450+130] = 2.0E0*I_ESP_Hx4y_Fy2z_a;
    abcd[iGrid*450+131] = 2.0E0*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*450+132] = 2.0E0*I_ESP_Hx2y2z_Fy2z_a;
    abcd[iGrid*450+133] = 2.0E0*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*450+134] = 2.0E0*I_ESP_Hx4z_Fy2z_a;
    abcd[iGrid*450+135] = 2.0E0*I_ESP_H5x_F3z_a-4*I_ESP_F3x_F3z;
    abcd[iGrid*450+136] = 2.0E0*I_ESP_H4xy_F3z_a-3*I_ESP_F2xy_F3z;
    abcd[iGrid*450+137] = 2.0E0*I_ESP_H4xz_F3z_a-3*I_ESP_F2xz_F3z;
    abcd[iGrid*450+138] = 2.0E0*I_ESP_H3x2y_F3z_a-2*I_ESP_Fx2y_F3z;
    abcd[iGrid*450+139] = 2.0E0*I_ESP_H3xyz_F3z_a-2*I_ESP_Fxyz_F3z;
    abcd[iGrid*450+140] = 2.0E0*I_ESP_H3x2z_F3z_a-2*I_ESP_Fx2z_F3z;
    abcd[iGrid*450+141] = 2.0E0*I_ESP_H2x3y_F3z_a-1*I_ESP_F3y_F3z;
    abcd[iGrid*450+142] = 2.0E0*I_ESP_H2x2yz_F3z_a-1*I_ESP_F2yz_F3z;
    abcd[iGrid*450+143] = 2.0E0*I_ESP_H2xy2z_F3z_a-1*I_ESP_Fy2z_F3z;
    abcd[iGrid*450+144] = 2.0E0*I_ESP_H2x3z_F3z_a-1*I_ESP_F3z_F3z;
    abcd[iGrid*450+145] = 2.0E0*I_ESP_Hx4y_F3z_a;
    abcd[iGrid*450+146] = 2.0E0*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*450+147] = 2.0E0*I_ESP_Hx2y2z_F3z_a;
    abcd[iGrid*450+148] = 2.0E0*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*450+149] = 2.0E0*I_ESP_Hx4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*450+150] = 2.0E0*I_ESP_H4xy_F3x_a;
    abcd[iGrid*450+151] = 2.0E0*I_ESP_H3x2y_F3x_a-1*I_ESP_F3x_F3x;
    abcd[iGrid*450+152] = 2.0E0*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*450+153] = 2.0E0*I_ESP_H2x3y_F3x_a-2*I_ESP_F2xy_F3x;
    abcd[iGrid*450+154] = 2.0E0*I_ESP_H2x2yz_F3x_a-1*I_ESP_F2xz_F3x;
    abcd[iGrid*450+155] = 2.0E0*I_ESP_H2xy2z_F3x_a;
    abcd[iGrid*450+156] = 2.0E0*I_ESP_Hx4y_F3x_a-3*I_ESP_Fx2y_F3x;
    abcd[iGrid*450+157] = 2.0E0*I_ESP_Hx3yz_F3x_a-2*I_ESP_Fxyz_F3x;
    abcd[iGrid*450+158] = 2.0E0*I_ESP_Hx2y2z_F3x_a-1*I_ESP_Fx2z_F3x;
    abcd[iGrid*450+159] = 2.0E0*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*450+160] = 2.0E0*I_ESP_H5y_F3x_a-4*I_ESP_F3y_F3x;
    abcd[iGrid*450+161] = 2.0E0*I_ESP_H4yz_F3x_a-3*I_ESP_F2yz_F3x;
    abcd[iGrid*450+162] = 2.0E0*I_ESP_H3y2z_F3x_a-2*I_ESP_Fy2z_F3x;
    abcd[iGrid*450+163] = 2.0E0*I_ESP_H2y3z_F3x_a-1*I_ESP_F3z_F3x;
    abcd[iGrid*450+164] = 2.0E0*I_ESP_Hy4z_F3x_a;
    abcd[iGrid*450+165] = 2.0E0*I_ESP_H4xy_F2xy_a;
    abcd[iGrid*450+166] = 2.0E0*I_ESP_H3x2y_F2xy_a-1*I_ESP_F3x_F2xy;
    abcd[iGrid*450+167] = 2.0E0*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*450+168] = 2.0E0*I_ESP_H2x3y_F2xy_a-2*I_ESP_F2xy_F2xy;
    abcd[iGrid*450+169] = 2.0E0*I_ESP_H2x2yz_F2xy_a-1*I_ESP_F2xz_F2xy;
    abcd[iGrid*450+170] = 2.0E0*I_ESP_H2xy2z_F2xy_a;
    abcd[iGrid*450+171] = 2.0E0*I_ESP_Hx4y_F2xy_a-3*I_ESP_Fx2y_F2xy;
    abcd[iGrid*450+172] = 2.0E0*I_ESP_Hx3yz_F2xy_a-2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*450+173] = 2.0E0*I_ESP_Hx2y2z_F2xy_a-1*I_ESP_Fx2z_F2xy;
    abcd[iGrid*450+174] = 2.0E0*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*450+175] = 2.0E0*I_ESP_H5y_F2xy_a-4*I_ESP_F3y_F2xy;
    abcd[iGrid*450+176] = 2.0E0*I_ESP_H4yz_F2xy_a-3*I_ESP_F2yz_F2xy;
    abcd[iGrid*450+177] = 2.0E0*I_ESP_H3y2z_F2xy_a-2*I_ESP_Fy2z_F2xy;
    abcd[iGrid*450+178] = 2.0E0*I_ESP_H2y3z_F2xy_a-1*I_ESP_F3z_F2xy;
    abcd[iGrid*450+179] = 2.0E0*I_ESP_Hy4z_F2xy_a;
    abcd[iGrid*450+180] = 2.0E0*I_ESP_H4xy_F2xz_a;
    abcd[iGrid*450+181] = 2.0E0*I_ESP_H3x2y_F2xz_a-1*I_ESP_F3x_F2xz;
    abcd[iGrid*450+182] = 2.0E0*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*450+183] = 2.0E0*I_ESP_H2x3y_F2xz_a-2*I_ESP_F2xy_F2xz;
    abcd[iGrid*450+184] = 2.0E0*I_ESP_H2x2yz_F2xz_a-1*I_ESP_F2xz_F2xz;
    abcd[iGrid*450+185] = 2.0E0*I_ESP_H2xy2z_F2xz_a;
    abcd[iGrid*450+186] = 2.0E0*I_ESP_Hx4y_F2xz_a-3*I_ESP_Fx2y_F2xz;
    abcd[iGrid*450+187] = 2.0E0*I_ESP_Hx3yz_F2xz_a-2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*450+188] = 2.0E0*I_ESP_Hx2y2z_F2xz_a-1*I_ESP_Fx2z_F2xz;
    abcd[iGrid*450+189] = 2.0E0*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*450+190] = 2.0E0*I_ESP_H5y_F2xz_a-4*I_ESP_F3y_F2xz;
    abcd[iGrid*450+191] = 2.0E0*I_ESP_H4yz_F2xz_a-3*I_ESP_F2yz_F2xz;
    abcd[iGrid*450+192] = 2.0E0*I_ESP_H3y2z_F2xz_a-2*I_ESP_Fy2z_F2xz;
    abcd[iGrid*450+193] = 2.0E0*I_ESP_H2y3z_F2xz_a-1*I_ESP_F3z_F2xz;
    abcd[iGrid*450+194] = 2.0E0*I_ESP_Hy4z_F2xz_a;
    abcd[iGrid*450+195] = 2.0E0*I_ESP_H4xy_Fx2y_a;
    abcd[iGrid*450+196] = 2.0E0*I_ESP_H3x2y_Fx2y_a-1*I_ESP_F3x_Fx2y;
    abcd[iGrid*450+197] = 2.0E0*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*450+198] = 2.0E0*I_ESP_H2x3y_Fx2y_a-2*I_ESP_F2xy_Fx2y;
    abcd[iGrid*450+199] = 2.0E0*I_ESP_H2x2yz_Fx2y_a-1*I_ESP_F2xz_Fx2y;
    abcd[iGrid*450+200] = 2.0E0*I_ESP_H2xy2z_Fx2y_a;
    abcd[iGrid*450+201] = 2.0E0*I_ESP_Hx4y_Fx2y_a-3*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*450+202] = 2.0E0*I_ESP_Hx3yz_Fx2y_a-2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*450+203] = 2.0E0*I_ESP_Hx2y2z_Fx2y_a-1*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*450+204] = 2.0E0*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*450+205] = 2.0E0*I_ESP_H5y_Fx2y_a-4*I_ESP_F3y_Fx2y;
    abcd[iGrid*450+206] = 2.0E0*I_ESP_H4yz_Fx2y_a-3*I_ESP_F2yz_Fx2y;
    abcd[iGrid*450+207] = 2.0E0*I_ESP_H3y2z_Fx2y_a-2*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*450+208] = 2.0E0*I_ESP_H2y3z_Fx2y_a-1*I_ESP_F3z_Fx2y;
    abcd[iGrid*450+209] = 2.0E0*I_ESP_Hy4z_Fx2y_a;
    abcd[iGrid*450+210] = 2.0E0*I_ESP_H4xy_Fxyz_a;
    abcd[iGrid*450+211] = 2.0E0*I_ESP_H3x2y_Fxyz_a-1*I_ESP_F3x_Fxyz;
    abcd[iGrid*450+212] = 2.0E0*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*450+213] = 2.0E0*I_ESP_H2x3y_Fxyz_a-2*I_ESP_F2xy_Fxyz;
    abcd[iGrid*450+214] = 2.0E0*I_ESP_H2x2yz_Fxyz_a-1*I_ESP_F2xz_Fxyz;
    abcd[iGrid*450+215] = 2.0E0*I_ESP_H2xy2z_Fxyz_a;
    abcd[iGrid*450+216] = 2.0E0*I_ESP_Hx4y_Fxyz_a-3*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*450+217] = 2.0E0*I_ESP_Hx3yz_Fxyz_a-2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*450+218] = 2.0E0*I_ESP_Hx2y2z_Fxyz_a-1*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*450+219] = 2.0E0*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*450+220] = 2.0E0*I_ESP_H5y_Fxyz_a-4*I_ESP_F3y_Fxyz;
    abcd[iGrid*450+221] = 2.0E0*I_ESP_H4yz_Fxyz_a-3*I_ESP_F2yz_Fxyz;
    abcd[iGrid*450+222] = 2.0E0*I_ESP_H3y2z_Fxyz_a-2*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*450+223] = 2.0E0*I_ESP_H2y3z_Fxyz_a-1*I_ESP_F3z_Fxyz;
    abcd[iGrid*450+224] = 2.0E0*I_ESP_Hy4z_Fxyz_a;
    abcd[iGrid*450+225] = 2.0E0*I_ESP_H4xy_Fx2z_a;
    abcd[iGrid*450+226] = 2.0E0*I_ESP_H3x2y_Fx2z_a-1*I_ESP_F3x_Fx2z;
    abcd[iGrid*450+227] = 2.0E0*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*450+228] = 2.0E0*I_ESP_H2x3y_Fx2z_a-2*I_ESP_F2xy_Fx2z;
    abcd[iGrid*450+229] = 2.0E0*I_ESP_H2x2yz_Fx2z_a-1*I_ESP_F2xz_Fx2z;
    abcd[iGrid*450+230] = 2.0E0*I_ESP_H2xy2z_Fx2z_a;
    abcd[iGrid*450+231] = 2.0E0*I_ESP_Hx4y_Fx2z_a-3*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*450+232] = 2.0E0*I_ESP_Hx3yz_Fx2z_a-2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*450+233] = 2.0E0*I_ESP_Hx2y2z_Fx2z_a-1*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*450+234] = 2.0E0*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*450+235] = 2.0E0*I_ESP_H5y_Fx2z_a-4*I_ESP_F3y_Fx2z;
    abcd[iGrid*450+236] = 2.0E0*I_ESP_H4yz_Fx2z_a-3*I_ESP_F2yz_Fx2z;
    abcd[iGrid*450+237] = 2.0E0*I_ESP_H3y2z_Fx2z_a-2*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*450+238] = 2.0E0*I_ESP_H2y3z_Fx2z_a-1*I_ESP_F3z_Fx2z;
    abcd[iGrid*450+239] = 2.0E0*I_ESP_Hy4z_Fx2z_a;
    abcd[iGrid*450+240] = 2.0E0*I_ESP_H4xy_F3y_a;
    abcd[iGrid*450+241] = 2.0E0*I_ESP_H3x2y_F3y_a-1*I_ESP_F3x_F3y;
    abcd[iGrid*450+242] = 2.0E0*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*450+243] = 2.0E0*I_ESP_H2x3y_F3y_a-2*I_ESP_F2xy_F3y;
    abcd[iGrid*450+244] = 2.0E0*I_ESP_H2x2yz_F3y_a-1*I_ESP_F2xz_F3y;
    abcd[iGrid*450+245] = 2.0E0*I_ESP_H2xy2z_F3y_a;
    abcd[iGrid*450+246] = 2.0E0*I_ESP_Hx4y_F3y_a-3*I_ESP_Fx2y_F3y;
    abcd[iGrid*450+247] = 2.0E0*I_ESP_Hx3yz_F3y_a-2*I_ESP_Fxyz_F3y;
    abcd[iGrid*450+248] = 2.0E0*I_ESP_Hx2y2z_F3y_a-1*I_ESP_Fx2z_F3y;
    abcd[iGrid*450+249] = 2.0E0*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*450+250] = 2.0E0*I_ESP_H5y_F3y_a-4*I_ESP_F3y_F3y;
    abcd[iGrid*450+251] = 2.0E0*I_ESP_H4yz_F3y_a-3*I_ESP_F2yz_F3y;
    abcd[iGrid*450+252] = 2.0E0*I_ESP_H3y2z_F3y_a-2*I_ESP_Fy2z_F3y;
    abcd[iGrid*450+253] = 2.0E0*I_ESP_H2y3z_F3y_a-1*I_ESP_F3z_F3y;
    abcd[iGrid*450+254] = 2.0E0*I_ESP_Hy4z_F3y_a;
    abcd[iGrid*450+255] = 2.0E0*I_ESP_H4xy_F2yz_a;
    abcd[iGrid*450+256] = 2.0E0*I_ESP_H3x2y_F2yz_a-1*I_ESP_F3x_F2yz;
    abcd[iGrid*450+257] = 2.0E0*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*450+258] = 2.0E0*I_ESP_H2x3y_F2yz_a-2*I_ESP_F2xy_F2yz;
    abcd[iGrid*450+259] = 2.0E0*I_ESP_H2x2yz_F2yz_a-1*I_ESP_F2xz_F2yz;
    abcd[iGrid*450+260] = 2.0E0*I_ESP_H2xy2z_F2yz_a;
    abcd[iGrid*450+261] = 2.0E0*I_ESP_Hx4y_F2yz_a-3*I_ESP_Fx2y_F2yz;
    abcd[iGrid*450+262] = 2.0E0*I_ESP_Hx3yz_F2yz_a-2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*450+263] = 2.0E0*I_ESP_Hx2y2z_F2yz_a-1*I_ESP_Fx2z_F2yz;
    abcd[iGrid*450+264] = 2.0E0*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*450+265] = 2.0E0*I_ESP_H5y_F2yz_a-4*I_ESP_F3y_F2yz;
    abcd[iGrid*450+266] = 2.0E0*I_ESP_H4yz_F2yz_a-3*I_ESP_F2yz_F2yz;
    abcd[iGrid*450+267] = 2.0E0*I_ESP_H3y2z_F2yz_a-2*I_ESP_Fy2z_F2yz;
    abcd[iGrid*450+268] = 2.0E0*I_ESP_H2y3z_F2yz_a-1*I_ESP_F3z_F2yz;
    abcd[iGrid*450+269] = 2.0E0*I_ESP_Hy4z_F2yz_a;
    abcd[iGrid*450+270] = 2.0E0*I_ESP_H4xy_Fy2z_a;
    abcd[iGrid*450+271] = 2.0E0*I_ESP_H3x2y_Fy2z_a-1*I_ESP_F3x_Fy2z;
    abcd[iGrid*450+272] = 2.0E0*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*450+273] = 2.0E0*I_ESP_H2x3y_Fy2z_a-2*I_ESP_F2xy_Fy2z;
    abcd[iGrid*450+274] = 2.0E0*I_ESP_H2x2yz_Fy2z_a-1*I_ESP_F2xz_Fy2z;
    abcd[iGrid*450+275] = 2.0E0*I_ESP_H2xy2z_Fy2z_a;
    abcd[iGrid*450+276] = 2.0E0*I_ESP_Hx4y_Fy2z_a-3*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*450+277] = 2.0E0*I_ESP_Hx3yz_Fy2z_a-2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*450+278] = 2.0E0*I_ESP_Hx2y2z_Fy2z_a-1*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*450+279] = 2.0E0*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*450+280] = 2.0E0*I_ESP_H5y_Fy2z_a-4*I_ESP_F3y_Fy2z;
    abcd[iGrid*450+281] = 2.0E0*I_ESP_H4yz_Fy2z_a-3*I_ESP_F2yz_Fy2z;
    abcd[iGrid*450+282] = 2.0E0*I_ESP_H3y2z_Fy2z_a-2*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*450+283] = 2.0E0*I_ESP_H2y3z_Fy2z_a-1*I_ESP_F3z_Fy2z;
    abcd[iGrid*450+284] = 2.0E0*I_ESP_Hy4z_Fy2z_a;
    abcd[iGrid*450+285] = 2.0E0*I_ESP_H4xy_F3z_a;
    abcd[iGrid*450+286] = 2.0E0*I_ESP_H3x2y_F3z_a-1*I_ESP_F3x_F3z;
    abcd[iGrid*450+287] = 2.0E0*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*450+288] = 2.0E0*I_ESP_H2x3y_F3z_a-2*I_ESP_F2xy_F3z;
    abcd[iGrid*450+289] = 2.0E0*I_ESP_H2x2yz_F3z_a-1*I_ESP_F2xz_F3z;
    abcd[iGrid*450+290] = 2.0E0*I_ESP_H2xy2z_F3z_a;
    abcd[iGrid*450+291] = 2.0E0*I_ESP_Hx4y_F3z_a-3*I_ESP_Fx2y_F3z;
    abcd[iGrid*450+292] = 2.0E0*I_ESP_Hx3yz_F3z_a-2*I_ESP_Fxyz_F3z;
    abcd[iGrid*450+293] = 2.0E0*I_ESP_Hx2y2z_F3z_a-1*I_ESP_Fx2z_F3z;
    abcd[iGrid*450+294] = 2.0E0*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*450+295] = 2.0E0*I_ESP_H5y_F3z_a-4*I_ESP_F3y_F3z;
    abcd[iGrid*450+296] = 2.0E0*I_ESP_H4yz_F3z_a-3*I_ESP_F2yz_F3z;
    abcd[iGrid*450+297] = 2.0E0*I_ESP_H3y2z_F3z_a-2*I_ESP_Fy2z_F3z;
    abcd[iGrid*450+298] = 2.0E0*I_ESP_H2y3z_F3z_a-1*I_ESP_F3z_F3z;
    abcd[iGrid*450+299] = 2.0E0*I_ESP_Hy4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*450+300] = 2.0E0*I_ESP_H4xz_F3x_a;
    abcd[iGrid*450+301] = 2.0E0*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*450+302] = 2.0E0*I_ESP_H3x2z_F3x_a-1*I_ESP_F3x_F3x;
    abcd[iGrid*450+303] = 2.0E0*I_ESP_H2x2yz_F3x_a;
    abcd[iGrid*450+304] = 2.0E0*I_ESP_H2xy2z_F3x_a-1*I_ESP_F2xy_F3x;
    abcd[iGrid*450+305] = 2.0E0*I_ESP_H2x3z_F3x_a-2*I_ESP_F2xz_F3x;
    abcd[iGrid*450+306] = 2.0E0*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*450+307] = 2.0E0*I_ESP_Hx2y2z_F3x_a-1*I_ESP_Fx2y_F3x;
    abcd[iGrid*450+308] = 2.0E0*I_ESP_Hxy3z_F3x_a-2*I_ESP_Fxyz_F3x;
    abcd[iGrid*450+309] = 2.0E0*I_ESP_Hx4z_F3x_a-3*I_ESP_Fx2z_F3x;
    abcd[iGrid*450+310] = 2.0E0*I_ESP_H4yz_F3x_a;
    abcd[iGrid*450+311] = 2.0E0*I_ESP_H3y2z_F3x_a-1*I_ESP_F3y_F3x;
    abcd[iGrid*450+312] = 2.0E0*I_ESP_H2y3z_F3x_a-2*I_ESP_F2yz_F3x;
    abcd[iGrid*450+313] = 2.0E0*I_ESP_Hy4z_F3x_a-3*I_ESP_Fy2z_F3x;
    abcd[iGrid*450+314] = 2.0E0*I_ESP_H5z_F3x_a-4*I_ESP_F3z_F3x;
    abcd[iGrid*450+315] = 2.0E0*I_ESP_H4xz_F2xy_a;
    abcd[iGrid*450+316] = 2.0E0*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*450+317] = 2.0E0*I_ESP_H3x2z_F2xy_a-1*I_ESP_F3x_F2xy;
    abcd[iGrid*450+318] = 2.0E0*I_ESP_H2x2yz_F2xy_a;
    abcd[iGrid*450+319] = 2.0E0*I_ESP_H2xy2z_F2xy_a-1*I_ESP_F2xy_F2xy;
    abcd[iGrid*450+320] = 2.0E0*I_ESP_H2x3z_F2xy_a-2*I_ESP_F2xz_F2xy;
    abcd[iGrid*450+321] = 2.0E0*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*450+322] = 2.0E0*I_ESP_Hx2y2z_F2xy_a-1*I_ESP_Fx2y_F2xy;
    abcd[iGrid*450+323] = 2.0E0*I_ESP_Hxy3z_F2xy_a-2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*450+324] = 2.0E0*I_ESP_Hx4z_F2xy_a-3*I_ESP_Fx2z_F2xy;
    abcd[iGrid*450+325] = 2.0E0*I_ESP_H4yz_F2xy_a;
    abcd[iGrid*450+326] = 2.0E0*I_ESP_H3y2z_F2xy_a-1*I_ESP_F3y_F2xy;
    abcd[iGrid*450+327] = 2.0E0*I_ESP_H2y3z_F2xy_a-2*I_ESP_F2yz_F2xy;
    abcd[iGrid*450+328] = 2.0E0*I_ESP_Hy4z_F2xy_a-3*I_ESP_Fy2z_F2xy;
    abcd[iGrid*450+329] = 2.0E0*I_ESP_H5z_F2xy_a-4*I_ESP_F3z_F2xy;
    abcd[iGrid*450+330] = 2.0E0*I_ESP_H4xz_F2xz_a;
    abcd[iGrid*450+331] = 2.0E0*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*450+332] = 2.0E0*I_ESP_H3x2z_F2xz_a-1*I_ESP_F3x_F2xz;
    abcd[iGrid*450+333] = 2.0E0*I_ESP_H2x2yz_F2xz_a;
    abcd[iGrid*450+334] = 2.0E0*I_ESP_H2xy2z_F2xz_a-1*I_ESP_F2xy_F2xz;
    abcd[iGrid*450+335] = 2.0E0*I_ESP_H2x3z_F2xz_a-2*I_ESP_F2xz_F2xz;
    abcd[iGrid*450+336] = 2.0E0*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*450+337] = 2.0E0*I_ESP_Hx2y2z_F2xz_a-1*I_ESP_Fx2y_F2xz;
    abcd[iGrid*450+338] = 2.0E0*I_ESP_Hxy3z_F2xz_a-2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*450+339] = 2.0E0*I_ESP_Hx4z_F2xz_a-3*I_ESP_Fx2z_F2xz;
    abcd[iGrid*450+340] = 2.0E0*I_ESP_H4yz_F2xz_a;
    abcd[iGrid*450+341] = 2.0E0*I_ESP_H3y2z_F2xz_a-1*I_ESP_F3y_F2xz;
    abcd[iGrid*450+342] = 2.0E0*I_ESP_H2y3z_F2xz_a-2*I_ESP_F2yz_F2xz;
    abcd[iGrid*450+343] = 2.0E0*I_ESP_Hy4z_F2xz_a-3*I_ESP_Fy2z_F2xz;
    abcd[iGrid*450+344] = 2.0E0*I_ESP_H5z_F2xz_a-4*I_ESP_F3z_F2xz;
    abcd[iGrid*450+345] = 2.0E0*I_ESP_H4xz_Fx2y_a;
    abcd[iGrid*450+346] = 2.0E0*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*450+347] = 2.0E0*I_ESP_H3x2z_Fx2y_a-1*I_ESP_F3x_Fx2y;
    abcd[iGrid*450+348] = 2.0E0*I_ESP_H2x2yz_Fx2y_a;
    abcd[iGrid*450+349] = 2.0E0*I_ESP_H2xy2z_Fx2y_a-1*I_ESP_F2xy_Fx2y;
    abcd[iGrid*450+350] = 2.0E0*I_ESP_H2x3z_Fx2y_a-2*I_ESP_F2xz_Fx2y;
    abcd[iGrid*450+351] = 2.0E0*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*450+352] = 2.0E0*I_ESP_Hx2y2z_Fx2y_a-1*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*450+353] = 2.0E0*I_ESP_Hxy3z_Fx2y_a-2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*450+354] = 2.0E0*I_ESP_Hx4z_Fx2y_a-3*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*450+355] = 2.0E0*I_ESP_H4yz_Fx2y_a;
    abcd[iGrid*450+356] = 2.0E0*I_ESP_H3y2z_Fx2y_a-1*I_ESP_F3y_Fx2y;
    abcd[iGrid*450+357] = 2.0E0*I_ESP_H2y3z_Fx2y_a-2*I_ESP_F2yz_Fx2y;
    abcd[iGrid*450+358] = 2.0E0*I_ESP_Hy4z_Fx2y_a-3*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*450+359] = 2.0E0*I_ESP_H5z_Fx2y_a-4*I_ESP_F3z_Fx2y;
    abcd[iGrid*450+360] = 2.0E0*I_ESP_H4xz_Fxyz_a;
    abcd[iGrid*450+361] = 2.0E0*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*450+362] = 2.0E0*I_ESP_H3x2z_Fxyz_a-1*I_ESP_F3x_Fxyz;
    abcd[iGrid*450+363] = 2.0E0*I_ESP_H2x2yz_Fxyz_a;
    abcd[iGrid*450+364] = 2.0E0*I_ESP_H2xy2z_Fxyz_a-1*I_ESP_F2xy_Fxyz;
    abcd[iGrid*450+365] = 2.0E0*I_ESP_H2x3z_Fxyz_a-2*I_ESP_F2xz_Fxyz;
    abcd[iGrid*450+366] = 2.0E0*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*450+367] = 2.0E0*I_ESP_Hx2y2z_Fxyz_a-1*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*450+368] = 2.0E0*I_ESP_Hxy3z_Fxyz_a-2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*450+369] = 2.0E0*I_ESP_Hx4z_Fxyz_a-3*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*450+370] = 2.0E0*I_ESP_H4yz_Fxyz_a;
    abcd[iGrid*450+371] = 2.0E0*I_ESP_H3y2z_Fxyz_a-1*I_ESP_F3y_Fxyz;
    abcd[iGrid*450+372] = 2.0E0*I_ESP_H2y3z_Fxyz_a-2*I_ESP_F2yz_Fxyz;
    abcd[iGrid*450+373] = 2.0E0*I_ESP_Hy4z_Fxyz_a-3*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*450+374] = 2.0E0*I_ESP_H5z_Fxyz_a-4*I_ESP_F3z_Fxyz;
    abcd[iGrid*450+375] = 2.0E0*I_ESP_H4xz_Fx2z_a;
    abcd[iGrid*450+376] = 2.0E0*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*450+377] = 2.0E0*I_ESP_H3x2z_Fx2z_a-1*I_ESP_F3x_Fx2z;
    abcd[iGrid*450+378] = 2.0E0*I_ESP_H2x2yz_Fx2z_a;
    abcd[iGrid*450+379] = 2.0E0*I_ESP_H2xy2z_Fx2z_a-1*I_ESP_F2xy_Fx2z;
    abcd[iGrid*450+380] = 2.0E0*I_ESP_H2x3z_Fx2z_a-2*I_ESP_F2xz_Fx2z;
    abcd[iGrid*450+381] = 2.0E0*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*450+382] = 2.0E0*I_ESP_Hx2y2z_Fx2z_a-1*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*450+383] = 2.0E0*I_ESP_Hxy3z_Fx2z_a-2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*450+384] = 2.0E0*I_ESP_Hx4z_Fx2z_a-3*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*450+385] = 2.0E0*I_ESP_H4yz_Fx2z_a;
    abcd[iGrid*450+386] = 2.0E0*I_ESP_H3y2z_Fx2z_a-1*I_ESP_F3y_Fx2z;
    abcd[iGrid*450+387] = 2.0E0*I_ESP_H2y3z_Fx2z_a-2*I_ESP_F2yz_Fx2z;
    abcd[iGrid*450+388] = 2.0E0*I_ESP_Hy4z_Fx2z_a-3*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*450+389] = 2.0E0*I_ESP_H5z_Fx2z_a-4*I_ESP_F3z_Fx2z;
    abcd[iGrid*450+390] = 2.0E0*I_ESP_H4xz_F3y_a;
    abcd[iGrid*450+391] = 2.0E0*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*450+392] = 2.0E0*I_ESP_H3x2z_F3y_a-1*I_ESP_F3x_F3y;
    abcd[iGrid*450+393] = 2.0E0*I_ESP_H2x2yz_F3y_a;
    abcd[iGrid*450+394] = 2.0E0*I_ESP_H2xy2z_F3y_a-1*I_ESP_F2xy_F3y;
    abcd[iGrid*450+395] = 2.0E0*I_ESP_H2x3z_F3y_a-2*I_ESP_F2xz_F3y;
    abcd[iGrid*450+396] = 2.0E0*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*450+397] = 2.0E0*I_ESP_Hx2y2z_F3y_a-1*I_ESP_Fx2y_F3y;
    abcd[iGrid*450+398] = 2.0E0*I_ESP_Hxy3z_F3y_a-2*I_ESP_Fxyz_F3y;
    abcd[iGrid*450+399] = 2.0E0*I_ESP_Hx4z_F3y_a-3*I_ESP_Fx2z_F3y;
    abcd[iGrid*450+400] = 2.0E0*I_ESP_H4yz_F3y_a;
    abcd[iGrid*450+401] = 2.0E0*I_ESP_H3y2z_F3y_a-1*I_ESP_F3y_F3y;
    abcd[iGrid*450+402] = 2.0E0*I_ESP_H2y3z_F3y_a-2*I_ESP_F2yz_F3y;
    abcd[iGrid*450+403] = 2.0E0*I_ESP_Hy4z_F3y_a-3*I_ESP_Fy2z_F3y;
    abcd[iGrid*450+404] = 2.0E0*I_ESP_H5z_F3y_a-4*I_ESP_F3z_F3y;
    abcd[iGrid*450+405] = 2.0E0*I_ESP_H4xz_F2yz_a;
    abcd[iGrid*450+406] = 2.0E0*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*450+407] = 2.0E0*I_ESP_H3x2z_F2yz_a-1*I_ESP_F3x_F2yz;
    abcd[iGrid*450+408] = 2.0E0*I_ESP_H2x2yz_F2yz_a;
    abcd[iGrid*450+409] = 2.0E0*I_ESP_H2xy2z_F2yz_a-1*I_ESP_F2xy_F2yz;
    abcd[iGrid*450+410] = 2.0E0*I_ESP_H2x3z_F2yz_a-2*I_ESP_F2xz_F2yz;
    abcd[iGrid*450+411] = 2.0E0*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*450+412] = 2.0E0*I_ESP_Hx2y2z_F2yz_a-1*I_ESP_Fx2y_F2yz;
    abcd[iGrid*450+413] = 2.0E0*I_ESP_Hxy3z_F2yz_a-2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*450+414] = 2.0E0*I_ESP_Hx4z_F2yz_a-3*I_ESP_Fx2z_F2yz;
    abcd[iGrid*450+415] = 2.0E0*I_ESP_H4yz_F2yz_a;
    abcd[iGrid*450+416] = 2.0E0*I_ESP_H3y2z_F2yz_a-1*I_ESP_F3y_F2yz;
    abcd[iGrid*450+417] = 2.0E0*I_ESP_H2y3z_F2yz_a-2*I_ESP_F2yz_F2yz;
    abcd[iGrid*450+418] = 2.0E0*I_ESP_Hy4z_F2yz_a-3*I_ESP_Fy2z_F2yz;
    abcd[iGrid*450+419] = 2.0E0*I_ESP_H5z_F2yz_a-4*I_ESP_F3z_F2yz;
    abcd[iGrid*450+420] = 2.0E0*I_ESP_H4xz_Fy2z_a;
    abcd[iGrid*450+421] = 2.0E0*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*450+422] = 2.0E0*I_ESP_H3x2z_Fy2z_a-1*I_ESP_F3x_Fy2z;
    abcd[iGrid*450+423] = 2.0E0*I_ESP_H2x2yz_Fy2z_a;
    abcd[iGrid*450+424] = 2.0E0*I_ESP_H2xy2z_Fy2z_a-1*I_ESP_F2xy_Fy2z;
    abcd[iGrid*450+425] = 2.0E0*I_ESP_H2x3z_Fy2z_a-2*I_ESP_F2xz_Fy2z;
    abcd[iGrid*450+426] = 2.0E0*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*450+427] = 2.0E0*I_ESP_Hx2y2z_Fy2z_a-1*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*450+428] = 2.0E0*I_ESP_Hxy3z_Fy2z_a-2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*450+429] = 2.0E0*I_ESP_Hx4z_Fy2z_a-3*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*450+430] = 2.0E0*I_ESP_H4yz_Fy2z_a;
    abcd[iGrid*450+431] = 2.0E0*I_ESP_H3y2z_Fy2z_a-1*I_ESP_F3y_Fy2z;
    abcd[iGrid*450+432] = 2.0E0*I_ESP_H2y3z_Fy2z_a-2*I_ESP_F2yz_Fy2z;
    abcd[iGrid*450+433] = 2.0E0*I_ESP_Hy4z_Fy2z_a-3*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*450+434] = 2.0E0*I_ESP_H5z_Fy2z_a-4*I_ESP_F3z_Fy2z;
    abcd[iGrid*450+435] = 2.0E0*I_ESP_H4xz_F3z_a;
    abcd[iGrid*450+436] = 2.0E0*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*450+437] = 2.0E0*I_ESP_H3x2z_F3z_a-1*I_ESP_F3x_F3z;
    abcd[iGrid*450+438] = 2.0E0*I_ESP_H2x2yz_F3z_a;
    abcd[iGrid*450+439] = 2.0E0*I_ESP_H2xy2z_F3z_a-1*I_ESP_F2xy_F3z;
    abcd[iGrid*450+440] = 2.0E0*I_ESP_H2x3z_F3z_a-2*I_ESP_F2xz_F3z;
    abcd[iGrid*450+441] = 2.0E0*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*450+442] = 2.0E0*I_ESP_Hx2y2z_F3z_a-1*I_ESP_Fx2y_F3z;
    abcd[iGrid*450+443] = 2.0E0*I_ESP_Hxy3z_F3z_a-2*I_ESP_Fxyz_F3z;
    abcd[iGrid*450+444] = 2.0E0*I_ESP_Hx4z_F3z_a-3*I_ESP_Fx2z_F3z;
    abcd[iGrid*450+445] = 2.0E0*I_ESP_H4yz_F3z_a;
    abcd[iGrid*450+446] = 2.0E0*I_ESP_H3y2z_F3z_a-1*I_ESP_F3y_F3z;
    abcd[iGrid*450+447] = 2.0E0*I_ESP_H2y3z_F3z_a-2*I_ESP_F2yz_F3z;
    abcd[iGrid*450+448] = 2.0E0*I_ESP_Hy4z_F3z_a-3*I_ESP_Fy2z_F3z;
    abcd[iGrid*450+449] = 2.0E0*I_ESP_H5z_F3z_a-4*I_ESP_F3z_F3z;
  }
}
