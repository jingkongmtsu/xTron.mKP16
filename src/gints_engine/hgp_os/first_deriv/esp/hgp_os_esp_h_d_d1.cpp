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
// BRA1 as redundant position, total RHS integrals evaluated as: 3057
// BRA2 as redundant position, total RHS integrals evaluated as: 2838
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

void hgp_os_esp_h_d_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
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
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 7 integrals are omitted 
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
    Double I_ESP_H5y_Px = I_ESP_Ix5y_S+ABX*I_ESP_H5y_S;
    Double I_ESP_H4yz_Px = I_ESP_Ix4yz_S+ABX*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Px = I_ESP_Ix3y2z_S+ABX*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Px = I_ESP_Ix2y3z_S+ABX*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Px = I_ESP_Ixy4z_S+ABX*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Px = I_ESP_Ix5z_S+ABX*I_ESP_H5z_S;
    Double I_ESP_H4xy_Py = I_ESP_I4x2y_S+ABY*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Py = I_ESP_I4xyz_S+ABY*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Py = I_ESP_I3x3y_S+ABY*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Py = I_ESP_I3x2yz_S+ABY*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Py = I_ESP_I3xy2z_S+ABY*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Py = I_ESP_I2x4y_S+ABY*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Py = I_ESP_I2x3yz_S+ABY*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Py = I_ESP_I2x2y2z_S+ABY*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Py = I_ESP_I2xy3z_S+ABY*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Py = I_ESP_Ix5y_S+ABY*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Py = I_ESP_Ix4yz_S+ABY*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Py = I_ESP_Ix3y2z_S+ABY*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Py = I_ESP_Ix2y3z_S+ABY*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Py = I_ESP_Ixy4z_S+ABY*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Py = I_ESP_I6y_S+ABY*I_ESP_H5y_S;
    Double I_ESP_H4yz_Py = I_ESP_I5yz_S+ABY*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Py = I_ESP_I4y2z_S+ABY*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Py = I_ESP_I3y3z_S+ABY*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Py = I_ESP_I2y4z_S+ABY*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Py = I_ESP_Iy5z_S+ABY*I_ESP_H5z_S;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_G4x_Dxy = I_ESP_H4xy_Px+ABY*I_ESP_G4x_Px;
    Double I_ESP_G3xy_Dxy = I_ESP_H3x2y_Px+ABY*I_ESP_G3xy_Px;
    Double I_ESP_G3xz_Dxy = I_ESP_H3xyz_Px+ABY*I_ESP_G3xz_Px;
    Double I_ESP_G2x2y_Dxy = I_ESP_H2x3y_Px+ABY*I_ESP_G2x2y_Px;
    Double I_ESP_G2xyz_Dxy = I_ESP_H2x2yz_Px+ABY*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_Dxy = I_ESP_H2xy2z_Px+ABY*I_ESP_G2x2z_Px;
    Double I_ESP_Gx3y_Dxy = I_ESP_Hx4y_Px+ABY*I_ESP_Gx3y_Px;
    Double I_ESP_Gx2yz_Dxy = I_ESP_Hx3yz_Px+ABY*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_Dxy = I_ESP_Hx2y2z_Px+ABY*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_Dxy = I_ESP_Hxy3z_Px+ABY*I_ESP_Gx3z_Px;
    Double I_ESP_G4y_Dxy = I_ESP_H5y_Px+ABY*I_ESP_G4y_Px;
    Double I_ESP_G3yz_Dxy = I_ESP_H4yz_Px+ABY*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_Dxy = I_ESP_H3y2z_Px+ABY*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_Dxy = I_ESP_H2y3z_Px+ABY*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_Dxy = I_ESP_Hy4z_Px+ABY*I_ESP_G4z_Px;
    Double I_ESP_G4x_Dxz = I_ESP_H4xz_Px+ABZ*I_ESP_G4x_Px;
    Double I_ESP_G3xy_Dxz = I_ESP_H3xyz_Px+ABZ*I_ESP_G3xy_Px;
    Double I_ESP_G3xz_Dxz = I_ESP_H3x2z_Px+ABZ*I_ESP_G3xz_Px;
    Double I_ESP_G2x2y_Dxz = I_ESP_H2x2yz_Px+ABZ*I_ESP_G2x2y_Px;
    Double I_ESP_G2xyz_Dxz = I_ESP_H2xy2z_Px+ABZ*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_Dxz = I_ESP_H2x3z_Px+ABZ*I_ESP_G2x2z_Px;
    Double I_ESP_Gx3y_Dxz = I_ESP_Hx3yz_Px+ABZ*I_ESP_Gx3y_Px;
    Double I_ESP_Gx2yz_Dxz = I_ESP_Hx2y2z_Px+ABZ*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_Dxz = I_ESP_Hxy3z_Px+ABZ*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_Dxz = I_ESP_Hx4z_Px+ABZ*I_ESP_Gx3z_Px;
    Double I_ESP_G4y_Dxz = I_ESP_H4yz_Px+ABZ*I_ESP_G4y_Px;
    Double I_ESP_G3yz_Dxz = I_ESP_H3y2z_Px+ABZ*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_Dxz = I_ESP_H2y3z_Px+ABZ*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_Dxz = I_ESP_Hy4z_Px+ABZ*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_Dxz = I_ESP_H5z_Px+ABZ*I_ESP_G4z_Px;
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
    Double I_ESP_G4x_Dyz = I_ESP_H4xz_Py+ABZ*I_ESP_G4x_Py;
    Double I_ESP_G3xy_Dyz = I_ESP_H3xyz_Py+ABZ*I_ESP_G3xy_Py;
    Double I_ESP_G3xz_Dyz = I_ESP_H3x2z_Py+ABZ*I_ESP_G3xz_Py;
    Double I_ESP_G2x2y_Dyz = I_ESP_H2x2yz_Py+ABZ*I_ESP_G2x2y_Py;
    Double I_ESP_G2xyz_Dyz = I_ESP_H2xy2z_Py+ABZ*I_ESP_G2xyz_Py;
    Double I_ESP_G2x2z_Dyz = I_ESP_H2x3z_Py+ABZ*I_ESP_G2x2z_Py;
    Double I_ESP_Gx3y_Dyz = I_ESP_Hx3yz_Py+ABZ*I_ESP_Gx3y_Py;
    Double I_ESP_Gx2yz_Dyz = I_ESP_Hx2y2z_Py+ABZ*I_ESP_Gx2yz_Py;
    Double I_ESP_Gxy2z_Dyz = I_ESP_Hxy3z_Py+ABZ*I_ESP_Gxy2z_Py;
    Double I_ESP_Gx3z_Dyz = I_ESP_Hx4z_Py+ABZ*I_ESP_Gx3z_Py;
    Double I_ESP_G4y_Dyz = I_ESP_H4yz_Py+ABZ*I_ESP_G4y_Py;
    Double I_ESP_G3yz_Dyz = I_ESP_H3y2z_Py+ABZ*I_ESP_G3yz_Py;
    Double I_ESP_G2y2z_Dyz = I_ESP_H2y3z_Py+ABZ*I_ESP_G2y2z_Py;
    Double I_ESP_Gy3z_Dyz = I_ESP_Hy4z_Py+ABZ*I_ESP_Gy3z_Py;
    Double I_ESP_G4z_Dyz = I_ESP_H5z_Py+ABZ*I_ESP_G4z_Py;
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
     * shell quartet name: SQ_ESP_K_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 9 integrals are omitted 
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
    Double I_ESP_K7y_Px_a = I_ESP_Lx7y_S_a+ABX*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Px_a = I_ESP_Lx6yz_S_a+ABX*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Px_a = I_ESP_Lx5y2z_S_a+ABX*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Px_a = I_ESP_Lx4y3z_S_a+ABX*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Px_a = I_ESP_Lx3y4z_S_a+ABX*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Px_a = I_ESP_Lx2y5z_S_a+ABX*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Px_a = I_ESP_Lxy6z_S_a+ABX*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Px_a = I_ESP_Lx7z_S_a+ABX*I_ESP_K7z_S_a;
    Double I_ESP_K6xy_Py_a = I_ESP_L6x2y_S_a+ABY*I_ESP_K6xy_S_a;
    Double I_ESP_K6xz_Py_a = I_ESP_L6xyz_S_a+ABY*I_ESP_K6xz_S_a;
    Double I_ESP_K5x2y_Py_a = I_ESP_L5x3y_S_a+ABY*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Py_a = I_ESP_L5x2yz_S_a+ABY*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Py_a = I_ESP_L5xy2z_S_a+ABY*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x3y_Py_a = I_ESP_L4x4y_S_a+ABY*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Py_a = I_ESP_L4x3yz_S_a+ABY*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Py_a = I_ESP_L4x2y2z_S_a+ABY*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Py_a = I_ESP_L4xy3z_S_a+ABY*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x4y_Py_a = I_ESP_L3x5y_S_a+ABY*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Py_a = I_ESP_L3x4yz_S_a+ABY*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Py_a = I_ESP_L3x3y2z_S_a+ABY*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Py_a = I_ESP_L3x2y3z_S_a+ABY*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Py_a = I_ESP_L3xy4z_S_a+ABY*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x5y_Py_a = I_ESP_L2x6y_S_a+ABY*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Py_a = I_ESP_L2x5yz_S_a+ABY*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Py_a = I_ESP_L2x4y2z_S_a+ABY*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Py_a = I_ESP_L2x3y3z_S_a+ABY*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Py_a = I_ESP_L2x2y4z_S_a+ABY*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Py_a = I_ESP_L2xy5z_S_a+ABY*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx6y_Py_a = I_ESP_Lx7y_S_a+ABY*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Py_a = I_ESP_Lx6yz_S_a+ABY*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Py_a = I_ESP_Lx5y2z_S_a+ABY*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Py_a = I_ESP_Lx4y3z_S_a+ABY*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Py_a = I_ESP_Lx3y4z_S_a+ABY*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Py_a = I_ESP_Lx2y5z_S_a+ABY*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Py_a = I_ESP_Lxy6z_S_a+ABY*I_ESP_Kx6z_S_a;
    Double I_ESP_K7y_Py_a = I_ESP_L8y_S_a+ABY*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Py_a = I_ESP_L7yz_S_a+ABY*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Py_a = I_ESP_L6y2z_S_a+ABY*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Py_a = I_ESP_L5y3z_S_a+ABY*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Py_a = I_ESP_L4y4z_S_a+ABY*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Py_a = I_ESP_L3y5z_S_a+ABY*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Py_a = I_ESP_L2y6z_S_a+ABY*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Py_a = I_ESP_Ly7z_S_a+ABY*I_ESP_K7z_S_a;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_I6x_Dxy_a = I_ESP_K6xy_Px_a+ABY*I_ESP_I6x_Px_a;
    Double I_ESP_I5xy_Dxy_a = I_ESP_K5x2y_Px_a+ABY*I_ESP_I5xy_Px_a;
    Double I_ESP_I5xz_Dxy_a = I_ESP_K5xyz_Px_a+ABY*I_ESP_I5xz_Px_a;
    Double I_ESP_I4x2y_Dxy_a = I_ESP_K4x3y_Px_a+ABY*I_ESP_I4x2y_Px_a;
    Double I_ESP_I4xyz_Dxy_a = I_ESP_K4x2yz_Px_a+ABY*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_Dxy_a = I_ESP_K4xy2z_Px_a+ABY*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x3y_Dxy_a = I_ESP_K3x4y_Px_a+ABY*I_ESP_I3x3y_Px_a;
    Double I_ESP_I3x2yz_Dxy_a = I_ESP_K3x3yz_Px_a+ABY*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_Dxy_a = I_ESP_K3x2y2z_Px_a+ABY*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_Dxy_a = I_ESP_K3xy3z_Px_a+ABY*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x4y_Dxy_a = I_ESP_K2x5y_Px_a+ABY*I_ESP_I2x4y_Px_a;
    Double I_ESP_I2x3yz_Dxy_a = I_ESP_K2x4yz_Px_a+ABY*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_Dxy_a = I_ESP_K2x3y2z_Px_a+ABY*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_Dxy_a = I_ESP_K2x2y3z_Px_a+ABY*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_Dxy_a = I_ESP_K2xy4z_Px_a+ABY*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix5y_Dxy_a = I_ESP_Kx6y_Px_a+ABY*I_ESP_Ix5y_Px_a;
    Double I_ESP_Ix4yz_Dxy_a = I_ESP_Kx5yz_Px_a+ABY*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_Dxy_a = I_ESP_Kx4y2z_Px_a+ABY*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_Dxy_a = I_ESP_Kx3y3z_Px_a+ABY*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_Dxy_a = I_ESP_Kx2y4z_Px_a+ABY*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_Dxy_a = I_ESP_Kxy5z_Px_a+ABY*I_ESP_Ix5z_Px_a;
    Double I_ESP_I6y_Dxy_a = I_ESP_K7y_Px_a+ABY*I_ESP_I6y_Px_a;
    Double I_ESP_I5yz_Dxy_a = I_ESP_K6yz_Px_a+ABY*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_Dxy_a = I_ESP_K5y2z_Px_a+ABY*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_Dxy_a = I_ESP_K4y3z_Px_a+ABY*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_Dxy_a = I_ESP_K3y4z_Px_a+ABY*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_Dxy_a = I_ESP_K2y5z_Px_a+ABY*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_Dxy_a = I_ESP_Ky6z_Px_a+ABY*I_ESP_I6z_Px_a;
    Double I_ESP_I6x_Dxz_a = I_ESP_K6xz_Px_a+ABZ*I_ESP_I6x_Px_a;
    Double I_ESP_I5xy_Dxz_a = I_ESP_K5xyz_Px_a+ABZ*I_ESP_I5xy_Px_a;
    Double I_ESP_I5xz_Dxz_a = I_ESP_K5x2z_Px_a+ABZ*I_ESP_I5xz_Px_a;
    Double I_ESP_I4x2y_Dxz_a = I_ESP_K4x2yz_Px_a+ABZ*I_ESP_I4x2y_Px_a;
    Double I_ESP_I4xyz_Dxz_a = I_ESP_K4xy2z_Px_a+ABZ*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_Dxz_a = I_ESP_K4x3z_Px_a+ABZ*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x3y_Dxz_a = I_ESP_K3x3yz_Px_a+ABZ*I_ESP_I3x3y_Px_a;
    Double I_ESP_I3x2yz_Dxz_a = I_ESP_K3x2y2z_Px_a+ABZ*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_Dxz_a = I_ESP_K3xy3z_Px_a+ABZ*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_Dxz_a = I_ESP_K3x4z_Px_a+ABZ*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x4y_Dxz_a = I_ESP_K2x4yz_Px_a+ABZ*I_ESP_I2x4y_Px_a;
    Double I_ESP_I2x3yz_Dxz_a = I_ESP_K2x3y2z_Px_a+ABZ*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_Dxz_a = I_ESP_K2x2y3z_Px_a+ABZ*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_Dxz_a = I_ESP_K2xy4z_Px_a+ABZ*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_Dxz_a = I_ESP_K2x5z_Px_a+ABZ*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix5y_Dxz_a = I_ESP_Kx5yz_Px_a+ABZ*I_ESP_Ix5y_Px_a;
    Double I_ESP_Ix4yz_Dxz_a = I_ESP_Kx4y2z_Px_a+ABZ*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_Dxz_a = I_ESP_Kx3y3z_Px_a+ABZ*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_Dxz_a = I_ESP_Kx2y4z_Px_a+ABZ*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_Dxz_a = I_ESP_Kxy5z_Px_a+ABZ*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_Dxz_a = I_ESP_Kx6z_Px_a+ABZ*I_ESP_Ix5z_Px_a;
    Double I_ESP_I6y_Dxz_a = I_ESP_K6yz_Px_a+ABZ*I_ESP_I6y_Px_a;
    Double I_ESP_I5yz_Dxz_a = I_ESP_K5y2z_Px_a+ABZ*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_Dxz_a = I_ESP_K4y3z_Px_a+ABZ*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_Dxz_a = I_ESP_K3y4z_Px_a+ABZ*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_Dxz_a = I_ESP_K2y5z_Px_a+ABZ*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_Dxz_a = I_ESP_Ky6z_Px_a+ABZ*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_Dxz_a = I_ESP_K7z_Px_a+ABZ*I_ESP_I6z_Px_a;
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
    Double I_ESP_I6x_Dyz_a = I_ESP_K6xz_Py_a+ABZ*I_ESP_I6x_Py_a;
    Double I_ESP_I5xy_Dyz_a = I_ESP_K5xyz_Py_a+ABZ*I_ESP_I5xy_Py_a;
    Double I_ESP_I5xz_Dyz_a = I_ESP_K5x2z_Py_a+ABZ*I_ESP_I5xz_Py_a;
    Double I_ESP_I4x2y_Dyz_a = I_ESP_K4x2yz_Py_a+ABZ*I_ESP_I4x2y_Py_a;
    Double I_ESP_I4xyz_Dyz_a = I_ESP_K4xy2z_Py_a+ABZ*I_ESP_I4xyz_Py_a;
    Double I_ESP_I4x2z_Dyz_a = I_ESP_K4x3z_Py_a+ABZ*I_ESP_I4x2z_Py_a;
    Double I_ESP_I3x3y_Dyz_a = I_ESP_K3x3yz_Py_a+ABZ*I_ESP_I3x3y_Py_a;
    Double I_ESP_I3x2yz_Dyz_a = I_ESP_K3x2y2z_Py_a+ABZ*I_ESP_I3x2yz_Py_a;
    Double I_ESP_I3xy2z_Dyz_a = I_ESP_K3xy3z_Py_a+ABZ*I_ESP_I3xy2z_Py_a;
    Double I_ESP_I3x3z_Dyz_a = I_ESP_K3x4z_Py_a+ABZ*I_ESP_I3x3z_Py_a;
    Double I_ESP_I2x4y_Dyz_a = I_ESP_K2x4yz_Py_a+ABZ*I_ESP_I2x4y_Py_a;
    Double I_ESP_I2x3yz_Dyz_a = I_ESP_K2x3y2z_Py_a+ABZ*I_ESP_I2x3yz_Py_a;
    Double I_ESP_I2x2y2z_Dyz_a = I_ESP_K2x2y3z_Py_a+ABZ*I_ESP_I2x2y2z_Py_a;
    Double I_ESP_I2xy3z_Dyz_a = I_ESP_K2xy4z_Py_a+ABZ*I_ESP_I2xy3z_Py_a;
    Double I_ESP_I2x4z_Dyz_a = I_ESP_K2x5z_Py_a+ABZ*I_ESP_I2x4z_Py_a;
    Double I_ESP_Ix5y_Dyz_a = I_ESP_Kx5yz_Py_a+ABZ*I_ESP_Ix5y_Py_a;
    Double I_ESP_Ix4yz_Dyz_a = I_ESP_Kx4y2z_Py_a+ABZ*I_ESP_Ix4yz_Py_a;
    Double I_ESP_Ix3y2z_Dyz_a = I_ESP_Kx3y3z_Py_a+ABZ*I_ESP_Ix3y2z_Py_a;
    Double I_ESP_Ix2y3z_Dyz_a = I_ESP_Kx2y4z_Py_a+ABZ*I_ESP_Ix2y3z_Py_a;
    Double I_ESP_Ixy4z_Dyz_a = I_ESP_Kxy5z_Py_a+ABZ*I_ESP_Ixy4z_Py_a;
    Double I_ESP_Ix5z_Dyz_a = I_ESP_Kx6z_Py_a+ABZ*I_ESP_Ix5z_Py_a;
    Double I_ESP_I6y_Dyz_a = I_ESP_K6yz_Py_a+ABZ*I_ESP_I6y_Py_a;
    Double I_ESP_I5yz_Dyz_a = I_ESP_K5y2z_Py_a+ABZ*I_ESP_I5yz_Py_a;
    Double I_ESP_I4y2z_Dyz_a = I_ESP_K4y3z_Py_a+ABZ*I_ESP_I4y2z_Py_a;
    Double I_ESP_I3y3z_Dyz_a = I_ESP_K3y4z_Py_a+ABZ*I_ESP_I3y3z_Py_a;
    Double I_ESP_I2y4z_Dyz_a = I_ESP_K2y5z_Py_a+ABZ*I_ESP_I2y4z_Py_a;
    Double I_ESP_Iy5z_Dyz_a = I_ESP_Ky6z_Py_a+ABZ*I_ESP_Iy5z_Py_a;
    Double I_ESP_I6z_Dyz_a = I_ESP_K7z_Py_a+ABZ*I_ESP_I6z_Py_a;
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
     * shell quartet name: SQ_ESP_H_D_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_a
     * RHS shell quartet name: SQ_ESP_G_D
     ************************************************************/
    abcd[iGrid*378+0] = 2.0E0*I_ESP_I6x_D2x_a-5*I_ESP_G4x_D2x;
    abcd[iGrid*378+1] = 2.0E0*I_ESP_I5xy_D2x_a-4*I_ESP_G3xy_D2x;
    abcd[iGrid*378+2] = 2.0E0*I_ESP_I5xz_D2x_a-4*I_ESP_G3xz_D2x;
    abcd[iGrid*378+3] = 2.0E0*I_ESP_I4x2y_D2x_a-3*I_ESP_G2x2y_D2x;
    abcd[iGrid*378+4] = 2.0E0*I_ESP_I4xyz_D2x_a-3*I_ESP_G2xyz_D2x;
    abcd[iGrid*378+5] = 2.0E0*I_ESP_I4x2z_D2x_a-3*I_ESP_G2x2z_D2x;
    abcd[iGrid*378+6] = 2.0E0*I_ESP_I3x3y_D2x_a-2*I_ESP_Gx3y_D2x;
    abcd[iGrid*378+7] = 2.0E0*I_ESP_I3x2yz_D2x_a-2*I_ESP_Gx2yz_D2x;
    abcd[iGrid*378+8] = 2.0E0*I_ESP_I3xy2z_D2x_a-2*I_ESP_Gxy2z_D2x;
    abcd[iGrid*378+9] = 2.0E0*I_ESP_I3x3z_D2x_a-2*I_ESP_Gx3z_D2x;
    abcd[iGrid*378+10] = 2.0E0*I_ESP_I2x4y_D2x_a-1*I_ESP_G4y_D2x;
    abcd[iGrid*378+11] = 2.0E0*I_ESP_I2x3yz_D2x_a-1*I_ESP_G3yz_D2x;
    abcd[iGrid*378+12] = 2.0E0*I_ESP_I2x2y2z_D2x_a-1*I_ESP_G2y2z_D2x;
    abcd[iGrid*378+13] = 2.0E0*I_ESP_I2xy3z_D2x_a-1*I_ESP_Gy3z_D2x;
    abcd[iGrid*378+14] = 2.0E0*I_ESP_I2x4z_D2x_a-1*I_ESP_G4z_D2x;
    abcd[iGrid*378+15] = 2.0E0*I_ESP_Ix5y_D2x_a;
    abcd[iGrid*378+16] = 2.0E0*I_ESP_Ix4yz_D2x_a;
    abcd[iGrid*378+17] = 2.0E0*I_ESP_Ix3y2z_D2x_a;
    abcd[iGrid*378+18] = 2.0E0*I_ESP_Ix2y3z_D2x_a;
    abcd[iGrid*378+19] = 2.0E0*I_ESP_Ixy4z_D2x_a;
    abcd[iGrid*378+20] = 2.0E0*I_ESP_Ix5z_D2x_a;
    abcd[iGrid*378+21] = 2.0E0*I_ESP_I6x_Dxy_a-5*I_ESP_G4x_Dxy;
    abcd[iGrid*378+22] = 2.0E0*I_ESP_I5xy_Dxy_a-4*I_ESP_G3xy_Dxy;
    abcd[iGrid*378+23] = 2.0E0*I_ESP_I5xz_Dxy_a-4*I_ESP_G3xz_Dxy;
    abcd[iGrid*378+24] = 2.0E0*I_ESP_I4x2y_Dxy_a-3*I_ESP_G2x2y_Dxy;
    abcd[iGrid*378+25] = 2.0E0*I_ESP_I4xyz_Dxy_a-3*I_ESP_G2xyz_Dxy;
    abcd[iGrid*378+26] = 2.0E0*I_ESP_I4x2z_Dxy_a-3*I_ESP_G2x2z_Dxy;
    abcd[iGrid*378+27] = 2.0E0*I_ESP_I3x3y_Dxy_a-2*I_ESP_Gx3y_Dxy;
    abcd[iGrid*378+28] = 2.0E0*I_ESP_I3x2yz_Dxy_a-2*I_ESP_Gx2yz_Dxy;
    abcd[iGrid*378+29] = 2.0E0*I_ESP_I3xy2z_Dxy_a-2*I_ESP_Gxy2z_Dxy;
    abcd[iGrid*378+30] = 2.0E0*I_ESP_I3x3z_Dxy_a-2*I_ESP_Gx3z_Dxy;
    abcd[iGrid*378+31] = 2.0E0*I_ESP_I2x4y_Dxy_a-1*I_ESP_G4y_Dxy;
    abcd[iGrid*378+32] = 2.0E0*I_ESP_I2x3yz_Dxy_a-1*I_ESP_G3yz_Dxy;
    abcd[iGrid*378+33] = 2.0E0*I_ESP_I2x2y2z_Dxy_a-1*I_ESP_G2y2z_Dxy;
    abcd[iGrid*378+34] = 2.0E0*I_ESP_I2xy3z_Dxy_a-1*I_ESP_Gy3z_Dxy;
    abcd[iGrid*378+35] = 2.0E0*I_ESP_I2x4z_Dxy_a-1*I_ESP_G4z_Dxy;
    abcd[iGrid*378+36] = 2.0E0*I_ESP_Ix5y_Dxy_a;
    abcd[iGrid*378+37] = 2.0E0*I_ESP_Ix4yz_Dxy_a;
    abcd[iGrid*378+38] = 2.0E0*I_ESP_Ix3y2z_Dxy_a;
    abcd[iGrid*378+39] = 2.0E0*I_ESP_Ix2y3z_Dxy_a;
    abcd[iGrid*378+40] = 2.0E0*I_ESP_Ixy4z_Dxy_a;
    abcd[iGrid*378+41] = 2.0E0*I_ESP_Ix5z_Dxy_a;
    abcd[iGrid*378+42] = 2.0E0*I_ESP_I6x_Dxz_a-5*I_ESP_G4x_Dxz;
    abcd[iGrid*378+43] = 2.0E0*I_ESP_I5xy_Dxz_a-4*I_ESP_G3xy_Dxz;
    abcd[iGrid*378+44] = 2.0E0*I_ESP_I5xz_Dxz_a-4*I_ESP_G3xz_Dxz;
    abcd[iGrid*378+45] = 2.0E0*I_ESP_I4x2y_Dxz_a-3*I_ESP_G2x2y_Dxz;
    abcd[iGrid*378+46] = 2.0E0*I_ESP_I4xyz_Dxz_a-3*I_ESP_G2xyz_Dxz;
    abcd[iGrid*378+47] = 2.0E0*I_ESP_I4x2z_Dxz_a-3*I_ESP_G2x2z_Dxz;
    abcd[iGrid*378+48] = 2.0E0*I_ESP_I3x3y_Dxz_a-2*I_ESP_Gx3y_Dxz;
    abcd[iGrid*378+49] = 2.0E0*I_ESP_I3x2yz_Dxz_a-2*I_ESP_Gx2yz_Dxz;
    abcd[iGrid*378+50] = 2.0E0*I_ESP_I3xy2z_Dxz_a-2*I_ESP_Gxy2z_Dxz;
    abcd[iGrid*378+51] = 2.0E0*I_ESP_I3x3z_Dxz_a-2*I_ESP_Gx3z_Dxz;
    abcd[iGrid*378+52] = 2.0E0*I_ESP_I2x4y_Dxz_a-1*I_ESP_G4y_Dxz;
    abcd[iGrid*378+53] = 2.0E0*I_ESP_I2x3yz_Dxz_a-1*I_ESP_G3yz_Dxz;
    abcd[iGrid*378+54] = 2.0E0*I_ESP_I2x2y2z_Dxz_a-1*I_ESP_G2y2z_Dxz;
    abcd[iGrid*378+55] = 2.0E0*I_ESP_I2xy3z_Dxz_a-1*I_ESP_Gy3z_Dxz;
    abcd[iGrid*378+56] = 2.0E0*I_ESP_I2x4z_Dxz_a-1*I_ESP_G4z_Dxz;
    abcd[iGrid*378+57] = 2.0E0*I_ESP_Ix5y_Dxz_a;
    abcd[iGrid*378+58] = 2.0E0*I_ESP_Ix4yz_Dxz_a;
    abcd[iGrid*378+59] = 2.0E0*I_ESP_Ix3y2z_Dxz_a;
    abcd[iGrid*378+60] = 2.0E0*I_ESP_Ix2y3z_Dxz_a;
    abcd[iGrid*378+61] = 2.0E0*I_ESP_Ixy4z_Dxz_a;
    abcd[iGrid*378+62] = 2.0E0*I_ESP_Ix5z_Dxz_a;
    abcd[iGrid*378+63] = 2.0E0*I_ESP_I6x_D2y_a-5*I_ESP_G4x_D2y;
    abcd[iGrid*378+64] = 2.0E0*I_ESP_I5xy_D2y_a-4*I_ESP_G3xy_D2y;
    abcd[iGrid*378+65] = 2.0E0*I_ESP_I5xz_D2y_a-4*I_ESP_G3xz_D2y;
    abcd[iGrid*378+66] = 2.0E0*I_ESP_I4x2y_D2y_a-3*I_ESP_G2x2y_D2y;
    abcd[iGrid*378+67] = 2.0E0*I_ESP_I4xyz_D2y_a-3*I_ESP_G2xyz_D2y;
    abcd[iGrid*378+68] = 2.0E0*I_ESP_I4x2z_D2y_a-3*I_ESP_G2x2z_D2y;
    abcd[iGrid*378+69] = 2.0E0*I_ESP_I3x3y_D2y_a-2*I_ESP_Gx3y_D2y;
    abcd[iGrid*378+70] = 2.0E0*I_ESP_I3x2yz_D2y_a-2*I_ESP_Gx2yz_D2y;
    abcd[iGrid*378+71] = 2.0E0*I_ESP_I3xy2z_D2y_a-2*I_ESP_Gxy2z_D2y;
    abcd[iGrid*378+72] = 2.0E0*I_ESP_I3x3z_D2y_a-2*I_ESP_Gx3z_D2y;
    abcd[iGrid*378+73] = 2.0E0*I_ESP_I2x4y_D2y_a-1*I_ESP_G4y_D2y;
    abcd[iGrid*378+74] = 2.0E0*I_ESP_I2x3yz_D2y_a-1*I_ESP_G3yz_D2y;
    abcd[iGrid*378+75] = 2.0E0*I_ESP_I2x2y2z_D2y_a-1*I_ESP_G2y2z_D2y;
    abcd[iGrid*378+76] = 2.0E0*I_ESP_I2xy3z_D2y_a-1*I_ESP_Gy3z_D2y;
    abcd[iGrid*378+77] = 2.0E0*I_ESP_I2x4z_D2y_a-1*I_ESP_G4z_D2y;
    abcd[iGrid*378+78] = 2.0E0*I_ESP_Ix5y_D2y_a;
    abcd[iGrid*378+79] = 2.0E0*I_ESP_Ix4yz_D2y_a;
    abcd[iGrid*378+80] = 2.0E0*I_ESP_Ix3y2z_D2y_a;
    abcd[iGrid*378+81] = 2.0E0*I_ESP_Ix2y3z_D2y_a;
    abcd[iGrid*378+82] = 2.0E0*I_ESP_Ixy4z_D2y_a;
    abcd[iGrid*378+83] = 2.0E0*I_ESP_Ix5z_D2y_a;
    abcd[iGrid*378+84] = 2.0E0*I_ESP_I6x_Dyz_a-5*I_ESP_G4x_Dyz;
    abcd[iGrid*378+85] = 2.0E0*I_ESP_I5xy_Dyz_a-4*I_ESP_G3xy_Dyz;
    abcd[iGrid*378+86] = 2.0E0*I_ESP_I5xz_Dyz_a-4*I_ESP_G3xz_Dyz;
    abcd[iGrid*378+87] = 2.0E0*I_ESP_I4x2y_Dyz_a-3*I_ESP_G2x2y_Dyz;
    abcd[iGrid*378+88] = 2.0E0*I_ESP_I4xyz_Dyz_a-3*I_ESP_G2xyz_Dyz;
    abcd[iGrid*378+89] = 2.0E0*I_ESP_I4x2z_Dyz_a-3*I_ESP_G2x2z_Dyz;
    abcd[iGrid*378+90] = 2.0E0*I_ESP_I3x3y_Dyz_a-2*I_ESP_Gx3y_Dyz;
    abcd[iGrid*378+91] = 2.0E0*I_ESP_I3x2yz_Dyz_a-2*I_ESP_Gx2yz_Dyz;
    abcd[iGrid*378+92] = 2.0E0*I_ESP_I3xy2z_Dyz_a-2*I_ESP_Gxy2z_Dyz;
    abcd[iGrid*378+93] = 2.0E0*I_ESP_I3x3z_Dyz_a-2*I_ESP_Gx3z_Dyz;
    abcd[iGrid*378+94] = 2.0E0*I_ESP_I2x4y_Dyz_a-1*I_ESP_G4y_Dyz;
    abcd[iGrid*378+95] = 2.0E0*I_ESP_I2x3yz_Dyz_a-1*I_ESP_G3yz_Dyz;
    abcd[iGrid*378+96] = 2.0E0*I_ESP_I2x2y2z_Dyz_a-1*I_ESP_G2y2z_Dyz;
    abcd[iGrid*378+97] = 2.0E0*I_ESP_I2xy3z_Dyz_a-1*I_ESP_Gy3z_Dyz;
    abcd[iGrid*378+98] = 2.0E0*I_ESP_I2x4z_Dyz_a-1*I_ESP_G4z_Dyz;
    abcd[iGrid*378+99] = 2.0E0*I_ESP_Ix5y_Dyz_a;
    abcd[iGrid*378+100] = 2.0E0*I_ESP_Ix4yz_Dyz_a;
    abcd[iGrid*378+101] = 2.0E0*I_ESP_Ix3y2z_Dyz_a;
    abcd[iGrid*378+102] = 2.0E0*I_ESP_Ix2y3z_Dyz_a;
    abcd[iGrid*378+103] = 2.0E0*I_ESP_Ixy4z_Dyz_a;
    abcd[iGrid*378+104] = 2.0E0*I_ESP_Ix5z_Dyz_a;
    abcd[iGrid*378+105] = 2.0E0*I_ESP_I6x_D2z_a-5*I_ESP_G4x_D2z;
    abcd[iGrid*378+106] = 2.0E0*I_ESP_I5xy_D2z_a-4*I_ESP_G3xy_D2z;
    abcd[iGrid*378+107] = 2.0E0*I_ESP_I5xz_D2z_a-4*I_ESP_G3xz_D2z;
    abcd[iGrid*378+108] = 2.0E0*I_ESP_I4x2y_D2z_a-3*I_ESP_G2x2y_D2z;
    abcd[iGrid*378+109] = 2.0E0*I_ESP_I4xyz_D2z_a-3*I_ESP_G2xyz_D2z;
    abcd[iGrid*378+110] = 2.0E0*I_ESP_I4x2z_D2z_a-3*I_ESP_G2x2z_D2z;
    abcd[iGrid*378+111] = 2.0E0*I_ESP_I3x3y_D2z_a-2*I_ESP_Gx3y_D2z;
    abcd[iGrid*378+112] = 2.0E0*I_ESP_I3x2yz_D2z_a-2*I_ESP_Gx2yz_D2z;
    abcd[iGrid*378+113] = 2.0E0*I_ESP_I3xy2z_D2z_a-2*I_ESP_Gxy2z_D2z;
    abcd[iGrid*378+114] = 2.0E0*I_ESP_I3x3z_D2z_a-2*I_ESP_Gx3z_D2z;
    abcd[iGrid*378+115] = 2.0E0*I_ESP_I2x4y_D2z_a-1*I_ESP_G4y_D2z;
    abcd[iGrid*378+116] = 2.0E0*I_ESP_I2x3yz_D2z_a-1*I_ESP_G3yz_D2z;
    abcd[iGrid*378+117] = 2.0E0*I_ESP_I2x2y2z_D2z_a-1*I_ESP_G2y2z_D2z;
    abcd[iGrid*378+118] = 2.0E0*I_ESP_I2xy3z_D2z_a-1*I_ESP_Gy3z_D2z;
    abcd[iGrid*378+119] = 2.0E0*I_ESP_I2x4z_D2z_a-1*I_ESP_G4z_D2z;
    abcd[iGrid*378+120] = 2.0E0*I_ESP_Ix5y_D2z_a;
    abcd[iGrid*378+121] = 2.0E0*I_ESP_Ix4yz_D2z_a;
    abcd[iGrid*378+122] = 2.0E0*I_ESP_Ix3y2z_D2z_a;
    abcd[iGrid*378+123] = 2.0E0*I_ESP_Ix2y3z_D2z_a;
    abcd[iGrid*378+124] = 2.0E0*I_ESP_Ixy4z_D2z_a;
    abcd[iGrid*378+125] = 2.0E0*I_ESP_Ix5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_a
     * RHS shell quartet name: SQ_ESP_G_D
     ************************************************************/
    abcd[iGrid*378+126] = 2.0E0*I_ESP_I5xy_D2x_a;
    abcd[iGrid*378+127] = 2.0E0*I_ESP_I4x2y_D2x_a-1*I_ESP_G4x_D2x;
    abcd[iGrid*378+128] = 2.0E0*I_ESP_I4xyz_D2x_a;
    abcd[iGrid*378+129] = 2.0E0*I_ESP_I3x3y_D2x_a-2*I_ESP_G3xy_D2x;
    abcd[iGrid*378+130] = 2.0E0*I_ESP_I3x2yz_D2x_a-1*I_ESP_G3xz_D2x;
    abcd[iGrid*378+131] = 2.0E0*I_ESP_I3xy2z_D2x_a;
    abcd[iGrid*378+132] = 2.0E0*I_ESP_I2x4y_D2x_a-3*I_ESP_G2x2y_D2x;
    abcd[iGrid*378+133] = 2.0E0*I_ESP_I2x3yz_D2x_a-2*I_ESP_G2xyz_D2x;
    abcd[iGrid*378+134] = 2.0E0*I_ESP_I2x2y2z_D2x_a-1*I_ESP_G2x2z_D2x;
    abcd[iGrid*378+135] = 2.0E0*I_ESP_I2xy3z_D2x_a;
    abcd[iGrid*378+136] = 2.0E0*I_ESP_Ix5y_D2x_a-4*I_ESP_Gx3y_D2x;
    abcd[iGrid*378+137] = 2.0E0*I_ESP_Ix4yz_D2x_a-3*I_ESP_Gx2yz_D2x;
    abcd[iGrid*378+138] = 2.0E0*I_ESP_Ix3y2z_D2x_a-2*I_ESP_Gxy2z_D2x;
    abcd[iGrid*378+139] = 2.0E0*I_ESP_Ix2y3z_D2x_a-1*I_ESP_Gx3z_D2x;
    abcd[iGrid*378+140] = 2.0E0*I_ESP_Ixy4z_D2x_a;
    abcd[iGrid*378+141] = 2.0E0*I_ESP_I6y_D2x_a-5*I_ESP_G4y_D2x;
    abcd[iGrid*378+142] = 2.0E0*I_ESP_I5yz_D2x_a-4*I_ESP_G3yz_D2x;
    abcd[iGrid*378+143] = 2.0E0*I_ESP_I4y2z_D2x_a-3*I_ESP_G2y2z_D2x;
    abcd[iGrid*378+144] = 2.0E0*I_ESP_I3y3z_D2x_a-2*I_ESP_Gy3z_D2x;
    abcd[iGrid*378+145] = 2.0E0*I_ESP_I2y4z_D2x_a-1*I_ESP_G4z_D2x;
    abcd[iGrid*378+146] = 2.0E0*I_ESP_Iy5z_D2x_a;
    abcd[iGrid*378+147] = 2.0E0*I_ESP_I5xy_Dxy_a;
    abcd[iGrid*378+148] = 2.0E0*I_ESP_I4x2y_Dxy_a-1*I_ESP_G4x_Dxy;
    abcd[iGrid*378+149] = 2.0E0*I_ESP_I4xyz_Dxy_a;
    abcd[iGrid*378+150] = 2.0E0*I_ESP_I3x3y_Dxy_a-2*I_ESP_G3xy_Dxy;
    abcd[iGrid*378+151] = 2.0E0*I_ESP_I3x2yz_Dxy_a-1*I_ESP_G3xz_Dxy;
    abcd[iGrid*378+152] = 2.0E0*I_ESP_I3xy2z_Dxy_a;
    abcd[iGrid*378+153] = 2.0E0*I_ESP_I2x4y_Dxy_a-3*I_ESP_G2x2y_Dxy;
    abcd[iGrid*378+154] = 2.0E0*I_ESP_I2x3yz_Dxy_a-2*I_ESP_G2xyz_Dxy;
    abcd[iGrid*378+155] = 2.0E0*I_ESP_I2x2y2z_Dxy_a-1*I_ESP_G2x2z_Dxy;
    abcd[iGrid*378+156] = 2.0E0*I_ESP_I2xy3z_Dxy_a;
    abcd[iGrid*378+157] = 2.0E0*I_ESP_Ix5y_Dxy_a-4*I_ESP_Gx3y_Dxy;
    abcd[iGrid*378+158] = 2.0E0*I_ESP_Ix4yz_Dxy_a-3*I_ESP_Gx2yz_Dxy;
    abcd[iGrid*378+159] = 2.0E0*I_ESP_Ix3y2z_Dxy_a-2*I_ESP_Gxy2z_Dxy;
    abcd[iGrid*378+160] = 2.0E0*I_ESP_Ix2y3z_Dxy_a-1*I_ESP_Gx3z_Dxy;
    abcd[iGrid*378+161] = 2.0E0*I_ESP_Ixy4z_Dxy_a;
    abcd[iGrid*378+162] = 2.0E0*I_ESP_I6y_Dxy_a-5*I_ESP_G4y_Dxy;
    abcd[iGrid*378+163] = 2.0E0*I_ESP_I5yz_Dxy_a-4*I_ESP_G3yz_Dxy;
    abcd[iGrid*378+164] = 2.0E0*I_ESP_I4y2z_Dxy_a-3*I_ESP_G2y2z_Dxy;
    abcd[iGrid*378+165] = 2.0E0*I_ESP_I3y3z_Dxy_a-2*I_ESP_Gy3z_Dxy;
    abcd[iGrid*378+166] = 2.0E0*I_ESP_I2y4z_Dxy_a-1*I_ESP_G4z_Dxy;
    abcd[iGrid*378+167] = 2.0E0*I_ESP_Iy5z_Dxy_a;
    abcd[iGrid*378+168] = 2.0E0*I_ESP_I5xy_Dxz_a;
    abcd[iGrid*378+169] = 2.0E0*I_ESP_I4x2y_Dxz_a-1*I_ESP_G4x_Dxz;
    abcd[iGrid*378+170] = 2.0E0*I_ESP_I4xyz_Dxz_a;
    abcd[iGrid*378+171] = 2.0E0*I_ESP_I3x3y_Dxz_a-2*I_ESP_G3xy_Dxz;
    abcd[iGrid*378+172] = 2.0E0*I_ESP_I3x2yz_Dxz_a-1*I_ESP_G3xz_Dxz;
    abcd[iGrid*378+173] = 2.0E0*I_ESP_I3xy2z_Dxz_a;
    abcd[iGrid*378+174] = 2.0E0*I_ESP_I2x4y_Dxz_a-3*I_ESP_G2x2y_Dxz;
    abcd[iGrid*378+175] = 2.0E0*I_ESP_I2x3yz_Dxz_a-2*I_ESP_G2xyz_Dxz;
    abcd[iGrid*378+176] = 2.0E0*I_ESP_I2x2y2z_Dxz_a-1*I_ESP_G2x2z_Dxz;
    abcd[iGrid*378+177] = 2.0E0*I_ESP_I2xy3z_Dxz_a;
    abcd[iGrid*378+178] = 2.0E0*I_ESP_Ix5y_Dxz_a-4*I_ESP_Gx3y_Dxz;
    abcd[iGrid*378+179] = 2.0E0*I_ESP_Ix4yz_Dxz_a-3*I_ESP_Gx2yz_Dxz;
    abcd[iGrid*378+180] = 2.0E0*I_ESP_Ix3y2z_Dxz_a-2*I_ESP_Gxy2z_Dxz;
    abcd[iGrid*378+181] = 2.0E0*I_ESP_Ix2y3z_Dxz_a-1*I_ESP_Gx3z_Dxz;
    abcd[iGrid*378+182] = 2.0E0*I_ESP_Ixy4z_Dxz_a;
    abcd[iGrid*378+183] = 2.0E0*I_ESP_I6y_Dxz_a-5*I_ESP_G4y_Dxz;
    abcd[iGrid*378+184] = 2.0E0*I_ESP_I5yz_Dxz_a-4*I_ESP_G3yz_Dxz;
    abcd[iGrid*378+185] = 2.0E0*I_ESP_I4y2z_Dxz_a-3*I_ESP_G2y2z_Dxz;
    abcd[iGrid*378+186] = 2.0E0*I_ESP_I3y3z_Dxz_a-2*I_ESP_Gy3z_Dxz;
    abcd[iGrid*378+187] = 2.0E0*I_ESP_I2y4z_Dxz_a-1*I_ESP_G4z_Dxz;
    abcd[iGrid*378+188] = 2.0E0*I_ESP_Iy5z_Dxz_a;
    abcd[iGrid*378+189] = 2.0E0*I_ESP_I5xy_D2y_a;
    abcd[iGrid*378+190] = 2.0E0*I_ESP_I4x2y_D2y_a-1*I_ESP_G4x_D2y;
    abcd[iGrid*378+191] = 2.0E0*I_ESP_I4xyz_D2y_a;
    abcd[iGrid*378+192] = 2.0E0*I_ESP_I3x3y_D2y_a-2*I_ESP_G3xy_D2y;
    abcd[iGrid*378+193] = 2.0E0*I_ESP_I3x2yz_D2y_a-1*I_ESP_G3xz_D2y;
    abcd[iGrid*378+194] = 2.0E0*I_ESP_I3xy2z_D2y_a;
    abcd[iGrid*378+195] = 2.0E0*I_ESP_I2x4y_D2y_a-3*I_ESP_G2x2y_D2y;
    abcd[iGrid*378+196] = 2.0E0*I_ESP_I2x3yz_D2y_a-2*I_ESP_G2xyz_D2y;
    abcd[iGrid*378+197] = 2.0E0*I_ESP_I2x2y2z_D2y_a-1*I_ESP_G2x2z_D2y;
    abcd[iGrid*378+198] = 2.0E0*I_ESP_I2xy3z_D2y_a;
    abcd[iGrid*378+199] = 2.0E0*I_ESP_Ix5y_D2y_a-4*I_ESP_Gx3y_D2y;
    abcd[iGrid*378+200] = 2.0E0*I_ESP_Ix4yz_D2y_a-3*I_ESP_Gx2yz_D2y;
    abcd[iGrid*378+201] = 2.0E0*I_ESP_Ix3y2z_D2y_a-2*I_ESP_Gxy2z_D2y;
    abcd[iGrid*378+202] = 2.0E0*I_ESP_Ix2y3z_D2y_a-1*I_ESP_Gx3z_D2y;
    abcd[iGrid*378+203] = 2.0E0*I_ESP_Ixy4z_D2y_a;
    abcd[iGrid*378+204] = 2.0E0*I_ESP_I6y_D2y_a-5*I_ESP_G4y_D2y;
    abcd[iGrid*378+205] = 2.0E0*I_ESP_I5yz_D2y_a-4*I_ESP_G3yz_D2y;
    abcd[iGrid*378+206] = 2.0E0*I_ESP_I4y2z_D2y_a-3*I_ESP_G2y2z_D2y;
    abcd[iGrid*378+207] = 2.0E0*I_ESP_I3y3z_D2y_a-2*I_ESP_Gy3z_D2y;
    abcd[iGrid*378+208] = 2.0E0*I_ESP_I2y4z_D2y_a-1*I_ESP_G4z_D2y;
    abcd[iGrid*378+209] = 2.0E0*I_ESP_Iy5z_D2y_a;
    abcd[iGrid*378+210] = 2.0E0*I_ESP_I5xy_Dyz_a;
    abcd[iGrid*378+211] = 2.0E0*I_ESP_I4x2y_Dyz_a-1*I_ESP_G4x_Dyz;
    abcd[iGrid*378+212] = 2.0E0*I_ESP_I4xyz_Dyz_a;
    abcd[iGrid*378+213] = 2.0E0*I_ESP_I3x3y_Dyz_a-2*I_ESP_G3xy_Dyz;
    abcd[iGrid*378+214] = 2.0E0*I_ESP_I3x2yz_Dyz_a-1*I_ESP_G3xz_Dyz;
    abcd[iGrid*378+215] = 2.0E0*I_ESP_I3xy2z_Dyz_a;
    abcd[iGrid*378+216] = 2.0E0*I_ESP_I2x4y_Dyz_a-3*I_ESP_G2x2y_Dyz;
    abcd[iGrid*378+217] = 2.0E0*I_ESP_I2x3yz_Dyz_a-2*I_ESP_G2xyz_Dyz;
    abcd[iGrid*378+218] = 2.0E0*I_ESP_I2x2y2z_Dyz_a-1*I_ESP_G2x2z_Dyz;
    abcd[iGrid*378+219] = 2.0E0*I_ESP_I2xy3z_Dyz_a;
    abcd[iGrid*378+220] = 2.0E0*I_ESP_Ix5y_Dyz_a-4*I_ESP_Gx3y_Dyz;
    abcd[iGrid*378+221] = 2.0E0*I_ESP_Ix4yz_Dyz_a-3*I_ESP_Gx2yz_Dyz;
    abcd[iGrid*378+222] = 2.0E0*I_ESP_Ix3y2z_Dyz_a-2*I_ESP_Gxy2z_Dyz;
    abcd[iGrid*378+223] = 2.0E0*I_ESP_Ix2y3z_Dyz_a-1*I_ESP_Gx3z_Dyz;
    abcd[iGrid*378+224] = 2.0E0*I_ESP_Ixy4z_Dyz_a;
    abcd[iGrid*378+225] = 2.0E0*I_ESP_I6y_Dyz_a-5*I_ESP_G4y_Dyz;
    abcd[iGrid*378+226] = 2.0E0*I_ESP_I5yz_Dyz_a-4*I_ESP_G3yz_Dyz;
    abcd[iGrid*378+227] = 2.0E0*I_ESP_I4y2z_Dyz_a-3*I_ESP_G2y2z_Dyz;
    abcd[iGrid*378+228] = 2.0E0*I_ESP_I3y3z_Dyz_a-2*I_ESP_Gy3z_Dyz;
    abcd[iGrid*378+229] = 2.0E0*I_ESP_I2y4z_Dyz_a-1*I_ESP_G4z_Dyz;
    abcd[iGrid*378+230] = 2.0E0*I_ESP_Iy5z_Dyz_a;
    abcd[iGrid*378+231] = 2.0E0*I_ESP_I5xy_D2z_a;
    abcd[iGrid*378+232] = 2.0E0*I_ESP_I4x2y_D2z_a-1*I_ESP_G4x_D2z;
    abcd[iGrid*378+233] = 2.0E0*I_ESP_I4xyz_D2z_a;
    abcd[iGrid*378+234] = 2.0E0*I_ESP_I3x3y_D2z_a-2*I_ESP_G3xy_D2z;
    abcd[iGrid*378+235] = 2.0E0*I_ESP_I3x2yz_D2z_a-1*I_ESP_G3xz_D2z;
    abcd[iGrid*378+236] = 2.0E0*I_ESP_I3xy2z_D2z_a;
    abcd[iGrid*378+237] = 2.0E0*I_ESP_I2x4y_D2z_a-3*I_ESP_G2x2y_D2z;
    abcd[iGrid*378+238] = 2.0E0*I_ESP_I2x3yz_D2z_a-2*I_ESP_G2xyz_D2z;
    abcd[iGrid*378+239] = 2.0E0*I_ESP_I2x2y2z_D2z_a-1*I_ESP_G2x2z_D2z;
    abcd[iGrid*378+240] = 2.0E0*I_ESP_I2xy3z_D2z_a;
    abcd[iGrid*378+241] = 2.0E0*I_ESP_Ix5y_D2z_a-4*I_ESP_Gx3y_D2z;
    abcd[iGrid*378+242] = 2.0E0*I_ESP_Ix4yz_D2z_a-3*I_ESP_Gx2yz_D2z;
    abcd[iGrid*378+243] = 2.0E0*I_ESP_Ix3y2z_D2z_a-2*I_ESP_Gxy2z_D2z;
    abcd[iGrid*378+244] = 2.0E0*I_ESP_Ix2y3z_D2z_a-1*I_ESP_Gx3z_D2z;
    abcd[iGrid*378+245] = 2.0E0*I_ESP_Ixy4z_D2z_a;
    abcd[iGrid*378+246] = 2.0E0*I_ESP_I6y_D2z_a-5*I_ESP_G4y_D2z;
    abcd[iGrid*378+247] = 2.0E0*I_ESP_I5yz_D2z_a-4*I_ESP_G3yz_D2z;
    abcd[iGrid*378+248] = 2.0E0*I_ESP_I4y2z_D2z_a-3*I_ESP_G2y2z_D2z;
    abcd[iGrid*378+249] = 2.0E0*I_ESP_I3y3z_D2z_a-2*I_ESP_Gy3z_D2z;
    abcd[iGrid*378+250] = 2.0E0*I_ESP_I2y4z_D2z_a-1*I_ESP_G4z_D2z;
    abcd[iGrid*378+251] = 2.0E0*I_ESP_Iy5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_a
     * RHS shell quartet name: SQ_ESP_G_D
     ************************************************************/
    abcd[iGrid*378+252] = 2.0E0*I_ESP_I5xz_D2x_a;
    abcd[iGrid*378+253] = 2.0E0*I_ESP_I4xyz_D2x_a;
    abcd[iGrid*378+254] = 2.0E0*I_ESP_I4x2z_D2x_a-1*I_ESP_G4x_D2x;
    abcd[iGrid*378+255] = 2.0E0*I_ESP_I3x2yz_D2x_a;
    abcd[iGrid*378+256] = 2.0E0*I_ESP_I3xy2z_D2x_a-1*I_ESP_G3xy_D2x;
    abcd[iGrid*378+257] = 2.0E0*I_ESP_I3x3z_D2x_a-2*I_ESP_G3xz_D2x;
    abcd[iGrid*378+258] = 2.0E0*I_ESP_I2x3yz_D2x_a;
    abcd[iGrid*378+259] = 2.0E0*I_ESP_I2x2y2z_D2x_a-1*I_ESP_G2x2y_D2x;
    abcd[iGrid*378+260] = 2.0E0*I_ESP_I2xy3z_D2x_a-2*I_ESP_G2xyz_D2x;
    abcd[iGrid*378+261] = 2.0E0*I_ESP_I2x4z_D2x_a-3*I_ESP_G2x2z_D2x;
    abcd[iGrid*378+262] = 2.0E0*I_ESP_Ix4yz_D2x_a;
    abcd[iGrid*378+263] = 2.0E0*I_ESP_Ix3y2z_D2x_a-1*I_ESP_Gx3y_D2x;
    abcd[iGrid*378+264] = 2.0E0*I_ESP_Ix2y3z_D2x_a-2*I_ESP_Gx2yz_D2x;
    abcd[iGrid*378+265] = 2.0E0*I_ESP_Ixy4z_D2x_a-3*I_ESP_Gxy2z_D2x;
    abcd[iGrid*378+266] = 2.0E0*I_ESP_Ix5z_D2x_a-4*I_ESP_Gx3z_D2x;
    abcd[iGrid*378+267] = 2.0E0*I_ESP_I5yz_D2x_a;
    abcd[iGrid*378+268] = 2.0E0*I_ESP_I4y2z_D2x_a-1*I_ESP_G4y_D2x;
    abcd[iGrid*378+269] = 2.0E0*I_ESP_I3y3z_D2x_a-2*I_ESP_G3yz_D2x;
    abcd[iGrid*378+270] = 2.0E0*I_ESP_I2y4z_D2x_a-3*I_ESP_G2y2z_D2x;
    abcd[iGrid*378+271] = 2.0E0*I_ESP_Iy5z_D2x_a-4*I_ESP_Gy3z_D2x;
    abcd[iGrid*378+272] = 2.0E0*I_ESP_I6z_D2x_a-5*I_ESP_G4z_D2x;
    abcd[iGrid*378+273] = 2.0E0*I_ESP_I5xz_Dxy_a;
    abcd[iGrid*378+274] = 2.0E0*I_ESP_I4xyz_Dxy_a;
    abcd[iGrid*378+275] = 2.0E0*I_ESP_I4x2z_Dxy_a-1*I_ESP_G4x_Dxy;
    abcd[iGrid*378+276] = 2.0E0*I_ESP_I3x2yz_Dxy_a;
    abcd[iGrid*378+277] = 2.0E0*I_ESP_I3xy2z_Dxy_a-1*I_ESP_G3xy_Dxy;
    abcd[iGrid*378+278] = 2.0E0*I_ESP_I3x3z_Dxy_a-2*I_ESP_G3xz_Dxy;
    abcd[iGrid*378+279] = 2.0E0*I_ESP_I2x3yz_Dxy_a;
    abcd[iGrid*378+280] = 2.0E0*I_ESP_I2x2y2z_Dxy_a-1*I_ESP_G2x2y_Dxy;
    abcd[iGrid*378+281] = 2.0E0*I_ESP_I2xy3z_Dxy_a-2*I_ESP_G2xyz_Dxy;
    abcd[iGrid*378+282] = 2.0E0*I_ESP_I2x4z_Dxy_a-3*I_ESP_G2x2z_Dxy;
    abcd[iGrid*378+283] = 2.0E0*I_ESP_Ix4yz_Dxy_a;
    abcd[iGrid*378+284] = 2.0E0*I_ESP_Ix3y2z_Dxy_a-1*I_ESP_Gx3y_Dxy;
    abcd[iGrid*378+285] = 2.0E0*I_ESP_Ix2y3z_Dxy_a-2*I_ESP_Gx2yz_Dxy;
    abcd[iGrid*378+286] = 2.0E0*I_ESP_Ixy4z_Dxy_a-3*I_ESP_Gxy2z_Dxy;
    abcd[iGrid*378+287] = 2.0E0*I_ESP_Ix5z_Dxy_a-4*I_ESP_Gx3z_Dxy;
    abcd[iGrid*378+288] = 2.0E0*I_ESP_I5yz_Dxy_a;
    abcd[iGrid*378+289] = 2.0E0*I_ESP_I4y2z_Dxy_a-1*I_ESP_G4y_Dxy;
    abcd[iGrid*378+290] = 2.0E0*I_ESP_I3y3z_Dxy_a-2*I_ESP_G3yz_Dxy;
    abcd[iGrid*378+291] = 2.0E0*I_ESP_I2y4z_Dxy_a-3*I_ESP_G2y2z_Dxy;
    abcd[iGrid*378+292] = 2.0E0*I_ESP_Iy5z_Dxy_a-4*I_ESP_Gy3z_Dxy;
    abcd[iGrid*378+293] = 2.0E0*I_ESP_I6z_Dxy_a-5*I_ESP_G4z_Dxy;
    abcd[iGrid*378+294] = 2.0E0*I_ESP_I5xz_Dxz_a;
    abcd[iGrid*378+295] = 2.0E0*I_ESP_I4xyz_Dxz_a;
    abcd[iGrid*378+296] = 2.0E0*I_ESP_I4x2z_Dxz_a-1*I_ESP_G4x_Dxz;
    abcd[iGrid*378+297] = 2.0E0*I_ESP_I3x2yz_Dxz_a;
    abcd[iGrid*378+298] = 2.0E0*I_ESP_I3xy2z_Dxz_a-1*I_ESP_G3xy_Dxz;
    abcd[iGrid*378+299] = 2.0E0*I_ESP_I3x3z_Dxz_a-2*I_ESP_G3xz_Dxz;
    abcd[iGrid*378+300] = 2.0E0*I_ESP_I2x3yz_Dxz_a;
    abcd[iGrid*378+301] = 2.0E0*I_ESP_I2x2y2z_Dxz_a-1*I_ESP_G2x2y_Dxz;
    abcd[iGrid*378+302] = 2.0E0*I_ESP_I2xy3z_Dxz_a-2*I_ESP_G2xyz_Dxz;
    abcd[iGrid*378+303] = 2.0E0*I_ESP_I2x4z_Dxz_a-3*I_ESP_G2x2z_Dxz;
    abcd[iGrid*378+304] = 2.0E0*I_ESP_Ix4yz_Dxz_a;
    abcd[iGrid*378+305] = 2.0E0*I_ESP_Ix3y2z_Dxz_a-1*I_ESP_Gx3y_Dxz;
    abcd[iGrid*378+306] = 2.0E0*I_ESP_Ix2y3z_Dxz_a-2*I_ESP_Gx2yz_Dxz;
    abcd[iGrid*378+307] = 2.0E0*I_ESP_Ixy4z_Dxz_a-3*I_ESP_Gxy2z_Dxz;
    abcd[iGrid*378+308] = 2.0E0*I_ESP_Ix5z_Dxz_a-4*I_ESP_Gx3z_Dxz;
    abcd[iGrid*378+309] = 2.0E0*I_ESP_I5yz_Dxz_a;
    abcd[iGrid*378+310] = 2.0E0*I_ESP_I4y2z_Dxz_a-1*I_ESP_G4y_Dxz;
    abcd[iGrid*378+311] = 2.0E0*I_ESP_I3y3z_Dxz_a-2*I_ESP_G3yz_Dxz;
    abcd[iGrid*378+312] = 2.0E0*I_ESP_I2y4z_Dxz_a-3*I_ESP_G2y2z_Dxz;
    abcd[iGrid*378+313] = 2.0E0*I_ESP_Iy5z_Dxz_a-4*I_ESP_Gy3z_Dxz;
    abcd[iGrid*378+314] = 2.0E0*I_ESP_I6z_Dxz_a-5*I_ESP_G4z_Dxz;
    abcd[iGrid*378+315] = 2.0E0*I_ESP_I5xz_D2y_a;
    abcd[iGrid*378+316] = 2.0E0*I_ESP_I4xyz_D2y_a;
    abcd[iGrid*378+317] = 2.0E0*I_ESP_I4x2z_D2y_a-1*I_ESP_G4x_D2y;
    abcd[iGrid*378+318] = 2.0E0*I_ESP_I3x2yz_D2y_a;
    abcd[iGrid*378+319] = 2.0E0*I_ESP_I3xy2z_D2y_a-1*I_ESP_G3xy_D2y;
    abcd[iGrid*378+320] = 2.0E0*I_ESP_I3x3z_D2y_a-2*I_ESP_G3xz_D2y;
    abcd[iGrid*378+321] = 2.0E0*I_ESP_I2x3yz_D2y_a;
    abcd[iGrid*378+322] = 2.0E0*I_ESP_I2x2y2z_D2y_a-1*I_ESP_G2x2y_D2y;
    abcd[iGrid*378+323] = 2.0E0*I_ESP_I2xy3z_D2y_a-2*I_ESP_G2xyz_D2y;
    abcd[iGrid*378+324] = 2.0E0*I_ESP_I2x4z_D2y_a-3*I_ESP_G2x2z_D2y;
    abcd[iGrid*378+325] = 2.0E0*I_ESP_Ix4yz_D2y_a;
    abcd[iGrid*378+326] = 2.0E0*I_ESP_Ix3y2z_D2y_a-1*I_ESP_Gx3y_D2y;
    abcd[iGrid*378+327] = 2.0E0*I_ESP_Ix2y3z_D2y_a-2*I_ESP_Gx2yz_D2y;
    abcd[iGrid*378+328] = 2.0E0*I_ESP_Ixy4z_D2y_a-3*I_ESP_Gxy2z_D2y;
    abcd[iGrid*378+329] = 2.0E0*I_ESP_Ix5z_D2y_a-4*I_ESP_Gx3z_D2y;
    abcd[iGrid*378+330] = 2.0E0*I_ESP_I5yz_D2y_a;
    abcd[iGrid*378+331] = 2.0E0*I_ESP_I4y2z_D2y_a-1*I_ESP_G4y_D2y;
    abcd[iGrid*378+332] = 2.0E0*I_ESP_I3y3z_D2y_a-2*I_ESP_G3yz_D2y;
    abcd[iGrid*378+333] = 2.0E0*I_ESP_I2y4z_D2y_a-3*I_ESP_G2y2z_D2y;
    abcd[iGrid*378+334] = 2.0E0*I_ESP_Iy5z_D2y_a-4*I_ESP_Gy3z_D2y;
    abcd[iGrid*378+335] = 2.0E0*I_ESP_I6z_D2y_a-5*I_ESP_G4z_D2y;
    abcd[iGrid*378+336] = 2.0E0*I_ESP_I5xz_Dyz_a;
    abcd[iGrid*378+337] = 2.0E0*I_ESP_I4xyz_Dyz_a;
    abcd[iGrid*378+338] = 2.0E0*I_ESP_I4x2z_Dyz_a-1*I_ESP_G4x_Dyz;
    abcd[iGrid*378+339] = 2.0E0*I_ESP_I3x2yz_Dyz_a;
    abcd[iGrid*378+340] = 2.0E0*I_ESP_I3xy2z_Dyz_a-1*I_ESP_G3xy_Dyz;
    abcd[iGrid*378+341] = 2.0E0*I_ESP_I3x3z_Dyz_a-2*I_ESP_G3xz_Dyz;
    abcd[iGrid*378+342] = 2.0E0*I_ESP_I2x3yz_Dyz_a;
    abcd[iGrid*378+343] = 2.0E0*I_ESP_I2x2y2z_Dyz_a-1*I_ESP_G2x2y_Dyz;
    abcd[iGrid*378+344] = 2.0E0*I_ESP_I2xy3z_Dyz_a-2*I_ESP_G2xyz_Dyz;
    abcd[iGrid*378+345] = 2.0E0*I_ESP_I2x4z_Dyz_a-3*I_ESP_G2x2z_Dyz;
    abcd[iGrid*378+346] = 2.0E0*I_ESP_Ix4yz_Dyz_a;
    abcd[iGrid*378+347] = 2.0E0*I_ESP_Ix3y2z_Dyz_a-1*I_ESP_Gx3y_Dyz;
    abcd[iGrid*378+348] = 2.0E0*I_ESP_Ix2y3z_Dyz_a-2*I_ESP_Gx2yz_Dyz;
    abcd[iGrid*378+349] = 2.0E0*I_ESP_Ixy4z_Dyz_a-3*I_ESP_Gxy2z_Dyz;
    abcd[iGrid*378+350] = 2.0E0*I_ESP_Ix5z_Dyz_a-4*I_ESP_Gx3z_Dyz;
    abcd[iGrid*378+351] = 2.0E0*I_ESP_I5yz_Dyz_a;
    abcd[iGrid*378+352] = 2.0E0*I_ESP_I4y2z_Dyz_a-1*I_ESP_G4y_Dyz;
    abcd[iGrid*378+353] = 2.0E0*I_ESP_I3y3z_Dyz_a-2*I_ESP_G3yz_Dyz;
    abcd[iGrid*378+354] = 2.0E0*I_ESP_I2y4z_Dyz_a-3*I_ESP_G2y2z_Dyz;
    abcd[iGrid*378+355] = 2.0E0*I_ESP_Iy5z_Dyz_a-4*I_ESP_Gy3z_Dyz;
    abcd[iGrid*378+356] = 2.0E0*I_ESP_I6z_Dyz_a-5*I_ESP_G4z_Dyz;
    abcd[iGrid*378+357] = 2.0E0*I_ESP_I5xz_D2z_a;
    abcd[iGrid*378+358] = 2.0E0*I_ESP_I4xyz_D2z_a;
    abcd[iGrid*378+359] = 2.0E0*I_ESP_I4x2z_D2z_a-1*I_ESP_G4x_D2z;
    abcd[iGrid*378+360] = 2.0E0*I_ESP_I3x2yz_D2z_a;
    abcd[iGrid*378+361] = 2.0E0*I_ESP_I3xy2z_D2z_a-1*I_ESP_G3xy_D2z;
    abcd[iGrid*378+362] = 2.0E0*I_ESP_I3x3z_D2z_a-2*I_ESP_G3xz_D2z;
    abcd[iGrid*378+363] = 2.0E0*I_ESP_I2x3yz_D2z_a;
    abcd[iGrid*378+364] = 2.0E0*I_ESP_I2x2y2z_D2z_a-1*I_ESP_G2x2y_D2z;
    abcd[iGrid*378+365] = 2.0E0*I_ESP_I2xy3z_D2z_a-2*I_ESP_G2xyz_D2z;
    abcd[iGrid*378+366] = 2.0E0*I_ESP_I2x4z_D2z_a-3*I_ESP_G2x2z_D2z;
    abcd[iGrid*378+367] = 2.0E0*I_ESP_Ix4yz_D2z_a;
    abcd[iGrid*378+368] = 2.0E0*I_ESP_Ix3y2z_D2z_a-1*I_ESP_Gx3y_D2z;
    abcd[iGrid*378+369] = 2.0E0*I_ESP_Ix2y3z_D2z_a-2*I_ESP_Gx2yz_D2z;
    abcd[iGrid*378+370] = 2.0E0*I_ESP_Ixy4z_D2z_a-3*I_ESP_Gxy2z_D2z;
    abcd[iGrid*378+371] = 2.0E0*I_ESP_Ix5z_D2z_a-4*I_ESP_Gx3z_D2z;
    abcd[iGrid*378+372] = 2.0E0*I_ESP_I5yz_D2z_a;
    abcd[iGrid*378+373] = 2.0E0*I_ESP_I4y2z_D2z_a-1*I_ESP_G4y_D2z;
    abcd[iGrid*378+374] = 2.0E0*I_ESP_I3y3z_D2z_a-2*I_ESP_G3yz_D2z;
    abcd[iGrid*378+375] = 2.0E0*I_ESP_I2y4z_D2z_a-3*I_ESP_G2y2z_D2z;
    abcd[iGrid*378+376] = 2.0E0*I_ESP_Iy5z_D2z_a-4*I_ESP_Gy3z_D2z;
    abcd[iGrid*378+377] = 2.0E0*I_ESP_I6z_D2z_a-5*I_ESP_G4z_D2z;
  }
}
