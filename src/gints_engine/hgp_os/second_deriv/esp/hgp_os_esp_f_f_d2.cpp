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
// BRA1 as redundant position, total RHS integrals evaluated as: 0
// BRA2 as redundant position, total RHS integrals evaluated as: 0
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: NOT AVIALABLE
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

void hgp_os_esp_f_f_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
    Double I_ESP_H5x_S_aa = 0.0E0;
    Double I_ESP_H4xy_S_aa = 0.0E0;
    Double I_ESP_H4xz_S_aa = 0.0E0;
    Double I_ESP_H3x2y_S_aa = 0.0E0;
    Double I_ESP_H3xyz_S_aa = 0.0E0;
    Double I_ESP_H3x2z_S_aa = 0.0E0;
    Double I_ESP_H2x3y_S_aa = 0.0E0;
    Double I_ESP_H2x2yz_S_aa = 0.0E0;
    Double I_ESP_H2xy2z_S_aa = 0.0E0;
    Double I_ESP_H2x3z_S_aa = 0.0E0;
    Double I_ESP_Hx4y_S_aa = 0.0E0;
    Double I_ESP_Hx3yz_S_aa = 0.0E0;
    Double I_ESP_Hx2y2z_S_aa = 0.0E0;
    Double I_ESP_Hxy3z_S_aa = 0.0E0;
    Double I_ESP_Hx4z_S_aa = 0.0E0;
    Double I_ESP_H5y_S_aa = 0.0E0;
    Double I_ESP_H4yz_S_aa = 0.0E0;
    Double I_ESP_H3y2z_S_aa = 0.0E0;
    Double I_ESP_H2y3z_S_aa = 0.0E0;
    Double I_ESP_Hy4z_S_aa = 0.0E0;
    Double I_ESP_H5z_S_aa = 0.0E0;
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
    Double I_ESP_F3x_S_a = 0.0E0;
    Double I_ESP_F2xy_S_a = 0.0E0;
    Double I_ESP_F2xz_S_a = 0.0E0;
    Double I_ESP_Fx2y_S_a = 0.0E0;
    Double I_ESP_Fxyz_S_a = 0.0E0;
    Double I_ESP_Fx2z_S_a = 0.0E0;
    Double I_ESP_F3y_S_a = 0.0E0;
    Double I_ESP_F2yz_S_a = 0.0E0;
    Double I_ESP_Fy2z_S_a = 0.0E0;
    Double I_ESP_F3z_S_a = 0.0E0;
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
    Double I_ESP_Px_S = 0.0E0;
    Double I_ESP_Py_S = 0.0E0;
    Double I_ESP_Pz_S = 0.0E0;

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
         * shell quartet name: SQ_ESP_H_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_aa_coefs = alpha*alpha;
        I_ESP_H5x_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5z_S_vrr;

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
         * shell quartet name: SQ_ESP_F_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_F_S_a_coefs = alpha;
        I_ESP_F3x_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3z_S_vrr;

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

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_Px_S += I_ESP_Px_S_vrr;
        I_ESP_Py_S += I_ESP_Py_S_vrr;
        I_ESP_Pz_S += I_ESP_Pz_S_vrr;
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
     * shell quartet name: SQ_ESP_P_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_D_S
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    Double I_ESP_Px_Px = I_ESP_D2x_S+ABX*I_ESP_Px_S;
    Double I_ESP_Py_Px = I_ESP_Dxy_S+ABX*I_ESP_Py_S;
    Double I_ESP_Pz_Px = I_ESP_Dxz_S+ABX*I_ESP_Pz_S;
    Double I_ESP_Px_Py = I_ESP_Dxy_S+ABY*I_ESP_Px_S;
    Double I_ESP_Py_Py = I_ESP_D2y_S+ABY*I_ESP_Py_S;
    Double I_ESP_Pz_Py = I_ESP_Dyz_S+ABY*I_ESP_Pz_S;
    Double I_ESP_Px_Pz = I_ESP_Dxz_S+ABZ*I_ESP_Px_S;
    Double I_ESP_Py_Pz = I_ESP_Dyz_S+ABZ*I_ESP_Py_S;
    Double I_ESP_Pz_Pz = I_ESP_D2z_S+ABZ*I_ESP_Pz_S;

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
     * shell quartet name: SQ_ESP_P_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_D_P
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    Double I_ESP_Px_D2x = I_ESP_D2x_Px+ABX*I_ESP_Px_Px;
    Double I_ESP_Py_D2x = I_ESP_Dxy_Px+ABX*I_ESP_Py_Px;
    Double I_ESP_Pz_D2x = I_ESP_Dxz_Px+ABX*I_ESP_Pz_Px;
    Double I_ESP_Px_Dxy = I_ESP_Dxy_Px+ABY*I_ESP_Px_Px;
    Double I_ESP_Py_Dxy = I_ESP_D2y_Px+ABY*I_ESP_Py_Px;
    Double I_ESP_Pz_Dxy = I_ESP_Dyz_Px+ABY*I_ESP_Pz_Px;
    Double I_ESP_Px_D2y = I_ESP_Dxy_Py+ABY*I_ESP_Px_Py;
    Double I_ESP_Py_D2y = I_ESP_D2y_Py+ABY*I_ESP_Py_Py;
    Double I_ESP_Pz_D2y = I_ESP_Dyz_Py+ABY*I_ESP_Pz_Py;
    Double I_ESP_Px_D2z = I_ESP_Dxz_Pz+ABZ*I_ESP_Px_Pz;
    Double I_ESP_Py_D2z = I_ESP_Dyz_Pz+ABZ*I_ESP_Py_Pz;
    Double I_ESP_Pz_D2z = I_ESP_D2z_Pz+ABZ*I_ESP_Pz_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 10 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_S
     * RHS shell quartet name: SQ_ESP_F_S
     ************************************************************/
    Double I_ESP_F3x_Px = I_ESP_G4x_S+ABX*I_ESP_F3x_S;
    Double I_ESP_F2xy_Px = I_ESP_G3xy_S+ABX*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Px = I_ESP_G3xz_S+ABX*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Px = I_ESP_G2x2y_S+ABX*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Px = I_ESP_G2xyz_S+ABX*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Px = I_ESP_G2x2z_S+ABX*I_ESP_Fx2z_S;
    Double I_ESP_F2yz_Px = I_ESP_Gx2yz_S+ABX*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Px = I_ESP_Gxy2z_S+ABX*I_ESP_Fy2z_S;
    Double I_ESP_F2xy_Py = I_ESP_G2x2y_S+ABY*I_ESP_F2xy_S;
    Double I_ESP_Fx2y_Py = I_ESP_Gx3y_S+ABY*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Py = I_ESP_Gx2yz_S+ABY*I_ESP_Fxyz_S;
    Double I_ESP_F3y_Py = I_ESP_G4y_S+ABY*I_ESP_F3y_S;
    Double I_ESP_F2yz_Py = I_ESP_G3yz_S+ABY*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Py = I_ESP_G2y2z_S+ABY*I_ESP_Fy2z_S;
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
     * totally 15 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_P
     * RHS shell quartet name: SQ_ESP_D_P
     ************************************************************/
    Double I_ESP_D2x_D2x = I_ESP_F3x_Px+ABX*I_ESP_D2x_Px;
    Double I_ESP_Dxy_D2x = I_ESP_F2xy_Px+ABX*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_D2x = I_ESP_F2xz_Px+ABX*I_ESP_Dxz_Px;
    Double I_ESP_D2y_D2x = I_ESP_Fx2y_Px+ABX*I_ESP_D2y_Px;
    Double I_ESP_Dyz_D2x = I_ESP_Fxyz_Px+ABX*I_ESP_Dyz_Px;
    Double I_ESP_D2z_D2x = I_ESP_Fx2z_Px+ABX*I_ESP_D2z_Px;
    Double I_ESP_Dxz_Dxy = I_ESP_Fxyz_Px+ABY*I_ESP_Dxz_Px;
    Double I_ESP_Dyz_Dxy = I_ESP_F2yz_Px+ABY*I_ESP_Dyz_Px;
    Double I_ESP_D2z_Dxy = I_ESP_Fy2z_Px+ABY*I_ESP_D2z_Px;
    Double I_ESP_D2x_D2y = I_ESP_F2xy_Py+ABY*I_ESP_D2x_Py;
    Double I_ESP_Dxy_D2y = I_ESP_Fx2y_Py+ABY*I_ESP_Dxy_Py;
    Double I_ESP_Dxz_D2y = I_ESP_Fxyz_Py+ABY*I_ESP_Dxz_Py;
    Double I_ESP_D2y_D2y = I_ESP_F3y_Py+ABY*I_ESP_D2y_Py;
    Double I_ESP_Dyz_D2y = I_ESP_F2yz_Py+ABY*I_ESP_Dyz_Py;
    Double I_ESP_D2z_D2y = I_ESP_Fy2z_Py+ABY*I_ESP_D2z_Py;
    Double I_ESP_D2x_D2z = I_ESP_F2xz_Pz+ABZ*I_ESP_D2x_Pz;
    Double I_ESP_Dxy_D2z = I_ESP_Fxyz_Pz+ABZ*I_ESP_Dxy_Pz;
    Double I_ESP_Dxz_D2z = I_ESP_Fx2z_Pz+ABZ*I_ESP_Dxz_Pz;
    Double I_ESP_D2y_D2z = I_ESP_F2yz_Pz+ABZ*I_ESP_D2y_Pz;
    Double I_ESP_Dyz_D2z = I_ESP_Fy2z_Pz+ABZ*I_ESP_Dyz_Pz;
    Double I_ESP_D2z_D2z = I_ESP_F3z_Pz+ABZ*I_ESP_D2z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_P_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_D_D
     * RHS shell quartet name: SQ_ESP_P_D
     ************************************************************/
    Double I_ESP_Px_F3x = I_ESP_D2x_D2x+ABX*I_ESP_Px_D2x;
    Double I_ESP_Py_F3x = I_ESP_Dxy_D2x+ABX*I_ESP_Py_D2x;
    Double I_ESP_Pz_F3x = I_ESP_Dxz_D2x+ABX*I_ESP_Pz_D2x;
    Double I_ESP_Px_F2xy = I_ESP_Dxy_D2x+ABY*I_ESP_Px_D2x;
    Double I_ESP_Py_F2xy = I_ESP_D2y_D2x+ABY*I_ESP_Py_D2x;
    Double I_ESP_Pz_F2xy = I_ESP_Dyz_D2x+ABY*I_ESP_Pz_D2x;
    Double I_ESP_Px_F2xz = I_ESP_Dxz_D2x+ABZ*I_ESP_Px_D2x;
    Double I_ESP_Py_F2xz = I_ESP_Dyz_D2x+ABZ*I_ESP_Py_D2x;
    Double I_ESP_Pz_F2xz = I_ESP_D2z_D2x+ABZ*I_ESP_Pz_D2x;
    Double I_ESP_Px_Fx2y = I_ESP_D2x_D2y+ABX*I_ESP_Px_D2y;
    Double I_ESP_Py_Fx2y = I_ESP_Dxy_D2y+ABX*I_ESP_Py_D2y;
    Double I_ESP_Pz_Fx2y = I_ESP_Dxz_D2y+ABX*I_ESP_Pz_D2y;
    Double I_ESP_Px_Fxyz = I_ESP_Dxz_Dxy+ABZ*I_ESP_Px_Dxy;
    Double I_ESP_Py_Fxyz = I_ESP_Dyz_Dxy+ABZ*I_ESP_Py_Dxy;
    Double I_ESP_Pz_Fxyz = I_ESP_D2z_Dxy+ABZ*I_ESP_Pz_Dxy;
    Double I_ESP_Px_Fx2z = I_ESP_D2x_D2z+ABX*I_ESP_Px_D2z;
    Double I_ESP_Py_Fx2z = I_ESP_Dxy_D2z+ABX*I_ESP_Py_D2z;
    Double I_ESP_Pz_Fx2z = I_ESP_Dxz_D2z+ABX*I_ESP_Pz_D2z;
    Double I_ESP_Px_F3y = I_ESP_Dxy_D2y+ABY*I_ESP_Px_D2y;
    Double I_ESP_Py_F3y = I_ESP_D2y_D2y+ABY*I_ESP_Py_D2y;
    Double I_ESP_Pz_F3y = I_ESP_Dyz_D2y+ABY*I_ESP_Pz_D2y;
    Double I_ESP_Px_F2yz = I_ESP_Dxz_D2y+ABZ*I_ESP_Px_D2y;
    Double I_ESP_Py_F2yz = I_ESP_Dyz_D2y+ABZ*I_ESP_Py_D2y;
    Double I_ESP_Pz_F2yz = I_ESP_D2z_D2y+ABZ*I_ESP_Pz_D2y;
    Double I_ESP_Px_Fy2z = I_ESP_Dxy_D2z+ABY*I_ESP_Px_D2z;
    Double I_ESP_Py_Fy2z = I_ESP_D2y_D2z+ABY*I_ESP_Py_D2z;
    Double I_ESP_Pz_Fy2z = I_ESP_Dyz_D2z+ABY*I_ESP_Pz_D2z;
    Double I_ESP_Px_F3z = I_ESP_Dxz_D2z+ABZ*I_ESP_Px_D2z;
    Double I_ESP_Py_F3z = I_ESP_Dyz_D2z+ABZ*I_ESP_Py_D2z;
    Double I_ESP_Pz_F3z = I_ESP_D2z_D2z+ABZ*I_ESP_Pz_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     ************************************************************/
    Double I_ESP_F3x_Px_a = I_ESP_G4x_S_a+ABX*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Px_a = I_ESP_G3xy_S_a+ABX*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Px_a = I_ESP_G3xz_S_a+ABX*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Px_a = I_ESP_G2x2y_S_a+ABX*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Px_a = I_ESP_G2xyz_S_a+ABX*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Px_a = I_ESP_G2x2z_S_a+ABX*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Px_a = I_ESP_Gx3y_S_a+ABX*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Px_a = I_ESP_Gx2yz_S_a+ABX*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Px_a = I_ESP_Gxy2z_S_a+ABX*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Px_a = I_ESP_Gx3z_S_a+ABX*I_ESP_F3z_S_a;
    Double I_ESP_F3x_Py_a = I_ESP_G3xy_S_a+ABY*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Py_a = I_ESP_G2x2y_S_a+ABY*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Py_a = I_ESP_G2xyz_S_a+ABY*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Py_a = I_ESP_Gx3y_S_a+ABY*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Py_a = I_ESP_Gx2yz_S_a+ABY*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Py_a = I_ESP_Gxy2z_S_a+ABY*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Py_a = I_ESP_G4y_S_a+ABY*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Py_a = I_ESP_G3yz_S_a+ABY*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Py_a = I_ESP_G2y2z_S_a+ABY*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Py_a = I_ESP_Gy3z_S_a+ABY*I_ESP_F3z_S_a;
    Double I_ESP_F3x_Pz_a = I_ESP_G3xz_S_a+ABZ*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Pz_a = I_ESP_G2xyz_S_a+ABZ*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Pz_a = I_ESP_G2x2z_S_a+ABZ*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Pz_a = I_ESP_Gx2yz_S_a+ABZ*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Pz_a = I_ESP_Gxy2z_S_a+ABZ*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Pz_a = I_ESP_Gx3z_S_a+ABZ*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Pz_a = I_ESP_G3yz_S_a+ABZ*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Pz_a = I_ESP_G2y2z_S_a+ABZ*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Pz_a = I_ESP_Gy3z_S_a+ABZ*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Pz_a = I_ESP_G4z_S_a+ABZ*I_ESP_F3z_S_a;

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
     * shell quartet name: SQ_ESP_F_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     ************************************************************/
    Double I_ESP_F3x_D2x_a = I_ESP_G4x_Px_a+ABX*I_ESP_F3x_Px_a;
    Double I_ESP_F2xy_D2x_a = I_ESP_G3xy_Px_a+ABX*I_ESP_F2xy_Px_a;
    Double I_ESP_F2xz_D2x_a = I_ESP_G3xz_Px_a+ABX*I_ESP_F2xz_Px_a;
    Double I_ESP_Fx2y_D2x_a = I_ESP_G2x2y_Px_a+ABX*I_ESP_Fx2y_Px_a;
    Double I_ESP_Fxyz_D2x_a = I_ESP_G2xyz_Px_a+ABX*I_ESP_Fxyz_Px_a;
    Double I_ESP_Fx2z_D2x_a = I_ESP_G2x2z_Px_a+ABX*I_ESP_Fx2z_Px_a;
    Double I_ESP_F3y_D2x_a = I_ESP_Gx3y_Px_a+ABX*I_ESP_F3y_Px_a;
    Double I_ESP_F2yz_D2x_a = I_ESP_Gx2yz_Px_a+ABX*I_ESP_F2yz_Px_a;
    Double I_ESP_Fy2z_D2x_a = I_ESP_Gxy2z_Px_a+ABX*I_ESP_Fy2z_Px_a;
    Double I_ESP_F3z_D2x_a = I_ESP_Gx3z_Px_a+ABX*I_ESP_F3z_Px_a;
    Double I_ESP_F3x_Dxy_a = I_ESP_G3xy_Px_a+ABY*I_ESP_F3x_Px_a;
    Double I_ESP_F2xy_Dxy_a = I_ESP_G2x2y_Px_a+ABY*I_ESP_F2xy_Px_a;
    Double I_ESP_F2xz_Dxy_a = I_ESP_G2xyz_Px_a+ABY*I_ESP_F2xz_Px_a;
    Double I_ESP_Fx2y_Dxy_a = I_ESP_Gx3y_Px_a+ABY*I_ESP_Fx2y_Px_a;
    Double I_ESP_Fxyz_Dxy_a = I_ESP_Gx2yz_Px_a+ABY*I_ESP_Fxyz_Px_a;
    Double I_ESP_Fx2z_Dxy_a = I_ESP_Gxy2z_Px_a+ABY*I_ESP_Fx2z_Px_a;
    Double I_ESP_F3y_Dxy_a = I_ESP_G4y_Px_a+ABY*I_ESP_F3y_Px_a;
    Double I_ESP_F2yz_Dxy_a = I_ESP_G3yz_Px_a+ABY*I_ESP_F2yz_Px_a;
    Double I_ESP_Fy2z_Dxy_a = I_ESP_G2y2z_Px_a+ABY*I_ESP_Fy2z_Px_a;
    Double I_ESP_F3z_Dxy_a = I_ESP_Gy3z_Px_a+ABY*I_ESP_F3z_Px_a;
    Double I_ESP_F3x_D2y_a = I_ESP_G3xy_Py_a+ABY*I_ESP_F3x_Py_a;
    Double I_ESP_F2xy_D2y_a = I_ESP_G2x2y_Py_a+ABY*I_ESP_F2xy_Py_a;
    Double I_ESP_F2xz_D2y_a = I_ESP_G2xyz_Py_a+ABY*I_ESP_F2xz_Py_a;
    Double I_ESP_Fx2y_D2y_a = I_ESP_Gx3y_Py_a+ABY*I_ESP_Fx2y_Py_a;
    Double I_ESP_Fxyz_D2y_a = I_ESP_Gx2yz_Py_a+ABY*I_ESP_Fxyz_Py_a;
    Double I_ESP_Fx2z_D2y_a = I_ESP_Gxy2z_Py_a+ABY*I_ESP_Fx2z_Py_a;
    Double I_ESP_F3y_D2y_a = I_ESP_G4y_Py_a+ABY*I_ESP_F3y_Py_a;
    Double I_ESP_F2yz_D2y_a = I_ESP_G3yz_Py_a+ABY*I_ESP_F2yz_Py_a;
    Double I_ESP_Fy2z_D2y_a = I_ESP_G2y2z_Py_a+ABY*I_ESP_Fy2z_Py_a;
    Double I_ESP_F3z_D2y_a = I_ESP_Gy3z_Py_a+ABY*I_ESP_F3z_Py_a;
    Double I_ESP_F3x_D2z_a = I_ESP_G3xz_Pz_a+ABZ*I_ESP_F3x_Pz_a;
    Double I_ESP_F2xy_D2z_a = I_ESP_G2xyz_Pz_a+ABZ*I_ESP_F2xy_Pz_a;
    Double I_ESP_F2xz_D2z_a = I_ESP_G2x2z_Pz_a+ABZ*I_ESP_F2xz_Pz_a;
    Double I_ESP_Fx2y_D2z_a = I_ESP_Gx2yz_Pz_a+ABZ*I_ESP_Fx2y_Pz_a;
    Double I_ESP_Fxyz_D2z_a = I_ESP_Gxy2z_Pz_a+ABZ*I_ESP_Fxyz_Pz_a;
    Double I_ESP_Fx2z_D2z_a = I_ESP_Gx3z_Pz_a+ABZ*I_ESP_Fx2z_Pz_a;
    Double I_ESP_F3y_D2z_a = I_ESP_G3yz_Pz_a+ABZ*I_ESP_F3y_Pz_a;
    Double I_ESP_F2yz_D2z_a = I_ESP_G2y2z_Pz_a+ABZ*I_ESP_F2yz_Pz_a;
    Double I_ESP_Fy2z_D2z_a = I_ESP_Gy3z_Pz_a+ABZ*I_ESP_Fy2z_Pz_a;
    Double I_ESP_F3z_D2z_a = I_ESP_G4z_Pz_a+ABZ*I_ESP_F3z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 14 integrals are omitted 
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
    Double I_ESP_H4yz_Px_a = I_ESP_Ix4yz_S_a+ABX*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Px_a = I_ESP_Ix3y2z_S_a+ABX*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Px_a = I_ESP_Ix2y3z_S_a+ABX*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Px_a = I_ESP_Ixy4z_S_a+ABX*I_ESP_Hy4z_S_a;
    Double I_ESP_H4xy_Py_a = I_ESP_I4x2y_S_a+ABY*I_ESP_H4xy_S_a;
    Double I_ESP_H3x2y_Py_a = I_ESP_I3x3y_S_a+ABY*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Py_a = I_ESP_I3x2yz_S_a+ABY*I_ESP_H3xyz_S_a;
    Double I_ESP_H2x3y_Py_a = I_ESP_I2x4y_S_a+ABY*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Py_a = I_ESP_I2x3yz_S_a+ABY*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Py_a = I_ESP_I2x2y2z_S_a+ABY*I_ESP_H2xy2z_S_a;
    Double I_ESP_Hx4y_Py_a = I_ESP_Ix5y_S_a+ABY*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Py_a = I_ESP_Ix4yz_S_a+ABY*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Py_a = I_ESP_Ix3y2z_S_a+ABY*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Py_a = I_ESP_Ix2y3z_S_a+ABY*I_ESP_Hxy3z_S_a;
    Double I_ESP_H5y_Py_a = I_ESP_I6y_S_a+ABY*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Py_a = I_ESP_I5yz_S_a+ABY*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Py_a = I_ESP_I4y2z_S_a+ABY*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Py_a = I_ESP_I3y3z_S_a+ABY*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Py_a = I_ESP_I2y4z_S_a+ABY*I_ESP_Hy4z_S_a;
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
     * totally 35 integrals are omitted 
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
    Double I_ESP_G3xz_Dxy_a = I_ESP_H3xyz_Px_a+ABY*I_ESP_G3xz_Px_a;
    Double I_ESP_G2xyz_Dxy_a = I_ESP_H2x2yz_Px_a+ABY*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_Dxy_a = I_ESP_H2xy2z_Px_a+ABY*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx2yz_Dxy_a = I_ESP_Hx3yz_Px_a+ABY*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_Dxy_a = I_ESP_Hx2y2z_Px_a+ABY*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_Dxy_a = I_ESP_Hxy3z_Px_a+ABY*I_ESP_Gx3z_Px_a;
    Double I_ESP_G3yz_Dxy_a = I_ESP_H4yz_Px_a+ABY*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_Dxy_a = I_ESP_H3y2z_Px_a+ABY*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_Dxy_a = I_ESP_H2y3z_Px_a+ABY*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_Dxy_a = I_ESP_Hy4z_Px_a+ABY*I_ESP_G4z_Px_a;
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
     * shell quartet name: SQ_ESP_F_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_D_a
     * RHS shell quartet name: SQ_ESP_F_D_a
     ************************************************************/
    Double I_ESP_F3x_F3x_a = I_ESP_G4x_D2x_a+ABX*I_ESP_F3x_D2x_a;
    Double I_ESP_F2xy_F3x_a = I_ESP_G3xy_D2x_a+ABX*I_ESP_F2xy_D2x_a;
    Double I_ESP_F2xz_F3x_a = I_ESP_G3xz_D2x_a+ABX*I_ESP_F2xz_D2x_a;
    Double I_ESP_Fx2y_F3x_a = I_ESP_G2x2y_D2x_a+ABX*I_ESP_Fx2y_D2x_a;
    Double I_ESP_Fxyz_F3x_a = I_ESP_G2xyz_D2x_a+ABX*I_ESP_Fxyz_D2x_a;
    Double I_ESP_Fx2z_F3x_a = I_ESP_G2x2z_D2x_a+ABX*I_ESP_Fx2z_D2x_a;
    Double I_ESP_F3y_F3x_a = I_ESP_Gx3y_D2x_a+ABX*I_ESP_F3y_D2x_a;
    Double I_ESP_F2yz_F3x_a = I_ESP_Gx2yz_D2x_a+ABX*I_ESP_F2yz_D2x_a;
    Double I_ESP_Fy2z_F3x_a = I_ESP_Gxy2z_D2x_a+ABX*I_ESP_Fy2z_D2x_a;
    Double I_ESP_F3z_F3x_a = I_ESP_Gx3z_D2x_a+ABX*I_ESP_F3z_D2x_a;
    Double I_ESP_F3x_F2xy_a = I_ESP_G3xy_D2x_a+ABY*I_ESP_F3x_D2x_a;
    Double I_ESP_F2xy_F2xy_a = I_ESP_G2x2y_D2x_a+ABY*I_ESP_F2xy_D2x_a;
    Double I_ESP_F2xz_F2xy_a = I_ESP_G2xyz_D2x_a+ABY*I_ESP_F2xz_D2x_a;
    Double I_ESP_Fx2y_F2xy_a = I_ESP_Gx3y_D2x_a+ABY*I_ESP_Fx2y_D2x_a;
    Double I_ESP_Fxyz_F2xy_a = I_ESP_Gx2yz_D2x_a+ABY*I_ESP_Fxyz_D2x_a;
    Double I_ESP_Fx2z_F2xy_a = I_ESP_Gxy2z_D2x_a+ABY*I_ESP_Fx2z_D2x_a;
    Double I_ESP_F3y_F2xy_a = I_ESP_G4y_D2x_a+ABY*I_ESP_F3y_D2x_a;
    Double I_ESP_F2yz_F2xy_a = I_ESP_G3yz_D2x_a+ABY*I_ESP_F2yz_D2x_a;
    Double I_ESP_Fy2z_F2xy_a = I_ESP_G2y2z_D2x_a+ABY*I_ESP_Fy2z_D2x_a;
    Double I_ESP_F3z_F2xy_a = I_ESP_Gy3z_D2x_a+ABY*I_ESP_F3z_D2x_a;
    Double I_ESP_F3x_F2xz_a = I_ESP_G3xz_D2x_a+ABZ*I_ESP_F3x_D2x_a;
    Double I_ESP_F2xy_F2xz_a = I_ESP_G2xyz_D2x_a+ABZ*I_ESP_F2xy_D2x_a;
    Double I_ESP_F2xz_F2xz_a = I_ESP_G2x2z_D2x_a+ABZ*I_ESP_F2xz_D2x_a;
    Double I_ESP_Fx2y_F2xz_a = I_ESP_Gx2yz_D2x_a+ABZ*I_ESP_Fx2y_D2x_a;
    Double I_ESP_Fxyz_F2xz_a = I_ESP_Gxy2z_D2x_a+ABZ*I_ESP_Fxyz_D2x_a;
    Double I_ESP_Fx2z_F2xz_a = I_ESP_Gx3z_D2x_a+ABZ*I_ESP_Fx2z_D2x_a;
    Double I_ESP_F3y_F2xz_a = I_ESP_G3yz_D2x_a+ABZ*I_ESP_F3y_D2x_a;
    Double I_ESP_F2yz_F2xz_a = I_ESP_G2y2z_D2x_a+ABZ*I_ESP_F2yz_D2x_a;
    Double I_ESP_Fy2z_F2xz_a = I_ESP_Gy3z_D2x_a+ABZ*I_ESP_Fy2z_D2x_a;
    Double I_ESP_F3z_F2xz_a = I_ESP_G4z_D2x_a+ABZ*I_ESP_F3z_D2x_a;
    Double I_ESP_F3x_Fx2y_a = I_ESP_G4x_D2y_a+ABX*I_ESP_F3x_D2y_a;
    Double I_ESP_F2xy_Fx2y_a = I_ESP_G3xy_D2y_a+ABX*I_ESP_F2xy_D2y_a;
    Double I_ESP_F2xz_Fx2y_a = I_ESP_G3xz_D2y_a+ABX*I_ESP_F2xz_D2y_a;
    Double I_ESP_Fx2y_Fx2y_a = I_ESP_G2x2y_D2y_a+ABX*I_ESP_Fx2y_D2y_a;
    Double I_ESP_Fxyz_Fx2y_a = I_ESP_G2xyz_D2y_a+ABX*I_ESP_Fxyz_D2y_a;
    Double I_ESP_Fx2z_Fx2y_a = I_ESP_G2x2z_D2y_a+ABX*I_ESP_Fx2z_D2y_a;
    Double I_ESP_F3y_Fx2y_a = I_ESP_Gx3y_D2y_a+ABX*I_ESP_F3y_D2y_a;
    Double I_ESP_F2yz_Fx2y_a = I_ESP_Gx2yz_D2y_a+ABX*I_ESP_F2yz_D2y_a;
    Double I_ESP_Fy2z_Fx2y_a = I_ESP_Gxy2z_D2y_a+ABX*I_ESP_Fy2z_D2y_a;
    Double I_ESP_F3z_Fx2y_a = I_ESP_Gx3z_D2y_a+ABX*I_ESP_F3z_D2y_a;
    Double I_ESP_F3x_Fxyz_a = I_ESP_G3xz_Dxy_a+ABZ*I_ESP_F3x_Dxy_a;
    Double I_ESP_F2xy_Fxyz_a = I_ESP_G2xyz_Dxy_a+ABZ*I_ESP_F2xy_Dxy_a;
    Double I_ESP_F2xz_Fxyz_a = I_ESP_G2x2z_Dxy_a+ABZ*I_ESP_F2xz_Dxy_a;
    Double I_ESP_Fx2y_Fxyz_a = I_ESP_Gx2yz_Dxy_a+ABZ*I_ESP_Fx2y_Dxy_a;
    Double I_ESP_Fxyz_Fxyz_a = I_ESP_Gxy2z_Dxy_a+ABZ*I_ESP_Fxyz_Dxy_a;
    Double I_ESP_Fx2z_Fxyz_a = I_ESP_Gx3z_Dxy_a+ABZ*I_ESP_Fx2z_Dxy_a;
    Double I_ESP_F3y_Fxyz_a = I_ESP_G3yz_Dxy_a+ABZ*I_ESP_F3y_Dxy_a;
    Double I_ESP_F2yz_Fxyz_a = I_ESP_G2y2z_Dxy_a+ABZ*I_ESP_F2yz_Dxy_a;
    Double I_ESP_Fy2z_Fxyz_a = I_ESP_Gy3z_Dxy_a+ABZ*I_ESP_Fy2z_Dxy_a;
    Double I_ESP_F3z_Fxyz_a = I_ESP_G4z_Dxy_a+ABZ*I_ESP_F3z_Dxy_a;
    Double I_ESP_F3x_Fx2z_a = I_ESP_G4x_D2z_a+ABX*I_ESP_F3x_D2z_a;
    Double I_ESP_F2xy_Fx2z_a = I_ESP_G3xy_D2z_a+ABX*I_ESP_F2xy_D2z_a;
    Double I_ESP_F2xz_Fx2z_a = I_ESP_G3xz_D2z_a+ABX*I_ESP_F2xz_D2z_a;
    Double I_ESP_Fx2y_Fx2z_a = I_ESP_G2x2y_D2z_a+ABX*I_ESP_Fx2y_D2z_a;
    Double I_ESP_Fxyz_Fx2z_a = I_ESP_G2xyz_D2z_a+ABX*I_ESP_Fxyz_D2z_a;
    Double I_ESP_Fx2z_Fx2z_a = I_ESP_G2x2z_D2z_a+ABX*I_ESP_Fx2z_D2z_a;
    Double I_ESP_F3y_Fx2z_a = I_ESP_Gx3y_D2z_a+ABX*I_ESP_F3y_D2z_a;
    Double I_ESP_F2yz_Fx2z_a = I_ESP_Gx2yz_D2z_a+ABX*I_ESP_F2yz_D2z_a;
    Double I_ESP_Fy2z_Fx2z_a = I_ESP_Gxy2z_D2z_a+ABX*I_ESP_Fy2z_D2z_a;
    Double I_ESP_F3z_Fx2z_a = I_ESP_Gx3z_D2z_a+ABX*I_ESP_F3z_D2z_a;
    Double I_ESP_F3x_F3y_a = I_ESP_G3xy_D2y_a+ABY*I_ESP_F3x_D2y_a;
    Double I_ESP_F2xy_F3y_a = I_ESP_G2x2y_D2y_a+ABY*I_ESP_F2xy_D2y_a;
    Double I_ESP_F2xz_F3y_a = I_ESP_G2xyz_D2y_a+ABY*I_ESP_F2xz_D2y_a;
    Double I_ESP_Fx2y_F3y_a = I_ESP_Gx3y_D2y_a+ABY*I_ESP_Fx2y_D2y_a;
    Double I_ESP_Fxyz_F3y_a = I_ESP_Gx2yz_D2y_a+ABY*I_ESP_Fxyz_D2y_a;
    Double I_ESP_Fx2z_F3y_a = I_ESP_Gxy2z_D2y_a+ABY*I_ESP_Fx2z_D2y_a;
    Double I_ESP_F3y_F3y_a = I_ESP_G4y_D2y_a+ABY*I_ESP_F3y_D2y_a;
    Double I_ESP_F2yz_F3y_a = I_ESP_G3yz_D2y_a+ABY*I_ESP_F2yz_D2y_a;
    Double I_ESP_Fy2z_F3y_a = I_ESP_G2y2z_D2y_a+ABY*I_ESP_Fy2z_D2y_a;
    Double I_ESP_F3z_F3y_a = I_ESP_Gy3z_D2y_a+ABY*I_ESP_F3z_D2y_a;
    Double I_ESP_F3x_F2yz_a = I_ESP_G3xz_D2y_a+ABZ*I_ESP_F3x_D2y_a;
    Double I_ESP_F2xy_F2yz_a = I_ESP_G2xyz_D2y_a+ABZ*I_ESP_F2xy_D2y_a;
    Double I_ESP_F2xz_F2yz_a = I_ESP_G2x2z_D2y_a+ABZ*I_ESP_F2xz_D2y_a;
    Double I_ESP_Fx2y_F2yz_a = I_ESP_Gx2yz_D2y_a+ABZ*I_ESP_Fx2y_D2y_a;
    Double I_ESP_Fxyz_F2yz_a = I_ESP_Gxy2z_D2y_a+ABZ*I_ESP_Fxyz_D2y_a;
    Double I_ESP_Fx2z_F2yz_a = I_ESP_Gx3z_D2y_a+ABZ*I_ESP_Fx2z_D2y_a;
    Double I_ESP_F3y_F2yz_a = I_ESP_G3yz_D2y_a+ABZ*I_ESP_F3y_D2y_a;
    Double I_ESP_F2yz_F2yz_a = I_ESP_G2y2z_D2y_a+ABZ*I_ESP_F2yz_D2y_a;
    Double I_ESP_Fy2z_F2yz_a = I_ESP_Gy3z_D2y_a+ABZ*I_ESP_Fy2z_D2y_a;
    Double I_ESP_F3z_F2yz_a = I_ESP_G4z_D2y_a+ABZ*I_ESP_F3z_D2y_a;
    Double I_ESP_F3x_Fy2z_a = I_ESP_G3xy_D2z_a+ABY*I_ESP_F3x_D2z_a;
    Double I_ESP_F2xy_Fy2z_a = I_ESP_G2x2y_D2z_a+ABY*I_ESP_F2xy_D2z_a;
    Double I_ESP_F2xz_Fy2z_a = I_ESP_G2xyz_D2z_a+ABY*I_ESP_F2xz_D2z_a;
    Double I_ESP_Fx2y_Fy2z_a = I_ESP_Gx3y_D2z_a+ABY*I_ESP_Fx2y_D2z_a;
    Double I_ESP_Fxyz_Fy2z_a = I_ESP_Gx2yz_D2z_a+ABY*I_ESP_Fxyz_D2z_a;
    Double I_ESP_Fx2z_Fy2z_a = I_ESP_Gxy2z_D2z_a+ABY*I_ESP_Fx2z_D2z_a;
    Double I_ESP_F3y_Fy2z_a = I_ESP_G4y_D2z_a+ABY*I_ESP_F3y_D2z_a;
    Double I_ESP_F2yz_Fy2z_a = I_ESP_G3yz_D2z_a+ABY*I_ESP_F2yz_D2z_a;
    Double I_ESP_Fy2z_Fy2z_a = I_ESP_G2y2z_D2z_a+ABY*I_ESP_Fy2z_D2z_a;
    Double I_ESP_F3z_Fy2z_a = I_ESP_Gy3z_D2z_a+ABY*I_ESP_F3z_D2z_a;
    Double I_ESP_F3x_F3z_a = I_ESP_G3xz_D2z_a+ABZ*I_ESP_F3x_D2z_a;
    Double I_ESP_F2xy_F3z_a = I_ESP_G2xyz_D2z_a+ABZ*I_ESP_F2xy_D2z_a;
    Double I_ESP_F2xz_F3z_a = I_ESP_G2x2z_D2z_a+ABZ*I_ESP_F2xz_D2z_a;
    Double I_ESP_Fx2y_F3z_a = I_ESP_Gx2yz_D2z_a+ABZ*I_ESP_Fx2y_D2z_a;
    Double I_ESP_Fxyz_F3z_a = I_ESP_Gxy2z_D2z_a+ABZ*I_ESP_Fxyz_D2z_a;
    Double I_ESP_Fx2z_F3z_a = I_ESP_Gx3z_D2z_a+ABZ*I_ESP_Fx2z_D2z_a;
    Double I_ESP_F3y_F3z_a = I_ESP_G3yz_D2z_a+ABZ*I_ESP_F3y_D2z_a;
    Double I_ESP_F2yz_F3z_a = I_ESP_G2y2z_D2z_a+ABZ*I_ESP_F2yz_D2z_a;
    Double I_ESP_Fy2z_F3z_a = I_ESP_Gy3z_D2z_a+ABZ*I_ESP_Fy2z_D2z_a;
    Double I_ESP_F3z_F3z_a = I_ESP_G4z_D2z_a+ABZ*I_ESP_F3z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_aa
     * RHS shell quartet name: SQ_ESP_H_S_aa
     ************************************************************/
    Double I_ESP_H5x_Px_aa = I_ESP_I6x_S_aa+ABX*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Px_aa = I_ESP_I5xy_S_aa+ABX*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Px_aa = I_ESP_I5xz_S_aa+ABX*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Px_aa = I_ESP_I4x2y_S_aa+ABX*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Px_aa = I_ESP_I4xyz_S_aa+ABX*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Px_aa = I_ESP_I4x2z_S_aa+ABX*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Px_aa = I_ESP_I3x3y_S_aa+ABX*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Px_aa = I_ESP_I3x2yz_S_aa+ABX*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Px_aa = I_ESP_I3xy2z_S_aa+ABX*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Px_aa = I_ESP_I3x3z_S_aa+ABX*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Px_aa = I_ESP_I2x4y_S_aa+ABX*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Px_aa = I_ESP_I2x3yz_S_aa+ABX*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Px_aa = I_ESP_I2x2y2z_S_aa+ABX*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Px_aa = I_ESP_I2xy3z_S_aa+ABX*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Px_aa = I_ESP_I2x4z_S_aa+ABX*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Px_aa = I_ESP_Ix5y_S_aa+ABX*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Px_aa = I_ESP_Ix4yz_S_aa+ABX*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Px_aa = I_ESP_Ix3y2z_S_aa+ABX*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Px_aa = I_ESP_Ix2y3z_S_aa+ABX*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Px_aa = I_ESP_Ixy4z_S_aa+ABX*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Px_aa = I_ESP_Ix5z_S_aa+ABX*I_ESP_H5z_S_aa;
    Double I_ESP_H5x_Py_aa = I_ESP_I5xy_S_aa+ABY*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Py_aa = I_ESP_I4x2y_S_aa+ABY*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Py_aa = I_ESP_I4xyz_S_aa+ABY*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Py_aa = I_ESP_I3x3y_S_aa+ABY*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Py_aa = I_ESP_I3x2yz_S_aa+ABY*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Py_aa = I_ESP_I3xy2z_S_aa+ABY*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Py_aa = I_ESP_I2x4y_S_aa+ABY*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Py_aa = I_ESP_I2x3yz_S_aa+ABY*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Py_aa = I_ESP_I2x2y2z_S_aa+ABY*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Py_aa = I_ESP_I2xy3z_S_aa+ABY*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Py_aa = I_ESP_Ix5y_S_aa+ABY*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Py_aa = I_ESP_Ix4yz_S_aa+ABY*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Py_aa = I_ESP_Ix3y2z_S_aa+ABY*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Py_aa = I_ESP_Ix2y3z_S_aa+ABY*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Py_aa = I_ESP_Ixy4z_S_aa+ABY*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Py_aa = I_ESP_I6y_S_aa+ABY*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Py_aa = I_ESP_I5yz_S_aa+ABY*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Py_aa = I_ESP_I4y2z_S_aa+ABY*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Py_aa = I_ESP_I3y3z_S_aa+ABY*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Py_aa = I_ESP_I2y4z_S_aa+ABY*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Py_aa = I_ESP_Iy5z_S_aa+ABY*I_ESP_H5z_S_aa;
    Double I_ESP_H5x_Pz_aa = I_ESP_I5xz_S_aa+ABZ*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Pz_aa = I_ESP_I4xyz_S_aa+ABZ*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Pz_aa = I_ESP_I4x2z_S_aa+ABZ*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Pz_aa = I_ESP_I3x2yz_S_aa+ABZ*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Pz_aa = I_ESP_I3xy2z_S_aa+ABZ*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Pz_aa = I_ESP_I3x3z_S_aa+ABZ*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Pz_aa = I_ESP_I2x3yz_S_aa+ABZ*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Pz_aa = I_ESP_I2x2y2z_S_aa+ABZ*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Pz_aa = I_ESP_I2xy3z_S_aa+ABZ*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Pz_aa = I_ESP_I2x4z_S_aa+ABZ*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Pz_aa = I_ESP_Ix4yz_S_aa+ABZ*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Pz_aa = I_ESP_Ix3y2z_S_aa+ABZ*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Pz_aa = I_ESP_Ix2y3z_S_aa+ABZ*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Pz_aa = I_ESP_Ixy4z_S_aa+ABZ*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Pz_aa = I_ESP_Ix5z_S_aa+ABZ*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Pz_aa = I_ESP_I5yz_S_aa+ABZ*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Pz_aa = I_ESP_I4y2z_S_aa+ABZ*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Pz_aa = I_ESP_I3y3z_S_aa+ABZ*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Pz_aa = I_ESP_I2y4z_S_aa+ABZ*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Pz_aa = I_ESP_Iy5z_S_aa+ABZ*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Pz_aa = I_ESP_I6z_S_aa+ABZ*I_ESP_H5z_S_aa;

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
     * shell quartet name: SQ_ESP_H_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_aa
     * RHS shell quartet name: SQ_ESP_H_P_aa
     ************************************************************/
    Double I_ESP_H5x_D2x_aa = I_ESP_I6x_Px_aa+ABX*I_ESP_H5x_Px_aa;
    Double I_ESP_H4xy_D2x_aa = I_ESP_I5xy_Px_aa+ABX*I_ESP_H4xy_Px_aa;
    Double I_ESP_H4xz_D2x_aa = I_ESP_I5xz_Px_aa+ABX*I_ESP_H4xz_Px_aa;
    Double I_ESP_H3x2y_D2x_aa = I_ESP_I4x2y_Px_aa+ABX*I_ESP_H3x2y_Px_aa;
    Double I_ESP_H3xyz_D2x_aa = I_ESP_I4xyz_Px_aa+ABX*I_ESP_H3xyz_Px_aa;
    Double I_ESP_H3x2z_D2x_aa = I_ESP_I4x2z_Px_aa+ABX*I_ESP_H3x2z_Px_aa;
    Double I_ESP_H2x3y_D2x_aa = I_ESP_I3x3y_Px_aa+ABX*I_ESP_H2x3y_Px_aa;
    Double I_ESP_H2x2yz_D2x_aa = I_ESP_I3x2yz_Px_aa+ABX*I_ESP_H2x2yz_Px_aa;
    Double I_ESP_H2xy2z_D2x_aa = I_ESP_I3xy2z_Px_aa+ABX*I_ESP_H2xy2z_Px_aa;
    Double I_ESP_H2x3z_D2x_aa = I_ESP_I3x3z_Px_aa+ABX*I_ESP_H2x3z_Px_aa;
    Double I_ESP_Hx4y_D2x_aa = I_ESP_I2x4y_Px_aa+ABX*I_ESP_Hx4y_Px_aa;
    Double I_ESP_Hx3yz_D2x_aa = I_ESP_I2x3yz_Px_aa+ABX*I_ESP_Hx3yz_Px_aa;
    Double I_ESP_Hx2y2z_D2x_aa = I_ESP_I2x2y2z_Px_aa+ABX*I_ESP_Hx2y2z_Px_aa;
    Double I_ESP_Hxy3z_D2x_aa = I_ESP_I2xy3z_Px_aa+ABX*I_ESP_Hxy3z_Px_aa;
    Double I_ESP_Hx4z_D2x_aa = I_ESP_I2x4z_Px_aa+ABX*I_ESP_Hx4z_Px_aa;
    Double I_ESP_H5y_D2x_aa = I_ESP_Ix5y_Px_aa+ABX*I_ESP_H5y_Px_aa;
    Double I_ESP_H4yz_D2x_aa = I_ESP_Ix4yz_Px_aa+ABX*I_ESP_H4yz_Px_aa;
    Double I_ESP_H3y2z_D2x_aa = I_ESP_Ix3y2z_Px_aa+ABX*I_ESP_H3y2z_Px_aa;
    Double I_ESP_H2y3z_D2x_aa = I_ESP_Ix2y3z_Px_aa+ABX*I_ESP_H2y3z_Px_aa;
    Double I_ESP_Hy4z_D2x_aa = I_ESP_Ixy4z_Px_aa+ABX*I_ESP_Hy4z_Px_aa;
    Double I_ESP_H5z_D2x_aa = I_ESP_Ix5z_Px_aa+ABX*I_ESP_H5z_Px_aa;
    Double I_ESP_H5x_Dxy_aa = I_ESP_I5xy_Px_aa+ABY*I_ESP_H5x_Px_aa;
    Double I_ESP_H4xy_Dxy_aa = I_ESP_I4x2y_Px_aa+ABY*I_ESP_H4xy_Px_aa;
    Double I_ESP_H4xz_Dxy_aa = I_ESP_I4xyz_Px_aa+ABY*I_ESP_H4xz_Px_aa;
    Double I_ESP_H3x2y_Dxy_aa = I_ESP_I3x3y_Px_aa+ABY*I_ESP_H3x2y_Px_aa;
    Double I_ESP_H3xyz_Dxy_aa = I_ESP_I3x2yz_Px_aa+ABY*I_ESP_H3xyz_Px_aa;
    Double I_ESP_H3x2z_Dxy_aa = I_ESP_I3xy2z_Px_aa+ABY*I_ESP_H3x2z_Px_aa;
    Double I_ESP_H2x3y_Dxy_aa = I_ESP_I2x4y_Px_aa+ABY*I_ESP_H2x3y_Px_aa;
    Double I_ESP_H2x2yz_Dxy_aa = I_ESP_I2x3yz_Px_aa+ABY*I_ESP_H2x2yz_Px_aa;
    Double I_ESP_H2xy2z_Dxy_aa = I_ESP_I2x2y2z_Px_aa+ABY*I_ESP_H2xy2z_Px_aa;
    Double I_ESP_H2x3z_Dxy_aa = I_ESP_I2xy3z_Px_aa+ABY*I_ESP_H2x3z_Px_aa;
    Double I_ESP_Hx4y_Dxy_aa = I_ESP_Ix5y_Px_aa+ABY*I_ESP_Hx4y_Px_aa;
    Double I_ESP_Hx3yz_Dxy_aa = I_ESP_Ix4yz_Px_aa+ABY*I_ESP_Hx3yz_Px_aa;
    Double I_ESP_Hx2y2z_Dxy_aa = I_ESP_Ix3y2z_Px_aa+ABY*I_ESP_Hx2y2z_Px_aa;
    Double I_ESP_Hxy3z_Dxy_aa = I_ESP_Ix2y3z_Px_aa+ABY*I_ESP_Hxy3z_Px_aa;
    Double I_ESP_Hx4z_Dxy_aa = I_ESP_Ixy4z_Px_aa+ABY*I_ESP_Hx4z_Px_aa;
    Double I_ESP_H5y_Dxy_aa = I_ESP_I6y_Px_aa+ABY*I_ESP_H5y_Px_aa;
    Double I_ESP_H4yz_Dxy_aa = I_ESP_I5yz_Px_aa+ABY*I_ESP_H4yz_Px_aa;
    Double I_ESP_H3y2z_Dxy_aa = I_ESP_I4y2z_Px_aa+ABY*I_ESP_H3y2z_Px_aa;
    Double I_ESP_H2y3z_Dxy_aa = I_ESP_I3y3z_Px_aa+ABY*I_ESP_H2y3z_Px_aa;
    Double I_ESP_Hy4z_Dxy_aa = I_ESP_I2y4z_Px_aa+ABY*I_ESP_Hy4z_Px_aa;
    Double I_ESP_H5z_Dxy_aa = I_ESP_Iy5z_Px_aa+ABY*I_ESP_H5z_Px_aa;
    Double I_ESP_H5x_D2y_aa = I_ESP_I5xy_Py_aa+ABY*I_ESP_H5x_Py_aa;
    Double I_ESP_H4xy_D2y_aa = I_ESP_I4x2y_Py_aa+ABY*I_ESP_H4xy_Py_aa;
    Double I_ESP_H4xz_D2y_aa = I_ESP_I4xyz_Py_aa+ABY*I_ESP_H4xz_Py_aa;
    Double I_ESP_H3x2y_D2y_aa = I_ESP_I3x3y_Py_aa+ABY*I_ESP_H3x2y_Py_aa;
    Double I_ESP_H3xyz_D2y_aa = I_ESP_I3x2yz_Py_aa+ABY*I_ESP_H3xyz_Py_aa;
    Double I_ESP_H3x2z_D2y_aa = I_ESP_I3xy2z_Py_aa+ABY*I_ESP_H3x2z_Py_aa;
    Double I_ESP_H2x3y_D2y_aa = I_ESP_I2x4y_Py_aa+ABY*I_ESP_H2x3y_Py_aa;
    Double I_ESP_H2x2yz_D2y_aa = I_ESP_I2x3yz_Py_aa+ABY*I_ESP_H2x2yz_Py_aa;
    Double I_ESP_H2xy2z_D2y_aa = I_ESP_I2x2y2z_Py_aa+ABY*I_ESP_H2xy2z_Py_aa;
    Double I_ESP_H2x3z_D2y_aa = I_ESP_I2xy3z_Py_aa+ABY*I_ESP_H2x3z_Py_aa;
    Double I_ESP_Hx4y_D2y_aa = I_ESP_Ix5y_Py_aa+ABY*I_ESP_Hx4y_Py_aa;
    Double I_ESP_Hx3yz_D2y_aa = I_ESP_Ix4yz_Py_aa+ABY*I_ESP_Hx3yz_Py_aa;
    Double I_ESP_Hx2y2z_D2y_aa = I_ESP_Ix3y2z_Py_aa+ABY*I_ESP_Hx2y2z_Py_aa;
    Double I_ESP_Hxy3z_D2y_aa = I_ESP_Ix2y3z_Py_aa+ABY*I_ESP_Hxy3z_Py_aa;
    Double I_ESP_Hx4z_D2y_aa = I_ESP_Ixy4z_Py_aa+ABY*I_ESP_Hx4z_Py_aa;
    Double I_ESP_H5y_D2y_aa = I_ESP_I6y_Py_aa+ABY*I_ESP_H5y_Py_aa;
    Double I_ESP_H4yz_D2y_aa = I_ESP_I5yz_Py_aa+ABY*I_ESP_H4yz_Py_aa;
    Double I_ESP_H3y2z_D2y_aa = I_ESP_I4y2z_Py_aa+ABY*I_ESP_H3y2z_Py_aa;
    Double I_ESP_H2y3z_D2y_aa = I_ESP_I3y3z_Py_aa+ABY*I_ESP_H2y3z_Py_aa;
    Double I_ESP_Hy4z_D2y_aa = I_ESP_I2y4z_Py_aa+ABY*I_ESP_Hy4z_Py_aa;
    Double I_ESP_H5z_D2y_aa = I_ESP_Iy5z_Py_aa+ABY*I_ESP_H5z_Py_aa;
    Double I_ESP_H5x_D2z_aa = I_ESP_I5xz_Pz_aa+ABZ*I_ESP_H5x_Pz_aa;
    Double I_ESP_H4xy_D2z_aa = I_ESP_I4xyz_Pz_aa+ABZ*I_ESP_H4xy_Pz_aa;
    Double I_ESP_H4xz_D2z_aa = I_ESP_I4x2z_Pz_aa+ABZ*I_ESP_H4xz_Pz_aa;
    Double I_ESP_H3x2y_D2z_aa = I_ESP_I3x2yz_Pz_aa+ABZ*I_ESP_H3x2y_Pz_aa;
    Double I_ESP_H3xyz_D2z_aa = I_ESP_I3xy2z_Pz_aa+ABZ*I_ESP_H3xyz_Pz_aa;
    Double I_ESP_H3x2z_D2z_aa = I_ESP_I3x3z_Pz_aa+ABZ*I_ESP_H3x2z_Pz_aa;
    Double I_ESP_H2x3y_D2z_aa = I_ESP_I2x3yz_Pz_aa+ABZ*I_ESP_H2x3y_Pz_aa;
    Double I_ESP_H2x2yz_D2z_aa = I_ESP_I2x2y2z_Pz_aa+ABZ*I_ESP_H2x2yz_Pz_aa;
    Double I_ESP_H2xy2z_D2z_aa = I_ESP_I2xy3z_Pz_aa+ABZ*I_ESP_H2xy2z_Pz_aa;
    Double I_ESP_H2x3z_D2z_aa = I_ESP_I2x4z_Pz_aa+ABZ*I_ESP_H2x3z_Pz_aa;
    Double I_ESP_Hx4y_D2z_aa = I_ESP_Ix4yz_Pz_aa+ABZ*I_ESP_Hx4y_Pz_aa;
    Double I_ESP_Hx3yz_D2z_aa = I_ESP_Ix3y2z_Pz_aa+ABZ*I_ESP_Hx3yz_Pz_aa;
    Double I_ESP_Hx2y2z_D2z_aa = I_ESP_Ix2y3z_Pz_aa+ABZ*I_ESP_Hx2y2z_Pz_aa;
    Double I_ESP_Hxy3z_D2z_aa = I_ESP_Ixy4z_Pz_aa+ABZ*I_ESP_Hxy3z_Pz_aa;
    Double I_ESP_Hx4z_D2z_aa = I_ESP_Ix5z_Pz_aa+ABZ*I_ESP_Hx4z_Pz_aa;
    Double I_ESP_H5y_D2z_aa = I_ESP_I5yz_Pz_aa+ABZ*I_ESP_H5y_Pz_aa;
    Double I_ESP_H4yz_D2z_aa = I_ESP_I4y2z_Pz_aa+ABZ*I_ESP_H4yz_Pz_aa;
    Double I_ESP_H3y2z_D2z_aa = I_ESP_I3y3z_Pz_aa+ABZ*I_ESP_H3y2z_Pz_aa;
    Double I_ESP_H2y3z_D2z_aa = I_ESP_I2y4z_Pz_aa+ABZ*I_ESP_H2y3z_Pz_aa;
    Double I_ESP_Hy4z_D2z_aa = I_ESP_Iy5z_Pz_aa+ABZ*I_ESP_Hy4z_Pz_aa;
    Double I_ESP_H5z_D2z_aa = I_ESP_I6z_Pz_aa+ABZ*I_ESP_H5z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 18 integrals are omitted 
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
    Double I_ESP_K6yz_Px_aa = I_ESP_Lx6yz_S_aa+ABX*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Px_aa = I_ESP_Lx5y2z_S_aa+ABX*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Px_aa = I_ESP_Lx4y3z_S_aa+ABX*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Px_aa = I_ESP_Lx3y4z_S_aa+ABX*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Px_aa = I_ESP_Lx2y5z_S_aa+ABX*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Px_aa = I_ESP_Lxy6z_S_aa+ABX*I_ESP_Ky6z_S_aa;
    Double I_ESP_K6xy_Py_aa = I_ESP_L6x2y_S_aa+ABY*I_ESP_K6xy_S_aa;
    Double I_ESP_K5x2y_Py_aa = I_ESP_L5x3y_S_aa+ABY*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Py_aa = I_ESP_L5x2yz_S_aa+ABY*I_ESP_K5xyz_S_aa;
    Double I_ESP_K4x3y_Py_aa = I_ESP_L4x4y_S_aa+ABY*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Py_aa = I_ESP_L4x3yz_S_aa+ABY*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Py_aa = I_ESP_L4x2y2z_S_aa+ABY*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K3x4y_Py_aa = I_ESP_L3x5y_S_aa+ABY*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Py_aa = I_ESP_L3x4yz_S_aa+ABY*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Py_aa = I_ESP_L3x3y2z_S_aa+ABY*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Py_aa = I_ESP_L3x2y3z_S_aa+ABY*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K2x5y_Py_aa = I_ESP_L2x6y_S_aa+ABY*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Py_aa = I_ESP_L2x5yz_S_aa+ABY*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Py_aa = I_ESP_L2x4y2z_S_aa+ABY*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Py_aa = I_ESP_L2x3y3z_S_aa+ABY*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Py_aa = I_ESP_L2x2y4z_S_aa+ABY*I_ESP_K2xy4z_S_aa;
    Double I_ESP_Kx6y_Py_aa = I_ESP_Lx7y_S_aa+ABY*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Py_aa = I_ESP_Lx6yz_S_aa+ABY*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Py_aa = I_ESP_Lx5y2z_S_aa+ABY*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Py_aa = I_ESP_Lx4y3z_S_aa+ABY*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Py_aa = I_ESP_Lx3y4z_S_aa+ABY*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Py_aa = I_ESP_Lx2y5z_S_aa+ABY*I_ESP_Kxy5z_S_aa;
    Double I_ESP_K7y_Py_aa = I_ESP_L8y_S_aa+ABY*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Py_aa = I_ESP_L7yz_S_aa+ABY*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Py_aa = I_ESP_L6y2z_S_aa+ABY*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Py_aa = I_ESP_L5y3z_S_aa+ABY*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Py_aa = I_ESP_L4y4z_S_aa+ABY*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Py_aa = I_ESP_L3y5z_S_aa+ABY*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Py_aa = I_ESP_L2y6z_S_aa+ABY*I_ESP_Ky6z_S_aa;
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
     * totally 63 integrals are omitted 
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
    Double I_ESP_I5xz_Dxy_aa = I_ESP_K5xyz_Px_aa+ABY*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4xyz_Dxy_aa = I_ESP_K4x2yz_Px_aa+ABY*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_Dxy_aa = I_ESP_K4xy2z_Px_aa+ABY*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x2yz_Dxy_aa = I_ESP_K3x3yz_Px_aa+ABY*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_Dxy_aa = I_ESP_K3x2y2z_Px_aa+ABY*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_Dxy_aa = I_ESP_K3xy3z_Px_aa+ABY*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x3yz_Dxy_aa = I_ESP_K2x4yz_Px_aa+ABY*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_Dxy_aa = I_ESP_K2x3y2z_Px_aa+ABY*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_Dxy_aa = I_ESP_K2x2y3z_Px_aa+ABY*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_Dxy_aa = I_ESP_K2xy4z_Px_aa+ABY*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix4yz_Dxy_aa = I_ESP_Kx5yz_Px_aa+ABY*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_Dxy_aa = I_ESP_Kx4y2z_Px_aa+ABY*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_Dxy_aa = I_ESP_Kx3y3z_Px_aa+ABY*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_Dxy_aa = I_ESP_Kx2y4z_Px_aa+ABY*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_Dxy_aa = I_ESP_Kxy5z_Px_aa+ABY*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I5yz_Dxy_aa = I_ESP_K6yz_Px_aa+ABY*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_Dxy_aa = I_ESP_K5y2z_Px_aa+ABY*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_Dxy_aa = I_ESP_K4y3z_Px_aa+ABY*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_Dxy_aa = I_ESP_K3y4z_Px_aa+ABY*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_Dxy_aa = I_ESP_K2y5z_Px_aa+ABY*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_Dxy_aa = I_ESP_Ky6z_Px_aa+ABY*I_ESP_I6z_Px_aa;
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
     * shell quartet name: SQ_ESP_H_F_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_aa
     ************************************************************/
    Double I_ESP_H5x_F3x_aa = I_ESP_I6x_D2x_aa+ABX*I_ESP_H5x_D2x_aa;
    Double I_ESP_H4xy_F3x_aa = I_ESP_I5xy_D2x_aa+ABX*I_ESP_H4xy_D2x_aa;
    Double I_ESP_H4xz_F3x_aa = I_ESP_I5xz_D2x_aa+ABX*I_ESP_H4xz_D2x_aa;
    Double I_ESP_H3x2y_F3x_aa = I_ESP_I4x2y_D2x_aa+ABX*I_ESP_H3x2y_D2x_aa;
    Double I_ESP_H3xyz_F3x_aa = I_ESP_I4xyz_D2x_aa+ABX*I_ESP_H3xyz_D2x_aa;
    Double I_ESP_H3x2z_F3x_aa = I_ESP_I4x2z_D2x_aa+ABX*I_ESP_H3x2z_D2x_aa;
    Double I_ESP_H2x3y_F3x_aa = I_ESP_I3x3y_D2x_aa+ABX*I_ESP_H2x3y_D2x_aa;
    Double I_ESP_H2x2yz_F3x_aa = I_ESP_I3x2yz_D2x_aa+ABX*I_ESP_H2x2yz_D2x_aa;
    Double I_ESP_H2xy2z_F3x_aa = I_ESP_I3xy2z_D2x_aa+ABX*I_ESP_H2xy2z_D2x_aa;
    Double I_ESP_H2x3z_F3x_aa = I_ESP_I3x3z_D2x_aa+ABX*I_ESP_H2x3z_D2x_aa;
    Double I_ESP_Hx4y_F3x_aa = I_ESP_I2x4y_D2x_aa+ABX*I_ESP_Hx4y_D2x_aa;
    Double I_ESP_Hx3yz_F3x_aa = I_ESP_I2x3yz_D2x_aa+ABX*I_ESP_Hx3yz_D2x_aa;
    Double I_ESP_Hx2y2z_F3x_aa = I_ESP_I2x2y2z_D2x_aa+ABX*I_ESP_Hx2y2z_D2x_aa;
    Double I_ESP_Hxy3z_F3x_aa = I_ESP_I2xy3z_D2x_aa+ABX*I_ESP_Hxy3z_D2x_aa;
    Double I_ESP_Hx4z_F3x_aa = I_ESP_I2x4z_D2x_aa+ABX*I_ESP_Hx4z_D2x_aa;
    Double I_ESP_H5y_F3x_aa = I_ESP_Ix5y_D2x_aa+ABX*I_ESP_H5y_D2x_aa;
    Double I_ESP_H4yz_F3x_aa = I_ESP_Ix4yz_D2x_aa+ABX*I_ESP_H4yz_D2x_aa;
    Double I_ESP_H3y2z_F3x_aa = I_ESP_Ix3y2z_D2x_aa+ABX*I_ESP_H3y2z_D2x_aa;
    Double I_ESP_H2y3z_F3x_aa = I_ESP_Ix2y3z_D2x_aa+ABX*I_ESP_H2y3z_D2x_aa;
    Double I_ESP_Hy4z_F3x_aa = I_ESP_Ixy4z_D2x_aa+ABX*I_ESP_Hy4z_D2x_aa;
    Double I_ESP_H5z_F3x_aa = I_ESP_Ix5z_D2x_aa+ABX*I_ESP_H5z_D2x_aa;
    Double I_ESP_H5x_F2xy_aa = I_ESP_I5xy_D2x_aa+ABY*I_ESP_H5x_D2x_aa;
    Double I_ESP_H4xy_F2xy_aa = I_ESP_I4x2y_D2x_aa+ABY*I_ESP_H4xy_D2x_aa;
    Double I_ESP_H4xz_F2xy_aa = I_ESP_I4xyz_D2x_aa+ABY*I_ESP_H4xz_D2x_aa;
    Double I_ESP_H3x2y_F2xy_aa = I_ESP_I3x3y_D2x_aa+ABY*I_ESP_H3x2y_D2x_aa;
    Double I_ESP_H3xyz_F2xy_aa = I_ESP_I3x2yz_D2x_aa+ABY*I_ESP_H3xyz_D2x_aa;
    Double I_ESP_H3x2z_F2xy_aa = I_ESP_I3xy2z_D2x_aa+ABY*I_ESP_H3x2z_D2x_aa;
    Double I_ESP_H2x3y_F2xy_aa = I_ESP_I2x4y_D2x_aa+ABY*I_ESP_H2x3y_D2x_aa;
    Double I_ESP_H2x2yz_F2xy_aa = I_ESP_I2x3yz_D2x_aa+ABY*I_ESP_H2x2yz_D2x_aa;
    Double I_ESP_H2xy2z_F2xy_aa = I_ESP_I2x2y2z_D2x_aa+ABY*I_ESP_H2xy2z_D2x_aa;
    Double I_ESP_H2x3z_F2xy_aa = I_ESP_I2xy3z_D2x_aa+ABY*I_ESP_H2x3z_D2x_aa;
    Double I_ESP_Hx4y_F2xy_aa = I_ESP_Ix5y_D2x_aa+ABY*I_ESP_Hx4y_D2x_aa;
    Double I_ESP_Hx3yz_F2xy_aa = I_ESP_Ix4yz_D2x_aa+ABY*I_ESP_Hx3yz_D2x_aa;
    Double I_ESP_Hx2y2z_F2xy_aa = I_ESP_Ix3y2z_D2x_aa+ABY*I_ESP_Hx2y2z_D2x_aa;
    Double I_ESP_Hxy3z_F2xy_aa = I_ESP_Ix2y3z_D2x_aa+ABY*I_ESP_Hxy3z_D2x_aa;
    Double I_ESP_Hx4z_F2xy_aa = I_ESP_Ixy4z_D2x_aa+ABY*I_ESP_Hx4z_D2x_aa;
    Double I_ESP_H5y_F2xy_aa = I_ESP_I6y_D2x_aa+ABY*I_ESP_H5y_D2x_aa;
    Double I_ESP_H4yz_F2xy_aa = I_ESP_I5yz_D2x_aa+ABY*I_ESP_H4yz_D2x_aa;
    Double I_ESP_H3y2z_F2xy_aa = I_ESP_I4y2z_D2x_aa+ABY*I_ESP_H3y2z_D2x_aa;
    Double I_ESP_H2y3z_F2xy_aa = I_ESP_I3y3z_D2x_aa+ABY*I_ESP_H2y3z_D2x_aa;
    Double I_ESP_Hy4z_F2xy_aa = I_ESP_I2y4z_D2x_aa+ABY*I_ESP_Hy4z_D2x_aa;
    Double I_ESP_H5z_F2xy_aa = I_ESP_Iy5z_D2x_aa+ABY*I_ESP_H5z_D2x_aa;
    Double I_ESP_H5x_F2xz_aa = I_ESP_I5xz_D2x_aa+ABZ*I_ESP_H5x_D2x_aa;
    Double I_ESP_H4xy_F2xz_aa = I_ESP_I4xyz_D2x_aa+ABZ*I_ESP_H4xy_D2x_aa;
    Double I_ESP_H4xz_F2xz_aa = I_ESP_I4x2z_D2x_aa+ABZ*I_ESP_H4xz_D2x_aa;
    Double I_ESP_H3x2y_F2xz_aa = I_ESP_I3x2yz_D2x_aa+ABZ*I_ESP_H3x2y_D2x_aa;
    Double I_ESP_H3xyz_F2xz_aa = I_ESP_I3xy2z_D2x_aa+ABZ*I_ESP_H3xyz_D2x_aa;
    Double I_ESP_H3x2z_F2xz_aa = I_ESP_I3x3z_D2x_aa+ABZ*I_ESP_H3x2z_D2x_aa;
    Double I_ESP_H2x3y_F2xz_aa = I_ESP_I2x3yz_D2x_aa+ABZ*I_ESP_H2x3y_D2x_aa;
    Double I_ESP_H2x2yz_F2xz_aa = I_ESP_I2x2y2z_D2x_aa+ABZ*I_ESP_H2x2yz_D2x_aa;
    Double I_ESP_H2xy2z_F2xz_aa = I_ESP_I2xy3z_D2x_aa+ABZ*I_ESP_H2xy2z_D2x_aa;
    Double I_ESP_H2x3z_F2xz_aa = I_ESP_I2x4z_D2x_aa+ABZ*I_ESP_H2x3z_D2x_aa;
    Double I_ESP_Hx4y_F2xz_aa = I_ESP_Ix4yz_D2x_aa+ABZ*I_ESP_Hx4y_D2x_aa;
    Double I_ESP_Hx3yz_F2xz_aa = I_ESP_Ix3y2z_D2x_aa+ABZ*I_ESP_Hx3yz_D2x_aa;
    Double I_ESP_Hx2y2z_F2xz_aa = I_ESP_Ix2y3z_D2x_aa+ABZ*I_ESP_Hx2y2z_D2x_aa;
    Double I_ESP_Hxy3z_F2xz_aa = I_ESP_Ixy4z_D2x_aa+ABZ*I_ESP_Hxy3z_D2x_aa;
    Double I_ESP_Hx4z_F2xz_aa = I_ESP_Ix5z_D2x_aa+ABZ*I_ESP_Hx4z_D2x_aa;
    Double I_ESP_H5y_F2xz_aa = I_ESP_I5yz_D2x_aa+ABZ*I_ESP_H5y_D2x_aa;
    Double I_ESP_H4yz_F2xz_aa = I_ESP_I4y2z_D2x_aa+ABZ*I_ESP_H4yz_D2x_aa;
    Double I_ESP_H3y2z_F2xz_aa = I_ESP_I3y3z_D2x_aa+ABZ*I_ESP_H3y2z_D2x_aa;
    Double I_ESP_H2y3z_F2xz_aa = I_ESP_I2y4z_D2x_aa+ABZ*I_ESP_H2y3z_D2x_aa;
    Double I_ESP_Hy4z_F2xz_aa = I_ESP_Iy5z_D2x_aa+ABZ*I_ESP_Hy4z_D2x_aa;
    Double I_ESP_H5z_F2xz_aa = I_ESP_I6z_D2x_aa+ABZ*I_ESP_H5z_D2x_aa;
    Double I_ESP_H5x_Fx2y_aa = I_ESP_I6x_D2y_aa+ABX*I_ESP_H5x_D2y_aa;
    Double I_ESP_H4xy_Fx2y_aa = I_ESP_I5xy_D2y_aa+ABX*I_ESP_H4xy_D2y_aa;
    Double I_ESP_H4xz_Fx2y_aa = I_ESP_I5xz_D2y_aa+ABX*I_ESP_H4xz_D2y_aa;
    Double I_ESP_H3x2y_Fx2y_aa = I_ESP_I4x2y_D2y_aa+ABX*I_ESP_H3x2y_D2y_aa;
    Double I_ESP_H3xyz_Fx2y_aa = I_ESP_I4xyz_D2y_aa+ABX*I_ESP_H3xyz_D2y_aa;
    Double I_ESP_H3x2z_Fx2y_aa = I_ESP_I4x2z_D2y_aa+ABX*I_ESP_H3x2z_D2y_aa;
    Double I_ESP_H2x3y_Fx2y_aa = I_ESP_I3x3y_D2y_aa+ABX*I_ESP_H2x3y_D2y_aa;
    Double I_ESP_H2x2yz_Fx2y_aa = I_ESP_I3x2yz_D2y_aa+ABX*I_ESP_H2x2yz_D2y_aa;
    Double I_ESP_H2xy2z_Fx2y_aa = I_ESP_I3xy2z_D2y_aa+ABX*I_ESP_H2xy2z_D2y_aa;
    Double I_ESP_H2x3z_Fx2y_aa = I_ESP_I3x3z_D2y_aa+ABX*I_ESP_H2x3z_D2y_aa;
    Double I_ESP_Hx4y_Fx2y_aa = I_ESP_I2x4y_D2y_aa+ABX*I_ESP_Hx4y_D2y_aa;
    Double I_ESP_Hx3yz_Fx2y_aa = I_ESP_I2x3yz_D2y_aa+ABX*I_ESP_Hx3yz_D2y_aa;
    Double I_ESP_Hx2y2z_Fx2y_aa = I_ESP_I2x2y2z_D2y_aa+ABX*I_ESP_Hx2y2z_D2y_aa;
    Double I_ESP_Hxy3z_Fx2y_aa = I_ESP_I2xy3z_D2y_aa+ABX*I_ESP_Hxy3z_D2y_aa;
    Double I_ESP_Hx4z_Fx2y_aa = I_ESP_I2x4z_D2y_aa+ABX*I_ESP_Hx4z_D2y_aa;
    Double I_ESP_H5y_Fx2y_aa = I_ESP_Ix5y_D2y_aa+ABX*I_ESP_H5y_D2y_aa;
    Double I_ESP_H4yz_Fx2y_aa = I_ESP_Ix4yz_D2y_aa+ABX*I_ESP_H4yz_D2y_aa;
    Double I_ESP_H3y2z_Fx2y_aa = I_ESP_Ix3y2z_D2y_aa+ABX*I_ESP_H3y2z_D2y_aa;
    Double I_ESP_H2y3z_Fx2y_aa = I_ESP_Ix2y3z_D2y_aa+ABX*I_ESP_H2y3z_D2y_aa;
    Double I_ESP_Hy4z_Fx2y_aa = I_ESP_Ixy4z_D2y_aa+ABX*I_ESP_Hy4z_D2y_aa;
    Double I_ESP_H5z_Fx2y_aa = I_ESP_Ix5z_D2y_aa+ABX*I_ESP_H5z_D2y_aa;
    Double I_ESP_H5x_Fxyz_aa = I_ESP_I5xz_Dxy_aa+ABZ*I_ESP_H5x_Dxy_aa;
    Double I_ESP_H4xy_Fxyz_aa = I_ESP_I4xyz_Dxy_aa+ABZ*I_ESP_H4xy_Dxy_aa;
    Double I_ESP_H4xz_Fxyz_aa = I_ESP_I4x2z_Dxy_aa+ABZ*I_ESP_H4xz_Dxy_aa;
    Double I_ESP_H3x2y_Fxyz_aa = I_ESP_I3x2yz_Dxy_aa+ABZ*I_ESP_H3x2y_Dxy_aa;
    Double I_ESP_H3xyz_Fxyz_aa = I_ESP_I3xy2z_Dxy_aa+ABZ*I_ESP_H3xyz_Dxy_aa;
    Double I_ESP_H3x2z_Fxyz_aa = I_ESP_I3x3z_Dxy_aa+ABZ*I_ESP_H3x2z_Dxy_aa;
    Double I_ESP_H2x3y_Fxyz_aa = I_ESP_I2x3yz_Dxy_aa+ABZ*I_ESP_H2x3y_Dxy_aa;
    Double I_ESP_H2x2yz_Fxyz_aa = I_ESP_I2x2y2z_Dxy_aa+ABZ*I_ESP_H2x2yz_Dxy_aa;
    Double I_ESP_H2xy2z_Fxyz_aa = I_ESP_I2xy3z_Dxy_aa+ABZ*I_ESP_H2xy2z_Dxy_aa;
    Double I_ESP_H2x3z_Fxyz_aa = I_ESP_I2x4z_Dxy_aa+ABZ*I_ESP_H2x3z_Dxy_aa;
    Double I_ESP_Hx4y_Fxyz_aa = I_ESP_Ix4yz_Dxy_aa+ABZ*I_ESP_Hx4y_Dxy_aa;
    Double I_ESP_Hx3yz_Fxyz_aa = I_ESP_Ix3y2z_Dxy_aa+ABZ*I_ESP_Hx3yz_Dxy_aa;
    Double I_ESP_Hx2y2z_Fxyz_aa = I_ESP_Ix2y3z_Dxy_aa+ABZ*I_ESP_Hx2y2z_Dxy_aa;
    Double I_ESP_Hxy3z_Fxyz_aa = I_ESP_Ixy4z_Dxy_aa+ABZ*I_ESP_Hxy3z_Dxy_aa;
    Double I_ESP_Hx4z_Fxyz_aa = I_ESP_Ix5z_Dxy_aa+ABZ*I_ESP_Hx4z_Dxy_aa;
    Double I_ESP_H5y_Fxyz_aa = I_ESP_I5yz_Dxy_aa+ABZ*I_ESP_H5y_Dxy_aa;
    Double I_ESP_H4yz_Fxyz_aa = I_ESP_I4y2z_Dxy_aa+ABZ*I_ESP_H4yz_Dxy_aa;
    Double I_ESP_H3y2z_Fxyz_aa = I_ESP_I3y3z_Dxy_aa+ABZ*I_ESP_H3y2z_Dxy_aa;
    Double I_ESP_H2y3z_Fxyz_aa = I_ESP_I2y4z_Dxy_aa+ABZ*I_ESP_H2y3z_Dxy_aa;
    Double I_ESP_Hy4z_Fxyz_aa = I_ESP_Iy5z_Dxy_aa+ABZ*I_ESP_Hy4z_Dxy_aa;
    Double I_ESP_H5z_Fxyz_aa = I_ESP_I6z_Dxy_aa+ABZ*I_ESP_H5z_Dxy_aa;
    Double I_ESP_H5x_Fx2z_aa = I_ESP_I6x_D2z_aa+ABX*I_ESP_H5x_D2z_aa;
    Double I_ESP_H4xy_Fx2z_aa = I_ESP_I5xy_D2z_aa+ABX*I_ESP_H4xy_D2z_aa;
    Double I_ESP_H4xz_Fx2z_aa = I_ESP_I5xz_D2z_aa+ABX*I_ESP_H4xz_D2z_aa;
    Double I_ESP_H3x2y_Fx2z_aa = I_ESP_I4x2y_D2z_aa+ABX*I_ESP_H3x2y_D2z_aa;
    Double I_ESP_H3xyz_Fx2z_aa = I_ESP_I4xyz_D2z_aa+ABX*I_ESP_H3xyz_D2z_aa;
    Double I_ESP_H3x2z_Fx2z_aa = I_ESP_I4x2z_D2z_aa+ABX*I_ESP_H3x2z_D2z_aa;
    Double I_ESP_H2x3y_Fx2z_aa = I_ESP_I3x3y_D2z_aa+ABX*I_ESP_H2x3y_D2z_aa;
    Double I_ESP_H2x2yz_Fx2z_aa = I_ESP_I3x2yz_D2z_aa+ABX*I_ESP_H2x2yz_D2z_aa;
    Double I_ESP_H2xy2z_Fx2z_aa = I_ESP_I3xy2z_D2z_aa+ABX*I_ESP_H2xy2z_D2z_aa;
    Double I_ESP_H2x3z_Fx2z_aa = I_ESP_I3x3z_D2z_aa+ABX*I_ESP_H2x3z_D2z_aa;
    Double I_ESP_Hx4y_Fx2z_aa = I_ESP_I2x4y_D2z_aa+ABX*I_ESP_Hx4y_D2z_aa;
    Double I_ESP_Hx3yz_Fx2z_aa = I_ESP_I2x3yz_D2z_aa+ABX*I_ESP_Hx3yz_D2z_aa;
    Double I_ESP_Hx2y2z_Fx2z_aa = I_ESP_I2x2y2z_D2z_aa+ABX*I_ESP_Hx2y2z_D2z_aa;
    Double I_ESP_Hxy3z_Fx2z_aa = I_ESP_I2xy3z_D2z_aa+ABX*I_ESP_Hxy3z_D2z_aa;
    Double I_ESP_Hx4z_Fx2z_aa = I_ESP_I2x4z_D2z_aa+ABX*I_ESP_Hx4z_D2z_aa;
    Double I_ESP_H5y_Fx2z_aa = I_ESP_Ix5y_D2z_aa+ABX*I_ESP_H5y_D2z_aa;
    Double I_ESP_H4yz_Fx2z_aa = I_ESP_Ix4yz_D2z_aa+ABX*I_ESP_H4yz_D2z_aa;
    Double I_ESP_H3y2z_Fx2z_aa = I_ESP_Ix3y2z_D2z_aa+ABX*I_ESP_H3y2z_D2z_aa;
    Double I_ESP_H2y3z_Fx2z_aa = I_ESP_Ix2y3z_D2z_aa+ABX*I_ESP_H2y3z_D2z_aa;
    Double I_ESP_Hy4z_Fx2z_aa = I_ESP_Ixy4z_D2z_aa+ABX*I_ESP_Hy4z_D2z_aa;
    Double I_ESP_H5z_Fx2z_aa = I_ESP_Ix5z_D2z_aa+ABX*I_ESP_H5z_D2z_aa;
    Double I_ESP_H5x_F3y_aa = I_ESP_I5xy_D2y_aa+ABY*I_ESP_H5x_D2y_aa;
    Double I_ESP_H4xy_F3y_aa = I_ESP_I4x2y_D2y_aa+ABY*I_ESP_H4xy_D2y_aa;
    Double I_ESP_H4xz_F3y_aa = I_ESP_I4xyz_D2y_aa+ABY*I_ESP_H4xz_D2y_aa;
    Double I_ESP_H3x2y_F3y_aa = I_ESP_I3x3y_D2y_aa+ABY*I_ESP_H3x2y_D2y_aa;
    Double I_ESP_H3xyz_F3y_aa = I_ESP_I3x2yz_D2y_aa+ABY*I_ESP_H3xyz_D2y_aa;
    Double I_ESP_H3x2z_F3y_aa = I_ESP_I3xy2z_D2y_aa+ABY*I_ESP_H3x2z_D2y_aa;
    Double I_ESP_H2x3y_F3y_aa = I_ESP_I2x4y_D2y_aa+ABY*I_ESP_H2x3y_D2y_aa;
    Double I_ESP_H2x2yz_F3y_aa = I_ESP_I2x3yz_D2y_aa+ABY*I_ESP_H2x2yz_D2y_aa;
    Double I_ESP_H2xy2z_F3y_aa = I_ESP_I2x2y2z_D2y_aa+ABY*I_ESP_H2xy2z_D2y_aa;
    Double I_ESP_H2x3z_F3y_aa = I_ESP_I2xy3z_D2y_aa+ABY*I_ESP_H2x3z_D2y_aa;
    Double I_ESP_Hx4y_F3y_aa = I_ESP_Ix5y_D2y_aa+ABY*I_ESP_Hx4y_D2y_aa;
    Double I_ESP_Hx3yz_F3y_aa = I_ESP_Ix4yz_D2y_aa+ABY*I_ESP_Hx3yz_D2y_aa;
    Double I_ESP_Hx2y2z_F3y_aa = I_ESP_Ix3y2z_D2y_aa+ABY*I_ESP_Hx2y2z_D2y_aa;
    Double I_ESP_Hxy3z_F3y_aa = I_ESP_Ix2y3z_D2y_aa+ABY*I_ESP_Hxy3z_D2y_aa;
    Double I_ESP_Hx4z_F3y_aa = I_ESP_Ixy4z_D2y_aa+ABY*I_ESP_Hx4z_D2y_aa;
    Double I_ESP_H5y_F3y_aa = I_ESP_I6y_D2y_aa+ABY*I_ESP_H5y_D2y_aa;
    Double I_ESP_H4yz_F3y_aa = I_ESP_I5yz_D2y_aa+ABY*I_ESP_H4yz_D2y_aa;
    Double I_ESP_H3y2z_F3y_aa = I_ESP_I4y2z_D2y_aa+ABY*I_ESP_H3y2z_D2y_aa;
    Double I_ESP_H2y3z_F3y_aa = I_ESP_I3y3z_D2y_aa+ABY*I_ESP_H2y3z_D2y_aa;
    Double I_ESP_Hy4z_F3y_aa = I_ESP_I2y4z_D2y_aa+ABY*I_ESP_Hy4z_D2y_aa;
    Double I_ESP_H5z_F3y_aa = I_ESP_Iy5z_D2y_aa+ABY*I_ESP_H5z_D2y_aa;
    Double I_ESP_H5x_F2yz_aa = I_ESP_I5xz_D2y_aa+ABZ*I_ESP_H5x_D2y_aa;
    Double I_ESP_H4xy_F2yz_aa = I_ESP_I4xyz_D2y_aa+ABZ*I_ESP_H4xy_D2y_aa;
    Double I_ESP_H4xz_F2yz_aa = I_ESP_I4x2z_D2y_aa+ABZ*I_ESP_H4xz_D2y_aa;
    Double I_ESP_H3x2y_F2yz_aa = I_ESP_I3x2yz_D2y_aa+ABZ*I_ESP_H3x2y_D2y_aa;
    Double I_ESP_H3xyz_F2yz_aa = I_ESP_I3xy2z_D2y_aa+ABZ*I_ESP_H3xyz_D2y_aa;
    Double I_ESP_H3x2z_F2yz_aa = I_ESP_I3x3z_D2y_aa+ABZ*I_ESP_H3x2z_D2y_aa;
    Double I_ESP_H2x3y_F2yz_aa = I_ESP_I2x3yz_D2y_aa+ABZ*I_ESP_H2x3y_D2y_aa;
    Double I_ESP_H2x2yz_F2yz_aa = I_ESP_I2x2y2z_D2y_aa+ABZ*I_ESP_H2x2yz_D2y_aa;
    Double I_ESP_H2xy2z_F2yz_aa = I_ESP_I2xy3z_D2y_aa+ABZ*I_ESP_H2xy2z_D2y_aa;
    Double I_ESP_H2x3z_F2yz_aa = I_ESP_I2x4z_D2y_aa+ABZ*I_ESP_H2x3z_D2y_aa;
    Double I_ESP_Hx4y_F2yz_aa = I_ESP_Ix4yz_D2y_aa+ABZ*I_ESP_Hx4y_D2y_aa;
    Double I_ESP_Hx3yz_F2yz_aa = I_ESP_Ix3y2z_D2y_aa+ABZ*I_ESP_Hx3yz_D2y_aa;
    Double I_ESP_Hx2y2z_F2yz_aa = I_ESP_Ix2y3z_D2y_aa+ABZ*I_ESP_Hx2y2z_D2y_aa;
    Double I_ESP_Hxy3z_F2yz_aa = I_ESP_Ixy4z_D2y_aa+ABZ*I_ESP_Hxy3z_D2y_aa;
    Double I_ESP_Hx4z_F2yz_aa = I_ESP_Ix5z_D2y_aa+ABZ*I_ESP_Hx4z_D2y_aa;
    Double I_ESP_H5y_F2yz_aa = I_ESP_I5yz_D2y_aa+ABZ*I_ESP_H5y_D2y_aa;
    Double I_ESP_H4yz_F2yz_aa = I_ESP_I4y2z_D2y_aa+ABZ*I_ESP_H4yz_D2y_aa;
    Double I_ESP_H3y2z_F2yz_aa = I_ESP_I3y3z_D2y_aa+ABZ*I_ESP_H3y2z_D2y_aa;
    Double I_ESP_H2y3z_F2yz_aa = I_ESP_I2y4z_D2y_aa+ABZ*I_ESP_H2y3z_D2y_aa;
    Double I_ESP_Hy4z_F2yz_aa = I_ESP_Iy5z_D2y_aa+ABZ*I_ESP_Hy4z_D2y_aa;
    Double I_ESP_H5z_F2yz_aa = I_ESP_I6z_D2y_aa+ABZ*I_ESP_H5z_D2y_aa;
    Double I_ESP_H5x_Fy2z_aa = I_ESP_I5xy_D2z_aa+ABY*I_ESP_H5x_D2z_aa;
    Double I_ESP_H4xy_Fy2z_aa = I_ESP_I4x2y_D2z_aa+ABY*I_ESP_H4xy_D2z_aa;
    Double I_ESP_H4xz_Fy2z_aa = I_ESP_I4xyz_D2z_aa+ABY*I_ESP_H4xz_D2z_aa;
    Double I_ESP_H3x2y_Fy2z_aa = I_ESP_I3x3y_D2z_aa+ABY*I_ESP_H3x2y_D2z_aa;
    Double I_ESP_H3xyz_Fy2z_aa = I_ESP_I3x2yz_D2z_aa+ABY*I_ESP_H3xyz_D2z_aa;
    Double I_ESP_H3x2z_Fy2z_aa = I_ESP_I3xy2z_D2z_aa+ABY*I_ESP_H3x2z_D2z_aa;
    Double I_ESP_H2x3y_Fy2z_aa = I_ESP_I2x4y_D2z_aa+ABY*I_ESP_H2x3y_D2z_aa;
    Double I_ESP_H2x2yz_Fy2z_aa = I_ESP_I2x3yz_D2z_aa+ABY*I_ESP_H2x2yz_D2z_aa;
    Double I_ESP_H2xy2z_Fy2z_aa = I_ESP_I2x2y2z_D2z_aa+ABY*I_ESP_H2xy2z_D2z_aa;
    Double I_ESP_H2x3z_Fy2z_aa = I_ESP_I2xy3z_D2z_aa+ABY*I_ESP_H2x3z_D2z_aa;
    Double I_ESP_Hx4y_Fy2z_aa = I_ESP_Ix5y_D2z_aa+ABY*I_ESP_Hx4y_D2z_aa;
    Double I_ESP_Hx3yz_Fy2z_aa = I_ESP_Ix4yz_D2z_aa+ABY*I_ESP_Hx3yz_D2z_aa;
    Double I_ESP_Hx2y2z_Fy2z_aa = I_ESP_Ix3y2z_D2z_aa+ABY*I_ESP_Hx2y2z_D2z_aa;
    Double I_ESP_Hxy3z_Fy2z_aa = I_ESP_Ix2y3z_D2z_aa+ABY*I_ESP_Hxy3z_D2z_aa;
    Double I_ESP_Hx4z_Fy2z_aa = I_ESP_Ixy4z_D2z_aa+ABY*I_ESP_Hx4z_D2z_aa;
    Double I_ESP_H5y_Fy2z_aa = I_ESP_I6y_D2z_aa+ABY*I_ESP_H5y_D2z_aa;
    Double I_ESP_H4yz_Fy2z_aa = I_ESP_I5yz_D2z_aa+ABY*I_ESP_H4yz_D2z_aa;
    Double I_ESP_H3y2z_Fy2z_aa = I_ESP_I4y2z_D2z_aa+ABY*I_ESP_H3y2z_D2z_aa;
    Double I_ESP_H2y3z_Fy2z_aa = I_ESP_I3y3z_D2z_aa+ABY*I_ESP_H2y3z_D2z_aa;
    Double I_ESP_Hy4z_Fy2z_aa = I_ESP_I2y4z_D2z_aa+ABY*I_ESP_Hy4z_D2z_aa;
    Double I_ESP_H5z_Fy2z_aa = I_ESP_Iy5z_D2z_aa+ABY*I_ESP_H5z_D2z_aa;
    Double I_ESP_H5x_F3z_aa = I_ESP_I5xz_D2z_aa+ABZ*I_ESP_H5x_D2z_aa;
    Double I_ESP_H4xy_F3z_aa = I_ESP_I4xyz_D2z_aa+ABZ*I_ESP_H4xy_D2z_aa;
    Double I_ESP_H4xz_F3z_aa = I_ESP_I4x2z_D2z_aa+ABZ*I_ESP_H4xz_D2z_aa;
    Double I_ESP_H3x2y_F3z_aa = I_ESP_I3x2yz_D2z_aa+ABZ*I_ESP_H3x2y_D2z_aa;
    Double I_ESP_H3xyz_F3z_aa = I_ESP_I3xy2z_D2z_aa+ABZ*I_ESP_H3xyz_D2z_aa;
    Double I_ESP_H3x2z_F3z_aa = I_ESP_I3x3z_D2z_aa+ABZ*I_ESP_H3x2z_D2z_aa;
    Double I_ESP_H2x3y_F3z_aa = I_ESP_I2x3yz_D2z_aa+ABZ*I_ESP_H2x3y_D2z_aa;
    Double I_ESP_H2x2yz_F3z_aa = I_ESP_I2x2y2z_D2z_aa+ABZ*I_ESP_H2x2yz_D2z_aa;
    Double I_ESP_H2xy2z_F3z_aa = I_ESP_I2xy3z_D2z_aa+ABZ*I_ESP_H2xy2z_D2z_aa;
    Double I_ESP_H2x3z_F3z_aa = I_ESP_I2x4z_D2z_aa+ABZ*I_ESP_H2x3z_D2z_aa;
    Double I_ESP_Hx4y_F3z_aa = I_ESP_Ix4yz_D2z_aa+ABZ*I_ESP_Hx4y_D2z_aa;
    Double I_ESP_Hx3yz_F3z_aa = I_ESP_Ix3y2z_D2z_aa+ABZ*I_ESP_Hx3yz_D2z_aa;
    Double I_ESP_Hx2y2z_F3z_aa = I_ESP_Ix2y3z_D2z_aa+ABZ*I_ESP_Hx2y2z_D2z_aa;
    Double I_ESP_Hxy3z_F3z_aa = I_ESP_Ixy4z_D2z_aa+ABZ*I_ESP_Hxy3z_D2z_aa;
    Double I_ESP_Hx4z_F3z_aa = I_ESP_Ix5z_D2z_aa+ABZ*I_ESP_Hx4z_D2z_aa;
    Double I_ESP_H5y_F3z_aa = I_ESP_I5yz_D2z_aa+ABZ*I_ESP_H5y_D2z_aa;
    Double I_ESP_H4yz_F3z_aa = I_ESP_I4y2z_D2z_aa+ABZ*I_ESP_H4yz_D2z_aa;
    Double I_ESP_H3y2z_F3z_aa = I_ESP_I3y3z_D2z_aa+ABZ*I_ESP_H3y2z_D2z_aa;
    Double I_ESP_H2y3z_F3z_aa = I_ESP_I2y4z_D2z_aa+ABZ*I_ESP_H2y3z_D2z_aa;
    Double I_ESP_Hy4z_F3z_aa = I_ESP_Iy5z_D2z_aa+ABZ*I_ESP_Hy4z_D2z_aa;
    Double I_ESP_H5z_F3z_aa = I_ESP_I6z_D2z_aa+ABZ*I_ESP_H5z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+0] = 4.0E0*I_ESP_H5x_F3x_aa-2.0E0*3*I_ESP_F3x_F3x_a-2.0E0*4*I_ESP_F3x_F3x_a+3*2*I_ESP_Px_F3x;
    abcd[iGrid*600+1] = 4.0E0*I_ESP_H4xy_F3x_aa-2.0E0*2*I_ESP_F2xy_F3x_a-2.0E0*3*I_ESP_F2xy_F3x_a+2*1*I_ESP_Py_F3x;
    abcd[iGrid*600+2] = 4.0E0*I_ESP_H4xz_F3x_aa-2.0E0*2*I_ESP_F2xz_F3x_a-2.0E0*3*I_ESP_F2xz_F3x_a+2*1*I_ESP_Pz_F3x;
    abcd[iGrid*600+3] = 4.0E0*I_ESP_H3x2y_F3x_aa-2.0E0*1*I_ESP_Fx2y_F3x_a-2.0E0*2*I_ESP_Fx2y_F3x_a;
    abcd[iGrid*600+4] = 4.0E0*I_ESP_H3xyz_F3x_aa-2.0E0*1*I_ESP_Fxyz_F3x_a-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+5] = 4.0E0*I_ESP_H3x2z_F3x_aa-2.0E0*1*I_ESP_Fx2z_F3x_a-2.0E0*2*I_ESP_Fx2z_F3x_a;
    abcd[iGrid*600+6] = 4.0E0*I_ESP_H2x3y_F3x_aa-2.0E0*1*I_ESP_F3y_F3x_a;
    abcd[iGrid*600+7] = 4.0E0*I_ESP_H2x2yz_F3x_aa-2.0E0*1*I_ESP_F2yz_F3x_a;
    abcd[iGrid*600+8] = 4.0E0*I_ESP_H2xy2z_F3x_aa-2.0E0*1*I_ESP_Fy2z_F3x_a;
    abcd[iGrid*600+9] = 4.0E0*I_ESP_H2x3z_F3x_aa-2.0E0*1*I_ESP_F3z_F3x_a;
    abcd[iGrid*600+10] = 4.0E0*I_ESP_H5x_F2xy_aa-2.0E0*3*I_ESP_F3x_F2xy_a-2.0E0*4*I_ESP_F3x_F2xy_a+3*2*I_ESP_Px_F2xy;
    abcd[iGrid*600+11] = 4.0E0*I_ESP_H4xy_F2xy_aa-2.0E0*2*I_ESP_F2xy_F2xy_a-2.0E0*3*I_ESP_F2xy_F2xy_a+2*1*I_ESP_Py_F2xy;
    abcd[iGrid*600+12] = 4.0E0*I_ESP_H4xz_F2xy_aa-2.0E0*2*I_ESP_F2xz_F2xy_a-2.0E0*3*I_ESP_F2xz_F2xy_a+2*1*I_ESP_Pz_F2xy;
    abcd[iGrid*600+13] = 4.0E0*I_ESP_H3x2y_F2xy_aa-2.0E0*1*I_ESP_Fx2y_F2xy_a-2.0E0*2*I_ESP_Fx2y_F2xy_a;
    abcd[iGrid*600+14] = 4.0E0*I_ESP_H3xyz_F2xy_aa-2.0E0*1*I_ESP_Fxyz_F2xy_a-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+15] = 4.0E0*I_ESP_H3x2z_F2xy_aa-2.0E0*1*I_ESP_Fx2z_F2xy_a-2.0E0*2*I_ESP_Fx2z_F2xy_a;
    abcd[iGrid*600+16] = 4.0E0*I_ESP_H2x3y_F2xy_aa-2.0E0*1*I_ESP_F3y_F2xy_a;
    abcd[iGrid*600+17] = 4.0E0*I_ESP_H2x2yz_F2xy_aa-2.0E0*1*I_ESP_F2yz_F2xy_a;
    abcd[iGrid*600+18] = 4.0E0*I_ESP_H2xy2z_F2xy_aa-2.0E0*1*I_ESP_Fy2z_F2xy_a;
    abcd[iGrid*600+19] = 4.0E0*I_ESP_H2x3z_F2xy_aa-2.0E0*1*I_ESP_F3z_F2xy_a;
    abcd[iGrid*600+20] = 4.0E0*I_ESP_H5x_F2xz_aa-2.0E0*3*I_ESP_F3x_F2xz_a-2.0E0*4*I_ESP_F3x_F2xz_a+3*2*I_ESP_Px_F2xz;
    abcd[iGrid*600+21] = 4.0E0*I_ESP_H4xy_F2xz_aa-2.0E0*2*I_ESP_F2xy_F2xz_a-2.0E0*3*I_ESP_F2xy_F2xz_a+2*1*I_ESP_Py_F2xz;
    abcd[iGrid*600+22] = 4.0E0*I_ESP_H4xz_F2xz_aa-2.0E0*2*I_ESP_F2xz_F2xz_a-2.0E0*3*I_ESP_F2xz_F2xz_a+2*1*I_ESP_Pz_F2xz;
    abcd[iGrid*600+23] = 4.0E0*I_ESP_H3x2y_F2xz_aa-2.0E0*1*I_ESP_Fx2y_F2xz_a-2.0E0*2*I_ESP_Fx2y_F2xz_a;
    abcd[iGrid*600+24] = 4.0E0*I_ESP_H3xyz_F2xz_aa-2.0E0*1*I_ESP_Fxyz_F2xz_a-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+25] = 4.0E0*I_ESP_H3x2z_F2xz_aa-2.0E0*1*I_ESP_Fx2z_F2xz_a-2.0E0*2*I_ESP_Fx2z_F2xz_a;
    abcd[iGrid*600+26] = 4.0E0*I_ESP_H2x3y_F2xz_aa-2.0E0*1*I_ESP_F3y_F2xz_a;
    abcd[iGrid*600+27] = 4.0E0*I_ESP_H2x2yz_F2xz_aa-2.0E0*1*I_ESP_F2yz_F2xz_a;
    abcd[iGrid*600+28] = 4.0E0*I_ESP_H2xy2z_F2xz_aa-2.0E0*1*I_ESP_Fy2z_F2xz_a;
    abcd[iGrid*600+29] = 4.0E0*I_ESP_H2x3z_F2xz_aa-2.0E0*1*I_ESP_F3z_F2xz_a;
    abcd[iGrid*600+30] = 4.0E0*I_ESP_H5x_Fx2y_aa-2.0E0*3*I_ESP_F3x_Fx2y_a-2.0E0*4*I_ESP_F3x_Fx2y_a+3*2*I_ESP_Px_Fx2y;
    abcd[iGrid*600+31] = 4.0E0*I_ESP_H4xy_Fx2y_aa-2.0E0*2*I_ESP_F2xy_Fx2y_a-2.0E0*3*I_ESP_F2xy_Fx2y_a+2*1*I_ESP_Py_Fx2y;
    abcd[iGrid*600+32] = 4.0E0*I_ESP_H4xz_Fx2y_aa-2.0E0*2*I_ESP_F2xz_Fx2y_a-2.0E0*3*I_ESP_F2xz_Fx2y_a+2*1*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+33] = 4.0E0*I_ESP_H3x2y_Fx2y_aa-2.0E0*1*I_ESP_Fx2y_Fx2y_a-2.0E0*2*I_ESP_Fx2y_Fx2y_a;
    abcd[iGrid*600+34] = 4.0E0*I_ESP_H3xyz_Fx2y_aa-2.0E0*1*I_ESP_Fxyz_Fx2y_a-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+35] = 4.0E0*I_ESP_H3x2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2z_Fx2y_a-2.0E0*2*I_ESP_Fx2z_Fx2y_a;
    abcd[iGrid*600+36] = 4.0E0*I_ESP_H2x3y_Fx2y_aa-2.0E0*1*I_ESP_F3y_Fx2y_a;
    abcd[iGrid*600+37] = 4.0E0*I_ESP_H2x2yz_Fx2y_aa-2.0E0*1*I_ESP_F2yz_Fx2y_a;
    abcd[iGrid*600+38] = 4.0E0*I_ESP_H2xy2z_Fx2y_aa-2.0E0*1*I_ESP_Fy2z_Fx2y_a;
    abcd[iGrid*600+39] = 4.0E0*I_ESP_H2x3z_Fx2y_aa-2.0E0*1*I_ESP_F3z_Fx2y_a;
    abcd[iGrid*600+40] = 4.0E0*I_ESP_H5x_Fxyz_aa-2.0E0*3*I_ESP_F3x_Fxyz_a-2.0E0*4*I_ESP_F3x_Fxyz_a+3*2*I_ESP_Px_Fxyz;
    abcd[iGrid*600+41] = 4.0E0*I_ESP_H4xy_Fxyz_aa-2.0E0*2*I_ESP_F2xy_Fxyz_a-2.0E0*3*I_ESP_F2xy_Fxyz_a+2*1*I_ESP_Py_Fxyz;
    abcd[iGrid*600+42] = 4.0E0*I_ESP_H4xz_Fxyz_aa-2.0E0*2*I_ESP_F2xz_Fxyz_a-2.0E0*3*I_ESP_F2xz_Fxyz_a+2*1*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+43] = 4.0E0*I_ESP_H3x2y_Fxyz_aa-2.0E0*1*I_ESP_Fx2y_Fxyz_a-2.0E0*2*I_ESP_Fx2y_Fxyz_a;
    abcd[iGrid*600+44] = 4.0E0*I_ESP_H3xyz_Fxyz_aa-2.0E0*1*I_ESP_Fxyz_Fxyz_a-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+45] = 4.0E0*I_ESP_H3x2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2z_Fxyz_a-2.0E0*2*I_ESP_Fx2z_Fxyz_a;
    abcd[iGrid*600+46] = 4.0E0*I_ESP_H2x3y_Fxyz_aa-2.0E0*1*I_ESP_F3y_Fxyz_a;
    abcd[iGrid*600+47] = 4.0E0*I_ESP_H2x2yz_Fxyz_aa-2.0E0*1*I_ESP_F2yz_Fxyz_a;
    abcd[iGrid*600+48] = 4.0E0*I_ESP_H2xy2z_Fxyz_aa-2.0E0*1*I_ESP_Fy2z_Fxyz_a;
    abcd[iGrid*600+49] = 4.0E0*I_ESP_H2x3z_Fxyz_aa-2.0E0*1*I_ESP_F3z_Fxyz_a;
    abcd[iGrid*600+50] = 4.0E0*I_ESP_H5x_Fx2z_aa-2.0E0*3*I_ESP_F3x_Fx2z_a-2.0E0*4*I_ESP_F3x_Fx2z_a+3*2*I_ESP_Px_Fx2z;
    abcd[iGrid*600+51] = 4.0E0*I_ESP_H4xy_Fx2z_aa-2.0E0*2*I_ESP_F2xy_Fx2z_a-2.0E0*3*I_ESP_F2xy_Fx2z_a+2*1*I_ESP_Py_Fx2z;
    abcd[iGrid*600+52] = 4.0E0*I_ESP_H4xz_Fx2z_aa-2.0E0*2*I_ESP_F2xz_Fx2z_a-2.0E0*3*I_ESP_F2xz_Fx2z_a+2*1*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+53] = 4.0E0*I_ESP_H3x2y_Fx2z_aa-2.0E0*1*I_ESP_Fx2y_Fx2z_a-2.0E0*2*I_ESP_Fx2y_Fx2z_a;
    abcd[iGrid*600+54] = 4.0E0*I_ESP_H3xyz_Fx2z_aa-2.0E0*1*I_ESP_Fxyz_Fx2z_a-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+55] = 4.0E0*I_ESP_H3x2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2z_Fx2z_a-2.0E0*2*I_ESP_Fx2z_Fx2z_a;
    abcd[iGrid*600+56] = 4.0E0*I_ESP_H2x3y_Fx2z_aa-2.0E0*1*I_ESP_F3y_Fx2z_a;
    abcd[iGrid*600+57] = 4.0E0*I_ESP_H2x2yz_Fx2z_aa-2.0E0*1*I_ESP_F2yz_Fx2z_a;
    abcd[iGrid*600+58] = 4.0E0*I_ESP_H2xy2z_Fx2z_aa-2.0E0*1*I_ESP_Fy2z_Fx2z_a;
    abcd[iGrid*600+59] = 4.0E0*I_ESP_H2x3z_Fx2z_aa-2.0E0*1*I_ESP_F3z_Fx2z_a;
    abcd[iGrid*600+60] = 4.0E0*I_ESP_H5x_F3y_aa-2.0E0*3*I_ESP_F3x_F3y_a-2.0E0*4*I_ESP_F3x_F3y_a+3*2*I_ESP_Px_F3y;
    abcd[iGrid*600+61] = 4.0E0*I_ESP_H4xy_F3y_aa-2.0E0*2*I_ESP_F2xy_F3y_a-2.0E0*3*I_ESP_F2xy_F3y_a+2*1*I_ESP_Py_F3y;
    abcd[iGrid*600+62] = 4.0E0*I_ESP_H4xz_F3y_aa-2.0E0*2*I_ESP_F2xz_F3y_a-2.0E0*3*I_ESP_F2xz_F3y_a+2*1*I_ESP_Pz_F3y;
    abcd[iGrid*600+63] = 4.0E0*I_ESP_H3x2y_F3y_aa-2.0E0*1*I_ESP_Fx2y_F3y_a-2.0E0*2*I_ESP_Fx2y_F3y_a;
    abcd[iGrid*600+64] = 4.0E0*I_ESP_H3xyz_F3y_aa-2.0E0*1*I_ESP_Fxyz_F3y_a-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+65] = 4.0E0*I_ESP_H3x2z_F3y_aa-2.0E0*1*I_ESP_Fx2z_F3y_a-2.0E0*2*I_ESP_Fx2z_F3y_a;
    abcd[iGrid*600+66] = 4.0E0*I_ESP_H2x3y_F3y_aa-2.0E0*1*I_ESP_F3y_F3y_a;
    abcd[iGrid*600+67] = 4.0E0*I_ESP_H2x2yz_F3y_aa-2.0E0*1*I_ESP_F2yz_F3y_a;
    abcd[iGrid*600+68] = 4.0E0*I_ESP_H2xy2z_F3y_aa-2.0E0*1*I_ESP_Fy2z_F3y_a;
    abcd[iGrid*600+69] = 4.0E0*I_ESP_H2x3z_F3y_aa-2.0E0*1*I_ESP_F3z_F3y_a;
    abcd[iGrid*600+70] = 4.0E0*I_ESP_H5x_F2yz_aa-2.0E0*3*I_ESP_F3x_F2yz_a-2.0E0*4*I_ESP_F3x_F2yz_a+3*2*I_ESP_Px_F2yz;
    abcd[iGrid*600+71] = 4.0E0*I_ESP_H4xy_F2yz_aa-2.0E0*2*I_ESP_F2xy_F2yz_a-2.0E0*3*I_ESP_F2xy_F2yz_a+2*1*I_ESP_Py_F2yz;
    abcd[iGrid*600+72] = 4.0E0*I_ESP_H4xz_F2yz_aa-2.0E0*2*I_ESP_F2xz_F2yz_a-2.0E0*3*I_ESP_F2xz_F2yz_a+2*1*I_ESP_Pz_F2yz;
    abcd[iGrid*600+73] = 4.0E0*I_ESP_H3x2y_F2yz_aa-2.0E0*1*I_ESP_Fx2y_F2yz_a-2.0E0*2*I_ESP_Fx2y_F2yz_a;
    abcd[iGrid*600+74] = 4.0E0*I_ESP_H3xyz_F2yz_aa-2.0E0*1*I_ESP_Fxyz_F2yz_a-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+75] = 4.0E0*I_ESP_H3x2z_F2yz_aa-2.0E0*1*I_ESP_Fx2z_F2yz_a-2.0E0*2*I_ESP_Fx2z_F2yz_a;
    abcd[iGrid*600+76] = 4.0E0*I_ESP_H2x3y_F2yz_aa-2.0E0*1*I_ESP_F3y_F2yz_a;
    abcd[iGrid*600+77] = 4.0E0*I_ESP_H2x2yz_F2yz_aa-2.0E0*1*I_ESP_F2yz_F2yz_a;
    abcd[iGrid*600+78] = 4.0E0*I_ESP_H2xy2z_F2yz_aa-2.0E0*1*I_ESP_Fy2z_F2yz_a;
    abcd[iGrid*600+79] = 4.0E0*I_ESP_H2x3z_F2yz_aa-2.0E0*1*I_ESP_F3z_F2yz_a;
    abcd[iGrid*600+80] = 4.0E0*I_ESP_H5x_Fy2z_aa-2.0E0*3*I_ESP_F3x_Fy2z_a-2.0E0*4*I_ESP_F3x_Fy2z_a+3*2*I_ESP_Px_Fy2z;
    abcd[iGrid*600+81] = 4.0E0*I_ESP_H4xy_Fy2z_aa-2.0E0*2*I_ESP_F2xy_Fy2z_a-2.0E0*3*I_ESP_F2xy_Fy2z_a+2*1*I_ESP_Py_Fy2z;
    abcd[iGrid*600+82] = 4.0E0*I_ESP_H4xz_Fy2z_aa-2.0E0*2*I_ESP_F2xz_Fy2z_a-2.0E0*3*I_ESP_F2xz_Fy2z_a+2*1*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+83] = 4.0E0*I_ESP_H3x2y_Fy2z_aa-2.0E0*1*I_ESP_Fx2y_Fy2z_a-2.0E0*2*I_ESP_Fx2y_Fy2z_a;
    abcd[iGrid*600+84] = 4.0E0*I_ESP_H3xyz_Fy2z_aa-2.0E0*1*I_ESP_Fxyz_Fy2z_a-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+85] = 4.0E0*I_ESP_H3x2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2z_Fy2z_a-2.0E0*2*I_ESP_Fx2z_Fy2z_a;
    abcd[iGrid*600+86] = 4.0E0*I_ESP_H2x3y_Fy2z_aa-2.0E0*1*I_ESP_F3y_Fy2z_a;
    abcd[iGrid*600+87] = 4.0E0*I_ESP_H2x2yz_Fy2z_aa-2.0E0*1*I_ESP_F2yz_Fy2z_a;
    abcd[iGrid*600+88] = 4.0E0*I_ESP_H2xy2z_Fy2z_aa-2.0E0*1*I_ESP_Fy2z_Fy2z_a;
    abcd[iGrid*600+89] = 4.0E0*I_ESP_H2x3z_Fy2z_aa-2.0E0*1*I_ESP_F3z_Fy2z_a;
    abcd[iGrid*600+90] = 4.0E0*I_ESP_H5x_F3z_aa-2.0E0*3*I_ESP_F3x_F3z_a-2.0E0*4*I_ESP_F3x_F3z_a+3*2*I_ESP_Px_F3z;
    abcd[iGrid*600+91] = 4.0E0*I_ESP_H4xy_F3z_aa-2.0E0*2*I_ESP_F2xy_F3z_a-2.0E0*3*I_ESP_F2xy_F3z_a+2*1*I_ESP_Py_F3z;
    abcd[iGrid*600+92] = 4.0E0*I_ESP_H4xz_F3z_aa-2.0E0*2*I_ESP_F2xz_F3z_a-2.0E0*3*I_ESP_F2xz_F3z_a+2*1*I_ESP_Pz_F3z;
    abcd[iGrid*600+93] = 4.0E0*I_ESP_H3x2y_F3z_aa-2.0E0*1*I_ESP_Fx2y_F3z_a-2.0E0*2*I_ESP_Fx2y_F3z_a;
    abcd[iGrid*600+94] = 4.0E0*I_ESP_H3xyz_F3z_aa-2.0E0*1*I_ESP_Fxyz_F3z_a-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+95] = 4.0E0*I_ESP_H3x2z_F3z_aa-2.0E0*1*I_ESP_Fx2z_F3z_a-2.0E0*2*I_ESP_Fx2z_F3z_a;
    abcd[iGrid*600+96] = 4.0E0*I_ESP_H2x3y_F3z_aa-2.0E0*1*I_ESP_F3y_F3z_a;
    abcd[iGrid*600+97] = 4.0E0*I_ESP_H2x2yz_F3z_aa-2.0E0*1*I_ESP_F2yz_F3z_a;
    abcd[iGrid*600+98] = 4.0E0*I_ESP_H2xy2z_F3z_aa-2.0E0*1*I_ESP_Fy2z_F3z_a;
    abcd[iGrid*600+99] = 4.0E0*I_ESP_H2x3z_F3z_aa-2.0E0*1*I_ESP_F3z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+100] = 4.0E0*I_ESP_H4xy_F3x_aa-2.0E0*3*I_ESP_F2xy_F3x_a;
    abcd[iGrid*600+101] = 4.0E0*I_ESP_H3x2y_F3x_aa-2.0E0*1*I_ESP_F3x_F3x_a-2.0E0*2*I_ESP_Fx2y_F3x_a+2*1*I_ESP_Px_F3x;
    abcd[iGrid*600+102] = 4.0E0*I_ESP_H3xyz_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+103] = 4.0E0*I_ESP_H2x3y_F3x_aa-2.0E0*2*I_ESP_F2xy_F3x_a-2.0E0*1*I_ESP_F3y_F3x_a+2*I_ESP_Py_F3x;
    abcd[iGrid*600+104] = 4.0E0*I_ESP_H2x2yz_F3x_aa-2.0E0*1*I_ESP_F2xz_F3x_a-2.0E0*1*I_ESP_F2yz_F3x_a+1*I_ESP_Pz_F3x;
    abcd[iGrid*600+105] = 4.0E0*I_ESP_H2xy2z_F3x_aa-2.0E0*1*I_ESP_Fy2z_F3x_a;
    abcd[iGrid*600+106] = 4.0E0*I_ESP_Hx4y_F3x_aa-2.0E0*3*I_ESP_Fx2y_F3x_a;
    abcd[iGrid*600+107] = 4.0E0*I_ESP_Hx3yz_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+108] = 4.0E0*I_ESP_Hx2y2z_F3x_aa-2.0E0*1*I_ESP_Fx2z_F3x_a;
    abcd[iGrid*600+109] = 4.0E0*I_ESP_Hxy3z_F3x_aa;
    abcd[iGrid*600+110] = 4.0E0*I_ESP_H4xy_F2xy_aa-2.0E0*3*I_ESP_F2xy_F2xy_a;
    abcd[iGrid*600+111] = 4.0E0*I_ESP_H3x2y_F2xy_aa-2.0E0*1*I_ESP_F3x_F2xy_a-2.0E0*2*I_ESP_Fx2y_F2xy_a+2*1*I_ESP_Px_F2xy;
    abcd[iGrid*600+112] = 4.0E0*I_ESP_H3xyz_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+113] = 4.0E0*I_ESP_H2x3y_F2xy_aa-2.0E0*2*I_ESP_F2xy_F2xy_a-2.0E0*1*I_ESP_F3y_F2xy_a+2*I_ESP_Py_F2xy;
    abcd[iGrid*600+114] = 4.0E0*I_ESP_H2x2yz_F2xy_aa-2.0E0*1*I_ESP_F2xz_F2xy_a-2.0E0*1*I_ESP_F2yz_F2xy_a+1*I_ESP_Pz_F2xy;
    abcd[iGrid*600+115] = 4.0E0*I_ESP_H2xy2z_F2xy_aa-2.0E0*1*I_ESP_Fy2z_F2xy_a;
    abcd[iGrid*600+116] = 4.0E0*I_ESP_Hx4y_F2xy_aa-2.0E0*3*I_ESP_Fx2y_F2xy_a;
    abcd[iGrid*600+117] = 4.0E0*I_ESP_Hx3yz_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+118] = 4.0E0*I_ESP_Hx2y2z_F2xy_aa-2.0E0*1*I_ESP_Fx2z_F2xy_a;
    abcd[iGrid*600+119] = 4.0E0*I_ESP_Hxy3z_F2xy_aa;
    abcd[iGrid*600+120] = 4.0E0*I_ESP_H4xy_F2xz_aa-2.0E0*3*I_ESP_F2xy_F2xz_a;
    abcd[iGrid*600+121] = 4.0E0*I_ESP_H3x2y_F2xz_aa-2.0E0*1*I_ESP_F3x_F2xz_a-2.0E0*2*I_ESP_Fx2y_F2xz_a+2*1*I_ESP_Px_F2xz;
    abcd[iGrid*600+122] = 4.0E0*I_ESP_H3xyz_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+123] = 4.0E0*I_ESP_H2x3y_F2xz_aa-2.0E0*2*I_ESP_F2xy_F2xz_a-2.0E0*1*I_ESP_F3y_F2xz_a+2*I_ESP_Py_F2xz;
    abcd[iGrid*600+124] = 4.0E0*I_ESP_H2x2yz_F2xz_aa-2.0E0*1*I_ESP_F2xz_F2xz_a-2.0E0*1*I_ESP_F2yz_F2xz_a+1*I_ESP_Pz_F2xz;
    abcd[iGrid*600+125] = 4.0E0*I_ESP_H2xy2z_F2xz_aa-2.0E0*1*I_ESP_Fy2z_F2xz_a;
    abcd[iGrid*600+126] = 4.0E0*I_ESP_Hx4y_F2xz_aa-2.0E0*3*I_ESP_Fx2y_F2xz_a;
    abcd[iGrid*600+127] = 4.0E0*I_ESP_Hx3yz_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+128] = 4.0E0*I_ESP_Hx2y2z_F2xz_aa-2.0E0*1*I_ESP_Fx2z_F2xz_a;
    abcd[iGrid*600+129] = 4.0E0*I_ESP_Hxy3z_F2xz_aa;
    abcd[iGrid*600+130] = 4.0E0*I_ESP_H4xy_Fx2y_aa-2.0E0*3*I_ESP_F2xy_Fx2y_a;
    abcd[iGrid*600+131] = 4.0E0*I_ESP_H3x2y_Fx2y_aa-2.0E0*1*I_ESP_F3x_Fx2y_a-2.0E0*2*I_ESP_Fx2y_Fx2y_a+2*1*I_ESP_Px_Fx2y;
    abcd[iGrid*600+132] = 4.0E0*I_ESP_H3xyz_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+133] = 4.0E0*I_ESP_H2x3y_Fx2y_aa-2.0E0*2*I_ESP_F2xy_Fx2y_a-2.0E0*1*I_ESP_F3y_Fx2y_a+2*I_ESP_Py_Fx2y;
    abcd[iGrid*600+134] = 4.0E0*I_ESP_H2x2yz_Fx2y_aa-2.0E0*1*I_ESP_F2xz_Fx2y_a-2.0E0*1*I_ESP_F2yz_Fx2y_a+1*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+135] = 4.0E0*I_ESP_H2xy2z_Fx2y_aa-2.0E0*1*I_ESP_Fy2z_Fx2y_a;
    abcd[iGrid*600+136] = 4.0E0*I_ESP_Hx4y_Fx2y_aa-2.0E0*3*I_ESP_Fx2y_Fx2y_a;
    abcd[iGrid*600+137] = 4.0E0*I_ESP_Hx3yz_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+138] = 4.0E0*I_ESP_Hx2y2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2z_Fx2y_a;
    abcd[iGrid*600+139] = 4.0E0*I_ESP_Hxy3z_Fx2y_aa;
    abcd[iGrid*600+140] = 4.0E0*I_ESP_H4xy_Fxyz_aa-2.0E0*3*I_ESP_F2xy_Fxyz_a;
    abcd[iGrid*600+141] = 4.0E0*I_ESP_H3x2y_Fxyz_aa-2.0E0*1*I_ESP_F3x_Fxyz_a-2.0E0*2*I_ESP_Fx2y_Fxyz_a+2*1*I_ESP_Px_Fxyz;
    abcd[iGrid*600+142] = 4.0E0*I_ESP_H3xyz_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+143] = 4.0E0*I_ESP_H2x3y_Fxyz_aa-2.0E0*2*I_ESP_F2xy_Fxyz_a-2.0E0*1*I_ESP_F3y_Fxyz_a+2*I_ESP_Py_Fxyz;
    abcd[iGrid*600+144] = 4.0E0*I_ESP_H2x2yz_Fxyz_aa-2.0E0*1*I_ESP_F2xz_Fxyz_a-2.0E0*1*I_ESP_F2yz_Fxyz_a+1*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+145] = 4.0E0*I_ESP_H2xy2z_Fxyz_aa-2.0E0*1*I_ESP_Fy2z_Fxyz_a;
    abcd[iGrid*600+146] = 4.0E0*I_ESP_Hx4y_Fxyz_aa-2.0E0*3*I_ESP_Fx2y_Fxyz_a;
    abcd[iGrid*600+147] = 4.0E0*I_ESP_Hx3yz_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+148] = 4.0E0*I_ESP_Hx2y2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2z_Fxyz_a;
    abcd[iGrid*600+149] = 4.0E0*I_ESP_Hxy3z_Fxyz_aa;
    abcd[iGrid*600+150] = 4.0E0*I_ESP_H4xy_Fx2z_aa-2.0E0*3*I_ESP_F2xy_Fx2z_a;
    abcd[iGrid*600+151] = 4.0E0*I_ESP_H3x2y_Fx2z_aa-2.0E0*1*I_ESP_F3x_Fx2z_a-2.0E0*2*I_ESP_Fx2y_Fx2z_a+2*1*I_ESP_Px_Fx2z;
    abcd[iGrid*600+152] = 4.0E0*I_ESP_H3xyz_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+153] = 4.0E0*I_ESP_H2x3y_Fx2z_aa-2.0E0*2*I_ESP_F2xy_Fx2z_a-2.0E0*1*I_ESP_F3y_Fx2z_a+2*I_ESP_Py_Fx2z;
    abcd[iGrid*600+154] = 4.0E0*I_ESP_H2x2yz_Fx2z_aa-2.0E0*1*I_ESP_F2xz_Fx2z_a-2.0E0*1*I_ESP_F2yz_Fx2z_a+1*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+155] = 4.0E0*I_ESP_H2xy2z_Fx2z_aa-2.0E0*1*I_ESP_Fy2z_Fx2z_a;
    abcd[iGrid*600+156] = 4.0E0*I_ESP_Hx4y_Fx2z_aa-2.0E0*3*I_ESP_Fx2y_Fx2z_a;
    abcd[iGrid*600+157] = 4.0E0*I_ESP_Hx3yz_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+158] = 4.0E0*I_ESP_Hx2y2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2z_Fx2z_a;
    abcd[iGrid*600+159] = 4.0E0*I_ESP_Hxy3z_Fx2z_aa;
    abcd[iGrid*600+160] = 4.0E0*I_ESP_H4xy_F3y_aa-2.0E0*3*I_ESP_F2xy_F3y_a;
    abcd[iGrid*600+161] = 4.0E0*I_ESP_H3x2y_F3y_aa-2.0E0*1*I_ESP_F3x_F3y_a-2.0E0*2*I_ESP_Fx2y_F3y_a+2*1*I_ESP_Px_F3y;
    abcd[iGrid*600+162] = 4.0E0*I_ESP_H3xyz_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+163] = 4.0E0*I_ESP_H2x3y_F3y_aa-2.0E0*2*I_ESP_F2xy_F3y_a-2.0E0*1*I_ESP_F3y_F3y_a+2*I_ESP_Py_F3y;
    abcd[iGrid*600+164] = 4.0E0*I_ESP_H2x2yz_F3y_aa-2.0E0*1*I_ESP_F2xz_F3y_a-2.0E0*1*I_ESP_F2yz_F3y_a+1*I_ESP_Pz_F3y;
    abcd[iGrid*600+165] = 4.0E0*I_ESP_H2xy2z_F3y_aa-2.0E0*1*I_ESP_Fy2z_F3y_a;
    abcd[iGrid*600+166] = 4.0E0*I_ESP_Hx4y_F3y_aa-2.0E0*3*I_ESP_Fx2y_F3y_a;
    abcd[iGrid*600+167] = 4.0E0*I_ESP_Hx3yz_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+168] = 4.0E0*I_ESP_Hx2y2z_F3y_aa-2.0E0*1*I_ESP_Fx2z_F3y_a;
    abcd[iGrid*600+169] = 4.0E0*I_ESP_Hxy3z_F3y_aa;
    abcd[iGrid*600+170] = 4.0E0*I_ESP_H4xy_F2yz_aa-2.0E0*3*I_ESP_F2xy_F2yz_a;
    abcd[iGrid*600+171] = 4.0E0*I_ESP_H3x2y_F2yz_aa-2.0E0*1*I_ESP_F3x_F2yz_a-2.0E0*2*I_ESP_Fx2y_F2yz_a+2*1*I_ESP_Px_F2yz;
    abcd[iGrid*600+172] = 4.0E0*I_ESP_H3xyz_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+173] = 4.0E0*I_ESP_H2x3y_F2yz_aa-2.0E0*2*I_ESP_F2xy_F2yz_a-2.0E0*1*I_ESP_F3y_F2yz_a+2*I_ESP_Py_F2yz;
    abcd[iGrid*600+174] = 4.0E0*I_ESP_H2x2yz_F2yz_aa-2.0E0*1*I_ESP_F2xz_F2yz_a-2.0E0*1*I_ESP_F2yz_F2yz_a+1*I_ESP_Pz_F2yz;
    abcd[iGrid*600+175] = 4.0E0*I_ESP_H2xy2z_F2yz_aa-2.0E0*1*I_ESP_Fy2z_F2yz_a;
    abcd[iGrid*600+176] = 4.0E0*I_ESP_Hx4y_F2yz_aa-2.0E0*3*I_ESP_Fx2y_F2yz_a;
    abcd[iGrid*600+177] = 4.0E0*I_ESP_Hx3yz_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+178] = 4.0E0*I_ESP_Hx2y2z_F2yz_aa-2.0E0*1*I_ESP_Fx2z_F2yz_a;
    abcd[iGrid*600+179] = 4.0E0*I_ESP_Hxy3z_F2yz_aa;
    abcd[iGrid*600+180] = 4.0E0*I_ESP_H4xy_Fy2z_aa-2.0E0*3*I_ESP_F2xy_Fy2z_a;
    abcd[iGrid*600+181] = 4.0E0*I_ESP_H3x2y_Fy2z_aa-2.0E0*1*I_ESP_F3x_Fy2z_a-2.0E0*2*I_ESP_Fx2y_Fy2z_a+2*1*I_ESP_Px_Fy2z;
    abcd[iGrid*600+182] = 4.0E0*I_ESP_H3xyz_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+183] = 4.0E0*I_ESP_H2x3y_Fy2z_aa-2.0E0*2*I_ESP_F2xy_Fy2z_a-2.0E0*1*I_ESP_F3y_Fy2z_a+2*I_ESP_Py_Fy2z;
    abcd[iGrid*600+184] = 4.0E0*I_ESP_H2x2yz_Fy2z_aa-2.0E0*1*I_ESP_F2xz_Fy2z_a-2.0E0*1*I_ESP_F2yz_Fy2z_a+1*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+185] = 4.0E0*I_ESP_H2xy2z_Fy2z_aa-2.0E0*1*I_ESP_Fy2z_Fy2z_a;
    abcd[iGrid*600+186] = 4.0E0*I_ESP_Hx4y_Fy2z_aa-2.0E0*3*I_ESP_Fx2y_Fy2z_a;
    abcd[iGrid*600+187] = 4.0E0*I_ESP_Hx3yz_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+188] = 4.0E0*I_ESP_Hx2y2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2z_Fy2z_a;
    abcd[iGrid*600+189] = 4.0E0*I_ESP_Hxy3z_Fy2z_aa;
    abcd[iGrid*600+190] = 4.0E0*I_ESP_H4xy_F3z_aa-2.0E0*3*I_ESP_F2xy_F3z_a;
    abcd[iGrid*600+191] = 4.0E0*I_ESP_H3x2y_F3z_aa-2.0E0*1*I_ESP_F3x_F3z_a-2.0E0*2*I_ESP_Fx2y_F3z_a+2*1*I_ESP_Px_F3z;
    abcd[iGrid*600+192] = 4.0E0*I_ESP_H3xyz_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+193] = 4.0E0*I_ESP_H2x3y_F3z_aa-2.0E0*2*I_ESP_F2xy_F3z_a-2.0E0*1*I_ESP_F3y_F3z_a+2*I_ESP_Py_F3z;
    abcd[iGrid*600+194] = 4.0E0*I_ESP_H2x2yz_F3z_aa-2.0E0*1*I_ESP_F2xz_F3z_a-2.0E0*1*I_ESP_F2yz_F3z_a+1*I_ESP_Pz_F3z;
    abcd[iGrid*600+195] = 4.0E0*I_ESP_H2xy2z_F3z_aa-2.0E0*1*I_ESP_Fy2z_F3z_a;
    abcd[iGrid*600+196] = 4.0E0*I_ESP_Hx4y_F3z_aa-2.0E0*3*I_ESP_Fx2y_F3z_a;
    abcd[iGrid*600+197] = 4.0E0*I_ESP_Hx3yz_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+198] = 4.0E0*I_ESP_Hx2y2z_F3z_aa-2.0E0*1*I_ESP_Fx2z_F3z_a;
    abcd[iGrid*600+199] = 4.0E0*I_ESP_Hxy3z_F3z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+200] = 4.0E0*I_ESP_H4xz_F3x_aa-2.0E0*3*I_ESP_F2xz_F3x_a;
    abcd[iGrid*600+201] = 4.0E0*I_ESP_H3xyz_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+202] = 4.0E0*I_ESP_H3x2z_F3x_aa-2.0E0*1*I_ESP_F3x_F3x_a-2.0E0*2*I_ESP_Fx2z_F3x_a+2*1*I_ESP_Px_F3x;
    abcd[iGrid*600+203] = 4.0E0*I_ESP_H2x2yz_F3x_aa-2.0E0*1*I_ESP_F2yz_F3x_a;
    abcd[iGrid*600+204] = 4.0E0*I_ESP_H2xy2z_F3x_aa-2.0E0*1*I_ESP_F2xy_F3x_a-2.0E0*1*I_ESP_Fy2z_F3x_a+1*I_ESP_Py_F3x;
    abcd[iGrid*600+205] = 4.0E0*I_ESP_H2x3z_F3x_aa-2.0E0*2*I_ESP_F2xz_F3x_a-2.0E0*1*I_ESP_F3z_F3x_a+2*I_ESP_Pz_F3x;
    abcd[iGrid*600+206] = 4.0E0*I_ESP_Hx3yz_F3x_aa;
    abcd[iGrid*600+207] = 4.0E0*I_ESP_Hx2y2z_F3x_aa-2.0E0*1*I_ESP_Fx2y_F3x_a;
    abcd[iGrid*600+208] = 4.0E0*I_ESP_Hxy3z_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+209] = 4.0E0*I_ESP_Hx4z_F3x_aa-2.0E0*3*I_ESP_Fx2z_F3x_a;
    abcd[iGrid*600+210] = 4.0E0*I_ESP_H4xz_F2xy_aa-2.0E0*3*I_ESP_F2xz_F2xy_a;
    abcd[iGrid*600+211] = 4.0E0*I_ESP_H3xyz_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+212] = 4.0E0*I_ESP_H3x2z_F2xy_aa-2.0E0*1*I_ESP_F3x_F2xy_a-2.0E0*2*I_ESP_Fx2z_F2xy_a+2*1*I_ESP_Px_F2xy;
    abcd[iGrid*600+213] = 4.0E0*I_ESP_H2x2yz_F2xy_aa-2.0E0*1*I_ESP_F2yz_F2xy_a;
    abcd[iGrid*600+214] = 4.0E0*I_ESP_H2xy2z_F2xy_aa-2.0E0*1*I_ESP_F2xy_F2xy_a-2.0E0*1*I_ESP_Fy2z_F2xy_a+1*I_ESP_Py_F2xy;
    abcd[iGrid*600+215] = 4.0E0*I_ESP_H2x3z_F2xy_aa-2.0E0*2*I_ESP_F2xz_F2xy_a-2.0E0*1*I_ESP_F3z_F2xy_a+2*I_ESP_Pz_F2xy;
    abcd[iGrid*600+216] = 4.0E0*I_ESP_Hx3yz_F2xy_aa;
    abcd[iGrid*600+217] = 4.0E0*I_ESP_Hx2y2z_F2xy_aa-2.0E0*1*I_ESP_Fx2y_F2xy_a;
    abcd[iGrid*600+218] = 4.0E0*I_ESP_Hxy3z_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+219] = 4.0E0*I_ESP_Hx4z_F2xy_aa-2.0E0*3*I_ESP_Fx2z_F2xy_a;
    abcd[iGrid*600+220] = 4.0E0*I_ESP_H4xz_F2xz_aa-2.0E0*3*I_ESP_F2xz_F2xz_a;
    abcd[iGrid*600+221] = 4.0E0*I_ESP_H3xyz_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+222] = 4.0E0*I_ESP_H3x2z_F2xz_aa-2.0E0*1*I_ESP_F3x_F2xz_a-2.0E0*2*I_ESP_Fx2z_F2xz_a+2*1*I_ESP_Px_F2xz;
    abcd[iGrid*600+223] = 4.0E0*I_ESP_H2x2yz_F2xz_aa-2.0E0*1*I_ESP_F2yz_F2xz_a;
    abcd[iGrid*600+224] = 4.0E0*I_ESP_H2xy2z_F2xz_aa-2.0E0*1*I_ESP_F2xy_F2xz_a-2.0E0*1*I_ESP_Fy2z_F2xz_a+1*I_ESP_Py_F2xz;
    abcd[iGrid*600+225] = 4.0E0*I_ESP_H2x3z_F2xz_aa-2.0E0*2*I_ESP_F2xz_F2xz_a-2.0E0*1*I_ESP_F3z_F2xz_a+2*I_ESP_Pz_F2xz;
    abcd[iGrid*600+226] = 4.0E0*I_ESP_Hx3yz_F2xz_aa;
    abcd[iGrid*600+227] = 4.0E0*I_ESP_Hx2y2z_F2xz_aa-2.0E0*1*I_ESP_Fx2y_F2xz_a;
    abcd[iGrid*600+228] = 4.0E0*I_ESP_Hxy3z_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+229] = 4.0E0*I_ESP_Hx4z_F2xz_aa-2.0E0*3*I_ESP_Fx2z_F2xz_a;
    abcd[iGrid*600+230] = 4.0E0*I_ESP_H4xz_Fx2y_aa-2.0E0*3*I_ESP_F2xz_Fx2y_a;
    abcd[iGrid*600+231] = 4.0E0*I_ESP_H3xyz_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+232] = 4.0E0*I_ESP_H3x2z_Fx2y_aa-2.0E0*1*I_ESP_F3x_Fx2y_a-2.0E0*2*I_ESP_Fx2z_Fx2y_a+2*1*I_ESP_Px_Fx2y;
    abcd[iGrid*600+233] = 4.0E0*I_ESP_H2x2yz_Fx2y_aa-2.0E0*1*I_ESP_F2yz_Fx2y_a;
    abcd[iGrid*600+234] = 4.0E0*I_ESP_H2xy2z_Fx2y_aa-2.0E0*1*I_ESP_F2xy_Fx2y_a-2.0E0*1*I_ESP_Fy2z_Fx2y_a+1*I_ESP_Py_Fx2y;
    abcd[iGrid*600+235] = 4.0E0*I_ESP_H2x3z_Fx2y_aa-2.0E0*2*I_ESP_F2xz_Fx2y_a-2.0E0*1*I_ESP_F3z_Fx2y_a+2*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+236] = 4.0E0*I_ESP_Hx3yz_Fx2y_aa;
    abcd[iGrid*600+237] = 4.0E0*I_ESP_Hx2y2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2y_Fx2y_a;
    abcd[iGrid*600+238] = 4.0E0*I_ESP_Hxy3z_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+239] = 4.0E0*I_ESP_Hx4z_Fx2y_aa-2.0E0*3*I_ESP_Fx2z_Fx2y_a;
    abcd[iGrid*600+240] = 4.0E0*I_ESP_H4xz_Fxyz_aa-2.0E0*3*I_ESP_F2xz_Fxyz_a;
    abcd[iGrid*600+241] = 4.0E0*I_ESP_H3xyz_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+242] = 4.0E0*I_ESP_H3x2z_Fxyz_aa-2.0E0*1*I_ESP_F3x_Fxyz_a-2.0E0*2*I_ESP_Fx2z_Fxyz_a+2*1*I_ESP_Px_Fxyz;
    abcd[iGrid*600+243] = 4.0E0*I_ESP_H2x2yz_Fxyz_aa-2.0E0*1*I_ESP_F2yz_Fxyz_a;
    abcd[iGrid*600+244] = 4.0E0*I_ESP_H2xy2z_Fxyz_aa-2.0E0*1*I_ESP_F2xy_Fxyz_a-2.0E0*1*I_ESP_Fy2z_Fxyz_a+1*I_ESP_Py_Fxyz;
    abcd[iGrid*600+245] = 4.0E0*I_ESP_H2x3z_Fxyz_aa-2.0E0*2*I_ESP_F2xz_Fxyz_a-2.0E0*1*I_ESP_F3z_Fxyz_a+2*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+246] = 4.0E0*I_ESP_Hx3yz_Fxyz_aa;
    abcd[iGrid*600+247] = 4.0E0*I_ESP_Hx2y2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2y_Fxyz_a;
    abcd[iGrid*600+248] = 4.0E0*I_ESP_Hxy3z_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+249] = 4.0E0*I_ESP_Hx4z_Fxyz_aa-2.0E0*3*I_ESP_Fx2z_Fxyz_a;
    abcd[iGrid*600+250] = 4.0E0*I_ESP_H4xz_Fx2z_aa-2.0E0*3*I_ESP_F2xz_Fx2z_a;
    abcd[iGrid*600+251] = 4.0E0*I_ESP_H3xyz_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+252] = 4.0E0*I_ESP_H3x2z_Fx2z_aa-2.0E0*1*I_ESP_F3x_Fx2z_a-2.0E0*2*I_ESP_Fx2z_Fx2z_a+2*1*I_ESP_Px_Fx2z;
    abcd[iGrid*600+253] = 4.0E0*I_ESP_H2x2yz_Fx2z_aa-2.0E0*1*I_ESP_F2yz_Fx2z_a;
    abcd[iGrid*600+254] = 4.0E0*I_ESP_H2xy2z_Fx2z_aa-2.0E0*1*I_ESP_F2xy_Fx2z_a-2.0E0*1*I_ESP_Fy2z_Fx2z_a+1*I_ESP_Py_Fx2z;
    abcd[iGrid*600+255] = 4.0E0*I_ESP_H2x3z_Fx2z_aa-2.0E0*2*I_ESP_F2xz_Fx2z_a-2.0E0*1*I_ESP_F3z_Fx2z_a+2*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+256] = 4.0E0*I_ESP_Hx3yz_Fx2z_aa;
    abcd[iGrid*600+257] = 4.0E0*I_ESP_Hx2y2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2y_Fx2z_a;
    abcd[iGrid*600+258] = 4.0E0*I_ESP_Hxy3z_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+259] = 4.0E0*I_ESP_Hx4z_Fx2z_aa-2.0E0*3*I_ESP_Fx2z_Fx2z_a;
    abcd[iGrid*600+260] = 4.0E0*I_ESP_H4xz_F3y_aa-2.0E0*3*I_ESP_F2xz_F3y_a;
    abcd[iGrid*600+261] = 4.0E0*I_ESP_H3xyz_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+262] = 4.0E0*I_ESP_H3x2z_F3y_aa-2.0E0*1*I_ESP_F3x_F3y_a-2.0E0*2*I_ESP_Fx2z_F3y_a+2*1*I_ESP_Px_F3y;
    abcd[iGrid*600+263] = 4.0E0*I_ESP_H2x2yz_F3y_aa-2.0E0*1*I_ESP_F2yz_F3y_a;
    abcd[iGrid*600+264] = 4.0E0*I_ESP_H2xy2z_F3y_aa-2.0E0*1*I_ESP_F2xy_F3y_a-2.0E0*1*I_ESP_Fy2z_F3y_a+1*I_ESP_Py_F3y;
    abcd[iGrid*600+265] = 4.0E0*I_ESP_H2x3z_F3y_aa-2.0E0*2*I_ESP_F2xz_F3y_a-2.0E0*1*I_ESP_F3z_F3y_a+2*I_ESP_Pz_F3y;
    abcd[iGrid*600+266] = 4.0E0*I_ESP_Hx3yz_F3y_aa;
    abcd[iGrid*600+267] = 4.0E0*I_ESP_Hx2y2z_F3y_aa-2.0E0*1*I_ESP_Fx2y_F3y_a;
    abcd[iGrid*600+268] = 4.0E0*I_ESP_Hxy3z_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+269] = 4.0E0*I_ESP_Hx4z_F3y_aa-2.0E0*3*I_ESP_Fx2z_F3y_a;
    abcd[iGrid*600+270] = 4.0E0*I_ESP_H4xz_F2yz_aa-2.0E0*3*I_ESP_F2xz_F2yz_a;
    abcd[iGrid*600+271] = 4.0E0*I_ESP_H3xyz_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+272] = 4.0E0*I_ESP_H3x2z_F2yz_aa-2.0E0*1*I_ESP_F3x_F2yz_a-2.0E0*2*I_ESP_Fx2z_F2yz_a+2*1*I_ESP_Px_F2yz;
    abcd[iGrid*600+273] = 4.0E0*I_ESP_H2x2yz_F2yz_aa-2.0E0*1*I_ESP_F2yz_F2yz_a;
    abcd[iGrid*600+274] = 4.0E0*I_ESP_H2xy2z_F2yz_aa-2.0E0*1*I_ESP_F2xy_F2yz_a-2.0E0*1*I_ESP_Fy2z_F2yz_a+1*I_ESP_Py_F2yz;
    abcd[iGrid*600+275] = 4.0E0*I_ESP_H2x3z_F2yz_aa-2.0E0*2*I_ESP_F2xz_F2yz_a-2.0E0*1*I_ESP_F3z_F2yz_a+2*I_ESP_Pz_F2yz;
    abcd[iGrid*600+276] = 4.0E0*I_ESP_Hx3yz_F2yz_aa;
    abcd[iGrid*600+277] = 4.0E0*I_ESP_Hx2y2z_F2yz_aa-2.0E0*1*I_ESP_Fx2y_F2yz_a;
    abcd[iGrid*600+278] = 4.0E0*I_ESP_Hxy3z_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+279] = 4.0E0*I_ESP_Hx4z_F2yz_aa-2.0E0*3*I_ESP_Fx2z_F2yz_a;
    abcd[iGrid*600+280] = 4.0E0*I_ESP_H4xz_Fy2z_aa-2.0E0*3*I_ESP_F2xz_Fy2z_a;
    abcd[iGrid*600+281] = 4.0E0*I_ESP_H3xyz_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+282] = 4.0E0*I_ESP_H3x2z_Fy2z_aa-2.0E0*1*I_ESP_F3x_Fy2z_a-2.0E0*2*I_ESP_Fx2z_Fy2z_a+2*1*I_ESP_Px_Fy2z;
    abcd[iGrid*600+283] = 4.0E0*I_ESP_H2x2yz_Fy2z_aa-2.0E0*1*I_ESP_F2yz_Fy2z_a;
    abcd[iGrid*600+284] = 4.0E0*I_ESP_H2xy2z_Fy2z_aa-2.0E0*1*I_ESP_F2xy_Fy2z_a-2.0E0*1*I_ESP_Fy2z_Fy2z_a+1*I_ESP_Py_Fy2z;
    abcd[iGrid*600+285] = 4.0E0*I_ESP_H2x3z_Fy2z_aa-2.0E0*2*I_ESP_F2xz_Fy2z_a-2.0E0*1*I_ESP_F3z_Fy2z_a+2*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+286] = 4.0E0*I_ESP_Hx3yz_Fy2z_aa;
    abcd[iGrid*600+287] = 4.0E0*I_ESP_Hx2y2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2y_Fy2z_a;
    abcd[iGrid*600+288] = 4.0E0*I_ESP_Hxy3z_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+289] = 4.0E0*I_ESP_Hx4z_Fy2z_aa-2.0E0*3*I_ESP_Fx2z_Fy2z_a;
    abcd[iGrid*600+290] = 4.0E0*I_ESP_H4xz_F3z_aa-2.0E0*3*I_ESP_F2xz_F3z_a;
    abcd[iGrid*600+291] = 4.0E0*I_ESP_H3xyz_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+292] = 4.0E0*I_ESP_H3x2z_F3z_aa-2.0E0*1*I_ESP_F3x_F3z_a-2.0E0*2*I_ESP_Fx2z_F3z_a+2*1*I_ESP_Px_F3z;
    abcd[iGrid*600+293] = 4.0E0*I_ESP_H2x2yz_F3z_aa-2.0E0*1*I_ESP_F2yz_F3z_a;
    abcd[iGrid*600+294] = 4.0E0*I_ESP_H2xy2z_F3z_aa-2.0E0*1*I_ESP_F2xy_F3z_a-2.0E0*1*I_ESP_Fy2z_F3z_a+1*I_ESP_Py_F3z;
    abcd[iGrid*600+295] = 4.0E0*I_ESP_H2x3z_F3z_aa-2.0E0*2*I_ESP_F2xz_F3z_a-2.0E0*1*I_ESP_F3z_F3z_a+2*I_ESP_Pz_F3z;
    abcd[iGrid*600+296] = 4.0E0*I_ESP_Hx3yz_F3z_aa;
    abcd[iGrid*600+297] = 4.0E0*I_ESP_Hx2y2z_F3z_aa-2.0E0*1*I_ESP_Fx2y_F3z_a;
    abcd[iGrid*600+298] = 4.0E0*I_ESP_Hxy3z_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+299] = 4.0E0*I_ESP_Hx4z_F3z_aa-2.0E0*3*I_ESP_Fx2z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+300] = 4.0E0*I_ESP_H3x2y_F3x_aa-2.0E0*1*I_ESP_F3x_F3x_a;
    abcd[iGrid*600+301] = 4.0E0*I_ESP_H2x3y_F3x_aa-2.0E0*1*I_ESP_F2xy_F3x_a-2.0E0*2*I_ESP_F2xy_F3x_a;
    abcd[iGrid*600+302] = 4.0E0*I_ESP_H2x2yz_F3x_aa-2.0E0*1*I_ESP_F2xz_F3x_a;
    abcd[iGrid*600+303] = 4.0E0*I_ESP_Hx4y_F3x_aa-2.0E0*2*I_ESP_Fx2y_F3x_a-2.0E0*3*I_ESP_Fx2y_F3x_a+2*1*I_ESP_Px_F3x;
    abcd[iGrid*600+304] = 4.0E0*I_ESP_Hx3yz_F3x_aa-2.0E0*1*I_ESP_Fxyz_F3x_a-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+305] = 4.0E0*I_ESP_Hx2y2z_F3x_aa-2.0E0*1*I_ESP_Fx2z_F3x_a;
    abcd[iGrid*600+306] = 4.0E0*I_ESP_H5y_F3x_aa-2.0E0*3*I_ESP_F3y_F3x_a-2.0E0*4*I_ESP_F3y_F3x_a+3*2*I_ESP_Py_F3x;
    abcd[iGrid*600+307] = 4.0E0*I_ESP_H4yz_F3x_aa-2.0E0*2*I_ESP_F2yz_F3x_a-2.0E0*3*I_ESP_F2yz_F3x_a+2*1*I_ESP_Pz_F3x;
    abcd[iGrid*600+308] = 4.0E0*I_ESP_H3y2z_F3x_aa-2.0E0*1*I_ESP_Fy2z_F3x_a-2.0E0*2*I_ESP_Fy2z_F3x_a;
    abcd[iGrid*600+309] = 4.0E0*I_ESP_H2y3z_F3x_aa-2.0E0*1*I_ESP_F3z_F3x_a;
    abcd[iGrid*600+310] = 4.0E0*I_ESP_H3x2y_F2xy_aa-2.0E0*1*I_ESP_F3x_F2xy_a;
    abcd[iGrid*600+311] = 4.0E0*I_ESP_H2x3y_F2xy_aa-2.0E0*1*I_ESP_F2xy_F2xy_a-2.0E0*2*I_ESP_F2xy_F2xy_a;
    abcd[iGrid*600+312] = 4.0E0*I_ESP_H2x2yz_F2xy_aa-2.0E0*1*I_ESP_F2xz_F2xy_a;
    abcd[iGrid*600+313] = 4.0E0*I_ESP_Hx4y_F2xy_aa-2.0E0*2*I_ESP_Fx2y_F2xy_a-2.0E0*3*I_ESP_Fx2y_F2xy_a+2*1*I_ESP_Px_F2xy;
    abcd[iGrid*600+314] = 4.0E0*I_ESP_Hx3yz_F2xy_aa-2.0E0*1*I_ESP_Fxyz_F2xy_a-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+315] = 4.0E0*I_ESP_Hx2y2z_F2xy_aa-2.0E0*1*I_ESP_Fx2z_F2xy_a;
    abcd[iGrid*600+316] = 4.0E0*I_ESP_H5y_F2xy_aa-2.0E0*3*I_ESP_F3y_F2xy_a-2.0E0*4*I_ESP_F3y_F2xy_a+3*2*I_ESP_Py_F2xy;
    abcd[iGrid*600+317] = 4.0E0*I_ESP_H4yz_F2xy_aa-2.0E0*2*I_ESP_F2yz_F2xy_a-2.0E0*3*I_ESP_F2yz_F2xy_a+2*1*I_ESP_Pz_F2xy;
    abcd[iGrid*600+318] = 4.0E0*I_ESP_H3y2z_F2xy_aa-2.0E0*1*I_ESP_Fy2z_F2xy_a-2.0E0*2*I_ESP_Fy2z_F2xy_a;
    abcd[iGrid*600+319] = 4.0E0*I_ESP_H2y3z_F2xy_aa-2.0E0*1*I_ESP_F3z_F2xy_a;
    abcd[iGrid*600+320] = 4.0E0*I_ESP_H3x2y_F2xz_aa-2.0E0*1*I_ESP_F3x_F2xz_a;
    abcd[iGrid*600+321] = 4.0E0*I_ESP_H2x3y_F2xz_aa-2.0E0*1*I_ESP_F2xy_F2xz_a-2.0E0*2*I_ESP_F2xy_F2xz_a;
    abcd[iGrid*600+322] = 4.0E0*I_ESP_H2x2yz_F2xz_aa-2.0E0*1*I_ESP_F2xz_F2xz_a;
    abcd[iGrid*600+323] = 4.0E0*I_ESP_Hx4y_F2xz_aa-2.0E0*2*I_ESP_Fx2y_F2xz_a-2.0E0*3*I_ESP_Fx2y_F2xz_a+2*1*I_ESP_Px_F2xz;
    abcd[iGrid*600+324] = 4.0E0*I_ESP_Hx3yz_F2xz_aa-2.0E0*1*I_ESP_Fxyz_F2xz_a-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+325] = 4.0E0*I_ESP_Hx2y2z_F2xz_aa-2.0E0*1*I_ESP_Fx2z_F2xz_a;
    abcd[iGrid*600+326] = 4.0E0*I_ESP_H5y_F2xz_aa-2.0E0*3*I_ESP_F3y_F2xz_a-2.0E0*4*I_ESP_F3y_F2xz_a+3*2*I_ESP_Py_F2xz;
    abcd[iGrid*600+327] = 4.0E0*I_ESP_H4yz_F2xz_aa-2.0E0*2*I_ESP_F2yz_F2xz_a-2.0E0*3*I_ESP_F2yz_F2xz_a+2*1*I_ESP_Pz_F2xz;
    abcd[iGrid*600+328] = 4.0E0*I_ESP_H3y2z_F2xz_aa-2.0E0*1*I_ESP_Fy2z_F2xz_a-2.0E0*2*I_ESP_Fy2z_F2xz_a;
    abcd[iGrid*600+329] = 4.0E0*I_ESP_H2y3z_F2xz_aa-2.0E0*1*I_ESP_F3z_F2xz_a;
    abcd[iGrid*600+330] = 4.0E0*I_ESP_H3x2y_Fx2y_aa-2.0E0*1*I_ESP_F3x_Fx2y_a;
    abcd[iGrid*600+331] = 4.0E0*I_ESP_H2x3y_Fx2y_aa-2.0E0*1*I_ESP_F2xy_Fx2y_a-2.0E0*2*I_ESP_F2xy_Fx2y_a;
    abcd[iGrid*600+332] = 4.0E0*I_ESP_H2x2yz_Fx2y_aa-2.0E0*1*I_ESP_F2xz_Fx2y_a;
    abcd[iGrid*600+333] = 4.0E0*I_ESP_Hx4y_Fx2y_aa-2.0E0*2*I_ESP_Fx2y_Fx2y_a-2.0E0*3*I_ESP_Fx2y_Fx2y_a+2*1*I_ESP_Px_Fx2y;
    abcd[iGrid*600+334] = 4.0E0*I_ESP_Hx3yz_Fx2y_aa-2.0E0*1*I_ESP_Fxyz_Fx2y_a-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+335] = 4.0E0*I_ESP_Hx2y2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2z_Fx2y_a;
    abcd[iGrid*600+336] = 4.0E0*I_ESP_H5y_Fx2y_aa-2.0E0*3*I_ESP_F3y_Fx2y_a-2.0E0*4*I_ESP_F3y_Fx2y_a+3*2*I_ESP_Py_Fx2y;
    abcd[iGrid*600+337] = 4.0E0*I_ESP_H4yz_Fx2y_aa-2.0E0*2*I_ESP_F2yz_Fx2y_a-2.0E0*3*I_ESP_F2yz_Fx2y_a+2*1*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+338] = 4.0E0*I_ESP_H3y2z_Fx2y_aa-2.0E0*1*I_ESP_Fy2z_Fx2y_a-2.0E0*2*I_ESP_Fy2z_Fx2y_a;
    abcd[iGrid*600+339] = 4.0E0*I_ESP_H2y3z_Fx2y_aa-2.0E0*1*I_ESP_F3z_Fx2y_a;
    abcd[iGrid*600+340] = 4.0E0*I_ESP_H3x2y_Fxyz_aa-2.0E0*1*I_ESP_F3x_Fxyz_a;
    abcd[iGrid*600+341] = 4.0E0*I_ESP_H2x3y_Fxyz_aa-2.0E0*1*I_ESP_F2xy_Fxyz_a-2.0E0*2*I_ESP_F2xy_Fxyz_a;
    abcd[iGrid*600+342] = 4.0E0*I_ESP_H2x2yz_Fxyz_aa-2.0E0*1*I_ESP_F2xz_Fxyz_a;
    abcd[iGrid*600+343] = 4.0E0*I_ESP_Hx4y_Fxyz_aa-2.0E0*2*I_ESP_Fx2y_Fxyz_a-2.0E0*3*I_ESP_Fx2y_Fxyz_a+2*1*I_ESP_Px_Fxyz;
    abcd[iGrid*600+344] = 4.0E0*I_ESP_Hx3yz_Fxyz_aa-2.0E0*1*I_ESP_Fxyz_Fxyz_a-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+345] = 4.0E0*I_ESP_Hx2y2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2z_Fxyz_a;
    abcd[iGrid*600+346] = 4.0E0*I_ESP_H5y_Fxyz_aa-2.0E0*3*I_ESP_F3y_Fxyz_a-2.0E0*4*I_ESP_F3y_Fxyz_a+3*2*I_ESP_Py_Fxyz;
    abcd[iGrid*600+347] = 4.0E0*I_ESP_H4yz_Fxyz_aa-2.0E0*2*I_ESP_F2yz_Fxyz_a-2.0E0*3*I_ESP_F2yz_Fxyz_a+2*1*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+348] = 4.0E0*I_ESP_H3y2z_Fxyz_aa-2.0E0*1*I_ESP_Fy2z_Fxyz_a-2.0E0*2*I_ESP_Fy2z_Fxyz_a;
    abcd[iGrid*600+349] = 4.0E0*I_ESP_H2y3z_Fxyz_aa-2.0E0*1*I_ESP_F3z_Fxyz_a;
    abcd[iGrid*600+350] = 4.0E0*I_ESP_H3x2y_Fx2z_aa-2.0E0*1*I_ESP_F3x_Fx2z_a;
    abcd[iGrid*600+351] = 4.0E0*I_ESP_H2x3y_Fx2z_aa-2.0E0*1*I_ESP_F2xy_Fx2z_a-2.0E0*2*I_ESP_F2xy_Fx2z_a;
    abcd[iGrid*600+352] = 4.0E0*I_ESP_H2x2yz_Fx2z_aa-2.0E0*1*I_ESP_F2xz_Fx2z_a;
    abcd[iGrid*600+353] = 4.0E0*I_ESP_Hx4y_Fx2z_aa-2.0E0*2*I_ESP_Fx2y_Fx2z_a-2.0E0*3*I_ESP_Fx2y_Fx2z_a+2*1*I_ESP_Px_Fx2z;
    abcd[iGrid*600+354] = 4.0E0*I_ESP_Hx3yz_Fx2z_aa-2.0E0*1*I_ESP_Fxyz_Fx2z_a-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+355] = 4.0E0*I_ESP_Hx2y2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2z_Fx2z_a;
    abcd[iGrid*600+356] = 4.0E0*I_ESP_H5y_Fx2z_aa-2.0E0*3*I_ESP_F3y_Fx2z_a-2.0E0*4*I_ESP_F3y_Fx2z_a+3*2*I_ESP_Py_Fx2z;
    abcd[iGrid*600+357] = 4.0E0*I_ESP_H4yz_Fx2z_aa-2.0E0*2*I_ESP_F2yz_Fx2z_a-2.0E0*3*I_ESP_F2yz_Fx2z_a+2*1*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+358] = 4.0E0*I_ESP_H3y2z_Fx2z_aa-2.0E0*1*I_ESP_Fy2z_Fx2z_a-2.0E0*2*I_ESP_Fy2z_Fx2z_a;
    abcd[iGrid*600+359] = 4.0E0*I_ESP_H2y3z_Fx2z_aa-2.0E0*1*I_ESP_F3z_Fx2z_a;
    abcd[iGrid*600+360] = 4.0E0*I_ESP_H3x2y_F3y_aa-2.0E0*1*I_ESP_F3x_F3y_a;
    abcd[iGrid*600+361] = 4.0E0*I_ESP_H2x3y_F3y_aa-2.0E0*1*I_ESP_F2xy_F3y_a-2.0E0*2*I_ESP_F2xy_F3y_a;
    abcd[iGrid*600+362] = 4.0E0*I_ESP_H2x2yz_F3y_aa-2.0E0*1*I_ESP_F2xz_F3y_a;
    abcd[iGrid*600+363] = 4.0E0*I_ESP_Hx4y_F3y_aa-2.0E0*2*I_ESP_Fx2y_F3y_a-2.0E0*3*I_ESP_Fx2y_F3y_a+2*1*I_ESP_Px_F3y;
    abcd[iGrid*600+364] = 4.0E0*I_ESP_Hx3yz_F3y_aa-2.0E0*1*I_ESP_Fxyz_F3y_a-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+365] = 4.0E0*I_ESP_Hx2y2z_F3y_aa-2.0E0*1*I_ESP_Fx2z_F3y_a;
    abcd[iGrid*600+366] = 4.0E0*I_ESP_H5y_F3y_aa-2.0E0*3*I_ESP_F3y_F3y_a-2.0E0*4*I_ESP_F3y_F3y_a+3*2*I_ESP_Py_F3y;
    abcd[iGrid*600+367] = 4.0E0*I_ESP_H4yz_F3y_aa-2.0E0*2*I_ESP_F2yz_F3y_a-2.0E0*3*I_ESP_F2yz_F3y_a+2*1*I_ESP_Pz_F3y;
    abcd[iGrid*600+368] = 4.0E0*I_ESP_H3y2z_F3y_aa-2.0E0*1*I_ESP_Fy2z_F3y_a-2.0E0*2*I_ESP_Fy2z_F3y_a;
    abcd[iGrid*600+369] = 4.0E0*I_ESP_H2y3z_F3y_aa-2.0E0*1*I_ESP_F3z_F3y_a;
    abcd[iGrid*600+370] = 4.0E0*I_ESP_H3x2y_F2yz_aa-2.0E0*1*I_ESP_F3x_F2yz_a;
    abcd[iGrid*600+371] = 4.0E0*I_ESP_H2x3y_F2yz_aa-2.0E0*1*I_ESP_F2xy_F2yz_a-2.0E0*2*I_ESP_F2xy_F2yz_a;
    abcd[iGrid*600+372] = 4.0E0*I_ESP_H2x2yz_F2yz_aa-2.0E0*1*I_ESP_F2xz_F2yz_a;
    abcd[iGrid*600+373] = 4.0E0*I_ESP_Hx4y_F2yz_aa-2.0E0*2*I_ESP_Fx2y_F2yz_a-2.0E0*3*I_ESP_Fx2y_F2yz_a+2*1*I_ESP_Px_F2yz;
    abcd[iGrid*600+374] = 4.0E0*I_ESP_Hx3yz_F2yz_aa-2.0E0*1*I_ESP_Fxyz_F2yz_a-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+375] = 4.0E0*I_ESP_Hx2y2z_F2yz_aa-2.0E0*1*I_ESP_Fx2z_F2yz_a;
    abcd[iGrid*600+376] = 4.0E0*I_ESP_H5y_F2yz_aa-2.0E0*3*I_ESP_F3y_F2yz_a-2.0E0*4*I_ESP_F3y_F2yz_a+3*2*I_ESP_Py_F2yz;
    abcd[iGrid*600+377] = 4.0E0*I_ESP_H4yz_F2yz_aa-2.0E0*2*I_ESP_F2yz_F2yz_a-2.0E0*3*I_ESP_F2yz_F2yz_a+2*1*I_ESP_Pz_F2yz;
    abcd[iGrid*600+378] = 4.0E0*I_ESP_H3y2z_F2yz_aa-2.0E0*1*I_ESP_Fy2z_F2yz_a-2.0E0*2*I_ESP_Fy2z_F2yz_a;
    abcd[iGrid*600+379] = 4.0E0*I_ESP_H2y3z_F2yz_aa-2.0E0*1*I_ESP_F3z_F2yz_a;
    abcd[iGrid*600+380] = 4.0E0*I_ESP_H3x2y_Fy2z_aa-2.0E0*1*I_ESP_F3x_Fy2z_a;
    abcd[iGrid*600+381] = 4.0E0*I_ESP_H2x3y_Fy2z_aa-2.0E0*1*I_ESP_F2xy_Fy2z_a-2.0E0*2*I_ESP_F2xy_Fy2z_a;
    abcd[iGrid*600+382] = 4.0E0*I_ESP_H2x2yz_Fy2z_aa-2.0E0*1*I_ESP_F2xz_Fy2z_a;
    abcd[iGrid*600+383] = 4.0E0*I_ESP_Hx4y_Fy2z_aa-2.0E0*2*I_ESP_Fx2y_Fy2z_a-2.0E0*3*I_ESP_Fx2y_Fy2z_a+2*1*I_ESP_Px_Fy2z;
    abcd[iGrid*600+384] = 4.0E0*I_ESP_Hx3yz_Fy2z_aa-2.0E0*1*I_ESP_Fxyz_Fy2z_a-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+385] = 4.0E0*I_ESP_Hx2y2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2z_Fy2z_a;
    abcd[iGrid*600+386] = 4.0E0*I_ESP_H5y_Fy2z_aa-2.0E0*3*I_ESP_F3y_Fy2z_a-2.0E0*4*I_ESP_F3y_Fy2z_a+3*2*I_ESP_Py_Fy2z;
    abcd[iGrid*600+387] = 4.0E0*I_ESP_H4yz_Fy2z_aa-2.0E0*2*I_ESP_F2yz_Fy2z_a-2.0E0*3*I_ESP_F2yz_Fy2z_a+2*1*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+388] = 4.0E0*I_ESP_H3y2z_Fy2z_aa-2.0E0*1*I_ESP_Fy2z_Fy2z_a-2.0E0*2*I_ESP_Fy2z_Fy2z_a;
    abcd[iGrid*600+389] = 4.0E0*I_ESP_H2y3z_Fy2z_aa-2.0E0*1*I_ESP_F3z_Fy2z_a;
    abcd[iGrid*600+390] = 4.0E0*I_ESP_H3x2y_F3z_aa-2.0E0*1*I_ESP_F3x_F3z_a;
    abcd[iGrid*600+391] = 4.0E0*I_ESP_H2x3y_F3z_aa-2.0E0*1*I_ESP_F2xy_F3z_a-2.0E0*2*I_ESP_F2xy_F3z_a;
    abcd[iGrid*600+392] = 4.0E0*I_ESP_H2x2yz_F3z_aa-2.0E0*1*I_ESP_F2xz_F3z_a;
    abcd[iGrid*600+393] = 4.0E0*I_ESP_Hx4y_F3z_aa-2.0E0*2*I_ESP_Fx2y_F3z_a-2.0E0*3*I_ESP_Fx2y_F3z_a+2*1*I_ESP_Px_F3z;
    abcd[iGrid*600+394] = 4.0E0*I_ESP_Hx3yz_F3z_aa-2.0E0*1*I_ESP_Fxyz_F3z_a-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+395] = 4.0E0*I_ESP_Hx2y2z_F3z_aa-2.0E0*1*I_ESP_Fx2z_F3z_a;
    abcd[iGrid*600+396] = 4.0E0*I_ESP_H5y_F3z_aa-2.0E0*3*I_ESP_F3y_F3z_a-2.0E0*4*I_ESP_F3y_F3z_a+3*2*I_ESP_Py_F3z;
    abcd[iGrid*600+397] = 4.0E0*I_ESP_H4yz_F3z_aa-2.0E0*2*I_ESP_F2yz_F3z_a-2.0E0*3*I_ESP_F2yz_F3z_a+2*1*I_ESP_Pz_F3z;
    abcd[iGrid*600+398] = 4.0E0*I_ESP_H3y2z_F3z_aa-2.0E0*1*I_ESP_Fy2z_F3z_a-2.0E0*2*I_ESP_Fy2z_F3z_a;
    abcd[iGrid*600+399] = 4.0E0*I_ESP_H2y3z_F3z_aa-2.0E0*1*I_ESP_F3z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+400] = 4.0E0*I_ESP_H3xyz_F3x_aa;
    abcd[iGrid*600+401] = 4.0E0*I_ESP_H2x2yz_F3x_aa-2.0E0*1*I_ESP_F2xz_F3x_a;
    abcd[iGrid*600+402] = 4.0E0*I_ESP_H2xy2z_F3x_aa-2.0E0*1*I_ESP_F2xy_F3x_a;
    abcd[iGrid*600+403] = 4.0E0*I_ESP_Hx3yz_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+404] = 4.0E0*I_ESP_Hx2y2z_F3x_aa-2.0E0*1*I_ESP_Fx2y_F3x_a-2.0E0*1*I_ESP_Fx2z_F3x_a+1*I_ESP_Px_F3x;
    abcd[iGrid*600+405] = 4.0E0*I_ESP_Hxy3z_F3x_aa-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+406] = 4.0E0*I_ESP_H4yz_F3x_aa-2.0E0*3*I_ESP_F2yz_F3x_a;
    abcd[iGrid*600+407] = 4.0E0*I_ESP_H3y2z_F3x_aa-2.0E0*1*I_ESP_F3y_F3x_a-2.0E0*2*I_ESP_Fy2z_F3x_a+2*1*I_ESP_Py_F3x;
    abcd[iGrid*600+408] = 4.0E0*I_ESP_H2y3z_F3x_aa-2.0E0*2*I_ESP_F2yz_F3x_a-2.0E0*1*I_ESP_F3z_F3x_a+2*I_ESP_Pz_F3x;
    abcd[iGrid*600+409] = 4.0E0*I_ESP_Hy4z_F3x_aa-2.0E0*3*I_ESP_Fy2z_F3x_a;
    abcd[iGrid*600+410] = 4.0E0*I_ESP_H3xyz_F2xy_aa;
    abcd[iGrid*600+411] = 4.0E0*I_ESP_H2x2yz_F2xy_aa-2.0E0*1*I_ESP_F2xz_F2xy_a;
    abcd[iGrid*600+412] = 4.0E0*I_ESP_H2xy2z_F2xy_aa-2.0E0*1*I_ESP_F2xy_F2xy_a;
    abcd[iGrid*600+413] = 4.0E0*I_ESP_Hx3yz_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+414] = 4.0E0*I_ESP_Hx2y2z_F2xy_aa-2.0E0*1*I_ESP_Fx2y_F2xy_a-2.0E0*1*I_ESP_Fx2z_F2xy_a+1*I_ESP_Px_F2xy;
    abcd[iGrid*600+415] = 4.0E0*I_ESP_Hxy3z_F2xy_aa-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+416] = 4.0E0*I_ESP_H4yz_F2xy_aa-2.0E0*3*I_ESP_F2yz_F2xy_a;
    abcd[iGrid*600+417] = 4.0E0*I_ESP_H3y2z_F2xy_aa-2.0E0*1*I_ESP_F3y_F2xy_a-2.0E0*2*I_ESP_Fy2z_F2xy_a+2*1*I_ESP_Py_F2xy;
    abcd[iGrid*600+418] = 4.0E0*I_ESP_H2y3z_F2xy_aa-2.0E0*2*I_ESP_F2yz_F2xy_a-2.0E0*1*I_ESP_F3z_F2xy_a+2*I_ESP_Pz_F2xy;
    abcd[iGrid*600+419] = 4.0E0*I_ESP_Hy4z_F2xy_aa-2.0E0*3*I_ESP_Fy2z_F2xy_a;
    abcd[iGrid*600+420] = 4.0E0*I_ESP_H3xyz_F2xz_aa;
    abcd[iGrid*600+421] = 4.0E0*I_ESP_H2x2yz_F2xz_aa-2.0E0*1*I_ESP_F2xz_F2xz_a;
    abcd[iGrid*600+422] = 4.0E0*I_ESP_H2xy2z_F2xz_aa-2.0E0*1*I_ESP_F2xy_F2xz_a;
    abcd[iGrid*600+423] = 4.0E0*I_ESP_Hx3yz_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+424] = 4.0E0*I_ESP_Hx2y2z_F2xz_aa-2.0E0*1*I_ESP_Fx2y_F2xz_a-2.0E0*1*I_ESP_Fx2z_F2xz_a+1*I_ESP_Px_F2xz;
    abcd[iGrid*600+425] = 4.0E0*I_ESP_Hxy3z_F2xz_aa-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+426] = 4.0E0*I_ESP_H4yz_F2xz_aa-2.0E0*3*I_ESP_F2yz_F2xz_a;
    abcd[iGrid*600+427] = 4.0E0*I_ESP_H3y2z_F2xz_aa-2.0E0*1*I_ESP_F3y_F2xz_a-2.0E0*2*I_ESP_Fy2z_F2xz_a+2*1*I_ESP_Py_F2xz;
    abcd[iGrid*600+428] = 4.0E0*I_ESP_H2y3z_F2xz_aa-2.0E0*2*I_ESP_F2yz_F2xz_a-2.0E0*1*I_ESP_F3z_F2xz_a+2*I_ESP_Pz_F2xz;
    abcd[iGrid*600+429] = 4.0E0*I_ESP_Hy4z_F2xz_aa-2.0E0*3*I_ESP_Fy2z_F2xz_a;
    abcd[iGrid*600+430] = 4.0E0*I_ESP_H3xyz_Fx2y_aa;
    abcd[iGrid*600+431] = 4.0E0*I_ESP_H2x2yz_Fx2y_aa-2.0E0*1*I_ESP_F2xz_Fx2y_a;
    abcd[iGrid*600+432] = 4.0E0*I_ESP_H2xy2z_Fx2y_aa-2.0E0*1*I_ESP_F2xy_Fx2y_a;
    abcd[iGrid*600+433] = 4.0E0*I_ESP_Hx3yz_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+434] = 4.0E0*I_ESP_Hx2y2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2y_Fx2y_a-2.0E0*1*I_ESP_Fx2z_Fx2y_a+1*I_ESP_Px_Fx2y;
    abcd[iGrid*600+435] = 4.0E0*I_ESP_Hxy3z_Fx2y_aa-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+436] = 4.0E0*I_ESP_H4yz_Fx2y_aa-2.0E0*3*I_ESP_F2yz_Fx2y_a;
    abcd[iGrid*600+437] = 4.0E0*I_ESP_H3y2z_Fx2y_aa-2.0E0*1*I_ESP_F3y_Fx2y_a-2.0E0*2*I_ESP_Fy2z_Fx2y_a+2*1*I_ESP_Py_Fx2y;
    abcd[iGrid*600+438] = 4.0E0*I_ESP_H2y3z_Fx2y_aa-2.0E0*2*I_ESP_F2yz_Fx2y_a-2.0E0*1*I_ESP_F3z_Fx2y_a+2*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+439] = 4.0E0*I_ESP_Hy4z_Fx2y_aa-2.0E0*3*I_ESP_Fy2z_Fx2y_a;
    abcd[iGrid*600+440] = 4.0E0*I_ESP_H3xyz_Fxyz_aa;
    abcd[iGrid*600+441] = 4.0E0*I_ESP_H2x2yz_Fxyz_aa-2.0E0*1*I_ESP_F2xz_Fxyz_a;
    abcd[iGrid*600+442] = 4.0E0*I_ESP_H2xy2z_Fxyz_aa-2.0E0*1*I_ESP_F2xy_Fxyz_a;
    abcd[iGrid*600+443] = 4.0E0*I_ESP_Hx3yz_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+444] = 4.0E0*I_ESP_Hx2y2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2y_Fxyz_a-2.0E0*1*I_ESP_Fx2z_Fxyz_a+1*I_ESP_Px_Fxyz;
    abcd[iGrid*600+445] = 4.0E0*I_ESP_Hxy3z_Fxyz_aa-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+446] = 4.0E0*I_ESP_H4yz_Fxyz_aa-2.0E0*3*I_ESP_F2yz_Fxyz_a;
    abcd[iGrid*600+447] = 4.0E0*I_ESP_H3y2z_Fxyz_aa-2.0E0*1*I_ESP_F3y_Fxyz_a-2.0E0*2*I_ESP_Fy2z_Fxyz_a+2*1*I_ESP_Py_Fxyz;
    abcd[iGrid*600+448] = 4.0E0*I_ESP_H2y3z_Fxyz_aa-2.0E0*2*I_ESP_F2yz_Fxyz_a-2.0E0*1*I_ESP_F3z_Fxyz_a+2*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+449] = 4.0E0*I_ESP_Hy4z_Fxyz_aa-2.0E0*3*I_ESP_Fy2z_Fxyz_a;
    abcd[iGrid*600+450] = 4.0E0*I_ESP_H3xyz_Fx2z_aa;
    abcd[iGrid*600+451] = 4.0E0*I_ESP_H2x2yz_Fx2z_aa-2.0E0*1*I_ESP_F2xz_Fx2z_a;
    abcd[iGrid*600+452] = 4.0E0*I_ESP_H2xy2z_Fx2z_aa-2.0E0*1*I_ESP_F2xy_Fx2z_a;
    abcd[iGrid*600+453] = 4.0E0*I_ESP_Hx3yz_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+454] = 4.0E0*I_ESP_Hx2y2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2y_Fx2z_a-2.0E0*1*I_ESP_Fx2z_Fx2z_a+1*I_ESP_Px_Fx2z;
    abcd[iGrid*600+455] = 4.0E0*I_ESP_Hxy3z_Fx2z_aa-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+456] = 4.0E0*I_ESP_H4yz_Fx2z_aa-2.0E0*3*I_ESP_F2yz_Fx2z_a;
    abcd[iGrid*600+457] = 4.0E0*I_ESP_H3y2z_Fx2z_aa-2.0E0*1*I_ESP_F3y_Fx2z_a-2.0E0*2*I_ESP_Fy2z_Fx2z_a+2*1*I_ESP_Py_Fx2z;
    abcd[iGrid*600+458] = 4.0E0*I_ESP_H2y3z_Fx2z_aa-2.0E0*2*I_ESP_F2yz_Fx2z_a-2.0E0*1*I_ESP_F3z_Fx2z_a+2*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+459] = 4.0E0*I_ESP_Hy4z_Fx2z_aa-2.0E0*3*I_ESP_Fy2z_Fx2z_a;
    abcd[iGrid*600+460] = 4.0E0*I_ESP_H3xyz_F3y_aa;
    abcd[iGrid*600+461] = 4.0E0*I_ESP_H2x2yz_F3y_aa-2.0E0*1*I_ESP_F2xz_F3y_a;
    abcd[iGrid*600+462] = 4.0E0*I_ESP_H2xy2z_F3y_aa-2.0E0*1*I_ESP_F2xy_F3y_a;
    abcd[iGrid*600+463] = 4.0E0*I_ESP_Hx3yz_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+464] = 4.0E0*I_ESP_Hx2y2z_F3y_aa-2.0E0*1*I_ESP_Fx2y_F3y_a-2.0E0*1*I_ESP_Fx2z_F3y_a+1*I_ESP_Px_F3y;
    abcd[iGrid*600+465] = 4.0E0*I_ESP_Hxy3z_F3y_aa-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+466] = 4.0E0*I_ESP_H4yz_F3y_aa-2.0E0*3*I_ESP_F2yz_F3y_a;
    abcd[iGrid*600+467] = 4.0E0*I_ESP_H3y2z_F3y_aa-2.0E0*1*I_ESP_F3y_F3y_a-2.0E0*2*I_ESP_Fy2z_F3y_a+2*1*I_ESP_Py_F3y;
    abcd[iGrid*600+468] = 4.0E0*I_ESP_H2y3z_F3y_aa-2.0E0*2*I_ESP_F2yz_F3y_a-2.0E0*1*I_ESP_F3z_F3y_a+2*I_ESP_Pz_F3y;
    abcd[iGrid*600+469] = 4.0E0*I_ESP_Hy4z_F3y_aa-2.0E0*3*I_ESP_Fy2z_F3y_a;
    abcd[iGrid*600+470] = 4.0E0*I_ESP_H3xyz_F2yz_aa;
    abcd[iGrid*600+471] = 4.0E0*I_ESP_H2x2yz_F2yz_aa-2.0E0*1*I_ESP_F2xz_F2yz_a;
    abcd[iGrid*600+472] = 4.0E0*I_ESP_H2xy2z_F2yz_aa-2.0E0*1*I_ESP_F2xy_F2yz_a;
    abcd[iGrid*600+473] = 4.0E0*I_ESP_Hx3yz_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+474] = 4.0E0*I_ESP_Hx2y2z_F2yz_aa-2.0E0*1*I_ESP_Fx2y_F2yz_a-2.0E0*1*I_ESP_Fx2z_F2yz_a+1*I_ESP_Px_F2yz;
    abcd[iGrid*600+475] = 4.0E0*I_ESP_Hxy3z_F2yz_aa-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+476] = 4.0E0*I_ESP_H4yz_F2yz_aa-2.0E0*3*I_ESP_F2yz_F2yz_a;
    abcd[iGrid*600+477] = 4.0E0*I_ESP_H3y2z_F2yz_aa-2.0E0*1*I_ESP_F3y_F2yz_a-2.0E0*2*I_ESP_Fy2z_F2yz_a+2*1*I_ESP_Py_F2yz;
    abcd[iGrid*600+478] = 4.0E0*I_ESP_H2y3z_F2yz_aa-2.0E0*2*I_ESP_F2yz_F2yz_a-2.0E0*1*I_ESP_F3z_F2yz_a+2*I_ESP_Pz_F2yz;
    abcd[iGrid*600+479] = 4.0E0*I_ESP_Hy4z_F2yz_aa-2.0E0*3*I_ESP_Fy2z_F2yz_a;
    abcd[iGrid*600+480] = 4.0E0*I_ESP_H3xyz_Fy2z_aa;
    abcd[iGrid*600+481] = 4.0E0*I_ESP_H2x2yz_Fy2z_aa-2.0E0*1*I_ESP_F2xz_Fy2z_a;
    abcd[iGrid*600+482] = 4.0E0*I_ESP_H2xy2z_Fy2z_aa-2.0E0*1*I_ESP_F2xy_Fy2z_a;
    abcd[iGrid*600+483] = 4.0E0*I_ESP_Hx3yz_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+484] = 4.0E0*I_ESP_Hx2y2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2y_Fy2z_a-2.0E0*1*I_ESP_Fx2z_Fy2z_a+1*I_ESP_Px_Fy2z;
    abcd[iGrid*600+485] = 4.0E0*I_ESP_Hxy3z_Fy2z_aa-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+486] = 4.0E0*I_ESP_H4yz_Fy2z_aa-2.0E0*3*I_ESP_F2yz_Fy2z_a;
    abcd[iGrid*600+487] = 4.0E0*I_ESP_H3y2z_Fy2z_aa-2.0E0*1*I_ESP_F3y_Fy2z_a-2.0E0*2*I_ESP_Fy2z_Fy2z_a+2*1*I_ESP_Py_Fy2z;
    abcd[iGrid*600+488] = 4.0E0*I_ESP_H2y3z_Fy2z_aa-2.0E0*2*I_ESP_F2yz_Fy2z_a-2.0E0*1*I_ESP_F3z_Fy2z_a+2*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+489] = 4.0E0*I_ESP_Hy4z_Fy2z_aa-2.0E0*3*I_ESP_Fy2z_Fy2z_a;
    abcd[iGrid*600+490] = 4.0E0*I_ESP_H3xyz_F3z_aa;
    abcd[iGrid*600+491] = 4.0E0*I_ESP_H2x2yz_F3z_aa-2.0E0*1*I_ESP_F2xz_F3z_a;
    abcd[iGrid*600+492] = 4.0E0*I_ESP_H2xy2z_F3z_aa-2.0E0*1*I_ESP_F2xy_F3z_a;
    abcd[iGrid*600+493] = 4.0E0*I_ESP_Hx3yz_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+494] = 4.0E0*I_ESP_Hx2y2z_F3z_aa-2.0E0*1*I_ESP_Fx2y_F3z_a-2.0E0*1*I_ESP_Fx2z_F3z_a+1*I_ESP_Px_F3z;
    abcd[iGrid*600+495] = 4.0E0*I_ESP_Hxy3z_F3z_aa-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+496] = 4.0E0*I_ESP_H4yz_F3z_aa-2.0E0*3*I_ESP_F2yz_F3z_a;
    abcd[iGrid*600+497] = 4.0E0*I_ESP_H3y2z_F3z_aa-2.0E0*1*I_ESP_F3y_F3z_a-2.0E0*2*I_ESP_Fy2z_F3z_a+2*1*I_ESP_Py_F3z;
    abcd[iGrid*600+498] = 4.0E0*I_ESP_H2y3z_F3z_aa-2.0E0*2*I_ESP_F2yz_F3z_a-2.0E0*1*I_ESP_F3z_F3z_a+2*I_ESP_Pz_F3z;
    abcd[iGrid*600+499] = 4.0E0*I_ESP_Hy4z_F3z_aa-2.0E0*3*I_ESP_Fy2z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_aa
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_F_F_a
     * RHS shell quartet name: SQ_ESP_P_F
     ************************************************************/
    abcd[iGrid*600+500] = 4.0E0*I_ESP_H3x2z_F3x_aa-2.0E0*1*I_ESP_F3x_F3x_a;
    abcd[iGrid*600+501] = 4.0E0*I_ESP_H2xy2z_F3x_aa-2.0E0*1*I_ESP_F2xy_F3x_a;
    abcd[iGrid*600+502] = 4.0E0*I_ESP_H2x3z_F3x_aa-2.0E0*1*I_ESP_F2xz_F3x_a-2.0E0*2*I_ESP_F2xz_F3x_a;
    abcd[iGrid*600+503] = 4.0E0*I_ESP_Hx2y2z_F3x_aa-2.0E0*1*I_ESP_Fx2y_F3x_a;
    abcd[iGrid*600+504] = 4.0E0*I_ESP_Hxy3z_F3x_aa-2.0E0*1*I_ESP_Fxyz_F3x_a-2.0E0*2*I_ESP_Fxyz_F3x_a;
    abcd[iGrid*600+505] = 4.0E0*I_ESP_Hx4z_F3x_aa-2.0E0*2*I_ESP_Fx2z_F3x_a-2.0E0*3*I_ESP_Fx2z_F3x_a+2*1*I_ESP_Px_F3x;
    abcd[iGrid*600+506] = 4.0E0*I_ESP_H3y2z_F3x_aa-2.0E0*1*I_ESP_F3y_F3x_a;
    abcd[iGrid*600+507] = 4.0E0*I_ESP_H2y3z_F3x_aa-2.0E0*1*I_ESP_F2yz_F3x_a-2.0E0*2*I_ESP_F2yz_F3x_a;
    abcd[iGrid*600+508] = 4.0E0*I_ESP_Hy4z_F3x_aa-2.0E0*2*I_ESP_Fy2z_F3x_a-2.0E0*3*I_ESP_Fy2z_F3x_a+2*1*I_ESP_Py_F3x;
    abcd[iGrid*600+509] = 4.0E0*I_ESP_H5z_F3x_aa-2.0E0*3*I_ESP_F3z_F3x_a-2.0E0*4*I_ESP_F3z_F3x_a+3*2*I_ESP_Pz_F3x;
    abcd[iGrid*600+510] = 4.0E0*I_ESP_H3x2z_F2xy_aa-2.0E0*1*I_ESP_F3x_F2xy_a;
    abcd[iGrid*600+511] = 4.0E0*I_ESP_H2xy2z_F2xy_aa-2.0E0*1*I_ESP_F2xy_F2xy_a;
    abcd[iGrid*600+512] = 4.0E0*I_ESP_H2x3z_F2xy_aa-2.0E0*1*I_ESP_F2xz_F2xy_a-2.0E0*2*I_ESP_F2xz_F2xy_a;
    abcd[iGrid*600+513] = 4.0E0*I_ESP_Hx2y2z_F2xy_aa-2.0E0*1*I_ESP_Fx2y_F2xy_a;
    abcd[iGrid*600+514] = 4.0E0*I_ESP_Hxy3z_F2xy_aa-2.0E0*1*I_ESP_Fxyz_F2xy_a-2.0E0*2*I_ESP_Fxyz_F2xy_a;
    abcd[iGrid*600+515] = 4.0E0*I_ESP_Hx4z_F2xy_aa-2.0E0*2*I_ESP_Fx2z_F2xy_a-2.0E0*3*I_ESP_Fx2z_F2xy_a+2*1*I_ESP_Px_F2xy;
    abcd[iGrid*600+516] = 4.0E0*I_ESP_H3y2z_F2xy_aa-2.0E0*1*I_ESP_F3y_F2xy_a;
    abcd[iGrid*600+517] = 4.0E0*I_ESP_H2y3z_F2xy_aa-2.0E0*1*I_ESP_F2yz_F2xy_a-2.0E0*2*I_ESP_F2yz_F2xy_a;
    abcd[iGrid*600+518] = 4.0E0*I_ESP_Hy4z_F2xy_aa-2.0E0*2*I_ESP_Fy2z_F2xy_a-2.0E0*3*I_ESP_Fy2z_F2xy_a+2*1*I_ESP_Py_F2xy;
    abcd[iGrid*600+519] = 4.0E0*I_ESP_H5z_F2xy_aa-2.0E0*3*I_ESP_F3z_F2xy_a-2.0E0*4*I_ESP_F3z_F2xy_a+3*2*I_ESP_Pz_F2xy;
    abcd[iGrid*600+520] = 4.0E0*I_ESP_H3x2z_F2xz_aa-2.0E0*1*I_ESP_F3x_F2xz_a;
    abcd[iGrid*600+521] = 4.0E0*I_ESP_H2xy2z_F2xz_aa-2.0E0*1*I_ESP_F2xy_F2xz_a;
    abcd[iGrid*600+522] = 4.0E0*I_ESP_H2x3z_F2xz_aa-2.0E0*1*I_ESP_F2xz_F2xz_a-2.0E0*2*I_ESP_F2xz_F2xz_a;
    abcd[iGrid*600+523] = 4.0E0*I_ESP_Hx2y2z_F2xz_aa-2.0E0*1*I_ESP_Fx2y_F2xz_a;
    abcd[iGrid*600+524] = 4.0E0*I_ESP_Hxy3z_F2xz_aa-2.0E0*1*I_ESP_Fxyz_F2xz_a-2.0E0*2*I_ESP_Fxyz_F2xz_a;
    abcd[iGrid*600+525] = 4.0E0*I_ESP_Hx4z_F2xz_aa-2.0E0*2*I_ESP_Fx2z_F2xz_a-2.0E0*3*I_ESP_Fx2z_F2xz_a+2*1*I_ESP_Px_F2xz;
    abcd[iGrid*600+526] = 4.0E0*I_ESP_H3y2z_F2xz_aa-2.0E0*1*I_ESP_F3y_F2xz_a;
    abcd[iGrid*600+527] = 4.0E0*I_ESP_H2y3z_F2xz_aa-2.0E0*1*I_ESP_F2yz_F2xz_a-2.0E0*2*I_ESP_F2yz_F2xz_a;
    abcd[iGrid*600+528] = 4.0E0*I_ESP_Hy4z_F2xz_aa-2.0E0*2*I_ESP_Fy2z_F2xz_a-2.0E0*3*I_ESP_Fy2z_F2xz_a+2*1*I_ESP_Py_F2xz;
    abcd[iGrid*600+529] = 4.0E0*I_ESP_H5z_F2xz_aa-2.0E0*3*I_ESP_F3z_F2xz_a-2.0E0*4*I_ESP_F3z_F2xz_a+3*2*I_ESP_Pz_F2xz;
    abcd[iGrid*600+530] = 4.0E0*I_ESP_H3x2z_Fx2y_aa-2.0E0*1*I_ESP_F3x_Fx2y_a;
    abcd[iGrid*600+531] = 4.0E0*I_ESP_H2xy2z_Fx2y_aa-2.0E0*1*I_ESP_F2xy_Fx2y_a;
    abcd[iGrid*600+532] = 4.0E0*I_ESP_H2x3z_Fx2y_aa-2.0E0*1*I_ESP_F2xz_Fx2y_a-2.0E0*2*I_ESP_F2xz_Fx2y_a;
    abcd[iGrid*600+533] = 4.0E0*I_ESP_Hx2y2z_Fx2y_aa-2.0E0*1*I_ESP_Fx2y_Fx2y_a;
    abcd[iGrid*600+534] = 4.0E0*I_ESP_Hxy3z_Fx2y_aa-2.0E0*1*I_ESP_Fxyz_Fx2y_a-2.0E0*2*I_ESP_Fxyz_Fx2y_a;
    abcd[iGrid*600+535] = 4.0E0*I_ESP_Hx4z_Fx2y_aa-2.0E0*2*I_ESP_Fx2z_Fx2y_a-2.0E0*3*I_ESP_Fx2z_Fx2y_a+2*1*I_ESP_Px_Fx2y;
    abcd[iGrid*600+536] = 4.0E0*I_ESP_H3y2z_Fx2y_aa-2.0E0*1*I_ESP_F3y_Fx2y_a;
    abcd[iGrid*600+537] = 4.0E0*I_ESP_H2y3z_Fx2y_aa-2.0E0*1*I_ESP_F2yz_Fx2y_a-2.0E0*2*I_ESP_F2yz_Fx2y_a;
    abcd[iGrid*600+538] = 4.0E0*I_ESP_Hy4z_Fx2y_aa-2.0E0*2*I_ESP_Fy2z_Fx2y_a-2.0E0*3*I_ESP_Fy2z_Fx2y_a+2*1*I_ESP_Py_Fx2y;
    abcd[iGrid*600+539] = 4.0E0*I_ESP_H5z_Fx2y_aa-2.0E0*3*I_ESP_F3z_Fx2y_a-2.0E0*4*I_ESP_F3z_Fx2y_a+3*2*I_ESP_Pz_Fx2y;
    abcd[iGrid*600+540] = 4.0E0*I_ESP_H3x2z_Fxyz_aa-2.0E0*1*I_ESP_F3x_Fxyz_a;
    abcd[iGrid*600+541] = 4.0E0*I_ESP_H2xy2z_Fxyz_aa-2.0E0*1*I_ESP_F2xy_Fxyz_a;
    abcd[iGrid*600+542] = 4.0E0*I_ESP_H2x3z_Fxyz_aa-2.0E0*1*I_ESP_F2xz_Fxyz_a-2.0E0*2*I_ESP_F2xz_Fxyz_a;
    abcd[iGrid*600+543] = 4.0E0*I_ESP_Hx2y2z_Fxyz_aa-2.0E0*1*I_ESP_Fx2y_Fxyz_a;
    abcd[iGrid*600+544] = 4.0E0*I_ESP_Hxy3z_Fxyz_aa-2.0E0*1*I_ESP_Fxyz_Fxyz_a-2.0E0*2*I_ESP_Fxyz_Fxyz_a;
    abcd[iGrid*600+545] = 4.0E0*I_ESP_Hx4z_Fxyz_aa-2.0E0*2*I_ESP_Fx2z_Fxyz_a-2.0E0*3*I_ESP_Fx2z_Fxyz_a+2*1*I_ESP_Px_Fxyz;
    abcd[iGrid*600+546] = 4.0E0*I_ESP_H3y2z_Fxyz_aa-2.0E0*1*I_ESP_F3y_Fxyz_a;
    abcd[iGrid*600+547] = 4.0E0*I_ESP_H2y3z_Fxyz_aa-2.0E0*1*I_ESP_F2yz_Fxyz_a-2.0E0*2*I_ESP_F2yz_Fxyz_a;
    abcd[iGrid*600+548] = 4.0E0*I_ESP_Hy4z_Fxyz_aa-2.0E0*2*I_ESP_Fy2z_Fxyz_a-2.0E0*3*I_ESP_Fy2z_Fxyz_a+2*1*I_ESP_Py_Fxyz;
    abcd[iGrid*600+549] = 4.0E0*I_ESP_H5z_Fxyz_aa-2.0E0*3*I_ESP_F3z_Fxyz_a-2.0E0*4*I_ESP_F3z_Fxyz_a+3*2*I_ESP_Pz_Fxyz;
    abcd[iGrid*600+550] = 4.0E0*I_ESP_H3x2z_Fx2z_aa-2.0E0*1*I_ESP_F3x_Fx2z_a;
    abcd[iGrid*600+551] = 4.0E0*I_ESP_H2xy2z_Fx2z_aa-2.0E0*1*I_ESP_F2xy_Fx2z_a;
    abcd[iGrid*600+552] = 4.0E0*I_ESP_H2x3z_Fx2z_aa-2.0E0*1*I_ESP_F2xz_Fx2z_a-2.0E0*2*I_ESP_F2xz_Fx2z_a;
    abcd[iGrid*600+553] = 4.0E0*I_ESP_Hx2y2z_Fx2z_aa-2.0E0*1*I_ESP_Fx2y_Fx2z_a;
    abcd[iGrid*600+554] = 4.0E0*I_ESP_Hxy3z_Fx2z_aa-2.0E0*1*I_ESP_Fxyz_Fx2z_a-2.0E0*2*I_ESP_Fxyz_Fx2z_a;
    abcd[iGrid*600+555] = 4.0E0*I_ESP_Hx4z_Fx2z_aa-2.0E0*2*I_ESP_Fx2z_Fx2z_a-2.0E0*3*I_ESP_Fx2z_Fx2z_a+2*1*I_ESP_Px_Fx2z;
    abcd[iGrid*600+556] = 4.0E0*I_ESP_H3y2z_Fx2z_aa-2.0E0*1*I_ESP_F3y_Fx2z_a;
    abcd[iGrid*600+557] = 4.0E0*I_ESP_H2y3z_Fx2z_aa-2.0E0*1*I_ESP_F2yz_Fx2z_a-2.0E0*2*I_ESP_F2yz_Fx2z_a;
    abcd[iGrid*600+558] = 4.0E0*I_ESP_Hy4z_Fx2z_aa-2.0E0*2*I_ESP_Fy2z_Fx2z_a-2.0E0*3*I_ESP_Fy2z_Fx2z_a+2*1*I_ESP_Py_Fx2z;
    abcd[iGrid*600+559] = 4.0E0*I_ESP_H5z_Fx2z_aa-2.0E0*3*I_ESP_F3z_Fx2z_a-2.0E0*4*I_ESP_F3z_Fx2z_a+3*2*I_ESP_Pz_Fx2z;
    abcd[iGrid*600+560] = 4.0E0*I_ESP_H3x2z_F3y_aa-2.0E0*1*I_ESP_F3x_F3y_a;
    abcd[iGrid*600+561] = 4.0E0*I_ESP_H2xy2z_F3y_aa-2.0E0*1*I_ESP_F2xy_F3y_a;
    abcd[iGrid*600+562] = 4.0E0*I_ESP_H2x3z_F3y_aa-2.0E0*1*I_ESP_F2xz_F3y_a-2.0E0*2*I_ESP_F2xz_F3y_a;
    abcd[iGrid*600+563] = 4.0E0*I_ESP_Hx2y2z_F3y_aa-2.0E0*1*I_ESP_Fx2y_F3y_a;
    abcd[iGrid*600+564] = 4.0E0*I_ESP_Hxy3z_F3y_aa-2.0E0*1*I_ESP_Fxyz_F3y_a-2.0E0*2*I_ESP_Fxyz_F3y_a;
    abcd[iGrid*600+565] = 4.0E0*I_ESP_Hx4z_F3y_aa-2.0E0*2*I_ESP_Fx2z_F3y_a-2.0E0*3*I_ESP_Fx2z_F3y_a+2*1*I_ESP_Px_F3y;
    abcd[iGrid*600+566] = 4.0E0*I_ESP_H3y2z_F3y_aa-2.0E0*1*I_ESP_F3y_F3y_a;
    abcd[iGrid*600+567] = 4.0E0*I_ESP_H2y3z_F3y_aa-2.0E0*1*I_ESP_F2yz_F3y_a-2.0E0*2*I_ESP_F2yz_F3y_a;
    abcd[iGrid*600+568] = 4.0E0*I_ESP_Hy4z_F3y_aa-2.0E0*2*I_ESP_Fy2z_F3y_a-2.0E0*3*I_ESP_Fy2z_F3y_a+2*1*I_ESP_Py_F3y;
    abcd[iGrid*600+569] = 4.0E0*I_ESP_H5z_F3y_aa-2.0E0*3*I_ESP_F3z_F3y_a-2.0E0*4*I_ESP_F3z_F3y_a+3*2*I_ESP_Pz_F3y;
    abcd[iGrid*600+570] = 4.0E0*I_ESP_H3x2z_F2yz_aa-2.0E0*1*I_ESP_F3x_F2yz_a;
    abcd[iGrid*600+571] = 4.0E0*I_ESP_H2xy2z_F2yz_aa-2.0E0*1*I_ESP_F2xy_F2yz_a;
    abcd[iGrid*600+572] = 4.0E0*I_ESP_H2x3z_F2yz_aa-2.0E0*1*I_ESP_F2xz_F2yz_a-2.0E0*2*I_ESP_F2xz_F2yz_a;
    abcd[iGrid*600+573] = 4.0E0*I_ESP_Hx2y2z_F2yz_aa-2.0E0*1*I_ESP_Fx2y_F2yz_a;
    abcd[iGrid*600+574] = 4.0E0*I_ESP_Hxy3z_F2yz_aa-2.0E0*1*I_ESP_Fxyz_F2yz_a-2.0E0*2*I_ESP_Fxyz_F2yz_a;
    abcd[iGrid*600+575] = 4.0E0*I_ESP_Hx4z_F2yz_aa-2.0E0*2*I_ESP_Fx2z_F2yz_a-2.0E0*3*I_ESP_Fx2z_F2yz_a+2*1*I_ESP_Px_F2yz;
    abcd[iGrid*600+576] = 4.0E0*I_ESP_H3y2z_F2yz_aa-2.0E0*1*I_ESP_F3y_F2yz_a;
    abcd[iGrid*600+577] = 4.0E0*I_ESP_H2y3z_F2yz_aa-2.0E0*1*I_ESP_F2yz_F2yz_a-2.0E0*2*I_ESP_F2yz_F2yz_a;
    abcd[iGrid*600+578] = 4.0E0*I_ESP_Hy4z_F2yz_aa-2.0E0*2*I_ESP_Fy2z_F2yz_a-2.0E0*3*I_ESP_Fy2z_F2yz_a+2*1*I_ESP_Py_F2yz;
    abcd[iGrid*600+579] = 4.0E0*I_ESP_H5z_F2yz_aa-2.0E0*3*I_ESP_F3z_F2yz_a-2.0E0*4*I_ESP_F3z_F2yz_a+3*2*I_ESP_Pz_F2yz;
    abcd[iGrid*600+580] = 4.0E0*I_ESP_H3x2z_Fy2z_aa-2.0E0*1*I_ESP_F3x_Fy2z_a;
    abcd[iGrid*600+581] = 4.0E0*I_ESP_H2xy2z_Fy2z_aa-2.0E0*1*I_ESP_F2xy_Fy2z_a;
    abcd[iGrid*600+582] = 4.0E0*I_ESP_H2x3z_Fy2z_aa-2.0E0*1*I_ESP_F2xz_Fy2z_a-2.0E0*2*I_ESP_F2xz_Fy2z_a;
    abcd[iGrid*600+583] = 4.0E0*I_ESP_Hx2y2z_Fy2z_aa-2.0E0*1*I_ESP_Fx2y_Fy2z_a;
    abcd[iGrid*600+584] = 4.0E0*I_ESP_Hxy3z_Fy2z_aa-2.0E0*1*I_ESP_Fxyz_Fy2z_a-2.0E0*2*I_ESP_Fxyz_Fy2z_a;
    abcd[iGrid*600+585] = 4.0E0*I_ESP_Hx4z_Fy2z_aa-2.0E0*2*I_ESP_Fx2z_Fy2z_a-2.0E0*3*I_ESP_Fx2z_Fy2z_a+2*1*I_ESP_Px_Fy2z;
    abcd[iGrid*600+586] = 4.0E0*I_ESP_H3y2z_Fy2z_aa-2.0E0*1*I_ESP_F3y_Fy2z_a;
    abcd[iGrid*600+587] = 4.0E0*I_ESP_H2y3z_Fy2z_aa-2.0E0*1*I_ESP_F2yz_Fy2z_a-2.0E0*2*I_ESP_F2yz_Fy2z_a;
    abcd[iGrid*600+588] = 4.0E0*I_ESP_Hy4z_Fy2z_aa-2.0E0*2*I_ESP_Fy2z_Fy2z_a-2.0E0*3*I_ESP_Fy2z_Fy2z_a+2*1*I_ESP_Py_Fy2z;
    abcd[iGrid*600+589] = 4.0E0*I_ESP_H5z_Fy2z_aa-2.0E0*3*I_ESP_F3z_Fy2z_a-2.0E0*4*I_ESP_F3z_Fy2z_a+3*2*I_ESP_Pz_Fy2z;
    abcd[iGrid*600+590] = 4.0E0*I_ESP_H3x2z_F3z_aa-2.0E0*1*I_ESP_F3x_F3z_a;
    abcd[iGrid*600+591] = 4.0E0*I_ESP_H2xy2z_F3z_aa-2.0E0*1*I_ESP_F2xy_F3z_a;
    abcd[iGrid*600+592] = 4.0E0*I_ESP_H2x3z_F3z_aa-2.0E0*1*I_ESP_F2xz_F3z_a-2.0E0*2*I_ESP_F2xz_F3z_a;
    abcd[iGrid*600+593] = 4.0E0*I_ESP_Hx2y2z_F3z_aa-2.0E0*1*I_ESP_Fx2y_F3z_a;
    abcd[iGrid*600+594] = 4.0E0*I_ESP_Hxy3z_F3z_aa-2.0E0*1*I_ESP_Fxyz_F3z_a-2.0E0*2*I_ESP_Fxyz_F3z_a;
    abcd[iGrid*600+595] = 4.0E0*I_ESP_Hx4z_F3z_aa-2.0E0*2*I_ESP_Fx2z_F3z_a-2.0E0*3*I_ESP_Fx2z_F3z_a+2*1*I_ESP_Px_F3z;
    abcd[iGrid*600+596] = 4.0E0*I_ESP_H3y2z_F3z_aa-2.0E0*1*I_ESP_F3y_F3z_a;
    abcd[iGrid*600+597] = 4.0E0*I_ESP_H2y3z_F3z_aa-2.0E0*1*I_ESP_F2yz_F3z_a-2.0E0*2*I_ESP_F2yz_F3z_a;
    abcd[iGrid*600+598] = 4.0E0*I_ESP_Hy4z_F3z_aa-2.0E0*2*I_ESP_Fy2z_F3z_a-2.0E0*3*I_ESP_Fy2z_F3z_a+2*1*I_ESP_Py_F3z;
    abcd[iGrid*600+599] = 4.0E0*I_ESP_H5z_F3z_aa-2.0E0*3*I_ESP_F3z_F3z_a-2.0E0*4*I_ESP_F3z_F3z_a+3*2*I_ESP_Pz_F3z;
  }
}
