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
// BRA1 as redundant position, total RHS integrals evaluated as: 3563
// BRA2 as redundant position, total RHS integrals evaluated as: 2652
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

void hgp_os_twobodyoverlap_h_sp_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_K7x_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_C5_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_L8x_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xy_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6xyz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x2yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5xy2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x3yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x2y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4xy3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x4yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x3y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x2y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3xy4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x5yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x4y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x3y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x2y4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2xy5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx6yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx5y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx4y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx3y4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx2y5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lxy6z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4y4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3y5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2y6z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ly7z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_C1005_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C1005 = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
    Double fbra  = ifac[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double I_TWOBODYOVERLAP_S_S_vrr = fbra;
    if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;


    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_S_vrr = PAX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_S_vrr = PAY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_S_vrr = PAZ*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_S_vrr = PAY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_S_vrr = PAZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_S_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_S_vrr = PAX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_S_vrr = PAY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_S_vrr = PAY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_S_vrr = PAZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_S_vrr = PAX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_S_vrr = PAY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_S_vrr = PAY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_S_vrr = PAX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_S_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_S_vrr = PAX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_S_vrr = PAY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_S_vrr = PAY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_S_vrr = PAZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_S_vrr = PAX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_S_vrr = PAY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_S_vrr = PAY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_S_vrr = PAX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_S_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_S_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_S_vrr = PAX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_S_vrr = PAY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_S_vrr = PAY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_S_vrr = PAZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_K7x_S_vrr = PAX*I_TWOBODYOVERLAP_I6x_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_S_vrr = PAY*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_S_vrr = PAY*I_TWOBODYOVERLAP_I5xy_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_S_vrr = PAX*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_S_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_S_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_S_vrr = PAX*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_S_vrr = PAY*I_TWOBODYOVERLAP_I6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_S_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_S_vrr = PAY*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_S_vrr = PAZ*I_TWOBODYOVERLAP_I6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_L8x_S_vrr = PAX*I_TWOBODYOVERLAP_K7x_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xy_S_vrr = PAY*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2y_S_vrr = PAY*I_TWOBODYOVERLAP_K6xy_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L6xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3y_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_L5x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L5xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4y_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5y_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_K3xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L3xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6y_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6y_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_K2xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L2xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7y_S_vrr = PAX*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx5y2z_S_vrr = PAX*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx2y5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Lxy6z_S_vrr = PAY*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7z_S_vrr = PAX*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8y_S_vrr = PAY*I_TWOBODYOVERLAP_K7y_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_L6y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L5y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_L4y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3y5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2y5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2y6z_S_vrr = PAY*I_TWOBODYOVERLAP_Ky6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Ly7z_S_vrr = PAY*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8z_S_vrr = PAZ*I_TWOBODYOVERLAP_K7z_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs = ic2*alpha*alpha;
    I_TWOBODYOVERLAP_K7x_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_C5_aa += SQ_TWOBODYOVERLAP_K_S_C5_aa_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C5_a_coefs = ic2*alpha;
    I_TWOBODYOVERLAP_H5x_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C5_a += SQ_TWOBODYOVERLAP_H_S_C5_a_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C5_coefs = ic2;
    I_TWOBODYOVERLAP_F3x_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C5 += SQ_TWOBODYOVERLAP_F_S_C5_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S_C1005_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs = ic2_1*alpha*alpha;
    I_TWOBODYOVERLAP_L8x_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L8x_S_vrr;
    I_TWOBODYOVERLAP_L7xy_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L7xy_S_vrr;
    I_TWOBODYOVERLAP_L7xz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L7xz_S_vrr;
    I_TWOBODYOVERLAP_L6x2y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    I_TWOBODYOVERLAP_L6xyz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L6xyz_S_vrr;
    I_TWOBODYOVERLAP_L6x2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    I_TWOBODYOVERLAP_L5x2yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L5x2yz_S_vrr;
    I_TWOBODYOVERLAP_L5xy2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L5xy2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    I_TWOBODYOVERLAP_L4x3yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4x3yz_S_vrr;
    I_TWOBODYOVERLAP_L4x2y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4x2y2z_S_vrr;
    I_TWOBODYOVERLAP_L4xy3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4xy3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    I_TWOBODYOVERLAP_L3x4yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3x4yz_S_vrr;
    I_TWOBODYOVERLAP_L3x3y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3x3y2z_S_vrr;
    I_TWOBODYOVERLAP_L3x2y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3x2y3z_S_vrr;
    I_TWOBODYOVERLAP_L3xy4z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3xy4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    I_TWOBODYOVERLAP_L2x5yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x5yz_S_vrr;
    I_TWOBODYOVERLAP_L2x4y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x4y2z_S_vrr;
    I_TWOBODYOVERLAP_L2x3y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x3y3z_S_vrr;
    I_TWOBODYOVERLAP_L2x2y4z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x2y4z_S_vrr;
    I_TWOBODYOVERLAP_L2xy5z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2xy5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    I_TWOBODYOVERLAP_Lx6yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx6yz_S_vrr;
    I_TWOBODYOVERLAP_Lx5y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx5y2z_S_vrr;
    I_TWOBODYOVERLAP_Lx4y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx4y3z_S_vrr;
    I_TWOBODYOVERLAP_Lx3y4z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx3y4z_S_vrr;
    I_TWOBODYOVERLAP_Lx2y5z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx2y5z_S_vrr;
    I_TWOBODYOVERLAP_Lxy6z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lxy6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    I_TWOBODYOVERLAP_L8y_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L8y_S_vrr;
    I_TWOBODYOVERLAP_L7yz_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L7yz_S_vrr;
    I_TWOBODYOVERLAP_L6y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    I_TWOBODYOVERLAP_L5y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    I_TWOBODYOVERLAP_L4y4z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    I_TWOBODYOVERLAP_L3y5z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    I_TWOBODYOVERLAP_L2y6z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    I_TWOBODYOVERLAP_Ly7z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Ly7z_S_vrr;
    I_TWOBODYOVERLAP_L8z_S_C1005_aa += SQ_TWOBODYOVERLAP_L_S_C1005_aa_coefs*I_TWOBODYOVERLAP_L8z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_C1005_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs = ic2_1*alpha*alpha;
    I_TWOBODYOVERLAP_K7x_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_C1005_aa += SQ_TWOBODYOVERLAP_K_S_C1005_aa_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_I6x_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_H5x_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C1005_a += SQ_TWOBODYOVERLAP_H_S_C1005_a_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_G4x_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_F3x_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C1005 += SQ_TWOBODYOVERLAP_F_S_C1005_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1005
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1005
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px_C1005 = I_TWOBODYOVERLAP_G4x_S_C1005+ABX*I_TWOBODYOVERLAP_F3x_S_C1005;
  Double I_TWOBODYOVERLAP_F2xy_Px_C1005 = I_TWOBODYOVERLAP_G3xy_S_C1005+ABX*I_TWOBODYOVERLAP_F2xy_S_C1005;
  Double I_TWOBODYOVERLAP_F2xz_Px_C1005 = I_TWOBODYOVERLAP_G3xz_S_C1005+ABX*I_TWOBODYOVERLAP_F2xz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2y_Px_C1005 = I_TWOBODYOVERLAP_G2x2y_S_C1005+ABX*I_TWOBODYOVERLAP_Fx2y_S_C1005;
  Double I_TWOBODYOVERLAP_Fxyz_Px_C1005 = I_TWOBODYOVERLAP_G2xyz_S_C1005+ABX*I_TWOBODYOVERLAP_Fxyz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2z_Px_C1005 = I_TWOBODYOVERLAP_G2x2z_S_C1005+ABX*I_TWOBODYOVERLAP_Fx2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3y_Px_C1005 = I_TWOBODYOVERLAP_Gx3y_S_C1005+ABX*I_TWOBODYOVERLAP_F3y_S_C1005;
  Double I_TWOBODYOVERLAP_F2yz_Px_C1005 = I_TWOBODYOVERLAP_Gx2yz_S_C1005+ABX*I_TWOBODYOVERLAP_F2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Fy2z_Px_C1005 = I_TWOBODYOVERLAP_Gxy2z_S_C1005+ABX*I_TWOBODYOVERLAP_Fy2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3z_Px_C1005 = I_TWOBODYOVERLAP_Gx3z_S_C1005+ABX*I_TWOBODYOVERLAP_F3z_S_C1005;
  Double I_TWOBODYOVERLAP_F3x_Py_C1005 = I_TWOBODYOVERLAP_G3xy_S_C1005+ABY*I_TWOBODYOVERLAP_F3x_S_C1005;
  Double I_TWOBODYOVERLAP_F2xy_Py_C1005 = I_TWOBODYOVERLAP_G2x2y_S_C1005+ABY*I_TWOBODYOVERLAP_F2xy_S_C1005;
  Double I_TWOBODYOVERLAP_F2xz_Py_C1005 = I_TWOBODYOVERLAP_G2xyz_S_C1005+ABY*I_TWOBODYOVERLAP_F2xz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2y_Py_C1005 = I_TWOBODYOVERLAP_Gx3y_S_C1005+ABY*I_TWOBODYOVERLAP_Fx2y_S_C1005;
  Double I_TWOBODYOVERLAP_Fxyz_Py_C1005 = I_TWOBODYOVERLAP_Gx2yz_S_C1005+ABY*I_TWOBODYOVERLAP_Fxyz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2z_Py_C1005 = I_TWOBODYOVERLAP_Gxy2z_S_C1005+ABY*I_TWOBODYOVERLAP_Fx2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3y_Py_C1005 = I_TWOBODYOVERLAP_G4y_S_C1005+ABY*I_TWOBODYOVERLAP_F3y_S_C1005;
  Double I_TWOBODYOVERLAP_F2yz_Py_C1005 = I_TWOBODYOVERLAP_G3yz_S_C1005+ABY*I_TWOBODYOVERLAP_F2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Fy2z_Py_C1005 = I_TWOBODYOVERLAP_G2y2z_S_C1005+ABY*I_TWOBODYOVERLAP_Fy2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3z_Py_C1005 = I_TWOBODYOVERLAP_Gy3z_S_C1005+ABY*I_TWOBODYOVERLAP_F3z_S_C1005;
  Double I_TWOBODYOVERLAP_F3x_Pz_C1005 = I_TWOBODYOVERLAP_G3xz_S_C1005+ABZ*I_TWOBODYOVERLAP_F3x_S_C1005;
  Double I_TWOBODYOVERLAP_F2xy_Pz_C1005 = I_TWOBODYOVERLAP_G2xyz_S_C1005+ABZ*I_TWOBODYOVERLAP_F2xy_S_C1005;
  Double I_TWOBODYOVERLAP_F2xz_Pz_C1005 = I_TWOBODYOVERLAP_G2x2z_S_C1005+ABZ*I_TWOBODYOVERLAP_F2xz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2y_Pz_C1005 = I_TWOBODYOVERLAP_Gx2yz_S_C1005+ABZ*I_TWOBODYOVERLAP_Fx2y_S_C1005;
  Double I_TWOBODYOVERLAP_Fxyz_Pz_C1005 = I_TWOBODYOVERLAP_Gxy2z_S_C1005+ABZ*I_TWOBODYOVERLAP_Fxyz_S_C1005;
  Double I_TWOBODYOVERLAP_Fx2z_Pz_C1005 = I_TWOBODYOVERLAP_Gx3z_S_C1005+ABZ*I_TWOBODYOVERLAP_Fx2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3y_Pz_C1005 = I_TWOBODYOVERLAP_G3yz_S_C1005+ABZ*I_TWOBODYOVERLAP_F3y_S_C1005;
  Double I_TWOBODYOVERLAP_F2yz_Pz_C1005 = I_TWOBODYOVERLAP_G2y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_F2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Fy2z_Pz_C1005 = I_TWOBODYOVERLAP_Gy3z_S_C1005+ABZ*I_TWOBODYOVERLAP_Fy2z_S_C1005;
  Double I_TWOBODYOVERLAP_F3z_Pz_C1005 = I_TWOBODYOVERLAP_G4z_S_C1005+ABZ*I_TWOBODYOVERLAP_F3z_S_C1005;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px_C1005_a = I_TWOBODYOVERLAP_I6x_S_C1005_a+ABX*I_TWOBODYOVERLAP_H5x_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xy_Px_C1005_a = I_TWOBODYOVERLAP_I5xy_S_C1005_a+ABX*I_TWOBODYOVERLAP_H4xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xz_Px_C1005_a = I_TWOBODYOVERLAP_I5xz_S_C1005_a+ABX*I_TWOBODYOVERLAP_H4xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2y_Px_C1005_a = I_TWOBODYOVERLAP_I4x2y_S_C1005_a+ABX*I_TWOBODYOVERLAP_H3x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3xyz_Px_C1005_a = I_TWOBODYOVERLAP_I4xyz_S_C1005_a+ABX*I_TWOBODYOVERLAP_H3xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2z_Px_C1005_a = I_TWOBODYOVERLAP_I4x2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H3x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3y_Px_C1005_a = I_TWOBODYOVERLAP_I3x3y_S_C1005_a+ABX*I_TWOBODYOVERLAP_H2x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a = I_TWOBODYOVERLAP_I3x2yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_H2x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a = I_TWOBODYOVERLAP_I3xy2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H2xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3z_Px_C1005_a = I_TWOBODYOVERLAP_I3x3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H2x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4y_Px_C1005_a = I_TWOBODYOVERLAP_I2x4y_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hx4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a = I_TWOBODYOVERLAP_I2x3yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hx3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a = I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a = I_TWOBODYOVERLAP_I2xy3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hxy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4z_Px_C1005_a = I_TWOBODYOVERLAP_I2x4z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hx4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5y_Px_C1005_a = I_TWOBODYOVERLAP_Ix5y_S_C1005_a+ABX*I_TWOBODYOVERLAP_H5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4yz_Px_C1005_a = I_TWOBODYOVERLAP_Ix4yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_H4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3y2z_Px_C1005_a = I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2y3z_Px_C1005_a = I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hy4z_Px_C1005_a = I_TWOBODYOVERLAP_Ixy4z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Hy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5z_Px_C1005_a = I_TWOBODYOVERLAP_Ix5z_S_C1005_a+ABX*I_TWOBODYOVERLAP_H5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5x_Py_C1005_a = I_TWOBODYOVERLAP_I5xy_S_C1005_a+ABY*I_TWOBODYOVERLAP_H5x_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xy_Py_C1005_a = I_TWOBODYOVERLAP_I4x2y_S_C1005_a+ABY*I_TWOBODYOVERLAP_H4xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xz_Py_C1005_a = I_TWOBODYOVERLAP_I4xyz_S_C1005_a+ABY*I_TWOBODYOVERLAP_H4xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2y_Py_C1005_a = I_TWOBODYOVERLAP_I3x3y_S_C1005_a+ABY*I_TWOBODYOVERLAP_H3x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3xyz_Py_C1005_a = I_TWOBODYOVERLAP_I3x2yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_H3xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2z_Py_C1005_a = I_TWOBODYOVERLAP_I3xy2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H3x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3y_Py_C1005_a = I_TWOBODYOVERLAP_I2x4y_S_C1005_a+ABY*I_TWOBODYOVERLAP_H2x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a = I_TWOBODYOVERLAP_I2x3yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_H2x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a = I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H2xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3z_Py_C1005_a = I_TWOBODYOVERLAP_I2xy3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H2x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4y_Py_C1005_a = I_TWOBODYOVERLAP_Ix5y_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hx4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a = I_TWOBODYOVERLAP_Ix4yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hx3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a = I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a = I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hxy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4z_Py_C1005_a = I_TWOBODYOVERLAP_Ixy4z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hx4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5y_Py_C1005_a = I_TWOBODYOVERLAP_I6y_S_C1005_a+ABY*I_TWOBODYOVERLAP_H5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4yz_Py_C1005_a = I_TWOBODYOVERLAP_I5yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_H4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3y2z_Py_C1005_a = I_TWOBODYOVERLAP_I4y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2y3z_Py_C1005_a = I_TWOBODYOVERLAP_I3y3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hy4z_Py_C1005_a = I_TWOBODYOVERLAP_I2y4z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Hy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5z_Py_C1005_a = I_TWOBODYOVERLAP_Iy5z_S_C1005_a+ABY*I_TWOBODYOVERLAP_H5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5x_Pz_C1005_a = I_TWOBODYOVERLAP_I5xz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H5x_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xy_Pz_C1005_a = I_TWOBODYOVERLAP_I4xyz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H4xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4xz_Pz_C1005_a = I_TWOBODYOVERLAP_I4x2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H4xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a = I_TWOBODYOVERLAP_I3x2yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H3x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a = I_TWOBODYOVERLAP_I3xy2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H3xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a = I_TWOBODYOVERLAP_I3x3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H3x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a = I_TWOBODYOVERLAP_I2x3yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H2x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a = I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a = I_TWOBODYOVERLAP_I2xy3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a = I_TWOBODYOVERLAP_I2x4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H2x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a = I_TWOBODYOVERLAP_Ix4yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hx4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a = I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a = I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a = I_TWOBODYOVERLAP_Ixy4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a = I_TWOBODYOVERLAP_Ix5z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hx4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5y_Pz_C1005_a = I_TWOBODYOVERLAP_I5yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_H4yz_Pz_C1005_a = I_TWOBODYOVERLAP_I4y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a = I_TWOBODYOVERLAP_I3y3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a = I_TWOBODYOVERLAP_I2y4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a = I_TWOBODYOVERLAP_Iy5z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Hy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_H5z_Pz_C1005_a = I_TWOBODYOVERLAP_I6z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_H5z_S_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C1005_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_Px_C1005_aa = I_TWOBODYOVERLAP_L8x_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K7x_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xy_Px_C1005_aa = I_TWOBODYOVERLAP_L7xy_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K6xy_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xz_Px_C1005_aa = I_TWOBODYOVERLAP_L7xz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K6xz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Px_C1005_aa = I_TWOBODYOVERLAP_L6x2y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K5x2y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Px_C1005_aa = I_TWOBODYOVERLAP_L6xyz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K5xyz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Px_C1005_aa = I_TWOBODYOVERLAP_L6x2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K5x2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Px_C1005_aa = I_TWOBODYOVERLAP_L5x3y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K4x3y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa = I_TWOBODYOVERLAP_L5x2yz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K4x2yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa = I_TWOBODYOVERLAP_L5xy2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K4xy2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Px_C1005_aa = I_TWOBODYOVERLAP_L5x3z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K4x3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Px_C1005_aa = I_TWOBODYOVERLAP_L4x4y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3x4y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa = I_TWOBODYOVERLAP_L4x3yz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3x3yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa = I_TWOBODYOVERLAP_L4x2y2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3x2y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa = I_TWOBODYOVERLAP_L4xy3z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3xy3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Px_C1005_aa = I_TWOBODYOVERLAP_L4x4z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3x4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Px_C1005_aa = I_TWOBODYOVERLAP_L3x5y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2x5y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa = I_TWOBODYOVERLAP_L3x4yz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2x4yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa = I_TWOBODYOVERLAP_L3x3y2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2x3y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa = I_TWOBODYOVERLAP_L3x2y3z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2x2y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa = I_TWOBODYOVERLAP_L3xy4z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2xy4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Px_C1005_aa = I_TWOBODYOVERLAP_L3x5z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2x5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Px_C1005_aa = I_TWOBODYOVERLAP_L2x6y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx6y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Px_C1005_aa = I_TWOBODYOVERLAP_L2x5yz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx5yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa = I_TWOBODYOVERLAP_L2x4y2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx4y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa = I_TWOBODYOVERLAP_L2x3y3z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx3y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa = I_TWOBODYOVERLAP_L2x2y4z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx2y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Px_C1005_aa = I_TWOBODYOVERLAP_L2xy5z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kxy5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Px_C1005_aa = I_TWOBODYOVERLAP_L2x6z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Kx6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7y_Px_C1005_aa = I_TWOBODYOVERLAP_Lx7y_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K7y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6yz_Px_C1005_aa = I_TWOBODYOVERLAP_Lx6yz_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K6yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Px_C1005_aa = I_TWOBODYOVERLAP_Lx5y2z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K5y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Px_C1005_aa = I_TWOBODYOVERLAP_Lx4y3z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K4y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Px_C1005_aa = I_TWOBODYOVERLAP_Lx3y4z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K3y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Px_C1005_aa = I_TWOBODYOVERLAP_Lx2y5z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K2y5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Px_C1005_aa = I_TWOBODYOVERLAP_Lxy6z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_Ky6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7z_Px_C1005_aa = I_TWOBODYOVERLAP_Lx7z_S_C1005_aa+ABX*I_TWOBODYOVERLAP_K7z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7x_Py_C1005_aa = I_TWOBODYOVERLAP_L7xy_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K7x_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xy_Py_C1005_aa = I_TWOBODYOVERLAP_L6x2y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K6xy_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xz_Py_C1005_aa = I_TWOBODYOVERLAP_L6xyz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K6xz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Py_C1005_aa = I_TWOBODYOVERLAP_L5x3y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K5x2y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Py_C1005_aa = I_TWOBODYOVERLAP_L5x2yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K5xyz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Py_C1005_aa = I_TWOBODYOVERLAP_L5xy2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K5x2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Py_C1005_aa = I_TWOBODYOVERLAP_L4x4y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K4x3y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa = I_TWOBODYOVERLAP_L4x3yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K4x2yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa = I_TWOBODYOVERLAP_L4x2y2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K4xy2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Py_C1005_aa = I_TWOBODYOVERLAP_L4xy3z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K4x3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Py_C1005_aa = I_TWOBODYOVERLAP_L3x5y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3x4y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa = I_TWOBODYOVERLAP_L3x4yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3x3yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa = I_TWOBODYOVERLAP_L3x3y2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3x2y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa = I_TWOBODYOVERLAP_L3x2y3z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3xy3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Py_C1005_aa = I_TWOBODYOVERLAP_L3xy4z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3x4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Py_C1005_aa = I_TWOBODYOVERLAP_L2x6y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2x5y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa = I_TWOBODYOVERLAP_L2x5yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2x4yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa = I_TWOBODYOVERLAP_L2x4y2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2x3y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa = I_TWOBODYOVERLAP_L2x3y3z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2x2y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa = I_TWOBODYOVERLAP_L2x2y4z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2xy4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Py_C1005_aa = I_TWOBODYOVERLAP_L2xy5z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2x5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Py_C1005_aa = I_TWOBODYOVERLAP_Lx7y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx6y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Py_C1005_aa = I_TWOBODYOVERLAP_Lx6yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx5yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa = I_TWOBODYOVERLAP_Lx5y2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx4y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa = I_TWOBODYOVERLAP_Lx4y3z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx3y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa = I_TWOBODYOVERLAP_Lx3y4z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx2y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Py_C1005_aa = I_TWOBODYOVERLAP_Lx2y5z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kxy5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Py_C1005_aa = I_TWOBODYOVERLAP_Lxy6z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Kx6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7y_Py_C1005_aa = I_TWOBODYOVERLAP_L8y_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K7y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6yz_Py_C1005_aa = I_TWOBODYOVERLAP_L7yz_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K6yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Py_C1005_aa = I_TWOBODYOVERLAP_L6y2z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K5y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Py_C1005_aa = I_TWOBODYOVERLAP_L5y3z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K4y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Py_C1005_aa = I_TWOBODYOVERLAP_L4y4z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K3y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Py_C1005_aa = I_TWOBODYOVERLAP_L3y5z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K2y5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Py_C1005_aa = I_TWOBODYOVERLAP_L2y6z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_Ky6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7z_Py_C1005_aa = I_TWOBODYOVERLAP_Ly7z_S_C1005_aa+ABY*I_TWOBODYOVERLAP_K7z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7x_Pz_C1005_aa = I_TWOBODYOVERLAP_L7xz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K7x_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xy_Pz_C1005_aa = I_TWOBODYOVERLAP_L6xyz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K6xy_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6xz_Pz_C1005_aa = I_TWOBODYOVERLAP_L6x2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K6xz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Pz_C1005_aa = I_TWOBODYOVERLAP_L5x2yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K5x2y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Pz_C1005_aa = I_TWOBODYOVERLAP_L5xy2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K5xyz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Pz_C1005_aa = I_TWOBODYOVERLAP_L5x3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K5x2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Pz_C1005_aa = I_TWOBODYOVERLAP_L4x3yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K4x3y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa = I_TWOBODYOVERLAP_L4x2y2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K4x2yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa = I_TWOBODYOVERLAP_L4xy3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K4xy2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Pz_C1005_aa = I_TWOBODYOVERLAP_L4x4z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K4x3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Pz_C1005_aa = I_TWOBODYOVERLAP_L3x4yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3x4y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa = I_TWOBODYOVERLAP_L3x3y2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3x3yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa = I_TWOBODYOVERLAP_L3x2y3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3x2y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa = I_TWOBODYOVERLAP_L3xy4z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3xy3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Pz_C1005_aa = I_TWOBODYOVERLAP_L3x5z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3x4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Pz_C1005_aa = I_TWOBODYOVERLAP_L2x5yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2x5y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa = I_TWOBODYOVERLAP_L2x4y2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2x4yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa = I_TWOBODYOVERLAP_L2x3y3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2x3y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa = I_TWOBODYOVERLAP_L2x2y4z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2x2y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa = I_TWOBODYOVERLAP_L2xy5z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2xy4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Pz_C1005_aa = I_TWOBODYOVERLAP_L2x6z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2x5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx6yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx6y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx5y2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx5yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx4y3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx4y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx3y4z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx3y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx2y5z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx2y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Pz_C1005_aa = I_TWOBODYOVERLAP_Lxy6z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kxy5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Pz_C1005_aa = I_TWOBODYOVERLAP_Lx7z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Kx6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7y_Pz_C1005_aa = I_TWOBODYOVERLAP_L7yz_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K7y_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K6yz_Pz_C1005_aa = I_TWOBODYOVERLAP_L6y2z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K6yz_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Pz_C1005_aa = I_TWOBODYOVERLAP_L5y3z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K5y2z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Pz_C1005_aa = I_TWOBODYOVERLAP_L4y4z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K4y3z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Pz_C1005_aa = I_TWOBODYOVERLAP_L3y5z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K3y4z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Pz_C1005_aa = I_TWOBODYOVERLAP_L2y6z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K2y5z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Pz_C1005_aa = I_TWOBODYOVERLAP_Ly7z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_Ky6z_S_C1005_aa;
  Double I_TWOBODYOVERLAP_K7z_Pz_C1005_aa = I_TWOBODYOVERLAP_L8z_S_C1005_aa+ABZ*I_TWOBODYOVERLAP_K7z_S_C1005_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[0] = 4.0E0*I_TWOBODYOVERLAP_K7x_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H5x_S_C5_a-2.0E0*6*I_TWOBODYOVERLAP_H5x_S_C5_a+5*4*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[1] = 4.0E0*I_TWOBODYOVERLAP_K6xy_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xy_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_H4xy_S_C5_a+4*3*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[2] = 4.0E0*I_TWOBODYOVERLAP_K6xz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xz_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_H4xz_S_C5_a+4*3*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[3] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_S_C5_a+3*2*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[4] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3xyz_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_S_C5_a+3*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[5] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_S_C5_a+3*2*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[6] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3y_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_S_C5_a+2*1*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[7] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a+2*1*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[8] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a+2*1*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[9] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_S_C5_a+2*1*I_TWOBODYOVERLAP_F3z_S_C5;
  abcd[10] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_S_C5_a;
  abcd[11] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[12] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a;
  abcd[13] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[14] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_S_C5_a;
  abcd[15] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_S_C5_a;
  abcd[16] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_S_C5_a;
  abcd[17] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_S_C5_a;
  abcd[18] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_S_C5_a;
  abcd[19] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_S_C5_a;
  abcd[20] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[21] = 4.0E0*I_TWOBODYOVERLAP_K7x_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5x_Px_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5x_Px_C1005_a+5*4*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[22] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xy_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Px_C1005_a+4*3*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[23] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xz_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Px_C1005_a+4*3*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[24] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[25] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[26] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[27] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[28] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[29] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[30] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[31] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a;
  abcd[32] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[33] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a;
  abcd[34] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[35] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a;
  abcd[36] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Px_C1005_a;
  abcd[37] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Px_C1005_a;
  abcd[38] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a;
  abcd[39] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a;
  abcd[40] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a;
  abcd[41] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Px_C1005_a;
  abcd[42] = 4.0E0*I_TWOBODYOVERLAP_K7x_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5x_Py_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5x_Py_C1005_a+5*4*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[43] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xy_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Py_C1005_a+4*3*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[44] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xz_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Py_C1005_a+4*3*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[45] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[46] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[47] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[48] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[49] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[50] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[51] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[52] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a;
  abcd[53] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[54] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a;
  abcd[55] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[56] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a;
  abcd[57] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Py_C1005_a;
  abcd[58] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Py_C1005_a;
  abcd[59] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a;
  abcd[60] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a;
  abcd[61] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a;
  abcd[62] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Py_C1005_a;
  abcd[63] = 4.0E0*I_TWOBODYOVERLAP_K7x_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5x_Pz_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5x_Pz_C1005_a+5*4*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[64] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[65] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[66] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[67] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[68] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[69] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[70] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[71] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[72] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Pz_C1005;
  abcd[73] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a;
  abcd[74] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[75] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a;
  abcd[76] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[77] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a;
  abcd[78] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Pz_C1005_a;
  abcd[79] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a;
  abcd[80] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a;
  abcd[81] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a;
  abcd[82] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a;
  abcd[83] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[84] = 4.0E0*I_TWOBODYOVERLAP_K6xy_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xy_S_C5_a;
  abcd[85] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_S_C5_a+4*1*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[86] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[87] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xy_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_S_C5_a+3*2*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[88] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a+3*1*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[89] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a;
  abcd[90] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_S_C5_a+2*3*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[91] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a+2*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[92] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a+2*1*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[93] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[94] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H5y_S_C5_a+4*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[95] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H4yz_S_C5_a+3*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[96] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_S_C5_a+2*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[97] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_S_C5_a+1*I_TWOBODYOVERLAP_F3z_S_C5;
  abcd[98] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_S_C5_a;
  abcd[99] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_S_C5_a;
  abcd[100] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[101] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a;
  abcd[102] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[103] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_S_C5_a;
  abcd[104] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_S_C5_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[105] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Px_C1005_a;
  abcd[106] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[107] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[108] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[109] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a+3*1*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[110] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a;
  abcd[111] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[112] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[113] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[114] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[115] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5y_Px_C1005_a+4*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[116] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Px_C1005_a+3*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[117] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a+2*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[118] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a+1*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[119] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a;
  abcd[120] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a;
  abcd[121] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[122] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a;
  abcd[123] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[124] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a;
  abcd[125] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Px_C1005_aa;
  abcd[126] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Py_C1005_a;
  abcd[127] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[128] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[129] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[130] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a+3*1*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[131] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a;
  abcd[132] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[133] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[134] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[135] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[136] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5y_Py_C1005_a+4*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[137] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Py_C1005_a+3*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[138] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a+2*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[139] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a+1*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[140] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a;
  abcd[141] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a;
  abcd[142] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[143] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a;
  abcd[144] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[145] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a;
  abcd[146] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Py_C1005_aa;
  abcd[147] = 4.0E0*I_TWOBODYOVERLAP_K6xy_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a;
  abcd[148] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[149] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[150] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[151] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a+3*1*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[152] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a;
  abcd[153] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[154] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[155] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[156] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[157] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5y_Pz_C1005_a+4*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[158] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a+3*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[159] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a+2*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[160] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a+1*I_TWOBODYOVERLAP_F3z_Pz_C1005;
  abcd[161] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a;
  abcd[162] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a;
  abcd[163] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[164] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a;
  abcd[165] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[166] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a;
  abcd[167] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Pz_C1005_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[168] = 4.0E0*I_TWOBODYOVERLAP_K6xz_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xz_S_C5_a;
  abcd[169] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[170] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_S_C5_a+4*1*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[171] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a;
  abcd[172] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a+3*1*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[173] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xz_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_S_C5_a+3*2*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[174] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[175] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a+2*1*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[176] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a+2*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[177] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_S_C5_a+2*3*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[178] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_S_C5_a;
  abcd[179] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_S_C5_a+1*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[180] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_S_C5_a+2*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[181] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_S_C5_a+3*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[182] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_S_C5_a+4*I_TWOBODYOVERLAP_F3z_S_C5;
  abcd[183] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_S_C5_aa;
  abcd[184] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_S_C5_a;
  abcd[185] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[186] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a;
  abcd[187] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[188] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[189] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Px_C1005_a;
  abcd[190] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[191] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[192] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a;
  abcd[193] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a+3*1*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[194] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[195] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[196] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[197] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[198] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[199] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Px_C1005_a;
  abcd[200] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a+1*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[201] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a+2*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[202] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a+3*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[203] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Px_C1005_a+4*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[204] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Px_C1005_aa;
  abcd[205] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a;
  abcd[206] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[207] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a;
  abcd[208] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[209] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a;
  abcd[210] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Py_C1005_a;
  abcd[211] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[212] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[213] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a;
  abcd[214] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a+3*1*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[215] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[216] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[217] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[218] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[219] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[220] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Py_C1005_a;
  abcd[221] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a+1*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[222] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a+2*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[223] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a+3*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[224] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Py_C1005_a+4*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[225] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Py_C1005_aa;
  abcd[226] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a;
  abcd[227] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[228] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a;
  abcd[229] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[230] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a;
  abcd[231] = 4.0E0*I_TWOBODYOVERLAP_K6xz_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a;
  abcd[232] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[233] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a+4*1*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[234] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a;
  abcd[235] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a+3*1*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[236] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[237] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[238] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[239] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[240] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a+2*3*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[241] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a;
  abcd[242] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a+1*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[243] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a+2*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[244] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a+3*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[245] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Pz_C1005_a+4*I_TWOBODYOVERLAP_F3z_Pz_C1005;
  abcd[246] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Pz_C1005_aa;
  abcd[247] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a;
  abcd[248] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[249] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a;
  abcd[250] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[251] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[252] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_S_C5_a;
  abcd[253] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H4xy_S_C5_a;
  abcd[254] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_S_C5_a;
  abcd[255] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2y_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_S_C5_a+2*1*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[256] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[257] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_S_C5_a;
  abcd[258] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_S_C5_a+3*2*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[259] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a+2*1*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[260] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2xy2z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_S_C5_a;
  abcd[261] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_S_C5_a;
  abcd[262] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4y_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_S_C5_a+4*3*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[263] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx3yz_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_S_C5_a+3*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[264] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a+2*1*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[265] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[266] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_S_C5_a;
  abcd[267] = 4.0E0*I_TWOBODYOVERLAP_K7y_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H5y_S_C5_a-2.0E0*6*I_TWOBODYOVERLAP_H5y_S_C5_a+5*4*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[268] = 4.0E0*I_TWOBODYOVERLAP_K6yz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H4yz_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_H4yz_S_C5_a+4*3*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[269] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_S_C5_a+3*2*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[270] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2y3z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_S_C5_a+2*1*I_TWOBODYOVERLAP_F3z_S_C5;
  abcd[271] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_S_C5_a;
  abcd[272] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[273] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Px_C1005_a;
  abcd[274] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Px_C1005_a;
  abcd[275] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Px_C1005_a;
  abcd[276] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[277] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[278] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a;
  abcd[279] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[280] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[281] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a;
  abcd[282] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a;
  abcd[283] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[284] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[285] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[286] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[287] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a;
  abcd[288] = 4.0E0*I_TWOBODYOVERLAP_K7y_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5y_Px_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5y_Px_C1005_a+5*4*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[289] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4yz_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Px_C1005_a+4*3*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[290] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[291] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[292] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a;
  abcd[293] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Px_C1005_a;
  abcd[294] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Py_C1005_a;
  abcd[295] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Py_C1005_a;
  abcd[296] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Py_C1005_a;
  abcd[297] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[298] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[299] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a;
  abcd[300] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[301] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[302] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a;
  abcd[303] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a;
  abcd[304] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[305] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[306] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[307] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[308] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a;
  abcd[309] = 4.0E0*I_TWOBODYOVERLAP_K7y_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5y_Py_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5y_Py_C1005_a+5*4*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[310] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4yz_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Py_C1005_a+4*3*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[311] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[312] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[313] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a;
  abcd[314] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Py_C1005_a;
  abcd[315] = 4.0E0*I_TWOBODYOVERLAP_K5x2y_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Pz_C1005_a;
  abcd[316] = 4.0E0*I_TWOBODYOVERLAP_K4x3y_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a;
  abcd[317] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a;
  abcd[318] = 4.0E0*I_TWOBODYOVERLAP_K3x4y_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[319] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[320] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a;
  abcd[321] = 4.0E0*I_TWOBODYOVERLAP_K2x5y_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[322] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[323] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a;
  abcd[324] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a;
  abcd[325] = 4.0E0*I_TWOBODYOVERLAP_Kx6y_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[326] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[327] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[328] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[329] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a;
  abcd[330] = 4.0E0*I_TWOBODYOVERLAP_K7y_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5y_Pz_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5y_Pz_C1005_a+5*4*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[331] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[332] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[333] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3z_Pz_C1005;
  abcd[334] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a;
  abcd[335] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[336] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_S_C5_aa;
  abcd[337] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_S_C5_a;
  abcd[338] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_S_C5_a;
  abcd[339] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[340] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_S_C5_a+1*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[341] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[342] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_S_C5_a;
  abcd[343] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_S_C5_a+2*1*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[344] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_S_C5_a+2*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[345] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a;
  abcd[346] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[347] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a+3*1*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[348] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_S_C5_a+2*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[349] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_S_C5_a+3*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[350] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_S_C5_a;
  abcd[351] = 4.0E0*I_TWOBODYOVERLAP_K6yz_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H4yz_S_C5_a;
  abcd[352] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_S_C5_a+4*1*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[353] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H4yz_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_S_C5_a+3*2*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[354] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_S_C5_a+2*3*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[355] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_S_C5_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_S_C5_a+4*I_TWOBODYOVERLAP_F3z_S_C5;
  abcd[356] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[357] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Px_C1005_aa;
  abcd[358] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Px_C1005_a;
  abcd[359] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Px_C1005_a;
  abcd[360] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[361] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a+1*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[362] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[363] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a;
  abcd[364] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[365] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a+2*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[366] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a;
  abcd[367] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[368] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a+3*1*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[369] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[370] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a+3*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[371] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a;
  abcd[372] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Px_C1005_a;
  abcd[373] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a+4*1*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[374] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[375] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a+2*3*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[376] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Px_C1005_a+4*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[377] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a;
  abcd[378] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Py_C1005_aa;
  abcd[379] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Py_C1005_a;
  abcd[380] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Py_C1005_a;
  abcd[381] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[382] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a+1*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[383] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[384] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a;
  abcd[385] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[386] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a+2*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[387] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a;
  abcd[388] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[389] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a+3*1*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[390] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[391] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a+3*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[392] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a;
  abcd[393] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Py_C1005_a;
  abcd[394] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a+4*1*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[395] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[396] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a+2*3*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[397] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Py_C1005_a+4*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[398] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a;
  abcd[399] = 4.0E0*I_TWOBODYOVERLAP_K5xyz_Pz_C1005_aa;
  abcd[400] = 4.0E0*I_TWOBODYOVERLAP_K4x2yz_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a;
  abcd[401] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a;
  abcd[402] = 4.0E0*I_TWOBODYOVERLAP_K3x3yz_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[403] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a+1*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[404] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[405] = 4.0E0*I_TWOBODYOVERLAP_K2x4yz_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a;
  abcd[406] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[407] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a+2*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[408] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a;
  abcd[409] = 4.0E0*I_TWOBODYOVERLAP_Kx5yz_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[410] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a+3*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[411] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a+2*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[412] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a+3*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[413] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a;
  abcd[414] = 4.0E0*I_TWOBODYOVERLAP_K6yz_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a;
  abcd[415] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a+4*1*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[416] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[417] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a+2*3*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[418] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a-2.0E0*1*I_TWOBODYOVERLAP_H5z_Pz_C1005_a+4*I_TWOBODYOVERLAP_F3z_Pz_C1005;
  abcd[419] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C5_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C5
   ************************************************************/
  abcd[420] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_S_C5_a;
  abcd[421] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_S_C5_a;
  abcd[422] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H4xz_S_C5_a;
  abcd[423] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_S_C5_a;
  abcd[424] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_S_C5_a;
  abcd[425] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_S_C5_a+2*1*I_TWOBODYOVERLAP_F3x_S_C5;
  abcd[426] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_S_C5_a;
  abcd[427] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x2yz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_S_C5_a;
  abcd[428] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_S_C5_a+2*1*I_TWOBODYOVERLAP_F2xy_S_C5;
  abcd[429] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_S_C5_a+3*2*I_TWOBODYOVERLAP_F2xz_S_C5;
  abcd[430] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_S_C5_a;
  abcd[431] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_S_C5_a;
  abcd[432] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_S_C5_a+2*1*I_TWOBODYOVERLAP_Fx2y_S_C5;
  abcd[433] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_Hxy3z_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_S_C5_a+3*2*I_TWOBODYOVERLAP_Fxyz_S_C5;
  abcd[434] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4z_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_S_C5_a+4*3*I_TWOBODYOVERLAP_Fx2z_S_C5;
  abcd[435] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_S_C5_a;
  abcd[436] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_S_C5_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_S_C5_a-2.0E0*2*I_TWOBODYOVERLAP_H4yz_S_C5_a;
  abcd[437] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_S_C5_aa-2.0E0*2*I_TWOBODYOVERLAP_H3y2z_S_C5_a-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_S_C5_a+2*1*I_TWOBODYOVERLAP_F3y_S_C5;
  abcd[438] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_S_C5_aa-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_S_C5_a-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_S_C5_a+3*2*I_TWOBODYOVERLAP_F2yz_S_C5;
  abcd[439] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_S_C5_aa-2.0E0*4*I_TWOBODYOVERLAP_Hy4z_S_C5_a-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_S_C5_a+4*3*I_TWOBODYOVERLAP_Fy2z_S_C5;
  abcd[440] = 4.0E0*I_TWOBODYOVERLAP_K7z_S_C5_aa-2.0E0*5*I_TWOBODYOVERLAP_H5z_S_C5_a-2.0E0*6*I_TWOBODYOVERLAP_H5z_S_C5_a+5*4*I_TWOBODYOVERLAP_F3z_S_C5;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_C1005_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1005
   ************************************************************/
  abcd[441] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Px_C1005_a;
  abcd[442] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Px_C1005_a;
  abcd[443] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Px_C1005_a;
  abcd[444] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Px_C1005_a;
  abcd[445] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Px_C1005_a;
  abcd[446] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Px_C1005;
  abcd[447] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Px_C1005_a;
  abcd[448] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Px_C1005_a;
  abcd[449] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Px_C1005;
  abcd[450] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Px_C1005;
  abcd[451] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Px_C1005_a;
  abcd[452] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Px_C1005_a;
  abcd[453] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Px_C1005;
  abcd[454] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Px_C1005;
  abcd[455] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Px_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2z_Px_C1005;
  abcd[456] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Px_C1005_a;
  abcd[457] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Px_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Px_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Px_C1005_a;
  abcd[458] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Px_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Px_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Px_C1005;
  abcd[459] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Px_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Px_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Px_C1005;
  abcd[460] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Px_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Px_C1005_a+4*3*I_TWOBODYOVERLAP_Fy2z_Px_C1005;
  abcd[461] = 4.0E0*I_TWOBODYOVERLAP_K7z_Px_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5z_Px_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5z_Px_C1005_a+5*4*I_TWOBODYOVERLAP_F3z_Px_C1005;
  abcd[462] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Py_C1005_a;
  abcd[463] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Py_C1005_a;
  abcd[464] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Py_C1005_a;
  abcd[465] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Py_C1005_a;
  abcd[466] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Py_C1005_a;
  abcd[467] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Py_C1005;
  abcd[468] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Py_C1005_a;
  abcd[469] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Py_C1005_a;
  abcd[470] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Py_C1005;
  abcd[471] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Py_C1005;
  abcd[472] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Py_C1005_a;
  abcd[473] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Py_C1005_a;
  abcd[474] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Py_C1005;
  abcd[475] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Py_C1005;
  abcd[476] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Py_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2z_Py_C1005;
  abcd[477] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Py_C1005_a;
  abcd[478] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Py_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Py_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Py_C1005_a;
  abcd[479] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Py_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Py_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Py_C1005;
  abcd[480] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Py_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Py_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Py_C1005;
  abcd[481] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Py_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Py_C1005_a+4*3*I_TWOBODYOVERLAP_Fy2z_Py_C1005;
  abcd[482] = 4.0E0*I_TWOBODYOVERLAP_K7z_Py_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5z_Py_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5z_Py_C1005_a+5*4*I_TWOBODYOVERLAP_F3z_Py_C1005;
  abcd[483] = 4.0E0*I_TWOBODYOVERLAP_K5x2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5x_Pz_C1005_a;
  abcd[484] = 4.0E0*I_TWOBODYOVERLAP_K4xy2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xy_Pz_C1005_a;
  abcd[485] = 4.0E0*I_TWOBODYOVERLAP_K4x3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4xz_Pz_C1005_a;
  abcd[486] = 4.0E0*I_TWOBODYOVERLAP_K3x2y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3x2y_Pz_C1005_a;
  abcd[487] = 4.0E0*I_TWOBODYOVERLAP_K3xy3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H3xyz_Pz_C1005_a;
  abcd[488] = 4.0E0*I_TWOBODYOVERLAP_K3x4z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3x2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3x_Pz_C1005;
  abcd[489] = 4.0E0*I_TWOBODYOVERLAP_K2x3y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x3y_Pz_C1005_a;
  abcd[490] = 4.0E0*I_TWOBODYOVERLAP_K2x2y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H2x2yz_Pz_C1005_a;
  abcd[491] = 4.0E0*I_TWOBODYOVERLAP_K2xy4z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H2xy2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F2xy_Pz_C1005;
  abcd[492] = 4.0E0*I_TWOBODYOVERLAP_K2x5z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2x3z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2xz_Pz_C1005;
  abcd[493] = 4.0E0*I_TWOBODYOVERLAP_Kx4y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx4y_Pz_C1005_a;
  abcd[494] = 4.0E0*I_TWOBODYOVERLAP_Kx3y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_Hx3yz_Pz_C1005_a;
  abcd[495] = 4.0E0*I_TWOBODYOVERLAP_Kx2y4z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1005;
  abcd[496] = 4.0E0*I_TWOBODYOVERLAP_Kxy5z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_Hxy3z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1005;
  abcd[497] = 4.0E0*I_TWOBODYOVERLAP_Kx6z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hx4z_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_Fx2z_Pz_C1005;
  abcd[498] = 4.0E0*I_TWOBODYOVERLAP_K5y2z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H5y_Pz_C1005_a;
  abcd[499] = 4.0E0*I_TWOBODYOVERLAP_K4y3z_Pz_C1005_aa-2.0E0*1*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a-2.0E0*2*I_TWOBODYOVERLAP_H4yz_Pz_C1005_a;
  abcd[500] = 4.0E0*I_TWOBODYOVERLAP_K3y4z_Pz_C1005_aa-2.0E0*2*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a-2.0E0*3*I_TWOBODYOVERLAP_H3y2z_Pz_C1005_a+2*1*I_TWOBODYOVERLAP_F3y_Pz_C1005;
  abcd[501] = 4.0E0*I_TWOBODYOVERLAP_K2y5z_Pz_C1005_aa-2.0E0*3*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a-2.0E0*4*I_TWOBODYOVERLAP_H2y3z_Pz_C1005_a+3*2*I_TWOBODYOVERLAP_F2yz_Pz_C1005;
  abcd[502] = 4.0E0*I_TWOBODYOVERLAP_Ky6z_Pz_C1005_aa-2.0E0*4*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a-2.0E0*5*I_TWOBODYOVERLAP_Hy4z_Pz_C1005_a+4*3*I_TWOBODYOVERLAP_Fy2z_Pz_C1005;
  abcd[503] = 4.0E0*I_TWOBODYOVERLAP_K7z_Pz_C1005_aa-2.0E0*5*I_TWOBODYOVERLAP_H5z_Pz_C1005_a-2.0E0*6*I_TWOBODYOVERLAP_H5z_Pz_C1005_a+5*4*I_TWOBODYOVERLAP_F3z_Pz_C1005;
}
