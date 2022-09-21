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
// BRA1  BRA2
// X  X
// X  Y
// X  Z
// Y  X
// Y  Y
// Y  Z
// Z  X
// Z  Y
// Z  Z
// ####
// BRA2  BRA2
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_nai_g_f_d2(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
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
  Double I_NAI_G4x_S = 0.0E0;
  Double I_NAI_G3xy_S = 0.0E0;
  Double I_NAI_G3xz_S = 0.0E0;
  Double I_NAI_G2x2y_S = 0.0E0;
  Double I_NAI_G2xyz_S = 0.0E0;
  Double I_NAI_G2x2z_S = 0.0E0;
  Double I_NAI_Gx3y_S = 0.0E0;
  Double I_NAI_Gx2yz_S = 0.0E0;
  Double I_NAI_Gxy2z_S = 0.0E0;
  Double I_NAI_Gx3z_S = 0.0E0;
  Double I_NAI_G4y_S = 0.0E0;
  Double I_NAI_G3yz_S = 0.0E0;
  Double I_NAI_G2y2z_S = 0.0E0;
  Double I_NAI_Gy3z_S = 0.0E0;
  Double I_NAI_G4z_S = 0.0E0;
  Double I_NAI_K7x_S_a = 0.0E0;
  Double I_NAI_K6xy_S_a = 0.0E0;
  Double I_NAI_K6xz_S_a = 0.0E0;
  Double I_NAI_K5x2y_S_a = 0.0E0;
  Double I_NAI_K5xyz_S_a = 0.0E0;
  Double I_NAI_K5x2z_S_a = 0.0E0;
  Double I_NAI_K4x3y_S_a = 0.0E0;
  Double I_NAI_K4x2yz_S_a = 0.0E0;
  Double I_NAI_K4xy2z_S_a = 0.0E0;
  Double I_NAI_K4x3z_S_a = 0.0E0;
  Double I_NAI_K3x4y_S_a = 0.0E0;
  Double I_NAI_K3x3yz_S_a = 0.0E0;
  Double I_NAI_K3x2y2z_S_a = 0.0E0;
  Double I_NAI_K3xy3z_S_a = 0.0E0;
  Double I_NAI_K3x4z_S_a = 0.0E0;
  Double I_NAI_K2x5y_S_a = 0.0E0;
  Double I_NAI_K2x4yz_S_a = 0.0E0;
  Double I_NAI_K2x3y2z_S_a = 0.0E0;
  Double I_NAI_K2x2y3z_S_a = 0.0E0;
  Double I_NAI_K2xy4z_S_a = 0.0E0;
  Double I_NAI_K2x5z_S_a = 0.0E0;
  Double I_NAI_Kx6y_S_a = 0.0E0;
  Double I_NAI_Kx5yz_S_a = 0.0E0;
  Double I_NAI_Kx4y2z_S_a = 0.0E0;
  Double I_NAI_Kx3y3z_S_a = 0.0E0;
  Double I_NAI_Kx2y4z_S_a = 0.0E0;
  Double I_NAI_Kxy5z_S_a = 0.0E0;
  Double I_NAI_Kx6z_S_a = 0.0E0;
  Double I_NAI_K7y_S_a = 0.0E0;
  Double I_NAI_K6yz_S_a = 0.0E0;
  Double I_NAI_K5y2z_S_a = 0.0E0;
  Double I_NAI_K4y3z_S_a = 0.0E0;
  Double I_NAI_K3y4z_S_a = 0.0E0;
  Double I_NAI_K2y5z_S_a = 0.0E0;
  Double I_NAI_Ky6z_S_a = 0.0E0;
  Double I_NAI_K7z_S_a = 0.0E0;
  Double I_NAI_I6x_S_a = 0.0E0;
  Double I_NAI_I5xy_S_a = 0.0E0;
  Double I_NAI_I5xz_S_a = 0.0E0;
  Double I_NAI_I4x2y_S_a = 0.0E0;
  Double I_NAI_I4xyz_S_a = 0.0E0;
  Double I_NAI_I4x2z_S_a = 0.0E0;
  Double I_NAI_I3x3y_S_a = 0.0E0;
  Double I_NAI_I3x2yz_S_a = 0.0E0;
  Double I_NAI_I3xy2z_S_a = 0.0E0;
  Double I_NAI_I3x3z_S_a = 0.0E0;
  Double I_NAI_I2x4y_S_a = 0.0E0;
  Double I_NAI_I2x3yz_S_a = 0.0E0;
  Double I_NAI_I2x2y2z_S_a = 0.0E0;
  Double I_NAI_I2xy3z_S_a = 0.0E0;
  Double I_NAI_I2x4z_S_a = 0.0E0;
  Double I_NAI_Ix5y_S_a = 0.0E0;
  Double I_NAI_Ix4yz_S_a = 0.0E0;
  Double I_NAI_Ix3y2z_S_a = 0.0E0;
  Double I_NAI_Ix2y3z_S_a = 0.0E0;
  Double I_NAI_Ixy4z_S_a = 0.0E0;
  Double I_NAI_Ix5z_S_a = 0.0E0;
  Double I_NAI_I6y_S_a = 0.0E0;
  Double I_NAI_I5yz_S_a = 0.0E0;
  Double I_NAI_I4y2z_S_a = 0.0E0;
  Double I_NAI_I3y3z_S_a = 0.0E0;
  Double I_NAI_I2y4z_S_a = 0.0E0;
  Double I_NAI_Iy5z_S_a = 0.0E0;
  Double I_NAI_I6z_S_a = 0.0E0;
  Double I_NAI_H5x_S_a = 0.0E0;
  Double I_NAI_H4xy_S_a = 0.0E0;
  Double I_NAI_H4xz_S_a = 0.0E0;
  Double I_NAI_H3x2y_S_a = 0.0E0;
  Double I_NAI_H3xyz_S_a = 0.0E0;
  Double I_NAI_H3x2z_S_a = 0.0E0;
  Double I_NAI_H2x3y_S_a = 0.0E0;
  Double I_NAI_H2x2yz_S_a = 0.0E0;
  Double I_NAI_H2xy2z_S_a = 0.0E0;
  Double I_NAI_H2x3z_S_a = 0.0E0;
  Double I_NAI_Hx4y_S_a = 0.0E0;
  Double I_NAI_Hx3yz_S_a = 0.0E0;
  Double I_NAI_Hx2y2z_S_a = 0.0E0;
  Double I_NAI_Hxy3z_S_a = 0.0E0;
  Double I_NAI_Hx4z_S_a = 0.0E0;
  Double I_NAI_H5y_S_a = 0.0E0;
  Double I_NAI_H4yz_S_a = 0.0E0;
  Double I_NAI_H3y2z_S_a = 0.0E0;
  Double I_NAI_H2y3z_S_a = 0.0E0;
  Double I_NAI_Hy4z_S_a = 0.0E0;
  Double I_NAI_H5z_S_a = 0.0E0;
  Double I_NAI_F3x_S = 0.0E0;
  Double I_NAI_F2xy_S = 0.0E0;
  Double I_NAI_F2xz_S = 0.0E0;
  Double I_NAI_Fx2y_S = 0.0E0;
  Double I_NAI_Fxyz_S = 0.0E0;
  Double I_NAI_Fx2z_S = 0.0E0;
  Double I_NAI_F3y_S = 0.0E0;
  Double I_NAI_F2yz_S = 0.0E0;
  Double I_NAI_Fy2z_S = 0.0E0;
  Double I_NAI_F3z_S = 0.0E0;
  Double I_NAI_M9x_S_aa = 0.0E0;
  Double I_NAI_M8xy_S_aa = 0.0E0;
  Double I_NAI_M8xz_S_aa = 0.0E0;
  Double I_NAI_M7x2y_S_aa = 0.0E0;
  Double I_NAI_M7xyz_S_aa = 0.0E0;
  Double I_NAI_M7x2z_S_aa = 0.0E0;
  Double I_NAI_M6x3y_S_aa = 0.0E0;
  Double I_NAI_M6x2yz_S_aa = 0.0E0;
  Double I_NAI_M6xy2z_S_aa = 0.0E0;
  Double I_NAI_M6x3z_S_aa = 0.0E0;
  Double I_NAI_M5x4y_S_aa = 0.0E0;
  Double I_NAI_M5x3yz_S_aa = 0.0E0;
  Double I_NAI_M5x2y2z_S_aa = 0.0E0;
  Double I_NAI_M5xy3z_S_aa = 0.0E0;
  Double I_NAI_M5x4z_S_aa = 0.0E0;
  Double I_NAI_M4x5y_S_aa = 0.0E0;
  Double I_NAI_M4x4yz_S_aa = 0.0E0;
  Double I_NAI_M4x3y2z_S_aa = 0.0E0;
  Double I_NAI_M4x2y3z_S_aa = 0.0E0;
  Double I_NAI_M4xy4z_S_aa = 0.0E0;
  Double I_NAI_M4x5z_S_aa = 0.0E0;
  Double I_NAI_M3x6y_S_aa = 0.0E0;
  Double I_NAI_M3x5yz_S_aa = 0.0E0;
  Double I_NAI_M3x4y2z_S_aa = 0.0E0;
  Double I_NAI_M3x3y3z_S_aa = 0.0E0;
  Double I_NAI_M3x2y4z_S_aa = 0.0E0;
  Double I_NAI_M3xy5z_S_aa = 0.0E0;
  Double I_NAI_M3x6z_S_aa = 0.0E0;
  Double I_NAI_M2x7y_S_aa = 0.0E0;
  Double I_NAI_M2x6yz_S_aa = 0.0E0;
  Double I_NAI_M2x5y2z_S_aa = 0.0E0;
  Double I_NAI_M2x4y3z_S_aa = 0.0E0;
  Double I_NAI_M2x3y4z_S_aa = 0.0E0;
  Double I_NAI_M2x2y5z_S_aa = 0.0E0;
  Double I_NAI_M2xy6z_S_aa = 0.0E0;
  Double I_NAI_M2x7z_S_aa = 0.0E0;
  Double I_NAI_Mx8y_S_aa = 0.0E0;
  Double I_NAI_Mx7yz_S_aa = 0.0E0;
  Double I_NAI_Mx6y2z_S_aa = 0.0E0;
  Double I_NAI_Mx5y3z_S_aa = 0.0E0;
  Double I_NAI_Mx4y4z_S_aa = 0.0E0;
  Double I_NAI_Mx3y5z_S_aa = 0.0E0;
  Double I_NAI_Mx2y6z_S_aa = 0.0E0;
  Double I_NAI_Mxy7z_S_aa = 0.0E0;
  Double I_NAI_Mx8z_S_aa = 0.0E0;
  Double I_NAI_M9y_S_aa = 0.0E0;
  Double I_NAI_M8yz_S_aa = 0.0E0;
  Double I_NAI_M7y2z_S_aa = 0.0E0;
  Double I_NAI_M6y3z_S_aa = 0.0E0;
  Double I_NAI_M5y4z_S_aa = 0.0E0;
  Double I_NAI_M4y5z_S_aa = 0.0E0;
  Double I_NAI_M3y6z_S_aa = 0.0E0;
  Double I_NAI_M2y7z_S_aa = 0.0E0;
  Double I_NAI_My8z_S_aa = 0.0E0;
  Double I_NAI_M9z_S_aa = 0.0E0;
  Double I_NAI_L8x_S_aa = 0.0E0;
  Double I_NAI_L7xy_S_aa = 0.0E0;
  Double I_NAI_L7xz_S_aa = 0.0E0;
  Double I_NAI_L6x2y_S_aa = 0.0E0;
  Double I_NAI_L6xyz_S_aa = 0.0E0;
  Double I_NAI_L6x2z_S_aa = 0.0E0;
  Double I_NAI_L5x3y_S_aa = 0.0E0;
  Double I_NAI_L5x2yz_S_aa = 0.0E0;
  Double I_NAI_L5xy2z_S_aa = 0.0E0;
  Double I_NAI_L5x3z_S_aa = 0.0E0;
  Double I_NAI_L4x4y_S_aa = 0.0E0;
  Double I_NAI_L4x3yz_S_aa = 0.0E0;
  Double I_NAI_L4x2y2z_S_aa = 0.0E0;
  Double I_NAI_L4xy3z_S_aa = 0.0E0;
  Double I_NAI_L4x4z_S_aa = 0.0E0;
  Double I_NAI_L3x5y_S_aa = 0.0E0;
  Double I_NAI_L3x4yz_S_aa = 0.0E0;
  Double I_NAI_L3x3y2z_S_aa = 0.0E0;
  Double I_NAI_L3x2y3z_S_aa = 0.0E0;
  Double I_NAI_L3xy4z_S_aa = 0.0E0;
  Double I_NAI_L3x5z_S_aa = 0.0E0;
  Double I_NAI_L2x6y_S_aa = 0.0E0;
  Double I_NAI_L2x5yz_S_aa = 0.0E0;
  Double I_NAI_L2x4y2z_S_aa = 0.0E0;
  Double I_NAI_L2x3y3z_S_aa = 0.0E0;
  Double I_NAI_L2x2y4z_S_aa = 0.0E0;
  Double I_NAI_L2xy5z_S_aa = 0.0E0;
  Double I_NAI_L2x6z_S_aa = 0.0E0;
  Double I_NAI_Lx7y_S_aa = 0.0E0;
  Double I_NAI_Lx6yz_S_aa = 0.0E0;
  Double I_NAI_Lx5y2z_S_aa = 0.0E0;
  Double I_NAI_Lx4y3z_S_aa = 0.0E0;
  Double I_NAI_Lx3y4z_S_aa = 0.0E0;
  Double I_NAI_Lx2y5z_S_aa = 0.0E0;
  Double I_NAI_Lxy6z_S_aa = 0.0E0;
  Double I_NAI_Lx7z_S_aa = 0.0E0;
  Double I_NAI_L8y_S_aa = 0.0E0;
  Double I_NAI_L7yz_S_aa = 0.0E0;
  Double I_NAI_L6y2z_S_aa = 0.0E0;
  Double I_NAI_L5y3z_S_aa = 0.0E0;
  Double I_NAI_L4y4z_S_aa = 0.0E0;
  Double I_NAI_L3y5z_S_aa = 0.0E0;
  Double I_NAI_L2y6z_S_aa = 0.0E0;
  Double I_NAI_Ly7z_S_aa = 0.0E0;
  Double I_NAI_L8z_S_aa = 0.0E0;
  Double I_NAI_K7x_S_aa = 0.0E0;
  Double I_NAI_K6xy_S_aa = 0.0E0;
  Double I_NAI_K6xz_S_aa = 0.0E0;
  Double I_NAI_K5x2y_S_aa = 0.0E0;
  Double I_NAI_K5xyz_S_aa = 0.0E0;
  Double I_NAI_K5x2z_S_aa = 0.0E0;
  Double I_NAI_K4x3y_S_aa = 0.0E0;
  Double I_NAI_K4x2yz_S_aa = 0.0E0;
  Double I_NAI_K4xy2z_S_aa = 0.0E0;
  Double I_NAI_K4x3z_S_aa = 0.0E0;
  Double I_NAI_K3x4y_S_aa = 0.0E0;
  Double I_NAI_K3x3yz_S_aa = 0.0E0;
  Double I_NAI_K3x2y2z_S_aa = 0.0E0;
  Double I_NAI_K3xy3z_S_aa = 0.0E0;
  Double I_NAI_K3x4z_S_aa = 0.0E0;
  Double I_NAI_K2x5y_S_aa = 0.0E0;
  Double I_NAI_K2x4yz_S_aa = 0.0E0;
  Double I_NAI_K2x3y2z_S_aa = 0.0E0;
  Double I_NAI_K2x2y3z_S_aa = 0.0E0;
  Double I_NAI_K2xy4z_S_aa = 0.0E0;
  Double I_NAI_K2x5z_S_aa = 0.0E0;
  Double I_NAI_Kx6y_S_aa = 0.0E0;
  Double I_NAI_Kx5yz_S_aa = 0.0E0;
  Double I_NAI_Kx4y2z_S_aa = 0.0E0;
  Double I_NAI_Kx3y3z_S_aa = 0.0E0;
  Double I_NAI_Kx2y4z_S_aa = 0.0E0;
  Double I_NAI_Kxy5z_S_aa = 0.0E0;
  Double I_NAI_Kx6z_S_aa = 0.0E0;
  Double I_NAI_K7y_S_aa = 0.0E0;
  Double I_NAI_K6yz_S_aa = 0.0E0;
  Double I_NAI_K5y2z_S_aa = 0.0E0;
  Double I_NAI_K4y3z_S_aa = 0.0E0;
  Double I_NAI_K3y4z_S_aa = 0.0E0;
  Double I_NAI_K2y5z_S_aa = 0.0E0;
  Double I_NAI_Ky6z_S_aa = 0.0E0;
  Double I_NAI_K7z_S_aa = 0.0E0;
  Double I_NAI_I6x_S_aa = 0.0E0;
  Double I_NAI_I5xy_S_aa = 0.0E0;
  Double I_NAI_I5xz_S_aa = 0.0E0;
  Double I_NAI_I4x2y_S_aa = 0.0E0;
  Double I_NAI_I4xyz_S_aa = 0.0E0;
  Double I_NAI_I4x2z_S_aa = 0.0E0;
  Double I_NAI_I3x3y_S_aa = 0.0E0;
  Double I_NAI_I3x2yz_S_aa = 0.0E0;
  Double I_NAI_I3xy2z_S_aa = 0.0E0;
  Double I_NAI_I3x3z_S_aa = 0.0E0;
  Double I_NAI_I2x4y_S_aa = 0.0E0;
  Double I_NAI_I2x3yz_S_aa = 0.0E0;
  Double I_NAI_I2x2y2z_S_aa = 0.0E0;
  Double I_NAI_I2xy3z_S_aa = 0.0E0;
  Double I_NAI_I2x4z_S_aa = 0.0E0;
  Double I_NAI_Ix5y_S_aa = 0.0E0;
  Double I_NAI_Ix4yz_S_aa = 0.0E0;
  Double I_NAI_Ix3y2z_S_aa = 0.0E0;
  Double I_NAI_Ix2y3z_S_aa = 0.0E0;
  Double I_NAI_Ixy4z_S_aa = 0.0E0;
  Double I_NAI_Ix5z_S_aa = 0.0E0;
  Double I_NAI_I6y_S_aa = 0.0E0;
  Double I_NAI_I5yz_S_aa = 0.0E0;
  Double I_NAI_I4y2z_S_aa = 0.0E0;
  Double I_NAI_I3y3z_S_aa = 0.0E0;
  Double I_NAI_I2y4z_S_aa = 0.0E0;
  Double I_NAI_Iy5z_S_aa = 0.0E0;
  Double I_NAI_I6z_S_aa = 0.0E0;
  Double I_NAI_G4x_S_a = 0.0E0;
  Double I_NAI_G3xy_S_a = 0.0E0;
  Double I_NAI_G3xz_S_a = 0.0E0;
  Double I_NAI_G2x2y_S_a = 0.0E0;
  Double I_NAI_G2xyz_S_a = 0.0E0;
  Double I_NAI_G2x2z_S_a = 0.0E0;
  Double I_NAI_Gx3y_S_a = 0.0E0;
  Double I_NAI_Gx2yz_S_a = 0.0E0;
  Double I_NAI_Gxy2z_S_a = 0.0E0;
  Double I_NAI_Gx3z_S_a = 0.0E0;
  Double I_NAI_G4y_S_a = 0.0E0;
  Double I_NAI_G3yz_S_a = 0.0E0;
  Double I_NAI_G2y2z_S_a = 0.0E0;
  Double I_NAI_Gy3z_S_a = 0.0E0;
  Double I_NAI_G4z_S_a = 0.0E0;
  Double I_NAI_D2x_S = 0.0E0;
  Double I_NAI_Dxy_S = 0.0E0;
  Double I_NAI_Dxz_S = 0.0E0;
  Double I_NAI_D2y_S = 0.0E0;
  Double I_NAI_Dyz_S = 0.0E0;
  Double I_NAI_D2z_S = 0.0E0;
  Double I_NAI_K7x_S_b = 0.0E0;
  Double I_NAI_K6xy_S_b = 0.0E0;
  Double I_NAI_K6xz_S_b = 0.0E0;
  Double I_NAI_K5x2y_S_b = 0.0E0;
  Double I_NAI_K5xyz_S_b = 0.0E0;
  Double I_NAI_K5x2z_S_b = 0.0E0;
  Double I_NAI_K4x3y_S_b = 0.0E0;
  Double I_NAI_K4x2yz_S_b = 0.0E0;
  Double I_NAI_K4xy2z_S_b = 0.0E0;
  Double I_NAI_K4x3z_S_b = 0.0E0;
  Double I_NAI_K3x4y_S_b = 0.0E0;
  Double I_NAI_K3x3yz_S_b = 0.0E0;
  Double I_NAI_K3x2y2z_S_b = 0.0E0;
  Double I_NAI_K3xy3z_S_b = 0.0E0;
  Double I_NAI_K3x4z_S_b = 0.0E0;
  Double I_NAI_K2x5y_S_b = 0.0E0;
  Double I_NAI_K2x4yz_S_b = 0.0E0;
  Double I_NAI_K2x3y2z_S_b = 0.0E0;
  Double I_NAI_K2x2y3z_S_b = 0.0E0;
  Double I_NAI_K2xy4z_S_b = 0.0E0;
  Double I_NAI_K2x5z_S_b = 0.0E0;
  Double I_NAI_Kx6y_S_b = 0.0E0;
  Double I_NAI_Kx5yz_S_b = 0.0E0;
  Double I_NAI_Kx4y2z_S_b = 0.0E0;
  Double I_NAI_Kx3y3z_S_b = 0.0E0;
  Double I_NAI_Kx2y4z_S_b = 0.0E0;
  Double I_NAI_Kxy5z_S_b = 0.0E0;
  Double I_NAI_Kx6z_S_b = 0.0E0;
  Double I_NAI_K7y_S_b = 0.0E0;
  Double I_NAI_K6yz_S_b = 0.0E0;
  Double I_NAI_K5y2z_S_b = 0.0E0;
  Double I_NAI_K4y3z_S_b = 0.0E0;
  Double I_NAI_K3y4z_S_b = 0.0E0;
  Double I_NAI_K2y5z_S_b = 0.0E0;
  Double I_NAI_Ky6z_S_b = 0.0E0;
  Double I_NAI_K7z_S_b = 0.0E0;
  Double I_NAI_I6x_S_b = 0.0E0;
  Double I_NAI_I5xy_S_b = 0.0E0;
  Double I_NAI_I5xz_S_b = 0.0E0;
  Double I_NAI_I4x2y_S_b = 0.0E0;
  Double I_NAI_I4xyz_S_b = 0.0E0;
  Double I_NAI_I4x2z_S_b = 0.0E0;
  Double I_NAI_I3x3y_S_b = 0.0E0;
  Double I_NAI_I3x2yz_S_b = 0.0E0;
  Double I_NAI_I3xy2z_S_b = 0.0E0;
  Double I_NAI_I3x3z_S_b = 0.0E0;
  Double I_NAI_I2x4y_S_b = 0.0E0;
  Double I_NAI_I2x3yz_S_b = 0.0E0;
  Double I_NAI_I2x2y2z_S_b = 0.0E0;
  Double I_NAI_I2xy3z_S_b = 0.0E0;
  Double I_NAI_I2x4z_S_b = 0.0E0;
  Double I_NAI_Ix5y_S_b = 0.0E0;
  Double I_NAI_Ix4yz_S_b = 0.0E0;
  Double I_NAI_Ix3y2z_S_b = 0.0E0;
  Double I_NAI_Ix2y3z_S_b = 0.0E0;
  Double I_NAI_Ixy4z_S_b = 0.0E0;
  Double I_NAI_Ix5z_S_b = 0.0E0;
  Double I_NAI_I6y_S_b = 0.0E0;
  Double I_NAI_I5yz_S_b = 0.0E0;
  Double I_NAI_I4y2z_S_b = 0.0E0;
  Double I_NAI_I3y3z_S_b = 0.0E0;
  Double I_NAI_I2y4z_S_b = 0.0E0;
  Double I_NAI_Iy5z_S_b = 0.0E0;
  Double I_NAI_I6z_S_b = 0.0E0;
  Double I_NAI_H5x_S_b = 0.0E0;
  Double I_NAI_H4xy_S_b = 0.0E0;
  Double I_NAI_H4xz_S_b = 0.0E0;
  Double I_NAI_H3x2y_S_b = 0.0E0;
  Double I_NAI_H3xyz_S_b = 0.0E0;
  Double I_NAI_H3x2z_S_b = 0.0E0;
  Double I_NAI_H2x3y_S_b = 0.0E0;
  Double I_NAI_H2x2yz_S_b = 0.0E0;
  Double I_NAI_H2xy2z_S_b = 0.0E0;
  Double I_NAI_H2x3z_S_b = 0.0E0;
  Double I_NAI_Hx4y_S_b = 0.0E0;
  Double I_NAI_Hx3yz_S_b = 0.0E0;
  Double I_NAI_Hx2y2z_S_b = 0.0E0;
  Double I_NAI_Hxy3z_S_b = 0.0E0;
  Double I_NAI_Hx4z_S_b = 0.0E0;
  Double I_NAI_H5y_S_b = 0.0E0;
  Double I_NAI_H4yz_S_b = 0.0E0;
  Double I_NAI_H3y2z_S_b = 0.0E0;
  Double I_NAI_H2y3z_S_b = 0.0E0;
  Double I_NAI_Hy4z_S_b = 0.0E0;
  Double I_NAI_H5z_S_b = 0.0E0;
  Double I_NAI_G4x_S_b = 0.0E0;
  Double I_NAI_G3xy_S_b = 0.0E0;
  Double I_NAI_G3xz_S_b = 0.0E0;
  Double I_NAI_G2x2y_S_b = 0.0E0;
  Double I_NAI_G2xyz_S_b = 0.0E0;
  Double I_NAI_G2x2z_S_b = 0.0E0;
  Double I_NAI_Gx3y_S_b = 0.0E0;
  Double I_NAI_Gx2yz_S_b = 0.0E0;
  Double I_NAI_Gxy2z_S_b = 0.0E0;
  Double I_NAI_Gx3z_S_b = 0.0E0;
  Double I_NAI_G4y_S_b = 0.0E0;
  Double I_NAI_G3yz_S_b = 0.0E0;
  Double I_NAI_G2y2z_S_b = 0.0E0;
  Double I_NAI_Gy3z_S_b = 0.0E0;
  Double I_NAI_G4z_S_b = 0.0E0;
  Double I_NAI_M9x_S_ab = 0.0E0;
  Double I_NAI_M8xy_S_ab = 0.0E0;
  Double I_NAI_M8xz_S_ab = 0.0E0;
  Double I_NAI_M7x2y_S_ab = 0.0E0;
  Double I_NAI_M7xyz_S_ab = 0.0E0;
  Double I_NAI_M7x2z_S_ab = 0.0E0;
  Double I_NAI_M6x3y_S_ab = 0.0E0;
  Double I_NAI_M6x2yz_S_ab = 0.0E0;
  Double I_NAI_M6xy2z_S_ab = 0.0E0;
  Double I_NAI_M6x3z_S_ab = 0.0E0;
  Double I_NAI_M5x4y_S_ab = 0.0E0;
  Double I_NAI_M5x3yz_S_ab = 0.0E0;
  Double I_NAI_M5x2y2z_S_ab = 0.0E0;
  Double I_NAI_M5xy3z_S_ab = 0.0E0;
  Double I_NAI_M5x4z_S_ab = 0.0E0;
  Double I_NAI_M4x5y_S_ab = 0.0E0;
  Double I_NAI_M4x4yz_S_ab = 0.0E0;
  Double I_NAI_M4x3y2z_S_ab = 0.0E0;
  Double I_NAI_M4x2y3z_S_ab = 0.0E0;
  Double I_NAI_M4xy4z_S_ab = 0.0E0;
  Double I_NAI_M4x5z_S_ab = 0.0E0;
  Double I_NAI_M3x6y_S_ab = 0.0E0;
  Double I_NAI_M3x5yz_S_ab = 0.0E0;
  Double I_NAI_M3x4y2z_S_ab = 0.0E0;
  Double I_NAI_M3x3y3z_S_ab = 0.0E0;
  Double I_NAI_M3x2y4z_S_ab = 0.0E0;
  Double I_NAI_M3xy5z_S_ab = 0.0E0;
  Double I_NAI_M3x6z_S_ab = 0.0E0;
  Double I_NAI_M2x7y_S_ab = 0.0E0;
  Double I_NAI_M2x6yz_S_ab = 0.0E0;
  Double I_NAI_M2x5y2z_S_ab = 0.0E0;
  Double I_NAI_M2x4y3z_S_ab = 0.0E0;
  Double I_NAI_M2x3y4z_S_ab = 0.0E0;
  Double I_NAI_M2x2y5z_S_ab = 0.0E0;
  Double I_NAI_M2xy6z_S_ab = 0.0E0;
  Double I_NAI_M2x7z_S_ab = 0.0E0;
  Double I_NAI_Mx8y_S_ab = 0.0E0;
  Double I_NAI_Mx7yz_S_ab = 0.0E0;
  Double I_NAI_Mx6y2z_S_ab = 0.0E0;
  Double I_NAI_Mx5y3z_S_ab = 0.0E0;
  Double I_NAI_Mx4y4z_S_ab = 0.0E0;
  Double I_NAI_Mx3y5z_S_ab = 0.0E0;
  Double I_NAI_Mx2y6z_S_ab = 0.0E0;
  Double I_NAI_Mxy7z_S_ab = 0.0E0;
  Double I_NAI_Mx8z_S_ab = 0.0E0;
  Double I_NAI_M9y_S_ab = 0.0E0;
  Double I_NAI_M8yz_S_ab = 0.0E0;
  Double I_NAI_M7y2z_S_ab = 0.0E0;
  Double I_NAI_M6y3z_S_ab = 0.0E0;
  Double I_NAI_M5y4z_S_ab = 0.0E0;
  Double I_NAI_M4y5z_S_ab = 0.0E0;
  Double I_NAI_M3y6z_S_ab = 0.0E0;
  Double I_NAI_M2y7z_S_ab = 0.0E0;
  Double I_NAI_My8z_S_ab = 0.0E0;
  Double I_NAI_M9z_S_ab = 0.0E0;
  Double I_NAI_L8x_S_ab = 0.0E0;
  Double I_NAI_L7xy_S_ab = 0.0E0;
  Double I_NAI_L7xz_S_ab = 0.0E0;
  Double I_NAI_L6x2y_S_ab = 0.0E0;
  Double I_NAI_L6xyz_S_ab = 0.0E0;
  Double I_NAI_L6x2z_S_ab = 0.0E0;
  Double I_NAI_L5x3y_S_ab = 0.0E0;
  Double I_NAI_L5x2yz_S_ab = 0.0E0;
  Double I_NAI_L5xy2z_S_ab = 0.0E0;
  Double I_NAI_L5x3z_S_ab = 0.0E0;
  Double I_NAI_L4x4y_S_ab = 0.0E0;
  Double I_NAI_L4x3yz_S_ab = 0.0E0;
  Double I_NAI_L4x2y2z_S_ab = 0.0E0;
  Double I_NAI_L4xy3z_S_ab = 0.0E0;
  Double I_NAI_L4x4z_S_ab = 0.0E0;
  Double I_NAI_L3x5y_S_ab = 0.0E0;
  Double I_NAI_L3x4yz_S_ab = 0.0E0;
  Double I_NAI_L3x3y2z_S_ab = 0.0E0;
  Double I_NAI_L3x2y3z_S_ab = 0.0E0;
  Double I_NAI_L3xy4z_S_ab = 0.0E0;
  Double I_NAI_L3x5z_S_ab = 0.0E0;
  Double I_NAI_L2x6y_S_ab = 0.0E0;
  Double I_NAI_L2x5yz_S_ab = 0.0E0;
  Double I_NAI_L2x4y2z_S_ab = 0.0E0;
  Double I_NAI_L2x3y3z_S_ab = 0.0E0;
  Double I_NAI_L2x2y4z_S_ab = 0.0E0;
  Double I_NAI_L2xy5z_S_ab = 0.0E0;
  Double I_NAI_L2x6z_S_ab = 0.0E0;
  Double I_NAI_Lx7y_S_ab = 0.0E0;
  Double I_NAI_Lx6yz_S_ab = 0.0E0;
  Double I_NAI_Lx5y2z_S_ab = 0.0E0;
  Double I_NAI_Lx4y3z_S_ab = 0.0E0;
  Double I_NAI_Lx3y4z_S_ab = 0.0E0;
  Double I_NAI_Lx2y5z_S_ab = 0.0E0;
  Double I_NAI_Lxy6z_S_ab = 0.0E0;
  Double I_NAI_Lx7z_S_ab = 0.0E0;
  Double I_NAI_L8y_S_ab = 0.0E0;
  Double I_NAI_L7yz_S_ab = 0.0E0;
  Double I_NAI_L6y2z_S_ab = 0.0E0;
  Double I_NAI_L5y3z_S_ab = 0.0E0;
  Double I_NAI_L4y4z_S_ab = 0.0E0;
  Double I_NAI_L3y5z_S_ab = 0.0E0;
  Double I_NAI_L2y6z_S_ab = 0.0E0;
  Double I_NAI_Ly7z_S_ab = 0.0E0;
  Double I_NAI_L8z_S_ab = 0.0E0;
  Double I_NAI_K7x_S_ab = 0.0E0;
  Double I_NAI_K6xy_S_ab = 0.0E0;
  Double I_NAI_K6xz_S_ab = 0.0E0;
  Double I_NAI_K5x2y_S_ab = 0.0E0;
  Double I_NAI_K5xyz_S_ab = 0.0E0;
  Double I_NAI_K5x2z_S_ab = 0.0E0;
  Double I_NAI_K4x3y_S_ab = 0.0E0;
  Double I_NAI_K4x2yz_S_ab = 0.0E0;
  Double I_NAI_K4xy2z_S_ab = 0.0E0;
  Double I_NAI_K4x3z_S_ab = 0.0E0;
  Double I_NAI_K3x4y_S_ab = 0.0E0;
  Double I_NAI_K3x3yz_S_ab = 0.0E0;
  Double I_NAI_K3x2y2z_S_ab = 0.0E0;
  Double I_NAI_K3xy3z_S_ab = 0.0E0;
  Double I_NAI_K3x4z_S_ab = 0.0E0;
  Double I_NAI_K2x5y_S_ab = 0.0E0;
  Double I_NAI_K2x4yz_S_ab = 0.0E0;
  Double I_NAI_K2x3y2z_S_ab = 0.0E0;
  Double I_NAI_K2x2y3z_S_ab = 0.0E0;
  Double I_NAI_K2xy4z_S_ab = 0.0E0;
  Double I_NAI_K2x5z_S_ab = 0.0E0;
  Double I_NAI_Kx6y_S_ab = 0.0E0;
  Double I_NAI_Kx5yz_S_ab = 0.0E0;
  Double I_NAI_Kx4y2z_S_ab = 0.0E0;
  Double I_NAI_Kx3y3z_S_ab = 0.0E0;
  Double I_NAI_Kx2y4z_S_ab = 0.0E0;
  Double I_NAI_Kxy5z_S_ab = 0.0E0;
  Double I_NAI_Kx6z_S_ab = 0.0E0;
  Double I_NAI_K7y_S_ab = 0.0E0;
  Double I_NAI_K6yz_S_ab = 0.0E0;
  Double I_NAI_K5y2z_S_ab = 0.0E0;
  Double I_NAI_K4y3z_S_ab = 0.0E0;
  Double I_NAI_K3y4z_S_ab = 0.0E0;
  Double I_NAI_K2y5z_S_ab = 0.0E0;
  Double I_NAI_Ky6z_S_ab = 0.0E0;
  Double I_NAI_K7z_S_ab = 0.0E0;
  Double I_NAI_I6x_S_ab = 0.0E0;
  Double I_NAI_I5xy_S_ab = 0.0E0;
  Double I_NAI_I5xz_S_ab = 0.0E0;
  Double I_NAI_I4x2y_S_ab = 0.0E0;
  Double I_NAI_I4xyz_S_ab = 0.0E0;
  Double I_NAI_I4x2z_S_ab = 0.0E0;
  Double I_NAI_I3x3y_S_ab = 0.0E0;
  Double I_NAI_I3x2yz_S_ab = 0.0E0;
  Double I_NAI_I3xy2z_S_ab = 0.0E0;
  Double I_NAI_I3x3z_S_ab = 0.0E0;
  Double I_NAI_I2x4y_S_ab = 0.0E0;
  Double I_NAI_I2x3yz_S_ab = 0.0E0;
  Double I_NAI_I2x2y2z_S_ab = 0.0E0;
  Double I_NAI_I2xy3z_S_ab = 0.0E0;
  Double I_NAI_I2x4z_S_ab = 0.0E0;
  Double I_NAI_Ix5y_S_ab = 0.0E0;
  Double I_NAI_Ix4yz_S_ab = 0.0E0;
  Double I_NAI_Ix3y2z_S_ab = 0.0E0;
  Double I_NAI_Ix2y3z_S_ab = 0.0E0;
  Double I_NAI_Ixy4z_S_ab = 0.0E0;
  Double I_NAI_Ix5z_S_ab = 0.0E0;
  Double I_NAI_I6y_S_ab = 0.0E0;
  Double I_NAI_I5yz_S_ab = 0.0E0;
  Double I_NAI_I4y2z_S_ab = 0.0E0;
  Double I_NAI_I3y3z_S_ab = 0.0E0;
  Double I_NAI_I2y4z_S_ab = 0.0E0;
  Double I_NAI_Iy5z_S_ab = 0.0E0;
  Double I_NAI_I6z_S_ab = 0.0E0;
  Double I_NAI_H5x_S_ab = 0.0E0;
  Double I_NAI_H4xy_S_ab = 0.0E0;
  Double I_NAI_H4xz_S_ab = 0.0E0;
  Double I_NAI_H3x2y_S_ab = 0.0E0;
  Double I_NAI_H3xyz_S_ab = 0.0E0;
  Double I_NAI_H3x2z_S_ab = 0.0E0;
  Double I_NAI_H2x3y_S_ab = 0.0E0;
  Double I_NAI_H2x2yz_S_ab = 0.0E0;
  Double I_NAI_H2xy2z_S_ab = 0.0E0;
  Double I_NAI_H2x3z_S_ab = 0.0E0;
  Double I_NAI_Hx4y_S_ab = 0.0E0;
  Double I_NAI_Hx3yz_S_ab = 0.0E0;
  Double I_NAI_Hx2y2z_S_ab = 0.0E0;
  Double I_NAI_Hxy3z_S_ab = 0.0E0;
  Double I_NAI_Hx4z_S_ab = 0.0E0;
  Double I_NAI_H5y_S_ab = 0.0E0;
  Double I_NAI_H4yz_S_ab = 0.0E0;
  Double I_NAI_H3y2z_S_ab = 0.0E0;
  Double I_NAI_H2y3z_S_ab = 0.0E0;
  Double I_NAI_Hy4z_S_ab = 0.0E0;
  Double I_NAI_H5z_S_ab = 0.0E0;
  Double I_NAI_F3x_S_b = 0.0E0;
  Double I_NAI_F2xy_S_b = 0.0E0;
  Double I_NAI_F2xz_S_b = 0.0E0;
  Double I_NAI_Fx2y_S_b = 0.0E0;
  Double I_NAI_Fxyz_S_b = 0.0E0;
  Double I_NAI_Fx2z_S_b = 0.0E0;
  Double I_NAI_F3y_S_b = 0.0E0;
  Double I_NAI_F2yz_S_b = 0.0E0;
  Double I_NAI_Fy2z_S_b = 0.0E0;
  Double I_NAI_F3z_S_b = 0.0E0;
  Double I_NAI_M9x_S_bb = 0.0E0;
  Double I_NAI_M8xy_S_bb = 0.0E0;
  Double I_NAI_M8xz_S_bb = 0.0E0;
  Double I_NAI_M7x2y_S_bb = 0.0E0;
  Double I_NAI_M7xyz_S_bb = 0.0E0;
  Double I_NAI_M7x2z_S_bb = 0.0E0;
  Double I_NAI_M6x3y_S_bb = 0.0E0;
  Double I_NAI_M6x2yz_S_bb = 0.0E0;
  Double I_NAI_M6xy2z_S_bb = 0.0E0;
  Double I_NAI_M6x3z_S_bb = 0.0E0;
  Double I_NAI_M5x4y_S_bb = 0.0E0;
  Double I_NAI_M5x3yz_S_bb = 0.0E0;
  Double I_NAI_M5x2y2z_S_bb = 0.0E0;
  Double I_NAI_M5xy3z_S_bb = 0.0E0;
  Double I_NAI_M5x4z_S_bb = 0.0E0;
  Double I_NAI_M4x5y_S_bb = 0.0E0;
  Double I_NAI_M4x4yz_S_bb = 0.0E0;
  Double I_NAI_M4x3y2z_S_bb = 0.0E0;
  Double I_NAI_M4x2y3z_S_bb = 0.0E0;
  Double I_NAI_M4xy4z_S_bb = 0.0E0;
  Double I_NAI_M4x5z_S_bb = 0.0E0;
  Double I_NAI_M3x6y_S_bb = 0.0E0;
  Double I_NAI_M3x5yz_S_bb = 0.0E0;
  Double I_NAI_M3x4y2z_S_bb = 0.0E0;
  Double I_NAI_M3x3y3z_S_bb = 0.0E0;
  Double I_NAI_M3x2y4z_S_bb = 0.0E0;
  Double I_NAI_M3xy5z_S_bb = 0.0E0;
  Double I_NAI_M3x6z_S_bb = 0.0E0;
  Double I_NAI_M2x7y_S_bb = 0.0E0;
  Double I_NAI_M2x6yz_S_bb = 0.0E0;
  Double I_NAI_M2x5y2z_S_bb = 0.0E0;
  Double I_NAI_M2x4y3z_S_bb = 0.0E0;
  Double I_NAI_M2x3y4z_S_bb = 0.0E0;
  Double I_NAI_M2x2y5z_S_bb = 0.0E0;
  Double I_NAI_M2xy6z_S_bb = 0.0E0;
  Double I_NAI_M2x7z_S_bb = 0.0E0;
  Double I_NAI_Mx8y_S_bb = 0.0E0;
  Double I_NAI_Mx7yz_S_bb = 0.0E0;
  Double I_NAI_Mx6y2z_S_bb = 0.0E0;
  Double I_NAI_Mx5y3z_S_bb = 0.0E0;
  Double I_NAI_Mx4y4z_S_bb = 0.0E0;
  Double I_NAI_Mx3y5z_S_bb = 0.0E0;
  Double I_NAI_Mx2y6z_S_bb = 0.0E0;
  Double I_NAI_Mxy7z_S_bb = 0.0E0;
  Double I_NAI_Mx8z_S_bb = 0.0E0;
  Double I_NAI_M9y_S_bb = 0.0E0;
  Double I_NAI_M8yz_S_bb = 0.0E0;
  Double I_NAI_M7y2z_S_bb = 0.0E0;
  Double I_NAI_M6y3z_S_bb = 0.0E0;
  Double I_NAI_M5y4z_S_bb = 0.0E0;
  Double I_NAI_M4y5z_S_bb = 0.0E0;
  Double I_NAI_M3y6z_S_bb = 0.0E0;
  Double I_NAI_M2y7z_S_bb = 0.0E0;
  Double I_NAI_My8z_S_bb = 0.0E0;
  Double I_NAI_M9z_S_bb = 0.0E0;
  Double I_NAI_L8x_S_bb = 0.0E0;
  Double I_NAI_L7xy_S_bb = 0.0E0;
  Double I_NAI_L7xz_S_bb = 0.0E0;
  Double I_NAI_L6x2y_S_bb = 0.0E0;
  Double I_NAI_L6xyz_S_bb = 0.0E0;
  Double I_NAI_L6x2z_S_bb = 0.0E0;
  Double I_NAI_L5x3y_S_bb = 0.0E0;
  Double I_NAI_L5x2yz_S_bb = 0.0E0;
  Double I_NAI_L5xy2z_S_bb = 0.0E0;
  Double I_NAI_L5x3z_S_bb = 0.0E0;
  Double I_NAI_L4x4y_S_bb = 0.0E0;
  Double I_NAI_L4x3yz_S_bb = 0.0E0;
  Double I_NAI_L4x2y2z_S_bb = 0.0E0;
  Double I_NAI_L4xy3z_S_bb = 0.0E0;
  Double I_NAI_L4x4z_S_bb = 0.0E0;
  Double I_NAI_L3x5y_S_bb = 0.0E0;
  Double I_NAI_L3x4yz_S_bb = 0.0E0;
  Double I_NAI_L3x3y2z_S_bb = 0.0E0;
  Double I_NAI_L3x2y3z_S_bb = 0.0E0;
  Double I_NAI_L3xy4z_S_bb = 0.0E0;
  Double I_NAI_L3x5z_S_bb = 0.0E0;
  Double I_NAI_L2x6y_S_bb = 0.0E0;
  Double I_NAI_L2x5yz_S_bb = 0.0E0;
  Double I_NAI_L2x4y2z_S_bb = 0.0E0;
  Double I_NAI_L2x3y3z_S_bb = 0.0E0;
  Double I_NAI_L2x2y4z_S_bb = 0.0E0;
  Double I_NAI_L2xy5z_S_bb = 0.0E0;
  Double I_NAI_L2x6z_S_bb = 0.0E0;
  Double I_NAI_Lx7y_S_bb = 0.0E0;
  Double I_NAI_Lx6yz_S_bb = 0.0E0;
  Double I_NAI_Lx5y2z_S_bb = 0.0E0;
  Double I_NAI_Lx4y3z_S_bb = 0.0E0;
  Double I_NAI_Lx3y4z_S_bb = 0.0E0;
  Double I_NAI_Lx2y5z_S_bb = 0.0E0;
  Double I_NAI_Lxy6z_S_bb = 0.0E0;
  Double I_NAI_Lx7z_S_bb = 0.0E0;
  Double I_NAI_L8y_S_bb = 0.0E0;
  Double I_NAI_L7yz_S_bb = 0.0E0;
  Double I_NAI_L6y2z_S_bb = 0.0E0;
  Double I_NAI_L5y3z_S_bb = 0.0E0;
  Double I_NAI_L4y4z_S_bb = 0.0E0;
  Double I_NAI_L3y5z_S_bb = 0.0E0;
  Double I_NAI_L2y6z_S_bb = 0.0E0;
  Double I_NAI_Ly7z_S_bb = 0.0E0;
  Double I_NAI_L8z_S_bb = 0.0E0;
  Double I_NAI_K7x_S_bb = 0.0E0;
  Double I_NAI_K6xy_S_bb = 0.0E0;
  Double I_NAI_K6xz_S_bb = 0.0E0;
  Double I_NAI_K5x2y_S_bb = 0.0E0;
  Double I_NAI_K5xyz_S_bb = 0.0E0;
  Double I_NAI_K5x2z_S_bb = 0.0E0;
  Double I_NAI_K4x3y_S_bb = 0.0E0;
  Double I_NAI_K4x2yz_S_bb = 0.0E0;
  Double I_NAI_K4xy2z_S_bb = 0.0E0;
  Double I_NAI_K4x3z_S_bb = 0.0E0;
  Double I_NAI_K3x4y_S_bb = 0.0E0;
  Double I_NAI_K3x3yz_S_bb = 0.0E0;
  Double I_NAI_K3x2y2z_S_bb = 0.0E0;
  Double I_NAI_K3xy3z_S_bb = 0.0E0;
  Double I_NAI_K3x4z_S_bb = 0.0E0;
  Double I_NAI_K2x5y_S_bb = 0.0E0;
  Double I_NAI_K2x4yz_S_bb = 0.0E0;
  Double I_NAI_K2x3y2z_S_bb = 0.0E0;
  Double I_NAI_K2x2y3z_S_bb = 0.0E0;
  Double I_NAI_K2xy4z_S_bb = 0.0E0;
  Double I_NAI_K2x5z_S_bb = 0.0E0;
  Double I_NAI_Kx6y_S_bb = 0.0E0;
  Double I_NAI_Kx5yz_S_bb = 0.0E0;
  Double I_NAI_Kx4y2z_S_bb = 0.0E0;
  Double I_NAI_Kx3y3z_S_bb = 0.0E0;
  Double I_NAI_Kx2y4z_S_bb = 0.0E0;
  Double I_NAI_Kxy5z_S_bb = 0.0E0;
  Double I_NAI_Kx6z_S_bb = 0.0E0;
  Double I_NAI_K7y_S_bb = 0.0E0;
  Double I_NAI_K6yz_S_bb = 0.0E0;
  Double I_NAI_K5y2z_S_bb = 0.0E0;
  Double I_NAI_K4y3z_S_bb = 0.0E0;
  Double I_NAI_K3y4z_S_bb = 0.0E0;
  Double I_NAI_K2y5z_S_bb = 0.0E0;
  Double I_NAI_Ky6z_S_bb = 0.0E0;
  Double I_NAI_K7z_S_bb = 0.0E0;
  Double I_NAI_I6x_S_bb = 0.0E0;
  Double I_NAI_I5xy_S_bb = 0.0E0;
  Double I_NAI_I5xz_S_bb = 0.0E0;
  Double I_NAI_I4x2y_S_bb = 0.0E0;
  Double I_NAI_I4xyz_S_bb = 0.0E0;
  Double I_NAI_I4x2z_S_bb = 0.0E0;
  Double I_NAI_I3x3y_S_bb = 0.0E0;
  Double I_NAI_I3x2yz_S_bb = 0.0E0;
  Double I_NAI_I3xy2z_S_bb = 0.0E0;
  Double I_NAI_I3x3z_S_bb = 0.0E0;
  Double I_NAI_I2x4y_S_bb = 0.0E0;
  Double I_NAI_I2x3yz_S_bb = 0.0E0;
  Double I_NAI_I2x2y2z_S_bb = 0.0E0;
  Double I_NAI_I2xy3z_S_bb = 0.0E0;
  Double I_NAI_I2x4z_S_bb = 0.0E0;
  Double I_NAI_Ix5y_S_bb = 0.0E0;
  Double I_NAI_Ix4yz_S_bb = 0.0E0;
  Double I_NAI_Ix3y2z_S_bb = 0.0E0;
  Double I_NAI_Ix2y3z_S_bb = 0.0E0;
  Double I_NAI_Ixy4z_S_bb = 0.0E0;
  Double I_NAI_Ix5z_S_bb = 0.0E0;
  Double I_NAI_I6y_S_bb = 0.0E0;
  Double I_NAI_I5yz_S_bb = 0.0E0;
  Double I_NAI_I4y2z_S_bb = 0.0E0;
  Double I_NAI_I3y3z_S_bb = 0.0E0;
  Double I_NAI_I2y4z_S_bb = 0.0E0;
  Double I_NAI_Iy5z_S_bb = 0.0E0;
  Double I_NAI_I6z_S_bb = 0.0E0;
  Double I_NAI_H5x_S_bb = 0.0E0;
  Double I_NAI_H4xy_S_bb = 0.0E0;
  Double I_NAI_H4xz_S_bb = 0.0E0;
  Double I_NAI_H3x2y_S_bb = 0.0E0;
  Double I_NAI_H3xyz_S_bb = 0.0E0;
  Double I_NAI_H3x2z_S_bb = 0.0E0;
  Double I_NAI_H2x3y_S_bb = 0.0E0;
  Double I_NAI_H2x2yz_S_bb = 0.0E0;
  Double I_NAI_H2xy2z_S_bb = 0.0E0;
  Double I_NAI_H2x3z_S_bb = 0.0E0;
  Double I_NAI_Hx4y_S_bb = 0.0E0;
  Double I_NAI_Hx3yz_S_bb = 0.0E0;
  Double I_NAI_Hx2y2z_S_bb = 0.0E0;
  Double I_NAI_Hxy3z_S_bb = 0.0E0;
  Double I_NAI_Hx4z_S_bb = 0.0E0;
  Double I_NAI_H5y_S_bb = 0.0E0;
  Double I_NAI_H4yz_S_bb = 0.0E0;
  Double I_NAI_H3y2z_S_bb = 0.0E0;
  Double I_NAI_H2y3z_S_bb = 0.0E0;
  Double I_NAI_Hy4z_S_bb = 0.0E0;
  Double I_NAI_H5z_S_bb = 0.0E0;
  Double I_NAI_G4x_S_bb = 0.0E0;
  Double I_NAI_G3xy_S_bb = 0.0E0;
  Double I_NAI_G3xz_S_bb = 0.0E0;
  Double I_NAI_G2x2y_S_bb = 0.0E0;
  Double I_NAI_G2xyz_S_bb = 0.0E0;
  Double I_NAI_G2x2z_S_bb = 0.0E0;
  Double I_NAI_Gx3y_S_bb = 0.0E0;
  Double I_NAI_Gx2yz_S_bb = 0.0E0;
  Double I_NAI_Gxy2z_S_bb = 0.0E0;
  Double I_NAI_Gx3z_S_bb = 0.0E0;
  Double I_NAI_G4y_S_bb = 0.0E0;
  Double I_NAI_G3yz_S_bb = 0.0E0;
  Double I_NAI_G2y2z_S_bb = 0.0E0;
  Double I_NAI_Gy3z_S_bb = 0.0E0;
  Double I_NAI_G4z_S_bb = 0.0E0;

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
      Double I_NAI_S_S_M9_vrr  = 0.0E0;

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
      double I_NAI_S_S_M9_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER53;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER51*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER49*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER47*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M9_vrr;
        I_NAI_S_S_M9_vrr = ONEOVER19*I_NAI_S_S_M9_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M9_vrr  = f*I_NAI_S_S_M9_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_NAI_S_S_M8_vrr  = ONEOVER17*(u2*I_NAI_S_S_M9_vrr+f);
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
        I_NAI_S_S_M9_vrr_d = oneO2u*(17.0E0*I_NAI_S_S_M8_vrr_d-f);

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
        I_NAI_S_S_M9_vrr = static_cast<Double>(I_NAI_S_S_M9_vrr_d);

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
        I_NAI_S_S_M9_vrr = oneO2u*(17.0E0*I_NAI_S_S_M8_vrr-f);

#endif

      }


      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M8
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M8
       * RHS shell quartet name: SQ_NAI_S_S_M9
       ************************************************************/
      Double I_NAI_Px_S_M8_vrr = PAX*I_NAI_S_S_M8_vrr-PNX*I_NAI_S_S_M9_vrr;
      Double I_NAI_Py_S_M8_vrr = PAY*I_NAI_S_S_M8_vrr-PNY*I_NAI_S_S_M9_vrr;
      Double I_NAI_Pz_S_M8_vrr = PAZ*I_NAI_S_S_M8_vrr-PNZ*I_NAI_S_S_M9_vrr;

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
       * shell quartet name: SQ_NAI_D_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M7
       * RHS shell quartet name: SQ_NAI_P_S_M8
       * RHS shell quartet name: SQ_NAI_S_S_M7
       * RHS shell quartet name: SQ_NAI_S_S_M8
       ************************************************************/
      Double I_NAI_D2x_S_M7_vrr = PAX*I_NAI_Px_S_M7_vrr-PNX*I_NAI_Px_S_M8_vrr+oned2z*I_NAI_S_S_M7_vrr-oned2z*I_NAI_S_S_M8_vrr;
      Double I_NAI_D2y_S_M7_vrr = PAY*I_NAI_Py_S_M7_vrr-PNY*I_NAI_Py_S_M8_vrr+oned2z*I_NAI_S_S_M7_vrr-oned2z*I_NAI_S_S_M8_vrr;
      Double I_NAI_D2z_S_M7_vrr = PAZ*I_NAI_Pz_S_M7_vrr-PNZ*I_NAI_Pz_S_M8_vrr+oned2z*I_NAI_S_S_M7_vrr-oned2z*I_NAI_S_S_M8_vrr;

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
       * shell quartet name: SQ_NAI_F_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M6
       * RHS shell quartet name: SQ_NAI_D_S_M7
       * RHS shell quartet name: SQ_NAI_P_S_M6
       * RHS shell quartet name: SQ_NAI_P_S_M7
       ************************************************************/
      Double I_NAI_F3x_S_M6_vrr = PAX*I_NAI_D2x_S_M6_vrr-PNX*I_NAI_D2x_S_M7_vrr+2*oned2z*I_NAI_Px_S_M6_vrr-2*oned2z*I_NAI_Px_S_M7_vrr;
      Double I_NAI_F3y_S_M6_vrr = PAY*I_NAI_D2y_S_M6_vrr-PNY*I_NAI_D2y_S_M7_vrr+2*oned2z*I_NAI_Py_S_M6_vrr-2*oned2z*I_NAI_Py_S_M7_vrr;
      Double I_NAI_F3z_S_M6_vrr = PAZ*I_NAI_D2z_S_M6_vrr-PNZ*I_NAI_D2z_S_M7_vrr+2*oned2z*I_NAI_Pz_S_M6_vrr-2*oned2z*I_NAI_Pz_S_M7_vrr;

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
       * shell quartet name: SQ_NAI_G_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 11 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M5
       * RHS shell quartet name: SQ_NAI_F_S_M6
       * RHS shell quartet name: SQ_NAI_D_S_M5
       * RHS shell quartet name: SQ_NAI_D_S_M6
       ************************************************************/
      Double I_NAI_G4x_S_M5_vrr = PAX*I_NAI_F3x_S_M5_vrr-PNX*I_NAI_F3x_S_M6_vrr+3*oned2z*I_NAI_D2x_S_M5_vrr-3*oned2z*I_NAI_D2x_S_M6_vrr;
      Double I_NAI_G3xy_S_M5_vrr = PAY*I_NAI_F3x_S_M5_vrr-PNY*I_NAI_F3x_S_M6_vrr;
      Double I_NAI_G4y_S_M5_vrr = PAY*I_NAI_F3y_S_M5_vrr-PNY*I_NAI_F3y_S_M6_vrr+3*oned2z*I_NAI_D2y_S_M5_vrr-3*oned2z*I_NAI_D2y_S_M6_vrr;
      Double I_NAI_G4z_S_M5_vrr = PAZ*I_NAI_F3z_S_M5_vrr-PNZ*I_NAI_F3z_S_M6_vrr+3*oned2z*I_NAI_D2z_S_M5_vrr-3*oned2z*I_NAI_D2z_S_M6_vrr;

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
       * shell quartet name: SQ_NAI_H_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 11 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M4
       * RHS shell quartet name: SQ_NAI_G_S_M5
       * RHS shell quartet name: SQ_NAI_F_S_M4
       * RHS shell quartet name: SQ_NAI_F_S_M5
       ************************************************************/
      Double I_NAI_H5x_S_M4_vrr = PAX*I_NAI_G4x_S_M4_vrr-PNX*I_NAI_G4x_S_M5_vrr+4*oned2z*I_NAI_F3x_S_M4_vrr-4*oned2z*I_NAI_F3x_S_M5_vrr;
      Double I_NAI_H4xy_S_M4_vrr = PAY*I_NAI_G4x_S_M4_vrr-PNY*I_NAI_G4x_S_M5_vrr;
      Double I_NAI_H4xz_S_M4_vrr = PAZ*I_NAI_G4x_S_M4_vrr-PNZ*I_NAI_G4x_S_M5_vrr;
      Double I_NAI_H3x2y_S_M4_vrr = PAY*I_NAI_G3xy_S_M4_vrr-PNY*I_NAI_G3xy_S_M5_vrr+oned2z*I_NAI_F3x_S_M4_vrr-oned2z*I_NAI_F3x_S_M5_vrr;
      Double I_NAI_Hx4y_S_M4_vrr = PAX*I_NAI_G4y_S_M4_vrr-PNX*I_NAI_G4y_S_M5_vrr;
      Double I_NAI_Hx4z_S_M4_vrr = PAX*I_NAI_G4z_S_M4_vrr-PNX*I_NAI_G4z_S_M5_vrr;
      Double I_NAI_H5y_S_M4_vrr = PAY*I_NAI_G4y_S_M4_vrr-PNY*I_NAI_G4y_S_M5_vrr+4*oned2z*I_NAI_F3y_S_M4_vrr-4*oned2z*I_NAI_F3y_S_M5_vrr;
      Double I_NAI_H4yz_S_M4_vrr = PAZ*I_NAI_G4y_S_M4_vrr-PNZ*I_NAI_G4y_S_M5_vrr;
      Double I_NAI_Hy4z_S_M4_vrr = PAY*I_NAI_G4z_S_M4_vrr-PNY*I_NAI_G4z_S_M5_vrr;
      Double I_NAI_H5z_S_M4_vrr = PAZ*I_NAI_G4z_S_M4_vrr-PNZ*I_NAI_G4z_S_M5_vrr+4*oned2z*I_NAI_F3z_S_M4_vrr-4*oned2z*I_NAI_F3z_S_M5_vrr;

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
       * shell quartet name: SQ_NAI_I_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 12 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S_M3
       * RHS shell quartet name: SQ_NAI_H_S_M4
       * RHS shell quartet name: SQ_NAI_G_S_M3
       * RHS shell quartet name: SQ_NAI_G_S_M4
       ************************************************************/
      Double I_NAI_I6x_S_M3_vrr = PAX*I_NAI_H5x_S_M3_vrr-PNX*I_NAI_H5x_S_M4_vrr+5*oned2z*I_NAI_G4x_S_M3_vrr-5*oned2z*I_NAI_G4x_S_M4_vrr;
      Double I_NAI_I5xy_S_M3_vrr = PAY*I_NAI_H5x_S_M3_vrr-PNY*I_NAI_H5x_S_M4_vrr;
      Double I_NAI_I5xz_S_M3_vrr = PAZ*I_NAI_H5x_S_M3_vrr-PNZ*I_NAI_H5x_S_M4_vrr;
      Double I_NAI_I4x2y_S_M3_vrr = PAY*I_NAI_H4xy_S_M3_vrr-PNY*I_NAI_H4xy_S_M4_vrr+oned2z*I_NAI_G4x_S_M3_vrr-oned2z*I_NAI_G4x_S_M4_vrr;
      Double I_NAI_I4x2z_S_M3_vrr = PAZ*I_NAI_H4xz_S_M3_vrr-PNZ*I_NAI_H4xz_S_M4_vrr+oned2z*I_NAI_G4x_S_M3_vrr-oned2z*I_NAI_G4x_S_M4_vrr;
      Double I_NAI_I3x3y_S_M3_vrr = PAY*I_NAI_H3x2y_S_M3_vrr-PNY*I_NAI_H3x2y_S_M4_vrr+2*oned2z*I_NAI_G3xy_S_M3_vrr-2*oned2z*I_NAI_G3xy_S_M4_vrr;
      Double I_NAI_I2x4y_S_M3_vrr = PAX*I_NAI_Hx4y_S_M3_vrr-PNX*I_NAI_Hx4y_S_M4_vrr+oned2z*I_NAI_G4y_S_M3_vrr-oned2z*I_NAI_G4y_S_M4_vrr;
      Double I_NAI_I2x4z_S_M3_vrr = PAX*I_NAI_Hx4z_S_M3_vrr-PNX*I_NAI_Hx4z_S_M4_vrr+oned2z*I_NAI_G4z_S_M3_vrr-oned2z*I_NAI_G4z_S_M4_vrr;
      Double I_NAI_Ix5y_S_M3_vrr = PAX*I_NAI_H5y_S_M3_vrr-PNX*I_NAI_H5y_S_M4_vrr;
      Double I_NAI_Ix5z_S_M3_vrr = PAX*I_NAI_H5z_S_M3_vrr-PNX*I_NAI_H5z_S_M4_vrr;
      Double I_NAI_I6y_S_M3_vrr = PAY*I_NAI_H5y_S_M3_vrr-PNY*I_NAI_H5y_S_M4_vrr+5*oned2z*I_NAI_G4y_S_M3_vrr-5*oned2z*I_NAI_G4y_S_M4_vrr;
      Double I_NAI_I5yz_S_M3_vrr = PAZ*I_NAI_H5y_S_M3_vrr-PNZ*I_NAI_H5y_S_M4_vrr;
      Double I_NAI_I4y2z_S_M3_vrr = PAZ*I_NAI_H4yz_S_M3_vrr-PNZ*I_NAI_H4yz_S_M4_vrr+oned2z*I_NAI_G4y_S_M3_vrr-oned2z*I_NAI_G4y_S_M4_vrr;
      Double I_NAI_I2y4z_S_M3_vrr = PAY*I_NAI_Hy4z_S_M3_vrr-PNY*I_NAI_Hy4z_S_M4_vrr+oned2z*I_NAI_G4z_S_M3_vrr-oned2z*I_NAI_G4z_S_M4_vrr;
      Double I_NAI_Iy5z_S_M3_vrr = PAY*I_NAI_H5z_S_M3_vrr-PNY*I_NAI_H5z_S_M4_vrr;
      Double I_NAI_I6z_S_M3_vrr = PAZ*I_NAI_H5z_S_M3_vrr-PNZ*I_NAI_H5z_S_M4_vrr+5*oned2z*I_NAI_G4z_S_M3_vrr-5*oned2z*I_NAI_G4z_S_M4_vrr;

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
       * shell quartet name: SQ_NAI_K_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 14 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_I_S_M2
       * RHS shell quartet name: SQ_NAI_I_S_M3
       * RHS shell quartet name: SQ_NAI_H_S_M2
       * RHS shell quartet name: SQ_NAI_H_S_M3
       ************************************************************/
      Double I_NAI_K7x_S_M2_vrr = PAX*I_NAI_I6x_S_M2_vrr-PNX*I_NAI_I6x_S_M3_vrr+6*oned2z*I_NAI_H5x_S_M2_vrr-6*oned2z*I_NAI_H5x_S_M3_vrr;
      Double I_NAI_K6xy_S_M2_vrr = PAY*I_NAI_I6x_S_M2_vrr-PNY*I_NAI_I6x_S_M3_vrr;
      Double I_NAI_K6xz_S_M2_vrr = PAZ*I_NAI_I6x_S_M2_vrr-PNZ*I_NAI_I6x_S_M3_vrr;
      Double I_NAI_K5x2y_S_M2_vrr = PAY*I_NAI_I5xy_S_M2_vrr-PNY*I_NAI_I5xy_S_M3_vrr+oned2z*I_NAI_H5x_S_M2_vrr-oned2z*I_NAI_H5x_S_M3_vrr;
      Double I_NAI_K5x2z_S_M2_vrr = PAZ*I_NAI_I5xz_S_M2_vrr-PNZ*I_NAI_I5xz_S_M3_vrr+oned2z*I_NAI_H5x_S_M2_vrr-oned2z*I_NAI_H5x_S_M3_vrr;
      Double I_NAI_K4x3y_S_M2_vrr = PAY*I_NAI_I4x2y_S_M2_vrr-PNY*I_NAI_I4x2y_S_M3_vrr+2*oned2z*I_NAI_H4xy_S_M2_vrr-2*oned2z*I_NAI_H4xy_S_M3_vrr;
      Double I_NAI_K4x3z_S_M2_vrr = PAZ*I_NAI_I4x2z_S_M2_vrr-PNZ*I_NAI_I4x2z_S_M3_vrr+2*oned2z*I_NAI_H4xz_S_M2_vrr-2*oned2z*I_NAI_H4xz_S_M3_vrr;
      Double I_NAI_K3x4y_S_M2_vrr = PAX*I_NAI_I2x4y_S_M2_vrr-PNX*I_NAI_I2x4y_S_M3_vrr+2*oned2z*I_NAI_Hx4y_S_M2_vrr-2*oned2z*I_NAI_Hx4y_S_M3_vrr;
      Double I_NAI_K3x3yz_S_M2_vrr = PAZ*I_NAI_I3x3y_S_M2_vrr-PNZ*I_NAI_I3x3y_S_M3_vrr;
      Double I_NAI_K3x4z_S_M2_vrr = PAX*I_NAI_I2x4z_S_M2_vrr-PNX*I_NAI_I2x4z_S_M3_vrr+2*oned2z*I_NAI_Hx4z_S_M2_vrr-2*oned2z*I_NAI_Hx4z_S_M3_vrr;
      Double I_NAI_K2x5y_S_M2_vrr = PAX*I_NAI_Ix5y_S_M2_vrr-PNX*I_NAI_Ix5y_S_M3_vrr+oned2z*I_NAI_H5y_S_M2_vrr-oned2z*I_NAI_H5y_S_M3_vrr;
      Double I_NAI_K2x5z_S_M2_vrr = PAX*I_NAI_Ix5z_S_M2_vrr-PNX*I_NAI_Ix5z_S_M3_vrr+oned2z*I_NAI_H5z_S_M2_vrr-oned2z*I_NAI_H5z_S_M3_vrr;
      Double I_NAI_Kx6y_S_M2_vrr = PAX*I_NAI_I6y_S_M2_vrr-PNX*I_NAI_I6y_S_M3_vrr;
      Double I_NAI_Kx6z_S_M2_vrr = PAX*I_NAI_I6z_S_M2_vrr-PNX*I_NAI_I6z_S_M3_vrr;
      Double I_NAI_K7y_S_M2_vrr = PAY*I_NAI_I6y_S_M2_vrr-PNY*I_NAI_I6y_S_M3_vrr+6*oned2z*I_NAI_H5y_S_M2_vrr-6*oned2z*I_NAI_H5y_S_M3_vrr;
      Double I_NAI_K6yz_S_M2_vrr = PAZ*I_NAI_I6y_S_M2_vrr-PNZ*I_NAI_I6y_S_M3_vrr;
      Double I_NAI_K5y2z_S_M2_vrr = PAZ*I_NAI_I5yz_S_M2_vrr-PNZ*I_NAI_I5yz_S_M3_vrr+oned2z*I_NAI_H5y_S_M2_vrr-oned2z*I_NAI_H5y_S_M3_vrr;
      Double I_NAI_K4y3z_S_M2_vrr = PAZ*I_NAI_I4y2z_S_M2_vrr-PNZ*I_NAI_I4y2z_S_M3_vrr+2*oned2z*I_NAI_H4yz_S_M2_vrr-2*oned2z*I_NAI_H4yz_S_M3_vrr;
      Double I_NAI_K3y4z_S_M2_vrr = PAY*I_NAI_I2y4z_S_M2_vrr-PNY*I_NAI_I2y4z_S_M3_vrr+2*oned2z*I_NAI_Hy4z_S_M2_vrr-2*oned2z*I_NAI_Hy4z_S_M3_vrr;
      Double I_NAI_K2y5z_S_M2_vrr = PAY*I_NAI_Iy5z_S_M2_vrr-PNY*I_NAI_Iy5z_S_M3_vrr+oned2z*I_NAI_H5z_S_M2_vrr-oned2z*I_NAI_H5z_S_M3_vrr;
      Double I_NAI_Ky6z_S_M2_vrr = PAY*I_NAI_I6z_S_M2_vrr-PNY*I_NAI_I6z_S_M3_vrr;
      Double I_NAI_K7z_S_M2_vrr = PAZ*I_NAI_I6z_S_M2_vrr-PNZ*I_NAI_I6z_S_M3_vrr+6*oned2z*I_NAI_H5z_S_M2_vrr-6*oned2z*I_NAI_H5z_S_M3_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_D2x_S_M1_vrr = PAX*I_NAI_Px_S_M1_vrr-PNX*I_NAI_Px_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_Dxy_S_M1_vrr = PAY*I_NAI_Px_S_M1_vrr-PNY*I_NAI_Px_S_M2_vrr;
      Double I_NAI_D2y_S_M1_vrr = PAY*I_NAI_Py_S_M1_vrr-PNY*I_NAI_Py_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_D2z_S_M1_vrr = PAZ*I_NAI_Pz_S_M1_vrr-PNZ*I_NAI_Pz_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       ************************************************************/
      Double I_NAI_F3x_S_M1_vrr = PAX*I_NAI_D2x_S_M1_vrr-PNX*I_NAI_D2x_S_M2_vrr+2*oned2z*I_NAI_Px_S_M1_vrr-2*oned2z*I_NAI_Px_S_M2_vrr;
      Double I_NAI_F2xy_S_M1_vrr = PAY*I_NAI_D2x_S_M1_vrr-PNY*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_F2xz_S_M1_vrr = PAZ*I_NAI_D2x_S_M1_vrr-PNZ*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Fx2y_S_M1_vrr = PAX*I_NAI_D2y_S_M1_vrr-PNX*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Fx2z_S_M1_vrr = PAX*I_NAI_D2z_S_M1_vrr-PNX*I_NAI_D2z_S_M2_vrr;
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
       * shell quartet name: SQ_NAI_L_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 11 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_K_S_M1
       * RHS shell quartet name: SQ_NAI_K_S_M2
       * RHS shell quartet name: SQ_NAI_I_S_M1
       * RHS shell quartet name: SQ_NAI_I_S_M2
       ************************************************************/
      Double I_NAI_L8x_S_M1_vrr = PAX*I_NAI_K7x_S_M1_vrr-PNX*I_NAI_K7x_S_M2_vrr+7*oned2z*I_NAI_I6x_S_M1_vrr-7*oned2z*I_NAI_I6x_S_M2_vrr;
      Double I_NAI_L7xy_S_M1_vrr = PAY*I_NAI_K7x_S_M1_vrr-PNY*I_NAI_K7x_S_M2_vrr;
      Double I_NAI_L7xz_S_M1_vrr = PAZ*I_NAI_K7x_S_M1_vrr-PNZ*I_NAI_K7x_S_M2_vrr;
      Double I_NAI_L6x2y_S_M1_vrr = PAY*I_NAI_K6xy_S_M1_vrr-PNY*I_NAI_K6xy_S_M2_vrr+oned2z*I_NAI_I6x_S_M1_vrr-oned2z*I_NAI_I6x_S_M2_vrr;
      Double I_NAI_L6x2z_S_M1_vrr = PAZ*I_NAI_K6xz_S_M1_vrr-PNZ*I_NAI_K6xz_S_M2_vrr+oned2z*I_NAI_I6x_S_M1_vrr-oned2z*I_NAI_I6x_S_M2_vrr;
      Double I_NAI_L5x3y_S_M1_vrr = PAY*I_NAI_K5x2y_S_M1_vrr-PNY*I_NAI_K5x2y_S_M2_vrr+2*oned2z*I_NAI_I5xy_S_M1_vrr-2*oned2z*I_NAI_I5xy_S_M2_vrr;
      Double I_NAI_L5x2yz_S_M1_vrr = PAZ*I_NAI_K5x2y_S_M1_vrr-PNZ*I_NAI_K5x2y_S_M2_vrr;
      Double I_NAI_L5x3z_S_M1_vrr = PAZ*I_NAI_K5x2z_S_M1_vrr-PNZ*I_NAI_K5x2z_S_M2_vrr+2*oned2z*I_NAI_I5xz_S_M1_vrr-2*oned2z*I_NAI_I5xz_S_M2_vrr;
      Double I_NAI_L4x4y_S_M1_vrr = PAY*I_NAI_K4x3y_S_M1_vrr-PNY*I_NAI_K4x3y_S_M2_vrr+3*oned2z*I_NAI_I4x2y_S_M1_vrr-3*oned2z*I_NAI_I4x2y_S_M2_vrr;
      Double I_NAI_L4x3yz_S_M1_vrr = PAZ*I_NAI_K4x3y_S_M1_vrr-PNZ*I_NAI_K4x3y_S_M2_vrr;
      Double I_NAI_L4xy3z_S_M1_vrr = PAY*I_NAI_K4x3z_S_M1_vrr-PNY*I_NAI_K4x3z_S_M2_vrr;
      Double I_NAI_L4x4z_S_M1_vrr = PAZ*I_NAI_K4x3z_S_M1_vrr-PNZ*I_NAI_K4x3z_S_M2_vrr+3*oned2z*I_NAI_I4x2z_S_M1_vrr-3*oned2z*I_NAI_I4x2z_S_M2_vrr;
      Double I_NAI_L3x5y_S_M1_vrr = PAX*I_NAI_K2x5y_S_M1_vrr-PNX*I_NAI_K2x5y_S_M2_vrr+2*oned2z*I_NAI_Ix5y_S_M1_vrr-2*oned2z*I_NAI_Ix5y_S_M2_vrr;
      Double I_NAI_L3x4yz_S_M1_vrr = PAZ*I_NAI_K3x4y_S_M1_vrr-PNZ*I_NAI_K3x4y_S_M2_vrr;
      Double I_NAI_L3x3y2z_S_M1_vrr = PAZ*I_NAI_K3x3yz_S_M1_vrr-PNZ*I_NAI_K3x3yz_S_M2_vrr+oned2z*I_NAI_I3x3y_S_M1_vrr-oned2z*I_NAI_I3x3y_S_M2_vrr;
      Double I_NAI_L3xy4z_S_M1_vrr = PAY*I_NAI_K3x4z_S_M1_vrr-PNY*I_NAI_K3x4z_S_M2_vrr;
      Double I_NAI_L3x5z_S_M1_vrr = PAX*I_NAI_K2x5z_S_M1_vrr-PNX*I_NAI_K2x5z_S_M2_vrr+2*oned2z*I_NAI_Ix5z_S_M1_vrr-2*oned2z*I_NAI_Ix5z_S_M2_vrr;
      Double I_NAI_L2x6y_S_M1_vrr = PAX*I_NAI_Kx6y_S_M1_vrr-PNX*I_NAI_Kx6y_S_M2_vrr+oned2z*I_NAI_I6y_S_M1_vrr-oned2z*I_NAI_I6y_S_M2_vrr;
      Double I_NAI_L2x5yz_S_M1_vrr = PAZ*I_NAI_K2x5y_S_M1_vrr-PNZ*I_NAI_K2x5y_S_M2_vrr;
      Double I_NAI_L2xy5z_S_M1_vrr = PAY*I_NAI_K2x5z_S_M1_vrr-PNY*I_NAI_K2x5z_S_M2_vrr;
      Double I_NAI_L2x6z_S_M1_vrr = PAX*I_NAI_Kx6z_S_M1_vrr-PNX*I_NAI_Kx6z_S_M2_vrr+oned2z*I_NAI_I6z_S_M1_vrr-oned2z*I_NAI_I6z_S_M2_vrr;
      Double I_NAI_Lx7y_S_M1_vrr = PAX*I_NAI_K7y_S_M1_vrr-PNX*I_NAI_K7y_S_M2_vrr;
      Double I_NAI_Lx4y3z_S_M1_vrr = PAX*I_NAI_K4y3z_S_M1_vrr-PNX*I_NAI_K4y3z_S_M2_vrr;
      Double I_NAI_Lx3y4z_S_M1_vrr = PAX*I_NAI_K3y4z_S_M1_vrr-PNX*I_NAI_K3y4z_S_M2_vrr;
      Double I_NAI_Lx7z_S_M1_vrr = PAX*I_NAI_K7z_S_M1_vrr-PNX*I_NAI_K7z_S_M2_vrr;
      Double I_NAI_L8y_S_M1_vrr = PAY*I_NAI_K7y_S_M1_vrr-PNY*I_NAI_K7y_S_M2_vrr+7*oned2z*I_NAI_I6y_S_M1_vrr-7*oned2z*I_NAI_I6y_S_M2_vrr;
      Double I_NAI_L7yz_S_M1_vrr = PAZ*I_NAI_K7y_S_M1_vrr-PNZ*I_NAI_K7y_S_M2_vrr;
      Double I_NAI_L6y2z_S_M1_vrr = PAZ*I_NAI_K6yz_S_M1_vrr-PNZ*I_NAI_K6yz_S_M2_vrr+oned2z*I_NAI_I6y_S_M1_vrr-oned2z*I_NAI_I6y_S_M2_vrr;
      Double I_NAI_L5y3z_S_M1_vrr = PAZ*I_NAI_K5y2z_S_M1_vrr-PNZ*I_NAI_K5y2z_S_M2_vrr+2*oned2z*I_NAI_I5yz_S_M1_vrr-2*oned2z*I_NAI_I5yz_S_M2_vrr;
      Double I_NAI_L4y4z_S_M1_vrr = PAZ*I_NAI_K4y3z_S_M1_vrr-PNZ*I_NAI_K4y3z_S_M2_vrr+3*oned2z*I_NAI_I4y2z_S_M1_vrr-3*oned2z*I_NAI_I4y2z_S_M2_vrr;
      Double I_NAI_L3y5z_S_M1_vrr = PAY*I_NAI_K2y5z_S_M1_vrr-PNY*I_NAI_K2y5z_S_M2_vrr+2*oned2z*I_NAI_Iy5z_S_M1_vrr-2*oned2z*I_NAI_Iy5z_S_M2_vrr;
      Double I_NAI_L2y6z_S_M1_vrr = PAY*I_NAI_Ky6z_S_M1_vrr-PNY*I_NAI_Ky6z_S_M2_vrr+oned2z*I_NAI_I6z_S_M1_vrr-oned2z*I_NAI_I6z_S_M2_vrr;
      Double I_NAI_Ly7z_S_M1_vrr = PAY*I_NAI_K7z_S_M1_vrr-PNY*I_NAI_K7z_S_M2_vrr;
      Double I_NAI_L8z_S_M1_vrr = PAZ*I_NAI_K7z_S_M1_vrr-PNZ*I_NAI_K7z_S_M2_vrr+7*oned2z*I_NAI_I6z_S_M1_vrr-7*oned2z*I_NAI_I6z_S_M2_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dxy_S_vrr = PAY*I_NAI_Px_S_vrr-PNY*I_NAI_Px_S_M1_vrr;
      Double I_NAI_Dxz_S_vrr = PAZ*I_NAI_Px_S_vrr-PNZ*I_NAI_Px_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dyz_S_vrr = PAZ*I_NAI_Py_S_vrr-PNZ*I_NAI_Py_S_M1_vrr;
      Double I_NAI_D2z_S_vrr = PAZ*I_NAI_Pz_S_vrr-PNZ*I_NAI_Pz_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       ************************************************************/
      Double I_NAI_F3x_S_vrr = PAX*I_NAI_D2x_S_vrr-PNX*I_NAI_D2x_S_M1_vrr+2*oned2z*I_NAI_Px_S_vrr-2*oned2z*I_NAI_Px_S_M1_vrr;
      Double I_NAI_F2xy_S_vrr = PAY*I_NAI_D2x_S_vrr-PNY*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F2xz_S_vrr = PAZ*I_NAI_D2x_S_vrr-PNZ*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Fx2y_S_vrr = PAX*I_NAI_D2y_S_vrr-PNX*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fxyz_S_vrr = PAZ*I_NAI_Dxy_S_vrr-PNZ*I_NAI_Dxy_S_M1_vrr;
      Double I_NAI_Fx2z_S_vrr = PAX*I_NAI_D2z_S_vrr-PNX*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3y_S_vrr = PAY*I_NAI_D2y_S_vrr-PNY*I_NAI_D2y_S_M1_vrr+2*oned2z*I_NAI_Py_S_vrr-2*oned2z*I_NAI_Py_S_M1_vrr;
      Double I_NAI_F2yz_S_vrr = PAZ*I_NAI_D2y_S_vrr-PNZ*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fy2z_S_vrr = PAY*I_NAI_D2z_S_vrr-PNY*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3z_S_vrr = PAZ*I_NAI_D2z_S_vrr-PNZ*I_NAI_D2z_S_M1_vrr+2*oned2z*I_NAI_Pz_S_vrr-2*oned2z*I_NAI_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       ************************************************************/
      Double I_NAI_G4x_S_vrr = PAX*I_NAI_F3x_S_vrr-PNX*I_NAI_F3x_S_M1_vrr+3*oned2z*I_NAI_D2x_S_vrr-3*oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G3xy_S_vrr = PAY*I_NAI_F3x_S_vrr-PNY*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G3xz_S_vrr = PAZ*I_NAI_F3x_S_vrr-PNZ*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G2x2y_S_vrr = PAY*I_NAI_F2xy_S_vrr-PNY*I_NAI_F2xy_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G2xyz_S_vrr = PAZ*I_NAI_F2xy_S_vrr-PNZ*I_NAI_F2xy_S_M1_vrr;
      Double I_NAI_G2x2z_S_vrr = PAZ*I_NAI_F2xz_S_vrr-PNZ*I_NAI_F2xz_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Gx3y_S_vrr = PAX*I_NAI_F3y_S_vrr-PNX*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_Gx2yz_S_vrr = PAZ*I_NAI_Fx2y_S_vrr-PNZ*I_NAI_Fx2y_S_M1_vrr;
      Double I_NAI_Gxy2z_S_vrr = PAY*I_NAI_Fx2z_S_vrr-PNY*I_NAI_Fx2z_S_M1_vrr;
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
       * shell quartet name: SQ_NAI_M_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_L_S
       * RHS shell quartet name: SQ_NAI_L_S_M1
       * RHS shell quartet name: SQ_NAI_K_S
       * RHS shell quartet name: SQ_NAI_K_S_M1
       ************************************************************/
      Double I_NAI_M9x_S_vrr = PAX*I_NAI_L8x_S_vrr-PNX*I_NAI_L8x_S_M1_vrr+8*oned2z*I_NAI_K7x_S_vrr-8*oned2z*I_NAI_K7x_S_M1_vrr;
      Double I_NAI_M8xy_S_vrr = PAY*I_NAI_L8x_S_vrr-PNY*I_NAI_L8x_S_M1_vrr;
      Double I_NAI_M8xz_S_vrr = PAZ*I_NAI_L8x_S_vrr-PNZ*I_NAI_L8x_S_M1_vrr;
      Double I_NAI_M7x2y_S_vrr = PAY*I_NAI_L7xy_S_vrr-PNY*I_NAI_L7xy_S_M1_vrr+oned2z*I_NAI_K7x_S_vrr-oned2z*I_NAI_K7x_S_M1_vrr;
      Double I_NAI_M7xyz_S_vrr = PAZ*I_NAI_L7xy_S_vrr-PNZ*I_NAI_L7xy_S_M1_vrr;
      Double I_NAI_M7x2z_S_vrr = PAZ*I_NAI_L7xz_S_vrr-PNZ*I_NAI_L7xz_S_M1_vrr+oned2z*I_NAI_K7x_S_vrr-oned2z*I_NAI_K7x_S_M1_vrr;
      Double I_NAI_M6x3y_S_vrr = PAY*I_NAI_L6x2y_S_vrr-PNY*I_NAI_L6x2y_S_M1_vrr+2*oned2z*I_NAI_K6xy_S_vrr-2*oned2z*I_NAI_K6xy_S_M1_vrr;
      Double I_NAI_M6x2yz_S_vrr = PAZ*I_NAI_L6x2y_S_vrr-PNZ*I_NAI_L6x2y_S_M1_vrr;
      Double I_NAI_M6xy2z_S_vrr = PAY*I_NAI_L6x2z_S_vrr-PNY*I_NAI_L6x2z_S_M1_vrr;
      Double I_NAI_M6x3z_S_vrr = PAZ*I_NAI_L6x2z_S_vrr-PNZ*I_NAI_L6x2z_S_M1_vrr+2*oned2z*I_NAI_K6xz_S_vrr-2*oned2z*I_NAI_K6xz_S_M1_vrr;
      Double I_NAI_M5x4y_S_vrr = PAY*I_NAI_L5x3y_S_vrr-PNY*I_NAI_L5x3y_S_M1_vrr+3*oned2z*I_NAI_K5x2y_S_vrr-3*oned2z*I_NAI_K5x2y_S_M1_vrr;
      Double I_NAI_M5x3yz_S_vrr = PAZ*I_NAI_L5x3y_S_vrr-PNZ*I_NAI_L5x3y_S_M1_vrr;
      Double I_NAI_M5x2y2z_S_vrr = PAZ*I_NAI_L5x2yz_S_vrr-PNZ*I_NAI_L5x2yz_S_M1_vrr+oned2z*I_NAI_K5x2y_S_vrr-oned2z*I_NAI_K5x2y_S_M1_vrr;
      Double I_NAI_M5xy3z_S_vrr = PAY*I_NAI_L5x3z_S_vrr-PNY*I_NAI_L5x3z_S_M1_vrr;
      Double I_NAI_M5x4z_S_vrr = PAZ*I_NAI_L5x3z_S_vrr-PNZ*I_NAI_L5x3z_S_M1_vrr+3*oned2z*I_NAI_K5x2z_S_vrr-3*oned2z*I_NAI_K5x2z_S_M1_vrr;
      Double I_NAI_M4x5y_S_vrr = PAX*I_NAI_L3x5y_S_vrr-PNX*I_NAI_L3x5y_S_M1_vrr+3*oned2z*I_NAI_K2x5y_S_vrr-3*oned2z*I_NAI_K2x5y_S_M1_vrr;
      Double I_NAI_M4x4yz_S_vrr = PAZ*I_NAI_L4x4y_S_vrr-PNZ*I_NAI_L4x4y_S_M1_vrr;
      Double I_NAI_M4x3y2z_S_vrr = PAZ*I_NAI_L4x3yz_S_vrr-PNZ*I_NAI_L4x3yz_S_M1_vrr+oned2z*I_NAI_K4x3y_S_vrr-oned2z*I_NAI_K4x3y_S_M1_vrr;
      Double I_NAI_M4x2y3z_S_vrr = PAY*I_NAI_L4xy3z_S_vrr-PNY*I_NAI_L4xy3z_S_M1_vrr+oned2z*I_NAI_K4x3z_S_vrr-oned2z*I_NAI_K4x3z_S_M1_vrr;
      Double I_NAI_M4xy4z_S_vrr = PAY*I_NAI_L4x4z_S_vrr-PNY*I_NAI_L4x4z_S_M1_vrr;
      Double I_NAI_M4x5z_S_vrr = PAX*I_NAI_L3x5z_S_vrr-PNX*I_NAI_L3x5z_S_M1_vrr+3*oned2z*I_NAI_K2x5z_S_vrr-3*oned2z*I_NAI_K2x5z_S_M1_vrr;
      Double I_NAI_M3x6y_S_vrr = PAX*I_NAI_L2x6y_S_vrr-PNX*I_NAI_L2x6y_S_M1_vrr+2*oned2z*I_NAI_Kx6y_S_vrr-2*oned2z*I_NAI_Kx6y_S_M1_vrr;
      Double I_NAI_M3x5yz_S_vrr = PAZ*I_NAI_L3x5y_S_vrr-PNZ*I_NAI_L3x5y_S_M1_vrr;
      Double I_NAI_M3x4y2z_S_vrr = PAZ*I_NAI_L3x4yz_S_vrr-PNZ*I_NAI_L3x4yz_S_M1_vrr+oned2z*I_NAI_K3x4y_S_vrr-oned2z*I_NAI_K3x4y_S_M1_vrr;
      Double I_NAI_M3x3y3z_S_vrr = PAZ*I_NAI_L3x3y2z_S_vrr-PNZ*I_NAI_L3x3y2z_S_M1_vrr+2*oned2z*I_NAI_K3x3yz_S_vrr-2*oned2z*I_NAI_K3x3yz_S_M1_vrr;
      Double I_NAI_M3x2y4z_S_vrr = PAY*I_NAI_L3xy4z_S_vrr-PNY*I_NAI_L3xy4z_S_M1_vrr+oned2z*I_NAI_K3x4z_S_vrr-oned2z*I_NAI_K3x4z_S_M1_vrr;
      Double I_NAI_M3xy5z_S_vrr = PAY*I_NAI_L3x5z_S_vrr-PNY*I_NAI_L3x5z_S_M1_vrr;
      Double I_NAI_M3x6z_S_vrr = PAX*I_NAI_L2x6z_S_vrr-PNX*I_NAI_L2x6z_S_M1_vrr+2*oned2z*I_NAI_Kx6z_S_vrr-2*oned2z*I_NAI_Kx6z_S_M1_vrr;
      Double I_NAI_M2x7y_S_vrr = PAX*I_NAI_Lx7y_S_vrr-PNX*I_NAI_Lx7y_S_M1_vrr+oned2z*I_NAI_K7y_S_vrr-oned2z*I_NAI_K7y_S_M1_vrr;
      Double I_NAI_M2x6yz_S_vrr = PAZ*I_NAI_L2x6y_S_vrr-PNZ*I_NAI_L2x6y_S_M1_vrr;
      Double I_NAI_M2x5y2z_S_vrr = PAZ*I_NAI_L2x5yz_S_vrr-PNZ*I_NAI_L2x5yz_S_M1_vrr+oned2z*I_NAI_K2x5y_S_vrr-oned2z*I_NAI_K2x5y_S_M1_vrr;
      Double I_NAI_M2x4y3z_S_vrr = PAX*I_NAI_Lx4y3z_S_vrr-PNX*I_NAI_Lx4y3z_S_M1_vrr+oned2z*I_NAI_K4y3z_S_vrr-oned2z*I_NAI_K4y3z_S_M1_vrr;
      Double I_NAI_M2x3y4z_S_vrr = PAX*I_NAI_Lx3y4z_S_vrr-PNX*I_NAI_Lx3y4z_S_M1_vrr+oned2z*I_NAI_K3y4z_S_vrr-oned2z*I_NAI_K3y4z_S_M1_vrr;
      Double I_NAI_M2x2y5z_S_vrr = PAY*I_NAI_L2xy5z_S_vrr-PNY*I_NAI_L2xy5z_S_M1_vrr+oned2z*I_NAI_K2x5z_S_vrr-oned2z*I_NAI_K2x5z_S_M1_vrr;
      Double I_NAI_M2xy6z_S_vrr = PAY*I_NAI_L2x6z_S_vrr-PNY*I_NAI_L2x6z_S_M1_vrr;
      Double I_NAI_M2x7z_S_vrr = PAX*I_NAI_Lx7z_S_vrr-PNX*I_NAI_Lx7z_S_M1_vrr+oned2z*I_NAI_K7z_S_vrr-oned2z*I_NAI_K7z_S_M1_vrr;
      Double I_NAI_Mx8y_S_vrr = PAX*I_NAI_L8y_S_vrr-PNX*I_NAI_L8y_S_M1_vrr;
      Double I_NAI_Mx7yz_S_vrr = PAZ*I_NAI_Lx7y_S_vrr-PNZ*I_NAI_Lx7y_S_M1_vrr;
      Double I_NAI_Mx6y2z_S_vrr = PAX*I_NAI_L6y2z_S_vrr-PNX*I_NAI_L6y2z_S_M1_vrr;
      Double I_NAI_Mx5y3z_S_vrr = PAX*I_NAI_L5y3z_S_vrr-PNX*I_NAI_L5y3z_S_M1_vrr;
      Double I_NAI_Mx4y4z_S_vrr = PAX*I_NAI_L4y4z_S_vrr-PNX*I_NAI_L4y4z_S_M1_vrr;
      Double I_NAI_Mx3y5z_S_vrr = PAX*I_NAI_L3y5z_S_vrr-PNX*I_NAI_L3y5z_S_M1_vrr;
      Double I_NAI_Mx2y6z_S_vrr = PAX*I_NAI_L2y6z_S_vrr-PNX*I_NAI_L2y6z_S_M1_vrr;
      Double I_NAI_Mxy7z_S_vrr = PAY*I_NAI_Lx7z_S_vrr-PNY*I_NAI_Lx7z_S_M1_vrr;
      Double I_NAI_Mx8z_S_vrr = PAX*I_NAI_L8z_S_vrr-PNX*I_NAI_L8z_S_M1_vrr;
      Double I_NAI_M9y_S_vrr = PAY*I_NAI_L8y_S_vrr-PNY*I_NAI_L8y_S_M1_vrr+8*oned2z*I_NAI_K7y_S_vrr-8*oned2z*I_NAI_K7y_S_M1_vrr;
      Double I_NAI_M8yz_S_vrr = PAZ*I_NAI_L8y_S_vrr-PNZ*I_NAI_L8y_S_M1_vrr;
      Double I_NAI_M7y2z_S_vrr = PAZ*I_NAI_L7yz_S_vrr-PNZ*I_NAI_L7yz_S_M1_vrr+oned2z*I_NAI_K7y_S_vrr-oned2z*I_NAI_K7y_S_M1_vrr;
      Double I_NAI_M6y3z_S_vrr = PAZ*I_NAI_L6y2z_S_vrr-PNZ*I_NAI_L6y2z_S_M1_vrr+2*oned2z*I_NAI_K6yz_S_vrr-2*oned2z*I_NAI_K6yz_S_M1_vrr;
      Double I_NAI_M5y4z_S_vrr = PAZ*I_NAI_L5y3z_S_vrr-PNZ*I_NAI_L5y3z_S_M1_vrr+3*oned2z*I_NAI_K5y2z_S_vrr-3*oned2z*I_NAI_K5y2z_S_M1_vrr;
      Double I_NAI_M4y5z_S_vrr = PAY*I_NAI_L3y5z_S_vrr-PNY*I_NAI_L3y5z_S_M1_vrr+3*oned2z*I_NAI_K2y5z_S_vrr-3*oned2z*I_NAI_K2y5z_S_M1_vrr;
      Double I_NAI_M3y6z_S_vrr = PAY*I_NAI_L2y6z_S_vrr-PNY*I_NAI_L2y6z_S_M1_vrr+2*oned2z*I_NAI_Ky6z_S_vrr-2*oned2z*I_NAI_Ky6z_S_M1_vrr;
      Double I_NAI_M2y7z_S_vrr = PAY*I_NAI_Ly7z_S_vrr-PNY*I_NAI_Ly7z_S_M1_vrr+oned2z*I_NAI_K7z_S_vrr-oned2z*I_NAI_K7z_S_M1_vrr;
      Double I_NAI_My8z_S_vrr = PAY*I_NAI_L8z_S_vrr-PNY*I_NAI_L8z_S_M1_vrr;
      Double I_NAI_M9z_S_vrr = PAZ*I_NAI_L8z_S_vrr-PNZ*I_NAI_L8z_S_M1_vrr+8*oned2z*I_NAI_K7z_S_vrr-8*oned2z*I_NAI_K7z_S_M1_vrr;

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

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_G4x_S += I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S += I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S += I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S += I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S += I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S += I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S += I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S += I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S += I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S += I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S += I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S += I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S += I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S += I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S += I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_a_coefs = alpha;
      I_NAI_K7x_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_a_coefs = alpha;
      I_NAI_I6x_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_a_coefs = alpha;
      I_NAI_H5x_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_F3x_S += I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S += I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S += I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S += I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S += I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S += I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S += I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S += I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S += I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S += I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_M_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_M_S_aa_coefs = alpha*alpha;
      I_NAI_M9x_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M9x_S_vrr;
      I_NAI_M8xy_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M8xy_S_vrr;
      I_NAI_M8xz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M8xz_S_vrr;
      I_NAI_M7x2y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M7x2y_S_vrr;
      I_NAI_M7xyz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M7xyz_S_vrr;
      I_NAI_M7x2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M7x2z_S_vrr;
      I_NAI_M6x3y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M6x3y_S_vrr;
      I_NAI_M6x2yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M6x2yz_S_vrr;
      I_NAI_M6xy2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M6xy2z_S_vrr;
      I_NAI_M6x3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M6x3z_S_vrr;
      I_NAI_M5x4y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5x4y_S_vrr;
      I_NAI_M5x3yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5x3yz_S_vrr;
      I_NAI_M5x2y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5x2y2z_S_vrr;
      I_NAI_M5xy3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5xy3z_S_vrr;
      I_NAI_M5x4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5x4z_S_vrr;
      I_NAI_M4x5y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4x5y_S_vrr;
      I_NAI_M4x4yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4x4yz_S_vrr;
      I_NAI_M4x3y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4x3y2z_S_vrr;
      I_NAI_M4x2y3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4x2y3z_S_vrr;
      I_NAI_M4xy4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4xy4z_S_vrr;
      I_NAI_M4x5z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4x5z_S_vrr;
      I_NAI_M3x6y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x6y_S_vrr;
      I_NAI_M3x5yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x5yz_S_vrr;
      I_NAI_M3x4y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x4y2z_S_vrr;
      I_NAI_M3x3y3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x3y3z_S_vrr;
      I_NAI_M3x2y4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x2y4z_S_vrr;
      I_NAI_M3xy5z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3xy5z_S_vrr;
      I_NAI_M3x6z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3x6z_S_vrr;
      I_NAI_M2x7y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x7y_S_vrr;
      I_NAI_M2x6yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x6yz_S_vrr;
      I_NAI_M2x5y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x5y2z_S_vrr;
      I_NAI_M2x4y3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x4y3z_S_vrr;
      I_NAI_M2x3y4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x3y4z_S_vrr;
      I_NAI_M2x2y5z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x2y5z_S_vrr;
      I_NAI_M2xy6z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2xy6z_S_vrr;
      I_NAI_M2x7z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2x7z_S_vrr;
      I_NAI_Mx8y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx8y_S_vrr;
      I_NAI_Mx7yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx7yz_S_vrr;
      I_NAI_Mx6y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx6y2z_S_vrr;
      I_NAI_Mx5y3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx5y3z_S_vrr;
      I_NAI_Mx4y4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx4y4z_S_vrr;
      I_NAI_Mx3y5z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx3y5z_S_vrr;
      I_NAI_Mx2y6z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx2y6z_S_vrr;
      I_NAI_Mxy7z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mxy7z_S_vrr;
      I_NAI_Mx8z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_Mx8z_S_vrr;
      I_NAI_M9y_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M9y_S_vrr;
      I_NAI_M8yz_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M8yz_S_vrr;
      I_NAI_M7y2z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M7y2z_S_vrr;
      I_NAI_M6y3z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M6y3z_S_vrr;
      I_NAI_M5y4z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M5y4z_S_vrr;
      I_NAI_M4y5z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M4y5z_S_vrr;
      I_NAI_M3y6z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M3y6z_S_vrr;
      I_NAI_M2y7z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M2y7z_S_vrr;
      I_NAI_My8z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_My8z_S_vrr;
      I_NAI_M9z_S_aa += SQ_NAI_M_S_aa_coefs*I_NAI_M9z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_L_S_aa_coefs = alpha*alpha;
      I_NAI_L8x_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S_aa += SQ_NAI_L_S_aa_coefs*I_NAI_L8z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_aa_coefs = alpha*alpha;
      I_NAI_K7x_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_aa += SQ_NAI_K_S_aa_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_aa_coefs = alpha*alpha;
      I_NAI_I6x_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_aa += SQ_NAI_I_S_aa_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_a_coefs = alpha;
      I_NAI_G4x_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_D2x_S += I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S += I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S += I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S += I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S += I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S += I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_b_coefs = beta;
      I_NAI_K7x_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_b_coefs = beta;
      I_NAI_I6x_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_b_coefs = beta;
      I_NAI_H5x_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_b_coefs = beta;
      I_NAI_G4x_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_M_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_M_S_ab_coefs = alpha*beta;
      I_NAI_M9x_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M9x_S_vrr;
      I_NAI_M8xy_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M8xy_S_vrr;
      I_NAI_M8xz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M8xz_S_vrr;
      I_NAI_M7x2y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M7x2y_S_vrr;
      I_NAI_M7xyz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M7xyz_S_vrr;
      I_NAI_M7x2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M7x2z_S_vrr;
      I_NAI_M6x3y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M6x3y_S_vrr;
      I_NAI_M6x2yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M6x2yz_S_vrr;
      I_NAI_M6xy2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M6xy2z_S_vrr;
      I_NAI_M6x3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M6x3z_S_vrr;
      I_NAI_M5x4y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5x4y_S_vrr;
      I_NAI_M5x3yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5x3yz_S_vrr;
      I_NAI_M5x2y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5x2y2z_S_vrr;
      I_NAI_M5xy3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5xy3z_S_vrr;
      I_NAI_M5x4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5x4z_S_vrr;
      I_NAI_M4x5y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4x5y_S_vrr;
      I_NAI_M4x4yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4x4yz_S_vrr;
      I_NAI_M4x3y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4x3y2z_S_vrr;
      I_NAI_M4x2y3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4x2y3z_S_vrr;
      I_NAI_M4xy4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4xy4z_S_vrr;
      I_NAI_M4x5z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4x5z_S_vrr;
      I_NAI_M3x6y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x6y_S_vrr;
      I_NAI_M3x5yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x5yz_S_vrr;
      I_NAI_M3x4y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x4y2z_S_vrr;
      I_NAI_M3x3y3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x3y3z_S_vrr;
      I_NAI_M3x2y4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x2y4z_S_vrr;
      I_NAI_M3xy5z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3xy5z_S_vrr;
      I_NAI_M3x6z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3x6z_S_vrr;
      I_NAI_M2x7y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x7y_S_vrr;
      I_NAI_M2x6yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x6yz_S_vrr;
      I_NAI_M2x5y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x5y2z_S_vrr;
      I_NAI_M2x4y3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x4y3z_S_vrr;
      I_NAI_M2x3y4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x3y4z_S_vrr;
      I_NAI_M2x2y5z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x2y5z_S_vrr;
      I_NAI_M2xy6z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2xy6z_S_vrr;
      I_NAI_M2x7z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2x7z_S_vrr;
      I_NAI_Mx8y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx8y_S_vrr;
      I_NAI_Mx7yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx7yz_S_vrr;
      I_NAI_Mx6y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx6y2z_S_vrr;
      I_NAI_Mx5y3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx5y3z_S_vrr;
      I_NAI_Mx4y4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx4y4z_S_vrr;
      I_NAI_Mx3y5z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx3y5z_S_vrr;
      I_NAI_Mx2y6z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx2y6z_S_vrr;
      I_NAI_Mxy7z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mxy7z_S_vrr;
      I_NAI_Mx8z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_Mx8z_S_vrr;
      I_NAI_M9y_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M9y_S_vrr;
      I_NAI_M8yz_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M8yz_S_vrr;
      I_NAI_M7y2z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M7y2z_S_vrr;
      I_NAI_M6y3z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M6y3z_S_vrr;
      I_NAI_M5y4z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M5y4z_S_vrr;
      I_NAI_M4y5z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M4y5z_S_vrr;
      I_NAI_M3y6z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M3y6z_S_vrr;
      I_NAI_M2y7z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M2y7z_S_vrr;
      I_NAI_My8z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_My8z_S_vrr;
      I_NAI_M9z_S_ab += SQ_NAI_M_S_ab_coefs*I_NAI_M9z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_L_S_ab_coefs = alpha*beta;
      I_NAI_L8x_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S_ab += SQ_NAI_L_S_ab_coefs*I_NAI_L8z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_ab_coefs = alpha*beta;
      I_NAI_K7x_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_ab += SQ_NAI_K_S_ab_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_ab_coefs = alpha*beta;
      I_NAI_I6x_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_ab += SQ_NAI_I_S_ab_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_ab_coefs = alpha*beta;
      I_NAI_H5x_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_ab += SQ_NAI_H_S_ab_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_b_coefs = beta;
      I_NAI_F3x_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_M_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_M_S_bb_coefs = beta*beta;
      I_NAI_M9x_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M9x_S_vrr;
      I_NAI_M8xy_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M8xy_S_vrr;
      I_NAI_M8xz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M8xz_S_vrr;
      I_NAI_M7x2y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M7x2y_S_vrr;
      I_NAI_M7xyz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M7xyz_S_vrr;
      I_NAI_M7x2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M7x2z_S_vrr;
      I_NAI_M6x3y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M6x3y_S_vrr;
      I_NAI_M6x2yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M6x2yz_S_vrr;
      I_NAI_M6xy2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M6xy2z_S_vrr;
      I_NAI_M6x3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M6x3z_S_vrr;
      I_NAI_M5x4y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5x4y_S_vrr;
      I_NAI_M5x3yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5x3yz_S_vrr;
      I_NAI_M5x2y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5x2y2z_S_vrr;
      I_NAI_M5xy3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5xy3z_S_vrr;
      I_NAI_M5x4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5x4z_S_vrr;
      I_NAI_M4x5y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4x5y_S_vrr;
      I_NAI_M4x4yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4x4yz_S_vrr;
      I_NAI_M4x3y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4x3y2z_S_vrr;
      I_NAI_M4x2y3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4x2y3z_S_vrr;
      I_NAI_M4xy4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4xy4z_S_vrr;
      I_NAI_M4x5z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4x5z_S_vrr;
      I_NAI_M3x6y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x6y_S_vrr;
      I_NAI_M3x5yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x5yz_S_vrr;
      I_NAI_M3x4y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x4y2z_S_vrr;
      I_NAI_M3x3y3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x3y3z_S_vrr;
      I_NAI_M3x2y4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x2y4z_S_vrr;
      I_NAI_M3xy5z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3xy5z_S_vrr;
      I_NAI_M3x6z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3x6z_S_vrr;
      I_NAI_M2x7y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x7y_S_vrr;
      I_NAI_M2x6yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x6yz_S_vrr;
      I_NAI_M2x5y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x5y2z_S_vrr;
      I_NAI_M2x4y3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x4y3z_S_vrr;
      I_NAI_M2x3y4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x3y4z_S_vrr;
      I_NAI_M2x2y5z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x2y5z_S_vrr;
      I_NAI_M2xy6z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2xy6z_S_vrr;
      I_NAI_M2x7z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2x7z_S_vrr;
      I_NAI_Mx8y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx8y_S_vrr;
      I_NAI_Mx7yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx7yz_S_vrr;
      I_NAI_Mx6y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx6y2z_S_vrr;
      I_NAI_Mx5y3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx5y3z_S_vrr;
      I_NAI_Mx4y4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx4y4z_S_vrr;
      I_NAI_Mx3y5z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx3y5z_S_vrr;
      I_NAI_Mx2y6z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx2y6z_S_vrr;
      I_NAI_Mxy7z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mxy7z_S_vrr;
      I_NAI_Mx8z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_Mx8z_S_vrr;
      I_NAI_M9y_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M9y_S_vrr;
      I_NAI_M8yz_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M8yz_S_vrr;
      I_NAI_M7y2z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M7y2z_S_vrr;
      I_NAI_M6y3z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M6y3z_S_vrr;
      I_NAI_M5y4z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M5y4z_S_vrr;
      I_NAI_M4y5z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M4y5z_S_vrr;
      I_NAI_M3y6z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M3y6z_S_vrr;
      I_NAI_M2y7z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M2y7z_S_vrr;
      I_NAI_My8z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_My8z_S_vrr;
      I_NAI_M9z_S_bb += SQ_NAI_M_S_bb_coefs*I_NAI_M9z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_L_S_bb_coefs = beta*beta;
      I_NAI_L8x_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S_bb += SQ_NAI_L_S_bb_coefs*I_NAI_L8z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_bb_coefs = beta*beta;
      I_NAI_K7x_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_bb += SQ_NAI_K_S_bb_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_bb_coefs = beta*beta;
      I_NAI_I6x_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_bb += SQ_NAI_I_S_bb_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_bb_coefs = beta*beta;
      I_NAI_H5x_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_bb += SQ_NAI_H_S_bb_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_bb_coefs = beta*beta;
      I_NAI_G4x_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_bb += SQ_NAI_G_S_bb_coefs*I_NAI_G4z_S_vrr;
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
   * shell quartet name: SQ_NAI_D_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  Double I_NAI_D2x_Px = I_NAI_F3x_S+ABX*I_NAI_D2x_S;
  Double I_NAI_Dxy_Px = I_NAI_F2xy_S+ABX*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Px = I_NAI_F2xz_S+ABX*I_NAI_Dxz_S;
  Double I_NAI_D2y_Px = I_NAI_Fx2y_S+ABX*I_NAI_D2y_S;
  Double I_NAI_Dyz_Px = I_NAI_Fxyz_S+ABX*I_NAI_Dyz_S;
  Double I_NAI_D2z_Px = I_NAI_Fx2z_S+ABX*I_NAI_D2z_S;
  Double I_NAI_D2x_Py = I_NAI_F2xy_S+ABY*I_NAI_D2x_S;
  Double I_NAI_Dxy_Py = I_NAI_Fx2y_S+ABY*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Py = I_NAI_Fxyz_S+ABY*I_NAI_Dxz_S;
  Double I_NAI_D2y_Py = I_NAI_F3y_S+ABY*I_NAI_D2y_S;
  Double I_NAI_Dyz_Py = I_NAI_F2yz_S+ABY*I_NAI_Dyz_S;
  Double I_NAI_D2z_Py = I_NAI_Fy2z_S+ABY*I_NAI_D2z_S;
  Double I_NAI_D2x_Pz = I_NAI_F2xz_S+ABZ*I_NAI_D2x_S;
  Double I_NAI_Dxy_Pz = I_NAI_Fxyz_S+ABZ*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Pz = I_NAI_Fx2z_S+ABZ*I_NAI_Dxz_S;
  Double I_NAI_D2y_Pz = I_NAI_F2yz_S+ABZ*I_NAI_D2y_S;
  Double I_NAI_Dyz_Pz = I_NAI_Fy2z_S+ABZ*I_NAI_Dyz_S;
  Double I_NAI_D2z_Pz = I_NAI_F3z_S+ABZ*I_NAI_D2z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S
   * RHS shell quartet name: SQ_NAI_F_S
   ************************************************************/
  Double I_NAI_F3x_Px = I_NAI_G4x_S+ABX*I_NAI_F3x_S;
  Double I_NAI_F2xy_Px = I_NAI_G3xy_S+ABX*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Px = I_NAI_G3xz_S+ABX*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Px = I_NAI_G2x2y_S+ABX*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Px = I_NAI_G2xyz_S+ABX*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Px = I_NAI_G2x2z_S+ABX*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Px = I_NAI_Gx3y_S+ABX*I_NAI_F3y_S;
  Double I_NAI_F2yz_Px = I_NAI_Gx2yz_S+ABX*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Px = I_NAI_Gxy2z_S+ABX*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Px = I_NAI_Gx3z_S+ABX*I_NAI_F3z_S;
  Double I_NAI_F3x_Py = I_NAI_G3xy_S+ABY*I_NAI_F3x_S;
  Double I_NAI_F2xy_Py = I_NAI_G2x2y_S+ABY*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Py = I_NAI_G2xyz_S+ABY*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Py = I_NAI_Gx3y_S+ABY*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Py = I_NAI_Gx2yz_S+ABY*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Py = I_NAI_Gxy2z_S+ABY*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Py = I_NAI_G4y_S+ABY*I_NAI_F3y_S;
  Double I_NAI_F2yz_Py = I_NAI_G3yz_S+ABY*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Py = I_NAI_G2y2z_S+ABY*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Py = I_NAI_Gy3z_S+ABY*I_NAI_F3z_S;
  Double I_NAI_F3x_Pz = I_NAI_G3xz_S+ABZ*I_NAI_F3x_S;
  Double I_NAI_F2xy_Pz = I_NAI_G2xyz_S+ABZ*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Pz = I_NAI_G2x2z_S+ABZ*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Pz = I_NAI_Gx2yz_S+ABZ*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Pz = I_NAI_Gxy2z_S+ABZ*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Pz = I_NAI_Gx3z_S+ABZ*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Pz = I_NAI_G3yz_S+ABZ*I_NAI_F3y_S;
  Double I_NAI_F2yz_Pz = I_NAI_G2y2z_S+ABZ*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Pz = I_NAI_Gy3z_S+ABZ*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Pz = I_NAI_G4z_S+ABZ*I_NAI_F3z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P
   * RHS shell quartet name: SQ_NAI_D_P
   ************************************************************/
  Double I_NAI_D2x_D2x = I_NAI_F3x_Px+ABX*I_NAI_D2x_Px;
  Double I_NAI_Dxy_D2x = I_NAI_F2xy_Px+ABX*I_NAI_Dxy_Px;
  Double I_NAI_Dxz_D2x = I_NAI_F2xz_Px+ABX*I_NAI_Dxz_Px;
  Double I_NAI_D2y_D2x = I_NAI_Fx2y_Px+ABX*I_NAI_D2y_Px;
  Double I_NAI_Dyz_D2x = I_NAI_Fxyz_Px+ABX*I_NAI_Dyz_Px;
  Double I_NAI_D2z_D2x = I_NAI_Fx2z_Px+ABX*I_NAI_D2z_Px;
  Double I_NAI_D2x_Dxy = I_NAI_F2xy_Px+ABY*I_NAI_D2x_Px;
  Double I_NAI_Dxy_Dxy = I_NAI_Fx2y_Px+ABY*I_NAI_Dxy_Px;
  Double I_NAI_Dxz_Dxy = I_NAI_Fxyz_Px+ABY*I_NAI_Dxz_Px;
  Double I_NAI_D2y_Dxy = I_NAI_F3y_Px+ABY*I_NAI_D2y_Px;
  Double I_NAI_Dyz_Dxy = I_NAI_F2yz_Px+ABY*I_NAI_Dyz_Px;
  Double I_NAI_D2z_Dxy = I_NAI_Fy2z_Px+ABY*I_NAI_D2z_Px;
  Double I_NAI_D2x_D2y = I_NAI_F2xy_Py+ABY*I_NAI_D2x_Py;
  Double I_NAI_Dxy_D2y = I_NAI_Fx2y_Py+ABY*I_NAI_Dxy_Py;
  Double I_NAI_Dxz_D2y = I_NAI_Fxyz_Py+ABY*I_NAI_Dxz_Py;
  Double I_NAI_D2y_D2y = I_NAI_F3y_Py+ABY*I_NAI_D2y_Py;
  Double I_NAI_Dyz_D2y = I_NAI_F2yz_Py+ABY*I_NAI_Dyz_Py;
  Double I_NAI_D2z_D2y = I_NAI_Fy2z_Py+ABY*I_NAI_D2z_Py;
  Double I_NAI_D2x_D2z = I_NAI_F2xz_Pz+ABZ*I_NAI_D2x_Pz;
  Double I_NAI_Dxy_D2z = I_NAI_Fxyz_Pz+ABZ*I_NAI_Dxy_Pz;
  Double I_NAI_Dxz_D2z = I_NAI_Fx2z_Pz+ABZ*I_NAI_Dxz_Pz;
  Double I_NAI_D2y_D2z = I_NAI_F2yz_Pz+ABZ*I_NAI_D2y_Pz;
  Double I_NAI_Dyz_D2z = I_NAI_Fy2z_Pz+ABZ*I_NAI_Dyz_Pz;
  Double I_NAI_D2z_D2z = I_NAI_F3z_Pz+ABZ*I_NAI_D2z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  Double I_NAI_G4x_Px = I_NAI_H5x_S+ABX*I_NAI_G4x_S;
  Double I_NAI_G3xy_Px = I_NAI_H4xy_S+ABX*I_NAI_G3xy_S;
  Double I_NAI_G3xz_Px = I_NAI_H4xz_S+ABX*I_NAI_G3xz_S;
  Double I_NAI_G2x2y_Px = I_NAI_H3x2y_S+ABX*I_NAI_G2x2y_S;
  Double I_NAI_G2xyz_Px = I_NAI_H3xyz_S+ABX*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Px = I_NAI_H3x2z_S+ABX*I_NAI_G2x2z_S;
  Double I_NAI_Gx3y_Px = I_NAI_H2x3y_S+ABX*I_NAI_Gx3y_S;
  Double I_NAI_Gx2yz_Px = I_NAI_H2x2yz_S+ABX*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Px = I_NAI_H2xy2z_S+ABX*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Px = I_NAI_H2x3z_S+ABX*I_NAI_Gx3z_S;
  Double I_NAI_G4y_Px = I_NAI_Hx4y_S+ABX*I_NAI_G4y_S;
  Double I_NAI_G3yz_Px = I_NAI_Hx3yz_S+ABX*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Px = I_NAI_Hx2y2z_S+ABX*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Px = I_NAI_Hxy3z_S+ABX*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Px = I_NAI_Hx4z_S+ABX*I_NAI_G4z_S;
  Double I_NAI_G4x_Py = I_NAI_H4xy_S+ABY*I_NAI_G4x_S;
  Double I_NAI_G3xy_Py = I_NAI_H3x2y_S+ABY*I_NAI_G3xy_S;
  Double I_NAI_G3xz_Py = I_NAI_H3xyz_S+ABY*I_NAI_G3xz_S;
  Double I_NAI_G2x2y_Py = I_NAI_H2x3y_S+ABY*I_NAI_G2x2y_S;
  Double I_NAI_G2xyz_Py = I_NAI_H2x2yz_S+ABY*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Py = I_NAI_H2xy2z_S+ABY*I_NAI_G2x2z_S;
  Double I_NAI_Gx3y_Py = I_NAI_Hx4y_S+ABY*I_NAI_Gx3y_S;
  Double I_NAI_Gx2yz_Py = I_NAI_Hx3yz_S+ABY*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Py = I_NAI_Hx2y2z_S+ABY*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Py = I_NAI_Hxy3z_S+ABY*I_NAI_Gx3z_S;
  Double I_NAI_G4y_Py = I_NAI_H5y_S+ABY*I_NAI_G4y_S;
  Double I_NAI_G3yz_Py = I_NAI_H4yz_S+ABY*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Py = I_NAI_H3y2z_S+ABY*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Py = I_NAI_H2y3z_S+ABY*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Py = I_NAI_Hy4z_S+ABY*I_NAI_G4z_S;
  Double I_NAI_G4x_Pz = I_NAI_H4xz_S+ABZ*I_NAI_G4x_S;
  Double I_NAI_G3xy_Pz = I_NAI_H3xyz_S+ABZ*I_NAI_G3xy_S;
  Double I_NAI_G3xz_Pz = I_NAI_H3x2z_S+ABZ*I_NAI_G3xz_S;
  Double I_NAI_G2x2y_Pz = I_NAI_H2x2yz_S+ABZ*I_NAI_G2x2y_S;
  Double I_NAI_G2xyz_Pz = I_NAI_H2xy2z_S+ABZ*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Pz = I_NAI_H2x3z_S+ABZ*I_NAI_G2x2z_S;
  Double I_NAI_Gx3y_Pz = I_NAI_Hx3yz_S+ABZ*I_NAI_Gx3y_S;
  Double I_NAI_Gx2yz_Pz = I_NAI_Hx2y2z_S+ABZ*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Pz = I_NAI_Hxy3z_S+ABZ*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Pz = I_NAI_Hx4z_S+ABZ*I_NAI_Gx3z_S;
  Double I_NAI_G4y_Pz = I_NAI_H4yz_S+ABZ*I_NAI_G4y_S;
  Double I_NAI_G3yz_Pz = I_NAI_H3y2z_S+ABZ*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Pz = I_NAI_H2y3z_S+ABZ*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Pz = I_NAI_Hy4z_S+ABZ*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Pz = I_NAI_H5z_S+ABZ*I_NAI_G4z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  Double I_NAI_F3x_D2x = I_NAI_G4x_Px+ABX*I_NAI_F3x_Px;
  Double I_NAI_F2xy_D2x = I_NAI_G3xy_Px+ABX*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_D2x = I_NAI_G3xz_Px+ABX*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_D2x = I_NAI_G2x2y_Px+ABX*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_D2x = I_NAI_G2xyz_Px+ABX*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_D2x = I_NAI_G2x2z_Px+ABX*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_D2x = I_NAI_Gx3y_Px+ABX*I_NAI_F3y_Px;
  Double I_NAI_F2yz_D2x = I_NAI_Gx2yz_Px+ABX*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_D2x = I_NAI_Gxy2z_Px+ABX*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_D2x = I_NAI_Gx3z_Px+ABX*I_NAI_F3z_Px;
  Double I_NAI_F3x_Dxy = I_NAI_G3xy_Px+ABY*I_NAI_F3x_Px;
  Double I_NAI_F2xy_Dxy = I_NAI_G2x2y_Px+ABY*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_Dxy = I_NAI_G2xyz_Px+ABY*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_Dxy = I_NAI_Gx3y_Px+ABY*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_Dxy = I_NAI_Gx2yz_Px+ABY*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_Dxy = I_NAI_Gxy2z_Px+ABY*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_Dxy = I_NAI_G4y_Px+ABY*I_NAI_F3y_Px;
  Double I_NAI_F2yz_Dxy = I_NAI_G3yz_Px+ABY*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_Dxy = I_NAI_G2y2z_Px+ABY*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_Dxy = I_NAI_Gy3z_Px+ABY*I_NAI_F3z_Px;
  Double I_NAI_F3x_Dxz = I_NAI_G3xz_Px+ABZ*I_NAI_F3x_Px;
  Double I_NAI_F2xy_Dxz = I_NAI_G2xyz_Px+ABZ*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_Dxz = I_NAI_G2x2z_Px+ABZ*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_Dxz = I_NAI_Gx2yz_Px+ABZ*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_Dxz = I_NAI_Gxy2z_Px+ABZ*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_Dxz = I_NAI_Gx3z_Px+ABZ*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_Dxz = I_NAI_G3yz_Px+ABZ*I_NAI_F3y_Px;
  Double I_NAI_F2yz_Dxz = I_NAI_G2y2z_Px+ABZ*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_Dxz = I_NAI_Gy3z_Px+ABZ*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_Dxz = I_NAI_G4z_Px+ABZ*I_NAI_F3z_Px;
  Double I_NAI_F3x_D2y = I_NAI_G3xy_Py+ABY*I_NAI_F3x_Py;
  Double I_NAI_F2xy_D2y = I_NAI_G2x2y_Py+ABY*I_NAI_F2xy_Py;
  Double I_NAI_F2xz_D2y = I_NAI_G2xyz_Py+ABY*I_NAI_F2xz_Py;
  Double I_NAI_Fx2y_D2y = I_NAI_Gx3y_Py+ABY*I_NAI_Fx2y_Py;
  Double I_NAI_Fxyz_D2y = I_NAI_Gx2yz_Py+ABY*I_NAI_Fxyz_Py;
  Double I_NAI_Fx2z_D2y = I_NAI_Gxy2z_Py+ABY*I_NAI_Fx2z_Py;
  Double I_NAI_F3y_D2y = I_NAI_G4y_Py+ABY*I_NAI_F3y_Py;
  Double I_NAI_F2yz_D2y = I_NAI_G3yz_Py+ABY*I_NAI_F2yz_Py;
  Double I_NAI_Fy2z_D2y = I_NAI_G2y2z_Py+ABY*I_NAI_Fy2z_Py;
  Double I_NAI_F3z_D2y = I_NAI_Gy3z_Py+ABY*I_NAI_F3z_Py;
  Double I_NAI_F3x_Dyz = I_NAI_G3xz_Py+ABZ*I_NAI_F3x_Py;
  Double I_NAI_F2xy_Dyz = I_NAI_G2xyz_Py+ABZ*I_NAI_F2xy_Py;
  Double I_NAI_F2xz_Dyz = I_NAI_G2x2z_Py+ABZ*I_NAI_F2xz_Py;
  Double I_NAI_Fx2y_Dyz = I_NAI_Gx2yz_Py+ABZ*I_NAI_Fx2y_Py;
  Double I_NAI_Fxyz_Dyz = I_NAI_Gxy2z_Py+ABZ*I_NAI_Fxyz_Py;
  Double I_NAI_Fx2z_Dyz = I_NAI_Gx3z_Py+ABZ*I_NAI_Fx2z_Py;
  Double I_NAI_F3y_Dyz = I_NAI_G3yz_Py+ABZ*I_NAI_F3y_Py;
  Double I_NAI_F2yz_Dyz = I_NAI_G2y2z_Py+ABZ*I_NAI_F2yz_Py;
  Double I_NAI_Fy2z_Dyz = I_NAI_Gy3z_Py+ABZ*I_NAI_Fy2z_Py;
  Double I_NAI_F3z_Dyz = I_NAI_G4z_Py+ABZ*I_NAI_F3z_Py;
  Double I_NAI_F3x_D2z = I_NAI_G3xz_Pz+ABZ*I_NAI_F3x_Pz;
  Double I_NAI_F2xy_D2z = I_NAI_G2xyz_Pz+ABZ*I_NAI_F2xy_Pz;
  Double I_NAI_F2xz_D2z = I_NAI_G2x2z_Pz+ABZ*I_NAI_F2xz_Pz;
  Double I_NAI_Fx2y_D2z = I_NAI_Gx2yz_Pz+ABZ*I_NAI_Fx2y_Pz;
  Double I_NAI_Fxyz_D2z = I_NAI_Gxy2z_Pz+ABZ*I_NAI_Fxyz_Pz;
  Double I_NAI_Fx2z_D2z = I_NAI_Gx3z_Pz+ABZ*I_NAI_Fx2z_Pz;
  Double I_NAI_F3y_D2z = I_NAI_G3yz_Pz+ABZ*I_NAI_F3y_Pz;
  Double I_NAI_F2yz_D2z = I_NAI_G2y2z_Pz+ABZ*I_NAI_F2yz_Pz;
  Double I_NAI_Fy2z_D2z = I_NAI_Gy3z_Pz+ABZ*I_NAI_Fy2z_Pz;
  Double I_NAI_F3z_D2z = I_NAI_G4z_Pz+ABZ*I_NAI_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_D_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D
   * RHS shell quartet name: SQ_NAI_D_D
   ************************************************************/
  Double I_NAI_D2x_F3x = I_NAI_F3x_D2x+ABX*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F3x = I_NAI_F2xy_D2x+ABX*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F3x = I_NAI_F2xz_D2x+ABX*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F3x = I_NAI_Fx2y_D2x+ABX*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F3x = I_NAI_Fxyz_D2x+ABX*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F3x = I_NAI_Fx2z_D2x+ABX*I_NAI_D2z_D2x;
  Double I_NAI_D2x_F2xy = I_NAI_F2xy_D2x+ABY*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F2xy = I_NAI_Fx2y_D2x+ABY*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F2xy = I_NAI_Fxyz_D2x+ABY*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F2xy = I_NAI_F3y_D2x+ABY*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F2xy = I_NAI_F2yz_D2x+ABY*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F2xy = I_NAI_Fy2z_D2x+ABY*I_NAI_D2z_D2x;
  Double I_NAI_D2x_F2xz = I_NAI_F2xz_D2x+ABZ*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F2xz = I_NAI_Fxyz_D2x+ABZ*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F2xz = I_NAI_Fx2z_D2x+ABZ*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F2xz = I_NAI_F2yz_D2x+ABZ*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F2xz = I_NAI_Fy2z_D2x+ABZ*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F2xz = I_NAI_F3z_D2x+ABZ*I_NAI_D2z_D2x;
  Double I_NAI_D2x_Fx2y = I_NAI_F3x_D2y+ABX*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_Fx2y = I_NAI_F2xy_D2y+ABX*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_Fx2y = I_NAI_F2xz_D2y+ABX*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_Fx2y = I_NAI_Fx2y_D2y+ABX*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_Fx2y = I_NAI_Fxyz_D2y+ABX*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_Fx2y = I_NAI_Fx2z_D2y+ABX*I_NAI_D2z_D2y;
  Double I_NAI_D2x_Fxyz = I_NAI_F2xz_Dxy+ABZ*I_NAI_D2x_Dxy;
  Double I_NAI_Dxy_Fxyz = I_NAI_Fxyz_Dxy+ABZ*I_NAI_Dxy_Dxy;
  Double I_NAI_Dxz_Fxyz = I_NAI_Fx2z_Dxy+ABZ*I_NAI_Dxz_Dxy;
  Double I_NAI_D2y_Fxyz = I_NAI_F2yz_Dxy+ABZ*I_NAI_D2y_Dxy;
  Double I_NAI_Dyz_Fxyz = I_NAI_Fy2z_Dxy+ABZ*I_NAI_Dyz_Dxy;
  Double I_NAI_D2z_Fxyz = I_NAI_F3z_Dxy+ABZ*I_NAI_D2z_Dxy;
  Double I_NAI_D2x_Fx2z = I_NAI_F3x_D2z+ABX*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_Fx2z = I_NAI_F2xy_D2z+ABX*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_Fx2z = I_NAI_F2xz_D2z+ABX*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_Fx2z = I_NAI_Fx2y_D2z+ABX*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_Fx2z = I_NAI_Fxyz_D2z+ABX*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_Fx2z = I_NAI_Fx2z_D2z+ABX*I_NAI_D2z_D2z;
  Double I_NAI_D2x_F3y = I_NAI_F2xy_D2y+ABY*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_F3y = I_NAI_Fx2y_D2y+ABY*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_F3y = I_NAI_Fxyz_D2y+ABY*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_F3y = I_NAI_F3y_D2y+ABY*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_F3y = I_NAI_F2yz_D2y+ABY*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_F3y = I_NAI_Fy2z_D2y+ABY*I_NAI_D2z_D2y;
  Double I_NAI_D2x_F2yz = I_NAI_F2xz_D2y+ABZ*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_F2yz = I_NAI_Fxyz_D2y+ABZ*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_F2yz = I_NAI_Fx2z_D2y+ABZ*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_F2yz = I_NAI_F2yz_D2y+ABZ*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_F2yz = I_NAI_Fy2z_D2y+ABZ*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_F2yz = I_NAI_F3z_D2y+ABZ*I_NAI_D2z_D2y;
  Double I_NAI_D2x_Fy2z = I_NAI_F2xy_D2z+ABY*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_Fy2z = I_NAI_Fx2y_D2z+ABY*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_Fy2z = I_NAI_Fxyz_D2z+ABY*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_Fy2z = I_NAI_F3y_D2z+ABY*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_Fy2z = I_NAI_F2yz_D2z+ABY*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_Fy2z = I_NAI_Fy2z_D2z+ABY*I_NAI_D2z_D2z;
  Double I_NAI_D2x_F3z = I_NAI_F2xz_D2z+ABZ*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_F3z = I_NAI_Fxyz_D2z+ABZ*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_F3z = I_NAI_Fx2z_D2z+ABZ*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_F3z = I_NAI_F2yz_D2z+ABZ*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_F3z = I_NAI_Fy2z_D2z+ABZ*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_F3z = I_NAI_F3z_D2z+ABZ*I_NAI_D2z_D2z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_a
   * RHS shell quartet name: SQ_NAI_G_S_a
   ************************************************************/
  Double I_NAI_G4x_Px_a = I_NAI_H5x_S_a+ABX*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Px_a = I_NAI_H4xy_S_a+ABX*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Px_a = I_NAI_H4xz_S_a+ABX*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Px_a = I_NAI_H3x2y_S_a+ABX*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Px_a = I_NAI_H3xyz_S_a+ABX*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Px_a = I_NAI_H3x2z_S_a+ABX*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Px_a = I_NAI_H2x3y_S_a+ABX*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Px_a = I_NAI_H2x2yz_S_a+ABX*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Px_a = I_NAI_H2xy2z_S_a+ABX*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Px_a = I_NAI_H2x3z_S_a+ABX*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Px_a = I_NAI_Hx4y_S_a+ABX*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Px_a = I_NAI_Hx3yz_S_a+ABX*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Px_a = I_NAI_Hx2y2z_S_a+ABX*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Px_a = I_NAI_Hxy3z_S_a+ABX*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Px_a = I_NAI_Hx4z_S_a+ABX*I_NAI_G4z_S_a;
  Double I_NAI_G4x_Py_a = I_NAI_H4xy_S_a+ABY*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Py_a = I_NAI_H3x2y_S_a+ABY*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Py_a = I_NAI_H3xyz_S_a+ABY*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Py_a = I_NAI_H2x3y_S_a+ABY*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Py_a = I_NAI_H2x2yz_S_a+ABY*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Py_a = I_NAI_H2xy2z_S_a+ABY*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Py_a = I_NAI_Hx4y_S_a+ABY*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Py_a = I_NAI_Hx3yz_S_a+ABY*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Py_a = I_NAI_Hx2y2z_S_a+ABY*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Py_a = I_NAI_Hxy3z_S_a+ABY*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Py_a = I_NAI_H5y_S_a+ABY*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Py_a = I_NAI_H4yz_S_a+ABY*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Py_a = I_NAI_H3y2z_S_a+ABY*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Py_a = I_NAI_H2y3z_S_a+ABY*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Py_a = I_NAI_Hy4z_S_a+ABY*I_NAI_G4z_S_a;
  Double I_NAI_G4x_Pz_a = I_NAI_H4xz_S_a+ABZ*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Pz_a = I_NAI_H3xyz_S_a+ABZ*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Pz_a = I_NAI_H3x2z_S_a+ABZ*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Pz_a = I_NAI_H2x2yz_S_a+ABZ*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Pz_a = I_NAI_H2xy2z_S_a+ABZ*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Pz_a = I_NAI_H2x3z_S_a+ABZ*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Pz_a = I_NAI_Hx3yz_S_a+ABZ*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Pz_a = I_NAI_Hx2y2z_S_a+ABZ*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Pz_a = I_NAI_Hxy3z_S_a+ABZ*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Pz_a = I_NAI_Hx4z_S_a+ABZ*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Pz_a = I_NAI_H4yz_S_a+ABZ*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Pz_a = I_NAI_H3y2z_S_a+ABZ*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Pz_a = I_NAI_H2y3z_S_a+ABZ*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Pz_a = I_NAI_Hy4z_S_a+ABZ*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Pz_a = I_NAI_H5z_S_a+ABZ*I_NAI_G4z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_H_S_a
   ************************************************************/
  Double I_NAI_H5x_Px_a = I_NAI_I6x_S_a+ABX*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Px_a = I_NAI_I5xy_S_a+ABX*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Px_a = I_NAI_I5xz_S_a+ABX*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Px_a = I_NAI_I4x2y_S_a+ABX*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Px_a = I_NAI_I4xyz_S_a+ABX*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Px_a = I_NAI_I4x2z_S_a+ABX*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Px_a = I_NAI_I3x3y_S_a+ABX*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Px_a = I_NAI_I3x2yz_S_a+ABX*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Px_a = I_NAI_I3xy2z_S_a+ABX*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Px_a = I_NAI_I3x3z_S_a+ABX*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Px_a = I_NAI_I2x4y_S_a+ABX*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Px_a = I_NAI_I2x3yz_S_a+ABX*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Px_a = I_NAI_I2x2y2z_S_a+ABX*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Px_a = I_NAI_I2xy3z_S_a+ABX*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Px_a = I_NAI_I2x4z_S_a+ABX*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Px_a = I_NAI_Ix5y_S_a+ABX*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Px_a = I_NAI_Ix4yz_S_a+ABX*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Px_a = I_NAI_Ix3y2z_S_a+ABX*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Px_a = I_NAI_Ix2y3z_S_a+ABX*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Px_a = I_NAI_Ixy4z_S_a+ABX*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Px_a = I_NAI_Ix5z_S_a+ABX*I_NAI_H5z_S_a;
  Double I_NAI_H5x_Py_a = I_NAI_I5xy_S_a+ABY*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Py_a = I_NAI_I4x2y_S_a+ABY*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Py_a = I_NAI_I4xyz_S_a+ABY*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Py_a = I_NAI_I3x3y_S_a+ABY*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Py_a = I_NAI_I3x2yz_S_a+ABY*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Py_a = I_NAI_I3xy2z_S_a+ABY*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Py_a = I_NAI_I2x4y_S_a+ABY*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Py_a = I_NAI_I2x3yz_S_a+ABY*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Py_a = I_NAI_I2x2y2z_S_a+ABY*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Py_a = I_NAI_I2xy3z_S_a+ABY*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Py_a = I_NAI_Ix5y_S_a+ABY*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Py_a = I_NAI_Ix4yz_S_a+ABY*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Py_a = I_NAI_Ix3y2z_S_a+ABY*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Py_a = I_NAI_Ix2y3z_S_a+ABY*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Py_a = I_NAI_Ixy4z_S_a+ABY*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Py_a = I_NAI_I6y_S_a+ABY*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Py_a = I_NAI_I5yz_S_a+ABY*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Py_a = I_NAI_I4y2z_S_a+ABY*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Py_a = I_NAI_I3y3z_S_a+ABY*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Py_a = I_NAI_I2y4z_S_a+ABY*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Py_a = I_NAI_Iy5z_S_a+ABY*I_NAI_H5z_S_a;
  Double I_NAI_H5x_Pz_a = I_NAI_I5xz_S_a+ABZ*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Pz_a = I_NAI_I4xyz_S_a+ABZ*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Pz_a = I_NAI_I4x2z_S_a+ABZ*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Pz_a = I_NAI_I3x2yz_S_a+ABZ*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Pz_a = I_NAI_I3xy2z_S_a+ABZ*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Pz_a = I_NAI_I3x3z_S_a+ABZ*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Pz_a = I_NAI_I2x3yz_S_a+ABZ*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Pz_a = I_NAI_I2x2y2z_S_a+ABZ*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Pz_a = I_NAI_I2xy3z_S_a+ABZ*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Pz_a = I_NAI_I2x4z_S_a+ABZ*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Pz_a = I_NAI_Ix4yz_S_a+ABZ*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Pz_a = I_NAI_Ix3y2z_S_a+ABZ*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Pz_a = I_NAI_Ix2y3z_S_a+ABZ*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Pz_a = I_NAI_Ixy4z_S_a+ABZ*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Pz_a = I_NAI_Ix5z_S_a+ABZ*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Pz_a = I_NAI_I5yz_S_a+ABZ*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Pz_a = I_NAI_I4y2z_S_a+ABZ*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Pz_a = I_NAI_I3y3z_S_a+ABZ*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Pz_a = I_NAI_I2y4z_S_a+ABZ*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Pz_a = I_NAI_Iy5z_S_a+ABZ*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Pz_a = I_NAI_I6z_S_a+ABZ*I_NAI_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_G_P_a
   ************************************************************/
  Double I_NAI_G4x_D2x_a = I_NAI_H5x_Px_a+ABX*I_NAI_G4x_Px_a;
  Double I_NAI_G3xy_D2x_a = I_NAI_H4xy_Px_a+ABX*I_NAI_G3xy_Px_a;
  Double I_NAI_G3xz_D2x_a = I_NAI_H4xz_Px_a+ABX*I_NAI_G3xz_Px_a;
  Double I_NAI_G2x2y_D2x_a = I_NAI_H3x2y_Px_a+ABX*I_NAI_G2x2y_Px_a;
  Double I_NAI_G2xyz_D2x_a = I_NAI_H3xyz_Px_a+ABX*I_NAI_G2xyz_Px_a;
  Double I_NAI_G2x2z_D2x_a = I_NAI_H3x2z_Px_a+ABX*I_NAI_G2x2z_Px_a;
  Double I_NAI_Gx3y_D2x_a = I_NAI_H2x3y_Px_a+ABX*I_NAI_Gx3y_Px_a;
  Double I_NAI_Gx2yz_D2x_a = I_NAI_H2x2yz_Px_a+ABX*I_NAI_Gx2yz_Px_a;
  Double I_NAI_Gxy2z_D2x_a = I_NAI_H2xy2z_Px_a+ABX*I_NAI_Gxy2z_Px_a;
  Double I_NAI_Gx3z_D2x_a = I_NAI_H2x3z_Px_a+ABX*I_NAI_Gx3z_Px_a;
  Double I_NAI_G4y_D2x_a = I_NAI_Hx4y_Px_a+ABX*I_NAI_G4y_Px_a;
  Double I_NAI_G3yz_D2x_a = I_NAI_Hx3yz_Px_a+ABX*I_NAI_G3yz_Px_a;
  Double I_NAI_G2y2z_D2x_a = I_NAI_Hx2y2z_Px_a+ABX*I_NAI_G2y2z_Px_a;
  Double I_NAI_Gy3z_D2x_a = I_NAI_Hxy3z_Px_a+ABX*I_NAI_Gy3z_Px_a;
  Double I_NAI_G4z_D2x_a = I_NAI_Hx4z_Px_a+ABX*I_NAI_G4z_Px_a;
  Double I_NAI_G4x_Dxy_a = I_NAI_H4xy_Px_a+ABY*I_NAI_G4x_Px_a;
  Double I_NAI_G3xy_Dxy_a = I_NAI_H3x2y_Px_a+ABY*I_NAI_G3xy_Px_a;
  Double I_NAI_G3xz_Dxy_a = I_NAI_H3xyz_Px_a+ABY*I_NAI_G3xz_Px_a;
  Double I_NAI_G2x2y_Dxy_a = I_NAI_H2x3y_Px_a+ABY*I_NAI_G2x2y_Px_a;
  Double I_NAI_G2xyz_Dxy_a = I_NAI_H2x2yz_Px_a+ABY*I_NAI_G2xyz_Px_a;
  Double I_NAI_G2x2z_Dxy_a = I_NAI_H2xy2z_Px_a+ABY*I_NAI_G2x2z_Px_a;
  Double I_NAI_Gx3y_Dxy_a = I_NAI_Hx4y_Px_a+ABY*I_NAI_Gx3y_Px_a;
  Double I_NAI_Gx2yz_Dxy_a = I_NAI_Hx3yz_Px_a+ABY*I_NAI_Gx2yz_Px_a;
  Double I_NAI_Gxy2z_Dxy_a = I_NAI_Hx2y2z_Px_a+ABY*I_NAI_Gxy2z_Px_a;
  Double I_NAI_Gx3z_Dxy_a = I_NAI_Hxy3z_Px_a+ABY*I_NAI_Gx3z_Px_a;
  Double I_NAI_G4y_Dxy_a = I_NAI_H5y_Px_a+ABY*I_NAI_G4y_Px_a;
  Double I_NAI_G3yz_Dxy_a = I_NAI_H4yz_Px_a+ABY*I_NAI_G3yz_Px_a;
  Double I_NAI_G2y2z_Dxy_a = I_NAI_H3y2z_Px_a+ABY*I_NAI_G2y2z_Px_a;
  Double I_NAI_Gy3z_Dxy_a = I_NAI_H2y3z_Px_a+ABY*I_NAI_Gy3z_Px_a;
  Double I_NAI_G4z_Dxy_a = I_NAI_Hy4z_Px_a+ABY*I_NAI_G4z_Px_a;
  Double I_NAI_G4x_D2y_a = I_NAI_H4xy_Py_a+ABY*I_NAI_G4x_Py_a;
  Double I_NAI_G3xy_D2y_a = I_NAI_H3x2y_Py_a+ABY*I_NAI_G3xy_Py_a;
  Double I_NAI_G3xz_D2y_a = I_NAI_H3xyz_Py_a+ABY*I_NAI_G3xz_Py_a;
  Double I_NAI_G2x2y_D2y_a = I_NAI_H2x3y_Py_a+ABY*I_NAI_G2x2y_Py_a;
  Double I_NAI_G2xyz_D2y_a = I_NAI_H2x2yz_Py_a+ABY*I_NAI_G2xyz_Py_a;
  Double I_NAI_G2x2z_D2y_a = I_NAI_H2xy2z_Py_a+ABY*I_NAI_G2x2z_Py_a;
  Double I_NAI_Gx3y_D2y_a = I_NAI_Hx4y_Py_a+ABY*I_NAI_Gx3y_Py_a;
  Double I_NAI_Gx2yz_D2y_a = I_NAI_Hx3yz_Py_a+ABY*I_NAI_Gx2yz_Py_a;
  Double I_NAI_Gxy2z_D2y_a = I_NAI_Hx2y2z_Py_a+ABY*I_NAI_Gxy2z_Py_a;
  Double I_NAI_Gx3z_D2y_a = I_NAI_Hxy3z_Py_a+ABY*I_NAI_Gx3z_Py_a;
  Double I_NAI_G4y_D2y_a = I_NAI_H5y_Py_a+ABY*I_NAI_G4y_Py_a;
  Double I_NAI_G3yz_D2y_a = I_NAI_H4yz_Py_a+ABY*I_NAI_G3yz_Py_a;
  Double I_NAI_G2y2z_D2y_a = I_NAI_H3y2z_Py_a+ABY*I_NAI_G2y2z_Py_a;
  Double I_NAI_Gy3z_D2y_a = I_NAI_H2y3z_Py_a+ABY*I_NAI_Gy3z_Py_a;
  Double I_NAI_G4z_D2y_a = I_NAI_Hy4z_Py_a+ABY*I_NAI_G4z_Py_a;
  Double I_NAI_G4x_D2z_a = I_NAI_H4xz_Pz_a+ABZ*I_NAI_G4x_Pz_a;
  Double I_NAI_G3xy_D2z_a = I_NAI_H3xyz_Pz_a+ABZ*I_NAI_G3xy_Pz_a;
  Double I_NAI_G3xz_D2z_a = I_NAI_H3x2z_Pz_a+ABZ*I_NAI_G3xz_Pz_a;
  Double I_NAI_G2x2y_D2z_a = I_NAI_H2x2yz_Pz_a+ABZ*I_NAI_G2x2y_Pz_a;
  Double I_NAI_G2xyz_D2z_a = I_NAI_H2xy2z_Pz_a+ABZ*I_NAI_G2xyz_Pz_a;
  Double I_NAI_G2x2z_D2z_a = I_NAI_H2x3z_Pz_a+ABZ*I_NAI_G2x2z_Pz_a;
  Double I_NAI_Gx3y_D2z_a = I_NAI_Hx3yz_Pz_a+ABZ*I_NAI_Gx3y_Pz_a;
  Double I_NAI_Gx2yz_D2z_a = I_NAI_Hx2y2z_Pz_a+ABZ*I_NAI_Gx2yz_Pz_a;
  Double I_NAI_Gxy2z_D2z_a = I_NAI_Hxy3z_Pz_a+ABZ*I_NAI_Gxy2z_Pz_a;
  Double I_NAI_Gx3z_D2z_a = I_NAI_Hx4z_Pz_a+ABZ*I_NAI_Gx3z_Pz_a;
  Double I_NAI_G4y_D2z_a = I_NAI_H4yz_Pz_a+ABZ*I_NAI_G4y_Pz_a;
  Double I_NAI_G3yz_D2z_a = I_NAI_H3y2z_Pz_a+ABZ*I_NAI_G3yz_Pz_a;
  Double I_NAI_G2y2z_D2z_a = I_NAI_H2y3z_Pz_a+ABZ*I_NAI_G2y2z_Pz_a;
  Double I_NAI_Gy3z_D2z_a = I_NAI_Hy4z_Pz_a+ABZ*I_NAI_Gy3z_Pz_a;
  Double I_NAI_G4z_D2z_a = I_NAI_H5z_Pz_a+ABZ*I_NAI_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 8 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_a
   * RHS shell quartet name: SQ_NAI_I_S_a
   ************************************************************/
  Double I_NAI_I6x_Px_a = I_NAI_K7x_S_a+ABX*I_NAI_I6x_S_a;
  Double I_NAI_I5xy_Px_a = I_NAI_K6xy_S_a+ABX*I_NAI_I5xy_S_a;
  Double I_NAI_I5xz_Px_a = I_NAI_K6xz_S_a+ABX*I_NAI_I5xz_S_a;
  Double I_NAI_I4x2y_Px_a = I_NAI_K5x2y_S_a+ABX*I_NAI_I4x2y_S_a;
  Double I_NAI_I4xyz_Px_a = I_NAI_K5xyz_S_a+ABX*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Px_a = I_NAI_K5x2z_S_a+ABX*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x3y_Px_a = I_NAI_K4x3y_S_a+ABX*I_NAI_I3x3y_S_a;
  Double I_NAI_I3x2yz_Px_a = I_NAI_K4x2yz_S_a+ABX*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Px_a = I_NAI_K4xy2z_S_a+ABX*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Px_a = I_NAI_K4x3z_S_a+ABX*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x4y_Px_a = I_NAI_K3x4y_S_a+ABX*I_NAI_I2x4y_S_a;
  Double I_NAI_I2x3yz_Px_a = I_NAI_K3x3yz_S_a+ABX*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Px_a = I_NAI_K3x2y2z_S_a+ABX*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Px_a = I_NAI_K3xy3z_S_a+ABX*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Px_a = I_NAI_K3x4z_S_a+ABX*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix5y_Px_a = I_NAI_K2x5y_S_a+ABX*I_NAI_Ix5y_S_a;
  Double I_NAI_Ix4yz_Px_a = I_NAI_K2x4yz_S_a+ABX*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Px_a = I_NAI_K2x3y2z_S_a+ABX*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Px_a = I_NAI_K2x2y3z_S_a+ABX*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Px_a = I_NAI_K2xy4z_S_a+ABX*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Px_a = I_NAI_K2x5z_S_a+ABX*I_NAI_Ix5z_S_a;
  Double I_NAI_I6y_Px_a = I_NAI_Kx6y_S_a+ABX*I_NAI_I6y_S_a;
  Double I_NAI_I5yz_Px_a = I_NAI_Kx5yz_S_a+ABX*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Px_a = I_NAI_Kx4y2z_S_a+ABX*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Px_a = I_NAI_Kx3y3z_S_a+ABX*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Px_a = I_NAI_Kx2y4z_S_a+ABX*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Px_a = I_NAI_Kxy5z_S_a+ABX*I_NAI_Iy5z_S_a;
  Double I_NAI_I6z_Px_a = I_NAI_Kx6z_S_a+ABX*I_NAI_I6z_S_a;
  Double I_NAI_I5xy_Py_a = I_NAI_K5x2y_S_a+ABY*I_NAI_I5xy_S_a;
  Double I_NAI_I5xz_Py_a = I_NAI_K5xyz_S_a+ABY*I_NAI_I5xz_S_a;
  Double I_NAI_I4x2y_Py_a = I_NAI_K4x3y_S_a+ABY*I_NAI_I4x2y_S_a;
  Double I_NAI_I4xyz_Py_a = I_NAI_K4x2yz_S_a+ABY*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Py_a = I_NAI_K4xy2z_S_a+ABY*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x3y_Py_a = I_NAI_K3x4y_S_a+ABY*I_NAI_I3x3y_S_a;
  Double I_NAI_I3x2yz_Py_a = I_NAI_K3x3yz_S_a+ABY*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Py_a = I_NAI_K3x2y2z_S_a+ABY*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Py_a = I_NAI_K3xy3z_S_a+ABY*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x4y_Py_a = I_NAI_K2x5y_S_a+ABY*I_NAI_I2x4y_S_a;
  Double I_NAI_I2x3yz_Py_a = I_NAI_K2x4yz_S_a+ABY*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Py_a = I_NAI_K2x3y2z_S_a+ABY*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Py_a = I_NAI_K2x2y3z_S_a+ABY*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Py_a = I_NAI_K2xy4z_S_a+ABY*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix5y_Py_a = I_NAI_Kx6y_S_a+ABY*I_NAI_Ix5y_S_a;
  Double I_NAI_Ix4yz_Py_a = I_NAI_Kx5yz_S_a+ABY*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Py_a = I_NAI_Kx4y2z_S_a+ABY*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Py_a = I_NAI_Kx3y3z_S_a+ABY*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Py_a = I_NAI_Kx2y4z_S_a+ABY*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Py_a = I_NAI_Kxy5z_S_a+ABY*I_NAI_Ix5z_S_a;
  Double I_NAI_I6y_Py_a = I_NAI_K7y_S_a+ABY*I_NAI_I6y_S_a;
  Double I_NAI_I5yz_Py_a = I_NAI_K6yz_S_a+ABY*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Py_a = I_NAI_K5y2z_S_a+ABY*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Py_a = I_NAI_K4y3z_S_a+ABY*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Py_a = I_NAI_K3y4z_S_a+ABY*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Py_a = I_NAI_K2y5z_S_a+ABY*I_NAI_Iy5z_S_a;
  Double I_NAI_I6z_Py_a = I_NAI_Ky6z_S_a+ABY*I_NAI_I6z_S_a;
  Double I_NAI_I5xz_Pz_a = I_NAI_K5x2z_S_a+ABZ*I_NAI_I5xz_S_a;
  Double I_NAI_I4xyz_Pz_a = I_NAI_K4xy2z_S_a+ABZ*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Pz_a = I_NAI_K4x3z_S_a+ABZ*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x2yz_Pz_a = I_NAI_K3x2y2z_S_a+ABZ*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Pz_a = I_NAI_K3xy3z_S_a+ABZ*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Pz_a = I_NAI_K3x4z_S_a+ABZ*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x3yz_Pz_a = I_NAI_K2x3y2z_S_a+ABZ*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Pz_a = I_NAI_K2x2y3z_S_a+ABZ*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Pz_a = I_NAI_K2xy4z_S_a+ABZ*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Pz_a = I_NAI_K2x5z_S_a+ABZ*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix4yz_Pz_a = I_NAI_Kx4y2z_S_a+ABZ*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Pz_a = I_NAI_Kx3y3z_S_a+ABZ*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Pz_a = I_NAI_Kx2y4z_S_a+ABZ*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Pz_a = I_NAI_Kxy5z_S_a+ABZ*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Pz_a = I_NAI_Kx6z_S_a+ABZ*I_NAI_Ix5z_S_a;
  Double I_NAI_I5yz_Pz_a = I_NAI_K5y2z_S_a+ABZ*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Pz_a = I_NAI_K4y3z_S_a+ABZ*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Pz_a = I_NAI_K3y4z_S_a+ABZ*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Pz_a = I_NAI_K2y5z_S_a+ABZ*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Pz_a = I_NAI_Ky6z_S_a+ABZ*I_NAI_Iy5z_S_a;
  Double I_NAI_I6z_Pz_a = I_NAI_K7z_S_a+ABZ*I_NAI_I6z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   ************************************************************/
  Double I_NAI_H5x_D2x_a = I_NAI_I6x_Px_a+ABX*I_NAI_H5x_Px_a;
  Double I_NAI_H4xy_D2x_a = I_NAI_I5xy_Px_a+ABX*I_NAI_H4xy_Px_a;
  Double I_NAI_H4xz_D2x_a = I_NAI_I5xz_Px_a+ABX*I_NAI_H4xz_Px_a;
  Double I_NAI_H3x2y_D2x_a = I_NAI_I4x2y_Px_a+ABX*I_NAI_H3x2y_Px_a;
  Double I_NAI_H3xyz_D2x_a = I_NAI_I4xyz_Px_a+ABX*I_NAI_H3xyz_Px_a;
  Double I_NAI_H3x2z_D2x_a = I_NAI_I4x2z_Px_a+ABX*I_NAI_H3x2z_Px_a;
  Double I_NAI_H2x3y_D2x_a = I_NAI_I3x3y_Px_a+ABX*I_NAI_H2x3y_Px_a;
  Double I_NAI_H2x2yz_D2x_a = I_NAI_I3x2yz_Px_a+ABX*I_NAI_H2x2yz_Px_a;
  Double I_NAI_H2xy2z_D2x_a = I_NAI_I3xy2z_Px_a+ABX*I_NAI_H2xy2z_Px_a;
  Double I_NAI_H2x3z_D2x_a = I_NAI_I3x3z_Px_a+ABX*I_NAI_H2x3z_Px_a;
  Double I_NAI_Hx4y_D2x_a = I_NAI_I2x4y_Px_a+ABX*I_NAI_Hx4y_Px_a;
  Double I_NAI_Hx3yz_D2x_a = I_NAI_I2x3yz_Px_a+ABX*I_NAI_Hx3yz_Px_a;
  Double I_NAI_Hx2y2z_D2x_a = I_NAI_I2x2y2z_Px_a+ABX*I_NAI_Hx2y2z_Px_a;
  Double I_NAI_Hxy3z_D2x_a = I_NAI_I2xy3z_Px_a+ABX*I_NAI_Hxy3z_Px_a;
  Double I_NAI_Hx4z_D2x_a = I_NAI_I2x4z_Px_a+ABX*I_NAI_Hx4z_Px_a;
  Double I_NAI_H5y_D2x_a = I_NAI_Ix5y_Px_a+ABX*I_NAI_H5y_Px_a;
  Double I_NAI_H4yz_D2x_a = I_NAI_Ix4yz_Px_a+ABX*I_NAI_H4yz_Px_a;
  Double I_NAI_H3y2z_D2x_a = I_NAI_Ix3y2z_Px_a+ABX*I_NAI_H3y2z_Px_a;
  Double I_NAI_H2y3z_D2x_a = I_NAI_Ix2y3z_Px_a+ABX*I_NAI_H2y3z_Px_a;
  Double I_NAI_Hy4z_D2x_a = I_NAI_Ixy4z_Px_a+ABX*I_NAI_Hy4z_Px_a;
  Double I_NAI_H5z_D2x_a = I_NAI_Ix5z_Px_a+ABX*I_NAI_H5z_Px_a;
  Double I_NAI_H5x_Dxy_a = I_NAI_I5xy_Px_a+ABY*I_NAI_H5x_Px_a;
  Double I_NAI_H4xy_Dxy_a = I_NAI_I4x2y_Px_a+ABY*I_NAI_H4xy_Px_a;
  Double I_NAI_H4xz_Dxy_a = I_NAI_I4xyz_Px_a+ABY*I_NAI_H4xz_Px_a;
  Double I_NAI_H3x2y_Dxy_a = I_NAI_I3x3y_Px_a+ABY*I_NAI_H3x2y_Px_a;
  Double I_NAI_H3xyz_Dxy_a = I_NAI_I3x2yz_Px_a+ABY*I_NAI_H3xyz_Px_a;
  Double I_NAI_H3x2z_Dxy_a = I_NAI_I3xy2z_Px_a+ABY*I_NAI_H3x2z_Px_a;
  Double I_NAI_H2x3y_Dxy_a = I_NAI_I2x4y_Px_a+ABY*I_NAI_H2x3y_Px_a;
  Double I_NAI_H2x2yz_Dxy_a = I_NAI_I2x3yz_Px_a+ABY*I_NAI_H2x2yz_Px_a;
  Double I_NAI_H2xy2z_Dxy_a = I_NAI_I2x2y2z_Px_a+ABY*I_NAI_H2xy2z_Px_a;
  Double I_NAI_H2x3z_Dxy_a = I_NAI_I2xy3z_Px_a+ABY*I_NAI_H2x3z_Px_a;
  Double I_NAI_Hx4y_Dxy_a = I_NAI_Ix5y_Px_a+ABY*I_NAI_Hx4y_Px_a;
  Double I_NAI_Hx3yz_Dxy_a = I_NAI_Ix4yz_Px_a+ABY*I_NAI_Hx3yz_Px_a;
  Double I_NAI_Hx2y2z_Dxy_a = I_NAI_Ix3y2z_Px_a+ABY*I_NAI_Hx2y2z_Px_a;
  Double I_NAI_Hxy3z_Dxy_a = I_NAI_Ix2y3z_Px_a+ABY*I_NAI_Hxy3z_Px_a;
  Double I_NAI_Hx4z_Dxy_a = I_NAI_Ixy4z_Px_a+ABY*I_NAI_Hx4z_Px_a;
  Double I_NAI_H5y_Dxy_a = I_NAI_I6y_Px_a+ABY*I_NAI_H5y_Px_a;
  Double I_NAI_H4yz_Dxy_a = I_NAI_I5yz_Px_a+ABY*I_NAI_H4yz_Px_a;
  Double I_NAI_H3y2z_Dxy_a = I_NAI_I4y2z_Px_a+ABY*I_NAI_H3y2z_Px_a;
  Double I_NAI_H2y3z_Dxy_a = I_NAI_I3y3z_Px_a+ABY*I_NAI_H2y3z_Px_a;
  Double I_NAI_Hy4z_Dxy_a = I_NAI_I2y4z_Px_a+ABY*I_NAI_Hy4z_Px_a;
  Double I_NAI_H5z_Dxy_a = I_NAI_Iy5z_Px_a+ABY*I_NAI_H5z_Px_a;
  Double I_NAI_H5x_Dxz_a = I_NAI_I5xz_Px_a+ABZ*I_NAI_H5x_Px_a;
  Double I_NAI_H4xy_Dxz_a = I_NAI_I4xyz_Px_a+ABZ*I_NAI_H4xy_Px_a;
  Double I_NAI_H4xz_Dxz_a = I_NAI_I4x2z_Px_a+ABZ*I_NAI_H4xz_Px_a;
  Double I_NAI_H3x2y_Dxz_a = I_NAI_I3x2yz_Px_a+ABZ*I_NAI_H3x2y_Px_a;
  Double I_NAI_H3xyz_Dxz_a = I_NAI_I3xy2z_Px_a+ABZ*I_NAI_H3xyz_Px_a;
  Double I_NAI_H3x2z_Dxz_a = I_NAI_I3x3z_Px_a+ABZ*I_NAI_H3x2z_Px_a;
  Double I_NAI_H2x3y_Dxz_a = I_NAI_I2x3yz_Px_a+ABZ*I_NAI_H2x3y_Px_a;
  Double I_NAI_H2x2yz_Dxz_a = I_NAI_I2x2y2z_Px_a+ABZ*I_NAI_H2x2yz_Px_a;
  Double I_NAI_H2xy2z_Dxz_a = I_NAI_I2xy3z_Px_a+ABZ*I_NAI_H2xy2z_Px_a;
  Double I_NAI_H2x3z_Dxz_a = I_NAI_I2x4z_Px_a+ABZ*I_NAI_H2x3z_Px_a;
  Double I_NAI_Hx4y_Dxz_a = I_NAI_Ix4yz_Px_a+ABZ*I_NAI_Hx4y_Px_a;
  Double I_NAI_Hx3yz_Dxz_a = I_NAI_Ix3y2z_Px_a+ABZ*I_NAI_Hx3yz_Px_a;
  Double I_NAI_Hx2y2z_Dxz_a = I_NAI_Ix2y3z_Px_a+ABZ*I_NAI_Hx2y2z_Px_a;
  Double I_NAI_Hxy3z_Dxz_a = I_NAI_Ixy4z_Px_a+ABZ*I_NAI_Hxy3z_Px_a;
  Double I_NAI_Hx4z_Dxz_a = I_NAI_Ix5z_Px_a+ABZ*I_NAI_Hx4z_Px_a;
  Double I_NAI_H5y_Dxz_a = I_NAI_I5yz_Px_a+ABZ*I_NAI_H5y_Px_a;
  Double I_NAI_H4yz_Dxz_a = I_NAI_I4y2z_Px_a+ABZ*I_NAI_H4yz_Px_a;
  Double I_NAI_H3y2z_Dxz_a = I_NAI_I3y3z_Px_a+ABZ*I_NAI_H3y2z_Px_a;
  Double I_NAI_H2y3z_Dxz_a = I_NAI_I2y4z_Px_a+ABZ*I_NAI_H2y3z_Px_a;
  Double I_NAI_Hy4z_Dxz_a = I_NAI_Iy5z_Px_a+ABZ*I_NAI_Hy4z_Px_a;
  Double I_NAI_H5z_Dxz_a = I_NAI_I6z_Px_a+ABZ*I_NAI_H5z_Px_a;
  Double I_NAI_H5x_D2y_a = I_NAI_I5xy_Py_a+ABY*I_NAI_H5x_Py_a;
  Double I_NAI_H4xy_D2y_a = I_NAI_I4x2y_Py_a+ABY*I_NAI_H4xy_Py_a;
  Double I_NAI_H4xz_D2y_a = I_NAI_I4xyz_Py_a+ABY*I_NAI_H4xz_Py_a;
  Double I_NAI_H3x2y_D2y_a = I_NAI_I3x3y_Py_a+ABY*I_NAI_H3x2y_Py_a;
  Double I_NAI_H3xyz_D2y_a = I_NAI_I3x2yz_Py_a+ABY*I_NAI_H3xyz_Py_a;
  Double I_NAI_H3x2z_D2y_a = I_NAI_I3xy2z_Py_a+ABY*I_NAI_H3x2z_Py_a;
  Double I_NAI_H2x3y_D2y_a = I_NAI_I2x4y_Py_a+ABY*I_NAI_H2x3y_Py_a;
  Double I_NAI_H2x2yz_D2y_a = I_NAI_I2x3yz_Py_a+ABY*I_NAI_H2x2yz_Py_a;
  Double I_NAI_H2xy2z_D2y_a = I_NAI_I2x2y2z_Py_a+ABY*I_NAI_H2xy2z_Py_a;
  Double I_NAI_H2x3z_D2y_a = I_NAI_I2xy3z_Py_a+ABY*I_NAI_H2x3z_Py_a;
  Double I_NAI_Hx4y_D2y_a = I_NAI_Ix5y_Py_a+ABY*I_NAI_Hx4y_Py_a;
  Double I_NAI_Hx3yz_D2y_a = I_NAI_Ix4yz_Py_a+ABY*I_NAI_Hx3yz_Py_a;
  Double I_NAI_Hx2y2z_D2y_a = I_NAI_Ix3y2z_Py_a+ABY*I_NAI_Hx2y2z_Py_a;
  Double I_NAI_Hxy3z_D2y_a = I_NAI_Ix2y3z_Py_a+ABY*I_NAI_Hxy3z_Py_a;
  Double I_NAI_Hx4z_D2y_a = I_NAI_Ixy4z_Py_a+ABY*I_NAI_Hx4z_Py_a;
  Double I_NAI_H5y_D2y_a = I_NAI_I6y_Py_a+ABY*I_NAI_H5y_Py_a;
  Double I_NAI_H4yz_D2y_a = I_NAI_I5yz_Py_a+ABY*I_NAI_H4yz_Py_a;
  Double I_NAI_H3y2z_D2y_a = I_NAI_I4y2z_Py_a+ABY*I_NAI_H3y2z_Py_a;
  Double I_NAI_H2y3z_D2y_a = I_NAI_I3y3z_Py_a+ABY*I_NAI_H2y3z_Py_a;
  Double I_NAI_Hy4z_D2y_a = I_NAI_I2y4z_Py_a+ABY*I_NAI_Hy4z_Py_a;
  Double I_NAI_H5z_D2y_a = I_NAI_Iy5z_Py_a+ABY*I_NAI_H5z_Py_a;
  Double I_NAI_H5x_Dyz_a = I_NAI_I5xz_Py_a+ABZ*I_NAI_H5x_Py_a;
  Double I_NAI_H4xy_Dyz_a = I_NAI_I4xyz_Py_a+ABZ*I_NAI_H4xy_Py_a;
  Double I_NAI_H4xz_Dyz_a = I_NAI_I4x2z_Py_a+ABZ*I_NAI_H4xz_Py_a;
  Double I_NAI_H3x2y_Dyz_a = I_NAI_I3x2yz_Py_a+ABZ*I_NAI_H3x2y_Py_a;
  Double I_NAI_H3xyz_Dyz_a = I_NAI_I3xy2z_Py_a+ABZ*I_NAI_H3xyz_Py_a;
  Double I_NAI_H3x2z_Dyz_a = I_NAI_I3x3z_Py_a+ABZ*I_NAI_H3x2z_Py_a;
  Double I_NAI_H2x3y_Dyz_a = I_NAI_I2x3yz_Py_a+ABZ*I_NAI_H2x3y_Py_a;
  Double I_NAI_H2x2yz_Dyz_a = I_NAI_I2x2y2z_Py_a+ABZ*I_NAI_H2x2yz_Py_a;
  Double I_NAI_H2xy2z_Dyz_a = I_NAI_I2xy3z_Py_a+ABZ*I_NAI_H2xy2z_Py_a;
  Double I_NAI_H2x3z_Dyz_a = I_NAI_I2x4z_Py_a+ABZ*I_NAI_H2x3z_Py_a;
  Double I_NAI_Hx4y_Dyz_a = I_NAI_Ix4yz_Py_a+ABZ*I_NAI_Hx4y_Py_a;
  Double I_NAI_Hx3yz_Dyz_a = I_NAI_Ix3y2z_Py_a+ABZ*I_NAI_Hx3yz_Py_a;
  Double I_NAI_Hx2y2z_Dyz_a = I_NAI_Ix2y3z_Py_a+ABZ*I_NAI_Hx2y2z_Py_a;
  Double I_NAI_Hxy3z_Dyz_a = I_NAI_Ixy4z_Py_a+ABZ*I_NAI_Hxy3z_Py_a;
  Double I_NAI_Hx4z_Dyz_a = I_NAI_Ix5z_Py_a+ABZ*I_NAI_Hx4z_Py_a;
  Double I_NAI_H5y_Dyz_a = I_NAI_I5yz_Py_a+ABZ*I_NAI_H5y_Py_a;
  Double I_NAI_H4yz_Dyz_a = I_NAI_I4y2z_Py_a+ABZ*I_NAI_H4yz_Py_a;
  Double I_NAI_H3y2z_Dyz_a = I_NAI_I3y3z_Py_a+ABZ*I_NAI_H3y2z_Py_a;
  Double I_NAI_H2y3z_Dyz_a = I_NAI_I2y4z_Py_a+ABZ*I_NAI_H2y3z_Py_a;
  Double I_NAI_Hy4z_Dyz_a = I_NAI_Iy5z_Py_a+ABZ*I_NAI_Hy4z_Py_a;
  Double I_NAI_H5z_Dyz_a = I_NAI_I6z_Py_a+ABZ*I_NAI_H5z_Py_a;
  Double I_NAI_H5x_D2z_a = I_NAI_I5xz_Pz_a+ABZ*I_NAI_H5x_Pz_a;
  Double I_NAI_H4xy_D2z_a = I_NAI_I4xyz_Pz_a+ABZ*I_NAI_H4xy_Pz_a;
  Double I_NAI_H4xz_D2z_a = I_NAI_I4x2z_Pz_a+ABZ*I_NAI_H4xz_Pz_a;
  Double I_NAI_H3x2y_D2z_a = I_NAI_I3x2yz_Pz_a+ABZ*I_NAI_H3x2y_Pz_a;
  Double I_NAI_H3xyz_D2z_a = I_NAI_I3xy2z_Pz_a+ABZ*I_NAI_H3xyz_Pz_a;
  Double I_NAI_H3x2z_D2z_a = I_NAI_I3x3z_Pz_a+ABZ*I_NAI_H3x2z_Pz_a;
  Double I_NAI_H2x3y_D2z_a = I_NAI_I2x3yz_Pz_a+ABZ*I_NAI_H2x3y_Pz_a;
  Double I_NAI_H2x2yz_D2z_a = I_NAI_I2x2y2z_Pz_a+ABZ*I_NAI_H2x2yz_Pz_a;
  Double I_NAI_H2xy2z_D2z_a = I_NAI_I2xy3z_Pz_a+ABZ*I_NAI_H2xy2z_Pz_a;
  Double I_NAI_H2x3z_D2z_a = I_NAI_I2x4z_Pz_a+ABZ*I_NAI_H2x3z_Pz_a;
  Double I_NAI_Hx4y_D2z_a = I_NAI_Ix4yz_Pz_a+ABZ*I_NAI_Hx4y_Pz_a;
  Double I_NAI_Hx3yz_D2z_a = I_NAI_Ix3y2z_Pz_a+ABZ*I_NAI_Hx3yz_Pz_a;
  Double I_NAI_Hx2y2z_D2z_a = I_NAI_Ix2y3z_Pz_a+ABZ*I_NAI_Hx2y2z_Pz_a;
  Double I_NAI_Hxy3z_D2z_a = I_NAI_Ixy4z_Pz_a+ABZ*I_NAI_Hxy3z_Pz_a;
  Double I_NAI_Hx4z_D2z_a = I_NAI_Ix5z_Pz_a+ABZ*I_NAI_Hx4z_Pz_a;
  Double I_NAI_H5y_D2z_a = I_NAI_I5yz_Pz_a+ABZ*I_NAI_H5y_Pz_a;
  Double I_NAI_H4yz_D2z_a = I_NAI_I4y2z_Pz_a+ABZ*I_NAI_H4yz_Pz_a;
  Double I_NAI_H3y2z_D2z_a = I_NAI_I3y3z_Pz_a+ABZ*I_NAI_H3y2z_Pz_a;
  Double I_NAI_H2y3z_D2z_a = I_NAI_I2y4z_Pz_a+ABZ*I_NAI_H2y3z_Pz_a;
  Double I_NAI_Hy4z_D2z_a = I_NAI_Iy5z_Pz_a+ABZ*I_NAI_Hy4z_Pz_a;
  Double I_NAI_H5z_D2z_a = I_NAI_I6z_Pz_a+ABZ*I_NAI_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_G_D_a
   ************************************************************/
  Double I_NAI_G4x_F3x_a = I_NAI_H5x_D2x_a+ABX*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F3x_a = I_NAI_H4xy_D2x_a+ABX*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F3x_a = I_NAI_H4xz_D2x_a+ABX*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F3x_a = I_NAI_H3x2y_D2x_a+ABX*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F3x_a = I_NAI_H3xyz_D2x_a+ABX*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F3x_a = I_NAI_H3x2z_D2x_a+ABX*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F3x_a = I_NAI_H2x3y_D2x_a+ABX*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F3x_a = I_NAI_H2x2yz_D2x_a+ABX*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F3x_a = I_NAI_H2xy2z_D2x_a+ABX*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F3x_a = I_NAI_H2x3z_D2x_a+ABX*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F3x_a = I_NAI_Hx4y_D2x_a+ABX*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F3x_a = I_NAI_Hx3yz_D2x_a+ABX*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F3x_a = I_NAI_Hx2y2z_D2x_a+ABX*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F3x_a = I_NAI_Hxy3z_D2x_a+ABX*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F3x_a = I_NAI_Hx4z_D2x_a+ABX*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_F2xy_a = I_NAI_H4xy_D2x_a+ABY*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F2xy_a = I_NAI_H3x2y_D2x_a+ABY*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F2xy_a = I_NAI_H3xyz_D2x_a+ABY*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F2xy_a = I_NAI_H2x3y_D2x_a+ABY*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F2xy_a = I_NAI_H2x2yz_D2x_a+ABY*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F2xy_a = I_NAI_H2xy2z_D2x_a+ABY*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F2xy_a = I_NAI_Hx4y_D2x_a+ABY*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F2xy_a = I_NAI_Hx3yz_D2x_a+ABY*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F2xy_a = I_NAI_Hx2y2z_D2x_a+ABY*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F2xy_a = I_NAI_Hxy3z_D2x_a+ABY*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F2xy_a = I_NAI_H5y_D2x_a+ABY*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F2xy_a = I_NAI_H4yz_D2x_a+ABY*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F2xy_a = I_NAI_H3y2z_D2x_a+ABY*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F2xy_a = I_NAI_H2y3z_D2x_a+ABY*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F2xy_a = I_NAI_Hy4z_D2x_a+ABY*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_F2xz_a = I_NAI_H4xz_D2x_a+ABZ*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F2xz_a = I_NAI_H3xyz_D2x_a+ABZ*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F2xz_a = I_NAI_H3x2z_D2x_a+ABZ*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F2xz_a = I_NAI_H2x2yz_D2x_a+ABZ*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F2xz_a = I_NAI_H2xy2z_D2x_a+ABZ*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F2xz_a = I_NAI_H2x3z_D2x_a+ABZ*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F2xz_a = I_NAI_Hx3yz_D2x_a+ABZ*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F2xz_a = I_NAI_Hx2y2z_D2x_a+ABZ*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F2xz_a = I_NAI_Hxy3z_D2x_a+ABZ*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F2xz_a = I_NAI_Hx4z_D2x_a+ABZ*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F2xz_a = I_NAI_H4yz_D2x_a+ABZ*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F2xz_a = I_NAI_H3y2z_D2x_a+ABZ*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F2xz_a = I_NAI_H2y3z_D2x_a+ABZ*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F2xz_a = I_NAI_Hy4z_D2x_a+ABZ*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F2xz_a = I_NAI_H5z_D2x_a+ABZ*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_Fx2y_a = I_NAI_H5x_D2y_a+ABX*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_Fx2y_a = I_NAI_H4xy_D2y_a+ABX*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_Fx2y_a = I_NAI_H4xz_D2y_a+ABX*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_Fx2y_a = I_NAI_H3x2y_D2y_a+ABX*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_Fx2y_a = I_NAI_H3xyz_D2y_a+ABX*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_Fx2y_a = I_NAI_H3x2z_D2y_a+ABX*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_Fx2y_a = I_NAI_H2x3y_D2y_a+ABX*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_Fx2y_a = I_NAI_H2x2yz_D2y_a+ABX*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_Fx2y_a = I_NAI_H2xy2z_D2y_a+ABX*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_Fx2y_a = I_NAI_H2x3z_D2y_a+ABX*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_Fx2y_a = I_NAI_Hx4y_D2y_a+ABX*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_Fx2y_a = I_NAI_Hx3yz_D2y_a+ABX*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_Fx2y_a = I_NAI_Hx2y2z_D2y_a+ABX*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_Fx2y_a = I_NAI_Hxy3z_D2y_a+ABX*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_Fx2y_a = I_NAI_Hx4z_D2y_a+ABX*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_Fxyz_a = I_NAI_H4xz_Dxy_a+ABZ*I_NAI_G4x_Dxy_a;
  Double I_NAI_G3xy_Fxyz_a = I_NAI_H3xyz_Dxy_a+ABZ*I_NAI_G3xy_Dxy_a;
  Double I_NAI_G3xz_Fxyz_a = I_NAI_H3x2z_Dxy_a+ABZ*I_NAI_G3xz_Dxy_a;
  Double I_NAI_G2x2y_Fxyz_a = I_NAI_H2x2yz_Dxy_a+ABZ*I_NAI_G2x2y_Dxy_a;
  Double I_NAI_G2xyz_Fxyz_a = I_NAI_H2xy2z_Dxy_a+ABZ*I_NAI_G2xyz_Dxy_a;
  Double I_NAI_G2x2z_Fxyz_a = I_NAI_H2x3z_Dxy_a+ABZ*I_NAI_G2x2z_Dxy_a;
  Double I_NAI_Gx3y_Fxyz_a = I_NAI_Hx3yz_Dxy_a+ABZ*I_NAI_Gx3y_Dxy_a;
  Double I_NAI_Gx2yz_Fxyz_a = I_NAI_Hx2y2z_Dxy_a+ABZ*I_NAI_Gx2yz_Dxy_a;
  Double I_NAI_Gxy2z_Fxyz_a = I_NAI_Hxy3z_Dxy_a+ABZ*I_NAI_Gxy2z_Dxy_a;
  Double I_NAI_Gx3z_Fxyz_a = I_NAI_Hx4z_Dxy_a+ABZ*I_NAI_Gx3z_Dxy_a;
  Double I_NAI_G4y_Fxyz_a = I_NAI_H4yz_Dxy_a+ABZ*I_NAI_G4y_Dxy_a;
  Double I_NAI_G3yz_Fxyz_a = I_NAI_H3y2z_Dxy_a+ABZ*I_NAI_G3yz_Dxy_a;
  Double I_NAI_G2y2z_Fxyz_a = I_NAI_H2y3z_Dxy_a+ABZ*I_NAI_G2y2z_Dxy_a;
  Double I_NAI_Gy3z_Fxyz_a = I_NAI_Hy4z_Dxy_a+ABZ*I_NAI_Gy3z_Dxy_a;
  Double I_NAI_G4z_Fxyz_a = I_NAI_H5z_Dxy_a+ABZ*I_NAI_G4z_Dxy_a;
  Double I_NAI_G4x_Fx2z_a = I_NAI_H5x_D2z_a+ABX*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_Fx2z_a = I_NAI_H4xy_D2z_a+ABX*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_Fx2z_a = I_NAI_H4xz_D2z_a+ABX*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_Fx2z_a = I_NAI_H3x2y_D2z_a+ABX*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_Fx2z_a = I_NAI_H3xyz_D2z_a+ABX*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_Fx2z_a = I_NAI_H3x2z_D2z_a+ABX*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_Fx2z_a = I_NAI_H2x3y_D2z_a+ABX*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_Fx2z_a = I_NAI_H2x2yz_D2z_a+ABX*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_Fx2z_a = I_NAI_H2xy2z_D2z_a+ABX*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_Fx2z_a = I_NAI_H2x3z_D2z_a+ABX*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_Fx2z_a = I_NAI_Hx4y_D2z_a+ABX*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_Fx2z_a = I_NAI_Hx3yz_D2z_a+ABX*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_Fx2z_a = I_NAI_Hx2y2z_D2z_a+ABX*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_Fx2z_a = I_NAI_Hxy3z_D2z_a+ABX*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_Fx2z_a = I_NAI_Hx4z_D2z_a+ABX*I_NAI_G4z_D2z_a;
  Double I_NAI_G4x_F3y_a = I_NAI_H4xy_D2y_a+ABY*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_F3y_a = I_NAI_H3x2y_D2y_a+ABY*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_F3y_a = I_NAI_H3xyz_D2y_a+ABY*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_F3y_a = I_NAI_H2x3y_D2y_a+ABY*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_F3y_a = I_NAI_H2x2yz_D2y_a+ABY*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_F3y_a = I_NAI_H2xy2z_D2y_a+ABY*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_F3y_a = I_NAI_Hx4y_D2y_a+ABY*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_F3y_a = I_NAI_Hx3yz_D2y_a+ABY*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_F3y_a = I_NAI_Hx2y2z_D2y_a+ABY*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_F3y_a = I_NAI_Hxy3z_D2y_a+ABY*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_F3y_a = I_NAI_H5y_D2y_a+ABY*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_F3y_a = I_NAI_H4yz_D2y_a+ABY*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_F3y_a = I_NAI_H3y2z_D2y_a+ABY*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_F3y_a = I_NAI_H2y3z_D2y_a+ABY*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_F3y_a = I_NAI_Hy4z_D2y_a+ABY*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_F2yz_a = I_NAI_H4xz_D2y_a+ABZ*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_F2yz_a = I_NAI_H3xyz_D2y_a+ABZ*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_F2yz_a = I_NAI_H3x2z_D2y_a+ABZ*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_F2yz_a = I_NAI_H2x2yz_D2y_a+ABZ*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_F2yz_a = I_NAI_H2xy2z_D2y_a+ABZ*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_F2yz_a = I_NAI_H2x3z_D2y_a+ABZ*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_F2yz_a = I_NAI_Hx3yz_D2y_a+ABZ*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_F2yz_a = I_NAI_Hx2y2z_D2y_a+ABZ*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_F2yz_a = I_NAI_Hxy3z_D2y_a+ABZ*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_F2yz_a = I_NAI_Hx4z_D2y_a+ABZ*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_F2yz_a = I_NAI_H4yz_D2y_a+ABZ*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_F2yz_a = I_NAI_H3y2z_D2y_a+ABZ*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_F2yz_a = I_NAI_H2y3z_D2y_a+ABZ*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_F2yz_a = I_NAI_Hy4z_D2y_a+ABZ*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_F2yz_a = I_NAI_H5z_D2y_a+ABZ*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_Fy2z_a = I_NAI_H4xy_D2z_a+ABY*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_Fy2z_a = I_NAI_H3x2y_D2z_a+ABY*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_Fy2z_a = I_NAI_H3xyz_D2z_a+ABY*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_Fy2z_a = I_NAI_H2x3y_D2z_a+ABY*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_Fy2z_a = I_NAI_H2x2yz_D2z_a+ABY*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_Fy2z_a = I_NAI_H2xy2z_D2z_a+ABY*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_Fy2z_a = I_NAI_Hx4y_D2z_a+ABY*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_Fy2z_a = I_NAI_Hx3yz_D2z_a+ABY*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_Fy2z_a = I_NAI_Hx2y2z_D2z_a+ABY*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_Fy2z_a = I_NAI_Hxy3z_D2z_a+ABY*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_Fy2z_a = I_NAI_H5y_D2z_a+ABY*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_Fy2z_a = I_NAI_H4yz_D2z_a+ABY*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_Fy2z_a = I_NAI_H3y2z_D2z_a+ABY*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_Fy2z_a = I_NAI_H2y3z_D2z_a+ABY*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_Fy2z_a = I_NAI_Hy4z_D2z_a+ABY*I_NAI_G4z_D2z_a;
  Double I_NAI_G4x_F3z_a = I_NAI_H4xz_D2z_a+ABZ*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_F3z_a = I_NAI_H3xyz_D2z_a+ABZ*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_F3z_a = I_NAI_H3x2z_D2z_a+ABZ*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_F3z_a = I_NAI_H2x2yz_D2z_a+ABZ*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_F3z_a = I_NAI_H2xy2z_D2z_a+ABZ*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_F3z_a = I_NAI_H2x3z_D2z_a+ABZ*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_F3z_a = I_NAI_Hx3yz_D2z_a+ABZ*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_F3z_a = I_NAI_Hx2y2z_D2z_a+ABZ*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_F3z_a = I_NAI_Hxy3z_D2z_a+ABZ*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_F3z_a = I_NAI_Hx4z_D2z_a+ABZ*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_F3z_a = I_NAI_H4yz_D2z_a+ABZ*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_F3z_a = I_NAI_H3y2z_D2z_a+ABZ*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_F3z_a = I_NAI_H2y3z_D2z_a+ABZ*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_F3z_a = I_NAI_Hy4z_D2z_a+ABZ*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_F3z_a = I_NAI_H5z_D2z_a+ABZ*I_NAI_G4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_b
   * RHS shell quartet name: SQ_NAI_F_S_b
   ************************************************************/
  Double I_NAI_F3x_Px_b = I_NAI_G4x_S_b+ABX*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Px_b = I_NAI_G3xy_S_b+ABX*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Px_b = I_NAI_G3xz_S_b+ABX*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Px_b = I_NAI_G2x2y_S_b+ABX*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Px_b = I_NAI_G2xyz_S_b+ABX*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Px_b = I_NAI_G2x2z_S_b+ABX*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Px_b = I_NAI_Gx3y_S_b+ABX*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Px_b = I_NAI_Gx2yz_S_b+ABX*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Px_b = I_NAI_Gxy2z_S_b+ABX*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Px_b = I_NAI_Gx3z_S_b+ABX*I_NAI_F3z_S_b;
  Double I_NAI_F3x_Py_b = I_NAI_G3xy_S_b+ABY*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Py_b = I_NAI_G2x2y_S_b+ABY*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Py_b = I_NAI_G2xyz_S_b+ABY*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Py_b = I_NAI_Gx3y_S_b+ABY*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Py_b = I_NAI_Gx2yz_S_b+ABY*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Py_b = I_NAI_Gxy2z_S_b+ABY*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Py_b = I_NAI_G4y_S_b+ABY*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Py_b = I_NAI_G3yz_S_b+ABY*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Py_b = I_NAI_G2y2z_S_b+ABY*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Py_b = I_NAI_Gy3z_S_b+ABY*I_NAI_F3z_S_b;
  Double I_NAI_F3x_Pz_b = I_NAI_G3xz_S_b+ABZ*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Pz_b = I_NAI_G2xyz_S_b+ABZ*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Pz_b = I_NAI_G2x2z_S_b+ABZ*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Pz_b = I_NAI_Gx2yz_S_b+ABZ*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Pz_b = I_NAI_Gxy2z_S_b+ABZ*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Pz_b = I_NAI_Gx3z_S_b+ABZ*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Pz_b = I_NAI_G3yz_S_b+ABZ*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Pz_b = I_NAI_G2y2z_S_b+ABZ*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Pz_b = I_NAI_Gy3z_S_b+ABZ*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Pz_b = I_NAI_G4z_S_b+ABZ*I_NAI_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_b
   * RHS shell quartet name: SQ_NAI_G_S_b
   ************************************************************/
  Double I_NAI_G4x_Px_b = I_NAI_H5x_S_b+ABX*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Px_b = I_NAI_H4xy_S_b+ABX*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Px_b = I_NAI_H4xz_S_b+ABX*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Px_b = I_NAI_H3x2y_S_b+ABX*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Px_b = I_NAI_H3xyz_S_b+ABX*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Px_b = I_NAI_H3x2z_S_b+ABX*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Px_b = I_NAI_H2x3y_S_b+ABX*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Px_b = I_NAI_H2x2yz_S_b+ABX*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Px_b = I_NAI_H2xy2z_S_b+ABX*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Px_b = I_NAI_H2x3z_S_b+ABX*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Px_b = I_NAI_Hx4y_S_b+ABX*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Px_b = I_NAI_Hx3yz_S_b+ABX*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Px_b = I_NAI_Hx2y2z_S_b+ABX*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Px_b = I_NAI_Hxy3z_S_b+ABX*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Px_b = I_NAI_Hx4z_S_b+ABX*I_NAI_G4z_S_b;
  Double I_NAI_G4x_Py_b = I_NAI_H4xy_S_b+ABY*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Py_b = I_NAI_H3x2y_S_b+ABY*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Py_b = I_NAI_H3xyz_S_b+ABY*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Py_b = I_NAI_H2x3y_S_b+ABY*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Py_b = I_NAI_H2x2yz_S_b+ABY*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Py_b = I_NAI_H2xy2z_S_b+ABY*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Py_b = I_NAI_Hx4y_S_b+ABY*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Py_b = I_NAI_Hx3yz_S_b+ABY*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Py_b = I_NAI_Hx2y2z_S_b+ABY*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Py_b = I_NAI_Hxy3z_S_b+ABY*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Py_b = I_NAI_H5y_S_b+ABY*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Py_b = I_NAI_H4yz_S_b+ABY*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Py_b = I_NAI_H3y2z_S_b+ABY*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Py_b = I_NAI_H2y3z_S_b+ABY*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Py_b = I_NAI_Hy4z_S_b+ABY*I_NAI_G4z_S_b;
  Double I_NAI_G4x_Pz_b = I_NAI_H4xz_S_b+ABZ*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Pz_b = I_NAI_H3xyz_S_b+ABZ*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Pz_b = I_NAI_H3x2z_S_b+ABZ*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Pz_b = I_NAI_H2x2yz_S_b+ABZ*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Pz_b = I_NAI_H2xy2z_S_b+ABZ*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Pz_b = I_NAI_H2x3z_S_b+ABZ*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Pz_b = I_NAI_Hx3yz_S_b+ABZ*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Pz_b = I_NAI_Hx2y2z_S_b+ABZ*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Pz_b = I_NAI_Hxy3z_S_b+ABZ*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Pz_b = I_NAI_Hx4z_S_b+ABZ*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Pz_b = I_NAI_H4yz_S_b+ABZ*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Pz_b = I_NAI_H3y2z_S_b+ABZ*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Pz_b = I_NAI_H2y3z_S_b+ABZ*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Pz_b = I_NAI_Hy4z_S_b+ABZ*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Pz_b = I_NAI_H5z_S_b+ABZ*I_NAI_G4z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_b
   * RHS shell quartet name: SQ_NAI_F_P_b
   ************************************************************/
  Double I_NAI_F3x_D2x_b = I_NAI_G4x_Px_b+ABX*I_NAI_F3x_Px_b;
  Double I_NAI_F2xy_D2x_b = I_NAI_G3xy_Px_b+ABX*I_NAI_F2xy_Px_b;
  Double I_NAI_F2xz_D2x_b = I_NAI_G3xz_Px_b+ABX*I_NAI_F2xz_Px_b;
  Double I_NAI_Fx2y_D2x_b = I_NAI_G2x2y_Px_b+ABX*I_NAI_Fx2y_Px_b;
  Double I_NAI_Fxyz_D2x_b = I_NAI_G2xyz_Px_b+ABX*I_NAI_Fxyz_Px_b;
  Double I_NAI_Fx2z_D2x_b = I_NAI_G2x2z_Px_b+ABX*I_NAI_Fx2z_Px_b;
  Double I_NAI_F3y_D2x_b = I_NAI_Gx3y_Px_b+ABX*I_NAI_F3y_Px_b;
  Double I_NAI_F2yz_D2x_b = I_NAI_Gx2yz_Px_b+ABX*I_NAI_F2yz_Px_b;
  Double I_NAI_Fy2z_D2x_b = I_NAI_Gxy2z_Px_b+ABX*I_NAI_Fy2z_Px_b;
  Double I_NAI_F3z_D2x_b = I_NAI_Gx3z_Px_b+ABX*I_NAI_F3z_Px_b;
  Double I_NAI_F3x_D2y_b = I_NAI_G3xy_Py_b+ABY*I_NAI_F3x_Py_b;
  Double I_NAI_F2xy_D2y_b = I_NAI_G2x2y_Py_b+ABY*I_NAI_F2xy_Py_b;
  Double I_NAI_F2xz_D2y_b = I_NAI_G2xyz_Py_b+ABY*I_NAI_F2xz_Py_b;
  Double I_NAI_Fx2y_D2y_b = I_NAI_Gx3y_Py_b+ABY*I_NAI_Fx2y_Py_b;
  Double I_NAI_Fxyz_D2y_b = I_NAI_Gx2yz_Py_b+ABY*I_NAI_Fxyz_Py_b;
  Double I_NAI_Fx2z_D2y_b = I_NAI_Gxy2z_Py_b+ABY*I_NAI_Fx2z_Py_b;
  Double I_NAI_F3y_D2y_b = I_NAI_G4y_Py_b+ABY*I_NAI_F3y_Py_b;
  Double I_NAI_F2yz_D2y_b = I_NAI_G3yz_Py_b+ABY*I_NAI_F2yz_Py_b;
  Double I_NAI_Fy2z_D2y_b = I_NAI_G2y2z_Py_b+ABY*I_NAI_Fy2z_Py_b;
  Double I_NAI_F3z_D2y_b = I_NAI_Gy3z_Py_b+ABY*I_NAI_F3z_Py_b;
  Double I_NAI_F3x_D2z_b = I_NAI_G3xz_Pz_b+ABZ*I_NAI_F3x_Pz_b;
  Double I_NAI_F2xy_D2z_b = I_NAI_G2xyz_Pz_b+ABZ*I_NAI_F2xy_Pz_b;
  Double I_NAI_F2xz_D2z_b = I_NAI_G2x2z_Pz_b+ABZ*I_NAI_F2xz_Pz_b;
  Double I_NAI_Fx2y_D2z_b = I_NAI_Gx2yz_Pz_b+ABZ*I_NAI_Fx2y_Pz_b;
  Double I_NAI_Fxyz_D2z_b = I_NAI_Gxy2z_Pz_b+ABZ*I_NAI_Fxyz_Pz_b;
  Double I_NAI_Fx2z_D2z_b = I_NAI_Gx3z_Pz_b+ABZ*I_NAI_Fx2z_Pz_b;
  Double I_NAI_F3y_D2z_b = I_NAI_G3yz_Pz_b+ABZ*I_NAI_F3y_Pz_b;
  Double I_NAI_F2yz_D2z_b = I_NAI_G2y2z_Pz_b+ABZ*I_NAI_F2yz_Pz_b;
  Double I_NAI_Fy2z_D2z_b = I_NAI_Gy3z_Pz_b+ABZ*I_NAI_Fy2z_Pz_b;
  Double I_NAI_F3z_D2z_b = I_NAI_G4z_Pz_b+ABZ*I_NAI_F3z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_b
   * RHS shell quartet name: SQ_NAI_H_S_b
   ************************************************************/
  Double I_NAI_H5x_Px_b = I_NAI_I6x_S_b+ABX*I_NAI_H5x_S_b;
  Double I_NAI_H4xy_Px_b = I_NAI_I5xy_S_b+ABX*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Px_b = I_NAI_I5xz_S_b+ABX*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Px_b = I_NAI_I4x2y_S_b+ABX*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Px_b = I_NAI_I4xyz_S_b+ABX*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Px_b = I_NAI_I4x2z_S_b+ABX*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Px_b = I_NAI_I3x3y_S_b+ABX*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Px_b = I_NAI_I3x2yz_S_b+ABX*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Px_b = I_NAI_I3xy2z_S_b+ABX*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Px_b = I_NAI_I3x3z_S_b+ABX*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Px_b = I_NAI_I2x4y_S_b+ABX*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Px_b = I_NAI_I2x3yz_S_b+ABX*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Px_b = I_NAI_I2x2y2z_S_b+ABX*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Px_b = I_NAI_I2xy3z_S_b+ABX*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Px_b = I_NAI_I2x4z_S_b+ABX*I_NAI_Hx4z_S_b;
  Double I_NAI_H5y_Px_b = I_NAI_Ix5y_S_b+ABX*I_NAI_H5y_S_b;
  Double I_NAI_H4yz_Px_b = I_NAI_Ix4yz_S_b+ABX*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Px_b = I_NAI_Ix3y2z_S_b+ABX*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Px_b = I_NAI_Ix2y3z_S_b+ABX*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Px_b = I_NAI_Ixy4z_S_b+ABX*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Px_b = I_NAI_Ix5z_S_b+ABX*I_NAI_H5z_S_b;
  Double I_NAI_H5x_Py_b = I_NAI_I5xy_S_b+ABY*I_NAI_H5x_S_b;
  Double I_NAI_H4xy_Py_b = I_NAI_I4x2y_S_b+ABY*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Py_b = I_NAI_I4xyz_S_b+ABY*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Py_b = I_NAI_I3x3y_S_b+ABY*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Py_b = I_NAI_I3x2yz_S_b+ABY*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Py_b = I_NAI_I3xy2z_S_b+ABY*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Py_b = I_NAI_I2x4y_S_b+ABY*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Py_b = I_NAI_I2x3yz_S_b+ABY*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Py_b = I_NAI_I2x2y2z_S_b+ABY*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Py_b = I_NAI_I2xy3z_S_b+ABY*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Py_b = I_NAI_Ix5y_S_b+ABY*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Py_b = I_NAI_Ix4yz_S_b+ABY*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Py_b = I_NAI_Ix3y2z_S_b+ABY*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Py_b = I_NAI_Ix2y3z_S_b+ABY*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Py_b = I_NAI_Ixy4z_S_b+ABY*I_NAI_Hx4z_S_b;
  Double I_NAI_H5y_Py_b = I_NAI_I6y_S_b+ABY*I_NAI_H5y_S_b;
  Double I_NAI_H4yz_Py_b = I_NAI_I5yz_S_b+ABY*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Py_b = I_NAI_I4y2z_S_b+ABY*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Py_b = I_NAI_I3y3z_S_b+ABY*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Py_b = I_NAI_I2y4z_S_b+ABY*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Py_b = I_NAI_Iy5z_S_b+ABY*I_NAI_H5z_S_b;
  Double I_NAI_H5x_Pz_b = I_NAI_I5xz_S_b+ABZ*I_NAI_H5x_S_b;
  Double I_NAI_H4xy_Pz_b = I_NAI_I4xyz_S_b+ABZ*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Pz_b = I_NAI_I4x2z_S_b+ABZ*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Pz_b = I_NAI_I3x2yz_S_b+ABZ*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Pz_b = I_NAI_I3xy2z_S_b+ABZ*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Pz_b = I_NAI_I3x3z_S_b+ABZ*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Pz_b = I_NAI_I2x3yz_S_b+ABZ*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Pz_b = I_NAI_I2x2y2z_S_b+ABZ*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Pz_b = I_NAI_I2xy3z_S_b+ABZ*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Pz_b = I_NAI_I2x4z_S_b+ABZ*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Pz_b = I_NAI_Ix4yz_S_b+ABZ*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Pz_b = I_NAI_Ix3y2z_S_b+ABZ*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Pz_b = I_NAI_Ix2y3z_S_b+ABZ*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Pz_b = I_NAI_Ixy4z_S_b+ABZ*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Pz_b = I_NAI_Ix5z_S_b+ABZ*I_NAI_Hx4z_S_b;
  Double I_NAI_H5y_Pz_b = I_NAI_I5yz_S_b+ABZ*I_NAI_H5y_S_b;
  Double I_NAI_H4yz_Pz_b = I_NAI_I4y2z_S_b+ABZ*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Pz_b = I_NAI_I3y3z_S_b+ABZ*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Pz_b = I_NAI_I2y4z_S_b+ABZ*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Pz_b = I_NAI_Iy5z_S_b+ABZ*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Pz_b = I_NAI_I6z_S_b+ABZ*I_NAI_H5z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_G_P_b
   ************************************************************/
  Double I_NAI_G4x_D2x_b = I_NAI_H5x_Px_b+ABX*I_NAI_G4x_Px_b;
  Double I_NAI_G3xy_D2x_b = I_NAI_H4xy_Px_b+ABX*I_NAI_G3xy_Px_b;
  Double I_NAI_G3xz_D2x_b = I_NAI_H4xz_Px_b+ABX*I_NAI_G3xz_Px_b;
  Double I_NAI_G2x2y_D2x_b = I_NAI_H3x2y_Px_b+ABX*I_NAI_G2x2y_Px_b;
  Double I_NAI_G2xyz_D2x_b = I_NAI_H3xyz_Px_b+ABX*I_NAI_G2xyz_Px_b;
  Double I_NAI_G2x2z_D2x_b = I_NAI_H3x2z_Px_b+ABX*I_NAI_G2x2z_Px_b;
  Double I_NAI_Gx3y_D2x_b = I_NAI_H2x3y_Px_b+ABX*I_NAI_Gx3y_Px_b;
  Double I_NAI_Gx2yz_D2x_b = I_NAI_H2x2yz_Px_b+ABX*I_NAI_Gx2yz_Px_b;
  Double I_NAI_Gxy2z_D2x_b = I_NAI_H2xy2z_Px_b+ABX*I_NAI_Gxy2z_Px_b;
  Double I_NAI_Gx3z_D2x_b = I_NAI_H2x3z_Px_b+ABX*I_NAI_Gx3z_Px_b;
  Double I_NAI_G4y_D2x_b = I_NAI_Hx4y_Px_b+ABX*I_NAI_G4y_Px_b;
  Double I_NAI_G3yz_D2x_b = I_NAI_Hx3yz_Px_b+ABX*I_NAI_G3yz_Px_b;
  Double I_NAI_G2y2z_D2x_b = I_NAI_Hx2y2z_Px_b+ABX*I_NAI_G2y2z_Px_b;
  Double I_NAI_Gy3z_D2x_b = I_NAI_Hxy3z_Px_b+ABX*I_NAI_Gy3z_Px_b;
  Double I_NAI_G4z_D2x_b = I_NAI_Hx4z_Px_b+ABX*I_NAI_G4z_Px_b;
  Double I_NAI_G4x_Dxy_b = I_NAI_H4xy_Px_b+ABY*I_NAI_G4x_Px_b;
  Double I_NAI_G3xy_Dxy_b = I_NAI_H3x2y_Px_b+ABY*I_NAI_G3xy_Px_b;
  Double I_NAI_G3xz_Dxy_b = I_NAI_H3xyz_Px_b+ABY*I_NAI_G3xz_Px_b;
  Double I_NAI_G2x2y_Dxy_b = I_NAI_H2x3y_Px_b+ABY*I_NAI_G2x2y_Px_b;
  Double I_NAI_G2xyz_Dxy_b = I_NAI_H2x2yz_Px_b+ABY*I_NAI_G2xyz_Px_b;
  Double I_NAI_G2x2z_Dxy_b = I_NAI_H2xy2z_Px_b+ABY*I_NAI_G2x2z_Px_b;
  Double I_NAI_Gx3y_Dxy_b = I_NAI_Hx4y_Px_b+ABY*I_NAI_Gx3y_Px_b;
  Double I_NAI_Gx2yz_Dxy_b = I_NAI_Hx3yz_Px_b+ABY*I_NAI_Gx2yz_Px_b;
  Double I_NAI_Gxy2z_Dxy_b = I_NAI_Hx2y2z_Px_b+ABY*I_NAI_Gxy2z_Px_b;
  Double I_NAI_Gx3z_Dxy_b = I_NAI_Hxy3z_Px_b+ABY*I_NAI_Gx3z_Px_b;
  Double I_NAI_G4y_Dxy_b = I_NAI_H5y_Px_b+ABY*I_NAI_G4y_Px_b;
  Double I_NAI_G3yz_Dxy_b = I_NAI_H4yz_Px_b+ABY*I_NAI_G3yz_Px_b;
  Double I_NAI_G2y2z_Dxy_b = I_NAI_H3y2z_Px_b+ABY*I_NAI_G2y2z_Px_b;
  Double I_NAI_Gy3z_Dxy_b = I_NAI_H2y3z_Px_b+ABY*I_NAI_Gy3z_Px_b;
  Double I_NAI_G4z_Dxy_b = I_NAI_Hy4z_Px_b+ABY*I_NAI_G4z_Px_b;
  Double I_NAI_G4x_D2y_b = I_NAI_H4xy_Py_b+ABY*I_NAI_G4x_Py_b;
  Double I_NAI_G3xy_D2y_b = I_NAI_H3x2y_Py_b+ABY*I_NAI_G3xy_Py_b;
  Double I_NAI_G3xz_D2y_b = I_NAI_H3xyz_Py_b+ABY*I_NAI_G3xz_Py_b;
  Double I_NAI_G2x2y_D2y_b = I_NAI_H2x3y_Py_b+ABY*I_NAI_G2x2y_Py_b;
  Double I_NAI_G2xyz_D2y_b = I_NAI_H2x2yz_Py_b+ABY*I_NAI_G2xyz_Py_b;
  Double I_NAI_G2x2z_D2y_b = I_NAI_H2xy2z_Py_b+ABY*I_NAI_G2x2z_Py_b;
  Double I_NAI_Gx3y_D2y_b = I_NAI_Hx4y_Py_b+ABY*I_NAI_Gx3y_Py_b;
  Double I_NAI_Gx2yz_D2y_b = I_NAI_Hx3yz_Py_b+ABY*I_NAI_Gx2yz_Py_b;
  Double I_NAI_Gxy2z_D2y_b = I_NAI_Hx2y2z_Py_b+ABY*I_NAI_Gxy2z_Py_b;
  Double I_NAI_Gx3z_D2y_b = I_NAI_Hxy3z_Py_b+ABY*I_NAI_Gx3z_Py_b;
  Double I_NAI_G4y_D2y_b = I_NAI_H5y_Py_b+ABY*I_NAI_G4y_Py_b;
  Double I_NAI_G3yz_D2y_b = I_NAI_H4yz_Py_b+ABY*I_NAI_G3yz_Py_b;
  Double I_NAI_G2y2z_D2y_b = I_NAI_H3y2z_Py_b+ABY*I_NAI_G2y2z_Py_b;
  Double I_NAI_Gy3z_D2y_b = I_NAI_H2y3z_Py_b+ABY*I_NAI_Gy3z_Py_b;
  Double I_NAI_G4z_D2y_b = I_NAI_Hy4z_Py_b+ABY*I_NAI_G4z_Py_b;
  Double I_NAI_G4x_D2z_b = I_NAI_H4xz_Pz_b+ABZ*I_NAI_G4x_Pz_b;
  Double I_NAI_G3xy_D2z_b = I_NAI_H3xyz_Pz_b+ABZ*I_NAI_G3xy_Pz_b;
  Double I_NAI_G3xz_D2z_b = I_NAI_H3x2z_Pz_b+ABZ*I_NAI_G3xz_Pz_b;
  Double I_NAI_G2x2y_D2z_b = I_NAI_H2x2yz_Pz_b+ABZ*I_NAI_G2x2y_Pz_b;
  Double I_NAI_G2xyz_D2z_b = I_NAI_H2xy2z_Pz_b+ABZ*I_NAI_G2xyz_Pz_b;
  Double I_NAI_G2x2z_D2z_b = I_NAI_H2x3z_Pz_b+ABZ*I_NAI_G2x2z_Pz_b;
  Double I_NAI_Gx3y_D2z_b = I_NAI_Hx3yz_Pz_b+ABZ*I_NAI_Gx3y_Pz_b;
  Double I_NAI_Gx2yz_D2z_b = I_NAI_Hx2y2z_Pz_b+ABZ*I_NAI_Gx2yz_Pz_b;
  Double I_NAI_Gxy2z_D2z_b = I_NAI_Hxy3z_Pz_b+ABZ*I_NAI_Gxy2z_Pz_b;
  Double I_NAI_Gx3z_D2z_b = I_NAI_Hx4z_Pz_b+ABZ*I_NAI_Gx3z_Pz_b;
  Double I_NAI_G4y_D2z_b = I_NAI_H4yz_Pz_b+ABZ*I_NAI_G4y_Pz_b;
  Double I_NAI_G3yz_D2z_b = I_NAI_H3y2z_Pz_b+ABZ*I_NAI_G3yz_Pz_b;
  Double I_NAI_G2y2z_D2z_b = I_NAI_H2y3z_Pz_b+ABZ*I_NAI_G2y2z_Pz_b;
  Double I_NAI_Gy3z_D2z_b = I_NAI_Hy4z_Pz_b+ABZ*I_NAI_Gy3z_Pz_b;
  Double I_NAI_G4z_D2z_b = I_NAI_H5z_Pz_b+ABZ*I_NAI_G4z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_F_D_b
   ************************************************************/
  Double I_NAI_F3x_F3x_b = I_NAI_G4x_D2x_b+ABX*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F3x_b = I_NAI_G3xy_D2x_b+ABX*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F3x_b = I_NAI_G3xz_D2x_b+ABX*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F3x_b = I_NAI_G2x2y_D2x_b+ABX*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F3x_b = I_NAI_G2xyz_D2x_b+ABX*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F3x_b = I_NAI_G2x2z_D2x_b+ABX*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F3x_b = I_NAI_Gx3y_D2x_b+ABX*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F3x_b = I_NAI_Gx2yz_D2x_b+ABX*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F3x_b = I_NAI_Gxy2z_D2x_b+ABX*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F3x_b = I_NAI_Gx3z_D2x_b+ABX*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_F2xy_b = I_NAI_G3xy_D2x_b+ABY*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F2xy_b = I_NAI_G2x2y_D2x_b+ABY*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F2xy_b = I_NAI_G2xyz_D2x_b+ABY*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F2xy_b = I_NAI_Gx3y_D2x_b+ABY*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F2xy_b = I_NAI_Gx2yz_D2x_b+ABY*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F2xy_b = I_NAI_Gxy2z_D2x_b+ABY*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F2xy_b = I_NAI_G4y_D2x_b+ABY*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F2xy_b = I_NAI_G3yz_D2x_b+ABY*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F2xy_b = I_NAI_G2y2z_D2x_b+ABY*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F2xy_b = I_NAI_Gy3z_D2x_b+ABY*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_F2xz_b = I_NAI_G3xz_D2x_b+ABZ*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F2xz_b = I_NAI_G2xyz_D2x_b+ABZ*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F2xz_b = I_NAI_G2x2z_D2x_b+ABZ*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F2xz_b = I_NAI_Gx2yz_D2x_b+ABZ*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F2xz_b = I_NAI_Gxy2z_D2x_b+ABZ*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F2xz_b = I_NAI_Gx3z_D2x_b+ABZ*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F2xz_b = I_NAI_G3yz_D2x_b+ABZ*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F2xz_b = I_NAI_G2y2z_D2x_b+ABZ*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F2xz_b = I_NAI_Gy3z_D2x_b+ABZ*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F2xz_b = I_NAI_G4z_D2x_b+ABZ*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_Fx2y_b = I_NAI_G4x_D2y_b+ABX*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_Fx2y_b = I_NAI_G3xy_D2y_b+ABX*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_Fx2y_b = I_NAI_G3xz_D2y_b+ABX*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_Fx2y_b = I_NAI_G2x2y_D2y_b+ABX*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_Fx2y_b = I_NAI_G2xyz_D2y_b+ABX*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_Fx2y_b = I_NAI_G2x2z_D2y_b+ABX*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_Fx2y_b = I_NAI_Gx3y_D2y_b+ABX*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_Fx2y_b = I_NAI_Gx2yz_D2y_b+ABX*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_Fx2y_b = I_NAI_Gxy2z_D2y_b+ABX*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_Fx2y_b = I_NAI_Gx3z_D2y_b+ABX*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_Fx2z_b = I_NAI_G4x_D2z_b+ABX*I_NAI_F3x_D2z_b;
  Double I_NAI_F2xy_Fx2z_b = I_NAI_G3xy_D2z_b+ABX*I_NAI_F2xy_D2z_b;
  Double I_NAI_F2xz_Fx2z_b = I_NAI_G3xz_D2z_b+ABX*I_NAI_F2xz_D2z_b;
  Double I_NAI_Fx2y_Fx2z_b = I_NAI_G2x2y_D2z_b+ABX*I_NAI_Fx2y_D2z_b;
  Double I_NAI_Fxyz_Fx2z_b = I_NAI_G2xyz_D2z_b+ABX*I_NAI_Fxyz_D2z_b;
  Double I_NAI_Fx2z_Fx2z_b = I_NAI_G2x2z_D2z_b+ABX*I_NAI_Fx2z_D2z_b;
  Double I_NAI_F3y_Fx2z_b = I_NAI_Gx3y_D2z_b+ABX*I_NAI_F3y_D2z_b;
  Double I_NAI_F2yz_Fx2z_b = I_NAI_Gx2yz_D2z_b+ABX*I_NAI_F2yz_D2z_b;
  Double I_NAI_Fy2z_Fx2z_b = I_NAI_Gxy2z_D2z_b+ABX*I_NAI_Fy2z_D2z_b;
  Double I_NAI_F3z_Fx2z_b = I_NAI_Gx3z_D2z_b+ABX*I_NAI_F3z_D2z_b;
  Double I_NAI_F3x_F3y_b = I_NAI_G3xy_D2y_b+ABY*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_F3y_b = I_NAI_G2x2y_D2y_b+ABY*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_F3y_b = I_NAI_G2xyz_D2y_b+ABY*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_F3y_b = I_NAI_Gx3y_D2y_b+ABY*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_F3y_b = I_NAI_Gx2yz_D2y_b+ABY*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_F3y_b = I_NAI_Gxy2z_D2y_b+ABY*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_F3y_b = I_NAI_G4y_D2y_b+ABY*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_F3y_b = I_NAI_G3yz_D2y_b+ABY*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_F3y_b = I_NAI_G2y2z_D2y_b+ABY*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_F3y_b = I_NAI_Gy3z_D2y_b+ABY*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_F2yz_b = I_NAI_G3xz_D2y_b+ABZ*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_F2yz_b = I_NAI_G2xyz_D2y_b+ABZ*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_F2yz_b = I_NAI_G2x2z_D2y_b+ABZ*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_F2yz_b = I_NAI_Gx2yz_D2y_b+ABZ*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_F2yz_b = I_NAI_Gxy2z_D2y_b+ABZ*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_F2yz_b = I_NAI_Gx3z_D2y_b+ABZ*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_F2yz_b = I_NAI_G3yz_D2y_b+ABZ*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_F2yz_b = I_NAI_G2y2z_D2y_b+ABZ*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_F2yz_b = I_NAI_Gy3z_D2y_b+ABZ*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_F2yz_b = I_NAI_G4z_D2y_b+ABZ*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_F3z_b = I_NAI_G3xz_D2z_b+ABZ*I_NAI_F3x_D2z_b;
  Double I_NAI_F2xy_F3z_b = I_NAI_G2xyz_D2z_b+ABZ*I_NAI_F2xy_D2z_b;
  Double I_NAI_F2xz_F3z_b = I_NAI_G2x2z_D2z_b+ABZ*I_NAI_F2xz_D2z_b;
  Double I_NAI_Fx2y_F3z_b = I_NAI_Gx2yz_D2z_b+ABZ*I_NAI_Fx2y_D2z_b;
  Double I_NAI_Fxyz_F3z_b = I_NAI_Gxy2z_D2z_b+ABZ*I_NAI_Fxyz_D2z_b;
  Double I_NAI_Fx2z_F3z_b = I_NAI_Gx3z_D2z_b+ABZ*I_NAI_Fx2z_D2z_b;
  Double I_NAI_F3y_F3z_b = I_NAI_G3yz_D2z_b+ABZ*I_NAI_F3y_D2z_b;
  Double I_NAI_F2yz_F3z_b = I_NAI_G2y2z_D2z_b+ABZ*I_NAI_F2yz_D2z_b;
  Double I_NAI_Fy2z_F3z_b = I_NAI_Gy3z_D2z_b+ABZ*I_NAI_Fy2z_D2z_b;
  Double I_NAI_F3z_F3z_b = I_NAI_G4z_D2z_b+ABZ*I_NAI_F3z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 16 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_b
   * RHS shell quartet name: SQ_NAI_I_S_b
   ************************************************************/
  Double I_NAI_I6x_Px_b = I_NAI_K7x_S_b+ABX*I_NAI_I6x_S_b;
  Double I_NAI_I5xy_Px_b = I_NAI_K6xy_S_b+ABX*I_NAI_I5xy_S_b;
  Double I_NAI_I5xz_Px_b = I_NAI_K6xz_S_b+ABX*I_NAI_I5xz_S_b;
  Double I_NAI_I4x2y_Px_b = I_NAI_K5x2y_S_b+ABX*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Px_b = I_NAI_K5xyz_S_b+ABX*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Px_b = I_NAI_K5x2z_S_b+ABX*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x3y_Px_b = I_NAI_K4x3y_S_b+ABX*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Px_b = I_NAI_K4x2yz_S_b+ABX*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Px_b = I_NAI_K4xy2z_S_b+ABX*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Px_b = I_NAI_K4x3z_S_b+ABX*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x4y_Px_b = I_NAI_K3x4y_S_b+ABX*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Px_b = I_NAI_K3x3yz_S_b+ABX*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Px_b = I_NAI_K3x2y2z_S_b+ABX*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Px_b = I_NAI_K3xy3z_S_b+ABX*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Px_b = I_NAI_K3x4z_S_b+ABX*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix5y_Px_b = I_NAI_K2x5y_S_b+ABX*I_NAI_Ix5y_S_b;
  Double I_NAI_Ix4yz_Px_b = I_NAI_K2x4yz_S_b+ABX*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Px_b = I_NAI_K2x3y2z_S_b+ABX*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Px_b = I_NAI_K2x2y3z_S_b+ABX*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Px_b = I_NAI_K2xy4z_S_b+ABX*I_NAI_Ixy4z_S_b;
  Double I_NAI_Ix5z_Px_b = I_NAI_K2x5z_S_b+ABX*I_NAI_Ix5z_S_b;
  Double I_NAI_I5yz_Px_b = I_NAI_Kx5yz_S_b+ABX*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Px_b = I_NAI_Kx4y2z_S_b+ABX*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Px_b = I_NAI_Kx3y3z_S_b+ABX*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Px_b = I_NAI_Kx2y4z_S_b+ABX*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Px_b = I_NAI_Kxy5z_S_b+ABX*I_NAI_Iy5z_S_b;
  Double I_NAI_I5xy_Py_b = I_NAI_K5x2y_S_b+ABY*I_NAI_I5xy_S_b;
  Double I_NAI_I4x2y_Py_b = I_NAI_K4x3y_S_b+ABY*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Py_b = I_NAI_K4x2yz_S_b+ABY*I_NAI_I4xyz_S_b;
  Double I_NAI_I3x3y_Py_b = I_NAI_K3x4y_S_b+ABY*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Py_b = I_NAI_K3x3yz_S_b+ABY*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Py_b = I_NAI_K3x2y2z_S_b+ABY*I_NAI_I3xy2z_S_b;
  Double I_NAI_I2x4y_Py_b = I_NAI_K2x5y_S_b+ABY*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Py_b = I_NAI_K2x4yz_S_b+ABY*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Py_b = I_NAI_K2x3y2z_S_b+ABY*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Py_b = I_NAI_K2x2y3z_S_b+ABY*I_NAI_I2xy3z_S_b;
  Double I_NAI_Ix5y_Py_b = I_NAI_Kx6y_S_b+ABY*I_NAI_Ix5y_S_b;
  Double I_NAI_Ix4yz_Py_b = I_NAI_Kx5yz_S_b+ABY*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Py_b = I_NAI_Kx4y2z_S_b+ABY*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Py_b = I_NAI_Kx3y3z_S_b+ABY*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Py_b = I_NAI_Kx2y4z_S_b+ABY*I_NAI_Ixy4z_S_b;
  Double I_NAI_I6y_Py_b = I_NAI_K7y_S_b+ABY*I_NAI_I6y_S_b;
  Double I_NAI_I5yz_Py_b = I_NAI_K6yz_S_b+ABY*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Py_b = I_NAI_K5y2z_S_b+ABY*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Py_b = I_NAI_K4y3z_S_b+ABY*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Py_b = I_NAI_K3y4z_S_b+ABY*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Py_b = I_NAI_K2y5z_S_b+ABY*I_NAI_Iy5z_S_b;
  Double I_NAI_I5xz_Pz_b = I_NAI_K5x2z_S_b+ABZ*I_NAI_I5xz_S_b;
  Double I_NAI_I4xyz_Pz_b = I_NAI_K4xy2z_S_b+ABZ*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Pz_b = I_NAI_K4x3z_S_b+ABZ*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x2yz_Pz_b = I_NAI_K3x2y2z_S_b+ABZ*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Pz_b = I_NAI_K3xy3z_S_b+ABZ*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Pz_b = I_NAI_K3x4z_S_b+ABZ*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x3yz_Pz_b = I_NAI_K2x3y2z_S_b+ABZ*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Pz_b = I_NAI_K2x2y3z_S_b+ABZ*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Pz_b = I_NAI_K2xy4z_S_b+ABZ*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Pz_b = I_NAI_K2x5z_S_b+ABZ*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix4yz_Pz_b = I_NAI_Kx4y2z_S_b+ABZ*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Pz_b = I_NAI_Kx3y3z_S_b+ABZ*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Pz_b = I_NAI_Kx2y4z_S_b+ABZ*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Pz_b = I_NAI_Kxy5z_S_b+ABZ*I_NAI_Ixy4z_S_b;
  Double I_NAI_Ix5z_Pz_b = I_NAI_Kx6z_S_b+ABZ*I_NAI_Ix5z_S_b;
  Double I_NAI_I5yz_Pz_b = I_NAI_K5y2z_S_b+ABZ*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Pz_b = I_NAI_K4y3z_S_b+ABZ*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Pz_b = I_NAI_K3y4z_S_b+ABZ*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Pz_b = I_NAI_K2y5z_S_b+ABZ*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Pz_b = I_NAI_Ky6z_S_b+ABZ*I_NAI_Iy5z_S_b;
  Double I_NAI_I6z_Pz_b = I_NAI_K7z_S_b+ABZ*I_NAI_I6z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 48 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  Double I_NAI_H5x_D2x_b = I_NAI_I6x_Px_b+ABX*I_NAI_H5x_Px_b;
  Double I_NAI_H4xy_D2x_b = I_NAI_I5xy_Px_b+ABX*I_NAI_H4xy_Px_b;
  Double I_NAI_H4xz_D2x_b = I_NAI_I5xz_Px_b+ABX*I_NAI_H4xz_Px_b;
  Double I_NAI_H3x2y_D2x_b = I_NAI_I4x2y_Px_b+ABX*I_NAI_H3x2y_Px_b;
  Double I_NAI_H3xyz_D2x_b = I_NAI_I4xyz_Px_b+ABX*I_NAI_H3xyz_Px_b;
  Double I_NAI_H3x2z_D2x_b = I_NAI_I4x2z_Px_b+ABX*I_NAI_H3x2z_Px_b;
  Double I_NAI_H2x3y_D2x_b = I_NAI_I3x3y_Px_b+ABX*I_NAI_H2x3y_Px_b;
  Double I_NAI_H2x2yz_D2x_b = I_NAI_I3x2yz_Px_b+ABX*I_NAI_H2x2yz_Px_b;
  Double I_NAI_H2xy2z_D2x_b = I_NAI_I3xy2z_Px_b+ABX*I_NAI_H2xy2z_Px_b;
  Double I_NAI_H2x3z_D2x_b = I_NAI_I3x3z_Px_b+ABX*I_NAI_H2x3z_Px_b;
  Double I_NAI_Hx4y_D2x_b = I_NAI_I2x4y_Px_b+ABX*I_NAI_Hx4y_Px_b;
  Double I_NAI_Hx3yz_D2x_b = I_NAI_I2x3yz_Px_b+ABX*I_NAI_Hx3yz_Px_b;
  Double I_NAI_Hx2y2z_D2x_b = I_NAI_I2x2y2z_Px_b+ABX*I_NAI_Hx2y2z_Px_b;
  Double I_NAI_Hxy3z_D2x_b = I_NAI_I2xy3z_Px_b+ABX*I_NAI_Hxy3z_Px_b;
  Double I_NAI_Hx4z_D2x_b = I_NAI_I2x4z_Px_b+ABX*I_NAI_Hx4z_Px_b;
  Double I_NAI_H5y_D2x_b = I_NAI_Ix5y_Px_b+ABX*I_NAI_H5y_Px_b;
  Double I_NAI_H4yz_D2x_b = I_NAI_Ix4yz_Px_b+ABX*I_NAI_H4yz_Px_b;
  Double I_NAI_H3y2z_D2x_b = I_NAI_Ix3y2z_Px_b+ABX*I_NAI_H3y2z_Px_b;
  Double I_NAI_H2y3z_D2x_b = I_NAI_Ix2y3z_Px_b+ABX*I_NAI_H2y3z_Px_b;
  Double I_NAI_Hy4z_D2x_b = I_NAI_Ixy4z_Px_b+ABX*I_NAI_Hy4z_Px_b;
  Double I_NAI_H5z_D2x_b = I_NAI_Ix5z_Px_b+ABX*I_NAI_H5z_Px_b;
  Double I_NAI_H4xz_Dxy_b = I_NAI_I4xyz_Px_b+ABY*I_NAI_H4xz_Px_b;
  Double I_NAI_H3xyz_Dxy_b = I_NAI_I3x2yz_Px_b+ABY*I_NAI_H3xyz_Px_b;
  Double I_NAI_H3x2z_Dxy_b = I_NAI_I3xy2z_Px_b+ABY*I_NAI_H3x2z_Px_b;
  Double I_NAI_H2x2yz_Dxy_b = I_NAI_I2x3yz_Px_b+ABY*I_NAI_H2x2yz_Px_b;
  Double I_NAI_H2xy2z_Dxy_b = I_NAI_I2x2y2z_Px_b+ABY*I_NAI_H2xy2z_Px_b;
  Double I_NAI_H2x3z_Dxy_b = I_NAI_I2xy3z_Px_b+ABY*I_NAI_H2x3z_Px_b;
  Double I_NAI_Hx3yz_Dxy_b = I_NAI_Ix4yz_Px_b+ABY*I_NAI_Hx3yz_Px_b;
  Double I_NAI_Hx2y2z_Dxy_b = I_NAI_Ix3y2z_Px_b+ABY*I_NAI_Hx2y2z_Px_b;
  Double I_NAI_Hxy3z_Dxy_b = I_NAI_Ix2y3z_Px_b+ABY*I_NAI_Hxy3z_Px_b;
  Double I_NAI_Hx4z_Dxy_b = I_NAI_Ixy4z_Px_b+ABY*I_NAI_Hx4z_Px_b;
  Double I_NAI_H4yz_Dxy_b = I_NAI_I5yz_Px_b+ABY*I_NAI_H4yz_Px_b;
  Double I_NAI_H3y2z_Dxy_b = I_NAI_I4y2z_Px_b+ABY*I_NAI_H3y2z_Px_b;
  Double I_NAI_H2y3z_Dxy_b = I_NAI_I3y3z_Px_b+ABY*I_NAI_H2y3z_Px_b;
  Double I_NAI_Hy4z_Dxy_b = I_NAI_I2y4z_Px_b+ABY*I_NAI_Hy4z_Px_b;
  Double I_NAI_H5z_Dxy_b = I_NAI_Iy5z_Px_b+ABY*I_NAI_H5z_Px_b;
  Double I_NAI_H5x_D2y_b = I_NAI_I5xy_Py_b+ABY*I_NAI_H5x_Py_b;
  Double I_NAI_H4xy_D2y_b = I_NAI_I4x2y_Py_b+ABY*I_NAI_H4xy_Py_b;
  Double I_NAI_H4xz_D2y_b = I_NAI_I4xyz_Py_b+ABY*I_NAI_H4xz_Py_b;
  Double I_NAI_H3x2y_D2y_b = I_NAI_I3x3y_Py_b+ABY*I_NAI_H3x2y_Py_b;
  Double I_NAI_H3xyz_D2y_b = I_NAI_I3x2yz_Py_b+ABY*I_NAI_H3xyz_Py_b;
  Double I_NAI_H3x2z_D2y_b = I_NAI_I3xy2z_Py_b+ABY*I_NAI_H3x2z_Py_b;
  Double I_NAI_H2x3y_D2y_b = I_NAI_I2x4y_Py_b+ABY*I_NAI_H2x3y_Py_b;
  Double I_NAI_H2x2yz_D2y_b = I_NAI_I2x3yz_Py_b+ABY*I_NAI_H2x2yz_Py_b;
  Double I_NAI_H2xy2z_D2y_b = I_NAI_I2x2y2z_Py_b+ABY*I_NAI_H2xy2z_Py_b;
  Double I_NAI_H2x3z_D2y_b = I_NAI_I2xy3z_Py_b+ABY*I_NAI_H2x3z_Py_b;
  Double I_NAI_Hx4y_D2y_b = I_NAI_Ix5y_Py_b+ABY*I_NAI_Hx4y_Py_b;
  Double I_NAI_Hx3yz_D2y_b = I_NAI_Ix4yz_Py_b+ABY*I_NAI_Hx3yz_Py_b;
  Double I_NAI_Hx2y2z_D2y_b = I_NAI_Ix3y2z_Py_b+ABY*I_NAI_Hx2y2z_Py_b;
  Double I_NAI_Hxy3z_D2y_b = I_NAI_Ix2y3z_Py_b+ABY*I_NAI_Hxy3z_Py_b;
  Double I_NAI_Hx4z_D2y_b = I_NAI_Ixy4z_Py_b+ABY*I_NAI_Hx4z_Py_b;
  Double I_NAI_H5y_D2y_b = I_NAI_I6y_Py_b+ABY*I_NAI_H5y_Py_b;
  Double I_NAI_H4yz_D2y_b = I_NAI_I5yz_Py_b+ABY*I_NAI_H4yz_Py_b;
  Double I_NAI_H3y2z_D2y_b = I_NAI_I4y2z_Py_b+ABY*I_NAI_H3y2z_Py_b;
  Double I_NAI_H2y3z_D2y_b = I_NAI_I3y3z_Py_b+ABY*I_NAI_H2y3z_Py_b;
  Double I_NAI_Hy4z_D2y_b = I_NAI_I2y4z_Py_b+ABY*I_NAI_Hy4z_Py_b;
  Double I_NAI_H5z_D2y_b = I_NAI_Iy5z_Py_b+ABY*I_NAI_H5z_Py_b;
  Double I_NAI_H5x_D2z_b = I_NAI_I5xz_Pz_b+ABZ*I_NAI_H5x_Pz_b;
  Double I_NAI_H4xy_D2z_b = I_NAI_I4xyz_Pz_b+ABZ*I_NAI_H4xy_Pz_b;
  Double I_NAI_H4xz_D2z_b = I_NAI_I4x2z_Pz_b+ABZ*I_NAI_H4xz_Pz_b;
  Double I_NAI_H3x2y_D2z_b = I_NAI_I3x2yz_Pz_b+ABZ*I_NAI_H3x2y_Pz_b;
  Double I_NAI_H3xyz_D2z_b = I_NAI_I3xy2z_Pz_b+ABZ*I_NAI_H3xyz_Pz_b;
  Double I_NAI_H3x2z_D2z_b = I_NAI_I3x3z_Pz_b+ABZ*I_NAI_H3x2z_Pz_b;
  Double I_NAI_H2x3y_D2z_b = I_NAI_I2x3yz_Pz_b+ABZ*I_NAI_H2x3y_Pz_b;
  Double I_NAI_H2x2yz_D2z_b = I_NAI_I2x2y2z_Pz_b+ABZ*I_NAI_H2x2yz_Pz_b;
  Double I_NAI_H2xy2z_D2z_b = I_NAI_I2xy3z_Pz_b+ABZ*I_NAI_H2xy2z_Pz_b;
  Double I_NAI_H2x3z_D2z_b = I_NAI_I2x4z_Pz_b+ABZ*I_NAI_H2x3z_Pz_b;
  Double I_NAI_Hx4y_D2z_b = I_NAI_Ix4yz_Pz_b+ABZ*I_NAI_Hx4y_Pz_b;
  Double I_NAI_Hx3yz_D2z_b = I_NAI_Ix3y2z_Pz_b+ABZ*I_NAI_Hx3yz_Pz_b;
  Double I_NAI_Hx2y2z_D2z_b = I_NAI_Ix2y3z_Pz_b+ABZ*I_NAI_Hx2y2z_Pz_b;
  Double I_NAI_Hxy3z_D2z_b = I_NAI_Ixy4z_Pz_b+ABZ*I_NAI_Hxy3z_Pz_b;
  Double I_NAI_Hx4z_D2z_b = I_NAI_Ix5z_Pz_b+ABZ*I_NAI_Hx4z_Pz_b;
  Double I_NAI_H5y_D2z_b = I_NAI_I5yz_Pz_b+ABZ*I_NAI_H5y_Pz_b;
  Double I_NAI_H4yz_D2z_b = I_NAI_I4y2z_Pz_b+ABZ*I_NAI_H4yz_Pz_b;
  Double I_NAI_H3y2z_D2z_b = I_NAI_I3y3z_Pz_b+ABZ*I_NAI_H3y2z_Pz_b;
  Double I_NAI_H2y3z_D2z_b = I_NAI_I2y4z_Pz_b+ABZ*I_NAI_H2y3z_Pz_b;
  Double I_NAI_Hy4z_D2z_b = I_NAI_Iy5z_Pz_b+ABZ*I_NAI_Hy4z_Pz_b;
  Double I_NAI_H5z_D2z_b = I_NAI_I6z_Pz_b+ABZ*I_NAI_H5z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_b
   * RHS shell quartet name: SQ_NAI_G_D_b
   ************************************************************/
  Double I_NAI_G4x_F3x_b = I_NAI_H5x_D2x_b+ABX*I_NAI_G4x_D2x_b;
  Double I_NAI_G3xy_F3x_b = I_NAI_H4xy_D2x_b+ABX*I_NAI_G3xy_D2x_b;
  Double I_NAI_G3xz_F3x_b = I_NAI_H4xz_D2x_b+ABX*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2x2y_F3x_b = I_NAI_H3x2y_D2x_b+ABX*I_NAI_G2x2y_D2x_b;
  Double I_NAI_G2xyz_F3x_b = I_NAI_H3xyz_D2x_b+ABX*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F3x_b = I_NAI_H3x2z_D2x_b+ABX*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx3y_F3x_b = I_NAI_H2x3y_D2x_b+ABX*I_NAI_Gx3y_D2x_b;
  Double I_NAI_Gx2yz_F3x_b = I_NAI_H2x2yz_D2x_b+ABX*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F3x_b = I_NAI_H2xy2z_D2x_b+ABX*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F3x_b = I_NAI_H2x3z_D2x_b+ABX*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G4y_F3x_b = I_NAI_Hx4y_D2x_b+ABX*I_NAI_G4y_D2x_b;
  Double I_NAI_G3yz_F3x_b = I_NAI_Hx3yz_D2x_b+ABX*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F3x_b = I_NAI_Hx2y2z_D2x_b+ABX*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F3x_b = I_NAI_Hxy3z_D2x_b+ABX*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F3x_b = I_NAI_Hx4z_D2x_b+ABX*I_NAI_G4z_D2x_b;
  Double I_NAI_G4x_F2xy_b = I_NAI_H4xy_D2x_b+ABY*I_NAI_G4x_D2x_b;
  Double I_NAI_G3xy_F2xy_b = I_NAI_H3x2y_D2x_b+ABY*I_NAI_G3xy_D2x_b;
  Double I_NAI_G3xz_F2xy_b = I_NAI_H3xyz_D2x_b+ABY*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2x2y_F2xy_b = I_NAI_H2x3y_D2x_b+ABY*I_NAI_G2x2y_D2x_b;
  Double I_NAI_G2xyz_F2xy_b = I_NAI_H2x2yz_D2x_b+ABY*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F2xy_b = I_NAI_H2xy2z_D2x_b+ABY*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx3y_F2xy_b = I_NAI_Hx4y_D2x_b+ABY*I_NAI_Gx3y_D2x_b;
  Double I_NAI_Gx2yz_F2xy_b = I_NAI_Hx3yz_D2x_b+ABY*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F2xy_b = I_NAI_Hx2y2z_D2x_b+ABY*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F2xy_b = I_NAI_Hxy3z_D2x_b+ABY*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G4y_F2xy_b = I_NAI_H5y_D2x_b+ABY*I_NAI_G4y_D2x_b;
  Double I_NAI_G3yz_F2xy_b = I_NAI_H4yz_D2x_b+ABY*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F2xy_b = I_NAI_H3y2z_D2x_b+ABY*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F2xy_b = I_NAI_H2y3z_D2x_b+ABY*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F2xy_b = I_NAI_Hy4z_D2x_b+ABY*I_NAI_G4z_D2x_b;
  Double I_NAI_G4x_F2xz_b = I_NAI_H4xz_D2x_b+ABZ*I_NAI_G4x_D2x_b;
  Double I_NAI_G3xy_F2xz_b = I_NAI_H3xyz_D2x_b+ABZ*I_NAI_G3xy_D2x_b;
  Double I_NAI_G3xz_F2xz_b = I_NAI_H3x2z_D2x_b+ABZ*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2x2y_F2xz_b = I_NAI_H2x2yz_D2x_b+ABZ*I_NAI_G2x2y_D2x_b;
  Double I_NAI_G2xyz_F2xz_b = I_NAI_H2xy2z_D2x_b+ABZ*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F2xz_b = I_NAI_H2x3z_D2x_b+ABZ*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx3y_F2xz_b = I_NAI_Hx3yz_D2x_b+ABZ*I_NAI_Gx3y_D2x_b;
  Double I_NAI_Gx2yz_F2xz_b = I_NAI_Hx2y2z_D2x_b+ABZ*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F2xz_b = I_NAI_Hxy3z_D2x_b+ABZ*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F2xz_b = I_NAI_Hx4z_D2x_b+ABZ*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G4y_F2xz_b = I_NAI_H4yz_D2x_b+ABZ*I_NAI_G4y_D2x_b;
  Double I_NAI_G3yz_F2xz_b = I_NAI_H3y2z_D2x_b+ABZ*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F2xz_b = I_NAI_H2y3z_D2x_b+ABZ*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F2xz_b = I_NAI_Hy4z_D2x_b+ABZ*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F2xz_b = I_NAI_H5z_D2x_b+ABZ*I_NAI_G4z_D2x_b;
  Double I_NAI_G4x_Fx2y_b = I_NAI_H5x_D2y_b+ABX*I_NAI_G4x_D2y_b;
  Double I_NAI_G3xy_Fx2y_b = I_NAI_H4xy_D2y_b+ABX*I_NAI_G3xy_D2y_b;
  Double I_NAI_G3xz_Fx2y_b = I_NAI_H4xz_D2y_b+ABX*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2x2y_Fx2y_b = I_NAI_H3x2y_D2y_b+ABX*I_NAI_G2x2y_D2y_b;
  Double I_NAI_G2xyz_Fx2y_b = I_NAI_H3xyz_D2y_b+ABX*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_Fx2y_b = I_NAI_H3x2z_D2y_b+ABX*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx3y_Fx2y_b = I_NAI_H2x3y_D2y_b+ABX*I_NAI_Gx3y_D2y_b;
  Double I_NAI_Gx2yz_Fx2y_b = I_NAI_H2x2yz_D2y_b+ABX*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_Fx2y_b = I_NAI_H2xy2z_D2y_b+ABX*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_Fx2y_b = I_NAI_H2x3z_D2y_b+ABX*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G4y_Fx2y_b = I_NAI_Hx4y_D2y_b+ABX*I_NAI_G4y_D2y_b;
  Double I_NAI_G3yz_Fx2y_b = I_NAI_Hx3yz_D2y_b+ABX*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_Fx2y_b = I_NAI_Hx2y2z_D2y_b+ABX*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_Fx2y_b = I_NAI_Hxy3z_D2y_b+ABX*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_Fx2y_b = I_NAI_Hx4z_D2y_b+ABX*I_NAI_G4z_D2y_b;
  Double I_NAI_G4x_Fxyz_b = I_NAI_H4xz_Dxy_b+ABZ*I_NAI_G4x_Dxy_b;
  Double I_NAI_G3xy_Fxyz_b = I_NAI_H3xyz_Dxy_b+ABZ*I_NAI_G3xy_Dxy_b;
  Double I_NAI_G3xz_Fxyz_b = I_NAI_H3x2z_Dxy_b+ABZ*I_NAI_G3xz_Dxy_b;
  Double I_NAI_G2x2y_Fxyz_b = I_NAI_H2x2yz_Dxy_b+ABZ*I_NAI_G2x2y_Dxy_b;
  Double I_NAI_G2xyz_Fxyz_b = I_NAI_H2xy2z_Dxy_b+ABZ*I_NAI_G2xyz_Dxy_b;
  Double I_NAI_G2x2z_Fxyz_b = I_NAI_H2x3z_Dxy_b+ABZ*I_NAI_G2x2z_Dxy_b;
  Double I_NAI_Gx3y_Fxyz_b = I_NAI_Hx3yz_Dxy_b+ABZ*I_NAI_Gx3y_Dxy_b;
  Double I_NAI_Gx2yz_Fxyz_b = I_NAI_Hx2y2z_Dxy_b+ABZ*I_NAI_Gx2yz_Dxy_b;
  Double I_NAI_Gxy2z_Fxyz_b = I_NAI_Hxy3z_Dxy_b+ABZ*I_NAI_Gxy2z_Dxy_b;
  Double I_NAI_Gx3z_Fxyz_b = I_NAI_Hx4z_Dxy_b+ABZ*I_NAI_Gx3z_Dxy_b;
  Double I_NAI_G4y_Fxyz_b = I_NAI_H4yz_Dxy_b+ABZ*I_NAI_G4y_Dxy_b;
  Double I_NAI_G3yz_Fxyz_b = I_NAI_H3y2z_Dxy_b+ABZ*I_NAI_G3yz_Dxy_b;
  Double I_NAI_G2y2z_Fxyz_b = I_NAI_H2y3z_Dxy_b+ABZ*I_NAI_G2y2z_Dxy_b;
  Double I_NAI_Gy3z_Fxyz_b = I_NAI_Hy4z_Dxy_b+ABZ*I_NAI_Gy3z_Dxy_b;
  Double I_NAI_G4z_Fxyz_b = I_NAI_H5z_Dxy_b+ABZ*I_NAI_G4z_Dxy_b;
  Double I_NAI_G4x_Fx2z_b = I_NAI_H5x_D2z_b+ABX*I_NAI_G4x_D2z_b;
  Double I_NAI_G3xy_Fx2z_b = I_NAI_H4xy_D2z_b+ABX*I_NAI_G3xy_D2z_b;
  Double I_NAI_G3xz_Fx2z_b = I_NAI_H4xz_D2z_b+ABX*I_NAI_G3xz_D2z_b;
  Double I_NAI_G2x2y_Fx2z_b = I_NAI_H3x2y_D2z_b+ABX*I_NAI_G2x2y_D2z_b;
  Double I_NAI_G2xyz_Fx2z_b = I_NAI_H3xyz_D2z_b+ABX*I_NAI_G2xyz_D2z_b;
  Double I_NAI_G2x2z_Fx2z_b = I_NAI_H3x2z_D2z_b+ABX*I_NAI_G2x2z_D2z_b;
  Double I_NAI_Gx3y_Fx2z_b = I_NAI_H2x3y_D2z_b+ABX*I_NAI_Gx3y_D2z_b;
  Double I_NAI_Gx2yz_Fx2z_b = I_NAI_H2x2yz_D2z_b+ABX*I_NAI_Gx2yz_D2z_b;
  Double I_NAI_Gxy2z_Fx2z_b = I_NAI_H2xy2z_D2z_b+ABX*I_NAI_Gxy2z_D2z_b;
  Double I_NAI_Gx3z_Fx2z_b = I_NAI_H2x3z_D2z_b+ABX*I_NAI_Gx3z_D2z_b;
  Double I_NAI_G4y_Fx2z_b = I_NAI_Hx4y_D2z_b+ABX*I_NAI_G4y_D2z_b;
  Double I_NAI_G3yz_Fx2z_b = I_NAI_Hx3yz_D2z_b+ABX*I_NAI_G3yz_D2z_b;
  Double I_NAI_G2y2z_Fx2z_b = I_NAI_Hx2y2z_D2z_b+ABX*I_NAI_G2y2z_D2z_b;
  Double I_NAI_Gy3z_Fx2z_b = I_NAI_Hxy3z_D2z_b+ABX*I_NAI_Gy3z_D2z_b;
  Double I_NAI_G4z_Fx2z_b = I_NAI_Hx4z_D2z_b+ABX*I_NAI_G4z_D2z_b;
  Double I_NAI_G4x_F3y_b = I_NAI_H4xy_D2y_b+ABY*I_NAI_G4x_D2y_b;
  Double I_NAI_G3xy_F3y_b = I_NAI_H3x2y_D2y_b+ABY*I_NAI_G3xy_D2y_b;
  Double I_NAI_G3xz_F3y_b = I_NAI_H3xyz_D2y_b+ABY*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2x2y_F3y_b = I_NAI_H2x3y_D2y_b+ABY*I_NAI_G2x2y_D2y_b;
  Double I_NAI_G2xyz_F3y_b = I_NAI_H2x2yz_D2y_b+ABY*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_F3y_b = I_NAI_H2xy2z_D2y_b+ABY*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx3y_F3y_b = I_NAI_Hx4y_D2y_b+ABY*I_NAI_Gx3y_D2y_b;
  Double I_NAI_Gx2yz_F3y_b = I_NAI_Hx3yz_D2y_b+ABY*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_F3y_b = I_NAI_Hx2y2z_D2y_b+ABY*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_F3y_b = I_NAI_Hxy3z_D2y_b+ABY*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G4y_F3y_b = I_NAI_H5y_D2y_b+ABY*I_NAI_G4y_D2y_b;
  Double I_NAI_G3yz_F3y_b = I_NAI_H4yz_D2y_b+ABY*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_F3y_b = I_NAI_H3y2z_D2y_b+ABY*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_F3y_b = I_NAI_H2y3z_D2y_b+ABY*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_F3y_b = I_NAI_Hy4z_D2y_b+ABY*I_NAI_G4z_D2y_b;
  Double I_NAI_G4x_F2yz_b = I_NAI_H4xz_D2y_b+ABZ*I_NAI_G4x_D2y_b;
  Double I_NAI_G3xy_F2yz_b = I_NAI_H3xyz_D2y_b+ABZ*I_NAI_G3xy_D2y_b;
  Double I_NAI_G3xz_F2yz_b = I_NAI_H3x2z_D2y_b+ABZ*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2x2y_F2yz_b = I_NAI_H2x2yz_D2y_b+ABZ*I_NAI_G2x2y_D2y_b;
  Double I_NAI_G2xyz_F2yz_b = I_NAI_H2xy2z_D2y_b+ABZ*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_F2yz_b = I_NAI_H2x3z_D2y_b+ABZ*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx3y_F2yz_b = I_NAI_Hx3yz_D2y_b+ABZ*I_NAI_Gx3y_D2y_b;
  Double I_NAI_Gx2yz_F2yz_b = I_NAI_Hx2y2z_D2y_b+ABZ*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_F2yz_b = I_NAI_Hxy3z_D2y_b+ABZ*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_F2yz_b = I_NAI_Hx4z_D2y_b+ABZ*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G4y_F2yz_b = I_NAI_H4yz_D2y_b+ABZ*I_NAI_G4y_D2y_b;
  Double I_NAI_G3yz_F2yz_b = I_NAI_H3y2z_D2y_b+ABZ*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_F2yz_b = I_NAI_H2y3z_D2y_b+ABZ*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_F2yz_b = I_NAI_Hy4z_D2y_b+ABZ*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_F2yz_b = I_NAI_H5z_D2y_b+ABZ*I_NAI_G4z_D2y_b;
  Double I_NAI_G4x_Fy2z_b = I_NAI_H4xy_D2z_b+ABY*I_NAI_G4x_D2z_b;
  Double I_NAI_G3xy_Fy2z_b = I_NAI_H3x2y_D2z_b+ABY*I_NAI_G3xy_D2z_b;
  Double I_NAI_G3xz_Fy2z_b = I_NAI_H3xyz_D2z_b+ABY*I_NAI_G3xz_D2z_b;
  Double I_NAI_G2x2y_Fy2z_b = I_NAI_H2x3y_D2z_b+ABY*I_NAI_G2x2y_D2z_b;
  Double I_NAI_G2xyz_Fy2z_b = I_NAI_H2x2yz_D2z_b+ABY*I_NAI_G2xyz_D2z_b;
  Double I_NAI_G2x2z_Fy2z_b = I_NAI_H2xy2z_D2z_b+ABY*I_NAI_G2x2z_D2z_b;
  Double I_NAI_Gx3y_Fy2z_b = I_NAI_Hx4y_D2z_b+ABY*I_NAI_Gx3y_D2z_b;
  Double I_NAI_Gx2yz_Fy2z_b = I_NAI_Hx3yz_D2z_b+ABY*I_NAI_Gx2yz_D2z_b;
  Double I_NAI_Gxy2z_Fy2z_b = I_NAI_Hx2y2z_D2z_b+ABY*I_NAI_Gxy2z_D2z_b;
  Double I_NAI_Gx3z_Fy2z_b = I_NAI_Hxy3z_D2z_b+ABY*I_NAI_Gx3z_D2z_b;
  Double I_NAI_G4y_Fy2z_b = I_NAI_H5y_D2z_b+ABY*I_NAI_G4y_D2z_b;
  Double I_NAI_G3yz_Fy2z_b = I_NAI_H4yz_D2z_b+ABY*I_NAI_G3yz_D2z_b;
  Double I_NAI_G2y2z_Fy2z_b = I_NAI_H3y2z_D2z_b+ABY*I_NAI_G2y2z_D2z_b;
  Double I_NAI_Gy3z_Fy2z_b = I_NAI_H2y3z_D2z_b+ABY*I_NAI_Gy3z_D2z_b;
  Double I_NAI_G4z_Fy2z_b = I_NAI_Hy4z_D2z_b+ABY*I_NAI_G4z_D2z_b;
  Double I_NAI_G4x_F3z_b = I_NAI_H4xz_D2z_b+ABZ*I_NAI_G4x_D2z_b;
  Double I_NAI_G3xy_F3z_b = I_NAI_H3xyz_D2z_b+ABZ*I_NAI_G3xy_D2z_b;
  Double I_NAI_G3xz_F3z_b = I_NAI_H3x2z_D2z_b+ABZ*I_NAI_G3xz_D2z_b;
  Double I_NAI_G2x2y_F3z_b = I_NAI_H2x2yz_D2z_b+ABZ*I_NAI_G2x2y_D2z_b;
  Double I_NAI_G2xyz_F3z_b = I_NAI_H2xy2z_D2z_b+ABZ*I_NAI_G2xyz_D2z_b;
  Double I_NAI_G2x2z_F3z_b = I_NAI_H2x3z_D2z_b+ABZ*I_NAI_G2x2z_D2z_b;
  Double I_NAI_Gx3y_F3z_b = I_NAI_Hx3yz_D2z_b+ABZ*I_NAI_Gx3y_D2z_b;
  Double I_NAI_Gx2yz_F3z_b = I_NAI_Hx2y2z_D2z_b+ABZ*I_NAI_Gx2yz_D2z_b;
  Double I_NAI_Gxy2z_F3z_b = I_NAI_Hxy3z_D2z_b+ABZ*I_NAI_Gxy2z_D2z_b;
  Double I_NAI_Gx3z_F3z_b = I_NAI_Hx4z_D2z_b+ABZ*I_NAI_Gx3z_D2z_b;
  Double I_NAI_G4y_F3z_b = I_NAI_H4yz_D2z_b+ABZ*I_NAI_G4y_D2z_b;
  Double I_NAI_G3yz_F3z_b = I_NAI_H3y2z_D2z_b+ABZ*I_NAI_G3yz_D2z_b;
  Double I_NAI_G2y2z_F3z_b = I_NAI_H2y3z_D2z_b+ABZ*I_NAI_G2y2z_D2z_b;
  Double I_NAI_Gy3z_F3z_b = I_NAI_Hy4z_D2z_b+ABZ*I_NAI_Gy3z_D2z_b;
  Double I_NAI_G4z_F3z_b = I_NAI_H5z_D2z_b+ABZ*I_NAI_G4z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_G_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_F_F_b
   ************************************************************/
  Double I_NAI_F3x_G4x_b = I_NAI_G4x_F3x_b+ABX*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G4x_b = I_NAI_G3xy_F3x_b+ABX*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G4x_b = I_NAI_G3xz_F3x_b+ABX*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G4x_b = I_NAI_G2x2y_F3x_b+ABX*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G4x_b = I_NAI_G2xyz_F3x_b+ABX*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G4x_b = I_NAI_G2x2z_F3x_b+ABX*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G4x_b = I_NAI_Gx3y_F3x_b+ABX*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G4x_b = I_NAI_Gx2yz_F3x_b+ABX*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G4x_b = I_NAI_Gxy2z_F3x_b+ABX*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G4x_b = I_NAI_Gx3z_F3x_b+ABX*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G3xy_b = I_NAI_G3xy_F3x_b+ABY*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G3xy_b = I_NAI_G2x2y_F3x_b+ABY*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G3xy_b = I_NAI_G2xyz_F3x_b+ABY*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G3xy_b = I_NAI_Gx3y_F3x_b+ABY*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G3xy_b = I_NAI_Gx2yz_F3x_b+ABY*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G3xy_b = I_NAI_Gxy2z_F3x_b+ABY*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G3xy_b = I_NAI_G4y_F3x_b+ABY*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G3xy_b = I_NAI_G3yz_F3x_b+ABY*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G3xy_b = I_NAI_G2y2z_F3x_b+ABY*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G3xy_b = I_NAI_Gy3z_F3x_b+ABY*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G3xz_b = I_NAI_G3xz_F3x_b+ABZ*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G3xz_b = I_NAI_G2xyz_F3x_b+ABZ*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G3xz_b = I_NAI_G2x2z_F3x_b+ABZ*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G3xz_b = I_NAI_Gx2yz_F3x_b+ABZ*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G3xz_b = I_NAI_Gxy2z_F3x_b+ABZ*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G3xz_b = I_NAI_Gx3z_F3x_b+ABZ*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G3xz_b = I_NAI_G3yz_F3x_b+ABZ*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G3xz_b = I_NAI_G2y2z_F3x_b+ABZ*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G3xz_b = I_NAI_Gy3z_F3x_b+ABZ*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G3xz_b = I_NAI_G4z_F3x_b+ABZ*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G2x2y_b = I_NAI_G3xy_F2xy_b+ABY*I_NAI_F3x_F2xy_b;
  Double I_NAI_F2xy_G2x2y_b = I_NAI_G2x2y_F2xy_b+ABY*I_NAI_F2xy_F2xy_b;
  Double I_NAI_F2xz_G2x2y_b = I_NAI_G2xyz_F2xy_b+ABY*I_NAI_F2xz_F2xy_b;
  Double I_NAI_Fx2y_G2x2y_b = I_NAI_Gx3y_F2xy_b+ABY*I_NAI_Fx2y_F2xy_b;
  Double I_NAI_Fxyz_G2x2y_b = I_NAI_Gx2yz_F2xy_b+ABY*I_NAI_Fxyz_F2xy_b;
  Double I_NAI_Fx2z_G2x2y_b = I_NAI_Gxy2z_F2xy_b+ABY*I_NAI_Fx2z_F2xy_b;
  Double I_NAI_F3y_G2x2y_b = I_NAI_G4y_F2xy_b+ABY*I_NAI_F3y_F2xy_b;
  Double I_NAI_F2yz_G2x2y_b = I_NAI_G3yz_F2xy_b+ABY*I_NAI_F2yz_F2xy_b;
  Double I_NAI_Fy2z_G2x2y_b = I_NAI_G2y2z_F2xy_b+ABY*I_NAI_Fy2z_F2xy_b;
  Double I_NAI_F3z_G2x2y_b = I_NAI_Gy3z_F2xy_b+ABY*I_NAI_F3z_F2xy_b;
  Double I_NAI_F3x_G2xyz_b = I_NAI_G3xz_F2xy_b+ABZ*I_NAI_F3x_F2xy_b;
  Double I_NAI_F2xy_G2xyz_b = I_NAI_G2xyz_F2xy_b+ABZ*I_NAI_F2xy_F2xy_b;
  Double I_NAI_F2xz_G2xyz_b = I_NAI_G2x2z_F2xy_b+ABZ*I_NAI_F2xz_F2xy_b;
  Double I_NAI_Fx2y_G2xyz_b = I_NAI_Gx2yz_F2xy_b+ABZ*I_NAI_Fx2y_F2xy_b;
  Double I_NAI_Fxyz_G2xyz_b = I_NAI_Gxy2z_F2xy_b+ABZ*I_NAI_Fxyz_F2xy_b;
  Double I_NAI_Fx2z_G2xyz_b = I_NAI_Gx3z_F2xy_b+ABZ*I_NAI_Fx2z_F2xy_b;
  Double I_NAI_F3y_G2xyz_b = I_NAI_G3yz_F2xy_b+ABZ*I_NAI_F3y_F2xy_b;
  Double I_NAI_F2yz_G2xyz_b = I_NAI_G2y2z_F2xy_b+ABZ*I_NAI_F2yz_F2xy_b;
  Double I_NAI_Fy2z_G2xyz_b = I_NAI_Gy3z_F2xy_b+ABZ*I_NAI_Fy2z_F2xy_b;
  Double I_NAI_F3z_G2xyz_b = I_NAI_G4z_F2xy_b+ABZ*I_NAI_F3z_F2xy_b;
  Double I_NAI_F3x_G2x2z_b = I_NAI_G3xz_F2xz_b+ABZ*I_NAI_F3x_F2xz_b;
  Double I_NAI_F2xy_G2x2z_b = I_NAI_G2xyz_F2xz_b+ABZ*I_NAI_F2xy_F2xz_b;
  Double I_NAI_F2xz_G2x2z_b = I_NAI_G2x2z_F2xz_b+ABZ*I_NAI_F2xz_F2xz_b;
  Double I_NAI_Fx2y_G2x2z_b = I_NAI_Gx2yz_F2xz_b+ABZ*I_NAI_Fx2y_F2xz_b;
  Double I_NAI_Fxyz_G2x2z_b = I_NAI_Gxy2z_F2xz_b+ABZ*I_NAI_Fxyz_F2xz_b;
  Double I_NAI_Fx2z_G2x2z_b = I_NAI_Gx3z_F2xz_b+ABZ*I_NAI_Fx2z_F2xz_b;
  Double I_NAI_F3y_G2x2z_b = I_NAI_G3yz_F2xz_b+ABZ*I_NAI_F3y_F2xz_b;
  Double I_NAI_F2yz_G2x2z_b = I_NAI_G2y2z_F2xz_b+ABZ*I_NAI_F2yz_F2xz_b;
  Double I_NAI_Fy2z_G2x2z_b = I_NAI_Gy3z_F2xz_b+ABZ*I_NAI_Fy2z_F2xz_b;
  Double I_NAI_F3z_G2x2z_b = I_NAI_G4z_F2xz_b+ABZ*I_NAI_F3z_F2xz_b;
  Double I_NAI_F3x_Gx3y_b = I_NAI_G4x_F3y_b+ABX*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_Gx3y_b = I_NAI_G3xy_F3y_b+ABX*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_Gx3y_b = I_NAI_G3xz_F3y_b+ABX*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_Gx3y_b = I_NAI_G2x2y_F3y_b+ABX*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_Gx3y_b = I_NAI_G2xyz_F3y_b+ABX*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_Gx3y_b = I_NAI_G2x2z_F3y_b+ABX*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_Gx3y_b = I_NAI_Gx3y_F3y_b+ABX*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_Gx3y_b = I_NAI_Gx2yz_F3y_b+ABX*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_Gx3y_b = I_NAI_Gxy2z_F3y_b+ABX*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_Gx3y_b = I_NAI_Gx3z_F3y_b+ABX*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_Gx2yz_b = I_NAI_G3xz_Fx2y_b+ABZ*I_NAI_F3x_Fx2y_b;
  Double I_NAI_F2xy_Gx2yz_b = I_NAI_G2xyz_Fx2y_b+ABZ*I_NAI_F2xy_Fx2y_b;
  Double I_NAI_F2xz_Gx2yz_b = I_NAI_G2x2z_Fx2y_b+ABZ*I_NAI_F2xz_Fx2y_b;
  Double I_NAI_Fx2y_Gx2yz_b = I_NAI_Gx2yz_Fx2y_b+ABZ*I_NAI_Fx2y_Fx2y_b;
  Double I_NAI_Fxyz_Gx2yz_b = I_NAI_Gxy2z_Fx2y_b+ABZ*I_NAI_Fxyz_Fx2y_b;
  Double I_NAI_Fx2z_Gx2yz_b = I_NAI_Gx3z_Fx2y_b+ABZ*I_NAI_Fx2z_Fx2y_b;
  Double I_NAI_F3y_Gx2yz_b = I_NAI_G3yz_Fx2y_b+ABZ*I_NAI_F3y_Fx2y_b;
  Double I_NAI_F2yz_Gx2yz_b = I_NAI_G2y2z_Fx2y_b+ABZ*I_NAI_F2yz_Fx2y_b;
  Double I_NAI_Fy2z_Gx2yz_b = I_NAI_Gy3z_Fx2y_b+ABZ*I_NAI_Fy2z_Fx2y_b;
  Double I_NAI_F3z_Gx2yz_b = I_NAI_G4z_Fx2y_b+ABZ*I_NAI_F3z_Fx2y_b;
  Double I_NAI_F3x_Gxy2z_b = I_NAI_G3xy_Fx2z_b+ABY*I_NAI_F3x_Fx2z_b;
  Double I_NAI_F2xy_Gxy2z_b = I_NAI_G2x2y_Fx2z_b+ABY*I_NAI_F2xy_Fx2z_b;
  Double I_NAI_F2xz_Gxy2z_b = I_NAI_G2xyz_Fx2z_b+ABY*I_NAI_F2xz_Fx2z_b;
  Double I_NAI_Fx2y_Gxy2z_b = I_NAI_Gx3y_Fx2z_b+ABY*I_NAI_Fx2y_Fx2z_b;
  Double I_NAI_Fxyz_Gxy2z_b = I_NAI_Gx2yz_Fx2z_b+ABY*I_NAI_Fxyz_Fx2z_b;
  Double I_NAI_Fx2z_Gxy2z_b = I_NAI_Gxy2z_Fx2z_b+ABY*I_NAI_Fx2z_Fx2z_b;
  Double I_NAI_F3y_Gxy2z_b = I_NAI_G4y_Fx2z_b+ABY*I_NAI_F3y_Fx2z_b;
  Double I_NAI_F2yz_Gxy2z_b = I_NAI_G3yz_Fx2z_b+ABY*I_NAI_F2yz_Fx2z_b;
  Double I_NAI_Fy2z_Gxy2z_b = I_NAI_G2y2z_Fx2z_b+ABY*I_NAI_Fy2z_Fx2z_b;
  Double I_NAI_F3z_Gxy2z_b = I_NAI_Gy3z_Fx2z_b+ABY*I_NAI_F3z_Fx2z_b;
  Double I_NAI_F3x_Gx3z_b = I_NAI_G4x_F3z_b+ABX*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_Gx3z_b = I_NAI_G3xy_F3z_b+ABX*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_Gx3z_b = I_NAI_G3xz_F3z_b+ABX*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_Gx3z_b = I_NAI_G2x2y_F3z_b+ABX*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_Gx3z_b = I_NAI_G2xyz_F3z_b+ABX*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_Gx3z_b = I_NAI_G2x2z_F3z_b+ABX*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_Gx3z_b = I_NAI_Gx3y_F3z_b+ABX*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_Gx3z_b = I_NAI_Gx2yz_F3z_b+ABX*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_Gx3z_b = I_NAI_Gxy2z_F3z_b+ABX*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_Gx3z_b = I_NAI_Gx3z_F3z_b+ABX*I_NAI_F3z_F3z_b;
  Double I_NAI_F3x_G4y_b = I_NAI_G3xy_F3y_b+ABY*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_G4y_b = I_NAI_G2x2y_F3y_b+ABY*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_G4y_b = I_NAI_G2xyz_F3y_b+ABY*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_G4y_b = I_NAI_Gx3y_F3y_b+ABY*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_G4y_b = I_NAI_Gx2yz_F3y_b+ABY*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_G4y_b = I_NAI_Gxy2z_F3y_b+ABY*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_G4y_b = I_NAI_G4y_F3y_b+ABY*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_G4y_b = I_NAI_G3yz_F3y_b+ABY*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_G4y_b = I_NAI_G2y2z_F3y_b+ABY*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_G4y_b = I_NAI_Gy3z_F3y_b+ABY*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_G3yz_b = I_NAI_G3xz_F3y_b+ABZ*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_G3yz_b = I_NAI_G2xyz_F3y_b+ABZ*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_G3yz_b = I_NAI_G2x2z_F3y_b+ABZ*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_G3yz_b = I_NAI_Gx2yz_F3y_b+ABZ*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_G3yz_b = I_NAI_Gxy2z_F3y_b+ABZ*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_G3yz_b = I_NAI_Gx3z_F3y_b+ABZ*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_G3yz_b = I_NAI_G3yz_F3y_b+ABZ*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_G3yz_b = I_NAI_G2y2z_F3y_b+ABZ*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_G3yz_b = I_NAI_Gy3z_F3y_b+ABZ*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_G3yz_b = I_NAI_G4z_F3y_b+ABZ*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_G2y2z_b = I_NAI_G3xz_F2yz_b+ABZ*I_NAI_F3x_F2yz_b;
  Double I_NAI_F2xy_G2y2z_b = I_NAI_G2xyz_F2yz_b+ABZ*I_NAI_F2xy_F2yz_b;
  Double I_NAI_F2xz_G2y2z_b = I_NAI_G2x2z_F2yz_b+ABZ*I_NAI_F2xz_F2yz_b;
  Double I_NAI_Fx2y_G2y2z_b = I_NAI_Gx2yz_F2yz_b+ABZ*I_NAI_Fx2y_F2yz_b;
  Double I_NAI_Fxyz_G2y2z_b = I_NAI_Gxy2z_F2yz_b+ABZ*I_NAI_Fxyz_F2yz_b;
  Double I_NAI_Fx2z_G2y2z_b = I_NAI_Gx3z_F2yz_b+ABZ*I_NAI_Fx2z_F2yz_b;
  Double I_NAI_F3y_G2y2z_b = I_NAI_G3yz_F2yz_b+ABZ*I_NAI_F3y_F2yz_b;
  Double I_NAI_F2yz_G2y2z_b = I_NAI_G2y2z_F2yz_b+ABZ*I_NAI_F2yz_F2yz_b;
  Double I_NAI_Fy2z_G2y2z_b = I_NAI_Gy3z_F2yz_b+ABZ*I_NAI_Fy2z_F2yz_b;
  Double I_NAI_F3z_G2y2z_b = I_NAI_G4z_F2yz_b+ABZ*I_NAI_F3z_F2yz_b;
  Double I_NAI_F3x_Gy3z_b = I_NAI_G3xy_F3z_b+ABY*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_Gy3z_b = I_NAI_G2x2y_F3z_b+ABY*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_Gy3z_b = I_NAI_G2xyz_F3z_b+ABY*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_Gy3z_b = I_NAI_Gx3y_F3z_b+ABY*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_Gy3z_b = I_NAI_Gx2yz_F3z_b+ABY*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_Gy3z_b = I_NAI_Gxy2z_F3z_b+ABY*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_Gy3z_b = I_NAI_G4y_F3z_b+ABY*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_Gy3z_b = I_NAI_G3yz_F3z_b+ABY*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_Gy3z_b = I_NAI_G2y2z_F3z_b+ABY*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_Gy3z_b = I_NAI_Gy3z_F3z_b+ABY*I_NAI_F3z_F3z_b;
  Double I_NAI_F3x_G4z_b = I_NAI_G3xz_F3z_b+ABZ*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_G4z_b = I_NAI_G2xyz_F3z_b+ABZ*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_G4z_b = I_NAI_G2x2z_F3z_b+ABZ*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_G4z_b = I_NAI_Gx2yz_F3z_b+ABZ*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_G4z_b = I_NAI_Gxy2z_F3z_b+ABZ*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_G4z_b = I_NAI_Gx3z_F3z_b+ABZ*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_G4z_b = I_NAI_G3yz_F3z_b+ABZ*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_G4z_b = I_NAI_G2y2z_F3z_b+ABZ*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_G4z_b = I_NAI_Gy3z_F3z_b+ABZ*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_G4z_b = I_NAI_G4z_F3z_b+ABZ*I_NAI_F3z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_aa
   * RHS shell quartet name: SQ_NAI_I_S_aa
   ************************************************************/
  Double I_NAI_I6x_Px_aa = I_NAI_K7x_S_aa+ABX*I_NAI_I6x_S_aa;
  Double I_NAI_I5xy_Px_aa = I_NAI_K6xy_S_aa+ABX*I_NAI_I5xy_S_aa;
  Double I_NAI_I5xz_Px_aa = I_NAI_K6xz_S_aa+ABX*I_NAI_I5xz_S_aa;
  Double I_NAI_I4x2y_Px_aa = I_NAI_K5x2y_S_aa+ABX*I_NAI_I4x2y_S_aa;
  Double I_NAI_I4xyz_Px_aa = I_NAI_K5xyz_S_aa+ABX*I_NAI_I4xyz_S_aa;
  Double I_NAI_I4x2z_Px_aa = I_NAI_K5x2z_S_aa+ABX*I_NAI_I4x2z_S_aa;
  Double I_NAI_I3x3y_Px_aa = I_NAI_K4x3y_S_aa+ABX*I_NAI_I3x3y_S_aa;
  Double I_NAI_I3x2yz_Px_aa = I_NAI_K4x2yz_S_aa+ABX*I_NAI_I3x2yz_S_aa;
  Double I_NAI_I3xy2z_Px_aa = I_NAI_K4xy2z_S_aa+ABX*I_NAI_I3xy2z_S_aa;
  Double I_NAI_I3x3z_Px_aa = I_NAI_K4x3z_S_aa+ABX*I_NAI_I3x3z_S_aa;
  Double I_NAI_I2x4y_Px_aa = I_NAI_K3x4y_S_aa+ABX*I_NAI_I2x4y_S_aa;
  Double I_NAI_I2x3yz_Px_aa = I_NAI_K3x3yz_S_aa+ABX*I_NAI_I2x3yz_S_aa;
  Double I_NAI_I2x2y2z_Px_aa = I_NAI_K3x2y2z_S_aa+ABX*I_NAI_I2x2y2z_S_aa;
  Double I_NAI_I2xy3z_Px_aa = I_NAI_K3xy3z_S_aa+ABX*I_NAI_I2xy3z_S_aa;
  Double I_NAI_I2x4z_Px_aa = I_NAI_K3x4z_S_aa+ABX*I_NAI_I2x4z_S_aa;
  Double I_NAI_Ix5y_Px_aa = I_NAI_K2x5y_S_aa+ABX*I_NAI_Ix5y_S_aa;
  Double I_NAI_Ix4yz_Px_aa = I_NAI_K2x4yz_S_aa+ABX*I_NAI_Ix4yz_S_aa;
  Double I_NAI_Ix3y2z_Px_aa = I_NAI_K2x3y2z_S_aa+ABX*I_NAI_Ix3y2z_S_aa;
  Double I_NAI_Ix2y3z_Px_aa = I_NAI_K2x2y3z_S_aa+ABX*I_NAI_Ix2y3z_S_aa;
  Double I_NAI_Ixy4z_Px_aa = I_NAI_K2xy4z_S_aa+ABX*I_NAI_Ixy4z_S_aa;
  Double I_NAI_Ix5z_Px_aa = I_NAI_K2x5z_S_aa+ABX*I_NAI_Ix5z_S_aa;
  Double I_NAI_I6y_Px_aa = I_NAI_Kx6y_S_aa+ABX*I_NAI_I6y_S_aa;
  Double I_NAI_I5yz_Px_aa = I_NAI_Kx5yz_S_aa+ABX*I_NAI_I5yz_S_aa;
  Double I_NAI_I4y2z_Px_aa = I_NAI_Kx4y2z_S_aa+ABX*I_NAI_I4y2z_S_aa;
  Double I_NAI_I3y3z_Px_aa = I_NAI_Kx3y3z_S_aa+ABX*I_NAI_I3y3z_S_aa;
  Double I_NAI_I2y4z_Px_aa = I_NAI_Kx2y4z_S_aa+ABX*I_NAI_I2y4z_S_aa;
  Double I_NAI_Iy5z_Px_aa = I_NAI_Kxy5z_S_aa+ABX*I_NAI_Iy5z_S_aa;
  Double I_NAI_I6z_Px_aa = I_NAI_Kx6z_S_aa+ABX*I_NAI_I6z_S_aa;
  Double I_NAI_I6x_Py_aa = I_NAI_K6xy_S_aa+ABY*I_NAI_I6x_S_aa;
  Double I_NAI_I5xy_Py_aa = I_NAI_K5x2y_S_aa+ABY*I_NAI_I5xy_S_aa;
  Double I_NAI_I5xz_Py_aa = I_NAI_K5xyz_S_aa+ABY*I_NAI_I5xz_S_aa;
  Double I_NAI_I4x2y_Py_aa = I_NAI_K4x3y_S_aa+ABY*I_NAI_I4x2y_S_aa;
  Double I_NAI_I4xyz_Py_aa = I_NAI_K4x2yz_S_aa+ABY*I_NAI_I4xyz_S_aa;
  Double I_NAI_I4x2z_Py_aa = I_NAI_K4xy2z_S_aa+ABY*I_NAI_I4x2z_S_aa;
  Double I_NAI_I3x3y_Py_aa = I_NAI_K3x4y_S_aa+ABY*I_NAI_I3x3y_S_aa;
  Double I_NAI_I3x2yz_Py_aa = I_NAI_K3x3yz_S_aa+ABY*I_NAI_I3x2yz_S_aa;
  Double I_NAI_I3xy2z_Py_aa = I_NAI_K3x2y2z_S_aa+ABY*I_NAI_I3xy2z_S_aa;
  Double I_NAI_I3x3z_Py_aa = I_NAI_K3xy3z_S_aa+ABY*I_NAI_I3x3z_S_aa;
  Double I_NAI_I2x4y_Py_aa = I_NAI_K2x5y_S_aa+ABY*I_NAI_I2x4y_S_aa;
  Double I_NAI_I2x3yz_Py_aa = I_NAI_K2x4yz_S_aa+ABY*I_NAI_I2x3yz_S_aa;
  Double I_NAI_I2x2y2z_Py_aa = I_NAI_K2x3y2z_S_aa+ABY*I_NAI_I2x2y2z_S_aa;
  Double I_NAI_I2xy3z_Py_aa = I_NAI_K2x2y3z_S_aa+ABY*I_NAI_I2xy3z_S_aa;
  Double I_NAI_I2x4z_Py_aa = I_NAI_K2xy4z_S_aa+ABY*I_NAI_I2x4z_S_aa;
  Double I_NAI_Ix5y_Py_aa = I_NAI_Kx6y_S_aa+ABY*I_NAI_Ix5y_S_aa;
  Double I_NAI_Ix4yz_Py_aa = I_NAI_Kx5yz_S_aa+ABY*I_NAI_Ix4yz_S_aa;
  Double I_NAI_Ix3y2z_Py_aa = I_NAI_Kx4y2z_S_aa+ABY*I_NAI_Ix3y2z_S_aa;
  Double I_NAI_Ix2y3z_Py_aa = I_NAI_Kx3y3z_S_aa+ABY*I_NAI_Ix2y3z_S_aa;
  Double I_NAI_Ixy4z_Py_aa = I_NAI_Kx2y4z_S_aa+ABY*I_NAI_Ixy4z_S_aa;
  Double I_NAI_Ix5z_Py_aa = I_NAI_Kxy5z_S_aa+ABY*I_NAI_Ix5z_S_aa;
  Double I_NAI_I6y_Py_aa = I_NAI_K7y_S_aa+ABY*I_NAI_I6y_S_aa;
  Double I_NAI_I5yz_Py_aa = I_NAI_K6yz_S_aa+ABY*I_NAI_I5yz_S_aa;
  Double I_NAI_I4y2z_Py_aa = I_NAI_K5y2z_S_aa+ABY*I_NAI_I4y2z_S_aa;
  Double I_NAI_I3y3z_Py_aa = I_NAI_K4y3z_S_aa+ABY*I_NAI_I3y3z_S_aa;
  Double I_NAI_I2y4z_Py_aa = I_NAI_K3y4z_S_aa+ABY*I_NAI_I2y4z_S_aa;
  Double I_NAI_Iy5z_Py_aa = I_NAI_K2y5z_S_aa+ABY*I_NAI_Iy5z_S_aa;
  Double I_NAI_I6z_Py_aa = I_NAI_Ky6z_S_aa+ABY*I_NAI_I6z_S_aa;
  Double I_NAI_I6x_Pz_aa = I_NAI_K6xz_S_aa+ABZ*I_NAI_I6x_S_aa;
  Double I_NAI_I5xy_Pz_aa = I_NAI_K5xyz_S_aa+ABZ*I_NAI_I5xy_S_aa;
  Double I_NAI_I5xz_Pz_aa = I_NAI_K5x2z_S_aa+ABZ*I_NAI_I5xz_S_aa;
  Double I_NAI_I4x2y_Pz_aa = I_NAI_K4x2yz_S_aa+ABZ*I_NAI_I4x2y_S_aa;
  Double I_NAI_I4xyz_Pz_aa = I_NAI_K4xy2z_S_aa+ABZ*I_NAI_I4xyz_S_aa;
  Double I_NAI_I4x2z_Pz_aa = I_NAI_K4x3z_S_aa+ABZ*I_NAI_I4x2z_S_aa;
  Double I_NAI_I3x3y_Pz_aa = I_NAI_K3x3yz_S_aa+ABZ*I_NAI_I3x3y_S_aa;
  Double I_NAI_I3x2yz_Pz_aa = I_NAI_K3x2y2z_S_aa+ABZ*I_NAI_I3x2yz_S_aa;
  Double I_NAI_I3xy2z_Pz_aa = I_NAI_K3xy3z_S_aa+ABZ*I_NAI_I3xy2z_S_aa;
  Double I_NAI_I3x3z_Pz_aa = I_NAI_K3x4z_S_aa+ABZ*I_NAI_I3x3z_S_aa;
  Double I_NAI_I2x4y_Pz_aa = I_NAI_K2x4yz_S_aa+ABZ*I_NAI_I2x4y_S_aa;
  Double I_NAI_I2x3yz_Pz_aa = I_NAI_K2x3y2z_S_aa+ABZ*I_NAI_I2x3yz_S_aa;
  Double I_NAI_I2x2y2z_Pz_aa = I_NAI_K2x2y3z_S_aa+ABZ*I_NAI_I2x2y2z_S_aa;
  Double I_NAI_I2xy3z_Pz_aa = I_NAI_K2xy4z_S_aa+ABZ*I_NAI_I2xy3z_S_aa;
  Double I_NAI_I2x4z_Pz_aa = I_NAI_K2x5z_S_aa+ABZ*I_NAI_I2x4z_S_aa;
  Double I_NAI_Ix5y_Pz_aa = I_NAI_Kx5yz_S_aa+ABZ*I_NAI_Ix5y_S_aa;
  Double I_NAI_Ix4yz_Pz_aa = I_NAI_Kx4y2z_S_aa+ABZ*I_NAI_Ix4yz_S_aa;
  Double I_NAI_Ix3y2z_Pz_aa = I_NAI_Kx3y3z_S_aa+ABZ*I_NAI_Ix3y2z_S_aa;
  Double I_NAI_Ix2y3z_Pz_aa = I_NAI_Kx2y4z_S_aa+ABZ*I_NAI_Ix2y3z_S_aa;
  Double I_NAI_Ixy4z_Pz_aa = I_NAI_Kxy5z_S_aa+ABZ*I_NAI_Ixy4z_S_aa;
  Double I_NAI_Ix5z_Pz_aa = I_NAI_Kx6z_S_aa+ABZ*I_NAI_Ix5z_S_aa;
  Double I_NAI_I6y_Pz_aa = I_NAI_K6yz_S_aa+ABZ*I_NAI_I6y_S_aa;
  Double I_NAI_I5yz_Pz_aa = I_NAI_K5y2z_S_aa+ABZ*I_NAI_I5yz_S_aa;
  Double I_NAI_I4y2z_Pz_aa = I_NAI_K4y3z_S_aa+ABZ*I_NAI_I4y2z_S_aa;
  Double I_NAI_I3y3z_Pz_aa = I_NAI_K3y4z_S_aa+ABZ*I_NAI_I3y3z_S_aa;
  Double I_NAI_I2y4z_Pz_aa = I_NAI_K2y5z_S_aa+ABZ*I_NAI_I2y4z_S_aa;
  Double I_NAI_Iy5z_Pz_aa = I_NAI_Ky6z_S_aa+ABZ*I_NAI_Iy5z_S_aa;
  Double I_NAI_I6z_Pz_aa = I_NAI_K7z_S_aa+ABZ*I_NAI_I6z_S_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_K_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S_aa
   * RHS shell quartet name: SQ_NAI_K_S_aa
   ************************************************************/
  Double I_NAI_K7x_Px_aa = I_NAI_L8x_S_aa+ABX*I_NAI_K7x_S_aa;
  Double I_NAI_K6xy_Px_aa = I_NAI_L7xy_S_aa+ABX*I_NAI_K6xy_S_aa;
  Double I_NAI_K6xz_Px_aa = I_NAI_L7xz_S_aa+ABX*I_NAI_K6xz_S_aa;
  Double I_NAI_K5x2y_Px_aa = I_NAI_L6x2y_S_aa+ABX*I_NAI_K5x2y_S_aa;
  Double I_NAI_K5xyz_Px_aa = I_NAI_L6xyz_S_aa+ABX*I_NAI_K5xyz_S_aa;
  Double I_NAI_K5x2z_Px_aa = I_NAI_L6x2z_S_aa+ABX*I_NAI_K5x2z_S_aa;
  Double I_NAI_K4x3y_Px_aa = I_NAI_L5x3y_S_aa+ABX*I_NAI_K4x3y_S_aa;
  Double I_NAI_K4x2yz_Px_aa = I_NAI_L5x2yz_S_aa+ABX*I_NAI_K4x2yz_S_aa;
  Double I_NAI_K4xy2z_Px_aa = I_NAI_L5xy2z_S_aa+ABX*I_NAI_K4xy2z_S_aa;
  Double I_NAI_K4x3z_Px_aa = I_NAI_L5x3z_S_aa+ABX*I_NAI_K4x3z_S_aa;
  Double I_NAI_K3x4y_Px_aa = I_NAI_L4x4y_S_aa+ABX*I_NAI_K3x4y_S_aa;
  Double I_NAI_K3x3yz_Px_aa = I_NAI_L4x3yz_S_aa+ABX*I_NAI_K3x3yz_S_aa;
  Double I_NAI_K3x2y2z_Px_aa = I_NAI_L4x2y2z_S_aa+ABX*I_NAI_K3x2y2z_S_aa;
  Double I_NAI_K3xy3z_Px_aa = I_NAI_L4xy3z_S_aa+ABX*I_NAI_K3xy3z_S_aa;
  Double I_NAI_K3x4z_Px_aa = I_NAI_L4x4z_S_aa+ABX*I_NAI_K3x4z_S_aa;
  Double I_NAI_K2x5y_Px_aa = I_NAI_L3x5y_S_aa+ABX*I_NAI_K2x5y_S_aa;
  Double I_NAI_K2x4yz_Px_aa = I_NAI_L3x4yz_S_aa+ABX*I_NAI_K2x4yz_S_aa;
  Double I_NAI_K2x3y2z_Px_aa = I_NAI_L3x3y2z_S_aa+ABX*I_NAI_K2x3y2z_S_aa;
  Double I_NAI_K2x2y3z_Px_aa = I_NAI_L3x2y3z_S_aa+ABX*I_NAI_K2x2y3z_S_aa;
  Double I_NAI_K2xy4z_Px_aa = I_NAI_L3xy4z_S_aa+ABX*I_NAI_K2xy4z_S_aa;
  Double I_NAI_K2x5z_Px_aa = I_NAI_L3x5z_S_aa+ABX*I_NAI_K2x5z_S_aa;
  Double I_NAI_Kx6y_Px_aa = I_NAI_L2x6y_S_aa+ABX*I_NAI_Kx6y_S_aa;
  Double I_NAI_Kx5yz_Px_aa = I_NAI_L2x5yz_S_aa+ABX*I_NAI_Kx5yz_S_aa;
  Double I_NAI_Kx4y2z_Px_aa = I_NAI_L2x4y2z_S_aa+ABX*I_NAI_Kx4y2z_S_aa;
  Double I_NAI_Kx3y3z_Px_aa = I_NAI_L2x3y3z_S_aa+ABX*I_NAI_Kx3y3z_S_aa;
  Double I_NAI_Kx2y4z_Px_aa = I_NAI_L2x2y4z_S_aa+ABX*I_NAI_Kx2y4z_S_aa;
  Double I_NAI_Kxy5z_Px_aa = I_NAI_L2xy5z_S_aa+ABX*I_NAI_Kxy5z_S_aa;
  Double I_NAI_Kx6z_Px_aa = I_NAI_L2x6z_S_aa+ABX*I_NAI_Kx6z_S_aa;
  Double I_NAI_K7y_Px_aa = I_NAI_Lx7y_S_aa+ABX*I_NAI_K7y_S_aa;
  Double I_NAI_K6yz_Px_aa = I_NAI_Lx6yz_S_aa+ABX*I_NAI_K6yz_S_aa;
  Double I_NAI_K5y2z_Px_aa = I_NAI_Lx5y2z_S_aa+ABX*I_NAI_K5y2z_S_aa;
  Double I_NAI_K4y3z_Px_aa = I_NAI_Lx4y3z_S_aa+ABX*I_NAI_K4y3z_S_aa;
  Double I_NAI_K3y4z_Px_aa = I_NAI_Lx3y4z_S_aa+ABX*I_NAI_K3y4z_S_aa;
  Double I_NAI_K2y5z_Px_aa = I_NAI_Lx2y5z_S_aa+ABX*I_NAI_K2y5z_S_aa;
  Double I_NAI_Ky6z_Px_aa = I_NAI_Lxy6z_S_aa+ABX*I_NAI_Ky6z_S_aa;
  Double I_NAI_K7z_Px_aa = I_NAI_Lx7z_S_aa+ABX*I_NAI_K7z_S_aa;
  Double I_NAI_K7x_Py_aa = I_NAI_L7xy_S_aa+ABY*I_NAI_K7x_S_aa;
  Double I_NAI_K6xy_Py_aa = I_NAI_L6x2y_S_aa+ABY*I_NAI_K6xy_S_aa;
  Double I_NAI_K6xz_Py_aa = I_NAI_L6xyz_S_aa+ABY*I_NAI_K6xz_S_aa;
  Double I_NAI_K5x2y_Py_aa = I_NAI_L5x3y_S_aa+ABY*I_NAI_K5x2y_S_aa;
  Double I_NAI_K5xyz_Py_aa = I_NAI_L5x2yz_S_aa+ABY*I_NAI_K5xyz_S_aa;
  Double I_NAI_K5x2z_Py_aa = I_NAI_L5xy2z_S_aa+ABY*I_NAI_K5x2z_S_aa;
  Double I_NAI_K4x3y_Py_aa = I_NAI_L4x4y_S_aa+ABY*I_NAI_K4x3y_S_aa;
  Double I_NAI_K4x2yz_Py_aa = I_NAI_L4x3yz_S_aa+ABY*I_NAI_K4x2yz_S_aa;
  Double I_NAI_K4xy2z_Py_aa = I_NAI_L4x2y2z_S_aa+ABY*I_NAI_K4xy2z_S_aa;
  Double I_NAI_K4x3z_Py_aa = I_NAI_L4xy3z_S_aa+ABY*I_NAI_K4x3z_S_aa;
  Double I_NAI_K3x4y_Py_aa = I_NAI_L3x5y_S_aa+ABY*I_NAI_K3x4y_S_aa;
  Double I_NAI_K3x3yz_Py_aa = I_NAI_L3x4yz_S_aa+ABY*I_NAI_K3x3yz_S_aa;
  Double I_NAI_K3x2y2z_Py_aa = I_NAI_L3x3y2z_S_aa+ABY*I_NAI_K3x2y2z_S_aa;
  Double I_NAI_K3xy3z_Py_aa = I_NAI_L3x2y3z_S_aa+ABY*I_NAI_K3xy3z_S_aa;
  Double I_NAI_K3x4z_Py_aa = I_NAI_L3xy4z_S_aa+ABY*I_NAI_K3x4z_S_aa;
  Double I_NAI_K2x5y_Py_aa = I_NAI_L2x6y_S_aa+ABY*I_NAI_K2x5y_S_aa;
  Double I_NAI_K2x4yz_Py_aa = I_NAI_L2x5yz_S_aa+ABY*I_NAI_K2x4yz_S_aa;
  Double I_NAI_K2x3y2z_Py_aa = I_NAI_L2x4y2z_S_aa+ABY*I_NAI_K2x3y2z_S_aa;
  Double I_NAI_K2x2y3z_Py_aa = I_NAI_L2x3y3z_S_aa+ABY*I_NAI_K2x2y3z_S_aa;
  Double I_NAI_K2xy4z_Py_aa = I_NAI_L2x2y4z_S_aa+ABY*I_NAI_K2xy4z_S_aa;
  Double I_NAI_K2x5z_Py_aa = I_NAI_L2xy5z_S_aa+ABY*I_NAI_K2x5z_S_aa;
  Double I_NAI_Kx6y_Py_aa = I_NAI_Lx7y_S_aa+ABY*I_NAI_Kx6y_S_aa;
  Double I_NAI_Kx5yz_Py_aa = I_NAI_Lx6yz_S_aa+ABY*I_NAI_Kx5yz_S_aa;
  Double I_NAI_Kx4y2z_Py_aa = I_NAI_Lx5y2z_S_aa+ABY*I_NAI_Kx4y2z_S_aa;
  Double I_NAI_Kx3y3z_Py_aa = I_NAI_Lx4y3z_S_aa+ABY*I_NAI_Kx3y3z_S_aa;
  Double I_NAI_Kx2y4z_Py_aa = I_NAI_Lx3y4z_S_aa+ABY*I_NAI_Kx2y4z_S_aa;
  Double I_NAI_Kxy5z_Py_aa = I_NAI_Lx2y5z_S_aa+ABY*I_NAI_Kxy5z_S_aa;
  Double I_NAI_Kx6z_Py_aa = I_NAI_Lxy6z_S_aa+ABY*I_NAI_Kx6z_S_aa;
  Double I_NAI_K7y_Py_aa = I_NAI_L8y_S_aa+ABY*I_NAI_K7y_S_aa;
  Double I_NAI_K6yz_Py_aa = I_NAI_L7yz_S_aa+ABY*I_NAI_K6yz_S_aa;
  Double I_NAI_K5y2z_Py_aa = I_NAI_L6y2z_S_aa+ABY*I_NAI_K5y2z_S_aa;
  Double I_NAI_K4y3z_Py_aa = I_NAI_L5y3z_S_aa+ABY*I_NAI_K4y3z_S_aa;
  Double I_NAI_K3y4z_Py_aa = I_NAI_L4y4z_S_aa+ABY*I_NAI_K3y4z_S_aa;
  Double I_NAI_K2y5z_Py_aa = I_NAI_L3y5z_S_aa+ABY*I_NAI_K2y5z_S_aa;
  Double I_NAI_Ky6z_Py_aa = I_NAI_L2y6z_S_aa+ABY*I_NAI_Ky6z_S_aa;
  Double I_NAI_K7z_Py_aa = I_NAI_Ly7z_S_aa+ABY*I_NAI_K7z_S_aa;
  Double I_NAI_K7x_Pz_aa = I_NAI_L7xz_S_aa+ABZ*I_NAI_K7x_S_aa;
  Double I_NAI_K6xy_Pz_aa = I_NAI_L6xyz_S_aa+ABZ*I_NAI_K6xy_S_aa;
  Double I_NAI_K6xz_Pz_aa = I_NAI_L6x2z_S_aa+ABZ*I_NAI_K6xz_S_aa;
  Double I_NAI_K5x2y_Pz_aa = I_NAI_L5x2yz_S_aa+ABZ*I_NAI_K5x2y_S_aa;
  Double I_NAI_K5xyz_Pz_aa = I_NAI_L5xy2z_S_aa+ABZ*I_NAI_K5xyz_S_aa;
  Double I_NAI_K5x2z_Pz_aa = I_NAI_L5x3z_S_aa+ABZ*I_NAI_K5x2z_S_aa;
  Double I_NAI_K4x3y_Pz_aa = I_NAI_L4x3yz_S_aa+ABZ*I_NAI_K4x3y_S_aa;
  Double I_NAI_K4x2yz_Pz_aa = I_NAI_L4x2y2z_S_aa+ABZ*I_NAI_K4x2yz_S_aa;
  Double I_NAI_K4xy2z_Pz_aa = I_NAI_L4xy3z_S_aa+ABZ*I_NAI_K4xy2z_S_aa;
  Double I_NAI_K4x3z_Pz_aa = I_NAI_L4x4z_S_aa+ABZ*I_NAI_K4x3z_S_aa;
  Double I_NAI_K3x4y_Pz_aa = I_NAI_L3x4yz_S_aa+ABZ*I_NAI_K3x4y_S_aa;
  Double I_NAI_K3x3yz_Pz_aa = I_NAI_L3x3y2z_S_aa+ABZ*I_NAI_K3x3yz_S_aa;
  Double I_NAI_K3x2y2z_Pz_aa = I_NAI_L3x2y3z_S_aa+ABZ*I_NAI_K3x2y2z_S_aa;
  Double I_NAI_K3xy3z_Pz_aa = I_NAI_L3xy4z_S_aa+ABZ*I_NAI_K3xy3z_S_aa;
  Double I_NAI_K3x4z_Pz_aa = I_NAI_L3x5z_S_aa+ABZ*I_NAI_K3x4z_S_aa;
  Double I_NAI_K2x5y_Pz_aa = I_NAI_L2x5yz_S_aa+ABZ*I_NAI_K2x5y_S_aa;
  Double I_NAI_K2x4yz_Pz_aa = I_NAI_L2x4y2z_S_aa+ABZ*I_NAI_K2x4yz_S_aa;
  Double I_NAI_K2x3y2z_Pz_aa = I_NAI_L2x3y3z_S_aa+ABZ*I_NAI_K2x3y2z_S_aa;
  Double I_NAI_K2x2y3z_Pz_aa = I_NAI_L2x2y4z_S_aa+ABZ*I_NAI_K2x2y3z_S_aa;
  Double I_NAI_K2xy4z_Pz_aa = I_NAI_L2xy5z_S_aa+ABZ*I_NAI_K2xy4z_S_aa;
  Double I_NAI_K2x5z_Pz_aa = I_NAI_L2x6z_S_aa+ABZ*I_NAI_K2x5z_S_aa;
  Double I_NAI_Kx6y_Pz_aa = I_NAI_Lx6yz_S_aa+ABZ*I_NAI_Kx6y_S_aa;
  Double I_NAI_Kx5yz_Pz_aa = I_NAI_Lx5y2z_S_aa+ABZ*I_NAI_Kx5yz_S_aa;
  Double I_NAI_Kx4y2z_Pz_aa = I_NAI_Lx4y3z_S_aa+ABZ*I_NAI_Kx4y2z_S_aa;
  Double I_NAI_Kx3y3z_Pz_aa = I_NAI_Lx3y4z_S_aa+ABZ*I_NAI_Kx3y3z_S_aa;
  Double I_NAI_Kx2y4z_Pz_aa = I_NAI_Lx2y5z_S_aa+ABZ*I_NAI_Kx2y4z_S_aa;
  Double I_NAI_Kxy5z_Pz_aa = I_NAI_Lxy6z_S_aa+ABZ*I_NAI_Kxy5z_S_aa;
  Double I_NAI_Kx6z_Pz_aa = I_NAI_Lx7z_S_aa+ABZ*I_NAI_Kx6z_S_aa;
  Double I_NAI_K7y_Pz_aa = I_NAI_L7yz_S_aa+ABZ*I_NAI_K7y_S_aa;
  Double I_NAI_K6yz_Pz_aa = I_NAI_L6y2z_S_aa+ABZ*I_NAI_K6yz_S_aa;
  Double I_NAI_K5y2z_Pz_aa = I_NAI_L5y3z_S_aa+ABZ*I_NAI_K5y2z_S_aa;
  Double I_NAI_K4y3z_Pz_aa = I_NAI_L4y4z_S_aa+ABZ*I_NAI_K4y3z_S_aa;
  Double I_NAI_K3y4z_Pz_aa = I_NAI_L3y5z_S_aa+ABZ*I_NAI_K3y4z_S_aa;
  Double I_NAI_K2y5z_Pz_aa = I_NAI_L2y6z_S_aa+ABZ*I_NAI_K2y5z_S_aa;
  Double I_NAI_Ky6z_Pz_aa = I_NAI_Ly7z_S_aa+ABZ*I_NAI_Ky6z_S_aa;
  Double I_NAI_K7z_Pz_aa = I_NAI_L8z_S_aa+ABZ*I_NAI_K7z_S_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 56 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_I_P_aa
   ************************************************************/
  Double I_NAI_I6x_D2x_aa = I_NAI_K7x_Px_aa+ABX*I_NAI_I6x_Px_aa;
  Double I_NAI_I5xy_D2x_aa = I_NAI_K6xy_Px_aa+ABX*I_NAI_I5xy_Px_aa;
  Double I_NAI_I5xz_D2x_aa = I_NAI_K6xz_Px_aa+ABX*I_NAI_I5xz_Px_aa;
  Double I_NAI_I4x2y_D2x_aa = I_NAI_K5x2y_Px_aa+ABX*I_NAI_I4x2y_Px_aa;
  Double I_NAI_I4xyz_D2x_aa = I_NAI_K5xyz_Px_aa+ABX*I_NAI_I4xyz_Px_aa;
  Double I_NAI_I4x2z_D2x_aa = I_NAI_K5x2z_Px_aa+ABX*I_NAI_I4x2z_Px_aa;
  Double I_NAI_I3x3y_D2x_aa = I_NAI_K4x3y_Px_aa+ABX*I_NAI_I3x3y_Px_aa;
  Double I_NAI_I3x2yz_D2x_aa = I_NAI_K4x2yz_Px_aa+ABX*I_NAI_I3x2yz_Px_aa;
  Double I_NAI_I3xy2z_D2x_aa = I_NAI_K4xy2z_Px_aa+ABX*I_NAI_I3xy2z_Px_aa;
  Double I_NAI_I3x3z_D2x_aa = I_NAI_K4x3z_Px_aa+ABX*I_NAI_I3x3z_Px_aa;
  Double I_NAI_I2x4y_D2x_aa = I_NAI_K3x4y_Px_aa+ABX*I_NAI_I2x4y_Px_aa;
  Double I_NAI_I2x3yz_D2x_aa = I_NAI_K3x3yz_Px_aa+ABX*I_NAI_I2x3yz_Px_aa;
  Double I_NAI_I2x2y2z_D2x_aa = I_NAI_K3x2y2z_Px_aa+ABX*I_NAI_I2x2y2z_Px_aa;
  Double I_NAI_I2xy3z_D2x_aa = I_NAI_K3xy3z_Px_aa+ABX*I_NAI_I2xy3z_Px_aa;
  Double I_NAI_I2x4z_D2x_aa = I_NAI_K3x4z_Px_aa+ABX*I_NAI_I2x4z_Px_aa;
  Double I_NAI_Ix5y_D2x_aa = I_NAI_K2x5y_Px_aa+ABX*I_NAI_Ix5y_Px_aa;
  Double I_NAI_Ix4yz_D2x_aa = I_NAI_K2x4yz_Px_aa+ABX*I_NAI_Ix4yz_Px_aa;
  Double I_NAI_Ix3y2z_D2x_aa = I_NAI_K2x3y2z_Px_aa+ABX*I_NAI_Ix3y2z_Px_aa;
  Double I_NAI_Ix2y3z_D2x_aa = I_NAI_K2x2y3z_Px_aa+ABX*I_NAI_Ix2y3z_Px_aa;
  Double I_NAI_Ixy4z_D2x_aa = I_NAI_K2xy4z_Px_aa+ABX*I_NAI_Ixy4z_Px_aa;
  Double I_NAI_Ix5z_D2x_aa = I_NAI_K2x5z_Px_aa+ABX*I_NAI_Ix5z_Px_aa;
  Double I_NAI_I6y_D2x_aa = I_NAI_Kx6y_Px_aa+ABX*I_NAI_I6y_Px_aa;
  Double I_NAI_I5yz_D2x_aa = I_NAI_Kx5yz_Px_aa+ABX*I_NAI_I5yz_Px_aa;
  Double I_NAI_I4y2z_D2x_aa = I_NAI_Kx4y2z_Px_aa+ABX*I_NAI_I4y2z_Px_aa;
  Double I_NAI_I3y3z_D2x_aa = I_NAI_Kx3y3z_Px_aa+ABX*I_NAI_I3y3z_Px_aa;
  Double I_NAI_I2y4z_D2x_aa = I_NAI_Kx2y4z_Px_aa+ABX*I_NAI_I2y4z_Px_aa;
  Double I_NAI_Iy5z_D2x_aa = I_NAI_Kxy5z_Px_aa+ABX*I_NAI_Iy5z_Px_aa;
  Double I_NAI_I6z_D2x_aa = I_NAI_Kx6z_Px_aa+ABX*I_NAI_I6z_Px_aa;
  Double I_NAI_I6x_Dxy_aa = I_NAI_K6xy_Px_aa+ABY*I_NAI_I6x_Px_aa;
  Double I_NAI_I5xy_Dxy_aa = I_NAI_K5x2y_Px_aa+ABY*I_NAI_I5xy_Px_aa;
  Double I_NAI_I5xz_Dxy_aa = I_NAI_K5xyz_Px_aa+ABY*I_NAI_I5xz_Px_aa;
  Double I_NAI_I4x2y_Dxy_aa = I_NAI_K4x3y_Px_aa+ABY*I_NAI_I4x2y_Px_aa;
  Double I_NAI_I4xyz_Dxy_aa = I_NAI_K4x2yz_Px_aa+ABY*I_NAI_I4xyz_Px_aa;
  Double I_NAI_I4x2z_Dxy_aa = I_NAI_K4xy2z_Px_aa+ABY*I_NAI_I4x2z_Px_aa;
  Double I_NAI_I3x3y_Dxy_aa = I_NAI_K3x4y_Px_aa+ABY*I_NAI_I3x3y_Px_aa;
  Double I_NAI_I3x2yz_Dxy_aa = I_NAI_K3x3yz_Px_aa+ABY*I_NAI_I3x2yz_Px_aa;
  Double I_NAI_I3xy2z_Dxy_aa = I_NAI_K3x2y2z_Px_aa+ABY*I_NAI_I3xy2z_Px_aa;
  Double I_NAI_I3x3z_Dxy_aa = I_NAI_K3xy3z_Px_aa+ABY*I_NAI_I3x3z_Px_aa;
  Double I_NAI_I2x4y_Dxy_aa = I_NAI_K2x5y_Px_aa+ABY*I_NAI_I2x4y_Px_aa;
  Double I_NAI_I2x3yz_Dxy_aa = I_NAI_K2x4yz_Px_aa+ABY*I_NAI_I2x3yz_Px_aa;
  Double I_NAI_I2x2y2z_Dxy_aa = I_NAI_K2x3y2z_Px_aa+ABY*I_NAI_I2x2y2z_Px_aa;
  Double I_NAI_I2xy3z_Dxy_aa = I_NAI_K2x2y3z_Px_aa+ABY*I_NAI_I2xy3z_Px_aa;
  Double I_NAI_I2x4z_Dxy_aa = I_NAI_K2xy4z_Px_aa+ABY*I_NAI_I2x4z_Px_aa;
  Double I_NAI_Ix5y_Dxy_aa = I_NAI_Kx6y_Px_aa+ABY*I_NAI_Ix5y_Px_aa;
  Double I_NAI_Ix4yz_Dxy_aa = I_NAI_Kx5yz_Px_aa+ABY*I_NAI_Ix4yz_Px_aa;
  Double I_NAI_Ix3y2z_Dxy_aa = I_NAI_Kx4y2z_Px_aa+ABY*I_NAI_Ix3y2z_Px_aa;
  Double I_NAI_Ix2y3z_Dxy_aa = I_NAI_Kx3y3z_Px_aa+ABY*I_NAI_Ix2y3z_Px_aa;
  Double I_NAI_Ixy4z_Dxy_aa = I_NAI_Kx2y4z_Px_aa+ABY*I_NAI_Ixy4z_Px_aa;
  Double I_NAI_Ix5z_Dxy_aa = I_NAI_Kxy5z_Px_aa+ABY*I_NAI_Ix5z_Px_aa;
  Double I_NAI_I6y_Dxy_aa = I_NAI_K7y_Px_aa+ABY*I_NAI_I6y_Px_aa;
  Double I_NAI_I5yz_Dxy_aa = I_NAI_K6yz_Px_aa+ABY*I_NAI_I5yz_Px_aa;
  Double I_NAI_I4y2z_Dxy_aa = I_NAI_K5y2z_Px_aa+ABY*I_NAI_I4y2z_Px_aa;
  Double I_NAI_I3y3z_Dxy_aa = I_NAI_K4y3z_Px_aa+ABY*I_NAI_I3y3z_Px_aa;
  Double I_NAI_I2y4z_Dxy_aa = I_NAI_K3y4z_Px_aa+ABY*I_NAI_I2y4z_Px_aa;
  Double I_NAI_Iy5z_Dxy_aa = I_NAI_K2y5z_Px_aa+ABY*I_NAI_Iy5z_Px_aa;
  Double I_NAI_I6z_Dxy_aa = I_NAI_Ky6z_Px_aa+ABY*I_NAI_I6z_Px_aa;
  Double I_NAI_I6x_D2y_aa = I_NAI_K6xy_Py_aa+ABY*I_NAI_I6x_Py_aa;
  Double I_NAI_I5xy_D2y_aa = I_NAI_K5x2y_Py_aa+ABY*I_NAI_I5xy_Py_aa;
  Double I_NAI_I5xz_D2y_aa = I_NAI_K5xyz_Py_aa+ABY*I_NAI_I5xz_Py_aa;
  Double I_NAI_I4x2y_D2y_aa = I_NAI_K4x3y_Py_aa+ABY*I_NAI_I4x2y_Py_aa;
  Double I_NAI_I4xyz_D2y_aa = I_NAI_K4x2yz_Py_aa+ABY*I_NAI_I4xyz_Py_aa;
  Double I_NAI_I4x2z_D2y_aa = I_NAI_K4xy2z_Py_aa+ABY*I_NAI_I4x2z_Py_aa;
  Double I_NAI_I3x3y_D2y_aa = I_NAI_K3x4y_Py_aa+ABY*I_NAI_I3x3y_Py_aa;
  Double I_NAI_I3x2yz_D2y_aa = I_NAI_K3x3yz_Py_aa+ABY*I_NAI_I3x2yz_Py_aa;
  Double I_NAI_I3xy2z_D2y_aa = I_NAI_K3x2y2z_Py_aa+ABY*I_NAI_I3xy2z_Py_aa;
  Double I_NAI_I3x3z_D2y_aa = I_NAI_K3xy3z_Py_aa+ABY*I_NAI_I3x3z_Py_aa;
  Double I_NAI_I2x4y_D2y_aa = I_NAI_K2x5y_Py_aa+ABY*I_NAI_I2x4y_Py_aa;
  Double I_NAI_I2x3yz_D2y_aa = I_NAI_K2x4yz_Py_aa+ABY*I_NAI_I2x3yz_Py_aa;
  Double I_NAI_I2x2y2z_D2y_aa = I_NAI_K2x3y2z_Py_aa+ABY*I_NAI_I2x2y2z_Py_aa;
  Double I_NAI_I2xy3z_D2y_aa = I_NAI_K2x2y3z_Py_aa+ABY*I_NAI_I2xy3z_Py_aa;
  Double I_NAI_I2x4z_D2y_aa = I_NAI_K2xy4z_Py_aa+ABY*I_NAI_I2x4z_Py_aa;
  Double I_NAI_Ix5y_D2y_aa = I_NAI_Kx6y_Py_aa+ABY*I_NAI_Ix5y_Py_aa;
  Double I_NAI_Ix4yz_D2y_aa = I_NAI_Kx5yz_Py_aa+ABY*I_NAI_Ix4yz_Py_aa;
  Double I_NAI_Ix3y2z_D2y_aa = I_NAI_Kx4y2z_Py_aa+ABY*I_NAI_Ix3y2z_Py_aa;
  Double I_NAI_Ix2y3z_D2y_aa = I_NAI_Kx3y3z_Py_aa+ABY*I_NAI_Ix2y3z_Py_aa;
  Double I_NAI_Ixy4z_D2y_aa = I_NAI_Kx2y4z_Py_aa+ABY*I_NAI_Ixy4z_Py_aa;
  Double I_NAI_Ix5z_D2y_aa = I_NAI_Kxy5z_Py_aa+ABY*I_NAI_Ix5z_Py_aa;
  Double I_NAI_I6y_D2y_aa = I_NAI_K7y_Py_aa+ABY*I_NAI_I6y_Py_aa;
  Double I_NAI_I5yz_D2y_aa = I_NAI_K6yz_Py_aa+ABY*I_NAI_I5yz_Py_aa;
  Double I_NAI_I4y2z_D2y_aa = I_NAI_K5y2z_Py_aa+ABY*I_NAI_I4y2z_Py_aa;
  Double I_NAI_I3y3z_D2y_aa = I_NAI_K4y3z_Py_aa+ABY*I_NAI_I3y3z_Py_aa;
  Double I_NAI_I2y4z_D2y_aa = I_NAI_K3y4z_Py_aa+ABY*I_NAI_I2y4z_Py_aa;
  Double I_NAI_Iy5z_D2y_aa = I_NAI_K2y5z_Py_aa+ABY*I_NAI_Iy5z_Py_aa;
  Double I_NAI_I6z_D2y_aa = I_NAI_Ky6z_Py_aa+ABY*I_NAI_I6z_Py_aa;
  Double I_NAI_I6x_D2z_aa = I_NAI_K6xz_Pz_aa+ABZ*I_NAI_I6x_Pz_aa;
  Double I_NAI_I5xy_D2z_aa = I_NAI_K5xyz_Pz_aa+ABZ*I_NAI_I5xy_Pz_aa;
  Double I_NAI_I5xz_D2z_aa = I_NAI_K5x2z_Pz_aa+ABZ*I_NAI_I5xz_Pz_aa;
  Double I_NAI_I4x2y_D2z_aa = I_NAI_K4x2yz_Pz_aa+ABZ*I_NAI_I4x2y_Pz_aa;
  Double I_NAI_I4xyz_D2z_aa = I_NAI_K4xy2z_Pz_aa+ABZ*I_NAI_I4xyz_Pz_aa;
  Double I_NAI_I4x2z_D2z_aa = I_NAI_K4x3z_Pz_aa+ABZ*I_NAI_I4x2z_Pz_aa;
  Double I_NAI_I3x3y_D2z_aa = I_NAI_K3x3yz_Pz_aa+ABZ*I_NAI_I3x3y_Pz_aa;
  Double I_NAI_I3x2yz_D2z_aa = I_NAI_K3x2y2z_Pz_aa+ABZ*I_NAI_I3x2yz_Pz_aa;
  Double I_NAI_I3xy2z_D2z_aa = I_NAI_K3xy3z_Pz_aa+ABZ*I_NAI_I3xy2z_Pz_aa;
  Double I_NAI_I3x3z_D2z_aa = I_NAI_K3x4z_Pz_aa+ABZ*I_NAI_I3x3z_Pz_aa;
  Double I_NAI_I2x4y_D2z_aa = I_NAI_K2x4yz_Pz_aa+ABZ*I_NAI_I2x4y_Pz_aa;
  Double I_NAI_I2x3yz_D2z_aa = I_NAI_K2x3y2z_Pz_aa+ABZ*I_NAI_I2x3yz_Pz_aa;
  Double I_NAI_I2x2y2z_D2z_aa = I_NAI_K2x2y3z_Pz_aa+ABZ*I_NAI_I2x2y2z_Pz_aa;
  Double I_NAI_I2xy3z_D2z_aa = I_NAI_K2xy4z_Pz_aa+ABZ*I_NAI_I2xy3z_Pz_aa;
  Double I_NAI_I2x4z_D2z_aa = I_NAI_K2x5z_Pz_aa+ABZ*I_NAI_I2x4z_Pz_aa;
  Double I_NAI_Ix5y_D2z_aa = I_NAI_Kx5yz_Pz_aa+ABZ*I_NAI_Ix5y_Pz_aa;
  Double I_NAI_Ix4yz_D2z_aa = I_NAI_Kx4y2z_Pz_aa+ABZ*I_NAI_Ix4yz_Pz_aa;
  Double I_NAI_Ix3y2z_D2z_aa = I_NAI_Kx3y3z_Pz_aa+ABZ*I_NAI_Ix3y2z_Pz_aa;
  Double I_NAI_Ix2y3z_D2z_aa = I_NAI_Kx2y4z_Pz_aa+ABZ*I_NAI_Ix2y3z_Pz_aa;
  Double I_NAI_Ixy4z_D2z_aa = I_NAI_Kxy5z_Pz_aa+ABZ*I_NAI_Ixy4z_Pz_aa;
  Double I_NAI_Ix5z_D2z_aa = I_NAI_Kx6z_Pz_aa+ABZ*I_NAI_Ix5z_Pz_aa;
  Double I_NAI_I6y_D2z_aa = I_NAI_K6yz_Pz_aa+ABZ*I_NAI_I6y_Pz_aa;
  Double I_NAI_I5yz_D2z_aa = I_NAI_K5y2z_Pz_aa+ABZ*I_NAI_I5yz_Pz_aa;
  Double I_NAI_I4y2z_D2z_aa = I_NAI_K4y3z_Pz_aa+ABZ*I_NAI_I4y2z_Pz_aa;
  Double I_NAI_I3y3z_D2z_aa = I_NAI_K3y4z_Pz_aa+ABZ*I_NAI_I3y3z_Pz_aa;
  Double I_NAI_I2y4z_D2z_aa = I_NAI_K2y5z_Pz_aa+ABZ*I_NAI_I2y4z_Pz_aa;
  Double I_NAI_Iy5z_D2z_aa = I_NAI_Ky6z_Pz_aa+ABZ*I_NAI_Iy5z_Pz_aa;
  Double I_NAI_I6z_D2z_aa = I_NAI_K7z_Pz_aa+ABZ*I_NAI_I6z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_L_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_S_aa
   * RHS shell quartet name: SQ_NAI_L_S_aa
   ************************************************************/
  Double I_NAI_L8x_Px_aa = I_NAI_M9x_S_aa+ABX*I_NAI_L8x_S_aa;
  Double I_NAI_L7xy_Px_aa = I_NAI_M8xy_S_aa+ABX*I_NAI_L7xy_S_aa;
  Double I_NAI_L7xz_Px_aa = I_NAI_M8xz_S_aa+ABX*I_NAI_L7xz_S_aa;
  Double I_NAI_L6x2y_Px_aa = I_NAI_M7x2y_S_aa+ABX*I_NAI_L6x2y_S_aa;
  Double I_NAI_L6xyz_Px_aa = I_NAI_M7xyz_S_aa+ABX*I_NAI_L6xyz_S_aa;
  Double I_NAI_L6x2z_Px_aa = I_NAI_M7x2z_S_aa+ABX*I_NAI_L6x2z_S_aa;
  Double I_NAI_L5x3y_Px_aa = I_NAI_M6x3y_S_aa+ABX*I_NAI_L5x3y_S_aa;
  Double I_NAI_L5x2yz_Px_aa = I_NAI_M6x2yz_S_aa+ABX*I_NAI_L5x2yz_S_aa;
  Double I_NAI_L5xy2z_Px_aa = I_NAI_M6xy2z_S_aa+ABX*I_NAI_L5xy2z_S_aa;
  Double I_NAI_L5x3z_Px_aa = I_NAI_M6x3z_S_aa+ABX*I_NAI_L5x3z_S_aa;
  Double I_NAI_L4x4y_Px_aa = I_NAI_M5x4y_S_aa+ABX*I_NAI_L4x4y_S_aa;
  Double I_NAI_L4x3yz_Px_aa = I_NAI_M5x3yz_S_aa+ABX*I_NAI_L4x3yz_S_aa;
  Double I_NAI_L4x2y2z_Px_aa = I_NAI_M5x2y2z_S_aa+ABX*I_NAI_L4x2y2z_S_aa;
  Double I_NAI_L4xy3z_Px_aa = I_NAI_M5xy3z_S_aa+ABX*I_NAI_L4xy3z_S_aa;
  Double I_NAI_L4x4z_Px_aa = I_NAI_M5x4z_S_aa+ABX*I_NAI_L4x4z_S_aa;
  Double I_NAI_L3x5y_Px_aa = I_NAI_M4x5y_S_aa+ABX*I_NAI_L3x5y_S_aa;
  Double I_NAI_L3x4yz_Px_aa = I_NAI_M4x4yz_S_aa+ABX*I_NAI_L3x4yz_S_aa;
  Double I_NAI_L3x3y2z_Px_aa = I_NAI_M4x3y2z_S_aa+ABX*I_NAI_L3x3y2z_S_aa;
  Double I_NAI_L3x2y3z_Px_aa = I_NAI_M4x2y3z_S_aa+ABX*I_NAI_L3x2y3z_S_aa;
  Double I_NAI_L3xy4z_Px_aa = I_NAI_M4xy4z_S_aa+ABX*I_NAI_L3xy4z_S_aa;
  Double I_NAI_L3x5z_Px_aa = I_NAI_M4x5z_S_aa+ABX*I_NAI_L3x5z_S_aa;
  Double I_NAI_L2x6y_Px_aa = I_NAI_M3x6y_S_aa+ABX*I_NAI_L2x6y_S_aa;
  Double I_NAI_L2x5yz_Px_aa = I_NAI_M3x5yz_S_aa+ABX*I_NAI_L2x5yz_S_aa;
  Double I_NAI_L2x4y2z_Px_aa = I_NAI_M3x4y2z_S_aa+ABX*I_NAI_L2x4y2z_S_aa;
  Double I_NAI_L2x3y3z_Px_aa = I_NAI_M3x3y3z_S_aa+ABX*I_NAI_L2x3y3z_S_aa;
  Double I_NAI_L2x2y4z_Px_aa = I_NAI_M3x2y4z_S_aa+ABX*I_NAI_L2x2y4z_S_aa;
  Double I_NAI_L2xy5z_Px_aa = I_NAI_M3xy5z_S_aa+ABX*I_NAI_L2xy5z_S_aa;
  Double I_NAI_L2x6z_Px_aa = I_NAI_M3x6z_S_aa+ABX*I_NAI_L2x6z_S_aa;
  Double I_NAI_Lx7y_Px_aa = I_NAI_M2x7y_S_aa+ABX*I_NAI_Lx7y_S_aa;
  Double I_NAI_Lx6yz_Px_aa = I_NAI_M2x6yz_S_aa+ABX*I_NAI_Lx6yz_S_aa;
  Double I_NAI_Lx5y2z_Px_aa = I_NAI_M2x5y2z_S_aa+ABX*I_NAI_Lx5y2z_S_aa;
  Double I_NAI_Lx4y3z_Px_aa = I_NAI_M2x4y3z_S_aa+ABX*I_NAI_Lx4y3z_S_aa;
  Double I_NAI_Lx3y4z_Px_aa = I_NAI_M2x3y4z_S_aa+ABX*I_NAI_Lx3y4z_S_aa;
  Double I_NAI_Lx2y5z_Px_aa = I_NAI_M2x2y5z_S_aa+ABX*I_NAI_Lx2y5z_S_aa;
  Double I_NAI_Lxy6z_Px_aa = I_NAI_M2xy6z_S_aa+ABX*I_NAI_Lxy6z_S_aa;
  Double I_NAI_Lx7z_Px_aa = I_NAI_M2x7z_S_aa+ABX*I_NAI_Lx7z_S_aa;
  Double I_NAI_L7yz_Px_aa = I_NAI_Mx7yz_S_aa+ABX*I_NAI_L7yz_S_aa;
  Double I_NAI_L6y2z_Px_aa = I_NAI_Mx6y2z_S_aa+ABX*I_NAI_L6y2z_S_aa;
  Double I_NAI_L5y3z_Px_aa = I_NAI_Mx5y3z_S_aa+ABX*I_NAI_L5y3z_S_aa;
  Double I_NAI_L4y4z_Px_aa = I_NAI_Mx4y4z_S_aa+ABX*I_NAI_L4y4z_S_aa;
  Double I_NAI_L3y5z_Px_aa = I_NAI_Mx3y5z_S_aa+ABX*I_NAI_L3y5z_S_aa;
  Double I_NAI_L2y6z_Px_aa = I_NAI_Mx2y6z_S_aa+ABX*I_NAI_L2y6z_S_aa;
  Double I_NAI_Ly7z_Px_aa = I_NAI_Mxy7z_S_aa+ABX*I_NAI_Ly7z_S_aa;
  Double I_NAI_L7xy_Py_aa = I_NAI_M7x2y_S_aa+ABY*I_NAI_L7xy_S_aa;
  Double I_NAI_L6x2y_Py_aa = I_NAI_M6x3y_S_aa+ABY*I_NAI_L6x2y_S_aa;
  Double I_NAI_L6xyz_Py_aa = I_NAI_M6x2yz_S_aa+ABY*I_NAI_L6xyz_S_aa;
  Double I_NAI_L5x3y_Py_aa = I_NAI_M5x4y_S_aa+ABY*I_NAI_L5x3y_S_aa;
  Double I_NAI_L5x2yz_Py_aa = I_NAI_M5x3yz_S_aa+ABY*I_NAI_L5x2yz_S_aa;
  Double I_NAI_L5xy2z_Py_aa = I_NAI_M5x2y2z_S_aa+ABY*I_NAI_L5xy2z_S_aa;
  Double I_NAI_L4x4y_Py_aa = I_NAI_M4x5y_S_aa+ABY*I_NAI_L4x4y_S_aa;
  Double I_NAI_L4x3yz_Py_aa = I_NAI_M4x4yz_S_aa+ABY*I_NAI_L4x3yz_S_aa;
  Double I_NAI_L4x2y2z_Py_aa = I_NAI_M4x3y2z_S_aa+ABY*I_NAI_L4x2y2z_S_aa;
  Double I_NAI_L4xy3z_Py_aa = I_NAI_M4x2y3z_S_aa+ABY*I_NAI_L4xy3z_S_aa;
  Double I_NAI_L3x5y_Py_aa = I_NAI_M3x6y_S_aa+ABY*I_NAI_L3x5y_S_aa;
  Double I_NAI_L3x4yz_Py_aa = I_NAI_M3x5yz_S_aa+ABY*I_NAI_L3x4yz_S_aa;
  Double I_NAI_L3x3y2z_Py_aa = I_NAI_M3x4y2z_S_aa+ABY*I_NAI_L3x3y2z_S_aa;
  Double I_NAI_L3x2y3z_Py_aa = I_NAI_M3x3y3z_S_aa+ABY*I_NAI_L3x2y3z_S_aa;
  Double I_NAI_L3xy4z_Py_aa = I_NAI_M3x2y4z_S_aa+ABY*I_NAI_L3xy4z_S_aa;
  Double I_NAI_L2x6y_Py_aa = I_NAI_M2x7y_S_aa+ABY*I_NAI_L2x6y_S_aa;
  Double I_NAI_L2x5yz_Py_aa = I_NAI_M2x6yz_S_aa+ABY*I_NAI_L2x5yz_S_aa;
  Double I_NAI_L2x4y2z_Py_aa = I_NAI_M2x5y2z_S_aa+ABY*I_NAI_L2x4y2z_S_aa;
  Double I_NAI_L2x3y3z_Py_aa = I_NAI_M2x4y3z_S_aa+ABY*I_NAI_L2x3y3z_S_aa;
  Double I_NAI_L2x2y4z_Py_aa = I_NAI_M2x3y4z_S_aa+ABY*I_NAI_L2x2y4z_S_aa;
  Double I_NAI_L2xy5z_Py_aa = I_NAI_M2x2y5z_S_aa+ABY*I_NAI_L2xy5z_S_aa;
  Double I_NAI_Lx7y_Py_aa = I_NAI_Mx8y_S_aa+ABY*I_NAI_Lx7y_S_aa;
  Double I_NAI_Lx6yz_Py_aa = I_NAI_Mx7yz_S_aa+ABY*I_NAI_Lx6yz_S_aa;
  Double I_NAI_Lx5y2z_Py_aa = I_NAI_Mx6y2z_S_aa+ABY*I_NAI_Lx5y2z_S_aa;
  Double I_NAI_Lx4y3z_Py_aa = I_NAI_Mx5y3z_S_aa+ABY*I_NAI_Lx4y3z_S_aa;
  Double I_NAI_Lx3y4z_Py_aa = I_NAI_Mx4y4z_S_aa+ABY*I_NAI_Lx3y4z_S_aa;
  Double I_NAI_Lx2y5z_Py_aa = I_NAI_Mx3y5z_S_aa+ABY*I_NAI_Lx2y5z_S_aa;
  Double I_NAI_Lxy6z_Py_aa = I_NAI_Mx2y6z_S_aa+ABY*I_NAI_Lxy6z_S_aa;
  Double I_NAI_L8y_Py_aa = I_NAI_M9y_S_aa+ABY*I_NAI_L8y_S_aa;
  Double I_NAI_L7yz_Py_aa = I_NAI_M8yz_S_aa+ABY*I_NAI_L7yz_S_aa;
  Double I_NAI_L6y2z_Py_aa = I_NAI_M7y2z_S_aa+ABY*I_NAI_L6y2z_S_aa;
  Double I_NAI_L5y3z_Py_aa = I_NAI_M6y3z_S_aa+ABY*I_NAI_L5y3z_S_aa;
  Double I_NAI_L4y4z_Py_aa = I_NAI_M5y4z_S_aa+ABY*I_NAI_L4y4z_S_aa;
  Double I_NAI_L3y5z_Py_aa = I_NAI_M4y5z_S_aa+ABY*I_NAI_L3y5z_S_aa;
  Double I_NAI_L2y6z_Py_aa = I_NAI_M3y6z_S_aa+ABY*I_NAI_L2y6z_S_aa;
  Double I_NAI_Ly7z_Py_aa = I_NAI_M2y7z_S_aa+ABY*I_NAI_Ly7z_S_aa;
  Double I_NAI_L7xz_Pz_aa = I_NAI_M7x2z_S_aa+ABZ*I_NAI_L7xz_S_aa;
  Double I_NAI_L6xyz_Pz_aa = I_NAI_M6xy2z_S_aa+ABZ*I_NAI_L6xyz_S_aa;
  Double I_NAI_L6x2z_Pz_aa = I_NAI_M6x3z_S_aa+ABZ*I_NAI_L6x2z_S_aa;
  Double I_NAI_L5x2yz_Pz_aa = I_NAI_M5x2y2z_S_aa+ABZ*I_NAI_L5x2yz_S_aa;
  Double I_NAI_L5xy2z_Pz_aa = I_NAI_M5xy3z_S_aa+ABZ*I_NAI_L5xy2z_S_aa;
  Double I_NAI_L5x3z_Pz_aa = I_NAI_M5x4z_S_aa+ABZ*I_NAI_L5x3z_S_aa;
  Double I_NAI_L4x3yz_Pz_aa = I_NAI_M4x3y2z_S_aa+ABZ*I_NAI_L4x3yz_S_aa;
  Double I_NAI_L4x2y2z_Pz_aa = I_NAI_M4x2y3z_S_aa+ABZ*I_NAI_L4x2y2z_S_aa;
  Double I_NAI_L4xy3z_Pz_aa = I_NAI_M4xy4z_S_aa+ABZ*I_NAI_L4xy3z_S_aa;
  Double I_NAI_L4x4z_Pz_aa = I_NAI_M4x5z_S_aa+ABZ*I_NAI_L4x4z_S_aa;
  Double I_NAI_L3x4yz_Pz_aa = I_NAI_M3x4y2z_S_aa+ABZ*I_NAI_L3x4yz_S_aa;
  Double I_NAI_L3x3y2z_Pz_aa = I_NAI_M3x3y3z_S_aa+ABZ*I_NAI_L3x3y2z_S_aa;
  Double I_NAI_L3x2y3z_Pz_aa = I_NAI_M3x2y4z_S_aa+ABZ*I_NAI_L3x2y3z_S_aa;
  Double I_NAI_L3xy4z_Pz_aa = I_NAI_M3xy5z_S_aa+ABZ*I_NAI_L3xy4z_S_aa;
  Double I_NAI_L3x5z_Pz_aa = I_NAI_M3x6z_S_aa+ABZ*I_NAI_L3x5z_S_aa;
  Double I_NAI_L2x5yz_Pz_aa = I_NAI_M2x5y2z_S_aa+ABZ*I_NAI_L2x5yz_S_aa;
  Double I_NAI_L2x4y2z_Pz_aa = I_NAI_M2x4y3z_S_aa+ABZ*I_NAI_L2x4y2z_S_aa;
  Double I_NAI_L2x3y3z_Pz_aa = I_NAI_M2x3y4z_S_aa+ABZ*I_NAI_L2x3y3z_S_aa;
  Double I_NAI_L2x2y4z_Pz_aa = I_NAI_M2x2y5z_S_aa+ABZ*I_NAI_L2x2y4z_S_aa;
  Double I_NAI_L2xy5z_Pz_aa = I_NAI_M2xy6z_S_aa+ABZ*I_NAI_L2xy5z_S_aa;
  Double I_NAI_L2x6z_Pz_aa = I_NAI_M2x7z_S_aa+ABZ*I_NAI_L2x6z_S_aa;
  Double I_NAI_Lx6yz_Pz_aa = I_NAI_Mx6y2z_S_aa+ABZ*I_NAI_Lx6yz_S_aa;
  Double I_NAI_Lx5y2z_Pz_aa = I_NAI_Mx5y3z_S_aa+ABZ*I_NAI_Lx5y2z_S_aa;
  Double I_NAI_Lx4y3z_Pz_aa = I_NAI_Mx4y4z_S_aa+ABZ*I_NAI_Lx4y3z_S_aa;
  Double I_NAI_Lx3y4z_Pz_aa = I_NAI_Mx3y5z_S_aa+ABZ*I_NAI_Lx3y4z_S_aa;
  Double I_NAI_Lx2y5z_Pz_aa = I_NAI_Mx2y6z_S_aa+ABZ*I_NAI_Lx2y5z_S_aa;
  Double I_NAI_Lxy6z_Pz_aa = I_NAI_Mxy7z_S_aa+ABZ*I_NAI_Lxy6z_S_aa;
  Double I_NAI_Lx7z_Pz_aa = I_NAI_Mx8z_S_aa+ABZ*I_NAI_Lx7z_S_aa;
  Double I_NAI_L7yz_Pz_aa = I_NAI_M7y2z_S_aa+ABZ*I_NAI_L7yz_S_aa;
  Double I_NAI_L6y2z_Pz_aa = I_NAI_M6y3z_S_aa+ABZ*I_NAI_L6y2z_S_aa;
  Double I_NAI_L5y3z_Pz_aa = I_NAI_M5y4z_S_aa+ABZ*I_NAI_L5y3z_S_aa;
  Double I_NAI_L4y4z_Pz_aa = I_NAI_M4y5z_S_aa+ABZ*I_NAI_L4y4z_S_aa;
  Double I_NAI_L3y5z_Pz_aa = I_NAI_M3y6z_S_aa+ABZ*I_NAI_L3y5z_S_aa;
  Double I_NAI_L2y6z_Pz_aa = I_NAI_M2y7z_S_aa+ABZ*I_NAI_L2y6z_S_aa;
  Double I_NAI_Ly7z_Pz_aa = I_NAI_My8z_S_aa+ABZ*I_NAI_Ly7z_S_aa;
  Double I_NAI_L8z_Pz_aa = I_NAI_M9z_S_aa+ABZ*I_NAI_L8z_S_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_K_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 80 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_P_aa
   * RHS shell quartet name: SQ_NAI_K_P_aa
   ************************************************************/
  Double I_NAI_K7x_D2x_aa = I_NAI_L8x_Px_aa+ABX*I_NAI_K7x_Px_aa;
  Double I_NAI_K6xy_D2x_aa = I_NAI_L7xy_Px_aa+ABX*I_NAI_K6xy_Px_aa;
  Double I_NAI_K6xz_D2x_aa = I_NAI_L7xz_Px_aa+ABX*I_NAI_K6xz_Px_aa;
  Double I_NAI_K5x2y_D2x_aa = I_NAI_L6x2y_Px_aa+ABX*I_NAI_K5x2y_Px_aa;
  Double I_NAI_K5xyz_D2x_aa = I_NAI_L6xyz_Px_aa+ABX*I_NAI_K5xyz_Px_aa;
  Double I_NAI_K5x2z_D2x_aa = I_NAI_L6x2z_Px_aa+ABX*I_NAI_K5x2z_Px_aa;
  Double I_NAI_K4x3y_D2x_aa = I_NAI_L5x3y_Px_aa+ABX*I_NAI_K4x3y_Px_aa;
  Double I_NAI_K4x2yz_D2x_aa = I_NAI_L5x2yz_Px_aa+ABX*I_NAI_K4x2yz_Px_aa;
  Double I_NAI_K4xy2z_D2x_aa = I_NAI_L5xy2z_Px_aa+ABX*I_NAI_K4xy2z_Px_aa;
  Double I_NAI_K4x3z_D2x_aa = I_NAI_L5x3z_Px_aa+ABX*I_NAI_K4x3z_Px_aa;
  Double I_NAI_K3x4y_D2x_aa = I_NAI_L4x4y_Px_aa+ABX*I_NAI_K3x4y_Px_aa;
  Double I_NAI_K3x3yz_D2x_aa = I_NAI_L4x3yz_Px_aa+ABX*I_NAI_K3x3yz_Px_aa;
  Double I_NAI_K3x2y2z_D2x_aa = I_NAI_L4x2y2z_Px_aa+ABX*I_NAI_K3x2y2z_Px_aa;
  Double I_NAI_K3xy3z_D2x_aa = I_NAI_L4xy3z_Px_aa+ABX*I_NAI_K3xy3z_Px_aa;
  Double I_NAI_K3x4z_D2x_aa = I_NAI_L4x4z_Px_aa+ABX*I_NAI_K3x4z_Px_aa;
  Double I_NAI_K2x5y_D2x_aa = I_NAI_L3x5y_Px_aa+ABX*I_NAI_K2x5y_Px_aa;
  Double I_NAI_K2x4yz_D2x_aa = I_NAI_L3x4yz_Px_aa+ABX*I_NAI_K2x4yz_Px_aa;
  Double I_NAI_K2x3y2z_D2x_aa = I_NAI_L3x3y2z_Px_aa+ABX*I_NAI_K2x3y2z_Px_aa;
  Double I_NAI_K2x2y3z_D2x_aa = I_NAI_L3x2y3z_Px_aa+ABX*I_NAI_K2x2y3z_Px_aa;
  Double I_NAI_K2xy4z_D2x_aa = I_NAI_L3xy4z_Px_aa+ABX*I_NAI_K2xy4z_Px_aa;
  Double I_NAI_K2x5z_D2x_aa = I_NAI_L3x5z_Px_aa+ABX*I_NAI_K2x5z_Px_aa;
  Double I_NAI_Kx6y_D2x_aa = I_NAI_L2x6y_Px_aa+ABX*I_NAI_Kx6y_Px_aa;
  Double I_NAI_Kx5yz_D2x_aa = I_NAI_L2x5yz_Px_aa+ABX*I_NAI_Kx5yz_Px_aa;
  Double I_NAI_Kx4y2z_D2x_aa = I_NAI_L2x4y2z_Px_aa+ABX*I_NAI_Kx4y2z_Px_aa;
  Double I_NAI_Kx3y3z_D2x_aa = I_NAI_L2x3y3z_Px_aa+ABX*I_NAI_Kx3y3z_Px_aa;
  Double I_NAI_Kx2y4z_D2x_aa = I_NAI_L2x2y4z_Px_aa+ABX*I_NAI_Kx2y4z_Px_aa;
  Double I_NAI_Kxy5z_D2x_aa = I_NAI_L2xy5z_Px_aa+ABX*I_NAI_Kxy5z_Px_aa;
  Double I_NAI_Kx6z_D2x_aa = I_NAI_L2x6z_Px_aa+ABX*I_NAI_Kx6z_Px_aa;
  Double I_NAI_K7y_D2x_aa = I_NAI_Lx7y_Px_aa+ABX*I_NAI_K7y_Px_aa;
  Double I_NAI_K6yz_D2x_aa = I_NAI_Lx6yz_Px_aa+ABX*I_NAI_K6yz_Px_aa;
  Double I_NAI_K5y2z_D2x_aa = I_NAI_Lx5y2z_Px_aa+ABX*I_NAI_K5y2z_Px_aa;
  Double I_NAI_K4y3z_D2x_aa = I_NAI_Lx4y3z_Px_aa+ABX*I_NAI_K4y3z_Px_aa;
  Double I_NAI_K3y4z_D2x_aa = I_NAI_Lx3y4z_Px_aa+ABX*I_NAI_K3y4z_Px_aa;
  Double I_NAI_K2y5z_D2x_aa = I_NAI_Lx2y5z_Px_aa+ABX*I_NAI_K2y5z_Px_aa;
  Double I_NAI_Ky6z_D2x_aa = I_NAI_Lxy6z_Px_aa+ABX*I_NAI_Ky6z_Px_aa;
  Double I_NAI_K7z_D2x_aa = I_NAI_Lx7z_Px_aa+ABX*I_NAI_K7z_Px_aa;
  Double I_NAI_K6xz_Dxy_aa = I_NAI_L6xyz_Px_aa+ABY*I_NAI_K6xz_Px_aa;
  Double I_NAI_K5xyz_Dxy_aa = I_NAI_L5x2yz_Px_aa+ABY*I_NAI_K5xyz_Px_aa;
  Double I_NAI_K5x2z_Dxy_aa = I_NAI_L5xy2z_Px_aa+ABY*I_NAI_K5x2z_Px_aa;
  Double I_NAI_K4x2yz_Dxy_aa = I_NAI_L4x3yz_Px_aa+ABY*I_NAI_K4x2yz_Px_aa;
  Double I_NAI_K4xy2z_Dxy_aa = I_NAI_L4x2y2z_Px_aa+ABY*I_NAI_K4xy2z_Px_aa;
  Double I_NAI_K4x3z_Dxy_aa = I_NAI_L4xy3z_Px_aa+ABY*I_NAI_K4x3z_Px_aa;
  Double I_NAI_K3x3yz_Dxy_aa = I_NAI_L3x4yz_Px_aa+ABY*I_NAI_K3x3yz_Px_aa;
  Double I_NAI_K3x2y2z_Dxy_aa = I_NAI_L3x3y2z_Px_aa+ABY*I_NAI_K3x2y2z_Px_aa;
  Double I_NAI_K3xy3z_Dxy_aa = I_NAI_L3x2y3z_Px_aa+ABY*I_NAI_K3xy3z_Px_aa;
  Double I_NAI_K3x4z_Dxy_aa = I_NAI_L3xy4z_Px_aa+ABY*I_NAI_K3x4z_Px_aa;
  Double I_NAI_K2x4yz_Dxy_aa = I_NAI_L2x5yz_Px_aa+ABY*I_NAI_K2x4yz_Px_aa;
  Double I_NAI_K2x3y2z_Dxy_aa = I_NAI_L2x4y2z_Px_aa+ABY*I_NAI_K2x3y2z_Px_aa;
  Double I_NAI_K2x2y3z_Dxy_aa = I_NAI_L2x3y3z_Px_aa+ABY*I_NAI_K2x2y3z_Px_aa;
  Double I_NAI_K2xy4z_Dxy_aa = I_NAI_L2x2y4z_Px_aa+ABY*I_NAI_K2xy4z_Px_aa;
  Double I_NAI_K2x5z_Dxy_aa = I_NAI_L2xy5z_Px_aa+ABY*I_NAI_K2x5z_Px_aa;
  Double I_NAI_Kx5yz_Dxy_aa = I_NAI_Lx6yz_Px_aa+ABY*I_NAI_Kx5yz_Px_aa;
  Double I_NAI_Kx4y2z_Dxy_aa = I_NAI_Lx5y2z_Px_aa+ABY*I_NAI_Kx4y2z_Px_aa;
  Double I_NAI_Kx3y3z_Dxy_aa = I_NAI_Lx4y3z_Px_aa+ABY*I_NAI_Kx3y3z_Px_aa;
  Double I_NAI_Kx2y4z_Dxy_aa = I_NAI_Lx3y4z_Px_aa+ABY*I_NAI_Kx2y4z_Px_aa;
  Double I_NAI_Kxy5z_Dxy_aa = I_NAI_Lx2y5z_Px_aa+ABY*I_NAI_Kxy5z_Px_aa;
  Double I_NAI_Kx6z_Dxy_aa = I_NAI_Lxy6z_Px_aa+ABY*I_NAI_Kx6z_Px_aa;
  Double I_NAI_K6yz_Dxy_aa = I_NAI_L7yz_Px_aa+ABY*I_NAI_K6yz_Px_aa;
  Double I_NAI_K5y2z_Dxy_aa = I_NAI_L6y2z_Px_aa+ABY*I_NAI_K5y2z_Px_aa;
  Double I_NAI_K4y3z_Dxy_aa = I_NAI_L5y3z_Px_aa+ABY*I_NAI_K4y3z_Px_aa;
  Double I_NAI_K3y4z_Dxy_aa = I_NAI_L4y4z_Px_aa+ABY*I_NAI_K3y4z_Px_aa;
  Double I_NAI_K2y5z_Dxy_aa = I_NAI_L3y5z_Px_aa+ABY*I_NAI_K2y5z_Px_aa;
  Double I_NAI_Ky6z_Dxy_aa = I_NAI_L2y6z_Px_aa+ABY*I_NAI_Ky6z_Px_aa;
  Double I_NAI_K7z_Dxy_aa = I_NAI_Ly7z_Px_aa+ABY*I_NAI_K7z_Px_aa;
  Double I_NAI_K7x_D2y_aa = I_NAI_L7xy_Py_aa+ABY*I_NAI_K7x_Py_aa;
  Double I_NAI_K6xy_D2y_aa = I_NAI_L6x2y_Py_aa+ABY*I_NAI_K6xy_Py_aa;
  Double I_NAI_K6xz_D2y_aa = I_NAI_L6xyz_Py_aa+ABY*I_NAI_K6xz_Py_aa;
  Double I_NAI_K5x2y_D2y_aa = I_NAI_L5x3y_Py_aa+ABY*I_NAI_K5x2y_Py_aa;
  Double I_NAI_K5xyz_D2y_aa = I_NAI_L5x2yz_Py_aa+ABY*I_NAI_K5xyz_Py_aa;
  Double I_NAI_K5x2z_D2y_aa = I_NAI_L5xy2z_Py_aa+ABY*I_NAI_K5x2z_Py_aa;
  Double I_NAI_K4x3y_D2y_aa = I_NAI_L4x4y_Py_aa+ABY*I_NAI_K4x3y_Py_aa;
  Double I_NAI_K4x2yz_D2y_aa = I_NAI_L4x3yz_Py_aa+ABY*I_NAI_K4x2yz_Py_aa;
  Double I_NAI_K4xy2z_D2y_aa = I_NAI_L4x2y2z_Py_aa+ABY*I_NAI_K4xy2z_Py_aa;
  Double I_NAI_K4x3z_D2y_aa = I_NAI_L4xy3z_Py_aa+ABY*I_NAI_K4x3z_Py_aa;
  Double I_NAI_K3x4y_D2y_aa = I_NAI_L3x5y_Py_aa+ABY*I_NAI_K3x4y_Py_aa;
  Double I_NAI_K3x3yz_D2y_aa = I_NAI_L3x4yz_Py_aa+ABY*I_NAI_K3x3yz_Py_aa;
  Double I_NAI_K3x2y2z_D2y_aa = I_NAI_L3x3y2z_Py_aa+ABY*I_NAI_K3x2y2z_Py_aa;
  Double I_NAI_K3xy3z_D2y_aa = I_NAI_L3x2y3z_Py_aa+ABY*I_NAI_K3xy3z_Py_aa;
  Double I_NAI_K3x4z_D2y_aa = I_NAI_L3xy4z_Py_aa+ABY*I_NAI_K3x4z_Py_aa;
  Double I_NAI_K2x5y_D2y_aa = I_NAI_L2x6y_Py_aa+ABY*I_NAI_K2x5y_Py_aa;
  Double I_NAI_K2x4yz_D2y_aa = I_NAI_L2x5yz_Py_aa+ABY*I_NAI_K2x4yz_Py_aa;
  Double I_NAI_K2x3y2z_D2y_aa = I_NAI_L2x4y2z_Py_aa+ABY*I_NAI_K2x3y2z_Py_aa;
  Double I_NAI_K2x2y3z_D2y_aa = I_NAI_L2x3y3z_Py_aa+ABY*I_NAI_K2x2y3z_Py_aa;
  Double I_NAI_K2xy4z_D2y_aa = I_NAI_L2x2y4z_Py_aa+ABY*I_NAI_K2xy4z_Py_aa;
  Double I_NAI_K2x5z_D2y_aa = I_NAI_L2xy5z_Py_aa+ABY*I_NAI_K2x5z_Py_aa;
  Double I_NAI_Kx6y_D2y_aa = I_NAI_Lx7y_Py_aa+ABY*I_NAI_Kx6y_Py_aa;
  Double I_NAI_Kx5yz_D2y_aa = I_NAI_Lx6yz_Py_aa+ABY*I_NAI_Kx5yz_Py_aa;
  Double I_NAI_Kx4y2z_D2y_aa = I_NAI_Lx5y2z_Py_aa+ABY*I_NAI_Kx4y2z_Py_aa;
  Double I_NAI_Kx3y3z_D2y_aa = I_NAI_Lx4y3z_Py_aa+ABY*I_NAI_Kx3y3z_Py_aa;
  Double I_NAI_Kx2y4z_D2y_aa = I_NAI_Lx3y4z_Py_aa+ABY*I_NAI_Kx2y4z_Py_aa;
  Double I_NAI_Kxy5z_D2y_aa = I_NAI_Lx2y5z_Py_aa+ABY*I_NAI_Kxy5z_Py_aa;
  Double I_NAI_Kx6z_D2y_aa = I_NAI_Lxy6z_Py_aa+ABY*I_NAI_Kx6z_Py_aa;
  Double I_NAI_K7y_D2y_aa = I_NAI_L8y_Py_aa+ABY*I_NAI_K7y_Py_aa;
  Double I_NAI_K6yz_D2y_aa = I_NAI_L7yz_Py_aa+ABY*I_NAI_K6yz_Py_aa;
  Double I_NAI_K5y2z_D2y_aa = I_NAI_L6y2z_Py_aa+ABY*I_NAI_K5y2z_Py_aa;
  Double I_NAI_K4y3z_D2y_aa = I_NAI_L5y3z_Py_aa+ABY*I_NAI_K4y3z_Py_aa;
  Double I_NAI_K3y4z_D2y_aa = I_NAI_L4y4z_Py_aa+ABY*I_NAI_K3y4z_Py_aa;
  Double I_NAI_K2y5z_D2y_aa = I_NAI_L3y5z_Py_aa+ABY*I_NAI_K2y5z_Py_aa;
  Double I_NAI_Ky6z_D2y_aa = I_NAI_L2y6z_Py_aa+ABY*I_NAI_Ky6z_Py_aa;
  Double I_NAI_K7z_D2y_aa = I_NAI_Ly7z_Py_aa+ABY*I_NAI_K7z_Py_aa;
  Double I_NAI_K7x_D2z_aa = I_NAI_L7xz_Pz_aa+ABZ*I_NAI_K7x_Pz_aa;
  Double I_NAI_K6xy_D2z_aa = I_NAI_L6xyz_Pz_aa+ABZ*I_NAI_K6xy_Pz_aa;
  Double I_NAI_K6xz_D2z_aa = I_NAI_L6x2z_Pz_aa+ABZ*I_NAI_K6xz_Pz_aa;
  Double I_NAI_K5x2y_D2z_aa = I_NAI_L5x2yz_Pz_aa+ABZ*I_NAI_K5x2y_Pz_aa;
  Double I_NAI_K5xyz_D2z_aa = I_NAI_L5xy2z_Pz_aa+ABZ*I_NAI_K5xyz_Pz_aa;
  Double I_NAI_K5x2z_D2z_aa = I_NAI_L5x3z_Pz_aa+ABZ*I_NAI_K5x2z_Pz_aa;
  Double I_NAI_K4x3y_D2z_aa = I_NAI_L4x3yz_Pz_aa+ABZ*I_NAI_K4x3y_Pz_aa;
  Double I_NAI_K4x2yz_D2z_aa = I_NAI_L4x2y2z_Pz_aa+ABZ*I_NAI_K4x2yz_Pz_aa;
  Double I_NAI_K4xy2z_D2z_aa = I_NAI_L4xy3z_Pz_aa+ABZ*I_NAI_K4xy2z_Pz_aa;
  Double I_NAI_K4x3z_D2z_aa = I_NAI_L4x4z_Pz_aa+ABZ*I_NAI_K4x3z_Pz_aa;
  Double I_NAI_K3x4y_D2z_aa = I_NAI_L3x4yz_Pz_aa+ABZ*I_NAI_K3x4y_Pz_aa;
  Double I_NAI_K3x3yz_D2z_aa = I_NAI_L3x3y2z_Pz_aa+ABZ*I_NAI_K3x3yz_Pz_aa;
  Double I_NAI_K3x2y2z_D2z_aa = I_NAI_L3x2y3z_Pz_aa+ABZ*I_NAI_K3x2y2z_Pz_aa;
  Double I_NAI_K3xy3z_D2z_aa = I_NAI_L3xy4z_Pz_aa+ABZ*I_NAI_K3xy3z_Pz_aa;
  Double I_NAI_K3x4z_D2z_aa = I_NAI_L3x5z_Pz_aa+ABZ*I_NAI_K3x4z_Pz_aa;
  Double I_NAI_K2x5y_D2z_aa = I_NAI_L2x5yz_Pz_aa+ABZ*I_NAI_K2x5y_Pz_aa;
  Double I_NAI_K2x4yz_D2z_aa = I_NAI_L2x4y2z_Pz_aa+ABZ*I_NAI_K2x4yz_Pz_aa;
  Double I_NAI_K2x3y2z_D2z_aa = I_NAI_L2x3y3z_Pz_aa+ABZ*I_NAI_K2x3y2z_Pz_aa;
  Double I_NAI_K2x2y3z_D2z_aa = I_NAI_L2x2y4z_Pz_aa+ABZ*I_NAI_K2x2y3z_Pz_aa;
  Double I_NAI_K2xy4z_D2z_aa = I_NAI_L2xy5z_Pz_aa+ABZ*I_NAI_K2xy4z_Pz_aa;
  Double I_NAI_K2x5z_D2z_aa = I_NAI_L2x6z_Pz_aa+ABZ*I_NAI_K2x5z_Pz_aa;
  Double I_NAI_Kx6y_D2z_aa = I_NAI_Lx6yz_Pz_aa+ABZ*I_NAI_Kx6y_Pz_aa;
  Double I_NAI_Kx5yz_D2z_aa = I_NAI_Lx5y2z_Pz_aa+ABZ*I_NAI_Kx5yz_Pz_aa;
  Double I_NAI_Kx4y2z_D2z_aa = I_NAI_Lx4y3z_Pz_aa+ABZ*I_NAI_Kx4y2z_Pz_aa;
  Double I_NAI_Kx3y3z_D2z_aa = I_NAI_Lx3y4z_Pz_aa+ABZ*I_NAI_Kx3y3z_Pz_aa;
  Double I_NAI_Kx2y4z_D2z_aa = I_NAI_Lx2y5z_Pz_aa+ABZ*I_NAI_Kx2y4z_Pz_aa;
  Double I_NAI_Kxy5z_D2z_aa = I_NAI_Lxy6z_Pz_aa+ABZ*I_NAI_Kxy5z_Pz_aa;
  Double I_NAI_Kx6z_D2z_aa = I_NAI_Lx7z_Pz_aa+ABZ*I_NAI_Kx6z_Pz_aa;
  Double I_NAI_K7y_D2z_aa = I_NAI_L7yz_Pz_aa+ABZ*I_NAI_K7y_Pz_aa;
  Double I_NAI_K6yz_D2z_aa = I_NAI_L6y2z_Pz_aa+ABZ*I_NAI_K6yz_Pz_aa;
  Double I_NAI_K5y2z_D2z_aa = I_NAI_L5y3z_Pz_aa+ABZ*I_NAI_K5y2z_Pz_aa;
  Double I_NAI_K4y3z_D2z_aa = I_NAI_L4y4z_Pz_aa+ABZ*I_NAI_K4y3z_Pz_aa;
  Double I_NAI_K3y4z_D2z_aa = I_NAI_L3y5z_Pz_aa+ABZ*I_NAI_K3y4z_Pz_aa;
  Double I_NAI_K2y5z_D2z_aa = I_NAI_L2y6z_Pz_aa+ABZ*I_NAI_K2y5z_Pz_aa;
  Double I_NAI_Ky6z_D2z_aa = I_NAI_Ly7z_Pz_aa+ABZ*I_NAI_Ky6z_Pz_aa;
  Double I_NAI_K7z_D2z_aa = I_NAI_L8z_Pz_aa+ABZ*I_NAI_K7z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_I_F_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_D_aa
   * RHS shell quartet name: SQ_NAI_I_D_aa
   ************************************************************/
  Double I_NAI_I6x_F3x_aa = I_NAI_K7x_D2x_aa+ABX*I_NAI_I6x_D2x_aa;
  Double I_NAI_I5xy_F3x_aa = I_NAI_K6xy_D2x_aa+ABX*I_NAI_I5xy_D2x_aa;
  Double I_NAI_I5xz_F3x_aa = I_NAI_K6xz_D2x_aa+ABX*I_NAI_I5xz_D2x_aa;
  Double I_NAI_I4x2y_F3x_aa = I_NAI_K5x2y_D2x_aa+ABX*I_NAI_I4x2y_D2x_aa;
  Double I_NAI_I4xyz_F3x_aa = I_NAI_K5xyz_D2x_aa+ABX*I_NAI_I4xyz_D2x_aa;
  Double I_NAI_I4x2z_F3x_aa = I_NAI_K5x2z_D2x_aa+ABX*I_NAI_I4x2z_D2x_aa;
  Double I_NAI_I3x3y_F3x_aa = I_NAI_K4x3y_D2x_aa+ABX*I_NAI_I3x3y_D2x_aa;
  Double I_NAI_I3x2yz_F3x_aa = I_NAI_K4x2yz_D2x_aa+ABX*I_NAI_I3x2yz_D2x_aa;
  Double I_NAI_I3xy2z_F3x_aa = I_NAI_K4xy2z_D2x_aa+ABX*I_NAI_I3xy2z_D2x_aa;
  Double I_NAI_I3x3z_F3x_aa = I_NAI_K4x3z_D2x_aa+ABX*I_NAI_I3x3z_D2x_aa;
  Double I_NAI_I2x4y_F3x_aa = I_NAI_K3x4y_D2x_aa+ABX*I_NAI_I2x4y_D2x_aa;
  Double I_NAI_I2x3yz_F3x_aa = I_NAI_K3x3yz_D2x_aa+ABX*I_NAI_I2x3yz_D2x_aa;
  Double I_NAI_I2x2y2z_F3x_aa = I_NAI_K3x2y2z_D2x_aa+ABX*I_NAI_I2x2y2z_D2x_aa;
  Double I_NAI_I2xy3z_F3x_aa = I_NAI_K3xy3z_D2x_aa+ABX*I_NAI_I2xy3z_D2x_aa;
  Double I_NAI_I2x4z_F3x_aa = I_NAI_K3x4z_D2x_aa+ABX*I_NAI_I2x4z_D2x_aa;
  Double I_NAI_Ix5y_F3x_aa = I_NAI_K2x5y_D2x_aa+ABX*I_NAI_Ix5y_D2x_aa;
  Double I_NAI_Ix4yz_F3x_aa = I_NAI_K2x4yz_D2x_aa+ABX*I_NAI_Ix4yz_D2x_aa;
  Double I_NAI_Ix3y2z_F3x_aa = I_NAI_K2x3y2z_D2x_aa+ABX*I_NAI_Ix3y2z_D2x_aa;
  Double I_NAI_Ix2y3z_F3x_aa = I_NAI_K2x2y3z_D2x_aa+ABX*I_NAI_Ix2y3z_D2x_aa;
  Double I_NAI_Ixy4z_F3x_aa = I_NAI_K2xy4z_D2x_aa+ABX*I_NAI_Ixy4z_D2x_aa;
  Double I_NAI_Ix5z_F3x_aa = I_NAI_K2x5z_D2x_aa+ABX*I_NAI_Ix5z_D2x_aa;
  Double I_NAI_I6y_F3x_aa = I_NAI_Kx6y_D2x_aa+ABX*I_NAI_I6y_D2x_aa;
  Double I_NAI_I5yz_F3x_aa = I_NAI_Kx5yz_D2x_aa+ABX*I_NAI_I5yz_D2x_aa;
  Double I_NAI_I4y2z_F3x_aa = I_NAI_Kx4y2z_D2x_aa+ABX*I_NAI_I4y2z_D2x_aa;
  Double I_NAI_I3y3z_F3x_aa = I_NAI_Kx3y3z_D2x_aa+ABX*I_NAI_I3y3z_D2x_aa;
  Double I_NAI_I2y4z_F3x_aa = I_NAI_Kx2y4z_D2x_aa+ABX*I_NAI_I2y4z_D2x_aa;
  Double I_NAI_Iy5z_F3x_aa = I_NAI_Kxy5z_D2x_aa+ABX*I_NAI_Iy5z_D2x_aa;
  Double I_NAI_I6z_F3x_aa = I_NAI_Kx6z_D2x_aa+ABX*I_NAI_I6z_D2x_aa;
  Double I_NAI_I6x_F2xy_aa = I_NAI_K6xy_D2x_aa+ABY*I_NAI_I6x_D2x_aa;
  Double I_NAI_I5xy_F2xy_aa = I_NAI_K5x2y_D2x_aa+ABY*I_NAI_I5xy_D2x_aa;
  Double I_NAI_I5xz_F2xy_aa = I_NAI_K5xyz_D2x_aa+ABY*I_NAI_I5xz_D2x_aa;
  Double I_NAI_I4x2y_F2xy_aa = I_NAI_K4x3y_D2x_aa+ABY*I_NAI_I4x2y_D2x_aa;
  Double I_NAI_I4xyz_F2xy_aa = I_NAI_K4x2yz_D2x_aa+ABY*I_NAI_I4xyz_D2x_aa;
  Double I_NAI_I4x2z_F2xy_aa = I_NAI_K4xy2z_D2x_aa+ABY*I_NAI_I4x2z_D2x_aa;
  Double I_NAI_I3x3y_F2xy_aa = I_NAI_K3x4y_D2x_aa+ABY*I_NAI_I3x3y_D2x_aa;
  Double I_NAI_I3x2yz_F2xy_aa = I_NAI_K3x3yz_D2x_aa+ABY*I_NAI_I3x2yz_D2x_aa;
  Double I_NAI_I3xy2z_F2xy_aa = I_NAI_K3x2y2z_D2x_aa+ABY*I_NAI_I3xy2z_D2x_aa;
  Double I_NAI_I3x3z_F2xy_aa = I_NAI_K3xy3z_D2x_aa+ABY*I_NAI_I3x3z_D2x_aa;
  Double I_NAI_I2x4y_F2xy_aa = I_NAI_K2x5y_D2x_aa+ABY*I_NAI_I2x4y_D2x_aa;
  Double I_NAI_I2x3yz_F2xy_aa = I_NAI_K2x4yz_D2x_aa+ABY*I_NAI_I2x3yz_D2x_aa;
  Double I_NAI_I2x2y2z_F2xy_aa = I_NAI_K2x3y2z_D2x_aa+ABY*I_NAI_I2x2y2z_D2x_aa;
  Double I_NAI_I2xy3z_F2xy_aa = I_NAI_K2x2y3z_D2x_aa+ABY*I_NAI_I2xy3z_D2x_aa;
  Double I_NAI_I2x4z_F2xy_aa = I_NAI_K2xy4z_D2x_aa+ABY*I_NAI_I2x4z_D2x_aa;
  Double I_NAI_Ix5y_F2xy_aa = I_NAI_Kx6y_D2x_aa+ABY*I_NAI_Ix5y_D2x_aa;
  Double I_NAI_Ix4yz_F2xy_aa = I_NAI_Kx5yz_D2x_aa+ABY*I_NAI_Ix4yz_D2x_aa;
  Double I_NAI_Ix3y2z_F2xy_aa = I_NAI_Kx4y2z_D2x_aa+ABY*I_NAI_Ix3y2z_D2x_aa;
  Double I_NAI_Ix2y3z_F2xy_aa = I_NAI_Kx3y3z_D2x_aa+ABY*I_NAI_Ix2y3z_D2x_aa;
  Double I_NAI_Ixy4z_F2xy_aa = I_NAI_Kx2y4z_D2x_aa+ABY*I_NAI_Ixy4z_D2x_aa;
  Double I_NAI_Ix5z_F2xy_aa = I_NAI_Kxy5z_D2x_aa+ABY*I_NAI_Ix5z_D2x_aa;
  Double I_NAI_I6y_F2xy_aa = I_NAI_K7y_D2x_aa+ABY*I_NAI_I6y_D2x_aa;
  Double I_NAI_I5yz_F2xy_aa = I_NAI_K6yz_D2x_aa+ABY*I_NAI_I5yz_D2x_aa;
  Double I_NAI_I4y2z_F2xy_aa = I_NAI_K5y2z_D2x_aa+ABY*I_NAI_I4y2z_D2x_aa;
  Double I_NAI_I3y3z_F2xy_aa = I_NAI_K4y3z_D2x_aa+ABY*I_NAI_I3y3z_D2x_aa;
  Double I_NAI_I2y4z_F2xy_aa = I_NAI_K3y4z_D2x_aa+ABY*I_NAI_I2y4z_D2x_aa;
  Double I_NAI_Iy5z_F2xy_aa = I_NAI_K2y5z_D2x_aa+ABY*I_NAI_Iy5z_D2x_aa;
  Double I_NAI_I6z_F2xy_aa = I_NAI_Ky6z_D2x_aa+ABY*I_NAI_I6z_D2x_aa;
  Double I_NAI_I6x_F2xz_aa = I_NAI_K6xz_D2x_aa+ABZ*I_NAI_I6x_D2x_aa;
  Double I_NAI_I5xy_F2xz_aa = I_NAI_K5xyz_D2x_aa+ABZ*I_NAI_I5xy_D2x_aa;
  Double I_NAI_I5xz_F2xz_aa = I_NAI_K5x2z_D2x_aa+ABZ*I_NAI_I5xz_D2x_aa;
  Double I_NAI_I4x2y_F2xz_aa = I_NAI_K4x2yz_D2x_aa+ABZ*I_NAI_I4x2y_D2x_aa;
  Double I_NAI_I4xyz_F2xz_aa = I_NAI_K4xy2z_D2x_aa+ABZ*I_NAI_I4xyz_D2x_aa;
  Double I_NAI_I4x2z_F2xz_aa = I_NAI_K4x3z_D2x_aa+ABZ*I_NAI_I4x2z_D2x_aa;
  Double I_NAI_I3x3y_F2xz_aa = I_NAI_K3x3yz_D2x_aa+ABZ*I_NAI_I3x3y_D2x_aa;
  Double I_NAI_I3x2yz_F2xz_aa = I_NAI_K3x2y2z_D2x_aa+ABZ*I_NAI_I3x2yz_D2x_aa;
  Double I_NAI_I3xy2z_F2xz_aa = I_NAI_K3xy3z_D2x_aa+ABZ*I_NAI_I3xy2z_D2x_aa;
  Double I_NAI_I3x3z_F2xz_aa = I_NAI_K3x4z_D2x_aa+ABZ*I_NAI_I3x3z_D2x_aa;
  Double I_NAI_I2x4y_F2xz_aa = I_NAI_K2x4yz_D2x_aa+ABZ*I_NAI_I2x4y_D2x_aa;
  Double I_NAI_I2x3yz_F2xz_aa = I_NAI_K2x3y2z_D2x_aa+ABZ*I_NAI_I2x3yz_D2x_aa;
  Double I_NAI_I2x2y2z_F2xz_aa = I_NAI_K2x2y3z_D2x_aa+ABZ*I_NAI_I2x2y2z_D2x_aa;
  Double I_NAI_I2xy3z_F2xz_aa = I_NAI_K2xy4z_D2x_aa+ABZ*I_NAI_I2xy3z_D2x_aa;
  Double I_NAI_I2x4z_F2xz_aa = I_NAI_K2x5z_D2x_aa+ABZ*I_NAI_I2x4z_D2x_aa;
  Double I_NAI_Ix5y_F2xz_aa = I_NAI_Kx5yz_D2x_aa+ABZ*I_NAI_Ix5y_D2x_aa;
  Double I_NAI_Ix4yz_F2xz_aa = I_NAI_Kx4y2z_D2x_aa+ABZ*I_NAI_Ix4yz_D2x_aa;
  Double I_NAI_Ix3y2z_F2xz_aa = I_NAI_Kx3y3z_D2x_aa+ABZ*I_NAI_Ix3y2z_D2x_aa;
  Double I_NAI_Ix2y3z_F2xz_aa = I_NAI_Kx2y4z_D2x_aa+ABZ*I_NAI_Ix2y3z_D2x_aa;
  Double I_NAI_Ixy4z_F2xz_aa = I_NAI_Kxy5z_D2x_aa+ABZ*I_NAI_Ixy4z_D2x_aa;
  Double I_NAI_Ix5z_F2xz_aa = I_NAI_Kx6z_D2x_aa+ABZ*I_NAI_Ix5z_D2x_aa;
  Double I_NAI_I6y_F2xz_aa = I_NAI_K6yz_D2x_aa+ABZ*I_NAI_I6y_D2x_aa;
  Double I_NAI_I5yz_F2xz_aa = I_NAI_K5y2z_D2x_aa+ABZ*I_NAI_I5yz_D2x_aa;
  Double I_NAI_I4y2z_F2xz_aa = I_NAI_K4y3z_D2x_aa+ABZ*I_NAI_I4y2z_D2x_aa;
  Double I_NAI_I3y3z_F2xz_aa = I_NAI_K3y4z_D2x_aa+ABZ*I_NAI_I3y3z_D2x_aa;
  Double I_NAI_I2y4z_F2xz_aa = I_NAI_K2y5z_D2x_aa+ABZ*I_NAI_I2y4z_D2x_aa;
  Double I_NAI_Iy5z_F2xz_aa = I_NAI_Ky6z_D2x_aa+ABZ*I_NAI_Iy5z_D2x_aa;
  Double I_NAI_I6z_F2xz_aa = I_NAI_K7z_D2x_aa+ABZ*I_NAI_I6z_D2x_aa;
  Double I_NAI_I6x_Fx2y_aa = I_NAI_K7x_D2y_aa+ABX*I_NAI_I6x_D2y_aa;
  Double I_NAI_I5xy_Fx2y_aa = I_NAI_K6xy_D2y_aa+ABX*I_NAI_I5xy_D2y_aa;
  Double I_NAI_I5xz_Fx2y_aa = I_NAI_K6xz_D2y_aa+ABX*I_NAI_I5xz_D2y_aa;
  Double I_NAI_I4x2y_Fx2y_aa = I_NAI_K5x2y_D2y_aa+ABX*I_NAI_I4x2y_D2y_aa;
  Double I_NAI_I4xyz_Fx2y_aa = I_NAI_K5xyz_D2y_aa+ABX*I_NAI_I4xyz_D2y_aa;
  Double I_NAI_I4x2z_Fx2y_aa = I_NAI_K5x2z_D2y_aa+ABX*I_NAI_I4x2z_D2y_aa;
  Double I_NAI_I3x3y_Fx2y_aa = I_NAI_K4x3y_D2y_aa+ABX*I_NAI_I3x3y_D2y_aa;
  Double I_NAI_I3x2yz_Fx2y_aa = I_NAI_K4x2yz_D2y_aa+ABX*I_NAI_I3x2yz_D2y_aa;
  Double I_NAI_I3xy2z_Fx2y_aa = I_NAI_K4xy2z_D2y_aa+ABX*I_NAI_I3xy2z_D2y_aa;
  Double I_NAI_I3x3z_Fx2y_aa = I_NAI_K4x3z_D2y_aa+ABX*I_NAI_I3x3z_D2y_aa;
  Double I_NAI_I2x4y_Fx2y_aa = I_NAI_K3x4y_D2y_aa+ABX*I_NAI_I2x4y_D2y_aa;
  Double I_NAI_I2x3yz_Fx2y_aa = I_NAI_K3x3yz_D2y_aa+ABX*I_NAI_I2x3yz_D2y_aa;
  Double I_NAI_I2x2y2z_Fx2y_aa = I_NAI_K3x2y2z_D2y_aa+ABX*I_NAI_I2x2y2z_D2y_aa;
  Double I_NAI_I2xy3z_Fx2y_aa = I_NAI_K3xy3z_D2y_aa+ABX*I_NAI_I2xy3z_D2y_aa;
  Double I_NAI_I2x4z_Fx2y_aa = I_NAI_K3x4z_D2y_aa+ABX*I_NAI_I2x4z_D2y_aa;
  Double I_NAI_Ix5y_Fx2y_aa = I_NAI_K2x5y_D2y_aa+ABX*I_NAI_Ix5y_D2y_aa;
  Double I_NAI_Ix4yz_Fx2y_aa = I_NAI_K2x4yz_D2y_aa+ABX*I_NAI_Ix4yz_D2y_aa;
  Double I_NAI_Ix3y2z_Fx2y_aa = I_NAI_K2x3y2z_D2y_aa+ABX*I_NAI_Ix3y2z_D2y_aa;
  Double I_NAI_Ix2y3z_Fx2y_aa = I_NAI_K2x2y3z_D2y_aa+ABX*I_NAI_Ix2y3z_D2y_aa;
  Double I_NAI_Ixy4z_Fx2y_aa = I_NAI_K2xy4z_D2y_aa+ABX*I_NAI_Ixy4z_D2y_aa;
  Double I_NAI_Ix5z_Fx2y_aa = I_NAI_K2x5z_D2y_aa+ABX*I_NAI_Ix5z_D2y_aa;
  Double I_NAI_I6y_Fx2y_aa = I_NAI_Kx6y_D2y_aa+ABX*I_NAI_I6y_D2y_aa;
  Double I_NAI_I5yz_Fx2y_aa = I_NAI_Kx5yz_D2y_aa+ABX*I_NAI_I5yz_D2y_aa;
  Double I_NAI_I4y2z_Fx2y_aa = I_NAI_Kx4y2z_D2y_aa+ABX*I_NAI_I4y2z_D2y_aa;
  Double I_NAI_I3y3z_Fx2y_aa = I_NAI_Kx3y3z_D2y_aa+ABX*I_NAI_I3y3z_D2y_aa;
  Double I_NAI_I2y4z_Fx2y_aa = I_NAI_Kx2y4z_D2y_aa+ABX*I_NAI_I2y4z_D2y_aa;
  Double I_NAI_Iy5z_Fx2y_aa = I_NAI_Kxy5z_D2y_aa+ABX*I_NAI_Iy5z_D2y_aa;
  Double I_NAI_I6z_Fx2y_aa = I_NAI_Kx6z_D2y_aa+ABX*I_NAI_I6z_D2y_aa;
  Double I_NAI_I6x_Fxyz_aa = I_NAI_K6xz_Dxy_aa+ABZ*I_NAI_I6x_Dxy_aa;
  Double I_NAI_I5xy_Fxyz_aa = I_NAI_K5xyz_Dxy_aa+ABZ*I_NAI_I5xy_Dxy_aa;
  Double I_NAI_I5xz_Fxyz_aa = I_NAI_K5x2z_Dxy_aa+ABZ*I_NAI_I5xz_Dxy_aa;
  Double I_NAI_I4x2y_Fxyz_aa = I_NAI_K4x2yz_Dxy_aa+ABZ*I_NAI_I4x2y_Dxy_aa;
  Double I_NAI_I4xyz_Fxyz_aa = I_NAI_K4xy2z_Dxy_aa+ABZ*I_NAI_I4xyz_Dxy_aa;
  Double I_NAI_I4x2z_Fxyz_aa = I_NAI_K4x3z_Dxy_aa+ABZ*I_NAI_I4x2z_Dxy_aa;
  Double I_NAI_I3x3y_Fxyz_aa = I_NAI_K3x3yz_Dxy_aa+ABZ*I_NAI_I3x3y_Dxy_aa;
  Double I_NAI_I3x2yz_Fxyz_aa = I_NAI_K3x2y2z_Dxy_aa+ABZ*I_NAI_I3x2yz_Dxy_aa;
  Double I_NAI_I3xy2z_Fxyz_aa = I_NAI_K3xy3z_Dxy_aa+ABZ*I_NAI_I3xy2z_Dxy_aa;
  Double I_NAI_I3x3z_Fxyz_aa = I_NAI_K3x4z_Dxy_aa+ABZ*I_NAI_I3x3z_Dxy_aa;
  Double I_NAI_I2x4y_Fxyz_aa = I_NAI_K2x4yz_Dxy_aa+ABZ*I_NAI_I2x4y_Dxy_aa;
  Double I_NAI_I2x3yz_Fxyz_aa = I_NAI_K2x3y2z_Dxy_aa+ABZ*I_NAI_I2x3yz_Dxy_aa;
  Double I_NAI_I2x2y2z_Fxyz_aa = I_NAI_K2x2y3z_Dxy_aa+ABZ*I_NAI_I2x2y2z_Dxy_aa;
  Double I_NAI_I2xy3z_Fxyz_aa = I_NAI_K2xy4z_Dxy_aa+ABZ*I_NAI_I2xy3z_Dxy_aa;
  Double I_NAI_I2x4z_Fxyz_aa = I_NAI_K2x5z_Dxy_aa+ABZ*I_NAI_I2x4z_Dxy_aa;
  Double I_NAI_Ix5y_Fxyz_aa = I_NAI_Kx5yz_Dxy_aa+ABZ*I_NAI_Ix5y_Dxy_aa;
  Double I_NAI_Ix4yz_Fxyz_aa = I_NAI_Kx4y2z_Dxy_aa+ABZ*I_NAI_Ix4yz_Dxy_aa;
  Double I_NAI_Ix3y2z_Fxyz_aa = I_NAI_Kx3y3z_Dxy_aa+ABZ*I_NAI_Ix3y2z_Dxy_aa;
  Double I_NAI_Ix2y3z_Fxyz_aa = I_NAI_Kx2y4z_Dxy_aa+ABZ*I_NAI_Ix2y3z_Dxy_aa;
  Double I_NAI_Ixy4z_Fxyz_aa = I_NAI_Kxy5z_Dxy_aa+ABZ*I_NAI_Ixy4z_Dxy_aa;
  Double I_NAI_Ix5z_Fxyz_aa = I_NAI_Kx6z_Dxy_aa+ABZ*I_NAI_Ix5z_Dxy_aa;
  Double I_NAI_I6y_Fxyz_aa = I_NAI_K6yz_Dxy_aa+ABZ*I_NAI_I6y_Dxy_aa;
  Double I_NAI_I5yz_Fxyz_aa = I_NAI_K5y2z_Dxy_aa+ABZ*I_NAI_I5yz_Dxy_aa;
  Double I_NAI_I4y2z_Fxyz_aa = I_NAI_K4y3z_Dxy_aa+ABZ*I_NAI_I4y2z_Dxy_aa;
  Double I_NAI_I3y3z_Fxyz_aa = I_NAI_K3y4z_Dxy_aa+ABZ*I_NAI_I3y3z_Dxy_aa;
  Double I_NAI_I2y4z_Fxyz_aa = I_NAI_K2y5z_Dxy_aa+ABZ*I_NAI_I2y4z_Dxy_aa;
  Double I_NAI_Iy5z_Fxyz_aa = I_NAI_Ky6z_Dxy_aa+ABZ*I_NAI_Iy5z_Dxy_aa;
  Double I_NAI_I6z_Fxyz_aa = I_NAI_K7z_Dxy_aa+ABZ*I_NAI_I6z_Dxy_aa;
  Double I_NAI_I6x_Fx2z_aa = I_NAI_K7x_D2z_aa+ABX*I_NAI_I6x_D2z_aa;
  Double I_NAI_I5xy_Fx2z_aa = I_NAI_K6xy_D2z_aa+ABX*I_NAI_I5xy_D2z_aa;
  Double I_NAI_I5xz_Fx2z_aa = I_NAI_K6xz_D2z_aa+ABX*I_NAI_I5xz_D2z_aa;
  Double I_NAI_I4x2y_Fx2z_aa = I_NAI_K5x2y_D2z_aa+ABX*I_NAI_I4x2y_D2z_aa;
  Double I_NAI_I4xyz_Fx2z_aa = I_NAI_K5xyz_D2z_aa+ABX*I_NAI_I4xyz_D2z_aa;
  Double I_NAI_I4x2z_Fx2z_aa = I_NAI_K5x2z_D2z_aa+ABX*I_NAI_I4x2z_D2z_aa;
  Double I_NAI_I3x3y_Fx2z_aa = I_NAI_K4x3y_D2z_aa+ABX*I_NAI_I3x3y_D2z_aa;
  Double I_NAI_I3x2yz_Fx2z_aa = I_NAI_K4x2yz_D2z_aa+ABX*I_NAI_I3x2yz_D2z_aa;
  Double I_NAI_I3xy2z_Fx2z_aa = I_NAI_K4xy2z_D2z_aa+ABX*I_NAI_I3xy2z_D2z_aa;
  Double I_NAI_I3x3z_Fx2z_aa = I_NAI_K4x3z_D2z_aa+ABX*I_NAI_I3x3z_D2z_aa;
  Double I_NAI_I2x4y_Fx2z_aa = I_NAI_K3x4y_D2z_aa+ABX*I_NAI_I2x4y_D2z_aa;
  Double I_NAI_I2x3yz_Fx2z_aa = I_NAI_K3x3yz_D2z_aa+ABX*I_NAI_I2x3yz_D2z_aa;
  Double I_NAI_I2x2y2z_Fx2z_aa = I_NAI_K3x2y2z_D2z_aa+ABX*I_NAI_I2x2y2z_D2z_aa;
  Double I_NAI_I2xy3z_Fx2z_aa = I_NAI_K3xy3z_D2z_aa+ABX*I_NAI_I2xy3z_D2z_aa;
  Double I_NAI_I2x4z_Fx2z_aa = I_NAI_K3x4z_D2z_aa+ABX*I_NAI_I2x4z_D2z_aa;
  Double I_NAI_Ix5y_Fx2z_aa = I_NAI_K2x5y_D2z_aa+ABX*I_NAI_Ix5y_D2z_aa;
  Double I_NAI_Ix4yz_Fx2z_aa = I_NAI_K2x4yz_D2z_aa+ABX*I_NAI_Ix4yz_D2z_aa;
  Double I_NAI_Ix3y2z_Fx2z_aa = I_NAI_K2x3y2z_D2z_aa+ABX*I_NAI_Ix3y2z_D2z_aa;
  Double I_NAI_Ix2y3z_Fx2z_aa = I_NAI_K2x2y3z_D2z_aa+ABX*I_NAI_Ix2y3z_D2z_aa;
  Double I_NAI_Ixy4z_Fx2z_aa = I_NAI_K2xy4z_D2z_aa+ABX*I_NAI_Ixy4z_D2z_aa;
  Double I_NAI_Ix5z_Fx2z_aa = I_NAI_K2x5z_D2z_aa+ABX*I_NAI_Ix5z_D2z_aa;
  Double I_NAI_I6y_Fx2z_aa = I_NAI_Kx6y_D2z_aa+ABX*I_NAI_I6y_D2z_aa;
  Double I_NAI_I5yz_Fx2z_aa = I_NAI_Kx5yz_D2z_aa+ABX*I_NAI_I5yz_D2z_aa;
  Double I_NAI_I4y2z_Fx2z_aa = I_NAI_Kx4y2z_D2z_aa+ABX*I_NAI_I4y2z_D2z_aa;
  Double I_NAI_I3y3z_Fx2z_aa = I_NAI_Kx3y3z_D2z_aa+ABX*I_NAI_I3y3z_D2z_aa;
  Double I_NAI_I2y4z_Fx2z_aa = I_NAI_Kx2y4z_D2z_aa+ABX*I_NAI_I2y4z_D2z_aa;
  Double I_NAI_Iy5z_Fx2z_aa = I_NAI_Kxy5z_D2z_aa+ABX*I_NAI_Iy5z_D2z_aa;
  Double I_NAI_I6z_Fx2z_aa = I_NAI_Kx6z_D2z_aa+ABX*I_NAI_I6z_D2z_aa;
  Double I_NAI_I6x_F3y_aa = I_NAI_K6xy_D2y_aa+ABY*I_NAI_I6x_D2y_aa;
  Double I_NAI_I5xy_F3y_aa = I_NAI_K5x2y_D2y_aa+ABY*I_NAI_I5xy_D2y_aa;
  Double I_NAI_I5xz_F3y_aa = I_NAI_K5xyz_D2y_aa+ABY*I_NAI_I5xz_D2y_aa;
  Double I_NAI_I4x2y_F3y_aa = I_NAI_K4x3y_D2y_aa+ABY*I_NAI_I4x2y_D2y_aa;
  Double I_NAI_I4xyz_F3y_aa = I_NAI_K4x2yz_D2y_aa+ABY*I_NAI_I4xyz_D2y_aa;
  Double I_NAI_I4x2z_F3y_aa = I_NAI_K4xy2z_D2y_aa+ABY*I_NAI_I4x2z_D2y_aa;
  Double I_NAI_I3x3y_F3y_aa = I_NAI_K3x4y_D2y_aa+ABY*I_NAI_I3x3y_D2y_aa;
  Double I_NAI_I3x2yz_F3y_aa = I_NAI_K3x3yz_D2y_aa+ABY*I_NAI_I3x2yz_D2y_aa;
  Double I_NAI_I3xy2z_F3y_aa = I_NAI_K3x2y2z_D2y_aa+ABY*I_NAI_I3xy2z_D2y_aa;
  Double I_NAI_I3x3z_F3y_aa = I_NAI_K3xy3z_D2y_aa+ABY*I_NAI_I3x3z_D2y_aa;
  Double I_NAI_I2x4y_F3y_aa = I_NAI_K2x5y_D2y_aa+ABY*I_NAI_I2x4y_D2y_aa;
  Double I_NAI_I2x3yz_F3y_aa = I_NAI_K2x4yz_D2y_aa+ABY*I_NAI_I2x3yz_D2y_aa;
  Double I_NAI_I2x2y2z_F3y_aa = I_NAI_K2x3y2z_D2y_aa+ABY*I_NAI_I2x2y2z_D2y_aa;
  Double I_NAI_I2xy3z_F3y_aa = I_NAI_K2x2y3z_D2y_aa+ABY*I_NAI_I2xy3z_D2y_aa;
  Double I_NAI_I2x4z_F3y_aa = I_NAI_K2xy4z_D2y_aa+ABY*I_NAI_I2x4z_D2y_aa;
  Double I_NAI_Ix5y_F3y_aa = I_NAI_Kx6y_D2y_aa+ABY*I_NAI_Ix5y_D2y_aa;
  Double I_NAI_Ix4yz_F3y_aa = I_NAI_Kx5yz_D2y_aa+ABY*I_NAI_Ix4yz_D2y_aa;
  Double I_NAI_Ix3y2z_F3y_aa = I_NAI_Kx4y2z_D2y_aa+ABY*I_NAI_Ix3y2z_D2y_aa;
  Double I_NAI_Ix2y3z_F3y_aa = I_NAI_Kx3y3z_D2y_aa+ABY*I_NAI_Ix2y3z_D2y_aa;
  Double I_NAI_Ixy4z_F3y_aa = I_NAI_Kx2y4z_D2y_aa+ABY*I_NAI_Ixy4z_D2y_aa;
  Double I_NAI_Ix5z_F3y_aa = I_NAI_Kxy5z_D2y_aa+ABY*I_NAI_Ix5z_D2y_aa;
  Double I_NAI_I6y_F3y_aa = I_NAI_K7y_D2y_aa+ABY*I_NAI_I6y_D2y_aa;
  Double I_NAI_I5yz_F3y_aa = I_NAI_K6yz_D2y_aa+ABY*I_NAI_I5yz_D2y_aa;
  Double I_NAI_I4y2z_F3y_aa = I_NAI_K5y2z_D2y_aa+ABY*I_NAI_I4y2z_D2y_aa;
  Double I_NAI_I3y3z_F3y_aa = I_NAI_K4y3z_D2y_aa+ABY*I_NAI_I3y3z_D2y_aa;
  Double I_NAI_I2y4z_F3y_aa = I_NAI_K3y4z_D2y_aa+ABY*I_NAI_I2y4z_D2y_aa;
  Double I_NAI_Iy5z_F3y_aa = I_NAI_K2y5z_D2y_aa+ABY*I_NAI_Iy5z_D2y_aa;
  Double I_NAI_I6z_F3y_aa = I_NAI_Ky6z_D2y_aa+ABY*I_NAI_I6z_D2y_aa;
  Double I_NAI_I6x_F2yz_aa = I_NAI_K6xz_D2y_aa+ABZ*I_NAI_I6x_D2y_aa;
  Double I_NAI_I5xy_F2yz_aa = I_NAI_K5xyz_D2y_aa+ABZ*I_NAI_I5xy_D2y_aa;
  Double I_NAI_I5xz_F2yz_aa = I_NAI_K5x2z_D2y_aa+ABZ*I_NAI_I5xz_D2y_aa;
  Double I_NAI_I4x2y_F2yz_aa = I_NAI_K4x2yz_D2y_aa+ABZ*I_NAI_I4x2y_D2y_aa;
  Double I_NAI_I4xyz_F2yz_aa = I_NAI_K4xy2z_D2y_aa+ABZ*I_NAI_I4xyz_D2y_aa;
  Double I_NAI_I4x2z_F2yz_aa = I_NAI_K4x3z_D2y_aa+ABZ*I_NAI_I4x2z_D2y_aa;
  Double I_NAI_I3x3y_F2yz_aa = I_NAI_K3x3yz_D2y_aa+ABZ*I_NAI_I3x3y_D2y_aa;
  Double I_NAI_I3x2yz_F2yz_aa = I_NAI_K3x2y2z_D2y_aa+ABZ*I_NAI_I3x2yz_D2y_aa;
  Double I_NAI_I3xy2z_F2yz_aa = I_NAI_K3xy3z_D2y_aa+ABZ*I_NAI_I3xy2z_D2y_aa;
  Double I_NAI_I3x3z_F2yz_aa = I_NAI_K3x4z_D2y_aa+ABZ*I_NAI_I3x3z_D2y_aa;
  Double I_NAI_I2x4y_F2yz_aa = I_NAI_K2x4yz_D2y_aa+ABZ*I_NAI_I2x4y_D2y_aa;
  Double I_NAI_I2x3yz_F2yz_aa = I_NAI_K2x3y2z_D2y_aa+ABZ*I_NAI_I2x3yz_D2y_aa;
  Double I_NAI_I2x2y2z_F2yz_aa = I_NAI_K2x2y3z_D2y_aa+ABZ*I_NAI_I2x2y2z_D2y_aa;
  Double I_NAI_I2xy3z_F2yz_aa = I_NAI_K2xy4z_D2y_aa+ABZ*I_NAI_I2xy3z_D2y_aa;
  Double I_NAI_I2x4z_F2yz_aa = I_NAI_K2x5z_D2y_aa+ABZ*I_NAI_I2x4z_D2y_aa;
  Double I_NAI_Ix5y_F2yz_aa = I_NAI_Kx5yz_D2y_aa+ABZ*I_NAI_Ix5y_D2y_aa;
  Double I_NAI_Ix4yz_F2yz_aa = I_NAI_Kx4y2z_D2y_aa+ABZ*I_NAI_Ix4yz_D2y_aa;
  Double I_NAI_Ix3y2z_F2yz_aa = I_NAI_Kx3y3z_D2y_aa+ABZ*I_NAI_Ix3y2z_D2y_aa;
  Double I_NAI_Ix2y3z_F2yz_aa = I_NAI_Kx2y4z_D2y_aa+ABZ*I_NAI_Ix2y3z_D2y_aa;
  Double I_NAI_Ixy4z_F2yz_aa = I_NAI_Kxy5z_D2y_aa+ABZ*I_NAI_Ixy4z_D2y_aa;
  Double I_NAI_Ix5z_F2yz_aa = I_NAI_Kx6z_D2y_aa+ABZ*I_NAI_Ix5z_D2y_aa;
  Double I_NAI_I6y_F2yz_aa = I_NAI_K6yz_D2y_aa+ABZ*I_NAI_I6y_D2y_aa;
  Double I_NAI_I5yz_F2yz_aa = I_NAI_K5y2z_D2y_aa+ABZ*I_NAI_I5yz_D2y_aa;
  Double I_NAI_I4y2z_F2yz_aa = I_NAI_K4y3z_D2y_aa+ABZ*I_NAI_I4y2z_D2y_aa;
  Double I_NAI_I3y3z_F2yz_aa = I_NAI_K3y4z_D2y_aa+ABZ*I_NAI_I3y3z_D2y_aa;
  Double I_NAI_I2y4z_F2yz_aa = I_NAI_K2y5z_D2y_aa+ABZ*I_NAI_I2y4z_D2y_aa;
  Double I_NAI_Iy5z_F2yz_aa = I_NAI_Ky6z_D2y_aa+ABZ*I_NAI_Iy5z_D2y_aa;
  Double I_NAI_I6z_F2yz_aa = I_NAI_K7z_D2y_aa+ABZ*I_NAI_I6z_D2y_aa;
  Double I_NAI_I6x_Fy2z_aa = I_NAI_K6xy_D2z_aa+ABY*I_NAI_I6x_D2z_aa;
  Double I_NAI_I5xy_Fy2z_aa = I_NAI_K5x2y_D2z_aa+ABY*I_NAI_I5xy_D2z_aa;
  Double I_NAI_I5xz_Fy2z_aa = I_NAI_K5xyz_D2z_aa+ABY*I_NAI_I5xz_D2z_aa;
  Double I_NAI_I4x2y_Fy2z_aa = I_NAI_K4x3y_D2z_aa+ABY*I_NAI_I4x2y_D2z_aa;
  Double I_NAI_I4xyz_Fy2z_aa = I_NAI_K4x2yz_D2z_aa+ABY*I_NAI_I4xyz_D2z_aa;
  Double I_NAI_I4x2z_Fy2z_aa = I_NAI_K4xy2z_D2z_aa+ABY*I_NAI_I4x2z_D2z_aa;
  Double I_NAI_I3x3y_Fy2z_aa = I_NAI_K3x4y_D2z_aa+ABY*I_NAI_I3x3y_D2z_aa;
  Double I_NAI_I3x2yz_Fy2z_aa = I_NAI_K3x3yz_D2z_aa+ABY*I_NAI_I3x2yz_D2z_aa;
  Double I_NAI_I3xy2z_Fy2z_aa = I_NAI_K3x2y2z_D2z_aa+ABY*I_NAI_I3xy2z_D2z_aa;
  Double I_NAI_I3x3z_Fy2z_aa = I_NAI_K3xy3z_D2z_aa+ABY*I_NAI_I3x3z_D2z_aa;
  Double I_NAI_I2x4y_Fy2z_aa = I_NAI_K2x5y_D2z_aa+ABY*I_NAI_I2x4y_D2z_aa;
  Double I_NAI_I2x3yz_Fy2z_aa = I_NAI_K2x4yz_D2z_aa+ABY*I_NAI_I2x3yz_D2z_aa;
  Double I_NAI_I2x2y2z_Fy2z_aa = I_NAI_K2x3y2z_D2z_aa+ABY*I_NAI_I2x2y2z_D2z_aa;
  Double I_NAI_I2xy3z_Fy2z_aa = I_NAI_K2x2y3z_D2z_aa+ABY*I_NAI_I2xy3z_D2z_aa;
  Double I_NAI_I2x4z_Fy2z_aa = I_NAI_K2xy4z_D2z_aa+ABY*I_NAI_I2x4z_D2z_aa;
  Double I_NAI_Ix5y_Fy2z_aa = I_NAI_Kx6y_D2z_aa+ABY*I_NAI_Ix5y_D2z_aa;
  Double I_NAI_Ix4yz_Fy2z_aa = I_NAI_Kx5yz_D2z_aa+ABY*I_NAI_Ix4yz_D2z_aa;
  Double I_NAI_Ix3y2z_Fy2z_aa = I_NAI_Kx4y2z_D2z_aa+ABY*I_NAI_Ix3y2z_D2z_aa;
  Double I_NAI_Ix2y3z_Fy2z_aa = I_NAI_Kx3y3z_D2z_aa+ABY*I_NAI_Ix2y3z_D2z_aa;
  Double I_NAI_Ixy4z_Fy2z_aa = I_NAI_Kx2y4z_D2z_aa+ABY*I_NAI_Ixy4z_D2z_aa;
  Double I_NAI_Ix5z_Fy2z_aa = I_NAI_Kxy5z_D2z_aa+ABY*I_NAI_Ix5z_D2z_aa;
  Double I_NAI_I6y_Fy2z_aa = I_NAI_K7y_D2z_aa+ABY*I_NAI_I6y_D2z_aa;
  Double I_NAI_I5yz_Fy2z_aa = I_NAI_K6yz_D2z_aa+ABY*I_NAI_I5yz_D2z_aa;
  Double I_NAI_I4y2z_Fy2z_aa = I_NAI_K5y2z_D2z_aa+ABY*I_NAI_I4y2z_D2z_aa;
  Double I_NAI_I3y3z_Fy2z_aa = I_NAI_K4y3z_D2z_aa+ABY*I_NAI_I3y3z_D2z_aa;
  Double I_NAI_I2y4z_Fy2z_aa = I_NAI_K3y4z_D2z_aa+ABY*I_NAI_I2y4z_D2z_aa;
  Double I_NAI_Iy5z_Fy2z_aa = I_NAI_K2y5z_D2z_aa+ABY*I_NAI_Iy5z_D2z_aa;
  Double I_NAI_I6z_Fy2z_aa = I_NAI_Ky6z_D2z_aa+ABY*I_NAI_I6z_D2z_aa;
  Double I_NAI_I6x_F3z_aa = I_NAI_K6xz_D2z_aa+ABZ*I_NAI_I6x_D2z_aa;
  Double I_NAI_I5xy_F3z_aa = I_NAI_K5xyz_D2z_aa+ABZ*I_NAI_I5xy_D2z_aa;
  Double I_NAI_I5xz_F3z_aa = I_NAI_K5x2z_D2z_aa+ABZ*I_NAI_I5xz_D2z_aa;
  Double I_NAI_I4x2y_F3z_aa = I_NAI_K4x2yz_D2z_aa+ABZ*I_NAI_I4x2y_D2z_aa;
  Double I_NAI_I4xyz_F3z_aa = I_NAI_K4xy2z_D2z_aa+ABZ*I_NAI_I4xyz_D2z_aa;
  Double I_NAI_I4x2z_F3z_aa = I_NAI_K4x3z_D2z_aa+ABZ*I_NAI_I4x2z_D2z_aa;
  Double I_NAI_I3x3y_F3z_aa = I_NAI_K3x3yz_D2z_aa+ABZ*I_NAI_I3x3y_D2z_aa;
  Double I_NAI_I3x2yz_F3z_aa = I_NAI_K3x2y2z_D2z_aa+ABZ*I_NAI_I3x2yz_D2z_aa;
  Double I_NAI_I3xy2z_F3z_aa = I_NAI_K3xy3z_D2z_aa+ABZ*I_NAI_I3xy2z_D2z_aa;
  Double I_NAI_I3x3z_F3z_aa = I_NAI_K3x4z_D2z_aa+ABZ*I_NAI_I3x3z_D2z_aa;
  Double I_NAI_I2x4y_F3z_aa = I_NAI_K2x4yz_D2z_aa+ABZ*I_NAI_I2x4y_D2z_aa;
  Double I_NAI_I2x3yz_F3z_aa = I_NAI_K2x3y2z_D2z_aa+ABZ*I_NAI_I2x3yz_D2z_aa;
  Double I_NAI_I2x2y2z_F3z_aa = I_NAI_K2x2y3z_D2z_aa+ABZ*I_NAI_I2x2y2z_D2z_aa;
  Double I_NAI_I2xy3z_F3z_aa = I_NAI_K2xy4z_D2z_aa+ABZ*I_NAI_I2xy3z_D2z_aa;
  Double I_NAI_I2x4z_F3z_aa = I_NAI_K2x5z_D2z_aa+ABZ*I_NAI_I2x4z_D2z_aa;
  Double I_NAI_Ix5y_F3z_aa = I_NAI_Kx5yz_D2z_aa+ABZ*I_NAI_Ix5y_D2z_aa;
  Double I_NAI_Ix4yz_F3z_aa = I_NAI_Kx4y2z_D2z_aa+ABZ*I_NAI_Ix4yz_D2z_aa;
  Double I_NAI_Ix3y2z_F3z_aa = I_NAI_Kx3y3z_D2z_aa+ABZ*I_NAI_Ix3y2z_D2z_aa;
  Double I_NAI_Ix2y3z_F3z_aa = I_NAI_Kx2y4z_D2z_aa+ABZ*I_NAI_Ix2y3z_D2z_aa;
  Double I_NAI_Ixy4z_F3z_aa = I_NAI_Kxy5z_D2z_aa+ABZ*I_NAI_Ixy4z_D2z_aa;
  Double I_NAI_Ix5z_F3z_aa = I_NAI_Kx6z_D2z_aa+ABZ*I_NAI_Ix5z_D2z_aa;
  Double I_NAI_I6y_F3z_aa = I_NAI_K6yz_D2z_aa+ABZ*I_NAI_I6y_D2z_aa;
  Double I_NAI_I5yz_F3z_aa = I_NAI_K5y2z_D2z_aa+ABZ*I_NAI_I5yz_D2z_aa;
  Double I_NAI_I4y2z_F3z_aa = I_NAI_K4y3z_D2z_aa+ABZ*I_NAI_I4y2z_D2z_aa;
  Double I_NAI_I3y3z_F3z_aa = I_NAI_K3y4z_D2z_aa+ABZ*I_NAI_I3y3z_D2z_aa;
  Double I_NAI_I2y4z_F3z_aa = I_NAI_K2y5z_D2z_aa+ABZ*I_NAI_I2y4z_D2z_aa;
  Double I_NAI_Iy5z_F3z_aa = I_NAI_Ky6z_D2z_aa+ABZ*I_NAI_Iy5z_D2z_aa;
  Double I_NAI_I6z_F3z_aa = I_NAI_K7z_D2z_aa+ABZ*I_NAI_I6z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_ab
   * RHS shell quartet name: SQ_NAI_H_S_ab
   ************************************************************/
  Double I_NAI_H5x_Px_ab = I_NAI_I6x_S_ab+ABX*I_NAI_H5x_S_ab;
  Double I_NAI_H4xy_Px_ab = I_NAI_I5xy_S_ab+ABX*I_NAI_H4xy_S_ab;
  Double I_NAI_H4xz_Px_ab = I_NAI_I5xz_S_ab+ABX*I_NAI_H4xz_S_ab;
  Double I_NAI_H3x2y_Px_ab = I_NAI_I4x2y_S_ab+ABX*I_NAI_H3x2y_S_ab;
  Double I_NAI_H3xyz_Px_ab = I_NAI_I4xyz_S_ab+ABX*I_NAI_H3xyz_S_ab;
  Double I_NAI_H3x2z_Px_ab = I_NAI_I4x2z_S_ab+ABX*I_NAI_H3x2z_S_ab;
  Double I_NAI_H2x3y_Px_ab = I_NAI_I3x3y_S_ab+ABX*I_NAI_H2x3y_S_ab;
  Double I_NAI_H2x2yz_Px_ab = I_NAI_I3x2yz_S_ab+ABX*I_NAI_H2x2yz_S_ab;
  Double I_NAI_H2xy2z_Px_ab = I_NAI_I3xy2z_S_ab+ABX*I_NAI_H2xy2z_S_ab;
  Double I_NAI_H2x3z_Px_ab = I_NAI_I3x3z_S_ab+ABX*I_NAI_H2x3z_S_ab;
  Double I_NAI_Hx4y_Px_ab = I_NAI_I2x4y_S_ab+ABX*I_NAI_Hx4y_S_ab;
  Double I_NAI_Hx3yz_Px_ab = I_NAI_I2x3yz_S_ab+ABX*I_NAI_Hx3yz_S_ab;
  Double I_NAI_Hx2y2z_Px_ab = I_NAI_I2x2y2z_S_ab+ABX*I_NAI_Hx2y2z_S_ab;
  Double I_NAI_Hxy3z_Px_ab = I_NAI_I2xy3z_S_ab+ABX*I_NAI_Hxy3z_S_ab;
  Double I_NAI_Hx4z_Px_ab = I_NAI_I2x4z_S_ab+ABX*I_NAI_Hx4z_S_ab;
  Double I_NAI_H5y_Px_ab = I_NAI_Ix5y_S_ab+ABX*I_NAI_H5y_S_ab;
  Double I_NAI_H4yz_Px_ab = I_NAI_Ix4yz_S_ab+ABX*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Px_ab = I_NAI_Ix3y2z_S_ab+ABX*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Px_ab = I_NAI_Ix2y3z_S_ab+ABX*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Px_ab = I_NAI_Ixy4z_S_ab+ABX*I_NAI_Hy4z_S_ab;
  Double I_NAI_H5z_Px_ab = I_NAI_Ix5z_S_ab+ABX*I_NAI_H5z_S_ab;
  Double I_NAI_H5x_Py_ab = I_NAI_I5xy_S_ab+ABY*I_NAI_H5x_S_ab;
  Double I_NAI_H4xy_Py_ab = I_NAI_I4x2y_S_ab+ABY*I_NAI_H4xy_S_ab;
  Double I_NAI_H4xz_Py_ab = I_NAI_I4xyz_S_ab+ABY*I_NAI_H4xz_S_ab;
  Double I_NAI_H3x2y_Py_ab = I_NAI_I3x3y_S_ab+ABY*I_NAI_H3x2y_S_ab;
  Double I_NAI_H3xyz_Py_ab = I_NAI_I3x2yz_S_ab+ABY*I_NAI_H3xyz_S_ab;
  Double I_NAI_H3x2z_Py_ab = I_NAI_I3xy2z_S_ab+ABY*I_NAI_H3x2z_S_ab;
  Double I_NAI_H2x3y_Py_ab = I_NAI_I2x4y_S_ab+ABY*I_NAI_H2x3y_S_ab;
  Double I_NAI_H2x2yz_Py_ab = I_NAI_I2x3yz_S_ab+ABY*I_NAI_H2x2yz_S_ab;
  Double I_NAI_H2xy2z_Py_ab = I_NAI_I2x2y2z_S_ab+ABY*I_NAI_H2xy2z_S_ab;
  Double I_NAI_H2x3z_Py_ab = I_NAI_I2xy3z_S_ab+ABY*I_NAI_H2x3z_S_ab;
  Double I_NAI_Hx4y_Py_ab = I_NAI_Ix5y_S_ab+ABY*I_NAI_Hx4y_S_ab;
  Double I_NAI_Hx3yz_Py_ab = I_NAI_Ix4yz_S_ab+ABY*I_NAI_Hx3yz_S_ab;
  Double I_NAI_Hx2y2z_Py_ab = I_NAI_Ix3y2z_S_ab+ABY*I_NAI_Hx2y2z_S_ab;
  Double I_NAI_Hxy3z_Py_ab = I_NAI_Ix2y3z_S_ab+ABY*I_NAI_Hxy3z_S_ab;
  Double I_NAI_Hx4z_Py_ab = I_NAI_Ixy4z_S_ab+ABY*I_NAI_Hx4z_S_ab;
  Double I_NAI_H5y_Py_ab = I_NAI_I6y_S_ab+ABY*I_NAI_H5y_S_ab;
  Double I_NAI_H4yz_Py_ab = I_NAI_I5yz_S_ab+ABY*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Py_ab = I_NAI_I4y2z_S_ab+ABY*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Py_ab = I_NAI_I3y3z_S_ab+ABY*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Py_ab = I_NAI_I2y4z_S_ab+ABY*I_NAI_Hy4z_S_ab;
  Double I_NAI_H5z_Py_ab = I_NAI_Iy5z_S_ab+ABY*I_NAI_H5z_S_ab;
  Double I_NAI_H5x_Pz_ab = I_NAI_I5xz_S_ab+ABZ*I_NAI_H5x_S_ab;
  Double I_NAI_H4xy_Pz_ab = I_NAI_I4xyz_S_ab+ABZ*I_NAI_H4xy_S_ab;
  Double I_NAI_H4xz_Pz_ab = I_NAI_I4x2z_S_ab+ABZ*I_NAI_H4xz_S_ab;
  Double I_NAI_H3x2y_Pz_ab = I_NAI_I3x2yz_S_ab+ABZ*I_NAI_H3x2y_S_ab;
  Double I_NAI_H3xyz_Pz_ab = I_NAI_I3xy2z_S_ab+ABZ*I_NAI_H3xyz_S_ab;
  Double I_NAI_H3x2z_Pz_ab = I_NAI_I3x3z_S_ab+ABZ*I_NAI_H3x2z_S_ab;
  Double I_NAI_H2x3y_Pz_ab = I_NAI_I2x3yz_S_ab+ABZ*I_NAI_H2x3y_S_ab;
  Double I_NAI_H2x2yz_Pz_ab = I_NAI_I2x2y2z_S_ab+ABZ*I_NAI_H2x2yz_S_ab;
  Double I_NAI_H2xy2z_Pz_ab = I_NAI_I2xy3z_S_ab+ABZ*I_NAI_H2xy2z_S_ab;
  Double I_NAI_H2x3z_Pz_ab = I_NAI_I2x4z_S_ab+ABZ*I_NAI_H2x3z_S_ab;
  Double I_NAI_Hx4y_Pz_ab = I_NAI_Ix4yz_S_ab+ABZ*I_NAI_Hx4y_S_ab;
  Double I_NAI_Hx3yz_Pz_ab = I_NAI_Ix3y2z_S_ab+ABZ*I_NAI_Hx3yz_S_ab;
  Double I_NAI_Hx2y2z_Pz_ab = I_NAI_Ix2y3z_S_ab+ABZ*I_NAI_Hx2y2z_S_ab;
  Double I_NAI_Hxy3z_Pz_ab = I_NAI_Ixy4z_S_ab+ABZ*I_NAI_Hxy3z_S_ab;
  Double I_NAI_Hx4z_Pz_ab = I_NAI_Ix5z_S_ab+ABZ*I_NAI_Hx4z_S_ab;
  Double I_NAI_H5y_Pz_ab = I_NAI_I5yz_S_ab+ABZ*I_NAI_H5y_S_ab;
  Double I_NAI_H4yz_Pz_ab = I_NAI_I4y2z_S_ab+ABZ*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Pz_ab = I_NAI_I3y3z_S_ab+ABZ*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Pz_ab = I_NAI_I2y4z_S_ab+ABZ*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Pz_ab = I_NAI_Iy5z_S_ab+ABZ*I_NAI_Hy4z_S_ab;
  Double I_NAI_H5z_Pz_ab = I_NAI_I6z_S_ab+ABZ*I_NAI_H5z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_ab
   * RHS shell quartet name: SQ_NAI_I_S_ab
   ************************************************************/
  Double I_NAI_I6x_Px_ab = I_NAI_K7x_S_ab+ABX*I_NAI_I6x_S_ab;
  Double I_NAI_I5xy_Px_ab = I_NAI_K6xy_S_ab+ABX*I_NAI_I5xy_S_ab;
  Double I_NAI_I5xz_Px_ab = I_NAI_K6xz_S_ab+ABX*I_NAI_I5xz_S_ab;
  Double I_NAI_I4x2y_Px_ab = I_NAI_K5x2y_S_ab+ABX*I_NAI_I4x2y_S_ab;
  Double I_NAI_I4xyz_Px_ab = I_NAI_K5xyz_S_ab+ABX*I_NAI_I4xyz_S_ab;
  Double I_NAI_I4x2z_Px_ab = I_NAI_K5x2z_S_ab+ABX*I_NAI_I4x2z_S_ab;
  Double I_NAI_I3x3y_Px_ab = I_NAI_K4x3y_S_ab+ABX*I_NAI_I3x3y_S_ab;
  Double I_NAI_I3x2yz_Px_ab = I_NAI_K4x2yz_S_ab+ABX*I_NAI_I3x2yz_S_ab;
  Double I_NAI_I3xy2z_Px_ab = I_NAI_K4xy2z_S_ab+ABX*I_NAI_I3xy2z_S_ab;
  Double I_NAI_I3x3z_Px_ab = I_NAI_K4x3z_S_ab+ABX*I_NAI_I3x3z_S_ab;
  Double I_NAI_I2x4y_Px_ab = I_NAI_K3x4y_S_ab+ABX*I_NAI_I2x4y_S_ab;
  Double I_NAI_I2x3yz_Px_ab = I_NAI_K3x3yz_S_ab+ABX*I_NAI_I2x3yz_S_ab;
  Double I_NAI_I2x2y2z_Px_ab = I_NAI_K3x2y2z_S_ab+ABX*I_NAI_I2x2y2z_S_ab;
  Double I_NAI_I2xy3z_Px_ab = I_NAI_K3xy3z_S_ab+ABX*I_NAI_I2xy3z_S_ab;
  Double I_NAI_I2x4z_Px_ab = I_NAI_K3x4z_S_ab+ABX*I_NAI_I2x4z_S_ab;
  Double I_NAI_Ix5y_Px_ab = I_NAI_K2x5y_S_ab+ABX*I_NAI_Ix5y_S_ab;
  Double I_NAI_Ix4yz_Px_ab = I_NAI_K2x4yz_S_ab+ABX*I_NAI_Ix4yz_S_ab;
  Double I_NAI_Ix3y2z_Px_ab = I_NAI_K2x3y2z_S_ab+ABX*I_NAI_Ix3y2z_S_ab;
  Double I_NAI_Ix2y3z_Px_ab = I_NAI_K2x2y3z_S_ab+ABX*I_NAI_Ix2y3z_S_ab;
  Double I_NAI_Ixy4z_Px_ab = I_NAI_K2xy4z_S_ab+ABX*I_NAI_Ixy4z_S_ab;
  Double I_NAI_Ix5z_Px_ab = I_NAI_K2x5z_S_ab+ABX*I_NAI_Ix5z_S_ab;
  Double I_NAI_I6y_Px_ab = I_NAI_Kx6y_S_ab+ABX*I_NAI_I6y_S_ab;
  Double I_NAI_I5yz_Px_ab = I_NAI_Kx5yz_S_ab+ABX*I_NAI_I5yz_S_ab;
  Double I_NAI_I4y2z_Px_ab = I_NAI_Kx4y2z_S_ab+ABX*I_NAI_I4y2z_S_ab;
  Double I_NAI_I3y3z_Px_ab = I_NAI_Kx3y3z_S_ab+ABX*I_NAI_I3y3z_S_ab;
  Double I_NAI_I2y4z_Px_ab = I_NAI_Kx2y4z_S_ab+ABX*I_NAI_I2y4z_S_ab;
  Double I_NAI_Iy5z_Px_ab = I_NAI_Kxy5z_S_ab+ABX*I_NAI_Iy5z_S_ab;
  Double I_NAI_I6z_Px_ab = I_NAI_Kx6z_S_ab+ABX*I_NAI_I6z_S_ab;
  Double I_NAI_I6x_Py_ab = I_NAI_K6xy_S_ab+ABY*I_NAI_I6x_S_ab;
  Double I_NAI_I5xy_Py_ab = I_NAI_K5x2y_S_ab+ABY*I_NAI_I5xy_S_ab;
  Double I_NAI_I5xz_Py_ab = I_NAI_K5xyz_S_ab+ABY*I_NAI_I5xz_S_ab;
  Double I_NAI_I4x2y_Py_ab = I_NAI_K4x3y_S_ab+ABY*I_NAI_I4x2y_S_ab;
  Double I_NAI_I4xyz_Py_ab = I_NAI_K4x2yz_S_ab+ABY*I_NAI_I4xyz_S_ab;
  Double I_NAI_I4x2z_Py_ab = I_NAI_K4xy2z_S_ab+ABY*I_NAI_I4x2z_S_ab;
  Double I_NAI_I3x3y_Py_ab = I_NAI_K3x4y_S_ab+ABY*I_NAI_I3x3y_S_ab;
  Double I_NAI_I3x2yz_Py_ab = I_NAI_K3x3yz_S_ab+ABY*I_NAI_I3x2yz_S_ab;
  Double I_NAI_I3xy2z_Py_ab = I_NAI_K3x2y2z_S_ab+ABY*I_NAI_I3xy2z_S_ab;
  Double I_NAI_I3x3z_Py_ab = I_NAI_K3xy3z_S_ab+ABY*I_NAI_I3x3z_S_ab;
  Double I_NAI_I2x4y_Py_ab = I_NAI_K2x5y_S_ab+ABY*I_NAI_I2x4y_S_ab;
  Double I_NAI_I2x3yz_Py_ab = I_NAI_K2x4yz_S_ab+ABY*I_NAI_I2x3yz_S_ab;
  Double I_NAI_I2x2y2z_Py_ab = I_NAI_K2x3y2z_S_ab+ABY*I_NAI_I2x2y2z_S_ab;
  Double I_NAI_I2xy3z_Py_ab = I_NAI_K2x2y3z_S_ab+ABY*I_NAI_I2xy3z_S_ab;
  Double I_NAI_I2x4z_Py_ab = I_NAI_K2xy4z_S_ab+ABY*I_NAI_I2x4z_S_ab;
  Double I_NAI_Ix5y_Py_ab = I_NAI_Kx6y_S_ab+ABY*I_NAI_Ix5y_S_ab;
  Double I_NAI_Ix4yz_Py_ab = I_NAI_Kx5yz_S_ab+ABY*I_NAI_Ix4yz_S_ab;
  Double I_NAI_Ix3y2z_Py_ab = I_NAI_Kx4y2z_S_ab+ABY*I_NAI_Ix3y2z_S_ab;
  Double I_NAI_Ix2y3z_Py_ab = I_NAI_Kx3y3z_S_ab+ABY*I_NAI_Ix2y3z_S_ab;
  Double I_NAI_Ixy4z_Py_ab = I_NAI_Kx2y4z_S_ab+ABY*I_NAI_Ixy4z_S_ab;
  Double I_NAI_Ix5z_Py_ab = I_NAI_Kxy5z_S_ab+ABY*I_NAI_Ix5z_S_ab;
  Double I_NAI_I6y_Py_ab = I_NAI_K7y_S_ab+ABY*I_NAI_I6y_S_ab;
  Double I_NAI_I5yz_Py_ab = I_NAI_K6yz_S_ab+ABY*I_NAI_I5yz_S_ab;
  Double I_NAI_I4y2z_Py_ab = I_NAI_K5y2z_S_ab+ABY*I_NAI_I4y2z_S_ab;
  Double I_NAI_I3y3z_Py_ab = I_NAI_K4y3z_S_ab+ABY*I_NAI_I3y3z_S_ab;
  Double I_NAI_I2y4z_Py_ab = I_NAI_K3y4z_S_ab+ABY*I_NAI_I2y4z_S_ab;
  Double I_NAI_Iy5z_Py_ab = I_NAI_K2y5z_S_ab+ABY*I_NAI_Iy5z_S_ab;
  Double I_NAI_I6z_Py_ab = I_NAI_Ky6z_S_ab+ABY*I_NAI_I6z_S_ab;
  Double I_NAI_I6x_Pz_ab = I_NAI_K6xz_S_ab+ABZ*I_NAI_I6x_S_ab;
  Double I_NAI_I5xy_Pz_ab = I_NAI_K5xyz_S_ab+ABZ*I_NAI_I5xy_S_ab;
  Double I_NAI_I5xz_Pz_ab = I_NAI_K5x2z_S_ab+ABZ*I_NAI_I5xz_S_ab;
  Double I_NAI_I4x2y_Pz_ab = I_NAI_K4x2yz_S_ab+ABZ*I_NAI_I4x2y_S_ab;
  Double I_NAI_I4xyz_Pz_ab = I_NAI_K4xy2z_S_ab+ABZ*I_NAI_I4xyz_S_ab;
  Double I_NAI_I4x2z_Pz_ab = I_NAI_K4x3z_S_ab+ABZ*I_NAI_I4x2z_S_ab;
  Double I_NAI_I3x3y_Pz_ab = I_NAI_K3x3yz_S_ab+ABZ*I_NAI_I3x3y_S_ab;
  Double I_NAI_I3x2yz_Pz_ab = I_NAI_K3x2y2z_S_ab+ABZ*I_NAI_I3x2yz_S_ab;
  Double I_NAI_I3xy2z_Pz_ab = I_NAI_K3xy3z_S_ab+ABZ*I_NAI_I3xy2z_S_ab;
  Double I_NAI_I3x3z_Pz_ab = I_NAI_K3x4z_S_ab+ABZ*I_NAI_I3x3z_S_ab;
  Double I_NAI_I2x4y_Pz_ab = I_NAI_K2x4yz_S_ab+ABZ*I_NAI_I2x4y_S_ab;
  Double I_NAI_I2x3yz_Pz_ab = I_NAI_K2x3y2z_S_ab+ABZ*I_NAI_I2x3yz_S_ab;
  Double I_NAI_I2x2y2z_Pz_ab = I_NAI_K2x2y3z_S_ab+ABZ*I_NAI_I2x2y2z_S_ab;
  Double I_NAI_I2xy3z_Pz_ab = I_NAI_K2xy4z_S_ab+ABZ*I_NAI_I2xy3z_S_ab;
  Double I_NAI_I2x4z_Pz_ab = I_NAI_K2x5z_S_ab+ABZ*I_NAI_I2x4z_S_ab;
  Double I_NAI_Ix5y_Pz_ab = I_NAI_Kx5yz_S_ab+ABZ*I_NAI_Ix5y_S_ab;
  Double I_NAI_Ix4yz_Pz_ab = I_NAI_Kx4y2z_S_ab+ABZ*I_NAI_Ix4yz_S_ab;
  Double I_NAI_Ix3y2z_Pz_ab = I_NAI_Kx3y3z_S_ab+ABZ*I_NAI_Ix3y2z_S_ab;
  Double I_NAI_Ix2y3z_Pz_ab = I_NAI_Kx2y4z_S_ab+ABZ*I_NAI_Ix2y3z_S_ab;
  Double I_NAI_Ixy4z_Pz_ab = I_NAI_Kxy5z_S_ab+ABZ*I_NAI_Ixy4z_S_ab;
  Double I_NAI_Ix5z_Pz_ab = I_NAI_Kx6z_S_ab+ABZ*I_NAI_Ix5z_S_ab;
  Double I_NAI_I6y_Pz_ab = I_NAI_K6yz_S_ab+ABZ*I_NAI_I6y_S_ab;
  Double I_NAI_I5yz_Pz_ab = I_NAI_K5y2z_S_ab+ABZ*I_NAI_I5yz_S_ab;
  Double I_NAI_I4y2z_Pz_ab = I_NAI_K4y3z_S_ab+ABZ*I_NAI_I4y2z_S_ab;
  Double I_NAI_I3y3z_Pz_ab = I_NAI_K3y4z_S_ab+ABZ*I_NAI_I3y3z_S_ab;
  Double I_NAI_I2y4z_Pz_ab = I_NAI_K2y5z_S_ab+ABZ*I_NAI_I2y4z_S_ab;
  Double I_NAI_Iy5z_Pz_ab = I_NAI_Ky6z_S_ab+ABZ*I_NAI_Iy5z_S_ab;
  Double I_NAI_I6z_Pz_ab = I_NAI_K7z_S_ab+ABZ*I_NAI_I6z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_ab
   * RHS shell quartet name: SQ_NAI_H_P_ab
   ************************************************************/
  Double I_NAI_H5x_D2x_ab = I_NAI_I6x_Px_ab+ABX*I_NAI_H5x_Px_ab;
  Double I_NAI_H4xy_D2x_ab = I_NAI_I5xy_Px_ab+ABX*I_NAI_H4xy_Px_ab;
  Double I_NAI_H4xz_D2x_ab = I_NAI_I5xz_Px_ab+ABX*I_NAI_H4xz_Px_ab;
  Double I_NAI_H3x2y_D2x_ab = I_NAI_I4x2y_Px_ab+ABX*I_NAI_H3x2y_Px_ab;
  Double I_NAI_H3xyz_D2x_ab = I_NAI_I4xyz_Px_ab+ABX*I_NAI_H3xyz_Px_ab;
  Double I_NAI_H3x2z_D2x_ab = I_NAI_I4x2z_Px_ab+ABX*I_NAI_H3x2z_Px_ab;
  Double I_NAI_H2x3y_D2x_ab = I_NAI_I3x3y_Px_ab+ABX*I_NAI_H2x3y_Px_ab;
  Double I_NAI_H2x2yz_D2x_ab = I_NAI_I3x2yz_Px_ab+ABX*I_NAI_H2x2yz_Px_ab;
  Double I_NAI_H2xy2z_D2x_ab = I_NAI_I3xy2z_Px_ab+ABX*I_NAI_H2xy2z_Px_ab;
  Double I_NAI_H2x3z_D2x_ab = I_NAI_I3x3z_Px_ab+ABX*I_NAI_H2x3z_Px_ab;
  Double I_NAI_Hx4y_D2x_ab = I_NAI_I2x4y_Px_ab+ABX*I_NAI_Hx4y_Px_ab;
  Double I_NAI_Hx3yz_D2x_ab = I_NAI_I2x3yz_Px_ab+ABX*I_NAI_Hx3yz_Px_ab;
  Double I_NAI_Hx2y2z_D2x_ab = I_NAI_I2x2y2z_Px_ab+ABX*I_NAI_Hx2y2z_Px_ab;
  Double I_NAI_Hxy3z_D2x_ab = I_NAI_I2xy3z_Px_ab+ABX*I_NAI_Hxy3z_Px_ab;
  Double I_NAI_Hx4z_D2x_ab = I_NAI_I2x4z_Px_ab+ABX*I_NAI_Hx4z_Px_ab;
  Double I_NAI_H5y_D2x_ab = I_NAI_Ix5y_Px_ab+ABX*I_NAI_H5y_Px_ab;
  Double I_NAI_H4yz_D2x_ab = I_NAI_Ix4yz_Px_ab+ABX*I_NAI_H4yz_Px_ab;
  Double I_NAI_H3y2z_D2x_ab = I_NAI_Ix3y2z_Px_ab+ABX*I_NAI_H3y2z_Px_ab;
  Double I_NAI_H2y3z_D2x_ab = I_NAI_Ix2y3z_Px_ab+ABX*I_NAI_H2y3z_Px_ab;
  Double I_NAI_Hy4z_D2x_ab = I_NAI_Ixy4z_Px_ab+ABX*I_NAI_Hy4z_Px_ab;
  Double I_NAI_H5z_D2x_ab = I_NAI_Ix5z_Px_ab+ABX*I_NAI_H5z_Px_ab;
  Double I_NAI_H5x_D2y_ab = I_NAI_I5xy_Py_ab+ABY*I_NAI_H5x_Py_ab;
  Double I_NAI_H4xy_D2y_ab = I_NAI_I4x2y_Py_ab+ABY*I_NAI_H4xy_Py_ab;
  Double I_NAI_H4xz_D2y_ab = I_NAI_I4xyz_Py_ab+ABY*I_NAI_H4xz_Py_ab;
  Double I_NAI_H3x2y_D2y_ab = I_NAI_I3x3y_Py_ab+ABY*I_NAI_H3x2y_Py_ab;
  Double I_NAI_H3xyz_D2y_ab = I_NAI_I3x2yz_Py_ab+ABY*I_NAI_H3xyz_Py_ab;
  Double I_NAI_H3x2z_D2y_ab = I_NAI_I3xy2z_Py_ab+ABY*I_NAI_H3x2z_Py_ab;
  Double I_NAI_H2x3y_D2y_ab = I_NAI_I2x4y_Py_ab+ABY*I_NAI_H2x3y_Py_ab;
  Double I_NAI_H2x2yz_D2y_ab = I_NAI_I2x3yz_Py_ab+ABY*I_NAI_H2x2yz_Py_ab;
  Double I_NAI_H2xy2z_D2y_ab = I_NAI_I2x2y2z_Py_ab+ABY*I_NAI_H2xy2z_Py_ab;
  Double I_NAI_H2x3z_D2y_ab = I_NAI_I2xy3z_Py_ab+ABY*I_NAI_H2x3z_Py_ab;
  Double I_NAI_Hx4y_D2y_ab = I_NAI_Ix5y_Py_ab+ABY*I_NAI_Hx4y_Py_ab;
  Double I_NAI_Hx3yz_D2y_ab = I_NAI_Ix4yz_Py_ab+ABY*I_NAI_Hx3yz_Py_ab;
  Double I_NAI_Hx2y2z_D2y_ab = I_NAI_Ix3y2z_Py_ab+ABY*I_NAI_Hx2y2z_Py_ab;
  Double I_NAI_Hxy3z_D2y_ab = I_NAI_Ix2y3z_Py_ab+ABY*I_NAI_Hxy3z_Py_ab;
  Double I_NAI_Hx4z_D2y_ab = I_NAI_Ixy4z_Py_ab+ABY*I_NAI_Hx4z_Py_ab;
  Double I_NAI_H5y_D2y_ab = I_NAI_I6y_Py_ab+ABY*I_NAI_H5y_Py_ab;
  Double I_NAI_H4yz_D2y_ab = I_NAI_I5yz_Py_ab+ABY*I_NAI_H4yz_Py_ab;
  Double I_NAI_H3y2z_D2y_ab = I_NAI_I4y2z_Py_ab+ABY*I_NAI_H3y2z_Py_ab;
  Double I_NAI_H2y3z_D2y_ab = I_NAI_I3y3z_Py_ab+ABY*I_NAI_H2y3z_Py_ab;
  Double I_NAI_Hy4z_D2y_ab = I_NAI_I2y4z_Py_ab+ABY*I_NAI_Hy4z_Py_ab;
  Double I_NAI_H5z_D2y_ab = I_NAI_Iy5z_Py_ab+ABY*I_NAI_H5z_Py_ab;
  Double I_NAI_H5x_D2z_ab = I_NAI_I5xz_Pz_ab+ABZ*I_NAI_H5x_Pz_ab;
  Double I_NAI_H4xy_D2z_ab = I_NAI_I4xyz_Pz_ab+ABZ*I_NAI_H4xy_Pz_ab;
  Double I_NAI_H4xz_D2z_ab = I_NAI_I4x2z_Pz_ab+ABZ*I_NAI_H4xz_Pz_ab;
  Double I_NAI_H3x2y_D2z_ab = I_NAI_I3x2yz_Pz_ab+ABZ*I_NAI_H3x2y_Pz_ab;
  Double I_NAI_H3xyz_D2z_ab = I_NAI_I3xy2z_Pz_ab+ABZ*I_NAI_H3xyz_Pz_ab;
  Double I_NAI_H3x2z_D2z_ab = I_NAI_I3x3z_Pz_ab+ABZ*I_NAI_H3x2z_Pz_ab;
  Double I_NAI_H2x3y_D2z_ab = I_NAI_I2x3yz_Pz_ab+ABZ*I_NAI_H2x3y_Pz_ab;
  Double I_NAI_H2x2yz_D2z_ab = I_NAI_I2x2y2z_Pz_ab+ABZ*I_NAI_H2x2yz_Pz_ab;
  Double I_NAI_H2xy2z_D2z_ab = I_NAI_I2xy3z_Pz_ab+ABZ*I_NAI_H2xy2z_Pz_ab;
  Double I_NAI_H2x3z_D2z_ab = I_NAI_I2x4z_Pz_ab+ABZ*I_NAI_H2x3z_Pz_ab;
  Double I_NAI_Hx4y_D2z_ab = I_NAI_Ix4yz_Pz_ab+ABZ*I_NAI_Hx4y_Pz_ab;
  Double I_NAI_Hx3yz_D2z_ab = I_NAI_Ix3y2z_Pz_ab+ABZ*I_NAI_Hx3yz_Pz_ab;
  Double I_NAI_Hx2y2z_D2z_ab = I_NAI_Ix2y3z_Pz_ab+ABZ*I_NAI_Hx2y2z_Pz_ab;
  Double I_NAI_Hxy3z_D2z_ab = I_NAI_Ixy4z_Pz_ab+ABZ*I_NAI_Hxy3z_Pz_ab;
  Double I_NAI_Hx4z_D2z_ab = I_NAI_Ix5z_Pz_ab+ABZ*I_NAI_Hx4z_Pz_ab;
  Double I_NAI_H5y_D2z_ab = I_NAI_I5yz_Pz_ab+ABZ*I_NAI_H5y_Pz_ab;
  Double I_NAI_H4yz_D2z_ab = I_NAI_I4y2z_Pz_ab+ABZ*I_NAI_H4yz_Pz_ab;
  Double I_NAI_H3y2z_D2z_ab = I_NAI_I3y3z_Pz_ab+ABZ*I_NAI_H3y2z_Pz_ab;
  Double I_NAI_H2y3z_D2z_ab = I_NAI_I2y4z_Pz_ab+ABZ*I_NAI_H2y3z_Pz_ab;
  Double I_NAI_Hy4z_D2z_ab = I_NAI_Iy5z_Pz_ab+ABZ*I_NAI_Hy4z_Pz_ab;
  Double I_NAI_H5z_D2z_ab = I_NAI_I6z_Pz_ab+ABZ*I_NAI_H5z_Pz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_K_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S_ab
   * RHS shell quartet name: SQ_NAI_K_S_ab
   ************************************************************/
  Double I_NAI_K7x_Px_ab = I_NAI_L8x_S_ab+ABX*I_NAI_K7x_S_ab;
  Double I_NAI_K6xy_Px_ab = I_NAI_L7xy_S_ab+ABX*I_NAI_K6xy_S_ab;
  Double I_NAI_K6xz_Px_ab = I_NAI_L7xz_S_ab+ABX*I_NAI_K6xz_S_ab;
  Double I_NAI_K5x2y_Px_ab = I_NAI_L6x2y_S_ab+ABX*I_NAI_K5x2y_S_ab;
  Double I_NAI_K5xyz_Px_ab = I_NAI_L6xyz_S_ab+ABX*I_NAI_K5xyz_S_ab;
  Double I_NAI_K5x2z_Px_ab = I_NAI_L6x2z_S_ab+ABX*I_NAI_K5x2z_S_ab;
  Double I_NAI_K4x3y_Px_ab = I_NAI_L5x3y_S_ab+ABX*I_NAI_K4x3y_S_ab;
  Double I_NAI_K4x2yz_Px_ab = I_NAI_L5x2yz_S_ab+ABX*I_NAI_K4x2yz_S_ab;
  Double I_NAI_K4xy2z_Px_ab = I_NAI_L5xy2z_S_ab+ABX*I_NAI_K4xy2z_S_ab;
  Double I_NAI_K4x3z_Px_ab = I_NAI_L5x3z_S_ab+ABX*I_NAI_K4x3z_S_ab;
  Double I_NAI_K3x4y_Px_ab = I_NAI_L4x4y_S_ab+ABX*I_NAI_K3x4y_S_ab;
  Double I_NAI_K3x3yz_Px_ab = I_NAI_L4x3yz_S_ab+ABX*I_NAI_K3x3yz_S_ab;
  Double I_NAI_K3x2y2z_Px_ab = I_NAI_L4x2y2z_S_ab+ABX*I_NAI_K3x2y2z_S_ab;
  Double I_NAI_K3xy3z_Px_ab = I_NAI_L4xy3z_S_ab+ABX*I_NAI_K3xy3z_S_ab;
  Double I_NAI_K3x4z_Px_ab = I_NAI_L4x4z_S_ab+ABX*I_NAI_K3x4z_S_ab;
  Double I_NAI_K2x5y_Px_ab = I_NAI_L3x5y_S_ab+ABX*I_NAI_K2x5y_S_ab;
  Double I_NAI_K2x4yz_Px_ab = I_NAI_L3x4yz_S_ab+ABX*I_NAI_K2x4yz_S_ab;
  Double I_NAI_K2x3y2z_Px_ab = I_NAI_L3x3y2z_S_ab+ABX*I_NAI_K2x3y2z_S_ab;
  Double I_NAI_K2x2y3z_Px_ab = I_NAI_L3x2y3z_S_ab+ABX*I_NAI_K2x2y3z_S_ab;
  Double I_NAI_K2xy4z_Px_ab = I_NAI_L3xy4z_S_ab+ABX*I_NAI_K2xy4z_S_ab;
  Double I_NAI_K2x5z_Px_ab = I_NAI_L3x5z_S_ab+ABX*I_NAI_K2x5z_S_ab;
  Double I_NAI_Kx6y_Px_ab = I_NAI_L2x6y_S_ab+ABX*I_NAI_Kx6y_S_ab;
  Double I_NAI_Kx5yz_Px_ab = I_NAI_L2x5yz_S_ab+ABX*I_NAI_Kx5yz_S_ab;
  Double I_NAI_Kx4y2z_Px_ab = I_NAI_L2x4y2z_S_ab+ABX*I_NAI_Kx4y2z_S_ab;
  Double I_NAI_Kx3y3z_Px_ab = I_NAI_L2x3y3z_S_ab+ABX*I_NAI_Kx3y3z_S_ab;
  Double I_NAI_Kx2y4z_Px_ab = I_NAI_L2x2y4z_S_ab+ABX*I_NAI_Kx2y4z_S_ab;
  Double I_NAI_Kxy5z_Px_ab = I_NAI_L2xy5z_S_ab+ABX*I_NAI_Kxy5z_S_ab;
  Double I_NAI_Kx6z_Px_ab = I_NAI_L2x6z_S_ab+ABX*I_NAI_Kx6z_S_ab;
  Double I_NAI_K7y_Px_ab = I_NAI_Lx7y_S_ab+ABX*I_NAI_K7y_S_ab;
  Double I_NAI_K6yz_Px_ab = I_NAI_Lx6yz_S_ab+ABX*I_NAI_K6yz_S_ab;
  Double I_NAI_K5y2z_Px_ab = I_NAI_Lx5y2z_S_ab+ABX*I_NAI_K5y2z_S_ab;
  Double I_NAI_K4y3z_Px_ab = I_NAI_Lx4y3z_S_ab+ABX*I_NAI_K4y3z_S_ab;
  Double I_NAI_K3y4z_Px_ab = I_NAI_Lx3y4z_S_ab+ABX*I_NAI_K3y4z_S_ab;
  Double I_NAI_K2y5z_Px_ab = I_NAI_Lx2y5z_S_ab+ABX*I_NAI_K2y5z_S_ab;
  Double I_NAI_Ky6z_Px_ab = I_NAI_Lxy6z_S_ab+ABX*I_NAI_Ky6z_S_ab;
  Double I_NAI_K7z_Px_ab = I_NAI_Lx7z_S_ab+ABX*I_NAI_K7z_S_ab;
  Double I_NAI_K6xy_Py_ab = I_NAI_L6x2y_S_ab+ABY*I_NAI_K6xy_S_ab;
  Double I_NAI_K6xz_Py_ab = I_NAI_L6xyz_S_ab+ABY*I_NAI_K6xz_S_ab;
  Double I_NAI_K5x2y_Py_ab = I_NAI_L5x3y_S_ab+ABY*I_NAI_K5x2y_S_ab;
  Double I_NAI_K5xyz_Py_ab = I_NAI_L5x2yz_S_ab+ABY*I_NAI_K5xyz_S_ab;
  Double I_NAI_K5x2z_Py_ab = I_NAI_L5xy2z_S_ab+ABY*I_NAI_K5x2z_S_ab;
  Double I_NAI_K4x3y_Py_ab = I_NAI_L4x4y_S_ab+ABY*I_NAI_K4x3y_S_ab;
  Double I_NAI_K4x2yz_Py_ab = I_NAI_L4x3yz_S_ab+ABY*I_NAI_K4x2yz_S_ab;
  Double I_NAI_K4xy2z_Py_ab = I_NAI_L4x2y2z_S_ab+ABY*I_NAI_K4xy2z_S_ab;
  Double I_NAI_K4x3z_Py_ab = I_NAI_L4xy3z_S_ab+ABY*I_NAI_K4x3z_S_ab;
  Double I_NAI_K3x4y_Py_ab = I_NAI_L3x5y_S_ab+ABY*I_NAI_K3x4y_S_ab;
  Double I_NAI_K3x3yz_Py_ab = I_NAI_L3x4yz_S_ab+ABY*I_NAI_K3x3yz_S_ab;
  Double I_NAI_K3x2y2z_Py_ab = I_NAI_L3x3y2z_S_ab+ABY*I_NAI_K3x2y2z_S_ab;
  Double I_NAI_K3xy3z_Py_ab = I_NAI_L3x2y3z_S_ab+ABY*I_NAI_K3xy3z_S_ab;
  Double I_NAI_K3x4z_Py_ab = I_NAI_L3xy4z_S_ab+ABY*I_NAI_K3x4z_S_ab;
  Double I_NAI_K2x5y_Py_ab = I_NAI_L2x6y_S_ab+ABY*I_NAI_K2x5y_S_ab;
  Double I_NAI_K2x4yz_Py_ab = I_NAI_L2x5yz_S_ab+ABY*I_NAI_K2x4yz_S_ab;
  Double I_NAI_K2x3y2z_Py_ab = I_NAI_L2x4y2z_S_ab+ABY*I_NAI_K2x3y2z_S_ab;
  Double I_NAI_K2x2y3z_Py_ab = I_NAI_L2x3y3z_S_ab+ABY*I_NAI_K2x2y3z_S_ab;
  Double I_NAI_K2xy4z_Py_ab = I_NAI_L2x2y4z_S_ab+ABY*I_NAI_K2xy4z_S_ab;
  Double I_NAI_K2x5z_Py_ab = I_NAI_L2xy5z_S_ab+ABY*I_NAI_K2x5z_S_ab;
  Double I_NAI_Kx6y_Py_ab = I_NAI_Lx7y_S_ab+ABY*I_NAI_Kx6y_S_ab;
  Double I_NAI_Kx5yz_Py_ab = I_NAI_Lx6yz_S_ab+ABY*I_NAI_Kx5yz_S_ab;
  Double I_NAI_Kx4y2z_Py_ab = I_NAI_Lx5y2z_S_ab+ABY*I_NAI_Kx4y2z_S_ab;
  Double I_NAI_Kx3y3z_Py_ab = I_NAI_Lx4y3z_S_ab+ABY*I_NAI_Kx3y3z_S_ab;
  Double I_NAI_Kx2y4z_Py_ab = I_NAI_Lx3y4z_S_ab+ABY*I_NAI_Kx2y4z_S_ab;
  Double I_NAI_Kxy5z_Py_ab = I_NAI_Lx2y5z_S_ab+ABY*I_NAI_Kxy5z_S_ab;
  Double I_NAI_Kx6z_Py_ab = I_NAI_Lxy6z_S_ab+ABY*I_NAI_Kx6z_S_ab;
  Double I_NAI_K7y_Py_ab = I_NAI_L8y_S_ab+ABY*I_NAI_K7y_S_ab;
  Double I_NAI_K6yz_Py_ab = I_NAI_L7yz_S_ab+ABY*I_NAI_K6yz_S_ab;
  Double I_NAI_K5y2z_Py_ab = I_NAI_L6y2z_S_ab+ABY*I_NAI_K5y2z_S_ab;
  Double I_NAI_K4y3z_Py_ab = I_NAI_L5y3z_S_ab+ABY*I_NAI_K4y3z_S_ab;
  Double I_NAI_K3y4z_Py_ab = I_NAI_L4y4z_S_ab+ABY*I_NAI_K3y4z_S_ab;
  Double I_NAI_K2y5z_Py_ab = I_NAI_L3y5z_S_ab+ABY*I_NAI_K2y5z_S_ab;
  Double I_NAI_Ky6z_Py_ab = I_NAI_L2y6z_S_ab+ABY*I_NAI_Ky6z_S_ab;
  Double I_NAI_K7z_Py_ab = I_NAI_Ly7z_S_ab+ABY*I_NAI_K7z_S_ab;
  Double I_NAI_K6xy_Pz_ab = I_NAI_L6xyz_S_ab+ABZ*I_NAI_K6xy_S_ab;
  Double I_NAI_K6xz_Pz_ab = I_NAI_L6x2z_S_ab+ABZ*I_NAI_K6xz_S_ab;
  Double I_NAI_K5x2y_Pz_ab = I_NAI_L5x2yz_S_ab+ABZ*I_NAI_K5x2y_S_ab;
  Double I_NAI_K5xyz_Pz_ab = I_NAI_L5xy2z_S_ab+ABZ*I_NAI_K5xyz_S_ab;
  Double I_NAI_K5x2z_Pz_ab = I_NAI_L5x3z_S_ab+ABZ*I_NAI_K5x2z_S_ab;
  Double I_NAI_K4x3y_Pz_ab = I_NAI_L4x3yz_S_ab+ABZ*I_NAI_K4x3y_S_ab;
  Double I_NAI_K4x2yz_Pz_ab = I_NAI_L4x2y2z_S_ab+ABZ*I_NAI_K4x2yz_S_ab;
  Double I_NAI_K4xy2z_Pz_ab = I_NAI_L4xy3z_S_ab+ABZ*I_NAI_K4xy2z_S_ab;
  Double I_NAI_K4x3z_Pz_ab = I_NAI_L4x4z_S_ab+ABZ*I_NAI_K4x3z_S_ab;
  Double I_NAI_K3x4y_Pz_ab = I_NAI_L3x4yz_S_ab+ABZ*I_NAI_K3x4y_S_ab;
  Double I_NAI_K3x3yz_Pz_ab = I_NAI_L3x3y2z_S_ab+ABZ*I_NAI_K3x3yz_S_ab;
  Double I_NAI_K3x2y2z_Pz_ab = I_NAI_L3x2y3z_S_ab+ABZ*I_NAI_K3x2y2z_S_ab;
  Double I_NAI_K3xy3z_Pz_ab = I_NAI_L3xy4z_S_ab+ABZ*I_NAI_K3xy3z_S_ab;
  Double I_NAI_K3x4z_Pz_ab = I_NAI_L3x5z_S_ab+ABZ*I_NAI_K3x4z_S_ab;
  Double I_NAI_K2x5y_Pz_ab = I_NAI_L2x5yz_S_ab+ABZ*I_NAI_K2x5y_S_ab;
  Double I_NAI_K2x4yz_Pz_ab = I_NAI_L2x4y2z_S_ab+ABZ*I_NAI_K2x4yz_S_ab;
  Double I_NAI_K2x3y2z_Pz_ab = I_NAI_L2x3y3z_S_ab+ABZ*I_NAI_K2x3y2z_S_ab;
  Double I_NAI_K2x2y3z_Pz_ab = I_NAI_L2x2y4z_S_ab+ABZ*I_NAI_K2x2y3z_S_ab;
  Double I_NAI_K2xy4z_Pz_ab = I_NAI_L2xy5z_S_ab+ABZ*I_NAI_K2xy4z_S_ab;
  Double I_NAI_K2x5z_Pz_ab = I_NAI_L2x6z_S_ab+ABZ*I_NAI_K2x5z_S_ab;
  Double I_NAI_Kx6y_Pz_ab = I_NAI_Lx6yz_S_ab+ABZ*I_NAI_Kx6y_S_ab;
  Double I_NAI_Kx5yz_Pz_ab = I_NAI_Lx5y2z_S_ab+ABZ*I_NAI_Kx5yz_S_ab;
  Double I_NAI_Kx4y2z_Pz_ab = I_NAI_Lx4y3z_S_ab+ABZ*I_NAI_Kx4y2z_S_ab;
  Double I_NAI_Kx3y3z_Pz_ab = I_NAI_Lx3y4z_S_ab+ABZ*I_NAI_Kx3y3z_S_ab;
  Double I_NAI_Kx2y4z_Pz_ab = I_NAI_Lx2y5z_S_ab+ABZ*I_NAI_Kx2y4z_S_ab;
  Double I_NAI_Kxy5z_Pz_ab = I_NAI_Lxy6z_S_ab+ABZ*I_NAI_Kxy5z_S_ab;
  Double I_NAI_Kx6z_Pz_ab = I_NAI_Lx7z_S_ab+ABZ*I_NAI_Kx6z_S_ab;
  Double I_NAI_K6yz_Pz_ab = I_NAI_L6y2z_S_ab+ABZ*I_NAI_K6yz_S_ab;
  Double I_NAI_K5y2z_Pz_ab = I_NAI_L5y3z_S_ab+ABZ*I_NAI_K5y2z_S_ab;
  Double I_NAI_K4y3z_Pz_ab = I_NAI_L4y4z_S_ab+ABZ*I_NAI_K4y3z_S_ab;
  Double I_NAI_K3y4z_Pz_ab = I_NAI_L3y5z_S_ab+ABZ*I_NAI_K3y4z_S_ab;
  Double I_NAI_K2y5z_Pz_ab = I_NAI_L2y6z_S_ab+ABZ*I_NAI_K2y5z_S_ab;
  Double I_NAI_Ky6z_Pz_ab = I_NAI_Ly7z_S_ab+ABZ*I_NAI_Ky6z_S_ab;
  Double I_NAI_K7z_Pz_ab = I_NAI_L8z_S_ab+ABZ*I_NAI_K7z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 84 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_ab
   * RHS shell quartet name: SQ_NAI_I_P_ab
   ************************************************************/
  Double I_NAI_I6x_D2x_ab = I_NAI_K7x_Px_ab+ABX*I_NAI_I6x_Px_ab;
  Double I_NAI_I5xy_D2x_ab = I_NAI_K6xy_Px_ab+ABX*I_NAI_I5xy_Px_ab;
  Double I_NAI_I5xz_D2x_ab = I_NAI_K6xz_Px_ab+ABX*I_NAI_I5xz_Px_ab;
  Double I_NAI_I4x2y_D2x_ab = I_NAI_K5x2y_Px_ab+ABX*I_NAI_I4x2y_Px_ab;
  Double I_NAI_I4xyz_D2x_ab = I_NAI_K5xyz_Px_ab+ABX*I_NAI_I4xyz_Px_ab;
  Double I_NAI_I4x2z_D2x_ab = I_NAI_K5x2z_Px_ab+ABX*I_NAI_I4x2z_Px_ab;
  Double I_NAI_I3x3y_D2x_ab = I_NAI_K4x3y_Px_ab+ABX*I_NAI_I3x3y_Px_ab;
  Double I_NAI_I3x2yz_D2x_ab = I_NAI_K4x2yz_Px_ab+ABX*I_NAI_I3x2yz_Px_ab;
  Double I_NAI_I3xy2z_D2x_ab = I_NAI_K4xy2z_Px_ab+ABX*I_NAI_I3xy2z_Px_ab;
  Double I_NAI_I3x3z_D2x_ab = I_NAI_K4x3z_Px_ab+ABX*I_NAI_I3x3z_Px_ab;
  Double I_NAI_I2x4y_D2x_ab = I_NAI_K3x4y_Px_ab+ABX*I_NAI_I2x4y_Px_ab;
  Double I_NAI_I2x3yz_D2x_ab = I_NAI_K3x3yz_Px_ab+ABX*I_NAI_I2x3yz_Px_ab;
  Double I_NAI_I2x2y2z_D2x_ab = I_NAI_K3x2y2z_Px_ab+ABX*I_NAI_I2x2y2z_Px_ab;
  Double I_NAI_I2xy3z_D2x_ab = I_NAI_K3xy3z_Px_ab+ABX*I_NAI_I2xy3z_Px_ab;
  Double I_NAI_I2x4z_D2x_ab = I_NAI_K3x4z_Px_ab+ABX*I_NAI_I2x4z_Px_ab;
  Double I_NAI_Ix5y_D2x_ab = I_NAI_K2x5y_Px_ab+ABX*I_NAI_Ix5y_Px_ab;
  Double I_NAI_Ix4yz_D2x_ab = I_NAI_K2x4yz_Px_ab+ABX*I_NAI_Ix4yz_Px_ab;
  Double I_NAI_Ix3y2z_D2x_ab = I_NAI_K2x3y2z_Px_ab+ABX*I_NAI_Ix3y2z_Px_ab;
  Double I_NAI_Ix2y3z_D2x_ab = I_NAI_K2x2y3z_Px_ab+ABX*I_NAI_Ix2y3z_Px_ab;
  Double I_NAI_Ixy4z_D2x_ab = I_NAI_K2xy4z_Px_ab+ABX*I_NAI_Ixy4z_Px_ab;
  Double I_NAI_Ix5z_D2x_ab = I_NAI_K2x5z_Px_ab+ABX*I_NAI_Ix5z_Px_ab;
  Double I_NAI_I6y_D2x_ab = I_NAI_Kx6y_Px_ab+ABX*I_NAI_I6y_Px_ab;
  Double I_NAI_I5yz_D2x_ab = I_NAI_Kx5yz_Px_ab+ABX*I_NAI_I5yz_Px_ab;
  Double I_NAI_I4y2z_D2x_ab = I_NAI_Kx4y2z_Px_ab+ABX*I_NAI_I4y2z_Px_ab;
  Double I_NAI_I3y3z_D2x_ab = I_NAI_Kx3y3z_Px_ab+ABX*I_NAI_I3y3z_Px_ab;
  Double I_NAI_I2y4z_D2x_ab = I_NAI_Kx2y4z_Px_ab+ABX*I_NAI_I2y4z_Px_ab;
  Double I_NAI_Iy5z_D2x_ab = I_NAI_Kxy5z_Px_ab+ABX*I_NAI_Iy5z_Px_ab;
  Double I_NAI_I6z_D2x_ab = I_NAI_Kx6z_Px_ab+ABX*I_NAI_I6z_Px_ab;
  Double I_NAI_I6x_D2y_ab = I_NAI_K6xy_Py_ab+ABY*I_NAI_I6x_Py_ab;
  Double I_NAI_I5xy_D2y_ab = I_NAI_K5x2y_Py_ab+ABY*I_NAI_I5xy_Py_ab;
  Double I_NAI_I5xz_D2y_ab = I_NAI_K5xyz_Py_ab+ABY*I_NAI_I5xz_Py_ab;
  Double I_NAI_I4x2y_D2y_ab = I_NAI_K4x3y_Py_ab+ABY*I_NAI_I4x2y_Py_ab;
  Double I_NAI_I4xyz_D2y_ab = I_NAI_K4x2yz_Py_ab+ABY*I_NAI_I4xyz_Py_ab;
  Double I_NAI_I4x2z_D2y_ab = I_NAI_K4xy2z_Py_ab+ABY*I_NAI_I4x2z_Py_ab;
  Double I_NAI_I3x3y_D2y_ab = I_NAI_K3x4y_Py_ab+ABY*I_NAI_I3x3y_Py_ab;
  Double I_NAI_I3x2yz_D2y_ab = I_NAI_K3x3yz_Py_ab+ABY*I_NAI_I3x2yz_Py_ab;
  Double I_NAI_I3xy2z_D2y_ab = I_NAI_K3x2y2z_Py_ab+ABY*I_NAI_I3xy2z_Py_ab;
  Double I_NAI_I3x3z_D2y_ab = I_NAI_K3xy3z_Py_ab+ABY*I_NAI_I3x3z_Py_ab;
  Double I_NAI_I2x4y_D2y_ab = I_NAI_K2x5y_Py_ab+ABY*I_NAI_I2x4y_Py_ab;
  Double I_NAI_I2x3yz_D2y_ab = I_NAI_K2x4yz_Py_ab+ABY*I_NAI_I2x3yz_Py_ab;
  Double I_NAI_I2x2y2z_D2y_ab = I_NAI_K2x3y2z_Py_ab+ABY*I_NAI_I2x2y2z_Py_ab;
  Double I_NAI_I2xy3z_D2y_ab = I_NAI_K2x2y3z_Py_ab+ABY*I_NAI_I2xy3z_Py_ab;
  Double I_NAI_I2x4z_D2y_ab = I_NAI_K2xy4z_Py_ab+ABY*I_NAI_I2x4z_Py_ab;
  Double I_NAI_Ix5y_D2y_ab = I_NAI_Kx6y_Py_ab+ABY*I_NAI_Ix5y_Py_ab;
  Double I_NAI_Ix4yz_D2y_ab = I_NAI_Kx5yz_Py_ab+ABY*I_NAI_Ix4yz_Py_ab;
  Double I_NAI_Ix3y2z_D2y_ab = I_NAI_Kx4y2z_Py_ab+ABY*I_NAI_Ix3y2z_Py_ab;
  Double I_NAI_Ix2y3z_D2y_ab = I_NAI_Kx3y3z_Py_ab+ABY*I_NAI_Ix2y3z_Py_ab;
  Double I_NAI_Ixy4z_D2y_ab = I_NAI_Kx2y4z_Py_ab+ABY*I_NAI_Ixy4z_Py_ab;
  Double I_NAI_Ix5z_D2y_ab = I_NAI_Kxy5z_Py_ab+ABY*I_NAI_Ix5z_Py_ab;
  Double I_NAI_I6y_D2y_ab = I_NAI_K7y_Py_ab+ABY*I_NAI_I6y_Py_ab;
  Double I_NAI_I5yz_D2y_ab = I_NAI_K6yz_Py_ab+ABY*I_NAI_I5yz_Py_ab;
  Double I_NAI_I4y2z_D2y_ab = I_NAI_K5y2z_Py_ab+ABY*I_NAI_I4y2z_Py_ab;
  Double I_NAI_I3y3z_D2y_ab = I_NAI_K4y3z_Py_ab+ABY*I_NAI_I3y3z_Py_ab;
  Double I_NAI_I2y4z_D2y_ab = I_NAI_K3y4z_Py_ab+ABY*I_NAI_I2y4z_Py_ab;
  Double I_NAI_Iy5z_D2y_ab = I_NAI_K2y5z_Py_ab+ABY*I_NAI_Iy5z_Py_ab;
  Double I_NAI_I6z_D2y_ab = I_NAI_Ky6z_Py_ab+ABY*I_NAI_I6z_Py_ab;
  Double I_NAI_I6x_D2z_ab = I_NAI_K6xz_Pz_ab+ABZ*I_NAI_I6x_Pz_ab;
  Double I_NAI_I5xy_D2z_ab = I_NAI_K5xyz_Pz_ab+ABZ*I_NAI_I5xy_Pz_ab;
  Double I_NAI_I5xz_D2z_ab = I_NAI_K5x2z_Pz_ab+ABZ*I_NAI_I5xz_Pz_ab;
  Double I_NAI_I4x2y_D2z_ab = I_NAI_K4x2yz_Pz_ab+ABZ*I_NAI_I4x2y_Pz_ab;
  Double I_NAI_I4xyz_D2z_ab = I_NAI_K4xy2z_Pz_ab+ABZ*I_NAI_I4xyz_Pz_ab;
  Double I_NAI_I4x2z_D2z_ab = I_NAI_K4x3z_Pz_ab+ABZ*I_NAI_I4x2z_Pz_ab;
  Double I_NAI_I3x3y_D2z_ab = I_NAI_K3x3yz_Pz_ab+ABZ*I_NAI_I3x3y_Pz_ab;
  Double I_NAI_I3x2yz_D2z_ab = I_NAI_K3x2y2z_Pz_ab+ABZ*I_NAI_I3x2yz_Pz_ab;
  Double I_NAI_I3xy2z_D2z_ab = I_NAI_K3xy3z_Pz_ab+ABZ*I_NAI_I3xy2z_Pz_ab;
  Double I_NAI_I3x3z_D2z_ab = I_NAI_K3x4z_Pz_ab+ABZ*I_NAI_I3x3z_Pz_ab;
  Double I_NAI_I2x4y_D2z_ab = I_NAI_K2x4yz_Pz_ab+ABZ*I_NAI_I2x4y_Pz_ab;
  Double I_NAI_I2x3yz_D2z_ab = I_NAI_K2x3y2z_Pz_ab+ABZ*I_NAI_I2x3yz_Pz_ab;
  Double I_NAI_I2x2y2z_D2z_ab = I_NAI_K2x2y3z_Pz_ab+ABZ*I_NAI_I2x2y2z_Pz_ab;
  Double I_NAI_I2xy3z_D2z_ab = I_NAI_K2xy4z_Pz_ab+ABZ*I_NAI_I2xy3z_Pz_ab;
  Double I_NAI_I2x4z_D2z_ab = I_NAI_K2x5z_Pz_ab+ABZ*I_NAI_I2x4z_Pz_ab;
  Double I_NAI_Ix5y_D2z_ab = I_NAI_Kx5yz_Pz_ab+ABZ*I_NAI_Ix5y_Pz_ab;
  Double I_NAI_Ix4yz_D2z_ab = I_NAI_Kx4y2z_Pz_ab+ABZ*I_NAI_Ix4yz_Pz_ab;
  Double I_NAI_Ix3y2z_D2z_ab = I_NAI_Kx3y3z_Pz_ab+ABZ*I_NAI_Ix3y2z_Pz_ab;
  Double I_NAI_Ix2y3z_D2z_ab = I_NAI_Kx2y4z_Pz_ab+ABZ*I_NAI_Ix2y3z_Pz_ab;
  Double I_NAI_Ixy4z_D2z_ab = I_NAI_Kxy5z_Pz_ab+ABZ*I_NAI_Ixy4z_Pz_ab;
  Double I_NAI_Ix5z_D2z_ab = I_NAI_Kx6z_Pz_ab+ABZ*I_NAI_Ix5z_Pz_ab;
  Double I_NAI_I6y_D2z_ab = I_NAI_K6yz_Pz_ab+ABZ*I_NAI_I6y_Pz_ab;
  Double I_NAI_I5yz_D2z_ab = I_NAI_K5y2z_Pz_ab+ABZ*I_NAI_I5yz_Pz_ab;
  Double I_NAI_I4y2z_D2z_ab = I_NAI_K4y3z_Pz_ab+ABZ*I_NAI_I4y2z_Pz_ab;
  Double I_NAI_I3y3z_D2z_ab = I_NAI_K3y4z_Pz_ab+ABZ*I_NAI_I3y3z_Pz_ab;
  Double I_NAI_I2y4z_D2z_ab = I_NAI_K2y5z_Pz_ab+ABZ*I_NAI_I2y4z_Pz_ab;
  Double I_NAI_Iy5z_D2z_ab = I_NAI_Ky6z_Pz_ab+ABZ*I_NAI_Iy5z_Pz_ab;
  Double I_NAI_I6z_D2z_ab = I_NAI_K7z_Pz_ab+ABZ*I_NAI_I6z_Pz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 42 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_H_D_ab
   ************************************************************/
  Double I_NAI_H5x_F3x_ab = I_NAI_I6x_D2x_ab+ABX*I_NAI_H5x_D2x_ab;
  Double I_NAI_H4xy_F3x_ab = I_NAI_I5xy_D2x_ab+ABX*I_NAI_H4xy_D2x_ab;
  Double I_NAI_H4xz_F3x_ab = I_NAI_I5xz_D2x_ab+ABX*I_NAI_H4xz_D2x_ab;
  Double I_NAI_H3x2y_F3x_ab = I_NAI_I4x2y_D2x_ab+ABX*I_NAI_H3x2y_D2x_ab;
  Double I_NAI_H3xyz_F3x_ab = I_NAI_I4xyz_D2x_ab+ABX*I_NAI_H3xyz_D2x_ab;
  Double I_NAI_H3x2z_F3x_ab = I_NAI_I4x2z_D2x_ab+ABX*I_NAI_H3x2z_D2x_ab;
  Double I_NAI_H2x3y_F3x_ab = I_NAI_I3x3y_D2x_ab+ABX*I_NAI_H2x3y_D2x_ab;
  Double I_NAI_H2x2yz_F3x_ab = I_NAI_I3x2yz_D2x_ab+ABX*I_NAI_H2x2yz_D2x_ab;
  Double I_NAI_H2xy2z_F3x_ab = I_NAI_I3xy2z_D2x_ab+ABX*I_NAI_H2xy2z_D2x_ab;
  Double I_NAI_H2x3z_F3x_ab = I_NAI_I3x3z_D2x_ab+ABX*I_NAI_H2x3z_D2x_ab;
  Double I_NAI_Hx4y_F3x_ab = I_NAI_I2x4y_D2x_ab+ABX*I_NAI_Hx4y_D2x_ab;
  Double I_NAI_Hx3yz_F3x_ab = I_NAI_I2x3yz_D2x_ab+ABX*I_NAI_Hx3yz_D2x_ab;
  Double I_NAI_Hx2y2z_F3x_ab = I_NAI_I2x2y2z_D2x_ab+ABX*I_NAI_Hx2y2z_D2x_ab;
  Double I_NAI_Hxy3z_F3x_ab = I_NAI_I2xy3z_D2x_ab+ABX*I_NAI_Hxy3z_D2x_ab;
  Double I_NAI_Hx4z_F3x_ab = I_NAI_I2x4z_D2x_ab+ABX*I_NAI_Hx4z_D2x_ab;
  Double I_NAI_H5y_F3x_ab = I_NAI_Ix5y_D2x_ab+ABX*I_NAI_H5y_D2x_ab;
  Double I_NAI_H4yz_F3x_ab = I_NAI_Ix4yz_D2x_ab+ABX*I_NAI_H4yz_D2x_ab;
  Double I_NAI_H3y2z_F3x_ab = I_NAI_Ix3y2z_D2x_ab+ABX*I_NAI_H3y2z_D2x_ab;
  Double I_NAI_H2y3z_F3x_ab = I_NAI_Ix2y3z_D2x_ab+ABX*I_NAI_H2y3z_D2x_ab;
  Double I_NAI_Hy4z_F3x_ab = I_NAI_Ixy4z_D2x_ab+ABX*I_NAI_Hy4z_D2x_ab;
  Double I_NAI_H5z_F3x_ab = I_NAI_Ix5z_D2x_ab+ABX*I_NAI_H5z_D2x_ab;
  Double I_NAI_H5x_F2xy_ab = I_NAI_I5xy_D2x_ab+ABY*I_NAI_H5x_D2x_ab;
  Double I_NAI_H4xy_F2xy_ab = I_NAI_I4x2y_D2x_ab+ABY*I_NAI_H4xy_D2x_ab;
  Double I_NAI_H4xz_F2xy_ab = I_NAI_I4xyz_D2x_ab+ABY*I_NAI_H4xz_D2x_ab;
  Double I_NAI_H3x2y_F2xy_ab = I_NAI_I3x3y_D2x_ab+ABY*I_NAI_H3x2y_D2x_ab;
  Double I_NAI_H3xyz_F2xy_ab = I_NAI_I3x2yz_D2x_ab+ABY*I_NAI_H3xyz_D2x_ab;
  Double I_NAI_H3x2z_F2xy_ab = I_NAI_I3xy2z_D2x_ab+ABY*I_NAI_H3x2z_D2x_ab;
  Double I_NAI_H2x3y_F2xy_ab = I_NAI_I2x4y_D2x_ab+ABY*I_NAI_H2x3y_D2x_ab;
  Double I_NAI_H2x2yz_F2xy_ab = I_NAI_I2x3yz_D2x_ab+ABY*I_NAI_H2x2yz_D2x_ab;
  Double I_NAI_H2xy2z_F2xy_ab = I_NAI_I2x2y2z_D2x_ab+ABY*I_NAI_H2xy2z_D2x_ab;
  Double I_NAI_H2x3z_F2xy_ab = I_NAI_I2xy3z_D2x_ab+ABY*I_NAI_H2x3z_D2x_ab;
  Double I_NAI_Hx4y_F2xy_ab = I_NAI_Ix5y_D2x_ab+ABY*I_NAI_Hx4y_D2x_ab;
  Double I_NAI_Hx3yz_F2xy_ab = I_NAI_Ix4yz_D2x_ab+ABY*I_NAI_Hx3yz_D2x_ab;
  Double I_NAI_Hx2y2z_F2xy_ab = I_NAI_Ix3y2z_D2x_ab+ABY*I_NAI_Hx2y2z_D2x_ab;
  Double I_NAI_Hxy3z_F2xy_ab = I_NAI_Ix2y3z_D2x_ab+ABY*I_NAI_Hxy3z_D2x_ab;
  Double I_NAI_Hx4z_F2xy_ab = I_NAI_Ixy4z_D2x_ab+ABY*I_NAI_Hx4z_D2x_ab;
  Double I_NAI_H5y_F2xy_ab = I_NAI_I6y_D2x_ab+ABY*I_NAI_H5y_D2x_ab;
  Double I_NAI_H4yz_F2xy_ab = I_NAI_I5yz_D2x_ab+ABY*I_NAI_H4yz_D2x_ab;
  Double I_NAI_H3y2z_F2xy_ab = I_NAI_I4y2z_D2x_ab+ABY*I_NAI_H3y2z_D2x_ab;
  Double I_NAI_H2y3z_F2xy_ab = I_NAI_I3y3z_D2x_ab+ABY*I_NAI_H2y3z_D2x_ab;
  Double I_NAI_Hy4z_F2xy_ab = I_NAI_I2y4z_D2x_ab+ABY*I_NAI_Hy4z_D2x_ab;
  Double I_NAI_H5z_F2xy_ab = I_NAI_Iy5z_D2x_ab+ABY*I_NAI_H5z_D2x_ab;
  Double I_NAI_H5x_F2xz_ab = I_NAI_I5xz_D2x_ab+ABZ*I_NAI_H5x_D2x_ab;
  Double I_NAI_H4xy_F2xz_ab = I_NAI_I4xyz_D2x_ab+ABZ*I_NAI_H4xy_D2x_ab;
  Double I_NAI_H4xz_F2xz_ab = I_NAI_I4x2z_D2x_ab+ABZ*I_NAI_H4xz_D2x_ab;
  Double I_NAI_H3x2y_F2xz_ab = I_NAI_I3x2yz_D2x_ab+ABZ*I_NAI_H3x2y_D2x_ab;
  Double I_NAI_H3xyz_F2xz_ab = I_NAI_I3xy2z_D2x_ab+ABZ*I_NAI_H3xyz_D2x_ab;
  Double I_NAI_H3x2z_F2xz_ab = I_NAI_I3x3z_D2x_ab+ABZ*I_NAI_H3x2z_D2x_ab;
  Double I_NAI_H2x3y_F2xz_ab = I_NAI_I2x3yz_D2x_ab+ABZ*I_NAI_H2x3y_D2x_ab;
  Double I_NAI_H2x2yz_F2xz_ab = I_NAI_I2x2y2z_D2x_ab+ABZ*I_NAI_H2x2yz_D2x_ab;
  Double I_NAI_H2xy2z_F2xz_ab = I_NAI_I2xy3z_D2x_ab+ABZ*I_NAI_H2xy2z_D2x_ab;
  Double I_NAI_H2x3z_F2xz_ab = I_NAI_I2x4z_D2x_ab+ABZ*I_NAI_H2x3z_D2x_ab;
  Double I_NAI_Hx4y_F2xz_ab = I_NAI_Ix4yz_D2x_ab+ABZ*I_NAI_Hx4y_D2x_ab;
  Double I_NAI_Hx3yz_F2xz_ab = I_NAI_Ix3y2z_D2x_ab+ABZ*I_NAI_Hx3yz_D2x_ab;
  Double I_NAI_Hx2y2z_F2xz_ab = I_NAI_Ix2y3z_D2x_ab+ABZ*I_NAI_Hx2y2z_D2x_ab;
  Double I_NAI_Hxy3z_F2xz_ab = I_NAI_Ixy4z_D2x_ab+ABZ*I_NAI_Hxy3z_D2x_ab;
  Double I_NAI_Hx4z_F2xz_ab = I_NAI_Ix5z_D2x_ab+ABZ*I_NAI_Hx4z_D2x_ab;
  Double I_NAI_H5y_F2xz_ab = I_NAI_I5yz_D2x_ab+ABZ*I_NAI_H5y_D2x_ab;
  Double I_NAI_H4yz_F2xz_ab = I_NAI_I4y2z_D2x_ab+ABZ*I_NAI_H4yz_D2x_ab;
  Double I_NAI_H3y2z_F2xz_ab = I_NAI_I3y3z_D2x_ab+ABZ*I_NAI_H3y2z_D2x_ab;
  Double I_NAI_H2y3z_F2xz_ab = I_NAI_I2y4z_D2x_ab+ABZ*I_NAI_H2y3z_D2x_ab;
  Double I_NAI_Hy4z_F2xz_ab = I_NAI_Iy5z_D2x_ab+ABZ*I_NAI_Hy4z_D2x_ab;
  Double I_NAI_H5z_F2xz_ab = I_NAI_I6z_D2x_ab+ABZ*I_NAI_H5z_D2x_ab;
  Double I_NAI_H5x_Fx2y_ab = I_NAI_I6x_D2y_ab+ABX*I_NAI_H5x_D2y_ab;
  Double I_NAI_H4xy_Fx2y_ab = I_NAI_I5xy_D2y_ab+ABX*I_NAI_H4xy_D2y_ab;
  Double I_NAI_H4xz_Fx2y_ab = I_NAI_I5xz_D2y_ab+ABX*I_NAI_H4xz_D2y_ab;
  Double I_NAI_H3x2y_Fx2y_ab = I_NAI_I4x2y_D2y_ab+ABX*I_NAI_H3x2y_D2y_ab;
  Double I_NAI_H3xyz_Fx2y_ab = I_NAI_I4xyz_D2y_ab+ABX*I_NAI_H3xyz_D2y_ab;
  Double I_NAI_H3x2z_Fx2y_ab = I_NAI_I4x2z_D2y_ab+ABX*I_NAI_H3x2z_D2y_ab;
  Double I_NAI_H2x3y_Fx2y_ab = I_NAI_I3x3y_D2y_ab+ABX*I_NAI_H2x3y_D2y_ab;
  Double I_NAI_H2x2yz_Fx2y_ab = I_NAI_I3x2yz_D2y_ab+ABX*I_NAI_H2x2yz_D2y_ab;
  Double I_NAI_H2xy2z_Fx2y_ab = I_NAI_I3xy2z_D2y_ab+ABX*I_NAI_H2xy2z_D2y_ab;
  Double I_NAI_H2x3z_Fx2y_ab = I_NAI_I3x3z_D2y_ab+ABX*I_NAI_H2x3z_D2y_ab;
  Double I_NAI_Hx4y_Fx2y_ab = I_NAI_I2x4y_D2y_ab+ABX*I_NAI_Hx4y_D2y_ab;
  Double I_NAI_Hx3yz_Fx2y_ab = I_NAI_I2x3yz_D2y_ab+ABX*I_NAI_Hx3yz_D2y_ab;
  Double I_NAI_Hx2y2z_Fx2y_ab = I_NAI_I2x2y2z_D2y_ab+ABX*I_NAI_Hx2y2z_D2y_ab;
  Double I_NAI_Hxy3z_Fx2y_ab = I_NAI_I2xy3z_D2y_ab+ABX*I_NAI_Hxy3z_D2y_ab;
  Double I_NAI_Hx4z_Fx2y_ab = I_NAI_I2x4z_D2y_ab+ABX*I_NAI_Hx4z_D2y_ab;
  Double I_NAI_H5y_Fx2y_ab = I_NAI_Ix5y_D2y_ab+ABX*I_NAI_H5y_D2y_ab;
  Double I_NAI_H4yz_Fx2y_ab = I_NAI_Ix4yz_D2y_ab+ABX*I_NAI_H4yz_D2y_ab;
  Double I_NAI_H3y2z_Fx2y_ab = I_NAI_Ix3y2z_D2y_ab+ABX*I_NAI_H3y2z_D2y_ab;
  Double I_NAI_H2y3z_Fx2y_ab = I_NAI_Ix2y3z_D2y_ab+ABX*I_NAI_H2y3z_D2y_ab;
  Double I_NAI_Hy4z_Fx2y_ab = I_NAI_Ixy4z_D2y_ab+ABX*I_NAI_Hy4z_D2y_ab;
  Double I_NAI_H5z_Fx2y_ab = I_NAI_Ix5z_D2y_ab+ABX*I_NAI_H5z_D2y_ab;
  Double I_NAI_H5x_Fx2z_ab = I_NAI_I6x_D2z_ab+ABX*I_NAI_H5x_D2z_ab;
  Double I_NAI_H4xy_Fx2z_ab = I_NAI_I5xy_D2z_ab+ABX*I_NAI_H4xy_D2z_ab;
  Double I_NAI_H4xz_Fx2z_ab = I_NAI_I5xz_D2z_ab+ABX*I_NAI_H4xz_D2z_ab;
  Double I_NAI_H3x2y_Fx2z_ab = I_NAI_I4x2y_D2z_ab+ABX*I_NAI_H3x2y_D2z_ab;
  Double I_NAI_H3xyz_Fx2z_ab = I_NAI_I4xyz_D2z_ab+ABX*I_NAI_H3xyz_D2z_ab;
  Double I_NAI_H3x2z_Fx2z_ab = I_NAI_I4x2z_D2z_ab+ABX*I_NAI_H3x2z_D2z_ab;
  Double I_NAI_H2x3y_Fx2z_ab = I_NAI_I3x3y_D2z_ab+ABX*I_NAI_H2x3y_D2z_ab;
  Double I_NAI_H2x2yz_Fx2z_ab = I_NAI_I3x2yz_D2z_ab+ABX*I_NAI_H2x2yz_D2z_ab;
  Double I_NAI_H2xy2z_Fx2z_ab = I_NAI_I3xy2z_D2z_ab+ABX*I_NAI_H2xy2z_D2z_ab;
  Double I_NAI_H2x3z_Fx2z_ab = I_NAI_I3x3z_D2z_ab+ABX*I_NAI_H2x3z_D2z_ab;
  Double I_NAI_Hx4y_Fx2z_ab = I_NAI_I2x4y_D2z_ab+ABX*I_NAI_Hx4y_D2z_ab;
  Double I_NAI_Hx3yz_Fx2z_ab = I_NAI_I2x3yz_D2z_ab+ABX*I_NAI_Hx3yz_D2z_ab;
  Double I_NAI_Hx2y2z_Fx2z_ab = I_NAI_I2x2y2z_D2z_ab+ABX*I_NAI_Hx2y2z_D2z_ab;
  Double I_NAI_Hxy3z_Fx2z_ab = I_NAI_I2xy3z_D2z_ab+ABX*I_NAI_Hxy3z_D2z_ab;
  Double I_NAI_Hx4z_Fx2z_ab = I_NAI_I2x4z_D2z_ab+ABX*I_NAI_Hx4z_D2z_ab;
  Double I_NAI_H5y_Fx2z_ab = I_NAI_Ix5y_D2z_ab+ABX*I_NAI_H5y_D2z_ab;
  Double I_NAI_H4yz_Fx2z_ab = I_NAI_Ix4yz_D2z_ab+ABX*I_NAI_H4yz_D2z_ab;
  Double I_NAI_H3y2z_Fx2z_ab = I_NAI_Ix3y2z_D2z_ab+ABX*I_NAI_H3y2z_D2z_ab;
  Double I_NAI_H2y3z_Fx2z_ab = I_NAI_Ix2y3z_D2z_ab+ABX*I_NAI_H2y3z_D2z_ab;
  Double I_NAI_Hy4z_Fx2z_ab = I_NAI_Ixy4z_D2z_ab+ABX*I_NAI_Hy4z_D2z_ab;
  Double I_NAI_H5z_Fx2z_ab = I_NAI_Ix5z_D2z_ab+ABX*I_NAI_H5z_D2z_ab;
  Double I_NAI_H5x_F3y_ab = I_NAI_I5xy_D2y_ab+ABY*I_NAI_H5x_D2y_ab;
  Double I_NAI_H4xy_F3y_ab = I_NAI_I4x2y_D2y_ab+ABY*I_NAI_H4xy_D2y_ab;
  Double I_NAI_H4xz_F3y_ab = I_NAI_I4xyz_D2y_ab+ABY*I_NAI_H4xz_D2y_ab;
  Double I_NAI_H3x2y_F3y_ab = I_NAI_I3x3y_D2y_ab+ABY*I_NAI_H3x2y_D2y_ab;
  Double I_NAI_H3xyz_F3y_ab = I_NAI_I3x2yz_D2y_ab+ABY*I_NAI_H3xyz_D2y_ab;
  Double I_NAI_H3x2z_F3y_ab = I_NAI_I3xy2z_D2y_ab+ABY*I_NAI_H3x2z_D2y_ab;
  Double I_NAI_H2x3y_F3y_ab = I_NAI_I2x4y_D2y_ab+ABY*I_NAI_H2x3y_D2y_ab;
  Double I_NAI_H2x2yz_F3y_ab = I_NAI_I2x3yz_D2y_ab+ABY*I_NAI_H2x2yz_D2y_ab;
  Double I_NAI_H2xy2z_F3y_ab = I_NAI_I2x2y2z_D2y_ab+ABY*I_NAI_H2xy2z_D2y_ab;
  Double I_NAI_H2x3z_F3y_ab = I_NAI_I2xy3z_D2y_ab+ABY*I_NAI_H2x3z_D2y_ab;
  Double I_NAI_Hx4y_F3y_ab = I_NAI_Ix5y_D2y_ab+ABY*I_NAI_Hx4y_D2y_ab;
  Double I_NAI_Hx3yz_F3y_ab = I_NAI_Ix4yz_D2y_ab+ABY*I_NAI_Hx3yz_D2y_ab;
  Double I_NAI_Hx2y2z_F3y_ab = I_NAI_Ix3y2z_D2y_ab+ABY*I_NAI_Hx2y2z_D2y_ab;
  Double I_NAI_Hxy3z_F3y_ab = I_NAI_Ix2y3z_D2y_ab+ABY*I_NAI_Hxy3z_D2y_ab;
  Double I_NAI_Hx4z_F3y_ab = I_NAI_Ixy4z_D2y_ab+ABY*I_NAI_Hx4z_D2y_ab;
  Double I_NAI_H5y_F3y_ab = I_NAI_I6y_D2y_ab+ABY*I_NAI_H5y_D2y_ab;
  Double I_NAI_H4yz_F3y_ab = I_NAI_I5yz_D2y_ab+ABY*I_NAI_H4yz_D2y_ab;
  Double I_NAI_H3y2z_F3y_ab = I_NAI_I4y2z_D2y_ab+ABY*I_NAI_H3y2z_D2y_ab;
  Double I_NAI_H2y3z_F3y_ab = I_NAI_I3y3z_D2y_ab+ABY*I_NAI_H2y3z_D2y_ab;
  Double I_NAI_Hy4z_F3y_ab = I_NAI_I2y4z_D2y_ab+ABY*I_NAI_Hy4z_D2y_ab;
  Double I_NAI_H5z_F3y_ab = I_NAI_Iy5z_D2y_ab+ABY*I_NAI_H5z_D2y_ab;
  Double I_NAI_H5x_F2yz_ab = I_NAI_I5xz_D2y_ab+ABZ*I_NAI_H5x_D2y_ab;
  Double I_NAI_H4xy_F2yz_ab = I_NAI_I4xyz_D2y_ab+ABZ*I_NAI_H4xy_D2y_ab;
  Double I_NAI_H4xz_F2yz_ab = I_NAI_I4x2z_D2y_ab+ABZ*I_NAI_H4xz_D2y_ab;
  Double I_NAI_H3x2y_F2yz_ab = I_NAI_I3x2yz_D2y_ab+ABZ*I_NAI_H3x2y_D2y_ab;
  Double I_NAI_H3xyz_F2yz_ab = I_NAI_I3xy2z_D2y_ab+ABZ*I_NAI_H3xyz_D2y_ab;
  Double I_NAI_H3x2z_F2yz_ab = I_NAI_I3x3z_D2y_ab+ABZ*I_NAI_H3x2z_D2y_ab;
  Double I_NAI_H2x3y_F2yz_ab = I_NAI_I2x3yz_D2y_ab+ABZ*I_NAI_H2x3y_D2y_ab;
  Double I_NAI_H2x2yz_F2yz_ab = I_NAI_I2x2y2z_D2y_ab+ABZ*I_NAI_H2x2yz_D2y_ab;
  Double I_NAI_H2xy2z_F2yz_ab = I_NAI_I2xy3z_D2y_ab+ABZ*I_NAI_H2xy2z_D2y_ab;
  Double I_NAI_H2x3z_F2yz_ab = I_NAI_I2x4z_D2y_ab+ABZ*I_NAI_H2x3z_D2y_ab;
  Double I_NAI_Hx4y_F2yz_ab = I_NAI_Ix4yz_D2y_ab+ABZ*I_NAI_Hx4y_D2y_ab;
  Double I_NAI_Hx3yz_F2yz_ab = I_NAI_Ix3y2z_D2y_ab+ABZ*I_NAI_Hx3yz_D2y_ab;
  Double I_NAI_Hx2y2z_F2yz_ab = I_NAI_Ix2y3z_D2y_ab+ABZ*I_NAI_Hx2y2z_D2y_ab;
  Double I_NAI_Hxy3z_F2yz_ab = I_NAI_Ixy4z_D2y_ab+ABZ*I_NAI_Hxy3z_D2y_ab;
  Double I_NAI_Hx4z_F2yz_ab = I_NAI_Ix5z_D2y_ab+ABZ*I_NAI_Hx4z_D2y_ab;
  Double I_NAI_H5y_F2yz_ab = I_NAI_I5yz_D2y_ab+ABZ*I_NAI_H5y_D2y_ab;
  Double I_NAI_H4yz_F2yz_ab = I_NAI_I4y2z_D2y_ab+ABZ*I_NAI_H4yz_D2y_ab;
  Double I_NAI_H3y2z_F2yz_ab = I_NAI_I3y3z_D2y_ab+ABZ*I_NAI_H3y2z_D2y_ab;
  Double I_NAI_H2y3z_F2yz_ab = I_NAI_I2y4z_D2y_ab+ABZ*I_NAI_H2y3z_D2y_ab;
  Double I_NAI_Hy4z_F2yz_ab = I_NAI_Iy5z_D2y_ab+ABZ*I_NAI_Hy4z_D2y_ab;
  Double I_NAI_H5z_F2yz_ab = I_NAI_I6z_D2y_ab+ABZ*I_NAI_H5z_D2y_ab;
  Double I_NAI_H5x_F3z_ab = I_NAI_I5xz_D2z_ab+ABZ*I_NAI_H5x_D2z_ab;
  Double I_NAI_H4xy_F3z_ab = I_NAI_I4xyz_D2z_ab+ABZ*I_NAI_H4xy_D2z_ab;
  Double I_NAI_H4xz_F3z_ab = I_NAI_I4x2z_D2z_ab+ABZ*I_NAI_H4xz_D2z_ab;
  Double I_NAI_H3x2y_F3z_ab = I_NAI_I3x2yz_D2z_ab+ABZ*I_NAI_H3x2y_D2z_ab;
  Double I_NAI_H3xyz_F3z_ab = I_NAI_I3xy2z_D2z_ab+ABZ*I_NAI_H3xyz_D2z_ab;
  Double I_NAI_H3x2z_F3z_ab = I_NAI_I3x3z_D2z_ab+ABZ*I_NAI_H3x2z_D2z_ab;
  Double I_NAI_H2x3y_F3z_ab = I_NAI_I2x3yz_D2z_ab+ABZ*I_NAI_H2x3y_D2z_ab;
  Double I_NAI_H2x2yz_F3z_ab = I_NAI_I2x2y2z_D2z_ab+ABZ*I_NAI_H2x2yz_D2z_ab;
  Double I_NAI_H2xy2z_F3z_ab = I_NAI_I2xy3z_D2z_ab+ABZ*I_NAI_H2xy2z_D2z_ab;
  Double I_NAI_H2x3z_F3z_ab = I_NAI_I2x4z_D2z_ab+ABZ*I_NAI_H2x3z_D2z_ab;
  Double I_NAI_Hx4y_F3z_ab = I_NAI_Ix4yz_D2z_ab+ABZ*I_NAI_Hx4y_D2z_ab;
  Double I_NAI_Hx3yz_F3z_ab = I_NAI_Ix3y2z_D2z_ab+ABZ*I_NAI_Hx3yz_D2z_ab;
  Double I_NAI_Hx2y2z_F3z_ab = I_NAI_Ix2y3z_D2z_ab+ABZ*I_NAI_Hx2y2z_D2z_ab;
  Double I_NAI_Hxy3z_F3z_ab = I_NAI_Ixy4z_D2z_ab+ABZ*I_NAI_Hxy3z_D2z_ab;
  Double I_NAI_Hx4z_F3z_ab = I_NAI_Ix5z_D2z_ab+ABZ*I_NAI_Hx4z_D2z_ab;
  Double I_NAI_H5y_F3z_ab = I_NAI_I5yz_D2z_ab+ABZ*I_NAI_H5y_D2z_ab;
  Double I_NAI_H4yz_F3z_ab = I_NAI_I4y2z_D2z_ab+ABZ*I_NAI_H4yz_D2z_ab;
  Double I_NAI_H3y2z_F3z_ab = I_NAI_I3y3z_D2z_ab+ABZ*I_NAI_H3y2z_D2z_ab;
  Double I_NAI_H2y3z_F3z_ab = I_NAI_I2y4z_D2z_ab+ABZ*I_NAI_H2y3z_D2z_ab;
  Double I_NAI_Hy4z_F3z_ab = I_NAI_Iy5z_D2z_ab+ABZ*I_NAI_Hy4z_D2z_ab;
  Double I_NAI_H5z_F3z_ab = I_NAI_I6z_D2z_ab+ABZ*I_NAI_H5z_D2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_L_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_S_ab
   * RHS shell quartet name: SQ_NAI_L_S_ab
   ************************************************************/
  Double I_NAI_L8x_Px_ab = I_NAI_M9x_S_ab+ABX*I_NAI_L8x_S_ab;
  Double I_NAI_L7xy_Px_ab = I_NAI_M8xy_S_ab+ABX*I_NAI_L7xy_S_ab;
  Double I_NAI_L7xz_Px_ab = I_NAI_M8xz_S_ab+ABX*I_NAI_L7xz_S_ab;
  Double I_NAI_L6x2y_Px_ab = I_NAI_M7x2y_S_ab+ABX*I_NAI_L6x2y_S_ab;
  Double I_NAI_L6xyz_Px_ab = I_NAI_M7xyz_S_ab+ABX*I_NAI_L6xyz_S_ab;
  Double I_NAI_L6x2z_Px_ab = I_NAI_M7x2z_S_ab+ABX*I_NAI_L6x2z_S_ab;
  Double I_NAI_L5x3y_Px_ab = I_NAI_M6x3y_S_ab+ABX*I_NAI_L5x3y_S_ab;
  Double I_NAI_L5x2yz_Px_ab = I_NAI_M6x2yz_S_ab+ABX*I_NAI_L5x2yz_S_ab;
  Double I_NAI_L5xy2z_Px_ab = I_NAI_M6xy2z_S_ab+ABX*I_NAI_L5xy2z_S_ab;
  Double I_NAI_L5x3z_Px_ab = I_NAI_M6x3z_S_ab+ABX*I_NAI_L5x3z_S_ab;
  Double I_NAI_L4x4y_Px_ab = I_NAI_M5x4y_S_ab+ABX*I_NAI_L4x4y_S_ab;
  Double I_NAI_L4x3yz_Px_ab = I_NAI_M5x3yz_S_ab+ABX*I_NAI_L4x3yz_S_ab;
  Double I_NAI_L4x2y2z_Px_ab = I_NAI_M5x2y2z_S_ab+ABX*I_NAI_L4x2y2z_S_ab;
  Double I_NAI_L4xy3z_Px_ab = I_NAI_M5xy3z_S_ab+ABX*I_NAI_L4xy3z_S_ab;
  Double I_NAI_L4x4z_Px_ab = I_NAI_M5x4z_S_ab+ABX*I_NAI_L4x4z_S_ab;
  Double I_NAI_L3x5y_Px_ab = I_NAI_M4x5y_S_ab+ABX*I_NAI_L3x5y_S_ab;
  Double I_NAI_L3x4yz_Px_ab = I_NAI_M4x4yz_S_ab+ABX*I_NAI_L3x4yz_S_ab;
  Double I_NAI_L3x3y2z_Px_ab = I_NAI_M4x3y2z_S_ab+ABX*I_NAI_L3x3y2z_S_ab;
  Double I_NAI_L3x2y3z_Px_ab = I_NAI_M4x2y3z_S_ab+ABX*I_NAI_L3x2y3z_S_ab;
  Double I_NAI_L3xy4z_Px_ab = I_NAI_M4xy4z_S_ab+ABX*I_NAI_L3xy4z_S_ab;
  Double I_NAI_L3x5z_Px_ab = I_NAI_M4x5z_S_ab+ABX*I_NAI_L3x5z_S_ab;
  Double I_NAI_L2x6y_Px_ab = I_NAI_M3x6y_S_ab+ABX*I_NAI_L2x6y_S_ab;
  Double I_NAI_L2x5yz_Px_ab = I_NAI_M3x5yz_S_ab+ABX*I_NAI_L2x5yz_S_ab;
  Double I_NAI_L2x4y2z_Px_ab = I_NAI_M3x4y2z_S_ab+ABX*I_NAI_L2x4y2z_S_ab;
  Double I_NAI_L2x3y3z_Px_ab = I_NAI_M3x3y3z_S_ab+ABX*I_NAI_L2x3y3z_S_ab;
  Double I_NAI_L2x2y4z_Px_ab = I_NAI_M3x2y4z_S_ab+ABX*I_NAI_L2x2y4z_S_ab;
  Double I_NAI_L2xy5z_Px_ab = I_NAI_M3xy5z_S_ab+ABX*I_NAI_L2xy5z_S_ab;
  Double I_NAI_L2x6z_Px_ab = I_NAI_M3x6z_S_ab+ABX*I_NAI_L2x6z_S_ab;
  Double I_NAI_Lx7y_Px_ab = I_NAI_M2x7y_S_ab+ABX*I_NAI_Lx7y_S_ab;
  Double I_NAI_Lx6yz_Px_ab = I_NAI_M2x6yz_S_ab+ABX*I_NAI_Lx6yz_S_ab;
  Double I_NAI_Lx5y2z_Px_ab = I_NAI_M2x5y2z_S_ab+ABX*I_NAI_Lx5y2z_S_ab;
  Double I_NAI_Lx4y3z_Px_ab = I_NAI_M2x4y3z_S_ab+ABX*I_NAI_Lx4y3z_S_ab;
  Double I_NAI_Lx3y4z_Px_ab = I_NAI_M2x3y4z_S_ab+ABX*I_NAI_Lx3y4z_S_ab;
  Double I_NAI_Lx2y5z_Px_ab = I_NAI_M2x2y5z_S_ab+ABX*I_NAI_Lx2y5z_S_ab;
  Double I_NAI_Lxy6z_Px_ab = I_NAI_M2xy6z_S_ab+ABX*I_NAI_Lxy6z_S_ab;
  Double I_NAI_Lx7z_Px_ab = I_NAI_M2x7z_S_ab+ABX*I_NAI_Lx7z_S_ab;
  Double I_NAI_L6x2y_Py_ab = I_NAI_M6x3y_S_ab+ABY*I_NAI_L6x2y_S_ab;
  Double I_NAI_L6xyz_Py_ab = I_NAI_M6x2yz_S_ab+ABY*I_NAI_L6xyz_S_ab;
  Double I_NAI_L5x3y_Py_ab = I_NAI_M5x4y_S_ab+ABY*I_NAI_L5x3y_S_ab;
  Double I_NAI_L5x2yz_Py_ab = I_NAI_M5x3yz_S_ab+ABY*I_NAI_L5x2yz_S_ab;
  Double I_NAI_L5xy2z_Py_ab = I_NAI_M5x2y2z_S_ab+ABY*I_NAI_L5xy2z_S_ab;
  Double I_NAI_L4x4y_Py_ab = I_NAI_M4x5y_S_ab+ABY*I_NAI_L4x4y_S_ab;
  Double I_NAI_L4x3yz_Py_ab = I_NAI_M4x4yz_S_ab+ABY*I_NAI_L4x3yz_S_ab;
  Double I_NAI_L4x2y2z_Py_ab = I_NAI_M4x3y2z_S_ab+ABY*I_NAI_L4x2y2z_S_ab;
  Double I_NAI_L4xy3z_Py_ab = I_NAI_M4x2y3z_S_ab+ABY*I_NAI_L4xy3z_S_ab;
  Double I_NAI_L3x5y_Py_ab = I_NAI_M3x6y_S_ab+ABY*I_NAI_L3x5y_S_ab;
  Double I_NAI_L3x4yz_Py_ab = I_NAI_M3x5yz_S_ab+ABY*I_NAI_L3x4yz_S_ab;
  Double I_NAI_L3x3y2z_Py_ab = I_NAI_M3x4y2z_S_ab+ABY*I_NAI_L3x3y2z_S_ab;
  Double I_NAI_L3x2y3z_Py_ab = I_NAI_M3x3y3z_S_ab+ABY*I_NAI_L3x2y3z_S_ab;
  Double I_NAI_L3xy4z_Py_ab = I_NAI_M3x2y4z_S_ab+ABY*I_NAI_L3xy4z_S_ab;
  Double I_NAI_L2x6y_Py_ab = I_NAI_M2x7y_S_ab+ABY*I_NAI_L2x6y_S_ab;
  Double I_NAI_L2x5yz_Py_ab = I_NAI_M2x6yz_S_ab+ABY*I_NAI_L2x5yz_S_ab;
  Double I_NAI_L2x4y2z_Py_ab = I_NAI_M2x5y2z_S_ab+ABY*I_NAI_L2x4y2z_S_ab;
  Double I_NAI_L2x3y3z_Py_ab = I_NAI_M2x4y3z_S_ab+ABY*I_NAI_L2x3y3z_S_ab;
  Double I_NAI_L2x2y4z_Py_ab = I_NAI_M2x3y4z_S_ab+ABY*I_NAI_L2x2y4z_S_ab;
  Double I_NAI_L2xy5z_Py_ab = I_NAI_M2x2y5z_S_ab+ABY*I_NAI_L2xy5z_S_ab;
  Double I_NAI_Lx7y_Py_ab = I_NAI_Mx8y_S_ab+ABY*I_NAI_Lx7y_S_ab;
  Double I_NAI_Lx6yz_Py_ab = I_NAI_Mx7yz_S_ab+ABY*I_NAI_Lx6yz_S_ab;
  Double I_NAI_Lx5y2z_Py_ab = I_NAI_Mx6y2z_S_ab+ABY*I_NAI_Lx5y2z_S_ab;
  Double I_NAI_Lx4y3z_Py_ab = I_NAI_Mx5y3z_S_ab+ABY*I_NAI_Lx4y3z_S_ab;
  Double I_NAI_Lx3y4z_Py_ab = I_NAI_Mx4y4z_S_ab+ABY*I_NAI_Lx3y4z_S_ab;
  Double I_NAI_Lx2y5z_Py_ab = I_NAI_Mx3y5z_S_ab+ABY*I_NAI_Lx2y5z_S_ab;
  Double I_NAI_Lxy6z_Py_ab = I_NAI_Mx2y6z_S_ab+ABY*I_NAI_Lxy6z_S_ab;
  Double I_NAI_L8y_Py_ab = I_NAI_M9y_S_ab+ABY*I_NAI_L8y_S_ab;
  Double I_NAI_L7yz_Py_ab = I_NAI_M8yz_S_ab+ABY*I_NAI_L7yz_S_ab;
  Double I_NAI_L6y2z_Py_ab = I_NAI_M7y2z_S_ab+ABY*I_NAI_L6y2z_S_ab;
  Double I_NAI_L5y3z_Py_ab = I_NAI_M6y3z_S_ab+ABY*I_NAI_L5y3z_S_ab;
  Double I_NAI_L4y4z_Py_ab = I_NAI_M5y4z_S_ab+ABY*I_NAI_L4y4z_S_ab;
  Double I_NAI_L3y5z_Py_ab = I_NAI_M4y5z_S_ab+ABY*I_NAI_L3y5z_S_ab;
  Double I_NAI_L2y6z_Py_ab = I_NAI_M3y6z_S_ab+ABY*I_NAI_L2y6z_S_ab;
  Double I_NAI_Ly7z_Py_ab = I_NAI_M2y7z_S_ab+ABY*I_NAI_Ly7z_S_ab;
  Double I_NAI_L6xyz_Pz_ab = I_NAI_M6xy2z_S_ab+ABZ*I_NAI_L6xyz_S_ab;
  Double I_NAI_L6x2z_Pz_ab = I_NAI_M6x3z_S_ab+ABZ*I_NAI_L6x2z_S_ab;
  Double I_NAI_L5x2yz_Pz_ab = I_NAI_M5x2y2z_S_ab+ABZ*I_NAI_L5x2yz_S_ab;
  Double I_NAI_L5xy2z_Pz_ab = I_NAI_M5xy3z_S_ab+ABZ*I_NAI_L5xy2z_S_ab;
  Double I_NAI_L5x3z_Pz_ab = I_NAI_M5x4z_S_ab+ABZ*I_NAI_L5x3z_S_ab;
  Double I_NAI_L4x3yz_Pz_ab = I_NAI_M4x3y2z_S_ab+ABZ*I_NAI_L4x3yz_S_ab;
  Double I_NAI_L4x2y2z_Pz_ab = I_NAI_M4x2y3z_S_ab+ABZ*I_NAI_L4x2y2z_S_ab;
  Double I_NAI_L4xy3z_Pz_ab = I_NAI_M4xy4z_S_ab+ABZ*I_NAI_L4xy3z_S_ab;
  Double I_NAI_L4x4z_Pz_ab = I_NAI_M4x5z_S_ab+ABZ*I_NAI_L4x4z_S_ab;
  Double I_NAI_L3x4yz_Pz_ab = I_NAI_M3x4y2z_S_ab+ABZ*I_NAI_L3x4yz_S_ab;
  Double I_NAI_L3x3y2z_Pz_ab = I_NAI_M3x3y3z_S_ab+ABZ*I_NAI_L3x3y2z_S_ab;
  Double I_NAI_L3x2y3z_Pz_ab = I_NAI_M3x2y4z_S_ab+ABZ*I_NAI_L3x2y3z_S_ab;
  Double I_NAI_L3xy4z_Pz_ab = I_NAI_M3xy5z_S_ab+ABZ*I_NAI_L3xy4z_S_ab;
  Double I_NAI_L3x5z_Pz_ab = I_NAI_M3x6z_S_ab+ABZ*I_NAI_L3x5z_S_ab;
  Double I_NAI_L2x5yz_Pz_ab = I_NAI_M2x5y2z_S_ab+ABZ*I_NAI_L2x5yz_S_ab;
  Double I_NAI_L2x4y2z_Pz_ab = I_NAI_M2x4y3z_S_ab+ABZ*I_NAI_L2x4y2z_S_ab;
  Double I_NAI_L2x3y3z_Pz_ab = I_NAI_M2x3y4z_S_ab+ABZ*I_NAI_L2x3y3z_S_ab;
  Double I_NAI_L2x2y4z_Pz_ab = I_NAI_M2x2y5z_S_ab+ABZ*I_NAI_L2x2y4z_S_ab;
  Double I_NAI_L2xy5z_Pz_ab = I_NAI_M2xy6z_S_ab+ABZ*I_NAI_L2xy5z_S_ab;
  Double I_NAI_L2x6z_Pz_ab = I_NAI_M2x7z_S_ab+ABZ*I_NAI_L2x6z_S_ab;
  Double I_NAI_Lx6yz_Pz_ab = I_NAI_Mx6y2z_S_ab+ABZ*I_NAI_Lx6yz_S_ab;
  Double I_NAI_Lx5y2z_Pz_ab = I_NAI_Mx5y3z_S_ab+ABZ*I_NAI_Lx5y2z_S_ab;
  Double I_NAI_Lx4y3z_Pz_ab = I_NAI_Mx4y4z_S_ab+ABZ*I_NAI_Lx4y3z_S_ab;
  Double I_NAI_Lx3y4z_Pz_ab = I_NAI_Mx3y5z_S_ab+ABZ*I_NAI_Lx3y4z_S_ab;
  Double I_NAI_Lx2y5z_Pz_ab = I_NAI_Mx2y6z_S_ab+ABZ*I_NAI_Lx2y5z_S_ab;
  Double I_NAI_Lxy6z_Pz_ab = I_NAI_Mxy7z_S_ab+ABZ*I_NAI_Lxy6z_S_ab;
  Double I_NAI_Lx7z_Pz_ab = I_NAI_Mx8z_S_ab+ABZ*I_NAI_Lx7z_S_ab;
  Double I_NAI_L6y2z_Pz_ab = I_NAI_M6y3z_S_ab+ABZ*I_NAI_L6y2z_S_ab;
  Double I_NAI_L5y3z_Pz_ab = I_NAI_M5y4z_S_ab+ABZ*I_NAI_L5y3z_S_ab;
  Double I_NAI_L4y4z_Pz_ab = I_NAI_M4y5z_S_ab+ABZ*I_NAI_L4y4z_S_ab;
  Double I_NAI_L3y5z_Pz_ab = I_NAI_M3y6z_S_ab+ABZ*I_NAI_L3y5z_S_ab;
  Double I_NAI_L2y6z_Pz_ab = I_NAI_M2y7z_S_ab+ABZ*I_NAI_L2y6z_S_ab;
  Double I_NAI_Ly7z_Pz_ab = I_NAI_My8z_S_ab+ABZ*I_NAI_Ly7z_S_ab;
  Double I_NAI_L8z_Pz_ab = I_NAI_M9z_S_ab+ABZ*I_NAI_L8z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_K_D_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 111 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_P_ab
   * RHS shell quartet name: SQ_NAI_K_P_ab
   ************************************************************/
  Double I_NAI_K7x_D2x_ab = I_NAI_L8x_Px_ab+ABX*I_NAI_K7x_Px_ab;
  Double I_NAI_K6xy_D2x_ab = I_NAI_L7xy_Px_ab+ABX*I_NAI_K6xy_Px_ab;
  Double I_NAI_K6xz_D2x_ab = I_NAI_L7xz_Px_ab+ABX*I_NAI_K6xz_Px_ab;
  Double I_NAI_K5x2y_D2x_ab = I_NAI_L6x2y_Px_ab+ABX*I_NAI_K5x2y_Px_ab;
  Double I_NAI_K5xyz_D2x_ab = I_NAI_L6xyz_Px_ab+ABX*I_NAI_K5xyz_Px_ab;
  Double I_NAI_K5x2z_D2x_ab = I_NAI_L6x2z_Px_ab+ABX*I_NAI_K5x2z_Px_ab;
  Double I_NAI_K4x3y_D2x_ab = I_NAI_L5x3y_Px_ab+ABX*I_NAI_K4x3y_Px_ab;
  Double I_NAI_K4x2yz_D2x_ab = I_NAI_L5x2yz_Px_ab+ABX*I_NAI_K4x2yz_Px_ab;
  Double I_NAI_K4xy2z_D2x_ab = I_NAI_L5xy2z_Px_ab+ABX*I_NAI_K4xy2z_Px_ab;
  Double I_NAI_K4x3z_D2x_ab = I_NAI_L5x3z_Px_ab+ABX*I_NAI_K4x3z_Px_ab;
  Double I_NAI_K3x4y_D2x_ab = I_NAI_L4x4y_Px_ab+ABX*I_NAI_K3x4y_Px_ab;
  Double I_NAI_K3x3yz_D2x_ab = I_NAI_L4x3yz_Px_ab+ABX*I_NAI_K3x3yz_Px_ab;
  Double I_NAI_K3x2y2z_D2x_ab = I_NAI_L4x2y2z_Px_ab+ABX*I_NAI_K3x2y2z_Px_ab;
  Double I_NAI_K3xy3z_D2x_ab = I_NAI_L4xy3z_Px_ab+ABX*I_NAI_K3xy3z_Px_ab;
  Double I_NAI_K3x4z_D2x_ab = I_NAI_L4x4z_Px_ab+ABX*I_NAI_K3x4z_Px_ab;
  Double I_NAI_K2x5y_D2x_ab = I_NAI_L3x5y_Px_ab+ABX*I_NAI_K2x5y_Px_ab;
  Double I_NAI_K2x4yz_D2x_ab = I_NAI_L3x4yz_Px_ab+ABX*I_NAI_K2x4yz_Px_ab;
  Double I_NAI_K2x3y2z_D2x_ab = I_NAI_L3x3y2z_Px_ab+ABX*I_NAI_K2x3y2z_Px_ab;
  Double I_NAI_K2x2y3z_D2x_ab = I_NAI_L3x2y3z_Px_ab+ABX*I_NAI_K2x2y3z_Px_ab;
  Double I_NAI_K2xy4z_D2x_ab = I_NAI_L3xy4z_Px_ab+ABX*I_NAI_K2xy4z_Px_ab;
  Double I_NAI_K2x5z_D2x_ab = I_NAI_L3x5z_Px_ab+ABX*I_NAI_K2x5z_Px_ab;
  Double I_NAI_Kx6y_D2x_ab = I_NAI_L2x6y_Px_ab+ABX*I_NAI_Kx6y_Px_ab;
  Double I_NAI_Kx5yz_D2x_ab = I_NAI_L2x5yz_Px_ab+ABX*I_NAI_Kx5yz_Px_ab;
  Double I_NAI_Kx4y2z_D2x_ab = I_NAI_L2x4y2z_Px_ab+ABX*I_NAI_Kx4y2z_Px_ab;
  Double I_NAI_Kx3y3z_D2x_ab = I_NAI_L2x3y3z_Px_ab+ABX*I_NAI_Kx3y3z_Px_ab;
  Double I_NAI_Kx2y4z_D2x_ab = I_NAI_L2x2y4z_Px_ab+ABX*I_NAI_Kx2y4z_Px_ab;
  Double I_NAI_Kxy5z_D2x_ab = I_NAI_L2xy5z_Px_ab+ABX*I_NAI_Kxy5z_Px_ab;
  Double I_NAI_Kx6z_D2x_ab = I_NAI_L2x6z_Px_ab+ABX*I_NAI_Kx6z_Px_ab;
  Double I_NAI_K7y_D2x_ab = I_NAI_Lx7y_Px_ab+ABX*I_NAI_K7y_Px_ab;
  Double I_NAI_K6yz_D2x_ab = I_NAI_Lx6yz_Px_ab+ABX*I_NAI_K6yz_Px_ab;
  Double I_NAI_K5y2z_D2x_ab = I_NAI_Lx5y2z_Px_ab+ABX*I_NAI_K5y2z_Px_ab;
  Double I_NAI_K4y3z_D2x_ab = I_NAI_Lx4y3z_Px_ab+ABX*I_NAI_K4y3z_Px_ab;
  Double I_NAI_K3y4z_D2x_ab = I_NAI_Lx3y4z_Px_ab+ABX*I_NAI_K3y4z_Px_ab;
  Double I_NAI_K2y5z_D2x_ab = I_NAI_Lx2y5z_Px_ab+ABX*I_NAI_K2y5z_Px_ab;
  Double I_NAI_Ky6z_D2x_ab = I_NAI_Lxy6z_Px_ab+ABX*I_NAI_Ky6z_Px_ab;
  Double I_NAI_K7z_D2x_ab = I_NAI_Lx7z_Px_ab+ABX*I_NAI_K7z_Px_ab;
  Double I_NAI_K6xy_D2y_ab = I_NAI_L6x2y_Py_ab+ABY*I_NAI_K6xy_Py_ab;
  Double I_NAI_K6xz_D2y_ab = I_NAI_L6xyz_Py_ab+ABY*I_NAI_K6xz_Py_ab;
  Double I_NAI_K5x2y_D2y_ab = I_NAI_L5x3y_Py_ab+ABY*I_NAI_K5x2y_Py_ab;
  Double I_NAI_K5xyz_D2y_ab = I_NAI_L5x2yz_Py_ab+ABY*I_NAI_K5xyz_Py_ab;
  Double I_NAI_K5x2z_D2y_ab = I_NAI_L5xy2z_Py_ab+ABY*I_NAI_K5x2z_Py_ab;
  Double I_NAI_K4x3y_D2y_ab = I_NAI_L4x4y_Py_ab+ABY*I_NAI_K4x3y_Py_ab;
  Double I_NAI_K4x2yz_D2y_ab = I_NAI_L4x3yz_Py_ab+ABY*I_NAI_K4x2yz_Py_ab;
  Double I_NAI_K4xy2z_D2y_ab = I_NAI_L4x2y2z_Py_ab+ABY*I_NAI_K4xy2z_Py_ab;
  Double I_NAI_K4x3z_D2y_ab = I_NAI_L4xy3z_Py_ab+ABY*I_NAI_K4x3z_Py_ab;
  Double I_NAI_K3x4y_D2y_ab = I_NAI_L3x5y_Py_ab+ABY*I_NAI_K3x4y_Py_ab;
  Double I_NAI_K3x3yz_D2y_ab = I_NAI_L3x4yz_Py_ab+ABY*I_NAI_K3x3yz_Py_ab;
  Double I_NAI_K3x2y2z_D2y_ab = I_NAI_L3x3y2z_Py_ab+ABY*I_NAI_K3x2y2z_Py_ab;
  Double I_NAI_K3xy3z_D2y_ab = I_NAI_L3x2y3z_Py_ab+ABY*I_NAI_K3xy3z_Py_ab;
  Double I_NAI_K3x4z_D2y_ab = I_NAI_L3xy4z_Py_ab+ABY*I_NAI_K3x4z_Py_ab;
  Double I_NAI_K2x5y_D2y_ab = I_NAI_L2x6y_Py_ab+ABY*I_NAI_K2x5y_Py_ab;
  Double I_NAI_K2x4yz_D2y_ab = I_NAI_L2x5yz_Py_ab+ABY*I_NAI_K2x4yz_Py_ab;
  Double I_NAI_K2x3y2z_D2y_ab = I_NAI_L2x4y2z_Py_ab+ABY*I_NAI_K2x3y2z_Py_ab;
  Double I_NAI_K2x2y3z_D2y_ab = I_NAI_L2x3y3z_Py_ab+ABY*I_NAI_K2x2y3z_Py_ab;
  Double I_NAI_K2xy4z_D2y_ab = I_NAI_L2x2y4z_Py_ab+ABY*I_NAI_K2xy4z_Py_ab;
  Double I_NAI_K2x5z_D2y_ab = I_NAI_L2xy5z_Py_ab+ABY*I_NAI_K2x5z_Py_ab;
  Double I_NAI_Kx6y_D2y_ab = I_NAI_Lx7y_Py_ab+ABY*I_NAI_Kx6y_Py_ab;
  Double I_NAI_Kx5yz_D2y_ab = I_NAI_Lx6yz_Py_ab+ABY*I_NAI_Kx5yz_Py_ab;
  Double I_NAI_Kx4y2z_D2y_ab = I_NAI_Lx5y2z_Py_ab+ABY*I_NAI_Kx4y2z_Py_ab;
  Double I_NAI_Kx3y3z_D2y_ab = I_NAI_Lx4y3z_Py_ab+ABY*I_NAI_Kx3y3z_Py_ab;
  Double I_NAI_Kx2y4z_D2y_ab = I_NAI_Lx3y4z_Py_ab+ABY*I_NAI_Kx2y4z_Py_ab;
  Double I_NAI_Kxy5z_D2y_ab = I_NAI_Lx2y5z_Py_ab+ABY*I_NAI_Kxy5z_Py_ab;
  Double I_NAI_Kx6z_D2y_ab = I_NAI_Lxy6z_Py_ab+ABY*I_NAI_Kx6z_Py_ab;
  Double I_NAI_K7y_D2y_ab = I_NAI_L8y_Py_ab+ABY*I_NAI_K7y_Py_ab;
  Double I_NAI_K6yz_D2y_ab = I_NAI_L7yz_Py_ab+ABY*I_NAI_K6yz_Py_ab;
  Double I_NAI_K5y2z_D2y_ab = I_NAI_L6y2z_Py_ab+ABY*I_NAI_K5y2z_Py_ab;
  Double I_NAI_K4y3z_D2y_ab = I_NAI_L5y3z_Py_ab+ABY*I_NAI_K4y3z_Py_ab;
  Double I_NAI_K3y4z_D2y_ab = I_NAI_L4y4z_Py_ab+ABY*I_NAI_K3y4z_Py_ab;
  Double I_NAI_K2y5z_D2y_ab = I_NAI_L3y5z_Py_ab+ABY*I_NAI_K2y5z_Py_ab;
  Double I_NAI_Ky6z_D2y_ab = I_NAI_L2y6z_Py_ab+ABY*I_NAI_Ky6z_Py_ab;
  Double I_NAI_K7z_D2y_ab = I_NAI_Ly7z_Py_ab+ABY*I_NAI_K7z_Py_ab;
  Double I_NAI_K6xy_D2z_ab = I_NAI_L6xyz_Pz_ab+ABZ*I_NAI_K6xy_Pz_ab;
  Double I_NAI_K6xz_D2z_ab = I_NAI_L6x2z_Pz_ab+ABZ*I_NAI_K6xz_Pz_ab;
  Double I_NAI_K5x2y_D2z_ab = I_NAI_L5x2yz_Pz_ab+ABZ*I_NAI_K5x2y_Pz_ab;
  Double I_NAI_K5xyz_D2z_ab = I_NAI_L5xy2z_Pz_ab+ABZ*I_NAI_K5xyz_Pz_ab;
  Double I_NAI_K5x2z_D2z_ab = I_NAI_L5x3z_Pz_ab+ABZ*I_NAI_K5x2z_Pz_ab;
  Double I_NAI_K4x3y_D2z_ab = I_NAI_L4x3yz_Pz_ab+ABZ*I_NAI_K4x3y_Pz_ab;
  Double I_NAI_K4x2yz_D2z_ab = I_NAI_L4x2y2z_Pz_ab+ABZ*I_NAI_K4x2yz_Pz_ab;
  Double I_NAI_K4xy2z_D2z_ab = I_NAI_L4xy3z_Pz_ab+ABZ*I_NAI_K4xy2z_Pz_ab;
  Double I_NAI_K4x3z_D2z_ab = I_NAI_L4x4z_Pz_ab+ABZ*I_NAI_K4x3z_Pz_ab;
  Double I_NAI_K3x4y_D2z_ab = I_NAI_L3x4yz_Pz_ab+ABZ*I_NAI_K3x4y_Pz_ab;
  Double I_NAI_K3x3yz_D2z_ab = I_NAI_L3x3y2z_Pz_ab+ABZ*I_NAI_K3x3yz_Pz_ab;
  Double I_NAI_K3x2y2z_D2z_ab = I_NAI_L3x2y3z_Pz_ab+ABZ*I_NAI_K3x2y2z_Pz_ab;
  Double I_NAI_K3xy3z_D2z_ab = I_NAI_L3xy4z_Pz_ab+ABZ*I_NAI_K3xy3z_Pz_ab;
  Double I_NAI_K3x4z_D2z_ab = I_NAI_L3x5z_Pz_ab+ABZ*I_NAI_K3x4z_Pz_ab;
  Double I_NAI_K2x5y_D2z_ab = I_NAI_L2x5yz_Pz_ab+ABZ*I_NAI_K2x5y_Pz_ab;
  Double I_NAI_K2x4yz_D2z_ab = I_NAI_L2x4y2z_Pz_ab+ABZ*I_NAI_K2x4yz_Pz_ab;
  Double I_NAI_K2x3y2z_D2z_ab = I_NAI_L2x3y3z_Pz_ab+ABZ*I_NAI_K2x3y2z_Pz_ab;
  Double I_NAI_K2x2y3z_D2z_ab = I_NAI_L2x2y4z_Pz_ab+ABZ*I_NAI_K2x2y3z_Pz_ab;
  Double I_NAI_K2xy4z_D2z_ab = I_NAI_L2xy5z_Pz_ab+ABZ*I_NAI_K2xy4z_Pz_ab;
  Double I_NAI_K2x5z_D2z_ab = I_NAI_L2x6z_Pz_ab+ABZ*I_NAI_K2x5z_Pz_ab;
  Double I_NAI_Kx6y_D2z_ab = I_NAI_Lx6yz_Pz_ab+ABZ*I_NAI_Kx6y_Pz_ab;
  Double I_NAI_Kx5yz_D2z_ab = I_NAI_Lx5y2z_Pz_ab+ABZ*I_NAI_Kx5yz_Pz_ab;
  Double I_NAI_Kx4y2z_D2z_ab = I_NAI_Lx4y3z_Pz_ab+ABZ*I_NAI_Kx4y2z_Pz_ab;
  Double I_NAI_Kx3y3z_D2z_ab = I_NAI_Lx3y4z_Pz_ab+ABZ*I_NAI_Kx3y3z_Pz_ab;
  Double I_NAI_Kx2y4z_D2z_ab = I_NAI_Lx2y5z_Pz_ab+ABZ*I_NAI_Kx2y4z_Pz_ab;
  Double I_NAI_Kxy5z_D2z_ab = I_NAI_Lxy6z_Pz_ab+ABZ*I_NAI_Kxy5z_Pz_ab;
  Double I_NAI_Kx6z_D2z_ab = I_NAI_Lx7z_Pz_ab+ABZ*I_NAI_Kx6z_Pz_ab;
  Double I_NAI_K6yz_D2z_ab = I_NAI_L6y2z_Pz_ab+ABZ*I_NAI_K6yz_Pz_ab;
  Double I_NAI_K5y2z_D2z_ab = I_NAI_L5y3z_Pz_ab+ABZ*I_NAI_K5y2z_Pz_ab;
  Double I_NAI_K4y3z_D2z_ab = I_NAI_L4y4z_Pz_ab+ABZ*I_NAI_K4y3z_Pz_ab;
  Double I_NAI_K3y4z_D2z_ab = I_NAI_L3y5z_Pz_ab+ABZ*I_NAI_K3y4z_Pz_ab;
  Double I_NAI_K2y5z_D2z_ab = I_NAI_L2y6z_Pz_ab+ABZ*I_NAI_K2y5z_Pz_ab;
  Double I_NAI_Ky6z_D2z_ab = I_NAI_Ly7z_Pz_ab+ABZ*I_NAI_Ky6z_Pz_ab;
  Double I_NAI_K7z_D2z_ab = I_NAI_L8z_Pz_ab+ABZ*I_NAI_K7z_Pz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_I_F_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 85 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_D_ab
   * RHS shell quartet name: SQ_NAI_I_D_ab
   ************************************************************/
  Double I_NAI_I6x_F3x_ab = I_NAI_K7x_D2x_ab+ABX*I_NAI_I6x_D2x_ab;
  Double I_NAI_I5xy_F3x_ab = I_NAI_K6xy_D2x_ab+ABX*I_NAI_I5xy_D2x_ab;
  Double I_NAI_I5xz_F3x_ab = I_NAI_K6xz_D2x_ab+ABX*I_NAI_I5xz_D2x_ab;
  Double I_NAI_I4x2y_F3x_ab = I_NAI_K5x2y_D2x_ab+ABX*I_NAI_I4x2y_D2x_ab;
  Double I_NAI_I4xyz_F3x_ab = I_NAI_K5xyz_D2x_ab+ABX*I_NAI_I4xyz_D2x_ab;
  Double I_NAI_I4x2z_F3x_ab = I_NAI_K5x2z_D2x_ab+ABX*I_NAI_I4x2z_D2x_ab;
  Double I_NAI_I3x3y_F3x_ab = I_NAI_K4x3y_D2x_ab+ABX*I_NAI_I3x3y_D2x_ab;
  Double I_NAI_I3x2yz_F3x_ab = I_NAI_K4x2yz_D2x_ab+ABX*I_NAI_I3x2yz_D2x_ab;
  Double I_NAI_I3xy2z_F3x_ab = I_NAI_K4xy2z_D2x_ab+ABX*I_NAI_I3xy2z_D2x_ab;
  Double I_NAI_I3x3z_F3x_ab = I_NAI_K4x3z_D2x_ab+ABX*I_NAI_I3x3z_D2x_ab;
  Double I_NAI_I2x4y_F3x_ab = I_NAI_K3x4y_D2x_ab+ABX*I_NAI_I2x4y_D2x_ab;
  Double I_NAI_I2x3yz_F3x_ab = I_NAI_K3x3yz_D2x_ab+ABX*I_NAI_I2x3yz_D2x_ab;
  Double I_NAI_I2x2y2z_F3x_ab = I_NAI_K3x2y2z_D2x_ab+ABX*I_NAI_I2x2y2z_D2x_ab;
  Double I_NAI_I2xy3z_F3x_ab = I_NAI_K3xy3z_D2x_ab+ABX*I_NAI_I2xy3z_D2x_ab;
  Double I_NAI_I2x4z_F3x_ab = I_NAI_K3x4z_D2x_ab+ABX*I_NAI_I2x4z_D2x_ab;
  Double I_NAI_Ix5y_F3x_ab = I_NAI_K2x5y_D2x_ab+ABX*I_NAI_Ix5y_D2x_ab;
  Double I_NAI_Ix4yz_F3x_ab = I_NAI_K2x4yz_D2x_ab+ABX*I_NAI_Ix4yz_D2x_ab;
  Double I_NAI_Ix3y2z_F3x_ab = I_NAI_K2x3y2z_D2x_ab+ABX*I_NAI_Ix3y2z_D2x_ab;
  Double I_NAI_Ix2y3z_F3x_ab = I_NAI_K2x2y3z_D2x_ab+ABX*I_NAI_Ix2y3z_D2x_ab;
  Double I_NAI_Ixy4z_F3x_ab = I_NAI_K2xy4z_D2x_ab+ABX*I_NAI_Ixy4z_D2x_ab;
  Double I_NAI_Ix5z_F3x_ab = I_NAI_K2x5z_D2x_ab+ABX*I_NAI_Ix5z_D2x_ab;
  Double I_NAI_I6y_F3x_ab = I_NAI_Kx6y_D2x_ab+ABX*I_NAI_I6y_D2x_ab;
  Double I_NAI_I5yz_F3x_ab = I_NAI_Kx5yz_D2x_ab+ABX*I_NAI_I5yz_D2x_ab;
  Double I_NAI_I4y2z_F3x_ab = I_NAI_Kx4y2z_D2x_ab+ABX*I_NAI_I4y2z_D2x_ab;
  Double I_NAI_I3y3z_F3x_ab = I_NAI_Kx3y3z_D2x_ab+ABX*I_NAI_I3y3z_D2x_ab;
  Double I_NAI_I2y4z_F3x_ab = I_NAI_Kx2y4z_D2x_ab+ABX*I_NAI_I2y4z_D2x_ab;
  Double I_NAI_Iy5z_F3x_ab = I_NAI_Kxy5z_D2x_ab+ABX*I_NAI_Iy5z_D2x_ab;
  Double I_NAI_I6z_F3x_ab = I_NAI_Kx6z_D2x_ab+ABX*I_NAI_I6z_D2x_ab;
  Double I_NAI_I5xy_F2xy_ab = I_NAI_K5x2y_D2x_ab+ABY*I_NAI_I5xy_D2x_ab;
  Double I_NAI_I5xz_F2xy_ab = I_NAI_K5xyz_D2x_ab+ABY*I_NAI_I5xz_D2x_ab;
  Double I_NAI_I4x2y_F2xy_ab = I_NAI_K4x3y_D2x_ab+ABY*I_NAI_I4x2y_D2x_ab;
  Double I_NAI_I4xyz_F2xy_ab = I_NAI_K4x2yz_D2x_ab+ABY*I_NAI_I4xyz_D2x_ab;
  Double I_NAI_I4x2z_F2xy_ab = I_NAI_K4xy2z_D2x_ab+ABY*I_NAI_I4x2z_D2x_ab;
  Double I_NAI_I3x3y_F2xy_ab = I_NAI_K3x4y_D2x_ab+ABY*I_NAI_I3x3y_D2x_ab;
  Double I_NAI_I3x2yz_F2xy_ab = I_NAI_K3x3yz_D2x_ab+ABY*I_NAI_I3x2yz_D2x_ab;
  Double I_NAI_I3xy2z_F2xy_ab = I_NAI_K3x2y2z_D2x_ab+ABY*I_NAI_I3xy2z_D2x_ab;
  Double I_NAI_I3x3z_F2xy_ab = I_NAI_K3xy3z_D2x_ab+ABY*I_NAI_I3x3z_D2x_ab;
  Double I_NAI_I2x4y_F2xy_ab = I_NAI_K2x5y_D2x_ab+ABY*I_NAI_I2x4y_D2x_ab;
  Double I_NAI_I2x3yz_F2xy_ab = I_NAI_K2x4yz_D2x_ab+ABY*I_NAI_I2x3yz_D2x_ab;
  Double I_NAI_I2x2y2z_F2xy_ab = I_NAI_K2x3y2z_D2x_ab+ABY*I_NAI_I2x2y2z_D2x_ab;
  Double I_NAI_I2xy3z_F2xy_ab = I_NAI_K2x2y3z_D2x_ab+ABY*I_NAI_I2xy3z_D2x_ab;
  Double I_NAI_I2x4z_F2xy_ab = I_NAI_K2xy4z_D2x_ab+ABY*I_NAI_I2x4z_D2x_ab;
  Double I_NAI_Ix5y_F2xy_ab = I_NAI_Kx6y_D2x_ab+ABY*I_NAI_Ix5y_D2x_ab;
  Double I_NAI_Ix4yz_F2xy_ab = I_NAI_Kx5yz_D2x_ab+ABY*I_NAI_Ix4yz_D2x_ab;
  Double I_NAI_Ix3y2z_F2xy_ab = I_NAI_Kx4y2z_D2x_ab+ABY*I_NAI_Ix3y2z_D2x_ab;
  Double I_NAI_Ix2y3z_F2xy_ab = I_NAI_Kx3y3z_D2x_ab+ABY*I_NAI_Ix2y3z_D2x_ab;
  Double I_NAI_Ixy4z_F2xy_ab = I_NAI_Kx2y4z_D2x_ab+ABY*I_NAI_Ixy4z_D2x_ab;
  Double I_NAI_Ix5z_F2xy_ab = I_NAI_Kxy5z_D2x_ab+ABY*I_NAI_Ix5z_D2x_ab;
  Double I_NAI_I6y_F2xy_ab = I_NAI_K7y_D2x_ab+ABY*I_NAI_I6y_D2x_ab;
  Double I_NAI_I5yz_F2xy_ab = I_NAI_K6yz_D2x_ab+ABY*I_NAI_I5yz_D2x_ab;
  Double I_NAI_I4y2z_F2xy_ab = I_NAI_K5y2z_D2x_ab+ABY*I_NAI_I4y2z_D2x_ab;
  Double I_NAI_I3y3z_F2xy_ab = I_NAI_K4y3z_D2x_ab+ABY*I_NAI_I3y3z_D2x_ab;
  Double I_NAI_I2y4z_F2xy_ab = I_NAI_K3y4z_D2x_ab+ABY*I_NAI_I2y4z_D2x_ab;
  Double I_NAI_Iy5z_F2xy_ab = I_NAI_K2y5z_D2x_ab+ABY*I_NAI_Iy5z_D2x_ab;
  Double I_NAI_I6z_F2xy_ab = I_NAI_Ky6z_D2x_ab+ABY*I_NAI_I6z_D2x_ab;
  Double I_NAI_I5xz_F2xz_ab = I_NAI_K5x2z_D2x_ab+ABZ*I_NAI_I5xz_D2x_ab;
  Double I_NAI_I4xyz_F2xz_ab = I_NAI_K4xy2z_D2x_ab+ABZ*I_NAI_I4xyz_D2x_ab;
  Double I_NAI_I4x2z_F2xz_ab = I_NAI_K4x3z_D2x_ab+ABZ*I_NAI_I4x2z_D2x_ab;
  Double I_NAI_I3x2yz_F2xz_ab = I_NAI_K3x2y2z_D2x_ab+ABZ*I_NAI_I3x2yz_D2x_ab;
  Double I_NAI_I3xy2z_F2xz_ab = I_NAI_K3xy3z_D2x_ab+ABZ*I_NAI_I3xy2z_D2x_ab;
  Double I_NAI_I3x3z_F2xz_ab = I_NAI_K3x4z_D2x_ab+ABZ*I_NAI_I3x3z_D2x_ab;
  Double I_NAI_I2x3yz_F2xz_ab = I_NAI_K2x3y2z_D2x_ab+ABZ*I_NAI_I2x3yz_D2x_ab;
  Double I_NAI_I2x2y2z_F2xz_ab = I_NAI_K2x2y3z_D2x_ab+ABZ*I_NAI_I2x2y2z_D2x_ab;
  Double I_NAI_I2xy3z_F2xz_ab = I_NAI_K2xy4z_D2x_ab+ABZ*I_NAI_I2xy3z_D2x_ab;
  Double I_NAI_I2x4z_F2xz_ab = I_NAI_K2x5z_D2x_ab+ABZ*I_NAI_I2x4z_D2x_ab;
  Double I_NAI_Ix4yz_F2xz_ab = I_NAI_Kx4y2z_D2x_ab+ABZ*I_NAI_Ix4yz_D2x_ab;
  Double I_NAI_Ix3y2z_F2xz_ab = I_NAI_Kx3y3z_D2x_ab+ABZ*I_NAI_Ix3y2z_D2x_ab;
  Double I_NAI_Ix2y3z_F2xz_ab = I_NAI_Kx2y4z_D2x_ab+ABZ*I_NAI_Ix2y3z_D2x_ab;
  Double I_NAI_Ixy4z_F2xz_ab = I_NAI_Kxy5z_D2x_ab+ABZ*I_NAI_Ixy4z_D2x_ab;
  Double I_NAI_Ix5z_F2xz_ab = I_NAI_Kx6z_D2x_ab+ABZ*I_NAI_Ix5z_D2x_ab;
  Double I_NAI_I5yz_F2xz_ab = I_NAI_K5y2z_D2x_ab+ABZ*I_NAI_I5yz_D2x_ab;
  Double I_NAI_I4y2z_F2xz_ab = I_NAI_K4y3z_D2x_ab+ABZ*I_NAI_I4y2z_D2x_ab;
  Double I_NAI_I3y3z_F2xz_ab = I_NAI_K3y4z_D2x_ab+ABZ*I_NAI_I3y3z_D2x_ab;
  Double I_NAI_I2y4z_F2xz_ab = I_NAI_K2y5z_D2x_ab+ABZ*I_NAI_I2y4z_D2x_ab;
  Double I_NAI_Iy5z_F2xz_ab = I_NAI_Ky6z_D2x_ab+ABZ*I_NAI_Iy5z_D2x_ab;
  Double I_NAI_I6z_F2xz_ab = I_NAI_K7z_D2x_ab+ABZ*I_NAI_I6z_D2x_ab;
  Double I_NAI_I5xz_Fx2y_ab = I_NAI_K6xz_D2y_ab+ABX*I_NAI_I5xz_D2y_ab;
  Double I_NAI_I4xyz_Fx2y_ab = I_NAI_K5xyz_D2y_ab+ABX*I_NAI_I4xyz_D2y_ab;
  Double I_NAI_I4x2z_Fx2y_ab = I_NAI_K5x2z_D2y_ab+ABX*I_NAI_I4x2z_D2y_ab;
  Double I_NAI_I3x2yz_Fx2y_ab = I_NAI_K4x2yz_D2y_ab+ABX*I_NAI_I3x2yz_D2y_ab;
  Double I_NAI_I3xy2z_Fx2y_ab = I_NAI_K4xy2z_D2y_ab+ABX*I_NAI_I3xy2z_D2y_ab;
  Double I_NAI_I3x3z_Fx2y_ab = I_NAI_K4x3z_D2y_ab+ABX*I_NAI_I3x3z_D2y_ab;
  Double I_NAI_I2x3yz_Fx2y_ab = I_NAI_K3x3yz_D2y_ab+ABX*I_NAI_I2x3yz_D2y_ab;
  Double I_NAI_I2x2y2z_Fx2y_ab = I_NAI_K3x2y2z_D2y_ab+ABX*I_NAI_I2x2y2z_D2y_ab;
  Double I_NAI_I2xy3z_Fx2y_ab = I_NAI_K3xy3z_D2y_ab+ABX*I_NAI_I2xy3z_D2y_ab;
  Double I_NAI_I2x4z_Fx2y_ab = I_NAI_K3x4z_D2y_ab+ABX*I_NAI_I2x4z_D2y_ab;
  Double I_NAI_Ix4yz_Fx2y_ab = I_NAI_K2x4yz_D2y_ab+ABX*I_NAI_Ix4yz_D2y_ab;
  Double I_NAI_Ix3y2z_Fx2y_ab = I_NAI_K2x3y2z_D2y_ab+ABX*I_NAI_Ix3y2z_D2y_ab;
  Double I_NAI_Ix2y3z_Fx2y_ab = I_NAI_K2x2y3z_D2y_ab+ABX*I_NAI_Ix2y3z_D2y_ab;
  Double I_NAI_Ixy4z_Fx2y_ab = I_NAI_K2xy4z_D2y_ab+ABX*I_NAI_Ixy4z_D2y_ab;
  Double I_NAI_Ix5z_Fx2y_ab = I_NAI_K2x5z_D2y_ab+ABX*I_NAI_Ix5z_D2y_ab;
  Double I_NAI_I5yz_Fx2y_ab = I_NAI_Kx5yz_D2y_ab+ABX*I_NAI_I5yz_D2y_ab;
  Double I_NAI_I4y2z_Fx2y_ab = I_NAI_Kx4y2z_D2y_ab+ABX*I_NAI_I4y2z_D2y_ab;
  Double I_NAI_I3y3z_Fx2y_ab = I_NAI_Kx3y3z_D2y_ab+ABX*I_NAI_I3y3z_D2y_ab;
  Double I_NAI_I2y4z_Fx2y_ab = I_NAI_Kx2y4z_D2y_ab+ABX*I_NAI_I2y4z_D2y_ab;
  Double I_NAI_Iy5z_Fx2y_ab = I_NAI_Kxy5z_D2y_ab+ABX*I_NAI_Iy5z_D2y_ab;
  Double I_NAI_I6z_Fx2y_ab = I_NAI_Kx6z_D2y_ab+ABX*I_NAI_I6z_D2y_ab;
  Double I_NAI_I5xy_Fx2z_ab = I_NAI_K6xy_D2z_ab+ABX*I_NAI_I5xy_D2z_ab;
  Double I_NAI_I4x2y_Fx2z_ab = I_NAI_K5x2y_D2z_ab+ABX*I_NAI_I4x2y_D2z_ab;
  Double I_NAI_I4xyz_Fx2z_ab = I_NAI_K5xyz_D2z_ab+ABX*I_NAI_I4xyz_D2z_ab;
  Double I_NAI_I3x3y_Fx2z_ab = I_NAI_K4x3y_D2z_ab+ABX*I_NAI_I3x3y_D2z_ab;
  Double I_NAI_I3x2yz_Fx2z_ab = I_NAI_K4x2yz_D2z_ab+ABX*I_NAI_I3x2yz_D2z_ab;
  Double I_NAI_I3xy2z_Fx2z_ab = I_NAI_K4xy2z_D2z_ab+ABX*I_NAI_I3xy2z_D2z_ab;
  Double I_NAI_I2x4y_Fx2z_ab = I_NAI_K3x4y_D2z_ab+ABX*I_NAI_I2x4y_D2z_ab;
  Double I_NAI_I2x3yz_Fx2z_ab = I_NAI_K3x3yz_D2z_ab+ABX*I_NAI_I2x3yz_D2z_ab;
  Double I_NAI_I2x2y2z_Fx2z_ab = I_NAI_K3x2y2z_D2z_ab+ABX*I_NAI_I2x2y2z_D2z_ab;
  Double I_NAI_I2xy3z_Fx2z_ab = I_NAI_K3xy3z_D2z_ab+ABX*I_NAI_I2xy3z_D2z_ab;
  Double I_NAI_Ix5y_Fx2z_ab = I_NAI_K2x5y_D2z_ab+ABX*I_NAI_Ix5y_D2z_ab;
  Double I_NAI_Ix4yz_Fx2z_ab = I_NAI_K2x4yz_D2z_ab+ABX*I_NAI_Ix4yz_D2z_ab;
  Double I_NAI_Ix3y2z_Fx2z_ab = I_NAI_K2x3y2z_D2z_ab+ABX*I_NAI_Ix3y2z_D2z_ab;
  Double I_NAI_Ix2y3z_Fx2z_ab = I_NAI_K2x2y3z_D2z_ab+ABX*I_NAI_Ix2y3z_D2z_ab;
  Double I_NAI_Ixy4z_Fx2z_ab = I_NAI_K2xy4z_D2z_ab+ABX*I_NAI_Ixy4z_D2z_ab;
  Double I_NAI_I6y_Fx2z_ab = I_NAI_Kx6y_D2z_ab+ABX*I_NAI_I6y_D2z_ab;
  Double I_NAI_I5yz_Fx2z_ab = I_NAI_Kx5yz_D2z_ab+ABX*I_NAI_I5yz_D2z_ab;
  Double I_NAI_I4y2z_Fx2z_ab = I_NAI_Kx4y2z_D2z_ab+ABX*I_NAI_I4y2z_D2z_ab;
  Double I_NAI_I3y3z_Fx2z_ab = I_NAI_Kx3y3z_D2z_ab+ABX*I_NAI_I3y3z_D2z_ab;
  Double I_NAI_I2y4z_Fx2z_ab = I_NAI_Kx2y4z_D2z_ab+ABX*I_NAI_I2y4z_D2z_ab;
  Double I_NAI_Iy5z_Fx2z_ab = I_NAI_Kxy5z_D2z_ab+ABX*I_NAI_Iy5z_D2z_ab;
  Double I_NAI_I6x_F3y_ab = I_NAI_K6xy_D2y_ab+ABY*I_NAI_I6x_D2y_ab;
  Double I_NAI_I5xy_F3y_ab = I_NAI_K5x2y_D2y_ab+ABY*I_NAI_I5xy_D2y_ab;
  Double I_NAI_I5xz_F3y_ab = I_NAI_K5xyz_D2y_ab+ABY*I_NAI_I5xz_D2y_ab;
  Double I_NAI_I4x2y_F3y_ab = I_NAI_K4x3y_D2y_ab+ABY*I_NAI_I4x2y_D2y_ab;
  Double I_NAI_I4xyz_F3y_ab = I_NAI_K4x2yz_D2y_ab+ABY*I_NAI_I4xyz_D2y_ab;
  Double I_NAI_I4x2z_F3y_ab = I_NAI_K4xy2z_D2y_ab+ABY*I_NAI_I4x2z_D2y_ab;
  Double I_NAI_I3x3y_F3y_ab = I_NAI_K3x4y_D2y_ab+ABY*I_NAI_I3x3y_D2y_ab;
  Double I_NAI_I3x2yz_F3y_ab = I_NAI_K3x3yz_D2y_ab+ABY*I_NAI_I3x2yz_D2y_ab;
  Double I_NAI_I3xy2z_F3y_ab = I_NAI_K3x2y2z_D2y_ab+ABY*I_NAI_I3xy2z_D2y_ab;
  Double I_NAI_I3x3z_F3y_ab = I_NAI_K3xy3z_D2y_ab+ABY*I_NAI_I3x3z_D2y_ab;
  Double I_NAI_I2x4y_F3y_ab = I_NAI_K2x5y_D2y_ab+ABY*I_NAI_I2x4y_D2y_ab;
  Double I_NAI_I2x3yz_F3y_ab = I_NAI_K2x4yz_D2y_ab+ABY*I_NAI_I2x3yz_D2y_ab;
  Double I_NAI_I2x2y2z_F3y_ab = I_NAI_K2x3y2z_D2y_ab+ABY*I_NAI_I2x2y2z_D2y_ab;
  Double I_NAI_I2xy3z_F3y_ab = I_NAI_K2x2y3z_D2y_ab+ABY*I_NAI_I2xy3z_D2y_ab;
  Double I_NAI_I2x4z_F3y_ab = I_NAI_K2xy4z_D2y_ab+ABY*I_NAI_I2x4z_D2y_ab;
  Double I_NAI_Ix5y_F3y_ab = I_NAI_Kx6y_D2y_ab+ABY*I_NAI_Ix5y_D2y_ab;
  Double I_NAI_Ix4yz_F3y_ab = I_NAI_Kx5yz_D2y_ab+ABY*I_NAI_Ix4yz_D2y_ab;
  Double I_NAI_Ix3y2z_F3y_ab = I_NAI_Kx4y2z_D2y_ab+ABY*I_NAI_Ix3y2z_D2y_ab;
  Double I_NAI_Ix2y3z_F3y_ab = I_NAI_Kx3y3z_D2y_ab+ABY*I_NAI_Ix2y3z_D2y_ab;
  Double I_NAI_Ixy4z_F3y_ab = I_NAI_Kx2y4z_D2y_ab+ABY*I_NAI_Ixy4z_D2y_ab;
  Double I_NAI_Ix5z_F3y_ab = I_NAI_Kxy5z_D2y_ab+ABY*I_NAI_Ix5z_D2y_ab;
  Double I_NAI_I6y_F3y_ab = I_NAI_K7y_D2y_ab+ABY*I_NAI_I6y_D2y_ab;
  Double I_NAI_I5yz_F3y_ab = I_NAI_K6yz_D2y_ab+ABY*I_NAI_I5yz_D2y_ab;
  Double I_NAI_I4y2z_F3y_ab = I_NAI_K5y2z_D2y_ab+ABY*I_NAI_I4y2z_D2y_ab;
  Double I_NAI_I3y3z_F3y_ab = I_NAI_K4y3z_D2y_ab+ABY*I_NAI_I3y3z_D2y_ab;
  Double I_NAI_I2y4z_F3y_ab = I_NAI_K3y4z_D2y_ab+ABY*I_NAI_I2y4z_D2y_ab;
  Double I_NAI_Iy5z_F3y_ab = I_NAI_K2y5z_D2y_ab+ABY*I_NAI_Iy5z_D2y_ab;
  Double I_NAI_I6z_F3y_ab = I_NAI_Ky6z_D2y_ab+ABY*I_NAI_I6z_D2y_ab;
  Double I_NAI_I5xz_F2yz_ab = I_NAI_K5x2z_D2y_ab+ABZ*I_NAI_I5xz_D2y_ab;
  Double I_NAI_I4xyz_F2yz_ab = I_NAI_K4xy2z_D2y_ab+ABZ*I_NAI_I4xyz_D2y_ab;
  Double I_NAI_I4x2z_F2yz_ab = I_NAI_K4x3z_D2y_ab+ABZ*I_NAI_I4x2z_D2y_ab;
  Double I_NAI_I3x2yz_F2yz_ab = I_NAI_K3x2y2z_D2y_ab+ABZ*I_NAI_I3x2yz_D2y_ab;
  Double I_NAI_I3xy2z_F2yz_ab = I_NAI_K3xy3z_D2y_ab+ABZ*I_NAI_I3xy2z_D2y_ab;
  Double I_NAI_I3x3z_F2yz_ab = I_NAI_K3x4z_D2y_ab+ABZ*I_NAI_I3x3z_D2y_ab;
  Double I_NAI_I2x3yz_F2yz_ab = I_NAI_K2x3y2z_D2y_ab+ABZ*I_NAI_I2x3yz_D2y_ab;
  Double I_NAI_I2x2y2z_F2yz_ab = I_NAI_K2x2y3z_D2y_ab+ABZ*I_NAI_I2x2y2z_D2y_ab;
  Double I_NAI_I2xy3z_F2yz_ab = I_NAI_K2xy4z_D2y_ab+ABZ*I_NAI_I2xy3z_D2y_ab;
  Double I_NAI_I2x4z_F2yz_ab = I_NAI_K2x5z_D2y_ab+ABZ*I_NAI_I2x4z_D2y_ab;
  Double I_NAI_Ix4yz_F2yz_ab = I_NAI_Kx4y2z_D2y_ab+ABZ*I_NAI_Ix4yz_D2y_ab;
  Double I_NAI_Ix3y2z_F2yz_ab = I_NAI_Kx3y3z_D2y_ab+ABZ*I_NAI_Ix3y2z_D2y_ab;
  Double I_NAI_Ix2y3z_F2yz_ab = I_NAI_Kx2y4z_D2y_ab+ABZ*I_NAI_Ix2y3z_D2y_ab;
  Double I_NAI_Ixy4z_F2yz_ab = I_NAI_Kxy5z_D2y_ab+ABZ*I_NAI_Ixy4z_D2y_ab;
  Double I_NAI_Ix5z_F2yz_ab = I_NAI_Kx6z_D2y_ab+ABZ*I_NAI_Ix5z_D2y_ab;
  Double I_NAI_I5yz_F2yz_ab = I_NAI_K5y2z_D2y_ab+ABZ*I_NAI_I5yz_D2y_ab;
  Double I_NAI_I4y2z_F2yz_ab = I_NAI_K4y3z_D2y_ab+ABZ*I_NAI_I4y2z_D2y_ab;
  Double I_NAI_I3y3z_F2yz_ab = I_NAI_K3y4z_D2y_ab+ABZ*I_NAI_I3y3z_D2y_ab;
  Double I_NAI_I2y4z_F2yz_ab = I_NAI_K2y5z_D2y_ab+ABZ*I_NAI_I2y4z_D2y_ab;
  Double I_NAI_Iy5z_F2yz_ab = I_NAI_Ky6z_D2y_ab+ABZ*I_NAI_Iy5z_D2y_ab;
  Double I_NAI_I6z_F2yz_ab = I_NAI_K7z_D2y_ab+ABZ*I_NAI_I6z_D2y_ab;
  Double I_NAI_I6x_F3z_ab = I_NAI_K6xz_D2z_ab+ABZ*I_NAI_I6x_D2z_ab;
  Double I_NAI_I5xy_F3z_ab = I_NAI_K5xyz_D2z_ab+ABZ*I_NAI_I5xy_D2z_ab;
  Double I_NAI_I5xz_F3z_ab = I_NAI_K5x2z_D2z_ab+ABZ*I_NAI_I5xz_D2z_ab;
  Double I_NAI_I4x2y_F3z_ab = I_NAI_K4x2yz_D2z_ab+ABZ*I_NAI_I4x2y_D2z_ab;
  Double I_NAI_I4xyz_F3z_ab = I_NAI_K4xy2z_D2z_ab+ABZ*I_NAI_I4xyz_D2z_ab;
  Double I_NAI_I4x2z_F3z_ab = I_NAI_K4x3z_D2z_ab+ABZ*I_NAI_I4x2z_D2z_ab;
  Double I_NAI_I3x3y_F3z_ab = I_NAI_K3x3yz_D2z_ab+ABZ*I_NAI_I3x3y_D2z_ab;
  Double I_NAI_I3x2yz_F3z_ab = I_NAI_K3x2y2z_D2z_ab+ABZ*I_NAI_I3x2yz_D2z_ab;
  Double I_NAI_I3xy2z_F3z_ab = I_NAI_K3xy3z_D2z_ab+ABZ*I_NAI_I3xy2z_D2z_ab;
  Double I_NAI_I3x3z_F3z_ab = I_NAI_K3x4z_D2z_ab+ABZ*I_NAI_I3x3z_D2z_ab;
  Double I_NAI_I2x4y_F3z_ab = I_NAI_K2x4yz_D2z_ab+ABZ*I_NAI_I2x4y_D2z_ab;
  Double I_NAI_I2x3yz_F3z_ab = I_NAI_K2x3y2z_D2z_ab+ABZ*I_NAI_I2x3yz_D2z_ab;
  Double I_NAI_I2x2y2z_F3z_ab = I_NAI_K2x2y3z_D2z_ab+ABZ*I_NAI_I2x2y2z_D2z_ab;
  Double I_NAI_I2xy3z_F3z_ab = I_NAI_K2xy4z_D2z_ab+ABZ*I_NAI_I2xy3z_D2z_ab;
  Double I_NAI_I2x4z_F3z_ab = I_NAI_K2x5z_D2z_ab+ABZ*I_NAI_I2x4z_D2z_ab;
  Double I_NAI_Ix5y_F3z_ab = I_NAI_Kx5yz_D2z_ab+ABZ*I_NAI_Ix5y_D2z_ab;
  Double I_NAI_Ix4yz_F3z_ab = I_NAI_Kx4y2z_D2z_ab+ABZ*I_NAI_Ix4yz_D2z_ab;
  Double I_NAI_Ix3y2z_F3z_ab = I_NAI_Kx3y3z_D2z_ab+ABZ*I_NAI_Ix3y2z_D2z_ab;
  Double I_NAI_Ix2y3z_F3z_ab = I_NAI_Kx2y4z_D2z_ab+ABZ*I_NAI_Ix2y3z_D2z_ab;
  Double I_NAI_Ixy4z_F3z_ab = I_NAI_Kxy5z_D2z_ab+ABZ*I_NAI_Ixy4z_D2z_ab;
  Double I_NAI_Ix5z_F3z_ab = I_NAI_Kx6z_D2z_ab+ABZ*I_NAI_Ix5z_D2z_ab;
  Double I_NAI_I6y_F3z_ab = I_NAI_K6yz_D2z_ab+ABZ*I_NAI_I6y_D2z_ab;
  Double I_NAI_I5yz_F3z_ab = I_NAI_K5y2z_D2z_ab+ABZ*I_NAI_I5yz_D2z_ab;
  Double I_NAI_I4y2z_F3z_ab = I_NAI_K4y3z_D2z_ab+ABZ*I_NAI_I4y2z_D2z_ab;
  Double I_NAI_I3y3z_F3z_ab = I_NAI_K3y4z_D2z_ab+ABZ*I_NAI_I3y3z_D2z_ab;
  Double I_NAI_I2y4z_F3z_ab = I_NAI_K2y5z_D2z_ab+ABZ*I_NAI_I2y4z_D2z_ab;
  Double I_NAI_Iy5z_F3z_ab = I_NAI_Ky6z_D2z_ab+ABZ*I_NAI_Iy5z_D2z_ab;
  Double I_NAI_I6z_F3z_ab = I_NAI_K7z_D2z_ab+ABZ*I_NAI_I6z_D2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_ab
   * RHS shell quartet name: SQ_NAI_H_F_ab
   ************************************************************/
  Double I_NAI_H5x_G4x_ab = I_NAI_I6x_F3x_ab+ABX*I_NAI_H5x_F3x_ab;
  Double I_NAI_H4xy_G4x_ab = I_NAI_I5xy_F3x_ab+ABX*I_NAI_H4xy_F3x_ab;
  Double I_NAI_H4xz_G4x_ab = I_NAI_I5xz_F3x_ab+ABX*I_NAI_H4xz_F3x_ab;
  Double I_NAI_H3x2y_G4x_ab = I_NAI_I4x2y_F3x_ab+ABX*I_NAI_H3x2y_F3x_ab;
  Double I_NAI_H3xyz_G4x_ab = I_NAI_I4xyz_F3x_ab+ABX*I_NAI_H3xyz_F3x_ab;
  Double I_NAI_H3x2z_G4x_ab = I_NAI_I4x2z_F3x_ab+ABX*I_NAI_H3x2z_F3x_ab;
  Double I_NAI_H2x3y_G4x_ab = I_NAI_I3x3y_F3x_ab+ABX*I_NAI_H2x3y_F3x_ab;
  Double I_NAI_H2x2yz_G4x_ab = I_NAI_I3x2yz_F3x_ab+ABX*I_NAI_H2x2yz_F3x_ab;
  Double I_NAI_H2xy2z_G4x_ab = I_NAI_I3xy2z_F3x_ab+ABX*I_NAI_H2xy2z_F3x_ab;
  Double I_NAI_H2x3z_G4x_ab = I_NAI_I3x3z_F3x_ab+ABX*I_NAI_H2x3z_F3x_ab;
  Double I_NAI_Hx4y_G4x_ab = I_NAI_I2x4y_F3x_ab+ABX*I_NAI_Hx4y_F3x_ab;
  Double I_NAI_Hx3yz_G4x_ab = I_NAI_I2x3yz_F3x_ab+ABX*I_NAI_Hx3yz_F3x_ab;
  Double I_NAI_Hx2y2z_G4x_ab = I_NAI_I2x2y2z_F3x_ab+ABX*I_NAI_Hx2y2z_F3x_ab;
  Double I_NAI_Hxy3z_G4x_ab = I_NAI_I2xy3z_F3x_ab+ABX*I_NAI_Hxy3z_F3x_ab;
  Double I_NAI_Hx4z_G4x_ab = I_NAI_I2x4z_F3x_ab+ABX*I_NAI_Hx4z_F3x_ab;
  Double I_NAI_H5y_G4x_ab = I_NAI_Ix5y_F3x_ab+ABX*I_NAI_H5y_F3x_ab;
  Double I_NAI_H4yz_G4x_ab = I_NAI_Ix4yz_F3x_ab+ABX*I_NAI_H4yz_F3x_ab;
  Double I_NAI_H3y2z_G4x_ab = I_NAI_Ix3y2z_F3x_ab+ABX*I_NAI_H3y2z_F3x_ab;
  Double I_NAI_H2y3z_G4x_ab = I_NAI_Ix2y3z_F3x_ab+ABX*I_NAI_H2y3z_F3x_ab;
  Double I_NAI_Hy4z_G4x_ab = I_NAI_Ixy4z_F3x_ab+ABX*I_NAI_Hy4z_F3x_ab;
  Double I_NAI_H5z_G4x_ab = I_NAI_Ix5z_F3x_ab+ABX*I_NAI_H5z_F3x_ab;
  Double I_NAI_H5x_G3xy_ab = I_NAI_I5xy_F3x_ab+ABY*I_NAI_H5x_F3x_ab;
  Double I_NAI_H4xy_G3xy_ab = I_NAI_I4x2y_F3x_ab+ABY*I_NAI_H4xy_F3x_ab;
  Double I_NAI_H4xz_G3xy_ab = I_NAI_I4xyz_F3x_ab+ABY*I_NAI_H4xz_F3x_ab;
  Double I_NAI_H3x2y_G3xy_ab = I_NAI_I3x3y_F3x_ab+ABY*I_NAI_H3x2y_F3x_ab;
  Double I_NAI_H3xyz_G3xy_ab = I_NAI_I3x2yz_F3x_ab+ABY*I_NAI_H3xyz_F3x_ab;
  Double I_NAI_H3x2z_G3xy_ab = I_NAI_I3xy2z_F3x_ab+ABY*I_NAI_H3x2z_F3x_ab;
  Double I_NAI_H2x3y_G3xy_ab = I_NAI_I2x4y_F3x_ab+ABY*I_NAI_H2x3y_F3x_ab;
  Double I_NAI_H2x2yz_G3xy_ab = I_NAI_I2x3yz_F3x_ab+ABY*I_NAI_H2x2yz_F3x_ab;
  Double I_NAI_H2xy2z_G3xy_ab = I_NAI_I2x2y2z_F3x_ab+ABY*I_NAI_H2xy2z_F3x_ab;
  Double I_NAI_H2x3z_G3xy_ab = I_NAI_I2xy3z_F3x_ab+ABY*I_NAI_H2x3z_F3x_ab;
  Double I_NAI_Hx4y_G3xy_ab = I_NAI_Ix5y_F3x_ab+ABY*I_NAI_Hx4y_F3x_ab;
  Double I_NAI_Hx3yz_G3xy_ab = I_NAI_Ix4yz_F3x_ab+ABY*I_NAI_Hx3yz_F3x_ab;
  Double I_NAI_Hx2y2z_G3xy_ab = I_NAI_Ix3y2z_F3x_ab+ABY*I_NAI_Hx2y2z_F3x_ab;
  Double I_NAI_Hxy3z_G3xy_ab = I_NAI_Ix2y3z_F3x_ab+ABY*I_NAI_Hxy3z_F3x_ab;
  Double I_NAI_Hx4z_G3xy_ab = I_NAI_Ixy4z_F3x_ab+ABY*I_NAI_Hx4z_F3x_ab;
  Double I_NAI_H5y_G3xy_ab = I_NAI_I6y_F3x_ab+ABY*I_NAI_H5y_F3x_ab;
  Double I_NAI_H4yz_G3xy_ab = I_NAI_I5yz_F3x_ab+ABY*I_NAI_H4yz_F3x_ab;
  Double I_NAI_H3y2z_G3xy_ab = I_NAI_I4y2z_F3x_ab+ABY*I_NAI_H3y2z_F3x_ab;
  Double I_NAI_H2y3z_G3xy_ab = I_NAI_I3y3z_F3x_ab+ABY*I_NAI_H2y3z_F3x_ab;
  Double I_NAI_Hy4z_G3xy_ab = I_NAI_I2y4z_F3x_ab+ABY*I_NAI_Hy4z_F3x_ab;
  Double I_NAI_H5z_G3xy_ab = I_NAI_Iy5z_F3x_ab+ABY*I_NAI_H5z_F3x_ab;
  Double I_NAI_H5x_G3xz_ab = I_NAI_I5xz_F3x_ab+ABZ*I_NAI_H5x_F3x_ab;
  Double I_NAI_H4xy_G3xz_ab = I_NAI_I4xyz_F3x_ab+ABZ*I_NAI_H4xy_F3x_ab;
  Double I_NAI_H4xz_G3xz_ab = I_NAI_I4x2z_F3x_ab+ABZ*I_NAI_H4xz_F3x_ab;
  Double I_NAI_H3x2y_G3xz_ab = I_NAI_I3x2yz_F3x_ab+ABZ*I_NAI_H3x2y_F3x_ab;
  Double I_NAI_H3xyz_G3xz_ab = I_NAI_I3xy2z_F3x_ab+ABZ*I_NAI_H3xyz_F3x_ab;
  Double I_NAI_H3x2z_G3xz_ab = I_NAI_I3x3z_F3x_ab+ABZ*I_NAI_H3x2z_F3x_ab;
  Double I_NAI_H2x3y_G3xz_ab = I_NAI_I2x3yz_F3x_ab+ABZ*I_NAI_H2x3y_F3x_ab;
  Double I_NAI_H2x2yz_G3xz_ab = I_NAI_I2x2y2z_F3x_ab+ABZ*I_NAI_H2x2yz_F3x_ab;
  Double I_NAI_H2xy2z_G3xz_ab = I_NAI_I2xy3z_F3x_ab+ABZ*I_NAI_H2xy2z_F3x_ab;
  Double I_NAI_H2x3z_G3xz_ab = I_NAI_I2x4z_F3x_ab+ABZ*I_NAI_H2x3z_F3x_ab;
  Double I_NAI_Hx4y_G3xz_ab = I_NAI_Ix4yz_F3x_ab+ABZ*I_NAI_Hx4y_F3x_ab;
  Double I_NAI_Hx3yz_G3xz_ab = I_NAI_Ix3y2z_F3x_ab+ABZ*I_NAI_Hx3yz_F3x_ab;
  Double I_NAI_Hx2y2z_G3xz_ab = I_NAI_Ix2y3z_F3x_ab+ABZ*I_NAI_Hx2y2z_F3x_ab;
  Double I_NAI_Hxy3z_G3xz_ab = I_NAI_Ixy4z_F3x_ab+ABZ*I_NAI_Hxy3z_F3x_ab;
  Double I_NAI_Hx4z_G3xz_ab = I_NAI_Ix5z_F3x_ab+ABZ*I_NAI_Hx4z_F3x_ab;
  Double I_NAI_H5y_G3xz_ab = I_NAI_I5yz_F3x_ab+ABZ*I_NAI_H5y_F3x_ab;
  Double I_NAI_H4yz_G3xz_ab = I_NAI_I4y2z_F3x_ab+ABZ*I_NAI_H4yz_F3x_ab;
  Double I_NAI_H3y2z_G3xz_ab = I_NAI_I3y3z_F3x_ab+ABZ*I_NAI_H3y2z_F3x_ab;
  Double I_NAI_H2y3z_G3xz_ab = I_NAI_I2y4z_F3x_ab+ABZ*I_NAI_H2y3z_F3x_ab;
  Double I_NAI_Hy4z_G3xz_ab = I_NAI_Iy5z_F3x_ab+ABZ*I_NAI_Hy4z_F3x_ab;
  Double I_NAI_H5z_G3xz_ab = I_NAI_I6z_F3x_ab+ABZ*I_NAI_H5z_F3x_ab;
  Double I_NAI_H5x_G2x2y_ab = I_NAI_I5xy_F2xy_ab+ABY*I_NAI_H5x_F2xy_ab;
  Double I_NAI_H4xy_G2x2y_ab = I_NAI_I4x2y_F2xy_ab+ABY*I_NAI_H4xy_F2xy_ab;
  Double I_NAI_H4xz_G2x2y_ab = I_NAI_I4xyz_F2xy_ab+ABY*I_NAI_H4xz_F2xy_ab;
  Double I_NAI_H3x2y_G2x2y_ab = I_NAI_I3x3y_F2xy_ab+ABY*I_NAI_H3x2y_F2xy_ab;
  Double I_NAI_H3xyz_G2x2y_ab = I_NAI_I3x2yz_F2xy_ab+ABY*I_NAI_H3xyz_F2xy_ab;
  Double I_NAI_H3x2z_G2x2y_ab = I_NAI_I3xy2z_F2xy_ab+ABY*I_NAI_H3x2z_F2xy_ab;
  Double I_NAI_H2x3y_G2x2y_ab = I_NAI_I2x4y_F2xy_ab+ABY*I_NAI_H2x3y_F2xy_ab;
  Double I_NAI_H2x2yz_G2x2y_ab = I_NAI_I2x3yz_F2xy_ab+ABY*I_NAI_H2x2yz_F2xy_ab;
  Double I_NAI_H2xy2z_G2x2y_ab = I_NAI_I2x2y2z_F2xy_ab+ABY*I_NAI_H2xy2z_F2xy_ab;
  Double I_NAI_H2x3z_G2x2y_ab = I_NAI_I2xy3z_F2xy_ab+ABY*I_NAI_H2x3z_F2xy_ab;
  Double I_NAI_Hx4y_G2x2y_ab = I_NAI_Ix5y_F2xy_ab+ABY*I_NAI_Hx4y_F2xy_ab;
  Double I_NAI_Hx3yz_G2x2y_ab = I_NAI_Ix4yz_F2xy_ab+ABY*I_NAI_Hx3yz_F2xy_ab;
  Double I_NAI_Hx2y2z_G2x2y_ab = I_NAI_Ix3y2z_F2xy_ab+ABY*I_NAI_Hx2y2z_F2xy_ab;
  Double I_NAI_Hxy3z_G2x2y_ab = I_NAI_Ix2y3z_F2xy_ab+ABY*I_NAI_Hxy3z_F2xy_ab;
  Double I_NAI_Hx4z_G2x2y_ab = I_NAI_Ixy4z_F2xy_ab+ABY*I_NAI_Hx4z_F2xy_ab;
  Double I_NAI_H5y_G2x2y_ab = I_NAI_I6y_F2xy_ab+ABY*I_NAI_H5y_F2xy_ab;
  Double I_NAI_H4yz_G2x2y_ab = I_NAI_I5yz_F2xy_ab+ABY*I_NAI_H4yz_F2xy_ab;
  Double I_NAI_H3y2z_G2x2y_ab = I_NAI_I4y2z_F2xy_ab+ABY*I_NAI_H3y2z_F2xy_ab;
  Double I_NAI_H2y3z_G2x2y_ab = I_NAI_I3y3z_F2xy_ab+ABY*I_NAI_H2y3z_F2xy_ab;
  Double I_NAI_Hy4z_G2x2y_ab = I_NAI_I2y4z_F2xy_ab+ABY*I_NAI_Hy4z_F2xy_ab;
  Double I_NAI_H5z_G2x2y_ab = I_NAI_Iy5z_F2xy_ab+ABY*I_NAI_H5z_F2xy_ab;
  Double I_NAI_H5x_G2xyz_ab = I_NAI_I5xz_F2xy_ab+ABZ*I_NAI_H5x_F2xy_ab;
  Double I_NAI_H4xy_G2xyz_ab = I_NAI_I4xyz_F2xy_ab+ABZ*I_NAI_H4xy_F2xy_ab;
  Double I_NAI_H4xz_G2xyz_ab = I_NAI_I4x2z_F2xy_ab+ABZ*I_NAI_H4xz_F2xy_ab;
  Double I_NAI_H3x2y_G2xyz_ab = I_NAI_I3x2yz_F2xy_ab+ABZ*I_NAI_H3x2y_F2xy_ab;
  Double I_NAI_H3xyz_G2xyz_ab = I_NAI_I3xy2z_F2xy_ab+ABZ*I_NAI_H3xyz_F2xy_ab;
  Double I_NAI_H3x2z_G2xyz_ab = I_NAI_I3x3z_F2xy_ab+ABZ*I_NAI_H3x2z_F2xy_ab;
  Double I_NAI_H2x3y_G2xyz_ab = I_NAI_I2x3yz_F2xy_ab+ABZ*I_NAI_H2x3y_F2xy_ab;
  Double I_NAI_H2x2yz_G2xyz_ab = I_NAI_I2x2y2z_F2xy_ab+ABZ*I_NAI_H2x2yz_F2xy_ab;
  Double I_NAI_H2xy2z_G2xyz_ab = I_NAI_I2xy3z_F2xy_ab+ABZ*I_NAI_H2xy2z_F2xy_ab;
  Double I_NAI_H2x3z_G2xyz_ab = I_NAI_I2x4z_F2xy_ab+ABZ*I_NAI_H2x3z_F2xy_ab;
  Double I_NAI_Hx4y_G2xyz_ab = I_NAI_Ix4yz_F2xy_ab+ABZ*I_NAI_Hx4y_F2xy_ab;
  Double I_NAI_Hx3yz_G2xyz_ab = I_NAI_Ix3y2z_F2xy_ab+ABZ*I_NAI_Hx3yz_F2xy_ab;
  Double I_NAI_Hx2y2z_G2xyz_ab = I_NAI_Ix2y3z_F2xy_ab+ABZ*I_NAI_Hx2y2z_F2xy_ab;
  Double I_NAI_Hxy3z_G2xyz_ab = I_NAI_Ixy4z_F2xy_ab+ABZ*I_NAI_Hxy3z_F2xy_ab;
  Double I_NAI_Hx4z_G2xyz_ab = I_NAI_Ix5z_F2xy_ab+ABZ*I_NAI_Hx4z_F2xy_ab;
  Double I_NAI_H5y_G2xyz_ab = I_NAI_I5yz_F2xy_ab+ABZ*I_NAI_H5y_F2xy_ab;
  Double I_NAI_H4yz_G2xyz_ab = I_NAI_I4y2z_F2xy_ab+ABZ*I_NAI_H4yz_F2xy_ab;
  Double I_NAI_H3y2z_G2xyz_ab = I_NAI_I3y3z_F2xy_ab+ABZ*I_NAI_H3y2z_F2xy_ab;
  Double I_NAI_H2y3z_G2xyz_ab = I_NAI_I2y4z_F2xy_ab+ABZ*I_NAI_H2y3z_F2xy_ab;
  Double I_NAI_Hy4z_G2xyz_ab = I_NAI_Iy5z_F2xy_ab+ABZ*I_NAI_Hy4z_F2xy_ab;
  Double I_NAI_H5z_G2xyz_ab = I_NAI_I6z_F2xy_ab+ABZ*I_NAI_H5z_F2xy_ab;
  Double I_NAI_H5x_G2x2z_ab = I_NAI_I5xz_F2xz_ab+ABZ*I_NAI_H5x_F2xz_ab;
  Double I_NAI_H4xy_G2x2z_ab = I_NAI_I4xyz_F2xz_ab+ABZ*I_NAI_H4xy_F2xz_ab;
  Double I_NAI_H4xz_G2x2z_ab = I_NAI_I4x2z_F2xz_ab+ABZ*I_NAI_H4xz_F2xz_ab;
  Double I_NAI_H3x2y_G2x2z_ab = I_NAI_I3x2yz_F2xz_ab+ABZ*I_NAI_H3x2y_F2xz_ab;
  Double I_NAI_H3xyz_G2x2z_ab = I_NAI_I3xy2z_F2xz_ab+ABZ*I_NAI_H3xyz_F2xz_ab;
  Double I_NAI_H3x2z_G2x2z_ab = I_NAI_I3x3z_F2xz_ab+ABZ*I_NAI_H3x2z_F2xz_ab;
  Double I_NAI_H2x3y_G2x2z_ab = I_NAI_I2x3yz_F2xz_ab+ABZ*I_NAI_H2x3y_F2xz_ab;
  Double I_NAI_H2x2yz_G2x2z_ab = I_NAI_I2x2y2z_F2xz_ab+ABZ*I_NAI_H2x2yz_F2xz_ab;
  Double I_NAI_H2xy2z_G2x2z_ab = I_NAI_I2xy3z_F2xz_ab+ABZ*I_NAI_H2xy2z_F2xz_ab;
  Double I_NAI_H2x3z_G2x2z_ab = I_NAI_I2x4z_F2xz_ab+ABZ*I_NAI_H2x3z_F2xz_ab;
  Double I_NAI_Hx4y_G2x2z_ab = I_NAI_Ix4yz_F2xz_ab+ABZ*I_NAI_Hx4y_F2xz_ab;
  Double I_NAI_Hx3yz_G2x2z_ab = I_NAI_Ix3y2z_F2xz_ab+ABZ*I_NAI_Hx3yz_F2xz_ab;
  Double I_NAI_Hx2y2z_G2x2z_ab = I_NAI_Ix2y3z_F2xz_ab+ABZ*I_NAI_Hx2y2z_F2xz_ab;
  Double I_NAI_Hxy3z_G2x2z_ab = I_NAI_Ixy4z_F2xz_ab+ABZ*I_NAI_Hxy3z_F2xz_ab;
  Double I_NAI_Hx4z_G2x2z_ab = I_NAI_Ix5z_F2xz_ab+ABZ*I_NAI_Hx4z_F2xz_ab;
  Double I_NAI_H5y_G2x2z_ab = I_NAI_I5yz_F2xz_ab+ABZ*I_NAI_H5y_F2xz_ab;
  Double I_NAI_H4yz_G2x2z_ab = I_NAI_I4y2z_F2xz_ab+ABZ*I_NAI_H4yz_F2xz_ab;
  Double I_NAI_H3y2z_G2x2z_ab = I_NAI_I3y3z_F2xz_ab+ABZ*I_NAI_H3y2z_F2xz_ab;
  Double I_NAI_H2y3z_G2x2z_ab = I_NAI_I2y4z_F2xz_ab+ABZ*I_NAI_H2y3z_F2xz_ab;
  Double I_NAI_Hy4z_G2x2z_ab = I_NAI_Iy5z_F2xz_ab+ABZ*I_NAI_Hy4z_F2xz_ab;
  Double I_NAI_H5z_G2x2z_ab = I_NAI_I6z_F2xz_ab+ABZ*I_NAI_H5z_F2xz_ab;
  Double I_NAI_H5x_Gx3y_ab = I_NAI_I6x_F3y_ab+ABX*I_NAI_H5x_F3y_ab;
  Double I_NAI_H4xy_Gx3y_ab = I_NAI_I5xy_F3y_ab+ABX*I_NAI_H4xy_F3y_ab;
  Double I_NAI_H4xz_Gx3y_ab = I_NAI_I5xz_F3y_ab+ABX*I_NAI_H4xz_F3y_ab;
  Double I_NAI_H3x2y_Gx3y_ab = I_NAI_I4x2y_F3y_ab+ABX*I_NAI_H3x2y_F3y_ab;
  Double I_NAI_H3xyz_Gx3y_ab = I_NAI_I4xyz_F3y_ab+ABX*I_NAI_H3xyz_F3y_ab;
  Double I_NAI_H3x2z_Gx3y_ab = I_NAI_I4x2z_F3y_ab+ABX*I_NAI_H3x2z_F3y_ab;
  Double I_NAI_H2x3y_Gx3y_ab = I_NAI_I3x3y_F3y_ab+ABX*I_NAI_H2x3y_F3y_ab;
  Double I_NAI_H2x2yz_Gx3y_ab = I_NAI_I3x2yz_F3y_ab+ABX*I_NAI_H2x2yz_F3y_ab;
  Double I_NAI_H2xy2z_Gx3y_ab = I_NAI_I3xy2z_F3y_ab+ABX*I_NAI_H2xy2z_F3y_ab;
  Double I_NAI_H2x3z_Gx3y_ab = I_NAI_I3x3z_F3y_ab+ABX*I_NAI_H2x3z_F3y_ab;
  Double I_NAI_Hx4y_Gx3y_ab = I_NAI_I2x4y_F3y_ab+ABX*I_NAI_Hx4y_F3y_ab;
  Double I_NAI_Hx3yz_Gx3y_ab = I_NAI_I2x3yz_F3y_ab+ABX*I_NAI_Hx3yz_F3y_ab;
  Double I_NAI_Hx2y2z_Gx3y_ab = I_NAI_I2x2y2z_F3y_ab+ABX*I_NAI_Hx2y2z_F3y_ab;
  Double I_NAI_Hxy3z_Gx3y_ab = I_NAI_I2xy3z_F3y_ab+ABX*I_NAI_Hxy3z_F3y_ab;
  Double I_NAI_Hx4z_Gx3y_ab = I_NAI_I2x4z_F3y_ab+ABX*I_NAI_Hx4z_F3y_ab;
  Double I_NAI_H5y_Gx3y_ab = I_NAI_Ix5y_F3y_ab+ABX*I_NAI_H5y_F3y_ab;
  Double I_NAI_H4yz_Gx3y_ab = I_NAI_Ix4yz_F3y_ab+ABX*I_NAI_H4yz_F3y_ab;
  Double I_NAI_H3y2z_Gx3y_ab = I_NAI_Ix3y2z_F3y_ab+ABX*I_NAI_H3y2z_F3y_ab;
  Double I_NAI_H2y3z_Gx3y_ab = I_NAI_Ix2y3z_F3y_ab+ABX*I_NAI_H2y3z_F3y_ab;
  Double I_NAI_Hy4z_Gx3y_ab = I_NAI_Ixy4z_F3y_ab+ABX*I_NAI_Hy4z_F3y_ab;
  Double I_NAI_H5z_Gx3y_ab = I_NAI_Ix5z_F3y_ab+ABX*I_NAI_H5z_F3y_ab;
  Double I_NAI_H5x_Gx2yz_ab = I_NAI_I5xz_Fx2y_ab+ABZ*I_NAI_H5x_Fx2y_ab;
  Double I_NAI_H4xy_Gx2yz_ab = I_NAI_I4xyz_Fx2y_ab+ABZ*I_NAI_H4xy_Fx2y_ab;
  Double I_NAI_H4xz_Gx2yz_ab = I_NAI_I4x2z_Fx2y_ab+ABZ*I_NAI_H4xz_Fx2y_ab;
  Double I_NAI_H3x2y_Gx2yz_ab = I_NAI_I3x2yz_Fx2y_ab+ABZ*I_NAI_H3x2y_Fx2y_ab;
  Double I_NAI_H3xyz_Gx2yz_ab = I_NAI_I3xy2z_Fx2y_ab+ABZ*I_NAI_H3xyz_Fx2y_ab;
  Double I_NAI_H3x2z_Gx2yz_ab = I_NAI_I3x3z_Fx2y_ab+ABZ*I_NAI_H3x2z_Fx2y_ab;
  Double I_NAI_H2x3y_Gx2yz_ab = I_NAI_I2x3yz_Fx2y_ab+ABZ*I_NAI_H2x3y_Fx2y_ab;
  Double I_NAI_H2x2yz_Gx2yz_ab = I_NAI_I2x2y2z_Fx2y_ab+ABZ*I_NAI_H2x2yz_Fx2y_ab;
  Double I_NAI_H2xy2z_Gx2yz_ab = I_NAI_I2xy3z_Fx2y_ab+ABZ*I_NAI_H2xy2z_Fx2y_ab;
  Double I_NAI_H2x3z_Gx2yz_ab = I_NAI_I2x4z_Fx2y_ab+ABZ*I_NAI_H2x3z_Fx2y_ab;
  Double I_NAI_Hx4y_Gx2yz_ab = I_NAI_Ix4yz_Fx2y_ab+ABZ*I_NAI_Hx4y_Fx2y_ab;
  Double I_NAI_Hx3yz_Gx2yz_ab = I_NAI_Ix3y2z_Fx2y_ab+ABZ*I_NAI_Hx3yz_Fx2y_ab;
  Double I_NAI_Hx2y2z_Gx2yz_ab = I_NAI_Ix2y3z_Fx2y_ab+ABZ*I_NAI_Hx2y2z_Fx2y_ab;
  Double I_NAI_Hxy3z_Gx2yz_ab = I_NAI_Ixy4z_Fx2y_ab+ABZ*I_NAI_Hxy3z_Fx2y_ab;
  Double I_NAI_Hx4z_Gx2yz_ab = I_NAI_Ix5z_Fx2y_ab+ABZ*I_NAI_Hx4z_Fx2y_ab;
  Double I_NAI_H5y_Gx2yz_ab = I_NAI_I5yz_Fx2y_ab+ABZ*I_NAI_H5y_Fx2y_ab;
  Double I_NAI_H4yz_Gx2yz_ab = I_NAI_I4y2z_Fx2y_ab+ABZ*I_NAI_H4yz_Fx2y_ab;
  Double I_NAI_H3y2z_Gx2yz_ab = I_NAI_I3y3z_Fx2y_ab+ABZ*I_NAI_H3y2z_Fx2y_ab;
  Double I_NAI_H2y3z_Gx2yz_ab = I_NAI_I2y4z_Fx2y_ab+ABZ*I_NAI_H2y3z_Fx2y_ab;
  Double I_NAI_Hy4z_Gx2yz_ab = I_NAI_Iy5z_Fx2y_ab+ABZ*I_NAI_Hy4z_Fx2y_ab;
  Double I_NAI_H5z_Gx2yz_ab = I_NAI_I6z_Fx2y_ab+ABZ*I_NAI_H5z_Fx2y_ab;
  Double I_NAI_H5x_Gxy2z_ab = I_NAI_I5xy_Fx2z_ab+ABY*I_NAI_H5x_Fx2z_ab;
  Double I_NAI_H4xy_Gxy2z_ab = I_NAI_I4x2y_Fx2z_ab+ABY*I_NAI_H4xy_Fx2z_ab;
  Double I_NAI_H4xz_Gxy2z_ab = I_NAI_I4xyz_Fx2z_ab+ABY*I_NAI_H4xz_Fx2z_ab;
  Double I_NAI_H3x2y_Gxy2z_ab = I_NAI_I3x3y_Fx2z_ab+ABY*I_NAI_H3x2y_Fx2z_ab;
  Double I_NAI_H3xyz_Gxy2z_ab = I_NAI_I3x2yz_Fx2z_ab+ABY*I_NAI_H3xyz_Fx2z_ab;
  Double I_NAI_H3x2z_Gxy2z_ab = I_NAI_I3xy2z_Fx2z_ab+ABY*I_NAI_H3x2z_Fx2z_ab;
  Double I_NAI_H2x3y_Gxy2z_ab = I_NAI_I2x4y_Fx2z_ab+ABY*I_NAI_H2x3y_Fx2z_ab;
  Double I_NAI_H2x2yz_Gxy2z_ab = I_NAI_I2x3yz_Fx2z_ab+ABY*I_NAI_H2x2yz_Fx2z_ab;
  Double I_NAI_H2xy2z_Gxy2z_ab = I_NAI_I2x2y2z_Fx2z_ab+ABY*I_NAI_H2xy2z_Fx2z_ab;
  Double I_NAI_H2x3z_Gxy2z_ab = I_NAI_I2xy3z_Fx2z_ab+ABY*I_NAI_H2x3z_Fx2z_ab;
  Double I_NAI_Hx4y_Gxy2z_ab = I_NAI_Ix5y_Fx2z_ab+ABY*I_NAI_Hx4y_Fx2z_ab;
  Double I_NAI_Hx3yz_Gxy2z_ab = I_NAI_Ix4yz_Fx2z_ab+ABY*I_NAI_Hx3yz_Fx2z_ab;
  Double I_NAI_Hx2y2z_Gxy2z_ab = I_NAI_Ix3y2z_Fx2z_ab+ABY*I_NAI_Hx2y2z_Fx2z_ab;
  Double I_NAI_Hxy3z_Gxy2z_ab = I_NAI_Ix2y3z_Fx2z_ab+ABY*I_NAI_Hxy3z_Fx2z_ab;
  Double I_NAI_Hx4z_Gxy2z_ab = I_NAI_Ixy4z_Fx2z_ab+ABY*I_NAI_Hx4z_Fx2z_ab;
  Double I_NAI_H5y_Gxy2z_ab = I_NAI_I6y_Fx2z_ab+ABY*I_NAI_H5y_Fx2z_ab;
  Double I_NAI_H4yz_Gxy2z_ab = I_NAI_I5yz_Fx2z_ab+ABY*I_NAI_H4yz_Fx2z_ab;
  Double I_NAI_H3y2z_Gxy2z_ab = I_NAI_I4y2z_Fx2z_ab+ABY*I_NAI_H3y2z_Fx2z_ab;
  Double I_NAI_H2y3z_Gxy2z_ab = I_NAI_I3y3z_Fx2z_ab+ABY*I_NAI_H2y3z_Fx2z_ab;
  Double I_NAI_Hy4z_Gxy2z_ab = I_NAI_I2y4z_Fx2z_ab+ABY*I_NAI_Hy4z_Fx2z_ab;
  Double I_NAI_H5z_Gxy2z_ab = I_NAI_Iy5z_Fx2z_ab+ABY*I_NAI_H5z_Fx2z_ab;
  Double I_NAI_H5x_Gx3z_ab = I_NAI_I6x_F3z_ab+ABX*I_NAI_H5x_F3z_ab;
  Double I_NAI_H4xy_Gx3z_ab = I_NAI_I5xy_F3z_ab+ABX*I_NAI_H4xy_F3z_ab;
  Double I_NAI_H4xz_Gx3z_ab = I_NAI_I5xz_F3z_ab+ABX*I_NAI_H4xz_F3z_ab;
  Double I_NAI_H3x2y_Gx3z_ab = I_NAI_I4x2y_F3z_ab+ABX*I_NAI_H3x2y_F3z_ab;
  Double I_NAI_H3xyz_Gx3z_ab = I_NAI_I4xyz_F3z_ab+ABX*I_NAI_H3xyz_F3z_ab;
  Double I_NAI_H3x2z_Gx3z_ab = I_NAI_I4x2z_F3z_ab+ABX*I_NAI_H3x2z_F3z_ab;
  Double I_NAI_H2x3y_Gx3z_ab = I_NAI_I3x3y_F3z_ab+ABX*I_NAI_H2x3y_F3z_ab;
  Double I_NAI_H2x2yz_Gx3z_ab = I_NAI_I3x2yz_F3z_ab+ABX*I_NAI_H2x2yz_F3z_ab;
  Double I_NAI_H2xy2z_Gx3z_ab = I_NAI_I3xy2z_F3z_ab+ABX*I_NAI_H2xy2z_F3z_ab;
  Double I_NAI_H2x3z_Gx3z_ab = I_NAI_I3x3z_F3z_ab+ABX*I_NAI_H2x3z_F3z_ab;
  Double I_NAI_Hx4y_Gx3z_ab = I_NAI_I2x4y_F3z_ab+ABX*I_NAI_Hx4y_F3z_ab;
  Double I_NAI_Hx3yz_Gx3z_ab = I_NAI_I2x3yz_F3z_ab+ABX*I_NAI_Hx3yz_F3z_ab;
  Double I_NAI_Hx2y2z_Gx3z_ab = I_NAI_I2x2y2z_F3z_ab+ABX*I_NAI_Hx2y2z_F3z_ab;
  Double I_NAI_Hxy3z_Gx3z_ab = I_NAI_I2xy3z_F3z_ab+ABX*I_NAI_Hxy3z_F3z_ab;
  Double I_NAI_Hx4z_Gx3z_ab = I_NAI_I2x4z_F3z_ab+ABX*I_NAI_Hx4z_F3z_ab;
  Double I_NAI_H5y_Gx3z_ab = I_NAI_Ix5y_F3z_ab+ABX*I_NAI_H5y_F3z_ab;
  Double I_NAI_H4yz_Gx3z_ab = I_NAI_Ix4yz_F3z_ab+ABX*I_NAI_H4yz_F3z_ab;
  Double I_NAI_H3y2z_Gx3z_ab = I_NAI_Ix3y2z_F3z_ab+ABX*I_NAI_H3y2z_F3z_ab;
  Double I_NAI_H2y3z_Gx3z_ab = I_NAI_Ix2y3z_F3z_ab+ABX*I_NAI_H2y3z_F3z_ab;
  Double I_NAI_Hy4z_Gx3z_ab = I_NAI_Ixy4z_F3z_ab+ABX*I_NAI_Hy4z_F3z_ab;
  Double I_NAI_H5z_Gx3z_ab = I_NAI_Ix5z_F3z_ab+ABX*I_NAI_H5z_F3z_ab;
  Double I_NAI_H5x_G4y_ab = I_NAI_I5xy_F3y_ab+ABY*I_NAI_H5x_F3y_ab;
  Double I_NAI_H4xy_G4y_ab = I_NAI_I4x2y_F3y_ab+ABY*I_NAI_H4xy_F3y_ab;
  Double I_NAI_H4xz_G4y_ab = I_NAI_I4xyz_F3y_ab+ABY*I_NAI_H4xz_F3y_ab;
  Double I_NAI_H3x2y_G4y_ab = I_NAI_I3x3y_F3y_ab+ABY*I_NAI_H3x2y_F3y_ab;
  Double I_NAI_H3xyz_G4y_ab = I_NAI_I3x2yz_F3y_ab+ABY*I_NAI_H3xyz_F3y_ab;
  Double I_NAI_H3x2z_G4y_ab = I_NAI_I3xy2z_F3y_ab+ABY*I_NAI_H3x2z_F3y_ab;
  Double I_NAI_H2x3y_G4y_ab = I_NAI_I2x4y_F3y_ab+ABY*I_NAI_H2x3y_F3y_ab;
  Double I_NAI_H2x2yz_G4y_ab = I_NAI_I2x3yz_F3y_ab+ABY*I_NAI_H2x2yz_F3y_ab;
  Double I_NAI_H2xy2z_G4y_ab = I_NAI_I2x2y2z_F3y_ab+ABY*I_NAI_H2xy2z_F3y_ab;
  Double I_NAI_H2x3z_G4y_ab = I_NAI_I2xy3z_F3y_ab+ABY*I_NAI_H2x3z_F3y_ab;
  Double I_NAI_Hx4y_G4y_ab = I_NAI_Ix5y_F3y_ab+ABY*I_NAI_Hx4y_F3y_ab;
  Double I_NAI_Hx3yz_G4y_ab = I_NAI_Ix4yz_F3y_ab+ABY*I_NAI_Hx3yz_F3y_ab;
  Double I_NAI_Hx2y2z_G4y_ab = I_NAI_Ix3y2z_F3y_ab+ABY*I_NAI_Hx2y2z_F3y_ab;
  Double I_NAI_Hxy3z_G4y_ab = I_NAI_Ix2y3z_F3y_ab+ABY*I_NAI_Hxy3z_F3y_ab;
  Double I_NAI_Hx4z_G4y_ab = I_NAI_Ixy4z_F3y_ab+ABY*I_NAI_Hx4z_F3y_ab;
  Double I_NAI_H5y_G4y_ab = I_NAI_I6y_F3y_ab+ABY*I_NAI_H5y_F3y_ab;
  Double I_NAI_H4yz_G4y_ab = I_NAI_I5yz_F3y_ab+ABY*I_NAI_H4yz_F3y_ab;
  Double I_NAI_H3y2z_G4y_ab = I_NAI_I4y2z_F3y_ab+ABY*I_NAI_H3y2z_F3y_ab;
  Double I_NAI_H2y3z_G4y_ab = I_NAI_I3y3z_F3y_ab+ABY*I_NAI_H2y3z_F3y_ab;
  Double I_NAI_Hy4z_G4y_ab = I_NAI_I2y4z_F3y_ab+ABY*I_NAI_Hy4z_F3y_ab;
  Double I_NAI_H5z_G4y_ab = I_NAI_Iy5z_F3y_ab+ABY*I_NAI_H5z_F3y_ab;
  Double I_NAI_H5x_G3yz_ab = I_NAI_I5xz_F3y_ab+ABZ*I_NAI_H5x_F3y_ab;
  Double I_NAI_H4xy_G3yz_ab = I_NAI_I4xyz_F3y_ab+ABZ*I_NAI_H4xy_F3y_ab;
  Double I_NAI_H4xz_G3yz_ab = I_NAI_I4x2z_F3y_ab+ABZ*I_NAI_H4xz_F3y_ab;
  Double I_NAI_H3x2y_G3yz_ab = I_NAI_I3x2yz_F3y_ab+ABZ*I_NAI_H3x2y_F3y_ab;
  Double I_NAI_H3xyz_G3yz_ab = I_NAI_I3xy2z_F3y_ab+ABZ*I_NAI_H3xyz_F3y_ab;
  Double I_NAI_H3x2z_G3yz_ab = I_NAI_I3x3z_F3y_ab+ABZ*I_NAI_H3x2z_F3y_ab;
  Double I_NAI_H2x3y_G3yz_ab = I_NAI_I2x3yz_F3y_ab+ABZ*I_NAI_H2x3y_F3y_ab;
  Double I_NAI_H2x2yz_G3yz_ab = I_NAI_I2x2y2z_F3y_ab+ABZ*I_NAI_H2x2yz_F3y_ab;
  Double I_NAI_H2xy2z_G3yz_ab = I_NAI_I2xy3z_F3y_ab+ABZ*I_NAI_H2xy2z_F3y_ab;
  Double I_NAI_H2x3z_G3yz_ab = I_NAI_I2x4z_F3y_ab+ABZ*I_NAI_H2x3z_F3y_ab;
  Double I_NAI_Hx4y_G3yz_ab = I_NAI_Ix4yz_F3y_ab+ABZ*I_NAI_Hx4y_F3y_ab;
  Double I_NAI_Hx3yz_G3yz_ab = I_NAI_Ix3y2z_F3y_ab+ABZ*I_NAI_Hx3yz_F3y_ab;
  Double I_NAI_Hx2y2z_G3yz_ab = I_NAI_Ix2y3z_F3y_ab+ABZ*I_NAI_Hx2y2z_F3y_ab;
  Double I_NAI_Hxy3z_G3yz_ab = I_NAI_Ixy4z_F3y_ab+ABZ*I_NAI_Hxy3z_F3y_ab;
  Double I_NAI_Hx4z_G3yz_ab = I_NAI_Ix5z_F3y_ab+ABZ*I_NAI_Hx4z_F3y_ab;
  Double I_NAI_H5y_G3yz_ab = I_NAI_I5yz_F3y_ab+ABZ*I_NAI_H5y_F3y_ab;
  Double I_NAI_H4yz_G3yz_ab = I_NAI_I4y2z_F3y_ab+ABZ*I_NAI_H4yz_F3y_ab;
  Double I_NAI_H3y2z_G3yz_ab = I_NAI_I3y3z_F3y_ab+ABZ*I_NAI_H3y2z_F3y_ab;
  Double I_NAI_H2y3z_G3yz_ab = I_NAI_I2y4z_F3y_ab+ABZ*I_NAI_H2y3z_F3y_ab;
  Double I_NAI_Hy4z_G3yz_ab = I_NAI_Iy5z_F3y_ab+ABZ*I_NAI_Hy4z_F3y_ab;
  Double I_NAI_H5z_G3yz_ab = I_NAI_I6z_F3y_ab+ABZ*I_NAI_H5z_F3y_ab;
  Double I_NAI_H5x_G2y2z_ab = I_NAI_I5xz_F2yz_ab+ABZ*I_NAI_H5x_F2yz_ab;
  Double I_NAI_H4xy_G2y2z_ab = I_NAI_I4xyz_F2yz_ab+ABZ*I_NAI_H4xy_F2yz_ab;
  Double I_NAI_H4xz_G2y2z_ab = I_NAI_I4x2z_F2yz_ab+ABZ*I_NAI_H4xz_F2yz_ab;
  Double I_NAI_H3x2y_G2y2z_ab = I_NAI_I3x2yz_F2yz_ab+ABZ*I_NAI_H3x2y_F2yz_ab;
  Double I_NAI_H3xyz_G2y2z_ab = I_NAI_I3xy2z_F2yz_ab+ABZ*I_NAI_H3xyz_F2yz_ab;
  Double I_NAI_H3x2z_G2y2z_ab = I_NAI_I3x3z_F2yz_ab+ABZ*I_NAI_H3x2z_F2yz_ab;
  Double I_NAI_H2x3y_G2y2z_ab = I_NAI_I2x3yz_F2yz_ab+ABZ*I_NAI_H2x3y_F2yz_ab;
  Double I_NAI_H2x2yz_G2y2z_ab = I_NAI_I2x2y2z_F2yz_ab+ABZ*I_NAI_H2x2yz_F2yz_ab;
  Double I_NAI_H2xy2z_G2y2z_ab = I_NAI_I2xy3z_F2yz_ab+ABZ*I_NAI_H2xy2z_F2yz_ab;
  Double I_NAI_H2x3z_G2y2z_ab = I_NAI_I2x4z_F2yz_ab+ABZ*I_NAI_H2x3z_F2yz_ab;
  Double I_NAI_Hx4y_G2y2z_ab = I_NAI_Ix4yz_F2yz_ab+ABZ*I_NAI_Hx4y_F2yz_ab;
  Double I_NAI_Hx3yz_G2y2z_ab = I_NAI_Ix3y2z_F2yz_ab+ABZ*I_NAI_Hx3yz_F2yz_ab;
  Double I_NAI_Hx2y2z_G2y2z_ab = I_NAI_Ix2y3z_F2yz_ab+ABZ*I_NAI_Hx2y2z_F2yz_ab;
  Double I_NAI_Hxy3z_G2y2z_ab = I_NAI_Ixy4z_F2yz_ab+ABZ*I_NAI_Hxy3z_F2yz_ab;
  Double I_NAI_Hx4z_G2y2z_ab = I_NAI_Ix5z_F2yz_ab+ABZ*I_NAI_Hx4z_F2yz_ab;
  Double I_NAI_H5y_G2y2z_ab = I_NAI_I5yz_F2yz_ab+ABZ*I_NAI_H5y_F2yz_ab;
  Double I_NAI_H4yz_G2y2z_ab = I_NAI_I4y2z_F2yz_ab+ABZ*I_NAI_H4yz_F2yz_ab;
  Double I_NAI_H3y2z_G2y2z_ab = I_NAI_I3y3z_F2yz_ab+ABZ*I_NAI_H3y2z_F2yz_ab;
  Double I_NAI_H2y3z_G2y2z_ab = I_NAI_I2y4z_F2yz_ab+ABZ*I_NAI_H2y3z_F2yz_ab;
  Double I_NAI_Hy4z_G2y2z_ab = I_NAI_Iy5z_F2yz_ab+ABZ*I_NAI_Hy4z_F2yz_ab;
  Double I_NAI_H5z_G2y2z_ab = I_NAI_I6z_F2yz_ab+ABZ*I_NAI_H5z_F2yz_ab;
  Double I_NAI_H5x_Gy3z_ab = I_NAI_I5xy_F3z_ab+ABY*I_NAI_H5x_F3z_ab;
  Double I_NAI_H4xy_Gy3z_ab = I_NAI_I4x2y_F3z_ab+ABY*I_NAI_H4xy_F3z_ab;
  Double I_NAI_H4xz_Gy3z_ab = I_NAI_I4xyz_F3z_ab+ABY*I_NAI_H4xz_F3z_ab;
  Double I_NAI_H3x2y_Gy3z_ab = I_NAI_I3x3y_F3z_ab+ABY*I_NAI_H3x2y_F3z_ab;
  Double I_NAI_H3xyz_Gy3z_ab = I_NAI_I3x2yz_F3z_ab+ABY*I_NAI_H3xyz_F3z_ab;
  Double I_NAI_H3x2z_Gy3z_ab = I_NAI_I3xy2z_F3z_ab+ABY*I_NAI_H3x2z_F3z_ab;
  Double I_NAI_H2x3y_Gy3z_ab = I_NAI_I2x4y_F3z_ab+ABY*I_NAI_H2x3y_F3z_ab;
  Double I_NAI_H2x2yz_Gy3z_ab = I_NAI_I2x3yz_F3z_ab+ABY*I_NAI_H2x2yz_F3z_ab;
  Double I_NAI_H2xy2z_Gy3z_ab = I_NAI_I2x2y2z_F3z_ab+ABY*I_NAI_H2xy2z_F3z_ab;
  Double I_NAI_H2x3z_Gy3z_ab = I_NAI_I2xy3z_F3z_ab+ABY*I_NAI_H2x3z_F3z_ab;
  Double I_NAI_Hx4y_Gy3z_ab = I_NAI_Ix5y_F3z_ab+ABY*I_NAI_Hx4y_F3z_ab;
  Double I_NAI_Hx3yz_Gy3z_ab = I_NAI_Ix4yz_F3z_ab+ABY*I_NAI_Hx3yz_F3z_ab;
  Double I_NAI_Hx2y2z_Gy3z_ab = I_NAI_Ix3y2z_F3z_ab+ABY*I_NAI_Hx2y2z_F3z_ab;
  Double I_NAI_Hxy3z_Gy3z_ab = I_NAI_Ix2y3z_F3z_ab+ABY*I_NAI_Hxy3z_F3z_ab;
  Double I_NAI_Hx4z_Gy3z_ab = I_NAI_Ixy4z_F3z_ab+ABY*I_NAI_Hx4z_F3z_ab;
  Double I_NAI_H5y_Gy3z_ab = I_NAI_I6y_F3z_ab+ABY*I_NAI_H5y_F3z_ab;
  Double I_NAI_H4yz_Gy3z_ab = I_NAI_I5yz_F3z_ab+ABY*I_NAI_H4yz_F3z_ab;
  Double I_NAI_H3y2z_Gy3z_ab = I_NAI_I4y2z_F3z_ab+ABY*I_NAI_H3y2z_F3z_ab;
  Double I_NAI_H2y3z_Gy3z_ab = I_NAI_I3y3z_F3z_ab+ABY*I_NAI_H2y3z_F3z_ab;
  Double I_NAI_Hy4z_Gy3z_ab = I_NAI_I2y4z_F3z_ab+ABY*I_NAI_Hy4z_F3z_ab;
  Double I_NAI_H5z_Gy3z_ab = I_NAI_Iy5z_F3z_ab+ABY*I_NAI_H5z_F3z_ab;
  Double I_NAI_H5x_G4z_ab = I_NAI_I5xz_F3z_ab+ABZ*I_NAI_H5x_F3z_ab;
  Double I_NAI_H4xy_G4z_ab = I_NAI_I4xyz_F3z_ab+ABZ*I_NAI_H4xy_F3z_ab;
  Double I_NAI_H4xz_G4z_ab = I_NAI_I4x2z_F3z_ab+ABZ*I_NAI_H4xz_F3z_ab;
  Double I_NAI_H3x2y_G4z_ab = I_NAI_I3x2yz_F3z_ab+ABZ*I_NAI_H3x2y_F3z_ab;
  Double I_NAI_H3xyz_G4z_ab = I_NAI_I3xy2z_F3z_ab+ABZ*I_NAI_H3xyz_F3z_ab;
  Double I_NAI_H3x2z_G4z_ab = I_NAI_I3x3z_F3z_ab+ABZ*I_NAI_H3x2z_F3z_ab;
  Double I_NAI_H2x3y_G4z_ab = I_NAI_I2x3yz_F3z_ab+ABZ*I_NAI_H2x3y_F3z_ab;
  Double I_NAI_H2x2yz_G4z_ab = I_NAI_I2x2y2z_F3z_ab+ABZ*I_NAI_H2x2yz_F3z_ab;
  Double I_NAI_H2xy2z_G4z_ab = I_NAI_I2xy3z_F3z_ab+ABZ*I_NAI_H2xy2z_F3z_ab;
  Double I_NAI_H2x3z_G4z_ab = I_NAI_I2x4z_F3z_ab+ABZ*I_NAI_H2x3z_F3z_ab;
  Double I_NAI_Hx4y_G4z_ab = I_NAI_Ix4yz_F3z_ab+ABZ*I_NAI_Hx4y_F3z_ab;
  Double I_NAI_Hx3yz_G4z_ab = I_NAI_Ix3y2z_F3z_ab+ABZ*I_NAI_Hx3yz_F3z_ab;
  Double I_NAI_Hx2y2z_G4z_ab = I_NAI_Ix2y3z_F3z_ab+ABZ*I_NAI_Hx2y2z_F3z_ab;
  Double I_NAI_Hxy3z_G4z_ab = I_NAI_Ixy4z_F3z_ab+ABZ*I_NAI_Hxy3z_F3z_ab;
  Double I_NAI_Hx4z_G4z_ab = I_NAI_Ix5z_F3z_ab+ABZ*I_NAI_Hx4z_F3z_ab;
  Double I_NAI_H5y_G4z_ab = I_NAI_I5yz_F3z_ab+ABZ*I_NAI_H5y_F3z_ab;
  Double I_NAI_H4yz_G4z_ab = I_NAI_I4y2z_F3z_ab+ABZ*I_NAI_H4yz_F3z_ab;
  Double I_NAI_H3y2z_G4z_ab = I_NAI_I3y3z_F3z_ab+ABZ*I_NAI_H3y2z_F3z_ab;
  Double I_NAI_H2y3z_G4z_ab = I_NAI_I2y4z_F3z_ab+ABZ*I_NAI_H2y3z_F3z_ab;
  Double I_NAI_Hy4z_G4z_ab = I_NAI_Iy5z_F3z_ab+ABZ*I_NAI_Hy4z_F3z_ab;
  Double I_NAI_H5z_G4z_ab = I_NAI_I6z_F3z_ab+ABZ*I_NAI_H5z_F3z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_bb
   * RHS shell quartet name: SQ_NAI_G_S_bb
   ************************************************************/
  Double I_NAI_G4x_Px_bb = I_NAI_H5x_S_bb+ABX*I_NAI_G4x_S_bb;
  Double I_NAI_G3xy_Px_bb = I_NAI_H4xy_S_bb+ABX*I_NAI_G3xy_S_bb;
  Double I_NAI_G3xz_Px_bb = I_NAI_H4xz_S_bb+ABX*I_NAI_G3xz_S_bb;
  Double I_NAI_G2x2y_Px_bb = I_NAI_H3x2y_S_bb+ABX*I_NAI_G2x2y_S_bb;
  Double I_NAI_G2xyz_Px_bb = I_NAI_H3xyz_S_bb+ABX*I_NAI_G2xyz_S_bb;
  Double I_NAI_G2x2z_Px_bb = I_NAI_H3x2z_S_bb+ABX*I_NAI_G2x2z_S_bb;
  Double I_NAI_Gx3y_Px_bb = I_NAI_H2x3y_S_bb+ABX*I_NAI_Gx3y_S_bb;
  Double I_NAI_Gx2yz_Px_bb = I_NAI_H2x2yz_S_bb+ABX*I_NAI_Gx2yz_S_bb;
  Double I_NAI_Gxy2z_Px_bb = I_NAI_H2xy2z_S_bb+ABX*I_NAI_Gxy2z_S_bb;
  Double I_NAI_Gx3z_Px_bb = I_NAI_H2x3z_S_bb+ABX*I_NAI_Gx3z_S_bb;
  Double I_NAI_G4y_Px_bb = I_NAI_Hx4y_S_bb+ABX*I_NAI_G4y_S_bb;
  Double I_NAI_G3yz_Px_bb = I_NAI_Hx3yz_S_bb+ABX*I_NAI_G3yz_S_bb;
  Double I_NAI_G2y2z_Px_bb = I_NAI_Hx2y2z_S_bb+ABX*I_NAI_G2y2z_S_bb;
  Double I_NAI_Gy3z_Px_bb = I_NAI_Hxy3z_S_bb+ABX*I_NAI_Gy3z_S_bb;
  Double I_NAI_G4z_Px_bb = I_NAI_Hx4z_S_bb+ABX*I_NAI_G4z_S_bb;
  Double I_NAI_G4x_Py_bb = I_NAI_H4xy_S_bb+ABY*I_NAI_G4x_S_bb;
  Double I_NAI_G3xy_Py_bb = I_NAI_H3x2y_S_bb+ABY*I_NAI_G3xy_S_bb;
  Double I_NAI_G3xz_Py_bb = I_NAI_H3xyz_S_bb+ABY*I_NAI_G3xz_S_bb;
  Double I_NAI_G2x2y_Py_bb = I_NAI_H2x3y_S_bb+ABY*I_NAI_G2x2y_S_bb;
  Double I_NAI_G2xyz_Py_bb = I_NAI_H2x2yz_S_bb+ABY*I_NAI_G2xyz_S_bb;
  Double I_NAI_G2x2z_Py_bb = I_NAI_H2xy2z_S_bb+ABY*I_NAI_G2x2z_S_bb;
  Double I_NAI_Gx3y_Py_bb = I_NAI_Hx4y_S_bb+ABY*I_NAI_Gx3y_S_bb;
  Double I_NAI_Gx2yz_Py_bb = I_NAI_Hx3yz_S_bb+ABY*I_NAI_Gx2yz_S_bb;
  Double I_NAI_Gxy2z_Py_bb = I_NAI_Hx2y2z_S_bb+ABY*I_NAI_Gxy2z_S_bb;
  Double I_NAI_Gx3z_Py_bb = I_NAI_Hxy3z_S_bb+ABY*I_NAI_Gx3z_S_bb;
  Double I_NAI_G4y_Py_bb = I_NAI_H5y_S_bb+ABY*I_NAI_G4y_S_bb;
  Double I_NAI_G3yz_Py_bb = I_NAI_H4yz_S_bb+ABY*I_NAI_G3yz_S_bb;
  Double I_NAI_G2y2z_Py_bb = I_NAI_H3y2z_S_bb+ABY*I_NAI_G2y2z_S_bb;
  Double I_NAI_Gy3z_Py_bb = I_NAI_H2y3z_S_bb+ABY*I_NAI_Gy3z_S_bb;
  Double I_NAI_G4z_Py_bb = I_NAI_Hy4z_S_bb+ABY*I_NAI_G4z_S_bb;
  Double I_NAI_G4x_Pz_bb = I_NAI_H4xz_S_bb+ABZ*I_NAI_G4x_S_bb;
  Double I_NAI_G3xy_Pz_bb = I_NAI_H3xyz_S_bb+ABZ*I_NAI_G3xy_S_bb;
  Double I_NAI_G3xz_Pz_bb = I_NAI_H3x2z_S_bb+ABZ*I_NAI_G3xz_S_bb;
  Double I_NAI_G2x2y_Pz_bb = I_NAI_H2x2yz_S_bb+ABZ*I_NAI_G2x2y_S_bb;
  Double I_NAI_G2xyz_Pz_bb = I_NAI_H2xy2z_S_bb+ABZ*I_NAI_G2xyz_S_bb;
  Double I_NAI_G2x2z_Pz_bb = I_NAI_H2x3z_S_bb+ABZ*I_NAI_G2x2z_S_bb;
  Double I_NAI_Gx3y_Pz_bb = I_NAI_Hx3yz_S_bb+ABZ*I_NAI_Gx3y_S_bb;
  Double I_NAI_Gx2yz_Pz_bb = I_NAI_Hx2y2z_S_bb+ABZ*I_NAI_Gx2yz_S_bb;
  Double I_NAI_Gxy2z_Pz_bb = I_NAI_Hxy3z_S_bb+ABZ*I_NAI_Gxy2z_S_bb;
  Double I_NAI_Gx3z_Pz_bb = I_NAI_Hx4z_S_bb+ABZ*I_NAI_Gx3z_S_bb;
  Double I_NAI_G4y_Pz_bb = I_NAI_H4yz_S_bb+ABZ*I_NAI_G4y_S_bb;
  Double I_NAI_G3yz_Pz_bb = I_NAI_H3y2z_S_bb+ABZ*I_NAI_G3yz_S_bb;
  Double I_NAI_G2y2z_Pz_bb = I_NAI_H2y3z_S_bb+ABZ*I_NAI_G2y2z_S_bb;
  Double I_NAI_Gy3z_Pz_bb = I_NAI_Hy4z_S_bb+ABZ*I_NAI_Gy3z_S_bb;
  Double I_NAI_G4z_Pz_bb = I_NAI_H5z_S_bb+ABZ*I_NAI_G4z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_bb
   * RHS shell quartet name: SQ_NAI_H_S_bb
   ************************************************************/
  Double I_NAI_H5x_Px_bb = I_NAI_I6x_S_bb+ABX*I_NAI_H5x_S_bb;
  Double I_NAI_H4xy_Px_bb = I_NAI_I5xy_S_bb+ABX*I_NAI_H4xy_S_bb;
  Double I_NAI_H4xz_Px_bb = I_NAI_I5xz_S_bb+ABX*I_NAI_H4xz_S_bb;
  Double I_NAI_H3x2y_Px_bb = I_NAI_I4x2y_S_bb+ABX*I_NAI_H3x2y_S_bb;
  Double I_NAI_H3xyz_Px_bb = I_NAI_I4xyz_S_bb+ABX*I_NAI_H3xyz_S_bb;
  Double I_NAI_H3x2z_Px_bb = I_NAI_I4x2z_S_bb+ABX*I_NAI_H3x2z_S_bb;
  Double I_NAI_H2x3y_Px_bb = I_NAI_I3x3y_S_bb+ABX*I_NAI_H2x3y_S_bb;
  Double I_NAI_H2x2yz_Px_bb = I_NAI_I3x2yz_S_bb+ABX*I_NAI_H2x2yz_S_bb;
  Double I_NAI_H2xy2z_Px_bb = I_NAI_I3xy2z_S_bb+ABX*I_NAI_H2xy2z_S_bb;
  Double I_NAI_H2x3z_Px_bb = I_NAI_I3x3z_S_bb+ABX*I_NAI_H2x3z_S_bb;
  Double I_NAI_Hx4y_Px_bb = I_NAI_I2x4y_S_bb+ABX*I_NAI_Hx4y_S_bb;
  Double I_NAI_Hx3yz_Px_bb = I_NAI_I2x3yz_S_bb+ABX*I_NAI_Hx3yz_S_bb;
  Double I_NAI_Hx2y2z_Px_bb = I_NAI_I2x2y2z_S_bb+ABX*I_NAI_Hx2y2z_S_bb;
  Double I_NAI_Hxy3z_Px_bb = I_NAI_I2xy3z_S_bb+ABX*I_NAI_Hxy3z_S_bb;
  Double I_NAI_Hx4z_Px_bb = I_NAI_I2x4z_S_bb+ABX*I_NAI_Hx4z_S_bb;
  Double I_NAI_H5y_Px_bb = I_NAI_Ix5y_S_bb+ABX*I_NAI_H5y_S_bb;
  Double I_NAI_H4yz_Px_bb = I_NAI_Ix4yz_S_bb+ABX*I_NAI_H4yz_S_bb;
  Double I_NAI_H3y2z_Px_bb = I_NAI_Ix3y2z_S_bb+ABX*I_NAI_H3y2z_S_bb;
  Double I_NAI_H2y3z_Px_bb = I_NAI_Ix2y3z_S_bb+ABX*I_NAI_H2y3z_S_bb;
  Double I_NAI_Hy4z_Px_bb = I_NAI_Ixy4z_S_bb+ABX*I_NAI_Hy4z_S_bb;
  Double I_NAI_H5z_Px_bb = I_NAI_Ix5z_S_bb+ABX*I_NAI_H5z_S_bb;
  Double I_NAI_H5x_Py_bb = I_NAI_I5xy_S_bb+ABY*I_NAI_H5x_S_bb;
  Double I_NAI_H4xy_Py_bb = I_NAI_I4x2y_S_bb+ABY*I_NAI_H4xy_S_bb;
  Double I_NAI_H4xz_Py_bb = I_NAI_I4xyz_S_bb+ABY*I_NAI_H4xz_S_bb;
  Double I_NAI_H3x2y_Py_bb = I_NAI_I3x3y_S_bb+ABY*I_NAI_H3x2y_S_bb;
  Double I_NAI_H3xyz_Py_bb = I_NAI_I3x2yz_S_bb+ABY*I_NAI_H3xyz_S_bb;
  Double I_NAI_H3x2z_Py_bb = I_NAI_I3xy2z_S_bb+ABY*I_NAI_H3x2z_S_bb;
  Double I_NAI_H2x3y_Py_bb = I_NAI_I2x4y_S_bb+ABY*I_NAI_H2x3y_S_bb;
  Double I_NAI_H2x2yz_Py_bb = I_NAI_I2x3yz_S_bb+ABY*I_NAI_H2x2yz_S_bb;
  Double I_NAI_H2xy2z_Py_bb = I_NAI_I2x2y2z_S_bb+ABY*I_NAI_H2xy2z_S_bb;
  Double I_NAI_H2x3z_Py_bb = I_NAI_I2xy3z_S_bb+ABY*I_NAI_H2x3z_S_bb;
  Double I_NAI_Hx4y_Py_bb = I_NAI_Ix5y_S_bb+ABY*I_NAI_Hx4y_S_bb;
  Double I_NAI_Hx3yz_Py_bb = I_NAI_Ix4yz_S_bb+ABY*I_NAI_Hx3yz_S_bb;
  Double I_NAI_Hx2y2z_Py_bb = I_NAI_Ix3y2z_S_bb+ABY*I_NAI_Hx2y2z_S_bb;
  Double I_NAI_Hxy3z_Py_bb = I_NAI_Ix2y3z_S_bb+ABY*I_NAI_Hxy3z_S_bb;
  Double I_NAI_Hx4z_Py_bb = I_NAI_Ixy4z_S_bb+ABY*I_NAI_Hx4z_S_bb;
  Double I_NAI_H5y_Py_bb = I_NAI_I6y_S_bb+ABY*I_NAI_H5y_S_bb;
  Double I_NAI_H4yz_Py_bb = I_NAI_I5yz_S_bb+ABY*I_NAI_H4yz_S_bb;
  Double I_NAI_H3y2z_Py_bb = I_NAI_I4y2z_S_bb+ABY*I_NAI_H3y2z_S_bb;
  Double I_NAI_H2y3z_Py_bb = I_NAI_I3y3z_S_bb+ABY*I_NAI_H2y3z_S_bb;
  Double I_NAI_Hy4z_Py_bb = I_NAI_I2y4z_S_bb+ABY*I_NAI_Hy4z_S_bb;
  Double I_NAI_H5z_Py_bb = I_NAI_Iy5z_S_bb+ABY*I_NAI_H5z_S_bb;
  Double I_NAI_H5x_Pz_bb = I_NAI_I5xz_S_bb+ABZ*I_NAI_H5x_S_bb;
  Double I_NAI_H4xy_Pz_bb = I_NAI_I4xyz_S_bb+ABZ*I_NAI_H4xy_S_bb;
  Double I_NAI_H4xz_Pz_bb = I_NAI_I4x2z_S_bb+ABZ*I_NAI_H4xz_S_bb;
  Double I_NAI_H3x2y_Pz_bb = I_NAI_I3x2yz_S_bb+ABZ*I_NAI_H3x2y_S_bb;
  Double I_NAI_H3xyz_Pz_bb = I_NAI_I3xy2z_S_bb+ABZ*I_NAI_H3xyz_S_bb;
  Double I_NAI_H3x2z_Pz_bb = I_NAI_I3x3z_S_bb+ABZ*I_NAI_H3x2z_S_bb;
  Double I_NAI_H2x3y_Pz_bb = I_NAI_I2x3yz_S_bb+ABZ*I_NAI_H2x3y_S_bb;
  Double I_NAI_H2x2yz_Pz_bb = I_NAI_I2x2y2z_S_bb+ABZ*I_NAI_H2x2yz_S_bb;
  Double I_NAI_H2xy2z_Pz_bb = I_NAI_I2xy3z_S_bb+ABZ*I_NAI_H2xy2z_S_bb;
  Double I_NAI_H2x3z_Pz_bb = I_NAI_I2x4z_S_bb+ABZ*I_NAI_H2x3z_S_bb;
  Double I_NAI_Hx4y_Pz_bb = I_NAI_Ix4yz_S_bb+ABZ*I_NAI_Hx4y_S_bb;
  Double I_NAI_Hx3yz_Pz_bb = I_NAI_Ix3y2z_S_bb+ABZ*I_NAI_Hx3yz_S_bb;
  Double I_NAI_Hx2y2z_Pz_bb = I_NAI_Ix2y3z_S_bb+ABZ*I_NAI_Hx2y2z_S_bb;
  Double I_NAI_Hxy3z_Pz_bb = I_NAI_Ixy4z_S_bb+ABZ*I_NAI_Hxy3z_S_bb;
  Double I_NAI_Hx4z_Pz_bb = I_NAI_Ix5z_S_bb+ABZ*I_NAI_Hx4z_S_bb;
  Double I_NAI_H5y_Pz_bb = I_NAI_I5yz_S_bb+ABZ*I_NAI_H5y_S_bb;
  Double I_NAI_H4yz_Pz_bb = I_NAI_I4y2z_S_bb+ABZ*I_NAI_H4yz_S_bb;
  Double I_NAI_H3y2z_Pz_bb = I_NAI_I3y3z_S_bb+ABZ*I_NAI_H3y2z_S_bb;
  Double I_NAI_H2y3z_Pz_bb = I_NAI_I2y4z_S_bb+ABZ*I_NAI_H2y3z_S_bb;
  Double I_NAI_Hy4z_Pz_bb = I_NAI_Iy5z_S_bb+ABZ*I_NAI_Hy4z_S_bb;
  Double I_NAI_H5z_Pz_bb = I_NAI_I6z_S_bb+ABZ*I_NAI_H5z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_bb
   * RHS shell quartet name: SQ_NAI_G_P_bb
   ************************************************************/
  Double I_NAI_G4x_D2x_bb = I_NAI_H5x_Px_bb+ABX*I_NAI_G4x_Px_bb;
  Double I_NAI_G3xy_D2x_bb = I_NAI_H4xy_Px_bb+ABX*I_NAI_G3xy_Px_bb;
  Double I_NAI_G3xz_D2x_bb = I_NAI_H4xz_Px_bb+ABX*I_NAI_G3xz_Px_bb;
  Double I_NAI_G2x2y_D2x_bb = I_NAI_H3x2y_Px_bb+ABX*I_NAI_G2x2y_Px_bb;
  Double I_NAI_G2xyz_D2x_bb = I_NAI_H3xyz_Px_bb+ABX*I_NAI_G2xyz_Px_bb;
  Double I_NAI_G2x2z_D2x_bb = I_NAI_H3x2z_Px_bb+ABX*I_NAI_G2x2z_Px_bb;
  Double I_NAI_Gx3y_D2x_bb = I_NAI_H2x3y_Px_bb+ABX*I_NAI_Gx3y_Px_bb;
  Double I_NAI_Gx2yz_D2x_bb = I_NAI_H2x2yz_Px_bb+ABX*I_NAI_Gx2yz_Px_bb;
  Double I_NAI_Gxy2z_D2x_bb = I_NAI_H2xy2z_Px_bb+ABX*I_NAI_Gxy2z_Px_bb;
  Double I_NAI_Gx3z_D2x_bb = I_NAI_H2x3z_Px_bb+ABX*I_NAI_Gx3z_Px_bb;
  Double I_NAI_G4y_D2x_bb = I_NAI_Hx4y_Px_bb+ABX*I_NAI_G4y_Px_bb;
  Double I_NAI_G3yz_D2x_bb = I_NAI_Hx3yz_Px_bb+ABX*I_NAI_G3yz_Px_bb;
  Double I_NAI_G2y2z_D2x_bb = I_NAI_Hx2y2z_Px_bb+ABX*I_NAI_G2y2z_Px_bb;
  Double I_NAI_Gy3z_D2x_bb = I_NAI_Hxy3z_Px_bb+ABX*I_NAI_Gy3z_Px_bb;
  Double I_NAI_G4z_D2x_bb = I_NAI_Hx4z_Px_bb+ABX*I_NAI_G4z_Px_bb;
  Double I_NAI_G4x_D2y_bb = I_NAI_H4xy_Py_bb+ABY*I_NAI_G4x_Py_bb;
  Double I_NAI_G3xy_D2y_bb = I_NAI_H3x2y_Py_bb+ABY*I_NAI_G3xy_Py_bb;
  Double I_NAI_G3xz_D2y_bb = I_NAI_H3xyz_Py_bb+ABY*I_NAI_G3xz_Py_bb;
  Double I_NAI_G2x2y_D2y_bb = I_NAI_H2x3y_Py_bb+ABY*I_NAI_G2x2y_Py_bb;
  Double I_NAI_G2xyz_D2y_bb = I_NAI_H2x2yz_Py_bb+ABY*I_NAI_G2xyz_Py_bb;
  Double I_NAI_G2x2z_D2y_bb = I_NAI_H2xy2z_Py_bb+ABY*I_NAI_G2x2z_Py_bb;
  Double I_NAI_Gx3y_D2y_bb = I_NAI_Hx4y_Py_bb+ABY*I_NAI_Gx3y_Py_bb;
  Double I_NAI_Gx2yz_D2y_bb = I_NAI_Hx3yz_Py_bb+ABY*I_NAI_Gx2yz_Py_bb;
  Double I_NAI_Gxy2z_D2y_bb = I_NAI_Hx2y2z_Py_bb+ABY*I_NAI_Gxy2z_Py_bb;
  Double I_NAI_Gx3z_D2y_bb = I_NAI_Hxy3z_Py_bb+ABY*I_NAI_Gx3z_Py_bb;
  Double I_NAI_G4y_D2y_bb = I_NAI_H5y_Py_bb+ABY*I_NAI_G4y_Py_bb;
  Double I_NAI_G3yz_D2y_bb = I_NAI_H4yz_Py_bb+ABY*I_NAI_G3yz_Py_bb;
  Double I_NAI_G2y2z_D2y_bb = I_NAI_H3y2z_Py_bb+ABY*I_NAI_G2y2z_Py_bb;
  Double I_NAI_Gy3z_D2y_bb = I_NAI_H2y3z_Py_bb+ABY*I_NAI_Gy3z_Py_bb;
  Double I_NAI_G4z_D2y_bb = I_NAI_Hy4z_Py_bb+ABY*I_NAI_G4z_Py_bb;
  Double I_NAI_G4x_D2z_bb = I_NAI_H4xz_Pz_bb+ABZ*I_NAI_G4x_Pz_bb;
  Double I_NAI_G3xy_D2z_bb = I_NAI_H3xyz_Pz_bb+ABZ*I_NAI_G3xy_Pz_bb;
  Double I_NAI_G3xz_D2z_bb = I_NAI_H3x2z_Pz_bb+ABZ*I_NAI_G3xz_Pz_bb;
  Double I_NAI_G2x2y_D2z_bb = I_NAI_H2x2yz_Pz_bb+ABZ*I_NAI_G2x2y_Pz_bb;
  Double I_NAI_G2xyz_D2z_bb = I_NAI_H2xy2z_Pz_bb+ABZ*I_NAI_G2xyz_Pz_bb;
  Double I_NAI_G2x2z_D2z_bb = I_NAI_H2x3z_Pz_bb+ABZ*I_NAI_G2x2z_Pz_bb;
  Double I_NAI_Gx3y_D2z_bb = I_NAI_Hx3yz_Pz_bb+ABZ*I_NAI_Gx3y_Pz_bb;
  Double I_NAI_Gx2yz_D2z_bb = I_NAI_Hx2y2z_Pz_bb+ABZ*I_NAI_Gx2yz_Pz_bb;
  Double I_NAI_Gxy2z_D2z_bb = I_NAI_Hxy3z_Pz_bb+ABZ*I_NAI_Gxy2z_Pz_bb;
  Double I_NAI_Gx3z_D2z_bb = I_NAI_Hx4z_Pz_bb+ABZ*I_NAI_Gx3z_Pz_bb;
  Double I_NAI_G4y_D2z_bb = I_NAI_H4yz_Pz_bb+ABZ*I_NAI_G4y_Pz_bb;
  Double I_NAI_G3yz_D2z_bb = I_NAI_H3y2z_Pz_bb+ABZ*I_NAI_G3yz_Pz_bb;
  Double I_NAI_G2y2z_D2z_bb = I_NAI_H2y3z_Pz_bb+ABZ*I_NAI_G2y2z_Pz_bb;
  Double I_NAI_Gy3z_D2z_bb = I_NAI_Hy4z_Pz_bb+ABZ*I_NAI_Gy3z_Pz_bb;
  Double I_NAI_G4z_D2z_bb = I_NAI_H5z_Pz_bb+ABZ*I_NAI_G4z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_bb
   * RHS shell quartet name: SQ_NAI_I_S_bb
   ************************************************************/
  Double I_NAI_I6x_Px_bb = I_NAI_K7x_S_bb+ABX*I_NAI_I6x_S_bb;
  Double I_NAI_I5xy_Px_bb = I_NAI_K6xy_S_bb+ABX*I_NAI_I5xy_S_bb;
  Double I_NAI_I5xz_Px_bb = I_NAI_K6xz_S_bb+ABX*I_NAI_I5xz_S_bb;
  Double I_NAI_I4x2y_Px_bb = I_NAI_K5x2y_S_bb+ABX*I_NAI_I4x2y_S_bb;
  Double I_NAI_I4xyz_Px_bb = I_NAI_K5xyz_S_bb+ABX*I_NAI_I4xyz_S_bb;
  Double I_NAI_I4x2z_Px_bb = I_NAI_K5x2z_S_bb+ABX*I_NAI_I4x2z_S_bb;
  Double I_NAI_I3x3y_Px_bb = I_NAI_K4x3y_S_bb+ABX*I_NAI_I3x3y_S_bb;
  Double I_NAI_I3x2yz_Px_bb = I_NAI_K4x2yz_S_bb+ABX*I_NAI_I3x2yz_S_bb;
  Double I_NAI_I3xy2z_Px_bb = I_NAI_K4xy2z_S_bb+ABX*I_NAI_I3xy2z_S_bb;
  Double I_NAI_I3x3z_Px_bb = I_NAI_K4x3z_S_bb+ABX*I_NAI_I3x3z_S_bb;
  Double I_NAI_I2x4y_Px_bb = I_NAI_K3x4y_S_bb+ABX*I_NAI_I2x4y_S_bb;
  Double I_NAI_I2x3yz_Px_bb = I_NAI_K3x3yz_S_bb+ABX*I_NAI_I2x3yz_S_bb;
  Double I_NAI_I2x2y2z_Px_bb = I_NAI_K3x2y2z_S_bb+ABX*I_NAI_I2x2y2z_S_bb;
  Double I_NAI_I2xy3z_Px_bb = I_NAI_K3xy3z_S_bb+ABX*I_NAI_I2xy3z_S_bb;
  Double I_NAI_I2x4z_Px_bb = I_NAI_K3x4z_S_bb+ABX*I_NAI_I2x4z_S_bb;
  Double I_NAI_Ix5y_Px_bb = I_NAI_K2x5y_S_bb+ABX*I_NAI_Ix5y_S_bb;
  Double I_NAI_Ix4yz_Px_bb = I_NAI_K2x4yz_S_bb+ABX*I_NAI_Ix4yz_S_bb;
  Double I_NAI_Ix3y2z_Px_bb = I_NAI_K2x3y2z_S_bb+ABX*I_NAI_Ix3y2z_S_bb;
  Double I_NAI_Ix2y3z_Px_bb = I_NAI_K2x2y3z_S_bb+ABX*I_NAI_Ix2y3z_S_bb;
  Double I_NAI_Ixy4z_Px_bb = I_NAI_K2xy4z_S_bb+ABX*I_NAI_Ixy4z_S_bb;
  Double I_NAI_Ix5z_Px_bb = I_NAI_K2x5z_S_bb+ABX*I_NAI_Ix5z_S_bb;
  Double I_NAI_I6y_Px_bb = I_NAI_Kx6y_S_bb+ABX*I_NAI_I6y_S_bb;
  Double I_NAI_I5yz_Px_bb = I_NAI_Kx5yz_S_bb+ABX*I_NAI_I5yz_S_bb;
  Double I_NAI_I4y2z_Px_bb = I_NAI_Kx4y2z_S_bb+ABX*I_NAI_I4y2z_S_bb;
  Double I_NAI_I3y3z_Px_bb = I_NAI_Kx3y3z_S_bb+ABX*I_NAI_I3y3z_S_bb;
  Double I_NAI_I2y4z_Px_bb = I_NAI_Kx2y4z_S_bb+ABX*I_NAI_I2y4z_S_bb;
  Double I_NAI_Iy5z_Px_bb = I_NAI_Kxy5z_S_bb+ABX*I_NAI_Iy5z_S_bb;
  Double I_NAI_I6z_Px_bb = I_NAI_Kx6z_S_bb+ABX*I_NAI_I6z_S_bb;
  Double I_NAI_I6x_Py_bb = I_NAI_K6xy_S_bb+ABY*I_NAI_I6x_S_bb;
  Double I_NAI_I5xy_Py_bb = I_NAI_K5x2y_S_bb+ABY*I_NAI_I5xy_S_bb;
  Double I_NAI_I5xz_Py_bb = I_NAI_K5xyz_S_bb+ABY*I_NAI_I5xz_S_bb;
  Double I_NAI_I4x2y_Py_bb = I_NAI_K4x3y_S_bb+ABY*I_NAI_I4x2y_S_bb;
  Double I_NAI_I4xyz_Py_bb = I_NAI_K4x2yz_S_bb+ABY*I_NAI_I4xyz_S_bb;
  Double I_NAI_I4x2z_Py_bb = I_NAI_K4xy2z_S_bb+ABY*I_NAI_I4x2z_S_bb;
  Double I_NAI_I3x3y_Py_bb = I_NAI_K3x4y_S_bb+ABY*I_NAI_I3x3y_S_bb;
  Double I_NAI_I3x2yz_Py_bb = I_NAI_K3x3yz_S_bb+ABY*I_NAI_I3x2yz_S_bb;
  Double I_NAI_I3xy2z_Py_bb = I_NAI_K3x2y2z_S_bb+ABY*I_NAI_I3xy2z_S_bb;
  Double I_NAI_I3x3z_Py_bb = I_NAI_K3xy3z_S_bb+ABY*I_NAI_I3x3z_S_bb;
  Double I_NAI_I2x4y_Py_bb = I_NAI_K2x5y_S_bb+ABY*I_NAI_I2x4y_S_bb;
  Double I_NAI_I2x3yz_Py_bb = I_NAI_K2x4yz_S_bb+ABY*I_NAI_I2x3yz_S_bb;
  Double I_NAI_I2x2y2z_Py_bb = I_NAI_K2x3y2z_S_bb+ABY*I_NAI_I2x2y2z_S_bb;
  Double I_NAI_I2xy3z_Py_bb = I_NAI_K2x2y3z_S_bb+ABY*I_NAI_I2xy3z_S_bb;
  Double I_NAI_I2x4z_Py_bb = I_NAI_K2xy4z_S_bb+ABY*I_NAI_I2x4z_S_bb;
  Double I_NAI_Ix5y_Py_bb = I_NAI_Kx6y_S_bb+ABY*I_NAI_Ix5y_S_bb;
  Double I_NAI_Ix4yz_Py_bb = I_NAI_Kx5yz_S_bb+ABY*I_NAI_Ix4yz_S_bb;
  Double I_NAI_Ix3y2z_Py_bb = I_NAI_Kx4y2z_S_bb+ABY*I_NAI_Ix3y2z_S_bb;
  Double I_NAI_Ix2y3z_Py_bb = I_NAI_Kx3y3z_S_bb+ABY*I_NAI_Ix2y3z_S_bb;
  Double I_NAI_Ixy4z_Py_bb = I_NAI_Kx2y4z_S_bb+ABY*I_NAI_Ixy4z_S_bb;
  Double I_NAI_Ix5z_Py_bb = I_NAI_Kxy5z_S_bb+ABY*I_NAI_Ix5z_S_bb;
  Double I_NAI_I6y_Py_bb = I_NAI_K7y_S_bb+ABY*I_NAI_I6y_S_bb;
  Double I_NAI_I5yz_Py_bb = I_NAI_K6yz_S_bb+ABY*I_NAI_I5yz_S_bb;
  Double I_NAI_I4y2z_Py_bb = I_NAI_K5y2z_S_bb+ABY*I_NAI_I4y2z_S_bb;
  Double I_NAI_I3y3z_Py_bb = I_NAI_K4y3z_S_bb+ABY*I_NAI_I3y3z_S_bb;
  Double I_NAI_I2y4z_Py_bb = I_NAI_K3y4z_S_bb+ABY*I_NAI_I2y4z_S_bb;
  Double I_NAI_Iy5z_Py_bb = I_NAI_K2y5z_S_bb+ABY*I_NAI_Iy5z_S_bb;
  Double I_NAI_I6z_Py_bb = I_NAI_Ky6z_S_bb+ABY*I_NAI_I6z_S_bb;
  Double I_NAI_I6x_Pz_bb = I_NAI_K6xz_S_bb+ABZ*I_NAI_I6x_S_bb;
  Double I_NAI_I5xy_Pz_bb = I_NAI_K5xyz_S_bb+ABZ*I_NAI_I5xy_S_bb;
  Double I_NAI_I5xz_Pz_bb = I_NAI_K5x2z_S_bb+ABZ*I_NAI_I5xz_S_bb;
  Double I_NAI_I4x2y_Pz_bb = I_NAI_K4x2yz_S_bb+ABZ*I_NAI_I4x2y_S_bb;
  Double I_NAI_I4xyz_Pz_bb = I_NAI_K4xy2z_S_bb+ABZ*I_NAI_I4xyz_S_bb;
  Double I_NAI_I4x2z_Pz_bb = I_NAI_K4x3z_S_bb+ABZ*I_NAI_I4x2z_S_bb;
  Double I_NAI_I3x3y_Pz_bb = I_NAI_K3x3yz_S_bb+ABZ*I_NAI_I3x3y_S_bb;
  Double I_NAI_I3x2yz_Pz_bb = I_NAI_K3x2y2z_S_bb+ABZ*I_NAI_I3x2yz_S_bb;
  Double I_NAI_I3xy2z_Pz_bb = I_NAI_K3xy3z_S_bb+ABZ*I_NAI_I3xy2z_S_bb;
  Double I_NAI_I3x3z_Pz_bb = I_NAI_K3x4z_S_bb+ABZ*I_NAI_I3x3z_S_bb;
  Double I_NAI_I2x4y_Pz_bb = I_NAI_K2x4yz_S_bb+ABZ*I_NAI_I2x4y_S_bb;
  Double I_NAI_I2x3yz_Pz_bb = I_NAI_K2x3y2z_S_bb+ABZ*I_NAI_I2x3yz_S_bb;
  Double I_NAI_I2x2y2z_Pz_bb = I_NAI_K2x2y3z_S_bb+ABZ*I_NAI_I2x2y2z_S_bb;
  Double I_NAI_I2xy3z_Pz_bb = I_NAI_K2xy4z_S_bb+ABZ*I_NAI_I2xy3z_S_bb;
  Double I_NAI_I2x4z_Pz_bb = I_NAI_K2x5z_S_bb+ABZ*I_NAI_I2x4z_S_bb;
  Double I_NAI_Ix5y_Pz_bb = I_NAI_Kx5yz_S_bb+ABZ*I_NAI_Ix5y_S_bb;
  Double I_NAI_Ix4yz_Pz_bb = I_NAI_Kx4y2z_S_bb+ABZ*I_NAI_Ix4yz_S_bb;
  Double I_NAI_Ix3y2z_Pz_bb = I_NAI_Kx3y3z_S_bb+ABZ*I_NAI_Ix3y2z_S_bb;
  Double I_NAI_Ix2y3z_Pz_bb = I_NAI_Kx2y4z_S_bb+ABZ*I_NAI_Ix2y3z_S_bb;
  Double I_NAI_Ixy4z_Pz_bb = I_NAI_Kxy5z_S_bb+ABZ*I_NAI_Ixy4z_S_bb;
  Double I_NAI_Ix5z_Pz_bb = I_NAI_Kx6z_S_bb+ABZ*I_NAI_Ix5z_S_bb;
  Double I_NAI_I6y_Pz_bb = I_NAI_K6yz_S_bb+ABZ*I_NAI_I6y_S_bb;
  Double I_NAI_I5yz_Pz_bb = I_NAI_K5y2z_S_bb+ABZ*I_NAI_I5yz_S_bb;
  Double I_NAI_I4y2z_Pz_bb = I_NAI_K4y3z_S_bb+ABZ*I_NAI_I4y2z_S_bb;
  Double I_NAI_I3y3z_Pz_bb = I_NAI_K3y4z_S_bb+ABZ*I_NAI_I3y3z_S_bb;
  Double I_NAI_I2y4z_Pz_bb = I_NAI_K2y5z_S_bb+ABZ*I_NAI_I2y4z_S_bb;
  Double I_NAI_Iy5z_Pz_bb = I_NAI_Ky6z_S_bb+ABZ*I_NAI_Iy5z_S_bb;
  Double I_NAI_I6z_Pz_bb = I_NAI_K7z_S_bb+ABZ*I_NAI_I6z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_bb
   * RHS shell quartet name: SQ_NAI_H_P_bb
   ************************************************************/
  Double I_NAI_H5x_D2x_bb = I_NAI_I6x_Px_bb+ABX*I_NAI_H5x_Px_bb;
  Double I_NAI_H4xy_D2x_bb = I_NAI_I5xy_Px_bb+ABX*I_NAI_H4xy_Px_bb;
  Double I_NAI_H4xz_D2x_bb = I_NAI_I5xz_Px_bb+ABX*I_NAI_H4xz_Px_bb;
  Double I_NAI_H3x2y_D2x_bb = I_NAI_I4x2y_Px_bb+ABX*I_NAI_H3x2y_Px_bb;
  Double I_NAI_H3xyz_D2x_bb = I_NAI_I4xyz_Px_bb+ABX*I_NAI_H3xyz_Px_bb;
  Double I_NAI_H3x2z_D2x_bb = I_NAI_I4x2z_Px_bb+ABX*I_NAI_H3x2z_Px_bb;
  Double I_NAI_H2x3y_D2x_bb = I_NAI_I3x3y_Px_bb+ABX*I_NAI_H2x3y_Px_bb;
  Double I_NAI_H2x2yz_D2x_bb = I_NAI_I3x2yz_Px_bb+ABX*I_NAI_H2x2yz_Px_bb;
  Double I_NAI_H2xy2z_D2x_bb = I_NAI_I3xy2z_Px_bb+ABX*I_NAI_H2xy2z_Px_bb;
  Double I_NAI_H2x3z_D2x_bb = I_NAI_I3x3z_Px_bb+ABX*I_NAI_H2x3z_Px_bb;
  Double I_NAI_Hx4y_D2x_bb = I_NAI_I2x4y_Px_bb+ABX*I_NAI_Hx4y_Px_bb;
  Double I_NAI_Hx3yz_D2x_bb = I_NAI_I2x3yz_Px_bb+ABX*I_NAI_Hx3yz_Px_bb;
  Double I_NAI_Hx2y2z_D2x_bb = I_NAI_I2x2y2z_Px_bb+ABX*I_NAI_Hx2y2z_Px_bb;
  Double I_NAI_Hxy3z_D2x_bb = I_NAI_I2xy3z_Px_bb+ABX*I_NAI_Hxy3z_Px_bb;
  Double I_NAI_Hx4z_D2x_bb = I_NAI_I2x4z_Px_bb+ABX*I_NAI_Hx4z_Px_bb;
  Double I_NAI_H5y_D2x_bb = I_NAI_Ix5y_Px_bb+ABX*I_NAI_H5y_Px_bb;
  Double I_NAI_H4yz_D2x_bb = I_NAI_Ix4yz_Px_bb+ABX*I_NAI_H4yz_Px_bb;
  Double I_NAI_H3y2z_D2x_bb = I_NAI_Ix3y2z_Px_bb+ABX*I_NAI_H3y2z_Px_bb;
  Double I_NAI_H2y3z_D2x_bb = I_NAI_Ix2y3z_Px_bb+ABX*I_NAI_H2y3z_Px_bb;
  Double I_NAI_Hy4z_D2x_bb = I_NAI_Ixy4z_Px_bb+ABX*I_NAI_Hy4z_Px_bb;
  Double I_NAI_H5z_D2x_bb = I_NAI_Ix5z_Px_bb+ABX*I_NAI_H5z_Px_bb;
  Double I_NAI_H5x_D2y_bb = I_NAI_I5xy_Py_bb+ABY*I_NAI_H5x_Py_bb;
  Double I_NAI_H4xy_D2y_bb = I_NAI_I4x2y_Py_bb+ABY*I_NAI_H4xy_Py_bb;
  Double I_NAI_H4xz_D2y_bb = I_NAI_I4xyz_Py_bb+ABY*I_NAI_H4xz_Py_bb;
  Double I_NAI_H3x2y_D2y_bb = I_NAI_I3x3y_Py_bb+ABY*I_NAI_H3x2y_Py_bb;
  Double I_NAI_H3xyz_D2y_bb = I_NAI_I3x2yz_Py_bb+ABY*I_NAI_H3xyz_Py_bb;
  Double I_NAI_H3x2z_D2y_bb = I_NAI_I3xy2z_Py_bb+ABY*I_NAI_H3x2z_Py_bb;
  Double I_NAI_H2x3y_D2y_bb = I_NAI_I2x4y_Py_bb+ABY*I_NAI_H2x3y_Py_bb;
  Double I_NAI_H2x2yz_D2y_bb = I_NAI_I2x3yz_Py_bb+ABY*I_NAI_H2x2yz_Py_bb;
  Double I_NAI_H2xy2z_D2y_bb = I_NAI_I2x2y2z_Py_bb+ABY*I_NAI_H2xy2z_Py_bb;
  Double I_NAI_H2x3z_D2y_bb = I_NAI_I2xy3z_Py_bb+ABY*I_NAI_H2x3z_Py_bb;
  Double I_NAI_Hx4y_D2y_bb = I_NAI_Ix5y_Py_bb+ABY*I_NAI_Hx4y_Py_bb;
  Double I_NAI_Hx3yz_D2y_bb = I_NAI_Ix4yz_Py_bb+ABY*I_NAI_Hx3yz_Py_bb;
  Double I_NAI_Hx2y2z_D2y_bb = I_NAI_Ix3y2z_Py_bb+ABY*I_NAI_Hx2y2z_Py_bb;
  Double I_NAI_Hxy3z_D2y_bb = I_NAI_Ix2y3z_Py_bb+ABY*I_NAI_Hxy3z_Py_bb;
  Double I_NAI_Hx4z_D2y_bb = I_NAI_Ixy4z_Py_bb+ABY*I_NAI_Hx4z_Py_bb;
  Double I_NAI_H5y_D2y_bb = I_NAI_I6y_Py_bb+ABY*I_NAI_H5y_Py_bb;
  Double I_NAI_H4yz_D2y_bb = I_NAI_I5yz_Py_bb+ABY*I_NAI_H4yz_Py_bb;
  Double I_NAI_H3y2z_D2y_bb = I_NAI_I4y2z_Py_bb+ABY*I_NAI_H3y2z_Py_bb;
  Double I_NAI_H2y3z_D2y_bb = I_NAI_I3y3z_Py_bb+ABY*I_NAI_H2y3z_Py_bb;
  Double I_NAI_Hy4z_D2y_bb = I_NAI_I2y4z_Py_bb+ABY*I_NAI_Hy4z_Py_bb;
  Double I_NAI_H5z_D2y_bb = I_NAI_Iy5z_Py_bb+ABY*I_NAI_H5z_Py_bb;
  Double I_NAI_H5x_D2z_bb = I_NAI_I5xz_Pz_bb+ABZ*I_NAI_H5x_Pz_bb;
  Double I_NAI_H4xy_D2z_bb = I_NAI_I4xyz_Pz_bb+ABZ*I_NAI_H4xy_Pz_bb;
  Double I_NAI_H4xz_D2z_bb = I_NAI_I4x2z_Pz_bb+ABZ*I_NAI_H4xz_Pz_bb;
  Double I_NAI_H3x2y_D2z_bb = I_NAI_I3x2yz_Pz_bb+ABZ*I_NAI_H3x2y_Pz_bb;
  Double I_NAI_H3xyz_D2z_bb = I_NAI_I3xy2z_Pz_bb+ABZ*I_NAI_H3xyz_Pz_bb;
  Double I_NAI_H3x2z_D2z_bb = I_NAI_I3x3z_Pz_bb+ABZ*I_NAI_H3x2z_Pz_bb;
  Double I_NAI_H2x3y_D2z_bb = I_NAI_I2x3yz_Pz_bb+ABZ*I_NAI_H2x3y_Pz_bb;
  Double I_NAI_H2x2yz_D2z_bb = I_NAI_I2x2y2z_Pz_bb+ABZ*I_NAI_H2x2yz_Pz_bb;
  Double I_NAI_H2xy2z_D2z_bb = I_NAI_I2xy3z_Pz_bb+ABZ*I_NAI_H2xy2z_Pz_bb;
  Double I_NAI_H2x3z_D2z_bb = I_NAI_I2x4z_Pz_bb+ABZ*I_NAI_H2x3z_Pz_bb;
  Double I_NAI_Hx4y_D2z_bb = I_NAI_Ix4yz_Pz_bb+ABZ*I_NAI_Hx4y_Pz_bb;
  Double I_NAI_Hx3yz_D2z_bb = I_NAI_Ix3y2z_Pz_bb+ABZ*I_NAI_Hx3yz_Pz_bb;
  Double I_NAI_Hx2y2z_D2z_bb = I_NAI_Ix2y3z_Pz_bb+ABZ*I_NAI_Hx2y2z_Pz_bb;
  Double I_NAI_Hxy3z_D2z_bb = I_NAI_Ixy4z_Pz_bb+ABZ*I_NAI_Hxy3z_Pz_bb;
  Double I_NAI_Hx4z_D2z_bb = I_NAI_Ix5z_Pz_bb+ABZ*I_NAI_Hx4z_Pz_bb;
  Double I_NAI_H5y_D2z_bb = I_NAI_I5yz_Pz_bb+ABZ*I_NAI_H5y_Pz_bb;
  Double I_NAI_H4yz_D2z_bb = I_NAI_I4y2z_Pz_bb+ABZ*I_NAI_H4yz_Pz_bb;
  Double I_NAI_H3y2z_D2z_bb = I_NAI_I3y3z_Pz_bb+ABZ*I_NAI_H3y2z_Pz_bb;
  Double I_NAI_H2y3z_D2z_bb = I_NAI_I2y4z_Pz_bb+ABZ*I_NAI_H2y3z_Pz_bb;
  Double I_NAI_Hy4z_D2z_bb = I_NAI_Iy5z_Pz_bb+ABZ*I_NAI_Hy4z_Pz_bb;
  Double I_NAI_H5z_D2z_bb = I_NAI_I6z_Pz_bb+ABZ*I_NAI_H5z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 60 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_bb
   * RHS shell quartet name: SQ_NAI_G_D_bb
   ************************************************************/
  Double I_NAI_G4x_F3x_bb = I_NAI_H5x_D2x_bb+ABX*I_NAI_G4x_D2x_bb;
  Double I_NAI_G3xy_F3x_bb = I_NAI_H4xy_D2x_bb+ABX*I_NAI_G3xy_D2x_bb;
  Double I_NAI_G3xz_F3x_bb = I_NAI_H4xz_D2x_bb+ABX*I_NAI_G3xz_D2x_bb;
  Double I_NAI_G2x2y_F3x_bb = I_NAI_H3x2y_D2x_bb+ABX*I_NAI_G2x2y_D2x_bb;
  Double I_NAI_G2xyz_F3x_bb = I_NAI_H3xyz_D2x_bb+ABX*I_NAI_G2xyz_D2x_bb;
  Double I_NAI_G2x2z_F3x_bb = I_NAI_H3x2z_D2x_bb+ABX*I_NAI_G2x2z_D2x_bb;
  Double I_NAI_Gx3y_F3x_bb = I_NAI_H2x3y_D2x_bb+ABX*I_NAI_Gx3y_D2x_bb;
  Double I_NAI_Gx2yz_F3x_bb = I_NAI_H2x2yz_D2x_bb+ABX*I_NAI_Gx2yz_D2x_bb;
  Double I_NAI_Gxy2z_F3x_bb = I_NAI_H2xy2z_D2x_bb+ABX*I_NAI_Gxy2z_D2x_bb;
  Double I_NAI_Gx3z_F3x_bb = I_NAI_H2x3z_D2x_bb+ABX*I_NAI_Gx3z_D2x_bb;
  Double I_NAI_G4y_F3x_bb = I_NAI_Hx4y_D2x_bb+ABX*I_NAI_G4y_D2x_bb;
  Double I_NAI_G3yz_F3x_bb = I_NAI_Hx3yz_D2x_bb+ABX*I_NAI_G3yz_D2x_bb;
  Double I_NAI_G2y2z_F3x_bb = I_NAI_Hx2y2z_D2x_bb+ABX*I_NAI_G2y2z_D2x_bb;
  Double I_NAI_Gy3z_F3x_bb = I_NAI_Hxy3z_D2x_bb+ABX*I_NAI_Gy3z_D2x_bb;
  Double I_NAI_G4z_F3x_bb = I_NAI_Hx4z_D2x_bb+ABX*I_NAI_G4z_D2x_bb;
  Double I_NAI_G4x_F2xy_bb = I_NAI_H4xy_D2x_bb+ABY*I_NAI_G4x_D2x_bb;
  Double I_NAI_G3xy_F2xy_bb = I_NAI_H3x2y_D2x_bb+ABY*I_NAI_G3xy_D2x_bb;
  Double I_NAI_G3xz_F2xy_bb = I_NAI_H3xyz_D2x_bb+ABY*I_NAI_G3xz_D2x_bb;
  Double I_NAI_G2x2y_F2xy_bb = I_NAI_H2x3y_D2x_bb+ABY*I_NAI_G2x2y_D2x_bb;
  Double I_NAI_G2xyz_F2xy_bb = I_NAI_H2x2yz_D2x_bb+ABY*I_NAI_G2xyz_D2x_bb;
  Double I_NAI_G2x2z_F2xy_bb = I_NAI_H2xy2z_D2x_bb+ABY*I_NAI_G2x2z_D2x_bb;
  Double I_NAI_Gx3y_F2xy_bb = I_NAI_Hx4y_D2x_bb+ABY*I_NAI_Gx3y_D2x_bb;
  Double I_NAI_Gx2yz_F2xy_bb = I_NAI_Hx3yz_D2x_bb+ABY*I_NAI_Gx2yz_D2x_bb;
  Double I_NAI_Gxy2z_F2xy_bb = I_NAI_Hx2y2z_D2x_bb+ABY*I_NAI_Gxy2z_D2x_bb;
  Double I_NAI_Gx3z_F2xy_bb = I_NAI_Hxy3z_D2x_bb+ABY*I_NAI_Gx3z_D2x_bb;
  Double I_NAI_G4y_F2xy_bb = I_NAI_H5y_D2x_bb+ABY*I_NAI_G4y_D2x_bb;
  Double I_NAI_G3yz_F2xy_bb = I_NAI_H4yz_D2x_bb+ABY*I_NAI_G3yz_D2x_bb;
  Double I_NAI_G2y2z_F2xy_bb = I_NAI_H3y2z_D2x_bb+ABY*I_NAI_G2y2z_D2x_bb;
  Double I_NAI_Gy3z_F2xy_bb = I_NAI_H2y3z_D2x_bb+ABY*I_NAI_Gy3z_D2x_bb;
  Double I_NAI_G4z_F2xy_bb = I_NAI_Hy4z_D2x_bb+ABY*I_NAI_G4z_D2x_bb;
  Double I_NAI_G4x_F2xz_bb = I_NAI_H4xz_D2x_bb+ABZ*I_NAI_G4x_D2x_bb;
  Double I_NAI_G3xy_F2xz_bb = I_NAI_H3xyz_D2x_bb+ABZ*I_NAI_G3xy_D2x_bb;
  Double I_NAI_G3xz_F2xz_bb = I_NAI_H3x2z_D2x_bb+ABZ*I_NAI_G3xz_D2x_bb;
  Double I_NAI_G2x2y_F2xz_bb = I_NAI_H2x2yz_D2x_bb+ABZ*I_NAI_G2x2y_D2x_bb;
  Double I_NAI_G2xyz_F2xz_bb = I_NAI_H2xy2z_D2x_bb+ABZ*I_NAI_G2xyz_D2x_bb;
  Double I_NAI_G2x2z_F2xz_bb = I_NAI_H2x3z_D2x_bb+ABZ*I_NAI_G2x2z_D2x_bb;
  Double I_NAI_Gx3y_F2xz_bb = I_NAI_Hx3yz_D2x_bb+ABZ*I_NAI_Gx3y_D2x_bb;
  Double I_NAI_Gx2yz_F2xz_bb = I_NAI_Hx2y2z_D2x_bb+ABZ*I_NAI_Gx2yz_D2x_bb;
  Double I_NAI_Gxy2z_F2xz_bb = I_NAI_Hxy3z_D2x_bb+ABZ*I_NAI_Gxy2z_D2x_bb;
  Double I_NAI_Gx3z_F2xz_bb = I_NAI_Hx4z_D2x_bb+ABZ*I_NAI_Gx3z_D2x_bb;
  Double I_NAI_G4y_F2xz_bb = I_NAI_H4yz_D2x_bb+ABZ*I_NAI_G4y_D2x_bb;
  Double I_NAI_G3yz_F2xz_bb = I_NAI_H3y2z_D2x_bb+ABZ*I_NAI_G3yz_D2x_bb;
  Double I_NAI_G2y2z_F2xz_bb = I_NAI_H2y3z_D2x_bb+ABZ*I_NAI_G2y2z_D2x_bb;
  Double I_NAI_Gy3z_F2xz_bb = I_NAI_Hy4z_D2x_bb+ABZ*I_NAI_Gy3z_D2x_bb;
  Double I_NAI_G4z_F2xz_bb = I_NAI_H5z_D2x_bb+ABZ*I_NAI_G4z_D2x_bb;
  Double I_NAI_G4x_F3y_bb = I_NAI_H4xy_D2y_bb+ABY*I_NAI_G4x_D2y_bb;
  Double I_NAI_G3xy_F3y_bb = I_NAI_H3x2y_D2y_bb+ABY*I_NAI_G3xy_D2y_bb;
  Double I_NAI_G3xz_F3y_bb = I_NAI_H3xyz_D2y_bb+ABY*I_NAI_G3xz_D2y_bb;
  Double I_NAI_G2x2y_F3y_bb = I_NAI_H2x3y_D2y_bb+ABY*I_NAI_G2x2y_D2y_bb;
  Double I_NAI_G2xyz_F3y_bb = I_NAI_H2x2yz_D2y_bb+ABY*I_NAI_G2xyz_D2y_bb;
  Double I_NAI_G2x2z_F3y_bb = I_NAI_H2xy2z_D2y_bb+ABY*I_NAI_G2x2z_D2y_bb;
  Double I_NAI_Gx3y_F3y_bb = I_NAI_Hx4y_D2y_bb+ABY*I_NAI_Gx3y_D2y_bb;
  Double I_NAI_Gx2yz_F3y_bb = I_NAI_Hx3yz_D2y_bb+ABY*I_NAI_Gx2yz_D2y_bb;
  Double I_NAI_Gxy2z_F3y_bb = I_NAI_Hx2y2z_D2y_bb+ABY*I_NAI_Gxy2z_D2y_bb;
  Double I_NAI_Gx3z_F3y_bb = I_NAI_Hxy3z_D2y_bb+ABY*I_NAI_Gx3z_D2y_bb;
  Double I_NAI_G4y_F3y_bb = I_NAI_H5y_D2y_bb+ABY*I_NAI_G4y_D2y_bb;
  Double I_NAI_G3yz_F3y_bb = I_NAI_H4yz_D2y_bb+ABY*I_NAI_G3yz_D2y_bb;
  Double I_NAI_G2y2z_F3y_bb = I_NAI_H3y2z_D2y_bb+ABY*I_NAI_G2y2z_D2y_bb;
  Double I_NAI_Gy3z_F3y_bb = I_NAI_H2y3z_D2y_bb+ABY*I_NAI_Gy3z_D2y_bb;
  Double I_NAI_G4z_F3y_bb = I_NAI_Hy4z_D2y_bb+ABY*I_NAI_G4z_D2y_bb;
  Double I_NAI_G4x_F2yz_bb = I_NAI_H4xz_D2y_bb+ABZ*I_NAI_G4x_D2y_bb;
  Double I_NAI_G3xy_F2yz_bb = I_NAI_H3xyz_D2y_bb+ABZ*I_NAI_G3xy_D2y_bb;
  Double I_NAI_G3xz_F2yz_bb = I_NAI_H3x2z_D2y_bb+ABZ*I_NAI_G3xz_D2y_bb;
  Double I_NAI_G2x2y_F2yz_bb = I_NAI_H2x2yz_D2y_bb+ABZ*I_NAI_G2x2y_D2y_bb;
  Double I_NAI_G2xyz_F2yz_bb = I_NAI_H2xy2z_D2y_bb+ABZ*I_NAI_G2xyz_D2y_bb;
  Double I_NAI_G2x2z_F2yz_bb = I_NAI_H2x3z_D2y_bb+ABZ*I_NAI_G2x2z_D2y_bb;
  Double I_NAI_Gx3y_F2yz_bb = I_NAI_Hx3yz_D2y_bb+ABZ*I_NAI_Gx3y_D2y_bb;
  Double I_NAI_Gx2yz_F2yz_bb = I_NAI_Hx2y2z_D2y_bb+ABZ*I_NAI_Gx2yz_D2y_bb;
  Double I_NAI_Gxy2z_F2yz_bb = I_NAI_Hxy3z_D2y_bb+ABZ*I_NAI_Gxy2z_D2y_bb;
  Double I_NAI_Gx3z_F2yz_bb = I_NAI_Hx4z_D2y_bb+ABZ*I_NAI_Gx3z_D2y_bb;
  Double I_NAI_G4y_F2yz_bb = I_NAI_H4yz_D2y_bb+ABZ*I_NAI_G4y_D2y_bb;
  Double I_NAI_G3yz_F2yz_bb = I_NAI_H3y2z_D2y_bb+ABZ*I_NAI_G3yz_D2y_bb;
  Double I_NAI_G2y2z_F2yz_bb = I_NAI_H2y3z_D2y_bb+ABZ*I_NAI_G2y2z_D2y_bb;
  Double I_NAI_Gy3z_F2yz_bb = I_NAI_Hy4z_D2y_bb+ABZ*I_NAI_Gy3z_D2y_bb;
  Double I_NAI_G4z_F2yz_bb = I_NAI_H5z_D2y_bb+ABZ*I_NAI_G4z_D2y_bb;
  Double I_NAI_G4x_F3z_bb = I_NAI_H4xz_D2z_bb+ABZ*I_NAI_G4x_D2z_bb;
  Double I_NAI_G3xy_F3z_bb = I_NAI_H3xyz_D2z_bb+ABZ*I_NAI_G3xy_D2z_bb;
  Double I_NAI_G3xz_F3z_bb = I_NAI_H3x2z_D2z_bb+ABZ*I_NAI_G3xz_D2z_bb;
  Double I_NAI_G2x2y_F3z_bb = I_NAI_H2x2yz_D2z_bb+ABZ*I_NAI_G2x2y_D2z_bb;
  Double I_NAI_G2xyz_F3z_bb = I_NAI_H2xy2z_D2z_bb+ABZ*I_NAI_G2xyz_D2z_bb;
  Double I_NAI_G2x2z_F3z_bb = I_NAI_H2x3z_D2z_bb+ABZ*I_NAI_G2x2z_D2z_bb;
  Double I_NAI_Gx3y_F3z_bb = I_NAI_Hx3yz_D2z_bb+ABZ*I_NAI_Gx3y_D2z_bb;
  Double I_NAI_Gx2yz_F3z_bb = I_NAI_Hx2y2z_D2z_bb+ABZ*I_NAI_Gx2yz_D2z_bb;
  Double I_NAI_Gxy2z_F3z_bb = I_NAI_Hxy3z_D2z_bb+ABZ*I_NAI_Gxy2z_D2z_bb;
  Double I_NAI_Gx3z_F3z_bb = I_NAI_Hx4z_D2z_bb+ABZ*I_NAI_Gx3z_D2z_bb;
  Double I_NAI_G4y_F3z_bb = I_NAI_H4yz_D2z_bb+ABZ*I_NAI_G4y_D2z_bb;
  Double I_NAI_G3yz_F3z_bb = I_NAI_H3y2z_D2z_bb+ABZ*I_NAI_G3yz_D2z_bb;
  Double I_NAI_G2y2z_F3z_bb = I_NAI_H2y3z_D2z_bb+ABZ*I_NAI_G2y2z_D2z_bb;
  Double I_NAI_Gy3z_F3z_bb = I_NAI_Hy4z_D2z_bb+ABZ*I_NAI_Gy3z_D2z_bb;
  Double I_NAI_G4z_F3z_bb = I_NAI_H5z_D2z_bb+ABZ*I_NAI_G4z_D2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_K_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 13 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S_bb
   * RHS shell quartet name: SQ_NAI_K_S_bb
   ************************************************************/
  Double I_NAI_K7x_Px_bb = I_NAI_L8x_S_bb+ABX*I_NAI_K7x_S_bb;
  Double I_NAI_K6xy_Px_bb = I_NAI_L7xy_S_bb+ABX*I_NAI_K6xy_S_bb;
  Double I_NAI_K6xz_Px_bb = I_NAI_L7xz_S_bb+ABX*I_NAI_K6xz_S_bb;
  Double I_NAI_K5x2y_Px_bb = I_NAI_L6x2y_S_bb+ABX*I_NAI_K5x2y_S_bb;
  Double I_NAI_K5xyz_Px_bb = I_NAI_L6xyz_S_bb+ABX*I_NAI_K5xyz_S_bb;
  Double I_NAI_K5x2z_Px_bb = I_NAI_L6x2z_S_bb+ABX*I_NAI_K5x2z_S_bb;
  Double I_NAI_K4x3y_Px_bb = I_NAI_L5x3y_S_bb+ABX*I_NAI_K4x3y_S_bb;
  Double I_NAI_K4x2yz_Px_bb = I_NAI_L5x2yz_S_bb+ABX*I_NAI_K4x2yz_S_bb;
  Double I_NAI_K4xy2z_Px_bb = I_NAI_L5xy2z_S_bb+ABX*I_NAI_K4xy2z_S_bb;
  Double I_NAI_K4x3z_Px_bb = I_NAI_L5x3z_S_bb+ABX*I_NAI_K4x3z_S_bb;
  Double I_NAI_K3x4y_Px_bb = I_NAI_L4x4y_S_bb+ABX*I_NAI_K3x4y_S_bb;
  Double I_NAI_K3x3yz_Px_bb = I_NAI_L4x3yz_S_bb+ABX*I_NAI_K3x3yz_S_bb;
  Double I_NAI_K3x2y2z_Px_bb = I_NAI_L4x2y2z_S_bb+ABX*I_NAI_K3x2y2z_S_bb;
  Double I_NAI_K3xy3z_Px_bb = I_NAI_L4xy3z_S_bb+ABX*I_NAI_K3xy3z_S_bb;
  Double I_NAI_K3x4z_Px_bb = I_NAI_L4x4z_S_bb+ABX*I_NAI_K3x4z_S_bb;
  Double I_NAI_K2x5y_Px_bb = I_NAI_L3x5y_S_bb+ABX*I_NAI_K2x5y_S_bb;
  Double I_NAI_K2x4yz_Px_bb = I_NAI_L3x4yz_S_bb+ABX*I_NAI_K2x4yz_S_bb;
  Double I_NAI_K2x3y2z_Px_bb = I_NAI_L3x3y2z_S_bb+ABX*I_NAI_K2x3y2z_S_bb;
  Double I_NAI_K2x2y3z_Px_bb = I_NAI_L3x2y3z_S_bb+ABX*I_NAI_K2x2y3z_S_bb;
  Double I_NAI_K2xy4z_Px_bb = I_NAI_L3xy4z_S_bb+ABX*I_NAI_K2xy4z_S_bb;
  Double I_NAI_K2x5z_Px_bb = I_NAI_L3x5z_S_bb+ABX*I_NAI_K2x5z_S_bb;
  Double I_NAI_Kx6y_Px_bb = I_NAI_L2x6y_S_bb+ABX*I_NAI_Kx6y_S_bb;
  Double I_NAI_Kx5yz_Px_bb = I_NAI_L2x5yz_S_bb+ABX*I_NAI_Kx5yz_S_bb;
  Double I_NAI_Kx4y2z_Px_bb = I_NAI_L2x4y2z_S_bb+ABX*I_NAI_Kx4y2z_S_bb;
  Double I_NAI_Kx3y3z_Px_bb = I_NAI_L2x3y3z_S_bb+ABX*I_NAI_Kx3y3z_S_bb;
  Double I_NAI_Kx2y4z_Px_bb = I_NAI_L2x2y4z_S_bb+ABX*I_NAI_Kx2y4z_S_bb;
  Double I_NAI_Kxy5z_Px_bb = I_NAI_L2xy5z_S_bb+ABX*I_NAI_Kxy5z_S_bb;
  Double I_NAI_Kx6z_Px_bb = I_NAI_L2x6z_S_bb+ABX*I_NAI_Kx6z_S_bb;
  Double I_NAI_K6yz_Px_bb = I_NAI_Lx6yz_S_bb+ABX*I_NAI_K6yz_S_bb;
  Double I_NAI_K5y2z_Px_bb = I_NAI_Lx5y2z_S_bb+ABX*I_NAI_K5y2z_S_bb;
  Double I_NAI_K4y3z_Px_bb = I_NAI_Lx4y3z_S_bb+ABX*I_NAI_K4y3z_S_bb;
  Double I_NAI_K3y4z_Px_bb = I_NAI_Lx3y4z_S_bb+ABX*I_NAI_K3y4z_S_bb;
  Double I_NAI_K2y5z_Px_bb = I_NAI_Lx2y5z_S_bb+ABX*I_NAI_K2y5z_S_bb;
  Double I_NAI_Ky6z_Px_bb = I_NAI_Lxy6z_S_bb+ABX*I_NAI_Ky6z_S_bb;
  Double I_NAI_K6xy_Py_bb = I_NAI_L6x2y_S_bb+ABY*I_NAI_K6xy_S_bb;
  Double I_NAI_K5x2y_Py_bb = I_NAI_L5x3y_S_bb+ABY*I_NAI_K5x2y_S_bb;
  Double I_NAI_K5xyz_Py_bb = I_NAI_L5x2yz_S_bb+ABY*I_NAI_K5xyz_S_bb;
  Double I_NAI_K5x2z_Py_bb = I_NAI_L5xy2z_S_bb+ABY*I_NAI_K5x2z_S_bb;
  Double I_NAI_K4x3y_Py_bb = I_NAI_L4x4y_S_bb+ABY*I_NAI_K4x3y_S_bb;
  Double I_NAI_K4x2yz_Py_bb = I_NAI_L4x3yz_S_bb+ABY*I_NAI_K4x2yz_S_bb;
  Double I_NAI_K4xy2z_Py_bb = I_NAI_L4x2y2z_S_bb+ABY*I_NAI_K4xy2z_S_bb;
  Double I_NAI_K4x3z_Py_bb = I_NAI_L4xy3z_S_bb+ABY*I_NAI_K4x3z_S_bb;
  Double I_NAI_K3x4y_Py_bb = I_NAI_L3x5y_S_bb+ABY*I_NAI_K3x4y_S_bb;
  Double I_NAI_K3x3yz_Py_bb = I_NAI_L3x4yz_S_bb+ABY*I_NAI_K3x3yz_S_bb;
  Double I_NAI_K3x2y2z_Py_bb = I_NAI_L3x3y2z_S_bb+ABY*I_NAI_K3x2y2z_S_bb;
  Double I_NAI_K3xy3z_Py_bb = I_NAI_L3x2y3z_S_bb+ABY*I_NAI_K3xy3z_S_bb;
  Double I_NAI_K3x4z_Py_bb = I_NAI_L3xy4z_S_bb+ABY*I_NAI_K3x4z_S_bb;
  Double I_NAI_K2x5y_Py_bb = I_NAI_L2x6y_S_bb+ABY*I_NAI_K2x5y_S_bb;
  Double I_NAI_K2x4yz_Py_bb = I_NAI_L2x5yz_S_bb+ABY*I_NAI_K2x4yz_S_bb;
  Double I_NAI_K2x3y2z_Py_bb = I_NAI_L2x4y2z_S_bb+ABY*I_NAI_K2x3y2z_S_bb;
  Double I_NAI_K2x2y3z_Py_bb = I_NAI_L2x3y3z_S_bb+ABY*I_NAI_K2x2y3z_S_bb;
  Double I_NAI_K2xy4z_Py_bb = I_NAI_L2x2y4z_S_bb+ABY*I_NAI_K2xy4z_S_bb;
  Double I_NAI_K2x5z_Py_bb = I_NAI_L2xy5z_S_bb+ABY*I_NAI_K2x5z_S_bb;
  Double I_NAI_Kx6y_Py_bb = I_NAI_Lx7y_S_bb+ABY*I_NAI_Kx6y_S_bb;
  Double I_NAI_Kx5yz_Py_bb = I_NAI_Lx6yz_S_bb+ABY*I_NAI_Kx5yz_S_bb;
  Double I_NAI_Kx4y2z_Py_bb = I_NAI_Lx5y2z_S_bb+ABY*I_NAI_Kx4y2z_S_bb;
  Double I_NAI_Kx3y3z_Py_bb = I_NAI_Lx4y3z_S_bb+ABY*I_NAI_Kx3y3z_S_bb;
  Double I_NAI_Kx2y4z_Py_bb = I_NAI_Lx3y4z_S_bb+ABY*I_NAI_Kx2y4z_S_bb;
  Double I_NAI_Kxy5z_Py_bb = I_NAI_Lx2y5z_S_bb+ABY*I_NAI_Kxy5z_S_bb;
  Double I_NAI_Kx6z_Py_bb = I_NAI_Lxy6z_S_bb+ABY*I_NAI_Kx6z_S_bb;
  Double I_NAI_K7y_Py_bb = I_NAI_L8y_S_bb+ABY*I_NAI_K7y_S_bb;
  Double I_NAI_K6yz_Py_bb = I_NAI_L7yz_S_bb+ABY*I_NAI_K6yz_S_bb;
  Double I_NAI_K5y2z_Py_bb = I_NAI_L6y2z_S_bb+ABY*I_NAI_K5y2z_S_bb;
  Double I_NAI_K4y3z_Py_bb = I_NAI_L5y3z_S_bb+ABY*I_NAI_K4y3z_S_bb;
  Double I_NAI_K3y4z_Py_bb = I_NAI_L4y4z_S_bb+ABY*I_NAI_K3y4z_S_bb;
  Double I_NAI_K2y5z_Py_bb = I_NAI_L3y5z_S_bb+ABY*I_NAI_K2y5z_S_bb;
  Double I_NAI_Ky6z_Py_bb = I_NAI_L2y6z_S_bb+ABY*I_NAI_Ky6z_S_bb;
  Double I_NAI_K6xz_Pz_bb = I_NAI_L6x2z_S_bb+ABZ*I_NAI_K6xz_S_bb;
  Double I_NAI_K5xyz_Pz_bb = I_NAI_L5xy2z_S_bb+ABZ*I_NAI_K5xyz_S_bb;
  Double I_NAI_K5x2z_Pz_bb = I_NAI_L5x3z_S_bb+ABZ*I_NAI_K5x2z_S_bb;
  Double I_NAI_K4x2yz_Pz_bb = I_NAI_L4x2y2z_S_bb+ABZ*I_NAI_K4x2yz_S_bb;
  Double I_NAI_K4xy2z_Pz_bb = I_NAI_L4xy3z_S_bb+ABZ*I_NAI_K4xy2z_S_bb;
  Double I_NAI_K4x3z_Pz_bb = I_NAI_L4x4z_S_bb+ABZ*I_NAI_K4x3z_S_bb;
  Double I_NAI_K3x3yz_Pz_bb = I_NAI_L3x3y2z_S_bb+ABZ*I_NAI_K3x3yz_S_bb;
  Double I_NAI_K3x2y2z_Pz_bb = I_NAI_L3x2y3z_S_bb+ABZ*I_NAI_K3x2y2z_S_bb;
  Double I_NAI_K3xy3z_Pz_bb = I_NAI_L3xy4z_S_bb+ABZ*I_NAI_K3xy3z_S_bb;
  Double I_NAI_K3x4z_Pz_bb = I_NAI_L3x5z_S_bb+ABZ*I_NAI_K3x4z_S_bb;
  Double I_NAI_K2x4yz_Pz_bb = I_NAI_L2x4y2z_S_bb+ABZ*I_NAI_K2x4yz_S_bb;
  Double I_NAI_K2x3y2z_Pz_bb = I_NAI_L2x3y3z_S_bb+ABZ*I_NAI_K2x3y2z_S_bb;
  Double I_NAI_K2x2y3z_Pz_bb = I_NAI_L2x2y4z_S_bb+ABZ*I_NAI_K2x2y3z_S_bb;
  Double I_NAI_K2xy4z_Pz_bb = I_NAI_L2xy5z_S_bb+ABZ*I_NAI_K2xy4z_S_bb;
  Double I_NAI_K2x5z_Pz_bb = I_NAI_L2x6z_S_bb+ABZ*I_NAI_K2x5z_S_bb;
  Double I_NAI_Kx5yz_Pz_bb = I_NAI_Lx5y2z_S_bb+ABZ*I_NAI_Kx5yz_S_bb;
  Double I_NAI_Kx4y2z_Pz_bb = I_NAI_Lx4y3z_S_bb+ABZ*I_NAI_Kx4y2z_S_bb;
  Double I_NAI_Kx3y3z_Pz_bb = I_NAI_Lx3y4z_S_bb+ABZ*I_NAI_Kx3y3z_S_bb;
  Double I_NAI_Kx2y4z_Pz_bb = I_NAI_Lx2y5z_S_bb+ABZ*I_NAI_Kx2y4z_S_bb;
  Double I_NAI_Kxy5z_Pz_bb = I_NAI_Lxy6z_S_bb+ABZ*I_NAI_Kxy5z_S_bb;
  Double I_NAI_Kx6z_Pz_bb = I_NAI_Lx7z_S_bb+ABZ*I_NAI_Kx6z_S_bb;
  Double I_NAI_K6yz_Pz_bb = I_NAI_L6y2z_S_bb+ABZ*I_NAI_K6yz_S_bb;
  Double I_NAI_K5y2z_Pz_bb = I_NAI_L5y3z_S_bb+ABZ*I_NAI_K5y2z_S_bb;
  Double I_NAI_K4y3z_Pz_bb = I_NAI_L4y4z_S_bb+ABZ*I_NAI_K4y3z_S_bb;
  Double I_NAI_K3y4z_Pz_bb = I_NAI_L3y5z_S_bb+ABZ*I_NAI_K3y4z_S_bb;
  Double I_NAI_K2y5z_Pz_bb = I_NAI_L2y6z_S_bb+ABZ*I_NAI_K2y5z_S_bb;
  Double I_NAI_Ky6z_Pz_bb = I_NAI_Ly7z_S_bb+ABZ*I_NAI_Ky6z_S_bb;
  Double I_NAI_K7z_Pz_bb = I_NAI_L8z_S_bb+ABZ*I_NAI_K7z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 84 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_bb
   * RHS shell quartet name: SQ_NAI_I_P_bb
   ************************************************************/
  Double I_NAI_I6x_D2x_bb = I_NAI_K7x_Px_bb+ABX*I_NAI_I6x_Px_bb;
  Double I_NAI_I5xy_D2x_bb = I_NAI_K6xy_Px_bb+ABX*I_NAI_I5xy_Px_bb;
  Double I_NAI_I5xz_D2x_bb = I_NAI_K6xz_Px_bb+ABX*I_NAI_I5xz_Px_bb;
  Double I_NAI_I4x2y_D2x_bb = I_NAI_K5x2y_Px_bb+ABX*I_NAI_I4x2y_Px_bb;
  Double I_NAI_I4xyz_D2x_bb = I_NAI_K5xyz_Px_bb+ABX*I_NAI_I4xyz_Px_bb;
  Double I_NAI_I4x2z_D2x_bb = I_NAI_K5x2z_Px_bb+ABX*I_NAI_I4x2z_Px_bb;
  Double I_NAI_I3x3y_D2x_bb = I_NAI_K4x3y_Px_bb+ABX*I_NAI_I3x3y_Px_bb;
  Double I_NAI_I3x2yz_D2x_bb = I_NAI_K4x2yz_Px_bb+ABX*I_NAI_I3x2yz_Px_bb;
  Double I_NAI_I3xy2z_D2x_bb = I_NAI_K4xy2z_Px_bb+ABX*I_NAI_I3xy2z_Px_bb;
  Double I_NAI_I3x3z_D2x_bb = I_NAI_K4x3z_Px_bb+ABX*I_NAI_I3x3z_Px_bb;
  Double I_NAI_I2x4y_D2x_bb = I_NAI_K3x4y_Px_bb+ABX*I_NAI_I2x4y_Px_bb;
  Double I_NAI_I2x3yz_D2x_bb = I_NAI_K3x3yz_Px_bb+ABX*I_NAI_I2x3yz_Px_bb;
  Double I_NAI_I2x2y2z_D2x_bb = I_NAI_K3x2y2z_Px_bb+ABX*I_NAI_I2x2y2z_Px_bb;
  Double I_NAI_I2xy3z_D2x_bb = I_NAI_K3xy3z_Px_bb+ABX*I_NAI_I2xy3z_Px_bb;
  Double I_NAI_I2x4z_D2x_bb = I_NAI_K3x4z_Px_bb+ABX*I_NAI_I2x4z_Px_bb;
  Double I_NAI_Ix5y_D2x_bb = I_NAI_K2x5y_Px_bb+ABX*I_NAI_Ix5y_Px_bb;
  Double I_NAI_Ix4yz_D2x_bb = I_NAI_K2x4yz_Px_bb+ABX*I_NAI_Ix4yz_Px_bb;
  Double I_NAI_Ix3y2z_D2x_bb = I_NAI_K2x3y2z_Px_bb+ABX*I_NAI_Ix3y2z_Px_bb;
  Double I_NAI_Ix2y3z_D2x_bb = I_NAI_K2x2y3z_Px_bb+ABX*I_NAI_Ix2y3z_Px_bb;
  Double I_NAI_Ixy4z_D2x_bb = I_NAI_K2xy4z_Px_bb+ABX*I_NAI_Ixy4z_Px_bb;
  Double I_NAI_Ix5z_D2x_bb = I_NAI_K2x5z_Px_bb+ABX*I_NAI_Ix5z_Px_bb;
  Double I_NAI_I6y_D2x_bb = I_NAI_Kx6y_Px_bb+ABX*I_NAI_I6y_Px_bb;
  Double I_NAI_I5yz_D2x_bb = I_NAI_Kx5yz_Px_bb+ABX*I_NAI_I5yz_Px_bb;
  Double I_NAI_I4y2z_D2x_bb = I_NAI_Kx4y2z_Px_bb+ABX*I_NAI_I4y2z_Px_bb;
  Double I_NAI_I3y3z_D2x_bb = I_NAI_Kx3y3z_Px_bb+ABX*I_NAI_I3y3z_Px_bb;
  Double I_NAI_I2y4z_D2x_bb = I_NAI_Kx2y4z_Px_bb+ABX*I_NAI_I2y4z_Px_bb;
  Double I_NAI_Iy5z_D2x_bb = I_NAI_Kxy5z_Px_bb+ABX*I_NAI_Iy5z_Px_bb;
  Double I_NAI_I6z_D2x_bb = I_NAI_Kx6z_Px_bb+ABX*I_NAI_I6z_Px_bb;
  Double I_NAI_I6x_D2y_bb = I_NAI_K6xy_Py_bb+ABY*I_NAI_I6x_Py_bb;
  Double I_NAI_I5xy_D2y_bb = I_NAI_K5x2y_Py_bb+ABY*I_NAI_I5xy_Py_bb;
  Double I_NAI_I5xz_D2y_bb = I_NAI_K5xyz_Py_bb+ABY*I_NAI_I5xz_Py_bb;
  Double I_NAI_I4x2y_D2y_bb = I_NAI_K4x3y_Py_bb+ABY*I_NAI_I4x2y_Py_bb;
  Double I_NAI_I4xyz_D2y_bb = I_NAI_K4x2yz_Py_bb+ABY*I_NAI_I4xyz_Py_bb;
  Double I_NAI_I4x2z_D2y_bb = I_NAI_K4xy2z_Py_bb+ABY*I_NAI_I4x2z_Py_bb;
  Double I_NAI_I3x3y_D2y_bb = I_NAI_K3x4y_Py_bb+ABY*I_NAI_I3x3y_Py_bb;
  Double I_NAI_I3x2yz_D2y_bb = I_NAI_K3x3yz_Py_bb+ABY*I_NAI_I3x2yz_Py_bb;
  Double I_NAI_I3xy2z_D2y_bb = I_NAI_K3x2y2z_Py_bb+ABY*I_NAI_I3xy2z_Py_bb;
  Double I_NAI_I3x3z_D2y_bb = I_NAI_K3xy3z_Py_bb+ABY*I_NAI_I3x3z_Py_bb;
  Double I_NAI_I2x4y_D2y_bb = I_NAI_K2x5y_Py_bb+ABY*I_NAI_I2x4y_Py_bb;
  Double I_NAI_I2x3yz_D2y_bb = I_NAI_K2x4yz_Py_bb+ABY*I_NAI_I2x3yz_Py_bb;
  Double I_NAI_I2x2y2z_D2y_bb = I_NAI_K2x3y2z_Py_bb+ABY*I_NAI_I2x2y2z_Py_bb;
  Double I_NAI_I2xy3z_D2y_bb = I_NAI_K2x2y3z_Py_bb+ABY*I_NAI_I2xy3z_Py_bb;
  Double I_NAI_I2x4z_D2y_bb = I_NAI_K2xy4z_Py_bb+ABY*I_NAI_I2x4z_Py_bb;
  Double I_NAI_Ix5y_D2y_bb = I_NAI_Kx6y_Py_bb+ABY*I_NAI_Ix5y_Py_bb;
  Double I_NAI_Ix4yz_D2y_bb = I_NAI_Kx5yz_Py_bb+ABY*I_NAI_Ix4yz_Py_bb;
  Double I_NAI_Ix3y2z_D2y_bb = I_NAI_Kx4y2z_Py_bb+ABY*I_NAI_Ix3y2z_Py_bb;
  Double I_NAI_Ix2y3z_D2y_bb = I_NAI_Kx3y3z_Py_bb+ABY*I_NAI_Ix2y3z_Py_bb;
  Double I_NAI_Ixy4z_D2y_bb = I_NAI_Kx2y4z_Py_bb+ABY*I_NAI_Ixy4z_Py_bb;
  Double I_NAI_Ix5z_D2y_bb = I_NAI_Kxy5z_Py_bb+ABY*I_NAI_Ix5z_Py_bb;
  Double I_NAI_I6y_D2y_bb = I_NAI_K7y_Py_bb+ABY*I_NAI_I6y_Py_bb;
  Double I_NAI_I5yz_D2y_bb = I_NAI_K6yz_Py_bb+ABY*I_NAI_I5yz_Py_bb;
  Double I_NAI_I4y2z_D2y_bb = I_NAI_K5y2z_Py_bb+ABY*I_NAI_I4y2z_Py_bb;
  Double I_NAI_I3y3z_D2y_bb = I_NAI_K4y3z_Py_bb+ABY*I_NAI_I3y3z_Py_bb;
  Double I_NAI_I2y4z_D2y_bb = I_NAI_K3y4z_Py_bb+ABY*I_NAI_I2y4z_Py_bb;
  Double I_NAI_Iy5z_D2y_bb = I_NAI_K2y5z_Py_bb+ABY*I_NAI_Iy5z_Py_bb;
  Double I_NAI_I6z_D2y_bb = I_NAI_Ky6z_Py_bb+ABY*I_NAI_I6z_Py_bb;
  Double I_NAI_I6x_D2z_bb = I_NAI_K6xz_Pz_bb+ABZ*I_NAI_I6x_Pz_bb;
  Double I_NAI_I5xy_D2z_bb = I_NAI_K5xyz_Pz_bb+ABZ*I_NAI_I5xy_Pz_bb;
  Double I_NAI_I5xz_D2z_bb = I_NAI_K5x2z_Pz_bb+ABZ*I_NAI_I5xz_Pz_bb;
  Double I_NAI_I4x2y_D2z_bb = I_NAI_K4x2yz_Pz_bb+ABZ*I_NAI_I4x2y_Pz_bb;
  Double I_NAI_I4xyz_D2z_bb = I_NAI_K4xy2z_Pz_bb+ABZ*I_NAI_I4xyz_Pz_bb;
  Double I_NAI_I4x2z_D2z_bb = I_NAI_K4x3z_Pz_bb+ABZ*I_NAI_I4x2z_Pz_bb;
  Double I_NAI_I3x3y_D2z_bb = I_NAI_K3x3yz_Pz_bb+ABZ*I_NAI_I3x3y_Pz_bb;
  Double I_NAI_I3x2yz_D2z_bb = I_NAI_K3x2y2z_Pz_bb+ABZ*I_NAI_I3x2yz_Pz_bb;
  Double I_NAI_I3xy2z_D2z_bb = I_NAI_K3xy3z_Pz_bb+ABZ*I_NAI_I3xy2z_Pz_bb;
  Double I_NAI_I3x3z_D2z_bb = I_NAI_K3x4z_Pz_bb+ABZ*I_NAI_I3x3z_Pz_bb;
  Double I_NAI_I2x4y_D2z_bb = I_NAI_K2x4yz_Pz_bb+ABZ*I_NAI_I2x4y_Pz_bb;
  Double I_NAI_I2x3yz_D2z_bb = I_NAI_K2x3y2z_Pz_bb+ABZ*I_NAI_I2x3yz_Pz_bb;
  Double I_NAI_I2x2y2z_D2z_bb = I_NAI_K2x2y3z_Pz_bb+ABZ*I_NAI_I2x2y2z_Pz_bb;
  Double I_NAI_I2xy3z_D2z_bb = I_NAI_K2xy4z_Pz_bb+ABZ*I_NAI_I2xy3z_Pz_bb;
  Double I_NAI_I2x4z_D2z_bb = I_NAI_K2x5z_Pz_bb+ABZ*I_NAI_I2x4z_Pz_bb;
  Double I_NAI_Ix5y_D2z_bb = I_NAI_Kx5yz_Pz_bb+ABZ*I_NAI_Ix5y_Pz_bb;
  Double I_NAI_Ix4yz_D2z_bb = I_NAI_Kx4y2z_Pz_bb+ABZ*I_NAI_Ix4yz_Pz_bb;
  Double I_NAI_Ix3y2z_D2z_bb = I_NAI_Kx3y3z_Pz_bb+ABZ*I_NAI_Ix3y2z_Pz_bb;
  Double I_NAI_Ix2y3z_D2z_bb = I_NAI_Kx2y4z_Pz_bb+ABZ*I_NAI_Ix2y3z_Pz_bb;
  Double I_NAI_Ixy4z_D2z_bb = I_NAI_Kxy5z_Pz_bb+ABZ*I_NAI_Ixy4z_Pz_bb;
  Double I_NAI_Ix5z_D2z_bb = I_NAI_Kx6z_Pz_bb+ABZ*I_NAI_Ix5z_Pz_bb;
  Double I_NAI_I6y_D2z_bb = I_NAI_K6yz_Pz_bb+ABZ*I_NAI_I6y_Pz_bb;
  Double I_NAI_I5yz_D2z_bb = I_NAI_K5y2z_Pz_bb+ABZ*I_NAI_I5yz_Pz_bb;
  Double I_NAI_I4y2z_D2z_bb = I_NAI_K4y3z_Pz_bb+ABZ*I_NAI_I4y2z_Pz_bb;
  Double I_NAI_I3y3z_D2z_bb = I_NAI_K3y4z_Pz_bb+ABZ*I_NAI_I3y3z_Pz_bb;
  Double I_NAI_I2y4z_D2z_bb = I_NAI_K2y5z_Pz_bb+ABZ*I_NAI_I2y4z_Pz_bb;
  Double I_NAI_Iy5z_D2z_bb = I_NAI_Ky6z_Pz_bb+ABZ*I_NAI_Iy5z_Pz_bb;
  Double I_NAI_I6z_D2z_bb = I_NAI_K7z_Pz_bb+ABZ*I_NAI_I6z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 87 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_bb
   * RHS shell quartet name: SQ_NAI_H_D_bb
   ************************************************************/
  Double I_NAI_H5x_F3x_bb = I_NAI_I6x_D2x_bb+ABX*I_NAI_H5x_D2x_bb;
  Double I_NAI_H4xy_F3x_bb = I_NAI_I5xy_D2x_bb+ABX*I_NAI_H4xy_D2x_bb;
  Double I_NAI_H4xz_F3x_bb = I_NAI_I5xz_D2x_bb+ABX*I_NAI_H4xz_D2x_bb;
  Double I_NAI_H3x2y_F3x_bb = I_NAI_I4x2y_D2x_bb+ABX*I_NAI_H3x2y_D2x_bb;
  Double I_NAI_H3xyz_F3x_bb = I_NAI_I4xyz_D2x_bb+ABX*I_NAI_H3xyz_D2x_bb;
  Double I_NAI_H3x2z_F3x_bb = I_NAI_I4x2z_D2x_bb+ABX*I_NAI_H3x2z_D2x_bb;
  Double I_NAI_H2x3y_F3x_bb = I_NAI_I3x3y_D2x_bb+ABX*I_NAI_H2x3y_D2x_bb;
  Double I_NAI_H2x2yz_F3x_bb = I_NAI_I3x2yz_D2x_bb+ABX*I_NAI_H2x2yz_D2x_bb;
  Double I_NAI_H2xy2z_F3x_bb = I_NAI_I3xy2z_D2x_bb+ABX*I_NAI_H2xy2z_D2x_bb;
  Double I_NAI_H2x3z_F3x_bb = I_NAI_I3x3z_D2x_bb+ABX*I_NAI_H2x3z_D2x_bb;
  Double I_NAI_Hx4y_F3x_bb = I_NAI_I2x4y_D2x_bb+ABX*I_NAI_Hx4y_D2x_bb;
  Double I_NAI_Hx3yz_F3x_bb = I_NAI_I2x3yz_D2x_bb+ABX*I_NAI_Hx3yz_D2x_bb;
  Double I_NAI_Hx2y2z_F3x_bb = I_NAI_I2x2y2z_D2x_bb+ABX*I_NAI_Hx2y2z_D2x_bb;
  Double I_NAI_Hxy3z_F3x_bb = I_NAI_I2xy3z_D2x_bb+ABX*I_NAI_Hxy3z_D2x_bb;
  Double I_NAI_Hx4z_F3x_bb = I_NAI_I2x4z_D2x_bb+ABX*I_NAI_Hx4z_D2x_bb;
  Double I_NAI_H5y_F3x_bb = I_NAI_Ix5y_D2x_bb+ABX*I_NAI_H5y_D2x_bb;
  Double I_NAI_H4yz_F3x_bb = I_NAI_Ix4yz_D2x_bb+ABX*I_NAI_H4yz_D2x_bb;
  Double I_NAI_H3y2z_F3x_bb = I_NAI_Ix3y2z_D2x_bb+ABX*I_NAI_H3y2z_D2x_bb;
  Double I_NAI_H2y3z_F3x_bb = I_NAI_Ix2y3z_D2x_bb+ABX*I_NAI_H2y3z_D2x_bb;
  Double I_NAI_Hy4z_F3x_bb = I_NAI_Ixy4z_D2x_bb+ABX*I_NAI_Hy4z_D2x_bb;
  Double I_NAI_H5z_F3x_bb = I_NAI_Ix5z_D2x_bb+ABX*I_NAI_H5z_D2x_bb;
  Double I_NAI_H4xy_F2xy_bb = I_NAI_I4x2y_D2x_bb+ABY*I_NAI_H4xy_D2x_bb;
  Double I_NAI_H4xz_F2xy_bb = I_NAI_I4xyz_D2x_bb+ABY*I_NAI_H4xz_D2x_bb;
  Double I_NAI_H3x2y_F2xy_bb = I_NAI_I3x3y_D2x_bb+ABY*I_NAI_H3x2y_D2x_bb;
  Double I_NAI_H3xyz_F2xy_bb = I_NAI_I3x2yz_D2x_bb+ABY*I_NAI_H3xyz_D2x_bb;
  Double I_NAI_H3x2z_F2xy_bb = I_NAI_I3xy2z_D2x_bb+ABY*I_NAI_H3x2z_D2x_bb;
  Double I_NAI_H2x3y_F2xy_bb = I_NAI_I2x4y_D2x_bb+ABY*I_NAI_H2x3y_D2x_bb;
  Double I_NAI_H2x2yz_F2xy_bb = I_NAI_I2x3yz_D2x_bb+ABY*I_NAI_H2x2yz_D2x_bb;
  Double I_NAI_H2xy2z_F2xy_bb = I_NAI_I2x2y2z_D2x_bb+ABY*I_NAI_H2xy2z_D2x_bb;
  Double I_NAI_H2x3z_F2xy_bb = I_NAI_I2xy3z_D2x_bb+ABY*I_NAI_H2x3z_D2x_bb;
  Double I_NAI_Hx4y_F2xy_bb = I_NAI_Ix5y_D2x_bb+ABY*I_NAI_Hx4y_D2x_bb;
  Double I_NAI_Hx3yz_F2xy_bb = I_NAI_Ix4yz_D2x_bb+ABY*I_NAI_Hx3yz_D2x_bb;
  Double I_NAI_Hx2y2z_F2xy_bb = I_NAI_Ix3y2z_D2x_bb+ABY*I_NAI_Hx2y2z_D2x_bb;
  Double I_NAI_Hxy3z_F2xy_bb = I_NAI_Ix2y3z_D2x_bb+ABY*I_NAI_Hxy3z_D2x_bb;
  Double I_NAI_Hx4z_F2xy_bb = I_NAI_Ixy4z_D2x_bb+ABY*I_NAI_Hx4z_D2x_bb;
  Double I_NAI_H5y_F2xy_bb = I_NAI_I6y_D2x_bb+ABY*I_NAI_H5y_D2x_bb;
  Double I_NAI_H4yz_F2xy_bb = I_NAI_I5yz_D2x_bb+ABY*I_NAI_H4yz_D2x_bb;
  Double I_NAI_H3y2z_F2xy_bb = I_NAI_I4y2z_D2x_bb+ABY*I_NAI_H3y2z_D2x_bb;
  Double I_NAI_H2y3z_F2xy_bb = I_NAI_I3y3z_D2x_bb+ABY*I_NAI_H2y3z_D2x_bb;
  Double I_NAI_Hy4z_F2xy_bb = I_NAI_I2y4z_D2x_bb+ABY*I_NAI_Hy4z_D2x_bb;
  Double I_NAI_H5z_F2xy_bb = I_NAI_Iy5z_D2x_bb+ABY*I_NAI_H5z_D2x_bb;
  Double I_NAI_H4xy_F2xz_bb = I_NAI_I4xyz_D2x_bb+ABZ*I_NAI_H4xy_D2x_bb;
  Double I_NAI_H4xz_F2xz_bb = I_NAI_I4x2z_D2x_bb+ABZ*I_NAI_H4xz_D2x_bb;
  Double I_NAI_H3x2y_F2xz_bb = I_NAI_I3x2yz_D2x_bb+ABZ*I_NAI_H3x2y_D2x_bb;
  Double I_NAI_H3xyz_F2xz_bb = I_NAI_I3xy2z_D2x_bb+ABZ*I_NAI_H3xyz_D2x_bb;
  Double I_NAI_H3x2z_F2xz_bb = I_NAI_I3x3z_D2x_bb+ABZ*I_NAI_H3x2z_D2x_bb;
  Double I_NAI_H2x3y_F2xz_bb = I_NAI_I2x3yz_D2x_bb+ABZ*I_NAI_H2x3y_D2x_bb;
  Double I_NAI_H2x2yz_F2xz_bb = I_NAI_I2x2y2z_D2x_bb+ABZ*I_NAI_H2x2yz_D2x_bb;
  Double I_NAI_H2xy2z_F2xz_bb = I_NAI_I2xy3z_D2x_bb+ABZ*I_NAI_H2xy2z_D2x_bb;
  Double I_NAI_H2x3z_F2xz_bb = I_NAI_I2x4z_D2x_bb+ABZ*I_NAI_H2x3z_D2x_bb;
  Double I_NAI_Hx4y_F2xz_bb = I_NAI_Ix4yz_D2x_bb+ABZ*I_NAI_Hx4y_D2x_bb;
  Double I_NAI_Hx3yz_F2xz_bb = I_NAI_Ix3y2z_D2x_bb+ABZ*I_NAI_Hx3yz_D2x_bb;
  Double I_NAI_Hx2y2z_F2xz_bb = I_NAI_Ix2y3z_D2x_bb+ABZ*I_NAI_Hx2y2z_D2x_bb;
  Double I_NAI_Hxy3z_F2xz_bb = I_NAI_Ixy4z_D2x_bb+ABZ*I_NAI_Hxy3z_D2x_bb;
  Double I_NAI_Hx4z_F2xz_bb = I_NAI_Ix5z_D2x_bb+ABZ*I_NAI_Hx4z_D2x_bb;
  Double I_NAI_H5y_F2xz_bb = I_NAI_I5yz_D2x_bb+ABZ*I_NAI_H5y_D2x_bb;
  Double I_NAI_H4yz_F2xz_bb = I_NAI_I4y2z_D2x_bb+ABZ*I_NAI_H4yz_D2x_bb;
  Double I_NAI_H3y2z_F2xz_bb = I_NAI_I3y3z_D2x_bb+ABZ*I_NAI_H3y2z_D2x_bb;
  Double I_NAI_H2y3z_F2xz_bb = I_NAI_I2y4z_D2x_bb+ABZ*I_NAI_H2y3z_D2x_bb;
  Double I_NAI_Hy4z_F2xz_bb = I_NAI_Iy5z_D2x_bb+ABZ*I_NAI_Hy4z_D2x_bb;
  Double I_NAI_H5z_F2xz_bb = I_NAI_I6z_D2x_bb+ABZ*I_NAI_H5z_D2x_bb;
  Double I_NAI_H5x_F3y_bb = I_NAI_I5xy_D2y_bb+ABY*I_NAI_H5x_D2y_bb;
  Double I_NAI_H4xy_F3y_bb = I_NAI_I4x2y_D2y_bb+ABY*I_NAI_H4xy_D2y_bb;
  Double I_NAI_H4xz_F3y_bb = I_NAI_I4xyz_D2y_bb+ABY*I_NAI_H4xz_D2y_bb;
  Double I_NAI_H3x2y_F3y_bb = I_NAI_I3x3y_D2y_bb+ABY*I_NAI_H3x2y_D2y_bb;
  Double I_NAI_H3xyz_F3y_bb = I_NAI_I3x2yz_D2y_bb+ABY*I_NAI_H3xyz_D2y_bb;
  Double I_NAI_H3x2z_F3y_bb = I_NAI_I3xy2z_D2y_bb+ABY*I_NAI_H3x2z_D2y_bb;
  Double I_NAI_H2x3y_F3y_bb = I_NAI_I2x4y_D2y_bb+ABY*I_NAI_H2x3y_D2y_bb;
  Double I_NAI_H2x2yz_F3y_bb = I_NAI_I2x3yz_D2y_bb+ABY*I_NAI_H2x2yz_D2y_bb;
  Double I_NAI_H2xy2z_F3y_bb = I_NAI_I2x2y2z_D2y_bb+ABY*I_NAI_H2xy2z_D2y_bb;
  Double I_NAI_H2x3z_F3y_bb = I_NAI_I2xy3z_D2y_bb+ABY*I_NAI_H2x3z_D2y_bb;
  Double I_NAI_Hx4y_F3y_bb = I_NAI_Ix5y_D2y_bb+ABY*I_NAI_Hx4y_D2y_bb;
  Double I_NAI_Hx3yz_F3y_bb = I_NAI_Ix4yz_D2y_bb+ABY*I_NAI_Hx3yz_D2y_bb;
  Double I_NAI_Hx2y2z_F3y_bb = I_NAI_Ix3y2z_D2y_bb+ABY*I_NAI_Hx2y2z_D2y_bb;
  Double I_NAI_Hxy3z_F3y_bb = I_NAI_Ix2y3z_D2y_bb+ABY*I_NAI_Hxy3z_D2y_bb;
  Double I_NAI_Hx4z_F3y_bb = I_NAI_Ixy4z_D2y_bb+ABY*I_NAI_Hx4z_D2y_bb;
  Double I_NAI_H5y_F3y_bb = I_NAI_I6y_D2y_bb+ABY*I_NAI_H5y_D2y_bb;
  Double I_NAI_H4yz_F3y_bb = I_NAI_I5yz_D2y_bb+ABY*I_NAI_H4yz_D2y_bb;
  Double I_NAI_H3y2z_F3y_bb = I_NAI_I4y2z_D2y_bb+ABY*I_NAI_H3y2z_D2y_bb;
  Double I_NAI_H2y3z_F3y_bb = I_NAI_I3y3z_D2y_bb+ABY*I_NAI_H2y3z_D2y_bb;
  Double I_NAI_Hy4z_F3y_bb = I_NAI_I2y4z_D2y_bb+ABY*I_NAI_Hy4z_D2y_bb;
  Double I_NAI_H5z_F3y_bb = I_NAI_Iy5z_D2y_bb+ABY*I_NAI_H5z_D2y_bb;
  Double I_NAI_H5x_F2yz_bb = I_NAI_I5xz_D2y_bb+ABZ*I_NAI_H5x_D2y_bb;
  Double I_NAI_H4xy_F2yz_bb = I_NAI_I4xyz_D2y_bb+ABZ*I_NAI_H4xy_D2y_bb;
  Double I_NAI_H4xz_F2yz_bb = I_NAI_I4x2z_D2y_bb+ABZ*I_NAI_H4xz_D2y_bb;
  Double I_NAI_H3x2y_F2yz_bb = I_NAI_I3x2yz_D2y_bb+ABZ*I_NAI_H3x2y_D2y_bb;
  Double I_NAI_H3xyz_F2yz_bb = I_NAI_I3xy2z_D2y_bb+ABZ*I_NAI_H3xyz_D2y_bb;
  Double I_NAI_H3x2z_F2yz_bb = I_NAI_I3x3z_D2y_bb+ABZ*I_NAI_H3x2z_D2y_bb;
  Double I_NAI_H2x3y_F2yz_bb = I_NAI_I2x3yz_D2y_bb+ABZ*I_NAI_H2x3y_D2y_bb;
  Double I_NAI_H2x2yz_F2yz_bb = I_NAI_I2x2y2z_D2y_bb+ABZ*I_NAI_H2x2yz_D2y_bb;
  Double I_NAI_H2xy2z_F2yz_bb = I_NAI_I2xy3z_D2y_bb+ABZ*I_NAI_H2xy2z_D2y_bb;
  Double I_NAI_H2x3z_F2yz_bb = I_NAI_I2x4z_D2y_bb+ABZ*I_NAI_H2x3z_D2y_bb;
  Double I_NAI_Hx4y_F2yz_bb = I_NAI_Ix4yz_D2y_bb+ABZ*I_NAI_Hx4y_D2y_bb;
  Double I_NAI_Hx3yz_F2yz_bb = I_NAI_Ix3y2z_D2y_bb+ABZ*I_NAI_Hx3yz_D2y_bb;
  Double I_NAI_Hx2y2z_F2yz_bb = I_NAI_Ix2y3z_D2y_bb+ABZ*I_NAI_Hx2y2z_D2y_bb;
  Double I_NAI_Hxy3z_F2yz_bb = I_NAI_Ixy4z_D2y_bb+ABZ*I_NAI_Hxy3z_D2y_bb;
  Double I_NAI_Hx4z_F2yz_bb = I_NAI_Ix5z_D2y_bb+ABZ*I_NAI_Hx4z_D2y_bb;
  Double I_NAI_H4yz_F2yz_bb = I_NAI_I4y2z_D2y_bb+ABZ*I_NAI_H4yz_D2y_bb;
  Double I_NAI_H3y2z_F2yz_bb = I_NAI_I3y3z_D2y_bb+ABZ*I_NAI_H3y2z_D2y_bb;
  Double I_NAI_H2y3z_F2yz_bb = I_NAI_I2y4z_D2y_bb+ABZ*I_NAI_H2y3z_D2y_bb;
  Double I_NAI_Hy4z_F2yz_bb = I_NAI_Iy5z_D2y_bb+ABZ*I_NAI_Hy4z_D2y_bb;
  Double I_NAI_H5z_F2yz_bb = I_NAI_I6z_D2y_bb+ABZ*I_NAI_H5z_D2y_bb;
  Double I_NAI_H5x_F3z_bb = I_NAI_I5xz_D2z_bb+ABZ*I_NAI_H5x_D2z_bb;
  Double I_NAI_H4xy_F3z_bb = I_NAI_I4xyz_D2z_bb+ABZ*I_NAI_H4xy_D2z_bb;
  Double I_NAI_H4xz_F3z_bb = I_NAI_I4x2z_D2z_bb+ABZ*I_NAI_H4xz_D2z_bb;
  Double I_NAI_H3x2y_F3z_bb = I_NAI_I3x2yz_D2z_bb+ABZ*I_NAI_H3x2y_D2z_bb;
  Double I_NAI_H3xyz_F3z_bb = I_NAI_I3xy2z_D2z_bb+ABZ*I_NAI_H3xyz_D2z_bb;
  Double I_NAI_H3x2z_F3z_bb = I_NAI_I3x3z_D2z_bb+ABZ*I_NAI_H3x2z_D2z_bb;
  Double I_NAI_H2x3y_F3z_bb = I_NAI_I2x3yz_D2z_bb+ABZ*I_NAI_H2x3y_D2z_bb;
  Double I_NAI_H2x2yz_F3z_bb = I_NAI_I2x2y2z_D2z_bb+ABZ*I_NAI_H2x2yz_D2z_bb;
  Double I_NAI_H2xy2z_F3z_bb = I_NAI_I2xy3z_D2z_bb+ABZ*I_NAI_H2xy2z_D2z_bb;
  Double I_NAI_H2x3z_F3z_bb = I_NAI_I2x4z_D2z_bb+ABZ*I_NAI_H2x3z_D2z_bb;
  Double I_NAI_Hx4y_F3z_bb = I_NAI_Ix4yz_D2z_bb+ABZ*I_NAI_Hx4y_D2z_bb;
  Double I_NAI_Hx3yz_F3z_bb = I_NAI_Ix3y2z_D2z_bb+ABZ*I_NAI_Hx3yz_D2z_bb;
  Double I_NAI_Hx2y2z_F3z_bb = I_NAI_Ix2y3z_D2z_bb+ABZ*I_NAI_Hx2y2z_D2z_bb;
  Double I_NAI_Hxy3z_F3z_bb = I_NAI_Ixy4z_D2z_bb+ABZ*I_NAI_Hxy3z_D2z_bb;
  Double I_NAI_Hx4z_F3z_bb = I_NAI_Ix5z_D2z_bb+ABZ*I_NAI_Hx4z_D2z_bb;
  Double I_NAI_H5y_F3z_bb = I_NAI_I5yz_D2z_bb+ABZ*I_NAI_H5y_D2z_bb;
  Double I_NAI_H4yz_F3z_bb = I_NAI_I4y2z_D2z_bb+ABZ*I_NAI_H4yz_D2z_bb;
  Double I_NAI_H3y2z_F3z_bb = I_NAI_I3y3z_D2z_bb+ABZ*I_NAI_H3y2z_D2z_bb;
  Double I_NAI_H2y3z_F3z_bb = I_NAI_I2y4z_D2z_bb+ABZ*I_NAI_H2y3z_D2z_bb;
  Double I_NAI_Hy4z_F3z_bb = I_NAI_Iy5z_D2z_bb+ABZ*I_NAI_Hy4z_D2z_bb;
  Double I_NAI_H5z_F3z_bb = I_NAI_I6z_D2z_bb+ABZ*I_NAI_H5z_D2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_G_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_G_F_bb
   ************************************************************/
  Double I_NAI_G4x_G4x_bb = I_NAI_H5x_F3x_bb+ABX*I_NAI_G4x_F3x_bb;
  Double I_NAI_G3xy_G4x_bb = I_NAI_H4xy_F3x_bb+ABX*I_NAI_G3xy_F3x_bb;
  Double I_NAI_G3xz_G4x_bb = I_NAI_H4xz_F3x_bb+ABX*I_NAI_G3xz_F3x_bb;
  Double I_NAI_G2x2y_G4x_bb = I_NAI_H3x2y_F3x_bb+ABX*I_NAI_G2x2y_F3x_bb;
  Double I_NAI_G2xyz_G4x_bb = I_NAI_H3xyz_F3x_bb+ABX*I_NAI_G2xyz_F3x_bb;
  Double I_NAI_G2x2z_G4x_bb = I_NAI_H3x2z_F3x_bb+ABX*I_NAI_G2x2z_F3x_bb;
  Double I_NAI_Gx3y_G4x_bb = I_NAI_H2x3y_F3x_bb+ABX*I_NAI_Gx3y_F3x_bb;
  Double I_NAI_Gx2yz_G4x_bb = I_NAI_H2x2yz_F3x_bb+ABX*I_NAI_Gx2yz_F3x_bb;
  Double I_NAI_Gxy2z_G4x_bb = I_NAI_H2xy2z_F3x_bb+ABX*I_NAI_Gxy2z_F3x_bb;
  Double I_NAI_Gx3z_G4x_bb = I_NAI_H2x3z_F3x_bb+ABX*I_NAI_Gx3z_F3x_bb;
  Double I_NAI_G4y_G4x_bb = I_NAI_Hx4y_F3x_bb+ABX*I_NAI_G4y_F3x_bb;
  Double I_NAI_G3yz_G4x_bb = I_NAI_Hx3yz_F3x_bb+ABX*I_NAI_G3yz_F3x_bb;
  Double I_NAI_G2y2z_G4x_bb = I_NAI_Hx2y2z_F3x_bb+ABX*I_NAI_G2y2z_F3x_bb;
  Double I_NAI_Gy3z_G4x_bb = I_NAI_Hxy3z_F3x_bb+ABX*I_NAI_Gy3z_F3x_bb;
  Double I_NAI_G4z_G4x_bb = I_NAI_Hx4z_F3x_bb+ABX*I_NAI_G4z_F3x_bb;
  Double I_NAI_G4x_G3xy_bb = I_NAI_H4xy_F3x_bb+ABY*I_NAI_G4x_F3x_bb;
  Double I_NAI_G3xy_G3xy_bb = I_NAI_H3x2y_F3x_bb+ABY*I_NAI_G3xy_F3x_bb;
  Double I_NAI_G3xz_G3xy_bb = I_NAI_H3xyz_F3x_bb+ABY*I_NAI_G3xz_F3x_bb;
  Double I_NAI_G2x2y_G3xy_bb = I_NAI_H2x3y_F3x_bb+ABY*I_NAI_G2x2y_F3x_bb;
  Double I_NAI_G2xyz_G3xy_bb = I_NAI_H2x2yz_F3x_bb+ABY*I_NAI_G2xyz_F3x_bb;
  Double I_NAI_G2x2z_G3xy_bb = I_NAI_H2xy2z_F3x_bb+ABY*I_NAI_G2x2z_F3x_bb;
  Double I_NAI_Gx3y_G3xy_bb = I_NAI_Hx4y_F3x_bb+ABY*I_NAI_Gx3y_F3x_bb;
  Double I_NAI_Gx2yz_G3xy_bb = I_NAI_Hx3yz_F3x_bb+ABY*I_NAI_Gx2yz_F3x_bb;
  Double I_NAI_Gxy2z_G3xy_bb = I_NAI_Hx2y2z_F3x_bb+ABY*I_NAI_Gxy2z_F3x_bb;
  Double I_NAI_Gx3z_G3xy_bb = I_NAI_Hxy3z_F3x_bb+ABY*I_NAI_Gx3z_F3x_bb;
  Double I_NAI_G4y_G3xy_bb = I_NAI_H5y_F3x_bb+ABY*I_NAI_G4y_F3x_bb;
  Double I_NAI_G3yz_G3xy_bb = I_NAI_H4yz_F3x_bb+ABY*I_NAI_G3yz_F3x_bb;
  Double I_NAI_G2y2z_G3xy_bb = I_NAI_H3y2z_F3x_bb+ABY*I_NAI_G2y2z_F3x_bb;
  Double I_NAI_Gy3z_G3xy_bb = I_NAI_H2y3z_F3x_bb+ABY*I_NAI_Gy3z_F3x_bb;
  Double I_NAI_G4z_G3xy_bb = I_NAI_Hy4z_F3x_bb+ABY*I_NAI_G4z_F3x_bb;
  Double I_NAI_G4x_G3xz_bb = I_NAI_H4xz_F3x_bb+ABZ*I_NAI_G4x_F3x_bb;
  Double I_NAI_G3xy_G3xz_bb = I_NAI_H3xyz_F3x_bb+ABZ*I_NAI_G3xy_F3x_bb;
  Double I_NAI_G3xz_G3xz_bb = I_NAI_H3x2z_F3x_bb+ABZ*I_NAI_G3xz_F3x_bb;
  Double I_NAI_G2x2y_G3xz_bb = I_NAI_H2x2yz_F3x_bb+ABZ*I_NAI_G2x2y_F3x_bb;
  Double I_NAI_G2xyz_G3xz_bb = I_NAI_H2xy2z_F3x_bb+ABZ*I_NAI_G2xyz_F3x_bb;
  Double I_NAI_G2x2z_G3xz_bb = I_NAI_H2x3z_F3x_bb+ABZ*I_NAI_G2x2z_F3x_bb;
  Double I_NAI_Gx3y_G3xz_bb = I_NAI_Hx3yz_F3x_bb+ABZ*I_NAI_Gx3y_F3x_bb;
  Double I_NAI_Gx2yz_G3xz_bb = I_NAI_Hx2y2z_F3x_bb+ABZ*I_NAI_Gx2yz_F3x_bb;
  Double I_NAI_Gxy2z_G3xz_bb = I_NAI_Hxy3z_F3x_bb+ABZ*I_NAI_Gxy2z_F3x_bb;
  Double I_NAI_Gx3z_G3xz_bb = I_NAI_Hx4z_F3x_bb+ABZ*I_NAI_Gx3z_F3x_bb;
  Double I_NAI_G4y_G3xz_bb = I_NAI_H4yz_F3x_bb+ABZ*I_NAI_G4y_F3x_bb;
  Double I_NAI_G3yz_G3xz_bb = I_NAI_H3y2z_F3x_bb+ABZ*I_NAI_G3yz_F3x_bb;
  Double I_NAI_G2y2z_G3xz_bb = I_NAI_H2y3z_F3x_bb+ABZ*I_NAI_G2y2z_F3x_bb;
  Double I_NAI_Gy3z_G3xz_bb = I_NAI_Hy4z_F3x_bb+ABZ*I_NAI_Gy3z_F3x_bb;
  Double I_NAI_G4z_G3xz_bb = I_NAI_H5z_F3x_bb+ABZ*I_NAI_G4z_F3x_bb;
  Double I_NAI_G4x_G2x2y_bb = I_NAI_H4xy_F2xy_bb+ABY*I_NAI_G4x_F2xy_bb;
  Double I_NAI_G3xy_G2x2y_bb = I_NAI_H3x2y_F2xy_bb+ABY*I_NAI_G3xy_F2xy_bb;
  Double I_NAI_G3xz_G2x2y_bb = I_NAI_H3xyz_F2xy_bb+ABY*I_NAI_G3xz_F2xy_bb;
  Double I_NAI_G2x2y_G2x2y_bb = I_NAI_H2x3y_F2xy_bb+ABY*I_NAI_G2x2y_F2xy_bb;
  Double I_NAI_G2xyz_G2x2y_bb = I_NAI_H2x2yz_F2xy_bb+ABY*I_NAI_G2xyz_F2xy_bb;
  Double I_NAI_G2x2z_G2x2y_bb = I_NAI_H2xy2z_F2xy_bb+ABY*I_NAI_G2x2z_F2xy_bb;
  Double I_NAI_Gx3y_G2x2y_bb = I_NAI_Hx4y_F2xy_bb+ABY*I_NAI_Gx3y_F2xy_bb;
  Double I_NAI_Gx2yz_G2x2y_bb = I_NAI_Hx3yz_F2xy_bb+ABY*I_NAI_Gx2yz_F2xy_bb;
  Double I_NAI_Gxy2z_G2x2y_bb = I_NAI_Hx2y2z_F2xy_bb+ABY*I_NAI_Gxy2z_F2xy_bb;
  Double I_NAI_Gx3z_G2x2y_bb = I_NAI_Hxy3z_F2xy_bb+ABY*I_NAI_Gx3z_F2xy_bb;
  Double I_NAI_G4y_G2x2y_bb = I_NAI_H5y_F2xy_bb+ABY*I_NAI_G4y_F2xy_bb;
  Double I_NAI_G3yz_G2x2y_bb = I_NAI_H4yz_F2xy_bb+ABY*I_NAI_G3yz_F2xy_bb;
  Double I_NAI_G2y2z_G2x2y_bb = I_NAI_H3y2z_F2xy_bb+ABY*I_NAI_G2y2z_F2xy_bb;
  Double I_NAI_Gy3z_G2x2y_bb = I_NAI_H2y3z_F2xy_bb+ABY*I_NAI_Gy3z_F2xy_bb;
  Double I_NAI_G4z_G2x2y_bb = I_NAI_Hy4z_F2xy_bb+ABY*I_NAI_G4z_F2xy_bb;
  Double I_NAI_G4x_G2x2z_bb = I_NAI_H4xz_F2xz_bb+ABZ*I_NAI_G4x_F2xz_bb;
  Double I_NAI_G3xy_G2x2z_bb = I_NAI_H3xyz_F2xz_bb+ABZ*I_NAI_G3xy_F2xz_bb;
  Double I_NAI_G3xz_G2x2z_bb = I_NAI_H3x2z_F2xz_bb+ABZ*I_NAI_G3xz_F2xz_bb;
  Double I_NAI_G2x2y_G2x2z_bb = I_NAI_H2x2yz_F2xz_bb+ABZ*I_NAI_G2x2y_F2xz_bb;
  Double I_NAI_G2xyz_G2x2z_bb = I_NAI_H2xy2z_F2xz_bb+ABZ*I_NAI_G2xyz_F2xz_bb;
  Double I_NAI_G2x2z_G2x2z_bb = I_NAI_H2x3z_F2xz_bb+ABZ*I_NAI_G2x2z_F2xz_bb;
  Double I_NAI_Gx3y_G2x2z_bb = I_NAI_Hx3yz_F2xz_bb+ABZ*I_NAI_Gx3y_F2xz_bb;
  Double I_NAI_Gx2yz_G2x2z_bb = I_NAI_Hx2y2z_F2xz_bb+ABZ*I_NAI_Gx2yz_F2xz_bb;
  Double I_NAI_Gxy2z_G2x2z_bb = I_NAI_Hxy3z_F2xz_bb+ABZ*I_NAI_Gxy2z_F2xz_bb;
  Double I_NAI_Gx3z_G2x2z_bb = I_NAI_Hx4z_F2xz_bb+ABZ*I_NAI_Gx3z_F2xz_bb;
  Double I_NAI_G4y_G2x2z_bb = I_NAI_H4yz_F2xz_bb+ABZ*I_NAI_G4y_F2xz_bb;
  Double I_NAI_G3yz_G2x2z_bb = I_NAI_H3y2z_F2xz_bb+ABZ*I_NAI_G3yz_F2xz_bb;
  Double I_NAI_G2y2z_G2x2z_bb = I_NAI_H2y3z_F2xz_bb+ABZ*I_NAI_G2y2z_F2xz_bb;
  Double I_NAI_Gy3z_G2x2z_bb = I_NAI_Hy4z_F2xz_bb+ABZ*I_NAI_Gy3z_F2xz_bb;
  Double I_NAI_G4z_G2x2z_bb = I_NAI_H5z_F2xz_bb+ABZ*I_NAI_G4z_F2xz_bb;
  Double I_NAI_G4x_Gx3y_bb = I_NAI_H5x_F3y_bb+ABX*I_NAI_G4x_F3y_bb;
  Double I_NAI_G3xy_Gx3y_bb = I_NAI_H4xy_F3y_bb+ABX*I_NAI_G3xy_F3y_bb;
  Double I_NAI_G3xz_Gx3y_bb = I_NAI_H4xz_F3y_bb+ABX*I_NAI_G3xz_F3y_bb;
  Double I_NAI_G2x2y_Gx3y_bb = I_NAI_H3x2y_F3y_bb+ABX*I_NAI_G2x2y_F3y_bb;
  Double I_NAI_G2xyz_Gx3y_bb = I_NAI_H3xyz_F3y_bb+ABX*I_NAI_G2xyz_F3y_bb;
  Double I_NAI_G2x2z_Gx3y_bb = I_NAI_H3x2z_F3y_bb+ABX*I_NAI_G2x2z_F3y_bb;
  Double I_NAI_Gx3y_Gx3y_bb = I_NAI_H2x3y_F3y_bb+ABX*I_NAI_Gx3y_F3y_bb;
  Double I_NAI_Gx2yz_Gx3y_bb = I_NAI_H2x2yz_F3y_bb+ABX*I_NAI_Gx2yz_F3y_bb;
  Double I_NAI_Gxy2z_Gx3y_bb = I_NAI_H2xy2z_F3y_bb+ABX*I_NAI_Gxy2z_F3y_bb;
  Double I_NAI_Gx3z_Gx3y_bb = I_NAI_H2x3z_F3y_bb+ABX*I_NAI_Gx3z_F3y_bb;
  Double I_NAI_G4y_Gx3y_bb = I_NAI_Hx4y_F3y_bb+ABX*I_NAI_G4y_F3y_bb;
  Double I_NAI_G3yz_Gx3y_bb = I_NAI_Hx3yz_F3y_bb+ABX*I_NAI_G3yz_F3y_bb;
  Double I_NAI_G2y2z_Gx3y_bb = I_NAI_Hx2y2z_F3y_bb+ABX*I_NAI_G2y2z_F3y_bb;
  Double I_NAI_Gy3z_Gx3y_bb = I_NAI_Hxy3z_F3y_bb+ABX*I_NAI_Gy3z_F3y_bb;
  Double I_NAI_G4z_Gx3y_bb = I_NAI_Hx4z_F3y_bb+ABX*I_NAI_G4z_F3y_bb;
  Double I_NAI_G4x_Gx3z_bb = I_NAI_H5x_F3z_bb+ABX*I_NAI_G4x_F3z_bb;
  Double I_NAI_G3xy_Gx3z_bb = I_NAI_H4xy_F3z_bb+ABX*I_NAI_G3xy_F3z_bb;
  Double I_NAI_G3xz_Gx3z_bb = I_NAI_H4xz_F3z_bb+ABX*I_NAI_G3xz_F3z_bb;
  Double I_NAI_G2x2y_Gx3z_bb = I_NAI_H3x2y_F3z_bb+ABX*I_NAI_G2x2y_F3z_bb;
  Double I_NAI_G2xyz_Gx3z_bb = I_NAI_H3xyz_F3z_bb+ABX*I_NAI_G2xyz_F3z_bb;
  Double I_NAI_G2x2z_Gx3z_bb = I_NAI_H3x2z_F3z_bb+ABX*I_NAI_G2x2z_F3z_bb;
  Double I_NAI_Gx3y_Gx3z_bb = I_NAI_H2x3y_F3z_bb+ABX*I_NAI_Gx3y_F3z_bb;
  Double I_NAI_Gx2yz_Gx3z_bb = I_NAI_H2x2yz_F3z_bb+ABX*I_NAI_Gx2yz_F3z_bb;
  Double I_NAI_Gxy2z_Gx3z_bb = I_NAI_H2xy2z_F3z_bb+ABX*I_NAI_Gxy2z_F3z_bb;
  Double I_NAI_Gx3z_Gx3z_bb = I_NAI_H2x3z_F3z_bb+ABX*I_NAI_Gx3z_F3z_bb;
  Double I_NAI_G4y_Gx3z_bb = I_NAI_Hx4y_F3z_bb+ABX*I_NAI_G4y_F3z_bb;
  Double I_NAI_G3yz_Gx3z_bb = I_NAI_Hx3yz_F3z_bb+ABX*I_NAI_G3yz_F3z_bb;
  Double I_NAI_G2y2z_Gx3z_bb = I_NAI_Hx2y2z_F3z_bb+ABX*I_NAI_G2y2z_F3z_bb;
  Double I_NAI_Gy3z_Gx3z_bb = I_NAI_Hxy3z_F3z_bb+ABX*I_NAI_Gy3z_F3z_bb;
  Double I_NAI_G4z_Gx3z_bb = I_NAI_Hx4z_F3z_bb+ABX*I_NAI_G4z_F3z_bb;
  Double I_NAI_G4x_G4y_bb = I_NAI_H4xy_F3y_bb+ABY*I_NAI_G4x_F3y_bb;
  Double I_NAI_G3xy_G4y_bb = I_NAI_H3x2y_F3y_bb+ABY*I_NAI_G3xy_F3y_bb;
  Double I_NAI_G3xz_G4y_bb = I_NAI_H3xyz_F3y_bb+ABY*I_NAI_G3xz_F3y_bb;
  Double I_NAI_G2x2y_G4y_bb = I_NAI_H2x3y_F3y_bb+ABY*I_NAI_G2x2y_F3y_bb;
  Double I_NAI_G2xyz_G4y_bb = I_NAI_H2x2yz_F3y_bb+ABY*I_NAI_G2xyz_F3y_bb;
  Double I_NAI_G2x2z_G4y_bb = I_NAI_H2xy2z_F3y_bb+ABY*I_NAI_G2x2z_F3y_bb;
  Double I_NAI_Gx3y_G4y_bb = I_NAI_Hx4y_F3y_bb+ABY*I_NAI_Gx3y_F3y_bb;
  Double I_NAI_Gx2yz_G4y_bb = I_NAI_Hx3yz_F3y_bb+ABY*I_NAI_Gx2yz_F3y_bb;
  Double I_NAI_Gxy2z_G4y_bb = I_NAI_Hx2y2z_F3y_bb+ABY*I_NAI_Gxy2z_F3y_bb;
  Double I_NAI_Gx3z_G4y_bb = I_NAI_Hxy3z_F3y_bb+ABY*I_NAI_Gx3z_F3y_bb;
  Double I_NAI_G4y_G4y_bb = I_NAI_H5y_F3y_bb+ABY*I_NAI_G4y_F3y_bb;
  Double I_NAI_G3yz_G4y_bb = I_NAI_H4yz_F3y_bb+ABY*I_NAI_G3yz_F3y_bb;
  Double I_NAI_G2y2z_G4y_bb = I_NAI_H3y2z_F3y_bb+ABY*I_NAI_G2y2z_F3y_bb;
  Double I_NAI_Gy3z_G4y_bb = I_NAI_H2y3z_F3y_bb+ABY*I_NAI_Gy3z_F3y_bb;
  Double I_NAI_G4z_G4y_bb = I_NAI_Hy4z_F3y_bb+ABY*I_NAI_G4z_F3y_bb;
  Double I_NAI_G4x_G3yz_bb = I_NAI_H4xz_F3y_bb+ABZ*I_NAI_G4x_F3y_bb;
  Double I_NAI_G3xy_G3yz_bb = I_NAI_H3xyz_F3y_bb+ABZ*I_NAI_G3xy_F3y_bb;
  Double I_NAI_G3xz_G3yz_bb = I_NAI_H3x2z_F3y_bb+ABZ*I_NAI_G3xz_F3y_bb;
  Double I_NAI_G2x2y_G3yz_bb = I_NAI_H2x2yz_F3y_bb+ABZ*I_NAI_G2x2y_F3y_bb;
  Double I_NAI_G2xyz_G3yz_bb = I_NAI_H2xy2z_F3y_bb+ABZ*I_NAI_G2xyz_F3y_bb;
  Double I_NAI_G2x2z_G3yz_bb = I_NAI_H2x3z_F3y_bb+ABZ*I_NAI_G2x2z_F3y_bb;
  Double I_NAI_Gx3y_G3yz_bb = I_NAI_Hx3yz_F3y_bb+ABZ*I_NAI_Gx3y_F3y_bb;
  Double I_NAI_Gx2yz_G3yz_bb = I_NAI_Hx2y2z_F3y_bb+ABZ*I_NAI_Gx2yz_F3y_bb;
  Double I_NAI_Gxy2z_G3yz_bb = I_NAI_Hxy3z_F3y_bb+ABZ*I_NAI_Gxy2z_F3y_bb;
  Double I_NAI_Gx3z_G3yz_bb = I_NAI_Hx4z_F3y_bb+ABZ*I_NAI_Gx3z_F3y_bb;
  Double I_NAI_G4y_G3yz_bb = I_NAI_H4yz_F3y_bb+ABZ*I_NAI_G4y_F3y_bb;
  Double I_NAI_G3yz_G3yz_bb = I_NAI_H3y2z_F3y_bb+ABZ*I_NAI_G3yz_F3y_bb;
  Double I_NAI_G2y2z_G3yz_bb = I_NAI_H2y3z_F3y_bb+ABZ*I_NAI_G2y2z_F3y_bb;
  Double I_NAI_Gy3z_G3yz_bb = I_NAI_Hy4z_F3y_bb+ABZ*I_NAI_Gy3z_F3y_bb;
  Double I_NAI_G4z_G3yz_bb = I_NAI_H5z_F3y_bb+ABZ*I_NAI_G4z_F3y_bb;
  Double I_NAI_G4x_G2y2z_bb = I_NAI_H4xz_F2yz_bb+ABZ*I_NAI_G4x_F2yz_bb;
  Double I_NAI_G3xy_G2y2z_bb = I_NAI_H3xyz_F2yz_bb+ABZ*I_NAI_G3xy_F2yz_bb;
  Double I_NAI_G3xz_G2y2z_bb = I_NAI_H3x2z_F2yz_bb+ABZ*I_NAI_G3xz_F2yz_bb;
  Double I_NAI_G2x2y_G2y2z_bb = I_NAI_H2x2yz_F2yz_bb+ABZ*I_NAI_G2x2y_F2yz_bb;
  Double I_NAI_G2xyz_G2y2z_bb = I_NAI_H2xy2z_F2yz_bb+ABZ*I_NAI_G2xyz_F2yz_bb;
  Double I_NAI_G2x2z_G2y2z_bb = I_NAI_H2x3z_F2yz_bb+ABZ*I_NAI_G2x2z_F2yz_bb;
  Double I_NAI_Gx3y_G2y2z_bb = I_NAI_Hx3yz_F2yz_bb+ABZ*I_NAI_Gx3y_F2yz_bb;
  Double I_NAI_Gx2yz_G2y2z_bb = I_NAI_Hx2y2z_F2yz_bb+ABZ*I_NAI_Gx2yz_F2yz_bb;
  Double I_NAI_Gxy2z_G2y2z_bb = I_NAI_Hxy3z_F2yz_bb+ABZ*I_NAI_Gxy2z_F2yz_bb;
  Double I_NAI_Gx3z_G2y2z_bb = I_NAI_Hx4z_F2yz_bb+ABZ*I_NAI_Gx3z_F2yz_bb;
  Double I_NAI_G4y_G2y2z_bb = I_NAI_H4yz_F2yz_bb+ABZ*I_NAI_G4y_F2yz_bb;
  Double I_NAI_G3yz_G2y2z_bb = I_NAI_H3y2z_F2yz_bb+ABZ*I_NAI_G3yz_F2yz_bb;
  Double I_NAI_G2y2z_G2y2z_bb = I_NAI_H2y3z_F2yz_bb+ABZ*I_NAI_G2y2z_F2yz_bb;
  Double I_NAI_Gy3z_G2y2z_bb = I_NAI_Hy4z_F2yz_bb+ABZ*I_NAI_Gy3z_F2yz_bb;
  Double I_NAI_G4z_G2y2z_bb = I_NAI_H5z_F2yz_bb+ABZ*I_NAI_G4z_F2yz_bb;
  Double I_NAI_G4x_Gy3z_bb = I_NAI_H4xy_F3z_bb+ABY*I_NAI_G4x_F3z_bb;
  Double I_NAI_G3xy_Gy3z_bb = I_NAI_H3x2y_F3z_bb+ABY*I_NAI_G3xy_F3z_bb;
  Double I_NAI_G3xz_Gy3z_bb = I_NAI_H3xyz_F3z_bb+ABY*I_NAI_G3xz_F3z_bb;
  Double I_NAI_G2x2y_Gy3z_bb = I_NAI_H2x3y_F3z_bb+ABY*I_NAI_G2x2y_F3z_bb;
  Double I_NAI_G2xyz_Gy3z_bb = I_NAI_H2x2yz_F3z_bb+ABY*I_NAI_G2xyz_F3z_bb;
  Double I_NAI_G2x2z_Gy3z_bb = I_NAI_H2xy2z_F3z_bb+ABY*I_NAI_G2x2z_F3z_bb;
  Double I_NAI_Gx3y_Gy3z_bb = I_NAI_Hx4y_F3z_bb+ABY*I_NAI_Gx3y_F3z_bb;
  Double I_NAI_Gx2yz_Gy3z_bb = I_NAI_Hx3yz_F3z_bb+ABY*I_NAI_Gx2yz_F3z_bb;
  Double I_NAI_Gxy2z_Gy3z_bb = I_NAI_Hx2y2z_F3z_bb+ABY*I_NAI_Gxy2z_F3z_bb;
  Double I_NAI_Gx3z_Gy3z_bb = I_NAI_Hxy3z_F3z_bb+ABY*I_NAI_Gx3z_F3z_bb;
  Double I_NAI_G4y_Gy3z_bb = I_NAI_H5y_F3z_bb+ABY*I_NAI_G4y_F3z_bb;
  Double I_NAI_G3yz_Gy3z_bb = I_NAI_H4yz_F3z_bb+ABY*I_NAI_G3yz_F3z_bb;
  Double I_NAI_G2y2z_Gy3z_bb = I_NAI_H3y2z_F3z_bb+ABY*I_NAI_G2y2z_F3z_bb;
  Double I_NAI_Gy3z_Gy3z_bb = I_NAI_H2y3z_F3z_bb+ABY*I_NAI_Gy3z_F3z_bb;
  Double I_NAI_G4z_Gy3z_bb = I_NAI_Hy4z_F3z_bb+ABY*I_NAI_G4z_F3z_bb;
  Double I_NAI_G4x_G4z_bb = I_NAI_H4xz_F3z_bb+ABZ*I_NAI_G4x_F3z_bb;
  Double I_NAI_G3xy_G4z_bb = I_NAI_H3xyz_F3z_bb+ABZ*I_NAI_G3xy_F3z_bb;
  Double I_NAI_G3xz_G4z_bb = I_NAI_H3x2z_F3z_bb+ABZ*I_NAI_G3xz_F3z_bb;
  Double I_NAI_G2x2y_G4z_bb = I_NAI_H2x2yz_F3z_bb+ABZ*I_NAI_G2x2y_F3z_bb;
  Double I_NAI_G2xyz_G4z_bb = I_NAI_H2xy2z_F3z_bb+ABZ*I_NAI_G2xyz_F3z_bb;
  Double I_NAI_G2x2z_G4z_bb = I_NAI_H2x3z_F3z_bb+ABZ*I_NAI_G2x2z_F3z_bb;
  Double I_NAI_Gx3y_G4z_bb = I_NAI_Hx3yz_F3z_bb+ABZ*I_NAI_Gx3y_F3z_bb;
  Double I_NAI_Gx2yz_G4z_bb = I_NAI_Hx2y2z_F3z_bb+ABZ*I_NAI_Gx2yz_F3z_bb;
  Double I_NAI_Gxy2z_G4z_bb = I_NAI_Hxy3z_F3z_bb+ABZ*I_NAI_Gxy2z_F3z_bb;
  Double I_NAI_Gx3z_G4z_bb = I_NAI_Hx4z_F3z_bb+ABZ*I_NAI_Gx3z_F3z_bb;
  Double I_NAI_G4y_G4z_bb = I_NAI_H4yz_F3z_bb+ABZ*I_NAI_G4y_F3z_bb;
  Double I_NAI_G3yz_G4z_bb = I_NAI_H3y2z_F3z_bb+ABZ*I_NAI_G3yz_F3z_bb;
  Double I_NAI_G2y2z_G4z_bb = I_NAI_H2y3z_F3z_bb+ABZ*I_NAI_G2y2z_F3z_bb;
  Double I_NAI_Gy3z_G4z_bb = I_NAI_Hy4z_F3z_bb+ABZ*I_NAI_Gy3z_F3z_bb;
  Double I_NAI_G4z_G4z_bb = I_NAI_H5z_F3z_bb+ABZ*I_NAI_G4z_F3z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_L_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 40 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_S_bb
   * RHS shell quartet name: SQ_NAI_L_S_bb
   ************************************************************/
  Double I_NAI_L8x_Px_bb = I_NAI_M9x_S_bb+ABX*I_NAI_L8x_S_bb;
  Double I_NAI_L7xy_Px_bb = I_NAI_M8xy_S_bb+ABX*I_NAI_L7xy_S_bb;
  Double I_NAI_L7xz_Px_bb = I_NAI_M8xz_S_bb+ABX*I_NAI_L7xz_S_bb;
  Double I_NAI_L6x2y_Px_bb = I_NAI_M7x2y_S_bb+ABX*I_NAI_L6x2y_S_bb;
  Double I_NAI_L6xyz_Px_bb = I_NAI_M7xyz_S_bb+ABX*I_NAI_L6xyz_S_bb;
  Double I_NAI_L6x2z_Px_bb = I_NAI_M7x2z_S_bb+ABX*I_NAI_L6x2z_S_bb;
  Double I_NAI_L5x3y_Px_bb = I_NAI_M6x3y_S_bb+ABX*I_NAI_L5x3y_S_bb;
  Double I_NAI_L5x2yz_Px_bb = I_NAI_M6x2yz_S_bb+ABX*I_NAI_L5x2yz_S_bb;
  Double I_NAI_L5xy2z_Px_bb = I_NAI_M6xy2z_S_bb+ABX*I_NAI_L5xy2z_S_bb;
  Double I_NAI_L5x3z_Px_bb = I_NAI_M6x3z_S_bb+ABX*I_NAI_L5x3z_S_bb;
  Double I_NAI_L4x4y_Px_bb = I_NAI_M5x4y_S_bb+ABX*I_NAI_L4x4y_S_bb;
  Double I_NAI_L4x3yz_Px_bb = I_NAI_M5x3yz_S_bb+ABX*I_NAI_L4x3yz_S_bb;
  Double I_NAI_L4x2y2z_Px_bb = I_NAI_M5x2y2z_S_bb+ABX*I_NAI_L4x2y2z_S_bb;
  Double I_NAI_L4xy3z_Px_bb = I_NAI_M5xy3z_S_bb+ABX*I_NAI_L4xy3z_S_bb;
  Double I_NAI_L4x4z_Px_bb = I_NAI_M5x4z_S_bb+ABX*I_NAI_L4x4z_S_bb;
  Double I_NAI_L3x5y_Px_bb = I_NAI_M4x5y_S_bb+ABX*I_NAI_L3x5y_S_bb;
  Double I_NAI_L3x4yz_Px_bb = I_NAI_M4x4yz_S_bb+ABX*I_NAI_L3x4yz_S_bb;
  Double I_NAI_L3x3y2z_Px_bb = I_NAI_M4x3y2z_S_bb+ABX*I_NAI_L3x3y2z_S_bb;
  Double I_NAI_L3x2y3z_Px_bb = I_NAI_M4x2y3z_S_bb+ABX*I_NAI_L3x2y3z_S_bb;
  Double I_NAI_L3xy4z_Px_bb = I_NAI_M4xy4z_S_bb+ABX*I_NAI_L3xy4z_S_bb;
  Double I_NAI_L3x5z_Px_bb = I_NAI_M4x5z_S_bb+ABX*I_NAI_L3x5z_S_bb;
  Double I_NAI_L2x6y_Px_bb = I_NAI_M3x6y_S_bb+ABX*I_NAI_L2x6y_S_bb;
  Double I_NAI_L2x5yz_Px_bb = I_NAI_M3x5yz_S_bb+ABX*I_NAI_L2x5yz_S_bb;
  Double I_NAI_L2x4y2z_Px_bb = I_NAI_M3x4y2z_S_bb+ABX*I_NAI_L2x4y2z_S_bb;
  Double I_NAI_L2x3y3z_Px_bb = I_NAI_M3x3y3z_S_bb+ABX*I_NAI_L2x3y3z_S_bb;
  Double I_NAI_L2x2y4z_Px_bb = I_NAI_M3x2y4z_S_bb+ABX*I_NAI_L2x2y4z_S_bb;
  Double I_NAI_L2xy5z_Px_bb = I_NAI_M3xy5z_S_bb+ABX*I_NAI_L2xy5z_S_bb;
  Double I_NAI_L2x6z_Px_bb = I_NAI_M3x6z_S_bb+ABX*I_NAI_L2x6z_S_bb;
  Double I_NAI_Lx6yz_Px_bb = I_NAI_M2x6yz_S_bb+ABX*I_NAI_Lx6yz_S_bb;
  Double I_NAI_Lx5y2z_Px_bb = I_NAI_M2x5y2z_S_bb+ABX*I_NAI_Lx5y2z_S_bb;
  Double I_NAI_Lx4y3z_Px_bb = I_NAI_M2x4y3z_S_bb+ABX*I_NAI_Lx4y3z_S_bb;
  Double I_NAI_Lx3y4z_Px_bb = I_NAI_M2x3y4z_S_bb+ABX*I_NAI_Lx3y4z_S_bb;
  Double I_NAI_Lx2y5z_Px_bb = I_NAI_M2x2y5z_S_bb+ABX*I_NAI_Lx2y5z_S_bb;
  Double I_NAI_Lxy6z_Px_bb = I_NAI_M2xy6z_S_bb+ABX*I_NAI_Lxy6z_S_bb;
  Double I_NAI_L6x2y_Py_bb = I_NAI_M6x3y_S_bb+ABY*I_NAI_L6x2y_S_bb;
  Double I_NAI_L5x3y_Py_bb = I_NAI_M5x4y_S_bb+ABY*I_NAI_L5x3y_S_bb;
  Double I_NAI_L5x2yz_Py_bb = I_NAI_M5x3yz_S_bb+ABY*I_NAI_L5x2yz_S_bb;
  Double I_NAI_L5xy2z_Py_bb = I_NAI_M5x2y2z_S_bb+ABY*I_NAI_L5xy2z_S_bb;
  Double I_NAI_L4x4y_Py_bb = I_NAI_M4x5y_S_bb+ABY*I_NAI_L4x4y_S_bb;
  Double I_NAI_L4x3yz_Py_bb = I_NAI_M4x4yz_S_bb+ABY*I_NAI_L4x3yz_S_bb;
  Double I_NAI_L4x2y2z_Py_bb = I_NAI_M4x3y2z_S_bb+ABY*I_NAI_L4x2y2z_S_bb;
  Double I_NAI_L4xy3z_Py_bb = I_NAI_M4x2y3z_S_bb+ABY*I_NAI_L4xy3z_S_bb;
  Double I_NAI_L3x5y_Py_bb = I_NAI_M3x6y_S_bb+ABY*I_NAI_L3x5y_S_bb;
  Double I_NAI_L3x4yz_Py_bb = I_NAI_M3x5yz_S_bb+ABY*I_NAI_L3x4yz_S_bb;
  Double I_NAI_L3x3y2z_Py_bb = I_NAI_M3x4y2z_S_bb+ABY*I_NAI_L3x3y2z_S_bb;
  Double I_NAI_L3x2y3z_Py_bb = I_NAI_M3x3y3z_S_bb+ABY*I_NAI_L3x2y3z_S_bb;
  Double I_NAI_L3xy4z_Py_bb = I_NAI_M3x2y4z_S_bb+ABY*I_NAI_L3xy4z_S_bb;
  Double I_NAI_L2x6y_Py_bb = I_NAI_M2x7y_S_bb+ABY*I_NAI_L2x6y_S_bb;
  Double I_NAI_L2x5yz_Py_bb = I_NAI_M2x6yz_S_bb+ABY*I_NAI_L2x5yz_S_bb;
  Double I_NAI_L2x4y2z_Py_bb = I_NAI_M2x5y2z_S_bb+ABY*I_NAI_L2x4y2z_S_bb;
  Double I_NAI_L2x3y3z_Py_bb = I_NAI_M2x4y3z_S_bb+ABY*I_NAI_L2x3y3z_S_bb;
  Double I_NAI_L2x2y4z_Py_bb = I_NAI_M2x3y4z_S_bb+ABY*I_NAI_L2x2y4z_S_bb;
  Double I_NAI_L2xy5z_Py_bb = I_NAI_M2x2y5z_S_bb+ABY*I_NAI_L2xy5z_S_bb;
  Double I_NAI_Lx7y_Py_bb = I_NAI_Mx8y_S_bb+ABY*I_NAI_Lx7y_S_bb;
  Double I_NAI_Lx6yz_Py_bb = I_NAI_Mx7yz_S_bb+ABY*I_NAI_Lx6yz_S_bb;
  Double I_NAI_Lx5y2z_Py_bb = I_NAI_Mx6y2z_S_bb+ABY*I_NAI_Lx5y2z_S_bb;
  Double I_NAI_Lx4y3z_Py_bb = I_NAI_Mx5y3z_S_bb+ABY*I_NAI_Lx4y3z_S_bb;
  Double I_NAI_Lx3y4z_Py_bb = I_NAI_Mx4y4z_S_bb+ABY*I_NAI_Lx3y4z_S_bb;
  Double I_NAI_Lx2y5z_Py_bb = I_NAI_Mx3y5z_S_bb+ABY*I_NAI_Lx2y5z_S_bb;
  Double I_NAI_Lxy6z_Py_bb = I_NAI_Mx2y6z_S_bb+ABY*I_NAI_Lxy6z_S_bb;
  Double I_NAI_L8y_Py_bb = I_NAI_M9y_S_bb+ABY*I_NAI_L8y_S_bb;
  Double I_NAI_L7yz_Py_bb = I_NAI_M8yz_S_bb+ABY*I_NAI_L7yz_S_bb;
  Double I_NAI_L6y2z_Py_bb = I_NAI_M7y2z_S_bb+ABY*I_NAI_L6y2z_S_bb;
  Double I_NAI_L5y3z_Py_bb = I_NAI_M6y3z_S_bb+ABY*I_NAI_L5y3z_S_bb;
  Double I_NAI_L4y4z_Py_bb = I_NAI_M5y4z_S_bb+ABY*I_NAI_L4y4z_S_bb;
  Double I_NAI_L3y5z_Py_bb = I_NAI_M4y5z_S_bb+ABY*I_NAI_L3y5z_S_bb;
  Double I_NAI_L2y6z_Py_bb = I_NAI_M3y6z_S_bb+ABY*I_NAI_L2y6z_S_bb;
  Double I_NAI_L6x2z_Pz_bb = I_NAI_M6x3z_S_bb+ABZ*I_NAI_L6x2z_S_bb;
  Double I_NAI_L5xy2z_Pz_bb = I_NAI_M5xy3z_S_bb+ABZ*I_NAI_L5xy2z_S_bb;
  Double I_NAI_L5x3z_Pz_bb = I_NAI_M5x4z_S_bb+ABZ*I_NAI_L5x3z_S_bb;
  Double I_NAI_L4x2y2z_Pz_bb = I_NAI_M4x2y3z_S_bb+ABZ*I_NAI_L4x2y2z_S_bb;
  Double I_NAI_L4xy3z_Pz_bb = I_NAI_M4xy4z_S_bb+ABZ*I_NAI_L4xy3z_S_bb;
  Double I_NAI_L4x4z_Pz_bb = I_NAI_M4x5z_S_bb+ABZ*I_NAI_L4x4z_S_bb;
  Double I_NAI_L3x3y2z_Pz_bb = I_NAI_M3x3y3z_S_bb+ABZ*I_NAI_L3x3y2z_S_bb;
  Double I_NAI_L3x2y3z_Pz_bb = I_NAI_M3x2y4z_S_bb+ABZ*I_NAI_L3x2y3z_S_bb;
  Double I_NAI_L3xy4z_Pz_bb = I_NAI_M3xy5z_S_bb+ABZ*I_NAI_L3xy4z_S_bb;
  Double I_NAI_L3x5z_Pz_bb = I_NAI_M3x6z_S_bb+ABZ*I_NAI_L3x5z_S_bb;
  Double I_NAI_L2x4y2z_Pz_bb = I_NAI_M2x4y3z_S_bb+ABZ*I_NAI_L2x4y2z_S_bb;
  Double I_NAI_L2x3y3z_Pz_bb = I_NAI_M2x3y4z_S_bb+ABZ*I_NAI_L2x3y3z_S_bb;
  Double I_NAI_L2x2y4z_Pz_bb = I_NAI_M2x2y5z_S_bb+ABZ*I_NAI_L2x2y4z_S_bb;
  Double I_NAI_L2xy5z_Pz_bb = I_NAI_M2xy6z_S_bb+ABZ*I_NAI_L2xy5z_S_bb;
  Double I_NAI_L2x6z_Pz_bb = I_NAI_M2x7z_S_bb+ABZ*I_NAI_L2x6z_S_bb;
  Double I_NAI_Lx5y2z_Pz_bb = I_NAI_Mx5y3z_S_bb+ABZ*I_NAI_Lx5y2z_S_bb;
  Double I_NAI_Lx4y3z_Pz_bb = I_NAI_Mx4y4z_S_bb+ABZ*I_NAI_Lx4y3z_S_bb;
  Double I_NAI_Lx3y4z_Pz_bb = I_NAI_Mx3y5z_S_bb+ABZ*I_NAI_Lx3y4z_S_bb;
  Double I_NAI_Lx2y5z_Pz_bb = I_NAI_Mx2y6z_S_bb+ABZ*I_NAI_Lx2y5z_S_bb;
  Double I_NAI_Lxy6z_Pz_bb = I_NAI_Mxy7z_S_bb+ABZ*I_NAI_Lxy6z_S_bb;
  Double I_NAI_Lx7z_Pz_bb = I_NAI_Mx8z_S_bb+ABZ*I_NAI_Lx7z_S_bb;
  Double I_NAI_L6y2z_Pz_bb = I_NAI_M6y3z_S_bb+ABZ*I_NAI_L6y2z_S_bb;
  Double I_NAI_L5y3z_Pz_bb = I_NAI_M5y4z_S_bb+ABZ*I_NAI_L5y3z_S_bb;
  Double I_NAI_L4y4z_Pz_bb = I_NAI_M4y5z_S_bb+ABZ*I_NAI_L4y4z_S_bb;
  Double I_NAI_L3y5z_Pz_bb = I_NAI_M3y6z_S_bb+ABZ*I_NAI_L3y5z_S_bb;
  Double I_NAI_L2y6z_Pz_bb = I_NAI_M2y7z_S_bb+ABZ*I_NAI_L2y6z_S_bb;
  Double I_NAI_Ly7z_Pz_bb = I_NAI_My8z_S_bb+ABZ*I_NAI_Ly7z_S_bb;
  Double I_NAI_L8z_Pz_bb = I_NAI_M9z_S_bb+ABZ*I_NAI_L8z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_K_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 121 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_P_bb
   * RHS shell quartet name: SQ_NAI_K_P_bb
   ************************************************************/
  Double I_NAI_K7x_D2x_bb = I_NAI_L8x_Px_bb+ABX*I_NAI_K7x_Px_bb;
  Double I_NAI_K6xy_D2x_bb = I_NAI_L7xy_Px_bb+ABX*I_NAI_K6xy_Px_bb;
  Double I_NAI_K6xz_D2x_bb = I_NAI_L7xz_Px_bb+ABX*I_NAI_K6xz_Px_bb;
  Double I_NAI_K5x2y_D2x_bb = I_NAI_L6x2y_Px_bb+ABX*I_NAI_K5x2y_Px_bb;
  Double I_NAI_K5xyz_D2x_bb = I_NAI_L6xyz_Px_bb+ABX*I_NAI_K5xyz_Px_bb;
  Double I_NAI_K5x2z_D2x_bb = I_NAI_L6x2z_Px_bb+ABX*I_NAI_K5x2z_Px_bb;
  Double I_NAI_K4x3y_D2x_bb = I_NAI_L5x3y_Px_bb+ABX*I_NAI_K4x3y_Px_bb;
  Double I_NAI_K4x2yz_D2x_bb = I_NAI_L5x2yz_Px_bb+ABX*I_NAI_K4x2yz_Px_bb;
  Double I_NAI_K4xy2z_D2x_bb = I_NAI_L5xy2z_Px_bb+ABX*I_NAI_K4xy2z_Px_bb;
  Double I_NAI_K4x3z_D2x_bb = I_NAI_L5x3z_Px_bb+ABX*I_NAI_K4x3z_Px_bb;
  Double I_NAI_K3x4y_D2x_bb = I_NAI_L4x4y_Px_bb+ABX*I_NAI_K3x4y_Px_bb;
  Double I_NAI_K3x3yz_D2x_bb = I_NAI_L4x3yz_Px_bb+ABX*I_NAI_K3x3yz_Px_bb;
  Double I_NAI_K3x2y2z_D2x_bb = I_NAI_L4x2y2z_Px_bb+ABX*I_NAI_K3x2y2z_Px_bb;
  Double I_NAI_K3xy3z_D2x_bb = I_NAI_L4xy3z_Px_bb+ABX*I_NAI_K3xy3z_Px_bb;
  Double I_NAI_K3x4z_D2x_bb = I_NAI_L4x4z_Px_bb+ABX*I_NAI_K3x4z_Px_bb;
  Double I_NAI_K2x5y_D2x_bb = I_NAI_L3x5y_Px_bb+ABX*I_NAI_K2x5y_Px_bb;
  Double I_NAI_K2x4yz_D2x_bb = I_NAI_L3x4yz_Px_bb+ABX*I_NAI_K2x4yz_Px_bb;
  Double I_NAI_K2x3y2z_D2x_bb = I_NAI_L3x3y2z_Px_bb+ABX*I_NAI_K2x3y2z_Px_bb;
  Double I_NAI_K2x2y3z_D2x_bb = I_NAI_L3x2y3z_Px_bb+ABX*I_NAI_K2x2y3z_Px_bb;
  Double I_NAI_K2xy4z_D2x_bb = I_NAI_L3xy4z_Px_bb+ABX*I_NAI_K2xy4z_Px_bb;
  Double I_NAI_K2x5z_D2x_bb = I_NAI_L3x5z_Px_bb+ABX*I_NAI_K2x5z_Px_bb;
  Double I_NAI_Kx6y_D2x_bb = I_NAI_L2x6y_Px_bb+ABX*I_NAI_Kx6y_Px_bb;
  Double I_NAI_Kx5yz_D2x_bb = I_NAI_L2x5yz_Px_bb+ABX*I_NAI_Kx5yz_Px_bb;
  Double I_NAI_Kx4y2z_D2x_bb = I_NAI_L2x4y2z_Px_bb+ABX*I_NAI_Kx4y2z_Px_bb;
  Double I_NAI_Kx3y3z_D2x_bb = I_NAI_L2x3y3z_Px_bb+ABX*I_NAI_Kx3y3z_Px_bb;
  Double I_NAI_Kx2y4z_D2x_bb = I_NAI_L2x2y4z_Px_bb+ABX*I_NAI_Kx2y4z_Px_bb;
  Double I_NAI_Kxy5z_D2x_bb = I_NAI_L2xy5z_Px_bb+ABX*I_NAI_Kxy5z_Px_bb;
  Double I_NAI_Kx6z_D2x_bb = I_NAI_L2x6z_Px_bb+ABX*I_NAI_Kx6z_Px_bb;
  Double I_NAI_K6yz_D2x_bb = I_NAI_Lx6yz_Px_bb+ABX*I_NAI_K6yz_Px_bb;
  Double I_NAI_K5y2z_D2x_bb = I_NAI_Lx5y2z_Px_bb+ABX*I_NAI_K5y2z_Px_bb;
  Double I_NAI_K4y3z_D2x_bb = I_NAI_Lx4y3z_Px_bb+ABX*I_NAI_K4y3z_Px_bb;
  Double I_NAI_K3y4z_D2x_bb = I_NAI_Lx3y4z_Px_bb+ABX*I_NAI_K3y4z_Px_bb;
  Double I_NAI_K2y5z_D2x_bb = I_NAI_Lx2y5z_Px_bb+ABX*I_NAI_K2y5z_Px_bb;
  Double I_NAI_Ky6z_D2x_bb = I_NAI_Lxy6z_Px_bb+ABX*I_NAI_Ky6z_Px_bb;
  Double I_NAI_K6xy_D2y_bb = I_NAI_L6x2y_Py_bb+ABY*I_NAI_K6xy_Py_bb;
  Double I_NAI_K5x2y_D2y_bb = I_NAI_L5x3y_Py_bb+ABY*I_NAI_K5x2y_Py_bb;
  Double I_NAI_K5xyz_D2y_bb = I_NAI_L5x2yz_Py_bb+ABY*I_NAI_K5xyz_Py_bb;
  Double I_NAI_K5x2z_D2y_bb = I_NAI_L5xy2z_Py_bb+ABY*I_NAI_K5x2z_Py_bb;
  Double I_NAI_K4x3y_D2y_bb = I_NAI_L4x4y_Py_bb+ABY*I_NAI_K4x3y_Py_bb;
  Double I_NAI_K4x2yz_D2y_bb = I_NAI_L4x3yz_Py_bb+ABY*I_NAI_K4x2yz_Py_bb;
  Double I_NAI_K4xy2z_D2y_bb = I_NAI_L4x2y2z_Py_bb+ABY*I_NAI_K4xy2z_Py_bb;
  Double I_NAI_K4x3z_D2y_bb = I_NAI_L4xy3z_Py_bb+ABY*I_NAI_K4x3z_Py_bb;
  Double I_NAI_K3x4y_D2y_bb = I_NAI_L3x5y_Py_bb+ABY*I_NAI_K3x4y_Py_bb;
  Double I_NAI_K3x3yz_D2y_bb = I_NAI_L3x4yz_Py_bb+ABY*I_NAI_K3x3yz_Py_bb;
  Double I_NAI_K3x2y2z_D2y_bb = I_NAI_L3x3y2z_Py_bb+ABY*I_NAI_K3x2y2z_Py_bb;
  Double I_NAI_K3xy3z_D2y_bb = I_NAI_L3x2y3z_Py_bb+ABY*I_NAI_K3xy3z_Py_bb;
  Double I_NAI_K3x4z_D2y_bb = I_NAI_L3xy4z_Py_bb+ABY*I_NAI_K3x4z_Py_bb;
  Double I_NAI_K2x5y_D2y_bb = I_NAI_L2x6y_Py_bb+ABY*I_NAI_K2x5y_Py_bb;
  Double I_NAI_K2x4yz_D2y_bb = I_NAI_L2x5yz_Py_bb+ABY*I_NAI_K2x4yz_Py_bb;
  Double I_NAI_K2x3y2z_D2y_bb = I_NAI_L2x4y2z_Py_bb+ABY*I_NAI_K2x3y2z_Py_bb;
  Double I_NAI_K2x2y3z_D2y_bb = I_NAI_L2x3y3z_Py_bb+ABY*I_NAI_K2x2y3z_Py_bb;
  Double I_NAI_K2xy4z_D2y_bb = I_NAI_L2x2y4z_Py_bb+ABY*I_NAI_K2xy4z_Py_bb;
  Double I_NAI_K2x5z_D2y_bb = I_NAI_L2xy5z_Py_bb+ABY*I_NAI_K2x5z_Py_bb;
  Double I_NAI_Kx6y_D2y_bb = I_NAI_Lx7y_Py_bb+ABY*I_NAI_Kx6y_Py_bb;
  Double I_NAI_Kx5yz_D2y_bb = I_NAI_Lx6yz_Py_bb+ABY*I_NAI_Kx5yz_Py_bb;
  Double I_NAI_Kx4y2z_D2y_bb = I_NAI_Lx5y2z_Py_bb+ABY*I_NAI_Kx4y2z_Py_bb;
  Double I_NAI_Kx3y3z_D2y_bb = I_NAI_Lx4y3z_Py_bb+ABY*I_NAI_Kx3y3z_Py_bb;
  Double I_NAI_Kx2y4z_D2y_bb = I_NAI_Lx3y4z_Py_bb+ABY*I_NAI_Kx2y4z_Py_bb;
  Double I_NAI_Kxy5z_D2y_bb = I_NAI_Lx2y5z_Py_bb+ABY*I_NAI_Kxy5z_Py_bb;
  Double I_NAI_Kx6z_D2y_bb = I_NAI_Lxy6z_Py_bb+ABY*I_NAI_Kx6z_Py_bb;
  Double I_NAI_K7y_D2y_bb = I_NAI_L8y_Py_bb+ABY*I_NAI_K7y_Py_bb;
  Double I_NAI_K6yz_D2y_bb = I_NAI_L7yz_Py_bb+ABY*I_NAI_K6yz_Py_bb;
  Double I_NAI_K5y2z_D2y_bb = I_NAI_L6y2z_Py_bb+ABY*I_NAI_K5y2z_Py_bb;
  Double I_NAI_K4y3z_D2y_bb = I_NAI_L5y3z_Py_bb+ABY*I_NAI_K4y3z_Py_bb;
  Double I_NAI_K3y4z_D2y_bb = I_NAI_L4y4z_Py_bb+ABY*I_NAI_K3y4z_Py_bb;
  Double I_NAI_K2y5z_D2y_bb = I_NAI_L3y5z_Py_bb+ABY*I_NAI_K2y5z_Py_bb;
  Double I_NAI_Ky6z_D2y_bb = I_NAI_L2y6z_Py_bb+ABY*I_NAI_Ky6z_Py_bb;
  Double I_NAI_K6xz_D2z_bb = I_NAI_L6x2z_Pz_bb+ABZ*I_NAI_K6xz_Pz_bb;
  Double I_NAI_K5xyz_D2z_bb = I_NAI_L5xy2z_Pz_bb+ABZ*I_NAI_K5xyz_Pz_bb;
  Double I_NAI_K5x2z_D2z_bb = I_NAI_L5x3z_Pz_bb+ABZ*I_NAI_K5x2z_Pz_bb;
  Double I_NAI_K4x2yz_D2z_bb = I_NAI_L4x2y2z_Pz_bb+ABZ*I_NAI_K4x2yz_Pz_bb;
  Double I_NAI_K4xy2z_D2z_bb = I_NAI_L4xy3z_Pz_bb+ABZ*I_NAI_K4xy2z_Pz_bb;
  Double I_NAI_K4x3z_D2z_bb = I_NAI_L4x4z_Pz_bb+ABZ*I_NAI_K4x3z_Pz_bb;
  Double I_NAI_K3x3yz_D2z_bb = I_NAI_L3x3y2z_Pz_bb+ABZ*I_NAI_K3x3yz_Pz_bb;
  Double I_NAI_K3x2y2z_D2z_bb = I_NAI_L3x2y3z_Pz_bb+ABZ*I_NAI_K3x2y2z_Pz_bb;
  Double I_NAI_K3xy3z_D2z_bb = I_NAI_L3xy4z_Pz_bb+ABZ*I_NAI_K3xy3z_Pz_bb;
  Double I_NAI_K3x4z_D2z_bb = I_NAI_L3x5z_Pz_bb+ABZ*I_NAI_K3x4z_Pz_bb;
  Double I_NAI_K2x4yz_D2z_bb = I_NAI_L2x4y2z_Pz_bb+ABZ*I_NAI_K2x4yz_Pz_bb;
  Double I_NAI_K2x3y2z_D2z_bb = I_NAI_L2x3y3z_Pz_bb+ABZ*I_NAI_K2x3y2z_Pz_bb;
  Double I_NAI_K2x2y3z_D2z_bb = I_NAI_L2x2y4z_Pz_bb+ABZ*I_NAI_K2x2y3z_Pz_bb;
  Double I_NAI_K2xy4z_D2z_bb = I_NAI_L2xy5z_Pz_bb+ABZ*I_NAI_K2xy4z_Pz_bb;
  Double I_NAI_K2x5z_D2z_bb = I_NAI_L2x6z_Pz_bb+ABZ*I_NAI_K2x5z_Pz_bb;
  Double I_NAI_Kx5yz_D2z_bb = I_NAI_Lx5y2z_Pz_bb+ABZ*I_NAI_Kx5yz_Pz_bb;
  Double I_NAI_Kx4y2z_D2z_bb = I_NAI_Lx4y3z_Pz_bb+ABZ*I_NAI_Kx4y2z_Pz_bb;
  Double I_NAI_Kx3y3z_D2z_bb = I_NAI_Lx3y4z_Pz_bb+ABZ*I_NAI_Kx3y3z_Pz_bb;
  Double I_NAI_Kx2y4z_D2z_bb = I_NAI_Lx2y5z_Pz_bb+ABZ*I_NAI_Kx2y4z_Pz_bb;
  Double I_NAI_Kxy5z_D2z_bb = I_NAI_Lxy6z_Pz_bb+ABZ*I_NAI_Kxy5z_Pz_bb;
  Double I_NAI_Kx6z_D2z_bb = I_NAI_Lx7z_Pz_bb+ABZ*I_NAI_Kx6z_Pz_bb;
  Double I_NAI_K6yz_D2z_bb = I_NAI_L6y2z_Pz_bb+ABZ*I_NAI_K6yz_Pz_bb;
  Double I_NAI_K5y2z_D2z_bb = I_NAI_L5y3z_Pz_bb+ABZ*I_NAI_K5y2z_Pz_bb;
  Double I_NAI_K4y3z_D2z_bb = I_NAI_L4y4z_Pz_bb+ABZ*I_NAI_K4y3z_Pz_bb;
  Double I_NAI_K3y4z_D2z_bb = I_NAI_L3y5z_Pz_bb+ABZ*I_NAI_K3y4z_Pz_bb;
  Double I_NAI_K2y5z_D2z_bb = I_NAI_L2y6z_Pz_bb+ABZ*I_NAI_K2y5z_Pz_bb;
  Double I_NAI_Ky6z_D2z_bb = I_NAI_Ly7z_Pz_bb+ABZ*I_NAI_Ky6z_Pz_bb;
  Double I_NAI_K7z_D2z_bb = I_NAI_L8z_Pz_bb+ABZ*I_NAI_K7z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_I_F_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 151 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_D_bb
   * RHS shell quartet name: SQ_NAI_I_D_bb
   ************************************************************/
  Double I_NAI_I6x_F3x_bb = I_NAI_K7x_D2x_bb+ABX*I_NAI_I6x_D2x_bb;
  Double I_NAI_I5xy_F3x_bb = I_NAI_K6xy_D2x_bb+ABX*I_NAI_I5xy_D2x_bb;
  Double I_NAI_I5xz_F3x_bb = I_NAI_K6xz_D2x_bb+ABX*I_NAI_I5xz_D2x_bb;
  Double I_NAI_I4x2y_F3x_bb = I_NAI_K5x2y_D2x_bb+ABX*I_NAI_I4x2y_D2x_bb;
  Double I_NAI_I4xyz_F3x_bb = I_NAI_K5xyz_D2x_bb+ABX*I_NAI_I4xyz_D2x_bb;
  Double I_NAI_I4x2z_F3x_bb = I_NAI_K5x2z_D2x_bb+ABX*I_NAI_I4x2z_D2x_bb;
  Double I_NAI_I3x3y_F3x_bb = I_NAI_K4x3y_D2x_bb+ABX*I_NAI_I3x3y_D2x_bb;
  Double I_NAI_I3x2yz_F3x_bb = I_NAI_K4x2yz_D2x_bb+ABX*I_NAI_I3x2yz_D2x_bb;
  Double I_NAI_I3xy2z_F3x_bb = I_NAI_K4xy2z_D2x_bb+ABX*I_NAI_I3xy2z_D2x_bb;
  Double I_NAI_I3x3z_F3x_bb = I_NAI_K4x3z_D2x_bb+ABX*I_NAI_I3x3z_D2x_bb;
  Double I_NAI_I2x4y_F3x_bb = I_NAI_K3x4y_D2x_bb+ABX*I_NAI_I2x4y_D2x_bb;
  Double I_NAI_I2x3yz_F3x_bb = I_NAI_K3x3yz_D2x_bb+ABX*I_NAI_I2x3yz_D2x_bb;
  Double I_NAI_I2x2y2z_F3x_bb = I_NAI_K3x2y2z_D2x_bb+ABX*I_NAI_I2x2y2z_D2x_bb;
  Double I_NAI_I2xy3z_F3x_bb = I_NAI_K3xy3z_D2x_bb+ABX*I_NAI_I2xy3z_D2x_bb;
  Double I_NAI_I2x4z_F3x_bb = I_NAI_K3x4z_D2x_bb+ABX*I_NAI_I2x4z_D2x_bb;
  Double I_NAI_Ix5y_F3x_bb = I_NAI_K2x5y_D2x_bb+ABX*I_NAI_Ix5y_D2x_bb;
  Double I_NAI_Ix4yz_F3x_bb = I_NAI_K2x4yz_D2x_bb+ABX*I_NAI_Ix4yz_D2x_bb;
  Double I_NAI_Ix3y2z_F3x_bb = I_NAI_K2x3y2z_D2x_bb+ABX*I_NAI_Ix3y2z_D2x_bb;
  Double I_NAI_Ix2y3z_F3x_bb = I_NAI_K2x2y3z_D2x_bb+ABX*I_NAI_Ix2y3z_D2x_bb;
  Double I_NAI_Ixy4z_F3x_bb = I_NAI_K2xy4z_D2x_bb+ABX*I_NAI_Ixy4z_D2x_bb;
  Double I_NAI_Ix5z_F3x_bb = I_NAI_K2x5z_D2x_bb+ABX*I_NAI_Ix5z_D2x_bb;
  Double I_NAI_I6y_F3x_bb = I_NAI_Kx6y_D2x_bb+ABX*I_NAI_I6y_D2x_bb;
  Double I_NAI_I5yz_F3x_bb = I_NAI_Kx5yz_D2x_bb+ABX*I_NAI_I5yz_D2x_bb;
  Double I_NAI_I4y2z_F3x_bb = I_NAI_Kx4y2z_D2x_bb+ABX*I_NAI_I4y2z_D2x_bb;
  Double I_NAI_I3y3z_F3x_bb = I_NAI_Kx3y3z_D2x_bb+ABX*I_NAI_I3y3z_D2x_bb;
  Double I_NAI_I2y4z_F3x_bb = I_NAI_Kx2y4z_D2x_bb+ABX*I_NAI_I2y4z_D2x_bb;
  Double I_NAI_Iy5z_F3x_bb = I_NAI_Kxy5z_D2x_bb+ABX*I_NAI_Iy5z_D2x_bb;
  Double I_NAI_I6z_F3x_bb = I_NAI_Kx6z_D2x_bb+ABX*I_NAI_I6z_D2x_bb;
  Double I_NAI_I4xyz_F2xy_bb = I_NAI_K4x2yz_D2x_bb+ABY*I_NAI_I4xyz_D2x_bb;
  Double I_NAI_I3x2yz_F2xy_bb = I_NAI_K3x3yz_D2x_bb+ABY*I_NAI_I3x2yz_D2x_bb;
  Double I_NAI_I3xy2z_F2xy_bb = I_NAI_K3x2y2z_D2x_bb+ABY*I_NAI_I3xy2z_D2x_bb;
  Double I_NAI_I2x3yz_F2xy_bb = I_NAI_K2x4yz_D2x_bb+ABY*I_NAI_I2x3yz_D2x_bb;
  Double I_NAI_I2x2y2z_F2xy_bb = I_NAI_K2x3y2z_D2x_bb+ABY*I_NAI_I2x2y2z_D2x_bb;
  Double I_NAI_I2xy3z_F2xy_bb = I_NAI_K2x2y3z_D2x_bb+ABY*I_NAI_I2xy3z_D2x_bb;
  Double I_NAI_Ix4yz_F2xy_bb = I_NAI_Kx5yz_D2x_bb+ABY*I_NAI_Ix4yz_D2x_bb;
  Double I_NAI_Ix3y2z_F2xy_bb = I_NAI_Kx4y2z_D2x_bb+ABY*I_NAI_Ix3y2z_D2x_bb;
  Double I_NAI_Ix2y3z_F2xy_bb = I_NAI_Kx3y3z_D2x_bb+ABY*I_NAI_Ix2y3z_D2x_bb;
  Double I_NAI_Ixy4z_F2xy_bb = I_NAI_Kx2y4z_D2x_bb+ABY*I_NAI_Ixy4z_D2x_bb;
  Double I_NAI_I5yz_F2xy_bb = I_NAI_K6yz_D2x_bb+ABY*I_NAI_I5yz_D2x_bb;
  Double I_NAI_I4y2z_F2xy_bb = I_NAI_K5y2z_D2x_bb+ABY*I_NAI_I4y2z_D2x_bb;
  Double I_NAI_I3y3z_F2xy_bb = I_NAI_K4y3z_D2x_bb+ABY*I_NAI_I3y3z_D2x_bb;
  Double I_NAI_I2y4z_F2xy_bb = I_NAI_K3y4z_D2x_bb+ABY*I_NAI_I2y4z_D2x_bb;
  Double I_NAI_Iy5z_F2xy_bb = I_NAI_K2y5z_D2x_bb+ABY*I_NAI_Iy5z_D2x_bb;
  Double I_NAI_I4xyz_F2xz_bb = I_NAI_K4xy2z_D2x_bb+ABZ*I_NAI_I4xyz_D2x_bb;
  Double I_NAI_I3x2yz_F2xz_bb = I_NAI_K3x2y2z_D2x_bb+ABZ*I_NAI_I3x2yz_D2x_bb;
  Double I_NAI_I3xy2z_F2xz_bb = I_NAI_K3xy3z_D2x_bb+ABZ*I_NAI_I3xy2z_D2x_bb;
  Double I_NAI_I2x3yz_F2xz_bb = I_NAI_K2x3y2z_D2x_bb+ABZ*I_NAI_I2x3yz_D2x_bb;
  Double I_NAI_I2x2y2z_F2xz_bb = I_NAI_K2x2y3z_D2x_bb+ABZ*I_NAI_I2x2y2z_D2x_bb;
  Double I_NAI_I2xy3z_F2xz_bb = I_NAI_K2xy4z_D2x_bb+ABZ*I_NAI_I2xy3z_D2x_bb;
  Double I_NAI_Ix4yz_F2xz_bb = I_NAI_Kx4y2z_D2x_bb+ABZ*I_NAI_Ix4yz_D2x_bb;
  Double I_NAI_Ix3y2z_F2xz_bb = I_NAI_Kx3y3z_D2x_bb+ABZ*I_NAI_Ix3y2z_D2x_bb;
  Double I_NAI_Ix2y3z_F2xz_bb = I_NAI_Kx2y4z_D2x_bb+ABZ*I_NAI_Ix2y3z_D2x_bb;
  Double I_NAI_Ixy4z_F2xz_bb = I_NAI_Kxy5z_D2x_bb+ABZ*I_NAI_Ixy4z_D2x_bb;
  Double I_NAI_I5yz_F2xz_bb = I_NAI_K5y2z_D2x_bb+ABZ*I_NAI_I5yz_D2x_bb;
  Double I_NAI_I4y2z_F2xz_bb = I_NAI_K4y3z_D2x_bb+ABZ*I_NAI_I4y2z_D2x_bb;
  Double I_NAI_I3y3z_F2xz_bb = I_NAI_K3y4z_D2x_bb+ABZ*I_NAI_I3y3z_D2x_bb;
  Double I_NAI_I2y4z_F2xz_bb = I_NAI_K2y5z_D2x_bb+ABZ*I_NAI_I2y4z_D2x_bb;
  Double I_NAI_Iy5z_F2xz_bb = I_NAI_Ky6z_D2x_bb+ABZ*I_NAI_Iy5z_D2x_bb;
  Double I_NAI_I6x_F3y_bb = I_NAI_K6xy_D2y_bb+ABY*I_NAI_I6x_D2y_bb;
  Double I_NAI_I5xy_F3y_bb = I_NAI_K5x2y_D2y_bb+ABY*I_NAI_I5xy_D2y_bb;
  Double I_NAI_I5xz_F3y_bb = I_NAI_K5xyz_D2y_bb+ABY*I_NAI_I5xz_D2y_bb;
  Double I_NAI_I4x2y_F3y_bb = I_NAI_K4x3y_D2y_bb+ABY*I_NAI_I4x2y_D2y_bb;
  Double I_NAI_I4xyz_F3y_bb = I_NAI_K4x2yz_D2y_bb+ABY*I_NAI_I4xyz_D2y_bb;
  Double I_NAI_I4x2z_F3y_bb = I_NAI_K4xy2z_D2y_bb+ABY*I_NAI_I4x2z_D2y_bb;
  Double I_NAI_I3x3y_F3y_bb = I_NAI_K3x4y_D2y_bb+ABY*I_NAI_I3x3y_D2y_bb;
  Double I_NAI_I3x2yz_F3y_bb = I_NAI_K3x3yz_D2y_bb+ABY*I_NAI_I3x2yz_D2y_bb;
  Double I_NAI_I3xy2z_F3y_bb = I_NAI_K3x2y2z_D2y_bb+ABY*I_NAI_I3xy2z_D2y_bb;
  Double I_NAI_I3x3z_F3y_bb = I_NAI_K3xy3z_D2y_bb+ABY*I_NAI_I3x3z_D2y_bb;
  Double I_NAI_I2x4y_F3y_bb = I_NAI_K2x5y_D2y_bb+ABY*I_NAI_I2x4y_D2y_bb;
  Double I_NAI_I2x3yz_F3y_bb = I_NAI_K2x4yz_D2y_bb+ABY*I_NAI_I2x3yz_D2y_bb;
  Double I_NAI_I2x2y2z_F3y_bb = I_NAI_K2x3y2z_D2y_bb+ABY*I_NAI_I2x2y2z_D2y_bb;
  Double I_NAI_I2xy3z_F3y_bb = I_NAI_K2x2y3z_D2y_bb+ABY*I_NAI_I2xy3z_D2y_bb;
  Double I_NAI_I2x4z_F3y_bb = I_NAI_K2xy4z_D2y_bb+ABY*I_NAI_I2x4z_D2y_bb;
  Double I_NAI_Ix5y_F3y_bb = I_NAI_Kx6y_D2y_bb+ABY*I_NAI_Ix5y_D2y_bb;
  Double I_NAI_Ix4yz_F3y_bb = I_NAI_Kx5yz_D2y_bb+ABY*I_NAI_Ix4yz_D2y_bb;
  Double I_NAI_Ix3y2z_F3y_bb = I_NAI_Kx4y2z_D2y_bb+ABY*I_NAI_Ix3y2z_D2y_bb;
  Double I_NAI_Ix2y3z_F3y_bb = I_NAI_Kx3y3z_D2y_bb+ABY*I_NAI_Ix2y3z_D2y_bb;
  Double I_NAI_Ixy4z_F3y_bb = I_NAI_Kx2y4z_D2y_bb+ABY*I_NAI_Ixy4z_D2y_bb;
  Double I_NAI_Ix5z_F3y_bb = I_NAI_Kxy5z_D2y_bb+ABY*I_NAI_Ix5z_D2y_bb;
  Double I_NAI_I6y_F3y_bb = I_NAI_K7y_D2y_bb+ABY*I_NAI_I6y_D2y_bb;
  Double I_NAI_I5yz_F3y_bb = I_NAI_K6yz_D2y_bb+ABY*I_NAI_I5yz_D2y_bb;
  Double I_NAI_I4y2z_F3y_bb = I_NAI_K5y2z_D2y_bb+ABY*I_NAI_I4y2z_D2y_bb;
  Double I_NAI_I3y3z_F3y_bb = I_NAI_K4y3z_D2y_bb+ABY*I_NAI_I3y3z_D2y_bb;
  Double I_NAI_I2y4z_F3y_bb = I_NAI_K3y4z_D2y_bb+ABY*I_NAI_I2y4z_D2y_bb;
  Double I_NAI_Iy5z_F3y_bb = I_NAI_K2y5z_D2y_bb+ABY*I_NAI_Iy5z_D2y_bb;
  Double I_NAI_I6z_F3y_bb = I_NAI_Ky6z_D2y_bb+ABY*I_NAI_I6z_D2y_bb;
  Double I_NAI_I5xz_F2yz_bb = I_NAI_K5x2z_D2y_bb+ABZ*I_NAI_I5xz_D2y_bb;
  Double I_NAI_I4xyz_F2yz_bb = I_NAI_K4xy2z_D2y_bb+ABZ*I_NAI_I4xyz_D2y_bb;
  Double I_NAI_I4x2z_F2yz_bb = I_NAI_K4x3z_D2y_bb+ABZ*I_NAI_I4x2z_D2y_bb;
  Double I_NAI_I3x2yz_F2yz_bb = I_NAI_K3x2y2z_D2y_bb+ABZ*I_NAI_I3x2yz_D2y_bb;
  Double I_NAI_I3xy2z_F2yz_bb = I_NAI_K3xy3z_D2y_bb+ABZ*I_NAI_I3xy2z_D2y_bb;
  Double I_NAI_I3x3z_F2yz_bb = I_NAI_K3x4z_D2y_bb+ABZ*I_NAI_I3x3z_D2y_bb;
  Double I_NAI_I2x3yz_F2yz_bb = I_NAI_K2x3y2z_D2y_bb+ABZ*I_NAI_I2x3yz_D2y_bb;
  Double I_NAI_I2x2y2z_F2yz_bb = I_NAI_K2x2y3z_D2y_bb+ABZ*I_NAI_I2x2y2z_D2y_bb;
  Double I_NAI_I2xy3z_F2yz_bb = I_NAI_K2xy4z_D2y_bb+ABZ*I_NAI_I2xy3z_D2y_bb;
  Double I_NAI_I2x4z_F2yz_bb = I_NAI_K2x5z_D2y_bb+ABZ*I_NAI_I2x4z_D2y_bb;
  Double I_NAI_Ix4yz_F2yz_bb = I_NAI_Kx4y2z_D2y_bb+ABZ*I_NAI_Ix4yz_D2y_bb;
  Double I_NAI_Ix3y2z_F2yz_bb = I_NAI_Kx3y3z_D2y_bb+ABZ*I_NAI_Ix3y2z_D2y_bb;
  Double I_NAI_Ix2y3z_F2yz_bb = I_NAI_Kx2y4z_D2y_bb+ABZ*I_NAI_Ix2y3z_D2y_bb;
  Double I_NAI_Ixy4z_F2yz_bb = I_NAI_Kxy5z_D2y_bb+ABZ*I_NAI_Ixy4z_D2y_bb;
  Double I_NAI_Ix5z_F2yz_bb = I_NAI_Kx6z_D2y_bb+ABZ*I_NAI_Ix5z_D2y_bb;
  Double I_NAI_I6x_F3z_bb = I_NAI_K6xz_D2z_bb+ABZ*I_NAI_I6x_D2z_bb;
  Double I_NAI_I5xy_F3z_bb = I_NAI_K5xyz_D2z_bb+ABZ*I_NAI_I5xy_D2z_bb;
  Double I_NAI_I5xz_F3z_bb = I_NAI_K5x2z_D2z_bb+ABZ*I_NAI_I5xz_D2z_bb;
  Double I_NAI_I4x2y_F3z_bb = I_NAI_K4x2yz_D2z_bb+ABZ*I_NAI_I4x2y_D2z_bb;
  Double I_NAI_I4xyz_F3z_bb = I_NAI_K4xy2z_D2z_bb+ABZ*I_NAI_I4xyz_D2z_bb;
  Double I_NAI_I4x2z_F3z_bb = I_NAI_K4x3z_D2z_bb+ABZ*I_NAI_I4x2z_D2z_bb;
  Double I_NAI_I3x3y_F3z_bb = I_NAI_K3x3yz_D2z_bb+ABZ*I_NAI_I3x3y_D2z_bb;
  Double I_NAI_I3x2yz_F3z_bb = I_NAI_K3x2y2z_D2z_bb+ABZ*I_NAI_I3x2yz_D2z_bb;
  Double I_NAI_I3xy2z_F3z_bb = I_NAI_K3xy3z_D2z_bb+ABZ*I_NAI_I3xy2z_D2z_bb;
  Double I_NAI_I3x3z_F3z_bb = I_NAI_K3x4z_D2z_bb+ABZ*I_NAI_I3x3z_D2z_bb;
  Double I_NAI_I2x4y_F3z_bb = I_NAI_K2x4yz_D2z_bb+ABZ*I_NAI_I2x4y_D2z_bb;
  Double I_NAI_I2x3yz_F3z_bb = I_NAI_K2x3y2z_D2z_bb+ABZ*I_NAI_I2x3yz_D2z_bb;
  Double I_NAI_I2x2y2z_F3z_bb = I_NAI_K2x2y3z_D2z_bb+ABZ*I_NAI_I2x2y2z_D2z_bb;
  Double I_NAI_I2xy3z_F3z_bb = I_NAI_K2xy4z_D2z_bb+ABZ*I_NAI_I2xy3z_D2z_bb;
  Double I_NAI_I2x4z_F3z_bb = I_NAI_K2x5z_D2z_bb+ABZ*I_NAI_I2x4z_D2z_bb;
  Double I_NAI_Ix5y_F3z_bb = I_NAI_Kx5yz_D2z_bb+ABZ*I_NAI_Ix5y_D2z_bb;
  Double I_NAI_Ix4yz_F3z_bb = I_NAI_Kx4y2z_D2z_bb+ABZ*I_NAI_Ix4yz_D2z_bb;
  Double I_NAI_Ix3y2z_F3z_bb = I_NAI_Kx3y3z_D2z_bb+ABZ*I_NAI_Ix3y2z_D2z_bb;
  Double I_NAI_Ix2y3z_F3z_bb = I_NAI_Kx2y4z_D2z_bb+ABZ*I_NAI_Ix2y3z_D2z_bb;
  Double I_NAI_Ixy4z_F3z_bb = I_NAI_Kxy5z_D2z_bb+ABZ*I_NAI_Ixy4z_D2z_bb;
  Double I_NAI_Ix5z_F3z_bb = I_NAI_Kx6z_D2z_bb+ABZ*I_NAI_Ix5z_D2z_bb;
  Double I_NAI_I6y_F3z_bb = I_NAI_K6yz_D2z_bb+ABZ*I_NAI_I6y_D2z_bb;
  Double I_NAI_I5yz_F3z_bb = I_NAI_K5y2z_D2z_bb+ABZ*I_NAI_I5yz_D2z_bb;
  Double I_NAI_I4y2z_F3z_bb = I_NAI_K4y3z_D2z_bb+ABZ*I_NAI_I4y2z_D2z_bb;
  Double I_NAI_I3y3z_F3z_bb = I_NAI_K3y4z_D2z_bb+ABZ*I_NAI_I3y3z_D2z_bb;
  Double I_NAI_I2y4z_F3z_bb = I_NAI_K2y5z_D2z_bb+ABZ*I_NAI_I2y4z_D2z_bb;
  Double I_NAI_Iy5z_F3z_bb = I_NAI_Ky6z_D2z_bb+ABZ*I_NAI_Iy5z_D2z_bb;
  Double I_NAI_I6z_F3z_bb = I_NAI_K7z_D2z_bb+ABZ*I_NAI_I6z_D2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 102 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_bb
   * RHS shell quartet name: SQ_NAI_H_F_bb
   ************************************************************/
  Double I_NAI_H5x_G4x_bb = I_NAI_I6x_F3x_bb+ABX*I_NAI_H5x_F3x_bb;
  Double I_NAI_H4xy_G4x_bb = I_NAI_I5xy_F3x_bb+ABX*I_NAI_H4xy_F3x_bb;
  Double I_NAI_H4xz_G4x_bb = I_NAI_I5xz_F3x_bb+ABX*I_NAI_H4xz_F3x_bb;
  Double I_NAI_H3x2y_G4x_bb = I_NAI_I4x2y_F3x_bb+ABX*I_NAI_H3x2y_F3x_bb;
  Double I_NAI_H3xyz_G4x_bb = I_NAI_I4xyz_F3x_bb+ABX*I_NAI_H3xyz_F3x_bb;
  Double I_NAI_H3x2z_G4x_bb = I_NAI_I4x2z_F3x_bb+ABX*I_NAI_H3x2z_F3x_bb;
  Double I_NAI_H2x3y_G4x_bb = I_NAI_I3x3y_F3x_bb+ABX*I_NAI_H2x3y_F3x_bb;
  Double I_NAI_H2x2yz_G4x_bb = I_NAI_I3x2yz_F3x_bb+ABX*I_NAI_H2x2yz_F3x_bb;
  Double I_NAI_H2xy2z_G4x_bb = I_NAI_I3xy2z_F3x_bb+ABX*I_NAI_H2xy2z_F3x_bb;
  Double I_NAI_H2x3z_G4x_bb = I_NAI_I3x3z_F3x_bb+ABX*I_NAI_H2x3z_F3x_bb;
  Double I_NAI_Hx4y_G4x_bb = I_NAI_I2x4y_F3x_bb+ABX*I_NAI_Hx4y_F3x_bb;
  Double I_NAI_Hx3yz_G4x_bb = I_NAI_I2x3yz_F3x_bb+ABX*I_NAI_Hx3yz_F3x_bb;
  Double I_NAI_Hx2y2z_G4x_bb = I_NAI_I2x2y2z_F3x_bb+ABX*I_NAI_Hx2y2z_F3x_bb;
  Double I_NAI_Hxy3z_G4x_bb = I_NAI_I2xy3z_F3x_bb+ABX*I_NAI_Hxy3z_F3x_bb;
  Double I_NAI_Hx4z_G4x_bb = I_NAI_I2x4z_F3x_bb+ABX*I_NAI_Hx4z_F3x_bb;
  Double I_NAI_H5y_G4x_bb = I_NAI_Ix5y_F3x_bb+ABX*I_NAI_H5y_F3x_bb;
  Double I_NAI_H4yz_G4x_bb = I_NAI_Ix4yz_F3x_bb+ABX*I_NAI_H4yz_F3x_bb;
  Double I_NAI_H3y2z_G4x_bb = I_NAI_Ix3y2z_F3x_bb+ABX*I_NAI_H3y2z_F3x_bb;
  Double I_NAI_H2y3z_G4x_bb = I_NAI_Ix2y3z_F3x_bb+ABX*I_NAI_H2y3z_F3x_bb;
  Double I_NAI_Hy4z_G4x_bb = I_NAI_Ixy4z_F3x_bb+ABX*I_NAI_Hy4z_F3x_bb;
  Double I_NAI_H5z_G4x_bb = I_NAI_Ix5z_F3x_bb+ABX*I_NAI_H5z_F3x_bb;
  Double I_NAI_H4xy_G3xy_bb = I_NAI_I4x2y_F3x_bb+ABY*I_NAI_H4xy_F3x_bb;
  Double I_NAI_H4xz_G3xy_bb = I_NAI_I4xyz_F3x_bb+ABY*I_NAI_H4xz_F3x_bb;
  Double I_NAI_H3x2y_G3xy_bb = I_NAI_I3x3y_F3x_bb+ABY*I_NAI_H3x2y_F3x_bb;
  Double I_NAI_H3xyz_G3xy_bb = I_NAI_I3x2yz_F3x_bb+ABY*I_NAI_H3xyz_F3x_bb;
  Double I_NAI_H3x2z_G3xy_bb = I_NAI_I3xy2z_F3x_bb+ABY*I_NAI_H3x2z_F3x_bb;
  Double I_NAI_H2x3y_G3xy_bb = I_NAI_I2x4y_F3x_bb+ABY*I_NAI_H2x3y_F3x_bb;
  Double I_NAI_H2x2yz_G3xy_bb = I_NAI_I2x3yz_F3x_bb+ABY*I_NAI_H2x2yz_F3x_bb;
  Double I_NAI_H2xy2z_G3xy_bb = I_NAI_I2x2y2z_F3x_bb+ABY*I_NAI_H2xy2z_F3x_bb;
  Double I_NAI_H2x3z_G3xy_bb = I_NAI_I2xy3z_F3x_bb+ABY*I_NAI_H2x3z_F3x_bb;
  Double I_NAI_Hx4y_G3xy_bb = I_NAI_Ix5y_F3x_bb+ABY*I_NAI_Hx4y_F3x_bb;
  Double I_NAI_Hx3yz_G3xy_bb = I_NAI_Ix4yz_F3x_bb+ABY*I_NAI_Hx3yz_F3x_bb;
  Double I_NAI_Hx2y2z_G3xy_bb = I_NAI_Ix3y2z_F3x_bb+ABY*I_NAI_Hx2y2z_F3x_bb;
  Double I_NAI_Hxy3z_G3xy_bb = I_NAI_Ix2y3z_F3x_bb+ABY*I_NAI_Hxy3z_F3x_bb;
  Double I_NAI_Hx4z_G3xy_bb = I_NAI_Ixy4z_F3x_bb+ABY*I_NAI_Hx4z_F3x_bb;
  Double I_NAI_H5y_G3xy_bb = I_NAI_I6y_F3x_bb+ABY*I_NAI_H5y_F3x_bb;
  Double I_NAI_H4yz_G3xy_bb = I_NAI_I5yz_F3x_bb+ABY*I_NAI_H4yz_F3x_bb;
  Double I_NAI_H3y2z_G3xy_bb = I_NAI_I4y2z_F3x_bb+ABY*I_NAI_H3y2z_F3x_bb;
  Double I_NAI_H2y3z_G3xy_bb = I_NAI_I3y3z_F3x_bb+ABY*I_NAI_H2y3z_F3x_bb;
  Double I_NAI_Hy4z_G3xy_bb = I_NAI_I2y4z_F3x_bb+ABY*I_NAI_Hy4z_F3x_bb;
  Double I_NAI_H5z_G3xy_bb = I_NAI_Iy5z_F3x_bb+ABY*I_NAI_H5z_F3x_bb;
  Double I_NAI_H4xz_G3xz_bb = I_NAI_I4x2z_F3x_bb+ABZ*I_NAI_H4xz_F3x_bb;
  Double I_NAI_H3xyz_G3xz_bb = I_NAI_I3xy2z_F3x_bb+ABZ*I_NAI_H3xyz_F3x_bb;
  Double I_NAI_H3x2z_G3xz_bb = I_NAI_I3x3z_F3x_bb+ABZ*I_NAI_H3x2z_F3x_bb;
  Double I_NAI_H2x2yz_G3xz_bb = I_NAI_I2x2y2z_F3x_bb+ABZ*I_NAI_H2x2yz_F3x_bb;
  Double I_NAI_H2xy2z_G3xz_bb = I_NAI_I2xy3z_F3x_bb+ABZ*I_NAI_H2xy2z_F3x_bb;
  Double I_NAI_H2x3z_G3xz_bb = I_NAI_I2x4z_F3x_bb+ABZ*I_NAI_H2x3z_F3x_bb;
  Double I_NAI_Hx3yz_G3xz_bb = I_NAI_Ix3y2z_F3x_bb+ABZ*I_NAI_Hx3yz_F3x_bb;
  Double I_NAI_Hx2y2z_G3xz_bb = I_NAI_Ix2y3z_F3x_bb+ABZ*I_NAI_Hx2y2z_F3x_bb;
  Double I_NAI_Hxy3z_G3xz_bb = I_NAI_Ixy4z_F3x_bb+ABZ*I_NAI_Hxy3z_F3x_bb;
  Double I_NAI_Hx4z_G3xz_bb = I_NAI_Ix5z_F3x_bb+ABZ*I_NAI_Hx4z_F3x_bb;
  Double I_NAI_H4yz_G3xz_bb = I_NAI_I4y2z_F3x_bb+ABZ*I_NAI_H4yz_F3x_bb;
  Double I_NAI_H3y2z_G3xz_bb = I_NAI_I3y3z_F3x_bb+ABZ*I_NAI_H3y2z_F3x_bb;
  Double I_NAI_H2y3z_G3xz_bb = I_NAI_I2y4z_F3x_bb+ABZ*I_NAI_H2y3z_F3x_bb;
  Double I_NAI_Hy4z_G3xz_bb = I_NAI_Iy5z_F3x_bb+ABZ*I_NAI_Hy4z_F3x_bb;
  Double I_NAI_H5z_G3xz_bb = I_NAI_I6z_F3x_bb+ABZ*I_NAI_H5z_F3x_bb;
  Double I_NAI_H4xz_G2x2y_bb = I_NAI_I4xyz_F2xy_bb+ABY*I_NAI_H4xz_F2xy_bb;
  Double I_NAI_H3xyz_G2x2y_bb = I_NAI_I3x2yz_F2xy_bb+ABY*I_NAI_H3xyz_F2xy_bb;
  Double I_NAI_H3x2z_G2x2y_bb = I_NAI_I3xy2z_F2xy_bb+ABY*I_NAI_H3x2z_F2xy_bb;
  Double I_NAI_H2x2yz_G2x2y_bb = I_NAI_I2x3yz_F2xy_bb+ABY*I_NAI_H2x2yz_F2xy_bb;
  Double I_NAI_H2xy2z_G2x2y_bb = I_NAI_I2x2y2z_F2xy_bb+ABY*I_NAI_H2xy2z_F2xy_bb;
  Double I_NAI_H2x3z_G2x2y_bb = I_NAI_I2xy3z_F2xy_bb+ABY*I_NAI_H2x3z_F2xy_bb;
  Double I_NAI_Hx3yz_G2x2y_bb = I_NAI_Ix4yz_F2xy_bb+ABY*I_NAI_Hx3yz_F2xy_bb;
  Double I_NAI_Hx2y2z_G2x2y_bb = I_NAI_Ix3y2z_F2xy_bb+ABY*I_NAI_Hx2y2z_F2xy_bb;
  Double I_NAI_Hxy3z_G2x2y_bb = I_NAI_Ix2y3z_F2xy_bb+ABY*I_NAI_Hxy3z_F2xy_bb;
  Double I_NAI_Hx4z_G2x2y_bb = I_NAI_Ixy4z_F2xy_bb+ABY*I_NAI_Hx4z_F2xy_bb;
  Double I_NAI_H4yz_G2x2y_bb = I_NAI_I5yz_F2xy_bb+ABY*I_NAI_H4yz_F2xy_bb;
  Double I_NAI_H3y2z_G2x2y_bb = I_NAI_I4y2z_F2xy_bb+ABY*I_NAI_H3y2z_F2xy_bb;
  Double I_NAI_H2y3z_G2x2y_bb = I_NAI_I3y3z_F2xy_bb+ABY*I_NAI_H2y3z_F2xy_bb;
  Double I_NAI_Hy4z_G2x2y_bb = I_NAI_I2y4z_F2xy_bb+ABY*I_NAI_Hy4z_F2xy_bb;
  Double I_NAI_H5z_G2x2y_bb = I_NAI_Iy5z_F2xy_bb+ABY*I_NAI_H5z_F2xy_bb;
  Double I_NAI_H4xy_G2x2z_bb = I_NAI_I4xyz_F2xz_bb+ABZ*I_NAI_H4xy_F2xz_bb;
  Double I_NAI_H3x2y_G2x2z_bb = I_NAI_I3x2yz_F2xz_bb+ABZ*I_NAI_H3x2y_F2xz_bb;
  Double I_NAI_H3xyz_G2x2z_bb = I_NAI_I3xy2z_F2xz_bb+ABZ*I_NAI_H3xyz_F2xz_bb;
  Double I_NAI_H2x3y_G2x2z_bb = I_NAI_I2x3yz_F2xz_bb+ABZ*I_NAI_H2x3y_F2xz_bb;
  Double I_NAI_H2x2yz_G2x2z_bb = I_NAI_I2x2y2z_F2xz_bb+ABZ*I_NAI_H2x2yz_F2xz_bb;
  Double I_NAI_H2xy2z_G2x2z_bb = I_NAI_I2xy3z_F2xz_bb+ABZ*I_NAI_H2xy2z_F2xz_bb;
  Double I_NAI_Hx4y_G2x2z_bb = I_NAI_Ix4yz_F2xz_bb+ABZ*I_NAI_Hx4y_F2xz_bb;
  Double I_NAI_Hx3yz_G2x2z_bb = I_NAI_Ix3y2z_F2xz_bb+ABZ*I_NAI_Hx3yz_F2xz_bb;
  Double I_NAI_Hx2y2z_G2x2z_bb = I_NAI_Ix2y3z_F2xz_bb+ABZ*I_NAI_Hx2y2z_F2xz_bb;
  Double I_NAI_Hxy3z_G2x2z_bb = I_NAI_Ixy4z_F2xz_bb+ABZ*I_NAI_Hxy3z_F2xz_bb;
  Double I_NAI_H5y_G2x2z_bb = I_NAI_I5yz_F2xz_bb+ABZ*I_NAI_H5y_F2xz_bb;
  Double I_NAI_H4yz_G2x2z_bb = I_NAI_I4y2z_F2xz_bb+ABZ*I_NAI_H4yz_F2xz_bb;
  Double I_NAI_H3y2z_G2x2z_bb = I_NAI_I3y3z_F2xz_bb+ABZ*I_NAI_H3y2z_F2xz_bb;
  Double I_NAI_H2y3z_G2x2z_bb = I_NAI_I2y4z_F2xz_bb+ABZ*I_NAI_H2y3z_F2xz_bb;
  Double I_NAI_Hy4z_G2x2z_bb = I_NAI_Iy5z_F2xz_bb+ABZ*I_NAI_Hy4z_F2xz_bb;
  Double I_NAI_H5x_Gx3y_bb = I_NAI_I6x_F3y_bb+ABX*I_NAI_H5x_F3y_bb;
  Double I_NAI_H4xy_Gx3y_bb = I_NAI_I5xy_F3y_bb+ABX*I_NAI_H4xy_F3y_bb;
  Double I_NAI_H4xz_Gx3y_bb = I_NAI_I5xz_F3y_bb+ABX*I_NAI_H4xz_F3y_bb;
  Double I_NAI_H3x2y_Gx3y_bb = I_NAI_I4x2y_F3y_bb+ABX*I_NAI_H3x2y_F3y_bb;
  Double I_NAI_H3xyz_Gx3y_bb = I_NAI_I4xyz_F3y_bb+ABX*I_NAI_H3xyz_F3y_bb;
  Double I_NAI_H3x2z_Gx3y_bb = I_NAI_I4x2z_F3y_bb+ABX*I_NAI_H3x2z_F3y_bb;
  Double I_NAI_H2x3y_Gx3y_bb = I_NAI_I3x3y_F3y_bb+ABX*I_NAI_H2x3y_F3y_bb;
  Double I_NAI_H2x2yz_Gx3y_bb = I_NAI_I3x2yz_F3y_bb+ABX*I_NAI_H2x2yz_F3y_bb;
  Double I_NAI_H2xy2z_Gx3y_bb = I_NAI_I3xy2z_F3y_bb+ABX*I_NAI_H2xy2z_F3y_bb;
  Double I_NAI_H2x3z_Gx3y_bb = I_NAI_I3x3z_F3y_bb+ABX*I_NAI_H2x3z_F3y_bb;
  Double I_NAI_Hx4y_Gx3y_bb = I_NAI_I2x4y_F3y_bb+ABX*I_NAI_Hx4y_F3y_bb;
  Double I_NAI_Hx3yz_Gx3y_bb = I_NAI_I2x3yz_F3y_bb+ABX*I_NAI_Hx3yz_F3y_bb;
  Double I_NAI_Hx2y2z_Gx3y_bb = I_NAI_I2x2y2z_F3y_bb+ABX*I_NAI_Hx2y2z_F3y_bb;
  Double I_NAI_Hxy3z_Gx3y_bb = I_NAI_I2xy3z_F3y_bb+ABX*I_NAI_Hxy3z_F3y_bb;
  Double I_NAI_Hx4z_Gx3y_bb = I_NAI_I2x4z_F3y_bb+ABX*I_NAI_Hx4z_F3y_bb;
  Double I_NAI_H4yz_Gx3y_bb = I_NAI_Ix4yz_F3y_bb+ABX*I_NAI_H4yz_F3y_bb;
  Double I_NAI_H3y2z_Gx3y_bb = I_NAI_Ix3y2z_F3y_bb+ABX*I_NAI_H3y2z_F3y_bb;
  Double I_NAI_H2y3z_Gx3y_bb = I_NAI_Ix2y3z_F3y_bb+ABX*I_NAI_H2y3z_F3y_bb;
  Double I_NAI_Hy4z_Gx3y_bb = I_NAI_Ixy4z_F3y_bb+ABX*I_NAI_Hy4z_F3y_bb;
  Double I_NAI_H5z_Gx3y_bb = I_NAI_Ix5z_F3y_bb+ABX*I_NAI_H5z_F3y_bb;
  Double I_NAI_H5x_Gx3z_bb = I_NAI_I6x_F3z_bb+ABX*I_NAI_H5x_F3z_bb;
  Double I_NAI_H4xy_Gx3z_bb = I_NAI_I5xy_F3z_bb+ABX*I_NAI_H4xy_F3z_bb;
  Double I_NAI_H4xz_Gx3z_bb = I_NAI_I5xz_F3z_bb+ABX*I_NAI_H4xz_F3z_bb;
  Double I_NAI_H3x2y_Gx3z_bb = I_NAI_I4x2y_F3z_bb+ABX*I_NAI_H3x2y_F3z_bb;
  Double I_NAI_H3xyz_Gx3z_bb = I_NAI_I4xyz_F3z_bb+ABX*I_NAI_H3xyz_F3z_bb;
  Double I_NAI_H3x2z_Gx3z_bb = I_NAI_I4x2z_F3z_bb+ABX*I_NAI_H3x2z_F3z_bb;
  Double I_NAI_H2x3y_Gx3z_bb = I_NAI_I3x3y_F3z_bb+ABX*I_NAI_H2x3y_F3z_bb;
  Double I_NAI_H2x2yz_Gx3z_bb = I_NAI_I3x2yz_F3z_bb+ABX*I_NAI_H2x2yz_F3z_bb;
  Double I_NAI_H2xy2z_Gx3z_bb = I_NAI_I3xy2z_F3z_bb+ABX*I_NAI_H2xy2z_F3z_bb;
  Double I_NAI_H2x3z_Gx3z_bb = I_NAI_I3x3z_F3z_bb+ABX*I_NAI_H2x3z_F3z_bb;
  Double I_NAI_Hx4y_Gx3z_bb = I_NAI_I2x4y_F3z_bb+ABX*I_NAI_Hx4y_F3z_bb;
  Double I_NAI_Hx3yz_Gx3z_bb = I_NAI_I2x3yz_F3z_bb+ABX*I_NAI_Hx3yz_F3z_bb;
  Double I_NAI_Hx2y2z_Gx3z_bb = I_NAI_I2x2y2z_F3z_bb+ABX*I_NAI_Hx2y2z_F3z_bb;
  Double I_NAI_Hxy3z_Gx3z_bb = I_NAI_I2xy3z_F3z_bb+ABX*I_NAI_Hxy3z_F3z_bb;
  Double I_NAI_Hx4z_Gx3z_bb = I_NAI_I2x4z_F3z_bb+ABX*I_NAI_Hx4z_F3z_bb;
  Double I_NAI_H5y_Gx3z_bb = I_NAI_Ix5y_F3z_bb+ABX*I_NAI_H5y_F3z_bb;
  Double I_NAI_H4yz_Gx3z_bb = I_NAI_Ix4yz_F3z_bb+ABX*I_NAI_H4yz_F3z_bb;
  Double I_NAI_H3y2z_Gx3z_bb = I_NAI_Ix3y2z_F3z_bb+ABX*I_NAI_H3y2z_F3z_bb;
  Double I_NAI_H2y3z_Gx3z_bb = I_NAI_Ix2y3z_F3z_bb+ABX*I_NAI_H2y3z_F3z_bb;
  Double I_NAI_Hy4z_Gx3z_bb = I_NAI_Ixy4z_F3z_bb+ABX*I_NAI_Hy4z_F3z_bb;
  Double I_NAI_H5x_G4y_bb = I_NAI_I5xy_F3y_bb+ABY*I_NAI_H5x_F3y_bb;
  Double I_NAI_H4xy_G4y_bb = I_NAI_I4x2y_F3y_bb+ABY*I_NAI_H4xy_F3y_bb;
  Double I_NAI_H4xz_G4y_bb = I_NAI_I4xyz_F3y_bb+ABY*I_NAI_H4xz_F3y_bb;
  Double I_NAI_H3x2y_G4y_bb = I_NAI_I3x3y_F3y_bb+ABY*I_NAI_H3x2y_F3y_bb;
  Double I_NAI_H3xyz_G4y_bb = I_NAI_I3x2yz_F3y_bb+ABY*I_NAI_H3xyz_F3y_bb;
  Double I_NAI_H3x2z_G4y_bb = I_NAI_I3xy2z_F3y_bb+ABY*I_NAI_H3x2z_F3y_bb;
  Double I_NAI_H2x3y_G4y_bb = I_NAI_I2x4y_F3y_bb+ABY*I_NAI_H2x3y_F3y_bb;
  Double I_NAI_H2x2yz_G4y_bb = I_NAI_I2x3yz_F3y_bb+ABY*I_NAI_H2x2yz_F3y_bb;
  Double I_NAI_H2xy2z_G4y_bb = I_NAI_I2x2y2z_F3y_bb+ABY*I_NAI_H2xy2z_F3y_bb;
  Double I_NAI_H2x3z_G4y_bb = I_NAI_I2xy3z_F3y_bb+ABY*I_NAI_H2x3z_F3y_bb;
  Double I_NAI_Hx4y_G4y_bb = I_NAI_Ix5y_F3y_bb+ABY*I_NAI_Hx4y_F3y_bb;
  Double I_NAI_Hx3yz_G4y_bb = I_NAI_Ix4yz_F3y_bb+ABY*I_NAI_Hx3yz_F3y_bb;
  Double I_NAI_Hx2y2z_G4y_bb = I_NAI_Ix3y2z_F3y_bb+ABY*I_NAI_Hx2y2z_F3y_bb;
  Double I_NAI_Hxy3z_G4y_bb = I_NAI_Ix2y3z_F3y_bb+ABY*I_NAI_Hxy3z_F3y_bb;
  Double I_NAI_Hx4z_G4y_bb = I_NAI_Ixy4z_F3y_bb+ABY*I_NAI_Hx4z_F3y_bb;
  Double I_NAI_H5y_G4y_bb = I_NAI_I6y_F3y_bb+ABY*I_NAI_H5y_F3y_bb;
  Double I_NAI_H4yz_G4y_bb = I_NAI_I5yz_F3y_bb+ABY*I_NAI_H4yz_F3y_bb;
  Double I_NAI_H3y2z_G4y_bb = I_NAI_I4y2z_F3y_bb+ABY*I_NAI_H3y2z_F3y_bb;
  Double I_NAI_H2y3z_G4y_bb = I_NAI_I3y3z_F3y_bb+ABY*I_NAI_H2y3z_F3y_bb;
  Double I_NAI_Hy4z_G4y_bb = I_NAI_I2y4z_F3y_bb+ABY*I_NAI_Hy4z_F3y_bb;
  Double I_NAI_H5z_G4y_bb = I_NAI_Iy5z_F3y_bb+ABY*I_NAI_H5z_F3y_bb;
  Double I_NAI_H4xz_G3yz_bb = I_NAI_I4x2z_F3y_bb+ABZ*I_NAI_H4xz_F3y_bb;
  Double I_NAI_H3xyz_G3yz_bb = I_NAI_I3xy2z_F3y_bb+ABZ*I_NAI_H3xyz_F3y_bb;
  Double I_NAI_H3x2z_G3yz_bb = I_NAI_I3x3z_F3y_bb+ABZ*I_NAI_H3x2z_F3y_bb;
  Double I_NAI_H2x2yz_G3yz_bb = I_NAI_I2x2y2z_F3y_bb+ABZ*I_NAI_H2x2yz_F3y_bb;
  Double I_NAI_H2xy2z_G3yz_bb = I_NAI_I2xy3z_F3y_bb+ABZ*I_NAI_H2xy2z_F3y_bb;
  Double I_NAI_H2x3z_G3yz_bb = I_NAI_I2x4z_F3y_bb+ABZ*I_NAI_H2x3z_F3y_bb;
  Double I_NAI_Hx3yz_G3yz_bb = I_NAI_Ix3y2z_F3y_bb+ABZ*I_NAI_Hx3yz_F3y_bb;
  Double I_NAI_Hx2y2z_G3yz_bb = I_NAI_Ix2y3z_F3y_bb+ABZ*I_NAI_Hx2y2z_F3y_bb;
  Double I_NAI_Hxy3z_G3yz_bb = I_NAI_Ixy4z_F3y_bb+ABZ*I_NAI_Hxy3z_F3y_bb;
  Double I_NAI_Hx4z_G3yz_bb = I_NAI_Ix5z_F3y_bb+ABZ*I_NAI_Hx4z_F3y_bb;
  Double I_NAI_H4yz_G3yz_bb = I_NAI_I4y2z_F3y_bb+ABZ*I_NAI_H4yz_F3y_bb;
  Double I_NAI_H3y2z_G3yz_bb = I_NAI_I3y3z_F3y_bb+ABZ*I_NAI_H3y2z_F3y_bb;
  Double I_NAI_H2y3z_G3yz_bb = I_NAI_I2y4z_F3y_bb+ABZ*I_NAI_H2y3z_F3y_bb;
  Double I_NAI_Hy4z_G3yz_bb = I_NAI_Iy5z_F3y_bb+ABZ*I_NAI_Hy4z_F3y_bb;
  Double I_NAI_H5z_G3yz_bb = I_NAI_I6z_F3y_bb+ABZ*I_NAI_H5z_F3y_bb;
  Double I_NAI_H5x_G2y2z_bb = I_NAI_I5xz_F2yz_bb+ABZ*I_NAI_H5x_F2yz_bb;
  Double I_NAI_H4xy_G2y2z_bb = I_NAI_I4xyz_F2yz_bb+ABZ*I_NAI_H4xy_F2yz_bb;
  Double I_NAI_H4xz_G2y2z_bb = I_NAI_I4x2z_F2yz_bb+ABZ*I_NAI_H4xz_F2yz_bb;
  Double I_NAI_H3x2y_G2y2z_bb = I_NAI_I3x2yz_F2yz_bb+ABZ*I_NAI_H3x2y_F2yz_bb;
  Double I_NAI_H3xyz_G2y2z_bb = I_NAI_I3xy2z_F2yz_bb+ABZ*I_NAI_H3xyz_F2yz_bb;
  Double I_NAI_H3x2z_G2y2z_bb = I_NAI_I3x3z_F2yz_bb+ABZ*I_NAI_H3x2z_F2yz_bb;
  Double I_NAI_H2x3y_G2y2z_bb = I_NAI_I2x3yz_F2yz_bb+ABZ*I_NAI_H2x3y_F2yz_bb;
  Double I_NAI_H2x2yz_G2y2z_bb = I_NAI_I2x2y2z_F2yz_bb+ABZ*I_NAI_H2x2yz_F2yz_bb;
  Double I_NAI_H2xy2z_G2y2z_bb = I_NAI_I2xy3z_F2yz_bb+ABZ*I_NAI_H2xy2z_F2yz_bb;
  Double I_NAI_H2x3z_G2y2z_bb = I_NAI_I2x4z_F2yz_bb+ABZ*I_NAI_H2x3z_F2yz_bb;
  Double I_NAI_Hx4y_G2y2z_bb = I_NAI_Ix4yz_F2yz_bb+ABZ*I_NAI_Hx4y_F2yz_bb;
  Double I_NAI_Hx3yz_G2y2z_bb = I_NAI_Ix3y2z_F2yz_bb+ABZ*I_NAI_Hx3yz_F2yz_bb;
  Double I_NAI_Hx2y2z_G2y2z_bb = I_NAI_Ix2y3z_F2yz_bb+ABZ*I_NAI_Hx2y2z_F2yz_bb;
  Double I_NAI_Hxy3z_G2y2z_bb = I_NAI_Ixy4z_F2yz_bb+ABZ*I_NAI_Hxy3z_F2yz_bb;
  Double I_NAI_Hx4z_G2y2z_bb = I_NAI_Ix5z_F2yz_bb+ABZ*I_NAI_Hx4z_F2yz_bb;
  Double I_NAI_H4xy_Gy3z_bb = I_NAI_I4x2y_F3z_bb+ABY*I_NAI_H4xy_F3z_bb;
  Double I_NAI_H3x2y_Gy3z_bb = I_NAI_I3x3y_F3z_bb+ABY*I_NAI_H3x2y_F3z_bb;
  Double I_NAI_H3xyz_Gy3z_bb = I_NAI_I3x2yz_F3z_bb+ABY*I_NAI_H3xyz_F3z_bb;
  Double I_NAI_H2x3y_Gy3z_bb = I_NAI_I2x4y_F3z_bb+ABY*I_NAI_H2x3y_F3z_bb;
  Double I_NAI_H2x2yz_Gy3z_bb = I_NAI_I2x3yz_F3z_bb+ABY*I_NAI_H2x2yz_F3z_bb;
  Double I_NAI_H2xy2z_Gy3z_bb = I_NAI_I2x2y2z_F3z_bb+ABY*I_NAI_H2xy2z_F3z_bb;
  Double I_NAI_Hx4y_Gy3z_bb = I_NAI_Ix5y_F3z_bb+ABY*I_NAI_Hx4y_F3z_bb;
  Double I_NAI_Hx3yz_Gy3z_bb = I_NAI_Ix4yz_F3z_bb+ABY*I_NAI_Hx3yz_F3z_bb;
  Double I_NAI_Hx2y2z_Gy3z_bb = I_NAI_Ix3y2z_F3z_bb+ABY*I_NAI_Hx2y2z_F3z_bb;
  Double I_NAI_Hxy3z_Gy3z_bb = I_NAI_Ix2y3z_F3z_bb+ABY*I_NAI_Hxy3z_F3z_bb;
  Double I_NAI_H5y_Gy3z_bb = I_NAI_I6y_F3z_bb+ABY*I_NAI_H5y_F3z_bb;
  Double I_NAI_H4yz_Gy3z_bb = I_NAI_I5yz_F3z_bb+ABY*I_NAI_H4yz_F3z_bb;
  Double I_NAI_H3y2z_Gy3z_bb = I_NAI_I4y2z_F3z_bb+ABY*I_NAI_H3y2z_F3z_bb;
  Double I_NAI_H2y3z_Gy3z_bb = I_NAI_I3y3z_F3z_bb+ABY*I_NAI_H2y3z_F3z_bb;
  Double I_NAI_Hy4z_Gy3z_bb = I_NAI_I2y4z_F3z_bb+ABY*I_NAI_Hy4z_F3z_bb;
  Double I_NAI_H5x_G4z_bb = I_NAI_I5xz_F3z_bb+ABZ*I_NAI_H5x_F3z_bb;
  Double I_NAI_H4xy_G4z_bb = I_NAI_I4xyz_F3z_bb+ABZ*I_NAI_H4xy_F3z_bb;
  Double I_NAI_H4xz_G4z_bb = I_NAI_I4x2z_F3z_bb+ABZ*I_NAI_H4xz_F3z_bb;
  Double I_NAI_H3x2y_G4z_bb = I_NAI_I3x2yz_F3z_bb+ABZ*I_NAI_H3x2y_F3z_bb;
  Double I_NAI_H3xyz_G4z_bb = I_NAI_I3xy2z_F3z_bb+ABZ*I_NAI_H3xyz_F3z_bb;
  Double I_NAI_H3x2z_G4z_bb = I_NAI_I3x3z_F3z_bb+ABZ*I_NAI_H3x2z_F3z_bb;
  Double I_NAI_H2x3y_G4z_bb = I_NAI_I2x3yz_F3z_bb+ABZ*I_NAI_H2x3y_F3z_bb;
  Double I_NAI_H2x2yz_G4z_bb = I_NAI_I2x2y2z_F3z_bb+ABZ*I_NAI_H2x2yz_F3z_bb;
  Double I_NAI_H2xy2z_G4z_bb = I_NAI_I2xy3z_F3z_bb+ABZ*I_NAI_H2xy2z_F3z_bb;
  Double I_NAI_H2x3z_G4z_bb = I_NAI_I2x4z_F3z_bb+ABZ*I_NAI_H2x3z_F3z_bb;
  Double I_NAI_Hx4y_G4z_bb = I_NAI_Ix4yz_F3z_bb+ABZ*I_NAI_Hx4y_F3z_bb;
  Double I_NAI_Hx3yz_G4z_bb = I_NAI_Ix3y2z_F3z_bb+ABZ*I_NAI_Hx3yz_F3z_bb;
  Double I_NAI_Hx2y2z_G4z_bb = I_NAI_Ix2y3z_F3z_bb+ABZ*I_NAI_Hx2y2z_F3z_bb;
  Double I_NAI_Hxy3z_G4z_bb = I_NAI_Ixy4z_F3z_bb+ABZ*I_NAI_Hxy3z_F3z_bb;
  Double I_NAI_Hx4z_G4z_bb = I_NAI_Ix5z_F3z_bb+ABZ*I_NAI_Hx4z_F3z_bb;
  Double I_NAI_H5y_G4z_bb = I_NAI_I5yz_F3z_bb+ABZ*I_NAI_H5y_F3z_bb;
  Double I_NAI_H4yz_G4z_bb = I_NAI_I4y2z_F3z_bb+ABZ*I_NAI_H4yz_F3z_bb;
  Double I_NAI_H3y2z_G4z_bb = I_NAI_I3y3z_F3z_bb+ABZ*I_NAI_H3y2z_F3z_bb;
  Double I_NAI_H2y3z_G4z_bb = I_NAI_I2y4z_F3z_bb+ABZ*I_NAI_H2y3z_F3z_bb;
  Double I_NAI_Hy4z_G4z_bb = I_NAI_Iy5z_F3z_bb+ABZ*I_NAI_Hy4z_F3z_bb;
  Double I_NAI_H5z_G4z_bb = I_NAI_I6z_F3z_bb+ABZ*I_NAI_H5z_F3z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_H_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_bb
   * RHS shell quartet name: SQ_NAI_G_G_bb
   ************************************************************/
  Double I_NAI_G4x_H5x_bb = I_NAI_H5x_G4x_bb+ABX*I_NAI_G4x_G4x_bb;
  Double I_NAI_G3xy_H5x_bb = I_NAI_H4xy_G4x_bb+ABX*I_NAI_G3xy_G4x_bb;
  Double I_NAI_G3xz_H5x_bb = I_NAI_H4xz_G4x_bb+ABX*I_NAI_G3xz_G4x_bb;
  Double I_NAI_G2x2y_H5x_bb = I_NAI_H3x2y_G4x_bb+ABX*I_NAI_G2x2y_G4x_bb;
  Double I_NAI_G2xyz_H5x_bb = I_NAI_H3xyz_G4x_bb+ABX*I_NAI_G2xyz_G4x_bb;
  Double I_NAI_G2x2z_H5x_bb = I_NAI_H3x2z_G4x_bb+ABX*I_NAI_G2x2z_G4x_bb;
  Double I_NAI_Gx3y_H5x_bb = I_NAI_H2x3y_G4x_bb+ABX*I_NAI_Gx3y_G4x_bb;
  Double I_NAI_Gx2yz_H5x_bb = I_NAI_H2x2yz_G4x_bb+ABX*I_NAI_Gx2yz_G4x_bb;
  Double I_NAI_Gxy2z_H5x_bb = I_NAI_H2xy2z_G4x_bb+ABX*I_NAI_Gxy2z_G4x_bb;
  Double I_NAI_Gx3z_H5x_bb = I_NAI_H2x3z_G4x_bb+ABX*I_NAI_Gx3z_G4x_bb;
  Double I_NAI_G4y_H5x_bb = I_NAI_Hx4y_G4x_bb+ABX*I_NAI_G4y_G4x_bb;
  Double I_NAI_G3yz_H5x_bb = I_NAI_Hx3yz_G4x_bb+ABX*I_NAI_G3yz_G4x_bb;
  Double I_NAI_G2y2z_H5x_bb = I_NAI_Hx2y2z_G4x_bb+ABX*I_NAI_G2y2z_G4x_bb;
  Double I_NAI_Gy3z_H5x_bb = I_NAI_Hxy3z_G4x_bb+ABX*I_NAI_Gy3z_G4x_bb;
  Double I_NAI_G4z_H5x_bb = I_NAI_Hx4z_G4x_bb+ABX*I_NAI_G4z_G4x_bb;
  Double I_NAI_G4x_H4xy_bb = I_NAI_H4xy_G4x_bb+ABY*I_NAI_G4x_G4x_bb;
  Double I_NAI_G3xy_H4xy_bb = I_NAI_H3x2y_G4x_bb+ABY*I_NAI_G3xy_G4x_bb;
  Double I_NAI_G3xz_H4xy_bb = I_NAI_H3xyz_G4x_bb+ABY*I_NAI_G3xz_G4x_bb;
  Double I_NAI_G2x2y_H4xy_bb = I_NAI_H2x3y_G4x_bb+ABY*I_NAI_G2x2y_G4x_bb;
  Double I_NAI_G2xyz_H4xy_bb = I_NAI_H2x2yz_G4x_bb+ABY*I_NAI_G2xyz_G4x_bb;
  Double I_NAI_G2x2z_H4xy_bb = I_NAI_H2xy2z_G4x_bb+ABY*I_NAI_G2x2z_G4x_bb;
  Double I_NAI_Gx3y_H4xy_bb = I_NAI_Hx4y_G4x_bb+ABY*I_NAI_Gx3y_G4x_bb;
  Double I_NAI_Gx2yz_H4xy_bb = I_NAI_Hx3yz_G4x_bb+ABY*I_NAI_Gx2yz_G4x_bb;
  Double I_NAI_Gxy2z_H4xy_bb = I_NAI_Hx2y2z_G4x_bb+ABY*I_NAI_Gxy2z_G4x_bb;
  Double I_NAI_Gx3z_H4xy_bb = I_NAI_Hxy3z_G4x_bb+ABY*I_NAI_Gx3z_G4x_bb;
  Double I_NAI_G4y_H4xy_bb = I_NAI_H5y_G4x_bb+ABY*I_NAI_G4y_G4x_bb;
  Double I_NAI_G3yz_H4xy_bb = I_NAI_H4yz_G4x_bb+ABY*I_NAI_G3yz_G4x_bb;
  Double I_NAI_G2y2z_H4xy_bb = I_NAI_H3y2z_G4x_bb+ABY*I_NAI_G2y2z_G4x_bb;
  Double I_NAI_Gy3z_H4xy_bb = I_NAI_H2y3z_G4x_bb+ABY*I_NAI_Gy3z_G4x_bb;
  Double I_NAI_G4z_H4xy_bb = I_NAI_Hy4z_G4x_bb+ABY*I_NAI_G4z_G4x_bb;
  Double I_NAI_G4x_H4xz_bb = I_NAI_H4xz_G4x_bb+ABZ*I_NAI_G4x_G4x_bb;
  Double I_NAI_G3xy_H4xz_bb = I_NAI_H3xyz_G4x_bb+ABZ*I_NAI_G3xy_G4x_bb;
  Double I_NAI_G3xz_H4xz_bb = I_NAI_H3x2z_G4x_bb+ABZ*I_NAI_G3xz_G4x_bb;
  Double I_NAI_G2x2y_H4xz_bb = I_NAI_H2x2yz_G4x_bb+ABZ*I_NAI_G2x2y_G4x_bb;
  Double I_NAI_G2xyz_H4xz_bb = I_NAI_H2xy2z_G4x_bb+ABZ*I_NAI_G2xyz_G4x_bb;
  Double I_NAI_G2x2z_H4xz_bb = I_NAI_H2x3z_G4x_bb+ABZ*I_NAI_G2x2z_G4x_bb;
  Double I_NAI_Gx3y_H4xz_bb = I_NAI_Hx3yz_G4x_bb+ABZ*I_NAI_Gx3y_G4x_bb;
  Double I_NAI_Gx2yz_H4xz_bb = I_NAI_Hx2y2z_G4x_bb+ABZ*I_NAI_Gx2yz_G4x_bb;
  Double I_NAI_Gxy2z_H4xz_bb = I_NAI_Hxy3z_G4x_bb+ABZ*I_NAI_Gxy2z_G4x_bb;
  Double I_NAI_Gx3z_H4xz_bb = I_NAI_Hx4z_G4x_bb+ABZ*I_NAI_Gx3z_G4x_bb;
  Double I_NAI_G4y_H4xz_bb = I_NAI_H4yz_G4x_bb+ABZ*I_NAI_G4y_G4x_bb;
  Double I_NAI_G3yz_H4xz_bb = I_NAI_H3y2z_G4x_bb+ABZ*I_NAI_G3yz_G4x_bb;
  Double I_NAI_G2y2z_H4xz_bb = I_NAI_H2y3z_G4x_bb+ABZ*I_NAI_G2y2z_G4x_bb;
  Double I_NAI_Gy3z_H4xz_bb = I_NAI_Hy4z_G4x_bb+ABZ*I_NAI_Gy3z_G4x_bb;
  Double I_NAI_G4z_H4xz_bb = I_NAI_H5z_G4x_bb+ABZ*I_NAI_G4z_G4x_bb;
  Double I_NAI_G4x_H3x2y_bb = I_NAI_H4xy_G3xy_bb+ABY*I_NAI_G4x_G3xy_bb;
  Double I_NAI_G3xy_H3x2y_bb = I_NAI_H3x2y_G3xy_bb+ABY*I_NAI_G3xy_G3xy_bb;
  Double I_NAI_G3xz_H3x2y_bb = I_NAI_H3xyz_G3xy_bb+ABY*I_NAI_G3xz_G3xy_bb;
  Double I_NAI_G2x2y_H3x2y_bb = I_NAI_H2x3y_G3xy_bb+ABY*I_NAI_G2x2y_G3xy_bb;
  Double I_NAI_G2xyz_H3x2y_bb = I_NAI_H2x2yz_G3xy_bb+ABY*I_NAI_G2xyz_G3xy_bb;
  Double I_NAI_G2x2z_H3x2y_bb = I_NAI_H2xy2z_G3xy_bb+ABY*I_NAI_G2x2z_G3xy_bb;
  Double I_NAI_Gx3y_H3x2y_bb = I_NAI_Hx4y_G3xy_bb+ABY*I_NAI_Gx3y_G3xy_bb;
  Double I_NAI_Gx2yz_H3x2y_bb = I_NAI_Hx3yz_G3xy_bb+ABY*I_NAI_Gx2yz_G3xy_bb;
  Double I_NAI_Gxy2z_H3x2y_bb = I_NAI_Hx2y2z_G3xy_bb+ABY*I_NAI_Gxy2z_G3xy_bb;
  Double I_NAI_Gx3z_H3x2y_bb = I_NAI_Hxy3z_G3xy_bb+ABY*I_NAI_Gx3z_G3xy_bb;
  Double I_NAI_G4y_H3x2y_bb = I_NAI_H5y_G3xy_bb+ABY*I_NAI_G4y_G3xy_bb;
  Double I_NAI_G3yz_H3x2y_bb = I_NAI_H4yz_G3xy_bb+ABY*I_NAI_G3yz_G3xy_bb;
  Double I_NAI_G2y2z_H3x2y_bb = I_NAI_H3y2z_G3xy_bb+ABY*I_NAI_G2y2z_G3xy_bb;
  Double I_NAI_Gy3z_H3x2y_bb = I_NAI_H2y3z_G3xy_bb+ABY*I_NAI_Gy3z_G3xy_bb;
  Double I_NAI_G4z_H3x2y_bb = I_NAI_Hy4z_G3xy_bb+ABY*I_NAI_G4z_G3xy_bb;
  Double I_NAI_G4x_H3xyz_bb = I_NAI_H4xz_G3xy_bb+ABZ*I_NAI_G4x_G3xy_bb;
  Double I_NAI_G3xy_H3xyz_bb = I_NAI_H3xyz_G3xy_bb+ABZ*I_NAI_G3xy_G3xy_bb;
  Double I_NAI_G3xz_H3xyz_bb = I_NAI_H3x2z_G3xy_bb+ABZ*I_NAI_G3xz_G3xy_bb;
  Double I_NAI_G2x2y_H3xyz_bb = I_NAI_H2x2yz_G3xy_bb+ABZ*I_NAI_G2x2y_G3xy_bb;
  Double I_NAI_G2xyz_H3xyz_bb = I_NAI_H2xy2z_G3xy_bb+ABZ*I_NAI_G2xyz_G3xy_bb;
  Double I_NAI_G2x2z_H3xyz_bb = I_NAI_H2x3z_G3xy_bb+ABZ*I_NAI_G2x2z_G3xy_bb;
  Double I_NAI_Gx3y_H3xyz_bb = I_NAI_Hx3yz_G3xy_bb+ABZ*I_NAI_Gx3y_G3xy_bb;
  Double I_NAI_Gx2yz_H3xyz_bb = I_NAI_Hx2y2z_G3xy_bb+ABZ*I_NAI_Gx2yz_G3xy_bb;
  Double I_NAI_Gxy2z_H3xyz_bb = I_NAI_Hxy3z_G3xy_bb+ABZ*I_NAI_Gxy2z_G3xy_bb;
  Double I_NAI_Gx3z_H3xyz_bb = I_NAI_Hx4z_G3xy_bb+ABZ*I_NAI_Gx3z_G3xy_bb;
  Double I_NAI_G4y_H3xyz_bb = I_NAI_H4yz_G3xy_bb+ABZ*I_NAI_G4y_G3xy_bb;
  Double I_NAI_G3yz_H3xyz_bb = I_NAI_H3y2z_G3xy_bb+ABZ*I_NAI_G3yz_G3xy_bb;
  Double I_NAI_G2y2z_H3xyz_bb = I_NAI_H2y3z_G3xy_bb+ABZ*I_NAI_G2y2z_G3xy_bb;
  Double I_NAI_Gy3z_H3xyz_bb = I_NAI_Hy4z_G3xy_bb+ABZ*I_NAI_Gy3z_G3xy_bb;
  Double I_NAI_G4z_H3xyz_bb = I_NAI_H5z_G3xy_bb+ABZ*I_NAI_G4z_G3xy_bb;
  Double I_NAI_G4x_H3x2z_bb = I_NAI_H4xz_G3xz_bb+ABZ*I_NAI_G4x_G3xz_bb;
  Double I_NAI_G3xy_H3x2z_bb = I_NAI_H3xyz_G3xz_bb+ABZ*I_NAI_G3xy_G3xz_bb;
  Double I_NAI_G3xz_H3x2z_bb = I_NAI_H3x2z_G3xz_bb+ABZ*I_NAI_G3xz_G3xz_bb;
  Double I_NAI_G2x2y_H3x2z_bb = I_NAI_H2x2yz_G3xz_bb+ABZ*I_NAI_G2x2y_G3xz_bb;
  Double I_NAI_G2xyz_H3x2z_bb = I_NAI_H2xy2z_G3xz_bb+ABZ*I_NAI_G2xyz_G3xz_bb;
  Double I_NAI_G2x2z_H3x2z_bb = I_NAI_H2x3z_G3xz_bb+ABZ*I_NAI_G2x2z_G3xz_bb;
  Double I_NAI_Gx3y_H3x2z_bb = I_NAI_Hx3yz_G3xz_bb+ABZ*I_NAI_Gx3y_G3xz_bb;
  Double I_NAI_Gx2yz_H3x2z_bb = I_NAI_Hx2y2z_G3xz_bb+ABZ*I_NAI_Gx2yz_G3xz_bb;
  Double I_NAI_Gxy2z_H3x2z_bb = I_NAI_Hxy3z_G3xz_bb+ABZ*I_NAI_Gxy2z_G3xz_bb;
  Double I_NAI_Gx3z_H3x2z_bb = I_NAI_Hx4z_G3xz_bb+ABZ*I_NAI_Gx3z_G3xz_bb;
  Double I_NAI_G4y_H3x2z_bb = I_NAI_H4yz_G3xz_bb+ABZ*I_NAI_G4y_G3xz_bb;
  Double I_NAI_G3yz_H3x2z_bb = I_NAI_H3y2z_G3xz_bb+ABZ*I_NAI_G3yz_G3xz_bb;
  Double I_NAI_G2y2z_H3x2z_bb = I_NAI_H2y3z_G3xz_bb+ABZ*I_NAI_G2y2z_G3xz_bb;
  Double I_NAI_Gy3z_H3x2z_bb = I_NAI_Hy4z_G3xz_bb+ABZ*I_NAI_Gy3z_G3xz_bb;
  Double I_NAI_G4z_H3x2z_bb = I_NAI_H5z_G3xz_bb+ABZ*I_NAI_G4z_G3xz_bb;
  Double I_NAI_G4x_H2x3y_bb = I_NAI_H5x_Gx3y_bb+ABX*I_NAI_G4x_Gx3y_bb;
  Double I_NAI_G3xy_H2x3y_bb = I_NAI_H4xy_Gx3y_bb+ABX*I_NAI_G3xy_Gx3y_bb;
  Double I_NAI_G3xz_H2x3y_bb = I_NAI_H4xz_Gx3y_bb+ABX*I_NAI_G3xz_Gx3y_bb;
  Double I_NAI_G2x2y_H2x3y_bb = I_NAI_H3x2y_Gx3y_bb+ABX*I_NAI_G2x2y_Gx3y_bb;
  Double I_NAI_G2xyz_H2x3y_bb = I_NAI_H3xyz_Gx3y_bb+ABX*I_NAI_G2xyz_Gx3y_bb;
  Double I_NAI_G2x2z_H2x3y_bb = I_NAI_H3x2z_Gx3y_bb+ABX*I_NAI_G2x2z_Gx3y_bb;
  Double I_NAI_Gx3y_H2x3y_bb = I_NAI_H2x3y_Gx3y_bb+ABX*I_NAI_Gx3y_Gx3y_bb;
  Double I_NAI_Gx2yz_H2x3y_bb = I_NAI_H2x2yz_Gx3y_bb+ABX*I_NAI_Gx2yz_Gx3y_bb;
  Double I_NAI_Gxy2z_H2x3y_bb = I_NAI_H2xy2z_Gx3y_bb+ABX*I_NAI_Gxy2z_Gx3y_bb;
  Double I_NAI_Gx3z_H2x3y_bb = I_NAI_H2x3z_Gx3y_bb+ABX*I_NAI_Gx3z_Gx3y_bb;
  Double I_NAI_G4y_H2x3y_bb = I_NAI_Hx4y_Gx3y_bb+ABX*I_NAI_G4y_Gx3y_bb;
  Double I_NAI_G3yz_H2x3y_bb = I_NAI_Hx3yz_Gx3y_bb+ABX*I_NAI_G3yz_Gx3y_bb;
  Double I_NAI_G2y2z_H2x3y_bb = I_NAI_Hx2y2z_Gx3y_bb+ABX*I_NAI_G2y2z_Gx3y_bb;
  Double I_NAI_Gy3z_H2x3y_bb = I_NAI_Hxy3z_Gx3y_bb+ABX*I_NAI_Gy3z_Gx3y_bb;
  Double I_NAI_G4z_H2x3y_bb = I_NAI_Hx4z_Gx3y_bb+ABX*I_NAI_G4z_Gx3y_bb;
  Double I_NAI_G4x_H2x2yz_bb = I_NAI_H4xz_G2x2y_bb+ABZ*I_NAI_G4x_G2x2y_bb;
  Double I_NAI_G3xy_H2x2yz_bb = I_NAI_H3xyz_G2x2y_bb+ABZ*I_NAI_G3xy_G2x2y_bb;
  Double I_NAI_G3xz_H2x2yz_bb = I_NAI_H3x2z_G2x2y_bb+ABZ*I_NAI_G3xz_G2x2y_bb;
  Double I_NAI_G2x2y_H2x2yz_bb = I_NAI_H2x2yz_G2x2y_bb+ABZ*I_NAI_G2x2y_G2x2y_bb;
  Double I_NAI_G2xyz_H2x2yz_bb = I_NAI_H2xy2z_G2x2y_bb+ABZ*I_NAI_G2xyz_G2x2y_bb;
  Double I_NAI_G2x2z_H2x2yz_bb = I_NAI_H2x3z_G2x2y_bb+ABZ*I_NAI_G2x2z_G2x2y_bb;
  Double I_NAI_Gx3y_H2x2yz_bb = I_NAI_Hx3yz_G2x2y_bb+ABZ*I_NAI_Gx3y_G2x2y_bb;
  Double I_NAI_Gx2yz_H2x2yz_bb = I_NAI_Hx2y2z_G2x2y_bb+ABZ*I_NAI_Gx2yz_G2x2y_bb;
  Double I_NAI_Gxy2z_H2x2yz_bb = I_NAI_Hxy3z_G2x2y_bb+ABZ*I_NAI_Gxy2z_G2x2y_bb;
  Double I_NAI_Gx3z_H2x2yz_bb = I_NAI_Hx4z_G2x2y_bb+ABZ*I_NAI_Gx3z_G2x2y_bb;
  Double I_NAI_G4y_H2x2yz_bb = I_NAI_H4yz_G2x2y_bb+ABZ*I_NAI_G4y_G2x2y_bb;
  Double I_NAI_G3yz_H2x2yz_bb = I_NAI_H3y2z_G2x2y_bb+ABZ*I_NAI_G3yz_G2x2y_bb;
  Double I_NAI_G2y2z_H2x2yz_bb = I_NAI_H2y3z_G2x2y_bb+ABZ*I_NAI_G2y2z_G2x2y_bb;
  Double I_NAI_Gy3z_H2x2yz_bb = I_NAI_Hy4z_G2x2y_bb+ABZ*I_NAI_Gy3z_G2x2y_bb;
  Double I_NAI_G4z_H2x2yz_bb = I_NAI_H5z_G2x2y_bb+ABZ*I_NAI_G4z_G2x2y_bb;
  Double I_NAI_G4x_H2xy2z_bb = I_NAI_H4xy_G2x2z_bb+ABY*I_NAI_G4x_G2x2z_bb;
  Double I_NAI_G3xy_H2xy2z_bb = I_NAI_H3x2y_G2x2z_bb+ABY*I_NAI_G3xy_G2x2z_bb;
  Double I_NAI_G3xz_H2xy2z_bb = I_NAI_H3xyz_G2x2z_bb+ABY*I_NAI_G3xz_G2x2z_bb;
  Double I_NAI_G2x2y_H2xy2z_bb = I_NAI_H2x3y_G2x2z_bb+ABY*I_NAI_G2x2y_G2x2z_bb;
  Double I_NAI_G2xyz_H2xy2z_bb = I_NAI_H2x2yz_G2x2z_bb+ABY*I_NAI_G2xyz_G2x2z_bb;
  Double I_NAI_G2x2z_H2xy2z_bb = I_NAI_H2xy2z_G2x2z_bb+ABY*I_NAI_G2x2z_G2x2z_bb;
  Double I_NAI_Gx3y_H2xy2z_bb = I_NAI_Hx4y_G2x2z_bb+ABY*I_NAI_Gx3y_G2x2z_bb;
  Double I_NAI_Gx2yz_H2xy2z_bb = I_NAI_Hx3yz_G2x2z_bb+ABY*I_NAI_Gx2yz_G2x2z_bb;
  Double I_NAI_Gxy2z_H2xy2z_bb = I_NAI_Hx2y2z_G2x2z_bb+ABY*I_NAI_Gxy2z_G2x2z_bb;
  Double I_NAI_Gx3z_H2xy2z_bb = I_NAI_Hxy3z_G2x2z_bb+ABY*I_NAI_Gx3z_G2x2z_bb;
  Double I_NAI_G4y_H2xy2z_bb = I_NAI_H5y_G2x2z_bb+ABY*I_NAI_G4y_G2x2z_bb;
  Double I_NAI_G3yz_H2xy2z_bb = I_NAI_H4yz_G2x2z_bb+ABY*I_NAI_G3yz_G2x2z_bb;
  Double I_NAI_G2y2z_H2xy2z_bb = I_NAI_H3y2z_G2x2z_bb+ABY*I_NAI_G2y2z_G2x2z_bb;
  Double I_NAI_Gy3z_H2xy2z_bb = I_NAI_H2y3z_G2x2z_bb+ABY*I_NAI_Gy3z_G2x2z_bb;
  Double I_NAI_G4z_H2xy2z_bb = I_NAI_Hy4z_G2x2z_bb+ABY*I_NAI_G4z_G2x2z_bb;
  Double I_NAI_G4x_H2x3z_bb = I_NAI_H5x_Gx3z_bb+ABX*I_NAI_G4x_Gx3z_bb;
  Double I_NAI_G3xy_H2x3z_bb = I_NAI_H4xy_Gx3z_bb+ABX*I_NAI_G3xy_Gx3z_bb;
  Double I_NAI_G3xz_H2x3z_bb = I_NAI_H4xz_Gx3z_bb+ABX*I_NAI_G3xz_Gx3z_bb;
  Double I_NAI_G2x2y_H2x3z_bb = I_NAI_H3x2y_Gx3z_bb+ABX*I_NAI_G2x2y_Gx3z_bb;
  Double I_NAI_G2xyz_H2x3z_bb = I_NAI_H3xyz_Gx3z_bb+ABX*I_NAI_G2xyz_Gx3z_bb;
  Double I_NAI_G2x2z_H2x3z_bb = I_NAI_H3x2z_Gx3z_bb+ABX*I_NAI_G2x2z_Gx3z_bb;
  Double I_NAI_Gx3y_H2x3z_bb = I_NAI_H2x3y_Gx3z_bb+ABX*I_NAI_Gx3y_Gx3z_bb;
  Double I_NAI_Gx2yz_H2x3z_bb = I_NAI_H2x2yz_Gx3z_bb+ABX*I_NAI_Gx2yz_Gx3z_bb;
  Double I_NAI_Gxy2z_H2x3z_bb = I_NAI_H2xy2z_Gx3z_bb+ABX*I_NAI_Gxy2z_Gx3z_bb;
  Double I_NAI_Gx3z_H2x3z_bb = I_NAI_H2x3z_Gx3z_bb+ABX*I_NAI_Gx3z_Gx3z_bb;
  Double I_NAI_G4y_H2x3z_bb = I_NAI_Hx4y_Gx3z_bb+ABX*I_NAI_G4y_Gx3z_bb;
  Double I_NAI_G3yz_H2x3z_bb = I_NAI_Hx3yz_Gx3z_bb+ABX*I_NAI_G3yz_Gx3z_bb;
  Double I_NAI_G2y2z_H2x3z_bb = I_NAI_Hx2y2z_Gx3z_bb+ABX*I_NAI_G2y2z_Gx3z_bb;
  Double I_NAI_Gy3z_H2x3z_bb = I_NAI_Hxy3z_Gx3z_bb+ABX*I_NAI_Gy3z_Gx3z_bb;
  Double I_NAI_G4z_H2x3z_bb = I_NAI_Hx4z_Gx3z_bb+ABX*I_NAI_G4z_Gx3z_bb;
  Double I_NAI_G4x_Hx4y_bb = I_NAI_H5x_G4y_bb+ABX*I_NAI_G4x_G4y_bb;
  Double I_NAI_G3xy_Hx4y_bb = I_NAI_H4xy_G4y_bb+ABX*I_NAI_G3xy_G4y_bb;
  Double I_NAI_G3xz_Hx4y_bb = I_NAI_H4xz_G4y_bb+ABX*I_NAI_G3xz_G4y_bb;
  Double I_NAI_G2x2y_Hx4y_bb = I_NAI_H3x2y_G4y_bb+ABX*I_NAI_G2x2y_G4y_bb;
  Double I_NAI_G2xyz_Hx4y_bb = I_NAI_H3xyz_G4y_bb+ABX*I_NAI_G2xyz_G4y_bb;
  Double I_NAI_G2x2z_Hx4y_bb = I_NAI_H3x2z_G4y_bb+ABX*I_NAI_G2x2z_G4y_bb;
  Double I_NAI_Gx3y_Hx4y_bb = I_NAI_H2x3y_G4y_bb+ABX*I_NAI_Gx3y_G4y_bb;
  Double I_NAI_Gx2yz_Hx4y_bb = I_NAI_H2x2yz_G4y_bb+ABX*I_NAI_Gx2yz_G4y_bb;
  Double I_NAI_Gxy2z_Hx4y_bb = I_NAI_H2xy2z_G4y_bb+ABX*I_NAI_Gxy2z_G4y_bb;
  Double I_NAI_Gx3z_Hx4y_bb = I_NAI_H2x3z_G4y_bb+ABX*I_NAI_Gx3z_G4y_bb;
  Double I_NAI_G4y_Hx4y_bb = I_NAI_Hx4y_G4y_bb+ABX*I_NAI_G4y_G4y_bb;
  Double I_NAI_G3yz_Hx4y_bb = I_NAI_Hx3yz_G4y_bb+ABX*I_NAI_G3yz_G4y_bb;
  Double I_NAI_G2y2z_Hx4y_bb = I_NAI_Hx2y2z_G4y_bb+ABX*I_NAI_G2y2z_G4y_bb;
  Double I_NAI_Gy3z_Hx4y_bb = I_NAI_Hxy3z_G4y_bb+ABX*I_NAI_Gy3z_G4y_bb;
  Double I_NAI_G4z_Hx4y_bb = I_NAI_Hx4z_G4y_bb+ABX*I_NAI_G4z_G4y_bb;
  Double I_NAI_G4x_Hx3yz_bb = I_NAI_H4xz_Gx3y_bb+ABZ*I_NAI_G4x_Gx3y_bb;
  Double I_NAI_G3xy_Hx3yz_bb = I_NAI_H3xyz_Gx3y_bb+ABZ*I_NAI_G3xy_Gx3y_bb;
  Double I_NAI_G3xz_Hx3yz_bb = I_NAI_H3x2z_Gx3y_bb+ABZ*I_NAI_G3xz_Gx3y_bb;
  Double I_NAI_G2x2y_Hx3yz_bb = I_NAI_H2x2yz_Gx3y_bb+ABZ*I_NAI_G2x2y_Gx3y_bb;
  Double I_NAI_G2xyz_Hx3yz_bb = I_NAI_H2xy2z_Gx3y_bb+ABZ*I_NAI_G2xyz_Gx3y_bb;
  Double I_NAI_G2x2z_Hx3yz_bb = I_NAI_H2x3z_Gx3y_bb+ABZ*I_NAI_G2x2z_Gx3y_bb;
  Double I_NAI_Gx3y_Hx3yz_bb = I_NAI_Hx3yz_Gx3y_bb+ABZ*I_NAI_Gx3y_Gx3y_bb;
  Double I_NAI_Gx2yz_Hx3yz_bb = I_NAI_Hx2y2z_Gx3y_bb+ABZ*I_NAI_Gx2yz_Gx3y_bb;
  Double I_NAI_Gxy2z_Hx3yz_bb = I_NAI_Hxy3z_Gx3y_bb+ABZ*I_NAI_Gxy2z_Gx3y_bb;
  Double I_NAI_Gx3z_Hx3yz_bb = I_NAI_Hx4z_Gx3y_bb+ABZ*I_NAI_Gx3z_Gx3y_bb;
  Double I_NAI_G4y_Hx3yz_bb = I_NAI_H4yz_Gx3y_bb+ABZ*I_NAI_G4y_Gx3y_bb;
  Double I_NAI_G3yz_Hx3yz_bb = I_NAI_H3y2z_Gx3y_bb+ABZ*I_NAI_G3yz_Gx3y_bb;
  Double I_NAI_G2y2z_Hx3yz_bb = I_NAI_H2y3z_Gx3y_bb+ABZ*I_NAI_G2y2z_Gx3y_bb;
  Double I_NAI_Gy3z_Hx3yz_bb = I_NAI_Hy4z_Gx3y_bb+ABZ*I_NAI_Gy3z_Gx3y_bb;
  Double I_NAI_G4z_Hx3yz_bb = I_NAI_H5z_Gx3y_bb+ABZ*I_NAI_G4z_Gx3y_bb;
  Double I_NAI_G4x_Hx2y2z_bb = I_NAI_H5x_G2y2z_bb+ABX*I_NAI_G4x_G2y2z_bb;
  Double I_NAI_G3xy_Hx2y2z_bb = I_NAI_H4xy_G2y2z_bb+ABX*I_NAI_G3xy_G2y2z_bb;
  Double I_NAI_G3xz_Hx2y2z_bb = I_NAI_H4xz_G2y2z_bb+ABX*I_NAI_G3xz_G2y2z_bb;
  Double I_NAI_G2x2y_Hx2y2z_bb = I_NAI_H3x2y_G2y2z_bb+ABX*I_NAI_G2x2y_G2y2z_bb;
  Double I_NAI_G2xyz_Hx2y2z_bb = I_NAI_H3xyz_G2y2z_bb+ABX*I_NAI_G2xyz_G2y2z_bb;
  Double I_NAI_G2x2z_Hx2y2z_bb = I_NAI_H3x2z_G2y2z_bb+ABX*I_NAI_G2x2z_G2y2z_bb;
  Double I_NAI_Gx3y_Hx2y2z_bb = I_NAI_H2x3y_G2y2z_bb+ABX*I_NAI_Gx3y_G2y2z_bb;
  Double I_NAI_Gx2yz_Hx2y2z_bb = I_NAI_H2x2yz_G2y2z_bb+ABX*I_NAI_Gx2yz_G2y2z_bb;
  Double I_NAI_Gxy2z_Hx2y2z_bb = I_NAI_H2xy2z_G2y2z_bb+ABX*I_NAI_Gxy2z_G2y2z_bb;
  Double I_NAI_Gx3z_Hx2y2z_bb = I_NAI_H2x3z_G2y2z_bb+ABX*I_NAI_Gx3z_G2y2z_bb;
  Double I_NAI_G4y_Hx2y2z_bb = I_NAI_Hx4y_G2y2z_bb+ABX*I_NAI_G4y_G2y2z_bb;
  Double I_NAI_G3yz_Hx2y2z_bb = I_NAI_Hx3yz_G2y2z_bb+ABX*I_NAI_G3yz_G2y2z_bb;
  Double I_NAI_G2y2z_Hx2y2z_bb = I_NAI_Hx2y2z_G2y2z_bb+ABX*I_NAI_G2y2z_G2y2z_bb;
  Double I_NAI_Gy3z_Hx2y2z_bb = I_NAI_Hxy3z_G2y2z_bb+ABX*I_NAI_Gy3z_G2y2z_bb;
  Double I_NAI_G4z_Hx2y2z_bb = I_NAI_Hx4z_G2y2z_bb+ABX*I_NAI_G4z_G2y2z_bb;
  Double I_NAI_G4x_Hxy3z_bb = I_NAI_H4xy_Gx3z_bb+ABY*I_NAI_G4x_Gx3z_bb;
  Double I_NAI_G3xy_Hxy3z_bb = I_NAI_H3x2y_Gx3z_bb+ABY*I_NAI_G3xy_Gx3z_bb;
  Double I_NAI_G3xz_Hxy3z_bb = I_NAI_H3xyz_Gx3z_bb+ABY*I_NAI_G3xz_Gx3z_bb;
  Double I_NAI_G2x2y_Hxy3z_bb = I_NAI_H2x3y_Gx3z_bb+ABY*I_NAI_G2x2y_Gx3z_bb;
  Double I_NAI_G2xyz_Hxy3z_bb = I_NAI_H2x2yz_Gx3z_bb+ABY*I_NAI_G2xyz_Gx3z_bb;
  Double I_NAI_G2x2z_Hxy3z_bb = I_NAI_H2xy2z_Gx3z_bb+ABY*I_NAI_G2x2z_Gx3z_bb;
  Double I_NAI_Gx3y_Hxy3z_bb = I_NAI_Hx4y_Gx3z_bb+ABY*I_NAI_Gx3y_Gx3z_bb;
  Double I_NAI_Gx2yz_Hxy3z_bb = I_NAI_Hx3yz_Gx3z_bb+ABY*I_NAI_Gx2yz_Gx3z_bb;
  Double I_NAI_Gxy2z_Hxy3z_bb = I_NAI_Hx2y2z_Gx3z_bb+ABY*I_NAI_Gxy2z_Gx3z_bb;
  Double I_NAI_Gx3z_Hxy3z_bb = I_NAI_Hxy3z_Gx3z_bb+ABY*I_NAI_Gx3z_Gx3z_bb;
  Double I_NAI_G4y_Hxy3z_bb = I_NAI_H5y_Gx3z_bb+ABY*I_NAI_G4y_Gx3z_bb;
  Double I_NAI_G3yz_Hxy3z_bb = I_NAI_H4yz_Gx3z_bb+ABY*I_NAI_G3yz_Gx3z_bb;
  Double I_NAI_G2y2z_Hxy3z_bb = I_NAI_H3y2z_Gx3z_bb+ABY*I_NAI_G2y2z_Gx3z_bb;
  Double I_NAI_Gy3z_Hxy3z_bb = I_NAI_H2y3z_Gx3z_bb+ABY*I_NAI_Gy3z_Gx3z_bb;
  Double I_NAI_G4z_Hxy3z_bb = I_NAI_Hy4z_Gx3z_bb+ABY*I_NAI_G4z_Gx3z_bb;
  Double I_NAI_G4x_Hx4z_bb = I_NAI_H5x_G4z_bb+ABX*I_NAI_G4x_G4z_bb;
  Double I_NAI_G3xy_Hx4z_bb = I_NAI_H4xy_G4z_bb+ABX*I_NAI_G3xy_G4z_bb;
  Double I_NAI_G3xz_Hx4z_bb = I_NAI_H4xz_G4z_bb+ABX*I_NAI_G3xz_G4z_bb;
  Double I_NAI_G2x2y_Hx4z_bb = I_NAI_H3x2y_G4z_bb+ABX*I_NAI_G2x2y_G4z_bb;
  Double I_NAI_G2xyz_Hx4z_bb = I_NAI_H3xyz_G4z_bb+ABX*I_NAI_G2xyz_G4z_bb;
  Double I_NAI_G2x2z_Hx4z_bb = I_NAI_H3x2z_G4z_bb+ABX*I_NAI_G2x2z_G4z_bb;
  Double I_NAI_Gx3y_Hx4z_bb = I_NAI_H2x3y_G4z_bb+ABX*I_NAI_Gx3y_G4z_bb;
  Double I_NAI_Gx2yz_Hx4z_bb = I_NAI_H2x2yz_G4z_bb+ABX*I_NAI_Gx2yz_G4z_bb;
  Double I_NAI_Gxy2z_Hx4z_bb = I_NAI_H2xy2z_G4z_bb+ABX*I_NAI_Gxy2z_G4z_bb;
  Double I_NAI_Gx3z_Hx4z_bb = I_NAI_H2x3z_G4z_bb+ABX*I_NAI_Gx3z_G4z_bb;
  Double I_NAI_G4y_Hx4z_bb = I_NAI_Hx4y_G4z_bb+ABX*I_NAI_G4y_G4z_bb;
  Double I_NAI_G3yz_Hx4z_bb = I_NAI_Hx3yz_G4z_bb+ABX*I_NAI_G3yz_G4z_bb;
  Double I_NAI_G2y2z_Hx4z_bb = I_NAI_Hx2y2z_G4z_bb+ABX*I_NAI_G2y2z_G4z_bb;
  Double I_NAI_Gy3z_Hx4z_bb = I_NAI_Hxy3z_G4z_bb+ABX*I_NAI_Gy3z_G4z_bb;
  Double I_NAI_G4z_Hx4z_bb = I_NAI_Hx4z_G4z_bb+ABX*I_NAI_G4z_G4z_bb;
  Double I_NAI_G4x_H5y_bb = I_NAI_H4xy_G4y_bb+ABY*I_NAI_G4x_G4y_bb;
  Double I_NAI_G3xy_H5y_bb = I_NAI_H3x2y_G4y_bb+ABY*I_NAI_G3xy_G4y_bb;
  Double I_NAI_G3xz_H5y_bb = I_NAI_H3xyz_G4y_bb+ABY*I_NAI_G3xz_G4y_bb;
  Double I_NAI_G2x2y_H5y_bb = I_NAI_H2x3y_G4y_bb+ABY*I_NAI_G2x2y_G4y_bb;
  Double I_NAI_G2xyz_H5y_bb = I_NAI_H2x2yz_G4y_bb+ABY*I_NAI_G2xyz_G4y_bb;
  Double I_NAI_G2x2z_H5y_bb = I_NAI_H2xy2z_G4y_bb+ABY*I_NAI_G2x2z_G4y_bb;
  Double I_NAI_Gx3y_H5y_bb = I_NAI_Hx4y_G4y_bb+ABY*I_NAI_Gx3y_G4y_bb;
  Double I_NAI_Gx2yz_H5y_bb = I_NAI_Hx3yz_G4y_bb+ABY*I_NAI_Gx2yz_G4y_bb;
  Double I_NAI_Gxy2z_H5y_bb = I_NAI_Hx2y2z_G4y_bb+ABY*I_NAI_Gxy2z_G4y_bb;
  Double I_NAI_Gx3z_H5y_bb = I_NAI_Hxy3z_G4y_bb+ABY*I_NAI_Gx3z_G4y_bb;
  Double I_NAI_G4y_H5y_bb = I_NAI_H5y_G4y_bb+ABY*I_NAI_G4y_G4y_bb;
  Double I_NAI_G3yz_H5y_bb = I_NAI_H4yz_G4y_bb+ABY*I_NAI_G3yz_G4y_bb;
  Double I_NAI_G2y2z_H5y_bb = I_NAI_H3y2z_G4y_bb+ABY*I_NAI_G2y2z_G4y_bb;
  Double I_NAI_Gy3z_H5y_bb = I_NAI_H2y3z_G4y_bb+ABY*I_NAI_Gy3z_G4y_bb;
  Double I_NAI_G4z_H5y_bb = I_NAI_Hy4z_G4y_bb+ABY*I_NAI_G4z_G4y_bb;
  Double I_NAI_G4x_H4yz_bb = I_NAI_H4xz_G4y_bb+ABZ*I_NAI_G4x_G4y_bb;
  Double I_NAI_G3xy_H4yz_bb = I_NAI_H3xyz_G4y_bb+ABZ*I_NAI_G3xy_G4y_bb;
  Double I_NAI_G3xz_H4yz_bb = I_NAI_H3x2z_G4y_bb+ABZ*I_NAI_G3xz_G4y_bb;
  Double I_NAI_G2x2y_H4yz_bb = I_NAI_H2x2yz_G4y_bb+ABZ*I_NAI_G2x2y_G4y_bb;
  Double I_NAI_G2xyz_H4yz_bb = I_NAI_H2xy2z_G4y_bb+ABZ*I_NAI_G2xyz_G4y_bb;
  Double I_NAI_G2x2z_H4yz_bb = I_NAI_H2x3z_G4y_bb+ABZ*I_NAI_G2x2z_G4y_bb;
  Double I_NAI_Gx3y_H4yz_bb = I_NAI_Hx3yz_G4y_bb+ABZ*I_NAI_Gx3y_G4y_bb;
  Double I_NAI_Gx2yz_H4yz_bb = I_NAI_Hx2y2z_G4y_bb+ABZ*I_NAI_Gx2yz_G4y_bb;
  Double I_NAI_Gxy2z_H4yz_bb = I_NAI_Hxy3z_G4y_bb+ABZ*I_NAI_Gxy2z_G4y_bb;
  Double I_NAI_Gx3z_H4yz_bb = I_NAI_Hx4z_G4y_bb+ABZ*I_NAI_Gx3z_G4y_bb;
  Double I_NAI_G4y_H4yz_bb = I_NAI_H4yz_G4y_bb+ABZ*I_NAI_G4y_G4y_bb;
  Double I_NAI_G3yz_H4yz_bb = I_NAI_H3y2z_G4y_bb+ABZ*I_NAI_G3yz_G4y_bb;
  Double I_NAI_G2y2z_H4yz_bb = I_NAI_H2y3z_G4y_bb+ABZ*I_NAI_G2y2z_G4y_bb;
  Double I_NAI_Gy3z_H4yz_bb = I_NAI_Hy4z_G4y_bb+ABZ*I_NAI_Gy3z_G4y_bb;
  Double I_NAI_G4z_H4yz_bb = I_NAI_H5z_G4y_bb+ABZ*I_NAI_G4z_G4y_bb;
  Double I_NAI_G4x_H3y2z_bb = I_NAI_H4xz_G3yz_bb+ABZ*I_NAI_G4x_G3yz_bb;
  Double I_NAI_G3xy_H3y2z_bb = I_NAI_H3xyz_G3yz_bb+ABZ*I_NAI_G3xy_G3yz_bb;
  Double I_NAI_G3xz_H3y2z_bb = I_NAI_H3x2z_G3yz_bb+ABZ*I_NAI_G3xz_G3yz_bb;
  Double I_NAI_G2x2y_H3y2z_bb = I_NAI_H2x2yz_G3yz_bb+ABZ*I_NAI_G2x2y_G3yz_bb;
  Double I_NAI_G2xyz_H3y2z_bb = I_NAI_H2xy2z_G3yz_bb+ABZ*I_NAI_G2xyz_G3yz_bb;
  Double I_NAI_G2x2z_H3y2z_bb = I_NAI_H2x3z_G3yz_bb+ABZ*I_NAI_G2x2z_G3yz_bb;
  Double I_NAI_Gx3y_H3y2z_bb = I_NAI_Hx3yz_G3yz_bb+ABZ*I_NAI_Gx3y_G3yz_bb;
  Double I_NAI_Gx2yz_H3y2z_bb = I_NAI_Hx2y2z_G3yz_bb+ABZ*I_NAI_Gx2yz_G3yz_bb;
  Double I_NAI_Gxy2z_H3y2z_bb = I_NAI_Hxy3z_G3yz_bb+ABZ*I_NAI_Gxy2z_G3yz_bb;
  Double I_NAI_Gx3z_H3y2z_bb = I_NAI_Hx4z_G3yz_bb+ABZ*I_NAI_Gx3z_G3yz_bb;
  Double I_NAI_G4y_H3y2z_bb = I_NAI_H4yz_G3yz_bb+ABZ*I_NAI_G4y_G3yz_bb;
  Double I_NAI_G3yz_H3y2z_bb = I_NAI_H3y2z_G3yz_bb+ABZ*I_NAI_G3yz_G3yz_bb;
  Double I_NAI_G2y2z_H3y2z_bb = I_NAI_H2y3z_G3yz_bb+ABZ*I_NAI_G2y2z_G3yz_bb;
  Double I_NAI_Gy3z_H3y2z_bb = I_NAI_Hy4z_G3yz_bb+ABZ*I_NAI_Gy3z_G3yz_bb;
  Double I_NAI_G4z_H3y2z_bb = I_NAI_H5z_G3yz_bb+ABZ*I_NAI_G4z_G3yz_bb;
  Double I_NAI_G4x_H2y3z_bb = I_NAI_H4xy_Gy3z_bb+ABY*I_NAI_G4x_Gy3z_bb;
  Double I_NAI_G3xy_H2y3z_bb = I_NAI_H3x2y_Gy3z_bb+ABY*I_NAI_G3xy_Gy3z_bb;
  Double I_NAI_G3xz_H2y3z_bb = I_NAI_H3xyz_Gy3z_bb+ABY*I_NAI_G3xz_Gy3z_bb;
  Double I_NAI_G2x2y_H2y3z_bb = I_NAI_H2x3y_Gy3z_bb+ABY*I_NAI_G2x2y_Gy3z_bb;
  Double I_NAI_G2xyz_H2y3z_bb = I_NAI_H2x2yz_Gy3z_bb+ABY*I_NAI_G2xyz_Gy3z_bb;
  Double I_NAI_G2x2z_H2y3z_bb = I_NAI_H2xy2z_Gy3z_bb+ABY*I_NAI_G2x2z_Gy3z_bb;
  Double I_NAI_Gx3y_H2y3z_bb = I_NAI_Hx4y_Gy3z_bb+ABY*I_NAI_Gx3y_Gy3z_bb;
  Double I_NAI_Gx2yz_H2y3z_bb = I_NAI_Hx3yz_Gy3z_bb+ABY*I_NAI_Gx2yz_Gy3z_bb;
  Double I_NAI_Gxy2z_H2y3z_bb = I_NAI_Hx2y2z_Gy3z_bb+ABY*I_NAI_Gxy2z_Gy3z_bb;
  Double I_NAI_Gx3z_H2y3z_bb = I_NAI_Hxy3z_Gy3z_bb+ABY*I_NAI_Gx3z_Gy3z_bb;
  Double I_NAI_G4y_H2y3z_bb = I_NAI_H5y_Gy3z_bb+ABY*I_NAI_G4y_Gy3z_bb;
  Double I_NAI_G3yz_H2y3z_bb = I_NAI_H4yz_Gy3z_bb+ABY*I_NAI_G3yz_Gy3z_bb;
  Double I_NAI_G2y2z_H2y3z_bb = I_NAI_H3y2z_Gy3z_bb+ABY*I_NAI_G2y2z_Gy3z_bb;
  Double I_NAI_Gy3z_H2y3z_bb = I_NAI_H2y3z_Gy3z_bb+ABY*I_NAI_Gy3z_Gy3z_bb;
  Double I_NAI_G4z_H2y3z_bb = I_NAI_Hy4z_Gy3z_bb+ABY*I_NAI_G4z_Gy3z_bb;
  Double I_NAI_G4x_Hy4z_bb = I_NAI_H4xy_G4z_bb+ABY*I_NAI_G4x_G4z_bb;
  Double I_NAI_G3xy_Hy4z_bb = I_NAI_H3x2y_G4z_bb+ABY*I_NAI_G3xy_G4z_bb;
  Double I_NAI_G3xz_Hy4z_bb = I_NAI_H3xyz_G4z_bb+ABY*I_NAI_G3xz_G4z_bb;
  Double I_NAI_G2x2y_Hy4z_bb = I_NAI_H2x3y_G4z_bb+ABY*I_NAI_G2x2y_G4z_bb;
  Double I_NAI_G2xyz_Hy4z_bb = I_NAI_H2x2yz_G4z_bb+ABY*I_NAI_G2xyz_G4z_bb;
  Double I_NAI_G2x2z_Hy4z_bb = I_NAI_H2xy2z_G4z_bb+ABY*I_NAI_G2x2z_G4z_bb;
  Double I_NAI_Gx3y_Hy4z_bb = I_NAI_Hx4y_G4z_bb+ABY*I_NAI_Gx3y_G4z_bb;
  Double I_NAI_Gx2yz_Hy4z_bb = I_NAI_Hx3yz_G4z_bb+ABY*I_NAI_Gx2yz_G4z_bb;
  Double I_NAI_Gxy2z_Hy4z_bb = I_NAI_Hx2y2z_G4z_bb+ABY*I_NAI_Gxy2z_G4z_bb;
  Double I_NAI_Gx3z_Hy4z_bb = I_NAI_Hxy3z_G4z_bb+ABY*I_NAI_Gx3z_G4z_bb;
  Double I_NAI_G4y_Hy4z_bb = I_NAI_H5y_G4z_bb+ABY*I_NAI_G4y_G4z_bb;
  Double I_NAI_G3yz_Hy4z_bb = I_NAI_H4yz_G4z_bb+ABY*I_NAI_G3yz_G4z_bb;
  Double I_NAI_G2y2z_Hy4z_bb = I_NAI_H3y2z_G4z_bb+ABY*I_NAI_G2y2z_G4z_bb;
  Double I_NAI_Gy3z_Hy4z_bb = I_NAI_H2y3z_G4z_bb+ABY*I_NAI_Gy3z_G4z_bb;
  Double I_NAI_G4z_Hy4z_bb = I_NAI_Hy4z_G4z_bb+ABY*I_NAI_G4z_G4z_bb;
  Double I_NAI_G4x_H5z_bb = I_NAI_H4xz_G4z_bb+ABZ*I_NAI_G4x_G4z_bb;
  Double I_NAI_G3xy_H5z_bb = I_NAI_H3xyz_G4z_bb+ABZ*I_NAI_G3xy_G4z_bb;
  Double I_NAI_G3xz_H5z_bb = I_NAI_H3x2z_G4z_bb+ABZ*I_NAI_G3xz_G4z_bb;
  Double I_NAI_G2x2y_H5z_bb = I_NAI_H2x2yz_G4z_bb+ABZ*I_NAI_G2x2y_G4z_bb;
  Double I_NAI_G2xyz_H5z_bb = I_NAI_H2xy2z_G4z_bb+ABZ*I_NAI_G2xyz_G4z_bb;
  Double I_NAI_G2x2z_H5z_bb = I_NAI_H2x3z_G4z_bb+ABZ*I_NAI_G2x2z_G4z_bb;
  Double I_NAI_Gx3y_H5z_bb = I_NAI_Hx3yz_G4z_bb+ABZ*I_NAI_Gx3y_G4z_bb;
  Double I_NAI_Gx2yz_H5z_bb = I_NAI_Hx2y2z_G4z_bb+ABZ*I_NAI_Gx2yz_G4z_bb;
  Double I_NAI_Gxy2z_H5z_bb = I_NAI_Hxy3z_G4z_bb+ABZ*I_NAI_Gxy2z_G4z_bb;
  Double I_NAI_Gx3z_H5z_bb = I_NAI_Hx4z_G4z_bb+ABZ*I_NAI_Gx3z_G4z_bb;
  Double I_NAI_G4y_H5z_bb = I_NAI_H4yz_G4z_bb+ABZ*I_NAI_G4y_G4z_bb;
  Double I_NAI_G3yz_H5z_bb = I_NAI_H3y2z_G4z_bb+ABZ*I_NAI_G3yz_G4z_bb;
  Double I_NAI_G2y2z_H5z_bb = I_NAI_H2y3z_G4z_bb+ABZ*I_NAI_G2y2z_G4z_bb;
  Double I_NAI_Gy3z_H5z_bb = I_NAI_Hy4z_G4z_bb+ABZ*I_NAI_Gy3z_G4z_bb;
  Double I_NAI_G4z_H5z_bb = I_NAI_H5z_G4z_bb+ABZ*I_NAI_G4z_G4z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[0] = 4.0E0*I_NAI_I6x_F3x_aa-2.0E0*4*I_NAI_G4x_F3x_a-2.0E0*5*I_NAI_G4x_F3x_a+4*3*I_NAI_D2x_F3x;
  abcd[1] = 4.0E0*I_NAI_I5xy_F3x_aa-2.0E0*3*I_NAI_G3xy_F3x_a-2.0E0*4*I_NAI_G3xy_F3x_a+3*2*I_NAI_Dxy_F3x;
  abcd[2] = 4.0E0*I_NAI_I5xz_F3x_aa-2.0E0*3*I_NAI_G3xz_F3x_a-2.0E0*4*I_NAI_G3xz_F3x_a+3*2*I_NAI_Dxz_F3x;
  abcd[3] = 4.0E0*I_NAI_I4x2y_F3x_aa-2.0E0*2*I_NAI_G2x2y_F3x_a-2.0E0*3*I_NAI_G2x2y_F3x_a+2*1*I_NAI_D2y_F3x;
  abcd[4] = 4.0E0*I_NAI_I4xyz_F3x_aa-2.0E0*2*I_NAI_G2xyz_F3x_a-2.0E0*3*I_NAI_G2xyz_F3x_a+2*1*I_NAI_Dyz_F3x;
  abcd[5] = 4.0E0*I_NAI_I4x2z_F3x_aa-2.0E0*2*I_NAI_G2x2z_F3x_a-2.0E0*3*I_NAI_G2x2z_F3x_a+2*1*I_NAI_D2z_F3x;
  abcd[6] = 4.0E0*I_NAI_I3x3y_F3x_aa-2.0E0*1*I_NAI_Gx3y_F3x_a-2.0E0*2*I_NAI_Gx3y_F3x_a;
  abcd[7] = 4.0E0*I_NAI_I3x2yz_F3x_aa-2.0E0*1*I_NAI_Gx2yz_F3x_a-2.0E0*2*I_NAI_Gx2yz_F3x_a;
  abcd[8] = 4.0E0*I_NAI_I3xy2z_F3x_aa-2.0E0*1*I_NAI_Gxy2z_F3x_a-2.0E0*2*I_NAI_Gxy2z_F3x_a;
  abcd[9] = 4.0E0*I_NAI_I3x3z_F3x_aa-2.0E0*1*I_NAI_Gx3z_F3x_a-2.0E0*2*I_NAI_Gx3z_F3x_a;
  abcd[10] = 4.0E0*I_NAI_I2x4y_F3x_aa-2.0E0*1*I_NAI_G4y_F3x_a;
  abcd[11] = 4.0E0*I_NAI_I2x3yz_F3x_aa-2.0E0*1*I_NAI_G3yz_F3x_a;
  abcd[12] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2y2z_F3x_a;
  abcd[13] = 4.0E0*I_NAI_I2xy3z_F3x_aa-2.0E0*1*I_NAI_Gy3z_F3x_a;
  abcd[14] = 4.0E0*I_NAI_I2x4z_F3x_aa-2.0E0*1*I_NAI_G4z_F3x_a;
  abcd[15] = 4.0E0*I_NAI_I6x_F2xy_aa-2.0E0*4*I_NAI_G4x_F2xy_a-2.0E0*5*I_NAI_G4x_F2xy_a+4*3*I_NAI_D2x_F2xy;
  abcd[16] = 4.0E0*I_NAI_I5xy_F2xy_aa-2.0E0*3*I_NAI_G3xy_F2xy_a-2.0E0*4*I_NAI_G3xy_F2xy_a+3*2*I_NAI_Dxy_F2xy;
  abcd[17] = 4.0E0*I_NAI_I5xz_F2xy_aa-2.0E0*3*I_NAI_G3xz_F2xy_a-2.0E0*4*I_NAI_G3xz_F2xy_a+3*2*I_NAI_Dxz_F2xy;
  abcd[18] = 4.0E0*I_NAI_I4x2y_F2xy_aa-2.0E0*2*I_NAI_G2x2y_F2xy_a-2.0E0*3*I_NAI_G2x2y_F2xy_a+2*1*I_NAI_D2y_F2xy;
  abcd[19] = 4.0E0*I_NAI_I4xyz_F2xy_aa-2.0E0*2*I_NAI_G2xyz_F2xy_a-2.0E0*3*I_NAI_G2xyz_F2xy_a+2*1*I_NAI_Dyz_F2xy;
  abcd[20] = 4.0E0*I_NAI_I4x2z_F2xy_aa-2.0E0*2*I_NAI_G2x2z_F2xy_a-2.0E0*3*I_NAI_G2x2z_F2xy_a+2*1*I_NAI_D2z_F2xy;
  abcd[21] = 4.0E0*I_NAI_I3x3y_F2xy_aa-2.0E0*1*I_NAI_Gx3y_F2xy_a-2.0E0*2*I_NAI_Gx3y_F2xy_a;
  abcd[22] = 4.0E0*I_NAI_I3x2yz_F2xy_aa-2.0E0*1*I_NAI_Gx2yz_F2xy_a-2.0E0*2*I_NAI_Gx2yz_F2xy_a;
  abcd[23] = 4.0E0*I_NAI_I3xy2z_F2xy_aa-2.0E0*1*I_NAI_Gxy2z_F2xy_a-2.0E0*2*I_NAI_Gxy2z_F2xy_a;
  abcd[24] = 4.0E0*I_NAI_I3x3z_F2xy_aa-2.0E0*1*I_NAI_Gx3z_F2xy_a-2.0E0*2*I_NAI_Gx3z_F2xy_a;
  abcd[25] = 4.0E0*I_NAI_I2x4y_F2xy_aa-2.0E0*1*I_NAI_G4y_F2xy_a;
  abcd[26] = 4.0E0*I_NAI_I2x3yz_F2xy_aa-2.0E0*1*I_NAI_G3yz_F2xy_a;
  abcd[27] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2y2z_F2xy_a;
  abcd[28] = 4.0E0*I_NAI_I2xy3z_F2xy_aa-2.0E0*1*I_NAI_Gy3z_F2xy_a;
  abcd[29] = 4.0E0*I_NAI_I2x4z_F2xy_aa-2.0E0*1*I_NAI_G4z_F2xy_a;
  abcd[30] = 4.0E0*I_NAI_I6x_F2xz_aa-2.0E0*4*I_NAI_G4x_F2xz_a-2.0E0*5*I_NAI_G4x_F2xz_a+4*3*I_NAI_D2x_F2xz;
  abcd[31] = 4.0E0*I_NAI_I5xy_F2xz_aa-2.0E0*3*I_NAI_G3xy_F2xz_a-2.0E0*4*I_NAI_G3xy_F2xz_a+3*2*I_NAI_Dxy_F2xz;
  abcd[32] = 4.0E0*I_NAI_I5xz_F2xz_aa-2.0E0*3*I_NAI_G3xz_F2xz_a-2.0E0*4*I_NAI_G3xz_F2xz_a+3*2*I_NAI_Dxz_F2xz;
  abcd[33] = 4.0E0*I_NAI_I4x2y_F2xz_aa-2.0E0*2*I_NAI_G2x2y_F2xz_a-2.0E0*3*I_NAI_G2x2y_F2xz_a+2*1*I_NAI_D2y_F2xz;
  abcd[34] = 4.0E0*I_NAI_I4xyz_F2xz_aa-2.0E0*2*I_NAI_G2xyz_F2xz_a-2.0E0*3*I_NAI_G2xyz_F2xz_a+2*1*I_NAI_Dyz_F2xz;
  abcd[35] = 4.0E0*I_NAI_I4x2z_F2xz_aa-2.0E0*2*I_NAI_G2x2z_F2xz_a-2.0E0*3*I_NAI_G2x2z_F2xz_a+2*1*I_NAI_D2z_F2xz;
  abcd[36] = 4.0E0*I_NAI_I3x3y_F2xz_aa-2.0E0*1*I_NAI_Gx3y_F2xz_a-2.0E0*2*I_NAI_Gx3y_F2xz_a;
  abcd[37] = 4.0E0*I_NAI_I3x2yz_F2xz_aa-2.0E0*1*I_NAI_Gx2yz_F2xz_a-2.0E0*2*I_NAI_Gx2yz_F2xz_a;
  abcd[38] = 4.0E0*I_NAI_I3xy2z_F2xz_aa-2.0E0*1*I_NAI_Gxy2z_F2xz_a-2.0E0*2*I_NAI_Gxy2z_F2xz_a;
  abcd[39] = 4.0E0*I_NAI_I3x3z_F2xz_aa-2.0E0*1*I_NAI_Gx3z_F2xz_a-2.0E0*2*I_NAI_Gx3z_F2xz_a;
  abcd[40] = 4.0E0*I_NAI_I2x4y_F2xz_aa-2.0E0*1*I_NAI_G4y_F2xz_a;
  abcd[41] = 4.0E0*I_NAI_I2x3yz_F2xz_aa-2.0E0*1*I_NAI_G3yz_F2xz_a;
  abcd[42] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2y2z_F2xz_a;
  abcd[43] = 4.0E0*I_NAI_I2xy3z_F2xz_aa-2.0E0*1*I_NAI_Gy3z_F2xz_a;
  abcd[44] = 4.0E0*I_NAI_I2x4z_F2xz_aa-2.0E0*1*I_NAI_G4z_F2xz_a;
  abcd[45] = 4.0E0*I_NAI_I6x_Fx2y_aa-2.0E0*4*I_NAI_G4x_Fx2y_a-2.0E0*5*I_NAI_G4x_Fx2y_a+4*3*I_NAI_D2x_Fx2y;
  abcd[46] = 4.0E0*I_NAI_I5xy_Fx2y_aa-2.0E0*3*I_NAI_G3xy_Fx2y_a-2.0E0*4*I_NAI_G3xy_Fx2y_a+3*2*I_NAI_Dxy_Fx2y;
  abcd[47] = 4.0E0*I_NAI_I5xz_Fx2y_aa-2.0E0*3*I_NAI_G3xz_Fx2y_a-2.0E0*4*I_NAI_G3xz_Fx2y_a+3*2*I_NAI_Dxz_Fx2y;
  abcd[48] = 4.0E0*I_NAI_I4x2y_Fx2y_aa-2.0E0*2*I_NAI_G2x2y_Fx2y_a-2.0E0*3*I_NAI_G2x2y_Fx2y_a+2*1*I_NAI_D2y_Fx2y;
  abcd[49] = 4.0E0*I_NAI_I4xyz_Fx2y_aa-2.0E0*2*I_NAI_G2xyz_Fx2y_a-2.0E0*3*I_NAI_G2xyz_Fx2y_a+2*1*I_NAI_Dyz_Fx2y;
  abcd[50] = 4.0E0*I_NAI_I4x2z_Fx2y_aa-2.0E0*2*I_NAI_G2x2z_Fx2y_a-2.0E0*3*I_NAI_G2x2z_Fx2y_a+2*1*I_NAI_D2z_Fx2y;
  abcd[51] = 4.0E0*I_NAI_I3x3y_Fx2y_aa-2.0E0*1*I_NAI_Gx3y_Fx2y_a-2.0E0*2*I_NAI_Gx3y_Fx2y_a;
  abcd[52] = 4.0E0*I_NAI_I3x2yz_Fx2y_aa-2.0E0*1*I_NAI_Gx2yz_Fx2y_a-2.0E0*2*I_NAI_Gx2yz_Fx2y_a;
  abcd[53] = 4.0E0*I_NAI_I3xy2z_Fx2y_aa-2.0E0*1*I_NAI_Gxy2z_Fx2y_a-2.0E0*2*I_NAI_Gxy2z_Fx2y_a;
  abcd[54] = 4.0E0*I_NAI_I3x3z_Fx2y_aa-2.0E0*1*I_NAI_Gx3z_Fx2y_a-2.0E0*2*I_NAI_Gx3z_Fx2y_a;
  abcd[55] = 4.0E0*I_NAI_I2x4y_Fx2y_aa-2.0E0*1*I_NAI_G4y_Fx2y_a;
  abcd[56] = 4.0E0*I_NAI_I2x3yz_Fx2y_aa-2.0E0*1*I_NAI_G3yz_Fx2y_a;
  abcd[57] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2y2z_Fx2y_a;
  abcd[58] = 4.0E0*I_NAI_I2xy3z_Fx2y_aa-2.0E0*1*I_NAI_Gy3z_Fx2y_a;
  abcd[59] = 4.0E0*I_NAI_I2x4z_Fx2y_aa-2.0E0*1*I_NAI_G4z_Fx2y_a;
  abcd[60] = 4.0E0*I_NAI_I6x_Fxyz_aa-2.0E0*4*I_NAI_G4x_Fxyz_a-2.0E0*5*I_NAI_G4x_Fxyz_a+4*3*I_NAI_D2x_Fxyz;
  abcd[61] = 4.0E0*I_NAI_I5xy_Fxyz_aa-2.0E0*3*I_NAI_G3xy_Fxyz_a-2.0E0*4*I_NAI_G3xy_Fxyz_a+3*2*I_NAI_Dxy_Fxyz;
  abcd[62] = 4.0E0*I_NAI_I5xz_Fxyz_aa-2.0E0*3*I_NAI_G3xz_Fxyz_a-2.0E0*4*I_NAI_G3xz_Fxyz_a+3*2*I_NAI_Dxz_Fxyz;
  abcd[63] = 4.0E0*I_NAI_I4x2y_Fxyz_aa-2.0E0*2*I_NAI_G2x2y_Fxyz_a-2.0E0*3*I_NAI_G2x2y_Fxyz_a+2*1*I_NAI_D2y_Fxyz;
  abcd[64] = 4.0E0*I_NAI_I4xyz_Fxyz_aa-2.0E0*2*I_NAI_G2xyz_Fxyz_a-2.0E0*3*I_NAI_G2xyz_Fxyz_a+2*1*I_NAI_Dyz_Fxyz;
  abcd[65] = 4.0E0*I_NAI_I4x2z_Fxyz_aa-2.0E0*2*I_NAI_G2x2z_Fxyz_a-2.0E0*3*I_NAI_G2x2z_Fxyz_a+2*1*I_NAI_D2z_Fxyz;
  abcd[66] = 4.0E0*I_NAI_I3x3y_Fxyz_aa-2.0E0*1*I_NAI_Gx3y_Fxyz_a-2.0E0*2*I_NAI_Gx3y_Fxyz_a;
  abcd[67] = 4.0E0*I_NAI_I3x2yz_Fxyz_aa-2.0E0*1*I_NAI_Gx2yz_Fxyz_a-2.0E0*2*I_NAI_Gx2yz_Fxyz_a;
  abcd[68] = 4.0E0*I_NAI_I3xy2z_Fxyz_aa-2.0E0*1*I_NAI_Gxy2z_Fxyz_a-2.0E0*2*I_NAI_Gxy2z_Fxyz_a;
  abcd[69] = 4.0E0*I_NAI_I3x3z_Fxyz_aa-2.0E0*1*I_NAI_Gx3z_Fxyz_a-2.0E0*2*I_NAI_Gx3z_Fxyz_a;
  abcd[70] = 4.0E0*I_NAI_I2x4y_Fxyz_aa-2.0E0*1*I_NAI_G4y_Fxyz_a;
  abcd[71] = 4.0E0*I_NAI_I2x3yz_Fxyz_aa-2.0E0*1*I_NAI_G3yz_Fxyz_a;
  abcd[72] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2y2z_Fxyz_a;
  abcd[73] = 4.0E0*I_NAI_I2xy3z_Fxyz_aa-2.0E0*1*I_NAI_Gy3z_Fxyz_a;
  abcd[74] = 4.0E0*I_NAI_I2x4z_Fxyz_aa-2.0E0*1*I_NAI_G4z_Fxyz_a;
  abcd[75] = 4.0E0*I_NAI_I6x_Fx2z_aa-2.0E0*4*I_NAI_G4x_Fx2z_a-2.0E0*5*I_NAI_G4x_Fx2z_a+4*3*I_NAI_D2x_Fx2z;
  abcd[76] = 4.0E0*I_NAI_I5xy_Fx2z_aa-2.0E0*3*I_NAI_G3xy_Fx2z_a-2.0E0*4*I_NAI_G3xy_Fx2z_a+3*2*I_NAI_Dxy_Fx2z;
  abcd[77] = 4.0E0*I_NAI_I5xz_Fx2z_aa-2.0E0*3*I_NAI_G3xz_Fx2z_a-2.0E0*4*I_NAI_G3xz_Fx2z_a+3*2*I_NAI_Dxz_Fx2z;
  abcd[78] = 4.0E0*I_NAI_I4x2y_Fx2z_aa-2.0E0*2*I_NAI_G2x2y_Fx2z_a-2.0E0*3*I_NAI_G2x2y_Fx2z_a+2*1*I_NAI_D2y_Fx2z;
  abcd[79] = 4.0E0*I_NAI_I4xyz_Fx2z_aa-2.0E0*2*I_NAI_G2xyz_Fx2z_a-2.0E0*3*I_NAI_G2xyz_Fx2z_a+2*1*I_NAI_Dyz_Fx2z;
  abcd[80] = 4.0E0*I_NAI_I4x2z_Fx2z_aa-2.0E0*2*I_NAI_G2x2z_Fx2z_a-2.0E0*3*I_NAI_G2x2z_Fx2z_a+2*1*I_NAI_D2z_Fx2z;
  abcd[81] = 4.0E0*I_NAI_I3x3y_Fx2z_aa-2.0E0*1*I_NAI_Gx3y_Fx2z_a-2.0E0*2*I_NAI_Gx3y_Fx2z_a;
  abcd[82] = 4.0E0*I_NAI_I3x2yz_Fx2z_aa-2.0E0*1*I_NAI_Gx2yz_Fx2z_a-2.0E0*2*I_NAI_Gx2yz_Fx2z_a;
  abcd[83] = 4.0E0*I_NAI_I3xy2z_Fx2z_aa-2.0E0*1*I_NAI_Gxy2z_Fx2z_a-2.0E0*2*I_NAI_Gxy2z_Fx2z_a;
  abcd[84] = 4.0E0*I_NAI_I3x3z_Fx2z_aa-2.0E0*1*I_NAI_Gx3z_Fx2z_a-2.0E0*2*I_NAI_Gx3z_Fx2z_a;
  abcd[85] = 4.0E0*I_NAI_I2x4y_Fx2z_aa-2.0E0*1*I_NAI_G4y_Fx2z_a;
  abcd[86] = 4.0E0*I_NAI_I2x3yz_Fx2z_aa-2.0E0*1*I_NAI_G3yz_Fx2z_a;
  abcd[87] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2y2z_Fx2z_a;
  abcd[88] = 4.0E0*I_NAI_I2xy3z_Fx2z_aa-2.0E0*1*I_NAI_Gy3z_Fx2z_a;
  abcd[89] = 4.0E0*I_NAI_I2x4z_Fx2z_aa-2.0E0*1*I_NAI_G4z_Fx2z_a;
  abcd[90] = 4.0E0*I_NAI_I6x_F3y_aa-2.0E0*4*I_NAI_G4x_F3y_a-2.0E0*5*I_NAI_G4x_F3y_a+4*3*I_NAI_D2x_F3y;
  abcd[91] = 4.0E0*I_NAI_I5xy_F3y_aa-2.0E0*3*I_NAI_G3xy_F3y_a-2.0E0*4*I_NAI_G3xy_F3y_a+3*2*I_NAI_Dxy_F3y;
  abcd[92] = 4.0E0*I_NAI_I5xz_F3y_aa-2.0E0*3*I_NAI_G3xz_F3y_a-2.0E0*4*I_NAI_G3xz_F3y_a+3*2*I_NAI_Dxz_F3y;
  abcd[93] = 4.0E0*I_NAI_I4x2y_F3y_aa-2.0E0*2*I_NAI_G2x2y_F3y_a-2.0E0*3*I_NAI_G2x2y_F3y_a+2*1*I_NAI_D2y_F3y;
  abcd[94] = 4.0E0*I_NAI_I4xyz_F3y_aa-2.0E0*2*I_NAI_G2xyz_F3y_a-2.0E0*3*I_NAI_G2xyz_F3y_a+2*1*I_NAI_Dyz_F3y;
  abcd[95] = 4.0E0*I_NAI_I4x2z_F3y_aa-2.0E0*2*I_NAI_G2x2z_F3y_a-2.0E0*3*I_NAI_G2x2z_F3y_a+2*1*I_NAI_D2z_F3y;
  abcd[96] = 4.0E0*I_NAI_I3x3y_F3y_aa-2.0E0*1*I_NAI_Gx3y_F3y_a-2.0E0*2*I_NAI_Gx3y_F3y_a;
  abcd[97] = 4.0E0*I_NAI_I3x2yz_F3y_aa-2.0E0*1*I_NAI_Gx2yz_F3y_a-2.0E0*2*I_NAI_Gx2yz_F3y_a;
  abcd[98] = 4.0E0*I_NAI_I3xy2z_F3y_aa-2.0E0*1*I_NAI_Gxy2z_F3y_a-2.0E0*2*I_NAI_Gxy2z_F3y_a;
  abcd[99] = 4.0E0*I_NAI_I3x3z_F3y_aa-2.0E0*1*I_NAI_Gx3z_F3y_a-2.0E0*2*I_NAI_Gx3z_F3y_a;
  abcd[100] = 4.0E0*I_NAI_I2x4y_F3y_aa-2.0E0*1*I_NAI_G4y_F3y_a;
  abcd[101] = 4.0E0*I_NAI_I2x3yz_F3y_aa-2.0E0*1*I_NAI_G3yz_F3y_a;
  abcd[102] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2y2z_F3y_a;
  abcd[103] = 4.0E0*I_NAI_I2xy3z_F3y_aa-2.0E0*1*I_NAI_Gy3z_F3y_a;
  abcd[104] = 4.0E0*I_NAI_I2x4z_F3y_aa-2.0E0*1*I_NAI_G4z_F3y_a;
  abcd[105] = 4.0E0*I_NAI_I6x_F2yz_aa-2.0E0*4*I_NAI_G4x_F2yz_a-2.0E0*5*I_NAI_G4x_F2yz_a+4*3*I_NAI_D2x_F2yz;
  abcd[106] = 4.0E0*I_NAI_I5xy_F2yz_aa-2.0E0*3*I_NAI_G3xy_F2yz_a-2.0E0*4*I_NAI_G3xy_F2yz_a+3*2*I_NAI_Dxy_F2yz;
  abcd[107] = 4.0E0*I_NAI_I5xz_F2yz_aa-2.0E0*3*I_NAI_G3xz_F2yz_a-2.0E0*4*I_NAI_G3xz_F2yz_a+3*2*I_NAI_Dxz_F2yz;
  abcd[108] = 4.0E0*I_NAI_I4x2y_F2yz_aa-2.0E0*2*I_NAI_G2x2y_F2yz_a-2.0E0*3*I_NAI_G2x2y_F2yz_a+2*1*I_NAI_D2y_F2yz;
  abcd[109] = 4.0E0*I_NAI_I4xyz_F2yz_aa-2.0E0*2*I_NAI_G2xyz_F2yz_a-2.0E0*3*I_NAI_G2xyz_F2yz_a+2*1*I_NAI_Dyz_F2yz;
  abcd[110] = 4.0E0*I_NAI_I4x2z_F2yz_aa-2.0E0*2*I_NAI_G2x2z_F2yz_a-2.0E0*3*I_NAI_G2x2z_F2yz_a+2*1*I_NAI_D2z_F2yz;
  abcd[111] = 4.0E0*I_NAI_I3x3y_F2yz_aa-2.0E0*1*I_NAI_Gx3y_F2yz_a-2.0E0*2*I_NAI_Gx3y_F2yz_a;
  abcd[112] = 4.0E0*I_NAI_I3x2yz_F2yz_aa-2.0E0*1*I_NAI_Gx2yz_F2yz_a-2.0E0*2*I_NAI_Gx2yz_F2yz_a;
  abcd[113] = 4.0E0*I_NAI_I3xy2z_F2yz_aa-2.0E0*1*I_NAI_Gxy2z_F2yz_a-2.0E0*2*I_NAI_Gxy2z_F2yz_a;
  abcd[114] = 4.0E0*I_NAI_I3x3z_F2yz_aa-2.0E0*1*I_NAI_Gx3z_F2yz_a-2.0E0*2*I_NAI_Gx3z_F2yz_a;
  abcd[115] = 4.0E0*I_NAI_I2x4y_F2yz_aa-2.0E0*1*I_NAI_G4y_F2yz_a;
  abcd[116] = 4.0E0*I_NAI_I2x3yz_F2yz_aa-2.0E0*1*I_NAI_G3yz_F2yz_a;
  abcd[117] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2y2z_F2yz_a;
  abcd[118] = 4.0E0*I_NAI_I2xy3z_F2yz_aa-2.0E0*1*I_NAI_Gy3z_F2yz_a;
  abcd[119] = 4.0E0*I_NAI_I2x4z_F2yz_aa-2.0E0*1*I_NAI_G4z_F2yz_a;
  abcd[120] = 4.0E0*I_NAI_I6x_Fy2z_aa-2.0E0*4*I_NAI_G4x_Fy2z_a-2.0E0*5*I_NAI_G4x_Fy2z_a+4*3*I_NAI_D2x_Fy2z;
  abcd[121] = 4.0E0*I_NAI_I5xy_Fy2z_aa-2.0E0*3*I_NAI_G3xy_Fy2z_a-2.0E0*4*I_NAI_G3xy_Fy2z_a+3*2*I_NAI_Dxy_Fy2z;
  abcd[122] = 4.0E0*I_NAI_I5xz_Fy2z_aa-2.0E0*3*I_NAI_G3xz_Fy2z_a-2.0E0*4*I_NAI_G3xz_Fy2z_a+3*2*I_NAI_Dxz_Fy2z;
  abcd[123] = 4.0E0*I_NAI_I4x2y_Fy2z_aa-2.0E0*2*I_NAI_G2x2y_Fy2z_a-2.0E0*3*I_NAI_G2x2y_Fy2z_a+2*1*I_NAI_D2y_Fy2z;
  abcd[124] = 4.0E0*I_NAI_I4xyz_Fy2z_aa-2.0E0*2*I_NAI_G2xyz_Fy2z_a-2.0E0*3*I_NAI_G2xyz_Fy2z_a+2*1*I_NAI_Dyz_Fy2z;
  abcd[125] = 4.0E0*I_NAI_I4x2z_Fy2z_aa-2.0E0*2*I_NAI_G2x2z_Fy2z_a-2.0E0*3*I_NAI_G2x2z_Fy2z_a+2*1*I_NAI_D2z_Fy2z;
  abcd[126] = 4.0E0*I_NAI_I3x3y_Fy2z_aa-2.0E0*1*I_NAI_Gx3y_Fy2z_a-2.0E0*2*I_NAI_Gx3y_Fy2z_a;
  abcd[127] = 4.0E0*I_NAI_I3x2yz_Fy2z_aa-2.0E0*1*I_NAI_Gx2yz_Fy2z_a-2.0E0*2*I_NAI_Gx2yz_Fy2z_a;
  abcd[128] = 4.0E0*I_NAI_I3xy2z_Fy2z_aa-2.0E0*1*I_NAI_Gxy2z_Fy2z_a-2.0E0*2*I_NAI_Gxy2z_Fy2z_a;
  abcd[129] = 4.0E0*I_NAI_I3x3z_Fy2z_aa-2.0E0*1*I_NAI_Gx3z_Fy2z_a-2.0E0*2*I_NAI_Gx3z_Fy2z_a;
  abcd[130] = 4.0E0*I_NAI_I2x4y_Fy2z_aa-2.0E0*1*I_NAI_G4y_Fy2z_a;
  abcd[131] = 4.0E0*I_NAI_I2x3yz_Fy2z_aa-2.0E0*1*I_NAI_G3yz_Fy2z_a;
  abcd[132] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2y2z_Fy2z_a;
  abcd[133] = 4.0E0*I_NAI_I2xy3z_Fy2z_aa-2.0E0*1*I_NAI_Gy3z_Fy2z_a;
  abcd[134] = 4.0E0*I_NAI_I2x4z_Fy2z_aa-2.0E0*1*I_NAI_G4z_Fy2z_a;
  abcd[135] = 4.0E0*I_NAI_I6x_F3z_aa-2.0E0*4*I_NAI_G4x_F3z_a-2.0E0*5*I_NAI_G4x_F3z_a+4*3*I_NAI_D2x_F3z;
  abcd[136] = 4.0E0*I_NAI_I5xy_F3z_aa-2.0E0*3*I_NAI_G3xy_F3z_a-2.0E0*4*I_NAI_G3xy_F3z_a+3*2*I_NAI_Dxy_F3z;
  abcd[137] = 4.0E0*I_NAI_I5xz_F3z_aa-2.0E0*3*I_NAI_G3xz_F3z_a-2.0E0*4*I_NAI_G3xz_F3z_a+3*2*I_NAI_Dxz_F3z;
  abcd[138] = 4.0E0*I_NAI_I4x2y_F3z_aa-2.0E0*2*I_NAI_G2x2y_F3z_a-2.0E0*3*I_NAI_G2x2y_F3z_a+2*1*I_NAI_D2y_F3z;
  abcd[139] = 4.0E0*I_NAI_I4xyz_F3z_aa-2.0E0*2*I_NAI_G2xyz_F3z_a-2.0E0*3*I_NAI_G2xyz_F3z_a+2*1*I_NAI_Dyz_F3z;
  abcd[140] = 4.0E0*I_NAI_I4x2z_F3z_aa-2.0E0*2*I_NAI_G2x2z_F3z_a-2.0E0*3*I_NAI_G2x2z_F3z_a+2*1*I_NAI_D2z_F3z;
  abcd[141] = 4.0E0*I_NAI_I3x3y_F3z_aa-2.0E0*1*I_NAI_Gx3y_F3z_a-2.0E0*2*I_NAI_Gx3y_F3z_a;
  abcd[142] = 4.0E0*I_NAI_I3x2yz_F3z_aa-2.0E0*1*I_NAI_Gx2yz_F3z_a-2.0E0*2*I_NAI_Gx2yz_F3z_a;
  abcd[143] = 4.0E0*I_NAI_I3xy2z_F3z_aa-2.0E0*1*I_NAI_Gxy2z_F3z_a-2.0E0*2*I_NAI_Gxy2z_F3z_a;
  abcd[144] = 4.0E0*I_NAI_I3x3z_F3z_aa-2.0E0*1*I_NAI_Gx3z_F3z_a-2.0E0*2*I_NAI_Gx3z_F3z_a;
  abcd[145] = 4.0E0*I_NAI_I2x4y_F3z_aa-2.0E0*1*I_NAI_G4y_F3z_a;
  abcd[146] = 4.0E0*I_NAI_I2x3yz_F3z_aa-2.0E0*1*I_NAI_G3yz_F3z_a;
  abcd[147] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2y2z_F3z_a;
  abcd[148] = 4.0E0*I_NAI_I2xy3z_F3z_aa-2.0E0*1*I_NAI_Gy3z_F3z_a;
  abcd[149] = 4.0E0*I_NAI_I2x4z_F3z_aa-2.0E0*1*I_NAI_G4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[150] = 4.0E0*I_NAI_I5xy_F3x_aa-2.0E0*4*I_NAI_G3xy_F3x_a;
  abcd[151] = 4.0E0*I_NAI_I4x2y_F3x_aa-2.0E0*1*I_NAI_G4x_F3x_a-2.0E0*3*I_NAI_G2x2y_F3x_a+3*1*I_NAI_D2x_F3x;
  abcd[152] = 4.0E0*I_NAI_I4xyz_F3x_aa-2.0E0*3*I_NAI_G2xyz_F3x_a;
  abcd[153] = 4.0E0*I_NAI_I3x3y_F3x_aa-2.0E0*2*I_NAI_G3xy_F3x_a-2.0E0*2*I_NAI_Gx3y_F3x_a+2*2*I_NAI_Dxy_F3x;
  abcd[154] = 4.0E0*I_NAI_I3x2yz_F3x_aa-2.0E0*1*I_NAI_G3xz_F3x_a-2.0E0*2*I_NAI_Gx2yz_F3x_a+2*1*I_NAI_Dxz_F3x;
  abcd[155] = 4.0E0*I_NAI_I3xy2z_F3x_aa-2.0E0*2*I_NAI_Gxy2z_F3x_a;
  abcd[156] = 4.0E0*I_NAI_I2x4y_F3x_aa-2.0E0*3*I_NAI_G2x2y_F3x_a-2.0E0*1*I_NAI_G4y_F3x_a+3*I_NAI_D2y_F3x;
  abcd[157] = 4.0E0*I_NAI_I2x3yz_F3x_aa-2.0E0*2*I_NAI_G2xyz_F3x_a-2.0E0*1*I_NAI_G3yz_F3x_a+2*I_NAI_Dyz_F3x;
  abcd[158] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2x2z_F3x_a-2.0E0*1*I_NAI_G2y2z_F3x_a+1*I_NAI_D2z_F3x;
  abcd[159] = 4.0E0*I_NAI_I2xy3z_F3x_aa-2.0E0*1*I_NAI_Gy3z_F3x_a;
  abcd[160] = 4.0E0*I_NAI_Ix5y_F3x_aa-2.0E0*4*I_NAI_Gx3y_F3x_a;
  abcd[161] = 4.0E0*I_NAI_Ix4yz_F3x_aa-2.0E0*3*I_NAI_Gx2yz_F3x_a;
  abcd[162] = 4.0E0*I_NAI_Ix3y2z_F3x_aa-2.0E0*2*I_NAI_Gxy2z_F3x_a;
  abcd[163] = 4.0E0*I_NAI_Ix2y3z_F3x_aa-2.0E0*1*I_NAI_Gx3z_F3x_a;
  abcd[164] = 4.0E0*I_NAI_Ixy4z_F3x_aa;
  abcd[165] = 4.0E0*I_NAI_I5xy_F2xy_aa-2.0E0*4*I_NAI_G3xy_F2xy_a;
  abcd[166] = 4.0E0*I_NAI_I4x2y_F2xy_aa-2.0E0*1*I_NAI_G4x_F2xy_a-2.0E0*3*I_NAI_G2x2y_F2xy_a+3*1*I_NAI_D2x_F2xy;
  abcd[167] = 4.0E0*I_NAI_I4xyz_F2xy_aa-2.0E0*3*I_NAI_G2xyz_F2xy_a;
  abcd[168] = 4.0E0*I_NAI_I3x3y_F2xy_aa-2.0E0*2*I_NAI_G3xy_F2xy_a-2.0E0*2*I_NAI_Gx3y_F2xy_a+2*2*I_NAI_Dxy_F2xy;
  abcd[169] = 4.0E0*I_NAI_I3x2yz_F2xy_aa-2.0E0*1*I_NAI_G3xz_F2xy_a-2.0E0*2*I_NAI_Gx2yz_F2xy_a+2*1*I_NAI_Dxz_F2xy;
  abcd[170] = 4.0E0*I_NAI_I3xy2z_F2xy_aa-2.0E0*2*I_NAI_Gxy2z_F2xy_a;
  abcd[171] = 4.0E0*I_NAI_I2x4y_F2xy_aa-2.0E0*3*I_NAI_G2x2y_F2xy_a-2.0E0*1*I_NAI_G4y_F2xy_a+3*I_NAI_D2y_F2xy;
  abcd[172] = 4.0E0*I_NAI_I2x3yz_F2xy_aa-2.0E0*2*I_NAI_G2xyz_F2xy_a-2.0E0*1*I_NAI_G3yz_F2xy_a+2*I_NAI_Dyz_F2xy;
  abcd[173] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2x2z_F2xy_a-2.0E0*1*I_NAI_G2y2z_F2xy_a+1*I_NAI_D2z_F2xy;
  abcd[174] = 4.0E0*I_NAI_I2xy3z_F2xy_aa-2.0E0*1*I_NAI_Gy3z_F2xy_a;
  abcd[175] = 4.0E0*I_NAI_Ix5y_F2xy_aa-2.0E0*4*I_NAI_Gx3y_F2xy_a;
  abcd[176] = 4.0E0*I_NAI_Ix4yz_F2xy_aa-2.0E0*3*I_NAI_Gx2yz_F2xy_a;
  abcd[177] = 4.0E0*I_NAI_Ix3y2z_F2xy_aa-2.0E0*2*I_NAI_Gxy2z_F2xy_a;
  abcd[178] = 4.0E0*I_NAI_Ix2y3z_F2xy_aa-2.0E0*1*I_NAI_Gx3z_F2xy_a;
  abcd[179] = 4.0E0*I_NAI_Ixy4z_F2xy_aa;
  abcd[180] = 4.0E0*I_NAI_I5xy_F2xz_aa-2.0E0*4*I_NAI_G3xy_F2xz_a;
  abcd[181] = 4.0E0*I_NAI_I4x2y_F2xz_aa-2.0E0*1*I_NAI_G4x_F2xz_a-2.0E0*3*I_NAI_G2x2y_F2xz_a+3*1*I_NAI_D2x_F2xz;
  abcd[182] = 4.0E0*I_NAI_I4xyz_F2xz_aa-2.0E0*3*I_NAI_G2xyz_F2xz_a;
  abcd[183] = 4.0E0*I_NAI_I3x3y_F2xz_aa-2.0E0*2*I_NAI_G3xy_F2xz_a-2.0E0*2*I_NAI_Gx3y_F2xz_a+2*2*I_NAI_Dxy_F2xz;
  abcd[184] = 4.0E0*I_NAI_I3x2yz_F2xz_aa-2.0E0*1*I_NAI_G3xz_F2xz_a-2.0E0*2*I_NAI_Gx2yz_F2xz_a+2*1*I_NAI_Dxz_F2xz;
  abcd[185] = 4.0E0*I_NAI_I3xy2z_F2xz_aa-2.0E0*2*I_NAI_Gxy2z_F2xz_a;
  abcd[186] = 4.0E0*I_NAI_I2x4y_F2xz_aa-2.0E0*3*I_NAI_G2x2y_F2xz_a-2.0E0*1*I_NAI_G4y_F2xz_a+3*I_NAI_D2y_F2xz;
  abcd[187] = 4.0E0*I_NAI_I2x3yz_F2xz_aa-2.0E0*2*I_NAI_G2xyz_F2xz_a-2.0E0*1*I_NAI_G3yz_F2xz_a+2*I_NAI_Dyz_F2xz;
  abcd[188] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2x2z_F2xz_a-2.0E0*1*I_NAI_G2y2z_F2xz_a+1*I_NAI_D2z_F2xz;
  abcd[189] = 4.0E0*I_NAI_I2xy3z_F2xz_aa-2.0E0*1*I_NAI_Gy3z_F2xz_a;
  abcd[190] = 4.0E0*I_NAI_Ix5y_F2xz_aa-2.0E0*4*I_NAI_Gx3y_F2xz_a;
  abcd[191] = 4.0E0*I_NAI_Ix4yz_F2xz_aa-2.0E0*3*I_NAI_Gx2yz_F2xz_a;
  abcd[192] = 4.0E0*I_NAI_Ix3y2z_F2xz_aa-2.0E0*2*I_NAI_Gxy2z_F2xz_a;
  abcd[193] = 4.0E0*I_NAI_Ix2y3z_F2xz_aa-2.0E0*1*I_NAI_Gx3z_F2xz_a;
  abcd[194] = 4.0E0*I_NAI_Ixy4z_F2xz_aa;
  abcd[195] = 4.0E0*I_NAI_I5xy_Fx2y_aa-2.0E0*4*I_NAI_G3xy_Fx2y_a;
  abcd[196] = 4.0E0*I_NAI_I4x2y_Fx2y_aa-2.0E0*1*I_NAI_G4x_Fx2y_a-2.0E0*3*I_NAI_G2x2y_Fx2y_a+3*1*I_NAI_D2x_Fx2y;
  abcd[197] = 4.0E0*I_NAI_I4xyz_Fx2y_aa-2.0E0*3*I_NAI_G2xyz_Fx2y_a;
  abcd[198] = 4.0E0*I_NAI_I3x3y_Fx2y_aa-2.0E0*2*I_NAI_G3xy_Fx2y_a-2.0E0*2*I_NAI_Gx3y_Fx2y_a+2*2*I_NAI_Dxy_Fx2y;
  abcd[199] = 4.0E0*I_NAI_I3x2yz_Fx2y_aa-2.0E0*1*I_NAI_G3xz_Fx2y_a-2.0E0*2*I_NAI_Gx2yz_Fx2y_a+2*1*I_NAI_Dxz_Fx2y;
  abcd[200] = 4.0E0*I_NAI_I3xy2z_Fx2y_aa-2.0E0*2*I_NAI_Gxy2z_Fx2y_a;
  abcd[201] = 4.0E0*I_NAI_I2x4y_Fx2y_aa-2.0E0*3*I_NAI_G2x2y_Fx2y_a-2.0E0*1*I_NAI_G4y_Fx2y_a+3*I_NAI_D2y_Fx2y;
  abcd[202] = 4.0E0*I_NAI_I2x3yz_Fx2y_aa-2.0E0*2*I_NAI_G2xyz_Fx2y_a-2.0E0*1*I_NAI_G3yz_Fx2y_a+2*I_NAI_Dyz_Fx2y;
  abcd[203] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2x2z_Fx2y_a-2.0E0*1*I_NAI_G2y2z_Fx2y_a+1*I_NAI_D2z_Fx2y;
  abcd[204] = 4.0E0*I_NAI_I2xy3z_Fx2y_aa-2.0E0*1*I_NAI_Gy3z_Fx2y_a;
  abcd[205] = 4.0E0*I_NAI_Ix5y_Fx2y_aa-2.0E0*4*I_NAI_Gx3y_Fx2y_a;
  abcd[206] = 4.0E0*I_NAI_Ix4yz_Fx2y_aa-2.0E0*3*I_NAI_Gx2yz_Fx2y_a;
  abcd[207] = 4.0E0*I_NAI_Ix3y2z_Fx2y_aa-2.0E0*2*I_NAI_Gxy2z_Fx2y_a;
  abcd[208] = 4.0E0*I_NAI_Ix2y3z_Fx2y_aa-2.0E0*1*I_NAI_Gx3z_Fx2y_a;
  abcd[209] = 4.0E0*I_NAI_Ixy4z_Fx2y_aa;
  abcd[210] = 4.0E0*I_NAI_I5xy_Fxyz_aa-2.0E0*4*I_NAI_G3xy_Fxyz_a;
  abcd[211] = 4.0E0*I_NAI_I4x2y_Fxyz_aa-2.0E0*1*I_NAI_G4x_Fxyz_a-2.0E0*3*I_NAI_G2x2y_Fxyz_a+3*1*I_NAI_D2x_Fxyz;
  abcd[212] = 4.0E0*I_NAI_I4xyz_Fxyz_aa-2.0E0*3*I_NAI_G2xyz_Fxyz_a;
  abcd[213] = 4.0E0*I_NAI_I3x3y_Fxyz_aa-2.0E0*2*I_NAI_G3xy_Fxyz_a-2.0E0*2*I_NAI_Gx3y_Fxyz_a+2*2*I_NAI_Dxy_Fxyz;
  abcd[214] = 4.0E0*I_NAI_I3x2yz_Fxyz_aa-2.0E0*1*I_NAI_G3xz_Fxyz_a-2.0E0*2*I_NAI_Gx2yz_Fxyz_a+2*1*I_NAI_Dxz_Fxyz;
  abcd[215] = 4.0E0*I_NAI_I3xy2z_Fxyz_aa-2.0E0*2*I_NAI_Gxy2z_Fxyz_a;
  abcd[216] = 4.0E0*I_NAI_I2x4y_Fxyz_aa-2.0E0*3*I_NAI_G2x2y_Fxyz_a-2.0E0*1*I_NAI_G4y_Fxyz_a+3*I_NAI_D2y_Fxyz;
  abcd[217] = 4.0E0*I_NAI_I2x3yz_Fxyz_aa-2.0E0*2*I_NAI_G2xyz_Fxyz_a-2.0E0*1*I_NAI_G3yz_Fxyz_a+2*I_NAI_Dyz_Fxyz;
  abcd[218] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2x2z_Fxyz_a-2.0E0*1*I_NAI_G2y2z_Fxyz_a+1*I_NAI_D2z_Fxyz;
  abcd[219] = 4.0E0*I_NAI_I2xy3z_Fxyz_aa-2.0E0*1*I_NAI_Gy3z_Fxyz_a;
  abcd[220] = 4.0E0*I_NAI_Ix5y_Fxyz_aa-2.0E0*4*I_NAI_Gx3y_Fxyz_a;
  abcd[221] = 4.0E0*I_NAI_Ix4yz_Fxyz_aa-2.0E0*3*I_NAI_Gx2yz_Fxyz_a;
  abcd[222] = 4.0E0*I_NAI_Ix3y2z_Fxyz_aa-2.0E0*2*I_NAI_Gxy2z_Fxyz_a;
  abcd[223] = 4.0E0*I_NAI_Ix2y3z_Fxyz_aa-2.0E0*1*I_NAI_Gx3z_Fxyz_a;
  abcd[224] = 4.0E0*I_NAI_Ixy4z_Fxyz_aa;
  abcd[225] = 4.0E0*I_NAI_I5xy_Fx2z_aa-2.0E0*4*I_NAI_G3xy_Fx2z_a;
  abcd[226] = 4.0E0*I_NAI_I4x2y_Fx2z_aa-2.0E0*1*I_NAI_G4x_Fx2z_a-2.0E0*3*I_NAI_G2x2y_Fx2z_a+3*1*I_NAI_D2x_Fx2z;
  abcd[227] = 4.0E0*I_NAI_I4xyz_Fx2z_aa-2.0E0*3*I_NAI_G2xyz_Fx2z_a;
  abcd[228] = 4.0E0*I_NAI_I3x3y_Fx2z_aa-2.0E0*2*I_NAI_G3xy_Fx2z_a-2.0E0*2*I_NAI_Gx3y_Fx2z_a+2*2*I_NAI_Dxy_Fx2z;
  abcd[229] = 4.0E0*I_NAI_I3x2yz_Fx2z_aa-2.0E0*1*I_NAI_G3xz_Fx2z_a-2.0E0*2*I_NAI_Gx2yz_Fx2z_a+2*1*I_NAI_Dxz_Fx2z;
  abcd[230] = 4.0E0*I_NAI_I3xy2z_Fx2z_aa-2.0E0*2*I_NAI_Gxy2z_Fx2z_a;
  abcd[231] = 4.0E0*I_NAI_I2x4y_Fx2z_aa-2.0E0*3*I_NAI_G2x2y_Fx2z_a-2.0E0*1*I_NAI_G4y_Fx2z_a+3*I_NAI_D2y_Fx2z;
  abcd[232] = 4.0E0*I_NAI_I2x3yz_Fx2z_aa-2.0E0*2*I_NAI_G2xyz_Fx2z_a-2.0E0*1*I_NAI_G3yz_Fx2z_a+2*I_NAI_Dyz_Fx2z;
  abcd[233] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2x2z_Fx2z_a-2.0E0*1*I_NAI_G2y2z_Fx2z_a+1*I_NAI_D2z_Fx2z;
  abcd[234] = 4.0E0*I_NAI_I2xy3z_Fx2z_aa-2.0E0*1*I_NAI_Gy3z_Fx2z_a;
  abcd[235] = 4.0E0*I_NAI_Ix5y_Fx2z_aa-2.0E0*4*I_NAI_Gx3y_Fx2z_a;
  abcd[236] = 4.0E0*I_NAI_Ix4yz_Fx2z_aa-2.0E0*3*I_NAI_Gx2yz_Fx2z_a;
  abcd[237] = 4.0E0*I_NAI_Ix3y2z_Fx2z_aa-2.0E0*2*I_NAI_Gxy2z_Fx2z_a;
  abcd[238] = 4.0E0*I_NAI_Ix2y3z_Fx2z_aa-2.0E0*1*I_NAI_Gx3z_Fx2z_a;
  abcd[239] = 4.0E0*I_NAI_Ixy4z_Fx2z_aa;
  abcd[240] = 4.0E0*I_NAI_I5xy_F3y_aa-2.0E0*4*I_NAI_G3xy_F3y_a;
  abcd[241] = 4.0E0*I_NAI_I4x2y_F3y_aa-2.0E0*1*I_NAI_G4x_F3y_a-2.0E0*3*I_NAI_G2x2y_F3y_a+3*1*I_NAI_D2x_F3y;
  abcd[242] = 4.0E0*I_NAI_I4xyz_F3y_aa-2.0E0*3*I_NAI_G2xyz_F3y_a;
  abcd[243] = 4.0E0*I_NAI_I3x3y_F3y_aa-2.0E0*2*I_NAI_G3xy_F3y_a-2.0E0*2*I_NAI_Gx3y_F3y_a+2*2*I_NAI_Dxy_F3y;
  abcd[244] = 4.0E0*I_NAI_I3x2yz_F3y_aa-2.0E0*1*I_NAI_G3xz_F3y_a-2.0E0*2*I_NAI_Gx2yz_F3y_a+2*1*I_NAI_Dxz_F3y;
  abcd[245] = 4.0E0*I_NAI_I3xy2z_F3y_aa-2.0E0*2*I_NAI_Gxy2z_F3y_a;
  abcd[246] = 4.0E0*I_NAI_I2x4y_F3y_aa-2.0E0*3*I_NAI_G2x2y_F3y_a-2.0E0*1*I_NAI_G4y_F3y_a+3*I_NAI_D2y_F3y;
  abcd[247] = 4.0E0*I_NAI_I2x3yz_F3y_aa-2.0E0*2*I_NAI_G2xyz_F3y_a-2.0E0*1*I_NAI_G3yz_F3y_a+2*I_NAI_Dyz_F3y;
  abcd[248] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2x2z_F3y_a-2.0E0*1*I_NAI_G2y2z_F3y_a+1*I_NAI_D2z_F3y;
  abcd[249] = 4.0E0*I_NAI_I2xy3z_F3y_aa-2.0E0*1*I_NAI_Gy3z_F3y_a;
  abcd[250] = 4.0E0*I_NAI_Ix5y_F3y_aa-2.0E0*4*I_NAI_Gx3y_F3y_a;
  abcd[251] = 4.0E0*I_NAI_Ix4yz_F3y_aa-2.0E0*3*I_NAI_Gx2yz_F3y_a;
  abcd[252] = 4.0E0*I_NAI_Ix3y2z_F3y_aa-2.0E0*2*I_NAI_Gxy2z_F3y_a;
  abcd[253] = 4.0E0*I_NAI_Ix2y3z_F3y_aa-2.0E0*1*I_NAI_Gx3z_F3y_a;
  abcd[254] = 4.0E0*I_NAI_Ixy4z_F3y_aa;
  abcd[255] = 4.0E0*I_NAI_I5xy_F2yz_aa-2.0E0*4*I_NAI_G3xy_F2yz_a;
  abcd[256] = 4.0E0*I_NAI_I4x2y_F2yz_aa-2.0E0*1*I_NAI_G4x_F2yz_a-2.0E0*3*I_NAI_G2x2y_F2yz_a+3*1*I_NAI_D2x_F2yz;
  abcd[257] = 4.0E0*I_NAI_I4xyz_F2yz_aa-2.0E0*3*I_NAI_G2xyz_F2yz_a;
  abcd[258] = 4.0E0*I_NAI_I3x3y_F2yz_aa-2.0E0*2*I_NAI_G3xy_F2yz_a-2.0E0*2*I_NAI_Gx3y_F2yz_a+2*2*I_NAI_Dxy_F2yz;
  abcd[259] = 4.0E0*I_NAI_I3x2yz_F2yz_aa-2.0E0*1*I_NAI_G3xz_F2yz_a-2.0E0*2*I_NAI_Gx2yz_F2yz_a+2*1*I_NAI_Dxz_F2yz;
  abcd[260] = 4.0E0*I_NAI_I3xy2z_F2yz_aa-2.0E0*2*I_NAI_Gxy2z_F2yz_a;
  abcd[261] = 4.0E0*I_NAI_I2x4y_F2yz_aa-2.0E0*3*I_NAI_G2x2y_F2yz_a-2.0E0*1*I_NAI_G4y_F2yz_a+3*I_NAI_D2y_F2yz;
  abcd[262] = 4.0E0*I_NAI_I2x3yz_F2yz_aa-2.0E0*2*I_NAI_G2xyz_F2yz_a-2.0E0*1*I_NAI_G3yz_F2yz_a+2*I_NAI_Dyz_F2yz;
  abcd[263] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2x2z_F2yz_a-2.0E0*1*I_NAI_G2y2z_F2yz_a+1*I_NAI_D2z_F2yz;
  abcd[264] = 4.0E0*I_NAI_I2xy3z_F2yz_aa-2.0E0*1*I_NAI_Gy3z_F2yz_a;
  abcd[265] = 4.0E0*I_NAI_Ix5y_F2yz_aa-2.0E0*4*I_NAI_Gx3y_F2yz_a;
  abcd[266] = 4.0E0*I_NAI_Ix4yz_F2yz_aa-2.0E0*3*I_NAI_Gx2yz_F2yz_a;
  abcd[267] = 4.0E0*I_NAI_Ix3y2z_F2yz_aa-2.0E0*2*I_NAI_Gxy2z_F2yz_a;
  abcd[268] = 4.0E0*I_NAI_Ix2y3z_F2yz_aa-2.0E0*1*I_NAI_Gx3z_F2yz_a;
  abcd[269] = 4.0E0*I_NAI_Ixy4z_F2yz_aa;
  abcd[270] = 4.0E0*I_NAI_I5xy_Fy2z_aa-2.0E0*4*I_NAI_G3xy_Fy2z_a;
  abcd[271] = 4.0E0*I_NAI_I4x2y_Fy2z_aa-2.0E0*1*I_NAI_G4x_Fy2z_a-2.0E0*3*I_NAI_G2x2y_Fy2z_a+3*1*I_NAI_D2x_Fy2z;
  abcd[272] = 4.0E0*I_NAI_I4xyz_Fy2z_aa-2.0E0*3*I_NAI_G2xyz_Fy2z_a;
  abcd[273] = 4.0E0*I_NAI_I3x3y_Fy2z_aa-2.0E0*2*I_NAI_G3xy_Fy2z_a-2.0E0*2*I_NAI_Gx3y_Fy2z_a+2*2*I_NAI_Dxy_Fy2z;
  abcd[274] = 4.0E0*I_NAI_I3x2yz_Fy2z_aa-2.0E0*1*I_NAI_G3xz_Fy2z_a-2.0E0*2*I_NAI_Gx2yz_Fy2z_a+2*1*I_NAI_Dxz_Fy2z;
  abcd[275] = 4.0E0*I_NAI_I3xy2z_Fy2z_aa-2.0E0*2*I_NAI_Gxy2z_Fy2z_a;
  abcd[276] = 4.0E0*I_NAI_I2x4y_Fy2z_aa-2.0E0*3*I_NAI_G2x2y_Fy2z_a-2.0E0*1*I_NAI_G4y_Fy2z_a+3*I_NAI_D2y_Fy2z;
  abcd[277] = 4.0E0*I_NAI_I2x3yz_Fy2z_aa-2.0E0*2*I_NAI_G2xyz_Fy2z_a-2.0E0*1*I_NAI_G3yz_Fy2z_a+2*I_NAI_Dyz_Fy2z;
  abcd[278] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2x2z_Fy2z_a-2.0E0*1*I_NAI_G2y2z_Fy2z_a+1*I_NAI_D2z_Fy2z;
  abcd[279] = 4.0E0*I_NAI_I2xy3z_Fy2z_aa-2.0E0*1*I_NAI_Gy3z_Fy2z_a;
  abcd[280] = 4.0E0*I_NAI_Ix5y_Fy2z_aa-2.0E0*4*I_NAI_Gx3y_Fy2z_a;
  abcd[281] = 4.0E0*I_NAI_Ix4yz_Fy2z_aa-2.0E0*3*I_NAI_Gx2yz_Fy2z_a;
  abcd[282] = 4.0E0*I_NAI_Ix3y2z_Fy2z_aa-2.0E0*2*I_NAI_Gxy2z_Fy2z_a;
  abcd[283] = 4.0E0*I_NAI_Ix2y3z_Fy2z_aa-2.0E0*1*I_NAI_Gx3z_Fy2z_a;
  abcd[284] = 4.0E0*I_NAI_Ixy4z_Fy2z_aa;
  abcd[285] = 4.0E0*I_NAI_I5xy_F3z_aa-2.0E0*4*I_NAI_G3xy_F3z_a;
  abcd[286] = 4.0E0*I_NAI_I4x2y_F3z_aa-2.0E0*1*I_NAI_G4x_F3z_a-2.0E0*3*I_NAI_G2x2y_F3z_a+3*1*I_NAI_D2x_F3z;
  abcd[287] = 4.0E0*I_NAI_I4xyz_F3z_aa-2.0E0*3*I_NAI_G2xyz_F3z_a;
  abcd[288] = 4.0E0*I_NAI_I3x3y_F3z_aa-2.0E0*2*I_NAI_G3xy_F3z_a-2.0E0*2*I_NAI_Gx3y_F3z_a+2*2*I_NAI_Dxy_F3z;
  abcd[289] = 4.0E0*I_NAI_I3x2yz_F3z_aa-2.0E0*1*I_NAI_G3xz_F3z_a-2.0E0*2*I_NAI_Gx2yz_F3z_a+2*1*I_NAI_Dxz_F3z;
  abcd[290] = 4.0E0*I_NAI_I3xy2z_F3z_aa-2.0E0*2*I_NAI_Gxy2z_F3z_a;
  abcd[291] = 4.0E0*I_NAI_I2x4y_F3z_aa-2.0E0*3*I_NAI_G2x2y_F3z_a-2.0E0*1*I_NAI_G4y_F3z_a+3*I_NAI_D2y_F3z;
  abcd[292] = 4.0E0*I_NAI_I2x3yz_F3z_aa-2.0E0*2*I_NAI_G2xyz_F3z_a-2.0E0*1*I_NAI_G3yz_F3z_a+2*I_NAI_Dyz_F3z;
  abcd[293] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2x2z_F3z_a-2.0E0*1*I_NAI_G2y2z_F3z_a+1*I_NAI_D2z_F3z;
  abcd[294] = 4.0E0*I_NAI_I2xy3z_F3z_aa-2.0E0*1*I_NAI_Gy3z_F3z_a;
  abcd[295] = 4.0E0*I_NAI_Ix5y_F3z_aa-2.0E0*4*I_NAI_Gx3y_F3z_a;
  abcd[296] = 4.0E0*I_NAI_Ix4yz_F3z_aa-2.0E0*3*I_NAI_Gx2yz_F3z_a;
  abcd[297] = 4.0E0*I_NAI_Ix3y2z_F3z_aa-2.0E0*2*I_NAI_Gxy2z_F3z_a;
  abcd[298] = 4.0E0*I_NAI_Ix2y3z_F3z_aa-2.0E0*1*I_NAI_Gx3z_F3z_a;
  abcd[299] = 4.0E0*I_NAI_Ixy4z_F3z_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[300] = 4.0E0*I_NAI_I5xz_F3x_aa-2.0E0*4*I_NAI_G3xz_F3x_a;
  abcd[301] = 4.0E0*I_NAI_I4xyz_F3x_aa-2.0E0*3*I_NAI_G2xyz_F3x_a;
  abcd[302] = 4.0E0*I_NAI_I4x2z_F3x_aa-2.0E0*1*I_NAI_G4x_F3x_a-2.0E0*3*I_NAI_G2x2z_F3x_a+3*1*I_NAI_D2x_F3x;
  abcd[303] = 4.0E0*I_NAI_I3x2yz_F3x_aa-2.0E0*2*I_NAI_Gx2yz_F3x_a;
  abcd[304] = 4.0E0*I_NAI_I3xy2z_F3x_aa-2.0E0*1*I_NAI_G3xy_F3x_a-2.0E0*2*I_NAI_Gxy2z_F3x_a+2*1*I_NAI_Dxy_F3x;
  abcd[305] = 4.0E0*I_NAI_I3x3z_F3x_aa-2.0E0*2*I_NAI_G3xz_F3x_a-2.0E0*2*I_NAI_Gx3z_F3x_a+2*2*I_NAI_Dxz_F3x;
  abcd[306] = 4.0E0*I_NAI_I2x3yz_F3x_aa-2.0E0*1*I_NAI_G3yz_F3x_a;
  abcd[307] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2x2y_F3x_a-2.0E0*1*I_NAI_G2y2z_F3x_a+1*I_NAI_D2y_F3x;
  abcd[308] = 4.0E0*I_NAI_I2xy3z_F3x_aa-2.0E0*2*I_NAI_G2xyz_F3x_a-2.0E0*1*I_NAI_Gy3z_F3x_a+2*I_NAI_Dyz_F3x;
  abcd[309] = 4.0E0*I_NAI_I2x4z_F3x_aa-2.0E0*3*I_NAI_G2x2z_F3x_a-2.0E0*1*I_NAI_G4z_F3x_a+3*I_NAI_D2z_F3x;
  abcd[310] = 4.0E0*I_NAI_Ix4yz_F3x_aa;
  abcd[311] = 4.0E0*I_NAI_Ix3y2z_F3x_aa-2.0E0*1*I_NAI_Gx3y_F3x_a;
  abcd[312] = 4.0E0*I_NAI_Ix2y3z_F3x_aa-2.0E0*2*I_NAI_Gx2yz_F3x_a;
  abcd[313] = 4.0E0*I_NAI_Ixy4z_F3x_aa-2.0E0*3*I_NAI_Gxy2z_F3x_a;
  abcd[314] = 4.0E0*I_NAI_Ix5z_F3x_aa-2.0E0*4*I_NAI_Gx3z_F3x_a;
  abcd[315] = 4.0E0*I_NAI_I5xz_F2xy_aa-2.0E0*4*I_NAI_G3xz_F2xy_a;
  abcd[316] = 4.0E0*I_NAI_I4xyz_F2xy_aa-2.0E0*3*I_NAI_G2xyz_F2xy_a;
  abcd[317] = 4.0E0*I_NAI_I4x2z_F2xy_aa-2.0E0*1*I_NAI_G4x_F2xy_a-2.0E0*3*I_NAI_G2x2z_F2xy_a+3*1*I_NAI_D2x_F2xy;
  abcd[318] = 4.0E0*I_NAI_I3x2yz_F2xy_aa-2.0E0*2*I_NAI_Gx2yz_F2xy_a;
  abcd[319] = 4.0E0*I_NAI_I3xy2z_F2xy_aa-2.0E0*1*I_NAI_G3xy_F2xy_a-2.0E0*2*I_NAI_Gxy2z_F2xy_a+2*1*I_NAI_Dxy_F2xy;
  abcd[320] = 4.0E0*I_NAI_I3x3z_F2xy_aa-2.0E0*2*I_NAI_G3xz_F2xy_a-2.0E0*2*I_NAI_Gx3z_F2xy_a+2*2*I_NAI_Dxz_F2xy;
  abcd[321] = 4.0E0*I_NAI_I2x3yz_F2xy_aa-2.0E0*1*I_NAI_G3yz_F2xy_a;
  abcd[322] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2x2y_F2xy_a-2.0E0*1*I_NAI_G2y2z_F2xy_a+1*I_NAI_D2y_F2xy;
  abcd[323] = 4.0E0*I_NAI_I2xy3z_F2xy_aa-2.0E0*2*I_NAI_G2xyz_F2xy_a-2.0E0*1*I_NAI_Gy3z_F2xy_a+2*I_NAI_Dyz_F2xy;
  abcd[324] = 4.0E0*I_NAI_I2x4z_F2xy_aa-2.0E0*3*I_NAI_G2x2z_F2xy_a-2.0E0*1*I_NAI_G4z_F2xy_a+3*I_NAI_D2z_F2xy;
  abcd[325] = 4.0E0*I_NAI_Ix4yz_F2xy_aa;
  abcd[326] = 4.0E0*I_NAI_Ix3y2z_F2xy_aa-2.0E0*1*I_NAI_Gx3y_F2xy_a;
  abcd[327] = 4.0E0*I_NAI_Ix2y3z_F2xy_aa-2.0E0*2*I_NAI_Gx2yz_F2xy_a;
  abcd[328] = 4.0E0*I_NAI_Ixy4z_F2xy_aa-2.0E0*3*I_NAI_Gxy2z_F2xy_a;
  abcd[329] = 4.0E0*I_NAI_Ix5z_F2xy_aa-2.0E0*4*I_NAI_Gx3z_F2xy_a;
  abcd[330] = 4.0E0*I_NAI_I5xz_F2xz_aa-2.0E0*4*I_NAI_G3xz_F2xz_a;
  abcd[331] = 4.0E0*I_NAI_I4xyz_F2xz_aa-2.0E0*3*I_NAI_G2xyz_F2xz_a;
  abcd[332] = 4.0E0*I_NAI_I4x2z_F2xz_aa-2.0E0*1*I_NAI_G4x_F2xz_a-2.0E0*3*I_NAI_G2x2z_F2xz_a+3*1*I_NAI_D2x_F2xz;
  abcd[333] = 4.0E0*I_NAI_I3x2yz_F2xz_aa-2.0E0*2*I_NAI_Gx2yz_F2xz_a;
  abcd[334] = 4.0E0*I_NAI_I3xy2z_F2xz_aa-2.0E0*1*I_NAI_G3xy_F2xz_a-2.0E0*2*I_NAI_Gxy2z_F2xz_a+2*1*I_NAI_Dxy_F2xz;
  abcd[335] = 4.0E0*I_NAI_I3x3z_F2xz_aa-2.0E0*2*I_NAI_G3xz_F2xz_a-2.0E0*2*I_NAI_Gx3z_F2xz_a+2*2*I_NAI_Dxz_F2xz;
  abcd[336] = 4.0E0*I_NAI_I2x3yz_F2xz_aa-2.0E0*1*I_NAI_G3yz_F2xz_a;
  abcd[337] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2x2y_F2xz_a-2.0E0*1*I_NAI_G2y2z_F2xz_a+1*I_NAI_D2y_F2xz;
  abcd[338] = 4.0E0*I_NAI_I2xy3z_F2xz_aa-2.0E0*2*I_NAI_G2xyz_F2xz_a-2.0E0*1*I_NAI_Gy3z_F2xz_a+2*I_NAI_Dyz_F2xz;
  abcd[339] = 4.0E0*I_NAI_I2x4z_F2xz_aa-2.0E0*3*I_NAI_G2x2z_F2xz_a-2.0E0*1*I_NAI_G4z_F2xz_a+3*I_NAI_D2z_F2xz;
  abcd[340] = 4.0E0*I_NAI_Ix4yz_F2xz_aa;
  abcd[341] = 4.0E0*I_NAI_Ix3y2z_F2xz_aa-2.0E0*1*I_NAI_Gx3y_F2xz_a;
  abcd[342] = 4.0E0*I_NAI_Ix2y3z_F2xz_aa-2.0E0*2*I_NAI_Gx2yz_F2xz_a;
  abcd[343] = 4.0E0*I_NAI_Ixy4z_F2xz_aa-2.0E0*3*I_NAI_Gxy2z_F2xz_a;
  abcd[344] = 4.0E0*I_NAI_Ix5z_F2xz_aa-2.0E0*4*I_NAI_Gx3z_F2xz_a;
  abcd[345] = 4.0E0*I_NAI_I5xz_Fx2y_aa-2.0E0*4*I_NAI_G3xz_Fx2y_a;
  abcd[346] = 4.0E0*I_NAI_I4xyz_Fx2y_aa-2.0E0*3*I_NAI_G2xyz_Fx2y_a;
  abcd[347] = 4.0E0*I_NAI_I4x2z_Fx2y_aa-2.0E0*1*I_NAI_G4x_Fx2y_a-2.0E0*3*I_NAI_G2x2z_Fx2y_a+3*1*I_NAI_D2x_Fx2y;
  abcd[348] = 4.0E0*I_NAI_I3x2yz_Fx2y_aa-2.0E0*2*I_NAI_Gx2yz_Fx2y_a;
  abcd[349] = 4.0E0*I_NAI_I3xy2z_Fx2y_aa-2.0E0*1*I_NAI_G3xy_Fx2y_a-2.0E0*2*I_NAI_Gxy2z_Fx2y_a+2*1*I_NAI_Dxy_Fx2y;
  abcd[350] = 4.0E0*I_NAI_I3x3z_Fx2y_aa-2.0E0*2*I_NAI_G3xz_Fx2y_a-2.0E0*2*I_NAI_Gx3z_Fx2y_a+2*2*I_NAI_Dxz_Fx2y;
  abcd[351] = 4.0E0*I_NAI_I2x3yz_Fx2y_aa-2.0E0*1*I_NAI_G3yz_Fx2y_a;
  abcd[352] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2x2y_Fx2y_a-2.0E0*1*I_NAI_G2y2z_Fx2y_a+1*I_NAI_D2y_Fx2y;
  abcd[353] = 4.0E0*I_NAI_I2xy3z_Fx2y_aa-2.0E0*2*I_NAI_G2xyz_Fx2y_a-2.0E0*1*I_NAI_Gy3z_Fx2y_a+2*I_NAI_Dyz_Fx2y;
  abcd[354] = 4.0E0*I_NAI_I2x4z_Fx2y_aa-2.0E0*3*I_NAI_G2x2z_Fx2y_a-2.0E0*1*I_NAI_G4z_Fx2y_a+3*I_NAI_D2z_Fx2y;
  abcd[355] = 4.0E0*I_NAI_Ix4yz_Fx2y_aa;
  abcd[356] = 4.0E0*I_NAI_Ix3y2z_Fx2y_aa-2.0E0*1*I_NAI_Gx3y_Fx2y_a;
  abcd[357] = 4.0E0*I_NAI_Ix2y3z_Fx2y_aa-2.0E0*2*I_NAI_Gx2yz_Fx2y_a;
  abcd[358] = 4.0E0*I_NAI_Ixy4z_Fx2y_aa-2.0E0*3*I_NAI_Gxy2z_Fx2y_a;
  abcd[359] = 4.0E0*I_NAI_Ix5z_Fx2y_aa-2.0E0*4*I_NAI_Gx3z_Fx2y_a;
  abcd[360] = 4.0E0*I_NAI_I5xz_Fxyz_aa-2.0E0*4*I_NAI_G3xz_Fxyz_a;
  abcd[361] = 4.0E0*I_NAI_I4xyz_Fxyz_aa-2.0E0*3*I_NAI_G2xyz_Fxyz_a;
  abcd[362] = 4.0E0*I_NAI_I4x2z_Fxyz_aa-2.0E0*1*I_NAI_G4x_Fxyz_a-2.0E0*3*I_NAI_G2x2z_Fxyz_a+3*1*I_NAI_D2x_Fxyz;
  abcd[363] = 4.0E0*I_NAI_I3x2yz_Fxyz_aa-2.0E0*2*I_NAI_Gx2yz_Fxyz_a;
  abcd[364] = 4.0E0*I_NAI_I3xy2z_Fxyz_aa-2.0E0*1*I_NAI_G3xy_Fxyz_a-2.0E0*2*I_NAI_Gxy2z_Fxyz_a+2*1*I_NAI_Dxy_Fxyz;
  abcd[365] = 4.0E0*I_NAI_I3x3z_Fxyz_aa-2.0E0*2*I_NAI_G3xz_Fxyz_a-2.0E0*2*I_NAI_Gx3z_Fxyz_a+2*2*I_NAI_Dxz_Fxyz;
  abcd[366] = 4.0E0*I_NAI_I2x3yz_Fxyz_aa-2.0E0*1*I_NAI_G3yz_Fxyz_a;
  abcd[367] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2x2y_Fxyz_a-2.0E0*1*I_NAI_G2y2z_Fxyz_a+1*I_NAI_D2y_Fxyz;
  abcd[368] = 4.0E0*I_NAI_I2xy3z_Fxyz_aa-2.0E0*2*I_NAI_G2xyz_Fxyz_a-2.0E0*1*I_NAI_Gy3z_Fxyz_a+2*I_NAI_Dyz_Fxyz;
  abcd[369] = 4.0E0*I_NAI_I2x4z_Fxyz_aa-2.0E0*3*I_NAI_G2x2z_Fxyz_a-2.0E0*1*I_NAI_G4z_Fxyz_a+3*I_NAI_D2z_Fxyz;
  abcd[370] = 4.0E0*I_NAI_Ix4yz_Fxyz_aa;
  abcd[371] = 4.0E0*I_NAI_Ix3y2z_Fxyz_aa-2.0E0*1*I_NAI_Gx3y_Fxyz_a;
  abcd[372] = 4.0E0*I_NAI_Ix2y3z_Fxyz_aa-2.0E0*2*I_NAI_Gx2yz_Fxyz_a;
  abcd[373] = 4.0E0*I_NAI_Ixy4z_Fxyz_aa-2.0E0*3*I_NAI_Gxy2z_Fxyz_a;
  abcd[374] = 4.0E0*I_NAI_Ix5z_Fxyz_aa-2.0E0*4*I_NAI_Gx3z_Fxyz_a;
  abcd[375] = 4.0E0*I_NAI_I5xz_Fx2z_aa-2.0E0*4*I_NAI_G3xz_Fx2z_a;
  abcd[376] = 4.0E0*I_NAI_I4xyz_Fx2z_aa-2.0E0*3*I_NAI_G2xyz_Fx2z_a;
  abcd[377] = 4.0E0*I_NAI_I4x2z_Fx2z_aa-2.0E0*1*I_NAI_G4x_Fx2z_a-2.0E0*3*I_NAI_G2x2z_Fx2z_a+3*1*I_NAI_D2x_Fx2z;
  abcd[378] = 4.0E0*I_NAI_I3x2yz_Fx2z_aa-2.0E0*2*I_NAI_Gx2yz_Fx2z_a;
  abcd[379] = 4.0E0*I_NAI_I3xy2z_Fx2z_aa-2.0E0*1*I_NAI_G3xy_Fx2z_a-2.0E0*2*I_NAI_Gxy2z_Fx2z_a+2*1*I_NAI_Dxy_Fx2z;
  abcd[380] = 4.0E0*I_NAI_I3x3z_Fx2z_aa-2.0E0*2*I_NAI_G3xz_Fx2z_a-2.0E0*2*I_NAI_Gx3z_Fx2z_a+2*2*I_NAI_Dxz_Fx2z;
  abcd[381] = 4.0E0*I_NAI_I2x3yz_Fx2z_aa-2.0E0*1*I_NAI_G3yz_Fx2z_a;
  abcd[382] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2x2y_Fx2z_a-2.0E0*1*I_NAI_G2y2z_Fx2z_a+1*I_NAI_D2y_Fx2z;
  abcd[383] = 4.0E0*I_NAI_I2xy3z_Fx2z_aa-2.0E0*2*I_NAI_G2xyz_Fx2z_a-2.0E0*1*I_NAI_Gy3z_Fx2z_a+2*I_NAI_Dyz_Fx2z;
  abcd[384] = 4.0E0*I_NAI_I2x4z_Fx2z_aa-2.0E0*3*I_NAI_G2x2z_Fx2z_a-2.0E0*1*I_NAI_G4z_Fx2z_a+3*I_NAI_D2z_Fx2z;
  abcd[385] = 4.0E0*I_NAI_Ix4yz_Fx2z_aa;
  abcd[386] = 4.0E0*I_NAI_Ix3y2z_Fx2z_aa-2.0E0*1*I_NAI_Gx3y_Fx2z_a;
  abcd[387] = 4.0E0*I_NAI_Ix2y3z_Fx2z_aa-2.0E0*2*I_NAI_Gx2yz_Fx2z_a;
  abcd[388] = 4.0E0*I_NAI_Ixy4z_Fx2z_aa-2.0E0*3*I_NAI_Gxy2z_Fx2z_a;
  abcd[389] = 4.0E0*I_NAI_Ix5z_Fx2z_aa-2.0E0*4*I_NAI_Gx3z_Fx2z_a;
  abcd[390] = 4.0E0*I_NAI_I5xz_F3y_aa-2.0E0*4*I_NAI_G3xz_F3y_a;
  abcd[391] = 4.0E0*I_NAI_I4xyz_F3y_aa-2.0E0*3*I_NAI_G2xyz_F3y_a;
  abcd[392] = 4.0E0*I_NAI_I4x2z_F3y_aa-2.0E0*1*I_NAI_G4x_F3y_a-2.0E0*3*I_NAI_G2x2z_F3y_a+3*1*I_NAI_D2x_F3y;
  abcd[393] = 4.0E0*I_NAI_I3x2yz_F3y_aa-2.0E0*2*I_NAI_Gx2yz_F3y_a;
  abcd[394] = 4.0E0*I_NAI_I3xy2z_F3y_aa-2.0E0*1*I_NAI_G3xy_F3y_a-2.0E0*2*I_NAI_Gxy2z_F3y_a+2*1*I_NAI_Dxy_F3y;
  abcd[395] = 4.0E0*I_NAI_I3x3z_F3y_aa-2.0E0*2*I_NAI_G3xz_F3y_a-2.0E0*2*I_NAI_Gx3z_F3y_a+2*2*I_NAI_Dxz_F3y;
  abcd[396] = 4.0E0*I_NAI_I2x3yz_F3y_aa-2.0E0*1*I_NAI_G3yz_F3y_a;
  abcd[397] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2x2y_F3y_a-2.0E0*1*I_NAI_G2y2z_F3y_a+1*I_NAI_D2y_F3y;
  abcd[398] = 4.0E0*I_NAI_I2xy3z_F3y_aa-2.0E0*2*I_NAI_G2xyz_F3y_a-2.0E0*1*I_NAI_Gy3z_F3y_a+2*I_NAI_Dyz_F3y;
  abcd[399] = 4.0E0*I_NAI_I2x4z_F3y_aa-2.0E0*3*I_NAI_G2x2z_F3y_a-2.0E0*1*I_NAI_G4z_F3y_a+3*I_NAI_D2z_F3y;
  abcd[400] = 4.0E0*I_NAI_Ix4yz_F3y_aa;
  abcd[401] = 4.0E0*I_NAI_Ix3y2z_F3y_aa-2.0E0*1*I_NAI_Gx3y_F3y_a;
  abcd[402] = 4.0E0*I_NAI_Ix2y3z_F3y_aa-2.0E0*2*I_NAI_Gx2yz_F3y_a;
  abcd[403] = 4.0E0*I_NAI_Ixy4z_F3y_aa-2.0E0*3*I_NAI_Gxy2z_F3y_a;
  abcd[404] = 4.0E0*I_NAI_Ix5z_F3y_aa-2.0E0*4*I_NAI_Gx3z_F3y_a;
  abcd[405] = 4.0E0*I_NAI_I5xz_F2yz_aa-2.0E0*4*I_NAI_G3xz_F2yz_a;
  abcd[406] = 4.0E0*I_NAI_I4xyz_F2yz_aa-2.0E0*3*I_NAI_G2xyz_F2yz_a;
  abcd[407] = 4.0E0*I_NAI_I4x2z_F2yz_aa-2.0E0*1*I_NAI_G4x_F2yz_a-2.0E0*3*I_NAI_G2x2z_F2yz_a+3*1*I_NAI_D2x_F2yz;
  abcd[408] = 4.0E0*I_NAI_I3x2yz_F2yz_aa-2.0E0*2*I_NAI_Gx2yz_F2yz_a;
  abcd[409] = 4.0E0*I_NAI_I3xy2z_F2yz_aa-2.0E0*1*I_NAI_G3xy_F2yz_a-2.0E0*2*I_NAI_Gxy2z_F2yz_a+2*1*I_NAI_Dxy_F2yz;
  abcd[410] = 4.0E0*I_NAI_I3x3z_F2yz_aa-2.0E0*2*I_NAI_G3xz_F2yz_a-2.0E0*2*I_NAI_Gx3z_F2yz_a+2*2*I_NAI_Dxz_F2yz;
  abcd[411] = 4.0E0*I_NAI_I2x3yz_F2yz_aa-2.0E0*1*I_NAI_G3yz_F2yz_a;
  abcd[412] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2x2y_F2yz_a-2.0E0*1*I_NAI_G2y2z_F2yz_a+1*I_NAI_D2y_F2yz;
  abcd[413] = 4.0E0*I_NAI_I2xy3z_F2yz_aa-2.0E0*2*I_NAI_G2xyz_F2yz_a-2.0E0*1*I_NAI_Gy3z_F2yz_a+2*I_NAI_Dyz_F2yz;
  abcd[414] = 4.0E0*I_NAI_I2x4z_F2yz_aa-2.0E0*3*I_NAI_G2x2z_F2yz_a-2.0E0*1*I_NAI_G4z_F2yz_a+3*I_NAI_D2z_F2yz;
  abcd[415] = 4.0E0*I_NAI_Ix4yz_F2yz_aa;
  abcd[416] = 4.0E0*I_NAI_Ix3y2z_F2yz_aa-2.0E0*1*I_NAI_Gx3y_F2yz_a;
  abcd[417] = 4.0E0*I_NAI_Ix2y3z_F2yz_aa-2.0E0*2*I_NAI_Gx2yz_F2yz_a;
  abcd[418] = 4.0E0*I_NAI_Ixy4z_F2yz_aa-2.0E0*3*I_NAI_Gxy2z_F2yz_a;
  abcd[419] = 4.0E0*I_NAI_Ix5z_F2yz_aa-2.0E0*4*I_NAI_Gx3z_F2yz_a;
  abcd[420] = 4.0E0*I_NAI_I5xz_Fy2z_aa-2.0E0*4*I_NAI_G3xz_Fy2z_a;
  abcd[421] = 4.0E0*I_NAI_I4xyz_Fy2z_aa-2.0E0*3*I_NAI_G2xyz_Fy2z_a;
  abcd[422] = 4.0E0*I_NAI_I4x2z_Fy2z_aa-2.0E0*1*I_NAI_G4x_Fy2z_a-2.0E0*3*I_NAI_G2x2z_Fy2z_a+3*1*I_NAI_D2x_Fy2z;
  abcd[423] = 4.0E0*I_NAI_I3x2yz_Fy2z_aa-2.0E0*2*I_NAI_Gx2yz_Fy2z_a;
  abcd[424] = 4.0E0*I_NAI_I3xy2z_Fy2z_aa-2.0E0*1*I_NAI_G3xy_Fy2z_a-2.0E0*2*I_NAI_Gxy2z_Fy2z_a+2*1*I_NAI_Dxy_Fy2z;
  abcd[425] = 4.0E0*I_NAI_I3x3z_Fy2z_aa-2.0E0*2*I_NAI_G3xz_Fy2z_a-2.0E0*2*I_NAI_Gx3z_Fy2z_a+2*2*I_NAI_Dxz_Fy2z;
  abcd[426] = 4.0E0*I_NAI_I2x3yz_Fy2z_aa-2.0E0*1*I_NAI_G3yz_Fy2z_a;
  abcd[427] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2x2y_Fy2z_a-2.0E0*1*I_NAI_G2y2z_Fy2z_a+1*I_NAI_D2y_Fy2z;
  abcd[428] = 4.0E0*I_NAI_I2xy3z_Fy2z_aa-2.0E0*2*I_NAI_G2xyz_Fy2z_a-2.0E0*1*I_NAI_Gy3z_Fy2z_a+2*I_NAI_Dyz_Fy2z;
  abcd[429] = 4.0E0*I_NAI_I2x4z_Fy2z_aa-2.0E0*3*I_NAI_G2x2z_Fy2z_a-2.0E0*1*I_NAI_G4z_Fy2z_a+3*I_NAI_D2z_Fy2z;
  abcd[430] = 4.0E0*I_NAI_Ix4yz_Fy2z_aa;
  abcd[431] = 4.0E0*I_NAI_Ix3y2z_Fy2z_aa-2.0E0*1*I_NAI_Gx3y_Fy2z_a;
  abcd[432] = 4.0E0*I_NAI_Ix2y3z_Fy2z_aa-2.0E0*2*I_NAI_Gx2yz_Fy2z_a;
  abcd[433] = 4.0E0*I_NAI_Ixy4z_Fy2z_aa-2.0E0*3*I_NAI_Gxy2z_Fy2z_a;
  abcd[434] = 4.0E0*I_NAI_Ix5z_Fy2z_aa-2.0E0*4*I_NAI_Gx3z_Fy2z_a;
  abcd[435] = 4.0E0*I_NAI_I5xz_F3z_aa-2.0E0*4*I_NAI_G3xz_F3z_a;
  abcd[436] = 4.0E0*I_NAI_I4xyz_F3z_aa-2.0E0*3*I_NAI_G2xyz_F3z_a;
  abcd[437] = 4.0E0*I_NAI_I4x2z_F3z_aa-2.0E0*1*I_NAI_G4x_F3z_a-2.0E0*3*I_NAI_G2x2z_F3z_a+3*1*I_NAI_D2x_F3z;
  abcd[438] = 4.0E0*I_NAI_I3x2yz_F3z_aa-2.0E0*2*I_NAI_Gx2yz_F3z_a;
  abcd[439] = 4.0E0*I_NAI_I3xy2z_F3z_aa-2.0E0*1*I_NAI_G3xy_F3z_a-2.0E0*2*I_NAI_Gxy2z_F3z_a+2*1*I_NAI_Dxy_F3z;
  abcd[440] = 4.0E0*I_NAI_I3x3z_F3z_aa-2.0E0*2*I_NAI_G3xz_F3z_a-2.0E0*2*I_NAI_Gx3z_F3z_a+2*2*I_NAI_Dxz_F3z;
  abcd[441] = 4.0E0*I_NAI_I2x3yz_F3z_aa-2.0E0*1*I_NAI_G3yz_F3z_a;
  abcd[442] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2x2y_F3z_a-2.0E0*1*I_NAI_G2y2z_F3z_a+1*I_NAI_D2y_F3z;
  abcd[443] = 4.0E0*I_NAI_I2xy3z_F3z_aa-2.0E0*2*I_NAI_G2xyz_F3z_a-2.0E0*1*I_NAI_Gy3z_F3z_a+2*I_NAI_Dyz_F3z;
  abcd[444] = 4.0E0*I_NAI_I2x4z_F3z_aa-2.0E0*3*I_NAI_G2x2z_F3z_a-2.0E0*1*I_NAI_G4z_F3z_a+3*I_NAI_D2z_F3z;
  abcd[445] = 4.0E0*I_NAI_Ix4yz_F3z_aa;
  abcd[446] = 4.0E0*I_NAI_Ix3y2z_F3z_aa-2.0E0*1*I_NAI_Gx3y_F3z_a;
  abcd[447] = 4.0E0*I_NAI_Ix2y3z_F3z_aa-2.0E0*2*I_NAI_Gx2yz_F3z_a;
  abcd[448] = 4.0E0*I_NAI_Ixy4z_F3z_aa-2.0E0*3*I_NAI_Gxy2z_F3z_a;
  abcd[449] = 4.0E0*I_NAI_Ix5z_F3z_aa-2.0E0*4*I_NAI_Gx3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[450] = 4.0E0*I_NAI_I4x2y_F3x_aa-2.0E0*1*I_NAI_G4x_F3x_a;
  abcd[451] = 4.0E0*I_NAI_I3x3y_F3x_aa-2.0E0*1*I_NAI_G3xy_F3x_a-2.0E0*2*I_NAI_G3xy_F3x_a;
  abcd[452] = 4.0E0*I_NAI_I3x2yz_F3x_aa-2.0E0*1*I_NAI_G3xz_F3x_a;
  abcd[453] = 4.0E0*I_NAI_I2x4y_F3x_aa-2.0E0*2*I_NAI_G2x2y_F3x_a-2.0E0*3*I_NAI_G2x2y_F3x_a+2*1*I_NAI_D2x_F3x;
  abcd[454] = 4.0E0*I_NAI_I2x3yz_F3x_aa-2.0E0*1*I_NAI_G2xyz_F3x_a-2.0E0*2*I_NAI_G2xyz_F3x_a;
  abcd[455] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2x2z_F3x_a;
  abcd[456] = 4.0E0*I_NAI_Ix5y_F3x_aa-2.0E0*3*I_NAI_Gx3y_F3x_a-2.0E0*4*I_NAI_Gx3y_F3x_a+3*2*I_NAI_Dxy_F3x;
  abcd[457] = 4.0E0*I_NAI_Ix4yz_F3x_aa-2.0E0*2*I_NAI_Gx2yz_F3x_a-2.0E0*3*I_NAI_Gx2yz_F3x_a+2*1*I_NAI_Dxz_F3x;
  abcd[458] = 4.0E0*I_NAI_Ix3y2z_F3x_aa-2.0E0*1*I_NAI_Gxy2z_F3x_a-2.0E0*2*I_NAI_Gxy2z_F3x_a;
  abcd[459] = 4.0E0*I_NAI_Ix2y3z_F3x_aa-2.0E0*1*I_NAI_Gx3z_F3x_a;
  abcd[460] = 4.0E0*I_NAI_I6y_F3x_aa-2.0E0*4*I_NAI_G4y_F3x_a-2.0E0*5*I_NAI_G4y_F3x_a+4*3*I_NAI_D2y_F3x;
  abcd[461] = 4.0E0*I_NAI_I5yz_F3x_aa-2.0E0*3*I_NAI_G3yz_F3x_a-2.0E0*4*I_NAI_G3yz_F3x_a+3*2*I_NAI_Dyz_F3x;
  abcd[462] = 4.0E0*I_NAI_I4y2z_F3x_aa-2.0E0*2*I_NAI_G2y2z_F3x_a-2.0E0*3*I_NAI_G2y2z_F3x_a+2*1*I_NAI_D2z_F3x;
  abcd[463] = 4.0E0*I_NAI_I3y3z_F3x_aa-2.0E0*1*I_NAI_Gy3z_F3x_a-2.0E0*2*I_NAI_Gy3z_F3x_a;
  abcd[464] = 4.0E0*I_NAI_I2y4z_F3x_aa-2.0E0*1*I_NAI_G4z_F3x_a;
  abcd[465] = 4.0E0*I_NAI_I4x2y_F2xy_aa-2.0E0*1*I_NAI_G4x_F2xy_a;
  abcd[466] = 4.0E0*I_NAI_I3x3y_F2xy_aa-2.0E0*1*I_NAI_G3xy_F2xy_a-2.0E0*2*I_NAI_G3xy_F2xy_a;
  abcd[467] = 4.0E0*I_NAI_I3x2yz_F2xy_aa-2.0E0*1*I_NAI_G3xz_F2xy_a;
  abcd[468] = 4.0E0*I_NAI_I2x4y_F2xy_aa-2.0E0*2*I_NAI_G2x2y_F2xy_a-2.0E0*3*I_NAI_G2x2y_F2xy_a+2*1*I_NAI_D2x_F2xy;
  abcd[469] = 4.0E0*I_NAI_I2x3yz_F2xy_aa-2.0E0*1*I_NAI_G2xyz_F2xy_a-2.0E0*2*I_NAI_G2xyz_F2xy_a;
  abcd[470] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2x2z_F2xy_a;
  abcd[471] = 4.0E0*I_NAI_Ix5y_F2xy_aa-2.0E0*3*I_NAI_Gx3y_F2xy_a-2.0E0*4*I_NAI_Gx3y_F2xy_a+3*2*I_NAI_Dxy_F2xy;
  abcd[472] = 4.0E0*I_NAI_Ix4yz_F2xy_aa-2.0E0*2*I_NAI_Gx2yz_F2xy_a-2.0E0*3*I_NAI_Gx2yz_F2xy_a+2*1*I_NAI_Dxz_F2xy;
  abcd[473] = 4.0E0*I_NAI_Ix3y2z_F2xy_aa-2.0E0*1*I_NAI_Gxy2z_F2xy_a-2.0E0*2*I_NAI_Gxy2z_F2xy_a;
  abcd[474] = 4.0E0*I_NAI_Ix2y3z_F2xy_aa-2.0E0*1*I_NAI_Gx3z_F2xy_a;
  abcd[475] = 4.0E0*I_NAI_I6y_F2xy_aa-2.0E0*4*I_NAI_G4y_F2xy_a-2.0E0*5*I_NAI_G4y_F2xy_a+4*3*I_NAI_D2y_F2xy;
  abcd[476] = 4.0E0*I_NAI_I5yz_F2xy_aa-2.0E0*3*I_NAI_G3yz_F2xy_a-2.0E0*4*I_NAI_G3yz_F2xy_a+3*2*I_NAI_Dyz_F2xy;
  abcd[477] = 4.0E0*I_NAI_I4y2z_F2xy_aa-2.0E0*2*I_NAI_G2y2z_F2xy_a-2.0E0*3*I_NAI_G2y2z_F2xy_a+2*1*I_NAI_D2z_F2xy;
  abcd[478] = 4.0E0*I_NAI_I3y3z_F2xy_aa-2.0E0*1*I_NAI_Gy3z_F2xy_a-2.0E0*2*I_NAI_Gy3z_F2xy_a;
  abcd[479] = 4.0E0*I_NAI_I2y4z_F2xy_aa-2.0E0*1*I_NAI_G4z_F2xy_a;
  abcd[480] = 4.0E0*I_NAI_I4x2y_F2xz_aa-2.0E0*1*I_NAI_G4x_F2xz_a;
  abcd[481] = 4.0E0*I_NAI_I3x3y_F2xz_aa-2.0E0*1*I_NAI_G3xy_F2xz_a-2.0E0*2*I_NAI_G3xy_F2xz_a;
  abcd[482] = 4.0E0*I_NAI_I3x2yz_F2xz_aa-2.0E0*1*I_NAI_G3xz_F2xz_a;
  abcd[483] = 4.0E0*I_NAI_I2x4y_F2xz_aa-2.0E0*2*I_NAI_G2x2y_F2xz_a-2.0E0*3*I_NAI_G2x2y_F2xz_a+2*1*I_NAI_D2x_F2xz;
  abcd[484] = 4.0E0*I_NAI_I2x3yz_F2xz_aa-2.0E0*1*I_NAI_G2xyz_F2xz_a-2.0E0*2*I_NAI_G2xyz_F2xz_a;
  abcd[485] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2x2z_F2xz_a;
  abcd[486] = 4.0E0*I_NAI_Ix5y_F2xz_aa-2.0E0*3*I_NAI_Gx3y_F2xz_a-2.0E0*4*I_NAI_Gx3y_F2xz_a+3*2*I_NAI_Dxy_F2xz;
  abcd[487] = 4.0E0*I_NAI_Ix4yz_F2xz_aa-2.0E0*2*I_NAI_Gx2yz_F2xz_a-2.0E0*3*I_NAI_Gx2yz_F2xz_a+2*1*I_NAI_Dxz_F2xz;
  abcd[488] = 4.0E0*I_NAI_Ix3y2z_F2xz_aa-2.0E0*1*I_NAI_Gxy2z_F2xz_a-2.0E0*2*I_NAI_Gxy2z_F2xz_a;
  abcd[489] = 4.0E0*I_NAI_Ix2y3z_F2xz_aa-2.0E0*1*I_NAI_Gx3z_F2xz_a;
  abcd[490] = 4.0E0*I_NAI_I6y_F2xz_aa-2.0E0*4*I_NAI_G4y_F2xz_a-2.0E0*5*I_NAI_G4y_F2xz_a+4*3*I_NAI_D2y_F2xz;
  abcd[491] = 4.0E0*I_NAI_I5yz_F2xz_aa-2.0E0*3*I_NAI_G3yz_F2xz_a-2.0E0*4*I_NAI_G3yz_F2xz_a+3*2*I_NAI_Dyz_F2xz;
  abcd[492] = 4.0E0*I_NAI_I4y2z_F2xz_aa-2.0E0*2*I_NAI_G2y2z_F2xz_a-2.0E0*3*I_NAI_G2y2z_F2xz_a+2*1*I_NAI_D2z_F2xz;
  abcd[493] = 4.0E0*I_NAI_I3y3z_F2xz_aa-2.0E0*1*I_NAI_Gy3z_F2xz_a-2.0E0*2*I_NAI_Gy3z_F2xz_a;
  abcd[494] = 4.0E0*I_NAI_I2y4z_F2xz_aa-2.0E0*1*I_NAI_G4z_F2xz_a;
  abcd[495] = 4.0E0*I_NAI_I4x2y_Fx2y_aa-2.0E0*1*I_NAI_G4x_Fx2y_a;
  abcd[496] = 4.0E0*I_NAI_I3x3y_Fx2y_aa-2.0E0*1*I_NAI_G3xy_Fx2y_a-2.0E0*2*I_NAI_G3xy_Fx2y_a;
  abcd[497] = 4.0E0*I_NAI_I3x2yz_Fx2y_aa-2.0E0*1*I_NAI_G3xz_Fx2y_a;
  abcd[498] = 4.0E0*I_NAI_I2x4y_Fx2y_aa-2.0E0*2*I_NAI_G2x2y_Fx2y_a-2.0E0*3*I_NAI_G2x2y_Fx2y_a+2*1*I_NAI_D2x_Fx2y;
  abcd[499] = 4.0E0*I_NAI_I2x3yz_Fx2y_aa-2.0E0*1*I_NAI_G2xyz_Fx2y_a-2.0E0*2*I_NAI_G2xyz_Fx2y_a;
  abcd[500] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2x2z_Fx2y_a;
  abcd[501] = 4.0E0*I_NAI_Ix5y_Fx2y_aa-2.0E0*3*I_NAI_Gx3y_Fx2y_a-2.0E0*4*I_NAI_Gx3y_Fx2y_a+3*2*I_NAI_Dxy_Fx2y;
  abcd[502] = 4.0E0*I_NAI_Ix4yz_Fx2y_aa-2.0E0*2*I_NAI_Gx2yz_Fx2y_a-2.0E0*3*I_NAI_Gx2yz_Fx2y_a+2*1*I_NAI_Dxz_Fx2y;
  abcd[503] = 4.0E0*I_NAI_Ix3y2z_Fx2y_aa-2.0E0*1*I_NAI_Gxy2z_Fx2y_a-2.0E0*2*I_NAI_Gxy2z_Fx2y_a;
  abcd[504] = 4.0E0*I_NAI_Ix2y3z_Fx2y_aa-2.0E0*1*I_NAI_Gx3z_Fx2y_a;
  abcd[505] = 4.0E0*I_NAI_I6y_Fx2y_aa-2.0E0*4*I_NAI_G4y_Fx2y_a-2.0E0*5*I_NAI_G4y_Fx2y_a+4*3*I_NAI_D2y_Fx2y;
  abcd[506] = 4.0E0*I_NAI_I5yz_Fx2y_aa-2.0E0*3*I_NAI_G3yz_Fx2y_a-2.0E0*4*I_NAI_G3yz_Fx2y_a+3*2*I_NAI_Dyz_Fx2y;
  abcd[507] = 4.0E0*I_NAI_I4y2z_Fx2y_aa-2.0E0*2*I_NAI_G2y2z_Fx2y_a-2.0E0*3*I_NAI_G2y2z_Fx2y_a+2*1*I_NAI_D2z_Fx2y;
  abcd[508] = 4.0E0*I_NAI_I3y3z_Fx2y_aa-2.0E0*1*I_NAI_Gy3z_Fx2y_a-2.0E0*2*I_NAI_Gy3z_Fx2y_a;
  abcd[509] = 4.0E0*I_NAI_I2y4z_Fx2y_aa-2.0E0*1*I_NAI_G4z_Fx2y_a;
  abcd[510] = 4.0E0*I_NAI_I4x2y_Fxyz_aa-2.0E0*1*I_NAI_G4x_Fxyz_a;
  abcd[511] = 4.0E0*I_NAI_I3x3y_Fxyz_aa-2.0E0*1*I_NAI_G3xy_Fxyz_a-2.0E0*2*I_NAI_G3xy_Fxyz_a;
  abcd[512] = 4.0E0*I_NAI_I3x2yz_Fxyz_aa-2.0E0*1*I_NAI_G3xz_Fxyz_a;
  abcd[513] = 4.0E0*I_NAI_I2x4y_Fxyz_aa-2.0E0*2*I_NAI_G2x2y_Fxyz_a-2.0E0*3*I_NAI_G2x2y_Fxyz_a+2*1*I_NAI_D2x_Fxyz;
  abcd[514] = 4.0E0*I_NAI_I2x3yz_Fxyz_aa-2.0E0*1*I_NAI_G2xyz_Fxyz_a-2.0E0*2*I_NAI_G2xyz_Fxyz_a;
  abcd[515] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2x2z_Fxyz_a;
  abcd[516] = 4.0E0*I_NAI_Ix5y_Fxyz_aa-2.0E0*3*I_NAI_Gx3y_Fxyz_a-2.0E0*4*I_NAI_Gx3y_Fxyz_a+3*2*I_NAI_Dxy_Fxyz;
  abcd[517] = 4.0E0*I_NAI_Ix4yz_Fxyz_aa-2.0E0*2*I_NAI_Gx2yz_Fxyz_a-2.0E0*3*I_NAI_Gx2yz_Fxyz_a+2*1*I_NAI_Dxz_Fxyz;
  abcd[518] = 4.0E0*I_NAI_Ix3y2z_Fxyz_aa-2.0E0*1*I_NAI_Gxy2z_Fxyz_a-2.0E0*2*I_NAI_Gxy2z_Fxyz_a;
  abcd[519] = 4.0E0*I_NAI_Ix2y3z_Fxyz_aa-2.0E0*1*I_NAI_Gx3z_Fxyz_a;
  abcd[520] = 4.0E0*I_NAI_I6y_Fxyz_aa-2.0E0*4*I_NAI_G4y_Fxyz_a-2.0E0*5*I_NAI_G4y_Fxyz_a+4*3*I_NAI_D2y_Fxyz;
  abcd[521] = 4.0E0*I_NAI_I5yz_Fxyz_aa-2.0E0*3*I_NAI_G3yz_Fxyz_a-2.0E0*4*I_NAI_G3yz_Fxyz_a+3*2*I_NAI_Dyz_Fxyz;
  abcd[522] = 4.0E0*I_NAI_I4y2z_Fxyz_aa-2.0E0*2*I_NAI_G2y2z_Fxyz_a-2.0E0*3*I_NAI_G2y2z_Fxyz_a+2*1*I_NAI_D2z_Fxyz;
  abcd[523] = 4.0E0*I_NAI_I3y3z_Fxyz_aa-2.0E0*1*I_NAI_Gy3z_Fxyz_a-2.0E0*2*I_NAI_Gy3z_Fxyz_a;
  abcd[524] = 4.0E0*I_NAI_I2y4z_Fxyz_aa-2.0E0*1*I_NAI_G4z_Fxyz_a;
  abcd[525] = 4.0E0*I_NAI_I4x2y_Fx2z_aa-2.0E0*1*I_NAI_G4x_Fx2z_a;
  abcd[526] = 4.0E0*I_NAI_I3x3y_Fx2z_aa-2.0E0*1*I_NAI_G3xy_Fx2z_a-2.0E0*2*I_NAI_G3xy_Fx2z_a;
  abcd[527] = 4.0E0*I_NAI_I3x2yz_Fx2z_aa-2.0E0*1*I_NAI_G3xz_Fx2z_a;
  abcd[528] = 4.0E0*I_NAI_I2x4y_Fx2z_aa-2.0E0*2*I_NAI_G2x2y_Fx2z_a-2.0E0*3*I_NAI_G2x2y_Fx2z_a+2*1*I_NAI_D2x_Fx2z;
  abcd[529] = 4.0E0*I_NAI_I2x3yz_Fx2z_aa-2.0E0*1*I_NAI_G2xyz_Fx2z_a-2.0E0*2*I_NAI_G2xyz_Fx2z_a;
  abcd[530] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2x2z_Fx2z_a;
  abcd[531] = 4.0E0*I_NAI_Ix5y_Fx2z_aa-2.0E0*3*I_NAI_Gx3y_Fx2z_a-2.0E0*4*I_NAI_Gx3y_Fx2z_a+3*2*I_NAI_Dxy_Fx2z;
  abcd[532] = 4.0E0*I_NAI_Ix4yz_Fx2z_aa-2.0E0*2*I_NAI_Gx2yz_Fx2z_a-2.0E0*3*I_NAI_Gx2yz_Fx2z_a+2*1*I_NAI_Dxz_Fx2z;
  abcd[533] = 4.0E0*I_NAI_Ix3y2z_Fx2z_aa-2.0E0*1*I_NAI_Gxy2z_Fx2z_a-2.0E0*2*I_NAI_Gxy2z_Fx2z_a;
  abcd[534] = 4.0E0*I_NAI_Ix2y3z_Fx2z_aa-2.0E0*1*I_NAI_Gx3z_Fx2z_a;
  abcd[535] = 4.0E0*I_NAI_I6y_Fx2z_aa-2.0E0*4*I_NAI_G4y_Fx2z_a-2.0E0*5*I_NAI_G4y_Fx2z_a+4*3*I_NAI_D2y_Fx2z;
  abcd[536] = 4.0E0*I_NAI_I5yz_Fx2z_aa-2.0E0*3*I_NAI_G3yz_Fx2z_a-2.0E0*4*I_NAI_G3yz_Fx2z_a+3*2*I_NAI_Dyz_Fx2z;
  abcd[537] = 4.0E0*I_NAI_I4y2z_Fx2z_aa-2.0E0*2*I_NAI_G2y2z_Fx2z_a-2.0E0*3*I_NAI_G2y2z_Fx2z_a+2*1*I_NAI_D2z_Fx2z;
  abcd[538] = 4.0E0*I_NAI_I3y3z_Fx2z_aa-2.0E0*1*I_NAI_Gy3z_Fx2z_a-2.0E0*2*I_NAI_Gy3z_Fx2z_a;
  abcd[539] = 4.0E0*I_NAI_I2y4z_Fx2z_aa-2.0E0*1*I_NAI_G4z_Fx2z_a;
  abcd[540] = 4.0E0*I_NAI_I4x2y_F3y_aa-2.0E0*1*I_NAI_G4x_F3y_a;
  abcd[541] = 4.0E0*I_NAI_I3x3y_F3y_aa-2.0E0*1*I_NAI_G3xy_F3y_a-2.0E0*2*I_NAI_G3xy_F3y_a;
  abcd[542] = 4.0E0*I_NAI_I3x2yz_F3y_aa-2.0E0*1*I_NAI_G3xz_F3y_a;
  abcd[543] = 4.0E0*I_NAI_I2x4y_F3y_aa-2.0E0*2*I_NAI_G2x2y_F3y_a-2.0E0*3*I_NAI_G2x2y_F3y_a+2*1*I_NAI_D2x_F3y;
  abcd[544] = 4.0E0*I_NAI_I2x3yz_F3y_aa-2.0E0*1*I_NAI_G2xyz_F3y_a-2.0E0*2*I_NAI_G2xyz_F3y_a;
  abcd[545] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2x2z_F3y_a;
  abcd[546] = 4.0E0*I_NAI_Ix5y_F3y_aa-2.0E0*3*I_NAI_Gx3y_F3y_a-2.0E0*4*I_NAI_Gx3y_F3y_a+3*2*I_NAI_Dxy_F3y;
  abcd[547] = 4.0E0*I_NAI_Ix4yz_F3y_aa-2.0E0*2*I_NAI_Gx2yz_F3y_a-2.0E0*3*I_NAI_Gx2yz_F3y_a+2*1*I_NAI_Dxz_F3y;
  abcd[548] = 4.0E0*I_NAI_Ix3y2z_F3y_aa-2.0E0*1*I_NAI_Gxy2z_F3y_a-2.0E0*2*I_NAI_Gxy2z_F3y_a;
  abcd[549] = 4.0E0*I_NAI_Ix2y3z_F3y_aa-2.0E0*1*I_NAI_Gx3z_F3y_a;
  abcd[550] = 4.0E0*I_NAI_I6y_F3y_aa-2.0E0*4*I_NAI_G4y_F3y_a-2.0E0*5*I_NAI_G4y_F3y_a+4*3*I_NAI_D2y_F3y;
  abcd[551] = 4.0E0*I_NAI_I5yz_F3y_aa-2.0E0*3*I_NAI_G3yz_F3y_a-2.0E0*4*I_NAI_G3yz_F3y_a+3*2*I_NAI_Dyz_F3y;
  abcd[552] = 4.0E0*I_NAI_I4y2z_F3y_aa-2.0E0*2*I_NAI_G2y2z_F3y_a-2.0E0*3*I_NAI_G2y2z_F3y_a+2*1*I_NAI_D2z_F3y;
  abcd[553] = 4.0E0*I_NAI_I3y3z_F3y_aa-2.0E0*1*I_NAI_Gy3z_F3y_a-2.0E0*2*I_NAI_Gy3z_F3y_a;
  abcd[554] = 4.0E0*I_NAI_I2y4z_F3y_aa-2.0E0*1*I_NAI_G4z_F3y_a;
  abcd[555] = 4.0E0*I_NAI_I4x2y_F2yz_aa-2.0E0*1*I_NAI_G4x_F2yz_a;
  abcd[556] = 4.0E0*I_NAI_I3x3y_F2yz_aa-2.0E0*1*I_NAI_G3xy_F2yz_a-2.0E0*2*I_NAI_G3xy_F2yz_a;
  abcd[557] = 4.0E0*I_NAI_I3x2yz_F2yz_aa-2.0E0*1*I_NAI_G3xz_F2yz_a;
  abcd[558] = 4.0E0*I_NAI_I2x4y_F2yz_aa-2.0E0*2*I_NAI_G2x2y_F2yz_a-2.0E0*3*I_NAI_G2x2y_F2yz_a+2*1*I_NAI_D2x_F2yz;
  abcd[559] = 4.0E0*I_NAI_I2x3yz_F2yz_aa-2.0E0*1*I_NAI_G2xyz_F2yz_a-2.0E0*2*I_NAI_G2xyz_F2yz_a;
  abcd[560] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2x2z_F2yz_a;
  abcd[561] = 4.0E0*I_NAI_Ix5y_F2yz_aa-2.0E0*3*I_NAI_Gx3y_F2yz_a-2.0E0*4*I_NAI_Gx3y_F2yz_a+3*2*I_NAI_Dxy_F2yz;
  abcd[562] = 4.0E0*I_NAI_Ix4yz_F2yz_aa-2.0E0*2*I_NAI_Gx2yz_F2yz_a-2.0E0*3*I_NAI_Gx2yz_F2yz_a+2*1*I_NAI_Dxz_F2yz;
  abcd[563] = 4.0E0*I_NAI_Ix3y2z_F2yz_aa-2.0E0*1*I_NAI_Gxy2z_F2yz_a-2.0E0*2*I_NAI_Gxy2z_F2yz_a;
  abcd[564] = 4.0E0*I_NAI_Ix2y3z_F2yz_aa-2.0E0*1*I_NAI_Gx3z_F2yz_a;
  abcd[565] = 4.0E0*I_NAI_I6y_F2yz_aa-2.0E0*4*I_NAI_G4y_F2yz_a-2.0E0*5*I_NAI_G4y_F2yz_a+4*3*I_NAI_D2y_F2yz;
  abcd[566] = 4.0E0*I_NAI_I5yz_F2yz_aa-2.0E0*3*I_NAI_G3yz_F2yz_a-2.0E0*4*I_NAI_G3yz_F2yz_a+3*2*I_NAI_Dyz_F2yz;
  abcd[567] = 4.0E0*I_NAI_I4y2z_F2yz_aa-2.0E0*2*I_NAI_G2y2z_F2yz_a-2.0E0*3*I_NAI_G2y2z_F2yz_a+2*1*I_NAI_D2z_F2yz;
  abcd[568] = 4.0E0*I_NAI_I3y3z_F2yz_aa-2.0E0*1*I_NAI_Gy3z_F2yz_a-2.0E0*2*I_NAI_Gy3z_F2yz_a;
  abcd[569] = 4.0E0*I_NAI_I2y4z_F2yz_aa-2.0E0*1*I_NAI_G4z_F2yz_a;
  abcd[570] = 4.0E0*I_NAI_I4x2y_Fy2z_aa-2.0E0*1*I_NAI_G4x_Fy2z_a;
  abcd[571] = 4.0E0*I_NAI_I3x3y_Fy2z_aa-2.0E0*1*I_NAI_G3xy_Fy2z_a-2.0E0*2*I_NAI_G3xy_Fy2z_a;
  abcd[572] = 4.0E0*I_NAI_I3x2yz_Fy2z_aa-2.0E0*1*I_NAI_G3xz_Fy2z_a;
  abcd[573] = 4.0E0*I_NAI_I2x4y_Fy2z_aa-2.0E0*2*I_NAI_G2x2y_Fy2z_a-2.0E0*3*I_NAI_G2x2y_Fy2z_a+2*1*I_NAI_D2x_Fy2z;
  abcd[574] = 4.0E0*I_NAI_I2x3yz_Fy2z_aa-2.0E0*1*I_NAI_G2xyz_Fy2z_a-2.0E0*2*I_NAI_G2xyz_Fy2z_a;
  abcd[575] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2x2z_Fy2z_a;
  abcd[576] = 4.0E0*I_NAI_Ix5y_Fy2z_aa-2.0E0*3*I_NAI_Gx3y_Fy2z_a-2.0E0*4*I_NAI_Gx3y_Fy2z_a+3*2*I_NAI_Dxy_Fy2z;
  abcd[577] = 4.0E0*I_NAI_Ix4yz_Fy2z_aa-2.0E0*2*I_NAI_Gx2yz_Fy2z_a-2.0E0*3*I_NAI_Gx2yz_Fy2z_a+2*1*I_NAI_Dxz_Fy2z;
  abcd[578] = 4.0E0*I_NAI_Ix3y2z_Fy2z_aa-2.0E0*1*I_NAI_Gxy2z_Fy2z_a-2.0E0*2*I_NAI_Gxy2z_Fy2z_a;
  abcd[579] = 4.0E0*I_NAI_Ix2y3z_Fy2z_aa-2.0E0*1*I_NAI_Gx3z_Fy2z_a;
  abcd[580] = 4.0E0*I_NAI_I6y_Fy2z_aa-2.0E0*4*I_NAI_G4y_Fy2z_a-2.0E0*5*I_NAI_G4y_Fy2z_a+4*3*I_NAI_D2y_Fy2z;
  abcd[581] = 4.0E0*I_NAI_I5yz_Fy2z_aa-2.0E0*3*I_NAI_G3yz_Fy2z_a-2.0E0*4*I_NAI_G3yz_Fy2z_a+3*2*I_NAI_Dyz_Fy2z;
  abcd[582] = 4.0E0*I_NAI_I4y2z_Fy2z_aa-2.0E0*2*I_NAI_G2y2z_Fy2z_a-2.0E0*3*I_NAI_G2y2z_Fy2z_a+2*1*I_NAI_D2z_Fy2z;
  abcd[583] = 4.0E0*I_NAI_I3y3z_Fy2z_aa-2.0E0*1*I_NAI_Gy3z_Fy2z_a-2.0E0*2*I_NAI_Gy3z_Fy2z_a;
  abcd[584] = 4.0E0*I_NAI_I2y4z_Fy2z_aa-2.0E0*1*I_NAI_G4z_Fy2z_a;
  abcd[585] = 4.0E0*I_NAI_I4x2y_F3z_aa-2.0E0*1*I_NAI_G4x_F3z_a;
  abcd[586] = 4.0E0*I_NAI_I3x3y_F3z_aa-2.0E0*1*I_NAI_G3xy_F3z_a-2.0E0*2*I_NAI_G3xy_F3z_a;
  abcd[587] = 4.0E0*I_NAI_I3x2yz_F3z_aa-2.0E0*1*I_NAI_G3xz_F3z_a;
  abcd[588] = 4.0E0*I_NAI_I2x4y_F3z_aa-2.0E0*2*I_NAI_G2x2y_F3z_a-2.0E0*3*I_NAI_G2x2y_F3z_a+2*1*I_NAI_D2x_F3z;
  abcd[589] = 4.0E0*I_NAI_I2x3yz_F3z_aa-2.0E0*1*I_NAI_G2xyz_F3z_a-2.0E0*2*I_NAI_G2xyz_F3z_a;
  abcd[590] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2x2z_F3z_a;
  abcd[591] = 4.0E0*I_NAI_Ix5y_F3z_aa-2.0E0*3*I_NAI_Gx3y_F3z_a-2.0E0*4*I_NAI_Gx3y_F3z_a+3*2*I_NAI_Dxy_F3z;
  abcd[592] = 4.0E0*I_NAI_Ix4yz_F3z_aa-2.0E0*2*I_NAI_Gx2yz_F3z_a-2.0E0*3*I_NAI_Gx2yz_F3z_a+2*1*I_NAI_Dxz_F3z;
  abcd[593] = 4.0E0*I_NAI_Ix3y2z_F3z_aa-2.0E0*1*I_NAI_Gxy2z_F3z_a-2.0E0*2*I_NAI_Gxy2z_F3z_a;
  abcd[594] = 4.0E0*I_NAI_Ix2y3z_F3z_aa-2.0E0*1*I_NAI_Gx3z_F3z_a;
  abcd[595] = 4.0E0*I_NAI_I6y_F3z_aa-2.0E0*4*I_NAI_G4y_F3z_a-2.0E0*5*I_NAI_G4y_F3z_a+4*3*I_NAI_D2y_F3z;
  abcd[596] = 4.0E0*I_NAI_I5yz_F3z_aa-2.0E0*3*I_NAI_G3yz_F3z_a-2.0E0*4*I_NAI_G3yz_F3z_a+3*2*I_NAI_Dyz_F3z;
  abcd[597] = 4.0E0*I_NAI_I4y2z_F3z_aa-2.0E0*2*I_NAI_G2y2z_F3z_a-2.0E0*3*I_NAI_G2y2z_F3z_a+2*1*I_NAI_D2z_F3z;
  abcd[598] = 4.0E0*I_NAI_I3y3z_F3z_aa-2.0E0*1*I_NAI_Gy3z_F3z_a-2.0E0*2*I_NAI_Gy3z_F3z_a;
  abcd[599] = 4.0E0*I_NAI_I2y4z_F3z_aa-2.0E0*1*I_NAI_G4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[600] = 4.0E0*I_NAI_I4xyz_F3x_aa;
  abcd[601] = 4.0E0*I_NAI_I3x2yz_F3x_aa-2.0E0*1*I_NAI_G3xz_F3x_a;
  abcd[602] = 4.0E0*I_NAI_I3xy2z_F3x_aa-2.0E0*1*I_NAI_G3xy_F3x_a;
  abcd[603] = 4.0E0*I_NAI_I2x3yz_F3x_aa-2.0E0*2*I_NAI_G2xyz_F3x_a;
  abcd[604] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2x2y_F3x_a-2.0E0*1*I_NAI_G2x2z_F3x_a+1*I_NAI_D2x_F3x;
  abcd[605] = 4.0E0*I_NAI_I2xy3z_F3x_aa-2.0E0*2*I_NAI_G2xyz_F3x_a;
  abcd[606] = 4.0E0*I_NAI_Ix4yz_F3x_aa-2.0E0*3*I_NAI_Gx2yz_F3x_a;
  abcd[607] = 4.0E0*I_NAI_Ix3y2z_F3x_aa-2.0E0*1*I_NAI_Gx3y_F3x_a-2.0E0*2*I_NAI_Gxy2z_F3x_a+2*1*I_NAI_Dxy_F3x;
  abcd[608] = 4.0E0*I_NAI_Ix2y3z_F3x_aa-2.0E0*2*I_NAI_Gx2yz_F3x_a-2.0E0*1*I_NAI_Gx3z_F3x_a+2*I_NAI_Dxz_F3x;
  abcd[609] = 4.0E0*I_NAI_Ixy4z_F3x_aa-2.0E0*3*I_NAI_Gxy2z_F3x_a;
  abcd[610] = 4.0E0*I_NAI_I5yz_F3x_aa-2.0E0*4*I_NAI_G3yz_F3x_a;
  abcd[611] = 4.0E0*I_NAI_I4y2z_F3x_aa-2.0E0*1*I_NAI_G4y_F3x_a-2.0E0*3*I_NAI_G2y2z_F3x_a+3*1*I_NAI_D2y_F3x;
  abcd[612] = 4.0E0*I_NAI_I3y3z_F3x_aa-2.0E0*2*I_NAI_G3yz_F3x_a-2.0E0*2*I_NAI_Gy3z_F3x_a+2*2*I_NAI_Dyz_F3x;
  abcd[613] = 4.0E0*I_NAI_I2y4z_F3x_aa-2.0E0*3*I_NAI_G2y2z_F3x_a-2.0E0*1*I_NAI_G4z_F3x_a+3*I_NAI_D2z_F3x;
  abcd[614] = 4.0E0*I_NAI_Iy5z_F3x_aa-2.0E0*4*I_NAI_Gy3z_F3x_a;
  abcd[615] = 4.0E0*I_NAI_I4xyz_F2xy_aa;
  abcd[616] = 4.0E0*I_NAI_I3x2yz_F2xy_aa-2.0E0*1*I_NAI_G3xz_F2xy_a;
  abcd[617] = 4.0E0*I_NAI_I3xy2z_F2xy_aa-2.0E0*1*I_NAI_G3xy_F2xy_a;
  abcd[618] = 4.0E0*I_NAI_I2x3yz_F2xy_aa-2.0E0*2*I_NAI_G2xyz_F2xy_a;
  abcd[619] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2x2y_F2xy_a-2.0E0*1*I_NAI_G2x2z_F2xy_a+1*I_NAI_D2x_F2xy;
  abcd[620] = 4.0E0*I_NAI_I2xy3z_F2xy_aa-2.0E0*2*I_NAI_G2xyz_F2xy_a;
  abcd[621] = 4.0E0*I_NAI_Ix4yz_F2xy_aa-2.0E0*3*I_NAI_Gx2yz_F2xy_a;
  abcd[622] = 4.0E0*I_NAI_Ix3y2z_F2xy_aa-2.0E0*1*I_NAI_Gx3y_F2xy_a-2.0E0*2*I_NAI_Gxy2z_F2xy_a+2*1*I_NAI_Dxy_F2xy;
  abcd[623] = 4.0E0*I_NAI_Ix2y3z_F2xy_aa-2.0E0*2*I_NAI_Gx2yz_F2xy_a-2.0E0*1*I_NAI_Gx3z_F2xy_a+2*I_NAI_Dxz_F2xy;
  abcd[624] = 4.0E0*I_NAI_Ixy4z_F2xy_aa-2.0E0*3*I_NAI_Gxy2z_F2xy_a;
  abcd[625] = 4.0E0*I_NAI_I5yz_F2xy_aa-2.0E0*4*I_NAI_G3yz_F2xy_a;
  abcd[626] = 4.0E0*I_NAI_I4y2z_F2xy_aa-2.0E0*1*I_NAI_G4y_F2xy_a-2.0E0*3*I_NAI_G2y2z_F2xy_a+3*1*I_NAI_D2y_F2xy;
  abcd[627] = 4.0E0*I_NAI_I3y3z_F2xy_aa-2.0E0*2*I_NAI_G3yz_F2xy_a-2.0E0*2*I_NAI_Gy3z_F2xy_a+2*2*I_NAI_Dyz_F2xy;
  abcd[628] = 4.0E0*I_NAI_I2y4z_F2xy_aa-2.0E0*3*I_NAI_G2y2z_F2xy_a-2.0E0*1*I_NAI_G4z_F2xy_a+3*I_NAI_D2z_F2xy;
  abcd[629] = 4.0E0*I_NAI_Iy5z_F2xy_aa-2.0E0*4*I_NAI_Gy3z_F2xy_a;
  abcd[630] = 4.0E0*I_NAI_I4xyz_F2xz_aa;
  abcd[631] = 4.0E0*I_NAI_I3x2yz_F2xz_aa-2.0E0*1*I_NAI_G3xz_F2xz_a;
  abcd[632] = 4.0E0*I_NAI_I3xy2z_F2xz_aa-2.0E0*1*I_NAI_G3xy_F2xz_a;
  abcd[633] = 4.0E0*I_NAI_I2x3yz_F2xz_aa-2.0E0*2*I_NAI_G2xyz_F2xz_a;
  abcd[634] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2x2y_F2xz_a-2.0E0*1*I_NAI_G2x2z_F2xz_a+1*I_NAI_D2x_F2xz;
  abcd[635] = 4.0E0*I_NAI_I2xy3z_F2xz_aa-2.0E0*2*I_NAI_G2xyz_F2xz_a;
  abcd[636] = 4.0E0*I_NAI_Ix4yz_F2xz_aa-2.0E0*3*I_NAI_Gx2yz_F2xz_a;
  abcd[637] = 4.0E0*I_NAI_Ix3y2z_F2xz_aa-2.0E0*1*I_NAI_Gx3y_F2xz_a-2.0E0*2*I_NAI_Gxy2z_F2xz_a+2*1*I_NAI_Dxy_F2xz;
  abcd[638] = 4.0E0*I_NAI_Ix2y3z_F2xz_aa-2.0E0*2*I_NAI_Gx2yz_F2xz_a-2.0E0*1*I_NAI_Gx3z_F2xz_a+2*I_NAI_Dxz_F2xz;
  abcd[639] = 4.0E0*I_NAI_Ixy4z_F2xz_aa-2.0E0*3*I_NAI_Gxy2z_F2xz_a;
  abcd[640] = 4.0E0*I_NAI_I5yz_F2xz_aa-2.0E0*4*I_NAI_G3yz_F2xz_a;
  abcd[641] = 4.0E0*I_NAI_I4y2z_F2xz_aa-2.0E0*1*I_NAI_G4y_F2xz_a-2.0E0*3*I_NAI_G2y2z_F2xz_a+3*1*I_NAI_D2y_F2xz;
  abcd[642] = 4.0E0*I_NAI_I3y3z_F2xz_aa-2.0E0*2*I_NAI_G3yz_F2xz_a-2.0E0*2*I_NAI_Gy3z_F2xz_a+2*2*I_NAI_Dyz_F2xz;
  abcd[643] = 4.0E0*I_NAI_I2y4z_F2xz_aa-2.0E0*3*I_NAI_G2y2z_F2xz_a-2.0E0*1*I_NAI_G4z_F2xz_a+3*I_NAI_D2z_F2xz;
  abcd[644] = 4.0E0*I_NAI_Iy5z_F2xz_aa-2.0E0*4*I_NAI_Gy3z_F2xz_a;
  abcd[645] = 4.0E0*I_NAI_I4xyz_Fx2y_aa;
  abcd[646] = 4.0E0*I_NAI_I3x2yz_Fx2y_aa-2.0E0*1*I_NAI_G3xz_Fx2y_a;
  abcd[647] = 4.0E0*I_NAI_I3xy2z_Fx2y_aa-2.0E0*1*I_NAI_G3xy_Fx2y_a;
  abcd[648] = 4.0E0*I_NAI_I2x3yz_Fx2y_aa-2.0E0*2*I_NAI_G2xyz_Fx2y_a;
  abcd[649] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2x2y_Fx2y_a-2.0E0*1*I_NAI_G2x2z_Fx2y_a+1*I_NAI_D2x_Fx2y;
  abcd[650] = 4.0E0*I_NAI_I2xy3z_Fx2y_aa-2.0E0*2*I_NAI_G2xyz_Fx2y_a;
  abcd[651] = 4.0E0*I_NAI_Ix4yz_Fx2y_aa-2.0E0*3*I_NAI_Gx2yz_Fx2y_a;
  abcd[652] = 4.0E0*I_NAI_Ix3y2z_Fx2y_aa-2.0E0*1*I_NAI_Gx3y_Fx2y_a-2.0E0*2*I_NAI_Gxy2z_Fx2y_a+2*1*I_NAI_Dxy_Fx2y;
  abcd[653] = 4.0E0*I_NAI_Ix2y3z_Fx2y_aa-2.0E0*2*I_NAI_Gx2yz_Fx2y_a-2.0E0*1*I_NAI_Gx3z_Fx2y_a+2*I_NAI_Dxz_Fx2y;
  abcd[654] = 4.0E0*I_NAI_Ixy4z_Fx2y_aa-2.0E0*3*I_NAI_Gxy2z_Fx2y_a;
  abcd[655] = 4.0E0*I_NAI_I5yz_Fx2y_aa-2.0E0*4*I_NAI_G3yz_Fx2y_a;
  abcd[656] = 4.0E0*I_NAI_I4y2z_Fx2y_aa-2.0E0*1*I_NAI_G4y_Fx2y_a-2.0E0*3*I_NAI_G2y2z_Fx2y_a+3*1*I_NAI_D2y_Fx2y;
  abcd[657] = 4.0E0*I_NAI_I3y3z_Fx2y_aa-2.0E0*2*I_NAI_G3yz_Fx2y_a-2.0E0*2*I_NAI_Gy3z_Fx2y_a+2*2*I_NAI_Dyz_Fx2y;
  abcd[658] = 4.0E0*I_NAI_I2y4z_Fx2y_aa-2.0E0*3*I_NAI_G2y2z_Fx2y_a-2.0E0*1*I_NAI_G4z_Fx2y_a+3*I_NAI_D2z_Fx2y;
  abcd[659] = 4.0E0*I_NAI_Iy5z_Fx2y_aa-2.0E0*4*I_NAI_Gy3z_Fx2y_a;
  abcd[660] = 4.0E0*I_NAI_I4xyz_Fxyz_aa;
  abcd[661] = 4.0E0*I_NAI_I3x2yz_Fxyz_aa-2.0E0*1*I_NAI_G3xz_Fxyz_a;
  abcd[662] = 4.0E0*I_NAI_I3xy2z_Fxyz_aa-2.0E0*1*I_NAI_G3xy_Fxyz_a;
  abcd[663] = 4.0E0*I_NAI_I2x3yz_Fxyz_aa-2.0E0*2*I_NAI_G2xyz_Fxyz_a;
  abcd[664] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2x2y_Fxyz_a-2.0E0*1*I_NAI_G2x2z_Fxyz_a+1*I_NAI_D2x_Fxyz;
  abcd[665] = 4.0E0*I_NAI_I2xy3z_Fxyz_aa-2.0E0*2*I_NAI_G2xyz_Fxyz_a;
  abcd[666] = 4.0E0*I_NAI_Ix4yz_Fxyz_aa-2.0E0*3*I_NAI_Gx2yz_Fxyz_a;
  abcd[667] = 4.0E0*I_NAI_Ix3y2z_Fxyz_aa-2.0E0*1*I_NAI_Gx3y_Fxyz_a-2.0E0*2*I_NAI_Gxy2z_Fxyz_a+2*1*I_NAI_Dxy_Fxyz;
  abcd[668] = 4.0E0*I_NAI_Ix2y3z_Fxyz_aa-2.0E0*2*I_NAI_Gx2yz_Fxyz_a-2.0E0*1*I_NAI_Gx3z_Fxyz_a+2*I_NAI_Dxz_Fxyz;
  abcd[669] = 4.0E0*I_NAI_Ixy4z_Fxyz_aa-2.0E0*3*I_NAI_Gxy2z_Fxyz_a;
  abcd[670] = 4.0E0*I_NAI_I5yz_Fxyz_aa-2.0E0*4*I_NAI_G3yz_Fxyz_a;
  abcd[671] = 4.0E0*I_NAI_I4y2z_Fxyz_aa-2.0E0*1*I_NAI_G4y_Fxyz_a-2.0E0*3*I_NAI_G2y2z_Fxyz_a+3*1*I_NAI_D2y_Fxyz;
  abcd[672] = 4.0E0*I_NAI_I3y3z_Fxyz_aa-2.0E0*2*I_NAI_G3yz_Fxyz_a-2.0E0*2*I_NAI_Gy3z_Fxyz_a+2*2*I_NAI_Dyz_Fxyz;
  abcd[673] = 4.0E0*I_NAI_I2y4z_Fxyz_aa-2.0E0*3*I_NAI_G2y2z_Fxyz_a-2.0E0*1*I_NAI_G4z_Fxyz_a+3*I_NAI_D2z_Fxyz;
  abcd[674] = 4.0E0*I_NAI_Iy5z_Fxyz_aa-2.0E0*4*I_NAI_Gy3z_Fxyz_a;
  abcd[675] = 4.0E0*I_NAI_I4xyz_Fx2z_aa;
  abcd[676] = 4.0E0*I_NAI_I3x2yz_Fx2z_aa-2.0E0*1*I_NAI_G3xz_Fx2z_a;
  abcd[677] = 4.0E0*I_NAI_I3xy2z_Fx2z_aa-2.0E0*1*I_NAI_G3xy_Fx2z_a;
  abcd[678] = 4.0E0*I_NAI_I2x3yz_Fx2z_aa-2.0E0*2*I_NAI_G2xyz_Fx2z_a;
  abcd[679] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2x2y_Fx2z_a-2.0E0*1*I_NAI_G2x2z_Fx2z_a+1*I_NAI_D2x_Fx2z;
  abcd[680] = 4.0E0*I_NAI_I2xy3z_Fx2z_aa-2.0E0*2*I_NAI_G2xyz_Fx2z_a;
  abcd[681] = 4.0E0*I_NAI_Ix4yz_Fx2z_aa-2.0E0*3*I_NAI_Gx2yz_Fx2z_a;
  abcd[682] = 4.0E0*I_NAI_Ix3y2z_Fx2z_aa-2.0E0*1*I_NAI_Gx3y_Fx2z_a-2.0E0*2*I_NAI_Gxy2z_Fx2z_a+2*1*I_NAI_Dxy_Fx2z;
  abcd[683] = 4.0E0*I_NAI_Ix2y3z_Fx2z_aa-2.0E0*2*I_NAI_Gx2yz_Fx2z_a-2.0E0*1*I_NAI_Gx3z_Fx2z_a+2*I_NAI_Dxz_Fx2z;
  abcd[684] = 4.0E0*I_NAI_Ixy4z_Fx2z_aa-2.0E0*3*I_NAI_Gxy2z_Fx2z_a;
  abcd[685] = 4.0E0*I_NAI_I5yz_Fx2z_aa-2.0E0*4*I_NAI_G3yz_Fx2z_a;
  abcd[686] = 4.0E0*I_NAI_I4y2z_Fx2z_aa-2.0E0*1*I_NAI_G4y_Fx2z_a-2.0E0*3*I_NAI_G2y2z_Fx2z_a+3*1*I_NAI_D2y_Fx2z;
  abcd[687] = 4.0E0*I_NAI_I3y3z_Fx2z_aa-2.0E0*2*I_NAI_G3yz_Fx2z_a-2.0E0*2*I_NAI_Gy3z_Fx2z_a+2*2*I_NAI_Dyz_Fx2z;
  abcd[688] = 4.0E0*I_NAI_I2y4z_Fx2z_aa-2.0E0*3*I_NAI_G2y2z_Fx2z_a-2.0E0*1*I_NAI_G4z_Fx2z_a+3*I_NAI_D2z_Fx2z;
  abcd[689] = 4.0E0*I_NAI_Iy5z_Fx2z_aa-2.0E0*4*I_NAI_Gy3z_Fx2z_a;
  abcd[690] = 4.0E0*I_NAI_I4xyz_F3y_aa;
  abcd[691] = 4.0E0*I_NAI_I3x2yz_F3y_aa-2.0E0*1*I_NAI_G3xz_F3y_a;
  abcd[692] = 4.0E0*I_NAI_I3xy2z_F3y_aa-2.0E0*1*I_NAI_G3xy_F3y_a;
  abcd[693] = 4.0E0*I_NAI_I2x3yz_F3y_aa-2.0E0*2*I_NAI_G2xyz_F3y_a;
  abcd[694] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2x2y_F3y_a-2.0E0*1*I_NAI_G2x2z_F3y_a+1*I_NAI_D2x_F3y;
  abcd[695] = 4.0E0*I_NAI_I2xy3z_F3y_aa-2.0E0*2*I_NAI_G2xyz_F3y_a;
  abcd[696] = 4.0E0*I_NAI_Ix4yz_F3y_aa-2.0E0*3*I_NAI_Gx2yz_F3y_a;
  abcd[697] = 4.0E0*I_NAI_Ix3y2z_F3y_aa-2.0E0*1*I_NAI_Gx3y_F3y_a-2.0E0*2*I_NAI_Gxy2z_F3y_a+2*1*I_NAI_Dxy_F3y;
  abcd[698] = 4.0E0*I_NAI_Ix2y3z_F3y_aa-2.0E0*2*I_NAI_Gx2yz_F3y_a-2.0E0*1*I_NAI_Gx3z_F3y_a+2*I_NAI_Dxz_F3y;
  abcd[699] = 4.0E0*I_NAI_Ixy4z_F3y_aa-2.0E0*3*I_NAI_Gxy2z_F3y_a;
  abcd[700] = 4.0E0*I_NAI_I5yz_F3y_aa-2.0E0*4*I_NAI_G3yz_F3y_a;
  abcd[701] = 4.0E0*I_NAI_I4y2z_F3y_aa-2.0E0*1*I_NAI_G4y_F3y_a-2.0E0*3*I_NAI_G2y2z_F3y_a+3*1*I_NAI_D2y_F3y;
  abcd[702] = 4.0E0*I_NAI_I3y3z_F3y_aa-2.0E0*2*I_NAI_G3yz_F3y_a-2.0E0*2*I_NAI_Gy3z_F3y_a+2*2*I_NAI_Dyz_F3y;
  abcd[703] = 4.0E0*I_NAI_I2y4z_F3y_aa-2.0E0*3*I_NAI_G2y2z_F3y_a-2.0E0*1*I_NAI_G4z_F3y_a+3*I_NAI_D2z_F3y;
  abcd[704] = 4.0E0*I_NAI_Iy5z_F3y_aa-2.0E0*4*I_NAI_Gy3z_F3y_a;
  abcd[705] = 4.0E0*I_NAI_I4xyz_F2yz_aa;
  abcd[706] = 4.0E0*I_NAI_I3x2yz_F2yz_aa-2.0E0*1*I_NAI_G3xz_F2yz_a;
  abcd[707] = 4.0E0*I_NAI_I3xy2z_F2yz_aa-2.0E0*1*I_NAI_G3xy_F2yz_a;
  abcd[708] = 4.0E0*I_NAI_I2x3yz_F2yz_aa-2.0E0*2*I_NAI_G2xyz_F2yz_a;
  abcd[709] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2x2y_F2yz_a-2.0E0*1*I_NAI_G2x2z_F2yz_a+1*I_NAI_D2x_F2yz;
  abcd[710] = 4.0E0*I_NAI_I2xy3z_F2yz_aa-2.0E0*2*I_NAI_G2xyz_F2yz_a;
  abcd[711] = 4.0E0*I_NAI_Ix4yz_F2yz_aa-2.0E0*3*I_NAI_Gx2yz_F2yz_a;
  abcd[712] = 4.0E0*I_NAI_Ix3y2z_F2yz_aa-2.0E0*1*I_NAI_Gx3y_F2yz_a-2.0E0*2*I_NAI_Gxy2z_F2yz_a+2*1*I_NAI_Dxy_F2yz;
  abcd[713] = 4.0E0*I_NAI_Ix2y3z_F2yz_aa-2.0E0*2*I_NAI_Gx2yz_F2yz_a-2.0E0*1*I_NAI_Gx3z_F2yz_a+2*I_NAI_Dxz_F2yz;
  abcd[714] = 4.0E0*I_NAI_Ixy4z_F2yz_aa-2.0E0*3*I_NAI_Gxy2z_F2yz_a;
  abcd[715] = 4.0E0*I_NAI_I5yz_F2yz_aa-2.0E0*4*I_NAI_G3yz_F2yz_a;
  abcd[716] = 4.0E0*I_NAI_I4y2z_F2yz_aa-2.0E0*1*I_NAI_G4y_F2yz_a-2.0E0*3*I_NAI_G2y2z_F2yz_a+3*1*I_NAI_D2y_F2yz;
  abcd[717] = 4.0E0*I_NAI_I3y3z_F2yz_aa-2.0E0*2*I_NAI_G3yz_F2yz_a-2.0E0*2*I_NAI_Gy3z_F2yz_a+2*2*I_NAI_Dyz_F2yz;
  abcd[718] = 4.0E0*I_NAI_I2y4z_F2yz_aa-2.0E0*3*I_NAI_G2y2z_F2yz_a-2.0E0*1*I_NAI_G4z_F2yz_a+3*I_NAI_D2z_F2yz;
  abcd[719] = 4.0E0*I_NAI_Iy5z_F2yz_aa-2.0E0*4*I_NAI_Gy3z_F2yz_a;
  abcd[720] = 4.0E0*I_NAI_I4xyz_Fy2z_aa;
  abcd[721] = 4.0E0*I_NAI_I3x2yz_Fy2z_aa-2.0E0*1*I_NAI_G3xz_Fy2z_a;
  abcd[722] = 4.0E0*I_NAI_I3xy2z_Fy2z_aa-2.0E0*1*I_NAI_G3xy_Fy2z_a;
  abcd[723] = 4.0E0*I_NAI_I2x3yz_Fy2z_aa-2.0E0*2*I_NAI_G2xyz_Fy2z_a;
  abcd[724] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2x2y_Fy2z_a-2.0E0*1*I_NAI_G2x2z_Fy2z_a+1*I_NAI_D2x_Fy2z;
  abcd[725] = 4.0E0*I_NAI_I2xy3z_Fy2z_aa-2.0E0*2*I_NAI_G2xyz_Fy2z_a;
  abcd[726] = 4.0E0*I_NAI_Ix4yz_Fy2z_aa-2.0E0*3*I_NAI_Gx2yz_Fy2z_a;
  abcd[727] = 4.0E0*I_NAI_Ix3y2z_Fy2z_aa-2.0E0*1*I_NAI_Gx3y_Fy2z_a-2.0E0*2*I_NAI_Gxy2z_Fy2z_a+2*1*I_NAI_Dxy_Fy2z;
  abcd[728] = 4.0E0*I_NAI_Ix2y3z_Fy2z_aa-2.0E0*2*I_NAI_Gx2yz_Fy2z_a-2.0E0*1*I_NAI_Gx3z_Fy2z_a+2*I_NAI_Dxz_Fy2z;
  abcd[729] = 4.0E0*I_NAI_Ixy4z_Fy2z_aa-2.0E0*3*I_NAI_Gxy2z_Fy2z_a;
  abcd[730] = 4.0E0*I_NAI_I5yz_Fy2z_aa-2.0E0*4*I_NAI_G3yz_Fy2z_a;
  abcd[731] = 4.0E0*I_NAI_I4y2z_Fy2z_aa-2.0E0*1*I_NAI_G4y_Fy2z_a-2.0E0*3*I_NAI_G2y2z_Fy2z_a+3*1*I_NAI_D2y_Fy2z;
  abcd[732] = 4.0E0*I_NAI_I3y3z_Fy2z_aa-2.0E0*2*I_NAI_G3yz_Fy2z_a-2.0E0*2*I_NAI_Gy3z_Fy2z_a+2*2*I_NAI_Dyz_Fy2z;
  abcd[733] = 4.0E0*I_NAI_I2y4z_Fy2z_aa-2.0E0*3*I_NAI_G2y2z_Fy2z_a-2.0E0*1*I_NAI_G4z_Fy2z_a+3*I_NAI_D2z_Fy2z;
  abcd[734] = 4.0E0*I_NAI_Iy5z_Fy2z_aa-2.0E0*4*I_NAI_Gy3z_Fy2z_a;
  abcd[735] = 4.0E0*I_NAI_I4xyz_F3z_aa;
  abcd[736] = 4.0E0*I_NAI_I3x2yz_F3z_aa-2.0E0*1*I_NAI_G3xz_F3z_a;
  abcd[737] = 4.0E0*I_NAI_I3xy2z_F3z_aa-2.0E0*1*I_NAI_G3xy_F3z_a;
  abcd[738] = 4.0E0*I_NAI_I2x3yz_F3z_aa-2.0E0*2*I_NAI_G2xyz_F3z_a;
  abcd[739] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2x2y_F3z_a-2.0E0*1*I_NAI_G2x2z_F3z_a+1*I_NAI_D2x_F3z;
  abcd[740] = 4.0E0*I_NAI_I2xy3z_F3z_aa-2.0E0*2*I_NAI_G2xyz_F3z_a;
  abcd[741] = 4.0E0*I_NAI_Ix4yz_F3z_aa-2.0E0*3*I_NAI_Gx2yz_F3z_a;
  abcd[742] = 4.0E0*I_NAI_Ix3y2z_F3z_aa-2.0E0*1*I_NAI_Gx3y_F3z_a-2.0E0*2*I_NAI_Gxy2z_F3z_a+2*1*I_NAI_Dxy_F3z;
  abcd[743] = 4.0E0*I_NAI_Ix2y3z_F3z_aa-2.0E0*2*I_NAI_Gx2yz_F3z_a-2.0E0*1*I_NAI_Gx3z_F3z_a+2*I_NAI_Dxz_F3z;
  abcd[744] = 4.0E0*I_NAI_Ixy4z_F3z_aa-2.0E0*3*I_NAI_Gxy2z_F3z_a;
  abcd[745] = 4.0E0*I_NAI_I5yz_F3z_aa-2.0E0*4*I_NAI_G3yz_F3z_a;
  abcd[746] = 4.0E0*I_NAI_I4y2z_F3z_aa-2.0E0*1*I_NAI_G4y_F3z_a-2.0E0*3*I_NAI_G2y2z_F3z_a+3*1*I_NAI_D2y_F3z;
  abcd[747] = 4.0E0*I_NAI_I3y3z_F3z_aa-2.0E0*2*I_NAI_G3yz_F3z_a-2.0E0*2*I_NAI_Gy3z_F3z_a+2*2*I_NAI_Dyz_F3z;
  abcd[748] = 4.0E0*I_NAI_I2y4z_F3z_aa-2.0E0*3*I_NAI_G2y2z_F3z_a-2.0E0*1*I_NAI_G4z_F3z_a+3*I_NAI_D2z_F3z;
  abcd[749] = 4.0E0*I_NAI_Iy5z_F3z_aa-2.0E0*4*I_NAI_Gy3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_aa
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[750] = 4.0E0*I_NAI_I4x2z_F3x_aa-2.0E0*1*I_NAI_G4x_F3x_a;
  abcd[751] = 4.0E0*I_NAI_I3xy2z_F3x_aa-2.0E0*1*I_NAI_G3xy_F3x_a;
  abcd[752] = 4.0E0*I_NAI_I3x3z_F3x_aa-2.0E0*1*I_NAI_G3xz_F3x_a-2.0E0*2*I_NAI_G3xz_F3x_a;
  abcd[753] = 4.0E0*I_NAI_I2x2y2z_F3x_aa-2.0E0*1*I_NAI_G2x2y_F3x_a;
  abcd[754] = 4.0E0*I_NAI_I2xy3z_F3x_aa-2.0E0*1*I_NAI_G2xyz_F3x_a-2.0E0*2*I_NAI_G2xyz_F3x_a;
  abcd[755] = 4.0E0*I_NAI_I2x4z_F3x_aa-2.0E0*2*I_NAI_G2x2z_F3x_a-2.0E0*3*I_NAI_G2x2z_F3x_a+2*1*I_NAI_D2x_F3x;
  abcd[756] = 4.0E0*I_NAI_Ix3y2z_F3x_aa-2.0E0*1*I_NAI_Gx3y_F3x_a;
  abcd[757] = 4.0E0*I_NAI_Ix2y3z_F3x_aa-2.0E0*1*I_NAI_Gx2yz_F3x_a-2.0E0*2*I_NAI_Gx2yz_F3x_a;
  abcd[758] = 4.0E0*I_NAI_Ixy4z_F3x_aa-2.0E0*2*I_NAI_Gxy2z_F3x_a-2.0E0*3*I_NAI_Gxy2z_F3x_a+2*1*I_NAI_Dxy_F3x;
  abcd[759] = 4.0E0*I_NAI_Ix5z_F3x_aa-2.0E0*3*I_NAI_Gx3z_F3x_a-2.0E0*4*I_NAI_Gx3z_F3x_a+3*2*I_NAI_Dxz_F3x;
  abcd[760] = 4.0E0*I_NAI_I4y2z_F3x_aa-2.0E0*1*I_NAI_G4y_F3x_a;
  abcd[761] = 4.0E0*I_NAI_I3y3z_F3x_aa-2.0E0*1*I_NAI_G3yz_F3x_a-2.0E0*2*I_NAI_G3yz_F3x_a;
  abcd[762] = 4.0E0*I_NAI_I2y4z_F3x_aa-2.0E0*2*I_NAI_G2y2z_F3x_a-2.0E0*3*I_NAI_G2y2z_F3x_a+2*1*I_NAI_D2y_F3x;
  abcd[763] = 4.0E0*I_NAI_Iy5z_F3x_aa-2.0E0*3*I_NAI_Gy3z_F3x_a-2.0E0*4*I_NAI_Gy3z_F3x_a+3*2*I_NAI_Dyz_F3x;
  abcd[764] = 4.0E0*I_NAI_I6z_F3x_aa-2.0E0*4*I_NAI_G4z_F3x_a-2.0E0*5*I_NAI_G4z_F3x_a+4*3*I_NAI_D2z_F3x;
  abcd[765] = 4.0E0*I_NAI_I4x2z_F2xy_aa-2.0E0*1*I_NAI_G4x_F2xy_a;
  abcd[766] = 4.0E0*I_NAI_I3xy2z_F2xy_aa-2.0E0*1*I_NAI_G3xy_F2xy_a;
  abcd[767] = 4.0E0*I_NAI_I3x3z_F2xy_aa-2.0E0*1*I_NAI_G3xz_F2xy_a-2.0E0*2*I_NAI_G3xz_F2xy_a;
  abcd[768] = 4.0E0*I_NAI_I2x2y2z_F2xy_aa-2.0E0*1*I_NAI_G2x2y_F2xy_a;
  abcd[769] = 4.0E0*I_NAI_I2xy3z_F2xy_aa-2.0E0*1*I_NAI_G2xyz_F2xy_a-2.0E0*2*I_NAI_G2xyz_F2xy_a;
  abcd[770] = 4.0E0*I_NAI_I2x4z_F2xy_aa-2.0E0*2*I_NAI_G2x2z_F2xy_a-2.0E0*3*I_NAI_G2x2z_F2xy_a+2*1*I_NAI_D2x_F2xy;
  abcd[771] = 4.0E0*I_NAI_Ix3y2z_F2xy_aa-2.0E0*1*I_NAI_Gx3y_F2xy_a;
  abcd[772] = 4.0E0*I_NAI_Ix2y3z_F2xy_aa-2.0E0*1*I_NAI_Gx2yz_F2xy_a-2.0E0*2*I_NAI_Gx2yz_F2xy_a;
  abcd[773] = 4.0E0*I_NAI_Ixy4z_F2xy_aa-2.0E0*2*I_NAI_Gxy2z_F2xy_a-2.0E0*3*I_NAI_Gxy2z_F2xy_a+2*1*I_NAI_Dxy_F2xy;
  abcd[774] = 4.0E0*I_NAI_Ix5z_F2xy_aa-2.0E0*3*I_NAI_Gx3z_F2xy_a-2.0E0*4*I_NAI_Gx3z_F2xy_a+3*2*I_NAI_Dxz_F2xy;
  abcd[775] = 4.0E0*I_NAI_I4y2z_F2xy_aa-2.0E0*1*I_NAI_G4y_F2xy_a;
  abcd[776] = 4.0E0*I_NAI_I3y3z_F2xy_aa-2.0E0*1*I_NAI_G3yz_F2xy_a-2.0E0*2*I_NAI_G3yz_F2xy_a;
  abcd[777] = 4.0E0*I_NAI_I2y4z_F2xy_aa-2.0E0*2*I_NAI_G2y2z_F2xy_a-2.0E0*3*I_NAI_G2y2z_F2xy_a+2*1*I_NAI_D2y_F2xy;
  abcd[778] = 4.0E0*I_NAI_Iy5z_F2xy_aa-2.0E0*3*I_NAI_Gy3z_F2xy_a-2.0E0*4*I_NAI_Gy3z_F2xy_a+3*2*I_NAI_Dyz_F2xy;
  abcd[779] = 4.0E0*I_NAI_I6z_F2xy_aa-2.0E0*4*I_NAI_G4z_F2xy_a-2.0E0*5*I_NAI_G4z_F2xy_a+4*3*I_NAI_D2z_F2xy;
  abcd[780] = 4.0E0*I_NAI_I4x2z_F2xz_aa-2.0E0*1*I_NAI_G4x_F2xz_a;
  abcd[781] = 4.0E0*I_NAI_I3xy2z_F2xz_aa-2.0E0*1*I_NAI_G3xy_F2xz_a;
  abcd[782] = 4.0E0*I_NAI_I3x3z_F2xz_aa-2.0E0*1*I_NAI_G3xz_F2xz_a-2.0E0*2*I_NAI_G3xz_F2xz_a;
  abcd[783] = 4.0E0*I_NAI_I2x2y2z_F2xz_aa-2.0E0*1*I_NAI_G2x2y_F2xz_a;
  abcd[784] = 4.0E0*I_NAI_I2xy3z_F2xz_aa-2.0E0*1*I_NAI_G2xyz_F2xz_a-2.0E0*2*I_NAI_G2xyz_F2xz_a;
  abcd[785] = 4.0E0*I_NAI_I2x4z_F2xz_aa-2.0E0*2*I_NAI_G2x2z_F2xz_a-2.0E0*3*I_NAI_G2x2z_F2xz_a+2*1*I_NAI_D2x_F2xz;
  abcd[786] = 4.0E0*I_NAI_Ix3y2z_F2xz_aa-2.0E0*1*I_NAI_Gx3y_F2xz_a;
  abcd[787] = 4.0E0*I_NAI_Ix2y3z_F2xz_aa-2.0E0*1*I_NAI_Gx2yz_F2xz_a-2.0E0*2*I_NAI_Gx2yz_F2xz_a;
  abcd[788] = 4.0E0*I_NAI_Ixy4z_F2xz_aa-2.0E0*2*I_NAI_Gxy2z_F2xz_a-2.0E0*3*I_NAI_Gxy2z_F2xz_a+2*1*I_NAI_Dxy_F2xz;
  abcd[789] = 4.0E0*I_NAI_Ix5z_F2xz_aa-2.0E0*3*I_NAI_Gx3z_F2xz_a-2.0E0*4*I_NAI_Gx3z_F2xz_a+3*2*I_NAI_Dxz_F2xz;
  abcd[790] = 4.0E0*I_NAI_I4y2z_F2xz_aa-2.0E0*1*I_NAI_G4y_F2xz_a;
  abcd[791] = 4.0E0*I_NAI_I3y3z_F2xz_aa-2.0E0*1*I_NAI_G3yz_F2xz_a-2.0E0*2*I_NAI_G3yz_F2xz_a;
  abcd[792] = 4.0E0*I_NAI_I2y4z_F2xz_aa-2.0E0*2*I_NAI_G2y2z_F2xz_a-2.0E0*3*I_NAI_G2y2z_F2xz_a+2*1*I_NAI_D2y_F2xz;
  abcd[793] = 4.0E0*I_NAI_Iy5z_F2xz_aa-2.0E0*3*I_NAI_Gy3z_F2xz_a-2.0E0*4*I_NAI_Gy3z_F2xz_a+3*2*I_NAI_Dyz_F2xz;
  abcd[794] = 4.0E0*I_NAI_I6z_F2xz_aa-2.0E0*4*I_NAI_G4z_F2xz_a-2.0E0*5*I_NAI_G4z_F2xz_a+4*3*I_NAI_D2z_F2xz;
  abcd[795] = 4.0E0*I_NAI_I4x2z_Fx2y_aa-2.0E0*1*I_NAI_G4x_Fx2y_a;
  abcd[796] = 4.0E0*I_NAI_I3xy2z_Fx2y_aa-2.0E0*1*I_NAI_G3xy_Fx2y_a;
  abcd[797] = 4.0E0*I_NAI_I3x3z_Fx2y_aa-2.0E0*1*I_NAI_G3xz_Fx2y_a-2.0E0*2*I_NAI_G3xz_Fx2y_a;
  abcd[798] = 4.0E0*I_NAI_I2x2y2z_Fx2y_aa-2.0E0*1*I_NAI_G2x2y_Fx2y_a;
  abcd[799] = 4.0E0*I_NAI_I2xy3z_Fx2y_aa-2.0E0*1*I_NAI_G2xyz_Fx2y_a-2.0E0*2*I_NAI_G2xyz_Fx2y_a;
  abcd[800] = 4.0E0*I_NAI_I2x4z_Fx2y_aa-2.0E0*2*I_NAI_G2x2z_Fx2y_a-2.0E0*3*I_NAI_G2x2z_Fx2y_a+2*1*I_NAI_D2x_Fx2y;
  abcd[801] = 4.0E0*I_NAI_Ix3y2z_Fx2y_aa-2.0E0*1*I_NAI_Gx3y_Fx2y_a;
  abcd[802] = 4.0E0*I_NAI_Ix2y3z_Fx2y_aa-2.0E0*1*I_NAI_Gx2yz_Fx2y_a-2.0E0*2*I_NAI_Gx2yz_Fx2y_a;
  abcd[803] = 4.0E0*I_NAI_Ixy4z_Fx2y_aa-2.0E0*2*I_NAI_Gxy2z_Fx2y_a-2.0E0*3*I_NAI_Gxy2z_Fx2y_a+2*1*I_NAI_Dxy_Fx2y;
  abcd[804] = 4.0E0*I_NAI_Ix5z_Fx2y_aa-2.0E0*3*I_NAI_Gx3z_Fx2y_a-2.0E0*4*I_NAI_Gx3z_Fx2y_a+3*2*I_NAI_Dxz_Fx2y;
  abcd[805] = 4.0E0*I_NAI_I4y2z_Fx2y_aa-2.0E0*1*I_NAI_G4y_Fx2y_a;
  abcd[806] = 4.0E0*I_NAI_I3y3z_Fx2y_aa-2.0E0*1*I_NAI_G3yz_Fx2y_a-2.0E0*2*I_NAI_G3yz_Fx2y_a;
  abcd[807] = 4.0E0*I_NAI_I2y4z_Fx2y_aa-2.0E0*2*I_NAI_G2y2z_Fx2y_a-2.0E0*3*I_NAI_G2y2z_Fx2y_a+2*1*I_NAI_D2y_Fx2y;
  abcd[808] = 4.0E0*I_NAI_Iy5z_Fx2y_aa-2.0E0*3*I_NAI_Gy3z_Fx2y_a-2.0E0*4*I_NAI_Gy3z_Fx2y_a+3*2*I_NAI_Dyz_Fx2y;
  abcd[809] = 4.0E0*I_NAI_I6z_Fx2y_aa-2.0E0*4*I_NAI_G4z_Fx2y_a-2.0E0*5*I_NAI_G4z_Fx2y_a+4*3*I_NAI_D2z_Fx2y;
  abcd[810] = 4.0E0*I_NAI_I4x2z_Fxyz_aa-2.0E0*1*I_NAI_G4x_Fxyz_a;
  abcd[811] = 4.0E0*I_NAI_I3xy2z_Fxyz_aa-2.0E0*1*I_NAI_G3xy_Fxyz_a;
  abcd[812] = 4.0E0*I_NAI_I3x3z_Fxyz_aa-2.0E0*1*I_NAI_G3xz_Fxyz_a-2.0E0*2*I_NAI_G3xz_Fxyz_a;
  abcd[813] = 4.0E0*I_NAI_I2x2y2z_Fxyz_aa-2.0E0*1*I_NAI_G2x2y_Fxyz_a;
  abcd[814] = 4.0E0*I_NAI_I2xy3z_Fxyz_aa-2.0E0*1*I_NAI_G2xyz_Fxyz_a-2.0E0*2*I_NAI_G2xyz_Fxyz_a;
  abcd[815] = 4.0E0*I_NAI_I2x4z_Fxyz_aa-2.0E0*2*I_NAI_G2x2z_Fxyz_a-2.0E0*3*I_NAI_G2x2z_Fxyz_a+2*1*I_NAI_D2x_Fxyz;
  abcd[816] = 4.0E0*I_NAI_Ix3y2z_Fxyz_aa-2.0E0*1*I_NAI_Gx3y_Fxyz_a;
  abcd[817] = 4.0E0*I_NAI_Ix2y3z_Fxyz_aa-2.0E0*1*I_NAI_Gx2yz_Fxyz_a-2.0E0*2*I_NAI_Gx2yz_Fxyz_a;
  abcd[818] = 4.0E0*I_NAI_Ixy4z_Fxyz_aa-2.0E0*2*I_NAI_Gxy2z_Fxyz_a-2.0E0*3*I_NAI_Gxy2z_Fxyz_a+2*1*I_NAI_Dxy_Fxyz;
  abcd[819] = 4.0E0*I_NAI_Ix5z_Fxyz_aa-2.0E0*3*I_NAI_Gx3z_Fxyz_a-2.0E0*4*I_NAI_Gx3z_Fxyz_a+3*2*I_NAI_Dxz_Fxyz;
  abcd[820] = 4.0E0*I_NAI_I4y2z_Fxyz_aa-2.0E0*1*I_NAI_G4y_Fxyz_a;
  abcd[821] = 4.0E0*I_NAI_I3y3z_Fxyz_aa-2.0E0*1*I_NAI_G3yz_Fxyz_a-2.0E0*2*I_NAI_G3yz_Fxyz_a;
  abcd[822] = 4.0E0*I_NAI_I2y4z_Fxyz_aa-2.0E0*2*I_NAI_G2y2z_Fxyz_a-2.0E0*3*I_NAI_G2y2z_Fxyz_a+2*1*I_NAI_D2y_Fxyz;
  abcd[823] = 4.0E0*I_NAI_Iy5z_Fxyz_aa-2.0E0*3*I_NAI_Gy3z_Fxyz_a-2.0E0*4*I_NAI_Gy3z_Fxyz_a+3*2*I_NAI_Dyz_Fxyz;
  abcd[824] = 4.0E0*I_NAI_I6z_Fxyz_aa-2.0E0*4*I_NAI_G4z_Fxyz_a-2.0E0*5*I_NAI_G4z_Fxyz_a+4*3*I_NAI_D2z_Fxyz;
  abcd[825] = 4.0E0*I_NAI_I4x2z_Fx2z_aa-2.0E0*1*I_NAI_G4x_Fx2z_a;
  abcd[826] = 4.0E0*I_NAI_I3xy2z_Fx2z_aa-2.0E0*1*I_NAI_G3xy_Fx2z_a;
  abcd[827] = 4.0E0*I_NAI_I3x3z_Fx2z_aa-2.0E0*1*I_NAI_G3xz_Fx2z_a-2.0E0*2*I_NAI_G3xz_Fx2z_a;
  abcd[828] = 4.0E0*I_NAI_I2x2y2z_Fx2z_aa-2.0E0*1*I_NAI_G2x2y_Fx2z_a;
  abcd[829] = 4.0E0*I_NAI_I2xy3z_Fx2z_aa-2.0E0*1*I_NAI_G2xyz_Fx2z_a-2.0E0*2*I_NAI_G2xyz_Fx2z_a;
  abcd[830] = 4.0E0*I_NAI_I2x4z_Fx2z_aa-2.0E0*2*I_NAI_G2x2z_Fx2z_a-2.0E0*3*I_NAI_G2x2z_Fx2z_a+2*1*I_NAI_D2x_Fx2z;
  abcd[831] = 4.0E0*I_NAI_Ix3y2z_Fx2z_aa-2.0E0*1*I_NAI_Gx3y_Fx2z_a;
  abcd[832] = 4.0E0*I_NAI_Ix2y3z_Fx2z_aa-2.0E0*1*I_NAI_Gx2yz_Fx2z_a-2.0E0*2*I_NAI_Gx2yz_Fx2z_a;
  abcd[833] = 4.0E0*I_NAI_Ixy4z_Fx2z_aa-2.0E0*2*I_NAI_Gxy2z_Fx2z_a-2.0E0*3*I_NAI_Gxy2z_Fx2z_a+2*1*I_NAI_Dxy_Fx2z;
  abcd[834] = 4.0E0*I_NAI_Ix5z_Fx2z_aa-2.0E0*3*I_NAI_Gx3z_Fx2z_a-2.0E0*4*I_NAI_Gx3z_Fx2z_a+3*2*I_NAI_Dxz_Fx2z;
  abcd[835] = 4.0E0*I_NAI_I4y2z_Fx2z_aa-2.0E0*1*I_NAI_G4y_Fx2z_a;
  abcd[836] = 4.0E0*I_NAI_I3y3z_Fx2z_aa-2.0E0*1*I_NAI_G3yz_Fx2z_a-2.0E0*2*I_NAI_G3yz_Fx2z_a;
  abcd[837] = 4.0E0*I_NAI_I2y4z_Fx2z_aa-2.0E0*2*I_NAI_G2y2z_Fx2z_a-2.0E0*3*I_NAI_G2y2z_Fx2z_a+2*1*I_NAI_D2y_Fx2z;
  abcd[838] = 4.0E0*I_NAI_Iy5z_Fx2z_aa-2.0E0*3*I_NAI_Gy3z_Fx2z_a-2.0E0*4*I_NAI_Gy3z_Fx2z_a+3*2*I_NAI_Dyz_Fx2z;
  abcd[839] = 4.0E0*I_NAI_I6z_Fx2z_aa-2.0E0*4*I_NAI_G4z_Fx2z_a-2.0E0*5*I_NAI_G4z_Fx2z_a+4*3*I_NAI_D2z_Fx2z;
  abcd[840] = 4.0E0*I_NAI_I4x2z_F3y_aa-2.0E0*1*I_NAI_G4x_F3y_a;
  abcd[841] = 4.0E0*I_NAI_I3xy2z_F3y_aa-2.0E0*1*I_NAI_G3xy_F3y_a;
  abcd[842] = 4.0E0*I_NAI_I3x3z_F3y_aa-2.0E0*1*I_NAI_G3xz_F3y_a-2.0E0*2*I_NAI_G3xz_F3y_a;
  abcd[843] = 4.0E0*I_NAI_I2x2y2z_F3y_aa-2.0E0*1*I_NAI_G2x2y_F3y_a;
  abcd[844] = 4.0E0*I_NAI_I2xy3z_F3y_aa-2.0E0*1*I_NAI_G2xyz_F3y_a-2.0E0*2*I_NAI_G2xyz_F3y_a;
  abcd[845] = 4.0E0*I_NAI_I2x4z_F3y_aa-2.0E0*2*I_NAI_G2x2z_F3y_a-2.0E0*3*I_NAI_G2x2z_F3y_a+2*1*I_NAI_D2x_F3y;
  abcd[846] = 4.0E0*I_NAI_Ix3y2z_F3y_aa-2.0E0*1*I_NAI_Gx3y_F3y_a;
  abcd[847] = 4.0E0*I_NAI_Ix2y3z_F3y_aa-2.0E0*1*I_NAI_Gx2yz_F3y_a-2.0E0*2*I_NAI_Gx2yz_F3y_a;
  abcd[848] = 4.0E0*I_NAI_Ixy4z_F3y_aa-2.0E0*2*I_NAI_Gxy2z_F3y_a-2.0E0*3*I_NAI_Gxy2z_F3y_a+2*1*I_NAI_Dxy_F3y;
  abcd[849] = 4.0E0*I_NAI_Ix5z_F3y_aa-2.0E0*3*I_NAI_Gx3z_F3y_a-2.0E0*4*I_NAI_Gx3z_F3y_a+3*2*I_NAI_Dxz_F3y;
  abcd[850] = 4.0E0*I_NAI_I4y2z_F3y_aa-2.0E0*1*I_NAI_G4y_F3y_a;
  abcd[851] = 4.0E0*I_NAI_I3y3z_F3y_aa-2.0E0*1*I_NAI_G3yz_F3y_a-2.0E0*2*I_NAI_G3yz_F3y_a;
  abcd[852] = 4.0E0*I_NAI_I2y4z_F3y_aa-2.0E0*2*I_NAI_G2y2z_F3y_a-2.0E0*3*I_NAI_G2y2z_F3y_a+2*1*I_NAI_D2y_F3y;
  abcd[853] = 4.0E0*I_NAI_Iy5z_F3y_aa-2.0E0*3*I_NAI_Gy3z_F3y_a-2.0E0*4*I_NAI_Gy3z_F3y_a+3*2*I_NAI_Dyz_F3y;
  abcd[854] = 4.0E0*I_NAI_I6z_F3y_aa-2.0E0*4*I_NAI_G4z_F3y_a-2.0E0*5*I_NAI_G4z_F3y_a+4*3*I_NAI_D2z_F3y;
  abcd[855] = 4.0E0*I_NAI_I4x2z_F2yz_aa-2.0E0*1*I_NAI_G4x_F2yz_a;
  abcd[856] = 4.0E0*I_NAI_I3xy2z_F2yz_aa-2.0E0*1*I_NAI_G3xy_F2yz_a;
  abcd[857] = 4.0E0*I_NAI_I3x3z_F2yz_aa-2.0E0*1*I_NAI_G3xz_F2yz_a-2.0E0*2*I_NAI_G3xz_F2yz_a;
  abcd[858] = 4.0E0*I_NAI_I2x2y2z_F2yz_aa-2.0E0*1*I_NAI_G2x2y_F2yz_a;
  abcd[859] = 4.0E0*I_NAI_I2xy3z_F2yz_aa-2.0E0*1*I_NAI_G2xyz_F2yz_a-2.0E0*2*I_NAI_G2xyz_F2yz_a;
  abcd[860] = 4.0E0*I_NAI_I2x4z_F2yz_aa-2.0E0*2*I_NAI_G2x2z_F2yz_a-2.0E0*3*I_NAI_G2x2z_F2yz_a+2*1*I_NAI_D2x_F2yz;
  abcd[861] = 4.0E0*I_NAI_Ix3y2z_F2yz_aa-2.0E0*1*I_NAI_Gx3y_F2yz_a;
  abcd[862] = 4.0E0*I_NAI_Ix2y3z_F2yz_aa-2.0E0*1*I_NAI_Gx2yz_F2yz_a-2.0E0*2*I_NAI_Gx2yz_F2yz_a;
  abcd[863] = 4.0E0*I_NAI_Ixy4z_F2yz_aa-2.0E0*2*I_NAI_Gxy2z_F2yz_a-2.0E0*3*I_NAI_Gxy2z_F2yz_a+2*1*I_NAI_Dxy_F2yz;
  abcd[864] = 4.0E0*I_NAI_Ix5z_F2yz_aa-2.0E0*3*I_NAI_Gx3z_F2yz_a-2.0E0*4*I_NAI_Gx3z_F2yz_a+3*2*I_NAI_Dxz_F2yz;
  abcd[865] = 4.0E0*I_NAI_I4y2z_F2yz_aa-2.0E0*1*I_NAI_G4y_F2yz_a;
  abcd[866] = 4.0E0*I_NAI_I3y3z_F2yz_aa-2.0E0*1*I_NAI_G3yz_F2yz_a-2.0E0*2*I_NAI_G3yz_F2yz_a;
  abcd[867] = 4.0E0*I_NAI_I2y4z_F2yz_aa-2.0E0*2*I_NAI_G2y2z_F2yz_a-2.0E0*3*I_NAI_G2y2z_F2yz_a+2*1*I_NAI_D2y_F2yz;
  abcd[868] = 4.0E0*I_NAI_Iy5z_F2yz_aa-2.0E0*3*I_NAI_Gy3z_F2yz_a-2.0E0*4*I_NAI_Gy3z_F2yz_a+3*2*I_NAI_Dyz_F2yz;
  abcd[869] = 4.0E0*I_NAI_I6z_F2yz_aa-2.0E0*4*I_NAI_G4z_F2yz_a-2.0E0*5*I_NAI_G4z_F2yz_a+4*3*I_NAI_D2z_F2yz;
  abcd[870] = 4.0E0*I_NAI_I4x2z_Fy2z_aa-2.0E0*1*I_NAI_G4x_Fy2z_a;
  abcd[871] = 4.0E0*I_NAI_I3xy2z_Fy2z_aa-2.0E0*1*I_NAI_G3xy_Fy2z_a;
  abcd[872] = 4.0E0*I_NAI_I3x3z_Fy2z_aa-2.0E0*1*I_NAI_G3xz_Fy2z_a-2.0E0*2*I_NAI_G3xz_Fy2z_a;
  abcd[873] = 4.0E0*I_NAI_I2x2y2z_Fy2z_aa-2.0E0*1*I_NAI_G2x2y_Fy2z_a;
  abcd[874] = 4.0E0*I_NAI_I2xy3z_Fy2z_aa-2.0E0*1*I_NAI_G2xyz_Fy2z_a-2.0E0*2*I_NAI_G2xyz_Fy2z_a;
  abcd[875] = 4.0E0*I_NAI_I2x4z_Fy2z_aa-2.0E0*2*I_NAI_G2x2z_Fy2z_a-2.0E0*3*I_NAI_G2x2z_Fy2z_a+2*1*I_NAI_D2x_Fy2z;
  abcd[876] = 4.0E0*I_NAI_Ix3y2z_Fy2z_aa-2.0E0*1*I_NAI_Gx3y_Fy2z_a;
  abcd[877] = 4.0E0*I_NAI_Ix2y3z_Fy2z_aa-2.0E0*1*I_NAI_Gx2yz_Fy2z_a-2.0E0*2*I_NAI_Gx2yz_Fy2z_a;
  abcd[878] = 4.0E0*I_NAI_Ixy4z_Fy2z_aa-2.0E0*2*I_NAI_Gxy2z_Fy2z_a-2.0E0*3*I_NAI_Gxy2z_Fy2z_a+2*1*I_NAI_Dxy_Fy2z;
  abcd[879] = 4.0E0*I_NAI_Ix5z_Fy2z_aa-2.0E0*3*I_NAI_Gx3z_Fy2z_a-2.0E0*4*I_NAI_Gx3z_Fy2z_a+3*2*I_NAI_Dxz_Fy2z;
  abcd[880] = 4.0E0*I_NAI_I4y2z_Fy2z_aa-2.0E0*1*I_NAI_G4y_Fy2z_a;
  abcd[881] = 4.0E0*I_NAI_I3y3z_Fy2z_aa-2.0E0*1*I_NAI_G3yz_Fy2z_a-2.0E0*2*I_NAI_G3yz_Fy2z_a;
  abcd[882] = 4.0E0*I_NAI_I2y4z_Fy2z_aa-2.0E0*2*I_NAI_G2y2z_Fy2z_a-2.0E0*3*I_NAI_G2y2z_Fy2z_a+2*1*I_NAI_D2y_Fy2z;
  abcd[883] = 4.0E0*I_NAI_Iy5z_Fy2z_aa-2.0E0*3*I_NAI_Gy3z_Fy2z_a-2.0E0*4*I_NAI_Gy3z_Fy2z_a+3*2*I_NAI_Dyz_Fy2z;
  abcd[884] = 4.0E0*I_NAI_I6z_Fy2z_aa-2.0E0*4*I_NAI_G4z_Fy2z_a-2.0E0*5*I_NAI_G4z_Fy2z_a+4*3*I_NAI_D2z_Fy2z;
  abcd[885] = 4.0E0*I_NAI_I4x2z_F3z_aa-2.0E0*1*I_NAI_G4x_F3z_a;
  abcd[886] = 4.0E0*I_NAI_I3xy2z_F3z_aa-2.0E0*1*I_NAI_G3xy_F3z_a;
  abcd[887] = 4.0E0*I_NAI_I3x3z_F3z_aa-2.0E0*1*I_NAI_G3xz_F3z_a-2.0E0*2*I_NAI_G3xz_F3z_a;
  abcd[888] = 4.0E0*I_NAI_I2x2y2z_F3z_aa-2.0E0*1*I_NAI_G2x2y_F3z_a;
  abcd[889] = 4.0E0*I_NAI_I2xy3z_F3z_aa-2.0E0*1*I_NAI_G2xyz_F3z_a-2.0E0*2*I_NAI_G2xyz_F3z_a;
  abcd[890] = 4.0E0*I_NAI_I2x4z_F3z_aa-2.0E0*2*I_NAI_G2x2z_F3z_a-2.0E0*3*I_NAI_G2x2z_F3z_a+2*1*I_NAI_D2x_F3z;
  abcd[891] = 4.0E0*I_NAI_Ix3y2z_F3z_aa-2.0E0*1*I_NAI_Gx3y_F3z_a;
  abcd[892] = 4.0E0*I_NAI_Ix2y3z_F3z_aa-2.0E0*1*I_NAI_Gx2yz_F3z_a-2.0E0*2*I_NAI_Gx2yz_F3z_a;
  abcd[893] = 4.0E0*I_NAI_Ixy4z_F3z_aa-2.0E0*2*I_NAI_Gxy2z_F3z_a-2.0E0*3*I_NAI_Gxy2z_F3z_a+2*1*I_NAI_Dxy_F3z;
  abcd[894] = 4.0E0*I_NAI_Ix5z_F3z_aa-2.0E0*3*I_NAI_Gx3z_F3z_a-2.0E0*4*I_NAI_Gx3z_F3z_a+3*2*I_NAI_Dxz_F3z;
  abcd[895] = 4.0E0*I_NAI_I4y2z_F3z_aa-2.0E0*1*I_NAI_G4y_F3z_a;
  abcd[896] = 4.0E0*I_NAI_I3y3z_F3z_aa-2.0E0*1*I_NAI_G3yz_F3z_a-2.0E0*2*I_NAI_G3yz_F3z_a;
  abcd[897] = 4.0E0*I_NAI_I2y4z_F3z_aa-2.0E0*2*I_NAI_G2y2z_F3z_a-2.0E0*3*I_NAI_G2y2z_F3z_a+2*1*I_NAI_D2y_F3z;
  abcd[898] = 4.0E0*I_NAI_Iy5z_F3z_aa-2.0E0*3*I_NAI_Gy3z_F3z_a-2.0E0*4*I_NAI_Gy3z_F3z_a+3*2*I_NAI_Dyz_F3z;
  abcd[899] = 4.0E0*I_NAI_I6z_F3z_aa-2.0E0*4*I_NAI_G4z_F3z_a-2.0E0*5*I_NAI_G4z_F3z_a+4*3*I_NAI_D2z_F3z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[900] = 4.0E0*I_NAI_H5x_G4x_ab-2.0E0*3*I_NAI_H5x_D2x_a-2.0E0*4*I_NAI_F3x_G4x_b+4*3*I_NAI_F3x_D2x;
  abcd[901] = 4.0E0*I_NAI_H4xy_G4x_ab-2.0E0*3*I_NAI_H4xy_D2x_a-2.0E0*3*I_NAI_F2xy_G4x_b+3*3*I_NAI_F2xy_D2x;
  abcd[902] = 4.0E0*I_NAI_H4xz_G4x_ab-2.0E0*3*I_NAI_H4xz_D2x_a-2.0E0*3*I_NAI_F2xz_G4x_b+3*3*I_NAI_F2xz_D2x;
  abcd[903] = 4.0E0*I_NAI_H3x2y_G4x_ab-2.0E0*3*I_NAI_H3x2y_D2x_a-2.0E0*2*I_NAI_Fx2y_G4x_b+2*3*I_NAI_Fx2y_D2x;
  abcd[904] = 4.0E0*I_NAI_H3xyz_G4x_ab-2.0E0*3*I_NAI_H3xyz_D2x_a-2.0E0*2*I_NAI_Fxyz_G4x_b+2*3*I_NAI_Fxyz_D2x;
  abcd[905] = 4.0E0*I_NAI_H3x2z_G4x_ab-2.0E0*3*I_NAI_H3x2z_D2x_a-2.0E0*2*I_NAI_Fx2z_G4x_b+2*3*I_NAI_Fx2z_D2x;
  abcd[906] = 4.0E0*I_NAI_H2x3y_G4x_ab-2.0E0*3*I_NAI_H2x3y_D2x_a-2.0E0*1*I_NAI_F3y_G4x_b+3*I_NAI_F3y_D2x;
  abcd[907] = 4.0E0*I_NAI_H2x2yz_G4x_ab-2.0E0*3*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2yz_G4x_b+3*I_NAI_F2yz_D2x;
  abcd[908] = 4.0E0*I_NAI_H2xy2z_G4x_ab-2.0E0*3*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_Fy2z_G4x_b+3*I_NAI_Fy2z_D2x;
  abcd[909] = 4.0E0*I_NAI_H2x3z_G4x_ab-2.0E0*3*I_NAI_H2x3z_D2x_a-2.0E0*1*I_NAI_F3z_G4x_b+3*I_NAI_F3z_D2x;
  abcd[910] = 4.0E0*I_NAI_Hx4y_G4x_ab-2.0E0*3*I_NAI_Hx4y_D2x_a;
  abcd[911] = 4.0E0*I_NAI_Hx3yz_G4x_ab-2.0E0*3*I_NAI_Hx3yz_D2x_a;
  abcd[912] = 4.0E0*I_NAI_Hx2y2z_G4x_ab-2.0E0*3*I_NAI_Hx2y2z_D2x_a;
  abcd[913] = 4.0E0*I_NAI_Hxy3z_G4x_ab-2.0E0*3*I_NAI_Hxy3z_D2x_a;
  abcd[914] = 4.0E0*I_NAI_Hx4z_G4x_ab-2.0E0*3*I_NAI_Hx4z_D2x_a;
  abcd[915] = 4.0E0*I_NAI_H5x_G3xy_ab-2.0E0*2*I_NAI_H5x_Dxy_a-2.0E0*4*I_NAI_F3x_G3xy_b+4*2*I_NAI_F3x_Dxy;
  abcd[916] = 4.0E0*I_NAI_H4xy_G3xy_ab-2.0E0*2*I_NAI_H4xy_Dxy_a-2.0E0*3*I_NAI_F2xy_G3xy_b+3*2*I_NAI_F2xy_Dxy;
  abcd[917] = 4.0E0*I_NAI_H4xz_G3xy_ab-2.0E0*2*I_NAI_H4xz_Dxy_a-2.0E0*3*I_NAI_F2xz_G3xy_b+3*2*I_NAI_F2xz_Dxy;
  abcd[918] = 4.0E0*I_NAI_H3x2y_G3xy_ab-2.0E0*2*I_NAI_H3x2y_Dxy_a-2.0E0*2*I_NAI_Fx2y_G3xy_b+2*2*I_NAI_Fx2y_Dxy;
  abcd[919] = 4.0E0*I_NAI_H3xyz_G3xy_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a-2.0E0*2*I_NAI_Fxyz_G3xy_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[920] = 4.0E0*I_NAI_H3x2z_G3xy_ab-2.0E0*2*I_NAI_H3x2z_Dxy_a-2.0E0*2*I_NAI_Fx2z_G3xy_b+2*2*I_NAI_Fx2z_Dxy;
  abcd[921] = 4.0E0*I_NAI_H2x3y_G3xy_ab-2.0E0*2*I_NAI_H2x3y_Dxy_a-2.0E0*1*I_NAI_F3y_G3xy_b+2*I_NAI_F3y_Dxy;
  abcd[922] = 4.0E0*I_NAI_H2x2yz_G3xy_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2yz_G3xy_b+2*I_NAI_F2yz_Dxy;
  abcd[923] = 4.0E0*I_NAI_H2xy2z_G3xy_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_Fy2z_G3xy_b+2*I_NAI_Fy2z_Dxy;
  abcd[924] = 4.0E0*I_NAI_H2x3z_G3xy_ab-2.0E0*2*I_NAI_H2x3z_Dxy_a-2.0E0*1*I_NAI_F3z_G3xy_b+2*I_NAI_F3z_Dxy;
  abcd[925] = 4.0E0*I_NAI_Hx4y_G3xy_ab-2.0E0*2*I_NAI_Hx4y_Dxy_a;
  abcd[926] = 4.0E0*I_NAI_Hx3yz_G3xy_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a;
  abcd[927] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a;
  abcd[928] = 4.0E0*I_NAI_Hxy3z_G3xy_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a;
  abcd[929] = 4.0E0*I_NAI_Hx4z_G3xy_ab-2.0E0*2*I_NAI_Hx4z_Dxy_a;
  abcd[930] = 4.0E0*I_NAI_H5x_G3xz_ab-2.0E0*2*I_NAI_H5x_Dxz_a-2.0E0*4*I_NAI_F3x_G3xz_b+4*2*I_NAI_F3x_Dxz;
  abcd[931] = 4.0E0*I_NAI_H4xy_G3xz_ab-2.0E0*2*I_NAI_H4xy_Dxz_a-2.0E0*3*I_NAI_F2xy_G3xz_b+3*2*I_NAI_F2xy_Dxz;
  abcd[932] = 4.0E0*I_NAI_H4xz_G3xz_ab-2.0E0*2*I_NAI_H4xz_Dxz_a-2.0E0*3*I_NAI_F2xz_G3xz_b+3*2*I_NAI_F2xz_Dxz;
  abcd[933] = 4.0E0*I_NAI_H3x2y_G3xz_ab-2.0E0*2*I_NAI_H3x2y_Dxz_a-2.0E0*2*I_NAI_Fx2y_G3xz_b+2*2*I_NAI_Fx2y_Dxz;
  abcd[934] = 4.0E0*I_NAI_H3xyz_G3xz_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a-2.0E0*2*I_NAI_Fxyz_G3xz_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[935] = 4.0E0*I_NAI_H3x2z_G3xz_ab-2.0E0*2*I_NAI_H3x2z_Dxz_a-2.0E0*2*I_NAI_Fx2z_G3xz_b+2*2*I_NAI_Fx2z_Dxz;
  abcd[936] = 4.0E0*I_NAI_H2x3y_G3xz_ab-2.0E0*2*I_NAI_H2x3y_Dxz_a-2.0E0*1*I_NAI_F3y_G3xz_b+2*I_NAI_F3y_Dxz;
  abcd[937] = 4.0E0*I_NAI_H2x2yz_G3xz_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2yz_G3xz_b+2*I_NAI_F2yz_Dxz;
  abcd[938] = 4.0E0*I_NAI_H2xy2z_G3xz_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_Fy2z_G3xz_b+2*I_NAI_Fy2z_Dxz;
  abcd[939] = 4.0E0*I_NAI_H2x3z_G3xz_ab-2.0E0*2*I_NAI_H2x3z_Dxz_a-2.0E0*1*I_NAI_F3z_G3xz_b+2*I_NAI_F3z_Dxz;
  abcd[940] = 4.0E0*I_NAI_Hx4y_G3xz_ab-2.0E0*2*I_NAI_Hx4y_Dxz_a;
  abcd[941] = 4.0E0*I_NAI_Hx3yz_G3xz_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a;
  abcd[942] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a;
  abcd[943] = 4.0E0*I_NAI_Hxy3z_G3xz_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a;
  abcd[944] = 4.0E0*I_NAI_Hx4z_G3xz_ab-2.0E0*2*I_NAI_Hx4z_Dxz_a;
  abcd[945] = 4.0E0*I_NAI_H5x_G2x2y_ab-2.0E0*1*I_NAI_H5x_D2y_a-2.0E0*4*I_NAI_F3x_G2x2y_b+4*1*I_NAI_F3x_D2y;
  abcd[946] = 4.0E0*I_NAI_H4xy_G2x2y_ab-2.0E0*1*I_NAI_H4xy_D2y_a-2.0E0*3*I_NAI_F2xy_G2x2y_b+3*1*I_NAI_F2xy_D2y;
  abcd[947] = 4.0E0*I_NAI_H4xz_G2x2y_ab-2.0E0*1*I_NAI_H4xz_D2y_a-2.0E0*3*I_NAI_F2xz_G2x2y_b+3*1*I_NAI_F2xz_D2y;
  abcd[948] = 4.0E0*I_NAI_H3x2y_G2x2y_ab-2.0E0*1*I_NAI_H3x2y_D2y_a-2.0E0*2*I_NAI_Fx2y_G2x2y_b+2*1*I_NAI_Fx2y_D2y;
  abcd[949] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2y_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2y;
  abcd[950] = 4.0E0*I_NAI_H3x2z_G2x2y_ab-2.0E0*1*I_NAI_H3x2z_D2y_a-2.0E0*2*I_NAI_Fx2z_G2x2y_b+2*1*I_NAI_Fx2z_D2y;
  abcd[951] = 4.0E0*I_NAI_H2x3y_G2x2y_ab-2.0E0*1*I_NAI_H2x3y_D2y_a-2.0E0*1*I_NAI_F3y_G2x2y_b+1*I_NAI_F3y_D2y;
  abcd[952] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2yz_G2x2y_b+1*I_NAI_F2yz_D2y;
  abcd[953] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_Fy2z_G2x2y_b+1*I_NAI_Fy2z_D2y;
  abcd[954] = 4.0E0*I_NAI_H2x3z_G2x2y_ab-2.0E0*1*I_NAI_H2x3z_D2y_a-2.0E0*1*I_NAI_F3z_G2x2y_b+1*I_NAI_F3z_D2y;
  abcd[955] = 4.0E0*I_NAI_Hx4y_G2x2y_ab-2.0E0*1*I_NAI_Hx4y_D2y_a;
  abcd[956] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a;
  abcd[957] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a;
  abcd[958] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a;
  abcd[959] = 4.0E0*I_NAI_Hx4z_G2x2y_ab-2.0E0*1*I_NAI_Hx4z_D2y_a;
  abcd[960] = 4.0E0*I_NAI_H5x_G2xyz_ab-2.0E0*1*I_NAI_H5x_Dyz_a-2.0E0*4*I_NAI_F3x_G2xyz_b+4*1*I_NAI_F3x_Dyz;
  abcd[961] = 4.0E0*I_NAI_H4xy_G2xyz_ab-2.0E0*1*I_NAI_H4xy_Dyz_a-2.0E0*3*I_NAI_F2xy_G2xyz_b+3*1*I_NAI_F2xy_Dyz;
  abcd[962] = 4.0E0*I_NAI_H4xz_G2xyz_ab-2.0E0*1*I_NAI_H4xz_Dyz_a-2.0E0*3*I_NAI_F2xz_G2xyz_b+3*1*I_NAI_F2xz_Dyz;
  abcd[963] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*1*I_NAI_H3x2y_Dyz_a-2.0E0*2*I_NAI_Fx2y_G2xyz_b+2*1*I_NAI_Fx2y_Dyz;
  abcd[964] = 4.0E0*I_NAI_H3xyz_G2xyz_ab-2.0E0*1*I_NAI_H3xyz_Dyz_a-2.0E0*2*I_NAI_Fxyz_G2xyz_b+2*1*I_NAI_Fxyz_Dyz;
  abcd[965] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*1*I_NAI_H3x2z_Dyz_a-2.0E0*2*I_NAI_Fx2z_G2xyz_b+2*1*I_NAI_Fx2z_Dyz;
  abcd[966] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*1*I_NAI_H2x3y_Dyz_a-2.0E0*1*I_NAI_F3y_G2xyz_b+1*I_NAI_F3y_Dyz;
  abcd[967] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2yz_G2xyz_b+1*I_NAI_F2yz_Dyz;
  abcd[968] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_Fy2z_G2xyz_b+1*I_NAI_Fy2z_Dyz;
  abcd[969] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*1*I_NAI_H2x3z_Dyz_a-2.0E0*1*I_NAI_F3z_G2xyz_b+1*I_NAI_F3z_Dyz;
  abcd[970] = 4.0E0*I_NAI_Hx4y_G2xyz_ab-2.0E0*1*I_NAI_Hx4y_Dyz_a;
  abcd[971] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab-2.0E0*1*I_NAI_Hx3yz_Dyz_a;
  abcd[972] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Hx2y2z_Dyz_a;
  abcd[973] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab-2.0E0*1*I_NAI_Hxy3z_Dyz_a;
  abcd[974] = 4.0E0*I_NAI_Hx4z_G2xyz_ab-2.0E0*1*I_NAI_Hx4z_Dyz_a;
  abcd[975] = 4.0E0*I_NAI_H5x_G2x2z_ab-2.0E0*1*I_NAI_H5x_D2z_a-2.0E0*4*I_NAI_F3x_G2x2z_b+4*1*I_NAI_F3x_D2z;
  abcd[976] = 4.0E0*I_NAI_H4xy_G2x2z_ab-2.0E0*1*I_NAI_H4xy_D2z_a-2.0E0*3*I_NAI_F2xy_G2x2z_b+3*1*I_NAI_F2xy_D2z;
  abcd[977] = 4.0E0*I_NAI_H4xz_G2x2z_ab-2.0E0*1*I_NAI_H4xz_D2z_a-2.0E0*3*I_NAI_F2xz_G2x2z_b+3*1*I_NAI_F2xz_D2z;
  abcd[978] = 4.0E0*I_NAI_H3x2y_G2x2z_ab-2.0E0*1*I_NAI_H3x2y_D2z_a-2.0E0*2*I_NAI_Fx2y_G2x2z_b+2*1*I_NAI_Fx2y_D2z;
  abcd[979] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[980] = 4.0E0*I_NAI_H3x2z_G2x2z_ab-2.0E0*1*I_NAI_H3x2z_D2z_a-2.0E0*2*I_NAI_Fx2z_G2x2z_b+2*1*I_NAI_Fx2z_D2z;
  abcd[981] = 4.0E0*I_NAI_H2x3y_G2x2z_ab-2.0E0*1*I_NAI_H2x3y_D2z_a-2.0E0*1*I_NAI_F3y_G2x2z_b+1*I_NAI_F3y_D2z;
  abcd[982] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2yz_G2x2z_b+1*I_NAI_F2yz_D2z;
  abcd[983] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_Fy2z_G2x2z_b+1*I_NAI_Fy2z_D2z;
  abcd[984] = 4.0E0*I_NAI_H2x3z_G2x2z_ab-2.0E0*1*I_NAI_H2x3z_D2z_a-2.0E0*1*I_NAI_F3z_G2x2z_b+1*I_NAI_F3z_D2z;
  abcd[985] = 4.0E0*I_NAI_Hx4y_G2x2z_ab-2.0E0*1*I_NAI_Hx4y_D2z_a;
  abcd[986] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a;
  abcd[987] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a;
  abcd[988] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a;
  abcd[989] = 4.0E0*I_NAI_Hx4z_G2x2z_ab-2.0E0*1*I_NAI_Hx4z_D2z_a;
  abcd[990] = 4.0E0*I_NAI_H5x_Gx3y_ab-2.0E0*4*I_NAI_F3x_Gx3y_b;
  abcd[991] = 4.0E0*I_NAI_H4xy_Gx3y_ab-2.0E0*3*I_NAI_F2xy_Gx3y_b;
  abcd[992] = 4.0E0*I_NAI_H4xz_Gx3y_ab-2.0E0*3*I_NAI_F2xz_Gx3y_b;
  abcd[993] = 4.0E0*I_NAI_H3x2y_Gx3y_ab-2.0E0*2*I_NAI_Fx2y_Gx3y_b;
  abcd[994] = 4.0E0*I_NAI_H3xyz_Gx3y_ab-2.0E0*2*I_NAI_Fxyz_Gx3y_b;
  abcd[995] = 4.0E0*I_NAI_H3x2z_Gx3y_ab-2.0E0*2*I_NAI_Fx2z_Gx3y_b;
  abcd[996] = 4.0E0*I_NAI_H2x3y_Gx3y_ab-2.0E0*1*I_NAI_F3y_Gx3y_b;
  abcd[997] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab-2.0E0*1*I_NAI_F2yz_Gx3y_b;
  abcd[998] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab-2.0E0*1*I_NAI_Fy2z_Gx3y_b;
  abcd[999] = 4.0E0*I_NAI_H2x3z_Gx3y_ab-2.0E0*1*I_NAI_F3z_Gx3y_b;
  abcd[1000] = 4.0E0*I_NAI_Hx4y_Gx3y_ab;
  abcd[1001] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab;
  abcd[1002] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab;
  abcd[1003] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab;
  abcd[1004] = 4.0E0*I_NAI_Hx4z_Gx3y_ab;
  abcd[1005] = 4.0E0*I_NAI_H5x_Gx2yz_ab-2.0E0*4*I_NAI_F3x_Gx2yz_b;
  abcd[1006] = 4.0E0*I_NAI_H4xy_Gx2yz_ab-2.0E0*3*I_NAI_F2xy_Gx2yz_b;
  abcd[1007] = 4.0E0*I_NAI_H4xz_Gx2yz_ab-2.0E0*3*I_NAI_F2xz_Gx2yz_b;
  abcd[1008] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*2*I_NAI_Fx2y_Gx2yz_b;
  abcd[1009] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[1010] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*2*I_NAI_Fx2z_Gx2yz_b;
  abcd[1011] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*1*I_NAI_F3y_Gx2yz_b;
  abcd[1012] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_F2yz_Gx2yz_b;
  abcd[1013] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_Fy2z_Gx2yz_b;
  abcd[1014] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*1*I_NAI_F3z_Gx2yz_b;
  abcd[1015] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab;
  abcd[1016] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab;
  abcd[1017] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab;
  abcd[1018] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab;
  abcd[1019] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab;
  abcd[1020] = 4.0E0*I_NAI_H5x_Gxy2z_ab-2.0E0*4*I_NAI_F3x_Gxy2z_b;
  abcd[1021] = 4.0E0*I_NAI_H4xy_Gxy2z_ab-2.0E0*3*I_NAI_F2xy_Gxy2z_b;
  abcd[1022] = 4.0E0*I_NAI_H4xz_Gxy2z_ab-2.0E0*3*I_NAI_F2xz_Gxy2z_b;
  abcd[1023] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*2*I_NAI_Fx2y_Gxy2z_b;
  abcd[1024] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[1025] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*2*I_NAI_Fx2z_Gxy2z_b;
  abcd[1026] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*1*I_NAI_F3y_Gxy2z_b;
  abcd[1027] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_F2yz_Gxy2z_b;
  abcd[1028] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_Fy2z_Gxy2z_b;
  abcd[1029] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*1*I_NAI_F3z_Gxy2z_b;
  abcd[1030] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab;
  abcd[1031] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab;
  abcd[1032] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab;
  abcd[1033] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab;
  abcd[1034] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab;
  abcd[1035] = 4.0E0*I_NAI_H5x_Gx3z_ab-2.0E0*4*I_NAI_F3x_Gx3z_b;
  abcd[1036] = 4.0E0*I_NAI_H4xy_Gx3z_ab-2.0E0*3*I_NAI_F2xy_Gx3z_b;
  abcd[1037] = 4.0E0*I_NAI_H4xz_Gx3z_ab-2.0E0*3*I_NAI_F2xz_Gx3z_b;
  abcd[1038] = 4.0E0*I_NAI_H3x2y_Gx3z_ab-2.0E0*2*I_NAI_Fx2y_Gx3z_b;
  abcd[1039] = 4.0E0*I_NAI_H3xyz_Gx3z_ab-2.0E0*2*I_NAI_Fxyz_Gx3z_b;
  abcd[1040] = 4.0E0*I_NAI_H3x2z_Gx3z_ab-2.0E0*2*I_NAI_Fx2z_Gx3z_b;
  abcd[1041] = 4.0E0*I_NAI_H2x3y_Gx3z_ab-2.0E0*1*I_NAI_F3y_Gx3z_b;
  abcd[1042] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab-2.0E0*1*I_NAI_F2yz_Gx3z_b;
  abcd[1043] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab-2.0E0*1*I_NAI_Fy2z_Gx3z_b;
  abcd[1044] = 4.0E0*I_NAI_H2x3z_Gx3z_ab-2.0E0*1*I_NAI_F3z_Gx3z_b;
  abcd[1045] = 4.0E0*I_NAI_Hx4y_Gx3z_ab;
  abcd[1046] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab;
  abcd[1047] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab;
  abcd[1048] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab;
  abcd[1049] = 4.0E0*I_NAI_Hx4z_Gx3z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1050] = 4.0E0*I_NAI_H5x_G3xy_ab-2.0E0*4*I_NAI_F3x_G3xy_b;
  abcd[1051] = 4.0E0*I_NAI_H4xy_G3xy_ab-2.0E0*3*I_NAI_F2xy_G3xy_b;
  abcd[1052] = 4.0E0*I_NAI_H4xz_G3xy_ab-2.0E0*3*I_NAI_F2xz_G3xy_b;
  abcd[1053] = 4.0E0*I_NAI_H3x2y_G3xy_ab-2.0E0*2*I_NAI_Fx2y_G3xy_b;
  abcd[1054] = 4.0E0*I_NAI_H3xyz_G3xy_ab-2.0E0*2*I_NAI_Fxyz_G3xy_b;
  abcd[1055] = 4.0E0*I_NAI_H3x2z_G3xy_ab-2.0E0*2*I_NAI_Fx2z_G3xy_b;
  abcd[1056] = 4.0E0*I_NAI_H2x3y_G3xy_ab-2.0E0*1*I_NAI_F3y_G3xy_b;
  abcd[1057] = 4.0E0*I_NAI_H2x2yz_G3xy_ab-2.0E0*1*I_NAI_F2yz_G3xy_b;
  abcd[1058] = 4.0E0*I_NAI_H2xy2z_G3xy_ab-2.0E0*1*I_NAI_Fy2z_G3xy_b;
  abcd[1059] = 4.0E0*I_NAI_H2x3z_G3xy_ab-2.0E0*1*I_NAI_F3z_G3xy_b;
  abcd[1060] = 4.0E0*I_NAI_Hx4y_G3xy_ab;
  abcd[1061] = 4.0E0*I_NAI_Hx3yz_G3xy_ab;
  abcd[1062] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab;
  abcd[1063] = 4.0E0*I_NAI_Hxy3z_G3xy_ab;
  abcd[1064] = 4.0E0*I_NAI_Hx4z_G3xy_ab;
  abcd[1065] = 4.0E0*I_NAI_H5x_G2x2y_ab-2.0E0*1*I_NAI_H5x_D2x_a-2.0E0*4*I_NAI_F3x_G2x2y_b+4*1*I_NAI_F3x_D2x;
  abcd[1066] = 4.0E0*I_NAI_H4xy_G2x2y_ab-2.0E0*1*I_NAI_H4xy_D2x_a-2.0E0*3*I_NAI_F2xy_G2x2y_b+3*1*I_NAI_F2xy_D2x;
  abcd[1067] = 4.0E0*I_NAI_H4xz_G2x2y_ab-2.0E0*1*I_NAI_H4xz_D2x_a-2.0E0*3*I_NAI_F2xz_G2x2y_b+3*1*I_NAI_F2xz_D2x;
  abcd[1068] = 4.0E0*I_NAI_H3x2y_G2x2y_ab-2.0E0*1*I_NAI_H3x2y_D2x_a-2.0E0*2*I_NAI_Fx2y_G2x2y_b+2*1*I_NAI_Fx2y_D2x;
  abcd[1069] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2x;
  abcd[1070] = 4.0E0*I_NAI_H3x2z_G2x2y_ab-2.0E0*1*I_NAI_H3x2z_D2x_a-2.0E0*2*I_NAI_Fx2z_G2x2y_b+2*1*I_NAI_Fx2z_D2x;
  abcd[1071] = 4.0E0*I_NAI_H2x3y_G2x2y_ab-2.0E0*1*I_NAI_H2x3y_D2x_a-2.0E0*1*I_NAI_F3y_G2x2y_b+1*I_NAI_F3y_D2x;
  abcd[1072] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2yz_G2x2y_b+1*I_NAI_F2yz_D2x;
  abcd[1073] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_Fy2z_G2x2y_b+1*I_NAI_Fy2z_D2x;
  abcd[1074] = 4.0E0*I_NAI_H2x3z_G2x2y_ab-2.0E0*1*I_NAI_H2x3z_D2x_a-2.0E0*1*I_NAI_F3z_G2x2y_b+1*I_NAI_F3z_D2x;
  abcd[1075] = 4.0E0*I_NAI_Hx4y_G2x2y_ab-2.0E0*1*I_NAI_Hx4y_D2x_a;
  abcd[1076] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a;
  abcd[1077] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a;
  abcd[1078] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a;
  abcd[1079] = 4.0E0*I_NAI_Hx4z_G2x2y_ab-2.0E0*1*I_NAI_Hx4z_D2x_a;
  abcd[1080] = 4.0E0*I_NAI_H5x_G2xyz_ab-2.0E0*4*I_NAI_F3x_G2xyz_b;
  abcd[1081] = 4.0E0*I_NAI_H4xy_G2xyz_ab-2.0E0*3*I_NAI_F2xy_G2xyz_b;
  abcd[1082] = 4.0E0*I_NAI_H4xz_G2xyz_ab-2.0E0*3*I_NAI_F2xz_G2xyz_b;
  abcd[1083] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*2*I_NAI_Fx2y_G2xyz_b;
  abcd[1084] = 4.0E0*I_NAI_H3xyz_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[1085] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*2*I_NAI_Fx2z_G2xyz_b;
  abcd[1086] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*1*I_NAI_F3y_G2xyz_b;
  abcd[1087] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_F2yz_G2xyz_b;
  abcd[1088] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_Fy2z_G2xyz_b;
  abcd[1089] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*1*I_NAI_F3z_G2xyz_b;
  abcd[1090] = 4.0E0*I_NAI_Hx4y_G2xyz_ab;
  abcd[1091] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab;
  abcd[1092] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab;
  abcd[1093] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab;
  abcd[1094] = 4.0E0*I_NAI_Hx4z_G2xyz_ab;
  abcd[1095] = 4.0E0*I_NAI_H5x_Gx3y_ab-2.0E0*2*I_NAI_H5x_Dxy_a-2.0E0*4*I_NAI_F3x_Gx3y_b+4*2*I_NAI_F3x_Dxy;
  abcd[1096] = 4.0E0*I_NAI_H4xy_Gx3y_ab-2.0E0*2*I_NAI_H4xy_Dxy_a-2.0E0*3*I_NAI_F2xy_Gx3y_b+3*2*I_NAI_F2xy_Dxy;
  abcd[1097] = 4.0E0*I_NAI_H4xz_Gx3y_ab-2.0E0*2*I_NAI_H4xz_Dxy_a-2.0E0*3*I_NAI_F2xz_Gx3y_b+3*2*I_NAI_F2xz_Dxy;
  abcd[1098] = 4.0E0*I_NAI_H3x2y_Gx3y_ab-2.0E0*2*I_NAI_H3x2y_Dxy_a-2.0E0*2*I_NAI_Fx2y_Gx3y_b+2*2*I_NAI_Fx2y_Dxy;
  abcd[1099] = 4.0E0*I_NAI_H3xyz_Gx3y_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gx3y_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[1100] = 4.0E0*I_NAI_H3x2z_Gx3y_ab-2.0E0*2*I_NAI_H3x2z_Dxy_a-2.0E0*2*I_NAI_Fx2z_Gx3y_b+2*2*I_NAI_Fx2z_Dxy;
  abcd[1101] = 4.0E0*I_NAI_H2x3y_Gx3y_ab-2.0E0*2*I_NAI_H2x3y_Dxy_a-2.0E0*1*I_NAI_F3y_Gx3y_b+2*I_NAI_F3y_Dxy;
  abcd[1102] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2yz_Gx3y_b+2*I_NAI_F2yz_Dxy;
  abcd[1103] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_Fy2z_Gx3y_b+2*I_NAI_Fy2z_Dxy;
  abcd[1104] = 4.0E0*I_NAI_H2x3z_Gx3y_ab-2.0E0*2*I_NAI_H2x3z_Dxy_a-2.0E0*1*I_NAI_F3z_Gx3y_b+2*I_NAI_F3z_Dxy;
  abcd[1105] = 4.0E0*I_NAI_Hx4y_Gx3y_ab-2.0E0*2*I_NAI_Hx4y_Dxy_a;
  abcd[1106] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a;
  abcd[1107] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a;
  abcd[1108] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a;
  abcd[1109] = 4.0E0*I_NAI_Hx4z_Gx3y_ab-2.0E0*2*I_NAI_Hx4z_Dxy_a;
  abcd[1110] = 4.0E0*I_NAI_H5x_Gx2yz_ab-2.0E0*1*I_NAI_H5x_Dxz_a-2.0E0*4*I_NAI_F3x_Gx2yz_b+4*1*I_NAI_F3x_Dxz;
  abcd[1111] = 4.0E0*I_NAI_H4xy_Gx2yz_ab-2.0E0*1*I_NAI_H4xy_Dxz_a-2.0E0*3*I_NAI_F2xy_Gx2yz_b+3*1*I_NAI_F2xy_Dxz;
  abcd[1112] = 4.0E0*I_NAI_H4xz_Gx2yz_ab-2.0E0*1*I_NAI_H4xz_Dxz_a-2.0E0*3*I_NAI_F2xz_Gx2yz_b+3*1*I_NAI_F2xz_Dxz;
  abcd[1113] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*1*I_NAI_H3x2y_Dxz_a-2.0E0*2*I_NAI_Fx2y_Gx2yz_b+2*1*I_NAI_Fx2y_Dxz;
  abcd[1114] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab-2.0E0*1*I_NAI_H3xyz_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx2yz_b+2*1*I_NAI_Fxyz_Dxz;
  abcd[1115] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*1*I_NAI_H3x2z_Dxz_a-2.0E0*2*I_NAI_Fx2z_Gx2yz_b+2*1*I_NAI_Fx2z_Dxz;
  abcd[1116] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*1*I_NAI_H2x3y_Dxz_a-2.0E0*1*I_NAI_F3y_Gx2yz_b+1*I_NAI_F3y_Dxz;
  abcd[1117] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2yz_Gx2yz_b+1*I_NAI_F2yz_Dxz;
  abcd[1118] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_Fy2z_Gx2yz_b+1*I_NAI_Fy2z_Dxz;
  abcd[1119] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*1*I_NAI_H2x3z_Dxz_a-2.0E0*1*I_NAI_F3z_Gx2yz_b+1*I_NAI_F3z_Dxz;
  abcd[1120] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab-2.0E0*1*I_NAI_Hx4y_Dxz_a;
  abcd[1121] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab-2.0E0*1*I_NAI_Hx3yz_Dxz_a;
  abcd[1122] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Hx2y2z_Dxz_a;
  abcd[1123] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab-2.0E0*1*I_NAI_Hxy3z_Dxz_a;
  abcd[1124] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab-2.0E0*1*I_NAI_Hx4z_Dxz_a;
  abcd[1125] = 4.0E0*I_NAI_H5x_Gxy2z_ab-2.0E0*4*I_NAI_F3x_Gxy2z_b;
  abcd[1126] = 4.0E0*I_NAI_H4xy_Gxy2z_ab-2.0E0*3*I_NAI_F2xy_Gxy2z_b;
  abcd[1127] = 4.0E0*I_NAI_H4xz_Gxy2z_ab-2.0E0*3*I_NAI_F2xz_Gxy2z_b;
  abcd[1128] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*2*I_NAI_Fx2y_Gxy2z_b;
  abcd[1129] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[1130] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*2*I_NAI_Fx2z_Gxy2z_b;
  abcd[1131] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*1*I_NAI_F3y_Gxy2z_b;
  abcd[1132] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_F2yz_Gxy2z_b;
  abcd[1133] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_Fy2z_Gxy2z_b;
  abcd[1134] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*1*I_NAI_F3z_Gxy2z_b;
  abcd[1135] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab;
  abcd[1136] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab;
  abcd[1137] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab;
  abcd[1138] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab;
  abcd[1139] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab;
  abcd[1140] = 4.0E0*I_NAI_H5x_G4y_ab-2.0E0*3*I_NAI_H5x_D2y_a-2.0E0*4*I_NAI_F3x_G4y_b+4*3*I_NAI_F3x_D2y;
  abcd[1141] = 4.0E0*I_NAI_H4xy_G4y_ab-2.0E0*3*I_NAI_H4xy_D2y_a-2.0E0*3*I_NAI_F2xy_G4y_b+3*3*I_NAI_F2xy_D2y;
  abcd[1142] = 4.0E0*I_NAI_H4xz_G4y_ab-2.0E0*3*I_NAI_H4xz_D2y_a-2.0E0*3*I_NAI_F2xz_G4y_b+3*3*I_NAI_F2xz_D2y;
  abcd[1143] = 4.0E0*I_NAI_H3x2y_G4y_ab-2.0E0*3*I_NAI_H3x2y_D2y_a-2.0E0*2*I_NAI_Fx2y_G4y_b+2*3*I_NAI_Fx2y_D2y;
  abcd[1144] = 4.0E0*I_NAI_H3xyz_G4y_ab-2.0E0*3*I_NAI_H3xyz_D2y_a-2.0E0*2*I_NAI_Fxyz_G4y_b+2*3*I_NAI_Fxyz_D2y;
  abcd[1145] = 4.0E0*I_NAI_H3x2z_G4y_ab-2.0E0*3*I_NAI_H3x2z_D2y_a-2.0E0*2*I_NAI_Fx2z_G4y_b+2*3*I_NAI_Fx2z_D2y;
  abcd[1146] = 4.0E0*I_NAI_H2x3y_G4y_ab-2.0E0*3*I_NAI_H2x3y_D2y_a-2.0E0*1*I_NAI_F3y_G4y_b+3*I_NAI_F3y_D2y;
  abcd[1147] = 4.0E0*I_NAI_H2x2yz_G4y_ab-2.0E0*3*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2yz_G4y_b+3*I_NAI_F2yz_D2y;
  abcd[1148] = 4.0E0*I_NAI_H2xy2z_G4y_ab-2.0E0*3*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_Fy2z_G4y_b+3*I_NAI_Fy2z_D2y;
  abcd[1149] = 4.0E0*I_NAI_H2x3z_G4y_ab-2.0E0*3*I_NAI_H2x3z_D2y_a-2.0E0*1*I_NAI_F3z_G4y_b+3*I_NAI_F3z_D2y;
  abcd[1150] = 4.0E0*I_NAI_Hx4y_G4y_ab-2.0E0*3*I_NAI_Hx4y_D2y_a;
  abcd[1151] = 4.0E0*I_NAI_Hx3yz_G4y_ab-2.0E0*3*I_NAI_Hx3yz_D2y_a;
  abcd[1152] = 4.0E0*I_NAI_Hx2y2z_G4y_ab-2.0E0*3*I_NAI_Hx2y2z_D2y_a;
  abcd[1153] = 4.0E0*I_NAI_Hxy3z_G4y_ab-2.0E0*3*I_NAI_Hxy3z_D2y_a;
  abcd[1154] = 4.0E0*I_NAI_Hx4z_G4y_ab-2.0E0*3*I_NAI_Hx4z_D2y_a;
  abcd[1155] = 4.0E0*I_NAI_H5x_G3yz_ab-2.0E0*2*I_NAI_H5x_Dyz_a-2.0E0*4*I_NAI_F3x_G3yz_b+4*2*I_NAI_F3x_Dyz;
  abcd[1156] = 4.0E0*I_NAI_H4xy_G3yz_ab-2.0E0*2*I_NAI_H4xy_Dyz_a-2.0E0*3*I_NAI_F2xy_G3yz_b+3*2*I_NAI_F2xy_Dyz;
  abcd[1157] = 4.0E0*I_NAI_H4xz_G3yz_ab-2.0E0*2*I_NAI_H4xz_Dyz_a-2.0E0*3*I_NAI_F2xz_G3yz_b+3*2*I_NAI_F2xz_Dyz;
  abcd[1158] = 4.0E0*I_NAI_H3x2y_G3yz_ab-2.0E0*2*I_NAI_H3x2y_Dyz_a-2.0E0*2*I_NAI_Fx2y_G3yz_b+2*2*I_NAI_Fx2y_Dyz;
  abcd[1159] = 4.0E0*I_NAI_H3xyz_G3yz_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a-2.0E0*2*I_NAI_Fxyz_G3yz_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[1160] = 4.0E0*I_NAI_H3x2z_G3yz_ab-2.0E0*2*I_NAI_H3x2z_Dyz_a-2.0E0*2*I_NAI_Fx2z_G3yz_b+2*2*I_NAI_Fx2z_Dyz;
  abcd[1161] = 4.0E0*I_NAI_H2x3y_G3yz_ab-2.0E0*2*I_NAI_H2x3y_Dyz_a-2.0E0*1*I_NAI_F3y_G3yz_b+2*I_NAI_F3y_Dyz;
  abcd[1162] = 4.0E0*I_NAI_H2x2yz_G3yz_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2yz_G3yz_b+2*I_NAI_F2yz_Dyz;
  abcd[1163] = 4.0E0*I_NAI_H2xy2z_G3yz_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_Fy2z_G3yz_b+2*I_NAI_Fy2z_Dyz;
  abcd[1164] = 4.0E0*I_NAI_H2x3z_G3yz_ab-2.0E0*2*I_NAI_H2x3z_Dyz_a-2.0E0*1*I_NAI_F3z_G3yz_b+2*I_NAI_F3z_Dyz;
  abcd[1165] = 4.0E0*I_NAI_Hx4y_G3yz_ab-2.0E0*2*I_NAI_Hx4y_Dyz_a;
  abcd[1166] = 4.0E0*I_NAI_Hx3yz_G3yz_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a;
  abcd[1167] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a;
  abcd[1168] = 4.0E0*I_NAI_Hxy3z_G3yz_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a;
  abcd[1169] = 4.0E0*I_NAI_Hx4z_G3yz_ab-2.0E0*2*I_NAI_Hx4z_Dyz_a;
  abcd[1170] = 4.0E0*I_NAI_H5x_G2y2z_ab-2.0E0*1*I_NAI_H5x_D2z_a-2.0E0*4*I_NAI_F3x_G2y2z_b+4*1*I_NAI_F3x_D2z;
  abcd[1171] = 4.0E0*I_NAI_H4xy_G2y2z_ab-2.0E0*1*I_NAI_H4xy_D2z_a-2.0E0*3*I_NAI_F2xy_G2y2z_b+3*1*I_NAI_F2xy_D2z;
  abcd[1172] = 4.0E0*I_NAI_H4xz_G2y2z_ab-2.0E0*1*I_NAI_H4xz_D2z_a-2.0E0*3*I_NAI_F2xz_G2y2z_b+3*1*I_NAI_F2xz_D2z;
  abcd[1173] = 4.0E0*I_NAI_H3x2y_G2y2z_ab-2.0E0*1*I_NAI_H3x2y_D2z_a-2.0E0*2*I_NAI_Fx2y_G2y2z_b+2*1*I_NAI_Fx2y_D2z;
  abcd[1174] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[1175] = 4.0E0*I_NAI_H3x2z_G2y2z_ab-2.0E0*1*I_NAI_H3x2z_D2z_a-2.0E0*2*I_NAI_Fx2z_G2y2z_b+2*1*I_NAI_Fx2z_D2z;
  abcd[1176] = 4.0E0*I_NAI_H2x3y_G2y2z_ab-2.0E0*1*I_NAI_H2x3y_D2z_a-2.0E0*1*I_NAI_F3y_G2y2z_b+1*I_NAI_F3y_D2z;
  abcd[1177] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2yz_G2y2z_b+1*I_NAI_F2yz_D2z;
  abcd[1178] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_Fy2z_G2y2z_b+1*I_NAI_Fy2z_D2z;
  abcd[1179] = 4.0E0*I_NAI_H2x3z_G2y2z_ab-2.0E0*1*I_NAI_H2x3z_D2z_a-2.0E0*1*I_NAI_F3z_G2y2z_b+1*I_NAI_F3z_D2z;
  abcd[1180] = 4.0E0*I_NAI_Hx4y_G2y2z_ab-2.0E0*1*I_NAI_Hx4y_D2z_a;
  abcd[1181] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a;
  abcd[1182] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a;
  abcd[1183] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a;
  abcd[1184] = 4.0E0*I_NAI_Hx4z_G2y2z_ab-2.0E0*1*I_NAI_Hx4z_D2z_a;
  abcd[1185] = 4.0E0*I_NAI_H5x_Gy3z_ab-2.0E0*4*I_NAI_F3x_Gy3z_b;
  abcd[1186] = 4.0E0*I_NAI_H4xy_Gy3z_ab-2.0E0*3*I_NAI_F2xy_Gy3z_b;
  abcd[1187] = 4.0E0*I_NAI_H4xz_Gy3z_ab-2.0E0*3*I_NAI_F2xz_Gy3z_b;
  abcd[1188] = 4.0E0*I_NAI_H3x2y_Gy3z_ab-2.0E0*2*I_NAI_Fx2y_Gy3z_b;
  abcd[1189] = 4.0E0*I_NAI_H3xyz_Gy3z_ab-2.0E0*2*I_NAI_Fxyz_Gy3z_b;
  abcd[1190] = 4.0E0*I_NAI_H3x2z_Gy3z_ab-2.0E0*2*I_NAI_Fx2z_Gy3z_b;
  abcd[1191] = 4.0E0*I_NAI_H2x3y_Gy3z_ab-2.0E0*1*I_NAI_F3y_Gy3z_b;
  abcd[1192] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab-2.0E0*1*I_NAI_F2yz_Gy3z_b;
  abcd[1193] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab-2.0E0*1*I_NAI_Fy2z_Gy3z_b;
  abcd[1194] = 4.0E0*I_NAI_H2x3z_Gy3z_ab-2.0E0*1*I_NAI_F3z_Gy3z_b;
  abcd[1195] = 4.0E0*I_NAI_Hx4y_Gy3z_ab;
  abcd[1196] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab;
  abcd[1197] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab;
  abcd[1198] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab;
  abcd[1199] = 4.0E0*I_NAI_Hx4z_Gy3z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1200] = 4.0E0*I_NAI_H5x_G3xz_ab-2.0E0*4*I_NAI_F3x_G3xz_b;
  abcd[1201] = 4.0E0*I_NAI_H4xy_G3xz_ab-2.0E0*3*I_NAI_F2xy_G3xz_b;
  abcd[1202] = 4.0E0*I_NAI_H4xz_G3xz_ab-2.0E0*3*I_NAI_F2xz_G3xz_b;
  abcd[1203] = 4.0E0*I_NAI_H3x2y_G3xz_ab-2.0E0*2*I_NAI_Fx2y_G3xz_b;
  abcd[1204] = 4.0E0*I_NAI_H3xyz_G3xz_ab-2.0E0*2*I_NAI_Fxyz_G3xz_b;
  abcd[1205] = 4.0E0*I_NAI_H3x2z_G3xz_ab-2.0E0*2*I_NAI_Fx2z_G3xz_b;
  abcd[1206] = 4.0E0*I_NAI_H2x3y_G3xz_ab-2.0E0*1*I_NAI_F3y_G3xz_b;
  abcd[1207] = 4.0E0*I_NAI_H2x2yz_G3xz_ab-2.0E0*1*I_NAI_F2yz_G3xz_b;
  abcd[1208] = 4.0E0*I_NAI_H2xy2z_G3xz_ab-2.0E0*1*I_NAI_Fy2z_G3xz_b;
  abcd[1209] = 4.0E0*I_NAI_H2x3z_G3xz_ab-2.0E0*1*I_NAI_F3z_G3xz_b;
  abcd[1210] = 4.0E0*I_NAI_Hx4y_G3xz_ab;
  abcd[1211] = 4.0E0*I_NAI_Hx3yz_G3xz_ab;
  abcd[1212] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab;
  abcd[1213] = 4.0E0*I_NAI_Hxy3z_G3xz_ab;
  abcd[1214] = 4.0E0*I_NAI_Hx4z_G3xz_ab;
  abcd[1215] = 4.0E0*I_NAI_H5x_G2xyz_ab-2.0E0*4*I_NAI_F3x_G2xyz_b;
  abcd[1216] = 4.0E0*I_NAI_H4xy_G2xyz_ab-2.0E0*3*I_NAI_F2xy_G2xyz_b;
  abcd[1217] = 4.0E0*I_NAI_H4xz_G2xyz_ab-2.0E0*3*I_NAI_F2xz_G2xyz_b;
  abcd[1218] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*2*I_NAI_Fx2y_G2xyz_b;
  abcd[1219] = 4.0E0*I_NAI_H3xyz_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[1220] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*2*I_NAI_Fx2z_G2xyz_b;
  abcd[1221] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*1*I_NAI_F3y_G2xyz_b;
  abcd[1222] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_F2yz_G2xyz_b;
  abcd[1223] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_Fy2z_G2xyz_b;
  abcd[1224] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*1*I_NAI_F3z_G2xyz_b;
  abcd[1225] = 4.0E0*I_NAI_Hx4y_G2xyz_ab;
  abcd[1226] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab;
  abcd[1227] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab;
  abcd[1228] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab;
  abcd[1229] = 4.0E0*I_NAI_Hx4z_G2xyz_ab;
  abcd[1230] = 4.0E0*I_NAI_H5x_G2x2z_ab-2.0E0*1*I_NAI_H5x_D2x_a-2.0E0*4*I_NAI_F3x_G2x2z_b+4*1*I_NAI_F3x_D2x;
  abcd[1231] = 4.0E0*I_NAI_H4xy_G2x2z_ab-2.0E0*1*I_NAI_H4xy_D2x_a-2.0E0*3*I_NAI_F2xy_G2x2z_b+3*1*I_NAI_F2xy_D2x;
  abcd[1232] = 4.0E0*I_NAI_H4xz_G2x2z_ab-2.0E0*1*I_NAI_H4xz_D2x_a-2.0E0*3*I_NAI_F2xz_G2x2z_b+3*1*I_NAI_F2xz_D2x;
  abcd[1233] = 4.0E0*I_NAI_H3x2y_G2x2z_ab-2.0E0*1*I_NAI_H3x2y_D2x_a-2.0E0*2*I_NAI_Fx2y_G2x2z_b+2*1*I_NAI_Fx2y_D2x;
  abcd[1234] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2x;
  abcd[1235] = 4.0E0*I_NAI_H3x2z_G2x2z_ab-2.0E0*1*I_NAI_H3x2z_D2x_a-2.0E0*2*I_NAI_Fx2z_G2x2z_b+2*1*I_NAI_Fx2z_D2x;
  abcd[1236] = 4.0E0*I_NAI_H2x3y_G2x2z_ab-2.0E0*1*I_NAI_H2x3y_D2x_a-2.0E0*1*I_NAI_F3y_G2x2z_b+1*I_NAI_F3y_D2x;
  abcd[1237] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2yz_G2x2z_b+1*I_NAI_F2yz_D2x;
  abcd[1238] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_Fy2z_G2x2z_b+1*I_NAI_Fy2z_D2x;
  abcd[1239] = 4.0E0*I_NAI_H2x3z_G2x2z_ab-2.0E0*1*I_NAI_H2x3z_D2x_a-2.0E0*1*I_NAI_F3z_G2x2z_b+1*I_NAI_F3z_D2x;
  abcd[1240] = 4.0E0*I_NAI_Hx4y_G2x2z_ab-2.0E0*1*I_NAI_Hx4y_D2x_a;
  abcd[1241] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a;
  abcd[1242] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a;
  abcd[1243] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a;
  abcd[1244] = 4.0E0*I_NAI_Hx4z_G2x2z_ab-2.0E0*1*I_NAI_Hx4z_D2x_a;
  abcd[1245] = 4.0E0*I_NAI_H5x_Gx2yz_ab-2.0E0*4*I_NAI_F3x_Gx2yz_b;
  abcd[1246] = 4.0E0*I_NAI_H4xy_Gx2yz_ab-2.0E0*3*I_NAI_F2xy_Gx2yz_b;
  abcd[1247] = 4.0E0*I_NAI_H4xz_Gx2yz_ab-2.0E0*3*I_NAI_F2xz_Gx2yz_b;
  abcd[1248] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*2*I_NAI_Fx2y_Gx2yz_b;
  abcd[1249] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[1250] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*2*I_NAI_Fx2z_Gx2yz_b;
  abcd[1251] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*1*I_NAI_F3y_Gx2yz_b;
  abcd[1252] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_F2yz_Gx2yz_b;
  abcd[1253] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_Fy2z_Gx2yz_b;
  abcd[1254] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*1*I_NAI_F3z_Gx2yz_b;
  abcd[1255] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab;
  abcd[1256] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab;
  abcd[1257] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab;
  abcd[1258] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab;
  abcd[1259] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab;
  abcd[1260] = 4.0E0*I_NAI_H5x_Gxy2z_ab-2.0E0*1*I_NAI_H5x_Dxy_a-2.0E0*4*I_NAI_F3x_Gxy2z_b+4*1*I_NAI_F3x_Dxy;
  abcd[1261] = 4.0E0*I_NAI_H4xy_Gxy2z_ab-2.0E0*1*I_NAI_H4xy_Dxy_a-2.0E0*3*I_NAI_F2xy_Gxy2z_b+3*1*I_NAI_F2xy_Dxy;
  abcd[1262] = 4.0E0*I_NAI_H4xz_Gxy2z_ab-2.0E0*1*I_NAI_H4xz_Dxy_a-2.0E0*3*I_NAI_F2xz_Gxy2z_b+3*1*I_NAI_F2xz_Dxy;
  abcd[1263] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*1*I_NAI_H3x2y_Dxy_a-2.0E0*2*I_NAI_Fx2y_Gxy2z_b+2*1*I_NAI_Fx2y_Dxy;
  abcd[1264] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab-2.0E0*1*I_NAI_H3xyz_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gxy2z_b+2*1*I_NAI_Fxyz_Dxy;
  abcd[1265] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*1*I_NAI_H3x2z_Dxy_a-2.0E0*2*I_NAI_Fx2z_Gxy2z_b+2*1*I_NAI_Fx2z_Dxy;
  abcd[1266] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*1*I_NAI_H2x3y_Dxy_a-2.0E0*1*I_NAI_F3y_Gxy2z_b+1*I_NAI_F3y_Dxy;
  abcd[1267] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2yz_Gxy2z_b+1*I_NAI_F2yz_Dxy;
  abcd[1268] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_Fy2z_Gxy2z_b+1*I_NAI_Fy2z_Dxy;
  abcd[1269] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*1*I_NAI_H2x3z_Dxy_a-2.0E0*1*I_NAI_F3z_Gxy2z_b+1*I_NAI_F3z_Dxy;
  abcd[1270] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab-2.0E0*1*I_NAI_Hx4y_Dxy_a;
  abcd[1271] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab-2.0E0*1*I_NAI_Hx3yz_Dxy_a;
  abcd[1272] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Hx2y2z_Dxy_a;
  abcd[1273] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab-2.0E0*1*I_NAI_Hxy3z_Dxy_a;
  abcd[1274] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab-2.0E0*1*I_NAI_Hx4z_Dxy_a;
  abcd[1275] = 4.0E0*I_NAI_H5x_Gx3z_ab-2.0E0*2*I_NAI_H5x_Dxz_a-2.0E0*4*I_NAI_F3x_Gx3z_b+4*2*I_NAI_F3x_Dxz;
  abcd[1276] = 4.0E0*I_NAI_H4xy_Gx3z_ab-2.0E0*2*I_NAI_H4xy_Dxz_a-2.0E0*3*I_NAI_F2xy_Gx3z_b+3*2*I_NAI_F2xy_Dxz;
  abcd[1277] = 4.0E0*I_NAI_H4xz_Gx3z_ab-2.0E0*2*I_NAI_H4xz_Dxz_a-2.0E0*3*I_NAI_F2xz_Gx3z_b+3*2*I_NAI_F2xz_Dxz;
  abcd[1278] = 4.0E0*I_NAI_H3x2y_Gx3z_ab-2.0E0*2*I_NAI_H3x2y_Dxz_a-2.0E0*2*I_NAI_Fx2y_Gx3z_b+2*2*I_NAI_Fx2y_Dxz;
  abcd[1279] = 4.0E0*I_NAI_H3xyz_Gx3z_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx3z_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[1280] = 4.0E0*I_NAI_H3x2z_Gx3z_ab-2.0E0*2*I_NAI_H3x2z_Dxz_a-2.0E0*2*I_NAI_Fx2z_Gx3z_b+2*2*I_NAI_Fx2z_Dxz;
  abcd[1281] = 4.0E0*I_NAI_H2x3y_Gx3z_ab-2.0E0*2*I_NAI_H2x3y_Dxz_a-2.0E0*1*I_NAI_F3y_Gx3z_b+2*I_NAI_F3y_Dxz;
  abcd[1282] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2yz_Gx3z_b+2*I_NAI_F2yz_Dxz;
  abcd[1283] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_Fy2z_Gx3z_b+2*I_NAI_Fy2z_Dxz;
  abcd[1284] = 4.0E0*I_NAI_H2x3z_Gx3z_ab-2.0E0*2*I_NAI_H2x3z_Dxz_a-2.0E0*1*I_NAI_F3z_Gx3z_b+2*I_NAI_F3z_Dxz;
  abcd[1285] = 4.0E0*I_NAI_Hx4y_Gx3z_ab-2.0E0*2*I_NAI_Hx4y_Dxz_a;
  abcd[1286] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a;
  abcd[1287] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a;
  abcd[1288] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a;
  abcd[1289] = 4.0E0*I_NAI_Hx4z_Gx3z_ab-2.0E0*2*I_NAI_Hx4z_Dxz_a;
  abcd[1290] = 4.0E0*I_NAI_H5x_G3yz_ab-2.0E0*4*I_NAI_F3x_G3yz_b;
  abcd[1291] = 4.0E0*I_NAI_H4xy_G3yz_ab-2.0E0*3*I_NAI_F2xy_G3yz_b;
  abcd[1292] = 4.0E0*I_NAI_H4xz_G3yz_ab-2.0E0*3*I_NAI_F2xz_G3yz_b;
  abcd[1293] = 4.0E0*I_NAI_H3x2y_G3yz_ab-2.0E0*2*I_NAI_Fx2y_G3yz_b;
  abcd[1294] = 4.0E0*I_NAI_H3xyz_G3yz_ab-2.0E0*2*I_NAI_Fxyz_G3yz_b;
  abcd[1295] = 4.0E0*I_NAI_H3x2z_G3yz_ab-2.0E0*2*I_NAI_Fx2z_G3yz_b;
  abcd[1296] = 4.0E0*I_NAI_H2x3y_G3yz_ab-2.0E0*1*I_NAI_F3y_G3yz_b;
  abcd[1297] = 4.0E0*I_NAI_H2x2yz_G3yz_ab-2.0E0*1*I_NAI_F2yz_G3yz_b;
  abcd[1298] = 4.0E0*I_NAI_H2xy2z_G3yz_ab-2.0E0*1*I_NAI_Fy2z_G3yz_b;
  abcd[1299] = 4.0E0*I_NAI_H2x3z_G3yz_ab-2.0E0*1*I_NAI_F3z_G3yz_b;
  abcd[1300] = 4.0E0*I_NAI_Hx4y_G3yz_ab;
  abcd[1301] = 4.0E0*I_NAI_Hx3yz_G3yz_ab;
  abcd[1302] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab;
  abcd[1303] = 4.0E0*I_NAI_Hxy3z_G3yz_ab;
  abcd[1304] = 4.0E0*I_NAI_Hx4z_G3yz_ab;
  abcd[1305] = 4.0E0*I_NAI_H5x_G2y2z_ab-2.0E0*1*I_NAI_H5x_D2y_a-2.0E0*4*I_NAI_F3x_G2y2z_b+4*1*I_NAI_F3x_D2y;
  abcd[1306] = 4.0E0*I_NAI_H4xy_G2y2z_ab-2.0E0*1*I_NAI_H4xy_D2y_a-2.0E0*3*I_NAI_F2xy_G2y2z_b+3*1*I_NAI_F2xy_D2y;
  abcd[1307] = 4.0E0*I_NAI_H4xz_G2y2z_ab-2.0E0*1*I_NAI_H4xz_D2y_a-2.0E0*3*I_NAI_F2xz_G2y2z_b+3*1*I_NAI_F2xz_D2y;
  abcd[1308] = 4.0E0*I_NAI_H3x2y_G2y2z_ab-2.0E0*1*I_NAI_H3x2y_D2y_a-2.0E0*2*I_NAI_Fx2y_G2y2z_b+2*1*I_NAI_Fx2y_D2y;
  abcd[1309] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2y_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2y;
  abcd[1310] = 4.0E0*I_NAI_H3x2z_G2y2z_ab-2.0E0*1*I_NAI_H3x2z_D2y_a-2.0E0*2*I_NAI_Fx2z_G2y2z_b+2*1*I_NAI_Fx2z_D2y;
  abcd[1311] = 4.0E0*I_NAI_H2x3y_G2y2z_ab-2.0E0*1*I_NAI_H2x3y_D2y_a-2.0E0*1*I_NAI_F3y_G2y2z_b+1*I_NAI_F3y_D2y;
  abcd[1312] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2yz_G2y2z_b+1*I_NAI_F2yz_D2y;
  abcd[1313] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_Fy2z_G2y2z_b+1*I_NAI_Fy2z_D2y;
  abcd[1314] = 4.0E0*I_NAI_H2x3z_G2y2z_ab-2.0E0*1*I_NAI_H2x3z_D2y_a-2.0E0*1*I_NAI_F3z_G2y2z_b+1*I_NAI_F3z_D2y;
  abcd[1315] = 4.0E0*I_NAI_Hx4y_G2y2z_ab-2.0E0*1*I_NAI_Hx4y_D2y_a;
  abcd[1316] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a;
  abcd[1317] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a;
  abcd[1318] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a;
  abcd[1319] = 4.0E0*I_NAI_Hx4z_G2y2z_ab-2.0E0*1*I_NAI_Hx4z_D2y_a;
  abcd[1320] = 4.0E0*I_NAI_H5x_Gy3z_ab-2.0E0*2*I_NAI_H5x_Dyz_a-2.0E0*4*I_NAI_F3x_Gy3z_b+4*2*I_NAI_F3x_Dyz;
  abcd[1321] = 4.0E0*I_NAI_H4xy_Gy3z_ab-2.0E0*2*I_NAI_H4xy_Dyz_a-2.0E0*3*I_NAI_F2xy_Gy3z_b+3*2*I_NAI_F2xy_Dyz;
  abcd[1322] = 4.0E0*I_NAI_H4xz_Gy3z_ab-2.0E0*2*I_NAI_H4xz_Dyz_a-2.0E0*3*I_NAI_F2xz_Gy3z_b+3*2*I_NAI_F2xz_Dyz;
  abcd[1323] = 4.0E0*I_NAI_H3x2y_Gy3z_ab-2.0E0*2*I_NAI_H3x2y_Dyz_a-2.0E0*2*I_NAI_Fx2y_Gy3z_b+2*2*I_NAI_Fx2y_Dyz;
  abcd[1324] = 4.0E0*I_NAI_H3xyz_Gy3z_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a-2.0E0*2*I_NAI_Fxyz_Gy3z_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[1325] = 4.0E0*I_NAI_H3x2z_Gy3z_ab-2.0E0*2*I_NAI_H3x2z_Dyz_a-2.0E0*2*I_NAI_Fx2z_Gy3z_b+2*2*I_NAI_Fx2z_Dyz;
  abcd[1326] = 4.0E0*I_NAI_H2x3y_Gy3z_ab-2.0E0*2*I_NAI_H2x3y_Dyz_a-2.0E0*1*I_NAI_F3y_Gy3z_b+2*I_NAI_F3y_Dyz;
  abcd[1327] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2yz_Gy3z_b+2*I_NAI_F2yz_Dyz;
  abcd[1328] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_Fy2z_Gy3z_b+2*I_NAI_Fy2z_Dyz;
  abcd[1329] = 4.0E0*I_NAI_H2x3z_Gy3z_ab-2.0E0*2*I_NAI_H2x3z_Dyz_a-2.0E0*1*I_NAI_F3z_Gy3z_b+2*I_NAI_F3z_Dyz;
  abcd[1330] = 4.0E0*I_NAI_Hx4y_Gy3z_ab-2.0E0*2*I_NAI_Hx4y_Dyz_a;
  abcd[1331] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a;
  abcd[1332] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a;
  abcd[1333] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a;
  abcd[1334] = 4.0E0*I_NAI_Hx4z_Gy3z_ab-2.0E0*2*I_NAI_Hx4z_Dyz_a;
  abcd[1335] = 4.0E0*I_NAI_H5x_G4z_ab-2.0E0*3*I_NAI_H5x_D2z_a-2.0E0*4*I_NAI_F3x_G4z_b+4*3*I_NAI_F3x_D2z;
  abcd[1336] = 4.0E0*I_NAI_H4xy_G4z_ab-2.0E0*3*I_NAI_H4xy_D2z_a-2.0E0*3*I_NAI_F2xy_G4z_b+3*3*I_NAI_F2xy_D2z;
  abcd[1337] = 4.0E0*I_NAI_H4xz_G4z_ab-2.0E0*3*I_NAI_H4xz_D2z_a-2.0E0*3*I_NAI_F2xz_G4z_b+3*3*I_NAI_F2xz_D2z;
  abcd[1338] = 4.0E0*I_NAI_H3x2y_G4z_ab-2.0E0*3*I_NAI_H3x2y_D2z_a-2.0E0*2*I_NAI_Fx2y_G4z_b+2*3*I_NAI_Fx2y_D2z;
  abcd[1339] = 4.0E0*I_NAI_H3xyz_G4z_ab-2.0E0*3*I_NAI_H3xyz_D2z_a-2.0E0*2*I_NAI_Fxyz_G4z_b+2*3*I_NAI_Fxyz_D2z;
  abcd[1340] = 4.0E0*I_NAI_H3x2z_G4z_ab-2.0E0*3*I_NAI_H3x2z_D2z_a-2.0E0*2*I_NAI_Fx2z_G4z_b+2*3*I_NAI_Fx2z_D2z;
  abcd[1341] = 4.0E0*I_NAI_H2x3y_G4z_ab-2.0E0*3*I_NAI_H2x3y_D2z_a-2.0E0*1*I_NAI_F3y_G4z_b+3*I_NAI_F3y_D2z;
  abcd[1342] = 4.0E0*I_NAI_H2x2yz_G4z_ab-2.0E0*3*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2yz_G4z_b+3*I_NAI_F2yz_D2z;
  abcd[1343] = 4.0E0*I_NAI_H2xy2z_G4z_ab-2.0E0*3*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_Fy2z_G4z_b+3*I_NAI_Fy2z_D2z;
  abcd[1344] = 4.0E0*I_NAI_H2x3z_G4z_ab-2.0E0*3*I_NAI_H2x3z_D2z_a-2.0E0*1*I_NAI_F3z_G4z_b+3*I_NAI_F3z_D2z;
  abcd[1345] = 4.0E0*I_NAI_Hx4y_G4z_ab-2.0E0*3*I_NAI_Hx4y_D2z_a;
  abcd[1346] = 4.0E0*I_NAI_Hx3yz_G4z_ab-2.0E0*3*I_NAI_Hx3yz_D2z_a;
  abcd[1347] = 4.0E0*I_NAI_Hx2y2z_G4z_ab-2.0E0*3*I_NAI_Hx2y2z_D2z_a;
  abcd[1348] = 4.0E0*I_NAI_Hxy3z_G4z_ab-2.0E0*3*I_NAI_Hxy3z_D2z_a;
  abcd[1349] = 4.0E0*I_NAI_Hx4z_G4z_ab-2.0E0*3*I_NAI_Hx4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1350] = 4.0E0*I_NAI_H4xy_G4x_ab-2.0E0*3*I_NAI_H4xy_D2x_a;
  abcd[1351] = 4.0E0*I_NAI_H3x2y_G4x_ab-2.0E0*3*I_NAI_H3x2y_D2x_a-2.0E0*1*I_NAI_F3x_G4x_b+3*I_NAI_F3x_D2x;
  abcd[1352] = 4.0E0*I_NAI_H3xyz_G4x_ab-2.0E0*3*I_NAI_H3xyz_D2x_a;
  abcd[1353] = 4.0E0*I_NAI_H2x3y_G4x_ab-2.0E0*3*I_NAI_H2x3y_D2x_a-2.0E0*2*I_NAI_F2xy_G4x_b+2*3*I_NAI_F2xy_D2x;
  abcd[1354] = 4.0E0*I_NAI_H2x2yz_G4x_ab-2.0E0*3*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2xz_G4x_b+3*I_NAI_F2xz_D2x;
  abcd[1355] = 4.0E0*I_NAI_H2xy2z_G4x_ab-2.0E0*3*I_NAI_H2xy2z_D2x_a;
  abcd[1356] = 4.0E0*I_NAI_Hx4y_G4x_ab-2.0E0*3*I_NAI_Hx4y_D2x_a-2.0E0*3*I_NAI_Fx2y_G4x_b+3*3*I_NAI_Fx2y_D2x;
  abcd[1357] = 4.0E0*I_NAI_Hx3yz_G4x_ab-2.0E0*3*I_NAI_Hx3yz_D2x_a-2.0E0*2*I_NAI_Fxyz_G4x_b+2*3*I_NAI_Fxyz_D2x;
  abcd[1358] = 4.0E0*I_NAI_Hx2y2z_G4x_ab-2.0E0*3*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2z_G4x_b+3*I_NAI_Fx2z_D2x;
  abcd[1359] = 4.0E0*I_NAI_Hxy3z_G4x_ab-2.0E0*3*I_NAI_Hxy3z_D2x_a;
  abcd[1360] = 4.0E0*I_NAI_H5y_G4x_ab-2.0E0*3*I_NAI_H5y_D2x_a-2.0E0*4*I_NAI_F3y_G4x_b+4*3*I_NAI_F3y_D2x;
  abcd[1361] = 4.0E0*I_NAI_H4yz_G4x_ab-2.0E0*3*I_NAI_H4yz_D2x_a-2.0E0*3*I_NAI_F2yz_G4x_b+3*3*I_NAI_F2yz_D2x;
  abcd[1362] = 4.0E0*I_NAI_H3y2z_G4x_ab-2.0E0*3*I_NAI_H3y2z_D2x_a-2.0E0*2*I_NAI_Fy2z_G4x_b+2*3*I_NAI_Fy2z_D2x;
  abcd[1363] = 4.0E0*I_NAI_H2y3z_G4x_ab-2.0E0*3*I_NAI_H2y3z_D2x_a-2.0E0*1*I_NAI_F3z_G4x_b+3*I_NAI_F3z_D2x;
  abcd[1364] = 4.0E0*I_NAI_Hy4z_G4x_ab-2.0E0*3*I_NAI_Hy4z_D2x_a;
  abcd[1365] = 4.0E0*I_NAI_H4xy_G3xy_ab-2.0E0*2*I_NAI_H4xy_Dxy_a;
  abcd[1366] = 4.0E0*I_NAI_H3x2y_G3xy_ab-2.0E0*2*I_NAI_H3x2y_Dxy_a-2.0E0*1*I_NAI_F3x_G3xy_b+2*I_NAI_F3x_Dxy;
  abcd[1367] = 4.0E0*I_NAI_H3xyz_G3xy_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a;
  abcd[1368] = 4.0E0*I_NAI_H2x3y_G3xy_ab-2.0E0*2*I_NAI_H2x3y_Dxy_a-2.0E0*2*I_NAI_F2xy_G3xy_b+2*2*I_NAI_F2xy_Dxy;
  abcd[1369] = 4.0E0*I_NAI_H2x2yz_G3xy_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2xz_G3xy_b+2*I_NAI_F2xz_Dxy;
  abcd[1370] = 4.0E0*I_NAI_H2xy2z_G3xy_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a;
  abcd[1371] = 4.0E0*I_NAI_Hx4y_G3xy_ab-2.0E0*2*I_NAI_Hx4y_Dxy_a-2.0E0*3*I_NAI_Fx2y_G3xy_b+3*2*I_NAI_Fx2y_Dxy;
  abcd[1372] = 4.0E0*I_NAI_Hx3yz_G3xy_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a-2.0E0*2*I_NAI_Fxyz_G3xy_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[1373] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2z_G3xy_b+2*I_NAI_Fx2z_Dxy;
  abcd[1374] = 4.0E0*I_NAI_Hxy3z_G3xy_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a;
  abcd[1375] = 4.0E0*I_NAI_H5y_G3xy_ab-2.0E0*2*I_NAI_H5y_Dxy_a-2.0E0*4*I_NAI_F3y_G3xy_b+4*2*I_NAI_F3y_Dxy;
  abcd[1376] = 4.0E0*I_NAI_H4yz_G3xy_ab-2.0E0*2*I_NAI_H4yz_Dxy_a-2.0E0*3*I_NAI_F2yz_G3xy_b+3*2*I_NAI_F2yz_Dxy;
  abcd[1377] = 4.0E0*I_NAI_H3y2z_G3xy_ab-2.0E0*2*I_NAI_H3y2z_Dxy_a-2.0E0*2*I_NAI_Fy2z_G3xy_b+2*2*I_NAI_Fy2z_Dxy;
  abcd[1378] = 4.0E0*I_NAI_H2y3z_G3xy_ab-2.0E0*2*I_NAI_H2y3z_Dxy_a-2.0E0*1*I_NAI_F3z_G3xy_b+2*I_NAI_F3z_Dxy;
  abcd[1379] = 4.0E0*I_NAI_Hy4z_G3xy_ab-2.0E0*2*I_NAI_Hy4z_Dxy_a;
  abcd[1380] = 4.0E0*I_NAI_H4xy_G3xz_ab-2.0E0*2*I_NAI_H4xy_Dxz_a;
  abcd[1381] = 4.0E0*I_NAI_H3x2y_G3xz_ab-2.0E0*2*I_NAI_H3x2y_Dxz_a-2.0E0*1*I_NAI_F3x_G3xz_b+2*I_NAI_F3x_Dxz;
  abcd[1382] = 4.0E0*I_NAI_H3xyz_G3xz_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a;
  abcd[1383] = 4.0E0*I_NAI_H2x3y_G3xz_ab-2.0E0*2*I_NAI_H2x3y_Dxz_a-2.0E0*2*I_NAI_F2xy_G3xz_b+2*2*I_NAI_F2xy_Dxz;
  abcd[1384] = 4.0E0*I_NAI_H2x2yz_G3xz_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2xz_G3xz_b+2*I_NAI_F2xz_Dxz;
  abcd[1385] = 4.0E0*I_NAI_H2xy2z_G3xz_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a;
  abcd[1386] = 4.0E0*I_NAI_Hx4y_G3xz_ab-2.0E0*2*I_NAI_Hx4y_Dxz_a-2.0E0*3*I_NAI_Fx2y_G3xz_b+3*2*I_NAI_Fx2y_Dxz;
  abcd[1387] = 4.0E0*I_NAI_Hx3yz_G3xz_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a-2.0E0*2*I_NAI_Fxyz_G3xz_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[1388] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2z_G3xz_b+2*I_NAI_Fx2z_Dxz;
  abcd[1389] = 4.0E0*I_NAI_Hxy3z_G3xz_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a;
  abcd[1390] = 4.0E0*I_NAI_H5y_G3xz_ab-2.0E0*2*I_NAI_H5y_Dxz_a-2.0E0*4*I_NAI_F3y_G3xz_b+4*2*I_NAI_F3y_Dxz;
  abcd[1391] = 4.0E0*I_NAI_H4yz_G3xz_ab-2.0E0*2*I_NAI_H4yz_Dxz_a-2.0E0*3*I_NAI_F2yz_G3xz_b+3*2*I_NAI_F2yz_Dxz;
  abcd[1392] = 4.0E0*I_NAI_H3y2z_G3xz_ab-2.0E0*2*I_NAI_H3y2z_Dxz_a-2.0E0*2*I_NAI_Fy2z_G3xz_b+2*2*I_NAI_Fy2z_Dxz;
  abcd[1393] = 4.0E0*I_NAI_H2y3z_G3xz_ab-2.0E0*2*I_NAI_H2y3z_Dxz_a-2.0E0*1*I_NAI_F3z_G3xz_b+2*I_NAI_F3z_Dxz;
  abcd[1394] = 4.0E0*I_NAI_Hy4z_G3xz_ab-2.0E0*2*I_NAI_Hy4z_Dxz_a;
  abcd[1395] = 4.0E0*I_NAI_H4xy_G2x2y_ab-2.0E0*1*I_NAI_H4xy_D2y_a;
  abcd[1396] = 4.0E0*I_NAI_H3x2y_G2x2y_ab-2.0E0*1*I_NAI_H3x2y_D2y_a-2.0E0*1*I_NAI_F3x_G2x2y_b+1*I_NAI_F3x_D2y;
  abcd[1397] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2y_a;
  abcd[1398] = 4.0E0*I_NAI_H2x3y_G2x2y_ab-2.0E0*1*I_NAI_H2x3y_D2y_a-2.0E0*2*I_NAI_F2xy_G2x2y_b+2*1*I_NAI_F2xy_D2y;
  abcd[1399] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2xz_G2x2y_b+1*I_NAI_F2xz_D2y;
  abcd[1400] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a;
  abcd[1401] = 4.0E0*I_NAI_Hx4y_G2x2y_ab-2.0E0*1*I_NAI_Hx4y_D2y_a-2.0E0*3*I_NAI_Fx2y_G2x2y_b+3*1*I_NAI_Fx2y_D2y;
  abcd[1402] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2y;
  abcd[1403] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2z_G2x2y_b+1*I_NAI_Fx2z_D2y;
  abcd[1404] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a;
  abcd[1405] = 4.0E0*I_NAI_H5y_G2x2y_ab-2.0E0*1*I_NAI_H5y_D2y_a-2.0E0*4*I_NAI_F3y_G2x2y_b+4*1*I_NAI_F3y_D2y;
  abcd[1406] = 4.0E0*I_NAI_H4yz_G2x2y_ab-2.0E0*1*I_NAI_H4yz_D2y_a-2.0E0*3*I_NAI_F2yz_G2x2y_b+3*1*I_NAI_F2yz_D2y;
  abcd[1407] = 4.0E0*I_NAI_H3y2z_G2x2y_ab-2.0E0*1*I_NAI_H3y2z_D2y_a-2.0E0*2*I_NAI_Fy2z_G2x2y_b+2*1*I_NAI_Fy2z_D2y;
  abcd[1408] = 4.0E0*I_NAI_H2y3z_G2x2y_ab-2.0E0*1*I_NAI_H2y3z_D2y_a-2.0E0*1*I_NAI_F3z_G2x2y_b+1*I_NAI_F3z_D2y;
  abcd[1409] = 4.0E0*I_NAI_Hy4z_G2x2y_ab-2.0E0*1*I_NAI_Hy4z_D2y_a;
  abcd[1410] = 4.0E0*I_NAI_H4xy_G2xyz_ab-2.0E0*1*I_NAI_H4xy_Dyz_a;
  abcd[1411] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*1*I_NAI_H3x2y_Dyz_a-2.0E0*1*I_NAI_F3x_G2xyz_b+1*I_NAI_F3x_Dyz;
  abcd[1412] = 4.0E0*I_NAI_H3xyz_G2xyz_ab-2.0E0*1*I_NAI_H3xyz_Dyz_a;
  abcd[1413] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*1*I_NAI_H2x3y_Dyz_a-2.0E0*2*I_NAI_F2xy_G2xyz_b+2*1*I_NAI_F2xy_Dyz;
  abcd[1414] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2xz_G2xyz_b+1*I_NAI_F2xz_Dyz;
  abcd[1415] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_H2xy2z_Dyz_a;
  abcd[1416] = 4.0E0*I_NAI_Hx4y_G2xyz_ab-2.0E0*1*I_NAI_Hx4y_Dyz_a-2.0E0*3*I_NAI_Fx2y_G2xyz_b+3*1*I_NAI_Fx2y_Dyz;
  abcd[1417] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab-2.0E0*1*I_NAI_Hx3yz_Dyz_a-2.0E0*2*I_NAI_Fxyz_G2xyz_b+2*1*I_NAI_Fxyz_Dyz;
  abcd[1418] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2z_G2xyz_b+1*I_NAI_Fx2z_Dyz;
  abcd[1419] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab-2.0E0*1*I_NAI_Hxy3z_Dyz_a;
  abcd[1420] = 4.0E0*I_NAI_H5y_G2xyz_ab-2.0E0*1*I_NAI_H5y_Dyz_a-2.0E0*4*I_NAI_F3y_G2xyz_b+4*1*I_NAI_F3y_Dyz;
  abcd[1421] = 4.0E0*I_NAI_H4yz_G2xyz_ab-2.0E0*1*I_NAI_H4yz_Dyz_a-2.0E0*3*I_NAI_F2yz_G2xyz_b+3*1*I_NAI_F2yz_Dyz;
  abcd[1422] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*1*I_NAI_H3y2z_Dyz_a-2.0E0*2*I_NAI_Fy2z_G2xyz_b+2*1*I_NAI_Fy2z_Dyz;
  abcd[1423] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*1*I_NAI_H2y3z_Dyz_a-2.0E0*1*I_NAI_F3z_G2xyz_b+1*I_NAI_F3z_Dyz;
  abcd[1424] = 4.0E0*I_NAI_Hy4z_G2xyz_ab-2.0E0*1*I_NAI_Hy4z_Dyz_a;
  abcd[1425] = 4.0E0*I_NAI_H4xy_G2x2z_ab-2.0E0*1*I_NAI_H4xy_D2z_a;
  abcd[1426] = 4.0E0*I_NAI_H3x2y_G2x2z_ab-2.0E0*1*I_NAI_H3x2y_D2z_a-2.0E0*1*I_NAI_F3x_G2x2z_b+1*I_NAI_F3x_D2z;
  abcd[1427] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a;
  abcd[1428] = 4.0E0*I_NAI_H2x3y_G2x2z_ab-2.0E0*1*I_NAI_H2x3y_D2z_a-2.0E0*2*I_NAI_F2xy_G2x2z_b+2*1*I_NAI_F2xy_D2z;
  abcd[1429] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2xz_G2x2z_b+1*I_NAI_F2xz_D2z;
  abcd[1430] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a;
  abcd[1431] = 4.0E0*I_NAI_Hx4y_G2x2z_ab-2.0E0*1*I_NAI_Hx4y_D2z_a-2.0E0*3*I_NAI_Fx2y_G2x2z_b+3*1*I_NAI_Fx2y_D2z;
  abcd[1432] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[1433] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2z_G2x2z_b+1*I_NAI_Fx2z_D2z;
  abcd[1434] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a;
  abcd[1435] = 4.0E0*I_NAI_H5y_G2x2z_ab-2.0E0*1*I_NAI_H5y_D2z_a-2.0E0*4*I_NAI_F3y_G2x2z_b+4*1*I_NAI_F3y_D2z;
  abcd[1436] = 4.0E0*I_NAI_H4yz_G2x2z_ab-2.0E0*1*I_NAI_H4yz_D2z_a-2.0E0*3*I_NAI_F2yz_G2x2z_b+3*1*I_NAI_F2yz_D2z;
  abcd[1437] = 4.0E0*I_NAI_H3y2z_G2x2z_ab-2.0E0*1*I_NAI_H3y2z_D2z_a-2.0E0*2*I_NAI_Fy2z_G2x2z_b+2*1*I_NAI_Fy2z_D2z;
  abcd[1438] = 4.0E0*I_NAI_H2y3z_G2x2z_ab-2.0E0*1*I_NAI_H2y3z_D2z_a-2.0E0*1*I_NAI_F3z_G2x2z_b+1*I_NAI_F3z_D2z;
  abcd[1439] = 4.0E0*I_NAI_Hy4z_G2x2z_ab-2.0E0*1*I_NAI_Hy4z_D2z_a;
  abcd[1440] = 4.0E0*I_NAI_H4xy_Gx3y_ab;
  abcd[1441] = 4.0E0*I_NAI_H3x2y_Gx3y_ab-2.0E0*1*I_NAI_F3x_Gx3y_b;
  abcd[1442] = 4.0E0*I_NAI_H3xyz_Gx3y_ab;
  abcd[1443] = 4.0E0*I_NAI_H2x3y_Gx3y_ab-2.0E0*2*I_NAI_F2xy_Gx3y_b;
  abcd[1444] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab-2.0E0*1*I_NAI_F2xz_Gx3y_b;
  abcd[1445] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab;
  abcd[1446] = 4.0E0*I_NAI_Hx4y_Gx3y_ab-2.0E0*3*I_NAI_Fx2y_Gx3y_b;
  abcd[1447] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab-2.0E0*2*I_NAI_Fxyz_Gx3y_b;
  abcd[1448] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab-2.0E0*1*I_NAI_Fx2z_Gx3y_b;
  abcd[1449] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab;
  abcd[1450] = 4.0E0*I_NAI_H5y_Gx3y_ab-2.0E0*4*I_NAI_F3y_Gx3y_b;
  abcd[1451] = 4.0E0*I_NAI_H4yz_Gx3y_ab-2.0E0*3*I_NAI_F2yz_Gx3y_b;
  abcd[1452] = 4.0E0*I_NAI_H3y2z_Gx3y_ab-2.0E0*2*I_NAI_Fy2z_Gx3y_b;
  abcd[1453] = 4.0E0*I_NAI_H2y3z_Gx3y_ab-2.0E0*1*I_NAI_F3z_Gx3y_b;
  abcd[1454] = 4.0E0*I_NAI_Hy4z_Gx3y_ab;
  abcd[1455] = 4.0E0*I_NAI_H4xy_Gx2yz_ab;
  abcd[1456] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*1*I_NAI_F3x_Gx2yz_b;
  abcd[1457] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab;
  abcd[1458] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*2*I_NAI_F2xy_Gx2yz_b;
  abcd[1459] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_F2xz_Gx2yz_b;
  abcd[1460] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab;
  abcd[1461] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab-2.0E0*3*I_NAI_Fx2y_Gx2yz_b;
  abcd[1462] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[1463] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Fx2z_Gx2yz_b;
  abcd[1464] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab;
  abcd[1465] = 4.0E0*I_NAI_H5y_Gx2yz_ab-2.0E0*4*I_NAI_F3y_Gx2yz_b;
  abcd[1466] = 4.0E0*I_NAI_H4yz_Gx2yz_ab-2.0E0*3*I_NAI_F2yz_Gx2yz_b;
  abcd[1467] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*2*I_NAI_Fy2z_Gx2yz_b;
  abcd[1468] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*1*I_NAI_F3z_Gx2yz_b;
  abcd[1469] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab;
  abcd[1470] = 4.0E0*I_NAI_H4xy_Gxy2z_ab;
  abcd[1471] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*1*I_NAI_F3x_Gxy2z_b;
  abcd[1472] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab;
  abcd[1473] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*2*I_NAI_F2xy_Gxy2z_b;
  abcd[1474] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_F2xz_Gxy2z_b;
  abcd[1475] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab;
  abcd[1476] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab-2.0E0*3*I_NAI_Fx2y_Gxy2z_b;
  abcd[1477] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[1478] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Fx2z_Gxy2z_b;
  abcd[1479] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab;
  abcd[1480] = 4.0E0*I_NAI_H5y_Gxy2z_ab-2.0E0*4*I_NAI_F3y_Gxy2z_b;
  abcd[1481] = 4.0E0*I_NAI_H4yz_Gxy2z_ab-2.0E0*3*I_NAI_F2yz_Gxy2z_b;
  abcd[1482] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*2*I_NAI_Fy2z_Gxy2z_b;
  abcd[1483] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*1*I_NAI_F3z_Gxy2z_b;
  abcd[1484] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab;
  abcd[1485] = 4.0E0*I_NAI_H4xy_Gx3z_ab;
  abcd[1486] = 4.0E0*I_NAI_H3x2y_Gx3z_ab-2.0E0*1*I_NAI_F3x_Gx3z_b;
  abcd[1487] = 4.0E0*I_NAI_H3xyz_Gx3z_ab;
  abcd[1488] = 4.0E0*I_NAI_H2x3y_Gx3z_ab-2.0E0*2*I_NAI_F2xy_Gx3z_b;
  abcd[1489] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab-2.0E0*1*I_NAI_F2xz_Gx3z_b;
  abcd[1490] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab;
  abcd[1491] = 4.0E0*I_NAI_Hx4y_Gx3z_ab-2.0E0*3*I_NAI_Fx2y_Gx3z_b;
  abcd[1492] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab-2.0E0*2*I_NAI_Fxyz_Gx3z_b;
  abcd[1493] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab-2.0E0*1*I_NAI_Fx2z_Gx3z_b;
  abcd[1494] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab;
  abcd[1495] = 4.0E0*I_NAI_H5y_Gx3z_ab-2.0E0*4*I_NAI_F3y_Gx3z_b;
  abcd[1496] = 4.0E0*I_NAI_H4yz_Gx3z_ab-2.0E0*3*I_NAI_F2yz_Gx3z_b;
  abcd[1497] = 4.0E0*I_NAI_H3y2z_Gx3z_ab-2.0E0*2*I_NAI_Fy2z_Gx3z_b;
  abcd[1498] = 4.0E0*I_NAI_H2y3z_Gx3z_ab-2.0E0*1*I_NAI_F3z_Gx3z_b;
  abcd[1499] = 4.0E0*I_NAI_Hy4z_Gx3z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1500] = 4.0E0*I_NAI_H4xy_G3xy_ab;
  abcd[1501] = 4.0E0*I_NAI_H3x2y_G3xy_ab-2.0E0*1*I_NAI_F3x_G3xy_b;
  abcd[1502] = 4.0E0*I_NAI_H3xyz_G3xy_ab;
  abcd[1503] = 4.0E0*I_NAI_H2x3y_G3xy_ab-2.0E0*2*I_NAI_F2xy_G3xy_b;
  abcd[1504] = 4.0E0*I_NAI_H2x2yz_G3xy_ab-2.0E0*1*I_NAI_F2xz_G3xy_b;
  abcd[1505] = 4.0E0*I_NAI_H2xy2z_G3xy_ab;
  abcd[1506] = 4.0E0*I_NAI_Hx4y_G3xy_ab-2.0E0*3*I_NAI_Fx2y_G3xy_b;
  abcd[1507] = 4.0E0*I_NAI_Hx3yz_G3xy_ab-2.0E0*2*I_NAI_Fxyz_G3xy_b;
  abcd[1508] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab-2.0E0*1*I_NAI_Fx2z_G3xy_b;
  abcd[1509] = 4.0E0*I_NAI_Hxy3z_G3xy_ab;
  abcd[1510] = 4.0E0*I_NAI_H5y_G3xy_ab-2.0E0*4*I_NAI_F3y_G3xy_b;
  abcd[1511] = 4.0E0*I_NAI_H4yz_G3xy_ab-2.0E0*3*I_NAI_F2yz_G3xy_b;
  abcd[1512] = 4.0E0*I_NAI_H3y2z_G3xy_ab-2.0E0*2*I_NAI_Fy2z_G3xy_b;
  abcd[1513] = 4.0E0*I_NAI_H2y3z_G3xy_ab-2.0E0*1*I_NAI_F3z_G3xy_b;
  abcd[1514] = 4.0E0*I_NAI_Hy4z_G3xy_ab;
  abcd[1515] = 4.0E0*I_NAI_H4xy_G2x2y_ab-2.0E0*1*I_NAI_H4xy_D2x_a;
  abcd[1516] = 4.0E0*I_NAI_H3x2y_G2x2y_ab-2.0E0*1*I_NAI_H3x2y_D2x_a-2.0E0*1*I_NAI_F3x_G2x2y_b+1*I_NAI_F3x_D2x;
  abcd[1517] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2x_a;
  abcd[1518] = 4.0E0*I_NAI_H2x3y_G2x2y_ab-2.0E0*1*I_NAI_H2x3y_D2x_a-2.0E0*2*I_NAI_F2xy_G2x2y_b+2*1*I_NAI_F2xy_D2x;
  abcd[1519] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2xz_G2x2y_b+1*I_NAI_F2xz_D2x;
  abcd[1520] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a;
  abcd[1521] = 4.0E0*I_NAI_Hx4y_G2x2y_ab-2.0E0*1*I_NAI_Hx4y_D2x_a-2.0E0*3*I_NAI_Fx2y_G2x2y_b+3*1*I_NAI_Fx2y_D2x;
  abcd[1522] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2x;
  abcd[1523] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2z_G2x2y_b+1*I_NAI_Fx2z_D2x;
  abcd[1524] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a;
  abcd[1525] = 4.0E0*I_NAI_H5y_G2x2y_ab-2.0E0*1*I_NAI_H5y_D2x_a-2.0E0*4*I_NAI_F3y_G2x2y_b+4*1*I_NAI_F3y_D2x;
  abcd[1526] = 4.0E0*I_NAI_H4yz_G2x2y_ab-2.0E0*1*I_NAI_H4yz_D2x_a-2.0E0*3*I_NAI_F2yz_G2x2y_b+3*1*I_NAI_F2yz_D2x;
  abcd[1527] = 4.0E0*I_NAI_H3y2z_G2x2y_ab-2.0E0*1*I_NAI_H3y2z_D2x_a-2.0E0*2*I_NAI_Fy2z_G2x2y_b+2*1*I_NAI_Fy2z_D2x;
  abcd[1528] = 4.0E0*I_NAI_H2y3z_G2x2y_ab-2.0E0*1*I_NAI_H2y3z_D2x_a-2.0E0*1*I_NAI_F3z_G2x2y_b+1*I_NAI_F3z_D2x;
  abcd[1529] = 4.0E0*I_NAI_Hy4z_G2x2y_ab-2.0E0*1*I_NAI_Hy4z_D2x_a;
  abcd[1530] = 4.0E0*I_NAI_H4xy_G2xyz_ab;
  abcd[1531] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*1*I_NAI_F3x_G2xyz_b;
  abcd[1532] = 4.0E0*I_NAI_H3xyz_G2xyz_ab;
  abcd[1533] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*2*I_NAI_F2xy_G2xyz_b;
  abcd[1534] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_F2xz_G2xyz_b;
  abcd[1535] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab;
  abcd[1536] = 4.0E0*I_NAI_Hx4y_G2xyz_ab-2.0E0*3*I_NAI_Fx2y_G2xyz_b;
  abcd[1537] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[1538] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Fx2z_G2xyz_b;
  abcd[1539] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab;
  abcd[1540] = 4.0E0*I_NAI_H5y_G2xyz_ab-2.0E0*4*I_NAI_F3y_G2xyz_b;
  abcd[1541] = 4.0E0*I_NAI_H4yz_G2xyz_ab-2.0E0*3*I_NAI_F2yz_G2xyz_b;
  abcd[1542] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*2*I_NAI_Fy2z_G2xyz_b;
  abcd[1543] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*1*I_NAI_F3z_G2xyz_b;
  abcd[1544] = 4.0E0*I_NAI_Hy4z_G2xyz_ab;
  abcd[1545] = 4.0E0*I_NAI_H4xy_Gx3y_ab-2.0E0*2*I_NAI_H4xy_Dxy_a;
  abcd[1546] = 4.0E0*I_NAI_H3x2y_Gx3y_ab-2.0E0*2*I_NAI_H3x2y_Dxy_a-2.0E0*1*I_NAI_F3x_Gx3y_b+2*I_NAI_F3x_Dxy;
  abcd[1547] = 4.0E0*I_NAI_H3xyz_Gx3y_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a;
  abcd[1548] = 4.0E0*I_NAI_H2x3y_Gx3y_ab-2.0E0*2*I_NAI_H2x3y_Dxy_a-2.0E0*2*I_NAI_F2xy_Gx3y_b+2*2*I_NAI_F2xy_Dxy;
  abcd[1549] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2xz_Gx3y_b+2*I_NAI_F2xz_Dxy;
  abcd[1550] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a;
  abcd[1551] = 4.0E0*I_NAI_Hx4y_Gx3y_ab-2.0E0*2*I_NAI_Hx4y_Dxy_a-2.0E0*3*I_NAI_Fx2y_Gx3y_b+3*2*I_NAI_Fx2y_Dxy;
  abcd[1552] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gx3y_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[1553] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2z_Gx3y_b+2*I_NAI_Fx2z_Dxy;
  abcd[1554] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a;
  abcd[1555] = 4.0E0*I_NAI_H5y_Gx3y_ab-2.0E0*2*I_NAI_H5y_Dxy_a-2.0E0*4*I_NAI_F3y_Gx3y_b+4*2*I_NAI_F3y_Dxy;
  abcd[1556] = 4.0E0*I_NAI_H4yz_Gx3y_ab-2.0E0*2*I_NAI_H4yz_Dxy_a-2.0E0*3*I_NAI_F2yz_Gx3y_b+3*2*I_NAI_F2yz_Dxy;
  abcd[1557] = 4.0E0*I_NAI_H3y2z_Gx3y_ab-2.0E0*2*I_NAI_H3y2z_Dxy_a-2.0E0*2*I_NAI_Fy2z_Gx3y_b+2*2*I_NAI_Fy2z_Dxy;
  abcd[1558] = 4.0E0*I_NAI_H2y3z_Gx3y_ab-2.0E0*2*I_NAI_H2y3z_Dxy_a-2.0E0*1*I_NAI_F3z_Gx3y_b+2*I_NAI_F3z_Dxy;
  abcd[1559] = 4.0E0*I_NAI_Hy4z_Gx3y_ab-2.0E0*2*I_NAI_Hy4z_Dxy_a;
  abcd[1560] = 4.0E0*I_NAI_H4xy_Gx2yz_ab-2.0E0*1*I_NAI_H4xy_Dxz_a;
  abcd[1561] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*1*I_NAI_H3x2y_Dxz_a-2.0E0*1*I_NAI_F3x_Gx2yz_b+1*I_NAI_F3x_Dxz;
  abcd[1562] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab-2.0E0*1*I_NAI_H3xyz_Dxz_a;
  abcd[1563] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*1*I_NAI_H2x3y_Dxz_a-2.0E0*2*I_NAI_F2xy_Gx2yz_b+2*1*I_NAI_F2xy_Dxz;
  abcd[1564] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2xz_Gx2yz_b+1*I_NAI_F2xz_Dxz;
  abcd[1565] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_H2xy2z_Dxz_a;
  abcd[1566] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab-2.0E0*1*I_NAI_Hx4y_Dxz_a-2.0E0*3*I_NAI_Fx2y_Gx2yz_b+3*1*I_NAI_Fx2y_Dxz;
  abcd[1567] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab-2.0E0*1*I_NAI_Hx3yz_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx2yz_b+2*1*I_NAI_Fxyz_Dxz;
  abcd[1568] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2z_Gx2yz_b+1*I_NAI_Fx2z_Dxz;
  abcd[1569] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab-2.0E0*1*I_NAI_Hxy3z_Dxz_a;
  abcd[1570] = 4.0E0*I_NAI_H5y_Gx2yz_ab-2.0E0*1*I_NAI_H5y_Dxz_a-2.0E0*4*I_NAI_F3y_Gx2yz_b+4*1*I_NAI_F3y_Dxz;
  abcd[1571] = 4.0E0*I_NAI_H4yz_Gx2yz_ab-2.0E0*1*I_NAI_H4yz_Dxz_a-2.0E0*3*I_NAI_F2yz_Gx2yz_b+3*1*I_NAI_F2yz_Dxz;
  abcd[1572] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*1*I_NAI_H3y2z_Dxz_a-2.0E0*2*I_NAI_Fy2z_Gx2yz_b+2*1*I_NAI_Fy2z_Dxz;
  abcd[1573] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*1*I_NAI_H2y3z_Dxz_a-2.0E0*1*I_NAI_F3z_Gx2yz_b+1*I_NAI_F3z_Dxz;
  abcd[1574] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab-2.0E0*1*I_NAI_Hy4z_Dxz_a;
  abcd[1575] = 4.0E0*I_NAI_H4xy_Gxy2z_ab;
  abcd[1576] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*1*I_NAI_F3x_Gxy2z_b;
  abcd[1577] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab;
  abcd[1578] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*2*I_NAI_F2xy_Gxy2z_b;
  abcd[1579] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_F2xz_Gxy2z_b;
  abcd[1580] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab;
  abcd[1581] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab-2.0E0*3*I_NAI_Fx2y_Gxy2z_b;
  abcd[1582] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[1583] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Fx2z_Gxy2z_b;
  abcd[1584] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab;
  abcd[1585] = 4.0E0*I_NAI_H5y_Gxy2z_ab-2.0E0*4*I_NAI_F3y_Gxy2z_b;
  abcd[1586] = 4.0E0*I_NAI_H4yz_Gxy2z_ab-2.0E0*3*I_NAI_F2yz_Gxy2z_b;
  abcd[1587] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*2*I_NAI_Fy2z_Gxy2z_b;
  abcd[1588] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*1*I_NAI_F3z_Gxy2z_b;
  abcd[1589] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab;
  abcd[1590] = 4.0E0*I_NAI_H4xy_G4y_ab-2.0E0*3*I_NAI_H4xy_D2y_a;
  abcd[1591] = 4.0E0*I_NAI_H3x2y_G4y_ab-2.0E0*3*I_NAI_H3x2y_D2y_a-2.0E0*1*I_NAI_F3x_G4y_b+3*I_NAI_F3x_D2y;
  abcd[1592] = 4.0E0*I_NAI_H3xyz_G4y_ab-2.0E0*3*I_NAI_H3xyz_D2y_a;
  abcd[1593] = 4.0E0*I_NAI_H2x3y_G4y_ab-2.0E0*3*I_NAI_H2x3y_D2y_a-2.0E0*2*I_NAI_F2xy_G4y_b+2*3*I_NAI_F2xy_D2y;
  abcd[1594] = 4.0E0*I_NAI_H2x2yz_G4y_ab-2.0E0*3*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2xz_G4y_b+3*I_NAI_F2xz_D2y;
  abcd[1595] = 4.0E0*I_NAI_H2xy2z_G4y_ab-2.0E0*3*I_NAI_H2xy2z_D2y_a;
  abcd[1596] = 4.0E0*I_NAI_Hx4y_G4y_ab-2.0E0*3*I_NAI_Hx4y_D2y_a-2.0E0*3*I_NAI_Fx2y_G4y_b+3*3*I_NAI_Fx2y_D2y;
  abcd[1597] = 4.0E0*I_NAI_Hx3yz_G4y_ab-2.0E0*3*I_NAI_Hx3yz_D2y_a-2.0E0*2*I_NAI_Fxyz_G4y_b+2*3*I_NAI_Fxyz_D2y;
  abcd[1598] = 4.0E0*I_NAI_Hx2y2z_G4y_ab-2.0E0*3*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2z_G4y_b+3*I_NAI_Fx2z_D2y;
  abcd[1599] = 4.0E0*I_NAI_Hxy3z_G4y_ab-2.0E0*3*I_NAI_Hxy3z_D2y_a;
  abcd[1600] = 4.0E0*I_NAI_H5y_G4y_ab-2.0E0*3*I_NAI_H5y_D2y_a-2.0E0*4*I_NAI_F3y_G4y_b+4*3*I_NAI_F3y_D2y;
  abcd[1601] = 4.0E0*I_NAI_H4yz_G4y_ab-2.0E0*3*I_NAI_H4yz_D2y_a-2.0E0*3*I_NAI_F2yz_G4y_b+3*3*I_NAI_F2yz_D2y;
  abcd[1602] = 4.0E0*I_NAI_H3y2z_G4y_ab-2.0E0*3*I_NAI_H3y2z_D2y_a-2.0E0*2*I_NAI_Fy2z_G4y_b+2*3*I_NAI_Fy2z_D2y;
  abcd[1603] = 4.0E0*I_NAI_H2y3z_G4y_ab-2.0E0*3*I_NAI_H2y3z_D2y_a-2.0E0*1*I_NAI_F3z_G4y_b+3*I_NAI_F3z_D2y;
  abcd[1604] = 4.0E0*I_NAI_Hy4z_G4y_ab-2.0E0*3*I_NAI_Hy4z_D2y_a;
  abcd[1605] = 4.0E0*I_NAI_H4xy_G3yz_ab-2.0E0*2*I_NAI_H4xy_Dyz_a;
  abcd[1606] = 4.0E0*I_NAI_H3x2y_G3yz_ab-2.0E0*2*I_NAI_H3x2y_Dyz_a-2.0E0*1*I_NAI_F3x_G3yz_b+2*I_NAI_F3x_Dyz;
  abcd[1607] = 4.0E0*I_NAI_H3xyz_G3yz_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a;
  abcd[1608] = 4.0E0*I_NAI_H2x3y_G3yz_ab-2.0E0*2*I_NAI_H2x3y_Dyz_a-2.0E0*2*I_NAI_F2xy_G3yz_b+2*2*I_NAI_F2xy_Dyz;
  abcd[1609] = 4.0E0*I_NAI_H2x2yz_G3yz_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2xz_G3yz_b+2*I_NAI_F2xz_Dyz;
  abcd[1610] = 4.0E0*I_NAI_H2xy2z_G3yz_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a;
  abcd[1611] = 4.0E0*I_NAI_Hx4y_G3yz_ab-2.0E0*2*I_NAI_Hx4y_Dyz_a-2.0E0*3*I_NAI_Fx2y_G3yz_b+3*2*I_NAI_Fx2y_Dyz;
  abcd[1612] = 4.0E0*I_NAI_Hx3yz_G3yz_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a-2.0E0*2*I_NAI_Fxyz_G3yz_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[1613] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2z_G3yz_b+2*I_NAI_Fx2z_Dyz;
  abcd[1614] = 4.0E0*I_NAI_Hxy3z_G3yz_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a;
  abcd[1615] = 4.0E0*I_NAI_H5y_G3yz_ab-2.0E0*2*I_NAI_H5y_Dyz_a-2.0E0*4*I_NAI_F3y_G3yz_b+4*2*I_NAI_F3y_Dyz;
  abcd[1616] = 4.0E0*I_NAI_H4yz_G3yz_ab-2.0E0*2*I_NAI_H4yz_Dyz_a-2.0E0*3*I_NAI_F2yz_G3yz_b+3*2*I_NAI_F2yz_Dyz;
  abcd[1617] = 4.0E0*I_NAI_H3y2z_G3yz_ab-2.0E0*2*I_NAI_H3y2z_Dyz_a-2.0E0*2*I_NAI_Fy2z_G3yz_b+2*2*I_NAI_Fy2z_Dyz;
  abcd[1618] = 4.0E0*I_NAI_H2y3z_G3yz_ab-2.0E0*2*I_NAI_H2y3z_Dyz_a-2.0E0*1*I_NAI_F3z_G3yz_b+2*I_NAI_F3z_Dyz;
  abcd[1619] = 4.0E0*I_NAI_Hy4z_G3yz_ab-2.0E0*2*I_NAI_Hy4z_Dyz_a;
  abcd[1620] = 4.0E0*I_NAI_H4xy_G2y2z_ab-2.0E0*1*I_NAI_H4xy_D2z_a;
  abcd[1621] = 4.0E0*I_NAI_H3x2y_G2y2z_ab-2.0E0*1*I_NAI_H3x2y_D2z_a-2.0E0*1*I_NAI_F3x_G2y2z_b+1*I_NAI_F3x_D2z;
  abcd[1622] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a;
  abcd[1623] = 4.0E0*I_NAI_H2x3y_G2y2z_ab-2.0E0*1*I_NAI_H2x3y_D2z_a-2.0E0*2*I_NAI_F2xy_G2y2z_b+2*1*I_NAI_F2xy_D2z;
  abcd[1624] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2xz_G2y2z_b+1*I_NAI_F2xz_D2z;
  abcd[1625] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a;
  abcd[1626] = 4.0E0*I_NAI_Hx4y_G2y2z_ab-2.0E0*1*I_NAI_Hx4y_D2z_a-2.0E0*3*I_NAI_Fx2y_G2y2z_b+3*1*I_NAI_Fx2y_D2z;
  abcd[1627] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[1628] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2z_G2y2z_b+1*I_NAI_Fx2z_D2z;
  abcd[1629] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a;
  abcd[1630] = 4.0E0*I_NAI_H5y_G2y2z_ab-2.0E0*1*I_NAI_H5y_D2z_a-2.0E0*4*I_NAI_F3y_G2y2z_b+4*1*I_NAI_F3y_D2z;
  abcd[1631] = 4.0E0*I_NAI_H4yz_G2y2z_ab-2.0E0*1*I_NAI_H4yz_D2z_a-2.0E0*3*I_NAI_F2yz_G2y2z_b+3*1*I_NAI_F2yz_D2z;
  abcd[1632] = 4.0E0*I_NAI_H3y2z_G2y2z_ab-2.0E0*1*I_NAI_H3y2z_D2z_a-2.0E0*2*I_NAI_Fy2z_G2y2z_b+2*1*I_NAI_Fy2z_D2z;
  abcd[1633] = 4.0E0*I_NAI_H2y3z_G2y2z_ab-2.0E0*1*I_NAI_H2y3z_D2z_a-2.0E0*1*I_NAI_F3z_G2y2z_b+1*I_NAI_F3z_D2z;
  abcd[1634] = 4.0E0*I_NAI_Hy4z_G2y2z_ab-2.0E0*1*I_NAI_Hy4z_D2z_a;
  abcd[1635] = 4.0E0*I_NAI_H4xy_Gy3z_ab;
  abcd[1636] = 4.0E0*I_NAI_H3x2y_Gy3z_ab-2.0E0*1*I_NAI_F3x_Gy3z_b;
  abcd[1637] = 4.0E0*I_NAI_H3xyz_Gy3z_ab;
  abcd[1638] = 4.0E0*I_NAI_H2x3y_Gy3z_ab-2.0E0*2*I_NAI_F2xy_Gy3z_b;
  abcd[1639] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab-2.0E0*1*I_NAI_F2xz_Gy3z_b;
  abcd[1640] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab;
  abcd[1641] = 4.0E0*I_NAI_Hx4y_Gy3z_ab-2.0E0*3*I_NAI_Fx2y_Gy3z_b;
  abcd[1642] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab-2.0E0*2*I_NAI_Fxyz_Gy3z_b;
  abcd[1643] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab-2.0E0*1*I_NAI_Fx2z_Gy3z_b;
  abcd[1644] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab;
  abcd[1645] = 4.0E0*I_NAI_H5y_Gy3z_ab-2.0E0*4*I_NAI_F3y_Gy3z_b;
  abcd[1646] = 4.0E0*I_NAI_H4yz_Gy3z_ab-2.0E0*3*I_NAI_F2yz_Gy3z_b;
  abcd[1647] = 4.0E0*I_NAI_H3y2z_Gy3z_ab-2.0E0*2*I_NAI_Fy2z_Gy3z_b;
  abcd[1648] = 4.0E0*I_NAI_H2y3z_Gy3z_ab-2.0E0*1*I_NAI_F3z_Gy3z_b;
  abcd[1649] = 4.0E0*I_NAI_Hy4z_Gy3z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1650] = 4.0E0*I_NAI_H4xy_G3xz_ab;
  abcd[1651] = 4.0E0*I_NAI_H3x2y_G3xz_ab-2.0E0*1*I_NAI_F3x_G3xz_b;
  abcd[1652] = 4.0E0*I_NAI_H3xyz_G3xz_ab;
  abcd[1653] = 4.0E0*I_NAI_H2x3y_G3xz_ab-2.0E0*2*I_NAI_F2xy_G3xz_b;
  abcd[1654] = 4.0E0*I_NAI_H2x2yz_G3xz_ab-2.0E0*1*I_NAI_F2xz_G3xz_b;
  abcd[1655] = 4.0E0*I_NAI_H2xy2z_G3xz_ab;
  abcd[1656] = 4.0E0*I_NAI_Hx4y_G3xz_ab-2.0E0*3*I_NAI_Fx2y_G3xz_b;
  abcd[1657] = 4.0E0*I_NAI_Hx3yz_G3xz_ab-2.0E0*2*I_NAI_Fxyz_G3xz_b;
  abcd[1658] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab-2.0E0*1*I_NAI_Fx2z_G3xz_b;
  abcd[1659] = 4.0E0*I_NAI_Hxy3z_G3xz_ab;
  abcd[1660] = 4.0E0*I_NAI_H5y_G3xz_ab-2.0E0*4*I_NAI_F3y_G3xz_b;
  abcd[1661] = 4.0E0*I_NAI_H4yz_G3xz_ab-2.0E0*3*I_NAI_F2yz_G3xz_b;
  abcd[1662] = 4.0E0*I_NAI_H3y2z_G3xz_ab-2.0E0*2*I_NAI_Fy2z_G3xz_b;
  abcd[1663] = 4.0E0*I_NAI_H2y3z_G3xz_ab-2.0E0*1*I_NAI_F3z_G3xz_b;
  abcd[1664] = 4.0E0*I_NAI_Hy4z_G3xz_ab;
  abcd[1665] = 4.0E0*I_NAI_H4xy_G2xyz_ab;
  abcd[1666] = 4.0E0*I_NAI_H3x2y_G2xyz_ab-2.0E0*1*I_NAI_F3x_G2xyz_b;
  abcd[1667] = 4.0E0*I_NAI_H3xyz_G2xyz_ab;
  abcd[1668] = 4.0E0*I_NAI_H2x3y_G2xyz_ab-2.0E0*2*I_NAI_F2xy_G2xyz_b;
  abcd[1669] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_F2xz_G2xyz_b;
  abcd[1670] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab;
  abcd[1671] = 4.0E0*I_NAI_Hx4y_G2xyz_ab-2.0E0*3*I_NAI_Fx2y_G2xyz_b;
  abcd[1672] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[1673] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Fx2z_G2xyz_b;
  abcd[1674] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab;
  abcd[1675] = 4.0E0*I_NAI_H5y_G2xyz_ab-2.0E0*4*I_NAI_F3y_G2xyz_b;
  abcd[1676] = 4.0E0*I_NAI_H4yz_G2xyz_ab-2.0E0*3*I_NAI_F2yz_G2xyz_b;
  abcd[1677] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*2*I_NAI_Fy2z_G2xyz_b;
  abcd[1678] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*1*I_NAI_F3z_G2xyz_b;
  abcd[1679] = 4.0E0*I_NAI_Hy4z_G2xyz_ab;
  abcd[1680] = 4.0E0*I_NAI_H4xy_G2x2z_ab-2.0E0*1*I_NAI_H4xy_D2x_a;
  abcd[1681] = 4.0E0*I_NAI_H3x2y_G2x2z_ab-2.0E0*1*I_NAI_H3x2y_D2x_a-2.0E0*1*I_NAI_F3x_G2x2z_b+1*I_NAI_F3x_D2x;
  abcd[1682] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2x_a;
  abcd[1683] = 4.0E0*I_NAI_H2x3y_G2x2z_ab-2.0E0*1*I_NAI_H2x3y_D2x_a-2.0E0*2*I_NAI_F2xy_G2x2z_b+2*1*I_NAI_F2xy_D2x;
  abcd[1684] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a-2.0E0*1*I_NAI_F2xz_G2x2z_b+1*I_NAI_F2xz_D2x;
  abcd[1685] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a;
  abcd[1686] = 4.0E0*I_NAI_Hx4y_G2x2z_ab-2.0E0*1*I_NAI_Hx4y_D2x_a-2.0E0*3*I_NAI_Fx2y_G2x2z_b+3*1*I_NAI_Fx2y_D2x;
  abcd[1687] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2x;
  abcd[1688] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2z_G2x2z_b+1*I_NAI_Fx2z_D2x;
  abcd[1689] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a;
  abcd[1690] = 4.0E0*I_NAI_H5y_G2x2z_ab-2.0E0*1*I_NAI_H5y_D2x_a-2.0E0*4*I_NAI_F3y_G2x2z_b+4*1*I_NAI_F3y_D2x;
  abcd[1691] = 4.0E0*I_NAI_H4yz_G2x2z_ab-2.0E0*1*I_NAI_H4yz_D2x_a-2.0E0*3*I_NAI_F2yz_G2x2z_b+3*1*I_NAI_F2yz_D2x;
  abcd[1692] = 4.0E0*I_NAI_H3y2z_G2x2z_ab-2.0E0*1*I_NAI_H3y2z_D2x_a-2.0E0*2*I_NAI_Fy2z_G2x2z_b+2*1*I_NAI_Fy2z_D2x;
  abcd[1693] = 4.0E0*I_NAI_H2y3z_G2x2z_ab-2.0E0*1*I_NAI_H2y3z_D2x_a-2.0E0*1*I_NAI_F3z_G2x2z_b+1*I_NAI_F3z_D2x;
  abcd[1694] = 4.0E0*I_NAI_Hy4z_G2x2z_ab-2.0E0*1*I_NAI_Hy4z_D2x_a;
  abcd[1695] = 4.0E0*I_NAI_H4xy_Gx2yz_ab;
  abcd[1696] = 4.0E0*I_NAI_H3x2y_Gx2yz_ab-2.0E0*1*I_NAI_F3x_Gx2yz_b;
  abcd[1697] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab;
  abcd[1698] = 4.0E0*I_NAI_H2x3y_Gx2yz_ab-2.0E0*2*I_NAI_F2xy_Gx2yz_b;
  abcd[1699] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_F2xz_Gx2yz_b;
  abcd[1700] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab;
  abcd[1701] = 4.0E0*I_NAI_Hx4y_Gx2yz_ab-2.0E0*3*I_NAI_Fx2y_Gx2yz_b;
  abcd[1702] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[1703] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Fx2z_Gx2yz_b;
  abcd[1704] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab;
  abcd[1705] = 4.0E0*I_NAI_H5y_Gx2yz_ab-2.0E0*4*I_NAI_F3y_Gx2yz_b;
  abcd[1706] = 4.0E0*I_NAI_H4yz_Gx2yz_ab-2.0E0*3*I_NAI_F2yz_Gx2yz_b;
  abcd[1707] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*2*I_NAI_Fy2z_Gx2yz_b;
  abcd[1708] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*1*I_NAI_F3z_Gx2yz_b;
  abcd[1709] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab;
  abcd[1710] = 4.0E0*I_NAI_H4xy_Gxy2z_ab-2.0E0*1*I_NAI_H4xy_Dxy_a;
  abcd[1711] = 4.0E0*I_NAI_H3x2y_Gxy2z_ab-2.0E0*1*I_NAI_H3x2y_Dxy_a-2.0E0*1*I_NAI_F3x_Gxy2z_b+1*I_NAI_F3x_Dxy;
  abcd[1712] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab-2.0E0*1*I_NAI_H3xyz_Dxy_a;
  abcd[1713] = 4.0E0*I_NAI_H2x3y_Gxy2z_ab-2.0E0*1*I_NAI_H2x3y_Dxy_a-2.0E0*2*I_NAI_F2xy_Gxy2z_b+2*1*I_NAI_F2xy_Dxy;
  abcd[1714] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_H2x2yz_Dxy_a-2.0E0*1*I_NAI_F2xz_Gxy2z_b+1*I_NAI_F2xz_Dxy;
  abcd[1715] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_H2xy2z_Dxy_a;
  abcd[1716] = 4.0E0*I_NAI_Hx4y_Gxy2z_ab-2.0E0*1*I_NAI_Hx4y_Dxy_a-2.0E0*3*I_NAI_Fx2y_Gxy2z_b+3*1*I_NAI_Fx2y_Dxy;
  abcd[1717] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab-2.0E0*1*I_NAI_Hx3yz_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gxy2z_b+2*1*I_NAI_Fxyz_Dxy;
  abcd[1718] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2z_Gxy2z_b+1*I_NAI_Fx2z_Dxy;
  abcd[1719] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab-2.0E0*1*I_NAI_Hxy3z_Dxy_a;
  abcd[1720] = 4.0E0*I_NAI_H5y_Gxy2z_ab-2.0E0*1*I_NAI_H5y_Dxy_a-2.0E0*4*I_NAI_F3y_Gxy2z_b+4*1*I_NAI_F3y_Dxy;
  abcd[1721] = 4.0E0*I_NAI_H4yz_Gxy2z_ab-2.0E0*1*I_NAI_H4yz_Dxy_a-2.0E0*3*I_NAI_F2yz_Gxy2z_b+3*1*I_NAI_F2yz_Dxy;
  abcd[1722] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*1*I_NAI_H3y2z_Dxy_a-2.0E0*2*I_NAI_Fy2z_Gxy2z_b+2*1*I_NAI_Fy2z_Dxy;
  abcd[1723] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*1*I_NAI_H2y3z_Dxy_a-2.0E0*1*I_NAI_F3z_Gxy2z_b+1*I_NAI_F3z_Dxy;
  abcd[1724] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab-2.0E0*1*I_NAI_Hy4z_Dxy_a;
  abcd[1725] = 4.0E0*I_NAI_H4xy_Gx3z_ab-2.0E0*2*I_NAI_H4xy_Dxz_a;
  abcd[1726] = 4.0E0*I_NAI_H3x2y_Gx3z_ab-2.0E0*2*I_NAI_H3x2y_Dxz_a-2.0E0*1*I_NAI_F3x_Gx3z_b+2*I_NAI_F3x_Dxz;
  abcd[1727] = 4.0E0*I_NAI_H3xyz_Gx3z_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a;
  abcd[1728] = 4.0E0*I_NAI_H2x3y_Gx3z_ab-2.0E0*2*I_NAI_H2x3y_Dxz_a-2.0E0*2*I_NAI_F2xy_Gx3z_b+2*2*I_NAI_F2xy_Dxz;
  abcd[1729] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a-2.0E0*1*I_NAI_F2xz_Gx3z_b+2*I_NAI_F2xz_Dxz;
  abcd[1730] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a;
  abcd[1731] = 4.0E0*I_NAI_Hx4y_Gx3z_ab-2.0E0*2*I_NAI_Hx4y_Dxz_a-2.0E0*3*I_NAI_Fx2y_Gx3z_b+3*2*I_NAI_Fx2y_Dxz;
  abcd[1732] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx3z_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[1733] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2z_Gx3z_b+2*I_NAI_Fx2z_Dxz;
  abcd[1734] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a;
  abcd[1735] = 4.0E0*I_NAI_H5y_Gx3z_ab-2.0E0*2*I_NAI_H5y_Dxz_a-2.0E0*4*I_NAI_F3y_Gx3z_b+4*2*I_NAI_F3y_Dxz;
  abcd[1736] = 4.0E0*I_NAI_H4yz_Gx3z_ab-2.0E0*2*I_NAI_H4yz_Dxz_a-2.0E0*3*I_NAI_F2yz_Gx3z_b+3*2*I_NAI_F2yz_Dxz;
  abcd[1737] = 4.0E0*I_NAI_H3y2z_Gx3z_ab-2.0E0*2*I_NAI_H3y2z_Dxz_a-2.0E0*2*I_NAI_Fy2z_Gx3z_b+2*2*I_NAI_Fy2z_Dxz;
  abcd[1738] = 4.0E0*I_NAI_H2y3z_Gx3z_ab-2.0E0*2*I_NAI_H2y3z_Dxz_a-2.0E0*1*I_NAI_F3z_Gx3z_b+2*I_NAI_F3z_Dxz;
  abcd[1739] = 4.0E0*I_NAI_Hy4z_Gx3z_ab-2.0E0*2*I_NAI_Hy4z_Dxz_a;
  abcd[1740] = 4.0E0*I_NAI_H4xy_G3yz_ab;
  abcd[1741] = 4.0E0*I_NAI_H3x2y_G3yz_ab-2.0E0*1*I_NAI_F3x_G3yz_b;
  abcd[1742] = 4.0E0*I_NAI_H3xyz_G3yz_ab;
  abcd[1743] = 4.0E0*I_NAI_H2x3y_G3yz_ab-2.0E0*2*I_NAI_F2xy_G3yz_b;
  abcd[1744] = 4.0E0*I_NAI_H2x2yz_G3yz_ab-2.0E0*1*I_NAI_F2xz_G3yz_b;
  abcd[1745] = 4.0E0*I_NAI_H2xy2z_G3yz_ab;
  abcd[1746] = 4.0E0*I_NAI_Hx4y_G3yz_ab-2.0E0*3*I_NAI_Fx2y_G3yz_b;
  abcd[1747] = 4.0E0*I_NAI_Hx3yz_G3yz_ab-2.0E0*2*I_NAI_Fxyz_G3yz_b;
  abcd[1748] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab-2.0E0*1*I_NAI_Fx2z_G3yz_b;
  abcd[1749] = 4.0E0*I_NAI_Hxy3z_G3yz_ab;
  abcd[1750] = 4.0E0*I_NAI_H5y_G3yz_ab-2.0E0*4*I_NAI_F3y_G3yz_b;
  abcd[1751] = 4.0E0*I_NAI_H4yz_G3yz_ab-2.0E0*3*I_NAI_F2yz_G3yz_b;
  abcd[1752] = 4.0E0*I_NAI_H3y2z_G3yz_ab-2.0E0*2*I_NAI_Fy2z_G3yz_b;
  abcd[1753] = 4.0E0*I_NAI_H2y3z_G3yz_ab-2.0E0*1*I_NAI_F3z_G3yz_b;
  abcd[1754] = 4.0E0*I_NAI_Hy4z_G3yz_ab;
  abcd[1755] = 4.0E0*I_NAI_H4xy_G2y2z_ab-2.0E0*1*I_NAI_H4xy_D2y_a;
  abcd[1756] = 4.0E0*I_NAI_H3x2y_G2y2z_ab-2.0E0*1*I_NAI_H3x2y_D2y_a-2.0E0*1*I_NAI_F3x_G2y2z_b+1*I_NAI_F3x_D2y;
  abcd[1757] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2y_a;
  abcd[1758] = 4.0E0*I_NAI_H2x3y_G2y2z_ab-2.0E0*1*I_NAI_H2x3y_D2y_a-2.0E0*2*I_NAI_F2xy_G2y2z_b+2*1*I_NAI_F2xy_D2y;
  abcd[1759] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a-2.0E0*1*I_NAI_F2xz_G2y2z_b+1*I_NAI_F2xz_D2y;
  abcd[1760] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a;
  abcd[1761] = 4.0E0*I_NAI_Hx4y_G2y2z_ab-2.0E0*1*I_NAI_Hx4y_D2y_a-2.0E0*3*I_NAI_Fx2y_G2y2z_b+3*1*I_NAI_Fx2y_D2y;
  abcd[1762] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2y;
  abcd[1763] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2z_G2y2z_b+1*I_NAI_Fx2z_D2y;
  abcd[1764] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a;
  abcd[1765] = 4.0E0*I_NAI_H5y_G2y2z_ab-2.0E0*1*I_NAI_H5y_D2y_a-2.0E0*4*I_NAI_F3y_G2y2z_b+4*1*I_NAI_F3y_D2y;
  abcd[1766] = 4.0E0*I_NAI_H4yz_G2y2z_ab-2.0E0*1*I_NAI_H4yz_D2y_a-2.0E0*3*I_NAI_F2yz_G2y2z_b+3*1*I_NAI_F2yz_D2y;
  abcd[1767] = 4.0E0*I_NAI_H3y2z_G2y2z_ab-2.0E0*1*I_NAI_H3y2z_D2y_a-2.0E0*2*I_NAI_Fy2z_G2y2z_b+2*1*I_NAI_Fy2z_D2y;
  abcd[1768] = 4.0E0*I_NAI_H2y3z_G2y2z_ab-2.0E0*1*I_NAI_H2y3z_D2y_a-2.0E0*1*I_NAI_F3z_G2y2z_b+1*I_NAI_F3z_D2y;
  abcd[1769] = 4.0E0*I_NAI_Hy4z_G2y2z_ab-2.0E0*1*I_NAI_Hy4z_D2y_a;
  abcd[1770] = 4.0E0*I_NAI_H4xy_Gy3z_ab-2.0E0*2*I_NAI_H4xy_Dyz_a;
  abcd[1771] = 4.0E0*I_NAI_H3x2y_Gy3z_ab-2.0E0*2*I_NAI_H3x2y_Dyz_a-2.0E0*1*I_NAI_F3x_Gy3z_b+2*I_NAI_F3x_Dyz;
  abcd[1772] = 4.0E0*I_NAI_H3xyz_Gy3z_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a;
  abcd[1773] = 4.0E0*I_NAI_H2x3y_Gy3z_ab-2.0E0*2*I_NAI_H2x3y_Dyz_a-2.0E0*2*I_NAI_F2xy_Gy3z_b+2*2*I_NAI_F2xy_Dyz;
  abcd[1774] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a-2.0E0*1*I_NAI_F2xz_Gy3z_b+2*I_NAI_F2xz_Dyz;
  abcd[1775] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a;
  abcd[1776] = 4.0E0*I_NAI_Hx4y_Gy3z_ab-2.0E0*2*I_NAI_Hx4y_Dyz_a-2.0E0*3*I_NAI_Fx2y_Gy3z_b+3*2*I_NAI_Fx2y_Dyz;
  abcd[1777] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a-2.0E0*2*I_NAI_Fxyz_Gy3z_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[1778] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2z_Gy3z_b+2*I_NAI_Fx2z_Dyz;
  abcd[1779] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a;
  abcd[1780] = 4.0E0*I_NAI_H5y_Gy3z_ab-2.0E0*2*I_NAI_H5y_Dyz_a-2.0E0*4*I_NAI_F3y_Gy3z_b+4*2*I_NAI_F3y_Dyz;
  abcd[1781] = 4.0E0*I_NAI_H4yz_Gy3z_ab-2.0E0*2*I_NAI_H4yz_Dyz_a-2.0E0*3*I_NAI_F2yz_Gy3z_b+3*2*I_NAI_F2yz_Dyz;
  abcd[1782] = 4.0E0*I_NAI_H3y2z_Gy3z_ab-2.0E0*2*I_NAI_H3y2z_Dyz_a-2.0E0*2*I_NAI_Fy2z_Gy3z_b+2*2*I_NAI_Fy2z_Dyz;
  abcd[1783] = 4.0E0*I_NAI_H2y3z_Gy3z_ab-2.0E0*2*I_NAI_H2y3z_Dyz_a-2.0E0*1*I_NAI_F3z_Gy3z_b+2*I_NAI_F3z_Dyz;
  abcd[1784] = 4.0E0*I_NAI_Hy4z_Gy3z_ab-2.0E0*2*I_NAI_Hy4z_Dyz_a;
  abcd[1785] = 4.0E0*I_NAI_H4xy_G4z_ab-2.0E0*3*I_NAI_H4xy_D2z_a;
  abcd[1786] = 4.0E0*I_NAI_H3x2y_G4z_ab-2.0E0*3*I_NAI_H3x2y_D2z_a-2.0E0*1*I_NAI_F3x_G4z_b+3*I_NAI_F3x_D2z;
  abcd[1787] = 4.0E0*I_NAI_H3xyz_G4z_ab-2.0E0*3*I_NAI_H3xyz_D2z_a;
  abcd[1788] = 4.0E0*I_NAI_H2x3y_G4z_ab-2.0E0*3*I_NAI_H2x3y_D2z_a-2.0E0*2*I_NAI_F2xy_G4z_b+2*3*I_NAI_F2xy_D2z;
  abcd[1789] = 4.0E0*I_NAI_H2x2yz_G4z_ab-2.0E0*3*I_NAI_H2x2yz_D2z_a-2.0E0*1*I_NAI_F2xz_G4z_b+3*I_NAI_F2xz_D2z;
  abcd[1790] = 4.0E0*I_NAI_H2xy2z_G4z_ab-2.0E0*3*I_NAI_H2xy2z_D2z_a;
  abcd[1791] = 4.0E0*I_NAI_Hx4y_G4z_ab-2.0E0*3*I_NAI_Hx4y_D2z_a-2.0E0*3*I_NAI_Fx2y_G4z_b+3*3*I_NAI_Fx2y_D2z;
  abcd[1792] = 4.0E0*I_NAI_Hx3yz_G4z_ab-2.0E0*3*I_NAI_Hx3yz_D2z_a-2.0E0*2*I_NAI_Fxyz_G4z_b+2*3*I_NAI_Fxyz_D2z;
  abcd[1793] = 4.0E0*I_NAI_Hx2y2z_G4z_ab-2.0E0*3*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2z_G4z_b+3*I_NAI_Fx2z_D2z;
  abcd[1794] = 4.0E0*I_NAI_Hxy3z_G4z_ab-2.0E0*3*I_NAI_Hxy3z_D2z_a;
  abcd[1795] = 4.0E0*I_NAI_H5y_G4z_ab-2.0E0*3*I_NAI_H5y_D2z_a-2.0E0*4*I_NAI_F3y_G4z_b+4*3*I_NAI_F3y_D2z;
  abcd[1796] = 4.0E0*I_NAI_H4yz_G4z_ab-2.0E0*3*I_NAI_H4yz_D2z_a-2.0E0*3*I_NAI_F2yz_G4z_b+3*3*I_NAI_F2yz_D2z;
  abcd[1797] = 4.0E0*I_NAI_H3y2z_G4z_ab-2.0E0*3*I_NAI_H3y2z_D2z_a-2.0E0*2*I_NAI_Fy2z_G4z_b+2*3*I_NAI_Fy2z_D2z;
  abcd[1798] = 4.0E0*I_NAI_H2y3z_G4z_ab-2.0E0*3*I_NAI_H2y3z_D2z_a-2.0E0*1*I_NAI_F3z_G4z_b+3*I_NAI_F3z_D2z;
  abcd[1799] = 4.0E0*I_NAI_Hy4z_G4z_ab-2.0E0*3*I_NAI_Hy4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1800] = 4.0E0*I_NAI_H4xz_G4x_ab-2.0E0*3*I_NAI_H4xz_D2x_a;
  abcd[1801] = 4.0E0*I_NAI_H3xyz_G4x_ab-2.0E0*3*I_NAI_H3xyz_D2x_a;
  abcd[1802] = 4.0E0*I_NAI_H3x2z_G4x_ab-2.0E0*3*I_NAI_H3x2z_D2x_a-2.0E0*1*I_NAI_F3x_G4x_b+3*I_NAI_F3x_D2x;
  abcd[1803] = 4.0E0*I_NAI_H2x2yz_G4x_ab-2.0E0*3*I_NAI_H2x2yz_D2x_a;
  abcd[1804] = 4.0E0*I_NAI_H2xy2z_G4x_ab-2.0E0*3*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_F2xy_G4x_b+3*I_NAI_F2xy_D2x;
  abcd[1805] = 4.0E0*I_NAI_H2x3z_G4x_ab-2.0E0*3*I_NAI_H2x3z_D2x_a-2.0E0*2*I_NAI_F2xz_G4x_b+2*3*I_NAI_F2xz_D2x;
  abcd[1806] = 4.0E0*I_NAI_Hx3yz_G4x_ab-2.0E0*3*I_NAI_Hx3yz_D2x_a;
  abcd[1807] = 4.0E0*I_NAI_Hx2y2z_G4x_ab-2.0E0*3*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2y_G4x_b+3*I_NAI_Fx2y_D2x;
  abcd[1808] = 4.0E0*I_NAI_Hxy3z_G4x_ab-2.0E0*3*I_NAI_Hxy3z_D2x_a-2.0E0*2*I_NAI_Fxyz_G4x_b+2*3*I_NAI_Fxyz_D2x;
  abcd[1809] = 4.0E0*I_NAI_Hx4z_G4x_ab-2.0E0*3*I_NAI_Hx4z_D2x_a-2.0E0*3*I_NAI_Fx2z_G4x_b+3*3*I_NAI_Fx2z_D2x;
  abcd[1810] = 4.0E0*I_NAI_H4yz_G4x_ab-2.0E0*3*I_NAI_H4yz_D2x_a;
  abcd[1811] = 4.0E0*I_NAI_H3y2z_G4x_ab-2.0E0*3*I_NAI_H3y2z_D2x_a-2.0E0*1*I_NAI_F3y_G4x_b+3*I_NAI_F3y_D2x;
  abcd[1812] = 4.0E0*I_NAI_H2y3z_G4x_ab-2.0E0*3*I_NAI_H2y3z_D2x_a-2.0E0*2*I_NAI_F2yz_G4x_b+2*3*I_NAI_F2yz_D2x;
  abcd[1813] = 4.0E0*I_NAI_Hy4z_G4x_ab-2.0E0*3*I_NAI_Hy4z_D2x_a-2.0E0*3*I_NAI_Fy2z_G4x_b+3*3*I_NAI_Fy2z_D2x;
  abcd[1814] = 4.0E0*I_NAI_H5z_G4x_ab-2.0E0*3*I_NAI_H5z_D2x_a-2.0E0*4*I_NAI_F3z_G4x_b+4*3*I_NAI_F3z_D2x;
  abcd[1815] = 4.0E0*I_NAI_H4xz_G3xy_ab-2.0E0*2*I_NAI_H4xz_Dxy_a;
  abcd[1816] = 4.0E0*I_NAI_H3xyz_G3xy_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a;
  abcd[1817] = 4.0E0*I_NAI_H3x2z_G3xy_ab-2.0E0*2*I_NAI_H3x2z_Dxy_a-2.0E0*1*I_NAI_F3x_G3xy_b+2*I_NAI_F3x_Dxy;
  abcd[1818] = 4.0E0*I_NAI_H2x2yz_G3xy_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a;
  abcd[1819] = 4.0E0*I_NAI_H2xy2z_G3xy_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_F2xy_G3xy_b+2*I_NAI_F2xy_Dxy;
  abcd[1820] = 4.0E0*I_NAI_H2x3z_G3xy_ab-2.0E0*2*I_NAI_H2x3z_Dxy_a-2.0E0*2*I_NAI_F2xz_G3xy_b+2*2*I_NAI_F2xz_Dxy;
  abcd[1821] = 4.0E0*I_NAI_Hx3yz_G3xy_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a;
  abcd[1822] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2y_G3xy_b+2*I_NAI_Fx2y_Dxy;
  abcd[1823] = 4.0E0*I_NAI_Hxy3z_G3xy_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a-2.0E0*2*I_NAI_Fxyz_G3xy_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[1824] = 4.0E0*I_NAI_Hx4z_G3xy_ab-2.0E0*2*I_NAI_Hx4z_Dxy_a-2.0E0*3*I_NAI_Fx2z_G3xy_b+3*2*I_NAI_Fx2z_Dxy;
  abcd[1825] = 4.0E0*I_NAI_H4yz_G3xy_ab-2.0E0*2*I_NAI_H4yz_Dxy_a;
  abcd[1826] = 4.0E0*I_NAI_H3y2z_G3xy_ab-2.0E0*2*I_NAI_H3y2z_Dxy_a-2.0E0*1*I_NAI_F3y_G3xy_b+2*I_NAI_F3y_Dxy;
  abcd[1827] = 4.0E0*I_NAI_H2y3z_G3xy_ab-2.0E0*2*I_NAI_H2y3z_Dxy_a-2.0E0*2*I_NAI_F2yz_G3xy_b+2*2*I_NAI_F2yz_Dxy;
  abcd[1828] = 4.0E0*I_NAI_Hy4z_G3xy_ab-2.0E0*2*I_NAI_Hy4z_Dxy_a-2.0E0*3*I_NAI_Fy2z_G3xy_b+3*2*I_NAI_Fy2z_Dxy;
  abcd[1829] = 4.0E0*I_NAI_H5z_G3xy_ab-2.0E0*2*I_NAI_H5z_Dxy_a-2.0E0*4*I_NAI_F3z_G3xy_b+4*2*I_NAI_F3z_Dxy;
  abcd[1830] = 4.0E0*I_NAI_H4xz_G3xz_ab-2.0E0*2*I_NAI_H4xz_Dxz_a;
  abcd[1831] = 4.0E0*I_NAI_H3xyz_G3xz_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a;
  abcd[1832] = 4.0E0*I_NAI_H3x2z_G3xz_ab-2.0E0*2*I_NAI_H3x2z_Dxz_a-2.0E0*1*I_NAI_F3x_G3xz_b+2*I_NAI_F3x_Dxz;
  abcd[1833] = 4.0E0*I_NAI_H2x2yz_G3xz_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a;
  abcd[1834] = 4.0E0*I_NAI_H2xy2z_G3xz_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_F2xy_G3xz_b+2*I_NAI_F2xy_Dxz;
  abcd[1835] = 4.0E0*I_NAI_H2x3z_G3xz_ab-2.0E0*2*I_NAI_H2x3z_Dxz_a-2.0E0*2*I_NAI_F2xz_G3xz_b+2*2*I_NAI_F2xz_Dxz;
  abcd[1836] = 4.0E0*I_NAI_Hx3yz_G3xz_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a;
  abcd[1837] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2y_G3xz_b+2*I_NAI_Fx2y_Dxz;
  abcd[1838] = 4.0E0*I_NAI_Hxy3z_G3xz_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a-2.0E0*2*I_NAI_Fxyz_G3xz_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[1839] = 4.0E0*I_NAI_Hx4z_G3xz_ab-2.0E0*2*I_NAI_Hx4z_Dxz_a-2.0E0*3*I_NAI_Fx2z_G3xz_b+3*2*I_NAI_Fx2z_Dxz;
  abcd[1840] = 4.0E0*I_NAI_H4yz_G3xz_ab-2.0E0*2*I_NAI_H4yz_Dxz_a;
  abcd[1841] = 4.0E0*I_NAI_H3y2z_G3xz_ab-2.0E0*2*I_NAI_H3y2z_Dxz_a-2.0E0*1*I_NAI_F3y_G3xz_b+2*I_NAI_F3y_Dxz;
  abcd[1842] = 4.0E0*I_NAI_H2y3z_G3xz_ab-2.0E0*2*I_NAI_H2y3z_Dxz_a-2.0E0*2*I_NAI_F2yz_G3xz_b+2*2*I_NAI_F2yz_Dxz;
  abcd[1843] = 4.0E0*I_NAI_Hy4z_G3xz_ab-2.0E0*2*I_NAI_Hy4z_Dxz_a-2.0E0*3*I_NAI_Fy2z_G3xz_b+3*2*I_NAI_Fy2z_Dxz;
  abcd[1844] = 4.0E0*I_NAI_H5z_G3xz_ab-2.0E0*2*I_NAI_H5z_Dxz_a-2.0E0*4*I_NAI_F3z_G3xz_b+4*2*I_NAI_F3z_Dxz;
  abcd[1845] = 4.0E0*I_NAI_H4xz_G2x2y_ab-2.0E0*1*I_NAI_H4xz_D2y_a;
  abcd[1846] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2y_a;
  abcd[1847] = 4.0E0*I_NAI_H3x2z_G2x2y_ab-2.0E0*1*I_NAI_H3x2z_D2y_a-2.0E0*1*I_NAI_F3x_G2x2y_b+1*I_NAI_F3x_D2y;
  abcd[1848] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a;
  abcd[1849] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_F2xy_G2x2y_b+1*I_NAI_F2xy_D2y;
  abcd[1850] = 4.0E0*I_NAI_H2x3z_G2x2y_ab-2.0E0*1*I_NAI_H2x3z_D2y_a-2.0E0*2*I_NAI_F2xz_G2x2y_b+2*1*I_NAI_F2xz_D2y;
  abcd[1851] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a;
  abcd[1852] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2y_G2x2y_b+1*I_NAI_Fx2y_D2y;
  abcd[1853] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2y;
  abcd[1854] = 4.0E0*I_NAI_Hx4z_G2x2y_ab-2.0E0*1*I_NAI_Hx4z_D2y_a-2.0E0*3*I_NAI_Fx2z_G2x2y_b+3*1*I_NAI_Fx2z_D2y;
  abcd[1855] = 4.0E0*I_NAI_H4yz_G2x2y_ab-2.0E0*1*I_NAI_H4yz_D2y_a;
  abcd[1856] = 4.0E0*I_NAI_H3y2z_G2x2y_ab-2.0E0*1*I_NAI_H3y2z_D2y_a-2.0E0*1*I_NAI_F3y_G2x2y_b+1*I_NAI_F3y_D2y;
  abcd[1857] = 4.0E0*I_NAI_H2y3z_G2x2y_ab-2.0E0*1*I_NAI_H2y3z_D2y_a-2.0E0*2*I_NAI_F2yz_G2x2y_b+2*1*I_NAI_F2yz_D2y;
  abcd[1858] = 4.0E0*I_NAI_Hy4z_G2x2y_ab-2.0E0*1*I_NAI_Hy4z_D2y_a-2.0E0*3*I_NAI_Fy2z_G2x2y_b+3*1*I_NAI_Fy2z_D2y;
  abcd[1859] = 4.0E0*I_NAI_H5z_G2x2y_ab-2.0E0*1*I_NAI_H5z_D2y_a-2.0E0*4*I_NAI_F3z_G2x2y_b+4*1*I_NAI_F3z_D2y;
  abcd[1860] = 4.0E0*I_NAI_H4xz_G2xyz_ab-2.0E0*1*I_NAI_H4xz_Dyz_a;
  abcd[1861] = 4.0E0*I_NAI_H3xyz_G2xyz_ab-2.0E0*1*I_NAI_H3xyz_Dyz_a;
  abcd[1862] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*1*I_NAI_H3x2z_Dyz_a-2.0E0*1*I_NAI_F3x_G2xyz_b+1*I_NAI_F3x_Dyz;
  abcd[1863] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab-2.0E0*1*I_NAI_H2x2yz_Dyz_a;
  abcd[1864] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_F2xy_G2xyz_b+1*I_NAI_F2xy_Dyz;
  abcd[1865] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*1*I_NAI_H2x3z_Dyz_a-2.0E0*2*I_NAI_F2xz_G2xyz_b+2*1*I_NAI_F2xz_Dyz;
  abcd[1866] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab-2.0E0*1*I_NAI_Hx3yz_Dyz_a;
  abcd[1867] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2y_G2xyz_b+1*I_NAI_Fx2y_Dyz;
  abcd[1868] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab-2.0E0*1*I_NAI_Hxy3z_Dyz_a-2.0E0*2*I_NAI_Fxyz_G2xyz_b+2*1*I_NAI_Fxyz_Dyz;
  abcd[1869] = 4.0E0*I_NAI_Hx4z_G2xyz_ab-2.0E0*1*I_NAI_Hx4z_Dyz_a-2.0E0*3*I_NAI_Fx2z_G2xyz_b+3*1*I_NAI_Fx2z_Dyz;
  abcd[1870] = 4.0E0*I_NAI_H4yz_G2xyz_ab-2.0E0*1*I_NAI_H4yz_Dyz_a;
  abcd[1871] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*1*I_NAI_H3y2z_Dyz_a-2.0E0*1*I_NAI_F3y_G2xyz_b+1*I_NAI_F3y_Dyz;
  abcd[1872] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*1*I_NAI_H2y3z_Dyz_a-2.0E0*2*I_NAI_F2yz_G2xyz_b+2*1*I_NAI_F2yz_Dyz;
  abcd[1873] = 4.0E0*I_NAI_Hy4z_G2xyz_ab-2.0E0*1*I_NAI_Hy4z_Dyz_a-2.0E0*3*I_NAI_Fy2z_G2xyz_b+3*1*I_NAI_Fy2z_Dyz;
  abcd[1874] = 4.0E0*I_NAI_H5z_G2xyz_ab-2.0E0*1*I_NAI_H5z_Dyz_a-2.0E0*4*I_NAI_F3z_G2xyz_b+4*1*I_NAI_F3z_Dyz;
  abcd[1875] = 4.0E0*I_NAI_H4xz_G2x2z_ab-2.0E0*1*I_NAI_H4xz_D2z_a;
  abcd[1876] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a;
  abcd[1877] = 4.0E0*I_NAI_H3x2z_G2x2z_ab-2.0E0*1*I_NAI_H3x2z_D2z_a-2.0E0*1*I_NAI_F3x_G2x2z_b+1*I_NAI_F3x_D2z;
  abcd[1878] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a;
  abcd[1879] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_F2xy_G2x2z_b+1*I_NAI_F2xy_D2z;
  abcd[1880] = 4.0E0*I_NAI_H2x3z_G2x2z_ab-2.0E0*1*I_NAI_H2x3z_D2z_a-2.0E0*2*I_NAI_F2xz_G2x2z_b+2*1*I_NAI_F2xz_D2z;
  abcd[1881] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a;
  abcd[1882] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2y_G2x2z_b+1*I_NAI_Fx2y_D2z;
  abcd[1883] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[1884] = 4.0E0*I_NAI_Hx4z_G2x2z_ab-2.0E0*1*I_NAI_Hx4z_D2z_a-2.0E0*3*I_NAI_Fx2z_G2x2z_b+3*1*I_NAI_Fx2z_D2z;
  abcd[1885] = 4.0E0*I_NAI_H4yz_G2x2z_ab-2.0E0*1*I_NAI_H4yz_D2z_a;
  abcd[1886] = 4.0E0*I_NAI_H3y2z_G2x2z_ab-2.0E0*1*I_NAI_H3y2z_D2z_a-2.0E0*1*I_NAI_F3y_G2x2z_b+1*I_NAI_F3y_D2z;
  abcd[1887] = 4.0E0*I_NAI_H2y3z_G2x2z_ab-2.0E0*1*I_NAI_H2y3z_D2z_a-2.0E0*2*I_NAI_F2yz_G2x2z_b+2*1*I_NAI_F2yz_D2z;
  abcd[1888] = 4.0E0*I_NAI_Hy4z_G2x2z_ab-2.0E0*1*I_NAI_Hy4z_D2z_a-2.0E0*3*I_NAI_Fy2z_G2x2z_b+3*1*I_NAI_Fy2z_D2z;
  abcd[1889] = 4.0E0*I_NAI_H5z_G2x2z_ab-2.0E0*1*I_NAI_H5z_D2z_a-2.0E0*4*I_NAI_F3z_G2x2z_b+4*1*I_NAI_F3z_D2z;
  abcd[1890] = 4.0E0*I_NAI_H4xz_Gx3y_ab;
  abcd[1891] = 4.0E0*I_NAI_H3xyz_Gx3y_ab;
  abcd[1892] = 4.0E0*I_NAI_H3x2z_Gx3y_ab-2.0E0*1*I_NAI_F3x_Gx3y_b;
  abcd[1893] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab;
  abcd[1894] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab-2.0E0*1*I_NAI_F2xy_Gx3y_b;
  abcd[1895] = 4.0E0*I_NAI_H2x3z_Gx3y_ab-2.0E0*2*I_NAI_F2xz_Gx3y_b;
  abcd[1896] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab;
  abcd[1897] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab-2.0E0*1*I_NAI_Fx2y_Gx3y_b;
  abcd[1898] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab-2.0E0*2*I_NAI_Fxyz_Gx3y_b;
  abcd[1899] = 4.0E0*I_NAI_Hx4z_Gx3y_ab-2.0E0*3*I_NAI_Fx2z_Gx3y_b;
  abcd[1900] = 4.0E0*I_NAI_H4yz_Gx3y_ab;
  abcd[1901] = 4.0E0*I_NAI_H3y2z_Gx3y_ab-2.0E0*1*I_NAI_F3y_Gx3y_b;
  abcd[1902] = 4.0E0*I_NAI_H2y3z_Gx3y_ab-2.0E0*2*I_NAI_F2yz_Gx3y_b;
  abcd[1903] = 4.0E0*I_NAI_Hy4z_Gx3y_ab-2.0E0*3*I_NAI_Fy2z_Gx3y_b;
  abcd[1904] = 4.0E0*I_NAI_H5z_Gx3y_ab-2.0E0*4*I_NAI_F3z_Gx3y_b;
  abcd[1905] = 4.0E0*I_NAI_H4xz_Gx2yz_ab;
  abcd[1906] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab;
  abcd[1907] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*1*I_NAI_F3x_Gx2yz_b;
  abcd[1908] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab;
  abcd[1909] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_F2xy_Gx2yz_b;
  abcd[1910] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*2*I_NAI_F2xz_Gx2yz_b;
  abcd[1911] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab;
  abcd[1912] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Fx2y_Gx2yz_b;
  abcd[1913] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[1914] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab-2.0E0*3*I_NAI_Fx2z_Gx2yz_b;
  abcd[1915] = 4.0E0*I_NAI_H4yz_Gx2yz_ab;
  abcd[1916] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*1*I_NAI_F3y_Gx2yz_b;
  abcd[1917] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*2*I_NAI_F2yz_Gx2yz_b;
  abcd[1918] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab-2.0E0*3*I_NAI_Fy2z_Gx2yz_b;
  abcd[1919] = 4.0E0*I_NAI_H5z_Gx2yz_ab-2.0E0*4*I_NAI_F3z_Gx2yz_b;
  abcd[1920] = 4.0E0*I_NAI_H4xz_Gxy2z_ab;
  abcd[1921] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab;
  abcd[1922] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*1*I_NAI_F3x_Gxy2z_b;
  abcd[1923] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab;
  abcd[1924] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_F2xy_Gxy2z_b;
  abcd[1925] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*2*I_NAI_F2xz_Gxy2z_b;
  abcd[1926] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab;
  abcd[1927] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Fx2y_Gxy2z_b;
  abcd[1928] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[1929] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab-2.0E0*3*I_NAI_Fx2z_Gxy2z_b;
  abcd[1930] = 4.0E0*I_NAI_H4yz_Gxy2z_ab;
  abcd[1931] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*1*I_NAI_F3y_Gxy2z_b;
  abcd[1932] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*2*I_NAI_F2yz_Gxy2z_b;
  abcd[1933] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab-2.0E0*3*I_NAI_Fy2z_Gxy2z_b;
  abcd[1934] = 4.0E0*I_NAI_H5z_Gxy2z_ab-2.0E0*4*I_NAI_F3z_Gxy2z_b;
  abcd[1935] = 4.0E0*I_NAI_H4xz_Gx3z_ab;
  abcd[1936] = 4.0E0*I_NAI_H3xyz_Gx3z_ab;
  abcd[1937] = 4.0E0*I_NAI_H3x2z_Gx3z_ab-2.0E0*1*I_NAI_F3x_Gx3z_b;
  abcd[1938] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab;
  abcd[1939] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab-2.0E0*1*I_NAI_F2xy_Gx3z_b;
  abcd[1940] = 4.0E0*I_NAI_H2x3z_Gx3z_ab-2.0E0*2*I_NAI_F2xz_Gx3z_b;
  abcd[1941] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab;
  abcd[1942] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab-2.0E0*1*I_NAI_Fx2y_Gx3z_b;
  abcd[1943] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab-2.0E0*2*I_NAI_Fxyz_Gx3z_b;
  abcd[1944] = 4.0E0*I_NAI_Hx4z_Gx3z_ab-2.0E0*3*I_NAI_Fx2z_Gx3z_b;
  abcd[1945] = 4.0E0*I_NAI_H4yz_Gx3z_ab;
  abcd[1946] = 4.0E0*I_NAI_H3y2z_Gx3z_ab-2.0E0*1*I_NAI_F3y_Gx3z_b;
  abcd[1947] = 4.0E0*I_NAI_H2y3z_Gx3z_ab-2.0E0*2*I_NAI_F2yz_Gx3z_b;
  abcd[1948] = 4.0E0*I_NAI_Hy4z_Gx3z_ab-2.0E0*3*I_NAI_Fy2z_Gx3z_b;
  abcd[1949] = 4.0E0*I_NAI_H5z_Gx3z_ab-2.0E0*4*I_NAI_F3z_Gx3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[1950] = 4.0E0*I_NAI_H4xz_G3xy_ab;
  abcd[1951] = 4.0E0*I_NAI_H3xyz_G3xy_ab;
  abcd[1952] = 4.0E0*I_NAI_H3x2z_G3xy_ab-2.0E0*1*I_NAI_F3x_G3xy_b;
  abcd[1953] = 4.0E0*I_NAI_H2x2yz_G3xy_ab;
  abcd[1954] = 4.0E0*I_NAI_H2xy2z_G3xy_ab-2.0E0*1*I_NAI_F2xy_G3xy_b;
  abcd[1955] = 4.0E0*I_NAI_H2x3z_G3xy_ab-2.0E0*2*I_NAI_F2xz_G3xy_b;
  abcd[1956] = 4.0E0*I_NAI_Hx3yz_G3xy_ab;
  abcd[1957] = 4.0E0*I_NAI_Hx2y2z_G3xy_ab-2.0E0*1*I_NAI_Fx2y_G3xy_b;
  abcd[1958] = 4.0E0*I_NAI_Hxy3z_G3xy_ab-2.0E0*2*I_NAI_Fxyz_G3xy_b;
  abcd[1959] = 4.0E0*I_NAI_Hx4z_G3xy_ab-2.0E0*3*I_NAI_Fx2z_G3xy_b;
  abcd[1960] = 4.0E0*I_NAI_H4yz_G3xy_ab;
  abcd[1961] = 4.0E0*I_NAI_H3y2z_G3xy_ab-2.0E0*1*I_NAI_F3y_G3xy_b;
  abcd[1962] = 4.0E0*I_NAI_H2y3z_G3xy_ab-2.0E0*2*I_NAI_F2yz_G3xy_b;
  abcd[1963] = 4.0E0*I_NAI_Hy4z_G3xy_ab-2.0E0*3*I_NAI_Fy2z_G3xy_b;
  abcd[1964] = 4.0E0*I_NAI_H5z_G3xy_ab-2.0E0*4*I_NAI_F3z_G3xy_b;
  abcd[1965] = 4.0E0*I_NAI_H4xz_G2x2y_ab-2.0E0*1*I_NAI_H4xz_D2x_a;
  abcd[1966] = 4.0E0*I_NAI_H3xyz_G2x2y_ab-2.0E0*1*I_NAI_H3xyz_D2x_a;
  abcd[1967] = 4.0E0*I_NAI_H3x2z_G2x2y_ab-2.0E0*1*I_NAI_H3x2z_D2x_a-2.0E0*1*I_NAI_F3x_G2x2y_b+1*I_NAI_F3x_D2x;
  abcd[1968] = 4.0E0*I_NAI_H2x2yz_G2x2y_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a;
  abcd[1969] = 4.0E0*I_NAI_H2xy2z_G2x2y_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_F2xy_G2x2y_b+1*I_NAI_F2xy_D2x;
  abcd[1970] = 4.0E0*I_NAI_H2x3z_G2x2y_ab-2.0E0*1*I_NAI_H2x3z_D2x_a-2.0E0*2*I_NAI_F2xz_G2x2y_b+2*1*I_NAI_F2xz_D2x;
  abcd[1971] = 4.0E0*I_NAI_Hx3yz_G2x2y_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a;
  abcd[1972] = 4.0E0*I_NAI_Hx2y2z_G2x2y_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2y_G2x2y_b+1*I_NAI_Fx2y_D2x;
  abcd[1973] = 4.0E0*I_NAI_Hxy3z_G2x2y_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2y_b+2*1*I_NAI_Fxyz_D2x;
  abcd[1974] = 4.0E0*I_NAI_Hx4z_G2x2y_ab-2.0E0*1*I_NAI_Hx4z_D2x_a-2.0E0*3*I_NAI_Fx2z_G2x2y_b+3*1*I_NAI_Fx2z_D2x;
  abcd[1975] = 4.0E0*I_NAI_H4yz_G2x2y_ab-2.0E0*1*I_NAI_H4yz_D2x_a;
  abcd[1976] = 4.0E0*I_NAI_H3y2z_G2x2y_ab-2.0E0*1*I_NAI_H3y2z_D2x_a-2.0E0*1*I_NAI_F3y_G2x2y_b+1*I_NAI_F3y_D2x;
  abcd[1977] = 4.0E0*I_NAI_H2y3z_G2x2y_ab-2.0E0*1*I_NAI_H2y3z_D2x_a-2.0E0*2*I_NAI_F2yz_G2x2y_b+2*1*I_NAI_F2yz_D2x;
  abcd[1978] = 4.0E0*I_NAI_Hy4z_G2x2y_ab-2.0E0*1*I_NAI_Hy4z_D2x_a-2.0E0*3*I_NAI_Fy2z_G2x2y_b+3*1*I_NAI_Fy2z_D2x;
  abcd[1979] = 4.0E0*I_NAI_H5z_G2x2y_ab-2.0E0*1*I_NAI_H5z_D2x_a-2.0E0*4*I_NAI_F3z_G2x2y_b+4*1*I_NAI_F3z_D2x;
  abcd[1980] = 4.0E0*I_NAI_H4xz_G2xyz_ab;
  abcd[1981] = 4.0E0*I_NAI_H3xyz_G2xyz_ab;
  abcd[1982] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*1*I_NAI_F3x_G2xyz_b;
  abcd[1983] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab;
  abcd[1984] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_F2xy_G2xyz_b;
  abcd[1985] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*2*I_NAI_F2xz_G2xyz_b;
  abcd[1986] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab;
  abcd[1987] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Fx2y_G2xyz_b;
  abcd[1988] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[1989] = 4.0E0*I_NAI_Hx4z_G2xyz_ab-2.0E0*3*I_NAI_Fx2z_G2xyz_b;
  abcd[1990] = 4.0E0*I_NAI_H4yz_G2xyz_ab;
  abcd[1991] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*1*I_NAI_F3y_G2xyz_b;
  abcd[1992] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*2*I_NAI_F2yz_G2xyz_b;
  abcd[1993] = 4.0E0*I_NAI_Hy4z_G2xyz_ab-2.0E0*3*I_NAI_Fy2z_G2xyz_b;
  abcd[1994] = 4.0E0*I_NAI_H5z_G2xyz_ab-2.0E0*4*I_NAI_F3z_G2xyz_b;
  abcd[1995] = 4.0E0*I_NAI_H4xz_Gx3y_ab-2.0E0*2*I_NAI_H4xz_Dxy_a;
  abcd[1996] = 4.0E0*I_NAI_H3xyz_Gx3y_ab-2.0E0*2*I_NAI_H3xyz_Dxy_a;
  abcd[1997] = 4.0E0*I_NAI_H3x2z_Gx3y_ab-2.0E0*2*I_NAI_H3x2z_Dxy_a-2.0E0*1*I_NAI_F3x_Gx3y_b+2*I_NAI_F3x_Dxy;
  abcd[1998] = 4.0E0*I_NAI_H2x2yz_Gx3y_ab-2.0E0*2*I_NAI_H2x2yz_Dxy_a;
  abcd[1999] = 4.0E0*I_NAI_H2xy2z_Gx3y_ab-2.0E0*2*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_F2xy_Gx3y_b+2*I_NAI_F2xy_Dxy;
  abcd[2000] = 4.0E0*I_NAI_H2x3z_Gx3y_ab-2.0E0*2*I_NAI_H2x3z_Dxy_a-2.0E0*2*I_NAI_F2xz_Gx3y_b+2*2*I_NAI_F2xz_Dxy;
  abcd[2001] = 4.0E0*I_NAI_Hx3yz_Gx3y_ab-2.0E0*2*I_NAI_Hx3yz_Dxy_a;
  abcd[2002] = 4.0E0*I_NAI_Hx2y2z_Gx3y_ab-2.0E0*2*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2y_Gx3y_b+2*I_NAI_Fx2y_Dxy;
  abcd[2003] = 4.0E0*I_NAI_Hxy3z_Gx3y_ab-2.0E0*2*I_NAI_Hxy3z_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gx3y_b+2*2*I_NAI_Fxyz_Dxy;
  abcd[2004] = 4.0E0*I_NAI_Hx4z_Gx3y_ab-2.0E0*2*I_NAI_Hx4z_Dxy_a-2.0E0*3*I_NAI_Fx2z_Gx3y_b+3*2*I_NAI_Fx2z_Dxy;
  abcd[2005] = 4.0E0*I_NAI_H4yz_Gx3y_ab-2.0E0*2*I_NAI_H4yz_Dxy_a;
  abcd[2006] = 4.0E0*I_NAI_H3y2z_Gx3y_ab-2.0E0*2*I_NAI_H3y2z_Dxy_a-2.0E0*1*I_NAI_F3y_Gx3y_b+2*I_NAI_F3y_Dxy;
  abcd[2007] = 4.0E0*I_NAI_H2y3z_Gx3y_ab-2.0E0*2*I_NAI_H2y3z_Dxy_a-2.0E0*2*I_NAI_F2yz_Gx3y_b+2*2*I_NAI_F2yz_Dxy;
  abcd[2008] = 4.0E0*I_NAI_Hy4z_Gx3y_ab-2.0E0*2*I_NAI_Hy4z_Dxy_a-2.0E0*3*I_NAI_Fy2z_Gx3y_b+3*2*I_NAI_Fy2z_Dxy;
  abcd[2009] = 4.0E0*I_NAI_H5z_Gx3y_ab-2.0E0*2*I_NAI_H5z_Dxy_a-2.0E0*4*I_NAI_F3z_Gx3y_b+4*2*I_NAI_F3z_Dxy;
  abcd[2010] = 4.0E0*I_NAI_H4xz_Gx2yz_ab-2.0E0*1*I_NAI_H4xz_Dxz_a;
  abcd[2011] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab-2.0E0*1*I_NAI_H3xyz_Dxz_a;
  abcd[2012] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*1*I_NAI_H3x2z_Dxz_a-2.0E0*1*I_NAI_F3x_Gx2yz_b+1*I_NAI_F3x_Dxz;
  abcd[2013] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab-2.0E0*1*I_NAI_H2x2yz_Dxz_a;
  abcd[2014] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_F2xy_Gx2yz_b+1*I_NAI_F2xy_Dxz;
  abcd[2015] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*1*I_NAI_H2x3z_Dxz_a-2.0E0*2*I_NAI_F2xz_Gx2yz_b+2*1*I_NAI_F2xz_Dxz;
  abcd[2016] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab-2.0E0*1*I_NAI_Hx3yz_Dxz_a;
  abcd[2017] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2y_Gx2yz_b+1*I_NAI_Fx2y_Dxz;
  abcd[2018] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab-2.0E0*1*I_NAI_Hxy3z_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx2yz_b+2*1*I_NAI_Fxyz_Dxz;
  abcd[2019] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab-2.0E0*1*I_NAI_Hx4z_Dxz_a-2.0E0*3*I_NAI_Fx2z_Gx2yz_b+3*1*I_NAI_Fx2z_Dxz;
  abcd[2020] = 4.0E0*I_NAI_H4yz_Gx2yz_ab-2.0E0*1*I_NAI_H4yz_Dxz_a;
  abcd[2021] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*1*I_NAI_H3y2z_Dxz_a-2.0E0*1*I_NAI_F3y_Gx2yz_b+1*I_NAI_F3y_Dxz;
  abcd[2022] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*1*I_NAI_H2y3z_Dxz_a-2.0E0*2*I_NAI_F2yz_Gx2yz_b+2*1*I_NAI_F2yz_Dxz;
  abcd[2023] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab-2.0E0*1*I_NAI_Hy4z_Dxz_a-2.0E0*3*I_NAI_Fy2z_Gx2yz_b+3*1*I_NAI_Fy2z_Dxz;
  abcd[2024] = 4.0E0*I_NAI_H5z_Gx2yz_ab-2.0E0*1*I_NAI_H5z_Dxz_a-2.0E0*4*I_NAI_F3z_Gx2yz_b+4*1*I_NAI_F3z_Dxz;
  abcd[2025] = 4.0E0*I_NAI_H4xz_Gxy2z_ab;
  abcd[2026] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab;
  abcd[2027] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*1*I_NAI_F3x_Gxy2z_b;
  abcd[2028] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab;
  abcd[2029] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_F2xy_Gxy2z_b;
  abcd[2030] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*2*I_NAI_F2xz_Gxy2z_b;
  abcd[2031] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab;
  abcd[2032] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Fx2y_Gxy2z_b;
  abcd[2033] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab-2.0E0*2*I_NAI_Fxyz_Gxy2z_b;
  abcd[2034] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab-2.0E0*3*I_NAI_Fx2z_Gxy2z_b;
  abcd[2035] = 4.0E0*I_NAI_H4yz_Gxy2z_ab;
  abcd[2036] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*1*I_NAI_F3y_Gxy2z_b;
  abcd[2037] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*2*I_NAI_F2yz_Gxy2z_b;
  abcd[2038] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab-2.0E0*3*I_NAI_Fy2z_Gxy2z_b;
  abcd[2039] = 4.0E0*I_NAI_H5z_Gxy2z_ab-2.0E0*4*I_NAI_F3z_Gxy2z_b;
  abcd[2040] = 4.0E0*I_NAI_H4xz_G4y_ab-2.0E0*3*I_NAI_H4xz_D2y_a;
  abcd[2041] = 4.0E0*I_NAI_H3xyz_G4y_ab-2.0E0*3*I_NAI_H3xyz_D2y_a;
  abcd[2042] = 4.0E0*I_NAI_H3x2z_G4y_ab-2.0E0*3*I_NAI_H3x2z_D2y_a-2.0E0*1*I_NAI_F3x_G4y_b+3*I_NAI_F3x_D2y;
  abcd[2043] = 4.0E0*I_NAI_H2x2yz_G4y_ab-2.0E0*3*I_NAI_H2x2yz_D2y_a;
  abcd[2044] = 4.0E0*I_NAI_H2xy2z_G4y_ab-2.0E0*3*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_F2xy_G4y_b+3*I_NAI_F2xy_D2y;
  abcd[2045] = 4.0E0*I_NAI_H2x3z_G4y_ab-2.0E0*3*I_NAI_H2x3z_D2y_a-2.0E0*2*I_NAI_F2xz_G4y_b+2*3*I_NAI_F2xz_D2y;
  abcd[2046] = 4.0E0*I_NAI_Hx3yz_G4y_ab-2.0E0*3*I_NAI_Hx3yz_D2y_a;
  abcd[2047] = 4.0E0*I_NAI_Hx2y2z_G4y_ab-2.0E0*3*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2y_G4y_b+3*I_NAI_Fx2y_D2y;
  abcd[2048] = 4.0E0*I_NAI_Hxy3z_G4y_ab-2.0E0*3*I_NAI_Hxy3z_D2y_a-2.0E0*2*I_NAI_Fxyz_G4y_b+2*3*I_NAI_Fxyz_D2y;
  abcd[2049] = 4.0E0*I_NAI_Hx4z_G4y_ab-2.0E0*3*I_NAI_Hx4z_D2y_a-2.0E0*3*I_NAI_Fx2z_G4y_b+3*3*I_NAI_Fx2z_D2y;
  abcd[2050] = 4.0E0*I_NAI_H4yz_G4y_ab-2.0E0*3*I_NAI_H4yz_D2y_a;
  abcd[2051] = 4.0E0*I_NAI_H3y2z_G4y_ab-2.0E0*3*I_NAI_H3y2z_D2y_a-2.0E0*1*I_NAI_F3y_G4y_b+3*I_NAI_F3y_D2y;
  abcd[2052] = 4.0E0*I_NAI_H2y3z_G4y_ab-2.0E0*3*I_NAI_H2y3z_D2y_a-2.0E0*2*I_NAI_F2yz_G4y_b+2*3*I_NAI_F2yz_D2y;
  abcd[2053] = 4.0E0*I_NAI_Hy4z_G4y_ab-2.0E0*3*I_NAI_Hy4z_D2y_a-2.0E0*3*I_NAI_Fy2z_G4y_b+3*3*I_NAI_Fy2z_D2y;
  abcd[2054] = 4.0E0*I_NAI_H5z_G4y_ab-2.0E0*3*I_NAI_H5z_D2y_a-2.0E0*4*I_NAI_F3z_G4y_b+4*3*I_NAI_F3z_D2y;
  abcd[2055] = 4.0E0*I_NAI_H4xz_G3yz_ab-2.0E0*2*I_NAI_H4xz_Dyz_a;
  abcd[2056] = 4.0E0*I_NAI_H3xyz_G3yz_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a;
  abcd[2057] = 4.0E0*I_NAI_H3x2z_G3yz_ab-2.0E0*2*I_NAI_H3x2z_Dyz_a-2.0E0*1*I_NAI_F3x_G3yz_b+2*I_NAI_F3x_Dyz;
  abcd[2058] = 4.0E0*I_NAI_H2x2yz_G3yz_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a;
  abcd[2059] = 4.0E0*I_NAI_H2xy2z_G3yz_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_F2xy_G3yz_b+2*I_NAI_F2xy_Dyz;
  abcd[2060] = 4.0E0*I_NAI_H2x3z_G3yz_ab-2.0E0*2*I_NAI_H2x3z_Dyz_a-2.0E0*2*I_NAI_F2xz_G3yz_b+2*2*I_NAI_F2xz_Dyz;
  abcd[2061] = 4.0E0*I_NAI_Hx3yz_G3yz_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a;
  abcd[2062] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2y_G3yz_b+2*I_NAI_Fx2y_Dyz;
  abcd[2063] = 4.0E0*I_NAI_Hxy3z_G3yz_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a-2.0E0*2*I_NAI_Fxyz_G3yz_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[2064] = 4.0E0*I_NAI_Hx4z_G3yz_ab-2.0E0*2*I_NAI_Hx4z_Dyz_a-2.0E0*3*I_NAI_Fx2z_G3yz_b+3*2*I_NAI_Fx2z_Dyz;
  abcd[2065] = 4.0E0*I_NAI_H4yz_G3yz_ab-2.0E0*2*I_NAI_H4yz_Dyz_a;
  abcd[2066] = 4.0E0*I_NAI_H3y2z_G3yz_ab-2.0E0*2*I_NAI_H3y2z_Dyz_a-2.0E0*1*I_NAI_F3y_G3yz_b+2*I_NAI_F3y_Dyz;
  abcd[2067] = 4.0E0*I_NAI_H2y3z_G3yz_ab-2.0E0*2*I_NAI_H2y3z_Dyz_a-2.0E0*2*I_NAI_F2yz_G3yz_b+2*2*I_NAI_F2yz_Dyz;
  abcd[2068] = 4.0E0*I_NAI_Hy4z_G3yz_ab-2.0E0*2*I_NAI_Hy4z_Dyz_a-2.0E0*3*I_NAI_Fy2z_G3yz_b+3*2*I_NAI_Fy2z_Dyz;
  abcd[2069] = 4.0E0*I_NAI_H5z_G3yz_ab-2.0E0*2*I_NAI_H5z_Dyz_a-2.0E0*4*I_NAI_F3z_G3yz_b+4*2*I_NAI_F3z_Dyz;
  abcd[2070] = 4.0E0*I_NAI_H4xz_G2y2z_ab-2.0E0*1*I_NAI_H4xz_D2z_a;
  abcd[2071] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2z_a;
  abcd[2072] = 4.0E0*I_NAI_H3x2z_G2y2z_ab-2.0E0*1*I_NAI_H3x2z_D2z_a-2.0E0*1*I_NAI_F3x_G2y2z_b+1*I_NAI_F3x_D2z;
  abcd[2073] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2z_a;
  abcd[2074] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_F2xy_G2y2z_b+1*I_NAI_F2xy_D2z;
  abcd[2075] = 4.0E0*I_NAI_H2x3z_G2y2z_ab-2.0E0*1*I_NAI_H2x3z_D2z_a-2.0E0*2*I_NAI_F2xz_G2y2z_b+2*1*I_NAI_F2xz_D2z;
  abcd[2076] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2z_a;
  abcd[2077] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2y_G2y2z_b+1*I_NAI_Fx2y_D2z;
  abcd[2078] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2z_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2z;
  abcd[2079] = 4.0E0*I_NAI_Hx4z_G2y2z_ab-2.0E0*1*I_NAI_Hx4z_D2z_a-2.0E0*3*I_NAI_Fx2z_G2y2z_b+3*1*I_NAI_Fx2z_D2z;
  abcd[2080] = 4.0E0*I_NAI_H4yz_G2y2z_ab-2.0E0*1*I_NAI_H4yz_D2z_a;
  abcd[2081] = 4.0E0*I_NAI_H3y2z_G2y2z_ab-2.0E0*1*I_NAI_H3y2z_D2z_a-2.0E0*1*I_NAI_F3y_G2y2z_b+1*I_NAI_F3y_D2z;
  abcd[2082] = 4.0E0*I_NAI_H2y3z_G2y2z_ab-2.0E0*1*I_NAI_H2y3z_D2z_a-2.0E0*2*I_NAI_F2yz_G2y2z_b+2*1*I_NAI_F2yz_D2z;
  abcd[2083] = 4.0E0*I_NAI_Hy4z_G2y2z_ab-2.0E0*1*I_NAI_Hy4z_D2z_a-2.0E0*3*I_NAI_Fy2z_G2y2z_b+3*1*I_NAI_Fy2z_D2z;
  abcd[2084] = 4.0E0*I_NAI_H5z_G2y2z_ab-2.0E0*1*I_NAI_H5z_D2z_a-2.0E0*4*I_NAI_F3z_G2y2z_b+4*1*I_NAI_F3z_D2z;
  abcd[2085] = 4.0E0*I_NAI_H4xz_Gy3z_ab;
  abcd[2086] = 4.0E0*I_NAI_H3xyz_Gy3z_ab;
  abcd[2087] = 4.0E0*I_NAI_H3x2z_Gy3z_ab-2.0E0*1*I_NAI_F3x_Gy3z_b;
  abcd[2088] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab;
  abcd[2089] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab-2.0E0*1*I_NAI_F2xy_Gy3z_b;
  abcd[2090] = 4.0E0*I_NAI_H2x3z_Gy3z_ab-2.0E0*2*I_NAI_F2xz_Gy3z_b;
  abcd[2091] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab;
  abcd[2092] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab-2.0E0*1*I_NAI_Fx2y_Gy3z_b;
  abcd[2093] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab-2.0E0*2*I_NAI_Fxyz_Gy3z_b;
  abcd[2094] = 4.0E0*I_NAI_Hx4z_Gy3z_ab-2.0E0*3*I_NAI_Fx2z_Gy3z_b;
  abcd[2095] = 4.0E0*I_NAI_H4yz_Gy3z_ab;
  abcd[2096] = 4.0E0*I_NAI_H3y2z_Gy3z_ab-2.0E0*1*I_NAI_F3y_Gy3z_b;
  abcd[2097] = 4.0E0*I_NAI_H2y3z_Gy3z_ab-2.0E0*2*I_NAI_F2yz_Gy3z_b;
  abcd[2098] = 4.0E0*I_NAI_Hy4z_Gy3z_ab-2.0E0*3*I_NAI_Fy2z_Gy3z_b;
  abcd[2099] = 4.0E0*I_NAI_H5z_Gy3z_ab-2.0E0*4*I_NAI_F3z_Gy3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_ab
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[2100] = 4.0E0*I_NAI_H4xz_G3xz_ab;
  abcd[2101] = 4.0E0*I_NAI_H3xyz_G3xz_ab;
  abcd[2102] = 4.0E0*I_NAI_H3x2z_G3xz_ab-2.0E0*1*I_NAI_F3x_G3xz_b;
  abcd[2103] = 4.0E0*I_NAI_H2x2yz_G3xz_ab;
  abcd[2104] = 4.0E0*I_NAI_H2xy2z_G3xz_ab-2.0E0*1*I_NAI_F2xy_G3xz_b;
  abcd[2105] = 4.0E0*I_NAI_H2x3z_G3xz_ab-2.0E0*2*I_NAI_F2xz_G3xz_b;
  abcd[2106] = 4.0E0*I_NAI_Hx3yz_G3xz_ab;
  abcd[2107] = 4.0E0*I_NAI_Hx2y2z_G3xz_ab-2.0E0*1*I_NAI_Fx2y_G3xz_b;
  abcd[2108] = 4.0E0*I_NAI_Hxy3z_G3xz_ab-2.0E0*2*I_NAI_Fxyz_G3xz_b;
  abcd[2109] = 4.0E0*I_NAI_Hx4z_G3xz_ab-2.0E0*3*I_NAI_Fx2z_G3xz_b;
  abcd[2110] = 4.0E0*I_NAI_H4yz_G3xz_ab;
  abcd[2111] = 4.0E0*I_NAI_H3y2z_G3xz_ab-2.0E0*1*I_NAI_F3y_G3xz_b;
  abcd[2112] = 4.0E0*I_NAI_H2y3z_G3xz_ab-2.0E0*2*I_NAI_F2yz_G3xz_b;
  abcd[2113] = 4.0E0*I_NAI_Hy4z_G3xz_ab-2.0E0*3*I_NAI_Fy2z_G3xz_b;
  abcd[2114] = 4.0E0*I_NAI_H5z_G3xz_ab-2.0E0*4*I_NAI_F3z_G3xz_b;
  abcd[2115] = 4.0E0*I_NAI_H4xz_G2xyz_ab;
  abcd[2116] = 4.0E0*I_NAI_H3xyz_G2xyz_ab;
  abcd[2117] = 4.0E0*I_NAI_H3x2z_G2xyz_ab-2.0E0*1*I_NAI_F3x_G2xyz_b;
  abcd[2118] = 4.0E0*I_NAI_H2x2yz_G2xyz_ab;
  abcd[2119] = 4.0E0*I_NAI_H2xy2z_G2xyz_ab-2.0E0*1*I_NAI_F2xy_G2xyz_b;
  abcd[2120] = 4.0E0*I_NAI_H2x3z_G2xyz_ab-2.0E0*2*I_NAI_F2xz_G2xyz_b;
  abcd[2121] = 4.0E0*I_NAI_Hx3yz_G2xyz_ab;
  abcd[2122] = 4.0E0*I_NAI_Hx2y2z_G2xyz_ab-2.0E0*1*I_NAI_Fx2y_G2xyz_b;
  abcd[2123] = 4.0E0*I_NAI_Hxy3z_G2xyz_ab-2.0E0*2*I_NAI_Fxyz_G2xyz_b;
  abcd[2124] = 4.0E0*I_NAI_Hx4z_G2xyz_ab-2.0E0*3*I_NAI_Fx2z_G2xyz_b;
  abcd[2125] = 4.0E0*I_NAI_H4yz_G2xyz_ab;
  abcd[2126] = 4.0E0*I_NAI_H3y2z_G2xyz_ab-2.0E0*1*I_NAI_F3y_G2xyz_b;
  abcd[2127] = 4.0E0*I_NAI_H2y3z_G2xyz_ab-2.0E0*2*I_NAI_F2yz_G2xyz_b;
  abcd[2128] = 4.0E0*I_NAI_Hy4z_G2xyz_ab-2.0E0*3*I_NAI_Fy2z_G2xyz_b;
  abcd[2129] = 4.0E0*I_NAI_H5z_G2xyz_ab-2.0E0*4*I_NAI_F3z_G2xyz_b;
  abcd[2130] = 4.0E0*I_NAI_H4xz_G2x2z_ab-2.0E0*1*I_NAI_H4xz_D2x_a;
  abcd[2131] = 4.0E0*I_NAI_H3xyz_G2x2z_ab-2.0E0*1*I_NAI_H3xyz_D2x_a;
  abcd[2132] = 4.0E0*I_NAI_H3x2z_G2x2z_ab-2.0E0*1*I_NAI_H3x2z_D2x_a-2.0E0*1*I_NAI_F3x_G2x2z_b+1*I_NAI_F3x_D2x;
  abcd[2133] = 4.0E0*I_NAI_H2x2yz_G2x2z_ab-2.0E0*1*I_NAI_H2x2yz_D2x_a;
  abcd[2134] = 4.0E0*I_NAI_H2xy2z_G2x2z_ab-2.0E0*1*I_NAI_H2xy2z_D2x_a-2.0E0*1*I_NAI_F2xy_G2x2z_b+1*I_NAI_F2xy_D2x;
  abcd[2135] = 4.0E0*I_NAI_H2x3z_G2x2z_ab-2.0E0*1*I_NAI_H2x3z_D2x_a-2.0E0*2*I_NAI_F2xz_G2x2z_b+2*1*I_NAI_F2xz_D2x;
  abcd[2136] = 4.0E0*I_NAI_Hx3yz_G2x2z_ab-2.0E0*1*I_NAI_Hx3yz_D2x_a;
  abcd[2137] = 4.0E0*I_NAI_Hx2y2z_G2x2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2x_a-2.0E0*1*I_NAI_Fx2y_G2x2z_b+1*I_NAI_Fx2y_D2x;
  abcd[2138] = 4.0E0*I_NAI_Hxy3z_G2x2z_ab-2.0E0*1*I_NAI_Hxy3z_D2x_a-2.0E0*2*I_NAI_Fxyz_G2x2z_b+2*1*I_NAI_Fxyz_D2x;
  abcd[2139] = 4.0E0*I_NAI_Hx4z_G2x2z_ab-2.0E0*1*I_NAI_Hx4z_D2x_a-2.0E0*3*I_NAI_Fx2z_G2x2z_b+3*1*I_NAI_Fx2z_D2x;
  abcd[2140] = 4.0E0*I_NAI_H4yz_G2x2z_ab-2.0E0*1*I_NAI_H4yz_D2x_a;
  abcd[2141] = 4.0E0*I_NAI_H3y2z_G2x2z_ab-2.0E0*1*I_NAI_H3y2z_D2x_a-2.0E0*1*I_NAI_F3y_G2x2z_b+1*I_NAI_F3y_D2x;
  abcd[2142] = 4.0E0*I_NAI_H2y3z_G2x2z_ab-2.0E0*1*I_NAI_H2y3z_D2x_a-2.0E0*2*I_NAI_F2yz_G2x2z_b+2*1*I_NAI_F2yz_D2x;
  abcd[2143] = 4.0E0*I_NAI_Hy4z_G2x2z_ab-2.0E0*1*I_NAI_Hy4z_D2x_a-2.0E0*3*I_NAI_Fy2z_G2x2z_b+3*1*I_NAI_Fy2z_D2x;
  abcd[2144] = 4.0E0*I_NAI_H5z_G2x2z_ab-2.0E0*1*I_NAI_H5z_D2x_a-2.0E0*4*I_NAI_F3z_G2x2z_b+4*1*I_NAI_F3z_D2x;
  abcd[2145] = 4.0E0*I_NAI_H4xz_Gx2yz_ab;
  abcd[2146] = 4.0E0*I_NAI_H3xyz_Gx2yz_ab;
  abcd[2147] = 4.0E0*I_NAI_H3x2z_Gx2yz_ab-2.0E0*1*I_NAI_F3x_Gx2yz_b;
  abcd[2148] = 4.0E0*I_NAI_H2x2yz_Gx2yz_ab;
  abcd[2149] = 4.0E0*I_NAI_H2xy2z_Gx2yz_ab-2.0E0*1*I_NAI_F2xy_Gx2yz_b;
  abcd[2150] = 4.0E0*I_NAI_H2x3z_Gx2yz_ab-2.0E0*2*I_NAI_F2xz_Gx2yz_b;
  abcd[2151] = 4.0E0*I_NAI_Hx3yz_Gx2yz_ab;
  abcd[2152] = 4.0E0*I_NAI_Hx2y2z_Gx2yz_ab-2.0E0*1*I_NAI_Fx2y_Gx2yz_b;
  abcd[2153] = 4.0E0*I_NAI_Hxy3z_Gx2yz_ab-2.0E0*2*I_NAI_Fxyz_Gx2yz_b;
  abcd[2154] = 4.0E0*I_NAI_Hx4z_Gx2yz_ab-2.0E0*3*I_NAI_Fx2z_Gx2yz_b;
  abcd[2155] = 4.0E0*I_NAI_H4yz_Gx2yz_ab;
  abcd[2156] = 4.0E0*I_NAI_H3y2z_Gx2yz_ab-2.0E0*1*I_NAI_F3y_Gx2yz_b;
  abcd[2157] = 4.0E0*I_NAI_H2y3z_Gx2yz_ab-2.0E0*2*I_NAI_F2yz_Gx2yz_b;
  abcd[2158] = 4.0E0*I_NAI_Hy4z_Gx2yz_ab-2.0E0*3*I_NAI_Fy2z_Gx2yz_b;
  abcd[2159] = 4.0E0*I_NAI_H5z_Gx2yz_ab-2.0E0*4*I_NAI_F3z_Gx2yz_b;
  abcd[2160] = 4.0E0*I_NAI_H4xz_Gxy2z_ab-2.0E0*1*I_NAI_H4xz_Dxy_a;
  abcd[2161] = 4.0E0*I_NAI_H3xyz_Gxy2z_ab-2.0E0*1*I_NAI_H3xyz_Dxy_a;
  abcd[2162] = 4.0E0*I_NAI_H3x2z_Gxy2z_ab-2.0E0*1*I_NAI_H3x2z_Dxy_a-2.0E0*1*I_NAI_F3x_Gxy2z_b+1*I_NAI_F3x_Dxy;
  abcd[2163] = 4.0E0*I_NAI_H2x2yz_Gxy2z_ab-2.0E0*1*I_NAI_H2x2yz_Dxy_a;
  abcd[2164] = 4.0E0*I_NAI_H2xy2z_Gxy2z_ab-2.0E0*1*I_NAI_H2xy2z_Dxy_a-2.0E0*1*I_NAI_F2xy_Gxy2z_b+1*I_NAI_F2xy_Dxy;
  abcd[2165] = 4.0E0*I_NAI_H2x3z_Gxy2z_ab-2.0E0*1*I_NAI_H2x3z_Dxy_a-2.0E0*2*I_NAI_F2xz_Gxy2z_b+2*1*I_NAI_F2xz_Dxy;
  abcd[2166] = 4.0E0*I_NAI_Hx3yz_Gxy2z_ab-2.0E0*1*I_NAI_Hx3yz_Dxy_a;
  abcd[2167] = 4.0E0*I_NAI_Hx2y2z_Gxy2z_ab-2.0E0*1*I_NAI_Hx2y2z_Dxy_a-2.0E0*1*I_NAI_Fx2y_Gxy2z_b+1*I_NAI_Fx2y_Dxy;
  abcd[2168] = 4.0E0*I_NAI_Hxy3z_Gxy2z_ab-2.0E0*1*I_NAI_Hxy3z_Dxy_a-2.0E0*2*I_NAI_Fxyz_Gxy2z_b+2*1*I_NAI_Fxyz_Dxy;
  abcd[2169] = 4.0E0*I_NAI_Hx4z_Gxy2z_ab-2.0E0*1*I_NAI_Hx4z_Dxy_a-2.0E0*3*I_NAI_Fx2z_Gxy2z_b+3*1*I_NAI_Fx2z_Dxy;
  abcd[2170] = 4.0E0*I_NAI_H4yz_Gxy2z_ab-2.0E0*1*I_NAI_H4yz_Dxy_a;
  abcd[2171] = 4.0E0*I_NAI_H3y2z_Gxy2z_ab-2.0E0*1*I_NAI_H3y2z_Dxy_a-2.0E0*1*I_NAI_F3y_Gxy2z_b+1*I_NAI_F3y_Dxy;
  abcd[2172] = 4.0E0*I_NAI_H2y3z_Gxy2z_ab-2.0E0*1*I_NAI_H2y3z_Dxy_a-2.0E0*2*I_NAI_F2yz_Gxy2z_b+2*1*I_NAI_F2yz_Dxy;
  abcd[2173] = 4.0E0*I_NAI_Hy4z_Gxy2z_ab-2.0E0*1*I_NAI_Hy4z_Dxy_a-2.0E0*3*I_NAI_Fy2z_Gxy2z_b+3*1*I_NAI_Fy2z_Dxy;
  abcd[2174] = 4.0E0*I_NAI_H5z_Gxy2z_ab-2.0E0*1*I_NAI_H5z_Dxy_a-2.0E0*4*I_NAI_F3z_Gxy2z_b+4*1*I_NAI_F3z_Dxy;
  abcd[2175] = 4.0E0*I_NAI_H4xz_Gx3z_ab-2.0E0*2*I_NAI_H4xz_Dxz_a;
  abcd[2176] = 4.0E0*I_NAI_H3xyz_Gx3z_ab-2.0E0*2*I_NAI_H3xyz_Dxz_a;
  abcd[2177] = 4.0E0*I_NAI_H3x2z_Gx3z_ab-2.0E0*2*I_NAI_H3x2z_Dxz_a-2.0E0*1*I_NAI_F3x_Gx3z_b+2*I_NAI_F3x_Dxz;
  abcd[2178] = 4.0E0*I_NAI_H2x2yz_Gx3z_ab-2.0E0*2*I_NAI_H2x2yz_Dxz_a;
  abcd[2179] = 4.0E0*I_NAI_H2xy2z_Gx3z_ab-2.0E0*2*I_NAI_H2xy2z_Dxz_a-2.0E0*1*I_NAI_F2xy_Gx3z_b+2*I_NAI_F2xy_Dxz;
  abcd[2180] = 4.0E0*I_NAI_H2x3z_Gx3z_ab-2.0E0*2*I_NAI_H2x3z_Dxz_a-2.0E0*2*I_NAI_F2xz_Gx3z_b+2*2*I_NAI_F2xz_Dxz;
  abcd[2181] = 4.0E0*I_NAI_Hx3yz_Gx3z_ab-2.0E0*2*I_NAI_Hx3yz_Dxz_a;
  abcd[2182] = 4.0E0*I_NAI_Hx2y2z_Gx3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dxz_a-2.0E0*1*I_NAI_Fx2y_Gx3z_b+2*I_NAI_Fx2y_Dxz;
  abcd[2183] = 4.0E0*I_NAI_Hxy3z_Gx3z_ab-2.0E0*2*I_NAI_Hxy3z_Dxz_a-2.0E0*2*I_NAI_Fxyz_Gx3z_b+2*2*I_NAI_Fxyz_Dxz;
  abcd[2184] = 4.0E0*I_NAI_Hx4z_Gx3z_ab-2.0E0*2*I_NAI_Hx4z_Dxz_a-2.0E0*3*I_NAI_Fx2z_Gx3z_b+3*2*I_NAI_Fx2z_Dxz;
  abcd[2185] = 4.0E0*I_NAI_H4yz_Gx3z_ab-2.0E0*2*I_NAI_H4yz_Dxz_a;
  abcd[2186] = 4.0E0*I_NAI_H3y2z_Gx3z_ab-2.0E0*2*I_NAI_H3y2z_Dxz_a-2.0E0*1*I_NAI_F3y_Gx3z_b+2*I_NAI_F3y_Dxz;
  abcd[2187] = 4.0E0*I_NAI_H2y3z_Gx3z_ab-2.0E0*2*I_NAI_H2y3z_Dxz_a-2.0E0*2*I_NAI_F2yz_Gx3z_b+2*2*I_NAI_F2yz_Dxz;
  abcd[2188] = 4.0E0*I_NAI_Hy4z_Gx3z_ab-2.0E0*2*I_NAI_Hy4z_Dxz_a-2.0E0*3*I_NAI_Fy2z_Gx3z_b+3*2*I_NAI_Fy2z_Dxz;
  abcd[2189] = 4.0E0*I_NAI_H5z_Gx3z_ab-2.0E0*2*I_NAI_H5z_Dxz_a-2.0E0*4*I_NAI_F3z_Gx3z_b+4*2*I_NAI_F3z_Dxz;
  abcd[2190] = 4.0E0*I_NAI_H4xz_G3yz_ab;
  abcd[2191] = 4.0E0*I_NAI_H3xyz_G3yz_ab;
  abcd[2192] = 4.0E0*I_NAI_H3x2z_G3yz_ab-2.0E0*1*I_NAI_F3x_G3yz_b;
  abcd[2193] = 4.0E0*I_NAI_H2x2yz_G3yz_ab;
  abcd[2194] = 4.0E0*I_NAI_H2xy2z_G3yz_ab-2.0E0*1*I_NAI_F2xy_G3yz_b;
  abcd[2195] = 4.0E0*I_NAI_H2x3z_G3yz_ab-2.0E0*2*I_NAI_F2xz_G3yz_b;
  abcd[2196] = 4.0E0*I_NAI_Hx3yz_G3yz_ab;
  abcd[2197] = 4.0E0*I_NAI_Hx2y2z_G3yz_ab-2.0E0*1*I_NAI_Fx2y_G3yz_b;
  abcd[2198] = 4.0E0*I_NAI_Hxy3z_G3yz_ab-2.0E0*2*I_NAI_Fxyz_G3yz_b;
  abcd[2199] = 4.0E0*I_NAI_Hx4z_G3yz_ab-2.0E0*3*I_NAI_Fx2z_G3yz_b;
  abcd[2200] = 4.0E0*I_NAI_H4yz_G3yz_ab;
  abcd[2201] = 4.0E0*I_NAI_H3y2z_G3yz_ab-2.0E0*1*I_NAI_F3y_G3yz_b;
  abcd[2202] = 4.0E0*I_NAI_H2y3z_G3yz_ab-2.0E0*2*I_NAI_F2yz_G3yz_b;
  abcd[2203] = 4.0E0*I_NAI_Hy4z_G3yz_ab-2.0E0*3*I_NAI_Fy2z_G3yz_b;
  abcd[2204] = 4.0E0*I_NAI_H5z_G3yz_ab-2.0E0*4*I_NAI_F3z_G3yz_b;
  abcd[2205] = 4.0E0*I_NAI_H4xz_G2y2z_ab-2.0E0*1*I_NAI_H4xz_D2y_a;
  abcd[2206] = 4.0E0*I_NAI_H3xyz_G2y2z_ab-2.0E0*1*I_NAI_H3xyz_D2y_a;
  abcd[2207] = 4.0E0*I_NAI_H3x2z_G2y2z_ab-2.0E0*1*I_NAI_H3x2z_D2y_a-2.0E0*1*I_NAI_F3x_G2y2z_b+1*I_NAI_F3x_D2y;
  abcd[2208] = 4.0E0*I_NAI_H2x2yz_G2y2z_ab-2.0E0*1*I_NAI_H2x2yz_D2y_a;
  abcd[2209] = 4.0E0*I_NAI_H2xy2z_G2y2z_ab-2.0E0*1*I_NAI_H2xy2z_D2y_a-2.0E0*1*I_NAI_F2xy_G2y2z_b+1*I_NAI_F2xy_D2y;
  abcd[2210] = 4.0E0*I_NAI_H2x3z_G2y2z_ab-2.0E0*1*I_NAI_H2x3z_D2y_a-2.0E0*2*I_NAI_F2xz_G2y2z_b+2*1*I_NAI_F2xz_D2y;
  abcd[2211] = 4.0E0*I_NAI_Hx3yz_G2y2z_ab-2.0E0*1*I_NAI_Hx3yz_D2y_a;
  abcd[2212] = 4.0E0*I_NAI_Hx2y2z_G2y2z_ab-2.0E0*1*I_NAI_Hx2y2z_D2y_a-2.0E0*1*I_NAI_Fx2y_G2y2z_b+1*I_NAI_Fx2y_D2y;
  abcd[2213] = 4.0E0*I_NAI_Hxy3z_G2y2z_ab-2.0E0*1*I_NAI_Hxy3z_D2y_a-2.0E0*2*I_NAI_Fxyz_G2y2z_b+2*1*I_NAI_Fxyz_D2y;
  abcd[2214] = 4.0E0*I_NAI_Hx4z_G2y2z_ab-2.0E0*1*I_NAI_Hx4z_D2y_a-2.0E0*3*I_NAI_Fx2z_G2y2z_b+3*1*I_NAI_Fx2z_D2y;
  abcd[2215] = 4.0E0*I_NAI_H4yz_G2y2z_ab-2.0E0*1*I_NAI_H4yz_D2y_a;
  abcd[2216] = 4.0E0*I_NAI_H3y2z_G2y2z_ab-2.0E0*1*I_NAI_H3y2z_D2y_a-2.0E0*1*I_NAI_F3y_G2y2z_b+1*I_NAI_F3y_D2y;
  abcd[2217] = 4.0E0*I_NAI_H2y3z_G2y2z_ab-2.0E0*1*I_NAI_H2y3z_D2y_a-2.0E0*2*I_NAI_F2yz_G2y2z_b+2*1*I_NAI_F2yz_D2y;
  abcd[2218] = 4.0E0*I_NAI_Hy4z_G2y2z_ab-2.0E0*1*I_NAI_Hy4z_D2y_a-2.0E0*3*I_NAI_Fy2z_G2y2z_b+3*1*I_NAI_Fy2z_D2y;
  abcd[2219] = 4.0E0*I_NAI_H5z_G2y2z_ab-2.0E0*1*I_NAI_H5z_D2y_a-2.0E0*4*I_NAI_F3z_G2y2z_b+4*1*I_NAI_F3z_D2y;
  abcd[2220] = 4.0E0*I_NAI_H4xz_Gy3z_ab-2.0E0*2*I_NAI_H4xz_Dyz_a;
  abcd[2221] = 4.0E0*I_NAI_H3xyz_Gy3z_ab-2.0E0*2*I_NAI_H3xyz_Dyz_a;
  abcd[2222] = 4.0E0*I_NAI_H3x2z_Gy3z_ab-2.0E0*2*I_NAI_H3x2z_Dyz_a-2.0E0*1*I_NAI_F3x_Gy3z_b+2*I_NAI_F3x_Dyz;
  abcd[2223] = 4.0E0*I_NAI_H2x2yz_Gy3z_ab-2.0E0*2*I_NAI_H2x2yz_Dyz_a;
  abcd[2224] = 4.0E0*I_NAI_H2xy2z_Gy3z_ab-2.0E0*2*I_NAI_H2xy2z_Dyz_a-2.0E0*1*I_NAI_F2xy_Gy3z_b+2*I_NAI_F2xy_Dyz;
  abcd[2225] = 4.0E0*I_NAI_H2x3z_Gy3z_ab-2.0E0*2*I_NAI_H2x3z_Dyz_a-2.0E0*2*I_NAI_F2xz_Gy3z_b+2*2*I_NAI_F2xz_Dyz;
  abcd[2226] = 4.0E0*I_NAI_Hx3yz_Gy3z_ab-2.0E0*2*I_NAI_Hx3yz_Dyz_a;
  abcd[2227] = 4.0E0*I_NAI_Hx2y2z_Gy3z_ab-2.0E0*2*I_NAI_Hx2y2z_Dyz_a-2.0E0*1*I_NAI_Fx2y_Gy3z_b+2*I_NAI_Fx2y_Dyz;
  abcd[2228] = 4.0E0*I_NAI_Hxy3z_Gy3z_ab-2.0E0*2*I_NAI_Hxy3z_Dyz_a-2.0E0*2*I_NAI_Fxyz_Gy3z_b+2*2*I_NAI_Fxyz_Dyz;
  abcd[2229] = 4.0E0*I_NAI_Hx4z_Gy3z_ab-2.0E0*2*I_NAI_Hx4z_Dyz_a-2.0E0*3*I_NAI_Fx2z_Gy3z_b+3*2*I_NAI_Fx2z_Dyz;
  abcd[2230] = 4.0E0*I_NAI_H4yz_Gy3z_ab-2.0E0*2*I_NAI_H4yz_Dyz_a;
  abcd[2231] = 4.0E0*I_NAI_H3y2z_Gy3z_ab-2.0E0*2*I_NAI_H3y2z_Dyz_a-2.0E0*1*I_NAI_F3y_Gy3z_b+2*I_NAI_F3y_Dyz;
  abcd[2232] = 4.0E0*I_NAI_H2y3z_Gy3z_ab-2.0E0*2*I_NAI_H2y3z_Dyz_a-2.0E0*2*I_NAI_F2yz_Gy3z_b+2*2*I_NAI_F2yz_Dyz;
  abcd[2233] = 4.0E0*I_NAI_Hy4z_Gy3z_ab-2.0E0*2*I_NAI_Hy4z_Dyz_a-2.0E0*3*I_NAI_Fy2z_Gy3z_b+3*2*I_NAI_Fy2z_Dyz;
  abcd[2234] = 4.0E0*I_NAI_H5z_Gy3z_ab-2.0E0*2*I_NAI_H5z_Dyz_a-2.0E0*4*I_NAI_F3z_Gy3z_b+4*2*I_NAI_F3z_Dyz;
  abcd[2235] = 4.0E0*I_NAI_H4xz_G4z_ab-2.0E0*3*I_NAI_H4xz_D2z_a;
  abcd[2236] = 4.0E0*I_NAI_H3xyz_G4z_ab-2.0E0*3*I_NAI_H3xyz_D2z_a;
  abcd[2237] = 4.0E0*I_NAI_H3x2z_G4z_ab-2.0E0*3*I_NAI_H3x2z_D2z_a-2.0E0*1*I_NAI_F3x_G4z_b+3*I_NAI_F3x_D2z;
  abcd[2238] = 4.0E0*I_NAI_H2x2yz_G4z_ab-2.0E0*3*I_NAI_H2x2yz_D2z_a;
  abcd[2239] = 4.0E0*I_NAI_H2xy2z_G4z_ab-2.0E0*3*I_NAI_H2xy2z_D2z_a-2.0E0*1*I_NAI_F2xy_G4z_b+3*I_NAI_F2xy_D2z;
  abcd[2240] = 4.0E0*I_NAI_H2x3z_G4z_ab-2.0E0*3*I_NAI_H2x3z_D2z_a-2.0E0*2*I_NAI_F2xz_G4z_b+2*3*I_NAI_F2xz_D2z;
  abcd[2241] = 4.0E0*I_NAI_Hx3yz_G4z_ab-2.0E0*3*I_NAI_Hx3yz_D2z_a;
  abcd[2242] = 4.0E0*I_NAI_Hx2y2z_G4z_ab-2.0E0*3*I_NAI_Hx2y2z_D2z_a-2.0E0*1*I_NAI_Fx2y_G4z_b+3*I_NAI_Fx2y_D2z;
  abcd[2243] = 4.0E0*I_NAI_Hxy3z_G4z_ab-2.0E0*3*I_NAI_Hxy3z_D2z_a-2.0E0*2*I_NAI_Fxyz_G4z_b+2*3*I_NAI_Fxyz_D2z;
  abcd[2244] = 4.0E0*I_NAI_Hx4z_G4z_ab-2.0E0*3*I_NAI_Hx4z_D2z_a-2.0E0*3*I_NAI_Fx2z_G4z_b+3*3*I_NAI_Fx2z_D2z;
  abcd[2245] = 4.0E0*I_NAI_H4yz_G4z_ab-2.0E0*3*I_NAI_H4yz_D2z_a;
  abcd[2246] = 4.0E0*I_NAI_H3y2z_G4z_ab-2.0E0*3*I_NAI_H3y2z_D2z_a-2.0E0*1*I_NAI_F3y_G4z_b+3*I_NAI_F3y_D2z;
  abcd[2247] = 4.0E0*I_NAI_H2y3z_G4z_ab-2.0E0*3*I_NAI_H2y3z_D2z_a-2.0E0*2*I_NAI_F2yz_G4z_b+2*3*I_NAI_F2yz_D2z;
  abcd[2248] = 4.0E0*I_NAI_Hy4z_G4z_ab-2.0E0*3*I_NAI_Hy4z_D2z_a-2.0E0*3*I_NAI_Fy2z_G4z_b+3*3*I_NAI_Fy2z_D2z;
  abcd[2249] = 4.0E0*I_NAI_H5z_G4z_ab-2.0E0*3*I_NAI_H5z_D2z_a-2.0E0*4*I_NAI_F3z_G4z_b+4*3*I_NAI_F3z_D2z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[2250] = 4.0E0*I_NAI_G4x_H5x_bb-2.0E0*3*I_NAI_G4x_F3x_b-2.0E0*4*I_NAI_G4x_F3x_b+3*2*I_NAI_G4x_Px;
  abcd[2251] = 4.0E0*I_NAI_G3xy_H5x_bb-2.0E0*3*I_NAI_G3xy_F3x_b-2.0E0*4*I_NAI_G3xy_F3x_b+3*2*I_NAI_G3xy_Px;
  abcd[2252] = 4.0E0*I_NAI_G3xz_H5x_bb-2.0E0*3*I_NAI_G3xz_F3x_b-2.0E0*4*I_NAI_G3xz_F3x_b+3*2*I_NAI_G3xz_Px;
  abcd[2253] = 4.0E0*I_NAI_G2x2y_H5x_bb-2.0E0*3*I_NAI_G2x2y_F3x_b-2.0E0*4*I_NAI_G2x2y_F3x_b+3*2*I_NAI_G2x2y_Px;
  abcd[2254] = 4.0E0*I_NAI_G2xyz_H5x_bb-2.0E0*3*I_NAI_G2xyz_F3x_b-2.0E0*4*I_NAI_G2xyz_F3x_b+3*2*I_NAI_G2xyz_Px;
  abcd[2255] = 4.0E0*I_NAI_G2x2z_H5x_bb-2.0E0*3*I_NAI_G2x2z_F3x_b-2.0E0*4*I_NAI_G2x2z_F3x_b+3*2*I_NAI_G2x2z_Px;
  abcd[2256] = 4.0E0*I_NAI_Gx3y_H5x_bb-2.0E0*3*I_NAI_Gx3y_F3x_b-2.0E0*4*I_NAI_Gx3y_F3x_b+3*2*I_NAI_Gx3y_Px;
  abcd[2257] = 4.0E0*I_NAI_Gx2yz_H5x_bb-2.0E0*3*I_NAI_Gx2yz_F3x_b-2.0E0*4*I_NAI_Gx2yz_F3x_b+3*2*I_NAI_Gx2yz_Px;
  abcd[2258] = 4.0E0*I_NAI_Gxy2z_H5x_bb-2.0E0*3*I_NAI_Gxy2z_F3x_b-2.0E0*4*I_NAI_Gxy2z_F3x_b+3*2*I_NAI_Gxy2z_Px;
  abcd[2259] = 4.0E0*I_NAI_Gx3z_H5x_bb-2.0E0*3*I_NAI_Gx3z_F3x_b-2.0E0*4*I_NAI_Gx3z_F3x_b+3*2*I_NAI_Gx3z_Px;
  abcd[2260] = 4.0E0*I_NAI_G4y_H5x_bb-2.0E0*3*I_NAI_G4y_F3x_b-2.0E0*4*I_NAI_G4y_F3x_b+3*2*I_NAI_G4y_Px;
  abcd[2261] = 4.0E0*I_NAI_G3yz_H5x_bb-2.0E0*3*I_NAI_G3yz_F3x_b-2.0E0*4*I_NAI_G3yz_F3x_b+3*2*I_NAI_G3yz_Px;
  abcd[2262] = 4.0E0*I_NAI_G2y2z_H5x_bb-2.0E0*3*I_NAI_G2y2z_F3x_b-2.0E0*4*I_NAI_G2y2z_F3x_b+3*2*I_NAI_G2y2z_Px;
  abcd[2263] = 4.0E0*I_NAI_Gy3z_H5x_bb-2.0E0*3*I_NAI_Gy3z_F3x_b-2.0E0*4*I_NAI_Gy3z_F3x_b+3*2*I_NAI_Gy3z_Px;
  abcd[2264] = 4.0E0*I_NAI_G4z_H5x_bb-2.0E0*3*I_NAI_G4z_F3x_b-2.0E0*4*I_NAI_G4z_F3x_b+3*2*I_NAI_G4z_Px;
  abcd[2265] = 4.0E0*I_NAI_G4x_H4xy_bb-2.0E0*2*I_NAI_G4x_F2xy_b-2.0E0*3*I_NAI_G4x_F2xy_b+2*1*I_NAI_G4x_Py;
  abcd[2266] = 4.0E0*I_NAI_G3xy_H4xy_bb-2.0E0*2*I_NAI_G3xy_F2xy_b-2.0E0*3*I_NAI_G3xy_F2xy_b+2*1*I_NAI_G3xy_Py;
  abcd[2267] = 4.0E0*I_NAI_G3xz_H4xy_bb-2.0E0*2*I_NAI_G3xz_F2xy_b-2.0E0*3*I_NAI_G3xz_F2xy_b+2*1*I_NAI_G3xz_Py;
  abcd[2268] = 4.0E0*I_NAI_G2x2y_H4xy_bb-2.0E0*2*I_NAI_G2x2y_F2xy_b-2.0E0*3*I_NAI_G2x2y_F2xy_b+2*1*I_NAI_G2x2y_Py;
  abcd[2269] = 4.0E0*I_NAI_G2xyz_H4xy_bb-2.0E0*2*I_NAI_G2xyz_F2xy_b-2.0E0*3*I_NAI_G2xyz_F2xy_b+2*1*I_NAI_G2xyz_Py;
  abcd[2270] = 4.0E0*I_NAI_G2x2z_H4xy_bb-2.0E0*2*I_NAI_G2x2z_F2xy_b-2.0E0*3*I_NAI_G2x2z_F2xy_b+2*1*I_NAI_G2x2z_Py;
  abcd[2271] = 4.0E0*I_NAI_Gx3y_H4xy_bb-2.0E0*2*I_NAI_Gx3y_F2xy_b-2.0E0*3*I_NAI_Gx3y_F2xy_b+2*1*I_NAI_Gx3y_Py;
  abcd[2272] = 4.0E0*I_NAI_Gx2yz_H4xy_bb-2.0E0*2*I_NAI_Gx2yz_F2xy_b-2.0E0*3*I_NAI_Gx2yz_F2xy_b+2*1*I_NAI_Gx2yz_Py;
  abcd[2273] = 4.0E0*I_NAI_Gxy2z_H4xy_bb-2.0E0*2*I_NAI_Gxy2z_F2xy_b-2.0E0*3*I_NAI_Gxy2z_F2xy_b+2*1*I_NAI_Gxy2z_Py;
  abcd[2274] = 4.0E0*I_NAI_Gx3z_H4xy_bb-2.0E0*2*I_NAI_Gx3z_F2xy_b-2.0E0*3*I_NAI_Gx3z_F2xy_b+2*1*I_NAI_Gx3z_Py;
  abcd[2275] = 4.0E0*I_NAI_G4y_H4xy_bb-2.0E0*2*I_NAI_G4y_F2xy_b-2.0E0*3*I_NAI_G4y_F2xy_b+2*1*I_NAI_G4y_Py;
  abcd[2276] = 4.0E0*I_NAI_G3yz_H4xy_bb-2.0E0*2*I_NAI_G3yz_F2xy_b-2.0E0*3*I_NAI_G3yz_F2xy_b+2*1*I_NAI_G3yz_Py;
  abcd[2277] = 4.0E0*I_NAI_G2y2z_H4xy_bb-2.0E0*2*I_NAI_G2y2z_F2xy_b-2.0E0*3*I_NAI_G2y2z_F2xy_b+2*1*I_NAI_G2y2z_Py;
  abcd[2278] = 4.0E0*I_NAI_Gy3z_H4xy_bb-2.0E0*2*I_NAI_Gy3z_F2xy_b-2.0E0*3*I_NAI_Gy3z_F2xy_b+2*1*I_NAI_Gy3z_Py;
  abcd[2279] = 4.0E0*I_NAI_G4z_H4xy_bb-2.0E0*2*I_NAI_G4z_F2xy_b-2.0E0*3*I_NAI_G4z_F2xy_b+2*1*I_NAI_G4z_Py;
  abcd[2280] = 4.0E0*I_NAI_G4x_H4xz_bb-2.0E0*2*I_NAI_G4x_F2xz_b-2.0E0*3*I_NAI_G4x_F2xz_b+2*1*I_NAI_G4x_Pz;
  abcd[2281] = 4.0E0*I_NAI_G3xy_H4xz_bb-2.0E0*2*I_NAI_G3xy_F2xz_b-2.0E0*3*I_NAI_G3xy_F2xz_b+2*1*I_NAI_G3xy_Pz;
  abcd[2282] = 4.0E0*I_NAI_G3xz_H4xz_bb-2.0E0*2*I_NAI_G3xz_F2xz_b-2.0E0*3*I_NAI_G3xz_F2xz_b+2*1*I_NAI_G3xz_Pz;
  abcd[2283] = 4.0E0*I_NAI_G2x2y_H4xz_bb-2.0E0*2*I_NAI_G2x2y_F2xz_b-2.0E0*3*I_NAI_G2x2y_F2xz_b+2*1*I_NAI_G2x2y_Pz;
  abcd[2284] = 4.0E0*I_NAI_G2xyz_H4xz_bb-2.0E0*2*I_NAI_G2xyz_F2xz_b-2.0E0*3*I_NAI_G2xyz_F2xz_b+2*1*I_NAI_G2xyz_Pz;
  abcd[2285] = 4.0E0*I_NAI_G2x2z_H4xz_bb-2.0E0*2*I_NAI_G2x2z_F2xz_b-2.0E0*3*I_NAI_G2x2z_F2xz_b+2*1*I_NAI_G2x2z_Pz;
  abcd[2286] = 4.0E0*I_NAI_Gx3y_H4xz_bb-2.0E0*2*I_NAI_Gx3y_F2xz_b-2.0E0*3*I_NAI_Gx3y_F2xz_b+2*1*I_NAI_Gx3y_Pz;
  abcd[2287] = 4.0E0*I_NAI_Gx2yz_H4xz_bb-2.0E0*2*I_NAI_Gx2yz_F2xz_b-2.0E0*3*I_NAI_Gx2yz_F2xz_b+2*1*I_NAI_Gx2yz_Pz;
  abcd[2288] = 4.0E0*I_NAI_Gxy2z_H4xz_bb-2.0E0*2*I_NAI_Gxy2z_F2xz_b-2.0E0*3*I_NAI_Gxy2z_F2xz_b+2*1*I_NAI_Gxy2z_Pz;
  abcd[2289] = 4.0E0*I_NAI_Gx3z_H4xz_bb-2.0E0*2*I_NAI_Gx3z_F2xz_b-2.0E0*3*I_NAI_Gx3z_F2xz_b+2*1*I_NAI_Gx3z_Pz;
  abcd[2290] = 4.0E0*I_NAI_G4y_H4xz_bb-2.0E0*2*I_NAI_G4y_F2xz_b-2.0E0*3*I_NAI_G4y_F2xz_b+2*1*I_NAI_G4y_Pz;
  abcd[2291] = 4.0E0*I_NAI_G3yz_H4xz_bb-2.0E0*2*I_NAI_G3yz_F2xz_b-2.0E0*3*I_NAI_G3yz_F2xz_b+2*1*I_NAI_G3yz_Pz;
  abcd[2292] = 4.0E0*I_NAI_G2y2z_H4xz_bb-2.0E0*2*I_NAI_G2y2z_F2xz_b-2.0E0*3*I_NAI_G2y2z_F2xz_b+2*1*I_NAI_G2y2z_Pz;
  abcd[2293] = 4.0E0*I_NAI_Gy3z_H4xz_bb-2.0E0*2*I_NAI_Gy3z_F2xz_b-2.0E0*3*I_NAI_Gy3z_F2xz_b+2*1*I_NAI_Gy3z_Pz;
  abcd[2294] = 4.0E0*I_NAI_G4z_H4xz_bb-2.0E0*2*I_NAI_G4z_F2xz_b-2.0E0*3*I_NAI_G4z_F2xz_b+2*1*I_NAI_G4z_Pz;
  abcd[2295] = 4.0E0*I_NAI_G4x_H3x2y_bb-2.0E0*1*I_NAI_G4x_Fx2y_b-2.0E0*2*I_NAI_G4x_Fx2y_b;
  abcd[2296] = 4.0E0*I_NAI_G3xy_H3x2y_bb-2.0E0*1*I_NAI_G3xy_Fx2y_b-2.0E0*2*I_NAI_G3xy_Fx2y_b;
  abcd[2297] = 4.0E0*I_NAI_G3xz_H3x2y_bb-2.0E0*1*I_NAI_G3xz_Fx2y_b-2.0E0*2*I_NAI_G3xz_Fx2y_b;
  abcd[2298] = 4.0E0*I_NAI_G2x2y_H3x2y_bb-2.0E0*1*I_NAI_G2x2y_Fx2y_b-2.0E0*2*I_NAI_G2x2y_Fx2y_b;
  abcd[2299] = 4.0E0*I_NAI_G2xyz_H3x2y_bb-2.0E0*1*I_NAI_G2xyz_Fx2y_b-2.0E0*2*I_NAI_G2xyz_Fx2y_b;
  abcd[2300] = 4.0E0*I_NAI_G2x2z_H3x2y_bb-2.0E0*1*I_NAI_G2x2z_Fx2y_b-2.0E0*2*I_NAI_G2x2z_Fx2y_b;
  abcd[2301] = 4.0E0*I_NAI_Gx3y_H3x2y_bb-2.0E0*1*I_NAI_Gx3y_Fx2y_b-2.0E0*2*I_NAI_Gx3y_Fx2y_b;
  abcd[2302] = 4.0E0*I_NAI_Gx2yz_H3x2y_bb-2.0E0*1*I_NAI_Gx2yz_Fx2y_b-2.0E0*2*I_NAI_Gx2yz_Fx2y_b;
  abcd[2303] = 4.0E0*I_NAI_Gxy2z_H3x2y_bb-2.0E0*1*I_NAI_Gxy2z_Fx2y_b-2.0E0*2*I_NAI_Gxy2z_Fx2y_b;
  abcd[2304] = 4.0E0*I_NAI_Gx3z_H3x2y_bb-2.0E0*1*I_NAI_Gx3z_Fx2y_b-2.0E0*2*I_NAI_Gx3z_Fx2y_b;
  abcd[2305] = 4.0E0*I_NAI_G4y_H3x2y_bb-2.0E0*1*I_NAI_G4y_Fx2y_b-2.0E0*2*I_NAI_G4y_Fx2y_b;
  abcd[2306] = 4.0E0*I_NAI_G3yz_H3x2y_bb-2.0E0*1*I_NAI_G3yz_Fx2y_b-2.0E0*2*I_NAI_G3yz_Fx2y_b;
  abcd[2307] = 4.0E0*I_NAI_G2y2z_H3x2y_bb-2.0E0*1*I_NAI_G2y2z_Fx2y_b-2.0E0*2*I_NAI_G2y2z_Fx2y_b;
  abcd[2308] = 4.0E0*I_NAI_Gy3z_H3x2y_bb-2.0E0*1*I_NAI_Gy3z_Fx2y_b-2.0E0*2*I_NAI_Gy3z_Fx2y_b;
  abcd[2309] = 4.0E0*I_NAI_G4z_H3x2y_bb-2.0E0*1*I_NAI_G4z_Fx2y_b-2.0E0*2*I_NAI_G4z_Fx2y_b;
  abcd[2310] = 4.0E0*I_NAI_G4x_H3xyz_bb-2.0E0*1*I_NAI_G4x_Fxyz_b-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2311] = 4.0E0*I_NAI_G3xy_H3xyz_bb-2.0E0*1*I_NAI_G3xy_Fxyz_b-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2312] = 4.0E0*I_NAI_G3xz_H3xyz_bb-2.0E0*1*I_NAI_G3xz_Fxyz_b-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2313] = 4.0E0*I_NAI_G2x2y_H3xyz_bb-2.0E0*1*I_NAI_G2x2y_Fxyz_b-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2314] = 4.0E0*I_NAI_G2xyz_H3xyz_bb-2.0E0*1*I_NAI_G2xyz_Fxyz_b-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2315] = 4.0E0*I_NAI_G2x2z_H3xyz_bb-2.0E0*1*I_NAI_G2x2z_Fxyz_b-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2316] = 4.0E0*I_NAI_Gx3y_H3xyz_bb-2.0E0*1*I_NAI_Gx3y_Fxyz_b-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2317] = 4.0E0*I_NAI_Gx2yz_H3xyz_bb-2.0E0*1*I_NAI_Gx2yz_Fxyz_b-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2318] = 4.0E0*I_NAI_Gxy2z_H3xyz_bb-2.0E0*1*I_NAI_Gxy2z_Fxyz_b-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2319] = 4.0E0*I_NAI_Gx3z_H3xyz_bb-2.0E0*1*I_NAI_Gx3z_Fxyz_b-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2320] = 4.0E0*I_NAI_G4y_H3xyz_bb-2.0E0*1*I_NAI_G4y_Fxyz_b-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2321] = 4.0E0*I_NAI_G3yz_H3xyz_bb-2.0E0*1*I_NAI_G3yz_Fxyz_b-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2322] = 4.0E0*I_NAI_G2y2z_H3xyz_bb-2.0E0*1*I_NAI_G2y2z_Fxyz_b-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2323] = 4.0E0*I_NAI_Gy3z_H3xyz_bb-2.0E0*1*I_NAI_Gy3z_Fxyz_b-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2324] = 4.0E0*I_NAI_G4z_H3xyz_bb-2.0E0*1*I_NAI_G4z_Fxyz_b-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2325] = 4.0E0*I_NAI_G4x_H3x2z_bb-2.0E0*1*I_NAI_G4x_Fx2z_b-2.0E0*2*I_NAI_G4x_Fx2z_b;
  abcd[2326] = 4.0E0*I_NAI_G3xy_H3x2z_bb-2.0E0*1*I_NAI_G3xy_Fx2z_b-2.0E0*2*I_NAI_G3xy_Fx2z_b;
  abcd[2327] = 4.0E0*I_NAI_G3xz_H3x2z_bb-2.0E0*1*I_NAI_G3xz_Fx2z_b-2.0E0*2*I_NAI_G3xz_Fx2z_b;
  abcd[2328] = 4.0E0*I_NAI_G2x2y_H3x2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2z_b-2.0E0*2*I_NAI_G2x2y_Fx2z_b;
  abcd[2329] = 4.0E0*I_NAI_G2xyz_H3x2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2z_b-2.0E0*2*I_NAI_G2xyz_Fx2z_b;
  abcd[2330] = 4.0E0*I_NAI_G2x2z_H3x2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2z_b-2.0E0*2*I_NAI_G2x2z_Fx2z_b;
  abcd[2331] = 4.0E0*I_NAI_Gx3y_H3x2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2z_b-2.0E0*2*I_NAI_Gx3y_Fx2z_b;
  abcd[2332] = 4.0E0*I_NAI_Gx2yz_H3x2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2z_b-2.0E0*2*I_NAI_Gx2yz_Fx2z_b;
  abcd[2333] = 4.0E0*I_NAI_Gxy2z_H3x2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2z_b-2.0E0*2*I_NAI_Gxy2z_Fx2z_b;
  abcd[2334] = 4.0E0*I_NAI_Gx3z_H3x2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2z_b-2.0E0*2*I_NAI_Gx3z_Fx2z_b;
  abcd[2335] = 4.0E0*I_NAI_G4y_H3x2z_bb-2.0E0*1*I_NAI_G4y_Fx2z_b-2.0E0*2*I_NAI_G4y_Fx2z_b;
  abcd[2336] = 4.0E0*I_NAI_G3yz_H3x2z_bb-2.0E0*1*I_NAI_G3yz_Fx2z_b-2.0E0*2*I_NAI_G3yz_Fx2z_b;
  abcd[2337] = 4.0E0*I_NAI_G2y2z_H3x2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2z_b-2.0E0*2*I_NAI_G2y2z_Fx2z_b;
  abcd[2338] = 4.0E0*I_NAI_Gy3z_H3x2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2z_b-2.0E0*2*I_NAI_Gy3z_Fx2z_b;
  abcd[2339] = 4.0E0*I_NAI_G4z_H3x2z_bb-2.0E0*1*I_NAI_G4z_Fx2z_b-2.0E0*2*I_NAI_G4z_Fx2z_b;
  abcd[2340] = 4.0E0*I_NAI_G4x_H2x3y_bb-2.0E0*1*I_NAI_G4x_F3y_b;
  abcd[2341] = 4.0E0*I_NAI_G3xy_H2x3y_bb-2.0E0*1*I_NAI_G3xy_F3y_b;
  abcd[2342] = 4.0E0*I_NAI_G3xz_H2x3y_bb-2.0E0*1*I_NAI_G3xz_F3y_b;
  abcd[2343] = 4.0E0*I_NAI_G2x2y_H2x3y_bb-2.0E0*1*I_NAI_G2x2y_F3y_b;
  abcd[2344] = 4.0E0*I_NAI_G2xyz_H2x3y_bb-2.0E0*1*I_NAI_G2xyz_F3y_b;
  abcd[2345] = 4.0E0*I_NAI_G2x2z_H2x3y_bb-2.0E0*1*I_NAI_G2x2z_F3y_b;
  abcd[2346] = 4.0E0*I_NAI_Gx3y_H2x3y_bb-2.0E0*1*I_NAI_Gx3y_F3y_b;
  abcd[2347] = 4.0E0*I_NAI_Gx2yz_H2x3y_bb-2.0E0*1*I_NAI_Gx2yz_F3y_b;
  abcd[2348] = 4.0E0*I_NAI_Gxy2z_H2x3y_bb-2.0E0*1*I_NAI_Gxy2z_F3y_b;
  abcd[2349] = 4.0E0*I_NAI_Gx3z_H2x3y_bb-2.0E0*1*I_NAI_Gx3z_F3y_b;
  abcd[2350] = 4.0E0*I_NAI_G4y_H2x3y_bb-2.0E0*1*I_NAI_G4y_F3y_b;
  abcd[2351] = 4.0E0*I_NAI_G3yz_H2x3y_bb-2.0E0*1*I_NAI_G3yz_F3y_b;
  abcd[2352] = 4.0E0*I_NAI_G2y2z_H2x3y_bb-2.0E0*1*I_NAI_G2y2z_F3y_b;
  abcd[2353] = 4.0E0*I_NAI_Gy3z_H2x3y_bb-2.0E0*1*I_NAI_Gy3z_F3y_b;
  abcd[2354] = 4.0E0*I_NAI_G4z_H2x3y_bb-2.0E0*1*I_NAI_G4z_F3y_b;
  abcd[2355] = 4.0E0*I_NAI_G4x_H2x2yz_bb-2.0E0*1*I_NAI_G4x_F2yz_b;
  abcd[2356] = 4.0E0*I_NAI_G3xy_H2x2yz_bb-2.0E0*1*I_NAI_G3xy_F2yz_b;
  abcd[2357] = 4.0E0*I_NAI_G3xz_H2x2yz_bb-2.0E0*1*I_NAI_G3xz_F2yz_b;
  abcd[2358] = 4.0E0*I_NAI_G2x2y_H2x2yz_bb-2.0E0*1*I_NAI_G2x2y_F2yz_b;
  abcd[2359] = 4.0E0*I_NAI_G2xyz_H2x2yz_bb-2.0E0*1*I_NAI_G2xyz_F2yz_b;
  abcd[2360] = 4.0E0*I_NAI_G2x2z_H2x2yz_bb-2.0E0*1*I_NAI_G2x2z_F2yz_b;
  abcd[2361] = 4.0E0*I_NAI_Gx3y_H2x2yz_bb-2.0E0*1*I_NAI_Gx3y_F2yz_b;
  abcd[2362] = 4.0E0*I_NAI_Gx2yz_H2x2yz_bb-2.0E0*1*I_NAI_Gx2yz_F2yz_b;
  abcd[2363] = 4.0E0*I_NAI_Gxy2z_H2x2yz_bb-2.0E0*1*I_NAI_Gxy2z_F2yz_b;
  abcd[2364] = 4.0E0*I_NAI_Gx3z_H2x2yz_bb-2.0E0*1*I_NAI_Gx3z_F2yz_b;
  abcd[2365] = 4.0E0*I_NAI_G4y_H2x2yz_bb-2.0E0*1*I_NAI_G4y_F2yz_b;
  abcd[2366] = 4.0E0*I_NAI_G3yz_H2x2yz_bb-2.0E0*1*I_NAI_G3yz_F2yz_b;
  abcd[2367] = 4.0E0*I_NAI_G2y2z_H2x2yz_bb-2.0E0*1*I_NAI_G2y2z_F2yz_b;
  abcd[2368] = 4.0E0*I_NAI_Gy3z_H2x2yz_bb-2.0E0*1*I_NAI_Gy3z_F2yz_b;
  abcd[2369] = 4.0E0*I_NAI_G4z_H2x2yz_bb-2.0E0*1*I_NAI_G4z_F2yz_b;
  abcd[2370] = 4.0E0*I_NAI_G4x_H2xy2z_bb-2.0E0*1*I_NAI_G4x_Fy2z_b;
  abcd[2371] = 4.0E0*I_NAI_G3xy_H2xy2z_bb-2.0E0*1*I_NAI_G3xy_Fy2z_b;
  abcd[2372] = 4.0E0*I_NAI_G3xz_H2xy2z_bb-2.0E0*1*I_NAI_G3xz_Fy2z_b;
  abcd[2373] = 4.0E0*I_NAI_G2x2y_H2xy2z_bb-2.0E0*1*I_NAI_G2x2y_Fy2z_b;
  abcd[2374] = 4.0E0*I_NAI_G2xyz_H2xy2z_bb-2.0E0*1*I_NAI_G2xyz_Fy2z_b;
  abcd[2375] = 4.0E0*I_NAI_G2x2z_H2xy2z_bb-2.0E0*1*I_NAI_G2x2z_Fy2z_b;
  abcd[2376] = 4.0E0*I_NAI_Gx3y_H2xy2z_bb-2.0E0*1*I_NAI_Gx3y_Fy2z_b;
  abcd[2377] = 4.0E0*I_NAI_Gx2yz_H2xy2z_bb-2.0E0*1*I_NAI_Gx2yz_Fy2z_b;
  abcd[2378] = 4.0E0*I_NAI_Gxy2z_H2xy2z_bb-2.0E0*1*I_NAI_Gxy2z_Fy2z_b;
  abcd[2379] = 4.0E0*I_NAI_Gx3z_H2xy2z_bb-2.0E0*1*I_NAI_Gx3z_Fy2z_b;
  abcd[2380] = 4.0E0*I_NAI_G4y_H2xy2z_bb-2.0E0*1*I_NAI_G4y_Fy2z_b;
  abcd[2381] = 4.0E0*I_NAI_G3yz_H2xy2z_bb-2.0E0*1*I_NAI_G3yz_Fy2z_b;
  abcd[2382] = 4.0E0*I_NAI_G2y2z_H2xy2z_bb-2.0E0*1*I_NAI_G2y2z_Fy2z_b;
  abcd[2383] = 4.0E0*I_NAI_Gy3z_H2xy2z_bb-2.0E0*1*I_NAI_Gy3z_Fy2z_b;
  abcd[2384] = 4.0E0*I_NAI_G4z_H2xy2z_bb-2.0E0*1*I_NAI_G4z_Fy2z_b;
  abcd[2385] = 4.0E0*I_NAI_G4x_H2x3z_bb-2.0E0*1*I_NAI_G4x_F3z_b;
  abcd[2386] = 4.0E0*I_NAI_G3xy_H2x3z_bb-2.0E0*1*I_NAI_G3xy_F3z_b;
  abcd[2387] = 4.0E0*I_NAI_G3xz_H2x3z_bb-2.0E0*1*I_NAI_G3xz_F3z_b;
  abcd[2388] = 4.0E0*I_NAI_G2x2y_H2x3z_bb-2.0E0*1*I_NAI_G2x2y_F3z_b;
  abcd[2389] = 4.0E0*I_NAI_G2xyz_H2x3z_bb-2.0E0*1*I_NAI_G2xyz_F3z_b;
  abcd[2390] = 4.0E0*I_NAI_G2x2z_H2x3z_bb-2.0E0*1*I_NAI_G2x2z_F3z_b;
  abcd[2391] = 4.0E0*I_NAI_Gx3y_H2x3z_bb-2.0E0*1*I_NAI_Gx3y_F3z_b;
  abcd[2392] = 4.0E0*I_NAI_Gx2yz_H2x3z_bb-2.0E0*1*I_NAI_Gx2yz_F3z_b;
  abcd[2393] = 4.0E0*I_NAI_Gxy2z_H2x3z_bb-2.0E0*1*I_NAI_Gxy2z_F3z_b;
  abcd[2394] = 4.0E0*I_NAI_Gx3z_H2x3z_bb-2.0E0*1*I_NAI_Gx3z_F3z_b;
  abcd[2395] = 4.0E0*I_NAI_G4y_H2x3z_bb-2.0E0*1*I_NAI_G4y_F3z_b;
  abcd[2396] = 4.0E0*I_NAI_G3yz_H2x3z_bb-2.0E0*1*I_NAI_G3yz_F3z_b;
  abcd[2397] = 4.0E0*I_NAI_G2y2z_H2x3z_bb-2.0E0*1*I_NAI_G2y2z_F3z_b;
  abcd[2398] = 4.0E0*I_NAI_Gy3z_H2x3z_bb-2.0E0*1*I_NAI_Gy3z_F3z_b;
  abcd[2399] = 4.0E0*I_NAI_G4z_H2x3z_bb-2.0E0*1*I_NAI_G4z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[2400] = 4.0E0*I_NAI_G4x_H4xy_bb-2.0E0*3*I_NAI_G4x_F2xy_b;
  abcd[2401] = 4.0E0*I_NAI_G3xy_H4xy_bb-2.0E0*3*I_NAI_G3xy_F2xy_b;
  abcd[2402] = 4.0E0*I_NAI_G3xz_H4xy_bb-2.0E0*3*I_NAI_G3xz_F2xy_b;
  abcd[2403] = 4.0E0*I_NAI_G2x2y_H4xy_bb-2.0E0*3*I_NAI_G2x2y_F2xy_b;
  abcd[2404] = 4.0E0*I_NAI_G2xyz_H4xy_bb-2.0E0*3*I_NAI_G2xyz_F2xy_b;
  abcd[2405] = 4.0E0*I_NAI_G2x2z_H4xy_bb-2.0E0*3*I_NAI_G2x2z_F2xy_b;
  abcd[2406] = 4.0E0*I_NAI_Gx3y_H4xy_bb-2.0E0*3*I_NAI_Gx3y_F2xy_b;
  abcd[2407] = 4.0E0*I_NAI_Gx2yz_H4xy_bb-2.0E0*3*I_NAI_Gx2yz_F2xy_b;
  abcd[2408] = 4.0E0*I_NAI_Gxy2z_H4xy_bb-2.0E0*3*I_NAI_Gxy2z_F2xy_b;
  abcd[2409] = 4.0E0*I_NAI_Gx3z_H4xy_bb-2.0E0*3*I_NAI_Gx3z_F2xy_b;
  abcd[2410] = 4.0E0*I_NAI_G4y_H4xy_bb-2.0E0*3*I_NAI_G4y_F2xy_b;
  abcd[2411] = 4.0E0*I_NAI_G3yz_H4xy_bb-2.0E0*3*I_NAI_G3yz_F2xy_b;
  abcd[2412] = 4.0E0*I_NAI_G2y2z_H4xy_bb-2.0E0*3*I_NAI_G2y2z_F2xy_b;
  abcd[2413] = 4.0E0*I_NAI_Gy3z_H4xy_bb-2.0E0*3*I_NAI_Gy3z_F2xy_b;
  abcd[2414] = 4.0E0*I_NAI_G4z_H4xy_bb-2.0E0*3*I_NAI_G4z_F2xy_b;
  abcd[2415] = 4.0E0*I_NAI_G4x_H3x2y_bb-2.0E0*1*I_NAI_G4x_F3x_b-2.0E0*2*I_NAI_G4x_Fx2y_b+2*1*I_NAI_G4x_Px;
  abcd[2416] = 4.0E0*I_NAI_G3xy_H3x2y_bb-2.0E0*1*I_NAI_G3xy_F3x_b-2.0E0*2*I_NAI_G3xy_Fx2y_b+2*1*I_NAI_G3xy_Px;
  abcd[2417] = 4.0E0*I_NAI_G3xz_H3x2y_bb-2.0E0*1*I_NAI_G3xz_F3x_b-2.0E0*2*I_NAI_G3xz_Fx2y_b+2*1*I_NAI_G3xz_Px;
  abcd[2418] = 4.0E0*I_NAI_G2x2y_H3x2y_bb-2.0E0*1*I_NAI_G2x2y_F3x_b-2.0E0*2*I_NAI_G2x2y_Fx2y_b+2*1*I_NAI_G2x2y_Px;
  abcd[2419] = 4.0E0*I_NAI_G2xyz_H3x2y_bb-2.0E0*1*I_NAI_G2xyz_F3x_b-2.0E0*2*I_NAI_G2xyz_Fx2y_b+2*1*I_NAI_G2xyz_Px;
  abcd[2420] = 4.0E0*I_NAI_G2x2z_H3x2y_bb-2.0E0*1*I_NAI_G2x2z_F3x_b-2.0E0*2*I_NAI_G2x2z_Fx2y_b+2*1*I_NAI_G2x2z_Px;
  abcd[2421] = 4.0E0*I_NAI_Gx3y_H3x2y_bb-2.0E0*1*I_NAI_Gx3y_F3x_b-2.0E0*2*I_NAI_Gx3y_Fx2y_b+2*1*I_NAI_Gx3y_Px;
  abcd[2422] = 4.0E0*I_NAI_Gx2yz_H3x2y_bb-2.0E0*1*I_NAI_Gx2yz_F3x_b-2.0E0*2*I_NAI_Gx2yz_Fx2y_b+2*1*I_NAI_Gx2yz_Px;
  abcd[2423] = 4.0E0*I_NAI_Gxy2z_H3x2y_bb-2.0E0*1*I_NAI_Gxy2z_F3x_b-2.0E0*2*I_NAI_Gxy2z_Fx2y_b+2*1*I_NAI_Gxy2z_Px;
  abcd[2424] = 4.0E0*I_NAI_Gx3z_H3x2y_bb-2.0E0*1*I_NAI_Gx3z_F3x_b-2.0E0*2*I_NAI_Gx3z_Fx2y_b+2*1*I_NAI_Gx3z_Px;
  abcd[2425] = 4.0E0*I_NAI_G4y_H3x2y_bb-2.0E0*1*I_NAI_G4y_F3x_b-2.0E0*2*I_NAI_G4y_Fx2y_b+2*1*I_NAI_G4y_Px;
  abcd[2426] = 4.0E0*I_NAI_G3yz_H3x2y_bb-2.0E0*1*I_NAI_G3yz_F3x_b-2.0E0*2*I_NAI_G3yz_Fx2y_b+2*1*I_NAI_G3yz_Px;
  abcd[2427] = 4.0E0*I_NAI_G2y2z_H3x2y_bb-2.0E0*1*I_NAI_G2y2z_F3x_b-2.0E0*2*I_NAI_G2y2z_Fx2y_b+2*1*I_NAI_G2y2z_Px;
  abcd[2428] = 4.0E0*I_NAI_Gy3z_H3x2y_bb-2.0E0*1*I_NAI_Gy3z_F3x_b-2.0E0*2*I_NAI_Gy3z_Fx2y_b+2*1*I_NAI_Gy3z_Px;
  abcd[2429] = 4.0E0*I_NAI_G4z_H3x2y_bb-2.0E0*1*I_NAI_G4z_F3x_b-2.0E0*2*I_NAI_G4z_Fx2y_b+2*1*I_NAI_G4z_Px;
  abcd[2430] = 4.0E0*I_NAI_G4x_H3xyz_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2431] = 4.0E0*I_NAI_G3xy_H3xyz_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2432] = 4.0E0*I_NAI_G3xz_H3xyz_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2433] = 4.0E0*I_NAI_G2x2y_H3xyz_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2434] = 4.0E0*I_NAI_G2xyz_H3xyz_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2435] = 4.0E0*I_NAI_G2x2z_H3xyz_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2436] = 4.0E0*I_NAI_Gx3y_H3xyz_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2437] = 4.0E0*I_NAI_Gx2yz_H3xyz_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2438] = 4.0E0*I_NAI_Gxy2z_H3xyz_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2439] = 4.0E0*I_NAI_Gx3z_H3xyz_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2440] = 4.0E0*I_NAI_G4y_H3xyz_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2441] = 4.0E0*I_NAI_G3yz_H3xyz_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2442] = 4.0E0*I_NAI_G2y2z_H3xyz_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2443] = 4.0E0*I_NAI_Gy3z_H3xyz_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2444] = 4.0E0*I_NAI_G4z_H3xyz_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2445] = 4.0E0*I_NAI_G4x_H2x3y_bb-2.0E0*2*I_NAI_G4x_F2xy_b-2.0E0*1*I_NAI_G4x_F3y_b+2*I_NAI_G4x_Py;
  abcd[2446] = 4.0E0*I_NAI_G3xy_H2x3y_bb-2.0E0*2*I_NAI_G3xy_F2xy_b-2.0E0*1*I_NAI_G3xy_F3y_b+2*I_NAI_G3xy_Py;
  abcd[2447] = 4.0E0*I_NAI_G3xz_H2x3y_bb-2.0E0*2*I_NAI_G3xz_F2xy_b-2.0E0*1*I_NAI_G3xz_F3y_b+2*I_NAI_G3xz_Py;
  abcd[2448] = 4.0E0*I_NAI_G2x2y_H2x3y_bb-2.0E0*2*I_NAI_G2x2y_F2xy_b-2.0E0*1*I_NAI_G2x2y_F3y_b+2*I_NAI_G2x2y_Py;
  abcd[2449] = 4.0E0*I_NAI_G2xyz_H2x3y_bb-2.0E0*2*I_NAI_G2xyz_F2xy_b-2.0E0*1*I_NAI_G2xyz_F3y_b+2*I_NAI_G2xyz_Py;
  abcd[2450] = 4.0E0*I_NAI_G2x2z_H2x3y_bb-2.0E0*2*I_NAI_G2x2z_F2xy_b-2.0E0*1*I_NAI_G2x2z_F3y_b+2*I_NAI_G2x2z_Py;
  abcd[2451] = 4.0E0*I_NAI_Gx3y_H2x3y_bb-2.0E0*2*I_NAI_Gx3y_F2xy_b-2.0E0*1*I_NAI_Gx3y_F3y_b+2*I_NAI_Gx3y_Py;
  abcd[2452] = 4.0E0*I_NAI_Gx2yz_H2x3y_bb-2.0E0*2*I_NAI_Gx2yz_F2xy_b-2.0E0*1*I_NAI_Gx2yz_F3y_b+2*I_NAI_Gx2yz_Py;
  abcd[2453] = 4.0E0*I_NAI_Gxy2z_H2x3y_bb-2.0E0*2*I_NAI_Gxy2z_F2xy_b-2.0E0*1*I_NAI_Gxy2z_F3y_b+2*I_NAI_Gxy2z_Py;
  abcd[2454] = 4.0E0*I_NAI_Gx3z_H2x3y_bb-2.0E0*2*I_NAI_Gx3z_F2xy_b-2.0E0*1*I_NAI_Gx3z_F3y_b+2*I_NAI_Gx3z_Py;
  abcd[2455] = 4.0E0*I_NAI_G4y_H2x3y_bb-2.0E0*2*I_NAI_G4y_F2xy_b-2.0E0*1*I_NAI_G4y_F3y_b+2*I_NAI_G4y_Py;
  abcd[2456] = 4.0E0*I_NAI_G3yz_H2x3y_bb-2.0E0*2*I_NAI_G3yz_F2xy_b-2.0E0*1*I_NAI_G3yz_F3y_b+2*I_NAI_G3yz_Py;
  abcd[2457] = 4.0E0*I_NAI_G2y2z_H2x3y_bb-2.0E0*2*I_NAI_G2y2z_F2xy_b-2.0E0*1*I_NAI_G2y2z_F3y_b+2*I_NAI_G2y2z_Py;
  abcd[2458] = 4.0E0*I_NAI_Gy3z_H2x3y_bb-2.0E0*2*I_NAI_Gy3z_F2xy_b-2.0E0*1*I_NAI_Gy3z_F3y_b+2*I_NAI_Gy3z_Py;
  abcd[2459] = 4.0E0*I_NAI_G4z_H2x3y_bb-2.0E0*2*I_NAI_G4z_F2xy_b-2.0E0*1*I_NAI_G4z_F3y_b+2*I_NAI_G4z_Py;
  abcd[2460] = 4.0E0*I_NAI_G4x_H2x2yz_bb-2.0E0*1*I_NAI_G4x_F2xz_b-2.0E0*1*I_NAI_G4x_F2yz_b+1*I_NAI_G4x_Pz;
  abcd[2461] = 4.0E0*I_NAI_G3xy_H2x2yz_bb-2.0E0*1*I_NAI_G3xy_F2xz_b-2.0E0*1*I_NAI_G3xy_F2yz_b+1*I_NAI_G3xy_Pz;
  abcd[2462] = 4.0E0*I_NAI_G3xz_H2x2yz_bb-2.0E0*1*I_NAI_G3xz_F2xz_b-2.0E0*1*I_NAI_G3xz_F2yz_b+1*I_NAI_G3xz_Pz;
  abcd[2463] = 4.0E0*I_NAI_G2x2y_H2x2yz_bb-2.0E0*1*I_NAI_G2x2y_F2xz_b-2.0E0*1*I_NAI_G2x2y_F2yz_b+1*I_NAI_G2x2y_Pz;
  abcd[2464] = 4.0E0*I_NAI_G2xyz_H2x2yz_bb-2.0E0*1*I_NAI_G2xyz_F2xz_b-2.0E0*1*I_NAI_G2xyz_F2yz_b+1*I_NAI_G2xyz_Pz;
  abcd[2465] = 4.0E0*I_NAI_G2x2z_H2x2yz_bb-2.0E0*1*I_NAI_G2x2z_F2xz_b-2.0E0*1*I_NAI_G2x2z_F2yz_b+1*I_NAI_G2x2z_Pz;
  abcd[2466] = 4.0E0*I_NAI_Gx3y_H2x2yz_bb-2.0E0*1*I_NAI_Gx3y_F2xz_b-2.0E0*1*I_NAI_Gx3y_F2yz_b+1*I_NAI_Gx3y_Pz;
  abcd[2467] = 4.0E0*I_NAI_Gx2yz_H2x2yz_bb-2.0E0*1*I_NAI_Gx2yz_F2xz_b-2.0E0*1*I_NAI_Gx2yz_F2yz_b+1*I_NAI_Gx2yz_Pz;
  abcd[2468] = 4.0E0*I_NAI_Gxy2z_H2x2yz_bb-2.0E0*1*I_NAI_Gxy2z_F2xz_b-2.0E0*1*I_NAI_Gxy2z_F2yz_b+1*I_NAI_Gxy2z_Pz;
  abcd[2469] = 4.0E0*I_NAI_Gx3z_H2x2yz_bb-2.0E0*1*I_NAI_Gx3z_F2xz_b-2.0E0*1*I_NAI_Gx3z_F2yz_b+1*I_NAI_Gx3z_Pz;
  abcd[2470] = 4.0E0*I_NAI_G4y_H2x2yz_bb-2.0E0*1*I_NAI_G4y_F2xz_b-2.0E0*1*I_NAI_G4y_F2yz_b+1*I_NAI_G4y_Pz;
  abcd[2471] = 4.0E0*I_NAI_G3yz_H2x2yz_bb-2.0E0*1*I_NAI_G3yz_F2xz_b-2.0E0*1*I_NAI_G3yz_F2yz_b+1*I_NAI_G3yz_Pz;
  abcd[2472] = 4.0E0*I_NAI_G2y2z_H2x2yz_bb-2.0E0*1*I_NAI_G2y2z_F2xz_b-2.0E0*1*I_NAI_G2y2z_F2yz_b+1*I_NAI_G2y2z_Pz;
  abcd[2473] = 4.0E0*I_NAI_Gy3z_H2x2yz_bb-2.0E0*1*I_NAI_Gy3z_F2xz_b-2.0E0*1*I_NAI_Gy3z_F2yz_b+1*I_NAI_Gy3z_Pz;
  abcd[2474] = 4.0E0*I_NAI_G4z_H2x2yz_bb-2.0E0*1*I_NAI_G4z_F2xz_b-2.0E0*1*I_NAI_G4z_F2yz_b+1*I_NAI_G4z_Pz;
  abcd[2475] = 4.0E0*I_NAI_G4x_H2xy2z_bb-2.0E0*1*I_NAI_G4x_Fy2z_b;
  abcd[2476] = 4.0E0*I_NAI_G3xy_H2xy2z_bb-2.0E0*1*I_NAI_G3xy_Fy2z_b;
  abcd[2477] = 4.0E0*I_NAI_G3xz_H2xy2z_bb-2.0E0*1*I_NAI_G3xz_Fy2z_b;
  abcd[2478] = 4.0E0*I_NAI_G2x2y_H2xy2z_bb-2.0E0*1*I_NAI_G2x2y_Fy2z_b;
  abcd[2479] = 4.0E0*I_NAI_G2xyz_H2xy2z_bb-2.0E0*1*I_NAI_G2xyz_Fy2z_b;
  abcd[2480] = 4.0E0*I_NAI_G2x2z_H2xy2z_bb-2.0E0*1*I_NAI_G2x2z_Fy2z_b;
  abcd[2481] = 4.0E0*I_NAI_Gx3y_H2xy2z_bb-2.0E0*1*I_NAI_Gx3y_Fy2z_b;
  abcd[2482] = 4.0E0*I_NAI_Gx2yz_H2xy2z_bb-2.0E0*1*I_NAI_Gx2yz_Fy2z_b;
  abcd[2483] = 4.0E0*I_NAI_Gxy2z_H2xy2z_bb-2.0E0*1*I_NAI_Gxy2z_Fy2z_b;
  abcd[2484] = 4.0E0*I_NAI_Gx3z_H2xy2z_bb-2.0E0*1*I_NAI_Gx3z_Fy2z_b;
  abcd[2485] = 4.0E0*I_NAI_G4y_H2xy2z_bb-2.0E0*1*I_NAI_G4y_Fy2z_b;
  abcd[2486] = 4.0E0*I_NAI_G3yz_H2xy2z_bb-2.0E0*1*I_NAI_G3yz_Fy2z_b;
  abcd[2487] = 4.0E0*I_NAI_G2y2z_H2xy2z_bb-2.0E0*1*I_NAI_G2y2z_Fy2z_b;
  abcd[2488] = 4.0E0*I_NAI_Gy3z_H2xy2z_bb-2.0E0*1*I_NAI_Gy3z_Fy2z_b;
  abcd[2489] = 4.0E0*I_NAI_G4z_H2xy2z_bb-2.0E0*1*I_NAI_G4z_Fy2z_b;
  abcd[2490] = 4.0E0*I_NAI_G4x_Hx4y_bb-2.0E0*3*I_NAI_G4x_Fx2y_b;
  abcd[2491] = 4.0E0*I_NAI_G3xy_Hx4y_bb-2.0E0*3*I_NAI_G3xy_Fx2y_b;
  abcd[2492] = 4.0E0*I_NAI_G3xz_Hx4y_bb-2.0E0*3*I_NAI_G3xz_Fx2y_b;
  abcd[2493] = 4.0E0*I_NAI_G2x2y_Hx4y_bb-2.0E0*3*I_NAI_G2x2y_Fx2y_b;
  abcd[2494] = 4.0E0*I_NAI_G2xyz_Hx4y_bb-2.0E0*3*I_NAI_G2xyz_Fx2y_b;
  abcd[2495] = 4.0E0*I_NAI_G2x2z_Hx4y_bb-2.0E0*3*I_NAI_G2x2z_Fx2y_b;
  abcd[2496] = 4.0E0*I_NAI_Gx3y_Hx4y_bb-2.0E0*3*I_NAI_Gx3y_Fx2y_b;
  abcd[2497] = 4.0E0*I_NAI_Gx2yz_Hx4y_bb-2.0E0*3*I_NAI_Gx2yz_Fx2y_b;
  abcd[2498] = 4.0E0*I_NAI_Gxy2z_Hx4y_bb-2.0E0*3*I_NAI_Gxy2z_Fx2y_b;
  abcd[2499] = 4.0E0*I_NAI_Gx3z_Hx4y_bb-2.0E0*3*I_NAI_Gx3z_Fx2y_b;
  abcd[2500] = 4.0E0*I_NAI_G4y_Hx4y_bb-2.0E0*3*I_NAI_G4y_Fx2y_b;
  abcd[2501] = 4.0E0*I_NAI_G3yz_Hx4y_bb-2.0E0*3*I_NAI_G3yz_Fx2y_b;
  abcd[2502] = 4.0E0*I_NAI_G2y2z_Hx4y_bb-2.0E0*3*I_NAI_G2y2z_Fx2y_b;
  abcd[2503] = 4.0E0*I_NAI_Gy3z_Hx4y_bb-2.0E0*3*I_NAI_Gy3z_Fx2y_b;
  abcd[2504] = 4.0E0*I_NAI_G4z_Hx4y_bb-2.0E0*3*I_NAI_G4z_Fx2y_b;
  abcd[2505] = 4.0E0*I_NAI_G4x_Hx3yz_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2506] = 4.0E0*I_NAI_G3xy_Hx3yz_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2507] = 4.0E0*I_NAI_G3xz_Hx3yz_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2508] = 4.0E0*I_NAI_G2x2y_Hx3yz_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2509] = 4.0E0*I_NAI_G2xyz_Hx3yz_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2510] = 4.0E0*I_NAI_G2x2z_Hx3yz_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2511] = 4.0E0*I_NAI_Gx3y_Hx3yz_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2512] = 4.0E0*I_NAI_Gx2yz_Hx3yz_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2513] = 4.0E0*I_NAI_Gxy2z_Hx3yz_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2514] = 4.0E0*I_NAI_Gx3z_Hx3yz_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2515] = 4.0E0*I_NAI_G4y_Hx3yz_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2516] = 4.0E0*I_NAI_G3yz_Hx3yz_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2517] = 4.0E0*I_NAI_G2y2z_Hx3yz_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2518] = 4.0E0*I_NAI_Gy3z_Hx3yz_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2519] = 4.0E0*I_NAI_G4z_Hx3yz_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2520] = 4.0E0*I_NAI_G4x_Hx2y2z_bb-2.0E0*1*I_NAI_G4x_Fx2z_b;
  abcd[2521] = 4.0E0*I_NAI_G3xy_Hx2y2z_bb-2.0E0*1*I_NAI_G3xy_Fx2z_b;
  abcd[2522] = 4.0E0*I_NAI_G3xz_Hx2y2z_bb-2.0E0*1*I_NAI_G3xz_Fx2z_b;
  abcd[2523] = 4.0E0*I_NAI_G2x2y_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2z_b;
  abcd[2524] = 4.0E0*I_NAI_G2xyz_Hx2y2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2z_b;
  abcd[2525] = 4.0E0*I_NAI_G2x2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2z_b;
  abcd[2526] = 4.0E0*I_NAI_Gx3y_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2z_b;
  abcd[2527] = 4.0E0*I_NAI_Gx2yz_Hx2y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2z_b;
  abcd[2528] = 4.0E0*I_NAI_Gxy2z_Hx2y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2z_b;
  abcd[2529] = 4.0E0*I_NAI_Gx3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2z_b;
  abcd[2530] = 4.0E0*I_NAI_G4y_Hx2y2z_bb-2.0E0*1*I_NAI_G4y_Fx2z_b;
  abcd[2531] = 4.0E0*I_NAI_G3yz_Hx2y2z_bb-2.0E0*1*I_NAI_G3yz_Fx2z_b;
  abcd[2532] = 4.0E0*I_NAI_G2y2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2z_b;
  abcd[2533] = 4.0E0*I_NAI_Gy3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2z_b;
  abcd[2534] = 4.0E0*I_NAI_G4z_Hx2y2z_bb-2.0E0*1*I_NAI_G4z_Fx2z_b;
  abcd[2535] = 4.0E0*I_NAI_G4x_Hxy3z_bb;
  abcd[2536] = 4.0E0*I_NAI_G3xy_Hxy3z_bb;
  abcd[2537] = 4.0E0*I_NAI_G3xz_Hxy3z_bb;
  abcd[2538] = 4.0E0*I_NAI_G2x2y_Hxy3z_bb;
  abcd[2539] = 4.0E0*I_NAI_G2xyz_Hxy3z_bb;
  abcd[2540] = 4.0E0*I_NAI_G2x2z_Hxy3z_bb;
  abcd[2541] = 4.0E0*I_NAI_Gx3y_Hxy3z_bb;
  abcd[2542] = 4.0E0*I_NAI_Gx2yz_Hxy3z_bb;
  abcd[2543] = 4.0E0*I_NAI_Gxy2z_Hxy3z_bb;
  abcd[2544] = 4.0E0*I_NAI_Gx3z_Hxy3z_bb;
  abcd[2545] = 4.0E0*I_NAI_G4y_Hxy3z_bb;
  abcd[2546] = 4.0E0*I_NAI_G3yz_Hxy3z_bb;
  abcd[2547] = 4.0E0*I_NAI_G2y2z_Hxy3z_bb;
  abcd[2548] = 4.0E0*I_NAI_Gy3z_Hxy3z_bb;
  abcd[2549] = 4.0E0*I_NAI_G4z_Hxy3z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[2550] = 4.0E0*I_NAI_G4x_H4xz_bb-2.0E0*3*I_NAI_G4x_F2xz_b;
  abcd[2551] = 4.0E0*I_NAI_G3xy_H4xz_bb-2.0E0*3*I_NAI_G3xy_F2xz_b;
  abcd[2552] = 4.0E0*I_NAI_G3xz_H4xz_bb-2.0E0*3*I_NAI_G3xz_F2xz_b;
  abcd[2553] = 4.0E0*I_NAI_G2x2y_H4xz_bb-2.0E0*3*I_NAI_G2x2y_F2xz_b;
  abcd[2554] = 4.0E0*I_NAI_G2xyz_H4xz_bb-2.0E0*3*I_NAI_G2xyz_F2xz_b;
  abcd[2555] = 4.0E0*I_NAI_G2x2z_H4xz_bb-2.0E0*3*I_NAI_G2x2z_F2xz_b;
  abcd[2556] = 4.0E0*I_NAI_Gx3y_H4xz_bb-2.0E0*3*I_NAI_Gx3y_F2xz_b;
  abcd[2557] = 4.0E0*I_NAI_Gx2yz_H4xz_bb-2.0E0*3*I_NAI_Gx2yz_F2xz_b;
  abcd[2558] = 4.0E0*I_NAI_Gxy2z_H4xz_bb-2.0E0*3*I_NAI_Gxy2z_F2xz_b;
  abcd[2559] = 4.0E0*I_NAI_Gx3z_H4xz_bb-2.0E0*3*I_NAI_Gx3z_F2xz_b;
  abcd[2560] = 4.0E0*I_NAI_G4y_H4xz_bb-2.0E0*3*I_NAI_G4y_F2xz_b;
  abcd[2561] = 4.0E0*I_NAI_G3yz_H4xz_bb-2.0E0*3*I_NAI_G3yz_F2xz_b;
  abcd[2562] = 4.0E0*I_NAI_G2y2z_H4xz_bb-2.0E0*3*I_NAI_G2y2z_F2xz_b;
  abcd[2563] = 4.0E0*I_NAI_Gy3z_H4xz_bb-2.0E0*3*I_NAI_Gy3z_F2xz_b;
  abcd[2564] = 4.0E0*I_NAI_G4z_H4xz_bb-2.0E0*3*I_NAI_G4z_F2xz_b;
  abcd[2565] = 4.0E0*I_NAI_G4x_H3xyz_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2566] = 4.0E0*I_NAI_G3xy_H3xyz_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2567] = 4.0E0*I_NAI_G3xz_H3xyz_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2568] = 4.0E0*I_NAI_G2x2y_H3xyz_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2569] = 4.0E0*I_NAI_G2xyz_H3xyz_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2570] = 4.0E0*I_NAI_G2x2z_H3xyz_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2571] = 4.0E0*I_NAI_Gx3y_H3xyz_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2572] = 4.0E0*I_NAI_Gx2yz_H3xyz_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2573] = 4.0E0*I_NAI_Gxy2z_H3xyz_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2574] = 4.0E0*I_NAI_Gx3z_H3xyz_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2575] = 4.0E0*I_NAI_G4y_H3xyz_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2576] = 4.0E0*I_NAI_G3yz_H3xyz_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2577] = 4.0E0*I_NAI_G2y2z_H3xyz_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2578] = 4.0E0*I_NAI_Gy3z_H3xyz_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2579] = 4.0E0*I_NAI_G4z_H3xyz_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2580] = 4.0E0*I_NAI_G4x_H3x2z_bb-2.0E0*1*I_NAI_G4x_F3x_b-2.0E0*2*I_NAI_G4x_Fx2z_b+2*1*I_NAI_G4x_Px;
  abcd[2581] = 4.0E0*I_NAI_G3xy_H3x2z_bb-2.0E0*1*I_NAI_G3xy_F3x_b-2.0E0*2*I_NAI_G3xy_Fx2z_b+2*1*I_NAI_G3xy_Px;
  abcd[2582] = 4.0E0*I_NAI_G3xz_H3x2z_bb-2.0E0*1*I_NAI_G3xz_F3x_b-2.0E0*2*I_NAI_G3xz_Fx2z_b+2*1*I_NAI_G3xz_Px;
  abcd[2583] = 4.0E0*I_NAI_G2x2y_H3x2z_bb-2.0E0*1*I_NAI_G2x2y_F3x_b-2.0E0*2*I_NAI_G2x2y_Fx2z_b+2*1*I_NAI_G2x2y_Px;
  abcd[2584] = 4.0E0*I_NAI_G2xyz_H3x2z_bb-2.0E0*1*I_NAI_G2xyz_F3x_b-2.0E0*2*I_NAI_G2xyz_Fx2z_b+2*1*I_NAI_G2xyz_Px;
  abcd[2585] = 4.0E0*I_NAI_G2x2z_H3x2z_bb-2.0E0*1*I_NAI_G2x2z_F3x_b-2.0E0*2*I_NAI_G2x2z_Fx2z_b+2*1*I_NAI_G2x2z_Px;
  abcd[2586] = 4.0E0*I_NAI_Gx3y_H3x2z_bb-2.0E0*1*I_NAI_Gx3y_F3x_b-2.0E0*2*I_NAI_Gx3y_Fx2z_b+2*1*I_NAI_Gx3y_Px;
  abcd[2587] = 4.0E0*I_NAI_Gx2yz_H3x2z_bb-2.0E0*1*I_NAI_Gx2yz_F3x_b-2.0E0*2*I_NAI_Gx2yz_Fx2z_b+2*1*I_NAI_Gx2yz_Px;
  abcd[2588] = 4.0E0*I_NAI_Gxy2z_H3x2z_bb-2.0E0*1*I_NAI_Gxy2z_F3x_b-2.0E0*2*I_NAI_Gxy2z_Fx2z_b+2*1*I_NAI_Gxy2z_Px;
  abcd[2589] = 4.0E0*I_NAI_Gx3z_H3x2z_bb-2.0E0*1*I_NAI_Gx3z_F3x_b-2.0E0*2*I_NAI_Gx3z_Fx2z_b+2*1*I_NAI_Gx3z_Px;
  abcd[2590] = 4.0E0*I_NAI_G4y_H3x2z_bb-2.0E0*1*I_NAI_G4y_F3x_b-2.0E0*2*I_NAI_G4y_Fx2z_b+2*1*I_NAI_G4y_Px;
  abcd[2591] = 4.0E0*I_NAI_G3yz_H3x2z_bb-2.0E0*1*I_NAI_G3yz_F3x_b-2.0E0*2*I_NAI_G3yz_Fx2z_b+2*1*I_NAI_G3yz_Px;
  abcd[2592] = 4.0E0*I_NAI_G2y2z_H3x2z_bb-2.0E0*1*I_NAI_G2y2z_F3x_b-2.0E0*2*I_NAI_G2y2z_Fx2z_b+2*1*I_NAI_G2y2z_Px;
  abcd[2593] = 4.0E0*I_NAI_Gy3z_H3x2z_bb-2.0E0*1*I_NAI_Gy3z_F3x_b-2.0E0*2*I_NAI_Gy3z_Fx2z_b+2*1*I_NAI_Gy3z_Px;
  abcd[2594] = 4.0E0*I_NAI_G4z_H3x2z_bb-2.0E0*1*I_NAI_G4z_F3x_b-2.0E0*2*I_NAI_G4z_Fx2z_b+2*1*I_NAI_G4z_Px;
  abcd[2595] = 4.0E0*I_NAI_G4x_H2x2yz_bb-2.0E0*1*I_NAI_G4x_F2yz_b;
  abcd[2596] = 4.0E0*I_NAI_G3xy_H2x2yz_bb-2.0E0*1*I_NAI_G3xy_F2yz_b;
  abcd[2597] = 4.0E0*I_NAI_G3xz_H2x2yz_bb-2.0E0*1*I_NAI_G3xz_F2yz_b;
  abcd[2598] = 4.0E0*I_NAI_G2x2y_H2x2yz_bb-2.0E0*1*I_NAI_G2x2y_F2yz_b;
  abcd[2599] = 4.0E0*I_NAI_G2xyz_H2x2yz_bb-2.0E0*1*I_NAI_G2xyz_F2yz_b;
  abcd[2600] = 4.0E0*I_NAI_G2x2z_H2x2yz_bb-2.0E0*1*I_NAI_G2x2z_F2yz_b;
  abcd[2601] = 4.0E0*I_NAI_Gx3y_H2x2yz_bb-2.0E0*1*I_NAI_Gx3y_F2yz_b;
  abcd[2602] = 4.0E0*I_NAI_Gx2yz_H2x2yz_bb-2.0E0*1*I_NAI_Gx2yz_F2yz_b;
  abcd[2603] = 4.0E0*I_NAI_Gxy2z_H2x2yz_bb-2.0E0*1*I_NAI_Gxy2z_F2yz_b;
  abcd[2604] = 4.0E0*I_NAI_Gx3z_H2x2yz_bb-2.0E0*1*I_NAI_Gx3z_F2yz_b;
  abcd[2605] = 4.0E0*I_NAI_G4y_H2x2yz_bb-2.0E0*1*I_NAI_G4y_F2yz_b;
  abcd[2606] = 4.0E0*I_NAI_G3yz_H2x2yz_bb-2.0E0*1*I_NAI_G3yz_F2yz_b;
  abcd[2607] = 4.0E0*I_NAI_G2y2z_H2x2yz_bb-2.0E0*1*I_NAI_G2y2z_F2yz_b;
  abcd[2608] = 4.0E0*I_NAI_Gy3z_H2x2yz_bb-2.0E0*1*I_NAI_Gy3z_F2yz_b;
  abcd[2609] = 4.0E0*I_NAI_G4z_H2x2yz_bb-2.0E0*1*I_NAI_G4z_F2yz_b;
  abcd[2610] = 4.0E0*I_NAI_G4x_H2xy2z_bb-2.0E0*1*I_NAI_G4x_F2xy_b-2.0E0*1*I_NAI_G4x_Fy2z_b+1*I_NAI_G4x_Py;
  abcd[2611] = 4.0E0*I_NAI_G3xy_H2xy2z_bb-2.0E0*1*I_NAI_G3xy_F2xy_b-2.0E0*1*I_NAI_G3xy_Fy2z_b+1*I_NAI_G3xy_Py;
  abcd[2612] = 4.0E0*I_NAI_G3xz_H2xy2z_bb-2.0E0*1*I_NAI_G3xz_F2xy_b-2.0E0*1*I_NAI_G3xz_Fy2z_b+1*I_NAI_G3xz_Py;
  abcd[2613] = 4.0E0*I_NAI_G2x2y_H2xy2z_bb-2.0E0*1*I_NAI_G2x2y_F2xy_b-2.0E0*1*I_NAI_G2x2y_Fy2z_b+1*I_NAI_G2x2y_Py;
  abcd[2614] = 4.0E0*I_NAI_G2xyz_H2xy2z_bb-2.0E0*1*I_NAI_G2xyz_F2xy_b-2.0E0*1*I_NAI_G2xyz_Fy2z_b+1*I_NAI_G2xyz_Py;
  abcd[2615] = 4.0E0*I_NAI_G2x2z_H2xy2z_bb-2.0E0*1*I_NAI_G2x2z_F2xy_b-2.0E0*1*I_NAI_G2x2z_Fy2z_b+1*I_NAI_G2x2z_Py;
  abcd[2616] = 4.0E0*I_NAI_Gx3y_H2xy2z_bb-2.0E0*1*I_NAI_Gx3y_F2xy_b-2.0E0*1*I_NAI_Gx3y_Fy2z_b+1*I_NAI_Gx3y_Py;
  abcd[2617] = 4.0E0*I_NAI_Gx2yz_H2xy2z_bb-2.0E0*1*I_NAI_Gx2yz_F2xy_b-2.0E0*1*I_NAI_Gx2yz_Fy2z_b+1*I_NAI_Gx2yz_Py;
  abcd[2618] = 4.0E0*I_NAI_Gxy2z_H2xy2z_bb-2.0E0*1*I_NAI_Gxy2z_F2xy_b-2.0E0*1*I_NAI_Gxy2z_Fy2z_b+1*I_NAI_Gxy2z_Py;
  abcd[2619] = 4.0E0*I_NAI_Gx3z_H2xy2z_bb-2.0E0*1*I_NAI_Gx3z_F2xy_b-2.0E0*1*I_NAI_Gx3z_Fy2z_b+1*I_NAI_Gx3z_Py;
  abcd[2620] = 4.0E0*I_NAI_G4y_H2xy2z_bb-2.0E0*1*I_NAI_G4y_F2xy_b-2.0E0*1*I_NAI_G4y_Fy2z_b+1*I_NAI_G4y_Py;
  abcd[2621] = 4.0E0*I_NAI_G3yz_H2xy2z_bb-2.0E0*1*I_NAI_G3yz_F2xy_b-2.0E0*1*I_NAI_G3yz_Fy2z_b+1*I_NAI_G3yz_Py;
  abcd[2622] = 4.0E0*I_NAI_G2y2z_H2xy2z_bb-2.0E0*1*I_NAI_G2y2z_F2xy_b-2.0E0*1*I_NAI_G2y2z_Fy2z_b+1*I_NAI_G2y2z_Py;
  abcd[2623] = 4.0E0*I_NAI_Gy3z_H2xy2z_bb-2.0E0*1*I_NAI_Gy3z_F2xy_b-2.0E0*1*I_NAI_Gy3z_Fy2z_b+1*I_NAI_Gy3z_Py;
  abcd[2624] = 4.0E0*I_NAI_G4z_H2xy2z_bb-2.0E0*1*I_NAI_G4z_F2xy_b-2.0E0*1*I_NAI_G4z_Fy2z_b+1*I_NAI_G4z_Py;
  abcd[2625] = 4.0E0*I_NAI_G4x_H2x3z_bb-2.0E0*2*I_NAI_G4x_F2xz_b-2.0E0*1*I_NAI_G4x_F3z_b+2*I_NAI_G4x_Pz;
  abcd[2626] = 4.0E0*I_NAI_G3xy_H2x3z_bb-2.0E0*2*I_NAI_G3xy_F2xz_b-2.0E0*1*I_NAI_G3xy_F3z_b+2*I_NAI_G3xy_Pz;
  abcd[2627] = 4.0E0*I_NAI_G3xz_H2x3z_bb-2.0E0*2*I_NAI_G3xz_F2xz_b-2.0E0*1*I_NAI_G3xz_F3z_b+2*I_NAI_G3xz_Pz;
  abcd[2628] = 4.0E0*I_NAI_G2x2y_H2x3z_bb-2.0E0*2*I_NAI_G2x2y_F2xz_b-2.0E0*1*I_NAI_G2x2y_F3z_b+2*I_NAI_G2x2y_Pz;
  abcd[2629] = 4.0E0*I_NAI_G2xyz_H2x3z_bb-2.0E0*2*I_NAI_G2xyz_F2xz_b-2.0E0*1*I_NAI_G2xyz_F3z_b+2*I_NAI_G2xyz_Pz;
  abcd[2630] = 4.0E0*I_NAI_G2x2z_H2x3z_bb-2.0E0*2*I_NAI_G2x2z_F2xz_b-2.0E0*1*I_NAI_G2x2z_F3z_b+2*I_NAI_G2x2z_Pz;
  abcd[2631] = 4.0E0*I_NAI_Gx3y_H2x3z_bb-2.0E0*2*I_NAI_Gx3y_F2xz_b-2.0E0*1*I_NAI_Gx3y_F3z_b+2*I_NAI_Gx3y_Pz;
  abcd[2632] = 4.0E0*I_NAI_Gx2yz_H2x3z_bb-2.0E0*2*I_NAI_Gx2yz_F2xz_b-2.0E0*1*I_NAI_Gx2yz_F3z_b+2*I_NAI_Gx2yz_Pz;
  abcd[2633] = 4.0E0*I_NAI_Gxy2z_H2x3z_bb-2.0E0*2*I_NAI_Gxy2z_F2xz_b-2.0E0*1*I_NAI_Gxy2z_F3z_b+2*I_NAI_Gxy2z_Pz;
  abcd[2634] = 4.0E0*I_NAI_Gx3z_H2x3z_bb-2.0E0*2*I_NAI_Gx3z_F2xz_b-2.0E0*1*I_NAI_Gx3z_F3z_b+2*I_NAI_Gx3z_Pz;
  abcd[2635] = 4.0E0*I_NAI_G4y_H2x3z_bb-2.0E0*2*I_NAI_G4y_F2xz_b-2.0E0*1*I_NAI_G4y_F3z_b+2*I_NAI_G4y_Pz;
  abcd[2636] = 4.0E0*I_NAI_G3yz_H2x3z_bb-2.0E0*2*I_NAI_G3yz_F2xz_b-2.0E0*1*I_NAI_G3yz_F3z_b+2*I_NAI_G3yz_Pz;
  abcd[2637] = 4.0E0*I_NAI_G2y2z_H2x3z_bb-2.0E0*2*I_NAI_G2y2z_F2xz_b-2.0E0*1*I_NAI_G2y2z_F3z_b+2*I_NAI_G2y2z_Pz;
  abcd[2638] = 4.0E0*I_NAI_Gy3z_H2x3z_bb-2.0E0*2*I_NAI_Gy3z_F2xz_b-2.0E0*1*I_NAI_Gy3z_F3z_b+2*I_NAI_Gy3z_Pz;
  abcd[2639] = 4.0E0*I_NAI_G4z_H2x3z_bb-2.0E0*2*I_NAI_G4z_F2xz_b-2.0E0*1*I_NAI_G4z_F3z_b+2*I_NAI_G4z_Pz;
  abcd[2640] = 4.0E0*I_NAI_G4x_Hx3yz_bb;
  abcd[2641] = 4.0E0*I_NAI_G3xy_Hx3yz_bb;
  abcd[2642] = 4.0E0*I_NAI_G3xz_Hx3yz_bb;
  abcd[2643] = 4.0E0*I_NAI_G2x2y_Hx3yz_bb;
  abcd[2644] = 4.0E0*I_NAI_G2xyz_Hx3yz_bb;
  abcd[2645] = 4.0E0*I_NAI_G2x2z_Hx3yz_bb;
  abcd[2646] = 4.0E0*I_NAI_Gx3y_Hx3yz_bb;
  abcd[2647] = 4.0E0*I_NAI_Gx2yz_Hx3yz_bb;
  abcd[2648] = 4.0E0*I_NAI_Gxy2z_Hx3yz_bb;
  abcd[2649] = 4.0E0*I_NAI_Gx3z_Hx3yz_bb;
  abcd[2650] = 4.0E0*I_NAI_G4y_Hx3yz_bb;
  abcd[2651] = 4.0E0*I_NAI_G3yz_Hx3yz_bb;
  abcd[2652] = 4.0E0*I_NAI_G2y2z_Hx3yz_bb;
  abcd[2653] = 4.0E0*I_NAI_Gy3z_Hx3yz_bb;
  abcd[2654] = 4.0E0*I_NAI_G4z_Hx3yz_bb;
  abcd[2655] = 4.0E0*I_NAI_G4x_Hx2y2z_bb-2.0E0*1*I_NAI_G4x_Fx2y_b;
  abcd[2656] = 4.0E0*I_NAI_G3xy_Hx2y2z_bb-2.0E0*1*I_NAI_G3xy_Fx2y_b;
  abcd[2657] = 4.0E0*I_NAI_G3xz_Hx2y2z_bb-2.0E0*1*I_NAI_G3xz_Fx2y_b;
  abcd[2658] = 4.0E0*I_NAI_G2x2y_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2y_b;
  abcd[2659] = 4.0E0*I_NAI_G2xyz_Hx2y2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2y_b;
  abcd[2660] = 4.0E0*I_NAI_G2x2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2y_b;
  abcd[2661] = 4.0E0*I_NAI_Gx3y_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2y_b;
  abcd[2662] = 4.0E0*I_NAI_Gx2yz_Hx2y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2y_b;
  abcd[2663] = 4.0E0*I_NAI_Gxy2z_Hx2y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2y_b;
  abcd[2664] = 4.0E0*I_NAI_Gx3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2y_b;
  abcd[2665] = 4.0E0*I_NAI_G4y_Hx2y2z_bb-2.0E0*1*I_NAI_G4y_Fx2y_b;
  abcd[2666] = 4.0E0*I_NAI_G3yz_Hx2y2z_bb-2.0E0*1*I_NAI_G3yz_Fx2y_b;
  abcd[2667] = 4.0E0*I_NAI_G2y2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2y_b;
  abcd[2668] = 4.0E0*I_NAI_Gy3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2y_b;
  abcd[2669] = 4.0E0*I_NAI_G4z_Hx2y2z_bb-2.0E0*1*I_NAI_G4z_Fx2y_b;
  abcd[2670] = 4.0E0*I_NAI_G4x_Hxy3z_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2671] = 4.0E0*I_NAI_G3xy_Hxy3z_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2672] = 4.0E0*I_NAI_G3xz_Hxy3z_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2673] = 4.0E0*I_NAI_G2x2y_Hxy3z_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2674] = 4.0E0*I_NAI_G2xyz_Hxy3z_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2675] = 4.0E0*I_NAI_G2x2z_Hxy3z_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2676] = 4.0E0*I_NAI_Gx3y_Hxy3z_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2677] = 4.0E0*I_NAI_Gx2yz_Hxy3z_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2678] = 4.0E0*I_NAI_Gxy2z_Hxy3z_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2679] = 4.0E0*I_NAI_Gx3z_Hxy3z_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2680] = 4.0E0*I_NAI_G4y_Hxy3z_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2681] = 4.0E0*I_NAI_G3yz_Hxy3z_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2682] = 4.0E0*I_NAI_G2y2z_Hxy3z_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2683] = 4.0E0*I_NAI_Gy3z_Hxy3z_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2684] = 4.0E0*I_NAI_G4z_Hxy3z_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2685] = 4.0E0*I_NAI_G4x_Hx4z_bb-2.0E0*3*I_NAI_G4x_Fx2z_b;
  abcd[2686] = 4.0E0*I_NAI_G3xy_Hx4z_bb-2.0E0*3*I_NAI_G3xy_Fx2z_b;
  abcd[2687] = 4.0E0*I_NAI_G3xz_Hx4z_bb-2.0E0*3*I_NAI_G3xz_Fx2z_b;
  abcd[2688] = 4.0E0*I_NAI_G2x2y_Hx4z_bb-2.0E0*3*I_NAI_G2x2y_Fx2z_b;
  abcd[2689] = 4.0E0*I_NAI_G2xyz_Hx4z_bb-2.0E0*3*I_NAI_G2xyz_Fx2z_b;
  abcd[2690] = 4.0E0*I_NAI_G2x2z_Hx4z_bb-2.0E0*3*I_NAI_G2x2z_Fx2z_b;
  abcd[2691] = 4.0E0*I_NAI_Gx3y_Hx4z_bb-2.0E0*3*I_NAI_Gx3y_Fx2z_b;
  abcd[2692] = 4.0E0*I_NAI_Gx2yz_Hx4z_bb-2.0E0*3*I_NAI_Gx2yz_Fx2z_b;
  abcd[2693] = 4.0E0*I_NAI_Gxy2z_Hx4z_bb-2.0E0*3*I_NAI_Gxy2z_Fx2z_b;
  abcd[2694] = 4.0E0*I_NAI_Gx3z_Hx4z_bb-2.0E0*3*I_NAI_Gx3z_Fx2z_b;
  abcd[2695] = 4.0E0*I_NAI_G4y_Hx4z_bb-2.0E0*3*I_NAI_G4y_Fx2z_b;
  abcd[2696] = 4.0E0*I_NAI_G3yz_Hx4z_bb-2.0E0*3*I_NAI_G3yz_Fx2z_b;
  abcd[2697] = 4.0E0*I_NAI_G2y2z_Hx4z_bb-2.0E0*3*I_NAI_G2y2z_Fx2z_b;
  abcd[2698] = 4.0E0*I_NAI_Gy3z_Hx4z_bb-2.0E0*3*I_NAI_Gy3z_Fx2z_b;
  abcd[2699] = 4.0E0*I_NAI_G4z_Hx4z_bb-2.0E0*3*I_NAI_G4z_Fx2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[2700] = 4.0E0*I_NAI_G4x_H3x2y_bb-2.0E0*1*I_NAI_G4x_F3x_b;
  abcd[2701] = 4.0E0*I_NAI_G3xy_H3x2y_bb-2.0E0*1*I_NAI_G3xy_F3x_b;
  abcd[2702] = 4.0E0*I_NAI_G3xz_H3x2y_bb-2.0E0*1*I_NAI_G3xz_F3x_b;
  abcd[2703] = 4.0E0*I_NAI_G2x2y_H3x2y_bb-2.0E0*1*I_NAI_G2x2y_F3x_b;
  abcd[2704] = 4.0E0*I_NAI_G2xyz_H3x2y_bb-2.0E0*1*I_NAI_G2xyz_F3x_b;
  abcd[2705] = 4.0E0*I_NAI_G2x2z_H3x2y_bb-2.0E0*1*I_NAI_G2x2z_F3x_b;
  abcd[2706] = 4.0E0*I_NAI_Gx3y_H3x2y_bb-2.0E0*1*I_NAI_Gx3y_F3x_b;
  abcd[2707] = 4.0E0*I_NAI_Gx2yz_H3x2y_bb-2.0E0*1*I_NAI_Gx2yz_F3x_b;
  abcd[2708] = 4.0E0*I_NAI_Gxy2z_H3x2y_bb-2.0E0*1*I_NAI_Gxy2z_F3x_b;
  abcd[2709] = 4.0E0*I_NAI_Gx3z_H3x2y_bb-2.0E0*1*I_NAI_Gx3z_F3x_b;
  abcd[2710] = 4.0E0*I_NAI_G4y_H3x2y_bb-2.0E0*1*I_NAI_G4y_F3x_b;
  abcd[2711] = 4.0E0*I_NAI_G3yz_H3x2y_bb-2.0E0*1*I_NAI_G3yz_F3x_b;
  abcd[2712] = 4.0E0*I_NAI_G2y2z_H3x2y_bb-2.0E0*1*I_NAI_G2y2z_F3x_b;
  abcd[2713] = 4.0E0*I_NAI_Gy3z_H3x2y_bb-2.0E0*1*I_NAI_Gy3z_F3x_b;
  abcd[2714] = 4.0E0*I_NAI_G4z_H3x2y_bb-2.0E0*1*I_NAI_G4z_F3x_b;
  abcd[2715] = 4.0E0*I_NAI_G4x_H2x3y_bb-2.0E0*1*I_NAI_G4x_F2xy_b-2.0E0*2*I_NAI_G4x_F2xy_b;
  abcd[2716] = 4.0E0*I_NAI_G3xy_H2x3y_bb-2.0E0*1*I_NAI_G3xy_F2xy_b-2.0E0*2*I_NAI_G3xy_F2xy_b;
  abcd[2717] = 4.0E0*I_NAI_G3xz_H2x3y_bb-2.0E0*1*I_NAI_G3xz_F2xy_b-2.0E0*2*I_NAI_G3xz_F2xy_b;
  abcd[2718] = 4.0E0*I_NAI_G2x2y_H2x3y_bb-2.0E0*1*I_NAI_G2x2y_F2xy_b-2.0E0*2*I_NAI_G2x2y_F2xy_b;
  abcd[2719] = 4.0E0*I_NAI_G2xyz_H2x3y_bb-2.0E0*1*I_NAI_G2xyz_F2xy_b-2.0E0*2*I_NAI_G2xyz_F2xy_b;
  abcd[2720] = 4.0E0*I_NAI_G2x2z_H2x3y_bb-2.0E0*1*I_NAI_G2x2z_F2xy_b-2.0E0*2*I_NAI_G2x2z_F2xy_b;
  abcd[2721] = 4.0E0*I_NAI_Gx3y_H2x3y_bb-2.0E0*1*I_NAI_Gx3y_F2xy_b-2.0E0*2*I_NAI_Gx3y_F2xy_b;
  abcd[2722] = 4.0E0*I_NAI_Gx2yz_H2x3y_bb-2.0E0*1*I_NAI_Gx2yz_F2xy_b-2.0E0*2*I_NAI_Gx2yz_F2xy_b;
  abcd[2723] = 4.0E0*I_NAI_Gxy2z_H2x3y_bb-2.0E0*1*I_NAI_Gxy2z_F2xy_b-2.0E0*2*I_NAI_Gxy2z_F2xy_b;
  abcd[2724] = 4.0E0*I_NAI_Gx3z_H2x3y_bb-2.0E0*1*I_NAI_Gx3z_F2xy_b-2.0E0*2*I_NAI_Gx3z_F2xy_b;
  abcd[2725] = 4.0E0*I_NAI_G4y_H2x3y_bb-2.0E0*1*I_NAI_G4y_F2xy_b-2.0E0*2*I_NAI_G4y_F2xy_b;
  abcd[2726] = 4.0E0*I_NAI_G3yz_H2x3y_bb-2.0E0*1*I_NAI_G3yz_F2xy_b-2.0E0*2*I_NAI_G3yz_F2xy_b;
  abcd[2727] = 4.0E0*I_NAI_G2y2z_H2x3y_bb-2.0E0*1*I_NAI_G2y2z_F2xy_b-2.0E0*2*I_NAI_G2y2z_F2xy_b;
  abcd[2728] = 4.0E0*I_NAI_Gy3z_H2x3y_bb-2.0E0*1*I_NAI_Gy3z_F2xy_b-2.0E0*2*I_NAI_Gy3z_F2xy_b;
  abcd[2729] = 4.0E0*I_NAI_G4z_H2x3y_bb-2.0E0*1*I_NAI_G4z_F2xy_b-2.0E0*2*I_NAI_G4z_F2xy_b;
  abcd[2730] = 4.0E0*I_NAI_G4x_H2x2yz_bb-2.0E0*1*I_NAI_G4x_F2xz_b;
  abcd[2731] = 4.0E0*I_NAI_G3xy_H2x2yz_bb-2.0E0*1*I_NAI_G3xy_F2xz_b;
  abcd[2732] = 4.0E0*I_NAI_G3xz_H2x2yz_bb-2.0E0*1*I_NAI_G3xz_F2xz_b;
  abcd[2733] = 4.0E0*I_NAI_G2x2y_H2x2yz_bb-2.0E0*1*I_NAI_G2x2y_F2xz_b;
  abcd[2734] = 4.0E0*I_NAI_G2xyz_H2x2yz_bb-2.0E0*1*I_NAI_G2xyz_F2xz_b;
  abcd[2735] = 4.0E0*I_NAI_G2x2z_H2x2yz_bb-2.0E0*1*I_NAI_G2x2z_F2xz_b;
  abcd[2736] = 4.0E0*I_NAI_Gx3y_H2x2yz_bb-2.0E0*1*I_NAI_Gx3y_F2xz_b;
  abcd[2737] = 4.0E0*I_NAI_Gx2yz_H2x2yz_bb-2.0E0*1*I_NAI_Gx2yz_F2xz_b;
  abcd[2738] = 4.0E0*I_NAI_Gxy2z_H2x2yz_bb-2.0E0*1*I_NAI_Gxy2z_F2xz_b;
  abcd[2739] = 4.0E0*I_NAI_Gx3z_H2x2yz_bb-2.0E0*1*I_NAI_Gx3z_F2xz_b;
  abcd[2740] = 4.0E0*I_NAI_G4y_H2x2yz_bb-2.0E0*1*I_NAI_G4y_F2xz_b;
  abcd[2741] = 4.0E0*I_NAI_G3yz_H2x2yz_bb-2.0E0*1*I_NAI_G3yz_F2xz_b;
  abcd[2742] = 4.0E0*I_NAI_G2y2z_H2x2yz_bb-2.0E0*1*I_NAI_G2y2z_F2xz_b;
  abcd[2743] = 4.0E0*I_NAI_Gy3z_H2x2yz_bb-2.0E0*1*I_NAI_Gy3z_F2xz_b;
  abcd[2744] = 4.0E0*I_NAI_G4z_H2x2yz_bb-2.0E0*1*I_NAI_G4z_F2xz_b;
  abcd[2745] = 4.0E0*I_NAI_G4x_Hx4y_bb-2.0E0*2*I_NAI_G4x_Fx2y_b-2.0E0*3*I_NAI_G4x_Fx2y_b+2*1*I_NAI_G4x_Px;
  abcd[2746] = 4.0E0*I_NAI_G3xy_Hx4y_bb-2.0E0*2*I_NAI_G3xy_Fx2y_b-2.0E0*3*I_NAI_G3xy_Fx2y_b+2*1*I_NAI_G3xy_Px;
  abcd[2747] = 4.0E0*I_NAI_G3xz_Hx4y_bb-2.0E0*2*I_NAI_G3xz_Fx2y_b-2.0E0*3*I_NAI_G3xz_Fx2y_b+2*1*I_NAI_G3xz_Px;
  abcd[2748] = 4.0E0*I_NAI_G2x2y_Hx4y_bb-2.0E0*2*I_NAI_G2x2y_Fx2y_b-2.0E0*3*I_NAI_G2x2y_Fx2y_b+2*1*I_NAI_G2x2y_Px;
  abcd[2749] = 4.0E0*I_NAI_G2xyz_Hx4y_bb-2.0E0*2*I_NAI_G2xyz_Fx2y_b-2.0E0*3*I_NAI_G2xyz_Fx2y_b+2*1*I_NAI_G2xyz_Px;
  abcd[2750] = 4.0E0*I_NAI_G2x2z_Hx4y_bb-2.0E0*2*I_NAI_G2x2z_Fx2y_b-2.0E0*3*I_NAI_G2x2z_Fx2y_b+2*1*I_NAI_G2x2z_Px;
  abcd[2751] = 4.0E0*I_NAI_Gx3y_Hx4y_bb-2.0E0*2*I_NAI_Gx3y_Fx2y_b-2.0E0*3*I_NAI_Gx3y_Fx2y_b+2*1*I_NAI_Gx3y_Px;
  abcd[2752] = 4.0E0*I_NAI_Gx2yz_Hx4y_bb-2.0E0*2*I_NAI_Gx2yz_Fx2y_b-2.0E0*3*I_NAI_Gx2yz_Fx2y_b+2*1*I_NAI_Gx2yz_Px;
  abcd[2753] = 4.0E0*I_NAI_Gxy2z_Hx4y_bb-2.0E0*2*I_NAI_Gxy2z_Fx2y_b-2.0E0*3*I_NAI_Gxy2z_Fx2y_b+2*1*I_NAI_Gxy2z_Px;
  abcd[2754] = 4.0E0*I_NAI_Gx3z_Hx4y_bb-2.0E0*2*I_NAI_Gx3z_Fx2y_b-2.0E0*3*I_NAI_Gx3z_Fx2y_b+2*1*I_NAI_Gx3z_Px;
  abcd[2755] = 4.0E0*I_NAI_G4y_Hx4y_bb-2.0E0*2*I_NAI_G4y_Fx2y_b-2.0E0*3*I_NAI_G4y_Fx2y_b+2*1*I_NAI_G4y_Px;
  abcd[2756] = 4.0E0*I_NAI_G3yz_Hx4y_bb-2.0E0*2*I_NAI_G3yz_Fx2y_b-2.0E0*3*I_NAI_G3yz_Fx2y_b+2*1*I_NAI_G3yz_Px;
  abcd[2757] = 4.0E0*I_NAI_G2y2z_Hx4y_bb-2.0E0*2*I_NAI_G2y2z_Fx2y_b-2.0E0*3*I_NAI_G2y2z_Fx2y_b+2*1*I_NAI_G2y2z_Px;
  abcd[2758] = 4.0E0*I_NAI_Gy3z_Hx4y_bb-2.0E0*2*I_NAI_Gy3z_Fx2y_b-2.0E0*3*I_NAI_Gy3z_Fx2y_b+2*1*I_NAI_Gy3z_Px;
  abcd[2759] = 4.0E0*I_NAI_G4z_Hx4y_bb-2.0E0*2*I_NAI_G4z_Fx2y_b-2.0E0*3*I_NAI_G4z_Fx2y_b+2*1*I_NAI_G4z_Px;
  abcd[2760] = 4.0E0*I_NAI_G4x_Hx3yz_bb-2.0E0*1*I_NAI_G4x_Fxyz_b-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2761] = 4.0E0*I_NAI_G3xy_Hx3yz_bb-2.0E0*1*I_NAI_G3xy_Fxyz_b-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2762] = 4.0E0*I_NAI_G3xz_Hx3yz_bb-2.0E0*1*I_NAI_G3xz_Fxyz_b-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2763] = 4.0E0*I_NAI_G2x2y_Hx3yz_bb-2.0E0*1*I_NAI_G2x2y_Fxyz_b-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2764] = 4.0E0*I_NAI_G2xyz_Hx3yz_bb-2.0E0*1*I_NAI_G2xyz_Fxyz_b-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2765] = 4.0E0*I_NAI_G2x2z_Hx3yz_bb-2.0E0*1*I_NAI_G2x2z_Fxyz_b-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2766] = 4.0E0*I_NAI_Gx3y_Hx3yz_bb-2.0E0*1*I_NAI_Gx3y_Fxyz_b-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2767] = 4.0E0*I_NAI_Gx2yz_Hx3yz_bb-2.0E0*1*I_NAI_Gx2yz_Fxyz_b-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2768] = 4.0E0*I_NAI_Gxy2z_Hx3yz_bb-2.0E0*1*I_NAI_Gxy2z_Fxyz_b-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2769] = 4.0E0*I_NAI_Gx3z_Hx3yz_bb-2.0E0*1*I_NAI_Gx3z_Fxyz_b-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2770] = 4.0E0*I_NAI_G4y_Hx3yz_bb-2.0E0*1*I_NAI_G4y_Fxyz_b-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2771] = 4.0E0*I_NAI_G3yz_Hx3yz_bb-2.0E0*1*I_NAI_G3yz_Fxyz_b-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2772] = 4.0E0*I_NAI_G2y2z_Hx3yz_bb-2.0E0*1*I_NAI_G2y2z_Fxyz_b-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2773] = 4.0E0*I_NAI_Gy3z_Hx3yz_bb-2.0E0*1*I_NAI_Gy3z_Fxyz_b-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2774] = 4.0E0*I_NAI_G4z_Hx3yz_bb-2.0E0*1*I_NAI_G4z_Fxyz_b-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2775] = 4.0E0*I_NAI_G4x_Hx2y2z_bb-2.0E0*1*I_NAI_G4x_Fx2z_b;
  abcd[2776] = 4.0E0*I_NAI_G3xy_Hx2y2z_bb-2.0E0*1*I_NAI_G3xy_Fx2z_b;
  abcd[2777] = 4.0E0*I_NAI_G3xz_Hx2y2z_bb-2.0E0*1*I_NAI_G3xz_Fx2z_b;
  abcd[2778] = 4.0E0*I_NAI_G2x2y_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2z_b;
  abcd[2779] = 4.0E0*I_NAI_G2xyz_Hx2y2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2z_b;
  abcd[2780] = 4.0E0*I_NAI_G2x2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2z_b;
  abcd[2781] = 4.0E0*I_NAI_Gx3y_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2z_b;
  abcd[2782] = 4.0E0*I_NAI_Gx2yz_Hx2y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2z_b;
  abcd[2783] = 4.0E0*I_NAI_Gxy2z_Hx2y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2z_b;
  abcd[2784] = 4.0E0*I_NAI_Gx3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2z_b;
  abcd[2785] = 4.0E0*I_NAI_G4y_Hx2y2z_bb-2.0E0*1*I_NAI_G4y_Fx2z_b;
  abcd[2786] = 4.0E0*I_NAI_G3yz_Hx2y2z_bb-2.0E0*1*I_NAI_G3yz_Fx2z_b;
  abcd[2787] = 4.0E0*I_NAI_G2y2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2z_b;
  abcd[2788] = 4.0E0*I_NAI_Gy3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2z_b;
  abcd[2789] = 4.0E0*I_NAI_G4z_Hx2y2z_bb-2.0E0*1*I_NAI_G4z_Fx2z_b;
  abcd[2790] = 4.0E0*I_NAI_G4x_H5y_bb-2.0E0*3*I_NAI_G4x_F3y_b-2.0E0*4*I_NAI_G4x_F3y_b+3*2*I_NAI_G4x_Py;
  abcd[2791] = 4.0E0*I_NAI_G3xy_H5y_bb-2.0E0*3*I_NAI_G3xy_F3y_b-2.0E0*4*I_NAI_G3xy_F3y_b+3*2*I_NAI_G3xy_Py;
  abcd[2792] = 4.0E0*I_NAI_G3xz_H5y_bb-2.0E0*3*I_NAI_G3xz_F3y_b-2.0E0*4*I_NAI_G3xz_F3y_b+3*2*I_NAI_G3xz_Py;
  abcd[2793] = 4.0E0*I_NAI_G2x2y_H5y_bb-2.0E0*3*I_NAI_G2x2y_F3y_b-2.0E0*4*I_NAI_G2x2y_F3y_b+3*2*I_NAI_G2x2y_Py;
  abcd[2794] = 4.0E0*I_NAI_G2xyz_H5y_bb-2.0E0*3*I_NAI_G2xyz_F3y_b-2.0E0*4*I_NAI_G2xyz_F3y_b+3*2*I_NAI_G2xyz_Py;
  abcd[2795] = 4.0E0*I_NAI_G2x2z_H5y_bb-2.0E0*3*I_NAI_G2x2z_F3y_b-2.0E0*4*I_NAI_G2x2z_F3y_b+3*2*I_NAI_G2x2z_Py;
  abcd[2796] = 4.0E0*I_NAI_Gx3y_H5y_bb-2.0E0*3*I_NAI_Gx3y_F3y_b-2.0E0*4*I_NAI_Gx3y_F3y_b+3*2*I_NAI_Gx3y_Py;
  abcd[2797] = 4.0E0*I_NAI_Gx2yz_H5y_bb-2.0E0*3*I_NAI_Gx2yz_F3y_b-2.0E0*4*I_NAI_Gx2yz_F3y_b+3*2*I_NAI_Gx2yz_Py;
  abcd[2798] = 4.0E0*I_NAI_Gxy2z_H5y_bb-2.0E0*3*I_NAI_Gxy2z_F3y_b-2.0E0*4*I_NAI_Gxy2z_F3y_b+3*2*I_NAI_Gxy2z_Py;
  abcd[2799] = 4.0E0*I_NAI_Gx3z_H5y_bb-2.0E0*3*I_NAI_Gx3z_F3y_b-2.0E0*4*I_NAI_Gx3z_F3y_b+3*2*I_NAI_Gx3z_Py;
  abcd[2800] = 4.0E0*I_NAI_G4y_H5y_bb-2.0E0*3*I_NAI_G4y_F3y_b-2.0E0*4*I_NAI_G4y_F3y_b+3*2*I_NAI_G4y_Py;
  abcd[2801] = 4.0E0*I_NAI_G3yz_H5y_bb-2.0E0*3*I_NAI_G3yz_F3y_b-2.0E0*4*I_NAI_G3yz_F3y_b+3*2*I_NAI_G3yz_Py;
  abcd[2802] = 4.0E0*I_NAI_G2y2z_H5y_bb-2.0E0*3*I_NAI_G2y2z_F3y_b-2.0E0*4*I_NAI_G2y2z_F3y_b+3*2*I_NAI_G2y2z_Py;
  abcd[2803] = 4.0E0*I_NAI_Gy3z_H5y_bb-2.0E0*3*I_NAI_Gy3z_F3y_b-2.0E0*4*I_NAI_Gy3z_F3y_b+3*2*I_NAI_Gy3z_Py;
  abcd[2804] = 4.0E0*I_NAI_G4z_H5y_bb-2.0E0*3*I_NAI_G4z_F3y_b-2.0E0*4*I_NAI_G4z_F3y_b+3*2*I_NAI_G4z_Py;
  abcd[2805] = 4.0E0*I_NAI_G4x_H4yz_bb-2.0E0*2*I_NAI_G4x_F2yz_b-2.0E0*3*I_NAI_G4x_F2yz_b+2*1*I_NAI_G4x_Pz;
  abcd[2806] = 4.0E0*I_NAI_G3xy_H4yz_bb-2.0E0*2*I_NAI_G3xy_F2yz_b-2.0E0*3*I_NAI_G3xy_F2yz_b+2*1*I_NAI_G3xy_Pz;
  abcd[2807] = 4.0E0*I_NAI_G3xz_H4yz_bb-2.0E0*2*I_NAI_G3xz_F2yz_b-2.0E0*3*I_NAI_G3xz_F2yz_b+2*1*I_NAI_G3xz_Pz;
  abcd[2808] = 4.0E0*I_NAI_G2x2y_H4yz_bb-2.0E0*2*I_NAI_G2x2y_F2yz_b-2.0E0*3*I_NAI_G2x2y_F2yz_b+2*1*I_NAI_G2x2y_Pz;
  abcd[2809] = 4.0E0*I_NAI_G2xyz_H4yz_bb-2.0E0*2*I_NAI_G2xyz_F2yz_b-2.0E0*3*I_NAI_G2xyz_F2yz_b+2*1*I_NAI_G2xyz_Pz;
  abcd[2810] = 4.0E0*I_NAI_G2x2z_H4yz_bb-2.0E0*2*I_NAI_G2x2z_F2yz_b-2.0E0*3*I_NAI_G2x2z_F2yz_b+2*1*I_NAI_G2x2z_Pz;
  abcd[2811] = 4.0E0*I_NAI_Gx3y_H4yz_bb-2.0E0*2*I_NAI_Gx3y_F2yz_b-2.0E0*3*I_NAI_Gx3y_F2yz_b+2*1*I_NAI_Gx3y_Pz;
  abcd[2812] = 4.0E0*I_NAI_Gx2yz_H4yz_bb-2.0E0*2*I_NAI_Gx2yz_F2yz_b-2.0E0*3*I_NAI_Gx2yz_F2yz_b+2*1*I_NAI_Gx2yz_Pz;
  abcd[2813] = 4.0E0*I_NAI_Gxy2z_H4yz_bb-2.0E0*2*I_NAI_Gxy2z_F2yz_b-2.0E0*3*I_NAI_Gxy2z_F2yz_b+2*1*I_NAI_Gxy2z_Pz;
  abcd[2814] = 4.0E0*I_NAI_Gx3z_H4yz_bb-2.0E0*2*I_NAI_Gx3z_F2yz_b-2.0E0*3*I_NAI_Gx3z_F2yz_b+2*1*I_NAI_Gx3z_Pz;
  abcd[2815] = 4.0E0*I_NAI_G4y_H4yz_bb-2.0E0*2*I_NAI_G4y_F2yz_b-2.0E0*3*I_NAI_G4y_F2yz_b+2*1*I_NAI_G4y_Pz;
  abcd[2816] = 4.0E0*I_NAI_G3yz_H4yz_bb-2.0E0*2*I_NAI_G3yz_F2yz_b-2.0E0*3*I_NAI_G3yz_F2yz_b+2*1*I_NAI_G3yz_Pz;
  abcd[2817] = 4.0E0*I_NAI_G2y2z_H4yz_bb-2.0E0*2*I_NAI_G2y2z_F2yz_b-2.0E0*3*I_NAI_G2y2z_F2yz_b+2*1*I_NAI_G2y2z_Pz;
  abcd[2818] = 4.0E0*I_NAI_Gy3z_H4yz_bb-2.0E0*2*I_NAI_Gy3z_F2yz_b-2.0E0*3*I_NAI_Gy3z_F2yz_b+2*1*I_NAI_Gy3z_Pz;
  abcd[2819] = 4.0E0*I_NAI_G4z_H4yz_bb-2.0E0*2*I_NAI_G4z_F2yz_b-2.0E0*3*I_NAI_G4z_F2yz_b+2*1*I_NAI_G4z_Pz;
  abcd[2820] = 4.0E0*I_NAI_G4x_H3y2z_bb-2.0E0*1*I_NAI_G4x_Fy2z_b-2.0E0*2*I_NAI_G4x_Fy2z_b;
  abcd[2821] = 4.0E0*I_NAI_G3xy_H3y2z_bb-2.0E0*1*I_NAI_G3xy_Fy2z_b-2.0E0*2*I_NAI_G3xy_Fy2z_b;
  abcd[2822] = 4.0E0*I_NAI_G3xz_H3y2z_bb-2.0E0*1*I_NAI_G3xz_Fy2z_b-2.0E0*2*I_NAI_G3xz_Fy2z_b;
  abcd[2823] = 4.0E0*I_NAI_G2x2y_H3y2z_bb-2.0E0*1*I_NAI_G2x2y_Fy2z_b-2.0E0*2*I_NAI_G2x2y_Fy2z_b;
  abcd[2824] = 4.0E0*I_NAI_G2xyz_H3y2z_bb-2.0E0*1*I_NAI_G2xyz_Fy2z_b-2.0E0*2*I_NAI_G2xyz_Fy2z_b;
  abcd[2825] = 4.0E0*I_NAI_G2x2z_H3y2z_bb-2.0E0*1*I_NAI_G2x2z_Fy2z_b-2.0E0*2*I_NAI_G2x2z_Fy2z_b;
  abcd[2826] = 4.0E0*I_NAI_Gx3y_H3y2z_bb-2.0E0*1*I_NAI_Gx3y_Fy2z_b-2.0E0*2*I_NAI_Gx3y_Fy2z_b;
  abcd[2827] = 4.0E0*I_NAI_Gx2yz_H3y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fy2z_b-2.0E0*2*I_NAI_Gx2yz_Fy2z_b;
  abcd[2828] = 4.0E0*I_NAI_Gxy2z_H3y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fy2z_b-2.0E0*2*I_NAI_Gxy2z_Fy2z_b;
  abcd[2829] = 4.0E0*I_NAI_Gx3z_H3y2z_bb-2.0E0*1*I_NAI_Gx3z_Fy2z_b-2.0E0*2*I_NAI_Gx3z_Fy2z_b;
  abcd[2830] = 4.0E0*I_NAI_G4y_H3y2z_bb-2.0E0*1*I_NAI_G4y_Fy2z_b-2.0E0*2*I_NAI_G4y_Fy2z_b;
  abcd[2831] = 4.0E0*I_NAI_G3yz_H3y2z_bb-2.0E0*1*I_NAI_G3yz_Fy2z_b-2.0E0*2*I_NAI_G3yz_Fy2z_b;
  abcd[2832] = 4.0E0*I_NAI_G2y2z_H3y2z_bb-2.0E0*1*I_NAI_G2y2z_Fy2z_b-2.0E0*2*I_NAI_G2y2z_Fy2z_b;
  abcd[2833] = 4.0E0*I_NAI_Gy3z_H3y2z_bb-2.0E0*1*I_NAI_Gy3z_Fy2z_b-2.0E0*2*I_NAI_Gy3z_Fy2z_b;
  abcd[2834] = 4.0E0*I_NAI_G4z_H3y2z_bb-2.0E0*1*I_NAI_G4z_Fy2z_b-2.0E0*2*I_NAI_G4z_Fy2z_b;
  abcd[2835] = 4.0E0*I_NAI_G4x_H2y3z_bb-2.0E0*1*I_NAI_G4x_F3z_b;
  abcd[2836] = 4.0E0*I_NAI_G3xy_H2y3z_bb-2.0E0*1*I_NAI_G3xy_F3z_b;
  abcd[2837] = 4.0E0*I_NAI_G3xz_H2y3z_bb-2.0E0*1*I_NAI_G3xz_F3z_b;
  abcd[2838] = 4.0E0*I_NAI_G2x2y_H2y3z_bb-2.0E0*1*I_NAI_G2x2y_F3z_b;
  abcd[2839] = 4.0E0*I_NAI_G2xyz_H2y3z_bb-2.0E0*1*I_NAI_G2xyz_F3z_b;
  abcd[2840] = 4.0E0*I_NAI_G2x2z_H2y3z_bb-2.0E0*1*I_NAI_G2x2z_F3z_b;
  abcd[2841] = 4.0E0*I_NAI_Gx3y_H2y3z_bb-2.0E0*1*I_NAI_Gx3y_F3z_b;
  abcd[2842] = 4.0E0*I_NAI_Gx2yz_H2y3z_bb-2.0E0*1*I_NAI_Gx2yz_F3z_b;
  abcd[2843] = 4.0E0*I_NAI_Gxy2z_H2y3z_bb-2.0E0*1*I_NAI_Gxy2z_F3z_b;
  abcd[2844] = 4.0E0*I_NAI_Gx3z_H2y3z_bb-2.0E0*1*I_NAI_Gx3z_F3z_b;
  abcd[2845] = 4.0E0*I_NAI_G4y_H2y3z_bb-2.0E0*1*I_NAI_G4y_F3z_b;
  abcd[2846] = 4.0E0*I_NAI_G3yz_H2y3z_bb-2.0E0*1*I_NAI_G3yz_F3z_b;
  abcd[2847] = 4.0E0*I_NAI_G2y2z_H2y3z_bb-2.0E0*1*I_NAI_G2y2z_F3z_b;
  abcd[2848] = 4.0E0*I_NAI_Gy3z_H2y3z_bb-2.0E0*1*I_NAI_Gy3z_F3z_b;
  abcd[2849] = 4.0E0*I_NAI_G4z_H2y3z_bb-2.0E0*1*I_NAI_G4z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[2850] = 4.0E0*I_NAI_G4x_H3xyz_bb;
  abcd[2851] = 4.0E0*I_NAI_G3xy_H3xyz_bb;
  abcd[2852] = 4.0E0*I_NAI_G3xz_H3xyz_bb;
  abcd[2853] = 4.0E0*I_NAI_G2x2y_H3xyz_bb;
  abcd[2854] = 4.0E0*I_NAI_G2xyz_H3xyz_bb;
  abcd[2855] = 4.0E0*I_NAI_G2x2z_H3xyz_bb;
  abcd[2856] = 4.0E0*I_NAI_Gx3y_H3xyz_bb;
  abcd[2857] = 4.0E0*I_NAI_Gx2yz_H3xyz_bb;
  abcd[2858] = 4.0E0*I_NAI_Gxy2z_H3xyz_bb;
  abcd[2859] = 4.0E0*I_NAI_Gx3z_H3xyz_bb;
  abcd[2860] = 4.0E0*I_NAI_G4y_H3xyz_bb;
  abcd[2861] = 4.0E0*I_NAI_G3yz_H3xyz_bb;
  abcd[2862] = 4.0E0*I_NAI_G2y2z_H3xyz_bb;
  abcd[2863] = 4.0E0*I_NAI_Gy3z_H3xyz_bb;
  abcd[2864] = 4.0E0*I_NAI_G4z_H3xyz_bb;
  abcd[2865] = 4.0E0*I_NAI_G4x_H2x2yz_bb-2.0E0*1*I_NAI_G4x_F2xz_b;
  abcd[2866] = 4.0E0*I_NAI_G3xy_H2x2yz_bb-2.0E0*1*I_NAI_G3xy_F2xz_b;
  abcd[2867] = 4.0E0*I_NAI_G3xz_H2x2yz_bb-2.0E0*1*I_NAI_G3xz_F2xz_b;
  abcd[2868] = 4.0E0*I_NAI_G2x2y_H2x2yz_bb-2.0E0*1*I_NAI_G2x2y_F2xz_b;
  abcd[2869] = 4.0E0*I_NAI_G2xyz_H2x2yz_bb-2.0E0*1*I_NAI_G2xyz_F2xz_b;
  abcd[2870] = 4.0E0*I_NAI_G2x2z_H2x2yz_bb-2.0E0*1*I_NAI_G2x2z_F2xz_b;
  abcd[2871] = 4.0E0*I_NAI_Gx3y_H2x2yz_bb-2.0E0*1*I_NAI_Gx3y_F2xz_b;
  abcd[2872] = 4.0E0*I_NAI_Gx2yz_H2x2yz_bb-2.0E0*1*I_NAI_Gx2yz_F2xz_b;
  abcd[2873] = 4.0E0*I_NAI_Gxy2z_H2x2yz_bb-2.0E0*1*I_NAI_Gxy2z_F2xz_b;
  abcd[2874] = 4.0E0*I_NAI_Gx3z_H2x2yz_bb-2.0E0*1*I_NAI_Gx3z_F2xz_b;
  abcd[2875] = 4.0E0*I_NAI_G4y_H2x2yz_bb-2.0E0*1*I_NAI_G4y_F2xz_b;
  abcd[2876] = 4.0E0*I_NAI_G3yz_H2x2yz_bb-2.0E0*1*I_NAI_G3yz_F2xz_b;
  abcd[2877] = 4.0E0*I_NAI_G2y2z_H2x2yz_bb-2.0E0*1*I_NAI_G2y2z_F2xz_b;
  abcd[2878] = 4.0E0*I_NAI_Gy3z_H2x2yz_bb-2.0E0*1*I_NAI_Gy3z_F2xz_b;
  abcd[2879] = 4.0E0*I_NAI_G4z_H2x2yz_bb-2.0E0*1*I_NAI_G4z_F2xz_b;
  abcd[2880] = 4.0E0*I_NAI_G4x_H2xy2z_bb-2.0E0*1*I_NAI_G4x_F2xy_b;
  abcd[2881] = 4.0E0*I_NAI_G3xy_H2xy2z_bb-2.0E0*1*I_NAI_G3xy_F2xy_b;
  abcd[2882] = 4.0E0*I_NAI_G3xz_H2xy2z_bb-2.0E0*1*I_NAI_G3xz_F2xy_b;
  abcd[2883] = 4.0E0*I_NAI_G2x2y_H2xy2z_bb-2.0E0*1*I_NAI_G2x2y_F2xy_b;
  abcd[2884] = 4.0E0*I_NAI_G2xyz_H2xy2z_bb-2.0E0*1*I_NAI_G2xyz_F2xy_b;
  abcd[2885] = 4.0E0*I_NAI_G2x2z_H2xy2z_bb-2.0E0*1*I_NAI_G2x2z_F2xy_b;
  abcd[2886] = 4.0E0*I_NAI_Gx3y_H2xy2z_bb-2.0E0*1*I_NAI_Gx3y_F2xy_b;
  abcd[2887] = 4.0E0*I_NAI_Gx2yz_H2xy2z_bb-2.0E0*1*I_NAI_Gx2yz_F2xy_b;
  abcd[2888] = 4.0E0*I_NAI_Gxy2z_H2xy2z_bb-2.0E0*1*I_NAI_Gxy2z_F2xy_b;
  abcd[2889] = 4.0E0*I_NAI_Gx3z_H2xy2z_bb-2.0E0*1*I_NAI_Gx3z_F2xy_b;
  abcd[2890] = 4.0E0*I_NAI_G4y_H2xy2z_bb-2.0E0*1*I_NAI_G4y_F2xy_b;
  abcd[2891] = 4.0E0*I_NAI_G3yz_H2xy2z_bb-2.0E0*1*I_NAI_G3yz_F2xy_b;
  abcd[2892] = 4.0E0*I_NAI_G2y2z_H2xy2z_bb-2.0E0*1*I_NAI_G2y2z_F2xy_b;
  abcd[2893] = 4.0E0*I_NAI_Gy3z_H2xy2z_bb-2.0E0*1*I_NAI_Gy3z_F2xy_b;
  abcd[2894] = 4.0E0*I_NAI_G4z_H2xy2z_bb-2.0E0*1*I_NAI_G4z_F2xy_b;
  abcd[2895] = 4.0E0*I_NAI_G4x_Hx3yz_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2896] = 4.0E0*I_NAI_G3xy_Hx3yz_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2897] = 4.0E0*I_NAI_G3xz_Hx3yz_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2898] = 4.0E0*I_NAI_G2x2y_Hx3yz_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2899] = 4.0E0*I_NAI_G2xyz_Hx3yz_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2900] = 4.0E0*I_NAI_G2x2z_Hx3yz_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2901] = 4.0E0*I_NAI_Gx3y_Hx3yz_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2902] = 4.0E0*I_NAI_Gx2yz_Hx3yz_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2903] = 4.0E0*I_NAI_Gxy2z_Hx3yz_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2904] = 4.0E0*I_NAI_Gx3z_Hx3yz_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2905] = 4.0E0*I_NAI_G4y_Hx3yz_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2906] = 4.0E0*I_NAI_G3yz_Hx3yz_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2907] = 4.0E0*I_NAI_G2y2z_Hx3yz_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2908] = 4.0E0*I_NAI_Gy3z_Hx3yz_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2909] = 4.0E0*I_NAI_G4z_Hx3yz_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2910] = 4.0E0*I_NAI_G4x_Hx2y2z_bb-2.0E0*1*I_NAI_G4x_Fx2y_b-2.0E0*1*I_NAI_G4x_Fx2z_b+1*I_NAI_G4x_Px;
  abcd[2911] = 4.0E0*I_NAI_G3xy_Hx2y2z_bb-2.0E0*1*I_NAI_G3xy_Fx2y_b-2.0E0*1*I_NAI_G3xy_Fx2z_b+1*I_NAI_G3xy_Px;
  abcd[2912] = 4.0E0*I_NAI_G3xz_Hx2y2z_bb-2.0E0*1*I_NAI_G3xz_Fx2y_b-2.0E0*1*I_NAI_G3xz_Fx2z_b+1*I_NAI_G3xz_Px;
  abcd[2913] = 4.0E0*I_NAI_G2x2y_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2y_b-2.0E0*1*I_NAI_G2x2y_Fx2z_b+1*I_NAI_G2x2y_Px;
  abcd[2914] = 4.0E0*I_NAI_G2xyz_Hx2y2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2y_b-2.0E0*1*I_NAI_G2xyz_Fx2z_b+1*I_NAI_G2xyz_Px;
  abcd[2915] = 4.0E0*I_NAI_G2x2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2y_b-2.0E0*1*I_NAI_G2x2z_Fx2z_b+1*I_NAI_G2x2z_Px;
  abcd[2916] = 4.0E0*I_NAI_Gx3y_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2y_b-2.0E0*1*I_NAI_Gx3y_Fx2z_b+1*I_NAI_Gx3y_Px;
  abcd[2917] = 4.0E0*I_NAI_Gx2yz_Hx2y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2y_b-2.0E0*1*I_NAI_Gx2yz_Fx2z_b+1*I_NAI_Gx2yz_Px;
  abcd[2918] = 4.0E0*I_NAI_Gxy2z_Hx2y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2y_b-2.0E0*1*I_NAI_Gxy2z_Fx2z_b+1*I_NAI_Gxy2z_Px;
  abcd[2919] = 4.0E0*I_NAI_Gx3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2y_b-2.0E0*1*I_NAI_Gx3z_Fx2z_b+1*I_NAI_Gx3z_Px;
  abcd[2920] = 4.0E0*I_NAI_G4y_Hx2y2z_bb-2.0E0*1*I_NAI_G4y_Fx2y_b-2.0E0*1*I_NAI_G4y_Fx2z_b+1*I_NAI_G4y_Px;
  abcd[2921] = 4.0E0*I_NAI_G3yz_Hx2y2z_bb-2.0E0*1*I_NAI_G3yz_Fx2y_b-2.0E0*1*I_NAI_G3yz_Fx2z_b+1*I_NAI_G3yz_Px;
  abcd[2922] = 4.0E0*I_NAI_G2y2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2y_b-2.0E0*1*I_NAI_G2y2z_Fx2z_b+1*I_NAI_G2y2z_Px;
  abcd[2923] = 4.0E0*I_NAI_Gy3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2y_b-2.0E0*1*I_NAI_Gy3z_Fx2z_b+1*I_NAI_Gy3z_Px;
  abcd[2924] = 4.0E0*I_NAI_G4z_Hx2y2z_bb-2.0E0*1*I_NAI_G4z_Fx2y_b-2.0E0*1*I_NAI_G4z_Fx2z_b+1*I_NAI_G4z_Px;
  abcd[2925] = 4.0E0*I_NAI_G4x_Hxy3z_bb-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[2926] = 4.0E0*I_NAI_G3xy_Hxy3z_bb-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[2927] = 4.0E0*I_NAI_G3xz_Hxy3z_bb-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[2928] = 4.0E0*I_NAI_G2x2y_Hxy3z_bb-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[2929] = 4.0E0*I_NAI_G2xyz_Hxy3z_bb-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[2930] = 4.0E0*I_NAI_G2x2z_Hxy3z_bb-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[2931] = 4.0E0*I_NAI_Gx3y_Hxy3z_bb-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[2932] = 4.0E0*I_NAI_Gx2yz_Hxy3z_bb-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[2933] = 4.0E0*I_NAI_Gxy2z_Hxy3z_bb-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[2934] = 4.0E0*I_NAI_Gx3z_Hxy3z_bb-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[2935] = 4.0E0*I_NAI_G4y_Hxy3z_bb-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[2936] = 4.0E0*I_NAI_G3yz_Hxy3z_bb-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[2937] = 4.0E0*I_NAI_G2y2z_Hxy3z_bb-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[2938] = 4.0E0*I_NAI_Gy3z_Hxy3z_bb-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[2939] = 4.0E0*I_NAI_G4z_Hxy3z_bb-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[2940] = 4.0E0*I_NAI_G4x_H4yz_bb-2.0E0*3*I_NAI_G4x_F2yz_b;
  abcd[2941] = 4.0E0*I_NAI_G3xy_H4yz_bb-2.0E0*3*I_NAI_G3xy_F2yz_b;
  abcd[2942] = 4.0E0*I_NAI_G3xz_H4yz_bb-2.0E0*3*I_NAI_G3xz_F2yz_b;
  abcd[2943] = 4.0E0*I_NAI_G2x2y_H4yz_bb-2.0E0*3*I_NAI_G2x2y_F2yz_b;
  abcd[2944] = 4.0E0*I_NAI_G2xyz_H4yz_bb-2.0E0*3*I_NAI_G2xyz_F2yz_b;
  abcd[2945] = 4.0E0*I_NAI_G2x2z_H4yz_bb-2.0E0*3*I_NAI_G2x2z_F2yz_b;
  abcd[2946] = 4.0E0*I_NAI_Gx3y_H4yz_bb-2.0E0*3*I_NAI_Gx3y_F2yz_b;
  abcd[2947] = 4.0E0*I_NAI_Gx2yz_H4yz_bb-2.0E0*3*I_NAI_Gx2yz_F2yz_b;
  abcd[2948] = 4.0E0*I_NAI_Gxy2z_H4yz_bb-2.0E0*3*I_NAI_Gxy2z_F2yz_b;
  abcd[2949] = 4.0E0*I_NAI_Gx3z_H4yz_bb-2.0E0*3*I_NAI_Gx3z_F2yz_b;
  abcd[2950] = 4.0E0*I_NAI_G4y_H4yz_bb-2.0E0*3*I_NAI_G4y_F2yz_b;
  abcd[2951] = 4.0E0*I_NAI_G3yz_H4yz_bb-2.0E0*3*I_NAI_G3yz_F2yz_b;
  abcd[2952] = 4.0E0*I_NAI_G2y2z_H4yz_bb-2.0E0*3*I_NAI_G2y2z_F2yz_b;
  abcd[2953] = 4.0E0*I_NAI_Gy3z_H4yz_bb-2.0E0*3*I_NAI_Gy3z_F2yz_b;
  abcd[2954] = 4.0E0*I_NAI_G4z_H4yz_bb-2.0E0*3*I_NAI_G4z_F2yz_b;
  abcd[2955] = 4.0E0*I_NAI_G4x_H3y2z_bb-2.0E0*1*I_NAI_G4x_F3y_b-2.0E0*2*I_NAI_G4x_Fy2z_b+2*1*I_NAI_G4x_Py;
  abcd[2956] = 4.0E0*I_NAI_G3xy_H3y2z_bb-2.0E0*1*I_NAI_G3xy_F3y_b-2.0E0*2*I_NAI_G3xy_Fy2z_b+2*1*I_NAI_G3xy_Py;
  abcd[2957] = 4.0E0*I_NAI_G3xz_H3y2z_bb-2.0E0*1*I_NAI_G3xz_F3y_b-2.0E0*2*I_NAI_G3xz_Fy2z_b+2*1*I_NAI_G3xz_Py;
  abcd[2958] = 4.0E0*I_NAI_G2x2y_H3y2z_bb-2.0E0*1*I_NAI_G2x2y_F3y_b-2.0E0*2*I_NAI_G2x2y_Fy2z_b+2*1*I_NAI_G2x2y_Py;
  abcd[2959] = 4.0E0*I_NAI_G2xyz_H3y2z_bb-2.0E0*1*I_NAI_G2xyz_F3y_b-2.0E0*2*I_NAI_G2xyz_Fy2z_b+2*1*I_NAI_G2xyz_Py;
  abcd[2960] = 4.0E0*I_NAI_G2x2z_H3y2z_bb-2.0E0*1*I_NAI_G2x2z_F3y_b-2.0E0*2*I_NAI_G2x2z_Fy2z_b+2*1*I_NAI_G2x2z_Py;
  abcd[2961] = 4.0E0*I_NAI_Gx3y_H3y2z_bb-2.0E0*1*I_NAI_Gx3y_F3y_b-2.0E0*2*I_NAI_Gx3y_Fy2z_b+2*1*I_NAI_Gx3y_Py;
  abcd[2962] = 4.0E0*I_NAI_Gx2yz_H3y2z_bb-2.0E0*1*I_NAI_Gx2yz_F3y_b-2.0E0*2*I_NAI_Gx2yz_Fy2z_b+2*1*I_NAI_Gx2yz_Py;
  abcd[2963] = 4.0E0*I_NAI_Gxy2z_H3y2z_bb-2.0E0*1*I_NAI_Gxy2z_F3y_b-2.0E0*2*I_NAI_Gxy2z_Fy2z_b+2*1*I_NAI_Gxy2z_Py;
  abcd[2964] = 4.0E0*I_NAI_Gx3z_H3y2z_bb-2.0E0*1*I_NAI_Gx3z_F3y_b-2.0E0*2*I_NAI_Gx3z_Fy2z_b+2*1*I_NAI_Gx3z_Py;
  abcd[2965] = 4.0E0*I_NAI_G4y_H3y2z_bb-2.0E0*1*I_NAI_G4y_F3y_b-2.0E0*2*I_NAI_G4y_Fy2z_b+2*1*I_NAI_G4y_Py;
  abcd[2966] = 4.0E0*I_NAI_G3yz_H3y2z_bb-2.0E0*1*I_NAI_G3yz_F3y_b-2.0E0*2*I_NAI_G3yz_Fy2z_b+2*1*I_NAI_G3yz_Py;
  abcd[2967] = 4.0E0*I_NAI_G2y2z_H3y2z_bb-2.0E0*1*I_NAI_G2y2z_F3y_b-2.0E0*2*I_NAI_G2y2z_Fy2z_b+2*1*I_NAI_G2y2z_Py;
  abcd[2968] = 4.0E0*I_NAI_Gy3z_H3y2z_bb-2.0E0*1*I_NAI_Gy3z_F3y_b-2.0E0*2*I_NAI_Gy3z_Fy2z_b+2*1*I_NAI_Gy3z_Py;
  abcd[2969] = 4.0E0*I_NAI_G4z_H3y2z_bb-2.0E0*1*I_NAI_G4z_F3y_b-2.0E0*2*I_NAI_G4z_Fy2z_b+2*1*I_NAI_G4z_Py;
  abcd[2970] = 4.0E0*I_NAI_G4x_H2y3z_bb-2.0E0*2*I_NAI_G4x_F2yz_b-2.0E0*1*I_NAI_G4x_F3z_b+2*I_NAI_G4x_Pz;
  abcd[2971] = 4.0E0*I_NAI_G3xy_H2y3z_bb-2.0E0*2*I_NAI_G3xy_F2yz_b-2.0E0*1*I_NAI_G3xy_F3z_b+2*I_NAI_G3xy_Pz;
  abcd[2972] = 4.0E0*I_NAI_G3xz_H2y3z_bb-2.0E0*2*I_NAI_G3xz_F2yz_b-2.0E0*1*I_NAI_G3xz_F3z_b+2*I_NAI_G3xz_Pz;
  abcd[2973] = 4.0E0*I_NAI_G2x2y_H2y3z_bb-2.0E0*2*I_NAI_G2x2y_F2yz_b-2.0E0*1*I_NAI_G2x2y_F3z_b+2*I_NAI_G2x2y_Pz;
  abcd[2974] = 4.0E0*I_NAI_G2xyz_H2y3z_bb-2.0E0*2*I_NAI_G2xyz_F2yz_b-2.0E0*1*I_NAI_G2xyz_F3z_b+2*I_NAI_G2xyz_Pz;
  abcd[2975] = 4.0E0*I_NAI_G2x2z_H2y3z_bb-2.0E0*2*I_NAI_G2x2z_F2yz_b-2.0E0*1*I_NAI_G2x2z_F3z_b+2*I_NAI_G2x2z_Pz;
  abcd[2976] = 4.0E0*I_NAI_Gx3y_H2y3z_bb-2.0E0*2*I_NAI_Gx3y_F2yz_b-2.0E0*1*I_NAI_Gx3y_F3z_b+2*I_NAI_Gx3y_Pz;
  abcd[2977] = 4.0E0*I_NAI_Gx2yz_H2y3z_bb-2.0E0*2*I_NAI_Gx2yz_F2yz_b-2.0E0*1*I_NAI_Gx2yz_F3z_b+2*I_NAI_Gx2yz_Pz;
  abcd[2978] = 4.0E0*I_NAI_Gxy2z_H2y3z_bb-2.0E0*2*I_NAI_Gxy2z_F2yz_b-2.0E0*1*I_NAI_Gxy2z_F3z_b+2*I_NAI_Gxy2z_Pz;
  abcd[2979] = 4.0E0*I_NAI_Gx3z_H2y3z_bb-2.0E0*2*I_NAI_Gx3z_F2yz_b-2.0E0*1*I_NAI_Gx3z_F3z_b+2*I_NAI_Gx3z_Pz;
  abcd[2980] = 4.0E0*I_NAI_G4y_H2y3z_bb-2.0E0*2*I_NAI_G4y_F2yz_b-2.0E0*1*I_NAI_G4y_F3z_b+2*I_NAI_G4y_Pz;
  abcd[2981] = 4.0E0*I_NAI_G3yz_H2y3z_bb-2.0E0*2*I_NAI_G3yz_F2yz_b-2.0E0*1*I_NAI_G3yz_F3z_b+2*I_NAI_G3yz_Pz;
  abcd[2982] = 4.0E0*I_NAI_G2y2z_H2y3z_bb-2.0E0*2*I_NAI_G2y2z_F2yz_b-2.0E0*1*I_NAI_G2y2z_F3z_b+2*I_NAI_G2y2z_Pz;
  abcd[2983] = 4.0E0*I_NAI_Gy3z_H2y3z_bb-2.0E0*2*I_NAI_Gy3z_F2yz_b-2.0E0*1*I_NAI_Gy3z_F3z_b+2*I_NAI_Gy3z_Pz;
  abcd[2984] = 4.0E0*I_NAI_G4z_H2y3z_bb-2.0E0*2*I_NAI_G4z_F2yz_b-2.0E0*1*I_NAI_G4z_F3z_b+2*I_NAI_G4z_Pz;
  abcd[2985] = 4.0E0*I_NAI_G4x_Hy4z_bb-2.0E0*3*I_NAI_G4x_Fy2z_b;
  abcd[2986] = 4.0E0*I_NAI_G3xy_Hy4z_bb-2.0E0*3*I_NAI_G3xy_Fy2z_b;
  abcd[2987] = 4.0E0*I_NAI_G3xz_Hy4z_bb-2.0E0*3*I_NAI_G3xz_Fy2z_b;
  abcd[2988] = 4.0E0*I_NAI_G2x2y_Hy4z_bb-2.0E0*3*I_NAI_G2x2y_Fy2z_b;
  abcd[2989] = 4.0E0*I_NAI_G2xyz_Hy4z_bb-2.0E0*3*I_NAI_G2xyz_Fy2z_b;
  abcd[2990] = 4.0E0*I_NAI_G2x2z_Hy4z_bb-2.0E0*3*I_NAI_G2x2z_Fy2z_b;
  abcd[2991] = 4.0E0*I_NAI_Gx3y_Hy4z_bb-2.0E0*3*I_NAI_Gx3y_Fy2z_b;
  abcd[2992] = 4.0E0*I_NAI_Gx2yz_Hy4z_bb-2.0E0*3*I_NAI_Gx2yz_Fy2z_b;
  abcd[2993] = 4.0E0*I_NAI_Gxy2z_Hy4z_bb-2.0E0*3*I_NAI_Gxy2z_Fy2z_b;
  abcd[2994] = 4.0E0*I_NAI_Gx3z_Hy4z_bb-2.0E0*3*I_NAI_Gx3z_Fy2z_b;
  abcd[2995] = 4.0E0*I_NAI_G4y_Hy4z_bb-2.0E0*3*I_NAI_G4y_Fy2z_b;
  abcd[2996] = 4.0E0*I_NAI_G3yz_Hy4z_bb-2.0E0*3*I_NAI_G3yz_Fy2z_b;
  abcd[2997] = 4.0E0*I_NAI_G2y2z_Hy4z_bb-2.0E0*3*I_NAI_G2y2z_Fy2z_b;
  abcd[2998] = 4.0E0*I_NAI_Gy3z_Hy4z_bb-2.0E0*3*I_NAI_Gy3z_Fy2z_b;
  abcd[2999] = 4.0E0*I_NAI_G4z_Hy4z_bb-2.0E0*3*I_NAI_G4z_Fy2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_H_bb
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  abcd[3000] = 4.0E0*I_NAI_G4x_H3x2z_bb-2.0E0*1*I_NAI_G4x_F3x_b;
  abcd[3001] = 4.0E0*I_NAI_G3xy_H3x2z_bb-2.0E0*1*I_NAI_G3xy_F3x_b;
  abcd[3002] = 4.0E0*I_NAI_G3xz_H3x2z_bb-2.0E0*1*I_NAI_G3xz_F3x_b;
  abcd[3003] = 4.0E0*I_NAI_G2x2y_H3x2z_bb-2.0E0*1*I_NAI_G2x2y_F3x_b;
  abcd[3004] = 4.0E0*I_NAI_G2xyz_H3x2z_bb-2.0E0*1*I_NAI_G2xyz_F3x_b;
  abcd[3005] = 4.0E0*I_NAI_G2x2z_H3x2z_bb-2.0E0*1*I_NAI_G2x2z_F3x_b;
  abcd[3006] = 4.0E0*I_NAI_Gx3y_H3x2z_bb-2.0E0*1*I_NAI_Gx3y_F3x_b;
  abcd[3007] = 4.0E0*I_NAI_Gx2yz_H3x2z_bb-2.0E0*1*I_NAI_Gx2yz_F3x_b;
  abcd[3008] = 4.0E0*I_NAI_Gxy2z_H3x2z_bb-2.0E0*1*I_NAI_Gxy2z_F3x_b;
  abcd[3009] = 4.0E0*I_NAI_Gx3z_H3x2z_bb-2.0E0*1*I_NAI_Gx3z_F3x_b;
  abcd[3010] = 4.0E0*I_NAI_G4y_H3x2z_bb-2.0E0*1*I_NAI_G4y_F3x_b;
  abcd[3011] = 4.0E0*I_NAI_G3yz_H3x2z_bb-2.0E0*1*I_NAI_G3yz_F3x_b;
  abcd[3012] = 4.0E0*I_NAI_G2y2z_H3x2z_bb-2.0E0*1*I_NAI_G2y2z_F3x_b;
  abcd[3013] = 4.0E0*I_NAI_Gy3z_H3x2z_bb-2.0E0*1*I_NAI_Gy3z_F3x_b;
  abcd[3014] = 4.0E0*I_NAI_G4z_H3x2z_bb-2.0E0*1*I_NAI_G4z_F3x_b;
  abcd[3015] = 4.0E0*I_NAI_G4x_H2xy2z_bb-2.0E0*1*I_NAI_G4x_F2xy_b;
  abcd[3016] = 4.0E0*I_NAI_G3xy_H2xy2z_bb-2.0E0*1*I_NAI_G3xy_F2xy_b;
  abcd[3017] = 4.0E0*I_NAI_G3xz_H2xy2z_bb-2.0E0*1*I_NAI_G3xz_F2xy_b;
  abcd[3018] = 4.0E0*I_NAI_G2x2y_H2xy2z_bb-2.0E0*1*I_NAI_G2x2y_F2xy_b;
  abcd[3019] = 4.0E0*I_NAI_G2xyz_H2xy2z_bb-2.0E0*1*I_NAI_G2xyz_F2xy_b;
  abcd[3020] = 4.0E0*I_NAI_G2x2z_H2xy2z_bb-2.0E0*1*I_NAI_G2x2z_F2xy_b;
  abcd[3021] = 4.0E0*I_NAI_Gx3y_H2xy2z_bb-2.0E0*1*I_NAI_Gx3y_F2xy_b;
  abcd[3022] = 4.0E0*I_NAI_Gx2yz_H2xy2z_bb-2.0E0*1*I_NAI_Gx2yz_F2xy_b;
  abcd[3023] = 4.0E0*I_NAI_Gxy2z_H2xy2z_bb-2.0E0*1*I_NAI_Gxy2z_F2xy_b;
  abcd[3024] = 4.0E0*I_NAI_Gx3z_H2xy2z_bb-2.0E0*1*I_NAI_Gx3z_F2xy_b;
  abcd[3025] = 4.0E0*I_NAI_G4y_H2xy2z_bb-2.0E0*1*I_NAI_G4y_F2xy_b;
  abcd[3026] = 4.0E0*I_NAI_G3yz_H2xy2z_bb-2.0E0*1*I_NAI_G3yz_F2xy_b;
  abcd[3027] = 4.0E0*I_NAI_G2y2z_H2xy2z_bb-2.0E0*1*I_NAI_G2y2z_F2xy_b;
  abcd[3028] = 4.0E0*I_NAI_Gy3z_H2xy2z_bb-2.0E0*1*I_NAI_Gy3z_F2xy_b;
  abcd[3029] = 4.0E0*I_NAI_G4z_H2xy2z_bb-2.0E0*1*I_NAI_G4z_F2xy_b;
  abcd[3030] = 4.0E0*I_NAI_G4x_H2x3z_bb-2.0E0*1*I_NAI_G4x_F2xz_b-2.0E0*2*I_NAI_G4x_F2xz_b;
  abcd[3031] = 4.0E0*I_NAI_G3xy_H2x3z_bb-2.0E0*1*I_NAI_G3xy_F2xz_b-2.0E0*2*I_NAI_G3xy_F2xz_b;
  abcd[3032] = 4.0E0*I_NAI_G3xz_H2x3z_bb-2.0E0*1*I_NAI_G3xz_F2xz_b-2.0E0*2*I_NAI_G3xz_F2xz_b;
  abcd[3033] = 4.0E0*I_NAI_G2x2y_H2x3z_bb-2.0E0*1*I_NAI_G2x2y_F2xz_b-2.0E0*2*I_NAI_G2x2y_F2xz_b;
  abcd[3034] = 4.0E0*I_NAI_G2xyz_H2x3z_bb-2.0E0*1*I_NAI_G2xyz_F2xz_b-2.0E0*2*I_NAI_G2xyz_F2xz_b;
  abcd[3035] = 4.0E0*I_NAI_G2x2z_H2x3z_bb-2.0E0*1*I_NAI_G2x2z_F2xz_b-2.0E0*2*I_NAI_G2x2z_F2xz_b;
  abcd[3036] = 4.0E0*I_NAI_Gx3y_H2x3z_bb-2.0E0*1*I_NAI_Gx3y_F2xz_b-2.0E0*2*I_NAI_Gx3y_F2xz_b;
  abcd[3037] = 4.0E0*I_NAI_Gx2yz_H2x3z_bb-2.0E0*1*I_NAI_Gx2yz_F2xz_b-2.0E0*2*I_NAI_Gx2yz_F2xz_b;
  abcd[3038] = 4.0E0*I_NAI_Gxy2z_H2x3z_bb-2.0E0*1*I_NAI_Gxy2z_F2xz_b-2.0E0*2*I_NAI_Gxy2z_F2xz_b;
  abcd[3039] = 4.0E0*I_NAI_Gx3z_H2x3z_bb-2.0E0*1*I_NAI_Gx3z_F2xz_b-2.0E0*2*I_NAI_Gx3z_F2xz_b;
  abcd[3040] = 4.0E0*I_NAI_G4y_H2x3z_bb-2.0E0*1*I_NAI_G4y_F2xz_b-2.0E0*2*I_NAI_G4y_F2xz_b;
  abcd[3041] = 4.0E0*I_NAI_G3yz_H2x3z_bb-2.0E0*1*I_NAI_G3yz_F2xz_b-2.0E0*2*I_NAI_G3yz_F2xz_b;
  abcd[3042] = 4.0E0*I_NAI_G2y2z_H2x3z_bb-2.0E0*1*I_NAI_G2y2z_F2xz_b-2.0E0*2*I_NAI_G2y2z_F2xz_b;
  abcd[3043] = 4.0E0*I_NAI_Gy3z_H2x3z_bb-2.0E0*1*I_NAI_Gy3z_F2xz_b-2.0E0*2*I_NAI_Gy3z_F2xz_b;
  abcd[3044] = 4.0E0*I_NAI_G4z_H2x3z_bb-2.0E0*1*I_NAI_G4z_F2xz_b-2.0E0*2*I_NAI_G4z_F2xz_b;
  abcd[3045] = 4.0E0*I_NAI_G4x_Hx2y2z_bb-2.0E0*1*I_NAI_G4x_Fx2y_b;
  abcd[3046] = 4.0E0*I_NAI_G3xy_Hx2y2z_bb-2.0E0*1*I_NAI_G3xy_Fx2y_b;
  abcd[3047] = 4.0E0*I_NAI_G3xz_Hx2y2z_bb-2.0E0*1*I_NAI_G3xz_Fx2y_b;
  abcd[3048] = 4.0E0*I_NAI_G2x2y_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2y_Fx2y_b;
  abcd[3049] = 4.0E0*I_NAI_G2xyz_Hx2y2z_bb-2.0E0*1*I_NAI_G2xyz_Fx2y_b;
  abcd[3050] = 4.0E0*I_NAI_G2x2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2x2z_Fx2y_b;
  abcd[3051] = 4.0E0*I_NAI_Gx3y_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3y_Fx2y_b;
  abcd[3052] = 4.0E0*I_NAI_Gx2yz_Hx2y2z_bb-2.0E0*1*I_NAI_Gx2yz_Fx2y_b;
  abcd[3053] = 4.0E0*I_NAI_Gxy2z_Hx2y2z_bb-2.0E0*1*I_NAI_Gxy2z_Fx2y_b;
  abcd[3054] = 4.0E0*I_NAI_Gx3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gx3z_Fx2y_b;
  abcd[3055] = 4.0E0*I_NAI_G4y_Hx2y2z_bb-2.0E0*1*I_NAI_G4y_Fx2y_b;
  abcd[3056] = 4.0E0*I_NAI_G3yz_Hx2y2z_bb-2.0E0*1*I_NAI_G3yz_Fx2y_b;
  abcd[3057] = 4.0E0*I_NAI_G2y2z_Hx2y2z_bb-2.0E0*1*I_NAI_G2y2z_Fx2y_b;
  abcd[3058] = 4.0E0*I_NAI_Gy3z_Hx2y2z_bb-2.0E0*1*I_NAI_Gy3z_Fx2y_b;
  abcd[3059] = 4.0E0*I_NAI_G4z_Hx2y2z_bb-2.0E0*1*I_NAI_G4z_Fx2y_b;
  abcd[3060] = 4.0E0*I_NAI_G4x_Hxy3z_bb-2.0E0*1*I_NAI_G4x_Fxyz_b-2.0E0*2*I_NAI_G4x_Fxyz_b;
  abcd[3061] = 4.0E0*I_NAI_G3xy_Hxy3z_bb-2.0E0*1*I_NAI_G3xy_Fxyz_b-2.0E0*2*I_NAI_G3xy_Fxyz_b;
  abcd[3062] = 4.0E0*I_NAI_G3xz_Hxy3z_bb-2.0E0*1*I_NAI_G3xz_Fxyz_b-2.0E0*2*I_NAI_G3xz_Fxyz_b;
  abcd[3063] = 4.0E0*I_NAI_G2x2y_Hxy3z_bb-2.0E0*1*I_NAI_G2x2y_Fxyz_b-2.0E0*2*I_NAI_G2x2y_Fxyz_b;
  abcd[3064] = 4.0E0*I_NAI_G2xyz_Hxy3z_bb-2.0E0*1*I_NAI_G2xyz_Fxyz_b-2.0E0*2*I_NAI_G2xyz_Fxyz_b;
  abcd[3065] = 4.0E0*I_NAI_G2x2z_Hxy3z_bb-2.0E0*1*I_NAI_G2x2z_Fxyz_b-2.0E0*2*I_NAI_G2x2z_Fxyz_b;
  abcd[3066] = 4.0E0*I_NAI_Gx3y_Hxy3z_bb-2.0E0*1*I_NAI_Gx3y_Fxyz_b-2.0E0*2*I_NAI_Gx3y_Fxyz_b;
  abcd[3067] = 4.0E0*I_NAI_Gx2yz_Hxy3z_bb-2.0E0*1*I_NAI_Gx2yz_Fxyz_b-2.0E0*2*I_NAI_Gx2yz_Fxyz_b;
  abcd[3068] = 4.0E0*I_NAI_Gxy2z_Hxy3z_bb-2.0E0*1*I_NAI_Gxy2z_Fxyz_b-2.0E0*2*I_NAI_Gxy2z_Fxyz_b;
  abcd[3069] = 4.0E0*I_NAI_Gx3z_Hxy3z_bb-2.0E0*1*I_NAI_Gx3z_Fxyz_b-2.0E0*2*I_NAI_Gx3z_Fxyz_b;
  abcd[3070] = 4.0E0*I_NAI_G4y_Hxy3z_bb-2.0E0*1*I_NAI_G4y_Fxyz_b-2.0E0*2*I_NAI_G4y_Fxyz_b;
  abcd[3071] = 4.0E0*I_NAI_G3yz_Hxy3z_bb-2.0E0*1*I_NAI_G3yz_Fxyz_b-2.0E0*2*I_NAI_G3yz_Fxyz_b;
  abcd[3072] = 4.0E0*I_NAI_G2y2z_Hxy3z_bb-2.0E0*1*I_NAI_G2y2z_Fxyz_b-2.0E0*2*I_NAI_G2y2z_Fxyz_b;
  abcd[3073] = 4.0E0*I_NAI_Gy3z_Hxy3z_bb-2.0E0*1*I_NAI_Gy3z_Fxyz_b-2.0E0*2*I_NAI_Gy3z_Fxyz_b;
  abcd[3074] = 4.0E0*I_NAI_G4z_Hxy3z_bb-2.0E0*1*I_NAI_G4z_Fxyz_b-2.0E0*2*I_NAI_G4z_Fxyz_b;
  abcd[3075] = 4.0E0*I_NAI_G4x_Hx4z_bb-2.0E0*2*I_NAI_G4x_Fx2z_b-2.0E0*3*I_NAI_G4x_Fx2z_b+2*1*I_NAI_G4x_Px;
  abcd[3076] = 4.0E0*I_NAI_G3xy_Hx4z_bb-2.0E0*2*I_NAI_G3xy_Fx2z_b-2.0E0*3*I_NAI_G3xy_Fx2z_b+2*1*I_NAI_G3xy_Px;
  abcd[3077] = 4.0E0*I_NAI_G3xz_Hx4z_bb-2.0E0*2*I_NAI_G3xz_Fx2z_b-2.0E0*3*I_NAI_G3xz_Fx2z_b+2*1*I_NAI_G3xz_Px;
  abcd[3078] = 4.0E0*I_NAI_G2x2y_Hx4z_bb-2.0E0*2*I_NAI_G2x2y_Fx2z_b-2.0E0*3*I_NAI_G2x2y_Fx2z_b+2*1*I_NAI_G2x2y_Px;
  abcd[3079] = 4.0E0*I_NAI_G2xyz_Hx4z_bb-2.0E0*2*I_NAI_G2xyz_Fx2z_b-2.0E0*3*I_NAI_G2xyz_Fx2z_b+2*1*I_NAI_G2xyz_Px;
  abcd[3080] = 4.0E0*I_NAI_G2x2z_Hx4z_bb-2.0E0*2*I_NAI_G2x2z_Fx2z_b-2.0E0*3*I_NAI_G2x2z_Fx2z_b+2*1*I_NAI_G2x2z_Px;
  abcd[3081] = 4.0E0*I_NAI_Gx3y_Hx4z_bb-2.0E0*2*I_NAI_Gx3y_Fx2z_b-2.0E0*3*I_NAI_Gx3y_Fx2z_b+2*1*I_NAI_Gx3y_Px;
  abcd[3082] = 4.0E0*I_NAI_Gx2yz_Hx4z_bb-2.0E0*2*I_NAI_Gx2yz_Fx2z_b-2.0E0*3*I_NAI_Gx2yz_Fx2z_b+2*1*I_NAI_Gx2yz_Px;
  abcd[3083] = 4.0E0*I_NAI_Gxy2z_Hx4z_bb-2.0E0*2*I_NAI_Gxy2z_Fx2z_b-2.0E0*3*I_NAI_Gxy2z_Fx2z_b+2*1*I_NAI_Gxy2z_Px;
  abcd[3084] = 4.0E0*I_NAI_Gx3z_Hx4z_bb-2.0E0*2*I_NAI_Gx3z_Fx2z_b-2.0E0*3*I_NAI_Gx3z_Fx2z_b+2*1*I_NAI_Gx3z_Px;
  abcd[3085] = 4.0E0*I_NAI_G4y_Hx4z_bb-2.0E0*2*I_NAI_G4y_Fx2z_b-2.0E0*3*I_NAI_G4y_Fx2z_b+2*1*I_NAI_G4y_Px;
  abcd[3086] = 4.0E0*I_NAI_G3yz_Hx4z_bb-2.0E0*2*I_NAI_G3yz_Fx2z_b-2.0E0*3*I_NAI_G3yz_Fx2z_b+2*1*I_NAI_G3yz_Px;
  abcd[3087] = 4.0E0*I_NAI_G2y2z_Hx4z_bb-2.0E0*2*I_NAI_G2y2z_Fx2z_b-2.0E0*3*I_NAI_G2y2z_Fx2z_b+2*1*I_NAI_G2y2z_Px;
  abcd[3088] = 4.0E0*I_NAI_Gy3z_Hx4z_bb-2.0E0*2*I_NAI_Gy3z_Fx2z_b-2.0E0*3*I_NAI_Gy3z_Fx2z_b+2*1*I_NAI_Gy3z_Px;
  abcd[3089] = 4.0E0*I_NAI_G4z_Hx4z_bb-2.0E0*2*I_NAI_G4z_Fx2z_b-2.0E0*3*I_NAI_G4z_Fx2z_b+2*1*I_NAI_G4z_Px;
  abcd[3090] = 4.0E0*I_NAI_G4x_H3y2z_bb-2.0E0*1*I_NAI_G4x_F3y_b;
  abcd[3091] = 4.0E0*I_NAI_G3xy_H3y2z_bb-2.0E0*1*I_NAI_G3xy_F3y_b;
  abcd[3092] = 4.0E0*I_NAI_G3xz_H3y2z_bb-2.0E0*1*I_NAI_G3xz_F3y_b;
  abcd[3093] = 4.0E0*I_NAI_G2x2y_H3y2z_bb-2.0E0*1*I_NAI_G2x2y_F3y_b;
  abcd[3094] = 4.0E0*I_NAI_G2xyz_H3y2z_bb-2.0E0*1*I_NAI_G2xyz_F3y_b;
  abcd[3095] = 4.0E0*I_NAI_G2x2z_H3y2z_bb-2.0E0*1*I_NAI_G2x2z_F3y_b;
  abcd[3096] = 4.0E0*I_NAI_Gx3y_H3y2z_bb-2.0E0*1*I_NAI_Gx3y_F3y_b;
  abcd[3097] = 4.0E0*I_NAI_Gx2yz_H3y2z_bb-2.0E0*1*I_NAI_Gx2yz_F3y_b;
  abcd[3098] = 4.0E0*I_NAI_Gxy2z_H3y2z_bb-2.0E0*1*I_NAI_Gxy2z_F3y_b;
  abcd[3099] = 4.0E0*I_NAI_Gx3z_H3y2z_bb-2.0E0*1*I_NAI_Gx3z_F3y_b;
  abcd[3100] = 4.0E0*I_NAI_G4y_H3y2z_bb-2.0E0*1*I_NAI_G4y_F3y_b;
  abcd[3101] = 4.0E0*I_NAI_G3yz_H3y2z_bb-2.0E0*1*I_NAI_G3yz_F3y_b;
  abcd[3102] = 4.0E0*I_NAI_G2y2z_H3y2z_bb-2.0E0*1*I_NAI_G2y2z_F3y_b;
  abcd[3103] = 4.0E0*I_NAI_Gy3z_H3y2z_bb-2.0E0*1*I_NAI_Gy3z_F3y_b;
  abcd[3104] = 4.0E0*I_NAI_G4z_H3y2z_bb-2.0E0*1*I_NAI_G4z_F3y_b;
  abcd[3105] = 4.0E0*I_NAI_G4x_H2y3z_bb-2.0E0*1*I_NAI_G4x_F2yz_b-2.0E0*2*I_NAI_G4x_F2yz_b;
  abcd[3106] = 4.0E0*I_NAI_G3xy_H2y3z_bb-2.0E0*1*I_NAI_G3xy_F2yz_b-2.0E0*2*I_NAI_G3xy_F2yz_b;
  abcd[3107] = 4.0E0*I_NAI_G3xz_H2y3z_bb-2.0E0*1*I_NAI_G3xz_F2yz_b-2.0E0*2*I_NAI_G3xz_F2yz_b;
  abcd[3108] = 4.0E0*I_NAI_G2x2y_H2y3z_bb-2.0E0*1*I_NAI_G2x2y_F2yz_b-2.0E0*2*I_NAI_G2x2y_F2yz_b;
  abcd[3109] = 4.0E0*I_NAI_G2xyz_H2y3z_bb-2.0E0*1*I_NAI_G2xyz_F2yz_b-2.0E0*2*I_NAI_G2xyz_F2yz_b;
  abcd[3110] = 4.0E0*I_NAI_G2x2z_H2y3z_bb-2.0E0*1*I_NAI_G2x2z_F2yz_b-2.0E0*2*I_NAI_G2x2z_F2yz_b;
  abcd[3111] = 4.0E0*I_NAI_Gx3y_H2y3z_bb-2.0E0*1*I_NAI_Gx3y_F2yz_b-2.0E0*2*I_NAI_Gx3y_F2yz_b;
  abcd[3112] = 4.0E0*I_NAI_Gx2yz_H2y3z_bb-2.0E0*1*I_NAI_Gx2yz_F2yz_b-2.0E0*2*I_NAI_Gx2yz_F2yz_b;
  abcd[3113] = 4.0E0*I_NAI_Gxy2z_H2y3z_bb-2.0E0*1*I_NAI_Gxy2z_F2yz_b-2.0E0*2*I_NAI_Gxy2z_F2yz_b;
  abcd[3114] = 4.0E0*I_NAI_Gx3z_H2y3z_bb-2.0E0*1*I_NAI_Gx3z_F2yz_b-2.0E0*2*I_NAI_Gx3z_F2yz_b;
  abcd[3115] = 4.0E0*I_NAI_G4y_H2y3z_bb-2.0E0*1*I_NAI_G4y_F2yz_b-2.0E0*2*I_NAI_G4y_F2yz_b;
  abcd[3116] = 4.0E0*I_NAI_G3yz_H2y3z_bb-2.0E0*1*I_NAI_G3yz_F2yz_b-2.0E0*2*I_NAI_G3yz_F2yz_b;
  abcd[3117] = 4.0E0*I_NAI_G2y2z_H2y3z_bb-2.0E0*1*I_NAI_G2y2z_F2yz_b-2.0E0*2*I_NAI_G2y2z_F2yz_b;
  abcd[3118] = 4.0E0*I_NAI_Gy3z_H2y3z_bb-2.0E0*1*I_NAI_Gy3z_F2yz_b-2.0E0*2*I_NAI_Gy3z_F2yz_b;
  abcd[3119] = 4.0E0*I_NAI_G4z_H2y3z_bb-2.0E0*1*I_NAI_G4z_F2yz_b-2.0E0*2*I_NAI_G4z_F2yz_b;
  abcd[3120] = 4.0E0*I_NAI_G4x_Hy4z_bb-2.0E0*2*I_NAI_G4x_Fy2z_b-2.0E0*3*I_NAI_G4x_Fy2z_b+2*1*I_NAI_G4x_Py;
  abcd[3121] = 4.0E0*I_NAI_G3xy_Hy4z_bb-2.0E0*2*I_NAI_G3xy_Fy2z_b-2.0E0*3*I_NAI_G3xy_Fy2z_b+2*1*I_NAI_G3xy_Py;
  abcd[3122] = 4.0E0*I_NAI_G3xz_Hy4z_bb-2.0E0*2*I_NAI_G3xz_Fy2z_b-2.0E0*3*I_NAI_G3xz_Fy2z_b+2*1*I_NAI_G3xz_Py;
  abcd[3123] = 4.0E0*I_NAI_G2x2y_Hy4z_bb-2.0E0*2*I_NAI_G2x2y_Fy2z_b-2.0E0*3*I_NAI_G2x2y_Fy2z_b+2*1*I_NAI_G2x2y_Py;
  abcd[3124] = 4.0E0*I_NAI_G2xyz_Hy4z_bb-2.0E0*2*I_NAI_G2xyz_Fy2z_b-2.0E0*3*I_NAI_G2xyz_Fy2z_b+2*1*I_NAI_G2xyz_Py;
  abcd[3125] = 4.0E0*I_NAI_G2x2z_Hy4z_bb-2.0E0*2*I_NAI_G2x2z_Fy2z_b-2.0E0*3*I_NAI_G2x2z_Fy2z_b+2*1*I_NAI_G2x2z_Py;
  abcd[3126] = 4.0E0*I_NAI_Gx3y_Hy4z_bb-2.0E0*2*I_NAI_Gx3y_Fy2z_b-2.0E0*3*I_NAI_Gx3y_Fy2z_b+2*1*I_NAI_Gx3y_Py;
  abcd[3127] = 4.0E0*I_NAI_Gx2yz_Hy4z_bb-2.0E0*2*I_NAI_Gx2yz_Fy2z_b-2.0E0*3*I_NAI_Gx2yz_Fy2z_b+2*1*I_NAI_Gx2yz_Py;
  abcd[3128] = 4.0E0*I_NAI_Gxy2z_Hy4z_bb-2.0E0*2*I_NAI_Gxy2z_Fy2z_b-2.0E0*3*I_NAI_Gxy2z_Fy2z_b+2*1*I_NAI_Gxy2z_Py;
  abcd[3129] = 4.0E0*I_NAI_Gx3z_Hy4z_bb-2.0E0*2*I_NAI_Gx3z_Fy2z_b-2.0E0*3*I_NAI_Gx3z_Fy2z_b+2*1*I_NAI_Gx3z_Py;
  abcd[3130] = 4.0E0*I_NAI_G4y_Hy4z_bb-2.0E0*2*I_NAI_G4y_Fy2z_b-2.0E0*3*I_NAI_G4y_Fy2z_b+2*1*I_NAI_G4y_Py;
  abcd[3131] = 4.0E0*I_NAI_G3yz_Hy4z_bb-2.0E0*2*I_NAI_G3yz_Fy2z_b-2.0E0*3*I_NAI_G3yz_Fy2z_b+2*1*I_NAI_G3yz_Py;
  abcd[3132] = 4.0E0*I_NAI_G2y2z_Hy4z_bb-2.0E0*2*I_NAI_G2y2z_Fy2z_b-2.0E0*3*I_NAI_G2y2z_Fy2z_b+2*1*I_NAI_G2y2z_Py;
  abcd[3133] = 4.0E0*I_NAI_Gy3z_Hy4z_bb-2.0E0*2*I_NAI_Gy3z_Fy2z_b-2.0E0*3*I_NAI_Gy3z_Fy2z_b+2*1*I_NAI_Gy3z_Py;
  abcd[3134] = 4.0E0*I_NAI_G4z_Hy4z_bb-2.0E0*2*I_NAI_G4z_Fy2z_b-2.0E0*3*I_NAI_G4z_Fy2z_b+2*1*I_NAI_G4z_Py;
  abcd[3135] = 4.0E0*I_NAI_G4x_H5z_bb-2.0E0*3*I_NAI_G4x_F3z_b-2.0E0*4*I_NAI_G4x_F3z_b+3*2*I_NAI_G4x_Pz;
  abcd[3136] = 4.0E0*I_NAI_G3xy_H5z_bb-2.0E0*3*I_NAI_G3xy_F3z_b-2.0E0*4*I_NAI_G3xy_F3z_b+3*2*I_NAI_G3xy_Pz;
  abcd[3137] = 4.0E0*I_NAI_G3xz_H5z_bb-2.0E0*3*I_NAI_G3xz_F3z_b-2.0E0*4*I_NAI_G3xz_F3z_b+3*2*I_NAI_G3xz_Pz;
  abcd[3138] = 4.0E0*I_NAI_G2x2y_H5z_bb-2.0E0*3*I_NAI_G2x2y_F3z_b-2.0E0*4*I_NAI_G2x2y_F3z_b+3*2*I_NAI_G2x2y_Pz;
  abcd[3139] = 4.0E0*I_NAI_G2xyz_H5z_bb-2.0E0*3*I_NAI_G2xyz_F3z_b-2.0E0*4*I_NAI_G2xyz_F3z_b+3*2*I_NAI_G2xyz_Pz;
  abcd[3140] = 4.0E0*I_NAI_G2x2z_H5z_bb-2.0E0*3*I_NAI_G2x2z_F3z_b-2.0E0*4*I_NAI_G2x2z_F3z_b+3*2*I_NAI_G2x2z_Pz;
  abcd[3141] = 4.0E0*I_NAI_Gx3y_H5z_bb-2.0E0*3*I_NAI_Gx3y_F3z_b-2.0E0*4*I_NAI_Gx3y_F3z_b+3*2*I_NAI_Gx3y_Pz;
  abcd[3142] = 4.0E0*I_NAI_Gx2yz_H5z_bb-2.0E0*3*I_NAI_Gx2yz_F3z_b-2.0E0*4*I_NAI_Gx2yz_F3z_b+3*2*I_NAI_Gx2yz_Pz;
  abcd[3143] = 4.0E0*I_NAI_Gxy2z_H5z_bb-2.0E0*3*I_NAI_Gxy2z_F3z_b-2.0E0*4*I_NAI_Gxy2z_F3z_b+3*2*I_NAI_Gxy2z_Pz;
  abcd[3144] = 4.0E0*I_NAI_Gx3z_H5z_bb-2.0E0*3*I_NAI_Gx3z_F3z_b-2.0E0*4*I_NAI_Gx3z_F3z_b+3*2*I_NAI_Gx3z_Pz;
  abcd[3145] = 4.0E0*I_NAI_G4y_H5z_bb-2.0E0*3*I_NAI_G4y_F3z_b-2.0E0*4*I_NAI_G4y_F3z_b+3*2*I_NAI_G4y_Pz;
  abcd[3146] = 4.0E0*I_NAI_G3yz_H5z_bb-2.0E0*3*I_NAI_G3yz_F3z_b-2.0E0*4*I_NAI_G3yz_F3z_b+3*2*I_NAI_G3yz_Pz;
  abcd[3147] = 4.0E0*I_NAI_G2y2z_H5z_bb-2.0E0*3*I_NAI_G2y2z_F3z_b-2.0E0*4*I_NAI_G2y2z_F3z_b+3*2*I_NAI_G2y2z_Pz;
  abcd[3148] = 4.0E0*I_NAI_Gy3z_H5z_bb-2.0E0*3*I_NAI_Gy3z_F3z_b-2.0E0*4*I_NAI_Gy3z_F3z_b+3*2*I_NAI_Gy3z_Pz;
  abcd[3149] = 4.0E0*I_NAI_G4z_H5z_bb-2.0E0*3*I_NAI_G4z_F3z_b-2.0E0*4*I_NAI_G4z_F3z_b+3*2*I_NAI_G4z_Pz;
}
