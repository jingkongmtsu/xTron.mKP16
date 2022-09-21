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

void hgp_os_nai_h_p_d2(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dxy_S_vrr = PAY*I_NAI_Px_S_vrr-PNY*I_NAI_Px_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
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
   * totally 0 integrals are omitted 
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
  Double I_NAI_G4x_Dxz_b = I_NAI_H4xz_Px_b+ABZ*I_NAI_G4x_Px_b;
  Double I_NAI_G3xy_Dxz_b = I_NAI_H3xyz_Px_b+ABZ*I_NAI_G3xy_Px_b;
  Double I_NAI_G3xz_Dxz_b = I_NAI_H3x2z_Px_b+ABZ*I_NAI_G3xz_Px_b;
  Double I_NAI_G2x2y_Dxz_b = I_NAI_H2x2yz_Px_b+ABZ*I_NAI_G2x2y_Px_b;
  Double I_NAI_G2xyz_Dxz_b = I_NAI_H2xy2z_Px_b+ABZ*I_NAI_G2xyz_Px_b;
  Double I_NAI_G2x2z_Dxz_b = I_NAI_H2x3z_Px_b+ABZ*I_NAI_G2x2z_Px_b;
  Double I_NAI_Gx3y_Dxz_b = I_NAI_Hx3yz_Px_b+ABZ*I_NAI_Gx3y_Px_b;
  Double I_NAI_Gx2yz_Dxz_b = I_NAI_Hx2y2z_Px_b+ABZ*I_NAI_Gx2yz_Px_b;
  Double I_NAI_Gxy2z_Dxz_b = I_NAI_Hxy3z_Px_b+ABZ*I_NAI_Gxy2z_Px_b;
  Double I_NAI_Gx3z_Dxz_b = I_NAI_Hx4z_Px_b+ABZ*I_NAI_Gx3z_Px_b;
  Double I_NAI_G4y_Dxz_b = I_NAI_H4yz_Px_b+ABZ*I_NAI_G4y_Px_b;
  Double I_NAI_G3yz_Dxz_b = I_NAI_H3y2z_Px_b+ABZ*I_NAI_G3yz_Px_b;
  Double I_NAI_G2y2z_Dxz_b = I_NAI_H2y3z_Px_b+ABZ*I_NAI_G2y2z_Px_b;
  Double I_NAI_Gy3z_Dxz_b = I_NAI_Hy4z_Px_b+ABZ*I_NAI_Gy3z_Px_b;
  Double I_NAI_G4z_Dxz_b = I_NAI_H5z_Px_b+ABZ*I_NAI_G4z_Px_b;
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
  Double I_NAI_G4x_Dyz_b = I_NAI_H4xz_Py_b+ABZ*I_NAI_G4x_Py_b;
  Double I_NAI_G3xy_Dyz_b = I_NAI_H3xyz_Py_b+ABZ*I_NAI_G3xy_Py_b;
  Double I_NAI_G3xz_Dyz_b = I_NAI_H3x2z_Py_b+ABZ*I_NAI_G3xz_Py_b;
  Double I_NAI_G2x2y_Dyz_b = I_NAI_H2x2yz_Py_b+ABZ*I_NAI_G2x2y_Py_b;
  Double I_NAI_G2xyz_Dyz_b = I_NAI_H2xy2z_Py_b+ABZ*I_NAI_G2xyz_Py_b;
  Double I_NAI_G2x2z_Dyz_b = I_NAI_H2x3z_Py_b+ABZ*I_NAI_G2x2z_Py_b;
  Double I_NAI_Gx3y_Dyz_b = I_NAI_Hx3yz_Py_b+ABZ*I_NAI_Gx3y_Py_b;
  Double I_NAI_Gx2yz_Dyz_b = I_NAI_Hx2y2z_Py_b+ABZ*I_NAI_Gx2yz_Py_b;
  Double I_NAI_Gxy2z_Dyz_b = I_NAI_Hxy3z_Py_b+ABZ*I_NAI_Gxy2z_Py_b;
  Double I_NAI_Gx3z_Dyz_b = I_NAI_Hx4z_Py_b+ABZ*I_NAI_Gx3z_Py_b;
  Double I_NAI_G4y_Dyz_b = I_NAI_H4yz_Py_b+ABZ*I_NAI_G4y_Py_b;
  Double I_NAI_G3yz_Dyz_b = I_NAI_H3y2z_Py_b+ABZ*I_NAI_G3yz_Py_b;
  Double I_NAI_G2y2z_Dyz_b = I_NAI_H2y3z_Py_b+ABZ*I_NAI_G2y2z_Py_b;
  Double I_NAI_Gy3z_Dyz_b = I_NAI_Hy4z_Py_b+ABZ*I_NAI_Gy3z_Py_b;
  Double I_NAI_G4z_Dyz_b = I_NAI_H5z_Py_b+ABZ*I_NAI_G4z_Py_b;
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
   * shell quartet name: SQ_NAI_K_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 9 integrals are omitted 
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
  Double I_NAI_K6xz_Pz_ab = I_NAI_L6x2z_S_ab+ABZ*I_NAI_K6xz_S_ab;
  Double I_NAI_K5xyz_Pz_ab = I_NAI_L5xy2z_S_ab+ABZ*I_NAI_K5xyz_S_ab;
  Double I_NAI_K5x2z_Pz_ab = I_NAI_L5x3z_S_ab+ABZ*I_NAI_K5x2z_S_ab;
  Double I_NAI_K4x2yz_Pz_ab = I_NAI_L4x2y2z_S_ab+ABZ*I_NAI_K4x2yz_S_ab;
  Double I_NAI_K4xy2z_Pz_ab = I_NAI_L4xy3z_S_ab+ABZ*I_NAI_K4xy2z_S_ab;
  Double I_NAI_K4x3z_Pz_ab = I_NAI_L4x4z_S_ab+ABZ*I_NAI_K4x3z_S_ab;
  Double I_NAI_K3x3yz_Pz_ab = I_NAI_L3x3y2z_S_ab+ABZ*I_NAI_K3x3yz_S_ab;
  Double I_NAI_K3x2y2z_Pz_ab = I_NAI_L3x2y3z_S_ab+ABZ*I_NAI_K3x2y2z_S_ab;
  Double I_NAI_K3xy3z_Pz_ab = I_NAI_L3xy4z_S_ab+ABZ*I_NAI_K3xy3z_S_ab;
  Double I_NAI_K3x4z_Pz_ab = I_NAI_L3x5z_S_ab+ABZ*I_NAI_K3x4z_S_ab;
  Double I_NAI_K2x4yz_Pz_ab = I_NAI_L2x4y2z_S_ab+ABZ*I_NAI_K2x4yz_S_ab;
  Double I_NAI_K2x3y2z_Pz_ab = I_NAI_L2x3y3z_S_ab+ABZ*I_NAI_K2x3y2z_S_ab;
  Double I_NAI_K2x2y3z_Pz_ab = I_NAI_L2x2y4z_S_ab+ABZ*I_NAI_K2x2y3z_S_ab;
  Double I_NAI_K2xy4z_Pz_ab = I_NAI_L2xy5z_S_ab+ABZ*I_NAI_K2xy4z_S_ab;
  Double I_NAI_K2x5z_Pz_ab = I_NAI_L2x6z_S_ab+ABZ*I_NAI_K2x5z_S_ab;
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
   * totally 0 integrals are omitted 
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
  Double I_NAI_I6x_Dxy_ab = I_NAI_K6xy_Px_ab+ABY*I_NAI_I6x_Px_ab;
  Double I_NAI_I5xy_Dxy_ab = I_NAI_K5x2y_Px_ab+ABY*I_NAI_I5xy_Px_ab;
  Double I_NAI_I5xz_Dxy_ab = I_NAI_K5xyz_Px_ab+ABY*I_NAI_I5xz_Px_ab;
  Double I_NAI_I4x2y_Dxy_ab = I_NAI_K4x3y_Px_ab+ABY*I_NAI_I4x2y_Px_ab;
  Double I_NAI_I4xyz_Dxy_ab = I_NAI_K4x2yz_Px_ab+ABY*I_NAI_I4xyz_Px_ab;
  Double I_NAI_I4x2z_Dxy_ab = I_NAI_K4xy2z_Px_ab+ABY*I_NAI_I4x2z_Px_ab;
  Double I_NAI_I3x3y_Dxy_ab = I_NAI_K3x4y_Px_ab+ABY*I_NAI_I3x3y_Px_ab;
  Double I_NAI_I3x2yz_Dxy_ab = I_NAI_K3x3yz_Px_ab+ABY*I_NAI_I3x2yz_Px_ab;
  Double I_NAI_I3xy2z_Dxy_ab = I_NAI_K3x2y2z_Px_ab+ABY*I_NAI_I3xy2z_Px_ab;
  Double I_NAI_I3x3z_Dxy_ab = I_NAI_K3xy3z_Px_ab+ABY*I_NAI_I3x3z_Px_ab;
  Double I_NAI_I2x4y_Dxy_ab = I_NAI_K2x5y_Px_ab+ABY*I_NAI_I2x4y_Px_ab;
  Double I_NAI_I2x3yz_Dxy_ab = I_NAI_K2x4yz_Px_ab+ABY*I_NAI_I2x3yz_Px_ab;
  Double I_NAI_I2x2y2z_Dxy_ab = I_NAI_K2x3y2z_Px_ab+ABY*I_NAI_I2x2y2z_Px_ab;
  Double I_NAI_I2xy3z_Dxy_ab = I_NAI_K2x2y3z_Px_ab+ABY*I_NAI_I2xy3z_Px_ab;
  Double I_NAI_I2x4z_Dxy_ab = I_NAI_K2xy4z_Px_ab+ABY*I_NAI_I2x4z_Px_ab;
  Double I_NAI_Ix5y_Dxy_ab = I_NAI_Kx6y_Px_ab+ABY*I_NAI_Ix5y_Px_ab;
  Double I_NAI_Ix4yz_Dxy_ab = I_NAI_Kx5yz_Px_ab+ABY*I_NAI_Ix4yz_Px_ab;
  Double I_NAI_Ix3y2z_Dxy_ab = I_NAI_Kx4y2z_Px_ab+ABY*I_NAI_Ix3y2z_Px_ab;
  Double I_NAI_Ix2y3z_Dxy_ab = I_NAI_Kx3y3z_Px_ab+ABY*I_NAI_Ix2y3z_Px_ab;
  Double I_NAI_Ixy4z_Dxy_ab = I_NAI_Kx2y4z_Px_ab+ABY*I_NAI_Ixy4z_Px_ab;
  Double I_NAI_Ix5z_Dxy_ab = I_NAI_Kxy5z_Px_ab+ABY*I_NAI_Ix5z_Px_ab;
  Double I_NAI_I6y_Dxy_ab = I_NAI_K7y_Px_ab+ABY*I_NAI_I6y_Px_ab;
  Double I_NAI_I5yz_Dxy_ab = I_NAI_K6yz_Px_ab+ABY*I_NAI_I5yz_Px_ab;
  Double I_NAI_I4y2z_Dxy_ab = I_NAI_K5y2z_Px_ab+ABY*I_NAI_I4y2z_Px_ab;
  Double I_NAI_I3y3z_Dxy_ab = I_NAI_K4y3z_Px_ab+ABY*I_NAI_I3y3z_Px_ab;
  Double I_NAI_I2y4z_Dxy_ab = I_NAI_K3y4z_Px_ab+ABY*I_NAI_I2y4z_Px_ab;
  Double I_NAI_Iy5z_Dxy_ab = I_NAI_K2y5z_Px_ab+ABY*I_NAI_Iy5z_Px_ab;
  Double I_NAI_I6z_Dxy_ab = I_NAI_Ky6z_Px_ab+ABY*I_NAI_I6z_Px_ab;
  Double I_NAI_I6x_Dxz_ab = I_NAI_K6xz_Px_ab+ABZ*I_NAI_I6x_Px_ab;
  Double I_NAI_I5xy_Dxz_ab = I_NAI_K5xyz_Px_ab+ABZ*I_NAI_I5xy_Px_ab;
  Double I_NAI_I5xz_Dxz_ab = I_NAI_K5x2z_Px_ab+ABZ*I_NAI_I5xz_Px_ab;
  Double I_NAI_I4x2y_Dxz_ab = I_NAI_K4x2yz_Px_ab+ABZ*I_NAI_I4x2y_Px_ab;
  Double I_NAI_I4xyz_Dxz_ab = I_NAI_K4xy2z_Px_ab+ABZ*I_NAI_I4xyz_Px_ab;
  Double I_NAI_I4x2z_Dxz_ab = I_NAI_K4x3z_Px_ab+ABZ*I_NAI_I4x2z_Px_ab;
  Double I_NAI_I3x3y_Dxz_ab = I_NAI_K3x3yz_Px_ab+ABZ*I_NAI_I3x3y_Px_ab;
  Double I_NAI_I3x2yz_Dxz_ab = I_NAI_K3x2y2z_Px_ab+ABZ*I_NAI_I3x2yz_Px_ab;
  Double I_NAI_I3xy2z_Dxz_ab = I_NAI_K3xy3z_Px_ab+ABZ*I_NAI_I3xy2z_Px_ab;
  Double I_NAI_I3x3z_Dxz_ab = I_NAI_K3x4z_Px_ab+ABZ*I_NAI_I3x3z_Px_ab;
  Double I_NAI_I2x4y_Dxz_ab = I_NAI_K2x4yz_Px_ab+ABZ*I_NAI_I2x4y_Px_ab;
  Double I_NAI_I2x3yz_Dxz_ab = I_NAI_K2x3y2z_Px_ab+ABZ*I_NAI_I2x3yz_Px_ab;
  Double I_NAI_I2x2y2z_Dxz_ab = I_NAI_K2x2y3z_Px_ab+ABZ*I_NAI_I2x2y2z_Px_ab;
  Double I_NAI_I2xy3z_Dxz_ab = I_NAI_K2xy4z_Px_ab+ABZ*I_NAI_I2xy3z_Px_ab;
  Double I_NAI_I2x4z_Dxz_ab = I_NAI_K2x5z_Px_ab+ABZ*I_NAI_I2x4z_Px_ab;
  Double I_NAI_Ix5y_Dxz_ab = I_NAI_Kx5yz_Px_ab+ABZ*I_NAI_Ix5y_Px_ab;
  Double I_NAI_Ix4yz_Dxz_ab = I_NAI_Kx4y2z_Px_ab+ABZ*I_NAI_Ix4yz_Px_ab;
  Double I_NAI_Ix3y2z_Dxz_ab = I_NAI_Kx3y3z_Px_ab+ABZ*I_NAI_Ix3y2z_Px_ab;
  Double I_NAI_Ix2y3z_Dxz_ab = I_NAI_Kx2y4z_Px_ab+ABZ*I_NAI_Ix2y3z_Px_ab;
  Double I_NAI_Ixy4z_Dxz_ab = I_NAI_Kxy5z_Px_ab+ABZ*I_NAI_Ixy4z_Px_ab;
  Double I_NAI_Ix5z_Dxz_ab = I_NAI_Kx6z_Px_ab+ABZ*I_NAI_Ix5z_Px_ab;
  Double I_NAI_I6y_Dxz_ab = I_NAI_K6yz_Px_ab+ABZ*I_NAI_I6y_Px_ab;
  Double I_NAI_I5yz_Dxz_ab = I_NAI_K5y2z_Px_ab+ABZ*I_NAI_I5yz_Px_ab;
  Double I_NAI_I4y2z_Dxz_ab = I_NAI_K4y3z_Px_ab+ABZ*I_NAI_I4y2z_Px_ab;
  Double I_NAI_I3y3z_Dxz_ab = I_NAI_K3y4z_Px_ab+ABZ*I_NAI_I3y3z_Px_ab;
  Double I_NAI_I2y4z_Dxz_ab = I_NAI_K2y5z_Px_ab+ABZ*I_NAI_I2y4z_Px_ab;
  Double I_NAI_Iy5z_Dxz_ab = I_NAI_Ky6z_Px_ab+ABZ*I_NAI_Iy5z_Px_ab;
  Double I_NAI_I6z_Dxz_ab = I_NAI_K7z_Px_ab+ABZ*I_NAI_I6z_Px_ab;
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
  Double I_NAI_I6x_Dyz_ab = I_NAI_K6xz_Py_ab+ABZ*I_NAI_I6x_Py_ab;
  Double I_NAI_I5xy_Dyz_ab = I_NAI_K5xyz_Py_ab+ABZ*I_NAI_I5xy_Py_ab;
  Double I_NAI_I5xz_Dyz_ab = I_NAI_K5x2z_Py_ab+ABZ*I_NAI_I5xz_Py_ab;
  Double I_NAI_I4x2y_Dyz_ab = I_NAI_K4x2yz_Py_ab+ABZ*I_NAI_I4x2y_Py_ab;
  Double I_NAI_I4xyz_Dyz_ab = I_NAI_K4xy2z_Py_ab+ABZ*I_NAI_I4xyz_Py_ab;
  Double I_NAI_I4x2z_Dyz_ab = I_NAI_K4x3z_Py_ab+ABZ*I_NAI_I4x2z_Py_ab;
  Double I_NAI_I3x3y_Dyz_ab = I_NAI_K3x3yz_Py_ab+ABZ*I_NAI_I3x3y_Py_ab;
  Double I_NAI_I3x2yz_Dyz_ab = I_NAI_K3x2y2z_Py_ab+ABZ*I_NAI_I3x2yz_Py_ab;
  Double I_NAI_I3xy2z_Dyz_ab = I_NAI_K3xy3z_Py_ab+ABZ*I_NAI_I3xy2z_Py_ab;
  Double I_NAI_I3x3z_Dyz_ab = I_NAI_K3x4z_Py_ab+ABZ*I_NAI_I3x3z_Py_ab;
  Double I_NAI_I2x4y_Dyz_ab = I_NAI_K2x4yz_Py_ab+ABZ*I_NAI_I2x4y_Py_ab;
  Double I_NAI_I2x3yz_Dyz_ab = I_NAI_K2x3y2z_Py_ab+ABZ*I_NAI_I2x3yz_Py_ab;
  Double I_NAI_I2x2y2z_Dyz_ab = I_NAI_K2x2y3z_Py_ab+ABZ*I_NAI_I2x2y2z_Py_ab;
  Double I_NAI_I2xy3z_Dyz_ab = I_NAI_K2xy4z_Py_ab+ABZ*I_NAI_I2xy3z_Py_ab;
  Double I_NAI_I2x4z_Dyz_ab = I_NAI_K2x5z_Py_ab+ABZ*I_NAI_I2x4z_Py_ab;
  Double I_NAI_Ix5y_Dyz_ab = I_NAI_Kx5yz_Py_ab+ABZ*I_NAI_Ix5y_Py_ab;
  Double I_NAI_Ix4yz_Dyz_ab = I_NAI_Kx4y2z_Py_ab+ABZ*I_NAI_Ix4yz_Py_ab;
  Double I_NAI_Ix3y2z_Dyz_ab = I_NAI_Kx3y3z_Py_ab+ABZ*I_NAI_Ix3y2z_Py_ab;
  Double I_NAI_Ix2y3z_Dyz_ab = I_NAI_Kx2y4z_Py_ab+ABZ*I_NAI_Ix2y3z_Py_ab;
  Double I_NAI_Ixy4z_Dyz_ab = I_NAI_Kxy5z_Py_ab+ABZ*I_NAI_Ixy4z_Py_ab;
  Double I_NAI_Ix5z_Dyz_ab = I_NAI_Kx6z_Py_ab+ABZ*I_NAI_Ix5z_Py_ab;
  Double I_NAI_I6y_Dyz_ab = I_NAI_K6yz_Py_ab+ABZ*I_NAI_I6y_Py_ab;
  Double I_NAI_I5yz_Dyz_ab = I_NAI_K5y2z_Py_ab+ABZ*I_NAI_I5yz_Py_ab;
  Double I_NAI_I4y2z_Dyz_ab = I_NAI_K4y3z_Py_ab+ABZ*I_NAI_I4y2z_Py_ab;
  Double I_NAI_I3y3z_Dyz_ab = I_NAI_K3y4z_Py_ab+ABZ*I_NAI_I3y3z_Py_ab;
  Double I_NAI_I2y4z_Dyz_ab = I_NAI_K2y5z_Py_ab+ABZ*I_NAI_I2y4z_Py_ab;
  Double I_NAI_Iy5z_Dyz_ab = I_NAI_Ky6z_Py_ab+ABZ*I_NAI_Iy5z_Py_ab;
  Double I_NAI_I6z_Dyz_ab = I_NAI_K7z_Py_ab+ABZ*I_NAI_I6z_Py_ab;
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
   * totally 42 integrals are omitted 
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
  Double I_NAI_H5x_Dxy_bb = I_NAI_I5xy_Px_bb+ABY*I_NAI_H5x_Px_bb;
  Double I_NAI_H4xy_Dxy_bb = I_NAI_I4x2y_Px_bb+ABY*I_NAI_H4xy_Px_bb;
  Double I_NAI_H4xz_Dxy_bb = I_NAI_I4xyz_Px_bb+ABY*I_NAI_H4xz_Px_bb;
  Double I_NAI_H3x2y_Dxy_bb = I_NAI_I3x3y_Px_bb+ABY*I_NAI_H3x2y_Px_bb;
  Double I_NAI_H3xyz_Dxy_bb = I_NAI_I3x2yz_Px_bb+ABY*I_NAI_H3xyz_Px_bb;
  Double I_NAI_H3x2z_Dxy_bb = I_NAI_I3xy2z_Px_bb+ABY*I_NAI_H3x2z_Px_bb;
  Double I_NAI_H2x3y_Dxy_bb = I_NAI_I2x4y_Px_bb+ABY*I_NAI_H2x3y_Px_bb;
  Double I_NAI_H2x2yz_Dxy_bb = I_NAI_I2x3yz_Px_bb+ABY*I_NAI_H2x2yz_Px_bb;
  Double I_NAI_H2xy2z_Dxy_bb = I_NAI_I2x2y2z_Px_bb+ABY*I_NAI_H2xy2z_Px_bb;
  Double I_NAI_H2x3z_Dxy_bb = I_NAI_I2xy3z_Px_bb+ABY*I_NAI_H2x3z_Px_bb;
  Double I_NAI_Hx4y_Dxy_bb = I_NAI_Ix5y_Px_bb+ABY*I_NAI_Hx4y_Px_bb;
  Double I_NAI_Hx3yz_Dxy_bb = I_NAI_Ix4yz_Px_bb+ABY*I_NAI_Hx3yz_Px_bb;
  Double I_NAI_Hx2y2z_Dxy_bb = I_NAI_Ix3y2z_Px_bb+ABY*I_NAI_Hx2y2z_Px_bb;
  Double I_NAI_Hxy3z_Dxy_bb = I_NAI_Ix2y3z_Px_bb+ABY*I_NAI_Hxy3z_Px_bb;
  Double I_NAI_Hx4z_Dxy_bb = I_NAI_Ixy4z_Px_bb+ABY*I_NAI_Hx4z_Px_bb;
  Double I_NAI_H5y_Dxy_bb = I_NAI_I6y_Px_bb+ABY*I_NAI_H5y_Px_bb;
  Double I_NAI_H4yz_Dxy_bb = I_NAI_I5yz_Px_bb+ABY*I_NAI_H4yz_Px_bb;
  Double I_NAI_H3y2z_Dxy_bb = I_NAI_I4y2z_Px_bb+ABY*I_NAI_H3y2z_Px_bb;
  Double I_NAI_H2y3z_Dxy_bb = I_NAI_I3y3z_Px_bb+ABY*I_NAI_H2y3z_Px_bb;
  Double I_NAI_Hy4z_Dxy_bb = I_NAI_I2y4z_Px_bb+ABY*I_NAI_Hy4z_Px_bb;
  Double I_NAI_H5z_Dxy_bb = I_NAI_Iy5z_Px_bb+ABY*I_NAI_H5z_Px_bb;
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
   * shell quartet name: SQ_NAI_K_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
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
  Double I_NAI_K4x3y_Py_bb = I_NAI_L4x4y_S_bb+ABY*I_NAI_K4x3y_S_bb;
  Double I_NAI_K4x2yz_Py_bb = I_NAI_L4x3yz_S_bb+ABY*I_NAI_K4x2yz_S_bb;
  Double I_NAI_K4xy2z_Py_bb = I_NAI_L4x2y2z_S_bb+ABY*I_NAI_K4xy2z_S_bb;
  Double I_NAI_K3x4y_Py_bb = I_NAI_L3x5y_S_bb+ABY*I_NAI_K3x4y_S_bb;
  Double I_NAI_K3x3yz_Py_bb = I_NAI_L3x4yz_S_bb+ABY*I_NAI_K3x3yz_S_bb;
  Double I_NAI_K3x2y2z_Py_bb = I_NAI_L3x3y2z_S_bb+ABY*I_NAI_K3x2y2z_S_bb;
  Double I_NAI_K3xy3z_Py_bb = I_NAI_L3x2y3z_S_bb+ABY*I_NAI_K3xy3z_S_bb;
  Double I_NAI_K2x5y_Py_bb = I_NAI_L2x6y_S_bb+ABY*I_NAI_K2x5y_S_bb;
  Double I_NAI_K2x4yz_Py_bb = I_NAI_L2x5yz_S_bb+ABY*I_NAI_K2x4yz_S_bb;
  Double I_NAI_K2x3y2z_Py_bb = I_NAI_L2x4y2z_S_bb+ABY*I_NAI_K2x3y2z_S_bb;
  Double I_NAI_K2x2y3z_Py_bb = I_NAI_L2x3y3z_S_bb+ABY*I_NAI_K2x2y3z_S_bb;
  Double I_NAI_K2xy4z_Py_bb = I_NAI_L2x2y4z_S_bb+ABY*I_NAI_K2xy4z_S_bb;
  Double I_NAI_Kx6y_Py_bb = I_NAI_Lx7y_S_bb+ABY*I_NAI_Kx6y_S_bb;
  Double I_NAI_Kx5yz_Py_bb = I_NAI_Lx6yz_S_bb+ABY*I_NAI_Kx5yz_S_bb;
  Double I_NAI_Kx4y2z_Py_bb = I_NAI_Lx5y2z_S_bb+ABY*I_NAI_Kx4y2z_S_bb;
  Double I_NAI_Kx3y3z_Py_bb = I_NAI_Lx4y3z_S_bb+ABY*I_NAI_Kx3y3z_S_bb;
  Double I_NAI_Kx2y4z_Py_bb = I_NAI_Lx3y4z_S_bb+ABY*I_NAI_Kx2y4z_S_bb;
  Double I_NAI_Kxy5z_Py_bb = I_NAI_Lx2y5z_S_bb+ABY*I_NAI_Kxy5z_S_bb;
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
   * totally 63 integrals are omitted 
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
  Double I_NAI_I5xz_Dxy_bb = I_NAI_K5xyz_Px_bb+ABY*I_NAI_I5xz_Px_bb;
  Double I_NAI_I4xyz_Dxy_bb = I_NAI_K4x2yz_Px_bb+ABY*I_NAI_I4xyz_Px_bb;
  Double I_NAI_I4x2z_Dxy_bb = I_NAI_K4xy2z_Px_bb+ABY*I_NAI_I4x2z_Px_bb;
  Double I_NAI_I3x2yz_Dxy_bb = I_NAI_K3x3yz_Px_bb+ABY*I_NAI_I3x2yz_Px_bb;
  Double I_NAI_I3xy2z_Dxy_bb = I_NAI_K3x2y2z_Px_bb+ABY*I_NAI_I3xy2z_Px_bb;
  Double I_NAI_I3x3z_Dxy_bb = I_NAI_K3xy3z_Px_bb+ABY*I_NAI_I3x3z_Px_bb;
  Double I_NAI_I2x3yz_Dxy_bb = I_NAI_K2x4yz_Px_bb+ABY*I_NAI_I2x3yz_Px_bb;
  Double I_NAI_I2x2y2z_Dxy_bb = I_NAI_K2x3y2z_Px_bb+ABY*I_NAI_I2x2y2z_Px_bb;
  Double I_NAI_I2xy3z_Dxy_bb = I_NAI_K2x2y3z_Px_bb+ABY*I_NAI_I2xy3z_Px_bb;
  Double I_NAI_I2x4z_Dxy_bb = I_NAI_K2xy4z_Px_bb+ABY*I_NAI_I2x4z_Px_bb;
  Double I_NAI_Ix4yz_Dxy_bb = I_NAI_Kx5yz_Px_bb+ABY*I_NAI_Ix4yz_Px_bb;
  Double I_NAI_Ix3y2z_Dxy_bb = I_NAI_Kx4y2z_Px_bb+ABY*I_NAI_Ix3y2z_Px_bb;
  Double I_NAI_Ix2y3z_Dxy_bb = I_NAI_Kx3y3z_Px_bb+ABY*I_NAI_Ix2y3z_Px_bb;
  Double I_NAI_Ixy4z_Dxy_bb = I_NAI_Kx2y4z_Px_bb+ABY*I_NAI_Ixy4z_Px_bb;
  Double I_NAI_Ix5z_Dxy_bb = I_NAI_Kxy5z_Px_bb+ABY*I_NAI_Ix5z_Px_bb;
  Double I_NAI_I5yz_Dxy_bb = I_NAI_K6yz_Px_bb+ABY*I_NAI_I5yz_Px_bb;
  Double I_NAI_I4y2z_Dxy_bb = I_NAI_K5y2z_Px_bb+ABY*I_NAI_I4y2z_Px_bb;
  Double I_NAI_I3y3z_Dxy_bb = I_NAI_K4y3z_Px_bb+ABY*I_NAI_I3y3z_Px_bb;
  Double I_NAI_I2y4z_Dxy_bb = I_NAI_K3y4z_Px_bb+ABY*I_NAI_I2y4z_Px_bb;
  Double I_NAI_Iy5z_Dxy_bb = I_NAI_K2y5z_Px_bb+ABY*I_NAI_Iy5z_Px_bb;
  Double I_NAI_I6z_Dxy_bb = I_NAI_Ky6z_Px_bb+ABY*I_NAI_I6z_Px_bb;
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
   * totally 0 integrals are omitted 
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
  Double I_NAI_H5x_F2xy_bb = I_NAI_I5xy_D2x_bb+ABY*I_NAI_H5x_D2x_bb;
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
  Double I_NAI_H5x_F2xz_bb = I_NAI_I5xz_D2x_bb+ABZ*I_NAI_H5x_D2x_bb;
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
  Double I_NAI_H5x_Fx2y_bb = I_NAI_I6x_D2y_bb+ABX*I_NAI_H5x_D2y_bb;
  Double I_NAI_H4xy_Fx2y_bb = I_NAI_I5xy_D2y_bb+ABX*I_NAI_H4xy_D2y_bb;
  Double I_NAI_H4xz_Fx2y_bb = I_NAI_I5xz_D2y_bb+ABX*I_NAI_H4xz_D2y_bb;
  Double I_NAI_H3x2y_Fx2y_bb = I_NAI_I4x2y_D2y_bb+ABX*I_NAI_H3x2y_D2y_bb;
  Double I_NAI_H3xyz_Fx2y_bb = I_NAI_I4xyz_D2y_bb+ABX*I_NAI_H3xyz_D2y_bb;
  Double I_NAI_H3x2z_Fx2y_bb = I_NAI_I4x2z_D2y_bb+ABX*I_NAI_H3x2z_D2y_bb;
  Double I_NAI_H2x3y_Fx2y_bb = I_NAI_I3x3y_D2y_bb+ABX*I_NAI_H2x3y_D2y_bb;
  Double I_NAI_H2x2yz_Fx2y_bb = I_NAI_I3x2yz_D2y_bb+ABX*I_NAI_H2x2yz_D2y_bb;
  Double I_NAI_H2xy2z_Fx2y_bb = I_NAI_I3xy2z_D2y_bb+ABX*I_NAI_H2xy2z_D2y_bb;
  Double I_NAI_H2x3z_Fx2y_bb = I_NAI_I3x3z_D2y_bb+ABX*I_NAI_H2x3z_D2y_bb;
  Double I_NAI_Hx4y_Fx2y_bb = I_NAI_I2x4y_D2y_bb+ABX*I_NAI_Hx4y_D2y_bb;
  Double I_NAI_Hx3yz_Fx2y_bb = I_NAI_I2x3yz_D2y_bb+ABX*I_NAI_Hx3yz_D2y_bb;
  Double I_NAI_Hx2y2z_Fx2y_bb = I_NAI_I2x2y2z_D2y_bb+ABX*I_NAI_Hx2y2z_D2y_bb;
  Double I_NAI_Hxy3z_Fx2y_bb = I_NAI_I2xy3z_D2y_bb+ABX*I_NAI_Hxy3z_D2y_bb;
  Double I_NAI_Hx4z_Fx2y_bb = I_NAI_I2x4z_D2y_bb+ABX*I_NAI_Hx4z_D2y_bb;
  Double I_NAI_H5y_Fx2y_bb = I_NAI_Ix5y_D2y_bb+ABX*I_NAI_H5y_D2y_bb;
  Double I_NAI_H4yz_Fx2y_bb = I_NAI_Ix4yz_D2y_bb+ABX*I_NAI_H4yz_D2y_bb;
  Double I_NAI_H3y2z_Fx2y_bb = I_NAI_Ix3y2z_D2y_bb+ABX*I_NAI_H3y2z_D2y_bb;
  Double I_NAI_H2y3z_Fx2y_bb = I_NAI_Ix2y3z_D2y_bb+ABX*I_NAI_H2y3z_D2y_bb;
  Double I_NAI_Hy4z_Fx2y_bb = I_NAI_Ixy4z_D2y_bb+ABX*I_NAI_Hy4z_D2y_bb;
  Double I_NAI_H5z_Fx2y_bb = I_NAI_Ix5z_D2y_bb+ABX*I_NAI_H5z_D2y_bb;
  Double I_NAI_H5x_Fxyz_bb = I_NAI_I5xz_Dxy_bb+ABZ*I_NAI_H5x_Dxy_bb;
  Double I_NAI_H4xy_Fxyz_bb = I_NAI_I4xyz_Dxy_bb+ABZ*I_NAI_H4xy_Dxy_bb;
  Double I_NAI_H4xz_Fxyz_bb = I_NAI_I4x2z_Dxy_bb+ABZ*I_NAI_H4xz_Dxy_bb;
  Double I_NAI_H3x2y_Fxyz_bb = I_NAI_I3x2yz_Dxy_bb+ABZ*I_NAI_H3x2y_Dxy_bb;
  Double I_NAI_H3xyz_Fxyz_bb = I_NAI_I3xy2z_Dxy_bb+ABZ*I_NAI_H3xyz_Dxy_bb;
  Double I_NAI_H3x2z_Fxyz_bb = I_NAI_I3x3z_Dxy_bb+ABZ*I_NAI_H3x2z_Dxy_bb;
  Double I_NAI_H2x3y_Fxyz_bb = I_NAI_I2x3yz_Dxy_bb+ABZ*I_NAI_H2x3y_Dxy_bb;
  Double I_NAI_H2x2yz_Fxyz_bb = I_NAI_I2x2y2z_Dxy_bb+ABZ*I_NAI_H2x2yz_Dxy_bb;
  Double I_NAI_H2xy2z_Fxyz_bb = I_NAI_I2xy3z_Dxy_bb+ABZ*I_NAI_H2xy2z_Dxy_bb;
  Double I_NAI_H2x3z_Fxyz_bb = I_NAI_I2x4z_Dxy_bb+ABZ*I_NAI_H2x3z_Dxy_bb;
  Double I_NAI_Hx4y_Fxyz_bb = I_NAI_Ix4yz_Dxy_bb+ABZ*I_NAI_Hx4y_Dxy_bb;
  Double I_NAI_Hx3yz_Fxyz_bb = I_NAI_Ix3y2z_Dxy_bb+ABZ*I_NAI_Hx3yz_Dxy_bb;
  Double I_NAI_Hx2y2z_Fxyz_bb = I_NAI_Ix2y3z_Dxy_bb+ABZ*I_NAI_Hx2y2z_Dxy_bb;
  Double I_NAI_Hxy3z_Fxyz_bb = I_NAI_Ixy4z_Dxy_bb+ABZ*I_NAI_Hxy3z_Dxy_bb;
  Double I_NAI_Hx4z_Fxyz_bb = I_NAI_Ix5z_Dxy_bb+ABZ*I_NAI_Hx4z_Dxy_bb;
  Double I_NAI_H5y_Fxyz_bb = I_NAI_I5yz_Dxy_bb+ABZ*I_NAI_H5y_Dxy_bb;
  Double I_NAI_H4yz_Fxyz_bb = I_NAI_I4y2z_Dxy_bb+ABZ*I_NAI_H4yz_Dxy_bb;
  Double I_NAI_H3y2z_Fxyz_bb = I_NAI_I3y3z_Dxy_bb+ABZ*I_NAI_H3y2z_Dxy_bb;
  Double I_NAI_H2y3z_Fxyz_bb = I_NAI_I2y4z_Dxy_bb+ABZ*I_NAI_H2y3z_Dxy_bb;
  Double I_NAI_Hy4z_Fxyz_bb = I_NAI_Iy5z_Dxy_bb+ABZ*I_NAI_Hy4z_Dxy_bb;
  Double I_NAI_H5z_Fxyz_bb = I_NAI_I6z_Dxy_bb+ABZ*I_NAI_H5z_Dxy_bb;
  Double I_NAI_H5x_Fx2z_bb = I_NAI_I6x_D2z_bb+ABX*I_NAI_H5x_D2z_bb;
  Double I_NAI_H4xy_Fx2z_bb = I_NAI_I5xy_D2z_bb+ABX*I_NAI_H4xy_D2z_bb;
  Double I_NAI_H4xz_Fx2z_bb = I_NAI_I5xz_D2z_bb+ABX*I_NAI_H4xz_D2z_bb;
  Double I_NAI_H3x2y_Fx2z_bb = I_NAI_I4x2y_D2z_bb+ABX*I_NAI_H3x2y_D2z_bb;
  Double I_NAI_H3xyz_Fx2z_bb = I_NAI_I4xyz_D2z_bb+ABX*I_NAI_H3xyz_D2z_bb;
  Double I_NAI_H3x2z_Fx2z_bb = I_NAI_I4x2z_D2z_bb+ABX*I_NAI_H3x2z_D2z_bb;
  Double I_NAI_H2x3y_Fx2z_bb = I_NAI_I3x3y_D2z_bb+ABX*I_NAI_H2x3y_D2z_bb;
  Double I_NAI_H2x2yz_Fx2z_bb = I_NAI_I3x2yz_D2z_bb+ABX*I_NAI_H2x2yz_D2z_bb;
  Double I_NAI_H2xy2z_Fx2z_bb = I_NAI_I3xy2z_D2z_bb+ABX*I_NAI_H2xy2z_D2z_bb;
  Double I_NAI_H2x3z_Fx2z_bb = I_NAI_I3x3z_D2z_bb+ABX*I_NAI_H2x3z_D2z_bb;
  Double I_NAI_Hx4y_Fx2z_bb = I_NAI_I2x4y_D2z_bb+ABX*I_NAI_Hx4y_D2z_bb;
  Double I_NAI_Hx3yz_Fx2z_bb = I_NAI_I2x3yz_D2z_bb+ABX*I_NAI_Hx3yz_D2z_bb;
  Double I_NAI_Hx2y2z_Fx2z_bb = I_NAI_I2x2y2z_D2z_bb+ABX*I_NAI_Hx2y2z_D2z_bb;
  Double I_NAI_Hxy3z_Fx2z_bb = I_NAI_I2xy3z_D2z_bb+ABX*I_NAI_Hxy3z_D2z_bb;
  Double I_NAI_Hx4z_Fx2z_bb = I_NAI_I2x4z_D2z_bb+ABX*I_NAI_Hx4z_D2z_bb;
  Double I_NAI_H5y_Fx2z_bb = I_NAI_Ix5y_D2z_bb+ABX*I_NAI_H5y_D2z_bb;
  Double I_NAI_H4yz_Fx2z_bb = I_NAI_Ix4yz_D2z_bb+ABX*I_NAI_H4yz_D2z_bb;
  Double I_NAI_H3y2z_Fx2z_bb = I_NAI_Ix3y2z_D2z_bb+ABX*I_NAI_H3y2z_D2z_bb;
  Double I_NAI_H2y3z_Fx2z_bb = I_NAI_Ix2y3z_D2z_bb+ABX*I_NAI_H2y3z_D2z_bb;
  Double I_NAI_Hy4z_Fx2z_bb = I_NAI_Ixy4z_D2z_bb+ABX*I_NAI_Hy4z_D2z_bb;
  Double I_NAI_H5z_Fx2z_bb = I_NAI_Ix5z_D2z_bb+ABX*I_NAI_H5z_D2z_bb;
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
  Double I_NAI_H5y_F2yz_bb = I_NAI_I5yz_D2y_bb+ABZ*I_NAI_H5y_D2y_bb;
  Double I_NAI_H4yz_F2yz_bb = I_NAI_I4y2z_D2y_bb+ABZ*I_NAI_H4yz_D2y_bb;
  Double I_NAI_H3y2z_F2yz_bb = I_NAI_I3y3z_D2y_bb+ABZ*I_NAI_H3y2z_D2y_bb;
  Double I_NAI_H2y3z_F2yz_bb = I_NAI_I2y4z_D2y_bb+ABZ*I_NAI_H2y3z_D2y_bb;
  Double I_NAI_Hy4z_F2yz_bb = I_NAI_Iy5z_D2y_bb+ABZ*I_NAI_Hy4z_D2y_bb;
  Double I_NAI_H5z_F2yz_bb = I_NAI_I6z_D2y_bb+ABZ*I_NAI_H5z_D2y_bb;
  Double I_NAI_H5x_Fy2z_bb = I_NAI_I5xy_D2z_bb+ABY*I_NAI_H5x_D2z_bb;
  Double I_NAI_H4xy_Fy2z_bb = I_NAI_I4x2y_D2z_bb+ABY*I_NAI_H4xy_D2z_bb;
  Double I_NAI_H4xz_Fy2z_bb = I_NAI_I4xyz_D2z_bb+ABY*I_NAI_H4xz_D2z_bb;
  Double I_NAI_H3x2y_Fy2z_bb = I_NAI_I3x3y_D2z_bb+ABY*I_NAI_H3x2y_D2z_bb;
  Double I_NAI_H3xyz_Fy2z_bb = I_NAI_I3x2yz_D2z_bb+ABY*I_NAI_H3xyz_D2z_bb;
  Double I_NAI_H3x2z_Fy2z_bb = I_NAI_I3xy2z_D2z_bb+ABY*I_NAI_H3x2z_D2z_bb;
  Double I_NAI_H2x3y_Fy2z_bb = I_NAI_I2x4y_D2z_bb+ABY*I_NAI_H2x3y_D2z_bb;
  Double I_NAI_H2x2yz_Fy2z_bb = I_NAI_I2x3yz_D2z_bb+ABY*I_NAI_H2x2yz_D2z_bb;
  Double I_NAI_H2xy2z_Fy2z_bb = I_NAI_I2x2y2z_D2z_bb+ABY*I_NAI_H2xy2z_D2z_bb;
  Double I_NAI_H2x3z_Fy2z_bb = I_NAI_I2xy3z_D2z_bb+ABY*I_NAI_H2x3z_D2z_bb;
  Double I_NAI_Hx4y_Fy2z_bb = I_NAI_Ix5y_D2z_bb+ABY*I_NAI_Hx4y_D2z_bb;
  Double I_NAI_Hx3yz_Fy2z_bb = I_NAI_Ix4yz_D2z_bb+ABY*I_NAI_Hx3yz_D2z_bb;
  Double I_NAI_Hx2y2z_Fy2z_bb = I_NAI_Ix3y2z_D2z_bb+ABY*I_NAI_Hx2y2z_D2z_bb;
  Double I_NAI_Hxy3z_Fy2z_bb = I_NAI_Ix2y3z_D2z_bb+ABY*I_NAI_Hxy3z_D2z_bb;
  Double I_NAI_Hx4z_Fy2z_bb = I_NAI_Ixy4z_D2z_bb+ABY*I_NAI_Hx4z_D2z_bb;
  Double I_NAI_H5y_Fy2z_bb = I_NAI_I6y_D2z_bb+ABY*I_NAI_H5y_D2z_bb;
  Double I_NAI_H4yz_Fy2z_bb = I_NAI_I5yz_D2z_bb+ABY*I_NAI_H4yz_D2z_bb;
  Double I_NAI_H3y2z_Fy2z_bb = I_NAI_I4y2z_D2z_bb+ABY*I_NAI_H3y2z_D2z_bb;
  Double I_NAI_H2y3z_Fy2z_bb = I_NAI_I3y3z_D2z_bb+ABY*I_NAI_H2y3z_D2z_bb;
  Double I_NAI_Hy4z_Fy2z_bb = I_NAI_I2y4z_D2z_bb+ABY*I_NAI_Hy4z_D2z_bb;
  Double I_NAI_H5z_Fy2z_bb = I_NAI_Iy5z_D2z_bb+ABY*I_NAI_H5z_D2z_bb;
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
   * shell quartet name: SQ_NAI_H_P_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[0] = 4.0E0*I_NAI_K7x_Px_aa-2.0E0*5*I_NAI_H5x_Px_a-2.0E0*6*I_NAI_H5x_Px_a+5*4*I_NAI_F3x_Px;
  abcd[1] = 4.0E0*I_NAI_K6xy_Px_aa-2.0E0*4*I_NAI_H4xy_Px_a-2.0E0*5*I_NAI_H4xy_Px_a+4*3*I_NAI_F2xy_Px;
  abcd[2] = 4.0E0*I_NAI_K6xz_Px_aa-2.0E0*4*I_NAI_H4xz_Px_a-2.0E0*5*I_NAI_H4xz_Px_a+4*3*I_NAI_F2xz_Px;
  abcd[3] = 4.0E0*I_NAI_K5x2y_Px_aa-2.0E0*3*I_NAI_H3x2y_Px_a-2.0E0*4*I_NAI_H3x2y_Px_a+3*2*I_NAI_Fx2y_Px;
  abcd[4] = 4.0E0*I_NAI_K5xyz_Px_aa-2.0E0*3*I_NAI_H3xyz_Px_a-2.0E0*4*I_NAI_H3xyz_Px_a+3*2*I_NAI_Fxyz_Px;
  abcd[5] = 4.0E0*I_NAI_K5x2z_Px_aa-2.0E0*3*I_NAI_H3x2z_Px_a-2.0E0*4*I_NAI_H3x2z_Px_a+3*2*I_NAI_Fx2z_Px;
  abcd[6] = 4.0E0*I_NAI_K4x3y_Px_aa-2.0E0*2*I_NAI_H2x3y_Px_a-2.0E0*3*I_NAI_H2x3y_Px_a+2*1*I_NAI_F3y_Px;
  abcd[7] = 4.0E0*I_NAI_K4x2yz_Px_aa-2.0E0*2*I_NAI_H2x2yz_Px_a-2.0E0*3*I_NAI_H2x2yz_Px_a+2*1*I_NAI_F2yz_Px;
  abcd[8] = 4.0E0*I_NAI_K4xy2z_Px_aa-2.0E0*2*I_NAI_H2xy2z_Px_a-2.0E0*3*I_NAI_H2xy2z_Px_a+2*1*I_NAI_Fy2z_Px;
  abcd[9] = 4.0E0*I_NAI_K4x3z_Px_aa-2.0E0*2*I_NAI_H2x3z_Px_a-2.0E0*3*I_NAI_H2x3z_Px_a+2*1*I_NAI_F3z_Px;
  abcd[10] = 4.0E0*I_NAI_K3x4y_Px_aa-2.0E0*1*I_NAI_Hx4y_Px_a-2.0E0*2*I_NAI_Hx4y_Px_a;
  abcd[11] = 4.0E0*I_NAI_K3x3yz_Px_aa-2.0E0*1*I_NAI_Hx3yz_Px_a-2.0E0*2*I_NAI_Hx3yz_Px_a;
  abcd[12] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_Hx2y2z_Px_a-2.0E0*2*I_NAI_Hx2y2z_Px_a;
  abcd[13] = 4.0E0*I_NAI_K3xy3z_Px_aa-2.0E0*1*I_NAI_Hxy3z_Px_a-2.0E0*2*I_NAI_Hxy3z_Px_a;
  abcd[14] = 4.0E0*I_NAI_K3x4z_Px_aa-2.0E0*1*I_NAI_Hx4z_Px_a-2.0E0*2*I_NAI_Hx4z_Px_a;
  abcd[15] = 4.0E0*I_NAI_K2x5y_Px_aa-2.0E0*1*I_NAI_H5y_Px_a;
  abcd[16] = 4.0E0*I_NAI_K2x4yz_Px_aa-2.0E0*1*I_NAI_H4yz_Px_a;
  abcd[17] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*1*I_NAI_H3y2z_Px_a;
  abcd[18] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*1*I_NAI_H2y3z_Px_a;
  abcd[19] = 4.0E0*I_NAI_K2xy4z_Px_aa-2.0E0*1*I_NAI_Hy4z_Px_a;
  abcd[20] = 4.0E0*I_NAI_K2x5z_Px_aa-2.0E0*1*I_NAI_H5z_Px_a;
  abcd[21] = 4.0E0*I_NAI_K7x_Py_aa-2.0E0*5*I_NAI_H5x_Py_a-2.0E0*6*I_NAI_H5x_Py_a+5*4*I_NAI_F3x_Py;
  abcd[22] = 4.0E0*I_NAI_K6xy_Py_aa-2.0E0*4*I_NAI_H4xy_Py_a-2.0E0*5*I_NAI_H4xy_Py_a+4*3*I_NAI_F2xy_Py;
  abcd[23] = 4.0E0*I_NAI_K6xz_Py_aa-2.0E0*4*I_NAI_H4xz_Py_a-2.0E0*5*I_NAI_H4xz_Py_a+4*3*I_NAI_F2xz_Py;
  abcd[24] = 4.0E0*I_NAI_K5x2y_Py_aa-2.0E0*3*I_NAI_H3x2y_Py_a-2.0E0*4*I_NAI_H3x2y_Py_a+3*2*I_NAI_Fx2y_Py;
  abcd[25] = 4.0E0*I_NAI_K5xyz_Py_aa-2.0E0*3*I_NAI_H3xyz_Py_a-2.0E0*4*I_NAI_H3xyz_Py_a+3*2*I_NAI_Fxyz_Py;
  abcd[26] = 4.0E0*I_NAI_K5x2z_Py_aa-2.0E0*3*I_NAI_H3x2z_Py_a-2.0E0*4*I_NAI_H3x2z_Py_a+3*2*I_NAI_Fx2z_Py;
  abcd[27] = 4.0E0*I_NAI_K4x3y_Py_aa-2.0E0*2*I_NAI_H2x3y_Py_a-2.0E0*3*I_NAI_H2x3y_Py_a+2*1*I_NAI_F3y_Py;
  abcd[28] = 4.0E0*I_NAI_K4x2yz_Py_aa-2.0E0*2*I_NAI_H2x2yz_Py_a-2.0E0*3*I_NAI_H2x2yz_Py_a+2*1*I_NAI_F2yz_Py;
  abcd[29] = 4.0E0*I_NAI_K4xy2z_Py_aa-2.0E0*2*I_NAI_H2xy2z_Py_a-2.0E0*3*I_NAI_H2xy2z_Py_a+2*1*I_NAI_Fy2z_Py;
  abcd[30] = 4.0E0*I_NAI_K4x3z_Py_aa-2.0E0*2*I_NAI_H2x3z_Py_a-2.0E0*3*I_NAI_H2x3z_Py_a+2*1*I_NAI_F3z_Py;
  abcd[31] = 4.0E0*I_NAI_K3x4y_Py_aa-2.0E0*1*I_NAI_Hx4y_Py_a-2.0E0*2*I_NAI_Hx4y_Py_a;
  abcd[32] = 4.0E0*I_NAI_K3x3yz_Py_aa-2.0E0*1*I_NAI_Hx3yz_Py_a-2.0E0*2*I_NAI_Hx3yz_Py_a;
  abcd[33] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_Hx2y2z_Py_a-2.0E0*2*I_NAI_Hx2y2z_Py_a;
  abcd[34] = 4.0E0*I_NAI_K3xy3z_Py_aa-2.0E0*1*I_NAI_Hxy3z_Py_a-2.0E0*2*I_NAI_Hxy3z_Py_a;
  abcd[35] = 4.0E0*I_NAI_K3x4z_Py_aa-2.0E0*1*I_NAI_Hx4z_Py_a-2.0E0*2*I_NAI_Hx4z_Py_a;
  abcd[36] = 4.0E0*I_NAI_K2x5y_Py_aa-2.0E0*1*I_NAI_H5y_Py_a;
  abcd[37] = 4.0E0*I_NAI_K2x4yz_Py_aa-2.0E0*1*I_NAI_H4yz_Py_a;
  abcd[38] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*1*I_NAI_H3y2z_Py_a;
  abcd[39] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*1*I_NAI_H2y3z_Py_a;
  abcd[40] = 4.0E0*I_NAI_K2xy4z_Py_aa-2.0E0*1*I_NAI_Hy4z_Py_a;
  abcd[41] = 4.0E0*I_NAI_K2x5z_Py_aa-2.0E0*1*I_NAI_H5z_Py_a;
  abcd[42] = 4.0E0*I_NAI_K7x_Pz_aa-2.0E0*5*I_NAI_H5x_Pz_a-2.0E0*6*I_NAI_H5x_Pz_a+5*4*I_NAI_F3x_Pz;
  abcd[43] = 4.0E0*I_NAI_K6xy_Pz_aa-2.0E0*4*I_NAI_H4xy_Pz_a-2.0E0*5*I_NAI_H4xy_Pz_a+4*3*I_NAI_F2xy_Pz;
  abcd[44] = 4.0E0*I_NAI_K6xz_Pz_aa-2.0E0*4*I_NAI_H4xz_Pz_a-2.0E0*5*I_NAI_H4xz_Pz_a+4*3*I_NAI_F2xz_Pz;
  abcd[45] = 4.0E0*I_NAI_K5x2y_Pz_aa-2.0E0*3*I_NAI_H3x2y_Pz_a-2.0E0*4*I_NAI_H3x2y_Pz_a+3*2*I_NAI_Fx2y_Pz;
  abcd[46] = 4.0E0*I_NAI_K5xyz_Pz_aa-2.0E0*3*I_NAI_H3xyz_Pz_a-2.0E0*4*I_NAI_H3xyz_Pz_a+3*2*I_NAI_Fxyz_Pz;
  abcd[47] = 4.0E0*I_NAI_K5x2z_Pz_aa-2.0E0*3*I_NAI_H3x2z_Pz_a-2.0E0*4*I_NAI_H3x2z_Pz_a+3*2*I_NAI_Fx2z_Pz;
  abcd[48] = 4.0E0*I_NAI_K4x3y_Pz_aa-2.0E0*2*I_NAI_H2x3y_Pz_a-2.0E0*3*I_NAI_H2x3y_Pz_a+2*1*I_NAI_F3y_Pz;
  abcd[49] = 4.0E0*I_NAI_K4x2yz_Pz_aa-2.0E0*2*I_NAI_H2x2yz_Pz_a-2.0E0*3*I_NAI_H2x2yz_Pz_a+2*1*I_NAI_F2yz_Pz;
  abcd[50] = 4.0E0*I_NAI_K4xy2z_Pz_aa-2.0E0*2*I_NAI_H2xy2z_Pz_a-2.0E0*3*I_NAI_H2xy2z_Pz_a+2*1*I_NAI_Fy2z_Pz;
  abcd[51] = 4.0E0*I_NAI_K4x3z_Pz_aa-2.0E0*2*I_NAI_H2x3z_Pz_a-2.0E0*3*I_NAI_H2x3z_Pz_a+2*1*I_NAI_F3z_Pz;
  abcd[52] = 4.0E0*I_NAI_K3x4y_Pz_aa-2.0E0*1*I_NAI_Hx4y_Pz_a-2.0E0*2*I_NAI_Hx4y_Pz_a;
  abcd[53] = 4.0E0*I_NAI_K3x3yz_Pz_aa-2.0E0*1*I_NAI_Hx3yz_Pz_a-2.0E0*2*I_NAI_Hx3yz_Pz_a;
  abcd[54] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_Hx2y2z_Pz_a-2.0E0*2*I_NAI_Hx2y2z_Pz_a;
  abcd[55] = 4.0E0*I_NAI_K3xy3z_Pz_aa-2.0E0*1*I_NAI_Hxy3z_Pz_a-2.0E0*2*I_NAI_Hxy3z_Pz_a;
  abcd[56] = 4.0E0*I_NAI_K3x4z_Pz_aa-2.0E0*1*I_NAI_Hx4z_Pz_a-2.0E0*2*I_NAI_Hx4z_Pz_a;
  abcd[57] = 4.0E0*I_NAI_K2x5y_Pz_aa-2.0E0*1*I_NAI_H5y_Pz_a;
  abcd[58] = 4.0E0*I_NAI_K2x4yz_Pz_aa-2.0E0*1*I_NAI_H4yz_Pz_a;
  abcd[59] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*1*I_NAI_H3y2z_Pz_a;
  abcd[60] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*1*I_NAI_H2y3z_Pz_a;
  abcd[61] = 4.0E0*I_NAI_K2xy4z_Pz_aa-2.0E0*1*I_NAI_Hy4z_Pz_a;
  abcd[62] = 4.0E0*I_NAI_K2x5z_Pz_aa-2.0E0*1*I_NAI_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[63] = 4.0E0*I_NAI_K6xy_Px_aa-2.0E0*5*I_NAI_H4xy_Px_a;
  abcd[64] = 4.0E0*I_NAI_K5x2y_Px_aa-2.0E0*1*I_NAI_H5x_Px_a-2.0E0*4*I_NAI_H3x2y_Px_a+4*1*I_NAI_F3x_Px;
  abcd[65] = 4.0E0*I_NAI_K5xyz_Px_aa-2.0E0*4*I_NAI_H3xyz_Px_a;
  abcd[66] = 4.0E0*I_NAI_K4x3y_Px_aa-2.0E0*2*I_NAI_H4xy_Px_a-2.0E0*3*I_NAI_H2x3y_Px_a+3*2*I_NAI_F2xy_Px;
  abcd[67] = 4.0E0*I_NAI_K4x2yz_Px_aa-2.0E0*1*I_NAI_H4xz_Px_a-2.0E0*3*I_NAI_H2x2yz_Px_a+3*1*I_NAI_F2xz_Px;
  abcd[68] = 4.0E0*I_NAI_K4xy2z_Px_aa-2.0E0*3*I_NAI_H2xy2z_Px_a;
  abcd[69] = 4.0E0*I_NAI_K3x4y_Px_aa-2.0E0*3*I_NAI_H3x2y_Px_a-2.0E0*2*I_NAI_Hx4y_Px_a+2*3*I_NAI_Fx2y_Px;
  abcd[70] = 4.0E0*I_NAI_K3x3yz_Px_aa-2.0E0*2*I_NAI_H3xyz_Px_a-2.0E0*2*I_NAI_Hx3yz_Px_a+2*2*I_NAI_Fxyz_Px;
  abcd[71] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_H3x2z_Px_a-2.0E0*2*I_NAI_Hx2y2z_Px_a+2*1*I_NAI_Fx2z_Px;
  abcd[72] = 4.0E0*I_NAI_K3xy3z_Px_aa-2.0E0*2*I_NAI_Hxy3z_Px_a;
  abcd[73] = 4.0E0*I_NAI_K2x5y_Px_aa-2.0E0*4*I_NAI_H2x3y_Px_a-2.0E0*1*I_NAI_H5y_Px_a+4*I_NAI_F3y_Px;
  abcd[74] = 4.0E0*I_NAI_K2x4yz_Px_aa-2.0E0*3*I_NAI_H2x2yz_Px_a-2.0E0*1*I_NAI_H4yz_Px_a+3*I_NAI_F2yz_Px;
  abcd[75] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*2*I_NAI_H2xy2z_Px_a-2.0E0*1*I_NAI_H3y2z_Px_a+2*I_NAI_Fy2z_Px;
  abcd[76] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*1*I_NAI_H2x3z_Px_a-2.0E0*1*I_NAI_H2y3z_Px_a+1*I_NAI_F3z_Px;
  abcd[77] = 4.0E0*I_NAI_K2xy4z_Px_aa-2.0E0*1*I_NAI_Hy4z_Px_a;
  abcd[78] = 4.0E0*I_NAI_Kx6y_Px_aa-2.0E0*5*I_NAI_Hx4y_Px_a;
  abcd[79] = 4.0E0*I_NAI_Kx5yz_Px_aa-2.0E0*4*I_NAI_Hx3yz_Px_a;
  abcd[80] = 4.0E0*I_NAI_Kx4y2z_Px_aa-2.0E0*3*I_NAI_Hx2y2z_Px_a;
  abcd[81] = 4.0E0*I_NAI_Kx3y3z_Px_aa-2.0E0*2*I_NAI_Hxy3z_Px_a;
  abcd[82] = 4.0E0*I_NAI_Kx2y4z_Px_aa-2.0E0*1*I_NAI_Hx4z_Px_a;
  abcd[83] = 4.0E0*I_NAI_Kxy5z_Px_aa;
  abcd[84] = 4.0E0*I_NAI_K6xy_Py_aa-2.0E0*5*I_NAI_H4xy_Py_a;
  abcd[85] = 4.0E0*I_NAI_K5x2y_Py_aa-2.0E0*1*I_NAI_H5x_Py_a-2.0E0*4*I_NAI_H3x2y_Py_a+4*1*I_NAI_F3x_Py;
  abcd[86] = 4.0E0*I_NAI_K5xyz_Py_aa-2.0E0*4*I_NAI_H3xyz_Py_a;
  abcd[87] = 4.0E0*I_NAI_K4x3y_Py_aa-2.0E0*2*I_NAI_H4xy_Py_a-2.0E0*3*I_NAI_H2x3y_Py_a+3*2*I_NAI_F2xy_Py;
  abcd[88] = 4.0E0*I_NAI_K4x2yz_Py_aa-2.0E0*1*I_NAI_H4xz_Py_a-2.0E0*3*I_NAI_H2x2yz_Py_a+3*1*I_NAI_F2xz_Py;
  abcd[89] = 4.0E0*I_NAI_K4xy2z_Py_aa-2.0E0*3*I_NAI_H2xy2z_Py_a;
  abcd[90] = 4.0E0*I_NAI_K3x4y_Py_aa-2.0E0*3*I_NAI_H3x2y_Py_a-2.0E0*2*I_NAI_Hx4y_Py_a+2*3*I_NAI_Fx2y_Py;
  abcd[91] = 4.0E0*I_NAI_K3x3yz_Py_aa-2.0E0*2*I_NAI_H3xyz_Py_a-2.0E0*2*I_NAI_Hx3yz_Py_a+2*2*I_NAI_Fxyz_Py;
  abcd[92] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_H3x2z_Py_a-2.0E0*2*I_NAI_Hx2y2z_Py_a+2*1*I_NAI_Fx2z_Py;
  abcd[93] = 4.0E0*I_NAI_K3xy3z_Py_aa-2.0E0*2*I_NAI_Hxy3z_Py_a;
  abcd[94] = 4.0E0*I_NAI_K2x5y_Py_aa-2.0E0*4*I_NAI_H2x3y_Py_a-2.0E0*1*I_NAI_H5y_Py_a+4*I_NAI_F3y_Py;
  abcd[95] = 4.0E0*I_NAI_K2x4yz_Py_aa-2.0E0*3*I_NAI_H2x2yz_Py_a-2.0E0*1*I_NAI_H4yz_Py_a+3*I_NAI_F2yz_Py;
  abcd[96] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*2*I_NAI_H2xy2z_Py_a-2.0E0*1*I_NAI_H3y2z_Py_a+2*I_NAI_Fy2z_Py;
  abcd[97] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*1*I_NAI_H2x3z_Py_a-2.0E0*1*I_NAI_H2y3z_Py_a+1*I_NAI_F3z_Py;
  abcd[98] = 4.0E0*I_NAI_K2xy4z_Py_aa-2.0E0*1*I_NAI_Hy4z_Py_a;
  abcd[99] = 4.0E0*I_NAI_Kx6y_Py_aa-2.0E0*5*I_NAI_Hx4y_Py_a;
  abcd[100] = 4.0E0*I_NAI_Kx5yz_Py_aa-2.0E0*4*I_NAI_Hx3yz_Py_a;
  abcd[101] = 4.0E0*I_NAI_Kx4y2z_Py_aa-2.0E0*3*I_NAI_Hx2y2z_Py_a;
  abcd[102] = 4.0E0*I_NAI_Kx3y3z_Py_aa-2.0E0*2*I_NAI_Hxy3z_Py_a;
  abcd[103] = 4.0E0*I_NAI_Kx2y4z_Py_aa-2.0E0*1*I_NAI_Hx4z_Py_a;
  abcd[104] = 4.0E0*I_NAI_Kxy5z_Py_aa;
  abcd[105] = 4.0E0*I_NAI_K6xy_Pz_aa-2.0E0*5*I_NAI_H4xy_Pz_a;
  abcd[106] = 4.0E0*I_NAI_K5x2y_Pz_aa-2.0E0*1*I_NAI_H5x_Pz_a-2.0E0*4*I_NAI_H3x2y_Pz_a+4*1*I_NAI_F3x_Pz;
  abcd[107] = 4.0E0*I_NAI_K5xyz_Pz_aa-2.0E0*4*I_NAI_H3xyz_Pz_a;
  abcd[108] = 4.0E0*I_NAI_K4x3y_Pz_aa-2.0E0*2*I_NAI_H4xy_Pz_a-2.0E0*3*I_NAI_H2x3y_Pz_a+3*2*I_NAI_F2xy_Pz;
  abcd[109] = 4.0E0*I_NAI_K4x2yz_Pz_aa-2.0E0*1*I_NAI_H4xz_Pz_a-2.0E0*3*I_NAI_H2x2yz_Pz_a+3*1*I_NAI_F2xz_Pz;
  abcd[110] = 4.0E0*I_NAI_K4xy2z_Pz_aa-2.0E0*3*I_NAI_H2xy2z_Pz_a;
  abcd[111] = 4.0E0*I_NAI_K3x4y_Pz_aa-2.0E0*3*I_NAI_H3x2y_Pz_a-2.0E0*2*I_NAI_Hx4y_Pz_a+2*3*I_NAI_Fx2y_Pz;
  abcd[112] = 4.0E0*I_NAI_K3x3yz_Pz_aa-2.0E0*2*I_NAI_H3xyz_Pz_a-2.0E0*2*I_NAI_Hx3yz_Pz_a+2*2*I_NAI_Fxyz_Pz;
  abcd[113] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_H3x2z_Pz_a-2.0E0*2*I_NAI_Hx2y2z_Pz_a+2*1*I_NAI_Fx2z_Pz;
  abcd[114] = 4.0E0*I_NAI_K3xy3z_Pz_aa-2.0E0*2*I_NAI_Hxy3z_Pz_a;
  abcd[115] = 4.0E0*I_NAI_K2x5y_Pz_aa-2.0E0*4*I_NAI_H2x3y_Pz_a-2.0E0*1*I_NAI_H5y_Pz_a+4*I_NAI_F3y_Pz;
  abcd[116] = 4.0E0*I_NAI_K2x4yz_Pz_aa-2.0E0*3*I_NAI_H2x2yz_Pz_a-2.0E0*1*I_NAI_H4yz_Pz_a+3*I_NAI_F2yz_Pz;
  abcd[117] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*2*I_NAI_H2xy2z_Pz_a-2.0E0*1*I_NAI_H3y2z_Pz_a+2*I_NAI_Fy2z_Pz;
  abcd[118] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*1*I_NAI_H2x3z_Pz_a-2.0E0*1*I_NAI_H2y3z_Pz_a+1*I_NAI_F3z_Pz;
  abcd[119] = 4.0E0*I_NAI_K2xy4z_Pz_aa-2.0E0*1*I_NAI_Hy4z_Pz_a;
  abcd[120] = 4.0E0*I_NAI_Kx6y_Pz_aa-2.0E0*5*I_NAI_Hx4y_Pz_a;
  abcd[121] = 4.0E0*I_NAI_Kx5yz_Pz_aa-2.0E0*4*I_NAI_Hx3yz_Pz_a;
  abcd[122] = 4.0E0*I_NAI_Kx4y2z_Pz_aa-2.0E0*3*I_NAI_Hx2y2z_Pz_a;
  abcd[123] = 4.0E0*I_NAI_Kx3y3z_Pz_aa-2.0E0*2*I_NAI_Hxy3z_Pz_a;
  abcd[124] = 4.0E0*I_NAI_Kx2y4z_Pz_aa-2.0E0*1*I_NAI_Hx4z_Pz_a;
  abcd[125] = 4.0E0*I_NAI_Kxy5z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[126] = 4.0E0*I_NAI_K6xz_Px_aa-2.0E0*5*I_NAI_H4xz_Px_a;
  abcd[127] = 4.0E0*I_NAI_K5xyz_Px_aa-2.0E0*4*I_NAI_H3xyz_Px_a;
  abcd[128] = 4.0E0*I_NAI_K5x2z_Px_aa-2.0E0*1*I_NAI_H5x_Px_a-2.0E0*4*I_NAI_H3x2z_Px_a+4*1*I_NAI_F3x_Px;
  abcd[129] = 4.0E0*I_NAI_K4x2yz_Px_aa-2.0E0*3*I_NAI_H2x2yz_Px_a;
  abcd[130] = 4.0E0*I_NAI_K4xy2z_Px_aa-2.0E0*1*I_NAI_H4xy_Px_a-2.0E0*3*I_NAI_H2xy2z_Px_a+3*1*I_NAI_F2xy_Px;
  abcd[131] = 4.0E0*I_NAI_K4x3z_Px_aa-2.0E0*2*I_NAI_H4xz_Px_a-2.0E0*3*I_NAI_H2x3z_Px_a+3*2*I_NAI_F2xz_Px;
  abcd[132] = 4.0E0*I_NAI_K3x3yz_Px_aa-2.0E0*2*I_NAI_Hx3yz_Px_a;
  abcd[133] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_H3x2y_Px_a-2.0E0*2*I_NAI_Hx2y2z_Px_a+2*1*I_NAI_Fx2y_Px;
  abcd[134] = 4.0E0*I_NAI_K3xy3z_Px_aa-2.0E0*2*I_NAI_H3xyz_Px_a-2.0E0*2*I_NAI_Hxy3z_Px_a+2*2*I_NAI_Fxyz_Px;
  abcd[135] = 4.0E0*I_NAI_K3x4z_Px_aa-2.0E0*3*I_NAI_H3x2z_Px_a-2.0E0*2*I_NAI_Hx4z_Px_a+2*3*I_NAI_Fx2z_Px;
  abcd[136] = 4.0E0*I_NAI_K2x4yz_Px_aa-2.0E0*1*I_NAI_H4yz_Px_a;
  abcd[137] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*1*I_NAI_H2x3y_Px_a-2.0E0*1*I_NAI_H3y2z_Px_a+1*I_NAI_F3y_Px;
  abcd[138] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*2*I_NAI_H2x2yz_Px_a-2.0E0*1*I_NAI_H2y3z_Px_a+2*I_NAI_F2yz_Px;
  abcd[139] = 4.0E0*I_NAI_K2xy4z_Px_aa-2.0E0*3*I_NAI_H2xy2z_Px_a-2.0E0*1*I_NAI_Hy4z_Px_a+3*I_NAI_Fy2z_Px;
  abcd[140] = 4.0E0*I_NAI_K2x5z_Px_aa-2.0E0*4*I_NAI_H2x3z_Px_a-2.0E0*1*I_NAI_H5z_Px_a+4*I_NAI_F3z_Px;
  abcd[141] = 4.0E0*I_NAI_Kx5yz_Px_aa;
  abcd[142] = 4.0E0*I_NAI_Kx4y2z_Px_aa-2.0E0*1*I_NAI_Hx4y_Px_a;
  abcd[143] = 4.0E0*I_NAI_Kx3y3z_Px_aa-2.0E0*2*I_NAI_Hx3yz_Px_a;
  abcd[144] = 4.0E0*I_NAI_Kx2y4z_Px_aa-2.0E0*3*I_NAI_Hx2y2z_Px_a;
  abcd[145] = 4.0E0*I_NAI_Kxy5z_Px_aa-2.0E0*4*I_NAI_Hxy3z_Px_a;
  abcd[146] = 4.0E0*I_NAI_Kx6z_Px_aa-2.0E0*5*I_NAI_Hx4z_Px_a;
  abcd[147] = 4.0E0*I_NAI_K6xz_Py_aa-2.0E0*5*I_NAI_H4xz_Py_a;
  abcd[148] = 4.0E0*I_NAI_K5xyz_Py_aa-2.0E0*4*I_NAI_H3xyz_Py_a;
  abcd[149] = 4.0E0*I_NAI_K5x2z_Py_aa-2.0E0*1*I_NAI_H5x_Py_a-2.0E0*4*I_NAI_H3x2z_Py_a+4*1*I_NAI_F3x_Py;
  abcd[150] = 4.0E0*I_NAI_K4x2yz_Py_aa-2.0E0*3*I_NAI_H2x2yz_Py_a;
  abcd[151] = 4.0E0*I_NAI_K4xy2z_Py_aa-2.0E0*1*I_NAI_H4xy_Py_a-2.0E0*3*I_NAI_H2xy2z_Py_a+3*1*I_NAI_F2xy_Py;
  abcd[152] = 4.0E0*I_NAI_K4x3z_Py_aa-2.0E0*2*I_NAI_H4xz_Py_a-2.0E0*3*I_NAI_H2x3z_Py_a+3*2*I_NAI_F2xz_Py;
  abcd[153] = 4.0E0*I_NAI_K3x3yz_Py_aa-2.0E0*2*I_NAI_Hx3yz_Py_a;
  abcd[154] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_H3x2y_Py_a-2.0E0*2*I_NAI_Hx2y2z_Py_a+2*1*I_NAI_Fx2y_Py;
  abcd[155] = 4.0E0*I_NAI_K3xy3z_Py_aa-2.0E0*2*I_NAI_H3xyz_Py_a-2.0E0*2*I_NAI_Hxy3z_Py_a+2*2*I_NAI_Fxyz_Py;
  abcd[156] = 4.0E0*I_NAI_K3x4z_Py_aa-2.0E0*3*I_NAI_H3x2z_Py_a-2.0E0*2*I_NAI_Hx4z_Py_a+2*3*I_NAI_Fx2z_Py;
  abcd[157] = 4.0E0*I_NAI_K2x4yz_Py_aa-2.0E0*1*I_NAI_H4yz_Py_a;
  abcd[158] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*1*I_NAI_H2x3y_Py_a-2.0E0*1*I_NAI_H3y2z_Py_a+1*I_NAI_F3y_Py;
  abcd[159] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*2*I_NAI_H2x2yz_Py_a-2.0E0*1*I_NAI_H2y3z_Py_a+2*I_NAI_F2yz_Py;
  abcd[160] = 4.0E0*I_NAI_K2xy4z_Py_aa-2.0E0*3*I_NAI_H2xy2z_Py_a-2.0E0*1*I_NAI_Hy4z_Py_a+3*I_NAI_Fy2z_Py;
  abcd[161] = 4.0E0*I_NAI_K2x5z_Py_aa-2.0E0*4*I_NAI_H2x3z_Py_a-2.0E0*1*I_NAI_H5z_Py_a+4*I_NAI_F3z_Py;
  abcd[162] = 4.0E0*I_NAI_Kx5yz_Py_aa;
  abcd[163] = 4.0E0*I_NAI_Kx4y2z_Py_aa-2.0E0*1*I_NAI_Hx4y_Py_a;
  abcd[164] = 4.0E0*I_NAI_Kx3y3z_Py_aa-2.0E0*2*I_NAI_Hx3yz_Py_a;
  abcd[165] = 4.0E0*I_NAI_Kx2y4z_Py_aa-2.0E0*3*I_NAI_Hx2y2z_Py_a;
  abcd[166] = 4.0E0*I_NAI_Kxy5z_Py_aa-2.0E0*4*I_NAI_Hxy3z_Py_a;
  abcd[167] = 4.0E0*I_NAI_Kx6z_Py_aa-2.0E0*5*I_NAI_Hx4z_Py_a;
  abcd[168] = 4.0E0*I_NAI_K6xz_Pz_aa-2.0E0*5*I_NAI_H4xz_Pz_a;
  abcd[169] = 4.0E0*I_NAI_K5xyz_Pz_aa-2.0E0*4*I_NAI_H3xyz_Pz_a;
  abcd[170] = 4.0E0*I_NAI_K5x2z_Pz_aa-2.0E0*1*I_NAI_H5x_Pz_a-2.0E0*4*I_NAI_H3x2z_Pz_a+4*1*I_NAI_F3x_Pz;
  abcd[171] = 4.0E0*I_NAI_K4x2yz_Pz_aa-2.0E0*3*I_NAI_H2x2yz_Pz_a;
  abcd[172] = 4.0E0*I_NAI_K4xy2z_Pz_aa-2.0E0*1*I_NAI_H4xy_Pz_a-2.0E0*3*I_NAI_H2xy2z_Pz_a+3*1*I_NAI_F2xy_Pz;
  abcd[173] = 4.0E0*I_NAI_K4x3z_Pz_aa-2.0E0*2*I_NAI_H4xz_Pz_a-2.0E0*3*I_NAI_H2x3z_Pz_a+3*2*I_NAI_F2xz_Pz;
  abcd[174] = 4.0E0*I_NAI_K3x3yz_Pz_aa-2.0E0*2*I_NAI_Hx3yz_Pz_a;
  abcd[175] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_H3x2y_Pz_a-2.0E0*2*I_NAI_Hx2y2z_Pz_a+2*1*I_NAI_Fx2y_Pz;
  abcd[176] = 4.0E0*I_NAI_K3xy3z_Pz_aa-2.0E0*2*I_NAI_H3xyz_Pz_a-2.0E0*2*I_NAI_Hxy3z_Pz_a+2*2*I_NAI_Fxyz_Pz;
  abcd[177] = 4.0E0*I_NAI_K3x4z_Pz_aa-2.0E0*3*I_NAI_H3x2z_Pz_a-2.0E0*2*I_NAI_Hx4z_Pz_a+2*3*I_NAI_Fx2z_Pz;
  abcd[178] = 4.0E0*I_NAI_K2x4yz_Pz_aa-2.0E0*1*I_NAI_H4yz_Pz_a;
  abcd[179] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*1*I_NAI_H2x3y_Pz_a-2.0E0*1*I_NAI_H3y2z_Pz_a+1*I_NAI_F3y_Pz;
  abcd[180] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*2*I_NAI_H2x2yz_Pz_a-2.0E0*1*I_NAI_H2y3z_Pz_a+2*I_NAI_F2yz_Pz;
  abcd[181] = 4.0E0*I_NAI_K2xy4z_Pz_aa-2.0E0*3*I_NAI_H2xy2z_Pz_a-2.0E0*1*I_NAI_Hy4z_Pz_a+3*I_NAI_Fy2z_Pz;
  abcd[182] = 4.0E0*I_NAI_K2x5z_Pz_aa-2.0E0*4*I_NAI_H2x3z_Pz_a-2.0E0*1*I_NAI_H5z_Pz_a+4*I_NAI_F3z_Pz;
  abcd[183] = 4.0E0*I_NAI_Kx5yz_Pz_aa;
  abcd[184] = 4.0E0*I_NAI_Kx4y2z_Pz_aa-2.0E0*1*I_NAI_Hx4y_Pz_a;
  abcd[185] = 4.0E0*I_NAI_Kx3y3z_Pz_aa-2.0E0*2*I_NAI_Hx3yz_Pz_a;
  abcd[186] = 4.0E0*I_NAI_Kx2y4z_Pz_aa-2.0E0*3*I_NAI_Hx2y2z_Pz_a;
  abcd[187] = 4.0E0*I_NAI_Kxy5z_Pz_aa-2.0E0*4*I_NAI_Hxy3z_Pz_a;
  abcd[188] = 4.0E0*I_NAI_Kx6z_Pz_aa-2.0E0*5*I_NAI_Hx4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[189] = 4.0E0*I_NAI_K5x2y_Px_aa-2.0E0*1*I_NAI_H5x_Px_a;
  abcd[190] = 4.0E0*I_NAI_K4x3y_Px_aa-2.0E0*1*I_NAI_H4xy_Px_a-2.0E0*2*I_NAI_H4xy_Px_a;
  abcd[191] = 4.0E0*I_NAI_K4x2yz_Px_aa-2.0E0*1*I_NAI_H4xz_Px_a;
  abcd[192] = 4.0E0*I_NAI_K3x4y_Px_aa-2.0E0*2*I_NAI_H3x2y_Px_a-2.0E0*3*I_NAI_H3x2y_Px_a+2*1*I_NAI_F3x_Px;
  abcd[193] = 4.0E0*I_NAI_K3x3yz_Px_aa-2.0E0*1*I_NAI_H3xyz_Px_a-2.0E0*2*I_NAI_H3xyz_Px_a;
  abcd[194] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_H3x2z_Px_a;
  abcd[195] = 4.0E0*I_NAI_K2x5y_Px_aa-2.0E0*3*I_NAI_H2x3y_Px_a-2.0E0*4*I_NAI_H2x3y_Px_a+3*2*I_NAI_F2xy_Px;
  abcd[196] = 4.0E0*I_NAI_K2x4yz_Px_aa-2.0E0*2*I_NAI_H2x2yz_Px_a-2.0E0*3*I_NAI_H2x2yz_Px_a+2*1*I_NAI_F2xz_Px;
  abcd[197] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*1*I_NAI_H2xy2z_Px_a-2.0E0*2*I_NAI_H2xy2z_Px_a;
  abcd[198] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*1*I_NAI_H2x3z_Px_a;
  abcd[199] = 4.0E0*I_NAI_Kx6y_Px_aa-2.0E0*4*I_NAI_Hx4y_Px_a-2.0E0*5*I_NAI_Hx4y_Px_a+4*3*I_NAI_Fx2y_Px;
  abcd[200] = 4.0E0*I_NAI_Kx5yz_Px_aa-2.0E0*3*I_NAI_Hx3yz_Px_a-2.0E0*4*I_NAI_Hx3yz_Px_a+3*2*I_NAI_Fxyz_Px;
  abcd[201] = 4.0E0*I_NAI_Kx4y2z_Px_aa-2.0E0*2*I_NAI_Hx2y2z_Px_a-2.0E0*3*I_NAI_Hx2y2z_Px_a+2*1*I_NAI_Fx2z_Px;
  abcd[202] = 4.0E0*I_NAI_Kx3y3z_Px_aa-2.0E0*1*I_NAI_Hxy3z_Px_a-2.0E0*2*I_NAI_Hxy3z_Px_a;
  abcd[203] = 4.0E0*I_NAI_Kx2y4z_Px_aa-2.0E0*1*I_NAI_Hx4z_Px_a;
  abcd[204] = 4.0E0*I_NAI_K7y_Px_aa-2.0E0*5*I_NAI_H5y_Px_a-2.0E0*6*I_NAI_H5y_Px_a+5*4*I_NAI_F3y_Px;
  abcd[205] = 4.0E0*I_NAI_K6yz_Px_aa-2.0E0*4*I_NAI_H4yz_Px_a-2.0E0*5*I_NAI_H4yz_Px_a+4*3*I_NAI_F2yz_Px;
  abcd[206] = 4.0E0*I_NAI_K5y2z_Px_aa-2.0E0*3*I_NAI_H3y2z_Px_a-2.0E0*4*I_NAI_H3y2z_Px_a+3*2*I_NAI_Fy2z_Px;
  abcd[207] = 4.0E0*I_NAI_K4y3z_Px_aa-2.0E0*2*I_NAI_H2y3z_Px_a-2.0E0*3*I_NAI_H2y3z_Px_a+2*1*I_NAI_F3z_Px;
  abcd[208] = 4.0E0*I_NAI_K3y4z_Px_aa-2.0E0*1*I_NAI_Hy4z_Px_a-2.0E0*2*I_NAI_Hy4z_Px_a;
  abcd[209] = 4.0E0*I_NAI_K2y5z_Px_aa-2.0E0*1*I_NAI_H5z_Px_a;
  abcd[210] = 4.0E0*I_NAI_K5x2y_Py_aa-2.0E0*1*I_NAI_H5x_Py_a;
  abcd[211] = 4.0E0*I_NAI_K4x3y_Py_aa-2.0E0*1*I_NAI_H4xy_Py_a-2.0E0*2*I_NAI_H4xy_Py_a;
  abcd[212] = 4.0E0*I_NAI_K4x2yz_Py_aa-2.0E0*1*I_NAI_H4xz_Py_a;
  abcd[213] = 4.0E0*I_NAI_K3x4y_Py_aa-2.0E0*2*I_NAI_H3x2y_Py_a-2.0E0*3*I_NAI_H3x2y_Py_a+2*1*I_NAI_F3x_Py;
  abcd[214] = 4.0E0*I_NAI_K3x3yz_Py_aa-2.0E0*1*I_NAI_H3xyz_Py_a-2.0E0*2*I_NAI_H3xyz_Py_a;
  abcd[215] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_H3x2z_Py_a;
  abcd[216] = 4.0E0*I_NAI_K2x5y_Py_aa-2.0E0*3*I_NAI_H2x3y_Py_a-2.0E0*4*I_NAI_H2x3y_Py_a+3*2*I_NAI_F2xy_Py;
  abcd[217] = 4.0E0*I_NAI_K2x4yz_Py_aa-2.0E0*2*I_NAI_H2x2yz_Py_a-2.0E0*3*I_NAI_H2x2yz_Py_a+2*1*I_NAI_F2xz_Py;
  abcd[218] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*1*I_NAI_H2xy2z_Py_a-2.0E0*2*I_NAI_H2xy2z_Py_a;
  abcd[219] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*1*I_NAI_H2x3z_Py_a;
  abcd[220] = 4.0E0*I_NAI_Kx6y_Py_aa-2.0E0*4*I_NAI_Hx4y_Py_a-2.0E0*5*I_NAI_Hx4y_Py_a+4*3*I_NAI_Fx2y_Py;
  abcd[221] = 4.0E0*I_NAI_Kx5yz_Py_aa-2.0E0*3*I_NAI_Hx3yz_Py_a-2.0E0*4*I_NAI_Hx3yz_Py_a+3*2*I_NAI_Fxyz_Py;
  abcd[222] = 4.0E0*I_NAI_Kx4y2z_Py_aa-2.0E0*2*I_NAI_Hx2y2z_Py_a-2.0E0*3*I_NAI_Hx2y2z_Py_a+2*1*I_NAI_Fx2z_Py;
  abcd[223] = 4.0E0*I_NAI_Kx3y3z_Py_aa-2.0E0*1*I_NAI_Hxy3z_Py_a-2.0E0*2*I_NAI_Hxy3z_Py_a;
  abcd[224] = 4.0E0*I_NAI_Kx2y4z_Py_aa-2.0E0*1*I_NAI_Hx4z_Py_a;
  abcd[225] = 4.0E0*I_NAI_K7y_Py_aa-2.0E0*5*I_NAI_H5y_Py_a-2.0E0*6*I_NAI_H5y_Py_a+5*4*I_NAI_F3y_Py;
  abcd[226] = 4.0E0*I_NAI_K6yz_Py_aa-2.0E0*4*I_NAI_H4yz_Py_a-2.0E0*5*I_NAI_H4yz_Py_a+4*3*I_NAI_F2yz_Py;
  abcd[227] = 4.0E0*I_NAI_K5y2z_Py_aa-2.0E0*3*I_NAI_H3y2z_Py_a-2.0E0*4*I_NAI_H3y2z_Py_a+3*2*I_NAI_Fy2z_Py;
  abcd[228] = 4.0E0*I_NAI_K4y3z_Py_aa-2.0E0*2*I_NAI_H2y3z_Py_a-2.0E0*3*I_NAI_H2y3z_Py_a+2*1*I_NAI_F3z_Py;
  abcd[229] = 4.0E0*I_NAI_K3y4z_Py_aa-2.0E0*1*I_NAI_Hy4z_Py_a-2.0E0*2*I_NAI_Hy4z_Py_a;
  abcd[230] = 4.0E0*I_NAI_K2y5z_Py_aa-2.0E0*1*I_NAI_H5z_Py_a;
  abcd[231] = 4.0E0*I_NAI_K5x2y_Pz_aa-2.0E0*1*I_NAI_H5x_Pz_a;
  abcd[232] = 4.0E0*I_NAI_K4x3y_Pz_aa-2.0E0*1*I_NAI_H4xy_Pz_a-2.0E0*2*I_NAI_H4xy_Pz_a;
  abcd[233] = 4.0E0*I_NAI_K4x2yz_Pz_aa-2.0E0*1*I_NAI_H4xz_Pz_a;
  abcd[234] = 4.0E0*I_NAI_K3x4y_Pz_aa-2.0E0*2*I_NAI_H3x2y_Pz_a-2.0E0*3*I_NAI_H3x2y_Pz_a+2*1*I_NAI_F3x_Pz;
  abcd[235] = 4.0E0*I_NAI_K3x3yz_Pz_aa-2.0E0*1*I_NAI_H3xyz_Pz_a-2.0E0*2*I_NAI_H3xyz_Pz_a;
  abcd[236] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_H3x2z_Pz_a;
  abcd[237] = 4.0E0*I_NAI_K2x5y_Pz_aa-2.0E0*3*I_NAI_H2x3y_Pz_a-2.0E0*4*I_NAI_H2x3y_Pz_a+3*2*I_NAI_F2xy_Pz;
  abcd[238] = 4.0E0*I_NAI_K2x4yz_Pz_aa-2.0E0*2*I_NAI_H2x2yz_Pz_a-2.0E0*3*I_NAI_H2x2yz_Pz_a+2*1*I_NAI_F2xz_Pz;
  abcd[239] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*1*I_NAI_H2xy2z_Pz_a-2.0E0*2*I_NAI_H2xy2z_Pz_a;
  abcd[240] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*1*I_NAI_H2x3z_Pz_a;
  abcd[241] = 4.0E0*I_NAI_Kx6y_Pz_aa-2.0E0*4*I_NAI_Hx4y_Pz_a-2.0E0*5*I_NAI_Hx4y_Pz_a+4*3*I_NAI_Fx2y_Pz;
  abcd[242] = 4.0E0*I_NAI_Kx5yz_Pz_aa-2.0E0*3*I_NAI_Hx3yz_Pz_a-2.0E0*4*I_NAI_Hx3yz_Pz_a+3*2*I_NAI_Fxyz_Pz;
  abcd[243] = 4.0E0*I_NAI_Kx4y2z_Pz_aa-2.0E0*2*I_NAI_Hx2y2z_Pz_a-2.0E0*3*I_NAI_Hx2y2z_Pz_a+2*1*I_NAI_Fx2z_Pz;
  abcd[244] = 4.0E0*I_NAI_Kx3y3z_Pz_aa-2.0E0*1*I_NAI_Hxy3z_Pz_a-2.0E0*2*I_NAI_Hxy3z_Pz_a;
  abcd[245] = 4.0E0*I_NAI_Kx2y4z_Pz_aa-2.0E0*1*I_NAI_Hx4z_Pz_a;
  abcd[246] = 4.0E0*I_NAI_K7y_Pz_aa-2.0E0*5*I_NAI_H5y_Pz_a-2.0E0*6*I_NAI_H5y_Pz_a+5*4*I_NAI_F3y_Pz;
  abcd[247] = 4.0E0*I_NAI_K6yz_Pz_aa-2.0E0*4*I_NAI_H4yz_Pz_a-2.0E0*5*I_NAI_H4yz_Pz_a+4*3*I_NAI_F2yz_Pz;
  abcd[248] = 4.0E0*I_NAI_K5y2z_Pz_aa-2.0E0*3*I_NAI_H3y2z_Pz_a-2.0E0*4*I_NAI_H3y2z_Pz_a+3*2*I_NAI_Fy2z_Pz;
  abcd[249] = 4.0E0*I_NAI_K4y3z_Pz_aa-2.0E0*2*I_NAI_H2y3z_Pz_a-2.0E0*3*I_NAI_H2y3z_Pz_a+2*1*I_NAI_F3z_Pz;
  abcd[250] = 4.0E0*I_NAI_K3y4z_Pz_aa-2.0E0*1*I_NAI_Hy4z_Pz_a-2.0E0*2*I_NAI_Hy4z_Pz_a;
  abcd[251] = 4.0E0*I_NAI_K2y5z_Pz_aa-2.0E0*1*I_NAI_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[252] = 4.0E0*I_NAI_K5xyz_Px_aa;
  abcd[253] = 4.0E0*I_NAI_K4x2yz_Px_aa-2.0E0*1*I_NAI_H4xz_Px_a;
  abcd[254] = 4.0E0*I_NAI_K4xy2z_Px_aa-2.0E0*1*I_NAI_H4xy_Px_a;
  abcd[255] = 4.0E0*I_NAI_K3x3yz_Px_aa-2.0E0*2*I_NAI_H3xyz_Px_a;
  abcd[256] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_H3x2y_Px_a-2.0E0*1*I_NAI_H3x2z_Px_a+1*I_NAI_F3x_Px;
  abcd[257] = 4.0E0*I_NAI_K3xy3z_Px_aa-2.0E0*2*I_NAI_H3xyz_Px_a;
  abcd[258] = 4.0E0*I_NAI_K2x4yz_Px_aa-2.0E0*3*I_NAI_H2x2yz_Px_a;
  abcd[259] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*1*I_NAI_H2x3y_Px_a-2.0E0*2*I_NAI_H2xy2z_Px_a+2*1*I_NAI_F2xy_Px;
  abcd[260] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*2*I_NAI_H2x2yz_Px_a-2.0E0*1*I_NAI_H2x3z_Px_a+2*I_NAI_F2xz_Px;
  abcd[261] = 4.0E0*I_NAI_K2xy4z_Px_aa-2.0E0*3*I_NAI_H2xy2z_Px_a;
  abcd[262] = 4.0E0*I_NAI_Kx5yz_Px_aa-2.0E0*4*I_NAI_Hx3yz_Px_a;
  abcd[263] = 4.0E0*I_NAI_Kx4y2z_Px_aa-2.0E0*1*I_NAI_Hx4y_Px_a-2.0E0*3*I_NAI_Hx2y2z_Px_a+3*1*I_NAI_Fx2y_Px;
  abcd[264] = 4.0E0*I_NAI_Kx3y3z_Px_aa-2.0E0*2*I_NAI_Hx3yz_Px_a-2.0E0*2*I_NAI_Hxy3z_Px_a+2*2*I_NAI_Fxyz_Px;
  abcd[265] = 4.0E0*I_NAI_Kx2y4z_Px_aa-2.0E0*3*I_NAI_Hx2y2z_Px_a-2.0E0*1*I_NAI_Hx4z_Px_a+3*I_NAI_Fx2z_Px;
  abcd[266] = 4.0E0*I_NAI_Kxy5z_Px_aa-2.0E0*4*I_NAI_Hxy3z_Px_a;
  abcd[267] = 4.0E0*I_NAI_K6yz_Px_aa-2.0E0*5*I_NAI_H4yz_Px_a;
  abcd[268] = 4.0E0*I_NAI_K5y2z_Px_aa-2.0E0*1*I_NAI_H5y_Px_a-2.0E0*4*I_NAI_H3y2z_Px_a+4*1*I_NAI_F3y_Px;
  abcd[269] = 4.0E0*I_NAI_K4y3z_Px_aa-2.0E0*2*I_NAI_H4yz_Px_a-2.0E0*3*I_NAI_H2y3z_Px_a+3*2*I_NAI_F2yz_Px;
  abcd[270] = 4.0E0*I_NAI_K3y4z_Px_aa-2.0E0*3*I_NAI_H3y2z_Px_a-2.0E0*2*I_NAI_Hy4z_Px_a+2*3*I_NAI_Fy2z_Px;
  abcd[271] = 4.0E0*I_NAI_K2y5z_Px_aa-2.0E0*4*I_NAI_H2y3z_Px_a-2.0E0*1*I_NAI_H5z_Px_a+4*I_NAI_F3z_Px;
  abcd[272] = 4.0E0*I_NAI_Ky6z_Px_aa-2.0E0*5*I_NAI_Hy4z_Px_a;
  abcd[273] = 4.0E0*I_NAI_K5xyz_Py_aa;
  abcd[274] = 4.0E0*I_NAI_K4x2yz_Py_aa-2.0E0*1*I_NAI_H4xz_Py_a;
  abcd[275] = 4.0E0*I_NAI_K4xy2z_Py_aa-2.0E0*1*I_NAI_H4xy_Py_a;
  abcd[276] = 4.0E0*I_NAI_K3x3yz_Py_aa-2.0E0*2*I_NAI_H3xyz_Py_a;
  abcd[277] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_H3x2y_Py_a-2.0E0*1*I_NAI_H3x2z_Py_a+1*I_NAI_F3x_Py;
  abcd[278] = 4.0E0*I_NAI_K3xy3z_Py_aa-2.0E0*2*I_NAI_H3xyz_Py_a;
  abcd[279] = 4.0E0*I_NAI_K2x4yz_Py_aa-2.0E0*3*I_NAI_H2x2yz_Py_a;
  abcd[280] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*1*I_NAI_H2x3y_Py_a-2.0E0*2*I_NAI_H2xy2z_Py_a+2*1*I_NAI_F2xy_Py;
  abcd[281] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*2*I_NAI_H2x2yz_Py_a-2.0E0*1*I_NAI_H2x3z_Py_a+2*I_NAI_F2xz_Py;
  abcd[282] = 4.0E0*I_NAI_K2xy4z_Py_aa-2.0E0*3*I_NAI_H2xy2z_Py_a;
  abcd[283] = 4.0E0*I_NAI_Kx5yz_Py_aa-2.0E0*4*I_NAI_Hx3yz_Py_a;
  abcd[284] = 4.0E0*I_NAI_Kx4y2z_Py_aa-2.0E0*1*I_NAI_Hx4y_Py_a-2.0E0*3*I_NAI_Hx2y2z_Py_a+3*1*I_NAI_Fx2y_Py;
  abcd[285] = 4.0E0*I_NAI_Kx3y3z_Py_aa-2.0E0*2*I_NAI_Hx3yz_Py_a-2.0E0*2*I_NAI_Hxy3z_Py_a+2*2*I_NAI_Fxyz_Py;
  abcd[286] = 4.0E0*I_NAI_Kx2y4z_Py_aa-2.0E0*3*I_NAI_Hx2y2z_Py_a-2.0E0*1*I_NAI_Hx4z_Py_a+3*I_NAI_Fx2z_Py;
  abcd[287] = 4.0E0*I_NAI_Kxy5z_Py_aa-2.0E0*4*I_NAI_Hxy3z_Py_a;
  abcd[288] = 4.0E0*I_NAI_K6yz_Py_aa-2.0E0*5*I_NAI_H4yz_Py_a;
  abcd[289] = 4.0E0*I_NAI_K5y2z_Py_aa-2.0E0*1*I_NAI_H5y_Py_a-2.0E0*4*I_NAI_H3y2z_Py_a+4*1*I_NAI_F3y_Py;
  abcd[290] = 4.0E0*I_NAI_K4y3z_Py_aa-2.0E0*2*I_NAI_H4yz_Py_a-2.0E0*3*I_NAI_H2y3z_Py_a+3*2*I_NAI_F2yz_Py;
  abcd[291] = 4.0E0*I_NAI_K3y4z_Py_aa-2.0E0*3*I_NAI_H3y2z_Py_a-2.0E0*2*I_NAI_Hy4z_Py_a+2*3*I_NAI_Fy2z_Py;
  abcd[292] = 4.0E0*I_NAI_K2y5z_Py_aa-2.0E0*4*I_NAI_H2y3z_Py_a-2.0E0*1*I_NAI_H5z_Py_a+4*I_NAI_F3z_Py;
  abcd[293] = 4.0E0*I_NAI_Ky6z_Py_aa-2.0E0*5*I_NAI_Hy4z_Py_a;
  abcd[294] = 4.0E0*I_NAI_K5xyz_Pz_aa;
  abcd[295] = 4.0E0*I_NAI_K4x2yz_Pz_aa-2.0E0*1*I_NAI_H4xz_Pz_a;
  abcd[296] = 4.0E0*I_NAI_K4xy2z_Pz_aa-2.0E0*1*I_NAI_H4xy_Pz_a;
  abcd[297] = 4.0E0*I_NAI_K3x3yz_Pz_aa-2.0E0*2*I_NAI_H3xyz_Pz_a;
  abcd[298] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_H3x2y_Pz_a-2.0E0*1*I_NAI_H3x2z_Pz_a+1*I_NAI_F3x_Pz;
  abcd[299] = 4.0E0*I_NAI_K3xy3z_Pz_aa-2.0E0*2*I_NAI_H3xyz_Pz_a;
  abcd[300] = 4.0E0*I_NAI_K2x4yz_Pz_aa-2.0E0*3*I_NAI_H2x2yz_Pz_a;
  abcd[301] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*1*I_NAI_H2x3y_Pz_a-2.0E0*2*I_NAI_H2xy2z_Pz_a+2*1*I_NAI_F2xy_Pz;
  abcd[302] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*2*I_NAI_H2x2yz_Pz_a-2.0E0*1*I_NAI_H2x3z_Pz_a+2*I_NAI_F2xz_Pz;
  abcd[303] = 4.0E0*I_NAI_K2xy4z_Pz_aa-2.0E0*3*I_NAI_H2xy2z_Pz_a;
  abcd[304] = 4.0E0*I_NAI_Kx5yz_Pz_aa-2.0E0*4*I_NAI_Hx3yz_Pz_a;
  abcd[305] = 4.0E0*I_NAI_Kx4y2z_Pz_aa-2.0E0*1*I_NAI_Hx4y_Pz_a-2.0E0*3*I_NAI_Hx2y2z_Pz_a+3*1*I_NAI_Fx2y_Pz;
  abcd[306] = 4.0E0*I_NAI_Kx3y3z_Pz_aa-2.0E0*2*I_NAI_Hx3yz_Pz_a-2.0E0*2*I_NAI_Hxy3z_Pz_a+2*2*I_NAI_Fxyz_Pz;
  abcd[307] = 4.0E0*I_NAI_Kx2y4z_Pz_aa-2.0E0*3*I_NAI_Hx2y2z_Pz_a-2.0E0*1*I_NAI_Hx4z_Pz_a+3*I_NAI_Fx2z_Pz;
  abcd[308] = 4.0E0*I_NAI_Kxy5z_Pz_aa-2.0E0*4*I_NAI_Hxy3z_Pz_a;
  abcd[309] = 4.0E0*I_NAI_K6yz_Pz_aa-2.0E0*5*I_NAI_H4yz_Pz_a;
  abcd[310] = 4.0E0*I_NAI_K5y2z_Pz_aa-2.0E0*1*I_NAI_H5y_Pz_a-2.0E0*4*I_NAI_H3y2z_Pz_a+4*1*I_NAI_F3y_Pz;
  abcd[311] = 4.0E0*I_NAI_K4y3z_Pz_aa-2.0E0*2*I_NAI_H4yz_Pz_a-2.0E0*3*I_NAI_H2y3z_Pz_a+3*2*I_NAI_F2yz_Pz;
  abcd[312] = 4.0E0*I_NAI_K3y4z_Pz_aa-2.0E0*3*I_NAI_H3y2z_Pz_a-2.0E0*2*I_NAI_Hy4z_Pz_a+2*3*I_NAI_Fy2z_Pz;
  abcd[313] = 4.0E0*I_NAI_K2y5z_Pz_aa-2.0E0*4*I_NAI_H2y3z_Pz_a-2.0E0*1*I_NAI_H5z_Pz_a+4*I_NAI_F3z_Pz;
  abcd[314] = 4.0E0*I_NAI_Ky6z_Pz_aa-2.0E0*5*I_NAI_Hy4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_aa
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  abcd[315] = 4.0E0*I_NAI_K5x2z_Px_aa-2.0E0*1*I_NAI_H5x_Px_a;
  abcd[316] = 4.0E0*I_NAI_K4xy2z_Px_aa-2.0E0*1*I_NAI_H4xy_Px_a;
  abcd[317] = 4.0E0*I_NAI_K4x3z_Px_aa-2.0E0*1*I_NAI_H4xz_Px_a-2.0E0*2*I_NAI_H4xz_Px_a;
  abcd[318] = 4.0E0*I_NAI_K3x2y2z_Px_aa-2.0E0*1*I_NAI_H3x2y_Px_a;
  abcd[319] = 4.0E0*I_NAI_K3xy3z_Px_aa-2.0E0*1*I_NAI_H3xyz_Px_a-2.0E0*2*I_NAI_H3xyz_Px_a;
  abcd[320] = 4.0E0*I_NAI_K3x4z_Px_aa-2.0E0*2*I_NAI_H3x2z_Px_a-2.0E0*3*I_NAI_H3x2z_Px_a+2*1*I_NAI_F3x_Px;
  abcd[321] = 4.0E0*I_NAI_K2x3y2z_Px_aa-2.0E0*1*I_NAI_H2x3y_Px_a;
  abcd[322] = 4.0E0*I_NAI_K2x2y3z_Px_aa-2.0E0*1*I_NAI_H2x2yz_Px_a-2.0E0*2*I_NAI_H2x2yz_Px_a;
  abcd[323] = 4.0E0*I_NAI_K2xy4z_Px_aa-2.0E0*2*I_NAI_H2xy2z_Px_a-2.0E0*3*I_NAI_H2xy2z_Px_a+2*1*I_NAI_F2xy_Px;
  abcd[324] = 4.0E0*I_NAI_K2x5z_Px_aa-2.0E0*3*I_NAI_H2x3z_Px_a-2.0E0*4*I_NAI_H2x3z_Px_a+3*2*I_NAI_F2xz_Px;
  abcd[325] = 4.0E0*I_NAI_Kx4y2z_Px_aa-2.0E0*1*I_NAI_Hx4y_Px_a;
  abcd[326] = 4.0E0*I_NAI_Kx3y3z_Px_aa-2.0E0*1*I_NAI_Hx3yz_Px_a-2.0E0*2*I_NAI_Hx3yz_Px_a;
  abcd[327] = 4.0E0*I_NAI_Kx2y4z_Px_aa-2.0E0*2*I_NAI_Hx2y2z_Px_a-2.0E0*3*I_NAI_Hx2y2z_Px_a+2*1*I_NAI_Fx2y_Px;
  abcd[328] = 4.0E0*I_NAI_Kxy5z_Px_aa-2.0E0*3*I_NAI_Hxy3z_Px_a-2.0E0*4*I_NAI_Hxy3z_Px_a+3*2*I_NAI_Fxyz_Px;
  abcd[329] = 4.0E0*I_NAI_Kx6z_Px_aa-2.0E0*4*I_NAI_Hx4z_Px_a-2.0E0*5*I_NAI_Hx4z_Px_a+4*3*I_NAI_Fx2z_Px;
  abcd[330] = 4.0E0*I_NAI_K5y2z_Px_aa-2.0E0*1*I_NAI_H5y_Px_a;
  abcd[331] = 4.0E0*I_NAI_K4y3z_Px_aa-2.0E0*1*I_NAI_H4yz_Px_a-2.0E0*2*I_NAI_H4yz_Px_a;
  abcd[332] = 4.0E0*I_NAI_K3y4z_Px_aa-2.0E0*2*I_NAI_H3y2z_Px_a-2.0E0*3*I_NAI_H3y2z_Px_a+2*1*I_NAI_F3y_Px;
  abcd[333] = 4.0E0*I_NAI_K2y5z_Px_aa-2.0E0*3*I_NAI_H2y3z_Px_a-2.0E0*4*I_NAI_H2y3z_Px_a+3*2*I_NAI_F2yz_Px;
  abcd[334] = 4.0E0*I_NAI_Ky6z_Px_aa-2.0E0*4*I_NAI_Hy4z_Px_a-2.0E0*5*I_NAI_Hy4z_Px_a+4*3*I_NAI_Fy2z_Px;
  abcd[335] = 4.0E0*I_NAI_K7z_Px_aa-2.0E0*5*I_NAI_H5z_Px_a-2.0E0*6*I_NAI_H5z_Px_a+5*4*I_NAI_F3z_Px;
  abcd[336] = 4.0E0*I_NAI_K5x2z_Py_aa-2.0E0*1*I_NAI_H5x_Py_a;
  abcd[337] = 4.0E0*I_NAI_K4xy2z_Py_aa-2.0E0*1*I_NAI_H4xy_Py_a;
  abcd[338] = 4.0E0*I_NAI_K4x3z_Py_aa-2.0E0*1*I_NAI_H4xz_Py_a-2.0E0*2*I_NAI_H4xz_Py_a;
  abcd[339] = 4.0E0*I_NAI_K3x2y2z_Py_aa-2.0E0*1*I_NAI_H3x2y_Py_a;
  abcd[340] = 4.0E0*I_NAI_K3xy3z_Py_aa-2.0E0*1*I_NAI_H3xyz_Py_a-2.0E0*2*I_NAI_H3xyz_Py_a;
  abcd[341] = 4.0E0*I_NAI_K3x4z_Py_aa-2.0E0*2*I_NAI_H3x2z_Py_a-2.0E0*3*I_NAI_H3x2z_Py_a+2*1*I_NAI_F3x_Py;
  abcd[342] = 4.0E0*I_NAI_K2x3y2z_Py_aa-2.0E0*1*I_NAI_H2x3y_Py_a;
  abcd[343] = 4.0E0*I_NAI_K2x2y3z_Py_aa-2.0E0*1*I_NAI_H2x2yz_Py_a-2.0E0*2*I_NAI_H2x2yz_Py_a;
  abcd[344] = 4.0E0*I_NAI_K2xy4z_Py_aa-2.0E0*2*I_NAI_H2xy2z_Py_a-2.0E0*3*I_NAI_H2xy2z_Py_a+2*1*I_NAI_F2xy_Py;
  abcd[345] = 4.0E0*I_NAI_K2x5z_Py_aa-2.0E0*3*I_NAI_H2x3z_Py_a-2.0E0*4*I_NAI_H2x3z_Py_a+3*2*I_NAI_F2xz_Py;
  abcd[346] = 4.0E0*I_NAI_Kx4y2z_Py_aa-2.0E0*1*I_NAI_Hx4y_Py_a;
  abcd[347] = 4.0E0*I_NAI_Kx3y3z_Py_aa-2.0E0*1*I_NAI_Hx3yz_Py_a-2.0E0*2*I_NAI_Hx3yz_Py_a;
  abcd[348] = 4.0E0*I_NAI_Kx2y4z_Py_aa-2.0E0*2*I_NAI_Hx2y2z_Py_a-2.0E0*3*I_NAI_Hx2y2z_Py_a+2*1*I_NAI_Fx2y_Py;
  abcd[349] = 4.0E0*I_NAI_Kxy5z_Py_aa-2.0E0*3*I_NAI_Hxy3z_Py_a-2.0E0*4*I_NAI_Hxy3z_Py_a+3*2*I_NAI_Fxyz_Py;
  abcd[350] = 4.0E0*I_NAI_Kx6z_Py_aa-2.0E0*4*I_NAI_Hx4z_Py_a-2.0E0*5*I_NAI_Hx4z_Py_a+4*3*I_NAI_Fx2z_Py;
  abcd[351] = 4.0E0*I_NAI_K5y2z_Py_aa-2.0E0*1*I_NAI_H5y_Py_a;
  abcd[352] = 4.0E0*I_NAI_K4y3z_Py_aa-2.0E0*1*I_NAI_H4yz_Py_a-2.0E0*2*I_NAI_H4yz_Py_a;
  abcd[353] = 4.0E0*I_NAI_K3y4z_Py_aa-2.0E0*2*I_NAI_H3y2z_Py_a-2.0E0*3*I_NAI_H3y2z_Py_a+2*1*I_NAI_F3y_Py;
  abcd[354] = 4.0E0*I_NAI_K2y5z_Py_aa-2.0E0*3*I_NAI_H2y3z_Py_a-2.0E0*4*I_NAI_H2y3z_Py_a+3*2*I_NAI_F2yz_Py;
  abcd[355] = 4.0E0*I_NAI_Ky6z_Py_aa-2.0E0*4*I_NAI_Hy4z_Py_a-2.0E0*5*I_NAI_Hy4z_Py_a+4*3*I_NAI_Fy2z_Py;
  abcd[356] = 4.0E0*I_NAI_K7z_Py_aa-2.0E0*5*I_NAI_H5z_Py_a-2.0E0*6*I_NAI_H5z_Py_a+5*4*I_NAI_F3z_Py;
  abcd[357] = 4.0E0*I_NAI_K5x2z_Pz_aa-2.0E0*1*I_NAI_H5x_Pz_a;
  abcd[358] = 4.0E0*I_NAI_K4xy2z_Pz_aa-2.0E0*1*I_NAI_H4xy_Pz_a;
  abcd[359] = 4.0E0*I_NAI_K4x3z_Pz_aa-2.0E0*1*I_NAI_H4xz_Pz_a-2.0E0*2*I_NAI_H4xz_Pz_a;
  abcd[360] = 4.0E0*I_NAI_K3x2y2z_Pz_aa-2.0E0*1*I_NAI_H3x2y_Pz_a;
  abcd[361] = 4.0E0*I_NAI_K3xy3z_Pz_aa-2.0E0*1*I_NAI_H3xyz_Pz_a-2.0E0*2*I_NAI_H3xyz_Pz_a;
  abcd[362] = 4.0E0*I_NAI_K3x4z_Pz_aa-2.0E0*2*I_NAI_H3x2z_Pz_a-2.0E0*3*I_NAI_H3x2z_Pz_a+2*1*I_NAI_F3x_Pz;
  abcd[363] = 4.0E0*I_NAI_K2x3y2z_Pz_aa-2.0E0*1*I_NAI_H2x3y_Pz_a;
  abcd[364] = 4.0E0*I_NAI_K2x2y3z_Pz_aa-2.0E0*1*I_NAI_H2x2yz_Pz_a-2.0E0*2*I_NAI_H2x2yz_Pz_a;
  abcd[365] = 4.0E0*I_NAI_K2xy4z_Pz_aa-2.0E0*2*I_NAI_H2xy2z_Pz_a-2.0E0*3*I_NAI_H2xy2z_Pz_a+2*1*I_NAI_F2xy_Pz;
  abcd[366] = 4.0E0*I_NAI_K2x5z_Pz_aa-2.0E0*3*I_NAI_H2x3z_Pz_a-2.0E0*4*I_NAI_H2x3z_Pz_a+3*2*I_NAI_F2xz_Pz;
  abcd[367] = 4.0E0*I_NAI_Kx4y2z_Pz_aa-2.0E0*1*I_NAI_Hx4y_Pz_a;
  abcd[368] = 4.0E0*I_NAI_Kx3y3z_Pz_aa-2.0E0*1*I_NAI_Hx3yz_Pz_a-2.0E0*2*I_NAI_Hx3yz_Pz_a;
  abcd[369] = 4.0E0*I_NAI_Kx2y4z_Pz_aa-2.0E0*2*I_NAI_Hx2y2z_Pz_a-2.0E0*3*I_NAI_Hx2y2z_Pz_a+2*1*I_NAI_Fx2y_Pz;
  abcd[370] = 4.0E0*I_NAI_Kxy5z_Pz_aa-2.0E0*3*I_NAI_Hxy3z_Pz_a-2.0E0*4*I_NAI_Hxy3z_Pz_a+3*2*I_NAI_Fxyz_Pz;
  abcd[371] = 4.0E0*I_NAI_Kx6z_Pz_aa-2.0E0*4*I_NAI_Hx4z_Pz_a-2.0E0*5*I_NAI_Hx4z_Pz_a+4*3*I_NAI_Fx2z_Pz;
  abcd[372] = 4.0E0*I_NAI_K5y2z_Pz_aa-2.0E0*1*I_NAI_H5y_Pz_a;
  abcd[373] = 4.0E0*I_NAI_K4y3z_Pz_aa-2.0E0*1*I_NAI_H4yz_Pz_a-2.0E0*2*I_NAI_H4yz_Pz_a;
  abcd[374] = 4.0E0*I_NAI_K3y4z_Pz_aa-2.0E0*2*I_NAI_H3y2z_Pz_a-2.0E0*3*I_NAI_H3y2z_Pz_a+2*1*I_NAI_F3y_Pz;
  abcd[375] = 4.0E0*I_NAI_K2y5z_Pz_aa-2.0E0*3*I_NAI_H2y3z_Pz_a-2.0E0*4*I_NAI_H2y3z_Pz_a+3*2*I_NAI_F2yz_Pz;
  abcd[376] = 4.0E0*I_NAI_Ky6z_Pz_aa-2.0E0*4*I_NAI_Hy4z_Pz_a-2.0E0*5*I_NAI_Hy4z_Pz_a+4*3*I_NAI_Fy2z_Pz;
  abcd[377] = 4.0E0*I_NAI_K7z_Pz_aa-2.0E0*5*I_NAI_H5z_Pz_a-2.0E0*6*I_NAI_H5z_Pz_a+5*4*I_NAI_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[378] = 4.0E0*I_NAI_I6x_D2x_ab-2.0E0*1*I_NAI_I6x_S_a-2.0E0*5*I_NAI_G4x_D2x_b+5*1*I_NAI_G4x_S;
  abcd[379] = 4.0E0*I_NAI_I5xy_D2x_ab-2.0E0*1*I_NAI_I5xy_S_a-2.0E0*4*I_NAI_G3xy_D2x_b+4*1*I_NAI_G3xy_S;
  abcd[380] = 4.0E0*I_NAI_I5xz_D2x_ab-2.0E0*1*I_NAI_I5xz_S_a-2.0E0*4*I_NAI_G3xz_D2x_b+4*1*I_NAI_G3xz_S;
  abcd[381] = 4.0E0*I_NAI_I4x2y_D2x_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*3*I_NAI_G2x2y_D2x_b+3*1*I_NAI_G2x2y_S;
  abcd[382] = 4.0E0*I_NAI_I4xyz_D2x_ab-2.0E0*1*I_NAI_I4xyz_S_a-2.0E0*3*I_NAI_G2xyz_D2x_b+3*1*I_NAI_G2xyz_S;
  abcd[383] = 4.0E0*I_NAI_I4x2z_D2x_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*3*I_NAI_G2x2z_D2x_b+3*1*I_NAI_G2x2z_S;
  abcd[384] = 4.0E0*I_NAI_I3x3y_D2x_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_Gx3y_D2x_b+2*1*I_NAI_Gx3y_S;
  abcd[385] = 4.0E0*I_NAI_I3x2yz_D2x_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*2*I_NAI_Gx2yz_D2x_b+2*1*I_NAI_Gx2yz_S;
  abcd[386] = 4.0E0*I_NAI_I3xy2z_D2x_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2x_b+2*1*I_NAI_Gxy2z_S;
  abcd[387] = 4.0E0*I_NAI_I3x3z_D2x_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_Gx3z_D2x_b+2*1*I_NAI_Gx3z_S;
  abcd[388] = 4.0E0*I_NAI_I2x4y_D2x_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*1*I_NAI_G4y_D2x_b+1*I_NAI_G4y_S;
  abcd[389] = 4.0E0*I_NAI_I2x3yz_D2x_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*1*I_NAI_G3yz_D2x_b+1*I_NAI_G3yz_S;
  abcd[390] = 4.0E0*I_NAI_I2x2y2z_D2x_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2y2z_D2x_b+1*I_NAI_G2y2z_S;
  abcd[391] = 4.0E0*I_NAI_I2xy3z_D2x_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*1*I_NAI_Gy3z_D2x_b+1*I_NAI_Gy3z_S;
  abcd[392] = 4.0E0*I_NAI_I2x4z_D2x_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*1*I_NAI_G4z_D2x_b+1*I_NAI_G4z_S;
  abcd[393] = 4.0E0*I_NAI_Ix5y_D2x_ab-2.0E0*1*I_NAI_Ix5y_S_a;
  abcd[394] = 4.0E0*I_NAI_Ix4yz_D2x_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[395] = 4.0E0*I_NAI_Ix3y2z_D2x_ab-2.0E0*1*I_NAI_Ix3y2z_S_a;
  abcd[396] = 4.0E0*I_NAI_Ix2y3z_D2x_ab-2.0E0*1*I_NAI_Ix2y3z_S_a;
  abcd[397] = 4.0E0*I_NAI_Ixy4z_D2x_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[398] = 4.0E0*I_NAI_Ix5z_D2x_ab-2.0E0*1*I_NAI_Ix5z_S_a;
  abcd[399] = 4.0E0*I_NAI_I6x_Dxy_ab-2.0E0*5*I_NAI_G4x_Dxy_b;
  abcd[400] = 4.0E0*I_NAI_I5xy_Dxy_ab-2.0E0*4*I_NAI_G3xy_Dxy_b;
  abcd[401] = 4.0E0*I_NAI_I5xz_Dxy_ab-2.0E0*4*I_NAI_G3xz_Dxy_b;
  abcd[402] = 4.0E0*I_NAI_I4x2y_Dxy_ab-2.0E0*3*I_NAI_G2x2y_Dxy_b;
  abcd[403] = 4.0E0*I_NAI_I4xyz_Dxy_ab-2.0E0*3*I_NAI_G2xyz_Dxy_b;
  abcd[404] = 4.0E0*I_NAI_I4x2z_Dxy_ab-2.0E0*3*I_NAI_G2x2z_Dxy_b;
  abcd[405] = 4.0E0*I_NAI_I3x3y_Dxy_ab-2.0E0*2*I_NAI_Gx3y_Dxy_b;
  abcd[406] = 4.0E0*I_NAI_I3x2yz_Dxy_ab-2.0E0*2*I_NAI_Gx2yz_Dxy_b;
  abcd[407] = 4.0E0*I_NAI_I3xy2z_Dxy_ab-2.0E0*2*I_NAI_Gxy2z_Dxy_b;
  abcd[408] = 4.0E0*I_NAI_I3x3z_Dxy_ab-2.0E0*2*I_NAI_Gx3z_Dxy_b;
  abcd[409] = 4.0E0*I_NAI_I2x4y_Dxy_ab-2.0E0*1*I_NAI_G4y_Dxy_b;
  abcd[410] = 4.0E0*I_NAI_I2x3yz_Dxy_ab-2.0E0*1*I_NAI_G3yz_Dxy_b;
  abcd[411] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2y2z_Dxy_b;
  abcd[412] = 4.0E0*I_NAI_I2xy3z_Dxy_ab-2.0E0*1*I_NAI_Gy3z_Dxy_b;
  abcd[413] = 4.0E0*I_NAI_I2x4z_Dxy_ab-2.0E0*1*I_NAI_G4z_Dxy_b;
  abcd[414] = 4.0E0*I_NAI_Ix5y_Dxy_ab;
  abcd[415] = 4.0E0*I_NAI_Ix4yz_Dxy_ab;
  abcd[416] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab;
  abcd[417] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab;
  abcd[418] = 4.0E0*I_NAI_Ixy4z_Dxy_ab;
  abcd[419] = 4.0E0*I_NAI_Ix5z_Dxy_ab;
  abcd[420] = 4.0E0*I_NAI_I6x_Dxz_ab-2.0E0*5*I_NAI_G4x_Dxz_b;
  abcd[421] = 4.0E0*I_NAI_I5xy_Dxz_ab-2.0E0*4*I_NAI_G3xy_Dxz_b;
  abcd[422] = 4.0E0*I_NAI_I5xz_Dxz_ab-2.0E0*4*I_NAI_G3xz_Dxz_b;
  abcd[423] = 4.0E0*I_NAI_I4x2y_Dxz_ab-2.0E0*3*I_NAI_G2x2y_Dxz_b;
  abcd[424] = 4.0E0*I_NAI_I4xyz_Dxz_ab-2.0E0*3*I_NAI_G2xyz_Dxz_b;
  abcd[425] = 4.0E0*I_NAI_I4x2z_Dxz_ab-2.0E0*3*I_NAI_G2x2z_Dxz_b;
  abcd[426] = 4.0E0*I_NAI_I3x3y_Dxz_ab-2.0E0*2*I_NAI_Gx3y_Dxz_b;
  abcd[427] = 4.0E0*I_NAI_I3x2yz_Dxz_ab-2.0E0*2*I_NAI_Gx2yz_Dxz_b;
  abcd[428] = 4.0E0*I_NAI_I3xy2z_Dxz_ab-2.0E0*2*I_NAI_Gxy2z_Dxz_b;
  abcd[429] = 4.0E0*I_NAI_I3x3z_Dxz_ab-2.0E0*2*I_NAI_Gx3z_Dxz_b;
  abcd[430] = 4.0E0*I_NAI_I2x4y_Dxz_ab-2.0E0*1*I_NAI_G4y_Dxz_b;
  abcd[431] = 4.0E0*I_NAI_I2x3yz_Dxz_ab-2.0E0*1*I_NAI_G3yz_Dxz_b;
  abcd[432] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2y2z_Dxz_b;
  abcd[433] = 4.0E0*I_NAI_I2xy3z_Dxz_ab-2.0E0*1*I_NAI_Gy3z_Dxz_b;
  abcd[434] = 4.0E0*I_NAI_I2x4z_Dxz_ab-2.0E0*1*I_NAI_G4z_Dxz_b;
  abcd[435] = 4.0E0*I_NAI_Ix5y_Dxz_ab;
  abcd[436] = 4.0E0*I_NAI_Ix4yz_Dxz_ab;
  abcd[437] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab;
  abcd[438] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab;
  abcd[439] = 4.0E0*I_NAI_Ixy4z_Dxz_ab;
  abcd[440] = 4.0E0*I_NAI_Ix5z_Dxz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[441] = 4.0E0*I_NAI_I6x_Dxy_ab-2.0E0*5*I_NAI_G4x_Dxy_b;
  abcd[442] = 4.0E0*I_NAI_I5xy_Dxy_ab-2.0E0*4*I_NAI_G3xy_Dxy_b;
  abcd[443] = 4.0E0*I_NAI_I5xz_Dxy_ab-2.0E0*4*I_NAI_G3xz_Dxy_b;
  abcd[444] = 4.0E0*I_NAI_I4x2y_Dxy_ab-2.0E0*3*I_NAI_G2x2y_Dxy_b;
  abcd[445] = 4.0E0*I_NAI_I4xyz_Dxy_ab-2.0E0*3*I_NAI_G2xyz_Dxy_b;
  abcd[446] = 4.0E0*I_NAI_I4x2z_Dxy_ab-2.0E0*3*I_NAI_G2x2z_Dxy_b;
  abcd[447] = 4.0E0*I_NAI_I3x3y_Dxy_ab-2.0E0*2*I_NAI_Gx3y_Dxy_b;
  abcd[448] = 4.0E0*I_NAI_I3x2yz_Dxy_ab-2.0E0*2*I_NAI_Gx2yz_Dxy_b;
  abcd[449] = 4.0E0*I_NAI_I3xy2z_Dxy_ab-2.0E0*2*I_NAI_Gxy2z_Dxy_b;
  abcd[450] = 4.0E0*I_NAI_I3x3z_Dxy_ab-2.0E0*2*I_NAI_Gx3z_Dxy_b;
  abcd[451] = 4.0E0*I_NAI_I2x4y_Dxy_ab-2.0E0*1*I_NAI_G4y_Dxy_b;
  abcd[452] = 4.0E0*I_NAI_I2x3yz_Dxy_ab-2.0E0*1*I_NAI_G3yz_Dxy_b;
  abcd[453] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2y2z_Dxy_b;
  abcd[454] = 4.0E0*I_NAI_I2xy3z_Dxy_ab-2.0E0*1*I_NAI_Gy3z_Dxy_b;
  abcd[455] = 4.0E0*I_NAI_I2x4z_Dxy_ab-2.0E0*1*I_NAI_G4z_Dxy_b;
  abcd[456] = 4.0E0*I_NAI_Ix5y_Dxy_ab;
  abcd[457] = 4.0E0*I_NAI_Ix4yz_Dxy_ab;
  abcd[458] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab;
  abcd[459] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab;
  abcd[460] = 4.0E0*I_NAI_Ixy4z_Dxy_ab;
  abcd[461] = 4.0E0*I_NAI_Ix5z_Dxy_ab;
  abcd[462] = 4.0E0*I_NAI_I6x_D2y_ab-2.0E0*1*I_NAI_I6x_S_a-2.0E0*5*I_NAI_G4x_D2y_b+5*1*I_NAI_G4x_S;
  abcd[463] = 4.0E0*I_NAI_I5xy_D2y_ab-2.0E0*1*I_NAI_I5xy_S_a-2.0E0*4*I_NAI_G3xy_D2y_b+4*1*I_NAI_G3xy_S;
  abcd[464] = 4.0E0*I_NAI_I5xz_D2y_ab-2.0E0*1*I_NAI_I5xz_S_a-2.0E0*4*I_NAI_G3xz_D2y_b+4*1*I_NAI_G3xz_S;
  abcd[465] = 4.0E0*I_NAI_I4x2y_D2y_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*3*I_NAI_G2x2y_D2y_b+3*1*I_NAI_G2x2y_S;
  abcd[466] = 4.0E0*I_NAI_I4xyz_D2y_ab-2.0E0*1*I_NAI_I4xyz_S_a-2.0E0*3*I_NAI_G2xyz_D2y_b+3*1*I_NAI_G2xyz_S;
  abcd[467] = 4.0E0*I_NAI_I4x2z_D2y_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*3*I_NAI_G2x2z_D2y_b+3*1*I_NAI_G2x2z_S;
  abcd[468] = 4.0E0*I_NAI_I3x3y_D2y_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_Gx3y_D2y_b+2*1*I_NAI_Gx3y_S;
  abcd[469] = 4.0E0*I_NAI_I3x2yz_D2y_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*2*I_NAI_Gx2yz_D2y_b+2*1*I_NAI_Gx2yz_S;
  abcd[470] = 4.0E0*I_NAI_I3xy2z_D2y_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2y_b+2*1*I_NAI_Gxy2z_S;
  abcd[471] = 4.0E0*I_NAI_I3x3z_D2y_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_Gx3z_D2y_b+2*1*I_NAI_Gx3z_S;
  abcd[472] = 4.0E0*I_NAI_I2x4y_D2y_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*1*I_NAI_G4y_D2y_b+1*I_NAI_G4y_S;
  abcd[473] = 4.0E0*I_NAI_I2x3yz_D2y_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*1*I_NAI_G3yz_D2y_b+1*I_NAI_G3yz_S;
  abcd[474] = 4.0E0*I_NAI_I2x2y2z_D2y_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2y2z_D2y_b+1*I_NAI_G2y2z_S;
  abcd[475] = 4.0E0*I_NAI_I2xy3z_D2y_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*1*I_NAI_Gy3z_D2y_b+1*I_NAI_Gy3z_S;
  abcd[476] = 4.0E0*I_NAI_I2x4z_D2y_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*1*I_NAI_G4z_D2y_b+1*I_NAI_G4z_S;
  abcd[477] = 4.0E0*I_NAI_Ix5y_D2y_ab-2.0E0*1*I_NAI_Ix5y_S_a;
  abcd[478] = 4.0E0*I_NAI_Ix4yz_D2y_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[479] = 4.0E0*I_NAI_Ix3y2z_D2y_ab-2.0E0*1*I_NAI_Ix3y2z_S_a;
  abcd[480] = 4.0E0*I_NAI_Ix2y3z_D2y_ab-2.0E0*1*I_NAI_Ix2y3z_S_a;
  abcd[481] = 4.0E0*I_NAI_Ixy4z_D2y_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[482] = 4.0E0*I_NAI_Ix5z_D2y_ab-2.0E0*1*I_NAI_Ix5z_S_a;
  abcd[483] = 4.0E0*I_NAI_I6x_Dyz_ab-2.0E0*5*I_NAI_G4x_Dyz_b;
  abcd[484] = 4.0E0*I_NAI_I5xy_Dyz_ab-2.0E0*4*I_NAI_G3xy_Dyz_b;
  abcd[485] = 4.0E0*I_NAI_I5xz_Dyz_ab-2.0E0*4*I_NAI_G3xz_Dyz_b;
  abcd[486] = 4.0E0*I_NAI_I4x2y_Dyz_ab-2.0E0*3*I_NAI_G2x2y_Dyz_b;
  abcd[487] = 4.0E0*I_NAI_I4xyz_Dyz_ab-2.0E0*3*I_NAI_G2xyz_Dyz_b;
  abcd[488] = 4.0E0*I_NAI_I4x2z_Dyz_ab-2.0E0*3*I_NAI_G2x2z_Dyz_b;
  abcd[489] = 4.0E0*I_NAI_I3x3y_Dyz_ab-2.0E0*2*I_NAI_Gx3y_Dyz_b;
  abcd[490] = 4.0E0*I_NAI_I3x2yz_Dyz_ab-2.0E0*2*I_NAI_Gx2yz_Dyz_b;
  abcd[491] = 4.0E0*I_NAI_I3xy2z_Dyz_ab-2.0E0*2*I_NAI_Gxy2z_Dyz_b;
  abcd[492] = 4.0E0*I_NAI_I3x3z_Dyz_ab-2.0E0*2*I_NAI_Gx3z_Dyz_b;
  abcd[493] = 4.0E0*I_NAI_I2x4y_Dyz_ab-2.0E0*1*I_NAI_G4y_Dyz_b;
  abcd[494] = 4.0E0*I_NAI_I2x3yz_Dyz_ab-2.0E0*1*I_NAI_G3yz_Dyz_b;
  abcd[495] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2y2z_Dyz_b;
  abcd[496] = 4.0E0*I_NAI_I2xy3z_Dyz_ab-2.0E0*1*I_NAI_Gy3z_Dyz_b;
  abcd[497] = 4.0E0*I_NAI_I2x4z_Dyz_ab-2.0E0*1*I_NAI_G4z_Dyz_b;
  abcd[498] = 4.0E0*I_NAI_Ix5y_Dyz_ab;
  abcd[499] = 4.0E0*I_NAI_Ix4yz_Dyz_ab;
  abcd[500] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab;
  abcd[501] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab;
  abcd[502] = 4.0E0*I_NAI_Ixy4z_Dyz_ab;
  abcd[503] = 4.0E0*I_NAI_Ix5z_Dyz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[504] = 4.0E0*I_NAI_I6x_Dxz_ab-2.0E0*5*I_NAI_G4x_Dxz_b;
  abcd[505] = 4.0E0*I_NAI_I5xy_Dxz_ab-2.0E0*4*I_NAI_G3xy_Dxz_b;
  abcd[506] = 4.0E0*I_NAI_I5xz_Dxz_ab-2.0E0*4*I_NAI_G3xz_Dxz_b;
  abcd[507] = 4.0E0*I_NAI_I4x2y_Dxz_ab-2.0E0*3*I_NAI_G2x2y_Dxz_b;
  abcd[508] = 4.0E0*I_NAI_I4xyz_Dxz_ab-2.0E0*3*I_NAI_G2xyz_Dxz_b;
  abcd[509] = 4.0E0*I_NAI_I4x2z_Dxz_ab-2.0E0*3*I_NAI_G2x2z_Dxz_b;
  abcd[510] = 4.0E0*I_NAI_I3x3y_Dxz_ab-2.0E0*2*I_NAI_Gx3y_Dxz_b;
  abcd[511] = 4.0E0*I_NAI_I3x2yz_Dxz_ab-2.0E0*2*I_NAI_Gx2yz_Dxz_b;
  abcd[512] = 4.0E0*I_NAI_I3xy2z_Dxz_ab-2.0E0*2*I_NAI_Gxy2z_Dxz_b;
  abcd[513] = 4.0E0*I_NAI_I3x3z_Dxz_ab-2.0E0*2*I_NAI_Gx3z_Dxz_b;
  abcd[514] = 4.0E0*I_NAI_I2x4y_Dxz_ab-2.0E0*1*I_NAI_G4y_Dxz_b;
  abcd[515] = 4.0E0*I_NAI_I2x3yz_Dxz_ab-2.0E0*1*I_NAI_G3yz_Dxz_b;
  abcd[516] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2y2z_Dxz_b;
  abcd[517] = 4.0E0*I_NAI_I2xy3z_Dxz_ab-2.0E0*1*I_NAI_Gy3z_Dxz_b;
  abcd[518] = 4.0E0*I_NAI_I2x4z_Dxz_ab-2.0E0*1*I_NAI_G4z_Dxz_b;
  abcd[519] = 4.0E0*I_NAI_Ix5y_Dxz_ab;
  abcd[520] = 4.0E0*I_NAI_Ix4yz_Dxz_ab;
  abcd[521] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab;
  abcd[522] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab;
  abcd[523] = 4.0E0*I_NAI_Ixy4z_Dxz_ab;
  abcd[524] = 4.0E0*I_NAI_Ix5z_Dxz_ab;
  abcd[525] = 4.0E0*I_NAI_I6x_Dyz_ab-2.0E0*5*I_NAI_G4x_Dyz_b;
  abcd[526] = 4.0E0*I_NAI_I5xy_Dyz_ab-2.0E0*4*I_NAI_G3xy_Dyz_b;
  abcd[527] = 4.0E0*I_NAI_I5xz_Dyz_ab-2.0E0*4*I_NAI_G3xz_Dyz_b;
  abcd[528] = 4.0E0*I_NAI_I4x2y_Dyz_ab-2.0E0*3*I_NAI_G2x2y_Dyz_b;
  abcd[529] = 4.0E0*I_NAI_I4xyz_Dyz_ab-2.0E0*3*I_NAI_G2xyz_Dyz_b;
  abcd[530] = 4.0E0*I_NAI_I4x2z_Dyz_ab-2.0E0*3*I_NAI_G2x2z_Dyz_b;
  abcd[531] = 4.0E0*I_NAI_I3x3y_Dyz_ab-2.0E0*2*I_NAI_Gx3y_Dyz_b;
  abcd[532] = 4.0E0*I_NAI_I3x2yz_Dyz_ab-2.0E0*2*I_NAI_Gx2yz_Dyz_b;
  abcd[533] = 4.0E0*I_NAI_I3xy2z_Dyz_ab-2.0E0*2*I_NAI_Gxy2z_Dyz_b;
  abcd[534] = 4.0E0*I_NAI_I3x3z_Dyz_ab-2.0E0*2*I_NAI_Gx3z_Dyz_b;
  abcd[535] = 4.0E0*I_NAI_I2x4y_Dyz_ab-2.0E0*1*I_NAI_G4y_Dyz_b;
  abcd[536] = 4.0E0*I_NAI_I2x3yz_Dyz_ab-2.0E0*1*I_NAI_G3yz_Dyz_b;
  abcd[537] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2y2z_Dyz_b;
  abcd[538] = 4.0E0*I_NAI_I2xy3z_Dyz_ab-2.0E0*1*I_NAI_Gy3z_Dyz_b;
  abcd[539] = 4.0E0*I_NAI_I2x4z_Dyz_ab-2.0E0*1*I_NAI_G4z_Dyz_b;
  abcd[540] = 4.0E0*I_NAI_Ix5y_Dyz_ab;
  abcd[541] = 4.0E0*I_NAI_Ix4yz_Dyz_ab;
  abcd[542] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab;
  abcd[543] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab;
  abcd[544] = 4.0E0*I_NAI_Ixy4z_Dyz_ab;
  abcd[545] = 4.0E0*I_NAI_Ix5z_Dyz_ab;
  abcd[546] = 4.0E0*I_NAI_I6x_D2z_ab-2.0E0*1*I_NAI_I6x_S_a-2.0E0*5*I_NAI_G4x_D2z_b+5*1*I_NAI_G4x_S;
  abcd[547] = 4.0E0*I_NAI_I5xy_D2z_ab-2.0E0*1*I_NAI_I5xy_S_a-2.0E0*4*I_NAI_G3xy_D2z_b+4*1*I_NAI_G3xy_S;
  abcd[548] = 4.0E0*I_NAI_I5xz_D2z_ab-2.0E0*1*I_NAI_I5xz_S_a-2.0E0*4*I_NAI_G3xz_D2z_b+4*1*I_NAI_G3xz_S;
  abcd[549] = 4.0E0*I_NAI_I4x2y_D2z_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*3*I_NAI_G2x2y_D2z_b+3*1*I_NAI_G2x2y_S;
  abcd[550] = 4.0E0*I_NAI_I4xyz_D2z_ab-2.0E0*1*I_NAI_I4xyz_S_a-2.0E0*3*I_NAI_G2xyz_D2z_b+3*1*I_NAI_G2xyz_S;
  abcd[551] = 4.0E0*I_NAI_I4x2z_D2z_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*3*I_NAI_G2x2z_D2z_b+3*1*I_NAI_G2x2z_S;
  abcd[552] = 4.0E0*I_NAI_I3x3y_D2z_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_Gx3y_D2z_b+2*1*I_NAI_Gx3y_S;
  abcd[553] = 4.0E0*I_NAI_I3x2yz_D2z_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*2*I_NAI_Gx2yz_D2z_b+2*1*I_NAI_Gx2yz_S;
  abcd[554] = 4.0E0*I_NAI_I3xy2z_D2z_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2z_b+2*1*I_NAI_Gxy2z_S;
  abcd[555] = 4.0E0*I_NAI_I3x3z_D2z_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_Gx3z_D2z_b+2*1*I_NAI_Gx3z_S;
  abcd[556] = 4.0E0*I_NAI_I2x4y_D2z_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*1*I_NAI_G4y_D2z_b+1*I_NAI_G4y_S;
  abcd[557] = 4.0E0*I_NAI_I2x3yz_D2z_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*1*I_NAI_G3yz_D2z_b+1*I_NAI_G3yz_S;
  abcd[558] = 4.0E0*I_NAI_I2x2y2z_D2z_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2y2z_D2z_b+1*I_NAI_G2y2z_S;
  abcd[559] = 4.0E0*I_NAI_I2xy3z_D2z_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*1*I_NAI_Gy3z_D2z_b+1*I_NAI_Gy3z_S;
  abcd[560] = 4.0E0*I_NAI_I2x4z_D2z_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*1*I_NAI_G4z_D2z_b+1*I_NAI_G4z_S;
  abcd[561] = 4.0E0*I_NAI_Ix5y_D2z_ab-2.0E0*1*I_NAI_Ix5y_S_a;
  abcd[562] = 4.0E0*I_NAI_Ix4yz_D2z_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[563] = 4.0E0*I_NAI_Ix3y2z_D2z_ab-2.0E0*1*I_NAI_Ix3y2z_S_a;
  abcd[564] = 4.0E0*I_NAI_Ix2y3z_D2z_ab-2.0E0*1*I_NAI_Ix2y3z_S_a;
  abcd[565] = 4.0E0*I_NAI_Ixy4z_D2z_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[566] = 4.0E0*I_NAI_Ix5z_D2z_ab-2.0E0*1*I_NAI_Ix5z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[567] = 4.0E0*I_NAI_I5xy_D2x_ab-2.0E0*1*I_NAI_I5xy_S_a;
  abcd[568] = 4.0E0*I_NAI_I4x2y_D2x_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*1*I_NAI_G4x_D2x_b+1*I_NAI_G4x_S;
  abcd[569] = 4.0E0*I_NAI_I4xyz_D2x_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[570] = 4.0E0*I_NAI_I3x3y_D2x_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_G3xy_D2x_b+2*1*I_NAI_G3xy_S;
  abcd[571] = 4.0E0*I_NAI_I3x2yz_D2x_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*1*I_NAI_G3xz_D2x_b+1*I_NAI_G3xz_S;
  abcd[572] = 4.0E0*I_NAI_I3xy2z_D2x_ab-2.0E0*1*I_NAI_I3xy2z_S_a;
  abcd[573] = 4.0E0*I_NAI_I2x4y_D2x_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*3*I_NAI_G2x2y_D2x_b+3*1*I_NAI_G2x2y_S;
  abcd[574] = 4.0E0*I_NAI_I2x3yz_D2x_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*2*I_NAI_G2xyz_D2x_b+2*1*I_NAI_G2xyz_S;
  abcd[575] = 4.0E0*I_NAI_I2x2y2z_D2x_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2z_D2x_b+1*I_NAI_G2x2z_S;
  abcd[576] = 4.0E0*I_NAI_I2xy3z_D2x_ab-2.0E0*1*I_NAI_I2xy3z_S_a;
  abcd[577] = 4.0E0*I_NAI_Ix5y_D2x_ab-2.0E0*1*I_NAI_Ix5y_S_a-2.0E0*4*I_NAI_Gx3y_D2x_b+4*1*I_NAI_Gx3y_S;
  abcd[578] = 4.0E0*I_NAI_Ix4yz_D2x_ab-2.0E0*1*I_NAI_Ix4yz_S_a-2.0E0*3*I_NAI_Gx2yz_D2x_b+3*1*I_NAI_Gx2yz_S;
  abcd[579] = 4.0E0*I_NAI_Ix3y2z_D2x_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2x_b+2*1*I_NAI_Gxy2z_S;
  abcd[580] = 4.0E0*I_NAI_Ix2y3z_D2x_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*1*I_NAI_Gx3z_D2x_b+1*I_NAI_Gx3z_S;
  abcd[581] = 4.0E0*I_NAI_Ixy4z_D2x_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[582] = 4.0E0*I_NAI_I6y_D2x_ab-2.0E0*1*I_NAI_I6y_S_a-2.0E0*5*I_NAI_G4y_D2x_b+5*1*I_NAI_G4y_S;
  abcd[583] = 4.0E0*I_NAI_I5yz_D2x_ab-2.0E0*1*I_NAI_I5yz_S_a-2.0E0*4*I_NAI_G3yz_D2x_b+4*1*I_NAI_G3yz_S;
  abcd[584] = 4.0E0*I_NAI_I4y2z_D2x_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*3*I_NAI_G2y2z_D2x_b+3*1*I_NAI_G2y2z_S;
  abcd[585] = 4.0E0*I_NAI_I3y3z_D2x_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_Gy3z_D2x_b+2*1*I_NAI_Gy3z_S;
  abcd[586] = 4.0E0*I_NAI_I2y4z_D2x_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*1*I_NAI_G4z_D2x_b+1*I_NAI_G4z_S;
  abcd[587] = 4.0E0*I_NAI_Iy5z_D2x_ab-2.0E0*1*I_NAI_Iy5z_S_a;
  abcd[588] = 4.0E0*I_NAI_I5xy_Dxy_ab;
  abcd[589] = 4.0E0*I_NAI_I4x2y_Dxy_ab-2.0E0*1*I_NAI_G4x_Dxy_b;
  abcd[590] = 4.0E0*I_NAI_I4xyz_Dxy_ab;
  abcd[591] = 4.0E0*I_NAI_I3x3y_Dxy_ab-2.0E0*2*I_NAI_G3xy_Dxy_b;
  abcd[592] = 4.0E0*I_NAI_I3x2yz_Dxy_ab-2.0E0*1*I_NAI_G3xz_Dxy_b;
  abcd[593] = 4.0E0*I_NAI_I3xy2z_Dxy_ab;
  abcd[594] = 4.0E0*I_NAI_I2x4y_Dxy_ab-2.0E0*3*I_NAI_G2x2y_Dxy_b;
  abcd[595] = 4.0E0*I_NAI_I2x3yz_Dxy_ab-2.0E0*2*I_NAI_G2xyz_Dxy_b;
  abcd[596] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2x2z_Dxy_b;
  abcd[597] = 4.0E0*I_NAI_I2xy3z_Dxy_ab;
  abcd[598] = 4.0E0*I_NAI_Ix5y_Dxy_ab-2.0E0*4*I_NAI_Gx3y_Dxy_b;
  abcd[599] = 4.0E0*I_NAI_Ix4yz_Dxy_ab-2.0E0*3*I_NAI_Gx2yz_Dxy_b;
  abcd[600] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab-2.0E0*2*I_NAI_Gxy2z_Dxy_b;
  abcd[601] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab-2.0E0*1*I_NAI_Gx3z_Dxy_b;
  abcd[602] = 4.0E0*I_NAI_Ixy4z_Dxy_ab;
  abcd[603] = 4.0E0*I_NAI_I6y_Dxy_ab-2.0E0*5*I_NAI_G4y_Dxy_b;
  abcd[604] = 4.0E0*I_NAI_I5yz_Dxy_ab-2.0E0*4*I_NAI_G3yz_Dxy_b;
  abcd[605] = 4.0E0*I_NAI_I4y2z_Dxy_ab-2.0E0*3*I_NAI_G2y2z_Dxy_b;
  abcd[606] = 4.0E0*I_NAI_I3y3z_Dxy_ab-2.0E0*2*I_NAI_Gy3z_Dxy_b;
  abcd[607] = 4.0E0*I_NAI_I2y4z_Dxy_ab-2.0E0*1*I_NAI_G4z_Dxy_b;
  abcd[608] = 4.0E0*I_NAI_Iy5z_Dxy_ab;
  abcd[609] = 4.0E0*I_NAI_I5xy_Dxz_ab;
  abcd[610] = 4.0E0*I_NAI_I4x2y_Dxz_ab-2.0E0*1*I_NAI_G4x_Dxz_b;
  abcd[611] = 4.0E0*I_NAI_I4xyz_Dxz_ab;
  abcd[612] = 4.0E0*I_NAI_I3x3y_Dxz_ab-2.0E0*2*I_NAI_G3xy_Dxz_b;
  abcd[613] = 4.0E0*I_NAI_I3x2yz_Dxz_ab-2.0E0*1*I_NAI_G3xz_Dxz_b;
  abcd[614] = 4.0E0*I_NAI_I3xy2z_Dxz_ab;
  abcd[615] = 4.0E0*I_NAI_I2x4y_Dxz_ab-2.0E0*3*I_NAI_G2x2y_Dxz_b;
  abcd[616] = 4.0E0*I_NAI_I2x3yz_Dxz_ab-2.0E0*2*I_NAI_G2xyz_Dxz_b;
  abcd[617] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2x2z_Dxz_b;
  abcd[618] = 4.0E0*I_NAI_I2xy3z_Dxz_ab;
  abcd[619] = 4.0E0*I_NAI_Ix5y_Dxz_ab-2.0E0*4*I_NAI_Gx3y_Dxz_b;
  abcd[620] = 4.0E0*I_NAI_Ix4yz_Dxz_ab-2.0E0*3*I_NAI_Gx2yz_Dxz_b;
  abcd[621] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab-2.0E0*2*I_NAI_Gxy2z_Dxz_b;
  abcd[622] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab-2.0E0*1*I_NAI_Gx3z_Dxz_b;
  abcd[623] = 4.0E0*I_NAI_Ixy4z_Dxz_ab;
  abcd[624] = 4.0E0*I_NAI_I6y_Dxz_ab-2.0E0*5*I_NAI_G4y_Dxz_b;
  abcd[625] = 4.0E0*I_NAI_I5yz_Dxz_ab-2.0E0*4*I_NAI_G3yz_Dxz_b;
  abcd[626] = 4.0E0*I_NAI_I4y2z_Dxz_ab-2.0E0*3*I_NAI_G2y2z_Dxz_b;
  abcd[627] = 4.0E0*I_NAI_I3y3z_Dxz_ab-2.0E0*2*I_NAI_Gy3z_Dxz_b;
  abcd[628] = 4.0E0*I_NAI_I2y4z_Dxz_ab-2.0E0*1*I_NAI_G4z_Dxz_b;
  abcd[629] = 4.0E0*I_NAI_Iy5z_Dxz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[630] = 4.0E0*I_NAI_I5xy_Dxy_ab;
  abcd[631] = 4.0E0*I_NAI_I4x2y_Dxy_ab-2.0E0*1*I_NAI_G4x_Dxy_b;
  abcd[632] = 4.0E0*I_NAI_I4xyz_Dxy_ab;
  abcd[633] = 4.0E0*I_NAI_I3x3y_Dxy_ab-2.0E0*2*I_NAI_G3xy_Dxy_b;
  abcd[634] = 4.0E0*I_NAI_I3x2yz_Dxy_ab-2.0E0*1*I_NAI_G3xz_Dxy_b;
  abcd[635] = 4.0E0*I_NAI_I3xy2z_Dxy_ab;
  abcd[636] = 4.0E0*I_NAI_I2x4y_Dxy_ab-2.0E0*3*I_NAI_G2x2y_Dxy_b;
  abcd[637] = 4.0E0*I_NAI_I2x3yz_Dxy_ab-2.0E0*2*I_NAI_G2xyz_Dxy_b;
  abcd[638] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2x2z_Dxy_b;
  abcd[639] = 4.0E0*I_NAI_I2xy3z_Dxy_ab;
  abcd[640] = 4.0E0*I_NAI_Ix5y_Dxy_ab-2.0E0*4*I_NAI_Gx3y_Dxy_b;
  abcd[641] = 4.0E0*I_NAI_Ix4yz_Dxy_ab-2.0E0*3*I_NAI_Gx2yz_Dxy_b;
  abcd[642] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab-2.0E0*2*I_NAI_Gxy2z_Dxy_b;
  abcd[643] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab-2.0E0*1*I_NAI_Gx3z_Dxy_b;
  abcd[644] = 4.0E0*I_NAI_Ixy4z_Dxy_ab;
  abcd[645] = 4.0E0*I_NAI_I6y_Dxy_ab-2.0E0*5*I_NAI_G4y_Dxy_b;
  abcd[646] = 4.0E0*I_NAI_I5yz_Dxy_ab-2.0E0*4*I_NAI_G3yz_Dxy_b;
  abcd[647] = 4.0E0*I_NAI_I4y2z_Dxy_ab-2.0E0*3*I_NAI_G2y2z_Dxy_b;
  abcd[648] = 4.0E0*I_NAI_I3y3z_Dxy_ab-2.0E0*2*I_NAI_Gy3z_Dxy_b;
  abcd[649] = 4.0E0*I_NAI_I2y4z_Dxy_ab-2.0E0*1*I_NAI_G4z_Dxy_b;
  abcd[650] = 4.0E0*I_NAI_Iy5z_Dxy_ab;
  abcd[651] = 4.0E0*I_NAI_I5xy_D2y_ab-2.0E0*1*I_NAI_I5xy_S_a;
  abcd[652] = 4.0E0*I_NAI_I4x2y_D2y_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*1*I_NAI_G4x_D2y_b+1*I_NAI_G4x_S;
  abcd[653] = 4.0E0*I_NAI_I4xyz_D2y_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[654] = 4.0E0*I_NAI_I3x3y_D2y_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_G3xy_D2y_b+2*1*I_NAI_G3xy_S;
  abcd[655] = 4.0E0*I_NAI_I3x2yz_D2y_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*1*I_NAI_G3xz_D2y_b+1*I_NAI_G3xz_S;
  abcd[656] = 4.0E0*I_NAI_I3xy2z_D2y_ab-2.0E0*1*I_NAI_I3xy2z_S_a;
  abcd[657] = 4.0E0*I_NAI_I2x4y_D2y_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*3*I_NAI_G2x2y_D2y_b+3*1*I_NAI_G2x2y_S;
  abcd[658] = 4.0E0*I_NAI_I2x3yz_D2y_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*2*I_NAI_G2xyz_D2y_b+2*1*I_NAI_G2xyz_S;
  abcd[659] = 4.0E0*I_NAI_I2x2y2z_D2y_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2z_D2y_b+1*I_NAI_G2x2z_S;
  abcd[660] = 4.0E0*I_NAI_I2xy3z_D2y_ab-2.0E0*1*I_NAI_I2xy3z_S_a;
  abcd[661] = 4.0E0*I_NAI_Ix5y_D2y_ab-2.0E0*1*I_NAI_Ix5y_S_a-2.0E0*4*I_NAI_Gx3y_D2y_b+4*1*I_NAI_Gx3y_S;
  abcd[662] = 4.0E0*I_NAI_Ix4yz_D2y_ab-2.0E0*1*I_NAI_Ix4yz_S_a-2.0E0*3*I_NAI_Gx2yz_D2y_b+3*1*I_NAI_Gx2yz_S;
  abcd[663] = 4.0E0*I_NAI_Ix3y2z_D2y_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2y_b+2*1*I_NAI_Gxy2z_S;
  abcd[664] = 4.0E0*I_NAI_Ix2y3z_D2y_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*1*I_NAI_Gx3z_D2y_b+1*I_NAI_Gx3z_S;
  abcd[665] = 4.0E0*I_NAI_Ixy4z_D2y_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[666] = 4.0E0*I_NAI_I6y_D2y_ab-2.0E0*1*I_NAI_I6y_S_a-2.0E0*5*I_NAI_G4y_D2y_b+5*1*I_NAI_G4y_S;
  abcd[667] = 4.0E0*I_NAI_I5yz_D2y_ab-2.0E0*1*I_NAI_I5yz_S_a-2.0E0*4*I_NAI_G3yz_D2y_b+4*1*I_NAI_G3yz_S;
  abcd[668] = 4.0E0*I_NAI_I4y2z_D2y_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*3*I_NAI_G2y2z_D2y_b+3*1*I_NAI_G2y2z_S;
  abcd[669] = 4.0E0*I_NAI_I3y3z_D2y_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_Gy3z_D2y_b+2*1*I_NAI_Gy3z_S;
  abcd[670] = 4.0E0*I_NAI_I2y4z_D2y_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*1*I_NAI_G4z_D2y_b+1*I_NAI_G4z_S;
  abcd[671] = 4.0E0*I_NAI_Iy5z_D2y_ab-2.0E0*1*I_NAI_Iy5z_S_a;
  abcd[672] = 4.0E0*I_NAI_I5xy_Dyz_ab;
  abcd[673] = 4.0E0*I_NAI_I4x2y_Dyz_ab-2.0E0*1*I_NAI_G4x_Dyz_b;
  abcd[674] = 4.0E0*I_NAI_I4xyz_Dyz_ab;
  abcd[675] = 4.0E0*I_NAI_I3x3y_Dyz_ab-2.0E0*2*I_NAI_G3xy_Dyz_b;
  abcd[676] = 4.0E0*I_NAI_I3x2yz_Dyz_ab-2.0E0*1*I_NAI_G3xz_Dyz_b;
  abcd[677] = 4.0E0*I_NAI_I3xy2z_Dyz_ab;
  abcd[678] = 4.0E0*I_NAI_I2x4y_Dyz_ab-2.0E0*3*I_NAI_G2x2y_Dyz_b;
  abcd[679] = 4.0E0*I_NAI_I2x3yz_Dyz_ab-2.0E0*2*I_NAI_G2xyz_Dyz_b;
  abcd[680] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2x2z_Dyz_b;
  abcd[681] = 4.0E0*I_NAI_I2xy3z_Dyz_ab;
  abcd[682] = 4.0E0*I_NAI_Ix5y_Dyz_ab-2.0E0*4*I_NAI_Gx3y_Dyz_b;
  abcd[683] = 4.0E0*I_NAI_Ix4yz_Dyz_ab-2.0E0*3*I_NAI_Gx2yz_Dyz_b;
  abcd[684] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab-2.0E0*2*I_NAI_Gxy2z_Dyz_b;
  abcd[685] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab-2.0E0*1*I_NAI_Gx3z_Dyz_b;
  abcd[686] = 4.0E0*I_NAI_Ixy4z_Dyz_ab;
  abcd[687] = 4.0E0*I_NAI_I6y_Dyz_ab-2.0E0*5*I_NAI_G4y_Dyz_b;
  abcd[688] = 4.0E0*I_NAI_I5yz_Dyz_ab-2.0E0*4*I_NAI_G3yz_Dyz_b;
  abcd[689] = 4.0E0*I_NAI_I4y2z_Dyz_ab-2.0E0*3*I_NAI_G2y2z_Dyz_b;
  abcd[690] = 4.0E0*I_NAI_I3y3z_Dyz_ab-2.0E0*2*I_NAI_Gy3z_Dyz_b;
  abcd[691] = 4.0E0*I_NAI_I2y4z_Dyz_ab-2.0E0*1*I_NAI_G4z_Dyz_b;
  abcd[692] = 4.0E0*I_NAI_Iy5z_Dyz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[693] = 4.0E0*I_NAI_I5xy_Dxz_ab;
  abcd[694] = 4.0E0*I_NAI_I4x2y_Dxz_ab-2.0E0*1*I_NAI_G4x_Dxz_b;
  abcd[695] = 4.0E0*I_NAI_I4xyz_Dxz_ab;
  abcd[696] = 4.0E0*I_NAI_I3x3y_Dxz_ab-2.0E0*2*I_NAI_G3xy_Dxz_b;
  abcd[697] = 4.0E0*I_NAI_I3x2yz_Dxz_ab-2.0E0*1*I_NAI_G3xz_Dxz_b;
  abcd[698] = 4.0E0*I_NAI_I3xy2z_Dxz_ab;
  abcd[699] = 4.0E0*I_NAI_I2x4y_Dxz_ab-2.0E0*3*I_NAI_G2x2y_Dxz_b;
  abcd[700] = 4.0E0*I_NAI_I2x3yz_Dxz_ab-2.0E0*2*I_NAI_G2xyz_Dxz_b;
  abcd[701] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2x2z_Dxz_b;
  abcd[702] = 4.0E0*I_NAI_I2xy3z_Dxz_ab;
  abcd[703] = 4.0E0*I_NAI_Ix5y_Dxz_ab-2.0E0*4*I_NAI_Gx3y_Dxz_b;
  abcd[704] = 4.0E0*I_NAI_Ix4yz_Dxz_ab-2.0E0*3*I_NAI_Gx2yz_Dxz_b;
  abcd[705] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab-2.0E0*2*I_NAI_Gxy2z_Dxz_b;
  abcd[706] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab-2.0E0*1*I_NAI_Gx3z_Dxz_b;
  abcd[707] = 4.0E0*I_NAI_Ixy4z_Dxz_ab;
  abcd[708] = 4.0E0*I_NAI_I6y_Dxz_ab-2.0E0*5*I_NAI_G4y_Dxz_b;
  abcd[709] = 4.0E0*I_NAI_I5yz_Dxz_ab-2.0E0*4*I_NAI_G3yz_Dxz_b;
  abcd[710] = 4.0E0*I_NAI_I4y2z_Dxz_ab-2.0E0*3*I_NAI_G2y2z_Dxz_b;
  abcd[711] = 4.0E0*I_NAI_I3y3z_Dxz_ab-2.0E0*2*I_NAI_Gy3z_Dxz_b;
  abcd[712] = 4.0E0*I_NAI_I2y4z_Dxz_ab-2.0E0*1*I_NAI_G4z_Dxz_b;
  abcd[713] = 4.0E0*I_NAI_Iy5z_Dxz_ab;
  abcd[714] = 4.0E0*I_NAI_I5xy_Dyz_ab;
  abcd[715] = 4.0E0*I_NAI_I4x2y_Dyz_ab-2.0E0*1*I_NAI_G4x_Dyz_b;
  abcd[716] = 4.0E0*I_NAI_I4xyz_Dyz_ab;
  abcd[717] = 4.0E0*I_NAI_I3x3y_Dyz_ab-2.0E0*2*I_NAI_G3xy_Dyz_b;
  abcd[718] = 4.0E0*I_NAI_I3x2yz_Dyz_ab-2.0E0*1*I_NAI_G3xz_Dyz_b;
  abcd[719] = 4.0E0*I_NAI_I3xy2z_Dyz_ab;
  abcd[720] = 4.0E0*I_NAI_I2x4y_Dyz_ab-2.0E0*3*I_NAI_G2x2y_Dyz_b;
  abcd[721] = 4.0E0*I_NAI_I2x3yz_Dyz_ab-2.0E0*2*I_NAI_G2xyz_Dyz_b;
  abcd[722] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2x2z_Dyz_b;
  abcd[723] = 4.0E0*I_NAI_I2xy3z_Dyz_ab;
  abcd[724] = 4.0E0*I_NAI_Ix5y_Dyz_ab-2.0E0*4*I_NAI_Gx3y_Dyz_b;
  abcd[725] = 4.0E0*I_NAI_Ix4yz_Dyz_ab-2.0E0*3*I_NAI_Gx2yz_Dyz_b;
  abcd[726] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab-2.0E0*2*I_NAI_Gxy2z_Dyz_b;
  abcd[727] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab-2.0E0*1*I_NAI_Gx3z_Dyz_b;
  abcd[728] = 4.0E0*I_NAI_Ixy4z_Dyz_ab;
  abcd[729] = 4.0E0*I_NAI_I6y_Dyz_ab-2.0E0*5*I_NAI_G4y_Dyz_b;
  abcd[730] = 4.0E0*I_NAI_I5yz_Dyz_ab-2.0E0*4*I_NAI_G3yz_Dyz_b;
  abcd[731] = 4.0E0*I_NAI_I4y2z_Dyz_ab-2.0E0*3*I_NAI_G2y2z_Dyz_b;
  abcd[732] = 4.0E0*I_NAI_I3y3z_Dyz_ab-2.0E0*2*I_NAI_Gy3z_Dyz_b;
  abcd[733] = 4.0E0*I_NAI_I2y4z_Dyz_ab-2.0E0*1*I_NAI_G4z_Dyz_b;
  abcd[734] = 4.0E0*I_NAI_Iy5z_Dyz_ab;
  abcd[735] = 4.0E0*I_NAI_I5xy_D2z_ab-2.0E0*1*I_NAI_I5xy_S_a;
  abcd[736] = 4.0E0*I_NAI_I4x2y_D2z_ab-2.0E0*1*I_NAI_I4x2y_S_a-2.0E0*1*I_NAI_G4x_D2z_b+1*I_NAI_G4x_S;
  abcd[737] = 4.0E0*I_NAI_I4xyz_D2z_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[738] = 4.0E0*I_NAI_I3x3y_D2z_ab-2.0E0*1*I_NAI_I3x3y_S_a-2.0E0*2*I_NAI_G3xy_D2z_b+2*1*I_NAI_G3xy_S;
  abcd[739] = 4.0E0*I_NAI_I3x2yz_D2z_ab-2.0E0*1*I_NAI_I3x2yz_S_a-2.0E0*1*I_NAI_G3xz_D2z_b+1*I_NAI_G3xz_S;
  abcd[740] = 4.0E0*I_NAI_I3xy2z_D2z_ab-2.0E0*1*I_NAI_I3xy2z_S_a;
  abcd[741] = 4.0E0*I_NAI_I2x4y_D2z_ab-2.0E0*1*I_NAI_I2x4y_S_a-2.0E0*3*I_NAI_G2x2y_D2z_b+3*1*I_NAI_G2x2y_S;
  abcd[742] = 4.0E0*I_NAI_I2x3yz_D2z_ab-2.0E0*1*I_NAI_I2x3yz_S_a-2.0E0*2*I_NAI_G2xyz_D2z_b+2*1*I_NAI_G2xyz_S;
  abcd[743] = 4.0E0*I_NAI_I2x2y2z_D2z_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2z_D2z_b+1*I_NAI_G2x2z_S;
  abcd[744] = 4.0E0*I_NAI_I2xy3z_D2z_ab-2.0E0*1*I_NAI_I2xy3z_S_a;
  abcd[745] = 4.0E0*I_NAI_Ix5y_D2z_ab-2.0E0*1*I_NAI_Ix5y_S_a-2.0E0*4*I_NAI_Gx3y_D2z_b+4*1*I_NAI_Gx3y_S;
  abcd[746] = 4.0E0*I_NAI_Ix4yz_D2z_ab-2.0E0*1*I_NAI_Ix4yz_S_a-2.0E0*3*I_NAI_Gx2yz_D2z_b+3*1*I_NAI_Gx2yz_S;
  abcd[747] = 4.0E0*I_NAI_Ix3y2z_D2z_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*2*I_NAI_Gxy2z_D2z_b+2*1*I_NAI_Gxy2z_S;
  abcd[748] = 4.0E0*I_NAI_Ix2y3z_D2z_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*1*I_NAI_Gx3z_D2z_b+1*I_NAI_Gx3z_S;
  abcd[749] = 4.0E0*I_NAI_Ixy4z_D2z_ab-2.0E0*1*I_NAI_Ixy4z_S_a;
  abcd[750] = 4.0E0*I_NAI_I6y_D2z_ab-2.0E0*1*I_NAI_I6y_S_a-2.0E0*5*I_NAI_G4y_D2z_b+5*1*I_NAI_G4y_S;
  abcd[751] = 4.0E0*I_NAI_I5yz_D2z_ab-2.0E0*1*I_NAI_I5yz_S_a-2.0E0*4*I_NAI_G3yz_D2z_b+4*1*I_NAI_G3yz_S;
  abcd[752] = 4.0E0*I_NAI_I4y2z_D2z_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*3*I_NAI_G2y2z_D2z_b+3*1*I_NAI_G2y2z_S;
  abcd[753] = 4.0E0*I_NAI_I3y3z_D2z_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_Gy3z_D2z_b+2*1*I_NAI_Gy3z_S;
  abcd[754] = 4.0E0*I_NAI_I2y4z_D2z_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*1*I_NAI_G4z_D2z_b+1*I_NAI_G4z_S;
  abcd[755] = 4.0E0*I_NAI_Iy5z_D2z_ab-2.0E0*1*I_NAI_Iy5z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[756] = 4.0E0*I_NAI_I5xz_D2x_ab-2.0E0*1*I_NAI_I5xz_S_a;
  abcd[757] = 4.0E0*I_NAI_I4xyz_D2x_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[758] = 4.0E0*I_NAI_I4x2z_D2x_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*1*I_NAI_G4x_D2x_b+1*I_NAI_G4x_S;
  abcd[759] = 4.0E0*I_NAI_I3x2yz_D2x_ab-2.0E0*1*I_NAI_I3x2yz_S_a;
  abcd[760] = 4.0E0*I_NAI_I3xy2z_D2x_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*1*I_NAI_G3xy_D2x_b+1*I_NAI_G3xy_S;
  abcd[761] = 4.0E0*I_NAI_I3x3z_D2x_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_G3xz_D2x_b+2*1*I_NAI_G3xz_S;
  abcd[762] = 4.0E0*I_NAI_I2x3yz_D2x_ab-2.0E0*1*I_NAI_I2x3yz_S_a;
  abcd[763] = 4.0E0*I_NAI_I2x2y2z_D2x_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2y_D2x_b+1*I_NAI_G2x2y_S;
  abcd[764] = 4.0E0*I_NAI_I2xy3z_D2x_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*2*I_NAI_G2xyz_D2x_b+2*1*I_NAI_G2xyz_S;
  abcd[765] = 4.0E0*I_NAI_I2x4z_D2x_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*3*I_NAI_G2x2z_D2x_b+3*1*I_NAI_G2x2z_S;
  abcd[766] = 4.0E0*I_NAI_Ix4yz_D2x_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[767] = 4.0E0*I_NAI_Ix3y2z_D2x_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*1*I_NAI_Gx3y_D2x_b+1*I_NAI_Gx3y_S;
  abcd[768] = 4.0E0*I_NAI_Ix2y3z_D2x_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*2*I_NAI_Gx2yz_D2x_b+2*1*I_NAI_Gx2yz_S;
  abcd[769] = 4.0E0*I_NAI_Ixy4z_D2x_ab-2.0E0*1*I_NAI_Ixy4z_S_a-2.0E0*3*I_NAI_Gxy2z_D2x_b+3*1*I_NAI_Gxy2z_S;
  abcd[770] = 4.0E0*I_NAI_Ix5z_D2x_ab-2.0E0*1*I_NAI_Ix5z_S_a-2.0E0*4*I_NAI_Gx3z_D2x_b+4*1*I_NAI_Gx3z_S;
  abcd[771] = 4.0E0*I_NAI_I5yz_D2x_ab-2.0E0*1*I_NAI_I5yz_S_a;
  abcd[772] = 4.0E0*I_NAI_I4y2z_D2x_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*1*I_NAI_G4y_D2x_b+1*I_NAI_G4y_S;
  abcd[773] = 4.0E0*I_NAI_I3y3z_D2x_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_G3yz_D2x_b+2*1*I_NAI_G3yz_S;
  abcd[774] = 4.0E0*I_NAI_I2y4z_D2x_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*3*I_NAI_G2y2z_D2x_b+3*1*I_NAI_G2y2z_S;
  abcd[775] = 4.0E0*I_NAI_Iy5z_D2x_ab-2.0E0*1*I_NAI_Iy5z_S_a-2.0E0*4*I_NAI_Gy3z_D2x_b+4*1*I_NAI_Gy3z_S;
  abcd[776] = 4.0E0*I_NAI_I6z_D2x_ab-2.0E0*1*I_NAI_I6z_S_a-2.0E0*5*I_NAI_G4z_D2x_b+5*1*I_NAI_G4z_S;
  abcd[777] = 4.0E0*I_NAI_I5xz_Dxy_ab;
  abcd[778] = 4.0E0*I_NAI_I4xyz_Dxy_ab;
  abcd[779] = 4.0E0*I_NAI_I4x2z_Dxy_ab-2.0E0*1*I_NAI_G4x_Dxy_b;
  abcd[780] = 4.0E0*I_NAI_I3x2yz_Dxy_ab;
  abcd[781] = 4.0E0*I_NAI_I3xy2z_Dxy_ab-2.0E0*1*I_NAI_G3xy_Dxy_b;
  abcd[782] = 4.0E0*I_NAI_I3x3z_Dxy_ab-2.0E0*2*I_NAI_G3xz_Dxy_b;
  abcd[783] = 4.0E0*I_NAI_I2x3yz_Dxy_ab;
  abcd[784] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2x2y_Dxy_b;
  abcd[785] = 4.0E0*I_NAI_I2xy3z_Dxy_ab-2.0E0*2*I_NAI_G2xyz_Dxy_b;
  abcd[786] = 4.0E0*I_NAI_I2x4z_Dxy_ab-2.0E0*3*I_NAI_G2x2z_Dxy_b;
  abcd[787] = 4.0E0*I_NAI_Ix4yz_Dxy_ab;
  abcd[788] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab-2.0E0*1*I_NAI_Gx3y_Dxy_b;
  abcd[789] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab-2.0E0*2*I_NAI_Gx2yz_Dxy_b;
  abcd[790] = 4.0E0*I_NAI_Ixy4z_Dxy_ab-2.0E0*3*I_NAI_Gxy2z_Dxy_b;
  abcd[791] = 4.0E0*I_NAI_Ix5z_Dxy_ab-2.0E0*4*I_NAI_Gx3z_Dxy_b;
  abcd[792] = 4.0E0*I_NAI_I5yz_Dxy_ab;
  abcd[793] = 4.0E0*I_NAI_I4y2z_Dxy_ab-2.0E0*1*I_NAI_G4y_Dxy_b;
  abcd[794] = 4.0E0*I_NAI_I3y3z_Dxy_ab-2.0E0*2*I_NAI_G3yz_Dxy_b;
  abcd[795] = 4.0E0*I_NAI_I2y4z_Dxy_ab-2.0E0*3*I_NAI_G2y2z_Dxy_b;
  abcd[796] = 4.0E0*I_NAI_Iy5z_Dxy_ab-2.0E0*4*I_NAI_Gy3z_Dxy_b;
  abcd[797] = 4.0E0*I_NAI_I6z_Dxy_ab-2.0E0*5*I_NAI_G4z_Dxy_b;
  abcd[798] = 4.0E0*I_NAI_I5xz_Dxz_ab;
  abcd[799] = 4.0E0*I_NAI_I4xyz_Dxz_ab;
  abcd[800] = 4.0E0*I_NAI_I4x2z_Dxz_ab-2.0E0*1*I_NAI_G4x_Dxz_b;
  abcd[801] = 4.0E0*I_NAI_I3x2yz_Dxz_ab;
  abcd[802] = 4.0E0*I_NAI_I3xy2z_Dxz_ab-2.0E0*1*I_NAI_G3xy_Dxz_b;
  abcd[803] = 4.0E0*I_NAI_I3x3z_Dxz_ab-2.0E0*2*I_NAI_G3xz_Dxz_b;
  abcd[804] = 4.0E0*I_NAI_I2x3yz_Dxz_ab;
  abcd[805] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2x2y_Dxz_b;
  abcd[806] = 4.0E0*I_NAI_I2xy3z_Dxz_ab-2.0E0*2*I_NAI_G2xyz_Dxz_b;
  abcd[807] = 4.0E0*I_NAI_I2x4z_Dxz_ab-2.0E0*3*I_NAI_G2x2z_Dxz_b;
  abcd[808] = 4.0E0*I_NAI_Ix4yz_Dxz_ab;
  abcd[809] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab-2.0E0*1*I_NAI_Gx3y_Dxz_b;
  abcd[810] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab-2.0E0*2*I_NAI_Gx2yz_Dxz_b;
  abcd[811] = 4.0E0*I_NAI_Ixy4z_Dxz_ab-2.0E0*3*I_NAI_Gxy2z_Dxz_b;
  abcd[812] = 4.0E0*I_NAI_Ix5z_Dxz_ab-2.0E0*4*I_NAI_Gx3z_Dxz_b;
  abcd[813] = 4.0E0*I_NAI_I5yz_Dxz_ab;
  abcd[814] = 4.0E0*I_NAI_I4y2z_Dxz_ab-2.0E0*1*I_NAI_G4y_Dxz_b;
  abcd[815] = 4.0E0*I_NAI_I3y3z_Dxz_ab-2.0E0*2*I_NAI_G3yz_Dxz_b;
  abcd[816] = 4.0E0*I_NAI_I2y4z_Dxz_ab-2.0E0*3*I_NAI_G2y2z_Dxz_b;
  abcd[817] = 4.0E0*I_NAI_Iy5z_Dxz_ab-2.0E0*4*I_NAI_Gy3z_Dxz_b;
  abcd[818] = 4.0E0*I_NAI_I6z_Dxz_ab-2.0E0*5*I_NAI_G4z_Dxz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[819] = 4.0E0*I_NAI_I5xz_Dxy_ab;
  abcd[820] = 4.0E0*I_NAI_I4xyz_Dxy_ab;
  abcd[821] = 4.0E0*I_NAI_I4x2z_Dxy_ab-2.0E0*1*I_NAI_G4x_Dxy_b;
  abcd[822] = 4.0E0*I_NAI_I3x2yz_Dxy_ab;
  abcd[823] = 4.0E0*I_NAI_I3xy2z_Dxy_ab-2.0E0*1*I_NAI_G3xy_Dxy_b;
  abcd[824] = 4.0E0*I_NAI_I3x3z_Dxy_ab-2.0E0*2*I_NAI_G3xz_Dxy_b;
  abcd[825] = 4.0E0*I_NAI_I2x3yz_Dxy_ab;
  abcd[826] = 4.0E0*I_NAI_I2x2y2z_Dxy_ab-2.0E0*1*I_NAI_G2x2y_Dxy_b;
  abcd[827] = 4.0E0*I_NAI_I2xy3z_Dxy_ab-2.0E0*2*I_NAI_G2xyz_Dxy_b;
  abcd[828] = 4.0E0*I_NAI_I2x4z_Dxy_ab-2.0E0*3*I_NAI_G2x2z_Dxy_b;
  abcd[829] = 4.0E0*I_NAI_Ix4yz_Dxy_ab;
  abcd[830] = 4.0E0*I_NAI_Ix3y2z_Dxy_ab-2.0E0*1*I_NAI_Gx3y_Dxy_b;
  abcd[831] = 4.0E0*I_NAI_Ix2y3z_Dxy_ab-2.0E0*2*I_NAI_Gx2yz_Dxy_b;
  abcd[832] = 4.0E0*I_NAI_Ixy4z_Dxy_ab-2.0E0*3*I_NAI_Gxy2z_Dxy_b;
  abcd[833] = 4.0E0*I_NAI_Ix5z_Dxy_ab-2.0E0*4*I_NAI_Gx3z_Dxy_b;
  abcd[834] = 4.0E0*I_NAI_I5yz_Dxy_ab;
  abcd[835] = 4.0E0*I_NAI_I4y2z_Dxy_ab-2.0E0*1*I_NAI_G4y_Dxy_b;
  abcd[836] = 4.0E0*I_NAI_I3y3z_Dxy_ab-2.0E0*2*I_NAI_G3yz_Dxy_b;
  abcd[837] = 4.0E0*I_NAI_I2y4z_Dxy_ab-2.0E0*3*I_NAI_G2y2z_Dxy_b;
  abcd[838] = 4.0E0*I_NAI_Iy5z_Dxy_ab-2.0E0*4*I_NAI_Gy3z_Dxy_b;
  abcd[839] = 4.0E0*I_NAI_I6z_Dxy_ab-2.0E0*5*I_NAI_G4z_Dxy_b;
  abcd[840] = 4.0E0*I_NAI_I5xz_D2y_ab-2.0E0*1*I_NAI_I5xz_S_a;
  abcd[841] = 4.0E0*I_NAI_I4xyz_D2y_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[842] = 4.0E0*I_NAI_I4x2z_D2y_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*1*I_NAI_G4x_D2y_b+1*I_NAI_G4x_S;
  abcd[843] = 4.0E0*I_NAI_I3x2yz_D2y_ab-2.0E0*1*I_NAI_I3x2yz_S_a;
  abcd[844] = 4.0E0*I_NAI_I3xy2z_D2y_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*1*I_NAI_G3xy_D2y_b+1*I_NAI_G3xy_S;
  abcd[845] = 4.0E0*I_NAI_I3x3z_D2y_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_G3xz_D2y_b+2*1*I_NAI_G3xz_S;
  abcd[846] = 4.0E0*I_NAI_I2x3yz_D2y_ab-2.0E0*1*I_NAI_I2x3yz_S_a;
  abcd[847] = 4.0E0*I_NAI_I2x2y2z_D2y_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2y_D2y_b+1*I_NAI_G2x2y_S;
  abcd[848] = 4.0E0*I_NAI_I2xy3z_D2y_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*2*I_NAI_G2xyz_D2y_b+2*1*I_NAI_G2xyz_S;
  abcd[849] = 4.0E0*I_NAI_I2x4z_D2y_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*3*I_NAI_G2x2z_D2y_b+3*1*I_NAI_G2x2z_S;
  abcd[850] = 4.0E0*I_NAI_Ix4yz_D2y_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[851] = 4.0E0*I_NAI_Ix3y2z_D2y_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*1*I_NAI_Gx3y_D2y_b+1*I_NAI_Gx3y_S;
  abcd[852] = 4.0E0*I_NAI_Ix2y3z_D2y_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*2*I_NAI_Gx2yz_D2y_b+2*1*I_NAI_Gx2yz_S;
  abcd[853] = 4.0E0*I_NAI_Ixy4z_D2y_ab-2.0E0*1*I_NAI_Ixy4z_S_a-2.0E0*3*I_NAI_Gxy2z_D2y_b+3*1*I_NAI_Gxy2z_S;
  abcd[854] = 4.0E0*I_NAI_Ix5z_D2y_ab-2.0E0*1*I_NAI_Ix5z_S_a-2.0E0*4*I_NAI_Gx3z_D2y_b+4*1*I_NAI_Gx3z_S;
  abcd[855] = 4.0E0*I_NAI_I5yz_D2y_ab-2.0E0*1*I_NAI_I5yz_S_a;
  abcd[856] = 4.0E0*I_NAI_I4y2z_D2y_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*1*I_NAI_G4y_D2y_b+1*I_NAI_G4y_S;
  abcd[857] = 4.0E0*I_NAI_I3y3z_D2y_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_G3yz_D2y_b+2*1*I_NAI_G3yz_S;
  abcd[858] = 4.0E0*I_NAI_I2y4z_D2y_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*3*I_NAI_G2y2z_D2y_b+3*1*I_NAI_G2y2z_S;
  abcd[859] = 4.0E0*I_NAI_Iy5z_D2y_ab-2.0E0*1*I_NAI_Iy5z_S_a-2.0E0*4*I_NAI_Gy3z_D2y_b+4*1*I_NAI_Gy3z_S;
  abcd[860] = 4.0E0*I_NAI_I6z_D2y_ab-2.0E0*1*I_NAI_I6z_S_a-2.0E0*5*I_NAI_G4z_D2y_b+5*1*I_NAI_G4z_S;
  abcd[861] = 4.0E0*I_NAI_I5xz_Dyz_ab;
  abcd[862] = 4.0E0*I_NAI_I4xyz_Dyz_ab;
  abcd[863] = 4.0E0*I_NAI_I4x2z_Dyz_ab-2.0E0*1*I_NAI_G4x_Dyz_b;
  abcd[864] = 4.0E0*I_NAI_I3x2yz_Dyz_ab;
  abcd[865] = 4.0E0*I_NAI_I3xy2z_Dyz_ab-2.0E0*1*I_NAI_G3xy_Dyz_b;
  abcd[866] = 4.0E0*I_NAI_I3x3z_Dyz_ab-2.0E0*2*I_NAI_G3xz_Dyz_b;
  abcd[867] = 4.0E0*I_NAI_I2x3yz_Dyz_ab;
  abcd[868] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2x2y_Dyz_b;
  abcd[869] = 4.0E0*I_NAI_I2xy3z_Dyz_ab-2.0E0*2*I_NAI_G2xyz_Dyz_b;
  abcd[870] = 4.0E0*I_NAI_I2x4z_Dyz_ab-2.0E0*3*I_NAI_G2x2z_Dyz_b;
  abcd[871] = 4.0E0*I_NAI_Ix4yz_Dyz_ab;
  abcd[872] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab-2.0E0*1*I_NAI_Gx3y_Dyz_b;
  abcd[873] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab-2.0E0*2*I_NAI_Gx2yz_Dyz_b;
  abcd[874] = 4.0E0*I_NAI_Ixy4z_Dyz_ab-2.0E0*3*I_NAI_Gxy2z_Dyz_b;
  abcd[875] = 4.0E0*I_NAI_Ix5z_Dyz_ab-2.0E0*4*I_NAI_Gx3z_Dyz_b;
  abcd[876] = 4.0E0*I_NAI_I5yz_Dyz_ab;
  abcd[877] = 4.0E0*I_NAI_I4y2z_Dyz_ab-2.0E0*1*I_NAI_G4y_Dyz_b;
  abcd[878] = 4.0E0*I_NAI_I3y3z_Dyz_ab-2.0E0*2*I_NAI_G3yz_Dyz_b;
  abcd[879] = 4.0E0*I_NAI_I2y4z_Dyz_ab-2.0E0*3*I_NAI_G2y2z_Dyz_b;
  abcd[880] = 4.0E0*I_NAI_Iy5z_Dyz_ab-2.0E0*4*I_NAI_Gy3z_Dyz_b;
  abcd[881] = 4.0E0*I_NAI_I6z_Dyz_ab-2.0E0*5*I_NAI_G4z_Dyz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_ab
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  abcd[882] = 4.0E0*I_NAI_I5xz_Dxz_ab;
  abcd[883] = 4.0E0*I_NAI_I4xyz_Dxz_ab;
  abcd[884] = 4.0E0*I_NAI_I4x2z_Dxz_ab-2.0E0*1*I_NAI_G4x_Dxz_b;
  abcd[885] = 4.0E0*I_NAI_I3x2yz_Dxz_ab;
  abcd[886] = 4.0E0*I_NAI_I3xy2z_Dxz_ab-2.0E0*1*I_NAI_G3xy_Dxz_b;
  abcd[887] = 4.0E0*I_NAI_I3x3z_Dxz_ab-2.0E0*2*I_NAI_G3xz_Dxz_b;
  abcd[888] = 4.0E0*I_NAI_I2x3yz_Dxz_ab;
  abcd[889] = 4.0E0*I_NAI_I2x2y2z_Dxz_ab-2.0E0*1*I_NAI_G2x2y_Dxz_b;
  abcd[890] = 4.0E0*I_NAI_I2xy3z_Dxz_ab-2.0E0*2*I_NAI_G2xyz_Dxz_b;
  abcd[891] = 4.0E0*I_NAI_I2x4z_Dxz_ab-2.0E0*3*I_NAI_G2x2z_Dxz_b;
  abcd[892] = 4.0E0*I_NAI_Ix4yz_Dxz_ab;
  abcd[893] = 4.0E0*I_NAI_Ix3y2z_Dxz_ab-2.0E0*1*I_NAI_Gx3y_Dxz_b;
  abcd[894] = 4.0E0*I_NAI_Ix2y3z_Dxz_ab-2.0E0*2*I_NAI_Gx2yz_Dxz_b;
  abcd[895] = 4.0E0*I_NAI_Ixy4z_Dxz_ab-2.0E0*3*I_NAI_Gxy2z_Dxz_b;
  abcd[896] = 4.0E0*I_NAI_Ix5z_Dxz_ab-2.0E0*4*I_NAI_Gx3z_Dxz_b;
  abcd[897] = 4.0E0*I_NAI_I5yz_Dxz_ab;
  abcd[898] = 4.0E0*I_NAI_I4y2z_Dxz_ab-2.0E0*1*I_NAI_G4y_Dxz_b;
  abcd[899] = 4.0E0*I_NAI_I3y3z_Dxz_ab-2.0E0*2*I_NAI_G3yz_Dxz_b;
  abcd[900] = 4.0E0*I_NAI_I2y4z_Dxz_ab-2.0E0*3*I_NAI_G2y2z_Dxz_b;
  abcd[901] = 4.0E0*I_NAI_Iy5z_Dxz_ab-2.0E0*4*I_NAI_Gy3z_Dxz_b;
  abcd[902] = 4.0E0*I_NAI_I6z_Dxz_ab-2.0E0*5*I_NAI_G4z_Dxz_b;
  abcd[903] = 4.0E0*I_NAI_I5xz_Dyz_ab;
  abcd[904] = 4.0E0*I_NAI_I4xyz_Dyz_ab;
  abcd[905] = 4.0E0*I_NAI_I4x2z_Dyz_ab-2.0E0*1*I_NAI_G4x_Dyz_b;
  abcd[906] = 4.0E0*I_NAI_I3x2yz_Dyz_ab;
  abcd[907] = 4.0E0*I_NAI_I3xy2z_Dyz_ab-2.0E0*1*I_NAI_G3xy_Dyz_b;
  abcd[908] = 4.0E0*I_NAI_I3x3z_Dyz_ab-2.0E0*2*I_NAI_G3xz_Dyz_b;
  abcd[909] = 4.0E0*I_NAI_I2x3yz_Dyz_ab;
  abcd[910] = 4.0E0*I_NAI_I2x2y2z_Dyz_ab-2.0E0*1*I_NAI_G2x2y_Dyz_b;
  abcd[911] = 4.0E0*I_NAI_I2xy3z_Dyz_ab-2.0E0*2*I_NAI_G2xyz_Dyz_b;
  abcd[912] = 4.0E0*I_NAI_I2x4z_Dyz_ab-2.0E0*3*I_NAI_G2x2z_Dyz_b;
  abcd[913] = 4.0E0*I_NAI_Ix4yz_Dyz_ab;
  abcd[914] = 4.0E0*I_NAI_Ix3y2z_Dyz_ab-2.0E0*1*I_NAI_Gx3y_Dyz_b;
  abcd[915] = 4.0E0*I_NAI_Ix2y3z_Dyz_ab-2.0E0*2*I_NAI_Gx2yz_Dyz_b;
  abcd[916] = 4.0E0*I_NAI_Ixy4z_Dyz_ab-2.0E0*3*I_NAI_Gxy2z_Dyz_b;
  abcd[917] = 4.0E0*I_NAI_Ix5z_Dyz_ab-2.0E0*4*I_NAI_Gx3z_Dyz_b;
  abcd[918] = 4.0E0*I_NAI_I5yz_Dyz_ab;
  abcd[919] = 4.0E0*I_NAI_I4y2z_Dyz_ab-2.0E0*1*I_NAI_G4y_Dyz_b;
  abcd[920] = 4.0E0*I_NAI_I3y3z_Dyz_ab-2.0E0*2*I_NAI_G3yz_Dyz_b;
  abcd[921] = 4.0E0*I_NAI_I2y4z_Dyz_ab-2.0E0*3*I_NAI_G2y2z_Dyz_b;
  abcd[922] = 4.0E0*I_NAI_Iy5z_Dyz_ab-2.0E0*4*I_NAI_Gy3z_Dyz_b;
  abcd[923] = 4.0E0*I_NAI_I6z_Dyz_ab-2.0E0*5*I_NAI_G4z_Dyz_b;
  abcd[924] = 4.0E0*I_NAI_I5xz_D2z_ab-2.0E0*1*I_NAI_I5xz_S_a;
  abcd[925] = 4.0E0*I_NAI_I4xyz_D2z_ab-2.0E0*1*I_NAI_I4xyz_S_a;
  abcd[926] = 4.0E0*I_NAI_I4x2z_D2z_ab-2.0E0*1*I_NAI_I4x2z_S_a-2.0E0*1*I_NAI_G4x_D2z_b+1*I_NAI_G4x_S;
  abcd[927] = 4.0E0*I_NAI_I3x2yz_D2z_ab-2.0E0*1*I_NAI_I3x2yz_S_a;
  abcd[928] = 4.0E0*I_NAI_I3xy2z_D2z_ab-2.0E0*1*I_NAI_I3xy2z_S_a-2.0E0*1*I_NAI_G3xy_D2z_b+1*I_NAI_G3xy_S;
  abcd[929] = 4.0E0*I_NAI_I3x3z_D2z_ab-2.0E0*1*I_NAI_I3x3z_S_a-2.0E0*2*I_NAI_G3xz_D2z_b+2*1*I_NAI_G3xz_S;
  abcd[930] = 4.0E0*I_NAI_I2x3yz_D2z_ab-2.0E0*1*I_NAI_I2x3yz_S_a;
  abcd[931] = 4.0E0*I_NAI_I2x2y2z_D2z_ab-2.0E0*1*I_NAI_I2x2y2z_S_a-2.0E0*1*I_NAI_G2x2y_D2z_b+1*I_NAI_G2x2y_S;
  abcd[932] = 4.0E0*I_NAI_I2xy3z_D2z_ab-2.0E0*1*I_NAI_I2xy3z_S_a-2.0E0*2*I_NAI_G2xyz_D2z_b+2*1*I_NAI_G2xyz_S;
  abcd[933] = 4.0E0*I_NAI_I2x4z_D2z_ab-2.0E0*1*I_NAI_I2x4z_S_a-2.0E0*3*I_NAI_G2x2z_D2z_b+3*1*I_NAI_G2x2z_S;
  abcd[934] = 4.0E0*I_NAI_Ix4yz_D2z_ab-2.0E0*1*I_NAI_Ix4yz_S_a;
  abcd[935] = 4.0E0*I_NAI_Ix3y2z_D2z_ab-2.0E0*1*I_NAI_Ix3y2z_S_a-2.0E0*1*I_NAI_Gx3y_D2z_b+1*I_NAI_Gx3y_S;
  abcd[936] = 4.0E0*I_NAI_Ix2y3z_D2z_ab-2.0E0*1*I_NAI_Ix2y3z_S_a-2.0E0*2*I_NAI_Gx2yz_D2z_b+2*1*I_NAI_Gx2yz_S;
  abcd[937] = 4.0E0*I_NAI_Ixy4z_D2z_ab-2.0E0*1*I_NAI_Ixy4z_S_a-2.0E0*3*I_NAI_Gxy2z_D2z_b+3*1*I_NAI_Gxy2z_S;
  abcd[938] = 4.0E0*I_NAI_Ix5z_D2z_ab-2.0E0*1*I_NAI_Ix5z_S_a-2.0E0*4*I_NAI_Gx3z_D2z_b+4*1*I_NAI_Gx3z_S;
  abcd[939] = 4.0E0*I_NAI_I5yz_D2z_ab-2.0E0*1*I_NAI_I5yz_S_a;
  abcd[940] = 4.0E0*I_NAI_I4y2z_D2z_ab-2.0E0*1*I_NAI_I4y2z_S_a-2.0E0*1*I_NAI_G4y_D2z_b+1*I_NAI_G4y_S;
  abcd[941] = 4.0E0*I_NAI_I3y3z_D2z_ab-2.0E0*1*I_NAI_I3y3z_S_a-2.0E0*2*I_NAI_G3yz_D2z_b+2*1*I_NAI_G3yz_S;
  abcd[942] = 4.0E0*I_NAI_I2y4z_D2z_ab-2.0E0*1*I_NAI_I2y4z_S_a-2.0E0*3*I_NAI_G2y2z_D2z_b+3*1*I_NAI_G2y2z_S;
  abcd[943] = 4.0E0*I_NAI_Iy5z_D2z_ab-2.0E0*1*I_NAI_Iy5z_S_a-2.0E0*4*I_NAI_Gy3z_D2z_b+4*1*I_NAI_Gy3z_S;
  abcd[944] = 4.0E0*I_NAI_I6z_D2z_ab-2.0E0*1*I_NAI_I6z_S_a-2.0E0*5*I_NAI_G4z_D2z_b+5*1*I_NAI_G4z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[945] = 4.0E0*I_NAI_H5x_F3x_bb-2.0E0*1*I_NAI_H5x_Px_b-2.0E0*2*I_NAI_H5x_Px_b;
  abcd[946] = 4.0E0*I_NAI_H4xy_F3x_bb-2.0E0*1*I_NAI_H4xy_Px_b-2.0E0*2*I_NAI_H4xy_Px_b;
  abcd[947] = 4.0E0*I_NAI_H4xz_F3x_bb-2.0E0*1*I_NAI_H4xz_Px_b-2.0E0*2*I_NAI_H4xz_Px_b;
  abcd[948] = 4.0E0*I_NAI_H3x2y_F3x_bb-2.0E0*1*I_NAI_H3x2y_Px_b-2.0E0*2*I_NAI_H3x2y_Px_b;
  abcd[949] = 4.0E0*I_NAI_H3xyz_F3x_bb-2.0E0*1*I_NAI_H3xyz_Px_b-2.0E0*2*I_NAI_H3xyz_Px_b;
  abcd[950] = 4.0E0*I_NAI_H3x2z_F3x_bb-2.0E0*1*I_NAI_H3x2z_Px_b-2.0E0*2*I_NAI_H3x2z_Px_b;
  abcd[951] = 4.0E0*I_NAI_H2x3y_F3x_bb-2.0E0*1*I_NAI_H2x3y_Px_b-2.0E0*2*I_NAI_H2x3y_Px_b;
  abcd[952] = 4.0E0*I_NAI_H2x2yz_F3x_bb-2.0E0*1*I_NAI_H2x2yz_Px_b-2.0E0*2*I_NAI_H2x2yz_Px_b;
  abcd[953] = 4.0E0*I_NAI_H2xy2z_F3x_bb-2.0E0*1*I_NAI_H2xy2z_Px_b-2.0E0*2*I_NAI_H2xy2z_Px_b;
  abcd[954] = 4.0E0*I_NAI_H2x3z_F3x_bb-2.0E0*1*I_NAI_H2x3z_Px_b-2.0E0*2*I_NAI_H2x3z_Px_b;
  abcd[955] = 4.0E0*I_NAI_Hx4y_F3x_bb-2.0E0*1*I_NAI_Hx4y_Px_b-2.0E0*2*I_NAI_Hx4y_Px_b;
  abcd[956] = 4.0E0*I_NAI_Hx3yz_F3x_bb-2.0E0*1*I_NAI_Hx3yz_Px_b-2.0E0*2*I_NAI_Hx3yz_Px_b;
  abcd[957] = 4.0E0*I_NAI_Hx2y2z_F3x_bb-2.0E0*1*I_NAI_Hx2y2z_Px_b-2.0E0*2*I_NAI_Hx2y2z_Px_b;
  abcd[958] = 4.0E0*I_NAI_Hxy3z_F3x_bb-2.0E0*1*I_NAI_Hxy3z_Px_b-2.0E0*2*I_NAI_Hxy3z_Px_b;
  abcd[959] = 4.0E0*I_NAI_Hx4z_F3x_bb-2.0E0*1*I_NAI_Hx4z_Px_b-2.0E0*2*I_NAI_Hx4z_Px_b;
  abcd[960] = 4.0E0*I_NAI_H5y_F3x_bb-2.0E0*1*I_NAI_H5y_Px_b-2.0E0*2*I_NAI_H5y_Px_b;
  abcd[961] = 4.0E0*I_NAI_H4yz_F3x_bb-2.0E0*1*I_NAI_H4yz_Px_b-2.0E0*2*I_NAI_H4yz_Px_b;
  abcd[962] = 4.0E0*I_NAI_H3y2z_F3x_bb-2.0E0*1*I_NAI_H3y2z_Px_b-2.0E0*2*I_NAI_H3y2z_Px_b;
  abcd[963] = 4.0E0*I_NAI_H2y3z_F3x_bb-2.0E0*1*I_NAI_H2y3z_Px_b-2.0E0*2*I_NAI_H2y3z_Px_b;
  abcd[964] = 4.0E0*I_NAI_Hy4z_F3x_bb-2.0E0*1*I_NAI_Hy4z_Px_b-2.0E0*2*I_NAI_Hy4z_Px_b;
  abcd[965] = 4.0E0*I_NAI_H5z_F3x_bb-2.0E0*1*I_NAI_H5z_Px_b-2.0E0*2*I_NAI_H5z_Px_b;
  abcd[966] = 4.0E0*I_NAI_H5x_F2xy_bb-2.0E0*1*I_NAI_H5x_Py_b;
  abcd[967] = 4.0E0*I_NAI_H4xy_F2xy_bb-2.0E0*1*I_NAI_H4xy_Py_b;
  abcd[968] = 4.0E0*I_NAI_H4xz_F2xy_bb-2.0E0*1*I_NAI_H4xz_Py_b;
  abcd[969] = 4.0E0*I_NAI_H3x2y_F2xy_bb-2.0E0*1*I_NAI_H3x2y_Py_b;
  abcd[970] = 4.0E0*I_NAI_H3xyz_F2xy_bb-2.0E0*1*I_NAI_H3xyz_Py_b;
  abcd[971] = 4.0E0*I_NAI_H3x2z_F2xy_bb-2.0E0*1*I_NAI_H3x2z_Py_b;
  abcd[972] = 4.0E0*I_NAI_H2x3y_F2xy_bb-2.0E0*1*I_NAI_H2x3y_Py_b;
  abcd[973] = 4.0E0*I_NAI_H2x2yz_F2xy_bb-2.0E0*1*I_NAI_H2x2yz_Py_b;
  abcd[974] = 4.0E0*I_NAI_H2xy2z_F2xy_bb-2.0E0*1*I_NAI_H2xy2z_Py_b;
  abcd[975] = 4.0E0*I_NAI_H2x3z_F2xy_bb-2.0E0*1*I_NAI_H2x3z_Py_b;
  abcd[976] = 4.0E0*I_NAI_Hx4y_F2xy_bb-2.0E0*1*I_NAI_Hx4y_Py_b;
  abcd[977] = 4.0E0*I_NAI_Hx3yz_F2xy_bb-2.0E0*1*I_NAI_Hx3yz_Py_b;
  abcd[978] = 4.0E0*I_NAI_Hx2y2z_F2xy_bb-2.0E0*1*I_NAI_Hx2y2z_Py_b;
  abcd[979] = 4.0E0*I_NAI_Hxy3z_F2xy_bb-2.0E0*1*I_NAI_Hxy3z_Py_b;
  abcd[980] = 4.0E0*I_NAI_Hx4z_F2xy_bb-2.0E0*1*I_NAI_Hx4z_Py_b;
  abcd[981] = 4.0E0*I_NAI_H5y_F2xy_bb-2.0E0*1*I_NAI_H5y_Py_b;
  abcd[982] = 4.0E0*I_NAI_H4yz_F2xy_bb-2.0E0*1*I_NAI_H4yz_Py_b;
  abcd[983] = 4.0E0*I_NAI_H3y2z_F2xy_bb-2.0E0*1*I_NAI_H3y2z_Py_b;
  abcd[984] = 4.0E0*I_NAI_H2y3z_F2xy_bb-2.0E0*1*I_NAI_H2y3z_Py_b;
  abcd[985] = 4.0E0*I_NAI_Hy4z_F2xy_bb-2.0E0*1*I_NAI_Hy4z_Py_b;
  abcd[986] = 4.0E0*I_NAI_H5z_F2xy_bb-2.0E0*1*I_NAI_H5z_Py_b;
  abcd[987] = 4.0E0*I_NAI_H5x_F2xz_bb-2.0E0*1*I_NAI_H5x_Pz_b;
  abcd[988] = 4.0E0*I_NAI_H4xy_F2xz_bb-2.0E0*1*I_NAI_H4xy_Pz_b;
  abcd[989] = 4.0E0*I_NAI_H4xz_F2xz_bb-2.0E0*1*I_NAI_H4xz_Pz_b;
  abcd[990] = 4.0E0*I_NAI_H3x2y_F2xz_bb-2.0E0*1*I_NAI_H3x2y_Pz_b;
  abcd[991] = 4.0E0*I_NAI_H3xyz_F2xz_bb-2.0E0*1*I_NAI_H3xyz_Pz_b;
  abcd[992] = 4.0E0*I_NAI_H3x2z_F2xz_bb-2.0E0*1*I_NAI_H3x2z_Pz_b;
  abcd[993] = 4.0E0*I_NAI_H2x3y_F2xz_bb-2.0E0*1*I_NAI_H2x3y_Pz_b;
  abcd[994] = 4.0E0*I_NAI_H2x2yz_F2xz_bb-2.0E0*1*I_NAI_H2x2yz_Pz_b;
  abcd[995] = 4.0E0*I_NAI_H2xy2z_F2xz_bb-2.0E0*1*I_NAI_H2xy2z_Pz_b;
  abcd[996] = 4.0E0*I_NAI_H2x3z_F2xz_bb-2.0E0*1*I_NAI_H2x3z_Pz_b;
  abcd[997] = 4.0E0*I_NAI_Hx4y_F2xz_bb-2.0E0*1*I_NAI_Hx4y_Pz_b;
  abcd[998] = 4.0E0*I_NAI_Hx3yz_F2xz_bb-2.0E0*1*I_NAI_Hx3yz_Pz_b;
  abcd[999] = 4.0E0*I_NAI_Hx2y2z_F2xz_bb-2.0E0*1*I_NAI_Hx2y2z_Pz_b;
  abcd[1000] = 4.0E0*I_NAI_Hxy3z_F2xz_bb-2.0E0*1*I_NAI_Hxy3z_Pz_b;
  abcd[1001] = 4.0E0*I_NAI_Hx4z_F2xz_bb-2.0E0*1*I_NAI_Hx4z_Pz_b;
  abcd[1002] = 4.0E0*I_NAI_H5y_F2xz_bb-2.0E0*1*I_NAI_H5y_Pz_b;
  abcd[1003] = 4.0E0*I_NAI_H4yz_F2xz_bb-2.0E0*1*I_NAI_H4yz_Pz_b;
  abcd[1004] = 4.0E0*I_NAI_H3y2z_F2xz_bb-2.0E0*1*I_NAI_H3y2z_Pz_b;
  abcd[1005] = 4.0E0*I_NAI_H2y3z_F2xz_bb-2.0E0*1*I_NAI_H2y3z_Pz_b;
  abcd[1006] = 4.0E0*I_NAI_Hy4z_F2xz_bb-2.0E0*1*I_NAI_Hy4z_Pz_b;
  abcd[1007] = 4.0E0*I_NAI_H5z_F2xz_bb-2.0E0*1*I_NAI_H5z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[1008] = 4.0E0*I_NAI_H5x_F2xy_bb-2.0E0*1*I_NAI_H5x_Py_b;
  abcd[1009] = 4.0E0*I_NAI_H4xy_F2xy_bb-2.0E0*1*I_NAI_H4xy_Py_b;
  abcd[1010] = 4.0E0*I_NAI_H4xz_F2xy_bb-2.0E0*1*I_NAI_H4xz_Py_b;
  abcd[1011] = 4.0E0*I_NAI_H3x2y_F2xy_bb-2.0E0*1*I_NAI_H3x2y_Py_b;
  abcd[1012] = 4.0E0*I_NAI_H3xyz_F2xy_bb-2.0E0*1*I_NAI_H3xyz_Py_b;
  abcd[1013] = 4.0E0*I_NAI_H3x2z_F2xy_bb-2.0E0*1*I_NAI_H3x2z_Py_b;
  abcd[1014] = 4.0E0*I_NAI_H2x3y_F2xy_bb-2.0E0*1*I_NAI_H2x3y_Py_b;
  abcd[1015] = 4.0E0*I_NAI_H2x2yz_F2xy_bb-2.0E0*1*I_NAI_H2x2yz_Py_b;
  abcd[1016] = 4.0E0*I_NAI_H2xy2z_F2xy_bb-2.0E0*1*I_NAI_H2xy2z_Py_b;
  abcd[1017] = 4.0E0*I_NAI_H2x3z_F2xy_bb-2.0E0*1*I_NAI_H2x3z_Py_b;
  abcd[1018] = 4.0E0*I_NAI_Hx4y_F2xy_bb-2.0E0*1*I_NAI_Hx4y_Py_b;
  abcd[1019] = 4.0E0*I_NAI_Hx3yz_F2xy_bb-2.0E0*1*I_NAI_Hx3yz_Py_b;
  abcd[1020] = 4.0E0*I_NAI_Hx2y2z_F2xy_bb-2.0E0*1*I_NAI_Hx2y2z_Py_b;
  abcd[1021] = 4.0E0*I_NAI_Hxy3z_F2xy_bb-2.0E0*1*I_NAI_Hxy3z_Py_b;
  abcd[1022] = 4.0E0*I_NAI_Hx4z_F2xy_bb-2.0E0*1*I_NAI_Hx4z_Py_b;
  abcd[1023] = 4.0E0*I_NAI_H5y_F2xy_bb-2.0E0*1*I_NAI_H5y_Py_b;
  abcd[1024] = 4.0E0*I_NAI_H4yz_F2xy_bb-2.0E0*1*I_NAI_H4yz_Py_b;
  abcd[1025] = 4.0E0*I_NAI_H3y2z_F2xy_bb-2.0E0*1*I_NAI_H3y2z_Py_b;
  abcd[1026] = 4.0E0*I_NAI_H2y3z_F2xy_bb-2.0E0*1*I_NAI_H2y3z_Py_b;
  abcd[1027] = 4.0E0*I_NAI_Hy4z_F2xy_bb-2.0E0*1*I_NAI_Hy4z_Py_b;
  abcd[1028] = 4.0E0*I_NAI_H5z_F2xy_bb-2.0E0*1*I_NAI_H5z_Py_b;
  abcd[1029] = 4.0E0*I_NAI_H5x_Fx2y_bb-2.0E0*1*I_NAI_H5x_Px_b;
  abcd[1030] = 4.0E0*I_NAI_H4xy_Fx2y_bb-2.0E0*1*I_NAI_H4xy_Px_b;
  abcd[1031] = 4.0E0*I_NAI_H4xz_Fx2y_bb-2.0E0*1*I_NAI_H4xz_Px_b;
  abcd[1032] = 4.0E0*I_NAI_H3x2y_Fx2y_bb-2.0E0*1*I_NAI_H3x2y_Px_b;
  abcd[1033] = 4.0E0*I_NAI_H3xyz_Fx2y_bb-2.0E0*1*I_NAI_H3xyz_Px_b;
  abcd[1034] = 4.0E0*I_NAI_H3x2z_Fx2y_bb-2.0E0*1*I_NAI_H3x2z_Px_b;
  abcd[1035] = 4.0E0*I_NAI_H2x3y_Fx2y_bb-2.0E0*1*I_NAI_H2x3y_Px_b;
  abcd[1036] = 4.0E0*I_NAI_H2x2yz_Fx2y_bb-2.0E0*1*I_NAI_H2x2yz_Px_b;
  abcd[1037] = 4.0E0*I_NAI_H2xy2z_Fx2y_bb-2.0E0*1*I_NAI_H2xy2z_Px_b;
  abcd[1038] = 4.0E0*I_NAI_H2x3z_Fx2y_bb-2.0E0*1*I_NAI_H2x3z_Px_b;
  abcd[1039] = 4.0E0*I_NAI_Hx4y_Fx2y_bb-2.0E0*1*I_NAI_Hx4y_Px_b;
  abcd[1040] = 4.0E0*I_NAI_Hx3yz_Fx2y_bb-2.0E0*1*I_NAI_Hx3yz_Px_b;
  abcd[1041] = 4.0E0*I_NAI_Hx2y2z_Fx2y_bb-2.0E0*1*I_NAI_Hx2y2z_Px_b;
  abcd[1042] = 4.0E0*I_NAI_Hxy3z_Fx2y_bb-2.0E0*1*I_NAI_Hxy3z_Px_b;
  abcd[1043] = 4.0E0*I_NAI_Hx4z_Fx2y_bb-2.0E0*1*I_NAI_Hx4z_Px_b;
  abcd[1044] = 4.0E0*I_NAI_H5y_Fx2y_bb-2.0E0*1*I_NAI_H5y_Px_b;
  abcd[1045] = 4.0E0*I_NAI_H4yz_Fx2y_bb-2.0E0*1*I_NAI_H4yz_Px_b;
  abcd[1046] = 4.0E0*I_NAI_H3y2z_Fx2y_bb-2.0E0*1*I_NAI_H3y2z_Px_b;
  abcd[1047] = 4.0E0*I_NAI_H2y3z_Fx2y_bb-2.0E0*1*I_NAI_H2y3z_Px_b;
  abcd[1048] = 4.0E0*I_NAI_Hy4z_Fx2y_bb-2.0E0*1*I_NAI_Hy4z_Px_b;
  abcd[1049] = 4.0E0*I_NAI_H5z_Fx2y_bb-2.0E0*1*I_NAI_H5z_Px_b;
  abcd[1050] = 4.0E0*I_NAI_H5x_Fxyz_bb;
  abcd[1051] = 4.0E0*I_NAI_H4xy_Fxyz_bb;
  abcd[1052] = 4.0E0*I_NAI_H4xz_Fxyz_bb;
  abcd[1053] = 4.0E0*I_NAI_H3x2y_Fxyz_bb;
  abcd[1054] = 4.0E0*I_NAI_H3xyz_Fxyz_bb;
  abcd[1055] = 4.0E0*I_NAI_H3x2z_Fxyz_bb;
  abcd[1056] = 4.0E0*I_NAI_H2x3y_Fxyz_bb;
  abcd[1057] = 4.0E0*I_NAI_H2x2yz_Fxyz_bb;
  abcd[1058] = 4.0E0*I_NAI_H2xy2z_Fxyz_bb;
  abcd[1059] = 4.0E0*I_NAI_H2x3z_Fxyz_bb;
  abcd[1060] = 4.0E0*I_NAI_Hx4y_Fxyz_bb;
  abcd[1061] = 4.0E0*I_NAI_Hx3yz_Fxyz_bb;
  abcd[1062] = 4.0E0*I_NAI_Hx2y2z_Fxyz_bb;
  abcd[1063] = 4.0E0*I_NAI_Hxy3z_Fxyz_bb;
  abcd[1064] = 4.0E0*I_NAI_Hx4z_Fxyz_bb;
  abcd[1065] = 4.0E0*I_NAI_H5y_Fxyz_bb;
  abcd[1066] = 4.0E0*I_NAI_H4yz_Fxyz_bb;
  abcd[1067] = 4.0E0*I_NAI_H3y2z_Fxyz_bb;
  abcd[1068] = 4.0E0*I_NAI_H2y3z_Fxyz_bb;
  abcd[1069] = 4.0E0*I_NAI_Hy4z_Fxyz_bb;
  abcd[1070] = 4.0E0*I_NAI_H5z_Fxyz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[1071] = 4.0E0*I_NAI_H5x_F2xz_bb-2.0E0*1*I_NAI_H5x_Pz_b;
  abcd[1072] = 4.0E0*I_NAI_H4xy_F2xz_bb-2.0E0*1*I_NAI_H4xy_Pz_b;
  abcd[1073] = 4.0E0*I_NAI_H4xz_F2xz_bb-2.0E0*1*I_NAI_H4xz_Pz_b;
  abcd[1074] = 4.0E0*I_NAI_H3x2y_F2xz_bb-2.0E0*1*I_NAI_H3x2y_Pz_b;
  abcd[1075] = 4.0E0*I_NAI_H3xyz_F2xz_bb-2.0E0*1*I_NAI_H3xyz_Pz_b;
  abcd[1076] = 4.0E0*I_NAI_H3x2z_F2xz_bb-2.0E0*1*I_NAI_H3x2z_Pz_b;
  abcd[1077] = 4.0E0*I_NAI_H2x3y_F2xz_bb-2.0E0*1*I_NAI_H2x3y_Pz_b;
  abcd[1078] = 4.0E0*I_NAI_H2x2yz_F2xz_bb-2.0E0*1*I_NAI_H2x2yz_Pz_b;
  abcd[1079] = 4.0E0*I_NAI_H2xy2z_F2xz_bb-2.0E0*1*I_NAI_H2xy2z_Pz_b;
  abcd[1080] = 4.0E0*I_NAI_H2x3z_F2xz_bb-2.0E0*1*I_NAI_H2x3z_Pz_b;
  abcd[1081] = 4.0E0*I_NAI_Hx4y_F2xz_bb-2.0E0*1*I_NAI_Hx4y_Pz_b;
  abcd[1082] = 4.0E0*I_NAI_Hx3yz_F2xz_bb-2.0E0*1*I_NAI_Hx3yz_Pz_b;
  abcd[1083] = 4.0E0*I_NAI_Hx2y2z_F2xz_bb-2.0E0*1*I_NAI_Hx2y2z_Pz_b;
  abcd[1084] = 4.0E0*I_NAI_Hxy3z_F2xz_bb-2.0E0*1*I_NAI_Hxy3z_Pz_b;
  abcd[1085] = 4.0E0*I_NAI_Hx4z_F2xz_bb-2.0E0*1*I_NAI_Hx4z_Pz_b;
  abcd[1086] = 4.0E0*I_NAI_H5y_F2xz_bb-2.0E0*1*I_NAI_H5y_Pz_b;
  abcd[1087] = 4.0E0*I_NAI_H4yz_F2xz_bb-2.0E0*1*I_NAI_H4yz_Pz_b;
  abcd[1088] = 4.0E0*I_NAI_H3y2z_F2xz_bb-2.0E0*1*I_NAI_H3y2z_Pz_b;
  abcd[1089] = 4.0E0*I_NAI_H2y3z_F2xz_bb-2.0E0*1*I_NAI_H2y3z_Pz_b;
  abcd[1090] = 4.0E0*I_NAI_Hy4z_F2xz_bb-2.0E0*1*I_NAI_Hy4z_Pz_b;
  abcd[1091] = 4.0E0*I_NAI_H5z_F2xz_bb-2.0E0*1*I_NAI_H5z_Pz_b;
  abcd[1092] = 4.0E0*I_NAI_H5x_Fxyz_bb;
  abcd[1093] = 4.0E0*I_NAI_H4xy_Fxyz_bb;
  abcd[1094] = 4.0E0*I_NAI_H4xz_Fxyz_bb;
  abcd[1095] = 4.0E0*I_NAI_H3x2y_Fxyz_bb;
  abcd[1096] = 4.0E0*I_NAI_H3xyz_Fxyz_bb;
  abcd[1097] = 4.0E0*I_NAI_H3x2z_Fxyz_bb;
  abcd[1098] = 4.0E0*I_NAI_H2x3y_Fxyz_bb;
  abcd[1099] = 4.0E0*I_NAI_H2x2yz_Fxyz_bb;
  abcd[1100] = 4.0E0*I_NAI_H2xy2z_Fxyz_bb;
  abcd[1101] = 4.0E0*I_NAI_H2x3z_Fxyz_bb;
  abcd[1102] = 4.0E0*I_NAI_Hx4y_Fxyz_bb;
  abcd[1103] = 4.0E0*I_NAI_Hx3yz_Fxyz_bb;
  abcd[1104] = 4.0E0*I_NAI_Hx2y2z_Fxyz_bb;
  abcd[1105] = 4.0E0*I_NAI_Hxy3z_Fxyz_bb;
  abcd[1106] = 4.0E0*I_NAI_Hx4z_Fxyz_bb;
  abcd[1107] = 4.0E0*I_NAI_H5y_Fxyz_bb;
  abcd[1108] = 4.0E0*I_NAI_H4yz_Fxyz_bb;
  abcd[1109] = 4.0E0*I_NAI_H3y2z_Fxyz_bb;
  abcd[1110] = 4.0E0*I_NAI_H2y3z_Fxyz_bb;
  abcd[1111] = 4.0E0*I_NAI_Hy4z_Fxyz_bb;
  abcd[1112] = 4.0E0*I_NAI_H5z_Fxyz_bb;
  abcd[1113] = 4.0E0*I_NAI_H5x_Fx2z_bb-2.0E0*1*I_NAI_H5x_Px_b;
  abcd[1114] = 4.0E0*I_NAI_H4xy_Fx2z_bb-2.0E0*1*I_NAI_H4xy_Px_b;
  abcd[1115] = 4.0E0*I_NAI_H4xz_Fx2z_bb-2.0E0*1*I_NAI_H4xz_Px_b;
  abcd[1116] = 4.0E0*I_NAI_H3x2y_Fx2z_bb-2.0E0*1*I_NAI_H3x2y_Px_b;
  abcd[1117] = 4.0E0*I_NAI_H3xyz_Fx2z_bb-2.0E0*1*I_NAI_H3xyz_Px_b;
  abcd[1118] = 4.0E0*I_NAI_H3x2z_Fx2z_bb-2.0E0*1*I_NAI_H3x2z_Px_b;
  abcd[1119] = 4.0E0*I_NAI_H2x3y_Fx2z_bb-2.0E0*1*I_NAI_H2x3y_Px_b;
  abcd[1120] = 4.0E0*I_NAI_H2x2yz_Fx2z_bb-2.0E0*1*I_NAI_H2x2yz_Px_b;
  abcd[1121] = 4.0E0*I_NAI_H2xy2z_Fx2z_bb-2.0E0*1*I_NAI_H2xy2z_Px_b;
  abcd[1122] = 4.0E0*I_NAI_H2x3z_Fx2z_bb-2.0E0*1*I_NAI_H2x3z_Px_b;
  abcd[1123] = 4.0E0*I_NAI_Hx4y_Fx2z_bb-2.0E0*1*I_NAI_Hx4y_Px_b;
  abcd[1124] = 4.0E0*I_NAI_Hx3yz_Fx2z_bb-2.0E0*1*I_NAI_Hx3yz_Px_b;
  abcd[1125] = 4.0E0*I_NAI_Hx2y2z_Fx2z_bb-2.0E0*1*I_NAI_Hx2y2z_Px_b;
  abcd[1126] = 4.0E0*I_NAI_Hxy3z_Fx2z_bb-2.0E0*1*I_NAI_Hxy3z_Px_b;
  abcd[1127] = 4.0E0*I_NAI_Hx4z_Fx2z_bb-2.0E0*1*I_NAI_Hx4z_Px_b;
  abcd[1128] = 4.0E0*I_NAI_H5y_Fx2z_bb-2.0E0*1*I_NAI_H5y_Px_b;
  abcd[1129] = 4.0E0*I_NAI_H4yz_Fx2z_bb-2.0E0*1*I_NAI_H4yz_Px_b;
  abcd[1130] = 4.0E0*I_NAI_H3y2z_Fx2z_bb-2.0E0*1*I_NAI_H3y2z_Px_b;
  abcd[1131] = 4.0E0*I_NAI_H2y3z_Fx2z_bb-2.0E0*1*I_NAI_H2y3z_Px_b;
  abcd[1132] = 4.0E0*I_NAI_Hy4z_Fx2z_bb-2.0E0*1*I_NAI_Hy4z_Px_b;
  abcd[1133] = 4.0E0*I_NAI_H5z_Fx2z_bb-2.0E0*1*I_NAI_H5z_Px_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[1134] = 4.0E0*I_NAI_H5x_Fx2y_bb-2.0E0*1*I_NAI_H5x_Px_b;
  abcd[1135] = 4.0E0*I_NAI_H4xy_Fx2y_bb-2.0E0*1*I_NAI_H4xy_Px_b;
  abcd[1136] = 4.0E0*I_NAI_H4xz_Fx2y_bb-2.0E0*1*I_NAI_H4xz_Px_b;
  abcd[1137] = 4.0E0*I_NAI_H3x2y_Fx2y_bb-2.0E0*1*I_NAI_H3x2y_Px_b;
  abcd[1138] = 4.0E0*I_NAI_H3xyz_Fx2y_bb-2.0E0*1*I_NAI_H3xyz_Px_b;
  abcd[1139] = 4.0E0*I_NAI_H3x2z_Fx2y_bb-2.0E0*1*I_NAI_H3x2z_Px_b;
  abcd[1140] = 4.0E0*I_NAI_H2x3y_Fx2y_bb-2.0E0*1*I_NAI_H2x3y_Px_b;
  abcd[1141] = 4.0E0*I_NAI_H2x2yz_Fx2y_bb-2.0E0*1*I_NAI_H2x2yz_Px_b;
  abcd[1142] = 4.0E0*I_NAI_H2xy2z_Fx2y_bb-2.0E0*1*I_NAI_H2xy2z_Px_b;
  abcd[1143] = 4.0E0*I_NAI_H2x3z_Fx2y_bb-2.0E0*1*I_NAI_H2x3z_Px_b;
  abcd[1144] = 4.0E0*I_NAI_Hx4y_Fx2y_bb-2.0E0*1*I_NAI_Hx4y_Px_b;
  abcd[1145] = 4.0E0*I_NAI_Hx3yz_Fx2y_bb-2.0E0*1*I_NAI_Hx3yz_Px_b;
  abcd[1146] = 4.0E0*I_NAI_Hx2y2z_Fx2y_bb-2.0E0*1*I_NAI_Hx2y2z_Px_b;
  abcd[1147] = 4.0E0*I_NAI_Hxy3z_Fx2y_bb-2.0E0*1*I_NAI_Hxy3z_Px_b;
  abcd[1148] = 4.0E0*I_NAI_Hx4z_Fx2y_bb-2.0E0*1*I_NAI_Hx4z_Px_b;
  abcd[1149] = 4.0E0*I_NAI_H5y_Fx2y_bb-2.0E0*1*I_NAI_H5y_Px_b;
  abcd[1150] = 4.0E0*I_NAI_H4yz_Fx2y_bb-2.0E0*1*I_NAI_H4yz_Px_b;
  abcd[1151] = 4.0E0*I_NAI_H3y2z_Fx2y_bb-2.0E0*1*I_NAI_H3y2z_Px_b;
  abcd[1152] = 4.0E0*I_NAI_H2y3z_Fx2y_bb-2.0E0*1*I_NAI_H2y3z_Px_b;
  abcd[1153] = 4.0E0*I_NAI_Hy4z_Fx2y_bb-2.0E0*1*I_NAI_Hy4z_Px_b;
  abcd[1154] = 4.0E0*I_NAI_H5z_Fx2y_bb-2.0E0*1*I_NAI_H5z_Px_b;
  abcd[1155] = 4.0E0*I_NAI_H5x_F3y_bb-2.0E0*1*I_NAI_H5x_Py_b-2.0E0*2*I_NAI_H5x_Py_b;
  abcd[1156] = 4.0E0*I_NAI_H4xy_F3y_bb-2.0E0*1*I_NAI_H4xy_Py_b-2.0E0*2*I_NAI_H4xy_Py_b;
  abcd[1157] = 4.0E0*I_NAI_H4xz_F3y_bb-2.0E0*1*I_NAI_H4xz_Py_b-2.0E0*2*I_NAI_H4xz_Py_b;
  abcd[1158] = 4.0E0*I_NAI_H3x2y_F3y_bb-2.0E0*1*I_NAI_H3x2y_Py_b-2.0E0*2*I_NAI_H3x2y_Py_b;
  abcd[1159] = 4.0E0*I_NAI_H3xyz_F3y_bb-2.0E0*1*I_NAI_H3xyz_Py_b-2.0E0*2*I_NAI_H3xyz_Py_b;
  abcd[1160] = 4.0E0*I_NAI_H3x2z_F3y_bb-2.0E0*1*I_NAI_H3x2z_Py_b-2.0E0*2*I_NAI_H3x2z_Py_b;
  abcd[1161] = 4.0E0*I_NAI_H2x3y_F3y_bb-2.0E0*1*I_NAI_H2x3y_Py_b-2.0E0*2*I_NAI_H2x3y_Py_b;
  abcd[1162] = 4.0E0*I_NAI_H2x2yz_F3y_bb-2.0E0*1*I_NAI_H2x2yz_Py_b-2.0E0*2*I_NAI_H2x2yz_Py_b;
  abcd[1163] = 4.0E0*I_NAI_H2xy2z_F3y_bb-2.0E0*1*I_NAI_H2xy2z_Py_b-2.0E0*2*I_NAI_H2xy2z_Py_b;
  abcd[1164] = 4.0E0*I_NAI_H2x3z_F3y_bb-2.0E0*1*I_NAI_H2x3z_Py_b-2.0E0*2*I_NAI_H2x3z_Py_b;
  abcd[1165] = 4.0E0*I_NAI_Hx4y_F3y_bb-2.0E0*1*I_NAI_Hx4y_Py_b-2.0E0*2*I_NAI_Hx4y_Py_b;
  abcd[1166] = 4.0E0*I_NAI_Hx3yz_F3y_bb-2.0E0*1*I_NAI_Hx3yz_Py_b-2.0E0*2*I_NAI_Hx3yz_Py_b;
  abcd[1167] = 4.0E0*I_NAI_Hx2y2z_F3y_bb-2.0E0*1*I_NAI_Hx2y2z_Py_b-2.0E0*2*I_NAI_Hx2y2z_Py_b;
  abcd[1168] = 4.0E0*I_NAI_Hxy3z_F3y_bb-2.0E0*1*I_NAI_Hxy3z_Py_b-2.0E0*2*I_NAI_Hxy3z_Py_b;
  abcd[1169] = 4.0E0*I_NAI_Hx4z_F3y_bb-2.0E0*1*I_NAI_Hx4z_Py_b-2.0E0*2*I_NAI_Hx4z_Py_b;
  abcd[1170] = 4.0E0*I_NAI_H5y_F3y_bb-2.0E0*1*I_NAI_H5y_Py_b-2.0E0*2*I_NAI_H5y_Py_b;
  abcd[1171] = 4.0E0*I_NAI_H4yz_F3y_bb-2.0E0*1*I_NAI_H4yz_Py_b-2.0E0*2*I_NAI_H4yz_Py_b;
  abcd[1172] = 4.0E0*I_NAI_H3y2z_F3y_bb-2.0E0*1*I_NAI_H3y2z_Py_b-2.0E0*2*I_NAI_H3y2z_Py_b;
  abcd[1173] = 4.0E0*I_NAI_H2y3z_F3y_bb-2.0E0*1*I_NAI_H2y3z_Py_b-2.0E0*2*I_NAI_H2y3z_Py_b;
  abcd[1174] = 4.0E0*I_NAI_Hy4z_F3y_bb-2.0E0*1*I_NAI_Hy4z_Py_b-2.0E0*2*I_NAI_Hy4z_Py_b;
  abcd[1175] = 4.0E0*I_NAI_H5z_F3y_bb-2.0E0*1*I_NAI_H5z_Py_b-2.0E0*2*I_NAI_H5z_Py_b;
  abcd[1176] = 4.0E0*I_NAI_H5x_F2yz_bb-2.0E0*1*I_NAI_H5x_Pz_b;
  abcd[1177] = 4.0E0*I_NAI_H4xy_F2yz_bb-2.0E0*1*I_NAI_H4xy_Pz_b;
  abcd[1178] = 4.0E0*I_NAI_H4xz_F2yz_bb-2.0E0*1*I_NAI_H4xz_Pz_b;
  abcd[1179] = 4.0E0*I_NAI_H3x2y_F2yz_bb-2.0E0*1*I_NAI_H3x2y_Pz_b;
  abcd[1180] = 4.0E0*I_NAI_H3xyz_F2yz_bb-2.0E0*1*I_NAI_H3xyz_Pz_b;
  abcd[1181] = 4.0E0*I_NAI_H3x2z_F2yz_bb-2.0E0*1*I_NAI_H3x2z_Pz_b;
  abcd[1182] = 4.0E0*I_NAI_H2x3y_F2yz_bb-2.0E0*1*I_NAI_H2x3y_Pz_b;
  abcd[1183] = 4.0E0*I_NAI_H2x2yz_F2yz_bb-2.0E0*1*I_NAI_H2x2yz_Pz_b;
  abcd[1184] = 4.0E0*I_NAI_H2xy2z_F2yz_bb-2.0E0*1*I_NAI_H2xy2z_Pz_b;
  abcd[1185] = 4.0E0*I_NAI_H2x3z_F2yz_bb-2.0E0*1*I_NAI_H2x3z_Pz_b;
  abcd[1186] = 4.0E0*I_NAI_Hx4y_F2yz_bb-2.0E0*1*I_NAI_Hx4y_Pz_b;
  abcd[1187] = 4.0E0*I_NAI_Hx3yz_F2yz_bb-2.0E0*1*I_NAI_Hx3yz_Pz_b;
  abcd[1188] = 4.0E0*I_NAI_Hx2y2z_F2yz_bb-2.0E0*1*I_NAI_Hx2y2z_Pz_b;
  abcd[1189] = 4.0E0*I_NAI_Hxy3z_F2yz_bb-2.0E0*1*I_NAI_Hxy3z_Pz_b;
  abcd[1190] = 4.0E0*I_NAI_Hx4z_F2yz_bb-2.0E0*1*I_NAI_Hx4z_Pz_b;
  abcd[1191] = 4.0E0*I_NAI_H5y_F2yz_bb-2.0E0*1*I_NAI_H5y_Pz_b;
  abcd[1192] = 4.0E0*I_NAI_H4yz_F2yz_bb-2.0E0*1*I_NAI_H4yz_Pz_b;
  abcd[1193] = 4.0E0*I_NAI_H3y2z_F2yz_bb-2.0E0*1*I_NAI_H3y2z_Pz_b;
  abcd[1194] = 4.0E0*I_NAI_H2y3z_F2yz_bb-2.0E0*1*I_NAI_H2y3z_Pz_b;
  abcd[1195] = 4.0E0*I_NAI_Hy4z_F2yz_bb-2.0E0*1*I_NAI_Hy4z_Pz_b;
  abcd[1196] = 4.0E0*I_NAI_H5z_F2yz_bb-2.0E0*1*I_NAI_H5z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[1197] = 4.0E0*I_NAI_H5x_Fxyz_bb;
  abcd[1198] = 4.0E0*I_NAI_H4xy_Fxyz_bb;
  abcd[1199] = 4.0E0*I_NAI_H4xz_Fxyz_bb;
  abcd[1200] = 4.0E0*I_NAI_H3x2y_Fxyz_bb;
  abcd[1201] = 4.0E0*I_NAI_H3xyz_Fxyz_bb;
  abcd[1202] = 4.0E0*I_NAI_H3x2z_Fxyz_bb;
  abcd[1203] = 4.0E0*I_NAI_H2x3y_Fxyz_bb;
  abcd[1204] = 4.0E0*I_NAI_H2x2yz_Fxyz_bb;
  abcd[1205] = 4.0E0*I_NAI_H2xy2z_Fxyz_bb;
  abcd[1206] = 4.0E0*I_NAI_H2x3z_Fxyz_bb;
  abcd[1207] = 4.0E0*I_NAI_Hx4y_Fxyz_bb;
  abcd[1208] = 4.0E0*I_NAI_Hx3yz_Fxyz_bb;
  abcd[1209] = 4.0E0*I_NAI_Hx2y2z_Fxyz_bb;
  abcd[1210] = 4.0E0*I_NAI_Hxy3z_Fxyz_bb;
  abcd[1211] = 4.0E0*I_NAI_Hx4z_Fxyz_bb;
  abcd[1212] = 4.0E0*I_NAI_H5y_Fxyz_bb;
  abcd[1213] = 4.0E0*I_NAI_H4yz_Fxyz_bb;
  abcd[1214] = 4.0E0*I_NAI_H3y2z_Fxyz_bb;
  abcd[1215] = 4.0E0*I_NAI_H2y3z_Fxyz_bb;
  abcd[1216] = 4.0E0*I_NAI_Hy4z_Fxyz_bb;
  abcd[1217] = 4.0E0*I_NAI_H5z_Fxyz_bb;
  abcd[1218] = 4.0E0*I_NAI_H5x_F2yz_bb-2.0E0*1*I_NAI_H5x_Pz_b;
  abcd[1219] = 4.0E0*I_NAI_H4xy_F2yz_bb-2.0E0*1*I_NAI_H4xy_Pz_b;
  abcd[1220] = 4.0E0*I_NAI_H4xz_F2yz_bb-2.0E0*1*I_NAI_H4xz_Pz_b;
  abcd[1221] = 4.0E0*I_NAI_H3x2y_F2yz_bb-2.0E0*1*I_NAI_H3x2y_Pz_b;
  abcd[1222] = 4.0E0*I_NAI_H3xyz_F2yz_bb-2.0E0*1*I_NAI_H3xyz_Pz_b;
  abcd[1223] = 4.0E0*I_NAI_H3x2z_F2yz_bb-2.0E0*1*I_NAI_H3x2z_Pz_b;
  abcd[1224] = 4.0E0*I_NAI_H2x3y_F2yz_bb-2.0E0*1*I_NAI_H2x3y_Pz_b;
  abcd[1225] = 4.0E0*I_NAI_H2x2yz_F2yz_bb-2.0E0*1*I_NAI_H2x2yz_Pz_b;
  abcd[1226] = 4.0E0*I_NAI_H2xy2z_F2yz_bb-2.0E0*1*I_NAI_H2xy2z_Pz_b;
  abcd[1227] = 4.0E0*I_NAI_H2x3z_F2yz_bb-2.0E0*1*I_NAI_H2x3z_Pz_b;
  abcd[1228] = 4.0E0*I_NAI_Hx4y_F2yz_bb-2.0E0*1*I_NAI_Hx4y_Pz_b;
  abcd[1229] = 4.0E0*I_NAI_Hx3yz_F2yz_bb-2.0E0*1*I_NAI_Hx3yz_Pz_b;
  abcd[1230] = 4.0E0*I_NAI_Hx2y2z_F2yz_bb-2.0E0*1*I_NAI_Hx2y2z_Pz_b;
  abcd[1231] = 4.0E0*I_NAI_Hxy3z_F2yz_bb-2.0E0*1*I_NAI_Hxy3z_Pz_b;
  abcd[1232] = 4.0E0*I_NAI_Hx4z_F2yz_bb-2.0E0*1*I_NAI_Hx4z_Pz_b;
  abcd[1233] = 4.0E0*I_NAI_H5y_F2yz_bb-2.0E0*1*I_NAI_H5y_Pz_b;
  abcd[1234] = 4.0E0*I_NAI_H4yz_F2yz_bb-2.0E0*1*I_NAI_H4yz_Pz_b;
  abcd[1235] = 4.0E0*I_NAI_H3y2z_F2yz_bb-2.0E0*1*I_NAI_H3y2z_Pz_b;
  abcd[1236] = 4.0E0*I_NAI_H2y3z_F2yz_bb-2.0E0*1*I_NAI_H2y3z_Pz_b;
  abcd[1237] = 4.0E0*I_NAI_Hy4z_F2yz_bb-2.0E0*1*I_NAI_Hy4z_Pz_b;
  abcd[1238] = 4.0E0*I_NAI_H5z_F2yz_bb-2.0E0*1*I_NAI_H5z_Pz_b;
  abcd[1239] = 4.0E0*I_NAI_H5x_Fy2z_bb-2.0E0*1*I_NAI_H5x_Py_b;
  abcd[1240] = 4.0E0*I_NAI_H4xy_Fy2z_bb-2.0E0*1*I_NAI_H4xy_Py_b;
  abcd[1241] = 4.0E0*I_NAI_H4xz_Fy2z_bb-2.0E0*1*I_NAI_H4xz_Py_b;
  abcd[1242] = 4.0E0*I_NAI_H3x2y_Fy2z_bb-2.0E0*1*I_NAI_H3x2y_Py_b;
  abcd[1243] = 4.0E0*I_NAI_H3xyz_Fy2z_bb-2.0E0*1*I_NAI_H3xyz_Py_b;
  abcd[1244] = 4.0E0*I_NAI_H3x2z_Fy2z_bb-2.0E0*1*I_NAI_H3x2z_Py_b;
  abcd[1245] = 4.0E0*I_NAI_H2x3y_Fy2z_bb-2.0E0*1*I_NAI_H2x3y_Py_b;
  abcd[1246] = 4.0E0*I_NAI_H2x2yz_Fy2z_bb-2.0E0*1*I_NAI_H2x2yz_Py_b;
  abcd[1247] = 4.0E0*I_NAI_H2xy2z_Fy2z_bb-2.0E0*1*I_NAI_H2xy2z_Py_b;
  abcd[1248] = 4.0E0*I_NAI_H2x3z_Fy2z_bb-2.0E0*1*I_NAI_H2x3z_Py_b;
  abcd[1249] = 4.0E0*I_NAI_Hx4y_Fy2z_bb-2.0E0*1*I_NAI_Hx4y_Py_b;
  abcd[1250] = 4.0E0*I_NAI_Hx3yz_Fy2z_bb-2.0E0*1*I_NAI_Hx3yz_Py_b;
  abcd[1251] = 4.0E0*I_NAI_Hx2y2z_Fy2z_bb-2.0E0*1*I_NAI_Hx2y2z_Py_b;
  abcd[1252] = 4.0E0*I_NAI_Hxy3z_Fy2z_bb-2.0E0*1*I_NAI_Hxy3z_Py_b;
  abcd[1253] = 4.0E0*I_NAI_Hx4z_Fy2z_bb-2.0E0*1*I_NAI_Hx4z_Py_b;
  abcd[1254] = 4.0E0*I_NAI_H5y_Fy2z_bb-2.0E0*1*I_NAI_H5y_Py_b;
  abcd[1255] = 4.0E0*I_NAI_H4yz_Fy2z_bb-2.0E0*1*I_NAI_H4yz_Py_b;
  abcd[1256] = 4.0E0*I_NAI_H3y2z_Fy2z_bb-2.0E0*1*I_NAI_H3y2z_Py_b;
  abcd[1257] = 4.0E0*I_NAI_H2y3z_Fy2z_bb-2.0E0*1*I_NAI_H2y3z_Py_b;
  abcd[1258] = 4.0E0*I_NAI_Hy4z_Fy2z_bb-2.0E0*1*I_NAI_Hy4z_Py_b;
  abcd[1259] = 4.0E0*I_NAI_H5z_Fy2z_bb-2.0E0*1*I_NAI_H5z_Py_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_bb
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  abcd[1260] = 4.0E0*I_NAI_H5x_Fx2z_bb-2.0E0*1*I_NAI_H5x_Px_b;
  abcd[1261] = 4.0E0*I_NAI_H4xy_Fx2z_bb-2.0E0*1*I_NAI_H4xy_Px_b;
  abcd[1262] = 4.0E0*I_NAI_H4xz_Fx2z_bb-2.0E0*1*I_NAI_H4xz_Px_b;
  abcd[1263] = 4.0E0*I_NAI_H3x2y_Fx2z_bb-2.0E0*1*I_NAI_H3x2y_Px_b;
  abcd[1264] = 4.0E0*I_NAI_H3xyz_Fx2z_bb-2.0E0*1*I_NAI_H3xyz_Px_b;
  abcd[1265] = 4.0E0*I_NAI_H3x2z_Fx2z_bb-2.0E0*1*I_NAI_H3x2z_Px_b;
  abcd[1266] = 4.0E0*I_NAI_H2x3y_Fx2z_bb-2.0E0*1*I_NAI_H2x3y_Px_b;
  abcd[1267] = 4.0E0*I_NAI_H2x2yz_Fx2z_bb-2.0E0*1*I_NAI_H2x2yz_Px_b;
  abcd[1268] = 4.0E0*I_NAI_H2xy2z_Fx2z_bb-2.0E0*1*I_NAI_H2xy2z_Px_b;
  abcd[1269] = 4.0E0*I_NAI_H2x3z_Fx2z_bb-2.0E0*1*I_NAI_H2x3z_Px_b;
  abcd[1270] = 4.0E0*I_NAI_Hx4y_Fx2z_bb-2.0E0*1*I_NAI_Hx4y_Px_b;
  abcd[1271] = 4.0E0*I_NAI_Hx3yz_Fx2z_bb-2.0E0*1*I_NAI_Hx3yz_Px_b;
  abcd[1272] = 4.0E0*I_NAI_Hx2y2z_Fx2z_bb-2.0E0*1*I_NAI_Hx2y2z_Px_b;
  abcd[1273] = 4.0E0*I_NAI_Hxy3z_Fx2z_bb-2.0E0*1*I_NAI_Hxy3z_Px_b;
  abcd[1274] = 4.0E0*I_NAI_Hx4z_Fx2z_bb-2.0E0*1*I_NAI_Hx4z_Px_b;
  abcd[1275] = 4.0E0*I_NAI_H5y_Fx2z_bb-2.0E0*1*I_NAI_H5y_Px_b;
  abcd[1276] = 4.0E0*I_NAI_H4yz_Fx2z_bb-2.0E0*1*I_NAI_H4yz_Px_b;
  abcd[1277] = 4.0E0*I_NAI_H3y2z_Fx2z_bb-2.0E0*1*I_NAI_H3y2z_Px_b;
  abcd[1278] = 4.0E0*I_NAI_H2y3z_Fx2z_bb-2.0E0*1*I_NAI_H2y3z_Px_b;
  abcd[1279] = 4.0E0*I_NAI_Hy4z_Fx2z_bb-2.0E0*1*I_NAI_Hy4z_Px_b;
  abcd[1280] = 4.0E0*I_NAI_H5z_Fx2z_bb-2.0E0*1*I_NAI_H5z_Px_b;
  abcd[1281] = 4.0E0*I_NAI_H5x_Fy2z_bb-2.0E0*1*I_NAI_H5x_Py_b;
  abcd[1282] = 4.0E0*I_NAI_H4xy_Fy2z_bb-2.0E0*1*I_NAI_H4xy_Py_b;
  abcd[1283] = 4.0E0*I_NAI_H4xz_Fy2z_bb-2.0E0*1*I_NAI_H4xz_Py_b;
  abcd[1284] = 4.0E0*I_NAI_H3x2y_Fy2z_bb-2.0E0*1*I_NAI_H3x2y_Py_b;
  abcd[1285] = 4.0E0*I_NAI_H3xyz_Fy2z_bb-2.0E0*1*I_NAI_H3xyz_Py_b;
  abcd[1286] = 4.0E0*I_NAI_H3x2z_Fy2z_bb-2.0E0*1*I_NAI_H3x2z_Py_b;
  abcd[1287] = 4.0E0*I_NAI_H2x3y_Fy2z_bb-2.0E0*1*I_NAI_H2x3y_Py_b;
  abcd[1288] = 4.0E0*I_NAI_H2x2yz_Fy2z_bb-2.0E0*1*I_NAI_H2x2yz_Py_b;
  abcd[1289] = 4.0E0*I_NAI_H2xy2z_Fy2z_bb-2.0E0*1*I_NAI_H2xy2z_Py_b;
  abcd[1290] = 4.0E0*I_NAI_H2x3z_Fy2z_bb-2.0E0*1*I_NAI_H2x3z_Py_b;
  abcd[1291] = 4.0E0*I_NAI_Hx4y_Fy2z_bb-2.0E0*1*I_NAI_Hx4y_Py_b;
  abcd[1292] = 4.0E0*I_NAI_Hx3yz_Fy2z_bb-2.0E0*1*I_NAI_Hx3yz_Py_b;
  abcd[1293] = 4.0E0*I_NAI_Hx2y2z_Fy2z_bb-2.0E0*1*I_NAI_Hx2y2z_Py_b;
  abcd[1294] = 4.0E0*I_NAI_Hxy3z_Fy2z_bb-2.0E0*1*I_NAI_Hxy3z_Py_b;
  abcd[1295] = 4.0E0*I_NAI_Hx4z_Fy2z_bb-2.0E0*1*I_NAI_Hx4z_Py_b;
  abcd[1296] = 4.0E0*I_NAI_H5y_Fy2z_bb-2.0E0*1*I_NAI_H5y_Py_b;
  abcd[1297] = 4.0E0*I_NAI_H4yz_Fy2z_bb-2.0E0*1*I_NAI_H4yz_Py_b;
  abcd[1298] = 4.0E0*I_NAI_H3y2z_Fy2z_bb-2.0E0*1*I_NAI_H3y2z_Py_b;
  abcd[1299] = 4.0E0*I_NAI_H2y3z_Fy2z_bb-2.0E0*1*I_NAI_H2y3z_Py_b;
  abcd[1300] = 4.0E0*I_NAI_Hy4z_Fy2z_bb-2.0E0*1*I_NAI_Hy4z_Py_b;
  abcd[1301] = 4.0E0*I_NAI_H5z_Fy2z_bb-2.0E0*1*I_NAI_H5z_Py_b;
  abcd[1302] = 4.0E0*I_NAI_H5x_F3z_bb-2.0E0*1*I_NAI_H5x_Pz_b-2.0E0*2*I_NAI_H5x_Pz_b;
  abcd[1303] = 4.0E0*I_NAI_H4xy_F3z_bb-2.0E0*1*I_NAI_H4xy_Pz_b-2.0E0*2*I_NAI_H4xy_Pz_b;
  abcd[1304] = 4.0E0*I_NAI_H4xz_F3z_bb-2.0E0*1*I_NAI_H4xz_Pz_b-2.0E0*2*I_NAI_H4xz_Pz_b;
  abcd[1305] = 4.0E0*I_NAI_H3x2y_F3z_bb-2.0E0*1*I_NAI_H3x2y_Pz_b-2.0E0*2*I_NAI_H3x2y_Pz_b;
  abcd[1306] = 4.0E0*I_NAI_H3xyz_F3z_bb-2.0E0*1*I_NAI_H3xyz_Pz_b-2.0E0*2*I_NAI_H3xyz_Pz_b;
  abcd[1307] = 4.0E0*I_NAI_H3x2z_F3z_bb-2.0E0*1*I_NAI_H3x2z_Pz_b-2.0E0*2*I_NAI_H3x2z_Pz_b;
  abcd[1308] = 4.0E0*I_NAI_H2x3y_F3z_bb-2.0E0*1*I_NAI_H2x3y_Pz_b-2.0E0*2*I_NAI_H2x3y_Pz_b;
  abcd[1309] = 4.0E0*I_NAI_H2x2yz_F3z_bb-2.0E0*1*I_NAI_H2x2yz_Pz_b-2.0E0*2*I_NAI_H2x2yz_Pz_b;
  abcd[1310] = 4.0E0*I_NAI_H2xy2z_F3z_bb-2.0E0*1*I_NAI_H2xy2z_Pz_b-2.0E0*2*I_NAI_H2xy2z_Pz_b;
  abcd[1311] = 4.0E0*I_NAI_H2x3z_F3z_bb-2.0E0*1*I_NAI_H2x3z_Pz_b-2.0E0*2*I_NAI_H2x3z_Pz_b;
  abcd[1312] = 4.0E0*I_NAI_Hx4y_F3z_bb-2.0E0*1*I_NAI_Hx4y_Pz_b-2.0E0*2*I_NAI_Hx4y_Pz_b;
  abcd[1313] = 4.0E0*I_NAI_Hx3yz_F3z_bb-2.0E0*1*I_NAI_Hx3yz_Pz_b-2.0E0*2*I_NAI_Hx3yz_Pz_b;
  abcd[1314] = 4.0E0*I_NAI_Hx2y2z_F3z_bb-2.0E0*1*I_NAI_Hx2y2z_Pz_b-2.0E0*2*I_NAI_Hx2y2z_Pz_b;
  abcd[1315] = 4.0E0*I_NAI_Hxy3z_F3z_bb-2.0E0*1*I_NAI_Hxy3z_Pz_b-2.0E0*2*I_NAI_Hxy3z_Pz_b;
  abcd[1316] = 4.0E0*I_NAI_Hx4z_F3z_bb-2.0E0*1*I_NAI_Hx4z_Pz_b-2.0E0*2*I_NAI_Hx4z_Pz_b;
  abcd[1317] = 4.0E0*I_NAI_H5y_F3z_bb-2.0E0*1*I_NAI_H5y_Pz_b-2.0E0*2*I_NAI_H5y_Pz_b;
  abcd[1318] = 4.0E0*I_NAI_H4yz_F3z_bb-2.0E0*1*I_NAI_H4yz_Pz_b-2.0E0*2*I_NAI_H4yz_Pz_b;
  abcd[1319] = 4.0E0*I_NAI_H3y2z_F3z_bb-2.0E0*1*I_NAI_H3y2z_Pz_b-2.0E0*2*I_NAI_H3y2z_Pz_b;
  abcd[1320] = 4.0E0*I_NAI_H2y3z_F3z_bb-2.0E0*1*I_NAI_H2y3z_Pz_b-2.0E0*2*I_NAI_H2y3z_Pz_b;
  abcd[1321] = 4.0E0*I_NAI_Hy4z_F3z_bb-2.0E0*1*I_NAI_Hy4z_Pz_b-2.0E0*2*I_NAI_Hy4z_Pz_b;
  abcd[1322] = 4.0E0*I_NAI_H5z_F3z_bb-2.0E0*1*I_NAI_H5z_Pz_b-2.0E0*2*I_NAI_H5z_Pz_b;
}
