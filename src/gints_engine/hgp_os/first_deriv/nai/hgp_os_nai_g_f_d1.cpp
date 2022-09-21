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
// BRA1
// X
// Y
// Z
// ####
// BRA2
// X
// Y
// Z
// ####

void hgp_os_nai_g_f_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
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
  Double I_NAI_L8x_S_a = 0.0E0;
  Double I_NAI_L7xy_S_a = 0.0E0;
  Double I_NAI_L7xz_S_a = 0.0E0;
  Double I_NAI_L6x2y_S_a = 0.0E0;
  Double I_NAI_L6xyz_S_a = 0.0E0;
  Double I_NAI_L6x2z_S_a = 0.0E0;
  Double I_NAI_L5x3y_S_a = 0.0E0;
  Double I_NAI_L5x2yz_S_a = 0.0E0;
  Double I_NAI_L5xy2z_S_a = 0.0E0;
  Double I_NAI_L5x3z_S_a = 0.0E0;
  Double I_NAI_L4x4y_S_a = 0.0E0;
  Double I_NAI_L4x3yz_S_a = 0.0E0;
  Double I_NAI_L4x2y2z_S_a = 0.0E0;
  Double I_NAI_L4xy3z_S_a = 0.0E0;
  Double I_NAI_L4x4z_S_a = 0.0E0;
  Double I_NAI_L3x5y_S_a = 0.0E0;
  Double I_NAI_L3x4yz_S_a = 0.0E0;
  Double I_NAI_L3x3y2z_S_a = 0.0E0;
  Double I_NAI_L3x2y3z_S_a = 0.0E0;
  Double I_NAI_L3xy4z_S_a = 0.0E0;
  Double I_NAI_L3x5z_S_a = 0.0E0;
  Double I_NAI_L2x6y_S_a = 0.0E0;
  Double I_NAI_L2x5yz_S_a = 0.0E0;
  Double I_NAI_L2x4y2z_S_a = 0.0E0;
  Double I_NAI_L2x3y3z_S_a = 0.0E0;
  Double I_NAI_L2x2y4z_S_a = 0.0E0;
  Double I_NAI_L2xy5z_S_a = 0.0E0;
  Double I_NAI_L2x6z_S_a = 0.0E0;
  Double I_NAI_Lx7y_S_a = 0.0E0;
  Double I_NAI_Lx6yz_S_a = 0.0E0;
  Double I_NAI_Lx5y2z_S_a = 0.0E0;
  Double I_NAI_Lx4y3z_S_a = 0.0E0;
  Double I_NAI_Lx3y4z_S_a = 0.0E0;
  Double I_NAI_Lx2y5z_S_a = 0.0E0;
  Double I_NAI_Lxy6z_S_a = 0.0E0;
  Double I_NAI_Lx7z_S_a = 0.0E0;
  Double I_NAI_L8y_S_a = 0.0E0;
  Double I_NAI_L7yz_S_a = 0.0E0;
  Double I_NAI_L6y2z_S_a = 0.0E0;
  Double I_NAI_L5y3z_S_a = 0.0E0;
  Double I_NAI_L4y4z_S_a = 0.0E0;
  Double I_NAI_L3y5z_S_a = 0.0E0;
  Double I_NAI_L2y6z_S_a = 0.0E0;
  Double I_NAI_Ly7z_S_a = 0.0E0;
  Double I_NAI_L8z_S_a = 0.0E0;
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
  Double I_NAI_L8x_S_b = 0.0E0;
  Double I_NAI_L7xy_S_b = 0.0E0;
  Double I_NAI_L7xz_S_b = 0.0E0;
  Double I_NAI_L6x2y_S_b = 0.0E0;
  Double I_NAI_L6xyz_S_b = 0.0E0;
  Double I_NAI_L6x2z_S_b = 0.0E0;
  Double I_NAI_L5x3y_S_b = 0.0E0;
  Double I_NAI_L5x2yz_S_b = 0.0E0;
  Double I_NAI_L5xy2z_S_b = 0.0E0;
  Double I_NAI_L5x3z_S_b = 0.0E0;
  Double I_NAI_L4x4y_S_b = 0.0E0;
  Double I_NAI_L4x3yz_S_b = 0.0E0;
  Double I_NAI_L4x2y2z_S_b = 0.0E0;
  Double I_NAI_L4xy3z_S_b = 0.0E0;
  Double I_NAI_L4x4z_S_b = 0.0E0;
  Double I_NAI_L3x5y_S_b = 0.0E0;
  Double I_NAI_L3x4yz_S_b = 0.0E0;
  Double I_NAI_L3x3y2z_S_b = 0.0E0;
  Double I_NAI_L3x2y3z_S_b = 0.0E0;
  Double I_NAI_L3xy4z_S_b = 0.0E0;
  Double I_NAI_L3x5z_S_b = 0.0E0;
  Double I_NAI_L2x6y_S_b = 0.0E0;
  Double I_NAI_L2x5yz_S_b = 0.0E0;
  Double I_NAI_L2x4y2z_S_b = 0.0E0;
  Double I_NAI_L2x3y3z_S_b = 0.0E0;
  Double I_NAI_L2x2y4z_S_b = 0.0E0;
  Double I_NAI_L2xy5z_S_b = 0.0E0;
  Double I_NAI_L2x6z_S_b = 0.0E0;
  Double I_NAI_Lx7y_S_b = 0.0E0;
  Double I_NAI_Lx6yz_S_b = 0.0E0;
  Double I_NAI_Lx5y2z_S_b = 0.0E0;
  Double I_NAI_Lx4y3z_S_b = 0.0E0;
  Double I_NAI_Lx3y4z_S_b = 0.0E0;
  Double I_NAI_Lx2y5z_S_b = 0.0E0;
  Double I_NAI_Lxy6z_S_b = 0.0E0;
  Double I_NAI_Lx7z_S_b = 0.0E0;
  Double I_NAI_L8y_S_b = 0.0E0;
  Double I_NAI_L7yz_S_b = 0.0E0;
  Double I_NAI_L6y2z_S_b = 0.0E0;
  Double I_NAI_L5y3z_S_b = 0.0E0;
  Double I_NAI_L4y4z_S_b = 0.0E0;
  Double I_NAI_L3y5z_S_b = 0.0E0;
  Double I_NAI_L2y6z_S_b = 0.0E0;
  Double I_NAI_Ly7z_S_b = 0.0E0;
  Double I_NAI_L8z_S_b = 0.0E0;
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
       * shell quartet name: SQ_NAI_L_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_L_S_a_coefs = alpha;
      I_NAI_L8x_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S_a += SQ_NAI_L_S_a_coefs*I_NAI_L8z_S_vrr;

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
       * shell quartet name: SQ_NAI_L_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_L_S_b_coefs = beta;
      I_NAI_L8x_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S_b += SQ_NAI_L_S_b_coefs*I_NAI_L8z_S_vrr;

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
   * totally 20 integrals are omitted 
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
   * shell quartet name: SQ_NAI_H_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 7 integrals are omitted 
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
  Double I_NAI_H4xz_Pz = I_NAI_I4x2z_S+ABZ*I_NAI_H4xz_S;
  Double I_NAI_H3xyz_Pz = I_NAI_I3xy2z_S+ABZ*I_NAI_H3xyz_S;
  Double I_NAI_H3x2z_Pz = I_NAI_I3x3z_S+ABZ*I_NAI_H3x2z_S;
  Double I_NAI_H2x2yz_Pz = I_NAI_I2x2y2z_S+ABZ*I_NAI_H2x2yz_S;
  Double I_NAI_H2xy2z_Pz = I_NAI_I2xy3z_S+ABZ*I_NAI_H2xy2z_S;
  Double I_NAI_H2x3z_Pz = I_NAI_I2x4z_S+ABZ*I_NAI_H2x3z_S;
  Double I_NAI_Hx3yz_Pz = I_NAI_Ix3y2z_S+ABZ*I_NAI_Hx3yz_S;
  Double I_NAI_Hx2y2z_Pz = I_NAI_Ix2y3z_S+ABZ*I_NAI_Hx2y2z_S;
  Double I_NAI_Hxy3z_Pz = I_NAI_Ixy4z_S+ABZ*I_NAI_Hxy3z_S;
  Double I_NAI_Hx4z_Pz = I_NAI_Ix5z_S+ABZ*I_NAI_Hx4z_S;
  Double I_NAI_H4yz_Pz = I_NAI_I4y2z_S+ABZ*I_NAI_H4yz_S;
  Double I_NAI_H3y2z_Pz = I_NAI_I3y3z_S+ABZ*I_NAI_H3y2z_S;
  Double I_NAI_H2y3z_Pz = I_NAI_I2y4z_S+ABZ*I_NAI_H2y3z_S;
  Double I_NAI_Hy4z_Pz = I_NAI_Iy5z_S+ABZ*I_NAI_Hy4z_S;
  Double I_NAI_H5z_Pz = I_NAI_I6z_S+ABZ*I_NAI_H5z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P
   * RHS shell quartet name: SQ_NAI_G_P
   ************************************************************/
  Double I_NAI_G4x_D2x = I_NAI_H5x_Px+ABX*I_NAI_G4x_Px;
  Double I_NAI_G3xy_D2x = I_NAI_H4xy_Px+ABX*I_NAI_G3xy_Px;
  Double I_NAI_G3xz_D2x = I_NAI_H4xz_Px+ABX*I_NAI_G3xz_Px;
  Double I_NAI_G2x2y_D2x = I_NAI_H3x2y_Px+ABX*I_NAI_G2x2y_Px;
  Double I_NAI_G2xyz_D2x = I_NAI_H3xyz_Px+ABX*I_NAI_G2xyz_Px;
  Double I_NAI_G2x2z_D2x = I_NAI_H3x2z_Px+ABX*I_NAI_G2x2z_Px;
  Double I_NAI_Gx3y_D2x = I_NAI_H2x3y_Px+ABX*I_NAI_Gx3y_Px;
  Double I_NAI_Gx2yz_D2x = I_NAI_H2x2yz_Px+ABX*I_NAI_Gx2yz_Px;
  Double I_NAI_Gxy2z_D2x = I_NAI_H2xy2z_Px+ABX*I_NAI_Gxy2z_Px;
  Double I_NAI_Gx3z_D2x = I_NAI_H2x3z_Px+ABX*I_NAI_Gx3z_Px;
  Double I_NAI_G4y_D2x = I_NAI_Hx4y_Px+ABX*I_NAI_G4y_Px;
  Double I_NAI_G3yz_D2x = I_NAI_Hx3yz_Px+ABX*I_NAI_G3yz_Px;
  Double I_NAI_G2y2z_D2x = I_NAI_Hx2y2z_Px+ABX*I_NAI_G2y2z_Px;
  Double I_NAI_Gy3z_D2x = I_NAI_Hxy3z_Px+ABX*I_NAI_Gy3z_Px;
  Double I_NAI_G4z_D2x = I_NAI_Hx4z_Px+ABX*I_NAI_G4z_Px;
  Double I_NAI_G4x_Dxy = I_NAI_H4xy_Px+ABY*I_NAI_G4x_Px;
  Double I_NAI_G3xy_Dxy = I_NAI_H3x2y_Px+ABY*I_NAI_G3xy_Px;
  Double I_NAI_G3xz_Dxy = I_NAI_H3xyz_Px+ABY*I_NAI_G3xz_Px;
  Double I_NAI_G2x2y_Dxy = I_NAI_H2x3y_Px+ABY*I_NAI_G2x2y_Px;
  Double I_NAI_G2xyz_Dxy = I_NAI_H2x2yz_Px+ABY*I_NAI_G2xyz_Px;
  Double I_NAI_G2x2z_Dxy = I_NAI_H2xy2z_Px+ABY*I_NAI_G2x2z_Px;
  Double I_NAI_Gx3y_Dxy = I_NAI_Hx4y_Px+ABY*I_NAI_Gx3y_Px;
  Double I_NAI_Gx2yz_Dxy = I_NAI_Hx3yz_Px+ABY*I_NAI_Gx2yz_Px;
  Double I_NAI_Gxy2z_Dxy = I_NAI_Hx2y2z_Px+ABY*I_NAI_Gxy2z_Px;
  Double I_NAI_Gx3z_Dxy = I_NAI_Hxy3z_Px+ABY*I_NAI_Gx3z_Px;
  Double I_NAI_G4y_Dxy = I_NAI_H5y_Px+ABY*I_NAI_G4y_Px;
  Double I_NAI_G3yz_Dxy = I_NAI_H4yz_Px+ABY*I_NAI_G3yz_Px;
  Double I_NAI_G2y2z_Dxy = I_NAI_H3y2z_Px+ABY*I_NAI_G2y2z_Px;
  Double I_NAI_Gy3z_Dxy = I_NAI_H2y3z_Px+ABY*I_NAI_Gy3z_Px;
  Double I_NAI_G4z_Dxy = I_NAI_Hy4z_Px+ABY*I_NAI_G4z_Px;
  Double I_NAI_G4x_Dxz = I_NAI_H4xz_Px+ABZ*I_NAI_G4x_Px;
  Double I_NAI_G3xy_Dxz = I_NAI_H3xyz_Px+ABZ*I_NAI_G3xy_Px;
  Double I_NAI_G3xz_Dxz = I_NAI_H3x2z_Px+ABZ*I_NAI_G3xz_Px;
  Double I_NAI_G2x2y_Dxz = I_NAI_H2x2yz_Px+ABZ*I_NAI_G2x2y_Px;
  Double I_NAI_G2xyz_Dxz = I_NAI_H2xy2z_Px+ABZ*I_NAI_G2xyz_Px;
  Double I_NAI_G2x2z_Dxz = I_NAI_H2x3z_Px+ABZ*I_NAI_G2x2z_Px;
  Double I_NAI_Gx3y_Dxz = I_NAI_Hx3yz_Px+ABZ*I_NAI_Gx3y_Px;
  Double I_NAI_Gx2yz_Dxz = I_NAI_Hx2y2z_Px+ABZ*I_NAI_Gx2yz_Px;
  Double I_NAI_Gxy2z_Dxz = I_NAI_Hxy3z_Px+ABZ*I_NAI_Gxy2z_Px;
  Double I_NAI_Gx3z_Dxz = I_NAI_Hx4z_Px+ABZ*I_NAI_Gx3z_Px;
  Double I_NAI_G4y_Dxz = I_NAI_H4yz_Px+ABZ*I_NAI_G4y_Px;
  Double I_NAI_G3yz_Dxz = I_NAI_H3y2z_Px+ABZ*I_NAI_G3yz_Px;
  Double I_NAI_G2y2z_Dxz = I_NAI_H2y3z_Px+ABZ*I_NAI_G2y2z_Px;
  Double I_NAI_Gy3z_Dxz = I_NAI_Hy4z_Px+ABZ*I_NAI_Gy3z_Px;
  Double I_NAI_G4z_Dxz = I_NAI_H5z_Px+ABZ*I_NAI_G4z_Px;
  Double I_NAI_G4x_D2y = I_NAI_H4xy_Py+ABY*I_NAI_G4x_Py;
  Double I_NAI_G3xy_D2y = I_NAI_H3x2y_Py+ABY*I_NAI_G3xy_Py;
  Double I_NAI_G3xz_D2y = I_NAI_H3xyz_Py+ABY*I_NAI_G3xz_Py;
  Double I_NAI_G2x2y_D2y = I_NAI_H2x3y_Py+ABY*I_NAI_G2x2y_Py;
  Double I_NAI_G2xyz_D2y = I_NAI_H2x2yz_Py+ABY*I_NAI_G2xyz_Py;
  Double I_NAI_G2x2z_D2y = I_NAI_H2xy2z_Py+ABY*I_NAI_G2x2z_Py;
  Double I_NAI_Gx3y_D2y = I_NAI_Hx4y_Py+ABY*I_NAI_Gx3y_Py;
  Double I_NAI_Gx2yz_D2y = I_NAI_Hx3yz_Py+ABY*I_NAI_Gx2yz_Py;
  Double I_NAI_Gxy2z_D2y = I_NAI_Hx2y2z_Py+ABY*I_NAI_Gxy2z_Py;
  Double I_NAI_Gx3z_D2y = I_NAI_Hxy3z_Py+ABY*I_NAI_Gx3z_Py;
  Double I_NAI_G4y_D2y = I_NAI_H5y_Py+ABY*I_NAI_G4y_Py;
  Double I_NAI_G3yz_D2y = I_NAI_H4yz_Py+ABY*I_NAI_G3yz_Py;
  Double I_NAI_G2y2z_D2y = I_NAI_H3y2z_Py+ABY*I_NAI_G2y2z_Py;
  Double I_NAI_Gy3z_D2y = I_NAI_H2y3z_Py+ABY*I_NAI_Gy3z_Py;
  Double I_NAI_G4z_D2y = I_NAI_Hy4z_Py+ABY*I_NAI_G4z_Py;
  Double I_NAI_G4x_Dyz = I_NAI_H4xz_Py+ABZ*I_NAI_G4x_Py;
  Double I_NAI_G3xy_Dyz = I_NAI_H3xyz_Py+ABZ*I_NAI_G3xy_Py;
  Double I_NAI_G3xz_Dyz = I_NAI_H3x2z_Py+ABZ*I_NAI_G3xz_Py;
  Double I_NAI_G2x2y_Dyz = I_NAI_H2x2yz_Py+ABZ*I_NAI_G2x2y_Py;
  Double I_NAI_G2xyz_Dyz = I_NAI_H2xy2z_Py+ABZ*I_NAI_G2xyz_Py;
  Double I_NAI_G2x2z_Dyz = I_NAI_H2x3z_Py+ABZ*I_NAI_G2x2z_Py;
  Double I_NAI_Gx3y_Dyz = I_NAI_Hx3yz_Py+ABZ*I_NAI_Gx3y_Py;
  Double I_NAI_Gx2yz_Dyz = I_NAI_Hx2y2z_Py+ABZ*I_NAI_Gx2yz_Py;
  Double I_NAI_Gxy2z_Dyz = I_NAI_Hxy3z_Py+ABZ*I_NAI_Gxy2z_Py;
  Double I_NAI_Gx3z_Dyz = I_NAI_Hx4z_Py+ABZ*I_NAI_Gx3z_Py;
  Double I_NAI_G4y_Dyz = I_NAI_H4yz_Py+ABZ*I_NAI_G4y_Py;
  Double I_NAI_G3yz_Dyz = I_NAI_H3y2z_Py+ABZ*I_NAI_G3yz_Py;
  Double I_NAI_G2y2z_Dyz = I_NAI_H2y3z_Py+ABZ*I_NAI_G2y2z_Py;
  Double I_NAI_Gy3z_Dyz = I_NAI_Hy4z_Py+ABZ*I_NAI_Gy3z_Py;
  Double I_NAI_G4z_Dyz = I_NAI_H5z_Py+ABZ*I_NAI_G4z_Py;
  Double I_NAI_G4x_D2z = I_NAI_H4xz_Pz+ABZ*I_NAI_G4x_Pz;
  Double I_NAI_G3xy_D2z = I_NAI_H3xyz_Pz+ABZ*I_NAI_G3xy_Pz;
  Double I_NAI_G3xz_D2z = I_NAI_H3x2z_Pz+ABZ*I_NAI_G3xz_Pz;
  Double I_NAI_G2x2y_D2z = I_NAI_H2x2yz_Pz+ABZ*I_NAI_G2x2y_Pz;
  Double I_NAI_G2xyz_D2z = I_NAI_H2xy2z_Pz+ABZ*I_NAI_G2xyz_Pz;
  Double I_NAI_G2x2z_D2z = I_NAI_H2x3z_Pz+ABZ*I_NAI_G2x2z_Pz;
  Double I_NAI_Gx3y_D2z = I_NAI_Hx3yz_Pz+ABZ*I_NAI_Gx3y_Pz;
  Double I_NAI_Gx2yz_D2z = I_NAI_Hx2y2z_Pz+ABZ*I_NAI_Gx2yz_Pz;
  Double I_NAI_Gxy2z_D2z = I_NAI_Hxy3z_Pz+ABZ*I_NAI_Gxy2z_Pz;
  Double I_NAI_Gx3z_D2z = I_NAI_Hx4z_Pz+ABZ*I_NAI_Gx3z_Pz;
  Double I_NAI_G4y_D2z = I_NAI_H4yz_Pz+ABZ*I_NAI_G4y_Pz;
  Double I_NAI_G3yz_D2z = I_NAI_H3y2z_Pz+ABZ*I_NAI_G3yz_Pz;
  Double I_NAI_G2y2z_D2z = I_NAI_H2y3z_Pz+ABZ*I_NAI_G2y2z_Pz;
  Double I_NAI_Gy3z_D2z = I_NAI_Hy4z_Pz+ABZ*I_NAI_Gy3z_Pz;
  Double I_NAI_G4z_D2z = I_NAI_H5z_Pz+ABZ*I_NAI_G4z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  Double I_NAI_F3x_F3x = I_NAI_G4x_D2x+ABX*I_NAI_F3x_D2x;
  Double I_NAI_F2xy_F3x = I_NAI_G3xy_D2x+ABX*I_NAI_F2xy_D2x;
  Double I_NAI_F2xz_F3x = I_NAI_G3xz_D2x+ABX*I_NAI_F2xz_D2x;
  Double I_NAI_Fx2y_F3x = I_NAI_G2x2y_D2x+ABX*I_NAI_Fx2y_D2x;
  Double I_NAI_Fxyz_F3x = I_NAI_G2xyz_D2x+ABX*I_NAI_Fxyz_D2x;
  Double I_NAI_Fx2z_F3x = I_NAI_G2x2z_D2x+ABX*I_NAI_Fx2z_D2x;
  Double I_NAI_F3y_F3x = I_NAI_Gx3y_D2x+ABX*I_NAI_F3y_D2x;
  Double I_NAI_F2yz_F3x = I_NAI_Gx2yz_D2x+ABX*I_NAI_F2yz_D2x;
  Double I_NAI_Fy2z_F3x = I_NAI_Gxy2z_D2x+ABX*I_NAI_Fy2z_D2x;
  Double I_NAI_F3z_F3x = I_NAI_Gx3z_D2x+ABX*I_NAI_F3z_D2x;
  Double I_NAI_F3x_F2xy = I_NAI_G3xy_D2x+ABY*I_NAI_F3x_D2x;
  Double I_NAI_F2xy_F2xy = I_NAI_G2x2y_D2x+ABY*I_NAI_F2xy_D2x;
  Double I_NAI_F2xz_F2xy = I_NAI_G2xyz_D2x+ABY*I_NAI_F2xz_D2x;
  Double I_NAI_Fx2y_F2xy = I_NAI_Gx3y_D2x+ABY*I_NAI_Fx2y_D2x;
  Double I_NAI_Fxyz_F2xy = I_NAI_Gx2yz_D2x+ABY*I_NAI_Fxyz_D2x;
  Double I_NAI_Fx2z_F2xy = I_NAI_Gxy2z_D2x+ABY*I_NAI_Fx2z_D2x;
  Double I_NAI_F3y_F2xy = I_NAI_G4y_D2x+ABY*I_NAI_F3y_D2x;
  Double I_NAI_F2yz_F2xy = I_NAI_G3yz_D2x+ABY*I_NAI_F2yz_D2x;
  Double I_NAI_Fy2z_F2xy = I_NAI_G2y2z_D2x+ABY*I_NAI_Fy2z_D2x;
  Double I_NAI_F3z_F2xy = I_NAI_Gy3z_D2x+ABY*I_NAI_F3z_D2x;
  Double I_NAI_F3x_F2xz = I_NAI_G3xz_D2x+ABZ*I_NAI_F3x_D2x;
  Double I_NAI_F2xy_F2xz = I_NAI_G2xyz_D2x+ABZ*I_NAI_F2xy_D2x;
  Double I_NAI_F2xz_F2xz = I_NAI_G2x2z_D2x+ABZ*I_NAI_F2xz_D2x;
  Double I_NAI_Fx2y_F2xz = I_NAI_Gx2yz_D2x+ABZ*I_NAI_Fx2y_D2x;
  Double I_NAI_Fxyz_F2xz = I_NAI_Gxy2z_D2x+ABZ*I_NAI_Fxyz_D2x;
  Double I_NAI_Fx2z_F2xz = I_NAI_Gx3z_D2x+ABZ*I_NAI_Fx2z_D2x;
  Double I_NAI_F3y_F2xz = I_NAI_G3yz_D2x+ABZ*I_NAI_F3y_D2x;
  Double I_NAI_F2yz_F2xz = I_NAI_G2y2z_D2x+ABZ*I_NAI_F2yz_D2x;
  Double I_NAI_Fy2z_F2xz = I_NAI_Gy3z_D2x+ABZ*I_NAI_Fy2z_D2x;
  Double I_NAI_F3z_F2xz = I_NAI_G4z_D2x+ABZ*I_NAI_F3z_D2x;
  Double I_NAI_F3x_Fx2y = I_NAI_G4x_D2y+ABX*I_NAI_F3x_D2y;
  Double I_NAI_F2xy_Fx2y = I_NAI_G3xy_D2y+ABX*I_NAI_F2xy_D2y;
  Double I_NAI_F2xz_Fx2y = I_NAI_G3xz_D2y+ABX*I_NAI_F2xz_D2y;
  Double I_NAI_Fx2y_Fx2y = I_NAI_G2x2y_D2y+ABX*I_NAI_Fx2y_D2y;
  Double I_NAI_Fxyz_Fx2y = I_NAI_G2xyz_D2y+ABX*I_NAI_Fxyz_D2y;
  Double I_NAI_Fx2z_Fx2y = I_NAI_G2x2z_D2y+ABX*I_NAI_Fx2z_D2y;
  Double I_NAI_F3y_Fx2y = I_NAI_Gx3y_D2y+ABX*I_NAI_F3y_D2y;
  Double I_NAI_F2yz_Fx2y = I_NAI_Gx2yz_D2y+ABX*I_NAI_F2yz_D2y;
  Double I_NAI_Fy2z_Fx2y = I_NAI_Gxy2z_D2y+ABX*I_NAI_Fy2z_D2y;
  Double I_NAI_F3z_Fx2y = I_NAI_Gx3z_D2y+ABX*I_NAI_F3z_D2y;
  Double I_NAI_F3x_Fxyz = I_NAI_G3xz_Dxy+ABZ*I_NAI_F3x_Dxy;
  Double I_NAI_F2xy_Fxyz = I_NAI_G2xyz_Dxy+ABZ*I_NAI_F2xy_Dxy;
  Double I_NAI_F2xz_Fxyz = I_NAI_G2x2z_Dxy+ABZ*I_NAI_F2xz_Dxy;
  Double I_NAI_Fx2y_Fxyz = I_NAI_Gx2yz_Dxy+ABZ*I_NAI_Fx2y_Dxy;
  Double I_NAI_Fxyz_Fxyz = I_NAI_Gxy2z_Dxy+ABZ*I_NAI_Fxyz_Dxy;
  Double I_NAI_Fx2z_Fxyz = I_NAI_Gx3z_Dxy+ABZ*I_NAI_Fx2z_Dxy;
  Double I_NAI_F3y_Fxyz = I_NAI_G3yz_Dxy+ABZ*I_NAI_F3y_Dxy;
  Double I_NAI_F2yz_Fxyz = I_NAI_G2y2z_Dxy+ABZ*I_NAI_F2yz_Dxy;
  Double I_NAI_Fy2z_Fxyz = I_NAI_Gy3z_Dxy+ABZ*I_NAI_Fy2z_Dxy;
  Double I_NAI_F3z_Fxyz = I_NAI_G4z_Dxy+ABZ*I_NAI_F3z_Dxy;
  Double I_NAI_F3x_Fx2z = I_NAI_G4x_D2z+ABX*I_NAI_F3x_D2z;
  Double I_NAI_F2xy_Fx2z = I_NAI_G3xy_D2z+ABX*I_NAI_F2xy_D2z;
  Double I_NAI_F2xz_Fx2z = I_NAI_G3xz_D2z+ABX*I_NAI_F2xz_D2z;
  Double I_NAI_Fx2y_Fx2z = I_NAI_G2x2y_D2z+ABX*I_NAI_Fx2y_D2z;
  Double I_NAI_Fxyz_Fx2z = I_NAI_G2xyz_D2z+ABX*I_NAI_Fxyz_D2z;
  Double I_NAI_Fx2z_Fx2z = I_NAI_G2x2z_D2z+ABX*I_NAI_Fx2z_D2z;
  Double I_NAI_F3y_Fx2z = I_NAI_Gx3y_D2z+ABX*I_NAI_F3y_D2z;
  Double I_NAI_F2yz_Fx2z = I_NAI_Gx2yz_D2z+ABX*I_NAI_F2yz_D2z;
  Double I_NAI_Fy2z_Fx2z = I_NAI_Gxy2z_D2z+ABX*I_NAI_Fy2z_D2z;
  Double I_NAI_F3z_Fx2z = I_NAI_Gx3z_D2z+ABX*I_NAI_F3z_D2z;
  Double I_NAI_F3x_F3y = I_NAI_G3xy_D2y+ABY*I_NAI_F3x_D2y;
  Double I_NAI_F2xy_F3y = I_NAI_G2x2y_D2y+ABY*I_NAI_F2xy_D2y;
  Double I_NAI_F2xz_F3y = I_NAI_G2xyz_D2y+ABY*I_NAI_F2xz_D2y;
  Double I_NAI_Fx2y_F3y = I_NAI_Gx3y_D2y+ABY*I_NAI_Fx2y_D2y;
  Double I_NAI_Fxyz_F3y = I_NAI_Gx2yz_D2y+ABY*I_NAI_Fxyz_D2y;
  Double I_NAI_Fx2z_F3y = I_NAI_Gxy2z_D2y+ABY*I_NAI_Fx2z_D2y;
  Double I_NAI_F3y_F3y = I_NAI_G4y_D2y+ABY*I_NAI_F3y_D2y;
  Double I_NAI_F2yz_F3y = I_NAI_G3yz_D2y+ABY*I_NAI_F2yz_D2y;
  Double I_NAI_Fy2z_F3y = I_NAI_G2y2z_D2y+ABY*I_NAI_Fy2z_D2y;
  Double I_NAI_F3z_F3y = I_NAI_Gy3z_D2y+ABY*I_NAI_F3z_D2y;
  Double I_NAI_F3x_F2yz = I_NAI_G3xz_D2y+ABZ*I_NAI_F3x_D2y;
  Double I_NAI_F2xy_F2yz = I_NAI_G2xyz_D2y+ABZ*I_NAI_F2xy_D2y;
  Double I_NAI_F2xz_F2yz = I_NAI_G2x2z_D2y+ABZ*I_NAI_F2xz_D2y;
  Double I_NAI_Fx2y_F2yz = I_NAI_Gx2yz_D2y+ABZ*I_NAI_Fx2y_D2y;
  Double I_NAI_Fxyz_F2yz = I_NAI_Gxy2z_D2y+ABZ*I_NAI_Fxyz_D2y;
  Double I_NAI_Fx2z_F2yz = I_NAI_Gx3z_D2y+ABZ*I_NAI_Fx2z_D2y;
  Double I_NAI_F3y_F2yz = I_NAI_G3yz_D2y+ABZ*I_NAI_F3y_D2y;
  Double I_NAI_F2yz_F2yz = I_NAI_G2y2z_D2y+ABZ*I_NAI_F2yz_D2y;
  Double I_NAI_Fy2z_F2yz = I_NAI_Gy3z_D2y+ABZ*I_NAI_Fy2z_D2y;
  Double I_NAI_F3z_F2yz = I_NAI_G4z_D2y+ABZ*I_NAI_F3z_D2y;
  Double I_NAI_F3x_Fy2z = I_NAI_G3xy_D2z+ABY*I_NAI_F3x_D2z;
  Double I_NAI_F2xy_Fy2z = I_NAI_G2x2y_D2z+ABY*I_NAI_F2xy_D2z;
  Double I_NAI_F2xz_Fy2z = I_NAI_G2xyz_D2z+ABY*I_NAI_F2xz_D2z;
  Double I_NAI_Fx2y_Fy2z = I_NAI_Gx3y_D2z+ABY*I_NAI_Fx2y_D2z;
  Double I_NAI_Fxyz_Fy2z = I_NAI_Gx2yz_D2z+ABY*I_NAI_Fxyz_D2z;
  Double I_NAI_Fx2z_Fy2z = I_NAI_Gxy2z_D2z+ABY*I_NAI_Fx2z_D2z;
  Double I_NAI_F3y_Fy2z = I_NAI_G4y_D2z+ABY*I_NAI_F3y_D2z;
  Double I_NAI_F2yz_Fy2z = I_NAI_G3yz_D2z+ABY*I_NAI_F2yz_D2z;
  Double I_NAI_Fy2z_Fy2z = I_NAI_G2y2z_D2z+ABY*I_NAI_Fy2z_D2z;
  Double I_NAI_F3z_Fy2z = I_NAI_Gy3z_D2z+ABY*I_NAI_F3z_D2z;
  Double I_NAI_F3x_F3z = I_NAI_G3xz_D2z+ABZ*I_NAI_F3x_D2z;
  Double I_NAI_F2xy_F3z = I_NAI_G2xyz_D2z+ABZ*I_NAI_F2xy_D2z;
  Double I_NAI_F2xz_F3z = I_NAI_G2x2z_D2z+ABZ*I_NAI_F2xz_D2z;
  Double I_NAI_Fx2y_F3z = I_NAI_Gx2yz_D2z+ABZ*I_NAI_Fx2y_D2z;
  Double I_NAI_Fxyz_F3z = I_NAI_Gxy2z_D2z+ABZ*I_NAI_Fxyz_D2z;
  Double I_NAI_Fx2z_F3z = I_NAI_Gx3z_D2z+ABZ*I_NAI_Fx2z_D2z;
  Double I_NAI_F3y_F3z = I_NAI_G3yz_D2z+ABZ*I_NAI_F3y_D2z;
  Double I_NAI_F2yz_F3z = I_NAI_G2y2z_D2z+ABZ*I_NAI_F2yz_D2z;
  Double I_NAI_Fy2z_F3z = I_NAI_Gy3z_D2z+ABZ*I_NAI_Fy2z_D2z;
  Double I_NAI_F3z_F3z = I_NAI_G4z_D2z+ABZ*I_NAI_F3z_D2z;

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
   * shell quartet name: SQ_NAI_I_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_NAI_I6x_Py_a = I_NAI_K6xy_S_a+ABY*I_NAI_I6x_S_a;
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
  Double I_NAI_I6x_Pz_a = I_NAI_K6xz_S_a+ABZ*I_NAI_I6x_S_a;
  Double I_NAI_I5xy_Pz_a = I_NAI_K5xyz_S_a+ABZ*I_NAI_I5xy_S_a;
  Double I_NAI_I5xz_Pz_a = I_NAI_K5x2z_S_a+ABZ*I_NAI_I5xz_S_a;
  Double I_NAI_I4x2y_Pz_a = I_NAI_K4x2yz_S_a+ABZ*I_NAI_I4x2y_S_a;
  Double I_NAI_I4xyz_Pz_a = I_NAI_K4xy2z_S_a+ABZ*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Pz_a = I_NAI_K4x3z_S_a+ABZ*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x3y_Pz_a = I_NAI_K3x3yz_S_a+ABZ*I_NAI_I3x3y_S_a;
  Double I_NAI_I3x2yz_Pz_a = I_NAI_K3x2y2z_S_a+ABZ*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Pz_a = I_NAI_K3xy3z_S_a+ABZ*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Pz_a = I_NAI_K3x4z_S_a+ABZ*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x4y_Pz_a = I_NAI_K2x4yz_S_a+ABZ*I_NAI_I2x4y_S_a;
  Double I_NAI_I2x3yz_Pz_a = I_NAI_K2x3y2z_S_a+ABZ*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Pz_a = I_NAI_K2x2y3z_S_a+ABZ*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Pz_a = I_NAI_K2xy4z_S_a+ABZ*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Pz_a = I_NAI_K2x5z_S_a+ABZ*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix5y_Pz_a = I_NAI_Kx5yz_S_a+ABZ*I_NAI_Ix5y_S_a;
  Double I_NAI_Ix4yz_Pz_a = I_NAI_Kx4y2z_S_a+ABZ*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Pz_a = I_NAI_Kx3y3z_S_a+ABZ*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Pz_a = I_NAI_Kx2y4z_S_a+ABZ*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Pz_a = I_NAI_Kxy5z_S_a+ABZ*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Pz_a = I_NAI_Kx6z_S_a+ABZ*I_NAI_Ix5z_S_a;
  Double I_NAI_I6y_Pz_a = I_NAI_K6yz_S_a+ABZ*I_NAI_I6y_S_a;
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
   * totally 42 integrals are omitted 
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
   * shell quartet name: SQ_NAI_K_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S_a
   * RHS shell quartet name: SQ_NAI_K_S_a
   ************************************************************/
  Double I_NAI_K7x_Px_a = I_NAI_L8x_S_a+ABX*I_NAI_K7x_S_a;
  Double I_NAI_K6xy_Px_a = I_NAI_L7xy_S_a+ABX*I_NAI_K6xy_S_a;
  Double I_NAI_K6xz_Px_a = I_NAI_L7xz_S_a+ABX*I_NAI_K6xz_S_a;
  Double I_NAI_K5x2y_Px_a = I_NAI_L6x2y_S_a+ABX*I_NAI_K5x2y_S_a;
  Double I_NAI_K5xyz_Px_a = I_NAI_L6xyz_S_a+ABX*I_NAI_K5xyz_S_a;
  Double I_NAI_K5x2z_Px_a = I_NAI_L6x2z_S_a+ABX*I_NAI_K5x2z_S_a;
  Double I_NAI_K4x3y_Px_a = I_NAI_L5x3y_S_a+ABX*I_NAI_K4x3y_S_a;
  Double I_NAI_K4x2yz_Px_a = I_NAI_L5x2yz_S_a+ABX*I_NAI_K4x2yz_S_a;
  Double I_NAI_K4xy2z_Px_a = I_NAI_L5xy2z_S_a+ABX*I_NAI_K4xy2z_S_a;
  Double I_NAI_K4x3z_Px_a = I_NAI_L5x3z_S_a+ABX*I_NAI_K4x3z_S_a;
  Double I_NAI_K3x4y_Px_a = I_NAI_L4x4y_S_a+ABX*I_NAI_K3x4y_S_a;
  Double I_NAI_K3x3yz_Px_a = I_NAI_L4x3yz_S_a+ABX*I_NAI_K3x3yz_S_a;
  Double I_NAI_K3x2y2z_Px_a = I_NAI_L4x2y2z_S_a+ABX*I_NAI_K3x2y2z_S_a;
  Double I_NAI_K3xy3z_Px_a = I_NAI_L4xy3z_S_a+ABX*I_NAI_K3xy3z_S_a;
  Double I_NAI_K3x4z_Px_a = I_NAI_L4x4z_S_a+ABX*I_NAI_K3x4z_S_a;
  Double I_NAI_K2x5y_Px_a = I_NAI_L3x5y_S_a+ABX*I_NAI_K2x5y_S_a;
  Double I_NAI_K2x4yz_Px_a = I_NAI_L3x4yz_S_a+ABX*I_NAI_K2x4yz_S_a;
  Double I_NAI_K2x3y2z_Px_a = I_NAI_L3x3y2z_S_a+ABX*I_NAI_K2x3y2z_S_a;
  Double I_NAI_K2x2y3z_Px_a = I_NAI_L3x2y3z_S_a+ABX*I_NAI_K2x2y3z_S_a;
  Double I_NAI_K2xy4z_Px_a = I_NAI_L3xy4z_S_a+ABX*I_NAI_K2xy4z_S_a;
  Double I_NAI_K2x5z_Px_a = I_NAI_L3x5z_S_a+ABX*I_NAI_K2x5z_S_a;
  Double I_NAI_Kx6y_Px_a = I_NAI_L2x6y_S_a+ABX*I_NAI_Kx6y_S_a;
  Double I_NAI_Kx5yz_Px_a = I_NAI_L2x5yz_S_a+ABX*I_NAI_Kx5yz_S_a;
  Double I_NAI_Kx4y2z_Px_a = I_NAI_L2x4y2z_S_a+ABX*I_NAI_Kx4y2z_S_a;
  Double I_NAI_Kx3y3z_Px_a = I_NAI_L2x3y3z_S_a+ABX*I_NAI_Kx3y3z_S_a;
  Double I_NAI_Kx2y4z_Px_a = I_NAI_L2x2y4z_S_a+ABX*I_NAI_Kx2y4z_S_a;
  Double I_NAI_Kxy5z_Px_a = I_NAI_L2xy5z_S_a+ABX*I_NAI_Kxy5z_S_a;
  Double I_NAI_Kx6z_Px_a = I_NAI_L2x6z_S_a+ABX*I_NAI_Kx6z_S_a;
  Double I_NAI_K6yz_Px_a = I_NAI_Lx6yz_S_a+ABX*I_NAI_K6yz_S_a;
  Double I_NAI_K5y2z_Px_a = I_NAI_Lx5y2z_S_a+ABX*I_NAI_K5y2z_S_a;
  Double I_NAI_K4y3z_Px_a = I_NAI_Lx4y3z_S_a+ABX*I_NAI_K4y3z_S_a;
  Double I_NAI_K3y4z_Px_a = I_NAI_Lx3y4z_S_a+ABX*I_NAI_K3y4z_S_a;
  Double I_NAI_K2y5z_Px_a = I_NAI_Lx2y5z_S_a+ABX*I_NAI_K2y5z_S_a;
  Double I_NAI_Ky6z_Px_a = I_NAI_Lxy6z_S_a+ABX*I_NAI_Ky6z_S_a;
  Double I_NAI_K6xy_Py_a = I_NAI_L6x2y_S_a+ABY*I_NAI_K6xy_S_a;
  Double I_NAI_K5x2y_Py_a = I_NAI_L5x3y_S_a+ABY*I_NAI_K5x2y_S_a;
  Double I_NAI_K5xyz_Py_a = I_NAI_L5x2yz_S_a+ABY*I_NAI_K5xyz_S_a;
  Double I_NAI_K4x3y_Py_a = I_NAI_L4x4y_S_a+ABY*I_NAI_K4x3y_S_a;
  Double I_NAI_K4x2yz_Py_a = I_NAI_L4x3yz_S_a+ABY*I_NAI_K4x2yz_S_a;
  Double I_NAI_K4xy2z_Py_a = I_NAI_L4x2y2z_S_a+ABY*I_NAI_K4xy2z_S_a;
  Double I_NAI_K3x4y_Py_a = I_NAI_L3x5y_S_a+ABY*I_NAI_K3x4y_S_a;
  Double I_NAI_K3x3yz_Py_a = I_NAI_L3x4yz_S_a+ABY*I_NAI_K3x3yz_S_a;
  Double I_NAI_K3x2y2z_Py_a = I_NAI_L3x3y2z_S_a+ABY*I_NAI_K3x2y2z_S_a;
  Double I_NAI_K3xy3z_Py_a = I_NAI_L3x2y3z_S_a+ABY*I_NAI_K3xy3z_S_a;
  Double I_NAI_K2x5y_Py_a = I_NAI_L2x6y_S_a+ABY*I_NAI_K2x5y_S_a;
  Double I_NAI_K2x4yz_Py_a = I_NAI_L2x5yz_S_a+ABY*I_NAI_K2x4yz_S_a;
  Double I_NAI_K2x3y2z_Py_a = I_NAI_L2x4y2z_S_a+ABY*I_NAI_K2x3y2z_S_a;
  Double I_NAI_K2x2y3z_Py_a = I_NAI_L2x3y3z_S_a+ABY*I_NAI_K2x2y3z_S_a;
  Double I_NAI_K2xy4z_Py_a = I_NAI_L2x2y4z_S_a+ABY*I_NAI_K2xy4z_S_a;
  Double I_NAI_Kx6y_Py_a = I_NAI_Lx7y_S_a+ABY*I_NAI_Kx6y_S_a;
  Double I_NAI_Kx5yz_Py_a = I_NAI_Lx6yz_S_a+ABY*I_NAI_Kx5yz_S_a;
  Double I_NAI_Kx4y2z_Py_a = I_NAI_Lx5y2z_S_a+ABY*I_NAI_Kx4y2z_S_a;
  Double I_NAI_Kx3y3z_Py_a = I_NAI_Lx4y3z_S_a+ABY*I_NAI_Kx3y3z_S_a;
  Double I_NAI_Kx2y4z_Py_a = I_NAI_Lx3y4z_S_a+ABY*I_NAI_Kx2y4z_S_a;
  Double I_NAI_Kxy5z_Py_a = I_NAI_Lx2y5z_S_a+ABY*I_NAI_Kxy5z_S_a;
  Double I_NAI_K7y_Py_a = I_NAI_L8y_S_a+ABY*I_NAI_K7y_S_a;
  Double I_NAI_K6yz_Py_a = I_NAI_L7yz_S_a+ABY*I_NAI_K6yz_S_a;
  Double I_NAI_K5y2z_Py_a = I_NAI_L6y2z_S_a+ABY*I_NAI_K5y2z_S_a;
  Double I_NAI_K4y3z_Py_a = I_NAI_L5y3z_S_a+ABY*I_NAI_K4y3z_S_a;
  Double I_NAI_K3y4z_Py_a = I_NAI_L4y4z_S_a+ABY*I_NAI_K3y4z_S_a;
  Double I_NAI_K2y5z_Py_a = I_NAI_L3y5z_S_a+ABY*I_NAI_K2y5z_S_a;
  Double I_NAI_Ky6z_Py_a = I_NAI_L2y6z_S_a+ABY*I_NAI_Ky6z_S_a;
  Double I_NAI_K6xz_Pz_a = I_NAI_L6x2z_S_a+ABZ*I_NAI_K6xz_S_a;
  Double I_NAI_K5xyz_Pz_a = I_NAI_L5xy2z_S_a+ABZ*I_NAI_K5xyz_S_a;
  Double I_NAI_K5x2z_Pz_a = I_NAI_L5x3z_S_a+ABZ*I_NAI_K5x2z_S_a;
  Double I_NAI_K4x2yz_Pz_a = I_NAI_L4x2y2z_S_a+ABZ*I_NAI_K4x2yz_S_a;
  Double I_NAI_K4xy2z_Pz_a = I_NAI_L4xy3z_S_a+ABZ*I_NAI_K4xy2z_S_a;
  Double I_NAI_K4x3z_Pz_a = I_NAI_L4x4z_S_a+ABZ*I_NAI_K4x3z_S_a;
  Double I_NAI_K3x3yz_Pz_a = I_NAI_L3x3y2z_S_a+ABZ*I_NAI_K3x3yz_S_a;
  Double I_NAI_K3x2y2z_Pz_a = I_NAI_L3x2y3z_S_a+ABZ*I_NAI_K3x2y2z_S_a;
  Double I_NAI_K3xy3z_Pz_a = I_NAI_L3xy4z_S_a+ABZ*I_NAI_K3xy3z_S_a;
  Double I_NAI_K3x4z_Pz_a = I_NAI_L3x5z_S_a+ABZ*I_NAI_K3x4z_S_a;
  Double I_NAI_K2x4yz_Pz_a = I_NAI_L2x4y2z_S_a+ABZ*I_NAI_K2x4yz_S_a;
  Double I_NAI_K2x3y2z_Pz_a = I_NAI_L2x3y3z_S_a+ABZ*I_NAI_K2x3y2z_S_a;
  Double I_NAI_K2x2y3z_Pz_a = I_NAI_L2x2y4z_S_a+ABZ*I_NAI_K2x2y3z_S_a;
  Double I_NAI_K2xy4z_Pz_a = I_NAI_L2xy5z_S_a+ABZ*I_NAI_K2xy4z_S_a;
  Double I_NAI_K2x5z_Pz_a = I_NAI_L2x6z_S_a+ABZ*I_NAI_K2x5z_S_a;
  Double I_NAI_Kx5yz_Pz_a = I_NAI_Lx5y2z_S_a+ABZ*I_NAI_Kx5yz_S_a;
  Double I_NAI_Kx4y2z_Pz_a = I_NAI_Lx4y3z_S_a+ABZ*I_NAI_Kx4y2z_S_a;
  Double I_NAI_Kx3y3z_Pz_a = I_NAI_Lx3y4z_S_a+ABZ*I_NAI_Kx3y3z_S_a;
  Double I_NAI_Kx2y4z_Pz_a = I_NAI_Lx2y5z_S_a+ABZ*I_NAI_Kx2y4z_S_a;
  Double I_NAI_Kxy5z_Pz_a = I_NAI_Lxy6z_S_a+ABZ*I_NAI_Kxy5z_S_a;
  Double I_NAI_Kx6z_Pz_a = I_NAI_Lx7z_S_a+ABZ*I_NAI_Kx6z_S_a;
  Double I_NAI_K6yz_Pz_a = I_NAI_L6y2z_S_a+ABZ*I_NAI_K6yz_S_a;
  Double I_NAI_K5y2z_Pz_a = I_NAI_L5y3z_S_a+ABZ*I_NAI_K5y2z_S_a;
  Double I_NAI_K4y3z_Pz_a = I_NAI_L4y4z_S_a+ABZ*I_NAI_K4y3z_S_a;
  Double I_NAI_K3y4z_Pz_a = I_NAI_L3y5z_S_a+ABZ*I_NAI_K3y4z_S_a;
  Double I_NAI_K2y5z_Pz_a = I_NAI_L2y6z_S_a+ABZ*I_NAI_K2y5z_S_a;
  Double I_NAI_Ky6z_Pz_a = I_NAI_Ly7z_S_a+ABZ*I_NAI_Ky6z_S_a;
  Double I_NAI_K7z_Pz_a = I_NAI_L8z_S_a+ABZ*I_NAI_K7z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_a
   * RHS shell quartet name: SQ_NAI_I_P_a
   ************************************************************/
  Double I_NAI_I6x_D2x_a = I_NAI_K7x_Px_a+ABX*I_NAI_I6x_Px_a;
  Double I_NAI_I5xy_D2x_a = I_NAI_K6xy_Px_a+ABX*I_NAI_I5xy_Px_a;
  Double I_NAI_I5xz_D2x_a = I_NAI_K6xz_Px_a+ABX*I_NAI_I5xz_Px_a;
  Double I_NAI_I4x2y_D2x_a = I_NAI_K5x2y_Px_a+ABX*I_NAI_I4x2y_Px_a;
  Double I_NAI_I4xyz_D2x_a = I_NAI_K5xyz_Px_a+ABX*I_NAI_I4xyz_Px_a;
  Double I_NAI_I4x2z_D2x_a = I_NAI_K5x2z_Px_a+ABX*I_NAI_I4x2z_Px_a;
  Double I_NAI_I3x3y_D2x_a = I_NAI_K4x3y_Px_a+ABX*I_NAI_I3x3y_Px_a;
  Double I_NAI_I3x2yz_D2x_a = I_NAI_K4x2yz_Px_a+ABX*I_NAI_I3x2yz_Px_a;
  Double I_NAI_I3xy2z_D2x_a = I_NAI_K4xy2z_Px_a+ABX*I_NAI_I3xy2z_Px_a;
  Double I_NAI_I3x3z_D2x_a = I_NAI_K4x3z_Px_a+ABX*I_NAI_I3x3z_Px_a;
  Double I_NAI_I2x4y_D2x_a = I_NAI_K3x4y_Px_a+ABX*I_NAI_I2x4y_Px_a;
  Double I_NAI_I2x3yz_D2x_a = I_NAI_K3x3yz_Px_a+ABX*I_NAI_I2x3yz_Px_a;
  Double I_NAI_I2x2y2z_D2x_a = I_NAI_K3x2y2z_Px_a+ABX*I_NAI_I2x2y2z_Px_a;
  Double I_NAI_I2xy3z_D2x_a = I_NAI_K3xy3z_Px_a+ABX*I_NAI_I2xy3z_Px_a;
  Double I_NAI_I2x4z_D2x_a = I_NAI_K3x4z_Px_a+ABX*I_NAI_I2x4z_Px_a;
  Double I_NAI_Ix5y_D2x_a = I_NAI_K2x5y_Px_a+ABX*I_NAI_Ix5y_Px_a;
  Double I_NAI_Ix4yz_D2x_a = I_NAI_K2x4yz_Px_a+ABX*I_NAI_Ix4yz_Px_a;
  Double I_NAI_Ix3y2z_D2x_a = I_NAI_K2x3y2z_Px_a+ABX*I_NAI_Ix3y2z_Px_a;
  Double I_NAI_Ix2y3z_D2x_a = I_NAI_K2x2y3z_Px_a+ABX*I_NAI_Ix2y3z_Px_a;
  Double I_NAI_Ixy4z_D2x_a = I_NAI_K2xy4z_Px_a+ABX*I_NAI_Ixy4z_Px_a;
  Double I_NAI_Ix5z_D2x_a = I_NAI_K2x5z_Px_a+ABX*I_NAI_Ix5z_Px_a;
  Double I_NAI_I6y_D2x_a = I_NAI_Kx6y_Px_a+ABX*I_NAI_I6y_Px_a;
  Double I_NAI_I5yz_D2x_a = I_NAI_Kx5yz_Px_a+ABX*I_NAI_I5yz_Px_a;
  Double I_NAI_I4y2z_D2x_a = I_NAI_Kx4y2z_Px_a+ABX*I_NAI_I4y2z_Px_a;
  Double I_NAI_I3y3z_D2x_a = I_NAI_Kx3y3z_Px_a+ABX*I_NAI_I3y3z_Px_a;
  Double I_NAI_I2y4z_D2x_a = I_NAI_Kx2y4z_Px_a+ABX*I_NAI_I2y4z_Px_a;
  Double I_NAI_Iy5z_D2x_a = I_NAI_Kxy5z_Px_a+ABX*I_NAI_Iy5z_Px_a;
  Double I_NAI_I6z_D2x_a = I_NAI_Kx6z_Px_a+ABX*I_NAI_I6z_Px_a;
  Double I_NAI_I5xz_Dxy_a = I_NAI_K5xyz_Px_a+ABY*I_NAI_I5xz_Px_a;
  Double I_NAI_I4xyz_Dxy_a = I_NAI_K4x2yz_Px_a+ABY*I_NAI_I4xyz_Px_a;
  Double I_NAI_I4x2z_Dxy_a = I_NAI_K4xy2z_Px_a+ABY*I_NAI_I4x2z_Px_a;
  Double I_NAI_I3x2yz_Dxy_a = I_NAI_K3x3yz_Px_a+ABY*I_NAI_I3x2yz_Px_a;
  Double I_NAI_I3xy2z_Dxy_a = I_NAI_K3x2y2z_Px_a+ABY*I_NAI_I3xy2z_Px_a;
  Double I_NAI_I3x3z_Dxy_a = I_NAI_K3xy3z_Px_a+ABY*I_NAI_I3x3z_Px_a;
  Double I_NAI_I2x3yz_Dxy_a = I_NAI_K2x4yz_Px_a+ABY*I_NAI_I2x3yz_Px_a;
  Double I_NAI_I2x2y2z_Dxy_a = I_NAI_K2x3y2z_Px_a+ABY*I_NAI_I2x2y2z_Px_a;
  Double I_NAI_I2xy3z_Dxy_a = I_NAI_K2x2y3z_Px_a+ABY*I_NAI_I2xy3z_Px_a;
  Double I_NAI_I2x4z_Dxy_a = I_NAI_K2xy4z_Px_a+ABY*I_NAI_I2x4z_Px_a;
  Double I_NAI_Ix4yz_Dxy_a = I_NAI_Kx5yz_Px_a+ABY*I_NAI_Ix4yz_Px_a;
  Double I_NAI_Ix3y2z_Dxy_a = I_NAI_Kx4y2z_Px_a+ABY*I_NAI_Ix3y2z_Px_a;
  Double I_NAI_Ix2y3z_Dxy_a = I_NAI_Kx3y3z_Px_a+ABY*I_NAI_Ix2y3z_Px_a;
  Double I_NAI_Ixy4z_Dxy_a = I_NAI_Kx2y4z_Px_a+ABY*I_NAI_Ixy4z_Px_a;
  Double I_NAI_Ix5z_Dxy_a = I_NAI_Kxy5z_Px_a+ABY*I_NAI_Ix5z_Px_a;
  Double I_NAI_I5yz_Dxy_a = I_NAI_K6yz_Px_a+ABY*I_NAI_I5yz_Px_a;
  Double I_NAI_I4y2z_Dxy_a = I_NAI_K5y2z_Px_a+ABY*I_NAI_I4y2z_Px_a;
  Double I_NAI_I3y3z_Dxy_a = I_NAI_K4y3z_Px_a+ABY*I_NAI_I3y3z_Px_a;
  Double I_NAI_I2y4z_Dxy_a = I_NAI_K3y4z_Px_a+ABY*I_NAI_I2y4z_Px_a;
  Double I_NAI_Iy5z_Dxy_a = I_NAI_K2y5z_Px_a+ABY*I_NAI_Iy5z_Px_a;
  Double I_NAI_I6z_Dxy_a = I_NAI_Ky6z_Px_a+ABY*I_NAI_I6z_Px_a;
  Double I_NAI_I6x_D2y_a = I_NAI_K6xy_Py_a+ABY*I_NAI_I6x_Py_a;
  Double I_NAI_I5xy_D2y_a = I_NAI_K5x2y_Py_a+ABY*I_NAI_I5xy_Py_a;
  Double I_NAI_I5xz_D2y_a = I_NAI_K5xyz_Py_a+ABY*I_NAI_I5xz_Py_a;
  Double I_NAI_I4x2y_D2y_a = I_NAI_K4x3y_Py_a+ABY*I_NAI_I4x2y_Py_a;
  Double I_NAI_I4xyz_D2y_a = I_NAI_K4x2yz_Py_a+ABY*I_NAI_I4xyz_Py_a;
  Double I_NAI_I4x2z_D2y_a = I_NAI_K4xy2z_Py_a+ABY*I_NAI_I4x2z_Py_a;
  Double I_NAI_I3x3y_D2y_a = I_NAI_K3x4y_Py_a+ABY*I_NAI_I3x3y_Py_a;
  Double I_NAI_I3x2yz_D2y_a = I_NAI_K3x3yz_Py_a+ABY*I_NAI_I3x2yz_Py_a;
  Double I_NAI_I3xy2z_D2y_a = I_NAI_K3x2y2z_Py_a+ABY*I_NAI_I3xy2z_Py_a;
  Double I_NAI_I3x3z_D2y_a = I_NAI_K3xy3z_Py_a+ABY*I_NAI_I3x3z_Py_a;
  Double I_NAI_I2x4y_D2y_a = I_NAI_K2x5y_Py_a+ABY*I_NAI_I2x4y_Py_a;
  Double I_NAI_I2x3yz_D2y_a = I_NAI_K2x4yz_Py_a+ABY*I_NAI_I2x3yz_Py_a;
  Double I_NAI_I2x2y2z_D2y_a = I_NAI_K2x3y2z_Py_a+ABY*I_NAI_I2x2y2z_Py_a;
  Double I_NAI_I2xy3z_D2y_a = I_NAI_K2x2y3z_Py_a+ABY*I_NAI_I2xy3z_Py_a;
  Double I_NAI_I2x4z_D2y_a = I_NAI_K2xy4z_Py_a+ABY*I_NAI_I2x4z_Py_a;
  Double I_NAI_Ix5y_D2y_a = I_NAI_Kx6y_Py_a+ABY*I_NAI_Ix5y_Py_a;
  Double I_NAI_Ix4yz_D2y_a = I_NAI_Kx5yz_Py_a+ABY*I_NAI_Ix4yz_Py_a;
  Double I_NAI_Ix3y2z_D2y_a = I_NAI_Kx4y2z_Py_a+ABY*I_NAI_Ix3y2z_Py_a;
  Double I_NAI_Ix2y3z_D2y_a = I_NAI_Kx3y3z_Py_a+ABY*I_NAI_Ix2y3z_Py_a;
  Double I_NAI_Ixy4z_D2y_a = I_NAI_Kx2y4z_Py_a+ABY*I_NAI_Ixy4z_Py_a;
  Double I_NAI_Ix5z_D2y_a = I_NAI_Kxy5z_Py_a+ABY*I_NAI_Ix5z_Py_a;
  Double I_NAI_I6y_D2y_a = I_NAI_K7y_Py_a+ABY*I_NAI_I6y_Py_a;
  Double I_NAI_I5yz_D2y_a = I_NAI_K6yz_Py_a+ABY*I_NAI_I5yz_Py_a;
  Double I_NAI_I4y2z_D2y_a = I_NAI_K5y2z_Py_a+ABY*I_NAI_I4y2z_Py_a;
  Double I_NAI_I3y3z_D2y_a = I_NAI_K4y3z_Py_a+ABY*I_NAI_I3y3z_Py_a;
  Double I_NAI_I2y4z_D2y_a = I_NAI_K3y4z_Py_a+ABY*I_NAI_I2y4z_Py_a;
  Double I_NAI_Iy5z_D2y_a = I_NAI_K2y5z_Py_a+ABY*I_NAI_Iy5z_Py_a;
  Double I_NAI_I6z_D2y_a = I_NAI_Ky6z_Py_a+ABY*I_NAI_I6z_Py_a;
  Double I_NAI_I6x_D2z_a = I_NAI_K6xz_Pz_a+ABZ*I_NAI_I6x_Pz_a;
  Double I_NAI_I5xy_D2z_a = I_NAI_K5xyz_Pz_a+ABZ*I_NAI_I5xy_Pz_a;
  Double I_NAI_I5xz_D2z_a = I_NAI_K5x2z_Pz_a+ABZ*I_NAI_I5xz_Pz_a;
  Double I_NAI_I4x2y_D2z_a = I_NAI_K4x2yz_Pz_a+ABZ*I_NAI_I4x2y_Pz_a;
  Double I_NAI_I4xyz_D2z_a = I_NAI_K4xy2z_Pz_a+ABZ*I_NAI_I4xyz_Pz_a;
  Double I_NAI_I4x2z_D2z_a = I_NAI_K4x3z_Pz_a+ABZ*I_NAI_I4x2z_Pz_a;
  Double I_NAI_I3x3y_D2z_a = I_NAI_K3x3yz_Pz_a+ABZ*I_NAI_I3x3y_Pz_a;
  Double I_NAI_I3x2yz_D2z_a = I_NAI_K3x2y2z_Pz_a+ABZ*I_NAI_I3x2yz_Pz_a;
  Double I_NAI_I3xy2z_D2z_a = I_NAI_K3xy3z_Pz_a+ABZ*I_NAI_I3xy2z_Pz_a;
  Double I_NAI_I3x3z_D2z_a = I_NAI_K3x4z_Pz_a+ABZ*I_NAI_I3x3z_Pz_a;
  Double I_NAI_I2x4y_D2z_a = I_NAI_K2x4yz_Pz_a+ABZ*I_NAI_I2x4y_Pz_a;
  Double I_NAI_I2x3yz_D2z_a = I_NAI_K2x3y2z_Pz_a+ABZ*I_NAI_I2x3yz_Pz_a;
  Double I_NAI_I2x2y2z_D2z_a = I_NAI_K2x2y3z_Pz_a+ABZ*I_NAI_I2x2y2z_Pz_a;
  Double I_NAI_I2xy3z_D2z_a = I_NAI_K2xy4z_Pz_a+ABZ*I_NAI_I2xy3z_Pz_a;
  Double I_NAI_I2x4z_D2z_a = I_NAI_K2x5z_Pz_a+ABZ*I_NAI_I2x4z_Pz_a;
  Double I_NAI_Ix5y_D2z_a = I_NAI_Kx5yz_Pz_a+ABZ*I_NAI_Ix5y_Pz_a;
  Double I_NAI_Ix4yz_D2z_a = I_NAI_Kx4y2z_Pz_a+ABZ*I_NAI_Ix4yz_Pz_a;
  Double I_NAI_Ix3y2z_D2z_a = I_NAI_Kx3y3z_Pz_a+ABZ*I_NAI_Ix3y2z_Pz_a;
  Double I_NAI_Ix2y3z_D2z_a = I_NAI_Kx2y4z_Pz_a+ABZ*I_NAI_Ix2y3z_Pz_a;
  Double I_NAI_Ixy4z_D2z_a = I_NAI_Kxy5z_Pz_a+ABZ*I_NAI_Ixy4z_Pz_a;
  Double I_NAI_Ix5z_D2z_a = I_NAI_Kx6z_Pz_a+ABZ*I_NAI_Ix5z_Pz_a;
  Double I_NAI_I6y_D2z_a = I_NAI_K6yz_Pz_a+ABZ*I_NAI_I6y_Pz_a;
  Double I_NAI_I5yz_D2z_a = I_NAI_K5y2z_Pz_a+ABZ*I_NAI_I5yz_Pz_a;
  Double I_NAI_I4y2z_D2z_a = I_NAI_K4y3z_Pz_a+ABZ*I_NAI_I4y2z_Pz_a;
  Double I_NAI_I3y3z_D2z_a = I_NAI_K3y4z_Pz_a+ABZ*I_NAI_I3y3z_Pz_a;
  Double I_NAI_I2y4z_D2z_a = I_NAI_K2y5z_Pz_a+ABZ*I_NAI_I2y4z_Pz_a;
  Double I_NAI_Iy5z_D2z_a = I_NAI_Ky6z_Pz_a+ABZ*I_NAI_Iy5z_Pz_a;
  Double I_NAI_I6z_D2z_a = I_NAI_K7z_Pz_a+ABZ*I_NAI_I6z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_a
   * RHS shell quartet name: SQ_NAI_H_D_a
   ************************************************************/
  Double I_NAI_H5x_F3x_a = I_NAI_I6x_D2x_a+ABX*I_NAI_H5x_D2x_a;
  Double I_NAI_H4xy_F3x_a = I_NAI_I5xy_D2x_a+ABX*I_NAI_H4xy_D2x_a;
  Double I_NAI_H4xz_F3x_a = I_NAI_I5xz_D2x_a+ABX*I_NAI_H4xz_D2x_a;
  Double I_NAI_H3x2y_F3x_a = I_NAI_I4x2y_D2x_a+ABX*I_NAI_H3x2y_D2x_a;
  Double I_NAI_H3xyz_F3x_a = I_NAI_I4xyz_D2x_a+ABX*I_NAI_H3xyz_D2x_a;
  Double I_NAI_H3x2z_F3x_a = I_NAI_I4x2z_D2x_a+ABX*I_NAI_H3x2z_D2x_a;
  Double I_NAI_H2x3y_F3x_a = I_NAI_I3x3y_D2x_a+ABX*I_NAI_H2x3y_D2x_a;
  Double I_NAI_H2x2yz_F3x_a = I_NAI_I3x2yz_D2x_a+ABX*I_NAI_H2x2yz_D2x_a;
  Double I_NAI_H2xy2z_F3x_a = I_NAI_I3xy2z_D2x_a+ABX*I_NAI_H2xy2z_D2x_a;
  Double I_NAI_H2x3z_F3x_a = I_NAI_I3x3z_D2x_a+ABX*I_NAI_H2x3z_D2x_a;
  Double I_NAI_Hx4y_F3x_a = I_NAI_I2x4y_D2x_a+ABX*I_NAI_Hx4y_D2x_a;
  Double I_NAI_Hx3yz_F3x_a = I_NAI_I2x3yz_D2x_a+ABX*I_NAI_Hx3yz_D2x_a;
  Double I_NAI_Hx2y2z_F3x_a = I_NAI_I2x2y2z_D2x_a+ABX*I_NAI_Hx2y2z_D2x_a;
  Double I_NAI_Hxy3z_F3x_a = I_NAI_I2xy3z_D2x_a+ABX*I_NAI_Hxy3z_D2x_a;
  Double I_NAI_Hx4z_F3x_a = I_NAI_I2x4z_D2x_a+ABX*I_NAI_Hx4z_D2x_a;
  Double I_NAI_H5y_F3x_a = I_NAI_Ix5y_D2x_a+ABX*I_NAI_H5y_D2x_a;
  Double I_NAI_H4yz_F3x_a = I_NAI_Ix4yz_D2x_a+ABX*I_NAI_H4yz_D2x_a;
  Double I_NAI_H3y2z_F3x_a = I_NAI_Ix3y2z_D2x_a+ABX*I_NAI_H3y2z_D2x_a;
  Double I_NAI_H2y3z_F3x_a = I_NAI_Ix2y3z_D2x_a+ABX*I_NAI_H2y3z_D2x_a;
  Double I_NAI_Hy4z_F3x_a = I_NAI_Ixy4z_D2x_a+ABX*I_NAI_Hy4z_D2x_a;
  Double I_NAI_H5z_F3x_a = I_NAI_Ix5z_D2x_a+ABX*I_NAI_H5z_D2x_a;
  Double I_NAI_H5x_F2xy_a = I_NAI_I5xy_D2x_a+ABY*I_NAI_H5x_D2x_a;
  Double I_NAI_H4xy_F2xy_a = I_NAI_I4x2y_D2x_a+ABY*I_NAI_H4xy_D2x_a;
  Double I_NAI_H4xz_F2xy_a = I_NAI_I4xyz_D2x_a+ABY*I_NAI_H4xz_D2x_a;
  Double I_NAI_H3x2y_F2xy_a = I_NAI_I3x3y_D2x_a+ABY*I_NAI_H3x2y_D2x_a;
  Double I_NAI_H3xyz_F2xy_a = I_NAI_I3x2yz_D2x_a+ABY*I_NAI_H3xyz_D2x_a;
  Double I_NAI_H3x2z_F2xy_a = I_NAI_I3xy2z_D2x_a+ABY*I_NAI_H3x2z_D2x_a;
  Double I_NAI_H2x3y_F2xy_a = I_NAI_I2x4y_D2x_a+ABY*I_NAI_H2x3y_D2x_a;
  Double I_NAI_H2x2yz_F2xy_a = I_NAI_I2x3yz_D2x_a+ABY*I_NAI_H2x2yz_D2x_a;
  Double I_NAI_H2xy2z_F2xy_a = I_NAI_I2x2y2z_D2x_a+ABY*I_NAI_H2xy2z_D2x_a;
  Double I_NAI_H2x3z_F2xy_a = I_NAI_I2xy3z_D2x_a+ABY*I_NAI_H2x3z_D2x_a;
  Double I_NAI_Hx4y_F2xy_a = I_NAI_Ix5y_D2x_a+ABY*I_NAI_Hx4y_D2x_a;
  Double I_NAI_Hx3yz_F2xy_a = I_NAI_Ix4yz_D2x_a+ABY*I_NAI_Hx3yz_D2x_a;
  Double I_NAI_Hx2y2z_F2xy_a = I_NAI_Ix3y2z_D2x_a+ABY*I_NAI_Hx2y2z_D2x_a;
  Double I_NAI_Hxy3z_F2xy_a = I_NAI_Ix2y3z_D2x_a+ABY*I_NAI_Hxy3z_D2x_a;
  Double I_NAI_Hx4z_F2xy_a = I_NAI_Ixy4z_D2x_a+ABY*I_NAI_Hx4z_D2x_a;
  Double I_NAI_H5y_F2xy_a = I_NAI_I6y_D2x_a+ABY*I_NAI_H5y_D2x_a;
  Double I_NAI_H4yz_F2xy_a = I_NAI_I5yz_D2x_a+ABY*I_NAI_H4yz_D2x_a;
  Double I_NAI_H3y2z_F2xy_a = I_NAI_I4y2z_D2x_a+ABY*I_NAI_H3y2z_D2x_a;
  Double I_NAI_H2y3z_F2xy_a = I_NAI_I3y3z_D2x_a+ABY*I_NAI_H2y3z_D2x_a;
  Double I_NAI_Hy4z_F2xy_a = I_NAI_I2y4z_D2x_a+ABY*I_NAI_Hy4z_D2x_a;
  Double I_NAI_H5z_F2xy_a = I_NAI_Iy5z_D2x_a+ABY*I_NAI_H5z_D2x_a;
  Double I_NAI_H5x_F2xz_a = I_NAI_I5xz_D2x_a+ABZ*I_NAI_H5x_D2x_a;
  Double I_NAI_H4xy_F2xz_a = I_NAI_I4xyz_D2x_a+ABZ*I_NAI_H4xy_D2x_a;
  Double I_NAI_H4xz_F2xz_a = I_NAI_I4x2z_D2x_a+ABZ*I_NAI_H4xz_D2x_a;
  Double I_NAI_H3x2y_F2xz_a = I_NAI_I3x2yz_D2x_a+ABZ*I_NAI_H3x2y_D2x_a;
  Double I_NAI_H3xyz_F2xz_a = I_NAI_I3xy2z_D2x_a+ABZ*I_NAI_H3xyz_D2x_a;
  Double I_NAI_H3x2z_F2xz_a = I_NAI_I3x3z_D2x_a+ABZ*I_NAI_H3x2z_D2x_a;
  Double I_NAI_H2x3y_F2xz_a = I_NAI_I2x3yz_D2x_a+ABZ*I_NAI_H2x3y_D2x_a;
  Double I_NAI_H2x2yz_F2xz_a = I_NAI_I2x2y2z_D2x_a+ABZ*I_NAI_H2x2yz_D2x_a;
  Double I_NAI_H2xy2z_F2xz_a = I_NAI_I2xy3z_D2x_a+ABZ*I_NAI_H2xy2z_D2x_a;
  Double I_NAI_H2x3z_F2xz_a = I_NAI_I2x4z_D2x_a+ABZ*I_NAI_H2x3z_D2x_a;
  Double I_NAI_Hx4y_F2xz_a = I_NAI_Ix4yz_D2x_a+ABZ*I_NAI_Hx4y_D2x_a;
  Double I_NAI_Hx3yz_F2xz_a = I_NAI_Ix3y2z_D2x_a+ABZ*I_NAI_Hx3yz_D2x_a;
  Double I_NAI_Hx2y2z_F2xz_a = I_NAI_Ix2y3z_D2x_a+ABZ*I_NAI_Hx2y2z_D2x_a;
  Double I_NAI_Hxy3z_F2xz_a = I_NAI_Ixy4z_D2x_a+ABZ*I_NAI_Hxy3z_D2x_a;
  Double I_NAI_Hx4z_F2xz_a = I_NAI_Ix5z_D2x_a+ABZ*I_NAI_Hx4z_D2x_a;
  Double I_NAI_H5y_F2xz_a = I_NAI_I5yz_D2x_a+ABZ*I_NAI_H5y_D2x_a;
  Double I_NAI_H4yz_F2xz_a = I_NAI_I4y2z_D2x_a+ABZ*I_NAI_H4yz_D2x_a;
  Double I_NAI_H3y2z_F2xz_a = I_NAI_I3y3z_D2x_a+ABZ*I_NAI_H3y2z_D2x_a;
  Double I_NAI_H2y3z_F2xz_a = I_NAI_I2y4z_D2x_a+ABZ*I_NAI_H2y3z_D2x_a;
  Double I_NAI_Hy4z_F2xz_a = I_NAI_Iy5z_D2x_a+ABZ*I_NAI_Hy4z_D2x_a;
  Double I_NAI_H5z_F2xz_a = I_NAI_I6z_D2x_a+ABZ*I_NAI_H5z_D2x_a;
  Double I_NAI_H5x_Fx2y_a = I_NAI_I6x_D2y_a+ABX*I_NAI_H5x_D2y_a;
  Double I_NAI_H4xy_Fx2y_a = I_NAI_I5xy_D2y_a+ABX*I_NAI_H4xy_D2y_a;
  Double I_NAI_H4xz_Fx2y_a = I_NAI_I5xz_D2y_a+ABX*I_NAI_H4xz_D2y_a;
  Double I_NAI_H3x2y_Fx2y_a = I_NAI_I4x2y_D2y_a+ABX*I_NAI_H3x2y_D2y_a;
  Double I_NAI_H3xyz_Fx2y_a = I_NAI_I4xyz_D2y_a+ABX*I_NAI_H3xyz_D2y_a;
  Double I_NAI_H3x2z_Fx2y_a = I_NAI_I4x2z_D2y_a+ABX*I_NAI_H3x2z_D2y_a;
  Double I_NAI_H2x3y_Fx2y_a = I_NAI_I3x3y_D2y_a+ABX*I_NAI_H2x3y_D2y_a;
  Double I_NAI_H2x2yz_Fx2y_a = I_NAI_I3x2yz_D2y_a+ABX*I_NAI_H2x2yz_D2y_a;
  Double I_NAI_H2xy2z_Fx2y_a = I_NAI_I3xy2z_D2y_a+ABX*I_NAI_H2xy2z_D2y_a;
  Double I_NAI_H2x3z_Fx2y_a = I_NAI_I3x3z_D2y_a+ABX*I_NAI_H2x3z_D2y_a;
  Double I_NAI_Hx4y_Fx2y_a = I_NAI_I2x4y_D2y_a+ABX*I_NAI_Hx4y_D2y_a;
  Double I_NAI_Hx3yz_Fx2y_a = I_NAI_I2x3yz_D2y_a+ABX*I_NAI_Hx3yz_D2y_a;
  Double I_NAI_Hx2y2z_Fx2y_a = I_NAI_I2x2y2z_D2y_a+ABX*I_NAI_Hx2y2z_D2y_a;
  Double I_NAI_Hxy3z_Fx2y_a = I_NAI_I2xy3z_D2y_a+ABX*I_NAI_Hxy3z_D2y_a;
  Double I_NAI_Hx4z_Fx2y_a = I_NAI_I2x4z_D2y_a+ABX*I_NAI_Hx4z_D2y_a;
  Double I_NAI_H5y_Fx2y_a = I_NAI_Ix5y_D2y_a+ABX*I_NAI_H5y_D2y_a;
  Double I_NAI_H4yz_Fx2y_a = I_NAI_Ix4yz_D2y_a+ABX*I_NAI_H4yz_D2y_a;
  Double I_NAI_H3y2z_Fx2y_a = I_NAI_Ix3y2z_D2y_a+ABX*I_NAI_H3y2z_D2y_a;
  Double I_NAI_H2y3z_Fx2y_a = I_NAI_Ix2y3z_D2y_a+ABX*I_NAI_H2y3z_D2y_a;
  Double I_NAI_Hy4z_Fx2y_a = I_NAI_Ixy4z_D2y_a+ABX*I_NAI_Hy4z_D2y_a;
  Double I_NAI_H5z_Fx2y_a = I_NAI_Ix5z_D2y_a+ABX*I_NAI_H5z_D2y_a;
  Double I_NAI_H5x_Fxyz_a = I_NAI_I5xz_Dxy_a+ABZ*I_NAI_H5x_Dxy_a;
  Double I_NAI_H4xy_Fxyz_a = I_NAI_I4xyz_Dxy_a+ABZ*I_NAI_H4xy_Dxy_a;
  Double I_NAI_H4xz_Fxyz_a = I_NAI_I4x2z_Dxy_a+ABZ*I_NAI_H4xz_Dxy_a;
  Double I_NAI_H3x2y_Fxyz_a = I_NAI_I3x2yz_Dxy_a+ABZ*I_NAI_H3x2y_Dxy_a;
  Double I_NAI_H3xyz_Fxyz_a = I_NAI_I3xy2z_Dxy_a+ABZ*I_NAI_H3xyz_Dxy_a;
  Double I_NAI_H3x2z_Fxyz_a = I_NAI_I3x3z_Dxy_a+ABZ*I_NAI_H3x2z_Dxy_a;
  Double I_NAI_H2x3y_Fxyz_a = I_NAI_I2x3yz_Dxy_a+ABZ*I_NAI_H2x3y_Dxy_a;
  Double I_NAI_H2x2yz_Fxyz_a = I_NAI_I2x2y2z_Dxy_a+ABZ*I_NAI_H2x2yz_Dxy_a;
  Double I_NAI_H2xy2z_Fxyz_a = I_NAI_I2xy3z_Dxy_a+ABZ*I_NAI_H2xy2z_Dxy_a;
  Double I_NAI_H2x3z_Fxyz_a = I_NAI_I2x4z_Dxy_a+ABZ*I_NAI_H2x3z_Dxy_a;
  Double I_NAI_Hx4y_Fxyz_a = I_NAI_Ix4yz_Dxy_a+ABZ*I_NAI_Hx4y_Dxy_a;
  Double I_NAI_Hx3yz_Fxyz_a = I_NAI_Ix3y2z_Dxy_a+ABZ*I_NAI_Hx3yz_Dxy_a;
  Double I_NAI_Hx2y2z_Fxyz_a = I_NAI_Ix2y3z_Dxy_a+ABZ*I_NAI_Hx2y2z_Dxy_a;
  Double I_NAI_Hxy3z_Fxyz_a = I_NAI_Ixy4z_Dxy_a+ABZ*I_NAI_Hxy3z_Dxy_a;
  Double I_NAI_Hx4z_Fxyz_a = I_NAI_Ix5z_Dxy_a+ABZ*I_NAI_Hx4z_Dxy_a;
  Double I_NAI_H5y_Fxyz_a = I_NAI_I5yz_Dxy_a+ABZ*I_NAI_H5y_Dxy_a;
  Double I_NAI_H4yz_Fxyz_a = I_NAI_I4y2z_Dxy_a+ABZ*I_NAI_H4yz_Dxy_a;
  Double I_NAI_H3y2z_Fxyz_a = I_NAI_I3y3z_Dxy_a+ABZ*I_NAI_H3y2z_Dxy_a;
  Double I_NAI_H2y3z_Fxyz_a = I_NAI_I2y4z_Dxy_a+ABZ*I_NAI_H2y3z_Dxy_a;
  Double I_NAI_Hy4z_Fxyz_a = I_NAI_Iy5z_Dxy_a+ABZ*I_NAI_Hy4z_Dxy_a;
  Double I_NAI_H5z_Fxyz_a = I_NAI_I6z_Dxy_a+ABZ*I_NAI_H5z_Dxy_a;
  Double I_NAI_H5x_Fx2z_a = I_NAI_I6x_D2z_a+ABX*I_NAI_H5x_D2z_a;
  Double I_NAI_H4xy_Fx2z_a = I_NAI_I5xy_D2z_a+ABX*I_NAI_H4xy_D2z_a;
  Double I_NAI_H4xz_Fx2z_a = I_NAI_I5xz_D2z_a+ABX*I_NAI_H4xz_D2z_a;
  Double I_NAI_H3x2y_Fx2z_a = I_NAI_I4x2y_D2z_a+ABX*I_NAI_H3x2y_D2z_a;
  Double I_NAI_H3xyz_Fx2z_a = I_NAI_I4xyz_D2z_a+ABX*I_NAI_H3xyz_D2z_a;
  Double I_NAI_H3x2z_Fx2z_a = I_NAI_I4x2z_D2z_a+ABX*I_NAI_H3x2z_D2z_a;
  Double I_NAI_H2x3y_Fx2z_a = I_NAI_I3x3y_D2z_a+ABX*I_NAI_H2x3y_D2z_a;
  Double I_NAI_H2x2yz_Fx2z_a = I_NAI_I3x2yz_D2z_a+ABX*I_NAI_H2x2yz_D2z_a;
  Double I_NAI_H2xy2z_Fx2z_a = I_NAI_I3xy2z_D2z_a+ABX*I_NAI_H2xy2z_D2z_a;
  Double I_NAI_H2x3z_Fx2z_a = I_NAI_I3x3z_D2z_a+ABX*I_NAI_H2x3z_D2z_a;
  Double I_NAI_Hx4y_Fx2z_a = I_NAI_I2x4y_D2z_a+ABX*I_NAI_Hx4y_D2z_a;
  Double I_NAI_Hx3yz_Fx2z_a = I_NAI_I2x3yz_D2z_a+ABX*I_NAI_Hx3yz_D2z_a;
  Double I_NAI_Hx2y2z_Fx2z_a = I_NAI_I2x2y2z_D2z_a+ABX*I_NAI_Hx2y2z_D2z_a;
  Double I_NAI_Hxy3z_Fx2z_a = I_NAI_I2xy3z_D2z_a+ABX*I_NAI_Hxy3z_D2z_a;
  Double I_NAI_Hx4z_Fx2z_a = I_NAI_I2x4z_D2z_a+ABX*I_NAI_Hx4z_D2z_a;
  Double I_NAI_H5y_Fx2z_a = I_NAI_Ix5y_D2z_a+ABX*I_NAI_H5y_D2z_a;
  Double I_NAI_H4yz_Fx2z_a = I_NAI_Ix4yz_D2z_a+ABX*I_NAI_H4yz_D2z_a;
  Double I_NAI_H3y2z_Fx2z_a = I_NAI_Ix3y2z_D2z_a+ABX*I_NAI_H3y2z_D2z_a;
  Double I_NAI_H2y3z_Fx2z_a = I_NAI_Ix2y3z_D2z_a+ABX*I_NAI_H2y3z_D2z_a;
  Double I_NAI_Hy4z_Fx2z_a = I_NAI_Ixy4z_D2z_a+ABX*I_NAI_Hy4z_D2z_a;
  Double I_NAI_H5z_Fx2z_a = I_NAI_Ix5z_D2z_a+ABX*I_NAI_H5z_D2z_a;
  Double I_NAI_H5x_F3y_a = I_NAI_I5xy_D2y_a+ABY*I_NAI_H5x_D2y_a;
  Double I_NAI_H4xy_F3y_a = I_NAI_I4x2y_D2y_a+ABY*I_NAI_H4xy_D2y_a;
  Double I_NAI_H4xz_F3y_a = I_NAI_I4xyz_D2y_a+ABY*I_NAI_H4xz_D2y_a;
  Double I_NAI_H3x2y_F3y_a = I_NAI_I3x3y_D2y_a+ABY*I_NAI_H3x2y_D2y_a;
  Double I_NAI_H3xyz_F3y_a = I_NAI_I3x2yz_D2y_a+ABY*I_NAI_H3xyz_D2y_a;
  Double I_NAI_H3x2z_F3y_a = I_NAI_I3xy2z_D2y_a+ABY*I_NAI_H3x2z_D2y_a;
  Double I_NAI_H2x3y_F3y_a = I_NAI_I2x4y_D2y_a+ABY*I_NAI_H2x3y_D2y_a;
  Double I_NAI_H2x2yz_F3y_a = I_NAI_I2x3yz_D2y_a+ABY*I_NAI_H2x2yz_D2y_a;
  Double I_NAI_H2xy2z_F3y_a = I_NAI_I2x2y2z_D2y_a+ABY*I_NAI_H2xy2z_D2y_a;
  Double I_NAI_H2x3z_F3y_a = I_NAI_I2xy3z_D2y_a+ABY*I_NAI_H2x3z_D2y_a;
  Double I_NAI_Hx4y_F3y_a = I_NAI_Ix5y_D2y_a+ABY*I_NAI_Hx4y_D2y_a;
  Double I_NAI_Hx3yz_F3y_a = I_NAI_Ix4yz_D2y_a+ABY*I_NAI_Hx3yz_D2y_a;
  Double I_NAI_Hx2y2z_F3y_a = I_NAI_Ix3y2z_D2y_a+ABY*I_NAI_Hx2y2z_D2y_a;
  Double I_NAI_Hxy3z_F3y_a = I_NAI_Ix2y3z_D2y_a+ABY*I_NAI_Hxy3z_D2y_a;
  Double I_NAI_Hx4z_F3y_a = I_NAI_Ixy4z_D2y_a+ABY*I_NAI_Hx4z_D2y_a;
  Double I_NAI_H5y_F3y_a = I_NAI_I6y_D2y_a+ABY*I_NAI_H5y_D2y_a;
  Double I_NAI_H4yz_F3y_a = I_NAI_I5yz_D2y_a+ABY*I_NAI_H4yz_D2y_a;
  Double I_NAI_H3y2z_F3y_a = I_NAI_I4y2z_D2y_a+ABY*I_NAI_H3y2z_D2y_a;
  Double I_NAI_H2y3z_F3y_a = I_NAI_I3y3z_D2y_a+ABY*I_NAI_H2y3z_D2y_a;
  Double I_NAI_Hy4z_F3y_a = I_NAI_I2y4z_D2y_a+ABY*I_NAI_Hy4z_D2y_a;
  Double I_NAI_H5z_F3y_a = I_NAI_Iy5z_D2y_a+ABY*I_NAI_H5z_D2y_a;
  Double I_NAI_H5x_F2yz_a = I_NAI_I5xz_D2y_a+ABZ*I_NAI_H5x_D2y_a;
  Double I_NAI_H4xy_F2yz_a = I_NAI_I4xyz_D2y_a+ABZ*I_NAI_H4xy_D2y_a;
  Double I_NAI_H4xz_F2yz_a = I_NAI_I4x2z_D2y_a+ABZ*I_NAI_H4xz_D2y_a;
  Double I_NAI_H3x2y_F2yz_a = I_NAI_I3x2yz_D2y_a+ABZ*I_NAI_H3x2y_D2y_a;
  Double I_NAI_H3xyz_F2yz_a = I_NAI_I3xy2z_D2y_a+ABZ*I_NAI_H3xyz_D2y_a;
  Double I_NAI_H3x2z_F2yz_a = I_NAI_I3x3z_D2y_a+ABZ*I_NAI_H3x2z_D2y_a;
  Double I_NAI_H2x3y_F2yz_a = I_NAI_I2x3yz_D2y_a+ABZ*I_NAI_H2x3y_D2y_a;
  Double I_NAI_H2x2yz_F2yz_a = I_NAI_I2x2y2z_D2y_a+ABZ*I_NAI_H2x2yz_D2y_a;
  Double I_NAI_H2xy2z_F2yz_a = I_NAI_I2xy3z_D2y_a+ABZ*I_NAI_H2xy2z_D2y_a;
  Double I_NAI_H2x3z_F2yz_a = I_NAI_I2x4z_D2y_a+ABZ*I_NAI_H2x3z_D2y_a;
  Double I_NAI_Hx4y_F2yz_a = I_NAI_Ix4yz_D2y_a+ABZ*I_NAI_Hx4y_D2y_a;
  Double I_NAI_Hx3yz_F2yz_a = I_NAI_Ix3y2z_D2y_a+ABZ*I_NAI_Hx3yz_D2y_a;
  Double I_NAI_Hx2y2z_F2yz_a = I_NAI_Ix2y3z_D2y_a+ABZ*I_NAI_Hx2y2z_D2y_a;
  Double I_NAI_Hxy3z_F2yz_a = I_NAI_Ixy4z_D2y_a+ABZ*I_NAI_Hxy3z_D2y_a;
  Double I_NAI_Hx4z_F2yz_a = I_NAI_Ix5z_D2y_a+ABZ*I_NAI_Hx4z_D2y_a;
  Double I_NAI_H5y_F2yz_a = I_NAI_I5yz_D2y_a+ABZ*I_NAI_H5y_D2y_a;
  Double I_NAI_H4yz_F2yz_a = I_NAI_I4y2z_D2y_a+ABZ*I_NAI_H4yz_D2y_a;
  Double I_NAI_H3y2z_F2yz_a = I_NAI_I3y3z_D2y_a+ABZ*I_NAI_H3y2z_D2y_a;
  Double I_NAI_H2y3z_F2yz_a = I_NAI_I2y4z_D2y_a+ABZ*I_NAI_H2y3z_D2y_a;
  Double I_NAI_Hy4z_F2yz_a = I_NAI_Iy5z_D2y_a+ABZ*I_NAI_Hy4z_D2y_a;
  Double I_NAI_H5z_F2yz_a = I_NAI_I6z_D2y_a+ABZ*I_NAI_H5z_D2y_a;
  Double I_NAI_H5x_Fy2z_a = I_NAI_I5xy_D2z_a+ABY*I_NAI_H5x_D2z_a;
  Double I_NAI_H4xy_Fy2z_a = I_NAI_I4x2y_D2z_a+ABY*I_NAI_H4xy_D2z_a;
  Double I_NAI_H4xz_Fy2z_a = I_NAI_I4xyz_D2z_a+ABY*I_NAI_H4xz_D2z_a;
  Double I_NAI_H3x2y_Fy2z_a = I_NAI_I3x3y_D2z_a+ABY*I_NAI_H3x2y_D2z_a;
  Double I_NAI_H3xyz_Fy2z_a = I_NAI_I3x2yz_D2z_a+ABY*I_NAI_H3xyz_D2z_a;
  Double I_NAI_H3x2z_Fy2z_a = I_NAI_I3xy2z_D2z_a+ABY*I_NAI_H3x2z_D2z_a;
  Double I_NAI_H2x3y_Fy2z_a = I_NAI_I2x4y_D2z_a+ABY*I_NAI_H2x3y_D2z_a;
  Double I_NAI_H2x2yz_Fy2z_a = I_NAI_I2x3yz_D2z_a+ABY*I_NAI_H2x2yz_D2z_a;
  Double I_NAI_H2xy2z_Fy2z_a = I_NAI_I2x2y2z_D2z_a+ABY*I_NAI_H2xy2z_D2z_a;
  Double I_NAI_H2x3z_Fy2z_a = I_NAI_I2xy3z_D2z_a+ABY*I_NAI_H2x3z_D2z_a;
  Double I_NAI_Hx4y_Fy2z_a = I_NAI_Ix5y_D2z_a+ABY*I_NAI_Hx4y_D2z_a;
  Double I_NAI_Hx3yz_Fy2z_a = I_NAI_Ix4yz_D2z_a+ABY*I_NAI_Hx3yz_D2z_a;
  Double I_NAI_Hx2y2z_Fy2z_a = I_NAI_Ix3y2z_D2z_a+ABY*I_NAI_Hx2y2z_D2z_a;
  Double I_NAI_Hxy3z_Fy2z_a = I_NAI_Ix2y3z_D2z_a+ABY*I_NAI_Hxy3z_D2z_a;
  Double I_NAI_Hx4z_Fy2z_a = I_NAI_Ixy4z_D2z_a+ABY*I_NAI_Hx4z_D2z_a;
  Double I_NAI_H5y_Fy2z_a = I_NAI_I6y_D2z_a+ABY*I_NAI_H5y_D2z_a;
  Double I_NAI_H4yz_Fy2z_a = I_NAI_I5yz_D2z_a+ABY*I_NAI_H4yz_D2z_a;
  Double I_NAI_H3y2z_Fy2z_a = I_NAI_I4y2z_D2z_a+ABY*I_NAI_H3y2z_D2z_a;
  Double I_NAI_H2y3z_Fy2z_a = I_NAI_I3y3z_D2z_a+ABY*I_NAI_H2y3z_D2z_a;
  Double I_NAI_Hy4z_Fy2z_a = I_NAI_I2y4z_D2z_a+ABY*I_NAI_Hy4z_D2z_a;
  Double I_NAI_H5z_Fy2z_a = I_NAI_Iy5z_D2z_a+ABY*I_NAI_H5z_D2z_a;
  Double I_NAI_H5x_F3z_a = I_NAI_I5xz_D2z_a+ABZ*I_NAI_H5x_D2z_a;
  Double I_NAI_H4xy_F3z_a = I_NAI_I4xyz_D2z_a+ABZ*I_NAI_H4xy_D2z_a;
  Double I_NAI_H4xz_F3z_a = I_NAI_I4x2z_D2z_a+ABZ*I_NAI_H4xz_D2z_a;
  Double I_NAI_H3x2y_F3z_a = I_NAI_I3x2yz_D2z_a+ABZ*I_NAI_H3x2y_D2z_a;
  Double I_NAI_H3xyz_F3z_a = I_NAI_I3xy2z_D2z_a+ABZ*I_NAI_H3xyz_D2z_a;
  Double I_NAI_H3x2z_F3z_a = I_NAI_I3x3z_D2z_a+ABZ*I_NAI_H3x2z_D2z_a;
  Double I_NAI_H2x3y_F3z_a = I_NAI_I2x3yz_D2z_a+ABZ*I_NAI_H2x3y_D2z_a;
  Double I_NAI_H2x2yz_F3z_a = I_NAI_I2x2y2z_D2z_a+ABZ*I_NAI_H2x2yz_D2z_a;
  Double I_NAI_H2xy2z_F3z_a = I_NAI_I2xy3z_D2z_a+ABZ*I_NAI_H2xy2z_D2z_a;
  Double I_NAI_H2x3z_F3z_a = I_NAI_I2x4z_D2z_a+ABZ*I_NAI_H2x3z_D2z_a;
  Double I_NAI_Hx4y_F3z_a = I_NAI_Ix4yz_D2z_a+ABZ*I_NAI_Hx4y_D2z_a;
  Double I_NAI_Hx3yz_F3z_a = I_NAI_Ix3y2z_D2z_a+ABZ*I_NAI_Hx3yz_D2z_a;
  Double I_NAI_Hx2y2z_F3z_a = I_NAI_Ix2y3z_D2z_a+ABZ*I_NAI_Hx2y2z_D2z_a;
  Double I_NAI_Hxy3z_F3z_a = I_NAI_Ixy4z_D2z_a+ABZ*I_NAI_Hxy3z_D2z_a;
  Double I_NAI_Hx4z_F3z_a = I_NAI_Ix5z_D2z_a+ABZ*I_NAI_Hx4z_D2z_a;
  Double I_NAI_H5y_F3z_a = I_NAI_I5yz_D2z_a+ABZ*I_NAI_H5y_D2z_a;
  Double I_NAI_H4yz_F3z_a = I_NAI_I4y2z_D2z_a+ABZ*I_NAI_H4yz_D2z_a;
  Double I_NAI_H3y2z_F3z_a = I_NAI_I3y3z_D2z_a+ABZ*I_NAI_H3y2z_D2z_a;
  Double I_NAI_H2y3z_F3z_a = I_NAI_I2y4z_D2z_a+ABZ*I_NAI_H2y3z_D2z_a;
  Double I_NAI_Hy4z_F3z_a = I_NAI_Iy5z_D2z_a+ABZ*I_NAI_Hy4z_D2z_a;
  Double I_NAI_H5z_F3z_a = I_NAI_I6z_D2z_a+ABZ*I_NAI_H5z_D2z_a;

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
   * totally 45 integrals are omitted 
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
   * shell quartet name: SQ_NAI_I_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
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
  Double I_NAI_I6y_Px_b = I_NAI_Kx6y_S_b+ABX*I_NAI_I6y_S_b;
  Double I_NAI_I5yz_Px_b = I_NAI_Kx5yz_S_b+ABX*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Px_b = I_NAI_Kx4y2z_S_b+ABX*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Px_b = I_NAI_Kx3y3z_S_b+ABX*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Px_b = I_NAI_Kx2y4z_S_b+ABX*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Px_b = I_NAI_Kxy5z_S_b+ABX*I_NAI_Iy5z_S_b;
  Double I_NAI_I6z_Px_b = I_NAI_Kx6z_S_b+ABX*I_NAI_I6z_S_b;
  Double I_NAI_I5xy_Py_b = I_NAI_K5x2y_S_b+ABY*I_NAI_I5xy_S_b;
  Double I_NAI_I5xz_Py_b = I_NAI_K5xyz_S_b+ABY*I_NAI_I5xz_S_b;
  Double I_NAI_I4x2y_Py_b = I_NAI_K4x3y_S_b+ABY*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Py_b = I_NAI_K4x2yz_S_b+ABY*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Py_b = I_NAI_K4xy2z_S_b+ABY*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x3y_Py_b = I_NAI_K3x4y_S_b+ABY*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Py_b = I_NAI_K3x3yz_S_b+ABY*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Py_b = I_NAI_K3x2y2z_S_b+ABY*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Py_b = I_NAI_K3xy3z_S_b+ABY*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x4y_Py_b = I_NAI_K2x5y_S_b+ABY*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Py_b = I_NAI_K2x4yz_S_b+ABY*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Py_b = I_NAI_K2x3y2z_S_b+ABY*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Py_b = I_NAI_K2x2y3z_S_b+ABY*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Py_b = I_NAI_K2xy4z_S_b+ABY*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix5y_Py_b = I_NAI_Kx6y_S_b+ABY*I_NAI_Ix5y_S_b;
  Double I_NAI_Ix4yz_Py_b = I_NAI_Kx5yz_S_b+ABY*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Py_b = I_NAI_Kx4y2z_S_b+ABY*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Py_b = I_NAI_Kx3y3z_S_b+ABY*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Py_b = I_NAI_Kx2y4z_S_b+ABY*I_NAI_Ixy4z_S_b;
  Double I_NAI_Ix5z_Py_b = I_NAI_Kxy5z_S_b+ABY*I_NAI_Ix5z_S_b;
  Double I_NAI_I6y_Py_b = I_NAI_K7y_S_b+ABY*I_NAI_I6y_S_b;
  Double I_NAI_I5yz_Py_b = I_NAI_K6yz_S_b+ABY*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Py_b = I_NAI_K5y2z_S_b+ABY*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Py_b = I_NAI_K4y3z_S_b+ABY*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Py_b = I_NAI_K3y4z_S_b+ABY*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Py_b = I_NAI_K2y5z_S_b+ABY*I_NAI_Iy5z_S_b;
  Double I_NAI_I6z_Py_b = I_NAI_Ky6z_S_b+ABY*I_NAI_I6z_S_b;
  Double I_NAI_I5xy_Pz_b = I_NAI_K5xyz_S_b+ABZ*I_NAI_I5xy_S_b;
  Double I_NAI_I5xz_Pz_b = I_NAI_K5x2z_S_b+ABZ*I_NAI_I5xz_S_b;
  Double I_NAI_I4x2y_Pz_b = I_NAI_K4x2yz_S_b+ABZ*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Pz_b = I_NAI_K4xy2z_S_b+ABZ*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Pz_b = I_NAI_K4x3z_S_b+ABZ*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x3y_Pz_b = I_NAI_K3x3yz_S_b+ABZ*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Pz_b = I_NAI_K3x2y2z_S_b+ABZ*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Pz_b = I_NAI_K3xy3z_S_b+ABZ*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Pz_b = I_NAI_K3x4z_S_b+ABZ*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x4y_Pz_b = I_NAI_K2x4yz_S_b+ABZ*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Pz_b = I_NAI_K2x3y2z_S_b+ABZ*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Pz_b = I_NAI_K2x2y3z_S_b+ABZ*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Pz_b = I_NAI_K2xy4z_S_b+ABZ*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Pz_b = I_NAI_K2x5z_S_b+ABZ*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix5y_Pz_b = I_NAI_Kx5yz_S_b+ABZ*I_NAI_Ix5y_S_b;
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
   * totally 63 integrals are omitted 
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
   * totally 30 integrals are omitted 
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
   * shell quartet name: SQ_NAI_K_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 27 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S_b
   * RHS shell quartet name: SQ_NAI_K_S_b
   ************************************************************/
  Double I_NAI_K7x_Px_b = I_NAI_L8x_S_b+ABX*I_NAI_K7x_S_b;
  Double I_NAI_K6xy_Px_b = I_NAI_L7xy_S_b+ABX*I_NAI_K6xy_S_b;
  Double I_NAI_K6xz_Px_b = I_NAI_L7xz_S_b+ABX*I_NAI_K6xz_S_b;
  Double I_NAI_K5x2y_Px_b = I_NAI_L6x2y_S_b+ABX*I_NAI_K5x2y_S_b;
  Double I_NAI_K5xyz_Px_b = I_NAI_L6xyz_S_b+ABX*I_NAI_K5xyz_S_b;
  Double I_NAI_K5x2z_Px_b = I_NAI_L6x2z_S_b+ABX*I_NAI_K5x2z_S_b;
  Double I_NAI_K4x3y_Px_b = I_NAI_L5x3y_S_b+ABX*I_NAI_K4x3y_S_b;
  Double I_NAI_K4x2yz_Px_b = I_NAI_L5x2yz_S_b+ABX*I_NAI_K4x2yz_S_b;
  Double I_NAI_K4xy2z_Px_b = I_NAI_L5xy2z_S_b+ABX*I_NAI_K4xy2z_S_b;
  Double I_NAI_K4x3z_Px_b = I_NAI_L5x3z_S_b+ABX*I_NAI_K4x3z_S_b;
  Double I_NAI_K3x4y_Px_b = I_NAI_L4x4y_S_b+ABX*I_NAI_K3x4y_S_b;
  Double I_NAI_K3x3yz_Px_b = I_NAI_L4x3yz_S_b+ABX*I_NAI_K3x3yz_S_b;
  Double I_NAI_K3x2y2z_Px_b = I_NAI_L4x2y2z_S_b+ABX*I_NAI_K3x2y2z_S_b;
  Double I_NAI_K3xy3z_Px_b = I_NAI_L4xy3z_S_b+ABX*I_NAI_K3xy3z_S_b;
  Double I_NAI_K3x4z_Px_b = I_NAI_L4x4z_S_b+ABX*I_NAI_K3x4z_S_b;
  Double I_NAI_K2x5y_Px_b = I_NAI_L3x5y_S_b+ABX*I_NAI_K2x5y_S_b;
  Double I_NAI_K2x4yz_Px_b = I_NAI_L3x4yz_S_b+ABX*I_NAI_K2x4yz_S_b;
  Double I_NAI_K2x3y2z_Px_b = I_NAI_L3x3y2z_S_b+ABX*I_NAI_K2x3y2z_S_b;
  Double I_NAI_K2x2y3z_Px_b = I_NAI_L3x2y3z_S_b+ABX*I_NAI_K2x2y3z_S_b;
  Double I_NAI_K2xy4z_Px_b = I_NAI_L3xy4z_S_b+ABX*I_NAI_K2xy4z_S_b;
  Double I_NAI_K2x5z_Px_b = I_NAI_L3x5z_S_b+ABX*I_NAI_K2x5z_S_b;
  Double I_NAI_Kx6y_Px_b = I_NAI_L2x6y_S_b+ABX*I_NAI_Kx6y_S_b;
  Double I_NAI_Kx5yz_Px_b = I_NAI_L2x5yz_S_b+ABX*I_NAI_Kx5yz_S_b;
  Double I_NAI_Kx4y2z_Px_b = I_NAI_L2x4y2z_S_b+ABX*I_NAI_Kx4y2z_S_b;
  Double I_NAI_Kx3y3z_Px_b = I_NAI_L2x3y3z_S_b+ABX*I_NAI_Kx3y3z_S_b;
  Double I_NAI_Kx2y4z_Px_b = I_NAI_L2x2y4z_S_b+ABX*I_NAI_Kx2y4z_S_b;
  Double I_NAI_Kxy5z_Px_b = I_NAI_L2xy5z_S_b+ABX*I_NAI_Kxy5z_S_b;
  Double I_NAI_Kx6z_Px_b = I_NAI_L2x6z_S_b+ABX*I_NAI_Kx6z_S_b;
  Double I_NAI_K5x2y_Py_b = I_NAI_L5x3y_S_b+ABY*I_NAI_K5x2y_S_b;
  Double I_NAI_K5xyz_Py_b = I_NAI_L5x2yz_S_b+ABY*I_NAI_K5xyz_S_b;
  Double I_NAI_K4x3y_Py_b = I_NAI_L4x4y_S_b+ABY*I_NAI_K4x3y_S_b;
  Double I_NAI_K4x2yz_Py_b = I_NAI_L4x3yz_S_b+ABY*I_NAI_K4x2yz_S_b;
  Double I_NAI_K4xy2z_Py_b = I_NAI_L4x2y2z_S_b+ABY*I_NAI_K4xy2z_S_b;
  Double I_NAI_K3x4y_Py_b = I_NAI_L3x5y_S_b+ABY*I_NAI_K3x4y_S_b;
  Double I_NAI_K3x3yz_Py_b = I_NAI_L3x4yz_S_b+ABY*I_NAI_K3x3yz_S_b;
  Double I_NAI_K3x2y2z_Py_b = I_NAI_L3x3y2z_S_b+ABY*I_NAI_K3x2y2z_S_b;
  Double I_NAI_K3xy3z_Py_b = I_NAI_L3x2y3z_S_b+ABY*I_NAI_K3xy3z_S_b;
  Double I_NAI_K2x5y_Py_b = I_NAI_L2x6y_S_b+ABY*I_NAI_K2x5y_S_b;
  Double I_NAI_K2x4yz_Py_b = I_NAI_L2x5yz_S_b+ABY*I_NAI_K2x4yz_S_b;
  Double I_NAI_K2x3y2z_Py_b = I_NAI_L2x4y2z_S_b+ABY*I_NAI_K2x3y2z_S_b;
  Double I_NAI_K2x2y3z_Py_b = I_NAI_L2x3y3z_S_b+ABY*I_NAI_K2x2y3z_S_b;
  Double I_NAI_K2xy4z_Py_b = I_NAI_L2x2y4z_S_b+ABY*I_NAI_K2xy4z_S_b;
  Double I_NAI_Kx6y_Py_b = I_NAI_Lx7y_S_b+ABY*I_NAI_Kx6y_S_b;
  Double I_NAI_Kx5yz_Py_b = I_NAI_Lx6yz_S_b+ABY*I_NAI_Kx5yz_S_b;
  Double I_NAI_Kx4y2z_Py_b = I_NAI_Lx5y2z_S_b+ABY*I_NAI_Kx4y2z_S_b;
  Double I_NAI_Kx3y3z_Py_b = I_NAI_Lx4y3z_S_b+ABY*I_NAI_Kx3y3z_S_b;
  Double I_NAI_Kx2y4z_Py_b = I_NAI_Lx3y4z_S_b+ABY*I_NAI_Kx2y4z_S_b;
  Double I_NAI_Kxy5z_Py_b = I_NAI_Lx2y5z_S_b+ABY*I_NAI_Kxy5z_S_b;
  Double I_NAI_K7y_Py_b = I_NAI_L8y_S_b+ABY*I_NAI_K7y_S_b;
  Double I_NAI_K6yz_Py_b = I_NAI_L7yz_S_b+ABY*I_NAI_K6yz_S_b;
  Double I_NAI_K5y2z_Py_b = I_NAI_L6y2z_S_b+ABY*I_NAI_K5y2z_S_b;
  Double I_NAI_K4y3z_Py_b = I_NAI_L5y3z_S_b+ABY*I_NAI_K4y3z_S_b;
  Double I_NAI_K3y4z_Py_b = I_NAI_L4y4z_S_b+ABY*I_NAI_K3y4z_S_b;
  Double I_NAI_K2y5z_Py_b = I_NAI_L3y5z_S_b+ABY*I_NAI_K2y5z_S_b;
  Double I_NAI_Ky6z_Py_b = I_NAI_L2y6z_S_b+ABY*I_NAI_Ky6z_S_b;
  Double I_NAI_K5xyz_Pz_b = I_NAI_L5xy2z_S_b+ABZ*I_NAI_K5xyz_S_b;
  Double I_NAI_K5x2z_Pz_b = I_NAI_L5x3z_S_b+ABZ*I_NAI_K5x2z_S_b;
  Double I_NAI_K4x2yz_Pz_b = I_NAI_L4x2y2z_S_b+ABZ*I_NAI_K4x2yz_S_b;
  Double I_NAI_K4xy2z_Pz_b = I_NAI_L4xy3z_S_b+ABZ*I_NAI_K4xy2z_S_b;
  Double I_NAI_K4x3z_Pz_b = I_NAI_L4x4z_S_b+ABZ*I_NAI_K4x3z_S_b;
  Double I_NAI_K3x3yz_Pz_b = I_NAI_L3x3y2z_S_b+ABZ*I_NAI_K3x3yz_S_b;
  Double I_NAI_K3x2y2z_Pz_b = I_NAI_L3x2y3z_S_b+ABZ*I_NAI_K3x2y2z_S_b;
  Double I_NAI_K3xy3z_Pz_b = I_NAI_L3xy4z_S_b+ABZ*I_NAI_K3xy3z_S_b;
  Double I_NAI_K3x4z_Pz_b = I_NAI_L3x5z_S_b+ABZ*I_NAI_K3x4z_S_b;
  Double I_NAI_K2x4yz_Pz_b = I_NAI_L2x4y2z_S_b+ABZ*I_NAI_K2x4yz_S_b;
  Double I_NAI_K2x3y2z_Pz_b = I_NAI_L2x3y3z_S_b+ABZ*I_NAI_K2x3y2z_S_b;
  Double I_NAI_K2x2y3z_Pz_b = I_NAI_L2x2y4z_S_b+ABZ*I_NAI_K2x2y3z_S_b;
  Double I_NAI_K2xy4z_Pz_b = I_NAI_L2xy5z_S_b+ABZ*I_NAI_K2xy4z_S_b;
  Double I_NAI_K2x5z_Pz_b = I_NAI_L2x6z_S_b+ABZ*I_NAI_K2x5z_S_b;
  Double I_NAI_Kx5yz_Pz_b = I_NAI_Lx5y2z_S_b+ABZ*I_NAI_Kx5yz_S_b;
  Double I_NAI_Kx4y2z_Pz_b = I_NAI_Lx4y3z_S_b+ABZ*I_NAI_Kx4y2z_S_b;
  Double I_NAI_Kx3y3z_Pz_b = I_NAI_Lx3y4z_S_b+ABZ*I_NAI_Kx3y3z_S_b;
  Double I_NAI_Kx2y4z_Pz_b = I_NAI_Lx2y5z_S_b+ABZ*I_NAI_Kx2y4z_S_b;
  Double I_NAI_Kxy5z_Pz_b = I_NAI_Lxy6z_S_b+ABZ*I_NAI_Kxy5z_S_b;
  Double I_NAI_Kx6z_Pz_b = I_NAI_Lx7z_S_b+ABZ*I_NAI_Kx6z_S_b;
  Double I_NAI_K5y2z_Pz_b = I_NAI_L5y3z_S_b+ABZ*I_NAI_K5y2z_S_b;
  Double I_NAI_K4y3z_Pz_b = I_NAI_L4y4z_S_b+ABZ*I_NAI_K4y3z_S_b;
  Double I_NAI_K3y4z_Pz_b = I_NAI_L3y5z_S_b+ABZ*I_NAI_K3y4z_S_b;
  Double I_NAI_K2y5z_Pz_b = I_NAI_L2y6z_S_b+ABZ*I_NAI_K2y5z_S_b;
  Double I_NAI_Ky6z_Pz_b = I_NAI_Ly7z_S_b+ABZ*I_NAI_Ky6z_S_b;
  Double I_NAI_K7z_Pz_b = I_NAI_L8z_S_b+ABZ*I_NAI_K7z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 87 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P_b
   * RHS shell quartet name: SQ_NAI_I_P_b
   ************************************************************/
  Double I_NAI_I6x_D2x_b = I_NAI_K7x_Px_b+ABX*I_NAI_I6x_Px_b;
  Double I_NAI_I5xy_D2x_b = I_NAI_K6xy_Px_b+ABX*I_NAI_I5xy_Px_b;
  Double I_NAI_I5xz_D2x_b = I_NAI_K6xz_Px_b+ABX*I_NAI_I5xz_Px_b;
  Double I_NAI_I4x2y_D2x_b = I_NAI_K5x2y_Px_b+ABX*I_NAI_I4x2y_Px_b;
  Double I_NAI_I4xyz_D2x_b = I_NAI_K5xyz_Px_b+ABX*I_NAI_I4xyz_Px_b;
  Double I_NAI_I4x2z_D2x_b = I_NAI_K5x2z_Px_b+ABX*I_NAI_I4x2z_Px_b;
  Double I_NAI_I3x3y_D2x_b = I_NAI_K4x3y_Px_b+ABX*I_NAI_I3x3y_Px_b;
  Double I_NAI_I3x2yz_D2x_b = I_NAI_K4x2yz_Px_b+ABX*I_NAI_I3x2yz_Px_b;
  Double I_NAI_I3xy2z_D2x_b = I_NAI_K4xy2z_Px_b+ABX*I_NAI_I3xy2z_Px_b;
  Double I_NAI_I3x3z_D2x_b = I_NAI_K4x3z_Px_b+ABX*I_NAI_I3x3z_Px_b;
  Double I_NAI_I2x4y_D2x_b = I_NAI_K3x4y_Px_b+ABX*I_NAI_I2x4y_Px_b;
  Double I_NAI_I2x3yz_D2x_b = I_NAI_K3x3yz_Px_b+ABX*I_NAI_I2x3yz_Px_b;
  Double I_NAI_I2x2y2z_D2x_b = I_NAI_K3x2y2z_Px_b+ABX*I_NAI_I2x2y2z_Px_b;
  Double I_NAI_I2xy3z_D2x_b = I_NAI_K3xy3z_Px_b+ABX*I_NAI_I2xy3z_Px_b;
  Double I_NAI_I2x4z_D2x_b = I_NAI_K3x4z_Px_b+ABX*I_NAI_I2x4z_Px_b;
  Double I_NAI_Ix5y_D2x_b = I_NAI_K2x5y_Px_b+ABX*I_NAI_Ix5y_Px_b;
  Double I_NAI_Ix4yz_D2x_b = I_NAI_K2x4yz_Px_b+ABX*I_NAI_Ix4yz_Px_b;
  Double I_NAI_Ix3y2z_D2x_b = I_NAI_K2x3y2z_Px_b+ABX*I_NAI_Ix3y2z_Px_b;
  Double I_NAI_Ix2y3z_D2x_b = I_NAI_K2x2y3z_Px_b+ABX*I_NAI_Ix2y3z_Px_b;
  Double I_NAI_Ixy4z_D2x_b = I_NAI_K2xy4z_Px_b+ABX*I_NAI_Ixy4z_Px_b;
  Double I_NAI_Ix5z_D2x_b = I_NAI_K2x5z_Px_b+ABX*I_NAI_Ix5z_Px_b;
  Double I_NAI_I6y_D2x_b = I_NAI_Kx6y_Px_b+ABX*I_NAI_I6y_Px_b;
  Double I_NAI_I5yz_D2x_b = I_NAI_Kx5yz_Px_b+ABX*I_NAI_I5yz_Px_b;
  Double I_NAI_I4y2z_D2x_b = I_NAI_Kx4y2z_Px_b+ABX*I_NAI_I4y2z_Px_b;
  Double I_NAI_I3y3z_D2x_b = I_NAI_Kx3y3z_Px_b+ABX*I_NAI_I3y3z_Px_b;
  Double I_NAI_I2y4z_D2x_b = I_NAI_Kx2y4z_Px_b+ABX*I_NAI_I2y4z_Px_b;
  Double I_NAI_Iy5z_D2x_b = I_NAI_Kxy5z_Px_b+ABX*I_NAI_Iy5z_Px_b;
  Double I_NAI_I6z_D2x_b = I_NAI_Kx6z_Px_b+ABX*I_NAI_I6z_Px_b;
  Double I_NAI_I5xy_D2y_b = I_NAI_K5x2y_Py_b+ABY*I_NAI_I5xy_Py_b;
  Double I_NAI_I5xz_D2y_b = I_NAI_K5xyz_Py_b+ABY*I_NAI_I5xz_Py_b;
  Double I_NAI_I4x2y_D2y_b = I_NAI_K4x3y_Py_b+ABY*I_NAI_I4x2y_Py_b;
  Double I_NAI_I4xyz_D2y_b = I_NAI_K4x2yz_Py_b+ABY*I_NAI_I4xyz_Py_b;
  Double I_NAI_I4x2z_D2y_b = I_NAI_K4xy2z_Py_b+ABY*I_NAI_I4x2z_Py_b;
  Double I_NAI_I3x3y_D2y_b = I_NAI_K3x4y_Py_b+ABY*I_NAI_I3x3y_Py_b;
  Double I_NAI_I3x2yz_D2y_b = I_NAI_K3x3yz_Py_b+ABY*I_NAI_I3x2yz_Py_b;
  Double I_NAI_I3xy2z_D2y_b = I_NAI_K3x2y2z_Py_b+ABY*I_NAI_I3xy2z_Py_b;
  Double I_NAI_I3x3z_D2y_b = I_NAI_K3xy3z_Py_b+ABY*I_NAI_I3x3z_Py_b;
  Double I_NAI_I2x4y_D2y_b = I_NAI_K2x5y_Py_b+ABY*I_NAI_I2x4y_Py_b;
  Double I_NAI_I2x3yz_D2y_b = I_NAI_K2x4yz_Py_b+ABY*I_NAI_I2x3yz_Py_b;
  Double I_NAI_I2x2y2z_D2y_b = I_NAI_K2x3y2z_Py_b+ABY*I_NAI_I2x2y2z_Py_b;
  Double I_NAI_I2xy3z_D2y_b = I_NAI_K2x2y3z_Py_b+ABY*I_NAI_I2xy3z_Py_b;
  Double I_NAI_I2x4z_D2y_b = I_NAI_K2xy4z_Py_b+ABY*I_NAI_I2x4z_Py_b;
  Double I_NAI_Ix5y_D2y_b = I_NAI_Kx6y_Py_b+ABY*I_NAI_Ix5y_Py_b;
  Double I_NAI_Ix4yz_D2y_b = I_NAI_Kx5yz_Py_b+ABY*I_NAI_Ix4yz_Py_b;
  Double I_NAI_Ix3y2z_D2y_b = I_NAI_Kx4y2z_Py_b+ABY*I_NAI_Ix3y2z_Py_b;
  Double I_NAI_Ix2y3z_D2y_b = I_NAI_Kx3y3z_Py_b+ABY*I_NAI_Ix2y3z_Py_b;
  Double I_NAI_Ixy4z_D2y_b = I_NAI_Kx2y4z_Py_b+ABY*I_NAI_Ixy4z_Py_b;
  Double I_NAI_Ix5z_D2y_b = I_NAI_Kxy5z_Py_b+ABY*I_NAI_Ix5z_Py_b;
  Double I_NAI_I6y_D2y_b = I_NAI_K7y_Py_b+ABY*I_NAI_I6y_Py_b;
  Double I_NAI_I5yz_D2y_b = I_NAI_K6yz_Py_b+ABY*I_NAI_I5yz_Py_b;
  Double I_NAI_I4y2z_D2y_b = I_NAI_K5y2z_Py_b+ABY*I_NAI_I4y2z_Py_b;
  Double I_NAI_I3y3z_D2y_b = I_NAI_K4y3z_Py_b+ABY*I_NAI_I3y3z_Py_b;
  Double I_NAI_I2y4z_D2y_b = I_NAI_K3y4z_Py_b+ABY*I_NAI_I2y4z_Py_b;
  Double I_NAI_Iy5z_D2y_b = I_NAI_K2y5z_Py_b+ABY*I_NAI_Iy5z_Py_b;
  Double I_NAI_I6z_D2y_b = I_NAI_Ky6z_Py_b+ABY*I_NAI_I6z_Py_b;
  Double I_NAI_I5xy_D2z_b = I_NAI_K5xyz_Pz_b+ABZ*I_NAI_I5xy_Pz_b;
  Double I_NAI_I5xz_D2z_b = I_NAI_K5x2z_Pz_b+ABZ*I_NAI_I5xz_Pz_b;
  Double I_NAI_I4x2y_D2z_b = I_NAI_K4x2yz_Pz_b+ABZ*I_NAI_I4x2y_Pz_b;
  Double I_NAI_I4xyz_D2z_b = I_NAI_K4xy2z_Pz_b+ABZ*I_NAI_I4xyz_Pz_b;
  Double I_NAI_I4x2z_D2z_b = I_NAI_K4x3z_Pz_b+ABZ*I_NAI_I4x2z_Pz_b;
  Double I_NAI_I3x3y_D2z_b = I_NAI_K3x3yz_Pz_b+ABZ*I_NAI_I3x3y_Pz_b;
  Double I_NAI_I3x2yz_D2z_b = I_NAI_K3x2y2z_Pz_b+ABZ*I_NAI_I3x2yz_Pz_b;
  Double I_NAI_I3xy2z_D2z_b = I_NAI_K3xy3z_Pz_b+ABZ*I_NAI_I3xy2z_Pz_b;
  Double I_NAI_I3x3z_D2z_b = I_NAI_K3x4z_Pz_b+ABZ*I_NAI_I3x3z_Pz_b;
  Double I_NAI_I2x4y_D2z_b = I_NAI_K2x4yz_Pz_b+ABZ*I_NAI_I2x4y_Pz_b;
  Double I_NAI_I2x3yz_D2z_b = I_NAI_K2x3y2z_Pz_b+ABZ*I_NAI_I2x3yz_Pz_b;
  Double I_NAI_I2x2y2z_D2z_b = I_NAI_K2x2y3z_Pz_b+ABZ*I_NAI_I2x2y2z_Pz_b;
  Double I_NAI_I2xy3z_D2z_b = I_NAI_K2xy4z_Pz_b+ABZ*I_NAI_I2xy3z_Pz_b;
  Double I_NAI_I2x4z_D2z_b = I_NAI_K2x5z_Pz_b+ABZ*I_NAI_I2x4z_Pz_b;
  Double I_NAI_Ix5y_D2z_b = I_NAI_Kx5yz_Pz_b+ABZ*I_NAI_Ix5y_Pz_b;
  Double I_NAI_Ix4yz_D2z_b = I_NAI_Kx4y2z_Pz_b+ABZ*I_NAI_Ix4yz_Pz_b;
  Double I_NAI_Ix3y2z_D2z_b = I_NAI_Kx3y3z_Pz_b+ABZ*I_NAI_Ix3y2z_Pz_b;
  Double I_NAI_Ix2y3z_D2z_b = I_NAI_Kx2y4z_Pz_b+ABZ*I_NAI_Ix2y3z_Pz_b;
  Double I_NAI_Ixy4z_D2z_b = I_NAI_Kxy5z_Pz_b+ABZ*I_NAI_Ixy4z_Pz_b;
  Double I_NAI_Ix5z_D2z_b = I_NAI_Kx6z_Pz_b+ABZ*I_NAI_Ix5z_Pz_b;
  Double I_NAI_I5yz_D2z_b = I_NAI_K5y2z_Pz_b+ABZ*I_NAI_I5yz_Pz_b;
  Double I_NAI_I4y2z_D2z_b = I_NAI_K4y3z_Pz_b+ABZ*I_NAI_I4y2z_Pz_b;
  Double I_NAI_I3y3z_D2z_b = I_NAI_K3y4z_Pz_b+ABZ*I_NAI_I3y3z_Pz_b;
  Double I_NAI_I2y4z_D2z_b = I_NAI_K2y5z_Pz_b+ABZ*I_NAI_I2y4z_Pz_b;
  Double I_NAI_Iy5z_D2z_b = I_NAI_Ky6z_Pz_b+ABZ*I_NAI_Iy5z_Pz_b;
  Double I_NAI_I6z_D2z_b = I_NAI_K7z_Pz_b+ABZ*I_NAI_I6z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 67 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D_b
   * RHS shell quartet name: SQ_NAI_H_D_b
   ************************************************************/
  Double I_NAI_H5x_F3x_b = I_NAI_I6x_D2x_b+ABX*I_NAI_H5x_D2x_b;
  Double I_NAI_H4xy_F3x_b = I_NAI_I5xy_D2x_b+ABX*I_NAI_H4xy_D2x_b;
  Double I_NAI_H4xz_F3x_b = I_NAI_I5xz_D2x_b+ABX*I_NAI_H4xz_D2x_b;
  Double I_NAI_H3x2y_F3x_b = I_NAI_I4x2y_D2x_b+ABX*I_NAI_H3x2y_D2x_b;
  Double I_NAI_H3xyz_F3x_b = I_NAI_I4xyz_D2x_b+ABX*I_NAI_H3xyz_D2x_b;
  Double I_NAI_H3x2z_F3x_b = I_NAI_I4x2z_D2x_b+ABX*I_NAI_H3x2z_D2x_b;
  Double I_NAI_H2x3y_F3x_b = I_NAI_I3x3y_D2x_b+ABX*I_NAI_H2x3y_D2x_b;
  Double I_NAI_H2x2yz_F3x_b = I_NAI_I3x2yz_D2x_b+ABX*I_NAI_H2x2yz_D2x_b;
  Double I_NAI_H2xy2z_F3x_b = I_NAI_I3xy2z_D2x_b+ABX*I_NAI_H2xy2z_D2x_b;
  Double I_NAI_H2x3z_F3x_b = I_NAI_I3x3z_D2x_b+ABX*I_NAI_H2x3z_D2x_b;
  Double I_NAI_Hx4y_F3x_b = I_NAI_I2x4y_D2x_b+ABX*I_NAI_Hx4y_D2x_b;
  Double I_NAI_Hx3yz_F3x_b = I_NAI_I2x3yz_D2x_b+ABX*I_NAI_Hx3yz_D2x_b;
  Double I_NAI_Hx2y2z_F3x_b = I_NAI_I2x2y2z_D2x_b+ABX*I_NAI_Hx2y2z_D2x_b;
  Double I_NAI_Hxy3z_F3x_b = I_NAI_I2xy3z_D2x_b+ABX*I_NAI_Hxy3z_D2x_b;
  Double I_NAI_Hx4z_F3x_b = I_NAI_I2x4z_D2x_b+ABX*I_NAI_Hx4z_D2x_b;
  Double I_NAI_H5y_F3x_b = I_NAI_Ix5y_D2x_b+ABX*I_NAI_H5y_D2x_b;
  Double I_NAI_H4yz_F3x_b = I_NAI_Ix4yz_D2x_b+ABX*I_NAI_H4yz_D2x_b;
  Double I_NAI_H3y2z_F3x_b = I_NAI_Ix3y2z_D2x_b+ABX*I_NAI_H3y2z_D2x_b;
  Double I_NAI_H2y3z_F3x_b = I_NAI_Ix2y3z_D2x_b+ABX*I_NAI_H2y3z_D2x_b;
  Double I_NAI_Hy4z_F3x_b = I_NAI_Ixy4z_D2x_b+ABX*I_NAI_Hy4z_D2x_b;
  Double I_NAI_H5z_F3x_b = I_NAI_Ix5z_D2x_b+ABX*I_NAI_H5z_D2x_b;
  Double I_NAI_H4xy_F2xy_b = I_NAI_I4x2y_D2x_b+ABY*I_NAI_H4xy_D2x_b;
  Double I_NAI_H4xz_F2xy_b = I_NAI_I4xyz_D2x_b+ABY*I_NAI_H4xz_D2x_b;
  Double I_NAI_H3x2y_F2xy_b = I_NAI_I3x3y_D2x_b+ABY*I_NAI_H3x2y_D2x_b;
  Double I_NAI_H3xyz_F2xy_b = I_NAI_I3x2yz_D2x_b+ABY*I_NAI_H3xyz_D2x_b;
  Double I_NAI_H3x2z_F2xy_b = I_NAI_I3xy2z_D2x_b+ABY*I_NAI_H3x2z_D2x_b;
  Double I_NAI_H2x3y_F2xy_b = I_NAI_I2x4y_D2x_b+ABY*I_NAI_H2x3y_D2x_b;
  Double I_NAI_H2x2yz_F2xy_b = I_NAI_I2x3yz_D2x_b+ABY*I_NAI_H2x2yz_D2x_b;
  Double I_NAI_H2xy2z_F2xy_b = I_NAI_I2x2y2z_D2x_b+ABY*I_NAI_H2xy2z_D2x_b;
  Double I_NAI_H2x3z_F2xy_b = I_NAI_I2xy3z_D2x_b+ABY*I_NAI_H2x3z_D2x_b;
  Double I_NAI_Hx4y_F2xy_b = I_NAI_Ix5y_D2x_b+ABY*I_NAI_Hx4y_D2x_b;
  Double I_NAI_Hx3yz_F2xy_b = I_NAI_Ix4yz_D2x_b+ABY*I_NAI_Hx3yz_D2x_b;
  Double I_NAI_Hx2y2z_F2xy_b = I_NAI_Ix3y2z_D2x_b+ABY*I_NAI_Hx2y2z_D2x_b;
  Double I_NAI_Hxy3z_F2xy_b = I_NAI_Ix2y3z_D2x_b+ABY*I_NAI_Hxy3z_D2x_b;
  Double I_NAI_Hx4z_F2xy_b = I_NAI_Ixy4z_D2x_b+ABY*I_NAI_Hx4z_D2x_b;
  Double I_NAI_H5y_F2xy_b = I_NAI_I6y_D2x_b+ABY*I_NAI_H5y_D2x_b;
  Double I_NAI_H4yz_F2xy_b = I_NAI_I5yz_D2x_b+ABY*I_NAI_H4yz_D2x_b;
  Double I_NAI_H3y2z_F2xy_b = I_NAI_I4y2z_D2x_b+ABY*I_NAI_H3y2z_D2x_b;
  Double I_NAI_H2y3z_F2xy_b = I_NAI_I3y3z_D2x_b+ABY*I_NAI_H2y3z_D2x_b;
  Double I_NAI_Hy4z_F2xy_b = I_NAI_I2y4z_D2x_b+ABY*I_NAI_Hy4z_D2x_b;
  Double I_NAI_H5z_F2xy_b = I_NAI_Iy5z_D2x_b+ABY*I_NAI_H5z_D2x_b;
  Double I_NAI_H4xz_F2xz_b = I_NAI_I4x2z_D2x_b+ABZ*I_NAI_H4xz_D2x_b;
  Double I_NAI_H3xyz_F2xz_b = I_NAI_I3xy2z_D2x_b+ABZ*I_NAI_H3xyz_D2x_b;
  Double I_NAI_H3x2z_F2xz_b = I_NAI_I3x3z_D2x_b+ABZ*I_NAI_H3x2z_D2x_b;
  Double I_NAI_H2x2yz_F2xz_b = I_NAI_I2x2y2z_D2x_b+ABZ*I_NAI_H2x2yz_D2x_b;
  Double I_NAI_H2xy2z_F2xz_b = I_NAI_I2xy3z_D2x_b+ABZ*I_NAI_H2xy2z_D2x_b;
  Double I_NAI_H2x3z_F2xz_b = I_NAI_I2x4z_D2x_b+ABZ*I_NAI_H2x3z_D2x_b;
  Double I_NAI_Hx3yz_F2xz_b = I_NAI_Ix3y2z_D2x_b+ABZ*I_NAI_Hx3yz_D2x_b;
  Double I_NAI_Hx2y2z_F2xz_b = I_NAI_Ix2y3z_D2x_b+ABZ*I_NAI_Hx2y2z_D2x_b;
  Double I_NAI_Hxy3z_F2xz_b = I_NAI_Ixy4z_D2x_b+ABZ*I_NAI_Hxy3z_D2x_b;
  Double I_NAI_Hx4z_F2xz_b = I_NAI_Ix5z_D2x_b+ABZ*I_NAI_Hx4z_D2x_b;
  Double I_NAI_H4yz_F2xz_b = I_NAI_I4y2z_D2x_b+ABZ*I_NAI_H4yz_D2x_b;
  Double I_NAI_H3y2z_F2xz_b = I_NAI_I3y3z_D2x_b+ABZ*I_NAI_H3y2z_D2x_b;
  Double I_NAI_H2y3z_F2xz_b = I_NAI_I2y4z_D2x_b+ABZ*I_NAI_H2y3z_D2x_b;
  Double I_NAI_Hy4z_F2xz_b = I_NAI_Iy5z_D2x_b+ABZ*I_NAI_Hy4z_D2x_b;
  Double I_NAI_H5z_F2xz_b = I_NAI_I6z_D2x_b+ABZ*I_NAI_H5z_D2x_b;
  Double I_NAI_H4xz_Fx2y_b = I_NAI_I5xz_D2y_b+ABX*I_NAI_H4xz_D2y_b;
  Double I_NAI_H3xyz_Fx2y_b = I_NAI_I4xyz_D2y_b+ABX*I_NAI_H3xyz_D2y_b;
  Double I_NAI_H3x2z_Fx2y_b = I_NAI_I4x2z_D2y_b+ABX*I_NAI_H3x2z_D2y_b;
  Double I_NAI_H2x2yz_Fx2y_b = I_NAI_I3x2yz_D2y_b+ABX*I_NAI_H2x2yz_D2y_b;
  Double I_NAI_H2xy2z_Fx2y_b = I_NAI_I3xy2z_D2y_b+ABX*I_NAI_H2xy2z_D2y_b;
  Double I_NAI_H2x3z_Fx2y_b = I_NAI_I3x3z_D2y_b+ABX*I_NAI_H2x3z_D2y_b;
  Double I_NAI_Hx3yz_Fx2y_b = I_NAI_I2x3yz_D2y_b+ABX*I_NAI_Hx3yz_D2y_b;
  Double I_NAI_Hx2y2z_Fx2y_b = I_NAI_I2x2y2z_D2y_b+ABX*I_NAI_Hx2y2z_D2y_b;
  Double I_NAI_Hxy3z_Fx2y_b = I_NAI_I2xy3z_D2y_b+ABX*I_NAI_Hxy3z_D2y_b;
  Double I_NAI_Hx4z_Fx2y_b = I_NAI_I2x4z_D2y_b+ABX*I_NAI_Hx4z_D2y_b;
  Double I_NAI_H4yz_Fx2y_b = I_NAI_Ix4yz_D2y_b+ABX*I_NAI_H4yz_D2y_b;
  Double I_NAI_H3y2z_Fx2y_b = I_NAI_Ix3y2z_D2y_b+ABX*I_NAI_H3y2z_D2y_b;
  Double I_NAI_H2y3z_Fx2y_b = I_NAI_Ix2y3z_D2y_b+ABX*I_NAI_H2y3z_D2y_b;
  Double I_NAI_Hy4z_Fx2y_b = I_NAI_Ixy4z_D2y_b+ABX*I_NAI_Hy4z_D2y_b;
  Double I_NAI_H5z_Fx2y_b = I_NAI_Ix5z_D2y_b+ABX*I_NAI_H5z_D2y_b;
  Double I_NAI_H4xy_Fx2z_b = I_NAI_I5xy_D2z_b+ABX*I_NAI_H4xy_D2z_b;
  Double I_NAI_H3x2y_Fx2z_b = I_NAI_I4x2y_D2z_b+ABX*I_NAI_H3x2y_D2z_b;
  Double I_NAI_H3xyz_Fx2z_b = I_NAI_I4xyz_D2z_b+ABX*I_NAI_H3xyz_D2z_b;
  Double I_NAI_H2x3y_Fx2z_b = I_NAI_I3x3y_D2z_b+ABX*I_NAI_H2x3y_D2z_b;
  Double I_NAI_H2x2yz_Fx2z_b = I_NAI_I3x2yz_D2z_b+ABX*I_NAI_H2x2yz_D2z_b;
  Double I_NAI_H2xy2z_Fx2z_b = I_NAI_I3xy2z_D2z_b+ABX*I_NAI_H2xy2z_D2z_b;
  Double I_NAI_Hx4y_Fx2z_b = I_NAI_I2x4y_D2z_b+ABX*I_NAI_Hx4y_D2z_b;
  Double I_NAI_Hx3yz_Fx2z_b = I_NAI_I2x3yz_D2z_b+ABX*I_NAI_Hx3yz_D2z_b;
  Double I_NAI_Hx2y2z_Fx2z_b = I_NAI_I2x2y2z_D2z_b+ABX*I_NAI_Hx2y2z_D2z_b;
  Double I_NAI_Hxy3z_Fx2z_b = I_NAI_I2xy3z_D2z_b+ABX*I_NAI_Hxy3z_D2z_b;
  Double I_NAI_H5y_Fx2z_b = I_NAI_Ix5y_D2z_b+ABX*I_NAI_H5y_D2z_b;
  Double I_NAI_H4yz_Fx2z_b = I_NAI_Ix4yz_D2z_b+ABX*I_NAI_H4yz_D2z_b;
  Double I_NAI_H3y2z_Fx2z_b = I_NAI_Ix3y2z_D2z_b+ABX*I_NAI_H3y2z_D2z_b;
  Double I_NAI_H2y3z_Fx2z_b = I_NAI_Ix2y3z_D2z_b+ABX*I_NAI_H2y3z_D2z_b;
  Double I_NAI_Hy4z_Fx2z_b = I_NAI_Ixy4z_D2z_b+ABX*I_NAI_Hy4z_D2z_b;
  Double I_NAI_H5x_F3y_b = I_NAI_I5xy_D2y_b+ABY*I_NAI_H5x_D2y_b;
  Double I_NAI_H4xy_F3y_b = I_NAI_I4x2y_D2y_b+ABY*I_NAI_H4xy_D2y_b;
  Double I_NAI_H4xz_F3y_b = I_NAI_I4xyz_D2y_b+ABY*I_NAI_H4xz_D2y_b;
  Double I_NAI_H3x2y_F3y_b = I_NAI_I3x3y_D2y_b+ABY*I_NAI_H3x2y_D2y_b;
  Double I_NAI_H3xyz_F3y_b = I_NAI_I3x2yz_D2y_b+ABY*I_NAI_H3xyz_D2y_b;
  Double I_NAI_H3x2z_F3y_b = I_NAI_I3xy2z_D2y_b+ABY*I_NAI_H3x2z_D2y_b;
  Double I_NAI_H2x3y_F3y_b = I_NAI_I2x4y_D2y_b+ABY*I_NAI_H2x3y_D2y_b;
  Double I_NAI_H2x2yz_F3y_b = I_NAI_I2x3yz_D2y_b+ABY*I_NAI_H2x2yz_D2y_b;
  Double I_NAI_H2xy2z_F3y_b = I_NAI_I2x2y2z_D2y_b+ABY*I_NAI_H2xy2z_D2y_b;
  Double I_NAI_H2x3z_F3y_b = I_NAI_I2xy3z_D2y_b+ABY*I_NAI_H2x3z_D2y_b;
  Double I_NAI_Hx4y_F3y_b = I_NAI_Ix5y_D2y_b+ABY*I_NAI_Hx4y_D2y_b;
  Double I_NAI_Hx3yz_F3y_b = I_NAI_Ix4yz_D2y_b+ABY*I_NAI_Hx3yz_D2y_b;
  Double I_NAI_Hx2y2z_F3y_b = I_NAI_Ix3y2z_D2y_b+ABY*I_NAI_Hx2y2z_D2y_b;
  Double I_NAI_Hxy3z_F3y_b = I_NAI_Ix2y3z_D2y_b+ABY*I_NAI_Hxy3z_D2y_b;
  Double I_NAI_Hx4z_F3y_b = I_NAI_Ixy4z_D2y_b+ABY*I_NAI_Hx4z_D2y_b;
  Double I_NAI_H5y_F3y_b = I_NAI_I6y_D2y_b+ABY*I_NAI_H5y_D2y_b;
  Double I_NAI_H4yz_F3y_b = I_NAI_I5yz_D2y_b+ABY*I_NAI_H4yz_D2y_b;
  Double I_NAI_H3y2z_F3y_b = I_NAI_I4y2z_D2y_b+ABY*I_NAI_H3y2z_D2y_b;
  Double I_NAI_H2y3z_F3y_b = I_NAI_I3y3z_D2y_b+ABY*I_NAI_H2y3z_D2y_b;
  Double I_NAI_Hy4z_F3y_b = I_NAI_I2y4z_D2y_b+ABY*I_NAI_Hy4z_D2y_b;
  Double I_NAI_H5z_F3y_b = I_NAI_Iy5z_D2y_b+ABY*I_NAI_H5z_D2y_b;
  Double I_NAI_H4xz_F2yz_b = I_NAI_I4x2z_D2y_b+ABZ*I_NAI_H4xz_D2y_b;
  Double I_NAI_H3xyz_F2yz_b = I_NAI_I3xy2z_D2y_b+ABZ*I_NAI_H3xyz_D2y_b;
  Double I_NAI_H3x2z_F2yz_b = I_NAI_I3x3z_D2y_b+ABZ*I_NAI_H3x2z_D2y_b;
  Double I_NAI_H2x2yz_F2yz_b = I_NAI_I2x2y2z_D2y_b+ABZ*I_NAI_H2x2yz_D2y_b;
  Double I_NAI_H2xy2z_F2yz_b = I_NAI_I2xy3z_D2y_b+ABZ*I_NAI_H2xy2z_D2y_b;
  Double I_NAI_H2x3z_F2yz_b = I_NAI_I2x4z_D2y_b+ABZ*I_NAI_H2x3z_D2y_b;
  Double I_NAI_Hx3yz_F2yz_b = I_NAI_Ix3y2z_D2y_b+ABZ*I_NAI_Hx3yz_D2y_b;
  Double I_NAI_Hx2y2z_F2yz_b = I_NAI_Ix2y3z_D2y_b+ABZ*I_NAI_Hx2y2z_D2y_b;
  Double I_NAI_Hxy3z_F2yz_b = I_NAI_Ixy4z_D2y_b+ABZ*I_NAI_Hxy3z_D2y_b;
  Double I_NAI_Hx4z_F2yz_b = I_NAI_Ix5z_D2y_b+ABZ*I_NAI_Hx4z_D2y_b;
  Double I_NAI_H4yz_F2yz_b = I_NAI_I4y2z_D2y_b+ABZ*I_NAI_H4yz_D2y_b;
  Double I_NAI_H3y2z_F2yz_b = I_NAI_I3y3z_D2y_b+ABZ*I_NAI_H3y2z_D2y_b;
  Double I_NAI_H2y3z_F2yz_b = I_NAI_I2y4z_D2y_b+ABZ*I_NAI_H2y3z_D2y_b;
  Double I_NAI_Hy4z_F2yz_b = I_NAI_Iy5z_D2y_b+ABZ*I_NAI_Hy4z_D2y_b;
  Double I_NAI_H5z_F2yz_b = I_NAI_I6z_D2y_b+ABZ*I_NAI_H5z_D2y_b;
  Double I_NAI_H5x_F3z_b = I_NAI_I5xz_D2z_b+ABZ*I_NAI_H5x_D2z_b;
  Double I_NAI_H4xy_F3z_b = I_NAI_I4xyz_D2z_b+ABZ*I_NAI_H4xy_D2z_b;
  Double I_NAI_H4xz_F3z_b = I_NAI_I4x2z_D2z_b+ABZ*I_NAI_H4xz_D2z_b;
  Double I_NAI_H3x2y_F3z_b = I_NAI_I3x2yz_D2z_b+ABZ*I_NAI_H3x2y_D2z_b;
  Double I_NAI_H3xyz_F3z_b = I_NAI_I3xy2z_D2z_b+ABZ*I_NAI_H3xyz_D2z_b;
  Double I_NAI_H3x2z_F3z_b = I_NAI_I3x3z_D2z_b+ABZ*I_NAI_H3x2z_D2z_b;
  Double I_NAI_H2x3y_F3z_b = I_NAI_I2x3yz_D2z_b+ABZ*I_NAI_H2x3y_D2z_b;
  Double I_NAI_H2x2yz_F3z_b = I_NAI_I2x2y2z_D2z_b+ABZ*I_NAI_H2x2yz_D2z_b;
  Double I_NAI_H2xy2z_F3z_b = I_NAI_I2xy3z_D2z_b+ABZ*I_NAI_H2xy2z_D2z_b;
  Double I_NAI_H2x3z_F3z_b = I_NAI_I2x4z_D2z_b+ABZ*I_NAI_H2x3z_D2z_b;
  Double I_NAI_Hx4y_F3z_b = I_NAI_Ix4yz_D2z_b+ABZ*I_NAI_Hx4y_D2z_b;
  Double I_NAI_Hx3yz_F3z_b = I_NAI_Ix3y2z_D2z_b+ABZ*I_NAI_Hx3yz_D2z_b;
  Double I_NAI_Hx2y2z_F3z_b = I_NAI_Ix2y3z_D2z_b+ABZ*I_NAI_Hx2y2z_D2z_b;
  Double I_NAI_Hxy3z_F3z_b = I_NAI_Ixy4z_D2z_b+ABZ*I_NAI_Hxy3z_D2z_b;
  Double I_NAI_Hx4z_F3z_b = I_NAI_Ix5z_D2z_b+ABZ*I_NAI_Hx4z_D2z_b;
  Double I_NAI_H5y_F3z_b = I_NAI_I5yz_D2z_b+ABZ*I_NAI_H5y_D2z_b;
  Double I_NAI_H4yz_F3z_b = I_NAI_I4y2z_D2z_b+ABZ*I_NAI_H4yz_D2z_b;
  Double I_NAI_H3y2z_F3z_b = I_NAI_I3y3z_D2z_b+ABZ*I_NAI_H3y2z_D2z_b;
  Double I_NAI_H2y3z_F3z_b = I_NAI_I2y4z_D2z_b+ABZ*I_NAI_H2y3z_D2z_b;
  Double I_NAI_Hy4z_F3z_b = I_NAI_Iy5z_D2z_b+ABZ*I_NAI_Hy4z_D2z_b;
  Double I_NAI_H5z_F3z_b = I_NAI_I6z_D2z_b+ABZ*I_NAI_H5z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_G_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_b
   * RHS shell quartet name: SQ_NAI_G_F_b
   ************************************************************/
  Double I_NAI_G4x_G4x_b = I_NAI_H5x_F3x_b+ABX*I_NAI_G4x_F3x_b;
  Double I_NAI_G3xy_G4x_b = I_NAI_H4xy_F3x_b+ABX*I_NAI_G3xy_F3x_b;
  Double I_NAI_G3xz_G4x_b = I_NAI_H4xz_F3x_b+ABX*I_NAI_G3xz_F3x_b;
  Double I_NAI_G2x2y_G4x_b = I_NAI_H3x2y_F3x_b+ABX*I_NAI_G2x2y_F3x_b;
  Double I_NAI_G2xyz_G4x_b = I_NAI_H3xyz_F3x_b+ABX*I_NAI_G2xyz_F3x_b;
  Double I_NAI_G2x2z_G4x_b = I_NAI_H3x2z_F3x_b+ABX*I_NAI_G2x2z_F3x_b;
  Double I_NAI_Gx3y_G4x_b = I_NAI_H2x3y_F3x_b+ABX*I_NAI_Gx3y_F3x_b;
  Double I_NAI_Gx2yz_G4x_b = I_NAI_H2x2yz_F3x_b+ABX*I_NAI_Gx2yz_F3x_b;
  Double I_NAI_Gxy2z_G4x_b = I_NAI_H2xy2z_F3x_b+ABX*I_NAI_Gxy2z_F3x_b;
  Double I_NAI_Gx3z_G4x_b = I_NAI_H2x3z_F3x_b+ABX*I_NAI_Gx3z_F3x_b;
  Double I_NAI_G4y_G4x_b = I_NAI_Hx4y_F3x_b+ABX*I_NAI_G4y_F3x_b;
  Double I_NAI_G3yz_G4x_b = I_NAI_Hx3yz_F3x_b+ABX*I_NAI_G3yz_F3x_b;
  Double I_NAI_G2y2z_G4x_b = I_NAI_Hx2y2z_F3x_b+ABX*I_NAI_G2y2z_F3x_b;
  Double I_NAI_Gy3z_G4x_b = I_NAI_Hxy3z_F3x_b+ABX*I_NAI_Gy3z_F3x_b;
  Double I_NAI_G4z_G4x_b = I_NAI_Hx4z_F3x_b+ABX*I_NAI_G4z_F3x_b;
  Double I_NAI_G4x_G3xy_b = I_NAI_H4xy_F3x_b+ABY*I_NAI_G4x_F3x_b;
  Double I_NAI_G3xy_G3xy_b = I_NAI_H3x2y_F3x_b+ABY*I_NAI_G3xy_F3x_b;
  Double I_NAI_G3xz_G3xy_b = I_NAI_H3xyz_F3x_b+ABY*I_NAI_G3xz_F3x_b;
  Double I_NAI_G2x2y_G3xy_b = I_NAI_H2x3y_F3x_b+ABY*I_NAI_G2x2y_F3x_b;
  Double I_NAI_G2xyz_G3xy_b = I_NAI_H2x2yz_F3x_b+ABY*I_NAI_G2xyz_F3x_b;
  Double I_NAI_G2x2z_G3xy_b = I_NAI_H2xy2z_F3x_b+ABY*I_NAI_G2x2z_F3x_b;
  Double I_NAI_Gx3y_G3xy_b = I_NAI_Hx4y_F3x_b+ABY*I_NAI_Gx3y_F3x_b;
  Double I_NAI_Gx2yz_G3xy_b = I_NAI_Hx3yz_F3x_b+ABY*I_NAI_Gx2yz_F3x_b;
  Double I_NAI_Gxy2z_G3xy_b = I_NAI_Hx2y2z_F3x_b+ABY*I_NAI_Gxy2z_F3x_b;
  Double I_NAI_Gx3z_G3xy_b = I_NAI_Hxy3z_F3x_b+ABY*I_NAI_Gx3z_F3x_b;
  Double I_NAI_G4y_G3xy_b = I_NAI_H5y_F3x_b+ABY*I_NAI_G4y_F3x_b;
  Double I_NAI_G3yz_G3xy_b = I_NAI_H4yz_F3x_b+ABY*I_NAI_G3yz_F3x_b;
  Double I_NAI_G2y2z_G3xy_b = I_NAI_H3y2z_F3x_b+ABY*I_NAI_G2y2z_F3x_b;
  Double I_NAI_Gy3z_G3xy_b = I_NAI_H2y3z_F3x_b+ABY*I_NAI_Gy3z_F3x_b;
  Double I_NAI_G4z_G3xy_b = I_NAI_Hy4z_F3x_b+ABY*I_NAI_G4z_F3x_b;
  Double I_NAI_G4x_G3xz_b = I_NAI_H4xz_F3x_b+ABZ*I_NAI_G4x_F3x_b;
  Double I_NAI_G3xy_G3xz_b = I_NAI_H3xyz_F3x_b+ABZ*I_NAI_G3xy_F3x_b;
  Double I_NAI_G3xz_G3xz_b = I_NAI_H3x2z_F3x_b+ABZ*I_NAI_G3xz_F3x_b;
  Double I_NAI_G2x2y_G3xz_b = I_NAI_H2x2yz_F3x_b+ABZ*I_NAI_G2x2y_F3x_b;
  Double I_NAI_G2xyz_G3xz_b = I_NAI_H2xy2z_F3x_b+ABZ*I_NAI_G2xyz_F3x_b;
  Double I_NAI_G2x2z_G3xz_b = I_NAI_H2x3z_F3x_b+ABZ*I_NAI_G2x2z_F3x_b;
  Double I_NAI_Gx3y_G3xz_b = I_NAI_Hx3yz_F3x_b+ABZ*I_NAI_Gx3y_F3x_b;
  Double I_NAI_Gx2yz_G3xz_b = I_NAI_Hx2y2z_F3x_b+ABZ*I_NAI_Gx2yz_F3x_b;
  Double I_NAI_Gxy2z_G3xz_b = I_NAI_Hxy3z_F3x_b+ABZ*I_NAI_Gxy2z_F3x_b;
  Double I_NAI_Gx3z_G3xz_b = I_NAI_Hx4z_F3x_b+ABZ*I_NAI_Gx3z_F3x_b;
  Double I_NAI_G4y_G3xz_b = I_NAI_H4yz_F3x_b+ABZ*I_NAI_G4y_F3x_b;
  Double I_NAI_G3yz_G3xz_b = I_NAI_H3y2z_F3x_b+ABZ*I_NAI_G3yz_F3x_b;
  Double I_NAI_G2y2z_G3xz_b = I_NAI_H2y3z_F3x_b+ABZ*I_NAI_G2y2z_F3x_b;
  Double I_NAI_Gy3z_G3xz_b = I_NAI_Hy4z_F3x_b+ABZ*I_NAI_Gy3z_F3x_b;
  Double I_NAI_G4z_G3xz_b = I_NAI_H5z_F3x_b+ABZ*I_NAI_G4z_F3x_b;
  Double I_NAI_G4x_G2x2y_b = I_NAI_H4xy_F2xy_b+ABY*I_NAI_G4x_F2xy_b;
  Double I_NAI_G3xy_G2x2y_b = I_NAI_H3x2y_F2xy_b+ABY*I_NAI_G3xy_F2xy_b;
  Double I_NAI_G3xz_G2x2y_b = I_NAI_H3xyz_F2xy_b+ABY*I_NAI_G3xz_F2xy_b;
  Double I_NAI_G2x2y_G2x2y_b = I_NAI_H2x3y_F2xy_b+ABY*I_NAI_G2x2y_F2xy_b;
  Double I_NAI_G2xyz_G2x2y_b = I_NAI_H2x2yz_F2xy_b+ABY*I_NAI_G2xyz_F2xy_b;
  Double I_NAI_G2x2z_G2x2y_b = I_NAI_H2xy2z_F2xy_b+ABY*I_NAI_G2x2z_F2xy_b;
  Double I_NAI_Gx3y_G2x2y_b = I_NAI_Hx4y_F2xy_b+ABY*I_NAI_Gx3y_F2xy_b;
  Double I_NAI_Gx2yz_G2x2y_b = I_NAI_Hx3yz_F2xy_b+ABY*I_NAI_Gx2yz_F2xy_b;
  Double I_NAI_Gxy2z_G2x2y_b = I_NAI_Hx2y2z_F2xy_b+ABY*I_NAI_Gxy2z_F2xy_b;
  Double I_NAI_Gx3z_G2x2y_b = I_NAI_Hxy3z_F2xy_b+ABY*I_NAI_Gx3z_F2xy_b;
  Double I_NAI_G4y_G2x2y_b = I_NAI_H5y_F2xy_b+ABY*I_NAI_G4y_F2xy_b;
  Double I_NAI_G3yz_G2x2y_b = I_NAI_H4yz_F2xy_b+ABY*I_NAI_G3yz_F2xy_b;
  Double I_NAI_G2y2z_G2x2y_b = I_NAI_H3y2z_F2xy_b+ABY*I_NAI_G2y2z_F2xy_b;
  Double I_NAI_Gy3z_G2x2y_b = I_NAI_H2y3z_F2xy_b+ABY*I_NAI_Gy3z_F2xy_b;
  Double I_NAI_G4z_G2x2y_b = I_NAI_Hy4z_F2xy_b+ABY*I_NAI_G4z_F2xy_b;
  Double I_NAI_G4x_G2xyz_b = I_NAI_H4xz_F2xy_b+ABZ*I_NAI_G4x_F2xy_b;
  Double I_NAI_G3xy_G2xyz_b = I_NAI_H3xyz_F2xy_b+ABZ*I_NAI_G3xy_F2xy_b;
  Double I_NAI_G3xz_G2xyz_b = I_NAI_H3x2z_F2xy_b+ABZ*I_NAI_G3xz_F2xy_b;
  Double I_NAI_G2x2y_G2xyz_b = I_NAI_H2x2yz_F2xy_b+ABZ*I_NAI_G2x2y_F2xy_b;
  Double I_NAI_G2xyz_G2xyz_b = I_NAI_H2xy2z_F2xy_b+ABZ*I_NAI_G2xyz_F2xy_b;
  Double I_NAI_G2x2z_G2xyz_b = I_NAI_H2x3z_F2xy_b+ABZ*I_NAI_G2x2z_F2xy_b;
  Double I_NAI_Gx3y_G2xyz_b = I_NAI_Hx3yz_F2xy_b+ABZ*I_NAI_Gx3y_F2xy_b;
  Double I_NAI_Gx2yz_G2xyz_b = I_NAI_Hx2y2z_F2xy_b+ABZ*I_NAI_Gx2yz_F2xy_b;
  Double I_NAI_Gxy2z_G2xyz_b = I_NAI_Hxy3z_F2xy_b+ABZ*I_NAI_Gxy2z_F2xy_b;
  Double I_NAI_Gx3z_G2xyz_b = I_NAI_Hx4z_F2xy_b+ABZ*I_NAI_Gx3z_F2xy_b;
  Double I_NAI_G4y_G2xyz_b = I_NAI_H4yz_F2xy_b+ABZ*I_NAI_G4y_F2xy_b;
  Double I_NAI_G3yz_G2xyz_b = I_NAI_H3y2z_F2xy_b+ABZ*I_NAI_G3yz_F2xy_b;
  Double I_NAI_G2y2z_G2xyz_b = I_NAI_H2y3z_F2xy_b+ABZ*I_NAI_G2y2z_F2xy_b;
  Double I_NAI_Gy3z_G2xyz_b = I_NAI_Hy4z_F2xy_b+ABZ*I_NAI_Gy3z_F2xy_b;
  Double I_NAI_G4z_G2xyz_b = I_NAI_H5z_F2xy_b+ABZ*I_NAI_G4z_F2xy_b;
  Double I_NAI_G4x_G2x2z_b = I_NAI_H4xz_F2xz_b+ABZ*I_NAI_G4x_F2xz_b;
  Double I_NAI_G3xy_G2x2z_b = I_NAI_H3xyz_F2xz_b+ABZ*I_NAI_G3xy_F2xz_b;
  Double I_NAI_G3xz_G2x2z_b = I_NAI_H3x2z_F2xz_b+ABZ*I_NAI_G3xz_F2xz_b;
  Double I_NAI_G2x2y_G2x2z_b = I_NAI_H2x2yz_F2xz_b+ABZ*I_NAI_G2x2y_F2xz_b;
  Double I_NAI_G2xyz_G2x2z_b = I_NAI_H2xy2z_F2xz_b+ABZ*I_NAI_G2xyz_F2xz_b;
  Double I_NAI_G2x2z_G2x2z_b = I_NAI_H2x3z_F2xz_b+ABZ*I_NAI_G2x2z_F2xz_b;
  Double I_NAI_Gx3y_G2x2z_b = I_NAI_Hx3yz_F2xz_b+ABZ*I_NAI_Gx3y_F2xz_b;
  Double I_NAI_Gx2yz_G2x2z_b = I_NAI_Hx2y2z_F2xz_b+ABZ*I_NAI_Gx2yz_F2xz_b;
  Double I_NAI_Gxy2z_G2x2z_b = I_NAI_Hxy3z_F2xz_b+ABZ*I_NAI_Gxy2z_F2xz_b;
  Double I_NAI_Gx3z_G2x2z_b = I_NAI_Hx4z_F2xz_b+ABZ*I_NAI_Gx3z_F2xz_b;
  Double I_NAI_G4y_G2x2z_b = I_NAI_H4yz_F2xz_b+ABZ*I_NAI_G4y_F2xz_b;
  Double I_NAI_G3yz_G2x2z_b = I_NAI_H3y2z_F2xz_b+ABZ*I_NAI_G3yz_F2xz_b;
  Double I_NAI_G2y2z_G2x2z_b = I_NAI_H2y3z_F2xz_b+ABZ*I_NAI_G2y2z_F2xz_b;
  Double I_NAI_Gy3z_G2x2z_b = I_NAI_Hy4z_F2xz_b+ABZ*I_NAI_Gy3z_F2xz_b;
  Double I_NAI_G4z_G2x2z_b = I_NAI_H5z_F2xz_b+ABZ*I_NAI_G4z_F2xz_b;
  Double I_NAI_G4x_Gx3y_b = I_NAI_H5x_F3y_b+ABX*I_NAI_G4x_F3y_b;
  Double I_NAI_G3xy_Gx3y_b = I_NAI_H4xy_F3y_b+ABX*I_NAI_G3xy_F3y_b;
  Double I_NAI_G3xz_Gx3y_b = I_NAI_H4xz_F3y_b+ABX*I_NAI_G3xz_F3y_b;
  Double I_NAI_G2x2y_Gx3y_b = I_NAI_H3x2y_F3y_b+ABX*I_NAI_G2x2y_F3y_b;
  Double I_NAI_G2xyz_Gx3y_b = I_NAI_H3xyz_F3y_b+ABX*I_NAI_G2xyz_F3y_b;
  Double I_NAI_G2x2z_Gx3y_b = I_NAI_H3x2z_F3y_b+ABX*I_NAI_G2x2z_F3y_b;
  Double I_NAI_Gx3y_Gx3y_b = I_NAI_H2x3y_F3y_b+ABX*I_NAI_Gx3y_F3y_b;
  Double I_NAI_Gx2yz_Gx3y_b = I_NAI_H2x2yz_F3y_b+ABX*I_NAI_Gx2yz_F3y_b;
  Double I_NAI_Gxy2z_Gx3y_b = I_NAI_H2xy2z_F3y_b+ABX*I_NAI_Gxy2z_F3y_b;
  Double I_NAI_Gx3z_Gx3y_b = I_NAI_H2x3z_F3y_b+ABX*I_NAI_Gx3z_F3y_b;
  Double I_NAI_G4y_Gx3y_b = I_NAI_Hx4y_F3y_b+ABX*I_NAI_G4y_F3y_b;
  Double I_NAI_G3yz_Gx3y_b = I_NAI_Hx3yz_F3y_b+ABX*I_NAI_G3yz_F3y_b;
  Double I_NAI_G2y2z_Gx3y_b = I_NAI_Hx2y2z_F3y_b+ABX*I_NAI_G2y2z_F3y_b;
  Double I_NAI_Gy3z_Gx3y_b = I_NAI_Hxy3z_F3y_b+ABX*I_NAI_Gy3z_F3y_b;
  Double I_NAI_G4z_Gx3y_b = I_NAI_Hx4z_F3y_b+ABX*I_NAI_G4z_F3y_b;
  Double I_NAI_G4x_Gx2yz_b = I_NAI_H4xz_Fx2y_b+ABZ*I_NAI_G4x_Fx2y_b;
  Double I_NAI_G3xy_Gx2yz_b = I_NAI_H3xyz_Fx2y_b+ABZ*I_NAI_G3xy_Fx2y_b;
  Double I_NAI_G3xz_Gx2yz_b = I_NAI_H3x2z_Fx2y_b+ABZ*I_NAI_G3xz_Fx2y_b;
  Double I_NAI_G2x2y_Gx2yz_b = I_NAI_H2x2yz_Fx2y_b+ABZ*I_NAI_G2x2y_Fx2y_b;
  Double I_NAI_G2xyz_Gx2yz_b = I_NAI_H2xy2z_Fx2y_b+ABZ*I_NAI_G2xyz_Fx2y_b;
  Double I_NAI_G2x2z_Gx2yz_b = I_NAI_H2x3z_Fx2y_b+ABZ*I_NAI_G2x2z_Fx2y_b;
  Double I_NAI_Gx3y_Gx2yz_b = I_NAI_Hx3yz_Fx2y_b+ABZ*I_NAI_Gx3y_Fx2y_b;
  Double I_NAI_Gx2yz_Gx2yz_b = I_NAI_Hx2y2z_Fx2y_b+ABZ*I_NAI_Gx2yz_Fx2y_b;
  Double I_NAI_Gxy2z_Gx2yz_b = I_NAI_Hxy3z_Fx2y_b+ABZ*I_NAI_Gxy2z_Fx2y_b;
  Double I_NAI_Gx3z_Gx2yz_b = I_NAI_Hx4z_Fx2y_b+ABZ*I_NAI_Gx3z_Fx2y_b;
  Double I_NAI_G4y_Gx2yz_b = I_NAI_H4yz_Fx2y_b+ABZ*I_NAI_G4y_Fx2y_b;
  Double I_NAI_G3yz_Gx2yz_b = I_NAI_H3y2z_Fx2y_b+ABZ*I_NAI_G3yz_Fx2y_b;
  Double I_NAI_G2y2z_Gx2yz_b = I_NAI_H2y3z_Fx2y_b+ABZ*I_NAI_G2y2z_Fx2y_b;
  Double I_NAI_Gy3z_Gx2yz_b = I_NAI_Hy4z_Fx2y_b+ABZ*I_NAI_Gy3z_Fx2y_b;
  Double I_NAI_G4z_Gx2yz_b = I_NAI_H5z_Fx2y_b+ABZ*I_NAI_G4z_Fx2y_b;
  Double I_NAI_G4x_Gxy2z_b = I_NAI_H4xy_Fx2z_b+ABY*I_NAI_G4x_Fx2z_b;
  Double I_NAI_G3xy_Gxy2z_b = I_NAI_H3x2y_Fx2z_b+ABY*I_NAI_G3xy_Fx2z_b;
  Double I_NAI_G3xz_Gxy2z_b = I_NAI_H3xyz_Fx2z_b+ABY*I_NAI_G3xz_Fx2z_b;
  Double I_NAI_G2x2y_Gxy2z_b = I_NAI_H2x3y_Fx2z_b+ABY*I_NAI_G2x2y_Fx2z_b;
  Double I_NAI_G2xyz_Gxy2z_b = I_NAI_H2x2yz_Fx2z_b+ABY*I_NAI_G2xyz_Fx2z_b;
  Double I_NAI_G2x2z_Gxy2z_b = I_NAI_H2xy2z_Fx2z_b+ABY*I_NAI_G2x2z_Fx2z_b;
  Double I_NAI_Gx3y_Gxy2z_b = I_NAI_Hx4y_Fx2z_b+ABY*I_NAI_Gx3y_Fx2z_b;
  Double I_NAI_Gx2yz_Gxy2z_b = I_NAI_Hx3yz_Fx2z_b+ABY*I_NAI_Gx2yz_Fx2z_b;
  Double I_NAI_Gxy2z_Gxy2z_b = I_NAI_Hx2y2z_Fx2z_b+ABY*I_NAI_Gxy2z_Fx2z_b;
  Double I_NAI_Gx3z_Gxy2z_b = I_NAI_Hxy3z_Fx2z_b+ABY*I_NAI_Gx3z_Fx2z_b;
  Double I_NAI_G4y_Gxy2z_b = I_NAI_H5y_Fx2z_b+ABY*I_NAI_G4y_Fx2z_b;
  Double I_NAI_G3yz_Gxy2z_b = I_NAI_H4yz_Fx2z_b+ABY*I_NAI_G3yz_Fx2z_b;
  Double I_NAI_G2y2z_Gxy2z_b = I_NAI_H3y2z_Fx2z_b+ABY*I_NAI_G2y2z_Fx2z_b;
  Double I_NAI_Gy3z_Gxy2z_b = I_NAI_H2y3z_Fx2z_b+ABY*I_NAI_Gy3z_Fx2z_b;
  Double I_NAI_G4z_Gxy2z_b = I_NAI_Hy4z_Fx2z_b+ABY*I_NAI_G4z_Fx2z_b;
  Double I_NAI_G4x_Gx3z_b = I_NAI_H5x_F3z_b+ABX*I_NAI_G4x_F3z_b;
  Double I_NAI_G3xy_Gx3z_b = I_NAI_H4xy_F3z_b+ABX*I_NAI_G3xy_F3z_b;
  Double I_NAI_G3xz_Gx3z_b = I_NAI_H4xz_F3z_b+ABX*I_NAI_G3xz_F3z_b;
  Double I_NAI_G2x2y_Gx3z_b = I_NAI_H3x2y_F3z_b+ABX*I_NAI_G2x2y_F3z_b;
  Double I_NAI_G2xyz_Gx3z_b = I_NAI_H3xyz_F3z_b+ABX*I_NAI_G2xyz_F3z_b;
  Double I_NAI_G2x2z_Gx3z_b = I_NAI_H3x2z_F3z_b+ABX*I_NAI_G2x2z_F3z_b;
  Double I_NAI_Gx3y_Gx3z_b = I_NAI_H2x3y_F3z_b+ABX*I_NAI_Gx3y_F3z_b;
  Double I_NAI_Gx2yz_Gx3z_b = I_NAI_H2x2yz_F3z_b+ABX*I_NAI_Gx2yz_F3z_b;
  Double I_NAI_Gxy2z_Gx3z_b = I_NAI_H2xy2z_F3z_b+ABX*I_NAI_Gxy2z_F3z_b;
  Double I_NAI_Gx3z_Gx3z_b = I_NAI_H2x3z_F3z_b+ABX*I_NAI_Gx3z_F3z_b;
  Double I_NAI_G4y_Gx3z_b = I_NAI_Hx4y_F3z_b+ABX*I_NAI_G4y_F3z_b;
  Double I_NAI_G3yz_Gx3z_b = I_NAI_Hx3yz_F3z_b+ABX*I_NAI_G3yz_F3z_b;
  Double I_NAI_G2y2z_Gx3z_b = I_NAI_Hx2y2z_F3z_b+ABX*I_NAI_G2y2z_F3z_b;
  Double I_NAI_Gy3z_Gx3z_b = I_NAI_Hxy3z_F3z_b+ABX*I_NAI_Gy3z_F3z_b;
  Double I_NAI_G4z_Gx3z_b = I_NAI_Hx4z_F3z_b+ABX*I_NAI_G4z_F3z_b;
  Double I_NAI_G4x_G4y_b = I_NAI_H4xy_F3y_b+ABY*I_NAI_G4x_F3y_b;
  Double I_NAI_G3xy_G4y_b = I_NAI_H3x2y_F3y_b+ABY*I_NAI_G3xy_F3y_b;
  Double I_NAI_G3xz_G4y_b = I_NAI_H3xyz_F3y_b+ABY*I_NAI_G3xz_F3y_b;
  Double I_NAI_G2x2y_G4y_b = I_NAI_H2x3y_F3y_b+ABY*I_NAI_G2x2y_F3y_b;
  Double I_NAI_G2xyz_G4y_b = I_NAI_H2x2yz_F3y_b+ABY*I_NAI_G2xyz_F3y_b;
  Double I_NAI_G2x2z_G4y_b = I_NAI_H2xy2z_F3y_b+ABY*I_NAI_G2x2z_F3y_b;
  Double I_NAI_Gx3y_G4y_b = I_NAI_Hx4y_F3y_b+ABY*I_NAI_Gx3y_F3y_b;
  Double I_NAI_Gx2yz_G4y_b = I_NAI_Hx3yz_F3y_b+ABY*I_NAI_Gx2yz_F3y_b;
  Double I_NAI_Gxy2z_G4y_b = I_NAI_Hx2y2z_F3y_b+ABY*I_NAI_Gxy2z_F3y_b;
  Double I_NAI_Gx3z_G4y_b = I_NAI_Hxy3z_F3y_b+ABY*I_NAI_Gx3z_F3y_b;
  Double I_NAI_G4y_G4y_b = I_NAI_H5y_F3y_b+ABY*I_NAI_G4y_F3y_b;
  Double I_NAI_G3yz_G4y_b = I_NAI_H4yz_F3y_b+ABY*I_NAI_G3yz_F3y_b;
  Double I_NAI_G2y2z_G4y_b = I_NAI_H3y2z_F3y_b+ABY*I_NAI_G2y2z_F3y_b;
  Double I_NAI_Gy3z_G4y_b = I_NAI_H2y3z_F3y_b+ABY*I_NAI_Gy3z_F3y_b;
  Double I_NAI_G4z_G4y_b = I_NAI_Hy4z_F3y_b+ABY*I_NAI_G4z_F3y_b;
  Double I_NAI_G4x_G3yz_b = I_NAI_H4xz_F3y_b+ABZ*I_NAI_G4x_F3y_b;
  Double I_NAI_G3xy_G3yz_b = I_NAI_H3xyz_F3y_b+ABZ*I_NAI_G3xy_F3y_b;
  Double I_NAI_G3xz_G3yz_b = I_NAI_H3x2z_F3y_b+ABZ*I_NAI_G3xz_F3y_b;
  Double I_NAI_G2x2y_G3yz_b = I_NAI_H2x2yz_F3y_b+ABZ*I_NAI_G2x2y_F3y_b;
  Double I_NAI_G2xyz_G3yz_b = I_NAI_H2xy2z_F3y_b+ABZ*I_NAI_G2xyz_F3y_b;
  Double I_NAI_G2x2z_G3yz_b = I_NAI_H2x3z_F3y_b+ABZ*I_NAI_G2x2z_F3y_b;
  Double I_NAI_Gx3y_G3yz_b = I_NAI_Hx3yz_F3y_b+ABZ*I_NAI_Gx3y_F3y_b;
  Double I_NAI_Gx2yz_G3yz_b = I_NAI_Hx2y2z_F3y_b+ABZ*I_NAI_Gx2yz_F3y_b;
  Double I_NAI_Gxy2z_G3yz_b = I_NAI_Hxy3z_F3y_b+ABZ*I_NAI_Gxy2z_F3y_b;
  Double I_NAI_Gx3z_G3yz_b = I_NAI_Hx4z_F3y_b+ABZ*I_NAI_Gx3z_F3y_b;
  Double I_NAI_G4y_G3yz_b = I_NAI_H4yz_F3y_b+ABZ*I_NAI_G4y_F3y_b;
  Double I_NAI_G3yz_G3yz_b = I_NAI_H3y2z_F3y_b+ABZ*I_NAI_G3yz_F3y_b;
  Double I_NAI_G2y2z_G3yz_b = I_NAI_H2y3z_F3y_b+ABZ*I_NAI_G2y2z_F3y_b;
  Double I_NAI_Gy3z_G3yz_b = I_NAI_Hy4z_F3y_b+ABZ*I_NAI_Gy3z_F3y_b;
  Double I_NAI_G4z_G3yz_b = I_NAI_H5z_F3y_b+ABZ*I_NAI_G4z_F3y_b;
  Double I_NAI_G4x_G2y2z_b = I_NAI_H4xz_F2yz_b+ABZ*I_NAI_G4x_F2yz_b;
  Double I_NAI_G3xy_G2y2z_b = I_NAI_H3xyz_F2yz_b+ABZ*I_NAI_G3xy_F2yz_b;
  Double I_NAI_G3xz_G2y2z_b = I_NAI_H3x2z_F2yz_b+ABZ*I_NAI_G3xz_F2yz_b;
  Double I_NAI_G2x2y_G2y2z_b = I_NAI_H2x2yz_F2yz_b+ABZ*I_NAI_G2x2y_F2yz_b;
  Double I_NAI_G2xyz_G2y2z_b = I_NAI_H2xy2z_F2yz_b+ABZ*I_NAI_G2xyz_F2yz_b;
  Double I_NAI_G2x2z_G2y2z_b = I_NAI_H2x3z_F2yz_b+ABZ*I_NAI_G2x2z_F2yz_b;
  Double I_NAI_Gx3y_G2y2z_b = I_NAI_Hx3yz_F2yz_b+ABZ*I_NAI_Gx3y_F2yz_b;
  Double I_NAI_Gx2yz_G2y2z_b = I_NAI_Hx2y2z_F2yz_b+ABZ*I_NAI_Gx2yz_F2yz_b;
  Double I_NAI_Gxy2z_G2y2z_b = I_NAI_Hxy3z_F2yz_b+ABZ*I_NAI_Gxy2z_F2yz_b;
  Double I_NAI_Gx3z_G2y2z_b = I_NAI_Hx4z_F2yz_b+ABZ*I_NAI_Gx3z_F2yz_b;
  Double I_NAI_G4y_G2y2z_b = I_NAI_H4yz_F2yz_b+ABZ*I_NAI_G4y_F2yz_b;
  Double I_NAI_G3yz_G2y2z_b = I_NAI_H3y2z_F2yz_b+ABZ*I_NAI_G3yz_F2yz_b;
  Double I_NAI_G2y2z_G2y2z_b = I_NAI_H2y3z_F2yz_b+ABZ*I_NAI_G2y2z_F2yz_b;
  Double I_NAI_Gy3z_G2y2z_b = I_NAI_Hy4z_F2yz_b+ABZ*I_NAI_Gy3z_F2yz_b;
  Double I_NAI_G4z_G2y2z_b = I_NAI_H5z_F2yz_b+ABZ*I_NAI_G4z_F2yz_b;
  Double I_NAI_G4x_Gy3z_b = I_NAI_H4xy_F3z_b+ABY*I_NAI_G4x_F3z_b;
  Double I_NAI_G3xy_Gy3z_b = I_NAI_H3x2y_F3z_b+ABY*I_NAI_G3xy_F3z_b;
  Double I_NAI_G3xz_Gy3z_b = I_NAI_H3xyz_F3z_b+ABY*I_NAI_G3xz_F3z_b;
  Double I_NAI_G2x2y_Gy3z_b = I_NAI_H2x3y_F3z_b+ABY*I_NAI_G2x2y_F3z_b;
  Double I_NAI_G2xyz_Gy3z_b = I_NAI_H2x2yz_F3z_b+ABY*I_NAI_G2xyz_F3z_b;
  Double I_NAI_G2x2z_Gy3z_b = I_NAI_H2xy2z_F3z_b+ABY*I_NAI_G2x2z_F3z_b;
  Double I_NAI_Gx3y_Gy3z_b = I_NAI_Hx4y_F3z_b+ABY*I_NAI_Gx3y_F3z_b;
  Double I_NAI_Gx2yz_Gy3z_b = I_NAI_Hx3yz_F3z_b+ABY*I_NAI_Gx2yz_F3z_b;
  Double I_NAI_Gxy2z_Gy3z_b = I_NAI_Hx2y2z_F3z_b+ABY*I_NAI_Gxy2z_F3z_b;
  Double I_NAI_Gx3z_Gy3z_b = I_NAI_Hxy3z_F3z_b+ABY*I_NAI_Gx3z_F3z_b;
  Double I_NAI_G4y_Gy3z_b = I_NAI_H5y_F3z_b+ABY*I_NAI_G4y_F3z_b;
  Double I_NAI_G3yz_Gy3z_b = I_NAI_H4yz_F3z_b+ABY*I_NAI_G3yz_F3z_b;
  Double I_NAI_G2y2z_Gy3z_b = I_NAI_H3y2z_F3z_b+ABY*I_NAI_G2y2z_F3z_b;
  Double I_NAI_Gy3z_Gy3z_b = I_NAI_H2y3z_F3z_b+ABY*I_NAI_Gy3z_F3z_b;
  Double I_NAI_G4z_Gy3z_b = I_NAI_Hy4z_F3z_b+ABY*I_NAI_G4z_F3z_b;
  Double I_NAI_G4x_G4z_b = I_NAI_H4xz_F3z_b+ABZ*I_NAI_G4x_F3z_b;
  Double I_NAI_G3xy_G4z_b = I_NAI_H3xyz_F3z_b+ABZ*I_NAI_G3xy_F3z_b;
  Double I_NAI_G3xz_G4z_b = I_NAI_H3x2z_F3z_b+ABZ*I_NAI_G3xz_F3z_b;
  Double I_NAI_G2x2y_G4z_b = I_NAI_H2x2yz_F3z_b+ABZ*I_NAI_G2x2y_F3z_b;
  Double I_NAI_G2xyz_G4z_b = I_NAI_H2xy2z_F3z_b+ABZ*I_NAI_G2xyz_F3z_b;
  Double I_NAI_G2x2z_G4z_b = I_NAI_H2x3z_F3z_b+ABZ*I_NAI_G2x2z_F3z_b;
  Double I_NAI_Gx3y_G4z_b = I_NAI_Hx3yz_F3z_b+ABZ*I_NAI_Gx3y_F3z_b;
  Double I_NAI_Gx2yz_G4z_b = I_NAI_Hx2y2z_F3z_b+ABZ*I_NAI_Gx2yz_F3z_b;
  Double I_NAI_Gxy2z_G4z_b = I_NAI_Hxy3z_F3z_b+ABZ*I_NAI_Gxy2z_F3z_b;
  Double I_NAI_Gx3z_G4z_b = I_NAI_Hx4z_F3z_b+ABZ*I_NAI_Gx3z_F3z_b;
  Double I_NAI_G4y_G4z_b = I_NAI_H4yz_F3z_b+ABZ*I_NAI_G4y_F3z_b;
  Double I_NAI_G3yz_G4z_b = I_NAI_H3y2z_F3z_b+ABZ*I_NAI_G3yz_F3z_b;
  Double I_NAI_G2y2z_G4z_b = I_NAI_H2y3z_F3z_b+ABZ*I_NAI_G2y2z_F3z_b;
  Double I_NAI_Gy3z_G4z_b = I_NAI_Hy4z_F3z_b+ABZ*I_NAI_Gy3z_F3z_b;
  Double I_NAI_G4z_G4z_b = I_NAI_H5z_F3z_b+ABZ*I_NAI_G4z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_a
   * RHS shell quartet name: SQ_NAI_F_F
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_H5x_F3x_a-4*I_NAI_F3x_F3x;
  abcd[1] = 2.0E0*I_NAI_H4xy_F3x_a-3*I_NAI_F2xy_F3x;
  abcd[2] = 2.0E0*I_NAI_H4xz_F3x_a-3*I_NAI_F2xz_F3x;
  abcd[3] = 2.0E0*I_NAI_H3x2y_F3x_a-2*I_NAI_Fx2y_F3x;
  abcd[4] = 2.0E0*I_NAI_H3xyz_F3x_a-2*I_NAI_Fxyz_F3x;
  abcd[5] = 2.0E0*I_NAI_H3x2z_F3x_a-2*I_NAI_Fx2z_F3x;
  abcd[6] = 2.0E0*I_NAI_H2x3y_F3x_a-1*I_NAI_F3y_F3x;
  abcd[7] = 2.0E0*I_NAI_H2x2yz_F3x_a-1*I_NAI_F2yz_F3x;
  abcd[8] = 2.0E0*I_NAI_H2xy2z_F3x_a-1*I_NAI_Fy2z_F3x;
  abcd[9] = 2.0E0*I_NAI_H2x3z_F3x_a-1*I_NAI_F3z_F3x;
  abcd[10] = 2.0E0*I_NAI_Hx4y_F3x_a;
  abcd[11] = 2.0E0*I_NAI_Hx3yz_F3x_a;
  abcd[12] = 2.0E0*I_NAI_Hx2y2z_F3x_a;
  abcd[13] = 2.0E0*I_NAI_Hxy3z_F3x_a;
  abcd[14] = 2.0E0*I_NAI_Hx4z_F3x_a;
  abcd[15] = 2.0E0*I_NAI_H5x_F2xy_a-4*I_NAI_F3x_F2xy;
  abcd[16] = 2.0E0*I_NAI_H4xy_F2xy_a-3*I_NAI_F2xy_F2xy;
  abcd[17] = 2.0E0*I_NAI_H4xz_F2xy_a-3*I_NAI_F2xz_F2xy;
  abcd[18] = 2.0E0*I_NAI_H3x2y_F2xy_a-2*I_NAI_Fx2y_F2xy;
  abcd[19] = 2.0E0*I_NAI_H3xyz_F2xy_a-2*I_NAI_Fxyz_F2xy;
  abcd[20] = 2.0E0*I_NAI_H3x2z_F2xy_a-2*I_NAI_Fx2z_F2xy;
  abcd[21] = 2.0E0*I_NAI_H2x3y_F2xy_a-1*I_NAI_F3y_F2xy;
  abcd[22] = 2.0E0*I_NAI_H2x2yz_F2xy_a-1*I_NAI_F2yz_F2xy;
  abcd[23] = 2.0E0*I_NAI_H2xy2z_F2xy_a-1*I_NAI_Fy2z_F2xy;
  abcd[24] = 2.0E0*I_NAI_H2x3z_F2xy_a-1*I_NAI_F3z_F2xy;
  abcd[25] = 2.0E0*I_NAI_Hx4y_F2xy_a;
  abcd[26] = 2.0E0*I_NAI_Hx3yz_F2xy_a;
  abcd[27] = 2.0E0*I_NAI_Hx2y2z_F2xy_a;
  abcd[28] = 2.0E0*I_NAI_Hxy3z_F2xy_a;
  abcd[29] = 2.0E0*I_NAI_Hx4z_F2xy_a;
  abcd[30] = 2.0E0*I_NAI_H5x_F2xz_a-4*I_NAI_F3x_F2xz;
  abcd[31] = 2.0E0*I_NAI_H4xy_F2xz_a-3*I_NAI_F2xy_F2xz;
  abcd[32] = 2.0E0*I_NAI_H4xz_F2xz_a-3*I_NAI_F2xz_F2xz;
  abcd[33] = 2.0E0*I_NAI_H3x2y_F2xz_a-2*I_NAI_Fx2y_F2xz;
  abcd[34] = 2.0E0*I_NAI_H3xyz_F2xz_a-2*I_NAI_Fxyz_F2xz;
  abcd[35] = 2.0E0*I_NAI_H3x2z_F2xz_a-2*I_NAI_Fx2z_F2xz;
  abcd[36] = 2.0E0*I_NAI_H2x3y_F2xz_a-1*I_NAI_F3y_F2xz;
  abcd[37] = 2.0E0*I_NAI_H2x2yz_F2xz_a-1*I_NAI_F2yz_F2xz;
  abcd[38] = 2.0E0*I_NAI_H2xy2z_F2xz_a-1*I_NAI_Fy2z_F2xz;
  abcd[39] = 2.0E0*I_NAI_H2x3z_F2xz_a-1*I_NAI_F3z_F2xz;
  abcd[40] = 2.0E0*I_NAI_Hx4y_F2xz_a;
  abcd[41] = 2.0E0*I_NAI_Hx3yz_F2xz_a;
  abcd[42] = 2.0E0*I_NAI_Hx2y2z_F2xz_a;
  abcd[43] = 2.0E0*I_NAI_Hxy3z_F2xz_a;
  abcd[44] = 2.0E0*I_NAI_Hx4z_F2xz_a;
  abcd[45] = 2.0E0*I_NAI_H5x_Fx2y_a-4*I_NAI_F3x_Fx2y;
  abcd[46] = 2.0E0*I_NAI_H4xy_Fx2y_a-3*I_NAI_F2xy_Fx2y;
  abcd[47] = 2.0E0*I_NAI_H4xz_Fx2y_a-3*I_NAI_F2xz_Fx2y;
  abcd[48] = 2.0E0*I_NAI_H3x2y_Fx2y_a-2*I_NAI_Fx2y_Fx2y;
  abcd[49] = 2.0E0*I_NAI_H3xyz_Fx2y_a-2*I_NAI_Fxyz_Fx2y;
  abcd[50] = 2.0E0*I_NAI_H3x2z_Fx2y_a-2*I_NAI_Fx2z_Fx2y;
  abcd[51] = 2.0E0*I_NAI_H2x3y_Fx2y_a-1*I_NAI_F3y_Fx2y;
  abcd[52] = 2.0E0*I_NAI_H2x2yz_Fx2y_a-1*I_NAI_F2yz_Fx2y;
  abcd[53] = 2.0E0*I_NAI_H2xy2z_Fx2y_a-1*I_NAI_Fy2z_Fx2y;
  abcd[54] = 2.0E0*I_NAI_H2x3z_Fx2y_a-1*I_NAI_F3z_Fx2y;
  abcd[55] = 2.0E0*I_NAI_Hx4y_Fx2y_a;
  abcd[56] = 2.0E0*I_NAI_Hx3yz_Fx2y_a;
  abcd[57] = 2.0E0*I_NAI_Hx2y2z_Fx2y_a;
  abcd[58] = 2.0E0*I_NAI_Hxy3z_Fx2y_a;
  abcd[59] = 2.0E0*I_NAI_Hx4z_Fx2y_a;
  abcd[60] = 2.0E0*I_NAI_H5x_Fxyz_a-4*I_NAI_F3x_Fxyz;
  abcd[61] = 2.0E0*I_NAI_H4xy_Fxyz_a-3*I_NAI_F2xy_Fxyz;
  abcd[62] = 2.0E0*I_NAI_H4xz_Fxyz_a-3*I_NAI_F2xz_Fxyz;
  abcd[63] = 2.0E0*I_NAI_H3x2y_Fxyz_a-2*I_NAI_Fx2y_Fxyz;
  abcd[64] = 2.0E0*I_NAI_H3xyz_Fxyz_a-2*I_NAI_Fxyz_Fxyz;
  abcd[65] = 2.0E0*I_NAI_H3x2z_Fxyz_a-2*I_NAI_Fx2z_Fxyz;
  abcd[66] = 2.0E0*I_NAI_H2x3y_Fxyz_a-1*I_NAI_F3y_Fxyz;
  abcd[67] = 2.0E0*I_NAI_H2x2yz_Fxyz_a-1*I_NAI_F2yz_Fxyz;
  abcd[68] = 2.0E0*I_NAI_H2xy2z_Fxyz_a-1*I_NAI_Fy2z_Fxyz;
  abcd[69] = 2.0E0*I_NAI_H2x3z_Fxyz_a-1*I_NAI_F3z_Fxyz;
  abcd[70] = 2.0E0*I_NAI_Hx4y_Fxyz_a;
  abcd[71] = 2.0E0*I_NAI_Hx3yz_Fxyz_a;
  abcd[72] = 2.0E0*I_NAI_Hx2y2z_Fxyz_a;
  abcd[73] = 2.0E0*I_NAI_Hxy3z_Fxyz_a;
  abcd[74] = 2.0E0*I_NAI_Hx4z_Fxyz_a;
  abcd[75] = 2.0E0*I_NAI_H5x_Fx2z_a-4*I_NAI_F3x_Fx2z;
  abcd[76] = 2.0E0*I_NAI_H4xy_Fx2z_a-3*I_NAI_F2xy_Fx2z;
  abcd[77] = 2.0E0*I_NAI_H4xz_Fx2z_a-3*I_NAI_F2xz_Fx2z;
  abcd[78] = 2.0E0*I_NAI_H3x2y_Fx2z_a-2*I_NAI_Fx2y_Fx2z;
  abcd[79] = 2.0E0*I_NAI_H3xyz_Fx2z_a-2*I_NAI_Fxyz_Fx2z;
  abcd[80] = 2.0E0*I_NAI_H3x2z_Fx2z_a-2*I_NAI_Fx2z_Fx2z;
  abcd[81] = 2.0E0*I_NAI_H2x3y_Fx2z_a-1*I_NAI_F3y_Fx2z;
  abcd[82] = 2.0E0*I_NAI_H2x2yz_Fx2z_a-1*I_NAI_F2yz_Fx2z;
  abcd[83] = 2.0E0*I_NAI_H2xy2z_Fx2z_a-1*I_NAI_Fy2z_Fx2z;
  abcd[84] = 2.0E0*I_NAI_H2x3z_Fx2z_a-1*I_NAI_F3z_Fx2z;
  abcd[85] = 2.0E0*I_NAI_Hx4y_Fx2z_a;
  abcd[86] = 2.0E0*I_NAI_Hx3yz_Fx2z_a;
  abcd[87] = 2.0E0*I_NAI_Hx2y2z_Fx2z_a;
  abcd[88] = 2.0E0*I_NAI_Hxy3z_Fx2z_a;
  abcd[89] = 2.0E0*I_NAI_Hx4z_Fx2z_a;
  abcd[90] = 2.0E0*I_NAI_H5x_F3y_a-4*I_NAI_F3x_F3y;
  abcd[91] = 2.0E0*I_NAI_H4xy_F3y_a-3*I_NAI_F2xy_F3y;
  abcd[92] = 2.0E0*I_NAI_H4xz_F3y_a-3*I_NAI_F2xz_F3y;
  abcd[93] = 2.0E0*I_NAI_H3x2y_F3y_a-2*I_NAI_Fx2y_F3y;
  abcd[94] = 2.0E0*I_NAI_H3xyz_F3y_a-2*I_NAI_Fxyz_F3y;
  abcd[95] = 2.0E0*I_NAI_H3x2z_F3y_a-2*I_NAI_Fx2z_F3y;
  abcd[96] = 2.0E0*I_NAI_H2x3y_F3y_a-1*I_NAI_F3y_F3y;
  abcd[97] = 2.0E0*I_NAI_H2x2yz_F3y_a-1*I_NAI_F2yz_F3y;
  abcd[98] = 2.0E0*I_NAI_H2xy2z_F3y_a-1*I_NAI_Fy2z_F3y;
  abcd[99] = 2.0E0*I_NAI_H2x3z_F3y_a-1*I_NAI_F3z_F3y;
  abcd[100] = 2.0E0*I_NAI_Hx4y_F3y_a;
  abcd[101] = 2.0E0*I_NAI_Hx3yz_F3y_a;
  abcd[102] = 2.0E0*I_NAI_Hx2y2z_F3y_a;
  abcd[103] = 2.0E0*I_NAI_Hxy3z_F3y_a;
  abcd[104] = 2.0E0*I_NAI_Hx4z_F3y_a;
  abcd[105] = 2.0E0*I_NAI_H5x_F2yz_a-4*I_NAI_F3x_F2yz;
  abcd[106] = 2.0E0*I_NAI_H4xy_F2yz_a-3*I_NAI_F2xy_F2yz;
  abcd[107] = 2.0E0*I_NAI_H4xz_F2yz_a-3*I_NAI_F2xz_F2yz;
  abcd[108] = 2.0E0*I_NAI_H3x2y_F2yz_a-2*I_NAI_Fx2y_F2yz;
  abcd[109] = 2.0E0*I_NAI_H3xyz_F2yz_a-2*I_NAI_Fxyz_F2yz;
  abcd[110] = 2.0E0*I_NAI_H3x2z_F2yz_a-2*I_NAI_Fx2z_F2yz;
  abcd[111] = 2.0E0*I_NAI_H2x3y_F2yz_a-1*I_NAI_F3y_F2yz;
  abcd[112] = 2.0E0*I_NAI_H2x2yz_F2yz_a-1*I_NAI_F2yz_F2yz;
  abcd[113] = 2.0E0*I_NAI_H2xy2z_F2yz_a-1*I_NAI_Fy2z_F2yz;
  abcd[114] = 2.0E0*I_NAI_H2x3z_F2yz_a-1*I_NAI_F3z_F2yz;
  abcd[115] = 2.0E0*I_NAI_Hx4y_F2yz_a;
  abcd[116] = 2.0E0*I_NAI_Hx3yz_F2yz_a;
  abcd[117] = 2.0E0*I_NAI_Hx2y2z_F2yz_a;
  abcd[118] = 2.0E0*I_NAI_Hxy3z_F2yz_a;
  abcd[119] = 2.0E0*I_NAI_Hx4z_F2yz_a;
  abcd[120] = 2.0E0*I_NAI_H5x_Fy2z_a-4*I_NAI_F3x_Fy2z;
  abcd[121] = 2.0E0*I_NAI_H4xy_Fy2z_a-3*I_NAI_F2xy_Fy2z;
  abcd[122] = 2.0E0*I_NAI_H4xz_Fy2z_a-3*I_NAI_F2xz_Fy2z;
  abcd[123] = 2.0E0*I_NAI_H3x2y_Fy2z_a-2*I_NAI_Fx2y_Fy2z;
  abcd[124] = 2.0E0*I_NAI_H3xyz_Fy2z_a-2*I_NAI_Fxyz_Fy2z;
  abcd[125] = 2.0E0*I_NAI_H3x2z_Fy2z_a-2*I_NAI_Fx2z_Fy2z;
  abcd[126] = 2.0E0*I_NAI_H2x3y_Fy2z_a-1*I_NAI_F3y_Fy2z;
  abcd[127] = 2.0E0*I_NAI_H2x2yz_Fy2z_a-1*I_NAI_F2yz_Fy2z;
  abcd[128] = 2.0E0*I_NAI_H2xy2z_Fy2z_a-1*I_NAI_Fy2z_Fy2z;
  abcd[129] = 2.0E0*I_NAI_H2x3z_Fy2z_a-1*I_NAI_F3z_Fy2z;
  abcd[130] = 2.0E0*I_NAI_Hx4y_Fy2z_a;
  abcd[131] = 2.0E0*I_NAI_Hx3yz_Fy2z_a;
  abcd[132] = 2.0E0*I_NAI_Hx2y2z_Fy2z_a;
  abcd[133] = 2.0E0*I_NAI_Hxy3z_Fy2z_a;
  abcd[134] = 2.0E0*I_NAI_Hx4z_Fy2z_a;
  abcd[135] = 2.0E0*I_NAI_H5x_F3z_a-4*I_NAI_F3x_F3z;
  abcd[136] = 2.0E0*I_NAI_H4xy_F3z_a-3*I_NAI_F2xy_F3z;
  abcd[137] = 2.0E0*I_NAI_H4xz_F3z_a-3*I_NAI_F2xz_F3z;
  abcd[138] = 2.0E0*I_NAI_H3x2y_F3z_a-2*I_NAI_Fx2y_F3z;
  abcd[139] = 2.0E0*I_NAI_H3xyz_F3z_a-2*I_NAI_Fxyz_F3z;
  abcd[140] = 2.0E0*I_NAI_H3x2z_F3z_a-2*I_NAI_Fx2z_F3z;
  abcd[141] = 2.0E0*I_NAI_H2x3y_F3z_a-1*I_NAI_F3y_F3z;
  abcd[142] = 2.0E0*I_NAI_H2x2yz_F3z_a-1*I_NAI_F2yz_F3z;
  abcd[143] = 2.0E0*I_NAI_H2xy2z_F3z_a-1*I_NAI_Fy2z_F3z;
  abcd[144] = 2.0E0*I_NAI_H2x3z_F3z_a-1*I_NAI_F3z_F3z;
  abcd[145] = 2.0E0*I_NAI_Hx4y_F3z_a;
  abcd[146] = 2.0E0*I_NAI_Hx3yz_F3z_a;
  abcd[147] = 2.0E0*I_NAI_Hx2y2z_F3z_a;
  abcd[148] = 2.0E0*I_NAI_Hxy3z_F3z_a;
  abcd[149] = 2.0E0*I_NAI_Hx4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_a
   * RHS shell quartet name: SQ_NAI_F_F
   ************************************************************/
  abcd[150] = 2.0E0*I_NAI_H4xy_F3x_a;
  abcd[151] = 2.0E0*I_NAI_H3x2y_F3x_a-1*I_NAI_F3x_F3x;
  abcd[152] = 2.0E0*I_NAI_H3xyz_F3x_a;
  abcd[153] = 2.0E0*I_NAI_H2x3y_F3x_a-2*I_NAI_F2xy_F3x;
  abcd[154] = 2.0E0*I_NAI_H2x2yz_F3x_a-1*I_NAI_F2xz_F3x;
  abcd[155] = 2.0E0*I_NAI_H2xy2z_F3x_a;
  abcd[156] = 2.0E0*I_NAI_Hx4y_F3x_a-3*I_NAI_Fx2y_F3x;
  abcd[157] = 2.0E0*I_NAI_Hx3yz_F3x_a-2*I_NAI_Fxyz_F3x;
  abcd[158] = 2.0E0*I_NAI_Hx2y2z_F3x_a-1*I_NAI_Fx2z_F3x;
  abcd[159] = 2.0E0*I_NAI_Hxy3z_F3x_a;
  abcd[160] = 2.0E0*I_NAI_H5y_F3x_a-4*I_NAI_F3y_F3x;
  abcd[161] = 2.0E0*I_NAI_H4yz_F3x_a-3*I_NAI_F2yz_F3x;
  abcd[162] = 2.0E0*I_NAI_H3y2z_F3x_a-2*I_NAI_Fy2z_F3x;
  abcd[163] = 2.0E0*I_NAI_H2y3z_F3x_a-1*I_NAI_F3z_F3x;
  abcd[164] = 2.0E0*I_NAI_Hy4z_F3x_a;
  abcd[165] = 2.0E0*I_NAI_H4xy_F2xy_a;
  abcd[166] = 2.0E0*I_NAI_H3x2y_F2xy_a-1*I_NAI_F3x_F2xy;
  abcd[167] = 2.0E0*I_NAI_H3xyz_F2xy_a;
  abcd[168] = 2.0E0*I_NAI_H2x3y_F2xy_a-2*I_NAI_F2xy_F2xy;
  abcd[169] = 2.0E0*I_NAI_H2x2yz_F2xy_a-1*I_NAI_F2xz_F2xy;
  abcd[170] = 2.0E0*I_NAI_H2xy2z_F2xy_a;
  abcd[171] = 2.0E0*I_NAI_Hx4y_F2xy_a-3*I_NAI_Fx2y_F2xy;
  abcd[172] = 2.0E0*I_NAI_Hx3yz_F2xy_a-2*I_NAI_Fxyz_F2xy;
  abcd[173] = 2.0E0*I_NAI_Hx2y2z_F2xy_a-1*I_NAI_Fx2z_F2xy;
  abcd[174] = 2.0E0*I_NAI_Hxy3z_F2xy_a;
  abcd[175] = 2.0E0*I_NAI_H5y_F2xy_a-4*I_NAI_F3y_F2xy;
  abcd[176] = 2.0E0*I_NAI_H4yz_F2xy_a-3*I_NAI_F2yz_F2xy;
  abcd[177] = 2.0E0*I_NAI_H3y2z_F2xy_a-2*I_NAI_Fy2z_F2xy;
  abcd[178] = 2.0E0*I_NAI_H2y3z_F2xy_a-1*I_NAI_F3z_F2xy;
  abcd[179] = 2.0E0*I_NAI_Hy4z_F2xy_a;
  abcd[180] = 2.0E0*I_NAI_H4xy_F2xz_a;
  abcd[181] = 2.0E0*I_NAI_H3x2y_F2xz_a-1*I_NAI_F3x_F2xz;
  abcd[182] = 2.0E0*I_NAI_H3xyz_F2xz_a;
  abcd[183] = 2.0E0*I_NAI_H2x3y_F2xz_a-2*I_NAI_F2xy_F2xz;
  abcd[184] = 2.0E0*I_NAI_H2x2yz_F2xz_a-1*I_NAI_F2xz_F2xz;
  abcd[185] = 2.0E0*I_NAI_H2xy2z_F2xz_a;
  abcd[186] = 2.0E0*I_NAI_Hx4y_F2xz_a-3*I_NAI_Fx2y_F2xz;
  abcd[187] = 2.0E0*I_NAI_Hx3yz_F2xz_a-2*I_NAI_Fxyz_F2xz;
  abcd[188] = 2.0E0*I_NAI_Hx2y2z_F2xz_a-1*I_NAI_Fx2z_F2xz;
  abcd[189] = 2.0E0*I_NAI_Hxy3z_F2xz_a;
  abcd[190] = 2.0E0*I_NAI_H5y_F2xz_a-4*I_NAI_F3y_F2xz;
  abcd[191] = 2.0E0*I_NAI_H4yz_F2xz_a-3*I_NAI_F2yz_F2xz;
  abcd[192] = 2.0E0*I_NAI_H3y2z_F2xz_a-2*I_NAI_Fy2z_F2xz;
  abcd[193] = 2.0E0*I_NAI_H2y3z_F2xz_a-1*I_NAI_F3z_F2xz;
  abcd[194] = 2.0E0*I_NAI_Hy4z_F2xz_a;
  abcd[195] = 2.0E0*I_NAI_H4xy_Fx2y_a;
  abcd[196] = 2.0E0*I_NAI_H3x2y_Fx2y_a-1*I_NAI_F3x_Fx2y;
  abcd[197] = 2.0E0*I_NAI_H3xyz_Fx2y_a;
  abcd[198] = 2.0E0*I_NAI_H2x3y_Fx2y_a-2*I_NAI_F2xy_Fx2y;
  abcd[199] = 2.0E0*I_NAI_H2x2yz_Fx2y_a-1*I_NAI_F2xz_Fx2y;
  abcd[200] = 2.0E0*I_NAI_H2xy2z_Fx2y_a;
  abcd[201] = 2.0E0*I_NAI_Hx4y_Fx2y_a-3*I_NAI_Fx2y_Fx2y;
  abcd[202] = 2.0E0*I_NAI_Hx3yz_Fx2y_a-2*I_NAI_Fxyz_Fx2y;
  abcd[203] = 2.0E0*I_NAI_Hx2y2z_Fx2y_a-1*I_NAI_Fx2z_Fx2y;
  abcd[204] = 2.0E0*I_NAI_Hxy3z_Fx2y_a;
  abcd[205] = 2.0E0*I_NAI_H5y_Fx2y_a-4*I_NAI_F3y_Fx2y;
  abcd[206] = 2.0E0*I_NAI_H4yz_Fx2y_a-3*I_NAI_F2yz_Fx2y;
  abcd[207] = 2.0E0*I_NAI_H3y2z_Fx2y_a-2*I_NAI_Fy2z_Fx2y;
  abcd[208] = 2.0E0*I_NAI_H2y3z_Fx2y_a-1*I_NAI_F3z_Fx2y;
  abcd[209] = 2.0E0*I_NAI_Hy4z_Fx2y_a;
  abcd[210] = 2.0E0*I_NAI_H4xy_Fxyz_a;
  abcd[211] = 2.0E0*I_NAI_H3x2y_Fxyz_a-1*I_NAI_F3x_Fxyz;
  abcd[212] = 2.0E0*I_NAI_H3xyz_Fxyz_a;
  abcd[213] = 2.0E0*I_NAI_H2x3y_Fxyz_a-2*I_NAI_F2xy_Fxyz;
  abcd[214] = 2.0E0*I_NAI_H2x2yz_Fxyz_a-1*I_NAI_F2xz_Fxyz;
  abcd[215] = 2.0E0*I_NAI_H2xy2z_Fxyz_a;
  abcd[216] = 2.0E0*I_NAI_Hx4y_Fxyz_a-3*I_NAI_Fx2y_Fxyz;
  abcd[217] = 2.0E0*I_NAI_Hx3yz_Fxyz_a-2*I_NAI_Fxyz_Fxyz;
  abcd[218] = 2.0E0*I_NAI_Hx2y2z_Fxyz_a-1*I_NAI_Fx2z_Fxyz;
  abcd[219] = 2.0E0*I_NAI_Hxy3z_Fxyz_a;
  abcd[220] = 2.0E0*I_NAI_H5y_Fxyz_a-4*I_NAI_F3y_Fxyz;
  abcd[221] = 2.0E0*I_NAI_H4yz_Fxyz_a-3*I_NAI_F2yz_Fxyz;
  abcd[222] = 2.0E0*I_NAI_H3y2z_Fxyz_a-2*I_NAI_Fy2z_Fxyz;
  abcd[223] = 2.0E0*I_NAI_H2y3z_Fxyz_a-1*I_NAI_F3z_Fxyz;
  abcd[224] = 2.0E0*I_NAI_Hy4z_Fxyz_a;
  abcd[225] = 2.0E0*I_NAI_H4xy_Fx2z_a;
  abcd[226] = 2.0E0*I_NAI_H3x2y_Fx2z_a-1*I_NAI_F3x_Fx2z;
  abcd[227] = 2.0E0*I_NAI_H3xyz_Fx2z_a;
  abcd[228] = 2.0E0*I_NAI_H2x3y_Fx2z_a-2*I_NAI_F2xy_Fx2z;
  abcd[229] = 2.0E0*I_NAI_H2x2yz_Fx2z_a-1*I_NAI_F2xz_Fx2z;
  abcd[230] = 2.0E0*I_NAI_H2xy2z_Fx2z_a;
  abcd[231] = 2.0E0*I_NAI_Hx4y_Fx2z_a-3*I_NAI_Fx2y_Fx2z;
  abcd[232] = 2.0E0*I_NAI_Hx3yz_Fx2z_a-2*I_NAI_Fxyz_Fx2z;
  abcd[233] = 2.0E0*I_NAI_Hx2y2z_Fx2z_a-1*I_NAI_Fx2z_Fx2z;
  abcd[234] = 2.0E0*I_NAI_Hxy3z_Fx2z_a;
  abcd[235] = 2.0E0*I_NAI_H5y_Fx2z_a-4*I_NAI_F3y_Fx2z;
  abcd[236] = 2.0E0*I_NAI_H4yz_Fx2z_a-3*I_NAI_F2yz_Fx2z;
  abcd[237] = 2.0E0*I_NAI_H3y2z_Fx2z_a-2*I_NAI_Fy2z_Fx2z;
  abcd[238] = 2.0E0*I_NAI_H2y3z_Fx2z_a-1*I_NAI_F3z_Fx2z;
  abcd[239] = 2.0E0*I_NAI_Hy4z_Fx2z_a;
  abcd[240] = 2.0E0*I_NAI_H4xy_F3y_a;
  abcd[241] = 2.0E0*I_NAI_H3x2y_F3y_a-1*I_NAI_F3x_F3y;
  abcd[242] = 2.0E0*I_NAI_H3xyz_F3y_a;
  abcd[243] = 2.0E0*I_NAI_H2x3y_F3y_a-2*I_NAI_F2xy_F3y;
  abcd[244] = 2.0E0*I_NAI_H2x2yz_F3y_a-1*I_NAI_F2xz_F3y;
  abcd[245] = 2.0E0*I_NAI_H2xy2z_F3y_a;
  abcd[246] = 2.0E0*I_NAI_Hx4y_F3y_a-3*I_NAI_Fx2y_F3y;
  abcd[247] = 2.0E0*I_NAI_Hx3yz_F3y_a-2*I_NAI_Fxyz_F3y;
  abcd[248] = 2.0E0*I_NAI_Hx2y2z_F3y_a-1*I_NAI_Fx2z_F3y;
  abcd[249] = 2.0E0*I_NAI_Hxy3z_F3y_a;
  abcd[250] = 2.0E0*I_NAI_H5y_F3y_a-4*I_NAI_F3y_F3y;
  abcd[251] = 2.0E0*I_NAI_H4yz_F3y_a-3*I_NAI_F2yz_F3y;
  abcd[252] = 2.0E0*I_NAI_H3y2z_F3y_a-2*I_NAI_Fy2z_F3y;
  abcd[253] = 2.0E0*I_NAI_H2y3z_F3y_a-1*I_NAI_F3z_F3y;
  abcd[254] = 2.0E0*I_NAI_Hy4z_F3y_a;
  abcd[255] = 2.0E0*I_NAI_H4xy_F2yz_a;
  abcd[256] = 2.0E0*I_NAI_H3x2y_F2yz_a-1*I_NAI_F3x_F2yz;
  abcd[257] = 2.0E0*I_NAI_H3xyz_F2yz_a;
  abcd[258] = 2.0E0*I_NAI_H2x3y_F2yz_a-2*I_NAI_F2xy_F2yz;
  abcd[259] = 2.0E0*I_NAI_H2x2yz_F2yz_a-1*I_NAI_F2xz_F2yz;
  abcd[260] = 2.0E0*I_NAI_H2xy2z_F2yz_a;
  abcd[261] = 2.0E0*I_NAI_Hx4y_F2yz_a-3*I_NAI_Fx2y_F2yz;
  abcd[262] = 2.0E0*I_NAI_Hx3yz_F2yz_a-2*I_NAI_Fxyz_F2yz;
  abcd[263] = 2.0E0*I_NAI_Hx2y2z_F2yz_a-1*I_NAI_Fx2z_F2yz;
  abcd[264] = 2.0E0*I_NAI_Hxy3z_F2yz_a;
  abcd[265] = 2.0E0*I_NAI_H5y_F2yz_a-4*I_NAI_F3y_F2yz;
  abcd[266] = 2.0E0*I_NAI_H4yz_F2yz_a-3*I_NAI_F2yz_F2yz;
  abcd[267] = 2.0E0*I_NAI_H3y2z_F2yz_a-2*I_NAI_Fy2z_F2yz;
  abcd[268] = 2.0E0*I_NAI_H2y3z_F2yz_a-1*I_NAI_F3z_F2yz;
  abcd[269] = 2.0E0*I_NAI_Hy4z_F2yz_a;
  abcd[270] = 2.0E0*I_NAI_H4xy_Fy2z_a;
  abcd[271] = 2.0E0*I_NAI_H3x2y_Fy2z_a-1*I_NAI_F3x_Fy2z;
  abcd[272] = 2.0E0*I_NAI_H3xyz_Fy2z_a;
  abcd[273] = 2.0E0*I_NAI_H2x3y_Fy2z_a-2*I_NAI_F2xy_Fy2z;
  abcd[274] = 2.0E0*I_NAI_H2x2yz_Fy2z_a-1*I_NAI_F2xz_Fy2z;
  abcd[275] = 2.0E0*I_NAI_H2xy2z_Fy2z_a;
  abcd[276] = 2.0E0*I_NAI_Hx4y_Fy2z_a-3*I_NAI_Fx2y_Fy2z;
  abcd[277] = 2.0E0*I_NAI_Hx3yz_Fy2z_a-2*I_NAI_Fxyz_Fy2z;
  abcd[278] = 2.0E0*I_NAI_Hx2y2z_Fy2z_a-1*I_NAI_Fx2z_Fy2z;
  abcd[279] = 2.0E0*I_NAI_Hxy3z_Fy2z_a;
  abcd[280] = 2.0E0*I_NAI_H5y_Fy2z_a-4*I_NAI_F3y_Fy2z;
  abcd[281] = 2.0E0*I_NAI_H4yz_Fy2z_a-3*I_NAI_F2yz_Fy2z;
  abcd[282] = 2.0E0*I_NAI_H3y2z_Fy2z_a-2*I_NAI_Fy2z_Fy2z;
  abcd[283] = 2.0E0*I_NAI_H2y3z_Fy2z_a-1*I_NAI_F3z_Fy2z;
  abcd[284] = 2.0E0*I_NAI_Hy4z_Fy2z_a;
  abcd[285] = 2.0E0*I_NAI_H4xy_F3z_a;
  abcd[286] = 2.0E0*I_NAI_H3x2y_F3z_a-1*I_NAI_F3x_F3z;
  abcd[287] = 2.0E0*I_NAI_H3xyz_F3z_a;
  abcd[288] = 2.0E0*I_NAI_H2x3y_F3z_a-2*I_NAI_F2xy_F3z;
  abcd[289] = 2.0E0*I_NAI_H2x2yz_F3z_a-1*I_NAI_F2xz_F3z;
  abcd[290] = 2.0E0*I_NAI_H2xy2z_F3z_a;
  abcd[291] = 2.0E0*I_NAI_Hx4y_F3z_a-3*I_NAI_Fx2y_F3z;
  abcd[292] = 2.0E0*I_NAI_Hx3yz_F3z_a-2*I_NAI_Fxyz_F3z;
  abcd[293] = 2.0E0*I_NAI_Hx2y2z_F3z_a-1*I_NAI_Fx2z_F3z;
  abcd[294] = 2.0E0*I_NAI_Hxy3z_F3z_a;
  abcd[295] = 2.0E0*I_NAI_H5y_F3z_a-4*I_NAI_F3y_F3z;
  abcd[296] = 2.0E0*I_NAI_H4yz_F3z_a-3*I_NAI_F2yz_F3z;
  abcd[297] = 2.0E0*I_NAI_H3y2z_F3z_a-2*I_NAI_Fy2z_F3z;
  abcd[298] = 2.0E0*I_NAI_H2y3z_F3z_a-1*I_NAI_F3z_F3z;
  abcd[299] = 2.0E0*I_NAI_Hy4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F_a
   * RHS shell quartet name: SQ_NAI_F_F
   ************************************************************/
  abcd[300] = 2.0E0*I_NAI_H4xz_F3x_a;
  abcd[301] = 2.0E0*I_NAI_H3xyz_F3x_a;
  abcd[302] = 2.0E0*I_NAI_H3x2z_F3x_a-1*I_NAI_F3x_F3x;
  abcd[303] = 2.0E0*I_NAI_H2x2yz_F3x_a;
  abcd[304] = 2.0E0*I_NAI_H2xy2z_F3x_a-1*I_NAI_F2xy_F3x;
  abcd[305] = 2.0E0*I_NAI_H2x3z_F3x_a-2*I_NAI_F2xz_F3x;
  abcd[306] = 2.0E0*I_NAI_Hx3yz_F3x_a;
  abcd[307] = 2.0E0*I_NAI_Hx2y2z_F3x_a-1*I_NAI_Fx2y_F3x;
  abcd[308] = 2.0E0*I_NAI_Hxy3z_F3x_a-2*I_NAI_Fxyz_F3x;
  abcd[309] = 2.0E0*I_NAI_Hx4z_F3x_a-3*I_NAI_Fx2z_F3x;
  abcd[310] = 2.0E0*I_NAI_H4yz_F3x_a;
  abcd[311] = 2.0E0*I_NAI_H3y2z_F3x_a-1*I_NAI_F3y_F3x;
  abcd[312] = 2.0E0*I_NAI_H2y3z_F3x_a-2*I_NAI_F2yz_F3x;
  abcd[313] = 2.0E0*I_NAI_Hy4z_F3x_a-3*I_NAI_Fy2z_F3x;
  abcd[314] = 2.0E0*I_NAI_H5z_F3x_a-4*I_NAI_F3z_F3x;
  abcd[315] = 2.0E0*I_NAI_H4xz_F2xy_a;
  abcd[316] = 2.0E0*I_NAI_H3xyz_F2xy_a;
  abcd[317] = 2.0E0*I_NAI_H3x2z_F2xy_a-1*I_NAI_F3x_F2xy;
  abcd[318] = 2.0E0*I_NAI_H2x2yz_F2xy_a;
  abcd[319] = 2.0E0*I_NAI_H2xy2z_F2xy_a-1*I_NAI_F2xy_F2xy;
  abcd[320] = 2.0E0*I_NAI_H2x3z_F2xy_a-2*I_NAI_F2xz_F2xy;
  abcd[321] = 2.0E0*I_NAI_Hx3yz_F2xy_a;
  abcd[322] = 2.0E0*I_NAI_Hx2y2z_F2xy_a-1*I_NAI_Fx2y_F2xy;
  abcd[323] = 2.0E0*I_NAI_Hxy3z_F2xy_a-2*I_NAI_Fxyz_F2xy;
  abcd[324] = 2.0E0*I_NAI_Hx4z_F2xy_a-3*I_NAI_Fx2z_F2xy;
  abcd[325] = 2.0E0*I_NAI_H4yz_F2xy_a;
  abcd[326] = 2.0E0*I_NAI_H3y2z_F2xy_a-1*I_NAI_F3y_F2xy;
  abcd[327] = 2.0E0*I_NAI_H2y3z_F2xy_a-2*I_NAI_F2yz_F2xy;
  abcd[328] = 2.0E0*I_NAI_Hy4z_F2xy_a-3*I_NAI_Fy2z_F2xy;
  abcd[329] = 2.0E0*I_NAI_H5z_F2xy_a-4*I_NAI_F3z_F2xy;
  abcd[330] = 2.0E0*I_NAI_H4xz_F2xz_a;
  abcd[331] = 2.0E0*I_NAI_H3xyz_F2xz_a;
  abcd[332] = 2.0E0*I_NAI_H3x2z_F2xz_a-1*I_NAI_F3x_F2xz;
  abcd[333] = 2.0E0*I_NAI_H2x2yz_F2xz_a;
  abcd[334] = 2.0E0*I_NAI_H2xy2z_F2xz_a-1*I_NAI_F2xy_F2xz;
  abcd[335] = 2.0E0*I_NAI_H2x3z_F2xz_a-2*I_NAI_F2xz_F2xz;
  abcd[336] = 2.0E0*I_NAI_Hx3yz_F2xz_a;
  abcd[337] = 2.0E0*I_NAI_Hx2y2z_F2xz_a-1*I_NAI_Fx2y_F2xz;
  abcd[338] = 2.0E0*I_NAI_Hxy3z_F2xz_a-2*I_NAI_Fxyz_F2xz;
  abcd[339] = 2.0E0*I_NAI_Hx4z_F2xz_a-3*I_NAI_Fx2z_F2xz;
  abcd[340] = 2.0E0*I_NAI_H4yz_F2xz_a;
  abcd[341] = 2.0E0*I_NAI_H3y2z_F2xz_a-1*I_NAI_F3y_F2xz;
  abcd[342] = 2.0E0*I_NAI_H2y3z_F2xz_a-2*I_NAI_F2yz_F2xz;
  abcd[343] = 2.0E0*I_NAI_Hy4z_F2xz_a-3*I_NAI_Fy2z_F2xz;
  abcd[344] = 2.0E0*I_NAI_H5z_F2xz_a-4*I_NAI_F3z_F2xz;
  abcd[345] = 2.0E0*I_NAI_H4xz_Fx2y_a;
  abcd[346] = 2.0E0*I_NAI_H3xyz_Fx2y_a;
  abcd[347] = 2.0E0*I_NAI_H3x2z_Fx2y_a-1*I_NAI_F3x_Fx2y;
  abcd[348] = 2.0E0*I_NAI_H2x2yz_Fx2y_a;
  abcd[349] = 2.0E0*I_NAI_H2xy2z_Fx2y_a-1*I_NAI_F2xy_Fx2y;
  abcd[350] = 2.0E0*I_NAI_H2x3z_Fx2y_a-2*I_NAI_F2xz_Fx2y;
  abcd[351] = 2.0E0*I_NAI_Hx3yz_Fx2y_a;
  abcd[352] = 2.0E0*I_NAI_Hx2y2z_Fx2y_a-1*I_NAI_Fx2y_Fx2y;
  abcd[353] = 2.0E0*I_NAI_Hxy3z_Fx2y_a-2*I_NAI_Fxyz_Fx2y;
  abcd[354] = 2.0E0*I_NAI_Hx4z_Fx2y_a-3*I_NAI_Fx2z_Fx2y;
  abcd[355] = 2.0E0*I_NAI_H4yz_Fx2y_a;
  abcd[356] = 2.0E0*I_NAI_H3y2z_Fx2y_a-1*I_NAI_F3y_Fx2y;
  abcd[357] = 2.0E0*I_NAI_H2y3z_Fx2y_a-2*I_NAI_F2yz_Fx2y;
  abcd[358] = 2.0E0*I_NAI_Hy4z_Fx2y_a-3*I_NAI_Fy2z_Fx2y;
  abcd[359] = 2.0E0*I_NAI_H5z_Fx2y_a-4*I_NAI_F3z_Fx2y;
  abcd[360] = 2.0E0*I_NAI_H4xz_Fxyz_a;
  abcd[361] = 2.0E0*I_NAI_H3xyz_Fxyz_a;
  abcd[362] = 2.0E0*I_NAI_H3x2z_Fxyz_a-1*I_NAI_F3x_Fxyz;
  abcd[363] = 2.0E0*I_NAI_H2x2yz_Fxyz_a;
  abcd[364] = 2.0E0*I_NAI_H2xy2z_Fxyz_a-1*I_NAI_F2xy_Fxyz;
  abcd[365] = 2.0E0*I_NAI_H2x3z_Fxyz_a-2*I_NAI_F2xz_Fxyz;
  abcd[366] = 2.0E0*I_NAI_Hx3yz_Fxyz_a;
  abcd[367] = 2.0E0*I_NAI_Hx2y2z_Fxyz_a-1*I_NAI_Fx2y_Fxyz;
  abcd[368] = 2.0E0*I_NAI_Hxy3z_Fxyz_a-2*I_NAI_Fxyz_Fxyz;
  abcd[369] = 2.0E0*I_NAI_Hx4z_Fxyz_a-3*I_NAI_Fx2z_Fxyz;
  abcd[370] = 2.0E0*I_NAI_H4yz_Fxyz_a;
  abcd[371] = 2.0E0*I_NAI_H3y2z_Fxyz_a-1*I_NAI_F3y_Fxyz;
  abcd[372] = 2.0E0*I_NAI_H2y3z_Fxyz_a-2*I_NAI_F2yz_Fxyz;
  abcd[373] = 2.0E0*I_NAI_Hy4z_Fxyz_a-3*I_NAI_Fy2z_Fxyz;
  abcd[374] = 2.0E0*I_NAI_H5z_Fxyz_a-4*I_NAI_F3z_Fxyz;
  abcd[375] = 2.0E0*I_NAI_H4xz_Fx2z_a;
  abcd[376] = 2.0E0*I_NAI_H3xyz_Fx2z_a;
  abcd[377] = 2.0E0*I_NAI_H3x2z_Fx2z_a-1*I_NAI_F3x_Fx2z;
  abcd[378] = 2.0E0*I_NAI_H2x2yz_Fx2z_a;
  abcd[379] = 2.0E0*I_NAI_H2xy2z_Fx2z_a-1*I_NAI_F2xy_Fx2z;
  abcd[380] = 2.0E0*I_NAI_H2x3z_Fx2z_a-2*I_NAI_F2xz_Fx2z;
  abcd[381] = 2.0E0*I_NAI_Hx3yz_Fx2z_a;
  abcd[382] = 2.0E0*I_NAI_Hx2y2z_Fx2z_a-1*I_NAI_Fx2y_Fx2z;
  abcd[383] = 2.0E0*I_NAI_Hxy3z_Fx2z_a-2*I_NAI_Fxyz_Fx2z;
  abcd[384] = 2.0E0*I_NAI_Hx4z_Fx2z_a-3*I_NAI_Fx2z_Fx2z;
  abcd[385] = 2.0E0*I_NAI_H4yz_Fx2z_a;
  abcd[386] = 2.0E0*I_NAI_H3y2z_Fx2z_a-1*I_NAI_F3y_Fx2z;
  abcd[387] = 2.0E0*I_NAI_H2y3z_Fx2z_a-2*I_NAI_F2yz_Fx2z;
  abcd[388] = 2.0E0*I_NAI_Hy4z_Fx2z_a-3*I_NAI_Fy2z_Fx2z;
  abcd[389] = 2.0E0*I_NAI_H5z_Fx2z_a-4*I_NAI_F3z_Fx2z;
  abcd[390] = 2.0E0*I_NAI_H4xz_F3y_a;
  abcd[391] = 2.0E0*I_NAI_H3xyz_F3y_a;
  abcd[392] = 2.0E0*I_NAI_H3x2z_F3y_a-1*I_NAI_F3x_F3y;
  abcd[393] = 2.0E0*I_NAI_H2x2yz_F3y_a;
  abcd[394] = 2.0E0*I_NAI_H2xy2z_F3y_a-1*I_NAI_F2xy_F3y;
  abcd[395] = 2.0E0*I_NAI_H2x3z_F3y_a-2*I_NAI_F2xz_F3y;
  abcd[396] = 2.0E0*I_NAI_Hx3yz_F3y_a;
  abcd[397] = 2.0E0*I_NAI_Hx2y2z_F3y_a-1*I_NAI_Fx2y_F3y;
  abcd[398] = 2.0E0*I_NAI_Hxy3z_F3y_a-2*I_NAI_Fxyz_F3y;
  abcd[399] = 2.0E0*I_NAI_Hx4z_F3y_a-3*I_NAI_Fx2z_F3y;
  abcd[400] = 2.0E0*I_NAI_H4yz_F3y_a;
  abcd[401] = 2.0E0*I_NAI_H3y2z_F3y_a-1*I_NAI_F3y_F3y;
  abcd[402] = 2.0E0*I_NAI_H2y3z_F3y_a-2*I_NAI_F2yz_F3y;
  abcd[403] = 2.0E0*I_NAI_Hy4z_F3y_a-3*I_NAI_Fy2z_F3y;
  abcd[404] = 2.0E0*I_NAI_H5z_F3y_a-4*I_NAI_F3z_F3y;
  abcd[405] = 2.0E0*I_NAI_H4xz_F2yz_a;
  abcd[406] = 2.0E0*I_NAI_H3xyz_F2yz_a;
  abcd[407] = 2.0E0*I_NAI_H3x2z_F2yz_a-1*I_NAI_F3x_F2yz;
  abcd[408] = 2.0E0*I_NAI_H2x2yz_F2yz_a;
  abcd[409] = 2.0E0*I_NAI_H2xy2z_F2yz_a-1*I_NAI_F2xy_F2yz;
  abcd[410] = 2.0E0*I_NAI_H2x3z_F2yz_a-2*I_NAI_F2xz_F2yz;
  abcd[411] = 2.0E0*I_NAI_Hx3yz_F2yz_a;
  abcd[412] = 2.0E0*I_NAI_Hx2y2z_F2yz_a-1*I_NAI_Fx2y_F2yz;
  abcd[413] = 2.0E0*I_NAI_Hxy3z_F2yz_a-2*I_NAI_Fxyz_F2yz;
  abcd[414] = 2.0E0*I_NAI_Hx4z_F2yz_a-3*I_NAI_Fx2z_F2yz;
  abcd[415] = 2.0E0*I_NAI_H4yz_F2yz_a;
  abcd[416] = 2.0E0*I_NAI_H3y2z_F2yz_a-1*I_NAI_F3y_F2yz;
  abcd[417] = 2.0E0*I_NAI_H2y3z_F2yz_a-2*I_NAI_F2yz_F2yz;
  abcd[418] = 2.0E0*I_NAI_Hy4z_F2yz_a-3*I_NAI_Fy2z_F2yz;
  abcd[419] = 2.0E0*I_NAI_H5z_F2yz_a-4*I_NAI_F3z_F2yz;
  abcd[420] = 2.0E0*I_NAI_H4xz_Fy2z_a;
  abcd[421] = 2.0E0*I_NAI_H3xyz_Fy2z_a;
  abcd[422] = 2.0E0*I_NAI_H3x2z_Fy2z_a-1*I_NAI_F3x_Fy2z;
  abcd[423] = 2.0E0*I_NAI_H2x2yz_Fy2z_a;
  abcd[424] = 2.0E0*I_NAI_H2xy2z_Fy2z_a-1*I_NAI_F2xy_Fy2z;
  abcd[425] = 2.0E0*I_NAI_H2x3z_Fy2z_a-2*I_NAI_F2xz_Fy2z;
  abcd[426] = 2.0E0*I_NAI_Hx3yz_Fy2z_a;
  abcd[427] = 2.0E0*I_NAI_Hx2y2z_Fy2z_a-1*I_NAI_Fx2y_Fy2z;
  abcd[428] = 2.0E0*I_NAI_Hxy3z_Fy2z_a-2*I_NAI_Fxyz_Fy2z;
  abcd[429] = 2.0E0*I_NAI_Hx4z_Fy2z_a-3*I_NAI_Fx2z_Fy2z;
  abcd[430] = 2.0E0*I_NAI_H4yz_Fy2z_a;
  abcd[431] = 2.0E0*I_NAI_H3y2z_Fy2z_a-1*I_NAI_F3y_Fy2z;
  abcd[432] = 2.0E0*I_NAI_H2y3z_Fy2z_a-2*I_NAI_F2yz_Fy2z;
  abcd[433] = 2.0E0*I_NAI_Hy4z_Fy2z_a-3*I_NAI_Fy2z_Fy2z;
  abcd[434] = 2.0E0*I_NAI_H5z_Fy2z_a-4*I_NAI_F3z_Fy2z;
  abcd[435] = 2.0E0*I_NAI_H4xz_F3z_a;
  abcd[436] = 2.0E0*I_NAI_H3xyz_F3z_a;
  abcd[437] = 2.0E0*I_NAI_H3x2z_F3z_a-1*I_NAI_F3x_F3z;
  abcd[438] = 2.0E0*I_NAI_H2x2yz_F3z_a;
  abcd[439] = 2.0E0*I_NAI_H2xy2z_F3z_a-1*I_NAI_F2xy_F3z;
  abcd[440] = 2.0E0*I_NAI_H2x3z_F3z_a-2*I_NAI_F2xz_F3z;
  abcd[441] = 2.0E0*I_NAI_Hx3yz_F3z_a;
  abcd[442] = 2.0E0*I_NAI_Hx2y2z_F3z_a-1*I_NAI_Fx2y_F3z;
  abcd[443] = 2.0E0*I_NAI_Hxy3z_F3z_a-2*I_NAI_Fxyz_F3z;
  abcd[444] = 2.0E0*I_NAI_Hx4z_F3z_a-3*I_NAI_Fx2z_F3z;
  abcd[445] = 2.0E0*I_NAI_H4yz_F3z_a;
  abcd[446] = 2.0E0*I_NAI_H3y2z_F3z_a-1*I_NAI_F3y_F3z;
  abcd[447] = 2.0E0*I_NAI_H2y3z_F3z_a-2*I_NAI_F2yz_F3z;
  abcd[448] = 2.0E0*I_NAI_Hy4z_F3z_a-3*I_NAI_Fy2z_F3z;
  abcd[449] = 2.0E0*I_NAI_H5z_F3z_a-4*I_NAI_F3z_F3z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_G_b
   * RHS shell quartet name: SQ_NAI_G_D
   ************************************************************/
  abcd[450] = 2.0E0*I_NAI_G4x_G4x_b-3*I_NAI_G4x_D2x;
  abcd[451] = 2.0E0*I_NAI_G3xy_G4x_b-3*I_NAI_G3xy_D2x;
  abcd[452] = 2.0E0*I_NAI_G3xz_G4x_b-3*I_NAI_G3xz_D2x;
  abcd[453] = 2.0E0*I_NAI_G2x2y_G4x_b-3*I_NAI_G2x2y_D2x;
  abcd[454] = 2.0E0*I_NAI_G2xyz_G4x_b-3*I_NAI_G2xyz_D2x;
  abcd[455] = 2.0E0*I_NAI_G2x2z_G4x_b-3*I_NAI_G2x2z_D2x;
  abcd[456] = 2.0E0*I_NAI_Gx3y_G4x_b-3*I_NAI_Gx3y_D2x;
  abcd[457] = 2.0E0*I_NAI_Gx2yz_G4x_b-3*I_NAI_Gx2yz_D2x;
  abcd[458] = 2.0E0*I_NAI_Gxy2z_G4x_b-3*I_NAI_Gxy2z_D2x;
  abcd[459] = 2.0E0*I_NAI_Gx3z_G4x_b-3*I_NAI_Gx3z_D2x;
  abcd[460] = 2.0E0*I_NAI_G4y_G4x_b-3*I_NAI_G4y_D2x;
  abcd[461] = 2.0E0*I_NAI_G3yz_G4x_b-3*I_NAI_G3yz_D2x;
  abcd[462] = 2.0E0*I_NAI_G2y2z_G4x_b-3*I_NAI_G2y2z_D2x;
  abcd[463] = 2.0E0*I_NAI_Gy3z_G4x_b-3*I_NAI_Gy3z_D2x;
  abcd[464] = 2.0E0*I_NAI_G4z_G4x_b-3*I_NAI_G4z_D2x;
  abcd[465] = 2.0E0*I_NAI_G4x_G3xy_b-2*I_NAI_G4x_Dxy;
  abcd[466] = 2.0E0*I_NAI_G3xy_G3xy_b-2*I_NAI_G3xy_Dxy;
  abcd[467] = 2.0E0*I_NAI_G3xz_G3xy_b-2*I_NAI_G3xz_Dxy;
  abcd[468] = 2.0E0*I_NAI_G2x2y_G3xy_b-2*I_NAI_G2x2y_Dxy;
  abcd[469] = 2.0E0*I_NAI_G2xyz_G3xy_b-2*I_NAI_G2xyz_Dxy;
  abcd[470] = 2.0E0*I_NAI_G2x2z_G3xy_b-2*I_NAI_G2x2z_Dxy;
  abcd[471] = 2.0E0*I_NAI_Gx3y_G3xy_b-2*I_NAI_Gx3y_Dxy;
  abcd[472] = 2.0E0*I_NAI_Gx2yz_G3xy_b-2*I_NAI_Gx2yz_Dxy;
  abcd[473] = 2.0E0*I_NAI_Gxy2z_G3xy_b-2*I_NAI_Gxy2z_Dxy;
  abcd[474] = 2.0E0*I_NAI_Gx3z_G3xy_b-2*I_NAI_Gx3z_Dxy;
  abcd[475] = 2.0E0*I_NAI_G4y_G3xy_b-2*I_NAI_G4y_Dxy;
  abcd[476] = 2.0E0*I_NAI_G3yz_G3xy_b-2*I_NAI_G3yz_Dxy;
  abcd[477] = 2.0E0*I_NAI_G2y2z_G3xy_b-2*I_NAI_G2y2z_Dxy;
  abcd[478] = 2.0E0*I_NAI_Gy3z_G3xy_b-2*I_NAI_Gy3z_Dxy;
  abcd[479] = 2.0E0*I_NAI_G4z_G3xy_b-2*I_NAI_G4z_Dxy;
  abcd[480] = 2.0E0*I_NAI_G4x_G3xz_b-2*I_NAI_G4x_Dxz;
  abcd[481] = 2.0E0*I_NAI_G3xy_G3xz_b-2*I_NAI_G3xy_Dxz;
  abcd[482] = 2.0E0*I_NAI_G3xz_G3xz_b-2*I_NAI_G3xz_Dxz;
  abcd[483] = 2.0E0*I_NAI_G2x2y_G3xz_b-2*I_NAI_G2x2y_Dxz;
  abcd[484] = 2.0E0*I_NAI_G2xyz_G3xz_b-2*I_NAI_G2xyz_Dxz;
  abcd[485] = 2.0E0*I_NAI_G2x2z_G3xz_b-2*I_NAI_G2x2z_Dxz;
  abcd[486] = 2.0E0*I_NAI_Gx3y_G3xz_b-2*I_NAI_Gx3y_Dxz;
  abcd[487] = 2.0E0*I_NAI_Gx2yz_G3xz_b-2*I_NAI_Gx2yz_Dxz;
  abcd[488] = 2.0E0*I_NAI_Gxy2z_G3xz_b-2*I_NAI_Gxy2z_Dxz;
  abcd[489] = 2.0E0*I_NAI_Gx3z_G3xz_b-2*I_NAI_Gx3z_Dxz;
  abcd[490] = 2.0E0*I_NAI_G4y_G3xz_b-2*I_NAI_G4y_Dxz;
  abcd[491] = 2.0E0*I_NAI_G3yz_G3xz_b-2*I_NAI_G3yz_Dxz;
  abcd[492] = 2.0E0*I_NAI_G2y2z_G3xz_b-2*I_NAI_G2y2z_Dxz;
  abcd[493] = 2.0E0*I_NAI_Gy3z_G3xz_b-2*I_NAI_Gy3z_Dxz;
  abcd[494] = 2.0E0*I_NAI_G4z_G3xz_b-2*I_NAI_G4z_Dxz;
  abcd[495] = 2.0E0*I_NAI_G4x_G2x2y_b-1*I_NAI_G4x_D2y;
  abcd[496] = 2.0E0*I_NAI_G3xy_G2x2y_b-1*I_NAI_G3xy_D2y;
  abcd[497] = 2.0E0*I_NAI_G3xz_G2x2y_b-1*I_NAI_G3xz_D2y;
  abcd[498] = 2.0E0*I_NAI_G2x2y_G2x2y_b-1*I_NAI_G2x2y_D2y;
  abcd[499] = 2.0E0*I_NAI_G2xyz_G2x2y_b-1*I_NAI_G2xyz_D2y;
  abcd[500] = 2.0E0*I_NAI_G2x2z_G2x2y_b-1*I_NAI_G2x2z_D2y;
  abcd[501] = 2.0E0*I_NAI_Gx3y_G2x2y_b-1*I_NAI_Gx3y_D2y;
  abcd[502] = 2.0E0*I_NAI_Gx2yz_G2x2y_b-1*I_NAI_Gx2yz_D2y;
  abcd[503] = 2.0E0*I_NAI_Gxy2z_G2x2y_b-1*I_NAI_Gxy2z_D2y;
  abcd[504] = 2.0E0*I_NAI_Gx3z_G2x2y_b-1*I_NAI_Gx3z_D2y;
  abcd[505] = 2.0E0*I_NAI_G4y_G2x2y_b-1*I_NAI_G4y_D2y;
  abcd[506] = 2.0E0*I_NAI_G3yz_G2x2y_b-1*I_NAI_G3yz_D2y;
  abcd[507] = 2.0E0*I_NAI_G2y2z_G2x2y_b-1*I_NAI_G2y2z_D2y;
  abcd[508] = 2.0E0*I_NAI_Gy3z_G2x2y_b-1*I_NAI_Gy3z_D2y;
  abcd[509] = 2.0E0*I_NAI_G4z_G2x2y_b-1*I_NAI_G4z_D2y;
  abcd[510] = 2.0E0*I_NAI_G4x_G2xyz_b-1*I_NAI_G4x_Dyz;
  abcd[511] = 2.0E0*I_NAI_G3xy_G2xyz_b-1*I_NAI_G3xy_Dyz;
  abcd[512] = 2.0E0*I_NAI_G3xz_G2xyz_b-1*I_NAI_G3xz_Dyz;
  abcd[513] = 2.0E0*I_NAI_G2x2y_G2xyz_b-1*I_NAI_G2x2y_Dyz;
  abcd[514] = 2.0E0*I_NAI_G2xyz_G2xyz_b-1*I_NAI_G2xyz_Dyz;
  abcd[515] = 2.0E0*I_NAI_G2x2z_G2xyz_b-1*I_NAI_G2x2z_Dyz;
  abcd[516] = 2.0E0*I_NAI_Gx3y_G2xyz_b-1*I_NAI_Gx3y_Dyz;
  abcd[517] = 2.0E0*I_NAI_Gx2yz_G2xyz_b-1*I_NAI_Gx2yz_Dyz;
  abcd[518] = 2.0E0*I_NAI_Gxy2z_G2xyz_b-1*I_NAI_Gxy2z_Dyz;
  abcd[519] = 2.0E0*I_NAI_Gx3z_G2xyz_b-1*I_NAI_Gx3z_Dyz;
  abcd[520] = 2.0E0*I_NAI_G4y_G2xyz_b-1*I_NAI_G4y_Dyz;
  abcd[521] = 2.0E0*I_NAI_G3yz_G2xyz_b-1*I_NAI_G3yz_Dyz;
  abcd[522] = 2.0E0*I_NAI_G2y2z_G2xyz_b-1*I_NAI_G2y2z_Dyz;
  abcd[523] = 2.0E0*I_NAI_Gy3z_G2xyz_b-1*I_NAI_Gy3z_Dyz;
  abcd[524] = 2.0E0*I_NAI_G4z_G2xyz_b-1*I_NAI_G4z_Dyz;
  abcd[525] = 2.0E0*I_NAI_G4x_G2x2z_b-1*I_NAI_G4x_D2z;
  abcd[526] = 2.0E0*I_NAI_G3xy_G2x2z_b-1*I_NAI_G3xy_D2z;
  abcd[527] = 2.0E0*I_NAI_G3xz_G2x2z_b-1*I_NAI_G3xz_D2z;
  abcd[528] = 2.0E0*I_NAI_G2x2y_G2x2z_b-1*I_NAI_G2x2y_D2z;
  abcd[529] = 2.0E0*I_NAI_G2xyz_G2x2z_b-1*I_NAI_G2xyz_D2z;
  abcd[530] = 2.0E0*I_NAI_G2x2z_G2x2z_b-1*I_NAI_G2x2z_D2z;
  abcd[531] = 2.0E0*I_NAI_Gx3y_G2x2z_b-1*I_NAI_Gx3y_D2z;
  abcd[532] = 2.0E0*I_NAI_Gx2yz_G2x2z_b-1*I_NAI_Gx2yz_D2z;
  abcd[533] = 2.0E0*I_NAI_Gxy2z_G2x2z_b-1*I_NAI_Gxy2z_D2z;
  abcd[534] = 2.0E0*I_NAI_Gx3z_G2x2z_b-1*I_NAI_Gx3z_D2z;
  abcd[535] = 2.0E0*I_NAI_G4y_G2x2z_b-1*I_NAI_G4y_D2z;
  abcd[536] = 2.0E0*I_NAI_G3yz_G2x2z_b-1*I_NAI_G3yz_D2z;
  abcd[537] = 2.0E0*I_NAI_G2y2z_G2x2z_b-1*I_NAI_G2y2z_D2z;
  abcd[538] = 2.0E0*I_NAI_Gy3z_G2x2z_b-1*I_NAI_Gy3z_D2z;
  abcd[539] = 2.0E0*I_NAI_G4z_G2x2z_b-1*I_NAI_G4z_D2z;
  abcd[540] = 2.0E0*I_NAI_G4x_Gx3y_b;
  abcd[541] = 2.0E0*I_NAI_G3xy_Gx3y_b;
  abcd[542] = 2.0E0*I_NAI_G3xz_Gx3y_b;
  abcd[543] = 2.0E0*I_NAI_G2x2y_Gx3y_b;
  abcd[544] = 2.0E0*I_NAI_G2xyz_Gx3y_b;
  abcd[545] = 2.0E0*I_NAI_G2x2z_Gx3y_b;
  abcd[546] = 2.0E0*I_NAI_Gx3y_Gx3y_b;
  abcd[547] = 2.0E0*I_NAI_Gx2yz_Gx3y_b;
  abcd[548] = 2.0E0*I_NAI_Gxy2z_Gx3y_b;
  abcd[549] = 2.0E0*I_NAI_Gx3z_Gx3y_b;
  abcd[550] = 2.0E0*I_NAI_G4y_Gx3y_b;
  abcd[551] = 2.0E0*I_NAI_G3yz_Gx3y_b;
  abcd[552] = 2.0E0*I_NAI_G2y2z_Gx3y_b;
  abcd[553] = 2.0E0*I_NAI_Gy3z_Gx3y_b;
  abcd[554] = 2.0E0*I_NAI_G4z_Gx3y_b;
  abcd[555] = 2.0E0*I_NAI_G4x_Gx2yz_b;
  abcd[556] = 2.0E0*I_NAI_G3xy_Gx2yz_b;
  abcd[557] = 2.0E0*I_NAI_G3xz_Gx2yz_b;
  abcd[558] = 2.0E0*I_NAI_G2x2y_Gx2yz_b;
  abcd[559] = 2.0E0*I_NAI_G2xyz_Gx2yz_b;
  abcd[560] = 2.0E0*I_NAI_G2x2z_Gx2yz_b;
  abcd[561] = 2.0E0*I_NAI_Gx3y_Gx2yz_b;
  abcd[562] = 2.0E0*I_NAI_Gx2yz_Gx2yz_b;
  abcd[563] = 2.0E0*I_NAI_Gxy2z_Gx2yz_b;
  abcd[564] = 2.0E0*I_NAI_Gx3z_Gx2yz_b;
  abcd[565] = 2.0E0*I_NAI_G4y_Gx2yz_b;
  abcd[566] = 2.0E0*I_NAI_G3yz_Gx2yz_b;
  abcd[567] = 2.0E0*I_NAI_G2y2z_Gx2yz_b;
  abcd[568] = 2.0E0*I_NAI_Gy3z_Gx2yz_b;
  abcd[569] = 2.0E0*I_NAI_G4z_Gx2yz_b;
  abcd[570] = 2.0E0*I_NAI_G4x_Gxy2z_b;
  abcd[571] = 2.0E0*I_NAI_G3xy_Gxy2z_b;
  abcd[572] = 2.0E0*I_NAI_G3xz_Gxy2z_b;
  abcd[573] = 2.0E0*I_NAI_G2x2y_Gxy2z_b;
  abcd[574] = 2.0E0*I_NAI_G2xyz_Gxy2z_b;
  abcd[575] = 2.0E0*I_NAI_G2x2z_Gxy2z_b;
  abcd[576] = 2.0E0*I_NAI_Gx3y_Gxy2z_b;
  abcd[577] = 2.0E0*I_NAI_Gx2yz_Gxy2z_b;
  abcd[578] = 2.0E0*I_NAI_Gxy2z_Gxy2z_b;
  abcd[579] = 2.0E0*I_NAI_Gx3z_Gxy2z_b;
  abcd[580] = 2.0E0*I_NAI_G4y_Gxy2z_b;
  abcd[581] = 2.0E0*I_NAI_G3yz_Gxy2z_b;
  abcd[582] = 2.0E0*I_NAI_G2y2z_Gxy2z_b;
  abcd[583] = 2.0E0*I_NAI_Gy3z_Gxy2z_b;
  abcd[584] = 2.0E0*I_NAI_G4z_Gxy2z_b;
  abcd[585] = 2.0E0*I_NAI_G4x_Gx3z_b;
  abcd[586] = 2.0E0*I_NAI_G3xy_Gx3z_b;
  abcd[587] = 2.0E0*I_NAI_G3xz_Gx3z_b;
  abcd[588] = 2.0E0*I_NAI_G2x2y_Gx3z_b;
  abcd[589] = 2.0E0*I_NAI_G2xyz_Gx3z_b;
  abcd[590] = 2.0E0*I_NAI_G2x2z_Gx3z_b;
  abcd[591] = 2.0E0*I_NAI_Gx3y_Gx3z_b;
  abcd[592] = 2.0E0*I_NAI_Gx2yz_Gx3z_b;
  abcd[593] = 2.0E0*I_NAI_Gxy2z_Gx3z_b;
  abcd[594] = 2.0E0*I_NAI_Gx3z_Gx3z_b;
  abcd[595] = 2.0E0*I_NAI_G4y_Gx3z_b;
  abcd[596] = 2.0E0*I_NAI_G3yz_Gx3z_b;
  abcd[597] = 2.0E0*I_NAI_G2y2z_Gx3z_b;
  abcd[598] = 2.0E0*I_NAI_Gy3z_Gx3z_b;
  abcd[599] = 2.0E0*I_NAI_G4z_Gx3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_G_b
   * RHS shell quartet name: SQ_NAI_G_D
   ************************************************************/
  abcd[600] = 2.0E0*I_NAI_G4x_G3xy_b;
  abcd[601] = 2.0E0*I_NAI_G3xy_G3xy_b;
  abcd[602] = 2.0E0*I_NAI_G3xz_G3xy_b;
  abcd[603] = 2.0E0*I_NAI_G2x2y_G3xy_b;
  abcd[604] = 2.0E0*I_NAI_G2xyz_G3xy_b;
  abcd[605] = 2.0E0*I_NAI_G2x2z_G3xy_b;
  abcd[606] = 2.0E0*I_NAI_Gx3y_G3xy_b;
  abcd[607] = 2.0E0*I_NAI_Gx2yz_G3xy_b;
  abcd[608] = 2.0E0*I_NAI_Gxy2z_G3xy_b;
  abcd[609] = 2.0E0*I_NAI_Gx3z_G3xy_b;
  abcd[610] = 2.0E0*I_NAI_G4y_G3xy_b;
  abcd[611] = 2.0E0*I_NAI_G3yz_G3xy_b;
  abcd[612] = 2.0E0*I_NAI_G2y2z_G3xy_b;
  abcd[613] = 2.0E0*I_NAI_Gy3z_G3xy_b;
  abcd[614] = 2.0E0*I_NAI_G4z_G3xy_b;
  abcd[615] = 2.0E0*I_NAI_G4x_G2x2y_b-1*I_NAI_G4x_D2x;
  abcd[616] = 2.0E0*I_NAI_G3xy_G2x2y_b-1*I_NAI_G3xy_D2x;
  abcd[617] = 2.0E0*I_NAI_G3xz_G2x2y_b-1*I_NAI_G3xz_D2x;
  abcd[618] = 2.0E0*I_NAI_G2x2y_G2x2y_b-1*I_NAI_G2x2y_D2x;
  abcd[619] = 2.0E0*I_NAI_G2xyz_G2x2y_b-1*I_NAI_G2xyz_D2x;
  abcd[620] = 2.0E0*I_NAI_G2x2z_G2x2y_b-1*I_NAI_G2x2z_D2x;
  abcd[621] = 2.0E0*I_NAI_Gx3y_G2x2y_b-1*I_NAI_Gx3y_D2x;
  abcd[622] = 2.0E0*I_NAI_Gx2yz_G2x2y_b-1*I_NAI_Gx2yz_D2x;
  abcd[623] = 2.0E0*I_NAI_Gxy2z_G2x2y_b-1*I_NAI_Gxy2z_D2x;
  abcd[624] = 2.0E0*I_NAI_Gx3z_G2x2y_b-1*I_NAI_Gx3z_D2x;
  abcd[625] = 2.0E0*I_NAI_G4y_G2x2y_b-1*I_NAI_G4y_D2x;
  abcd[626] = 2.0E0*I_NAI_G3yz_G2x2y_b-1*I_NAI_G3yz_D2x;
  abcd[627] = 2.0E0*I_NAI_G2y2z_G2x2y_b-1*I_NAI_G2y2z_D2x;
  abcd[628] = 2.0E0*I_NAI_Gy3z_G2x2y_b-1*I_NAI_Gy3z_D2x;
  abcd[629] = 2.0E0*I_NAI_G4z_G2x2y_b-1*I_NAI_G4z_D2x;
  abcd[630] = 2.0E0*I_NAI_G4x_G2xyz_b;
  abcd[631] = 2.0E0*I_NAI_G3xy_G2xyz_b;
  abcd[632] = 2.0E0*I_NAI_G3xz_G2xyz_b;
  abcd[633] = 2.0E0*I_NAI_G2x2y_G2xyz_b;
  abcd[634] = 2.0E0*I_NAI_G2xyz_G2xyz_b;
  abcd[635] = 2.0E0*I_NAI_G2x2z_G2xyz_b;
  abcd[636] = 2.0E0*I_NAI_Gx3y_G2xyz_b;
  abcd[637] = 2.0E0*I_NAI_Gx2yz_G2xyz_b;
  abcd[638] = 2.0E0*I_NAI_Gxy2z_G2xyz_b;
  abcd[639] = 2.0E0*I_NAI_Gx3z_G2xyz_b;
  abcd[640] = 2.0E0*I_NAI_G4y_G2xyz_b;
  abcd[641] = 2.0E0*I_NAI_G3yz_G2xyz_b;
  abcd[642] = 2.0E0*I_NAI_G2y2z_G2xyz_b;
  abcd[643] = 2.0E0*I_NAI_Gy3z_G2xyz_b;
  abcd[644] = 2.0E0*I_NAI_G4z_G2xyz_b;
  abcd[645] = 2.0E0*I_NAI_G4x_Gx3y_b-2*I_NAI_G4x_Dxy;
  abcd[646] = 2.0E0*I_NAI_G3xy_Gx3y_b-2*I_NAI_G3xy_Dxy;
  abcd[647] = 2.0E0*I_NAI_G3xz_Gx3y_b-2*I_NAI_G3xz_Dxy;
  abcd[648] = 2.0E0*I_NAI_G2x2y_Gx3y_b-2*I_NAI_G2x2y_Dxy;
  abcd[649] = 2.0E0*I_NAI_G2xyz_Gx3y_b-2*I_NAI_G2xyz_Dxy;
  abcd[650] = 2.0E0*I_NAI_G2x2z_Gx3y_b-2*I_NAI_G2x2z_Dxy;
  abcd[651] = 2.0E0*I_NAI_Gx3y_Gx3y_b-2*I_NAI_Gx3y_Dxy;
  abcd[652] = 2.0E0*I_NAI_Gx2yz_Gx3y_b-2*I_NAI_Gx2yz_Dxy;
  abcd[653] = 2.0E0*I_NAI_Gxy2z_Gx3y_b-2*I_NAI_Gxy2z_Dxy;
  abcd[654] = 2.0E0*I_NAI_Gx3z_Gx3y_b-2*I_NAI_Gx3z_Dxy;
  abcd[655] = 2.0E0*I_NAI_G4y_Gx3y_b-2*I_NAI_G4y_Dxy;
  abcd[656] = 2.0E0*I_NAI_G3yz_Gx3y_b-2*I_NAI_G3yz_Dxy;
  abcd[657] = 2.0E0*I_NAI_G2y2z_Gx3y_b-2*I_NAI_G2y2z_Dxy;
  abcd[658] = 2.0E0*I_NAI_Gy3z_Gx3y_b-2*I_NAI_Gy3z_Dxy;
  abcd[659] = 2.0E0*I_NAI_G4z_Gx3y_b-2*I_NAI_G4z_Dxy;
  abcd[660] = 2.0E0*I_NAI_G4x_Gx2yz_b-1*I_NAI_G4x_Dxz;
  abcd[661] = 2.0E0*I_NAI_G3xy_Gx2yz_b-1*I_NAI_G3xy_Dxz;
  abcd[662] = 2.0E0*I_NAI_G3xz_Gx2yz_b-1*I_NAI_G3xz_Dxz;
  abcd[663] = 2.0E0*I_NAI_G2x2y_Gx2yz_b-1*I_NAI_G2x2y_Dxz;
  abcd[664] = 2.0E0*I_NAI_G2xyz_Gx2yz_b-1*I_NAI_G2xyz_Dxz;
  abcd[665] = 2.0E0*I_NAI_G2x2z_Gx2yz_b-1*I_NAI_G2x2z_Dxz;
  abcd[666] = 2.0E0*I_NAI_Gx3y_Gx2yz_b-1*I_NAI_Gx3y_Dxz;
  abcd[667] = 2.0E0*I_NAI_Gx2yz_Gx2yz_b-1*I_NAI_Gx2yz_Dxz;
  abcd[668] = 2.0E0*I_NAI_Gxy2z_Gx2yz_b-1*I_NAI_Gxy2z_Dxz;
  abcd[669] = 2.0E0*I_NAI_Gx3z_Gx2yz_b-1*I_NAI_Gx3z_Dxz;
  abcd[670] = 2.0E0*I_NAI_G4y_Gx2yz_b-1*I_NAI_G4y_Dxz;
  abcd[671] = 2.0E0*I_NAI_G3yz_Gx2yz_b-1*I_NAI_G3yz_Dxz;
  abcd[672] = 2.0E0*I_NAI_G2y2z_Gx2yz_b-1*I_NAI_G2y2z_Dxz;
  abcd[673] = 2.0E0*I_NAI_Gy3z_Gx2yz_b-1*I_NAI_Gy3z_Dxz;
  abcd[674] = 2.0E0*I_NAI_G4z_Gx2yz_b-1*I_NAI_G4z_Dxz;
  abcd[675] = 2.0E0*I_NAI_G4x_Gxy2z_b;
  abcd[676] = 2.0E0*I_NAI_G3xy_Gxy2z_b;
  abcd[677] = 2.0E0*I_NAI_G3xz_Gxy2z_b;
  abcd[678] = 2.0E0*I_NAI_G2x2y_Gxy2z_b;
  abcd[679] = 2.0E0*I_NAI_G2xyz_Gxy2z_b;
  abcd[680] = 2.0E0*I_NAI_G2x2z_Gxy2z_b;
  abcd[681] = 2.0E0*I_NAI_Gx3y_Gxy2z_b;
  abcd[682] = 2.0E0*I_NAI_Gx2yz_Gxy2z_b;
  abcd[683] = 2.0E0*I_NAI_Gxy2z_Gxy2z_b;
  abcd[684] = 2.0E0*I_NAI_Gx3z_Gxy2z_b;
  abcd[685] = 2.0E0*I_NAI_G4y_Gxy2z_b;
  abcd[686] = 2.0E0*I_NAI_G3yz_Gxy2z_b;
  abcd[687] = 2.0E0*I_NAI_G2y2z_Gxy2z_b;
  abcd[688] = 2.0E0*I_NAI_Gy3z_Gxy2z_b;
  abcd[689] = 2.0E0*I_NAI_G4z_Gxy2z_b;
  abcd[690] = 2.0E0*I_NAI_G4x_G4y_b-3*I_NAI_G4x_D2y;
  abcd[691] = 2.0E0*I_NAI_G3xy_G4y_b-3*I_NAI_G3xy_D2y;
  abcd[692] = 2.0E0*I_NAI_G3xz_G4y_b-3*I_NAI_G3xz_D2y;
  abcd[693] = 2.0E0*I_NAI_G2x2y_G4y_b-3*I_NAI_G2x2y_D2y;
  abcd[694] = 2.0E0*I_NAI_G2xyz_G4y_b-3*I_NAI_G2xyz_D2y;
  abcd[695] = 2.0E0*I_NAI_G2x2z_G4y_b-3*I_NAI_G2x2z_D2y;
  abcd[696] = 2.0E0*I_NAI_Gx3y_G4y_b-3*I_NAI_Gx3y_D2y;
  abcd[697] = 2.0E0*I_NAI_Gx2yz_G4y_b-3*I_NAI_Gx2yz_D2y;
  abcd[698] = 2.0E0*I_NAI_Gxy2z_G4y_b-3*I_NAI_Gxy2z_D2y;
  abcd[699] = 2.0E0*I_NAI_Gx3z_G4y_b-3*I_NAI_Gx3z_D2y;
  abcd[700] = 2.0E0*I_NAI_G4y_G4y_b-3*I_NAI_G4y_D2y;
  abcd[701] = 2.0E0*I_NAI_G3yz_G4y_b-3*I_NAI_G3yz_D2y;
  abcd[702] = 2.0E0*I_NAI_G2y2z_G4y_b-3*I_NAI_G2y2z_D2y;
  abcd[703] = 2.0E0*I_NAI_Gy3z_G4y_b-3*I_NAI_Gy3z_D2y;
  abcd[704] = 2.0E0*I_NAI_G4z_G4y_b-3*I_NAI_G4z_D2y;
  abcd[705] = 2.0E0*I_NAI_G4x_G3yz_b-2*I_NAI_G4x_Dyz;
  abcd[706] = 2.0E0*I_NAI_G3xy_G3yz_b-2*I_NAI_G3xy_Dyz;
  abcd[707] = 2.0E0*I_NAI_G3xz_G3yz_b-2*I_NAI_G3xz_Dyz;
  abcd[708] = 2.0E0*I_NAI_G2x2y_G3yz_b-2*I_NAI_G2x2y_Dyz;
  abcd[709] = 2.0E0*I_NAI_G2xyz_G3yz_b-2*I_NAI_G2xyz_Dyz;
  abcd[710] = 2.0E0*I_NAI_G2x2z_G3yz_b-2*I_NAI_G2x2z_Dyz;
  abcd[711] = 2.0E0*I_NAI_Gx3y_G3yz_b-2*I_NAI_Gx3y_Dyz;
  abcd[712] = 2.0E0*I_NAI_Gx2yz_G3yz_b-2*I_NAI_Gx2yz_Dyz;
  abcd[713] = 2.0E0*I_NAI_Gxy2z_G3yz_b-2*I_NAI_Gxy2z_Dyz;
  abcd[714] = 2.0E0*I_NAI_Gx3z_G3yz_b-2*I_NAI_Gx3z_Dyz;
  abcd[715] = 2.0E0*I_NAI_G4y_G3yz_b-2*I_NAI_G4y_Dyz;
  abcd[716] = 2.0E0*I_NAI_G3yz_G3yz_b-2*I_NAI_G3yz_Dyz;
  abcd[717] = 2.0E0*I_NAI_G2y2z_G3yz_b-2*I_NAI_G2y2z_Dyz;
  abcd[718] = 2.0E0*I_NAI_Gy3z_G3yz_b-2*I_NAI_Gy3z_Dyz;
  abcd[719] = 2.0E0*I_NAI_G4z_G3yz_b-2*I_NAI_G4z_Dyz;
  abcd[720] = 2.0E0*I_NAI_G4x_G2y2z_b-1*I_NAI_G4x_D2z;
  abcd[721] = 2.0E0*I_NAI_G3xy_G2y2z_b-1*I_NAI_G3xy_D2z;
  abcd[722] = 2.0E0*I_NAI_G3xz_G2y2z_b-1*I_NAI_G3xz_D2z;
  abcd[723] = 2.0E0*I_NAI_G2x2y_G2y2z_b-1*I_NAI_G2x2y_D2z;
  abcd[724] = 2.0E0*I_NAI_G2xyz_G2y2z_b-1*I_NAI_G2xyz_D2z;
  abcd[725] = 2.0E0*I_NAI_G2x2z_G2y2z_b-1*I_NAI_G2x2z_D2z;
  abcd[726] = 2.0E0*I_NAI_Gx3y_G2y2z_b-1*I_NAI_Gx3y_D2z;
  abcd[727] = 2.0E0*I_NAI_Gx2yz_G2y2z_b-1*I_NAI_Gx2yz_D2z;
  abcd[728] = 2.0E0*I_NAI_Gxy2z_G2y2z_b-1*I_NAI_Gxy2z_D2z;
  abcd[729] = 2.0E0*I_NAI_Gx3z_G2y2z_b-1*I_NAI_Gx3z_D2z;
  abcd[730] = 2.0E0*I_NAI_G4y_G2y2z_b-1*I_NAI_G4y_D2z;
  abcd[731] = 2.0E0*I_NAI_G3yz_G2y2z_b-1*I_NAI_G3yz_D2z;
  abcd[732] = 2.0E0*I_NAI_G2y2z_G2y2z_b-1*I_NAI_G2y2z_D2z;
  abcd[733] = 2.0E0*I_NAI_Gy3z_G2y2z_b-1*I_NAI_Gy3z_D2z;
  abcd[734] = 2.0E0*I_NAI_G4z_G2y2z_b-1*I_NAI_G4z_D2z;
  abcd[735] = 2.0E0*I_NAI_G4x_Gy3z_b;
  abcd[736] = 2.0E0*I_NAI_G3xy_Gy3z_b;
  abcd[737] = 2.0E0*I_NAI_G3xz_Gy3z_b;
  abcd[738] = 2.0E0*I_NAI_G2x2y_Gy3z_b;
  abcd[739] = 2.0E0*I_NAI_G2xyz_Gy3z_b;
  abcd[740] = 2.0E0*I_NAI_G2x2z_Gy3z_b;
  abcd[741] = 2.0E0*I_NAI_Gx3y_Gy3z_b;
  abcd[742] = 2.0E0*I_NAI_Gx2yz_Gy3z_b;
  abcd[743] = 2.0E0*I_NAI_Gxy2z_Gy3z_b;
  abcd[744] = 2.0E0*I_NAI_Gx3z_Gy3z_b;
  abcd[745] = 2.0E0*I_NAI_G4y_Gy3z_b;
  abcd[746] = 2.0E0*I_NAI_G3yz_Gy3z_b;
  abcd[747] = 2.0E0*I_NAI_G2y2z_Gy3z_b;
  abcd[748] = 2.0E0*I_NAI_Gy3z_Gy3z_b;
  abcd[749] = 2.0E0*I_NAI_G4z_Gy3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_G_b
   * RHS shell quartet name: SQ_NAI_G_D
   ************************************************************/
  abcd[750] = 2.0E0*I_NAI_G4x_G3xz_b;
  abcd[751] = 2.0E0*I_NAI_G3xy_G3xz_b;
  abcd[752] = 2.0E0*I_NAI_G3xz_G3xz_b;
  abcd[753] = 2.0E0*I_NAI_G2x2y_G3xz_b;
  abcd[754] = 2.0E0*I_NAI_G2xyz_G3xz_b;
  abcd[755] = 2.0E0*I_NAI_G2x2z_G3xz_b;
  abcd[756] = 2.0E0*I_NAI_Gx3y_G3xz_b;
  abcd[757] = 2.0E0*I_NAI_Gx2yz_G3xz_b;
  abcd[758] = 2.0E0*I_NAI_Gxy2z_G3xz_b;
  abcd[759] = 2.0E0*I_NAI_Gx3z_G3xz_b;
  abcd[760] = 2.0E0*I_NAI_G4y_G3xz_b;
  abcd[761] = 2.0E0*I_NAI_G3yz_G3xz_b;
  abcd[762] = 2.0E0*I_NAI_G2y2z_G3xz_b;
  abcd[763] = 2.0E0*I_NAI_Gy3z_G3xz_b;
  abcd[764] = 2.0E0*I_NAI_G4z_G3xz_b;
  abcd[765] = 2.0E0*I_NAI_G4x_G2xyz_b;
  abcd[766] = 2.0E0*I_NAI_G3xy_G2xyz_b;
  abcd[767] = 2.0E0*I_NAI_G3xz_G2xyz_b;
  abcd[768] = 2.0E0*I_NAI_G2x2y_G2xyz_b;
  abcd[769] = 2.0E0*I_NAI_G2xyz_G2xyz_b;
  abcd[770] = 2.0E0*I_NAI_G2x2z_G2xyz_b;
  abcd[771] = 2.0E0*I_NAI_Gx3y_G2xyz_b;
  abcd[772] = 2.0E0*I_NAI_Gx2yz_G2xyz_b;
  abcd[773] = 2.0E0*I_NAI_Gxy2z_G2xyz_b;
  abcd[774] = 2.0E0*I_NAI_Gx3z_G2xyz_b;
  abcd[775] = 2.0E0*I_NAI_G4y_G2xyz_b;
  abcd[776] = 2.0E0*I_NAI_G3yz_G2xyz_b;
  abcd[777] = 2.0E0*I_NAI_G2y2z_G2xyz_b;
  abcd[778] = 2.0E0*I_NAI_Gy3z_G2xyz_b;
  abcd[779] = 2.0E0*I_NAI_G4z_G2xyz_b;
  abcd[780] = 2.0E0*I_NAI_G4x_G2x2z_b-1*I_NAI_G4x_D2x;
  abcd[781] = 2.0E0*I_NAI_G3xy_G2x2z_b-1*I_NAI_G3xy_D2x;
  abcd[782] = 2.0E0*I_NAI_G3xz_G2x2z_b-1*I_NAI_G3xz_D2x;
  abcd[783] = 2.0E0*I_NAI_G2x2y_G2x2z_b-1*I_NAI_G2x2y_D2x;
  abcd[784] = 2.0E0*I_NAI_G2xyz_G2x2z_b-1*I_NAI_G2xyz_D2x;
  abcd[785] = 2.0E0*I_NAI_G2x2z_G2x2z_b-1*I_NAI_G2x2z_D2x;
  abcd[786] = 2.0E0*I_NAI_Gx3y_G2x2z_b-1*I_NAI_Gx3y_D2x;
  abcd[787] = 2.0E0*I_NAI_Gx2yz_G2x2z_b-1*I_NAI_Gx2yz_D2x;
  abcd[788] = 2.0E0*I_NAI_Gxy2z_G2x2z_b-1*I_NAI_Gxy2z_D2x;
  abcd[789] = 2.0E0*I_NAI_Gx3z_G2x2z_b-1*I_NAI_Gx3z_D2x;
  abcd[790] = 2.0E0*I_NAI_G4y_G2x2z_b-1*I_NAI_G4y_D2x;
  abcd[791] = 2.0E0*I_NAI_G3yz_G2x2z_b-1*I_NAI_G3yz_D2x;
  abcd[792] = 2.0E0*I_NAI_G2y2z_G2x2z_b-1*I_NAI_G2y2z_D2x;
  abcd[793] = 2.0E0*I_NAI_Gy3z_G2x2z_b-1*I_NAI_Gy3z_D2x;
  abcd[794] = 2.0E0*I_NAI_G4z_G2x2z_b-1*I_NAI_G4z_D2x;
  abcd[795] = 2.0E0*I_NAI_G4x_Gx2yz_b;
  abcd[796] = 2.0E0*I_NAI_G3xy_Gx2yz_b;
  abcd[797] = 2.0E0*I_NAI_G3xz_Gx2yz_b;
  abcd[798] = 2.0E0*I_NAI_G2x2y_Gx2yz_b;
  abcd[799] = 2.0E0*I_NAI_G2xyz_Gx2yz_b;
  abcd[800] = 2.0E0*I_NAI_G2x2z_Gx2yz_b;
  abcd[801] = 2.0E0*I_NAI_Gx3y_Gx2yz_b;
  abcd[802] = 2.0E0*I_NAI_Gx2yz_Gx2yz_b;
  abcd[803] = 2.0E0*I_NAI_Gxy2z_Gx2yz_b;
  abcd[804] = 2.0E0*I_NAI_Gx3z_Gx2yz_b;
  abcd[805] = 2.0E0*I_NAI_G4y_Gx2yz_b;
  abcd[806] = 2.0E0*I_NAI_G3yz_Gx2yz_b;
  abcd[807] = 2.0E0*I_NAI_G2y2z_Gx2yz_b;
  abcd[808] = 2.0E0*I_NAI_Gy3z_Gx2yz_b;
  abcd[809] = 2.0E0*I_NAI_G4z_Gx2yz_b;
  abcd[810] = 2.0E0*I_NAI_G4x_Gxy2z_b-1*I_NAI_G4x_Dxy;
  abcd[811] = 2.0E0*I_NAI_G3xy_Gxy2z_b-1*I_NAI_G3xy_Dxy;
  abcd[812] = 2.0E0*I_NAI_G3xz_Gxy2z_b-1*I_NAI_G3xz_Dxy;
  abcd[813] = 2.0E0*I_NAI_G2x2y_Gxy2z_b-1*I_NAI_G2x2y_Dxy;
  abcd[814] = 2.0E0*I_NAI_G2xyz_Gxy2z_b-1*I_NAI_G2xyz_Dxy;
  abcd[815] = 2.0E0*I_NAI_G2x2z_Gxy2z_b-1*I_NAI_G2x2z_Dxy;
  abcd[816] = 2.0E0*I_NAI_Gx3y_Gxy2z_b-1*I_NAI_Gx3y_Dxy;
  abcd[817] = 2.0E0*I_NAI_Gx2yz_Gxy2z_b-1*I_NAI_Gx2yz_Dxy;
  abcd[818] = 2.0E0*I_NAI_Gxy2z_Gxy2z_b-1*I_NAI_Gxy2z_Dxy;
  abcd[819] = 2.0E0*I_NAI_Gx3z_Gxy2z_b-1*I_NAI_Gx3z_Dxy;
  abcd[820] = 2.0E0*I_NAI_G4y_Gxy2z_b-1*I_NAI_G4y_Dxy;
  abcd[821] = 2.0E0*I_NAI_G3yz_Gxy2z_b-1*I_NAI_G3yz_Dxy;
  abcd[822] = 2.0E0*I_NAI_G2y2z_Gxy2z_b-1*I_NAI_G2y2z_Dxy;
  abcd[823] = 2.0E0*I_NAI_Gy3z_Gxy2z_b-1*I_NAI_Gy3z_Dxy;
  abcd[824] = 2.0E0*I_NAI_G4z_Gxy2z_b-1*I_NAI_G4z_Dxy;
  abcd[825] = 2.0E0*I_NAI_G4x_Gx3z_b-2*I_NAI_G4x_Dxz;
  abcd[826] = 2.0E0*I_NAI_G3xy_Gx3z_b-2*I_NAI_G3xy_Dxz;
  abcd[827] = 2.0E0*I_NAI_G3xz_Gx3z_b-2*I_NAI_G3xz_Dxz;
  abcd[828] = 2.0E0*I_NAI_G2x2y_Gx3z_b-2*I_NAI_G2x2y_Dxz;
  abcd[829] = 2.0E0*I_NAI_G2xyz_Gx3z_b-2*I_NAI_G2xyz_Dxz;
  abcd[830] = 2.0E0*I_NAI_G2x2z_Gx3z_b-2*I_NAI_G2x2z_Dxz;
  abcd[831] = 2.0E0*I_NAI_Gx3y_Gx3z_b-2*I_NAI_Gx3y_Dxz;
  abcd[832] = 2.0E0*I_NAI_Gx2yz_Gx3z_b-2*I_NAI_Gx2yz_Dxz;
  abcd[833] = 2.0E0*I_NAI_Gxy2z_Gx3z_b-2*I_NAI_Gxy2z_Dxz;
  abcd[834] = 2.0E0*I_NAI_Gx3z_Gx3z_b-2*I_NAI_Gx3z_Dxz;
  abcd[835] = 2.0E0*I_NAI_G4y_Gx3z_b-2*I_NAI_G4y_Dxz;
  abcd[836] = 2.0E0*I_NAI_G3yz_Gx3z_b-2*I_NAI_G3yz_Dxz;
  abcd[837] = 2.0E0*I_NAI_G2y2z_Gx3z_b-2*I_NAI_G2y2z_Dxz;
  abcd[838] = 2.0E0*I_NAI_Gy3z_Gx3z_b-2*I_NAI_Gy3z_Dxz;
  abcd[839] = 2.0E0*I_NAI_G4z_Gx3z_b-2*I_NAI_G4z_Dxz;
  abcd[840] = 2.0E0*I_NAI_G4x_G3yz_b;
  abcd[841] = 2.0E0*I_NAI_G3xy_G3yz_b;
  abcd[842] = 2.0E0*I_NAI_G3xz_G3yz_b;
  abcd[843] = 2.0E0*I_NAI_G2x2y_G3yz_b;
  abcd[844] = 2.0E0*I_NAI_G2xyz_G3yz_b;
  abcd[845] = 2.0E0*I_NAI_G2x2z_G3yz_b;
  abcd[846] = 2.0E0*I_NAI_Gx3y_G3yz_b;
  abcd[847] = 2.0E0*I_NAI_Gx2yz_G3yz_b;
  abcd[848] = 2.0E0*I_NAI_Gxy2z_G3yz_b;
  abcd[849] = 2.0E0*I_NAI_Gx3z_G3yz_b;
  abcd[850] = 2.0E0*I_NAI_G4y_G3yz_b;
  abcd[851] = 2.0E0*I_NAI_G3yz_G3yz_b;
  abcd[852] = 2.0E0*I_NAI_G2y2z_G3yz_b;
  abcd[853] = 2.0E0*I_NAI_Gy3z_G3yz_b;
  abcd[854] = 2.0E0*I_NAI_G4z_G3yz_b;
  abcd[855] = 2.0E0*I_NAI_G4x_G2y2z_b-1*I_NAI_G4x_D2y;
  abcd[856] = 2.0E0*I_NAI_G3xy_G2y2z_b-1*I_NAI_G3xy_D2y;
  abcd[857] = 2.0E0*I_NAI_G3xz_G2y2z_b-1*I_NAI_G3xz_D2y;
  abcd[858] = 2.0E0*I_NAI_G2x2y_G2y2z_b-1*I_NAI_G2x2y_D2y;
  abcd[859] = 2.0E0*I_NAI_G2xyz_G2y2z_b-1*I_NAI_G2xyz_D2y;
  abcd[860] = 2.0E0*I_NAI_G2x2z_G2y2z_b-1*I_NAI_G2x2z_D2y;
  abcd[861] = 2.0E0*I_NAI_Gx3y_G2y2z_b-1*I_NAI_Gx3y_D2y;
  abcd[862] = 2.0E0*I_NAI_Gx2yz_G2y2z_b-1*I_NAI_Gx2yz_D2y;
  abcd[863] = 2.0E0*I_NAI_Gxy2z_G2y2z_b-1*I_NAI_Gxy2z_D2y;
  abcd[864] = 2.0E0*I_NAI_Gx3z_G2y2z_b-1*I_NAI_Gx3z_D2y;
  abcd[865] = 2.0E0*I_NAI_G4y_G2y2z_b-1*I_NAI_G4y_D2y;
  abcd[866] = 2.0E0*I_NAI_G3yz_G2y2z_b-1*I_NAI_G3yz_D2y;
  abcd[867] = 2.0E0*I_NAI_G2y2z_G2y2z_b-1*I_NAI_G2y2z_D2y;
  abcd[868] = 2.0E0*I_NAI_Gy3z_G2y2z_b-1*I_NAI_Gy3z_D2y;
  abcd[869] = 2.0E0*I_NAI_G4z_G2y2z_b-1*I_NAI_G4z_D2y;
  abcd[870] = 2.0E0*I_NAI_G4x_Gy3z_b-2*I_NAI_G4x_Dyz;
  abcd[871] = 2.0E0*I_NAI_G3xy_Gy3z_b-2*I_NAI_G3xy_Dyz;
  abcd[872] = 2.0E0*I_NAI_G3xz_Gy3z_b-2*I_NAI_G3xz_Dyz;
  abcd[873] = 2.0E0*I_NAI_G2x2y_Gy3z_b-2*I_NAI_G2x2y_Dyz;
  abcd[874] = 2.0E0*I_NAI_G2xyz_Gy3z_b-2*I_NAI_G2xyz_Dyz;
  abcd[875] = 2.0E0*I_NAI_G2x2z_Gy3z_b-2*I_NAI_G2x2z_Dyz;
  abcd[876] = 2.0E0*I_NAI_Gx3y_Gy3z_b-2*I_NAI_Gx3y_Dyz;
  abcd[877] = 2.0E0*I_NAI_Gx2yz_Gy3z_b-2*I_NAI_Gx2yz_Dyz;
  abcd[878] = 2.0E0*I_NAI_Gxy2z_Gy3z_b-2*I_NAI_Gxy2z_Dyz;
  abcd[879] = 2.0E0*I_NAI_Gx3z_Gy3z_b-2*I_NAI_Gx3z_Dyz;
  abcd[880] = 2.0E0*I_NAI_G4y_Gy3z_b-2*I_NAI_G4y_Dyz;
  abcd[881] = 2.0E0*I_NAI_G3yz_Gy3z_b-2*I_NAI_G3yz_Dyz;
  abcd[882] = 2.0E0*I_NAI_G2y2z_Gy3z_b-2*I_NAI_G2y2z_Dyz;
  abcd[883] = 2.0E0*I_NAI_Gy3z_Gy3z_b-2*I_NAI_Gy3z_Dyz;
  abcd[884] = 2.0E0*I_NAI_G4z_Gy3z_b-2*I_NAI_G4z_Dyz;
  abcd[885] = 2.0E0*I_NAI_G4x_G4z_b-3*I_NAI_G4x_D2z;
  abcd[886] = 2.0E0*I_NAI_G3xy_G4z_b-3*I_NAI_G3xy_D2z;
  abcd[887] = 2.0E0*I_NAI_G3xz_G4z_b-3*I_NAI_G3xz_D2z;
  abcd[888] = 2.0E0*I_NAI_G2x2y_G4z_b-3*I_NAI_G2x2y_D2z;
  abcd[889] = 2.0E0*I_NAI_G2xyz_G4z_b-3*I_NAI_G2xyz_D2z;
  abcd[890] = 2.0E0*I_NAI_G2x2z_G4z_b-3*I_NAI_G2x2z_D2z;
  abcd[891] = 2.0E0*I_NAI_Gx3y_G4z_b-3*I_NAI_Gx3y_D2z;
  abcd[892] = 2.0E0*I_NAI_Gx2yz_G4z_b-3*I_NAI_Gx2yz_D2z;
  abcd[893] = 2.0E0*I_NAI_Gxy2z_G4z_b-3*I_NAI_Gxy2z_D2z;
  abcd[894] = 2.0E0*I_NAI_Gx3z_G4z_b-3*I_NAI_Gx3z_D2z;
  abcd[895] = 2.0E0*I_NAI_G4y_G4z_b-3*I_NAI_G4y_D2z;
  abcd[896] = 2.0E0*I_NAI_G3yz_G4z_b-3*I_NAI_G3yz_D2z;
  abcd[897] = 2.0E0*I_NAI_G2y2z_G4z_b-3*I_NAI_G2y2z_D2z;
  abcd[898] = 2.0E0*I_NAI_Gy3z_G4z_b-3*I_NAI_Gy3z_D2z;
  abcd[899] = 2.0E0*I_NAI_G4z_G4z_b-3*I_NAI_G4z_D2z;
}
