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

void hgp_os_nai_h_f_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_K7x_S = 0.0E0;
  Double I_NAI_K6xy_S = 0.0E0;
  Double I_NAI_K6xz_S = 0.0E0;
  Double I_NAI_K5x2y_S = 0.0E0;
  Double I_NAI_K5xyz_S = 0.0E0;
  Double I_NAI_K5x2z_S = 0.0E0;
  Double I_NAI_K4x3y_S = 0.0E0;
  Double I_NAI_K4x2yz_S = 0.0E0;
  Double I_NAI_K4xy2z_S = 0.0E0;
  Double I_NAI_K4x3z_S = 0.0E0;
  Double I_NAI_K3x4y_S = 0.0E0;
  Double I_NAI_K3x3yz_S = 0.0E0;
  Double I_NAI_K3x2y2z_S = 0.0E0;
  Double I_NAI_K3xy3z_S = 0.0E0;
  Double I_NAI_K3x4z_S = 0.0E0;
  Double I_NAI_K2x5y_S = 0.0E0;
  Double I_NAI_K2x4yz_S = 0.0E0;
  Double I_NAI_K2x3y2z_S = 0.0E0;
  Double I_NAI_K2x2y3z_S = 0.0E0;
  Double I_NAI_K2xy4z_S = 0.0E0;
  Double I_NAI_K2x5z_S = 0.0E0;
  Double I_NAI_Kx6y_S = 0.0E0;
  Double I_NAI_Kx5yz_S = 0.0E0;
  Double I_NAI_Kx4y2z_S = 0.0E0;
  Double I_NAI_Kx3y3z_S = 0.0E0;
  Double I_NAI_Kx2y4z_S = 0.0E0;
  Double I_NAI_Kxy5z_S = 0.0E0;
  Double I_NAI_Kx6z_S = 0.0E0;
  Double I_NAI_K7y_S = 0.0E0;
  Double I_NAI_K6yz_S = 0.0E0;
  Double I_NAI_K5y2z_S = 0.0E0;
  Double I_NAI_K4y3z_S = 0.0E0;
  Double I_NAI_K3y4z_S = 0.0E0;
  Double I_NAI_K2y5z_S = 0.0E0;
  Double I_NAI_Ky6z_S = 0.0E0;
  Double I_NAI_K7z_S = 0.0E0;
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
  Double I_NAI_M9x_S_a = 0.0E0;
  Double I_NAI_M8xy_S_a = 0.0E0;
  Double I_NAI_M8xz_S_a = 0.0E0;
  Double I_NAI_M7x2y_S_a = 0.0E0;
  Double I_NAI_M7xyz_S_a = 0.0E0;
  Double I_NAI_M7x2z_S_a = 0.0E0;
  Double I_NAI_M6x3y_S_a = 0.0E0;
  Double I_NAI_M6x2yz_S_a = 0.0E0;
  Double I_NAI_M6xy2z_S_a = 0.0E0;
  Double I_NAI_M6x3z_S_a = 0.0E0;
  Double I_NAI_M5x4y_S_a = 0.0E0;
  Double I_NAI_M5x3yz_S_a = 0.0E0;
  Double I_NAI_M5x2y2z_S_a = 0.0E0;
  Double I_NAI_M5xy3z_S_a = 0.0E0;
  Double I_NAI_M5x4z_S_a = 0.0E0;
  Double I_NAI_M4x5y_S_a = 0.0E0;
  Double I_NAI_M4x4yz_S_a = 0.0E0;
  Double I_NAI_M4x3y2z_S_a = 0.0E0;
  Double I_NAI_M4x2y3z_S_a = 0.0E0;
  Double I_NAI_M4xy4z_S_a = 0.0E0;
  Double I_NAI_M4x5z_S_a = 0.0E0;
  Double I_NAI_M3x6y_S_a = 0.0E0;
  Double I_NAI_M3x5yz_S_a = 0.0E0;
  Double I_NAI_M3x4y2z_S_a = 0.0E0;
  Double I_NAI_M3x3y3z_S_a = 0.0E0;
  Double I_NAI_M3x2y4z_S_a = 0.0E0;
  Double I_NAI_M3xy5z_S_a = 0.0E0;
  Double I_NAI_M3x6z_S_a = 0.0E0;
  Double I_NAI_M2x7y_S_a = 0.0E0;
  Double I_NAI_M2x6yz_S_a = 0.0E0;
  Double I_NAI_M2x5y2z_S_a = 0.0E0;
  Double I_NAI_M2x4y3z_S_a = 0.0E0;
  Double I_NAI_M2x3y4z_S_a = 0.0E0;
  Double I_NAI_M2x2y5z_S_a = 0.0E0;
  Double I_NAI_M2xy6z_S_a = 0.0E0;
  Double I_NAI_M2x7z_S_a = 0.0E0;
  Double I_NAI_Mx8y_S_a = 0.0E0;
  Double I_NAI_Mx7yz_S_a = 0.0E0;
  Double I_NAI_Mx6y2z_S_a = 0.0E0;
  Double I_NAI_Mx5y3z_S_a = 0.0E0;
  Double I_NAI_Mx4y4z_S_a = 0.0E0;
  Double I_NAI_Mx3y5z_S_a = 0.0E0;
  Double I_NAI_Mx2y6z_S_a = 0.0E0;
  Double I_NAI_Mxy7z_S_a = 0.0E0;
  Double I_NAI_Mx8z_S_a = 0.0E0;
  Double I_NAI_M9y_S_a = 0.0E0;
  Double I_NAI_M8yz_S_a = 0.0E0;
  Double I_NAI_M7y2z_S_a = 0.0E0;
  Double I_NAI_M6y3z_S_a = 0.0E0;
  Double I_NAI_M5y4z_S_a = 0.0E0;
  Double I_NAI_M4y5z_S_a = 0.0E0;
  Double I_NAI_M3y6z_S_a = 0.0E0;
  Double I_NAI_M2y7z_S_a = 0.0E0;
  Double I_NAI_My8z_S_a = 0.0E0;
  Double I_NAI_M9z_S_a = 0.0E0;
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
  Double I_NAI_M9x_S_b = 0.0E0;
  Double I_NAI_M8xy_S_b = 0.0E0;
  Double I_NAI_M8xz_S_b = 0.0E0;
  Double I_NAI_M7x2y_S_b = 0.0E0;
  Double I_NAI_M7xyz_S_b = 0.0E0;
  Double I_NAI_M7x2z_S_b = 0.0E0;
  Double I_NAI_M6x3y_S_b = 0.0E0;
  Double I_NAI_M6x2yz_S_b = 0.0E0;
  Double I_NAI_M6xy2z_S_b = 0.0E0;
  Double I_NAI_M6x3z_S_b = 0.0E0;
  Double I_NAI_M5x4y_S_b = 0.0E0;
  Double I_NAI_M5x3yz_S_b = 0.0E0;
  Double I_NAI_M5x2y2z_S_b = 0.0E0;
  Double I_NAI_M5xy3z_S_b = 0.0E0;
  Double I_NAI_M5x4z_S_b = 0.0E0;
  Double I_NAI_M4x5y_S_b = 0.0E0;
  Double I_NAI_M4x4yz_S_b = 0.0E0;
  Double I_NAI_M4x3y2z_S_b = 0.0E0;
  Double I_NAI_M4x2y3z_S_b = 0.0E0;
  Double I_NAI_M4xy4z_S_b = 0.0E0;
  Double I_NAI_M4x5z_S_b = 0.0E0;
  Double I_NAI_M3x6y_S_b = 0.0E0;
  Double I_NAI_M3x5yz_S_b = 0.0E0;
  Double I_NAI_M3x4y2z_S_b = 0.0E0;
  Double I_NAI_M3x3y3z_S_b = 0.0E0;
  Double I_NAI_M3x2y4z_S_b = 0.0E0;
  Double I_NAI_M3xy5z_S_b = 0.0E0;
  Double I_NAI_M3x6z_S_b = 0.0E0;
  Double I_NAI_M2x7y_S_b = 0.0E0;
  Double I_NAI_M2x6yz_S_b = 0.0E0;
  Double I_NAI_M2x5y2z_S_b = 0.0E0;
  Double I_NAI_M2x4y3z_S_b = 0.0E0;
  Double I_NAI_M2x3y4z_S_b = 0.0E0;
  Double I_NAI_M2x2y5z_S_b = 0.0E0;
  Double I_NAI_M2xy6z_S_b = 0.0E0;
  Double I_NAI_M2x7z_S_b = 0.0E0;
  Double I_NAI_Mx8y_S_b = 0.0E0;
  Double I_NAI_Mx7yz_S_b = 0.0E0;
  Double I_NAI_Mx6y2z_S_b = 0.0E0;
  Double I_NAI_Mx5y3z_S_b = 0.0E0;
  Double I_NAI_Mx4y4z_S_b = 0.0E0;
  Double I_NAI_Mx3y5z_S_b = 0.0E0;
  Double I_NAI_Mx2y6z_S_b = 0.0E0;
  Double I_NAI_Mxy7z_S_b = 0.0E0;
  Double I_NAI_Mx8z_S_b = 0.0E0;
  Double I_NAI_M9y_S_b = 0.0E0;
  Double I_NAI_M8yz_S_b = 0.0E0;
  Double I_NAI_M7y2z_S_b = 0.0E0;
  Double I_NAI_M6y3z_S_b = 0.0E0;
  Double I_NAI_M5y4z_S_b = 0.0E0;
  Double I_NAI_M4y5z_S_b = 0.0E0;
  Double I_NAI_M3y6z_S_b = 0.0E0;
  Double I_NAI_M2y7z_S_b = 0.0E0;
  Double I_NAI_My8z_S_b = 0.0E0;
  Double I_NAI_M9z_S_b = 0.0E0;
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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_D2x_S_M1_vrr = PAX*I_NAI_Px_S_M1_vrr-PNX*I_NAI_Px_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_D2z_S_vrr = PAZ*I_NAI_Pz_S_vrr-PNZ*I_NAI_Pz_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       ************************************************************/
      Double I_NAI_F3x_S_vrr = PAX*I_NAI_D2x_S_vrr-PNX*I_NAI_D2x_S_M1_vrr+2*oned2z*I_NAI_Px_S_vrr-2*oned2z*I_NAI_Px_S_M1_vrr;
      Double I_NAI_F2xy_S_vrr = PAY*I_NAI_D2x_S_vrr-PNY*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F2xz_S_vrr = PAZ*I_NAI_D2x_S_vrr-PNZ*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Fx2y_S_vrr = PAX*I_NAI_D2y_S_vrr-PNX*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fx2z_S_vrr = PAX*I_NAI_D2z_S_vrr-PNX*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3y_S_vrr = PAY*I_NAI_D2y_S_vrr-PNY*I_NAI_D2y_S_M1_vrr+2*oned2z*I_NAI_Py_S_vrr-2*oned2z*I_NAI_Py_S_M1_vrr;
      Double I_NAI_F2yz_S_vrr = PAZ*I_NAI_D2y_S_vrr-PNZ*I_NAI_D2y_S_M1_vrr;
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
       * shell quartet name: SQ_NAI_K_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_K7x_S += I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S += I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S += I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S += I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S += I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S += I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S += I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S += I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S += I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S += I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S += I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S += I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S += I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S += I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S += I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S += I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S += I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S += I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S += I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S += I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S += I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S += I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S += I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S += I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S += I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S += I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S += I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S += I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S += I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S += I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S += I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S += I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S += I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S += I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S += I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S += I_NAI_K7z_S_vrr;

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
       * shell quartet name: SQ_NAI_M_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_M_S_a_coefs = alpha;
      I_NAI_M9x_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M9x_S_vrr;
      I_NAI_M8xy_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M8xy_S_vrr;
      I_NAI_M8xz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M8xz_S_vrr;
      I_NAI_M7x2y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M7x2y_S_vrr;
      I_NAI_M7xyz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M7xyz_S_vrr;
      I_NAI_M7x2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M7x2z_S_vrr;
      I_NAI_M6x3y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M6x3y_S_vrr;
      I_NAI_M6x2yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M6x2yz_S_vrr;
      I_NAI_M6xy2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M6xy2z_S_vrr;
      I_NAI_M6x3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M6x3z_S_vrr;
      I_NAI_M5x4y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5x4y_S_vrr;
      I_NAI_M5x3yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5x3yz_S_vrr;
      I_NAI_M5x2y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5x2y2z_S_vrr;
      I_NAI_M5xy3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5xy3z_S_vrr;
      I_NAI_M5x4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5x4z_S_vrr;
      I_NAI_M4x5y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4x5y_S_vrr;
      I_NAI_M4x4yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4x4yz_S_vrr;
      I_NAI_M4x3y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4x3y2z_S_vrr;
      I_NAI_M4x2y3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4x2y3z_S_vrr;
      I_NAI_M4xy4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4xy4z_S_vrr;
      I_NAI_M4x5z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4x5z_S_vrr;
      I_NAI_M3x6y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x6y_S_vrr;
      I_NAI_M3x5yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x5yz_S_vrr;
      I_NAI_M3x4y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x4y2z_S_vrr;
      I_NAI_M3x3y3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x3y3z_S_vrr;
      I_NAI_M3x2y4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x2y4z_S_vrr;
      I_NAI_M3xy5z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3xy5z_S_vrr;
      I_NAI_M3x6z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3x6z_S_vrr;
      I_NAI_M2x7y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x7y_S_vrr;
      I_NAI_M2x6yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x6yz_S_vrr;
      I_NAI_M2x5y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x5y2z_S_vrr;
      I_NAI_M2x4y3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x4y3z_S_vrr;
      I_NAI_M2x3y4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x3y4z_S_vrr;
      I_NAI_M2x2y5z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x2y5z_S_vrr;
      I_NAI_M2xy6z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2xy6z_S_vrr;
      I_NAI_M2x7z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2x7z_S_vrr;
      I_NAI_Mx8y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx8y_S_vrr;
      I_NAI_Mx7yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx7yz_S_vrr;
      I_NAI_Mx6y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx6y2z_S_vrr;
      I_NAI_Mx5y3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx5y3z_S_vrr;
      I_NAI_Mx4y4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx4y4z_S_vrr;
      I_NAI_Mx3y5z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx3y5z_S_vrr;
      I_NAI_Mx2y6z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx2y6z_S_vrr;
      I_NAI_Mxy7z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mxy7z_S_vrr;
      I_NAI_Mx8z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_Mx8z_S_vrr;
      I_NAI_M9y_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M9y_S_vrr;
      I_NAI_M8yz_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M8yz_S_vrr;
      I_NAI_M7y2z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M7y2z_S_vrr;
      I_NAI_M6y3z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M6y3z_S_vrr;
      I_NAI_M5y4z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M5y4z_S_vrr;
      I_NAI_M4y5z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M4y5z_S_vrr;
      I_NAI_M3y6z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M3y6z_S_vrr;
      I_NAI_M2y7z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M2y7z_S_vrr;
      I_NAI_My8z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_My8z_S_vrr;
      I_NAI_M9z_S_a += SQ_NAI_M_S_a_coefs*I_NAI_M9z_S_vrr;

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
       * shell quartet name: SQ_NAI_M_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_M_S_b_coefs = beta;
      I_NAI_M9x_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M9x_S_vrr;
      I_NAI_M8xy_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M8xy_S_vrr;
      I_NAI_M8xz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M8xz_S_vrr;
      I_NAI_M7x2y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M7x2y_S_vrr;
      I_NAI_M7xyz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M7xyz_S_vrr;
      I_NAI_M7x2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M7x2z_S_vrr;
      I_NAI_M6x3y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M6x3y_S_vrr;
      I_NAI_M6x2yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M6x2yz_S_vrr;
      I_NAI_M6xy2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M6xy2z_S_vrr;
      I_NAI_M6x3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M6x3z_S_vrr;
      I_NAI_M5x4y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5x4y_S_vrr;
      I_NAI_M5x3yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5x3yz_S_vrr;
      I_NAI_M5x2y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5x2y2z_S_vrr;
      I_NAI_M5xy3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5xy3z_S_vrr;
      I_NAI_M5x4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5x4z_S_vrr;
      I_NAI_M4x5y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4x5y_S_vrr;
      I_NAI_M4x4yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4x4yz_S_vrr;
      I_NAI_M4x3y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4x3y2z_S_vrr;
      I_NAI_M4x2y3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4x2y3z_S_vrr;
      I_NAI_M4xy4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4xy4z_S_vrr;
      I_NAI_M4x5z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4x5z_S_vrr;
      I_NAI_M3x6y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x6y_S_vrr;
      I_NAI_M3x5yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x5yz_S_vrr;
      I_NAI_M3x4y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x4y2z_S_vrr;
      I_NAI_M3x3y3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x3y3z_S_vrr;
      I_NAI_M3x2y4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x2y4z_S_vrr;
      I_NAI_M3xy5z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3xy5z_S_vrr;
      I_NAI_M3x6z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3x6z_S_vrr;
      I_NAI_M2x7y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x7y_S_vrr;
      I_NAI_M2x6yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x6yz_S_vrr;
      I_NAI_M2x5y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x5y2z_S_vrr;
      I_NAI_M2x4y3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x4y3z_S_vrr;
      I_NAI_M2x3y4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x3y4z_S_vrr;
      I_NAI_M2x2y5z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x2y5z_S_vrr;
      I_NAI_M2xy6z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2xy6z_S_vrr;
      I_NAI_M2x7z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2x7z_S_vrr;
      I_NAI_Mx8y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx8y_S_vrr;
      I_NAI_Mx7yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx7yz_S_vrr;
      I_NAI_Mx6y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx6y2z_S_vrr;
      I_NAI_Mx5y3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx5y3z_S_vrr;
      I_NAI_Mx4y4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx4y4z_S_vrr;
      I_NAI_Mx3y5z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx3y5z_S_vrr;
      I_NAI_Mx2y6z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx2y6z_S_vrr;
      I_NAI_Mxy7z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mxy7z_S_vrr;
      I_NAI_Mx8z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_Mx8z_S_vrr;
      I_NAI_M9y_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M9y_S_vrr;
      I_NAI_M8yz_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M8yz_S_vrr;
      I_NAI_M7y2z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M7y2z_S_vrr;
      I_NAI_M6y3z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M6y3z_S_vrr;
      I_NAI_M5y4z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M5y4z_S_vrr;
      I_NAI_M4y5z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M4y5z_S_vrr;
      I_NAI_M3y6z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M3y6z_S_vrr;
      I_NAI_M2y7z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M2y7z_S_vrr;
      I_NAI_My8z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_My8z_S_vrr;
      I_NAI_M9z_S_b += SQ_NAI_M_S_b_coefs*I_NAI_M9z_S_vrr;

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
   * shell quartet name: SQ_NAI_H_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_NAI_H5x_Py = I_NAI_I5xy_S+ABY*I_NAI_H5x_S;
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
  Double I_NAI_H5x_Pz = I_NAI_I5xz_S+ABZ*I_NAI_H5x_S;
  Double I_NAI_H4xy_Pz = I_NAI_I4xyz_S+ABZ*I_NAI_H4xy_S;
  Double I_NAI_H4xz_Pz = I_NAI_I4x2z_S+ABZ*I_NAI_H4xz_S;
  Double I_NAI_H3x2y_Pz = I_NAI_I3x2yz_S+ABZ*I_NAI_H3x2y_S;
  Double I_NAI_H3xyz_Pz = I_NAI_I3xy2z_S+ABZ*I_NAI_H3xyz_S;
  Double I_NAI_H3x2z_Pz = I_NAI_I3x3z_S+ABZ*I_NAI_H3x2z_S;
  Double I_NAI_H2x3y_Pz = I_NAI_I2x3yz_S+ABZ*I_NAI_H2x3y_S;
  Double I_NAI_H2x2yz_Pz = I_NAI_I2x2y2z_S+ABZ*I_NAI_H2x2yz_S;
  Double I_NAI_H2xy2z_Pz = I_NAI_I2xy3z_S+ABZ*I_NAI_H2xy2z_S;
  Double I_NAI_H2x3z_Pz = I_NAI_I2x4z_S+ABZ*I_NAI_H2x3z_S;
  Double I_NAI_Hx4y_Pz = I_NAI_Ix4yz_S+ABZ*I_NAI_Hx4y_S;
  Double I_NAI_Hx3yz_Pz = I_NAI_Ix3y2z_S+ABZ*I_NAI_Hx3yz_S;
  Double I_NAI_Hx2y2z_Pz = I_NAI_Ix2y3z_S+ABZ*I_NAI_Hx2y2z_S;
  Double I_NAI_Hxy3z_Pz = I_NAI_Ixy4z_S+ABZ*I_NAI_Hxy3z_S;
  Double I_NAI_Hx4z_Pz = I_NAI_Ix5z_S+ABZ*I_NAI_Hx4z_S;
  Double I_NAI_H5y_Pz = I_NAI_I5yz_S+ABZ*I_NAI_H5y_S;
  Double I_NAI_H4yz_Pz = I_NAI_I4y2z_S+ABZ*I_NAI_H4yz_S;
  Double I_NAI_H3y2z_Pz = I_NAI_I3y3z_S+ABZ*I_NAI_H3y2z_S;
  Double I_NAI_H2y3z_Pz = I_NAI_I2y4z_S+ABZ*I_NAI_H2y3z_S;
  Double I_NAI_Hy4z_Pz = I_NAI_Iy5z_S+ABZ*I_NAI_Hy4z_S;
  Double I_NAI_H5z_Pz = I_NAI_I6z_S+ABZ*I_NAI_H5z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
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
   * shell quartet name: SQ_NAI_I_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 8 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S
   * RHS shell quartet name: SQ_NAI_I_S
   ************************************************************/
  Double I_NAI_I6x_Px = I_NAI_K7x_S+ABX*I_NAI_I6x_S;
  Double I_NAI_I5xy_Px = I_NAI_K6xy_S+ABX*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Px = I_NAI_K6xz_S+ABX*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Px = I_NAI_K5x2y_S+ABX*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Px = I_NAI_K5xyz_S+ABX*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Px = I_NAI_K5x2z_S+ABX*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Px = I_NAI_K4x3y_S+ABX*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Px = I_NAI_K4x2yz_S+ABX*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Px = I_NAI_K4xy2z_S+ABX*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Px = I_NAI_K4x3z_S+ABX*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Px = I_NAI_K3x4y_S+ABX*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Px = I_NAI_K3x3yz_S+ABX*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Px = I_NAI_K3x2y2z_S+ABX*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Px = I_NAI_K3xy3z_S+ABX*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Px = I_NAI_K3x4z_S+ABX*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Px = I_NAI_K2x5y_S+ABX*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Px = I_NAI_K2x4yz_S+ABX*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Px = I_NAI_K2x3y2z_S+ABX*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Px = I_NAI_K2x2y3z_S+ABX*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Px = I_NAI_K2xy4z_S+ABX*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Px = I_NAI_K2x5z_S+ABX*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Px = I_NAI_Kx6y_S+ABX*I_NAI_I6y_S;
  Double I_NAI_I5yz_Px = I_NAI_Kx5yz_S+ABX*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Px = I_NAI_Kx4y2z_S+ABX*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Px = I_NAI_Kx3y3z_S+ABX*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Px = I_NAI_Kx2y4z_S+ABX*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Px = I_NAI_Kxy5z_S+ABX*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Px = I_NAI_Kx6z_S+ABX*I_NAI_I6z_S;
  Double I_NAI_I5xy_Py = I_NAI_K5x2y_S+ABY*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Py = I_NAI_K5xyz_S+ABY*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Py = I_NAI_K4x3y_S+ABY*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Py = I_NAI_K4x2yz_S+ABY*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Py = I_NAI_K4xy2z_S+ABY*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Py = I_NAI_K3x4y_S+ABY*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Py = I_NAI_K3x3yz_S+ABY*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Py = I_NAI_K3x2y2z_S+ABY*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Py = I_NAI_K3xy3z_S+ABY*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Py = I_NAI_K2x5y_S+ABY*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Py = I_NAI_K2x4yz_S+ABY*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Py = I_NAI_K2x3y2z_S+ABY*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Py = I_NAI_K2x2y3z_S+ABY*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Py = I_NAI_K2xy4z_S+ABY*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Py = I_NAI_Kx6y_S+ABY*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Py = I_NAI_Kx5yz_S+ABY*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Py = I_NAI_Kx4y2z_S+ABY*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Py = I_NAI_Kx3y3z_S+ABY*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Py = I_NAI_Kx2y4z_S+ABY*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Py = I_NAI_Kxy5z_S+ABY*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Py = I_NAI_K7y_S+ABY*I_NAI_I6y_S;
  Double I_NAI_I5yz_Py = I_NAI_K6yz_S+ABY*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Py = I_NAI_K5y2z_S+ABY*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Py = I_NAI_K4y3z_S+ABY*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Py = I_NAI_K3y4z_S+ABY*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Py = I_NAI_K2y5z_S+ABY*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Py = I_NAI_Ky6z_S+ABY*I_NAI_I6z_S;
  Double I_NAI_I5xz_Pz = I_NAI_K5x2z_S+ABZ*I_NAI_I5xz_S;
  Double I_NAI_I4xyz_Pz = I_NAI_K4xy2z_S+ABZ*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Pz = I_NAI_K4x3z_S+ABZ*I_NAI_I4x2z_S;
  Double I_NAI_I3x2yz_Pz = I_NAI_K3x2y2z_S+ABZ*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Pz = I_NAI_K3xy3z_S+ABZ*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Pz = I_NAI_K3x4z_S+ABZ*I_NAI_I3x3z_S;
  Double I_NAI_I2x3yz_Pz = I_NAI_K2x3y2z_S+ABZ*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Pz = I_NAI_K2x2y3z_S+ABZ*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Pz = I_NAI_K2xy4z_S+ABZ*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Pz = I_NAI_K2x5z_S+ABZ*I_NAI_I2x4z_S;
  Double I_NAI_Ix4yz_Pz = I_NAI_Kx4y2z_S+ABZ*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Pz = I_NAI_Kx3y3z_S+ABZ*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Pz = I_NAI_Kx2y4z_S+ABZ*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Pz = I_NAI_Kxy5z_S+ABZ*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Pz = I_NAI_Kx6z_S+ABZ*I_NAI_Ix5z_S;
  Double I_NAI_I5yz_Pz = I_NAI_K5y2z_S+ABZ*I_NAI_I5yz_S;
  Double I_NAI_I4y2z_Pz = I_NAI_K4y3z_S+ABZ*I_NAI_I4y2z_S;
  Double I_NAI_I3y3z_Pz = I_NAI_K3y4z_S+ABZ*I_NAI_I3y3z_S;
  Double I_NAI_I2y4z_Pz = I_NAI_K2y5z_S+ABZ*I_NAI_I2y4z_S;
  Double I_NAI_Iy5z_Pz = I_NAI_Ky6z_S+ABZ*I_NAI_Iy5z_S;
  Double I_NAI_I6z_Pz = I_NAI_K7z_S+ABZ*I_NAI_I6z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P
   * RHS shell quartet name: SQ_NAI_H_P
   ************************************************************/
  Double I_NAI_H5x_D2x = I_NAI_I6x_Px+ABX*I_NAI_H5x_Px;
  Double I_NAI_H4xy_D2x = I_NAI_I5xy_Px+ABX*I_NAI_H4xy_Px;
  Double I_NAI_H4xz_D2x = I_NAI_I5xz_Px+ABX*I_NAI_H4xz_Px;
  Double I_NAI_H3x2y_D2x = I_NAI_I4x2y_Px+ABX*I_NAI_H3x2y_Px;
  Double I_NAI_H3xyz_D2x = I_NAI_I4xyz_Px+ABX*I_NAI_H3xyz_Px;
  Double I_NAI_H3x2z_D2x = I_NAI_I4x2z_Px+ABX*I_NAI_H3x2z_Px;
  Double I_NAI_H2x3y_D2x = I_NAI_I3x3y_Px+ABX*I_NAI_H2x3y_Px;
  Double I_NAI_H2x2yz_D2x = I_NAI_I3x2yz_Px+ABX*I_NAI_H2x2yz_Px;
  Double I_NAI_H2xy2z_D2x = I_NAI_I3xy2z_Px+ABX*I_NAI_H2xy2z_Px;
  Double I_NAI_H2x3z_D2x = I_NAI_I3x3z_Px+ABX*I_NAI_H2x3z_Px;
  Double I_NAI_Hx4y_D2x = I_NAI_I2x4y_Px+ABX*I_NAI_Hx4y_Px;
  Double I_NAI_Hx3yz_D2x = I_NAI_I2x3yz_Px+ABX*I_NAI_Hx3yz_Px;
  Double I_NAI_Hx2y2z_D2x = I_NAI_I2x2y2z_Px+ABX*I_NAI_Hx2y2z_Px;
  Double I_NAI_Hxy3z_D2x = I_NAI_I2xy3z_Px+ABX*I_NAI_Hxy3z_Px;
  Double I_NAI_Hx4z_D2x = I_NAI_I2x4z_Px+ABX*I_NAI_Hx4z_Px;
  Double I_NAI_H5y_D2x = I_NAI_Ix5y_Px+ABX*I_NAI_H5y_Px;
  Double I_NAI_H4yz_D2x = I_NAI_Ix4yz_Px+ABX*I_NAI_H4yz_Px;
  Double I_NAI_H3y2z_D2x = I_NAI_Ix3y2z_Px+ABX*I_NAI_H3y2z_Px;
  Double I_NAI_H2y3z_D2x = I_NAI_Ix2y3z_Px+ABX*I_NAI_H2y3z_Px;
  Double I_NAI_Hy4z_D2x = I_NAI_Ixy4z_Px+ABX*I_NAI_Hy4z_Px;
  Double I_NAI_H5z_D2x = I_NAI_Ix5z_Px+ABX*I_NAI_H5z_Px;
  Double I_NAI_H5x_Dxy = I_NAI_I5xy_Px+ABY*I_NAI_H5x_Px;
  Double I_NAI_H4xy_Dxy = I_NAI_I4x2y_Px+ABY*I_NAI_H4xy_Px;
  Double I_NAI_H4xz_Dxy = I_NAI_I4xyz_Px+ABY*I_NAI_H4xz_Px;
  Double I_NAI_H3x2y_Dxy = I_NAI_I3x3y_Px+ABY*I_NAI_H3x2y_Px;
  Double I_NAI_H3xyz_Dxy = I_NAI_I3x2yz_Px+ABY*I_NAI_H3xyz_Px;
  Double I_NAI_H3x2z_Dxy = I_NAI_I3xy2z_Px+ABY*I_NAI_H3x2z_Px;
  Double I_NAI_H2x3y_Dxy = I_NAI_I2x4y_Px+ABY*I_NAI_H2x3y_Px;
  Double I_NAI_H2x2yz_Dxy = I_NAI_I2x3yz_Px+ABY*I_NAI_H2x2yz_Px;
  Double I_NAI_H2xy2z_Dxy = I_NAI_I2x2y2z_Px+ABY*I_NAI_H2xy2z_Px;
  Double I_NAI_H2x3z_Dxy = I_NAI_I2xy3z_Px+ABY*I_NAI_H2x3z_Px;
  Double I_NAI_Hx4y_Dxy = I_NAI_Ix5y_Px+ABY*I_NAI_Hx4y_Px;
  Double I_NAI_Hx3yz_Dxy = I_NAI_Ix4yz_Px+ABY*I_NAI_Hx3yz_Px;
  Double I_NAI_Hx2y2z_Dxy = I_NAI_Ix3y2z_Px+ABY*I_NAI_Hx2y2z_Px;
  Double I_NAI_Hxy3z_Dxy = I_NAI_Ix2y3z_Px+ABY*I_NAI_Hxy3z_Px;
  Double I_NAI_Hx4z_Dxy = I_NAI_Ixy4z_Px+ABY*I_NAI_Hx4z_Px;
  Double I_NAI_H5y_Dxy = I_NAI_I6y_Px+ABY*I_NAI_H5y_Px;
  Double I_NAI_H4yz_Dxy = I_NAI_I5yz_Px+ABY*I_NAI_H4yz_Px;
  Double I_NAI_H3y2z_Dxy = I_NAI_I4y2z_Px+ABY*I_NAI_H3y2z_Px;
  Double I_NAI_H2y3z_Dxy = I_NAI_I3y3z_Px+ABY*I_NAI_H2y3z_Px;
  Double I_NAI_Hy4z_Dxy = I_NAI_I2y4z_Px+ABY*I_NAI_Hy4z_Px;
  Double I_NAI_H5z_Dxy = I_NAI_Iy5z_Px+ABY*I_NAI_H5z_Px;
  Double I_NAI_H5x_Dxz = I_NAI_I5xz_Px+ABZ*I_NAI_H5x_Px;
  Double I_NAI_H4xy_Dxz = I_NAI_I4xyz_Px+ABZ*I_NAI_H4xy_Px;
  Double I_NAI_H4xz_Dxz = I_NAI_I4x2z_Px+ABZ*I_NAI_H4xz_Px;
  Double I_NAI_H3x2y_Dxz = I_NAI_I3x2yz_Px+ABZ*I_NAI_H3x2y_Px;
  Double I_NAI_H3xyz_Dxz = I_NAI_I3xy2z_Px+ABZ*I_NAI_H3xyz_Px;
  Double I_NAI_H3x2z_Dxz = I_NAI_I3x3z_Px+ABZ*I_NAI_H3x2z_Px;
  Double I_NAI_H2x3y_Dxz = I_NAI_I2x3yz_Px+ABZ*I_NAI_H2x3y_Px;
  Double I_NAI_H2x2yz_Dxz = I_NAI_I2x2y2z_Px+ABZ*I_NAI_H2x2yz_Px;
  Double I_NAI_H2xy2z_Dxz = I_NAI_I2xy3z_Px+ABZ*I_NAI_H2xy2z_Px;
  Double I_NAI_H2x3z_Dxz = I_NAI_I2x4z_Px+ABZ*I_NAI_H2x3z_Px;
  Double I_NAI_Hx4y_Dxz = I_NAI_Ix4yz_Px+ABZ*I_NAI_Hx4y_Px;
  Double I_NAI_Hx3yz_Dxz = I_NAI_Ix3y2z_Px+ABZ*I_NAI_Hx3yz_Px;
  Double I_NAI_Hx2y2z_Dxz = I_NAI_Ix2y3z_Px+ABZ*I_NAI_Hx2y2z_Px;
  Double I_NAI_Hxy3z_Dxz = I_NAI_Ixy4z_Px+ABZ*I_NAI_Hxy3z_Px;
  Double I_NAI_Hx4z_Dxz = I_NAI_Ix5z_Px+ABZ*I_NAI_Hx4z_Px;
  Double I_NAI_H5y_Dxz = I_NAI_I5yz_Px+ABZ*I_NAI_H5y_Px;
  Double I_NAI_H4yz_Dxz = I_NAI_I4y2z_Px+ABZ*I_NAI_H4yz_Px;
  Double I_NAI_H3y2z_Dxz = I_NAI_I3y3z_Px+ABZ*I_NAI_H3y2z_Px;
  Double I_NAI_H2y3z_Dxz = I_NAI_I2y4z_Px+ABZ*I_NAI_H2y3z_Px;
  Double I_NAI_Hy4z_Dxz = I_NAI_Iy5z_Px+ABZ*I_NAI_Hy4z_Px;
  Double I_NAI_H5z_Dxz = I_NAI_I6z_Px+ABZ*I_NAI_H5z_Px;
  Double I_NAI_H5x_D2y = I_NAI_I5xy_Py+ABY*I_NAI_H5x_Py;
  Double I_NAI_H4xy_D2y = I_NAI_I4x2y_Py+ABY*I_NAI_H4xy_Py;
  Double I_NAI_H4xz_D2y = I_NAI_I4xyz_Py+ABY*I_NAI_H4xz_Py;
  Double I_NAI_H3x2y_D2y = I_NAI_I3x3y_Py+ABY*I_NAI_H3x2y_Py;
  Double I_NAI_H3xyz_D2y = I_NAI_I3x2yz_Py+ABY*I_NAI_H3xyz_Py;
  Double I_NAI_H3x2z_D2y = I_NAI_I3xy2z_Py+ABY*I_NAI_H3x2z_Py;
  Double I_NAI_H2x3y_D2y = I_NAI_I2x4y_Py+ABY*I_NAI_H2x3y_Py;
  Double I_NAI_H2x2yz_D2y = I_NAI_I2x3yz_Py+ABY*I_NAI_H2x2yz_Py;
  Double I_NAI_H2xy2z_D2y = I_NAI_I2x2y2z_Py+ABY*I_NAI_H2xy2z_Py;
  Double I_NAI_H2x3z_D2y = I_NAI_I2xy3z_Py+ABY*I_NAI_H2x3z_Py;
  Double I_NAI_Hx4y_D2y = I_NAI_Ix5y_Py+ABY*I_NAI_Hx4y_Py;
  Double I_NAI_Hx3yz_D2y = I_NAI_Ix4yz_Py+ABY*I_NAI_Hx3yz_Py;
  Double I_NAI_Hx2y2z_D2y = I_NAI_Ix3y2z_Py+ABY*I_NAI_Hx2y2z_Py;
  Double I_NAI_Hxy3z_D2y = I_NAI_Ix2y3z_Py+ABY*I_NAI_Hxy3z_Py;
  Double I_NAI_Hx4z_D2y = I_NAI_Ixy4z_Py+ABY*I_NAI_Hx4z_Py;
  Double I_NAI_H5y_D2y = I_NAI_I6y_Py+ABY*I_NAI_H5y_Py;
  Double I_NAI_H4yz_D2y = I_NAI_I5yz_Py+ABY*I_NAI_H4yz_Py;
  Double I_NAI_H3y2z_D2y = I_NAI_I4y2z_Py+ABY*I_NAI_H3y2z_Py;
  Double I_NAI_H2y3z_D2y = I_NAI_I3y3z_Py+ABY*I_NAI_H2y3z_Py;
  Double I_NAI_Hy4z_D2y = I_NAI_I2y4z_Py+ABY*I_NAI_Hy4z_Py;
  Double I_NAI_H5z_D2y = I_NAI_Iy5z_Py+ABY*I_NAI_H5z_Py;
  Double I_NAI_H5x_Dyz = I_NAI_I5xz_Py+ABZ*I_NAI_H5x_Py;
  Double I_NAI_H4xy_Dyz = I_NAI_I4xyz_Py+ABZ*I_NAI_H4xy_Py;
  Double I_NAI_H4xz_Dyz = I_NAI_I4x2z_Py+ABZ*I_NAI_H4xz_Py;
  Double I_NAI_H3x2y_Dyz = I_NAI_I3x2yz_Py+ABZ*I_NAI_H3x2y_Py;
  Double I_NAI_H3xyz_Dyz = I_NAI_I3xy2z_Py+ABZ*I_NAI_H3xyz_Py;
  Double I_NAI_H3x2z_Dyz = I_NAI_I3x3z_Py+ABZ*I_NAI_H3x2z_Py;
  Double I_NAI_H2x3y_Dyz = I_NAI_I2x3yz_Py+ABZ*I_NAI_H2x3y_Py;
  Double I_NAI_H2x2yz_Dyz = I_NAI_I2x2y2z_Py+ABZ*I_NAI_H2x2yz_Py;
  Double I_NAI_H2xy2z_Dyz = I_NAI_I2xy3z_Py+ABZ*I_NAI_H2xy2z_Py;
  Double I_NAI_H2x3z_Dyz = I_NAI_I2x4z_Py+ABZ*I_NAI_H2x3z_Py;
  Double I_NAI_Hx4y_Dyz = I_NAI_Ix4yz_Py+ABZ*I_NAI_Hx4y_Py;
  Double I_NAI_Hx3yz_Dyz = I_NAI_Ix3y2z_Py+ABZ*I_NAI_Hx3yz_Py;
  Double I_NAI_Hx2y2z_Dyz = I_NAI_Ix2y3z_Py+ABZ*I_NAI_Hx2y2z_Py;
  Double I_NAI_Hxy3z_Dyz = I_NAI_Ixy4z_Py+ABZ*I_NAI_Hxy3z_Py;
  Double I_NAI_Hx4z_Dyz = I_NAI_Ix5z_Py+ABZ*I_NAI_Hx4z_Py;
  Double I_NAI_H5y_Dyz = I_NAI_I5yz_Py+ABZ*I_NAI_H5y_Py;
  Double I_NAI_H4yz_Dyz = I_NAI_I4y2z_Py+ABZ*I_NAI_H4yz_Py;
  Double I_NAI_H3y2z_Dyz = I_NAI_I3y3z_Py+ABZ*I_NAI_H3y2z_Py;
  Double I_NAI_H2y3z_Dyz = I_NAI_I2y4z_Py+ABZ*I_NAI_H2y3z_Py;
  Double I_NAI_Hy4z_Dyz = I_NAI_Iy5z_Py+ABZ*I_NAI_Hy4z_Py;
  Double I_NAI_H5z_Dyz = I_NAI_I6z_Py+ABZ*I_NAI_H5z_Py;
  Double I_NAI_H5x_D2z = I_NAI_I5xz_Pz+ABZ*I_NAI_H5x_Pz;
  Double I_NAI_H4xy_D2z = I_NAI_I4xyz_Pz+ABZ*I_NAI_H4xy_Pz;
  Double I_NAI_H4xz_D2z = I_NAI_I4x2z_Pz+ABZ*I_NAI_H4xz_Pz;
  Double I_NAI_H3x2y_D2z = I_NAI_I3x2yz_Pz+ABZ*I_NAI_H3x2y_Pz;
  Double I_NAI_H3xyz_D2z = I_NAI_I3xy2z_Pz+ABZ*I_NAI_H3xyz_Pz;
  Double I_NAI_H3x2z_D2z = I_NAI_I3x3z_Pz+ABZ*I_NAI_H3x2z_Pz;
  Double I_NAI_H2x3y_D2z = I_NAI_I2x3yz_Pz+ABZ*I_NAI_H2x3y_Pz;
  Double I_NAI_H2x2yz_D2z = I_NAI_I2x2y2z_Pz+ABZ*I_NAI_H2x2yz_Pz;
  Double I_NAI_H2xy2z_D2z = I_NAI_I2xy3z_Pz+ABZ*I_NAI_H2xy2z_Pz;
  Double I_NAI_H2x3z_D2z = I_NAI_I2x4z_Pz+ABZ*I_NAI_H2x3z_Pz;
  Double I_NAI_Hx4y_D2z = I_NAI_Ix4yz_Pz+ABZ*I_NAI_Hx4y_Pz;
  Double I_NAI_Hx3yz_D2z = I_NAI_Ix3y2z_Pz+ABZ*I_NAI_Hx3yz_Pz;
  Double I_NAI_Hx2y2z_D2z = I_NAI_Ix2y3z_Pz+ABZ*I_NAI_Hx2y2z_Pz;
  Double I_NAI_Hxy3z_D2z = I_NAI_Ixy4z_Pz+ABZ*I_NAI_Hxy3z_Pz;
  Double I_NAI_Hx4z_D2z = I_NAI_Ix5z_Pz+ABZ*I_NAI_Hx4z_Pz;
  Double I_NAI_H5y_D2z = I_NAI_I5yz_Pz+ABZ*I_NAI_H5y_Pz;
  Double I_NAI_H4yz_D2z = I_NAI_I4y2z_Pz+ABZ*I_NAI_H4yz_Pz;
  Double I_NAI_H3y2z_D2z = I_NAI_I3y3z_Pz+ABZ*I_NAI_H3y2z_Pz;
  Double I_NAI_H2y3z_D2z = I_NAI_I2y4z_Pz+ABZ*I_NAI_H2y3z_Pz;
  Double I_NAI_Hy4z_D2z = I_NAI_Iy5z_Pz+ABZ*I_NAI_Hy4z_Pz;
  Double I_NAI_H5z_D2z = I_NAI_I6z_Pz+ABZ*I_NAI_H5z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D
   * RHS shell quartet name: SQ_NAI_G_D
   ************************************************************/
  Double I_NAI_G4x_F3x = I_NAI_H5x_D2x+ABX*I_NAI_G4x_D2x;
  Double I_NAI_G3xy_F3x = I_NAI_H4xy_D2x+ABX*I_NAI_G3xy_D2x;
  Double I_NAI_G3xz_F3x = I_NAI_H4xz_D2x+ABX*I_NAI_G3xz_D2x;
  Double I_NAI_G2x2y_F3x = I_NAI_H3x2y_D2x+ABX*I_NAI_G2x2y_D2x;
  Double I_NAI_G2xyz_F3x = I_NAI_H3xyz_D2x+ABX*I_NAI_G2xyz_D2x;
  Double I_NAI_G2x2z_F3x = I_NAI_H3x2z_D2x+ABX*I_NAI_G2x2z_D2x;
  Double I_NAI_Gx3y_F3x = I_NAI_H2x3y_D2x+ABX*I_NAI_Gx3y_D2x;
  Double I_NAI_Gx2yz_F3x = I_NAI_H2x2yz_D2x+ABX*I_NAI_Gx2yz_D2x;
  Double I_NAI_Gxy2z_F3x = I_NAI_H2xy2z_D2x+ABX*I_NAI_Gxy2z_D2x;
  Double I_NAI_Gx3z_F3x = I_NAI_H2x3z_D2x+ABX*I_NAI_Gx3z_D2x;
  Double I_NAI_G4y_F3x = I_NAI_Hx4y_D2x+ABX*I_NAI_G4y_D2x;
  Double I_NAI_G3yz_F3x = I_NAI_Hx3yz_D2x+ABX*I_NAI_G3yz_D2x;
  Double I_NAI_G2y2z_F3x = I_NAI_Hx2y2z_D2x+ABX*I_NAI_G2y2z_D2x;
  Double I_NAI_Gy3z_F3x = I_NAI_Hxy3z_D2x+ABX*I_NAI_Gy3z_D2x;
  Double I_NAI_G4z_F3x = I_NAI_Hx4z_D2x+ABX*I_NAI_G4z_D2x;
  Double I_NAI_G4x_F2xy = I_NAI_H4xy_D2x+ABY*I_NAI_G4x_D2x;
  Double I_NAI_G3xy_F2xy = I_NAI_H3x2y_D2x+ABY*I_NAI_G3xy_D2x;
  Double I_NAI_G3xz_F2xy = I_NAI_H3xyz_D2x+ABY*I_NAI_G3xz_D2x;
  Double I_NAI_G2x2y_F2xy = I_NAI_H2x3y_D2x+ABY*I_NAI_G2x2y_D2x;
  Double I_NAI_G2xyz_F2xy = I_NAI_H2x2yz_D2x+ABY*I_NAI_G2xyz_D2x;
  Double I_NAI_G2x2z_F2xy = I_NAI_H2xy2z_D2x+ABY*I_NAI_G2x2z_D2x;
  Double I_NAI_Gx3y_F2xy = I_NAI_Hx4y_D2x+ABY*I_NAI_Gx3y_D2x;
  Double I_NAI_Gx2yz_F2xy = I_NAI_Hx3yz_D2x+ABY*I_NAI_Gx2yz_D2x;
  Double I_NAI_Gxy2z_F2xy = I_NAI_Hx2y2z_D2x+ABY*I_NAI_Gxy2z_D2x;
  Double I_NAI_Gx3z_F2xy = I_NAI_Hxy3z_D2x+ABY*I_NAI_Gx3z_D2x;
  Double I_NAI_G4y_F2xy = I_NAI_H5y_D2x+ABY*I_NAI_G4y_D2x;
  Double I_NAI_G3yz_F2xy = I_NAI_H4yz_D2x+ABY*I_NAI_G3yz_D2x;
  Double I_NAI_G2y2z_F2xy = I_NAI_H3y2z_D2x+ABY*I_NAI_G2y2z_D2x;
  Double I_NAI_Gy3z_F2xy = I_NAI_H2y3z_D2x+ABY*I_NAI_Gy3z_D2x;
  Double I_NAI_G4z_F2xy = I_NAI_Hy4z_D2x+ABY*I_NAI_G4z_D2x;
  Double I_NAI_G4x_F2xz = I_NAI_H4xz_D2x+ABZ*I_NAI_G4x_D2x;
  Double I_NAI_G3xy_F2xz = I_NAI_H3xyz_D2x+ABZ*I_NAI_G3xy_D2x;
  Double I_NAI_G3xz_F2xz = I_NAI_H3x2z_D2x+ABZ*I_NAI_G3xz_D2x;
  Double I_NAI_G2x2y_F2xz = I_NAI_H2x2yz_D2x+ABZ*I_NAI_G2x2y_D2x;
  Double I_NAI_G2xyz_F2xz = I_NAI_H2xy2z_D2x+ABZ*I_NAI_G2xyz_D2x;
  Double I_NAI_G2x2z_F2xz = I_NAI_H2x3z_D2x+ABZ*I_NAI_G2x2z_D2x;
  Double I_NAI_Gx3y_F2xz = I_NAI_Hx3yz_D2x+ABZ*I_NAI_Gx3y_D2x;
  Double I_NAI_Gx2yz_F2xz = I_NAI_Hx2y2z_D2x+ABZ*I_NAI_Gx2yz_D2x;
  Double I_NAI_Gxy2z_F2xz = I_NAI_Hxy3z_D2x+ABZ*I_NAI_Gxy2z_D2x;
  Double I_NAI_Gx3z_F2xz = I_NAI_Hx4z_D2x+ABZ*I_NAI_Gx3z_D2x;
  Double I_NAI_G4y_F2xz = I_NAI_H4yz_D2x+ABZ*I_NAI_G4y_D2x;
  Double I_NAI_G3yz_F2xz = I_NAI_H3y2z_D2x+ABZ*I_NAI_G3yz_D2x;
  Double I_NAI_G2y2z_F2xz = I_NAI_H2y3z_D2x+ABZ*I_NAI_G2y2z_D2x;
  Double I_NAI_Gy3z_F2xz = I_NAI_Hy4z_D2x+ABZ*I_NAI_Gy3z_D2x;
  Double I_NAI_G4z_F2xz = I_NAI_H5z_D2x+ABZ*I_NAI_G4z_D2x;
  Double I_NAI_G4x_Fx2y = I_NAI_H5x_D2y+ABX*I_NAI_G4x_D2y;
  Double I_NAI_G3xy_Fx2y = I_NAI_H4xy_D2y+ABX*I_NAI_G3xy_D2y;
  Double I_NAI_G3xz_Fx2y = I_NAI_H4xz_D2y+ABX*I_NAI_G3xz_D2y;
  Double I_NAI_G2x2y_Fx2y = I_NAI_H3x2y_D2y+ABX*I_NAI_G2x2y_D2y;
  Double I_NAI_G2xyz_Fx2y = I_NAI_H3xyz_D2y+ABX*I_NAI_G2xyz_D2y;
  Double I_NAI_G2x2z_Fx2y = I_NAI_H3x2z_D2y+ABX*I_NAI_G2x2z_D2y;
  Double I_NAI_Gx3y_Fx2y = I_NAI_H2x3y_D2y+ABX*I_NAI_Gx3y_D2y;
  Double I_NAI_Gx2yz_Fx2y = I_NAI_H2x2yz_D2y+ABX*I_NAI_Gx2yz_D2y;
  Double I_NAI_Gxy2z_Fx2y = I_NAI_H2xy2z_D2y+ABX*I_NAI_Gxy2z_D2y;
  Double I_NAI_Gx3z_Fx2y = I_NAI_H2x3z_D2y+ABX*I_NAI_Gx3z_D2y;
  Double I_NAI_G4y_Fx2y = I_NAI_Hx4y_D2y+ABX*I_NAI_G4y_D2y;
  Double I_NAI_G3yz_Fx2y = I_NAI_Hx3yz_D2y+ABX*I_NAI_G3yz_D2y;
  Double I_NAI_G2y2z_Fx2y = I_NAI_Hx2y2z_D2y+ABX*I_NAI_G2y2z_D2y;
  Double I_NAI_Gy3z_Fx2y = I_NAI_Hxy3z_D2y+ABX*I_NAI_Gy3z_D2y;
  Double I_NAI_G4z_Fx2y = I_NAI_Hx4z_D2y+ABX*I_NAI_G4z_D2y;
  Double I_NAI_G4x_Fxyz = I_NAI_H4xz_Dxy+ABZ*I_NAI_G4x_Dxy;
  Double I_NAI_G3xy_Fxyz = I_NAI_H3xyz_Dxy+ABZ*I_NAI_G3xy_Dxy;
  Double I_NAI_G3xz_Fxyz = I_NAI_H3x2z_Dxy+ABZ*I_NAI_G3xz_Dxy;
  Double I_NAI_G2x2y_Fxyz = I_NAI_H2x2yz_Dxy+ABZ*I_NAI_G2x2y_Dxy;
  Double I_NAI_G2xyz_Fxyz = I_NAI_H2xy2z_Dxy+ABZ*I_NAI_G2xyz_Dxy;
  Double I_NAI_G2x2z_Fxyz = I_NAI_H2x3z_Dxy+ABZ*I_NAI_G2x2z_Dxy;
  Double I_NAI_Gx3y_Fxyz = I_NAI_Hx3yz_Dxy+ABZ*I_NAI_Gx3y_Dxy;
  Double I_NAI_Gx2yz_Fxyz = I_NAI_Hx2y2z_Dxy+ABZ*I_NAI_Gx2yz_Dxy;
  Double I_NAI_Gxy2z_Fxyz = I_NAI_Hxy3z_Dxy+ABZ*I_NAI_Gxy2z_Dxy;
  Double I_NAI_Gx3z_Fxyz = I_NAI_Hx4z_Dxy+ABZ*I_NAI_Gx3z_Dxy;
  Double I_NAI_G4y_Fxyz = I_NAI_H4yz_Dxy+ABZ*I_NAI_G4y_Dxy;
  Double I_NAI_G3yz_Fxyz = I_NAI_H3y2z_Dxy+ABZ*I_NAI_G3yz_Dxy;
  Double I_NAI_G2y2z_Fxyz = I_NAI_H2y3z_Dxy+ABZ*I_NAI_G2y2z_Dxy;
  Double I_NAI_Gy3z_Fxyz = I_NAI_Hy4z_Dxy+ABZ*I_NAI_Gy3z_Dxy;
  Double I_NAI_G4z_Fxyz = I_NAI_H5z_Dxy+ABZ*I_NAI_G4z_Dxy;
  Double I_NAI_G4x_Fx2z = I_NAI_H5x_D2z+ABX*I_NAI_G4x_D2z;
  Double I_NAI_G3xy_Fx2z = I_NAI_H4xy_D2z+ABX*I_NAI_G3xy_D2z;
  Double I_NAI_G3xz_Fx2z = I_NAI_H4xz_D2z+ABX*I_NAI_G3xz_D2z;
  Double I_NAI_G2x2y_Fx2z = I_NAI_H3x2y_D2z+ABX*I_NAI_G2x2y_D2z;
  Double I_NAI_G2xyz_Fx2z = I_NAI_H3xyz_D2z+ABX*I_NAI_G2xyz_D2z;
  Double I_NAI_G2x2z_Fx2z = I_NAI_H3x2z_D2z+ABX*I_NAI_G2x2z_D2z;
  Double I_NAI_Gx3y_Fx2z = I_NAI_H2x3y_D2z+ABX*I_NAI_Gx3y_D2z;
  Double I_NAI_Gx2yz_Fx2z = I_NAI_H2x2yz_D2z+ABX*I_NAI_Gx2yz_D2z;
  Double I_NAI_Gxy2z_Fx2z = I_NAI_H2xy2z_D2z+ABX*I_NAI_Gxy2z_D2z;
  Double I_NAI_Gx3z_Fx2z = I_NAI_H2x3z_D2z+ABX*I_NAI_Gx3z_D2z;
  Double I_NAI_G4y_Fx2z = I_NAI_Hx4y_D2z+ABX*I_NAI_G4y_D2z;
  Double I_NAI_G3yz_Fx2z = I_NAI_Hx3yz_D2z+ABX*I_NAI_G3yz_D2z;
  Double I_NAI_G2y2z_Fx2z = I_NAI_Hx2y2z_D2z+ABX*I_NAI_G2y2z_D2z;
  Double I_NAI_Gy3z_Fx2z = I_NAI_Hxy3z_D2z+ABX*I_NAI_Gy3z_D2z;
  Double I_NAI_G4z_Fx2z = I_NAI_Hx4z_D2z+ABX*I_NAI_G4z_D2z;
  Double I_NAI_G4x_F3y = I_NAI_H4xy_D2y+ABY*I_NAI_G4x_D2y;
  Double I_NAI_G3xy_F3y = I_NAI_H3x2y_D2y+ABY*I_NAI_G3xy_D2y;
  Double I_NAI_G3xz_F3y = I_NAI_H3xyz_D2y+ABY*I_NAI_G3xz_D2y;
  Double I_NAI_G2x2y_F3y = I_NAI_H2x3y_D2y+ABY*I_NAI_G2x2y_D2y;
  Double I_NAI_G2xyz_F3y = I_NAI_H2x2yz_D2y+ABY*I_NAI_G2xyz_D2y;
  Double I_NAI_G2x2z_F3y = I_NAI_H2xy2z_D2y+ABY*I_NAI_G2x2z_D2y;
  Double I_NAI_Gx3y_F3y = I_NAI_Hx4y_D2y+ABY*I_NAI_Gx3y_D2y;
  Double I_NAI_Gx2yz_F3y = I_NAI_Hx3yz_D2y+ABY*I_NAI_Gx2yz_D2y;
  Double I_NAI_Gxy2z_F3y = I_NAI_Hx2y2z_D2y+ABY*I_NAI_Gxy2z_D2y;
  Double I_NAI_Gx3z_F3y = I_NAI_Hxy3z_D2y+ABY*I_NAI_Gx3z_D2y;
  Double I_NAI_G4y_F3y = I_NAI_H5y_D2y+ABY*I_NAI_G4y_D2y;
  Double I_NAI_G3yz_F3y = I_NAI_H4yz_D2y+ABY*I_NAI_G3yz_D2y;
  Double I_NAI_G2y2z_F3y = I_NAI_H3y2z_D2y+ABY*I_NAI_G2y2z_D2y;
  Double I_NAI_Gy3z_F3y = I_NAI_H2y3z_D2y+ABY*I_NAI_Gy3z_D2y;
  Double I_NAI_G4z_F3y = I_NAI_Hy4z_D2y+ABY*I_NAI_G4z_D2y;
  Double I_NAI_G4x_F2yz = I_NAI_H4xz_D2y+ABZ*I_NAI_G4x_D2y;
  Double I_NAI_G3xy_F2yz = I_NAI_H3xyz_D2y+ABZ*I_NAI_G3xy_D2y;
  Double I_NAI_G3xz_F2yz = I_NAI_H3x2z_D2y+ABZ*I_NAI_G3xz_D2y;
  Double I_NAI_G2x2y_F2yz = I_NAI_H2x2yz_D2y+ABZ*I_NAI_G2x2y_D2y;
  Double I_NAI_G2xyz_F2yz = I_NAI_H2xy2z_D2y+ABZ*I_NAI_G2xyz_D2y;
  Double I_NAI_G2x2z_F2yz = I_NAI_H2x3z_D2y+ABZ*I_NAI_G2x2z_D2y;
  Double I_NAI_Gx3y_F2yz = I_NAI_Hx3yz_D2y+ABZ*I_NAI_Gx3y_D2y;
  Double I_NAI_Gx2yz_F2yz = I_NAI_Hx2y2z_D2y+ABZ*I_NAI_Gx2yz_D2y;
  Double I_NAI_Gxy2z_F2yz = I_NAI_Hxy3z_D2y+ABZ*I_NAI_Gxy2z_D2y;
  Double I_NAI_Gx3z_F2yz = I_NAI_Hx4z_D2y+ABZ*I_NAI_Gx3z_D2y;
  Double I_NAI_G4y_F2yz = I_NAI_H4yz_D2y+ABZ*I_NAI_G4y_D2y;
  Double I_NAI_G3yz_F2yz = I_NAI_H3y2z_D2y+ABZ*I_NAI_G3yz_D2y;
  Double I_NAI_G2y2z_F2yz = I_NAI_H2y3z_D2y+ABZ*I_NAI_G2y2z_D2y;
  Double I_NAI_Gy3z_F2yz = I_NAI_Hy4z_D2y+ABZ*I_NAI_Gy3z_D2y;
  Double I_NAI_G4z_F2yz = I_NAI_H5z_D2y+ABZ*I_NAI_G4z_D2y;
  Double I_NAI_G4x_Fy2z = I_NAI_H4xy_D2z+ABY*I_NAI_G4x_D2z;
  Double I_NAI_G3xy_Fy2z = I_NAI_H3x2y_D2z+ABY*I_NAI_G3xy_D2z;
  Double I_NAI_G3xz_Fy2z = I_NAI_H3xyz_D2z+ABY*I_NAI_G3xz_D2z;
  Double I_NAI_G2x2y_Fy2z = I_NAI_H2x3y_D2z+ABY*I_NAI_G2x2y_D2z;
  Double I_NAI_G2xyz_Fy2z = I_NAI_H2x2yz_D2z+ABY*I_NAI_G2xyz_D2z;
  Double I_NAI_G2x2z_Fy2z = I_NAI_H2xy2z_D2z+ABY*I_NAI_G2x2z_D2z;
  Double I_NAI_Gx3y_Fy2z = I_NAI_Hx4y_D2z+ABY*I_NAI_Gx3y_D2z;
  Double I_NAI_Gx2yz_Fy2z = I_NAI_Hx3yz_D2z+ABY*I_NAI_Gx2yz_D2z;
  Double I_NAI_Gxy2z_Fy2z = I_NAI_Hx2y2z_D2z+ABY*I_NAI_Gxy2z_D2z;
  Double I_NAI_Gx3z_Fy2z = I_NAI_Hxy3z_D2z+ABY*I_NAI_Gx3z_D2z;
  Double I_NAI_G4y_Fy2z = I_NAI_H5y_D2z+ABY*I_NAI_G4y_D2z;
  Double I_NAI_G3yz_Fy2z = I_NAI_H4yz_D2z+ABY*I_NAI_G3yz_D2z;
  Double I_NAI_G2y2z_Fy2z = I_NAI_H3y2z_D2z+ABY*I_NAI_G2y2z_D2z;
  Double I_NAI_Gy3z_Fy2z = I_NAI_H2y3z_D2z+ABY*I_NAI_Gy3z_D2z;
  Double I_NAI_G4z_Fy2z = I_NAI_Hy4z_D2z+ABY*I_NAI_G4z_D2z;
  Double I_NAI_G4x_F3z = I_NAI_H4xz_D2z+ABZ*I_NAI_G4x_D2z;
  Double I_NAI_G3xy_F3z = I_NAI_H3xyz_D2z+ABZ*I_NAI_G3xy_D2z;
  Double I_NAI_G3xz_F3z = I_NAI_H3x2z_D2z+ABZ*I_NAI_G3xz_D2z;
  Double I_NAI_G2x2y_F3z = I_NAI_H2x2yz_D2z+ABZ*I_NAI_G2x2y_D2z;
  Double I_NAI_G2xyz_F3z = I_NAI_H2xy2z_D2z+ABZ*I_NAI_G2xyz_D2z;
  Double I_NAI_G2x2z_F3z = I_NAI_H2x3z_D2z+ABZ*I_NAI_G2x2z_D2z;
  Double I_NAI_Gx3y_F3z = I_NAI_Hx3yz_D2z+ABZ*I_NAI_Gx3y_D2z;
  Double I_NAI_Gx2yz_F3z = I_NAI_Hx2y2z_D2z+ABZ*I_NAI_Gx2yz_D2z;
  Double I_NAI_Gxy2z_F3z = I_NAI_Hxy3z_D2z+ABZ*I_NAI_Gxy2z_D2z;
  Double I_NAI_Gx3z_F3z = I_NAI_Hx4z_D2z+ABZ*I_NAI_Gx3z_D2z;
  Double I_NAI_G4y_F3z = I_NAI_H4yz_D2z+ABZ*I_NAI_G4y_D2z;
  Double I_NAI_G3yz_F3z = I_NAI_H3y2z_D2z+ABZ*I_NAI_G3yz_D2z;
  Double I_NAI_G2y2z_F3z = I_NAI_H2y3z_D2z+ABZ*I_NAI_G2y2z_D2z;
  Double I_NAI_Gy3z_F3z = I_NAI_Hy4z_D2z+ABZ*I_NAI_Gy3z_D2z;
  Double I_NAI_G4z_F3z = I_NAI_H5z_D2z+ABZ*I_NAI_G4z_D2z;

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
   * shell quartet name: SQ_NAI_K_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_NAI_K7y_Px_a = I_NAI_Lx7y_S_a+ABX*I_NAI_K7y_S_a;
  Double I_NAI_K6yz_Px_a = I_NAI_Lx6yz_S_a+ABX*I_NAI_K6yz_S_a;
  Double I_NAI_K5y2z_Px_a = I_NAI_Lx5y2z_S_a+ABX*I_NAI_K5y2z_S_a;
  Double I_NAI_K4y3z_Px_a = I_NAI_Lx4y3z_S_a+ABX*I_NAI_K4y3z_S_a;
  Double I_NAI_K3y4z_Px_a = I_NAI_Lx3y4z_S_a+ABX*I_NAI_K3y4z_S_a;
  Double I_NAI_K2y5z_Px_a = I_NAI_Lx2y5z_S_a+ABX*I_NAI_K2y5z_S_a;
  Double I_NAI_Ky6z_Px_a = I_NAI_Lxy6z_S_a+ABX*I_NAI_Ky6z_S_a;
  Double I_NAI_K7z_Px_a = I_NAI_Lx7z_S_a+ABX*I_NAI_K7z_S_a;
  Double I_NAI_K7x_Py_a = I_NAI_L7xy_S_a+ABY*I_NAI_K7x_S_a;
  Double I_NAI_K6xy_Py_a = I_NAI_L6x2y_S_a+ABY*I_NAI_K6xy_S_a;
  Double I_NAI_K6xz_Py_a = I_NAI_L6xyz_S_a+ABY*I_NAI_K6xz_S_a;
  Double I_NAI_K5x2y_Py_a = I_NAI_L5x3y_S_a+ABY*I_NAI_K5x2y_S_a;
  Double I_NAI_K5xyz_Py_a = I_NAI_L5x2yz_S_a+ABY*I_NAI_K5xyz_S_a;
  Double I_NAI_K5x2z_Py_a = I_NAI_L5xy2z_S_a+ABY*I_NAI_K5x2z_S_a;
  Double I_NAI_K4x3y_Py_a = I_NAI_L4x4y_S_a+ABY*I_NAI_K4x3y_S_a;
  Double I_NAI_K4x2yz_Py_a = I_NAI_L4x3yz_S_a+ABY*I_NAI_K4x2yz_S_a;
  Double I_NAI_K4xy2z_Py_a = I_NAI_L4x2y2z_S_a+ABY*I_NAI_K4xy2z_S_a;
  Double I_NAI_K4x3z_Py_a = I_NAI_L4xy3z_S_a+ABY*I_NAI_K4x3z_S_a;
  Double I_NAI_K3x4y_Py_a = I_NAI_L3x5y_S_a+ABY*I_NAI_K3x4y_S_a;
  Double I_NAI_K3x3yz_Py_a = I_NAI_L3x4yz_S_a+ABY*I_NAI_K3x3yz_S_a;
  Double I_NAI_K3x2y2z_Py_a = I_NAI_L3x3y2z_S_a+ABY*I_NAI_K3x2y2z_S_a;
  Double I_NAI_K3xy3z_Py_a = I_NAI_L3x2y3z_S_a+ABY*I_NAI_K3xy3z_S_a;
  Double I_NAI_K3x4z_Py_a = I_NAI_L3xy4z_S_a+ABY*I_NAI_K3x4z_S_a;
  Double I_NAI_K2x5y_Py_a = I_NAI_L2x6y_S_a+ABY*I_NAI_K2x5y_S_a;
  Double I_NAI_K2x4yz_Py_a = I_NAI_L2x5yz_S_a+ABY*I_NAI_K2x4yz_S_a;
  Double I_NAI_K2x3y2z_Py_a = I_NAI_L2x4y2z_S_a+ABY*I_NAI_K2x3y2z_S_a;
  Double I_NAI_K2x2y3z_Py_a = I_NAI_L2x3y3z_S_a+ABY*I_NAI_K2x2y3z_S_a;
  Double I_NAI_K2xy4z_Py_a = I_NAI_L2x2y4z_S_a+ABY*I_NAI_K2xy4z_S_a;
  Double I_NAI_K2x5z_Py_a = I_NAI_L2xy5z_S_a+ABY*I_NAI_K2x5z_S_a;
  Double I_NAI_Kx6y_Py_a = I_NAI_Lx7y_S_a+ABY*I_NAI_Kx6y_S_a;
  Double I_NAI_Kx5yz_Py_a = I_NAI_Lx6yz_S_a+ABY*I_NAI_Kx5yz_S_a;
  Double I_NAI_Kx4y2z_Py_a = I_NAI_Lx5y2z_S_a+ABY*I_NAI_Kx4y2z_S_a;
  Double I_NAI_Kx3y3z_Py_a = I_NAI_Lx4y3z_S_a+ABY*I_NAI_Kx3y3z_S_a;
  Double I_NAI_Kx2y4z_Py_a = I_NAI_Lx3y4z_S_a+ABY*I_NAI_Kx2y4z_S_a;
  Double I_NAI_Kxy5z_Py_a = I_NAI_Lx2y5z_S_a+ABY*I_NAI_Kxy5z_S_a;
  Double I_NAI_Kx6z_Py_a = I_NAI_Lxy6z_S_a+ABY*I_NAI_Kx6z_S_a;
  Double I_NAI_K7y_Py_a = I_NAI_L8y_S_a+ABY*I_NAI_K7y_S_a;
  Double I_NAI_K6yz_Py_a = I_NAI_L7yz_S_a+ABY*I_NAI_K6yz_S_a;
  Double I_NAI_K5y2z_Py_a = I_NAI_L6y2z_S_a+ABY*I_NAI_K5y2z_S_a;
  Double I_NAI_K4y3z_Py_a = I_NAI_L5y3z_S_a+ABY*I_NAI_K4y3z_S_a;
  Double I_NAI_K3y4z_Py_a = I_NAI_L4y4z_S_a+ABY*I_NAI_K3y4z_S_a;
  Double I_NAI_K2y5z_Py_a = I_NAI_L3y5z_S_a+ABY*I_NAI_K2y5z_S_a;
  Double I_NAI_Ky6z_Py_a = I_NAI_L2y6z_S_a+ABY*I_NAI_Ky6z_S_a;
  Double I_NAI_K7z_Py_a = I_NAI_Ly7z_S_a+ABY*I_NAI_K7z_S_a;
  Double I_NAI_K7x_Pz_a = I_NAI_L7xz_S_a+ABZ*I_NAI_K7x_S_a;
  Double I_NAI_K6xy_Pz_a = I_NAI_L6xyz_S_a+ABZ*I_NAI_K6xy_S_a;
  Double I_NAI_K6xz_Pz_a = I_NAI_L6x2z_S_a+ABZ*I_NAI_K6xz_S_a;
  Double I_NAI_K5x2y_Pz_a = I_NAI_L5x2yz_S_a+ABZ*I_NAI_K5x2y_S_a;
  Double I_NAI_K5xyz_Pz_a = I_NAI_L5xy2z_S_a+ABZ*I_NAI_K5xyz_S_a;
  Double I_NAI_K5x2z_Pz_a = I_NAI_L5x3z_S_a+ABZ*I_NAI_K5x2z_S_a;
  Double I_NAI_K4x3y_Pz_a = I_NAI_L4x3yz_S_a+ABZ*I_NAI_K4x3y_S_a;
  Double I_NAI_K4x2yz_Pz_a = I_NAI_L4x2y2z_S_a+ABZ*I_NAI_K4x2yz_S_a;
  Double I_NAI_K4xy2z_Pz_a = I_NAI_L4xy3z_S_a+ABZ*I_NAI_K4xy2z_S_a;
  Double I_NAI_K4x3z_Pz_a = I_NAI_L4x4z_S_a+ABZ*I_NAI_K4x3z_S_a;
  Double I_NAI_K3x4y_Pz_a = I_NAI_L3x4yz_S_a+ABZ*I_NAI_K3x4y_S_a;
  Double I_NAI_K3x3yz_Pz_a = I_NAI_L3x3y2z_S_a+ABZ*I_NAI_K3x3yz_S_a;
  Double I_NAI_K3x2y2z_Pz_a = I_NAI_L3x2y3z_S_a+ABZ*I_NAI_K3x2y2z_S_a;
  Double I_NAI_K3xy3z_Pz_a = I_NAI_L3xy4z_S_a+ABZ*I_NAI_K3xy3z_S_a;
  Double I_NAI_K3x4z_Pz_a = I_NAI_L3x5z_S_a+ABZ*I_NAI_K3x4z_S_a;
  Double I_NAI_K2x5y_Pz_a = I_NAI_L2x5yz_S_a+ABZ*I_NAI_K2x5y_S_a;
  Double I_NAI_K2x4yz_Pz_a = I_NAI_L2x4y2z_S_a+ABZ*I_NAI_K2x4yz_S_a;
  Double I_NAI_K2x3y2z_Pz_a = I_NAI_L2x3y3z_S_a+ABZ*I_NAI_K2x3y2z_S_a;
  Double I_NAI_K2x2y3z_Pz_a = I_NAI_L2x2y4z_S_a+ABZ*I_NAI_K2x2y3z_S_a;
  Double I_NAI_K2xy4z_Pz_a = I_NAI_L2xy5z_S_a+ABZ*I_NAI_K2xy4z_S_a;
  Double I_NAI_K2x5z_Pz_a = I_NAI_L2x6z_S_a+ABZ*I_NAI_K2x5z_S_a;
  Double I_NAI_Kx6y_Pz_a = I_NAI_Lx6yz_S_a+ABZ*I_NAI_Kx6y_S_a;
  Double I_NAI_Kx5yz_Pz_a = I_NAI_Lx5y2z_S_a+ABZ*I_NAI_Kx5yz_S_a;
  Double I_NAI_Kx4y2z_Pz_a = I_NAI_Lx4y3z_S_a+ABZ*I_NAI_Kx4y2z_S_a;
  Double I_NAI_Kx3y3z_Pz_a = I_NAI_Lx3y4z_S_a+ABZ*I_NAI_Kx3y3z_S_a;
  Double I_NAI_Kx2y4z_Pz_a = I_NAI_Lx2y5z_S_a+ABZ*I_NAI_Kx2y4z_S_a;
  Double I_NAI_Kxy5z_Pz_a = I_NAI_Lxy6z_S_a+ABZ*I_NAI_Kxy5z_S_a;
  Double I_NAI_Kx6z_Pz_a = I_NAI_Lx7z_S_a+ABZ*I_NAI_Kx6z_S_a;
  Double I_NAI_K7y_Pz_a = I_NAI_L7yz_S_a+ABZ*I_NAI_K7y_S_a;
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
   * totally 56 integrals are omitted 
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
  Double I_NAI_I6x_Dxy_a = I_NAI_K6xy_Px_a+ABY*I_NAI_I6x_Px_a;
  Double I_NAI_I5xy_Dxy_a = I_NAI_K5x2y_Px_a+ABY*I_NAI_I5xy_Px_a;
  Double I_NAI_I5xz_Dxy_a = I_NAI_K5xyz_Px_a+ABY*I_NAI_I5xz_Px_a;
  Double I_NAI_I4x2y_Dxy_a = I_NAI_K4x3y_Px_a+ABY*I_NAI_I4x2y_Px_a;
  Double I_NAI_I4xyz_Dxy_a = I_NAI_K4x2yz_Px_a+ABY*I_NAI_I4xyz_Px_a;
  Double I_NAI_I4x2z_Dxy_a = I_NAI_K4xy2z_Px_a+ABY*I_NAI_I4x2z_Px_a;
  Double I_NAI_I3x3y_Dxy_a = I_NAI_K3x4y_Px_a+ABY*I_NAI_I3x3y_Px_a;
  Double I_NAI_I3x2yz_Dxy_a = I_NAI_K3x3yz_Px_a+ABY*I_NAI_I3x2yz_Px_a;
  Double I_NAI_I3xy2z_Dxy_a = I_NAI_K3x2y2z_Px_a+ABY*I_NAI_I3xy2z_Px_a;
  Double I_NAI_I3x3z_Dxy_a = I_NAI_K3xy3z_Px_a+ABY*I_NAI_I3x3z_Px_a;
  Double I_NAI_I2x4y_Dxy_a = I_NAI_K2x5y_Px_a+ABY*I_NAI_I2x4y_Px_a;
  Double I_NAI_I2x3yz_Dxy_a = I_NAI_K2x4yz_Px_a+ABY*I_NAI_I2x3yz_Px_a;
  Double I_NAI_I2x2y2z_Dxy_a = I_NAI_K2x3y2z_Px_a+ABY*I_NAI_I2x2y2z_Px_a;
  Double I_NAI_I2xy3z_Dxy_a = I_NAI_K2x2y3z_Px_a+ABY*I_NAI_I2xy3z_Px_a;
  Double I_NAI_I2x4z_Dxy_a = I_NAI_K2xy4z_Px_a+ABY*I_NAI_I2x4z_Px_a;
  Double I_NAI_Ix5y_Dxy_a = I_NAI_Kx6y_Px_a+ABY*I_NAI_Ix5y_Px_a;
  Double I_NAI_Ix4yz_Dxy_a = I_NAI_Kx5yz_Px_a+ABY*I_NAI_Ix4yz_Px_a;
  Double I_NAI_Ix3y2z_Dxy_a = I_NAI_Kx4y2z_Px_a+ABY*I_NAI_Ix3y2z_Px_a;
  Double I_NAI_Ix2y3z_Dxy_a = I_NAI_Kx3y3z_Px_a+ABY*I_NAI_Ix2y3z_Px_a;
  Double I_NAI_Ixy4z_Dxy_a = I_NAI_Kx2y4z_Px_a+ABY*I_NAI_Ixy4z_Px_a;
  Double I_NAI_Ix5z_Dxy_a = I_NAI_Kxy5z_Px_a+ABY*I_NAI_Ix5z_Px_a;
  Double I_NAI_I6y_Dxy_a = I_NAI_K7y_Px_a+ABY*I_NAI_I6y_Px_a;
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
   * shell quartet name: SQ_NAI_L_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_S_a
   * RHS shell quartet name: SQ_NAI_L_S_a
   ************************************************************/
  Double I_NAI_L8x_Px_a = I_NAI_M9x_S_a+ABX*I_NAI_L8x_S_a;
  Double I_NAI_L7xy_Px_a = I_NAI_M8xy_S_a+ABX*I_NAI_L7xy_S_a;
  Double I_NAI_L7xz_Px_a = I_NAI_M8xz_S_a+ABX*I_NAI_L7xz_S_a;
  Double I_NAI_L6x2y_Px_a = I_NAI_M7x2y_S_a+ABX*I_NAI_L6x2y_S_a;
  Double I_NAI_L6xyz_Px_a = I_NAI_M7xyz_S_a+ABX*I_NAI_L6xyz_S_a;
  Double I_NAI_L6x2z_Px_a = I_NAI_M7x2z_S_a+ABX*I_NAI_L6x2z_S_a;
  Double I_NAI_L5x3y_Px_a = I_NAI_M6x3y_S_a+ABX*I_NAI_L5x3y_S_a;
  Double I_NAI_L5x2yz_Px_a = I_NAI_M6x2yz_S_a+ABX*I_NAI_L5x2yz_S_a;
  Double I_NAI_L5xy2z_Px_a = I_NAI_M6xy2z_S_a+ABX*I_NAI_L5xy2z_S_a;
  Double I_NAI_L5x3z_Px_a = I_NAI_M6x3z_S_a+ABX*I_NAI_L5x3z_S_a;
  Double I_NAI_L4x4y_Px_a = I_NAI_M5x4y_S_a+ABX*I_NAI_L4x4y_S_a;
  Double I_NAI_L4x3yz_Px_a = I_NAI_M5x3yz_S_a+ABX*I_NAI_L4x3yz_S_a;
  Double I_NAI_L4x2y2z_Px_a = I_NAI_M5x2y2z_S_a+ABX*I_NAI_L4x2y2z_S_a;
  Double I_NAI_L4xy3z_Px_a = I_NAI_M5xy3z_S_a+ABX*I_NAI_L4xy3z_S_a;
  Double I_NAI_L4x4z_Px_a = I_NAI_M5x4z_S_a+ABX*I_NAI_L4x4z_S_a;
  Double I_NAI_L3x5y_Px_a = I_NAI_M4x5y_S_a+ABX*I_NAI_L3x5y_S_a;
  Double I_NAI_L3x4yz_Px_a = I_NAI_M4x4yz_S_a+ABX*I_NAI_L3x4yz_S_a;
  Double I_NAI_L3x3y2z_Px_a = I_NAI_M4x3y2z_S_a+ABX*I_NAI_L3x3y2z_S_a;
  Double I_NAI_L3x2y3z_Px_a = I_NAI_M4x2y3z_S_a+ABX*I_NAI_L3x2y3z_S_a;
  Double I_NAI_L3xy4z_Px_a = I_NAI_M4xy4z_S_a+ABX*I_NAI_L3xy4z_S_a;
  Double I_NAI_L3x5z_Px_a = I_NAI_M4x5z_S_a+ABX*I_NAI_L3x5z_S_a;
  Double I_NAI_L2x6y_Px_a = I_NAI_M3x6y_S_a+ABX*I_NAI_L2x6y_S_a;
  Double I_NAI_L2x5yz_Px_a = I_NAI_M3x5yz_S_a+ABX*I_NAI_L2x5yz_S_a;
  Double I_NAI_L2x4y2z_Px_a = I_NAI_M3x4y2z_S_a+ABX*I_NAI_L2x4y2z_S_a;
  Double I_NAI_L2x3y3z_Px_a = I_NAI_M3x3y3z_S_a+ABX*I_NAI_L2x3y3z_S_a;
  Double I_NAI_L2x2y4z_Px_a = I_NAI_M3x2y4z_S_a+ABX*I_NAI_L2x2y4z_S_a;
  Double I_NAI_L2xy5z_Px_a = I_NAI_M3xy5z_S_a+ABX*I_NAI_L2xy5z_S_a;
  Double I_NAI_L2x6z_Px_a = I_NAI_M3x6z_S_a+ABX*I_NAI_L2x6z_S_a;
  Double I_NAI_Lx7y_Px_a = I_NAI_M2x7y_S_a+ABX*I_NAI_Lx7y_S_a;
  Double I_NAI_Lx6yz_Px_a = I_NAI_M2x6yz_S_a+ABX*I_NAI_Lx6yz_S_a;
  Double I_NAI_Lx5y2z_Px_a = I_NAI_M2x5y2z_S_a+ABX*I_NAI_Lx5y2z_S_a;
  Double I_NAI_Lx4y3z_Px_a = I_NAI_M2x4y3z_S_a+ABX*I_NAI_Lx4y3z_S_a;
  Double I_NAI_Lx3y4z_Px_a = I_NAI_M2x3y4z_S_a+ABX*I_NAI_Lx3y4z_S_a;
  Double I_NAI_Lx2y5z_Px_a = I_NAI_M2x2y5z_S_a+ABX*I_NAI_Lx2y5z_S_a;
  Double I_NAI_Lxy6z_Px_a = I_NAI_M2xy6z_S_a+ABX*I_NAI_Lxy6z_S_a;
  Double I_NAI_Lx7z_Px_a = I_NAI_M2x7z_S_a+ABX*I_NAI_Lx7z_S_a;
  Double I_NAI_L7yz_Px_a = I_NAI_Mx7yz_S_a+ABX*I_NAI_L7yz_S_a;
  Double I_NAI_L6y2z_Px_a = I_NAI_Mx6y2z_S_a+ABX*I_NAI_L6y2z_S_a;
  Double I_NAI_L5y3z_Px_a = I_NAI_Mx5y3z_S_a+ABX*I_NAI_L5y3z_S_a;
  Double I_NAI_L4y4z_Px_a = I_NAI_Mx4y4z_S_a+ABX*I_NAI_L4y4z_S_a;
  Double I_NAI_L3y5z_Px_a = I_NAI_Mx3y5z_S_a+ABX*I_NAI_L3y5z_S_a;
  Double I_NAI_L2y6z_Px_a = I_NAI_Mx2y6z_S_a+ABX*I_NAI_L2y6z_S_a;
  Double I_NAI_Ly7z_Px_a = I_NAI_Mxy7z_S_a+ABX*I_NAI_Ly7z_S_a;
  Double I_NAI_L7xy_Py_a = I_NAI_M7x2y_S_a+ABY*I_NAI_L7xy_S_a;
  Double I_NAI_L6x2y_Py_a = I_NAI_M6x3y_S_a+ABY*I_NAI_L6x2y_S_a;
  Double I_NAI_L6xyz_Py_a = I_NAI_M6x2yz_S_a+ABY*I_NAI_L6xyz_S_a;
  Double I_NAI_L5x3y_Py_a = I_NAI_M5x4y_S_a+ABY*I_NAI_L5x3y_S_a;
  Double I_NAI_L5x2yz_Py_a = I_NAI_M5x3yz_S_a+ABY*I_NAI_L5x2yz_S_a;
  Double I_NAI_L5xy2z_Py_a = I_NAI_M5x2y2z_S_a+ABY*I_NAI_L5xy2z_S_a;
  Double I_NAI_L4x4y_Py_a = I_NAI_M4x5y_S_a+ABY*I_NAI_L4x4y_S_a;
  Double I_NAI_L4x3yz_Py_a = I_NAI_M4x4yz_S_a+ABY*I_NAI_L4x3yz_S_a;
  Double I_NAI_L4x2y2z_Py_a = I_NAI_M4x3y2z_S_a+ABY*I_NAI_L4x2y2z_S_a;
  Double I_NAI_L4xy3z_Py_a = I_NAI_M4x2y3z_S_a+ABY*I_NAI_L4xy3z_S_a;
  Double I_NAI_L3x5y_Py_a = I_NAI_M3x6y_S_a+ABY*I_NAI_L3x5y_S_a;
  Double I_NAI_L3x4yz_Py_a = I_NAI_M3x5yz_S_a+ABY*I_NAI_L3x4yz_S_a;
  Double I_NAI_L3x3y2z_Py_a = I_NAI_M3x4y2z_S_a+ABY*I_NAI_L3x3y2z_S_a;
  Double I_NAI_L3x2y3z_Py_a = I_NAI_M3x3y3z_S_a+ABY*I_NAI_L3x2y3z_S_a;
  Double I_NAI_L3xy4z_Py_a = I_NAI_M3x2y4z_S_a+ABY*I_NAI_L3xy4z_S_a;
  Double I_NAI_L2x6y_Py_a = I_NAI_M2x7y_S_a+ABY*I_NAI_L2x6y_S_a;
  Double I_NAI_L2x5yz_Py_a = I_NAI_M2x6yz_S_a+ABY*I_NAI_L2x5yz_S_a;
  Double I_NAI_L2x4y2z_Py_a = I_NAI_M2x5y2z_S_a+ABY*I_NAI_L2x4y2z_S_a;
  Double I_NAI_L2x3y3z_Py_a = I_NAI_M2x4y3z_S_a+ABY*I_NAI_L2x3y3z_S_a;
  Double I_NAI_L2x2y4z_Py_a = I_NAI_M2x3y4z_S_a+ABY*I_NAI_L2x2y4z_S_a;
  Double I_NAI_L2xy5z_Py_a = I_NAI_M2x2y5z_S_a+ABY*I_NAI_L2xy5z_S_a;
  Double I_NAI_Lx7y_Py_a = I_NAI_Mx8y_S_a+ABY*I_NAI_Lx7y_S_a;
  Double I_NAI_Lx6yz_Py_a = I_NAI_Mx7yz_S_a+ABY*I_NAI_Lx6yz_S_a;
  Double I_NAI_Lx5y2z_Py_a = I_NAI_Mx6y2z_S_a+ABY*I_NAI_Lx5y2z_S_a;
  Double I_NAI_Lx4y3z_Py_a = I_NAI_Mx5y3z_S_a+ABY*I_NAI_Lx4y3z_S_a;
  Double I_NAI_Lx3y4z_Py_a = I_NAI_Mx4y4z_S_a+ABY*I_NAI_Lx3y4z_S_a;
  Double I_NAI_Lx2y5z_Py_a = I_NAI_Mx3y5z_S_a+ABY*I_NAI_Lx2y5z_S_a;
  Double I_NAI_Lxy6z_Py_a = I_NAI_Mx2y6z_S_a+ABY*I_NAI_Lxy6z_S_a;
  Double I_NAI_L8y_Py_a = I_NAI_M9y_S_a+ABY*I_NAI_L8y_S_a;
  Double I_NAI_L7yz_Py_a = I_NAI_M8yz_S_a+ABY*I_NAI_L7yz_S_a;
  Double I_NAI_L6y2z_Py_a = I_NAI_M7y2z_S_a+ABY*I_NAI_L6y2z_S_a;
  Double I_NAI_L5y3z_Py_a = I_NAI_M6y3z_S_a+ABY*I_NAI_L5y3z_S_a;
  Double I_NAI_L4y4z_Py_a = I_NAI_M5y4z_S_a+ABY*I_NAI_L4y4z_S_a;
  Double I_NAI_L3y5z_Py_a = I_NAI_M4y5z_S_a+ABY*I_NAI_L3y5z_S_a;
  Double I_NAI_L2y6z_Py_a = I_NAI_M3y6z_S_a+ABY*I_NAI_L2y6z_S_a;
  Double I_NAI_Ly7z_Py_a = I_NAI_M2y7z_S_a+ABY*I_NAI_Ly7z_S_a;
  Double I_NAI_L7xz_Pz_a = I_NAI_M7x2z_S_a+ABZ*I_NAI_L7xz_S_a;
  Double I_NAI_L6xyz_Pz_a = I_NAI_M6xy2z_S_a+ABZ*I_NAI_L6xyz_S_a;
  Double I_NAI_L6x2z_Pz_a = I_NAI_M6x3z_S_a+ABZ*I_NAI_L6x2z_S_a;
  Double I_NAI_L5x2yz_Pz_a = I_NAI_M5x2y2z_S_a+ABZ*I_NAI_L5x2yz_S_a;
  Double I_NAI_L5xy2z_Pz_a = I_NAI_M5xy3z_S_a+ABZ*I_NAI_L5xy2z_S_a;
  Double I_NAI_L5x3z_Pz_a = I_NAI_M5x4z_S_a+ABZ*I_NAI_L5x3z_S_a;
  Double I_NAI_L4x3yz_Pz_a = I_NAI_M4x3y2z_S_a+ABZ*I_NAI_L4x3yz_S_a;
  Double I_NAI_L4x2y2z_Pz_a = I_NAI_M4x2y3z_S_a+ABZ*I_NAI_L4x2y2z_S_a;
  Double I_NAI_L4xy3z_Pz_a = I_NAI_M4xy4z_S_a+ABZ*I_NAI_L4xy3z_S_a;
  Double I_NAI_L4x4z_Pz_a = I_NAI_M4x5z_S_a+ABZ*I_NAI_L4x4z_S_a;
  Double I_NAI_L3x4yz_Pz_a = I_NAI_M3x4y2z_S_a+ABZ*I_NAI_L3x4yz_S_a;
  Double I_NAI_L3x3y2z_Pz_a = I_NAI_M3x3y3z_S_a+ABZ*I_NAI_L3x3y2z_S_a;
  Double I_NAI_L3x2y3z_Pz_a = I_NAI_M3x2y4z_S_a+ABZ*I_NAI_L3x2y3z_S_a;
  Double I_NAI_L3xy4z_Pz_a = I_NAI_M3xy5z_S_a+ABZ*I_NAI_L3xy4z_S_a;
  Double I_NAI_L3x5z_Pz_a = I_NAI_M3x6z_S_a+ABZ*I_NAI_L3x5z_S_a;
  Double I_NAI_L2x5yz_Pz_a = I_NAI_M2x5y2z_S_a+ABZ*I_NAI_L2x5yz_S_a;
  Double I_NAI_L2x4y2z_Pz_a = I_NAI_M2x4y3z_S_a+ABZ*I_NAI_L2x4y2z_S_a;
  Double I_NAI_L2x3y3z_Pz_a = I_NAI_M2x3y4z_S_a+ABZ*I_NAI_L2x3y3z_S_a;
  Double I_NAI_L2x2y4z_Pz_a = I_NAI_M2x2y5z_S_a+ABZ*I_NAI_L2x2y4z_S_a;
  Double I_NAI_L2xy5z_Pz_a = I_NAI_M2xy6z_S_a+ABZ*I_NAI_L2xy5z_S_a;
  Double I_NAI_L2x6z_Pz_a = I_NAI_M2x7z_S_a+ABZ*I_NAI_L2x6z_S_a;
  Double I_NAI_Lx6yz_Pz_a = I_NAI_Mx6y2z_S_a+ABZ*I_NAI_Lx6yz_S_a;
  Double I_NAI_Lx5y2z_Pz_a = I_NAI_Mx5y3z_S_a+ABZ*I_NAI_Lx5y2z_S_a;
  Double I_NAI_Lx4y3z_Pz_a = I_NAI_Mx4y4z_S_a+ABZ*I_NAI_Lx4y3z_S_a;
  Double I_NAI_Lx3y4z_Pz_a = I_NAI_Mx3y5z_S_a+ABZ*I_NAI_Lx3y4z_S_a;
  Double I_NAI_Lx2y5z_Pz_a = I_NAI_Mx2y6z_S_a+ABZ*I_NAI_Lx2y5z_S_a;
  Double I_NAI_Lxy6z_Pz_a = I_NAI_Mxy7z_S_a+ABZ*I_NAI_Lxy6z_S_a;
  Double I_NAI_Lx7z_Pz_a = I_NAI_Mx8z_S_a+ABZ*I_NAI_Lx7z_S_a;
  Double I_NAI_L7yz_Pz_a = I_NAI_M7y2z_S_a+ABZ*I_NAI_L7yz_S_a;
  Double I_NAI_L6y2z_Pz_a = I_NAI_M6y3z_S_a+ABZ*I_NAI_L6y2z_S_a;
  Double I_NAI_L5y3z_Pz_a = I_NAI_M5y4z_S_a+ABZ*I_NAI_L5y3z_S_a;
  Double I_NAI_L4y4z_Pz_a = I_NAI_M4y5z_S_a+ABZ*I_NAI_L4y4z_S_a;
  Double I_NAI_L3y5z_Pz_a = I_NAI_M3y6z_S_a+ABZ*I_NAI_L3y5z_S_a;
  Double I_NAI_L2y6z_Pz_a = I_NAI_M2y7z_S_a+ABZ*I_NAI_L2y6z_S_a;
  Double I_NAI_Ly7z_Pz_a = I_NAI_My8z_S_a+ABZ*I_NAI_Ly7z_S_a;
  Double I_NAI_L8z_Pz_a = I_NAI_M9z_S_a+ABZ*I_NAI_L8z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_K_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 80 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_P_a
   * RHS shell quartet name: SQ_NAI_K_P_a
   ************************************************************/
  Double I_NAI_K7x_D2x_a = I_NAI_L8x_Px_a+ABX*I_NAI_K7x_Px_a;
  Double I_NAI_K6xy_D2x_a = I_NAI_L7xy_Px_a+ABX*I_NAI_K6xy_Px_a;
  Double I_NAI_K6xz_D2x_a = I_NAI_L7xz_Px_a+ABX*I_NAI_K6xz_Px_a;
  Double I_NAI_K5x2y_D2x_a = I_NAI_L6x2y_Px_a+ABX*I_NAI_K5x2y_Px_a;
  Double I_NAI_K5xyz_D2x_a = I_NAI_L6xyz_Px_a+ABX*I_NAI_K5xyz_Px_a;
  Double I_NAI_K5x2z_D2x_a = I_NAI_L6x2z_Px_a+ABX*I_NAI_K5x2z_Px_a;
  Double I_NAI_K4x3y_D2x_a = I_NAI_L5x3y_Px_a+ABX*I_NAI_K4x3y_Px_a;
  Double I_NAI_K4x2yz_D2x_a = I_NAI_L5x2yz_Px_a+ABX*I_NAI_K4x2yz_Px_a;
  Double I_NAI_K4xy2z_D2x_a = I_NAI_L5xy2z_Px_a+ABX*I_NAI_K4xy2z_Px_a;
  Double I_NAI_K4x3z_D2x_a = I_NAI_L5x3z_Px_a+ABX*I_NAI_K4x3z_Px_a;
  Double I_NAI_K3x4y_D2x_a = I_NAI_L4x4y_Px_a+ABX*I_NAI_K3x4y_Px_a;
  Double I_NAI_K3x3yz_D2x_a = I_NAI_L4x3yz_Px_a+ABX*I_NAI_K3x3yz_Px_a;
  Double I_NAI_K3x2y2z_D2x_a = I_NAI_L4x2y2z_Px_a+ABX*I_NAI_K3x2y2z_Px_a;
  Double I_NAI_K3xy3z_D2x_a = I_NAI_L4xy3z_Px_a+ABX*I_NAI_K3xy3z_Px_a;
  Double I_NAI_K3x4z_D2x_a = I_NAI_L4x4z_Px_a+ABX*I_NAI_K3x4z_Px_a;
  Double I_NAI_K2x5y_D2x_a = I_NAI_L3x5y_Px_a+ABX*I_NAI_K2x5y_Px_a;
  Double I_NAI_K2x4yz_D2x_a = I_NAI_L3x4yz_Px_a+ABX*I_NAI_K2x4yz_Px_a;
  Double I_NAI_K2x3y2z_D2x_a = I_NAI_L3x3y2z_Px_a+ABX*I_NAI_K2x3y2z_Px_a;
  Double I_NAI_K2x2y3z_D2x_a = I_NAI_L3x2y3z_Px_a+ABX*I_NAI_K2x2y3z_Px_a;
  Double I_NAI_K2xy4z_D2x_a = I_NAI_L3xy4z_Px_a+ABX*I_NAI_K2xy4z_Px_a;
  Double I_NAI_K2x5z_D2x_a = I_NAI_L3x5z_Px_a+ABX*I_NAI_K2x5z_Px_a;
  Double I_NAI_Kx6y_D2x_a = I_NAI_L2x6y_Px_a+ABX*I_NAI_Kx6y_Px_a;
  Double I_NAI_Kx5yz_D2x_a = I_NAI_L2x5yz_Px_a+ABX*I_NAI_Kx5yz_Px_a;
  Double I_NAI_Kx4y2z_D2x_a = I_NAI_L2x4y2z_Px_a+ABX*I_NAI_Kx4y2z_Px_a;
  Double I_NAI_Kx3y3z_D2x_a = I_NAI_L2x3y3z_Px_a+ABX*I_NAI_Kx3y3z_Px_a;
  Double I_NAI_Kx2y4z_D2x_a = I_NAI_L2x2y4z_Px_a+ABX*I_NAI_Kx2y4z_Px_a;
  Double I_NAI_Kxy5z_D2x_a = I_NAI_L2xy5z_Px_a+ABX*I_NAI_Kxy5z_Px_a;
  Double I_NAI_Kx6z_D2x_a = I_NAI_L2x6z_Px_a+ABX*I_NAI_Kx6z_Px_a;
  Double I_NAI_K7y_D2x_a = I_NAI_Lx7y_Px_a+ABX*I_NAI_K7y_Px_a;
  Double I_NAI_K6yz_D2x_a = I_NAI_Lx6yz_Px_a+ABX*I_NAI_K6yz_Px_a;
  Double I_NAI_K5y2z_D2x_a = I_NAI_Lx5y2z_Px_a+ABX*I_NAI_K5y2z_Px_a;
  Double I_NAI_K4y3z_D2x_a = I_NAI_Lx4y3z_Px_a+ABX*I_NAI_K4y3z_Px_a;
  Double I_NAI_K3y4z_D2x_a = I_NAI_Lx3y4z_Px_a+ABX*I_NAI_K3y4z_Px_a;
  Double I_NAI_K2y5z_D2x_a = I_NAI_Lx2y5z_Px_a+ABX*I_NAI_K2y5z_Px_a;
  Double I_NAI_Ky6z_D2x_a = I_NAI_Lxy6z_Px_a+ABX*I_NAI_Ky6z_Px_a;
  Double I_NAI_K7z_D2x_a = I_NAI_Lx7z_Px_a+ABX*I_NAI_K7z_Px_a;
  Double I_NAI_K6xz_Dxy_a = I_NAI_L6xyz_Px_a+ABY*I_NAI_K6xz_Px_a;
  Double I_NAI_K5xyz_Dxy_a = I_NAI_L5x2yz_Px_a+ABY*I_NAI_K5xyz_Px_a;
  Double I_NAI_K5x2z_Dxy_a = I_NAI_L5xy2z_Px_a+ABY*I_NAI_K5x2z_Px_a;
  Double I_NAI_K4x2yz_Dxy_a = I_NAI_L4x3yz_Px_a+ABY*I_NAI_K4x2yz_Px_a;
  Double I_NAI_K4xy2z_Dxy_a = I_NAI_L4x2y2z_Px_a+ABY*I_NAI_K4xy2z_Px_a;
  Double I_NAI_K4x3z_Dxy_a = I_NAI_L4xy3z_Px_a+ABY*I_NAI_K4x3z_Px_a;
  Double I_NAI_K3x3yz_Dxy_a = I_NAI_L3x4yz_Px_a+ABY*I_NAI_K3x3yz_Px_a;
  Double I_NAI_K3x2y2z_Dxy_a = I_NAI_L3x3y2z_Px_a+ABY*I_NAI_K3x2y2z_Px_a;
  Double I_NAI_K3xy3z_Dxy_a = I_NAI_L3x2y3z_Px_a+ABY*I_NAI_K3xy3z_Px_a;
  Double I_NAI_K3x4z_Dxy_a = I_NAI_L3xy4z_Px_a+ABY*I_NAI_K3x4z_Px_a;
  Double I_NAI_K2x4yz_Dxy_a = I_NAI_L2x5yz_Px_a+ABY*I_NAI_K2x4yz_Px_a;
  Double I_NAI_K2x3y2z_Dxy_a = I_NAI_L2x4y2z_Px_a+ABY*I_NAI_K2x3y2z_Px_a;
  Double I_NAI_K2x2y3z_Dxy_a = I_NAI_L2x3y3z_Px_a+ABY*I_NAI_K2x2y3z_Px_a;
  Double I_NAI_K2xy4z_Dxy_a = I_NAI_L2x2y4z_Px_a+ABY*I_NAI_K2xy4z_Px_a;
  Double I_NAI_K2x5z_Dxy_a = I_NAI_L2xy5z_Px_a+ABY*I_NAI_K2x5z_Px_a;
  Double I_NAI_Kx5yz_Dxy_a = I_NAI_Lx6yz_Px_a+ABY*I_NAI_Kx5yz_Px_a;
  Double I_NAI_Kx4y2z_Dxy_a = I_NAI_Lx5y2z_Px_a+ABY*I_NAI_Kx4y2z_Px_a;
  Double I_NAI_Kx3y3z_Dxy_a = I_NAI_Lx4y3z_Px_a+ABY*I_NAI_Kx3y3z_Px_a;
  Double I_NAI_Kx2y4z_Dxy_a = I_NAI_Lx3y4z_Px_a+ABY*I_NAI_Kx2y4z_Px_a;
  Double I_NAI_Kxy5z_Dxy_a = I_NAI_Lx2y5z_Px_a+ABY*I_NAI_Kxy5z_Px_a;
  Double I_NAI_Kx6z_Dxy_a = I_NAI_Lxy6z_Px_a+ABY*I_NAI_Kx6z_Px_a;
  Double I_NAI_K6yz_Dxy_a = I_NAI_L7yz_Px_a+ABY*I_NAI_K6yz_Px_a;
  Double I_NAI_K5y2z_Dxy_a = I_NAI_L6y2z_Px_a+ABY*I_NAI_K5y2z_Px_a;
  Double I_NAI_K4y3z_Dxy_a = I_NAI_L5y3z_Px_a+ABY*I_NAI_K4y3z_Px_a;
  Double I_NAI_K3y4z_Dxy_a = I_NAI_L4y4z_Px_a+ABY*I_NAI_K3y4z_Px_a;
  Double I_NAI_K2y5z_Dxy_a = I_NAI_L3y5z_Px_a+ABY*I_NAI_K2y5z_Px_a;
  Double I_NAI_Ky6z_Dxy_a = I_NAI_L2y6z_Px_a+ABY*I_NAI_Ky6z_Px_a;
  Double I_NAI_K7z_Dxy_a = I_NAI_Ly7z_Px_a+ABY*I_NAI_K7z_Px_a;
  Double I_NAI_K7x_D2y_a = I_NAI_L7xy_Py_a+ABY*I_NAI_K7x_Py_a;
  Double I_NAI_K6xy_D2y_a = I_NAI_L6x2y_Py_a+ABY*I_NAI_K6xy_Py_a;
  Double I_NAI_K6xz_D2y_a = I_NAI_L6xyz_Py_a+ABY*I_NAI_K6xz_Py_a;
  Double I_NAI_K5x2y_D2y_a = I_NAI_L5x3y_Py_a+ABY*I_NAI_K5x2y_Py_a;
  Double I_NAI_K5xyz_D2y_a = I_NAI_L5x2yz_Py_a+ABY*I_NAI_K5xyz_Py_a;
  Double I_NAI_K5x2z_D2y_a = I_NAI_L5xy2z_Py_a+ABY*I_NAI_K5x2z_Py_a;
  Double I_NAI_K4x3y_D2y_a = I_NAI_L4x4y_Py_a+ABY*I_NAI_K4x3y_Py_a;
  Double I_NAI_K4x2yz_D2y_a = I_NAI_L4x3yz_Py_a+ABY*I_NAI_K4x2yz_Py_a;
  Double I_NAI_K4xy2z_D2y_a = I_NAI_L4x2y2z_Py_a+ABY*I_NAI_K4xy2z_Py_a;
  Double I_NAI_K4x3z_D2y_a = I_NAI_L4xy3z_Py_a+ABY*I_NAI_K4x3z_Py_a;
  Double I_NAI_K3x4y_D2y_a = I_NAI_L3x5y_Py_a+ABY*I_NAI_K3x4y_Py_a;
  Double I_NAI_K3x3yz_D2y_a = I_NAI_L3x4yz_Py_a+ABY*I_NAI_K3x3yz_Py_a;
  Double I_NAI_K3x2y2z_D2y_a = I_NAI_L3x3y2z_Py_a+ABY*I_NAI_K3x2y2z_Py_a;
  Double I_NAI_K3xy3z_D2y_a = I_NAI_L3x2y3z_Py_a+ABY*I_NAI_K3xy3z_Py_a;
  Double I_NAI_K3x4z_D2y_a = I_NAI_L3xy4z_Py_a+ABY*I_NAI_K3x4z_Py_a;
  Double I_NAI_K2x5y_D2y_a = I_NAI_L2x6y_Py_a+ABY*I_NAI_K2x5y_Py_a;
  Double I_NAI_K2x4yz_D2y_a = I_NAI_L2x5yz_Py_a+ABY*I_NAI_K2x4yz_Py_a;
  Double I_NAI_K2x3y2z_D2y_a = I_NAI_L2x4y2z_Py_a+ABY*I_NAI_K2x3y2z_Py_a;
  Double I_NAI_K2x2y3z_D2y_a = I_NAI_L2x3y3z_Py_a+ABY*I_NAI_K2x2y3z_Py_a;
  Double I_NAI_K2xy4z_D2y_a = I_NAI_L2x2y4z_Py_a+ABY*I_NAI_K2xy4z_Py_a;
  Double I_NAI_K2x5z_D2y_a = I_NAI_L2xy5z_Py_a+ABY*I_NAI_K2x5z_Py_a;
  Double I_NAI_Kx6y_D2y_a = I_NAI_Lx7y_Py_a+ABY*I_NAI_Kx6y_Py_a;
  Double I_NAI_Kx5yz_D2y_a = I_NAI_Lx6yz_Py_a+ABY*I_NAI_Kx5yz_Py_a;
  Double I_NAI_Kx4y2z_D2y_a = I_NAI_Lx5y2z_Py_a+ABY*I_NAI_Kx4y2z_Py_a;
  Double I_NAI_Kx3y3z_D2y_a = I_NAI_Lx4y3z_Py_a+ABY*I_NAI_Kx3y3z_Py_a;
  Double I_NAI_Kx2y4z_D2y_a = I_NAI_Lx3y4z_Py_a+ABY*I_NAI_Kx2y4z_Py_a;
  Double I_NAI_Kxy5z_D2y_a = I_NAI_Lx2y5z_Py_a+ABY*I_NAI_Kxy5z_Py_a;
  Double I_NAI_Kx6z_D2y_a = I_NAI_Lxy6z_Py_a+ABY*I_NAI_Kx6z_Py_a;
  Double I_NAI_K7y_D2y_a = I_NAI_L8y_Py_a+ABY*I_NAI_K7y_Py_a;
  Double I_NAI_K6yz_D2y_a = I_NAI_L7yz_Py_a+ABY*I_NAI_K6yz_Py_a;
  Double I_NAI_K5y2z_D2y_a = I_NAI_L6y2z_Py_a+ABY*I_NAI_K5y2z_Py_a;
  Double I_NAI_K4y3z_D2y_a = I_NAI_L5y3z_Py_a+ABY*I_NAI_K4y3z_Py_a;
  Double I_NAI_K3y4z_D2y_a = I_NAI_L4y4z_Py_a+ABY*I_NAI_K3y4z_Py_a;
  Double I_NAI_K2y5z_D2y_a = I_NAI_L3y5z_Py_a+ABY*I_NAI_K2y5z_Py_a;
  Double I_NAI_Ky6z_D2y_a = I_NAI_L2y6z_Py_a+ABY*I_NAI_Ky6z_Py_a;
  Double I_NAI_K7z_D2y_a = I_NAI_Ly7z_Py_a+ABY*I_NAI_K7z_Py_a;
  Double I_NAI_K7x_D2z_a = I_NAI_L7xz_Pz_a+ABZ*I_NAI_K7x_Pz_a;
  Double I_NAI_K6xy_D2z_a = I_NAI_L6xyz_Pz_a+ABZ*I_NAI_K6xy_Pz_a;
  Double I_NAI_K6xz_D2z_a = I_NAI_L6x2z_Pz_a+ABZ*I_NAI_K6xz_Pz_a;
  Double I_NAI_K5x2y_D2z_a = I_NAI_L5x2yz_Pz_a+ABZ*I_NAI_K5x2y_Pz_a;
  Double I_NAI_K5xyz_D2z_a = I_NAI_L5xy2z_Pz_a+ABZ*I_NAI_K5xyz_Pz_a;
  Double I_NAI_K5x2z_D2z_a = I_NAI_L5x3z_Pz_a+ABZ*I_NAI_K5x2z_Pz_a;
  Double I_NAI_K4x3y_D2z_a = I_NAI_L4x3yz_Pz_a+ABZ*I_NAI_K4x3y_Pz_a;
  Double I_NAI_K4x2yz_D2z_a = I_NAI_L4x2y2z_Pz_a+ABZ*I_NAI_K4x2yz_Pz_a;
  Double I_NAI_K4xy2z_D2z_a = I_NAI_L4xy3z_Pz_a+ABZ*I_NAI_K4xy2z_Pz_a;
  Double I_NAI_K4x3z_D2z_a = I_NAI_L4x4z_Pz_a+ABZ*I_NAI_K4x3z_Pz_a;
  Double I_NAI_K3x4y_D2z_a = I_NAI_L3x4yz_Pz_a+ABZ*I_NAI_K3x4y_Pz_a;
  Double I_NAI_K3x3yz_D2z_a = I_NAI_L3x3y2z_Pz_a+ABZ*I_NAI_K3x3yz_Pz_a;
  Double I_NAI_K3x2y2z_D2z_a = I_NAI_L3x2y3z_Pz_a+ABZ*I_NAI_K3x2y2z_Pz_a;
  Double I_NAI_K3xy3z_D2z_a = I_NAI_L3xy4z_Pz_a+ABZ*I_NAI_K3xy3z_Pz_a;
  Double I_NAI_K3x4z_D2z_a = I_NAI_L3x5z_Pz_a+ABZ*I_NAI_K3x4z_Pz_a;
  Double I_NAI_K2x5y_D2z_a = I_NAI_L2x5yz_Pz_a+ABZ*I_NAI_K2x5y_Pz_a;
  Double I_NAI_K2x4yz_D2z_a = I_NAI_L2x4y2z_Pz_a+ABZ*I_NAI_K2x4yz_Pz_a;
  Double I_NAI_K2x3y2z_D2z_a = I_NAI_L2x3y3z_Pz_a+ABZ*I_NAI_K2x3y2z_Pz_a;
  Double I_NAI_K2x2y3z_D2z_a = I_NAI_L2x2y4z_Pz_a+ABZ*I_NAI_K2x2y3z_Pz_a;
  Double I_NAI_K2xy4z_D2z_a = I_NAI_L2xy5z_Pz_a+ABZ*I_NAI_K2xy4z_Pz_a;
  Double I_NAI_K2x5z_D2z_a = I_NAI_L2x6z_Pz_a+ABZ*I_NAI_K2x5z_Pz_a;
  Double I_NAI_Kx6y_D2z_a = I_NAI_Lx6yz_Pz_a+ABZ*I_NAI_Kx6y_Pz_a;
  Double I_NAI_Kx5yz_D2z_a = I_NAI_Lx5y2z_Pz_a+ABZ*I_NAI_Kx5yz_Pz_a;
  Double I_NAI_Kx4y2z_D2z_a = I_NAI_Lx4y3z_Pz_a+ABZ*I_NAI_Kx4y2z_Pz_a;
  Double I_NAI_Kx3y3z_D2z_a = I_NAI_Lx3y4z_Pz_a+ABZ*I_NAI_Kx3y3z_Pz_a;
  Double I_NAI_Kx2y4z_D2z_a = I_NAI_Lx2y5z_Pz_a+ABZ*I_NAI_Kx2y4z_Pz_a;
  Double I_NAI_Kxy5z_D2z_a = I_NAI_Lxy6z_Pz_a+ABZ*I_NAI_Kxy5z_Pz_a;
  Double I_NAI_Kx6z_D2z_a = I_NAI_Lx7z_Pz_a+ABZ*I_NAI_Kx6z_Pz_a;
  Double I_NAI_K7y_D2z_a = I_NAI_L7yz_Pz_a+ABZ*I_NAI_K7y_Pz_a;
  Double I_NAI_K6yz_D2z_a = I_NAI_L6y2z_Pz_a+ABZ*I_NAI_K6yz_Pz_a;
  Double I_NAI_K5y2z_D2z_a = I_NAI_L5y3z_Pz_a+ABZ*I_NAI_K5y2z_Pz_a;
  Double I_NAI_K4y3z_D2z_a = I_NAI_L4y4z_Pz_a+ABZ*I_NAI_K4y3z_Pz_a;
  Double I_NAI_K3y4z_D2z_a = I_NAI_L3y5z_Pz_a+ABZ*I_NAI_K3y4z_Pz_a;
  Double I_NAI_K2y5z_D2z_a = I_NAI_L2y6z_Pz_a+ABZ*I_NAI_K2y5z_Pz_a;
  Double I_NAI_Ky6z_D2z_a = I_NAI_Ly7z_Pz_a+ABZ*I_NAI_Ky6z_Pz_a;
  Double I_NAI_K7z_D2z_a = I_NAI_L8z_Pz_a+ABZ*I_NAI_K7z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_I_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_D_a
   * RHS shell quartet name: SQ_NAI_I_D_a
   ************************************************************/
  Double I_NAI_I6x_F3x_a = I_NAI_K7x_D2x_a+ABX*I_NAI_I6x_D2x_a;
  Double I_NAI_I5xy_F3x_a = I_NAI_K6xy_D2x_a+ABX*I_NAI_I5xy_D2x_a;
  Double I_NAI_I5xz_F3x_a = I_NAI_K6xz_D2x_a+ABX*I_NAI_I5xz_D2x_a;
  Double I_NAI_I4x2y_F3x_a = I_NAI_K5x2y_D2x_a+ABX*I_NAI_I4x2y_D2x_a;
  Double I_NAI_I4xyz_F3x_a = I_NAI_K5xyz_D2x_a+ABX*I_NAI_I4xyz_D2x_a;
  Double I_NAI_I4x2z_F3x_a = I_NAI_K5x2z_D2x_a+ABX*I_NAI_I4x2z_D2x_a;
  Double I_NAI_I3x3y_F3x_a = I_NAI_K4x3y_D2x_a+ABX*I_NAI_I3x3y_D2x_a;
  Double I_NAI_I3x2yz_F3x_a = I_NAI_K4x2yz_D2x_a+ABX*I_NAI_I3x2yz_D2x_a;
  Double I_NAI_I3xy2z_F3x_a = I_NAI_K4xy2z_D2x_a+ABX*I_NAI_I3xy2z_D2x_a;
  Double I_NAI_I3x3z_F3x_a = I_NAI_K4x3z_D2x_a+ABX*I_NAI_I3x3z_D2x_a;
  Double I_NAI_I2x4y_F3x_a = I_NAI_K3x4y_D2x_a+ABX*I_NAI_I2x4y_D2x_a;
  Double I_NAI_I2x3yz_F3x_a = I_NAI_K3x3yz_D2x_a+ABX*I_NAI_I2x3yz_D2x_a;
  Double I_NAI_I2x2y2z_F3x_a = I_NAI_K3x2y2z_D2x_a+ABX*I_NAI_I2x2y2z_D2x_a;
  Double I_NAI_I2xy3z_F3x_a = I_NAI_K3xy3z_D2x_a+ABX*I_NAI_I2xy3z_D2x_a;
  Double I_NAI_I2x4z_F3x_a = I_NAI_K3x4z_D2x_a+ABX*I_NAI_I2x4z_D2x_a;
  Double I_NAI_Ix5y_F3x_a = I_NAI_K2x5y_D2x_a+ABX*I_NAI_Ix5y_D2x_a;
  Double I_NAI_Ix4yz_F3x_a = I_NAI_K2x4yz_D2x_a+ABX*I_NAI_Ix4yz_D2x_a;
  Double I_NAI_Ix3y2z_F3x_a = I_NAI_K2x3y2z_D2x_a+ABX*I_NAI_Ix3y2z_D2x_a;
  Double I_NAI_Ix2y3z_F3x_a = I_NAI_K2x2y3z_D2x_a+ABX*I_NAI_Ix2y3z_D2x_a;
  Double I_NAI_Ixy4z_F3x_a = I_NAI_K2xy4z_D2x_a+ABX*I_NAI_Ixy4z_D2x_a;
  Double I_NAI_Ix5z_F3x_a = I_NAI_K2x5z_D2x_a+ABX*I_NAI_Ix5z_D2x_a;
  Double I_NAI_I6y_F3x_a = I_NAI_Kx6y_D2x_a+ABX*I_NAI_I6y_D2x_a;
  Double I_NAI_I5yz_F3x_a = I_NAI_Kx5yz_D2x_a+ABX*I_NAI_I5yz_D2x_a;
  Double I_NAI_I4y2z_F3x_a = I_NAI_Kx4y2z_D2x_a+ABX*I_NAI_I4y2z_D2x_a;
  Double I_NAI_I3y3z_F3x_a = I_NAI_Kx3y3z_D2x_a+ABX*I_NAI_I3y3z_D2x_a;
  Double I_NAI_I2y4z_F3x_a = I_NAI_Kx2y4z_D2x_a+ABX*I_NAI_I2y4z_D2x_a;
  Double I_NAI_Iy5z_F3x_a = I_NAI_Kxy5z_D2x_a+ABX*I_NAI_Iy5z_D2x_a;
  Double I_NAI_I6z_F3x_a = I_NAI_Kx6z_D2x_a+ABX*I_NAI_I6z_D2x_a;
  Double I_NAI_I6x_F2xy_a = I_NAI_K6xy_D2x_a+ABY*I_NAI_I6x_D2x_a;
  Double I_NAI_I5xy_F2xy_a = I_NAI_K5x2y_D2x_a+ABY*I_NAI_I5xy_D2x_a;
  Double I_NAI_I5xz_F2xy_a = I_NAI_K5xyz_D2x_a+ABY*I_NAI_I5xz_D2x_a;
  Double I_NAI_I4x2y_F2xy_a = I_NAI_K4x3y_D2x_a+ABY*I_NAI_I4x2y_D2x_a;
  Double I_NAI_I4xyz_F2xy_a = I_NAI_K4x2yz_D2x_a+ABY*I_NAI_I4xyz_D2x_a;
  Double I_NAI_I4x2z_F2xy_a = I_NAI_K4xy2z_D2x_a+ABY*I_NAI_I4x2z_D2x_a;
  Double I_NAI_I3x3y_F2xy_a = I_NAI_K3x4y_D2x_a+ABY*I_NAI_I3x3y_D2x_a;
  Double I_NAI_I3x2yz_F2xy_a = I_NAI_K3x3yz_D2x_a+ABY*I_NAI_I3x2yz_D2x_a;
  Double I_NAI_I3xy2z_F2xy_a = I_NAI_K3x2y2z_D2x_a+ABY*I_NAI_I3xy2z_D2x_a;
  Double I_NAI_I3x3z_F2xy_a = I_NAI_K3xy3z_D2x_a+ABY*I_NAI_I3x3z_D2x_a;
  Double I_NAI_I2x4y_F2xy_a = I_NAI_K2x5y_D2x_a+ABY*I_NAI_I2x4y_D2x_a;
  Double I_NAI_I2x3yz_F2xy_a = I_NAI_K2x4yz_D2x_a+ABY*I_NAI_I2x3yz_D2x_a;
  Double I_NAI_I2x2y2z_F2xy_a = I_NAI_K2x3y2z_D2x_a+ABY*I_NAI_I2x2y2z_D2x_a;
  Double I_NAI_I2xy3z_F2xy_a = I_NAI_K2x2y3z_D2x_a+ABY*I_NAI_I2xy3z_D2x_a;
  Double I_NAI_I2x4z_F2xy_a = I_NAI_K2xy4z_D2x_a+ABY*I_NAI_I2x4z_D2x_a;
  Double I_NAI_Ix5y_F2xy_a = I_NAI_Kx6y_D2x_a+ABY*I_NAI_Ix5y_D2x_a;
  Double I_NAI_Ix4yz_F2xy_a = I_NAI_Kx5yz_D2x_a+ABY*I_NAI_Ix4yz_D2x_a;
  Double I_NAI_Ix3y2z_F2xy_a = I_NAI_Kx4y2z_D2x_a+ABY*I_NAI_Ix3y2z_D2x_a;
  Double I_NAI_Ix2y3z_F2xy_a = I_NAI_Kx3y3z_D2x_a+ABY*I_NAI_Ix2y3z_D2x_a;
  Double I_NAI_Ixy4z_F2xy_a = I_NAI_Kx2y4z_D2x_a+ABY*I_NAI_Ixy4z_D2x_a;
  Double I_NAI_Ix5z_F2xy_a = I_NAI_Kxy5z_D2x_a+ABY*I_NAI_Ix5z_D2x_a;
  Double I_NAI_I6y_F2xy_a = I_NAI_K7y_D2x_a+ABY*I_NAI_I6y_D2x_a;
  Double I_NAI_I5yz_F2xy_a = I_NAI_K6yz_D2x_a+ABY*I_NAI_I5yz_D2x_a;
  Double I_NAI_I4y2z_F2xy_a = I_NAI_K5y2z_D2x_a+ABY*I_NAI_I4y2z_D2x_a;
  Double I_NAI_I3y3z_F2xy_a = I_NAI_K4y3z_D2x_a+ABY*I_NAI_I3y3z_D2x_a;
  Double I_NAI_I2y4z_F2xy_a = I_NAI_K3y4z_D2x_a+ABY*I_NAI_I2y4z_D2x_a;
  Double I_NAI_Iy5z_F2xy_a = I_NAI_K2y5z_D2x_a+ABY*I_NAI_Iy5z_D2x_a;
  Double I_NAI_I6z_F2xy_a = I_NAI_Ky6z_D2x_a+ABY*I_NAI_I6z_D2x_a;
  Double I_NAI_I6x_F2xz_a = I_NAI_K6xz_D2x_a+ABZ*I_NAI_I6x_D2x_a;
  Double I_NAI_I5xy_F2xz_a = I_NAI_K5xyz_D2x_a+ABZ*I_NAI_I5xy_D2x_a;
  Double I_NAI_I5xz_F2xz_a = I_NAI_K5x2z_D2x_a+ABZ*I_NAI_I5xz_D2x_a;
  Double I_NAI_I4x2y_F2xz_a = I_NAI_K4x2yz_D2x_a+ABZ*I_NAI_I4x2y_D2x_a;
  Double I_NAI_I4xyz_F2xz_a = I_NAI_K4xy2z_D2x_a+ABZ*I_NAI_I4xyz_D2x_a;
  Double I_NAI_I4x2z_F2xz_a = I_NAI_K4x3z_D2x_a+ABZ*I_NAI_I4x2z_D2x_a;
  Double I_NAI_I3x3y_F2xz_a = I_NAI_K3x3yz_D2x_a+ABZ*I_NAI_I3x3y_D2x_a;
  Double I_NAI_I3x2yz_F2xz_a = I_NAI_K3x2y2z_D2x_a+ABZ*I_NAI_I3x2yz_D2x_a;
  Double I_NAI_I3xy2z_F2xz_a = I_NAI_K3xy3z_D2x_a+ABZ*I_NAI_I3xy2z_D2x_a;
  Double I_NAI_I3x3z_F2xz_a = I_NAI_K3x4z_D2x_a+ABZ*I_NAI_I3x3z_D2x_a;
  Double I_NAI_I2x4y_F2xz_a = I_NAI_K2x4yz_D2x_a+ABZ*I_NAI_I2x4y_D2x_a;
  Double I_NAI_I2x3yz_F2xz_a = I_NAI_K2x3y2z_D2x_a+ABZ*I_NAI_I2x3yz_D2x_a;
  Double I_NAI_I2x2y2z_F2xz_a = I_NAI_K2x2y3z_D2x_a+ABZ*I_NAI_I2x2y2z_D2x_a;
  Double I_NAI_I2xy3z_F2xz_a = I_NAI_K2xy4z_D2x_a+ABZ*I_NAI_I2xy3z_D2x_a;
  Double I_NAI_I2x4z_F2xz_a = I_NAI_K2x5z_D2x_a+ABZ*I_NAI_I2x4z_D2x_a;
  Double I_NAI_Ix5y_F2xz_a = I_NAI_Kx5yz_D2x_a+ABZ*I_NAI_Ix5y_D2x_a;
  Double I_NAI_Ix4yz_F2xz_a = I_NAI_Kx4y2z_D2x_a+ABZ*I_NAI_Ix4yz_D2x_a;
  Double I_NAI_Ix3y2z_F2xz_a = I_NAI_Kx3y3z_D2x_a+ABZ*I_NAI_Ix3y2z_D2x_a;
  Double I_NAI_Ix2y3z_F2xz_a = I_NAI_Kx2y4z_D2x_a+ABZ*I_NAI_Ix2y3z_D2x_a;
  Double I_NAI_Ixy4z_F2xz_a = I_NAI_Kxy5z_D2x_a+ABZ*I_NAI_Ixy4z_D2x_a;
  Double I_NAI_Ix5z_F2xz_a = I_NAI_Kx6z_D2x_a+ABZ*I_NAI_Ix5z_D2x_a;
  Double I_NAI_I6y_F2xz_a = I_NAI_K6yz_D2x_a+ABZ*I_NAI_I6y_D2x_a;
  Double I_NAI_I5yz_F2xz_a = I_NAI_K5y2z_D2x_a+ABZ*I_NAI_I5yz_D2x_a;
  Double I_NAI_I4y2z_F2xz_a = I_NAI_K4y3z_D2x_a+ABZ*I_NAI_I4y2z_D2x_a;
  Double I_NAI_I3y3z_F2xz_a = I_NAI_K3y4z_D2x_a+ABZ*I_NAI_I3y3z_D2x_a;
  Double I_NAI_I2y4z_F2xz_a = I_NAI_K2y5z_D2x_a+ABZ*I_NAI_I2y4z_D2x_a;
  Double I_NAI_Iy5z_F2xz_a = I_NAI_Ky6z_D2x_a+ABZ*I_NAI_Iy5z_D2x_a;
  Double I_NAI_I6z_F2xz_a = I_NAI_K7z_D2x_a+ABZ*I_NAI_I6z_D2x_a;
  Double I_NAI_I6x_Fx2y_a = I_NAI_K7x_D2y_a+ABX*I_NAI_I6x_D2y_a;
  Double I_NAI_I5xy_Fx2y_a = I_NAI_K6xy_D2y_a+ABX*I_NAI_I5xy_D2y_a;
  Double I_NAI_I5xz_Fx2y_a = I_NAI_K6xz_D2y_a+ABX*I_NAI_I5xz_D2y_a;
  Double I_NAI_I4x2y_Fx2y_a = I_NAI_K5x2y_D2y_a+ABX*I_NAI_I4x2y_D2y_a;
  Double I_NAI_I4xyz_Fx2y_a = I_NAI_K5xyz_D2y_a+ABX*I_NAI_I4xyz_D2y_a;
  Double I_NAI_I4x2z_Fx2y_a = I_NAI_K5x2z_D2y_a+ABX*I_NAI_I4x2z_D2y_a;
  Double I_NAI_I3x3y_Fx2y_a = I_NAI_K4x3y_D2y_a+ABX*I_NAI_I3x3y_D2y_a;
  Double I_NAI_I3x2yz_Fx2y_a = I_NAI_K4x2yz_D2y_a+ABX*I_NAI_I3x2yz_D2y_a;
  Double I_NAI_I3xy2z_Fx2y_a = I_NAI_K4xy2z_D2y_a+ABX*I_NAI_I3xy2z_D2y_a;
  Double I_NAI_I3x3z_Fx2y_a = I_NAI_K4x3z_D2y_a+ABX*I_NAI_I3x3z_D2y_a;
  Double I_NAI_I2x4y_Fx2y_a = I_NAI_K3x4y_D2y_a+ABX*I_NAI_I2x4y_D2y_a;
  Double I_NAI_I2x3yz_Fx2y_a = I_NAI_K3x3yz_D2y_a+ABX*I_NAI_I2x3yz_D2y_a;
  Double I_NAI_I2x2y2z_Fx2y_a = I_NAI_K3x2y2z_D2y_a+ABX*I_NAI_I2x2y2z_D2y_a;
  Double I_NAI_I2xy3z_Fx2y_a = I_NAI_K3xy3z_D2y_a+ABX*I_NAI_I2xy3z_D2y_a;
  Double I_NAI_I2x4z_Fx2y_a = I_NAI_K3x4z_D2y_a+ABX*I_NAI_I2x4z_D2y_a;
  Double I_NAI_Ix5y_Fx2y_a = I_NAI_K2x5y_D2y_a+ABX*I_NAI_Ix5y_D2y_a;
  Double I_NAI_Ix4yz_Fx2y_a = I_NAI_K2x4yz_D2y_a+ABX*I_NAI_Ix4yz_D2y_a;
  Double I_NAI_Ix3y2z_Fx2y_a = I_NAI_K2x3y2z_D2y_a+ABX*I_NAI_Ix3y2z_D2y_a;
  Double I_NAI_Ix2y3z_Fx2y_a = I_NAI_K2x2y3z_D2y_a+ABX*I_NAI_Ix2y3z_D2y_a;
  Double I_NAI_Ixy4z_Fx2y_a = I_NAI_K2xy4z_D2y_a+ABX*I_NAI_Ixy4z_D2y_a;
  Double I_NAI_Ix5z_Fx2y_a = I_NAI_K2x5z_D2y_a+ABX*I_NAI_Ix5z_D2y_a;
  Double I_NAI_I6y_Fx2y_a = I_NAI_Kx6y_D2y_a+ABX*I_NAI_I6y_D2y_a;
  Double I_NAI_I5yz_Fx2y_a = I_NAI_Kx5yz_D2y_a+ABX*I_NAI_I5yz_D2y_a;
  Double I_NAI_I4y2z_Fx2y_a = I_NAI_Kx4y2z_D2y_a+ABX*I_NAI_I4y2z_D2y_a;
  Double I_NAI_I3y3z_Fx2y_a = I_NAI_Kx3y3z_D2y_a+ABX*I_NAI_I3y3z_D2y_a;
  Double I_NAI_I2y4z_Fx2y_a = I_NAI_Kx2y4z_D2y_a+ABX*I_NAI_I2y4z_D2y_a;
  Double I_NAI_Iy5z_Fx2y_a = I_NAI_Kxy5z_D2y_a+ABX*I_NAI_Iy5z_D2y_a;
  Double I_NAI_I6z_Fx2y_a = I_NAI_Kx6z_D2y_a+ABX*I_NAI_I6z_D2y_a;
  Double I_NAI_I6x_Fxyz_a = I_NAI_K6xz_Dxy_a+ABZ*I_NAI_I6x_Dxy_a;
  Double I_NAI_I5xy_Fxyz_a = I_NAI_K5xyz_Dxy_a+ABZ*I_NAI_I5xy_Dxy_a;
  Double I_NAI_I5xz_Fxyz_a = I_NAI_K5x2z_Dxy_a+ABZ*I_NAI_I5xz_Dxy_a;
  Double I_NAI_I4x2y_Fxyz_a = I_NAI_K4x2yz_Dxy_a+ABZ*I_NAI_I4x2y_Dxy_a;
  Double I_NAI_I4xyz_Fxyz_a = I_NAI_K4xy2z_Dxy_a+ABZ*I_NAI_I4xyz_Dxy_a;
  Double I_NAI_I4x2z_Fxyz_a = I_NAI_K4x3z_Dxy_a+ABZ*I_NAI_I4x2z_Dxy_a;
  Double I_NAI_I3x3y_Fxyz_a = I_NAI_K3x3yz_Dxy_a+ABZ*I_NAI_I3x3y_Dxy_a;
  Double I_NAI_I3x2yz_Fxyz_a = I_NAI_K3x2y2z_Dxy_a+ABZ*I_NAI_I3x2yz_Dxy_a;
  Double I_NAI_I3xy2z_Fxyz_a = I_NAI_K3xy3z_Dxy_a+ABZ*I_NAI_I3xy2z_Dxy_a;
  Double I_NAI_I3x3z_Fxyz_a = I_NAI_K3x4z_Dxy_a+ABZ*I_NAI_I3x3z_Dxy_a;
  Double I_NAI_I2x4y_Fxyz_a = I_NAI_K2x4yz_Dxy_a+ABZ*I_NAI_I2x4y_Dxy_a;
  Double I_NAI_I2x3yz_Fxyz_a = I_NAI_K2x3y2z_Dxy_a+ABZ*I_NAI_I2x3yz_Dxy_a;
  Double I_NAI_I2x2y2z_Fxyz_a = I_NAI_K2x2y3z_Dxy_a+ABZ*I_NAI_I2x2y2z_Dxy_a;
  Double I_NAI_I2xy3z_Fxyz_a = I_NAI_K2xy4z_Dxy_a+ABZ*I_NAI_I2xy3z_Dxy_a;
  Double I_NAI_I2x4z_Fxyz_a = I_NAI_K2x5z_Dxy_a+ABZ*I_NAI_I2x4z_Dxy_a;
  Double I_NAI_Ix5y_Fxyz_a = I_NAI_Kx5yz_Dxy_a+ABZ*I_NAI_Ix5y_Dxy_a;
  Double I_NAI_Ix4yz_Fxyz_a = I_NAI_Kx4y2z_Dxy_a+ABZ*I_NAI_Ix4yz_Dxy_a;
  Double I_NAI_Ix3y2z_Fxyz_a = I_NAI_Kx3y3z_Dxy_a+ABZ*I_NAI_Ix3y2z_Dxy_a;
  Double I_NAI_Ix2y3z_Fxyz_a = I_NAI_Kx2y4z_Dxy_a+ABZ*I_NAI_Ix2y3z_Dxy_a;
  Double I_NAI_Ixy4z_Fxyz_a = I_NAI_Kxy5z_Dxy_a+ABZ*I_NAI_Ixy4z_Dxy_a;
  Double I_NAI_Ix5z_Fxyz_a = I_NAI_Kx6z_Dxy_a+ABZ*I_NAI_Ix5z_Dxy_a;
  Double I_NAI_I6y_Fxyz_a = I_NAI_K6yz_Dxy_a+ABZ*I_NAI_I6y_Dxy_a;
  Double I_NAI_I5yz_Fxyz_a = I_NAI_K5y2z_Dxy_a+ABZ*I_NAI_I5yz_Dxy_a;
  Double I_NAI_I4y2z_Fxyz_a = I_NAI_K4y3z_Dxy_a+ABZ*I_NAI_I4y2z_Dxy_a;
  Double I_NAI_I3y3z_Fxyz_a = I_NAI_K3y4z_Dxy_a+ABZ*I_NAI_I3y3z_Dxy_a;
  Double I_NAI_I2y4z_Fxyz_a = I_NAI_K2y5z_Dxy_a+ABZ*I_NAI_I2y4z_Dxy_a;
  Double I_NAI_Iy5z_Fxyz_a = I_NAI_Ky6z_Dxy_a+ABZ*I_NAI_Iy5z_Dxy_a;
  Double I_NAI_I6z_Fxyz_a = I_NAI_K7z_Dxy_a+ABZ*I_NAI_I6z_Dxy_a;
  Double I_NAI_I6x_Fx2z_a = I_NAI_K7x_D2z_a+ABX*I_NAI_I6x_D2z_a;
  Double I_NAI_I5xy_Fx2z_a = I_NAI_K6xy_D2z_a+ABX*I_NAI_I5xy_D2z_a;
  Double I_NAI_I5xz_Fx2z_a = I_NAI_K6xz_D2z_a+ABX*I_NAI_I5xz_D2z_a;
  Double I_NAI_I4x2y_Fx2z_a = I_NAI_K5x2y_D2z_a+ABX*I_NAI_I4x2y_D2z_a;
  Double I_NAI_I4xyz_Fx2z_a = I_NAI_K5xyz_D2z_a+ABX*I_NAI_I4xyz_D2z_a;
  Double I_NAI_I4x2z_Fx2z_a = I_NAI_K5x2z_D2z_a+ABX*I_NAI_I4x2z_D2z_a;
  Double I_NAI_I3x3y_Fx2z_a = I_NAI_K4x3y_D2z_a+ABX*I_NAI_I3x3y_D2z_a;
  Double I_NAI_I3x2yz_Fx2z_a = I_NAI_K4x2yz_D2z_a+ABX*I_NAI_I3x2yz_D2z_a;
  Double I_NAI_I3xy2z_Fx2z_a = I_NAI_K4xy2z_D2z_a+ABX*I_NAI_I3xy2z_D2z_a;
  Double I_NAI_I3x3z_Fx2z_a = I_NAI_K4x3z_D2z_a+ABX*I_NAI_I3x3z_D2z_a;
  Double I_NAI_I2x4y_Fx2z_a = I_NAI_K3x4y_D2z_a+ABX*I_NAI_I2x4y_D2z_a;
  Double I_NAI_I2x3yz_Fx2z_a = I_NAI_K3x3yz_D2z_a+ABX*I_NAI_I2x3yz_D2z_a;
  Double I_NAI_I2x2y2z_Fx2z_a = I_NAI_K3x2y2z_D2z_a+ABX*I_NAI_I2x2y2z_D2z_a;
  Double I_NAI_I2xy3z_Fx2z_a = I_NAI_K3xy3z_D2z_a+ABX*I_NAI_I2xy3z_D2z_a;
  Double I_NAI_I2x4z_Fx2z_a = I_NAI_K3x4z_D2z_a+ABX*I_NAI_I2x4z_D2z_a;
  Double I_NAI_Ix5y_Fx2z_a = I_NAI_K2x5y_D2z_a+ABX*I_NAI_Ix5y_D2z_a;
  Double I_NAI_Ix4yz_Fx2z_a = I_NAI_K2x4yz_D2z_a+ABX*I_NAI_Ix4yz_D2z_a;
  Double I_NAI_Ix3y2z_Fx2z_a = I_NAI_K2x3y2z_D2z_a+ABX*I_NAI_Ix3y2z_D2z_a;
  Double I_NAI_Ix2y3z_Fx2z_a = I_NAI_K2x2y3z_D2z_a+ABX*I_NAI_Ix2y3z_D2z_a;
  Double I_NAI_Ixy4z_Fx2z_a = I_NAI_K2xy4z_D2z_a+ABX*I_NAI_Ixy4z_D2z_a;
  Double I_NAI_Ix5z_Fx2z_a = I_NAI_K2x5z_D2z_a+ABX*I_NAI_Ix5z_D2z_a;
  Double I_NAI_I6y_Fx2z_a = I_NAI_Kx6y_D2z_a+ABX*I_NAI_I6y_D2z_a;
  Double I_NAI_I5yz_Fx2z_a = I_NAI_Kx5yz_D2z_a+ABX*I_NAI_I5yz_D2z_a;
  Double I_NAI_I4y2z_Fx2z_a = I_NAI_Kx4y2z_D2z_a+ABX*I_NAI_I4y2z_D2z_a;
  Double I_NAI_I3y3z_Fx2z_a = I_NAI_Kx3y3z_D2z_a+ABX*I_NAI_I3y3z_D2z_a;
  Double I_NAI_I2y4z_Fx2z_a = I_NAI_Kx2y4z_D2z_a+ABX*I_NAI_I2y4z_D2z_a;
  Double I_NAI_Iy5z_Fx2z_a = I_NAI_Kxy5z_D2z_a+ABX*I_NAI_Iy5z_D2z_a;
  Double I_NAI_I6z_Fx2z_a = I_NAI_Kx6z_D2z_a+ABX*I_NAI_I6z_D2z_a;
  Double I_NAI_I6x_F3y_a = I_NAI_K6xy_D2y_a+ABY*I_NAI_I6x_D2y_a;
  Double I_NAI_I5xy_F3y_a = I_NAI_K5x2y_D2y_a+ABY*I_NAI_I5xy_D2y_a;
  Double I_NAI_I5xz_F3y_a = I_NAI_K5xyz_D2y_a+ABY*I_NAI_I5xz_D2y_a;
  Double I_NAI_I4x2y_F3y_a = I_NAI_K4x3y_D2y_a+ABY*I_NAI_I4x2y_D2y_a;
  Double I_NAI_I4xyz_F3y_a = I_NAI_K4x2yz_D2y_a+ABY*I_NAI_I4xyz_D2y_a;
  Double I_NAI_I4x2z_F3y_a = I_NAI_K4xy2z_D2y_a+ABY*I_NAI_I4x2z_D2y_a;
  Double I_NAI_I3x3y_F3y_a = I_NAI_K3x4y_D2y_a+ABY*I_NAI_I3x3y_D2y_a;
  Double I_NAI_I3x2yz_F3y_a = I_NAI_K3x3yz_D2y_a+ABY*I_NAI_I3x2yz_D2y_a;
  Double I_NAI_I3xy2z_F3y_a = I_NAI_K3x2y2z_D2y_a+ABY*I_NAI_I3xy2z_D2y_a;
  Double I_NAI_I3x3z_F3y_a = I_NAI_K3xy3z_D2y_a+ABY*I_NAI_I3x3z_D2y_a;
  Double I_NAI_I2x4y_F3y_a = I_NAI_K2x5y_D2y_a+ABY*I_NAI_I2x4y_D2y_a;
  Double I_NAI_I2x3yz_F3y_a = I_NAI_K2x4yz_D2y_a+ABY*I_NAI_I2x3yz_D2y_a;
  Double I_NAI_I2x2y2z_F3y_a = I_NAI_K2x3y2z_D2y_a+ABY*I_NAI_I2x2y2z_D2y_a;
  Double I_NAI_I2xy3z_F3y_a = I_NAI_K2x2y3z_D2y_a+ABY*I_NAI_I2xy3z_D2y_a;
  Double I_NAI_I2x4z_F3y_a = I_NAI_K2xy4z_D2y_a+ABY*I_NAI_I2x4z_D2y_a;
  Double I_NAI_Ix5y_F3y_a = I_NAI_Kx6y_D2y_a+ABY*I_NAI_Ix5y_D2y_a;
  Double I_NAI_Ix4yz_F3y_a = I_NAI_Kx5yz_D2y_a+ABY*I_NAI_Ix4yz_D2y_a;
  Double I_NAI_Ix3y2z_F3y_a = I_NAI_Kx4y2z_D2y_a+ABY*I_NAI_Ix3y2z_D2y_a;
  Double I_NAI_Ix2y3z_F3y_a = I_NAI_Kx3y3z_D2y_a+ABY*I_NAI_Ix2y3z_D2y_a;
  Double I_NAI_Ixy4z_F3y_a = I_NAI_Kx2y4z_D2y_a+ABY*I_NAI_Ixy4z_D2y_a;
  Double I_NAI_Ix5z_F3y_a = I_NAI_Kxy5z_D2y_a+ABY*I_NAI_Ix5z_D2y_a;
  Double I_NAI_I6y_F3y_a = I_NAI_K7y_D2y_a+ABY*I_NAI_I6y_D2y_a;
  Double I_NAI_I5yz_F3y_a = I_NAI_K6yz_D2y_a+ABY*I_NAI_I5yz_D2y_a;
  Double I_NAI_I4y2z_F3y_a = I_NAI_K5y2z_D2y_a+ABY*I_NAI_I4y2z_D2y_a;
  Double I_NAI_I3y3z_F3y_a = I_NAI_K4y3z_D2y_a+ABY*I_NAI_I3y3z_D2y_a;
  Double I_NAI_I2y4z_F3y_a = I_NAI_K3y4z_D2y_a+ABY*I_NAI_I2y4z_D2y_a;
  Double I_NAI_Iy5z_F3y_a = I_NAI_K2y5z_D2y_a+ABY*I_NAI_Iy5z_D2y_a;
  Double I_NAI_I6z_F3y_a = I_NAI_Ky6z_D2y_a+ABY*I_NAI_I6z_D2y_a;
  Double I_NAI_I6x_F2yz_a = I_NAI_K6xz_D2y_a+ABZ*I_NAI_I6x_D2y_a;
  Double I_NAI_I5xy_F2yz_a = I_NAI_K5xyz_D2y_a+ABZ*I_NAI_I5xy_D2y_a;
  Double I_NAI_I5xz_F2yz_a = I_NAI_K5x2z_D2y_a+ABZ*I_NAI_I5xz_D2y_a;
  Double I_NAI_I4x2y_F2yz_a = I_NAI_K4x2yz_D2y_a+ABZ*I_NAI_I4x2y_D2y_a;
  Double I_NAI_I4xyz_F2yz_a = I_NAI_K4xy2z_D2y_a+ABZ*I_NAI_I4xyz_D2y_a;
  Double I_NAI_I4x2z_F2yz_a = I_NAI_K4x3z_D2y_a+ABZ*I_NAI_I4x2z_D2y_a;
  Double I_NAI_I3x3y_F2yz_a = I_NAI_K3x3yz_D2y_a+ABZ*I_NAI_I3x3y_D2y_a;
  Double I_NAI_I3x2yz_F2yz_a = I_NAI_K3x2y2z_D2y_a+ABZ*I_NAI_I3x2yz_D2y_a;
  Double I_NAI_I3xy2z_F2yz_a = I_NAI_K3xy3z_D2y_a+ABZ*I_NAI_I3xy2z_D2y_a;
  Double I_NAI_I3x3z_F2yz_a = I_NAI_K3x4z_D2y_a+ABZ*I_NAI_I3x3z_D2y_a;
  Double I_NAI_I2x4y_F2yz_a = I_NAI_K2x4yz_D2y_a+ABZ*I_NAI_I2x4y_D2y_a;
  Double I_NAI_I2x3yz_F2yz_a = I_NAI_K2x3y2z_D2y_a+ABZ*I_NAI_I2x3yz_D2y_a;
  Double I_NAI_I2x2y2z_F2yz_a = I_NAI_K2x2y3z_D2y_a+ABZ*I_NAI_I2x2y2z_D2y_a;
  Double I_NAI_I2xy3z_F2yz_a = I_NAI_K2xy4z_D2y_a+ABZ*I_NAI_I2xy3z_D2y_a;
  Double I_NAI_I2x4z_F2yz_a = I_NAI_K2x5z_D2y_a+ABZ*I_NAI_I2x4z_D2y_a;
  Double I_NAI_Ix5y_F2yz_a = I_NAI_Kx5yz_D2y_a+ABZ*I_NAI_Ix5y_D2y_a;
  Double I_NAI_Ix4yz_F2yz_a = I_NAI_Kx4y2z_D2y_a+ABZ*I_NAI_Ix4yz_D2y_a;
  Double I_NAI_Ix3y2z_F2yz_a = I_NAI_Kx3y3z_D2y_a+ABZ*I_NAI_Ix3y2z_D2y_a;
  Double I_NAI_Ix2y3z_F2yz_a = I_NAI_Kx2y4z_D2y_a+ABZ*I_NAI_Ix2y3z_D2y_a;
  Double I_NAI_Ixy4z_F2yz_a = I_NAI_Kxy5z_D2y_a+ABZ*I_NAI_Ixy4z_D2y_a;
  Double I_NAI_Ix5z_F2yz_a = I_NAI_Kx6z_D2y_a+ABZ*I_NAI_Ix5z_D2y_a;
  Double I_NAI_I6y_F2yz_a = I_NAI_K6yz_D2y_a+ABZ*I_NAI_I6y_D2y_a;
  Double I_NAI_I5yz_F2yz_a = I_NAI_K5y2z_D2y_a+ABZ*I_NAI_I5yz_D2y_a;
  Double I_NAI_I4y2z_F2yz_a = I_NAI_K4y3z_D2y_a+ABZ*I_NAI_I4y2z_D2y_a;
  Double I_NAI_I3y3z_F2yz_a = I_NAI_K3y4z_D2y_a+ABZ*I_NAI_I3y3z_D2y_a;
  Double I_NAI_I2y4z_F2yz_a = I_NAI_K2y5z_D2y_a+ABZ*I_NAI_I2y4z_D2y_a;
  Double I_NAI_Iy5z_F2yz_a = I_NAI_Ky6z_D2y_a+ABZ*I_NAI_Iy5z_D2y_a;
  Double I_NAI_I6z_F2yz_a = I_NAI_K7z_D2y_a+ABZ*I_NAI_I6z_D2y_a;
  Double I_NAI_I6x_Fy2z_a = I_NAI_K6xy_D2z_a+ABY*I_NAI_I6x_D2z_a;
  Double I_NAI_I5xy_Fy2z_a = I_NAI_K5x2y_D2z_a+ABY*I_NAI_I5xy_D2z_a;
  Double I_NAI_I5xz_Fy2z_a = I_NAI_K5xyz_D2z_a+ABY*I_NAI_I5xz_D2z_a;
  Double I_NAI_I4x2y_Fy2z_a = I_NAI_K4x3y_D2z_a+ABY*I_NAI_I4x2y_D2z_a;
  Double I_NAI_I4xyz_Fy2z_a = I_NAI_K4x2yz_D2z_a+ABY*I_NAI_I4xyz_D2z_a;
  Double I_NAI_I4x2z_Fy2z_a = I_NAI_K4xy2z_D2z_a+ABY*I_NAI_I4x2z_D2z_a;
  Double I_NAI_I3x3y_Fy2z_a = I_NAI_K3x4y_D2z_a+ABY*I_NAI_I3x3y_D2z_a;
  Double I_NAI_I3x2yz_Fy2z_a = I_NAI_K3x3yz_D2z_a+ABY*I_NAI_I3x2yz_D2z_a;
  Double I_NAI_I3xy2z_Fy2z_a = I_NAI_K3x2y2z_D2z_a+ABY*I_NAI_I3xy2z_D2z_a;
  Double I_NAI_I3x3z_Fy2z_a = I_NAI_K3xy3z_D2z_a+ABY*I_NAI_I3x3z_D2z_a;
  Double I_NAI_I2x4y_Fy2z_a = I_NAI_K2x5y_D2z_a+ABY*I_NAI_I2x4y_D2z_a;
  Double I_NAI_I2x3yz_Fy2z_a = I_NAI_K2x4yz_D2z_a+ABY*I_NAI_I2x3yz_D2z_a;
  Double I_NAI_I2x2y2z_Fy2z_a = I_NAI_K2x3y2z_D2z_a+ABY*I_NAI_I2x2y2z_D2z_a;
  Double I_NAI_I2xy3z_Fy2z_a = I_NAI_K2x2y3z_D2z_a+ABY*I_NAI_I2xy3z_D2z_a;
  Double I_NAI_I2x4z_Fy2z_a = I_NAI_K2xy4z_D2z_a+ABY*I_NAI_I2x4z_D2z_a;
  Double I_NAI_Ix5y_Fy2z_a = I_NAI_Kx6y_D2z_a+ABY*I_NAI_Ix5y_D2z_a;
  Double I_NAI_Ix4yz_Fy2z_a = I_NAI_Kx5yz_D2z_a+ABY*I_NAI_Ix4yz_D2z_a;
  Double I_NAI_Ix3y2z_Fy2z_a = I_NAI_Kx4y2z_D2z_a+ABY*I_NAI_Ix3y2z_D2z_a;
  Double I_NAI_Ix2y3z_Fy2z_a = I_NAI_Kx3y3z_D2z_a+ABY*I_NAI_Ix2y3z_D2z_a;
  Double I_NAI_Ixy4z_Fy2z_a = I_NAI_Kx2y4z_D2z_a+ABY*I_NAI_Ixy4z_D2z_a;
  Double I_NAI_Ix5z_Fy2z_a = I_NAI_Kxy5z_D2z_a+ABY*I_NAI_Ix5z_D2z_a;
  Double I_NAI_I6y_Fy2z_a = I_NAI_K7y_D2z_a+ABY*I_NAI_I6y_D2z_a;
  Double I_NAI_I5yz_Fy2z_a = I_NAI_K6yz_D2z_a+ABY*I_NAI_I5yz_D2z_a;
  Double I_NAI_I4y2z_Fy2z_a = I_NAI_K5y2z_D2z_a+ABY*I_NAI_I4y2z_D2z_a;
  Double I_NAI_I3y3z_Fy2z_a = I_NAI_K4y3z_D2z_a+ABY*I_NAI_I3y3z_D2z_a;
  Double I_NAI_I2y4z_Fy2z_a = I_NAI_K3y4z_D2z_a+ABY*I_NAI_I2y4z_D2z_a;
  Double I_NAI_Iy5z_Fy2z_a = I_NAI_K2y5z_D2z_a+ABY*I_NAI_Iy5z_D2z_a;
  Double I_NAI_I6z_Fy2z_a = I_NAI_Ky6z_D2z_a+ABY*I_NAI_I6z_D2z_a;
  Double I_NAI_I6x_F3z_a = I_NAI_K6xz_D2z_a+ABZ*I_NAI_I6x_D2z_a;
  Double I_NAI_I5xy_F3z_a = I_NAI_K5xyz_D2z_a+ABZ*I_NAI_I5xy_D2z_a;
  Double I_NAI_I5xz_F3z_a = I_NAI_K5x2z_D2z_a+ABZ*I_NAI_I5xz_D2z_a;
  Double I_NAI_I4x2y_F3z_a = I_NAI_K4x2yz_D2z_a+ABZ*I_NAI_I4x2y_D2z_a;
  Double I_NAI_I4xyz_F3z_a = I_NAI_K4xy2z_D2z_a+ABZ*I_NAI_I4xyz_D2z_a;
  Double I_NAI_I4x2z_F3z_a = I_NAI_K4x3z_D2z_a+ABZ*I_NAI_I4x2z_D2z_a;
  Double I_NAI_I3x3y_F3z_a = I_NAI_K3x3yz_D2z_a+ABZ*I_NAI_I3x3y_D2z_a;
  Double I_NAI_I3x2yz_F3z_a = I_NAI_K3x2y2z_D2z_a+ABZ*I_NAI_I3x2yz_D2z_a;
  Double I_NAI_I3xy2z_F3z_a = I_NAI_K3xy3z_D2z_a+ABZ*I_NAI_I3xy2z_D2z_a;
  Double I_NAI_I3x3z_F3z_a = I_NAI_K3x4z_D2z_a+ABZ*I_NAI_I3x3z_D2z_a;
  Double I_NAI_I2x4y_F3z_a = I_NAI_K2x4yz_D2z_a+ABZ*I_NAI_I2x4y_D2z_a;
  Double I_NAI_I2x3yz_F3z_a = I_NAI_K2x3y2z_D2z_a+ABZ*I_NAI_I2x3yz_D2z_a;
  Double I_NAI_I2x2y2z_F3z_a = I_NAI_K2x2y3z_D2z_a+ABZ*I_NAI_I2x2y2z_D2z_a;
  Double I_NAI_I2xy3z_F3z_a = I_NAI_K2xy4z_D2z_a+ABZ*I_NAI_I2xy3z_D2z_a;
  Double I_NAI_I2x4z_F3z_a = I_NAI_K2x5z_D2z_a+ABZ*I_NAI_I2x4z_D2z_a;
  Double I_NAI_Ix5y_F3z_a = I_NAI_Kx5yz_D2z_a+ABZ*I_NAI_Ix5y_D2z_a;
  Double I_NAI_Ix4yz_F3z_a = I_NAI_Kx4y2z_D2z_a+ABZ*I_NAI_Ix4yz_D2z_a;
  Double I_NAI_Ix3y2z_F3z_a = I_NAI_Kx3y3z_D2z_a+ABZ*I_NAI_Ix3y2z_D2z_a;
  Double I_NAI_Ix2y3z_F3z_a = I_NAI_Kx2y4z_D2z_a+ABZ*I_NAI_Ix2y3z_D2z_a;
  Double I_NAI_Ixy4z_F3z_a = I_NAI_Kxy5z_D2z_a+ABZ*I_NAI_Ixy4z_D2z_a;
  Double I_NAI_Ix5z_F3z_a = I_NAI_Kx6z_D2z_a+ABZ*I_NAI_Ix5z_D2z_a;
  Double I_NAI_I6y_F3z_a = I_NAI_K6yz_D2z_a+ABZ*I_NAI_I6y_D2z_a;
  Double I_NAI_I5yz_F3z_a = I_NAI_K5y2z_D2z_a+ABZ*I_NAI_I5yz_D2z_a;
  Double I_NAI_I4y2z_F3z_a = I_NAI_K4y3z_D2z_a+ABZ*I_NAI_I4y2z_D2z_a;
  Double I_NAI_I3y3z_F3z_a = I_NAI_K3y4z_D2z_a+ABZ*I_NAI_I3y3z_D2z_a;
  Double I_NAI_I2y4z_F3z_a = I_NAI_K2y5z_D2z_a+ABZ*I_NAI_I2y4z_D2z_a;
  Double I_NAI_Iy5z_F3z_a = I_NAI_Ky6z_D2z_a+ABZ*I_NAI_Iy5z_D2z_a;
  Double I_NAI_I6z_F3z_a = I_NAI_K7z_D2z_a+ABZ*I_NAI_I6z_D2z_a;

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
   * shell quartet name: SQ_NAI_I_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_NAI_I6x_Py_b = I_NAI_K6xy_S_b+ABY*I_NAI_I6x_S_b;
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
  Double I_NAI_I6x_Pz_b = I_NAI_K6xz_S_b+ABZ*I_NAI_I6x_S_b;
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
  Double I_NAI_I6y_Pz_b = I_NAI_K6yz_S_b+ABZ*I_NAI_I6y_S_b;
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
   * shell quartet name: SQ_NAI_K_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
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
  Double I_NAI_K7y_Px_b = I_NAI_Lx7y_S_b+ABX*I_NAI_K7y_S_b;
  Double I_NAI_K6yz_Px_b = I_NAI_Lx6yz_S_b+ABX*I_NAI_K6yz_S_b;
  Double I_NAI_K5y2z_Px_b = I_NAI_Lx5y2z_S_b+ABX*I_NAI_K5y2z_S_b;
  Double I_NAI_K4y3z_Px_b = I_NAI_Lx4y3z_S_b+ABX*I_NAI_K4y3z_S_b;
  Double I_NAI_K3y4z_Px_b = I_NAI_Lx3y4z_S_b+ABX*I_NAI_K3y4z_S_b;
  Double I_NAI_K2y5z_Px_b = I_NAI_Lx2y5z_S_b+ABX*I_NAI_K2y5z_S_b;
  Double I_NAI_Ky6z_Px_b = I_NAI_Lxy6z_S_b+ABX*I_NAI_Ky6z_S_b;
  Double I_NAI_K7z_Px_b = I_NAI_Lx7z_S_b+ABX*I_NAI_K7z_S_b;
  Double I_NAI_K6xy_Py_b = I_NAI_L6x2y_S_b+ABY*I_NAI_K6xy_S_b;
  Double I_NAI_K6xz_Py_b = I_NAI_L6xyz_S_b+ABY*I_NAI_K6xz_S_b;
  Double I_NAI_K5x2y_Py_b = I_NAI_L5x3y_S_b+ABY*I_NAI_K5x2y_S_b;
  Double I_NAI_K5xyz_Py_b = I_NAI_L5x2yz_S_b+ABY*I_NAI_K5xyz_S_b;
  Double I_NAI_K5x2z_Py_b = I_NAI_L5xy2z_S_b+ABY*I_NAI_K5x2z_S_b;
  Double I_NAI_K4x3y_Py_b = I_NAI_L4x4y_S_b+ABY*I_NAI_K4x3y_S_b;
  Double I_NAI_K4x2yz_Py_b = I_NAI_L4x3yz_S_b+ABY*I_NAI_K4x2yz_S_b;
  Double I_NAI_K4xy2z_Py_b = I_NAI_L4x2y2z_S_b+ABY*I_NAI_K4xy2z_S_b;
  Double I_NAI_K4x3z_Py_b = I_NAI_L4xy3z_S_b+ABY*I_NAI_K4x3z_S_b;
  Double I_NAI_K3x4y_Py_b = I_NAI_L3x5y_S_b+ABY*I_NAI_K3x4y_S_b;
  Double I_NAI_K3x3yz_Py_b = I_NAI_L3x4yz_S_b+ABY*I_NAI_K3x3yz_S_b;
  Double I_NAI_K3x2y2z_Py_b = I_NAI_L3x3y2z_S_b+ABY*I_NAI_K3x2y2z_S_b;
  Double I_NAI_K3xy3z_Py_b = I_NAI_L3x2y3z_S_b+ABY*I_NAI_K3xy3z_S_b;
  Double I_NAI_K3x4z_Py_b = I_NAI_L3xy4z_S_b+ABY*I_NAI_K3x4z_S_b;
  Double I_NAI_K2x5y_Py_b = I_NAI_L2x6y_S_b+ABY*I_NAI_K2x5y_S_b;
  Double I_NAI_K2x4yz_Py_b = I_NAI_L2x5yz_S_b+ABY*I_NAI_K2x4yz_S_b;
  Double I_NAI_K2x3y2z_Py_b = I_NAI_L2x4y2z_S_b+ABY*I_NAI_K2x3y2z_S_b;
  Double I_NAI_K2x2y3z_Py_b = I_NAI_L2x3y3z_S_b+ABY*I_NAI_K2x2y3z_S_b;
  Double I_NAI_K2xy4z_Py_b = I_NAI_L2x2y4z_S_b+ABY*I_NAI_K2xy4z_S_b;
  Double I_NAI_K2x5z_Py_b = I_NAI_L2xy5z_S_b+ABY*I_NAI_K2x5z_S_b;
  Double I_NAI_Kx6y_Py_b = I_NAI_Lx7y_S_b+ABY*I_NAI_Kx6y_S_b;
  Double I_NAI_Kx5yz_Py_b = I_NAI_Lx6yz_S_b+ABY*I_NAI_Kx5yz_S_b;
  Double I_NAI_Kx4y2z_Py_b = I_NAI_Lx5y2z_S_b+ABY*I_NAI_Kx4y2z_S_b;
  Double I_NAI_Kx3y3z_Py_b = I_NAI_Lx4y3z_S_b+ABY*I_NAI_Kx3y3z_S_b;
  Double I_NAI_Kx2y4z_Py_b = I_NAI_Lx3y4z_S_b+ABY*I_NAI_Kx2y4z_S_b;
  Double I_NAI_Kxy5z_Py_b = I_NAI_Lx2y5z_S_b+ABY*I_NAI_Kxy5z_S_b;
  Double I_NAI_Kx6z_Py_b = I_NAI_Lxy6z_S_b+ABY*I_NAI_Kx6z_S_b;
  Double I_NAI_K7y_Py_b = I_NAI_L8y_S_b+ABY*I_NAI_K7y_S_b;
  Double I_NAI_K6yz_Py_b = I_NAI_L7yz_S_b+ABY*I_NAI_K6yz_S_b;
  Double I_NAI_K5y2z_Py_b = I_NAI_L6y2z_S_b+ABY*I_NAI_K5y2z_S_b;
  Double I_NAI_K4y3z_Py_b = I_NAI_L5y3z_S_b+ABY*I_NAI_K4y3z_S_b;
  Double I_NAI_K3y4z_Py_b = I_NAI_L4y4z_S_b+ABY*I_NAI_K3y4z_S_b;
  Double I_NAI_K2y5z_Py_b = I_NAI_L3y5z_S_b+ABY*I_NAI_K2y5z_S_b;
  Double I_NAI_Ky6z_Py_b = I_NAI_L2y6z_S_b+ABY*I_NAI_Ky6z_S_b;
  Double I_NAI_K7z_Py_b = I_NAI_Ly7z_S_b+ABY*I_NAI_K7z_S_b;
  Double I_NAI_K6xy_Pz_b = I_NAI_L6xyz_S_b+ABZ*I_NAI_K6xy_S_b;
  Double I_NAI_K6xz_Pz_b = I_NAI_L6x2z_S_b+ABZ*I_NAI_K6xz_S_b;
  Double I_NAI_K5x2y_Pz_b = I_NAI_L5x2yz_S_b+ABZ*I_NAI_K5x2y_S_b;
  Double I_NAI_K5xyz_Pz_b = I_NAI_L5xy2z_S_b+ABZ*I_NAI_K5xyz_S_b;
  Double I_NAI_K5x2z_Pz_b = I_NAI_L5x3z_S_b+ABZ*I_NAI_K5x2z_S_b;
  Double I_NAI_K4x3y_Pz_b = I_NAI_L4x3yz_S_b+ABZ*I_NAI_K4x3y_S_b;
  Double I_NAI_K4x2yz_Pz_b = I_NAI_L4x2y2z_S_b+ABZ*I_NAI_K4x2yz_S_b;
  Double I_NAI_K4xy2z_Pz_b = I_NAI_L4xy3z_S_b+ABZ*I_NAI_K4xy2z_S_b;
  Double I_NAI_K4x3z_Pz_b = I_NAI_L4x4z_S_b+ABZ*I_NAI_K4x3z_S_b;
  Double I_NAI_K3x4y_Pz_b = I_NAI_L3x4yz_S_b+ABZ*I_NAI_K3x4y_S_b;
  Double I_NAI_K3x3yz_Pz_b = I_NAI_L3x3y2z_S_b+ABZ*I_NAI_K3x3yz_S_b;
  Double I_NAI_K3x2y2z_Pz_b = I_NAI_L3x2y3z_S_b+ABZ*I_NAI_K3x2y2z_S_b;
  Double I_NAI_K3xy3z_Pz_b = I_NAI_L3xy4z_S_b+ABZ*I_NAI_K3xy3z_S_b;
  Double I_NAI_K3x4z_Pz_b = I_NAI_L3x5z_S_b+ABZ*I_NAI_K3x4z_S_b;
  Double I_NAI_K2x5y_Pz_b = I_NAI_L2x5yz_S_b+ABZ*I_NAI_K2x5y_S_b;
  Double I_NAI_K2x4yz_Pz_b = I_NAI_L2x4y2z_S_b+ABZ*I_NAI_K2x4yz_S_b;
  Double I_NAI_K2x3y2z_Pz_b = I_NAI_L2x3y3z_S_b+ABZ*I_NAI_K2x3y2z_S_b;
  Double I_NAI_K2x2y3z_Pz_b = I_NAI_L2x2y4z_S_b+ABZ*I_NAI_K2x2y3z_S_b;
  Double I_NAI_K2xy4z_Pz_b = I_NAI_L2xy5z_S_b+ABZ*I_NAI_K2xy4z_S_b;
  Double I_NAI_K2x5z_Pz_b = I_NAI_L2x6z_S_b+ABZ*I_NAI_K2x5z_S_b;
  Double I_NAI_Kx6y_Pz_b = I_NAI_Lx6yz_S_b+ABZ*I_NAI_Kx6y_S_b;
  Double I_NAI_Kx5yz_Pz_b = I_NAI_Lx5y2z_S_b+ABZ*I_NAI_Kx5yz_S_b;
  Double I_NAI_Kx4y2z_Pz_b = I_NAI_Lx4y3z_S_b+ABZ*I_NAI_Kx4y2z_S_b;
  Double I_NAI_Kx3y3z_Pz_b = I_NAI_Lx3y4z_S_b+ABZ*I_NAI_Kx3y3z_S_b;
  Double I_NAI_Kx2y4z_Pz_b = I_NAI_Lx2y5z_S_b+ABZ*I_NAI_Kx2y4z_S_b;
  Double I_NAI_Kxy5z_Pz_b = I_NAI_Lxy6z_S_b+ABZ*I_NAI_Kxy5z_S_b;
  Double I_NAI_Kx6z_Pz_b = I_NAI_Lx7z_S_b+ABZ*I_NAI_Kx6z_S_b;
  Double I_NAI_K6yz_Pz_b = I_NAI_L6y2z_S_b+ABZ*I_NAI_K6yz_S_b;
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
   * totally 84 integrals are omitted 
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
  Double I_NAI_I6x_D2y_b = I_NAI_K6xy_Py_b+ABY*I_NAI_I6x_Py_b;
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
  Double I_NAI_I6x_D2z_b = I_NAI_K6xz_Pz_b+ABZ*I_NAI_I6x_Pz_b;
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
  Double I_NAI_I6y_D2z_b = I_NAI_K6yz_Pz_b+ABZ*I_NAI_I6y_Pz_b;
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
   * totally 42 integrals are omitted 
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
  Double I_NAI_H5x_F2xy_b = I_NAI_I5xy_D2x_b+ABY*I_NAI_H5x_D2x_b;
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
  Double I_NAI_H5x_F2xz_b = I_NAI_I5xz_D2x_b+ABZ*I_NAI_H5x_D2x_b;
  Double I_NAI_H4xy_F2xz_b = I_NAI_I4xyz_D2x_b+ABZ*I_NAI_H4xy_D2x_b;
  Double I_NAI_H4xz_F2xz_b = I_NAI_I4x2z_D2x_b+ABZ*I_NAI_H4xz_D2x_b;
  Double I_NAI_H3x2y_F2xz_b = I_NAI_I3x2yz_D2x_b+ABZ*I_NAI_H3x2y_D2x_b;
  Double I_NAI_H3xyz_F2xz_b = I_NAI_I3xy2z_D2x_b+ABZ*I_NAI_H3xyz_D2x_b;
  Double I_NAI_H3x2z_F2xz_b = I_NAI_I3x3z_D2x_b+ABZ*I_NAI_H3x2z_D2x_b;
  Double I_NAI_H2x3y_F2xz_b = I_NAI_I2x3yz_D2x_b+ABZ*I_NAI_H2x3y_D2x_b;
  Double I_NAI_H2x2yz_F2xz_b = I_NAI_I2x2y2z_D2x_b+ABZ*I_NAI_H2x2yz_D2x_b;
  Double I_NAI_H2xy2z_F2xz_b = I_NAI_I2xy3z_D2x_b+ABZ*I_NAI_H2xy2z_D2x_b;
  Double I_NAI_H2x3z_F2xz_b = I_NAI_I2x4z_D2x_b+ABZ*I_NAI_H2x3z_D2x_b;
  Double I_NAI_Hx4y_F2xz_b = I_NAI_Ix4yz_D2x_b+ABZ*I_NAI_Hx4y_D2x_b;
  Double I_NAI_Hx3yz_F2xz_b = I_NAI_Ix3y2z_D2x_b+ABZ*I_NAI_Hx3yz_D2x_b;
  Double I_NAI_Hx2y2z_F2xz_b = I_NAI_Ix2y3z_D2x_b+ABZ*I_NAI_Hx2y2z_D2x_b;
  Double I_NAI_Hxy3z_F2xz_b = I_NAI_Ixy4z_D2x_b+ABZ*I_NAI_Hxy3z_D2x_b;
  Double I_NAI_Hx4z_F2xz_b = I_NAI_Ix5z_D2x_b+ABZ*I_NAI_Hx4z_D2x_b;
  Double I_NAI_H5y_F2xz_b = I_NAI_I5yz_D2x_b+ABZ*I_NAI_H5y_D2x_b;
  Double I_NAI_H4yz_F2xz_b = I_NAI_I4y2z_D2x_b+ABZ*I_NAI_H4yz_D2x_b;
  Double I_NAI_H3y2z_F2xz_b = I_NAI_I3y3z_D2x_b+ABZ*I_NAI_H3y2z_D2x_b;
  Double I_NAI_H2y3z_F2xz_b = I_NAI_I2y4z_D2x_b+ABZ*I_NAI_H2y3z_D2x_b;
  Double I_NAI_Hy4z_F2xz_b = I_NAI_Iy5z_D2x_b+ABZ*I_NAI_Hy4z_D2x_b;
  Double I_NAI_H5z_F2xz_b = I_NAI_I6z_D2x_b+ABZ*I_NAI_H5z_D2x_b;
  Double I_NAI_H5x_Fx2y_b = I_NAI_I6x_D2y_b+ABX*I_NAI_H5x_D2y_b;
  Double I_NAI_H4xy_Fx2y_b = I_NAI_I5xy_D2y_b+ABX*I_NAI_H4xy_D2y_b;
  Double I_NAI_H4xz_Fx2y_b = I_NAI_I5xz_D2y_b+ABX*I_NAI_H4xz_D2y_b;
  Double I_NAI_H3x2y_Fx2y_b = I_NAI_I4x2y_D2y_b+ABX*I_NAI_H3x2y_D2y_b;
  Double I_NAI_H3xyz_Fx2y_b = I_NAI_I4xyz_D2y_b+ABX*I_NAI_H3xyz_D2y_b;
  Double I_NAI_H3x2z_Fx2y_b = I_NAI_I4x2z_D2y_b+ABX*I_NAI_H3x2z_D2y_b;
  Double I_NAI_H2x3y_Fx2y_b = I_NAI_I3x3y_D2y_b+ABX*I_NAI_H2x3y_D2y_b;
  Double I_NAI_H2x2yz_Fx2y_b = I_NAI_I3x2yz_D2y_b+ABX*I_NAI_H2x2yz_D2y_b;
  Double I_NAI_H2xy2z_Fx2y_b = I_NAI_I3xy2z_D2y_b+ABX*I_NAI_H2xy2z_D2y_b;
  Double I_NAI_H2x3z_Fx2y_b = I_NAI_I3x3z_D2y_b+ABX*I_NAI_H2x3z_D2y_b;
  Double I_NAI_Hx4y_Fx2y_b = I_NAI_I2x4y_D2y_b+ABX*I_NAI_Hx4y_D2y_b;
  Double I_NAI_Hx3yz_Fx2y_b = I_NAI_I2x3yz_D2y_b+ABX*I_NAI_Hx3yz_D2y_b;
  Double I_NAI_Hx2y2z_Fx2y_b = I_NAI_I2x2y2z_D2y_b+ABX*I_NAI_Hx2y2z_D2y_b;
  Double I_NAI_Hxy3z_Fx2y_b = I_NAI_I2xy3z_D2y_b+ABX*I_NAI_Hxy3z_D2y_b;
  Double I_NAI_Hx4z_Fx2y_b = I_NAI_I2x4z_D2y_b+ABX*I_NAI_Hx4z_D2y_b;
  Double I_NAI_H5y_Fx2y_b = I_NAI_Ix5y_D2y_b+ABX*I_NAI_H5y_D2y_b;
  Double I_NAI_H4yz_Fx2y_b = I_NAI_Ix4yz_D2y_b+ABX*I_NAI_H4yz_D2y_b;
  Double I_NAI_H3y2z_Fx2y_b = I_NAI_Ix3y2z_D2y_b+ABX*I_NAI_H3y2z_D2y_b;
  Double I_NAI_H2y3z_Fx2y_b = I_NAI_Ix2y3z_D2y_b+ABX*I_NAI_H2y3z_D2y_b;
  Double I_NAI_Hy4z_Fx2y_b = I_NAI_Ixy4z_D2y_b+ABX*I_NAI_Hy4z_D2y_b;
  Double I_NAI_H5z_Fx2y_b = I_NAI_Ix5z_D2y_b+ABX*I_NAI_H5z_D2y_b;
  Double I_NAI_H5x_Fx2z_b = I_NAI_I6x_D2z_b+ABX*I_NAI_H5x_D2z_b;
  Double I_NAI_H4xy_Fx2z_b = I_NAI_I5xy_D2z_b+ABX*I_NAI_H4xy_D2z_b;
  Double I_NAI_H4xz_Fx2z_b = I_NAI_I5xz_D2z_b+ABX*I_NAI_H4xz_D2z_b;
  Double I_NAI_H3x2y_Fx2z_b = I_NAI_I4x2y_D2z_b+ABX*I_NAI_H3x2y_D2z_b;
  Double I_NAI_H3xyz_Fx2z_b = I_NAI_I4xyz_D2z_b+ABX*I_NAI_H3xyz_D2z_b;
  Double I_NAI_H3x2z_Fx2z_b = I_NAI_I4x2z_D2z_b+ABX*I_NAI_H3x2z_D2z_b;
  Double I_NAI_H2x3y_Fx2z_b = I_NAI_I3x3y_D2z_b+ABX*I_NAI_H2x3y_D2z_b;
  Double I_NAI_H2x2yz_Fx2z_b = I_NAI_I3x2yz_D2z_b+ABX*I_NAI_H2x2yz_D2z_b;
  Double I_NAI_H2xy2z_Fx2z_b = I_NAI_I3xy2z_D2z_b+ABX*I_NAI_H2xy2z_D2z_b;
  Double I_NAI_H2x3z_Fx2z_b = I_NAI_I3x3z_D2z_b+ABX*I_NAI_H2x3z_D2z_b;
  Double I_NAI_Hx4y_Fx2z_b = I_NAI_I2x4y_D2z_b+ABX*I_NAI_Hx4y_D2z_b;
  Double I_NAI_Hx3yz_Fx2z_b = I_NAI_I2x3yz_D2z_b+ABX*I_NAI_Hx3yz_D2z_b;
  Double I_NAI_Hx2y2z_Fx2z_b = I_NAI_I2x2y2z_D2z_b+ABX*I_NAI_Hx2y2z_D2z_b;
  Double I_NAI_Hxy3z_Fx2z_b = I_NAI_I2xy3z_D2z_b+ABX*I_NAI_Hxy3z_D2z_b;
  Double I_NAI_Hx4z_Fx2z_b = I_NAI_I2x4z_D2z_b+ABX*I_NAI_Hx4z_D2z_b;
  Double I_NAI_H5y_Fx2z_b = I_NAI_Ix5y_D2z_b+ABX*I_NAI_H5y_D2z_b;
  Double I_NAI_H4yz_Fx2z_b = I_NAI_Ix4yz_D2z_b+ABX*I_NAI_H4yz_D2z_b;
  Double I_NAI_H3y2z_Fx2z_b = I_NAI_Ix3y2z_D2z_b+ABX*I_NAI_H3y2z_D2z_b;
  Double I_NAI_H2y3z_Fx2z_b = I_NAI_Ix2y3z_D2z_b+ABX*I_NAI_H2y3z_D2z_b;
  Double I_NAI_Hy4z_Fx2z_b = I_NAI_Ixy4z_D2z_b+ABX*I_NAI_Hy4z_D2z_b;
  Double I_NAI_H5z_Fx2z_b = I_NAI_Ix5z_D2z_b+ABX*I_NAI_H5z_D2z_b;
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
  Double I_NAI_H5x_F2yz_b = I_NAI_I5xz_D2y_b+ABZ*I_NAI_H5x_D2y_b;
  Double I_NAI_H4xy_F2yz_b = I_NAI_I4xyz_D2y_b+ABZ*I_NAI_H4xy_D2y_b;
  Double I_NAI_H4xz_F2yz_b = I_NAI_I4x2z_D2y_b+ABZ*I_NAI_H4xz_D2y_b;
  Double I_NAI_H3x2y_F2yz_b = I_NAI_I3x2yz_D2y_b+ABZ*I_NAI_H3x2y_D2y_b;
  Double I_NAI_H3xyz_F2yz_b = I_NAI_I3xy2z_D2y_b+ABZ*I_NAI_H3xyz_D2y_b;
  Double I_NAI_H3x2z_F2yz_b = I_NAI_I3x3z_D2y_b+ABZ*I_NAI_H3x2z_D2y_b;
  Double I_NAI_H2x3y_F2yz_b = I_NAI_I2x3yz_D2y_b+ABZ*I_NAI_H2x3y_D2y_b;
  Double I_NAI_H2x2yz_F2yz_b = I_NAI_I2x2y2z_D2y_b+ABZ*I_NAI_H2x2yz_D2y_b;
  Double I_NAI_H2xy2z_F2yz_b = I_NAI_I2xy3z_D2y_b+ABZ*I_NAI_H2xy2z_D2y_b;
  Double I_NAI_H2x3z_F2yz_b = I_NAI_I2x4z_D2y_b+ABZ*I_NAI_H2x3z_D2y_b;
  Double I_NAI_Hx4y_F2yz_b = I_NAI_Ix4yz_D2y_b+ABZ*I_NAI_Hx4y_D2y_b;
  Double I_NAI_Hx3yz_F2yz_b = I_NAI_Ix3y2z_D2y_b+ABZ*I_NAI_Hx3yz_D2y_b;
  Double I_NAI_Hx2y2z_F2yz_b = I_NAI_Ix2y3z_D2y_b+ABZ*I_NAI_Hx2y2z_D2y_b;
  Double I_NAI_Hxy3z_F2yz_b = I_NAI_Ixy4z_D2y_b+ABZ*I_NAI_Hxy3z_D2y_b;
  Double I_NAI_Hx4z_F2yz_b = I_NAI_Ix5z_D2y_b+ABZ*I_NAI_Hx4z_D2y_b;
  Double I_NAI_H5y_F2yz_b = I_NAI_I5yz_D2y_b+ABZ*I_NAI_H5y_D2y_b;
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
   * shell quartet name: SQ_NAI_L_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_S_b
   * RHS shell quartet name: SQ_NAI_L_S_b
   ************************************************************/
  Double I_NAI_L8x_Px_b = I_NAI_M9x_S_b+ABX*I_NAI_L8x_S_b;
  Double I_NAI_L7xy_Px_b = I_NAI_M8xy_S_b+ABX*I_NAI_L7xy_S_b;
  Double I_NAI_L7xz_Px_b = I_NAI_M8xz_S_b+ABX*I_NAI_L7xz_S_b;
  Double I_NAI_L6x2y_Px_b = I_NAI_M7x2y_S_b+ABX*I_NAI_L6x2y_S_b;
  Double I_NAI_L6xyz_Px_b = I_NAI_M7xyz_S_b+ABX*I_NAI_L6xyz_S_b;
  Double I_NAI_L6x2z_Px_b = I_NAI_M7x2z_S_b+ABX*I_NAI_L6x2z_S_b;
  Double I_NAI_L5x3y_Px_b = I_NAI_M6x3y_S_b+ABX*I_NAI_L5x3y_S_b;
  Double I_NAI_L5x2yz_Px_b = I_NAI_M6x2yz_S_b+ABX*I_NAI_L5x2yz_S_b;
  Double I_NAI_L5xy2z_Px_b = I_NAI_M6xy2z_S_b+ABX*I_NAI_L5xy2z_S_b;
  Double I_NAI_L5x3z_Px_b = I_NAI_M6x3z_S_b+ABX*I_NAI_L5x3z_S_b;
  Double I_NAI_L4x4y_Px_b = I_NAI_M5x4y_S_b+ABX*I_NAI_L4x4y_S_b;
  Double I_NAI_L4x3yz_Px_b = I_NAI_M5x3yz_S_b+ABX*I_NAI_L4x3yz_S_b;
  Double I_NAI_L4x2y2z_Px_b = I_NAI_M5x2y2z_S_b+ABX*I_NAI_L4x2y2z_S_b;
  Double I_NAI_L4xy3z_Px_b = I_NAI_M5xy3z_S_b+ABX*I_NAI_L4xy3z_S_b;
  Double I_NAI_L4x4z_Px_b = I_NAI_M5x4z_S_b+ABX*I_NAI_L4x4z_S_b;
  Double I_NAI_L3x5y_Px_b = I_NAI_M4x5y_S_b+ABX*I_NAI_L3x5y_S_b;
  Double I_NAI_L3x4yz_Px_b = I_NAI_M4x4yz_S_b+ABX*I_NAI_L3x4yz_S_b;
  Double I_NAI_L3x3y2z_Px_b = I_NAI_M4x3y2z_S_b+ABX*I_NAI_L3x3y2z_S_b;
  Double I_NAI_L3x2y3z_Px_b = I_NAI_M4x2y3z_S_b+ABX*I_NAI_L3x2y3z_S_b;
  Double I_NAI_L3xy4z_Px_b = I_NAI_M4xy4z_S_b+ABX*I_NAI_L3xy4z_S_b;
  Double I_NAI_L3x5z_Px_b = I_NAI_M4x5z_S_b+ABX*I_NAI_L3x5z_S_b;
  Double I_NAI_L2x6y_Px_b = I_NAI_M3x6y_S_b+ABX*I_NAI_L2x6y_S_b;
  Double I_NAI_L2x5yz_Px_b = I_NAI_M3x5yz_S_b+ABX*I_NAI_L2x5yz_S_b;
  Double I_NAI_L2x4y2z_Px_b = I_NAI_M3x4y2z_S_b+ABX*I_NAI_L2x4y2z_S_b;
  Double I_NAI_L2x3y3z_Px_b = I_NAI_M3x3y3z_S_b+ABX*I_NAI_L2x3y3z_S_b;
  Double I_NAI_L2x2y4z_Px_b = I_NAI_M3x2y4z_S_b+ABX*I_NAI_L2x2y4z_S_b;
  Double I_NAI_L2xy5z_Px_b = I_NAI_M3xy5z_S_b+ABX*I_NAI_L2xy5z_S_b;
  Double I_NAI_L2x6z_Px_b = I_NAI_M3x6z_S_b+ABX*I_NAI_L2x6z_S_b;
  Double I_NAI_Lx7y_Px_b = I_NAI_M2x7y_S_b+ABX*I_NAI_Lx7y_S_b;
  Double I_NAI_Lx6yz_Px_b = I_NAI_M2x6yz_S_b+ABX*I_NAI_Lx6yz_S_b;
  Double I_NAI_Lx5y2z_Px_b = I_NAI_M2x5y2z_S_b+ABX*I_NAI_Lx5y2z_S_b;
  Double I_NAI_Lx4y3z_Px_b = I_NAI_M2x4y3z_S_b+ABX*I_NAI_Lx4y3z_S_b;
  Double I_NAI_Lx3y4z_Px_b = I_NAI_M2x3y4z_S_b+ABX*I_NAI_Lx3y4z_S_b;
  Double I_NAI_Lx2y5z_Px_b = I_NAI_M2x2y5z_S_b+ABX*I_NAI_Lx2y5z_S_b;
  Double I_NAI_Lxy6z_Px_b = I_NAI_M2xy6z_S_b+ABX*I_NAI_Lxy6z_S_b;
  Double I_NAI_Lx7z_Px_b = I_NAI_M2x7z_S_b+ABX*I_NAI_Lx7z_S_b;
  Double I_NAI_L6x2y_Py_b = I_NAI_M6x3y_S_b+ABY*I_NAI_L6x2y_S_b;
  Double I_NAI_L6xyz_Py_b = I_NAI_M6x2yz_S_b+ABY*I_NAI_L6xyz_S_b;
  Double I_NAI_L5x3y_Py_b = I_NAI_M5x4y_S_b+ABY*I_NAI_L5x3y_S_b;
  Double I_NAI_L5x2yz_Py_b = I_NAI_M5x3yz_S_b+ABY*I_NAI_L5x2yz_S_b;
  Double I_NAI_L5xy2z_Py_b = I_NAI_M5x2y2z_S_b+ABY*I_NAI_L5xy2z_S_b;
  Double I_NAI_L4x4y_Py_b = I_NAI_M4x5y_S_b+ABY*I_NAI_L4x4y_S_b;
  Double I_NAI_L4x3yz_Py_b = I_NAI_M4x4yz_S_b+ABY*I_NAI_L4x3yz_S_b;
  Double I_NAI_L4x2y2z_Py_b = I_NAI_M4x3y2z_S_b+ABY*I_NAI_L4x2y2z_S_b;
  Double I_NAI_L4xy3z_Py_b = I_NAI_M4x2y3z_S_b+ABY*I_NAI_L4xy3z_S_b;
  Double I_NAI_L3x5y_Py_b = I_NAI_M3x6y_S_b+ABY*I_NAI_L3x5y_S_b;
  Double I_NAI_L3x4yz_Py_b = I_NAI_M3x5yz_S_b+ABY*I_NAI_L3x4yz_S_b;
  Double I_NAI_L3x3y2z_Py_b = I_NAI_M3x4y2z_S_b+ABY*I_NAI_L3x3y2z_S_b;
  Double I_NAI_L3x2y3z_Py_b = I_NAI_M3x3y3z_S_b+ABY*I_NAI_L3x2y3z_S_b;
  Double I_NAI_L3xy4z_Py_b = I_NAI_M3x2y4z_S_b+ABY*I_NAI_L3xy4z_S_b;
  Double I_NAI_L2x6y_Py_b = I_NAI_M2x7y_S_b+ABY*I_NAI_L2x6y_S_b;
  Double I_NAI_L2x5yz_Py_b = I_NAI_M2x6yz_S_b+ABY*I_NAI_L2x5yz_S_b;
  Double I_NAI_L2x4y2z_Py_b = I_NAI_M2x5y2z_S_b+ABY*I_NAI_L2x4y2z_S_b;
  Double I_NAI_L2x3y3z_Py_b = I_NAI_M2x4y3z_S_b+ABY*I_NAI_L2x3y3z_S_b;
  Double I_NAI_L2x2y4z_Py_b = I_NAI_M2x3y4z_S_b+ABY*I_NAI_L2x2y4z_S_b;
  Double I_NAI_L2xy5z_Py_b = I_NAI_M2x2y5z_S_b+ABY*I_NAI_L2xy5z_S_b;
  Double I_NAI_Lx7y_Py_b = I_NAI_Mx8y_S_b+ABY*I_NAI_Lx7y_S_b;
  Double I_NAI_Lx6yz_Py_b = I_NAI_Mx7yz_S_b+ABY*I_NAI_Lx6yz_S_b;
  Double I_NAI_Lx5y2z_Py_b = I_NAI_Mx6y2z_S_b+ABY*I_NAI_Lx5y2z_S_b;
  Double I_NAI_Lx4y3z_Py_b = I_NAI_Mx5y3z_S_b+ABY*I_NAI_Lx4y3z_S_b;
  Double I_NAI_Lx3y4z_Py_b = I_NAI_Mx4y4z_S_b+ABY*I_NAI_Lx3y4z_S_b;
  Double I_NAI_Lx2y5z_Py_b = I_NAI_Mx3y5z_S_b+ABY*I_NAI_Lx2y5z_S_b;
  Double I_NAI_Lxy6z_Py_b = I_NAI_Mx2y6z_S_b+ABY*I_NAI_Lxy6z_S_b;
  Double I_NAI_L8y_Py_b = I_NAI_M9y_S_b+ABY*I_NAI_L8y_S_b;
  Double I_NAI_L7yz_Py_b = I_NAI_M8yz_S_b+ABY*I_NAI_L7yz_S_b;
  Double I_NAI_L6y2z_Py_b = I_NAI_M7y2z_S_b+ABY*I_NAI_L6y2z_S_b;
  Double I_NAI_L5y3z_Py_b = I_NAI_M6y3z_S_b+ABY*I_NAI_L5y3z_S_b;
  Double I_NAI_L4y4z_Py_b = I_NAI_M5y4z_S_b+ABY*I_NAI_L4y4z_S_b;
  Double I_NAI_L3y5z_Py_b = I_NAI_M4y5z_S_b+ABY*I_NAI_L3y5z_S_b;
  Double I_NAI_L2y6z_Py_b = I_NAI_M3y6z_S_b+ABY*I_NAI_L2y6z_S_b;
  Double I_NAI_Ly7z_Py_b = I_NAI_M2y7z_S_b+ABY*I_NAI_Ly7z_S_b;
  Double I_NAI_L6xyz_Pz_b = I_NAI_M6xy2z_S_b+ABZ*I_NAI_L6xyz_S_b;
  Double I_NAI_L6x2z_Pz_b = I_NAI_M6x3z_S_b+ABZ*I_NAI_L6x2z_S_b;
  Double I_NAI_L5x2yz_Pz_b = I_NAI_M5x2y2z_S_b+ABZ*I_NAI_L5x2yz_S_b;
  Double I_NAI_L5xy2z_Pz_b = I_NAI_M5xy3z_S_b+ABZ*I_NAI_L5xy2z_S_b;
  Double I_NAI_L5x3z_Pz_b = I_NAI_M5x4z_S_b+ABZ*I_NAI_L5x3z_S_b;
  Double I_NAI_L4x3yz_Pz_b = I_NAI_M4x3y2z_S_b+ABZ*I_NAI_L4x3yz_S_b;
  Double I_NAI_L4x2y2z_Pz_b = I_NAI_M4x2y3z_S_b+ABZ*I_NAI_L4x2y2z_S_b;
  Double I_NAI_L4xy3z_Pz_b = I_NAI_M4xy4z_S_b+ABZ*I_NAI_L4xy3z_S_b;
  Double I_NAI_L4x4z_Pz_b = I_NAI_M4x5z_S_b+ABZ*I_NAI_L4x4z_S_b;
  Double I_NAI_L3x4yz_Pz_b = I_NAI_M3x4y2z_S_b+ABZ*I_NAI_L3x4yz_S_b;
  Double I_NAI_L3x3y2z_Pz_b = I_NAI_M3x3y3z_S_b+ABZ*I_NAI_L3x3y2z_S_b;
  Double I_NAI_L3x2y3z_Pz_b = I_NAI_M3x2y4z_S_b+ABZ*I_NAI_L3x2y3z_S_b;
  Double I_NAI_L3xy4z_Pz_b = I_NAI_M3xy5z_S_b+ABZ*I_NAI_L3xy4z_S_b;
  Double I_NAI_L3x5z_Pz_b = I_NAI_M3x6z_S_b+ABZ*I_NAI_L3x5z_S_b;
  Double I_NAI_L2x5yz_Pz_b = I_NAI_M2x5y2z_S_b+ABZ*I_NAI_L2x5yz_S_b;
  Double I_NAI_L2x4y2z_Pz_b = I_NAI_M2x4y3z_S_b+ABZ*I_NAI_L2x4y2z_S_b;
  Double I_NAI_L2x3y3z_Pz_b = I_NAI_M2x3y4z_S_b+ABZ*I_NAI_L2x3y3z_S_b;
  Double I_NAI_L2x2y4z_Pz_b = I_NAI_M2x2y5z_S_b+ABZ*I_NAI_L2x2y4z_S_b;
  Double I_NAI_L2xy5z_Pz_b = I_NAI_M2xy6z_S_b+ABZ*I_NAI_L2xy5z_S_b;
  Double I_NAI_L2x6z_Pz_b = I_NAI_M2x7z_S_b+ABZ*I_NAI_L2x6z_S_b;
  Double I_NAI_Lx6yz_Pz_b = I_NAI_Mx6y2z_S_b+ABZ*I_NAI_Lx6yz_S_b;
  Double I_NAI_Lx5y2z_Pz_b = I_NAI_Mx5y3z_S_b+ABZ*I_NAI_Lx5y2z_S_b;
  Double I_NAI_Lx4y3z_Pz_b = I_NAI_Mx4y4z_S_b+ABZ*I_NAI_Lx4y3z_S_b;
  Double I_NAI_Lx3y4z_Pz_b = I_NAI_Mx3y5z_S_b+ABZ*I_NAI_Lx3y4z_S_b;
  Double I_NAI_Lx2y5z_Pz_b = I_NAI_Mx2y6z_S_b+ABZ*I_NAI_Lx2y5z_S_b;
  Double I_NAI_Lxy6z_Pz_b = I_NAI_Mxy7z_S_b+ABZ*I_NAI_Lxy6z_S_b;
  Double I_NAI_Lx7z_Pz_b = I_NAI_Mx8z_S_b+ABZ*I_NAI_Lx7z_S_b;
  Double I_NAI_L6y2z_Pz_b = I_NAI_M6y3z_S_b+ABZ*I_NAI_L6y2z_S_b;
  Double I_NAI_L5y3z_Pz_b = I_NAI_M5y4z_S_b+ABZ*I_NAI_L5y3z_S_b;
  Double I_NAI_L4y4z_Pz_b = I_NAI_M4y5z_S_b+ABZ*I_NAI_L4y4z_S_b;
  Double I_NAI_L3y5z_Pz_b = I_NAI_M3y6z_S_b+ABZ*I_NAI_L3y5z_S_b;
  Double I_NAI_L2y6z_Pz_b = I_NAI_M2y7z_S_b+ABZ*I_NAI_L2y6z_S_b;
  Double I_NAI_Ly7z_Pz_b = I_NAI_My8z_S_b+ABZ*I_NAI_Ly7z_S_b;
  Double I_NAI_L8z_Pz_b = I_NAI_M9z_S_b+ABZ*I_NAI_L8z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_K_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 111 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_P_b
   * RHS shell quartet name: SQ_NAI_K_P_b
   ************************************************************/
  Double I_NAI_K7x_D2x_b = I_NAI_L8x_Px_b+ABX*I_NAI_K7x_Px_b;
  Double I_NAI_K6xy_D2x_b = I_NAI_L7xy_Px_b+ABX*I_NAI_K6xy_Px_b;
  Double I_NAI_K6xz_D2x_b = I_NAI_L7xz_Px_b+ABX*I_NAI_K6xz_Px_b;
  Double I_NAI_K5x2y_D2x_b = I_NAI_L6x2y_Px_b+ABX*I_NAI_K5x2y_Px_b;
  Double I_NAI_K5xyz_D2x_b = I_NAI_L6xyz_Px_b+ABX*I_NAI_K5xyz_Px_b;
  Double I_NAI_K5x2z_D2x_b = I_NAI_L6x2z_Px_b+ABX*I_NAI_K5x2z_Px_b;
  Double I_NAI_K4x3y_D2x_b = I_NAI_L5x3y_Px_b+ABX*I_NAI_K4x3y_Px_b;
  Double I_NAI_K4x2yz_D2x_b = I_NAI_L5x2yz_Px_b+ABX*I_NAI_K4x2yz_Px_b;
  Double I_NAI_K4xy2z_D2x_b = I_NAI_L5xy2z_Px_b+ABX*I_NAI_K4xy2z_Px_b;
  Double I_NAI_K4x3z_D2x_b = I_NAI_L5x3z_Px_b+ABX*I_NAI_K4x3z_Px_b;
  Double I_NAI_K3x4y_D2x_b = I_NAI_L4x4y_Px_b+ABX*I_NAI_K3x4y_Px_b;
  Double I_NAI_K3x3yz_D2x_b = I_NAI_L4x3yz_Px_b+ABX*I_NAI_K3x3yz_Px_b;
  Double I_NAI_K3x2y2z_D2x_b = I_NAI_L4x2y2z_Px_b+ABX*I_NAI_K3x2y2z_Px_b;
  Double I_NAI_K3xy3z_D2x_b = I_NAI_L4xy3z_Px_b+ABX*I_NAI_K3xy3z_Px_b;
  Double I_NAI_K3x4z_D2x_b = I_NAI_L4x4z_Px_b+ABX*I_NAI_K3x4z_Px_b;
  Double I_NAI_K2x5y_D2x_b = I_NAI_L3x5y_Px_b+ABX*I_NAI_K2x5y_Px_b;
  Double I_NAI_K2x4yz_D2x_b = I_NAI_L3x4yz_Px_b+ABX*I_NAI_K2x4yz_Px_b;
  Double I_NAI_K2x3y2z_D2x_b = I_NAI_L3x3y2z_Px_b+ABX*I_NAI_K2x3y2z_Px_b;
  Double I_NAI_K2x2y3z_D2x_b = I_NAI_L3x2y3z_Px_b+ABX*I_NAI_K2x2y3z_Px_b;
  Double I_NAI_K2xy4z_D2x_b = I_NAI_L3xy4z_Px_b+ABX*I_NAI_K2xy4z_Px_b;
  Double I_NAI_K2x5z_D2x_b = I_NAI_L3x5z_Px_b+ABX*I_NAI_K2x5z_Px_b;
  Double I_NAI_Kx6y_D2x_b = I_NAI_L2x6y_Px_b+ABX*I_NAI_Kx6y_Px_b;
  Double I_NAI_Kx5yz_D2x_b = I_NAI_L2x5yz_Px_b+ABX*I_NAI_Kx5yz_Px_b;
  Double I_NAI_Kx4y2z_D2x_b = I_NAI_L2x4y2z_Px_b+ABX*I_NAI_Kx4y2z_Px_b;
  Double I_NAI_Kx3y3z_D2x_b = I_NAI_L2x3y3z_Px_b+ABX*I_NAI_Kx3y3z_Px_b;
  Double I_NAI_Kx2y4z_D2x_b = I_NAI_L2x2y4z_Px_b+ABX*I_NAI_Kx2y4z_Px_b;
  Double I_NAI_Kxy5z_D2x_b = I_NAI_L2xy5z_Px_b+ABX*I_NAI_Kxy5z_Px_b;
  Double I_NAI_Kx6z_D2x_b = I_NAI_L2x6z_Px_b+ABX*I_NAI_Kx6z_Px_b;
  Double I_NAI_K7y_D2x_b = I_NAI_Lx7y_Px_b+ABX*I_NAI_K7y_Px_b;
  Double I_NAI_K6yz_D2x_b = I_NAI_Lx6yz_Px_b+ABX*I_NAI_K6yz_Px_b;
  Double I_NAI_K5y2z_D2x_b = I_NAI_Lx5y2z_Px_b+ABX*I_NAI_K5y2z_Px_b;
  Double I_NAI_K4y3z_D2x_b = I_NAI_Lx4y3z_Px_b+ABX*I_NAI_K4y3z_Px_b;
  Double I_NAI_K3y4z_D2x_b = I_NAI_Lx3y4z_Px_b+ABX*I_NAI_K3y4z_Px_b;
  Double I_NAI_K2y5z_D2x_b = I_NAI_Lx2y5z_Px_b+ABX*I_NAI_K2y5z_Px_b;
  Double I_NAI_Ky6z_D2x_b = I_NAI_Lxy6z_Px_b+ABX*I_NAI_Ky6z_Px_b;
  Double I_NAI_K7z_D2x_b = I_NAI_Lx7z_Px_b+ABX*I_NAI_K7z_Px_b;
  Double I_NAI_K6xy_D2y_b = I_NAI_L6x2y_Py_b+ABY*I_NAI_K6xy_Py_b;
  Double I_NAI_K6xz_D2y_b = I_NAI_L6xyz_Py_b+ABY*I_NAI_K6xz_Py_b;
  Double I_NAI_K5x2y_D2y_b = I_NAI_L5x3y_Py_b+ABY*I_NAI_K5x2y_Py_b;
  Double I_NAI_K5xyz_D2y_b = I_NAI_L5x2yz_Py_b+ABY*I_NAI_K5xyz_Py_b;
  Double I_NAI_K5x2z_D2y_b = I_NAI_L5xy2z_Py_b+ABY*I_NAI_K5x2z_Py_b;
  Double I_NAI_K4x3y_D2y_b = I_NAI_L4x4y_Py_b+ABY*I_NAI_K4x3y_Py_b;
  Double I_NAI_K4x2yz_D2y_b = I_NAI_L4x3yz_Py_b+ABY*I_NAI_K4x2yz_Py_b;
  Double I_NAI_K4xy2z_D2y_b = I_NAI_L4x2y2z_Py_b+ABY*I_NAI_K4xy2z_Py_b;
  Double I_NAI_K4x3z_D2y_b = I_NAI_L4xy3z_Py_b+ABY*I_NAI_K4x3z_Py_b;
  Double I_NAI_K3x4y_D2y_b = I_NAI_L3x5y_Py_b+ABY*I_NAI_K3x4y_Py_b;
  Double I_NAI_K3x3yz_D2y_b = I_NAI_L3x4yz_Py_b+ABY*I_NAI_K3x3yz_Py_b;
  Double I_NAI_K3x2y2z_D2y_b = I_NAI_L3x3y2z_Py_b+ABY*I_NAI_K3x2y2z_Py_b;
  Double I_NAI_K3xy3z_D2y_b = I_NAI_L3x2y3z_Py_b+ABY*I_NAI_K3xy3z_Py_b;
  Double I_NAI_K3x4z_D2y_b = I_NAI_L3xy4z_Py_b+ABY*I_NAI_K3x4z_Py_b;
  Double I_NAI_K2x5y_D2y_b = I_NAI_L2x6y_Py_b+ABY*I_NAI_K2x5y_Py_b;
  Double I_NAI_K2x4yz_D2y_b = I_NAI_L2x5yz_Py_b+ABY*I_NAI_K2x4yz_Py_b;
  Double I_NAI_K2x3y2z_D2y_b = I_NAI_L2x4y2z_Py_b+ABY*I_NAI_K2x3y2z_Py_b;
  Double I_NAI_K2x2y3z_D2y_b = I_NAI_L2x3y3z_Py_b+ABY*I_NAI_K2x2y3z_Py_b;
  Double I_NAI_K2xy4z_D2y_b = I_NAI_L2x2y4z_Py_b+ABY*I_NAI_K2xy4z_Py_b;
  Double I_NAI_K2x5z_D2y_b = I_NAI_L2xy5z_Py_b+ABY*I_NAI_K2x5z_Py_b;
  Double I_NAI_Kx6y_D2y_b = I_NAI_Lx7y_Py_b+ABY*I_NAI_Kx6y_Py_b;
  Double I_NAI_Kx5yz_D2y_b = I_NAI_Lx6yz_Py_b+ABY*I_NAI_Kx5yz_Py_b;
  Double I_NAI_Kx4y2z_D2y_b = I_NAI_Lx5y2z_Py_b+ABY*I_NAI_Kx4y2z_Py_b;
  Double I_NAI_Kx3y3z_D2y_b = I_NAI_Lx4y3z_Py_b+ABY*I_NAI_Kx3y3z_Py_b;
  Double I_NAI_Kx2y4z_D2y_b = I_NAI_Lx3y4z_Py_b+ABY*I_NAI_Kx2y4z_Py_b;
  Double I_NAI_Kxy5z_D2y_b = I_NAI_Lx2y5z_Py_b+ABY*I_NAI_Kxy5z_Py_b;
  Double I_NAI_Kx6z_D2y_b = I_NAI_Lxy6z_Py_b+ABY*I_NAI_Kx6z_Py_b;
  Double I_NAI_K7y_D2y_b = I_NAI_L8y_Py_b+ABY*I_NAI_K7y_Py_b;
  Double I_NAI_K6yz_D2y_b = I_NAI_L7yz_Py_b+ABY*I_NAI_K6yz_Py_b;
  Double I_NAI_K5y2z_D2y_b = I_NAI_L6y2z_Py_b+ABY*I_NAI_K5y2z_Py_b;
  Double I_NAI_K4y3z_D2y_b = I_NAI_L5y3z_Py_b+ABY*I_NAI_K4y3z_Py_b;
  Double I_NAI_K3y4z_D2y_b = I_NAI_L4y4z_Py_b+ABY*I_NAI_K3y4z_Py_b;
  Double I_NAI_K2y5z_D2y_b = I_NAI_L3y5z_Py_b+ABY*I_NAI_K2y5z_Py_b;
  Double I_NAI_Ky6z_D2y_b = I_NAI_L2y6z_Py_b+ABY*I_NAI_Ky6z_Py_b;
  Double I_NAI_K7z_D2y_b = I_NAI_Ly7z_Py_b+ABY*I_NAI_K7z_Py_b;
  Double I_NAI_K6xy_D2z_b = I_NAI_L6xyz_Pz_b+ABZ*I_NAI_K6xy_Pz_b;
  Double I_NAI_K6xz_D2z_b = I_NAI_L6x2z_Pz_b+ABZ*I_NAI_K6xz_Pz_b;
  Double I_NAI_K5x2y_D2z_b = I_NAI_L5x2yz_Pz_b+ABZ*I_NAI_K5x2y_Pz_b;
  Double I_NAI_K5xyz_D2z_b = I_NAI_L5xy2z_Pz_b+ABZ*I_NAI_K5xyz_Pz_b;
  Double I_NAI_K5x2z_D2z_b = I_NAI_L5x3z_Pz_b+ABZ*I_NAI_K5x2z_Pz_b;
  Double I_NAI_K4x3y_D2z_b = I_NAI_L4x3yz_Pz_b+ABZ*I_NAI_K4x3y_Pz_b;
  Double I_NAI_K4x2yz_D2z_b = I_NAI_L4x2y2z_Pz_b+ABZ*I_NAI_K4x2yz_Pz_b;
  Double I_NAI_K4xy2z_D2z_b = I_NAI_L4xy3z_Pz_b+ABZ*I_NAI_K4xy2z_Pz_b;
  Double I_NAI_K4x3z_D2z_b = I_NAI_L4x4z_Pz_b+ABZ*I_NAI_K4x3z_Pz_b;
  Double I_NAI_K3x4y_D2z_b = I_NAI_L3x4yz_Pz_b+ABZ*I_NAI_K3x4y_Pz_b;
  Double I_NAI_K3x3yz_D2z_b = I_NAI_L3x3y2z_Pz_b+ABZ*I_NAI_K3x3yz_Pz_b;
  Double I_NAI_K3x2y2z_D2z_b = I_NAI_L3x2y3z_Pz_b+ABZ*I_NAI_K3x2y2z_Pz_b;
  Double I_NAI_K3xy3z_D2z_b = I_NAI_L3xy4z_Pz_b+ABZ*I_NAI_K3xy3z_Pz_b;
  Double I_NAI_K3x4z_D2z_b = I_NAI_L3x5z_Pz_b+ABZ*I_NAI_K3x4z_Pz_b;
  Double I_NAI_K2x5y_D2z_b = I_NAI_L2x5yz_Pz_b+ABZ*I_NAI_K2x5y_Pz_b;
  Double I_NAI_K2x4yz_D2z_b = I_NAI_L2x4y2z_Pz_b+ABZ*I_NAI_K2x4yz_Pz_b;
  Double I_NAI_K2x3y2z_D2z_b = I_NAI_L2x3y3z_Pz_b+ABZ*I_NAI_K2x3y2z_Pz_b;
  Double I_NAI_K2x2y3z_D2z_b = I_NAI_L2x2y4z_Pz_b+ABZ*I_NAI_K2x2y3z_Pz_b;
  Double I_NAI_K2xy4z_D2z_b = I_NAI_L2xy5z_Pz_b+ABZ*I_NAI_K2xy4z_Pz_b;
  Double I_NAI_K2x5z_D2z_b = I_NAI_L2x6z_Pz_b+ABZ*I_NAI_K2x5z_Pz_b;
  Double I_NAI_Kx6y_D2z_b = I_NAI_Lx6yz_Pz_b+ABZ*I_NAI_Kx6y_Pz_b;
  Double I_NAI_Kx5yz_D2z_b = I_NAI_Lx5y2z_Pz_b+ABZ*I_NAI_Kx5yz_Pz_b;
  Double I_NAI_Kx4y2z_D2z_b = I_NAI_Lx4y3z_Pz_b+ABZ*I_NAI_Kx4y2z_Pz_b;
  Double I_NAI_Kx3y3z_D2z_b = I_NAI_Lx3y4z_Pz_b+ABZ*I_NAI_Kx3y3z_Pz_b;
  Double I_NAI_Kx2y4z_D2z_b = I_NAI_Lx2y5z_Pz_b+ABZ*I_NAI_Kx2y4z_Pz_b;
  Double I_NAI_Kxy5z_D2z_b = I_NAI_Lxy6z_Pz_b+ABZ*I_NAI_Kxy5z_Pz_b;
  Double I_NAI_Kx6z_D2z_b = I_NAI_Lx7z_Pz_b+ABZ*I_NAI_Kx6z_Pz_b;
  Double I_NAI_K6yz_D2z_b = I_NAI_L6y2z_Pz_b+ABZ*I_NAI_K6yz_Pz_b;
  Double I_NAI_K5y2z_D2z_b = I_NAI_L5y3z_Pz_b+ABZ*I_NAI_K5y2z_Pz_b;
  Double I_NAI_K4y3z_D2z_b = I_NAI_L4y4z_Pz_b+ABZ*I_NAI_K4y3z_Pz_b;
  Double I_NAI_K3y4z_D2z_b = I_NAI_L3y5z_Pz_b+ABZ*I_NAI_K3y4z_Pz_b;
  Double I_NAI_K2y5z_D2z_b = I_NAI_L2y6z_Pz_b+ABZ*I_NAI_K2y5z_Pz_b;
  Double I_NAI_Ky6z_D2z_b = I_NAI_Ly7z_Pz_b+ABZ*I_NAI_Ky6z_Pz_b;
  Double I_NAI_K7z_D2z_b = I_NAI_L8z_Pz_b+ABZ*I_NAI_K7z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 85 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_D_b
   * RHS shell quartet name: SQ_NAI_I_D_b
   ************************************************************/
  Double I_NAI_I6x_F3x_b = I_NAI_K7x_D2x_b+ABX*I_NAI_I6x_D2x_b;
  Double I_NAI_I5xy_F3x_b = I_NAI_K6xy_D2x_b+ABX*I_NAI_I5xy_D2x_b;
  Double I_NAI_I5xz_F3x_b = I_NAI_K6xz_D2x_b+ABX*I_NAI_I5xz_D2x_b;
  Double I_NAI_I4x2y_F3x_b = I_NAI_K5x2y_D2x_b+ABX*I_NAI_I4x2y_D2x_b;
  Double I_NAI_I4xyz_F3x_b = I_NAI_K5xyz_D2x_b+ABX*I_NAI_I4xyz_D2x_b;
  Double I_NAI_I4x2z_F3x_b = I_NAI_K5x2z_D2x_b+ABX*I_NAI_I4x2z_D2x_b;
  Double I_NAI_I3x3y_F3x_b = I_NAI_K4x3y_D2x_b+ABX*I_NAI_I3x3y_D2x_b;
  Double I_NAI_I3x2yz_F3x_b = I_NAI_K4x2yz_D2x_b+ABX*I_NAI_I3x2yz_D2x_b;
  Double I_NAI_I3xy2z_F3x_b = I_NAI_K4xy2z_D2x_b+ABX*I_NAI_I3xy2z_D2x_b;
  Double I_NAI_I3x3z_F3x_b = I_NAI_K4x3z_D2x_b+ABX*I_NAI_I3x3z_D2x_b;
  Double I_NAI_I2x4y_F3x_b = I_NAI_K3x4y_D2x_b+ABX*I_NAI_I2x4y_D2x_b;
  Double I_NAI_I2x3yz_F3x_b = I_NAI_K3x3yz_D2x_b+ABX*I_NAI_I2x3yz_D2x_b;
  Double I_NAI_I2x2y2z_F3x_b = I_NAI_K3x2y2z_D2x_b+ABX*I_NAI_I2x2y2z_D2x_b;
  Double I_NAI_I2xy3z_F3x_b = I_NAI_K3xy3z_D2x_b+ABX*I_NAI_I2xy3z_D2x_b;
  Double I_NAI_I2x4z_F3x_b = I_NAI_K3x4z_D2x_b+ABX*I_NAI_I2x4z_D2x_b;
  Double I_NAI_Ix5y_F3x_b = I_NAI_K2x5y_D2x_b+ABX*I_NAI_Ix5y_D2x_b;
  Double I_NAI_Ix4yz_F3x_b = I_NAI_K2x4yz_D2x_b+ABX*I_NAI_Ix4yz_D2x_b;
  Double I_NAI_Ix3y2z_F3x_b = I_NAI_K2x3y2z_D2x_b+ABX*I_NAI_Ix3y2z_D2x_b;
  Double I_NAI_Ix2y3z_F3x_b = I_NAI_K2x2y3z_D2x_b+ABX*I_NAI_Ix2y3z_D2x_b;
  Double I_NAI_Ixy4z_F3x_b = I_NAI_K2xy4z_D2x_b+ABX*I_NAI_Ixy4z_D2x_b;
  Double I_NAI_Ix5z_F3x_b = I_NAI_K2x5z_D2x_b+ABX*I_NAI_Ix5z_D2x_b;
  Double I_NAI_I6y_F3x_b = I_NAI_Kx6y_D2x_b+ABX*I_NAI_I6y_D2x_b;
  Double I_NAI_I5yz_F3x_b = I_NAI_Kx5yz_D2x_b+ABX*I_NAI_I5yz_D2x_b;
  Double I_NAI_I4y2z_F3x_b = I_NAI_Kx4y2z_D2x_b+ABX*I_NAI_I4y2z_D2x_b;
  Double I_NAI_I3y3z_F3x_b = I_NAI_Kx3y3z_D2x_b+ABX*I_NAI_I3y3z_D2x_b;
  Double I_NAI_I2y4z_F3x_b = I_NAI_Kx2y4z_D2x_b+ABX*I_NAI_I2y4z_D2x_b;
  Double I_NAI_Iy5z_F3x_b = I_NAI_Kxy5z_D2x_b+ABX*I_NAI_Iy5z_D2x_b;
  Double I_NAI_I6z_F3x_b = I_NAI_Kx6z_D2x_b+ABX*I_NAI_I6z_D2x_b;
  Double I_NAI_I5xy_F2xy_b = I_NAI_K5x2y_D2x_b+ABY*I_NAI_I5xy_D2x_b;
  Double I_NAI_I5xz_F2xy_b = I_NAI_K5xyz_D2x_b+ABY*I_NAI_I5xz_D2x_b;
  Double I_NAI_I4x2y_F2xy_b = I_NAI_K4x3y_D2x_b+ABY*I_NAI_I4x2y_D2x_b;
  Double I_NAI_I4xyz_F2xy_b = I_NAI_K4x2yz_D2x_b+ABY*I_NAI_I4xyz_D2x_b;
  Double I_NAI_I4x2z_F2xy_b = I_NAI_K4xy2z_D2x_b+ABY*I_NAI_I4x2z_D2x_b;
  Double I_NAI_I3x3y_F2xy_b = I_NAI_K3x4y_D2x_b+ABY*I_NAI_I3x3y_D2x_b;
  Double I_NAI_I3x2yz_F2xy_b = I_NAI_K3x3yz_D2x_b+ABY*I_NAI_I3x2yz_D2x_b;
  Double I_NAI_I3xy2z_F2xy_b = I_NAI_K3x2y2z_D2x_b+ABY*I_NAI_I3xy2z_D2x_b;
  Double I_NAI_I3x3z_F2xy_b = I_NAI_K3xy3z_D2x_b+ABY*I_NAI_I3x3z_D2x_b;
  Double I_NAI_I2x4y_F2xy_b = I_NAI_K2x5y_D2x_b+ABY*I_NAI_I2x4y_D2x_b;
  Double I_NAI_I2x3yz_F2xy_b = I_NAI_K2x4yz_D2x_b+ABY*I_NAI_I2x3yz_D2x_b;
  Double I_NAI_I2x2y2z_F2xy_b = I_NAI_K2x3y2z_D2x_b+ABY*I_NAI_I2x2y2z_D2x_b;
  Double I_NAI_I2xy3z_F2xy_b = I_NAI_K2x2y3z_D2x_b+ABY*I_NAI_I2xy3z_D2x_b;
  Double I_NAI_I2x4z_F2xy_b = I_NAI_K2xy4z_D2x_b+ABY*I_NAI_I2x4z_D2x_b;
  Double I_NAI_Ix5y_F2xy_b = I_NAI_Kx6y_D2x_b+ABY*I_NAI_Ix5y_D2x_b;
  Double I_NAI_Ix4yz_F2xy_b = I_NAI_Kx5yz_D2x_b+ABY*I_NAI_Ix4yz_D2x_b;
  Double I_NAI_Ix3y2z_F2xy_b = I_NAI_Kx4y2z_D2x_b+ABY*I_NAI_Ix3y2z_D2x_b;
  Double I_NAI_Ix2y3z_F2xy_b = I_NAI_Kx3y3z_D2x_b+ABY*I_NAI_Ix2y3z_D2x_b;
  Double I_NAI_Ixy4z_F2xy_b = I_NAI_Kx2y4z_D2x_b+ABY*I_NAI_Ixy4z_D2x_b;
  Double I_NAI_Ix5z_F2xy_b = I_NAI_Kxy5z_D2x_b+ABY*I_NAI_Ix5z_D2x_b;
  Double I_NAI_I6y_F2xy_b = I_NAI_K7y_D2x_b+ABY*I_NAI_I6y_D2x_b;
  Double I_NAI_I5yz_F2xy_b = I_NAI_K6yz_D2x_b+ABY*I_NAI_I5yz_D2x_b;
  Double I_NAI_I4y2z_F2xy_b = I_NAI_K5y2z_D2x_b+ABY*I_NAI_I4y2z_D2x_b;
  Double I_NAI_I3y3z_F2xy_b = I_NAI_K4y3z_D2x_b+ABY*I_NAI_I3y3z_D2x_b;
  Double I_NAI_I2y4z_F2xy_b = I_NAI_K3y4z_D2x_b+ABY*I_NAI_I2y4z_D2x_b;
  Double I_NAI_Iy5z_F2xy_b = I_NAI_K2y5z_D2x_b+ABY*I_NAI_Iy5z_D2x_b;
  Double I_NAI_I6z_F2xy_b = I_NAI_Ky6z_D2x_b+ABY*I_NAI_I6z_D2x_b;
  Double I_NAI_I5xz_F2xz_b = I_NAI_K5x2z_D2x_b+ABZ*I_NAI_I5xz_D2x_b;
  Double I_NAI_I4xyz_F2xz_b = I_NAI_K4xy2z_D2x_b+ABZ*I_NAI_I4xyz_D2x_b;
  Double I_NAI_I4x2z_F2xz_b = I_NAI_K4x3z_D2x_b+ABZ*I_NAI_I4x2z_D2x_b;
  Double I_NAI_I3x2yz_F2xz_b = I_NAI_K3x2y2z_D2x_b+ABZ*I_NAI_I3x2yz_D2x_b;
  Double I_NAI_I3xy2z_F2xz_b = I_NAI_K3xy3z_D2x_b+ABZ*I_NAI_I3xy2z_D2x_b;
  Double I_NAI_I3x3z_F2xz_b = I_NAI_K3x4z_D2x_b+ABZ*I_NAI_I3x3z_D2x_b;
  Double I_NAI_I2x3yz_F2xz_b = I_NAI_K2x3y2z_D2x_b+ABZ*I_NAI_I2x3yz_D2x_b;
  Double I_NAI_I2x2y2z_F2xz_b = I_NAI_K2x2y3z_D2x_b+ABZ*I_NAI_I2x2y2z_D2x_b;
  Double I_NAI_I2xy3z_F2xz_b = I_NAI_K2xy4z_D2x_b+ABZ*I_NAI_I2xy3z_D2x_b;
  Double I_NAI_I2x4z_F2xz_b = I_NAI_K2x5z_D2x_b+ABZ*I_NAI_I2x4z_D2x_b;
  Double I_NAI_Ix4yz_F2xz_b = I_NAI_Kx4y2z_D2x_b+ABZ*I_NAI_Ix4yz_D2x_b;
  Double I_NAI_Ix3y2z_F2xz_b = I_NAI_Kx3y3z_D2x_b+ABZ*I_NAI_Ix3y2z_D2x_b;
  Double I_NAI_Ix2y3z_F2xz_b = I_NAI_Kx2y4z_D2x_b+ABZ*I_NAI_Ix2y3z_D2x_b;
  Double I_NAI_Ixy4z_F2xz_b = I_NAI_Kxy5z_D2x_b+ABZ*I_NAI_Ixy4z_D2x_b;
  Double I_NAI_Ix5z_F2xz_b = I_NAI_Kx6z_D2x_b+ABZ*I_NAI_Ix5z_D2x_b;
  Double I_NAI_I5yz_F2xz_b = I_NAI_K5y2z_D2x_b+ABZ*I_NAI_I5yz_D2x_b;
  Double I_NAI_I4y2z_F2xz_b = I_NAI_K4y3z_D2x_b+ABZ*I_NAI_I4y2z_D2x_b;
  Double I_NAI_I3y3z_F2xz_b = I_NAI_K3y4z_D2x_b+ABZ*I_NAI_I3y3z_D2x_b;
  Double I_NAI_I2y4z_F2xz_b = I_NAI_K2y5z_D2x_b+ABZ*I_NAI_I2y4z_D2x_b;
  Double I_NAI_Iy5z_F2xz_b = I_NAI_Ky6z_D2x_b+ABZ*I_NAI_Iy5z_D2x_b;
  Double I_NAI_I6z_F2xz_b = I_NAI_K7z_D2x_b+ABZ*I_NAI_I6z_D2x_b;
  Double I_NAI_I5xz_Fx2y_b = I_NAI_K6xz_D2y_b+ABX*I_NAI_I5xz_D2y_b;
  Double I_NAI_I4xyz_Fx2y_b = I_NAI_K5xyz_D2y_b+ABX*I_NAI_I4xyz_D2y_b;
  Double I_NAI_I4x2z_Fx2y_b = I_NAI_K5x2z_D2y_b+ABX*I_NAI_I4x2z_D2y_b;
  Double I_NAI_I3x2yz_Fx2y_b = I_NAI_K4x2yz_D2y_b+ABX*I_NAI_I3x2yz_D2y_b;
  Double I_NAI_I3xy2z_Fx2y_b = I_NAI_K4xy2z_D2y_b+ABX*I_NAI_I3xy2z_D2y_b;
  Double I_NAI_I3x3z_Fx2y_b = I_NAI_K4x3z_D2y_b+ABX*I_NAI_I3x3z_D2y_b;
  Double I_NAI_I2x3yz_Fx2y_b = I_NAI_K3x3yz_D2y_b+ABX*I_NAI_I2x3yz_D2y_b;
  Double I_NAI_I2x2y2z_Fx2y_b = I_NAI_K3x2y2z_D2y_b+ABX*I_NAI_I2x2y2z_D2y_b;
  Double I_NAI_I2xy3z_Fx2y_b = I_NAI_K3xy3z_D2y_b+ABX*I_NAI_I2xy3z_D2y_b;
  Double I_NAI_I2x4z_Fx2y_b = I_NAI_K3x4z_D2y_b+ABX*I_NAI_I2x4z_D2y_b;
  Double I_NAI_Ix4yz_Fx2y_b = I_NAI_K2x4yz_D2y_b+ABX*I_NAI_Ix4yz_D2y_b;
  Double I_NAI_Ix3y2z_Fx2y_b = I_NAI_K2x3y2z_D2y_b+ABX*I_NAI_Ix3y2z_D2y_b;
  Double I_NAI_Ix2y3z_Fx2y_b = I_NAI_K2x2y3z_D2y_b+ABX*I_NAI_Ix2y3z_D2y_b;
  Double I_NAI_Ixy4z_Fx2y_b = I_NAI_K2xy4z_D2y_b+ABX*I_NAI_Ixy4z_D2y_b;
  Double I_NAI_Ix5z_Fx2y_b = I_NAI_K2x5z_D2y_b+ABX*I_NAI_Ix5z_D2y_b;
  Double I_NAI_I5yz_Fx2y_b = I_NAI_Kx5yz_D2y_b+ABX*I_NAI_I5yz_D2y_b;
  Double I_NAI_I4y2z_Fx2y_b = I_NAI_Kx4y2z_D2y_b+ABX*I_NAI_I4y2z_D2y_b;
  Double I_NAI_I3y3z_Fx2y_b = I_NAI_Kx3y3z_D2y_b+ABX*I_NAI_I3y3z_D2y_b;
  Double I_NAI_I2y4z_Fx2y_b = I_NAI_Kx2y4z_D2y_b+ABX*I_NAI_I2y4z_D2y_b;
  Double I_NAI_Iy5z_Fx2y_b = I_NAI_Kxy5z_D2y_b+ABX*I_NAI_Iy5z_D2y_b;
  Double I_NAI_I6z_Fx2y_b = I_NAI_Kx6z_D2y_b+ABX*I_NAI_I6z_D2y_b;
  Double I_NAI_I5xy_Fx2z_b = I_NAI_K6xy_D2z_b+ABX*I_NAI_I5xy_D2z_b;
  Double I_NAI_I4x2y_Fx2z_b = I_NAI_K5x2y_D2z_b+ABX*I_NAI_I4x2y_D2z_b;
  Double I_NAI_I4xyz_Fx2z_b = I_NAI_K5xyz_D2z_b+ABX*I_NAI_I4xyz_D2z_b;
  Double I_NAI_I3x3y_Fx2z_b = I_NAI_K4x3y_D2z_b+ABX*I_NAI_I3x3y_D2z_b;
  Double I_NAI_I3x2yz_Fx2z_b = I_NAI_K4x2yz_D2z_b+ABX*I_NAI_I3x2yz_D2z_b;
  Double I_NAI_I3xy2z_Fx2z_b = I_NAI_K4xy2z_D2z_b+ABX*I_NAI_I3xy2z_D2z_b;
  Double I_NAI_I2x4y_Fx2z_b = I_NAI_K3x4y_D2z_b+ABX*I_NAI_I2x4y_D2z_b;
  Double I_NAI_I2x3yz_Fx2z_b = I_NAI_K3x3yz_D2z_b+ABX*I_NAI_I2x3yz_D2z_b;
  Double I_NAI_I2x2y2z_Fx2z_b = I_NAI_K3x2y2z_D2z_b+ABX*I_NAI_I2x2y2z_D2z_b;
  Double I_NAI_I2xy3z_Fx2z_b = I_NAI_K3xy3z_D2z_b+ABX*I_NAI_I2xy3z_D2z_b;
  Double I_NAI_Ix5y_Fx2z_b = I_NAI_K2x5y_D2z_b+ABX*I_NAI_Ix5y_D2z_b;
  Double I_NAI_Ix4yz_Fx2z_b = I_NAI_K2x4yz_D2z_b+ABX*I_NAI_Ix4yz_D2z_b;
  Double I_NAI_Ix3y2z_Fx2z_b = I_NAI_K2x3y2z_D2z_b+ABX*I_NAI_Ix3y2z_D2z_b;
  Double I_NAI_Ix2y3z_Fx2z_b = I_NAI_K2x2y3z_D2z_b+ABX*I_NAI_Ix2y3z_D2z_b;
  Double I_NAI_Ixy4z_Fx2z_b = I_NAI_K2xy4z_D2z_b+ABX*I_NAI_Ixy4z_D2z_b;
  Double I_NAI_I6y_Fx2z_b = I_NAI_Kx6y_D2z_b+ABX*I_NAI_I6y_D2z_b;
  Double I_NAI_I5yz_Fx2z_b = I_NAI_Kx5yz_D2z_b+ABX*I_NAI_I5yz_D2z_b;
  Double I_NAI_I4y2z_Fx2z_b = I_NAI_Kx4y2z_D2z_b+ABX*I_NAI_I4y2z_D2z_b;
  Double I_NAI_I3y3z_Fx2z_b = I_NAI_Kx3y3z_D2z_b+ABX*I_NAI_I3y3z_D2z_b;
  Double I_NAI_I2y4z_Fx2z_b = I_NAI_Kx2y4z_D2z_b+ABX*I_NAI_I2y4z_D2z_b;
  Double I_NAI_Iy5z_Fx2z_b = I_NAI_Kxy5z_D2z_b+ABX*I_NAI_Iy5z_D2z_b;
  Double I_NAI_I6x_F3y_b = I_NAI_K6xy_D2y_b+ABY*I_NAI_I6x_D2y_b;
  Double I_NAI_I5xy_F3y_b = I_NAI_K5x2y_D2y_b+ABY*I_NAI_I5xy_D2y_b;
  Double I_NAI_I5xz_F3y_b = I_NAI_K5xyz_D2y_b+ABY*I_NAI_I5xz_D2y_b;
  Double I_NAI_I4x2y_F3y_b = I_NAI_K4x3y_D2y_b+ABY*I_NAI_I4x2y_D2y_b;
  Double I_NAI_I4xyz_F3y_b = I_NAI_K4x2yz_D2y_b+ABY*I_NAI_I4xyz_D2y_b;
  Double I_NAI_I4x2z_F3y_b = I_NAI_K4xy2z_D2y_b+ABY*I_NAI_I4x2z_D2y_b;
  Double I_NAI_I3x3y_F3y_b = I_NAI_K3x4y_D2y_b+ABY*I_NAI_I3x3y_D2y_b;
  Double I_NAI_I3x2yz_F3y_b = I_NAI_K3x3yz_D2y_b+ABY*I_NAI_I3x2yz_D2y_b;
  Double I_NAI_I3xy2z_F3y_b = I_NAI_K3x2y2z_D2y_b+ABY*I_NAI_I3xy2z_D2y_b;
  Double I_NAI_I3x3z_F3y_b = I_NAI_K3xy3z_D2y_b+ABY*I_NAI_I3x3z_D2y_b;
  Double I_NAI_I2x4y_F3y_b = I_NAI_K2x5y_D2y_b+ABY*I_NAI_I2x4y_D2y_b;
  Double I_NAI_I2x3yz_F3y_b = I_NAI_K2x4yz_D2y_b+ABY*I_NAI_I2x3yz_D2y_b;
  Double I_NAI_I2x2y2z_F3y_b = I_NAI_K2x3y2z_D2y_b+ABY*I_NAI_I2x2y2z_D2y_b;
  Double I_NAI_I2xy3z_F3y_b = I_NAI_K2x2y3z_D2y_b+ABY*I_NAI_I2xy3z_D2y_b;
  Double I_NAI_I2x4z_F3y_b = I_NAI_K2xy4z_D2y_b+ABY*I_NAI_I2x4z_D2y_b;
  Double I_NAI_Ix5y_F3y_b = I_NAI_Kx6y_D2y_b+ABY*I_NAI_Ix5y_D2y_b;
  Double I_NAI_Ix4yz_F3y_b = I_NAI_Kx5yz_D2y_b+ABY*I_NAI_Ix4yz_D2y_b;
  Double I_NAI_Ix3y2z_F3y_b = I_NAI_Kx4y2z_D2y_b+ABY*I_NAI_Ix3y2z_D2y_b;
  Double I_NAI_Ix2y3z_F3y_b = I_NAI_Kx3y3z_D2y_b+ABY*I_NAI_Ix2y3z_D2y_b;
  Double I_NAI_Ixy4z_F3y_b = I_NAI_Kx2y4z_D2y_b+ABY*I_NAI_Ixy4z_D2y_b;
  Double I_NAI_Ix5z_F3y_b = I_NAI_Kxy5z_D2y_b+ABY*I_NAI_Ix5z_D2y_b;
  Double I_NAI_I6y_F3y_b = I_NAI_K7y_D2y_b+ABY*I_NAI_I6y_D2y_b;
  Double I_NAI_I5yz_F3y_b = I_NAI_K6yz_D2y_b+ABY*I_NAI_I5yz_D2y_b;
  Double I_NAI_I4y2z_F3y_b = I_NAI_K5y2z_D2y_b+ABY*I_NAI_I4y2z_D2y_b;
  Double I_NAI_I3y3z_F3y_b = I_NAI_K4y3z_D2y_b+ABY*I_NAI_I3y3z_D2y_b;
  Double I_NAI_I2y4z_F3y_b = I_NAI_K3y4z_D2y_b+ABY*I_NAI_I2y4z_D2y_b;
  Double I_NAI_Iy5z_F3y_b = I_NAI_K2y5z_D2y_b+ABY*I_NAI_Iy5z_D2y_b;
  Double I_NAI_I6z_F3y_b = I_NAI_Ky6z_D2y_b+ABY*I_NAI_I6z_D2y_b;
  Double I_NAI_I5xz_F2yz_b = I_NAI_K5x2z_D2y_b+ABZ*I_NAI_I5xz_D2y_b;
  Double I_NAI_I4xyz_F2yz_b = I_NAI_K4xy2z_D2y_b+ABZ*I_NAI_I4xyz_D2y_b;
  Double I_NAI_I4x2z_F2yz_b = I_NAI_K4x3z_D2y_b+ABZ*I_NAI_I4x2z_D2y_b;
  Double I_NAI_I3x2yz_F2yz_b = I_NAI_K3x2y2z_D2y_b+ABZ*I_NAI_I3x2yz_D2y_b;
  Double I_NAI_I3xy2z_F2yz_b = I_NAI_K3xy3z_D2y_b+ABZ*I_NAI_I3xy2z_D2y_b;
  Double I_NAI_I3x3z_F2yz_b = I_NAI_K3x4z_D2y_b+ABZ*I_NAI_I3x3z_D2y_b;
  Double I_NAI_I2x3yz_F2yz_b = I_NAI_K2x3y2z_D2y_b+ABZ*I_NAI_I2x3yz_D2y_b;
  Double I_NAI_I2x2y2z_F2yz_b = I_NAI_K2x2y3z_D2y_b+ABZ*I_NAI_I2x2y2z_D2y_b;
  Double I_NAI_I2xy3z_F2yz_b = I_NAI_K2xy4z_D2y_b+ABZ*I_NAI_I2xy3z_D2y_b;
  Double I_NAI_I2x4z_F2yz_b = I_NAI_K2x5z_D2y_b+ABZ*I_NAI_I2x4z_D2y_b;
  Double I_NAI_Ix4yz_F2yz_b = I_NAI_Kx4y2z_D2y_b+ABZ*I_NAI_Ix4yz_D2y_b;
  Double I_NAI_Ix3y2z_F2yz_b = I_NAI_Kx3y3z_D2y_b+ABZ*I_NAI_Ix3y2z_D2y_b;
  Double I_NAI_Ix2y3z_F2yz_b = I_NAI_Kx2y4z_D2y_b+ABZ*I_NAI_Ix2y3z_D2y_b;
  Double I_NAI_Ixy4z_F2yz_b = I_NAI_Kxy5z_D2y_b+ABZ*I_NAI_Ixy4z_D2y_b;
  Double I_NAI_Ix5z_F2yz_b = I_NAI_Kx6z_D2y_b+ABZ*I_NAI_Ix5z_D2y_b;
  Double I_NAI_I5yz_F2yz_b = I_NAI_K5y2z_D2y_b+ABZ*I_NAI_I5yz_D2y_b;
  Double I_NAI_I4y2z_F2yz_b = I_NAI_K4y3z_D2y_b+ABZ*I_NAI_I4y2z_D2y_b;
  Double I_NAI_I3y3z_F2yz_b = I_NAI_K3y4z_D2y_b+ABZ*I_NAI_I3y3z_D2y_b;
  Double I_NAI_I2y4z_F2yz_b = I_NAI_K2y5z_D2y_b+ABZ*I_NAI_I2y4z_D2y_b;
  Double I_NAI_Iy5z_F2yz_b = I_NAI_Ky6z_D2y_b+ABZ*I_NAI_Iy5z_D2y_b;
  Double I_NAI_I6z_F2yz_b = I_NAI_K7z_D2y_b+ABZ*I_NAI_I6z_D2y_b;
  Double I_NAI_I6x_F3z_b = I_NAI_K6xz_D2z_b+ABZ*I_NAI_I6x_D2z_b;
  Double I_NAI_I5xy_F3z_b = I_NAI_K5xyz_D2z_b+ABZ*I_NAI_I5xy_D2z_b;
  Double I_NAI_I5xz_F3z_b = I_NAI_K5x2z_D2z_b+ABZ*I_NAI_I5xz_D2z_b;
  Double I_NAI_I4x2y_F3z_b = I_NAI_K4x2yz_D2z_b+ABZ*I_NAI_I4x2y_D2z_b;
  Double I_NAI_I4xyz_F3z_b = I_NAI_K4xy2z_D2z_b+ABZ*I_NAI_I4xyz_D2z_b;
  Double I_NAI_I4x2z_F3z_b = I_NAI_K4x3z_D2z_b+ABZ*I_NAI_I4x2z_D2z_b;
  Double I_NAI_I3x3y_F3z_b = I_NAI_K3x3yz_D2z_b+ABZ*I_NAI_I3x3y_D2z_b;
  Double I_NAI_I3x2yz_F3z_b = I_NAI_K3x2y2z_D2z_b+ABZ*I_NAI_I3x2yz_D2z_b;
  Double I_NAI_I3xy2z_F3z_b = I_NAI_K3xy3z_D2z_b+ABZ*I_NAI_I3xy2z_D2z_b;
  Double I_NAI_I3x3z_F3z_b = I_NAI_K3x4z_D2z_b+ABZ*I_NAI_I3x3z_D2z_b;
  Double I_NAI_I2x4y_F3z_b = I_NAI_K2x4yz_D2z_b+ABZ*I_NAI_I2x4y_D2z_b;
  Double I_NAI_I2x3yz_F3z_b = I_NAI_K2x3y2z_D2z_b+ABZ*I_NAI_I2x3yz_D2z_b;
  Double I_NAI_I2x2y2z_F3z_b = I_NAI_K2x2y3z_D2z_b+ABZ*I_NAI_I2x2y2z_D2z_b;
  Double I_NAI_I2xy3z_F3z_b = I_NAI_K2xy4z_D2z_b+ABZ*I_NAI_I2xy3z_D2z_b;
  Double I_NAI_I2x4z_F3z_b = I_NAI_K2x5z_D2z_b+ABZ*I_NAI_I2x4z_D2z_b;
  Double I_NAI_Ix5y_F3z_b = I_NAI_Kx5yz_D2z_b+ABZ*I_NAI_Ix5y_D2z_b;
  Double I_NAI_Ix4yz_F3z_b = I_NAI_Kx4y2z_D2z_b+ABZ*I_NAI_Ix4yz_D2z_b;
  Double I_NAI_Ix3y2z_F3z_b = I_NAI_Kx3y3z_D2z_b+ABZ*I_NAI_Ix3y2z_D2z_b;
  Double I_NAI_Ix2y3z_F3z_b = I_NAI_Kx2y4z_D2z_b+ABZ*I_NAI_Ix2y3z_D2z_b;
  Double I_NAI_Ixy4z_F3z_b = I_NAI_Kxy5z_D2z_b+ABZ*I_NAI_Ixy4z_D2z_b;
  Double I_NAI_Ix5z_F3z_b = I_NAI_Kx6z_D2z_b+ABZ*I_NAI_Ix5z_D2z_b;
  Double I_NAI_I6y_F3z_b = I_NAI_K6yz_D2z_b+ABZ*I_NAI_I6y_D2z_b;
  Double I_NAI_I5yz_F3z_b = I_NAI_K5y2z_D2z_b+ABZ*I_NAI_I5yz_D2z_b;
  Double I_NAI_I4y2z_F3z_b = I_NAI_K4y3z_D2z_b+ABZ*I_NAI_I4y2z_D2z_b;
  Double I_NAI_I3y3z_F3z_b = I_NAI_K3y4z_D2z_b+ABZ*I_NAI_I3y3z_D2z_b;
  Double I_NAI_I2y4z_F3z_b = I_NAI_K2y5z_D2z_b+ABZ*I_NAI_I2y4z_D2z_b;
  Double I_NAI_Iy5z_F3z_b = I_NAI_Ky6z_D2z_b+ABZ*I_NAI_Iy5z_D2z_b;
  Double I_NAI_I6z_F3z_b = I_NAI_K7z_D2z_b+ABZ*I_NAI_I6z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_b
   * RHS shell quartet name: SQ_NAI_H_F_b
   ************************************************************/
  Double I_NAI_H5x_G4x_b = I_NAI_I6x_F3x_b+ABX*I_NAI_H5x_F3x_b;
  Double I_NAI_H4xy_G4x_b = I_NAI_I5xy_F3x_b+ABX*I_NAI_H4xy_F3x_b;
  Double I_NAI_H4xz_G4x_b = I_NAI_I5xz_F3x_b+ABX*I_NAI_H4xz_F3x_b;
  Double I_NAI_H3x2y_G4x_b = I_NAI_I4x2y_F3x_b+ABX*I_NAI_H3x2y_F3x_b;
  Double I_NAI_H3xyz_G4x_b = I_NAI_I4xyz_F3x_b+ABX*I_NAI_H3xyz_F3x_b;
  Double I_NAI_H3x2z_G4x_b = I_NAI_I4x2z_F3x_b+ABX*I_NAI_H3x2z_F3x_b;
  Double I_NAI_H2x3y_G4x_b = I_NAI_I3x3y_F3x_b+ABX*I_NAI_H2x3y_F3x_b;
  Double I_NAI_H2x2yz_G4x_b = I_NAI_I3x2yz_F3x_b+ABX*I_NAI_H2x2yz_F3x_b;
  Double I_NAI_H2xy2z_G4x_b = I_NAI_I3xy2z_F3x_b+ABX*I_NAI_H2xy2z_F3x_b;
  Double I_NAI_H2x3z_G4x_b = I_NAI_I3x3z_F3x_b+ABX*I_NAI_H2x3z_F3x_b;
  Double I_NAI_Hx4y_G4x_b = I_NAI_I2x4y_F3x_b+ABX*I_NAI_Hx4y_F3x_b;
  Double I_NAI_Hx3yz_G4x_b = I_NAI_I2x3yz_F3x_b+ABX*I_NAI_Hx3yz_F3x_b;
  Double I_NAI_Hx2y2z_G4x_b = I_NAI_I2x2y2z_F3x_b+ABX*I_NAI_Hx2y2z_F3x_b;
  Double I_NAI_Hxy3z_G4x_b = I_NAI_I2xy3z_F3x_b+ABX*I_NAI_Hxy3z_F3x_b;
  Double I_NAI_Hx4z_G4x_b = I_NAI_I2x4z_F3x_b+ABX*I_NAI_Hx4z_F3x_b;
  Double I_NAI_H5y_G4x_b = I_NAI_Ix5y_F3x_b+ABX*I_NAI_H5y_F3x_b;
  Double I_NAI_H4yz_G4x_b = I_NAI_Ix4yz_F3x_b+ABX*I_NAI_H4yz_F3x_b;
  Double I_NAI_H3y2z_G4x_b = I_NAI_Ix3y2z_F3x_b+ABX*I_NAI_H3y2z_F3x_b;
  Double I_NAI_H2y3z_G4x_b = I_NAI_Ix2y3z_F3x_b+ABX*I_NAI_H2y3z_F3x_b;
  Double I_NAI_Hy4z_G4x_b = I_NAI_Ixy4z_F3x_b+ABX*I_NAI_Hy4z_F3x_b;
  Double I_NAI_H5z_G4x_b = I_NAI_Ix5z_F3x_b+ABX*I_NAI_H5z_F3x_b;
  Double I_NAI_H5x_G3xy_b = I_NAI_I5xy_F3x_b+ABY*I_NAI_H5x_F3x_b;
  Double I_NAI_H4xy_G3xy_b = I_NAI_I4x2y_F3x_b+ABY*I_NAI_H4xy_F3x_b;
  Double I_NAI_H4xz_G3xy_b = I_NAI_I4xyz_F3x_b+ABY*I_NAI_H4xz_F3x_b;
  Double I_NAI_H3x2y_G3xy_b = I_NAI_I3x3y_F3x_b+ABY*I_NAI_H3x2y_F3x_b;
  Double I_NAI_H3xyz_G3xy_b = I_NAI_I3x2yz_F3x_b+ABY*I_NAI_H3xyz_F3x_b;
  Double I_NAI_H3x2z_G3xy_b = I_NAI_I3xy2z_F3x_b+ABY*I_NAI_H3x2z_F3x_b;
  Double I_NAI_H2x3y_G3xy_b = I_NAI_I2x4y_F3x_b+ABY*I_NAI_H2x3y_F3x_b;
  Double I_NAI_H2x2yz_G3xy_b = I_NAI_I2x3yz_F3x_b+ABY*I_NAI_H2x2yz_F3x_b;
  Double I_NAI_H2xy2z_G3xy_b = I_NAI_I2x2y2z_F3x_b+ABY*I_NAI_H2xy2z_F3x_b;
  Double I_NAI_H2x3z_G3xy_b = I_NAI_I2xy3z_F3x_b+ABY*I_NAI_H2x3z_F3x_b;
  Double I_NAI_Hx4y_G3xy_b = I_NAI_Ix5y_F3x_b+ABY*I_NAI_Hx4y_F3x_b;
  Double I_NAI_Hx3yz_G3xy_b = I_NAI_Ix4yz_F3x_b+ABY*I_NAI_Hx3yz_F3x_b;
  Double I_NAI_Hx2y2z_G3xy_b = I_NAI_Ix3y2z_F3x_b+ABY*I_NAI_Hx2y2z_F3x_b;
  Double I_NAI_Hxy3z_G3xy_b = I_NAI_Ix2y3z_F3x_b+ABY*I_NAI_Hxy3z_F3x_b;
  Double I_NAI_Hx4z_G3xy_b = I_NAI_Ixy4z_F3x_b+ABY*I_NAI_Hx4z_F3x_b;
  Double I_NAI_H5y_G3xy_b = I_NAI_I6y_F3x_b+ABY*I_NAI_H5y_F3x_b;
  Double I_NAI_H4yz_G3xy_b = I_NAI_I5yz_F3x_b+ABY*I_NAI_H4yz_F3x_b;
  Double I_NAI_H3y2z_G3xy_b = I_NAI_I4y2z_F3x_b+ABY*I_NAI_H3y2z_F3x_b;
  Double I_NAI_H2y3z_G3xy_b = I_NAI_I3y3z_F3x_b+ABY*I_NAI_H2y3z_F3x_b;
  Double I_NAI_Hy4z_G3xy_b = I_NAI_I2y4z_F3x_b+ABY*I_NAI_Hy4z_F3x_b;
  Double I_NAI_H5z_G3xy_b = I_NAI_Iy5z_F3x_b+ABY*I_NAI_H5z_F3x_b;
  Double I_NAI_H5x_G3xz_b = I_NAI_I5xz_F3x_b+ABZ*I_NAI_H5x_F3x_b;
  Double I_NAI_H4xy_G3xz_b = I_NAI_I4xyz_F3x_b+ABZ*I_NAI_H4xy_F3x_b;
  Double I_NAI_H4xz_G3xz_b = I_NAI_I4x2z_F3x_b+ABZ*I_NAI_H4xz_F3x_b;
  Double I_NAI_H3x2y_G3xz_b = I_NAI_I3x2yz_F3x_b+ABZ*I_NAI_H3x2y_F3x_b;
  Double I_NAI_H3xyz_G3xz_b = I_NAI_I3xy2z_F3x_b+ABZ*I_NAI_H3xyz_F3x_b;
  Double I_NAI_H3x2z_G3xz_b = I_NAI_I3x3z_F3x_b+ABZ*I_NAI_H3x2z_F3x_b;
  Double I_NAI_H2x3y_G3xz_b = I_NAI_I2x3yz_F3x_b+ABZ*I_NAI_H2x3y_F3x_b;
  Double I_NAI_H2x2yz_G3xz_b = I_NAI_I2x2y2z_F3x_b+ABZ*I_NAI_H2x2yz_F3x_b;
  Double I_NAI_H2xy2z_G3xz_b = I_NAI_I2xy3z_F3x_b+ABZ*I_NAI_H2xy2z_F3x_b;
  Double I_NAI_H2x3z_G3xz_b = I_NAI_I2x4z_F3x_b+ABZ*I_NAI_H2x3z_F3x_b;
  Double I_NAI_Hx4y_G3xz_b = I_NAI_Ix4yz_F3x_b+ABZ*I_NAI_Hx4y_F3x_b;
  Double I_NAI_Hx3yz_G3xz_b = I_NAI_Ix3y2z_F3x_b+ABZ*I_NAI_Hx3yz_F3x_b;
  Double I_NAI_Hx2y2z_G3xz_b = I_NAI_Ix2y3z_F3x_b+ABZ*I_NAI_Hx2y2z_F3x_b;
  Double I_NAI_Hxy3z_G3xz_b = I_NAI_Ixy4z_F3x_b+ABZ*I_NAI_Hxy3z_F3x_b;
  Double I_NAI_Hx4z_G3xz_b = I_NAI_Ix5z_F3x_b+ABZ*I_NAI_Hx4z_F3x_b;
  Double I_NAI_H5y_G3xz_b = I_NAI_I5yz_F3x_b+ABZ*I_NAI_H5y_F3x_b;
  Double I_NAI_H4yz_G3xz_b = I_NAI_I4y2z_F3x_b+ABZ*I_NAI_H4yz_F3x_b;
  Double I_NAI_H3y2z_G3xz_b = I_NAI_I3y3z_F3x_b+ABZ*I_NAI_H3y2z_F3x_b;
  Double I_NAI_H2y3z_G3xz_b = I_NAI_I2y4z_F3x_b+ABZ*I_NAI_H2y3z_F3x_b;
  Double I_NAI_Hy4z_G3xz_b = I_NAI_Iy5z_F3x_b+ABZ*I_NAI_Hy4z_F3x_b;
  Double I_NAI_H5z_G3xz_b = I_NAI_I6z_F3x_b+ABZ*I_NAI_H5z_F3x_b;
  Double I_NAI_H5x_G2x2y_b = I_NAI_I5xy_F2xy_b+ABY*I_NAI_H5x_F2xy_b;
  Double I_NAI_H4xy_G2x2y_b = I_NAI_I4x2y_F2xy_b+ABY*I_NAI_H4xy_F2xy_b;
  Double I_NAI_H4xz_G2x2y_b = I_NAI_I4xyz_F2xy_b+ABY*I_NAI_H4xz_F2xy_b;
  Double I_NAI_H3x2y_G2x2y_b = I_NAI_I3x3y_F2xy_b+ABY*I_NAI_H3x2y_F2xy_b;
  Double I_NAI_H3xyz_G2x2y_b = I_NAI_I3x2yz_F2xy_b+ABY*I_NAI_H3xyz_F2xy_b;
  Double I_NAI_H3x2z_G2x2y_b = I_NAI_I3xy2z_F2xy_b+ABY*I_NAI_H3x2z_F2xy_b;
  Double I_NAI_H2x3y_G2x2y_b = I_NAI_I2x4y_F2xy_b+ABY*I_NAI_H2x3y_F2xy_b;
  Double I_NAI_H2x2yz_G2x2y_b = I_NAI_I2x3yz_F2xy_b+ABY*I_NAI_H2x2yz_F2xy_b;
  Double I_NAI_H2xy2z_G2x2y_b = I_NAI_I2x2y2z_F2xy_b+ABY*I_NAI_H2xy2z_F2xy_b;
  Double I_NAI_H2x3z_G2x2y_b = I_NAI_I2xy3z_F2xy_b+ABY*I_NAI_H2x3z_F2xy_b;
  Double I_NAI_Hx4y_G2x2y_b = I_NAI_Ix5y_F2xy_b+ABY*I_NAI_Hx4y_F2xy_b;
  Double I_NAI_Hx3yz_G2x2y_b = I_NAI_Ix4yz_F2xy_b+ABY*I_NAI_Hx3yz_F2xy_b;
  Double I_NAI_Hx2y2z_G2x2y_b = I_NAI_Ix3y2z_F2xy_b+ABY*I_NAI_Hx2y2z_F2xy_b;
  Double I_NAI_Hxy3z_G2x2y_b = I_NAI_Ix2y3z_F2xy_b+ABY*I_NAI_Hxy3z_F2xy_b;
  Double I_NAI_Hx4z_G2x2y_b = I_NAI_Ixy4z_F2xy_b+ABY*I_NAI_Hx4z_F2xy_b;
  Double I_NAI_H5y_G2x2y_b = I_NAI_I6y_F2xy_b+ABY*I_NAI_H5y_F2xy_b;
  Double I_NAI_H4yz_G2x2y_b = I_NAI_I5yz_F2xy_b+ABY*I_NAI_H4yz_F2xy_b;
  Double I_NAI_H3y2z_G2x2y_b = I_NAI_I4y2z_F2xy_b+ABY*I_NAI_H3y2z_F2xy_b;
  Double I_NAI_H2y3z_G2x2y_b = I_NAI_I3y3z_F2xy_b+ABY*I_NAI_H2y3z_F2xy_b;
  Double I_NAI_Hy4z_G2x2y_b = I_NAI_I2y4z_F2xy_b+ABY*I_NAI_Hy4z_F2xy_b;
  Double I_NAI_H5z_G2x2y_b = I_NAI_Iy5z_F2xy_b+ABY*I_NAI_H5z_F2xy_b;
  Double I_NAI_H5x_G2xyz_b = I_NAI_I5xz_F2xy_b+ABZ*I_NAI_H5x_F2xy_b;
  Double I_NAI_H4xy_G2xyz_b = I_NAI_I4xyz_F2xy_b+ABZ*I_NAI_H4xy_F2xy_b;
  Double I_NAI_H4xz_G2xyz_b = I_NAI_I4x2z_F2xy_b+ABZ*I_NAI_H4xz_F2xy_b;
  Double I_NAI_H3x2y_G2xyz_b = I_NAI_I3x2yz_F2xy_b+ABZ*I_NAI_H3x2y_F2xy_b;
  Double I_NAI_H3xyz_G2xyz_b = I_NAI_I3xy2z_F2xy_b+ABZ*I_NAI_H3xyz_F2xy_b;
  Double I_NAI_H3x2z_G2xyz_b = I_NAI_I3x3z_F2xy_b+ABZ*I_NAI_H3x2z_F2xy_b;
  Double I_NAI_H2x3y_G2xyz_b = I_NAI_I2x3yz_F2xy_b+ABZ*I_NAI_H2x3y_F2xy_b;
  Double I_NAI_H2x2yz_G2xyz_b = I_NAI_I2x2y2z_F2xy_b+ABZ*I_NAI_H2x2yz_F2xy_b;
  Double I_NAI_H2xy2z_G2xyz_b = I_NAI_I2xy3z_F2xy_b+ABZ*I_NAI_H2xy2z_F2xy_b;
  Double I_NAI_H2x3z_G2xyz_b = I_NAI_I2x4z_F2xy_b+ABZ*I_NAI_H2x3z_F2xy_b;
  Double I_NAI_Hx4y_G2xyz_b = I_NAI_Ix4yz_F2xy_b+ABZ*I_NAI_Hx4y_F2xy_b;
  Double I_NAI_Hx3yz_G2xyz_b = I_NAI_Ix3y2z_F2xy_b+ABZ*I_NAI_Hx3yz_F2xy_b;
  Double I_NAI_Hx2y2z_G2xyz_b = I_NAI_Ix2y3z_F2xy_b+ABZ*I_NAI_Hx2y2z_F2xy_b;
  Double I_NAI_Hxy3z_G2xyz_b = I_NAI_Ixy4z_F2xy_b+ABZ*I_NAI_Hxy3z_F2xy_b;
  Double I_NAI_Hx4z_G2xyz_b = I_NAI_Ix5z_F2xy_b+ABZ*I_NAI_Hx4z_F2xy_b;
  Double I_NAI_H5y_G2xyz_b = I_NAI_I5yz_F2xy_b+ABZ*I_NAI_H5y_F2xy_b;
  Double I_NAI_H4yz_G2xyz_b = I_NAI_I4y2z_F2xy_b+ABZ*I_NAI_H4yz_F2xy_b;
  Double I_NAI_H3y2z_G2xyz_b = I_NAI_I3y3z_F2xy_b+ABZ*I_NAI_H3y2z_F2xy_b;
  Double I_NAI_H2y3z_G2xyz_b = I_NAI_I2y4z_F2xy_b+ABZ*I_NAI_H2y3z_F2xy_b;
  Double I_NAI_Hy4z_G2xyz_b = I_NAI_Iy5z_F2xy_b+ABZ*I_NAI_Hy4z_F2xy_b;
  Double I_NAI_H5z_G2xyz_b = I_NAI_I6z_F2xy_b+ABZ*I_NAI_H5z_F2xy_b;
  Double I_NAI_H5x_G2x2z_b = I_NAI_I5xz_F2xz_b+ABZ*I_NAI_H5x_F2xz_b;
  Double I_NAI_H4xy_G2x2z_b = I_NAI_I4xyz_F2xz_b+ABZ*I_NAI_H4xy_F2xz_b;
  Double I_NAI_H4xz_G2x2z_b = I_NAI_I4x2z_F2xz_b+ABZ*I_NAI_H4xz_F2xz_b;
  Double I_NAI_H3x2y_G2x2z_b = I_NAI_I3x2yz_F2xz_b+ABZ*I_NAI_H3x2y_F2xz_b;
  Double I_NAI_H3xyz_G2x2z_b = I_NAI_I3xy2z_F2xz_b+ABZ*I_NAI_H3xyz_F2xz_b;
  Double I_NAI_H3x2z_G2x2z_b = I_NAI_I3x3z_F2xz_b+ABZ*I_NAI_H3x2z_F2xz_b;
  Double I_NAI_H2x3y_G2x2z_b = I_NAI_I2x3yz_F2xz_b+ABZ*I_NAI_H2x3y_F2xz_b;
  Double I_NAI_H2x2yz_G2x2z_b = I_NAI_I2x2y2z_F2xz_b+ABZ*I_NAI_H2x2yz_F2xz_b;
  Double I_NAI_H2xy2z_G2x2z_b = I_NAI_I2xy3z_F2xz_b+ABZ*I_NAI_H2xy2z_F2xz_b;
  Double I_NAI_H2x3z_G2x2z_b = I_NAI_I2x4z_F2xz_b+ABZ*I_NAI_H2x3z_F2xz_b;
  Double I_NAI_Hx4y_G2x2z_b = I_NAI_Ix4yz_F2xz_b+ABZ*I_NAI_Hx4y_F2xz_b;
  Double I_NAI_Hx3yz_G2x2z_b = I_NAI_Ix3y2z_F2xz_b+ABZ*I_NAI_Hx3yz_F2xz_b;
  Double I_NAI_Hx2y2z_G2x2z_b = I_NAI_Ix2y3z_F2xz_b+ABZ*I_NAI_Hx2y2z_F2xz_b;
  Double I_NAI_Hxy3z_G2x2z_b = I_NAI_Ixy4z_F2xz_b+ABZ*I_NAI_Hxy3z_F2xz_b;
  Double I_NAI_Hx4z_G2x2z_b = I_NAI_Ix5z_F2xz_b+ABZ*I_NAI_Hx4z_F2xz_b;
  Double I_NAI_H5y_G2x2z_b = I_NAI_I5yz_F2xz_b+ABZ*I_NAI_H5y_F2xz_b;
  Double I_NAI_H4yz_G2x2z_b = I_NAI_I4y2z_F2xz_b+ABZ*I_NAI_H4yz_F2xz_b;
  Double I_NAI_H3y2z_G2x2z_b = I_NAI_I3y3z_F2xz_b+ABZ*I_NAI_H3y2z_F2xz_b;
  Double I_NAI_H2y3z_G2x2z_b = I_NAI_I2y4z_F2xz_b+ABZ*I_NAI_H2y3z_F2xz_b;
  Double I_NAI_Hy4z_G2x2z_b = I_NAI_Iy5z_F2xz_b+ABZ*I_NAI_Hy4z_F2xz_b;
  Double I_NAI_H5z_G2x2z_b = I_NAI_I6z_F2xz_b+ABZ*I_NAI_H5z_F2xz_b;
  Double I_NAI_H5x_Gx3y_b = I_NAI_I6x_F3y_b+ABX*I_NAI_H5x_F3y_b;
  Double I_NAI_H4xy_Gx3y_b = I_NAI_I5xy_F3y_b+ABX*I_NAI_H4xy_F3y_b;
  Double I_NAI_H4xz_Gx3y_b = I_NAI_I5xz_F3y_b+ABX*I_NAI_H4xz_F3y_b;
  Double I_NAI_H3x2y_Gx3y_b = I_NAI_I4x2y_F3y_b+ABX*I_NAI_H3x2y_F3y_b;
  Double I_NAI_H3xyz_Gx3y_b = I_NAI_I4xyz_F3y_b+ABX*I_NAI_H3xyz_F3y_b;
  Double I_NAI_H3x2z_Gx3y_b = I_NAI_I4x2z_F3y_b+ABX*I_NAI_H3x2z_F3y_b;
  Double I_NAI_H2x3y_Gx3y_b = I_NAI_I3x3y_F3y_b+ABX*I_NAI_H2x3y_F3y_b;
  Double I_NAI_H2x2yz_Gx3y_b = I_NAI_I3x2yz_F3y_b+ABX*I_NAI_H2x2yz_F3y_b;
  Double I_NAI_H2xy2z_Gx3y_b = I_NAI_I3xy2z_F3y_b+ABX*I_NAI_H2xy2z_F3y_b;
  Double I_NAI_H2x3z_Gx3y_b = I_NAI_I3x3z_F3y_b+ABX*I_NAI_H2x3z_F3y_b;
  Double I_NAI_Hx4y_Gx3y_b = I_NAI_I2x4y_F3y_b+ABX*I_NAI_Hx4y_F3y_b;
  Double I_NAI_Hx3yz_Gx3y_b = I_NAI_I2x3yz_F3y_b+ABX*I_NAI_Hx3yz_F3y_b;
  Double I_NAI_Hx2y2z_Gx3y_b = I_NAI_I2x2y2z_F3y_b+ABX*I_NAI_Hx2y2z_F3y_b;
  Double I_NAI_Hxy3z_Gx3y_b = I_NAI_I2xy3z_F3y_b+ABX*I_NAI_Hxy3z_F3y_b;
  Double I_NAI_Hx4z_Gx3y_b = I_NAI_I2x4z_F3y_b+ABX*I_NAI_Hx4z_F3y_b;
  Double I_NAI_H5y_Gx3y_b = I_NAI_Ix5y_F3y_b+ABX*I_NAI_H5y_F3y_b;
  Double I_NAI_H4yz_Gx3y_b = I_NAI_Ix4yz_F3y_b+ABX*I_NAI_H4yz_F3y_b;
  Double I_NAI_H3y2z_Gx3y_b = I_NAI_Ix3y2z_F3y_b+ABX*I_NAI_H3y2z_F3y_b;
  Double I_NAI_H2y3z_Gx3y_b = I_NAI_Ix2y3z_F3y_b+ABX*I_NAI_H2y3z_F3y_b;
  Double I_NAI_Hy4z_Gx3y_b = I_NAI_Ixy4z_F3y_b+ABX*I_NAI_Hy4z_F3y_b;
  Double I_NAI_H5z_Gx3y_b = I_NAI_Ix5z_F3y_b+ABX*I_NAI_H5z_F3y_b;
  Double I_NAI_H5x_Gx2yz_b = I_NAI_I5xz_Fx2y_b+ABZ*I_NAI_H5x_Fx2y_b;
  Double I_NAI_H4xy_Gx2yz_b = I_NAI_I4xyz_Fx2y_b+ABZ*I_NAI_H4xy_Fx2y_b;
  Double I_NAI_H4xz_Gx2yz_b = I_NAI_I4x2z_Fx2y_b+ABZ*I_NAI_H4xz_Fx2y_b;
  Double I_NAI_H3x2y_Gx2yz_b = I_NAI_I3x2yz_Fx2y_b+ABZ*I_NAI_H3x2y_Fx2y_b;
  Double I_NAI_H3xyz_Gx2yz_b = I_NAI_I3xy2z_Fx2y_b+ABZ*I_NAI_H3xyz_Fx2y_b;
  Double I_NAI_H3x2z_Gx2yz_b = I_NAI_I3x3z_Fx2y_b+ABZ*I_NAI_H3x2z_Fx2y_b;
  Double I_NAI_H2x3y_Gx2yz_b = I_NAI_I2x3yz_Fx2y_b+ABZ*I_NAI_H2x3y_Fx2y_b;
  Double I_NAI_H2x2yz_Gx2yz_b = I_NAI_I2x2y2z_Fx2y_b+ABZ*I_NAI_H2x2yz_Fx2y_b;
  Double I_NAI_H2xy2z_Gx2yz_b = I_NAI_I2xy3z_Fx2y_b+ABZ*I_NAI_H2xy2z_Fx2y_b;
  Double I_NAI_H2x3z_Gx2yz_b = I_NAI_I2x4z_Fx2y_b+ABZ*I_NAI_H2x3z_Fx2y_b;
  Double I_NAI_Hx4y_Gx2yz_b = I_NAI_Ix4yz_Fx2y_b+ABZ*I_NAI_Hx4y_Fx2y_b;
  Double I_NAI_Hx3yz_Gx2yz_b = I_NAI_Ix3y2z_Fx2y_b+ABZ*I_NAI_Hx3yz_Fx2y_b;
  Double I_NAI_Hx2y2z_Gx2yz_b = I_NAI_Ix2y3z_Fx2y_b+ABZ*I_NAI_Hx2y2z_Fx2y_b;
  Double I_NAI_Hxy3z_Gx2yz_b = I_NAI_Ixy4z_Fx2y_b+ABZ*I_NAI_Hxy3z_Fx2y_b;
  Double I_NAI_Hx4z_Gx2yz_b = I_NAI_Ix5z_Fx2y_b+ABZ*I_NAI_Hx4z_Fx2y_b;
  Double I_NAI_H5y_Gx2yz_b = I_NAI_I5yz_Fx2y_b+ABZ*I_NAI_H5y_Fx2y_b;
  Double I_NAI_H4yz_Gx2yz_b = I_NAI_I4y2z_Fx2y_b+ABZ*I_NAI_H4yz_Fx2y_b;
  Double I_NAI_H3y2z_Gx2yz_b = I_NAI_I3y3z_Fx2y_b+ABZ*I_NAI_H3y2z_Fx2y_b;
  Double I_NAI_H2y3z_Gx2yz_b = I_NAI_I2y4z_Fx2y_b+ABZ*I_NAI_H2y3z_Fx2y_b;
  Double I_NAI_Hy4z_Gx2yz_b = I_NAI_Iy5z_Fx2y_b+ABZ*I_NAI_Hy4z_Fx2y_b;
  Double I_NAI_H5z_Gx2yz_b = I_NAI_I6z_Fx2y_b+ABZ*I_NAI_H5z_Fx2y_b;
  Double I_NAI_H5x_Gxy2z_b = I_NAI_I5xy_Fx2z_b+ABY*I_NAI_H5x_Fx2z_b;
  Double I_NAI_H4xy_Gxy2z_b = I_NAI_I4x2y_Fx2z_b+ABY*I_NAI_H4xy_Fx2z_b;
  Double I_NAI_H4xz_Gxy2z_b = I_NAI_I4xyz_Fx2z_b+ABY*I_NAI_H4xz_Fx2z_b;
  Double I_NAI_H3x2y_Gxy2z_b = I_NAI_I3x3y_Fx2z_b+ABY*I_NAI_H3x2y_Fx2z_b;
  Double I_NAI_H3xyz_Gxy2z_b = I_NAI_I3x2yz_Fx2z_b+ABY*I_NAI_H3xyz_Fx2z_b;
  Double I_NAI_H3x2z_Gxy2z_b = I_NAI_I3xy2z_Fx2z_b+ABY*I_NAI_H3x2z_Fx2z_b;
  Double I_NAI_H2x3y_Gxy2z_b = I_NAI_I2x4y_Fx2z_b+ABY*I_NAI_H2x3y_Fx2z_b;
  Double I_NAI_H2x2yz_Gxy2z_b = I_NAI_I2x3yz_Fx2z_b+ABY*I_NAI_H2x2yz_Fx2z_b;
  Double I_NAI_H2xy2z_Gxy2z_b = I_NAI_I2x2y2z_Fx2z_b+ABY*I_NAI_H2xy2z_Fx2z_b;
  Double I_NAI_H2x3z_Gxy2z_b = I_NAI_I2xy3z_Fx2z_b+ABY*I_NAI_H2x3z_Fx2z_b;
  Double I_NAI_Hx4y_Gxy2z_b = I_NAI_Ix5y_Fx2z_b+ABY*I_NAI_Hx4y_Fx2z_b;
  Double I_NAI_Hx3yz_Gxy2z_b = I_NAI_Ix4yz_Fx2z_b+ABY*I_NAI_Hx3yz_Fx2z_b;
  Double I_NAI_Hx2y2z_Gxy2z_b = I_NAI_Ix3y2z_Fx2z_b+ABY*I_NAI_Hx2y2z_Fx2z_b;
  Double I_NAI_Hxy3z_Gxy2z_b = I_NAI_Ix2y3z_Fx2z_b+ABY*I_NAI_Hxy3z_Fx2z_b;
  Double I_NAI_Hx4z_Gxy2z_b = I_NAI_Ixy4z_Fx2z_b+ABY*I_NAI_Hx4z_Fx2z_b;
  Double I_NAI_H5y_Gxy2z_b = I_NAI_I6y_Fx2z_b+ABY*I_NAI_H5y_Fx2z_b;
  Double I_NAI_H4yz_Gxy2z_b = I_NAI_I5yz_Fx2z_b+ABY*I_NAI_H4yz_Fx2z_b;
  Double I_NAI_H3y2z_Gxy2z_b = I_NAI_I4y2z_Fx2z_b+ABY*I_NAI_H3y2z_Fx2z_b;
  Double I_NAI_H2y3z_Gxy2z_b = I_NAI_I3y3z_Fx2z_b+ABY*I_NAI_H2y3z_Fx2z_b;
  Double I_NAI_Hy4z_Gxy2z_b = I_NAI_I2y4z_Fx2z_b+ABY*I_NAI_Hy4z_Fx2z_b;
  Double I_NAI_H5z_Gxy2z_b = I_NAI_Iy5z_Fx2z_b+ABY*I_NAI_H5z_Fx2z_b;
  Double I_NAI_H5x_Gx3z_b = I_NAI_I6x_F3z_b+ABX*I_NAI_H5x_F3z_b;
  Double I_NAI_H4xy_Gx3z_b = I_NAI_I5xy_F3z_b+ABX*I_NAI_H4xy_F3z_b;
  Double I_NAI_H4xz_Gx3z_b = I_NAI_I5xz_F3z_b+ABX*I_NAI_H4xz_F3z_b;
  Double I_NAI_H3x2y_Gx3z_b = I_NAI_I4x2y_F3z_b+ABX*I_NAI_H3x2y_F3z_b;
  Double I_NAI_H3xyz_Gx3z_b = I_NAI_I4xyz_F3z_b+ABX*I_NAI_H3xyz_F3z_b;
  Double I_NAI_H3x2z_Gx3z_b = I_NAI_I4x2z_F3z_b+ABX*I_NAI_H3x2z_F3z_b;
  Double I_NAI_H2x3y_Gx3z_b = I_NAI_I3x3y_F3z_b+ABX*I_NAI_H2x3y_F3z_b;
  Double I_NAI_H2x2yz_Gx3z_b = I_NAI_I3x2yz_F3z_b+ABX*I_NAI_H2x2yz_F3z_b;
  Double I_NAI_H2xy2z_Gx3z_b = I_NAI_I3xy2z_F3z_b+ABX*I_NAI_H2xy2z_F3z_b;
  Double I_NAI_H2x3z_Gx3z_b = I_NAI_I3x3z_F3z_b+ABX*I_NAI_H2x3z_F3z_b;
  Double I_NAI_Hx4y_Gx3z_b = I_NAI_I2x4y_F3z_b+ABX*I_NAI_Hx4y_F3z_b;
  Double I_NAI_Hx3yz_Gx3z_b = I_NAI_I2x3yz_F3z_b+ABX*I_NAI_Hx3yz_F3z_b;
  Double I_NAI_Hx2y2z_Gx3z_b = I_NAI_I2x2y2z_F3z_b+ABX*I_NAI_Hx2y2z_F3z_b;
  Double I_NAI_Hxy3z_Gx3z_b = I_NAI_I2xy3z_F3z_b+ABX*I_NAI_Hxy3z_F3z_b;
  Double I_NAI_Hx4z_Gx3z_b = I_NAI_I2x4z_F3z_b+ABX*I_NAI_Hx4z_F3z_b;
  Double I_NAI_H5y_Gx3z_b = I_NAI_Ix5y_F3z_b+ABX*I_NAI_H5y_F3z_b;
  Double I_NAI_H4yz_Gx3z_b = I_NAI_Ix4yz_F3z_b+ABX*I_NAI_H4yz_F3z_b;
  Double I_NAI_H3y2z_Gx3z_b = I_NAI_Ix3y2z_F3z_b+ABX*I_NAI_H3y2z_F3z_b;
  Double I_NAI_H2y3z_Gx3z_b = I_NAI_Ix2y3z_F3z_b+ABX*I_NAI_H2y3z_F3z_b;
  Double I_NAI_Hy4z_Gx3z_b = I_NAI_Ixy4z_F3z_b+ABX*I_NAI_Hy4z_F3z_b;
  Double I_NAI_H5z_Gx3z_b = I_NAI_Ix5z_F3z_b+ABX*I_NAI_H5z_F3z_b;
  Double I_NAI_H5x_G4y_b = I_NAI_I5xy_F3y_b+ABY*I_NAI_H5x_F3y_b;
  Double I_NAI_H4xy_G4y_b = I_NAI_I4x2y_F3y_b+ABY*I_NAI_H4xy_F3y_b;
  Double I_NAI_H4xz_G4y_b = I_NAI_I4xyz_F3y_b+ABY*I_NAI_H4xz_F3y_b;
  Double I_NAI_H3x2y_G4y_b = I_NAI_I3x3y_F3y_b+ABY*I_NAI_H3x2y_F3y_b;
  Double I_NAI_H3xyz_G4y_b = I_NAI_I3x2yz_F3y_b+ABY*I_NAI_H3xyz_F3y_b;
  Double I_NAI_H3x2z_G4y_b = I_NAI_I3xy2z_F3y_b+ABY*I_NAI_H3x2z_F3y_b;
  Double I_NAI_H2x3y_G4y_b = I_NAI_I2x4y_F3y_b+ABY*I_NAI_H2x3y_F3y_b;
  Double I_NAI_H2x2yz_G4y_b = I_NAI_I2x3yz_F3y_b+ABY*I_NAI_H2x2yz_F3y_b;
  Double I_NAI_H2xy2z_G4y_b = I_NAI_I2x2y2z_F3y_b+ABY*I_NAI_H2xy2z_F3y_b;
  Double I_NAI_H2x3z_G4y_b = I_NAI_I2xy3z_F3y_b+ABY*I_NAI_H2x3z_F3y_b;
  Double I_NAI_Hx4y_G4y_b = I_NAI_Ix5y_F3y_b+ABY*I_NAI_Hx4y_F3y_b;
  Double I_NAI_Hx3yz_G4y_b = I_NAI_Ix4yz_F3y_b+ABY*I_NAI_Hx3yz_F3y_b;
  Double I_NAI_Hx2y2z_G4y_b = I_NAI_Ix3y2z_F3y_b+ABY*I_NAI_Hx2y2z_F3y_b;
  Double I_NAI_Hxy3z_G4y_b = I_NAI_Ix2y3z_F3y_b+ABY*I_NAI_Hxy3z_F3y_b;
  Double I_NAI_Hx4z_G4y_b = I_NAI_Ixy4z_F3y_b+ABY*I_NAI_Hx4z_F3y_b;
  Double I_NAI_H5y_G4y_b = I_NAI_I6y_F3y_b+ABY*I_NAI_H5y_F3y_b;
  Double I_NAI_H4yz_G4y_b = I_NAI_I5yz_F3y_b+ABY*I_NAI_H4yz_F3y_b;
  Double I_NAI_H3y2z_G4y_b = I_NAI_I4y2z_F3y_b+ABY*I_NAI_H3y2z_F3y_b;
  Double I_NAI_H2y3z_G4y_b = I_NAI_I3y3z_F3y_b+ABY*I_NAI_H2y3z_F3y_b;
  Double I_NAI_Hy4z_G4y_b = I_NAI_I2y4z_F3y_b+ABY*I_NAI_Hy4z_F3y_b;
  Double I_NAI_H5z_G4y_b = I_NAI_Iy5z_F3y_b+ABY*I_NAI_H5z_F3y_b;
  Double I_NAI_H5x_G3yz_b = I_NAI_I5xz_F3y_b+ABZ*I_NAI_H5x_F3y_b;
  Double I_NAI_H4xy_G3yz_b = I_NAI_I4xyz_F3y_b+ABZ*I_NAI_H4xy_F3y_b;
  Double I_NAI_H4xz_G3yz_b = I_NAI_I4x2z_F3y_b+ABZ*I_NAI_H4xz_F3y_b;
  Double I_NAI_H3x2y_G3yz_b = I_NAI_I3x2yz_F3y_b+ABZ*I_NAI_H3x2y_F3y_b;
  Double I_NAI_H3xyz_G3yz_b = I_NAI_I3xy2z_F3y_b+ABZ*I_NAI_H3xyz_F3y_b;
  Double I_NAI_H3x2z_G3yz_b = I_NAI_I3x3z_F3y_b+ABZ*I_NAI_H3x2z_F3y_b;
  Double I_NAI_H2x3y_G3yz_b = I_NAI_I2x3yz_F3y_b+ABZ*I_NAI_H2x3y_F3y_b;
  Double I_NAI_H2x2yz_G3yz_b = I_NAI_I2x2y2z_F3y_b+ABZ*I_NAI_H2x2yz_F3y_b;
  Double I_NAI_H2xy2z_G3yz_b = I_NAI_I2xy3z_F3y_b+ABZ*I_NAI_H2xy2z_F3y_b;
  Double I_NAI_H2x3z_G3yz_b = I_NAI_I2x4z_F3y_b+ABZ*I_NAI_H2x3z_F3y_b;
  Double I_NAI_Hx4y_G3yz_b = I_NAI_Ix4yz_F3y_b+ABZ*I_NAI_Hx4y_F3y_b;
  Double I_NAI_Hx3yz_G3yz_b = I_NAI_Ix3y2z_F3y_b+ABZ*I_NAI_Hx3yz_F3y_b;
  Double I_NAI_Hx2y2z_G3yz_b = I_NAI_Ix2y3z_F3y_b+ABZ*I_NAI_Hx2y2z_F3y_b;
  Double I_NAI_Hxy3z_G3yz_b = I_NAI_Ixy4z_F3y_b+ABZ*I_NAI_Hxy3z_F3y_b;
  Double I_NAI_Hx4z_G3yz_b = I_NAI_Ix5z_F3y_b+ABZ*I_NAI_Hx4z_F3y_b;
  Double I_NAI_H5y_G3yz_b = I_NAI_I5yz_F3y_b+ABZ*I_NAI_H5y_F3y_b;
  Double I_NAI_H4yz_G3yz_b = I_NAI_I4y2z_F3y_b+ABZ*I_NAI_H4yz_F3y_b;
  Double I_NAI_H3y2z_G3yz_b = I_NAI_I3y3z_F3y_b+ABZ*I_NAI_H3y2z_F3y_b;
  Double I_NAI_H2y3z_G3yz_b = I_NAI_I2y4z_F3y_b+ABZ*I_NAI_H2y3z_F3y_b;
  Double I_NAI_Hy4z_G3yz_b = I_NAI_Iy5z_F3y_b+ABZ*I_NAI_Hy4z_F3y_b;
  Double I_NAI_H5z_G3yz_b = I_NAI_I6z_F3y_b+ABZ*I_NAI_H5z_F3y_b;
  Double I_NAI_H5x_G2y2z_b = I_NAI_I5xz_F2yz_b+ABZ*I_NAI_H5x_F2yz_b;
  Double I_NAI_H4xy_G2y2z_b = I_NAI_I4xyz_F2yz_b+ABZ*I_NAI_H4xy_F2yz_b;
  Double I_NAI_H4xz_G2y2z_b = I_NAI_I4x2z_F2yz_b+ABZ*I_NAI_H4xz_F2yz_b;
  Double I_NAI_H3x2y_G2y2z_b = I_NAI_I3x2yz_F2yz_b+ABZ*I_NAI_H3x2y_F2yz_b;
  Double I_NAI_H3xyz_G2y2z_b = I_NAI_I3xy2z_F2yz_b+ABZ*I_NAI_H3xyz_F2yz_b;
  Double I_NAI_H3x2z_G2y2z_b = I_NAI_I3x3z_F2yz_b+ABZ*I_NAI_H3x2z_F2yz_b;
  Double I_NAI_H2x3y_G2y2z_b = I_NAI_I2x3yz_F2yz_b+ABZ*I_NAI_H2x3y_F2yz_b;
  Double I_NAI_H2x2yz_G2y2z_b = I_NAI_I2x2y2z_F2yz_b+ABZ*I_NAI_H2x2yz_F2yz_b;
  Double I_NAI_H2xy2z_G2y2z_b = I_NAI_I2xy3z_F2yz_b+ABZ*I_NAI_H2xy2z_F2yz_b;
  Double I_NAI_H2x3z_G2y2z_b = I_NAI_I2x4z_F2yz_b+ABZ*I_NAI_H2x3z_F2yz_b;
  Double I_NAI_Hx4y_G2y2z_b = I_NAI_Ix4yz_F2yz_b+ABZ*I_NAI_Hx4y_F2yz_b;
  Double I_NAI_Hx3yz_G2y2z_b = I_NAI_Ix3y2z_F2yz_b+ABZ*I_NAI_Hx3yz_F2yz_b;
  Double I_NAI_Hx2y2z_G2y2z_b = I_NAI_Ix2y3z_F2yz_b+ABZ*I_NAI_Hx2y2z_F2yz_b;
  Double I_NAI_Hxy3z_G2y2z_b = I_NAI_Ixy4z_F2yz_b+ABZ*I_NAI_Hxy3z_F2yz_b;
  Double I_NAI_Hx4z_G2y2z_b = I_NAI_Ix5z_F2yz_b+ABZ*I_NAI_Hx4z_F2yz_b;
  Double I_NAI_H5y_G2y2z_b = I_NAI_I5yz_F2yz_b+ABZ*I_NAI_H5y_F2yz_b;
  Double I_NAI_H4yz_G2y2z_b = I_NAI_I4y2z_F2yz_b+ABZ*I_NAI_H4yz_F2yz_b;
  Double I_NAI_H3y2z_G2y2z_b = I_NAI_I3y3z_F2yz_b+ABZ*I_NAI_H3y2z_F2yz_b;
  Double I_NAI_H2y3z_G2y2z_b = I_NAI_I2y4z_F2yz_b+ABZ*I_NAI_H2y3z_F2yz_b;
  Double I_NAI_Hy4z_G2y2z_b = I_NAI_Iy5z_F2yz_b+ABZ*I_NAI_Hy4z_F2yz_b;
  Double I_NAI_H5z_G2y2z_b = I_NAI_I6z_F2yz_b+ABZ*I_NAI_H5z_F2yz_b;
  Double I_NAI_H5x_Gy3z_b = I_NAI_I5xy_F3z_b+ABY*I_NAI_H5x_F3z_b;
  Double I_NAI_H4xy_Gy3z_b = I_NAI_I4x2y_F3z_b+ABY*I_NAI_H4xy_F3z_b;
  Double I_NAI_H4xz_Gy3z_b = I_NAI_I4xyz_F3z_b+ABY*I_NAI_H4xz_F3z_b;
  Double I_NAI_H3x2y_Gy3z_b = I_NAI_I3x3y_F3z_b+ABY*I_NAI_H3x2y_F3z_b;
  Double I_NAI_H3xyz_Gy3z_b = I_NAI_I3x2yz_F3z_b+ABY*I_NAI_H3xyz_F3z_b;
  Double I_NAI_H3x2z_Gy3z_b = I_NAI_I3xy2z_F3z_b+ABY*I_NAI_H3x2z_F3z_b;
  Double I_NAI_H2x3y_Gy3z_b = I_NAI_I2x4y_F3z_b+ABY*I_NAI_H2x3y_F3z_b;
  Double I_NAI_H2x2yz_Gy3z_b = I_NAI_I2x3yz_F3z_b+ABY*I_NAI_H2x2yz_F3z_b;
  Double I_NAI_H2xy2z_Gy3z_b = I_NAI_I2x2y2z_F3z_b+ABY*I_NAI_H2xy2z_F3z_b;
  Double I_NAI_H2x3z_Gy3z_b = I_NAI_I2xy3z_F3z_b+ABY*I_NAI_H2x3z_F3z_b;
  Double I_NAI_Hx4y_Gy3z_b = I_NAI_Ix5y_F3z_b+ABY*I_NAI_Hx4y_F3z_b;
  Double I_NAI_Hx3yz_Gy3z_b = I_NAI_Ix4yz_F3z_b+ABY*I_NAI_Hx3yz_F3z_b;
  Double I_NAI_Hx2y2z_Gy3z_b = I_NAI_Ix3y2z_F3z_b+ABY*I_NAI_Hx2y2z_F3z_b;
  Double I_NAI_Hxy3z_Gy3z_b = I_NAI_Ix2y3z_F3z_b+ABY*I_NAI_Hxy3z_F3z_b;
  Double I_NAI_Hx4z_Gy3z_b = I_NAI_Ixy4z_F3z_b+ABY*I_NAI_Hx4z_F3z_b;
  Double I_NAI_H5y_Gy3z_b = I_NAI_I6y_F3z_b+ABY*I_NAI_H5y_F3z_b;
  Double I_NAI_H4yz_Gy3z_b = I_NAI_I5yz_F3z_b+ABY*I_NAI_H4yz_F3z_b;
  Double I_NAI_H3y2z_Gy3z_b = I_NAI_I4y2z_F3z_b+ABY*I_NAI_H3y2z_F3z_b;
  Double I_NAI_H2y3z_Gy3z_b = I_NAI_I3y3z_F3z_b+ABY*I_NAI_H2y3z_F3z_b;
  Double I_NAI_Hy4z_Gy3z_b = I_NAI_I2y4z_F3z_b+ABY*I_NAI_Hy4z_F3z_b;
  Double I_NAI_H5z_Gy3z_b = I_NAI_Iy5z_F3z_b+ABY*I_NAI_H5z_F3z_b;
  Double I_NAI_H5x_G4z_b = I_NAI_I5xz_F3z_b+ABZ*I_NAI_H5x_F3z_b;
  Double I_NAI_H4xy_G4z_b = I_NAI_I4xyz_F3z_b+ABZ*I_NAI_H4xy_F3z_b;
  Double I_NAI_H4xz_G4z_b = I_NAI_I4x2z_F3z_b+ABZ*I_NAI_H4xz_F3z_b;
  Double I_NAI_H3x2y_G4z_b = I_NAI_I3x2yz_F3z_b+ABZ*I_NAI_H3x2y_F3z_b;
  Double I_NAI_H3xyz_G4z_b = I_NAI_I3xy2z_F3z_b+ABZ*I_NAI_H3xyz_F3z_b;
  Double I_NAI_H3x2z_G4z_b = I_NAI_I3x3z_F3z_b+ABZ*I_NAI_H3x2z_F3z_b;
  Double I_NAI_H2x3y_G4z_b = I_NAI_I2x3yz_F3z_b+ABZ*I_NAI_H2x3y_F3z_b;
  Double I_NAI_H2x2yz_G4z_b = I_NAI_I2x2y2z_F3z_b+ABZ*I_NAI_H2x2yz_F3z_b;
  Double I_NAI_H2xy2z_G4z_b = I_NAI_I2xy3z_F3z_b+ABZ*I_NAI_H2xy2z_F3z_b;
  Double I_NAI_H2x3z_G4z_b = I_NAI_I2x4z_F3z_b+ABZ*I_NAI_H2x3z_F3z_b;
  Double I_NAI_Hx4y_G4z_b = I_NAI_Ix4yz_F3z_b+ABZ*I_NAI_Hx4y_F3z_b;
  Double I_NAI_Hx3yz_G4z_b = I_NAI_Ix3y2z_F3z_b+ABZ*I_NAI_Hx3yz_F3z_b;
  Double I_NAI_Hx2y2z_G4z_b = I_NAI_Ix2y3z_F3z_b+ABZ*I_NAI_Hx2y2z_F3z_b;
  Double I_NAI_Hxy3z_G4z_b = I_NAI_Ixy4z_F3z_b+ABZ*I_NAI_Hxy3z_F3z_b;
  Double I_NAI_Hx4z_G4z_b = I_NAI_Ix5z_F3z_b+ABZ*I_NAI_Hx4z_F3z_b;
  Double I_NAI_H5y_G4z_b = I_NAI_I5yz_F3z_b+ABZ*I_NAI_H5y_F3z_b;
  Double I_NAI_H4yz_G4z_b = I_NAI_I4y2z_F3z_b+ABZ*I_NAI_H4yz_F3z_b;
  Double I_NAI_H3y2z_G4z_b = I_NAI_I3y3z_F3z_b+ABZ*I_NAI_H3y2z_F3z_b;
  Double I_NAI_H2y3z_G4z_b = I_NAI_I2y4z_F3z_b+ABZ*I_NAI_H2y3z_F3z_b;
  Double I_NAI_Hy4z_G4z_b = I_NAI_Iy5z_F3z_b+ABZ*I_NAI_Hy4z_F3z_b;
  Double I_NAI_H5z_G4z_b = I_NAI_I6z_F3z_b+ABZ*I_NAI_H5z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_a
   * RHS shell quartet name: SQ_NAI_G_F
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_I6x_F3x_a-5*I_NAI_G4x_F3x;
  abcd[1] = 2.0E0*I_NAI_I5xy_F3x_a-4*I_NAI_G3xy_F3x;
  abcd[2] = 2.0E0*I_NAI_I5xz_F3x_a-4*I_NAI_G3xz_F3x;
  abcd[3] = 2.0E0*I_NAI_I4x2y_F3x_a-3*I_NAI_G2x2y_F3x;
  abcd[4] = 2.0E0*I_NAI_I4xyz_F3x_a-3*I_NAI_G2xyz_F3x;
  abcd[5] = 2.0E0*I_NAI_I4x2z_F3x_a-3*I_NAI_G2x2z_F3x;
  abcd[6] = 2.0E0*I_NAI_I3x3y_F3x_a-2*I_NAI_Gx3y_F3x;
  abcd[7] = 2.0E0*I_NAI_I3x2yz_F3x_a-2*I_NAI_Gx2yz_F3x;
  abcd[8] = 2.0E0*I_NAI_I3xy2z_F3x_a-2*I_NAI_Gxy2z_F3x;
  abcd[9] = 2.0E0*I_NAI_I3x3z_F3x_a-2*I_NAI_Gx3z_F3x;
  abcd[10] = 2.0E0*I_NAI_I2x4y_F3x_a-1*I_NAI_G4y_F3x;
  abcd[11] = 2.0E0*I_NAI_I2x3yz_F3x_a-1*I_NAI_G3yz_F3x;
  abcd[12] = 2.0E0*I_NAI_I2x2y2z_F3x_a-1*I_NAI_G2y2z_F3x;
  abcd[13] = 2.0E0*I_NAI_I2xy3z_F3x_a-1*I_NAI_Gy3z_F3x;
  abcd[14] = 2.0E0*I_NAI_I2x4z_F3x_a-1*I_NAI_G4z_F3x;
  abcd[15] = 2.0E0*I_NAI_Ix5y_F3x_a;
  abcd[16] = 2.0E0*I_NAI_Ix4yz_F3x_a;
  abcd[17] = 2.0E0*I_NAI_Ix3y2z_F3x_a;
  abcd[18] = 2.0E0*I_NAI_Ix2y3z_F3x_a;
  abcd[19] = 2.0E0*I_NAI_Ixy4z_F3x_a;
  abcd[20] = 2.0E0*I_NAI_Ix5z_F3x_a;
  abcd[21] = 2.0E0*I_NAI_I6x_F2xy_a-5*I_NAI_G4x_F2xy;
  abcd[22] = 2.0E0*I_NAI_I5xy_F2xy_a-4*I_NAI_G3xy_F2xy;
  abcd[23] = 2.0E0*I_NAI_I5xz_F2xy_a-4*I_NAI_G3xz_F2xy;
  abcd[24] = 2.0E0*I_NAI_I4x2y_F2xy_a-3*I_NAI_G2x2y_F2xy;
  abcd[25] = 2.0E0*I_NAI_I4xyz_F2xy_a-3*I_NAI_G2xyz_F2xy;
  abcd[26] = 2.0E0*I_NAI_I4x2z_F2xy_a-3*I_NAI_G2x2z_F2xy;
  abcd[27] = 2.0E0*I_NAI_I3x3y_F2xy_a-2*I_NAI_Gx3y_F2xy;
  abcd[28] = 2.0E0*I_NAI_I3x2yz_F2xy_a-2*I_NAI_Gx2yz_F2xy;
  abcd[29] = 2.0E0*I_NAI_I3xy2z_F2xy_a-2*I_NAI_Gxy2z_F2xy;
  abcd[30] = 2.0E0*I_NAI_I3x3z_F2xy_a-2*I_NAI_Gx3z_F2xy;
  abcd[31] = 2.0E0*I_NAI_I2x4y_F2xy_a-1*I_NAI_G4y_F2xy;
  abcd[32] = 2.0E0*I_NAI_I2x3yz_F2xy_a-1*I_NAI_G3yz_F2xy;
  abcd[33] = 2.0E0*I_NAI_I2x2y2z_F2xy_a-1*I_NAI_G2y2z_F2xy;
  abcd[34] = 2.0E0*I_NAI_I2xy3z_F2xy_a-1*I_NAI_Gy3z_F2xy;
  abcd[35] = 2.0E0*I_NAI_I2x4z_F2xy_a-1*I_NAI_G4z_F2xy;
  abcd[36] = 2.0E0*I_NAI_Ix5y_F2xy_a;
  abcd[37] = 2.0E0*I_NAI_Ix4yz_F2xy_a;
  abcd[38] = 2.0E0*I_NAI_Ix3y2z_F2xy_a;
  abcd[39] = 2.0E0*I_NAI_Ix2y3z_F2xy_a;
  abcd[40] = 2.0E0*I_NAI_Ixy4z_F2xy_a;
  abcd[41] = 2.0E0*I_NAI_Ix5z_F2xy_a;
  abcd[42] = 2.0E0*I_NAI_I6x_F2xz_a-5*I_NAI_G4x_F2xz;
  abcd[43] = 2.0E0*I_NAI_I5xy_F2xz_a-4*I_NAI_G3xy_F2xz;
  abcd[44] = 2.0E0*I_NAI_I5xz_F2xz_a-4*I_NAI_G3xz_F2xz;
  abcd[45] = 2.0E0*I_NAI_I4x2y_F2xz_a-3*I_NAI_G2x2y_F2xz;
  abcd[46] = 2.0E0*I_NAI_I4xyz_F2xz_a-3*I_NAI_G2xyz_F2xz;
  abcd[47] = 2.0E0*I_NAI_I4x2z_F2xz_a-3*I_NAI_G2x2z_F2xz;
  abcd[48] = 2.0E0*I_NAI_I3x3y_F2xz_a-2*I_NAI_Gx3y_F2xz;
  abcd[49] = 2.0E0*I_NAI_I3x2yz_F2xz_a-2*I_NAI_Gx2yz_F2xz;
  abcd[50] = 2.0E0*I_NAI_I3xy2z_F2xz_a-2*I_NAI_Gxy2z_F2xz;
  abcd[51] = 2.0E0*I_NAI_I3x3z_F2xz_a-2*I_NAI_Gx3z_F2xz;
  abcd[52] = 2.0E0*I_NAI_I2x4y_F2xz_a-1*I_NAI_G4y_F2xz;
  abcd[53] = 2.0E0*I_NAI_I2x3yz_F2xz_a-1*I_NAI_G3yz_F2xz;
  abcd[54] = 2.0E0*I_NAI_I2x2y2z_F2xz_a-1*I_NAI_G2y2z_F2xz;
  abcd[55] = 2.0E0*I_NAI_I2xy3z_F2xz_a-1*I_NAI_Gy3z_F2xz;
  abcd[56] = 2.0E0*I_NAI_I2x4z_F2xz_a-1*I_NAI_G4z_F2xz;
  abcd[57] = 2.0E0*I_NAI_Ix5y_F2xz_a;
  abcd[58] = 2.0E0*I_NAI_Ix4yz_F2xz_a;
  abcd[59] = 2.0E0*I_NAI_Ix3y2z_F2xz_a;
  abcd[60] = 2.0E0*I_NAI_Ix2y3z_F2xz_a;
  abcd[61] = 2.0E0*I_NAI_Ixy4z_F2xz_a;
  abcd[62] = 2.0E0*I_NAI_Ix5z_F2xz_a;
  abcd[63] = 2.0E0*I_NAI_I6x_Fx2y_a-5*I_NAI_G4x_Fx2y;
  abcd[64] = 2.0E0*I_NAI_I5xy_Fx2y_a-4*I_NAI_G3xy_Fx2y;
  abcd[65] = 2.0E0*I_NAI_I5xz_Fx2y_a-4*I_NAI_G3xz_Fx2y;
  abcd[66] = 2.0E0*I_NAI_I4x2y_Fx2y_a-3*I_NAI_G2x2y_Fx2y;
  abcd[67] = 2.0E0*I_NAI_I4xyz_Fx2y_a-3*I_NAI_G2xyz_Fx2y;
  abcd[68] = 2.0E0*I_NAI_I4x2z_Fx2y_a-3*I_NAI_G2x2z_Fx2y;
  abcd[69] = 2.0E0*I_NAI_I3x3y_Fx2y_a-2*I_NAI_Gx3y_Fx2y;
  abcd[70] = 2.0E0*I_NAI_I3x2yz_Fx2y_a-2*I_NAI_Gx2yz_Fx2y;
  abcd[71] = 2.0E0*I_NAI_I3xy2z_Fx2y_a-2*I_NAI_Gxy2z_Fx2y;
  abcd[72] = 2.0E0*I_NAI_I3x3z_Fx2y_a-2*I_NAI_Gx3z_Fx2y;
  abcd[73] = 2.0E0*I_NAI_I2x4y_Fx2y_a-1*I_NAI_G4y_Fx2y;
  abcd[74] = 2.0E0*I_NAI_I2x3yz_Fx2y_a-1*I_NAI_G3yz_Fx2y;
  abcd[75] = 2.0E0*I_NAI_I2x2y2z_Fx2y_a-1*I_NAI_G2y2z_Fx2y;
  abcd[76] = 2.0E0*I_NAI_I2xy3z_Fx2y_a-1*I_NAI_Gy3z_Fx2y;
  abcd[77] = 2.0E0*I_NAI_I2x4z_Fx2y_a-1*I_NAI_G4z_Fx2y;
  abcd[78] = 2.0E0*I_NAI_Ix5y_Fx2y_a;
  abcd[79] = 2.0E0*I_NAI_Ix4yz_Fx2y_a;
  abcd[80] = 2.0E0*I_NAI_Ix3y2z_Fx2y_a;
  abcd[81] = 2.0E0*I_NAI_Ix2y3z_Fx2y_a;
  abcd[82] = 2.0E0*I_NAI_Ixy4z_Fx2y_a;
  abcd[83] = 2.0E0*I_NAI_Ix5z_Fx2y_a;
  abcd[84] = 2.0E0*I_NAI_I6x_Fxyz_a-5*I_NAI_G4x_Fxyz;
  abcd[85] = 2.0E0*I_NAI_I5xy_Fxyz_a-4*I_NAI_G3xy_Fxyz;
  abcd[86] = 2.0E0*I_NAI_I5xz_Fxyz_a-4*I_NAI_G3xz_Fxyz;
  abcd[87] = 2.0E0*I_NAI_I4x2y_Fxyz_a-3*I_NAI_G2x2y_Fxyz;
  abcd[88] = 2.0E0*I_NAI_I4xyz_Fxyz_a-3*I_NAI_G2xyz_Fxyz;
  abcd[89] = 2.0E0*I_NAI_I4x2z_Fxyz_a-3*I_NAI_G2x2z_Fxyz;
  abcd[90] = 2.0E0*I_NAI_I3x3y_Fxyz_a-2*I_NAI_Gx3y_Fxyz;
  abcd[91] = 2.0E0*I_NAI_I3x2yz_Fxyz_a-2*I_NAI_Gx2yz_Fxyz;
  abcd[92] = 2.0E0*I_NAI_I3xy2z_Fxyz_a-2*I_NAI_Gxy2z_Fxyz;
  abcd[93] = 2.0E0*I_NAI_I3x3z_Fxyz_a-2*I_NAI_Gx3z_Fxyz;
  abcd[94] = 2.0E0*I_NAI_I2x4y_Fxyz_a-1*I_NAI_G4y_Fxyz;
  abcd[95] = 2.0E0*I_NAI_I2x3yz_Fxyz_a-1*I_NAI_G3yz_Fxyz;
  abcd[96] = 2.0E0*I_NAI_I2x2y2z_Fxyz_a-1*I_NAI_G2y2z_Fxyz;
  abcd[97] = 2.0E0*I_NAI_I2xy3z_Fxyz_a-1*I_NAI_Gy3z_Fxyz;
  abcd[98] = 2.0E0*I_NAI_I2x4z_Fxyz_a-1*I_NAI_G4z_Fxyz;
  abcd[99] = 2.0E0*I_NAI_Ix5y_Fxyz_a;
  abcd[100] = 2.0E0*I_NAI_Ix4yz_Fxyz_a;
  abcd[101] = 2.0E0*I_NAI_Ix3y2z_Fxyz_a;
  abcd[102] = 2.0E0*I_NAI_Ix2y3z_Fxyz_a;
  abcd[103] = 2.0E0*I_NAI_Ixy4z_Fxyz_a;
  abcd[104] = 2.0E0*I_NAI_Ix5z_Fxyz_a;
  abcd[105] = 2.0E0*I_NAI_I6x_Fx2z_a-5*I_NAI_G4x_Fx2z;
  abcd[106] = 2.0E0*I_NAI_I5xy_Fx2z_a-4*I_NAI_G3xy_Fx2z;
  abcd[107] = 2.0E0*I_NAI_I5xz_Fx2z_a-4*I_NAI_G3xz_Fx2z;
  abcd[108] = 2.0E0*I_NAI_I4x2y_Fx2z_a-3*I_NAI_G2x2y_Fx2z;
  abcd[109] = 2.0E0*I_NAI_I4xyz_Fx2z_a-3*I_NAI_G2xyz_Fx2z;
  abcd[110] = 2.0E0*I_NAI_I4x2z_Fx2z_a-3*I_NAI_G2x2z_Fx2z;
  abcd[111] = 2.0E0*I_NAI_I3x3y_Fx2z_a-2*I_NAI_Gx3y_Fx2z;
  abcd[112] = 2.0E0*I_NAI_I3x2yz_Fx2z_a-2*I_NAI_Gx2yz_Fx2z;
  abcd[113] = 2.0E0*I_NAI_I3xy2z_Fx2z_a-2*I_NAI_Gxy2z_Fx2z;
  abcd[114] = 2.0E0*I_NAI_I3x3z_Fx2z_a-2*I_NAI_Gx3z_Fx2z;
  abcd[115] = 2.0E0*I_NAI_I2x4y_Fx2z_a-1*I_NAI_G4y_Fx2z;
  abcd[116] = 2.0E0*I_NAI_I2x3yz_Fx2z_a-1*I_NAI_G3yz_Fx2z;
  abcd[117] = 2.0E0*I_NAI_I2x2y2z_Fx2z_a-1*I_NAI_G2y2z_Fx2z;
  abcd[118] = 2.0E0*I_NAI_I2xy3z_Fx2z_a-1*I_NAI_Gy3z_Fx2z;
  abcd[119] = 2.0E0*I_NAI_I2x4z_Fx2z_a-1*I_NAI_G4z_Fx2z;
  abcd[120] = 2.0E0*I_NAI_Ix5y_Fx2z_a;
  abcd[121] = 2.0E0*I_NAI_Ix4yz_Fx2z_a;
  abcd[122] = 2.0E0*I_NAI_Ix3y2z_Fx2z_a;
  abcd[123] = 2.0E0*I_NAI_Ix2y3z_Fx2z_a;
  abcd[124] = 2.0E0*I_NAI_Ixy4z_Fx2z_a;
  abcd[125] = 2.0E0*I_NAI_Ix5z_Fx2z_a;
  abcd[126] = 2.0E0*I_NAI_I6x_F3y_a-5*I_NAI_G4x_F3y;
  abcd[127] = 2.0E0*I_NAI_I5xy_F3y_a-4*I_NAI_G3xy_F3y;
  abcd[128] = 2.0E0*I_NAI_I5xz_F3y_a-4*I_NAI_G3xz_F3y;
  abcd[129] = 2.0E0*I_NAI_I4x2y_F3y_a-3*I_NAI_G2x2y_F3y;
  abcd[130] = 2.0E0*I_NAI_I4xyz_F3y_a-3*I_NAI_G2xyz_F3y;
  abcd[131] = 2.0E0*I_NAI_I4x2z_F3y_a-3*I_NAI_G2x2z_F3y;
  abcd[132] = 2.0E0*I_NAI_I3x3y_F3y_a-2*I_NAI_Gx3y_F3y;
  abcd[133] = 2.0E0*I_NAI_I3x2yz_F3y_a-2*I_NAI_Gx2yz_F3y;
  abcd[134] = 2.0E0*I_NAI_I3xy2z_F3y_a-2*I_NAI_Gxy2z_F3y;
  abcd[135] = 2.0E0*I_NAI_I3x3z_F3y_a-2*I_NAI_Gx3z_F3y;
  abcd[136] = 2.0E0*I_NAI_I2x4y_F3y_a-1*I_NAI_G4y_F3y;
  abcd[137] = 2.0E0*I_NAI_I2x3yz_F3y_a-1*I_NAI_G3yz_F3y;
  abcd[138] = 2.0E0*I_NAI_I2x2y2z_F3y_a-1*I_NAI_G2y2z_F3y;
  abcd[139] = 2.0E0*I_NAI_I2xy3z_F3y_a-1*I_NAI_Gy3z_F3y;
  abcd[140] = 2.0E0*I_NAI_I2x4z_F3y_a-1*I_NAI_G4z_F3y;
  abcd[141] = 2.0E0*I_NAI_Ix5y_F3y_a;
  abcd[142] = 2.0E0*I_NAI_Ix4yz_F3y_a;
  abcd[143] = 2.0E0*I_NAI_Ix3y2z_F3y_a;
  abcd[144] = 2.0E0*I_NAI_Ix2y3z_F3y_a;
  abcd[145] = 2.0E0*I_NAI_Ixy4z_F3y_a;
  abcd[146] = 2.0E0*I_NAI_Ix5z_F3y_a;
  abcd[147] = 2.0E0*I_NAI_I6x_F2yz_a-5*I_NAI_G4x_F2yz;
  abcd[148] = 2.0E0*I_NAI_I5xy_F2yz_a-4*I_NAI_G3xy_F2yz;
  abcd[149] = 2.0E0*I_NAI_I5xz_F2yz_a-4*I_NAI_G3xz_F2yz;
  abcd[150] = 2.0E0*I_NAI_I4x2y_F2yz_a-3*I_NAI_G2x2y_F2yz;
  abcd[151] = 2.0E0*I_NAI_I4xyz_F2yz_a-3*I_NAI_G2xyz_F2yz;
  abcd[152] = 2.0E0*I_NAI_I4x2z_F2yz_a-3*I_NAI_G2x2z_F2yz;
  abcd[153] = 2.0E0*I_NAI_I3x3y_F2yz_a-2*I_NAI_Gx3y_F2yz;
  abcd[154] = 2.0E0*I_NAI_I3x2yz_F2yz_a-2*I_NAI_Gx2yz_F2yz;
  abcd[155] = 2.0E0*I_NAI_I3xy2z_F2yz_a-2*I_NAI_Gxy2z_F2yz;
  abcd[156] = 2.0E0*I_NAI_I3x3z_F2yz_a-2*I_NAI_Gx3z_F2yz;
  abcd[157] = 2.0E0*I_NAI_I2x4y_F2yz_a-1*I_NAI_G4y_F2yz;
  abcd[158] = 2.0E0*I_NAI_I2x3yz_F2yz_a-1*I_NAI_G3yz_F2yz;
  abcd[159] = 2.0E0*I_NAI_I2x2y2z_F2yz_a-1*I_NAI_G2y2z_F2yz;
  abcd[160] = 2.0E0*I_NAI_I2xy3z_F2yz_a-1*I_NAI_Gy3z_F2yz;
  abcd[161] = 2.0E0*I_NAI_I2x4z_F2yz_a-1*I_NAI_G4z_F2yz;
  abcd[162] = 2.0E0*I_NAI_Ix5y_F2yz_a;
  abcd[163] = 2.0E0*I_NAI_Ix4yz_F2yz_a;
  abcd[164] = 2.0E0*I_NAI_Ix3y2z_F2yz_a;
  abcd[165] = 2.0E0*I_NAI_Ix2y3z_F2yz_a;
  abcd[166] = 2.0E0*I_NAI_Ixy4z_F2yz_a;
  abcd[167] = 2.0E0*I_NAI_Ix5z_F2yz_a;
  abcd[168] = 2.0E0*I_NAI_I6x_Fy2z_a-5*I_NAI_G4x_Fy2z;
  abcd[169] = 2.0E0*I_NAI_I5xy_Fy2z_a-4*I_NAI_G3xy_Fy2z;
  abcd[170] = 2.0E0*I_NAI_I5xz_Fy2z_a-4*I_NAI_G3xz_Fy2z;
  abcd[171] = 2.0E0*I_NAI_I4x2y_Fy2z_a-3*I_NAI_G2x2y_Fy2z;
  abcd[172] = 2.0E0*I_NAI_I4xyz_Fy2z_a-3*I_NAI_G2xyz_Fy2z;
  abcd[173] = 2.0E0*I_NAI_I4x2z_Fy2z_a-3*I_NAI_G2x2z_Fy2z;
  abcd[174] = 2.0E0*I_NAI_I3x3y_Fy2z_a-2*I_NAI_Gx3y_Fy2z;
  abcd[175] = 2.0E0*I_NAI_I3x2yz_Fy2z_a-2*I_NAI_Gx2yz_Fy2z;
  abcd[176] = 2.0E0*I_NAI_I3xy2z_Fy2z_a-2*I_NAI_Gxy2z_Fy2z;
  abcd[177] = 2.0E0*I_NAI_I3x3z_Fy2z_a-2*I_NAI_Gx3z_Fy2z;
  abcd[178] = 2.0E0*I_NAI_I2x4y_Fy2z_a-1*I_NAI_G4y_Fy2z;
  abcd[179] = 2.0E0*I_NAI_I2x3yz_Fy2z_a-1*I_NAI_G3yz_Fy2z;
  abcd[180] = 2.0E0*I_NAI_I2x2y2z_Fy2z_a-1*I_NAI_G2y2z_Fy2z;
  abcd[181] = 2.0E0*I_NAI_I2xy3z_Fy2z_a-1*I_NAI_Gy3z_Fy2z;
  abcd[182] = 2.0E0*I_NAI_I2x4z_Fy2z_a-1*I_NAI_G4z_Fy2z;
  abcd[183] = 2.0E0*I_NAI_Ix5y_Fy2z_a;
  abcd[184] = 2.0E0*I_NAI_Ix4yz_Fy2z_a;
  abcd[185] = 2.0E0*I_NAI_Ix3y2z_Fy2z_a;
  abcd[186] = 2.0E0*I_NAI_Ix2y3z_Fy2z_a;
  abcd[187] = 2.0E0*I_NAI_Ixy4z_Fy2z_a;
  abcd[188] = 2.0E0*I_NAI_Ix5z_Fy2z_a;
  abcd[189] = 2.0E0*I_NAI_I6x_F3z_a-5*I_NAI_G4x_F3z;
  abcd[190] = 2.0E0*I_NAI_I5xy_F3z_a-4*I_NAI_G3xy_F3z;
  abcd[191] = 2.0E0*I_NAI_I5xz_F3z_a-4*I_NAI_G3xz_F3z;
  abcd[192] = 2.0E0*I_NAI_I4x2y_F3z_a-3*I_NAI_G2x2y_F3z;
  abcd[193] = 2.0E0*I_NAI_I4xyz_F3z_a-3*I_NAI_G2xyz_F3z;
  abcd[194] = 2.0E0*I_NAI_I4x2z_F3z_a-3*I_NAI_G2x2z_F3z;
  abcd[195] = 2.0E0*I_NAI_I3x3y_F3z_a-2*I_NAI_Gx3y_F3z;
  abcd[196] = 2.0E0*I_NAI_I3x2yz_F3z_a-2*I_NAI_Gx2yz_F3z;
  abcd[197] = 2.0E0*I_NAI_I3xy2z_F3z_a-2*I_NAI_Gxy2z_F3z;
  abcd[198] = 2.0E0*I_NAI_I3x3z_F3z_a-2*I_NAI_Gx3z_F3z;
  abcd[199] = 2.0E0*I_NAI_I2x4y_F3z_a-1*I_NAI_G4y_F3z;
  abcd[200] = 2.0E0*I_NAI_I2x3yz_F3z_a-1*I_NAI_G3yz_F3z;
  abcd[201] = 2.0E0*I_NAI_I2x2y2z_F3z_a-1*I_NAI_G2y2z_F3z;
  abcd[202] = 2.0E0*I_NAI_I2xy3z_F3z_a-1*I_NAI_Gy3z_F3z;
  abcd[203] = 2.0E0*I_NAI_I2x4z_F3z_a-1*I_NAI_G4z_F3z;
  abcd[204] = 2.0E0*I_NAI_Ix5y_F3z_a;
  abcd[205] = 2.0E0*I_NAI_Ix4yz_F3z_a;
  abcd[206] = 2.0E0*I_NAI_Ix3y2z_F3z_a;
  abcd[207] = 2.0E0*I_NAI_Ix2y3z_F3z_a;
  abcd[208] = 2.0E0*I_NAI_Ixy4z_F3z_a;
  abcd[209] = 2.0E0*I_NAI_Ix5z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_a
   * RHS shell quartet name: SQ_NAI_G_F
   ************************************************************/
  abcd[210] = 2.0E0*I_NAI_I5xy_F3x_a;
  abcd[211] = 2.0E0*I_NAI_I4x2y_F3x_a-1*I_NAI_G4x_F3x;
  abcd[212] = 2.0E0*I_NAI_I4xyz_F3x_a;
  abcd[213] = 2.0E0*I_NAI_I3x3y_F3x_a-2*I_NAI_G3xy_F3x;
  abcd[214] = 2.0E0*I_NAI_I3x2yz_F3x_a-1*I_NAI_G3xz_F3x;
  abcd[215] = 2.0E0*I_NAI_I3xy2z_F3x_a;
  abcd[216] = 2.0E0*I_NAI_I2x4y_F3x_a-3*I_NAI_G2x2y_F3x;
  abcd[217] = 2.0E0*I_NAI_I2x3yz_F3x_a-2*I_NAI_G2xyz_F3x;
  abcd[218] = 2.0E0*I_NAI_I2x2y2z_F3x_a-1*I_NAI_G2x2z_F3x;
  abcd[219] = 2.0E0*I_NAI_I2xy3z_F3x_a;
  abcd[220] = 2.0E0*I_NAI_Ix5y_F3x_a-4*I_NAI_Gx3y_F3x;
  abcd[221] = 2.0E0*I_NAI_Ix4yz_F3x_a-3*I_NAI_Gx2yz_F3x;
  abcd[222] = 2.0E0*I_NAI_Ix3y2z_F3x_a-2*I_NAI_Gxy2z_F3x;
  abcd[223] = 2.0E0*I_NAI_Ix2y3z_F3x_a-1*I_NAI_Gx3z_F3x;
  abcd[224] = 2.0E0*I_NAI_Ixy4z_F3x_a;
  abcd[225] = 2.0E0*I_NAI_I6y_F3x_a-5*I_NAI_G4y_F3x;
  abcd[226] = 2.0E0*I_NAI_I5yz_F3x_a-4*I_NAI_G3yz_F3x;
  abcd[227] = 2.0E0*I_NAI_I4y2z_F3x_a-3*I_NAI_G2y2z_F3x;
  abcd[228] = 2.0E0*I_NAI_I3y3z_F3x_a-2*I_NAI_Gy3z_F3x;
  abcd[229] = 2.0E0*I_NAI_I2y4z_F3x_a-1*I_NAI_G4z_F3x;
  abcd[230] = 2.0E0*I_NAI_Iy5z_F3x_a;
  abcd[231] = 2.0E0*I_NAI_I5xy_F2xy_a;
  abcd[232] = 2.0E0*I_NAI_I4x2y_F2xy_a-1*I_NAI_G4x_F2xy;
  abcd[233] = 2.0E0*I_NAI_I4xyz_F2xy_a;
  abcd[234] = 2.0E0*I_NAI_I3x3y_F2xy_a-2*I_NAI_G3xy_F2xy;
  abcd[235] = 2.0E0*I_NAI_I3x2yz_F2xy_a-1*I_NAI_G3xz_F2xy;
  abcd[236] = 2.0E0*I_NAI_I3xy2z_F2xy_a;
  abcd[237] = 2.0E0*I_NAI_I2x4y_F2xy_a-3*I_NAI_G2x2y_F2xy;
  abcd[238] = 2.0E0*I_NAI_I2x3yz_F2xy_a-2*I_NAI_G2xyz_F2xy;
  abcd[239] = 2.0E0*I_NAI_I2x2y2z_F2xy_a-1*I_NAI_G2x2z_F2xy;
  abcd[240] = 2.0E0*I_NAI_I2xy3z_F2xy_a;
  abcd[241] = 2.0E0*I_NAI_Ix5y_F2xy_a-4*I_NAI_Gx3y_F2xy;
  abcd[242] = 2.0E0*I_NAI_Ix4yz_F2xy_a-3*I_NAI_Gx2yz_F2xy;
  abcd[243] = 2.0E0*I_NAI_Ix3y2z_F2xy_a-2*I_NAI_Gxy2z_F2xy;
  abcd[244] = 2.0E0*I_NAI_Ix2y3z_F2xy_a-1*I_NAI_Gx3z_F2xy;
  abcd[245] = 2.0E0*I_NAI_Ixy4z_F2xy_a;
  abcd[246] = 2.0E0*I_NAI_I6y_F2xy_a-5*I_NAI_G4y_F2xy;
  abcd[247] = 2.0E0*I_NAI_I5yz_F2xy_a-4*I_NAI_G3yz_F2xy;
  abcd[248] = 2.0E0*I_NAI_I4y2z_F2xy_a-3*I_NAI_G2y2z_F2xy;
  abcd[249] = 2.0E0*I_NAI_I3y3z_F2xy_a-2*I_NAI_Gy3z_F2xy;
  abcd[250] = 2.0E0*I_NAI_I2y4z_F2xy_a-1*I_NAI_G4z_F2xy;
  abcd[251] = 2.0E0*I_NAI_Iy5z_F2xy_a;
  abcd[252] = 2.0E0*I_NAI_I5xy_F2xz_a;
  abcd[253] = 2.0E0*I_NAI_I4x2y_F2xz_a-1*I_NAI_G4x_F2xz;
  abcd[254] = 2.0E0*I_NAI_I4xyz_F2xz_a;
  abcd[255] = 2.0E0*I_NAI_I3x3y_F2xz_a-2*I_NAI_G3xy_F2xz;
  abcd[256] = 2.0E0*I_NAI_I3x2yz_F2xz_a-1*I_NAI_G3xz_F2xz;
  abcd[257] = 2.0E0*I_NAI_I3xy2z_F2xz_a;
  abcd[258] = 2.0E0*I_NAI_I2x4y_F2xz_a-3*I_NAI_G2x2y_F2xz;
  abcd[259] = 2.0E0*I_NAI_I2x3yz_F2xz_a-2*I_NAI_G2xyz_F2xz;
  abcd[260] = 2.0E0*I_NAI_I2x2y2z_F2xz_a-1*I_NAI_G2x2z_F2xz;
  abcd[261] = 2.0E0*I_NAI_I2xy3z_F2xz_a;
  abcd[262] = 2.0E0*I_NAI_Ix5y_F2xz_a-4*I_NAI_Gx3y_F2xz;
  abcd[263] = 2.0E0*I_NAI_Ix4yz_F2xz_a-3*I_NAI_Gx2yz_F2xz;
  abcd[264] = 2.0E0*I_NAI_Ix3y2z_F2xz_a-2*I_NAI_Gxy2z_F2xz;
  abcd[265] = 2.0E0*I_NAI_Ix2y3z_F2xz_a-1*I_NAI_Gx3z_F2xz;
  abcd[266] = 2.0E0*I_NAI_Ixy4z_F2xz_a;
  abcd[267] = 2.0E0*I_NAI_I6y_F2xz_a-5*I_NAI_G4y_F2xz;
  abcd[268] = 2.0E0*I_NAI_I5yz_F2xz_a-4*I_NAI_G3yz_F2xz;
  abcd[269] = 2.0E0*I_NAI_I4y2z_F2xz_a-3*I_NAI_G2y2z_F2xz;
  abcd[270] = 2.0E0*I_NAI_I3y3z_F2xz_a-2*I_NAI_Gy3z_F2xz;
  abcd[271] = 2.0E0*I_NAI_I2y4z_F2xz_a-1*I_NAI_G4z_F2xz;
  abcd[272] = 2.0E0*I_NAI_Iy5z_F2xz_a;
  abcd[273] = 2.0E0*I_NAI_I5xy_Fx2y_a;
  abcd[274] = 2.0E0*I_NAI_I4x2y_Fx2y_a-1*I_NAI_G4x_Fx2y;
  abcd[275] = 2.0E0*I_NAI_I4xyz_Fx2y_a;
  abcd[276] = 2.0E0*I_NAI_I3x3y_Fx2y_a-2*I_NAI_G3xy_Fx2y;
  abcd[277] = 2.0E0*I_NAI_I3x2yz_Fx2y_a-1*I_NAI_G3xz_Fx2y;
  abcd[278] = 2.0E0*I_NAI_I3xy2z_Fx2y_a;
  abcd[279] = 2.0E0*I_NAI_I2x4y_Fx2y_a-3*I_NAI_G2x2y_Fx2y;
  abcd[280] = 2.0E0*I_NAI_I2x3yz_Fx2y_a-2*I_NAI_G2xyz_Fx2y;
  abcd[281] = 2.0E0*I_NAI_I2x2y2z_Fx2y_a-1*I_NAI_G2x2z_Fx2y;
  abcd[282] = 2.0E0*I_NAI_I2xy3z_Fx2y_a;
  abcd[283] = 2.0E0*I_NAI_Ix5y_Fx2y_a-4*I_NAI_Gx3y_Fx2y;
  abcd[284] = 2.0E0*I_NAI_Ix4yz_Fx2y_a-3*I_NAI_Gx2yz_Fx2y;
  abcd[285] = 2.0E0*I_NAI_Ix3y2z_Fx2y_a-2*I_NAI_Gxy2z_Fx2y;
  abcd[286] = 2.0E0*I_NAI_Ix2y3z_Fx2y_a-1*I_NAI_Gx3z_Fx2y;
  abcd[287] = 2.0E0*I_NAI_Ixy4z_Fx2y_a;
  abcd[288] = 2.0E0*I_NAI_I6y_Fx2y_a-5*I_NAI_G4y_Fx2y;
  abcd[289] = 2.0E0*I_NAI_I5yz_Fx2y_a-4*I_NAI_G3yz_Fx2y;
  abcd[290] = 2.0E0*I_NAI_I4y2z_Fx2y_a-3*I_NAI_G2y2z_Fx2y;
  abcd[291] = 2.0E0*I_NAI_I3y3z_Fx2y_a-2*I_NAI_Gy3z_Fx2y;
  abcd[292] = 2.0E0*I_NAI_I2y4z_Fx2y_a-1*I_NAI_G4z_Fx2y;
  abcd[293] = 2.0E0*I_NAI_Iy5z_Fx2y_a;
  abcd[294] = 2.0E0*I_NAI_I5xy_Fxyz_a;
  abcd[295] = 2.0E0*I_NAI_I4x2y_Fxyz_a-1*I_NAI_G4x_Fxyz;
  abcd[296] = 2.0E0*I_NAI_I4xyz_Fxyz_a;
  abcd[297] = 2.0E0*I_NAI_I3x3y_Fxyz_a-2*I_NAI_G3xy_Fxyz;
  abcd[298] = 2.0E0*I_NAI_I3x2yz_Fxyz_a-1*I_NAI_G3xz_Fxyz;
  abcd[299] = 2.0E0*I_NAI_I3xy2z_Fxyz_a;
  abcd[300] = 2.0E0*I_NAI_I2x4y_Fxyz_a-3*I_NAI_G2x2y_Fxyz;
  abcd[301] = 2.0E0*I_NAI_I2x3yz_Fxyz_a-2*I_NAI_G2xyz_Fxyz;
  abcd[302] = 2.0E0*I_NAI_I2x2y2z_Fxyz_a-1*I_NAI_G2x2z_Fxyz;
  abcd[303] = 2.0E0*I_NAI_I2xy3z_Fxyz_a;
  abcd[304] = 2.0E0*I_NAI_Ix5y_Fxyz_a-4*I_NAI_Gx3y_Fxyz;
  abcd[305] = 2.0E0*I_NAI_Ix4yz_Fxyz_a-3*I_NAI_Gx2yz_Fxyz;
  abcd[306] = 2.0E0*I_NAI_Ix3y2z_Fxyz_a-2*I_NAI_Gxy2z_Fxyz;
  abcd[307] = 2.0E0*I_NAI_Ix2y3z_Fxyz_a-1*I_NAI_Gx3z_Fxyz;
  abcd[308] = 2.0E0*I_NAI_Ixy4z_Fxyz_a;
  abcd[309] = 2.0E0*I_NAI_I6y_Fxyz_a-5*I_NAI_G4y_Fxyz;
  abcd[310] = 2.0E0*I_NAI_I5yz_Fxyz_a-4*I_NAI_G3yz_Fxyz;
  abcd[311] = 2.0E0*I_NAI_I4y2z_Fxyz_a-3*I_NAI_G2y2z_Fxyz;
  abcd[312] = 2.0E0*I_NAI_I3y3z_Fxyz_a-2*I_NAI_Gy3z_Fxyz;
  abcd[313] = 2.0E0*I_NAI_I2y4z_Fxyz_a-1*I_NAI_G4z_Fxyz;
  abcd[314] = 2.0E0*I_NAI_Iy5z_Fxyz_a;
  abcd[315] = 2.0E0*I_NAI_I5xy_Fx2z_a;
  abcd[316] = 2.0E0*I_NAI_I4x2y_Fx2z_a-1*I_NAI_G4x_Fx2z;
  abcd[317] = 2.0E0*I_NAI_I4xyz_Fx2z_a;
  abcd[318] = 2.0E0*I_NAI_I3x3y_Fx2z_a-2*I_NAI_G3xy_Fx2z;
  abcd[319] = 2.0E0*I_NAI_I3x2yz_Fx2z_a-1*I_NAI_G3xz_Fx2z;
  abcd[320] = 2.0E0*I_NAI_I3xy2z_Fx2z_a;
  abcd[321] = 2.0E0*I_NAI_I2x4y_Fx2z_a-3*I_NAI_G2x2y_Fx2z;
  abcd[322] = 2.0E0*I_NAI_I2x3yz_Fx2z_a-2*I_NAI_G2xyz_Fx2z;
  abcd[323] = 2.0E0*I_NAI_I2x2y2z_Fx2z_a-1*I_NAI_G2x2z_Fx2z;
  abcd[324] = 2.0E0*I_NAI_I2xy3z_Fx2z_a;
  abcd[325] = 2.0E0*I_NAI_Ix5y_Fx2z_a-4*I_NAI_Gx3y_Fx2z;
  abcd[326] = 2.0E0*I_NAI_Ix4yz_Fx2z_a-3*I_NAI_Gx2yz_Fx2z;
  abcd[327] = 2.0E0*I_NAI_Ix3y2z_Fx2z_a-2*I_NAI_Gxy2z_Fx2z;
  abcd[328] = 2.0E0*I_NAI_Ix2y3z_Fx2z_a-1*I_NAI_Gx3z_Fx2z;
  abcd[329] = 2.0E0*I_NAI_Ixy4z_Fx2z_a;
  abcd[330] = 2.0E0*I_NAI_I6y_Fx2z_a-5*I_NAI_G4y_Fx2z;
  abcd[331] = 2.0E0*I_NAI_I5yz_Fx2z_a-4*I_NAI_G3yz_Fx2z;
  abcd[332] = 2.0E0*I_NAI_I4y2z_Fx2z_a-3*I_NAI_G2y2z_Fx2z;
  abcd[333] = 2.0E0*I_NAI_I3y3z_Fx2z_a-2*I_NAI_Gy3z_Fx2z;
  abcd[334] = 2.0E0*I_NAI_I2y4z_Fx2z_a-1*I_NAI_G4z_Fx2z;
  abcd[335] = 2.0E0*I_NAI_Iy5z_Fx2z_a;
  abcd[336] = 2.0E0*I_NAI_I5xy_F3y_a;
  abcd[337] = 2.0E0*I_NAI_I4x2y_F3y_a-1*I_NAI_G4x_F3y;
  abcd[338] = 2.0E0*I_NAI_I4xyz_F3y_a;
  abcd[339] = 2.0E0*I_NAI_I3x3y_F3y_a-2*I_NAI_G3xy_F3y;
  abcd[340] = 2.0E0*I_NAI_I3x2yz_F3y_a-1*I_NAI_G3xz_F3y;
  abcd[341] = 2.0E0*I_NAI_I3xy2z_F3y_a;
  abcd[342] = 2.0E0*I_NAI_I2x4y_F3y_a-3*I_NAI_G2x2y_F3y;
  abcd[343] = 2.0E0*I_NAI_I2x3yz_F3y_a-2*I_NAI_G2xyz_F3y;
  abcd[344] = 2.0E0*I_NAI_I2x2y2z_F3y_a-1*I_NAI_G2x2z_F3y;
  abcd[345] = 2.0E0*I_NAI_I2xy3z_F3y_a;
  abcd[346] = 2.0E0*I_NAI_Ix5y_F3y_a-4*I_NAI_Gx3y_F3y;
  abcd[347] = 2.0E0*I_NAI_Ix4yz_F3y_a-3*I_NAI_Gx2yz_F3y;
  abcd[348] = 2.0E0*I_NAI_Ix3y2z_F3y_a-2*I_NAI_Gxy2z_F3y;
  abcd[349] = 2.0E0*I_NAI_Ix2y3z_F3y_a-1*I_NAI_Gx3z_F3y;
  abcd[350] = 2.0E0*I_NAI_Ixy4z_F3y_a;
  abcd[351] = 2.0E0*I_NAI_I6y_F3y_a-5*I_NAI_G4y_F3y;
  abcd[352] = 2.0E0*I_NAI_I5yz_F3y_a-4*I_NAI_G3yz_F3y;
  abcd[353] = 2.0E0*I_NAI_I4y2z_F3y_a-3*I_NAI_G2y2z_F3y;
  abcd[354] = 2.0E0*I_NAI_I3y3z_F3y_a-2*I_NAI_Gy3z_F3y;
  abcd[355] = 2.0E0*I_NAI_I2y4z_F3y_a-1*I_NAI_G4z_F3y;
  abcd[356] = 2.0E0*I_NAI_Iy5z_F3y_a;
  abcd[357] = 2.0E0*I_NAI_I5xy_F2yz_a;
  abcd[358] = 2.0E0*I_NAI_I4x2y_F2yz_a-1*I_NAI_G4x_F2yz;
  abcd[359] = 2.0E0*I_NAI_I4xyz_F2yz_a;
  abcd[360] = 2.0E0*I_NAI_I3x3y_F2yz_a-2*I_NAI_G3xy_F2yz;
  abcd[361] = 2.0E0*I_NAI_I3x2yz_F2yz_a-1*I_NAI_G3xz_F2yz;
  abcd[362] = 2.0E0*I_NAI_I3xy2z_F2yz_a;
  abcd[363] = 2.0E0*I_NAI_I2x4y_F2yz_a-3*I_NAI_G2x2y_F2yz;
  abcd[364] = 2.0E0*I_NAI_I2x3yz_F2yz_a-2*I_NAI_G2xyz_F2yz;
  abcd[365] = 2.0E0*I_NAI_I2x2y2z_F2yz_a-1*I_NAI_G2x2z_F2yz;
  abcd[366] = 2.0E0*I_NAI_I2xy3z_F2yz_a;
  abcd[367] = 2.0E0*I_NAI_Ix5y_F2yz_a-4*I_NAI_Gx3y_F2yz;
  abcd[368] = 2.0E0*I_NAI_Ix4yz_F2yz_a-3*I_NAI_Gx2yz_F2yz;
  abcd[369] = 2.0E0*I_NAI_Ix3y2z_F2yz_a-2*I_NAI_Gxy2z_F2yz;
  abcd[370] = 2.0E0*I_NAI_Ix2y3z_F2yz_a-1*I_NAI_Gx3z_F2yz;
  abcd[371] = 2.0E0*I_NAI_Ixy4z_F2yz_a;
  abcd[372] = 2.0E0*I_NAI_I6y_F2yz_a-5*I_NAI_G4y_F2yz;
  abcd[373] = 2.0E0*I_NAI_I5yz_F2yz_a-4*I_NAI_G3yz_F2yz;
  abcd[374] = 2.0E0*I_NAI_I4y2z_F2yz_a-3*I_NAI_G2y2z_F2yz;
  abcd[375] = 2.0E0*I_NAI_I3y3z_F2yz_a-2*I_NAI_Gy3z_F2yz;
  abcd[376] = 2.0E0*I_NAI_I2y4z_F2yz_a-1*I_NAI_G4z_F2yz;
  abcd[377] = 2.0E0*I_NAI_Iy5z_F2yz_a;
  abcd[378] = 2.0E0*I_NAI_I5xy_Fy2z_a;
  abcd[379] = 2.0E0*I_NAI_I4x2y_Fy2z_a-1*I_NAI_G4x_Fy2z;
  abcd[380] = 2.0E0*I_NAI_I4xyz_Fy2z_a;
  abcd[381] = 2.0E0*I_NAI_I3x3y_Fy2z_a-2*I_NAI_G3xy_Fy2z;
  abcd[382] = 2.0E0*I_NAI_I3x2yz_Fy2z_a-1*I_NAI_G3xz_Fy2z;
  abcd[383] = 2.0E0*I_NAI_I3xy2z_Fy2z_a;
  abcd[384] = 2.0E0*I_NAI_I2x4y_Fy2z_a-3*I_NAI_G2x2y_Fy2z;
  abcd[385] = 2.0E0*I_NAI_I2x3yz_Fy2z_a-2*I_NAI_G2xyz_Fy2z;
  abcd[386] = 2.0E0*I_NAI_I2x2y2z_Fy2z_a-1*I_NAI_G2x2z_Fy2z;
  abcd[387] = 2.0E0*I_NAI_I2xy3z_Fy2z_a;
  abcd[388] = 2.0E0*I_NAI_Ix5y_Fy2z_a-4*I_NAI_Gx3y_Fy2z;
  abcd[389] = 2.0E0*I_NAI_Ix4yz_Fy2z_a-3*I_NAI_Gx2yz_Fy2z;
  abcd[390] = 2.0E0*I_NAI_Ix3y2z_Fy2z_a-2*I_NAI_Gxy2z_Fy2z;
  abcd[391] = 2.0E0*I_NAI_Ix2y3z_Fy2z_a-1*I_NAI_Gx3z_Fy2z;
  abcd[392] = 2.0E0*I_NAI_Ixy4z_Fy2z_a;
  abcd[393] = 2.0E0*I_NAI_I6y_Fy2z_a-5*I_NAI_G4y_Fy2z;
  abcd[394] = 2.0E0*I_NAI_I5yz_Fy2z_a-4*I_NAI_G3yz_Fy2z;
  abcd[395] = 2.0E0*I_NAI_I4y2z_Fy2z_a-3*I_NAI_G2y2z_Fy2z;
  abcd[396] = 2.0E0*I_NAI_I3y3z_Fy2z_a-2*I_NAI_Gy3z_Fy2z;
  abcd[397] = 2.0E0*I_NAI_I2y4z_Fy2z_a-1*I_NAI_G4z_Fy2z;
  abcd[398] = 2.0E0*I_NAI_Iy5z_Fy2z_a;
  abcd[399] = 2.0E0*I_NAI_I5xy_F3z_a;
  abcd[400] = 2.0E0*I_NAI_I4x2y_F3z_a-1*I_NAI_G4x_F3z;
  abcd[401] = 2.0E0*I_NAI_I4xyz_F3z_a;
  abcd[402] = 2.0E0*I_NAI_I3x3y_F3z_a-2*I_NAI_G3xy_F3z;
  abcd[403] = 2.0E0*I_NAI_I3x2yz_F3z_a-1*I_NAI_G3xz_F3z;
  abcd[404] = 2.0E0*I_NAI_I3xy2z_F3z_a;
  abcd[405] = 2.0E0*I_NAI_I2x4y_F3z_a-3*I_NAI_G2x2y_F3z;
  abcd[406] = 2.0E0*I_NAI_I2x3yz_F3z_a-2*I_NAI_G2xyz_F3z;
  abcd[407] = 2.0E0*I_NAI_I2x2y2z_F3z_a-1*I_NAI_G2x2z_F3z;
  abcd[408] = 2.0E0*I_NAI_I2xy3z_F3z_a;
  abcd[409] = 2.0E0*I_NAI_Ix5y_F3z_a-4*I_NAI_Gx3y_F3z;
  abcd[410] = 2.0E0*I_NAI_Ix4yz_F3z_a-3*I_NAI_Gx2yz_F3z;
  abcd[411] = 2.0E0*I_NAI_Ix3y2z_F3z_a-2*I_NAI_Gxy2z_F3z;
  abcd[412] = 2.0E0*I_NAI_Ix2y3z_F3z_a-1*I_NAI_Gx3z_F3z;
  abcd[413] = 2.0E0*I_NAI_Ixy4z_F3z_a;
  abcd[414] = 2.0E0*I_NAI_I6y_F3z_a-5*I_NAI_G4y_F3z;
  abcd[415] = 2.0E0*I_NAI_I5yz_F3z_a-4*I_NAI_G3yz_F3z;
  abcd[416] = 2.0E0*I_NAI_I4y2z_F3z_a-3*I_NAI_G2y2z_F3z;
  abcd[417] = 2.0E0*I_NAI_I3y3z_F3z_a-2*I_NAI_Gy3z_F3z;
  abcd[418] = 2.0E0*I_NAI_I2y4z_F3z_a-1*I_NAI_G4z_F3z;
  abcd[419] = 2.0E0*I_NAI_Iy5z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_F_a
   * RHS shell quartet name: SQ_NAI_G_F
   ************************************************************/
  abcd[420] = 2.0E0*I_NAI_I5xz_F3x_a;
  abcd[421] = 2.0E0*I_NAI_I4xyz_F3x_a;
  abcd[422] = 2.0E0*I_NAI_I4x2z_F3x_a-1*I_NAI_G4x_F3x;
  abcd[423] = 2.0E0*I_NAI_I3x2yz_F3x_a;
  abcd[424] = 2.0E0*I_NAI_I3xy2z_F3x_a-1*I_NAI_G3xy_F3x;
  abcd[425] = 2.0E0*I_NAI_I3x3z_F3x_a-2*I_NAI_G3xz_F3x;
  abcd[426] = 2.0E0*I_NAI_I2x3yz_F3x_a;
  abcd[427] = 2.0E0*I_NAI_I2x2y2z_F3x_a-1*I_NAI_G2x2y_F3x;
  abcd[428] = 2.0E0*I_NAI_I2xy3z_F3x_a-2*I_NAI_G2xyz_F3x;
  abcd[429] = 2.0E0*I_NAI_I2x4z_F3x_a-3*I_NAI_G2x2z_F3x;
  abcd[430] = 2.0E0*I_NAI_Ix4yz_F3x_a;
  abcd[431] = 2.0E0*I_NAI_Ix3y2z_F3x_a-1*I_NAI_Gx3y_F3x;
  abcd[432] = 2.0E0*I_NAI_Ix2y3z_F3x_a-2*I_NAI_Gx2yz_F3x;
  abcd[433] = 2.0E0*I_NAI_Ixy4z_F3x_a-3*I_NAI_Gxy2z_F3x;
  abcd[434] = 2.0E0*I_NAI_Ix5z_F3x_a-4*I_NAI_Gx3z_F3x;
  abcd[435] = 2.0E0*I_NAI_I5yz_F3x_a;
  abcd[436] = 2.0E0*I_NAI_I4y2z_F3x_a-1*I_NAI_G4y_F3x;
  abcd[437] = 2.0E0*I_NAI_I3y3z_F3x_a-2*I_NAI_G3yz_F3x;
  abcd[438] = 2.0E0*I_NAI_I2y4z_F3x_a-3*I_NAI_G2y2z_F3x;
  abcd[439] = 2.0E0*I_NAI_Iy5z_F3x_a-4*I_NAI_Gy3z_F3x;
  abcd[440] = 2.0E0*I_NAI_I6z_F3x_a-5*I_NAI_G4z_F3x;
  abcd[441] = 2.0E0*I_NAI_I5xz_F2xy_a;
  abcd[442] = 2.0E0*I_NAI_I4xyz_F2xy_a;
  abcd[443] = 2.0E0*I_NAI_I4x2z_F2xy_a-1*I_NAI_G4x_F2xy;
  abcd[444] = 2.0E0*I_NAI_I3x2yz_F2xy_a;
  abcd[445] = 2.0E0*I_NAI_I3xy2z_F2xy_a-1*I_NAI_G3xy_F2xy;
  abcd[446] = 2.0E0*I_NAI_I3x3z_F2xy_a-2*I_NAI_G3xz_F2xy;
  abcd[447] = 2.0E0*I_NAI_I2x3yz_F2xy_a;
  abcd[448] = 2.0E0*I_NAI_I2x2y2z_F2xy_a-1*I_NAI_G2x2y_F2xy;
  abcd[449] = 2.0E0*I_NAI_I2xy3z_F2xy_a-2*I_NAI_G2xyz_F2xy;
  abcd[450] = 2.0E0*I_NAI_I2x4z_F2xy_a-3*I_NAI_G2x2z_F2xy;
  abcd[451] = 2.0E0*I_NAI_Ix4yz_F2xy_a;
  abcd[452] = 2.0E0*I_NAI_Ix3y2z_F2xy_a-1*I_NAI_Gx3y_F2xy;
  abcd[453] = 2.0E0*I_NAI_Ix2y3z_F2xy_a-2*I_NAI_Gx2yz_F2xy;
  abcd[454] = 2.0E0*I_NAI_Ixy4z_F2xy_a-3*I_NAI_Gxy2z_F2xy;
  abcd[455] = 2.0E0*I_NAI_Ix5z_F2xy_a-4*I_NAI_Gx3z_F2xy;
  abcd[456] = 2.0E0*I_NAI_I5yz_F2xy_a;
  abcd[457] = 2.0E0*I_NAI_I4y2z_F2xy_a-1*I_NAI_G4y_F2xy;
  abcd[458] = 2.0E0*I_NAI_I3y3z_F2xy_a-2*I_NAI_G3yz_F2xy;
  abcd[459] = 2.0E0*I_NAI_I2y4z_F2xy_a-3*I_NAI_G2y2z_F2xy;
  abcd[460] = 2.0E0*I_NAI_Iy5z_F2xy_a-4*I_NAI_Gy3z_F2xy;
  abcd[461] = 2.0E0*I_NAI_I6z_F2xy_a-5*I_NAI_G4z_F2xy;
  abcd[462] = 2.0E0*I_NAI_I5xz_F2xz_a;
  abcd[463] = 2.0E0*I_NAI_I4xyz_F2xz_a;
  abcd[464] = 2.0E0*I_NAI_I4x2z_F2xz_a-1*I_NAI_G4x_F2xz;
  abcd[465] = 2.0E0*I_NAI_I3x2yz_F2xz_a;
  abcd[466] = 2.0E0*I_NAI_I3xy2z_F2xz_a-1*I_NAI_G3xy_F2xz;
  abcd[467] = 2.0E0*I_NAI_I3x3z_F2xz_a-2*I_NAI_G3xz_F2xz;
  abcd[468] = 2.0E0*I_NAI_I2x3yz_F2xz_a;
  abcd[469] = 2.0E0*I_NAI_I2x2y2z_F2xz_a-1*I_NAI_G2x2y_F2xz;
  abcd[470] = 2.0E0*I_NAI_I2xy3z_F2xz_a-2*I_NAI_G2xyz_F2xz;
  abcd[471] = 2.0E0*I_NAI_I2x4z_F2xz_a-3*I_NAI_G2x2z_F2xz;
  abcd[472] = 2.0E0*I_NAI_Ix4yz_F2xz_a;
  abcd[473] = 2.0E0*I_NAI_Ix3y2z_F2xz_a-1*I_NAI_Gx3y_F2xz;
  abcd[474] = 2.0E0*I_NAI_Ix2y3z_F2xz_a-2*I_NAI_Gx2yz_F2xz;
  abcd[475] = 2.0E0*I_NAI_Ixy4z_F2xz_a-3*I_NAI_Gxy2z_F2xz;
  abcd[476] = 2.0E0*I_NAI_Ix5z_F2xz_a-4*I_NAI_Gx3z_F2xz;
  abcd[477] = 2.0E0*I_NAI_I5yz_F2xz_a;
  abcd[478] = 2.0E0*I_NAI_I4y2z_F2xz_a-1*I_NAI_G4y_F2xz;
  abcd[479] = 2.0E0*I_NAI_I3y3z_F2xz_a-2*I_NAI_G3yz_F2xz;
  abcd[480] = 2.0E0*I_NAI_I2y4z_F2xz_a-3*I_NAI_G2y2z_F2xz;
  abcd[481] = 2.0E0*I_NAI_Iy5z_F2xz_a-4*I_NAI_Gy3z_F2xz;
  abcd[482] = 2.0E0*I_NAI_I6z_F2xz_a-5*I_NAI_G4z_F2xz;
  abcd[483] = 2.0E0*I_NAI_I5xz_Fx2y_a;
  abcd[484] = 2.0E0*I_NAI_I4xyz_Fx2y_a;
  abcd[485] = 2.0E0*I_NAI_I4x2z_Fx2y_a-1*I_NAI_G4x_Fx2y;
  abcd[486] = 2.0E0*I_NAI_I3x2yz_Fx2y_a;
  abcd[487] = 2.0E0*I_NAI_I3xy2z_Fx2y_a-1*I_NAI_G3xy_Fx2y;
  abcd[488] = 2.0E0*I_NAI_I3x3z_Fx2y_a-2*I_NAI_G3xz_Fx2y;
  abcd[489] = 2.0E0*I_NAI_I2x3yz_Fx2y_a;
  abcd[490] = 2.0E0*I_NAI_I2x2y2z_Fx2y_a-1*I_NAI_G2x2y_Fx2y;
  abcd[491] = 2.0E0*I_NAI_I2xy3z_Fx2y_a-2*I_NAI_G2xyz_Fx2y;
  abcd[492] = 2.0E0*I_NAI_I2x4z_Fx2y_a-3*I_NAI_G2x2z_Fx2y;
  abcd[493] = 2.0E0*I_NAI_Ix4yz_Fx2y_a;
  abcd[494] = 2.0E0*I_NAI_Ix3y2z_Fx2y_a-1*I_NAI_Gx3y_Fx2y;
  abcd[495] = 2.0E0*I_NAI_Ix2y3z_Fx2y_a-2*I_NAI_Gx2yz_Fx2y;
  abcd[496] = 2.0E0*I_NAI_Ixy4z_Fx2y_a-3*I_NAI_Gxy2z_Fx2y;
  abcd[497] = 2.0E0*I_NAI_Ix5z_Fx2y_a-4*I_NAI_Gx3z_Fx2y;
  abcd[498] = 2.0E0*I_NAI_I5yz_Fx2y_a;
  abcd[499] = 2.0E0*I_NAI_I4y2z_Fx2y_a-1*I_NAI_G4y_Fx2y;
  abcd[500] = 2.0E0*I_NAI_I3y3z_Fx2y_a-2*I_NAI_G3yz_Fx2y;
  abcd[501] = 2.0E0*I_NAI_I2y4z_Fx2y_a-3*I_NAI_G2y2z_Fx2y;
  abcd[502] = 2.0E0*I_NAI_Iy5z_Fx2y_a-4*I_NAI_Gy3z_Fx2y;
  abcd[503] = 2.0E0*I_NAI_I6z_Fx2y_a-5*I_NAI_G4z_Fx2y;
  abcd[504] = 2.0E0*I_NAI_I5xz_Fxyz_a;
  abcd[505] = 2.0E0*I_NAI_I4xyz_Fxyz_a;
  abcd[506] = 2.0E0*I_NAI_I4x2z_Fxyz_a-1*I_NAI_G4x_Fxyz;
  abcd[507] = 2.0E0*I_NAI_I3x2yz_Fxyz_a;
  abcd[508] = 2.0E0*I_NAI_I3xy2z_Fxyz_a-1*I_NAI_G3xy_Fxyz;
  abcd[509] = 2.0E0*I_NAI_I3x3z_Fxyz_a-2*I_NAI_G3xz_Fxyz;
  abcd[510] = 2.0E0*I_NAI_I2x3yz_Fxyz_a;
  abcd[511] = 2.0E0*I_NAI_I2x2y2z_Fxyz_a-1*I_NAI_G2x2y_Fxyz;
  abcd[512] = 2.0E0*I_NAI_I2xy3z_Fxyz_a-2*I_NAI_G2xyz_Fxyz;
  abcd[513] = 2.0E0*I_NAI_I2x4z_Fxyz_a-3*I_NAI_G2x2z_Fxyz;
  abcd[514] = 2.0E0*I_NAI_Ix4yz_Fxyz_a;
  abcd[515] = 2.0E0*I_NAI_Ix3y2z_Fxyz_a-1*I_NAI_Gx3y_Fxyz;
  abcd[516] = 2.0E0*I_NAI_Ix2y3z_Fxyz_a-2*I_NAI_Gx2yz_Fxyz;
  abcd[517] = 2.0E0*I_NAI_Ixy4z_Fxyz_a-3*I_NAI_Gxy2z_Fxyz;
  abcd[518] = 2.0E0*I_NAI_Ix5z_Fxyz_a-4*I_NAI_Gx3z_Fxyz;
  abcd[519] = 2.0E0*I_NAI_I5yz_Fxyz_a;
  abcd[520] = 2.0E0*I_NAI_I4y2z_Fxyz_a-1*I_NAI_G4y_Fxyz;
  abcd[521] = 2.0E0*I_NAI_I3y3z_Fxyz_a-2*I_NAI_G3yz_Fxyz;
  abcd[522] = 2.0E0*I_NAI_I2y4z_Fxyz_a-3*I_NAI_G2y2z_Fxyz;
  abcd[523] = 2.0E0*I_NAI_Iy5z_Fxyz_a-4*I_NAI_Gy3z_Fxyz;
  abcd[524] = 2.0E0*I_NAI_I6z_Fxyz_a-5*I_NAI_G4z_Fxyz;
  abcd[525] = 2.0E0*I_NAI_I5xz_Fx2z_a;
  abcd[526] = 2.0E0*I_NAI_I4xyz_Fx2z_a;
  abcd[527] = 2.0E0*I_NAI_I4x2z_Fx2z_a-1*I_NAI_G4x_Fx2z;
  abcd[528] = 2.0E0*I_NAI_I3x2yz_Fx2z_a;
  abcd[529] = 2.0E0*I_NAI_I3xy2z_Fx2z_a-1*I_NAI_G3xy_Fx2z;
  abcd[530] = 2.0E0*I_NAI_I3x3z_Fx2z_a-2*I_NAI_G3xz_Fx2z;
  abcd[531] = 2.0E0*I_NAI_I2x3yz_Fx2z_a;
  abcd[532] = 2.0E0*I_NAI_I2x2y2z_Fx2z_a-1*I_NAI_G2x2y_Fx2z;
  abcd[533] = 2.0E0*I_NAI_I2xy3z_Fx2z_a-2*I_NAI_G2xyz_Fx2z;
  abcd[534] = 2.0E0*I_NAI_I2x4z_Fx2z_a-3*I_NAI_G2x2z_Fx2z;
  abcd[535] = 2.0E0*I_NAI_Ix4yz_Fx2z_a;
  abcd[536] = 2.0E0*I_NAI_Ix3y2z_Fx2z_a-1*I_NAI_Gx3y_Fx2z;
  abcd[537] = 2.0E0*I_NAI_Ix2y3z_Fx2z_a-2*I_NAI_Gx2yz_Fx2z;
  abcd[538] = 2.0E0*I_NAI_Ixy4z_Fx2z_a-3*I_NAI_Gxy2z_Fx2z;
  abcd[539] = 2.0E0*I_NAI_Ix5z_Fx2z_a-4*I_NAI_Gx3z_Fx2z;
  abcd[540] = 2.0E0*I_NAI_I5yz_Fx2z_a;
  abcd[541] = 2.0E0*I_NAI_I4y2z_Fx2z_a-1*I_NAI_G4y_Fx2z;
  abcd[542] = 2.0E0*I_NAI_I3y3z_Fx2z_a-2*I_NAI_G3yz_Fx2z;
  abcd[543] = 2.0E0*I_NAI_I2y4z_Fx2z_a-3*I_NAI_G2y2z_Fx2z;
  abcd[544] = 2.0E0*I_NAI_Iy5z_Fx2z_a-4*I_NAI_Gy3z_Fx2z;
  abcd[545] = 2.0E0*I_NAI_I6z_Fx2z_a-5*I_NAI_G4z_Fx2z;
  abcd[546] = 2.0E0*I_NAI_I5xz_F3y_a;
  abcd[547] = 2.0E0*I_NAI_I4xyz_F3y_a;
  abcd[548] = 2.0E0*I_NAI_I4x2z_F3y_a-1*I_NAI_G4x_F3y;
  abcd[549] = 2.0E0*I_NAI_I3x2yz_F3y_a;
  abcd[550] = 2.0E0*I_NAI_I3xy2z_F3y_a-1*I_NAI_G3xy_F3y;
  abcd[551] = 2.0E0*I_NAI_I3x3z_F3y_a-2*I_NAI_G3xz_F3y;
  abcd[552] = 2.0E0*I_NAI_I2x3yz_F3y_a;
  abcd[553] = 2.0E0*I_NAI_I2x2y2z_F3y_a-1*I_NAI_G2x2y_F3y;
  abcd[554] = 2.0E0*I_NAI_I2xy3z_F3y_a-2*I_NAI_G2xyz_F3y;
  abcd[555] = 2.0E0*I_NAI_I2x4z_F3y_a-3*I_NAI_G2x2z_F3y;
  abcd[556] = 2.0E0*I_NAI_Ix4yz_F3y_a;
  abcd[557] = 2.0E0*I_NAI_Ix3y2z_F3y_a-1*I_NAI_Gx3y_F3y;
  abcd[558] = 2.0E0*I_NAI_Ix2y3z_F3y_a-2*I_NAI_Gx2yz_F3y;
  abcd[559] = 2.0E0*I_NAI_Ixy4z_F3y_a-3*I_NAI_Gxy2z_F3y;
  abcd[560] = 2.0E0*I_NAI_Ix5z_F3y_a-4*I_NAI_Gx3z_F3y;
  abcd[561] = 2.0E0*I_NAI_I5yz_F3y_a;
  abcd[562] = 2.0E0*I_NAI_I4y2z_F3y_a-1*I_NAI_G4y_F3y;
  abcd[563] = 2.0E0*I_NAI_I3y3z_F3y_a-2*I_NAI_G3yz_F3y;
  abcd[564] = 2.0E0*I_NAI_I2y4z_F3y_a-3*I_NAI_G2y2z_F3y;
  abcd[565] = 2.0E0*I_NAI_Iy5z_F3y_a-4*I_NAI_Gy3z_F3y;
  abcd[566] = 2.0E0*I_NAI_I6z_F3y_a-5*I_NAI_G4z_F3y;
  abcd[567] = 2.0E0*I_NAI_I5xz_F2yz_a;
  abcd[568] = 2.0E0*I_NAI_I4xyz_F2yz_a;
  abcd[569] = 2.0E0*I_NAI_I4x2z_F2yz_a-1*I_NAI_G4x_F2yz;
  abcd[570] = 2.0E0*I_NAI_I3x2yz_F2yz_a;
  abcd[571] = 2.0E0*I_NAI_I3xy2z_F2yz_a-1*I_NAI_G3xy_F2yz;
  abcd[572] = 2.0E0*I_NAI_I3x3z_F2yz_a-2*I_NAI_G3xz_F2yz;
  abcd[573] = 2.0E0*I_NAI_I2x3yz_F2yz_a;
  abcd[574] = 2.0E0*I_NAI_I2x2y2z_F2yz_a-1*I_NAI_G2x2y_F2yz;
  abcd[575] = 2.0E0*I_NAI_I2xy3z_F2yz_a-2*I_NAI_G2xyz_F2yz;
  abcd[576] = 2.0E0*I_NAI_I2x4z_F2yz_a-3*I_NAI_G2x2z_F2yz;
  abcd[577] = 2.0E0*I_NAI_Ix4yz_F2yz_a;
  abcd[578] = 2.0E0*I_NAI_Ix3y2z_F2yz_a-1*I_NAI_Gx3y_F2yz;
  abcd[579] = 2.0E0*I_NAI_Ix2y3z_F2yz_a-2*I_NAI_Gx2yz_F2yz;
  abcd[580] = 2.0E0*I_NAI_Ixy4z_F2yz_a-3*I_NAI_Gxy2z_F2yz;
  abcd[581] = 2.0E0*I_NAI_Ix5z_F2yz_a-4*I_NAI_Gx3z_F2yz;
  abcd[582] = 2.0E0*I_NAI_I5yz_F2yz_a;
  abcd[583] = 2.0E0*I_NAI_I4y2z_F2yz_a-1*I_NAI_G4y_F2yz;
  abcd[584] = 2.0E0*I_NAI_I3y3z_F2yz_a-2*I_NAI_G3yz_F2yz;
  abcd[585] = 2.0E0*I_NAI_I2y4z_F2yz_a-3*I_NAI_G2y2z_F2yz;
  abcd[586] = 2.0E0*I_NAI_Iy5z_F2yz_a-4*I_NAI_Gy3z_F2yz;
  abcd[587] = 2.0E0*I_NAI_I6z_F2yz_a-5*I_NAI_G4z_F2yz;
  abcd[588] = 2.0E0*I_NAI_I5xz_Fy2z_a;
  abcd[589] = 2.0E0*I_NAI_I4xyz_Fy2z_a;
  abcd[590] = 2.0E0*I_NAI_I4x2z_Fy2z_a-1*I_NAI_G4x_Fy2z;
  abcd[591] = 2.0E0*I_NAI_I3x2yz_Fy2z_a;
  abcd[592] = 2.0E0*I_NAI_I3xy2z_Fy2z_a-1*I_NAI_G3xy_Fy2z;
  abcd[593] = 2.0E0*I_NAI_I3x3z_Fy2z_a-2*I_NAI_G3xz_Fy2z;
  abcd[594] = 2.0E0*I_NAI_I2x3yz_Fy2z_a;
  abcd[595] = 2.0E0*I_NAI_I2x2y2z_Fy2z_a-1*I_NAI_G2x2y_Fy2z;
  abcd[596] = 2.0E0*I_NAI_I2xy3z_Fy2z_a-2*I_NAI_G2xyz_Fy2z;
  abcd[597] = 2.0E0*I_NAI_I2x4z_Fy2z_a-3*I_NAI_G2x2z_Fy2z;
  abcd[598] = 2.0E0*I_NAI_Ix4yz_Fy2z_a;
  abcd[599] = 2.0E0*I_NAI_Ix3y2z_Fy2z_a-1*I_NAI_Gx3y_Fy2z;
  abcd[600] = 2.0E0*I_NAI_Ix2y3z_Fy2z_a-2*I_NAI_Gx2yz_Fy2z;
  abcd[601] = 2.0E0*I_NAI_Ixy4z_Fy2z_a-3*I_NAI_Gxy2z_Fy2z;
  abcd[602] = 2.0E0*I_NAI_Ix5z_Fy2z_a-4*I_NAI_Gx3z_Fy2z;
  abcd[603] = 2.0E0*I_NAI_I5yz_Fy2z_a;
  abcd[604] = 2.0E0*I_NAI_I4y2z_Fy2z_a-1*I_NAI_G4y_Fy2z;
  abcd[605] = 2.0E0*I_NAI_I3y3z_Fy2z_a-2*I_NAI_G3yz_Fy2z;
  abcd[606] = 2.0E0*I_NAI_I2y4z_Fy2z_a-3*I_NAI_G2y2z_Fy2z;
  abcd[607] = 2.0E0*I_NAI_Iy5z_Fy2z_a-4*I_NAI_Gy3z_Fy2z;
  abcd[608] = 2.0E0*I_NAI_I6z_Fy2z_a-5*I_NAI_G4z_Fy2z;
  abcd[609] = 2.0E0*I_NAI_I5xz_F3z_a;
  abcd[610] = 2.0E0*I_NAI_I4xyz_F3z_a;
  abcd[611] = 2.0E0*I_NAI_I4x2z_F3z_a-1*I_NAI_G4x_F3z;
  abcd[612] = 2.0E0*I_NAI_I3x2yz_F3z_a;
  abcd[613] = 2.0E0*I_NAI_I3xy2z_F3z_a-1*I_NAI_G3xy_F3z;
  abcd[614] = 2.0E0*I_NAI_I3x3z_F3z_a-2*I_NAI_G3xz_F3z;
  abcd[615] = 2.0E0*I_NAI_I2x3yz_F3z_a;
  abcd[616] = 2.0E0*I_NAI_I2x2y2z_F3z_a-1*I_NAI_G2x2y_F3z;
  abcd[617] = 2.0E0*I_NAI_I2xy3z_F3z_a-2*I_NAI_G2xyz_F3z;
  abcd[618] = 2.0E0*I_NAI_I2x4z_F3z_a-3*I_NAI_G2x2z_F3z;
  abcd[619] = 2.0E0*I_NAI_Ix4yz_F3z_a;
  abcd[620] = 2.0E0*I_NAI_Ix3y2z_F3z_a-1*I_NAI_Gx3y_F3z;
  abcd[621] = 2.0E0*I_NAI_Ix2y3z_F3z_a-2*I_NAI_Gx2yz_F3z;
  abcd[622] = 2.0E0*I_NAI_Ixy4z_F3z_a-3*I_NAI_Gxy2z_F3z;
  abcd[623] = 2.0E0*I_NAI_Ix5z_F3z_a-4*I_NAI_Gx3z_F3z;
  abcd[624] = 2.0E0*I_NAI_I5yz_F3z_a;
  abcd[625] = 2.0E0*I_NAI_I4y2z_F3z_a-1*I_NAI_G4y_F3z;
  abcd[626] = 2.0E0*I_NAI_I3y3z_F3z_a-2*I_NAI_G3yz_F3z;
  abcd[627] = 2.0E0*I_NAI_I2y4z_F3z_a-3*I_NAI_G2y2z_F3z;
  abcd[628] = 2.0E0*I_NAI_Iy5z_F3z_a-4*I_NAI_Gy3z_F3z;
  abcd[629] = 2.0E0*I_NAI_I6z_F3z_a-5*I_NAI_G4z_F3z;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_b
   * RHS shell quartet name: SQ_NAI_H_D
   ************************************************************/
  abcd[630] = 2.0E0*I_NAI_H5x_G4x_b-3*I_NAI_H5x_D2x;
  abcd[631] = 2.0E0*I_NAI_H4xy_G4x_b-3*I_NAI_H4xy_D2x;
  abcd[632] = 2.0E0*I_NAI_H4xz_G4x_b-3*I_NAI_H4xz_D2x;
  abcd[633] = 2.0E0*I_NAI_H3x2y_G4x_b-3*I_NAI_H3x2y_D2x;
  abcd[634] = 2.0E0*I_NAI_H3xyz_G4x_b-3*I_NAI_H3xyz_D2x;
  abcd[635] = 2.0E0*I_NAI_H3x2z_G4x_b-3*I_NAI_H3x2z_D2x;
  abcd[636] = 2.0E0*I_NAI_H2x3y_G4x_b-3*I_NAI_H2x3y_D2x;
  abcd[637] = 2.0E0*I_NAI_H2x2yz_G4x_b-3*I_NAI_H2x2yz_D2x;
  abcd[638] = 2.0E0*I_NAI_H2xy2z_G4x_b-3*I_NAI_H2xy2z_D2x;
  abcd[639] = 2.0E0*I_NAI_H2x3z_G4x_b-3*I_NAI_H2x3z_D2x;
  abcd[640] = 2.0E0*I_NAI_Hx4y_G4x_b-3*I_NAI_Hx4y_D2x;
  abcd[641] = 2.0E0*I_NAI_Hx3yz_G4x_b-3*I_NAI_Hx3yz_D2x;
  abcd[642] = 2.0E0*I_NAI_Hx2y2z_G4x_b-3*I_NAI_Hx2y2z_D2x;
  abcd[643] = 2.0E0*I_NAI_Hxy3z_G4x_b-3*I_NAI_Hxy3z_D2x;
  abcd[644] = 2.0E0*I_NAI_Hx4z_G4x_b-3*I_NAI_Hx4z_D2x;
  abcd[645] = 2.0E0*I_NAI_H5y_G4x_b-3*I_NAI_H5y_D2x;
  abcd[646] = 2.0E0*I_NAI_H4yz_G4x_b-3*I_NAI_H4yz_D2x;
  abcd[647] = 2.0E0*I_NAI_H3y2z_G4x_b-3*I_NAI_H3y2z_D2x;
  abcd[648] = 2.0E0*I_NAI_H2y3z_G4x_b-3*I_NAI_H2y3z_D2x;
  abcd[649] = 2.0E0*I_NAI_Hy4z_G4x_b-3*I_NAI_Hy4z_D2x;
  abcd[650] = 2.0E0*I_NAI_H5z_G4x_b-3*I_NAI_H5z_D2x;
  abcd[651] = 2.0E0*I_NAI_H5x_G3xy_b-2*I_NAI_H5x_Dxy;
  abcd[652] = 2.0E0*I_NAI_H4xy_G3xy_b-2*I_NAI_H4xy_Dxy;
  abcd[653] = 2.0E0*I_NAI_H4xz_G3xy_b-2*I_NAI_H4xz_Dxy;
  abcd[654] = 2.0E0*I_NAI_H3x2y_G3xy_b-2*I_NAI_H3x2y_Dxy;
  abcd[655] = 2.0E0*I_NAI_H3xyz_G3xy_b-2*I_NAI_H3xyz_Dxy;
  abcd[656] = 2.0E0*I_NAI_H3x2z_G3xy_b-2*I_NAI_H3x2z_Dxy;
  abcd[657] = 2.0E0*I_NAI_H2x3y_G3xy_b-2*I_NAI_H2x3y_Dxy;
  abcd[658] = 2.0E0*I_NAI_H2x2yz_G3xy_b-2*I_NAI_H2x2yz_Dxy;
  abcd[659] = 2.0E0*I_NAI_H2xy2z_G3xy_b-2*I_NAI_H2xy2z_Dxy;
  abcd[660] = 2.0E0*I_NAI_H2x3z_G3xy_b-2*I_NAI_H2x3z_Dxy;
  abcd[661] = 2.0E0*I_NAI_Hx4y_G3xy_b-2*I_NAI_Hx4y_Dxy;
  abcd[662] = 2.0E0*I_NAI_Hx3yz_G3xy_b-2*I_NAI_Hx3yz_Dxy;
  abcd[663] = 2.0E0*I_NAI_Hx2y2z_G3xy_b-2*I_NAI_Hx2y2z_Dxy;
  abcd[664] = 2.0E0*I_NAI_Hxy3z_G3xy_b-2*I_NAI_Hxy3z_Dxy;
  abcd[665] = 2.0E0*I_NAI_Hx4z_G3xy_b-2*I_NAI_Hx4z_Dxy;
  abcd[666] = 2.0E0*I_NAI_H5y_G3xy_b-2*I_NAI_H5y_Dxy;
  abcd[667] = 2.0E0*I_NAI_H4yz_G3xy_b-2*I_NAI_H4yz_Dxy;
  abcd[668] = 2.0E0*I_NAI_H3y2z_G3xy_b-2*I_NAI_H3y2z_Dxy;
  abcd[669] = 2.0E0*I_NAI_H2y3z_G3xy_b-2*I_NAI_H2y3z_Dxy;
  abcd[670] = 2.0E0*I_NAI_Hy4z_G3xy_b-2*I_NAI_Hy4z_Dxy;
  abcd[671] = 2.0E0*I_NAI_H5z_G3xy_b-2*I_NAI_H5z_Dxy;
  abcd[672] = 2.0E0*I_NAI_H5x_G3xz_b-2*I_NAI_H5x_Dxz;
  abcd[673] = 2.0E0*I_NAI_H4xy_G3xz_b-2*I_NAI_H4xy_Dxz;
  abcd[674] = 2.0E0*I_NAI_H4xz_G3xz_b-2*I_NAI_H4xz_Dxz;
  abcd[675] = 2.0E0*I_NAI_H3x2y_G3xz_b-2*I_NAI_H3x2y_Dxz;
  abcd[676] = 2.0E0*I_NAI_H3xyz_G3xz_b-2*I_NAI_H3xyz_Dxz;
  abcd[677] = 2.0E0*I_NAI_H3x2z_G3xz_b-2*I_NAI_H3x2z_Dxz;
  abcd[678] = 2.0E0*I_NAI_H2x3y_G3xz_b-2*I_NAI_H2x3y_Dxz;
  abcd[679] = 2.0E0*I_NAI_H2x2yz_G3xz_b-2*I_NAI_H2x2yz_Dxz;
  abcd[680] = 2.0E0*I_NAI_H2xy2z_G3xz_b-2*I_NAI_H2xy2z_Dxz;
  abcd[681] = 2.0E0*I_NAI_H2x3z_G3xz_b-2*I_NAI_H2x3z_Dxz;
  abcd[682] = 2.0E0*I_NAI_Hx4y_G3xz_b-2*I_NAI_Hx4y_Dxz;
  abcd[683] = 2.0E0*I_NAI_Hx3yz_G3xz_b-2*I_NAI_Hx3yz_Dxz;
  abcd[684] = 2.0E0*I_NAI_Hx2y2z_G3xz_b-2*I_NAI_Hx2y2z_Dxz;
  abcd[685] = 2.0E0*I_NAI_Hxy3z_G3xz_b-2*I_NAI_Hxy3z_Dxz;
  abcd[686] = 2.0E0*I_NAI_Hx4z_G3xz_b-2*I_NAI_Hx4z_Dxz;
  abcd[687] = 2.0E0*I_NAI_H5y_G3xz_b-2*I_NAI_H5y_Dxz;
  abcd[688] = 2.0E0*I_NAI_H4yz_G3xz_b-2*I_NAI_H4yz_Dxz;
  abcd[689] = 2.0E0*I_NAI_H3y2z_G3xz_b-2*I_NAI_H3y2z_Dxz;
  abcd[690] = 2.0E0*I_NAI_H2y3z_G3xz_b-2*I_NAI_H2y3z_Dxz;
  abcd[691] = 2.0E0*I_NAI_Hy4z_G3xz_b-2*I_NAI_Hy4z_Dxz;
  abcd[692] = 2.0E0*I_NAI_H5z_G3xz_b-2*I_NAI_H5z_Dxz;
  abcd[693] = 2.0E0*I_NAI_H5x_G2x2y_b-1*I_NAI_H5x_D2y;
  abcd[694] = 2.0E0*I_NAI_H4xy_G2x2y_b-1*I_NAI_H4xy_D2y;
  abcd[695] = 2.0E0*I_NAI_H4xz_G2x2y_b-1*I_NAI_H4xz_D2y;
  abcd[696] = 2.0E0*I_NAI_H3x2y_G2x2y_b-1*I_NAI_H3x2y_D2y;
  abcd[697] = 2.0E0*I_NAI_H3xyz_G2x2y_b-1*I_NAI_H3xyz_D2y;
  abcd[698] = 2.0E0*I_NAI_H3x2z_G2x2y_b-1*I_NAI_H3x2z_D2y;
  abcd[699] = 2.0E0*I_NAI_H2x3y_G2x2y_b-1*I_NAI_H2x3y_D2y;
  abcd[700] = 2.0E0*I_NAI_H2x2yz_G2x2y_b-1*I_NAI_H2x2yz_D2y;
  abcd[701] = 2.0E0*I_NAI_H2xy2z_G2x2y_b-1*I_NAI_H2xy2z_D2y;
  abcd[702] = 2.0E0*I_NAI_H2x3z_G2x2y_b-1*I_NAI_H2x3z_D2y;
  abcd[703] = 2.0E0*I_NAI_Hx4y_G2x2y_b-1*I_NAI_Hx4y_D2y;
  abcd[704] = 2.0E0*I_NAI_Hx3yz_G2x2y_b-1*I_NAI_Hx3yz_D2y;
  abcd[705] = 2.0E0*I_NAI_Hx2y2z_G2x2y_b-1*I_NAI_Hx2y2z_D2y;
  abcd[706] = 2.0E0*I_NAI_Hxy3z_G2x2y_b-1*I_NAI_Hxy3z_D2y;
  abcd[707] = 2.0E0*I_NAI_Hx4z_G2x2y_b-1*I_NAI_Hx4z_D2y;
  abcd[708] = 2.0E0*I_NAI_H5y_G2x2y_b-1*I_NAI_H5y_D2y;
  abcd[709] = 2.0E0*I_NAI_H4yz_G2x2y_b-1*I_NAI_H4yz_D2y;
  abcd[710] = 2.0E0*I_NAI_H3y2z_G2x2y_b-1*I_NAI_H3y2z_D2y;
  abcd[711] = 2.0E0*I_NAI_H2y3z_G2x2y_b-1*I_NAI_H2y3z_D2y;
  abcd[712] = 2.0E0*I_NAI_Hy4z_G2x2y_b-1*I_NAI_Hy4z_D2y;
  abcd[713] = 2.0E0*I_NAI_H5z_G2x2y_b-1*I_NAI_H5z_D2y;
  abcd[714] = 2.0E0*I_NAI_H5x_G2xyz_b-1*I_NAI_H5x_Dyz;
  abcd[715] = 2.0E0*I_NAI_H4xy_G2xyz_b-1*I_NAI_H4xy_Dyz;
  abcd[716] = 2.0E0*I_NAI_H4xz_G2xyz_b-1*I_NAI_H4xz_Dyz;
  abcd[717] = 2.0E0*I_NAI_H3x2y_G2xyz_b-1*I_NAI_H3x2y_Dyz;
  abcd[718] = 2.0E0*I_NAI_H3xyz_G2xyz_b-1*I_NAI_H3xyz_Dyz;
  abcd[719] = 2.0E0*I_NAI_H3x2z_G2xyz_b-1*I_NAI_H3x2z_Dyz;
  abcd[720] = 2.0E0*I_NAI_H2x3y_G2xyz_b-1*I_NAI_H2x3y_Dyz;
  abcd[721] = 2.0E0*I_NAI_H2x2yz_G2xyz_b-1*I_NAI_H2x2yz_Dyz;
  abcd[722] = 2.0E0*I_NAI_H2xy2z_G2xyz_b-1*I_NAI_H2xy2z_Dyz;
  abcd[723] = 2.0E0*I_NAI_H2x3z_G2xyz_b-1*I_NAI_H2x3z_Dyz;
  abcd[724] = 2.0E0*I_NAI_Hx4y_G2xyz_b-1*I_NAI_Hx4y_Dyz;
  abcd[725] = 2.0E0*I_NAI_Hx3yz_G2xyz_b-1*I_NAI_Hx3yz_Dyz;
  abcd[726] = 2.0E0*I_NAI_Hx2y2z_G2xyz_b-1*I_NAI_Hx2y2z_Dyz;
  abcd[727] = 2.0E0*I_NAI_Hxy3z_G2xyz_b-1*I_NAI_Hxy3z_Dyz;
  abcd[728] = 2.0E0*I_NAI_Hx4z_G2xyz_b-1*I_NAI_Hx4z_Dyz;
  abcd[729] = 2.0E0*I_NAI_H5y_G2xyz_b-1*I_NAI_H5y_Dyz;
  abcd[730] = 2.0E0*I_NAI_H4yz_G2xyz_b-1*I_NAI_H4yz_Dyz;
  abcd[731] = 2.0E0*I_NAI_H3y2z_G2xyz_b-1*I_NAI_H3y2z_Dyz;
  abcd[732] = 2.0E0*I_NAI_H2y3z_G2xyz_b-1*I_NAI_H2y3z_Dyz;
  abcd[733] = 2.0E0*I_NAI_Hy4z_G2xyz_b-1*I_NAI_Hy4z_Dyz;
  abcd[734] = 2.0E0*I_NAI_H5z_G2xyz_b-1*I_NAI_H5z_Dyz;
  abcd[735] = 2.0E0*I_NAI_H5x_G2x2z_b-1*I_NAI_H5x_D2z;
  abcd[736] = 2.0E0*I_NAI_H4xy_G2x2z_b-1*I_NAI_H4xy_D2z;
  abcd[737] = 2.0E0*I_NAI_H4xz_G2x2z_b-1*I_NAI_H4xz_D2z;
  abcd[738] = 2.0E0*I_NAI_H3x2y_G2x2z_b-1*I_NAI_H3x2y_D2z;
  abcd[739] = 2.0E0*I_NAI_H3xyz_G2x2z_b-1*I_NAI_H3xyz_D2z;
  abcd[740] = 2.0E0*I_NAI_H3x2z_G2x2z_b-1*I_NAI_H3x2z_D2z;
  abcd[741] = 2.0E0*I_NAI_H2x3y_G2x2z_b-1*I_NAI_H2x3y_D2z;
  abcd[742] = 2.0E0*I_NAI_H2x2yz_G2x2z_b-1*I_NAI_H2x2yz_D2z;
  abcd[743] = 2.0E0*I_NAI_H2xy2z_G2x2z_b-1*I_NAI_H2xy2z_D2z;
  abcd[744] = 2.0E0*I_NAI_H2x3z_G2x2z_b-1*I_NAI_H2x3z_D2z;
  abcd[745] = 2.0E0*I_NAI_Hx4y_G2x2z_b-1*I_NAI_Hx4y_D2z;
  abcd[746] = 2.0E0*I_NAI_Hx3yz_G2x2z_b-1*I_NAI_Hx3yz_D2z;
  abcd[747] = 2.0E0*I_NAI_Hx2y2z_G2x2z_b-1*I_NAI_Hx2y2z_D2z;
  abcd[748] = 2.0E0*I_NAI_Hxy3z_G2x2z_b-1*I_NAI_Hxy3z_D2z;
  abcd[749] = 2.0E0*I_NAI_Hx4z_G2x2z_b-1*I_NAI_Hx4z_D2z;
  abcd[750] = 2.0E0*I_NAI_H5y_G2x2z_b-1*I_NAI_H5y_D2z;
  abcd[751] = 2.0E0*I_NAI_H4yz_G2x2z_b-1*I_NAI_H4yz_D2z;
  abcd[752] = 2.0E0*I_NAI_H3y2z_G2x2z_b-1*I_NAI_H3y2z_D2z;
  abcd[753] = 2.0E0*I_NAI_H2y3z_G2x2z_b-1*I_NAI_H2y3z_D2z;
  abcd[754] = 2.0E0*I_NAI_Hy4z_G2x2z_b-1*I_NAI_Hy4z_D2z;
  abcd[755] = 2.0E0*I_NAI_H5z_G2x2z_b-1*I_NAI_H5z_D2z;
  abcd[756] = 2.0E0*I_NAI_H5x_Gx3y_b;
  abcd[757] = 2.0E0*I_NAI_H4xy_Gx3y_b;
  abcd[758] = 2.0E0*I_NAI_H4xz_Gx3y_b;
  abcd[759] = 2.0E0*I_NAI_H3x2y_Gx3y_b;
  abcd[760] = 2.0E0*I_NAI_H3xyz_Gx3y_b;
  abcd[761] = 2.0E0*I_NAI_H3x2z_Gx3y_b;
  abcd[762] = 2.0E0*I_NAI_H2x3y_Gx3y_b;
  abcd[763] = 2.0E0*I_NAI_H2x2yz_Gx3y_b;
  abcd[764] = 2.0E0*I_NAI_H2xy2z_Gx3y_b;
  abcd[765] = 2.0E0*I_NAI_H2x3z_Gx3y_b;
  abcd[766] = 2.0E0*I_NAI_Hx4y_Gx3y_b;
  abcd[767] = 2.0E0*I_NAI_Hx3yz_Gx3y_b;
  abcd[768] = 2.0E0*I_NAI_Hx2y2z_Gx3y_b;
  abcd[769] = 2.0E0*I_NAI_Hxy3z_Gx3y_b;
  abcd[770] = 2.0E0*I_NAI_Hx4z_Gx3y_b;
  abcd[771] = 2.0E0*I_NAI_H5y_Gx3y_b;
  abcd[772] = 2.0E0*I_NAI_H4yz_Gx3y_b;
  abcd[773] = 2.0E0*I_NAI_H3y2z_Gx3y_b;
  abcd[774] = 2.0E0*I_NAI_H2y3z_Gx3y_b;
  abcd[775] = 2.0E0*I_NAI_Hy4z_Gx3y_b;
  abcd[776] = 2.0E0*I_NAI_H5z_Gx3y_b;
  abcd[777] = 2.0E0*I_NAI_H5x_Gx2yz_b;
  abcd[778] = 2.0E0*I_NAI_H4xy_Gx2yz_b;
  abcd[779] = 2.0E0*I_NAI_H4xz_Gx2yz_b;
  abcd[780] = 2.0E0*I_NAI_H3x2y_Gx2yz_b;
  abcd[781] = 2.0E0*I_NAI_H3xyz_Gx2yz_b;
  abcd[782] = 2.0E0*I_NAI_H3x2z_Gx2yz_b;
  abcd[783] = 2.0E0*I_NAI_H2x3y_Gx2yz_b;
  abcd[784] = 2.0E0*I_NAI_H2x2yz_Gx2yz_b;
  abcd[785] = 2.0E0*I_NAI_H2xy2z_Gx2yz_b;
  abcd[786] = 2.0E0*I_NAI_H2x3z_Gx2yz_b;
  abcd[787] = 2.0E0*I_NAI_Hx4y_Gx2yz_b;
  abcd[788] = 2.0E0*I_NAI_Hx3yz_Gx2yz_b;
  abcd[789] = 2.0E0*I_NAI_Hx2y2z_Gx2yz_b;
  abcd[790] = 2.0E0*I_NAI_Hxy3z_Gx2yz_b;
  abcd[791] = 2.0E0*I_NAI_Hx4z_Gx2yz_b;
  abcd[792] = 2.0E0*I_NAI_H5y_Gx2yz_b;
  abcd[793] = 2.0E0*I_NAI_H4yz_Gx2yz_b;
  abcd[794] = 2.0E0*I_NAI_H3y2z_Gx2yz_b;
  abcd[795] = 2.0E0*I_NAI_H2y3z_Gx2yz_b;
  abcd[796] = 2.0E0*I_NAI_Hy4z_Gx2yz_b;
  abcd[797] = 2.0E0*I_NAI_H5z_Gx2yz_b;
  abcd[798] = 2.0E0*I_NAI_H5x_Gxy2z_b;
  abcd[799] = 2.0E0*I_NAI_H4xy_Gxy2z_b;
  abcd[800] = 2.0E0*I_NAI_H4xz_Gxy2z_b;
  abcd[801] = 2.0E0*I_NAI_H3x2y_Gxy2z_b;
  abcd[802] = 2.0E0*I_NAI_H3xyz_Gxy2z_b;
  abcd[803] = 2.0E0*I_NAI_H3x2z_Gxy2z_b;
  abcd[804] = 2.0E0*I_NAI_H2x3y_Gxy2z_b;
  abcd[805] = 2.0E0*I_NAI_H2x2yz_Gxy2z_b;
  abcd[806] = 2.0E0*I_NAI_H2xy2z_Gxy2z_b;
  abcd[807] = 2.0E0*I_NAI_H2x3z_Gxy2z_b;
  abcd[808] = 2.0E0*I_NAI_Hx4y_Gxy2z_b;
  abcd[809] = 2.0E0*I_NAI_Hx3yz_Gxy2z_b;
  abcd[810] = 2.0E0*I_NAI_Hx2y2z_Gxy2z_b;
  abcd[811] = 2.0E0*I_NAI_Hxy3z_Gxy2z_b;
  abcd[812] = 2.0E0*I_NAI_Hx4z_Gxy2z_b;
  abcd[813] = 2.0E0*I_NAI_H5y_Gxy2z_b;
  abcd[814] = 2.0E0*I_NAI_H4yz_Gxy2z_b;
  abcd[815] = 2.0E0*I_NAI_H3y2z_Gxy2z_b;
  abcd[816] = 2.0E0*I_NAI_H2y3z_Gxy2z_b;
  abcd[817] = 2.0E0*I_NAI_Hy4z_Gxy2z_b;
  abcd[818] = 2.0E0*I_NAI_H5z_Gxy2z_b;
  abcd[819] = 2.0E0*I_NAI_H5x_Gx3z_b;
  abcd[820] = 2.0E0*I_NAI_H4xy_Gx3z_b;
  abcd[821] = 2.0E0*I_NAI_H4xz_Gx3z_b;
  abcd[822] = 2.0E0*I_NAI_H3x2y_Gx3z_b;
  abcd[823] = 2.0E0*I_NAI_H3xyz_Gx3z_b;
  abcd[824] = 2.0E0*I_NAI_H3x2z_Gx3z_b;
  abcd[825] = 2.0E0*I_NAI_H2x3y_Gx3z_b;
  abcd[826] = 2.0E0*I_NAI_H2x2yz_Gx3z_b;
  abcd[827] = 2.0E0*I_NAI_H2xy2z_Gx3z_b;
  abcd[828] = 2.0E0*I_NAI_H2x3z_Gx3z_b;
  abcd[829] = 2.0E0*I_NAI_Hx4y_Gx3z_b;
  abcd[830] = 2.0E0*I_NAI_Hx3yz_Gx3z_b;
  abcd[831] = 2.0E0*I_NAI_Hx2y2z_Gx3z_b;
  abcd[832] = 2.0E0*I_NAI_Hxy3z_Gx3z_b;
  abcd[833] = 2.0E0*I_NAI_Hx4z_Gx3z_b;
  abcd[834] = 2.0E0*I_NAI_H5y_Gx3z_b;
  abcd[835] = 2.0E0*I_NAI_H4yz_Gx3z_b;
  abcd[836] = 2.0E0*I_NAI_H3y2z_Gx3z_b;
  abcd[837] = 2.0E0*I_NAI_H2y3z_Gx3z_b;
  abcd[838] = 2.0E0*I_NAI_Hy4z_Gx3z_b;
  abcd[839] = 2.0E0*I_NAI_H5z_Gx3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_b
   * RHS shell quartet name: SQ_NAI_H_D
   ************************************************************/
  abcd[840] = 2.0E0*I_NAI_H5x_G3xy_b;
  abcd[841] = 2.0E0*I_NAI_H4xy_G3xy_b;
  abcd[842] = 2.0E0*I_NAI_H4xz_G3xy_b;
  abcd[843] = 2.0E0*I_NAI_H3x2y_G3xy_b;
  abcd[844] = 2.0E0*I_NAI_H3xyz_G3xy_b;
  abcd[845] = 2.0E0*I_NAI_H3x2z_G3xy_b;
  abcd[846] = 2.0E0*I_NAI_H2x3y_G3xy_b;
  abcd[847] = 2.0E0*I_NAI_H2x2yz_G3xy_b;
  abcd[848] = 2.0E0*I_NAI_H2xy2z_G3xy_b;
  abcd[849] = 2.0E0*I_NAI_H2x3z_G3xy_b;
  abcd[850] = 2.0E0*I_NAI_Hx4y_G3xy_b;
  abcd[851] = 2.0E0*I_NAI_Hx3yz_G3xy_b;
  abcd[852] = 2.0E0*I_NAI_Hx2y2z_G3xy_b;
  abcd[853] = 2.0E0*I_NAI_Hxy3z_G3xy_b;
  abcd[854] = 2.0E0*I_NAI_Hx4z_G3xy_b;
  abcd[855] = 2.0E0*I_NAI_H5y_G3xy_b;
  abcd[856] = 2.0E0*I_NAI_H4yz_G3xy_b;
  abcd[857] = 2.0E0*I_NAI_H3y2z_G3xy_b;
  abcd[858] = 2.0E0*I_NAI_H2y3z_G3xy_b;
  abcd[859] = 2.0E0*I_NAI_Hy4z_G3xy_b;
  abcd[860] = 2.0E0*I_NAI_H5z_G3xy_b;
  abcd[861] = 2.0E0*I_NAI_H5x_G2x2y_b-1*I_NAI_H5x_D2x;
  abcd[862] = 2.0E0*I_NAI_H4xy_G2x2y_b-1*I_NAI_H4xy_D2x;
  abcd[863] = 2.0E0*I_NAI_H4xz_G2x2y_b-1*I_NAI_H4xz_D2x;
  abcd[864] = 2.0E0*I_NAI_H3x2y_G2x2y_b-1*I_NAI_H3x2y_D2x;
  abcd[865] = 2.0E0*I_NAI_H3xyz_G2x2y_b-1*I_NAI_H3xyz_D2x;
  abcd[866] = 2.0E0*I_NAI_H3x2z_G2x2y_b-1*I_NAI_H3x2z_D2x;
  abcd[867] = 2.0E0*I_NAI_H2x3y_G2x2y_b-1*I_NAI_H2x3y_D2x;
  abcd[868] = 2.0E0*I_NAI_H2x2yz_G2x2y_b-1*I_NAI_H2x2yz_D2x;
  abcd[869] = 2.0E0*I_NAI_H2xy2z_G2x2y_b-1*I_NAI_H2xy2z_D2x;
  abcd[870] = 2.0E0*I_NAI_H2x3z_G2x2y_b-1*I_NAI_H2x3z_D2x;
  abcd[871] = 2.0E0*I_NAI_Hx4y_G2x2y_b-1*I_NAI_Hx4y_D2x;
  abcd[872] = 2.0E0*I_NAI_Hx3yz_G2x2y_b-1*I_NAI_Hx3yz_D2x;
  abcd[873] = 2.0E0*I_NAI_Hx2y2z_G2x2y_b-1*I_NAI_Hx2y2z_D2x;
  abcd[874] = 2.0E0*I_NAI_Hxy3z_G2x2y_b-1*I_NAI_Hxy3z_D2x;
  abcd[875] = 2.0E0*I_NAI_Hx4z_G2x2y_b-1*I_NAI_Hx4z_D2x;
  abcd[876] = 2.0E0*I_NAI_H5y_G2x2y_b-1*I_NAI_H5y_D2x;
  abcd[877] = 2.0E0*I_NAI_H4yz_G2x2y_b-1*I_NAI_H4yz_D2x;
  abcd[878] = 2.0E0*I_NAI_H3y2z_G2x2y_b-1*I_NAI_H3y2z_D2x;
  abcd[879] = 2.0E0*I_NAI_H2y3z_G2x2y_b-1*I_NAI_H2y3z_D2x;
  abcd[880] = 2.0E0*I_NAI_Hy4z_G2x2y_b-1*I_NAI_Hy4z_D2x;
  abcd[881] = 2.0E0*I_NAI_H5z_G2x2y_b-1*I_NAI_H5z_D2x;
  abcd[882] = 2.0E0*I_NAI_H5x_G2xyz_b;
  abcd[883] = 2.0E0*I_NAI_H4xy_G2xyz_b;
  abcd[884] = 2.0E0*I_NAI_H4xz_G2xyz_b;
  abcd[885] = 2.0E0*I_NAI_H3x2y_G2xyz_b;
  abcd[886] = 2.0E0*I_NAI_H3xyz_G2xyz_b;
  abcd[887] = 2.0E0*I_NAI_H3x2z_G2xyz_b;
  abcd[888] = 2.0E0*I_NAI_H2x3y_G2xyz_b;
  abcd[889] = 2.0E0*I_NAI_H2x2yz_G2xyz_b;
  abcd[890] = 2.0E0*I_NAI_H2xy2z_G2xyz_b;
  abcd[891] = 2.0E0*I_NAI_H2x3z_G2xyz_b;
  abcd[892] = 2.0E0*I_NAI_Hx4y_G2xyz_b;
  abcd[893] = 2.0E0*I_NAI_Hx3yz_G2xyz_b;
  abcd[894] = 2.0E0*I_NAI_Hx2y2z_G2xyz_b;
  abcd[895] = 2.0E0*I_NAI_Hxy3z_G2xyz_b;
  abcd[896] = 2.0E0*I_NAI_Hx4z_G2xyz_b;
  abcd[897] = 2.0E0*I_NAI_H5y_G2xyz_b;
  abcd[898] = 2.0E0*I_NAI_H4yz_G2xyz_b;
  abcd[899] = 2.0E0*I_NAI_H3y2z_G2xyz_b;
  abcd[900] = 2.0E0*I_NAI_H2y3z_G2xyz_b;
  abcd[901] = 2.0E0*I_NAI_Hy4z_G2xyz_b;
  abcd[902] = 2.0E0*I_NAI_H5z_G2xyz_b;
  abcd[903] = 2.0E0*I_NAI_H5x_Gx3y_b-2*I_NAI_H5x_Dxy;
  abcd[904] = 2.0E0*I_NAI_H4xy_Gx3y_b-2*I_NAI_H4xy_Dxy;
  abcd[905] = 2.0E0*I_NAI_H4xz_Gx3y_b-2*I_NAI_H4xz_Dxy;
  abcd[906] = 2.0E0*I_NAI_H3x2y_Gx3y_b-2*I_NAI_H3x2y_Dxy;
  abcd[907] = 2.0E0*I_NAI_H3xyz_Gx3y_b-2*I_NAI_H3xyz_Dxy;
  abcd[908] = 2.0E0*I_NAI_H3x2z_Gx3y_b-2*I_NAI_H3x2z_Dxy;
  abcd[909] = 2.0E0*I_NAI_H2x3y_Gx3y_b-2*I_NAI_H2x3y_Dxy;
  abcd[910] = 2.0E0*I_NAI_H2x2yz_Gx3y_b-2*I_NAI_H2x2yz_Dxy;
  abcd[911] = 2.0E0*I_NAI_H2xy2z_Gx3y_b-2*I_NAI_H2xy2z_Dxy;
  abcd[912] = 2.0E0*I_NAI_H2x3z_Gx3y_b-2*I_NAI_H2x3z_Dxy;
  abcd[913] = 2.0E0*I_NAI_Hx4y_Gx3y_b-2*I_NAI_Hx4y_Dxy;
  abcd[914] = 2.0E0*I_NAI_Hx3yz_Gx3y_b-2*I_NAI_Hx3yz_Dxy;
  abcd[915] = 2.0E0*I_NAI_Hx2y2z_Gx3y_b-2*I_NAI_Hx2y2z_Dxy;
  abcd[916] = 2.0E0*I_NAI_Hxy3z_Gx3y_b-2*I_NAI_Hxy3z_Dxy;
  abcd[917] = 2.0E0*I_NAI_Hx4z_Gx3y_b-2*I_NAI_Hx4z_Dxy;
  abcd[918] = 2.0E0*I_NAI_H5y_Gx3y_b-2*I_NAI_H5y_Dxy;
  abcd[919] = 2.0E0*I_NAI_H4yz_Gx3y_b-2*I_NAI_H4yz_Dxy;
  abcd[920] = 2.0E0*I_NAI_H3y2z_Gx3y_b-2*I_NAI_H3y2z_Dxy;
  abcd[921] = 2.0E0*I_NAI_H2y3z_Gx3y_b-2*I_NAI_H2y3z_Dxy;
  abcd[922] = 2.0E0*I_NAI_Hy4z_Gx3y_b-2*I_NAI_Hy4z_Dxy;
  abcd[923] = 2.0E0*I_NAI_H5z_Gx3y_b-2*I_NAI_H5z_Dxy;
  abcd[924] = 2.0E0*I_NAI_H5x_Gx2yz_b-1*I_NAI_H5x_Dxz;
  abcd[925] = 2.0E0*I_NAI_H4xy_Gx2yz_b-1*I_NAI_H4xy_Dxz;
  abcd[926] = 2.0E0*I_NAI_H4xz_Gx2yz_b-1*I_NAI_H4xz_Dxz;
  abcd[927] = 2.0E0*I_NAI_H3x2y_Gx2yz_b-1*I_NAI_H3x2y_Dxz;
  abcd[928] = 2.0E0*I_NAI_H3xyz_Gx2yz_b-1*I_NAI_H3xyz_Dxz;
  abcd[929] = 2.0E0*I_NAI_H3x2z_Gx2yz_b-1*I_NAI_H3x2z_Dxz;
  abcd[930] = 2.0E0*I_NAI_H2x3y_Gx2yz_b-1*I_NAI_H2x3y_Dxz;
  abcd[931] = 2.0E0*I_NAI_H2x2yz_Gx2yz_b-1*I_NAI_H2x2yz_Dxz;
  abcd[932] = 2.0E0*I_NAI_H2xy2z_Gx2yz_b-1*I_NAI_H2xy2z_Dxz;
  abcd[933] = 2.0E0*I_NAI_H2x3z_Gx2yz_b-1*I_NAI_H2x3z_Dxz;
  abcd[934] = 2.0E0*I_NAI_Hx4y_Gx2yz_b-1*I_NAI_Hx4y_Dxz;
  abcd[935] = 2.0E0*I_NAI_Hx3yz_Gx2yz_b-1*I_NAI_Hx3yz_Dxz;
  abcd[936] = 2.0E0*I_NAI_Hx2y2z_Gx2yz_b-1*I_NAI_Hx2y2z_Dxz;
  abcd[937] = 2.0E0*I_NAI_Hxy3z_Gx2yz_b-1*I_NAI_Hxy3z_Dxz;
  abcd[938] = 2.0E0*I_NAI_Hx4z_Gx2yz_b-1*I_NAI_Hx4z_Dxz;
  abcd[939] = 2.0E0*I_NAI_H5y_Gx2yz_b-1*I_NAI_H5y_Dxz;
  abcd[940] = 2.0E0*I_NAI_H4yz_Gx2yz_b-1*I_NAI_H4yz_Dxz;
  abcd[941] = 2.0E0*I_NAI_H3y2z_Gx2yz_b-1*I_NAI_H3y2z_Dxz;
  abcd[942] = 2.0E0*I_NAI_H2y3z_Gx2yz_b-1*I_NAI_H2y3z_Dxz;
  abcd[943] = 2.0E0*I_NAI_Hy4z_Gx2yz_b-1*I_NAI_Hy4z_Dxz;
  abcd[944] = 2.0E0*I_NAI_H5z_Gx2yz_b-1*I_NAI_H5z_Dxz;
  abcd[945] = 2.0E0*I_NAI_H5x_Gxy2z_b;
  abcd[946] = 2.0E0*I_NAI_H4xy_Gxy2z_b;
  abcd[947] = 2.0E0*I_NAI_H4xz_Gxy2z_b;
  abcd[948] = 2.0E0*I_NAI_H3x2y_Gxy2z_b;
  abcd[949] = 2.0E0*I_NAI_H3xyz_Gxy2z_b;
  abcd[950] = 2.0E0*I_NAI_H3x2z_Gxy2z_b;
  abcd[951] = 2.0E0*I_NAI_H2x3y_Gxy2z_b;
  abcd[952] = 2.0E0*I_NAI_H2x2yz_Gxy2z_b;
  abcd[953] = 2.0E0*I_NAI_H2xy2z_Gxy2z_b;
  abcd[954] = 2.0E0*I_NAI_H2x3z_Gxy2z_b;
  abcd[955] = 2.0E0*I_NAI_Hx4y_Gxy2z_b;
  abcd[956] = 2.0E0*I_NAI_Hx3yz_Gxy2z_b;
  abcd[957] = 2.0E0*I_NAI_Hx2y2z_Gxy2z_b;
  abcd[958] = 2.0E0*I_NAI_Hxy3z_Gxy2z_b;
  abcd[959] = 2.0E0*I_NAI_Hx4z_Gxy2z_b;
  abcd[960] = 2.0E0*I_NAI_H5y_Gxy2z_b;
  abcd[961] = 2.0E0*I_NAI_H4yz_Gxy2z_b;
  abcd[962] = 2.0E0*I_NAI_H3y2z_Gxy2z_b;
  abcd[963] = 2.0E0*I_NAI_H2y3z_Gxy2z_b;
  abcd[964] = 2.0E0*I_NAI_Hy4z_Gxy2z_b;
  abcd[965] = 2.0E0*I_NAI_H5z_Gxy2z_b;
  abcd[966] = 2.0E0*I_NAI_H5x_G4y_b-3*I_NAI_H5x_D2y;
  abcd[967] = 2.0E0*I_NAI_H4xy_G4y_b-3*I_NAI_H4xy_D2y;
  abcd[968] = 2.0E0*I_NAI_H4xz_G4y_b-3*I_NAI_H4xz_D2y;
  abcd[969] = 2.0E0*I_NAI_H3x2y_G4y_b-3*I_NAI_H3x2y_D2y;
  abcd[970] = 2.0E0*I_NAI_H3xyz_G4y_b-3*I_NAI_H3xyz_D2y;
  abcd[971] = 2.0E0*I_NAI_H3x2z_G4y_b-3*I_NAI_H3x2z_D2y;
  abcd[972] = 2.0E0*I_NAI_H2x3y_G4y_b-3*I_NAI_H2x3y_D2y;
  abcd[973] = 2.0E0*I_NAI_H2x2yz_G4y_b-3*I_NAI_H2x2yz_D2y;
  abcd[974] = 2.0E0*I_NAI_H2xy2z_G4y_b-3*I_NAI_H2xy2z_D2y;
  abcd[975] = 2.0E0*I_NAI_H2x3z_G4y_b-3*I_NAI_H2x3z_D2y;
  abcd[976] = 2.0E0*I_NAI_Hx4y_G4y_b-3*I_NAI_Hx4y_D2y;
  abcd[977] = 2.0E0*I_NAI_Hx3yz_G4y_b-3*I_NAI_Hx3yz_D2y;
  abcd[978] = 2.0E0*I_NAI_Hx2y2z_G4y_b-3*I_NAI_Hx2y2z_D2y;
  abcd[979] = 2.0E0*I_NAI_Hxy3z_G4y_b-3*I_NAI_Hxy3z_D2y;
  abcd[980] = 2.0E0*I_NAI_Hx4z_G4y_b-3*I_NAI_Hx4z_D2y;
  abcd[981] = 2.0E0*I_NAI_H5y_G4y_b-3*I_NAI_H5y_D2y;
  abcd[982] = 2.0E0*I_NAI_H4yz_G4y_b-3*I_NAI_H4yz_D2y;
  abcd[983] = 2.0E0*I_NAI_H3y2z_G4y_b-3*I_NAI_H3y2z_D2y;
  abcd[984] = 2.0E0*I_NAI_H2y3z_G4y_b-3*I_NAI_H2y3z_D2y;
  abcd[985] = 2.0E0*I_NAI_Hy4z_G4y_b-3*I_NAI_Hy4z_D2y;
  abcd[986] = 2.0E0*I_NAI_H5z_G4y_b-3*I_NAI_H5z_D2y;
  abcd[987] = 2.0E0*I_NAI_H5x_G3yz_b-2*I_NAI_H5x_Dyz;
  abcd[988] = 2.0E0*I_NAI_H4xy_G3yz_b-2*I_NAI_H4xy_Dyz;
  abcd[989] = 2.0E0*I_NAI_H4xz_G3yz_b-2*I_NAI_H4xz_Dyz;
  abcd[990] = 2.0E0*I_NAI_H3x2y_G3yz_b-2*I_NAI_H3x2y_Dyz;
  abcd[991] = 2.0E0*I_NAI_H3xyz_G3yz_b-2*I_NAI_H3xyz_Dyz;
  abcd[992] = 2.0E0*I_NAI_H3x2z_G3yz_b-2*I_NAI_H3x2z_Dyz;
  abcd[993] = 2.0E0*I_NAI_H2x3y_G3yz_b-2*I_NAI_H2x3y_Dyz;
  abcd[994] = 2.0E0*I_NAI_H2x2yz_G3yz_b-2*I_NAI_H2x2yz_Dyz;
  abcd[995] = 2.0E0*I_NAI_H2xy2z_G3yz_b-2*I_NAI_H2xy2z_Dyz;
  abcd[996] = 2.0E0*I_NAI_H2x3z_G3yz_b-2*I_NAI_H2x3z_Dyz;
  abcd[997] = 2.0E0*I_NAI_Hx4y_G3yz_b-2*I_NAI_Hx4y_Dyz;
  abcd[998] = 2.0E0*I_NAI_Hx3yz_G3yz_b-2*I_NAI_Hx3yz_Dyz;
  abcd[999] = 2.0E0*I_NAI_Hx2y2z_G3yz_b-2*I_NAI_Hx2y2z_Dyz;
  abcd[1000] = 2.0E0*I_NAI_Hxy3z_G3yz_b-2*I_NAI_Hxy3z_Dyz;
  abcd[1001] = 2.0E0*I_NAI_Hx4z_G3yz_b-2*I_NAI_Hx4z_Dyz;
  abcd[1002] = 2.0E0*I_NAI_H5y_G3yz_b-2*I_NAI_H5y_Dyz;
  abcd[1003] = 2.0E0*I_NAI_H4yz_G3yz_b-2*I_NAI_H4yz_Dyz;
  abcd[1004] = 2.0E0*I_NAI_H3y2z_G3yz_b-2*I_NAI_H3y2z_Dyz;
  abcd[1005] = 2.0E0*I_NAI_H2y3z_G3yz_b-2*I_NAI_H2y3z_Dyz;
  abcd[1006] = 2.0E0*I_NAI_Hy4z_G3yz_b-2*I_NAI_Hy4z_Dyz;
  abcd[1007] = 2.0E0*I_NAI_H5z_G3yz_b-2*I_NAI_H5z_Dyz;
  abcd[1008] = 2.0E0*I_NAI_H5x_G2y2z_b-1*I_NAI_H5x_D2z;
  abcd[1009] = 2.0E0*I_NAI_H4xy_G2y2z_b-1*I_NAI_H4xy_D2z;
  abcd[1010] = 2.0E0*I_NAI_H4xz_G2y2z_b-1*I_NAI_H4xz_D2z;
  abcd[1011] = 2.0E0*I_NAI_H3x2y_G2y2z_b-1*I_NAI_H3x2y_D2z;
  abcd[1012] = 2.0E0*I_NAI_H3xyz_G2y2z_b-1*I_NAI_H3xyz_D2z;
  abcd[1013] = 2.0E0*I_NAI_H3x2z_G2y2z_b-1*I_NAI_H3x2z_D2z;
  abcd[1014] = 2.0E0*I_NAI_H2x3y_G2y2z_b-1*I_NAI_H2x3y_D2z;
  abcd[1015] = 2.0E0*I_NAI_H2x2yz_G2y2z_b-1*I_NAI_H2x2yz_D2z;
  abcd[1016] = 2.0E0*I_NAI_H2xy2z_G2y2z_b-1*I_NAI_H2xy2z_D2z;
  abcd[1017] = 2.0E0*I_NAI_H2x3z_G2y2z_b-1*I_NAI_H2x3z_D2z;
  abcd[1018] = 2.0E0*I_NAI_Hx4y_G2y2z_b-1*I_NAI_Hx4y_D2z;
  abcd[1019] = 2.0E0*I_NAI_Hx3yz_G2y2z_b-1*I_NAI_Hx3yz_D2z;
  abcd[1020] = 2.0E0*I_NAI_Hx2y2z_G2y2z_b-1*I_NAI_Hx2y2z_D2z;
  abcd[1021] = 2.0E0*I_NAI_Hxy3z_G2y2z_b-1*I_NAI_Hxy3z_D2z;
  abcd[1022] = 2.0E0*I_NAI_Hx4z_G2y2z_b-1*I_NAI_Hx4z_D2z;
  abcd[1023] = 2.0E0*I_NAI_H5y_G2y2z_b-1*I_NAI_H5y_D2z;
  abcd[1024] = 2.0E0*I_NAI_H4yz_G2y2z_b-1*I_NAI_H4yz_D2z;
  abcd[1025] = 2.0E0*I_NAI_H3y2z_G2y2z_b-1*I_NAI_H3y2z_D2z;
  abcd[1026] = 2.0E0*I_NAI_H2y3z_G2y2z_b-1*I_NAI_H2y3z_D2z;
  abcd[1027] = 2.0E0*I_NAI_Hy4z_G2y2z_b-1*I_NAI_Hy4z_D2z;
  abcd[1028] = 2.0E0*I_NAI_H5z_G2y2z_b-1*I_NAI_H5z_D2z;
  abcd[1029] = 2.0E0*I_NAI_H5x_Gy3z_b;
  abcd[1030] = 2.0E0*I_NAI_H4xy_Gy3z_b;
  abcd[1031] = 2.0E0*I_NAI_H4xz_Gy3z_b;
  abcd[1032] = 2.0E0*I_NAI_H3x2y_Gy3z_b;
  abcd[1033] = 2.0E0*I_NAI_H3xyz_Gy3z_b;
  abcd[1034] = 2.0E0*I_NAI_H3x2z_Gy3z_b;
  abcd[1035] = 2.0E0*I_NAI_H2x3y_Gy3z_b;
  abcd[1036] = 2.0E0*I_NAI_H2x2yz_Gy3z_b;
  abcd[1037] = 2.0E0*I_NAI_H2xy2z_Gy3z_b;
  abcd[1038] = 2.0E0*I_NAI_H2x3z_Gy3z_b;
  abcd[1039] = 2.0E0*I_NAI_Hx4y_Gy3z_b;
  abcd[1040] = 2.0E0*I_NAI_Hx3yz_Gy3z_b;
  abcd[1041] = 2.0E0*I_NAI_Hx2y2z_Gy3z_b;
  abcd[1042] = 2.0E0*I_NAI_Hxy3z_Gy3z_b;
  abcd[1043] = 2.0E0*I_NAI_Hx4z_Gy3z_b;
  abcd[1044] = 2.0E0*I_NAI_H5y_Gy3z_b;
  abcd[1045] = 2.0E0*I_NAI_H4yz_Gy3z_b;
  abcd[1046] = 2.0E0*I_NAI_H3y2z_Gy3z_b;
  abcd[1047] = 2.0E0*I_NAI_H2y3z_Gy3z_b;
  abcd[1048] = 2.0E0*I_NAI_Hy4z_Gy3z_b;
  abcd[1049] = 2.0E0*I_NAI_H5z_Gy3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_G_b
   * RHS shell quartet name: SQ_NAI_H_D
   ************************************************************/
  abcd[1050] = 2.0E0*I_NAI_H5x_G3xz_b;
  abcd[1051] = 2.0E0*I_NAI_H4xy_G3xz_b;
  abcd[1052] = 2.0E0*I_NAI_H4xz_G3xz_b;
  abcd[1053] = 2.0E0*I_NAI_H3x2y_G3xz_b;
  abcd[1054] = 2.0E0*I_NAI_H3xyz_G3xz_b;
  abcd[1055] = 2.0E0*I_NAI_H3x2z_G3xz_b;
  abcd[1056] = 2.0E0*I_NAI_H2x3y_G3xz_b;
  abcd[1057] = 2.0E0*I_NAI_H2x2yz_G3xz_b;
  abcd[1058] = 2.0E0*I_NAI_H2xy2z_G3xz_b;
  abcd[1059] = 2.0E0*I_NAI_H2x3z_G3xz_b;
  abcd[1060] = 2.0E0*I_NAI_Hx4y_G3xz_b;
  abcd[1061] = 2.0E0*I_NAI_Hx3yz_G3xz_b;
  abcd[1062] = 2.0E0*I_NAI_Hx2y2z_G3xz_b;
  abcd[1063] = 2.0E0*I_NAI_Hxy3z_G3xz_b;
  abcd[1064] = 2.0E0*I_NAI_Hx4z_G3xz_b;
  abcd[1065] = 2.0E0*I_NAI_H5y_G3xz_b;
  abcd[1066] = 2.0E0*I_NAI_H4yz_G3xz_b;
  abcd[1067] = 2.0E0*I_NAI_H3y2z_G3xz_b;
  abcd[1068] = 2.0E0*I_NAI_H2y3z_G3xz_b;
  abcd[1069] = 2.0E0*I_NAI_Hy4z_G3xz_b;
  abcd[1070] = 2.0E0*I_NAI_H5z_G3xz_b;
  abcd[1071] = 2.0E0*I_NAI_H5x_G2xyz_b;
  abcd[1072] = 2.0E0*I_NAI_H4xy_G2xyz_b;
  abcd[1073] = 2.0E0*I_NAI_H4xz_G2xyz_b;
  abcd[1074] = 2.0E0*I_NAI_H3x2y_G2xyz_b;
  abcd[1075] = 2.0E0*I_NAI_H3xyz_G2xyz_b;
  abcd[1076] = 2.0E0*I_NAI_H3x2z_G2xyz_b;
  abcd[1077] = 2.0E0*I_NAI_H2x3y_G2xyz_b;
  abcd[1078] = 2.0E0*I_NAI_H2x2yz_G2xyz_b;
  abcd[1079] = 2.0E0*I_NAI_H2xy2z_G2xyz_b;
  abcd[1080] = 2.0E0*I_NAI_H2x3z_G2xyz_b;
  abcd[1081] = 2.0E0*I_NAI_Hx4y_G2xyz_b;
  abcd[1082] = 2.0E0*I_NAI_Hx3yz_G2xyz_b;
  abcd[1083] = 2.0E0*I_NAI_Hx2y2z_G2xyz_b;
  abcd[1084] = 2.0E0*I_NAI_Hxy3z_G2xyz_b;
  abcd[1085] = 2.0E0*I_NAI_Hx4z_G2xyz_b;
  abcd[1086] = 2.0E0*I_NAI_H5y_G2xyz_b;
  abcd[1087] = 2.0E0*I_NAI_H4yz_G2xyz_b;
  abcd[1088] = 2.0E0*I_NAI_H3y2z_G2xyz_b;
  abcd[1089] = 2.0E0*I_NAI_H2y3z_G2xyz_b;
  abcd[1090] = 2.0E0*I_NAI_Hy4z_G2xyz_b;
  abcd[1091] = 2.0E0*I_NAI_H5z_G2xyz_b;
  abcd[1092] = 2.0E0*I_NAI_H5x_G2x2z_b-1*I_NAI_H5x_D2x;
  abcd[1093] = 2.0E0*I_NAI_H4xy_G2x2z_b-1*I_NAI_H4xy_D2x;
  abcd[1094] = 2.0E0*I_NAI_H4xz_G2x2z_b-1*I_NAI_H4xz_D2x;
  abcd[1095] = 2.0E0*I_NAI_H3x2y_G2x2z_b-1*I_NAI_H3x2y_D2x;
  abcd[1096] = 2.0E0*I_NAI_H3xyz_G2x2z_b-1*I_NAI_H3xyz_D2x;
  abcd[1097] = 2.0E0*I_NAI_H3x2z_G2x2z_b-1*I_NAI_H3x2z_D2x;
  abcd[1098] = 2.0E0*I_NAI_H2x3y_G2x2z_b-1*I_NAI_H2x3y_D2x;
  abcd[1099] = 2.0E0*I_NAI_H2x2yz_G2x2z_b-1*I_NAI_H2x2yz_D2x;
  abcd[1100] = 2.0E0*I_NAI_H2xy2z_G2x2z_b-1*I_NAI_H2xy2z_D2x;
  abcd[1101] = 2.0E0*I_NAI_H2x3z_G2x2z_b-1*I_NAI_H2x3z_D2x;
  abcd[1102] = 2.0E0*I_NAI_Hx4y_G2x2z_b-1*I_NAI_Hx4y_D2x;
  abcd[1103] = 2.0E0*I_NAI_Hx3yz_G2x2z_b-1*I_NAI_Hx3yz_D2x;
  abcd[1104] = 2.0E0*I_NAI_Hx2y2z_G2x2z_b-1*I_NAI_Hx2y2z_D2x;
  abcd[1105] = 2.0E0*I_NAI_Hxy3z_G2x2z_b-1*I_NAI_Hxy3z_D2x;
  abcd[1106] = 2.0E0*I_NAI_Hx4z_G2x2z_b-1*I_NAI_Hx4z_D2x;
  abcd[1107] = 2.0E0*I_NAI_H5y_G2x2z_b-1*I_NAI_H5y_D2x;
  abcd[1108] = 2.0E0*I_NAI_H4yz_G2x2z_b-1*I_NAI_H4yz_D2x;
  abcd[1109] = 2.0E0*I_NAI_H3y2z_G2x2z_b-1*I_NAI_H3y2z_D2x;
  abcd[1110] = 2.0E0*I_NAI_H2y3z_G2x2z_b-1*I_NAI_H2y3z_D2x;
  abcd[1111] = 2.0E0*I_NAI_Hy4z_G2x2z_b-1*I_NAI_Hy4z_D2x;
  abcd[1112] = 2.0E0*I_NAI_H5z_G2x2z_b-1*I_NAI_H5z_D2x;
  abcd[1113] = 2.0E0*I_NAI_H5x_Gx2yz_b;
  abcd[1114] = 2.0E0*I_NAI_H4xy_Gx2yz_b;
  abcd[1115] = 2.0E0*I_NAI_H4xz_Gx2yz_b;
  abcd[1116] = 2.0E0*I_NAI_H3x2y_Gx2yz_b;
  abcd[1117] = 2.0E0*I_NAI_H3xyz_Gx2yz_b;
  abcd[1118] = 2.0E0*I_NAI_H3x2z_Gx2yz_b;
  abcd[1119] = 2.0E0*I_NAI_H2x3y_Gx2yz_b;
  abcd[1120] = 2.0E0*I_NAI_H2x2yz_Gx2yz_b;
  abcd[1121] = 2.0E0*I_NAI_H2xy2z_Gx2yz_b;
  abcd[1122] = 2.0E0*I_NAI_H2x3z_Gx2yz_b;
  abcd[1123] = 2.0E0*I_NAI_Hx4y_Gx2yz_b;
  abcd[1124] = 2.0E0*I_NAI_Hx3yz_Gx2yz_b;
  abcd[1125] = 2.0E0*I_NAI_Hx2y2z_Gx2yz_b;
  abcd[1126] = 2.0E0*I_NAI_Hxy3z_Gx2yz_b;
  abcd[1127] = 2.0E0*I_NAI_Hx4z_Gx2yz_b;
  abcd[1128] = 2.0E0*I_NAI_H5y_Gx2yz_b;
  abcd[1129] = 2.0E0*I_NAI_H4yz_Gx2yz_b;
  abcd[1130] = 2.0E0*I_NAI_H3y2z_Gx2yz_b;
  abcd[1131] = 2.0E0*I_NAI_H2y3z_Gx2yz_b;
  abcd[1132] = 2.0E0*I_NAI_Hy4z_Gx2yz_b;
  abcd[1133] = 2.0E0*I_NAI_H5z_Gx2yz_b;
  abcd[1134] = 2.0E0*I_NAI_H5x_Gxy2z_b-1*I_NAI_H5x_Dxy;
  abcd[1135] = 2.0E0*I_NAI_H4xy_Gxy2z_b-1*I_NAI_H4xy_Dxy;
  abcd[1136] = 2.0E0*I_NAI_H4xz_Gxy2z_b-1*I_NAI_H4xz_Dxy;
  abcd[1137] = 2.0E0*I_NAI_H3x2y_Gxy2z_b-1*I_NAI_H3x2y_Dxy;
  abcd[1138] = 2.0E0*I_NAI_H3xyz_Gxy2z_b-1*I_NAI_H3xyz_Dxy;
  abcd[1139] = 2.0E0*I_NAI_H3x2z_Gxy2z_b-1*I_NAI_H3x2z_Dxy;
  abcd[1140] = 2.0E0*I_NAI_H2x3y_Gxy2z_b-1*I_NAI_H2x3y_Dxy;
  abcd[1141] = 2.0E0*I_NAI_H2x2yz_Gxy2z_b-1*I_NAI_H2x2yz_Dxy;
  abcd[1142] = 2.0E0*I_NAI_H2xy2z_Gxy2z_b-1*I_NAI_H2xy2z_Dxy;
  abcd[1143] = 2.0E0*I_NAI_H2x3z_Gxy2z_b-1*I_NAI_H2x3z_Dxy;
  abcd[1144] = 2.0E0*I_NAI_Hx4y_Gxy2z_b-1*I_NAI_Hx4y_Dxy;
  abcd[1145] = 2.0E0*I_NAI_Hx3yz_Gxy2z_b-1*I_NAI_Hx3yz_Dxy;
  abcd[1146] = 2.0E0*I_NAI_Hx2y2z_Gxy2z_b-1*I_NAI_Hx2y2z_Dxy;
  abcd[1147] = 2.0E0*I_NAI_Hxy3z_Gxy2z_b-1*I_NAI_Hxy3z_Dxy;
  abcd[1148] = 2.0E0*I_NAI_Hx4z_Gxy2z_b-1*I_NAI_Hx4z_Dxy;
  abcd[1149] = 2.0E0*I_NAI_H5y_Gxy2z_b-1*I_NAI_H5y_Dxy;
  abcd[1150] = 2.0E0*I_NAI_H4yz_Gxy2z_b-1*I_NAI_H4yz_Dxy;
  abcd[1151] = 2.0E0*I_NAI_H3y2z_Gxy2z_b-1*I_NAI_H3y2z_Dxy;
  abcd[1152] = 2.0E0*I_NAI_H2y3z_Gxy2z_b-1*I_NAI_H2y3z_Dxy;
  abcd[1153] = 2.0E0*I_NAI_Hy4z_Gxy2z_b-1*I_NAI_Hy4z_Dxy;
  abcd[1154] = 2.0E0*I_NAI_H5z_Gxy2z_b-1*I_NAI_H5z_Dxy;
  abcd[1155] = 2.0E0*I_NAI_H5x_Gx3z_b-2*I_NAI_H5x_Dxz;
  abcd[1156] = 2.0E0*I_NAI_H4xy_Gx3z_b-2*I_NAI_H4xy_Dxz;
  abcd[1157] = 2.0E0*I_NAI_H4xz_Gx3z_b-2*I_NAI_H4xz_Dxz;
  abcd[1158] = 2.0E0*I_NAI_H3x2y_Gx3z_b-2*I_NAI_H3x2y_Dxz;
  abcd[1159] = 2.0E0*I_NAI_H3xyz_Gx3z_b-2*I_NAI_H3xyz_Dxz;
  abcd[1160] = 2.0E0*I_NAI_H3x2z_Gx3z_b-2*I_NAI_H3x2z_Dxz;
  abcd[1161] = 2.0E0*I_NAI_H2x3y_Gx3z_b-2*I_NAI_H2x3y_Dxz;
  abcd[1162] = 2.0E0*I_NAI_H2x2yz_Gx3z_b-2*I_NAI_H2x2yz_Dxz;
  abcd[1163] = 2.0E0*I_NAI_H2xy2z_Gx3z_b-2*I_NAI_H2xy2z_Dxz;
  abcd[1164] = 2.0E0*I_NAI_H2x3z_Gx3z_b-2*I_NAI_H2x3z_Dxz;
  abcd[1165] = 2.0E0*I_NAI_Hx4y_Gx3z_b-2*I_NAI_Hx4y_Dxz;
  abcd[1166] = 2.0E0*I_NAI_Hx3yz_Gx3z_b-2*I_NAI_Hx3yz_Dxz;
  abcd[1167] = 2.0E0*I_NAI_Hx2y2z_Gx3z_b-2*I_NAI_Hx2y2z_Dxz;
  abcd[1168] = 2.0E0*I_NAI_Hxy3z_Gx3z_b-2*I_NAI_Hxy3z_Dxz;
  abcd[1169] = 2.0E0*I_NAI_Hx4z_Gx3z_b-2*I_NAI_Hx4z_Dxz;
  abcd[1170] = 2.0E0*I_NAI_H5y_Gx3z_b-2*I_NAI_H5y_Dxz;
  abcd[1171] = 2.0E0*I_NAI_H4yz_Gx3z_b-2*I_NAI_H4yz_Dxz;
  abcd[1172] = 2.0E0*I_NAI_H3y2z_Gx3z_b-2*I_NAI_H3y2z_Dxz;
  abcd[1173] = 2.0E0*I_NAI_H2y3z_Gx3z_b-2*I_NAI_H2y3z_Dxz;
  abcd[1174] = 2.0E0*I_NAI_Hy4z_Gx3z_b-2*I_NAI_Hy4z_Dxz;
  abcd[1175] = 2.0E0*I_NAI_H5z_Gx3z_b-2*I_NAI_H5z_Dxz;
  abcd[1176] = 2.0E0*I_NAI_H5x_G3yz_b;
  abcd[1177] = 2.0E0*I_NAI_H4xy_G3yz_b;
  abcd[1178] = 2.0E0*I_NAI_H4xz_G3yz_b;
  abcd[1179] = 2.0E0*I_NAI_H3x2y_G3yz_b;
  abcd[1180] = 2.0E0*I_NAI_H3xyz_G3yz_b;
  abcd[1181] = 2.0E0*I_NAI_H3x2z_G3yz_b;
  abcd[1182] = 2.0E0*I_NAI_H2x3y_G3yz_b;
  abcd[1183] = 2.0E0*I_NAI_H2x2yz_G3yz_b;
  abcd[1184] = 2.0E0*I_NAI_H2xy2z_G3yz_b;
  abcd[1185] = 2.0E0*I_NAI_H2x3z_G3yz_b;
  abcd[1186] = 2.0E0*I_NAI_Hx4y_G3yz_b;
  abcd[1187] = 2.0E0*I_NAI_Hx3yz_G3yz_b;
  abcd[1188] = 2.0E0*I_NAI_Hx2y2z_G3yz_b;
  abcd[1189] = 2.0E0*I_NAI_Hxy3z_G3yz_b;
  abcd[1190] = 2.0E0*I_NAI_Hx4z_G3yz_b;
  abcd[1191] = 2.0E0*I_NAI_H5y_G3yz_b;
  abcd[1192] = 2.0E0*I_NAI_H4yz_G3yz_b;
  abcd[1193] = 2.0E0*I_NAI_H3y2z_G3yz_b;
  abcd[1194] = 2.0E0*I_NAI_H2y3z_G3yz_b;
  abcd[1195] = 2.0E0*I_NAI_Hy4z_G3yz_b;
  abcd[1196] = 2.0E0*I_NAI_H5z_G3yz_b;
  abcd[1197] = 2.0E0*I_NAI_H5x_G2y2z_b-1*I_NAI_H5x_D2y;
  abcd[1198] = 2.0E0*I_NAI_H4xy_G2y2z_b-1*I_NAI_H4xy_D2y;
  abcd[1199] = 2.0E0*I_NAI_H4xz_G2y2z_b-1*I_NAI_H4xz_D2y;
  abcd[1200] = 2.0E0*I_NAI_H3x2y_G2y2z_b-1*I_NAI_H3x2y_D2y;
  abcd[1201] = 2.0E0*I_NAI_H3xyz_G2y2z_b-1*I_NAI_H3xyz_D2y;
  abcd[1202] = 2.0E0*I_NAI_H3x2z_G2y2z_b-1*I_NAI_H3x2z_D2y;
  abcd[1203] = 2.0E0*I_NAI_H2x3y_G2y2z_b-1*I_NAI_H2x3y_D2y;
  abcd[1204] = 2.0E0*I_NAI_H2x2yz_G2y2z_b-1*I_NAI_H2x2yz_D2y;
  abcd[1205] = 2.0E0*I_NAI_H2xy2z_G2y2z_b-1*I_NAI_H2xy2z_D2y;
  abcd[1206] = 2.0E0*I_NAI_H2x3z_G2y2z_b-1*I_NAI_H2x3z_D2y;
  abcd[1207] = 2.0E0*I_NAI_Hx4y_G2y2z_b-1*I_NAI_Hx4y_D2y;
  abcd[1208] = 2.0E0*I_NAI_Hx3yz_G2y2z_b-1*I_NAI_Hx3yz_D2y;
  abcd[1209] = 2.0E0*I_NAI_Hx2y2z_G2y2z_b-1*I_NAI_Hx2y2z_D2y;
  abcd[1210] = 2.0E0*I_NAI_Hxy3z_G2y2z_b-1*I_NAI_Hxy3z_D2y;
  abcd[1211] = 2.0E0*I_NAI_Hx4z_G2y2z_b-1*I_NAI_Hx4z_D2y;
  abcd[1212] = 2.0E0*I_NAI_H5y_G2y2z_b-1*I_NAI_H5y_D2y;
  abcd[1213] = 2.0E0*I_NAI_H4yz_G2y2z_b-1*I_NAI_H4yz_D2y;
  abcd[1214] = 2.0E0*I_NAI_H3y2z_G2y2z_b-1*I_NAI_H3y2z_D2y;
  abcd[1215] = 2.0E0*I_NAI_H2y3z_G2y2z_b-1*I_NAI_H2y3z_D2y;
  abcd[1216] = 2.0E0*I_NAI_Hy4z_G2y2z_b-1*I_NAI_Hy4z_D2y;
  abcd[1217] = 2.0E0*I_NAI_H5z_G2y2z_b-1*I_NAI_H5z_D2y;
  abcd[1218] = 2.0E0*I_NAI_H5x_Gy3z_b-2*I_NAI_H5x_Dyz;
  abcd[1219] = 2.0E0*I_NAI_H4xy_Gy3z_b-2*I_NAI_H4xy_Dyz;
  abcd[1220] = 2.0E0*I_NAI_H4xz_Gy3z_b-2*I_NAI_H4xz_Dyz;
  abcd[1221] = 2.0E0*I_NAI_H3x2y_Gy3z_b-2*I_NAI_H3x2y_Dyz;
  abcd[1222] = 2.0E0*I_NAI_H3xyz_Gy3z_b-2*I_NAI_H3xyz_Dyz;
  abcd[1223] = 2.0E0*I_NAI_H3x2z_Gy3z_b-2*I_NAI_H3x2z_Dyz;
  abcd[1224] = 2.0E0*I_NAI_H2x3y_Gy3z_b-2*I_NAI_H2x3y_Dyz;
  abcd[1225] = 2.0E0*I_NAI_H2x2yz_Gy3z_b-2*I_NAI_H2x2yz_Dyz;
  abcd[1226] = 2.0E0*I_NAI_H2xy2z_Gy3z_b-2*I_NAI_H2xy2z_Dyz;
  abcd[1227] = 2.0E0*I_NAI_H2x3z_Gy3z_b-2*I_NAI_H2x3z_Dyz;
  abcd[1228] = 2.0E0*I_NAI_Hx4y_Gy3z_b-2*I_NAI_Hx4y_Dyz;
  abcd[1229] = 2.0E0*I_NAI_Hx3yz_Gy3z_b-2*I_NAI_Hx3yz_Dyz;
  abcd[1230] = 2.0E0*I_NAI_Hx2y2z_Gy3z_b-2*I_NAI_Hx2y2z_Dyz;
  abcd[1231] = 2.0E0*I_NAI_Hxy3z_Gy3z_b-2*I_NAI_Hxy3z_Dyz;
  abcd[1232] = 2.0E0*I_NAI_Hx4z_Gy3z_b-2*I_NAI_Hx4z_Dyz;
  abcd[1233] = 2.0E0*I_NAI_H5y_Gy3z_b-2*I_NAI_H5y_Dyz;
  abcd[1234] = 2.0E0*I_NAI_H4yz_Gy3z_b-2*I_NAI_H4yz_Dyz;
  abcd[1235] = 2.0E0*I_NAI_H3y2z_Gy3z_b-2*I_NAI_H3y2z_Dyz;
  abcd[1236] = 2.0E0*I_NAI_H2y3z_Gy3z_b-2*I_NAI_H2y3z_Dyz;
  abcd[1237] = 2.0E0*I_NAI_Hy4z_Gy3z_b-2*I_NAI_Hy4z_Dyz;
  abcd[1238] = 2.0E0*I_NAI_H5z_Gy3z_b-2*I_NAI_H5z_Dyz;
  abcd[1239] = 2.0E0*I_NAI_H5x_G4z_b-3*I_NAI_H5x_D2z;
  abcd[1240] = 2.0E0*I_NAI_H4xy_G4z_b-3*I_NAI_H4xy_D2z;
  abcd[1241] = 2.0E0*I_NAI_H4xz_G4z_b-3*I_NAI_H4xz_D2z;
  abcd[1242] = 2.0E0*I_NAI_H3x2y_G4z_b-3*I_NAI_H3x2y_D2z;
  abcd[1243] = 2.0E0*I_NAI_H3xyz_G4z_b-3*I_NAI_H3xyz_D2z;
  abcd[1244] = 2.0E0*I_NAI_H3x2z_G4z_b-3*I_NAI_H3x2z_D2z;
  abcd[1245] = 2.0E0*I_NAI_H2x3y_G4z_b-3*I_NAI_H2x3y_D2z;
  abcd[1246] = 2.0E0*I_NAI_H2x2yz_G4z_b-3*I_NAI_H2x2yz_D2z;
  abcd[1247] = 2.0E0*I_NAI_H2xy2z_G4z_b-3*I_NAI_H2xy2z_D2z;
  abcd[1248] = 2.0E0*I_NAI_H2x3z_G4z_b-3*I_NAI_H2x3z_D2z;
  abcd[1249] = 2.0E0*I_NAI_Hx4y_G4z_b-3*I_NAI_Hx4y_D2z;
  abcd[1250] = 2.0E0*I_NAI_Hx3yz_G4z_b-3*I_NAI_Hx3yz_D2z;
  abcd[1251] = 2.0E0*I_NAI_Hx2y2z_G4z_b-3*I_NAI_Hx2y2z_D2z;
  abcd[1252] = 2.0E0*I_NAI_Hxy3z_G4z_b-3*I_NAI_Hxy3z_D2z;
  abcd[1253] = 2.0E0*I_NAI_Hx4z_G4z_b-3*I_NAI_Hx4z_D2z;
  abcd[1254] = 2.0E0*I_NAI_H5y_G4z_b-3*I_NAI_H5y_D2z;
  abcd[1255] = 2.0E0*I_NAI_H4yz_G4z_b-3*I_NAI_H4yz_D2z;
  abcd[1256] = 2.0E0*I_NAI_H3y2z_G4z_b-3*I_NAI_H3y2z_D2z;
  abcd[1257] = 2.0E0*I_NAI_H2y3z_G4z_b-3*I_NAI_H2y3z_D2z;
  abcd[1258] = 2.0E0*I_NAI_Hy4z_G4z_b-3*I_NAI_Hy4z_D2z;
  abcd[1259] = 2.0E0*I_NAI_H5z_G4z_b-3*I_NAI_H5z_D2z;
}
