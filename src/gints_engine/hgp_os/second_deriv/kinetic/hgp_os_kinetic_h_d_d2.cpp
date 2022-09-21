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
// BRA1 as redundant position, total RHS integrals evaluated as: 7420
// BRA2 as redundant position, total RHS integrals evaluated as: 6660
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

void hgp_os_kinetic_h_d_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_K7x_D2x_aa = 0.0E0;
  Double I_KINETIC_K6xy_D2x_aa = 0.0E0;
  Double I_KINETIC_K6xz_D2x_aa = 0.0E0;
  Double I_KINETIC_K5x2y_D2x_aa = 0.0E0;
  Double I_KINETIC_K5xyz_D2x_aa = 0.0E0;
  Double I_KINETIC_K5x2z_D2x_aa = 0.0E0;
  Double I_KINETIC_K4x3y_D2x_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_D2x_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_D2x_aa = 0.0E0;
  Double I_KINETIC_K4x3z_D2x_aa = 0.0E0;
  Double I_KINETIC_K3x4y_D2x_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_D2x_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_D2x_aa = 0.0E0;
  Double I_KINETIC_K3x4z_D2x_aa = 0.0E0;
  Double I_KINETIC_K2x5y_D2x_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_D2x_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_D2x_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_D2x_aa = 0.0E0;
  Double I_KINETIC_K2x5z_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx6y_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_D2x_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_D2x_aa = 0.0E0;
  Double I_KINETIC_Kx6z_D2x_aa = 0.0E0;
  Double I_KINETIC_K7y_D2x_aa = 0.0E0;
  Double I_KINETIC_K6yz_D2x_aa = 0.0E0;
  Double I_KINETIC_K5y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_K4y3z_D2x_aa = 0.0E0;
  Double I_KINETIC_K3y4z_D2x_aa = 0.0E0;
  Double I_KINETIC_K2y5z_D2x_aa = 0.0E0;
  Double I_KINETIC_Ky6z_D2x_aa = 0.0E0;
  Double I_KINETIC_K7z_D2x_aa = 0.0E0;
  Double I_KINETIC_K7x_Dxy_aa = 0.0E0;
  Double I_KINETIC_K6xy_Dxy_aa = 0.0E0;
  Double I_KINETIC_K6xz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Dxy_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Dxy_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K7y_Dxy_aa = 0.0E0;
  Double I_KINETIC_K6yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K7z_Dxy_aa = 0.0E0;
  Double I_KINETIC_K7x_Dxz_aa = 0.0E0;
  Double I_KINETIC_K6xy_Dxz_aa = 0.0E0;
  Double I_KINETIC_K6xz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Dxz_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Dxz_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K7y_Dxz_aa = 0.0E0;
  Double I_KINETIC_K6yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K7z_Dxz_aa = 0.0E0;
  Double I_KINETIC_K7x_D2y_aa = 0.0E0;
  Double I_KINETIC_K6xy_D2y_aa = 0.0E0;
  Double I_KINETIC_K6xz_D2y_aa = 0.0E0;
  Double I_KINETIC_K5x2y_D2y_aa = 0.0E0;
  Double I_KINETIC_K5xyz_D2y_aa = 0.0E0;
  Double I_KINETIC_K5x2z_D2y_aa = 0.0E0;
  Double I_KINETIC_K4x3y_D2y_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_D2y_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_D2y_aa = 0.0E0;
  Double I_KINETIC_K4x3z_D2y_aa = 0.0E0;
  Double I_KINETIC_K3x4y_D2y_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_D2y_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_D2y_aa = 0.0E0;
  Double I_KINETIC_K3x4z_D2y_aa = 0.0E0;
  Double I_KINETIC_K2x5y_D2y_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_D2y_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_D2y_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_D2y_aa = 0.0E0;
  Double I_KINETIC_K2x5z_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx6y_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_D2y_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_D2y_aa = 0.0E0;
  Double I_KINETIC_Kx6z_D2y_aa = 0.0E0;
  Double I_KINETIC_K7y_D2y_aa = 0.0E0;
  Double I_KINETIC_K6yz_D2y_aa = 0.0E0;
  Double I_KINETIC_K5y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_K4y3z_D2y_aa = 0.0E0;
  Double I_KINETIC_K3y4z_D2y_aa = 0.0E0;
  Double I_KINETIC_K2y5z_D2y_aa = 0.0E0;
  Double I_KINETIC_Ky6z_D2y_aa = 0.0E0;
  Double I_KINETIC_K7z_D2y_aa = 0.0E0;
  Double I_KINETIC_K7x_Dyz_aa = 0.0E0;
  Double I_KINETIC_K6xy_Dyz_aa = 0.0E0;
  Double I_KINETIC_K6xz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Dyz_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Dyz_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K7y_Dyz_aa = 0.0E0;
  Double I_KINETIC_K6yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K7z_Dyz_aa = 0.0E0;
  Double I_KINETIC_K7x_D2z_aa = 0.0E0;
  Double I_KINETIC_K6xy_D2z_aa = 0.0E0;
  Double I_KINETIC_K6xz_D2z_aa = 0.0E0;
  Double I_KINETIC_K5x2y_D2z_aa = 0.0E0;
  Double I_KINETIC_K5xyz_D2z_aa = 0.0E0;
  Double I_KINETIC_K5x2z_D2z_aa = 0.0E0;
  Double I_KINETIC_K4x3y_D2z_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_D2z_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_D2z_aa = 0.0E0;
  Double I_KINETIC_K4x3z_D2z_aa = 0.0E0;
  Double I_KINETIC_K3x4y_D2z_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_D2z_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_D2z_aa = 0.0E0;
  Double I_KINETIC_K3x4z_D2z_aa = 0.0E0;
  Double I_KINETIC_K2x5y_D2z_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_D2z_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_D2z_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_D2z_aa = 0.0E0;
  Double I_KINETIC_K2x5z_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx6y_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_D2z_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_D2z_aa = 0.0E0;
  Double I_KINETIC_Kx6z_D2z_aa = 0.0E0;
  Double I_KINETIC_K7y_D2z_aa = 0.0E0;
  Double I_KINETIC_K6yz_D2z_aa = 0.0E0;
  Double I_KINETIC_K5y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_K4y3z_D2z_aa = 0.0E0;
  Double I_KINETIC_K3y4z_D2z_aa = 0.0E0;
  Double I_KINETIC_K2y5z_D2z_aa = 0.0E0;
  Double I_KINETIC_Ky6z_D2z_aa = 0.0E0;
  Double I_KINETIC_K7z_D2z_aa = 0.0E0;
  Double I_KINETIC_H5x_D2x_a = 0.0E0;
  Double I_KINETIC_H4xy_D2x_a = 0.0E0;
  Double I_KINETIC_H4xz_D2x_a = 0.0E0;
  Double I_KINETIC_H3x2y_D2x_a = 0.0E0;
  Double I_KINETIC_H3xyz_D2x_a = 0.0E0;
  Double I_KINETIC_H3x2z_D2x_a = 0.0E0;
  Double I_KINETIC_H2x3y_D2x_a = 0.0E0;
  Double I_KINETIC_H2x2yz_D2x_a = 0.0E0;
  Double I_KINETIC_H2xy2z_D2x_a = 0.0E0;
  Double I_KINETIC_H2x3z_D2x_a = 0.0E0;
  Double I_KINETIC_Hx4y_D2x_a = 0.0E0;
  Double I_KINETIC_Hx3yz_D2x_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2x_a = 0.0E0;
  Double I_KINETIC_Hxy3z_D2x_a = 0.0E0;
  Double I_KINETIC_Hx4z_D2x_a = 0.0E0;
  Double I_KINETIC_H5y_D2x_a = 0.0E0;
  Double I_KINETIC_H4yz_D2x_a = 0.0E0;
  Double I_KINETIC_H3y2z_D2x_a = 0.0E0;
  Double I_KINETIC_H2y3z_D2x_a = 0.0E0;
  Double I_KINETIC_Hy4z_D2x_a = 0.0E0;
  Double I_KINETIC_H5z_D2x_a = 0.0E0;
  Double I_KINETIC_H5x_Dxy_a = 0.0E0;
  Double I_KINETIC_H4xy_Dxy_a = 0.0E0;
  Double I_KINETIC_H4xz_Dxy_a = 0.0E0;
  Double I_KINETIC_H3x2y_Dxy_a = 0.0E0;
  Double I_KINETIC_H3xyz_Dxy_a = 0.0E0;
  Double I_KINETIC_H3x2z_Dxy_a = 0.0E0;
  Double I_KINETIC_H2x3y_Dxy_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Dxy_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Dxy_a = 0.0E0;
  Double I_KINETIC_H2x3z_Dxy_a = 0.0E0;
  Double I_KINETIC_Hx4y_Dxy_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Dxy_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dxy_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Dxy_a = 0.0E0;
  Double I_KINETIC_Hx4z_Dxy_a = 0.0E0;
  Double I_KINETIC_H5y_Dxy_a = 0.0E0;
  Double I_KINETIC_H4yz_Dxy_a = 0.0E0;
  Double I_KINETIC_H3y2z_Dxy_a = 0.0E0;
  Double I_KINETIC_H2y3z_Dxy_a = 0.0E0;
  Double I_KINETIC_Hy4z_Dxy_a = 0.0E0;
  Double I_KINETIC_H5z_Dxy_a = 0.0E0;
  Double I_KINETIC_H5x_Dxz_a = 0.0E0;
  Double I_KINETIC_H4xy_Dxz_a = 0.0E0;
  Double I_KINETIC_H4xz_Dxz_a = 0.0E0;
  Double I_KINETIC_H3x2y_Dxz_a = 0.0E0;
  Double I_KINETIC_H3xyz_Dxz_a = 0.0E0;
  Double I_KINETIC_H3x2z_Dxz_a = 0.0E0;
  Double I_KINETIC_H2x3y_Dxz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Dxz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Dxz_a = 0.0E0;
  Double I_KINETIC_H2x3z_Dxz_a = 0.0E0;
  Double I_KINETIC_Hx4y_Dxz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Dxz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dxz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Dxz_a = 0.0E0;
  Double I_KINETIC_Hx4z_Dxz_a = 0.0E0;
  Double I_KINETIC_H5y_Dxz_a = 0.0E0;
  Double I_KINETIC_H4yz_Dxz_a = 0.0E0;
  Double I_KINETIC_H3y2z_Dxz_a = 0.0E0;
  Double I_KINETIC_H2y3z_Dxz_a = 0.0E0;
  Double I_KINETIC_Hy4z_Dxz_a = 0.0E0;
  Double I_KINETIC_H5z_Dxz_a = 0.0E0;
  Double I_KINETIC_H5x_D2y_a = 0.0E0;
  Double I_KINETIC_H4xy_D2y_a = 0.0E0;
  Double I_KINETIC_H4xz_D2y_a = 0.0E0;
  Double I_KINETIC_H3x2y_D2y_a = 0.0E0;
  Double I_KINETIC_H3xyz_D2y_a = 0.0E0;
  Double I_KINETIC_H3x2z_D2y_a = 0.0E0;
  Double I_KINETIC_H2x3y_D2y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_D2y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_D2y_a = 0.0E0;
  Double I_KINETIC_H2x3z_D2y_a = 0.0E0;
  Double I_KINETIC_Hx4y_D2y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_D2y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_D2y_a = 0.0E0;
  Double I_KINETIC_Hx4z_D2y_a = 0.0E0;
  Double I_KINETIC_H5y_D2y_a = 0.0E0;
  Double I_KINETIC_H4yz_D2y_a = 0.0E0;
  Double I_KINETIC_H3y2z_D2y_a = 0.0E0;
  Double I_KINETIC_H2y3z_D2y_a = 0.0E0;
  Double I_KINETIC_Hy4z_D2y_a = 0.0E0;
  Double I_KINETIC_H5z_D2y_a = 0.0E0;
  Double I_KINETIC_H5x_Dyz_a = 0.0E0;
  Double I_KINETIC_H4xy_Dyz_a = 0.0E0;
  Double I_KINETIC_H4xz_Dyz_a = 0.0E0;
  Double I_KINETIC_H3x2y_Dyz_a = 0.0E0;
  Double I_KINETIC_H3xyz_Dyz_a = 0.0E0;
  Double I_KINETIC_H3x2z_Dyz_a = 0.0E0;
  Double I_KINETIC_H2x3y_Dyz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Dyz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Dyz_a = 0.0E0;
  Double I_KINETIC_H2x3z_Dyz_a = 0.0E0;
  Double I_KINETIC_Hx4y_Dyz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Dyz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dyz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Dyz_a = 0.0E0;
  Double I_KINETIC_Hx4z_Dyz_a = 0.0E0;
  Double I_KINETIC_H5y_Dyz_a = 0.0E0;
  Double I_KINETIC_H4yz_Dyz_a = 0.0E0;
  Double I_KINETIC_H3y2z_Dyz_a = 0.0E0;
  Double I_KINETIC_H2y3z_Dyz_a = 0.0E0;
  Double I_KINETIC_Hy4z_Dyz_a = 0.0E0;
  Double I_KINETIC_H5z_Dyz_a = 0.0E0;
  Double I_KINETIC_H5x_D2z_a = 0.0E0;
  Double I_KINETIC_H4xy_D2z_a = 0.0E0;
  Double I_KINETIC_H4xz_D2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_D2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_D2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_D2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_D2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_D2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_D2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_D2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_D2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_D2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_D2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_D2z_a = 0.0E0;
  Double I_KINETIC_H5y_D2z_a = 0.0E0;
  Double I_KINETIC_H4yz_D2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_D2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_D2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_D2z_a = 0.0E0;
  Double I_KINETIC_H5z_D2z_a = 0.0E0;
  Double I_KINETIC_F3x_D2x = 0.0E0;
  Double I_KINETIC_F2xy_D2x = 0.0E0;
  Double I_KINETIC_F2xz_D2x = 0.0E0;
  Double I_KINETIC_Fx2y_D2x = 0.0E0;
  Double I_KINETIC_Fxyz_D2x = 0.0E0;
  Double I_KINETIC_Fx2z_D2x = 0.0E0;
  Double I_KINETIC_F3y_D2x = 0.0E0;
  Double I_KINETIC_F2yz_D2x = 0.0E0;
  Double I_KINETIC_Fy2z_D2x = 0.0E0;
  Double I_KINETIC_F3z_D2x = 0.0E0;
  Double I_KINETIC_F3x_Dxy = 0.0E0;
  Double I_KINETIC_F2xy_Dxy = 0.0E0;
  Double I_KINETIC_F2xz_Dxy = 0.0E0;
  Double I_KINETIC_Fx2y_Dxy = 0.0E0;
  Double I_KINETIC_Fxyz_Dxy = 0.0E0;
  Double I_KINETIC_Fx2z_Dxy = 0.0E0;
  Double I_KINETIC_F3y_Dxy = 0.0E0;
  Double I_KINETIC_F2yz_Dxy = 0.0E0;
  Double I_KINETIC_Fy2z_Dxy = 0.0E0;
  Double I_KINETIC_F3z_Dxy = 0.0E0;
  Double I_KINETIC_F3x_Dxz = 0.0E0;
  Double I_KINETIC_F2xy_Dxz = 0.0E0;
  Double I_KINETIC_F2xz_Dxz = 0.0E0;
  Double I_KINETIC_Fx2y_Dxz = 0.0E0;
  Double I_KINETIC_Fxyz_Dxz = 0.0E0;
  Double I_KINETIC_Fx2z_Dxz = 0.0E0;
  Double I_KINETIC_F3y_Dxz = 0.0E0;
  Double I_KINETIC_F2yz_Dxz = 0.0E0;
  Double I_KINETIC_Fy2z_Dxz = 0.0E0;
  Double I_KINETIC_F3z_Dxz = 0.0E0;
  Double I_KINETIC_F3x_D2y = 0.0E0;
  Double I_KINETIC_F2xy_D2y = 0.0E0;
  Double I_KINETIC_F2xz_D2y = 0.0E0;
  Double I_KINETIC_Fx2y_D2y = 0.0E0;
  Double I_KINETIC_Fxyz_D2y = 0.0E0;
  Double I_KINETIC_Fx2z_D2y = 0.0E0;
  Double I_KINETIC_F3y_D2y = 0.0E0;
  Double I_KINETIC_F2yz_D2y = 0.0E0;
  Double I_KINETIC_Fy2z_D2y = 0.0E0;
  Double I_KINETIC_F3z_D2y = 0.0E0;
  Double I_KINETIC_F3x_Dyz = 0.0E0;
  Double I_KINETIC_F2xy_Dyz = 0.0E0;
  Double I_KINETIC_F2xz_Dyz = 0.0E0;
  Double I_KINETIC_Fx2y_Dyz = 0.0E0;
  Double I_KINETIC_Fxyz_Dyz = 0.0E0;
  Double I_KINETIC_Fx2z_Dyz = 0.0E0;
  Double I_KINETIC_F3y_Dyz = 0.0E0;
  Double I_KINETIC_F2yz_Dyz = 0.0E0;
  Double I_KINETIC_Fy2z_Dyz = 0.0E0;
  Double I_KINETIC_F3z_Dyz = 0.0E0;
  Double I_KINETIC_F3x_D2z = 0.0E0;
  Double I_KINETIC_F2xy_D2z = 0.0E0;
  Double I_KINETIC_F2xz_D2z = 0.0E0;
  Double I_KINETIC_Fx2y_D2z = 0.0E0;
  Double I_KINETIC_Fxyz_D2z = 0.0E0;
  Double I_KINETIC_Fx2z_D2z = 0.0E0;
  Double I_KINETIC_F3y_D2z = 0.0E0;
  Double I_KINETIC_F2yz_D2z = 0.0E0;
  Double I_KINETIC_Fy2z_D2z = 0.0E0;
  Double I_KINETIC_F3z_D2z = 0.0E0;

  Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double xi    = alpha*beta*onedz;
    Double twoxi = 2.0E0*xi;
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    Double adz   = alpha*onedz;
    Double bdz   = beta*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double PBX   = PX - B[0];
    Double PBY   = PY - B[1];
    Double PBZ   = PZ - B[2];
    Double I_KINETIC_S_S_vrr = ic2*fbra*xi*(3.0E0-twoxi*AB2);
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_KINETIC_S_S_vrr)<THRESHOLD_MATH) continue;


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
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_S_vrr = PAZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PBX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PBX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PBX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PBY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PBY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PBY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_Px_vrr = PBX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xy_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Px_vrr = PBX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Px_vrr = PBX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Py_vrr = PBY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Py_vrr = PBY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Py_vrr = PBY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xy_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PAX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PAY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PAX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PAX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PAY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PAX*I_TWOBODYOVERLAP_F3x_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PAY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PAX*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PAX*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PAY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_Px_vrr = PAX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Px_vrr = PAY*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Px_vrr = PAX*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Px_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Px_vrr = PAX*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Px_vrr = PAY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Px_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Px_vrr = PAY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_Py_vrr = PAX*I_TWOBODYOVERLAP_G4x_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Py_vrr = PAY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Py_vrr = PAX*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Py_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Py_vrr = PAX*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Py_vrr = PAY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Py_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Py_vrr = PAY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 7 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_S_vrr = PAX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_S_vrr = PAY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_S_vrr = PAY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_S_vrr = PAX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_S_vrr = PAX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_S_vrr = PAY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_S_vrr = PAY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_S_vrr = PAZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5x_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4xy_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4xz_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5y_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4yz_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5y_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 21 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_Px_vrr = PAX*I_TWOBODYOVERLAP_H5x_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Px_vrr = PAY*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Px_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Px_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Px_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Px_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Px_vrr = PAX*I_TWOBODYOVERLAP_H5y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Px_vrr = PAX*I_TWOBODYOVERLAP_H5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Px_vrr = PAY*I_TWOBODYOVERLAP_H5y_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Px_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Px_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Px_vrr = PAY*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6z_Px_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6x_Py_vrr = PAX*I_TWOBODYOVERLAP_H5x_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Py_vrr = PAY*I_TWOBODYOVERLAP_H5x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Py_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Py_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Py_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Py_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Py_vrr = PAX*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Py_vrr = PAX*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6y_Py_vrr = PAY*I_TWOBODYOVERLAP_H5y_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Py_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Py_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Py_vrr = PAY*I_TWOBODYOVERLAP_H5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Py_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6x_Pz_vrr = PAX*I_TWOBODYOVERLAP_H5x_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Pz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Pz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Pz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Pz_vrr = PAX*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Pz_vrr = PAX*I_TWOBODYOVERLAP_H5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I6y_Pz_vrr = PAY*I_TWOBODYOVERLAP_H5y_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Pz_vrr = PAY*I_TWOBODYOVERLAP_H5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I6z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_D2x_vrr = PBX*I_TWOBODYOVERLAP_I6x_Px_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_I5xy_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_I5xz_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_I4x2y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I4x2z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_I3x3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_I3x2yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I3x3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2x_vrr = PBX*I_TWOBODYOVERLAP_I2x4y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_I2x3yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I2xy3z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I2x4z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Ix5y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Ix5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2x_vrr = PBX*I_TWOBODYOVERLAP_I6y_Px_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_I5yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I4y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I3y3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I2y4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Iy5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2x_vrr = PBX*I_TWOBODYOVERLAP_I6z_Px_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I5xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I5xz_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I4x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I4x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I3x3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I3x2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I3x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I2x4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I2x3yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I2xy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I2x4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Ix5y_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I6y_Px_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I5yz_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I4y2z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I3y3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I2y4z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Iy5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I5xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I3x2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I2x3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I2xy3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5z_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I5yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I4y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I3y3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I2y4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Iy5z_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_I6z_Px_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6x_D2y_vrr = PBY*I_TWOBODYOVERLAP_I6x_Py_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_I5xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_I5xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_I4x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I4x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_I3x3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_I3x2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I3x3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_I2x4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_I2x3yz_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I2xy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I2x4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Ix5y_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Ix5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2y_vrr = PBY*I_TWOBODYOVERLAP_I6y_Py_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_I5yz_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_H4yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I4y2z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I3y3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I2y4z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Iy5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2y_vrr = PBY*I_TWOBODYOVERLAP_I6z_Py_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I5xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I3x2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I2x3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I2xy3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5z_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I5yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I4y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I3y3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I2y4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Iy5z_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_I6z_Py_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I6x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I5xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I5xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I4x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I4x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I3x3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I3x2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I3x3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I2x4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I2x3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I2xy3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I2x4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Ix5y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Ix5z_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I6y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I5yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I4y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I3y3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I2y4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_H2y3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Iy5z_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_Hy4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_I6z_Pz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_K7x_D2x_vrr = PAX*I_TWOBODYOVERLAP_I6x_D2x_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_TWOBODYOVERLAP_K6xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_I6x_D2x_vrr;
    Double I_TWOBODYOVERLAP_K6xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I6x_D2x_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_I5xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_D2x_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x4y_Px_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_D2x_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x4z_Px_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_D2x_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_D2x_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_D2x_vrr = PAX*I_TWOBODYOVERLAP_I6y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I4y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_D2x_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I3y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_D2x_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I2y4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_D2x_vrr = PAX*I_TWOBODYOVERLAP_I6z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_TWOBODYOVERLAP_K7y_D2x_vrr = PAY*I_TWOBODYOVERLAP_I6y_D2x_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K6yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I6y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_D2x_vrr = PAY*I_TWOBODYOVERLAP_I6z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K7z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_I6z_D2x_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_K7x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I6x_Dxy_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I6x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I6x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I5xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_Py_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_Py_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2xy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_Px_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I6y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_Py_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_I6z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Py_vrr;
    Double I_TWOBODYOVERLAP_K7y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I6y_Dxy_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I6y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_Px_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_I6z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_TWOBODYOVERLAP_K7z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_I6z_Dxy_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_K7x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I6x_Dxz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Pz_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I6x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I6x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I5xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_Px_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I3x2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_Px_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I2x3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I6y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_I6z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K7y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I6y_Dxz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I6y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_Px_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_I6z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_K7z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_I6z_Dxz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_TWOBODYOVERLAP_K7x_D2y_vrr = PAX*I_TWOBODYOVERLAP_I6x_D2y_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_K6xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_I6x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_TWOBODYOVERLAP_K6xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I6x_D2y_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_I5xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I4x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I4x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_D2y_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I3x3z_Py_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_D2y_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_D2y_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I2xy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x4z_Py_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_D2y_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_D2y_vrr = PAX*I_TWOBODYOVERLAP_I6y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_D2y_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_D2y_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_D2y_vrr = PAX*I_TWOBODYOVERLAP_I6z_D2y_vrr;
    Double I_TWOBODYOVERLAP_K7y_D2y_vrr = PAY*I_TWOBODYOVERLAP_I6y_D2y_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_TWOBODYOVERLAP_K6yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I6y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I2y4z_Py_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_D2y_vrr = PAY*I_TWOBODYOVERLAP_I6z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_I6z_Py_vrr;
    Double I_TWOBODYOVERLAP_K7z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_I6z_D2y_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_K7x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I6x_Dyz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I6x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Pz_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I6x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I5xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_Py_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I3x2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_Py_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I2x3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I2xy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I6y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_I6z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_K7y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I6y_Dyz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I6y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_Py_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_I6z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K7z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_I6z_Dyz_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_I6z_Py_vrr;
    Double I_TWOBODYOVERLAP_K7x_D2z_vrr = PAX*I_TWOBODYOVERLAP_I6x_D2z_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_K6xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_I6x_D2z_vrr;
    Double I_TWOBODYOVERLAP_K6xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I6x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I6x_Pz_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_I5xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I4x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I4x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_D2z_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I3x3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I3x2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_D2z_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_D2z_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_D2z_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_D2z_vrr = PAX*I_TWOBODYOVERLAP_I6y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_D2z_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_D2z_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_D2z_vrr = PAX*I_TWOBODYOVERLAP_I6z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K7y_D2z_vrr = PAY*I_TWOBODYOVERLAP_I6y_D2z_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_K6yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I6y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I6y_Pz_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I4y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_D2z_vrr = PAY*I_TWOBODYOVERLAP_I6z_D2z_vrr;
    Double I_TWOBODYOVERLAP_K7z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_I6z_D2z_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_I6z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_Px_S_vrr = PAX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_S_vrr = PAY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_S_vrr = PAZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dxy_S_vrr = PAY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_S_vrr = PAZ*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dyz_S_vrr = PAZ*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PBX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PBX*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_Dxz_Px_vrr = PBX*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PBX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Dyz_Px_vrr = PBX*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PBX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PBY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PBY*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_Dxz_Py_vrr = PBY*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PBY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Dyz_Py_vrr = PBY*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PBY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PBZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PBZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_Dxz_Pz_vrr = PBZ*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PBZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Dyz_Pz_vrr = PBZ*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PBZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_F3x_S_vrr = PAX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_S_vrr-2*bdz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_F2xy_S_vrr = PAY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_S_vrr = PAZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_S_vrr = PAX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_S_vrr = PAZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_S_vrr = PAX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_S_vrr = PAY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_F2yz_S_vrr = PAZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_S_vrr = PAY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_S_vrr = PAZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_KINETIC_F3x_Px_vrr = PBX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_F2xy_Px_vrr = PBX*I_KINETIC_F2xy_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_Px_vrr = PBX*I_KINETIC_F2xz_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_Fx2y_Px_vrr = PBX*I_KINETIC_Fx2y_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_KINETIC_Fxyz_Px_vrr = PBX*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_KINETIC_Fx2z_Px_vrr = PBX*I_KINETIC_Fx2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_KINETIC_F3y_Px_vrr = PBX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_F2yz_Px_vrr = PBX*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_Fy2z_Px_vrr = PBX*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_KINETIC_F3z_Px_vrr = PBX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_F3x_Py_vrr = PBY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_F2xy_Py_vrr = PBY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_Py_vrr = PBY*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_Fx2y_Py_vrr = PBY*I_KINETIC_Fx2y_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_KINETIC_Fxyz_Py_vrr = PBY*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_KINETIC_Fx2z_Py_vrr = PBY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_KINETIC_F3y_Py_vrr = PBY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_F2yz_Py_vrr = PBY*I_KINETIC_F2yz_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_Fy2z_Py_vrr = PBY*I_KINETIC_Fy2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_KINETIC_F3z_Py_vrr = PBY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_F3x_Pz_vrr = PBZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_F2xy_Pz_vrr = PBZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_Pz_vrr = PBZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_Fx2y_Pz_vrr = PBZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_KINETIC_Fxyz_Pz_vrr = PBZ*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_KINETIC_Fx2z_Pz_vrr = PBZ*I_KINETIC_Fx2z_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_KINETIC_F3y_Pz_vrr = PBZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_F2yz_Pz_vrr = PBZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_Fy2z_Pz_vrr = PBZ*I_KINETIC_Fy2z_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_KINETIC_F3z_Pz_vrr = PBZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_G4x_S_vrr = PAX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G3xy_S_vrr = PAY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_S_vrr = PAZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_S_vrr = PAY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G2xyz_S_vrr = PAZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_S_vrr = PAZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Gx3y_S_vrr = PAX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_S_vrr = PAZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_S_vrr = PAY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_S_vrr = PAX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_S_vrr = PAY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_G3yz_S_vrr = PAZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_S_vrr = PAZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Gy3z_S_vrr = PAY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_S_vrr = PAZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PBX*I_KINETIC_F3x_Px_vrr+3*oned2z*I_KINETIC_D2x_Px_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PBX*I_KINETIC_F2xy_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PBX*I_KINETIC_F2xz_Px_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2x_vrr = PBX*I_KINETIC_Fx2y_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2x_vrr = PBX*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dyz_Px_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2x_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2x_vrr = PBX*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PBX*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PBX*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2x_vrr = PBX*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PBX*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PBY*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PBY*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PBY*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PBY*I_KINETIC_Fx2y_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fxyz_Dxy_vrr = PBY*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PBY*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PBY*I_KINETIC_F3y_Px_vrr+3*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PBY*I_KINETIC_F2yz_Px_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_Fy2z_Dxy_vrr = PBY*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PBY*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PBZ*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_F2xy_Dxz_vrr = PBZ*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_KINETIC_F2xz_Dxz_vrr = PBZ*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_KINETIC_Fx2y_Dxz_vrr = PBZ*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxz_vrr;
    Double I_KINETIC_Fxyz_Dxz_vrr = PBZ*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxz_vrr;
    Double I_KINETIC_Fx2z_Dxz_vrr = PBZ*I_KINETIC_Fx2z_Px_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PBZ*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_F2yz_Dxz_vrr = PBZ*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_KINETIC_Fy2z_Dxz_vrr = PBZ*I_KINETIC_Fy2z_Px_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PBZ*I_KINETIC_F3z_Px_vrr+3*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PBY*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PBY*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PBY*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2y_vrr = PBY*I_KINETIC_Fx2y_Py_vrr+2*oned2z*I_KINETIC_Dxy_Py_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2y_vrr = PBY*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Dxz_Py_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2y_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PBY*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PBY*I_KINETIC_F3y_Py_vrr+3*oned2z*I_KINETIC_D2y_Py_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PBY*I_KINETIC_F2yz_Py_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2y_vrr = PBY*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_D2z_Py_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PBY*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PBZ*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_F2xy_Dyz_vrr = PBZ*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_KINETIC_F2xz_Dyz_vrr = PBZ*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_KINETIC_Fx2y_Dyz_vrr = PBZ*I_KINETIC_Fx2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dyz_vrr;
    Double I_KINETIC_Fxyz_Dyz_vrr = PBZ*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dyz_vrr;
    Double I_KINETIC_Fx2z_Dyz_vrr = PBZ*I_KINETIC_Fx2z_Py_vrr+2*oned2z*I_KINETIC_Dxz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PBZ*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_F2yz_Dyz_vrr = PBZ*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_KINETIC_Fy2z_Dyz_vrr = PBZ*I_KINETIC_Fy2z_Py_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PBZ*I_KINETIC_F3z_Py_vrr+3*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PBZ*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PBZ*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PBZ*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2z_vrr = PBZ*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2z_vrr = PBZ*I_KINETIC_Fxyz_Pz_vrr+oned2z*I_KINETIC_Dxy_Pz_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2z_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PBZ*I_KINETIC_Fx2z_Pz_vrr+2*oned2z*I_KINETIC_Dxz_Pz_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PBZ*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PBZ*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2z_vrr = PBZ*I_KINETIC_Fy2z_Pz_vrr+2*oned2z*I_KINETIC_Dyz_Pz_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PBZ*I_KINETIC_F3z_Pz_vrr+3*oned2z*I_KINETIC_D2z_Pz_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PBX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PBX*I_KINETIC_G3xy_S_vrr+3*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PBX*I_KINETIC_G3xz_S_vrr+3*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PBX*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_Px_vrr = PBX*I_KINETIC_G2xyz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PBX*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PBX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_Px_vrr = PBX*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_Px_vrr = PBX*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PBX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PBX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PBX*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PBX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PBX*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PBX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PBY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PBY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PBY*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PBY*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_Py_vrr = PBY*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PBY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PBY*I_KINETIC_Gx3y_S_vrr+3*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_Py_vrr = PBY*I_KINETIC_Gx2yz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_Py_vrr = PBY*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PBY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PBY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PBY*I_KINETIC_G3yz_S_vrr+3*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PBY*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PBY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PBY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PBZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PBZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PBZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PBZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_Pz_vrr = PBZ*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PBZ*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PBZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_Pz_vrr = PBZ*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_Pz_vrr = PBZ*I_KINETIC_Gxy2z_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PBZ*I_KINETIC_Gx3z_S_vrr+3*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PBZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PBZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PBZ*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PBZ*I_KINETIC_Gy3z_S_vrr+3*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PBZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_H5x_S_vrr = PAX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H4xy_S_vrr = PAY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_S_vrr = PAZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_S_vrr = PAY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H3xyz_S_vrr = PAZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_S_vrr = PAZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H2x3y_S_vrr = PAX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2x2yz_S_vrr = PAZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_S_vrr = PAY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_S_vrr = PAX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hx4y_S_vrr = PAX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_S_vrr = PAZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_S_vrr = PAX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_S_vrr = PAY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_S_vrr = PAX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_S_vrr = PAY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H4yz_S_vrr = PAZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_S_vrr = PAZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2y3z_S_vrr = PAY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hy4z_S_vrr = PAY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_S_vrr = PAZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_KINETIC_H5x_Px_vrr = PBX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_KINETIC_H4xy_Px_vrr = PBX*I_KINETIC_H4xy_S_vrr+4*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_Px_vrr = PBX*I_KINETIC_H4xz_S_vrr+4*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_Px_vrr = PBX*I_KINETIC_H3x2y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_KINETIC_H3xyz_Px_vrr = PBX*I_KINETIC_H3xyz_S_vrr+3*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_KINETIC_H3x2z_Px_vrr = PBX*I_KINETIC_H3x2z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_KINETIC_H2x3y_Px_vrr = PBX*I_KINETIC_H2x3y_S_vrr+2*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_KINETIC_H2x2yz_Px_vrr = PBX*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2xy2z_Px_vrr = PBX*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_KINETIC_H2x3z_Px_vrr = PBX*I_KINETIC_H2x3z_S_vrr+2*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_KINETIC_Hx4y_Px_vrr = PBX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx3yz_Px_vrr = PBX*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Px_vrr;
    Double I_KINETIC_Hx2y2z_Px_vrr = PBX*I_KINETIC_Hx2y2z_S_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr;
    Double I_KINETIC_Hxy3z_Px_vrr = PBX*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Px_vrr;
    Double I_KINETIC_Hx4z_Px_vrr = PBX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_Px_vrr = PBX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_KINETIC_H4yz_Px_vrr = PBX*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_Px_vrr = PBX*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_KINETIC_H2y3z_Px_vrr = PBX*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_KINETIC_Hy4z_Px_vrr = PBX*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_Px_vrr = PBX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_KINETIC_H5x_Py_vrr = PBY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_KINETIC_H4xy_Py_vrr = PBY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_Py_vrr = PBY*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_Py_vrr = PBY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_KINETIC_H3xyz_Py_vrr = PBY*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Py_vrr;
    Double I_KINETIC_H3x2z_Py_vrr = PBY*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_KINETIC_H2x3y_Py_vrr = PBY*I_KINETIC_H2x3y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_KINETIC_H2x2yz_Py_vrr = PBY*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2xy2z_Py_vrr = PBY*I_KINETIC_H2xy2z_S_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_KINETIC_H2x3z_Py_vrr = PBY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_KINETIC_Hx4y_Py_vrr = PBY*I_KINETIC_Hx4y_S_vrr+4*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx3yz_Py_vrr = PBY*I_KINETIC_Hx3yz_S_vrr+3*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Py_vrr;
    Double I_KINETIC_Hx2y2z_Py_vrr = PBY*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr;
    Double I_KINETIC_Hxy3z_Py_vrr = PBY*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Py_vrr;
    Double I_KINETIC_Hx4z_Py_vrr = PBY*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_Py_vrr = PBY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_KINETIC_H4yz_Py_vrr = PBY*I_KINETIC_H4yz_S_vrr+4*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_Py_vrr = PBY*I_KINETIC_H3y2z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_KINETIC_H2y3z_Py_vrr = PBY*I_KINETIC_H2y3z_S_vrr+2*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_KINETIC_Hy4z_Py_vrr = PBY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_Py_vrr = PBY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_KINETIC_H5x_Pz_vrr = PBZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_KINETIC_H4xy_Pz_vrr = PBZ*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_Pz_vrr = PBZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_Pz_vrr = PBZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_KINETIC_H3xyz_Pz_vrr = PBZ*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Pz_vrr;
    Double I_KINETIC_H3x2z_Pz_vrr = PBZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_KINETIC_H2x3y_Pz_vrr = PBZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_Pz_vrr = PBZ*I_KINETIC_H2x2yz_S_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2xy2z_Pz_vrr = PBZ*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr;
    Double I_KINETIC_H2x3z_Pz_vrr = PBZ*I_KINETIC_H2x3z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_KINETIC_Hx4y_Pz_vrr = PBZ*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx3yz_Pz_vrr = PBZ*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr;
    Double I_KINETIC_Hx2y2z_Pz_vrr = PBZ*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr;
    Double I_KINETIC_Hxy3z_Pz_vrr = PBZ*I_KINETIC_Hxy3z_S_vrr+3*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr;
    Double I_KINETIC_Hx4z_Pz_vrr = PBZ*I_KINETIC_Hx4z_S_vrr+4*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_Pz_vrr = PBZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_KINETIC_H4yz_Pz_vrr = PBZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_Pz_vrr = PBZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_KINETIC_H2y3z_Pz_vrr = PBZ*I_KINETIC_H2y3z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_KINETIC_Hy4z_Pz_vrr = PBZ*I_KINETIC_Hy4z_S_vrr+4*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_Pz_vrr = PBZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 7 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_I6x_S_vrr = PAX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I5xy_S_vrr = PAY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_S_vrr = PAZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_S_vrr = PAY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I4x2z_S_vrr = PAZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I3x3y_S_vrr = PAY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_I3x2yz_S_vrr = PAZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3x3z_S_vrr = PAZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_I2x4y_S_vrr = PAX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I2x3yz_S_vrr = PAZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2xy3z_S_vrr = PAY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_S_vrr = PAX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Ix5y_S_vrr = PAX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix5z_S_vrr = PAX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_S_vrr = PAY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I5yz_S_vrr = PAZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_S_vrr = PAZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I3y3z_S_vrr = PAZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_I2y4z_S_vrr = PAY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Iy5z_S_vrr = PAY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_S_vrr = PAZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_KINETIC_H5x_D2x_vrr = PBX*I_KINETIC_H5x_Px_vrr+5*oned2z*I_KINETIC_G4x_Px_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2x_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2x_vrr = PBX*I_KINETIC_H4xy_Px_vrr+4*oned2z*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2x_vrr = PBX*I_KINETIC_H4xz_Px_vrr+4*oned2z*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2x_vrr = PBX*I_KINETIC_H3x2y_Px_vrr+3*oned2z*I_KINETIC_G2x2y_Px_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2x_vrr = PBX*I_KINETIC_H3xyz_Px_vrr+3*oned2z*I_KINETIC_G2xyz_Px_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2x_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2x_vrr = PBX*I_KINETIC_H3x2z_Px_vrr+3*oned2z*I_KINETIC_G2x2z_Px_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2x_vrr = PBX*I_KINETIC_H2x3y_Px_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2x_vrr = PBX*I_KINETIC_H2x2yz_Px_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2x_vrr = PBX*I_KINETIC_H2xy2z_Px_vrr+2*oned2z*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2x_vrr = PBX*I_KINETIC_H2x3z_Px_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2x_vrr = PBX*I_KINETIC_Hx4y_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2x_vrr = PBX*I_KINETIC_Hx3yz_Px_vrr+oned2z*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2x_vrr = PBX*I_KINETIC_Hx2y2z_Px_vrr+oned2z*I_KINETIC_G2y2z_Px_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2x_vrr = PBX*I_KINETIC_Hxy3z_Px_vrr+oned2z*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2x_vrr = PBX*I_KINETIC_Hx4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2x_vrr = PBX*I_KINETIC_H5y_Px_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2x_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2x_vrr = PBX*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2x_vrr = PBX*I_KINETIC_H3y2z_Px_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2x_vrr = PBX*I_KINETIC_H2y3z_Px_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2x_vrr = PBX*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2x_vrr = PBX*I_KINETIC_H5z_Px_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2x_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_Dxy_vrr = PBY*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_KINETIC_H4xy_Dxy_vrr = PBY*I_KINETIC_H4xy_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_KINETIC_H4xz_Dxy_vrr = PBY*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_KINETIC_H3x2y_Dxy_vrr = PBY*I_KINETIC_H3x2y_Px_vrr+2*oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_KINETIC_H3xyz_Dxy_vrr = PBY*I_KINETIC_H3xyz_Px_vrr+oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dxy_vrr;
    Double I_KINETIC_H3x2z_Dxy_vrr = PBY*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr;
    Double I_KINETIC_H2x3y_Dxy_vrr = PBY*I_KINETIC_H2x3y_Px_vrr+3*oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_KINETIC_H2x2yz_Dxy_vrr = PBY*I_KINETIC_H2x2yz_Px_vrr+2*oned2z*I_KINETIC_G2xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr;
    Double I_KINETIC_H2xy2z_Dxy_vrr = PBY*I_KINETIC_H2xy2z_Px_vrr+oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dxy_vrr;
    Double I_KINETIC_H2x3z_Dxy_vrr = PBY*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr;
    Double I_KINETIC_Hx4y_Dxy_vrr = PBY*I_KINETIC_Hx4y_Px_vrr+4*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_KINETIC_Hx3yz_Dxy_vrr = PBY*I_KINETIC_Hx3yz_Px_vrr+3*oned2z*I_KINETIC_Gx2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dxy_vrr;
    Double I_KINETIC_Hx2y2z_Dxy_vrr = PBY*I_KINETIC_Hx2y2z_Px_vrr+2*oned2z*I_KINETIC_Gxy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dxy_vrr;
    Double I_KINETIC_Hxy3z_Dxy_vrr = PBY*I_KINETIC_Hxy3z_Px_vrr+oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dxy_vrr;
    Double I_KINETIC_Hx4z_Dxy_vrr = PBY*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_KINETIC_H5y_Dxy_vrr = PBY*I_KINETIC_H5y_Px_vrr+5*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_KINETIC_H4yz_Dxy_vrr = PBY*I_KINETIC_H4yz_Px_vrr+4*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_KINETIC_H3y2z_Dxy_vrr = PBY*I_KINETIC_H3y2z_Px_vrr+3*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr;
    Double I_KINETIC_H2y3z_Dxy_vrr = PBY*I_KINETIC_H2y3z_Px_vrr+2*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr;
    Double I_KINETIC_Hy4z_Dxy_vrr = PBY*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_KINETIC_H5z_Dxy_vrr = PBY*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_KINETIC_H5x_Dxz_vrr = PBZ*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_KINETIC_H4xy_Dxz_vrr = PBZ*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_KINETIC_H4xz_Dxz_vrr = PBZ*I_KINETIC_H4xz_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_KINETIC_H3x2y_Dxz_vrr = PBZ*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr;
    Double I_KINETIC_H3xyz_Dxz_vrr = PBZ*I_KINETIC_H3xyz_Px_vrr+oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dxz_vrr;
    Double I_KINETIC_H3x2z_Dxz_vrr = PBZ*I_KINETIC_H3x2z_Px_vrr+2*oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr;
    Double I_KINETIC_H2x3y_Dxz_vrr = PBZ*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr;
    Double I_KINETIC_H2x2yz_Dxz_vrr = PBZ*I_KINETIC_H2x2yz_Px_vrr+oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr;
    Double I_KINETIC_H2xy2z_Dxz_vrr = PBZ*I_KINETIC_H2xy2z_Px_vrr+2*oned2z*I_KINETIC_G2xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dxz_vrr;
    Double I_KINETIC_H2x3z_Dxz_vrr = PBZ*I_KINETIC_H2x3z_Px_vrr+3*oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_KINETIC_Hx4y_Dxz_vrr = PBZ*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_KINETIC_Hx3yz_Dxz_vrr = PBZ*I_KINETIC_Hx3yz_Px_vrr+oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dxz_vrr;
    Double I_KINETIC_Hx2y2z_Dxz_vrr = PBZ*I_KINETIC_Hx2y2z_Px_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dxz_vrr;
    Double I_KINETIC_Hxy3z_Dxz_vrr = PBZ*I_KINETIC_Hxy3z_Px_vrr+3*oned2z*I_KINETIC_Gxy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dxz_vrr;
    Double I_KINETIC_Hx4z_Dxz_vrr = PBZ*I_KINETIC_Hx4z_Px_vrr+4*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_KINETIC_H5y_Dxz_vrr = PBZ*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_KINETIC_H4yz_Dxz_vrr = PBZ*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_KINETIC_H3y2z_Dxz_vrr = PBZ*I_KINETIC_H3y2z_Px_vrr+2*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr;
    Double I_KINETIC_H2y3z_Dxz_vrr = PBZ*I_KINETIC_H2y3z_Px_vrr+3*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr;
    Double I_KINETIC_Hy4z_Dxz_vrr = PBZ*I_KINETIC_Hy4z_Px_vrr+4*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_KINETIC_H5z_Dxz_vrr = PBZ*I_KINETIC_H5z_Px_vrr+5*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_KINETIC_H5x_D2y_vrr = PBY*I_KINETIC_H5x_Py_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2y_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2y_vrr = PBY*I_KINETIC_H4xy_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2y_vrr = PBY*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2y_vrr = PBY*I_KINETIC_H3x2y_Py_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2y_vrr = PBY*I_KINETIC_H3xyz_Py_vrr+oned2z*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2y_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2y_vrr = PBY*I_KINETIC_H3x2z_Py_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2y_vrr = PBY*I_KINETIC_H2x3y_Py_vrr+3*oned2z*I_KINETIC_G2x2y_Py_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2y_vrr = PBY*I_KINETIC_H2x2yz_Py_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2y_vrr = PBY*I_KINETIC_H2xy2z_Py_vrr+oned2z*I_KINETIC_G2x2z_Py_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2y_vrr = PBY*I_KINETIC_H2x3z_Py_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2y_vrr = PBY*I_KINETIC_Hx4y_Py_vrr+4*oned2z*I_KINETIC_Gx3y_Py_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2y_vrr = PBY*I_KINETIC_Hx3yz_Py_vrr+3*oned2z*I_KINETIC_Gx2yz_Py_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2y_vrr = PBY*I_KINETIC_Hx2y2z_Py_vrr+2*oned2z*I_KINETIC_Gxy2z_Py_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2y_vrr = PBY*I_KINETIC_Hxy3z_Py_vrr+oned2z*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2y_vrr = PBY*I_KINETIC_Hx4z_Py_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2y_vrr = PBY*I_KINETIC_H5y_Py_vrr+5*oned2z*I_KINETIC_G4y_Py_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2y_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2y_vrr = PBY*I_KINETIC_H4yz_Py_vrr+4*oned2z*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2y_vrr = PBY*I_KINETIC_H3y2z_Py_vrr+3*oned2z*I_KINETIC_G2y2z_Py_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2y_vrr = PBY*I_KINETIC_H2y3z_Py_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2y_vrr = PBY*I_KINETIC_Hy4z_Py_vrr+oned2z*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2y_vrr = PBY*I_KINETIC_H5z_Py_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2y_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_Dyz_vrr = PBZ*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_KINETIC_H4xy_Dyz_vrr = PBZ*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_KINETIC_H4xz_Dyz_vrr = PBZ*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_KINETIC_H3x2y_Dyz_vrr = PBZ*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr;
    Double I_KINETIC_H3xyz_Dyz_vrr = PBZ*I_KINETIC_H3xyz_Py_vrr+oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dyz_vrr;
    Double I_KINETIC_H3x2z_Dyz_vrr = PBZ*I_KINETIC_H3x2z_Py_vrr+2*oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr;
    Double I_KINETIC_H2x3y_Dyz_vrr = PBZ*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr;
    Double I_KINETIC_H2x2yz_Dyz_vrr = PBZ*I_KINETIC_H2x2yz_Py_vrr+oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr;
    Double I_KINETIC_H2xy2z_Dyz_vrr = PBZ*I_KINETIC_H2xy2z_Py_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dyz_vrr;
    Double I_KINETIC_H2x3z_Dyz_vrr = PBZ*I_KINETIC_H2x3z_Py_vrr+3*oned2z*I_KINETIC_G2x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr;
    Double I_KINETIC_Hx4y_Dyz_vrr = PBZ*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_KINETIC_Hx3yz_Dyz_vrr = PBZ*I_KINETIC_Hx3yz_Py_vrr+oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dyz_vrr;
    Double I_KINETIC_Hx2y2z_Dyz_vrr = PBZ*I_KINETIC_Hx2y2z_Py_vrr+2*oned2z*I_KINETIC_Gx2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dyz_vrr;
    Double I_KINETIC_Hxy3z_Dyz_vrr = PBZ*I_KINETIC_Hxy3z_Py_vrr+3*oned2z*I_KINETIC_Gxy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dyz_vrr;
    Double I_KINETIC_Hx4z_Dyz_vrr = PBZ*I_KINETIC_Hx4z_Py_vrr+4*oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_KINETIC_H5y_Dyz_vrr = PBZ*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_KINETIC_H4yz_Dyz_vrr = PBZ*I_KINETIC_H4yz_Py_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_KINETIC_H3y2z_Dyz_vrr = PBZ*I_KINETIC_H3y2z_Py_vrr+2*oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr;
    Double I_KINETIC_H2y3z_Dyz_vrr = PBZ*I_KINETIC_H2y3z_Py_vrr+3*oned2z*I_KINETIC_G2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr;
    Double I_KINETIC_Hy4z_Dyz_vrr = PBZ*I_KINETIC_Hy4z_Py_vrr+4*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_KINETIC_H5z_Dyz_vrr = PBZ*I_KINETIC_H5z_Py_vrr+5*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_KINETIC_H5x_D2z_vrr = PBZ*I_KINETIC_H5x_Pz_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2z_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2z_vrr = PBZ*I_KINETIC_H4xy_Pz_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2z_vrr = PBZ*I_KINETIC_H4xz_Pz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2z_vrr = PBZ*I_KINETIC_H3x2y_Pz_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2z_vrr = PBZ*I_KINETIC_H3xyz_Pz_vrr+oned2z*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2z_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2z_vrr = PBZ*I_KINETIC_H3x2z_Pz_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2z_vrr = PBZ*I_KINETIC_H2x3y_Pz_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2z_vrr = PBZ*I_KINETIC_H2x2yz_Pz_vrr+oned2z*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2z_vrr = PBZ*I_KINETIC_H2xy2z_Pz_vrr+2*oned2z*I_KINETIC_G2xyz_Pz_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2z_vrr = PBZ*I_KINETIC_H2x3z_Pz_vrr+3*oned2z*I_KINETIC_G2x2z_Pz_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2z_vrr = PBZ*I_KINETIC_Hx4y_Pz_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2z_vrr = PBZ*I_KINETIC_Hx3yz_Pz_vrr+oned2z*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2z_vrr = PBZ*I_KINETIC_Hx2y2z_Pz_vrr+2*oned2z*I_KINETIC_Gx2yz_Pz_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2z_vrr = PBZ*I_KINETIC_Hxy3z_Pz_vrr+3*oned2z*I_KINETIC_Gxy2z_Pz_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2z_vrr = PBZ*I_KINETIC_Hx4z_Pz_vrr+4*oned2z*I_KINETIC_Gx3z_Pz_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2z_vrr = PBZ*I_KINETIC_H5y_Pz_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2z_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2z_vrr = PBZ*I_KINETIC_H4yz_Pz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2z_vrr = PBZ*I_KINETIC_H3y2z_Pz_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2z_vrr = PBZ*I_KINETIC_H2y3z_Pz_vrr+3*oned2z*I_KINETIC_G2y2z_Pz_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2z_vrr = PBZ*I_KINETIC_Hy4z_Pz_vrr+4*oned2z*I_KINETIC_Gy3z_Pz_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2z_vrr = PBZ*I_KINETIC_H5z_Pz_vrr+5*oned2z*I_KINETIC_G4z_Pz_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2z_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 21 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_I6x_Px_vrr = PAX*I_KINETIC_H5x_Px_vrr+5*oned2z*I_KINETIC_G4x_Px_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Px_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_I5xy_Px_vrr = PAY*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_KINETIC_I5xz_Px_vrr = PAZ*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Px_vrr;
    Double I_KINETIC_I4x2y_Px_vrr = PAY*I_KINETIC_H4xy_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Px_vrr-bdz*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_I4x2z_Px_vrr = PAZ*I_KINETIC_H4xz_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Px_vrr-bdz*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_I3x3y_Px_vrr = PAY*I_KINETIC_H3x2y_Px_vrr+2*oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Px_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_I3x2yz_Px_vrr = PAZ*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Px_vrr;
    Double I_KINETIC_I3x3z_Px_vrr = PAZ*I_KINETIC_H3x2z_Px_vrr+2*oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Px_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_I2x4y_Px_vrr = PAX*I_KINETIC_Hx4y_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Px_vrr-bdz*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_I2x3yz_Px_vrr = PAZ*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Px_vrr;
    Double I_KINETIC_I2xy3z_Px_vrr = PAY*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Px_vrr;
    Double I_KINETIC_I2x4z_Px_vrr = PAX*I_KINETIC_Hx4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Px_vrr-bdz*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_Ix5y_Px_vrr = PAX*I_KINETIC_H5y_Px_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_KINETIC_Ix5z_Px_vrr = PAX*I_KINETIC_H5z_Px_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_KINETIC_I6y_Px_vrr = PAY*I_KINETIC_H5y_Px_vrr+5*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Px_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_I5yz_Px_vrr = PAZ*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Px_vrr;
    Double I_KINETIC_I4y2z_Px_vrr = PAZ*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Px_vrr-bdz*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_I3y3z_Px_vrr = PAZ*I_KINETIC_H3y2z_Px_vrr+2*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Px_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_I2y4z_Px_vrr = PAY*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Px_vrr-bdz*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_Iy5z_Px_vrr = PAY*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Px_vrr;
    Double I_KINETIC_I6z_Px_vrr = PAZ*I_KINETIC_H5z_Px_vrr+5*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Px_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_I6x_Py_vrr = PAX*I_KINETIC_H5x_Py_vrr+5*oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Py_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_I5xy_Py_vrr = PAY*I_KINETIC_H5x_Py_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_KINETIC_I5xz_Py_vrr = PAZ*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Py_vrr;
    Double I_KINETIC_I4x2y_Py_vrr = PAY*I_KINETIC_H4xy_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Py_vrr-bdz*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_I4x2z_Py_vrr = PAZ*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Py_vrr-bdz*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_I3x3y_Py_vrr = PAY*I_KINETIC_H3x2y_Py_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Py_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_I3x2yz_Py_vrr = PAZ*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Py_vrr;
    Double I_KINETIC_I3x3z_Py_vrr = PAZ*I_KINETIC_H3x2z_Py_vrr+2*oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Py_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_I2x4y_Py_vrr = PAX*I_KINETIC_Hx4y_Py_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Py_vrr-bdz*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_I2x3yz_Py_vrr = PAZ*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Py_vrr;
    Double I_KINETIC_I2xy3z_Py_vrr = PAY*I_KINETIC_H2x3z_Py_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Py_vrr;
    Double I_KINETIC_I2x4z_Py_vrr = PAX*I_KINETIC_Hx4z_Py_vrr+oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Py_vrr-bdz*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_Ix5y_Py_vrr = PAX*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_KINETIC_Ix5z_Py_vrr = PAX*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Py_vrr;
    Double I_KINETIC_I6y_Py_vrr = PAY*I_KINETIC_H5y_Py_vrr+5*oned2z*I_KINETIC_G4y_Py_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Py_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_I5yz_Py_vrr = PAZ*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Py_vrr;
    Double I_KINETIC_I4y2z_Py_vrr = PAZ*I_KINETIC_H4yz_Py_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Py_vrr-bdz*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_I3y3z_Py_vrr = PAZ*I_KINETIC_H3y2z_Py_vrr+2*oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Py_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_I2y4z_Py_vrr = PAY*I_KINETIC_Hy4z_Py_vrr+oned2z*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Py_vrr-bdz*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_Iy5z_Py_vrr = PAY*I_KINETIC_H5z_Py_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Py_vrr;
    Double I_KINETIC_I6z_Py_vrr = PAZ*I_KINETIC_H5z_Py_vrr+5*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Py_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_I6x_Pz_vrr = PAX*I_KINETIC_H5x_Pz_vrr+5*oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Pz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_I5xy_Pz_vrr = PAY*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Pz_vrr;
    Double I_KINETIC_I5xz_Pz_vrr = PAZ*I_KINETIC_H5x_Pz_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Pz_vrr;
    Double I_KINETIC_I4x2y_Pz_vrr = PAY*I_KINETIC_H4xy_Pz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_I4x2z_Pz_vrr = PAZ*I_KINETIC_H4xz_Pz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_I3x3y_Pz_vrr = PAY*I_KINETIC_H3x2y_Pz_vrr+2*oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_I3x2yz_Pz_vrr = PAZ*I_KINETIC_H3x2y_Pz_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Pz_vrr;
    Double I_KINETIC_I3x3z_Pz_vrr = PAZ*I_KINETIC_H3x2z_Pz_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_I2x4y_Pz_vrr = PAX*I_KINETIC_Hx4y_Pz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_I2x3yz_Pz_vrr = PAZ*I_KINETIC_H2x3y_Pz_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Pz_vrr;
    Double I_KINETIC_I2xy3z_Pz_vrr = PAY*I_KINETIC_H2x3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Pz_vrr;
    Double I_KINETIC_I2x4z_Pz_vrr = PAX*I_KINETIC_Hx4z_Pz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_KINETIC_Ix5y_Pz_vrr = PAX*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Pz_vrr;
    Double I_KINETIC_Ix5z_Pz_vrr = PAX*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Pz_vrr;
    Double I_KINETIC_I6y_Pz_vrr = PAY*I_KINETIC_H5y_Pz_vrr+5*oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Pz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_I5yz_Pz_vrr = PAZ*I_KINETIC_H5y_Pz_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Pz_vrr;
    Double I_KINETIC_I4y2z_Pz_vrr = PAZ*I_KINETIC_H4yz_Pz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_I3y3z_Pz_vrr = PAZ*I_KINETIC_H3y2z_Pz_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_I2y4z_Pz_vrr = PAY*I_KINETIC_Hy4z_Pz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Pz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_KINETIC_Iy5z_Pz_vrr = PAY*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Pz_vrr;
    Double I_KINETIC_I6z_Pz_vrr = PAZ*I_KINETIC_H5z_Pz_vrr+5*oned2z*I_KINETIC_G4z_Pz_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Pz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_I_P
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_KINETIC_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_KINETIC_I6x_D2x_vrr = PBX*I_KINETIC_I6x_Px_vrr+6*oned2z*I_KINETIC_H5x_Px_vrr+oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2x_vrr-adz*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_KINETIC_I5xy_D2x_vrr = PBX*I_KINETIC_I5xy_Px_vrr+5*oned2z*I_KINETIC_H4xy_Px_vrr+oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2x_vrr-adz*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_D2x_vrr = PBX*I_KINETIC_I5xz_Px_vrr+5*oned2z*I_KINETIC_H4xz_Px_vrr+oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2x_vrr-adz*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_D2x_vrr = PBX*I_KINETIC_I4x2y_Px_vrr+4*oned2z*I_KINETIC_H3x2y_Px_vrr+oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_KINETIC_I4x2z_D2x_vrr = PBX*I_KINETIC_I4x2z_Px_vrr+4*oned2z*I_KINETIC_H3x2z_Px_vrr+oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_KINETIC_I3x3y_D2x_vrr = PBX*I_KINETIC_I3x3y_Px_vrr+3*oned2z*I_KINETIC_H2x3y_Px_vrr+oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2x_vrr-adz*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_KINETIC_I3x2yz_D2x_vrr = PBX*I_KINETIC_I3x2yz_Px_vrr+3*oned2z*I_KINETIC_H2x2yz_Px_vrr+oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3x3z_D2x_vrr = PBX*I_KINETIC_I3x3z_Px_vrr+3*oned2z*I_KINETIC_H2x3z_Px_vrr+oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2x_vrr-adz*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_KINETIC_I2x4y_D2x_vrr = PBX*I_KINETIC_I2x4y_Px_vrr+2*oned2z*I_KINETIC_Hx4y_Px_vrr+oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2x_vrr-adz*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_KINETIC_I2x3yz_D2x_vrr = PBX*I_KINETIC_I2x3yz_Px_vrr+2*oned2z*I_KINETIC_Hx3yz_Px_vrr+oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2x_vrr-adz*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2xy3z_D2x_vrr = PBX*I_KINETIC_I2xy3z_Px_vrr+2*oned2z*I_KINETIC_Hxy3z_Px_vrr+oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2x_vrr-adz*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_D2x_vrr = PBX*I_KINETIC_I2x4z_Px_vrr+2*oned2z*I_KINETIC_Hx4z_Px_vrr+oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2x_vrr-adz*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_KINETIC_Ix5y_D2x_vrr = PBX*I_KINETIC_Ix5y_Px_vrr+oned2z*I_KINETIC_H5y_Px_vrr+oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2x_vrr-adz*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix5z_D2x_vrr = PBX*I_KINETIC_Ix5z_Px_vrr+oned2z*I_KINETIC_H5z_Px_vrr+oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2x_vrr-adz*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_D2x_vrr = PBX*I_KINETIC_I6y_Px_vrr+oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2x_vrr-adz*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_KINETIC_I5yz_D2x_vrr = PBX*I_KINETIC_I5yz_Px_vrr+oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2x_vrr-adz*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_D2x_vrr = PBX*I_KINETIC_I4y2z_Px_vrr+oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_KINETIC_I3y3z_D2x_vrr = PBX*I_KINETIC_I3y3z_Px_vrr+oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2x_vrr-adz*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_KINETIC_I2y4z_D2x_vrr = PBX*I_KINETIC_I2y4z_Px_vrr+oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2x_vrr-adz*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_KINETIC_Iy5z_D2x_vrr = PBX*I_KINETIC_Iy5z_Px_vrr+oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2x_vrr-adz*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_D2x_vrr = PBX*I_KINETIC_I6z_Px_vrr+oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2x_vrr-adz*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_KINETIC_I6x_Dxy_vrr = PBY*I_KINETIC_I6x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dxy_vrr;
    Double I_KINETIC_I5xy_Dxy_vrr = PBY*I_KINETIC_I5xy_Px_vrr+oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dxy_vrr;
    Double I_KINETIC_I5xz_Dxy_vrr = PBY*I_KINETIC_I5xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dxy_vrr;
    Double I_KINETIC_I4x2y_Dxy_vrr = PBY*I_KINETIC_I4x2y_Px_vrr+2*oned2z*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dxy_vrr;
    Double I_KINETIC_I4x2z_Dxy_vrr = PBY*I_KINETIC_I4x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dxy_vrr;
    Double I_KINETIC_I3x3y_Dxy_vrr = PBY*I_KINETIC_I3x3y_Px_vrr+3*oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dxy_vrr;
    Double I_KINETIC_I3x2yz_Dxy_vrr = PBY*I_KINETIC_I3x2yz_Px_vrr+2*oned2z*I_KINETIC_H3xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dxy_vrr;
    Double I_KINETIC_I3x3z_Dxy_vrr = PBY*I_KINETIC_I3x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dxy_vrr;
    Double I_KINETIC_I2x4y_Dxy_vrr = PBY*I_KINETIC_I2x4y_Px_vrr+4*oned2z*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dxy_vrr;
    Double I_KINETIC_I2x3yz_Dxy_vrr = PBY*I_KINETIC_I2x3yz_Px_vrr+3*oned2z*I_KINETIC_H2x2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dxy_vrr;
    Double I_KINETIC_I2xy3z_Dxy_vrr = PBY*I_KINETIC_I2xy3z_Px_vrr+oned2z*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dxy_vrr;
    Double I_KINETIC_I2x4z_Dxy_vrr = PBY*I_KINETIC_I2x4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dxy_vrr;
    Double I_KINETIC_Ix5y_Dxy_vrr = PBY*I_KINETIC_Ix5y_Px_vrr+5*oned2z*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dxy_vrr;
    Double I_KINETIC_Ix5z_Dxy_vrr = PBY*I_KINETIC_Ix5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dxy_vrr;
    Double I_KINETIC_I6y_Dxy_vrr = PBY*I_KINETIC_I6y_Px_vrr+6*oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dxy_vrr;
    Double I_KINETIC_I5yz_Dxy_vrr = PBY*I_KINETIC_I5yz_Px_vrr+5*oned2z*I_KINETIC_H4yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dxy_vrr;
    Double I_KINETIC_I4y2z_Dxy_vrr = PBY*I_KINETIC_I4y2z_Px_vrr+4*oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dxy_vrr;
    Double I_KINETIC_I3y3z_Dxy_vrr = PBY*I_KINETIC_I3y3z_Px_vrr+3*oned2z*I_KINETIC_H2y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dxy_vrr;
    Double I_KINETIC_I2y4z_Dxy_vrr = PBY*I_KINETIC_I2y4z_Px_vrr+2*oned2z*I_KINETIC_Hy4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dxy_vrr;
    Double I_KINETIC_Iy5z_Dxy_vrr = PBY*I_KINETIC_Iy5z_Px_vrr+oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dxy_vrr;
    Double I_KINETIC_I6z_Dxy_vrr = PBY*I_KINETIC_I6z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dxy_vrr;
    Double I_KINETIC_I6x_Dxz_vrr = PBZ*I_KINETIC_I6x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dxz_vrr;
    Double I_KINETIC_I5xy_Dxz_vrr = PBZ*I_KINETIC_I5xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dxz_vrr;
    Double I_KINETIC_I5xz_Dxz_vrr = PBZ*I_KINETIC_I5xz_Px_vrr+oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dxz_vrr;
    Double I_KINETIC_I4x2y_Dxz_vrr = PBZ*I_KINETIC_I4x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dxz_vrr;
    Double I_KINETIC_I4x2z_Dxz_vrr = PBZ*I_KINETIC_I4x2z_Px_vrr+2*oned2z*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dxz_vrr;
    Double I_KINETIC_I3x3y_Dxz_vrr = PBZ*I_KINETIC_I3x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dxz_vrr;
    Double I_KINETIC_I3x2yz_Dxz_vrr = PBZ*I_KINETIC_I3x2yz_Px_vrr+oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dxz_vrr;
    Double I_KINETIC_I3x3z_Dxz_vrr = PBZ*I_KINETIC_I3x3z_Px_vrr+3*oned2z*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dxz_vrr;
    Double I_KINETIC_I2x4y_Dxz_vrr = PBZ*I_KINETIC_I2x4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dxz_vrr;
    Double I_KINETIC_I2x3yz_Dxz_vrr = PBZ*I_KINETIC_I2x3yz_Px_vrr+oned2z*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dxz_vrr;
    Double I_KINETIC_I2xy3z_Dxz_vrr = PBZ*I_KINETIC_I2xy3z_Px_vrr+3*oned2z*I_KINETIC_H2xy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dxz_vrr;
    Double I_KINETIC_I2x4z_Dxz_vrr = PBZ*I_KINETIC_I2x4z_Px_vrr+4*oned2z*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dxz_vrr;
    Double I_KINETIC_Ix5y_Dxz_vrr = PBZ*I_KINETIC_Ix5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dxz_vrr;
    Double I_KINETIC_Ix5z_Dxz_vrr = PBZ*I_KINETIC_Ix5z_Px_vrr+5*oned2z*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dxz_vrr;
    Double I_KINETIC_I6y_Dxz_vrr = PBZ*I_KINETIC_I6y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dxz_vrr;
    Double I_KINETIC_I5yz_Dxz_vrr = PBZ*I_KINETIC_I5yz_Px_vrr+oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dxz_vrr;
    Double I_KINETIC_I4y2z_Dxz_vrr = PBZ*I_KINETIC_I4y2z_Px_vrr+2*oned2z*I_KINETIC_H4yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dxz_vrr;
    Double I_KINETIC_I3y3z_Dxz_vrr = PBZ*I_KINETIC_I3y3z_Px_vrr+3*oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dxz_vrr;
    Double I_KINETIC_I2y4z_Dxz_vrr = PBZ*I_KINETIC_I2y4z_Px_vrr+4*oned2z*I_KINETIC_H2y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dxz_vrr;
    Double I_KINETIC_Iy5z_Dxz_vrr = PBZ*I_KINETIC_Iy5z_Px_vrr+5*oned2z*I_KINETIC_Hy4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dxz_vrr;
    Double I_KINETIC_I6z_Dxz_vrr = PBZ*I_KINETIC_I6z_Px_vrr+6*oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dxz_vrr;
    Double I_KINETIC_I6x_D2y_vrr = PBY*I_KINETIC_I6x_Py_vrr+oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2y_vrr-adz*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_KINETIC_I5xy_D2y_vrr = PBY*I_KINETIC_I5xy_Py_vrr+oned2z*I_KINETIC_H5x_Py_vrr+oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2y_vrr-adz*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_D2y_vrr = PBY*I_KINETIC_I5xz_Py_vrr+oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2y_vrr-adz*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_D2y_vrr = PBY*I_KINETIC_I4x2y_Py_vrr+2*oned2z*I_KINETIC_H4xy_Py_vrr+oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_KINETIC_I4x2z_D2y_vrr = PBY*I_KINETIC_I4x2z_Py_vrr+oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_KINETIC_I3x3y_D2y_vrr = PBY*I_KINETIC_I3x3y_Py_vrr+3*oned2z*I_KINETIC_H3x2y_Py_vrr+oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2y_vrr-adz*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_KINETIC_I3x2yz_D2y_vrr = PBY*I_KINETIC_I3x2yz_Py_vrr+2*oned2z*I_KINETIC_H3xyz_Py_vrr+oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3x3z_D2y_vrr = PBY*I_KINETIC_I3x3z_Py_vrr+oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2y_vrr-adz*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_KINETIC_I2x4y_D2y_vrr = PBY*I_KINETIC_I2x4y_Py_vrr+4*oned2z*I_KINETIC_H2x3y_Py_vrr+oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2y_vrr-adz*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_KINETIC_I2x3yz_D2y_vrr = PBY*I_KINETIC_I2x3yz_Py_vrr+3*oned2z*I_KINETIC_H2x2yz_Py_vrr+oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2xy3z_D2y_vrr = PBY*I_KINETIC_I2xy3z_Py_vrr+oned2z*I_KINETIC_H2x3z_Py_vrr+oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_D2y_vrr = PBY*I_KINETIC_I2x4z_Py_vrr+oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2y_vrr-adz*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_KINETIC_Ix5y_D2y_vrr = PBY*I_KINETIC_Ix5y_Py_vrr+5*oned2z*I_KINETIC_Hx4y_Py_vrr+oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2y_vrr-adz*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix5z_D2y_vrr = PBY*I_KINETIC_Ix5z_Py_vrr+oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2y_vrr-adz*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_D2y_vrr = PBY*I_KINETIC_I6y_Py_vrr+6*oned2z*I_KINETIC_H5y_Py_vrr+oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2y_vrr-adz*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_KINETIC_I5yz_D2y_vrr = PBY*I_KINETIC_I5yz_Py_vrr+5*oned2z*I_KINETIC_H4yz_Py_vrr+oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2y_vrr-adz*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_D2y_vrr = PBY*I_KINETIC_I4y2z_Py_vrr+4*oned2z*I_KINETIC_H3y2z_Py_vrr+oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_KINETIC_I3y3z_D2y_vrr = PBY*I_KINETIC_I3y3z_Py_vrr+3*oned2z*I_KINETIC_H2y3z_Py_vrr+oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2y_vrr-adz*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_KINETIC_I2y4z_D2y_vrr = PBY*I_KINETIC_I2y4z_Py_vrr+2*oned2z*I_KINETIC_Hy4z_Py_vrr+oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2y_vrr-adz*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_KINETIC_Iy5z_D2y_vrr = PBY*I_KINETIC_Iy5z_Py_vrr+oned2z*I_KINETIC_H5z_Py_vrr+oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2y_vrr-adz*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_D2y_vrr = PBY*I_KINETIC_I6z_Py_vrr+oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2y_vrr-adz*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_KINETIC_I6x_Dyz_vrr = PBZ*I_KINETIC_I6x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dyz_vrr;
    Double I_KINETIC_I5xy_Dyz_vrr = PBZ*I_KINETIC_I5xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dyz_vrr;
    Double I_KINETIC_I5xz_Dyz_vrr = PBZ*I_KINETIC_I5xz_Py_vrr+oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dyz_vrr;
    Double I_KINETIC_I4x2y_Dyz_vrr = PBZ*I_KINETIC_I4x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dyz_vrr;
    Double I_KINETIC_I4x2z_Dyz_vrr = PBZ*I_KINETIC_I4x2z_Py_vrr+2*oned2z*I_KINETIC_H4xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dyz_vrr;
    Double I_KINETIC_I3x3y_Dyz_vrr = PBZ*I_KINETIC_I3x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dyz_vrr;
    Double I_KINETIC_I3x2yz_Dyz_vrr = PBZ*I_KINETIC_I3x2yz_Py_vrr+oned2z*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dyz_vrr;
    Double I_KINETIC_I3x3z_Dyz_vrr = PBZ*I_KINETIC_I3x3z_Py_vrr+3*oned2z*I_KINETIC_H3x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dyz_vrr;
    Double I_KINETIC_I2x4y_Dyz_vrr = PBZ*I_KINETIC_I2x4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dyz_vrr;
    Double I_KINETIC_I2x3yz_Dyz_vrr = PBZ*I_KINETIC_I2x3yz_Py_vrr+oned2z*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dyz_vrr;
    Double I_KINETIC_I2xy3z_Dyz_vrr = PBZ*I_KINETIC_I2xy3z_Py_vrr+3*oned2z*I_KINETIC_H2xy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dyz_vrr;
    Double I_KINETIC_I2x4z_Dyz_vrr = PBZ*I_KINETIC_I2x4z_Py_vrr+4*oned2z*I_KINETIC_H2x3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dyz_vrr;
    Double I_KINETIC_Ix5y_Dyz_vrr = PBZ*I_KINETIC_Ix5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dyz_vrr;
    Double I_KINETIC_Ix5z_Dyz_vrr = PBZ*I_KINETIC_Ix5z_Py_vrr+5*oned2z*I_KINETIC_Hx4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dyz_vrr;
    Double I_KINETIC_I6y_Dyz_vrr = PBZ*I_KINETIC_I6y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dyz_vrr;
    Double I_KINETIC_I5yz_Dyz_vrr = PBZ*I_KINETIC_I5yz_Py_vrr+oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dyz_vrr;
    Double I_KINETIC_I4y2z_Dyz_vrr = PBZ*I_KINETIC_I4y2z_Py_vrr+2*oned2z*I_KINETIC_H4yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dyz_vrr;
    Double I_KINETIC_I3y3z_Dyz_vrr = PBZ*I_KINETIC_I3y3z_Py_vrr+3*oned2z*I_KINETIC_H3y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dyz_vrr;
    Double I_KINETIC_I2y4z_Dyz_vrr = PBZ*I_KINETIC_I2y4z_Py_vrr+4*oned2z*I_KINETIC_H2y3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dyz_vrr;
    Double I_KINETIC_Iy5z_Dyz_vrr = PBZ*I_KINETIC_Iy5z_Py_vrr+5*oned2z*I_KINETIC_Hy4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dyz_vrr;
    Double I_KINETIC_I6z_Dyz_vrr = PBZ*I_KINETIC_I6z_Py_vrr+6*oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dyz_vrr;
    Double I_KINETIC_I6x_D2z_vrr = PBZ*I_KINETIC_I6x_Pz_vrr+oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2z_vrr-adz*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_KINETIC_I5xy_D2z_vrr = PBZ*I_KINETIC_I5xy_Pz_vrr+oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2z_vrr-adz*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_D2z_vrr = PBZ*I_KINETIC_I5xz_Pz_vrr+oned2z*I_KINETIC_H5x_Pz_vrr+oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2z_vrr-adz*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_D2z_vrr = PBZ*I_KINETIC_I4x2y_Pz_vrr+oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_KINETIC_I4x2z_D2z_vrr = PBZ*I_KINETIC_I4x2z_Pz_vrr+2*oned2z*I_KINETIC_H4xz_Pz_vrr+oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_KINETIC_I3x3y_D2z_vrr = PBZ*I_KINETIC_I3x3y_Pz_vrr+oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2z_vrr-adz*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_KINETIC_I3x2yz_D2z_vrr = PBZ*I_KINETIC_I3x2yz_Pz_vrr+oned2z*I_KINETIC_H3x2y_Pz_vrr+oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3x3z_D2z_vrr = PBZ*I_KINETIC_I3x3z_Pz_vrr+3*oned2z*I_KINETIC_H3x2z_Pz_vrr+oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2z_vrr-adz*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_KINETIC_I2x4y_D2z_vrr = PBZ*I_KINETIC_I2x4y_Pz_vrr+oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2z_vrr-adz*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_KINETIC_I2x3yz_D2z_vrr = PBZ*I_KINETIC_I2x3yz_Pz_vrr+oned2z*I_KINETIC_H2x3y_Pz_vrr+oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2xy3z_D2z_vrr = PBZ*I_KINETIC_I2xy3z_Pz_vrr+3*oned2z*I_KINETIC_H2xy2z_Pz_vrr+oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_D2z_vrr = PBZ*I_KINETIC_I2x4z_Pz_vrr+4*oned2z*I_KINETIC_H2x3z_Pz_vrr+oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2z_vrr-adz*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_KINETIC_Ix5y_D2z_vrr = PBZ*I_KINETIC_Ix5y_Pz_vrr+oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2z_vrr-adz*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix5z_D2z_vrr = PBZ*I_KINETIC_Ix5z_Pz_vrr+5*oned2z*I_KINETIC_Hx4z_Pz_vrr+oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2z_vrr-adz*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_D2z_vrr = PBZ*I_KINETIC_I6y_Pz_vrr+oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2z_vrr-adz*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_KINETIC_I5yz_D2z_vrr = PBZ*I_KINETIC_I5yz_Pz_vrr+oned2z*I_KINETIC_H5y_Pz_vrr+oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2z_vrr-adz*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_D2z_vrr = PBZ*I_KINETIC_I4y2z_Pz_vrr+2*oned2z*I_KINETIC_H4yz_Pz_vrr+oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_KINETIC_I3y3z_D2z_vrr = PBZ*I_KINETIC_I3y3z_Pz_vrr+3*oned2z*I_KINETIC_H3y2z_Pz_vrr+oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2z_vrr-adz*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_KINETIC_I2y4z_D2z_vrr = PBZ*I_KINETIC_I2y4z_Pz_vrr+4*oned2z*I_KINETIC_H2y3z_Pz_vrr+oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2z_vrr-adz*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_KINETIC_Iy5z_D2z_vrr = PBZ*I_KINETIC_Iy5z_Pz_vrr+5*oned2z*I_KINETIC_Hy4z_Pz_vrr+oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2z_vrr-adz*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_D2z_vrr = PBZ*I_KINETIC_I6z_Pz_vrr+6*oned2z*I_KINETIC_H5z_Pz_vrr+oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2z_vrr-adz*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_K_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_I_D
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_KINETIC_I_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     ************************************************************/
    Double I_KINETIC_K7x_D2x_vrr = PAX*I_KINETIC_I6x_D2x_vrr+6*oned2z*I_KINETIC_H5x_D2x_vrr+2*oned2z*I_KINETIC_I6x_Px_vrr+twoxi*I_TWOBODYOVERLAP_K7x_D2x_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_K6xy_D2x_vrr = PAY*I_KINETIC_I6x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_D2x_vrr;
    Double I_KINETIC_K6xz_D2x_vrr = PAZ*I_KINETIC_I6x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_D2x_vrr;
    Double I_KINETIC_K5x2y_D2x_vrr = PAY*I_KINETIC_I5xy_D2x_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_K5xyz_D2x_vrr = PAZ*I_KINETIC_I5xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_D2x_vrr;
    Double I_KINETIC_K5x2z_D2x_vrr = PAZ*I_KINETIC_I5xz_D2x_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_K4x3y_D2x_vrr = PAY*I_KINETIC_I4x2y_D2x_vrr+2*oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_K4x2yz_D2x_vrr = PAZ*I_KINETIC_I4x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_D2x_vrr;
    Double I_KINETIC_K4xy2z_D2x_vrr = PAY*I_KINETIC_I4x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_D2x_vrr;
    Double I_KINETIC_K4x3z_D2x_vrr = PAZ*I_KINETIC_I4x2z_D2x_vrr+2*oned2z*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_K3x4y_D2x_vrr = PAX*I_KINETIC_I2x4y_D2x_vrr+2*oned2z*I_KINETIC_Hx4y_D2x_vrr+2*oned2z*I_KINETIC_I2x4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_K3x3yz_D2x_vrr = PAZ*I_KINETIC_I3x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_D2x_vrr;
    Double I_KINETIC_K3x2y2z_D2x_vrr = PAZ*I_KINETIC_I3x2yz_D2x_vrr+oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_KINETIC_K3xy3z_D2x_vrr = PAY*I_KINETIC_I3x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_D2x_vrr;
    Double I_KINETIC_K3x4z_D2x_vrr = PAX*I_KINETIC_I2x4z_D2x_vrr+2*oned2z*I_KINETIC_Hx4z_D2x_vrr+2*oned2z*I_KINETIC_I2x4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_K2x5y_D2x_vrr = PAX*I_KINETIC_Ix5y_D2x_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+2*oned2z*I_KINETIC_Ix5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_K2x4yz_D2x_vrr = PAZ*I_KINETIC_I2x4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_D2x_vrr;
    Double I_KINETIC_K2x3y2z_D2x_vrr = PAZ*I_KINETIC_I2x3yz_D2x_vrr+oned2z*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_KINETIC_K2x2y3z_D2x_vrr = PAY*I_KINETIC_I2xy3z_D2x_vrr+oned2z*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_KINETIC_K2xy4z_D2x_vrr = PAY*I_KINETIC_I2x4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_D2x_vrr;
    Double I_KINETIC_K2x5z_D2x_vrr = PAX*I_KINETIC_Ix5z_D2x_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+2*oned2z*I_KINETIC_Ix5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_Kx6y_D2x_vrr = PAX*I_KINETIC_I6y_D2x_vrr+2*oned2z*I_KINETIC_I6y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_D2x_vrr;
    Double I_KINETIC_Kx5yz_D2x_vrr = PAZ*I_KINETIC_Ix5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_D2x_vrr;
    Double I_KINETIC_Kx4y2z_D2x_vrr = PAX*I_KINETIC_I4y2z_D2x_vrr+2*oned2z*I_KINETIC_I4y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_D2x_vrr;
    Double I_KINETIC_Kx3y3z_D2x_vrr = PAX*I_KINETIC_I3y3z_D2x_vrr+2*oned2z*I_KINETIC_I3y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_D2x_vrr;
    Double I_KINETIC_Kx2y4z_D2x_vrr = PAX*I_KINETIC_I2y4z_D2x_vrr+2*oned2z*I_KINETIC_I2y4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_D2x_vrr;
    Double I_KINETIC_Kxy5z_D2x_vrr = PAY*I_KINETIC_Ix5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_D2x_vrr;
    Double I_KINETIC_Kx6z_D2x_vrr = PAX*I_KINETIC_I6z_D2x_vrr+2*oned2z*I_KINETIC_I6z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_D2x_vrr;
    Double I_KINETIC_K7y_D2x_vrr = PAY*I_KINETIC_I6y_D2x_vrr+6*oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K7y_D2x_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_K6yz_D2x_vrr = PAZ*I_KINETIC_I6y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_D2x_vrr;
    Double I_KINETIC_K5y2z_D2x_vrr = PAZ*I_KINETIC_I5yz_D2x_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_K4y3z_D2x_vrr = PAZ*I_KINETIC_I4y2z_D2x_vrr+2*oned2z*I_KINETIC_H4yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_K3y4z_D2x_vrr = PAY*I_KINETIC_I2y4z_D2x_vrr+2*oned2z*I_KINETIC_Hy4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_K2y5z_D2x_vrr = PAY*I_KINETIC_Iy5z_D2x_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_D2x_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_Ky6z_D2x_vrr = PAY*I_KINETIC_I6z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_D2x_vrr;
    Double I_KINETIC_K7z_D2x_vrr = PAZ*I_KINETIC_I6z_D2x_vrr+6*oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_K7z_D2x_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_K7x_Dxy_vrr = PAX*I_KINETIC_I6x_Dxy_vrr+6*oned2z*I_KINETIC_H5x_Dxy_vrr+oned2z*I_KINETIC_I6x_Py_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Dxy_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_KINETIC_K6xy_Dxy_vrr = PAY*I_KINETIC_I6x_Dxy_vrr+oned2z*I_KINETIC_I6x_Px_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Dxy_vrr;
    Double I_KINETIC_K6xz_Dxy_vrr = PAZ*I_KINETIC_I6x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Dxy_vrr;
    Double I_KINETIC_K5x2y_Dxy_vrr = PAY*I_KINETIC_I5xy_Dxy_vrr+oned2z*I_KINETIC_H5x_Dxy_vrr+oned2z*I_KINETIC_I5xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_KINETIC_K5xyz_Dxy_vrr = PAZ*I_KINETIC_I5xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Dxy_vrr;
    Double I_KINETIC_K5x2z_Dxy_vrr = PAZ*I_KINETIC_I5xz_Dxy_vrr+oned2z*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_KINETIC_K4x3y_Dxy_vrr = PAY*I_KINETIC_I4x2y_Dxy_vrr+2*oned2z*I_KINETIC_H4xy_Dxy_vrr+oned2z*I_KINETIC_I4x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_KINETIC_K4x2yz_Dxy_vrr = PAZ*I_KINETIC_I4x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Dxy_vrr;
    Double I_KINETIC_K4xy2z_Dxy_vrr = PAY*I_KINETIC_I4x2z_Dxy_vrr+oned2z*I_KINETIC_I4x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Dxy_vrr;
    Double I_KINETIC_K4x3z_Dxy_vrr = PAZ*I_KINETIC_I4x2z_Dxy_vrr+2*oned2z*I_KINETIC_H4xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_KINETIC_K3x4y_Dxy_vrr = PAX*I_KINETIC_I2x4y_Dxy_vrr+2*oned2z*I_KINETIC_Hx4y_Dxy_vrr+oned2z*I_KINETIC_I2x4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_KINETIC_K3x3yz_Dxy_vrr = PAZ*I_KINETIC_I3x3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Dxy_vrr;
    Double I_KINETIC_K3x2y2z_Dxy_vrr = PAZ*I_KINETIC_I3x2yz_Dxy_vrr+oned2z*I_KINETIC_H3x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_KINETIC_K3xy3z_Dxy_vrr = PAY*I_KINETIC_I3x3z_Dxy_vrr+oned2z*I_KINETIC_I3x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Dxy_vrr;
    Double I_KINETIC_K3x4z_Dxy_vrr = PAX*I_KINETIC_I2x4z_Dxy_vrr+2*oned2z*I_KINETIC_Hx4z_Dxy_vrr+oned2z*I_KINETIC_I2x4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_KINETIC_K2x5y_Dxy_vrr = PAX*I_KINETIC_Ix5y_Dxy_vrr+oned2z*I_KINETIC_H5y_Dxy_vrr+oned2z*I_KINETIC_Ix5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_KINETIC_K2x4yz_Dxy_vrr = PAZ*I_KINETIC_I2x4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Dxy_vrr;
    Double I_KINETIC_K2x3y2z_Dxy_vrr = PAZ*I_KINETIC_I2x3yz_Dxy_vrr+oned2z*I_KINETIC_H2x3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_KINETIC_K2x2y3z_Dxy_vrr = PAY*I_KINETIC_I2xy3z_Dxy_vrr+oned2z*I_KINETIC_H2x3z_Dxy_vrr+oned2z*I_KINETIC_I2xy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr;
    Double I_KINETIC_K2xy4z_Dxy_vrr = PAY*I_KINETIC_I2x4z_Dxy_vrr+oned2z*I_KINETIC_I2x4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Dxy_vrr;
    Double I_KINETIC_K2x5z_Dxy_vrr = PAX*I_KINETIC_Ix5z_Dxy_vrr+oned2z*I_KINETIC_H5z_Dxy_vrr+oned2z*I_KINETIC_Ix5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_KINETIC_Kx6y_Dxy_vrr = PAX*I_KINETIC_I6y_Dxy_vrr+oned2z*I_KINETIC_I6y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Dxy_vrr;
    Double I_KINETIC_Kx5yz_Dxy_vrr = PAZ*I_KINETIC_Ix5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Dxy_vrr;
    Double I_KINETIC_Kx4y2z_Dxy_vrr = PAX*I_KINETIC_I4y2z_Dxy_vrr+oned2z*I_KINETIC_I4y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Dxy_vrr;
    Double I_KINETIC_Kx3y3z_Dxy_vrr = PAX*I_KINETIC_I3y3z_Dxy_vrr+oned2z*I_KINETIC_I3y3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Dxy_vrr;
    Double I_KINETIC_Kx2y4z_Dxy_vrr = PAX*I_KINETIC_I2y4z_Dxy_vrr+oned2z*I_KINETIC_I2y4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Dxy_vrr;
    Double I_KINETIC_Kxy5z_Dxy_vrr = PAY*I_KINETIC_Ix5z_Dxy_vrr+oned2z*I_KINETIC_Ix5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Dxy_vrr;
    Double I_KINETIC_Kx6z_Dxy_vrr = PAX*I_KINETIC_I6z_Dxy_vrr+oned2z*I_KINETIC_I6z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Dxy_vrr;
    Double I_KINETIC_K7y_Dxy_vrr = PAY*I_KINETIC_I6y_Dxy_vrr+6*oned2z*I_KINETIC_H5y_Dxy_vrr+oned2z*I_KINETIC_I6y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Dxy_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_KINETIC_K6yz_Dxy_vrr = PAZ*I_KINETIC_I6y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Dxy_vrr;
    Double I_KINETIC_K5y2z_Dxy_vrr = PAZ*I_KINETIC_I5yz_Dxy_vrr+oned2z*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_KINETIC_K4y3z_Dxy_vrr = PAZ*I_KINETIC_I4y2z_Dxy_vrr+2*oned2z*I_KINETIC_H4yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_KINETIC_K3y4z_Dxy_vrr = PAY*I_KINETIC_I2y4z_Dxy_vrr+2*oned2z*I_KINETIC_Hy4z_Dxy_vrr+oned2z*I_KINETIC_I2y4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_KINETIC_K2y5z_Dxy_vrr = PAY*I_KINETIC_Iy5z_Dxy_vrr+oned2z*I_KINETIC_H5z_Dxy_vrr+oned2z*I_KINETIC_Iy5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_KINETIC_Ky6z_Dxy_vrr = PAY*I_KINETIC_I6z_Dxy_vrr+oned2z*I_KINETIC_I6z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Dxy_vrr;
    Double I_KINETIC_K7z_Dxy_vrr = PAZ*I_KINETIC_I6z_Dxy_vrr+6*oned2z*I_KINETIC_H5z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Dxy_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_KINETIC_K7x_Dxz_vrr = PAX*I_KINETIC_I6x_Dxz_vrr+6*oned2z*I_KINETIC_H5x_Dxz_vrr+oned2z*I_KINETIC_I6x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Dxz_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_KINETIC_K6xy_Dxz_vrr = PAY*I_KINETIC_I6x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Dxz_vrr;
    Double I_KINETIC_K6xz_Dxz_vrr = PAZ*I_KINETIC_I6x_Dxz_vrr+oned2z*I_KINETIC_I6x_Px_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Dxz_vrr;
    Double I_KINETIC_K5x2y_Dxz_vrr = PAY*I_KINETIC_I5xy_Dxz_vrr+oned2z*I_KINETIC_H5x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_KINETIC_K5xyz_Dxz_vrr = PAZ*I_KINETIC_I5xy_Dxz_vrr+oned2z*I_KINETIC_I5xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Dxz_vrr;
    Double I_KINETIC_K5x2z_Dxz_vrr = PAZ*I_KINETIC_I5xz_Dxz_vrr+oned2z*I_KINETIC_H5x_Dxz_vrr+oned2z*I_KINETIC_I5xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_KINETIC_K4x3y_Dxz_vrr = PAY*I_KINETIC_I4x2y_Dxz_vrr+2*oned2z*I_KINETIC_H4xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_KINETIC_K4x2yz_Dxz_vrr = PAZ*I_KINETIC_I4x2y_Dxz_vrr+oned2z*I_KINETIC_I4x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Dxz_vrr;
    Double I_KINETIC_K4xy2z_Dxz_vrr = PAY*I_KINETIC_I4x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Dxz_vrr;
    Double I_KINETIC_K4x3z_Dxz_vrr = PAZ*I_KINETIC_I4x2z_Dxz_vrr+2*oned2z*I_KINETIC_H4xz_Dxz_vrr+oned2z*I_KINETIC_I4x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_KINETIC_K3x4y_Dxz_vrr = PAX*I_KINETIC_I2x4y_Dxz_vrr+2*oned2z*I_KINETIC_Hx4y_Dxz_vrr+oned2z*I_KINETIC_I2x4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_KINETIC_K3x3yz_Dxz_vrr = PAZ*I_KINETIC_I3x3y_Dxz_vrr+oned2z*I_KINETIC_I3x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Dxz_vrr;
    Double I_KINETIC_K3x2y2z_Dxz_vrr = PAZ*I_KINETIC_I3x2yz_Dxz_vrr+oned2z*I_KINETIC_H3x2y_Dxz_vrr+oned2z*I_KINETIC_I3x2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr;
    Double I_KINETIC_K3xy3z_Dxz_vrr = PAY*I_KINETIC_I3x3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Dxz_vrr;
    Double I_KINETIC_K3x4z_Dxz_vrr = PAX*I_KINETIC_I2x4z_Dxz_vrr+2*oned2z*I_KINETIC_Hx4z_Dxz_vrr+oned2z*I_KINETIC_I2x4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_KINETIC_K2x5y_Dxz_vrr = PAX*I_KINETIC_Ix5y_Dxz_vrr+oned2z*I_KINETIC_H5y_Dxz_vrr+oned2z*I_KINETIC_Ix5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_KINETIC_K2x4yz_Dxz_vrr = PAZ*I_KINETIC_I2x4y_Dxz_vrr+oned2z*I_KINETIC_I2x4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Dxz_vrr;
    Double I_KINETIC_K2x3y2z_Dxz_vrr = PAZ*I_KINETIC_I2x3yz_Dxz_vrr+oned2z*I_KINETIC_H2x3y_Dxz_vrr+oned2z*I_KINETIC_I2x3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr;
    Double I_KINETIC_K2x2y3z_Dxz_vrr = PAY*I_KINETIC_I2xy3z_Dxz_vrr+oned2z*I_KINETIC_H2x3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_KINETIC_K2xy4z_Dxz_vrr = PAY*I_KINETIC_I2x4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Dxz_vrr;
    Double I_KINETIC_K2x5z_Dxz_vrr = PAX*I_KINETIC_Ix5z_Dxz_vrr+oned2z*I_KINETIC_H5z_Dxz_vrr+oned2z*I_KINETIC_Ix5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_KINETIC_Kx6y_Dxz_vrr = PAX*I_KINETIC_I6y_Dxz_vrr+oned2z*I_KINETIC_I6y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Dxz_vrr;
    Double I_KINETIC_Kx5yz_Dxz_vrr = PAZ*I_KINETIC_Ix5y_Dxz_vrr+oned2z*I_KINETIC_Ix5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Dxz_vrr;
    Double I_KINETIC_Kx4y2z_Dxz_vrr = PAX*I_KINETIC_I4y2z_Dxz_vrr+oned2z*I_KINETIC_I4y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Dxz_vrr;
    Double I_KINETIC_Kx3y3z_Dxz_vrr = PAX*I_KINETIC_I3y3z_Dxz_vrr+oned2z*I_KINETIC_I3y3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Dxz_vrr;
    Double I_KINETIC_Kx2y4z_Dxz_vrr = PAX*I_KINETIC_I2y4z_Dxz_vrr+oned2z*I_KINETIC_I2y4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Dxz_vrr;
    Double I_KINETIC_Kxy5z_Dxz_vrr = PAY*I_KINETIC_Ix5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Dxz_vrr;
    Double I_KINETIC_Kx6z_Dxz_vrr = PAX*I_KINETIC_I6z_Dxz_vrr+oned2z*I_KINETIC_I6z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Dxz_vrr;
    Double I_KINETIC_K7y_Dxz_vrr = PAY*I_KINETIC_I6y_Dxz_vrr+6*oned2z*I_KINETIC_H5y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Dxz_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_KINETIC_K6yz_Dxz_vrr = PAZ*I_KINETIC_I6y_Dxz_vrr+oned2z*I_KINETIC_I6y_Px_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Dxz_vrr;
    Double I_KINETIC_K5y2z_Dxz_vrr = PAZ*I_KINETIC_I5yz_Dxz_vrr+oned2z*I_KINETIC_H5y_Dxz_vrr+oned2z*I_KINETIC_I5yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_KINETIC_K4y3z_Dxz_vrr = PAZ*I_KINETIC_I4y2z_Dxz_vrr+2*oned2z*I_KINETIC_H4yz_Dxz_vrr+oned2z*I_KINETIC_I4y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_KINETIC_K3y4z_Dxz_vrr = PAY*I_KINETIC_I2y4z_Dxz_vrr+2*oned2z*I_KINETIC_Hy4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_KINETIC_K2y5z_Dxz_vrr = PAY*I_KINETIC_Iy5z_Dxz_vrr+oned2z*I_KINETIC_H5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_KINETIC_Ky6z_Dxz_vrr = PAY*I_KINETIC_I6z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Dxz_vrr;
    Double I_KINETIC_K7z_Dxz_vrr = PAZ*I_KINETIC_I6z_Dxz_vrr+6*oned2z*I_KINETIC_H5z_Dxz_vrr+oned2z*I_KINETIC_I6z_Px_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Dxz_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_KINETIC_K7x_D2y_vrr = PAX*I_KINETIC_I6x_D2y_vrr+6*oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K7x_D2y_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_KINETIC_K6xy_D2y_vrr = PAY*I_KINETIC_I6x_D2y_vrr+2*oned2z*I_KINETIC_I6x_Py_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_D2y_vrr;
    Double I_KINETIC_K6xz_D2y_vrr = PAZ*I_KINETIC_I6x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_D2y_vrr;
    Double I_KINETIC_K5x2y_D2y_vrr = PAY*I_KINETIC_I5xy_D2y_vrr+oned2z*I_KINETIC_H5x_D2y_vrr+2*oned2z*I_KINETIC_I5xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_KINETIC_K5xyz_D2y_vrr = PAZ*I_KINETIC_I5xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_D2y_vrr;
    Double I_KINETIC_K5x2z_D2y_vrr = PAZ*I_KINETIC_I5xz_D2y_vrr+oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_KINETIC_K4x3y_D2y_vrr = PAY*I_KINETIC_I4x2y_D2y_vrr+2*oned2z*I_KINETIC_H4xy_D2y_vrr+2*oned2z*I_KINETIC_I4x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_KINETIC_K4x2yz_D2y_vrr = PAZ*I_KINETIC_I4x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_D2y_vrr;
    Double I_KINETIC_K4xy2z_D2y_vrr = PAY*I_KINETIC_I4x2z_D2y_vrr+2*oned2z*I_KINETIC_I4x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_D2y_vrr;
    Double I_KINETIC_K4x3z_D2y_vrr = PAZ*I_KINETIC_I4x2z_D2y_vrr+2*oned2z*I_KINETIC_H4xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_KINETIC_K3x4y_D2y_vrr = PAX*I_KINETIC_I2x4y_D2y_vrr+2*oned2z*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_KINETIC_K3x3yz_D2y_vrr = PAZ*I_KINETIC_I3x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_D2y_vrr;
    Double I_KINETIC_K3x2y2z_D2y_vrr = PAZ*I_KINETIC_I3x2yz_D2y_vrr+oned2z*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_KINETIC_K3xy3z_D2y_vrr = PAY*I_KINETIC_I3x3z_D2y_vrr+2*oned2z*I_KINETIC_I3x3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_D2y_vrr;
    Double I_KINETIC_K3x4z_D2y_vrr = PAX*I_KINETIC_I2x4z_D2y_vrr+2*oned2z*I_KINETIC_Hx4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_KINETIC_K2x5y_D2y_vrr = PAX*I_KINETIC_Ix5y_D2y_vrr+oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_KINETIC_K2x4yz_D2y_vrr = PAZ*I_KINETIC_I2x4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_D2y_vrr;
    Double I_KINETIC_K2x3y2z_D2y_vrr = PAZ*I_KINETIC_I2x3yz_D2y_vrr+oned2z*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_KINETIC_K2x2y3z_D2y_vrr = PAY*I_KINETIC_I2xy3z_D2y_vrr+oned2z*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_I2xy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_KINETIC_K2xy4z_D2y_vrr = PAY*I_KINETIC_I2x4z_D2y_vrr+2*oned2z*I_KINETIC_I2x4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_D2y_vrr;
    Double I_KINETIC_K2x5z_D2y_vrr = PAX*I_KINETIC_Ix5z_D2y_vrr+oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_KINETIC_Kx6y_D2y_vrr = PAX*I_KINETIC_I6y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_D2y_vrr;
    Double I_KINETIC_Kx5yz_D2y_vrr = PAZ*I_KINETIC_Ix5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_D2y_vrr;
    Double I_KINETIC_Kx4y2z_D2y_vrr = PAX*I_KINETIC_I4y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_D2y_vrr;
    Double I_KINETIC_Kx3y3z_D2y_vrr = PAX*I_KINETIC_I3y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_D2y_vrr;
    Double I_KINETIC_Kx2y4z_D2y_vrr = PAX*I_KINETIC_I2y4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_D2y_vrr;
    Double I_KINETIC_Kxy5z_D2y_vrr = PAY*I_KINETIC_Ix5z_D2y_vrr+2*oned2z*I_KINETIC_Ix5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_D2y_vrr;
    Double I_KINETIC_Kx6z_D2y_vrr = PAX*I_KINETIC_I6z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_D2y_vrr;
    Double I_KINETIC_K7y_D2y_vrr = PAY*I_KINETIC_I6y_D2y_vrr+6*oned2z*I_KINETIC_H5y_D2y_vrr+2*oned2z*I_KINETIC_I6y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K7y_D2y_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_KINETIC_K6yz_D2y_vrr = PAZ*I_KINETIC_I6y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_D2y_vrr;
    Double I_KINETIC_K5y2z_D2y_vrr = PAZ*I_KINETIC_I5yz_D2y_vrr+oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_KINETIC_K4y3z_D2y_vrr = PAZ*I_KINETIC_I4y2z_D2y_vrr+2*oned2z*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_KINETIC_K3y4z_D2y_vrr = PAY*I_KINETIC_I2y4z_D2y_vrr+2*oned2z*I_KINETIC_Hy4z_D2y_vrr+2*oned2z*I_KINETIC_I2y4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_KINETIC_K2y5z_D2y_vrr = PAY*I_KINETIC_Iy5z_D2y_vrr+oned2z*I_KINETIC_H5z_D2y_vrr+2*oned2z*I_KINETIC_Iy5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_D2y_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_KINETIC_Ky6z_D2y_vrr = PAY*I_KINETIC_I6z_D2y_vrr+2*oned2z*I_KINETIC_I6z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_D2y_vrr;
    Double I_KINETIC_K7z_D2y_vrr = PAZ*I_KINETIC_I6z_D2y_vrr+6*oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_K7z_D2y_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_KINETIC_K7x_Dyz_vrr = PAX*I_KINETIC_I6x_Dyz_vrr+6*oned2z*I_KINETIC_H5x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Dyz_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_KINETIC_K6xy_Dyz_vrr = PAY*I_KINETIC_I6x_Dyz_vrr+oned2z*I_KINETIC_I6x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Dyz_vrr;
    Double I_KINETIC_K6xz_Dyz_vrr = PAZ*I_KINETIC_I6x_Dyz_vrr+oned2z*I_KINETIC_I6x_Py_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Dyz_vrr;
    Double I_KINETIC_K5x2y_Dyz_vrr = PAY*I_KINETIC_I5xy_Dyz_vrr+oned2z*I_KINETIC_H5x_Dyz_vrr+oned2z*I_KINETIC_I5xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_KINETIC_K5xyz_Dyz_vrr = PAZ*I_KINETIC_I5xy_Dyz_vrr+oned2z*I_KINETIC_I5xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Dyz_vrr;
    Double I_KINETIC_K5x2z_Dyz_vrr = PAZ*I_KINETIC_I5xz_Dyz_vrr+oned2z*I_KINETIC_H5x_Dyz_vrr+oned2z*I_KINETIC_I5xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_KINETIC_K4x3y_Dyz_vrr = PAY*I_KINETIC_I4x2y_Dyz_vrr+2*oned2z*I_KINETIC_H4xy_Dyz_vrr+oned2z*I_KINETIC_I4x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_KINETIC_K4x2yz_Dyz_vrr = PAZ*I_KINETIC_I4x2y_Dyz_vrr+oned2z*I_KINETIC_I4x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Dyz_vrr;
    Double I_KINETIC_K4xy2z_Dyz_vrr = PAY*I_KINETIC_I4x2z_Dyz_vrr+oned2z*I_KINETIC_I4x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Dyz_vrr;
    Double I_KINETIC_K4x3z_Dyz_vrr = PAZ*I_KINETIC_I4x2z_Dyz_vrr+2*oned2z*I_KINETIC_H4xz_Dyz_vrr+oned2z*I_KINETIC_I4x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_KINETIC_K3x4y_Dyz_vrr = PAX*I_KINETIC_I2x4y_Dyz_vrr+2*oned2z*I_KINETIC_Hx4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_KINETIC_K3x3yz_Dyz_vrr = PAZ*I_KINETIC_I3x3y_Dyz_vrr+oned2z*I_KINETIC_I3x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Dyz_vrr;
    Double I_KINETIC_K3x2y2z_Dyz_vrr = PAZ*I_KINETIC_I3x2yz_Dyz_vrr+oned2z*I_KINETIC_H3x2y_Dyz_vrr+oned2z*I_KINETIC_I3x2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr;
    Double I_KINETIC_K3xy3z_Dyz_vrr = PAY*I_KINETIC_I3x3z_Dyz_vrr+oned2z*I_KINETIC_I3x3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Dyz_vrr;
    Double I_KINETIC_K3x4z_Dyz_vrr = PAX*I_KINETIC_I2x4z_Dyz_vrr+2*oned2z*I_KINETIC_Hx4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_KINETIC_K2x5y_Dyz_vrr = PAX*I_KINETIC_Ix5y_Dyz_vrr+oned2z*I_KINETIC_H5y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_KINETIC_K2x4yz_Dyz_vrr = PAZ*I_KINETIC_I2x4y_Dyz_vrr+oned2z*I_KINETIC_I2x4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Dyz_vrr;
    Double I_KINETIC_K2x3y2z_Dyz_vrr = PAZ*I_KINETIC_I2x3yz_Dyz_vrr+oned2z*I_KINETIC_H2x3y_Dyz_vrr+oned2z*I_KINETIC_I2x3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr;
    Double I_KINETIC_K2x2y3z_Dyz_vrr = PAY*I_KINETIC_I2xy3z_Dyz_vrr+oned2z*I_KINETIC_H2x3z_Dyz_vrr+oned2z*I_KINETIC_I2xy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr;
    Double I_KINETIC_K2xy4z_Dyz_vrr = PAY*I_KINETIC_I2x4z_Dyz_vrr+oned2z*I_KINETIC_I2x4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Dyz_vrr;
    Double I_KINETIC_K2x5z_Dyz_vrr = PAX*I_KINETIC_Ix5z_Dyz_vrr+oned2z*I_KINETIC_H5z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_KINETIC_Kx6y_Dyz_vrr = PAX*I_KINETIC_I6y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Dyz_vrr;
    Double I_KINETIC_Kx5yz_Dyz_vrr = PAZ*I_KINETIC_Ix5y_Dyz_vrr+oned2z*I_KINETIC_Ix5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Dyz_vrr;
    Double I_KINETIC_Kx4y2z_Dyz_vrr = PAX*I_KINETIC_I4y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Dyz_vrr;
    Double I_KINETIC_Kx3y3z_Dyz_vrr = PAX*I_KINETIC_I3y3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Dyz_vrr;
    Double I_KINETIC_Kx2y4z_Dyz_vrr = PAX*I_KINETIC_I2y4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Dyz_vrr;
    Double I_KINETIC_Kxy5z_Dyz_vrr = PAY*I_KINETIC_Ix5z_Dyz_vrr+oned2z*I_KINETIC_Ix5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Dyz_vrr;
    Double I_KINETIC_Kx6z_Dyz_vrr = PAX*I_KINETIC_I6z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Dyz_vrr;
    Double I_KINETIC_K7y_Dyz_vrr = PAY*I_KINETIC_I6y_Dyz_vrr+6*oned2z*I_KINETIC_H5y_Dyz_vrr+oned2z*I_KINETIC_I6y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Dyz_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_KINETIC_K6yz_Dyz_vrr = PAZ*I_KINETIC_I6y_Dyz_vrr+oned2z*I_KINETIC_I6y_Py_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Dyz_vrr;
    Double I_KINETIC_K5y2z_Dyz_vrr = PAZ*I_KINETIC_I5yz_Dyz_vrr+oned2z*I_KINETIC_H5y_Dyz_vrr+oned2z*I_KINETIC_I5yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_KINETIC_K4y3z_Dyz_vrr = PAZ*I_KINETIC_I4y2z_Dyz_vrr+2*oned2z*I_KINETIC_H4yz_Dyz_vrr+oned2z*I_KINETIC_I4y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_KINETIC_K3y4z_Dyz_vrr = PAY*I_KINETIC_I2y4z_Dyz_vrr+2*oned2z*I_KINETIC_Hy4z_Dyz_vrr+oned2z*I_KINETIC_I2y4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_KINETIC_K2y5z_Dyz_vrr = PAY*I_KINETIC_Iy5z_Dyz_vrr+oned2z*I_KINETIC_H5z_Dyz_vrr+oned2z*I_KINETIC_Iy5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_KINETIC_Ky6z_Dyz_vrr = PAY*I_KINETIC_I6z_Dyz_vrr+oned2z*I_KINETIC_I6z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Dyz_vrr;
    Double I_KINETIC_K7z_Dyz_vrr = PAZ*I_KINETIC_I6z_Dyz_vrr+6*oned2z*I_KINETIC_H5z_Dyz_vrr+oned2z*I_KINETIC_I6z_Py_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Dyz_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_KINETIC_K7x_D2z_vrr = PAX*I_KINETIC_I6x_D2z_vrr+6*oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K7x_D2z_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_KINETIC_K6xy_D2z_vrr = PAY*I_KINETIC_I6x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_D2z_vrr;
    Double I_KINETIC_K6xz_D2z_vrr = PAZ*I_KINETIC_I6x_D2z_vrr+2*oned2z*I_KINETIC_I6x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_D2z_vrr;
    Double I_KINETIC_K5x2y_D2z_vrr = PAY*I_KINETIC_I5xy_D2z_vrr+oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_KINETIC_K5xyz_D2z_vrr = PAZ*I_KINETIC_I5xy_D2z_vrr+2*oned2z*I_KINETIC_I5xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_D2z_vrr;
    Double I_KINETIC_K5x2z_D2z_vrr = PAZ*I_KINETIC_I5xz_D2z_vrr+oned2z*I_KINETIC_H5x_D2z_vrr+2*oned2z*I_KINETIC_I5xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_KINETIC_K4x3y_D2z_vrr = PAY*I_KINETIC_I4x2y_D2z_vrr+2*oned2z*I_KINETIC_H4xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_KINETIC_K4x2yz_D2z_vrr = PAZ*I_KINETIC_I4x2y_D2z_vrr+2*oned2z*I_KINETIC_I4x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_D2z_vrr;
    Double I_KINETIC_K4xy2z_D2z_vrr = PAY*I_KINETIC_I4x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_D2z_vrr;
    Double I_KINETIC_K4x3z_D2z_vrr = PAZ*I_KINETIC_I4x2z_D2z_vrr+2*oned2z*I_KINETIC_H4xz_D2z_vrr+2*oned2z*I_KINETIC_I4x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_KINETIC_K3x4y_D2z_vrr = PAX*I_KINETIC_I2x4y_D2z_vrr+2*oned2z*I_KINETIC_Hx4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_KINETIC_K3x3yz_D2z_vrr = PAZ*I_KINETIC_I3x3y_D2z_vrr+2*oned2z*I_KINETIC_I3x3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_D2z_vrr;
    Double I_KINETIC_K3x2y2z_D2z_vrr = PAZ*I_KINETIC_I3x2yz_D2z_vrr+oned2z*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_I3x2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_D2z_vrr;
    Double I_KINETIC_K3xy3z_D2z_vrr = PAY*I_KINETIC_I3x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_D2z_vrr;
    Double I_KINETIC_K3x4z_D2z_vrr = PAX*I_KINETIC_I2x4z_D2z_vrr+2*oned2z*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_KINETIC_K2x5y_D2z_vrr = PAX*I_KINETIC_Ix5y_D2z_vrr+oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_KINETIC_K2x4yz_D2z_vrr = PAZ*I_KINETIC_I2x4y_D2z_vrr+2*oned2z*I_KINETIC_I2x4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_D2z_vrr;
    Double I_KINETIC_K2x3y2z_D2z_vrr = PAZ*I_KINETIC_I2x3yz_D2z_vrr+oned2z*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_I2x3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_D2z_vrr;
    Double I_KINETIC_K2x2y3z_D2z_vrr = PAY*I_KINETIC_I2xy3z_D2z_vrr+oned2z*I_KINETIC_H2x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_KINETIC_K2xy4z_D2z_vrr = PAY*I_KINETIC_I2x4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_D2z_vrr;
    Double I_KINETIC_K2x5z_D2z_vrr = PAX*I_KINETIC_Ix5z_D2z_vrr+oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_KINETIC_Kx6y_D2z_vrr = PAX*I_KINETIC_I6y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_D2z_vrr;
    Double I_KINETIC_Kx5yz_D2z_vrr = PAZ*I_KINETIC_Ix5y_D2z_vrr+2*oned2z*I_KINETIC_Ix5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_D2z_vrr;
    Double I_KINETIC_Kx4y2z_D2z_vrr = PAX*I_KINETIC_I4y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_D2z_vrr;
    Double I_KINETIC_Kx3y3z_D2z_vrr = PAX*I_KINETIC_I3y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_D2z_vrr;
    Double I_KINETIC_Kx2y4z_D2z_vrr = PAX*I_KINETIC_I2y4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_D2z_vrr;
    Double I_KINETIC_Kxy5z_D2z_vrr = PAY*I_KINETIC_Ix5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_D2z_vrr;
    Double I_KINETIC_Kx6z_D2z_vrr = PAX*I_KINETIC_I6z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_D2z_vrr;
    Double I_KINETIC_K7y_D2z_vrr = PAY*I_KINETIC_I6y_D2z_vrr+6*oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K7y_D2z_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_KINETIC_K6yz_D2z_vrr = PAZ*I_KINETIC_I6y_D2z_vrr+2*oned2z*I_KINETIC_I6y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_D2z_vrr;
    Double I_KINETIC_K5y2z_D2z_vrr = PAZ*I_KINETIC_I5yz_D2z_vrr+oned2z*I_KINETIC_H5y_D2z_vrr+2*oned2z*I_KINETIC_I5yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_KINETIC_K4y3z_D2z_vrr = PAZ*I_KINETIC_I4y2z_D2z_vrr+2*oned2z*I_KINETIC_H4yz_D2z_vrr+2*oned2z*I_KINETIC_I4y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_KINETIC_K3y4z_D2z_vrr = PAY*I_KINETIC_I2y4z_D2z_vrr+2*oned2z*I_KINETIC_Hy4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_KINETIC_K2y5z_D2z_vrr = PAY*I_KINETIC_Iy5z_D2z_vrr+oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_D2z_vrr-bdz*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_KINETIC_Ky6z_D2z_vrr = PAY*I_KINETIC_I6z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_D2z_vrr;
    Double I_KINETIC_K7z_D2z_vrr = PAZ*I_KINETIC_I6z_D2z_vrr+6*oned2z*I_KINETIC_H5z_D2z_vrr+2*oned2z*I_KINETIC_I6z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_K7z_D2z_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_K_D_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_K_D_aa_coefs = alpha*alpha;
    I_KINETIC_K7x_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_D2x_vrr;
    I_KINETIC_K6xy_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_D2x_vrr;
    I_KINETIC_K6xz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_D2x_vrr;
    I_KINETIC_K5x2y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_D2x_vrr;
    I_KINETIC_K5xyz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_D2x_vrr;
    I_KINETIC_K5x2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_D2x_vrr;
    I_KINETIC_K4x3y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_D2x_vrr;
    I_KINETIC_K4x2yz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_D2x_vrr;
    I_KINETIC_K4xy2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_D2x_vrr;
    I_KINETIC_K4x3z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_D2x_vrr;
    I_KINETIC_K3x4y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_D2x_vrr;
    I_KINETIC_K3x3yz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_D2x_vrr;
    I_KINETIC_K3x2y2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_D2x_vrr;
    I_KINETIC_K3xy3z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_D2x_vrr;
    I_KINETIC_K3x4z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_D2x_vrr;
    I_KINETIC_K2x5y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_D2x_vrr;
    I_KINETIC_K2x4yz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_D2x_vrr;
    I_KINETIC_K2x3y2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_D2x_vrr;
    I_KINETIC_K2x2y3z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_D2x_vrr;
    I_KINETIC_K2xy4z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_D2x_vrr;
    I_KINETIC_K2x5z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_D2x_vrr;
    I_KINETIC_Kx6y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_D2x_vrr;
    I_KINETIC_Kx5yz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_D2x_vrr;
    I_KINETIC_Kx4y2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_D2x_vrr;
    I_KINETIC_Kx3y3z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_D2x_vrr;
    I_KINETIC_Kx2y4z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_D2x_vrr;
    I_KINETIC_Kxy5z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_D2x_vrr;
    I_KINETIC_Kx6z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_D2x_vrr;
    I_KINETIC_K7y_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_D2x_vrr;
    I_KINETIC_K6yz_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_D2x_vrr;
    I_KINETIC_K5y2z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_D2x_vrr;
    I_KINETIC_K4y3z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_D2x_vrr;
    I_KINETIC_K3y4z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_D2x_vrr;
    I_KINETIC_K2y5z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_D2x_vrr;
    I_KINETIC_Ky6z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_D2x_vrr;
    I_KINETIC_K7z_D2x_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_D2x_vrr;
    I_KINETIC_K7x_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_Dxy_vrr;
    I_KINETIC_K6xy_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_Dxy_vrr;
    I_KINETIC_K6xz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_Dxy_vrr;
    I_KINETIC_K5x2y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_Dxy_vrr;
    I_KINETIC_K5xyz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_Dxy_vrr;
    I_KINETIC_K5x2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_Dxy_vrr;
    I_KINETIC_K4x3y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_Dxy_vrr;
    I_KINETIC_K4x2yz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_Dxy_vrr;
    I_KINETIC_K4xy2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_Dxy_vrr;
    I_KINETIC_K4x3z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_Dxy_vrr;
    I_KINETIC_K3x4y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_Dxy_vrr;
    I_KINETIC_K3x3yz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_Dxy_vrr;
    I_KINETIC_K3x2y2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_Dxy_vrr;
    I_KINETIC_K3xy3z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_Dxy_vrr;
    I_KINETIC_K3x4z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_Dxy_vrr;
    I_KINETIC_K2x5y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_Dxy_vrr;
    I_KINETIC_K2x4yz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_Dxy_vrr;
    I_KINETIC_K2x3y2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_Dxy_vrr;
    I_KINETIC_K2x2y3z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_Dxy_vrr;
    I_KINETIC_K2xy4z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_Dxy_vrr;
    I_KINETIC_K2x5z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_Dxy_vrr;
    I_KINETIC_Kx6y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_Dxy_vrr;
    I_KINETIC_Kx5yz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_Dxy_vrr;
    I_KINETIC_Kx4y2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_Dxy_vrr;
    I_KINETIC_Kx3y3z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_Dxy_vrr;
    I_KINETIC_Kx2y4z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_Dxy_vrr;
    I_KINETIC_Kxy5z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_Dxy_vrr;
    I_KINETIC_Kx6z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_Dxy_vrr;
    I_KINETIC_K7y_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_Dxy_vrr;
    I_KINETIC_K6yz_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_Dxy_vrr;
    I_KINETIC_K5y2z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_Dxy_vrr;
    I_KINETIC_K4y3z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_Dxy_vrr;
    I_KINETIC_K3y4z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_Dxy_vrr;
    I_KINETIC_K2y5z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_Dxy_vrr;
    I_KINETIC_Ky6z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_Dxy_vrr;
    I_KINETIC_K7z_Dxy_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_Dxy_vrr;
    I_KINETIC_K7x_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_Dxz_vrr;
    I_KINETIC_K6xy_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_Dxz_vrr;
    I_KINETIC_K6xz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_Dxz_vrr;
    I_KINETIC_K5x2y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_Dxz_vrr;
    I_KINETIC_K5xyz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_Dxz_vrr;
    I_KINETIC_K5x2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_Dxz_vrr;
    I_KINETIC_K4x3y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_Dxz_vrr;
    I_KINETIC_K4x2yz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_Dxz_vrr;
    I_KINETIC_K4xy2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_Dxz_vrr;
    I_KINETIC_K4x3z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_Dxz_vrr;
    I_KINETIC_K3x4y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_Dxz_vrr;
    I_KINETIC_K3x3yz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_Dxz_vrr;
    I_KINETIC_K3x2y2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_Dxz_vrr;
    I_KINETIC_K3xy3z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_Dxz_vrr;
    I_KINETIC_K3x4z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_Dxz_vrr;
    I_KINETIC_K2x5y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_Dxz_vrr;
    I_KINETIC_K2x4yz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_Dxz_vrr;
    I_KINETIC_K2x3y2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_Dxz_vrr;
    I_KINETIC_K2x2y3z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_Dxz_vrr;
    I_KINETIC_K2xy4z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_Dxz_vrr;
    I_KINETIC_K2x5z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_Dxz_vrr;
    I_KINETIC_Kx6y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_Dxz_vrr;
    I_KINETIC_Kx5yz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_Dxz_vrr;
    I_KINETIC_Kx4y2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_Dxz_vrr;
    I_KINETIC_Kx3y3z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_Dxz_vrr;
    I_KINETIC_Kx2y4z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_Dxz_vrr;
    I_KINETIC_Kxy5z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_Dxz_vrr;
    I_KINETIC_Kx6z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_Dxz_vrr;
    I_KINETIC_K7y_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_Dxz_vrr;
    I_KINETIC_K6yz_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_Dxz_vrr;
    I_KINETIC_K5y2z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_Dxz_vrr;
    I_KINETIC_K4y3z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_Dxz_vrr;
    I_KINETIC_K3y4z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_Dxz_vrr;
    I_KINETIC_K2y5z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_Dxz_vrr;
    I_KINETIC_Ky6z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_Dxz_vrr;
    I_KINETIC_K7z_Dxz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_Dxz_vrr;
    I_KINETIC_K7x_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_D2y_vrr;
    I_KINETIC_K6xy_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_D2y_vrr;
    I_KINETIC_K6xz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_D2y_vrr;
    I_KINETIC_K5x2y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_D2y_vrr;
    I_KINETIC_K5xyz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_D2y_vrr;
    I_KINETIC_K5x2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_D2y_vrr;
    I_KINETIC_K4x3y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_D2y_vrr;
    I_KINETIC_K4x2yz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_D2y_vrr;
    I_KINETIC_K4xy2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_D2y_vrr;
    I_KINETIC_K4x3z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_D2y_vrr;
    I_KINETIC_K3x4y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_D2y_vrr;
    I_KINETIC_K3x3yz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_D2y_vrr;
    I_KINETIC_K3x2y2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_D2y_vrr;
    I_KINETIC_K3xy3z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_D2y_vrr;
    I_KINETIC_K3x4z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_D2y_vrr;
    I_KINETIC_K2x5y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_D2y_vrr;
    I_KINETIC_K2x4yz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_D2y_vrr;
    I_KINETIC_K2x3y2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_D2y_vrr;
    I_KINETIC_K2x2y3z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_D2y_vrr;
    I_KINETIC_K2xy4z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_D2y_vrr;
    I_KINETIC_K2x5z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_D2y_vrr;
    I_KINETIC_Kx6y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_D2y_vrr;
    I_KINETIC_Kx5yz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_D2y_vrr;
    I_KINETIC_Kx4y2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_D2y_vrr;
    I_KINETIC_Kx3y3z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_D2y_vrr;
    I_KINETIC_Kx2y4z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_D2y_vrr;
    I_KINETIC_Kxy5z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_D2y_vrr;
    I_KINETIC_Kx6z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_D2y_vrr;
    I_KINETIC_K7y_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_D2y_vrr;
    I_KINETIC_K6yz_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_D2y_vrr;
    I_KINETIC_K5y2z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_D2y_vrr;
    I_KINETIC_K4y3z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_D2y_vrr;
    I_KINETIC_K3y4z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_D2y_vrr;
    I_KINETIC_K2y5z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_D2y_vrr;
    I_KINETIC_Ky6z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_D2y_vrr;
    I_KINETIC_K7z_D2y_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_D2y_vrr;
    I_KINETIC_K7x_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_Dyz_vrr;
    I_KINETIC_K6xy_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_Dyz_vrr;
    I_KINETIC_K6xz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_Dyz_vrr;
    I_KINETIC_K5x2y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_Dyz_vrr;
    I_KINETIC_K5xyz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_Dyz_vrr;
    I_KINETIC_K5x2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_Dyz_vrr;
    I_KINETIC_K4x3y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_Dyz_vrr;
    I_KINETIC_K4x2yz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_Dyz_vrr;
    I_KINETIC_K4xy2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_Dyz_vrr;
    I_KINETIC_K4x3z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_Dyz_vrr;
    I_KINETIC_K3x4y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_Dyz_vrr;
    I_KINETIC_K3x3yz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_Dyz_vrr;
    I_KINETIC_K3x2y2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_Dyz_vrr;
    I_KINETIC_K3xy3z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_Dyz_vrr;
    I_KINETIC_K3x4z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_Dyz_vrr;
    I_KINETIC_K2x5y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_Dyz_vrr;
    I_KINETIC_K2x4yz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_Dyz_vrr;
    I_KINETIC_K2x3y2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_Dyz_vrr;
    I_KINETIC_K2x2y3z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_Dyz_vrr;
    I_KINETIC_K2xy4z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_Dyz_vrr;
    I_KINETIC_K2x5z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_Dyz_vrr;
    I_KINETIC_Kx6y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_Dyz_vrr;
    I_KINETIC_Kx5yz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_Dyz_vrr;
    I_KINETIC_Kx4y2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_Dyz_vrr;
    I_KINETIC_Kx3y3z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_Dyz_vrr;
    I_KINETIC_Kx2y4z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_Dyz_vrr;
    I_KINETIC_Kxy5z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_Dyz_vrr;
    I_KINETIC_Kx6z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_Dyz_vrr;
    I_KINETIC_K7y_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_Dyz_vrr;
    I_KINETIC_K6yz_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_Dyz_vrr;
    I_KINETIC_K5y2z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_Dyz_vrr;
    I_KINETIC_K4y3z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_Dyz_vrr;
    I_KINETIC_K3y4z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_Dyz_vrr;
    I_KINETIC_K2y5z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_Dyz_vrr;
    I_KINETIC_Ky6z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_Dyz_vrr;
    I_KINETIC_K7z_Dyz_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_Dyz_vrr;
    I_KINETIC_K7x_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7x_D2z_vrr;
    I_KINETIC_K6xy_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xy_D2z_vrr;
    I_KINETIC_K6xz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6xz_D2z_vrr;
    I_KINETIC_K5x2y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2y_D2z_vrr;
    I_KINETIC_K5xyz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5xyz_D2z_vrr;
    I_KINETIC_K5x2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5x2z_D2z_vrr;
    I_KINETIC_K4x3y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3y_D2z_vrr;
    I_KINETIC_K4x2yz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x2yz_D2z_vrr;
    I_KINETIC_K4xy2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4xy2z_D2z_vrr;
    I_KINETIC_K4x3z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4x3z_D2z_vrr;
    I_KINETIC_K3x4y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4y_D2z_vrr;
    I_KINETIC_K3x3yz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x3yz_D2z_vrr;
    I_KINETIC_K3x2y2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x2y2z_D2z_vrr;
    I_KINETIC_K3xy3z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3xy3z_D2z_vrr;
    I_KINETIC_K3x4z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3x4z_D2z_vrr;
    I_KINETIC_K2x5y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5y_D2z_vrr;
    I_KINETIC_K2x4yz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x4yz_D2z_vrr;
    I_KINETIC_K2x3y2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x3y2z_D2z_vrr;
    I_KINETIC_K2x2y3z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x2y3z_D2z_vrr;
    I_KINETIC_K2xy4z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2xy4z_D2z_vrr;
    I_KINETIC_K2x5z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2x5z_D2z_vrr;
    I_KINETIC_Kx6y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6y_D2z_vrr;
    I_KINETIC_Kx5yz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx5yz_D2z_vrr;
    I_KINETIC_Kx4y2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx4y2z_D2z_vrr;
    I_KINETIC_Kx3y3z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx3y3z_D2z_vrr;
    I_KINETIC_Kx2y4z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx2y4z_D2z_vrr;
    I_KINETIC_Kxy5z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kxy5z_D2z_vrr;
    I_KINETIC_Kx6z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Kx6z_D2z_vrr;
    I_KINETIC_K7y_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7y_D2z_vrr;
    I_KINETIC_K6yz_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K6yz_D2z_vrr;
    I_KINETIC_K5y2z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K5y2z_D2z_vrr;
    I_KINETIC_K4y3z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K4y3z_D2z_vrr;
    I_KINETIC_K3y4z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K3y4z_D2z_vrr;
    I_KINETIC_K2y5z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K2y5z_D2z_vrr;
    I_KINETIC_Ky6z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_Ky6z_D2z_vrr;
    I_KINETIC_K7z_D2z_aa += SQ_KINETIC_K_D_aa_coefs*I_KINETIC_K7z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_D_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_D_a_coefs = alpha;
    I_KINETIC_H5x_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_D2x_vrr;
    I_KINETIC_H4xy_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_D2x_vrr;
    I_KINETIC_H4xz_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_D2x_vrr;
    I_KINETIC_H3x2y_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_D2x_vrr;
    I_KINETIC_H3xyz_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_D2x_vrr;
    I_KINETIC_H3x2z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_D2x_vrr;
    I_KINETIC_H2x3y_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_D2x_vrr;
    I_KINETIC_H2x2yz_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_D2x_vrr;
    I_KINETIC_H2xy2z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_D2x_vrr;
    I_KINETIC_H2x3z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_D2x_vrr;
    I_KINETIC_Hx4y_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_D2x_vrr;
    I_KINETIC_Hx3yz_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_D2x_vrr;
    I_KINETIC_Hx2y2z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_D2x_vrr;
    I_KINETIC_Hxy3z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_D2x_vrr;
    I_KINETIC_Hx4z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_D2x_vrr;
    I_KINETIC_H5y_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_D2x_vrr;
    I_KINETIC_H4yz_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_D2x_vrr;
    I_KINETIC_H3y2z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_D2x_vrr;
    I_KINETIC_H2y3z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_D2x_vrr;
    I_KINETIC_Hy4z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_D2x_vrr;
    I_KINETIC_H5z_D2x_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_D2x_vrr;
    I_KINETIC_H5x_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_Dxy_vrr;
    I_KINETIC_H4xy_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_Dxy_vrr;
    I_KINETIC_H4xz_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_Dxy_vrr;
    I_KINETIC_H3x2y_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_Dxy_vrr;
    I_KINETIC_H3xyz_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_Dxy_vrr;
    I_KINETIC_H3x2z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_Dxy_vrr;
    I_KINETIC_H2x3y_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_Dxy_vrr;
    I_KINETIC_H2x2yz_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_Dxy_vrr;
    I_KINETIC_H2xy2z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_Dxy_vrr;
    I_KINETIC_H2x3z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_Dxy_vrr;
    I_KINETIC_Hx4y_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_Dxy_vrr;
    I_KINETIC_Hx3yz_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_Dxy_vrr;
    I_KINETIC_Hx2y2z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_Dxy_vrr;
    I_KINETIC_Hxy3z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_Dxy_vrr;
    I_KINETIC_Hx4z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_Dxy_vrr;
    I_KINETIC_H5y_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_Dxy_vrr;
    I_KINETIC_H4yz_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_Dxy_vrr;
    I_KINETIC_H3y2z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_Dxy_vrr;
    I_KINETIC_H2y3z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_Dxy_vrr;
    I_KINETIC_Hy4z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_Dxy_vrr;
    I_KINETIC_H5z_Dxy_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_Dxy_vrr;
    I_KINETIC_H5x_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_Dxz_vrr;
    I_KINETIC_H4xy_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_Dxz_vrr;
    I_KINETIC_H4xz_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_Dxz_vrr;
    I_KINETIC_H3x2y_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_Dxz_vrr;
    I_KINETIC_H3xyz_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_Dxz_vrr;
    I_KINETIC_H3x2z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_Dxz_vrr;
    I_KINETIC_H2x3y_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_Dxz_vrr;
    I_KINETIC_H2x2yz_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_Dxz_vrr;
    I_KINETIC_H2xy2z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_Dxz_vrr;
    I_KINETIC_H2x3z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_Dxz_vrr;
    I_KINETIC_Hx4y_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_Dxz_vrr;
    I_KINETIC_Hx3yz_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_Dxz_vrr;
    I_KINETIC_Hx2y2z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_Dxz_vrr;
    I_KINETIC_Hxy3z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_Dxz_vrr;
    I_KINETIC_Hx4z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_Dxz_vrr;
    I_KINETIC_H5y_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_Dxz_vrr;
    I_KINETIC_H4yz_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_Dxz_vrr;
    I_KINETIC_H3y2z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_Dxz_vrr;
    I_KINETIC_H2y3z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_Dxz_vrr;
    I_KINETIC_Hy4z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_Dxz_vrr;
    I_KINETIC_H5z_Dxz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_Dxz_vrr;
    I_KINETIC_H5x_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_D2y_vrr;
    I_KINETIC_H4xy_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_D2y_vrr;
    I_KINETIC_H4xz_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_D2y_vrr;
    I_KINETIC_H3x2y_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_D2y_vrr;
    I_KINETIC_H3xyz_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_D2y_vrr;
    I_KINETIC_H3x2z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_D2y_vrr;
    I_KINETIC_H2x3y_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_D2y_vrr;
    I_KINETIC_H2x2yz_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_D2y_vrr;
    I_KINETIC_H2xy2z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_D2y_vrr;
    I_KINETIC_H2x3z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_D2y_vrr;
    I_KINETIC_Hx4y_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_D2y_vrr;
    I_KINETIC_Hx3yz_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_D2y_vrr;
    I_KINETIC_Hx2y2z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_D2y_vrr;
    I_KINETIC_Hxy3z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_D2y_vrr;
    I_KINETIC_Hx4z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_D2y_vrr;
    I_KINETIC_H5y_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_D2y_vrr;
    I_KINETIC_H4yz_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_D2y_vrr;
    I_KINETIC_H3y2z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_D2y_vrr;
    I_KINETIC_H2y3z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_D2y_vrr;
    I_KINETIC_Hy4z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_D2y_vrr;
    I_KINETIC_H5z_D2y_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_D2y_vrr;
    I_KINETIC_H5x_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_Dyz_vrr;
    I_KINETIC_H4xy_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_Dyz_vrr;
    I_KINETIC_H4xz_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_Dyz_vrr;
    I_KINETIC_H3x2y_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_Dyz_vrr;
    I_KINETIC_H3xyz_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_Dyz_vrr;
    I_KINETIC_H3x2z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_Dyz_vrr;
    I_KINETIC_H2x3y_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_Dyz_vrr;
    I_KINETIC_H2x2yz_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_Dyz_vrr;
    I_KINETIC_H2xy2z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_Dyz_vrr;
    I_KINETIC_H2x3z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_Dyz_vrr;
    I_KINETIC_Hx4y_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_Dyz_vrr;
    I_KINETIC_Hx3yz_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_Dyz_vrr;
    I_KINETIC_Hx2y2z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_Dyz_vrr;
    I_KINETIC_Hxy3z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_Dyz_vrr;
    I_KINETIC_Hx4z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_Dyz_vrr;
    I_KINETIC_H5y_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_Dyz_vrr;
    I_KINETIC_H4yz_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_Dyz_vrr;
    I_KINETIC_H3y2z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_Dyz_vrr;
    I_KINETIC_H2y3z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_Dyz_vrr;
    I_KINETIC_Hy4z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_Dyz_vrr;
    I_KINETIC_H5z_Dyz_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_Dyz_vrr;
    I_KINETIC_H5x_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5x_D2z_vrr;
    I_KINETIC_H4xy_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xy_D2z_vrr;
    I_KINETIC_H4xz_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4xz_D2z_vrr;
    I_KINETIC_H3x2y_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2y_D2z_vrr;
    I_KINETIC_H3xyz_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3xyz_D2z_vrr;
    I_KINETIC_H3x2z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3x2z_D2z_vrr;
    I_KINETIC_H2x3y_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3y_D2z_vrr;
    I_KINETIC_H2x2yz_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x2yz_D2z_vrr;
    I_KINETIC_H2xy2z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2xy2z_D2z_vrr;
    I_KINETIC_H2x3z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2x3z_D2z_vrr;
    I_KINETIC_Hx4y_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4y_D2z_vrr;
    I_KINETIC_Hx3yz_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx3yz_D2z_vrr;
    I_KINETIC_Hx2y2z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx2y2z_D2z_vrr;
    I_KINETIC_Hxy3z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hxy3z_D2z_vrr;
    I_KINETIC_Hx4z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hx4z_D2z_vrr;
    I_KINETIC_H5y_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5y_D2z_vrr;
    I_KINETIC_H4yz_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H4yz_D2z_vrr;
    I_KINETIC_H3y2z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H3y2z_D2z_vrr;
    I_KINETIC_H2y3z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H2y3z_D2z_vrr;
    I_KINETIC_Hy4z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_Hy4z_D2z_vrr;
    I_KINETIC_H5z_D2z_a += SQ_KINETIC_H_D_a_coefs*I_KINETIC_H5z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_F3x_D2x += I_KINETIC_F3x_D2x_vrr;
    I_KINETIC_F2xy_D2x += I_KINETIC_F2xy_D2x_vrr;
    I_KINETIC_F2xz_D2x += I_KINETIC_F2xz_D2x_vrr;
    I_KINETIC_Fx2y_D2x += I_KINETIC_Fx2y_D2x_vrr;
    I_KINETIC_Fxyz_D2x += I_KINETIC_Fxyz_D2x_vrr;
    I_KINETIC_Fx2z_D2x += I_KINETIC_Fx2z_D2x_vrr;
    I_KINETIC_F3y_D2x += I_KINETIC_F3y_D2x_vrr;
    I_KINETIC_F2yz_D2x += I_KINETIC_F2yz_D2x_vrr;
    I_KINETIC_Fy2z_D2x += I_KINETIC_Fy2z_D2x_vrr;
    I_KINETIC_F3z_D2x += I_KINETIC_F3z_D2x_vrr;
    I_KINETIC_F3x_Dxy += I_KINETIC_F3x_Dxy_vrr;
    I_KINETIC_F2xy_Dxy += I_KINETIC_F2xy_Dxy_vrr;
    I_KINETIC_F2xz_Dxy += I_KINETIC_F2xz_Dxy_vrr;
    I_KINETIC_Fx2y_Dxy += I_KINETIC_Fx2y_Dxy_vrr;
    I_KINETIC_Fxyz_Dxy += I_KINETIC_Fxyz_Dxy_vrr;
    I_KINETIC_Fx2z_Dxy += I_KINETIC_Fx2z_Dxy_vrr;
    I_KINETIC_F3y_Dxy += I_KINETIC_F3y_Dxy_vrr;
    I_KINETIC_F2yz_Dxy += I_KINETIC_F2yz_Dxy_vrr;
    I_KINETIC_Fy2z_Dxy += I_KINETIC_Fy2z_Dxy_vrr;
    I_KINETIC_F3z_Dxy += I_KINETIC_F3z_Dxy_vrr;
    I_KINETIC_F3x_Dxz += I_KINETIC_F3x_Dxz_vrr;
    I_KINETIC_F2xy_Dxz += I_KINETIC_F2xy_Dxz_vrr;
    I_KINETIC_F2xz_Dxz += I_KINETIC_F2xz_Dxz_vrr;
    I_KINETIC_Fx2y_Dxz += I_KINETIC_Fx2y_Dxz_vrr;
    I_KINETIC_Fxyz_Dxz += I_KINETIC_Fxyz_Dxz_vrr;
    I_KINETIC_Fx2z_Dxz += I_KINETIC_Fx2z_Dxz_vrr;
    I_KINETIC_F3y_Dxz += I_KINETIC_F3y_Dxz_vrr;
    I_KINETIC_F2yz_Dxz += I_KINETIC_F2yz_Dxz_vrr;
    I_KINETIC_Fy2z_Dxz += I_KINETIC_Fy2z_Dxz_vrr;
    I_KINETIC_F3z_Dxz += I_KINETIC_F3z_Dxz_vrr;
    I_KINETIC_F3x_D2y += I_KINETIC_F3x_D2y_vrr;
    I_KINETIC_F2xy_D2y += I_KINETIC_F2xy_D2y_vrr;
    I_KINETIC_F2xz_D2y += I_KINETIC_F2xz_D2y_vrr;
    I_KINETIC_Fx2y_D2y += I_KINETIC_Fx2y_D2y_vrr;
    I_KINETIC_Fxyz_D2y += I_KINETIC_Fxyz_D2y_vrr;
    I_KINETIC_Fx2z_D2y += I_KINETIC_Fx2z_D2y_vrr;
    I_KINETIC_F3y_D2y += I_KINETIC_F3y_D2y_vrr;
    I_KINETIC_F2yz_D2y += I_KINETIC_F2yz_D2y_vrr;
    I_KINETIC_Fy2z_D2y += I_KINETIC_Fy2z_D2y_vrr;
    I_KINETIC_F3z_D2y += I_KINETIC_F3z_D2y_vrr;
    I_KINETIC_F3x_Dyz += I_KINETIC_F3x_Dyz_vrr;
    I_KINETIC_F2xy_Dyz += I_KINETIC_F2xy_Dyz_vrr;
    I_KINETIC_F2xz_Dyz += I_KINETIC_F2xz_Dyz_vrr;
    I_KINETIC_Fx2y_Dyz += I_KINETIC_Fx2y_Dyz_vrr;
    I_KINETIC_Fxyz_Dyz += I_KINETIC_Fxyz_Dyz_vrr;
    I_KINETIC_Fx2z_Dyz += I_KINETIC_Fx2z_Dyz_vrr;
    I_KINETIC_F3y_Dyz += I_KINETIC_F3y_Dyz_vrr;
    I_KINETIC_F2yz_Dyz += I_KINETIC_F2yz_Dyz_vrr;
    I_KINETIC_Fy2z_Dyz += I_KINETIC_Fy2z_Dyz_vrr;
    I_KINETIC_F3z_Dyz += I_KINETIC_F3z_Dyz_vrr;
    I_KINETIC_F3x_D2z += I_KINETIC_F3x_D2z_vrr;
    I_KINETIC_F2xy_D2z += I_KINETIC_F2xy_D2z_vrr;
    I_KINETIC_F2xz_D2z += I_KINETIC_F2xz_D2z_vrr;
    I_KINETIC_Fx2y_D2z += I_KINETIC_Fx2y_D2z_vrr;
    I_KINETIC_Fxyz_D2z += I_KINETIC_Fxyz_D2z_vrr;
    I_KINETIC_Fx2z_D2z += I_KINETIC_Fx2z_D2z_vrr;
    I_KINETIC_F3y_D2z += I_KINETIC_F3y_D2z_vrr;
    I_KINETIC_F2yz_D2z += I_KINETIC_F2yz_D2z_vrr;
    I_KINETIC_Fy2z_D2z += I_KINETIC_Fy2z_D2z_vrr;
    I_KINETIC_F3z_D2z += I_KINETIC_F3z_D2z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_K7x_D2x_aa-2.0E0*5*I_KINETIC_H5x_D2x_a-2.0E0*6*I_KINETIC_H5x_D2x_a+5*4*I_KINETIC_F3x_D2x;
  abcd[1] = 4.0E0*I_KINETIC_K6xy_D2x_aa-2.0E0*4*I_KINETIC_H4xy_D2x_a-2.0E0*5*I_KINETIC_H4xy_D2x_a+4*3*I_KINETIC_F2xy_D2x;
  abcd[2] = 4.0E0*I_KINETIC_K6xz_D2x_aa-2.0E0*4*I_KINETIC_H4xz_D2x_a-2.0E0*5*I_KINETIC_H4xz_D2x_a+4*3*I_KINETIC_F2xz_D2x;
  abcd[3] = 4.0E0*I_KINETIC_K5x2y_D2x_aa-2.0E0*3*I_KINETIC_H3x2y_D2x_a-2.0E0*4*I_KINETIC_H3x2y_D2x_a+3*2*I_KINETIC_Fx2y_D2x;
  abcd[4] = 4.0E0*I_KINETIC_K5xyz_D2x_aa-2.0E0*3*I_KINETIC_H3xyz_D2x_a-2.0E0*4*I_KINETIC_H3xyz_D2x_a+3*2*I_KINETIC_Fxyz_D2x;
  abcd[5] = 4.0E0*I_KINETIC_K5x2z_D2x_aa-2.0E0*3*I_KINETIC_H3x2z_D2x_a-2.0E0*4*I_KINETIC_H3x2z_D2x_a+3*2*I_KINETIC_Fx2z_D2x;
  abcd[6] = 4.0E0*I_KINETIC_K4x3y_D2x_aa-2.0E0*2*I_KINETIC_H2x3y_D2x_a-2.0E0*3*I_KINETIC_H2x3y_D2x_a+2*1*I_KINETIC_F3y_D2x;
  abcd[7] = 4.0E0*I_KINETIC_K4x2yz_D2x_aa-2.0E0*2*I_KINETIC_H2x2yz_D2x_a-2.0E0*3*I_KINETIC_H2x2yz_D2x_a+2*1*I_KINETIC_F2yz_D2x;
  abcd[8] = 4.0E0*I_KINETIC_K4xy2z_D2x_aa-2.0E0*2*I_KINETIC_H2xy2z_D2x_a-2.0E0*3*I_KINETIC_H2xy2z_D2x_a+2*1*I_KINETIC_Fy2z_D2x;
  abcd[9] = 4.0E0*I_KINETIC_K4x3z_D2x_aa-2.0E0*2*I_KINETIC_H2x3z_D2x_a-2.0E0*3*I_KINETIC_H2x3z_D2x_a+2*1*I_KINETIC_F3z_D2x;
  abcd[10] = 4.0E0*I_KINETIC_K3x4y_D2x_aa-2.0E0*1*I_KINETIC_Hx4y_D2x_a-2.0E0*2*I_KINETIC_Hx4y_D2x_a;
  abcd[11] = 4.0E0*I_KINETIC_K3x3yz_D2x_aa-2.0E0*1*I_KINETIC_Hx3yz_D2x_a-2.0E0*2*I_KINETIC_Hx3yz_D2x_a;
  abcd[12] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_Hx2y2z_D2x_a-2.0E0*2*I_KINETIC_Hx2y2z_D2x_a;
  abcd[13] = 4.0E0*I_KINETIC_K3xy3z_D2x_aa-2.0E0*1*I_KINETIC_Hxy3z_D2x_a-2.0E0*2*I_KINETIC_Hxy3z_D2x_a;
  abcd[14] = 4.0E0*I_KINETIC_K3x4z_D2x_aa-2.0E0*1*I_KINETIC_Hx4z_D2x_a-2.0E0*2*I_KINETIC_Hx4z_D2x_a;
  abcd[15] = 4.0E0*I_KINETIC_K2x5y_D2x_aa-2.0E0*1*I_KINETIC_H5y_D2x_a;
  abcd[16] = 4.0E0*I_KINETIC_K2x4yz_D2x_aa-2.0E0*1*I_KINETIC_H4yz_D2x_a;
  abcd[17] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*1*I_KINETIC_H3y2z_D2x_a;
  abcd[18] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*1*I_KINETIC_H2y3z_D2x_a;
  abcd[19] = 4.0E0*I_KINETIC_K2xy4z_D2x_aa-2.0E0*1*I_KINETIC_Hy4z_D2x_a;
  abcd[20] = 4.0E0*I_KINETIC_K2x5z_D2x_aa-2.0E0*1*I_KINETIC_H5z_D2x_a;
  abcd[21] = 4.0E0*I_KINETIC_K7x_Dxy_aa-2.0E0*5*I_KINETIC_H5x_Dxy_a-2.0E0*6*I_KINETIC_H5x_Dxy_a+5*4*I_KINETIC_F3x_Dxy;
  abcd[22] = 4.0E0*I_KINETIC_K6xy_Dxy_aa-2.0E0*4*I_KINETIC_H4xy_Dxy_a-2.0E0*5*I_KINETIC_H4xy_Dxy_a+4*3*I_KINETIC_F2xy_Dxy;
  abcd[23] = 4.0E0*I_KINETIC_K6xz_Dxy_aa-2.0E0*4*I_KINETIC_H4xz_Dxy_a-2.0E0*5*I_KINETIC_H4xz_Dxy_a+4*3*I_KINETIC_F2xz_Dxy;
  abcd[24] = 4.0E0*I_KINETIC_K5x2y_Dxy_aa-2.0E0*3*I_KINETIC_H3x2y_Dxy_a-2.0E0*4*I_KINETIC_H3x2y_Dxy_a+3*2*I_KINETIC_Fx2y_Dxy;
  abcd[25] = 4.0E0*I_KINETIC_K5xyz_Dxy_aa-2.0E0*3*I_KINETIC_H3xyz_Dxy_a-2.0E0*4*I_KINETIC_H3xyz_Dxy_a+3*2*I_KINETIC_Fxyz_Dxy;
  abcd[26] = 4.0E0*I_KINETIC_K5x2z_Dxy_aa-2.0E0*3*I_KINETIC_H3x2z_Dxy_a-2.0E0*4*I_KINETIC_H3x2z_Dxy_a+3*2*I_KINETIC_Fx2z_Dxy;
  abcd[27] = 4.0E0*I_KINETIC_K4x3y_Dxy_aa-2.0E0*2*I_KINETIC_H2x3y_Dxy_a-2.0E0*3*I_KINETIC_H2x3y_Dxy_a+2*1*I_KINETIC_F3y_Dxy;
  abcd[28] = 4.0E0*I_KINETIC_K4x2yz_Dxy_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxy_a-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a+2*1*I_KINETIC_F2yz_Dxy;
  abcd[29] = 4.0E0*I_KINETIC_K4xy2z_Dxy_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxy_a-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a+2*1*I_KINETIC_Fy2z_Dxy;
  abcd[30] = 4.0E0*I_KINETIC_K4x3z_Dxy_aa-2.0E0*2*I_KINETIC_H2x3z_Dxy_a-2.0E0*3*I_KINETIC_H2x3z_Dxy_a+2*1*I_KINETIC_F3z_Dxy;
  abcd[31] = 4.0E0*I_KINETIC_K3x4y_Dxy_aa-2.0E0*1*I_KINETIC_Hx4y_Dxy_a-2.0E0*2*I_KINETIC_Hx4y_Dxy_a;
  abcd[32] = 4.0E0*I_KINETIC_K3x3yz_Dxy_aa-2.0E0*1*I_KINETIC_Hx3yz_Dxy_a-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a;
  abcd[33] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Hx2y2z_Dxy_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxy_a;
  abcd[34] = 4.0E0*I_KINETIC_K3xy3z_Dxy_aa-2.0E0*1*I_KINETIC_Hxy3z_Dxy_a-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a;
  abcd[35] = 4.0E0*I_KINETIC_K3x4z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4z_Dxy_a-2.0E0*2*I_KINETIC_Hx4z_Dxy_a;
  abcd[36] = 4.0E0*I_KINETIC_K2x5y_Dxy_aa-2.0E0*1*I_KINETIC_H5y_Dxy_a;
  abcd[37] = 4.0E0*I_KINETIC_K2x4yz_Dxy_aa-2.0E0*1*I_KINETIC_H4yz_Dxy_a;
  abcd[38] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3y2z_Dxy_a;
  abcd[39] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*1*I_KINETIC_H2y3z_Dxy_a;
  abcd[40] = 4.0E0*I_KINETIC_K2xy4z_Dxy_aa-2.0E0*1*I_KINETIC_Hy4z_Dxy_a;
  abcd[41] = 4.0E0*I_KINETIC_K2x5z_Dxy_aa-2.0E0*1*I_KINETIC_H5z_Dxy_a;
  abcd[42] = 4.0E0*I_KINETIC_K7x_Dxz_aa-2.0E0*5*I_KINETIC_H5x_Dxz_a-2.0E0*6*I_KINETIC_H5x_Dxz_a+5*4*I_KINETIC_F3x_Dxz;
  abcd[43] = 4.0E0*I_KINETIC_K6xy_Dxz_aa-2.0E0*4*I_KINETIC_H4xy_Dxz_a-2.0E0*5*I_KINETIC_H4xy_Dxz_a+4*3*I_KINETIC_F2xy_Dxz;
  abcd[44] = 4.0E0*I_KINETIC_K6xz_Dxz_aa-2.0E0*4*I_KINETIC_H4xz_Dxz_a-2.0E0*5*I_KINETIC_H4xz_Dxz_a+4*3*I_KINETIC_F2xz_Dxz;
  abcd[45] = 4.0E0*I_KINETIC_K5x2y_Dxz_aa-2.0E0*3*I_KINETIC_H3x2y_Dxz_a-2.0E0*4*I_KINETIC_H3x2y_Dxz_a+3*2*I_KINETIC_Fx2y_Dxz;
  abcd[46] = 4.0E0*I_KINETIC_K5xyz_Dxz_aa-2.0E0*3*I_KINETIC_H3xyz_Dxz_a-2.0E0*4*I_KINETIC_H3xyz_Dxz_a+3*2*I_KINETIC_Fxyz_Dxz;
  abcd[47] = 4.0E0*I_KINETIC_K5x2z_Dxz_aa-2.0E0*3*I_KINETIC_H3x2z_Dxz_a-2.0E0*4*I_KINETIC_H3x2z_Dxz_a+3*2*I_KINETIC_Fx2z_Dxz;
  abcd[48] = 4.0E0*I_KINETIC_K4x3y_Dxz_aa-2.0E0*2*I_KINETIC_H2x3y_Dxz_a-2.0E0*3*I_KINETIC_H2x3y_Dxz_a+2*1*I_KINETIC_F3y_Dxz;
  abcd[49] = 4.0E0*I_KINETIC_K4x2yz_Dxz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxz_a-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a+2*1*I_KINETIC_F2yz_Dxz;
  abcd[50] = 4.0E0*I_KINETIC_K4xy2z_Dxz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxz_a-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a+2*1*I_KINETIC_Fy2z_Dxz;
  abcd[51] = 4.0E0*I_KINETIC_K4x3z_Dxz_aa-2.0E0*2*I_KINETIC_H2x3z_Dxz_a-2.0E0*3*I_KINETIC_H2x3z_Dxz_a+2*1*I_KINETIC_F3z_Dxz;
  abcd[52] = 4.0E0*I_KINETIC_K3x4y_Dxz_aa-2.0E0*1*I_KINETIC_Hx4y_Dxz_a-2.0E0*2*I_KINETIC_Hx4y_Dxz_a;
  abcd[53] = 4.0E0*I_KINETIC_K3x3yz_Dxz_aa-2.0E0*1*I_KINETIC_Hx3yz_Dxz_a-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a;
  abcd[54] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Hx2y2z_Dxz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxz_a;
  abcd[55] = 4.0E0*I_KINETIC_K3xy3z_Dxz_aa-2.0E0*1*I_KINETIC_Hxy3z_Dxz_a-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a;
  abcd[56] = 4.0E0*I_KINETIC_K3x4z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4z_Dxz_a-2.0E0*2*I_KINETIC_Hx4z_Dxz_a;
  abcd[57] = 4.0E0*I_KINETIC_K2x5y_Dxz_aa-2.0E0*1*I_KINETIC_H5y_Dxz_a;
  abcd[58] = 4.0E0*I_KINETIC_K2x4yz_Dxz_aa-2.0E0*1*I_KINETIC_H4yz_Dxz_a;
  abcd[59] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3y2z_Dxz_a;
  abcd[60] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*1*I_KINETIC_H2y3z_Dxz_a;
  abcd[61] = 4.0E0*I_KINETIC_K2xy4z_Dxz_aa-2.0E0*1*I_KINETIC_Hy4z_Dxz_a;
  abcd[62] = 4.0E0*I_KINETIC_K2x5z_Dxz_aa-2.0E0*1*I_KINETIC_H5z_Dxz_a;
  abcd[63] = 4.0E0*I_KINETIC_K7x_D2y_aa-2.0E0*5*I_KINETIC_H5x_D2y_a-2.0E0*6*I_KINETIC_H5x_D2y_a+5*4*I_KINETIC_F3x_D2y;
  abcd[64] = 4.0E0*I_KINETIC_K6xy_D2y_aa-2.0E0*4*I_KINETIC_H4xy_D2y_a-2.0E0*5*I_KINETIC_H4xy_D2y_a+4*3*I_KINETIC_F2xy_D2y;
  abcd[65] = 4.0E0*I_KINETIC_K6xz_D2y_aa-2.0E0*4*I_KINETIC_H4xz_D2y_a-2.0E0*5*I_KINETIC_H4xz_D2y_a+4*3*I_KINETIC_F2xz_D2y;
  abcd[66] = 4.0E0*I_KINETIC_K5x2y_D2y_aa-2.0E0*3*I_KINETIC_H3x2y_D2y_a-2.0E0*4*I_KINETIC_H3x2y_D2y_a+3*2*I_KINETIC_Fx2y_D2y;
  abcd[67] = 4.0E0*I_KINETIC_K5xyz_D2y_aa-2.0E0*3*I_KINETIC_H3xyz_D2y_a-2.0E0*4*I_KINETIC_H3xyz_D2y_a+3*2*I_KINETIC_Fxyz_D2y;
  abcd[68] = 4.0E0*I_KINETIC_K5x2z_D2y_aa-2.0E0*3*I_KINETIC_H3x2z_D2y_a-2.0E0*4*I_KINETIC_H3x2z_D2y_a+3*2*I_KINETIC_Fx2z_D2y;
  abcd[69] = 4.0E0*I_KINETIC_K4x3y_D2y_aa-2.0E0*2*I_KINETIC_H2x3y_D2y_a-2.0E0*3*I_KINETIC_H2x3y_D2y_a+2*1*I_KINETIC_F3y_D2y;
  abcd[70] = 4.0E0*I_KINETIC_K4x2yz_D2y_aa-2.0E0*2*I_KINETIC_H2x2yz_D2y_a-2.0E0*3*I_KINETIC_H2x2yz_D2y_a+2*1*I_KINETIC_F2yz_D2y;
  abcd[71] = 4.0E0*I_KINETIC_K4xy2z_D2y_aa-2.0E0*2*I_KINETIC_H2xy2z_D2y_a-2.0E0*3*I_KINETIC_H2xy2z_D2y_a+2*1*I_KINETIC_Fy2z_D2y;
  abcd[72] = 4.0E0*I_KINETIC_K4x3z_D2y_aa-2.0E0*2*I_KINETIC_H2x3z_D2y_a-2.0E0*3*I_KINETIC_H2x3z_D2y_a+2*1*I_KINETIC_F3z_D2y;
  abcd[73] = 4.0E0*I_KINETIC_K3x4y_D2y_aa-2.0E0*1*I_KINETIC_Hx4y_D2y_a-2.0E0*2*I_KINETIC_Hx4y_D2y_a;
  abcd[74] = 4.0E0*I_KINETIC_K3x3yz_D2y_aa-2.0E0*1*I_KINETIC_Hx3yz_D2y_a-2.0E0*2*I_KINETIC_Hx3yz_D2y_a;
  abcd[75] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_Hx2y2z_D2y_a-2.0E0*2*I_KINETIC_Hx2y2z_D2y_a;
  abcd[76] = 4.0E0*I_KINETIC_K3xy3z_D2y_aa-2.0E0*1*I_KINETIC_Hxy3z_D2y_a-2.0E0*2*I_KINETIC_Hxy3z_D2y_a;
  abcd[77] = 4.0E0*I_KINETIC_K3x4z_D2y_aa-2.0E0*1*I_KINETIC_Hx4z_D2y_a-2.0E0*2*I_KINETIC_Hx4z_D2y_a;
  abcd[78] = 4.0E0*I_KINETIC_K2x5y_D2y_aa-2.0E0*1*I_KINETIC_H5y_D2y_a;
  abcd[79] = 4.0E0*I_KINETIC_K2x4yz_D2y_aa-2.0E0*1*I_KINETIC_H4yz_D2y_a;
  abcd[80] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*1*I_KINETIC_H3y2z_D2y_a;
  abcd[81] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*1*I_KINETIC_H2y3z_D2y_a;
  abcd[82] = 4.0E0*I_KINETIC_K2xy4z_D2y_aa-2.0E0*1*I_KINETIC_Hy4z_D2y_a;
  abcd[83] = 4.0E0*I_KINETIC_K2x5z_D2y_aa-2.0E0*1*I_KINETIC_H5z_D2y_a;
  abcd[84] = 4.0E0*I_KINETIC_K7x_Dyz_aa-2.0E0*5*I_KINETIC_H5x_Dyz_a-2.0E0*6*I_KINETIC_H5x_Dyz_a+5*4*I_KINETIC_F3x_Dyz;
  abcd[85] = 4.0E0*I_KINETIC_K6xy_Dyz_aa-2.0E0*4*I_KINETIC_H4xy_Dyz_a-2.0E0*5*I_KINETIC_H4xy_Dyz_a+4*3*I_KINETIC_F2xy_Dyz;
  abcd[86] = 4.0E0*I_KINETIC_K6xz_Dyz_aa-2.0E0*4*I_KINETIC_H4xz_Dyz_a-2.0E0*5*I_KINETIC_H4xz_Dyz_a+4*3*I_KINETIC_F2xz_Dyz;
  abcd[87] = 4.0E0*I_KINETIC_K5x2y_Dyz_aa-2.0E0*3*I_KINETIC_H3x2y_Dyz_a-2.0E0*4*I_KINETIC_H3x2y_Dyz_a+3*2*I_KINETIC_Fx2y_Dyz;
  abcd[88] = 4.0E0*I_KINETIC_K5xyz_Dyz_aa-2.0E0*3*I_KINETIC_H3xyz_Dyz_a-2.0E0*4*I_KINETIC_H3xyz_Dyz_a+3*2*I_KINETIC_Fxyz_Dyz;
  abcd[89] = 4.0E0*I_KINETIC_K5x2z_Dyz_aa-2.0E0*3*I_KINETIC_H3x2z_Dyz_a-2.0E0*4*I_KINETIC_H3x2z_Dyz_a+3*2*I_KINETIC_Fx2z_Dyz;
  abcd[90] = 4.0E0*I_KINETIC_K4x3y_Dyz_aa-2.0E0*2*I_KINETIC_H2x3y_Dyz_a-2.0E0*3*I_KINETIC_H2x3y_Dyz_a+2*1*I_KINETIC_F3y_Dyz;
  abcd[91] = 4.0E0*I_KINETIC_K4x2yz_Dyz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dyz_a-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a+2*1*I_KINETIC_F2yz_Dyz;
  abcd[92] = 4.0E0*I_KINETIC_K4xy2z_Dyz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dyz_a-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a+2*1*I_KINETIC_Fy2z_Dyz;
  abcd[93] = 4.0E0*I_KINETIC_K4x3z_Dyz_aa-2.0E0*2*I_KINETIC_H2x3z_Dyz_a-2.0E0*3*I_KINETIC_H2x3z_Dyz_a+2*1*I_KINETIC_F3z_Dyz;
  abcd[94] = 4.0E0*I_KINETIC_K3x4y_Dyz_aa-2.0E0*1*I_KINETIC_Hx4y_Dyz_a-2.0E0*2*I_KINETIC_Hx4y_Dyz_a;
  abcd[95] = 4.0E0*I_KINETIC_K3x3yz_Dyz_aa-2.0E0*1*I_KINETIC_Hx3yz_Dyz_a-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a;
  abcd[96] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Hx2y2z_Dyz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dyz_a;
  abcd[97] = 4.0E0*I_KINETIC_K3xy3z_Dyz_aa-2.0E0*1*I_KINETIC_Hxy3z_Dyz_a-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a;
  abcd[98] = 4.0E0*I_KINETIC_K3x4z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4z_Dyz_a-2.0E0*2*I_KINETIC_Hx4z_Dyz_a;
  abcd[99] = 4.0E0*I_KINETIC_K2x5y_Dyz_aa-2.0E0*1*I_KINETIC_H5y_Dyz_a;
  abcd[100] = 4.0E0*I_KINETIC_K2x4yz_Dyz_aa-2.0E0*1*I_KINETIC_H4yz_Dyz_a;
  abcd[101] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3y2z_Dyz_a;
  abcd[102] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*1*I_KINETIC_H2y3z_Dyz_a;
  abcd[103] = 4.0E0*I_KINETIC_K2xy4z_Dyz_aa-2.0E0*1*I_KINETIC_Hy4z_Dyz_a;
  abcd[104] = 4.0E0*I_KINETIC_K2x5z_Dyz_aa-2.0E0*1*I_KINETIC_H5z_Dyz_a;
  abcd[105] = 4.0E0*I_KINETIC_K7x_D2z_aa-2.0E0*5*I_KINETIC_H5x_D2z_a-2.0E0*6*I_KINETIC_H5x_D2z_a+5*4*I_KINETIC_F3x_D2z;
  abcd[106] = 4.0E0*I_KINETIC_K6xy_D2z_aa-2.0E0*4*I_KINETIC_H4xy_D2z_a-2.0E0*5*I_KINETIC_H4xy_D2z_a+4*3*I_KINETIC_F2xy_D2z;
  abcd[107] = 4.0E0*I_KINETIC_K6xz_D2z_aa-2.0E0*4*I_KINETIC_H4xz_D2z_a-2.0E0*5*I_KINETIC_H4xz_D2z_a+4*3*I_KINETIC_F2xz_D2z;
  abcd[108] = 4.0E0*I_KINETIC_K5x2y_D2z_aa-2.0E0*3*I_KINETIC_H3x2y_D2z_a-2.0E0*4*I_KINETIC_H3x2y_D2z_a+3*2*I_KINETIC_Fx2y_D2z;
  abcd[109] = 4.0E0*I_KINETIC_K5xyz_D2z_aa-2.0E0*3*I_KINETIC_H3xyz_D2z_a-2.0E0*4*I_KINETIC_H3xyz_D2z_a+3*2*I_KINETIC_Fxyz_D2z;
  abcd[110] = 4.0E0*I_KINETIC_K5x2z_D2z_aa-2.0E0*3*I_KINETIC_H3x2z_D2z_a-2.0E0*4*I_KINETIC_H3x2z_D2z_a+3*2*I_KINETIC_Fx2z_D2z;
  abcd[111] = 4.0E0*I_KINETIC_K4x3y_D2z_aa-2.0E0*2*I_KINETIC_H2x3y_D2z_a-2.0E0*3*I_KINETIC_H2x3y_D2z_a+2*1*I_KINETIC_F3y_D2z;
  abcd[112] = 4.0E0*I_KINETIC_K4x2yz_D2z_aa-2.0E0*2*I_KINETIC_H2x2yz_D2z_a-2.0E0*3*I_KINETIC_H2x2yz_D2z_a+2*1*I_KINETIC_F2yz_D2z;
  abcd[113] = 4.0E0*I_KINETIC_K4xy2z_D2z_aa-2.0E0*2*I_KINETIC_H2xy2z_D2z_a-2.0E0*3*I_KINETIC_H2xy2z_D2z_a+2*1*I_KINETIC_Fy2z_D2z;
  abcd[114] = 4.0E0*I_KINETIC_K4x3z_D2z_aa-2.0E0*2*I_KINETIC_H2x3z_D2z_a-2.0E0*3*I_KINETIC_H2x3z_D2z_a+2*1*I_KINETIC_F3z_D2z;
  abcd[115] = 4.0E0*I_KINETIC_K3x4y_D2z_aa-2.0E0*1*I_KINETIC_Hx4y_D2z_a-2.0E0*2*I_KINETIC_Hx4y_D2z_a;
  abcd[116] = 4.0E0*I_KINETIC_K3x3yz_D2z_aa-2.0E0*1*I_KINETIC_Hx3yz_D2z_a-2.0E0*2*I_KINETIC_Hx3yz_D2z_a;
  abcd[117] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_Hx2y2z_D2z_a-2.0E0*2*I_KINETIC_Hx2y2z_D2z_a;
  abcd[118] = 4.0E0*I_KINETIC_K3xy3z_D2z_aa-2.0E0*1*I_KINETIC_Hxy3z_D2z_a-2.0E0*2*I_KINETIC_Hxy3z_D2z_a;
  abcd[119] = 4.0E0*I_KINETIC_K3x4z_D2z_aa-2.0E0*1*I_KINETIC_Hx4z_D2z_a-2.0E0*2*I_KINETIC_Hx4z_D2z_a;
  abcd[120] = 4.0E0*I_KINETIC_K2x5y_D2z_aa-2.0E0*1*I_KINETIC_H5y_D2z_a;
  abcd[121] = 4.0E0*I_KINETIC_K2x4yz_D2z_aa-2.0E0*1*I_KINETIC_H4yz_D2z_a;
  abcd[122] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*1*I_KINETIC_H3y2z_D2z_a;
  abcd[123] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*1*I_KINETIC_H2y3z_D2z_a;
  abcd[124] = 4.0E0*I_KINETIC_K2xy4z_D2z_aa-2.0E0*1*I_KINETIC_Hy4z_D2z_a;
  abcd[125] = 4.0E0*I_KINETIC_K2x5z_D2z_aa-2.0E0*1*I_KINETIC_H5z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[126] = 4.0E0*I_KINETIC_K6xy_D2x_aa-2.0E0*5*I_KINETIC_H4xy_D2x_a;
  abcd[127] = 4.0E0*I_KINETIC_K5x2y_D2x_aa-2.0E0*1*I_KINETIC_H5x_D2x_a-2.0E0*4*I_KINETIC_H3x2y_D2x_a+4*1*I_KINETIC_F3x_D2x;
  abcd[128] = 4.0E0*I_KINETIC_K5xyz_D2x_aa-2.0E0*4*I_KINETIC_H3xyz_D2x_a;
  abcd[129] = 4.0E0*I_KINETIC_K4x3y_D2x_aa-2.0E0*2*I_KINETIC_H4xy_D2x_a-2.0E0*3*I_KINETIC_H2x3y_D2x_a+3*2*I_KINETIC_F2xy_D2x;
  abcd[130] = 4.0E0*I_KINETIC_K4x2yz_D2x_aa-2.0E0*1*I_KINETIC_H4xz_D2x_a-2.0E0*3*I_KINETIC_H2x2yz_D2x_a+3*1*I_KINETIC_F2xz_D2x;
  abcd[131] = 4.0E0*I_KINETIC_K4xy2z_D2x_aa-2.0E0*3*I_KINETIC_H2xy2z_D2x_a;
  abcd[132] = 4.0E0*I_KINETIC_K3x4y_D2x_aa-2.0E0*3*I_KINETIC_H3x2y_D2x_a-2.0E0*2*I_KINETIC_Hx4y_D2x_a+2*3*I_KINETIC_Fx2y_D2x;
  abcd[133] = 4.0E0*I_KINETIC_K3x3yz_D2x_aa-2.0E0*2*I_KINETIC_H3xyz_D2x_a-2.0E0*2*I_KINETIC_Hx3yz_D2x_a+2*2*I_KINETIC_Fxyz_D2x;
  abcd[134] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_H3x2z_D2x_a-2.0E0*2*I_KINETIC_Hx2y2z_D2x_a+2*1*I_KINETIC_Fx2z_D2x;
  abcd[135] = 4.0E0*I_KINETIC_K3xy3z_D2x_aa-2.0E0*2*I_KINETIC_Hxy3z_D2x_a;
  abcd[136] = 4.0E0*I_KINETIC_K2x5y_D2x_aa-2.0E0*4*I_KINETIC_H2x3y_D2x_a-2.0E0*1*I_KINETIC_H5y_D2x_a+4*I_KINETIC_F3y_D2x;
  abcd[137] = 4.0E0*I_KINETIC_K2x4yz_D2x_aa-2.0E0*3*I_KINETIC_H2x2yz_D2x_a-2.0E0*1*I_KINETIC_H4yz_D2x_a+3*I_KINETIC_F2yz_D2x;
  abcd[138] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*2*I_KINETIC_H2xy2z_D2x_a-2.0E0*1*I_KINETIC_H3y2z_D2x_a+2*I_KINETIC_Fy2z_D2x;
  abcd[139] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*1*I_KINETIC_H2x3z_D2x_a-2.0E0*1*I_KINETIC_H2y3z_D2x_a+1*I_KINETIC_F3z_D2x;
  abcd[140] = 4.0E0*I_KINETIC_K2xy4z_D2x_aa-2.0E0*1*I_KINETIC_Hy4z_D2x_a;
  abcd[141] = 4.0E0*I_KINETIC_Kx6y_D2x_aa-2.0E0*5*I_KINETIC_Hx4y_D2x_a;
  abcd[142] = 4.0E0*I_KINETIC_Kx5yz_D2x_aa-2.0E0*4*I_KINETIC_Hx3yz_D2x_a;
  abcd[143] = 4.0E0*I_KINETIC_Kx4y2z_D2x_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a;
  abcd[144] = 4.0E0*I_KINETIC_Kx3y3z_D2x_aa-2.0E0*2*I_KINETIC_Hxy3z_D2x_a;
  abcd[145] = 4.0E0*I_KINETIC_Kx2y4z_D2x_aa-2.0E0*1*I_KINETIC_Hx4z_D2x_a;
  abcd[146] = 4.0E0*I_KINETIC_Kxy5z_D2x_aa;
  abcd[147] = 4.0E0*I_KINETIC_K6xy_Dxy_aa-2.0E0*5*I_KINETIC_H4xy_Dxy_a;
  abcd[148] = 4.0E0*I_KINETIC_K5x2y_Dxy_aa-2.0E0*1*I_KINETIC_H5x_Dxy_a-2.0E0*4*I_KINETIC_H3x2y_Dxy_a+4*1*I_KINETIC_F3x_Dxy;
  abcd[149] = 4.0E0*I_KINETIC_K5xyz_Dxy_aa-2.0E0*4*I_KINETIC_H3xyz_Dxy_a;
  abcd[150] = 4.0E0*I_KINETIC_K4x3y_Dxy_aa-2.0E0*2*I_KINETIC_H4xy_Dxy_a-2.0E0*3*I_KINETIC_H2x3y_Dxy_a+3*2*I_KINETIC_F2xy_Dxy;
  abcd[151] = 4.0E0*I_KINETIC_K4x2yz_Dxy_aa-2.0E0*1*I_KINETIC_H4xz_Dxy_a-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a+3*1*I_KINETIC_F2xz_Dxy;
  abcd[152] = 4.0E0*I_KINETIC_K4xy2z_Dxy_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a;
  abcd[153] = 4.0E0*I_KINETIC_K3x4y_Dxy_aa-2.0E0*3*I_KINETIC_H3x2y_Dxy_a-2.0E0*2*I_KINETIC_Hx4y_Dxy_a+2*3*I_KINETIC_Fx2y_Dxy;
  abcd[154] = 4.0E0*I_KINETIC_K3x3yz_Dxy_aa-2.0E0*2*I_KINETIC_H3xyz_Dxy_a-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a+2*2*I_KINETIC_Fxyz_Dxy;
  abcd[155] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3x2z_Dxy_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxy_a+2*1*I_KINETIC_Fx2z_Dxy;
  abcd[156] = 4.0E0*I_KINETIC_K3xy3z_Dxy_aa-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a;
  abcd[157] = 4.0E0*I_KINETIC_K2x5y_Dxy_aa-2.0E0*4*I_KINETIC_H2x3y_Dxy_a-2.0E0*1*I_KINETIC_H5y_Dxy_a+4*I_KINETIC_F3y_Dxy;
  abcd[158] = 4.0E0*I_KINETIC_K2x4yz_Dxy_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a-2.0E0*1*I_KINETIC_H4yz_Dxy_a+3*I_KINETIC_F2yz_Dxy;
  abcd[159] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxy_a-2.0E0*1*I_KINETIC_H3y2z_Dxy_a+2*I_KINETIC_Fy2z_Dxy;
  abcd[160] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*1*I_KINETIC_H2x3z_Dxy_a-2.0E0*1*I_KINETIC_H2y3z_Dxy_a+1*I_KINETIC_F3z_Dxy;
  abcd[161] = 4.0E0*I_KINETIC_K2xy4z_Dxy_aa-2.0E0*1*I_KINETIC_Hy4z_Dxy_a;
  abcd[162] = 4.0E0*I_KINETIC_Kx6y_Dxy_aa-2.0E0*5*I_KINETIC_Hx4y_Dxy_a;
  abcd[163] = 4.0E0*I_KINETIC_Kx5yz_Dxy_aa-2.0E0*4*I_KINETIC_Hx3yz_Dxy_a;
  abcd[164] = 4.0E0*I_KINETIC_Kx4y2z_Dxy_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a;
  abcd[165] = 4.0E0*I_KINETIC_Kx3y3z_Dxy_aa-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a;
  abcd[166] = 4.0E0*I_KINETIC_Kx2y4z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4z_Dxy_a;
  abcd[167] = 4.0E0*I_KINETIC_Kxy5z_Dxy_aa;
  abcd[168] = 4.0E0*I_KINETIC_K6xy_Dxz_aa-2.0E0*5*I_KINETIC_H4xy_Dxz_a;
  abcd[169] = 4.0E0*I_KINETIC_K5x2y_Dxz_aa-2.0E0*1*I_KINETIC_H5x_Dxz_a-2.0E0*4*I_KINETIC_H3x2y_Dxz_a+4*1*I_KINETIC_F3x_Dxz;
  abcd[170] = 4.0E0*I_KINETIC_K5xyz_Dxz_aa-2.0E0*4*I_KINETIC_H3xyz_Dxz_a;
  abcd[171] = 4.0E0*I_KINETIC_K4x3y_Dxz_aa-2.0E0*2*I_KINETIC_H4xy_Dxz_a-2.0E0*3*I_KINETIC_H2x3y_Dxz_a+3*2*I_KINETIC_F2xy_Dxz;
  abcd[172] = 4.0E0*I_KINETIC_K4x2yz_Dxz_aa-2.0E0*1*I_KINETIC_H4xz_Dxz_a-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a+3*1*I_KINETIC_F2xz_Dxz;
  abcd[173] = 4.0E0*I_KINETIC_K4xy2z_Dxz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a;
  abcd[174] = 4.0E0*I_KINETIC_K3x4y_Dxz_aa-2.0E0*3*I_KINETIC_H3x2y_Dxz_a-2.0E0*2*I_KINETIC_Hx4y_Dxz_a+2*3*I_KINETIC_Fx2y_Dxz;
  abcd[175] = 4.0E0*I_KINETIC_K3x3yz_Dxz_aa-2.0E0*2*I_KINETIC_H3xyz_Dxz_a-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a+2*2*I_KINETIC_Fxyz_Dxz;
  abcd[176] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3x2z_Dxz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxz_a+2*1*I_KINETIC_Fx2z_Dxz;
  abcd[177] = 4.0E0*I_KINETIC_K3xy3z_Dxz_aa-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a;
  abcd[178] = 4.0E0*I_KINETIC_K2x5y_Dxz_aa-2.0E0*4*I_KINETIC_H2x3y_Dxz_a-2.0E0*1*I_KINETIC_H5y_Dxz_a+4*I_KINETIC_F3y_Dxz;
  abcd[179] = 4.0E0*I_KINETIC_K2x4yz_Dxz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a-2.0E0*1*I_KINETIC_H4yz_Dxz_a+3*I_KINETIC_F2yz_Dxz;
  abcd[180] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxz_a-2.0E0*1*I_KINETIC_H3y2z_Dxz_a+2*I_KINETIC_Fy2z_Dxz;
  abcd[181] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*1*I_KINETIC_H2x3z_Dxz_a-2.0E0*1*I_KINETIC_H2y3z_Dxz_a+1*I_KINETIC_F3z_Dxz;
  abcd[182] = 4.0E0*I_KINETIC_K2xy4z_Dxz_aa-2.0E0*1*I_KINETIC_Hy4z_Dxz_a;
  abcd[183] = 4.0E0*I_KINETIC_Kx6y_Dxz_aa-2.0E0*5*I_KINETIC_Hx4y_Dxz_a;
  abcd[184] = 4.0E0*I_KINETIC_Kx5yz_Dxz_aa-2.0E0*4*I_KINETIC_Hx3yz_Dxz_a;
  abcd[185] = 4.0E0*I_KINETIC_Kx4y2z_Dxz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a;
  abcd[186] = 4.0E0*I_KINETIC_Kx3y3z_Dxz_aa-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a;
  abcd[187] = 4.0E0*I_KINETIC_Kx2y4z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4z_Dxz_a;
  abcd[188] = 4.0E0*I_KINETIC_Kxy5z_Dxz_aa;
  abcd[189] = 4.0E0*I_KINETIC_K6xy_D2y_aa-2.0E0*5*I_KINETIC_H4xy_D2y_a;
  abcd[190] = 4.0E0*I_KINETIC_K5x2y_D2y_aa-2.0E0*1*I_KINETIC_H5x_D2y_a-2.0E0*4*I_KINETIC_H3x2y_D2y_a+4*1*I_KINETIC_F3x_D2y;
  abcd[191] = 4.0E0*I_KINETIC_K5xyz_D2y_aa-2.0E0*4*I_KINETIC_H3xyz_D2y_a;
  abcd[192] = 4.0E0*I_KINETIC_K4x3y_D2y_aa-2.0E0*2*I_KINETIC_H4xy_D2y_a-2.0E0*3*I_KINETIC_H2x3y_D2y_a+3*2*I_KINETIC_F2xy_D2y;
  abcd[193] = 4.0E0*I_KINETIC_K4x2yz_D2y_aa-2.0E0*1*I_KINETIC_H4xz_D2y_a-2.0E0*3*I_KINETIC_H2x2yz_D2y_a+3*1*I_KINETIC_F2xz_D2y;
  abcd[194] = 4.0E0*I_KINETIC_K4xy2z_D2y_aa-2.0E0*3*I_KINETIC_H2xy2z_D2y_a;
  abcd[195] = 4.0E0*I_KINETIC_K3x4y_D2y_aa-2.0E0*3*I_KINETIC_H3x2y_D2y_a-2.0E0*2*I_KINETIC_Hx4y_D2y_a+2*3*I_KINETIC_Fx2y_D2y;
  abcd[196] = 4.0E0*I_KINETIC_K3x3yz_D2y_aa-2.0E0*2*I_KINETIC_H3xyz_D2y_a-2.0E0*2*I_KINETIC_Hx3yz_D2y_a+2*2*I_KINETIC_Fxyz_D2y;
  abcd[197] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_H3x2z_D2y_a-2.0E0*2*I_KINETIC_Hx2y2z_D2y_a+2*1*I_KINETIC_Fx2z_D2y;
  abcd[198] = 4.0E0*I_KINETIC_K3xy3z_D2y_aa-2.0E0*2*I_KINETIC_Hxy3z_D2y_a;
  abcd[199] = 4.0E0*I_KINETIC_K2x5y_D2y_aa-2.0E0*4*I_KINETIC_H2x3y_D2y_a-2.0E0*1*I_KINETIC_H5y_D2y_a+4*I_KINETIC_F3y_D2y;
  abcd[200] = 4.0E0*I_KINETIC_K2x4yz_D2y_aa-2.0E0*3*I_KINETIC_H2x2yz_D2y_a-2.0E0*1*I_KINETIC_H4yz_D2y_a+3*I_KINETIC_F2yz_D2y;
  abcd[201] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*2*I_KINETIC_H2xy2z_D2y_a-2.0E0*1*I_KINETIC_H3y2z_D2y_a+2*I_KINETIC_Fy2z_D2y;
  abcd[202] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*1*I_KINETIC_H2x3z_D2y_a-2.0E0*1*I_KINETIC_H2y3z_D2y_a+1*I_KINETIC_F3z_D2y;
  abcd[203] = 4.0E0*I_KINETIC_K2xy4z_D2y_aa-2.0E0*1*I_KINETIC_Hy4z_D2y_a;
  abcd[204] = 4.0E0*I_KINETIC_Kx6y_D2y_aa-2.0E0*5*I_KINETIC_Hx4y_D2y_a;
  abcd[205] = 4.0E0*I_KINETIC_Kx5yz_D2y_aa-2.0E0*4*I_KINETIC_Hx3yz_D2y_a;
  abcd[206] = 4.0E0*I_KINETIC_Kx4y2z_D2y_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a;
  abcd[207] = 4.0E0*I_KINETIC_Kx3y3z_D2y_aa-2.0E0*2*I_KINETIC_Hxy3z_D2y_a;
  abcd[208] = 4.0E0*I_KINETIC_Kx2y4z_D2y_aa-2.0E0*1*I_KINETIC_Hx4z_D2y_a;
  abcd[209] = 4.0E0*I_KINETIC_Kxy5z_D2y_aa;
  abcd[210] = 4.0E0*I_KINETIC_K6xy_Dyz_aa-2.0E0*5*I_KINETIC_H4xy_Dyz_a;
  abcd[211] = 4.0E0*I_KINETIC_K5x2y_Dyz_aa-2.0E0*1*I_KINETIC_H5x_Dyz_a-2.0E0*4*I_KINETIC_H3x2y_Dyz_a+4*1*I_KINETIC_F3x_Dyz;
  abcd[212] = 4.0E0*I_KINETIC_K5xyz_Dyz_aa-2.0E0*4*I_KINETIC_H3xyz_Dyz_a;
  abcd[213] = 4.0E0*I_KINETIC_K4x3y_Dyz_aa-2.0E0*2*I_KINETIC_H4xy_Dyz_a-2.0E0*3*I_KINETIC_H2x3y_Dyz_a+3*2*I_KINETIC_F2xy_Dyz;
  abcd[214] = 4.0E0*I_KINETIC_K4x2yz_Dyz_aa-2.0E0*1*I_KINETIC_H4xz_Dyz_a-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a+3*1*I_KINETIC_F2xz_Dyz;
  abcd[215] = 4.0E0*I_KINETIC_K4xy2z_Dyz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a;
  abcd[216] = 4.0E0*I_KINETIC_K3x4y_Dyz_aa-2.0E0*3*I_KINETIC_H3x2y_Dyz_a-2.0E0*2*I_KINETIC_Hx4y_Dyz_a+2*3*I_KINETIC_Fx2y_Dyz;
  abcd[217] = 4.0E0*I_KINETIC_K3x3yz_Dyz_aa-2.0E0*2*I_KINETIC_H3xyz_Dyz_a-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a+2*2*I_KINETIC_Fxyz_Dyz;
  abcd[218] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3x2z_Dyz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dyz_a+2*1*I_KINETIC_Fx2z_Dyz;
  abcd[219] = 4.0E0*I_KINETIC_K3xy3z_Dyz_aa-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a;
  abcd[220] = 4.0E0*I_KINETIC_K2x5y_Dyz_aa-2.0E0*4*I_KINETIC_H2x3y_Dyz_a-2.0E0*1*I_KINETIC_H5y_Dyz_a+4*I_KINETIC_F3y_Dyz;
  abcd[221] = 4.0E0*I_KINETIC_K2x4yz_Dyz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a-2.0E0*1*I_KINETIC_H4yz_Dyz_a+3*I_KINETIC_F2yz_Dyz;
  abcd[222] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dyz_a-2.0E0*1*I_KINETIC_H3y2z_Dyz_a+2*I_KINETIC_Fy2z_Dyz;
  abcd[223] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*1*I_KINETIC_H2x3z_Dyz_a-2.0E0*1*I_KINETIC_H2y3z_Dyz_a+1*I_KINETIC_F3z_Dyz;
  abcd[224] = 4.0E0*I_KINETIC_K2xy4z_Dyz_aa-2.0E0*1*I_KINETIC_Hy4z_Dyz_a;
  abcd[225] = 4.0E0*I_KINETIC_Kx6y_Dyz_aa-2.0E0*5*I_KINETIC_Hx4y_Dyz_a;
  abcd[226] = 4.0E0*I_KINETIC_Kx5yz_Dyz_aa-2.0E0*4*I_KINETIC_Hx3yz_Dyz_a;
  abcd[227] = 4.0E0*I_KINETIC_Kx4y2z_Dyz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a;
  abcd[228] = 4.0E0*I_KINETIC_Kx3y3z_Dyz_aa-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a;
  abcd[229] = 4.0E0*I_KINETIC_Kx2y4z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4z_Dyz_a;
  abcd[230] = 4.0E0*I_KINETIC_Kxy5z_Dyz_aa;
  abcd[231] = 4.0E0*I_KINETIC_K6xy_D2z_aa-2.0E0*5*I_KINETIC_H4xy_D2z_a;
  abcd[232] = 4.0E0*I_KINETIC_K5x2y_D2z_aa-2.0E0*1*I_KINETIC_H5x_D2z_a-2.0E0*4*I_KINETIC_H3x2y_D2z_a+4*1*I_KINETIC_F3x_D2z;
  abcd[233] = 4.0E0*I_KINETIC_K5xyz_D2z_aa-2.0E0*4*I_KINETIC_H3xyz_D2z_a;
  abcd[234] = 4.0E0*I_KINETIC_K4x3y_D2z_aa-2.0E0*2*I_KINETIC_H4xy_D2z_a-2.0E0*3*I_KINETIC_H2x3y_D2z_a+3*2*I_KINETIC_F2xy_D2z;
  abcd[235] = 4.0E0*I_KINETIC_K4x2yz_D2z_aa-2.0E0*1*I_KINETIC_H4xz_D2z_a-2.0E0*3*I_KINETIC_H2x2yz_D2z_a+3*1*I_KINETIC_F2xz_D2z;
  abcd[236] = 4.0E0*I_KINETIC_K4xy2z_D2z_aa-2.0E0*3*I_KINETIC_H2xy2z_D2z_a;
  abcd[237] = 4.0E0*I_KINETIC_K3x4y_D2z_aa-2.0E0*3*I_KINETIC_H3x2y_D2z_a-2.0E0*2*I_KINETIC_Hx4y_D2z_a+2*3*I_KINETIC_Fx2y_D2z;
  abcd[238] = 4.0E0*I_KINETIC_K3x3yz_D2z_aa-2.0E0*2*I_KINETIC_H3xyz_D2z_a-2.0E0*2*I_KINETIC_Hx3yz_D2z_a+2*2*I_KINETIC_Fxyz_D2z;
  abcd[239] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_H3x2z_D2z_a-2.0E0*2*I_KINETIC_Hx2y2z_D2z_a+2*1*I_KINETIC_Fx2z_D2z;
  abcd[240] = 4.0E0*I_KINETIC_K3xy3z_D2z_aa-2.0E0*2*I_KINETIC_Hxy3z_D2z_a;
  abcd[241] = 4.0E0*I_KINETIC_K2x5y_D2z_aa-2.0E0*4*I_KINETIC_H2x3y_D2z_a-2.0E0*1*I_KINETIC_H5y_D2z_a+4*I_KINETIC_F3y_D2z;
  abcd[242] = 4.0E0*I_KINETIC_K2x4yz_D2z_aa-2.0E0*3*I_KINETIC_H2x2yz_D2z_a-2.0E0*1*I_KINETIC_H4yz_D2z_a+3*I_KINETIC_F2yz_D2z;
  abcd[243] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*2*I_KINETIC_H2xy2z_D2z_a-2.0E0*1*I_KINETIC_H3y2z_D2z_a+2*I_KINETIC_Fy2z_D2z;
  abcd[244] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*1*I_KINETIC_H2x3z_D2z_a-2.0E0*1*I_KINETIC_H2y3z_D2z_a+1*I_KINETIC_F3z_D2z;
  abcd[245] = 4.0E0*I_KINETIC_K2xy4z_D2z_aa-2.0E0*1*I_KINETIC_Hy4z_D2z_a;
  abcd[246] = 4.0E0*I_KINETIC_Kx6y_D2z_aa-2.0E0*5*I_KINETIC_Hx4y_D2z_a;
  abcd[247] = 4.0E0*I_KINETIC_Kx5yz_D2z_aa-2.0E0*4*I_KINETIC_Hx3yz_D2z_a;
  abcd[248] = 4.0E0*I_KINETIC_Kx4y2z_D2z_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a;
  abcd[249] = 4.0E0*I_KINETIC_Kx3y3z_D2z_aa-2.0E0*2*I_KINETIC_Hxy3z_D2z_a;
  abcd[250] = 4.0E0*I_KINETIC_Kx2y4z_D2z_aa-2.0E0*1*I_KINETIC_Hx4z_D2z_a;
  abcd[251] = 4.0E0*I_KINETIC_Kxy5z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[252] = 4.0E0*I_KINETIC_K6xz_D2x_aa-2.0E0*5*I_KINETIC_H4xz_D2x_a;
  abcd[253] = 4.0E0*I_KINETIC_K5xyz_D2x_aa-2.0E0*4*I_KINETIC_H3xyz_D2x_a;
  abcd[254] = 4.0E0*I_KINETIC_K5x2z_D2x_aa-2.0E0*1*I_KINETIC_H5x_D2x_a-2.0E0*4*I_KINETIC_H3x2z_D2x_a+4*1*I_KINETIC_F3x_D2x;
  abcd[255] = 4.0E0*I_KINETIC_K4x2yz_D2x_aa-2.0E0*3*I_KINETIC_H2x2yz_D2x_a;
  abcd[256] = 4.0E0*I_KINETIC_K4xy2z_D2x_aa-2.0E0*1*I_KINETIC_H4xy_D2x_a-2.0E0*3*I_KINETIC_H2xy2z_D2x_a+3*1*I_KINETIC_F2xy_D2x;
  abcd[257] = 4.0E0*I_KINETIC_K4x3z_D2x_aa-2.0E0*2*I_KINETIC_H4xz_D2x_a-2.0E0*3*I_KINETIC_H2x3z_D2x_a+3*2*I_KINETIC_F2xz_D2x;
  abcd[258] = 4.0E0*I_KINETIC_K3x3yz_D2x_aa-2.0E0*2*I_KINETIC_Hx3yz_D2x_a;
  abcd[259] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_H3x2y_D2x_a-2.0E0*2*I_KINETIC_Hx2y2z_D2x_a+2*1*I_KINETIC_Fx2y_D2x;
  abcd[260] = 4.0E0*I_KINETIC_K3xy3z_D2x_aa-2.0E0*2*I_KINETIC_H3xyz_D2x_a-2.0E0*2*I_KINETIC_Hxy3z_D2x_a+2*2*I_KINETIC_Fxyz_D2x;
  abcd[261] = 4.0E0*I_KINETIC_K3x4z_D2x_aa-2.0E0*3*I_KINETIC_H3x2z_D2x_a-2.0E0*2*I_KINETIC_Hx4z_D2x_a+2*3*I_KINETIC_Fx2z_D2x;
  abcd[262] = 4.0E0*I_KINETIC_K2x4yz_D2x_aa-2.0E0*1*I_KINETIC_H4yz_D2x_a;
  abcd[263] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*1*I_KINETIC_H2x3y_D2x_a-2.0E0*1*I_KINETIC_H3y2z_D2x_a+1*I_KINETIC_F3y_D2x;
  abcd[264] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*2*I_KINETIC_H2x2yz_D2x_a-2.0E0*1*I_KINETIC_H2y3z_D2x_a+2*I_KINETIC_F2yz_D2x;
  abcd[265] = 4.0E0*I_KINETIC_K2xy4z_D2x_aa-2.0E0*3*I_KINETIC_H2xy2z_D2x_a-2.0E0*1*I_KINETIC_Hy4z_D2x_a+3*I_KINETIC_Fy2z_D2x;
  abcd[266] = 4.0E0*I_KINETIC_K2x5z_D2x_aa-2.0E0*4*I_KINETIC_H2x3z_D2x_a-2.0E0*1*I_KINETIC_H5z_D2x_a+4*I_KINETIC_F3z_D2x;
  abcd[267] = 4.0E0*I_KINETIC_Kx5yz_D2x_aa;
  abcd[268] = 4.0E0*I_KINETIC_Kx4y2z_D2x_aa-2.0E0*1*I_KINETIC_Hx4y_D2x_a;
  abcd[269] = 4.0E0*I_KINETIC_Kx3y3z_D2x_aa-2.0E0*2*I_KINETIC_Hx3yz_D2x_a;
  abcd[270] = 4.0E0*I_KINETIC_Kx2y4z_D2x_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a;
  abcd[271] = 4.0E0*I_KINETIC_Kxy5z_D2x_aa-2.0E0*4*I_KINETIC_Hxy3z_D2x_a;
  abcd[272] = 4.0E0*I_KINETIC_Kx6z_D2x_aa-2.0E0*5*I_KINETIC_Hx4z_D2x_a;
  abcd[273] = 4.0E0*I_KINETIC_K6xz_Dxy_aa-2.0E0*5*I_KINETIC_H4xz_Dxy_a;
  abcd[274] = 4.0E0*I_KINETIC_K5xyz_Dxy_aa-2.0E0*4*I_KINETIC_H3xyz_Dxy_a;
  abcd[275] = 4.0E0*I_KINETIC_K5x2z_Dxy_aa-2.0E0*1*I_KINETIC_H5x_Dxy_a-2.0E0*4*I_KINETIC_H3x2z_Dxy_a+4*1*I_KINETIC_F3x_Dxy;
  abcd[276] = 4.0E0*I_KINETIC_K4x2yz_Dxy_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a;
  abcd[277] = 4.0E0*I_KINETIC_K4xy2z_Dxy_aa-2.0E0*1*I_KINETIC_H4xy_Dxy_a-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a+3*1*I_KINETIC_F2xy_Dxy;
  abcd[278] = 4.0E0*I_KINETIC_K4x3z_Dxy_aa-2.0E0*2*I_KINETIC_H4xz_Dxy_a-2.0E0*3*I_KINETIC_H2x3z_Dxy_a+3*2*I_KINETIC_F2xz_Dxy;
  abcd[279] = 4.0E0*I_KINETIC_K3x3yz_Dxy_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a;
  abcd[280] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3x2y_Dxy_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxy_a+2*1*I_KINETIC_Fx2y_Dxy;
  abcd[281] = 4.0E0*I_KINETIC_K3xy3z_Dxy_aa-2.0E0*2*I_KINETIC_H3xyz_Dxy_a-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a+2*2*I_KINETIC_Fxyz_Dxy;
  abcd[282] = 4.0E0*I_KINETIC_K3x4z_Dxy_aa-2.0E0*3*I_KINETIC_H3x2z_Dxy_a-2.0E0*2*I_KINETIC_Hx4z_Dxy_a+2*3*I_KINETIC_Fx2z_Dxy;
  abcd[283] = 4.0E0*I_KINETIC_K2x4yz_Dxy_aa-2.0E0*1*I_KINETIC_H4yz_Dxy_a;
  abcd[284] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*1*I_KINETIC_H2x3y_Dxy_a-2.0E0*1*I_KINETIC_H3y2z_Dxy_a+1*I_KINETIC_F3y_Dxy;
  abcd[285] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxy_a-2.0E0*1*I_KINETIC_H2y3z_Dxy_a+2*I_KINETIC_F2yz_Dxy;
  abcd[286] = 4.0E0*I_KINETIC_K2xy4z_Dxy_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a-2.0E0*1*I_KINETIC_Hy4z_Dxy_a+3*I_KINETIC_Fy2z_Dxy;
  abcd[287] = 4.0E0*I_KINETIC_K2x5z_Dxy_aa-2.0E0*4*I_KINETIC_H2x3z_Dxy_a-2.0E0*1*I_KINETIC_H5z_Dxy_a+4*I_KINETIC_F3z_Dxy;
  abcd[288] = 4.0E0*I_KINETIC_Kx5yz_Dxy_aa;
  abcd[289] = 4.0E0*I_KINETIC_Kx4y2z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4y_Dxy_a;
  abcd[290] = 4.0E0*I_KINETIC_Kx3y3z_Dxy_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a;
  abcd[291] = 4.0E0*I_KINETIC_Kx2y4z_Dxy_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a;
  abcd[292] = 4.0E0*I_KINETIC_Kxy5z_Dxy_aa-2.0E0*4*I_KINETIC_Hxy3z_Dxy_a;
  abcd[293] = 4.0E0*I_KINETIC_Kx6z_Dxy_aa-2.0E0*5*I_KINETIC_Hx4z_Dxy_a;
  abcd[294] = 4.0E0*I_KINETIC_K6xz_Dxz_aa-2.0E0*5*I_KINETIC_H4xz_Dxz_a;
  abcd[295] = 4.0E0*I_KINETIC_K5xyz_Dxz_aa-2.0E0*4*I_KINETIC_H3xyz_Dxz_a;
  abcd[296] = 4.0E0*I_KINETIC_K5x2z_Dxz_aa-2.0E0*1*I_KINETIC_H5x_Dxz_a-2.0E0*4*I_KINETIC_H3x2z_Dxz_a+4*1*I_KINETIC_F3x_Dxz;
  abcd[297] = 4.0E0*I_KINETIC_K4x2yz_Dxz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a;
  abcd[298] = 4.0E0*I_KINETIC_K4xy2z_Dxz_aa-2.0E0*1*I_KINETIC_H4xy_Dxz_a-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a+3*1*I_KINETIC_F2xy_Dxz;
  abcd[299] = 4.0E0*I_KINETIC_K4x3z_Dxz_aa-2.0E0*2*I_KINETIC_H4xz_Dxz_a-2.0E0*3*I_KINETIC_H2x3z_Dxz_a+3*2*I_KINETIC_F2xz_Dxz;
  abcd[300] = 4.0E0*I_KINETIC_K3x3yz_Dxz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a;
  abcd[301] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3x2y_Dxz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dxz_a+2*1*I_KINETIC_Fx2y_Dxz;
  abcd[302] = 4.0E0*I_KINETIC_K3xy3z_Dxz_aa-2.0E0*2*I_KINETIC_H3xyz_Dxz_a-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a+2*2*I_KINETIC_Fxyz_Dxz;
  abcd[303] = 4.0E0*I_KINETIC_K3x4z_Dxz_aa-2.0E0*3*I_KINETIC_H3x2z_Dxz_a-2.0E0*2*I_KINETIC_Hx4z_Dxz_a+2*3*I_KINETIC_Fx2z_Dxz;
  abcd[304] = 4.0E0*I_KINETIC_K2x4yz_Dxz_aa-2.0E0*1*I_KINETIC_H4yz_Dxz_a;
  abcd[305] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*1*I_KINETIC_H2x3y_Dxz_a-2.0E0*1*I_KINETIC_H3y2z_Dxz_a+1*I_KINETIC_F3y_Dxz;
  abcd[306] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxz_a-2.0E0*1*I_KINETIC_H2y3z_Dxz_a+2*I_KINETIC_F2yz_Dxz;
  abcd[307] = 4.0E0*I_KINETIC_K2xy4z_Dxz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a-2.0E0*1*I_KINETIC_Hy4z_Dxz_a+3*I_KINETIC_Fy2z_Dxz;
  abcd[308] = 4.0E0*I_KINETIC_K2x5z_Dxz_aa-2.0E0*4*I_KINETIC_H2x3z_Dxz_a-2.0E0*1*I_KINETIC_H5z_Dxz_a+4*I_KINETIC_F3z_Dxz;
  abcd[309] = 4.0E0*I_KINETIC_Kx5yz_Dxz_aa;
  abcd[310] = 4.0E0*I_KINETIC_Kx4y2z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4y_Dxz_a;
  abcd[311] = 4.0E0*I_KINETIC_Kx3y3z_Dxz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a;
  abcd[312] = 4.0E0*I_KINETIC_Kx2y4z_Dxz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a;
  abcd[313] = 4.0E0*I_KINETIC_Kxy5z_Dxz_aa-2.0E0*4*I_KINETIC_Hxy3z_Dxz_a;
  abcd[314] = 4.0E0*I_KINETIC_Kx6z_Dxz_aa-2.0E0*5*I_KINETIC_Hx4z_Dxz_a;
  abcd[315] = 4.0E0*I_KINETIC_K6xz_D2y_aa-2.0E0*5*I_KINETIC_H4xz_D2y_a;
  abcd[316] = 4.0E0*I_KINETIC_K5xyz_D2y_aa-2.0E0*4*I_KINETIC_H3xyz_D2y_a;
  abcd[317] = 4.0E0*I_KINETIC_K5x2z_D2y_aa-2.0E0*1*I_KINETIC_H5x_D2y_a-2.0E0*4*I_KINETIC_H3x2z_D2y_a+4*1*I_KINETIC_F3x_D2y;
  abcd[318] = 4.0E0*I_KINETIC_K4x2yz_D2y_aa-2.0E0*3*I_KINETIC_H2x2yz_D2y_a;
  abcd[319] = 4.0E0*I_KINETIC_K4xy2z_D2y_aa-2.0E0*1*I_KINETIC_H4xy_D2y_a-2.0E0*3*I_KINETIC_H2xy2z_D2y_a+3*1*I_KINETIC_F2xy_D2y;
  abcd[320] = 4.0E0*I_KINETIC_K4x3z_D2y_aa-2.0E0*2*I_KINETIC_H4xz_D2y_a-2.0E0*3*I_KINETIC_H2x3z_D2y_a+3*2*I_KINETIC_F2xz_D2y;
  abcd[321] = 4.0E0*I_KINETIC_K3x3yz_D2y_aa-2.0E0*2*I_KINETIC_Hx3yz_D2y_a;
  abcd[322] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_H3x2y_D2y_a-2.0E0*2*I_KINETIC_Hx2y2z_D2y_a+2*1*I_KINETIC_Fx2y_D2y;
  abcd[323] = 4.0E0*I_KINETIC_K3xy3z_D2y_aa-2.0E0*2*I_KINETIC_H3xyz_D2y_a-2.0E0*2*I_KINETIC_Hxy3z_D2y_a+2*2*I_KINETIC_Fxyz_D2y;
  abcd[324] = 4.0E0*I_KINETIC_K3x4z_D2y_aa-2.0E0*3*I_KINETIC_H3x2z_D2y_a-2.0E0*2*I_KINETIC_Hx4z_D2y_a+2*3*I_KINETIC_Fx2z_D2y;
  abcd[325] = 4.0E0*I_KINETIC_K2x4yz_D2y_aa-2.0E0*1*I_KINETIC_H4yz_D2y_a;
  abcd[326] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*1*I_KINETIC_H2x3y_D2y_a-2.0E0*1*I_KINETIC_H3y2z_D2y_a+1*I_KINETIC_F3y_D2y;
  abcd[327] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*2*I_KINETIC_H2x2yz_D2y_a-2.0E0*1*I_KINETIC_H2y3z_D2y_a+2*I_KINETIC_F2yz_D2y;
  abcd[328] = 4.0E0*I_KINETIC_K2xy4z_D2y_aa-2.0E0*3*I_KINETIC_H2xy2z_D2y_a-2.0E0*1*I_KINETIC_Hy4z_D2y_a+3*I_KINETIC_Fy2z_D2y;
  abcd[329] = 4.0E0*I_KINETIC_K2x5z_D2y_aa-2.0E0*4*I_KINETIC_H2x3z_D2y_a-2.0E0*1*I_KINETIC_H5z_D2y_a+4*I_KINETIC_F3z_D2y;
  abcd[330] = 4.0E0*I_KINETIC_Kx5yz_D2y_aa;
  abcd[331] = 4.0E0*I_KINETIC_Kx4y2z_D2y_aa-2.0E0*1*I_KINETIC_Hx4y_D2y_a;
  abcd[332] = 4.0E0*I_KINETIC_Kx3y3z_D2y_aa-2.0E0*2*I_KINETIC_Hx3yz_D2y_a;
  abcd[333] = 4.0E0*I_KINETIC_Kx2y4z_D2y_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a;
  abcd[334] = 4.0E0*I_KINETIC_Kxy5z_D2y_aa-2.0E0*4*I_KINETIC_Hxy3z_D2y_a;
  abcd[335] = 4.0E0*I_KINETIC_Kx6z_D2y_aa-2.0E0*5*I_KINETIC_Hx4z_D2y_a;
  abcd[336] = 4.0E0*I_KINETIC_K6xz_Dyz_aa-2.0E0*5*I_KINETIC_H4xz_Dyz_a;
  abcd[337] = 4.0E0*I_KINETIC_K5xyz_Dyz_aa-2.0E0*4*I_KINETIC_H3xyz_Dyz_a;
  abcd[338] = 4.0E0*I_KINETIC_K5x2z_Dyz_aa-2.0E0*1*I_KINETIC_H5x_Dyz_a-2.0E0*4*I_KINETIC_H3x2z_Dyz_a+4*1*I_KINETIC_F3x_Dyz;
  abcd[339] = 4.0E0*I_KINETIC_K4x2yz_Dyz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a;
  abcd[340] = 4.0E0*I_KINETIC_K4xy2z_Dyz_aa-2.0E0*1*I_KINETIC_H4xy_Dyz_a-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a+3*1*I_KINETIC_F2xy_Dyz;
  abcd[341] = 4.0E0*I_KINETIC_K4x3z_Dyz_aa-2.0E0*2*I_KINETIC_H4xz_Dyz_a-2.0E0*3*I_KINETIC_H2x3z_Dyz_a+3*2*I_KINETIC_F2xz_Dyz;
  abcd[342] = 4.0E0*I_KINETIC_K3x3yz_Dyz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a;
  abcd[343] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3x2y_Dyz_a-2.0E0*2*I_KINETIC_Hx2y2z_Dyz_a+2*1*I_KINETIC_Fx2y_Dyz;
  abcd[344] = 4.0E0*I_KINETIC_K3xy3z_Dyz_aa-2.0E0*2*I_KINETIC_H3xyz_Dyz_a-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a+2*2*I_KINETIC_Fxyz_Dyz;
  abcd[345] = 4.0E0*I_KINETIC_K3x4z_Dyz_aa-2.0E0*3*I_KINETIC_H3x2z_Dyz_a-2.0E0*2*I_KINETIC_Hx4z_Dyz_a+2*3*I_KINETIC_Fx2z_Dyz;
  abcd[346] = 4.0E0*I_KINETIC_K2x4yz_Dyz_aa-2.0E0*1*I_KINETIC_H4yz_Dyz_a;
  abcd[347] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*1*I_KINETIC_H2x3y_Dyz_a-2.0E0*1*I_KINETIC_H3y2z_Dyz_a+1*I_KINETIC_F3y_Dyz;
  abcd[348] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dyz_a-2.0E0*1*I_KINETIC_H2y3z_Dyz_a+2*I_KINETIC_F2yz_Dyz;
  abcd[349] = 4.0E0*I_KINETIC_K2xy4z_Dyz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a-2.0E0*1*I_KINETIC_Hy4z_Dyz_a+3*I_KINETIC_Fy2z_Dyz;
  abcd[350] = 4.0E0*I_KINETIC_K2x5z_Dyz_aa-2.0E0*4*I_KINETIC_H2x3z_Dyz_a-2.0E0*1*I_KINETIC_H5z_Dyz_a+4*I_KINETIC_F3z_Dyz;
  abcd[351] = 4.0E0*I_KINETIC_Kx5yz_Dyz_aa;
  abcd[352] = 4.0E0*I_KINETIC_Kx4y2z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4y_Dyz_a;
  abcd[353] = 4.0E0*I_KINETIC_Kx3y3z_Dyz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a;
  abcd[354] = 4.0E0*I_KINETIC_Kx2y4z_Dyz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a;
  abcd[355] = 4.0E0*I_KINETIC_Kxy5z_Dyz_aa-2.0E0*4*I_KINETIC_Hxy3z_Dyz_a;
  abcd[356] = 4.0E0*I_KINETIC_Kx6z_Dyz_aa-2.0E0*5*I_KINETIC_Hx4z_Dyz_a;
  abcd[357] = 4.0E0*I_KINETIC_K6xz_D2z_aa-2.0E0*5*I_KINETIC_H4xz_D2z_a;
  abcd[358] = 4.0E0*I_KINETIC_K5xyz_D2z_aa-2.0E0*4*I_KINETIC_H3xyz_D2z_a;
  abcd[359] = 4.0E0*I_KINETIC_K5x2z_D2z_aa-2.0E0*1*I_KINETIC_H5x_D2z_a-2.0E0*4*I_KINETIC_H3x2z_D2z_a+4*1*I_KINETIC_F3x_D2z;
  abcd[360] = 4.0E0*I_KINETIC_K4x2yz_D2z_aa-2.0E0*3*I_KINETIC_H2x2yz_D2z_a;
  abcd[361] = 4.0E0*I_KINETIC_K4xy2z_D2z_aa-2.0E0*1*I_KINETIC_H4xy_D2z_a-2.0E0*3*I_KINETIC_H2xy2z_D2z_a+3*1*I_KINETIC_F2xy_D2z;
  abcd[362] = 4.0E0*I_KINETIC_K4x3z_D2z_aa-2.0E0*2*I_KINETIC_H4xz_D2z_a-2.0E0*3*I_KINETIC_H2x3z_D2z_a+3*2*I_KINETIC_F2xz_D2z;
  abcd[363] = 4.0E0*I_KINETIC_K3x3yz_D2z_aa-2.0E0*2*I_KINETIC_Hx3yz_D2z_a;
  abcd[364] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_H3x2y_D2z_a-2.0E0*2*I_KINETIC_Hx2y2z_D2z_a+2*1*I_KINETIC_Fx2y_D2z;
  abcd[365] = 4.0E0*I_KINETIC_K3xy3z_D2z_aa-2.0E0*2*I_KINETIC_H3xyz_D2z_a-2.0E0*2*I_KINETIC_Hxy3z_D2z_a+2*2*I_KINETIC_Fxyz_D2z;
  abcd[366] = 4.0E0*I_KINETIC_K3x4z_D2z_aa-2.0E0*3*I_KINETIC_H3x2z_D2z_a-2.0E0*2*I_KINETIC_Hx4z_D2z_a+2*3*I_KINETIC_Fx2z_D2z;
  abcd[367] = 4.0E0*I_KINETIC_K2x4yz_D2z_aa-2.0E0*1*I_KINETIC_H4yz_D2z_a;
  abcd[368] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*1*I_KINETIC_H2x3y_D2z_a-2.0E0*1*I_KINETIC_H3y2z_D2z_a+1*I_KINETIC_F3y_D2z;
  abcd[369] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*2*I_KINETIC_H2x2yz_D2z_a-2.0E0*1*I_KINETIC_H2y3z_D2z_a+2*I_KINETIC_F2yz_D2z;
  abcd[370] = 4.0E0*I_KINETIC_K2xy4z_D2z_aa-2.0E0*3*I_KINETIC_H2xy2z_D2z_a-2.0E0*1*I_KINETIC_Hy4z_D2z_a+3*I_KINETIC_Fy2z_D2z;
  abcd[371] = 4.0E0*I_KINETIC_K2x5z_D2z_aa-2.0E0*4*I_KINETIC_H2x3z_D2z_a-2.0E0*1*I_KINETIC_H5z_D2z_a+4*I_KINETIC_F3z_D2z;
  abcd[372] = 4.0E0*I_KINETIC_Kx5yz_D2z_aa;
  abcd[373] = 4.0E0*I_KINETIC_Kx4y2z_D2z_aa-2.0E0*1*I_KINETIC_Hx4y_D2z_a;
  abcd[374] = 4.0E0*I_KINETIC_Kx3y3z_D2z_aa-2.0E0*2*I_KINETIC_Hx3yz_D2z_a;
  abcd[375] = 4.0E0*I_KINETIC_Kx2y4z_D2z_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a;
  abcd[376] = 4.0E0*I_KINETIC_Kxy5z_D2z_aa-2.0E0*4*I_KINETIC_Hxy3z_D2z_a;
  abcd[377] = 4.0E0*I_KINETIC_Kx6z_D2z_aa-2.0E0*5*I_KINETIC_Hx4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[378] = 4.0E0*I_KINETIC_K5x2y_D2x_aa-2.0E0*1*I_KINETIC_H5x_D2x_a;
  abcd[379] = 4.0E0*I_KINETIC_K4x3y_D2x_aa-2.0E0*1*I_KINETIC_H4xy_D2x_a-2.0E0*2*I_KINETIC_H4xy_D2x_a;
  abcd[380] = 4.0E0*I_KINETIC_K4x2yz_D2x_aa-2.0E0*1*I_KINETIC_H4xz_D2x_a;
  abcd[381] = 4.0E0*I_KINETIC_K3x4y_D2x_aa-2.0E0*2*I_KINETIC_H3x2y_D2x_a-2.0E0*3*I_KINETIC_H3x2y_D2x_a+2*1*I_KINETIC_F3x_D2x;
  abcd[382] = 4.0E0*I_KINETIC_K3x3yz_D2x_aa-2.0E0*1*I_KINETIC_H3xyz_D2x_a-2.0E0*2*I_KINETIC_H3xyz_D2x_a;
  abcd[383] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_H3x2z_D2x_a;
  abcd[384] = 4.0E0*I_KINETIC_K2x5y_D2x_aa-2.0E0*3*I_KINETIC_H2x3y_D2x_a-2.0E0*4*I_KINETIC_H2x3y_D2x_a+3*2*I_KINETIC_F2xy_D2x;
  abcd[385] = 4.0E0*I_KINETIC_K2x4yz_D2x_aa-2.0E0*2*I_KINETIC_H2x2yz_D2x_a-2.0E0*3*I_KINETIC_H2x2yz_D2x_a+2*1*I_KINETIC_F2xz_D2x;
  abcd[386] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*1*I_KINETIC_H2xy2z_D2x_a-2.0E0*2*I_KINETIC_H2xy2z_D2x_a;
  abcd[387] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*1*I_KINETIC_H2x3z_D2x_a;
  abcd[388] = 4.0E0*I_KINETIC_Kx6y_D2x_aa-2.0E0*4*I_KINETIC_Hx4y_D2x_a-2.0E0*5*I_KINETIC_Hx4y_D2x_a+4*3*I_KINETIC_Fx2y_D2x;
  abcd[389] = 4.0E0*I_KINETIC_Kx5yz_D2x_aa-2.0E0*3*I_KINETIC_Hx3yz_D2x_a-2.0E0*4*I_KINETIC_Hx3yz_D2x_a+3*2*I_KINETIC_Fxyz_D2x;
  abcd[390] = 4.0E0*I_KINETIC_Kx4y2z_D2x_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2x_a-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a+2*1*I_KINETIC_Fx2z_D2x;
  abcd[391] = 4.0E0*I_KINETIC_Kx3y3z_D2x_aa-2.0E0*1*I_KINETIC_Hxy3z_D2x_a-2.0E0*2*I_KINETIC_Hxy3z_D2x_a;
  abcd[392] = 4.0E0*I_KINETIC_Kx2y4z_D2x_aa-2.0E0*1*I_KINETIC_Hx4z_D2x_a;
  abcd[393] = 4.0E0*I_KINETIC_K7y_D2x_aa-2.0E0*5*I_KINETIC_H5y_D2x_a-2.0E0*6*I_KINETIC_H5y_D2x_a+5*4*I_KINETIC_F3y_D2x;
  abcd[394] = 4.0E0*I_KINETIC_K6yz_D2x_aa-2.0E0*4*I_KINETIC_H4yz_D2x_a-2.0E0*5*I_KINETIC_H4yz_D2x_a+4*3*I_KINETIC_F2yz_D2x;
  abcd[395] = 4.0E0*I_KINETIC_K5y2z_D2x_aa-2.0E0*3*I_KINETIC_H3y2z_D2x_a-2.0E0*4*I_KINETIC_H3y2z_D2x_a+3*2*I_KINETIC_Fy2z_D2x;
  abcd[396] = 4.0E0*I_KINETIC_K4y3z_D2x_aa-2.0E0*2*I_KINETIC_H2y3z_D2x_a-2.0E0*3*I_KINETIC_H2y3z_D2x_a+2*1*I_KINETIC_F3z_D2x;
  abcd[397] = 4.0E0*I_KINETIC_K3y4z_D2x_aa-2.0E0*1*I_KINETIC_Hy4z_D2x_a-2.0E0*2*I_KINETIC_Hy4z_D2x_a;
  abcd[398] = 4.0E0*I_KINETIC_K2y5z_D2x_aa-2.0E0*1*I_KINETIC_H5z_D2x_a;
  abcd[399] = 4.0E0*I_KINETIC_K5x2y_Dxy_aa-2.0E0*1*I_KINETIC_H5x_Dxy_a;
  abcd[400] = 4.0E0*I_KINETIC_K4x3y_Dxy_aa-2.0E0*1*I_KINETIC_H4xy_Dxy_a-2.0E0*2*I_KINETIC_H4xy_Dxy_a;
  abcd[401] = 4.0E0*I_KINETIC_K4x2yz_Dxy_aa-2.0E0*1*I_KINETIC_H4xz_Dxy_a;
  abcd[402] = 4.0E0*I_KINETIC_K3x4y_Dxy_aa-2.0E0*2*I_KINETIC_H3x2y_Dxy_a-2.0E0*3*I_KINETIC_H3x2y_Dxy_a+2*1*I_KINETIC_F3x_Dxy;
  abcd[403] = 4.0E0*I_KINETIC_K3x3yz_Dxy_aa-2.0E0*1*I_KINETIC_H3xyz_Dxy_a-2.0E0*2*I_KINETIC_H3xyz_Dxy_a;
  abcd[404] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3x2z_Dxy_a;
  abcd[405] = 4.0E0*I_KINETIC_K2x5y_Dxy_aa-2.0E0*3*I_KINETIC_H2x3y_Dxy_a-2.0E0*4*I_KINETIC_H2x3y_Dxy_a+3*2*I_KINETIC_F2xy_Dxy;
  abcd[406] = 4.0E0*I_KINETIC_K2x4yz_Dxy_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxy_a-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a+2*1*I_KINETIC_F2xz_Dxy;
  abcd[407] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*1*I_KINETIC_H2xy2z_Dxy_a-2.0E0*2*I_KINETIC_H2xy2z_Dxy_a;
  abcd[408] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*1*I_KINETIC_H2x3z_Dxy_a;
  abcd[409] = 4.0E0*I_KINETIC_Kx6y_Dxy_aa-2.0E0*4*I_KINETIC_Hx4y_Dxy_a-2.0E0*5*I_KINETIC_Hx4y_Dxy_a+4*3*I_KINETIC_Fx2y_Dxy;
  abcd[410] = 4.0E0*I_KINETIC_Kx5yz_Dxy_aa-2.0E0*3*I_KINETIC_Hx3yz_Dxy_a-2.0E0*4*I_KINETIC_Hx3yz_Dxy_a+3*2*I_KINETIC_Fxyz_Dxy;
  abcd[411] = 4.0E0*I_KINETIC_Kx4y2z_Dxy_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dxy_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a+2*1*I_KINETIC_Fx2z_Dxy;
  abcd[412] = 4.0E0*I_KINETIC_Kx3y3z_Dxy_aa-2.0E0*1*I_KINETIC_Hxy3z_Dxy_a-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a;
  abcd[413] = 4.0E0*I_KINETIC_Kx2y4z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4z_Dxy_a;
  abcd[414] = 4.0E0*I_KINETIC_K7y_Dxy_aa-2.0E0*5*I_KINETIC_H5y_Dxy_a-2.0E0*6*I_KINETIC_H5y_Dxy_a+5*4*I_KINETIC_F3y_Dxy;
  abcd[415] = 4.0E0*I_KINETIC_K6yz_Dxy_aa-2.0E0*4*I_KINETIC_H4yz_Dxy_a-2.0E0*5*I_KINETIC_H4yz_Dxy_a+4*3*I_KINETIC_F2yz_Dxy;
  abcd[416] = 4.0E0*I_KINETIC_K5y2z_Dxy_aa-2.0E0*3*I_KINETIC_H3y2z_Dxy_a-2.0E0*4*I_KINETIC_H3y2z_Dxy_a+3*2*I_KINETIC_Fy2z_Dxy;
  abcd[417] = 4.0E0*I_KINETIC_K4y3z_Dxy_aa-2.0E0*2*I_KINETIC_H2y3z_Dxy_a-2.0E0*3*I_KINETIC_H2y3z_Dxy_a+2*1*I_KINETIC_F3z_Dxy;
  abcd[418] = 4.0E0*I_KINETIC_K3y4z_Dxy_aa-2.0E0*1*I_KINETIC_Hy4z_Dxy_a-2.0E0*2*I_KINETIC_Hy4z_Dxy_a;
  abcd[419] = 4.0E0*I_KINETIC_K2y5z_Dxy_aa-2.0E0*1*I_KINETIC_H5z_Dxy_a;
  abcd[420] = 4.0E0*I_KINETIC_K5x2y_Dxz_aa-2.0E0*1*I_KINETIC_H5x_Dxz_a;
  abcd[421] = 4.0E0*I_KINETIC_K4x3y_Dxz_aa-2.0E0*1*I_KINETIC_H4xy_Dxz_a-2.0E0*2*I_KINETIC_H4xy_Dxz_a;
  abcd[422] = 4.0E0*I_KINETIC_K4x2yz_Dxz_aa-2.0E0*1*I_KINETIC_H4xz_Dxz_a;
  abcd[423] = 4.0E0*I_KINETIC_K3x4y_Dxz_aa-2.0E0*2*I_KINETIC_H3x2y_Dxz_a-2.0E0*3*I_KINETIC_H3x2y_Dxz_a+2*1*I_KINETIC_F3x_Dxz;
  abcd[424] = 4.0E0*I_KINETIC_K3x3yz_Dxz_aa-2.0E0*1*I_KINETIC_H3xyz_Dxz_a-2.0E0*2*I_KINETIC_H3xyz_Dxz_a;
  abcd[425] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3x2z_Dxz_a;
  abcd[426] = 4.0E0*I_KINETIC_K2x5y_Dxz_aa-2.0E0*3*I_KINETIC_H2x3y_Dxz_a-2.0E0*4*I_KINETIC_H2x3y_Dxz_a+3*2*I_KINETIC_F2xy_Dxz;
  abcd[427] = 4.0E0*I_KINETIC_K2x4yz_Dxz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxz_a-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a+2*1*I_KINETIC_F2xz_Dxz;
  abcd[428] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*1*I_KINETIC_H2xy2z_Dxz_a-2.0E0*2*I_KINETIC_H2xy2z_Dxz_a;
  abcd[429] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*1*I_KINETIC_H2x3z_Dxz_a;
  abcd[430] = 4.0E0*I_KINETIC_Kx6y_Dxz_aa-2.0E0*4*I_KINETIC_Hx4y_Dxz_a-2.0E0*5*I_KINETIC_Hx4y_Dxz_a+4*3*I_KINETIC_Fx2y_Dxz;
  abcd[431] = 4.0E0*I_KINETIC_Kx5yz_Dxz_aa-2.0E0*3*I_KINETIC_Hx3yz_Dxz_a-2.0E0*4*I_KINETIC_Hx3yz_Dxz_a+3*2*I_KINETIC_Fxyz_Dxz;
  abcd[432] = 4.0E0*I_KINETIC_Kx4y2z_Dxz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dxz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a+2*1*I_KINETIC_Fx2z_Dxz;
  abcd[433] = 4.0E0*I_KINETIC_Kx3y3z_Dxz_aa-2.0E0*1*I_KINETIC_Hxy3z_Dxz_a-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a;
  abcd[434] = 4.0E0*I_KINETIC_Kx2y4z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4z_Dxz_a;
  abcd[435] = 4.0E0*I_KINETIC_K7y_Dxz_aa-2.0E0*5*I_KINETIC_H5y_Dxz_a-2.0E0*6*I_KINETIC_H5y_Dxz_a+5*4*I_KINETIC_F3y_Dxz;
  abcd[436] = 4.0E0*I_KINETIC_K6yz_Dxz_aa-2.0E0*4*I_KINETIC_H4yz_Dxz_a-2.0E0*5*I_KINETIC_H4yz_Dxz_a+4*3*I_KINETIC_F2yz_Dxz;
  abcd[437] = 4.0E0*I_KINETIC_K5y2z_Dxz_aa-2.0E0*3*I_KINETIC_H3y2z_Dxz_a-2.0E0*4*I_KINETIC_H3y2z_Dxz_a+3*2*I_KINETIC_Fy2z_Dxz;
  abcd[438] = 4.0E0*I_KINETIC_K4y3z_Dxz_aa-2.0E0*2*I_KINETIC_H2y3z_Dxz_a-2.0E0*3*I_KINETIC_H2y3z_Dxz_a+2*1*I_KINETIC_F3z_Dxz;
  abcd[439] = 4.0E0*I_KINETIC_K3y4z_Dxz_aa-2.0E0*1*I_KINETIC_Hy4z_Dxz_a-2.0E0*2*I_KINETIC_Hy4z_Dxz_a;
  abcd[440] = 4.0E0*I_KINETIC_K2y5z_Dxz_aa-2.0E0*1*I_KINETIC_H5z_Dxz_a;
  abcd[441] = 4.0E0*I_KINETIC_K5x2y_D2y_aa-2.0E0*1*I_KINETIC_H5x_D2y_a;
  abcd[442] = 4.0E0*I_KINETIC_K4x3y_D2y_aa-2.0E0*1*I_KINETIC_H4xy_D2y_a-2.0E0*2*I_KINETIC_H4xy_D2y_a;
  abcd[443] = 4.0E0*I_KINETIC_K4x2yz_D2y_aa-2.0E0*1*I_KINETIC_H4xz_D2y_a;
  abcd[444] = 4.0E0*I_KINETIC_K3x4y_D2y_aa-2.0E0*2*I_KINETIC_H3x2y_D2y_a-2.0E0*3*I_KINETIC_H3x2y_D2y_a+2*1*I_KINETIC_F3x_D2y;
  abcd[445] = 4.0E0*I_KINETIC_K3x3yz_D2y_aa-2.0E0*1*I_KINETIC_H3xyz_D2y_a-2.0E0*2*I_KINETIC_H3xyz_D2y_a;
  abcd[446] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_H3x2z_D2y_a;
  abcd[447] = 4.0E0*I_KINETIC_K2x5y_D2y_aa-2.0E0*3*I_KINETIC_H2x3y_D2y_a-2.0E0*4*I_KINETIC_H2x3y_D2y_a+3*2*I_KINETIC_F2xy_D2y;
  abcd[448] = 4.0E0*I_KINETIC_K2x4yz_D2y_aa-2.0E0*2*I_KINETIC_H2x2yz_D2y_a-2.0E0*3*I_KINETIC_H2x2yz_D2y_a+2*1*I_KINETIC_F2xz_D2y;
  abcd[449] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*1*I_KINETIC_H2xy2z_D2y_a-2.0E0*2*I_KINETIC_H2xy2z_D2y_a;
  abcd[450] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*1*I_KINETIC_H2x3z_D2y_a;
  abcd[451] = 4.0E0*I_KINETIC_Kx6y_D2y_aa-2.0E0*4*I_KINETIC_Hx4y_D2y_a-2.0E0*5*I_KINETIC_Hx4y_D2y_a+4*3*I_KINETIC_Fx2y_D2y;
  abcd[452] = 4.0E0*I_KINETIC_Kx5yz_D2y_aa-2.0E0*3*I_KINETIC_Hx3yz_D2y_a-2.0E0*4*I_KINETIC_Hx3yz_D2y_a+3*2*I_KINETIC_Fxyz_D2y;
  abcd[453] = 4.0E0*I_KINETIC_Kx4y2z_D2y_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2y_a-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a+2*1*I_KINETIC_Fx2z_D2y;
  abcd[454] = 4.0E0*I_KINETIC_Kx3y3z_D2y_aa-2.0E0*1*I_KINETIC_Hxy3z_D2y_a-2.0E0*2*I_KINETIC_Hxy3z_D2y_a;
  abcd[455] = 4.0E0*I_KINETIC_Kx2y4z_D2y_aa-2.0E0*1*I_KINETIC_Hx4z_D2y_a;
  abcd[456] = 4.0E0*I_KINETIC_K7y_D2y_aa-2.0E0*5*I_KINETIC_H5y_D2y_a-2.0E0*6*I_KINETIC_H5y_D2y_a+5*4*I_KINETIC_F3y_D2y;
  abcd[457] = 4.0E0*I_KINETIC_K6yz_D2y_aa-2.0E0*4*I_KINETIC_H4yz_D2y_a-2.0E0*5*I_KINETIC_H4yz_D2y_a+4*3*I_KINETIC_F2yz_D2y;
  abcd[458] = 4.0E0*I_KINETIC_K5y2z_D2y_aa-2.0E0*3*I_KINETIC_H3y2z_D2y_a-2.0E0*4*I_KINETIC_H3y2z_D2y_a+3*2*I_KINETIC_Fy2z_D2y;
  abcd[459] = 4.0E0*I_KINETIC_K4y3z_D2y_aa-2.0E0*2*I_KINETIC_H2y3z_D2y_a-2.0E0*3*I_KINETIC_H2y3z_D2y_a+2*1*I_KINETIC_F3z_D2y;
  abcd[460] = 4.0E0*I_KINETIC_K3y4z_D2y_aa-2.0E0*1*I_KINETIC_Hy4z_D2y_a-2.0E0*2*I_KINETIC_Hy4z_D2y_a;
  abcd[461] = 4.0E0*I_KINETIC_K2y5z_D2y_aa-2.0E0*1*I_KINETIC_H5z_D2y_a;
  abcd[462] = 4.0E0*I_KINETIC_K5x2y_Dyz_aa-2.0E0*1*I_KINETIC_H5x_Dyz_a;
  abcd[463] = 4.0E0*I_KINETIC_K4x3y_Dyz_aa-2.0E0*1*I_KINETIC_H4xy_Dyz_a-2.0E0*2*I_KINETIC_H4xy_Dyz_a;
  abcd[464] = 4.0E0*I_KINETIC_K4x2yz_Dyz_aa-2.0E0*1*I_KINETIC_H4xz_Dyz_a;
  abcd[465] = 4.0E0*I_KINETIC_K3x4y_Dyz_aa-2.0E0*2*I_KINETIC_H3x2y_Dyz_a-2.0E0*3*I_KINETIC_H3x2y_Dyz_a+2*1*I_KINETIC_F3x_Dyz;
  abcd[466] = 4.0E0*I_KINETIC_K3x3yz_Dyz_aa-2.0E0*1*I_KINETIC_H3xyz_Dyz_a-2.0E0*2*I_KINETIC_H3xyz_Dyz_a;
  abcd[467] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3x2z_Dyz_a;
  abcd[468] = 4.0E0*I_KINETIC_K2x5y_Dyz_aa-2.0E0*3*I_KINETIC_H2x3y_Dyz_a-2.0E0*4*I_KINETIC_H2x3y_Dyz_a+3*2*I_KINETIC_F2xy_Dyz;
  abcd[469] = 4.0E0*I_KINETIC_K2x4yz_Dyz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dyz_a-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a+2*1*I_KINETIC_F2xz_Dyz;
  abcd[470] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*1*I_KINETIC_H2xy2z_Dyz_a-2.0E0*2*I_KINETIC_H2xy2z_Dyz_a;
  abcd[471] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*1*I_KINETIC_H2x3z_Dyz_a;
  abcd[472] = 4.0E0*I_KINETIC_Kx6y_Dyz_aa-2.0E0*4*I_KINETIC_Hx4y_Dyz_a-2.0E0*5*I_KINETIC_Hx4y_Dyz_a+4*3*I_KINETIC_Fx2y_Dyz;
  abcd[473] = 4.0E0*I_KINETIC_Kx5yz_Dyz_aa-2.0E0*3*I_KINETIC_Hx3yz_Dyz_a-2.0E0*4*I_KINETIC_Hx3yz_Dyz_a+3*2*I_KINETIC_Fxyz_Dyz;
  abcd[474] = 4.0E0*I_KINETIC_Kx4y2z_Dyz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dyz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a+2*1*I_KINETIC_Fx2z_Dyz;
  abcd[475] = 4.0E0*I_KINETIC_Kx3y3z_Dyz_aa-2.0E0*1*I_KINETIC_Hxy3z_Dyz_a-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a;
  abcd[476] = 4.0E0*I_KINETIC_Kx2y4z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4z_Dyz_a;
  abcd[477] = 4.0E0*I_KINETIC_K7y_Dyz_aa-2.0E0*5*I_KINETIC_H5y_Dyz_a-2.0E0*6*I_KINETIC_H5y_Dyz_a+5*4*I_KINETIC_F3y_Dyz;
  abcd[478] = 4.0E0*I_KINETIC_K6yz_Dyz_aa-2.0E0*4*I_KINETIC_H4yz_Dyz_a-2.0E0*5*I_KINETIC_H4yz_Dyz_a+4*3*I_KINETIC_F2yz_Dyz;
  abcd[479] = 4.0E0*I_KINETIC_K5y2z_Dyz_aa-2.0E0*3*I_KINETIC_H3y2z_Dyz_a-2.0E0*4*I_KINETIC_H3y2z_Dyz_a+3*2*I_KINETIC_Fy2z_Dyz;
  abcd[480] = 4.0E0*I_KINETIC_K4y3z_Dyz_aa-2.0E0*2*I_KINETIC_H2y3z_Dyz_a-2.0E0*3*I_KINETIC_H2y3z_Dyz_a+2*1*I_KINETIC_F3z_Dyz;
  abcd[481] = 4.0E0*I_KINETIC_K3y4z_Dyz_aa-2.0E0*1*I_KINETIC_Hy4z_Dyz_a-2.0E0*2*I_KINETIC_Hy4z_Dyz_a;
  abcd[482] = 4.0E0*I_KINETIC_K2y5z_Dyz_aa-2.0E0*1*I_KINETIC_H5z_Dyz_a;
  abcd[483] = 4.0E0*I_KINETIC_K5x2y_D2z_aa-2.0E0*1*I_KINETIC_H5x_D2z_a;
  abcd[484] = 4.0E0*I_KINETIC_K4x3y_D2z_aa-2.0E0*1*I_KINETIC_H4xy_D2z_a-2.0E0*2*I_KINETIC_H4xy_D2z_a;
  abcd[485] = 4.0E0*I_KINETIC_K4x2yz_D2z_aa-2.0E0*1*I_KINETIC_H4xz_D2z_a;
  abcd[486] = 4.0E0*I_KINETIC_K3x4y_D2z_aa-2.0E0*2*I_KINETIC_H3x2y_D2z_a-2.0E0*3*I_KINETIC_H3x2y_D2z_a+2*1*I_KINETIC_F3x_D2z;
  abcd[487] = 4.0E0*I_KINETIC_K3x3yz_D2z_aa-2.0E0*1*I_KINETIC_H3xyz_D2z_a-2.0E0*2*I_KINETIC_H3xyz_D2z_a;
  abcd[488] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_H3x2z_D2z_a;
  abcd[489] = 4.0E0*I_KINETIC_K2x5y_D2z_aa-2.0E0*3*I_KINETIC_H2x3y_D2z_a-2.0E0*4*I_KINETIC_H2x3y_D2z_a+3*2*I_KINETIC_F2xy_D2z;
  abcd[490] = 4.0E0*I_KINETIC_K2x4yz_D2z_aa-2.0E0*2*I_KINETIC_H2x2yz_D2z_a-2.0E0*3*I_KINETIC_H2x2yz_D2z_a+2*1*I_KINETIC_F2xz_D2z;
  abcd[491] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*1*I_KINETIC_H2xy2z_D2z_a-2.0E0*2*I_KINETIC_H2xy2z_D2z_a;
  abcd[492] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*1*I_KINETIC_H2x3z_D2z_a;
  abcd[493] = 4.0E0*I_KINETIC_Kx6y_D2z_aa-2.0E0*4*I_KINETIC_Hx4y_D2z_a-2.0E0*5*I_KINETIC_Hx4y_D2z_a+4*3*I_KINETIC_Fx2y_D2z;
  abcd[494] = 4.0E0*I_KINETIC_Kx5yz_D2z_aa-2.0E0*3*I_KINETIC_Hx3yz_D2z_a-2.0E0*4*I_KINETIC_Hx3yz_D2z_a+3*2*I_KINETIC_Fxyz_D2z;
  abcd[495] = 4.0E0*I_KINETIC_Kx4y2z_D2z_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2z_a-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a+2*1*I_KINETIC_Fx2z_D2z;
  abcd[496] = 4.0E0*I_KINETIC_Kx3y3z_D2z_aa-2.0E0*1*I_KINETIC_Hxy3z_D2z_a-2.0E0*2*I_KINETIC_Hxy3z_D2z_a;
  abcd[497] = 4.0E0*I_KINETIC_Kx2y4z_D2z_aa-2.0E0*1*I_KINETIC_Hx4z_D2z_a;
  abcd[498] = 4.0E0*I_KINETIC_K7y_D2z_aa-2.0E0*5*I_KINETIC_H5y_D2z_a-2.0E0*6*I_KINETIC_H5y_D2z_a+5*4*I_KINETIC_F3y_D2z;
  abcd[499] = 4.0E0*I_KINETIC_K6yz_D2z_aa-2.0E0*4*I_KINETIC_H4yz_D2z_a-2.0E0*5*I_KINETIC_H4yz_D2z_a+4*3*I_KINETIC_F2yz_D2z;
  abcd[500] = 4.0E0*I_KINETIC_K5y2z_D2z_aa-2.0E0*3*I_KINETIC_H3y2z_D2z_a-2.0E0*4*I_KINETIC_H3y2z_D2z_a+3*2*I_KINETIC_Fy2z_D2z;
  abcd[501] = 4.0E0*I_KINETIC_K4y3z_D2z_aa-2.0E0*2*I_KINETIC_H2y3z_D2z_a-2.0E0*3*I_KINETIC_H2y3z_D2z_a+2*1*I_KINETIC_F3z_D2z;
  abcd[502] = 4.0E0*I_KINETIC_K3y4z_D2z_aa-2.0E0*1*I_KINETIC_Hy4z_D2z_a-2.0E0*2*I_KINETIC_Hy4z_D2z_a;
  abcd[503] = 4.0E0*I_KINETIC_K2y5z_D2z_aa-2.0E0*1*I_KINETIC_H5z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[504] = 4.0E0*I_KINETIC_K5xyz_D2x_aa;
  abcd[505] = 4.0E0*I_KINETIC_K4x2yz_D2x_aa-2.0E0*1*I_KINETIC_H4xz_D2x_a;
  abcd[506] = 4.0E0*I_KINETIC_K4xy2z_D2x_aa-2.0E0*1*I_KINETIC_H4xy_D2x_a;
  abcd[507] = 4.0E0*I_KINETIC_K3x3yz_D2x_aa-2.0E0*2*I_KINETIC_H3xyz_D2x_a;
  abcd[508] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_H3x2y_D2x_a-2.0E0*1*I_KINETIC_H3x2z_D2x_a+1*I_KINETIC_F3x_D2x;
  abcd[509] = 4.0E0*I_KINETIC_K3xy3z_D2x_aa-2.0E0*2*I_KINETIC_H3xyz_D2x_a;
  abcd[510] = 4.0E0*I_KINETIC_K2x4yz_D2x_aa-2.0E0*3*I_KINETIC_H2x2yz_D2x_a;
  abcd[511] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*1*I_KINETIC_H2x3y_D2x_a-2.0E0*2*I_KINETIC_H2xy2z_D2x_a+2*1*I_KINETIC_F2xy_D2x;
  abcd[512] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*2*I_KINETIC_H2x2yz_D2x_a-2.0E0*1*I_KINETIC_H2x3z_D2x_a+2*I_KINETIC_F2xz_D2x;
  abcd[513] = 4.0E0*I_KINETIC_K2xy4z_D2x_aa-2.0E0*3*I_KINETIC_H2xy2z_D2x_a;
  abcd[514] = 4.0E0*I_KINETIC_Kx5yz_D2x_aa-2.0E0*4*I_KINETIC_Hx3yz_D2x_a;
  abcd[515] = 4.0E0*I_KINETIC_Kx4y2z_D2x_aa-2.0E0*1*I_KINETIC_Hx4y_D2x_a-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a+3*1*I_KINETIC_Fx2y_D2x;
  abcd[516] = 4.0E0*I_KINETIC_Kx3y3z_D2x_aa-2.0E0*2*I_KINETIC_Hx3yz_D2x_a-2.0E0*2*I_KINETIC_Hxy3z_D2x_a+2*2*I_KINETIC_Fxyz_D2x;
  abcd[517] = 4.0E0*I_KINETIC_Kx2y4z_D2x_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a-2.0E0*1*I_KINETIC_Hx4z_D2x_a+3*I_KINETIC_Fx2z_D2x;
  abcd[518] = 4.0E0*I_KINETIC_Kxy5z_D2x_aa-2.0E0*4*I_KINETIC_Hxy3z_D2x_a;
  abcd[519] = 4.0E0*I_KINETIC_K6yz_D2x_aa-2.0E0*5*I_KINETIC_H4yz_D2x_a;
  abcd[520] = 4.0E0*I_KINETIC_K5y2z_D2x_aa-2.0E0*1*I_KINETIC_H5y_D2x_a-2.0E0*4*I_KINETIC_H3y2z_D2x_a+4*1*I_KINETIC_F3y_D2x;
  abcd[521] = 4.0E0*I_KINETIC_K4y3z_D2x_aa-2.0E0*2*I_KINETIC_H4yz_D2x_a-2.0E0*3*I_KINETIC_H2y3z_D2x_a+3*2*I_KINETIC_F2yz_D2x;
  abcd[522] = 4.0E0*I_KINETIC_K3y4z_D2x_aa-2.0E0*3*I_KINETIC_H3y2z_D2x_a-2.0E0*2*I_KINETIC_Hy4z_D2x_a+2*3*I_KINETIC_Fy2z_D2x;
  abcd[523] = 4.0E0*I_KINETIC_K2y5z_D2x_aa-2.0E0*4*I_KINETIC_H2y3z_D2x_a-2.0E0*1*I_KINETIC_H5z_D2x_a+4*I_KINETIC_F3z_D2x;
  abcd[524] = 4.0E0*I_KINETIC_Ky6z_D2x_aa-2.0E0*5*I_KINETIC_Hy4z_D2x_a;
  abcd[525] = 4.0E0*I_KINETIC_K5xyz_Dxy_aa;
  abcd[526] = 4.0E0*I_KINETIC_K4x2yz_Dxy_aa-2.0E0*1*I_KINETIC_H4xz_Dxy_a;
  abcd[527] = 4.0E0*I_KINETIC_K4xy2z_Dxy_aa-2.0E0*1*I_KINETIC_H4xy_Dxy_a;
  abcd[528] = 4.0E0*I_KINETIC_K3x3yz_Dxy_aa-2.0E0*2*I_KINETIC_H3xyz_Dxy_a;
  abcd[529] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3x2y_Dxy_a-2.0E0*1*I_KINETIC_H3x2z_Dxy_a+1*I_KINETIC_F3x_Dxy;
  abcd[530] = 4.0E0*I_KINETIC_K3xy3z_Dxy_aa-2.0E0*2*I_KINETIC_H3xyz_Dxy_a;
  abcd[531] = 4.0E0*I_KINETIC_K2x4yz_Dxy_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxy_a;
  abcd[532] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*1*I_KINETIC_H2x3y_Dxy_a-2.0E0*2*I_KINETIC_H2xy2z_Dxy_a+2*1*I_KINETIC_F2xy_Dxy;
  abcd[533] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxy_a-2.0E0*1*I_KINETIC_H2x3z_Dxy_a+2*I_KINETIC_F2xz_Dxy;
  abcd[534] = 4.0E0*I_KINETIC_K2xy4z_Dxy_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a;
  abcd[535] = 4.0E0*I_KINETIC_Kx5yz_Dxy_aa-2.0E0*4*I_KINETIC_Hx3yz_Dxy_a;
  abcd[536] = 4.0E0*I_KINETIC_Kx4y2z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4y_Dxy_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a+3*1*I_KINETIC_Fx2y_Dxy;
  abcd[537] = 4.0E0*I_KINETIC_Kx3y3z_Dxy_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a-2.0E0*2*I_KINETIC_Hxy3z_Dxy_a+2*2*I_KINETIC_Fxyz_Dxy;
  abcd[538] = 4.0E0*I_KINETIC_Kx2y4z_Dxy_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a-2.0E0*1*I_KINETIC_Hx4z_Dxy_a+3*I_KINETIC_Fx2z_Dxy;
  abcd[539] = 4.0E0*I_KINETIC_Kxy5z_Dxy_aa-2.0E0*4*I_KINETIC_Hxy3z_Dxy_a;
  abcd[540] = 4.0E0*I_KINETIC_K6yz_Dxy_aa-2.0E0*5*I_KINETIC_H4yz_Dxy_a;
  abcd[541] = 4.0E0*I_KINETIC_K5y2z_Dxy_aa-2.0E0*1*I_KINETIC_H5y_Dxy_a-2.0E0*4*I_KINETIC_H3y2z_Dxy_a+4*1*I_KINETIC_F3y_Dxy;
  abcd[542] = 4.0E0*I_KINETIC_K4y3z_Dxy_aa-2.0E0*2*I_KINETIC_H4yz_Dxy_a-2.0E0*3*I_KINETIC_H2y3z_Dxy_a+3*2*I_KINETIC_F2yz_Dxy;
  abcd[543] = 4.0E0*I_KINETIC_K3y4z_Dxy_aa-2.0E0*3*I_KINETIC_H3y2z_Dxy_a-2.0E0*2*I_KINETIC_Hy4z_Dxy_a+2*3*I_KINETIC_Fy2z_Dxy;
  abcd[544] = 4.0E0*I_KINETIC_K2y5z_Dxy_aa-2.0E0*4*I_KINETIC_H2y3z_Dxy_a-2.0E0*1*I_KINETIC_H5z_Dxy_a+4*I_KINETIC_F3z_Dxy;
  abcd[545] = 4.0E0*I_KINETIC_Ky6z_Dxy_aa-2.0E0*5*I_KINETIC_Hy4z_Dxy_a;
  abcd[546] = 4.0E0*I_KINETIC_K5xyz_Dxz_aa;
  abcd[547] = 4.0E0*I_KINETIC_K4x2yz_Dxz_aa-2.0E0*1*I_KINETIC_H4xz_Dxz_a;
  abcd[548] = 4.0E0*I_KINETIC_K4xy2z_Dxz_aa-2.0E0*1*I_KINETIC_H4xy_Dxz_a;
  abcd[549] = 4.0E0*I_KINETIC_K3x3yz_Dxz_aa-2.0E0*2*I_KINETIC_H3xyz_Dxz_a;
  abcd[550] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3x2y_Dxz_a-2.0E0*1*I_KINETIC_H3x2z_Dxz_a+1*I_KINETIC_F3x_Dxz;
  abcd[551] = 4.0E0*I_KINETIC_K3xy3z_Dxz_aa-2.0E0*2*I_KINETIC_H3xyz_Dxz_a;
  abcd[552] = 4.0E0*I_KINETIC_K2x4yz_Dxz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dxz_a;
  abcd[553] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*1*I_KINETIC_H2x3y_Dxz_a-2.0E0*2*I_KINETIC_H2xy2z_Dxz_a+2*1*I_KINETIC_F2xy_Dxz;
  abcd[554] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dxz_a-2.0E0*1*I_KINETIC_H2x3z_Dxz_a+2*I_KINETIC_F2xz_Dxz;
  abcd[555] = 4.0E0*I_KINETIC_K2xy4z_Dxz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a;
  abcd[556] = 4.0E0*I_KINETIC_Kx5yz_Dxz_aa-2.0E0*4*I_KINETIC_Hx3yz_Dxz_a;
  abcd[557] = 4.0E0*I_KINETIC_Kx4y2z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4y_Dxz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a+3*1*I_KINETIC_Fx2y_Dxz;
  abcd[558] = 4.0E0*I_KINETIC_Kx3y3z_Dxz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a-2.0E0*2*I_KINETIC_Hxy3z_Dxz_a+2*2*I_KINETIC_Fxyz_Dxz;
  abcd[559] = 4.0E0*I_KINETIC_Kx2y4z_Dxz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a-2.0E0*1*I_KINETIC_Hx4z_Dxz_a+3*I_KINETIC_Fx2z_Dxz;
  abcd[560] = 4.0E0*I_KINETIC_Kxy5z_Dxz_aa-2.0E0*4*I_KINETIC_Hxy3z_Dxz_a;
  abcd[561] = 4.0E0*I_KINETIC_K6yz_Dxz_aa-2.0E0*5*I_KINETIC_H4yz_Dxz_a;
  abcd[562] = 4.0E0*I_KINETIC_K5y2z_Dxz_aa-2.0E0*1*I_KINETIC_H5y_Dxz_a-2.0E0*4*I_KINETIC_H3y2z_Dxz_a+4*1*I_KINETIC_F3y_Dxz;
  abcd[563] = 4.0E0*I_KINETIC_K4y3z_Dxz_aa-2.0E0*2*I_KINETIC_H4yz_Dxz_a-2.0E0*3*I_KINETIC_H2y3z_Dxz_a+3*2*I_KINETIC_F2yz_Dxz;
  abcd[564] = 4.0E0*I_KINETIC_K3y4z_Dxz_aa-2.0E0*3*I_KINETIC_H3y2z_Dxz_a-2.0E0*2*I_KINETIC_Hy4z_Dxz_a+2*3*I_KINETIC_Fy2z_Dxz;
  abcd[565] = 4.0E0*I_KINETIC_K2y5z_Dxz_aa-2.0E0*4*I_KINETIC_H2y3z_Dxz_a-2.0E0*1*I_KINETIC_H5z_Dxz_a+4*I_KINETIC_F3z_Dxz;
  abcd[566] = 4.0E0*I_KINETIC_Ky6z_Dxz_aa-2.0E0*5*I_KINETIC_Hy4z_Dxz_a;
  abcd[567] = 4.0E0*I_KINETIC_K5xyz_D2y_aa;
  abcd[568] = 4.0E0*I_KINETIC_K4x2yz_D2y_aa-2.0E0*1*I_KINETIC_H4xz_D2y_a;
  abcd[569] = 4.0E0*I_KINETIC_K4xy2z_D2y_aa-2.0E0*1*I_KINETIC_H4xy_D2y_a;
  abcd[570] = 4.0E0*I_KINETIC_K3x3yz_D2y_aa-2.0E0*2*I_KINETIC_H3xyz_D2y_a;
  abcd[571] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_H3x2y_D2y_a-2.0E0*1*I_KINETIC_H3x2z_D2y_a+1*I_KINETIC_F3x_D2y;
  abcd[572] = 4.0E0*I_KINETIC_K3xy3z_D2y_aa-2.0E0*2*I_KINETIC_H3xyz_D2y_a;
  abcd[573] = 4.0E0*I_KINETIC_K2x4yz_D2y_aa-2.0E0*3*I_KINETIC_H2x2yz_D2y_a;
  abcd[574] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*1*I_KINETIC_H2x3y_D2y_a-2.0E0*2*I_KINETIC_H2xy2z_D2y_a+2*1*I_KINETIC_F2xy_D2y;
  abcd[575] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*2*I_KINETIC_H2x2yz_D2y_a-2.0E0*1*I_KINETIC_H2x3z_D2y_a+2*I_KINETIC_F2xz_D2y;
  abcd[576] = 4.0E0*I_KINETIC_K2xy4z_D2y_aa-2.0E0*3*I_KINETIC_H2xy2z_D2y_a;
  abcd[577] = 4.0E0*I_KINETIC_Kx5yz_D2y_aa-2.0E0*4*I_KINETIC_Hx3yz_D2y_a;
  abcd[578] = 4.0E0*I_KINETIC_Kx4y2z_D2y_aa-2.0E0*1*I_KINETIC_Hx4y_D2y_a-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a+3*1*I_KINETIC_Fx2y_D2y;
  abcd[579] = 4.0E0*I_KINETIC_Kx3y3z_D2y_aa-2.0E0*2*I_KINETIC_Hx3yz_D2y_a-2.0E0*2*I_KINETIC_Hxy3z_D2y_a+2*2*I_KINETIC_Fxyz_D2y;
  abcd[580] = 4.0E0*I_KINETIC_Kx2y4z_D2y_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a-2.0E0*1*I_KINETIC_Hx4z_D2y_a+3*I_KINETIC_Fx2z_D2y;
  abcd[581] = 4.0E0*I_KINETIC_Kxy5z_D2y_aa-2.0E0*4*I_KINETIC_Hxy3z_D2y_a;
  abcd[582] = 4.0E0*I_KINETIC_K6yz_D2y_aa-2.0E0*5*I_KINETIC_H4yz_D2y_a;
  abcd[583] = 4.0E0*I_KINETIC_K5y2z_D2y_aa-2.0E0*1*I_KINETIC_H5y_D2y_a-2.0E0*4*I_KINETIC_H3y2z_D2y_a+4*1*I_KINETIC_F3y_D2y;
  abcd[584] = 4.0E0*I_KINETIC_K4y3z_D2y_aa-2.0E0*2*I_KINETIC_H4yz_D2y_a-2.0E0*3*I_KINETIC_H2y3z_D2y_a+3*2*I_KINETIC_F2yz_D2y;
  abcd[585] = 4.0E0*I_KINETIC_K3y4z_D2y_aa-2.0E0*3*I_KINETIC_H3y2z_D2y_a-2.0E0*2*I_KINETIC_Hy4z_D2y_a+2*3*I_KINETIC_Fy2z_D2y;
  abcd[586] = 4.0E0*I_KINETIC_K2y5z_D2y_aa-2.0E0*4*I_KINETIC_H2y3z_D2y_a-2.0E0*1*I_KINETIC_H5z_D2y_a+4*I_KINETIC_F3z_D2y;
  abcd[587] = 4.0E0*I_KINETIC_Ky6z_D2y_aa-2.0E0*5*I_KINETIC_Hy4z_D2y_a;
  abcd[588] = 4.0E0*I_KINETIC_K5xyz_Dyz_aa;
  abcd[589] = 4.0E0*I_KINETIC_K4x2yz_Dyz_aa-2.0E0*1*I_KINETIC_H4xz_Dyz_a;
  abcd[590] = 4.0E0*I_KINETIC_K4xy2z_Dyz_aa-2.0E0*1*I_KINETIC_H4xy_Dyz_a;
  abcd[591] = 4.0E0*I_KINETIC_K3x3yz_Dyz_aa-2.0E0*2*I_KINETIC_H3xyz_Dyz_a;
  abcd[592] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3x2y_Dyz_a-2.0E0*1*I_KINETIC_H3x2z_Dyz_a+1*I_KINETIC_F3x_Dyz;
  abcd[593] = 4.0E0*I_KINETIC_K3xy3z_Dyz_aa-2.0E0*2*I_KINETIC_H3xyz_Dyz_a;
  abcd[594] = 4.0E0*I_KINETIC_K2x4yz_Dyz_aa-2.0E0*3*I_KINETIC_H2x2yz_Dyz_a;
  abcd[595] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*1*I_KINETIC_H2x3y_Dyz_a-2.0E0*2*I_KINETIC_H2xy2z_Dyz_a+2*1*I_KINETIC_F2xy_Dyz;
  abcd[596] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*2*I_KINETIC_H2x2yz_Dyz_a-2.0E0*1*I_KINETIC_H2x3z_Dyz_a+2*I_KINETIC_F2xz_Dyz;
  abcd[597] = 4.0E0*I_KINETIC_K2xy4z_Dyz_aa-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a;
  abcd[598] = 4.0E0*I_KINETIC_Kx5yz_Dyz_aa-2.0E0*4*I_KINETIC_Hx3yz_Dyz_a;
  abcd[599] = 4.0E0*I_KINETIC_Kx4y2z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4y_Dyz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a+3*1*I_KINETIC_Fx2y_Dyz;
  abcd[600] = 4.0E0*I_KINETIC_Kx3y3z_Dyz_aa-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a-2.0E0*2*I_KINETIC_Hxy3z_Dyz_a+2*2*I_KINETIC_Fxyz_Dyz;
  abcd[601] = 4.0E0*I_KINETIC_Kx2y4z_Dyz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a-2.0E0*1*I_KINETIC_Hx4z_Dyz_a+3*I_KINETIC_Fx2z_Dyz;
  abcd[602] = 4.0E0*I_KINETIC_Kxy5z_Dyz_aa-2.0E0*4*I_KINETIC_Hxy3z_Dyz_a;
  abcd[603] = 4.0E0*I_KINETIC_K6yz_Dyz_aa-2.0E0*5*I_KINETIC_H4yz_Dyz_a;
  abcd[604] = 4.0E0*I_KINETIC_K5y2z_Dyz_aa-2.0E0*1*I_KINETIC_H5y_Dyz_a-2.0E0*4*I_KINETIC_H3y2z_Dyz_a+4*1*I_KINETIC_F3y_Dyz;
  abcd[605] = 4.0E0*I_KINETIC_K4y3z_Dyz_aa-2.0E0*2*I_KINETIC_H4yz_Dyz_a-2.0E0*3*I_KINETIC_H2y3z_Dyz_a+3*2*I_KINETIC_F2yz_Dyz;
  abcd[606] = 4.0E0*I_KINETIC_K3y4z_Dyz_aa-2.0E0*3*I_KINETIC_H3y2z_Dyz_a-2.0E0*2*I_KINETIC_Hy4z_Dyz_a+2*3*I_KINETIC_Fy2z_Dyz;
  abcd[607] = 4.0E0*I_KINETIC_K2y5z_Dyz_aa-2.0E0*4*I_KINETIC_H2y3z_Dyz_a-2.0E0*1*I_KINETIC_H5z_Dyz_a+4*I_KINETIC_F3z_Dyz;
  abcd[608] = 4.0E0*I_KINETIC_Ky6z_Dyz_aa-2.0E0*5*I_KINETIC_Hy4z_Dyz_a;
  abcd[609] = 4.0E0*I_KINETIC_K5xyz_D2z_aa;
  abcd[610] = 4.0E0*I_KINETIC_K4x2yz_D2z_aa-2.0E0*1*I_KINETIC_H4xz_D2z_a;
  abcd[611] = 4.0E0*I_KINETIC_K4xy2z_D2z_aa-2.0E0*1*I_KINETIC_H4xy_D2z_a;
  abcd[612] = 4.0E0*I_KINETIC_K3x3yz_D2z_aa-2.0E0*2*I_KINETIC_H3xyz_D2z_a;
  abcd[613] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_H3x2y_D2z_a-2.0E0*1*I_KINETIC_H3x2z_D2z_a+1*I_KINETIC_F3x_D2z;
  abcd[614] = 4.0E0*I_KINETIC_K3xy3z_D2z_aa-2.0E0*2*I_KINETIC_H3xyz_D2z_a;
  abcd[615] = 4.0E0*I_KINETIC_K2x4yz_D2z_aa-2.0E0*3*I_KINETIC_H2x2yz_D2z_a;
  abcd[616] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*1*I_KINETIC_H2x3y_D2z_a-2.0E0*2*I_KINETIC_H2xy2z_D2z_a+2*1*I_KINETIC_F2xy_D2z;
  abcd[617] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*2*I_KINETIC_H2x2yz_D2z_a-2.0E0*1*I_KINETIC_H2x3z_D2z_a+2*I_KINETIC_F2xz_D2z;
  abcd[618] = 4.0E0*I_KINETIC_K2xy4z_D2z_aa-2.0E0*3*I_KINETIC_H2xy2z_D2z_a;
  abcd[619] = 4.0E0*I_KINETIC_Kx5yz_D2z_aa-2.0E0*4*I_KINETIC_Hx3yz_D2z_a;
  abcd[620] = 4.0E0*I_KINETIC_Kx4y2z_D2z_aa-2.0E0*1*I_KINETIC_Hx4y_D2z_a-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a+3*1*I_KINETIC_Fx2y_D2z;
  abcd[621] = 4.0E0*I_KINETIC_Kx3y3z_D2z_aa-2.0E0*2*I_KINETIC_Hx3yz_D2z_a-2.0E0*2*I_KINETIC_Hxy3z_D2z_a+2*2*I_KINETIC_Fxyz_D2z;
  abcd[622] = 4.0E0*I_KINETIC_Kx2y4z_D2z_aa-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a-2.0E0*1*I_KINETIC_Hx4z_D2z_a+3*I_KINETIC_Fx2z_D2z;
  abcd[623] = 4.0E0*I_KINETIC_Kxy5z_D2z_aa-2.0E0*4*I_KINETIC_Hxy3z_D2z_a;
  abcd[624] = 4.0E0*I_KINETIC_K6yz_D2z_aa-2.0E0*5*I_KINETIC_H4yz_D2z_a;
  abcd[625] = 4.0E0*I_KINETIC_K5y2z_D2z_aa-2.0E0*1*I_KINETIC_H5y_D2z_a-2.0E0*4*I_KINETIC_H3y2z_D2z_a+4*1*I_KINETIC_F3y_D2z;
  abcd[626] = 4.0E0*I_KINETIC_K4y3z_D2z_aa-2.0E0*2*I_KINETIC_H4yz_D2z_a-2.0E0*3*I_KINETIC_H2y3z_D2z_a+3*2*I_KINETIC_F2yz_D2z;
  abcd[627] = 4.0E0*I_KINETIC_K3y4z_D2z_aa-2.0E0*3*I_KINETIC_H3y2z_D2z_a-2.0E0*2*I_KINETIC_Hy4z_D2z_a+2*3*I_KINETIC_Fy2z_D2z;
  abcd[628] = 4.0E0*I_KINETIC_K2y5z_D2z_aa-2.0E0*4*I_KINETIC_H2y3z_D2z_a-2.0E0*1*I_KINETIC_H5z_D2z_a+4*I_KINETIC_F3z_D2z;
  abcd[629] = 4.0E0*I_KINETIC_Ky6z_D2z_aa-2.0E0*5*I_KINETIC_Hy4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_D_aa
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_H_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D
   ************************************************************/
  abcd[630] = 4.0E0*I_KINETIC_K5x2z_D2x_aa-2.0E0*1*I_KINETIC_H5x_D2x_a;
  abcd[631] = 4.0E0*I_KINETIC_K4xy2z_D2x_aa-2.0E0*1*I_KINETIC_H4xy_D2x_a;
  abcd[632] = 4.0E0*I_KINETIC_K4x3z_D2x_aa-2.0E0*1*I_KINETIC_H4xz_D2x_a-2.0E0*2*I_KINETIC_H4xz_D2x_a;
  abcd[633] = 4.0E0*I_KINETIC_K3x2y2z_D2x_aa-2.0E0*1*I_KINETIC_H3x2y_D2x_a;
  abcd[634] = 4.0E0*I_KINETIC_K3xy3z_D2x_aa-2.0E0*1*I_KINETIC_H3xyz_D2x_a-2.0E0*2*I_KINETIC_H3xyz_D2x_a;
  abcd[635] = 4.0E0*I_KINETIC_K3x4z_D2x_aa-2.0E0*2*I_KINETIC_H3x2z_D2x_a-2.0E0*3*I_KINETIC_H3x2z_D2x_a+2*1*I_KINETIC_F3x_D2x;
  abcd[636] = 4.0E0*I_KINETIC_K2x3y2z_D2x_aa-2.0E0*1*I_KINETIC_H2x3y_D2x_a;
  abcd[637] = 4.0E0*I_KINETIC_K2x2y3z_D2x_aa-2.0E0*1*I_KINETIC_H2x2yz_D2x_a-2.0E0*2*I_KINETIC_H2x2yz_D2x_a;
  abcd[638] = 4.0E0*I_KINETIC_K2xy4z_D2x_aa-2.0E0*2*I_KINETIC_H2xy2z_D2x_a-2.0E0*3*I_KINETIC_H2xy2z_D2x_a+2*1*I_KINETIC_F2xy_D2x;
  abcd[639] = 4.0E0*I_KINETIC_K2x5z_D2x_aa-2.0E0*3*I_KINETIC_H2x3z_D2x_a-2.0E0*4*I_KINETIC_H2x3z_D2x_a+3*2*I_KINETIC_F2xz_D2x;
  abcd[640] = 4.0E0*I_KINETIC_Kx4y2z_D2x_aa-2.0E0*1*I_KINETIC_Hx4y_D2x_a;
  abcd[641] = 4.0E0*I_KINETIC_Kx3y3z_D2x_aa-2.0E0*1*I_KINETIC_Hx3yz_D2x_a-2.0E0*2*I_KINETIC_Hx3yz_D2x_a;
  abcd[642] = 4.0E0*I_KINETIC_Kx2y4z_D2x_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2x_a-2.0E0*3*I_KINETIC_Hx2y2z_D2x_a+2*1*I_KINETIC_Fx2y_D2x;
  abcd[643] = 4.0E0*I_KINETIC_Kxy5z_D2x_aa-2.0E0*3*I_KINETIC_Hxy3z_D2x_a-2.0E0*4*I_KINETIC_Hxy3z_D2x_a+3*2*I_KINETIC_Fxyz_D2x;
  abcd[644] = 4.0E0*I_KINETIC_Kx6z_D2x_aa-2.0E0*4*I_KINETIC_Hx4z_D2x_a-2.0E0*5*I_KINETIC_Hx4z_D2x_a+4*3*I_KINETIC_Fx2z_D2x;
  abcd[645] = 4.0E0*I_KINETIC_K5y2z_D2x_aa-2.0E0*1*I_KINETIC_H5y_D2x_a;
  abcd[646] = 4.0E0*I_KINETIC_K4y3z_D2x_aa-2.0E0*1*I_KINETIC_H4yz_D2x_a-2.0E0*2*I_KINETIC_H4yz_D2x_a;
  abcd[647] = 4.0E0*I_KINETIC_K3y4z_D2x_aa-2.0E0*2*I_KINETIC_H3y2z_D2x_a-2.0E0*3*I_KINETIC_H3y2z_D2x_a+2*1*I_KINETIC_F3y_D2x;
  abcd[648] = 4.0E0*I_KINETIC_K2y5z_D2x_aa-2.0E0*3*I_KINETIC_H2y3z_D2x_a-2.0E0*4*I_KINETIC_H2y3z_D2x_a+3*2*I_KINETIC_F2yz_D2x;
  abcd[649] = 4.0E0*I_KINETIC_Ky6z_D2x_aa-2.0E0*4*I_KINETIC_Hy4z_D2x_a-2.0E0*5*I_KINETIC_Hy4z_D2x_a+4*3*I_KINETIC_Fy2z_D2x;
  abcd[650] = 4.0E0*I_KINETIC_K7z_D2x_aa-2.0E0*5*I_KINETIC_H5z_D2x_a-2.0E0*6*I_KINETIC_H5z_D2x_a+5*4*I_KINETIC_F3z_D2x;
  abcd[651] = 4.0E0*I_KINETIC_K5x2z_Dxy_aa-2.0E0*1*I_KINETIC_H5x_Dxy_a;
  abcd[652] = 4.0E0*I_KINETIC_K4xy2z_Dxy_aa-2.0E0*1*I_KINETIC_H4xy_Dxy_a;
  abcd[653] = 4.0E0*I_KINETIC_K4x3z_Dxy_aa-2.0E0*1*I_KINETIC_H4xz_Dxy_a-2.0E0*2*I_KINETIC_H4xz_Dxy_a;
  abcd[654] = 4.0E0*I_KINETIC_K3x2y2z_Dxy_aa-2.0E0*1*I_KINETIC_H3x2y_Dxy_a;
  abcd[655] = 4.0E0*I_KINETIC_K3xy3z_Dxy_aa-2.0E0*1*I_KINETIC_H3xyz_Dxy_a-2.0E0*2*I_KINETIC_H3xyz_Dxy_a;
  abcd[656] = 4.0E0*I_KINETIC_K3x4z_Dxy_aa-2.0E0*2*I_KINETIC_H3x2z_Dxy_a-2.0E0*3*I_KINETIC_H3x2z_Dxy_a+2*1*I_KINETIC_F3x_Dxy;
  abcd[657] = 4.0E0*I_KINETIC_K2x3y2z_Dxy_aa-2.0E0*1*I_KINETIC_H2x3y_Dxy_a;
  abcd[658] = 4.0E0*I_KINETIC_K2x2y3z_Dxy_aa-2.0E0*1*I_KINETIC_H2x2yz_Dxy_a-2.0E0*2*I_KINETIC_H2x2yz_Dxy_a;
  abcd[659] = 4.0E0*I_KINETIC_K2xy4z_Dxy_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxy_a-2.0E0*3*I_KINETIC_H2xy2z_Dxy_a+2*1*I_KINETIC_F2xy_Dxy;
  abcd[660] = 4.0E0*I_KINETIC_K2x5z_Dxy_aa-2.0E0*3*I_KINETIC_H2x3z_Dxy_a-2.0E0*4*I_KINETIC_H2x3z_Dxy_a+3*2*I_KINETIC_F2xz_Dxy;
  abcd[661] = 4.0E0*I_KINETIC_Kx4y2z_Dxy_aa-2.0E0*1*I_KINETIC_Hx4y_Dxy_a;
  abcd[662] = 4.0E0*I_KINETIC_Kx3y3z_Dxy_aa-2.0E0*1*I_KINETIC_Hx3yz_Dxy_a-2.0E0*2*I_KINETIC_Hx3yz_Dxy_a;
  abcd[663] = 4.0E0*I_KINETIC_Kx2y4z_Dxy_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dxy_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxy_a+2*1*I_KINETIC_Fx2y_Dxy;
  abcd[664] = 4.0E0*I_KINETIC_Kxy5z_Dxy_aa-2.0E0*3*I_KINETIC_Hxy3z_Dxy_a-2.0E0*4*I_KINETIC_Hxy3z_Dxy_a+3*2*I_KINETIC_Fxyz_Dxy;
  abcd[665] = 4.0E0*I_KINETIC_Kx6z_Dxy_aa-2.0E0*4*I_KINETIC_Hx4z_Dxy_a-2.0E0*5*I_KINETIC_Hx4z_Dxy_a+4*3*I_KINETIC_Fx2z_Dxy;
  abcd[666] = 4.0E0*I_KINETIC_K5y2z_Dxy_aa-2.0E0*1*I_KINETIC_H5y_Dxy_a;
  abcd[667] = 4.0E0*I_KINETIC_K4y3z_Dxy_aa-2.0E0*1*I_KINETIC_H4yz_Dxy_a-2.0E0*2*I_KINETIC_H4yz_Dxy_a;
  abcd[668] = 4.0E0*I_KINETIC_K3y4z_Dxy_aa-2.0E0*2*I_KINETIC_H3y2z_Dxy_a-2.0E0*3*I_KINETIC_H3y2z_Dxy_a+2*1*I_KINETIC_F3y_Dxy;
  abcd[669] = 4.0E0*I_KINETIC_K2y5z_Dxy_aa-2.0E0*3*I_KINETIC_H2y3z_Dxy_a-2.0E0*4*I_KINETIC_H2y3z_Dxy_a+3*2*I_KINETIC_F2yz_Dxy;
  abcd[670] = 4.0E0*I_KINETIC_Ky6z_Dxy_aa-2.0E0*4*I_KINETIC_Hy4z_Dxy_a-2.0E0*5*I_KINETIC_Hy4z_Dxy_a+4*3*I_KINETIC_Fy2z_Dxy;
  abcd[671] = 4.0E0*I_KINETIC_K7z_Dxy_aa-2.0E0*5*I_KINETIC_H5z_Dxy_a-2.0E0*6*I_KINETIC_H5z_Dxy_a+5*4*I_KINETIC_F3z_Dxy;
  abcd[672] = 4.0E0*I_KINETIC_K5x2z_Dxz_aa-2.0E0*1*I_KINETIC_H5x_Dxz_a;
  abcd[673] = 4.0E0*I_KINETIC_K4xy2z_Dxz_aa-2.0E0*1*I_KINETIC_H4xy_Dxz_a;
  abcd[674] = 4.0E0*I_KINETIC_K4x3z_Dxz_aa-2.0E0*1*I_KINETIC_H4xz_Dxz_a-2.0E0*2*I_KINETIC_H4xz_Dxz_a;
  abcd[675] = 4.0E0*I_KINETIC_K3x2y2z_Dxz_aa-2.0E0*1*I_KINETIC_H3x2y_Dxz_a;
  abcd[676] = 4.0E0*I_KINETIC_K3xy3z_Dxz_aa-2.0E0*1*I_KINETIC_H3xyz_Dxz_a-2.0E0*2*I_KINETIC_H3xyz_Dxz_a;
  abcd[677] = 4.0E0*I_KINETIC_K3x4z_Dxz_aa-2.0E0*2*I_KINETIC_H3x2z_Dxz_a-2.0E0*3*I_KINETIC_H3x2z_Dxz_a+2*1*I_KINETIC_F3x_Dxz;
  abcd[678] = 4.0E0*I_KINETIC_K2x3y2z_Dxz_aa-2.0E0*1*I_KINETIC_H2x3y_Dxz_a;
  abcd[679] = 4.0E0*I_KINETIC_K2x2y3z_Dxz_aa-2.0E0*1*I_KINETIC_H2x2yz_Dxz_a-2.0E0*2*I_KINETIC_H2x2yz_Dxz_a;
  abcd[680] = 4.0E0*I_KINETIC_K2xy4z_Dxz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dxz_a-2.0E0*3*I_KINETIC_H2xy2z_Dxz_a+2*1*I_KINETIC_F2xy_Dxz;
  abcd[681] = 4.0E0*I_KINETIC_K2x5z_Dxz_aa-2.0E0*3*I_KINETIC_H2x3z_Dxz_a-2.0E0*4*I_KINETIC_H2x3z_Dxz_a+3*2*I_KINETIC_F2xz_Dxz;
  abcd[682] = 4.0E0*I_KINETIC_Kx4y2z_Dxz_aa-2.0E0*1*I_KINETIC_Hx4y_Dxz_a;
  abcd[683] = 4.0E0*I_KINETIC_Kx3y3z_Dxz_aa-2.0E0*1*I_KINETIC_Hx3yz_Dxz_a-2.0E0*2*I_KINETIC_Hx3yz_Dxz_a;
  abcd[684] = 4.0E0*I_KINETIC_Kx2y4z_Dxz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dxz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dxz_a+2*1*I_KINETIC_Fx2y_Dxz;
  abcd[685] = 4.0E0*I_KINETIC_Kxy5z_Dxz_aa-2.0E0*3*I_KINETIC_Hxy3z_Dxz_a-2.0E0*4*I_KINETIC_Hxy3z_Dxz_a+3*2*I_KINETIC_Fxyz_Dxz;
  abcd[686] = 4.0E0*I_KINETIC_Kx6z_Dxz_aa-2.0E0*4*I_KINETIC_Hx4z_Dxz_a-2.0E0*5*I_KINETIC_Hx4z_Dxz_a+4*3*I_KINETIC_Fx2z_Dxz;
  abcd[687] = 4.0E0*I_KINETIC_K5y2z_Dxz_aa-2.0E0*1*I_KINETIC_H5y_Dxz_a;
  abcd[688] = 4.0E0*I_KINETIC_K4y3z_Dxz_aa-2.0E0*1*I_KINETIC_H4yz_Dxz_a-2.0E0*2*I_KINETIC_H4yz_Dxz_a;
  abcd[689] = 4.0E0*I_KINETIC_K3y4z_Dxz_aa-2.0E0*2*I_KINETIC_H3y2z_Dxz_a-2.0E0*3*I_KINETIC_H3y2z_Dxz_a+2*1*I_KINETIC_F3y_Dxz;
  abcd[690] = 4.0E0*I_KINETIC_K2y5z_Dxz_aa-2.0E0*3*I_KINETIC_H2y3z_Dxz_a-2.0E0*4*I_KINETIC_H2y3z_Dxz_a+3*2*I_KINETIC_F2yz_Dxz;
  abcd[691] = 4.0E0*I_KINETIC_Ky6z_Dxz_aa-2.0E0*4*I_KINETIC_Hy4z_Dxz_a-2.0E0*5*I_KINETIC_Hy4z_Dxz_a+4*3*I_KINETIC_Fy2z_Dxz;
  abcd[692] = 4.0E0*I_KINETIC_K7z_Dxz_aa-2.0E0*5*I_KINETIC_H5z_Dxz_a-2.0E0*6*I_KINETIC_H5z_Dxz_a+5*4*I_KINETIC_F3z_Dxz;
  abcd[693] = 4.0E0*I_KINETIC_K5x2z_D2y_aa-2.0E0*1*I_KINETIC_H5x_D2y_a;
  abcd[694] = 4.0E0*I_KINETIC_K4xy2z_D2y_aa-2.0E0*1*I_KINETIC_H4xy_D2y_a;
  abcd[695] = 4.0E0*I_KINETIC_K4x3z_D2y_aa-2.0E0*1*I_KINETIC_H4xz_D2y_a-2.0E0*2*I_KINETIC_H4xz_D2y_a;
  abcd[696] = 4.0E0*I_KINETIC_K3x2y2z_D2y_aa-2.0E0*1*I_KINETIC_H3x2y_D2y_a;
  abcd[697] = 4.0E0*I_KINETIC_K3xy3z_D2y_aa-2.0E0*1*I_KINETIC_H3xyz_D2y_a-2.0E0*2*I_KINETIC_H3xyz_D2y_a;
  abcd[698] = 4.0E0*I_KINETIC_K3x4z_D2y_aa-2.0E0*2*I_KINETIC_H3x2z_D2y_a-2.0E0*3*I_KINETIC_H3x2z_D2y_a+2*1*I_KINETIC_F3x_D2y;
  abcd[699] = 4.0E0*I_KINETIC_K2x3y2z_D2y_aa-2.0E0*1*I_KINETIC_H2x3y_D2y_a;
  abcd[700] = 4.0E0*I_KINETIC_K2x2y3z_D2y_aa-2.0E0*1*I_KINETIC_H2x2yz_D2y_a-2.0E0*2*I_KINETIC_H2x2yz_D2y_a;
  abcd[701] = 4.0E0*I_KINETIC_K2xy4z_D2y_aa-2.0E0*2*I_KINETIC_H2xy2z_D2y_a-2.0E0*3*I_KINETIC_H2xy2z_D2y_a+2*1*I_KINETIC_F2xy_D2y;
  abcd[702] = 4.0E0*I_KINETIC_K2x5z_D2y_aa-2.0E0*3*I_KINETIC_H2x3z_D2y_a-2.0E0*4*I_KINETIC_H2x3z_D2y_a+3*2*I_KINETIC_F2xz_D2y;
  abcd[703] = 4.0E0*I_KINETIC_Kx4y2z_D2y_aa-2.0E0*1*I_KINETIC_Hx4y_D2y_a;
  abcd[704] = 4.0E0*I_KINETIC_Kx3y3z_D2y_aa-2.0E0*1*I_KINETIC_Hx3yz_D2y_a-2.0E0*2*I_KINETIC_Hx3yz_D2y_a;
  abcd[705] = 4.0E0*I_KINETIC_Kx2y4z_D2y_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2y_a-2.0E0*3*I_KINETIC_Hx2y2z_D2y_a+2*1*I_KINETIC_Fx2y_D2y;
  abcd[706] = 4.0E0*I_KINETIC_Kxy5z_D2y_aa-2.0E0*3*I_KINETIC_Hxy3z_D2y_a-2.0E0*4*I_KINETIC_Hxy3z_D2y_a+3*2*I_KINETIC_Fxyz_D2y;
  abcd[707] = 4.0E0*I_KINETIC_Kx6z_D2y_aa-2.0E0*4*I_KINETIC_Hx4z_D2y_a-2.0E0*5*I_KINETIC_Hx4z_D2y_a+4*3*I_KINETIC_Fx2z_D2y;
  abcd[708] = 4.0E0*I_KINETIC_K5y2z_D2y_aa-2.0E0*1*I_KINETIC_H5y_D2y_a;
  abcd[709] = 4.0E0*I_KINETIC_K4y3z_D2y_aa-2.0E0*1*I_KINETIC_H4yz_D2y_a-2.0E0*2*I_KINETIC_H4yz_D2y_a;
  abcd[710] = 4.0E0*I_KINETIC_K3y4z_D2y_aa-2.0E0*2*I_KINETIC_H3y2z_D2y_a-2.0E0*3*I_KINETIC_H3y2z_D2y_a+2*1*I_KINETIC_F3y_D2y;
  abcd[711] = 4.0E0*I_KINETIC_K2y5z_D2y_aa-2.0E0*3*I_KINETIC_H2y3z_D2y_a-2.0E0*4*I_KINETIC_H2y3z_D2y_a+3*2*I_KINETIC_F2yz_D2y;
  abcd[712] = 4.0E0*I_KINETIC_Ky6z_D2y_aa-2.0E0*4*I_KINETIC_Hy4z_D2y_a-2.0E0*5*I_KINETIC_Hy4z_D2y_a+4*3*I_KINETIC_Fy2z_D2y;
  abcd[713] = 4.0E0*I_KINETIC_K7z_D2y_aa-2.0E0*5*I_KINETIC_H5z_D2y_a-2.0E0*6*I_KINETIC_H5z_D2y_a+5*4*I_KINETIC_F3z_D2y;
  abcd[714] = 4.0E0*I_KINETIC_K5x2z_Dyz_aa-2.0E0*1*I_KINETIC_H5x_Dyz_a;
  abcd[715] = 4.0E0*I_KINETIC_K4xy2z_Dyz_aa-2.0E0*1*I_KINETIC_H4xy_Dyz_a;
  abcd[716] = 4.0E0*I_KINETIC_K4x3z_Dyz_aa-2.0E0*1*I_KINETIC_H4xz_Dyz_a-2.0E0*2*I_KINETIC_H4xz_Dyz_a;
  abcd[717] = 4.0E0*I_KINETIC_K3x2y2z_Dyz_aa-2.0E0*1*I_KINETIC_H3x2y_Dyz_a;
  abcd[718] = 4.0E0*I_KINETIC_K3xy3z_Dyz_aa-2.0E0*1*I_KINETIC_H3xyz_Dyz_a-2.0E0*2*I_KINETIC_H3xyz_Dyz_a;
  abcd[719] = 4.0E0*I_KINETIC_K3x4z_Dyz_aa-2.0E0*2*I_KINETIC_H3x2z_Dyz_a-2.0E0*3*I_KINETIC_H3x2z_Dyz_a+2*1*I_KINETIC_F3x_Dyz;
  abcd[720] = 4.0E0*I_KINETIC_K2x3y2z_Dyz_aa-2.0E0*1*I_KINETIC_H2x3y_Dyz_a;
  abcd[721] = 4.0E0*I_KINETIC_K2x2y3z_Dyz_aa-2.0E0*1*I_KINETIC_H2x2yz_Dyz_a-2.0E0*2*I_KINETIC_H2x2yz_Dyz_a;
  abcd[722] = 4.0E0*I_KINETIC_K2xy4z_Dyz_aa-2.0E0*2*I_KINETIC_H2xy2z_Dyz_a-2.0E0*3*I_KINETIC_H2xy2z_Dyz_a+2*1*I_KINETIC_F2xy_Dyz;
  abcd[723] = 4.0E0*I_KINETIC_K2x5z_Dyz_aa-2.0E0*3*I_KINETIC_H2x3z_Dyz_a-2.0E0*4*I_KINETIC_H2x3z_Dyz_a+3*2*I_KINETIC_F2xz_Dyz;
  abcd[724] = 4.0E0*I_KINETIC_Kx4y2z_Dyz_aa-2.0E0*1*I_KINETIC_Hx4y_Dyz_a;
  abcd[725] = 4.0E0*I_KINETIC_Kx3y3z_Dyz_aa-2.0E0*1*I_KINETIC_Hx3yz_Dyz_a-2.0E0*2*I_KINETIC_Hx3yz_Dyz_a;
  abcd[726] = 4.0E0*I_KINETIC_Kx2y4z_Dyz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Dyz_a-2.0E0*3*I_KINETIC_Hx2y2z_Dyz_a+2*1*I_KINETIC_Fx2y_Dyz;
  abcd[727] = 4.0E0*I_KINETIC_Kxy5z_Dyz_aa-2.0E0*3*I_KINETIC_Hxy3z_Dyz_a-2.0E0*4*I_KINETIC_Hxy3z_Dyz_a+3*2*I_KINETIC_Fxyz_Dyz;
  abcd[728] = 4.0E0*I_KINETIC_Kx6z_Dyz_aa-2.0E0*4*I_KINETIC_Hx4z_Dyz_a-2.0E0*5*I_KINETIC_Hx4z_Dyz_a+4*3*I_KINETIC_Fx2z_Dyz;
  abcd[729] = 4.0E0*I_KINETIC_K5y2z_Dyz_aa-2.0E0*1*I_KINETIC_H5y_Dyz_a;
  abcd[730] = 4.0E0*I_KINETIC_K4y3z_Dyz_aa-2.0E0*1*I_KINETIC_H4yz_Dyz_a-2.0E0*2*I_KINETIC_H4yz_Dyz_a;
  abcd[731] = 4.0E0*I_KINETIC_K3y4z_Dyz_aa-2.0E0*2*I_KINETIC_H3y2z_Dyz_a-2.0E0*3*I_KINETIC_H3y2z_Dyz_a+2*1*I_KINETIC_F3y_Dyz;
  abcd[732] = 4.0E0*I_KINETIC_K2y5z_Dyz_aa-2.0E0*3*I_KINETIC_H2y3z_Dyz_a-2.0E0*4*I_KINETIC_H2y3z_Dyz_a+3*2*I_KINETIC_F2yz_Dyz;
  abcd[733] = 4.0E0*I_KINETIC_Ky6z_Dyz_aa-2.0E0*4*I_KINETIC_Hy4z_Dyz_a-2.0E0*5*I_KINETIC_Hy4z_Dyz_a+4*3*I_KINETIC_Fy2z_Dyz;
  abcd[734] = 4.0E0*I_KINETIC_K7z_Dyz_aa-2.0E0*5*I_KINETIC_H5z_Dyz_a-2.0E0*6*I_KINETIC_H5z_Dyz_a+5*4*I_KINETIC_F3z_Dyz;
  abcd[735] = 4.0E0*I_KINETIC_K5x2z_D2z_aa-2.0E0*1*I_KINETIC_H5x_D2z_a;
  abcd[736] = 4.0E0*I_KINETIC_K4xy2z_D2z_aa-2.0E0*1*I_KINETIC_H4xy_D2z_a;
  abcd[737] = 4.0E0*I_KINETIC_K4x3z_D2z_aa-2.0E0*1*I_KINETIC_H4xz_D2z_a-2.0E0*2*I_KINETIC_H4xz_D2z_a;
  abcd[738] = 4.0E0*I_KINETIC_K3x2y2z_D2z_aa-2.0E0*1*I_KINETIC_H3x2y_D2z_a;
  abcd[739] = 4.0E0*I_KINETIC_K3xy3z_D2z_aa-2.0E0*1*I_KINETIC_H3xyz_D2z_a-2.0E0*2*I_KINETIC_H3xyz_D2z_a;
  abcd[740] = 4.0E0*I_KINETIC_K3x4z_D2z_aa-2.0E0*2*I_KINETIC_H3x2z_D2z_a-2.0E0*3*I_KINETIC_H3x2z_D2z_a+2*1*I_KINETIC_F3x_D2z;
  abcd[741] = 4.0E0*I_KINETIC_K2x3y2z_D2z_aa-2.0E0*1*I_KINETIC_H2x3y_D2z_a;
  abcd[742] = 4.0E0*I_KINETIC_K2x2y3z_D2z_aa-2.0E0*1*I_KINETIC_H2x2yz_D2z_a-2.0E0*2*I_KINETIC_H2x2yz_D2z_a;
  abcd[743] = 4.0E0*I_KINETIC_K2xy4z_D2z_aa-2.0E0*2*I_KINETIC_H2xy2z_D2z_a-2.0E0*3*I_KINETIC_H2xy2z_D2z_a+2*1*I_KINETIC_F2xy_D2z;
  abcd[744] = 4.0E0*I_KINETIC_K2x5z_D2z_aa-2.0E0*3*I_KINETIC_H2x3z_D2z_a-2.0E0*4*I_KINETIC_H2x3z_D2z_a+3*2*I_KINETIC_F2xz_D2z;
  abcd[745] = 4.0E0*I_KINETIC_Kx4y2z_D2z_aa-2.0E0*1*I_KINETIC_Hx4y_D2z_a;
  abcd[746] = 4.0E0*I_KINETIC_Kx3y3z_D2z_aa-2.0E0*1*I_KINETIC_Hx3yz_D2z_a-2.0E0*2*I_KINETIC_Hx3yz_D2z_a;
  abcd[747] = 4.0E0*I_KINETIC_Kx2y4z_D2z_aa-2.0E0*2*I_KINETIC_Hx2y2z_D2z_a-2.0E0*3*I_KINETIC_Hx2y2z_D2z_a+2*1*I_KINETIC_Fx2y_D2z;
  abcd[748] = 4.0E0*I_KINETIC_Kxy5z_D2z_aa-2.0E0*3*I_KINETIC_Hxy3z_D2z_a-2.0E0*4*I_KINETIC_Hxy3z_D2z_a+3*2*I_KINETIC_Fxyz_D2z;
  abcd[749] = 4.0E0*I_KINETIC_Kx6z_D2z_aa-2.0E0*4*I_KINETIC_Hx4z_D2z_a-2.0E0*5*I_KINETIC_Hx4z_D2z_a+4*3*I_KINETIC_Fx2z_D2z;
  abcd[750] = 4.0E0*I_KINETIC_K5y2z_D2z_aa-2.0E0*1*I_KINETIC_H5y_D2z_a;
  abcd[751] = 4.0E0*I_KINETIC_K4y3z_D2z_aa-2.0E0*1*I_KINETIC_H4yz_D2z_a-2.0E0*2*I_KINETIC_H4yz_D2z_a;
  abcd[752] = 4.0E0*I_KINETIC_K3y4z_D2z_aa-2.0E0*2*I_KINETIC_H3y2z_D2z_a-2.0E0*3*I_KINETIC_H3y2z_D2z_a+2*1*I_KINETIC_F3y_D2z;
  abcd[753] = 4.0E0*I_KINETIC_K2y5z_D2z_aa-2.0E0*3*I_KINETIC_H2y3z_D2z_a-2.0E0*4*I_KINETIC_H2y3z_D2z_a+3*2*I_KINETIC_F2yz_D2z;
  abcd[754] = 4.0E0*I_KINETIC_Ky6z_D2z_aa-2.0E0*4*I_KINETIC_Hy4z_D2z_a-2.0E0*5*I_KINETIC_Hy4z_D2z_a+4*3*I_KINETIC_Fy2z_D2z;
  abcd[755] = 4.0E0*I_KINETIC_K7z_D2z_aa-2.0E0*5*I_KINETIC_H5z_D2z_a-2.0E0*6*I_KINETIC_H5z_D2z_a+5*4*I_KINETIC_F3z_D2z;
}
