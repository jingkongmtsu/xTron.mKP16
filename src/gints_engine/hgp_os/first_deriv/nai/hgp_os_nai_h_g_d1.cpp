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

void hgp_os_nai_h_g_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_L8x_S = 0.0E0;
  Double I_NAI_L7xy_S = 0.0E0;
  Double I_NAI_L7xz_S = 0.0E0;
  Double I_NAI_L6x2y_S = 0.0E0;
  Double I_NAI_L6xyz_S = 0.0E0;
  Double I_NAI_L6x2z_S = 0.0E0;
  Double I_NAI_L5x3y_S = 0.0E0;
  Double I_NAI_L5x2yz_S = 0.0E0;
  Double I_NAI_L5xy2z_S = 0.0E0;
  Double I_NAI_L5x3z_S = 0.0E0;
  Double I_NAI_L4x4y_S = 0.0E0;
  Double I_NAI_L4x3yz_S = 0.0E0;
  Double I_NAI_L4x2y2z_S = 0.0E0;
  Double I_NAI_L4xy3z_S = 0.0E0;
  Double I_NAI_L4x4z_S = 0.0E0;
  Double I_NAI_L3x5y_S = 0.0E0;
  Double I_NAI_L3x4yz_S = 0.0E0;
  Double I_NAI_L3x3y2z_S = 0.0E0;
  Double I_NAI_L3x2y3z_S = 0.0E0;
  Double I_NAI_L3xy4z_S = 0.0E0;
  Double I_NAI_L3x5z_S = 0.0E0;
  Double I_NAI_L2x6y_S = 0.0E0;
  Double I_NAI_L2x5yz_S = 0.0E0;
  Double I_NAI_L2x4y2z_S = 0.0E0;
  Double I_NAI_L2x3y3z_S = 0.0E0;
  Double I_NAI_L2x2y4z_S = 0.0E0;
  Double I_NAI_L2xy5z_S = 0.0E0;
  Double I_NAI_L2x6z_S = 0.0E0;
  Double I_NAI_Lx7y_S = 0.0E0;
  Double I_NAI_Lx6yz_S = 0.0E0;
  Double I_NAI_Lx5y2z_S = 0.0E0;
  Double I_NAI_Lx4y3z_S = 0.0E0;
  Double I_NAI_Lx3y4z_S = 0.0E0;
  Double I_NAI_Lx2y5z_S = 0.0E0;
  Double I_NAI_Lxy6z_S = 0.0E0;
  Double I_NAI_Lx7z_S = 0.0E0;
  Double I_NAI_L8y_S = 0.0E0;
  Double I_NAI_L7yz_S = 0.0E0;
  Double I_NAI_L6y2z_S = 0.0E0;
  Double I_NAI_L5y3z_S = 0.0E0;
  Double I_NAI_L4y4z_S = 0.0E0;
  Double I_NAI_L3y5z_S = 0.0E0;
  Double I_NAI_L2y6z_S = 0.0E0;
  Double I_NAI_Ly7z_S = 0.0E0;
  Double I_NAI_L8z_S = 0.0E0;
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
  Double I_NAI_N10x_S_a = 0.0E0;
  Double I_NAI_N9xy_S_a = 0.0E0;
  Double I_NAI_N9xz_S_a = 0.0E0;
  Double I_NAI_N8x2y_S_a = 0.0E0;
  Double I_NAI_N8xyz_S_a = 0.0E0;
  Double I_NAI_N8x2z_S_a = 0.0E0;
  Double I_NAI_N7x3y_S_a = 0.0E0;
  Double I_NAI_N7x2yz_S_a = 0.0E0;
  Double I_NAI_N7xy2z_S_a = 0.0E0;
  Double I_NAI_N7x3z_S_a = 0.0E0;
  Double I_NAI_N6x4y_S_a = 0.0E0;
  Double I_NAI_N6x3yz_S_a = 0.0E0;
  Double I_NAI_N6x2y2z_S_a = 0.0E0;
  Double I_NAI_N6xy3z_S_a = 0.0E0;
  Double I_NAI_N6x4z_S_a = 0.0E0;
  Double I_NAI_N5x5y_S_a = 0.0E0;
  Double I_NAI_N5x4yz_S_a = 0.0E0;
  Double I_NAI_N5x3y2z_S_a = 0.0E0;
  Double I_NAI_N5x2y3z_S_a = 0.0E0;
  Double I_NAI_N5xy4z_S_a = 0.0E0;
  Double I_NAI_N5x5z_S_a = 0.0E0;
  Double I_NAI_N4x6y_S_a = 0.0E0;
  Double I_NAI_N4x5yz_S_a = 0.0E0;
  Double I_NAI_N4x4y2z_S_a = 0.0E0;
  Double I_NAI_N4x3y3z_S_a = 0.0E0;
  Double I_NAI_N4x2y4z_S_a = 0.0E0;
  Double I_NAI_N4xy5z_S_a = 0.0E0;
  Double I_NAI_N4x6z_S_a = 0.0E0;
  Double I_NAI_N3x7y_S_a = 0.0E0;
  Double I_NAI_N3x6yz_S_a = 0.0E0;
  Double I_NAI_N3x5y2z_S_a = 0.0E0;
  Double I_NAI_N3x4y3z_S_a = 0.0E0;
  Double I_NAI_N3x3y4z_S_a = 0.0E0;
  Double I_NAI_N3x2y5z_S_a = 0.0E0;
  Double I_NAI_N3xy6z_S_a = 0.0E0;
  Double I_NAI_N3x7z_S_a = 0.0E0;
  Double I_NAI_N2x8y_S_a = 0.0E0;
  Double I_NAI_N2x7yz_S_a = 0.0E0;
  Double I_NAI_N2x6y2z_S_a = 0.0E0;
  Double I_NAI_N2x5y3z_S_a = 0.0E0;
  Double I_NAI_N2x4y4z_S_a = 0.0E0;
  Double I_NAI_N2x3y5z_S_a = 0.0E0;
  Double I_NAI_N2x2y6z_S_a = 0.0E0;
  Double I_NAI_N2xy7z_S_a = 0.0E0;
  Double I_NAI_N2x8z_S_a = 0.0E0;
  Double I_NAI_Nx9y_S_a = 0.0E0;
  Double I_NAI_Nx8yz_S_a = 0.0E0;
  Double I_NAI_Nx7y2z_S_a = 0.0E0;
  Double I_NAI_Nx6y3z_S_a = 0.0E0;
  Double I_NAI_Nx5y4z_S_a = 0.0E0;
  Double I_NAI_Nx4y5z_S_a = 0.0E0;
  Double I_NAI_Nx3y6z_S_a = 0.0E0;
  Double I_NAI_Nx2y7z_S_a = 0.0E0;
  Double I_NAI_Nxy8z_S_a = 0.0E0;
  Double I_NAI_Nx9z_S_a = 0.0E0;
  Double I_NAI_N10y_S_a = 0.0E0;
  Double I_NAI_N9yz_S_a = 0.0E0;
  Double I_NAI_N8y2z_S_a = 0.0E0;
  Double I_NAI_N7y3z_S_a = 0.0E0;
  Double I_NAI_N6y4z_S_a = 0.0E0;
  Double I_NAI_N5y5z_S_a = 0.0E0;
  Double I_NAI_N4y6z_S_a = 0.0E0;
  Double I_NAI_N3y7z_S_a = 0.0E0;
  Double I_NAI_N2y8z_S_a = 0.0E0;
  Double I_NAI_Ny9z_S_a = 0.0E0;
  Double I_NAI_N10z_S_a = 0.0E0;
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
  Double I_NAI_N10x_S_b = 0.0E0;
  Double I_NAI_N9xy_S_b = 0.0E0;
  Double I_NAI_N9xz_S_b = 0.0E0;
  Double I_NAI_N8x2y_S_b = 0.0E0;
  Double I_NAI_N8xyz_S_b = 0.0E0;
  Double I_NAI_N8x2z_S_b = 0.0E0;
  Double I_NAI_N7x3y_S_b = 0.0E0;
  Double I_NAI_N7x2yz_S_b = 0.0E0;
  Double I_NAI_N7xy2z_S_b = 0.0E0;
  Double I_NAI_N7x3z_S_b = 0.0E0;
  Double I_NAI_N6x4y_S_b = 0.0E0;
  Double I_NAI_N6x3yz_S_b = 0.0E0;
  Double I_NAI_N6x2y2z_S_b = 0.0E0;
  Double I_NAI_N6xy3z_S_b = 0.0E0;
  Double I_NAI_N6x4z_S_b = 0.0E0;
  Double I_NAI_N5x5y_S_b = 0.0E0;
  Double I_NAI_N5x4yz_S_b = 0.0E0;
  Double I_NAI_N5x3y2z_S_b = 0.0E0;
  Double I_NAI_N5x2y3z_S_b = 0.0E0;
  Double I_NAI_N5xy4z_S_b = 0.0E0;
  Double I_NAI_N5x5z_S_b = 0.0E0;
  Double I_NAI_N4x6y_S_b = 0.0E0;
  Double I_NAI_N4x5yz_S_b = 0.0E0;
  Double I_NAI_N4x4y2z_S_b = 0.0E0;
  Double I_NAI_N4x3y3z_S_b = 0.0E0;
  Double I_NAI_N4x2y4z_S_b = 0.0E0;
  Double I_NAI_N4xy5z_S_b = 0.0E0;
  Double I_NAI_N4x6z_S_b = 0.0E0;
  Double I_NAI_N3x7y_S_b = 0.0E0;
  Double I_NAI_N3x6yz_S_b = 0.0E0;
  Double I_NAI_N3x5y2z_S_b = 0.0E0;
  Double I_NAI_N3x4y3z_S_b = 0.0E0;
  Double I_NAI_N3x3y4z_S_b = 0.0E0;
  Double I_NAI_N3x2y5z_S_b = 0.0E0;
  Double I_NAI_N3xy6z_S_b = 0.0E0;
  Double I_NAI_N3x7z_S_b = 0.0E0;
  Double I_NAI_N2x8y_S_b = 0.0E0;
  Double I_NAI_N2x7yz_S_b = 0.0E0;
  Double I_NAI_N2x6y2z_S_b = 0.0E0;
  Double I_NAI_N2x5y3z_S_b = 0.0E0;
  Double I_NAI_N2x4y4z_S_b = 0.0E0;
  Double I_NAI_N2x3y5z_S_b = 0.0E0;
  Double I_NAI_N2x2y6z_S_b = 0.0E0;
  Double I_NAI_N2xy7z_S_b = 0.0E0;
  Double I_NAI_N2x8z_S_b = 0.0E0;
  Double I_NAI_Nx9y_S_b = 0.0E0;
  Double I_NAI_Nx8yz_S_b = 0.0E0;
  Double I_NAI_Nx7y2z_S_b = 0.0E0;
  Double I_NAI_Nx6y3z_S_b = 0.0E0;
  Double I_NAI_Nx5y4z_S_b = 0.0E0;
  Double I_NAI_Nx4y5z_S_b = 0.0E0;
  Double I_NAI_Nx3y6z_S_b = 0.0E0;
  Double I_NAI_Nx2y7z_S_b = 0.0E0;
  Double I_NAI_Nxy8z_S_b = 0.0E0;
  Double I_NAI_Nx9z_S_b = 0.0E0;
  Double I_NAI_N10y_S_b = 0.0E0;
  Double I_NAI_N9yz_S_b = 0.0E0;
  Double I_NAI_N8y2z_S_b = 0.0E0;
  Double I_NAI_N7y3z_S_b = 0.0E0;
  Double I_NAI_N6y4z_S_b = 0.0E0;
  Double I_NAI_N5y5z_S_b = 0.0E0;
  Double I_NAI_N4y6z_S_b = 0.0E0;
  Double I_NAI_N3y7z_S_b = 0.0E0;
  Double I_NAI_N2y8z_S_b = 0.0E0;
  Double I_NAI_Ny9z_S_b = 0.0E0;
  Double I_NAI_N10z_S_b = 0.0E0;
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
      Double I_NAI_S_S_M10_vrr  = 0.0E0;

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
      double I_NAI_S_S_M10_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER55;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER53*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER51*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER49*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER47*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M10_vrr;
        I_NAI_S_S_M10_vrr = ONEOVER21*I_NAI_S_S_M10_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M10_vrr  = f*I_NAI_S_S_M10_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_NAI_S_S_M9_vrr  = ONEOVER19*(u2*I_NAI_S_S_M10_vrr+f);
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
        I_NAI_S_S_M10_vrr_d = oneO2u*(19.0E0*I_NAI_S_S_M9_vrr_d-f);

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
        I_NAI_S_S_M10_vrr = static_cast<Double>(I_NAI_S_S_M10_vrr_d);

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
        I_NAI_S_S_M10_vrr = oneO2u*(19.0E0*I_NAI_S_S_M9_vrr-f);

#endif

      }


      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M9
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M9
       * RHS shell quartet name: SQ_NAI_S_S_M10
       ************************************************************/
      Double I_NAI_Px_S_M9_vrr = PAX*I_NAI_S_S_M9_vrr-PNX*I_NAI_S_S_M10_vrr;
      Double I_NAI_Py_S_M9_vrr = PAY*I_NAI_S_S_M9_vrr-PNY*I_NAI_S_S_M10_vrr;
      Double I_NAI_Pz_S_M9_vrr = PAZ*I_NAI_S_S_M9_vrr-PNZ*I_NAI_S_S_M10_vrr;

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
       * shell quartet name: SQ_NAI_D_S_M8
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M8
       * RHS shell quartet name: SQ_NAI_P_S_M9
       * RHS shell quartet name: SQ_NAI_S_S_M8
       * RHS shell quartet name: SQ_NAI_S_S_M9
       ************************************************************/
      Double I_NAI_D2x_S_M8_vrr = PAX*I_NAI_Px_S_M8_vrr-PNX*I_NAI_Px_S_M9_vrr+oned2z*I_NAI_S_S_M8_vrr-oned2z*I_NAI_S_S_M9_vrr;
      Double I_NAI_D2y_S_M8_vrr = PAY*I_NAI_Py_S_M8_vrr-PNY*I_NAI_Py_S_M9_vrr+oned2z*I_NAI_S_S_M8_vrr-oned2z*I_NAI_S_S_M9_vrr;
      Double I_NAI_D2z_S_M8_vrr = PAZ*I_NAI_Pz_S_M8_vrr-PNZ*I_NAI_Pz_S_M9_vrr+oned2z*I_NAI_S_S_M8_vrr-oned2z*I_NAI_S_S_M9_vrr;

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
       * shell quartet name: SQ_NAI_F_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M7
       * RHS shell quartet name: SQ_NAI_D_S_M8
       * RHS shell quartet name: SQ_NAI_P_S_M7
       * RHS shell quartet name: SQ_NAI_P_S_M8
       ************************************************************/
      Double I_NAI_F3x_S_M7_vrr = PAX*I_NAI_D2x_S_M7_vrr-PNX*I_NAI_D2x_S_M8_vrr+2*oned2z*I_NAI_Px_S_M7_vrr-2*oned2z*I_NAI_Px_S_M8_vrr;
      Double I_NAI_F3y_S_M7_vrr = PAY*I_NAI_D2y_S_M7_vrr-PNY*I_NAI_D2y_S_M8_vrr+2*oned2z*I_NAI_Py_S_M7_vrr-2*oned2z*I_NAI_Py_S_M8_vrr;
      Double I_NAI_F3z_S_M7_vrr = PAZ*I_NAI_D2z_S_M7_vrr-PNZ*I_NAI_D2z_S_M8_vrr+2*oned2z*I_NAI_Pz_S_M7_vrr-2*oned2z*I_NAI_Pz_S_M8_vrr;

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
       * shell quartet name: SQ_NAI_G_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 12 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M6
       * RHS shell quartet name: SQ_NAI_F_S_M7
       * RHS shell quartet name: SQ_NAI_D_S_M6
       * RHS shell quartet name: SQ_NAI_D_S_M7
       ************************************************************/
      Double I_NAI_G4x_S_M6_vrr = PAX*I_NAI_F3x_S_M6_vrr-PNX*I_NAI_F3x_S_M7_vrr+3*oned2z*I_NAI_D2x_S_M6_vrr-3*oned2z*I_NAI_D2x_S_M7_vrr;
      Double I_NAI_G4y_S_M6_vrr = PAY*I_NAI_F3y_S_M6_vrr-PNY*I_NAI_F3y_S_M7_vrr+3*oned2z*I_NAI_D2y_S_M6_vrr-3*oned2z*I_NAI_D2y_S_M7_vrr;
      Double I_NAI_G4z_S_M6_vrr = PAZ*I_NAI_F3z_S_M6_vrr-PNZ*I_NAI_F3z_S_M7_vrr+3*oned2z*I_NAI_D2z_S_M6_vrr-3*oned2z*I_NAI_D2z_S_M7_vrr;

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
       * shell quartet name: SQ_NAI_H_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 13 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M5
       * RHS shell quartet name: SQ_NAI_G_S_M6
       * RHS shell quartet name: SQ_NAI_F_S_M5
       * RHS shell quartet name: SQ_NAI_F_S_M6
       ************************************************************/
      Double I_NAI_H5x_S_M5_vrr = PAX*I_NAI_G4x_S_M5_vrr-PNX*I_NAI_G4x_S_M6_vrr+4*oned2z*I_NAI_F3x_S_M5_vrr-4*oned2z*I_NAI_F3x_S_M6_vrr;
      Double I_NAI_H4xy_S_M5_vrr = PAY*I_NAI_G4x_S_M5_vrr-PNY*I_NAI_G4x_S_M6_vrr;
      Double I_NAI_H4xz_S_M5_vrr = PAZ*I_NAI_G4x_S_M5_vrr-PNZ*I_NAI_G4x_S_M6_vrr;
      Double I_NAI_Hx4y_S_M5_vrr = PAX*I_NAI_G4y_S_M5_vrr-PNX*I_NAI_G4y_S_M6_vrr;
      Double I_NAI_Hx4z_S_M5_vrr = PAX*I_NAI_G4z_S_M5_vrr-PNX*I_NAI_G4z_S_M6_vrr;
      Double I_NAI_H5y_S_M5_vrr = PAY*I_NAI_G4y_S_M5_vrr-PNY*I_NAI_G4y_S_M6_vrr+4*oned2z*I_NAI_F3y_S_M5_vrr-4*oned2z*I_NAI_F3y_S_M6_vrr;
      Double I_NAI_H4yz_S_M5_vrr = PAZ*I_NAI_G4y_S_M5_vrr-PNZ*I_NAI_G4y_S_M6_vrr;
      Double I_NAI_H5z_S_M5_vrr = PAZ*I_NAI_G4z_S_M5_vrr-PNZ*I_NAI_G4z_S_M6_vrr+4*oned2z*I_NAI_F3z_S_M5_vrr-4*oned2z*I_NAI_F3z_S_M6_vrr;

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
       * shell quartet name: SQ_NAI_I_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 14 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S_M4
       * RHS shell quartet name: SQ_NAI_H_S_M5
       * RHS shell quartet name: SQ_NAI_G_S_M4
       * RHS shell quartet name: SQ_NAI_G_S_M5
       ************************************************************/
      Double I_NAI_I6x_S_M4_vrr = PAX*I_NAI_H5x_S_M4_vrr-PNX*I_NAI_H5x_S_M5_vrr+5*oned2z*I_NAI_G4x_S_M4_vrr-5*oned2z*I_NAI_G4x_S_M5_vrr;
      Double I_NAI_I5xy_S_M4_vrr = PAY*I_NAI_H5x_S_M4_vrr-PNY*I_NAI_H5x_S_M5_vrr;
      Double I_NAI_I5xz_S_M4_vrr = PAZ*I_NAI_H5x_S_M4_vrr-PNZ*I_NAI_H5x_S_M5_vrr;
      Double I_NAI_I4x2y_S_M4_vrr = PAY*I_NAI_H4xy_S_M4_vrr-PNY*I_NAI_H4xy_S_M5_vrr+oned2z*I_NAI_G4x_S_M4_vrr-oned2z*I_NAI_G4x_S_M5_vrr;
      Double I_NAI_I4x2z_S_M4_vrr = PAZ*I_NAI_H4xz_S_M4_vrr-PNZ*I_NAI_H4xz_S_M5_vrr+oned2z*I_NAI_G4x_S_M4_vrr-oned2z*I_NAI_G4x_S_M5_vrr;
      Double I_NAI_I2x4y_S_M4_vrr = PAX*I_NAI_Hx4y_S_M4_vrr-PNX*I_NAI_Hx4y_S_M5_vrr+oned2z*I_NAI_G4y_S_M4_vrr-oned2z*I_NAI_G4y_S_M5_vrr;
      Double I_NAI_I2x4z_S_M4_vrr = PAX*I_NAI_Hx4z_S_M4_vrr-PNX*I_NAI_Hx4z_S_M5_vrr+oned2z*I_NAI_G4z_S_M4_vrr-oned2z*I_NAI_G4z_S_M5_vrr;
      Double I_NAI_Ix5y_S_M4_vrr = PAX*I_NAI_H5y_S_M4_vrr-PNX*I_NAI_H5y_S_M5_vrr;
      Double I_NAI_Ix5z_S_M4_vrr = PAX*I_NAI_H5z_S_M4_vrr-PNX*I_NAI_H5z_S_M5_vrr;
      Double I_NAI_I6y_S_M4_vrr = PAY*I_NAI_H5y_S_M4_vrr-PNY*I_NAI_H5y_S_M5_vrr+5*oned2z*I_NAI_G4y_S_M4_vrr-5*oned2z*I_NAI_G4y_S_M5_vrr;
      Double I_NAI_I5yz_S_M4_vrr = PAZ*I_NAI_H5y_S_M4_vrr-PNZ*I_NAI_H5y_S_M5_vrr;
      Double I_NAI_I4y2z_S_M4_vrr = PAZ*I_NAI_H4yz_S_M4_vrr-PNZ*I_NAI_H4yz_S_M5_vrr+oned2z*I_NAI_G4y_S_M4_vrr-oned2z*I_NAI_G4y_S_M5_vrr;
      Double I_NAI_Iy5z_S_M4_vrr = PAY*I_NAI_H5z_S_M4_vrr-PNY*I_NAI_H5z_S_M5_vrr;
      Double I_NAI_I6z_S_M4_vrr = PAZ*I_NAI_H5z_S_M4_vrr-PNZ*I_NAI_H5z_S_M5_vrr+5*oned2z*I_NAI_G4z_S_M4_vrr-5*oned2z*I_NAI_G4z_S_M5_vrr;

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
       * shell quartet name: SQ_NAI_K_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 16 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_I_S_M3
       * RHS shell quartet name: SQ_NAI_I_S_M4
       * RHS shell quartet name: SQ_NAI_H_S_M3
       * RHS shell quartet name: SQ_NAI_H_S_M4
       ************************************************************/
      Double I_NAI_K7x_S_M3_vrr = PAX*I_NAI_I6x_S_M3_vrr-PNX*I_NAI_I6x_S_M4_vrr+6*oned2z*I_NAI_H5x_S_M3_vrr-6*oned2z*I_NAI_H5x_S_M4_vrr;
      Double I_NAI_K6xy_S_M3_vrr = PAY*I_NAI_I6x_S_M3_vrr-PNY*I_NAI_I6x_S_M4_vrr;
      Double I_NAI_K6xz_S_M3_vrr = PAZ*I_NAI_I6x_S_M3_vrr-PNZ*I_NAI_I6x_S_M4_vrr;
      Double I_NAI_K5x2y_S_M3_vrr = PAY*I_NAI_I5xy_S_M3_vrr-PNY*I_NAI_I5xy_S_M4_vrr+oned2z*I_NAI_H5x_S_M3_vrr-oned2z*I_NAI_H5x_S_M4_vrr;
      Double I_NAI_K5x2z_S_M3_vrr = PAZ*I_NAI_I5xz_S_M3_vrr-PNZ*I_NAI_I5xz_S_M4_vrr+oned2z*I_NAI_H5x_S_M3_vrr-oned2z*I_NAI_H5x_S_M4_vrr;
      Double I_NAI_K4x3y_S_M3_vrr = PAY*I_NAI_I4x2y_S_M3_vrr-PNY*I_NAI_I4x2y_S_M4_vrr+2*oned2z*I_NAI_H4xy_S_M3_vrr-2*oned2z*I_NAI_H4xy_S_M4_vrr;
      Double I_NAI_K4x3z_S_M3_vrr = PAZ*I_NAI_I4x2z_S_M3_vrr-PNZ*I_NAI_I4x2z_S_M4_vrr+2*oned2z*I_NAI_H4xz_S_M3_vrr-2*oned2z*I_NAI_H4xz_S_M4_vrr;
      Double I_NAI_K3x4y_S_M3_vrr = PAX*I_NAI_I2x4y_S_M3_vrr-PNX*I_NAI_I2x4y_S_M4_vrr+2*oned2z*I_NAI_Hx4y_S_M3_vrr-2*oned2z*I_NAI_Hx4y_S_M4_vrr;
      Double I_NAI_K3x4z_S_M3_vrr = PAX*I_NAI_I2x4z_S_M3_vrr-PNX*I_NAI_I2x4z_S_M4_vrr+2*oned2z*I_NAI_Hx4z_S_M3_vrr-2*oned2z*I_NAI_Hx4z_S_M4_vrr;
      Double I_NAI_K2x5y_S_M3_vrr = PAX*I_NAI_Ix5y_S_M3_vrr-PNX*I_NAI_Ix5y_S_M4_vrr+oned2z*I_NAI_H5y_S_M3_vrr-oned2z*I_NAI_H5y_S_M4_vrr;
      Double I_NAI_K2x5z_S_M3_vrr = PAX*I_NAI_Ix5z_S_M3_vrr-PNX*I_NAI_Ix5z_S_M4_vrr+oned2z*I_NAI_H5z_S_M3_vrr-oned2z*I_NAI_H5z_S_M4_vrr;
      Double I_NAI_Kx6y_S_M3_vrr = PAX*I_NAI_I6y_S_M3_vrr-PNX*I_NAI_I6y_S_M4_vrr;
      Double I_NAI_Kx6z_S_M3_vrr = PAX*I_NAI_I6z_S_M3_vrr-PNX*I_NAI_I6z_S_M4_vrr;
      Double I_NAI_K7y_S_M3_vrr = PAY*I_NAI_I6y_S_M3_vrr-PNY*I_NAI_I6y_S_M4_vrr+6*oned2z*I_NAI_H5y_S_M3_vrr-6*oned2z*I_NAI_H5y_S_M4_vrr;
      Double I_NAI_K6yz_S_M3_vrr = PAZ*I_NAI_I6y_S_M3_vrr-PNZ*I_NAI_I6y_S_M4_vrr;
      Double I_NAI_K5y2z_S_M3_vrr = PAZ*I_NAI_I5yz_S_M3_vrr-PNZ*I_NAI_I5yz_S_M4_vrr+oned2z*I_NAI_H5y_S_M3_vrr-oned2z*I_NAI_H5y_S_M4_vrr;
      Double I_NAI_K4y3z_S_M3_vrr = PAZ*I_NAI_I4y2z_S_M3_vrr-PNZ*I_NAI_I4y2z_S_M4_vrr+2*oned2z*I_NAI_H4yz_S_M3_vrr-2*oned2z*I_NAI_H4yz_S_M4_vrr;
      Double I_NAI_K2y5z_S_M3_vrr = PAY*I_NAI_Iy5z_S_M3_vrr-PNY*I_NAI_Iy5z_S_M4_vrr+oned2z*I_NAI_H5z_S_M3_vrr-oned2z*I_NAI_H5z_S_M4_vrr;
      Double I_NAI_Ky6z_S_M3_vrr = PAY*I_NAI_I6z_S_M3_vrr-PNY*I_NAI_I6z_S_M4_vrr;
      Double I_NAI_K7z_S_M3_vrr = PAZ*I_NAI_I6z_S_M3_vrr-PNZ*I_NAI_I6z_S_M4_vrr+6*oned2z*I_NAI_H5z_S_M3_vrr-6*oned2z*I_NAI_H5z_S_M4_vrr;

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
       * shell quartet name: SQ_NAI_L_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 18 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_K_S_M2
       * RHS shell quartet name: SQ_NAI_K_S_M3
       * RHS shell quartet name: SQ_NAI_I_S_M2
       * RHS shell quartet name: SQ_NAI_I_S_M3
       ************************************************************/
      Double I_NAI_L8x_S_M2_vrr = PAX*I_NAI_K7x_S_M2_vrr-PNX*I_NAI_K7x_S_M3_vrr+7*oned2z*I_NAI_I6x_S_M2_vrr-7*oned2z*I_NAI_I6x_S_M3_vrr;
      Double I_NAI_L7xy_S_M2_vrr = PAY*I_NAI_K7x_S_M2_vrr-PNY*I_NAI_K7x_S_M3_vrr;
      Double I_NAI_L7xz_S_M2_vrr = PAZ*I_NAI_K7x_S_M2_vrr-PNZ*I_NAI_K7x_S_M3_vrr;
      Double I_NAI_L6x2y_S_M2_vrr = PAY*I_NAI_K6xy_S_M2_vrr-PNY*I_NAI_K6xy_S_M3_vrr+oned2z*I_NAI_I6x_S_M2_vrr-oned2z*I_NAI_I6x_S_M3_vrr;
      Double I_NAI_L6x2z_S_M2_vrr = PAZ*I_NAI_K6xz_S_M2_vrr-PNZ*I_NAI_K6xz_S_M3_vrr+oned2z*I_NAI_I6x_S_M2_vrr-oned2z*I_NAI_I6x_S_M3_vrr;
      Double I_NAI_L5x3y_S_M2_vrr = PAY*I_NAI_K5x2y_S_M2_vrr-PNY*I_NAI_K5x2y_S_M3_vrr+2*oned2z*I_NAI_I5xy_S_M2_vrr-2*oned2z*I_NAI_I5xy_S_M3_vrr;
      Double I_NAI_L5x3z_S_M2_vrr = PAZ*I_NAI_K5x2z_S_M2_vrr-PNZ*I_NAI_K5x2z_S_M3_vrr+2*oned2z*I_NAI_I5xz_S_M2_vrr-2*oned2z*I_NAI_I5xz_S_M3_vrr;
      Double I_NAI_L4x4y_S_M2_vrr = PAY*I_NAI_K4x3y_S_M2_vrr-PNY*I_NAI_K4x3y_S_M3_vrr+3*oned2z*I_NAI_I4x2y_S_M2_vrr-3*oned2z*I_NAI_I4x2y_S_M3_vrr;
      Double I_NAI_L4x3yz_S_M2_vrr = PAZ*I_NAI_K4x3y_S_M2_vrr-PNZ*I_NAI_K4x3y_S_M3_vrr;
      Double I_NAI_L4x4z_S_M2_vrr = PAZ*I_NAI_K4x3z_S_M2_vrr-PNZ*I_NAI_K4x3z_S_M3_vrr+3*oned2z*I_NAI_I4x2z_S_M2_vrr-3*oned2z*I_NAI_I4x2z_S_M3_vrr;
      Double I_NAI_L3x5y_S_M2_vrr = PAX*I_NAI_K2x5y_S_M2_vrr-PNX*I_NAI_K2x5y_S_M3_vrr+2*oned2z*I_NAI_Ix5y_S_M2_vrr-2*oned2z*I_NAI_Ix5y_S_M3_vrr;
      Double I_NAI_L3x4yz_S_M2_vrr = PAZ*I_NAI_K3x4y_S_M2_vrr-PNZ*I_NAI_K3x4y_S_M3_vrr;
      Double I_NAI_L3xy4z_S_M2_vrr = PAY*I_NAI_K3x4z_S_M2_vrr-PNY*I_NAI_K3x4z_S_M3_vrr;
      Double I_NAI_L3x5z_S_M2_vrr = PAX*I_NAI_K2x5z_S_M2_vrr-PNX*I_NAI_K2x5z_S_M3_vrr+2*oned2z*I_NAI_Ix5z_S_M2_vrr-2*oned2z*I_NAI_Ix5z_S_M3_vrr;
      Double I_NAI_L2x6y_S_M2_vrr = PAX*I_NAI_Kx6y_S_M2_vrr-PNX*I_NAI_Kx6y_S_M3_vrr+oned2z*I_NAI_I6y_S_M2_vrr-oned2z*I_NAI_I6y_S_M3_vrr;
      Double I_NAI_L2x6z_S_M2_vrr = PAX*I_NAI_Kx6z_S_M2_vrr-PNX*I_NAI_Kx6z_S_M3_vrr+oned2z*I_NAI_I6z_S_M2_vrr-oned2z*I_NAI_I6z_S_M3_vrr;
      Double I_NAI_Lx7y_S_M2_vrr = PAX*I_NAI_K7y_S_M2_vrr-PNX*I_NAI_K7y_S_M3_vrr;
      Double I_NAI_Lx7z_S_M2_vrr = PAX*I_NAI_K7z_S_M2_vrr-PNX*I_NAI_K7z_S_M3_vrr;
      Double I_NAI_L8y_S_M2_vrr = PAY*I_NAI_K7y_S_M2_vrr-PNY*I_NAI_K7y_S_M3_vrr+7*oned2z*I_NAI_I6y_S_M2_vrr-7*oned2z*I_NAI_I6y_S_M3_vrr;
      Double I_NAI_L7yz_S_M2_vrr = PAZ*I_NAI_K7y_S_M2_vrr-PNZ*I_NAI_K7y_S_M3_vrr;
      Double I_NAI_L6y2z_S_M2_vrr = PAZ*I_NAI_K6yz_S_M2_vrr-PNZ*I_NAI_K6yz_S_M3_vrr+oned2z*I_NAI_I6y_S_M2_vrr-oned2z*I_NAI_I6y_S_M3_vrr;
      Double I_NAI_L5y3z_S_M2_vrr = PAZ*I_NAI_K5y2z_S_M2_vrr-PNZ*I_NAI_K5y2z_S_M3_vrr+2*oned2z*I_NAI_I5yz_S_M2_vrr-2*oned2z*I_NAI_I5yz_S_M3_vrr;
      Double I_NAI_L4y4z_S_M2_vrr = PAZ*I_NAI_K4y3z_S_M2_vrr-PNZ*I_NAI_K4y3z_S_M3_vrr+3*oned2z*I_NAI_I4y2z_S_M2_vrr-3*oned2z*I_NAI_I4y2z_S_M3_vrr;
      Double I_NAI_L3y5z_S_M2_vrr = PAY*I_NAI_K2y5z_S_M2_vrr-PNY*I_NAI_K2y5z_S_M3_vrr+2*oned2z*I_NAI_Iy5z_S_M2_vrr-2*oned2z*I_NAI_Iy5z_S_M3_vrr;
      Double I_NAI_L2y6z_S_M2_vrr = PAY*I_NAI_Ky6z_S_M2_vrr-PNY*I_NAI_Ky6z_S_M3_vrr+oned2z*I_NAI_I6z_S_M2_vrr-oned2z*I_NAI_I6z_S_M3_vrr;
      Double I_NAI_Ly7z_S_M2_vrr = PAY*I_NAI_K7z_S_M2_vrr-PNY*I_NAI_K7z_S_M3_vrr;
      Double I_NAI_L8z_S_M2_vrr = PAZ*I_NAI_K7z_S_M2_vrr-PNZ*I_NAI_K7z_S_M3_vrr+7*oned2z*I_NAI_I6z_S_M2_vrr-7*oned2z*I_NAI_I6z_S_M3_vrr;

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
       * shell quartet name: SQ_NAI_M_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 13 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_L_S_M1
       * RHS shell quartet name: SQ_NAI_L_S_M2
       * RHS shell quartet name: SQ_NAI_K_S_M1
       * RHS shell quartet name: SQ_NAI_K_S_M2
       ************************************************************/
      Double I_NAI_M9x_S_M1_vrr = PAX*I_NAI_L8x_S_M1_vrr-PNX*I_NAI_L8x_S_M2_vrr+8*oned2z*I_NAI_K7x_S_M1_vrr-8*oned2z*I_NAI_K7x_S_M2_vrr;
      Double I_NAI_M8xy_S_M1_vrr = PAY*I_NAI_L8x_S_M1_vrr-PNY*I_NAI_L8x_S_M2_vrr;
      Double I_NAI_M8xz_S_M1_vrr = PAZ*I_NAI_L8x_S_M1_vrr-PNZ*I_NAI_L8x_S_M2_vrr;
      Double I_NAI_M7x2y_S_M1_vrr = PAY*I_NAI_L7xy_S_M1_vrr-PNY*I_NAI_L7xy_S_M2_vrr+oned2z*I_NAI_K7x_S_M1_vrr-oned2z*I_NAI_K7x_S_M2_vrr;
      Double I_NAI_M7x2z_S_M1_vrr = PAZ*I_NAI_L7xz_S_M1_vrr-PNZ*I_NAI_L7xz_S_M2_vrr+oned2z*I_NAI_K7x_S_M1_vrr-oned2z*I_NAI_K7x_S_M2_vrr;
      Double I_NAI_M6x3y_S_M1_vrr = PAY*I_NAI_L6x2y_S_M1_vrr-PNY*I_NAI_L6x2y_S_M2_vrr+2*oned2z*I_NAI_K6xy_S_M1_vrr-2*oned2z*I_NAI_K6xy_S_M2_vrr;
      Double I_NAI_M6x2yz_S_M1_vrr = PAZ*I_NAI_L6x2y_S_M1_vrr-PNZ*I_NAI_L6x2y_S_M2_vrr;
      Double I_NAI_M6x3z_S_M1_vrr = PAZ*I_NAI_L6x2z_S_M1_vrr-PNZ*I_NAI_L6x2z_S_M2_vrr+2*oned2z*I_NAI_K6xz_S_M1_vrr-2*oned2z*I_NAI_K6xz_S_M2_vrr;
      Double I_NAI_M5x4y_S_M1_vrr = PAY*I_NAI_L5x3y_S_M1_vrr-PNY*I_NAI_L5x3y_S_M2_vrr+3*oned2z*I_NAI_K5x2y_S_M1_vrr-3*oned2z*I_NAI_K5x2y_S_M2_vrr;
      Double I_NAI_M5x3yz_S_M1_vrr = PAZ*I_NAI_L5x3y_S_M1_vrr-PNZ*I_NAI_L5x3y_S_M2_vrr;
      Double I_NAI_M5xy3z_S_M1_vrr = PAY*I_NAI_L5x3z_S_M1_vrr-PNY*I_NAI_L5x3z_S_M2_vrr;
      Double I_NAI_M5x4z_S_M1_vrr = PAZ*I_NAI_L5x3z_S_M1_vrr-PNZ*I_NAI_L5x3z_S_M2_vrr+3*oned2z*I_NAI_K5x2z_S_M1_vrr-3*oned2z*I_NAI_K5x2z_S_M2_vrr;
      Double I_NAI_M4x5y_S_M1_vrr = PAX*I_NAI_L3x5y_S_M1_vrr-PNX*I_NAI_L3x5y_S_M2_vrr+3*oned2z*I_NAI_K2x5y_S_M1_vrr-3*oned2z*I_NAI_K2x5y_S_M2_vrr;
      Double I_NAI_M4x4yz_S_M1_vrr = PAZ*I_NAI_L4x4y_S_M1_vrr-PNZ*I_NAI_L4x4y_S_M2_vrr;
      Double I_NAI_M4x3y2z_S_M1_vrr = PAZ*I_NAI_L4x3yz_S_M1_vrr-PNZ*I_NAI_L4x3yz_S_M2_vrr+oned2z*I_NAI_K4x3y_S_M1_vrr-oned2z*I_NAI_K4x3y_S_M2_vrr;
      Double I_NAI_M4xy4z_S_M1_vrr = PAY*I_NAI_L4x4z_S_M1_vrr-PNY*I_NAI_L4x4z_S_M2_vrr;
      Double I_NAI_M4x5z_S_M1_vrr = PAX*I_NAI_L3x5z_S_M1_vrr-PNX*I_NAI_L3x5z_S_M2_vrr+3*oned2z*I_NAI_K2x5z_S_M1_vrr-3*oned2z*I_NAI_K2x5z_S_M2_vrr;
      Double I_NAI_M3x6y_S_M1_vrr = PAX*I_NAI_L2x6y_S_M1_vrr-PNX*I_NAI_L2x6y_S_M2_vrr+2*oned2z*I_NAI_Kx6y_S_M1_vrr-2*oned2z*I_NAI_Kx6y_S_M2_vrr;
      Double I_NAI_M3x5yz_S_M1_vrr = PAZ*I_NAI_L3x5y_S_M1_vrr-PNZ*I_NAI_L3x5y_S_M2_vrr;
      Double I_NAI_M3x4y2z_S_M1_vrr = PAZ*I_NAI_L3x4yz_S_M1_vrr-PNZ*I_NAI_L3x4yz_S_M2_vrr+oned2z*I_NAI_K3x4y_S_M1_vrr-oned2z*I_NAI_K3x4y_S_M2_vrr;
      Double I_NAI_M3x2y4z_S_M1_vrr = PAY*I_NAI_L3xy4z_S_M1_vrr-PNY*I_NAI_L3xy4z_S_M2_vrr+oned2z*I_NAI_K3x4z_S_M1_vrr-oned2z*I_NAI_K3x4z_S_M2_vrr;
      Double I_NAI_M3xy5z_S_M1_vrr = PAY*I_NAI_L3x5z_S_M1_vrr-PNY*I_NAI_L3x5z_S_M2_vrr;
      Double I_NAI_M3x6z_S_M1_vrr = PAX*I_NAI_L2x6z_S_M1_vrr-PNX*I_NAI_L2x6z_S_M2_vrr+2*oned2z*I_NAI_Kx6z_S_M1_vrr-2*oned2z*I_NAI_Kx6z_S_M2_vrr;
      Double I_NAI_M2x7y_S_M1_vrr = PAX*I_NAI_Lx7y_S_M1_vrr-PNX*I_NAI_Lx7y_S_M2_vrr+oned2z*I_NAI_K7y_S_M1_vrr-oned2z*I_NAI_K7y_S_M2_vrr;
      Double I_NAI_M2x6yz_S_M1_vrr = PAZ*I_NAI_L2x6y_S_M1_vrr-PNZ*I_NAI_L2x6y_S_M2_vrr;
      Double I_NAI_M2xy6z_S_M1_vrr = PAY*I_NAI_L2x6z_S_M1_vrr-PNY*I_NAI_L2x6z_S_M2_vrr;
      Double I_NAI_M2x7z_S_M1_vrr = PAX*I_NAI_Lx7z_S_M1_vrr-PNX*I_NAI_Lx7z_S_M2_vrr+oned2z*I_NAI_K7z_S_M1_vrr-oned2z*I_NAI_K7z_S_M2_vrr;
      Double I_NAI_Mx8y_S_M1_vrr = PAX*I_NAI_L8y_S_M1_vrr-PNX*I_NAI_L8y_S_M2_vrr;
      Double I_NAI_Mx5y3z_S_M1_vrr = PAX*I_NAI_L5y3z_S_M1_vrr-PNX*I_NAI_L5y3z_S_M2_vrr;
      Double I_NAI_Mx4y4z_S_M1_vrr = PAX*I_NAI_L4y4z_S_M1_vrr-PNX*I_NAI_L4y4z_S_M2_vrr;
      Double I_NAI_Mx3y5z_S_M1_vrr = PAX*I_NAI_L3y5z_S_M1_vrr-PNX*I_NAI_L3y5z_S_M2_vrr;
      Double I_NAI_Mx8z_S_M1_vrr = PAX*I_NAI_L8z_S_M1_vrr-PNX*I_NAI_L8z_S_M2_vrr;
      Double I_NAI_M9y_S_M1_vrr = PAY*I_NAI_L8y_S_M1_vrr-PNY*I_NAI_L8y_S_M2_vrr+8*oned2z*I_NAI_K7y_S_M1_vrr-8*oned2z*I_NAI_K7y_S_M2_vrr;
      Double I_NAI_M8yz_S_M1_vrr = PAZ*I_NAI_L8y_S_M1_vrr-PNZ*I_NAI_L8y_S_M2_vrr;
      Double I_NAI_M7y2z_S_M1_vrr = PAZ*I_NAI_L7yz_S_M1_vrr-PNZ*I_NAI_L7yz_S_M2_vrr+oned2z*I_NAI_K7y_S_M1_vrr-oned2z*I_NAI_K7y_S_M2_vrr;
      Double I_NAI_M6y3z_S_M1_vrr = PAZ*I_NAI_L6y2z_S_M1_vrr-PNZ*I_NAI_L6y2z_S_M2_vrr+2*oned2z*I_NAI_K6yz_S_M1_vrr-2*oned2z*I_NAI_K6yz_S_M2_vrr;
      Double I_NAI_M5y4z_S_M1_vrr = PAZ*I_NAI_L5y3z_S_M1_vrr-PNZ*I_NAI_L5y3z_S_M2_vrr+3*oned2z*I_NAI_K5y2z_S_M1_vrr-3*oned2z*I_NAI_K5y2z_S_M2_vrr;
      Double I_NAI_M4y5z_S_M1_vrr = PAY*I_NAI_L3y5z_S_M1_vrr-PNY*I_NAI_L3y5z_S_M2_vrr+3*oned2z*I_NAI_K2y5z_S_M1_vrr-3*oned2z*I_NAI_K2y5z_S_M2_vrr;
      Double I_NAI_M3y6z_S_M1_vrr = PAY*I_NAI_L2y6z_S_M1_vrr-PNY*I_NAI_L2y6z_S_M2_vrr+2*oned2z*I_NAI_Ky6z_S_M1_vrr-2*oned2z*I_NAI_Ky6z_S_M2_vrr;
      Double I_NAI_M2y7z_S_M1_vrr = PAY*I_NAI_Ly7z_S_M1_vrr-PNY*I_NAI_Ly7z_S_M2_vrr+oned2z*I_NAI_K7z_S_M1_vrr-oned2z*I_NAI_K7z_S_M2_vrr;
      Double I_NAI_My8z_S_M1_vrr = PAY*I_NAI_L8z_S_M1_vrr-PNY*I_NAI_L8z_S_M2_vrr;
      Double I_NAI_M9z_S_M1_vrr = PAZ*I_NAI_L8z_S_M1_vrr-PNZ*I_NAI_L8z_S_M2_vrr+8*oned2z*I_NAI_K7z_S_M1_vrr-8*oned2z*I_NAI_K7z_S_M2_vrr;

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
       * shell quartet name: SQ_NAI_N_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_M_S
       * RHS shell quartet name: SQ_NAI_M_S_M1
       * RHS shell quartet name: SQ_NAI_L_S
       * RHS shell quartet name: SQ_NAI_L_S_M1
       ************************************************************/
      Double I_NAI_N10x_S_vrr = PAX*I_NAI_M9x_S_vrr-PNX*I_NAI_M9x_S_M1_vrr+9*oned2z*I_NAI_L8x_S_vrr-9*oned2z*I_NAI_L8x_S_M1_vrr;
      Double I_NAI_N9xy_S_vrr = PAY*I_NAI_M9x_S_vrr-PNY*I_NAI_M9x_S_M1_vrr;
      Double I_NAI_N9xz_S_vrr = PAZ*I_NAI_M9x_S_vrr-PNZ*I_NAI_M9x_S_M1_vrr;
      Double I_NAI_N8x2y_S_vrr = PAY*I_NAI_M8xy_S_vrr-PNY*I_NAI_M8xy_S_M1_vrr+oned2z*I_NAI_L8x_S_vrr-oned2z*I_NAI_L8x_S_M1_vrr;
      Double I_NAI_N8xyz_S_vrr = PAZ*I_NAI_M8xy_S_vrr-PNZ*I_NAI_M8xy_S_M1_vrr;
      Double I_NAI_N8x2z_S_vrr = PAZ*I_NAI_M8xz_S_vrr-PNZ*I_NAI_M8xz_S_M1_vrr+oned2z*I_NAI_L8x_S_vrr-oned2z*I_NAI_L8x_S_M1_vrr;
      Double I_NAI_N7x3y_S_vrr = PAY*I_NAI_M7x2y_S_vrr-PNY*I_NAI_M7x2y_S_M1_vrr+2*oned2z*I_NAI_L7xy_S_vrr-2*oned2z*I_NAI_L7xy_S_M1_vrr;
      Double I_NAI_N7x2yz_S_vrr = PAZ*I_NAI_M7x2y_S_vrr-PNZ*I_NAI_M7x2y_S_M1_vrr;
      Double I_NAI_N7xy2z_S_vrr = PAY*I_NAI_M7x2z_S_vrr-PNY*I_NAI_M7x2z_S_M1_vrr;
      Double I_NAI_N7x3z_S_vrr = PAZ*I_NAI_M7x2z_S_vrr-PNZ*I_NAI_M7x2z_S_M1_vrr+2*oned2z*I_NAI_L7xz_S_vrr-2*oned2z*I_NAI_L7xz_S_M1_vrr;
      Double I_NAI_N6x4y_S_vrr = PAY*I_NAI_M6x3y_S_vrr-PNY*I_NAI_M6x3y_S_M1_vrr+3*oned2z*I_NAI_L6x2y_S_vrr-3*oned2z*I_NAI_L6x2y_S_M1_vrr;
      Double I_NAI_N6x3yz_S_vrr = PAZ*I_NAI_M6x3y_S_vrr-PNZ*I_NAI_M6x3y_S_M1_vrr;
      Double I_NAI_N6x2y2z_S_vrr = PAZ*I_NAI_M6x2yz_S_vrr-PNZ*I_NAI_M6x2yz_S_M1_vrr+oned2z*I_NAI_L6x2y_S_vrr-oned2z*I_NAI_L6x2y_S_M1_vrr;
      Double I_NAI_N6xy3z_S_vrr = PAY*I_NAI_M6x3z_S_vrr-PNY*I_NAI_M6x3z_S_M1_vrr;
      Double I_NAI_N6x4z_S_vrr = PAZ*I_NAI_M6x3z_S_vrr-PNZ*I_NAI_M6x3z_S_M1_vrr+3*oned2z*I_NAI_L6x2z_S_vrr-3*oned2z*I_NAI_L6x2z_S_M1_vrr;
      Double I_NAI_N5x5y_S_vrr = PAY*I_NAI_M5x4y_S_vrr-PNY*I_NAI_M5x4y_S_M1_vrr+4*oned2z*I_NAI_L5x3y_S_vrr-4*oned2z*I_NAI_L5x3y_S_M1_vrr;
      Double I_NAI_N5x4yz_S_vrr = PAZ*I_NAI_M5x4y_S_vrr-PNZ*I_NAI_M5x4y_S_M1_vrr;
      Double I_NAI_N5x3y2z_S_vrr = PAZ*I_NAI_M5x3yz_S_vrr-PNZ*I_NAI_M5x3yz_S_M1_vrr+oned2z*I_NAI_L5x3y_S_vrr-oned2z*I_NAI_L5x3y_S_M1_vrr;
      Double I_NAI_N5x2y3z_S_vrr = PAY*I_NAI_M5xy3z_S_vrr-PNY*I_NAI_M5xy3z_S_M1_vrr+oned2z*I_NAI_L5x3z_S_vrr-oned2z*I_NAI_L5x3z_S_M1_vrr;
      Double I_NAI_N5xy4z_S_vrr = PAY*I_NAI_M5x4z_S_vrr-PNY*I_NAI_M5x4z_S_M1_vrr;
      Double I_NAI_N5x5z_S_vrr = PAZ*I_NAI_M5x4z_S_vrr-PNZ*I_NAI_M5x4z_S_M1_vrr+4*oned2z*I_NAI_L5x3z_S_vrr-4*oned2z*I_NAI_L5x3z_S_M1_vrr;
      Double I_NAI_N4x6y_S_vrr = PAX*I_NAI_M3x6y_S_vrr-PNX*I_NAI_M3x6y_S_M1_vrr+3*oned2z*I_NAI_L2x6y_S_vrr-3*oned2z*I_NAI_L2x6y_S_M1_vrr;
      Double I_NAI_N4x5yz_S_vrr = PAZ*I_NAI_M4x5y_S_vrr-PNZ*I_NAI_M4x5y_S_M1_vrr;
      Double I_NAI_N4x4y2z_S_vrr = PAZ*I_NAI_M4x4yz_S_vrr-PNZ*I_NAI_M4x4yz_S_M1_vrr+oned2z*I_NAI_L4x4y_S_vrr-oned2z*I_NAI_L4x4y_S_M1_vrr;
      Double I_NAI_N4x3y3z_S_vrr = PAZ*I_NAI_M4x3y2z_S_vrr-PNZ*I_NAI_M4x3y2z_S_M1_vrr+2*oned2z*I_NAI_L4x3yz_S_vrr-2*oned2z*I_NAI_L4x3yz_S_M1_vrr;
      Double I_NAI_N4x2y4z_S_vrr = PAY*I_NAI_M4xy4z_S_vrr-PNY*I_NAI_M4xy4z_S_M1_vrr+oned2z*I_NAI_L4x4z_S_vrr-oned2z*I_NAI_L4x4z_S_M1_vrr;
      Double I_NAI_N4xy5z_S_vrr = PAY*I_NAI_M4x5z_S_vrr-PNY*I_NAI_M4x5z_S_M1_vrr;
      Double I_NAI_N4x6z_S_vrr = PAX*I_NAI_M3x6z_S_vrr-PNX*I_NAI_M3x6z_S_M1_vrr+3*oned2z*I_NAI_L2x6z_S_vrr-3*oned2z*I_NAI_L2x6z_S_M1_vrr;
      Double I_NAI_N3x7y_S_vrr = PAX*I_NAI_M2x7y_S_vrr-PNX*I_NAI_M2x7y_S_M1_vrr+2*oned2z*I_NAI_Lx7y_S_vrr-2*oned2z*I_NAI_Lx7y_S_M1_vrr;
      Double I_NAI_N3x6yz_S_vrr = PAZ*I_NAI_M3x6y_S_vrr-PNZ*I_NAI_M3x6y_S_M1_vrr;
      Double I_NAI_N3x5y2z_S_vrr = PAZ*I_NAI_M3x5yz_S_vrr-PNZ*I_NAI_M3x5yz_S_M1_vrr+oned2z*I_NAI_L3x5y_S_vrr-oned2z*I_NAI_L3x5y_S_M1_vrr;
      Double I_NAI_N3x4y3z_S_vrr = PAZ*I_NAI_M3x4y2z_S_vrr-PNZ*I_NAI_M3x4y2z_S_M1_vrr+2*oned2z*I_NAI_L3x4yz_S_vrr-2*oned2z*I_NAI_L3x4yz_S_M1_vrr;
      Double I_NAI_N3x3y4z_S_vrr = PAY*I_NAI_M3x2y4z_S_vrr-PNY*I_NAI_M3x2y4z_S_M1_vrr+2*oned2z*I_NAI_L3xy4z_S_vrr-2*oned2z*I_NAI_L3xy4z_S_M1_vrr;
      Double I_NAI_N3x2y5z_S_vrr = PAY*I_NAI_M3xy5z_S_vrr-PNY*I_NAI_M3xy5z_S_M1_vrr+oned2z*I_NAI_L3x5z_S_vrr-oned2z*I_NAI_L3x5z_S_M1_vrr;
      Double I_NAI_N3xy6z_S_vrr = PAY*I_NAI_M3x6z_S_vrr-PNY*I_NAI_M3x6z_S_M1_vrr;
      Double I_NAI_N3x7z_S_vrr = PAX*I_NAI_M2x7z_S_vrr-PNX*I_NAI_M2x7z_S_M1_vrr+2*oned2z*I_NAI_Lx7z_S_vrr-2*oned2z*I_NAI_Lx7z_S_M1_vrr;
      Double I_NAI_N2x8y_S_vrr = PAX*I_NAI_Mx8y_S_vrr-PNX*I_NAI_Mx8y_S_M1_vrr+oned2z*I_NAI_L8y_S_vrr-oned2z*I_NAI_L8y_S_M1_vrr;
      Double I_NAI_N2x7yz_S_vrr = PAZ*I_NAI_M2x7y_S_vrr-PNZ*I_NAI_M2x7y_S_M1_vrr;
      Double I_NAI_N2x6y2z_S_vrr = PAZ*I_NAI_M2x6yz_S_vrr-PNZ*I_NAI_M2x6yz_S_M1_vrr+oned2z*I_NAI_L2x6y_S_vrr-oned2z*I_NAI_L2x6y_S_M1_vrr;
      Double I_NAI_N2x5y3z_S_vrr = PAX*I_NAI_Mx5y3z_S_vrr-PNX*I_NAI_Mx5y3z_S_M1_vrr+oned2z*I_NAI_L5y3z_S_vrr-oned2z*I_NAI_L5y3z_S_M1_vrr;
      Double I_NAI_N2x4y4z_S_vrr = PAX*I_NAI_Mx4y4z_S_vrr-PNX*I_NAI_Mx4y4z_S_M1_vrr+oned2z*I_NAI_L4y4z_S_vrr-oned2z*I_NAI_L4y4z_S_M1_vrr;
      Double I_NAI_N2x3y5z_S_vrr = PAX*I_NAI_Mx3y5z_S_vrr-PNX*I_NAI_Mx3y5z_S_M1_vrr+oned2z*I_NAI_L3y5z_S_vrr-oned2z*I_NAI_L3y5z_S_M1_vrr;
      Double I_NAI_N2x2y6z_S_vrr = PAY*I_NAI_M2xy6z_S_vrr-PNY*I_NAI_M2xy6z_S_M1_vrr+oned2z*I_NAI_L2x6z_S_vrr-oned2z*I_NAI_L2x6z_S_M1_vrr;
      Double I_NAI_N2xy7z_S_vrr = PAY*I_NAI_M2x7z_S_vrr-PNY*I_NAI_M2x7z_S_M1_vrr;
      Double I_NAI_N2x8z_S_vrr = PAX*I_NAI_Mx8z_S_vrr-PNX*I_NAI_Mx8z_S_M1_vrr+oned2z*I_NAI_L8z_S_vrr-oned2z*I_NAI_L8z_S_M1_vrr;
      Double I_NAI_Nx9y_S_vrr = PAX*I_NAI_M9y_S_vrr-PNX*I_NAI_M9y_S_M1_vrr;
      Double I_NAI_Nx8yz_S_vrr = PAZ*I_NAI_Mx8y_S_vrr-PNZ*I_NAI_Mx8y_S_M1_vrr;
      Double I_NAI_Nx7y2z_S_vrr = PAX*I_NAI_M7y2z_S_vrr-PNX*I_NAI_M7y2z_S_M1_vrr;
      Double I_NAI_Nx6y3z_S_vrr = PAX*I_NAI_M6y3z_S_vrr-PNX*I_NAI_M6y3z_S_M1_vrr;
      Double I_NAI_Nx5y4z_S_vrr = PAX*I_NAI_M5y4z_S_vrr-PNX*I_NAI_M5y4z_S_M1_vrr;
      Double I_NAI_Nx4y5z_S_vrr = PAX*I_NAI_M4y5z_S_vrr-PNX*I_NAI_M4y5z_S_M1_vrr;
      Double I_NAI_Nx3y6z_S_vrr = PAX*I_NAI_M3y6z_S_vrr-PNX*I_NAI_M3y6z_S_M1_vrr;
      Double I_NAI_Nx2y7z_S_vrr = PAX*I_NAI_M2y7z_S_vrr-PNX*I_NAI_M2y7z_S_M1_vrr;
      Double I_NAI_Nxy8z_S_vrr = PAY*I_NAI_Mx8z_S_vrr-PNY*I_NAI_Mx8z_S_M1_vrr;
      Double I_NAI_Nx9z_S_vrr = PAX*I_NAI_M9z_S_vrr-PNX*I_NAI_M9z_S_M1_vrr;
      Double I_NAI_N10y_S_vrr = PAY*I_NAI_M9y_S_vrr-PNY*I_NAI_M9y_S_M1_vrr+9*oned2z*I_NAI_L8y_S_vrr-9*oned2z*I_NAI_L8y_S_M1_vrr;
      Double I_NAI_N9yz_S_vrr = PAZ*I_NAI_M9y_S_vrr-PNZ*I_NAI_M9y_S_M1_vrr;
      Double I_NAI_N8y2z_S_vrr = PAZ*I_NAI_M8yz_S_vrr-PNZ*I_NAI_M8yz_S_M1_vrr+oned2z*I_NAI_L8y_S_vrr-oned2z*I_NAI_L8y_S_M1_vrr;
      Double I_NAI_N7y3z_S_vrr = PAZ*I_NAI_M7y2z_S_vrr-PNZ*I_NAI_M7y2z_S_M1_vrr+2*oned2z*I_NAI_L7yz_S_vrr-2*oned2z*I_NAI_L7yz_S_M1_vrr;
      Double I_NAI_N6y4z_S_vrr = PAZ*I_NAI_M6y3z_S_vrr-PNZ*I_NAI_M6y3z_S_M1_vrr+3*oned2z*I_NAI_L6y2z_S_vrr-3*oned2z*I_NAI_L6y2z_S_M1_vrr;
      Double I_NAI_N5y5z_S_vrr = PAZ*I_NAI_M5y4z_S_vrr-PNZ*I_NAI_M5y4z_S_M1_vrr+4*oned2z*I_NAI_L5y3z_S_vrr-4*oned2z*I_NAI_L5y3z_S_M1_vrr;
      Double I_NAI_N4y6z_S_vrr = PAY*I_NAI_M3y6z_S_vrr-PNY*I_NAI_M3y6z_S_M1_vrr+3*oned2z*I_NAI_L2y6z_S_vrr-3*oned2z*I_NAI_L2y6z_S_M1_vrr;
      Double I_NAI_N3y7z_S_vrr = PAY*I_NAI_M2y7z_S_vrr-PNY*I_NAI_M2y7z_S_M1_vrr+2*oned2z*I_NAI_Ly7z_S_vrr-2*oned2z*I_NAI_Ly7z_S_M1_vrr;
      Double I_NAI_N2y8z_S_vrr = PAY*I_NAI_My8z_S_vrr-PNY*I_NAI_My8z_S_M1_vrr+oned2z*I_NAI_L8z_S_vrr-oned2z*I_NAI_L8z_S_M1_vrr;
      Double I_NAI_Ny9z_S_vrr = PAY*I_NAI_M9z_S_vrr-PNY*I_NAI_M9z_S_M1_vrr;
      Double I_NAI_N10z_S_vrr = PAZ*I_NAI_M9z_S_vrr-PNZ*I_NAI_M9z_S_M1_vrr+9*oned2z*I_NAI_L8z_S_vrr-9*oned2z*I_NAI_L8z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_L_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_L8x_S += I_NAI_L8x_S_vrr;
      I_NAI_L7xy_S += I_NAI_L7xy_S_vrr;
      I_NAI_L7xz_S += I_NAI_L7xz_S_vrr;
      I_NAI_L6x2y_S += I_NAI_L6x2y_S_vrr;
      I_NAI_L6xyz_S += I_NAI_L6xyz_S_vrr;
      I_NAI_L6x2z_S += I_NAI_L6x2z_S_vrr;
      I_NAI_L5x3y_S += I_NAI_L5x3y_S_vrr;
      I_NAI_L5x2yz_S += I_NAI_L5x2yz_S_vrr;
      I_NAI_L5xy2z_S += I_NAI_L5xy2z_S_vrr;
      I_NAI_L5x3z_S += I_NAI_L5x3z_S_vrr;
      I_NAI_L4x4y_S += I_NAI_L4x4y_S_vrr;
      I_NAI_L4x3yz_S += I_NAI_L4x3yz_S_vrr;
      I_NAI_L4x2y2z_S += I_NAI_L4x2y2z_S_vrr;
      I_NAI_L4xy3z_S += I_NAI_L4xy3z_S_vrr;
      I_NAI_L4x4z_S += I_NAI_L4x4z_S_vrr;
      I_NAI_L3x5y_S += I_NAI_L3x5y_S_vrr;
      I_NAI_L3x4yz_S += I_NAI_L3x4yz_S_vrr;
      I_NAI_L3x3y2z_S += I_NAI_L3x3y2z_S_vrr;
      I_NAI_L3x2y3z_S += I_NAI_L3x2y3z_S_vrr;
      I_NAI_L3xy4z_S += I_NAI_L3xy4z_S_vrr;
      I_NAI_L3x5z_S += I_NAI_L3x5z_S_vrr;
      I_NAI_L2x6y_S += I_NAI_L2x6y_S_vrr;
      I_NAI_L2x5yz_S += I_NAI_L2x5yz_S_vrr;
      I_NAI_L2x4y2z_S += I_NAI_L2x4y2z_S_vrr;
      I_NAI_L2x3y3z_S += I_NAI_L2x3y3z_S_vrr;
      I_NAI_L2x2y4z_S += I_NAI_L2x2y4z_S_vrr;
      I_NAI_L2xy5z_S += I_NAI_L2xy5z_S_vrr;
      I_NAI_L2x6z_S += I_NAI_L2x6z_S_vrr;
      I_NAI_Lx7y_S += I_NAI_Lx7y_S_vrr;
      I_NAI_Lx6yz_S += I_NAI_Lx6yz_S_vrr;
      I_NAI_Lx5y2z_S += I_NAI_Lx5y2z_S_vrr;
      I_NAI_Lx4y3z_S += I_NAI_Lx4y3z_S_vrr;
      I_NAI_Lx3y4z_S += I_NAI_Lx3y4z_S_vrr;
      I_NAI_Lx2y5z_S += I_NAI_Lx2y5z_S_vrr;
      I_NAI_Lxy6z_S += I_NAI_Lxy6z_S_vrr;
      I_NAI_Lx7z_S += I_NAI_Lx7z_S_vrr;
      I_NAI_L8y_S += I_NAI_L8y_S_vrr;
      I_NAI_L7yz_S += I_NAI_L7yz_S_vrr;
      I_NAI_L6y2z_S += I_NAI_L6y2z_S_vrr;
      I_NAI_L5y3z_S += I_NAI_L5y3z_S_vrr;
      I_NAI_L4y4z_S += I_NAI_L4y4z_S_vrr;
      I_NAI_L3y5z_S += I_NAI_L3y5z_S_vrr;
      I_NAI_L2y6z_S += I_NAI_L2y6z_S_vrr;
      I_NAI_Ly7z_S += I_NAI_Ly7z_S_vrr;
      I_NAI_L8z_S += I_NAI_L8z_S_vrr;

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
       * shell quartet name: SQ_NAI_N_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_N_S_a_coefs = alpha;
      I_NAI_N10x_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N10x_S_vrr;
      I_NAI_N9xy_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N9xy_S_vrr;
      I_NAI_N9xz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N9xz_S_vrr;
      I_NAI_N8x2y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N8x2y_S_vrr;
      I_NAI_N8xyz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N8xyz_S_vrr;
      I_NAI_N8x2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N8x2z_S_vrr;
      I_NAI_N7x3y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N7x3y_S_vrr;
      I_NAI_N7x2yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N7x2yz_S_vrr;
      I_NAI_N7xy2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N7xy2z_S_vrr;
      I_NAI_N7x3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N7x3z_S_vrr;
      I_NAI_N6x4y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6x4y_S_vrr;
      I_NAI_N6x3yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6x3yz_S_vrr;
      I_NAI_N6x2y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6x2y2z_S_vrr;
      I_NAI_N6xy3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6xy3z_S_vrr;
      I_NAI_N6x4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6x4z_S_vrr;
      I_NAI_N5x5y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5x5y_S_vrr;
      I_NAI_N5x4yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5x4yz_S_vrr;
      I_NAI_N5x3y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5x3y2z_S_vrr;
      I_NAI_N5x2y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5x2y3z_S_vrr;
      I_NAI_N5xy4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5xy4z_S_vrr;
      I_NAI_N5x5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5x5z_S_vrr;
      I_NAI_N4x6y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x6y_S_vrr;
      I_NAI_N4x5yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x5yz_S_vrr;
      I_NAI_N4x4y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x4y2z_S_vrr;
      I_NAI_N4x3y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x3y3z_S_vrr;
      I_NAI_N4x2y4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x2y4z_S_vrr;
      I_NAI_N4xy5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4xy5z_S_vrr;
      I_NAI_N4x6z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4x6z_S_vrr;
      I_NAI_N3x7y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x7y_S_vrr;
      I_NAI_N3x6yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x6yz_S_vrr;
      I_NAI_N3x5y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x5y2z_S_vrr;
      I_NAI_N3x4y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x4y3z_S_vrr;
      I_NAI_N3x3y4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x3y4z_S_vrr;
      I_NAI_N3x2y5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x2y5z_S_vrr;
      I_NAI_N3xy6z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3xy6z_S_vrr;
      I_NAI_N3x7z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3x7z_S_vrr;
      I_NAI_N2x8y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x8y_S_vrr;
      I_NAI_N2x7yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x7yz_S_vrr;
      I_NAI_N2x6y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x6y2z_S_vrr;
      I_NAI_N2x5y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x5y3z_S_vrr;
      I_NAI_N2x4y4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x4y4z_S_vrr;
      I_NAI_N2x3y5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x3y5z_S_vrr;
      I_NAI_N2x2y6z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x2y6z_S_vrr;
      I_NAI_N2xy7z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2xy7z_S_vrr;
      I_NAI_N2x8z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2x8z_S_vrr;
      I_NAI_Nx9y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx9y_S_vrr;
      I_NAI_Nx8yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx8yz_S_vrr;
      I_NAI_Nx7y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx7y2z_S_vrr;
      I_NAI_Nx6y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx6y3z_S_vrr;
      I_NAI_Nx5y4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx5y4z_S_vrr;
      I_NAI_Nx4y5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx4y5z_S_vrr;
      I_NAI_Nx3y6z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx3y6z_S_vrr;
      I_NAI_Nx2y7z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx2y7z_S_vrr;
      I_NAI_Nxy8z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nxy8z_S_vrr;
      I_NAI_Nx9z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Nx9z_S_vrr;
      I_NAI_N10y_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N10y_S_vrr;
      I_NAI_N9yz_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N9yz_S_vrr;
      I_NAI_N8y2z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N8y2z_S_vrr;
      I_NAI_N7y3z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N7y3z_S_vrr;
      I_NAI_N6y4z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N6y4z_S_vrr;
      I_NAI_N5y5z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N5y5z_S_vrr;
      I_NAI_N4y6z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N4y6z_S_vrr;
      I_NAI_N3y7z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N3y7z_S_vrr;
      I_NAI_N2y8z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N2y8z_S_vrr;
      I_NAI_Ny9z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_Ny9z_S_vrr;
      I_NAI_N10z_S_a += SQ_NAI_N_S_a_coefs*I_NAI_N10z_S_vrr;

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
       * shell quartet name: SQ_NAI_N_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_N_S_b_coefs = beta;
      I_NAI_N10x_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N10x_S_vrr;
      I_NAI_N9xy_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N9xy_S_vrr;
      I_NAI_N9xz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N9xz_S_vrr;
      I_NAI_N8x2y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N8x2y_S_vrr;
      I_NAI_N8xyz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N8xyz_S_vrr;
      I_NAI_N8x2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N8x2z_S_vrr;
      I_NAI_N7x3y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N7x3y_S_vrr;
      I_NAI_N7x2yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N7x2yz_S_vrr;
      I_NAI_N7xy2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N7xy2z_S_vrr;
      I_NAI_N7x3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N7x3z_S_vrr;
      I_NAI_N6x4y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6x4y_S_vrr;
      I_NAI_N6x3yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6x3yz_S_vrr;
      I_NAI_N6x2y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6x2y2z_S_vrr;
      I_NAI_N6xy3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6xy3z_S_vrr;
      I_NAI_N6x4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6x4z_S_vrr;
      I_NAI_N5x5y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5x5y_S_vrr;
      I_NAI_N5x4yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5x4yz_S_vrr;
      I_NAI_N5x3y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5x3y2z_S_vrr;
      I_NAI_N5x2y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5x2y3z_S_vrr;
      I_NAI_N5xy4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5xy4z_S_vrr;
      I_NAI_N5x5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5x5z_S_vrr;
      I_NAI_N4x6y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x6y_S_vrr;
      I_NAI_N4x5yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x5yz_S_vrr;
      I_NAI_N4x4y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x4y2z_S_vrr;
      I_NAI_N4x3y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x3y3z_S_vrr;
      I_NAI_N4x2y4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x2y4z_S_vrr;
      I_NAI_N4xy5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4xy5z_S_vrr;
      I_NAI_N4x6z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4x6z_S_vrr;
      I_NAI_N3x7y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x7y_S_vrr;
      I_NAI_N3x6yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x6yz_S_vrr;
      I_NAI_N3x5y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x5y2z_S_vrr;
      I_NAI_N3x4y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x4y3z_S_vrr;
      I_NAI_N3x3y4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x3y4z_S_vrr;
      I_NAI_N3x2y5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x2y5z_S_vrr;
      I_NAI_N3xy6z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3xy6z_S_vrr;
      I_NAI_N3x7z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3x7z_S_vrr;
      I_NAI_N2x8y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x8y_S_vrr;
      I_NAI_N2x7yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x7yz_S_vrr;
      I_NAI_N2x6y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x6y2z_S_vrr;
      I_NAI_N2x5y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x5y3z_S_vrr;
      I_NAI_N2x4y4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x4y4z_S_vrr;
      I_NAI_N2x3y5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x3y5z_S_vrr;
      I_NAI_N2x2y6z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x2y6z_S_vrr;
      I_NAI_N2xy7z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2xy7z_S_vrr;
      I_NAI_N2x8z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2x8z_S_vrr;
      I_NAI_Nx9y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx9y_S_vrr;
      I_NAI_Nx8yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx8yz_S_vrr;
      I_NAI_Nx7y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx7y2z_S_vrr;
      I_NAI_Nx6y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx6y3z_S_vrr;
      I_NAI_Nx5y4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx5y4z_S_vrr;
      I_NAI_Nx4y5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx4y5z_S_vrr;
      I_NAI_Nx3y6z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx3y6z_S_vrr;
      I_NAI_Nx2y7z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx2y7z_S_vrr;
      I_NAI_Nxy8z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nxy8z_S_vrr;
      I_NAI_Nx9z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Nx9z_S_vrr;
      I_NAI_N10y_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N10y_S_vrr;
      I_NAI_N9yz_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N9yz_S_vrr;
      I_NAI_N8y2z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N8y2z_S_vrr;
      I_NAI_N7y3z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N7y3z_S_vrr;
      I_NAI_N6y4z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N6y4z_S_vrr;
      I_NAI_N5y5z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N5y5z_S_vrr;
      I_NAI_N4y6z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N4y6z_S_vrr;
      I_NAI_N3y7z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N3y7z_S_vrr;
      I_NAI_N2y8z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N2y8z_S_vrr;
      I_NAI_Ny9z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_Ny9z_S_vrr;
      I_NAI_N10z_S_b += SQ_NAI_N_S_b_coefs*I_NAI_N10z_S_vrr;

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
   * totally 45 integrals are omitted 
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
   * totally 0 integrals are omitted 
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
  Double I_NAI_I6x_Py = I_NAI_K6xy_S+ABY*I_NAI_I6x_S;
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
  Double I_NAI_I6x_Pz = I_NAI_K6xz_S+ABZ*I_NAI_I6x_S;
  Double I_NAI_I5xy_Pz = I_NAI_K5xyz_S+ABZ*I_NAI_I5xy_S;
  Double I_NAI_I5xz_Pz = I_NAI_K5x2z_S+ABZ*I_NAI_I5xz_S;
  Double I_NAI_I4x2y_Pz = I_NAI_K4x2yz_S+ABZ*I_NAI_I4x2y_S;
  Double I_NAI_I4xyz_Pz = I_NAI_K4xy2z_S+ABZ*I_NAI_I4xyz_S;
  Double I_NAI_I4x2z_Pz = I_NAI_K4x3z_S+ABZ*I_NAI_I4x2z_S;
  Double I_NAI_I3x3y_Pz = I_NAI_K3x3yz_S+ABZ*I_NAI_I3x3y_S;
  Double I_NAI_I3x2yz_Pz = I_NAI_K3x2y2z_S+ABZ*I_NAI_I3x2yz_S;
  Double I_NAI_I3xy2z_Pz = I_NAI_K3xy3z_S+ABZ*I_NAI_I3xy2z_S;
  Double I_NAI_I3x3z_Pz = I_NAI_K3x4z_S+ABZ*I_NAI_I3x3z_S;
  Double I_NAI_I2x4y_Pz = I_NAI_K2x4yz_S+ABZ*I_NAI_I2x4y_S;
  Double I_NAI_I2x3yz_Pz = I_NAI_K2x3y2z_S+ABZ*I_NAI_I2x3yz_S;
  Double I_NAI_I2x2y2z_Pz = I_NAI_K2x2y3z_S+ABZ*I_NAI_I2x2y2z_S;
  Double I_NAI_I2xy3z_Pz = I_NAI_K2xy4z_S+ABZ*I_NAI_I2xy3z_S;
  Double I_NAI_I2x4z_Pz = I_NAI_K2x5z_S+ABZ*I_NAI_I2x4z_S;
  Double I_NAI_Ix5y_Pz = I_NAI_Kx5yz_S+ABZ*I_NAI_Ix5y_S;
  Double I_NAI_Ix4yz_Pz = I_NAI_Kx4y2z_S+ABZ*I_NAI_Ix4yz_S;
  Double I_NAI_Ix3y2z_Pz = I_NAI_Kx3y3z_S+ABZ*I_NAI_Ix3y2z_S;
  Double I_NAI_Ix2y3z_Pz = I_NAI_Kx2y4z_S+ABZ*I_NAI_Ix2y3z_S;
  Double I_NAI_Ixy4z_Pz = I_NAI_Kxy5z_S+ABZ*I_NAI_Ixy4z_S;
  Double I_NAI_Ix5z_Pz = I_NAI_Kx6z_S+ABZ*I_NAI_Ix5z_S;
  Double I_NAI_I6y_Pz = I_NAI_K6yz_S+ABZ*I_NAI_I6y_S;
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
   * totally 42 integrals are omitted 
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
   * totally 30 integrals are omitted 
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
   * shell quartet name: SQ_NAI_K_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_S
   * RHS shell quartet name: SQ_NAI_K_S
   ************************************************************/
  Double I_NAI_K7x_Px = I_NAI_L8x_S+ABX*I_NAI_K7x_S;
  Double I_NAI_K6xy_Px = I_NAI_L7xy_S+ABX*I_NAI_K6xy_S;
  Double I_NAI_K6xz_Px = I_NAI_L7xz_S+ABX*I_NAI_K6xz_S;
  Double I_NAI_K5x2y_Px = I_NAI_L6x2y_S+ABX*I_NAI_K5x2y_S;
  Double I_NAI_K5xyz_Px = I_NAI_L6xyz_S+ABX*I_NAI_K5xyz_S;
  Double I_NAI_K5x2z_Px = I_NAI_L6x2z_S+ABX*I_NAI_K5x2z_S;
  Double I_NAI_K4x3y_Px = I_NAI_L5x3y_S+ABX*I_NAI_K4x3y_S;
  Double I_NAI_K4x2yz_Px = I_NAI_L5x2yz_S+ABX*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Px = I_NAI_L5xy2z_S+ABX*I_NAI_K4xy2z_S;
  Double I_NAI_K4x3z_Px = I_NAI_L5x3z_S+ABX*I_NAI_K4x3z_S;
  Double I_NAI_K3x4y_Px = I_NAI_L4x4y_S+ABX*I_NAI_K3x4y_S;
  Double I_NAI_K3x3yz_Px = I_NAI_L4x3yz_S+ABX*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Px = I_NAI_L4x2y2z_S+ABX*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Px = I_NAI_L4xy3z_S+ABX*I_NAI_K3xy3z_S;
  Double I_NAI_K3x4z_Px = I_NAI_L4x4z_S+ABX*I_NAI_K3x4z_S;
  Double I_NAI_K2x5y_Px = I_NAI_L3x5y_S+ABX*I_NAI_K2x5y_S;
  Double I_NAI_K2x4yz_Px = I_NAI_L3x4yz_S+ABX*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Px = I_NAI_L3x3y2z_S+ABX*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Px = I_NAI_L3x2y3z_S+ABX*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Px = I_NAI_L3xy4z_S+ABX*I_NAI_K2xy4z_S;
  Double I_NAI_K2x5z_Px = I_NAI_L3x5z_S+ABX*I_NAI_K2x5z_S;
  Double I_NAI_Kx6y_Px = I_NAI_L2x6y_S+ABX*I_NAI_Kx6y_S;
  Double I_NAI_Kx5yz_Px = I_NAI_L2x5yz_S+ABX*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Px = I_NAI_L2x4y2z_S+ABX*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Px = I_NAI_L2x3y3z_S+ABX*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Px = I_NAI_L2x2y4z_S+ABX*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Px = I_NAI_L2xy5z_S+ABX*I_NAI_Kxy5z_S;
  Double I_NAI_Kx6z_Px = I_NAI_L2x6z_S+ABX*I_NAI_Kx6z_S;
  Double I_NAI_K6yz_Px = I_NAI_Lx6yz_S+ABX*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Px = I_NAI_Lx5y2z_S+ABX*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Px = I_NAI_Lx4y3z_S+ABX*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Px = I_NAI_Lx3y4z_S+ABX*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Px = I_NAI_Lx2y5z_S+ABX*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Px = I_NAI_Lxy6z_S+ABX*I_NAI_Ky6z_S;
  Double I_NAI_K6xy_Py = I_NAI_L6x2y_S+ABY*I_NAI_K6xy_S;
  Double I_NAI_K5x2y_Py = I_NAI_L5x3y_S+ABY*I_NAI_K5x2y_S;
  Double I_NAI_K5xyz_Py = I_NAI_L5x2yz_S+ABY*I_NAI_K5xyz_S;
  Double I_NAI_K4x3y_Py = I_NAI_L4x4y_S+ABY*I_NAI_K4x3y_S;
  Double I_NAI_K4x2yz_Py = I_NAI_L4x3yz_S+ABY*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Py = I_NAI_L4x2y2z_S+ABY*I_NAI_K4xy2z_S;
  Double I_NAI_K3x4y_Py = I_NAI_L3x5y_S+ABY*I_NAI_K3x4y_S;
  Double I_NAI_K3x3yz_Py = I_NAI_L3x4yz_S+ABY*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Py = I_NAI_L3x3y2z_S+ABY*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Py = I_NAI_L3x2y3z_S+ABY*I_NAI_K3xy3z_S;
  Double I_NAI_K2x5y_Py = I_NAI_L2x6y_S+ABY*I_NAI_K2x5y_S;
  Double I_NAI_K2x4yz_Py = I_NAI_L2x5yz_S+ABY*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Py = I_NAI_L2x4y2z_S+ABY*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Py = I_NAI_L2x3y3z_S+ABY*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Py = I_NAI_L2x2y4z_S+ABY*I_NAI_K2xy4z_S;
  Double I_NAI_Kx6y_Py = I_NAI_Lx7y_S+ABY*I_NAI_Kx6y_S;
  Double I_NAI_Kx5yz_Py = I_NAI_Lx6yz_S+ABY*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Py = I_NAI_Lx5y2z_S+ABY*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Py = I_NAI_Lx4y3z_S+ABY*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Py = I_NAI_Lx3y4z_S+ABY*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Py = I_NAI_Lx2y5z_S+ABY*I_NAI_Kxy5z_S;
  Double I_NAI_K7y_Py = I_NAI_L8y_S+ABY*I_NAI_K7y_S;
  Double I_NAI_K6yz_Py = I_NAI_L7yz_S+ABY*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Py = I_NAI_L6y2z_S+ABY*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Py = I_NAI_L5y3z_S+ABY*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Py = I_NAI_L4y4z_S+ABY*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Py = I_NAI_L3y5z_S+ABY*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Py = I_NAI_L2y6z_S+ABY*I_NAI_Ky6z_S;
  Double I_NAI_K6xz_Pz = I_NAI_L6x2z_S+ABZ*I_NAI_K6xz_S;
  Double I_NAI_K5xyz_Pz = I_NAI_L5xy2z_S+ABZ*I_NAI_K5xyz_S;
  Double I_NAI_K5x2z_Pz = I_NAI_L5x3z_S+ABZ*I_NAI_K5x2z_S;
  Double I_NAI_K4x2yz_Pz = I_NAI_L4x2y2z_S+ABZ*I_NAI_K4x2yz_S;
  Double I_NAI_K4xy2z_Pz = I_NAI_L4xy3z_S+ABZ*I_NAI_K4xy2z_S;
  Double I_NAI_K4x3z_Pz = I_NAI_L4x4z_S+ABZ*I_NAI_K4x3z_S;
  Double I_NAI_K3x3yz_Pz = I_NAI_L3x3y2z_S+ABZ*I_NAI_K3x3yz_S;
  Double I_NAI_K3x2y2z_Pz = I_NAI_L3x2y3z_S+ABZ*I_NAI_K3x2y2z_S;
  Double I_NAI_K3xy3z_Pz = I_NAI_L3xy4z_S+ABZ*I_NAI_K3xy3z_S;
  Double I_NAI_K3x4z_Pz = I_NAI_L3x5z_S+ABZ*I_NAI_K3x4z_S;
  Double I_NAI_K2x4yz_Pz = I_NAI_L2x4y2z_S+ABZ*I_NAI_K2x4yz_S;
  Double I_NAI_K2x3y2z_Pz = I_NAI_L2x3y3z_S+ABZ*I_NAI_K2x3y2z_S;
  Double I_NAI_K2x2y3z_Pz = I_NAI_L2x2y4z_S+ABZ*I_NAI_K2x2y3z_S;
  Double I_NAI_K2xy4z_Pz = I_NAI_L2xy5z_S+ABZ*I_NAI_K2xy4z_S;
  Double I_NAI_K2x5z_Pz = I_NAI_L2x6z_S+ABZ*I_NAI_K2x5z_S;
  Double I_NAI_Kx5yz_Pz = I_NAI_Lx5y2z_S+ABZ*I_NAI_Kx5yz_S;
  Double I_NAI_Kx4y2z_Pz = I_NAI_Lx4y3z_S+ABZ*I_NAI_Kx4y2z_S;
  Double I_NAI_Kx3y3z_Pz = I_NAI_Lx3y4z_S+ABZ*I_NAI_Kx3y3z_S;
  Double I_NAI_Kx2y4z_Pz = I_NAI_Lx2y5z_S+ABZ*I_NAI_Kx2y4z_S;
  Double I_NAI_Kxy5z_Pz = I_NAI_Lxy6z_S+ABZ*I_NAI_Kxy5z_S;
  Double I_NAI_Kx6z_Pz = I_NAI_Lx7z_S+ABZ*I_NAI_Kx6z_S;
  Double I_NAI_K6yz_Pz = I_NAI_L6y2z_S+ABZ*I_NAI_K6yz_S;
  Double I_NAI_K5y2z_Pz = I_NAI_L5y3z_S+ABZ*I_NAI_K5y2z_S;
  Double I_NAI_K4y3z_Pz = I_NAI_L4y4z_S+ABZ*I_NAI_K4y3z_S;
  Double I_NAI_K3y4z_Pz = I_NAI_L3y5z_S+ABZ*I_NAI_K3y4z_S;
  Double I_NAI_K2y5z_Pz = I_NAI_L2y6z_S+ABZ*I_NAI_K2y5z_S;
  Double I_NAI_Ky6z_Pz = I_NAI_Ly7z_S+ABZ*I_NAI_Ky6z_S;
  Double I_NAI_K7z_Pz = I_NAI_L8z_S+ABZ*I_NAI_K7z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_I_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_P
   * RHS shell quartet name: SQ_NAI_I_P
   ************************************************************/
  Double I_NAI_I6x_D2x = I_NAI_K7x_Px+ABX*I_NAI_I6x_Px;
  Double I_NAI_I5xy_D2x = I_NAI_K6xy_Px+ABX*I_NAI_I5xy_Px;
  Double I_NAI_I5xz_D2x = I_NAI_K6xz_Px+ABX*I_NAI_I5xz_Px;
  Double I_NAI_I4x2y_D2x = I_NAI_K5x2y_Px+ABX*I_NAI_I4x2y_Px;
  Double I_NAI_I4xyz_D2x = I_NAI_K5xyz_Px+ABX*I_NAI_I4xyz_Px;
  Double I_NAI_I4x2z_D2x = I_NAI_K5x2z_Px+ABX*I_NAI_I4x2z_Px;
  Double I_NAI_I3x3y_D2x = I_NAI_K4x3y_Px+ABX*I_NAI_I3x3y_Px;
  Double I_NAI_I3x2yz_D2x = I_NAI_K4x2yz_Px+ABX*I_NAI_I3x2yz_Px;
  Double I_NAI_I3xy2z_D2x = I_NAI_K4xy2z_Px+ABX*I_NAI_I3xy2z_Px;
  Double I_NAI_I3x3z_D2x = I_NAI_K4x3z_Px+ABX*I_NAI_I3x3z_Px;
  Double I_NAI_I2x4y_D2x = I_NAI_K3x4y_Px+ABX*I_NAI_I2x4y_Px;
  Double I_NAI_I2x3yz_D2x = I_NAI_K3x3yz_Px+ABX*I_NAI_I2x3yz_Px;
  Double I_NAI_I2x2y2z_D2x = I_NAI_K3x2y2z_Px+ABX*I_NAI_I2x2y2z_Px;
  Double I_NAI_I2xy3z_D2x = I_NAI_K3xy3z_Px+ABX*I_NAI_I2xy3z_Px;
  Double I_NAI_I2x4z_D2x = I_NAI_K3x4z_Px+ABX*I_NAI_I2x4z_Px;
  Double I_NAI_Ix5y_D2x = I_NAI_K2x5y_Px+ABX*I_NAI_Ix5y_Px;
  Double I_NAI_Ix4yz_D2x = I_NAI_K2x4yz_Px+ABX*I_NAI_Ix4yz_Px;
  Double I_NAI_Ix3y2z_D2x = I_NAI_K2x3y2z_Px+ABX*I_NAI_Ix3y2z_Px;
  Double I_NAI_Ix2y3z_D2x = I_NAI_K2x2y3z_Px+ABX*I_NAI_Ix2y3z_Px;
  Double I_NAI_Ixy4z_D2x = I_NAI_K2xy4z_Px+ABX*I_NAI_Ixy4z_Px;
  Double I_NAI_Ix5z_D2x = I_NAI_K2x5z_Px+ABX*I_NAI_Ix5z_Px;
  Double I_NAI_I6y_D2x = I_NAI_Kx6y_Px+ABX*I_NAI_I6y_Px;
  Double I_NAI_I5yz_D2x = I_NAI_Kx5yz_Px+ABX*I_NAI_I5yz_Px;
  Double I_NAI_I4y2z_D2x = I_NAI_Kx4y2z_Px+ABX*I_NAI_I4y2z_Px;
  Double I_NAI_I3y3z_D2x = I_NAI_Kx3y3z_Px+ABX*I_NAI_I3y3z_Px;
  Double I_NAI_I2y4z_D2x = I_NAI_Kx2y4z_Px+ABX*I_NAI_I2y4z_Px;
  Double I_NAI_Iy5z_D2x = I_NAI_Kxy5z_Px+ABX*I_NAI_Iy5z_Px;
  Double I_NAI_I6z_D2x = I_NAI_Kx6z_Px+ABX*I_NAI_I6z_Px;
  Double I_NAI_I5xz_Dxy = I_NAI_K5xyz_Px+ABY*I_NAI_I5xz_Px;
  Double I_NAI_I4xyz_Dxy = I_NAI_K4x2yz_Px+ABY*I_NAI_I4xyz_Px;
  Double I_NAI_I4x2z_Dxy = I_NAI_K4xy2z_Px+ABY*I_NAI_I4x2z_Px;
  Double I_NAI_I3x2yz_Dxy = I_NAI_K3x3yz_Px+ABY*I_NAI_I3x2yz_Px;
  Double I_NAI_I3xy2z_Dxy = I_NAI_K3x2y2z_Px+ABY*I_NAI_I3xy2z_Px;
  Double I_NAI_I3x3z_Dxy = I_NAI_K3xy3z_Px+ABY*I_NAI_I3x3z_Px;
  Double I_NAI_I2x3yz_Dxy = I_NAI_K2x4yz_Px+ABY*I_NAI_I2x3yz_Px;
  Double I_NAI_I2x2y2z_Dxy = I_NAI_K2x3y2z_Px+ABY*I_NAI_I2x2y2z_Px;
  Double I_NAI_I2xy3z_Dxy = I_NAI_K2x2y3z_Px+ABY*I_NAI_I2xy3z_Px;
  Double I_NAI_I2x4z_Dxy = I_NAI_K2xy4z_Px+ABY*I_NAI_I2x4z_Px;
  Double I_NAI_Ix4yz_Dxy = I_NAI_Kx5yz_Px+ABY*I_NAI_Ix4yz_Px;
  Double I_NAI_Ix3y2z_Dxy = I_NAI_Kx4y2z_Px+ABY*I_NAI_Ix3y2z_Px;
  Double I_NAI_Ix2y3z_Dxy = I_NAI_Kx3y3z_Px+ABY*I_NAI_Ix2y3z_Px;
  Double I_NAI_Ixy4z_Dxy = I_NAI_Kx2y4z_Px+ABY*I_NAI_Ixy4z_Px;
  Double I_NAI_Ix5z_Dxy = I_NAI_Kxy5z_Px+ABY*I_NAI_Ix5z_Px;
  Double I_NAI_I5yz_Dxy = I_NAI_K6yz_Px+ABY*I_NAI_I5yz_Px;
  Double I_NAI_I4y2z_Dxy = I_NAI_K5y2z_Px+ABY*I_NAI_I4y2z_Px;
  Double I_NAI_I3y3z_Dxy = I_NAI_K4y3z_Px+ABY*I_NAI_I3y3z_Px;
  Double I_NAI_I2y4z_Dxy = I_NAI_K3y4z_Px+ABY*I_NAI_I2y4z_Px;
  Double I_NAI_Iy5z_Dxy = I_NAI_K2y5z_Px+ABY*I_NAI_Iy5z_Px;
  Double I_NAI_I6z_Dxy = I_NAI_Ky6z_Px+ABY*I_NAI_I6z_Px;
  Double I_NAI_I6x_D2y = I_NAI_K6xy_Py+ABY*I_NAI_I6x_Py;
  Double I_NAI_I5xy_D2y = I_NAI_K5x2y_Py+ABY*I_NAI_I5xy_Py;
  Double I_NAI_I5xz_D2y = I_NAI_K5xyz_Py+ABY*I_NAI_I5xz_Py;
  Double I_NAI_I4x2y_D2y = I_NAI_K4x3y_Py+ABY*I_NAI_I4x2y_Py;
  Double I_NAI_I4xyz_D2y = I_NAI_K4x2yz_Py+ABY*I_NAI_I4xyz_Py;
  Double I_NAI_I4x2z_D2y = I_NAI_K4xy2z_Py+ABY*I_NAI_I4x2z_Py;
  Double I_NAI_I3x3y_D2y = I_NAI_K3x4y_Py+ABY*I_NAI_I3x3y_Py;
  Double I_NAI_I3x2yz_D2y = I_NAI_K3x3yz_Py+ABY*I_NAI_I3x2yz_Py;
  Double I_NAI_I3xy2z_D2y = I_NAI_K3x2y2z_Py+ABY*I_NAI_I3xy2z_Py;
  Double I_NAI_I3x3z_D2y = I_NAI_K3xy3z_Py+ABY*I_NAI_I3x3z_Py;
  Double I_NAI_I2x4y_D2y = I_NAI_K2x5y_Py+ABY*I_NAI_I2x4y_Py;
  Double I_NAI_I2x3yz_D2y = I_NAI_K2x4yz_Py+ABY*I_NAI_I2x3yz_Py;
  Double I_NAI_I2x2y2z_D2y = I_NAI_K2x3y2z_Py+ABY*I_NAI_I2x2y2z_Py;
  Double I_NAI_I2xy3z_D2y = I_NAI_K2x2y3z_Py+ABY*I_NAI_I2xy3z_Py;
  Double I_NAI_I2x4z_D2y = I_NAI_K2xy4z_Py+ABY*I_NAI_I2x4z_Py;
  Double I_NAI_Ix5y_D2y = I_NAI_Kx6y_Py+ABY*I_NAI_Ix5y_Py;
  Double I_NAI_Ix4yz_D2y = I_NAI_Kx5yz_Py+ABY*I_NAI_Ix4yz_Py;
  Double I_NAI_Ix3y2z_D2y = I_NAI_Kx4y2z_Py+ABY*I_NAI_Ix3y2z_Py;
  Double I_NAI_Ix2y3z_D2y = I_NAI_Kx3y3z_Py+ABY*I_NAI_Ix2y3z_Py;
  Double I_NAI_Ixy4z_D2y = I_NAI_Kx2y4z_Py+ABY*I_NAI_Ixy4z_Py;
  Double I_NAI_Ix5z_D2y = I_NAI_Kxy5z_Py+ABY*I_NAI_Ix5z_Py;
  Double I_NAI_I6y_D2y = I_NAI_K7y_Py+ABY*I_NAI_I6y_Py;
  Double I_NAI_I5yz_D2y = I_NAI_K6yz_Py+ABY*I_NAI_I5yz_Py;
  Double I_NAI_I4y2z_D2y = I_NAI_K5y2z_Py+ABY*I_NAI_I4y2z_Py;
  Double I_NAI_I3y3z_D2y = I_NAI_K4y3z_Py+ABY*I_NAI_I3y3z_Py;
  Double I_NAI_I2y4z_D2y = I_NAI_K3y4z_Py+ABY*I_NAI_I2y4z_Py;
  Double I_NAI_Iy5z_D2y = I_NAI_K2y5z_Py+ABY*I_NAI_Iy5z_Py;
  Double I_NAI_I6z_D2y = I_NAI_Ky6z_Py+ABY*I_NAI_I6z_Py;
  Double I_NAI_I6x_D2z = I_NAI_K6xz_Pz+ABZ*I_NAI_I6x_Pz;
  Double I_NAI_I5xy_D2z = I_NAI_K5xyz_Pz+ABZ*I_NAI_I5xy_Pz;
  Double I_NAI_I5xz_D2z = I_NAI_K5x2z_Pz+ABZ*I_NAI_I5xz_Pz;
  Double I_NAI_I4x2y_D2z = I_NAI_K4x2yz_Pz+ABZ*I_NAI_I4x2y_Pz;
  Double I_NAI_I4xyz_D2z = I_NAI_K4xy2z_Pz+ABZ*I_NAI_I4xyz_Pz;
  Double I_NAI_I4x2z_D2z = I_NAI_K4x3z_Pz+ABZ*I_NAI_I4x2z_Pz;
  Double I_NAI_I3x3y_D2z = I_NAI_K3x3yz_Pz+ABZ*I_NAI_I3x3y_Pz;
  Double I_NAI_I3x2yz_D2z = I_NAI_K3x2y2z_Pz+ABZ*I_NAI_I3x2yz_Pz;
  Double I_NAI_I3xy2z_D2z = I_NAI_K3xy3z_Pz+ABZ*I_NAI_I3xy2z_Pz;
  Double I_NAI_I3x3z_D2z = I_NAI_K3x4z_Pz+ABZ*I_NAI_I3x3z_Pz;
  Double I_NAI_I2x4y_D2z = I_NAI_K2x4yz_Pz+ABZ*I_NAI_I2x4y_Pz;
  Double I_NAI_I2x3yz_D2z = I_NAI_K2x3y2z_Pz+ABZ*I_NAI_I2x3yz_Pz;
  Double I_NAI_I2x2y2z_D2z = I_NAI_K2x2y3z_Pz+ABZ*I_NAI_I2x2y2z_Pz;
  Double I_NAI_I2xy3z_D2z = I_NAI_K2xy4z_Pz+ABZ*I_NAI_I2xy3z_Pz;
  Double I_NAI_I2x4z_D2z = I_NAI_K2x5z_Pz+ABZ*I_NAI_I2x4z_Pz;
  Double I_NAI_Ix5y_D2z = I_NAI_Kx5yz_Pz+ABZ*I_NAI_Ix5y_Pz;
  Double I_NAI_Ix4yz_D2z = I_NAI_Kx4y2z_Pz+ABZ*I_NAI_Ix4yz_Pz;
  Double I_NAI_Ix3y2z_D2z = I_NAI_Kx3y3z_Pz+ABZ*I_NAI_Ix3y2z_Pz;
  Double I_NAI_Ix2y3z_D2z = I_NAI_Kx2y4z_Pz+ABZ*I_NAI_Ix2y3z_Pz;
  Double I_NAI_Ixy4z_D2z = I_NAI_Kxy5z_Pz+ABZ*I_NAI_Ixy4z_Pz;
  Double I_NAI_Ix5z_D2z = I_NAI_Kx6z_Pz+ABZ*I_NAI_Ix5z_Pz;
  Double I_NAI_I6y_D2z = I_NAI_K6yz_Pz+ABZ*I_NAI_I6y_Pz;
  Double I_NAI_I5yz_D2z = I_NAI_K5y2z_Pz+ABZ*I_NAI_I5yz_Pz;
  Double I_NAI_I4y2z_D2z = I_NAI_K4y3z_Pz+ABZ*I_NAI_I4y2z_Pz;
  Double I_NAI_I3y3z_D2z = I_NAI_K3y4z_Pz+ABZ*I_NAI_I3y3z_Pz;
  Double I_NAI_I2y4z_D2z = I_NAI_K2y5z_Pz+ABZ*I_NAI_I2y4z_Pz;
  Double I_NAI_Iy5z_D2z = I_NAI_Ky6z_Pz+ABZ*I_NAI_Iy5z_Pz;
  Double I_NAI_I6z_D2z = I_NAI_K7z_Pz+ABZ*I_NAI_I6z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_H_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_D
   * RHS shell quartet name: SQ_NAI_H_D
   ************************************************************/
  Double I_NAI_H5x_F3x = I_NAI_I6x_D2x+ABX*I_NAI_H5x_D2x;
  Double I_NAI_H4xy_F3x = I_NAI_I5xy_D2x+ABX*I_NAI_H4xy_D2x;
  Double I_NAI_H4xz_F3x = I_NAI_I5xz_D2x+ABX*I_NAI_H4xz_D2x;
  Double I_NAI_H3x2y_F3x = I_NAI_I4x2y_D2x+ABX*I_NAI_H3x2y_D2x;
  Double I_NAI_H3xyz_F3x = I_NAI_I4xyz_D2x+ABX*I_NAI_H3xyz_D2x;
  Double I_NAI_H3x2z_F3x = I_NAI_I4x2z_D2x+ABX*I_NAI_H3x2z_D2x;
  Double I_NAI_H2x3y_F3x = I_NAI_I3x3y_D2x+ABX*I_NAI_H2x3y_D2x;
  Double I_NAI_H2x2yz_F3x = I_NAI_I3x2yz_D2x+ABX*I_NAI_H2x2yz_D2x;
  Double I_NAI_H2xy2z_F3x = I_NAI_I3xy2z_D2x+ABX*I_NAI_H2xy2z_D2x;
  Double I_NAI_H2x3z_F3x = I_NAI_I3x3z_D2x+ABX*I_NAI_H2x3z_D2x;
  Double I_NAI_Hx4y_F3x = I_NAI_I2x4y_D2x+ABX*I_NAI_Hx4y_D2x;
  Double I_NAI_Hx3yz_F3x = I_NAI_I2x3yz_D2x+ABX*I_NAI_Hx3yz_D2x;
  Double I_NAI_Hx2y2z_F3x = I_NAI_I2x2y2z_D2x+ABX*I_NAI_Hx2y2z_D2x;
  Double I_NAI_Hxy3z_F3x = I_NAI_I2xy3z_D2x+ABX*I_NAI_Hxy3z_D2x;
  Double I_NAI_Hx4z_F3x = I_NAI_I2x4z_D2x+ABX*I_NAI_Hx4z_D2x;
  Double I_NAI_H5y_F3x = I_NAI_Ix5y_D2x+ABX*I_NAI_H5y_D2x;
  Double I_NAI_H4yz_F3x = I_NAI_Ix4yz_D2x+ABX*I_NAI_H4yz_D2x;
  Double I_NAI_H3y2z_F3x = I_NAI_Ix3y2z_D2x+ABX*I_NAI_H3y2z_D2x;
  Double I_NAI_H2y3z_F3x = I_NAI_Ix2y3z_D2x+ABX*I_NAI_H2y3z_D2x;
  Double I_NAI_Hy4z_F3x = I_NAI_Ixy4z_D2x+ABX*I_NAI_Hy4z_D2x;
  Double I_NAI_H5z_F3x = I_NAI_Ix5z_D2x+ABX*I_NAI_H5z_D2x;
  Double I_NAI_H5x_F2xy = I_NAI_I5xy_D2x+ABY*I_NAI_H5x_D2x;
  Double I_NAI_H4xy_F2xy = I_NAI_I4x2y_D2x+ABY*I_NAI_H4xy_D2x;
  Double I_NAI_H4xz_F2xy = I_NAI_I4xyz_D2x+ABY*I_NAI_H4xz_D2x;
  Double I_NAI_H3x2y_F2xy = I_NAI_I3x3y_D2x+ABY*I_NAI_H3x2y_D2x;
  Double I_NAI_H3xyz_F2xy = I_NAI_I3x2yz_D2x+ABY*I_NAI_H3xyz_D2x;
  Double I_NAI_H3x2z_F2xy = I_NAI_I3xy2z_D2x+ABY*I_NAI_H3x2z_D2x;
  Double I_NAI_H2x3y_F2xy = I_NAI_I2x4y_D2x+ABY*I_NAI_H2x3y_D2x;
  Double I_NAI_H2x2yz_F2xy = I_NAI_I2x3yz_D2x+ABY*I_NAI_H2x2yz_D2x;
  Double I_NAI_H2xy2z_F2xy = I_NAI_I2x2y2z_D2x+ABY*I_NAI_H2xy2z_D2x;
  Double I_NAI_H2x3z_F2xy = I_NAI_I2xy3z_D2x+ABY*I_NAI_H2x3z_D2x;
  Double I_NAI_Hx4y_F2xy = I_NAI_Ix5y_D2x+ABY*I_NAI_Hx4y_D2x;
  Double I_NAI_Hx3yz_F2xy = I_NAI_Ix4yz_D2x+ABY*I_NAI_Hx3yz_D2x;
  Double I_NAI_Hx2y2z_F2xy = I_NAI_Ix3y2z_D2x+ABY*I_NAI_Hx2y2z_D2x;
  Double I_NAI_Hxy3z_F2xy = I_NAI_Ix2y3z_D2x+ABY*I_NAI_Hxy3z_D2x;
  Double I_NAI_Hx4z_F2xy = I_NAI_Ixy4z_D2x+ABY*I_NAI_Hx4z_D2x;
  Double I_NAI_H5y_F2xy = I_NAI_I6y_D2x+ABY*I_NAI_H5y_D2x;
  Double I_NAI_H4yz_F2xy = I_NAI_I5yz_D2x+ABY*I_NAI_H4yz_D2x;
  Double I_NAI_H3y2z_F2xy = I_NAI_I4y2z_D2x+ABY*I_NAI_H3y2z_D2x;
  Double I_NAI_H2y3z_F2xy = I_NAI_I3y3z_D2x+ABY*I_NAI_H2y3z_D2x;
  Double I_NAI_Hy4z_F2xy = I_NAI_I2y4z_D2x+ABY*I_NAI_Hy4z_D2x;
  Double I_NAI_H5z_F2xy = I_NAI_Iy5z_D2x+ABY*I_NAI_H5z_D2x;
  Double I_NAI_H5x_F2xz = I_NAI_I5xz_D2x+ABZ*I_NAI_H5x_D2x;
  Double I_NAI_H4xy_F2xz = I_NAI_I4xyz_D2x+ABZ*I_NAI_H4xy_D2x;
  Double I_NAI_H4xz_F2xz = I_NAI_I4x2z_D2x+ABZ*I_NAI_H4xz_D2x;
  Double I_NAI_H3x2y_F2xz = I_NAI_I3x2yz_D2x+ABZ*I_NAI_H3x2y_D2x;
  Double I_NAI_H3xyz_F2xz = I_NAI_I3xy2z_D2x+ABZ*I_NAI_H3xyz_D2x;
  Double I_NAI_H3x2z_F2xz = I_NAI_I3x3z_D2x+ABZ*I_NAI_H3x2z_D2x;
  Double I_NAI_H2x3y_F2xz = I_NAI_I2x3yz_D2x+ABZ*I_NAI_H2x3y_D2x;
  Double I_NAI_H2x2yz_F2xz = I_NAI_I2x2y2z_D2x+ABZ*I_NAI_H2x2yz_D2x;
  Double I_NAI_H2xy2z_F2xz = I_NAI_I2xy3z_D2x+ABZ*I_NAI_H2xy2z_D2x;
  Double I_NAI_H2x3z_F2xz = I_NAI_I2x4z_D2x+ABZ*I_NAI_H2x3z_D2x;
  Double I_NAI_Hx4y_F2xz = I_NAI_Ix4yz_D2x+ABZ*I_NAI_Hx4y_D2x;
  Double I_NAI_Hx3yz_F2xz = I_NAI_Ix3y2z_D2x+ABZ*I_NAI_Hx3yz_D2x;
  Double I_NAI_Hx2y2z_F2xz = I_NAI_Ix2y3z_D2x+ABZ*I_NAI_Hx2y2z_D2x;
  Double I_NAI_Hxy3z_F2xz = I_NAI_Ixy4z_D2x+ABZ*I_NAI_Hxy3z_D2x;
  Double I_NAI_Hx4z_F2xz = I_NAI_Ix5z_D2x+ABZ*I_NAI_Hx4z_D2x;
  Double I_NAI_H5y_F2xz = I_NAI_I5yz_D2x+ABZ*I_NAI_H5y_D2x;
  Double I_NAI_H4yz_F2xz = I_NAI_I4y2z_D2x+ABZ*I_NAI_H4yz_D2x;
  Double I_NAI_H3y2z_F2xz = I_NAI_I3y3z_D2x+ABZ*I_NAI_H3y2z_D2x;
  Double I_NAI_H2y3z_F2xz = I_NAI_I2y4z_D2x+ABZ*I_NAI_H2y3z_D2x;
  Double I_NAI_Hy4z_F2xz = I_NAI_Iy5z_D2x+ABZ*I_NAI_Hy4z_D2x;
  Double I_NAI_H5z_F2xz = I_NAI_I6z_D2x+ABZ*I_NAI_H5z_D2x;
  Double I_NAI_H5x_Fx2y = I_NAI_I6x_D2y+ABX*I_NAI_H5x_D2y;
  Double I_NAI_H4xy_Fx2y = I_NAI_I5xy_D2y+ABX*I_NAI_H4xy_D2y;
  Double I_NAI_H4xz_Fx2y = I_NAI_I5xz_D2y+ABX*I_NAI_H4xz_D2y;
  Double I_NAI_H3x2y_Fx2y = I_NAI_I4x2y_D2y+ABX*I_NAI_H3x2y_D2y;
  Double I_NAI_H3xyz_Fx2y = I_NAI_I4xyz_D2y+ABX*I_NAI_H3xyz_D2y;
  Double I_NAI_H3x2z_Fx2y = I_NAI_I4x2z_D2y+ABX*I_NAI_H3x2z_D2y;
  Double I_NAI_H2x3y_Fx2y = I_NAI_I3x3y_D2y+ABX*I_NAI_H2x3y_D2y;
  Double I_NAI_H2x2yz_Fx2y = I_NAI_I3x2yz_D2y+ABX*I_NAI_H2x2yz_D2y;
  Double I_NAI_H2xy2z_Fx2y = I_NAI_I3xy2z_D2y+ABX*I_NAI_H2xy2z_D2y;
  Double I_NAI_H2x3z_Fx2y = I_NAI_I3x3z_D2y+ABX*I_NAI_H2x3z_D2y;
  Double I_NAI_Hx4y_Fx2y = I_NAI_I2x4y_D2y+ABX*I_NAI_Hx4y_D2y;
  Double I_NAI_Hx3yz_Fx2y = I_NAI_I2x3yz_D2y+ABX*I_NAI_Hx3yz_D2y;
  Double I_NAI_Hx2y2z_Fx2y = I_NAI_I2x2y2z_D2y+ABX*I_NAI_Hx2y2z_D2y;
  Double I_NAI_Hxy3z_Fx2y = I_NAI_I2xy3z_D2y+ABX*I_NAI_Hxy3z_D2y;
  Double I_NAI_Hx4z_Fx2y = I_NAI_I2x4z_D2y+ABX*I_NAI_Hx4z_D2y;
  Double I_NAI_H5y_Fx2y = I_NAI_Ix5y_D2y+ABX*I_NAI_H5y_D2y;
  Double I_NAI_H4yz_Fx2y = I_NAI_Ix4yz_D2y+ABX*I_NAI_H4yz_D2y;
  Double I_NAI_H3y2z_Fx2y = I_NAI_Ix3y2z_D2y+ABX*I_NAI_H3y2z_D2y;
  Double I_NAI_H2y3z_Fx2y = I_NAI_Ix2y3z_D2y+ABX*I_NAI_H2y3z_D2y;
  Double I_NAI_Hy4z_Fx2y = I_NAI_Ixy4z_D2y+ABX*I_NAI_Hy4z_D2y;
  Double I_NAI_H5z_Fx2y = I_NAI_Ix5z_D2y+ABX*I_NAI_H5z_D2y;
  Double I_NAI_H5x_Fxyz = I_NAI_I5xz_Dxy+ABZ*I_NAI_H5x_Dxy;
  Double I_NAI_H4xy_Fxyz = I_NAI_I4xyz_Dxy+ABZ*I_NAI_H4xy_Dxy;
  Double I_NAI_H4xz_Fxyz = I_NAI_I4x2z_Dxy+ABZ*I_NAI_H4xz_Dxy;
  Double I_NAI_H3x2y_Fxyz = I_NAI_I3x2yz_Dxy+ABZ*I_NAI_H3x2y_Dxy;
  Double I_NAI_H3xyz_Fxyz = I_NAI_I3xy2z_Dxy+ABZ*I_NAI_H3xyz_Dxy;
  Double I_NAI_H3x2z_Fxyz = I_NAI_I3x3z_Dxy+ABZ*I_NAI_H3x2z_Dxy;
  Double I_NAI_H2x3y_Fxyz = I_NAI_I2x3yz_Dxy+ABZ*I_NAI_H2x3y_Dxy;
  Double I_NAI_H2x2yz_Fxyz = I_NAI_I2x2y2z_Dxy+ABZ*I_NAI_H2x2yz_Dxy;
  Double I_NAI_H2xy2z_Fxyz = I_NAI_I2xy3z_Dxy+ABZ*I_NAI_H2xy2z_Dxy;
  Double I_NAI_H2x3z_Fxyz = I_NAI_I2x4z_Dxy+ABZ*I_NAI_H2x3z_Dxy;
  Double I_NAI_Hx4y_Fxyz = I_NAI_Ix4yz_Dxy+ABZ*I_NAI_Hx4y_Dxy;
  Double I_NAI_Hx3yz_Fxyz = I_NAI_Ix3y2z_Dxy+ABZ*I_NAI_Hx3yz_Dxy;
  Double I_NAI_Hx2y2z_Fxyz = I_NAI_Ix2y3z_Dxy+ABZ*I_NAI_Hx2y2z_Dxy;
  Double I_NAI_Hxy3z_Fxyz = I_NAI_Ixy4z_Dxy+ABZ*I_NAI_Hxy3z_Dxy;
  Double I_NAI_Hx4z_Fxyz = I_NAI_Ix5z_Dxy+ABZ*I_NAI_Hx4z_Dxy;
  Double I_NAI_H5y_Fxyz = I_NAI_I5yz_Dxy+ABZ*I_NAI_H5y_Dxy;
  Double I_NAI_H4yz_Fxyz = I_NAI_I4y2z_Dxy+ABZ*I_NAI_H4yz_Dxy;
  Double I_NAI_H3y2z_Fxyz = I_NAI_I3y3z_Dxy+ABZ*I_NAI_H3y2z_Dxy;
  Double I_NAI_H2y3z_Fxyz = I_NAI_I2y4z_Dxy+ABZ*I_NAI_H2y3z_Dxy;
  Double I_NAI_Hy4z_Fxyz = I_NAI_Iy5z_Dxy+ABZ*I_NAI_Hy4z_Dxy;
  Double I_NAI_H5z_Fxyz = I_NAI_I6z_Dxy+ABZ*I_NAI_H5z_Dxy;
  Double I_NAI_H5x_Fx2z = I_NAI_I6x_D2z+ABX*I_NAI_H5x_D2z;
  Double I_NAI_H4xy_Fx2z = I_NAI_I5xy_D2z+ABX*I_NAI_H4xy_D2z;
  Double I_NAI_H4xz_Fx2z = I_NAI_I5xz_D2z+ABX*I_NAI_H4xz_D2z;
  Double I_NAI_H3x2y_Fx2z = I_NAI_I4x2y_D2z+ABX*I_NAI_H3x2y_D2z;
  Double I_NAI_H3xyz_Fx2z = I_NAI_I4xyz_D2z+ABX*I_NAI_H3xyz_D2z;
  Double I_NAI_H3x2z_Fx2z = I_NAI_I4x2z_D2z+ABX*I_NAI_H3x2z_D2z;
  Double I_NAI_H2x3y_Fx2z = I_NAI_I3x3y_D2z+ABX*I_NAI_H2x3y_D2z;
  Double I_NAI_H2x2yz_Fx2z = I_NAI_I3x2yz_D2z+ABX*I_NAI_H2x2yz_D2z;
  Double I_NAI_H2xy2z_Fx2z = I_NAI_I3xy2z_D2z+ABX*I_NAI_H2xy2z_D2z;
  Double I_NAI_H2x3z_Fx2z = I_NAI_I3x3z_D2z+ABX*I_NAI_H2x3z_D2z;
  Double I_NAI_Hx4y_Fx2z = I_NAI_I2x4y_D2z+ABX*I_NAI_Hx4y_D2z;
  Double I_NAI_Hx3yz_Fx2z = I_NAI_I2x3yz_D2z+ABX*I_NAI_Hx3yz_D2z;
  Double I_NAI_Hx2y2z_Fx2z = I_NAI_I2x2y2z_D2z+ABX*I_NAI_Hx2y2z_D2z;
  Double I_NAI_Hxy3z_Fx2z = I_NAI_I2xy3z_D2z+ABX*I_NAI_Hxy3z_D2z;
  Double I_NAI_Hx4z_Fx2z = I_NAI_I2x4z_D2z+ABX*I_NAI_Hx4z_D2z;
  Double I_NAI_H5y_Fx2z = I_NAI_Ix5y_D2z+ABX*I_NAI_H5y_D2z;
  Double I_NAI_H4yz_Fx2z = I_NAI_Ix4yz_D2z+ABX*I_NAI_H4yz_D2z;
  Double I_NAI_H3y2z_Fx2z = I_NAI_Ix3y2z_D2z+ABX*I_NAI_H3y2z_D2z;
  Double I_NAI_H2y3z_Fx2z = I_NAI_Ix2y3z_D2z+ABX*I_NAI_H2y3z_D2z;
  Double I_NAI_Hy4z_Fx2z = I_NAI_Ixy4z_D2z+ABX*I_NAI_Hy4z_D2z;
  Double I_NAI_H5z_Fx2z = I_NAI_Ix5z_D2z+ABX*I_NAI_H5z_D2z;
  Double I_NAI_H5x_F3y = I_NAI_I5xy_D2y+ABY*I_NAI_H5x_D2y;
  Double I_NAI_H4xy_F3y = I_NAI_I4x2y_D2y+ABY*I_NAI_H4xy_D2y;
  Double I_NAI_H4xz_F3y = I_NAI_I4xyz_D2y+ABY*I_NAI_H4xz_D2y;
  Double I_NAI_H3x2y_F3y = I_NAI_I3x3y_D2y+ABY*I_NAI_H3x2y_D2y;
  Double I_NAI_H3xyz_F3y = I_NAI_I3x2yz_D2y+ABY*I_NAI_H3xyz_D2y;
  Double I_NAI_H3x2z_F3y = I_NAI_I3xy2z_D2y+ABY*I_NAI_H3x2z_D2y;
  Double I_NAI_H2x3y_F3y = I_NAI_I2x4y_D2y+ABY*I_NAI_H2x3y_D2y;
  Double I_NAI_H2x2yz_F3y = I_NAI_I2x3yz_D2y+ABY*I_NAI_H2x2yz_D2y;
  Double I_NAI_H2xy2z_F3y = I_NAI_I2x2y2z_D2y+ABY*I_NAI_H2xy2z_D2y;
  Double I_NAI_H2x3z_F3y = I_NAI_I2xy3z_D2y+ABY*I_NAI_H2x3z_D2y;
  Double I_NAI_Hx4y_F3y = I_NAI_Ix5y_D2y+ABY*I_NAI_Hx4y_D2y;
  Double I_NAI_Hx3yz_F3y = I_NAI_Ix4yz_D2y+ABY*I_NAI_Hx3yz_D2y;
  Double I_NAI_Hx2y2z_F3y = I_NAI_Ix3y2z_D2y+ABY*I_NAI_Hx2y2z_D2y;
  Double I_NAI_Hxy3z_F3y = I_NAI_Ix2y3z_D2y+ABY*I_NAI_Hxy3z_D2y;
  Double I_NAI_Hx4z_F3y = I_NAI_Ixy4z_D2y+ABY*I_NAI_Hx4z_D2y;
  Double I_NAI_H5y_F3y = I_NAI_I6y_D2y+ABY*I_NAI_H5y_D2y;
  Double I_NAI_H4yz_F3y = I_NAI_I5yz_D2y+ABY*I_NAI_H4yz_D2y;
  Double I_NAI_H3y2z_F3y = I_NAI_I4y2z_D2y+ABY*I_NAI_H3y2z_D2y;
  Double I_NAI_H2y3z_F3y = I_NAI_I3y3z_D2y+ABY*I_NAI_H2y3z_D2y;
  Double I_NAI_Hy4z_F3y = I_NAI_I2y4z_D2y+ABY*I_NAI_Hy4z_D2y;
  Double I_NAI_H5z_F3y = I_NAI_Iy5z_D2y+ABY*I_NAI_H5z_D2y;
  Double I_NAI_H5x_F2yz = I_NAI_I5xz_D2y+ABZ*I_NAI_H5x_D2y;
  Double I_NAI_H4xy_F2yz = I_NAI_I4xyz_D2y+ABZ*I_NAI_H4xy_D2y;
  Double I_NAI_H4xz_F2yz = I_NAI_I4x2z_D2y+ABZ*I_NAI_H4xz_D2y;
  Double I_NAI_H3x2y_F2yz = I_NAI_I3x2yz_D2y+ABZ*I_NAI_H3x2y_D2y;
  Double I_NAI_H3xyz_F2yz = I_NAI_I3xy2z_D2y+ABZ*I_NAI_H3xyz_D2y;
  Double I_NAI_H3x2z_F2yz = I_NAI_I3x3z_D2y+ABZ*I_NAI_H3x2z_D2y;
  Double I_NAI_H2x3y_F2yz = I_NAI_I2x3yz_D2y+ABZ*I_NAI_H2x3y_D2y;
  Double I_NAI_H2x2yz_F2yz = I_NAI_I2x2y2z_D2y+ABZ*I_NAI_H2x2yz_D2y;
  Double I_NAI_H2xy2z_F2yz = I_NAI_I2xy3z_D2y+ABZ*I_NAI_H2xy2z_D2y;
  Double I_NAI_H2x3z_F2yz = I_NAI_I2x4z_D2y+ABZ*I_NAI_H2x3z_D2y;
  Double I_NAI_Hx4y_F2yz = I_NAI_Ix4yz_D2y+ABZ*I_NAI_Hx4y_D2y;
  Double I_NAI_Hx3yz_F2yz = I_NAI_Ix3y2z_D2y+ABZ*I_NAI_Hx3yz_D2y;
  Double I_NAI_Hx2y2z_F2yz = I_NAI_Ix2y3z_D2y+ABZ*I_NAI_Hx2y2z_D2y;
  Double I_NAI_Hxy3z_F2yz = I_NAI_Ixy4z_D2y+ABZ*I_NAI_Hxy3z_D2y;
  Double I_NAI_Hx4z_F2yz = I_NAI_Ix5z_D2y+ABZ*I_NAI_Hx4z_D2y;
  Double I_NAI_H5y_F2yz = I_NAI_I5yz_D2y+ABZ*I_NAI_H5y_D2y;
  Double I_NAI_H4yz_F2yz = I_NAI_I4y2z_D2y+ABZ*I_NAI_H4yz_D2y;
  Double I_NAI_H3y2z_F2yz = I_NAI_I3y3z_D2y+ABZ*I_NAI_H3y2z_D2y;
  Double I_NAI_H2y3z_F2yz = I_NAI_I2y4z_D2y+ABZ*I_NAI_H2y3z_D2y;
  Double I_NAI_Hy4z_F2yz = I_NAI_Iy5z_D2y+ABZ*I_NAI_Hy4z_D2y;
  Double I_NAI_H5z_F2yz = I_NAI_I6z_D2y+ABZ*I_NAI_H5z_D2y;
  Double I_NAI_H5x_Fy2z = I_NAI_I5xy_D2z+ABY*I_NAI_H5x_D2z;
  Double I_NAI_H4xy_Fy2z = I_NAI_I4x2y_D2z+ABY*I_NAI_H4xy_D2z;
  Double I_NAI_H4xz_Fy2z = I_NAI_I4xyz_D2z+ABY*I_NAI_H4xz_D2z;
  Double I_NAI_H3x2y_Fy2z = I_NAI_I3x3y_D2z+ABY*I_NAI_H3x2y_D2z;
  Double I_NAI_H3xyz_Fy2z = I_NAI_I3x2yz_D2z+ABY*I_NAI_H3xyz_D2z;
  Double I_NAI_H3x2z_Fy2z = I_NAI_I3xy2z_D2z+ABY*I_NAI_H3x2z_D2z;
  Double I_NAI_H2x3y_Fy2z = I_NAI_I2x4y_D2z+ABY*I_NAI_H2x3y_D2z;
  Double I_NAI_H2x2yz_Fy2z = I_NAI_I2x3yz_D2z+ABY*I_NAI_H2x2yz_D2z;
  Double I_NAI_H2xy2z_Fy2z = I_NAI_I2x2y2z_D2z+ABY*I_NAI_H2xy2z_D2z;
  Double I_NAI_H2x3z_Fy2z = I_NAI_I2xy3z_D2z+ABY*I_NAI_H2x3z_D2z;
  Double I_NAI_Hx4y_Fy2z = I_NAI_Ix5y_D2z+ABY*I_NAI_Hx4y_D2z;
  Double I_NAI_Hx3yz_Fy2z = I_NAI_Ix4yz_D2z+ABY*I_NAI_Hx3yz_D2z;
  Double I_NAI_Hx2y2z_Fy2z = I_NAI_Ix3y2z_D2z+ABY*I_NAI_Hx2y2z_D2z;
  Double I_NAI_Hxy3z_Fy2z = I_NAI_Ix2y3z_D2z+ABY*I_NAI_Hxy3z_D2z;
  Double I_NAI_Hx4z_Fy2z = I_NAI_Ixy4z_D2z+ABY*I_NAI_Hx4z_D2z;
  Double I_NAI_H5y_Fy2z = I_NAI_I6y_D2z+ABY*I_NAI_H5y_D2z;
  Double I_NAI_H4yz_Fy2z = I_NAI_I5yz_D2z+ABY*I_NAI_H4yz_D2z;
  Double I_NAI_H3y2z_Fy2z = I_NAI_I4y2z_D2z+ABY*I_NAI_H3y2z_D2z;
  Double I_NAI_H2y3z_Fy2z = I_NAI_I3y3z_D2z+ABY*I_NAI_H2y3z_D2z;
  Double I_NAI_Hy4z_Fy2z = I_NAI_I2y4z_D2z+ABY*I_NAI_Hy4z_D2z;
  Double I_NAI_H5z_Fy2z = I_NAI_Iy5z_D2z+ABY*I_NAI_H5z_D2z;
  Double I_NAI_H5x_F3z = I_NAI_I5xz_D2z+ABZ*I_NAI_H5x_D2z;
  Double I_NAI_H4xy_F3z = I_NAI_I4xyz_D2z+ABZ*I_NAI_H4xy_D2z;
  Double I_NAI_H4xz_F3z = I_NAI_I4x2z_D2z+ABZ*I_NAI_H4xz_D2z;
  Double I_NAI_H3x2y_F3z = I_NAI_I3x2yz_D2z+ABZ*I_NAI_H3x2y_D2z;
  Double I_NAI_H3xyz_F3z = I_NAI_I3xy2z_D2z+ABZ*I_NAI_H3xyz_D2z;
  Double I_NAI_H3x2z_F3z = I_NAI_I3x3z_D2z+ABZ*I_NAI_H3x2z_D2z;
  Double I_NAI_H2x3y_F3z = I_NAI_I2x3yz_D2z+ABZ*I_NAI_H2x3y_D2z;
  Double I_NAI_H2x2yz_F3z = I_NAI_I2x2y2z_D2z+ABZ*I_NAI_H2x2yz_D2z;
  Double I_NAI_H2xy2z_F3z = I_NAI_I2xy3z_D2z+ABZ*I_NAI_H2xy2z_D2z;
  Double I_NAI_H2x3z_F3z = I_NAI_I2x4z_D2z+ABZ*I_NAI_H2x3z_D2z;
  Double I_NAI_Hx4y_F3z = I_NAI_Ix4yz_D2z+ABZ*I_NAI_Hx4y_D2z;
  Double I_NAI_Hx3yz_F3z = I_NAI_Ix3y2z_D2z+ABZ*I_NAI_Hx3yz_D2z;
  Double I_NAI_Hx2y2z_F3z = I_NAI_Ix2y3z_D2z+ABZ*I_NAI_Hx2y2z_D2z;
  Double I_NAI_Hxy3z_F3z = I_NAI_Ixy4z_D2z+ABZ*I_NAI_Hxy3z_D2z;
  Double I_NAI_Hx4z_F3z = I_NAI_Ix5z_D2z+ABZ*I_NAI_Hx4z_D2z;
  Double I_NAI_H5y_F3z = I_NAI_I5yz_D2z+ABZ*I_NAI_H5y_D2z;
  Double I_NAI_H4yz_F3z = I_NAI_I4y2z_D2z+ABZ*I_NAI_H4yz_D2z;
  Double I_NAI_H3y2z_F3z = I_NAI_I3y3z_D2z+ABZ*I_NAI_H3y2z_D2z;
  Double I_NAI_H2y3z_F3z = I_NAI_I2y4z_D2z+ABZ*I_NAI_H2y3z_D2z;
  Double I_NAI_Hy4z_F3z = I_NAI_Iy5z_D2z+ABZ*I_NAI_Hy4z_D2z;
  Double I_NAI_H5z_F3z = I_NAI_I6z_D2z+ABZ*I_NAI_H5z_D2z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_G
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_F
   * RHS shell quartet name: SQ_NAI_G_F
   ************************************************************/
  Double I_NAI_G4x_G4x = I_NAI_H5x_F3x+ABX*I_NAI_G4x_F3x;
  Double I_NAI_G3xy_G4x = I_NAI_H4xy_F3x+ABX*I_NAI_G3xy_F3x;
  Double I_NAI_G3xz_G4x = I_NAI_H4xz_F3x+ABX*I_NAI_G3xz_F3x;
  Double I_NAI_G2x2y_G4x = I_NAI_H3x2y_F3x+ABX*I_NAI_G2x2y_F3x;
  Double I_NAI_G2xyz_G4x = I_NAI_H3xyz_F3x+ABX*I_NAI_G2xyz_F3x;
  Double I_NAI_G2x2z_G4x = I_NAI_H3x2z_F3x+ABX*I_NAI_G2x2z_F3x;
  Double I_NAI_Gx3y_G4x = I_NAI_H2x3y_F3x+ABX*I_NAI_Gx3y_F3x;
  Double I_NAI_Gx2yz_G4x = I_NAI_H2x2yz_F3x+ABX*I_NAI_Gx2yz_F3x;
  Double I_NAI_Gxy2z_G4x = I_NAI_H2xy2z_F3x+ABX*I_NAI_Gxy2z_F3x;
  Double I_NAI_Gx3z_G4x = I_NAI_H2x3z_F3x+ABX*I_NAI_Gx3z_F3x;
  Double I_NAI_G4y_G4x = I_NAI_Hx4y_F3x+ABX*I_NAI_G4y_F3x;
  Double I_NAI_G3yz_G4x = I_NAI_Hx3yz_F3x+ABX*I_NAI_G3yz_F3x;
  Double I_NAI_G2y2z_G4x = I_NAI_Hx2y2z_F3x+ABX*I_NAI_G2y2z_F3x;
  Double I_NAI_Gy3z_G4x = I_NAI_Hxy3z_F3x+ABX*I_NAI_Gy3z_F3x;
  Double I_NAI_G4z_G4x = I_NAI_Hx4z_F3x+ABX*I_NAI_G4z_F3x;
  Double I_NAI_G4x_G3xy = I_NAI_H4xy_F3x+ABY*I_NAI_G4x_F3x;
  Double I_NAI_G3xy_G3xy = I_NAI_H3x2y_F3x+ABY*I_NAI_G3xy_F3x;
  Double I_NAI_G3xz_G3xy = I_NAI_H3xyz_F3x+ABY*I_NAI_G3xz_F3x;
  Double I_NAI_G2x2y_G3xy = I_NAI_H2x3y_F3x+ABY*I_NAI_G2x2y_F3x;
  Double I_NAI_G2xyz_G3xy = I_NAI_H2x2yz_F3x+ABY*I_NAI_G2xyz_F3x;
  Double I_NAI_G2x2z_G3xy = I_NAI_H2xy2z_F3x+ABY*I_NAI_G2x2z_F3x;
  Double I_NAI_Gx3y_G3xy = I_NAI_Hx4y_F3x+ABY*I_NAI_Gx3y_F3x;
  Double I_NAI_Gx2yz_G3xy = I_NAI_Hx3yz_F3x+ABY*I_NAI_Gx2yz_F3x;
  Double I_NAI_Gxy2z_G3xy = I_NAI_Hx2y2z_F3x+ABY*I_NAI_Gxy2z_F3x;
  Double I_NAI_Gx3z_G3xy = I_NAI_Hxy3z_F3x+ABY*I_NAI_Gx3z_F3x;
  Double I_NAI_G4y_G3xy = I_NAI_H5y_F3x+ABY*I_NAI_G4y_F3x;
  Double I_NAI_G3yz_G3xy = I_NAI_H4yz_F3x+ABY*I_NAI_G3yz_F3x;
  Double I_NAI_G2y2z_G3xy = I_NAI_H3y2z_F3x+ABY*I_NAI_G2y2z_F3x;
  Double I_NAI_Gy3z_G3xy = I_NAI_H2y3z_F3x+ABY*I_NAI_Gy3z_F3x;
  Double I_NAI_G4z_G3xy = I_NAI_Hy4z_F3x+ABY*I_NAI_G4z_F3x;
  Double I_NAI_G4x_G3xz = I_NAI_H4xz_F3x+ABZ*I_NAI_G4x_F3x;
  Double I_NAI_G3xy_G3xz = I_NAI_H3xyz_F3x+ABZ*I_NAI_G3xy_F3x;
  Double I_NAI_G3xz_G3xz = I_NAI_H3x2z_F3x+ABZ*I_NAI_G3xz_F3x;
  Double I_NAI_G2x2y_G3xz = I_NAI_H2x2yz_F3x+ABZ*I_NAI_G2x2y_F3x;
  Double I_NAI_G2xyz_G3xz = I_NAI_H2xy2z_F3x+ABZ*I_NAI_G2xyz_F3x;
  Double I_NAI_G2x2z_G3xz = I_NAI_H2x3z_F3x+ABZ*I_NAI_G2x2z_F3x;
  Double I_NAI_Gx3y_G3xz = I_NAI_Hx3yz_F3x+ABZ*I_NAI_Gx3y_F3x;
  Double I_NAI_Gx2yz_G3xz = I_NAI_Hx2y2z_F3x+ABZ*I_NAI_Gx2yz_F3x;
  Double I_NAI_Gxy2z_G3xz = I_NAI_Hxy3z_F3x+ABZ*I_NAI_Gxy2z_F3x;
  Double I_NAI_Gx3z_G3xz = I_NAI_Hx4z_F3x+ABZ*I_NAI_Gx3z_F3x;
  Double I_NAI_G4y_G3xz = I_NAI_H4yz_F3x+ABZ*I_NAI_G4y_F3x;
  Double I_NAI_G3yz_G3xz = I_NAI_H3y2z_F3x+ABZ*I_NAI_G3yz_F3x;
  Double I_NAI_G2y2z_G3xz = I_NAI_H2y3z_F3x+ABZ*I_NAI_G2y2z_F3x;
  Double I_NAI_Gy3z_G3xz = I_NAI_Hy4z_F3x+ABZ*I_NAI_Gy3z_F3x;
  Double I_NAI_G4z_G3xz = I_NAI_H5z_F3x+ABZ*I_NAI_G4z_F3x;
  Double I_NAI_G4x_G2x2y = I_NAI_H4xy_F2xy+ABY*I_NAI_G4x_F2xy;
  Double I_NAI_G3xy_G2x2y = I_NAI_H3x2y_F2xy+ABY*I_NAI_G3xy_F2xy;
  Double I_NAI_G3xz_G2x2y = I_NAI_H3xyz_F2xy+ABY*I_NAI_G3xz_F2xy;
  Double I_NAI_G2x2y_G2x2y = I_NAI_H2x3y_F2xy+ABY*I_NAI_G2x2y_F2xy;
  Double I_NAI_G2xyz_G2x2y = I_NAI_H2x2yz_F2xy+ABY*I_NAI_G2xyz_F2xy;
  Double I_NAI_G2x2z_G2x2y = I_NAI_H2xy2z_F2xy+ABY*I_NAI_G2x2z_F2xy;
  Double I_NAI_Gx3y_G2x2y = I_NAI_Hx4y_F2xy+ABY*I_NAI_Gx3y_F2xy;
  Double I_NAI_Gx2yz_G2x2y = I_NAI_Hx3yz_F2xy+ABY*I_NAI_Gx2yz_F2xy;
  Double I_NAI_Gxy2z_G2x2y = I_NAI_Hx2y2z_F2xy+ABY*I_NAI_Gxy2z_F2xy;
  Double I_NAI_Gx3z_G2x2y = I_NAI_Hxy3z_F2xy+ABY*I_NAI_Gx3z_F2xy;
  Double I_NAI_G4y_G2x2y = I_NAI_H5y_F2xy+ABY*I_NAI_G4y_F2xy;
  Double I_NAI_G3yz_G2x2y = I_NAI_H4yz_F2xy+ABY*I_NAI_G3yz_F2xy;
  Double I_NAI_G2y2z_G2x2y = I_NAI_H3y2z_F2xy+ABY*I_NAI_G2y2z_F2xy;
  Double I_NAI_Gy3z_G2x2y = I_NAI_H2y3z_F2xy+ABY*I_NAI_Gy3z_F2xy;
  Double I_NAI_G4z_G2x2y = I_NAI_Hy4z_F2xy+ABY*I_NAI_G4z_F2xy;
  Double I_NAI_G4x_G2xyz = I_NAI_H4xz_F2xy+ABZ*I_NAI_G4x_F2xy;
  Double I_NAI_G3xy_G2xyz = I_NAI_H3xyz_F2xy+ABZ*I_NAI_G3xy_F2xy;
  Double I_NAI_G3xz_G2xyz = I_NAI_H3x2z_F2xy+ABZ*I_NAI_G3xz_F2xy;
  Double I_NAI_G2x2y_G2xyz = I_NAI_H2x2yz_F2xy+ABZ*I_NAI_G2x2y_F2xy;
  Double I_NAI_G2xyz_G2xyz = I_NAI_H2xy2z_F2xy+ABZ*I_NAI_G2xyz_F2xy;
  Double I_NAI_G2x2z_G2xyz = I_NAI_H2x3z_F2xy+ABZ*I_NAI_G2x2z_F2xy;
  Double I_NAI_Gx3y_G2xyz = I_NAI_Hx3yz_F2xy+ABZ*I_NAI_Gx3y_F2xy;
  Double I_NAI_Gx2yz_G2xyz = I_NAI_Hx2y2z_F2xy+ABZ*I_NAI_Gx2yz_F2xy;
  Double I_NAI_Gxy2z_G2xyz = I_NAI_Hxy3z_F2xy+ABZ*I_NAI_Gxy2z_F2xy;
  Double I_NAI_Gx3z_G2xyz = I_NAI_Hx4z_F2xy+ABZ*I_NAI_Gx3z_F2xy;
  Double I_NAI_G4y_G2xyz = I_NAI_H4yz_F2xy+ABZ*I_NAI_G4y_F2xy;
  Double I_NAI_G3yz_G2xyz = I_NAI_H3y2z_F2xy+ABZ*I_NAI_G3yz_F2xy;
  Double I_NAI_G2y2z_G2xyz = I_NAI_H2y3z_F2xy+ABZ*I_NAI_G2y2z_F2xy;
  Double I_NAI_Gy3z_G2xyz = I_NAI_Hy4z_F2xy+ABZ*I_NAI_Gy3z_F2xy;
  Double I_NAI_G4z_G2xyz = I_NAI_H5z_F2xy+ABZ*I_NAI_G4z_F2xy;
  Double I_NAI_G4x_G2x2z = I_NAI_H4xz_F2xz+ABZ*I_NAI_G4x_F2xz;
  Double I_NAI_G3xy_G2x2z = I_NAI_H3xyz_F2xz+ABZ*I_NAI_G3xy_F2xz;
  Double I_NAI_G3xz_G2x2z = I_NAI_H3x2z_F2xz+ABZ*I_NAI_G3xz_F2xz;
  Double I_NAI_G2x2y_G2x2z = I_NAI_H2x2yz_F2xz+ABZ*I_NAI_G2x2y_F2xz;
  Double I_NAI_G2xyz_G2x2z = I_NAI_H2xy2z_F2xz+ABZ*I_NAI_G2xyz_F2xz;
  Double I_NAI_G2x2z_G2x2z = I_NAI_H2x3z_F2xz+ABZ*I_NAI_G2x2z_F2xz;
  Double I_NAI_Gx3y_G2x2z = I_NAI_Hx3yz_F2xz+ABZ*I_NAI_Gx3y_F2xz;
  Double I_NAI_Gx2yz_G2x2z = I_NAI_Hx2y2z_F2xz+ABZ*I_NAI_Gx2yz_F2xz;
  Double I_NAI_Gxy2z_G2x2z = I_NAI_Hxy3z_F2xz+ABZ*I_NAI_Gxy2z_F2xz;
  Double I_NAI_Gx3z_G2x2z = I_NAI_Hx4z_F2xz+ABZ*I_NAI_Gx3z_F2xz;
  Double I_NAI_G4y_G2x2z = I_NAI_H4yz_F2xz+ABZ*I_NAI_G4y_F2xz;
  Double I_NAI_G3yz_G2x2z = I_NAI_H3y2z_F2xz+ABZ*I_NAI_G3yz_F2xz;
  Double I_NAI_G2y2z_G2x2z = I_NAI_H2y3z_F2xz+ABZ*I_NAI_G2y2z_F2xz;
  Double I_NAI_Gy3z_G2x2z = I_NAI_Hy4z_F2xz+ABZ*I_NAI_Gy3z_F2xz;
  Double I_NAI_G4z_G2x2z = I_NAI_H5z_F2xz+ABZ*I_NAI_G4z_F2xz;
  Double I_NAI_G4x_Gx3y = I_NAI_H5x_F3y+ABX*I_NAI_G4x_F3y;
  Double I_NAI_G3xy_Gx3y = I_NAI_H4xy_F3y+ABX*I_NAI_G3xy_F3y;
  Double I_NAI_G3xz_Gx3y = I_NAI_H4xz_F3y+ABX*I_NAI_G3xz_F3y;
  Double I_NAI_G2x2y_Gx3y = I_NAI_H3x2y_F3y+ABX*I_NAI_G2x2y_F3y;
  Double I_NAI_G2xyz_Gx3y = I_NAI_H3xyz_F3y+ABX*I_NAI_G2xyz_F3y;
  Double I_NAI_G2x2z_Gx3y = I_NAI_H3x2z_F3y+ABX*I_NAI_G2x2z_F3y;
  Double I_NAI_Gx3y_Gx3y = I_NAI_H2x3y_F3y+ABX*I_NAI_Gx3y_F3y;
  Double I_NAI_Gx2yz_Gx3y = I_NAI_H2x2yz_F3y+ABX*I_NAI_Gx2yz_F3y;
  Double I_NAI_Gxy2z_Gx3y = I_NAI_H2xy2z_F3y+ABX*I_NAI_Gxy2z_F3y;
  Double I_NAI_Gx3z_Gx3y = I_NAI_H2x3z_F3y+ABX*I_NAI_Gx3z_F3y;
  Double I_NAI_G4y_Gx3y = I_NAI_Hx4y_F3y+ABX*I_NAI_G4y_F3y;
  Double I_NAI_G3yz_Gx3y = I_NAI_Hx3yz_F3y+ABX*I_NAI_G3yz_F3y;
  Double I_NAI_G2y2z_Gx3y = I_NAI_Hx2y2z_F3y+ABX*I_NAI_G2y2z_F3y;
  Double I_NAI_Gy3z_Gx3y = I_NAI_Hxy3z_F3y+ABX*I_NAI_Gy3z_F3y;
  Double I_NAI_G4z_Gx3y = I_NAI_Hx4z_F3y+ABX*I_NAI_G4z_F3y;
  Double I_NAI_G4x_Gx2yz = I_NAI_H4xz_Fx2y+ABZ*I_NAI_G4x_Fx2y;
  Double I_NAI_G3xy_Gx2yz = I_NAI_H3xyz_Fx2y+ABZ*I_NAI_G3xy_Fx2y;
  Double I_NAI_G3xz_Gx2yz = I_NAI_H3x2z_Fx2y+ABZ*I_NAI_G3xz_Fx2y;
  Double I_NAI_G2x2y_Gx2yz = I_NAI_H2x2yz_Fx2y+ABZ*I_NAI_G2x2y_Fx2y;
  Double I_NAI_G2xyz_Gx2yz = I_NAI_H2xy2z_Fx2y+ABZ*I_NAI_G2xyz_Fx2y;
  Double I_NAI_G2x2z_Gx2yz = I_NAI_H2x3z_Fx2y+ABZ*I_NAI_G2x2z_Fx2y;
  Double I_NAI_Gx3y_Gx2yz = I_NAI_Hx3yz_Fx2y+ABZ*I_NAI_Gx3y_Fx2y;
  Double I_NAI_Gx2yz_Gx2yz = I_NAI_Hx2y2z_Fx2y+ABZ*I_NAI_Gx2yz_Fx2y;
  Double I_NAI_Gxy2z_Gx2yz = I_NAI_Hxy3z_Fx2y+ABZ*I_NAI_Gxy2z_Fx2y;
  Double I_NAI_Gx3z_Gx2yz = I_NAI_Hx4z_Fx2y+ABZ*I_NAI_Gx3z_Fx2y;
  Double I_NAI_G4y_Gx2yz = I_NAI_H4yz_Fx2y+ABZ*I_NAI_G4y_Fx2y;
  Double I_NAI_G3yz_Gx2yz = I_NAI_H3y2z_Fx2y+ABZ*I_NAI_G3yz_Fx2y;
  Double I_NAI_G2y2z_Gx2yz = I_NAI_H2y3z_Fx2y+ABZ*I_NAI_G2y2z_Fx2y;
  Double I_NAI_Gy3z_Gx2yz = I_NAI_Hy4z_Fx2y+ABZ*I_NAI_Gy3z_Fx2y;
  Double I_NAI_G4z_Gx2yz = I_NAI_H5z_Fx2y+ABZ*I_NAI_G4z_Fx2y;
  Double I_NAI_G4x_Gxy2z = I_NAI_H4xy_Fx2z+ABY*I_NAI_G4x_Fx2z;
  Double I_NAI_G3xy_Gxy2z = I_NAI_H3x2y_Fx2z+ABY*I_NAI_G3xy_Fx2z;
  Double I_NAI_G3xz_Gxy2z = I_NAI_H3xyz_Fx2z+ABY*I_NAI_G3xz_Fx2z;
  Double I_NAI_G2x2y_Gxy2z = I_NAI_H2x3y_Fx2z+ABY*I_NAI_G2x2y_Fx2z;
  Double I_NAI_G2xyz_Gxy2z = I_NAI_H2x2yz_Fx2z+ABY*I_NAI_G2xyz_Fx2z;
  Double I_NAI_G2x2z_Gxy2z = I_NAI_H2xy2z_Fx2z+ABY*I_NAI_G2x2z_Fx2z;
  Double I_NAI_Gx3y_Gxy2z = I_NAI_Hx4y_Fx2z+ABY*I_NAI_Gx3y_Fx2z;
  Double I_NAI_Gx2yz_Gxy2z = I_NAI_Hx3yz_Fx2z+ABY*I_NAI_Gx2yz_Fx2z;
  Double I_NAI_Gxy2z_Gxy2z = I_NAI_Hx2y2z_Fx2z+ABY*I_NAI_Gxy2z_Fx2z;
  Double I_NAI_Gx3z_Gxy2z = I_NAI_Hxy3z_Fx2z+ABY*I_NAI_Gx3z_Fx2z;
  Double I_NAI_G4y_Gxy2z = I_NAI_H5y_Fx2z+ABY*I_NAI_G4y_Fx2z;
  Double I_NAI_G3yz_Gxy2z = I_NAI_H4yz_Fx2z+ABY*I_NAI_G3yz_Fx2z;
  Double I_NAI_G2y2z_Gxy2z = I_NAI_H3y2z_Fx2z+ABY*I_NAI_G2y2z_Fx2z;
  Double I_NAI_Gy3z_Gxy2z = I_NAI_H2y3z_Fx2z+ABY*I_NAI_Gy3z_Fx2z;
  Double I_NAI_G4z_Gxy2z = I_NAI_Hy4z_Fx2z+ABY*I_NAI_G4z_Fx2z;
  Double I_NAI_G4x_Gx3z = I_NAI_H5x_F3z+ABX*I_NAI_G4x_F3z;
  Double I_NAI_G3xy_Gx3z = I_NAI_H4xy_F3z+ABX*I_NAI_G3xy_F3z;
  Double I_NAI_G3xz_Gx3z = I_NAI_H4xz_F3z+ABX*I_NAI_G3xz_F3z;
  Double I_NAI_G2x2y_Gx3z = I_NAI_H3x2y_F3z+ABX*I_NAI_G2x2y_F3z;
  Double I_NAI_G2xyz_Gx3z = I_NAI_H3xyz_F3z+ABX*I_NAI_G2xyz_F3z;
  Double I_NAI_G2x2z_Gx3z = I_NAI_H3x2z_F3z+ABX*I_NAI_G2x2z_F3z;
  Double I_NAI_Gx3y_Gx3z = I_NAI_H2x3y_F3z+ABX*I_NAI_Gx3y_F3z;
  Double I_NAI_Gx2yz_Gx3z = I_NAI_H2x2yz_F3z+ABX*I_NAI_Gx2yz_F3z;
  Double I_NAI_Gxy2z_Gx3z = I_NAI_H2xy2z_F3z+ABX*I_NAI_Gxy2z_F3z;
  Double I_NAI_Gx3z_Gx3z = I_NAI_H2x3z_F3z+ABX*I_NAI_Gx3z_F3z;
  Double I_NAI_G4y_Gx3z = I_NAI_Hx4y_F3z+ABX*I_NAI_G4y_F3z;
  Double I_NAI_G3yz_Gx3z = I_NAI_Hx3yz_F3z+ABX*I_NAI_G3yz_F3z;
  Double I_NAI_G2y2z_Gx3z = I_NAI_Hx2y2z_F3z+ABX*I_NAI_G2y2z_F3z;
  Double I_NAI_Gy3z_Gx3z = I_NAI_Hxy3z_F3z+ABX*I_NAI_Gy3z_F3z;
  Double I_NAI_G4z_Gx3z = I_NAI_Hx4z_F3z+ABX*I_NAI_G4z_F3z;
  Double I_NAI_G4x_G4y = I_NAI_H4xy_F3y+ABY*I_NAI_G4x_F3y;
  Double I_NAI_G3xy_G4y = I_NAI_H3x2y_F3y+ABY*I_NAI_G3xy_F3y;
  Double I_NAI_G3xz_G4y = I_NAI_H3xyz_F3y+ABY*I_NAI_G3xz_F3y;
  Double I_NAI_G2x2y_G4y = I_NAI_H2x3y_F3y+ABY*I_NAI_G2x2y_F3y;
  Double I_NAI_G2xyz_G4y = I_NAI_H2x2yz_F3y+ABY*I_NAI_G2xyz_F3y;
  Double I_NAI_G2x2z_G4y = I_NAI_H2xy2z_F3y+ABY*I_NAI_G2x2z_F3y;
  Double I_NAI_Gx3y_G4y = I_NAI_Hx4y_F3y+ABY*I_NAI_Gx3y_F3y;
  Double I_NAI_Gx2yz_G4y = I_NAI_Hx3yz_F3y+ABY*I_NAI_Gx2yz_F3y;
  Double I_NAI_Gxy2z_G4y = I_NAI_Hx2y2z_F3y+ABY*I_NAI_Gxy2z_F3y;
  Double I_NAI_Gx3z_G4y = I_NAI_Hxy3z_F3y+ABY*I_NAI_Gx3z_F3y;
  Double I_NAI_G4y_G4y = I_NAI_H5y_F3y+ABY*I_NAI_G4y_F3y;
  Double I_NAI_G3yz_G4y = I_NAI_H4yz_F3y+ABY*I_NAI_G3yz_F3y;
  Double I_NAI_G2y2z_G4y = I_NAI_H3y2z_F3y+ABY*I_NAI_G2y2z_F3y;
  Double I_NAI_Gy3z_G4y = I_NAI_H2y3z_F3y+ABY*I_NAI_Gy3z_F3y;
  Double I_NAI_G4z_G4y = I_NAI_Hy4z_F3y+ABY*I_NAI_G4z_F3y;
  Double I_NAI_G4x_G3yz = I_NAI_H4xz_F3y+ABZ*I_NAI_G4x_F3y;
  Double I_NAI_G3xy_G3yz = I_NAI_H3xyz_F3y+ABZ*I_NAI_G3xy_F3y;
  Double I_NAI_G3xz_G3yz = I_NAI_H3x2z_F3y+ABZ*I_NAI_G3xz_F3y;
  Double I_NAI_G2x2y_G3yz = I_NAI_H2x2yz_F3y+ABZ*I_NAI_G2x2y_F3y;
  Double I_NAI_G2xyz_G3yz = I_NAI_H2xy2z_F3y+ABZ*I_NAI_G2xyz_F3y;
  Double I_NAI_G2x2z_G3yz = I_NAI_H2x3z_F3y+ABZ*I_NAI_G2x2z_F3y;
  Double I_NAI_Gx3y_G3yz = I_NAI_Hx3yz_F3y+ABZ*I_NAI_Gx3y_F3y;
  Double I_NAI_Gx2yz_G3yz = I_NAI_Hx2y2z_F3y+ABZ*I_NAI_Gx2yz_F3y;
  Double I_NAI_Gxy2z_G3yz = I_NAI_Hxy3z_F3y+ABZ*I_NAI_Gxy2z_F3y;
  Double I_NAI_Gx3z_G3yz = I_NAI_Hx4z_F3y+ABZ*I_NAI_Gx3z_F3y;
  Double I_NAI_G4y_G3yz = I_NAI_H4yz_F3y+ABZ*I_NAI_G4y_F3y;
  Double I_NAI_G3yz_G3yz = I_NAI_H3y2z_F3y+ABZ*I_NAI_G3yz_F3y;
  Double I_NAI_G2y2z_G3yz = I_NAI_H2y3z_F3y+ABZ*I_NAI_G2y2z_F3y;
  Double I_NAI_Gy3z_G3yz = I_NAI_Hy4z_F3y+ABZ*I_NAI_Gy3z_F3y;
  Double I_NAI_G4z_G3yz = I_NAI_H5z_F3y+ABZ*I_NAI_G4z_F3y;
  Double I_NAI_G4x_G2y2z = I_NAI_H4xz_F2yz+ABZ*I_NAI_G4x_F2yz;
  Double I_NAI_G3xy_G2y2z = I_NAI_H3xyz_F2yz+ABZ*I_NAI_G3xy_F2yz;
  Double I_NAI_G3xz_G2y2z = I_NAI_H3x2z_F2yz+ABZ*I_NAI_G3xz_F2yz;
  Double I_NAI_G2x2y_G2y2z = I_NAI_H2x2yz_F2yz+ABZ*I_NAI_G2x2y_F2yz;
  Double I_NAI_G2xyz_G2y2z = I_NAI_H2xy2z_F2yz+ABZ*I_NAI_G2xyz_F2yz;
  Double I_NAI_G2x2z_G2y2z = I_NAI_H2x3z_F2yz+ABZ*I_NAI_G2x2z_F2yz;
  Double I_NAI_Gx3y_G2y2z = I_NAI_Hx3yz_F2yz+ABZ*I_NAI_Gx3y_F2yz;
  Double I_NAI_Gx2yz_G2y2z = I_NAI_Hx2y2z_F2yz+ABZ*I_NAI_Gx2yz_F2yz;
  Double I_NAI_Gxy2z_G2y2z = I_NAI_Hxy3z_F2yz+ABZ*I_NAI_Gxy2z_F2yz;
  Double I_NAI_Gx3z_G2y2z = I_NAI_Hx4z_F2yz+ABZ*I_NAI_Gx3z_F2yz;
  Double I_NAI_G4y_G2y2z = I_NAI_H4yz_F2yz+ABZ*I_NAI_G4y_F2yz;
  Double I_NAI_G3yz_G2y2z = I_NAI_H3y2z_F2yz+ABZ*I_NAI_G3yz_F2yz;
  Double I_NAI_G2y2z_G2y2z = I_NAI_H2y3z_F2yz+ABZ*I_NAI_G2y2z_F2yz;
  Double I_NAI_Gy3z_G2y2z = I_NAI_Hy4z_F2yz+ABZ*I_NAI_Gy3z_F2yz;
  Double I_NAI_G4z_G2y2z = I_NAI_H5z_F2yz+ABZ*I_NAI_G4z_F2yz;
  Double I_NAI_G4x_Gy3z = I_NAI_H4xy_F3z+ABY*I_NAI_G4x_F3z;
  Double I_NAI_G3xy_Gy3z = I_NAI_H3x2y_F3z+ABY*I_NAI_G3xy_F3z;
  Double I_NAI_G3xz_Gy3z = I_NAI_H3xyz_F3z+ABY*I_NAI_G3xz_F3z;
  Double I_NAI_G2x2y_Gy3z = I_NAI_H2x3y_F3z+ABY*I_NAI_G2x2y_F3z;
  Double I_NAI_G2xyz_Gy3z = I_NAI_H2x2yz_F3z+ABY*I_NAI_G2xyz_F3z;
  Double I_NAI_G2x2z_Gy3z = I_NAI_H2xy2z_F3z+ABY*I_NAI_G2x2z_F3z;
  Double I_NAI_Gx3y_Gy3z = I_NAI_Hx4y_F3z+ABY*I_NAI_Gx3y_F3z;
  Double I_NAI_Gx2yz_Gy3z = I_NAI_Hx3yz_F3z+ABY*I_NAI_Gx2yz_F3z;
  Double I_NAI_Gxy2z_Gy3z = I_NAI_Hx2y2z_F3z+ABY*I_NAI_Gxy2z_F3z;
  Double I_NAI_Gx3z_Gy3z = I_NAI_Hxy3z_F3z+ABY*I_NAI_Gx3z_F3z;
  Double I_NAI_G4y_Gy3z = I_NAI_H5y_F3z+ABY*I_NAI_G4y_F3z;
  Double I_NAI_G3yz_Gy3z = I_NAI_H4yz_F3z+ABY*I_NAI_G3yz_F3z;
  Double I_NAI_G2y2z_Gy3z = I_NAI_H3y2z_F3z+ABY*I_NAI_G2y2z_F3z;
  Double I_NAI_Gy3z_Gy3z = I_NAI_H2y3z_F3z+ABY*I_NAI_Gy3z_F3z;
  Double I_NAI_G4z_Gy3z = I_NAI_Hy4z_F3z+ABY*I_NAI_G4z_F3z;
  Double I_NAI_G4x_G4z = I_NAI_H4xz_F3z+ABZ*I_NAI_G4x_F3z;
  Double I_NAI_G3xy_G4z = I_NAI_H3xyz_F3z+ABZ*I_NAI_G3xy_F3z;
  Double I_NAI_G3xz_G4z = I_NAI_H3x2z_F3z+ABZ*I_NAI_G3xz_F3z;
  Double I_NAI_G2x2y_G4z = I_NAI_H2x2yz_F3z+ABZ*I_NAI_G2x2y_F3z;
  Double I_NAI_G2xyz_G4z = I_NAI_H2xy2z_F3z+ABZ*I_NAI_G2xyz_F3z;
  Double I_NAI_G2x2z_G4z = I_NAI_H2x3z_F3z+ABZ*I_NAI_G2x2z_F3z;
  Double I_NAI_Gx3y_G4z = I_NAI_Hx3yz_F3z+ABZ*I_NAI_Gx3y_F3z;
  Double I_NAI_Gx2yz_G4z = I_NAI_Hx2y2z_F3z+ABZ*I_NAI_Gx2yz_F3z;
  Double I_NAI_Gxy2z_G4z = I_NAI_Hxy3z_F3z+ABZ*I_NAI_Gxy2z_F3z;
  Double I_NAI_Gx3z_G4z = I_NAI_Hx4z_F3z+ABZ*I_NAI_Gx3z_F3z;
  Double I_NAI_G4y_G4z = I_NAI_H4yz_F3z+ABZ*I_NAI_G4y_F3z;
  Double I_NAI_G3yz_G4z = I_NAI_H3y2z_F3z+ABZ*I_NAI_G3yz_F3z;
  Double I_NAI_G2y2z_G4z = I_NAI_H2y3z_F3z+ABZ*I_NAI_G2y2z_F3z;
  Double I_NAI_Gy3z_G4z = I_NAI_Hy4z_F3z+ABZ*I_NAI_Gy3z_F3z;
  Double I_NAI_G4z_G4z = I_NAI_H5z_F3z+ABZ*I_NAI_G4z_F3z;

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
   * totally 84 integrals are omitted 
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
   * totally 3 integrals are omitted 
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
  Double I_NAI_L8y_Px_a = I_NAI_Mx8y_S_a+ABX*I_NAI_L8y_S_a;
  Double I_NAI_L7yz_Px_a = I_NAI_Mx7yz_S_a+ABX*I_NAI_L7yz_S_a;
  Double I_NAI_L6y2z_Px_a = I_NAI_Mx6y2z_S_a+ABX*I_NAI_L6y2z_S_a;
  Double I_NAI_L5y3z_Px_a = I_NAI_Mx5y3z_S_a+ABX*I_NAI_L5y3z_S_a;
  Double I_NAI_L4y4z_Px_a = I_NAI_Mx4y4z_S_a+ABX*I_NAI_L4y4z_S_a;
  Double I_NAI_L3y5z_Px_a = I_NAI_Mx3y5z_S_a+ABX*I_NAI_L3y5z_S_a;
  Double I_NAI_L2y6z_Px_a = I_NAI_Mx2y6z_S_a+ABX*I_NAI_L2y6z_S_a;
  Double I_NAI_Ly7z_Px_a = I_NAI_Mxy7z_S_a+ABX*I_NAI_Ly7z_S_a;
  Double I_NAI_L8z_Px_a = I_NAI_Mx8z_S_a+ABX*I_NAI_L8z_S_a;
  Double I_NAI_L7xy_Py_a = I_NAI_M7x2y_S_a+ABY*I_NAI_L7xy_S_a;
  Double I_NAI_L7xz_Py_a = I_NAI_M7xyz_S_a+ABY*I_NAI_L7xz_S_a;
  Double I_NAI_L6x2y_Py_a = I_NAI_M6x3y_S_a+ABY*I_NAI_L6x2y_S_a;
  Double I_NAI_L6xyz_Py_a = I_NAI_M6x2yz_S_a+ABY*I_NAI_L6xyz_S_a;
  Double I_NAI_L6x2z_Py_a = I_NAI_M6xy2z_S_a+ABY*I_NAI_L6x2z_S_a;
  Double I_NAI_L5x3y_Py_a = I_NAI_M5x4y_S_a+ABY*I_NAI_L5x3y_S_a;
  Double I_NAI_L5x2yz_Py_a = I_NAI_M5x3yz_S_a+ABY*I_NAI_L5x2yz_S_a;
  Double I_NAI_L5xy2z_Py_a = I_NAI_M5x2y2z_S_a+ABY*I_NAI_L5xy2z_S_a;
  Double I_NAI_L5x3z_Py_a = I_NAI_M5xy3z_S_a+ABY*I_NAI_L5x3z_S_a;
  Double I_NAI_L4x4y_Py_a = I_NAI_M4x5y_S_a+ABY*I_NAI_L4x4y_S_a;
  Double I_NAI_L4x3yz_Py_a = I_NAI_M4x4yz_S_a+ABY*I_NAI_L4x3yz_S_a;
  Double I_NAI_L4x2y2z_Py_a = I_NAI_M4x3y2z_S_a+ABY*I_NAI_L4x2y2z_S_a;
  Double I_NAI_L4xy3z_Py_a = I_NAI_M4x2y3z_S_a+ABY*I_NAI_L4xy3z_S_a;
  Double I_NAI_L4x4z_Py_a = I_NAI_M4xy4z_S_a+ABY*I_NAI_L4x4z_S_a;
  Double I_NAI_L3x5y_Py_a = I_NAI_M3x6y_S_a+ABY*I_NAI_L3x5y_S_a;
  Double I_NAI_L3x4yz_Py_a = I_NAI_M3x5yz_S_a+ABY*I_NAI_L3x4yz_S_a;
  Double I_NAI_L3x3y2z_Py_a = I_NAI_M3x4y2z_S_a+ABY*I_NAI_L3x3y2z_S_a;
  Double I_NAI_L3x2y3z_Py_a = I_NAI_M3x3y3z_S_a+ABY*I_NAI_L3x2y3z_S_a;
  Double I_NAI_L3xy4z_Py_a = I_NAI_M3x2y4z_S_a+ABY*I_NAI_L3xy4z_S_a;
  Double I_NAI_L3x5z_Py_a = I_NAI_M3xy5z_S_a+ABY*I_NAI_L3x5z_S_a;
  Double I_NAI_L2x6y_Py_a = I_NAI_M2x7y_S_a+ABY*I_NAI_L2x6y_S_a;
  Double I_NAI_L2x5yz_Py_a = I_NAI_M2x6yz_S_a+ABY*I_NAI_L2x5yz_S_a;
  Double I_NAI_L2x4y2z_Py_a = I_NAI_M2x5y2z_S_a+ABY*I_NAI_L2x4y2z_S_a;
  Double I_NAI_L2x3y3z_Py_a = I_NAI_M2x4y3z_S_a+ABY*I_NAI_L2x3y3z_S_a;
  Double I_NAI_L2x2y4z_Py_a = I_NAI_M2x3y4z_S_a+ABY*I_NAI_L2x2y4z_S_a;
  Double I_NAI_L2xy5z_Py_a = I_NAI_M2x2y5z_S_a+ABY*I_NAI_L2xy5z_S_a;
  Double I_NAI_L2x6z_Py_a = I_NAI_M2xy6z_S_a+ABY*I_NAI_L2x6z_S_a;
  Double I_NAI_Lx7y_Py_a = I_NAI_Mx8y_S_a+ABY*I_NAI_Lx7y_S_a;
  Double I_NAI_Lx6yz_Py_a = I_NAI_Mx7yz_S_a+ABY*I_NAI_Lx6yz_S_a;
  Double I_NAI_Lx5y2z_Py_a = I_NAI_Mx6y2z_S_a+ABY*I_NAI_Lx5y2z_S_a;
  Double I_NAI_Lx4y3z_Py_a = I_NAI_Mx5y3z_S_a+ABY*I_NAI_Lx4y3z_S_a;
  Double I_NAI_Lx3y4z_Py_a = I_NAI_Mx4y4z_S_a+ABY*I_NAI_Lx3y4z_S_a;
  Double I_NAI_Lx2y5z_Py_a = I_NAI_Mx3y5z_S_a+ABY*I_NAI_Lx2y5z_S_a;
  Double I_NAI_Lxy6z_Py_a = I_NAI_Mx2y6z_S_a+ABY*I_NAI_Lxy6z_S_a;
  Double I_NAI_Lx7z_Py_a = I_NAI_Mxy7z_S_a+ABY*I_NAI_Lx7z_S_a;
  Double I_NAI_L8y_Py_a = I_NAI_M9y_S_a+ABY*I_NAI_L8y_S_a;
  Double I_NAI_L7yz_Py_a = I_NAI_M8yz_S_a+ABY*I_NAI_L7yz_S_a;
  Double I_NAI_L6y2z_Py_a = I_NAI_M7y2z_S_a+ABY*I_NAI_L6y2z_S_a;
  Double I_NAI_L5y3z_Py_a = I_NAI_M6y3z_S_a+ABY*I_NAI_L5y3z_S_a;
  Double I_NAI_L4y4z_Py_a = I_NAI_M5y4z_S_a+ABY*I_NAI_L4y4z_S_a;
  Double I_NAI_L3y5z_Py_a = I_NAI_M4y5z_S_a+ABY*I_NAI_L3y5z_S_a;
  Double I_NAI_L2y6z_Py_a = I_NAI_M3y6z_S_a+ABY*I_NAI_L2y6z_S_a;
  Double I_NAI_Ly7z_Py_a = I_NAI_M2y7z_S_a+ABY*I_NAI_Ly7z_S_a;
  Double I_NAI_L8z_Py_a = I_NAI_My8z_S_a+ABY*I_NAI_L8z_S_a;
  Double I_NAI_L7xy_Pz_a = I_NAI_M7xyz_S_a+ABZ*I_NAI_L7xy_S_a;
  Double I_NAI_L7xz_Pz_a = I_NAI_M7x2z_S_a+ABZ*I_NAI_L7xz_S_a;
  Double I_NAI_L6x2y_Pz_a = I_NAI_M6x2yz_S_a+ABZ*I_NAI_L6x2y_S_a;
  Double I_NAI_L6xyz_Pz_a = I_NAI_M6xy2z_S_a+ABZ*I_NAI_L6xyz_S_a;
  Double I_NAI_L6x2z_Pz_a = I_NAI_M6x3z_S_a+ABZ*I_NAI_L6x2z_S_a;
  Double I_NAI_L5x3y_Pz_a = I_NAI_M5x3yz_S_a+ABZ*I_NAI_L5x3y_S_a;
  Double I_NAI_L5x2yz_Pz_a = I_NAI_M5x2y2z_S_a+ABZ*I_NAI_L5x2yz_S_a;
  Double I_NAI_L5xy2z_Pz_a = I_NAI_M5xy3z_S_a+ABZ*I_NAI_L5xy2z_S_a;
  Double I_NAI_L5x3z_Pz_a = I_NAI_M5x4z_S_a+ABZ*I_NAI_L5x3z_S_a;
  Double I_NAI_L4x4y_Pz_a = I_NAI_M4x4yz_S_a+ABZ*I_NAI_L4x4y_S_a;
  Double I_NAI_L4x3yz_Pz_a = I_NAI_M4x3y2z_S_a+ABZ*I_NAI_L4x3yz_S_a;
  Double I_NAI_L4x2y2z_Pz_a = I_NAI_M4x2y3z_S_a+ABZ*I_NAI_L4x2y2z_S_a;
  Double I_NAI_L4xy3z_Pz_a = I_NAI_M4xy4z_S_a+ABZ*I_NAI_L4xy3z_S_a;
  Double I_NAI_L4x4z_Pz_a = I_NAI_M4x5z_S_a+ABZ*I_NAI_L4x4z_S_a;
  Double I_NAI_L3x5y_Pz_a = I_NAI_M3x5yz_S_a+ABZ*I_NAI_L3x5y_S_a;
  Double I_NAI_L3x4yz_Pz_a = I_NAI_M3x4y2z_S_a+ABZ*I_NAI_L3x4yz_S_a;
  Double I_NAI_L3x3y2z_Pz_a = I_NAI_M3x3y3z_S_a+ABZ*I_NAI_L3x3y2z_S_a;
  Double I_NAI_L3x2y3z_Pz_a = I_NAI_M3x2y4z_S_a+ABZ*I_NAI_L3x2y3z_S_a;
  Double I_NAI_L3xy4z_Pz_a = I_NAI_M3xy5z_S_a+ABZ*I_NAI_L3xy4z_S_a;
  Double I_NAI_L3x5z_Pz_a = I_NAI_M3x6z_S_a+ABZ*I_NAI_L3x5z_S_a;
  Double I_NAI_L2x6y_Pz_a = I_NAI_M2x6yz_S_a+ABZ*I_NAI_L2x6y_S_a;
  Double I_NAI_L2x5yz_Pz_a = I_NAI_M2x5y2z_S_a+ABZ*I_NAI_L2x5yz_S_a;
  Double I_NAI_L2x4y2z_Pz_a = I_NAI_M2x4y3z_S_a+ABZ*I_NAI_L2x4y2z_S_a;
  Double I_NAI_L2x3y3z_Pz_a = I_NAI_M2x3y4z_S_a+ABZ*I_NAI_L2x3y3z_S_a;
  Double I_NAI_L2x2y4z_Pz_a = I_NAI_M2x2y5z_S_a+ABZ*I_NAI_L2x2y4z_S_a;
  Double I_NAI_L2xy5z_Pz_a = I_NAI_M2xy6z_S_a+ABZ*I_NAI_L2xy5z_S_a;
  Double I_NAI_L2x6z_Pz_a = I_NAI_M2x7z_S_a+ABZ*I_NAI_L2x6z_S_a;
  Double I_NAI_Lx7y_Pz_a = I_NAI_Mx7yz_S_a+ABZ*I_NAI_Lx7y_S_a;
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
   * totally 108 integrals are omitted 
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
   * totally 56 integrals are omitted 
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
   * shell quartet name: SQ_NAI_M_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 33 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_N_S_a
   * RHS shell quartet name: SQ_NAI_M_S_a
   ************************************************************/
  Double I_NAI_M9x_Px_a = I_NAI_N10x_S_a+ABX*I_NAI_M9x_S_a;
  Double I_NAI_M8xy_Px_a = I_NAI_N9xy_S_a+ABX*I_NAI_M8xy_S_a;
  Double I_NAI_M8xz_Px_a = I_NAI_N9xz_S_a+ABX*I_NAI_M8xz_S_a;
  Double I_NAI_M7x2y_Px_a = I_NAI_N8x2y_S_a+ABX*I_NAI_M7x2y_S_a;
  Double I_NAI_M7xyz_Px_a = I_NAI_N8xyz_S_a+ABX*I_NAI_M7xyz_S_a;
  Double I_NAI_M7x2z_Px_a = I_NAI_N8x2z_S_a+ABX*I_NAI_M7x2z_S_a;
  Double I_NAI_M6x3y_Px_a = I_NAI_N7x3y_S_a+ABX*I_NAI_M6x3y_S_a;
  Double I_NAI_M6x2yz_Px_a = I_NAI_N7x2yz_S_a+ABX*I_NAI_M6x2yz_S_a;
  Double I_NAI_M6xy2z_Px_a = I_NAI_N7xy2z_S_a+ABX*I_NAI_M6xy2z_S_a;
  Double I_NAI_M6x3z_Px_a = I_NAI_N7x3z_S_a+ABX*I_NAI_M6x3z_S_a;
  Double I_NAI_M5x4y_Px_a = I_NAI_N6x4y_S_a+ABX*I_NAI_M5x4y_S_a;
  Double I_NAI_M5x3yz_Px_a = I_NAI_N6x3yz_S_a+ABX*I_NAI_M5x3yz_S_a;
  Double I_NAI_M5x2y2z_Px_a = I_NAI_N6x2y2z_S_a+ABX*I_NAI_M5x2y2z_S_a;
  Double I_NAI_M5xy3z_Px_a = I_NAI_N6xy3z_S_a+ABX*I_NAI_M5xy3z_S_a;
  Double I_NAI_M5x4z_Px_a = I_NAI_N6x4z_S_a+ABX*I_NAI_M5x4z_S_a;
  Double I_NAI_M4x5y_Px_a = I_NAI_N5x5y_S_a+ABX*I_NAI_M4x5y_S_a;
  Double I_NAI_M4x4yz_Px_a = I_NAI_N5x4yz_S_a+ABX*I_NAI_M4x4yz_S_a;
  Double I_NAI_M4x3y2z_Px_a = I_NAI_N5x3y2z_S_a+ABX*I_NAI_M4x3y2z_S_a;
  Double I_NAI_M4x2y3z_Px_a = I_NAI_N5x2y3z_S_a+ABX*I_NAI_M4x2y3z_S_a;
  Double I_NAI_M4xy4z_Px_a = I_NAI_N5xy4z_S_a+ABX*I_NAI_M4xy4z_S_a;
  Double I_NAI_M4x5z_Px_a = I_NAI_N5x5z_S_a+ABX*I_NAI_M4x5z_S_a;
  Double I_NAI_M3x6y_Px_a = I_NAI_N4x6y_S_a+ABX*I_NAI_M3x6y_S_a;
  Double I_NAI_M3x5yz_Px_a = I_NAI_N4x5yz_S_a+ABX*I_NAI_M3x5yz_S_a;
  Double I_NAI_M3x4y2z_Px_a = I_NAI_N4x4y2z_S_a+ABX*I_NAI_M3x4y2z_S_a;
  Double I_NAI_M3x3y3z_Px_a = I_NAI_N4x3y3z_S_a+ABX*I_NAI_M3x3y3z_S_a;
  Double I_NAI_M3x2y4z_Px_a = I_NAI_N4x2y4z_S_a+ABX*I_NAI_M3x2y4z_S_a;
  Double I_NAI_M3xy5z_Px_a = I_NAI_N4xy5z_S_a+ABX*I_NAI_M3xy5z_S_a;
  Double I_NAI_M3x6z_Px_a = I_NAI_N4x6z_S_a+ABX*I_NAI_M3x6z_S_a;
  Double I_NAI_M2x7y_Px_a = I_NAI_N3x7y_S_a+ABX*I_NAI_M2x7y_S_a;
  Double I_NAI_M2x6yz_Px_a = I_NAI_N3x6yz_S_a+ABX*I_NAI_M2x6yz_S_a;
  Double I_NAI_M2x5y2z_Px_a = I_NAI_N3x5y2z_S_a+ABX*I_NAI_M2x5y2z_S_a;
  Double I_NAI_M2x4y3z_Px_a = I_NAI_N3x4y3z_S_a+ABX*I_NAI_M2x4y3z_S_a;
  Double I_NAI_M2x3y4z_Px_a = I_NAI_N3x3y4z_S_a+ABX*I_NAI_M2x3y4z_S_a;
  Double I_NAI_M2x2y5z_Px_a = I_NAI_N3x2y5z_S_a+ABX*I_NAI_M2x2y5z_S_a;
  Double I_NAI_M2xy6z_Px_a = I_NAI_N3xy6z_S_a+ABX*I_NAI_M2xy6z_S_a;
  Double I_NAI_M2x7z_Px_a = I_NAI_N3x7z_S_a+ABX*I_NAI_M2x7z_S_a;
  Double I_NAI_Mx8y_Px_a = I_NAI_N2x8y_S_a+ABX*I_NAI_Mx8y_S_a;
  Double I_NAI_Mx7yz_Px_a = I_NAI_N2x7yz_S_a+ABX*I_NAI_Mx7yz_S_a;
  Double I_NAI_Mx6y2z_Px_a = I_NAI_N2x6y2z_S_a+ABX*I_NAI_Mx6y2z_S_a;
  Double I_NAI_Mx5y3z_Px_a = I_NAI_N2x5y3z_S_a+ABX*I_NAI_Mx5y3z_S_a;
  Double I_NAI_Mx4y4z_Px_a = I_NAI_N2x4y4z_S_a+ABX*I_NAI_Mx4y4z_S_a;
  Double I_NAI_Mx3y5z_Px_a = I_NAI_N2x3y5z_S_a+ABX*I_NAI_Mx3y5z_S_a;
  Double I_NAI_Mx2y6z_Px_a = I_NAI_N2x2y6z_S_a+ABX*I_NAI_Mx2y6z_S_a;
  Double I_NAI_Mxy7z_Px_a = I_NAI_N2xy7z_S_a+ABX*I_NAI_Mxy7z_S_a;
  Double I_NAI_Mx8z_Px_a = I_NAI_N2x8z_S_a+ABX*I_NAI_Mx8z_S_a;
  Double I_NAI_M7x2y_Py_a = I_NAI_N7x3y_S_a+ABY*I_NAI_M7x2y_S_a;
  Double I_NAI_M7xyz_Py_a = I_NAI_N7x2yz_S_a+ABY*I_NAI_M7xyz_S_a;
  Double I_NAI_M6x3y_Py_a = I_NAI_N6x4y_S_a+ABY*I_NAI_M6x3y_S_a;
  Double I_NAI_M6x2yz_Py_a = I_NAI_N6x3yz_S_a+ABY*I_NAI_M6x2yz_S_a;
  Double I_NAI_M6xy2z_Py_a = I_NAI_N6x2y2z_S_a+ABY*I_NAI_M6xy2z_S_a;
  Double I_NAI_M5x4y_Py_a = I_NAI_N5x5y_S_a+ABY*I_NAI_M5x4y_S_a;
  Double I_NAI_M5x3yz_Py_a = I_NAI_N5x4yz_S_a+ABY*I_NAI_M5x3yz_S_a;
  Double I_NAI_M5x2y2z_Py_a = I_NAI_N5x3y2z_S_a+ABY*I_NAI_M5x2y2z_S_a;
  Double I_NAI_M5xy3z_Py_a = I_NAI_N5x2y3z_S_a+ABY*I_NAI_M5xy3z_S_a;
  Double I_NAI_M4x5y_Py_a = I_NAI_N4x6y_S_a+ABY*I_NAI_M4x5y_S_a;
  Double I_NAI_M4x4yz_Py_a = I_NAI_N4x5yz_S_a+ABY*I_NAI_M4x4yz_S_a;
  Double I_NAI_M4x3y2z_Py_a = I_NAI_N4x4y2z_S_a+ABY*I_NAI_M4x3y2z_S_a;
  Double I_NAI_M4x2y3z_Py_a = I_NAI_N4x3y3z_S_a+ABY*I_NAI_M4x2y3z_S_a;
  Double I_NAI_M4xy4z_Py_a = I_NAI_N4x2y4z_S_a+ABY*I_NAI_M4xy4z_S_a;
  Double I_NAI_M3x6y_Py_a = I_NAI_N3x7y_S_a+ABY*I_NAI_M3x6y_S_a;
  Double I_NAI_M3x5yz_Py_a = I_NAI_N3x6yz_S_a+ABY*I_NAI_M3x5yz_S_a;
  Double I_NAI_M3x4y2z_Py_a = I_NAI_N3x5y2z_S_a+ABY*I_NAI_M3x4y2z_S_a;
  Double I_NAI_M3x3y3z_Py_a = I_NAI_N3x4y3z_S_a+ABY*I_NAI_M3x3y3z_S_a;
  Double I_NAI_M3x2y4z_Py_a = I_NAI_N3x3y4z_S_a+ABY*I_NAI_M3x2y4z_S_a;
  Double I_NAI_M3xy5z_Py_a = I_NAI_N3x2y5z_S_a+ABY*I_NAI_M3xy5z_S_a;
  Double I_NAI_M2x7y_Py_a = I_NAI_N2x8y_S_a+ABY*I_NAI_M2x7y_S_a;
  Double I_NAI_M2x6yz_Py_a = I_NAI_N2x7yz_S_a+ABY*I_NAI_M2x6yz_S_a;
  Double I_NAI_M2x5y2z_Py_a = I_NAI_N2x6y2z_S_a+ABY*I_NAI_M2x5y2z_S_a;
  Double I_NAI_M2x4y3z_Py_a = I_NAI_N2x5y3z_S_a+ABY*I_NAI_M2x4y3z_S_a;
  Double I_NAI_M2x3y4z_Py_a = I_NAI_N2x4y4z_S_a+ABY*I_NAI_M2x3y4z_S_a;
  Double I_NAI_M2x2y5z_Py_a = I_NAI_N2x3y5z_S_a+ABY*I_NAI_M2x2y5z_S_a;
  Double I_NAI_M2xy6z_Py_a = I_NAI_N2x2y6z_S_a+ABY*I_NAI_M2xy6z_S_a;
  Double I_NAI_Mx8y_Py_a = I_NAI_Nx9y_S_a+ABY*I_NAI_Mx8y_S_a;
  Double I_NAI_Mx7yz_Py_a = I_NAI_Nx8yz_S_a+ABY*I_NAI_Mx7yz_S_a;
  Double I_NAI_Mx6y2z_Py_a = I_NAI_Nx7y2z_S_a+ABY*I_NAI_Mx6y2z_S_a;
  Double I_NAI_Mx5y3z_Py_a = I_NAI_Nx6y3z_S_a+ABY*I_NAI_Mx5y3z_S_a;
  Double I_NAI_Mx4y4z_Py_a = I_NAI_Nx5y4z_S_a+ABY*I_NAI_Mx4y4z_S_a;
  Double I_NAI_Mx3y5z_Py_a = I_NAI_Nx4y5z_S_a+ABY*I_NAI_Mx3y5z_S_a;
  Double I_NAI_Mx2y6z_Py_a = I_NAI_Nx3y6z_S_a+ABY*I_NAI_Mx2y6z_S_a;
  Double I_NAI_Mxy7z_Py_a = I_NAI_Nx2y7z_S_a+ABY*I_NAI_Mxy7z_S_a;
  Double I_NAI_M9y_Py_a = I_NAI_N10y_S_a+ABY*I_NAI_M9y_S_a;
  Double I_NAI_M8yz_Py_a = I_NAI_N9yz_S_a+ABY*I_NAI_M8yz_S_a;
  Double I_NAI_M7y2z_Py_a = I_NAI_N8y2z_S_a+ABY*I_NAI_M7y2z_S_a;
  Double I_NAI_M6y3z_Py_a = I_NAI_N7y3z_S_a+ABY*I_NAI_M6y3z_S_a;
  Double I_NAI_M5y4z_Py_a = I_NAI_N6y4z_S_a+ABY*I_NAI_M5y4z_S_a;
  Double I_NAI_M4y5z_Py_a = I_NAI_N5y5z_S_a+ABY*I_NAI_M4y5z_S_a;
  Double I_NAI_M3y6z_Py_a = I_NAI_N4y6z_S_a+ABY*I_NAI_M3y6z_S_a;
  Double I_NAI_M2y7z_Py_a = I_NAI_N3y7z_S_a+ABY*I_NAI_M2y7z_S_a;
  Double I_NAI_My8z_Py_a = I_NAI_N2y8z_S_a+ABY*I_NAI_My8z_S_a;
  Double I_NAI_M7xyz_Pz_a = I_NAI_N7xy2z_S_a+ABZ*I_NAI_M7xyz_S_a;
  Double I_NAI_M7x2z_Pz_a = I_NAI_N7x3z_S_a+ABZ*I_NAI_M7x2z_S_a;
  Double I_NAI_M6x2yz_Pz_a = I_NAI_N6x2y2z_S_a+ABZ*I_NAI_M6x2yz_S_a;
  Double I_NAI_M6xy2z_Pz_a = I_NAI_N6xy3z_S_a+ABZ*I_NAI_M6xy2z_S_a;
  Double I_NAI_M6x3z_Pz_a = I_NAI_N6x4z_S_a+ABZ*I_NAI_M6x3z_S_a;
  Double I_NAI_M5x3yz_Pz_a = I_NAI_N5x3y2z_S_a+ABZ*I_NAI_M5x3yz_S_a;
  Double I_NAI_M5x2y2z_Pz_a = I_NAI_N5x2y3z_S_a+ABZ*I_NAI_M5x2y2z_S_a;
  Double I_NAI_M5xy3z_Pz_a = I_NAI_N5xy4z_S_a+ABZ*I_NAI_M5xy3z_S_a;
  Double I_NAI_M5x4z_Pz_a = I_NAI_N5x5z_S_a+ABZ*I_NAI_M5x4z_S_a;
  Double I_NAI_M4x4yz_Pz_a = I_NAI_N4x4y2z_S_a+ABZ*I_NAI_M4x4yz_S_a;
  Double I_NAI_M4x3y2z_Pz_a = I_NAI_N4x3y3z_S_a+ABZ*I_NAI_M4x3y2z_S_a;
  Double I_NAI_M4x2y3z_Pz_a = I_NAI_N4x2y4z_S_a+ABZ*I_NAI_M4x2y3z_S_a;
  Double I_NAI_M4xy4z_Pz_a = I_NAI_N4xy5z_S_a+ABZ*I_NAI_M4xy4z_S_a;
  Double I_NAI_M4x5z_Pz_a = I_NAI_N4x6z_S_a+ABZ*I_NAI_M4x5z_S_a;
  Double I_NAI_M3x5yz_Pz_a = I_NAI_N3x5y2z_S_a+ABZ*I_NAI_M3x5yz_S_a;
  Double I_NAI_M3x4y2z_Pz_a = I_NAI_N3x4y3z_S_a+ABZ*I_NAI_M3x4y2z_S_a;
  Double I_NAI_M3x3y3z_Pz_a = I_NAI_N3x3y4z_S_a+ABZ*I_NAI_M3x3y3z_S_a;
  Double I_NAI_M3x2y4z_Pz_a = I_NAI_N3x2y5z_S_a+ABZ*I_NAI_M3x2y4z_S_a;
  Double I_NAI_M3xy5z_Pz_a = I_NAI_N3xy6z_S_a+ABZ*I_NAI_M3xy5z_S_a;
  Double I_NAI_M3x6z_Pz_a = I_NAI_N3x7z_S_a+ABZ*I_NAI_M3x6z_S_a;
  Double I_NAI_M2x6yz_Pz_a = I_NAI_N2x6y2z_S_a+ABZ*I_NAI_M2x6yz_S_a;
  Double I_NAI_M2x5y2z_Pz_a = I_NAI_N2x5y3z_S_a+ABZ*I_NAI_M2x5y2z_S_a;
  Double I_NAI_M2x4y3z_Pz_a = I_NAI_N2x4y4z_S_a+ABZ*I_NAI_M2x4y3z_S_a;
  Double I_NAI_M2x3y4z_Pz_a = I_NAI_N2x3y5z_S_a+ABZ*I_NAI_M2x3y4z_S_a;
  Double I_NAI_M2x2y5z_Pz_a = I_NAI_N2x2y6z_S_a+ABZ*I_NAI_M2x2y5z_S_a;
  Double I_NAI_M2xy6z_Pz_a = I_NAI_N2xy7z_S_a+ABZ*I_NAI_M2xy6z_S_a;
  Double I_NAI_M2x7z_Pz_a = I_NAI_N2x8z_S_a+ABZ*I_NAI_M2x7z_S_a;
  Double I_NAI_Mx7yz_Pz_a = I_NAI_Nx7y2z_S_a+ABZ*I_NAI_Mx7yz_S_a;
  Double I_NAI_Mx6y2z_Pz_a = I_NAI_Nx6y3z_S_a+ABZ*I_NAI_Mx6y2z_S_a;
  Double I_NAI_Mx5y3z_Pz_a = I_NAI_Nx5y4z_S_a+ABZ*I_NAI_Mx5y3z_S_a;
  Double I_NAI_Mx4y4z_Pz_a = I_NAI_Nx4y5z_S_a+ABZ*I_NAI_Mx4y4z_S_a;
  Double I_NAI_Mx3y5z_Pz_a = I_NAI_Nx3y6z_S_a+ABZ*I_NAI_Mx3y5z_S_a;
  Double I_NAI_Mx2y6z_Pz_a = I_NAI_Nx2y7z_S_a+ABZ*I_NAI_Mx2y6z_S_a;
  Double I_NAI_Mxy7z_Pz_a = I_NAI_Nxy8z_S_a+ABZ*I_NAI_Mxy7z_S_a;
  Double I_NAI_Mx8z_Pz_a = I_NAI_Nx9z_S_a+ABZ*I_NAI_Mx8z_S_a;
  Double I_NAI_M7y2z_Pz_a = I_NAI_N7y3z_S_a+ABZ*I_NAI_M7y2z_S_a;
  Double I_NAI_M6y3z_Pz_a = I_NAI_N6y4z_S_a+ABZ*I_NAI_M6y3z_S_a;
  Double I_NAI_M5y4z_Pz_a = I_NAI_N5y5z_S_a+ABZ*I_NAI_M5y4z_S_a;
  Double I_NAI_M4y5z_Pz_a = I_NAI_N4y6z_S_a+ABZ*I_NAI_M4y5z_S_a;
  Double I_NAI_M3y6z_Pz_a = I_NAI_N3y7z_S_a+ABZ*I_NAI_M3y6z_S_a;
  Double I_NAI_M2y7z_Pz_a = I_NAI_N2y8z_S_a+ABZ*I_NAI_M2y7z_S_a;
  Double I_NAI_My8z_Pz_a = I_NAI_Ny9z_S_a+ABZ*I_NAI_My8z_S_a;
  Double I_NAI_M9z_Pz_a = I_NAI_N10z_S_a+ABZ*I_NAI_M9z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_L_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 138 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_P_a
   * RHS shell quartet name: SQ_NAI_L_P_a
   ************************************************************/
  Double I_NAI_L8x_D2x_a = I_NAI_M9x_Px_a+ABX*I_NAI_L8x_Px_a;
  Double I_NAI_L7xy_D2x_a = I_NAI_M8xy_Px_a+ABX*I_NAI_L7xy_Px_a;
  Double I_NAI_L7xz_D2x_a = I_NAI_M8xz_Px_a+ABX*I_NAI_L7xz_Px_a;
  Double I_NAI_L6x2y_D2x_a = I_NAI_M7x2y_Px_a+ABX*I_NAI_L6x2y_Px_a;
  Double I_NAI_L6xyz_D2x_a = I_NAI_M7xyz_Px_a+ABX*I_NAI_L6xyz_Px_a;
  Double I_NAI_L6x2z_D2x_a = I_NAI_M7x2z_Px_a+ABX*I_NAI_L6x2z_Px_a;
  Double I_NAI_L5x3y_D2x_a = I_NAI_M6x3y_Px_a+ABX*I_NAI_L5x3y_Px_a;
  Double I_NAI_L5x2yz_D2x_a = I_NAI_M6x2yz_Px_a+ABX*I_NAI_L5x2yz_Px_a;
  Double I_NAI_L5xy2z_D2x_a = I_NAI_M6xy2z_Px_a+ABX*I_NAI_L5xy2z_Px_a;
  Double I_NAI_L5x3z_D2x_a = I_NAI_M6x3z_Px_a+ABX*I_NAI_L5x3z_Px_a;
  Double I_NAI_L4x4y_D2x_a = I_NAI_M5x4y_Px_a+ABX*I_NAI_L4x4y_Px_a;
  Double I_NAI_L4x3yz_D2x_a = I_NAI_M5x3yz_Px_a+ABX*I_NAI_L4x3yz_Px_a;
  Double I_NAI_L4x2y2z_D2x_a = I_NAI_M5x2y2z_Px_a+ABX*I_NAI_L4x2y2z_Px_a;
  Double I_NAI_L4xy3z_D2x_a = I_NAI_M5xy3z_Px_a+ABX*I_NAI_L4xy3z_Px_a;
  Double I_NAI_L4x4z_D2x_a = I_NAI_M5x4z_Px_a+ABX*I_NAI_L4x4z_Px_a;
  Double I_NAI_L3x5y_D2x_a = I_NAI_M4x5y_Px_a+ABX*I_NAI_L3x5y_Px_a;
  Double I_NAI_L3x4yz_D2x_a = I_NAI_M4x4yz_Px_a+ABX*I_NAI_L3x4yz_Px_a;
  Double I_NAI_L3x3y2z_D2x_a = I_NAI_M4x3y2z_Px_a+ABX*I_NAI_L3x3y2z_Px_a;
  Double I_NAI_L3x2y3z_D2x_a = I_NAI_M4x2y3z_Px_a+ABX*I_NAI_L3x2y3z_Px_a;
  Double I_NAI_L3xy4z_D2x_a = I_NAI_M4xy4z_Px_a+ABX*I_NAI_L3xy4z_Px_a;
  Double I_NAI_L3x5z_D2x_a = I_NAI_M4x5z_Px_a+ABX*I_NAI_L3x5z_Px_a;
  Double I_NAI_L2x6y_D2x_a = I_NAI_M3x6y_Px_a+ABX*I_NAI_L2x6y_Px_a;
  Double I_NAI_L2x5yz_D2x_a = I_NAI_M3x5yz_Px_a+ABX*I_NAI_L2x5yz_Px_a;
  Double I_NAI_L2x4y2z_D2x_a = I_NAI_M3x4y2z_Px_a+ABX*I_NAI_L2x4y2z_Px_a;
  Double I_NAI_L2x3y3z_D2x_a = I_NAI_M3x3y3z_Px_a+ABX*I_NAI_L2x3y3z_Px_a;
  Double I_NAI_L2x2y4z_D2x_a = I_NAI_M3x2y4z_Px_a+ABX*I_NAI_L2x2y4z_Px_a;
  Double I_NAI_L2xy5z_D2x_a = I_NAI_M3xy5z_Px_a+ABX*I_NAI_L2xy5z_Px_a;
  Double I_NAI_L2x6z_D2x_a = I_NAI_M3x6z_Px_a+ABX*I_NAI_L2x6z_Px_a;
  Double I_NAI_Lx7y_D2x_a = I_NAI_M2x7y_Px_a+ABX*I_NAI_Lx7y_Px_a;
  Double I_NAI_Lx6yz_D2x_a = I_NAI_M2x6yz_Px_a+ABX*I_NAI_Lx6yz_Px_a;
  Double I_NAI_Lx5y2z_D2x_a = I_NAI_M2x5y2z_Px_a+ABX*I_NAI_Lx5y2z_Px_a;
  Double I_NAI_Lx4y3z_D2x_a = I_NAI_M2x4y3z_Px_a+ABX*I_NAI_Lx4y3z_Px_a;
  Double I_NAI_Lx3y4z_D2x_a = I_NAI_M2x3y4z_Px_a+ABX*I_NAI_Lx3y4z_Px_a;
  Double I_NAI_Lx2y5z_D2x_a = I_NAI_M2x2y5z_Px_a+ABX*I_NAI_Lx2y5z_Px_a;
  Double I_NAI_Lxy6z_D2x_a = I_NAI_M2xy6z_Px_a+ABX*I_NAI_Lxy6z_Px_a;
  Double I_NAI_Lx7z_D2x_a = I_NAI_M2x7z_Px_a+ABX*I_NAI_Lx7z_Px_a;
  Double I_NAI_L8y_D2x_a = I_NAI_Mx8y_Px_a+ABX*I_NAI_L8y_Px_a;
  Double I_NAI_L7yz_D2x_a = I_NAI_Mx7yz_Px_a+ABX*I_NAI_L7yz_Px_a;
  Double I_NAI_L6y2z_D2x_a = I_NAI_Mx6y2z_Px_a+ABX*I_NAI_L6y2z_Px_a;
  Double I_NAI_L5y3z_D2x_a = I_NAI_Mx5y3z_Px_a+ABX*I_NAI_L5y3z_Px_a;
  Double I_NAI_L4y4z_D2x_a = I_NAI_Mx4y4z_Px_a+ABX*I_NAI_L4y4z_Px_a;
  Double I_NAI_L3y5z_D2x_a = I_NAI_Mx3y5z_Px_a+ABX*I_NAI_L3y5z_Px_a;
  Double I_NAI_L2y6z_D2x_a = I_NAI_Mx2y6z_Px_a+ABX*I_NAI_L2y6z_Px_a;
  Double I_NAI_Ly7z_D2x_a = I_NAI_Mxy7z_Px_a+ABX*I_NAI_Ly7z_Px_a;
  Double I_NAI_L8z_D2x_a = I_NAI_Mx8z_Px_a+ABX*I_NAI_L8z_Px_a;
  Double I_NAI_L7xy_D2y_a = I_NAI_M7x2y_Py_a+ABY*I_NAI_L7xy_Py_a;
  Double I_NAI_L7xz_D2y_a = I_NAI_M7xyz_Py_a+ABY*I_NAI_L7xz_Py_a;
  Double I_NAI_L6x2y_D2y_a = I_NAI_M6x3y_Py_a+ABY*I_NAI_L6x2y_Py_a;
  Double I_NAI_L6xyz_D2y_a = I_NAI_M6x2yz_Py_a+ABY*I_NAI_L6xyz_Py_a;
  Double I_NAI_L6x2z_D2y_a = I_NAI_M6xy2z_Py_a+ABY*I_NAI_L6x2z_Py_a;
  Double I_NAI_L5x3y_D2y_a = I_NAI_M5x4y_Py_a+ABY*I_NAI_L5x3y_Py_a;
  Double I_NAI_L5x2yz_D2y_a = I_NAI_M5x3yz_Py_a+ABY*I_NAI_L5x2yz_Py_a;
  Double I_NAI_L5xy2z_D2y_a = I_NAI_M5x2y2z_Py_a+ABY*I_NAI_L5xy2z_Py_a;
  Double I_NAI_L5x3z_D2y_a = I_NAI_M5xy3z_Py_a+ABY*I_NAI_L5x3z_Py_a;
  Double I_NAI_L4x4y_D2y_a = I_NAI_M4x5y_Py_a+ABY*I_NAI_L4x4y_Py_a;
  Double I_NAI_L4x3yz_D2y_a = I_NAI_M4x4yz_Py_a+ABY*I_NAI_L4x3yz_Py_a;
  Double I_NAI_L4x2y2z_D2y_a = I_NAI_M4x3y2z_Py_a+ABY*I_NAI_L4x2y2z_Py_a;
  Double I_NAI_L4xy3z_D2y_a = I_NAI_M4x2y3z_Py_a+ABY*I_NAI_L4xy3z_Py_a;
  Double I_NAI_L4x4z_D2y_a = I_NAI_M4xy4z_Py_a+ABY*I_NAI_L4x4z_Py_a;
  Double I_NAI_L3x5y_D2y_a = I_NAI_M3x6y_Py_a+ABY*I_NAI_L3x5y_Py_a;
  Double I_NAI_L3x4yz_D2y_a = I_NAI_M3x5yz_Py_a+ABY*I_NAI_L3x4yz_Py_a;
  Double I_NAI_L3x3y2z_D2y_a = I_NAI_M3x4y2z_Py_a+ABY*I_NAI_L3x3y2z_Py_a;
  Double I_NAI_L3x2y3z_D2y_a = I_NAI_M3x3y3z_Py_a+ABY*I_NAI_L3x2y3z_Py_a;
  Double I_NAI_L3xy4z_D2y_a = I_NAI_M3x2y4z_Py_a+ABY*I_NAI_L3xy4z_Py_a;
  Double I_NAI_L3x5z_D2y_a = I_NAI_M3xy5z_Py_a+ABY*I_NAI_L3x5z_Py_a;
  Double I_NAI_L2x6y_D2y_a = I_NAI_M2x7y_Py_a+ABY*I_NAI_L2x6y_Py_a;
  Double I_NAI_L2x5yz_D2y_a = I_NAI_M2x6yz_Py_a+ABY*I_NAI_L2x5yz_Py_a;
  Double I_NAI_L2x4y2z_D2y_a = I_NAI_M2x5y2z_Py_a+ABY*I_NAI_L2x4y2z_Py_a;
  Double I_NAI_L2x3y3z_D2y_a = I_NAI_M2x4y3z_Py_a+ABY*I_NAI_L2x3y3z_Py_a;
  Double I_NAI_L2x2y4z_D2y_a = I_NAI_M2x3y4z_Py_a+ABY*I_NAI_L2x2y4z_Py_a;
  Double I_NAI_L2xy5z_D2y_a = I_NAI_M2x2y5z_Py_a+ABY*I_NAI_L2xy5z_Py_a;
  Double I_NAI_L2x6z_D2y_a = I_NAI_M2xy6z_Py_a+ABY*I_NAI_L2x6z_Py_a;
  Double I_NAI_Lx7y_D2y_a = I_NAI_Mx8y_Py_a+ABY*I_NAI_Lx7y_Py_a;
  Double I_NAI_Lx6yz_D2y_a = I_NAI_Mx7yz_Py_a+ABY*I_NAI_Lx6yz_Py_a;
  Double I_NAI_Lx5y2z_D2y_a = I_NAI_Mx6y2z_Py_a+ABY*I_NAI_Lx5y2z_Py_a;
  Double I_NAI_Lx4y3z_D2y_a = I_NAI_Mx5y3z_Py_a+ABY*I_NAI_Lx4y3z_Py_a;
  Double I_NAI_Lx3y4z_D2y_a = I_NAI_Mx4y4z_Py_a+ABY*I_NAI_Lx3y4z_Py_a;
  Double I_NAI_Lx2y5z_D2y_a = I_NAI_Mx3y5z_Py_a+ABY*I_NAI_Lx2y5z_Py_a;
  Double I_NAI_Lxy6z_D2y_a = I_NAI_Mx2y6z_Py_a+ABY*I_NAI_Lxy6z_Py_a;
  Double I_NAI_Lx7z_D2y_a = I_NAI_Mxy7z_Py_a+ABY*I_NAI_Lx7z_Py_a;
  Double I_NAI_L8y_D2y_a = I_NAI_M9y_Py_a+ABY*I_NAI_L8y_Py_a;
  Double I_NAI_L7yz_D2y_a = I_NAI_M8yz_Py_a+ABY*I_NAI_L7yz_Py_a;
  Double I_NAI_L6y2z_D2y_a = I_NAI_M7y2z_Py_a+ABY*I_NAI_L6y2z_Py_a;
  Double I_NAI_L5y3z_D2y_a = I_NAI_M6y3z_Py_a+ABY*I_NAI_L5y3z_Py_a;
  Double I_NAI_L4y4z_D2y_a = I_NAI_M5y4z_Py_a+ABY*I_NAI_L4y4z_Py_a;
  Double I_NAI_L3y5z_D2y_a = I_NAI_M4y5z_Py_a+ABY*I_NAI_L3y5z_Py_a;
  Double I_NAI_L2y6z_D2y_a = I_NAI_M3y6z_Py_a+ABY*I_NAI_L2y6z_Py_a;
  Double I_NAI_Ly7z_D2y_a = I_NAI_M2y7z_Py_a+ABY*I_NAI_Ly7z_Py_a;
  Double I_NAI_L8z_D2y_a = I_NAI_My8z_Py_a+ABY*I_NAI_L8z_Py_a;
  Double I_NAI_L7xy_D2z_a = I_NAI_M7xyz_Pz_a+ABZ*I_NAI_L7xy_Pz_a;
  Double I_NAI_L7xz_D2z_a = I_NAI_M7x2z_Pz_a+ABZ*I_NAI_L7xz_Pz_a;
  Double I_NAI_L6x2y_D2z_a = I_NAI_M6x2yz_Pz_a+ABZ*I_NAI_L6x2y_Pz_a;
  Double I_NAI_L6xyz_D2z_a = I_NAI_M6xy2z_Pz_a+ABZ*I_NAI_L6xyz_Pz_a;
  Double I_NAI_L6x2z_D2z_a = I_NAI_M6x3z_Pz_a+ABZ*I_NAI_L6x2z_Pz_a;
  Double I_NAI_L5x3y_D2z_a = I_NAI_M5x3yz_Pz_a+ABZ*I_NAI_L5x3y_Pz_a;
  Double I_NAI_L5x2yz_D2z_a = I_NAI_M5x2y2z_Pz_a+ABZ*I_NAI_L5x2yz_Pz_a;
  Double I_NAI_L5xy2z_D2z_a = I_NAI_M5xy3z_Pz_a+ABZ*I_NAI_L5xy2z_Pz_a;
  Double I_NAI_L5x3z_D2z_a = I_NAI_M5x4z_Pz_a+ABZ*I_NAI_L5x3z_Pz_a;
  Double I_NAI_L4x4y_D2z_a = I_NAI_M4x4yz_Pz_a+ABZ*I_NAI_L4x4y_Pz_a;
  Double I_NAI_L4x3yz_D2z_a = I_NAI_M4x3y2z_Pz_a+ABZ*I_NAI_L4x3yz_Pz_a;
  Double I_NAI_L4x2y2z_D2z_a = I_NAI_M4x2y3z_Pz_a+ABZ*I_NAI_L4x2y2z_Pz_a;
  Double I_NAI_L4xy3z_D2z_a = I_NAI_M4xy4z_Pz_a+ABZ*I_NAI_L4xy3z_Pz_a;
  Double I_NAI_L4x4z_D2z_a = I_NAI_M4x5z_Pz_a+ABZ*I_NAI_L4x4z_Pz_a;
  Double I_NAI_L3x5y_D2z_a = I_NAI_M3x5yz_Pz_a+ABZ*I_NAI_L3x5y_Pz_a;
  Double I_NAI_L3x4yz_D2z_a = I_NAI_M3x4y2z_Pz_a+ABZ*I_NAI_L3x4yz_Pz_a;
  Double I_NAI_L3x3y2z_D2z_a = I_NAI_M3x3y3z_Pz_a+ABZ*I_NAI_L3x3y2z_Pz_a;
  Double I_NAI_L3x2y3z_D2z_a = I_NAI_M3x2y4z_Pz_a+ABZ*I_NAI_L3x2y3z_Pz_a;
  Double I_NAI_L3xy4z_D2z_a = I_NAI_M3xy5z_Pz_a+ABZ*I_NAI_L3xy4z_Pz_a;
  Double I_NAI_L3x5z_D2z_a = I_NAI_M3x6z_Pz_a+ABZ*I_NAI_L3x5z_Pz_a;
  Double I_NAI_L2x6y_D2z_a = I_NAI_M2x6yz_Pz_a+ABZ*I_NAI_L2x6y_Pz_a;
  Double I_NAI_L2x5yz_D2z_a = I_NAI_M2x5y2z_Pz_a+ABZ*I_NAI_L2x5yz_Pz_a;
  Double I_NAI_L2x4y2z_D2z_a = I_NAI_M2x4y3z_Pz_a+ABZ*I_NAI_L2x4y2z_Pz_a;
  Double I_NAI_L2x3y3z_D2z_a = I_NAI_M2x3y4z_Pz_a+ABZ*I_NAI_L2x3y3z_Pz_a;
  Double I_NAI_L2x2y4z_D2z_a = I_NAI_M2x2y5z_Pz_a+ABZ*I_NAI_L2x2y4z_Pz_a;
  Double I_NAI_L2xy5z_D2z_a = I_NAI_M2xy6z_Pz_a+ABZ*I_NAI_L2xy5z_Pz_a;
  Double I_NAI_L2x6z_D2z_a = I_NAI_M2x7z_Pz_a+ABZ*I_NAI_L2x6z_Pz_a;
  Double I_NAI_Lx7y_D2z_a = I_NAI_Mx7yz_Pz_a+ABZ*I_NAI_Lx7y_Pz_a;
  Double I_NAI_Lx6yz_D2z_a = I_NAI_Mx6y2z_Pz_a+ABZ*I_NAI_Lx6yz_Pz_a;
  Double I_NAI_Lx5y2z_D2z_a = I_NAI_Mx5y3z_Pz_a+ABZ*I_NAI_Lx5y2z_Pz_a;
  Double I_NAI_Lx4y3z_D2z_a = I_NAI_Mx4y4z_Pz_a+ABZ*I_NAI_Lx4y3z_Pz_a;
  Double I_NAI_Lx3y4z_D2z_a = I_NAI_Mx3y5z_Pz_a+ABZ*I_NAI_Lx3y4z_Pz_a;
  Double I_NAI_Lx2y5z_D2z_a = I_NAI_Mx2y6z_Pz_a+ABZ*I_NAI_Lx2y5z_Pz_a;
  Double I_NAI_Lxy6z_D2z_a = I_NAI_Mxy7z_Pz_a+ABZ*I_NAI_Lxy6z_Pz_a;
  Double I_NAI_Lx7z_D2z_a = I_NAI_Mx8z_Pz_a+ABZ*I_NAI_Lx7z_Pz_a;
  Double I_NAI_L7yz_D2z_a = I_NAI_M7y2z_Pz_a+ABZ*I_NAI_L7yz_Pz_a;
  Double I_NAI_L6y2z_D2z_a = I_NAI_M6y3z_Pz_a+ABZ*I_NAI_L6y2z_Pz_a;
  Double I_NAI_L5y3z_D2z_a = I_NAI_M5y4z_Pz_a+ABZ*I_NAI_L5y3z_Pz_a;
  Double I_NAI_L4y4z_D2z_a = I_NAI_M4y5z_Pz_a+ABZ*I_NAI_L4y4z_Pz_a;
  Double I_NAI_L3y5z_D2z_a = I_NAI_M3y6z_Pz_a+ABZ*I_NAI_L3y5z_Pz_a;
  Double I_NAI_L2y6z_D2z_a = I_NAI_M2y7z_Pz_a+ABZ*I_NAI_L2y6z_Pz_a;
  Double I_NAI_Ly7z_D2z_a = I_NAI_My8z_Pz_a+ABZ*I_NAI_Ly7z_Pz_a;
  Double I_NAI_L8z_D2z_a = I_NAI_M9z_Pz_a+ABZ*I_NAI_L8z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_K_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 105 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_D_a
   * RHS shell quartet name: SQ_NAI_K_D_a
   ************************************************************/
  Double I_NAI_K7x_F3x_a = I_NAI_L8x_D2x_a+ABX*I_NAI_K7x_D2x_a;
  Double I_NAI_K6xy_F3x_a = I_NAI_L7xy_D2x_a+ABX*I_NAI_K6xy_D2x_a;
  Double I_NAI_K6xz_F3x_a = I_NAI_L7xz_D2x_a+ABX*I_NAI_K6xz_D2x_a;
  Double I_NAI_K5x2y_F3x_a = I_NAI_L6x2y_D2x_a+ABX*I_NAI_K5x2y_D2x_a;
  Double I_NAI_K5xyz_F3x_a = I_NAI_L6xyz_D2x_a+ABX*I_NAI_K5xyz_D2x_a;
  Double I_NAI_K5x2z_F3x_a = I_NAI_L6x2z_D2x_a+ABX*I_NAI_K5x2z_D2x_a;
  Double I_NAI_K4x3y_F3x_a = I_NAI_L5x3y_D2x_a+ABX*I_NAI_K4x3y_D2x_a;
  Double I_NAI_K4x2yz_F3x_a = I_NAI_L5x2yz_D2x_a+ABX*I_NAI_K4x2yz_D2x_a;
  Double I_NAI_K4xy2z_F3x_a = I_NAI_L5xy2z_D2x_a+ABX*I_NAI_K4xy2z_D2x_a;
  Double I_NAI_K4x3z_F3x_a = I_NAI_L5x3z_D2x_a+ABX*I_NAI_K4x3z_D2x_a;
  Double I_NAI_K3x4y_F3x_a = I_NAI_L4x4y_D2x_a+ABX*I_NAI_K3x4y_D2x_a;
  Double I_NAI_K3x3yz_F3x_a = I_NAI_L4x3yz_D2x_a+ABX*I_NAI_K3x3yz_D2x_a;
  Double I_NAI_K3x2y2z_F3x_a = I_NAI_L4x2y2z_D2x_a+ABX*I_NAI_K3x2y2z_D2x_a;
  Double I_NAI_K3xy3z_F3x_a = I_NAI_L4xy3z_D2x_a+ABX*I_NAI_K3xy3z_D2x_a;
  Double I_NAI_K3x4z_F3x_a = I_NAI_L4x4z_D2x_a+ABX*I_NAI_K3x4z_D2x_a;
  Double I_NAI_K2x5y_F3x_a = I_NAI_L3x5y_D2x_a+ABX*I_NAI_K2x5y_D2x_a;
  Double I_NAI_K2x4yz_F3x_a = I_NAI_L3x4yz_D2x_a+ABX*I_NAI_K2x4yz_D2x_a;
  Double I_NAI_K2x3y2z_F3x_a = I_NAI_L3x3y2z_D2x_a+ABX*I_NAI_K2x3y2z_D2x_a;
  Double I_NAI_K2x2y3z_F3x_a = I_NAI_L3x2y3z_D2x_a+ABX*I_NAI_K2x2y3z_D2x_a;
  Double I_NAI_K2xy4z_F3x_a = I_NAI_L3xy4z_D2x_a+ABX*I_NAI_K2xy4z_D2x_a;
  Double I_NAI_K2x5z_F3x_a = I_NAI_L3x5z_D2x_a+ABX*I_NAI_K2x5z_D2x_a;
  Double I_NAI_Kx6y_F3x_a = I_NAI_L2x6y_D2x_a+ABX*I_NAI_Kx6y_D2x_a;
  Double I_NAI_Kx5yz_F3x_a = I_NAI_L2x5yz_D2x_a+ABX*I_NAI_Kx5yz_D2x_a;
  Double I_NAI_Kx4y2z_F3x_a = I_NAI_L2x4y2z_D2x_a+ABX*I_NAI_Kx4y2z_D2x_a;
  Double I_NAI_Kx3y3z_F3x_a = I_NAI_L2x3y3z_D2x_a+ABX*I_NAI_Kx3y3z_D2x_a;
  Double I_NAI_Kx2y4z_F3x_a = I_NAI_L2x2y4z_D2x_a+ABX*I_NAI_Kx2y4z_D2x_a;
  Double I_NAI_Kxy5z_F3x_a = I_NAI_L2xy5z_D2x_a+ABX*I_NAI_Kxy5z_D2x_a;
  Double I_NAI_Kx6z_F3x_a = I_NAI_L2x6z_D2x_a+ABX*I_NAI_Kx6z_D2x_a;
  Double I_NAI_K7y_F3x_a = I_NAI_Lx7y_D2x_a+ABX*I_NAI_K7y_D2x_a;
  Double I_NAI_K6yz_F3x_a = I_NAI_Lx6yz_D2x_a+ABX*I_NAI_K6yz_D2x_a;
  Double I_NAI_K5y2z_F3x_a = I_NAI_Lx5y2z_D2x_a+ABX*I_NAI_K5y2z_D2x_a;
  Double I_NAI_K4y3z_F3x_a = I_NAI_Lx4y3z_D2x_a+ABX*I_NAI_K4y3z_D2x_a;
  Double I_NAI_K3y4z_F3x_a = I_NAI_Lx3y4z_D2x_a+ABX*I_NAI_K3y4z_D2x_a;
  Double I_NAI_K2y5z_F3x_a = I_NAI_Lx2y5z_D2x_a+ABX*I_NAI_K2y5z_D2x_a;
  Double I_NAI_Ky6z_F3x_a = I_NAI_Lxy6z_D2x_a+ABX*I_NAI_Ky6z_D2x_a;
  Double I_NAI_K7z_F3x_a = I_NAI_Lx7z_D2x_a+ABX*I_NAI_K7z_D2x_a;
  Double I_NAI_K6xy_F2xy_a = I_NAI_L6x2y_D2x_a+ABY*I_NAI_K6xy_D2x_a;
  Double I_NAI_K6xz_F2xy_a = I_NAI_L6xyz_D2x_a+ABY*I_NAI_K6xz_D2x_a;
  Double I_NAI_K5x2y_F2xy_a = I_NAI_L5x3y_D2x_a+ABY*I_NAI_K5x2y_D2x_a;
  Double I_NAI_K5xyz_F2xy_a = I_NAI_L5x2yz_D2x_a+ABY*I_NAI_K5xyz_D2x_a;
  Double I_NAI_K5x2z_F2xy_a = I_NAI_L5xy2z_D2x_a+ABY*I_NAI_K5x2z_D2x_a;
  Double I_NAI_K4x3y_F2xy_a = I_NAI_L4x4y_D2x_a+ABY*I_NAI_K4x3y_D2x_a;
  Double I_NAI_K4x2yz_F2xy_a = I_NAI_L4x3yz_D2x_a+ABY*I_NAI_K4x2yz_D2x_a;
  Double I_NAI_K4xy2z_F2xy_a = I_NAI_L4x2y2z_D2x_a+ABY*I_NAI_K4xy2z_D2x_a;
  Double I_NAI_K4x3z_F2xy_a = I_NAI_L4xy3z_D2x_a+ABY*I_NAI_K4x3z_D2x_a;
  Double I_NAI_K3x4y_F2xy_a = I_NAI_L3x5y_D2x_a+ABY*I_NAI_K3x4y_D2x_a;
  Double I_NAI_K3x3yz_F2xy_a = I_NAI_L3x4yz_D2x_a+ABY*I_NAI_K3x3yz_D2x_a;
  Double I_NAI_K3x2y2z_F2xy_a = I_NAI_L3x3y2z_D2x_a+ABY*I_NAI_K3x2y2z_D2x_a;
  Double I_NAI_K3xy3z_F2xy_a = I_NAI_L3x2y3z_D2x_a+ABY*I_NAI_K3xy3z_D2x_a;
  Double I_NAI_K3x4z_F2xy_a = I_NAI_L3xy4z_D2x_a+ABY*I_NAI_K3x4z_D2x_a;
  Double I_NAI_K2x5y_F2xy_a = I_NAI_L2x6y_D2x_a+ABY*I_NAI_K2x5y_D2x_a;
  Double I_NAI_K2x4yz_F2xy_a = I_NAI_L2x5yz_D2x_a+ABY*I_NAI_K2x4yz_D2x_a;
  Double I_NAI_K2x3y2z_F2xy_a = I_NAI_L2x4y2z_D2x_a+ABY*I_NAI_K2x3y2z_D2x_a;
  Double I_NAI_K2x2y3z_F2xy_a = I_NAI_L2x3y3z_D2x_a+ABY*I_NAI_K2x2y3z_D2x_a;
  Double I_NAI_K2xy4z_F2xy_a = I_NAI_L2x2y4z_D2x_a+ABY*I_NAI_K2xy4z_D2x_a;
  Double I_NAI_K2x5z_F2xy_a = I_NAI_L2xy5z_D2x_a+ABY*I_NAI_K2x5z_D2x_a;
  Double I_NAI_Kx6y_F2xy_a = I_NAI_Lx7y_D2x_a+ABY*I_NAI_Kx6y_D2x_a;
  Double I_NAI_Kx5yz_F2xy_a = I_NAI_Lx6yz_D2x_a+ABY*I_NAI_Kx5yz_D2x_a;
  Double I_NAI_Kx4y2z_F2xy_a = I_NAI_Lx5y2z_D2x_a+ABY*I_NAI_Kx4y2z_D2x_a;
  Double I_NAI_Kx3y3z_F2xy_a = I_NAI_Lx4y3z_D2x_a+ABY*I_NAI_Kx3y3z_D2x_a;
  Double I_NAI_Kx2y4z_F2xy_a = I_NAI_Lx3y4z_D2x_a+ABY*I_NAI_Kx2y4z_D2x_a;
  Double I_NAI_Kxy5z_F2xy_a = I_NAI_Lx2y5z_D2x_a+ABY*I_NAI_Kxy5z_D2x_a;
  Double I_NAI_Kx6z_F2xy_a = I_NAI_Lxy6z_D2x_a+ABY*I_NAI_Kx6z_D2x_a;
  Double I_NAI_K7y_F2xy_a = I_NAI_L8y_D2x_a+ABY*I_NAI_K7y_D2x_a;
  Double I_NAI_K6yz_F2xy_a = I_NAI_L7yz_D2x_a+ABY*I_NAI_K6yz_D2x_a;
  Double I_NAI_K5y2z_F2xy_a = I_NAI_L6y2z_D2x_a+ABY*I_NAI_K5y2z_D2x_a;
  Double I_NAI_K4y3z_F2xy_a = I_NAI_L5y3z_D2x_a+ABY*I_NAI_K4y3z_D2x_a;
  Double I_NAI_K3y4z_F2xy_a = I_NAI_L4y4z_D2x_a+ABY*I_NAI_K3y4z_D2x_a;
  Double I_NAI_K2y5z_F2xy_a = I_NAI_L3y5z_D2x_a+ABY*I_NAI_K2y5z_D2x_a;
  Double I_NAI_Ky6z_F2xy_a = I_NAI_L2y6z_D2x_a+ABY*I_NAI_Ky6z_D2x_a;
  Double I_NAI_K7z_F2xy_a = I_NAI_Ly7z_D2x_a+ABY*I_NAI_K7z_D2x_a;
  Double I_NAI_K6xz_F2xz_a = I_NAI_L6x2z_D2x_a+ABZ*I_NAI_K6xz_D2x_a;
  Double I_NAI_K5xyz_F2xz_a = I_NAI_L5xy2z_D2x_a+ABZ*I_NAI_K5xyz_D2x_a;
  Double I_NAI_K5x2z_F2xz_a = I_NAI_L5x3z_D2x_a+ABZ*I_NAI_K5x2z_D2x_a;
  Double I_NAI_K4x2yz_F2xz_a = I_NAI_L4x2y2z_D2x_a+ABZ*I_NAI_K4x2yz_D2x_a;
  Double I_NAI_K4xy2z_F2xz_a = I_NAI_L4xy3z_D2x_a+ABZ*I_NAI_K4xy2z_D2x_a;
  Double I_NAI_K4x3z_F2xz_a = I_NAI_L4x4z_D2x_a+ABZ*I_NAI_K4x3z_D2x_a;
  Double I_NAI_K3x3yz_F2xz_a = I_NAI_L3x3y2z_D2x_a+ABZ*I_NAI_K3x3yz_D2x_a;
  Double I_NAI_K3x2y2z_F2xz_a = I_NAI_L3x2y3z_D2x_a+ABZ*I_NAI_K3x2y2z_D2x_a;
  Double I_NAI_K3xy3z_F2xz_a = I_NAI_L3xy4z_D2x_a+ABZ*I_NAI_K3xy3z_D2x_a;
  Double I_NAI_K3x4z_F2xz_a = I_NAI_L3x5z_D2x_a+ABZ*I_NAI_K3x4z_D2x_a;
  Double I_NAI_K2x4yz_F2xz_a = I_NAI_L2x4y2z_D2x_a+ABZ*I_NAI_K2x4yz_D2x_a;
  Double I_NAI_K2x3y2z_F2xz_a = I_NAI_L2x3y3z_D2x_a+ABZ*I_NAI_K2x3y2z_D2x_a;
  Double I_NAI_K2x2y3z_F2xz_a = I_NAI_L2x2y4z_D2x_a+ABZ*I_NAI_K2x2y3z_D2x_a;
  Double I_NAI_K2xy4z_F2xz_a = I_NAI_L2xy5z_D2x_a+ABZ*I_NAI_K2xy4z_D2x_a;
  Double I_NAI_K2x5z_F2xz_a = I_NAI_L2x6z_D2x_a+ABZ*I_NAI_K2x5z_D2x_a;
  Double I_NAI_Kx5yz_F2xz_a = I_NAI_Lx5y2z_D2x_a+ABZ*I_NAI_Kx5yz_D2x_a;
  Double I_NAI_Kx4y2z_F2xz_a = I_NAI_Lx4y3z_D2x_a+ABZ*I_NAI_Kx4y2z_D2x_a;
  Double I_NAI_Kx3y3z_F2xz_a = I_NAI_Lx3y4z_D2x_a+ABZ*I_NAI_Kx3y3z_D2x_a;
  Double I_NAI_Kx2y4z_F2xz_a = I_NAI_Lx2y5z_D2x_a+ABZ*I_NAI_Kx2y4z_D2x_a;
  Double I_NAI_Kxy5z_F2xz_a = I_NAI_Lxy6z_D2x_a+ABZ*I_NAI_Kxy5z_D2x_a;
  Double I_NAI_Kx6z_F2xz_a = I_NAI_Lx7z_D2x_a+ABZ*I_NAI_Kx6z_D2x_a;
  Double I_NAI_K6yz_F2xz_a = I_NAI_L6y2z_D2x_a+ABZ*I_NAI_K6yz_D2x_a;
  Double I_NAI_K5y2z_F2xz_a = I_NAI_L5y3z_D2x_a+ABZ*I_NAI_K5y2z_D2x_a;
  Double I_NAI_K4y3z_F2xz_a = I_NAI_L4y4z_D2x_a+ABZ*I_NAI_K4y3z_D2x_a;
  Double I_NAI_K3y4z_F2xz_a = I_NAI_L3y5z_D2x_a+ABZ*I_NAI_K3y4z_D2x_a;
  Double I_NAI_K2y5z_F2xz_a = I_NAI_L2y6z_D2x_a+ABZ*I_NAI_K2y5z_D2x_a;
  Double I_NAI_Ky6z_F2xz_a = I_NAI_Ly7z_D2x_a+ABZ*I_NAI_Ky6z_D2x_a;
  Double I_NAI_K7z_F2xz_a = I_NAI_L8z_D2x_a+ABZ*I_NAI_K7z_D2x_a;
  Double I_NAI_K6xz_Fx2y_a = I_NAI_L7xz_D2y_a+ABX*I_NAI_K6xz_D2y_a;
  Double I_NAI_K5xyz_Fx2y_a = I_NAI_L6xyz_D2y_a+ABX*I_NAI_K5xyz_D2y_a;
  Double I_NAI_K5x2z_Fx2y_a = I_NAI_L6x2z_D2y_a+ABX*I_NAI_K5x2z_D2y_a;
  Double I_NAI_K4x2yz_Fx2y_a = I_NAI_L5x2yz_D2y_a+ABX*I_NAI_K4x2yz_D2y_a;
  Double I_NAI_K4xy2z_Fx2y_a = I_NAI_L5xy2z_D2y_a+ABX*I_NAI_K4xy2z_D2y_a;
  Double I_NAI_K4x3z_Fx2y_a = I_NAI_L5x3z_D2y_a+ABX*I_NAI_K4x3z_D2y_a;
  Double I_NAI_K3x3yz_Fx2y_a = I_NAI_L4x3yz_D2y_a+ABX*I_NAI_K3x3yz_D2y_a;
  Double I_NAI_K3x2y2z_Fx2y_a = I_NAI_L4x2y2z_D2y_a+ABX*I_NAI_K3x2y2z_D2y_a;
  Double I_NAI_K3xy3z_Fx2y_a = I_NAI_L4xy3z_D2y_a+ABX*I_NAI_K3xy3z_D2y_a;
  Double I_NAI_K3x4z_Fx2y_a = I_NAI_L4x4z_D2y_a+ABX*I_NAI_K3x4z_D2y_a;
  Double I_NAI_K2x4yz_Fx2y_a = I_NAI_L3x4yz_D2y_a+ABX*I_NAI_K2x4yz_D2y_a;
  Double I_NAI_K2x3y2z_Fx2y_a = I_NAI_L3x3y2z_D2y_a+ABX*I_NAI_K2x3y2z_D2y_a;
  Double I_NAI_K2x2y3z_Fx2y_a = I_NAI_L3x2y3z_D2y_a+ABX*I_NAI_K2x2y3z_D2y_a;
  Double I_NAI_K2xy4z_Fx2y_a = I_NAI_L3xy4z_D2y_a+ABX*I_NAI_K2xy4z_D2y_a;
  Double I_NAI_K2x5z_Fx2y_a = I_NAI_L3x5z_D2y_a+ABX*I_NAI_K2x5z_D2y_a;
  Double I_NAI_Kx5yz_Fx2y_a = I_NAI_L2x5yz_D2y_a+ABX*I_NAI_Kx5yz_D2y_a;
  Double I_NAI_Kx4y2z_Fx2y_a = I_NAI_L2x4y2z_D2y_a+ABX*I_NAI_Kx4y2z_D2y_a;
  Double I_NAI_Kx3y3z_Fx2y_a = I_NAI_L2x3y3z_D2y_a+ABX*I_NAI_Kx3y3z_D2y_a;
  Double I_NAI_Kx2y4z_Fx2y_a = I_NAI_L2x2y4z_D2y_a+ABX*I_NAI_Kx2y4z_D2y_a;
  Double I_NAI_Kxy5z_Fx2y_a = I_NAI_L2xy5z_D2y_a+ABX*I_NAI_Kxy5z_D2y_a;
  Double I_NAI_Kx6z_Fx2y_a = I_NAI_L2x6z_D2y_a+ABX*I_NAI_Kx6z_D2y_a;
  Double I_NAI_K6yz_Fx2y_a = I_NAI_Lx6yz_D2y_a+ABX*I_NAI_K6yz_D2y_a;
  Double I_NAI_K5y2z_Fx2y_a = I_NAI_Lx5y2z_D2y_a+ABX*I_NAI_K5y2z_D2y_a;
  Double I_NAI_K4y3z_Fx2y_a = I_NAI_Lx4y3z_D2y_a+ABX*I_NAI_K4y3z_D2y_a;
  Double I_NAI_K3y4z_Fx2y_a = I_NAI_Lx3y4z_D2y_a+ABX*I_NAI_K3y4z_D2y_a;
  Double I_NAI_K2y5z_Fx2y_a = I_NAI_Lx2y5z_D2y_a+ABX*I_NAI_K2y5z_D2y_a;
  Double I_NAI_Ky6z_Fx2y_a = I_NAI_Lxy6z_D2y_a+ABX*I_NAI_Ky6z_D2y_a;
  Double I_NAI_K7z_Fx2y_a = I_NAI_Lx7z_D2y_a+ABX*I_NAI_K7z_D2y_a;
  Double I_NAI_K6xy_Fx2z_a = I_NAI_L7xy_D2z_a+ABX*I_NAI_K6xy_D2z_a;
  Double I_NAI_K5x2y_Fx2z_a = I_NAI_L6x2y_D2z_a+ABX*I_NAI_K5x2y_D2z_a;
  Double I_NAI_K5xyz_Fx2z_a = I_NAI_L6xyz_D2z_a+ABX*I_NAI_K5xyz_D2z_a;
  Double I_NAI_K4x3y_Fx2z_a = I_NAI_L5x3y_D2z_a+ABX*I_NAI_K4x3y_D2z_a;
  Double I_NAI_K4x2yz_Fx2z_a = I_NAI_L5x2yz_D2z_a+ABX*I_NAI_K4x2yz_D2z_a;
  Double I_NAI_K4xy2z_Fx2z_a = I_NAI_L5xy2z_D2z_a+ABX*I_NAI_K4xy2z_D2z_a;
  Double I_NAI_K3x4y_Fx2z_a = I_NAI_L4x4y_D2z_a+ABX*I_NAI_K3x4y_D2z_a;
  Double I_NAI_K3x3yz_Fx2z_a = I_NAI_L4x3yz_D2z_a+ABX*I_NAI_K3x3yz_D2z_a;
  Double I_NAI_K3x2y2z_Fx2z_a = I_NAI_L4x2y2z_D2z_a+ABX*I_NAI_K3x2y2z_D2z_a;
  Double I_NAI_K3xy3z_Fx2z_a = I_NAI_L4xy3z_D2z_a+ABX*I_NAI_K3xy3z_D2z_a;
  Double I_NAI_K2x5y_Fx2z_a = I_NAI_L3x5y_D2z_a+ABX*I_NAI_K2x5y_D2z_a;
  Double I_NAI_K2x4yz_Fx2z_a = I_NAI_L3x4yz_D2z_a+ABX*I_NAI_K2x4yz_D2z_a;
  Double I_NAI_K2x3y2z_Fx2z_a = I_NAI_L3x3y2z_D2z_a+ABX*I_NAI_K2x3y2z_D2z_a;
  Double I_NAI_K2x2y3z_Fx2z_a = I_NAI_L3x2y3z_D2z_a+ABX*I_NAI_K2x2y3z_D2z_a;
  Double I_NAI_K2xy4z_Fx2z_a = I_NAI_L3xy4z_D2z_a+ABX*I_NAI_K2xy4z_D2z_a;
  Double I_NAI_Kx6y_Fx2z_a = I_NAI_L2x6y_D2z_a+ABX*I_NAI_Kx6y_D2z_a;
  Double I_NAI_Kx5yz_Fx2z_a = I_NAI_L2x5yz_D2z_a+ABX*I_NAI_Kx5yz_D2z_a;
  Double I_NAI_Kx4y2z_Fx2z_a = I_NAI_L2x4y2z_D2z_a+ABX*I_NAI_Kx4y2z_D2z_a;
  Double I_NAI_Kx3y3z_Fx2z_a = I_NAI_L2x3y3z_D2z_a+ABX*I_NAI_Kx3y3z_D2z_a;
  Double I_NAI_Kx2y4z_Fx2z_a = I_NAI_L2x2y4z_D2z_a+ABX*I_NAI_Kx2y4z_D2z_a;
  Double I_NAI_Kxy5z_Fx2z_a = I_NAI_L2xy5z_D2z_a+ABX*I_NAI_Kxy5z_D2z_a;
  Double I_NAI_K7y_Fx2z_a = I_NAI_Lx7y_D2z_a+ABX*I_NAI_K7y_D2z_a;
  Double I_NAI_K6yz_Fx2z_a = I_NAI_Lx6yz_D2z_a+ABX*I_NAI_K6yz_D2z_a;
  Double I_NAI_K5y2z_Fx2z_a = I_NAI_Lx5y2z_D2z_a+ABX*I_NAI_K5y2z_D2z_a;
  Double I_NAI_K4y3z_Fx2z_a = I_NAI_Lx4y3z_D2z_a+ABX*I_NAI_K4y3z_D2z_a;
  Double I_NAI_K3y4z_Fx2z_a = I_NAI_Lx3y4z_D2z_a+ABX*I_NAI_K3y4z_D2z_a;
  Double I_NAI_K2y5z_Fx2z_a = I_NAI_Lx2y5z_D2z_a+ABX*I_NAI_K2y5z_D2z_a;
  Double I_NAI_Ky6z_Fx2z_a = I_NAI_Lxy6z_D2z_a+ABX*I_NAI_Ky6z_D2z_a;
  Double I_NAI_K7x_F3y_a = I_NAI_L7xy_D2y_a+ABY*I_NAI_K7x_D2y_a;
  Double I_NAI_K6xy_F3y_a = I_NAI_L6x2y_D2y_a+ABY*I_NAI_K6xy_D2y_a;
  Double I_NAI_K6xz_F3y_a = I_NAI_L6xyz_D2y_a+ABY*I_NAI_K6xz_D2y_a;
  Double I_NAI_K5x2y_F3y_a = I_NAI_L5x3y_D2y_a+ABY*I_NAI_K5x2y_D2y_a;
  Double I_NAI_K5xyz_F3y_a = I_NAI_L5x2yz_D2y_a+ABY*I_NAI_K5xyz_D2y_a;
  Double I_NAI_K5x2z_F3y_a = I_NAI_L5xy2z_D2y_a+ABY*I_NAI_K5x2z_D2y_a;
  Double I_NAI_K4x3y_F3y_a = I_NAI_L4x4y_D2y_a+ABY*I_NAI_K4x3y_D2y_a;
  Double I_NAI_K4x2yz_F3y_a = I_NAI_L4x3yz_D2y_a+ABY*I_NAI_K4x2yz_D2y_a;
  Double I_NAI_K4xy2z_F3y_a = I_NAI_L4x2y2z_D2y_a+ABY*I_NAI_K4xy2z_D2y_a;
  Double I_NAI_K4x3z_F3y_a = I_NAI_L4xy3z_D2y_a+ABY*I_NAI_K4x3z_D2y_a;
  Double I_NAI_K3x4y_F3y_a = I_NAI_L3x5y_D2y_a+ABY*I_NAI_K3x4y_D2y_a;
  Double I_NAI_K3x3yz_F3y_a = I_NAI_L3x4yz_D2y_a+ABY*I_NAI_K3x3yz_D2y_a;
  Double I_NAI_K3x2y2z_F3y_a = I_NAI_L3x3y2z_D2y_a+ABY*I_NAI_K3x2y2z_D2y_a;
  Double I_NAI_K3xy3z_F3y_a = I_NAI_L3x2y3z_D2y_a+ABY*I_NAI_K3xy3z_D2y_a;
  Double I_NAI_K3x4z_F3y_a = I_NAI_L3xy4z_D2y_a+ABY*I_NAI_K3x4z_D2y_a;
  Double I_NAI_K2x5y_F3y_a = I_NAI_L2x6y_D2y_a+ABY*I_NAI_K2x5y_D2y_a;
  Double I_NAI_K2x4yz_F3y_a = I_NAI_L2x5yz_D2y_a+ABY*I_NAI_K2x4yz_D2y_a;
  Double I_NAI_K2x3y2z_F3y_a = I_NAI_L2x4y2z_D2y_a+ABY*I_NAI_K2x3y2z_D2y_a;
  Double I_NAI_K2x2y3z_F3y_a = I_NAI_L2x3y3z_D2y_a+ABY*I_NAI_K2x2y3z_D2y_a;
  Double I_NAI_K2xy4z_F3y_a = I_NAI_L2x2y4z_D2y_a+ABY*I_NAI_K2xy4z_D2y_a;
  Double I_NAI_K2x5z_F3y_a = I_NAI_L2xy5z_D2y_a+ABY*I_NAI_K2x5z_D2y_a;
  Double I_NAI_Kx6y_F3y_a = I_NAI_Lx7y_D2y_a+ABY*I_NAI_Kx6y_D2y_a;
  Double I_NAI_Kx5yz_F3y_a = I_NAI_Lx6yz_D2y_a+ABY*I_NAI_Kx5yz_D2y_a;
  Double I_NAI_Kx4y2z_F3y_a = I_NAI_Lx5y2z_D2y_a+ABY*I_NAI_Kx4y2z_D2y_a;
  Double I_NAI_Kx3y3z_F3y_a = I_NAI_Lx4y3z_D2y_a+ABY*I_NAI_Kx3y3z_D2y_a;
  Double I_NAI_Kx2y4z_F3y_a = I_NAI_Lx3y4z_D2y_a+ABY*I_NAI_Kx2y4z_D2y_a;
  Double I_NAI_Kxy5z_F3y_a = I_NAI_Lx2y5z_D2y_a+ABY*I_NAI_Kxy5z_D2y_a;
  Double I_NAI_Kx6z_F3y_a = I_NAI_Lxy6z_D2y_a+ABY*I_NAI_Kx6z_D2y_a;
  Double I_NAI_K7y_F3y_a = I_NAI_L8y_D2y_a+ABY*I_NAI_K7y_D2y_a;
  Double I_NAI_K6yz_F3y_a = I_NAI_L7yz_D2y_a+ABY*I_NAI_K6yz_D2y_a;
  Double I_NAI_K5y2z_F3y_a = I_NAI_L6y2z_D2y_a+ABY*I_NAI_K5y2z_D2y_a;
  Double I_NAI_K4y3z_F3y_a = I_NAI_L5y3z_D2y_a+ABY*I_NAI_K4y3z_D2y_a;
  Double I_NAI_K3y4z_F3y_a = I_NAI_L4y4z_D2y_a+ABY*I_NAI_K3y4z_D2y_a;
  Double I_NAI_K2y5z_F3y_a = I_NAI_L3y5z_D2y_a+ABY*I_NAI_K2y5z_D2y_a;
  Double I_NAI_Ky6z_F3y_a = I_NAI_L2y6z_D2y_a+ABY*I_NAI_Ky6z_D2y_a;
  Double I_NAI_K7z_F3y_a = I_NAI_Ly7z_D2y_a+ABY*I_NAI_K7z_D2y_a;
  Double I_NAI_K6xz_F2yz_a = I_NAI_L6x2z_D2y_a+ABZ*I_NAI_K6xz_D2y_a;
  Double I_NAI_K5xyz_F2yz_a = I_NAI_L5xy2z_D2y_a+ABZ*I_NAI_K5xyz_D2y_a;
  Double I_NAI_K5x2z_F2yz_a = I_NAI_L5x3z_D2y_a+ABZ*I_NAI_K5x2z_D2y_a;
  Double I_NAI_K4x2yz_F2yz_a = I_NAI_L4x2y2z_D2y_a+ABZ*I_NAI_K4x2yz_D2y_a;
  Double I_NAI_K4xy2z_F2yz_a = I_NAI_L4xy3z_D2y_a+ABZ*I_NAI_K4xy2z_D2y_a;
  Double I_NAI_K4x3z_F2yz_a = I_NAI_L4x4z_D2y_a+ABZ*I_NAI_K4x3z_D2y_a;
  Double I_NAI_K3x3yz_F2yz_a = I_NAI_L3x3y2z_D2y_a+ABZ*I_NAI_K3x3yz_D2y_a;
  Double I_NAI_K3x2y2z_F2yz_a = I_NAI_L3x2y3z_D2y_a+ABZ*I_NAI_K3x2y2z_D2y_a;
  Double I_NAI_K3xy3z_F2yz_a = I_NAI_L3xy4z_D2y_a+ABZ*I_NAI_K3xy3z_D2y_a;
  Double I_NAI_K3x4z_F2yz_a = I_NAI_L3x5z_D2y_a+ABZ*I_NAI_K3x4z_D2y_a;
  Double I_NAI_K2x4yz_F2yz_a = I_NAI_L2x4y2z_D2y_a+ABZ*I_NAI_K2x4yz_D2y_a;
  Double I_NAI_K2x3y2z_F2yz_a = I_NAI_L2x3y3z_D2y_a+ABZ*I_NAI_K2x3y2z_D2y_a;
  Double I_NAI_K2x2y3z_F2yz_a = I_NAI_L2x2y4z_D2y_a+ABZ*I_NAI_K2x2y3z_D2y_a;
  Double I_NAI_K2xy4z_F2yz_a = I_NAI_L2xy5z_D2y_a+ABZ*I_NAI_K2xy4z_D2y_a;
  Double I_NAI_K2x5z_F2yz_a = I_NAI_L2x6z_D2y_a+ABZ*I_NAI_K2x5z_D2y_a;
  Double I_NAI_Kx5yz_F2yz_a = I_NAI_Lx5y2z_D2y_a+ABZ*I_NAI_Kx5yz_D2y_a;
  Double I_NAI_Kx4y2z_F2yz_a = I_NAI_Lx4y3z_D2y_a+ABZ*I_NAI_Kx4y2z_D2y_a;
  Double I_NAI_Kx3y3z_F2yz_a = I_NAI_Lx3y4z_D2y_a+ABZ*I_NAI_Kx3y3z_D2y_a;
  Double I_NAI_Kx2y4z_F2yz_a = I_NAI_Lx2y5z_D2y_a+ABZ*I_NAI_Kx2y4z_D2y_a;
  Double I_NAI_Kxy5z_F2yz_a = I_NAI_Lxy6z_D2y_a+ABZ*I_NAI_Kxy5z_D2y_a;
  Double I_NAI_Kx6z_F2yz_a = I_NAI_Lx7z_D2y_a+ABZ*I_NAI_Kx6z_D2y_a;
  Double I_NAI_K6yz_F2yz_a = I_NAI_L6y2z_D2y_a+ABZ*I_NAI_K6yz_D2y_a;
  Double I_NAI_K5y2z_F2yz_a = I_NAI_L5y3z_D2y_a+ABZ*I_NAI_K5y2z_D2y_a;
  Double I_NAI_K4y3z_F2yz_a = I_NAI_L4y4z_D2y_a+ABZ*I_NAI_K4y3z_D2y_a;
  Double I_NAI_K3y4z_F2yz_a = I_NAI_L3y5z_D2y_a+ABZ*I_NAI_K3y4z_D2y_a;
  Double I_NAI_K2y5z_F2yz_a = I_NAI_L2y6z_D2y_a+ABZ*I_NAI_K2y5z_D2y_a;
  Double I_NAI_Ky6z_F2yz_a = I_NAI_Ly7z_D2y_a+ABZ*I_NAI_Ky6z_D2y_a;
  Double I_NAI_K7z_F2yz_a = I_NAI_L8z_D2y_a+ABZ*I_NAI_K7z_D2y_a;
  Double I_NAI_K7x_F3z_a = I_NAI_L7xz_D2z_a+ABZ*I_NAI_K7x_D2z_a;
  Double I_NAI_K6xy_F3z_a = I_NAI_L6xyz_D2z_a+ABZ*I_NAI_K6xy_D2z_a;
  Double I_NAI_K6xz_F3z_a = I_NAI_L6x2z_D2z_a+ABZ*I_NAI_K6xz_D2z_a;
  Double I_NAI_K5x2y_F3z_a = I_NAI_L5x2yz_D2z_a+ABZ*I_NAI_K5x2y_D2z_a;
  Double I_NAI_K5xyz_F3z_a = I_NAI_L5xy2z_D2z_a+ABZ*I_NAI_K5xyz_D2z_a;
  Double I_NAI_K5x2z_F3z_a = I_NAI_L5x3z_D2z_a+ABZ*I_NAI_K5x2z_D2z_a;
  Double I_NAI_K4x3y_F3z_a = I_NAI_L4x3yz_D2z_a+ABZ*I_NAI_K4x3y_D2z_a;
  Double I_NAI_K4x2yz_F3z_a = I_NAI_L4x2y2z_D2z_a+ABZ*I_NAI_K4x2yz_D2z_a;
  Double I_NAI_K4xy2z_F3z_a = I_NAI_L4xy3z_D2z_a+ABZ*I_NAI_K4xy2z_D2z_a;
  Double I_NAI_K4x3z_F3z_a = I_NAI_L4x4z_D2z_a+ABZ*I_NAI_K4x3z_D2z_a;
  Double I_NAI_K3x4y_F3z_a = I_NAI_L3x4yz_D2z_a+ABZ*I_NAI_K3x4y_D2z_a;
  Double I_NAI_K3x3yz_F3z_a = I_NAI_L3x3y2z_D2z_a+ABZ*I_NAI_K3x3yz_D2z_a;
  Double I_NAI_K3x2y2z_F3z_a = I_NAI_L3x2y3z_D2z_a+ABZ*I_NAI_K3x2y2z_D2z_a;
  Double I_NAI_K3xy3z_F3z_a = I_NAI_L3xy4z_D2z_a+ABZ*I_NAI_K3xy3z_D2z_a;
  Double I_NAI_K3x4z_F3z_a = I_NAI_L3x5z_D2z_a+ABZ*I_NAI_K3x4z_D2z_a;
  Double I_NAI_K2x5y_F3z_a = I_NAI_L2x5yz_D2z_a+ABZ*I_NAI_K2x5y_D2z_a;
  Double I_NAI_K2x4yz_F3z_a = I_NAI_L2x4y2z_D2z_a+ABZ*I_NAI_K2x4yz_D2z_a;
  Double I_NAI_K2x3y2z_F3z_a = I_NAI_L2x3y3z_D2z_a+ABZ*I_NAI_K2x3y2z_D2z_a;
  Double I_NAI_K2x2y3z_F3z_a = I_NAI_L2x2y4z_D2z_a+ABZ*I_NAI_K2x2y3z_D2z_a;
  Double I_NAI_K2xy4z_F3z_a = I_NAI_L2xy5z_D2z_a+ABZ*I_NAI_K2xy4z_D2z_a;
  Double I_NAI_K2x5z_F3z_a = I_NAI_L2x6z_D2z_a+ABZ*I_NAI_K2x5z_D2z_a;
  Double I_NAI_Kx6y_F3z_a = I_NAI_Lx6yz_D2z_a+ABZ*I_NAI_Kx6y_D2z_a;
  Double I_NAI_Kx5yz_F3z_a = I_NAI_Lx5y2z_D2z_a+ABZ*I_NAI_Kx5yz_D2z_a;
  Double I_NAI_Kx4y2z_F3z_a = I_NAI_Lx4y3z_D2z_a+ABZ*I_NAI_Kx4y2z_D2z_a;
  Double I_NAI_Kx3y3z_F3z_a = I_NAI_Lx3y4z_D2z_a+ABZ*I_NAI_Kx3y3z_D2z_a;
  Double I_NAI_Kx2y4z_F3z_a = I_NAI_Lx2y5z_D2z_a+ABZ*I_NAI_Kx2y4z_D2z_a;
  Double I_NAI_Kxy5z_F3z_a = I_NAI_Lxy6z_D2z_a+ABZ*I_NAI_Kxy5z_D2z_a;
  Double I_NAI_Kx6z_F3z_a = I_NAI_Lx7z_D2z_a+ABZ*I_NAI_Kx6z_D2z_a;
  Double I_NAI_K7y_F3z_a = I_NAI_L7yz_D2z_a+ABZ*I_NAI_K7y_D2z_a;
  Double I_NAI_K6yz_F3z_a = I_NAI_L6y2z_D2z_a+ABZ*I_NAI_K6yz_D2z_a;
  Double I_NAI_K5y2z_F3z_a = I_NAI_L5y3z_D2z_a+ABZ*I_NAI_K5y2z_D2z_a;
  Double I_NAI_K4y3z_F3z_a = I_NAI_L4y4z_D2z_a+ABZ*I_NAI_K4y3z_D2z_a;
  Double I_NAI_K3y4z_F3z_a = I_NAI_L3y5z_D2z_a+ABZ*I_NAI_K3y4z_D2z_a;
  Double I_NAI_K2y5z_F3z_a = I_NAI_L2y6z_D2z_a+ABZ*I_NAI_K2y5z_D2z_a;
  Double I_NAI_Ky6z_F3z_a = I_NAI_Ly7z_D2z_a+ABZ*I_NAI_Ky6z_D2z_a;
  Double I_NAI_K7z_F3z_a = I_NAI_L8z_D2z_a+ABZ*I_NAI_K7z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_I_G_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_F_a
   * RHS shell quartet name: SQ_NAI_I_F_a
   ************************************************************/
  Double I_NAI_I6x_G4x_a = I_NAI_K7x_F3x_a+ABX*I_NAI_I6x_F3x_a;
  Double I_NAI_I5xy_G4x_a = I_NAI_K6xy_F3x_a+ABX*I_NAI_I5xy_F3x_a;
  Double I_NAI_I5xz_G4x_a = I_NAI_K6xz_F3x_a+ABX*I_NAI_I5xz_F3x_a;
  Double I_NAI_I4x2y_G4x_a = I_NAI_K5x2y_F3x_a+ABX*I_NAI_I4x2y_F3x_a;
  Double I_NAI_I4xyz_G4x_a = I_NAI_K5xyz_F3x_a+ABX*I_NAI_I4xyz_F3x_a;
  Double I_NAI_I4x2z_G4x_a = I_NAI_K5x2z_F3x_a+ABX*I_NAI_I4x2z_F3x_a;
  Double I_NAI_I3x3y_G4x_a = I_NAI_K4x3y_F3x_a+ABX*I_NAI_I3x3y_F3x_a;
  Double I_NAI_I3x2yz_G4x_a = I_NAI_K4x2yz_F3x_a+ABX*I_NAI_I3x2yz_F3x_a;
  Double I_NAI_I3xy2z_G4x_a = I_NAI_K4xy2z_F3x_a+ABX*I_NAI_I3xy2z_F3x_a;
  Double I_NAI_I3x3z_G4x_a = I_NAI_K4x3z_F3x_a+ABX*I_NAI_I3x3z_F3x_a;
  Double I_NAI_I2x4y_G4x_a = I_NAI_K3x4y_F3x_a+ABX*I_NAI_I2x4y_F3x_a;
  Double I_NAI_I2x3yz_G4x_a = I_NAI_K3x3yz_F3x_a+ABX*I_NAI_I2x3yz_F3x_a;
  Double I_NAI_I2x2y2z_G4x_a = I_NAI_K3x2y2z_F3x_a+ABX*I_NAI_I2x2y2z_F3x_a;
  Double I_NAI_I2xy3z_G4x_a = I_NAI_K3xy3z_F3x_a+ABX*I_NAI_I2xy3z_F3x_a;
  Double I_NAI_I2x4z_G4x_a = I_NAI_K3x4z_F3x_a+ABX*I_NAI_I2x4z_F3x_a;
  Double I_NAI_Ix5y_G4x_a = I_NAI_K2x5y_F3x_a+ABX*I_NAI_Ix5y_F3x_a;
  Double I_NAI_Ix4yz_G4x_a = I_NAI_K2x4yz_F3x_a+ABX*I_NAI_Ix4yz_F3x_a;
  Double I_NAI_Ix3y2z_G4x_a = I_NAI_K2x3y2z_F3x_a+ABX*I_NAI_Ix3y2z_F3x_a;
  Double I_NAI_Ix2y3z_G4x_a = I_NAI_K2x2y3z_F3x_a+ABX*I_NAI_Ix2y3z_F3x_a;
  Double I_NAI_Ixy4z_G4x_a = I_NAI_K2xy4z_F3x_a+ABX*I_NAI_Ixy4z_F3x_a;
  Double I_NAI_Ix5z_G4x_a = I_NAI_K2x5z_F3x_a+ABX*I_NAI_Ix5z_F3x_a;
  Double I_NAI_I6y_G4x_a = I_NAI_Kx6y_F3x_a+ABX*I_NAI_I6y_F3x_a;
  Double I_NAI_I5yz_G4x_a = I_NAI_Kx5yz_F3x_a+ABX*I_NAI_I5yz_F3x_a;
  Double I_NAI_I4y2z_G4x_a = I_NAI_Kx4y2z_F3x_a+ABX*I_NAI_I4y2z_F3x_a;
  Double I_NAI_I3y3z_G4x_a = I_NAI_Kx3y3z_F3x_a+ABX*I_NAI_I3y3z_F3x_a;
  Double I_NAI_I2y4z_G4x_a = I_NAI_Kx2y4z_F3x_a+ABX*I_NAI_I2y4z_F3x_a;
  Double I_NAI_Iy5z_G4x_a = I_NAI_Kxy5z_F3x_a+ABX*I_NAI_Iy5z_F3x_a;
  Double I_NAI_I6z_G4x_a = I_NAI_Kx6z_F3x_a+ABX*I_NAI_I6z_F3x_a;
  Double I_NAI_I6x_G3xy_a = I_NAI_K6xy_F3x_a+ABY*I_NAI_I6x_F3x_a;
  Double I_NAI_I5xy_G3xy_a = I_NAI_K5x2y_F3x_a+ABY*I_NAI_I5xy_F3x_a;
  Double I_NAI_I5xz_G3xy_a = I_NAI_K5xyz_F3x_a+ABY*I_NAI_I5xz_F3x_a;
  Double I_NAI_I4x2y_G3xy_a = I_NAI_K4x3y_F3x_a+ABY*I_NAI_I4x2y_F3x_a;
  Double I_NAI_I4xyz_G3xy_a = I_NAI_K4x2yz_F3x_a+ABY*I_NAI_I4xyz_F3x_a;
  Double I_NAI_I4x2z_G3xy_a = I_NAI_K4xy2z_F3x_a+ABY*I_NAI_I4x2z_F3x_a;
  Double I_NAI_I3x3y_G3xy_a = I_NAI_K3x4y_F3x_a+ABY*I_NAI_I3x3y_F3x_a;
  Double I_NAI_I3x2yz_G3xy_a = I_NAI_K3x3yz_F3x_a+ABY*I_NAI_I3x2yz_F3x_a;
  Double I_NAI_I3xy2z_G3xy_a = I_NAI_K3x2y2z_F3x_a+ABY*I_NAI_I3xy2z_F3x_a;
  Double I_NAI_I3x3z_G3xy_a = I_NAI_K3xy3z_F3x_a+ABY*I_NAI_I3x3z_F3x_a;
  Double I_NAI_I2x4y_G3xy_a = I_NAI_K2x5y_F3x_a+ABY*I_NAI_I2x4y_F3x_a;
  Double I_NAI_I2x3yz_G3xy_a = I_NAI_K2x4yz_F3x_a+ABY*I_NAI_I2x3yz_F3x_a;
  Double I_NAI_I2x2y2z_G3xy_a = I_NAI_K2x3y2z_F3x_a+ABY*I_NAI_I2x2y2z_F3x_a;
  Double I_NAI_I2xy3z_G3xy_a = I_NAI_K2x2y3z_F3x_a+ABY*I_NAI_I2xy3z_F3x_a;
  Double I_NAI_I2x4z_G3xy_a = I_NAI_K2xy4z_F3x_a+ABY*I_NAI_I2x4z_F3x_a;
  Double I_NAI_Ix5y_G3xy_a = I_NAI_Kx6y_F3x_a+ABY*I_NAI_Ix5y_F3x_a;
  Double I_NAI_Ix4yz_G3xy_a = I_NAI_Kx5yz_F3x_a+ABY*I_NAI_Ix4yz_F3x_a;
  Double I_NAI_Ix3y2z_G3xy_a = I_NAI_Kx4y2z_F3x_a+ABY*I_NAI_Ix3y2z_F3x_a;
  Double I_NAI_Ix2y3z_G3xy_a = I_NAI_Kx3y3z_F3x_a+ABY*I_NAI_Ix2y3z_F3x_a;
  Double I_NAI_Ixy4z_G3xy_a = I_NAI_Kx2y4z_F3x_a+ABY*I_NAI_Ixy4z_F3x_a;
  Double I_NAI_Ix5z_G3xy_a = I_NAI_Kxy5z_F3x_a+ABY*I_NAI_Ix5z_F3x_a;
  Double I_NAI_I6y_G3xy_a = I_NAI_K7y_F3x_a+ABY*I_NAI_I6y_F3x_a;
  Double I_NAI_I5yz_G3xy_a = I_NAI_K6yz_F3x_a+ABY*I_NAI_I5yz_F3x_a;
  Double I_NAI_I4y2z_G3xy_a = I_NAI_K5y2z_F3x_a+ABY*I_NAI_I4y2z_F3x_a;
  Double I_NAI_I3y3z_G3xy_a = I_NAI_K4y3z_F3x_a+ABY*I_NAI_I3y3z_F3x_a;
  Double I_NAI_I2y4z_G3xy_a = I_NAI_K3y4z_F3x_a+ABY*I_NAI_I2y4z_F3x_a;
  Double I_NAI_Iy5z_G3xy_a = I_NAI_K2y5z_F3x_a+ABY*I_NAI_Iy5z_F3x_a;
  Double I_NAI_I6z_G3xy_a = I_NAI_Ky6z_F3x_a+ABY*I_NAI_I6z_F3x_a;
  Double I_NAI_I6x_G3xz_a = I_NAI_K6xz_F3x_a+ABZ*I_NAI_I6x_F3x_a;
  Double I_NAI_I5xy_G3xz_a = I_NAI_K5xyz_F3x_a+ABZ*I_NAI_I5xy_F3x_a;
  Double I_NAI_I5xz_G3xz_a = I_NAI_K5x2z_F3x_a+ABZ*I_NAI_I5xz_F3x_a;
  Double I_NAI_I4x2y_G3xz_a = I_NAI_K4x2yz_F3x_a+ABZ*I_NAI_I4x2y_F3x_a;
  Double I_NAI_I4xyz_G3xz_a = I_NAI_K4xy2z_F3x_a+ABZ*I_NAI_I4xyz_F3x_a;
  Double I_NAI_I4x2z_G3xz_a = I_NAI_K4x3z_F3x_a+ABZ*I_NAI_I4x2z_F3x_a;
  Double I_NAI_I3x3y_G3xz_a = I_NAI_K3x3yz_F3x_a+ABZ*I_NAI_I3x3y_F3x_a;
  Double I_NAI_I3x2yz_G3xz_a = I_NAI_K3x2y2z_F3x_a+ABZ*I_NAI_I3x2yz_F3x_a;
  Double I_NAI_I3xy2z_G3xz_a = I_NAI_K3xy3z_F3x_a+ABZ*I_NAI_I3xy2z_F3x_a;
  Double I_NAI_I3x3z_G3xz_a = I_NAI_K3x4z_F3x_a+ABZ*I_NAI_I3x3z_F3x_a;
  Double I_NAI_I2x4y_G3xz_a = I_NAI_K2x4yz_F3x_a+ABZ*I_NAI_I2x4y_F3x_a;
  Double I_NAI_I2x3yz_G3xz_a = I_NAI_K2x3y2z_F3x_a+ABZ*I_NAI_I2x3yz_F3x_a;
  Double I_NAI_I2x2y2z_G3xz_a = I_NAI_K2x2y3z_F3x_a+ABZ*I_NAI_I2x2y2z_F3x_a;
  Double I_NAI_I2xy3z_G3xz_a = I_NAI_K2xy4z_F3x_a+ABZ*I_NAI_I2xy3z_F3x_a;
  Double I_NAI_I2x4z_G3xz_a = I_NAI_K2x5z_F3x_a+ABZ*I_NAI_I2x4z_F3x_a;
  Double I_NAI_Ix5y_G3xz_a = I_NAI_Kx5yz_F3x_a+ABZ*I_NAI_Ix5y_F3x_a;
  Double I_NAI_Ix4yz_G3xz_a = I_NAI_Kx4y2z_F3x_a+ABZ*I_NAI_Ix4yz_F3x_a;
  Double I_NAI_Ix3y2z_G3xz_a = I_NAI_Kx3y3z_F3x_a+ABZ*I_NAI_Ix3y2z_F3x_a;
  Double I_NAI_Ix2y3z_G3xz_a = I_NAI_Kx2y4z_F3x_a+ABZ*I_NAI_Ix2y3z_F3x_a;
  Double I_NAI_Ixy4z_G3xz_a = I_NAI_Kxy5z_F3x_a+ABZ*I_NAI_Ixy4z_F3x_a;
  Double I_NAI_Ix5z_G3xz_a = I_NAI_Kx6z_F3x_a+ABZ*I_NAI_Ix5z_F3x_a;
  Double I_NAI_I6y_G3xz_a = I_NAI_K6yz_F3x_a+ABZ*I_NAI_I6y_F3x_a;
  Double I_NAI_I5yz_G3xz_a = I_NAI_K5y2z_F3x_a+ABZ*I_NAI_I5yz_F3x_a;
  Double I_NAI_I4y2z_G3xz_a = I_NAI_K4y3z_F3x_a+ABZ*I_NAI_I4y2z_F3x_a;
  Double I_NAI_I3y3z_G3xz_a = I_NAI_K3y4z_F3x_a+ABZ*I_NAI_I3y3z_F3x_a;
  Double I_NAI_I2y4z_G3xz_a = I_NAI_K2y5z_F3x_a+ABZ*I_NAI_I2y4z_F3x_a;
  Double I_NAI_Iy5z_G3xz_a = I_NAI_Ky6z_F3x_a+ABZ*I_NAI_Iy5z_F3x_a;
  Double I_NAI_I6z_G3xz_a = I_NAI_K7z_F3x_a+ABZ*I_NAI_I6z_F3x_a;
  Double I_NAI_I6x_G2x2y_a = I_NAI_K6xy_F2xy_a+ABY*I_NAI_I6x_F2xy_a;
  Double I_NAI_I5xy_G2x2y_a = I_NAI_K5x2y_F2xy_a+ABY*I_NAI_I5xy_F2xy_a;
  Double I_NAI_I5xz_G2x2y_a = I_NAI_K5xyz_F2xy_a+ABY*I_NAI_I5xz_F2xy_a;
  Double I_NAI_I4x2y_G2x2y_a = I_NAI_K4x3y_F2xy_a+ABY*I_NAI_I4x2y_F2xy_a;
  Double I_NAI_I4xyz_G2x2y_a = I_NAI_K4x2yz_F2xy_a+ABY*I_NAI_I4xyz_F2xy_a;
  Double I_NAI_I4x2z_G2x2y_a = I_NAI_K4xy2z_F2xy_a+ABY*I_NAI_I4x2z_F2xy_a;
  Double I_NAI_I3x3y_G2x2y_a = I_NAI_K3x4y_F2xy_a+ABY*I_NAI_I3x3y_F2xy_a;
  Double I_NAI_I3x2yz_G2x2y_a = I_NAI_K3x3yz_F2xy_a+ABY*I_NAI_I3x2yz_F2xy_a;
  Double I_NAI_I3xy2z_G2x2y_a = I_NAI_K3x2y2z_F2xy_a+ABY*I_NAI_I3xy2z_F2xy_a;
  Double I_NAI_I3x3z_G2x2y_a = I_NAI_K3xy3z_F2xy_a+ABY*I_NAI_I3x3z_F2xy_a;
  Double I_NAI_I2x4y_G2x2y_a = I_NAI_K2x5y_F2xy_a+ABY*I_NAI_I2x4y_F2xy_a;
  Double I_NAI_I2x3yz_G2x2y_a = I_NAI_K2x4yz_F2xy_a+ABY*I_NAI_I2x3yz_F2xy_a;
  Double I_NAI_I2x2y2z_G2x2y_a = I_NAI_K2x3y2z_F2xy_a+ABY*I_NAI_I2x2y2z_F2xy_a;
  Double I_NAI_I2xy3z_G2x2y_a = I_NAI_K2x2y3z_F2xy_a+ABY*I_NAI_I2xy3z_F2xy_a;
  Double I_NAI_I2x4z_G2x2y_a = I_NAI_K2xy4z_F2xy_a+ABY*I_NAI_I2x4z_F2xy_a;
  Double I_NAI_Ix5y_G2x2y_a = I_NAI_Kx6y_F2xy_a+ABY*I_NAI_Ix5y_F2xy_a;
  Double I_NAI_Ix4yz_G2x2y_a = I_NAI_Kx5yz_F2xy_a+ABY*I_NAI_Ix4yz_F2xy_a;
  Double I_NAI_Ix3y2z_G2x2y_a = I_NAI_Kx4y2z_F2xy_a+ABY*I_NAI_Ix3y2z_F2xy_a;
  Double I_NAI_Ix2y3z_G2x2y_a = I_NAI_Kx3y3z_F2xy_a+ABY*I_NAI_Ix2y3z_F2xy_a;
  Double I_NAI_Ixy4z_G2x2y_a = I_NAI_Kx2y4z_F2xy_a+ABY*I_NAI_Ixy4z_F2xy_a;
  Double I_NAI_Ix5z_G2x2y_a = I_NAI_Kxy5z_F2xy_a+ABY*I_NAI_Ix5z_F2xy_a;
  Double I_NAI_I6y_G2x2y_a = I_NAI_K7y_F2xy_a+ABY*I_NAI_I6y_F2xy_a;
  Double I_NAI_I5yz_G2x2y_a = I_NAI_K6yz_F2xy_a+ABY*I_NAI_I5yz_F2xy_a;
  Double I_NAI_I4y2z_G2x2y_a = I_NAI_K5y2z_F2xy_a+ABY*I_NAI_I4y2z_F2xy_a;
  Double I_NAI_I3y3z_G2x2y_a = I_NAI_K4y3z_F2xy_a+ABY*I_NAI_I3y3z_F2xy_a;
  Double I_NAI_I2y4z_G2x2y_a = I_NAI_K3y4z_F2xy_a+ABY*I_NAI_I2y4z_F2xy_a;
  Double I_NAI_Iy5z_G2x2y_a = I_NAI_K2y5z_F2xy_a+ABY*I_NAI_Iy5z_F2xy_a;
  Double I_NAI_I6z_G2x2y_a = I_NAI_Ky6z_F2xy_a+ABY*I_NAI_I6z_F2xy_a;
  Double I_NAI_I6x_G2xyz_a = I_NAI_K6xz_F2xy_a+ABZ*I_NAI_I6x_F2xy_a;
  Double I_NAI_I5xy_G2xyz_a = I_NAI_K5xyz_F2xy_a+ABZ*I_NAI_I5xy_F2xy_a;
  Double I_NAI_I5xz_G2xyz_a = I_NAI_K5x2z_F2xy_a+ABZ*I_NAI_I5xz_F2xy_a;
  Double I_NAI_I4x2y_G2xyz_a = I_NAI_K4x2yz_F2xy_a+ABZ*I_NAI_I4x2y_F2xy_a;
  Double I_NAI_I4xyz_G2xyz_a = I_NAI_K4xy2z_F2xy_a+ABZ*I_NAI_I4xyz_F2xy_a;
  Double I_NAI_I4x2z_G2xyz_a = I_NAI_K4x3z_F2xy_a+ABZ*I_NAI_I4x2z_F2xy_a;
  Double I_NAI_I3x3y_G2xyz_a = I_NAI_K3x3yz_F2xy_a+ABZ*I_NAI_I3x3y_F2xy_a;
  Double I_NAI_I3x2yz_G2xyz_a = I_NAI_K3x2y2z_F2xy_a+ABZ*I_NAI_I3x2yz_F2xy_a;
  Double I_NAI_I3xy2z_G2xyz_a = I_NAI_K3xy3z_F2xy_a+ABZ*I_NAI_I3xy2z_F2xy_a;
  Double I_NAI_I3x3z_G2xyz_a = I_NAI_K3x4z_F2xy_a+ABZ*I_NAI_I3x3z_F2xy_a;
  Double I_NAI_I2x4y_G2xyz_a = I_NAI_K2x4yz_F2xy_a+ABZ*I_NAI_I2x4y_F2xy_a;
  Double I_NAI_I2x3yz_G2xyz_a = I_NAI_K2x3y2z_F2xy_a+ABZ*I_NAI_I2x3yz_F2xy_a;
  Double I_NAI_I2x2y2z_G2xyz_a = I_NAI_K2x2y3z_F2xy_a+ABZ*I_NAI_I2x2y2z_F2xy_a;
  Double I_NAI_I2xy3z_G2xyz_a = I_NAI_K2xy4z_F2xy_a+ABZ*I_NAI_I2xy3z_F2xy_a;
  Double I_NAI_I2x4z_G2xyz_a = I_NAI_K2x5z_F2xy_a+ABZ*I_NAI_I2x4z_F2xy_a;
  Double I_NAI_Ix5y_G2xyz_a = I_NAI_Kx5yz_F2xy_a+ABZ*I_NAI_Ix5y_F2xy_a;
  Double I_NAI_Ix4yz_G2xyz_a = I_NAI_Kx4y2z_F2xy_a+ABZ*I_NAI_Ix4yz_F2xy_a;
  Double I_NAI_Ix3y2z_G2xyz_a = I_NAI_Kx3y3z_F2xy_a+ABZ*I_NAI_Ix3y2z_F2xy_a;
  Double I_NAI_Ix2y3z_G2xyz_a = I_NAI_Kx2y4z_F2xy_a+ABZ*I_NAI_Ix2y3z_F2xy_a;
  Double I_NAI_Ixy4z_G2xyz_a = I_NAI_Kxy5z_F2xy_a+ABZ*I_NAI_Ixy4z_F2xy_a;
  Double I_NAI_Ix5z_G2xyz_a = I_NAI_Kx6z_F2xy_a+ABZ*I_NAI_Ix5z_F2xy_a;
  Double I_NAI_I6y_G2xyz_a = I_NAI_K6yz_F2xy_a+ABZ*I_NAI_I6y_F2xy_a;
  Double I_NAI_I5yz_G2xyz_a = I_NAI_K5y2z_F2xy_a+ABZ*I_NAI_I5yz_F2xy_a;
  Double I_NAI_I4y2z_G2xyz_a = I_NAI_K4y3z_F2xy_a+ABZ*I_NAI_I4y2z_F2xy_a;
  Double I_NAI_I3y3z_G2xyz_a = I_NAI_K3y4z_F2xy_a+ABZ*I_NAI_I3y3z_F2xy_a;
  Double I_NAI_I2y4z_G2xyz_a = I_NAI_K2y5z_F2xy_a+ABZ*I_NAI_I2y4z_F2xy_a;
  Double I_NAI_Iy5z_G2xyz_a = I_NAI_Ky6z_F2xy_a+ABZ*I_NAI_Iy5z_F2xy_a;
  Double I_NAI_I6z_G2xyz_a = I_NAI_K7z_F2xy_a+ABZ*I_NAI_I6z_F2xy_a;
  Double I_NAI_I6x_G2x2z_a = I_NAI_K6xz_F2xz_a+ABZ*I_NAI_I6x_F2xz_a;
  Double I_NAI_I5xy_G2x2z_a = I_NAI_K5xyz_F2xz_a+ABZ*I_NAI_I5xy_F2xz_a;
  Double I_NAI_I5xz_G2x2z_a = I_NAI_K5x2z_F2xz_a+ABZ*I_NAI_I5xz_F2xz_a;
  Double I_NAI_I4x2y_G2x2z_a = I_NAI_K4x2yz_F2xz_a+ABZ*I_NAI_I4x2y_F2xz_a;
  Double I_NAI_I4xyz_G2x2z_a = I_NAI_K4xy2z_F2xz_a+ABZ*I_NAI_I4xyz_F2xz_a;
  Double I_NAI_I4x2z_G2x2z_a = I_NAI_K4x3z_F2xz_a+ABZ*I_NAI_I4x2z_F2xz_a;
  Double I_NAI_I3x3y_G2x2z_a = I_NAI_K3x3yz_F2xz_a+ABZ*I_NAI_I3x3y_F2xz_a;
  Double I_NAI_I3x2yz_G2x2z_a = I_NAI_K3x2y2z_F2xz_a+ABZ*I_NAI_I3x2yz_F2xz_a;
  Double I_NAI_I3xy2z_G2x2z_a = I_NAI_K3xy3z_F2xz_a+ABZ*I_NAI_I3xy2z_F2xz_a;
  Double I_NAI_I3x3z_G2x2z_a = I_NAI_K3x4z_F2xz_a+ABZ*I_NAI_I3x3z_F2xz_a;
  Double I_NAI_I2x4y_G2x2z_a = I_NAI_K2x4yz_F2xz_a+ABZ*I_NAI_I2x4y_F2xz_a;
  Double I_NAI_I2x3yz_G2x2z_a = I_NAI_K2x3y2z_F2xz_a+ABZ*I_NAI_I2x3yz_F2xz_a;
  Double I_NAI_I2x2y2z_G2x2z_a = I_NAI_K2x2y3z_F2xz_a+ABZ*I_NAI_I2x2y2z_F2xz_a;
  Double I_NAI_I2xy3z_G2x2z_a = I_NAI_K2xy4z_F2xz_a+ABZ*I_NAI_I2xy3z_F2xz_a;
  Double I_NAI_I2x4z_G2x2z_a = I_NAI_K2x5z_F2xz_a+ABZ*I_NAI_I2x4z_F2xz_a;
  Double I_NAI_Ix5y_G2x2z_a = I_NAI_Kx5yz_F2xz_a+ABZ*I_NAI_Ix5y_F2xz_a;
  Double I_NAI_Ix4yz_G2x2z_a = I_NAI_Kx4y2z_F2xz_a+ABZ*I_NAI_Ix4yz_F2xz_a;
  Double I_NAI_Ix3y2z_G2x2z_a = I_NAI_Kx3y3z_F2xz_a+ABZ*I_NAI_Ix3y2z_F2xz_a;
  Double I_NAI_Ix2y3z_G2x2z_a = I_NAI_Kx2y4z_F2xz_a+ABZ*I_NAI_Ix2y3z_F2xz_a;
  Double I_NAI_Ixy4z_G2x2z_a = I_NAI_Kxy5z_F2xz_a+ABZ*I_NAI_Ixy4z_F2xz_a;
  Double I_NAI_Ix5z_G2x2z_a = I_NAI_Kx6z_F2xz_a+ABZ*I_NAI_Ix5z_F2xz_a;
  Double I_NAI_I6y_G2x2z_a = I_NAI_K6yz_F2xz_a+ABZ*I_NAI_I6y_F2xz_a;
  Double I_NAI_I5yz_G2x2z_a = I_NAI_K5y2z_F2xz_a+ABZ*I_NAI_I5yz_F2xz_a;
  Double I_NAI_I4y2z_G2x2z_a = I_NAI_K4y3z_F2xz_a+ABZ*I_NAI_I4y2z_F2xz_a;
  Double I_NAI_I3y3z_G2x2z_a = I_NAI_K3y4z_F2xz_a+ABZ*I_NAI_I3y3z_F2xz_a;
  Double I_NAI_I2y4z_G2x2z_a = I_NAI_K2y5z_F2xz_a+ABZ*I_NAI_I2y4z_F2xz_a;
  Double I_NAI_Iy5z_G2x2z_a = I_NAI_Ky6z_F2xz_a+ABZ*I_NAI_Iy5z_F2xz_a;
  Double I_NAI_I6z_G2x2z_a = I_NAI_K7z_F2xz_a+ABZ*I_NAI_I6z_F2xz_a;
  Double I_NAI_I6x_Gx3y_a = I_NAI_K7x_F3y_a+ABX*I_NAI_I6x_F3y_a;
  Double I_NAI_I5xy_Gx3y_a = I_NAI_K6xy_F3y_a+ABX*I_NAI_I5xy_F3y_a;
  Double I_NAI_I5xz_Gx3y_a = I_NAI_K6xz_F3y_a+ABX*I_NAI_I5xz_F3y_a;
  Double I_NAI_I4x2y_Gx3y_a = I_NAI_K5x2y_F3y_a+ABX*I_NAI_I4x2y_F3y_a;
  Double I_NAI_I4xyz_Gx3y_a = I_NAI_K5xyz_F3y_a+ABX*I_NAI_I4xyz_F3y_a;
  Double I_NAI_I4x2z_Gx3y_a = I_NAI_K5x2z_F3y_a+ABX*I_NAI_I4x2z_F3y_a;
  Double I_NAI_I3x3y_Gx3y_a = I_NAI_K4x3y_F3y_a+ABX*I_NAI_I3x3y_F3y_a;
  Double I_NAI_I3x2yz_Gx3y_a = I_NAI_K4x2yz_F3y_a+ABX*I_NAI_I3x2yz_F3y_a;
  Double I_NAI_I3xy2z_Gx3y_a = I_NAI_K4xy2z_F3y_a+ABX*I_NAI_I3xy2z_F3y_a;
  Double I_NAI_I3x3z_Gx3y_a = I_NAI_K4x3z_F3y_a+ABX*I_NAI_I3x3z_F3y_a;
  Double I_NAI_I2x4y_Gx3y_a = I_NAI_K3x4y_F3y_a+ABX*I_NAI_I2x4y_F3y_a;
  Double I_NAI_I2x3yz_Gx3y_a = I_NAI_K3x3yz_F3y_a+ABX*I_NAI_I2x3yz_F3y_a;
  Double I_NAI_I2x2y2z_Gx3y_a = I_NAI_K3x2y2z_F3y_a+ABX*I_NAI_I2x2y2z_F3y_a;
  Double I_NAI_I2xy3z_Gx3y_a = I_NAI_K3xy3z_F3y_a+ABX*I_NAI_I2xy3z_F3y_a;
  Double I_NAI_I2x4z_Gx3y_a = I_NAI_K3x4z_F3y_a+ABX*I_NAI_I2x4z_F3y_a;
  Double I_NAI_Ix5y_Gx3y_a = I_NAI_K2x5y_F3y_a+ABX*I_NAI_Ix5y_F3y_a;
  Double I_NAI_Ix4yz_Gx3y_a = I_NAI_K2x4yz_F3y_a+ABX*I_NAI_Ix4yz_F3y_a;
  Double I_NAI_Ix3y2z_Gx3y_a = I_NAI_K2x3y2z_F3y_a+ABX*I_NAI_Ix3y2z_F3y_a;
  Double I_NAI_Ix2y3z_Gx3y_a = I_NAI_K2x2y3z_F3y_a+ABX*I_NAI_Ix2y3z_F3y_a;
  Double I_NAI_Ixy4z_Gx3y_a = I_NAI_K2xy4z_F3y_a+ABX*I_NAI_Ixy4z_F3y_a;
  Double I_NAI_Ix5z_Gx3y_a = I_NAI_K2x5z_F3y_a+ABX*I_NAI_Ix5z_F3y_a;
  Double I_NAI_I6y_Gx3y_a = I_NAI_Kx6y_F3y_a+ABX*I_NAI_I6y_F3y_a;
  Double I_NAI_I5yz_Gx3y_a = I_NAI_Kx5yz_F3y_a+ABX*I_NAI_I5yz_F3y_a;
  Double I_NAI_I4y2z_Gx3y_a = I_NAI_Kx4y2z_F3y_a+ABX*I_NAI_I4y2z_F3y_a;
  Double I_NAI_I3y3z_Gx3y_a = I_NAI_Kx3y3z_F3y_a+ABX*I_NAI_I3y3z_F3y_a;
  Double I_NAI_I2y4z_Gx3y_a = I_NAI_Kx2y4z_F3y_a+ABX*I_NAI_I2y4z_F3y_a;
  Double I_NAI_Iy5z_Gx3y_a = I_NAI_Kxy5z_F3y_a+ABX*I_NAI_Iy5z_F3y_a;
  Double I_NAI_I6z_Gx3y_a = I_NAI_Kx6z_F3y_a+ABX*I_NAI_I6z_F3y_a;
  Double I_NAI_I6x_Gx2yz_a = I_NAI_K6xz_Fx2y_a+ABZ*I_NAI_I6x_Fx2y_a;
  Double I_NAI_I5xy_Gx2yz_a = I_NAI_K5xyz_Fx2y_a+ABZ*I_NAI_I5xy_Fx2y_a;
  Double I_NAI_I5xz_Gx2yz_a = I_NAI_K5x2z_Fx2y_a+ABZ*I_NAI_I5xz_Fx2y_a;
  Double I_NAI_I4x2y_Gx2yz_a = I_NAI_K4x2yz_Fx2y_a+ABZ*I_NAI_I4x2y_Fx2y_a;
  Double I_NAI_I4xyz_Gx2yz_a = I_NAI_K4xy2z_Fx2y_a+ABZ*I_NAI_I4xyz_Fx2y_a;
  Double I_NAI_I4x2z_Gx2yz_a = I_NAI_K4x3z_Fx2y_a+ABZ*I_NAI_I4x2z_Fx2y_a;
  Double I_NAI_I3x3y_Gx2yz_a = I_NAI_K3x3yz_Fx2y_a+ABZ*I_NAI_I3x3y_Fx2y_a;
  Double I_NAI_I3x2yz_Gx2yz_a = I_NAI_K3x2y2z_Fx2y_a+ABZ*I_NAI_I3x2yz_Fx2y_a;
  Double I_NAI_I3xy2z_Gx2yz_a = I_NAI_K3xy3z_Fx2y_a+ABZ*I_NAI_I3xy2z_Fx2y_a;
  Double I_NAI_I3x3z_Gx2yz_a = I_NAI_K3x4z_Fx2y_a+ABZ*I_NAI_I3x3z_Fx2y_a;
  Double I_NAI_I2x4y_Gx2yz_a = I_NAI_K2x4yz_Fx2y_a+ABZ*I_NAI_I2x4y_Fx2y_a;
  Double I_NAI_I2x3yz_Gx2yz_a = I_NAI_K2x3y2z_Fx2y_a+ABZ*I_NAI_I2x3yz_Fx2y_a;
  Double I_NAI_I2x2y2z_Gx2yz_a = I_NAI_K2x2y3z_Fx2y_a+ABZ*I_NAI_I2x2y2z_Fx2y_a;
  Double I_NAI_I2xy3z_Gx2yz_a = I_NAI_K2xy4z_Fx2y_a+ABZ*I_NAI_I2xy3z_Fx2y_a;
  Double I_NAI_I2x4z_Gx2yz_a = I_NAI_K2x5z_Fx2y_a+ABZ*I_NAI_I2x4z_Fx2y_a;
  Double I_NAI_Ix5y_Gx2yz_a = I_NAI_Kx5yz_Fx2y_a+ABZ*I_NAI_Ix5y_Fx2y_a;
  Double I_NAI_Ix4yz_Gx2yz_a = I_NAI_Kx4y2z_Fx2y_a+ABZ*I_NAI_Ix4yz_Fx2y_a;
  Double I_NAI_Ix3y2z_Gx2yz_a = I_NAI_Kx3y3z_Fx2y_a+ABZ*I_NAI_Ix3y2z_Fx2y_a;
  Double I_NAI_Ix2y3z_Gx2yz_a = I_NAI_Kx2y4z_Fx2y_a+ABZ*I_NAI_Ix2y3z_Fx2y_a;
  Double I_NAI_Ixy4z_Gx2yz_a = I_NAI_Kxy5z_Fx2y_a+ABZ*I_NAI_Ixy4z_Fx2y_a;
  Double I_NAI_Ix5z_Gx2yz_a = I_NAI_Kx6z_Fx2y_a+ABZ*I_NAI_Ix5z_Fx2y_a;
  Double I_NAI_I6y_Gx2yz_a = I_NAI_K6yz_Fx2y_a+ABZ*I_NAI_I6y_Fx2y_a;
  Double I_NAI_I5yz_Gx2yz_a = I_NAI_K5y2z_Fx2y_a+ABZ*I_NAI_I5yz_Fx2y_a;
  Double I_NAI_I4y2z_Gx2yz_a = I_NAI_K4y3z_Fx2y_a+ABZ*I_NAI_I4y2z_Fx2y_a;
  Double I_NAI_I3y3z_Gx2yz_a = I_NAI_K3y4z_Fx2y_a+ABZ*I_NAI_I3y3z_Fx2y_a;
  Double I_NAI_I2y4z_Gx2yz_a = I_NAI_K2y5z_Fx2y_a+ABZ*I_NAI_I2y4z_Fx2y_a;
  Double I_NAI_Iy5z_Gx2yz_a = I_NAI_Ky6z_Fx2y_a+ABZ*I_NAI_Iy5z_Fx2y_a;
  Double I_NAI_I6z_Gx2yz_a = I_NAI_K7z_Fx2y_a+ABZ*I_NAI_I6z_Fx2y_a;
  Double I_NAI_I6x_Gxy2z_a = I_NAI_K6xy_Fx2z_a+ABY*I_NAI_I6x_Fx2z_a;
  Double I_NAI_I5xy_Gxy2z_a = I_NAI_K5x2y_Fx2z_a+ABY*I_NAI_I5xy_Fx2z_a;
  Double I_NAI_I5xz_Gxy2z_a = I_NAI_K5xyz_Fx2z_a+ABY*I_NAI_I5xz_Fx2z_a;
  Double I_NAI_I4x2y_Gxy2z_a = I_NAI_K4x3y_Fx2z_a+ABY*I_NAI_I4x2y_Fx2z_a;
  Double I_NAI_I4xyz_Gxy2z_a = I_NAI_K4x2yz_Fx2z_a+ABY*I_NAI_I4xyz_Fx2z_a;
  Double I_NAI_I4x2z_Gxy2z_a = I_NAI_K4xy2z_Fx2z_a+ABY*I_NAI_I4x2z_Fx2z_a;
  Double I_NAI_I3x3y_Gxy2z_a = I_NAI_K3x4y_Fx2z_a+ABY*I_NAI_I3x3y_Fx2z_a;
  Double I_NAI_I3x2yz_Gxy2z_a = I_NAI_K3x3yz_Fx2z_a+ABY*I_NAI_I3x2yz_Fx2z_a;
  Double I_NAI_I3xy2z_Gxy2z_a = I_NAI_K3x2y2z_Fx2z_a+ABY*I_NAI_I3xy2z_Fx2z_a;
  Double I_NAI_I3x3z_Gxy2z_a = I_NAI_K3xy3z_Fx2z_a+ABY*I_NAI_I3x3z_Fx2z_a;
  Double I_NAI_I2x4y_Gxy2z_a = I_NAI_K2x5y_Fx2z_a+ABY*I_NAI_I2x4y_Fx2z_a;
  Double I_NAI_I2x3yz_Gxy2z_a = I_NAI_K2x4yz_Fx2z_a+ABY*I_NAI_I2x3yz_Fx2z_a;
  Double I_NAI_I2x2y2z_Gxy2z_a = I_NAI_K2x3y2z_Fx2z_a+ABY*I_NAI_I2x2y2z_Fx2z_a;
  Double I_NAI_I2xy3z_Gxy2z_a = I_NAI_K2x2y3z_Fx2z_a+ABY*I_NAI_I2xy3z_Fx2z_a;
  Double I_NAI_I2x4z_Gxy2z_a = I_NAI_K2xy4z_Fx2z_a+ABY*I_NAI_I2x4z_Fx2z_a;
  Double I_NAI_Ix5y_Gxy2z_a = I_NAI_Kx6y_Fx2z_a+ABY*I_NAI_Ix5y_Fx2z_a;
  Double I_NAI_Ix4yz_Gxy2z_a = I_NAI_Kx5yz_Fx2z_a+ABY*I_NAI_Ix4yz_Fx2z_a;
  Double I_NAI_Ix3y2z_Gxy2z_a = I_NAI_Kx4y2z_Fx2z_a+ABY*I_NAI_Ix3y2z_Fx2z_a;
  Double I_NAI_Ix2y3z_Gxy2z_a = I_NAI_Kx3y3z_Fx2z_a+ABY*I_NAI_Ix2y3z_Fx2z_a;
  Double I_NAI_Ixy4z_Gxy2z_a = I_NAI_Kx2y4z_Fx2z_a+ABY*I_NAI_Ixy4z_Fx2z_a;
  Double I_NAI_Ix5z_Gxy2z_a = I_NAI_Kxy5z_Fx2z_a+ABY*I_NAI_Ix5z_Fx2z_a;
  Double I_NAI_I6y_Gxy2z_a = I_NAI_K7y_Fx2z_a+ABY*I_NAI_I6y_Fx2z_a;
  Double I_NAI_I5yz_Gxy2z_a = I_NAI_K6yz_Fx2z_a+ABY*I_NAI_I5yz_Fx2z_a;
  Double I_NAI_I4y2z_Gxy2z_a = I_NAI_K5y2z_Fx2z_a+ABY*I_NAI_I4y2z_Fx2z_a;
  Double I_NAI_I3y3z_Gxy2z_a = I_NAI_K4y3z_Fx2z_a+ABY*I_NAI_I3y3z_Fx2z_a;
  Double I_NAI_I2y4z_Gxy2z_a = I_NAI_K3y4z_Fx2z_a+ABY*I_NAI_I2y4z_Fx2z_a;
  Double I_NAI_Iy5z_Gxy2z_a = I_NAI_K2y5z_Fx2z_a+ABY*I_NAI_Iy5z_Fx2z_a;
  Double I_NAI_I6z_Gxy2z_a = I_NAI_Ky6z_Fx2z_a+ABY*I_NAI_I6z_Fx2z_a;
  Double I_NAI_I6x_Gx3z_a = I_NAI_K7x_F3z_a+ABX*I_NAI_I6x_F3z_a;
  Double I_NAI_I5xy_Gx3z_a = I_NAI_K6xy_F3z_a+ABX*I_NAI_I5xy_F3z_a;
  Double I_NAI_I5xz_Gx3z_a = I_NAI_K6xz_F3z_a+ABX*I_NAI_I5xz_F3z_a;
  Double I_NAI_I4x2y_Gx3z_a = I_NAI_K5x2y_F3z_a+ABX*I_NAI_I4x2y_F3z_a;
  Double I_NAI_I4xyz_Gx3z_a = I_NAI_K5xyz_F3z_a+ABX*I_NAI_I4xyz_F3z_a;
  Double I_NAI_I4x2z_Gx3z_a = I_NAI_K5x2z_F3z_a+ABX*I_NAI_I4x2z_F3z_a;
  Double I_NAI_I3x3y_Gx3z_a = I_NAI_K4x3y_F3z_a+ABX*I_NAI_I3x3y_F3z_a;
  Double I_NAI_I3x2yz_Gx3z_a = I_NAI_K4x2yz_F3z_a+ABX*I_NAI_I3x2yz_F3z_a;
  Double I_NAI_I3xy2z_Gx3z_a = I_NAI_K4xy2z_F3z_a+ABX*I_NAI_I3xy2z_F3z_a;
  Double I_NAI_I3x3z_Gx3z_a = I_NAI_K4x3z_F3z_a+ABX*I_NAI_I3x3z_F3z_a;
  Double I_NAI_I2x4y_Gx3z_a = I_NAI_K3x4y_F3z_a+ABX*I_NAI_I2x4y_F3z_a;
  Double I_NAI_I2x3yz_Gx3z_a = I_NAI_K3x3yz_F3z_a+ABX*I_NAI_I2x3yz_F3z_a;
  Double I_NAI_I2x2y2z_Gx3z_a = I_NAI_K3x2y2z_F3z_a+ABX*I_NAI_I2x2y2z_F3z_a;
  Double I_NAI_I2xy3z_Gx3z_a = I_NAI_K3xy3z_F3z_a+ABX*I_NAI_I2xy3z_F3z_a;
  Double I_NAI_I2x4z_Gx3z_a = I_NAI_K3x4z_F3z_a+ABX*I_NAI_I2x4z_F3z_a;
  Double I_NAI_Ix5y_Gx3z_a = I_NAI_K2x5y_F3z_a+ABX*I_NAI_Ix5y_F3z_a;
  Double I_NAI_Ix4yz_Gx3z_a = I_NAI_K2x4yz_F3z_a+ABX*I_NAI_Ix4yz_F3z_a;
  Double I_NAI_Ix3y2z_Gx3z_a = I_NAI_K2x3y2z_F3z_a+ABX*I_NAI_Ix3y2z_F3z_a;
  Double I_NAI_Ix2y3z_Gx3z_a = I_NAI_K2x2y3z_F3z_a+ABX*I_NAI_Ix2y3z_F3z_a;
  Double I_NAI_Ixy4z_Gx3z_a = I_NAI_K2xy4z_F3z_a+ABX*I_NAI_Ixy4z_F3z_a;
  Double I_NAI_Ix5z_Gx3z_a = I_NAI_K2x5z_F3z_a+ABX*I_NAI_Ix5z_F3z_a;
  Double I_NAI_I6y_Gx3z_a = I_NAI_Kx6y_F3z_a+ABX*I_NAI_I6y_F3z_a;
  Double I_NAI_I5yz_Gx3z_a = I_NAI_Kx5yz_F3z_a+ABX*I_NAI_I5yz_F3z_a;
  Double I_NAI_I4y2z_Gx3z_a = I_NAI_Kx4y2z_F3z_a+ABX*I_NAI_I4y2z_F3z_a;
  Double I_NAI_I3y3z_Gx3z_a = I_NAI_Kx3y3z_F3z_a+ABX*I_NAI_I3y3z_F3z_a;
  Double I_NAI_I2y4z_Gx3z_a = I_NAI_Kx2y4z_F3z_a+ABX*I_NAI_I2y4z_F3z_a;
  Double I_NAI_Iy5z_Gx3z_a = I_NAI_Kxy5z_F3z_a+ABX*I_NAI_Iy5z_F3z_a;
  Double I_NAI_I6z_Gx3z_a = I_NAI_Kx6z_F3z_a+ABX*I_NAI_I6z_F3z_a;
  Double I_NAI_I6x_G4y_a = I_NAI_K6xy_F3y_a+ABY*I_NAI_I6x_F3y_a;
  Double I_NAI_I5xy_G4y_a = I_NAI_K5x2y_F3y_a+ABY*I_NAI_I5xy_F3y_a;
  Double I_NAI_I5xz_G4y_a = I_NAI_K5xyz_F3y_a+ABY*I_NAI_I5xz_F3y_a;
  Double I_NAI_I4x2y_G4y_a = I_NAI_K4x3y_F3y_a+ABY*I_NAI_I4x2y_F3y_a;
  Double I_NAI_I4xyz_G4y_a = I_NAI_K4x2yz_F3y_a+ABY*I_NAI_I4xyz_F3y_a;
  Double I_NAI_I4x2z_G4y_a = I_NAI_K4xy2z_F3y_a+ABY*I_NAI_I4x2z_F3y_a;
  Double I_NAI_I3x3y_G4y_a = I_NAI_K3x4y_F3y_a+ABY*I_NAI_I3x3y_F3y_a;
  Double I_NAI_I3x2yz_G4y_a = I_NAI_K3x3yz_F3y_a+ABY*I_NAI_I3x2yz_F3y_a;
  Double I_NAI_I3xy2z_G4y_a = I_NAI_K3x2y2z_F3y_a+ABY*I_NAI_I3xy2z_F3y_a;
  Double I_NAI_I3x3z_G4y_a = I_NAI_K3xy3z_F3y_a+ABY*I_NAI_I3x3z_F3y_a;
  Double I_NAI_I2x4y_G4y_a = I_NAI_K2x5y_F3y_a+ABY*I_NAI_I2x4y_F3y_a;
  Double I_NAI_I2x3yz_G4y_a = I_NAI_K2x4yz_F3y_a+ABY*I_NAI_I2x3yz_F3y_a;
  Double I_NAI_I2x2y2z_G4y_a = I_NAI_K2x3y2z_F3y_a+ABY*I_NAI_I2x2y2z_F3y_a;
  Double I_NAI_I2xy3z_G4y_a = I_NAI_K2x2y3z_F3y_a+ABY*I_NAI_I2xy3z_F3y_a;
  Double I_NAI_I2x4z_G4y_a = I_NAI_K2xy4z_F3y_a+ABY*I_NAI_I2x4z_F3y_a;
  Double I_NAI_Ix5y_G4y_a = I_NAI_Kx6y_F3y_a+ABY*I_NAI_Ix5y_F3y_a;
  Double I_NAI_Ix4yz_G4y_a = I_NAI_Kx5yz_F3y_a+ABY*I_NAI_Ix4yz_F3y_a;
  Double I_NAI_Ix3y2z_G4y_a = I_NAI_Kx4y2z_F3y_a+ABY*I_NAI_Ix3y2z_F3y_a;
  Double I_NAI_Ix2y3z_G4y_a = I_NAI_Kx3y3z_F3y_a+ABY*I_NAI_Ix2y3z_F3y_a;
  Double I_NAI_Ixy4z_G4y_a = I_NAI_Kx2y4z_F3y_a+ABY*I_NAI_Ixy4z_F3y_a;
  Double I_NAI_Ix5z_G4y_a = I_NAI_Kxy5z_F3y_a+ABY*I_NAI_Ix5z_F3y_a;
  Double I_NAI_I6y_G4y_a = I_NAI_K7y_F3y_a+ABY*I_NAI_I6y_F3y_a;
  Double I_NAI_I5yz_G4y_a = I_NAI_K6yz_F3y_a+ABY*I_NAI_I5yz_F3y_a;
  Double I_NAI_I4y2z_G4y_a = I_NAI_K5y2z_F3y_a+ABY*I_NAI_I4y2z_F3y_a;
  Double I_NAI_I3y3z_G4y_a = I_NAI_K4y3z_F3y_a+ABY*I_NAI_I3y3z_F3y_a;
  Double I_NAI_I2y4z_G4y_a = I_NAI_K3y4z_F3y_a+ABY*I_NAI_I2y4z_F3y_a;
  Double I_NAI_Iy5z_G4y_a = I_NAI_K2y5z_F3y_a+ABY*I_NAI_Iy5z_F3y_a;
  Double I_NAI_I6z_G4y_a = I_NAI_Ky6z_F3y_a+ABY*I_NAI_I6z_F3y_a;
  Double I_NAI_I6x_G3yz_a = I_NAI_K6xz_F3y_a+ABZ*I_NAI_I6x_F3y_a;
  Double I_NAI_I5xy_G3yz_a = I_NAI_K5xyz_F3y_a+ABZ*I_NAI_I5xy_F3y_a;
  Double I_NAI_I5xz_G3yz_a = I_NAI_K5x2z_F3y_a+ABZ*I_NAI_I5xz_F3y_a;
  Double I_NAI_I4x2y_G3yz_a = I_NAI_K4x2yz_F3y_a+ABZ*I_NAI_I4x2y_F3y_a;
  Double I_NAI_I4xyz_G3yz_a = I_NAI_K4xy2z_F3y_a+ABZ*I_NAI_I4xyz_F3y_a;
  Double I_NAI_I4x2z_G3yz_a = I_NAI_K4x3z_F3y_a+ABZ*I_NAI_I4x2z_F3y_a;
  Double I_NAI_I3x3y_G3yz_a = I_NAI_K3x3yz_F3y_a+ABZ*I_NAI_I3x3y_F3y_a;
  Double I_NAI_I3x2yz_G3yz_a = I_NAI_K3x2y2z_F3y_a+ABZ*I_NAI_I3x2yz_F3y_a;
  Double I_NAI_I3xy2z_G3yz_a = I_NAI_K3xy3z_F3y_a+ABZ*I_NAI_I3xy2z_F3y_a;
  Double I_NAI_I3x3z_G3yz_a = I_NAI_K3x4z_F3y_a+ABZ*I_NAI_I3x3z_F3y_a;
  Double I_NAI_I2x4y_G3yz_a = I_NAI_K2x4yz_F3y_a+ABZ*I_NAI_I2x4y_F3y_a;
  Double I_NAI_I2x3yz_G3yz_a = I_NAI_K2x3y2z_F3y_a+ABZ*I_NAI_I2x3yz_F3y_a;
  Double I_NAI_I2x2y2z_G3yz_a = I_NAI_K2x2y3z_F3y_a+ABZ*I_NAI_I2x2y2z_F3y_a;
  Double I_NAI_I2xy3z_G3yz_a = I_NAI_K2xy4z_F3y_a+ABZ*I_NAI_I2xy3z_F3y_a;
  Double I_NAI_I2x4z_G3yz_a = I_NAI_K2x5z_F3y_a+ABZ*I_NAI_I2x4z_F3y_a;
  Double I_NAI_Ix5y_G3yz_a = I_NAI_Kx5yz_F3y_a+ABZ*I_NAI_Ix5y_F3y_a;
  Double I_NAI_Ix4yz_G3yz_a = I_NAI_Kx4y2z_F3y_a+ABZ*I_NAI_Ix4yz_F3y_a;
  Double I_NAI_Ix3y2z_G3yz_a = I_NAI_Kx3y3z_F3y_a+ABZ*I_NAI_Ix3y2z_F3y_a;
  Double I_NAI_Ix2y3z_G3yz_a = I_NAI_Kx2y4z_F3y_a+ABZ*I_NAI_Ix2y3z_F3y_a;
  Double I_NAI_Ixy4z_G3yz_a = I_NAI_Kxy5z_F3y_a+ABZ*I_NAI_Ixy4z_F3y_a;
  Double I_NAI_Ix5z_G3yz_a = I_NAI_Kx6z_F3y_a+ABZ*I_NAI_Ix5z_F3y_a;
  Double I_NAI_I6y_G3yz_a = I_NAI_K6yz_F3y_a+ABZ*I_NAI_I6y_F3y_a;
  Double I_NAI_I5yz_G3yz_a = I_NAI_K5y2z_F3y_a+ABZ*I_NAI_I5yz_F3y_a;
  Double I_NAI_I4y2z_G3yz_a = I_NAI_K4y3z_F3y_a+ABZ*I_NAI_I4y2z_F3y_a;
  Double I_NAI_I3y3z_G3yz_a = I_NAI_K3y4z_F3y_a+ABZ*I_NAI_I3y3z_F3y_a;
  Double I_NAI_I2y4z_G3yz_a = I_NAI_K2y5z_F3y_a+ABZ*I_NAI_I2y4z_F3y_a;
  Double I_NAI_Iy5z_G3yz_a = I_NAI_Ky6z_F3y_a+ABZ*I_NAI_Iy5z_F3y_a;
  Double I_NAI_I6z_G3yz_a = I_NAI_K7z_F3y_a+ABZ*I_NAI_I6z_F3y_a;
  Double I_NAI_I6x_G2y2z_a = I_NAI_K6xz_F2yz_a+ABZ*I_NAI_I6x_F2yz_a;
  Double I_NAI_I5xy_G2y2z_a = I_NAI_K5xyz_F2yz_a+ABZ*I_NAI_I5xy_F2yz_a;
  Double I_NAI_I5xz_G2y2z_a = I_NAI_K5x2z_F2yz_a+ABZ*I_NAI_I5xz_F2yz_a;
  Double I_NAI_I4x2y_G2y2z_a = I_NAI_K4x2yz_F2yz_a+ABZ*I_NAI_I4x2y_F2yz_a;
  Double I_NAI_I4xyz_G2y2z_a = I_NAI_K4xy2z_F2yz_a+ABZ*I_NAI_I4xyz_F2yz_a;
  Double I_NAI_I4x2z_G2y2z_a = I_NAI_K4x3z_F2yz_a+ABZ*I_NAI_I4x2z_F2yz_a;
  Double I_NAI_I3x3y_G2y2z_a = I_NAI_K3x3yz_F2yz_a+ABZ*I_NAI_I3x3y_F2yz_a;
  Double I_NAI_I3x2yz_G2y2z_a = I_NAI_K3x2y2z_F2yz_a+ABZ*I_NAI_I3x2yz_F2yz_a;
  Double I_NAI_I3xy2z_G2y2z_a = I_NAI_K3xy3z_F2yz_a+ABZ*I_NAI_I3xy2z_F2yz_a;
  Double I_NAI_I3x3z_G2y2z_a = I_NAI_K3x4z_F2yz_a+ABZ*I_NAI_I3x3z_F2yz_a;
  Double I_NAI_I2x4y_G2y2z_a = I_NAI_K2x4yz_F2yz_a+ABZ*I_NAI_I2x4y_F2yz_a;
  Double I_NAI_I2x3yz_G2y2z_a = I_NAI_K2x3y2z_F2yz_a+ABZ*I_NAI_I2x3yz_F2yz_a;
  Double I_NAI_I2x2y2z_G2y2z_a = I_NAI_K2x2y3z_F2yz_a+ABZ*I_NAI_I2x2y2z_F2yz_a;
  Double I_NAI_I2xy3z_G2y2z_a = I_NAI_K2xy4z_F2yz_a+ABZ*I_NAI_I2xy3z_F2yz_a;
  Double I_NAI_I2x4z_G2y2z_a = I_NAI_K2x5z_F2yz_a+ABZ*I_NAI_I2x4z_F2yz_a;
  Double I_NAI_Ix5y_G2y2z_a = I_NAI_Kx5yz_F2yz_a+ABZ*I_NAI_Ix5y_F2yz_a;
  Double I_NAI_Ix4yz_G2y2z_a = I_NAI_Kx4y2z_F2yz_a+ABZ*I_NAI_Ix4yz_F2yz_a;
  Double I_NAI_Ix3y2z_G2y2z_a = I_NAI_Kx3y3z_F2yz_a+ABZ*I_NAI_Ix3y2z_F2yz_a;
  Double I_NAI_Ix2y3z_G2y2z_a = I_NAI_Kx2y4z_F2yz_a+ABZ*I_NAI_Ix2y3z_F2yz_a;
  Double I_NAI_Ixy4z_G2y2z_a = I_NAI_Kxy5z_F2yz_a+ABZ*I_NAI_Ixy4z_F2yz_a;
  Double I_NAI_Ix5z_G2y2z_a = I_NAI_Kx6z_F2yz_a+ABZ*I_NAI_Ix5z_F2yz_a;
  Double I_NAI_I6y_G2y2z_a = I_NAI_K6yz_F2yz_a+ABZ*I_NAI_I6y_F2yz_a;
  Double I_NAI_I5yz_G2y2z_a = I_NAI_K5y2z_F2yz_a+ABZ*I_NAI_I5yz_F2yz_a;
  Double I_NAI_I4y2z_G2y2z_a = I_NAI_K4y3z_F2yz_a+ABZ*I_NAI_I4y2z_F2yz_a;
  Double I_NAI_I3y3z_G2y2z_a = I_NAI_K3y4z_F2yz_a+ABZ*I_NAI_I3y3z_F2yz_a;
  Double I_NAI_I2y4z_G2y2z_a = I_NAI_K2y5z_F2yz_a+ABZ*I_NAI_I2y4z_F2yz_a;
  Double I_NAI_Iy5z_G2y2z_a = I_NAI_Ky6z_F2yz_a+ABZ*I_NAI_Iy5z_F2yz_a;
  Double I_NAI_I6z_G2y2z_a = I_NAI_K7z_F2yz_a+ABZ*I_NAI_I6z_F2yz_a;
  Double I_NAI_I6x_Gy3z_a = I_NAI_K6xy_F3z_a+ABY*I_NAI_I6x_F3z_a;
  Double I_NAI_I5xy_Gy3z_a = I_NAI_K5x2y_F3z_a+ABY*I_NAI_I5xy_F3z_a;
  Double I_NAI_I5xz_Gy3z_a = I_NAI_K5xyz_F3z_a+ABY*I_NAI_I5xz_F3z_a;
  Double I_NAI_I4x2y_Gy3z_a = I_NAI_K4x3y_F3z_a+ABY*I_NAI_I4x2y_F3z_a;
  Double I_NAI_I4xyz_Gy3z_a = I_NAI_K4x2yz_F3z_a+ABY*I_NAI_I4xyz_F3z_a;
  Double I_NAI_I4x2z_Gy3z_a = I_NAI_K4xy2z_F3z_a+ABY*I_NAI_I4x2z_F3z_a;
  Double I_NAI_I3x3y_Gy3z_a = I_NAI_K3x4y_F3z_a+ABY*I_NAI_I3x3y_F3z_a;
  Double I_NAI_I3x2yz_Gy3z_a = I_NAI_K3x3yz_F3z_a+ABY*I_NAI_I3x2yz_F3z_a;
  Double I_NAI_I3xy2z_Gy3z_a = I_NAI_K3x2y2z_F3z_a+ABY*I_NAI_I3xy2z_F3z_a;
  Double I_NAI_I3x3z_Gy3z_a = I_NAI_K3xy3z_F3z_a+ABY*I_NAI_I3x3z_F3z_a;
  Double I_NAI_I2x4y_Gy3z_a = I_NAI_K2x5y_F3z_a+ABY*I_NAI_I2x4y_F3z_a;
  Double I_NAI_I2x3yz_Gy3z_a = I_NAI_K2x4yz_F3z_a+ABY*I_NAI_I2x3yz_F3z_a;
  Double I_NAI_I2x2y2z_Gy3z_a = I_NAI_K2x3y2z_F3z_a+ABY*I_NAI_I2x2y2z_F3z_a;
  Double I_NAI_I2xy3z_Gy3z_a = I_NAI_K2x2y3z_F3z_a+ABY*I_NAI_I2xy3z_F3z_a;
  Double I_NAI_I2x4z_Gy3z_a = I_NAI_K2xy4z_F3z_a+ABY*I_NAI_I2x4z_F3z_a;
  Double I_NAI_Ix5y_Gy3z_a = I_NAI_Kx6y_F3z_a+ABY*I_NAI_Ix5y_F3z_a;
  Double I_NAI_Ix4yz_Gy3z_a = I_NAI_Kx5yz_F3z_a+ABY*I_NAI_Ix4yz_F3z_a;
  Double I_NAI_Ix3y2z_Gy3z_a = I_NAI_Kx4y2z_F3z_a+ABY*I_NAI_Ix3y2z_F3z_a;
  Double I_NAI_Ix2y3z_Gy3z_a = I_NAI_Kx3y3z_F3z_a+ABY*I_NAI_Ix2y3z_F3z_a;
  Double I_NAI_Ixy4z_Gy3z_a = I_NAI_Kx2y4z_F3z_a+ABY*I_NAI_Ixy4z_F3z_a;
  Double I_NAI_Ix5z_Gy3z_a = I_NAI_Kxy5z_F3z_a+ABY*I_NAI_Ix5z_F3z_a;
  Double I_NAI_I6y_Gy3z_a = I_NAI_K7y_F3z_a+ABY*I_NAI_I6y_F3z_a;
  Double I_NAI_I5yz_Gy3z_a = I_NAI_K6yz_F3z_a+ABY*I_NAI_I5yz_F3z_a;
  Double I_NAI_I4y2z_Gy3z_a = I_NAI_K5y2z_F3z_a+ABY*I_NAI_I4y2z_F3z_a;
  Double I_NAI_I3y3z_Gy3z_a = I_NAI_K4y3z_F3z_a+ABY*I_NAI_I3y3z_F3z_a;
  Double I_NAI_I2y4z_Gy3z_a = I_NAI_K3y4z_F3z_a+ABY*I_NAI_I2y4z_F3z_a;
  Double I_NAI_Iy5z_Gy3z_a = I_NAI_K2y5z_F3z_a+ABY*I_NAI_Iy5z_F3z_a;
  Double I_NAI_I6z_Gy3z_a = I_NAI_Ky6z_F3z_a+ABY*I_NAI_I6z_F3z_a;
  Double I_NAI_I6x_G4z_a = I_NAI_K6xz_F3z_a+ABZ*I_NAI_I6x_F3z_a;
  Double I_NAI_I5xy_G4z_a = I_NAI_K5xyz_F3z_a+ABZ*I_NAI_I5xy_F3z_a;
  Double I_NAI_I5xz_G4z_a = I_NAI_K5x2z_F3z_a+ABZ*I_NAI_I5xz_F3z_a;
  Double I_NAI_I4x2y_G4z_a = I_NAI_K4x2yz_F3z_a+ABZ*I_NAI_I4x2y_F3z_a;
  Double I_NAI_I4xyz_G4z_a = I_NAI_K4xy2z_F3z_a+ABZ*I_NAI_I4xyz_F3z_a;
  Double I_NAI_I4x2z_G4z_a = I_NAI_K4x3z_F3z_a+ABZ*I_NAI_I4x2z_F3z_a;
  Double I_NAI_I3x3y_G4z_a = I_NAI_K3x3yz_F3z_a+ABZ*I_NAI_I3x3y_F3z_a;
  Double I_NAI_I3x2yz_G4z_a = I_NAI_K3x2y2z_F3z_a+ABZ*I_NAI_I3x2yz_F3z_a;
  Double I_NAI_I3xy2z_G4z_a = I_NAI_K3xy3z_F3z_a+ABZ*I_NAI_I3xy2z_F3z_a;
  Double I_NAI_I3x3z_G4z_a = I_NAI_K3x4z_F3z_a+ABZ*I_NAI_I3x3z_F3z_a;
  Double I_NAI_I2x4y_G4z_a = I_NAI_K2x4yz_F3z_a+ABZ*I_NAI_I2x4y_F3z_a;
  Double I_NAI_I2x3yz_G4z_a = I_NAI_K2x3y2z_F3z_a+ABZ*I_NAI_I2x3yz_F3z_a;
  Double I_NAI_I2x2y2z_G4z_a = I_NAI_K2x2y3z_F3z_a+ABZ*I_NAI_I2x2y2z_F3z_a;
  Double I_NAI_I2xy3z_G4z_a = I_NAI_K2xy4z_F3z_a+ABZ*I_NAI_I2xy3z_F3z_a;
  Double I_NAI_I2x4z_G4z_a = I_NAI_K2x5z_F3z_a+ABZ*I_NAI_I2x4z_F3z_a;
  Double I_NAI_Ix5y_G4z_a = I_NAI_Kx5yz_F3z_a+ABZ*I_NAI_Ix5y_F3z_a;
  Double I_NAI_Ix4yz_G4z_a = I_NAI_Kx4y2z_F3z_a+ABZ*I_NAI_Ix4yz_F3z_a;
  Double I_NAI_Ix3y2z_G4z_a = I_NAI_Kx3y3z_F3z_a+ABZ*I_NAI_Ix3y2z_F3z_a;
  Double I_NAI_Ix2y3z_G4z_a = I_NAI_Kx2y4z_F3z_a+ABZ*I_NAI_Ix2y3z_F3z_a;
  Double I_NAI_Ixy4z_G4z_a = I_NAI_Kxy5z_F3z_a+ABZ*I_NAI_Ixy4z_F3z_a;
  Double I_NAI_Ix5z_G4z_a = I_NAI_Kx6z_F3z_a+ABZ*I_NAI_Ix5z_F3z_a;
  Double I_NAI_I6y_G4z_a = I_NAI_K6yz_F3z_a+ABZ*I_NAI_I6y_F3z_a;
  Double I_NAI_I5yz_G4z_a = I_NAI_K5y2z_F3z_a+ABZ*I_NAI_I5yz_F3z_a;
  Double I_NAI_I4y2z_G4z_a = I_NAI_K4y3z_F3z_a+ABZ*I_NAI_I4y2z_F3z_a;
  Double I_NAI_I3y3z_G4z_a = I_NAI_K3y4z_F3z_a+ABZ*I_NAI_I3y3z_F3z_a;
  Double I_NAI_I2y4z_G4z_a = I_NAI_K2y5z_F3z_a+ABZ*I_NAI_I2y4z_F3z_a;
  Double I_NAI_Iy5z_G4z_a = I_NAI_Ky6z_F3z_a+ABZ*I_NAI_Iy5z_F3z_a;
  Double I_NAI_I6z_G4z_a = I_NAI_K7z_F3z_a+ABZ*I_NAI_I6z_F3z_a;

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
   * totally 0 integrals are omitted 
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
  Double I_NAI_K7x_Py_b = I_NAI_L7xy_S_b+ABY*I_NAI_K7x_S_b;
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
  Double I_NAI_K7x_Pz_b = I_NAI_L7xz_S_b+ABZ*I_NAI_K7x_S_b;
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
  Double I_NAI_K7y_Pz_b = I_NAI_L7yz_S_b+ABZ*I_NAI_K7y_S_b;
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
   * totally 84 integrals are omitted 
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
   * totally 14 integrals are omitted 
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
  Double I_NAI_L7yz_Px_b = I_NAI_Mx7yz_S_b+ABX*I_NAI_L7yz_S_b;
  Double I_NAI_L6y2z_Px_b = I_NAI_Mx6y2z_S_b+ABX*I_NAI_L6y2z_S_b;
  Double I_NAI_L5y3z_Px_b = I_NAI_Mx5y3z_S_b+ABX*I_NAI_L5y3z_S_b;
  Double I_NAI_L4y4z_Px_b = I_NAI_Mx4y4z_S_b+ABX*I_NAI_L4y4z_S_b;
  Double I_NAI_L3y5z_Px_b = I_NAI_Mx3y5z_S_b+ABX*I_NAI_L3y5z_S_b;
  Double I_NAI_L2y6z_Px_b = I_NAI_Mx2y6z_S_b+ABX*I_NAI_L2y6z_S_b;
  Double I_NAI_Ly7z_Px_b = I_NAI_Mxy7z_S_b+ABX*I_NAI_Ly7z_S_b;
  Double I_NAI_L7xy_Py_b = I_NAI_M7x2y_S_b+ABY*I_NAI_L7xy_S_b;
  Double I_NAI_L6x2y_Py_b = I_NAI_M6x3y_S_b+ABY*I_NAI_L6x2y_S_b;
  Double I_NAI_L6xyz_Py_b = I_NAI_M6x2yz_S_b+ABY*I_NAI_L6xyz_S_b;
  Double I_NAI_L6x2z_Py_b = I_NAI_M6xy2z_S_b+ABY*I_NAI_L6x2z_S_b;
  Double I_NAI_L5x3y_Py_b = I_NAI_M5x4y_S_b+ABY*I_NAI_L5x3y_S_b;
  Double I_NAI_L5x2yz_Py_b = I_NAI_M5x3yz_S_b+ABY*I_NAI_L5x2yz_S_b;
  Double I_NAI_L5xy2z_Py_b = I_NAI_M5x2y2z_S_b+ABY*I_NAI_L5xy2z_S_b;
  Double I_NAI_L5x3z_Py_b = I_NAI_M5xy3z_S_b+ABY*I_NAI_L5x3z_S_b;
  Double I_NAI_L4x4y_Py_b = I_NAI_M4x5y_S_b+ABY*I_NAI_L4x4y_S_b;
  Double I_NAI_L4x3yz_Py_b = I_NAI_M4x4yz_S_b+ABY*I_NAI_L4x3yz_S_b;
  Double I_NAI_L4x2y2z_Py_b = I_NAI_M4x3y2z_S_b+ABY*I_NAI_L4x2y2z_S_b;
  Double I_NAI_L4xy3z_Py_b = I_NAI_M4x2y3z_S_b+ABY*I_NAI_L4xy3z_S_b;
  Double I_NAI_L4x4z_Py_b = I_NAI_M4xy4z_S_b+ABY*I_NAI_L4x4z_S_b;
  Double I_NAI_L3x5y_Py_b = I_NAI_M3x6y_S_b+ABY*I_NAI_L3x5y_S_b;
  Double I_NAI_L3x4yz_Py_b = I_NAI_M3x5yz_S_b+ABY*I_NAI_L3x4yz_S_b;
  Double I_NAI_L3x3y2z_Py_b = I_NAI_M3x4y2z_S_b+ABY*I_NAI_L3x3y2z_S_b;
  Double I_NAI_L3x2y3z_Py_b = I_NAI_M3x3y3z_S_b+ABY*I_NAI_L3x2y3z_S_b;
  Double I_NAI_L3xy4z_Py_b = I_NAI_M3x2y4z_S_b+ABY*I_NAI_L3xy4z_S_b;
  Double I_NAI_L3x5z_Py_b = I_NAI_M3xy5z_S_b+ABY*I_NAI_L3x5z_S_b;
  Double I_NAI_L2x6y_Py_b = I_NAI_M2x7y_S_b+ABY*I_NAI_L2x6y_S_b;
  Double I_NAI_L2x5yz_Py_b = I_NAI_M2x6yz_S_b+ABY*I_NAI_L2x5yz_S_b;
  Double I_NAI_L2x4y2z_Py_b = I_NAI_M2x5y2z_S_b+ABY*I_NAI_L2x4y2z_S_b;
  Double I_NAI_L2x3y3z_Py_b = I_NAI_M2x4y3z_S_b+ABY*I_NAI_L2x3y3z_S_b;
  Double I_NAI_L2x2y4z_Py_b = I_NAI_M2x3y4z_S_b+ABY*I_NAI_L2x2y4z_S_b;
  Double I_NAI_L2xy5z_Py_b = I_NAI_M2x2y5z_S_b+ABY*I_NAI_L2xy5z_S_b;
  Double I_NAI_L2x6z_Py_b = I_NAI_M2xy6z_S_b+ABY*I_NAI_L2x6z_S_b;
  Double I_NAI_Lx7y_Py_b = I_NAI_Mx8y_S_b+ABY*I_NAI_Lx7y_S_b;
  Double I_NAI_Lx6yz_Py_b = I_NAI_Mx7yz_S_b+ABY*I_NAI_Lx6yz_S_b;
  Double I_NAI_Lx5y2z_Py_b = I_NAI_Mx6y2z_S_b+ABY*I_NAI_Lx5y2z_S_b;
  Double I_NAI_Lx4y3z_Py_b = I_NAI_Mx5y3z_S_b+ABY*I_NAI_Lx4y3z_S_b;
  Double I_NAI_Lx3y4z_Py_b = I_NAI_Mx4y4z_S_b+ABY*I_NAI_Lx3y4z_S_b;
  Double I_NAI_Lx2y5z_Py_b = I_NAI_Mx3y5z_S_b+ABY*I_NAI_Lx2y5z_S_b;
  Double I_NAI_Lxy6z_Py_b = I_NAI_Mx2y6z_S_b+ABY*I_NAI_Lxy6z_S_b;
  Double I_NAI_Lx7z_Py_b = I_NAI_Mxy7z_S_b+ABY*I_NAI_Lx7z_S_b;
  Double I_NAI_L8y_Py_b = I_NAI_M9y_S_b+ABY*I_NAI_L8y_S_b;
  Double I_NAI_L7yz_Py_b = I_NAI_M8yz_S_b+ABY*I_NAI_L7yz_S_b;
  Double I_NAI_L6y2z_Py_b = I_NAI_M7y2z_S_b+ABY*I_NAI_L6y2z_S_b;
  Double I_NAI_L5y3z_Py_b = I_NAI_M6y3z_S_b+ABY*I_NAI_L5y3z_S_b;
  Double I_NAI_L4y4z_Py_b = I_NAI_M5y4z_S_b+ABY*I_NAI_L4y4z_S_b;
  Double I_NAI_L3y5z_Py_b = I_NAI_M4y5z_S_b+ABY*I_NAI_L3y5z_S_b;
  Double I_NAI_L2y6z_Py_b = I_NAI_M3y6z_S_b+ABY*I_NAI_L2y6z_S_b;
  Double I_NAI_Ly7z_Py_b = I_NAI_M2y7z_S_b+ABY*I_NAI_Ly7z_S_b;
  Double I_NAI_L7xz_Pz_b = I_NAI_M7x2z_S_b+ABZ*I_NAI_L7xz_S_b;
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
  Double I_NAI_L7yz_Pz_b = I_NAI_M7y2z_S_b+ABZ*I_NAI_L7yz_S_b;
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
   * totally 108 integrals are omitted 
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
  Double I_NAI_K7x_D2y_b = I_NAI_L7xy_Py_b+ABY*I_NAI_K7x_Py_b;
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
  Double I_NAI_K7x_D2z_b = I_NAI_L7xz_Pz_b+ABZ*I_NAI_K7x_Pz_b;
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
  Double I_NAI_K7y_D2z_b = I_NAI_L7yz_Pz_b+ABZ*I_NAI_K7y_Pz_b;
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
   * totally 115 integrals are omitted 
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
  Double I_NAI_I5xy_F2xz_b = I_NAI_K5xyz_D2x_b+ABZ*I_NAI_I5xy_D2x_b;
  Double I_NAI_I5xz_F2xz_b = I_NAI_K5x2z_D2x_b+ABZ*I_NAI_I5xz_D2x_b;
  Double I_NAI_I4x2y_F2xz_b = I_NAI_K4x2yz_D2x_b+ABZ*I_NAI_I4x2y_D2x_b;
  Double I_NAI_I4xyz_F2xz_b = I_NAI_K4xy2z_D2x_b+ABZ*I_NAI_I4xyz_D2x_b;
  Double I_NAI_I4x2z_F2xz_b = I_NAI_K4x3z_D2x_b+ABZ*I_NAI_I4x2z_D2x_b;
  Double I_NAI_I3x3y_F2xz_b = I_NAI_K3x3yz_D2x_b+ABZ*I_NAI_I3x3y_D2x_b;
  Double I_NAI_I3x2yz_F2xz_b = I_NAI_K3x2y2z_D2x_b+ABZ*I_NAI_I3x2yz_D2x_b;
  Double I_NAI_I3xy2z_F2xz_b = I_NAI_K3xy3z_D2x_b+ABZ*I_NAI_I3xy2z_D2x_b;
  Double I_NAI_I3x3z_F2xz_b = I_NAI_K3x4z_D2x_b+ABZ*I_NAI_I3x3z_D2x_b;
  Double I_NAI_I2x4y_F2xz_b = I_NAI_K2x4yz_D2x_b+ABZ*I_NAI_I2x4y_D2x_b;
  Double I_NAI_I2x3yz_F2xz_b = I_NAI_K2x3y2z_D2x_b+ABZ*I_NAI_I2x3yz_D2x_b;
  Double I_NAI_I2x2y2z_F2xz_b = I_NAI_K2x2y3z_D2x_b+ABZ*I_NAI_I2x2y2z_D2x_b;
  Double I_NAI_I2xy3z_F2xz_b = I_NAI_K2xy4z_D2x_b+ABZ*I_NAI_I2xy3z_D2x_b;
  Double I_NAI_I2x4z_F2xz_b = I_NAI_K2x5z_D2x_b+ABZ*I_NAI_I2x4z_D2x_b;
  Double I_NAI_Ix5y_F2xz_b = I_NAI_Kx5yz_D2x_b+ABZ*I_NAI_Ix5y_D2x_b;
  Double I_NAI_Ix4yz_F2xz_b = I_NAI_Kx4y2z_D2x_b+ABZ*I_NAI_Ix4yz_D2x_b;
  Double I_NAI_Ix3y2z_F2xz_b = I_NAI_Kx3y3z_D2x_b+ABZ*I_NAI_Ix3y2z_D2x_b;
  Double I_NAI_Ix2y3z_F2xz_b = I_NAI_Kx2y4z_D2x_b+ABZ*I_NAI_Ix2y3z_D2x_b;
  Double I_NAI_Ixy4z_F2xz_b = I_NAI_Kxy5z_D2x_b+ABZ*I_NAI_Ixy4z_D2x_b;
  Double I_NAI_Ix5z_F2xz_b = I_NAI_Kx6z_D2x_b+ABZ*I_NAI_Ix5z_D2x_b;
  Double I_NAI_I6y_F2xz_b = I_NAI_K6yz_D2x_b+ABZ*I_NAI_I6y_D2x_b;
  Double I_NAI_I5yz_F2xz_b = I_NAI_K5y2z_D2x_b+ABZ*I_NAI_I5yz_D2x_b;
  Double I_NAI_I4y2z_F2xz_b = I_NAI_K4y3z_D2x_b+ABZ*I_NAI_I4y2z_D2x_b;
  Double I_NAI_I3y3z_F2xz_b = I_NAI_K3y4z_D2x_b+ABZ*I_NAI_I3y3z_D2x_b;
  Double I_NAI_I2y4z_F2xz_b = I_NAI_K2y5z_D2x_b+ABZ*I_NAI_I2y4z_D2x_b;
  Double I_NAI_Iy5z_F2xz_b = I_NAI_Ky6z_D2x_b+ABZ*I_NAI_Iy5z_D2x_b;
  Double I_NAI_I6z_F2xz_b = I_NAI_K7z_D2x_b+ABZ*I_NAI_I6z_D2x_b;
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
  Double I_NAI_I6x_F2yz_b = I_NAI_K6xz_D2y_b+ABZ*I_NAI_I6x_D2y_b;
  Double I_NAI_I5xy_F2yz_b = I_NAI_K5xyz_D2y_b+ABZ*I_NAI_I5xy_D2y_b;
  Double I_NAI_I5xz_F2yz_b = I_NAI_K5x2z_D2y_b+ABZ*I_NAI_I5xz_D2y_b;
  Double I_NAI_I4x2y_F2yz_b = I_NAI_K4x2yz_D2y_b+ABZ*I_NAI_I4x2y_D2y_b;
  Double I_NAI_I4xyz_F2yz_b = I_NAI_K4xy2z_D2y_b+ABZ*I_NAI_I4xyz_D2y_b;
  Double I_NAI_I4x2z_F2yz_b = I_NAI_K4x3z_D2y_b+ABZ*I_NAI_I4x2z_D2y_b;
  Double I_NAI_I3x3y_F2yz_b = I_NAI_K3x3yz_D2y_b+ABZ*I_NAI_I3x3y_D2y_b;
  Double I_NAI_I3x2yz_F2yz_b = I_NAI_K3x2y2z_D2y_b+ABZ*I_NAI_I3x2yz_D2y_b;
  Double I_NAI_I3xy2z_F2yz_b = I_NAI_K3xy3z_D2y_b+ABZ*I_NAI_I3xy2z_D2y_b;
  Double I_NAI_I3x3z_F2yz_b = I_NAI_K3x4z_D2y_b+ABZ*I_NAI_I3x3z_D2y_b;
  Double I_NAI_I2x4y_F2yz_b = I_NAI_K2x4yz_D2y_b+ABZ*I_NAI_I2x4y_D2y_b;
  Double I_NAI_I2x3yz_F2yz_b = I_NAI_K2x3y2z_D2y_b+ABZ*I_NAI_I2x3yz_D2y_b;
  Double I_NAI_I2x2y2z_F2yz_b = I_NAI_K2x2y3z_D2y_b+ABZ*I_NAI_I2x2y2z_D2y_b;
  Double I_NAI_I2xy3z_F2yz_b = I_NAI_K2xy4z_D2y_b+ABZ*I_NAI_I2xy3z_D2y_b;
  Double I_NAI_I2x4z_F2yz_b = I_NAI_K2x5z_D2y_b+ABZ*I_NAI_I2x4z_D2y_b;
  Double I_NAI_Ix5y_F2yz_b = I_NAI_Kx5yz_D2y_b+ABZ*I_NAI_Ix5y_D2y_b;
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
   * totally 63 integrals are omitted 
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
   * shell quartet name: SQ_NAI_M_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 44 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_N_S_b
   * RHS shell quartet name: SQ_NAI_M_S_b
   ************************************************************/
  Double I_NAI_M9x_Px_b = I_NAI_N10x_S_b+ABX*I_NAI_M9x_S_b;
  Double I_NAI_M8xy_Px_b = I_NAI_N9xy_S_b+ABX*I_NAI_M8xy_S_b;
  Double I_NAI_M8xz_Px_b = I_NAI_N9xz_S_b+ABX*I_NAI_M8xz_S_b;
  Double I_NAI_M7x2y_Px_b = I_NAI_N8x2y_S_b+ABX*I_NAI_M7x2y_S_b;
  Double I_NAI_M7xyz_Px_b = I_NAI_N8xyz_S_b+ABX*I_NAI_M7xyz_S_b;
  Double I_NAI_M7x2z_Px_b = I_NAI_N8x2z_S_b+ABX*I_NAI_M7x2z_S_b;
  Double I_NAI_M6x3y_Px_b = I_NAI_N7x3y_S_b+ABX*I_NAI_M6x3y_S_b;
  Double I_NAI_M6x2yz_Px_b = I_NAI_N7x2yz_S_b+ABX*I_NAI_M6x2yz_S_b;
  Double I_NAI_M6xy2z_Px_b = I_NAI_N7xy2z_S_b+ABX*I_NAI_M6xy2z_S_b;
  Double I_NAI_M6x3z_Px_b = I_NAI_N7x3z_S_b+ABX*I_NAI_M6x3z_S_b;
  Double I_NAI_M5x4y_Px_b = I_NAI_N6x4y_S_b+ABX*I_NAI_M5x4y_S_b;
  Double I_NAI_M5x3yz_Px_b = I_NAI_N6x3yz_S_b+ABX*I_NAI_M5x3yz_S_b;
  Double I_NAI_M5x2y2z_Px_b = I_NAI_N6x2y2z_S_b+ABX*I_NAI_M5x2y2z_S_b;
  Double I_NAI_M5xy3z_Px_b = I_NAI_N6xy3z_S_b+ABX*I_NAI_M5xy3z_S_b;
  Double I_NAI_M5x4z_Px_b = I_NAI_N6x4z_S_b+ABX*I_NAI_M5x4z_S_b;
  Double I_NAI_M4x5y_Px_b = I_NAI_N5x5y_S_b+ABX*I_NAI_M4x5y_S_b;
  Double I_NAI_M4x4yz_Px_b = I_NAI_N5x4yz_S_b+ABX*I_NAI_M4x4yz_S_b;
  Double I_NAI_M4x3y2z_Px_b = I_NAI_N5x3y2z_S_b+ABX*I_NAI_M4x3y2z_S_b;
  Double I_NAI_M4x2y3z_Px_b = I_NAI_N5x2y3z_S_b+ABX*I_NAI_M4x2y3z_S_b;
  Double I_NAI_M4xy4z_Px_b = I_NAI_N5xy4z_S_b+ABX*I_NAI_M4xy4z_S_b;
  Double I_NAI_M4x5z_Px_b = I_NAI_N5x5z_S_b+ABX*I_NAI_M4x5z_S_b;
  Double I_NAI_M3x6y_Px_b = I_NAI_N4x6y_S_b+ABX*I_NAI_M3x6y_S_b;
  Double I_NAI_M3x5yz_Px_b = I_NAI_N4x5yz_S_b+ABX*I_NAI_M3x5yz_S_b;
  Double I_NAI_M3x4y2z_Px_b = I_NAI_N4x4y2z_S_b+ABX*I_NAI_M3x4y2z_S_b;
  Double I_NAI_M3x3y3z_Px_b = I_NAI_N4x3y3z_S_b+ABX*I_NAI_M3x3y3z_S_b;
  Double I_NAI_M3x2y4z_Px_b = I_NAI_N4x2y4z_S_b+ABX*I_NAI_M3x2y4z_S_b;
  Double I_NAI_M3xy5z_Px_b = I_NAI_N4xy5z_S_b+ABX*I_NAI_M3xy5z_S_b;
  Double I_NAI_M3x6z_Px_b = I_NAI_N4x6z_S_b+ABX*I_NAI_M3x6z_S_b;
  Double I_NAI_M2x7y_Px_b = I_NAI_N3x7y_S_b+ABX*I_NAI_M2x7y_S_b;
  Double I_NAI_M2x6yz_Px_b = I_NAI_N3x6yz_S_b+ABX*I_NAI_M2x6yz_S_b;
  Double I_NAI_M2x5y2z_Px_b = I_NAI_N3x5y2z_S_b+ABX*I_NAI_M2x5y2z_S_b;
  Double I_NAI_M2x4y3z_Px_b = I_NAI_N3x4y3z_S_b+ABX*I_NAI_M2x4y3z_S_b;
  Double I_NAI_M2x3y4z_Px_b = I_NAI_N3x3y4z_S_b+ABX*I_NAI_M2x3y4z_S_b;
  Double I_NAI_M2x2y5z_Px_b = I_NAI_N3x2y5z_S_b+ABX*I_NAI_M2x2y5z_S_b;
  Double I_NAI_M2xy6z_Px_b = I_NAI_N3xy6z_S_b+ABX*I_NAI_M2xy6z_S_b;
  Double I_NAI_M2x7z_Px_b = I_NAI_N3x7z_S_b+ABX*I_NAI_M2x7z_S_b;
  Double I_NAI_Mx7yz_Px_b = I_NAI_N2x7yz_S_b+ABX*I_NAI_Mx7yz_S_b;
  Double I_NAI_Mx6y2z_Px_b = I_NAI_N2x6y2z_S_b+ABX*I_NAI_Mx6y2z_S_b;
  Double I_NAI_Mx5y3z_Px_b = I_NAI_N2x5y3z_S_b+ABX*I_NAI_Mx5y3z_S_b;
  Double I_NAI_Mx4y4z_Px_b = I_NAI_N2x4y4z_S_b+ABX*I_NAI_Mx4y4z_S_b;
  Double I_NAI_Mx3y5z_Px_b = I_NAI_N2x3y5z_S_b+ABX*I_NAI_Mx3y5z_S_b;
  Double I_NAI_Mx2y6z_Px_b = I_NAI_N2x2y6z_S_b+ABX*I_NAI_Mx2y6z_S_b;
  Double I_NAI_Mxy7z_Px_b = I_NAI_N2xy7z_S_b+ABX*I_NAI_Mxy7z_S_b;
  Double I_NAI_M7x2y_Py_b = I_NAI_N7x3y_S_b+ABY*I_NAI_M7x2y_S_b;
  Double I_NAI_M6x3y_Py_b = I_NAI_N6x4y_S_b+ABY*I_NAI_M6x3y_S_b;
  Double I_NAI_M6x2yz_Py_b = I_NAI_N6x3yz_S_b+ABY*I_NAI_M6x2yz_S_b;
  Double I_NAI_M6xy2z_Py_b = I_NAI_N6x2y2z_S_b+ABY*I_NAI_M6xy2z_S_b;
  Double I_NAI_M5x4y_Py_b = I_NAI_N5x5y_S_b+ABY*I_NAI_M5x4y_S_b;
  Double I_NAI_M5x3yz_Py_b = I_NAI_N5x4yz_S_b+ABY*I_NAI_M5x3yz_S_b;
  Double I_NAI_M5x2y2z_Py_b = I_NAI_N5x3y2z_S_b+ABY*I_NAI_M5x2y2z_S_b;
  Double I_NAI_M5xy3z_Py_b = I_NAI_N5x2y3z_S_b+ABY*I_NAI_M5xy3z_S_b;
  Double I_NAI_M4x5y_Py_b = I_NAI_N4x6y_S_b+ABY*I_NAI_M4x5y_S_b;
  Double I_NAI_M4x4yz_Py_b = I_NAI_N4x5yz_S_b+ABY*I_NAI_M4x4yz_S_b;
  Double I_NAI_M4x3y2z_Py_b = I_NAI_N4x4y2z_S_b+ABY*I_NAI_M4x3y2z_S_b;
  Double I_NAI_M4x2y3z_Py_b = I_NAI_N4x3y3z_S_b+ABY*I_NAI_M4x2y3z_S_b;
  Double I_NAI_M4xy4z_Py_b = I_NAI_N4x2y4z_S_b+ABY*I_NAI_M4xy4z_S_b;
  Double I_NAI_M3x6y_Py_b = I_NAI_N3x7y_S_b+ABY*I_NAI_M3x6y_S_b;
  Double I_NAI_M3x5yz_Py_b = I_NAI_N3x6yz_S_b+ABY*I_NAI_M3x5yz_S_b;
  Double I_NAI_M3x4y2z_Py_b = I_NAI_N3x5y2z_S_b+ABY*I_NAI_M3x4y2z_S_b;
  Double I_NAI_M3x3y3z_Py_b = I_NAI_N3x4y3z_S_b+ABY*I_NAI_M3x3y3z_S_b;
  Double I_NAI_M3x2y4z_Py_b = I_NAI_N3x3y4z_S_b+ABY*I_NAI_M3x2y4z_S_b;
  Double I_NAI_M3xy5z_Py_b = I_NAI_N3x2y5z_S_b+ABY*I_NAI_M3xy5z_S_b;
  Double I_NAI_M2x7y_Py_b = I_NAI_N2x8y_S_b+ABY*I_NAI_M2x7y_S_b;
  Double I_NAI_M2x6yz_Py_b = I_NAI_N2x7yz_S_b+ABY*I_NAI_M2x6yz_S_b;
  Double I_NAI_M2x5y2z_Py_b = I_NAI_N2x6y2z_S_b+ABY*I_NAI_M2x5y2z_S_b;
  Double I_NAI_M2x4y3z_Py_b = I_NAI_N2x5y3z_S_b+ABY*I_NAI_M2x4y3z_S_b;
  Double I_NAI_M2x3y4z_Py_b = I_NAI_N2x4y4z_S_b+ABY*I_NAI_M2x3y4z_S_b;
  Double I_NAI_M2x2y5z_Py_b = I_NAI_N2x3y5z_S_b+ABY*I_NAI_M2x2y5z_S_b;
  Double I_NAI_M2xy6z_Py_b = I_NAI_N2x2y6z_S_b+ABY*I_NAI_M2xy6z_S_b;
  Double I_NAI_Mx8y_Py_b = I_NAI_Nx9y_S_b+ABY*I_NAI_Mx8y_S_b;
  Double I_NAI_Mx7yz_Py_b = I_NAI_Nx8yz_S_b+ABY*I_NAI_Mx7yz_S_b;
  Double I_NAI_Mx6y2z_Py_b = I_NAI_Nx7y2z_S_b+ABY*I_NAI_Mx6y2z_S_b;
  Double I_NAI_Mx5y3z_Py_b = I_NAI_Nx6y3z_S_b+ABY*I_NAI_Mx5y3z_S_b;
  Double I_NAI_Mx4y4z_Py_b = I_NAI_Nx5y4z_S_b+ABY*I_NAI_Mx4y4z_S_b;
  Double I_NAI_Mx3y5z_Py_b = I_NAI_Nx4y5z_S_b+ABY*I_NAI_Mx3y5z_S_b;
  Double I_NAI_Mx2y6z_Py_b = I_NAI_Nx3y6z_S_b+ABY*I_NAI_Mx2y6z_S_b;
  Double I_NAI_Mxy7z_Py_b = I_NAI_Nx2y7z_S_b+ABY*I_NAI_Mxy7z_S_b;
  Double I_NAI_M9y_Py_b = I_NAI_N10y_S_b+ABY*I_NAI_M9y_S_b;
  Double I_NAI_M8yz_Py_b = I_NAI_N9yz_S_b+ABY*I_NAI_M8yz_S_b;
  Double I_NAI_M7y2z_Py_b = I_NAI_N8y2z_S_b+ABY*I_NAI_M7y2z_S_b;
  Double I_NAI_M6y3z_Py_b = I_NAI_N7y3z_S_b+ABY*I_NAI_M6y3z_S_b;
  Double I_NAI_M5y4z_Py_b = I_NAI_N6y4z_S_b+ABY*I_NAI_M5y4z_S_b;
  Double I_NAI_M4y5z_Py_b = I_NAI_N5y5z_S_b+ABY*I_NAI_M4y5z_S_b;
  Double I_NAI_M3y6z_Py_b = I_NAI_N4y6z_S_b+ABY*I_NAI_M3y6z_S_b;
  Double I_NAI_M2y7z_Py_b = I_NAI_N3y7z_S_b+ABY*I_NAI_M2y7z_S_b;
  Double I_NAI_M7x2z_Pz_b = I_NAI_N7x3z_S_b+ABZ*I_NAI_M7x2z_S_b;
  Double I_NAI_M6xy2z_Pz_b = I_NAI_N6xy3z_S_b+ABZ*I_NAI_M6xy2z_S_b;
  Double I_NAI_M6x3z_Pz_b = I_NAI_N6x4z_S_b+ABZ*I_NAI_M6x3z_S_b;
  Double I_NAI_M5x2y2z_Pz_b = I_NAI_N5x2y3z_S_b+ABZ*I_NAI_M5x2y2z_S_b;
  Double I_NAI_M5xy3z_Pz_b = I_NAI_N5xy4z_S_b+ABZ*I_NAI_M5xy3z_S_b;
  Double I_NAI_M5x4z_Pz_b = I_NAI_N5x5z_S_b+ABZ*I_NAI_M5x4z_S_b;
  Double I_NAI_M4x3y2z_Pz_b = I_NAI_N4x3y3z_S_b+ABZ*I_NAI_M4x3y2z_S_b;
  Double I_NAI_M4x2y3z_Pz_b = I_NAI_N4x2y4z_S_b+ABZ*I_NAI_M4x2y3z_S_b;
  Double I_NAI_M4xy4z_Pz_b = I_NAI_N4xy5z_S_b+ABZ*I_NAI_M4xy4z_S_b;
  Double I_NAI_M4x5z_Pz_b = I_NAI_N4x6z_S_b+ABZ*I_NAI_M4x5z_S_b;
  Double I_NAI_M3x4y2z_Pz_b = I_NAI_N3x4y3z_S_b+ABZ*I_NAI_M3x4y2z_S_b;
  Double I_NAI_M3x3y3z_Pz_b = I_NAI_N3x3y4z_S_b+ABZ*I_NAI_M3x3y3z_S_b;
  Double I_NAI_M3x2y4z_Pz_b = I_NAI_N3x2y5z_S_b+ABZ*I_NAI_M3x2y4z_S_b;
  Double I_NAI_M3xy5z_Pz_b = I_NAI_N3xy6z_S_b+ABZ*I_NAI_M3xy5z_S_b;
  Double I_NAI_M3x6z_Pz_b = I_NAI_N3x7z_S_b+ABZ*I_NAI_M3x6z_S_b;
  Double I_NAI_M2x5y2z_Pz_b = I_NAI_N2x5y3z_S_b+ABZ*I_NAI_M2x5y2z_S_b;
  Double I_NAI_M2x4y3z_Pz_b = I_NAI_N2x4y4z_S_b+ABZ*I_NAI_M2x4y3z_S_b;
  Double I_NAI_M2x3y4z_Pz_b = I_NAI_N2x3y5z_S_b+ABZ*I_NAI_M2x3y4z_S_b;
  Double I_NAI_M2x2y5z_Pz_b = I_NAI_N2x2y6z_S_b+ABZ*I_NAI_M2x2y5z_S_b;
  Double I_NAI_M2xy6z_Pz_b = I_NAI_N2xy7z_S_b+ABZ*I_NAI_M2xy6z_S_b;
  Double I_NAI_M2x7z_Pz_b = I_NAI_N2x8z_S_b+ABZ*I_NAI_M2x7z_S_b;
  Double I_NAI_Mx6y2z_Pz_b = I_NAI_Nx6y3z_S_b+ABZ*I_NAI_Mx6y2z_S_b;
  Double I_NAI_Mx5y3z_Pz_b = I_NAI_Nx5y4z_S_b+ABZ*I_NAI_Mx5y3z_S_b;
  Double I_NAI_Mx4y4z_Pz_b = I_NAI_Nx4y5z_S_b+ABZ*I_NAI_Mx4y4z_S_b;
  Double I_NAI_Mx3y5z_Pz_b = I_NAI_Nx3y6z_S_b+ABZ*I_NAI_Mx3y5z_S_b;
  Double I_NAI_Mx2y6z_Pz_b = I_NAI_Nx2y7z_S_b+ABZ*I_NAI_Mx2y6z_S_b;
  Double I_NAI_Mxy7z_Pz_b = I_NAI_Nxy8z_S_b+ABZ*I_NAI_Mxy7z_S_b;
  Double I_NAI_Mx8z_Pz_b = I_NAI_Nx9z_S_b+ABZ*I_NAI_Mx8z_S_b;
  Double I_NAI_M7y2z_Pz_b = I_NAI_N7y3z_S_b+ABZ*I_NAI_M7y2z_S_b;
  Double I_NAI_M6y3z_Pz_b = I_NAI_N6y4z_S_b+ABZ*I_NAI_M6y3z_S_b;
  Double I_NAI_M5y4z_Pz_b = I_NAI_N5y5z_S_b+ABZ*I_NAI_M5y4z_S_b;
  Double I_NAI_M4y5z_Pz_b = I_NAI_N4y6z_S_b+ABZ*I_NAI_M4y5z_S_b;
  Double I_NAI_M3y6z_Pz_b = I_NAI_N3y7z_S_b+ABZ*I_NAI_M3y6z_S_b;
  Double I_NAI_M2y7z_Pz_b = I_NAI_N2y8z_S_b+ABZ*I_NAI_M2y7z_S_b;
  Double I_NAI_My8z_Pz_b = I_NAI_Ny9z_S_b+ABZ*I_NAI_My8z_S_b;
  Double I_NAI_M9z_Pz_b = I_NAI_N10z_S_b+ABZ*I_NAI_M9z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_L_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 149 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_M_P_b
   * RHS shell quartet name: SQ_NAI_L_P_b
   ************************************************************/
  Double I_NAI_L8x_D2x_b = I_NAI_M9x_Px_b+ABX*I_NAI_L8x_Px_b;
  Double I_NAI_L7xy_D2x_b = I_NAI_M8xy_Px_b+ABX*I_NAI_L7xy_Px_b;
  Double I_NAI_L7xz_D2x_b = I_NAI_M8xz_Px_b+ABX*I_NAI_L7xz_Px_b;
  Double I_NAI_L6x2y_D2x_b = I_NAI_M7x2y_Px_b+ABX*I_NAI_L6x2y_Px_b;
  Double I_NAI_L6xyz_D2x_b = I_NAI_M7xyz_Px_b+ABX*I_NAI_L6xyz_Px_b;
  Double I_NAI_L6x2z_D2x_b = I_NAI_M7x2z_Px_b+ABX*I_NAI_L6x2z_Px_b;
  Double I_NAI_L5x3y_D2x_b = I_NAI_M6x3y_Px_b+ABX*I_NAI_L5x3y_Px_b;
  Double I_NAI_L5x2yz_D2x_b = I_NAI_M6x2yz_Px_b+ABX*I_NAI_L5x2yz_Px_b;
  Double I_NAI_L5xy2z_D2x_b = I_NAI_M6xy2z_Px_b+ABX*I_NAI_L5xy2z_Px_b;
  Double I_NAI_L5x3z_D2x_b = I_NAI_M6x3z_Px_b+ABX*I_NAI_L5x3z_Px_b;
  Double I_NAI_L4x4y_D2x_b = I_NAI_M5x4y_Px_b+ABX*I_NAI_L4x4y_Px_b;
  Double I_NAI_L4x3yz_D2x_b = I_NAI_M5x3yz_Px_b+ABX*I_NAI_L4x3yz_Px_b;
  Double I_NAI_L4x2y2z_D2x_b = I_NAI_M5x2y2z_Px_b+ABX*I_NAI_L4x2y2z_Px_b;
  Double I_NAI_L4xy3z_D2x_b = I_NAI_M5xy3z_Px_b+ABX*I_NAI_L4xy3z_Px_b;
  Double I_NAI_L4x4z_D2x_b = I_NAI_M5x4z_Px_b+ABX*I_NAI_L4x4z_Px_b;
  Double I_NAI_L3x5y_D2x_b = I_NAI_M4x5y_Px_b+ABX*I_NAI_L3x5y_Px_b;
  Double I_NAI_L3x4yz_D2x_b = I_NAI_M4x4yz_Px_b+ABX*I_NAI_L3x4yz_Px_b;
  Double I_NAI_L3x3y2z_D2x_b = I_NAI_M4x3y2z_Px_b+ABX*I_NAI_L3x3y2z_Px_b;
  Double I_NAI_L3x2y3z_D2x_b = I_NAI_M4x2y3z_Px_b+ABX*I_NAI_L3x2y3z_Px_b;
  Double I_NAI_L3xy4z_D2x_b = I_NAI_M4xy4z_Px_b+ABX*I_NAI_L3xy4z_Px_b;
  Double I_NAI_L3x5z_D2x_b = I_NAI_M4x5z_Px_b+ABX*I_NAI_L3x5z_Px_b;
  Double I_NAI_L2x6y_D2x_b = I_NAI_M3x6y_Px_b+ABX*I_NAI_L2x6y_Px_b;
  Double I_NAI_L2x5yz_D2x_b = I_NAI_M3x5yz_Px_b+ABX*I_NAI_L2x5yz_Px_b;
  Double I_NAI_L2x4y2z_D2x_b = I_NAI_M3x4y2z_Px_b+ABX*I_NAI_L2x4y2z_Px_b;
  Double I_NAI_L2x3y3z_D2x_b = I_NAI_M3x3y3z_Px_b+ABX*I_NAI_L2x3y3z_Px_b;
  Double I_NAI_L2x2y4z_D2x_b = I_NAI_M3x2y4z_Px_b+ABX*I_NAI_L2x2y4z_Px_b;
  Double I_NAI_L2xy5z_D2x_b = I_NAI_M3xy5z_Px_b+ABX*I_NAI_L2xy5z_Px_b;
  Double I_NAI_L2x6z_D2x_b = I_NAI_M3x6z_Px_b+ABX*I_NAI_L2x6z_Px_b;
  Double I_NAI_Lx7y_D2x_b = I_NAI_M2x7y_Px_b+ABX*I_NAI_Lx7y_Px_b;
  Double I_NAI_Lx6yz_D2x_b = I_NAI_M2x6yz_Px_b+ABX*I_NAI_Lx6yz_Px_b;
  Double I_NAI_Lx5y2z_D2x_b = I_NAI_M2x5y2z_Px_b+ABX*I_NAI_Lx5y2z_Px_b;
  Double I_NAI_Lx4y3z_D2x_b = I_NAI_M2x4y3z_Px_b+ABX*I_NAI_Lx4y3z_Px_b;
  Double I_NAI_Lx3y4z_D2x_b = I_NAI_M2x3y4z_Px_b+ABX*I_NAI_Lx3y4z_Px_b;
  Double I_NAI_Lx2y5z_D2x_b = I_NAI_M2x2y5z_Px_b+ABX*I_NAI_Lx2y5z_Px_b;
  Double I_NAI_Lxy6z_D2x_b = I_NAI_M2xy6z_Px_b+ABX*I_NAI_Lxy6z_Px_b;
  Double I_NAI_Lx7z_D2x_b = I_NAI_M2x7z_Px_b+ABX*I_NAI_Lx7z_Px_b;
  Double I_NAI_L7yz_D2x_b = I_NAI_Mx7yz_Px_b+ABX*I_NAI_L7yz_Px_b;
  Double I_NAI_L6y2z_D2x_b = I_NAI_Mx6y2z_Px_b+ABX*I_NAI_L6y2z_Px_b;
  Double I_NAI_L5y3z_D2x_b = I_NAI_Mx5y3z_Px_b+ABX*I_NAI_L5y3z_Px_b;
  Double I_NAI_L4y4z_D2x_b = I_NAI_Mx4y4z_Px_b+ABX*I_NAI_L4y4z_Px_b;
  Double I_NAI_L3y5z_D2x_b = I_NAI_Mx3y5z_Px_b+ABX*I_NAI_L3y5z_Px_b;
  Double I_NAI_L2y6z_D2x_b = I_NAI_Mx2y6z_Px_b+ABX*I_NAI_L2y6z_Px_b;
  Double I_NAI_Ly7z_D2x_b = I_NAI_Mxy7z_Px_b+ABX*I_NAI_Ly7z_Px_b;
  Double I_NAI_L7xy_D2y_b = I_NAI_M7x2y_Py_b+ABY*I_NAI_L7xy_Py_b;
  Double I_NAI_L6x2y_D2y_b = I_NAI_M6x3y_Py_b+ABY*I_NAI_L6x2y_Py_b;
  Double I_NAI_L6xyz_D2y_b = I_NAI_M6x2yz_Py_b+ABY*I_NAI_L6xyz_Py_b;
  Double I_NAI_L6x2z_D2y_b = I_NAI_M6xy2z_Py_b+ABY*I_NAI_L6x2z_Py_b;
  Double I_NAI_L5x3y_D2y_b = I_NAI_M5x4y_Py_b+ABY*I_NAI_L5x3y_Py_b;
  Double I_NAI_L5x2yz_D2y_b = I_NAI_M5x3yz_Py_b+ABY*I_NAI_L5x2yz_Py_b;
  Double I_NAI_L5xy2z_D2y_b = I_NAI_M5x2y2z_Py_b+ABY*I_NAI_L5xy2z_Py_b;
  Double I_NAI_L5x3z_D2y_b = I_NAI_M5xy3z_Py_b+ABY*I_NAI_L5x3z_Py_b;
  Double I_NAI_L4x4y_D2y_b = I_NAI_M4x5y_Py_b+ABY*I_NAI_L4x4y_Py_b;
  Double I_NAI_L4x3yz_D2y_b = I_NAI_M4x4yz_Py_b+ABY*I_NAI_L4x3yz_Py_b;
  Double I_NAI_L4x2y2z_D2y_b = I_NAI_M4x3y2z_Py_b+ABY*I_NAI_L4x2y2z_Py_b;
  Double I_NAI_L4xy3z_D2y_b = I_NAI_M4x2y3z_Py_b+ABY*I_NAI_L4xy3z_Py_b;
  Double I_NAI_L4x4z_D2y_b = I_NAI_M4xy4z_Py_b+ABY*I_NAI_L4x4z_Py_b;
  Double I_NAI_L3x5y_D2y_b = I_NAI_M3x6y_Py_b+ABY*I_NAI_L3x5y_Py_b;
  Double I_NAI_L3x4yz_D2y_b = I_NAI_M3x5yz_Py_b+ABY*I_NAI_L3x4yz_Py_b;
  Double I_NAI_L3x3y2z_D2y_b = I_NAI_M3x4y2z_Py_b+ABY*I_NAI_L3x3y2z_Py_b;
  Double I_NAI_L3x2y3z_D2y_b = I_NAI_M3x3y3z_Py_b+ABY*I_NAI_L3x2y3z_Py_b;
  Double I_NAI_L3xy4z_D2y_b = I_NAI_M3x2y4z_Py_b+ABY*I_NAI_L3xy4z_Py_b;
  Double I_NAI_L3x5z_D2y_b = I_NAI_M3xy5z_Py_b+ABY*I_NAI_L3x5z_Py_b;
  Double I_NAI_L2x6y_D2y_b = I_NAI_M2x7y_Py_b+ABY*I_NAI_L2x6y_Py_b;
  Double I_NAI_L2x5yz_D2y_b = I_NAI_M2x6yz_Py_b+ABY*I_NAI_L2x5yz_Py_b;
  Double I_NAI_L2x4y2z_D2y_b = I_NAI_M2x5y2z_Py_b+ABY*I_NAI_L2x4y2z_Py_b;
  Double I_NAI_L2x3y3z_D2y_b = I_NAI_M2x4y3z_Py_b+ABY*I_NAI_L2x3y3z_Py_b;
  Double I_NAI_L2x2y4z_D2y_b = I_NAI_M2x3y4z_Py_b+ABY*I_NAI_L2x2y4z_Py_b;
  Double I_NAI_L2xy5z_D2y_b = I_NAI_M2x2y5z_Py_b+ABY*I_NAI_L2xy5z_Py_b;
  Double I_NAI_L2x6z_D2y_b = I_NAI_M2xy6z_Py_b+ABY*I_NAI_L2x6z_Py_b;
  Double I_NAI_Lx7y_D2y_b = I_NAI_Mx8y_Py_b+ABY*I_NAI_Lx7y_Py_b;
  Double I_NAI_Lx6yz_D2y_b = I_NAI_Mx7yz_Py_b+ABY*I_NAI_Lx6yz_Py_b;
  Double I_NAI_Lx5y2z_D2y_b = I_NAI_Mx6y2z_Py_b+ABY*I_NAI_Lx5y2z_Py_b;
  Double I_NAI_Lx4y3z_D2y_b = I_NAI_Mx5y3z_Py_b+ABY*I_NAI_Lx4y3z_Py_b;
  Double I_NAI_Lx3y4z_D2y_b = I_NAI_Mx4y4z_Py_b+ABY*I_NAI_Lx3y4z_Py_b;
  Double I_NAI_Lx2y5z_D2y_b = I_NAI_Mx3y5z_Py_b+ABY*I_NAI_Lx2y5z_Py_b;
  Double I_NAI_Lxy6z_D2y_b = I_NAI_Mx2y6z_Py_b+ABY*I_NAI_Lxy6z_Py_b;
  Double I_NAI_Lx7z_D2y_b = I_NAI_Mxy7z_Py_b+ABY*I_NAI_Lx7z_Py_b;
  Double I_NAI_L8y_D2y_b = I_NAI_M9y_Py_b+ABY*I_NAI_L8y_Py_b;
  Double I_NAI_L7yz_D2y_b = I_NAI_M8yz_Py_b+ABY*I_NAI_L7yz_Py_b;
  Double I_NAI_L6y2z_D2y_b = I_NAI_M7y2z_Py_b+ABY*I_NAI_L6y2z_Py_b;
  Double I_NAI_L5y3z_D2y_b = I_NAI_M6y3z_Py_b+ABY*I_NAI_L5y3z_Py_b;
  Double I_NAI_L4y4z_D2y_b = I_NAI_M5y4z_Py_b+ABY*I_NAI_L4y4z_Py_b;
  Double I_NAI_L3y5z_D2y_b = I_NAI_M4y5z_Py_b+ABY*I_NAI_L3y5z_Py_b;
  Double I_NAI_L2y6z_D2y_b = I_NAI_M3y6z_Py_b+ABY*I_NAI_L2y6z_Py_b;
  Double I_NAI_Ly7z_D2y_b = I_NAI_M2y7z_Py_b+ABY*I_NAI_Ly7z_Py_b;
  Double I_NAI_L7xz_D2z_b = I_NAI_M7x2z_Pz_b+ABZ*I_NAI_L7xz_Pz_b;
  Double I_NAI_L6xyz_D2z_b = I_NAI_M6xy2z_Pz_b+ABZ*I_NAI_L6xyz_Pz_b;
  Double I_NAI_L6x2z_D2z_b = I_NAI_M6x3z_Pz_b+ABZ*I_NAI_L6x2z_Pz_b;
  Double I_NAI_L5x2yz_D2z_b = I_NAI_M5x2y2z_Pz_b+ABZ*I_NAI_L5x2yz_Pz_b;
  Double I_NAI_L5xy2z_D2z_b = I_NAI_M5xy3z_Pz_b+ABZ*I_NAI_L5xy2z_Pz_b;
  Double I_NAI_L5x3z_D2z_b = I_NAI_M5x4z_Pz_b+ABZ*I_NAI_L5x3z_Pz_b;
  Double I_NAI_L4x3yz_D2z_b = I_NAI_M4x3y2z_Pz_b+ABZ*I_NAI_L4x3yz_Pz_b;
  Double I_NAI_L4x2y2z_D2z_b = I_NAI_M4x2y3z_Pz_b+ABZ*I_NAI_L4x2y2z_Pz_b;
  Double I_NAI_L4xy3z_D2z_b = I_NAI_M4xy4z_Pz_b+ABZ*I_NAI_L4xy3z_Pz_b;
  Double I_NAI_L4x4z_D2z_b = I_NAI_M4x5z_Pz_b+ABZ*I_NAI_L4x4z_Pz_b;
  Double I_NAI_L3x4yz_D2z_b = I_NAI_M3x4y2z_Pz_b+ABZ*I_NAI_L3x4yz_Pz_b;
  Double I_NAI_L3x3y2z_D2z_b = I_NAI_M3x3y3z_Pz_b+ABZ*I_NAI_L3x3y2z_Pz_b;
  Double I_NAI_L3x2y3z_D2z_b = I_NAI_M3x2y4z_Pz_b+ABZ*I_NAI_L3x2y3z_Pz_b;
  Double I_NAI_L3xy4z_D2z_b = I_NAI_M3xy5z_Pz_b+ABZ*I_NAI_L3xy4z_Pz_b;
  Double I_NAI_L3x5z_D2z_b = I_NAI_M3x6z_Pz_b+ABZ*I_NAI_L3x5z_Pz_b;
  Double I_NAI_L2x5yz_D2z_b = I_NAI_M2x5y2z_Pz_b+ABZ*I_NAI_L2x5yz_Pz_b;
  Double I_NAI_L2x4y2z_D2z_b = I_NAI_M2x4y3z_Pz_b+ABZ*I_NAI_L2x4y2z_Pz_b;
  Double I_NAI_L2x3y3z_D2z_b = I_NAI_M2x3y4z_Pz_b+ABZ*I_NAI_L2x3y3z_Pz_b;
  Double I_NAI_L2x2y4z_D2z_b = I_NAI_M2x2y5z_Pz_b+ABZ*I_NAI_L2x2y4z_Pz_b;
  Double I_NAI_L2xy5z_D2z_b = I_NAI_M2xy6z_Pz_b+ABZ*I_NAI_L2xy5z_Pz_b;
  Double I_NAI_L2x6z_D2z_b = I_NAI_M2x7z_Pz_b+ABZ*I_NAI_L2x6z_Pz_b;
  Double I_NAI_Lx6yz_D2z_b = I_NAI_Mx6y2z_Pz_b+ABZ*I_NAI_Lx6yz_Pz_b;
  Double I_NAI_Lx5y2z_D2z_b = I_NAI_Mx5y3z_Pz_b+ABZ*I_NAI_Lx5y2z_Pz_b;
  Double I_NAI_Lx4y3z_D2z_b = I_NAI_Mx4y4z_Pz_b+ABZ*I_NAI_Lx4y3z_Pz_b;
  Double I_NAI_Lx3y4z_D2z_b = I_NAI_Mx3y5z_Pz_b+ABZ*I_NAI_Lx3y4z_Pz_b;
  Double I_NAI_Lx2y5z_D2z_b = I_NAI_Mx2y6z_Pz_b+ABZ*I_NAI_Lx2y5z_Pz_b;
  Double I_NAI_Lxy6z_D2z_b = I_NAI_Mxy7z_Pz_b+ABZ*I_NAI_Lxy6z_Pz_b;
  Double I_NAI_Lx7z_D2z_b = I_NAI_Mx8z_Pz_b+ABZ*I_NAI_Lx7z_Pz_b;
  Double I_NAI_L7yz_D2z_b = I_NAI_M7y2z_Pz_b+ABZ*I_NAI_L7yz_Pz_b;
  Double I_NAI_L6y2z_D2z_b = I_NAI_M6y3z_Pz_b+ABZ*I_NAI_L6y2z_Pz_b;
  Double I_NAI_L5y3z_D2z_b = I_NAI_M5y4z_Pz_b+ABZ*I_NAI_L5y3z_Pz_b;
  Double I_NAI_L4y4z_D2z_b = I_NAI_M4y5z_Pz_b+ABZ*I_NAI_L4y4z_Pz_b;
  Double I_NAI_L3y5z_D2z_b = I_NAI_M3y6z_Pz_b+ABZ*I_NAI_L3y5z_Pz_b;
  Double I_NAI_L2y6z_D2z_b = I_NAI_M2y7z_Pz_b+ABZ*I_NAI_L2y6z_Pz_b;
  Double I_NAI_Ly7z_D2z_b = I_NAI_My8z_Pz_b+ABZ*I_NAI_Ly7z_Pz_b;
  Double I_NAI_L8z_D2z_b = I_NAI_M9z_Pz_b+ABZ*I_NAI_L8z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_K_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 189 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_L_D_b
   * RHS shell quartet name: SQ_NAI_K_D_b
   ************************************************************/
  Double I_NAI_K7x_F3x_b = I_NAI_L8x_D2x_b+ABX*I_NAI_K7x_D2x_b;
  Double I_NAI_K6xy_F3x_b = I_NAI_L7xy_D2x_b+ABX*I_NAI_K6xy_D2x_b;
  Double I_NAI_K6xz_F3x_b = I_NAI_L7xz_D2x_b+ABX*I_NAI_K6xz_D2x_b;
  Double I_NAI_K5x2y_F3x_b = I_NAI_L6x2y_D2x_b+ABX*I_NAI_K5x2y_D2x_b;
  Double I_NAI_K5xyz_F3x_b = I_NAI_L6xyz_D2x_b+ABX*I_NAI_K5xyz_D2x_b;
  Double I_NAI_K5x2z_F3x_b = I_NAI_L6x2z_D2x_b+ABX*I_NAI_K5x2z_D2x_b;
  Double I_NAI_K4x3y_F3x_b = I_NAI_L5x3y_D2x_b+ABX*I_NAI_K4x3y_D2x_b;
  Double I_NAI_K4x2yz_F3x_b = I_NAI_L5x2yz_D2x_b+ABX*I_NAI_K4x2yz_D2x_b;
  Double I_NAI_K4xy2z_F3x_b = I_NAI_L5xy2z_D2x_b+ABX*I_NAI_K4xy2z_D2x_b;
  Double I_NAI_K4x3z_F3x_b = I_NAI_L5x3z_D2x_b+ABX*I_NAI_K4x3z_D2x_b;
  Double I_NAI_K3x4y_F3x_b = I_NAI_L4x4y_D2x_b+ABX*I_NAI_K3x4y_D2x_b;
  Double I_NAI_K3x3yz_F3x_b = I_NAI_L4x3yz_D2x_b+ABX*I_NAI_K3x3yz_D2x_b;
  Double I_NAI_K3x2y2z_F3x_b = I_NAI_L4x2y2z_D2x_b+ABX*I_NAI_K3x2y2z_D2x_b;
  Double I_NAI_K3xy3z_F3x_b = I_NAI_L4xy3z_D2x_b+ABX*I_NAI_K3xy3z_D2x_b;
  Double I_NAI_K3x4z_F3x_b = I_NAI_L4x4z_D2x_b+ABX*I_NAI_K3x4z_D2x_b;
  Double I_NAI_K2x5y_F3x_b = I_NAI_L3x5y_D2x_b+ABX*I_NAI_K2x5y_D2x_b;
  Double I_NAI_K2x4yz_F3x_b = I_NAI_L3x4yz_D2x_b+ABX*I_NAI_K2x4yz_D2x_b;
  Double I_NAI_K2x3y2z_F3x_b = I_NAI_L3x3y2z_D2x_b+ABX*I_NAI_K2x3y2z_D2x_b;
  Double I_NAI_K2x2y3z_F3x_b = I_NAI_L3x2y3z_D2x_b+ABX*I_NAI_K2x2y3z_D2x_b;
  Double I_NAI_K2xy4z_F3x_b = I_NAI_L3xy4z_D2x_b+ABX*I_NAI_K2xy4z_D2x_b;
  Double I_NAI_K2x5z_F3x_b = I_NAI_L3x5z_D2x_b+ABX*I_NAI_K2x5z_D2x_b;
  Double I_NAI_Kx6y_F3x_b = I_NAI_L2x6y_D2x_b+ABX*I_NAI_Kx6y_D2x_b;
  Double I_NAI_Kx5yz_F3x_b = I_NAI_L2x5yz_D2x_b+ABX*I_NAI_Kx5yz_D2x_b;
  Double I_NAI_Kx4y2z_F3x_b = I_NAI_L2x4y2z_D2x_b+ABX*I_NAI_Kx4y2z_D2x_b;
  Double I_NAI_Kx3y3z_F3x_b = I_NAI_L2x3y3z_D2x_b+ABX*I_NAI_Kx3y3z_D2x_b;
  Double I_NAI_Kx2y4z_F3x_b = I_NAI_L2x2y4z_D2x_b+ABX*I_NAI_Kx2y4z_D2x_b;
  Double I_NAI_Kxy5z_F3x_b = I_NAI_L2xy5z_D2x_b+ABX*I_NAI_Kxy5z_D2x_b;
  Double I_NAI_Kx6z_F3x_b = I_NAI_L2x6z_D2x_b+ABX*I_NAI_Kx6z_D2x_b;
  Double I_NAI_K7y_F3x_b = I_NAI_Lx7y_D2x_b+ABX*I_NAI_K7y_D2x_b;
  Double I_NAI_K6yz_F3x_b = I_NAI_Lx6yz_D2x_b+ABX*I_NAI_K6yz_D2x_b;
  Double I_NAI_K5y2z_F3x_b = I_NAI_Lx5y2z_D2x_b+ABX*I_NAI_K5y2z_D2x_b;
  Double I_NAI_K4y3z_F3x_b = I_NAI_Lx4y3z_D2x_b+ABX*I_NAI_K4y3z_D2x_b;
  Double I_NAI_K3y4z_F3x_b = I_NAI_Lx3y4z_D2x_b+ABX*I_NAI_K3y4z_D2x_b;
  Double I_NAI_K2y5z_F3x_b = I_NAI_Lx2y5z_D2x_b+ABX*I_NAI_K2y5z_D2x_b;
  Double I_NAI_Ky6z_F3x_b = I_NAI_Lxy6z_D2x_b+ABX*I_NAI_Ky6z_D2x_b;
  Double I_NAI_K7z_F3x_b = I_NAI_Lx7z_D2x_b+ABX*I_NAI_K7z_D2x_b;
  Double I_NAI_K5xyz_F2xy_b = I_NAI_L5x2yz_D2x_b+ABY*I_NAI_K5xyz_D2x_b;
  Double I_NAI_K4x2yz_F2xy_b = I_NAI_L4x3yz_D2x_b+ABY*I_NAI_K4x2yz_D2x_b;
  Double I_NAI_K4xy2z_F2xy_b = I_NAI_L4x2y2z_D2x_b+ABY*I_NAI_K4xy2z_D2x_b;
  Double I_NAI_K3x3yz_F2xy_b = I_NAI_L3x4yz_D2x_b+ABY*I_NAI_K3x3yz_D2x_b;
  Double I_NAI_K3x2y2z_F2xy_b = I_NAI_L3x3y2z_D2x_b+ABY*I_NAI_K3x2y2z_D2x_b;
  Double I_NAI_K3xy3z_F2xy_b = I_NAI_L3x2y3z_D2x_b+ABY*I_NAI_K3xy3z_D2x_b;
  Double I_NAI_K2x4yz_F2xy_b = I_NAI_L2x5yz_D2x_b+ABY*I_NAI_K2x4yz_D2x_b;
  Double I_NAI_K2x3y2z_F2xy_b = I_NAI_L2x4y2z_D2x_b+ABY*I_NAI_K2x3y2z_D2x_b;
  Double I_NAI_K2x2y3z_F2xy_b = I_NAI_L2x3y3z_D2x_b+ABY*I_NAI_K2x2y3z_D2x_b;
  Double I_NAI_K2xy4z_F2xy_b = I_NAI_L2x2y4z_D2x_b+ABY*I_NAI_K2xy4z_D2x_b;
  Double I_NAI_Kx5yz_F2xy_b = I_NAI_Lx6yz_D2x_b+ABY*I_NAI_Kx5yz_D2x_b;
  Double I_NAI_Kx4y2z_F2xy_b = I_NAI_Lx5y2z_D2x_b+ABY*I_NAI_Kx4y2z_D2x_b;
  Double I_NAI_Kx3y3z_F2xy_b = I_NAI_Lx4y3z_D2x_b+ABY*I_NAI_Kx3y3z_D2x_b;
  Double I_NAI_Kx2y4z_F2xy_b = I_NAI_Lx3y4z_D2x_b+ABY*I_NAI_Kx2y4z_D2x_b;
  Double I_NAI_Kxy5z_F2xy_b = I_NAI_Lx2y5z_D2x_b+ABY*I_NAI_Kxy5z_D2x_b;
  Double I_NAI_K6yz_F2xy_b = I_NAI_L7yz_D2x_b+ABY*I_NAI_K6yz_D2x_b;
  Double I_NAI_K5y2z_F2xy_b = I_NAI_L6y2z_D2x_b+ABY*I_NAI_K5y2z_D2x_b;
  Double I_NAI_K4y3z_F2xy_b = I_NAI_L5y3z_D2x_b+ABY*I_NAI_K4y3z_D2x_b;
  Double I_NAI_K3y4z_F2xy_b = I_NAI_L4y4z_D2x_b+ABY*I_NAI_K3y4z_D2x_b;
  Double I_NAI_K2y5z_F2xy_b = I_NAI_L3y5z_D2x_b+ABY*I_NAI_K2y5z_D2x_b;
  Double I_NAI_Ky6z_F2xy_b = I_NAI_L2y6z_D2x_b+ABY*I_NAI_Ky6z_D2x_b;
  Double I_NAI_K5xyz_F2xz_b = I_NAI_L5xy2z_D2x_b+ABZ*I_NAI_K5xyz_D2x_b;
  Double I_NAI_K4x2yz_F2xz_b = I_NAI_L4x2y2z_D2x_b+ABZ*I_NAI_K4x2yz_D2x_b;
  Double I_NAI_K4xy2z_F2xz_b = I_NAI_L4xy3z_D2x_b+ABZ*I_NAI_K4xy2z_D2x_b;
  Double I_NAI_K3x3yz_F2xz_b = I_NAI_L3x3y2z_D2x_b+ABZ*I_NAI_K3x3yz_D2x_b;
  Double I_NAI_K3x2y2z_F2xz_b = I_NAI_L3x2y3z_D2x_b+ABZ*I_NAI_K3x2y2z_D2x_b;
  Double I_NAI_K3xy3z_F2xz_b = I_NAI_L3xy4z_D2x_b+ABZ*I_NAI_K3xy3z_D2x_b;
  Double I_NAI_K2x4yz_F2xz_b = I_NAI_L2x4y2z_D2x_b+ABZ*I_NAI_K2x4yz_D2x_b;
  Double I_NAI_K2x3y2z_F2xz_b = I_NAI_L2x3y3z_D2x_b+ABZ*I_NAI_K2x3y2z_D2x_b;
  Double I_NAI_K2x2y3z_F2xz_b = I_NAI_L2x2y4z_D2x_b+ABZ*I_NAI_K2x2y3z_D2x_b;
  Double I_NAI_K2xy4z_F2xz_b = I_NAI_L2xy5z_D2x_b+ABZ*I_NAI_K2xy4z_D2x_b;
  Double I_NAI_Kx5yz_F2xz_b = I_NAI_Lx5y2z_D2x_b+ABZ*I_NAI_Kx5yz_D2x_b;
  Double I_NAI_Kx4y2z_F2xz_b = I_NAI_Lx4y3z_D2x_b+ABZ*I_NAI_Kx4y2z_D2x_b;
  Double I_NAI_Kx3y3z_F2xz_b = I_NAI_Lx3y4z_D2x_b+ABZ*I_NAI_Kx3y3z_D2x_b;
  Double I_NAI_Kx2y4z_F2xz_b = I_NAI_Lx2y5z_D2x_b+ABZ*I_NAI_Kx2y4z_D2x_b;
  Double I_NAI_Kxy5z_F2xz_b = I_NAI_Lxy6z_D2x_b+ABZ*I_NAI_Kxy5z_D2x_b;
  Double I_NAI_K6yz_F2xz_b = I_NAI_L6y2z_D2x_b+ABZ*I_NAI_K6yz_D2x_b;
  Double I_NAI_K5y2z_F2xz_b = I_NAI_L5y3z_D2x_b+ABZ*I_NAI_K5y2z_D2x_b;
  Double I_NAI_K4y3z_F2xz_b = I_NAI_L4y4z_D2x_b+ABZ*I_NAI_K4y3z_D2x_b;
  Double I_NAI_K3y4z_F2xz_b = I_NAI_L3y5z_D2x_b+ABZ*I_NAI_K3y4z_D2x_b;
  Double I_NAI_K2y5z_F2xz_b = I_NAI_L2y6z_D2x_b+ABZ*I_NAI_K2y5z_D2x_b;
  Double I_NAI_Ky6z_F2xz_b = I_NAI_Ly7z_D2x_b+ABZ*I_NAI_Ky6z_D2x_b;
  Double I_NAI_K7x_F3y_b = I_NAI_L7xy_D2y_b+ABY*I_NAI_K7x_D2y_b;
  Double I_NAI_K6xy_F3y_b = I_NAI_L6x2y_D2y_b+ABY*I_NAI_K6xy_D2y_b;
  Double I_NAI_K6xz_F3y_b = I_NAI_L6xyz_D2y_b+ABY*I_NAI_K6xz_D2y_b;
  Double I_NAI_K5x2y_F3y_b = I_NAI_L5x3y_D2y_b+ABY*I_NAI_K5x2y_D2y_b;
  Double I_NAI_K5xyz_F3y_b = I_NAI_L5x2yz_D2y_b+ABY*I_NAI_K5xyz_D2y_b;
  Double I_NAI_K5x2z_F3y_b = I_NAI_L5xy2z_D2y_b+ABY*I_NAI_K5x2z_D2y_b;
  Double I_NAI_K4x3y_F3y_b = I_NAI_L4x4y_D2y_b+ABY*I_NAI_K4x3y_D2y_b;
  Double I_NAI_K4x2yz_F3y_b = I_NAI_L4x3yz_D2y_b+ABY*I_NAI_K4x2yz_D2y_b;
  Double I_NAI_K4xy2z_F3y_b = I_NAI_L4x2y2z_D2y_b+ABY*I_NAI_K4xy2z_D2y_b;
  Double I_NAI_K4x3z_F3y_b = I_NAI_L4xy3z_D2y_b+ABY*I_NAI_K4x3z_D2y_b;
  Double I_NAI_K3x4y_F3y_b = I_NAI_L3x5y_D2y_b+ABY*I_NAI_K3x4y_D2y_b;
  Double I_NAI_K3x3yz_F3y_b = I_NAI_L3x4yz_D2y_b+ABY*I_NAI_K3x3yz_D2y_b;
  Double I_NAI_K3x2y2z_F3y_b = I_NAI_L3x3y2z_D2y_b+ABY*I_NAI_K3x2y2z_D2y_b;
  Double I_NAI_K3xy3z_F3y_b = I_NAI_L3x2y3z_D2y_b+ABY*I_NAI_K3xy3z_D2y_b;
  Double I_NAI_K3x4z_F3y_b = I_NAI_L3xy4z_D2y_b+ABY*I_NAI_K3x4z_D2y_b;
  Double I_NAI_K2x5y_F3y_b = I_NAI_L2x6y_D2y_b+ABY*I_NAI_K2x5y_D2y_b;
  Double I_NAI_K2x4yz_F3y_b = I_NAI_L2x5yz_D2y_b+ABY*I_NAI_K2x4yz_D2y_b;
  Double I_NAI_K2x3y2z_F3y_b = I_NAI_L2x4y2z_D2y_b+ABY*I_NAI_K2x3y2z_D2y_b;
  Double I_NAI_K2x2y3z_F3y_b = I_NAI_L2x3y3z_D2y_b+ABY*I_NAI_K2x2y3z_D2y_b;
  Double I_NAI_K2xy4z_F3y_b = I_NAI_L2x2y4z_D2y_b+ABY*I_NAI_K2xy4z_D2y_b;
  Double I_NAI_K2x5z_F3y_b = I_NAI_L2xy5z_D2y_b+ABY*I_NAI_K2x5z_D2y_b;
  Double I_NAI_Kx6y_F3y_b = I_NAI_Lx7y_D2y_b+ABY*I_NAI_Kx6y_D2y_b;
  Double I_NAI_Kx5yz_F3y_b = I_NAI_Lx6yz_D2y_b+ABY*I_NAI_Kx5yz_D2y_b;
  Double I_NAI_Kx4y2z_F3y_b = I_NAI_Lx5y2z_D2y_b+ABY*I_NAI_Kx4y2z_D2y_b;
  Double I_NAI_Kx3y3z_F3y_b = I_NAI_Lx4y3z_D2y_b+ABY*I_NAI_Kx3y3z_D2y_b;
  Double I_NAI_Kx2y4z_F3y_b = I_NAI_Lx3y4z_D2y_b+ABY*I_NAI_Kx2y4z_D2y_b;
  Double I_NAI_Kxy5z_F3y_b = I_NAI_Lx2y5z_D2y_b+ABY*I_NAI_Kxy5z_D2y_b;
  Double I_NAI_Kx6z_F3y_b = I_NAI_Lxy6z_D2y_b+ABY*I_NAI_Kx6z_D2y_b;
  Double I_NAI_K7y_F3y_b = I_NAI_L8y_D2y_b+ABY*I_NAI_K7y_D2y_b;
  Double I_NAI_K6yz_F3y_b = I_NAI_L7yz_D2y_b+ABY*I_NAI_K6yz_D2y_b;
  Double I_NAI_K5y2z_F3y_b = I_NAI_L6y2z_D2y_b+ABY*I_NAI_K5y2z_D2y_b;
  Double I_NAI_K4y3z_F3y_b = I_NAI_L5y3z_D2y_b+ABY*I_NAI_K4y3z_D2y_b;
  Double I_NAI_K3y4z_F3y_b = I_NAI_L4y4z_D2y_b+ABY*I_NAI_K3y4z_D2y_b;
  Double I_NAI_K2y5z_F3y_b = I_NAI_L3y5z_D2y_b+ABY*I_NAI_K2y5z_D2y_b;
  Double I_NAI_Ky6z_F3y_b = I_NAI_L2y6z_D2y_b+ABY*I_NAI_Ky6z_D2y_b;
  Double I_NAI_K7z_F3y_b = I_NAI_Ly7z_D2y_b+ABY*I_NAI_K7z_D2y_b;
  Double I_NAI_K6xz_F2yz_b = I_NAI_L6x2z_D2y_b+ABZ*I_NAI_K6xz_D2y_b;
  Double I_NAI_K5xyz_F2yz_b = I_NAI_L5xy2z_D2y_b+ABZ*I_NAI_K5xyz_D2y_b;
  Double I_NAI_K5x2z_F2yz_b = I_NAI_L5x3z_D2y_b+ABZ*I_NAI_K5x2z_D2y_b;
  Double I_NAI_K4x2yz_F2yz_b = I_NAI_L4x2y2z_D2y_b+ABZ*I_NAI_K4x2yz_D2y_b;
  Double I_NAI_K4xy2z_F2yz_b = I_NAI_L4xy3z_D2y_b+ABZ*I_NAI_K4xy2z_D2y_b;
  Double I_NAI_K4x3z_F2yz_b = I_NAI_L4x4z_D2y_b+ABZ*I_NAI_K4x3z_D2y_b;
  Double I_NAI_K3x3yz_F2yz_b = I_NAI_L3x3y2z_D2y_b+ABZ*I_NAI_K3x3yz_D2y_b;
  Double I_NAI_K3x2y2z_F2yz_b = I_NAI_L3x2y3z_D2y_b+ABZ*I_NAI_K3x2y2z_D2y_b;
  Double I_NAI_K3xy3z_F2yz_b = I_NAI_L3xy4z_D2y_b+ABZ*I_NAI_K3xy3z_D2y_b;
  Double I_NAI_K3x4z_F2yz_b = I_NAI_L3x5z_D2y_b+ABZ*I_NAI_K3x4z_D2y_b;
  Double I_NAI_K2x4yz_F2yz_b = I_NAI_L2x4y2z_D2y_b+ABZ*I_NAI_K2x4yz_D2y_b;
  Double I_NAI_K2x3y2z_F2yz_b = I_NAI_L2x3y3z_D2y_b+ABZ*I_NAI_K2x3y2z_D2y_b;
  Double I_NAI_K2x2y3z_F2yz_b = I_NAI_L2x2y4z_D2y_b+ABZ*I_NAI_K2x2y3z_D2y_b;
  Double I_NAI_K2xy4z_F2yz_b = I_NAI_L2xy5z_D2y_b+ABZ*I_NAI_K2xy4z_D2y_b;
  Double I_NAI_K2x5z_F2yz_b = I_NAI_L2x6z_D2y_b+ABZ*I_NAI_K2x5z_D2y_b;
  Double I_NAI_Kx5yz_F2yz_b = I_NAI_Lx5y2z_D2y_b+ABZ*I_NAI_Kx5yz_D2y_b;
  Double I_NAI_Kx4y2z_F2yz_b = I_NAI_Lx4y3z_D2y_b+ABZ*I_NAI_Kx4y2z_D2y_b;
  Double I_NAI_Kx3y3z_F2yz_b = I_NAI_Lx3y4z_D2y_b+ABZ*I_NAI_Kx3y3z_D2y_b;
  Double I_NAI_Kx2y4z_F2yz_b = I_NAI_Lx2y5z_D2y_b+ABZ*I_NAI_Kx2y4z_D2y_b;
  Double I_NAI_Kxy5z_F2yz_b = I_NAI_Lxy6z_D2y_b+ABZ*I_NAI_Kxy5z_D2y_b;
  Double I_NAI_Kx6z_F2yz_b = I_NAI_Lx7z_D2y_b+ABZ*I_NAI_Kx6z_D2y_b;
  Double I_NAI_K7x_F3z_b = I_NAI_L7xz_D2z_b+ABZ*I_NAI_K7x_D2z_b;
  Double I_NAI_K6xy_F3z_b = I_NAI_L6xyz_D2z_b+ABZ*I_NAI_K6xy_D2z_b;
  Double I_NAI_K6xz_F3z_b = I_NAI_L6x2z_D2z_b+ABZ*I_NAI_K6xz_D2z_b;
  Double I_NAI_K5x2y_F3z_b = I_NAI_L5x2yz_D2z_b+ABZ*I_NAI_K5x2y_D2z_b;
  Double I_NAI_K5xyz_F3z_b = I_NAI_L5xy2z_D2z_b+ABZ*I_NAI_K5xyz_D2z_b;
  Double I_NAI_K5x2z_F3z_b = I_NAI_L5x3z_D2z_b+ABZ*I_NAI_K5x2z_D2z_b;
  Double I_NAI_K4x3y_F3z_b = I_NAI_L4x3yz_D2z_b+ABZ*I_NAI_K4x3y_D2z_b;
  Double I_NAI_K4x2yz_F3z_b = I_NAI_L4x2y2z_D2z_b+ABZ*I_NAI_K4x2yz_D2z_b;
  Double I_NAI_K4xy2z_F3z_b = I_NAI_L4xy3z_D2z_b+ABZ*I_NAI_K4xy2z_D2z_b;
  Double I_NAI_K4x3z_F3z_b = I_NAI_L4x4z_D2z_b+ABZ*I_NAI_K4x3z_D2z_b;
  Double I_NAI_K3x4y_F3z_b = I_NAI_L3x4yz_D2z_b+ABZ*I_NAI_K3x4y_D2z_b;
  Double I_NAI_K3x3yz_F3z_b = I_NAI_L3x3y2z_D2z_b+ABZ*I_NAI_K3x3yz_D2z_b;
  Double I_NAI_K3x2y2z_F3z_b = I_NAI_L3x2y3z_D2z_b+ABZ*I_NAI_K3x2y2z_D2z_b;
  Double I_NAI_K3xy3z_F3z_b = I_NAI_L3xy4z_D2z_b+ABZ*I_NAI_K3xy3z_D2z_b;
  Double I_NAI_K3x4z_F3z_b = I_NAI_L3x5z_D2z_b+ABZ*I_NAI_K3x4z_D2z_b;
  Double I_NAI_K2x5y_F3z_b = I_NAI_L2x5yz_D2z_b+ABZ*I_NAI_K2x5y_D2z_b;
  Double I_NAI_K2x4yz_F3z_b = I_NAI_L2x4y2z_D2z_b+ABZ*I_NAI_K2x4yz_D2z_b;
  Double I_NAI_K2x3y2z_F3z_b = I_NAI_L2x3y3z_D2z_b+ABZ*I_NAI_K2x3y2z_D2z_b;
  Double I_NAI_K2x2y3z_F3z_b = I_NAI_L2x2y4z_D2z_b+ABZ*I_NAI_K2x2y3z_D2z_b;
  Double I_NAI_K2xy4z_F3z_b = I_NAI_L2xy5z_D2z_b+ABZ*I_NAI_K2xy4z_D2z_b;
  Double I_NAI_K2x5z_F3z_b = I_NAI_L2x6z_D2z_b+ABZ*I_NAI_K2x5z_D2z_b;
  Double I_NAI_Kx6y_F3z_b = I_NAI_Lx6yz_D2z_b+ABZ*I_NAI_Kx6y_D2z_b;
  Double I_NAI_Kx5yz_F3z_b = I_NAI_Lx5y2z_D2z_b+ABZ*I_NAI_Kx5yz_D2z_b;
  Double I_NAI_Kx4y2z_F3z_b = I_NAI_Lx4y3z_D2z_b+ABZ*I_NAI_Kx4y2z_D2z_b;
  Double I_NAI_Kx3y3z_F3z_b = I_NAI_Lx3y4z_D2z_b+ABZ*I_NAI_Kx3y3z_D2z_b;
  Double I_NAI_Kx2y4z_F3z_b = I_NAI_Lx2y5z_D2z_b+ABZ*I_NAI_Kx2y4z_D2z_b;
  Double I_NAI_Kxy5z_F3z_b = I_NAI_Lxy6z_D2z_b+ABZ*I_NAI_Kxy5z_D2z_b;
  Double I_NAI_Kx6z_F3z_b = I_NAI_Lx7z_D2z_b+ABZ*I_NAI_Kx6z_D2z_b;
  Double I_NAI_K7y_F3z_b = I_NAI_L7yz_D2z_b+ABZ*I_NAI_K7y_D2z_b;
  Double I_NAI_K6yz_F3z_b = I_NAI_L6y2z_D2z_b+ABZ*I_NAI_K6yz_D2z_b;
  Double I_NAI_K5y2z_F3z_b = I_NAI_L5y3z_D2z_b+ABZ*I_NAI_K5y2z_D2z_b;
  Double I_NAI_K4y3z_F3z_b = I_NAI_L4y4z_D2z_b+ABZ*I_NAI_K4y3z_D2z_b;
  Double I_NAI_K3y4z_F3z_b = I_NAI_L3y5z_D2z_b+ABZ*I_NAI_K3y4z_D2z_b;
  Double I_NAI_K2y5z_F3z_b = I_NAI_L2y6z_D2z_b+ABZ*I_NAI_K2y5z_D2z_b;
  Double I_NAI_Ky6z_F3z_b = I_NAI_Ly7z_D2z_b+ABZ*I_NAI_Ky6z_D2z_b;
  Double I_NAI_K7z_F3z_b = I_NAI_L8z_D2z_b+ABZ*I_NAI_K7z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_G_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 129 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_F_b
   * RHS shell quartet name: SQ_NAI_I_F_b
   ************************************************************/
  Double I_NAI_I6x_G4x_b = I_NAI_K7x_F3x_b+ABX*I_NAI_I6x_F3x_b;
  Double I_NAI_I5xy_G4x_b = I_NAI_K6xy_F3x_b+ABX*I_NAI_I5xy_F3x_b;
  Double I_NAI_I5xz_G4x_b = I_NAI_K6xz_F3x_b+ABX*I_NAI_I5xz_F3x_b;
  Double I_NAI_I4x2y_G4x_b = I_NAI_K5x2y_F3x_b+ABX*I_NAI_I4x2y_F3x_b;
  Double I_NAI_I4xyz_G4x_b = I_NAI_K5xyz_F3x_b+ABX*I_NAI_I4xyz_F3x_b;
  Double I_NAI_I4x2z_G4x_b = I_NAI_K5x2z_F3x_b+ABX*I_NAI_I4x2z_F3x_b;
  Double I_NAI_I3x3y_G4x_b = I_NAI_K4x3y_F3x_b+ABX*I_NAI_I3x3y_F3x_b;
  Double I_NAI_I3x2yz_G4x_b = I_NAI_K4x2yz_F3x_b+ABX*I_NAI_I3x2yz_F3x_b;
  Double I_NAI_I3xy2z_G4x_b = I_NAI_K4xy2z_F3x_b+ABX*I_NAI_I3xy2z_F3x_b;
  Double I_NAI_I3x3z_G4x_b = I_NAI_K4x3z_F3x_b+ABX*I_NAI_I3x3z_F3x_b;
  Double I_NAI_I2x4y_G4x_b = I_NAI_K3x4y_F3x_b+ABX*I_NAI_I2x4y_F3x_b;
  Double I_NAI_I2x3yz_G4x_b = I_NAI_K3x3yz_F3x_b+ABX*I_NAI_I2x3yz_F3x_b;
  Double I_NAI_I2x2y2z_G4x_b = I_NAI_K3x2y2z_F3x_b+ABX*I_NAI_I2x2y2z_F3x_b;
  Double I_NAI_I2xy3z_G4x_b = I_NAI_K3xy3z_F3x_b+ABX*I_NAI_I2xy3z_F3x_b;
  Double I_NAI_I2x4z_G4x_b = I_NAI_K3x4z_F3x_b+ABX*I_NAI_I2x4z_F3x_b;
  Double I_NAI_Ix5y_G4x_b = I_NAI_K2x5y_F3x_b+ABX*I_NAI_Ix5y_F3x_b;
  Double I_NAI_Ix4yz_G4x_b = I_NAI_K2x4yz_F3x_b+ABX*I_NAI_Ix4yz_F3x_b;
  Double I_NAI_Ix3y2z_G4x_b = I_NAI_K2x3y2z_F3x_b+ABX*I_NAI_Ix3y2z_F3x_b;
  Double I_NAI_Ix2y3z_G4x_b = I_NAI_K2x2y3z_F3x_b+ABX*I_NAI_Ix2y3z_F3x_b;
  Double I_NAI_Ixy4z_G4x_b = I_NAI_K2xy4z_F3x_b+ABX*I_NAI_Ixy4z_F3x_b;
  Double I_NAI_Ix5z_G4x_b = I_NAI_K2x5z_F3x_b+ABX*I_NAI_Ix5z_F3x_b;
  Double I_NAI_I6y_G4x_b = I_NAI_Kx6y_F3x_b+ABX*I_NAI_I6y_F3x_b;
  Double I_NAI_I5yz_G4x_b = I_NAI_Kx5yz_F3x_b+ABX*I_NAI_I5yz_F3x_b;
  Double I_NAI_I4y2z_G4x_b = I_NAI_Kx4y2z_F3x_b+ABX*I_NAI_I4y2z_F3x_b;
  Double I_NAI_I3y3z_G4x_b = I_NAI_Kx3y3z_F3x_b+ABX*I_NAI_I3y3z_F3x_b;
  Double I_NAI_I2y4z_G4x_b = I_NAI_Kx2y4z_F3x_b+ABX*I_NAI_I2y4z_F3x_b;
  Double I_NAI_Iy5z_G4x_b = I_NAI_Kxy5z_F3x_b+ABX*I_NAI_Iy5z_F3x_b;
  Double I_NAI_I6z_G4x_b = I_NAI_Kx6z_F3x_b+ABX*I_NAI_I6z_F3x_b;
  Double I_NAI_I5xy_G3xy_b = I_NAI_K5x2y_F3x_b+ABY*I_NAI_I5xy_F3x_b;
  Double I_NAI_I5xz_G3xy_b = I_NAI_K5xyz_F3x_b+ABY*I_NAI_I5xz_F3x_b;
  Double I_NAI_I4x2y_G3xy_b = I_NAI_K4x3y_F3x_b+ABY*I_NAI_I4x2y_F3x_b;
  Double I_NAI_I4xyz_G3xy_b = I_NAI_K4x2yz_F3x_b+ABY*I_NAI_I4xyz_F3x_b;
  Double I_NAI_I4x2z_G3xy_b = I_NAI_K4xy2z_F3x_b+ABY*I_NAI_I4x2z_F3x_b;
  Double I_NAI_I3x3y_G3xy_b = I_NAI_K3x4y_F3x_b+ABY*I_NAI_I3x3y_F3x_b;
  Double I_NAI_I3x2yz_G3xy_b = I_NAI_K3x3yz_F3x_b+ABY*I_NAI_I3x2yz_F3x_b;
  Double I_NAI_I3xy2z_G3xy_b = I_NAI_K3x2y2z_F3x_b+ABY*I_NAI_I3xy2z_F3x_b;
  Double I_NAI_I3x3z_G3xy_b = I_NAI_K3xy3z_F3x_b+ABY*I_NAI_I3x3z_F3x_b;
  Double I_NAI_I2x4y_G3xy_b = I_NAI_K2x5y_F3x_b+ABY*I_NAI_I2x4y_F3x_b;
  Double I_NAI_I2x3yz_G3xy_b = I_NAI_K2x4yz_F3x_b+ABY*I_NAI_I2x3yz_F3x_b;
  Double I_NAI_I2x2y2z_G3xy_b = I_NAI_K2x3y2z_F3x_b+ABY*I_NAI_I2x2y2z_F3x_b;
  Double I_NAI_I2xy3z_G3xy_b = I_NAI_K2x2y3z_F3x_b+ABY*I_NAI_I2xy3z_F3x_b;
  Double I_NAI_I2x4z_G3xy_b = I_NAI_K2xy4z_F3x_b+ABY*I_NAI_I2x4z_F3x_b;
  Double I_NAI_Ix5y_G3xy_b = I_NAI_Kx6y_F3x_b+ABY*I_NAI_Ix5y_F3x_b;
  Double I_NAI_Ix4yz_G3xy_b = I_NAI_Kx5yz_F3x_b+ABY*I_NAI_Ix4yz_F3x_b;
  Double I_NAI_Ix3y2z_G3xy_b = I_NAI_Kx4y2z_F3x_b+ABY*I_NAI_Ix3y2z_F3x_b;
  Double I_NAI_Ix2y3z_G3xy_b = I_NAI_Kx3y3z_F3x_b+ABY*I_NAI_Ix2y3z_F3x_b;
  Double I_NAI_Ixy4z_G3xy_b = I_NAI_Kx2y4z_F3x_b+ABY*I_NAI_Ixy4z_F3x_b;
  Double I_NAI_Ix5z_G3xy_b = I_NAI_Kxy5z_F3x_b+ABY*I_NAI_Ix5z_F3x_b;
  Double I_NAI_I6y_G3xy_b = I_NAI_K7y_F3x_b+ABY*I_NAI_I6y_F3x_b;
  Double I_NAI_I5yz_G3xy_b = I_NAI_K6yz_F3x_b+ABY*I_NAI_I5yz_F3x_b;
  Double I_NAI_I4y2z_G3xy_b = I_NAI_K5y2z_F3x_b+ABY*I_NAI_I4y2z_F3x_b;
  Double I_NAI_I3y3z_G3xy_b = I_NAI_K4y3z_F3x_b+ABY*I_NAI_I3y3z_F3x_b;
  Double I_NAI_I2y4z_G3xy_b = I_NAI_K3y4z_F3x_b+ABY*I_NAI_I2y4z_F3x_b;
  Double I_NAI_Iy5z_G3xy_b = I_NAI_K2y5z_F3x_b+ABY*I_NAI_Iy5z_F3x_b;
  Double I_NAI_I6z_G3xy_b = I_NAI_Ky6z_F3x_b+ABY*I_NAI_I6z_F3x_b;
  Double I_NAI_I5xz_G3xz_b = I_NAI_K5x2z_F3x_b+ABZ*I_NAI_I5xz_F3x_b;
  Double I_NAI_I4xyz_G3xz_b = I_NAI_K4xy2z_F3x_b+ABZ*I_NAI_I4xyz_F3x_b;
  Double I_NAI_I4x2z_G3xz_b = I_NAI_K4x3z_F3x_b+ABZ*I_NAI_I4x2z_F3x_b;
  Double I_NAI_I3x2yz_G3xz_b = I_NAI_K3x2y2z_F3x_b+ABZ*I_NAI_I3x2yz_F3x_b;
  Double I_NAI_I3xy2z_G3xz_b = I_NAI_K3xy3z_F3x_b+ABZ*I_NAI_I3xy2z_F3x_b;
  Double I_NAI_I3x3z_G3xz_b = I_NAI_K3x4z_F3x_b+ABZ*I_NAI_I3x3z_F3x_b;
  Double I_NAI_I2x3yz_G3xz_b = I_NAI_K2x3y2z_F3x_b+ABZ*I_NAI_I2x3yz_F3x_b;
  Double I_NAI_I2x2y2z_G3xz_b = I_NAI_K2x2y3z_F3x_b+ABZ*I_NAI_I2x2y2z_F3x_b;
  Double I_NAI_I2xy3z_G3xz_b = I_NAI_K2xy4z_F3x_b+ABZ*I_NAI_I2xy3z_F3x_b;
  Double I_NAI_I2x4z_G3xz_b = I_NAI_K2x5z_F3x_b+ABZ*I_NAI_I2x4z_F3x_b;
  Double I_NAI_Ix4yz_G3xz_b = I_NAI_Kx4y2z_F3x_b+ABZ*I_NAI_Ix4yz_F3x_b;
  Double I_NAI_Ix3y2z_G3xz_b = I_NAI_Kx3y3z_F3x_b+ABZ*I_NAI_Ix3y2z_F3x_b;
  Double I_NAI_Ix2y3z_G3xz_b = I_NAI_Kx2y4z_F3x_b+ABZ*I_NAI_Ix2y3z_F3x_b;
  Double I_NAI_Ixy4z_G3xz_b = I_NAI_Kxy5z_F3x_b+ABZ*I_NAI_Ixy4z_F3x_b;
  Double I_NAI_Ix5z_G3xz_b = I_NAI_Kx6z_F3x_b+ABZ*I_NAI_Ix5z_F3x_b;
  Double I_NAI_I5yz_G3xz_b = I_NAI_K5y2z_F3x_b+ABZ*I_NAI_I5yz_F3x_b;
  Double I_NAI_I4y2z_G3xz_b = I_NAI_K4y3z_F3x_b+ABZ*I_NAI_I4y2z_F3x_b;
  Double I_NAI_I3y3z_G3xz_b = I_NAI_K3y4z_F3x_b+ABZ*I_NAI_I3y3z_F3x_b;
  Double I_NAI_I2y4z_G3xz_b = I_NAI_K2y5z_F3x_b+ABZ*I_NAI_I2y4z_F3x_b;
  Double I_NAI_Iy5z_G3xz_b = I_NAI_Ky6z_F3x_b+ABZ*I_NAI_Iy5z_F3x_b;
  Double I_NAI_I6z_G3xz_b = I_NAI_K7z_F3x_b+ABZ*I_NAI_I6z_F3x_b;
  Double I_NAI_I5xz_G2x2y_b = I_NAI_K5xyz_F2xy_b+ABY*I_NAI_I5xz_F2xy_b;
  Double I_NAI_I4xyz_G2x2y_b = I_NAI_K4x2yz_F2xy_b+ABY*I_NAI_I4xyz_F2xy_b;
  Double I_NAI_I4x2z_G2x2y_b = I_NAI_K4xy2z_F2xy_b+ABY*I_NAI_I4x2z_F2xy_b;
  Double I_NAI_I3x2yz_G2x2y_b = I_NAI_K3x3yz_F2xy_b+ABY*I_NAI_I3x2yz_F2xy_b;
  Double I_NAI_I3xy2z_G2x2y_b = I_NAI_K3x2y2z_F2xy_b+ABY*I_NAI_I3xy2z_F2xy_b;
  Double I_NAI_I3x3z_G2x2y_b = I_NAI_K3xy3z_F2xy_b+ABY*I_NAI_I3x3z_F2xy_b;
  Double I_NAI_I2x3yz_G2x2y_b = I_NAI_K2x4yz_F2xy_b+ABY*I_NAI_I2x3yz_F2xy_b;
  Double I_NAI_I2x2y2z_G2x2y_b = I_NAI_K2x3y2z_F2xy_b+ABY*I_NAI_I2x2y2z_F2xy_b;
  Double I_NAI_I2xy3z_G2x2y_b = I_NAI_K2x2y3z_F2xy_b+ABY*I_NAI_I2xy3z_F2xy_b;
  Double I_NAI_I2x4z_G2x2y_b = I_NAI_K2xy4z_F2xy_b+ABY*I_NAI_I2x4z_F2xy_b;
  Double I_NAI_Ix4yz_G2x2y_b = I_NAI_Kx5yz_F2xy_b+ABY*I_NAI_Ix4yz_F2xy_b;
  Double I_NAI_Ix3y2z_G2x2y_b = I_NAI_Kx4y2z_F2xy_b+ABY*I_NAI_Ix3y2z_F2xy_b;
  Double I_NAI_Ix2y3z_G2x2y_b = I_NAI_Kx3y3z_F2xy_b+ABY*I_NAI_Ix2y3z_F2xy_b;
  Double I_NAI_Ixy4z_G2x2y_b = I_NAI_Kx2y4z_F2xy_b+ABY*I_NAI_Ixy4z_F2xy_b;
  Double I_NAI_Ix5z_G2x2y_b = I_NAI_Kxy5z_F2xy_b+ABY*I_NAI_Ix5z_F2xy_b;
  Double I_NAI_I5yz_G2x2y_b = I_NAI_K6yz_F2xy_b+ABY*I_NAI_I5yz_F2xy_b;
  Double I_NAI_I4y2z_G2x2y_b = I_NAI_K5y2z_F2xy_b+ABY*I_NAI_I4y2z_F2xy_b;
  Double I_NAI_I3y3z_G2x2y_b = I_NAI_K4y3z_F2xy_b+ABY*I_NAI_I3y3z_F2xy_b;
  Double I_NAI_I2y4z_G2x2y_b = I_NAI_K3y4z_F2xy_b+ABY*I_NAI_I2y4z_F2xy_b;
  Double I_NAI_Iy5z_G2x2y_b = I_NAI_K2y5z_F2xy_b+ABY*I_NAI_Iy5z_F2xy_b;
  Double I_NAI_I6z_G2x2y_b = I_NAI_Ky6z_F2xy_b+ABY*I_NAI_I6z_F2xy_b;
  Double I_NAI_I5xy_G2x2z_b = I_NAI_K5xyz_F2xz_b+ABZ*I_NAI_I5xy_F2xz_b;
  Double I_NAI_I4x2y_G2x2z_b = I_NAI_K4x2yz_F2xz_b+ABZ*I_NAI_I4x2y_F2xz_b;
  Double I_NAI_I4xyz_G2x2z_b = I_NAI_K4xy2z_F2xz_b+ABZ*I_NAI_I4xyz_F2xz_b;
  Double I_NAI_I3x3y_G2x2z_b = I_NAI_K3x3yz_F2xz_b+ABZ*I_NAI_I3x3y_F2xz_b;
  Double I_NAI_I3x2yz_G2x2z_b = I_NAI_K3x2y2z_F2xz_b+ABZ*I_NAI_I3x2yz_F2xz_b;
  Double I_NAI_I3xy2z_G2x2z_b = I_NAI_K3xy3z_F2xz_b+ABZ*I_NAI_I3xy2z_F2xz_b;
  Double I_NAI_I2x4y_G2x2z_b = I_NAI_K2x4yz_F2xz_b+ABZ*I_NAI_I2x4y_F2xz_b;
  Double I_NAI_I2x3yz_G2x2z_b = I_NAI_K2x3y2z_F2xz_b+ABZ*I_NAI_I2x3yz_F2xz_b;
  Double I_NAI_I2x2y2z_G2x2z_b = I_NAI_K2x2y3z_F2xz_b+ABZ*I_NAI_I2x2y2z_F2xz_b;
  Double I_NAI_I2xy3z_G2x2z_b = I_NAI_K2xy4z_F2xz_b+ABZ*I_NAI_I2xy3z_F2xz_b;
  Double I_NAI_Ix5y_G2x2z_b = I_NAI_Kx5yz_F2xz_b+ABZ*I_NAI_Ix5y_F2xz_b;
  Double I_NAI_Ix4yz_G2x2z_b = I_NAI_Kx4y2z_F2xz_b+ABZ*I_NAI_Ix4yz_F2xz_b;
  Double I_NAI_Ix3y2z_G2x2z_b = I_NAI_Kx3y3z_F2xz_b+ABZ*I_NAI_Ix3y2z_F2xz_b;
  Double I_NAI_Ix2y3z_G2x2z_b = I_NAI_Kx2y4z_F2xz_b+ABZ*I_NAI_Ix2y3z_F2xz_b;
  Double I_NAI_Ixy4z_G2x2z_b = I_NAI_Kxy5z_F2xz_b+ABZ*I_NAI_Ixy4z_F2xz_b;
  Double I_NAI_I6y_G2x2z_b = I_NAI_K6yz_F2xz_b+ABZ*I_NAI_I6y_F2xz_b;
  Double I_NAI_I5yz_G2x2z_b = I_NAI_K5y2z_F2xz_b+ABZ*I_NAI_I5yz_F2xz_b;
  Double I_NAI_I4y2z_G2x2z_b = I_NAI_K4y3z_F2xz_b+ABZ*I_NAI_I4y2z_F2xz_b;
  Double I_NAI_I3y3z_G2x2z_b = I_NAI_K3y4z_F2xz_b+ABZ*I_NAI_I3y3z_F2xz_b;
  Double I_NAI_I2y4z_G2x2z_b = I_NAI_K2y5z_F2xz_b+ABZ*I_NAI_I2y4z_F2xz_b;
  Double I_NAI_Iy5z_G2x2z_b = I_NAI_Ky6z_F2xz_b+ABZ*I_NAI_Iy5z_F2xz_b;
  Double I_NAI_I6x_Gx3y_b = I_NAI_K7x_F3y_b+ABX*I_NAI_I6x_F3y_b;
  Double I_NAI_I5xy_Gx3y_b = I_NAI_K6xy_F3y_b+ABX*I_NAI_I5xy_F3y_b;
  Double I_NAI_I5xz_Gx3y_b = I_NAI_K6xz_F3y_b+ABX*I_NAI_I5xz_F3y_b;
  Double I_NAI_I4x2y_Gx3y_b = I_NAI_K5x2y_F3y_b+ABX*I_NAI_I4x2y_F3y_b;
  Double I_NAI_I4xyz_Gx3y_b = I_NAI_K5xyz_F3y_b+ABX*I_NAI_I4xyz_F3y_b;
  Double I_NAI_I4x2z_Gx3y_b = I_NAI_K5x2z_F3y_b+ABX*I_NAI_I4x2z_F3y_b;
  Double I_NAI_I3x3y_Gx3y_b = I_NAI_K4x3y_F3y_b+ABX*I_NAI_I3x3y_F3y_b;
  Double I_NAI_I3x2yz_Gx3y_b = I_NAI_K4x2yz_F3y_b+ABX*I_NAI_I3x2yz_F3y_b;
  Double I_NAI_I3xy2z_Gx3y_b = I_NAI_K4xy2z_F3y_b+ABX*I_NAI_I3xy2z_F3y_b;
  Double I_NAI_I3x3z_Gx3y_b = I_NAI_K4x3z_F3y_b+ABX*I_NAI_I3x3z_F3y_b;
  Double I_NAI_I2x4y_Gx3y_b = I_NAI_K3x4y_F3y_b+ABX*I_NAI_I2x4y_F3y_b;
  Double I_NAI_I2x3yz_Gx3y_b = I_NAI_K3x3yz_F3y_b+ABX*I_NAI_I2x3yz_F3y_b;
  Double I_NAI_I2x2y2z_Gx3y_b = I_NAI_K3x2y2z_F3y_b+ABX*I_NAI_I2x2y2z_F3y_b;
  Double I_NAI_I2xy3z_Gx3y_b = I_NAI_K3xy3z_F3y_b+ABX*I_NAI_I2xy3z_F3y_b;
  Double I_NAI_I2x4z_Gx3y_b = I_NAI_K3x4z_F3y_b+ABX*I_NAI_I2x4z_F3y_b;
  Double I_NAI_Ix5y_Gx3y_b = I_NAI_K2x5y_F3y_b+ABX*I_NAI_Ix5y_F3y_b;
  Double I_NAI_Ix4yz_Gx3y_b = I_NAI_K2x4yz_F3y_b+ABX*I_NAI_Ix4yz_F3y_b;
  Double I_NAI_Ix3y2z_Gx3y_b = I_NAI_K2x3y2z_F3y_b+ABX*I_NAI_Ix3y2z_F3y_b;
  Double I_NAI_Ix2y3z_Gx3y_b = I_NAI_K2x2y3z_F3y_b+ABX*I_NAI_Ix2y3z_F3y_b;
  Double I_NAI_Ixy4z_Gx3y_b = I_NAI_K2xy4z_F3y_b+ABX*I_NAI_Ixy4z_F3y_b;
  Double I_NAI_Ix5z_Gx3y_b = I_NAI_K2x5z_F3y_b+ABX*I_NAI_Ix5z_F3y_b;
  Double I_NAI_I5yz_Gx3y_b = I_NAI_Kx5yz_F3y_b+ABX*I_NAI_I5yz_F3y_b;
  Double I_NAI_I4y2z_Gx3y_b = I_NAI_Kx4y2z_F3y_b+ABX*I_NAI_I4y2z_F3y_b;
  Double I_NAI_I3y3z_Gx3y_b = I_NAI_Kx3y3z_F3y_b+ABX*I_NAI_I3y3z_F3y_b;
  Double I_NAI_I2y4z_Gx3y_b = I_NAI_Kx2y4z_F3y_b+ABX*I_NAI_I2y4z_F3y_b;
  Double I_NAI_Iy5z_Gx3y_b = I_NAI_Kxy5z_F3y_b+ABX*I_NAI_Iy5z_F3y_b;
  Double I_NAI_I6z_Gx3y_b = I_NAI_Kx6z_F3y_b+ABX*I_NAI_I6z_F3y_b;
  Double I_NAI_I6x_Gx3z_b = I_NAI_K7x_F3z_b+ABX*I_NAI_I6x_F3z_b;
  Double I_NAI_I5xy_Gx3z_b = I_NAI_K6xy_F3z_b+ABX*I_NAI_I5xy_F3z_b;
  Double I_NAI_I5xz_Gx3z_b = I_NAI_K6xz_F3z_b+ABX*I_NAI_I5xz_F3z_b;
  Double I_NAI_I4x2y_Gx3z_b = I_NAI_K5x2y_F3z_b+ABX*I_NAI_I4x2y_F3z_b;
  Double I_NAI_I4xyz_Gx3z_b = I_NAI_K5xyz_F3z_b+ABX*I_NAI_I4xyz_F3z_b;
  Double I_NAI_I4x2z_Gx3z_b = I_NAI_K5x2z_F3z_b+ABX*I_NAI_I4x2z_F3z_b;
  Double I_NAI_I3x3y_Gx3z_b = I_NAI_K4x3y_F3z_b+ABX*I_NAI_I3x3y_F3z_b;
  Double I_NAI_I3x2yz_Gx3z_b = I_NAI_K4x2yz_F3z_b+ABX*I_NAI_I3x2yz_F3z_b;
  Double I_NAI_I3xy2z_Gx3z_b = I_NAI_K4xy2z_F3z_b+ABX*I_NAI_I3xy2z_F3z_b;
  Double I_NAI_I3x3z_Gx3z_b = I_NAI_K4x3z_F3z_b+ABX*I_NAI_I3x3z_F3z_b;
  Double I_NAI_I2x4y_Gx3z_b = I_NAI_K3x4y_F3z_b+ABX*I_NAI_I2x4y_F3z_b;
  Double I_NAI_I2x3yz_Gx3z_b = I_NAI_K3x3yz_F3z_b+ABX*I_NAI_I2x3yz_F3z_b;
  Double I_NAI_I2x2y2z_Gx3z_b = I_NAI_K3x2y2z_F3z_b+ABX*I_NAI_I2x2y2z_F3z_b;
  Double I_NAI_I2xy3z_Gx3z_b = I_NAI_K3xy3z_F3z_b+ABX*I_NAI_I2xy3z_F3z_b;
  Double I_NAI_I2x4z_Gx3z_b = I_NAI_K3x4z_F3z_b+ABX*I_NAI_I2x4z_F3z_b;
  Double I_NAI_Ix5y_Gx3z_b = I_NAI_K2x5y_F3z_b+ABX*I_NAI_Ix5y_F3z_b;
  Double I_NAI_Ix4yz_Gx3z_b = I_NAI_K2x4yz_F3z_b+ABX*I_NAI_Ix4yz_F3z_b;
  Double I_NAI_Ix3y2z_Gx3z_b = I_NAI_K2x3y2z_F3z_b+ABX*I_NAI_Ix3y2z_F3z_b;
  Double I_NAI_Ix2y3z_Gx3z_b = I_NAI_K2x2y3z_F3z_b+ABX*I_NAI_Ix2y3z_F3z_b;
  Double I_NAI_Ixy4z_Gx3z_b = I_NAI_K2xy4z_F3z_b+ABX*I_NAI_Ixy4z_F3z_b;
  Double I_NAI_Ix5z_Gx3z_b = I_NAI_K2x5z_F3z_b+ABX*I_NAI_Ix5z_F3z_b;
  Double I_NAI_I6y_Gx3z_b = I_NAI_Kx6y_F3z_b+ABX*I_NAI_I6y_F3z_b;
  Double I_NAI_I5yz_Gx3z_b = I_NAI_Kx5yz_F3z_b+ABX*I_NAI_I5yz_F3z_b;
  Double I_NAI_I4y2z_Gx3z_b = I_NAI_Kx4y2z_F3z_b+ABX*I_NAI_I4y2z_F3z_b;
  Double I_NAI_I3y3z_Gx3z_b = I_NAI_Kx3y3z_F3z_b+ABX*I_NAI_I3y3z_F3z_b;
  Double I_NAI_I2y4z_Gx3z_b = I_NAI_Kx2y4z_F3z_b+ABX*I_NAI_I2y4z_F3z_b;
  Double I_NAI_Iy5z_Gx3z_b = I_NAI_Kxy5z_F3z_b+ABX*I_NAI_Iy5z_F3z_b;
  Double I_NAI_I6x_G4y_b = I_NAI_K6xy_F3y_b+ABY*I_NAI_I6x_F3y_b;
  Double I_NAI_I5xy_G4y_b = I_NAI_K5x2y_F3y_b+ABY*I_NAI_I5xy_F3y_b;
  Double I_NAI_I5xz_G4y_b = I_NAI_K5xyz_F3y_b+ABY*I_NAI_I5xz_F3y_b;
  Double I_NAI_I4x2y_G4y_b = I_NAI_K4x3y_F3y_b+ABY*I_NAI_I4x2y_F3y_b;
  Double I_NAI_I4xyz_G4y_b = I_NAI_K4x2yz_F3y_b+ABY*I_NAI_I4xyz_F3y_b;
  Double I_NAI_I4x2z_G4y_b = I_NAI_K4xy2z_F3y_b+ABY*I_NAI_I4x2z_F3y_b;
  Double I_NAI_I3x3y_G4y_b = I_NAI_K3x4y_F3y_b+ABY*I_NAI_I3x3y_F3y_b;
  Double I_NAI_I3x2yz_G4y_b = I_NAI_K3x3yz_F3y_b+ABY*I_NAI_I3x2yz_F3y_b;
  Double I_NAI_I3xy2z_G4y_b = I_NAI_K3x2y2z_F3y_b+ABY*I_NAI_I3xy2z_F3y_b;
  Double I_NAI_I3x3z_G4y_b = I_NAI_K3xy3z_F3y_b+ABY*I_NAI_I3x3z_F3y_b;
  Double I_NAI_I2x4y_G4y_b = I_NAI_K2x5y_F3y_b+ABY*I_NAI_I2x4y_F3y_b;
  Double I_NAI_I2x3yz_G4y_b = I_NAI_K2x4yz_F3y_b+ABY*I_NAI_I2x3yz_F3y_b;
  Double I_NAI_I2x2y2z_G4y_b = I_NAI_K2x3y2z_F3y_b+ABY*I_NAI_I2x2y2z_F3y_b;
  Double I_NAI_I2xy3z_G4y_b = I_NAI_K2x2y3z_F3y_b+ABY*I_NAI_I2xy3z_F3y_b;
  Double I_NAI_I2x4z_G4y_b = I_NAI_K2xy4z_F3y_b+ABY*I_NAI_I2x4z_F3y_b;
  Double I_NAI_Ix5y_G4y_b = I_NAI_Kx6y_F3y_b+ABY*I_NAI_Ix5y_F3y_b;
  Double I_NAI_Ix4yz_G4y_b = I_NAI_Kx5yz_F3y_b+ABY*I_NAI_Ix4yz_F3y_b;
  Double I_NAI_Ix3y2z_G4y_b = I_NAI_Kx4y2z_F3y_b+ABY*I_NAI_Ix3y2z_F3y_b;
  Double I_NAI_Ix2y3z_G4y_b = I_NAI_Kx3y3z_F3y_b+ABY*I_NAI_Ix2y3z_F3y_b;
  Double I_NAI_Ixy4z_G4y_b = I_NAI_Kx2y4z_F3y_b+ABY*I_NAI_Ixy4z_F3y_b;
  Double I_NAI_Ix5z_G4y_b = I_NAI_Kxy5z_F3y_b+ABY*I_NAI_Ix5z_F3y_b;
  Double I_NAI_I6y_G4y_b = I_NAI_K7y_F3y_b+ABY*I_NAI_I6y_F3y_b;
  Double I_NAI_I5yz_G4y_b = I_NAI_K6yz_F3y_b+ABY*I_NAI_I5yz_F3y_b;
  Double I_NAI_I4y2z_G4y_b = I_NAI_K5y2z_F3y_b+ABY*I_NAI_I4y2z_F3y_b;
  Double I_NAI_I3y3z_G4y_b = I_NAI_K4y3z_F3y_b+ABY*I_NAI_I3y3z_F3y_b;
  Double I_NAI_I2y4z_G4y_b = I_NAI_K3y4z_F3y_b+ABY*I_NAI_I2y4z_F3y_b;
  Double I_NAI_Iy5z_G4y_b = I_NAI_K2y5z_F3y_b+ABY*I_NAI_Iy5z_F3y_b;
  Double I_NAI_I6z_G4y_b = I_NAI_Ky6z_F3y_b+ABY*I_NAI_I6z_F3y_b;
  Double I_NAI_I5xz_G3yz_b = I_NAI_K5x2z_F3y_b+ABZ*I_NAI_I5xz_F3y_b;
  Double I_NAI_I4xyz_G3yz_b = I_NAI_K4xy2z_F3y_b+ABZ*I_NAI_I4xyz_F3y_b;
  Double I_NAI_I4x2z_G3yz_b = I_NAI_K4x3z_F3y_b+ABZ*I_NAI_I4x2z_F3y_b;
  Double I_NAI_I3x2yz_G3yz_b = I_NAI_K3x2y2z_F3y_b+ABZ*I_NAI_I3x2yz_F3y_b;
  Double I_NAI_I3xy2z_G3yz_b = I_NAI_K3xy3z_F3y_b+ABZ*I_NAI_I3xy2z_F3y_b;
  Double I_NAI_I3x3z_G3yz_b = I_NAI_K3x4z_F3y_b+ABZ*I_NAI_I3x3z_F3y_b;
  Double I_NAI_I2x3yz_G3yz_b = I_NAI_K2x3y2z_F3y_b+ABZ*I_NAI_I2x3yz_F3y_b;
  Double I_NAI_I2x2y2z_G3yz_b = I_NAI_K2x2y3z_F3y_b+ABZ*I_NAI_I2x2y2z_F3y_b;
  Double I_NAI_I2xy3z_G3yz_b = I_NAI_K2xy4z_F3y_b+ABZ*I_NAI_I2xy3z_F3y_b;
  Double I_NAI_I2x4z_G3yz_b = I_NAI_K2x5z_F3y_b+ABZ*I_NAI_I2x4z_F3y_b;
  Double I_NAI_Ix4yz_G3yz_b = I_NAI_Kx4y2z_F3y_b+ABZ*I_NAI_Ix4yz_F3y_b;
  Double I_NAI_Ix3y2z_G3yz_b = I_NAI_Kx3y3z_F3y_b+ABZ*I_NAI_Ix3y2z_F3y_b;
  Double I_NAI_Ix2y3z_G3yz_b = I_NAI_Kx2y4z_F3y_b+ABZ*I_NAI_Ix2y3z_F3y_b;
  Double I_NAI_Ixy4z_G3yz_b = I_NAI_Kxy5z_F3y_b+ABZ*I_NAI_Ixy4z_F3y_b;
  Double I_NAI_Ix5z_G3yz_b = I_NAI_Kx6z_F3y_b+ABZ*I_NAI_Ix5z_F3y_b;
  Double I_NAI_I5yz_G3yz_b = I_NAI_K5y2z_F3y_b+ABZ*I_NAI_I5yz_F3y_b;
  Double I_NAI_I4y2z_G3yz_b = I_NAI_K4y3z_F3y_b+ABZ*I_NAI_I4y2z_F3y_b;
  Double I_NAI_I3y3z_G3yz_b = I_NAI_K3y4z_F3y_b+ABZ*I_NAI_I3y3z_F3y_b;
  Double I_NAI_I2y4z_G3yz_b = I_NAI_K2y5z_F3y_b+ABZ*I_NAI_I2y4z_F3y_b;
  Double I_NAI_Iy5z_G3yz_b = I_NAI_Ky6z_F3y_b+ABZ*I_NAI_Iy5z_F3y_b;
  Double I_NAI_I6z_G3yz_b = I_NAI_K7z_F3y_b+ABZ*I_NAI_I6z_F3y_b;
  Double I_NAI_I6x_G2y2z_b = I_NAI_K6xz_F2yz_b+ABZ*I_NAI_I6x_F2yz_b;
  Double I_NAI_I5xy_G2y2z_b = I_NAI_K5xyz_F2yz_b+ABZ*I_NAI_I5xy_F2yz_b;
  Double I_NAI_I5xz_G2y2z_b = I_NAI_K5x2z_F2yz_b+ABZ*I_NAI_I5xz_F2yz_b;
  Double I_NAI_I4x2y_G2y2z_b = I_NAI_K4x2yz_F2yz_b+ABZ*I_NAI_I4x2y_F2yz_b;
  Double I_NAI_I4xyz_G2y2z_b = I_NAI_K4xy2z_F2yz_b+ABZ*I_NAI_I4xyz_F2yz_b;
  Double I_NAI_I4x2z_G2y2z_b = I_NAI_K4x3z_F2yz_b+ABZ*I_NAI_I4x2z_F2yz_b;
  Double I_NAI_I3x3y_G2y2z_b = I_NAI_K3x3yz_F2yz_b+ABZ*I_NAI_I3x3y_F2yz_b;
  Double I_NAI_I3x2yz_G2y2z_b = I_NAI_K3x2y2z_F2yz_b+ABZ*I_NAI_I3x2yz_F2yz_b;
  Double I_NAI_I3xy2z_G2y2z_b = I_NAI_K3xy3z_F2yz_b+ABZ*I_NAI_I3xy2z_F2yz_b;
  Double I_NAI_I3x3z_G2y2z_b = I_NAI_K3x4z_F2yz_b+ABZ*I_NAI_I3x3z_F2yz_b;
  Double I_NAI_I2x4y_G2y2z_b = I_NAI_K2x4yz_F2yz_b+ABZ*I_NAI_I2x4y_F2yz_b;
  Double I_NAI_I2x3yz_G2y2z_b = I_NAI_K2x3y2z_F2yz_b+ABZ*I_NAI_I2x3yz_F2yz_b;
  Double I_NAI_I2x2y2z_G2y2z_b = I_NAI_K2x2y3z_F2yz_b+ABZ*I_NAI_I2x2y2z_F2yz_b;
  Double I_NAI_I2xy3z_G2y2z_b = I_NAI_K2xy4z_F2yz_b+ABZ*I_NAI_I2xy3z_F2yz_b;
  Double I_NAI_I2x4z_G2y2z_b = I_NAI_K2x5z_F2yz_b+ABZ*I_NAI_I2x4z_F2yz_b;
  Double I_NAI_Ix5y_G2y2z_b = I_NAI_Kx5yz_F2yz_b+ABZ*I_NAI_Ix5y_F2yz_b;
  Double I_NAI_Ix4yz_G2y2z_b = I_NAI_Kx4y2z_F2yz_b+ABZ*I_NAI_Ix4yz_F2yz_b;
  Double I_NAI_Ix3y2z_G2y2z_b = I_NAI_Kx3y3z_F2yz_b+ABZ*I_NAI_Ix3y2z_F2yz_b;
  Double I_NAI_Ix2y3z_G2y2z_b = I_NAI_Kx2y4z_F2yz_b+ABZ*I_NAI_Ix2y3z_F2yz_b;
  Double I_NAI_Ixy4z_G2y2z_b = I_NAI_Kxy5z_F2yz_b+ABZ*I_NAI_Ixy4z_F2yz_b;
  Double I_NAI_Ix5z_G2y2z_b = I_NAI_Kx6z_F2yz_b+ABZ*I_NAI_Ix5z_F2yz_b;
  Double I_NAI_I5xy_Gy3z_b = I_NAI_K5x2y_F3z_b+ABY*I_NAI_I5xy_F3z_b;
  Double I_NAI_I4x2y_Gy3z_b = I_NAI_K4x3y_F3z_b+ABY*I_NAI_I4x2y_F3z_b;
  Double I_NAI_I4xyz_Gy3z_b = I_NAI_K4x2yz_F3z_b+ABY*I_NAI_I4xyz_F3z_b;
  Double I_NAI_I3x3y_Gy3z_b = I_NAI_K3x4y_F3z_b+ABY*I_NAI_I3x3y_F3z_b;
  Double I_NAI_I3x2yz_Gy3z_b = I_NAI_K3x3yz_F3z_b+ABY*I_NAI_I3x2yz_F3z_b;
  Double I_NAI_I3xy2z_Gy3z_b = I_NAI_K3x2y2z_F3z_b+ABY*I_NAI_I3xy2z_F3z_b;
  Double I_NAI_I2x4y_Gy3z_b = I_NAI_K2x5y_F3z_b+ABY*I_NAI_I2x4y_F3z_b;
  Double I_NAI_I2x3yz_Gy3z_b = I_NAI_K2x4yz_F3z_b+ABY*I_NAI_I2x3yz_F3z_b;
  Double I_NAI_I2x2y2z_Gy3z_b = I_NAI_K2x3y2z_F3z_b+ABY*I_NAI_I2x2y2z_F3z_b;
  Double I_NAI_I2xy3z_Gy3z_b = I_NAI_K2x2y3z_F3z_b+ABY*I_NAI_I2xy3z_F3z_b;
  Double I_NAI_Ix5y_Gy3z_b = I_NAI_Kx6y_F3z_b+ABY*I_NAI_Ix5y_F3z_b;
  Double I_NAI_Ix4yz_Gy3z_b = I_NAI_Kx5yz_F3z_b+ABY*I_NAI_Ix4yz_F3z_b;
  Double I_NAI_Ix3y2z_Gy3z_b = I_NAI_Kx4y2z_F3z_b+ABY*I_NAI_Ix3y2z_F3z_b;
  Double I_NAI_Ix2y3z_Gy3z_b = I_NAI_Kx3y3z_F3z_b+ABY*I_NAI_Ix2y3z_F3z_b;
  Double I_NAI_Ixy4z_Gy3z_b = I_NAI_Kx2y4z_F3z_b+ABY*I_NAI_Ixy4z_F3z_b;
  Double I_NAI_I6y_Gy3z_b = I_NAI_K7y_F3z_b+ABY*I_NAI_I6y_F3z_b;
  Double I_NAI_I5yz_Gy3z_b = I_NAI_K6yz_F3z_b+ABY*I_NAI_I5yz_F3z_b;
  Double I_NAI_I4y2z_Gy3z_b = I_NAI_K5y2z_F3z_b+ABY*I_NAI_I4y2z_F3z_b;
  Double I_NAI_I3y3z_Gy3z_b = I_NAI_K4y3z_F3z_b+ABY*I_NAI_I3y3z_F3z_b;
  Double I_NAI_I2y4z_Gy3z_b = I_NAI_K3y4z_F3z_b+ABY*I_NAI_I2y4z_F3z_b;
  Double I_NAI_Iy5z_Gy3z_b = I_NAI_K2y5z_F3z_b+ABY*I_NAI_Iy5z_F3z_b;
  Double I_NAI_I6x_G4z_b = I_NAI_K6xz_F3z_b+ABZ*I_NAI_I6x_F3z_b;
  Double I_NAI_I5xy_G4z_b = I_NAI_K5xyz_F3z_b+ABZ*I_NAI_I5xy_F3z_b;
  Double I_NAI_I5xz_G4z_b = I_NAI_K5x2z_F3z_b+ABZ*I_NAI_I5xz_F3z_b;
  Double I_NAI_I4x2y_G4z_b = I_NAI_K4x2yz_F3z_b+ABZ*I_NAI_I4x2y_F3z_b;
  Double I_NAI_I4xyz_G4z_b = I_NAI_K4xy2z_F3z_b+ABZ*I_NAI_I4xyz_F3z_b;
  Double I_NAI_I4x2z_G4z_b = I_NAI_K4x3z_F3z_b+ABZ*I_NAI_I4x2z_F3z_b;
  Double I_NAI_I3x3y_G4z_b = I_NAI_K3x3yz_F3z_b+ABZ*I_NAI_I3x3y_F3z_b;
  Double I_NAI_I3x2yz_G4z_b = I_NAI_K3x2y2z_F3z_b+ABZ*I_NAI_I3x2yz_F3z_b;
  Double I_NAI_I3xy2z_G4z_b = I_NAI_K3xy3z_F3z_b+ABZ*I_NAI_I3xy2z_F3z_b;
  Double I_NAI_I3x3z_G4z_b = I_NAI_K3x4z_F3z_b+ABZ*I_NAI_I3x3z_F3z_b;
  Double I_NAI_I2x4y_G4z_b = I_NAI_K2x4yz_F3z_b+ABZ*I_NAI_I2x4y_F3z_b;
  Double I_NAI_I2x3yz_G4z_b = I_NAI_K2x3y2z_F3z_b+ABZ*I_NAI_I2x3yz_F3z_b;
  Double I_NAI_I2x2y2z_G4z_b = I_NAI_K2x2y3z_F3z_b+ABZ*I_NAI_I2x2y2z_F3z_b;
  Double I_NAI_I2xy3z_G4z_b = I_NAI_K2xy4z_F3z_b+ABZ*I_NAI_I2xy3z_F3z_b;
  Double I_NAI_I2x4z_G4z_b = I_NAI_K2x5z_F3z_b+ABZ*I_NAI_I2x4z_F3z_b;
  Double I_NAI_Ix5y_G4z_b = I_NAI_Kx5yz_F3z_b+ABZ*I_NAI_Ix5y_F3z_b;
  Double I_NAI_Ix4yz_G4z_b = I_NAI_Kx4y2z_F3z_b+ABZ*I_NAI_Ix4yz_F3z_b;
  Double I_NAI_Ix3y2z_G4z_b = I_NAI_Kx3y3z_F3z_b+ABZ*I_NAI_Ix3y2z_F3z_b;
  Double I_NAI_Ix2y3z_G4z_b = I_NAI_Kx2y4z_F3z_b+ABZ*I_NAI_Ix2y3z_F3z_b;
  Double I_NAI_Ixy4z_G4z_b = I_NAI_Kxy5z_F3z_b+ABZ*I_NAI_Ixy4z_F3z_b;
  Double I_NAI_Ix5z_G4z_b = I_NAI_Kx6z_F3z_b+ABZ*I_NAI_Ix5z_F3z_b;
  Double I_NAI_I6y_G4z_b = I_NAI_K6yz_F3z_b+ABZ*I_NAI_I6y_F3z_b;
  Double I_NAI_I5yz_G4z_b = I_NAI_K5y2z_F3z_b+ABZ*I_NAI_I5yz_F3z_b;
  Double I_NAI_I4y2z_G4z_b = I_NAI_K4y3z_F3z_b+ABZ*I_NAI_I4y2z_F3z_b;
  Double I_NAI_I3y3z_G4z_b = I_NAI_K3y4z_F3z_b+ABZ*I_NAI_I3y3z_F3z_b;
  Double I_NAI_I2y4z_G4z_b = I_NAI_K2y5z_F3z_b+ABZ*I_NAI_I2y4z_F3z_b;
  Double I_NAI_Iy5z_G4z_b = I_NAI_Ky6z_F3z_b+ABZ*I_NAI_Iy5z_F3z_b;
  Double I_NAI_I6z_G4z_b = I_NAI_K7z_F3z_b+ABZ*I_NAI_I6z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_H_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_G_b
   * RHS shell quartet name: SQ_NAI_H_G_b
   ************************************************************/
  Double I_NAI_H5x_H5x_b = I_NAI_I6x_G4x_b+ABX*I_NAI_H5x_G4x_b;
  Double I_NAI_H4xy_H5x_b = I_NAI_I5xy_G4x_b+ABX*I_NAI_H4xy_G4x_b;
  Double I_NAI_H4xz_H5x_b = I_NAI_I5xz_G4x_b+ABX*I_NAI_H4xz_G4x_b;
  Double I_NAI_H3x2y_H5x_b = I_NAI_I4x2y_G4x_b+ABX*I_NAI_H3x2y_G4x_b;
  Double I_NAI_H3xyz_H5x_b = I_NAI_I4xyz_G4x_b+ABX*I_NAI_H3xyz_G4x_b;
  Double I_NAI_H3x2z_H5x_b = I_NAI_I4x2z_G4x_b+ABX*I_NAI_H3x2z_G4x_b;
  Double I_NAI_H2x3y_H5x_b = I_NAI_I3x3y_G4x_b+ABX*I_NAI_H2x3y_G4x_b;
  Double I_NAI_H2x2yz_H5x_b = I_NAI_I3x2yz_G4x_b+ABX*I_NAI_H2x2yz_G4x_b;
  Double I_NAI_H2xy2z_H5x_b = I_NAI_I3xy2z_G4x_b+ABX*I_NAI_H2xy2z_G4x_b;
  Double I_NAI_H2x3z_H5x_b = I_NAI_I3x3z_G4x_b+ABX*I_NAI_H2x3z_G4x_b;
  Double I_NAI_Hx4y_H5x_b = I_NAI_I2x4y_G4x_b+ABX*I_NAI_Hx4y_G4x_b;
  Double I_NAI_Hx3yz_H5x_b = I_NAI_I2x3yz_G4x_b+ABX*I_NAI_Hx3yz_G4x_b;
  Double I_NAI_Hx2y2z_H5x_b = I_NAI_I2x2y2z_G4x_b+ABX*I_NAI_Hx2y2z_G4x_b;
  Double I_NAI_Hxy3z_H5x_b = I_NAI_I2xy3z_G4x_b+ABX*I_NAI_Hxy3z_G4x_b;
  Double I_NAI_Hx4z_H5x_b = I_NAI_I2x4z_G4x_b+ABX*I_NAI_Hx4z_G4x_b;
  Double I_NAI_H5y_H5x_b = I_NAI_Ix5y_G4x_b+ABX*I_NAI_H5y_G4x_b;
  Double I_NAI_H4yz_H5x_b = I_NAI_Ix4yz_G4x_b+ABX*I_NAI_H4yz_G4x_b;
  Double I_NAI_H3y2z_H5x_b = I_NAI_Ix3y2z_G4x_b+ABX*I_NAI_H3y2z_G4x_b;
  Double I_NAI_H2y3z_H5x_b = I_NAI_Ix2y3z_G4x_b+ABX*I_NAI_H2y3z_G4x_b;
  Double I_NAI_Hy4z_H5x_b = I_NAI_Ixy4z_G4x_b+ABX*I_NAI_Hy4z_G4x_b;
  Double I_NAI_H5z_H5x_b = I_NAI_Ix5z_G4x_b+ABX*I_NAI_H5z_G4x_b;
  Double I_NAI_H5x_H4xy_b = I_NAI_I5xy_G4x_b+ABY*I_NAI_H5x_G4x_b;
  Double I_NAI_H4xy_H4xy_b = I_NAI_I4x2y_G4x_b+ABY*I_NAI_H4xy_G4x_b;
  Double I_NAI_H4xz_H4xy_b = I_NAI_I4xyz_G4x_b+ABY*I_NAI_H4xz_G4x_b;
  Double I_NAI_H3x2y_H4xy_b = I_NAI_I3x3y_G4x_b+ABY*I_NAI_H3x2y_G4x_b;
  Double I_NAI_H3xyz_H4xy_b = I_NAI_I3x2yz_G4x_b+ABY*I_NAI_H3xyz_G4x_b;
  Double I_NAI_H3x2z_H4xy_b = I_NAI_I3xy2z_G4x_b+ABY*I_NAI_H3x2z_G4x_b;
  Double I_NAI_H2x3y_H4xy_b = I_NAI_I2x4y_G4x_b+ABY*I_NAI_H2x3y_G4x_b;
  Double I_NAI_H2x2yz_H4xy_b = I_NAI_I2x3yz_G4x_b+ABY*I_NAI_H2x2yz_G4x_b;
  Double I_NAI_H2xy2z_H4xy_b = I_NAI_I2x2y2z_G4x_b+ABY*I_NAI_H2xy2z_G4x_b;
  Double I_NAI_H2x3z_H4xy_b = I_NAI_I2xy3z_G4x_b+ABY*I_NAI_H2x3z_G4x_b;
  Double I_NAI_Hx4y_H4xy_b = I_NAI_Ix5y_G4x_b+ABY*I_NAI_Hx4y_G4x_b;
  Double I_NAI_Hx3yz_H4xy_b = I_NAI_Ix4yz_G4x_b+ABY*I_NAI_Hx3yz_G4x_b;
  Double I_NAI_Hx2y2z_H4xy_b = I_NAI_Ix3y2z_G4x_b+ABY*I_NAI_Hx2y2z_G4x_b;
  Double I_NAI_Hxy3z_H4xy_b = I_NAI_Ix2y3z_G4x_b+ABY*I_NAI_Hxy3z_G4x_b;
  Double I_NAI_Hx4z_H4xy_b = I_NAI_Ixy4z_G4x_b+ABY*I_NAI_Hx4z_G4x_b;
  Double I_NAI_H5y_H4xy_b = I_NAI_I6y_G4x_b+ABY*I_NAI_H5y_G4x_b;
  Double I_NAI_H4yz_H4xy_b = I_NAI_I5yz_G4x_b+ABY*I_NAI_H4yz_G4x_b;
  Double I_NAI_H3y2z_H4xy_b = I_NAI_I4y2z_G4x_b+ABY*I_NAI_H3y2z_G4x_b;
  Double I_NAI_H2y3z_H4xy_b = I_NAI_I3y3z_G4x_b+ABY*I_NAI_H2y3z_G4x_b;
  Double I_NAI_Hy4z_H4xy_b = I_NAI_I2y4z_G4x_b+ABY*I_NAI_Hy4z_G4x_b;
  Double I_NAI_H5z_H4xy_b = I_NAI_Iy5z_G4x_b+ABY*I_NAI_H5z_G4x_b;
  Double I_NAI_H5x_H4xz_b = I_NAI_I5xz_G4x_b+ABZ*I_NAI_H5x_G4x_b;
  Double I_NAI_H4xy_H4xz_b = I_NAI_I4xyz_G4x_b+ABZ*I_NAI_H4xy_G4x_b;
  Double I_NAI_H4xz_H4xz_b = I_NAI_I4x2z_G4x_b+ABZ*I_NAI_H4xz_G4x_b;
  Double I_NAI_H3x2y_H4xz_b = I_NAI_I3x2yz_G4x_b+ABZ*I_NAI_H3x2y_G4x_b;
  Double I_NAI_H3xyz_H4xz_b = I_NAI_I3xy2z_G4x_b+ABZ*I_NAI_H3xyz_G4x_b;
  Double I_NAI_H3x2z_H4xz_b = I_NAI_I3x3z_G4x_b+ABZ*I_NAI_H3x2z_G4x_b;
  Double I_NAI_H2x3y_H4xz_b = I_NAI_I2x3yz_G4x_b+ABZ*I_NAI_H2x3y_G4x_b;
  Double I_NAI_H2x2yz_H4xz_b = I_NAI_I2x2y2z_G4x_b+ABZ*I_NAI_H2x2yz_G4x_b;
  Double I_NAI_H2xy2z_H4xz_b = I_NAI_I2xy3z_G4x_b+ABZ*I_NAI_H2xy2z_G4x_b;
  Double I_NAI_H2x3z_H4xz_b = I_NAI_I2x4z_G4x_b+ABZ*I_NAI_H2x3z_G4x_b;
  Double I_NAI_Hx4y_H4xz_b = I_NAI_Ix4yz_G4x_b+ABZ*I_NAI_Hx4y_G4x_b;
  Double I_NAI_Hx3yz_H4xz_b = I_NAI_Ix3y2z_G4x_b+ABZ*I_NAI_Hx3yz_G4x_b;
  Double I_NAI_Hx2y2z_H4xz_b = I_NAI_Ix2y3z_G4x_b+ABZ*I_NAI_Hx2y2z_G4x_b;
  Double I_NAI_Hxy3z_H4xz_b = I_NAI_Ixy4z_G4x_b+ABZ*I_NAI_Hxy3z_G4x_b;
  Double I_NAI_Hx4z_H4xz_b = I_NAI_Ix5z_G4x_b+ABZ*I_NAI_Hx4z_G4x_b;
  Double I_NAI_H5y_H4xz_b = I_NAI_I5yz_G4x_b+ABZ*I_NAI_H5y_G4x_b;
  Double I_NAI_H4yz_H4xz_b = I_NAI_I4y2z_G4x_b+ABZ*I_NAI_H4yz_G4x_b;
  Double I_NAI_H3y2z_H4xz_b = I_NAI_I3y3z_G4x_b+ABZ*I_NAI_H3y2z_G4x_b;
  Double I_NAI_H2y3z_H4xz_b = I_NAI_I2y4z_G4x_b+ABZ*I_NAI_H2y3z_G4x_b;
  Double I_NAI_Hy4z_H4xz_b = I_NAI_Iy5z_G4x_b+ABZ*I_NAI_Hy4z_G4x_b;
  Double I_NAI_H5z_H4xz_b = I_NAI_I6z_G4x_b+ABZ*I_NAI_H5z_G4x_b;
  Double I_NAI_H5x_H3x2y_b = I_NAI_I5xy_G3xy_b+ABY*I_NAI_H5x_G3xy_b;
  Double I_NAI_H4xy_H3x2y_b = I_NAI_I4x2y_G3xy_b+ABY*I_NAI_H4xy_G3xy_b;
  Double I_NAI_H4xz_H3x2y_b = I_NAI_I4xyz_G3xy_b+ABY*I_NAI_H4xz_G3xy_b;
  Double I_NAI_H3x2y_H3x2y_b = I_NAI_I3x3y_G3xy_b+ABY*I_NAI_H3x2y_G3xy_b;
  Double I_NAI_H3xyz_H3x2y_b = I_NAI_I3x2yz_G3xy_b+ABY*I_NAI_H3xyz_G3xy_b;
  Double I_NAI_H3x2z_H3x2y_b = I_NAI_I3xy2z_G3xy_b+ABY*I_NAI_H3x2z_G3xy_b;
  Double I_NAI_H2x3y_H3x2y_b = I_NAI_I2x4y_G3xy_b+ABY*I_NAI_H2x3y_G3xy_b;
  Double I_NAI_H2x2yz_H3x2y_b = I_NAI_I2x3yz_G3xy_b+ABY*I_NAI_H2x2yz_G3xy_b;
  Double I_NAI_H2xy2z_H3x2y_b = I_NAI_I2x2y2z_G3xy_b+ABY*I_NAI_H2xy2z_G3xy_b;
  Double I_NAI_H2x3z_H3x2y_b = I_NAI_I2xy3z_G3xy_b+ABY*I_NAI_H2x3z_G3xy_b;
  Double I_NAI_Hx4y_H3x2y_b = I_NAI_Ix5y_G3xy_b+ABY*I_NAI_Hx4y_G3xy_b;
  Double I_NAI_Hx3yz_H3x2y_b = I_NAI_Ix4yz_G3xy_b+ABY*I_NAI_Hx3yz_G3xy_b;
  Double I_NAI_Hx2y2z_H3x2y_b = I_NAI_Ix3y2z_G3xy_b+ABY*I_NAI_Hx2y2z_G3xy_b;
  Double I_NAI_Hxy3z_H3x2y_b = I_NAI_Ix2y3z_G3xy_b+ABY*I_NAI_Hxy3z_G3xy_b;
  Double I_NAI_Hx4z_H3x2y_b = I_NAI_Ixy4z_G3xy_b+ABY*I_NAI_Hx4z_G3xy_b;
  Double I_NAI_H5y_H3x2y_b = I_NAI_I6y_G3xy_b+ABY*I_NAI_H5y_G3xy_b;
  Double I_NAI_H4yz_H3x2y_b = I_NAI_I5yz_G3xy_b+ABY*I_NAI_H4yz_G3xy_b;
  Double I_NAI_H3y2z_H3x2y_b = I_NAI_I4y2z_G3xy_b+ABY*I_NAI_H3y2z_G3xy_b;
  Double I_NAI_H2y3z_H3x2y_b = I_NAI_I3y3z_G3xy_b+ABY*I_NAI_H2y3z_G3xy_b;
  Double I_NAI_Hy4z_H3x2y_b = I_NAI_I2y4z_G3xy_b+ABY*I_NAI_Hy4z_G3xy_b;
  Double I_NAI_H5z_H3x2y_b = I_NAI_Iy5z_G3xy_b+ABY*I_NAI_H5z_G3xy_b;
  Double I_NAI_H5x_H3xyz_b = I_NAI_I5xz_G3xy_b+ABZ*I_NAI_H5x_G3xy_b;
  Double I_NAI_H4xy_H3xyz_b = I_NAI_I4xyz_G3xy_b+ABZ*I_NAI_H4xy_G3xy_b;
  Double I_NAI_H4xz_H3xyz_b = I_NAI_I4x2z_G3xy_b+ABZ*I_NAI_H4xz_G3xy_b;
  Double I_NAI_H3x2y_H3xyz_b = I_NAI_I3x2yz_G3xy_b+ABZ*I_NAI_H3x2y_G3xy_b;
  Double I_NAI_H3xyz_H3xyz_b = I_NAI_I3xy2z_G3xy_b+ABZ*I_NAI_H3xyz_G3xy_b;
  Double I_NAI_H3x2z_H3xyz_b = I_NAI_I3x3z_G3xy_b+ABZ*I_NAI_H3x2z_G3xy_b;
  Double I_NAI_H2x3y_H3xyz_b = I_NAI_I2x3yz_G3xy_b+ABZ*I_NAI_H2x3y_G3xy_b;
  Double I_NAI_H2x2yz_H3xyz_b = I_NAI_I2x2y2z_G3xy_b+ABZ*I_NAI_H2x2yz_G3xy_b;
  Double I_NAI_H2xy2z_H3xyz_b = I_NAI_I2xy3z_G3xy_b+ABZ*I_NAI_H2xy2z_G3xy_b;
  Double I_NAI_H2x3z_H3xyz_b = I_NAI_I2x4z_G3xy_b+ABZ*I_NAI_H2x3z_G3xy_b;
  Double I_NAI_Hx4y_H3xyz_b = I_NAI_Ix4yz_G3xy_b+ABZ*I_NAI_Hx4y_G3xy_b;
  Double I_NAI_Hx3yz_H3xyz_b = I_NAI_Ix3y2z_G3xy_b+ABZ*I_NAI_Hx3yz_G3xy_b;
  Double I_NAI_Hx2y2z_H3xyz_b = I_NAI_Ix2y3z_G3xy_b+ABZ*I_NAI_Hx2y2z_G3xy_b;
  Double I_NAI_Hxy3z_H3xyz_b = I_NAI_Ixy4z_G3xy_b+ABZ*I_NAI_Hxy3z_G3xy_b;
  Double I_NAI_Hx4z_H3xyz_b = I_NAI_Ix5z_G3xy_b+ABZ*I_NAI_Hx4z_G3xy_b;
  Double I_NAI_H5y_H3xyz_b = I_NAI_I5yz_G3xy_b+ABZ*I_NAI_H5y_G3xy_b;
  Double I_NAI_H4yz_H3xyz_b = I_NAI_I4y2z_G3xy_b+ABZ*I_NAI_H4yz_G3xy_b;
  Double I_NAI_H3y2z_H3xyz_b = I_NAI_I3y3z_G3xy_b+ABZ*I_NAI_H3y2z_G3xy_b;
  Double I_NAI_H2y3z_H3xyz_b = I_NAI_I2y4z_G3xy_b+ABZ*I_NAI_H2y3z_G3xy_b;
  Double I_NAI_Hy4z_H3xyz_b = I_NAI_Iy5z_G3xy_b+ABZ*I_NAI_Hy4z_G3xy_b;
  Double I_NAI_H5z_H3xyz_b = I_NAI_I6z_G3xy_b+ABZ*I_NAI_H5z_G3xy_b;
  Double I_NAI_H5x_H3x2z_b = I_NAI_I5xz_G3xz_b+ABZ*I_NAI_H5x_G3xz_b;
  Double I_NAI_H4xy_H3x2z_b = I_NAI_I4xyz_G3xz_b+ABZ*I_NAI_H4xy_G3xz_b;
  Double I_NAI_H4xz_H3x2z_b = I_NAI_I4x2z_G3xz_b+ABZ*I_NAI_H4xz_G3xz_b;
  Double I_NAI_H3x2y_H3x2z_b = I_NAI_I3x2yz_G3xz_b+ABZ*I_NAI_H3x2y_G3xz_b;
  Double I_NAI_H3xyz_H3x2z_b = I_NAI_I3xy2z_G3xz_b+ABZ*I_NAI_H3xyz_G3xz_b;
  Double I_NAI_H3x2z_H3x2z_b = I_NAI_I3x3z_G3xz_b+ABZ*I_NAI_H3x2z_G3xz_b;
  Double I_NAI_H2x3y_H3x2z_b = I_NAI_I2x3yz_G3xz_b+ABZ*I_NAI_H2x3y_G3xz_b;
  Double I_NAI_H2x2yz_H3x2z_b = I_NAI_I2x2y2z_G3xz_b+ABZ*I_NAI_H2x2yz_G3xz_b;
  Double I_NAI_H2xy2z_H3x2z_b = I_NAI_I2xy3z_G3xz_b+ABZ*I_NAI_H2xy2z_G3xz_b;
  Double I_NAI_H2x3z_H3x2z_b = I_NAI_I2x4z_G3xz_b+ABZ*I_NAI_H2x3z_G3xz_b;
  Double I_NAI_Hx4y_H3x2z_b = I_NAI_Ix4yz_G3xz_b+ABZ*I_NAI_Hx4y_G3xz_b;
  Double I_NAI_Hx3yz_H3x2z_b = I_NAI_Ix3y2z_G3xz_b+ABZ*I_NAI_Hx3yz_G3xz_b;
  Double I_NAI_Hx2y2z_H3x2z_b = I_NAI_Ix2y3z_G3xz_b+ABZ*I_NAI_Hx2y2z_G3xz_b;
  Double I_NAI_Hxy3z_H3x2z_b = I_NAI_Ixy4z_G3xz_b+ABZ*I_NAI_Hxy3z_G3xz_b;
  Double I_NAI_Hx4z_H3x2z_b = I_NAI_Ix5z_G3xz_b+ABZ*I_NAI_Hx4z_G3xz_b;
  Double I_NAI_H5y_H3x2z_b = I_NAI_I5yz_G3xz_b+ABZ*I_NAI_H5y_G3xz_b;
  Double I_NAI_H4yz_H3x2z_b = I_NAI_I4y2z_G3xz_b+ABZ*I_NAI_H4yz_G3xz_b;
  Double I_NAI_H3y2z_H3x2z_b = I_NAI_I3y3z_G3xz_b+ABZ*I_NAI_H3y2z_G3xz_b;
  Double I_NAI_H2y3z_H3x2z_b = I_NAI_I2y4z_G3xz_b+ABZ*I_NAI_H2y3z_G3xz_b;
  Double I_NAI_Hy4z_H3x2z_b = I_NAI_Iy5z_G3xz_b+ABZ*I_NAI_Hy4z_G3xz_b;
  Double I_NAI_H5z_H3x2z_b = I_NAI_I6z_G3xz_b+ABZ*I_NAI_H5z_G3xz_b;
  Double I_NAI_H5x_H2x3y_b = I_NAI_I6x_Gx3y_b+ABX*I_NAI_H5x_Gx3y_b;
  Double I_NAI_H4xy_H2x3y_b = I_NAI_I5xy_Gx3y_b+ABX*I_NAI_H4xy_Gx3y_b;
  Double I_NAI_H4xz_H2x3y_b = I_NAI_I5xz_Gx3y_b+ABX*I_NAI_H4xz_Gx3y_b;
  Double I_NAI_H3x2y_H2x3y_b = I_NAI_I4x2y_Gx3y_b+ABX*I_NAI_H3x2y_Gx3y_b;
  Double I_NAI_H3xyz_H2x3y_b = I_NAI_I4xyz_Gx3y_b+ABX*I_NAI_H3xyz_Gx3y_b;
  Double I_NAI_H3x2z_H2x3y_b = I_NAI_I4x2z_Gx3y_b+ABX*I_NAI_H3x2z_Gx3y_b;
  Double I_NAI_H2x3y_H2x3y_b = I_NAI_I3x3y_Gx3y_b+ABX*I_NAI_H2x3y_Gx3y_b;
  Double I_NAI_H2x2yz_H2x3y_b = I_NAI_I3x2yz_Gx3y_b+ABX*I_NAI_H2x2yz_Gx3y_b;
  Double I_NAI_H2xy2z_H2x3y_b = I_NAI_I3xy2z_Gx3y_b+ABX*I_NAI_H2xy2z_Gx3y_b;
  Double I_NAI_H2x3z_H2x3y_b = I_NAI_I3x3z_Gx3y_b+ABX*I_NAI_H2x3z_Gx3y_b;
  Double I_NAI_Hx4y_H2x3y_b = I_NAI_I2x4y_Gx3y_b+ABX*I_NAI_Hx4y_Gx3y_b;
  Double I_NAI_Hx3yz_H2x3y_b = I_NAI_I2x3yz_Gx3y_b+ABX*I_NAI_Hx3yz_Gx3y_b;
  Double I_NAI_Hx2y2z_H2x3y_b = I_NAI_I2x2y2z_Gx3y_b+ABX*I_NAI_Hx2y2z_Gx3y_b;
  Double I_NAI_Hxy3z_H2x3y_b = I_NAI_I2xy3z_Gx3y_b+ABX*I_NAI_Hxy3z_Gx3y_b;
  Double I_NAI_Hx4z_H2x3y_b = I_NAI_I2x4z_Gx3y_b+ABX*I_NAI_Hx4z_Gx3y_b;
  Double I_NAI_H5y_H2x3y_b = I_NAI_Ix5y_Gx3y_b+ABX*I_NAI_H5y_Gx3y_b;
  Double I_NAI_H4yz_H2x3y_b = I_NAI_Ix4yz_Gx3y_b+ABX*I_NAI_H4yz_Gx3y_b;
  Double I_NAI_H3y2z_H2x3y_b = I_NAI_Ix3y2z_Gx3y_b+ABX*I_NAI_H3y2z_Gx3y_b;
  Double I_NAI_H2y3z_H2x3y_b = I_NAI_Ix2y3z_Gx3y_b+ABX*I_NAI_H2y3z_Gx3y_b;
  Double I_NAI_Hy4z_H2x3y_b = I_NAI_Ixy4z_Gx3y_b+ABX*I_NAI_Hy4z_Gx3y_b;
  Double I_NAI_H5z_H2x3y_b = I_NAI_Ix5z_Gx3y_b+ABX*I_NAI_H5z_Gx3y_b;
  Double I_NAI_H5x_H2x2yz_b = I_NAI_I5xz_G2x2y_b+ABZ*I_NAI_H5x_G2x2y_b;
  Double I_NAI_H4xy_H2x2yz_b = I_NAI_I4xyz_G2x2y_b+ABZ*I_NAI_H4xy_G2x2y_b;
  Double I_NAI_H4xz_H2x2yz_b = I_NAI_I4x2z_G2x2y_b+ABZ*I_NAI_H4xz_G2x2y_b;
  Double I_NAI_H3x2y_H2x2yz_b = I_NAI_I3x2yz_G2x2y_b+ABZ*I_NAI_H3x2y_G2x2y_b;
  Double I_NAI_H3xyz_H2x2yz_b = I_NAI_I3xy2z_G2x2y_b+ABZ*I_NAI_H3xyz_G2x2y_b;
  Double I_NAI_H3x2z_H2x2yz_b = I_NAI_I3x3z_G2x2y_b+ABZ*I_NAI_H3x2z_G2x2y_b;
  Double I_NAI_H2x3y_H2x2yz_b = I_NAI_I2x3yz_G2x2y_b+ABZ*I_NAI_H2x3y_G2x2y_b;
  Double I_NAI_H2x2yz_H2x2yz_b = I_NAI_I2x2y2z_G2x2y_b+ABZ*I_NAI_H2x2yz_G2x2y_b;
  Double I_NAI_H2xy2z_H2x2yz_b = I_NAI_I2xy3z_G2x2y_b+ABZ*I_NAI_H2xy2z_G2x2y_b;
  Double I_NAI_H2x3z_H2x2yz_b = I_NAI_I2x4z_G2x2y_b+ABZ*I_NAI_H2x3z_G2x2y_b;
  Double I_NAI_Hx4y_H2x2yz_b = I_NAI_Ix4yz_G2x2y_b+ABZ*I_NAI_Hx4y_G2x2y_b;
  Double I_NAI_Hx3yz_H2x2yz_b = I_NAI_Ix3y2z_G2x2y_b+ABZ*I_NAI_Hx3yz_G2x2y_b;
  Double I_NAI_Hx2y2z_H2x2yz_b = I_NAI_Ix2y3z_G2x2y_b+ABZ*I_NAI_Hx2y2z_G2x2y_b;
  Double I_NAI_Hxy3z_H2x2yz_b = I_NAI_Ixy4z_G2x2y_b+ABZ*I_NAI_Hxy3z_G2x2y_b;
  Double I_NAI_Hx4z_H2x2yz_b = I_NAI_Ix5z_G2x2y_b+ABZ*I_NAI_Hx4z_G2x2y_b;
  Double I_NAI_H5y_H2x2yz_b = I_NAI_I5yz_G2x2y_b+ABZ*I_NAI_H5y_G2x2y_b;
  Double I_NAI_H4yz_H2x2yz_b = I_NAI_I4y2z_G2x2y_b+ABZ*I_NAI_H4yz_G2x2y_b;
  Double I_NAI_H3y2z_H2x2yz_b = I_NAI_I3y3z_G2x2y_b+ABZ*I_NAI_H3y2z_G2x2y_b;
  Double I_NAI_H2y3z_H2x2yz_b = I_NAI_I2y4z_G2x2y_b+ABZ*I_NAI_H2y3z_G2x2y_b;
  Double I_NAI_Hy4z_H2x2yz_b = I_NAI_Iy5z_G2x2y_b+ABZ*I_NAI_Hy4z_G2x2y_b;
  Double I_NAI_H5z_H2x2yz_b = I_NAI_I6z_G2x2y_b+ABZ*I_NAI_H5z_G2x2y_b;
  Double I_NAI_H5x_H2xy2z_b = I_NAI_I5xy_G2x2z_b+ABY*I_NAI_H5x_G2x2z_b;
  Double I_NAI_H4xy_H2xy2z_b = I_NAI_I4x2y_G2x2z_b+ABY*I_NAI_H4xy_G2x2z_b;
  Double I_NAI_H4xz_H2xy2z_b = I_NAI_I4xyz_G2x2z_b+ABY*I_NAI_H4xz_G2x2z_b;
  Double I_NAI_H3x2y_H2xy2z_b = I_NAI_I3x3y_G2x2z_b+ABY*I_NAI_H3x2y_G2x2z_b;
  Double I_NAI_H3xyz_H2xy2z_b = I_NAI_I3x2yz_G2x2z_b+ABY*I_NAI_H3xyz_G2x2z_b;
  Double I_NAI_H3x2z_H2xy2z_b = I_NAI_I3xy2z_G2x2z_b+ABY*I_NAI_H3x2z_G2x2z_b;
  Double I_NAI_H2x3y_H2xy2z_b = I_NAI_I2x4y_G2x2z_b+ABY*I_NAI_H2x3y_G2x2z_b;
  Double I_NAI_H2x2yz_H2xy2z_b = I_NAI_I2x3yz_G2x2z_b+ABY*I_NAI_H2x2yz_G2x2z_b;
  Double I_NAI_H2xy2z_H2xy2z_b = I_NAI_I2x2y2z_G2x2z_b+ABY*I_NAI_H2xy2z_G2x2z_b;
  Double I_NAI_H2x3z_H2xy2z_b = I_NAI_I2xy3z_G2x2z_b+ABY*I_NAI_H2x3z_G2x2z_b;
  Double I_NAI_Hx4y_H2xy2z_b = I_NAI_Ix5y_G2x2z_b+ABY*I_NAI_Hx4y_G2x2z_b;
  Double I_NAI_Hx3yz_H2xy2z_b = I_NAI_Ix4yz_G2x2z_b+ABY*I_NAI_Hx3yz_G2x2z_b;
  Double I_NAI_Hx2y2z_H2xy2z_b = I_NAI_Ix3y2z_G2x2z_b+ABY*I_NAI_Hx2y2z_G2x2z_b;
  Double I_NAI_Hxy3z_H2xy2z_b = I_NAI_Ix2y3z_G2x2z_b+ABY*I_NAI_Hxy3z_G2x2z_b;
  Double I_NAI_Hx4z_H2xy2z_b = I_NAI_Ixy4z_G2x2z_b+ABY*I_NAI_Hx4z_G2x2z_b;
  Double I_NAI_H5y_H2xy2z_b = I_NAI_I6y_G2x2z_b+ABY*I_NAI_H5y_G2x2z_b;
  Double I_NAI_H4yz_H2xy2z_b = I_NAI_I5yz_G2x2z_b+ABY*I_NAI_H4yz_G2x2z_b;
  Double I_NAI_H3y2z_H2xy2z_b = I_NAI_I4y2z_G2x2z_b+ABY*I_NAI_H3y2z_G2x2z_b;
  Double I_NAI_H2y3z_H2xy2z_b = I_NAI_I3y3z_G2x2z_b+ABY*I_NAI_H2y3z_G2x2z_b;
  Double I_NAI_Hy4z_H2xy2z_b = I_NAI_I2y4z_G2x2z_b+ABY*I_NAI_Hy4z_G2x2z_b;
  Double I_NAI_H5z_H2xy2z_b = I_NAI_Iy5z_G2x2z_b+ABY*I_NAI_H5z_G2x2z_b;
  Double I_NAI_H5x_H2x3z_b = I_NAI_I6x_Gx3z_b+ABX*I_NAI_H5x_Gx3z_b;
  Double I_NAI_H4xy_H2x3z_b = I_NAI_I5xy_Gx3z_b+ABX*I_NAI_H4xy_Gx3z_b;
  Double I_NAI_H4xz_H2x3z_b = I_NAI_I5xz_Gx3z_b+ABX*I_NAI_H4xz_Gx3z_b;
  Double I_NAI_H3x2y_H2x3z_b = I_NAI_I4x2y_Gx3z_b+ABX*I_NAI_H3x2y_Gx3z_b;
  Double I_NAI_H3xyz_H2x3z_b = I_NAI_I4xyz_Gx3z_b+ABX*I_NAI_H3xyz_Gx3z_b;
  Double I_NAI_H3x2z_H2x3z_b = I_NAI_I4x2z_Gx3z_b+ABX*I_NAI_H3x2z_Gx3z_b;
  Double I_NAI_H2x3y_H2x3z_b = I_NAI_I3x3y_Gx3z_b+ABX*I_NAI_H2x3y_Gx3z_b;
  Double I_NAI_H2x2yz_H2x3z_b = I_NAI_I3x2yz_Gx3z_b+ABX*I_NAI_H2x2yz_Gx3z_b;
  Double I_NAI_H2xy2z_H2x3z_b = I_NAI_I3xy2z_Gx3z_b+ABX*I_NAI_H2xy2z_Gx3z_b;
  Double I_NAI_H2x3z_H2x3z_b = I_NAI_I3x3z_Gx3z_b+ABX*I_NAI_H2x3z_Gx3z_b;
  Double I_NAI_Hx4y_H2x3z_b = I_NAI_I2x4y_Gx3z_b+ABX*I_NAI_Hx4y_Gx3z_b;
  Double I_NAI_Hx3yz_H2x3z_b = I_NAI_I2x3yz_Gx3z_b+ABX*I_NAI_Hx3yz_Gx3z_b;
  Double I_NAI_Hx2y2z_H2x3z_b = I_NAI_I2x2y2z_Gx3z_b+ABX*I_NAI_Hx2y2z_Gx3z_b;
  Double I_NAI_Hxy3z_H2x3z_b = I_NAI_I2xy3z_Gx3z_b+ABX*I_NAI_Hxy3z_Gx3z_b;
  Double I_NAI_Hx4z_H2x3z_b = I_NAI_I2x4z_Gx3z_b+ABX*I_NAI_Hx4z_Gx3z_b;
  Double I_NAI_H5y_H2x3z_b = I_NAI_Ix5y_Gx3z_b+ABX*I_NAI_H5y_Gx3z_b;
  Double I_NAI_H4yz_H2x3z_b = I_NAI_Ix4yz_Gx3z_b+ABX*I_NAI_H4yz_Gx3z_b;
  Double I_NAI_H3y2z_H2x3z_b = I_NAI_Ix3y2z_Gx3z_b+ABX*I_NAI_H3y2z_Gx3z_b;
  Double I_NAI_H2y3z_H2x3z_b = I_NAI_Ix2y3z_Gx3z_b+ABX*I_NAI_H2y3z_Gx3z_b;
  Double I_NAI_Hy4z_H2x3z_b = I_NAI_Ixy4z_Gx3z_b+ABX*I_NAI_Hy4z_Gx3z_b;
  Double I_NAI_H5z_H2x3z_b = I_NAI_Ix5z_Gx3z_b+ABX*I_NAI_H5z_Gx3z_b;
  Double I_NAI_H5x_Hx4y_b = I_NAI_I6x_G4y_b+ABX*I_NAI_H5x_G4y_b;
  Double I_NAI_H4xy_Hx4y_b = I_NAI_I5xy_G4y_b+ABX*I_NAI_H4xy_G4y_b;
  Double I_NAI_H4xz_Hx4y_b = I_NAI_I5xz_G4y_b+ABX*I_NAI_H4xz_G4y_b;
  Double I_NAI_H3x2y_Hx4y_b = I_NAI_I4x2y_G4y_b+ABX*I_NAI_H3x2y_G4y_b;
  Double I_NAI_H3xyz_Hx4y_b = I_NAI_I4xyz_G4y_b+ABX*I_NAI_H3xyz_G4y_b;
  Double I_NAI_H3x2z_Hx4y_b = I_NAI_I4x2z_G4y_b+ABX*I_NAI_H3x2z_G4y_b;
  Double I_NAI_H2x3y_Hx4y_b = I_NAI_I3x3y_G4y_b+ABX*I_NAI_H2x3y_G4y_b;
  Double I_NAI_H2x2yz_Hx4y_b = I_NAI_I3x2yz_G4y_b+ABX*I_NAI_H2x2yz_G4y_b;
  Double I_NAI_H2xy2z_Hx4y_b = I_NAI_I3xy2z_G4y_b+ABX*I_NAI_H2xy2z_G4y_b;
  Double I_NAI_H2x3z_Hx4y_b = I_NAI_I3x3z_G4y_b+ABX*I_NAI_H2x3z_G4y_b;
  Double I_NAI_Hx4y_Hx4y_b = I_NAI_I2x4y_G4y_b+ABX*I_NAI_Hx4y_G4y_b;
  Double I_NAI_Hx3yz_Hx4y_b = I_NAI_I2x3yz_G4y_b+ABX*I_NAI_Hx3yz_G4y_b;
  Double I_NAI_Hx2y2z_Hx4y_b = I_NAI_I2x2y2z_G4y_b+ABX*I_NAI_Hx2y2z_G4y_b;
  Double I_NAI_Hxy3z_Hx4y_b = I_NAI_I2xy3z_G4y_b+ABX*I_NAI_Hxy3z_G4y_b;
  Double I_NAI_Hx4z_Hx4y_b = I_NAI_I2x4z_G4y_b+ABX*I_NAI_Hx4z_G4y_b;
  Double I_NAI_H5y_Hx4y_b = I_NAI_Ix5y_G4y_b+ABX*I_NAI_H5y_G4y_b;
  Double I_NAI_H4yz_Hx4y_b = I_NAI_Ix4yz_G4y_b+ABX*I_NAI_H4yz_G4y_b;
  Double I_NAI_H3y2z_Hx4y_b = I_NAI_Ix3y2z_G4y_b+ABX*I_NAI_H3y2z_G4y_b;
  Double I_NAI_H2y3z_Hx4y_b = I_NAI_Ix2y3z_G4y_b+ABX*I_NAI_H2y3z_G4y_b;
  Double I_NAI_Hy4z_Hx4y_b = I_NAI_Ixy4z_G4y_b+ABX*I_NAI_Hy4z_G4y_b;
  Double I_NAI_H5z_Hx4y_b = I_NAI_Ix5z_G4y_b+ABX*I_NAI_H5z_G4y_b;
  Double I_NAI_H5x_Hx3yz_b = I_NAI_I5xz_Gx3y_b+ABZ*I_NAI_H5x_Gx3y_b;
  Double I_NAI_H4xy_Hx3yz_b = I_NAI_I4xyz_Gx3y_b+ABZ*I_NAI_H4xy_Gx3y_b;
  Double I_NAI_H4xz_Hx3yz_b = I_NAI_I4x2z_Gx3y_b+ABZ*I_NAI_H4xz_Gx3y_b;
  Double I_NAI_H3x2y_Hx3yz_b = I_NAI_I3x2yz_Gx3y_b+ABZ*I_NAI_H3x2y_Gx3y_b;
  Double I_NAI_H3xyz_Hx3yz_b = I_NAI_I3xy2z_Gx3y_b+ABZ*I_NAI_H3xyz_Gx3y_b;
  Double I_NAI_H3x2z_Hx3yz_b = I_NAI_I3x3z_Gx3y_b+ABZ*I_NAI_H3x2z_Gx3y_b;
  Double I_NAI_H2x3y_Hx3yz_b = I_NAI_I2x3yz_Gx3y_b+ABZ*I_NAI_H2x3y_Gx3y_b;
  Double I_NAI_H2x2yz_Hx3yz_b = I_NAI_I2x2y2z_Gx3y_b+ABZ*I_NAI_H2x2yz_Gx3y_b;
  Double I_NAI_H2xy2z_Hx3yz_b = I_NAI_I2xy3z_Gx3y_b+ABZ*I_NAI_H2xy2z_Gx3y_b;
  Double I_NAI_H2x3z_Hx3yz_b = I_NAI_I2x4z_Gx3y_b+ABZ*I_NAI_H2x3z_Gx3y_b;
  Double I_NAI_Hx4y_Hx3yz_b = I_NAI_Ix4yz_Gx3y_b+ABZ*I_NAI_Hx4y_Gx3y_b;
  Double I_NAI_Hx3yz_Hx3yz_b = I_NAI_Ix3y2z_Gx3y_b+ABZ*I_NAI_Hx3yz_Gx3y_b;
  Double I_NAI_Hx2y2z_Hx3yz_b = I_NAI_Ix2y3z_Gx3y_b+ABZ*I_NAI_Hx2y2z_Gx3y_b;
  Double I_NAI_Hxy3z_Hx3yz_b = I_NAI_Ixy4z_Gx3y_b+ABZ*I_NAI_Hxy3z_Gx3y_b;
  Double I_NAI_Hx4z_Hx3yz_b = I_NAI_Ix5z_Gx3y_b+ABZ*I_NAI_Hx4z_Gx3y_b;
  Double I_NAI_H5y_Hx3yz_b = I_NAI_I5yz_Gx3y_b+ABZ*I_NAI_H5y_Gx3y_b;
  Double I_NAI_H4yz_Hx3yz_b = I_NAI_I4y2z_Gx3y_b+ABZ*I_NAI_H4yz_Gx3y_b;
  Double I_NAI_H3y2z_Hx3yz_b = I_NAI_I3y3z_Gx3y_b+ABZ*I_NAI_H3y2z_Gx3y_b;
  Double I_NAI_H2y3z_Hx3yz_b = I_NAI_I2y4z_Gx3y_b+ABZ*I_NAI_H2y3z_Gx3y_b;
  Double I_NAI_Hy4z_Hx3yz_b = I_NAI_Iy5z_Gx3y_b+ABZ*I_NAI_Hy4z_Gx3y_b;
  Double I_NAI_H5z_Hx3yz_b = I_NAI_I6z_Gx3y_b+ABZ*I_NAI_H5z_Gx3y_b;
  Double I_NAI_H5x_Hx2y2z_b = I_NAI_I6x_G2y2z_b+ABX*I_NAI_H5x_G2y2z_b;
  Double I_NAI_H4xy_Hx2y2z_b = I_NAI_I5xy_G2y2z_b+ABX*I_NAI_H4xy_G2y2z_b;
  Double I_NAI_H4xz_Hx2y2z_b = I_NAI_I5xz_G2y2z_b+ABX*I_NAI_H4xz_G2y2z_b;
  Double I_NAI_H3x2y_Hx2y2z_b = I_NAI_I4x2y_G2y2z_b+ABX*I_NAI_H3x2y_G2y2z_b;
  Double I_NAI_H3xyz_Hx2y2z_b = I_NAI_I4xyz_G2y2z_b+ABX*I_NAI_H3xyz_G2y2z_b;
  Double I_NAI_H3x2z_Hx2y2z_b = I_NAI_I4x2z_G2y2z_b+ABX*I_NAI_H3x2z_G2y2z_b;
  Double I_NAI_H2x3y_Hx2y2z_b = I_NAI_I3x3y_G2y2z_b+ABX*I_NAI_H2x3y_G2y2z_b;
  Double I_NAI_H2x2yz_Hx2y2z_b = I_NAI_I3x2yz_G2y2z_b+ABX*I_NAI_H2x2yz_G2y2z_b;
  Double I_NAI_H2xy2z_Hx2y2z_b = I_NAI_I3xy2z_G2y2z_b+ABX*I_NAI_H2xy2z_G2y2z_b;
  Double I_NAI_H2x3z_Hx2y2z_b = I_NAI_I3x3z_G2y2z_b+ABX*I_NAI_H2x3z_G2y2z_b;
  Double I_NAI_Hx4y_Hx2y2z_b = I_NAI_I2x4y_G2y2z_b+ABX*I_NAI_Hx4y_G2y2z_b;
  Double I_NAI_Hx3yz_Hx2y2z_b = I_NAI_I2x3yz_G2y2z_b+ABX*I_NAI_Hx3yz_G2y2z_b;
  Double I_NAI_Hx2y2z_Hx2y2z_b = I_NAI_I2x2y2z_G2y2z_b+ABX*I_NAI_Hx2y2z_G2y2z_b;
  Double I_NAI_Hxy3z_Hx2y2z_b = I_NAI_I2xy3z_G2y2z_b+ABX*I_NAI_Hxy3z_G2y2z_b;
  Double I_NAI_Hx4z_Hx2y2z_b = I_NAI_I2x4z_G2y2z_b+ABX*I_NAI_Hx4z_G2y2z_b;
  Double I_NAI_H5y_Hx2y2z_b = I_NAI_Ix5y_G2y2z_b+ABX*I_NAI_H5y_G2y2z_b;
  Double I_NAI_H4yz_Hx2y2z_b = I_NAI_Ix4yz_G2y2z_b+ABX*I_NAI_H4yz_G2y2z_b;
  Double I_NAI_H3y2z_Hx2y2z_b = I_NAI_Ix3y2z_G2y2z_b+ABX*I_NAI_H3y2z_G2y2z_b;
  Double I_NAI_H2y3z_Hx2y2z_b = I_NAI_Ix2y3z_G2y2z_b+ABX*I_NAI_H2y3z_G2y2z_b;
  Double I_NAI_Hy4z_Hx2y2z_b = I_NAI_Ixy4z_G2y2z_b+ABX*I_NAI_Hy4z_G2y2z_b;
  Double I_NAI_H5z_Hx2y2z_b = I_NAI_Ix5z_G2y2z_b+ABX*I_NAI_H5z_G2y2z_b;
  Double I_NAI_H5x_Hxy3z_b = I_NAI_I5xy_Gx3z_b+ABY*I_NAI_H5x_Gx3z_b;
  Double I_NAI_H4xy_Hxy3z_b = I_NAI_I4x2y_Gx3z_b+ABY*I_NAI_H4xy_Gx3z_b;
  Double I_NAI_H4xz_Hxy3z_b = I_NAI_I4xyz_Gx3z_b+ABY*I_NAI_H4xz_Gx3z_b;
  Double I_NAI_H3x2y_Hxy3z_b = I_NAI_I3x3y_Gx3z_b+ABY*I_NAI_H3x2y_Gx3z_b;
  Double I_NAI_H3xyz_Hxy3z_b = I_NAI_I3x2yz_Gx3z_b+ABY*I_NAI_H3xyz_Gx3z_b;
  Double I_NAI_H3x2z_Hxy3z_b = I_NAI_I3xy2z_Gx3z_b+ABY*I_NAI_H3x2z_Gx3z_b;
  Double I_NAI_H2x3y_Hxy3z_b = I_NAI_I2x4y_Gx3z_b+ABY*I_NAI_H2x3y_Gx3z_b;
  Double I_NAI_H2x2yz_Hxy3z_b = I_NAI_I2x3yz_Gx3z_b+ABY*I_NAI_H2x2yz_Gx3z_b;
  Double I_NAI_H2xy2z_Hxy3z_b = I_NAI_I2x2y2z_Gx3z_b+ABY*I_NAI_H2xy2z_Gx3z_b;
  Double I_NAI_H2x3z_Hxy3z_b = I_NAI_I2xy3z_Gx3z_b+ABY*I_NAI_H2x3z_Gx3z_b;
  Double I_NAI_Hx4y_Hxy3z_b = I_NAI_Ix5y_Gx3z_b+ABY*I_NAI_Hx4y_Gx3z_b;
  Double I_NAI_Hx3yz_Hxy3z_b = I_NAI_Ix4yz_Gx3z_b+ABY*I_NAI_Hx3yz_Gx3z_b;
  Double I_NAI_Hx2y2z_Hxy3z_b = I_NAI_Ix3y2z_Gx3z_b+ABY*I_NAI_Hx2y2z_Gx3z_b;
  Double I_NAI_Hxy3z_Hxy3z_b = I_NAI_Ix2y3z_Gx3z_b+ABY*I_NAI_Hxy3z_Gx3z_b;
  Double I_NAI_Hx4z_Hxy3z_b = I_NAI_Ixy4z_Gx3z_b+ABY*I_NAI_Hx4z_Gx3z_b;
  Double I_NAI_H5y_Hxy3z_b = I_NAI_I6y_Gx3z_b+ABY*I_NAI_H5y_Gx3z_b;
  Double I_NAI_H4yz_Hxy3z_b = I_NAI_I5yz_Gx3z_b+ABY*I_NAI_H4yz_Gx3z_b;
  Double I_NAI_H3y2z_Hxy3z_b = I_NAI_I4y2z_Gx3z_b+ABY*I_NAI_H3y2z_Gx3z_b;
  Double I_NAI_H2y3z_Hxy3z_b = I_NAI_I3y3z_Gx3z_b+ABY*I_NAI_H2y3z_Gx3z_b;
  Double I_NAI_Hy4z_Hxy3z_b = I_NAI_I2y4z_Gx3z_b+ABY*I_NAI_Hy4z_Gx3z_b;
  Double I_NAI_H5z_Hxy3z_b = I_NAI_Iy5z_Gx3z_b+ABY*I_NAI_H5z_Gx3z_b;
  Double I_NAI_H5x_Hx4z_b = I_NAI_I6x_G4z_b+ABX*I_NAI_H5x_G4z_b;
  Double I_NAI_H4xy_Hx4z_b = I_NAI_I5xy_G4z_b+ABX*I_NAI_H4xy_G4z_b;
  Double I_NAI_H4xz_Hx4z_b = I_NAI_I5xz_G4z_b+ABX*I_NAI_H4xz_G4z_b;
  Double I_NAI_H3x2y_Hx4z_b = I_NAI_I4x2y_G4z_b+ABX*I_NAI_H3x2y_G4z_b;
  Double I_NAI_H3xyz_Hx4z_b = I_NAI_I4xyz_G4z_b+ABX*I_NAI_H3xyz_G4z_b;
  Double I_NAI_H3x2z_Hx4z_b = I_NAI_I4x2z_G4z_b+ABX*I_NAI_H3x2z_G4z_b;
  Double I_NAI_H2x3y_Hx4z_b = I_NAI_I3x3y_G4z_b+ABX*I_NAI_H2x3y_G4z_b;
  Double I_NAI_H2x2yz_Hx4z_b = I_NAI_I3x2yz_G4z_b+ABX*I_NAI_H2x2yz_G4z_b;
  Double I_NAI_H2xy2z_Hx4z_b = I_NAI_I3xy2z_G4z_b+ABX*I_NAI_H2xy2z_G4z_b;
  Double I_NAI_H2x3z_Hx4z_b = I_NAI_I3x3z_G4z_b+ABX*I_NAI_H2x3z_G4z_b;
  Double I_NAI_Hx4y_Hx4z_b = I_NAI_I2x4y_G4z_b+ABX*I_NAI_Hx4y_G4z_b;
  Double I_NAI_Hx3yz_Hx4z_b = I_NAI_I2x3yz_G4z_b+ABX*I_NAI_Hx3yz_G4z_b;
  Double I_NAI_Hx2y2z_Hx4z_b = I_NAI_I2x2y2z_G4z_b+ABX*I_NAI_Hx2y2z_G4z_b;
  Double I_NAI_Hxy3z_Hx4z_b = I_NAI_I2xy3z_G4z_b+ABX*I_NAI_Hxy3z_G4z_b;
  Double I_NAI_Hx4z_Hx4z_b = I_NAI_I2x4z_G4z_b+ABX*I_NAI_Hx4z_G4z_b;
  Double I_NAI_H5y_Hx4z_b = I_NAI_Ix5y_G4z_b+ABX*I_NAI_H5y_G4z_b;
  Double I_NAI_H4yz_Hx4z_b = I_NAI_Ix4yz_G4z_b+ABX*I_NAI_H4yz_G4z_b;
  Double I_NAI_H3y2z_Hx4z_b = I_NAI_Ix3y2z_G4z_b+ABX*I_NAI_H3y2z_G4z_b;
  Double I_NAI_H2y3z_Hx4z_b = I_NAI_Ix2y3z_G4z_b+ABX*I_NAI_H2y3z_G4z_b;
  Double I_NAI_Hy4z_Hx4z_b = I_NAI_Ixy4z_G4z_b+ABX*I_NAI_Hy4z_G4z_b;
  Double I_NAI_H5z_Hx4z_b = I_NAI_Ix5z_G4z_b+ABX*I_NAI_H5z_G4z_b;
  Double I_NAI_H5x_H5y_b = I_NAI_I5xy_G4y_b+ABY*I_NAI_H5x_G4y_b;
  Double I_NAI_H4xy_H5y_b = I_NAI_I4x2y_G4y_b+ABY*I_NAI_H4xy_G4y_b;
  Double I_NAI_H4xz_H5y_b = I_NAI_I4xyz_G4y_b+ABY*I_NAI_H4xz_G4y_b;
  Double I_NAI_H3x2y_H5y_b = I_NAI_I3x3y_G4y_b+ABY*I_NAI_H3x2y_G4y_b;
  Double I_NAI_H3xyz_H5y_b = I_NAI_I3x2yz_G4y_b+ABY*I_NAI_H3xyz_G4y_b;
  Double I_NAI_H3x2z_H5y_b = I_NAI_I3xy2z_G4y_b+ABY*I_NAI_H3x2z_G4y_b;
  Double I_NAI_H2x3y_H5y_b = I_NAI_I2x4y_G4y_b+ABY*I_NAI_H2x3y_G4y_b;
  Double I_NAI_H2x2yz_H5y_b = I_NAI_I2x3yz_G4y_b+ABY*I_NAI_H2x2yz_G4y_b;
  Double I_NAI_H2xy2z_H5y_b = I_NAI_I2x2y2z_G4y_b+ABY*I_NAI_H2xy2z_G4y_b;
  Double I_NAI_H2x3z_H5y_b = I_NAI_I2xy3z_G4y_b+ABY*I_NAI_H2x3z_G4y_b;
  Double I_NAI_Hx4y_H5y_b = I_NAI_Ix5y_G4y_b+ABY*I_NAI_Hx4y_G4y_b;
  Double I_NAI_Hx3yz_H5y_b = I_NAI_Ix4yz_G4y_b+ABY*I_NAI_Hx3yz_G4y_b;
  Double I_NAI_Hx2y2z_H5y_b = I_NAI_Ix3y2z_G4y_b+ABY*I_NAI_Hx2y2z_G4y_b;
  Double I_NAI_Hxy3z_H5y_b = I_NAI_Ix2y3z_G4y_b+ABY*I_NAI_Hxy3z_G4y_b;
  Double I_NAI_Hx4z_H5y_b = I_NAI_Ixy4z_G4y_b+ABY*I_NAI_Hx4z_G4y_b;
  Double I_NAI_H5y_H5y_b = I_NAI_I6y_G4y_b+ABY*I_NAI_H5y_G4y_b;
  Double I_NAI_H4yz_H5y_b = I_NAI_I5yz_G4y_b+ABY*I_NAI_H4yz_G4y_b;
  Double I_NAI_H3y2z_H5y_b = I_NAI_I4y2z_G4y_b+ABY*I_NAI_H3y2z_G4y_b;
  Double I_NAI_H2y3z_H5y_b = I_NAI_I3y3z_G4y_b+ABY*I_NAI_H2y3z_G4y_b;
  Double I_NAI_Hy4z_H5y_b = I_NAI_I2y4z_G4y_b+ABY*I_NAI_Hy4z_G4y_b;
  Double I_NAI_H5z_H5y_b = I_NAI_Iy5z_G4y_b+ABY*I_NAI_H5z_G4y_b;
  Double I_NAI_H5x_H4yz_b = I_NAI_I5xz_G4y_b+ABZ*I_NAI_H5x_G4y_b;
  Double I_NAI_H4xy_H4yz_b = I_NAI_I4xyz_G4y_b+ABZ*I_NAI_H4xy_G4y_b;
  Double I_NAI_H4xz_H4yz_b = I_NAI_I4x2z_G4y_b+ABZ*I_NAI_H4xz_G4y_b;
  Double I_NAI_H3x2y_H4yz_b = I_NAI_I3x2yz_G4y_b+ABZ*I_NAI_H3x2y_G4y_b;
  Double I_NAI_H3xyz_H4yz_b = I_NAI_I3xy2z_G4y_b+ABZ*I_NAI_H3xyz_G4y_b;
  Double I_NAI_H3x2z_H4yz_b = I_NAI_I3x3z_G4y_b+ABZ*I_NAI_H3x2z_G4y_b;
  Double I_NAI_H2x3y_H4yz_b = I_NAI_I2x3yz_G4y_b+ABZ*I_NAI_H2x3y_G4y_b;
  Double I_NAI_H2x2yz_H4yz_b = I_NAI_I2x2y2z_G4y_b+ABZ*I_NAI_H2x2yz_G4y_b;
  Double I_NAI_H2xy2z_H4yz_b = I_NAI_I2xy3z_G4y_b+ABZ*I_NAI_H2xy2z_G4y_b;
  Double I_NAI_H2x3z_H4yz_b = I_NAI_I2x4z_G4y_b+ABZ*I_NAI_H2x3z_G4y_b;
  Double I_NAI_Hx4y_H4yz_b = I_NAI_Ix4yz_G4y_b+ABZ*I_NAI_Hx4y_G4y_b;
  Double I_NAI_Hx3yz_H4yz_b = I_NAI_Ix3y2z_G4y_b+ABZ*I_NAI_Hx3yz_G4y_b;
  Double I_NAI_Hx2y2z_H4yz_b = I_NAI_Ix2y3z_G4y_b+ABZ*I_NAI_Hx2y2z_G4y_b;
  Double I_NAI_Hxy3z_H4yz_b = I_NAI_Ixy4z_G4y_b+ABZ*I_NAI_Hxy3z_G4y_b;
  Double I_NAI_Hx4z_H4yz_b = I_NAI_Ix5z_G4y_b+ABZ*I_NAI_Hx4z_G4y_b;
  Double I_NAI_H5y_H4yz_b = I_NAI_I5yz_G4y_b+ABZ*I_NAI_H5y_G4y_b;
  Double I_NAI_H4yz_H4yz_b = I_NAI_I4y2z_G4y_b+ABZ*I_NAI_H4yz_G4y_b;
  Double I_NAI_H3y2z_H4yz_b = I_NAI_I3y3z_G4y_b+ABZ*I_NAI_H3y2z_G4y_b;
  Double I_NAI_H2y3z_H4yz_b = I_NAI_I2y4z_G4y_b+ABZ*I_NAI_H2y3z_G4y_b;
  Double I_NAI_Hy4z_H4yz_b = I_NAI_Iy5z_G4y_b+ABZ*I_NAI_Hy4z_G4y_b;
  Double I_NAI_H5z_H4yz_b = I_NAI_I6z_G4y_b+ABZ*I_NAI_H5z_G4y_b;
  Double I_NAI_H5x_H3y2z_b = I_NAI_I5xz_G3yz_b+ABZ*I_NAI_H5x_G3yz_b;
  Double I_NAI_H4xy_H3y2z_b = I_NAI_I4xyz_G3yz_b+ABZ*I_NAI_H4xy_G3yz_b;
  Double I_NAI_H4xz_H3y2z_b = I_NAI_I4x2z_G3yz_b+ABZ*I_NAI_H4xz_G3yz_b;
  Double I_NAI_H3x2y_H3y2z_b = I_NAI_I3x2yz_G3yz_b+ABZ*I_NAI_H3x2y_G3yz_b;
  Double I_NAI_H3xyz_H3y2z_b = I_NAI_I3xy2z_G3yz_b+ABZ*I_NAI_H3xyz_G3yz_b;
  Double I_NAI_H3x2z_H3y2z_b = I_NAI_I3x3z_G3yz_b+ABZ*I_NAI_H3x2z_G3yz_b;
  Double I_NAI_H2x3y_H3y2z_b = I_NAI_I2x3yz_G3yz_b+ABZ*I_NAI_H2x3y_G3yz_b;
  Double I_NAI_H2x2yz_H3y2z_b = I_NAI_I2x2y2z_G3yz_b+ABZ*I_NAI_H2x2yz_G3yz_b;
  Double I_NAI_H2xy2z_H3y2z_b = I_NAI_I2xy3z_G3yz_b+ABZ*I_NAI_H2xy2z_G3yz_b;
  Double I_NAI_H2x3z_H3y2z_b = I_NAI_I2x4z_G3yz_b+ABZ*I_NAI_H2x3z_G3yz_b;
  Double I_NAI_Hx4y_H3y2z_b = I_NAI_Ix4yz_G3yz_b+ABZ*I_NAI_Hx4y_G3yz_b;
  Double I_NAI_Hx3yz_H3y2z_b = I_NAI_Ix3y2z_G3yz_b+ABZ*I_NAI_Hx3yz_G3yz_b;
  Double I_NAI_Hx2y2z_H3y2z_b = I_NAI_Ix2y3z_G3yz_b+ABZ*I_NAI_Hx2y2z_G3yz_b;
  Double I_NAI_Hxy3z_H3y2z_b = I_NAI_Ixy4z_G3yz_b+ABZ*I_NAI_Hxy3z_G3yz_b;
  Double I_NAI_Hx4z_H3y2z_b = I_NAI_Ix5z_G3yz_b+ABZ*I_NAI_Hx4z_G3yz_b;
  Double I_NAI_H5y_H3y2z_b = I_NAI_I5yz_G3yz_b+ABZ*I_NAI_H5y_G3yz_b;
  Double I_NAI_H4yz_H3y2z_b = I_NAI_I4y2z_G3yz_b+ABZ*I_NAI_H4yz_G3yz_b;
  Double I_NAI_H3y2z_H3y2z_b = I_NAI_I3y3z_G3yz_b+ABZ*I_NAI_H3y2z_G3yz_b;
  Double I_NAI_H2y3z_H3y2z_b = I_NAI_I2y4z_G3yz_b+ABZ*I_NAI_H2y3z_G3yz_b;
  Double I_NAI_Hy4z_H3y2z_b = I_NAI_Iy5z_G3yz_b+ABZ*I_NAI_Hy4z_G3yz_b;
  Double I_NAI_H5z_H3y2z_b = I_NAI_I6z_G3yz_b+ABZ*I_NAI_H5z_G3yz_b;
  Double I_NAI_H5x_H2y3z_b = I_NAI_I5xy_Gy3z_b+ABY*I_NAI_H5x_Gy3z_b;
  Double I_NAI_H4xy_H2y3z_b = I_NAI_I4x2y_Gy3z_b+ABY*I_NAI_H4xy_Gy3z_b;
  Double I_NAI_H4xz_H2y3z_b = I_NAI_I4xyz_Gy3z_b+ABY*I_NAI_H4xz_Gy3z_b;
  Double I_NAI_H3x2y_H2y3z_b = I_NAI_I3x3y_Gy3z_b+ABY*I_NAI_H3x2y_Gy3z_b;
  Double I_NAI_H3xyz_H2y3z_b = I_NAI_I3x2yz_Gy3z_b+ABY*I_NAI_H3xyz_Gy3z_b;
  Double I_NAI_H3x2z_H2y3z_b = I_NAI_I3xy2z_Gy3z_b+ABY*I_NAI_H3x2z_Gy3z_b;
  Double I_NAI_H2x3y_H2y3z_b = I_NAI_I2x4y_Gy3z_b+ABY*I_NAI_H2x3y_Gy3z_b;
  Double I_NAI_H2x2yz_H2y3z_b = I_NAI_I2x3yz_Gy3z_b+ABY*I_NAI_H2x2yz_Gy3z_b;
  Double I_NAI_H2xy2z_H2y3z_b = I_NAI_I2x2y2z_Gy3z_b+ABY*I_NAI_H2xy2z_Gy3z_b;
  Double I_NAI_H2x3z_H2y3z_b = I_NAI_I2xy3z_Gy3z_b+ABY*I_NAI_H2x3z_Gy3z_b;
  Double I_NAI_Hx4y_H2y3z_b = I_NAI_Ix5y_Gy3z_b+ABY*I_NAI_Hx4y_Gy3z_b;
  Double I_NAI_Hx3yz_H2y3z_b = I_NAI_Ix4yz_Gy3z_b+ABY*I_NAI_Hx3yz_Gy3z_b;
  Double I_NAI_Hx2y2z_H2y3z_b = I_NAI_Ix3y2z_Gy3z_b+ABY*I_NAI_Hx2y2z_Gy3z_b;
  Double I_NAI_Hxy3z_H2y3z_b = I_NAI_Ix2y3z_Gy3z_b+ABY*I_NAI_Hxy3z_Gy3z_b;
  Double I_NAI_Hx4z_H2y3z_b = I_NAI_Ixy4z_Gy3z_b+ABY*I_NAI_Hx4z_Gy3z_b;
  Double I_NAI_H5y_H2y3z_b = I_NAI_I6y_Gy3z_b+ABY*I_NAI_H5y_Gy3z_b;
  Double I_NAI_H4yz_H2y3z_b = I_NAI_I5yz_Gy3z_b+ABY*I_NAI_H4yz_Gy3z_b;
  Double I_NAI_H3y2z_H2y3z_b = I_NAI_I4y2z_Gy3z_b+ABY*I_NAI_H3y2z_Gy3z_b;
  Double I_NAI_H2y3z_H2y3z_b = I_NAI_I3y3z_Gy3z_b+ABY*I_NAI_H2y3z_Gy3z_b;
  Double I_NAI_Hy4z_H2y3z_b = I_NAI_I2y4z_Gy3z_b+ABY*I_NAI_Hy4z_Gy3z_b;
  Double I_NAI_H5z_H2y3z_b = I_NAI_Iy5z_Gy3z_b+ABY*I_NAI_H5z_Gy3z_b;
  Double I_NAI_H5x_Hy4z_b = I_NAI_I5xy_G4z_b+ABY*I_NAI_H5x_G4z_b;
  Double I_NAI_H4xy_Hy4z_b = I_NAI_I4x2y_G4z_b+ABY*I_NAI_H4xy_G4z_b;
  Double I_NAI_H4xz_Hy4z_b = I_NAI_I4xyz_G4z_b+ABY*I_NAI_H4xz_G4z_b;
  Double I_NAI_H3x2y_Hy4z_b = I_NAI_I3x3y_G4z_b+ABY*I_NAI_H3x2y_G4z_b;
  Double I_NAI_H3xyz_Hy4z_b = I_NAI_I3x2yz_G4z_b+ABY*I_NAI_H3xyz_G4z_b;
  Double I_NAI_H3x2z_Hy4z_b = I_NAI_I3xy2z_G4z_b+ABY*I_NAI_H3x2z_G4z_b;
  Double I_NAI_H2x3y_Hy4z_b = I_NAI_I2x4y_G4z_b+ABY*I_NAI_H2x3y_G4z_b;
  Double I_NAI_H2x2yz_Hy4z_b = I_NAI_I2x3yz_G4z_b+ABY*I_NAI_H2x2yz_G4z_b;
  Double I_NAI_H2xy2z_Hy4z_b = I_NAI_I2x2y2z_G4z_b+ABY*I_NAI_H2xy2z_G4z_b;
  Double I_NAI_H2x3z_Hy4z_b = I_NAI_I2xy3z_G4z_b+ABY*I_NAI_H2x3z_G4z_b;
  Double I_NAI_Hx4y_Hy4z_b = I_NAI_Ix5y_G4z_b+ABY*I_NAI_Hx4y_G4z_b;
  Double I_NAI_Hx3yz_Hy4z_b = I_NAI_Ix4yz_G4z_b+ABY*I_NAI_Hx3yz_G4z_b;
  Double I_NAI_Hx2y2z_Hy4z_b = I_NAI_Ix3y2z_G4z_b+ABY*I_NAI_Hx2y2z_G4z_b;
  Double I_NAI_Hxy3z_Hy4z_b = I_NAI_Ix2y3z_G4z_b+ABY*I_NAI_Hxy3z_G4z_b;
  Double I_NAI_Hx4z_Hy4z_b = I_NAI_Ixy4z_G4z_b+ABY*I_NAI_Hx4z_G4z_b;
  Double I_NAI_H5y_Hy4z_b = I_NAI_I6y_G4z_b+ABY*I_NAI_H5y_G4z_b;
  Double I_NAI_H4yz_Hy4z_b = I_NAI_I5yz_G4z_b+ABY*I_NAI_H4yz_G4z_b;
  Double I_NAI_H3y2z_Hy4z_b = I_NAI_I4y2z_G4z_b+ABY*I_NAI_H3y2z_G4z_b;
  Double I_NAI_H2y3z_Hy4z_b = I_NAI_I3y3z_G4z_b+ABY*I_NAI_H2y3z_G4z_b;
  Double I_NAI_Hy4z_Hy4z_b = I_NAI_I2y4z_G4z_b+ABY*I_NAI_Hy4z_G4z_b;
  Double I_NAI_H5z_Hy4z_b = I_NAI_Iy5z_G4z_b+ABY*I_NAI_H5z_G4z_b;
  Double I_NAI_H5x_H5z_b = I_NAI_I5xz_G4z_b+ABZ*I_NAI_H5x_G4z_b;
  Double I_NAI_H4xy_H5z_b = I_NAI_I4xyz_G4z_b+ABZ*I_NAI_H4xy_G4z_b;
  Double I_NAI_H4xz_H5z_b = I_NAI_I4x2z_G4z_b+ABZ*I_NAI_H4xz_G4z_b;
  Double I_NAI_H3x2y_H5z_b = I_NAI_I3x2yz_G4z_b+ABZ*I_NAI_H3x2y_G4z_b;
  Double I_NAI_H3xyz_H5z_b = I_NAI_I3xy2z_G4z_b+ABZ*I_NAI_H3xyz_G4z_b;
  Double I_NAI_H3x2z_H5z_b = I_NAI_I3x3z_G4z_b+ABZ*I_NAI_H3x2z_G4z_b;
  Double I_NAI_H2x3y_H5z_b = I_NAI_I2x3yz_G4z_b+ABZ*I_NAI_H2x3y_G4z_b;
  Double I_NAI_H2x2yz_H5z_b = I_NAI_I2x2y2z_G4z_b+ABZ*I_NAI_H2x2yz_G4z_b;
  Double I_NAI_H2xy2z_H5z_b = I_NAI_I2xy3z_G4z_b+ABZ*I_NAI_H2xy2z_G4z_b;
  Double I_NAI_H2x3z_H5z_b = I_NAI_I2x4z_G4z_b+ABZ*I_NAI_H2x3z_G4z_b;
  Double I_NAI_Hx4y_H5z_b = I_NAI_Ix4yz_G4z_b+ABZ*I_NAI_Hx4y_G4z_b;
  Double I_NAI_Hx3yz_H5z_b = I_NAI_Ix3y2z_G4z_b+ABZ*I_NAI_Hx3yz_G4z_b;
  Double I_NAI_Hx2y2z_H5z_b = I_NAI_Ix2y3z_G4z_b+ABZ*I_NAI_Hx2y2z_G4z_b;
  Double I_NAI_Hxy3z_H5z_b = I_NAI_Ixy4z_G4z_b+ABZ*I_NAI_Hxy3z_G4z_b;
  Double I_NAI_Hx4z_H5z_b = I_NAI_Ix5z_G4z_b+ABZ*I_NAI_Hx4z_G4z_b;
  Double I_NAI_H5y_H5z_b = I_NAI_I5yz_G4z_b+ABZ*I_NAI_H5y_G4z_b;
  Double I_NAI_H4yz_H5z_b = I_NAI_I4y2z_G4z_b+ABZ*I_NAI_H4yz_G4z_b;
  Double I_NAI_H3y2z_H5z_b = I_NAI_I3y3z_G4z_b+ABZ*I_NAI_H3y2z_G4z_b;
  Double I_NAI_H2y3z_H5z_b = I_NAI_I2y4z_G4z_b+ABZ*I_NAI_H2y3z_G4z_b;
  Double I_NAI_Hy4z_H5z_b = I_NAI_Iy5z_G4z_b+ABZ*I_NAI_Hy4z_G4z_b;
  Double I_NAI_H5z_H5z_b = I_NAI_I6z_G4z_b+ABZ*I_NAI_H5z_G4z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_G_a
   * RHS shell quartet name: SQ_NAI_G_G
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_I6x_G4x_a-5*I_NAI_G4x_G4x;
  abcd[1] = 2.0E0*I_NAI_I5xy_G4x_a-4*I_NAI_G3xy_G4x;
  abcd[2] = 2.0E0*I_NAI_I5xz_G4x_a-4*I_NAI_G3xz_G4x;
  abcd[3] = 2.0E0*I_NAI_I4x2y_G4x_a-3*I_NAI_G2x2y_G4x;
  abcd[4] = 2.0E0*I_NAI_I4xyz_G4x_a-3*I_NAI_G2xyz_G4x;
  abcd[5] = 2.0E0*I_NAI_I4x2z_G4x_a-3*I_NAI_G2x2z_G4x;
  abcd[6] = 2.0E0*I_NAI_I3x3y_G4x_a-2*I_NAI_Gx3y_G4x;
  abcd[7] = 2.0E0*I_NAI_I3x2yz_G4x_a-2*I_NAI_Gx2yz_G4x;
  abcd[8] = 2.0E0*I_NAI_I3xy2z_G4x_a-2*I_NAI_Gxy2z_G4x;
  abcd[9] = 2.0E0*I_NAI_I3x3z_G4x_a-2*I_NAI_Gx3z_G4x;
  abcd[10] = 2.0E0*I_NAI_I2x4y_G4x_a-1*I_NAI_G4y_G4x;
  abcd[11] = 2.0E0*I_NAI_I2x3yz_G4x_a-1*I_NAI_G3yz_G4x;
  abcd[12] = 2.0E0*I_NAI_I2x2y2z_G4x_a-1*I_NAI_G2y2z_G4x;
  abcd[13] = 2.0E0*I_NAI_I2xy3z_G4x_a-1*I_NAI_Gy3z_G4x;
  abcd[14] = 2.0E0*I_NAI_I2x4z_G4x_a-1*I_NAI_G4z_G4x;
  abcd[15] = 2.0E0*I_NAI_Ix5y_G4x_a;
  abcd[16] = 2.0E0*I_NAI_Ix4yz_G4x_a;
  abcd[17] = 2.0E0*I_NAI_Ix3y2z_G4x_a;
  abcd[18] = 2.0E0*I_NAI_Ix2y3z_G4x_a;
  abcd[19] = 2.0E0*I_NAI_Ixy4z_G4x_a;
  abcd[20] = 2.0E0*I_NAI_Ix5z_G4x_a;
  abcd[21] = 2.0E0*I_NAI_I6x_G3xy_a-5*I_NAI_G4x_G3xy;
  abcd[22] = 2.0E0*I_NAI_I5xy_G3xy_a-4*I_NAI_G3xy_G3xy;
  abcd[23] = 2.0E0*I_NAI_I5xz_G3xy_a-4*I_NAI_G3xz_G3xy;
  abcd[24] = 2.0E0*I_NAI_I4x2y_G3xy_a-3*I_NAI_G2x2y_G3xy;
  abcd[25] = 2.0E0*I_NAI_I4xyz_G3xy_a-3*I_NAI_G2xyz_G3xy;
  abcd[26] = 2.0E0*I_NAI_I4x2z_G3xy_a-3*I_NAI_G2x2z_G3xy;
  abcd[27] = 2.0E0*I_NAI_I3x3y_G3xy_a-2*I_NAI_Gx3y_G3xy;
  abcd[28] = 2.0E0*I_NAI_I3x2yz_G3xy_a-2*I_NAI_Gx2yz_G3xy;
  abcd[29] = 2.0E0*I_NAI_I3xy2z_G3xy_a-2*I_NAI_Gxy2z_G3xy;
  abcd[30] = 2.0E0*I_NAI_I3x3z_G3xy_a-2*I_NAI_Gx3z_G3xy;
  abcd[31] = 2.0E0*I_NAI_I2x4y_G3xy_a-1*I_NAI_G4y_G3xy;
  abcd[32] = 2.0E0*I_NAI_I2x3yz_G3xy_a-1*I_NAI_G3yz_G3xy;
  abcd[33] = 2.0E0*I_NAI_I2x2y2z_G3xy_a-1*I_NAI_G2y2z_G3xy;
  abcd[34] = 2.0E0*I_NAI_I2xy3z_G3xy_a-1*I_NAI_Gy3z_G3xy;
  abcd[35] = 2.0E0*I_NAI_I2x4z_G3xy_a-1*I_NAI_G4z_G3xy;
  abcd[36] = 2.0E0*I_NAI_Ix5y_G3xy_a;
  abcd[37] = 2.0E0*I_NAI_Ix4yz_G3xy_a;
  abcd[38] = 2.0E0*I_NAI_Ix3y2z_G3xy_a;
  abcd[39] = 2.0E0*I_NAI_Ix2y3z_G3xy_a;
  abcd[40] = 2.0E0*I_NAI_Ixy4z_G3xy_a;
  abcd[41] = 2.0E0*I_NAI_Ix5z_G3xy_a;
  abcd[42] = 2.0E0*I_NAI_I6x_G3xz_a-5*I_NAI_G4x_G3xz;
  abcd[43] = 2.0E0*I_NAI_I5xy_G3xz_a-4*I_NAI_G3xy_G3xz;
  abcd[44] = 2.0E0*I_NAI_I5xz_G3xz_a-4*I_NAI_G3xz_G3xz;
  abcd[45] = 2.0E0*I_NAI_I4x2y_G3xz_a-3*I_NAI_G2x2y_G3xz;
  abcd[46] = 2.0E0*I_NAI_I4xyz_G3xz_a-3*I_NAI_G2xyz_G3xz;
  abcd[47] = 2.0E0*I_NAI_I4x2z_G3xz_a-3*I_NAI_G2x2z_G3xz;
  abcd[48] = 2.0E0*I_NAI_I3x3y_G3xz_a-2*I_NAI_Gx3y_G3xz;
  abcd[49] = 2.0E0*I_NAI_I3x2yz_G3xz_a-2*I_NAI_Gx2yz_G3xz;
  abcd[50] = 2.0E0*I_NAI_I3xy2z_G3xz_a-2*I_NAI_Gxy2z_G3xz;
  abcd[51] = 2.0E0*I_NAI_I3x3z_G3xz_a-2*I_NAI_Gx3z_G3xz;
  abcd[52] = 2.0E0*I_NAI_I2x4y_G3xz_a-1*I_NAI_G4y_G3xz;
  abcd[53] = 2.0E0*I_NAI_I2x3yz_G3xz_a-1*I_NAI_G3yz_G3xz;
  abcd[54] = 2.0E0*I_NAI_I2x2y2z_G3xz_a-1*I_NAI_G2y2z_G3xz;
  abcd[55] = 2.0E0*I_NAI_I2xy3z_G3xz_a-1*I_NAI_Gy3z_G3xz;
  abcd[56] = 2.0E0*I_NAI_I2x4z_G3xz_a-1*I_NAI_G4z_G3xz;
  abcd[57] = 2.0E0*I_NAI_Ix5y_G3xz_a;
  abcd[58] = 2.0E0*I_NAI_Ix4yz_G3xz_a;
  abcd[59] = 2.0E0*I_NAI_Ix3y2z_G3xz_a;
  abcd[60] = 2.0E0*I_NAI_Ix2y3z_G3xz_a;
  abcd[61] = 2.0E0*I_NAI_Ixy4z_G3xz_a;
  abcd[62] = 2.0E0*I_NAI_Ix5z_G3xz_a;
  abcd[63] = 2.0E0*I_NAI_I6x_G2x2y_a-5*I_NAI_G4x_G2x2y;
  abcd[64] = 2.0E0*I_NAI_I5xy_G2x2y_a-4*I_NAI_G3xy_G2x2y;
  abcd[65] = 2.0E0*I_NAI_I5xz_G2x2y_a-4*I_NAI_G3xz_G2x2y;
  abcd[66] = 2.0E0*I_NAI_I4x2y_G2x2y_a-3*I_NAI_G2x2y_G2x2y;
  abcd[67] = 2.0E0*I_NAI_I4xyz_G2x2y_a-3*I_NAI_G2xyz_G2x2y;
  abcd[68] = 2.0E0*I_NAI_I4x2z_G2x2y_a-3*I_NAI_G2x2z_G2x2y;
  abcd[69] = 2.0E0*I_NAI_I3x3y_G2x2y_a-2*I_NAI_Gx3y_G2x2y;
  abcd[70] = 2.0E0*I_NAI_I3x2yz_G2x2y_a-2*I_NAI_Gx2yz_G2x2y;
  abcd[71] = 2.0E0*I_NAI_I3xy2z_G2x2y_a-2*I_NAI_Gxy2z_G2x2y;
  abcd[72] = 2.0E0*I_NAI_I3x3z_G2x2y_a-2*I_NAI_Gx3z_G2x2y;
  abcd[73] = 2.0E0*I_NAI_I2x4y_G2x2y_a-1*I_NAI_G4y_G2x2y;
  abcd[74] = 2.0E0*I_NAI_I2x3yz_G2x2y_a-1*I_NAI_G3yz_G2x2y;
  abcd[75] = 2.0E0*I_NAI_I2x2y2z_G2x2y_a-1*I_NAI_G2y2z_G2x2y;
  abcd[76] = 2.0E0*I_NAI_I2xy3z_G2x2y_a-1*I_NAI_Gy3z_G2x2y;
  abcd[77] = 2.0E0*I_NAI_I2x4z_G2x2y_a-1*I_NAI_G4z_G2x2y;
  abcd[78] = 2.0E0*I_NAI_Ix5y_G2x2y_a;
  abcd[79] = 2.0E0*I_NAI_Ix4yz_G2x2y_a;
  abcd[80] = 2.0E0*I_NAI_Ix3y2z_G2x2y_a;
  abcd[81] = 2.0E0*I_NAI_Ix2y3z_G2x2y_a;
  abcd[82] = 2.0E0*I_NAI_Ixy4z_G2x2y_a;
  abcd[83] = 2.0E0*I_NAI_Ix5z_G2x2y_a;
  abcd[84] = 2.0E0*I_NAI_I6x_G2xyz_a-5*I_NAI_G4x_G2xyz;
  abcd[85] = 2.0E0*I_NAI_I5xy_G2xyz_a-4*I_NAI_G3xy_G2xyz;
  abcd[86] = 2.0E0*I_NAI_I5xz_G2xyz_a-4*I_NAI_G3xz_G2xyz;
  abcd[87] = 2.0E0*I_NAI_I4x2y_G2xyz_a-3*I_NAI_G2x2y_G2xyz;
  abcd[88] = 2.0E0*I_NAI_I4xyz_G2xyz_a-3*I_NAI_G2xyz_G2xyz;
  abcd[89] = 2.0E0*I_NAI_I4x2z_G2xyz_a-3*I_NAI_G2x2z_G2xyz;
  abcd[90] = 2.0E0*I_NAI_I3x3y_G2xyz_a-2*I_NAI_Gx3y_G2xyz;
  abcd[91] = 2.0E0*I_NAI_I3x2yz_G2xyz_a-2*I_NAI_Gx2yz_G2xyz;
  abcd[92] = 2.0E0*I_NAI_I3xy2z_G2xyz_a-2*I_NAI_Gxy2z_G2xyz;
  abcd[93] = 2.0E0*I_NAI_I3x3z_G2xyz_a-2*I_NAI_Gx3z_G2xyz;
  abcd[94] = 2.0E0*I_NAI_I2x4y_G2xyz_a-1*I_NAI_G4y_G2xyz;
  abcd[95] = 2.0E0*I_NAI_I2x3yz_G2xyz_a-1*I_NAI_G3yz_G2xyz;
  abcd[96] = 2.0E0*I_NAI_I2x2y2z_G2xyz_a-1*I_NAI_G2y2z_G2xyz;
  abcd[97] = 2.0E0*I_NAI_I2xy3z_G2xyz_a-1*I_NAI_Gy3z_G2xyz;
  abcd[98] = 2.0E0*I_NAI_I2x4z_G2xyz_a-1*I_NAI_G4z_G2xyz;
  abcd[99] = 2.0E0*I_NAI_Ix5y_G2xyz_a;
  abcd[100] = 2.0E0*I_NAI_Ix4yz_G2xyz_a;
  abcd[101] = 2.0E0*I_NAI_Ix3y2z_G2xyz_a;
  abcd[102] = 2.0E0*I_NAI_Ix2y3z_G2xyz_a;
  abcd[103] = 2.0E0*I_NAI_Ixy4z_G2xyz_a;
  abcd[104] = 2.0E0*I_NAI_Ix5z_G2xyz_a;
  abcd[105] = 2.0E0*I_NAI_I6x_G2x2z_a-5*I_NAI_G4x_G2x2z;
  abcd[106] = 2.0E0*I_NAI_I5xy_G2x2z_a-4*I_NAI_G3xy_G2x2z;
  abcd[107] = 2.0E0*I_NAI_I5xz_G2x2z_a-4*I_NAI_G3xz_G2x2z;
  abcd[108] = 2.0E0*I_NAI_I4x2y_G2x2z_a-3*I_NAI_G2x2y_G2x2z;
  abcd[109] = 2.0E0*I_NAI_I4xyz_G2x2z_a-3*I_NAI_G2xyz_G2x2z;
  abcd[110] = 2.0E0*I_NAI_I4x2z_G2x2z_a-3*I_NAI_G2x2z_G2x2z;
  abcd[111] = 2.0E0*I_NAI_I3x3y_G2x2z_a-2*I_NAI_Gx3y_G2x2z;
  abcd[112] = 2.0E0*I_NAI_I3x2yz_G2x2z_a-2*I_NAI_Gx2yz_G2x2z;
  abcd[113] = 2.0E0*I_NAI_I3xy2z_G2x2z_a-2*I_NAI_Gxy2z_G2x2z;
  abcd[114] = 2.0E0*I_NAI_I3x3z_G2x2z_a-2*I_NAI_Gx3z_G2x2z;
  abcd[115] = 2.0E0*I_NAI_I2x4y_G2x2z_a-1*I_NAI_G4y_G2x2z;
  abcd[116] = 2.0E0*I_NAI_I2x3yz_G2x2z_a-1*I_NAI_G3yz_G2x2z;
  abcd[117] = 2.0E0*I_NAI_I2x2y2z_G2x2z_a-1*I_NAI_G2y2z_G2x2z;
  abcd[118] = 2.0E0*I_NAI_I2xy3z_G2x2z_a-1*I_NAI_Gy3z_G2x2z;
  abcd[119] = 2.0E0*I_NAI_I2x4z_G2x2z_a-1*I_NAI_G4z_G2x2z;
  abcd[120] = 2.0E0*I_NAI_Ix5y_G2x2z_a;
  abcd[121] = 2.0E0*I_NAI_Ix4yz_G2x2z_a;
  abcd[122] = 2.0E0*I_NAI_Ix3y2z_G2x2z_a;
  abcd[123] = 2.0E0*I_NAI_Ix2y3z_G2x2z_a;
  abcd[124] = 2.0E0*I_NAI_Ixy4z_G2x2z_a;
  abcd[125] = 2.0E0*I_NAI_Ix5z_G2x2z_a;
  abcd[126] = 2.0E0*I_NAI_I6x_Gx3y_a-5*I_NAI_G4x_Gx3y;
  abcd[127] = 2.0E0*I_NAI_I5xy_Gx3y_a-4*I_NAI_G3xy_Gx3y;
  abcd[128] = 2.0E0*I_NAI_I5xz_Gx3y_a-4*I_NAI_G3xz_Gx3y;
  abcd[129] = 2.0E0*I_NAI_I4x2y_Gx3y_a-3*I_NAI_G2x2y_Gx3y;
  abcd[130] = 2.0E0*I_NAI_I4xyz_Gx3y_a-3*I_NAI_G2xyz_Gx3y;
  abcd[131] = 2.0E0*I_NAI_I4x2z_Gx3y_a-3*I_NAI_G2x2z_Gx3y;
  abcd[132] = 2.0E0*I_NAI_I3x3y_Gx3y_a-2*I_NAI_Gx3y_Gx3y;
  abcd[133] = 2.0E0*I_NAI_I3x2yz_Gx3y_a-2*I_NAI_Gx2yz_Gx3y;
  abcd[134] = 2.0E0*I_NAI_I3xy2z_Gx3y_a-2*I_NAI_Gxy2z_Gx3y;
  abcd[135] = 2.0E0*I_NAI_I3x3z_Gx3y_a-2*I_NAI_Gx3z_Gx3y;
  abcd[136] = 2.0E0*I_NAI_I2x4y_Gx3y_a-1*I_NAI_G4y_Gx3y;
  abcd[137] = 2.0E0*I_NAI_I2x3yz_Gx3y_a-1*I_NAI_G3yz_Gx3y;
  abcd[138] = 2.0E0*I_NAI_I2x2y2z_Gx3y_a-1*I_NAI_G2y2z_Gx3y;
  abcd[139] = 2.0E0*I_NAI_I2xy3z_Gx3y_a-1*I_NAI_Gy3z_Gx3y;
  abcd[140] = 2.0E0*I_NAI_I2x4z_Gx3y_a-1*I_NAI_G4z_Gx3y;
  abcd[141] = 2.0E0*I_NAI_Ix5y_Gx3y_a;
  abcd[142] = 2.0E0*I_NAI_Ix4yz_Gx3y_a;
  abcd[143] = 2.0E0*I_NAI_Ix3y2z_Gx3y_a;
  abcd[144] = 2.0E0*I_NAI_Ix2y3z_Gx3y_a;
  abcd[145] = 2.0E0*I_NAI_Ixy4z_Gx3y_a;
  abcd[146] = 2.0E0*I_NAI_Ix5z_Gx3y_a;
  abcd[147] = 2.0E0*I_NAI_I6x_Gx2yz_a-5*I_NAI_G4x_Gx2yz;
  abcd[148] = 2.0E0*I_NAI_I5xy_Gx2yz_a-4*I_NAI_G3xy_Gx2yz;
  abcd[149] = 2.0E0*I_NAI_I5xz_Gx2yz_a-4*I_NAI_G3xz_Gx2yz;
  abcd[150] = 2.0E0*I_NAI_I4x2y_Gx2yz_a-3*I_NAI_G2x2y_Gx2yz;
  abcd[151] = 2.0E0*I_NAI_I4xyz_Gx2yz_a-3*I_NAI_G2xyz_Gx2yz;
  abcd[152] = 2.0E0*I_NAI_I4x2z_Gx2yz_a-3*I_NAI_G2x2z_Gx2yz;
  abcd[153] = 2.0E0*I_NAI_I3x3y_Gx2yz_a-2*I_NAI_Gx3y_Gx2yz;
  abcd[154] = 2.0E0*I_NAI_I3x2yz_Gx2yz_a-2*I_NAI_Gx2yz_Gx2yz;
  abcd[155] = 2.0E0*I_NAI_I3xy2z_Gx2yz_a-2*I_NAI_Gxy2z_Gx2yz;
  abcd[156] = 2.0E0*I_NAI_I3x3z_Gx2yz_a-2*I_NAI_Gx3z_Gx2yz;
  abcd[157] = 2.0E0*I_NAI_I2x4y_Gx2yz_a-1*I_NAI_G4y_Gx2yz;
  abcd[158] = 2.0E0*I_NAI_I2x3yz_Gx2yz_a-1*I_NAI_G3yz_Gx2yz;
  abcd[159] = 2.0E0*I_NAI_I2x2y2z_Gx2yz_a-1*I_NAI_G2y2z_Gx2yz;
  abcd[160] = 2.0E0*I_NAI_I2xy3z_Gx2yz_a-1*I_NAI_Gy3z_Gx2yz;
  abcd[161] = 2.0E0*I_NAI_I2x4z_Gx2yz_a-1*I_NAI_G4z_Gx2yz;
  abcd[162] = 2.0E0*I_NAI_Ix5y_Gx2yz_a;
  abcd[163] = 2.0E0*I_NAI_Ix4yz_Gx2yz_a;
  abcd[164] = 2.0E0*I_NAI_Ix3y2z_Gx2yz_a;
  abcd[165] = 2.0E0*I_NAI_Ix2y3z_Gx2yz_a;
  abcd[166] = 2.0E0*I_NAI_Ixy4z_Gx2yz_a;
  abcd[167] = 2.0E0*I_NAI_Ix5z_Gx2yz_a;
  abcd[168] = 2.0E0*I_NAI_I6x_Gxy2z_a-5*I_NAI_G4x_Gxy2z;
  abcd[169] = 2.0E0*I_NAI_I5xy_Gxy2z_a-4*I_NAI_G3xy_Gxy2z;
  abcd[170] = 2.0E0*I_NAI_I5xz_Gxy2z_a-4*I_NAI_G3xz_Gxy2z;
  abcd[171] = 2.0E0*I_NAI_I4x2y_Gxy2z_a-3*I_NAI_G2x2y_Gxy2z;
  abcd[172] = 2.0E0*I_NAI_I4xyz_Gxy2z_a-3*I_NAI_G2xyz_Gxy2z;
  abcd[173] = 2.0E0*I_NAI_I4x2z_Gxy2z_a-3*I_NAI_G2x2z_Gxy2z;
  abcd[174] = 2.0E0*I_NAI_I3x3y_Gxy2z_a-2*I_NAI_Gx3y_Gxy2z;
  abcd[175] = 2.0E0*I_NAI_I3x2yz_Gxy2z_a-2*I_NAI_Gx2yz_Gxy2z;
  abcd[176] = 2.0E0*I_NAI_I3xy2z_Gxy2z_a-2*I_NAI_Gxy2z_Gxy2z;
  abcd[177] = 2.0E0*I_NAI_I3x3z_Gxy2z_a-2*I_NAI_Gx3z_Gxy2z;
  abcd[178] = 2.0E0*I_NAI_I2x4y_Gxy2z_a-1*I_NAI_G4y_Gxy2z;
  abcd[179] = 2.0E0*I_NAI_I2x3yz_Gxy2z_a-1*I_NAI_G3yz_Gxy2z;
  abcd[180] = 2.0E0*I_NAI_I2x2y2z_Gxy2z_a-1*I_NAI_G2y2z_Gxy2z;
  abcd[181] = 2.0E0*I_NAI_I2xy3z_Gxy2z_a-1*I_NAI_Gy3z_Gxy2z;
  abcd[182] = 2.0E0*I_NAI_I2x4z_Gxy2z_a-1*I_NAI_G4z_Gxy2z;
  abcd[183] = 2.0E0*I_NAI_Ix5y_Gxy2z_a;
  abcd[184] = 2.0E0*I_NAI_Ix4yz_Gxy2z_a;
  abcd[185] = 2.0E0*I_NAI_Ix3y2z_Gxy2z_a;
  abcd[186] = 2.0E0*I_NAI_Ix2y3z_Gxy2z_a;
  abcd[187] = 2.0E0*I_NAI_Ixy4z_Gxy2z_a;
  abcd[188] = 2.0E0*I_NAI_Ix5z_Gxy2z_a;
  abcd[189] = 2.0E0*I_NAI_I6x_Gx3z_a-5*I_NAI_G4x_Gx3z;
  abcd[190] = 2.0E0*I_NAI_I5xy_Gx3z_a-4*I_NAI_G3xy_Gx3z;
  abcd[191] = 2.0E0*I_NAI_I5xz_Gx3z_a-4*I_NAI_G3xz_Gx3z;
  abcd[192] = 2.0E0*I_NAI_I4x2y_Gx3z_a-3*I_NAI_G2x2y_Gx3z;
  abcd[193] = 2.0E0*I_NAI_I4xyz_Gx3z_a-3*I_NAI_G2xyz_Gx3z;
  abcd[194] = 2.0E0*I_NAI_I4x2z_Gx3z_a-3*I_NAI_G2x2z_Gx3z;
  abcd[195] = 2.0E0*I_NAI_I3x3y_Gx3z_a-2*I_NAI_Gx3y_Gx3z;
  abcd[196] = 2.0E0*I_NAI_I3x2yz_Gx3z_a-2*I_NAI_Gx2yz_Gx3z;
  abcd[197] = 2.0E0*I_NAI_I3xy2z_Gx3z_a-2*I_NAI_Gxy2z_Gx3z;
  abcd[198] = 2.0E0*I_NAI_I3x3z_Gx3z_a-2*I_NAI_Gx3z_Gx3z;
  abcd[199] = 2.0E0*I_NAI_I2x4y_Gx3z_a-1*I_NAI_G4y_Gx3z;
  abcd[200] = 2.0E0*I_NAI_I2x3yz_Gx3z_a-1*I_NAI_G3yz_Gx3z;
  abcd[201] = 2.0E0*I_NAI_I2x2y2z_Gx3z_a-1*I_NAI_G2y2z_Gx3z;
  abcd[202] = 2.0E0*I_NAI_I2xy3z_Gx3z_a-1*I_NAI_Gy3z_Gx3z;
  abcd[203] = 2.0E0*I_NAI_I2x4z_Gx3z_a-1*I_NAI_G4z_Gx3z;
  abcd[204] = 2.0E0*I_NAI_Ix5y_Gx3z_a;
  abcd[205] = 2.0E0*I_NAI_Ix4yz_Gx3z_a;
  abcd[206] = 2.0E0*I_NAI_Ix3y2z_Gx3z_a;
  abcd[207] = 2.0E0*I_NAI_Ix2y3z_Gx3z_a;
  abcd[208] = 2.0E0*I_NAI_Ixy4z_Gx3z_a;
  abcd[209] = 2.0E0*I_NAI_Ix5z_Gx3z_a;
  abcd[210] = 2.0E0*I_NAI_I6x_G4y_a-5*I_NAI_G4x_G4y;
  abcd[211] = 2.0E0*I_NAI_I5xy_G4y_a-4*I_NAI_G3xy_G4y;
  abcd[212] = 2.0E0*I_NAI_I5xz_G4y_a-4*I_NAI_G3xz_G4y;
  abcd[213] = 2.0E0*I_NAI_I4x2y_G4y_a-3*I_NAI_G2x2y_G4y;
  abcd[214] = 2.0E0*I_NAI_I4xyz_G4y_a-3*I_NAI_G2xyz_G4y;
  abcd[215] = 2.0E0*I_NAI_I4x2z_G4y_a-3*I_NAI_G2x2z_G4y;
  abcd[216] = 2.0E0*I_NAI_I3x3y_G4y_a-2*I_NAI_Gx3y_G4y;
  abcd[217] = 2.0E0*I_NAI_I3x2yz_G4y_a-2*I_NAI_Gx2yz_G4y;
  abcd[218] = 2.0E0*I_NAI_I3xy2z_G4y_a-2*I_NAI_Gxy2z_G4y;
  abcd[219] = 2.0E0*I_NAI_I3x3z_G4y_a-2*I_NAI_Gx3z_G4y;
  abcd[220] = 2.0E0*I_NAI_I2x4y_G4y_a-1*I_NAI_G4y_G4y;
  abcd[221] = 2.0E0*I_NAI_I2x3yz_G4y_a-1*I_NAI_G3yz_G4y;
  abcd[222] = 2.0E0*I_NAI_I2x2y2z_G4y_a-1*I_NAI_G2y2z_G4y;
  abcd[223] = 2.0E0*I_NAI_I2xy3z_G4y_a-1*I_NAI_Gy3z_G4y;
  abcd[224] = 2.0E0*I_NAI_I2x4z_G4y_a-1*I_NAI_G4z_G4y;
  abcd[225] = 2.0E0*I_NAI_Ix5y_G4y_a;
  abcd[226] = 2.0E0*I_NAI_Ix4yz_G4y_a;
  abcd[227] = 2.0E0*I_NAI_Ix3y2z_G4y_a;
  abcd[228] = 2.0E0*I_NAI_Ix2y3z_G4y_a;
  abcd[229] = 2.0E0*I_NAI_Ixy4z_G4y_a;
  abcd[230] = 2.0E0*I_NAI_Ix5z_G4y_a;
  abcd[231] = 2.0E0*I_NAI_I6x_G3yz_a-5*I_NAI_G4x_G3yz;
  abcd[232] = 2.0E0*I_NAI_I5xy_G3yz_a-4*I_NAI_G3xy_G3yz;
  abcd[233] = 2.0E0*I_NAI_I5xz_G3yz_a-4*I_NAI_G3xz_G3yz;
  abcd[234] = 2.0E0*I_NAI_I4x2y_G3yz_a-3*I_NAI_G2x2y_G3yz;
  abcd[235] = 2.0E0*I_NAI_I4xyz_G3yz_a-3*I_NAI_G2xyz_G3yz;
  abcd[236] = 2.0E0*I_NAI_I4x2z_G3yz_a-3*I_NAI_G2x2z_G3yz;
  abcd[237] = 2.0E0*I_NAI_I3x3y_G3yz_a-2*I_NAI_Gx3y_G3yz;
  abcd[238] = 2.0E0*I_NAI_I3x2yz_G3yz_a-2*I_NAI_Gx2yz_G3yz;
  abcd[239] = 2.0E0*I_NAI_I3xy2z_G3yz_a-2*I_NAI_Gxy2z_G3yz;
  abcd[240] = 2.0E0*I_NAI_I3x3z_G3yz_a-2*I_NAI_Gx3z_G3yz;
  abcd[241] = 2.0E0*I_NAI_I2x4y_G3yz_a-1*I_NAI_G4y_G3yz;
  abcd[242] = 2.0E0*I_NAI_I2x3yz_G3yz_a-1*I_NAI_G3yz_G3yz;
  abcd[243] = 2.0E0*I_NAI_I2x2y2z_G3yz_a-1*I_NAI_G2y2z_G3yz;
  abcd[244] = 2.0E0*I_NAI_I2xy3z_G3yz_a-1*I_NAI_Gy3z_G3yz;
  abcd[245] = 2.0E0*I_NAI_I2x4z_G3yz_a-1*I_NAI_G4z_G3yz;
  abcd[246] = 2.0E0*I_NAI_Ix5y_G3yz_a;
  abcd[247] = 2.0E0*I_NAI_Ix4yz_G3yz_a;
  abcd[248] = 2.0E0*I_NAI_Ix3y2z_G3yz_a;
  abcd[249] = 2.0E0*I_NAI_Ix2y3z_G3yz_a;
  abcd[250] = 2.0E0*I_NAI_Ixy4z_G3yz_a;
  abcd[251] = 2.0E0*I_NAI_Ix5z_G3yz_a;
  abcd[252] = 2.0E0*I_NAI_I6x_G2y2z_a-5*I_NAI_G4x_G2y2z;
  abcd[253] = 2.0E0*I_NAI_I5xy_G2y2z_a-4*I_NAI_G3xy_G2y2z;
  abcd[254] = 2.0E0*I_NAI_I5xz_G2y2z_a-4*I_NAI_G3xz_G2y2z;
  abcd[255] = 2.0E0*I_NAI_I4x2y_G2y2z_a-3*I_NAI_G2x2y_G2y2z;
  abcd[256] = 2.0E0*I_NAI_I4xyz_G2y2z_a-3*I_NAI_G2xyz_G2y2z;
  abcd[257] = 2.0E0*I_NAI_I4x2z_G2y2z_a-3*I_NAI_G2x2z_G2y2z;
  abcd[258] = 2.0E0*I_NAI_I3x3y_G2y2z_a-2*I_NAI_Gx3y_G2y2z;
  abcd[259] = 2.0E0*I_NAI_I3x2yz_G2y2z_a-2*I_NAI_Gx2yz_G2y2z;
  abcd[260] = 2.0E0*I_NAI_I3xy2z_G2y2z_a-2*I_NAI_Gxy2z_G2y2z;
  abcd[261] = 2.0E0*I_NAI_I3x3z_G2y2z_a-2*I_NAI_Gx3z_G2y2z;
  abcd[262] = 2.0E0*I_NAI_I2x4y_G2y2z_a-1*I_NAI_G4y_G2y2z;
  abcd[263] = 2.0E0*I_NAI_I2x3yz_G2y2z_a-1*I_NAI_G3yz_G2y2z;
  abcd[264] = 2.0E0*I_NAI_I2x2y2z_G2y2z_a-1*I_NAI_G2y2z_G2y2z;
  abcd[265] = 2.0E0*I_NAI_I2xy3z_G2y2z_a-1*I_NAI_Gy3z_G2y2z;
  abcd[266] = 2.0E0*I_NAI_I2x4z_G2y2z_a-1*I_NAI_G4z_G2y2z;
  abcd[267] = 2.0E0*I_NAI_Ix5y_G2y2z_a;
  abcd[268] = 2.0E0*I_NAI_Ix4yz_G2y2z_a;
  abcd[269] = 2.0E0*I_NAI_Ix3y2z_G2y2z_a;
  abcd[270] = 2.0E0*I_NAI_Ix2y3z_G2y2z_a;
  abcd[271] = 2.0E0*I_NAI_Ixy4z_G2y2z_a;
  abcd[272] = 2.0E0*I_NAI_Ix5z_G2y2z_a;
  abcd[273] = 2.0E0*I_NAI_I6x_Gy3z_a-5*I_NAI_G4x_Gy3z;
  abcd[274] = 2.0E0*I_NAI_I5xy_Gy3z_a-4*I_NAI_G3xy_Gy3z;
  abcd[275] = 2.0E0*I_NAI_I5xz_Gy3z_a-4*I_NAI_G3xz_Gy3z;
  abcd[276] = 2.0E0*I_NAI_I4x2y_Gy3z_a-3*I_NAI_G2x2y_Gy3z;
  abcd[277] = 2.0E0*I_NAI_I4xyz_Gy3z_a-3*I_NAI_G2xyz_Gy3z;
  abcd[278] = 2.0E0*I_NAI_I4x2z_Gy3z_a-3*I_NAI_G2x2z_Gy3z;
  abcd[279] = 2.0E0*I_NAI_I3x3y_Gy3z_a-2*I_NAI_Gx3y_Gy3z;
  abcd[280] = 2.0E0*I_NAI_I3x2yz_Gy3z_a-2*I_NAI_Gx2yz_Gy3z;
  abcd[281] = 2.0E0*I_NAI_I3xy2z_Gy3z_a-2*I_NAI_Gxy2z_Gy3z;
  abcd[282] = 2.0E0*I_NAI_I3x3z_Gy3z_a-2*I_NAI_Gx3z_Gy3z;
  abcd[283] = 2.0E0*I_NAI_I2x4y_Gy3z_a-1*I_NAI_G4y_Gy3z;
  abcd[284] = 2.0E0*I_NAI_I2x3yz_Gy3z_a-1*I_NAI_G3yz_Gy3z;
  abcd[285] = 2.0E0*I_NAI_I2x2y2z_Gy3z_a-1*I_NAI_G2y2z_Gy3z;
  abcd[286] = 2.0E0*I_NAI_I2xy3z_Gy3z_a-1*I_NAI_Gy3z_Gy3z;
  abcd[287] = 2.0E0*I_NAI_I2x4z_Gy3z_a-1*I_NAI_G4z_Gy3z;
  abcd[288] = 2.0E0*I_NAI_Ix5y_Gy3z_a;
  abcd[289] = 2.0E0*I_NAI_Ix4yz_Gy3z_a;
  abcd[290] = 2.0E0*I_NAI_Ix3y2z_Gy3z_a;
  abcd[291] = 2.0E0*I_NAI_Ix2y3z_Gy3z_a;
  abcd[292] = 2.0E0*I_NAI_Ixy4z_Gy3z_a;
  abcd[293] = 2.0E0*I_NAI_Ix5z_Gy3z_a;
  abcd[294] = 2.0E0*I_NAI_I6x_G4z_a-5*I_NAI_G4x_G4z;
  abcd[295] = 2.0E0*I_NAI_I5xy_G4z_a-4*I_NAI_G3xy_G4z;
  abcd[296] = 2.0E0*I_NAI_I5xz_G4z_a-4*I_NAI_G3xz_G4z;
  abcd[297] = 2.0E0*I_NAI_I4x2y_G4z_a-3*I_NAI_G2x2y_G4z;
  abcd[298] = 2.0E0*I_NAI_I4xyz_G4z_a-3*I_NAI_G2xyz_G4z;
  abcd[299] = 2.0E0*I_NAI_I4x2z_G4z_a-3*I_NAI_G2x2z_G4z;
  abcd[300] = 2.0E0*I_NAI_I3x3y_G4z_a-2*I_NAI_Gx3y_G4z;
  abcd[301] = 2.0E0*I_NAI_I3x2yz_G4z_a-2*I_NAI_Gx2yz_G4z;
  abcd[302] = 2.0E0*I_NAI_I3xy2z_G4z_a-2*I_NAI_Gxy2z_G4z;
  abcd[303] = 2.0E0*I_NAI_I3x3z_G4z_a-2*I_NAI_Gx3z_G4z;
  abcd[304] = 2.0E0*I_NAI_I2x4y_G4z_a-1*I_NAI_G4y_G4z;
  abcd[305] = 2.0E0*I_NAI_I2x3yz_G4z_a-1*I_NAI_G3yz_G4z;
  abcd[306] = 2.0E0*I_NAI_I2x2y2z_G4z_a-1*I_NAI_G2y2z_G4z;
  abcd[307] = 2.0E0*I_NAI_I2xy3z_G4z_a-1*I_NAI_Gy3z_G4z;
  abcd[308] = 2.0E0*I_NAI_I2x4z_G4z_a-1*I_NAI_G4z_G4z;
  abcd[309] = 2.0E0*I_NAI_Ix5y_G4z_a;
  abcd[310] = 2.0E0*I_NAI_Ix4yz_G4z_a;
  abcd[311] = 2.0E0*I_NAI_Ix3y2z_G4z_a;
  abcd[312] = 2.0E0*I_NAI_Ix2y3z_G4z_a;
  abcd[313] = 2.0E0*I_NAI_Ixy4z_G4z_a;
  abcd[314] = 2.0E0*I_NAI_Ix5z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_G_a
   * RHS shell quartet name: SQ_NAI_G_G
   ************************************************************/
  abcd[315] = 2.0E0*I_NAI_I5xy_G4x_a;
  abcd[316] = 2.0E0*I_NAI_I4x2y_G4x_a-1*I_NAI_G4x_G4x;
  abcd[317] = 2.0E0*I_NAI_I4xyz_G4x_a;
  abcd[318] = 2.0E0*I_NAI_I3x3y_G4x_a-2*I_NAI_G3xy_G4x;
  abcd[319] = 2.0E0*I_NAI_I3x2yz_G4x_a-1*I_NAI_G3xz_G4x;
  abcd[320] = 2.0E0*I_NAI_I3xy2z_G4x_a;
  abcd[321] = 2.0E0*I_NAI_I2x4y_G4x_a-3*I_NAI_G2x2y_G4x;
  abcd[322] = 2.0E0*I_NAI_I2x3yz_G4x_a-2*I_NAI_G2xyz_G4x;
  abcd[323] = 2.0E0*I_NAI_I2x2y2z_G4x_a-1*I_NAI_G2x2z_G4x;
  abcd[324] = 2.0E0*I_NAI_I2xy3z_G4x_a;
  abcd[325] = 2.0E0*I_NAI_Ix5y_G4x_a-4*I_NAI_Gx3y_G4x;
  abcd[326] = 2.0E0*I_NAI_Ix4yz_G4x_a-3*I_NAI_Gx2yz_G4x;
  abcd[327] = 2.0E0*I_NAI_Ix3y2z_G4x_a-2*I_NAI_Gxy2z_G4x;
  abcd[328] = 2.0E0*I_NAI_Ix2y3z_G4x_a-1*I_NAI_Gx3z_G4x;
  abcd[329] = 2.0E0*I_NAI_Ixy4z_G4x_a;
  abcd[330] = 2.0E0*I_NAI_I6y_G4x_a-5*I_NAI_G4y_G4x;
  abcd[331] = 2.0E0*I_NAI_I5yz_G4x_a-4*I_NAI_G3yz_G4x;
  abcd[332] = 2.0E0*I_NAI_I4y2z_G4x_a-3*I_NAI_G2y2z_G4x;
  abcd[333] = 2.0E0*I_NAI_I3y3z_G4x_a-2*I_NAI_Gy3z_G4x;
  abcd[334] = 2.0E0*I_NAI_I2y4z_G4x_a-1*I_NAI_G4z_G4x;
  abcd[335] = 2.0E0*I_NAI_Iy5z_G4x_a;
  abcd[336] = 2.0E0*I_NAI_I5xy_G3xy_a;
  abcd[337] = 2.0E0*I_NAI_I4x2y_G3xy_a-1*I_NAI_G4x_G3xy;
  abcd[338] = 2.0E0*I_NAI_I4xyz_G3xy_a;
  abcd[339] = 2.0E0*I_NAI_I3x3y_G3xy_a-2*I_NAI_G3xy_G3xy;
  abcd[340] = 2.0E0*I_NAI_I3x2yz_G3xy_a-1*I_NAI_G3xz_G3xy;
  abcd[341] = 2.0E0*I_NAI_I3xy2z_G3xy_a;
  abcd[342] = 2.0E0*I_NAI_I2x4y_G3xy_a-3*I_NAI_G2x2y_G3xy;
  abcd[343] = 2.0E0*I_NAI_I2x3yz_G3xy_a-2*I_NAI_G2xyz_G3xy;
  abcd[344] = 2.0E0*I_NAI_I2x2y2z_G3xy_a-1*I_NAI_G2x2z_G3xy;
  abcd[345] = 2.0E0*I_NAI_I2xy3z_G3xy_a;
  abcd[346] = 2.0E0*I_NAI_Ix5y_G3xy_a-4*I_NAI_Gx3y_G3xy;
  abcd[347] = 2.0E0*I_NAI_Ix4yz_G3xy_a-3*I_NAI_Gx2yz_G3xy;
  abcd[348] = 2.0E0*I_NAI_Ix3y2z_G3xy_a-2*I_NAI_Gxy2z_G3xy;
  abcd[349] = 2.0E0*I_NAI_Ix2y3z_G3xy_a-1*I_NAI_Gx3z_G3xy;
  abcd[350] = 2.0E0*I_NAI_Ixy4z_G3xy_a;
  abcd[351] = 2.0E0*I_NAI_I6y_G3xy_a-5*I_NAI_G4y_G3xy;
  abcd[352] = 2.0E0*I_NAI_I5yz_G3xy_a-4*I_NAI_G3yz_G3xy;
  abcd[353] = 2.0E0*I_NAI_I4y2z_G3xy_a-3*I_NAI_G2y2z_G3xy;
  abcd[354] = 2.0E0*I_NAI_I3y3z_G3xy_a-2*I_NAI_Gy3z_G3xy;
  abcd[355] = 2.0E0*I_NAI_I2y4z_G3xy_a-1*I_NAI_G4z_G3xy;
  abcd[356] = 2.0E0*I_NAI_Iy5z_G3xy_a;
  abcd[357] = 2.0E0*I_NAI_I5xy_G3xz_a;
  abcd[358] = 2.0E0*I_NAI_I4x2y_G3xz_a-1*I_NAI_G4x_G3xz;
  abcd[359] = 2.0E0*I_NAI_I4xyz_G3xz_a;
  abcd[360] = 2.0E0*I_NAI_I3x3y_G3xz_a-2*I_NAI_G3xy_G3xz;
  abcd[361] = 2.0E0*I_NAI_I3x2yz_G3xz_a-1*I_NAI_G3xz_G3xz;
  abcd[362] = 2.0E0*I_NAI_I3xy2z_G3xz_a;
  abcd[363] = 2.0E0*I_NAI_I2x4y_G3xz_a-3*I_NAI_G2x2y_G3xz;
  abcd[364] = 2.0E0*I_NAI_I2x3yz_G3xz_a-2*I_NAI_G2xyz_G3xz;
  abcd[365] = 2.0E0*I_NAI_I2x2y2z_G3xz_a-1*I_NAI_G2x2z_G3xz;
  abcd[366] = 2.0E0*I_NAI_I2xy3z_G3xz_a;
  abcd[367] = 2.0E0*I_NAI_Ix5y_G3xz_a-4*I_NAI_Gx3y_G3xz;
  abcd[368] = 2.0E0*I_NAI_Ix4yz_G3xz_a-3*I_NAI_Gx2yz_G3xz;
  abcd[369] = 2.0E0*I_NAI_Ix3y2z_G3xz_a-2*I_NAI_Gxy2z_G3xz;
  abcd[370] = 2.0E0*I_NAI_Ix2y3z_G3xz_a-1*I_NAI_Gx3z_G3xz;
  abcd[371] = 2.0E0*I_NAI_Ixy4z_G3xz_a;
  abcd[372] = 2.0E0*I_NAI_I6y_G3xz_a-5*I_NAI_G4y_G3xz;
  abcd[373] = 2.0E0*I_NAI_I5yz_G3xz_a-4*I_NAI_G3yz_G3xz;
  abcd[374] = 2.0E0*I_NAI_I4y2z_G3xz_a-3*I_NAI_G2y2z_G3xz;
  abcd[375] = 2.0E0*I_NAI_I3y3z_G3xz_a-2*I_NAI_Gy3z_G3xz;
  abcd[376] = 2.0E0*I_NAI_I2y4z_G3xz_a-1*I_NAI_G4z_G3xz;
  abcd[377] = 2.0E0*I_NAI_Iy5z_G3xz_a;
  abcd[378] = 2.0E0*I_NAI_I5xy_G2x2y_a;
  abcd[379] = 2.0E0*I_NAI_I4x2y_G2x2y_a-1*I_NAI_G4x_G2x2y;
  abcd[380] = 2.0E0*I_NAI_I4xyz_G2x2y_a;
  abcd[381] = 2.0E0*I_NAI_I3x3y_G2x2y_a-2*I_NAI_G3xy_G2x2y;
  abcd[382] = 2.0E0*I_NAI_I3x2yz_G2x2y_a-1*I_NAI_G3xz_G2x2y;
  abcd[383] = 2.0E0*I_NAI_I3xy2z_G2x2y_a;
  abcd[384] = 2.0E0*I_NAI_I2x4y_G2x2y_a-3*I_NAI_G2x2y_G2x2y;
  abcd[385] = 2.0E0*I_NAI_I2x3yz_G2x2y_a-2*I_NAI_G2xyz_G2x2y;
  abcd[386] = 2.0E0*I_NAI_I2x2y2z_G2x2y_a-1*I_NAI_G2x2z_G2x2y;
  abcd[387] = 2.0E0*I_NAI_I2xy3z_G2x2y_a;
  abcd[388] = 2.0E0*I_NAI_Ix5y_G2x2y_a-4*I_NAI_Gx3y_G2x2y;
  abcd[389] = 2.0E0*I_NAI_Ix4yz_G2x2y_a-3*I_NAI_Gx2yz_G2x2y;
  abcd[390] = 2.0E0*I_NAI_Ix3y2z_G2x2y_a-2*I_NAI_Gxy2z_G2x2y;
  abcd[391] = 2.0E0*I_NAI_Ix2y3z_G2x2y_a-1*I_NAI_Gx3z_G2x2y;
  abcd[392] = 2.0E0*I_NAI_Ixy4z_G2x2y_a;
  abcd[393] = 2.0E0*I_NAI_I6y_G2x2y_a-5*I_NAI_G4y_G2x2y;
  abcd[394] = 2.0E0*I_NAI_I5yz_G2x2y_a-4*I_NAI_G3yz_G2x2y;
  abcd[395] = 2.0E0*I_NAI_I4y2z_G2x2y_a-3*I_NAI_G2y2z_G2x2y;
  abcd[396] = 2.0E0*I_NAI_I3y3z_G2x2y_a-2*I_NAI_Gy3z_G2x2y;
  abcd[397] = 2.0E0*I_NAI_I2y4z_G2x2y_a-1*I_NAI_G4z_G2x2y;
  abcd[398] = 2.0E0*I_NAI_Iy5z_G2x2y_a;
  abcd[399] = 2.0E0*I_NAI_I5xy_G2xyz_a;
  abcd[400] = 2.0E0*I_NAI_I4x2y_G2xyz_a-1*I_NAI_G4x_G2xyz;
  abcd[401] = 2.0E0*I_NAI_I4xyz_G2xyz_a;
  abcd[402] = 2.0E0*I_NAI_I3x3y_G2xyz_a-2*I_NAI_G3xy_G2xyz;
  abcd[403] = 2.0E0*I_NAI_I3x2yz_G2xyz_a-1*I_NAI_G3xz_G2xyz;
  abcd[404] = 2.0E0*I_NAI_I3xy2z_G2xyz_a;
  abcd[405] = 2.0E0*I_NAI_I2x4y_G2xyz_a-3*I_NAI_G2x2y_G2xyz;
  abcd[406] = 2.0E0*I_NAI_I2x3yz_G2xyz_a-2*I_NAI_G2xyz_G2xyz;
  abcd[407] = 2.0E0*I_NAI_I2x2y2z_G2xyz_a-1*I_NAI_G2x2z_G2xyz;
  abcd[408] = 2.0E0*I_NAI_I2xy3z_G2xyz_a;
  abcd[409] = 2.0E0*I_NAI_Ix5y_G2xyz_a-4*I_NAI_Gx3y_G2xyz;
  abcd[410] = 2.0E0*I_NAI_Ix4yz_G2xyz_a-3*I_NAI_Gx2yz_G2xyz;
  abcd[411] = 2.0E0*I_NAI_Ix3y2z_G2xyz_a-2*I_NAI_Gxy2z_G2xyz;
  abcd[412] = 2.0E0*I_NAI_Ix2y3z_G2xyz_a-1*I_NAI_Gx3z_G2xyz;
  abcd[413] = 2.0E0*I_NAI_Ixy4z_G2xyz_a;
  abcd[414] = 2.0E0*I_NAI_I6y_G2xyz_a-5*I_NAI_G4y_G2xyz;
  abcd[415] = 2.0E0*I_NAI_I5yz_G2xyz_a-4*I_NAI_G3yz_G2xyz;
  abcd[416] = 2.0E0*I_NAI_I4y2z_G2xyz_a-3*I_NAI_G2y2z_G2xyz;
  abcd[417] = 2.0E0*I_NAI_I3y3z_G2xyz_a-2*I_NAI_Gy3z_G2xyz;
  abcd[418] = 2.0E0*I_NAI_I2y4z_G2xyz_a-1*I_NAI_G4z_G2xyz;
  abcd[419] = 2.0E0*I_NAI_Iy5z_G2xyz_a;
  abcd[420] = 2.0E0*I_NAI_I5xy_G2x2z_a;
  abcd[421] = 2.0E0*I_NAI_I4x2y_G2x2z_a-1*I_NAI_G4x_G2x2z;
  abcd[422] = 2.0E0*I_NAI_I4xyz_G2x2z_a;
  abcd[423] = 2.0E0*I_NAI_I3x3y_G2x2z_a-2*I_NAI_G3xy_G2x2z;
  abcd[424] = 2.0E0*I_NAI_I3x2yz_G2x2z_a-1*I_NAI_G3xz_G2x2z;
  abcd[425] = 2.0E0*I_NAI_I3xy2z_G2x2z_a;
  abcd[426] = 2.0E0*I_NAI_I2x4y_G2x2z_a-3*I_NAI_G2x2y_G2x2z;
  abcd[427] = 2.0E0*I_NAI_I2x3yz_G2x2z_a-2*I_NAI_G2xyz_G2x2z;
  abcd[428] = 2.0E0*I_NAI_I2x2y2z_G2x2z_a-1*I_NAI_G2x2z_G2x2z;
  abcd[429] = 2.0E0*I_NAI_I2xy3z_G2x2z_a;
  abcd[430] = 2.0E0*I_NAI_Ix5y_G2x2z_a-4*I_NAI_Gx3y_G2x2z;
  abcd[431] = 2.0E0*I_NAI_Ix4yz_G2x2z_a-3*I_NAI_Gx2yz_G2x2z;
  abcd[432] = 2.0E0*I_NAI_Ix3y2z_G2x2z_a-2*I_NAI_Gxy2z_G2x2z;
  abcd[433] = 2.0E0*I_NAI_Ix2y3z_G2x2z_a-1*I_NAI_Gx3z_G2x2z;
  abcd[434] = 2.0E0*I_NAI_Ixy4z_G2x2z_a;
  abcd[435] = 2.0E0*I_NAI_I6y_G2x2z_a-5*I_NAI_G4y_G2x2z;
  abcd[436] = 2.0E0*I_NAI_I5yz_G2x2z_a-4*I_NAI_G3yz_G2x2z;
  abcd[437] = 2.0E0*I_NAI_I4y2z_G2x2z_a-3*I_NAI_G2y2z_G2x2z;
  abcd[438] = 2.0E0*I_NAI_I3y3z_G2x2z_a-2*I_NAI_Gy3z_G2x2z;
  abcd[439] = 2.0E0*I_NAI_I2y4z_G2x2z_a-1*I_NAI_G4z_G2x2z;
  abcd[440] = 2.0E0*I_NAI_Iy5z_G2x2z_a;
  abcd[441] = 2.0E0*I_NAI_I5xy_Gx3y_a;
  abcd[442] = 2.0E0*I_NAI_I4x2y_Gx3y_a-1*I_NAI_G4x_Gx3y;
  abcd[443] = 2.0E0*I_NAI_I4xyz_Gx3y_a;
  abcd[444] = 2.0E0*I_NAI_I3x3y_Gx3y_a-2*I_NAI_G3xy_Gx3y;
  abcd[445] = 2.0E0*I_NAI_I3x2yz_Gx3y_a-1*I_NAI_G3xz_Gx3y;
  abcd[446] = 2.0E0*I_NAI_I3xy2z_Gx3y_a;
  abcd[447] = 2.0E0*I_NAI_I2x4y_Gx3y_a-3*I_NAI_G2x2y_Gx3y;
  abcd[448] = 2.0E0*I_NAI_I2x3yz_Gx3y_a-2*I_NAI_G2xyz_Gx3y;
  abcd[449] = 2.0E0*I_NAI_I2x2y2z_Gx3y_a-1*I_NAI_G2x2z_Gx3y;
  abcd[450] = 2.0E0*I_NAI_I2xy3z_Gx3y_a;
  abcd[451] = 2.0E0*I_NAI_Ix5y_Gx3y_a-4*I_NAI_Gx3y_Gx3y;
  abcd[452] = 2.0E0*I_NAI_Ix4yz_Gx3y_a-3*I_NAI_Gx2yz_Gx3y;
  abcd[453] = 2.0E0*I_NAI_Ix3y2z_Gx3y_a-2*I_NAI_Gxy2z_Gx3y;
  abcd[454] = 2.0E0*I_NAI_Ix2y3z_Gx3y_a-1*I_NAI_Gx3z_Gx3y;
  abcd[455] = 2.0E0*I_NAI_Ixy4z_Gx3y_a;
  abcd[456] = 2.0E0*I_NAI_I6y_Gx3y_a-5*I_NAI_G4y_Gx3y;
  abcd[457] = 2.0E0*I_NAI_I5yz_Gx3y_a-4*I_NAI_G3yz_Gx3y;
  abcd[458] = 2.0E0*I_NAI_I4y2z_Gx3y_a-3*I_NAI_G2y2z_Gx3y;
  abcd[459] = 2.0E0*I_NAI_I3y3z_Gx3y_a-2*I_NAI_Gy3z_Gx3y;
  abcd[460] = 2.0E0*I_NAI_I2y4z_Gx3y_a-1*I_NAI_G4z_Gx3y;
  abcd[461] = 2.0E0*I_NAI_Iy5z_Gx3y_a;
  abcd[462] = 2.0E0*I_NAI_I5xy_Gx2yz_a;
  abcd[463] = 2.0E0*I_NAI_I4x2y_Gx2yz_a-1*I_NAI_G4x_Gx2yz;
  abcd[464] = 2.0E0*I_NAI_I4xyz_Gx2yz_a;
  abcd[465] = 2.0E0*I_NAI_I3x3y_Gx2yz_a-2*I_NAI_G3xy_Gx2yz;
  abcd[466] = 2.0E0*I_NAI_I3x2yz_Gx2yz_a-1*I_NAI_G3xz_Gx2yz;
  abcd[467] = 2.0E0*I_NAI_I3xy2z_Gx2yz_a;
  abcd[468] = 2.0E0*I_NAI_I2x4y_Gx2yz_a-3*I_NAI_G2x2y_Gx2yz;
  abcd[469] = 2.0E0*I_NAI_I2x3yz_Gx2yz_a-2*I_NAI_G2xyz_Gx2yz;
  abcd[470] = 2.0E0*I_NAI_I2x2y2z_Gx2yz_a-1*I_NAI_G2x2z_Gx2yz;
  abcd[471] = 2.0E0*I_NAI_I2xy3z_Gx2yz_a;
  abcd[472] = 2.0E0*I_NAI_Ix5y_Gx2yz_a-4*I_NAI_Gx3y_Gx2yz;
  abcd[473] = 2.0E0*I_NAI_Ix4yz_Gx2yz_a-3*I_NAI_Gx2yz_Gx2yz;
  abcd[474] = 2.0E0*I_NAI_Ix3y2z_Gx2yz_a-2*I_NAI_Gxy2z_Gx2yz;
  abcd[475] = 2.0E0*I_NAI_Ix2y3z_Gx2yz_a-1*I_NAI_Gx3z_Gx2yz;
  abcd[476] = 2.0E0*I_NAI_Ixy4z_Gx2yz_a;
  abcd[477] = 2.0E0*I_NAI_I6y_Gx2yz_a-5*I_NAI_G4y_Gx2yz;
  abcd[478] = 2.0E0*I_NAI_I5yz_Gx2yz_a-4*I_NAI_G3yz_Gx2yz;
  abcd[479] = 2.0E0*I_NAI_I4y2z_Gx2yz_a-3*I_NAI_G2y2z_Gx2yz;
  abcd[480] = 2.0E0*I_NAI_I3y3z_Gx2yz_a-2*I_NAI_Gy3z_Gx2yz;
  abcd[481] = 2.0E0*I_NAI_I2y4z_Gx2yz_a-1*I_NAI_G4z_Gx2yz;
  abcd[482] = 2.0E0*I_NAI_Iy5z_Gx2yz_a;
  abcd[483] = 2.0E0*I_NAI_I5xy_Gxy2z_a;
  abcd[484] = 2.0E0*I_NAI_I4x2y_Gxy2z_a-1*I_NAI_G4x_Gxy2z;
  abcd[485] = 2.0E0*I_NAI_I4xyz_Gxy2z_a;
  abcd[486] = 2.0E0*I_NAI_I3x3y_Gxy2z_a-2*I_NAI_G3xy_Gxy2z;
  abcd[487] = 2.0E0*I_NAI_I3x2yz_Gxy2z_a-1*I_NAI_G3xz_Gxy2z;
  abcd[488] = 2.0E0*I_NAI_I3xy2z_Gxy2z_a;
  abcd[489] = 2.0E0*I_NAI_I2x4y_Gxy2z_a-3*I_NAI_G2x2y_Gxy2z;
  abcd[490] = 2.0E0*I_NAI_I2x3yz_Gxy2z_a-2*I_NAI_G2xyz_Gxy2z;
  abcd[491] = 2.0E0*I_NAI_I2x2y2z_Gxy2z_a-1*I_NAI_G2x2z_Gxy2z;
  abcd[492] = 2.0E0*I_NAI_I2xy3z_Gxy2z_a;
  abcd[493] = 2.0E0*I_NAI_Ix5y_Gxy2z_a-4*I_NAI_Gx3y_Gxy2z;
  abcd[494] = 2.0E0*I_NAI_Ix4yz_Gxy2z_a-3*I_NAI_Gx2yz_Gxy2z;
  abcd[495] = 2.0E0*I_NAI_Ix3y2z_Gxy2z_a-2*I_NAI_Gxy2z_Gxy2z;
  abcd[496] = 2.0E0*I_NAI_Ix2y3z_Gxy2z_a-1*I_NAI_Gx3z_Gxy2z;
  abcd[497] = 2.0E0*I_NAI_Ixy4z_Gxy2z_a;
  abcd[498] = 2.0E0*I_NAI_I6y_Gxy2z_a-5*I_NAI_G4y_Gxy2z;
  abcd[499] = 2.0E0*I_NAI_I5yz_Gxy2z_a-4*I_NAI_G3yz_Gxy2z;
  abcd[500] = 2.0E0*I_NAI_I4y2z_Gxy2z_a-3*I_NAI_G2y2z_Gxy2z;
  abcd[501] = 2.0E0*I_NAI_I3y3z_Gxy2z_a-2*I_NAI_Gy3z_Gxy2z;
  abcd[502] = 2.0E0*I_NAI_I2y4z_Gxy2z_a-1*I_NAI_G4z_Gxy2z;
  abcd[503] = 2.0E0*I_NAI_Iy5z_Gxy2z_a;
  abcd[504] = 2.0E0*I_NAI_I5xy_Gx3z_a;
  abcd[505] = 2.0E0*I_NAI_I4x2y_Gx3z_a-1*I_NAI_G4x_Gx3z;
  abcd[506] = 2.0E0*I_NAI_I4xyz_Gx3z_a;
  abcd[507] = 2.0E0*I_NAI_I3x3y_Gx3z_a-2*I_NAI_G3xy_Gx3z;
  abcd[508] = 2.0E0*I_NAI_I3x2yz_Gx3z_a-1*I_NAI_G3xz_Gx3z;
  abcd[509] = 2.0E0*I_NAI_I3xy2z_Gx3z_a;
  abcd[510] = 2.0E0*I_NAI_I2x4y_Gx3z_a-3*I_NAI_G2x2y_Gx3z;
  abcd[511] = 2.0E0*I_NAI_I2x3yz_Gx3z_a-2*I_NAI_G2xyz_Gx3z;
  abcd[512] = 2.0E0*I_NAI_I2x2y2z_Gx3z_a-1*I_NAI_G2x2z_Gx3z;
  abcd[513] = 2.0E0*I_NAI_I2xy3z_Gx3z_a;
  abcd[514] = 2.0E0*I_NAI_Ix5y_Gx3z_a-4*I_NAI_Gx3y_Gx3z;
  abcd[515] = 2.0E0*I_NAI_Ix4yz_Gx3z_a-3*I_NAI_Gx2yz_Gx3z;
  abcd[516] = 2.0E0*I_NAI_Ix3y2z_Gx3z_a-2*I_NAI_Gxy2z_Gx3z;
  abcd[517] = 2.0E0*I_NAI_Ix2y3z_Gx3z_a-1*I_NAI_Gx3z_Gx3z;
  abcd[518] = 2.0E0*I_NAI_Ixy4z_Gx3z_a;
  abcd[519] = 2.0E0*I_NAI_I6y_Gx3z_a-5*I_NAI_G4y_Gx3z;
  abcd[520] = 2.0E0*I_NAI_I5yz_Gx3z_a-4*I_NAI_G3yz_Gx3z;
  abcd[521] = 2.0E0*I_NAI_I4y2z_Gx3z_a-3*I_NAI_G2y2z_Gx3z;
  abcd[522] = 2.0E0*I_NAI_I3y3z_Gx3z_a-2*I_NAI_Gy3z_Gx3z;
  abcd[523] = 2.0E0*I_NAI_I2y4z_Gx3z_a-1*I_NAI_G4z_Gx3z;
  abcd[524] = 2.0E0*I_NAI_Iy5z_Gx3z_a;
  abcd[525] = 2.0E0*I_NAI_I5xy_G4y_a;
  abcd[526] = 2.0E0*I_NAI_I4x2y_G4y_a-1*I_NAI_G4x_G4y;
  abcd[527] = 2.0E0*I_NAI_I4xyz_G4y_a;
  abcd[528] = 2.0E0*I_NAI_I3x3y_G4y_a-2*I_NAI_G3xy_G4y;
  abcd[529] = 2.0E0*I_NAI_I3x2yz_G4y_a-1*I_NAI_G3xz_G4y;
  abcd[530] = 2.0E0*I_NAI_I3xy2z_G4y_a;
  abcd[531] = 2.0E0*I_NAI_I2x4y_G4y_a-3*I_NAI_G2x2y_G4y;
  abcd[532] = 2.0E0*I_NAI_I2x3yz_G4y_a-2*I_NAI_G2xyz_G4y;
  abcd[533] = 2.0E0*I_NAI_I2x2y2z_G4y_a-1*I_NAI_G2x2z_G4y;
  abcd[534] = 2.0E0*I_NAI_I2xy3z_G4y_a;
  abcd[535] = 2.0E0*I_NAI_Ix5y_G4y_a-4*I_NAI_Gx3y_G4y;
  abcd[536] = 2.0E0*I_NAI_Ix4yz_G4y_a-3*I_NAI_Gx2yz_G4y;
  abcd[537] = 2.0E0*I_NAI_Ix3y2z_G4y_a-2*I_NAI_Gxy2z_G4y;
  abcd[538] = 2.0E0*I_NAI_Ix2y3z_G4y_a-1*I_NAI_Gx3z_G4y;
  abcd[539] = 2.0E0*I_NAI_Ixy4z_G4y_a;
  abcd[540] = 2.0E0*I_NAI_I6y_G4y_a-5*I_NAI_G4y_G4y;
  abcd[541] = 2.0E0*I_NAI_I5yz_G4y_a-4*I_NAI_G3yz_G4y;
  abcd[542] = 2.0E0*I_NAI_I4y2z_G4y_a-3*I_NAI_G2y2z_G4y;
  abcd[543] = 2.0E0*I_NAI_I3y3z_G4y_a-2*I_NAI_Gy3z_G4y;
  abcd[544] = 2.0E0*I_NAI_I2y4z_G4y_a-1*I_NAI_G4z_G4y;
  abcd[545] = 2.0E0*I_NAI_Iy5z_G4y_a;
  abcd[546] = 2.0E0*I_NAI_I5xy_G3yz_a;
  abcd[547] = 2.0E0*I_NAI_I4x2y_G3yz_a-1*I_NAI_G4x_G3yz;
  abcd[548] = 2.0E0*I_NAI_I4xyz_G3yz_a;
  abcd[549] = 2.0E0*I_NAI_I3x3y_G3yz_a-2*I_NAI_G3xy_G3yz;
  abcd[550] = 2.0E0*I_NAI_I3x2yz_G3yz_a-1*I_NAI_G3xz_G3yz;
  abcd[551] = 2.0E0*I_NAI_I3xy2z_G3yz_a;
  abcd[552] = 2.0E0*I_NAI_I2x4y_G3yz_a-3*I_NAI_G2x2y_G3yz;
  abcd[553] = 2.0E0*I_NAI_I2x3yz_G3yz_a-2*I_NAI_G2xyz_G3yz;
  abcd[554] = 2.0E0*I_NAI_I2x2y2z_G3yz_a-1*I_NAI_G2x2z_G3yz;
  abcd[555] = 2.0E0*I_NAI_I2xy3z_G3yz_a;
  abcd[556] = 2.0E0*I_NAI_Ix5y_G3yz_a-4*I_NAI_Gx3y_G3yz;
  abcd[557] = 2.0E0*I_NAI_Ix4yz_G3yz_a-3*I_NAI_Gx2yz_G3yz;
  abcd[558] = 2.0E0*I_NAI_Ix3y2z_G3yz_a-2*I_NAI_Gxy2z_G3yz;
  abcd[559] = 2.0E0*I_NAI_Ix2y3z_G3yz_a-1*I_NAI_Gx3z_G3yz;
  abcd[560] = 2.0E0*I_NAI_Ixy4z_G3yz_a;
  abcd[561] = 2.0E0*I_NAI_I6y_G3yz_a-5*I_NAI_G4y_G3yz;
  abcd[562] = 2.0E0*I_NAI_I5yz_G3yz_a-4*I_NAI_G3yz_G3yz;
  abcd[563] = 2.0E0*I_NAI_I4y2z_G3yz_a-3*I_NAI_G2y2z_G3yz;
  abcd[564] = 2.0E0*I_NAI_I3y3z_G3yz_a-2*I_NAI_Gy3z_G3yz;
  abcd[565] = 2.0E0*I_NAI_I2y4z_G3yz_a-1*I_NAI_G4z_G3yz;
  abcd[566] = 2.0E0*I_NAI_Iy5z_G3yz_a;
  abcd[567] = 2.0E0*I_NAI_I5xy_G2y2z_a;
  abcd[568] = 2.0E0*I_NAI_I4x2y_G2y2z_a-1*I_NAI_G4x_G2y2z;
  abcd[569] = 2.0E0*I_NAI_I4xyz_G2y2z_a;
  abcd[570] = 2.0E0*I_NAI_I3x3y_G2y2z_a-2*I_NAI_G3xy_G2y2z;
  abcd[571] = 2.0E0*I_NAI_I3x2yz_G2y2z_a-1*I_NAI_G3xz_G2y2z;
  abcd[572] = 2.0E0*I_NAI_I3xy2z_G2y2z_a;
  abcd[573] = 2.0E0*I_NAI_I2x4y_G2y2z_a-3*I_NAI_G2x2y_G2y2z;
  abcd[574] = 2.0E0*I_NAI_I2x3yz_G2y2z_a-2*I_NAI_G2xyz_G2y2z;
  abcd[575] = 2.0E0*I_NAI_I2x2y2z_G2y2z_a-1*I_NAI_G2x2z_G2y2z;
  abcd[576] = 2.0E0*I_NAI_I2xy3z_G2y2z_a;
  abcd[577] = 2.0E0*I_NAI_Ix5y_G2y2z_a-4*I_NAI_Gx3y_G2y2z;
  abcd[578] = 2.0E0*I_NAI_Ix4yz_G2y2z_a-3*I_NAI_Gx2yz_G2y2z;
  abcd[579] = 2.0E0*I_NAI_Ix3y2z_G2y2z_a-2*I_NAI_Gxy2z_G2y2z;
  abcd[580] = 2.0E0*I_NAI_Ix2y3z_G2y2z_a-1*I_NAI_Gx3z_G2y2z;
  abcd[581] = 2.0E0*I_NAI_Ixy4z_G2y2z_a;
  abcd[582] = 2.0E0*I_NAI_I6y_G2y2z_a-5*I_NAI_G4y_G2y2z;
  abcd[583] = 2.0E0*I_NAI_I5yz_G2y2z_a-4*I_NAI_G3yz_G2y2z;
  abcd[584] = 2.0E0*I_NAI_I4y2z_G2y2z_a-3*I_NAI_G2y2z_G2y2z;
  abcd[585] = 2.0E0*I_NAI_I3y3z_G2y2z_a-2*I_NAI_Gy3z_G2y2z;
  abcd[586] = 2.0E0*I_NAI_I2y4z_G2y2z_a-1*I_NAI_G4z_G2y2z;
  abcd[587] = 2.0E0*I_NAI_Iy5z_G2y2z_a;
  abcd[588] = 2.0E0*I_NAI_I5xy_Gy3z_a;
  abcd[589] = 2.0E0*I_NAI_I4x2y_Gy3z_a-1*I_NAI_G4x_Gy3z;
  abcd[590] = 2.0E0*I_NAI_I4xyz_Gy3z_a;
  abcd[591] = 2.0E0*I_NAI_I3x3y_Gy3z_a-2*I_NAI_G3xy_Gy3z;
  abcd[592] = 2.0E0*I_NAI_I3x2yz_Gy3z_a-1*I_NAI_G3xz_Gy3z;
  abcd[593] = 2.0E0*I_NAI_I3xy2z_Gy3z_a;
  abcd[594] = 2.0E0*I_NAI_I2x4y_Gy3z_a-3*I_NAI_G2x2y_Gy3z;
  abcd[595] = 2.0E0*I_NAI_I2x3yz_Gy3z_a-2*I_NAI_G2xyz_Gy3z;
  abcd[596] = 2.0E0*I_NAI_I2x2y2z_Gy3z_a-1*I_NAI_G2x2z_Gy3z;
  abcd[597] = 2.0E0*I_NAI_I2xy3z_Gy3z_a;
  abcd[598] = 2.0E0*I_NAI_Ix5y_Gy3z_a-4*I_NAI_Gx3y_Gy3z;
  abcd[599] = 2.0E0*I_NAI_Ix4yz_Gy3z_a-3*I_NAI_Gx2yz_Gy3z;
  abcd[600] = 2.0E0*I_NAI_Ix3y2z_Gy3z_a-2*I_NAI_Gxy2z_Gy3z;
  abcd[601] = 2.0E0*I_NAI_Ix2y3z_Gy3z_a-1*I_NAI_Gx3z_Gy3z;
  abcd[602] = 2.0E0*I_NAI_Ixy4z_Gy3z_a;
  abcd[603] = 2.0E0*I_NAI_I6y_Gy3z_a-5*I_NAI_G4y_Gy3z;
  abcd[604] = 2.0E0*I_NAI_I5yz_Gy3z_a-4*I_NAI_G3yz_Gy3z;
  abcd[605] = 2.0E0*I_NAI_I4y2z_Gy3z_a-3*I_NAI_G2y2z_Gy3z;
  abcd[606] = 2.0E0*I_NAI_I3y3z_Gy3z_a-2*I_NAI_Gy3z_Gy3z;
  abcd[607] = 2.0E0*I_NAI_I2y4z_Gy3z_a-1*I_NAI_G4z_Gy3z;
  abcd[608] = 2.0E0*I_NAI_Iy5z_Gy3z_a;
  abcd[609] = 2.0E0*I_NAI_I5xy_G4z_a;
  abcd[610] = 2.0E0*I_NAI_I4x2y_G4z_a-1*I_NAI_G4x_G4z;
  abcd[611] = 2.0E0*I_NAI_I4xyz_G4z_a;
  abcd[612] = 2.0E0*I_NAI_I3x3y_G4z_a-2*I_NAI_G3xy_G4z;
  abcd[613] = 2.0E0*I_NAI_I3x2yz_G4z_a-1*I_NAI_G3xz_G4z;
  abcd[614] = 2.0E0*I_NAI_I3xy2z_G4z_a;
  abcd[615] = 2.0E0*I_NAI_I2x4y_G4z_a-3*I_NAI_G2x2y_G4z;
  abcd[616] = 2.0E0*I_NAI_I2x3yz_G4z_a-2*I_NAI_G2xyz_G4z;
  abcd[617] = 2.0E0*I_NAI_I2x2y2z_G4z_a-1*I_NAI_G2x2z_G4z;
  abcd[618] = 2.0E0*I_NAI_I2xy3z_G4z_a;
  abcd[619] = 2.0E0*I_NAI_Ix5y_G4z_a-4*I_NAI_Gx3y_G4z;
  abcd[620] = 2.0E0*I_NAI_Ix4yz_G4z_a-3*I_NAI_Gx2yz_G4z;
  abcd[621] = 2.0E0*I_NAI_Ix3y2z_G4z_a-2*I_NAI_Gxy2z_G4z;
  abcd[622] = 2.0E0*I_NAI_Ix2y3z_G4z_a-1*I_NAI_Gx3z_G4z;
  abcd[623] = 2.0E0*I_NAI_Ixy4z_G4z_a;
  abcd[624] = 2.0E0*I_NAI_I6y_G4z_a-5*I_NAI_G4y_G4z;
  abcd[625] = 2.0E0*I_NAI_I5yz_G4z_a-4*I_NAI_G3yz_G4z;
  abcd[626] = 2.0E0*I_NAI_I4y2z_G4z_a-3*I_NAI_G2y2z_G4z;
  abcd[627] = 2.0E0*I_NAI_I3y3z_G4z_a-2*I_NAI_Gy3z_G4z;
  abcd[628] = 2.0E0*I_NAI_I2y4z_G4z_a-1*I_NAI_G4z_G4z;
  abcd[629] = 2.0E0*I_NAI_Iy5z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_G_a
   * RHS shell quartet name: SQ_NAI_G_G
   ************************************************************/
  abcd[630] = 2.0E0*I_NAI_I5xz_G4x_a;
  abcd[631] = 2.0E0*I_NAI_I4xyz_G4x_a;
  abcd[632] = 2.0E0*I_NAI_I4x2z_G4x_a-1*I_NAI_G4x_G4x;
  abcd[633] = 2.0E0*I_NAI_I3x2yz_G4x_a;
  abcd[634] = 2.0E0*I_NAI_I3xy2z_G4x_a-1*I_NAI_G3xy_G4x;
  abcd[635] = 2.0E0*I_NAI_I3x3z_G4x_a-2*I_NAI_G3xz_G4x;
  abcd[636] = 2.0E0*I_NAI_I2x3yz_G4x_a;
  abcd[637] = 2.0E0*I_NAI_I2x2y2z_G4x_a-1*I_NAI_G2x2y_G4x;
  abcd[638] = 2.0E0*I_NAI_I2xy3z_G4x_a-2*I_NAI_G2xyz_G4x;
  abcd[639] = 2.0E0*I_NAI_I2x4z_G4x_a-3*I_NAI_G2x2z_G4x;
  abcd[640] = 2.0E0*I_NAI_Ix4yz_G4x_a;
  abcd[641] = 2.0E0*I_NAI_Ix3y2z_G4x_a-1*I_NAI_Gx3y_G4x;
  abcd[642] = 2.0E0*I_NAI_Ix2y3z_G4x_a-2*I_NAI_Gx2yz_G4x;
  abcd[643] = 2.0E0*I_NAI_Ixy4z_G4x_a-3*I_NAI_Gxy2z_G4x;
  abcd[644] = 2.0E0*I_NAI_Ix5z_G4x_a-4*I_NAI_Gx3z_G4x;
  abcd[645] = 2.0E0*I_NAI_I5yz_G4x_a;
  abcd[646] = 2.0E0*I_NAI_I4y2z_G4x_a-1*I_NAI_G4y_G4x;
  abcd[647] = 2.0E0*I_NAI_I3y3z_G4x_a-2*I_NAI_G3yz_G4x;
  abcd[648] = 2.0E0*I_NAI_I2y4z_G4x_a-3*I_NAI_G2y2z_G4x;
  abcd[649] = 2.0E0*I_NAI_Iy5z_G4x_a-4*I_NAI_Gy3z_G4x;
  abcd[650] = 2.0E0*I_NAI_I6z_G4x_a-5*I_NAI_G4z_G4x;
  abcd[651] = 2.0E0*I_NAI_I5xz_G3xy_a;
  abcd[652] = 2.0E0*I_NAI_I4xyz_G3xy_a;
  abcd[653] = 2.0E0*I_NAI_I4x2z_G3xy_a-1*I_NAI_G4x_G3xy;
  abcd[654] = 2.0E0*I_NAI_I3x2yz_G3xy_a;
  abcd[655] = 2.0E0*I_NAI_I3xy2z_G3xy_a-1*I_NAI_G3xy_G3xy;
  abcd[656] = 2.0E0*I_NAI_I3x3z_G3xy_a-2*I_NAI_G3xz_G3xy;
  abcd[657] = 2.0E0*I_NAI_I2x3yz_G3xy_a;
  abcd[658] = 2.0E0*I_NAI_I2x2y2z_G3xy_a-1*I_NAI_G2x2y_G3xy;
  abcd[659] = 2.0E0*I_NAI_I2xy3z_G3xy_a-2*I_NAI_G2xyz_G3xy;
  abcd[660] = 2.0E0*I_NAI_I2x4z_G3xy_a-3*I_NAI_G2x2z_G3xy;
  abcd[661] = 2.0E0*I_NAI_Ix4yz_G3xy_a;
  abcd[662] = 2.0E0*I_NAI_Ix3y2z_G3xy_a-1*I_NAI_Gx3y_G3xy;
  abcd[663] = 2.0E0*I_NAI_Ix2y3z_G3xy_a-2*I_NAI_Gx2yz_G3xy;
  abcd[664] = 2.0E0*I_NAI_Ixy4z_G3xy_a-3*I_NAI_Gxy2z_G3xy;
  abcd[665] = 2.0E0*I_NAI_Ix5z_G3xy_a-4*I_NAI_Gx3z_G3xy;
  abcd[666] = 2.0E0*I_NAI_I5yz_G3xy_a;
  abcd[667] = 2.0E0*I_NAI_I4y2z_G3xy_a-1*I_NAI_G4y_G3xy;
  abcd[668] = 2.0E0*I_NAI_I3y3z_G3xy_a-2*I_NAI_G3yz_G3xy;
  abcd[669] = 2.0E0*I_NAI_I2y4z_G3xy_a-3*I_NAI_G2y2z_G3xy;
  abcd[670] = 2.0E0*I_NAI_Iy5z_G3xy_a-4*I_NAI_Gy3z_G3xy;
  abcd[671] = 2.0E0*I_NAI_I6z_G3xy_a-5*I_NAI_G4z_G3xy;
  abcd[672] = 2.0E0*I_NAI_I5xz_G3xz_a;
  abcd[673] = 2.0E0*I_NAI_I4xyz_G3xz_a;
  abcd[674] = 2.0E0*I_NAI_I4x2z_G3xz_a-1*I_NAI_G4x_G3xz;
  abcd[675] = 2.0E0*I_NAI_I3x2yz_G3xz_a;
  abcd[676] = 2.0E0*I_NAI_I3xy2z_G3xz_a-1*I_NAI_G3xy_G3xz;
  abcd[677] = 2.0E0*I_NAI_I3x3z_G3xz_a-2*I_NAI_G3xz_G3xz;
  abcd[678] = 2.0E0*I_NAI_I2x3yz_G3xz_a;
  abcd[679] = 2.0E0*I_NAI_I2x2y2z_G3xz_a-1*I_NAI_G2x2y_G3xz;
  abcd[680] = 2.0E0*I_NAI_I2xy3z_G3xz_a-2*I_NAI_G2xyz_G3xz;
  abcd[681] = 2.0E0*I_NAI_I2x4z_G3xz_a-3*I_NAI_G2x2z_G3xz;
  abcd[682] = 2.0E0*I_NAI_Ix4yz_G3xz_a;
  abcd[683] = 2.0E0*I_NAI_Ix3y2z_G3xz_a-1*I_NAI_Gx3y_G3xz;
  abcd[684] = 2.0E0*I_NAI_Ix2y3z_G3xz_a-2*I_NAI_Gx2yz_G3xz;
  abcd[685] = 2.0E0*I_NAI_Ixy4z_G3xz_a-3*I_NAI_Gxy2z_G3xz;
  abcd[686] = 2.0E0*I_NAI_Ix5z_G3xz_a-4*I_NAI_Gx3z_G3xz;
  abcd[687] = 2.0E0*I_NAI_I5yz_G3xz_a;
  abcd[688] = 2.0E0*I_NAI_I4y2z_G3xz_a-1*I_NAI_G4y_G3xz;
  abcd[689] = 2.0E0*I_NAI_I3y3z_G3xz_a-2*I_NAI_G3yz_G3xz;
  abcd[690] = 2.0E0*I_NAI_I2y4z_G3xz_a-3*I_NAI_G2y2z_G3xz;
  abcd[691] = 2.0E0*I_NAI_Iy5z_G3xz_a-4*I_NAI_Gy3z_G3xz;
  abcd[692] = 2.0E0*I_NAI_I6z_G3xz_a-5*I_NAI_G4z_G3xz;
  abcd[693] = 2.0E0*I_NAI_I5xz_G2x2y_a;
  abcd[694] = 2.0E0*I_NAI_I4xyz_G2x2y_a;
  abcd[695] = 2.0E0*I_NAI_I4x2z_G2x2y_a-1*I_NAI_G4x_G2x2y;
  abcd[696] = 2.0E0*I_NAI_I3x2yz_G2x2y_a;
  abcd[697] = 2.0E0*I_NAI_I3xy2z_G2x2y_a-1*I_NAI_G3xy_G2x2y;
  abcd[698] = 2.0E0*I_NAI_I3x3z_G2x2y_a-2*I_NAI_G3xz_G2x2y;
  abcd[699] = 2.0E0*I_NAI_I2x3yz_G2x2y_a;
  abcd[700] = 2.0E0*I_NAI_I2x2y2z_G2x2y_a-1*I_NAI_G2x2y_G2x2y;
  abcd[701] = 2.0E0*I_NAI_I2xy3z_G2x2y_a-2*I_NAI_G2xyz_G2x2y;
  abcd[702] = 2.0E0*I_NAI_I2x4z_G2x2y_a-3*I_NAI_G2x2z_G2x2y;
  abcd[703] = 2.0E0*I_NAI_Ix4yz_G2x2y_a;
  abcd[704] = 2.0E0*I_NAI_Ix3y2z_G2x2y_a-1*I_NAI_Gx3y_G2x2y;
  abcd[705] = 2.0E0*I_NAI_Ix2y3z_G2x2y_a-2*I_NAI_Gx2yz_G2x2y;
  abcd[706] = 2.0E0*I_NAI_Ixy4z_G2x2y_a-3*I_NAI_Gxy2z_G2x2y;
  abcd[707] = 2.0E0*I_NAI_Ix5z_G2x2y_a-4*I_NAI_Gx3z_G2x2y;
  abcd[708] = 2.0E0*I_NAI_I5yz_G2x2y_a;
  abcd[709] = 2.0E0*I_NAI_I4y2z_G2x2y_a-1*I_NAI_G4y_G2x2y;
  abcd[710] = 2.0E0*I_NAI_I3y3z_G2x2y_a-2*I_NAI_G3yz_G2x2y;
  abcd[711] = 2.0E0*I_NAI_I2y4z_G2x2y_a-3*I_NAI_G2y2z_G2x2y;
  abcd[712] = 2.0E0*I_NAI_Iy5z_G2x2y_a-4*I_NAI_Gy3z_G2x2y;
  abcd[713] = 2.0E0*I_NAI_I6z_G2x2y_a-5*I_NAI_G4z_G2x2y;
  abcd[714] = 2.0E0*I_NAI_I5xz_G2xyz_a;
  abcd[715] = 2.0E0*I_NAI_I4xyz_G2xyz_a;
  abcd[716] = 2.0E0*I_NAI_I4x2z_G2xyz_a-1*I_NAI_G4x_G2xyz;
  abcd[717] = 2.0E0*I_NAI_I3x2yz_G2xyz_a;
  abcd[718] = 2.0E0*I_NAI_I3xy2z_G2xyz_a-1*I_NAI_G3xy_G2xyz;
  abcd[719] = 2.0E0*I_NAI_I3x3z_G2xyz_a-2*I_NAI_G3xz_G2xyz;
  abcd[720] = 2.0E0*I_NAI_I2x3yz_G2xyz_a;
  abcd[721] = 2.0E0*I_NAI_I2x2y2z_G2xyz_a-1*I_NAI_G2x2y_G2xyz;
  abcd[722] = 2.0E0*I_NAI_I2xy3z_G2xyz_a-2*I_NAI_G2xyz_G2xyz;
  abcd[723] = 2.0E0*I_NAI_I2x4z_G2xyz_a-3*I_NAI_G2x2z_G2xyz;
  abcd[724] = 2.0E0*I_NAI_Ix4yz_G2xyz_a;
  abcd[725] = 2.0E0*I_NAI_Ix3y2z_G2xyz_a-1*I_NAI_Gx3y_G2xyz;
  abcd[726] = 2.0E0*I_NAI_Ix2y3z_G2xyz_a-2*I_NAI_Gx2yz_G2xyz;
  abcd[727] = 2.0E0*I_NAI_Ixy4z_G2xyz_a-3*I_NAI_Gxy2z_G2xyz;
  abcd[728] = 2.0E0*I_NAI_Ix5z_G2xyz_a-4*I_NAI_Gx3z_G2xyz;
  abcd[729] = 2.0E0*I_NAI_I5yz_G2xyz_a;
  abcd[730] = 2.0E0*I_NAI_I4y2z_G2xyz_a-1*I_NAI_G4y_G2xyz;
  abcd[731] = 2.0E0*I_NAI_I3y3z_G2xyz_a-2*I_NAI_G3yz_G2xyz;
  abcd[732] = 2.0E0*I_NAI_I2y4z_G2xyz_a-3*I_NAI_G2y2z_G2xyz;
  abcd[733] = 2.0E0*I_NAI_Iy5z_G2xyz_a-4*I_NAI_Gy3z_G2xyz;
  abcd[734] = 2.0E0*I_NAI_I6z_G2xyz_a-5*I_NAI_G4z_G2xyz;
  abcd[735] = 2.0E0*I_NAI_I5xz_G2x2z_a;
  abcd[736] = 2.0E0*I_NAI_I4xyz_G2x2z_a;
  abcd[737] = 2.0E0*I_NAI_I4x2z_G2x2z_a-1*I_NAI_G4x_G2x2z;
  abcd[738] = 2.0E0*I_NAI_I3x2yz_G2x2z_a;
  abcd[739] = 2.0E0*I_NAI_I3xy2z_G2x2z_a-1*I_NAI_G3xy_G2x2z;
  abcd[740] = 2.0E0*I_NAI_I3x3z_G2x2z_a-2*I_NAI_G3xz_G2x2z;
  abcd[741] = 2.0E0*I_NAI_I2x3yz_G2x2z_a;
  abcd[742] = 2.0E0*I_NAI_I2x2y2z_G2x2z_a-1*I_NAI_G2x2y_G2x2z;
  abcd[743] = 2.0E0*I_NAI_I2xy3z_G2x2z_a-2*I_NAI_G2xyz_G2x2z;
  abcd[744] = 2.0E0*I_NAI_I2x4z_G2x2z_a-3*I_NAI_G2x2z_G2x2z;
  abcd[745] = 2.0E0*I_NAI_Ix4yz_G2x2z_a;
  abcd[746] = 2.0E0*I_NAI_Ix3y2z_G2x2z_a-1*I_NAI_Gx3y_G2x2z;
  abcd[747] = 2.0E0*I_NAI_Ix2y3z_G2x2z_a-2*I_NAI_Gx2yz_G2x2z;
  abcd[748] = 2.0E0*I_NAI_Ixy4z_G2x2z_a-3*I_NAI_Gxy2z_G2x2z;
  abcd[749] = 2.0E0*I_NAI_Ix5z_G2x2z_a-4*I_NAI_Gx3z_G2x2z;
  abcd[750] = 2.0E0*I_NAI_I5yz_G2x2z_a;
  abcd[751] = 2.0E0*I_NAI_I4y2z_G2x2z_a-1*I_NAI_G4y_G2x2z;
  abcd[752] = 2.0E0*I_NAI_I3y3z_G2x2z_a-2*I_NAI_G3yz_G2x2z;
  abcd[753] = 2.0E0*I_NAI_I2y4z_G2x2z_a-3*I_NAI_G2y2z_G2x2z;
  abcd[754] = 2.0E0*I_NAI_Iy5z_G2x2z_a-4*I_NAI_Gy3z_G2x2z;
  abcd[755] = 2.0E0*I_NAI_I6z_G2x2z_a-5*I_NAI_G4z_G2x2z;
  abcd[756] = 2.0E0*I_NAI_I5xz_Gx3y_a;
  abcd[757] = 2.0E0*I_NAI_I4xyz_Gx3y_a;
  abcd[758] = 2.0E0*I_NAI_I4x2z_Gx3y_a-1*I_NAI_G4x_Gx3y;
  abcd[759] = 2.0E0*I_NAI_I3x2yz_Gx3y_a;
  abcd[760] = 2.0E0*I_NAI_I3xy2z_Gx3y_a-1*I_NAI_G3xy_Gx3y;
  abcd[761] = 2.0E0*I_NAI_I3x3z_Gx3y_a-2*I_NAI_G3xz_Gx3y;
  abcd[762] = 2.0E0*I_NAI_I2x3yz_Gx3y_a;
  abcd[763] = 2.0E0*I_NAI_I2x2y2z_Gx3y_a-1*I_NAI_G2x2y_Gx3y;
  abcd[764] = 2.0E0*I_NAI_I2xy3z_Gx3y_a-2*I_NAI_G2xyz_Gx3y;
  abcd[765] = 2.0E0*I_NAI_I2x4z_Gx3y_a-3*I_NAI_G2x2z_Gx3y;
  abcd[766] = 2.0E0*I_NAI_Ix4yz_Gx3y_a;
  abcd[767] = 2.0E0*I_NAI_Ix3y2z_Gx3y_a-1*I_NAI_Gx3y_Gx3y;
  abcd[768] = 2.0E0*I_NAI_Ix2y3z_Gx3y_a-2*I_NAI_Gx2yz_Gx3y;
  abcd[769] = 2.0E0*I_NAI_Ixy4z_Gx3y_a-3*I_NAI_Gxy2z_Gx3y;
  abcd[770] = 2.0E0*I_NAI_Ix5z_Gx3y_a-4*I_NAI_Gx3z_Gx3y;
  abcd[771] = 2.0E0*I_NAI_I5yz_Gx3y_a;
  abcd[772] = 2.0E0*I_NAI_I4y2z_Gx3y_a-1*I_NAI_G4y_Gx3y;
  abcd[773] = 2.0E0*I_NAI_I3y3z_Gx3y_a-2*I_NAI_G3yz_Gx3y;
  abcd[774] = 2.0E0*I_NAI_I2y4z_Gx3y_a-3*I_NAI_G2y2z_Gx3y;
  abcd[775] = 2.0E0*I_NAI_Iy5z_Gx3y_a-4*I_NAI_Gy3z_Gx3y;
  abcd[776] = 2.0E0*I_NAI_I6z_Gx3y_a-5*I_NAI_G4z_Gx3y;
  abcd[777] = 2.0E0*I_NAI_I5xz_Gx2yz_a;
  abcd[778] = 2.0E0*I_NAI_I4xyz_Gx2yz_a;
  abcd[779] = 2.0E0*I_NAI_I4x2z_Gx2yz_a-1*I_NAI_G4x_Gx2yz;
  abcd[780] = 2.0E0*I_NAI_I3x2yz_Gx2yz_a;
  abcd[781] = 2.0E0*I_NAI_I3xy2z_Gx2yz_a-1*I_NAI_G3xy_Gx2yz;
  abcd[782] = 2.0E0*I_NAI_I3x3z_Gx2yz_a-2*I_NAI_G3xz_Gx2yz;
  abcd[783] = 2.0E0*I_NAI_I2x3yz_Gx2yz_a;
  abcd[784] = 2.0E0*I_NAI_I2x2y2z_Gx2yz_a-1*I_NAI_G2x2y_Gx2yz;
  abcd[785] = 2.0E0*I_NAI_I2xy3z_Gx2yz_a-2*I_NAI_G2xyz_Gx2yz;
  abcd[786] = 2.0E0*I_NAI_I2x4z_Gx2yz_a-3*I_NAI_G2x2z_Gx2yz;
  abcd[787] = 2.0E0*I_NAI_Ix4yz_Gx2yz_a;
  abcd[788] = 2.0E0*I_NAI_Ix3y2z_Gx2yz_a-1*I_NAI_Gx3y_Gx2yz;
  abcd[789] = 2.0E0*I_NAI_Ix2y3z_Gx2yz_a-2*I_NAI_Gx2yz_Gx2yz;
  abcd[790] = 2.0E0*I_NAI_Ixy4z_Gx2yz_a-3*I_NAI_Gxy2z_Gx2yz;
  abcd[791] = 2.0E0*I_NAI_Ix5z_Gx2yz_a-4*I_NAI_Gx3z_Gx2yz;
  abcd[792] = 2.0E0*I_NAI_I5yz_Gx2yz_a;
  abcd[793] = 2.0E0*I_NAI_I4y2z_Gx2yz_a-1*I_NAI_G4y_Gx2yz;
  abcd[794] = 2.0E0*I_NAI_I3y3z_Gx2yz_a-2*I_NAI_G3yz_Gx2yz;
  abcd[795] = 2.0E0*I_NAI_I2y4z_Gx2yz_a-3*I_NAI_G2y2z_Gx2yz;
  abcd[796] = 2.0E0*I_NAI_Iy5z_Gx2yz_a-4*I_NAI_Gy3z_Gx2yz;
  abcd[797] = 2.0E0*I_NAI_I6z_Gx2yz_a-5*I_NAI_G4z_Gx2yz;
  abcd[798] = 2.0E0*I_NAI_I5xz_Gxy2z_a;
  abcd[799] = 2.0E0*I_NAI_I4xyz_Gxy2z_a;
  abcd[800] = 2.0E0*I_NAI_I4x2z_Gxy2z_a-1*I_NAI_G4x_Gxy2z;
  abcd[801] = 2.0E0*I_NAI_I3x2yz_Gxy2z_a;
  abcd[802] = 2.0E0*I_NAI_I3xy2z_Gxy2z_a-1*I_NAI_G3xy_Gxy2z;
  abcd[803] = 2.0E0*I_NAI_I3x3z_Gxy2z_a-2*I_NAI_G3xz_Gxy2z;
  abcd[804] = 2.0E0*I_NAI_I2x3yz_Gxy2z_a;
  abcd[805] = 2.0E0*I_NAI_I2x2y2z_Gxy2z_a-1*I_NAI_G2x2y_Gxy2z;
  abcd[806] = 2.0E0*I_NAI_I2xy3z_Gxy2z_a-2*I_NAI_G2xyz_Gxy2z;
  abcd[807] = 2.0E0*I_NAI_I2x4z_Gxy2z_a-3*I_NAI_G2x2z_Gxy2z;
  abcd[808] = 2.0E0*I_NAI_Ix4yz_Gxy2z_a;
  abcd[809] = 2.0E0*I_NAI_Ix3y2z_Gxy2z_a-1*I_NAI_Gx3y_Gxy2z;
  abcd[810] = 2.0E0*I_NAI_Ix2y3z_Gxy2z_a-2*I_NAI_Gx2yz_Gxy2z;
  abcd[811] = 2.0E0*I_NAI_Ixy4z_Gxy2z_a-3*I_NAI_Gxy2z_Gxy2z;
  abcd[812] = 2.0E0*I_NAI_Ix5z_Gxy2z_a-4*I_NAI_Gx3z_Gxy2z;
  abcd[813] = 2.0E0*I_NAI_I5yz_Gxy2z_a;
  abcd[814] = 2.0E0*I_NAI_I4y2z_Gxy2z_a-1*I_NAI_G4y_Gxy2z;
  abcd[815] = 2.0E0*I_NAI_I3y3z_Gxy2z_a-2*I_NAI_G3yz_Gxy2z;
  abcd[816] = 2.0E0*I_NAI_I2y4z_Gxy2z_a-3*I_NAI_G2y2z_Gxy2z;
  abcd[817] = 2.0E0*I_NAI_Iy5z_Gxy2z_a-4*I_NAI_Gy3z_Gxy2z;
  abcd[818] = 2.0E0*I_NAI_I6z_Gxy2z_a-5*I_NAI_G4z_Gxy2z;
  abcd[819] = 2.0E0*I_NAI_I5xz_Gx3z_a;
  abcd[820] = 2.0E0*I_NAI_I4xyz_Gx3z_a;
  abcd[821] = 2.0E0*I_NAI_I4x2z_Gx3z_a-1*I_NAI_G4x_Gx3z;
  abcd[822] = 2.0E0*I_NAI_I3x2yz_Gx3z_a;
  abcd[823] = 2.0E0*I_NAI_I3xy2z_Gx3z_a-1*I_NAI_G3xy_Gx3z;
  abcd[824] = 2.0E0*I_NAI_I3x3z_Gx3z_a-2*I_NAI_G3xz_Gx3z;
  abcd[825] = 2.0E0*I_NAI_I2x3yz_Gx3z_a;
  abcd[826] = 2.0E0*I_NAI_I2x2y2z_Gx3z_a-1*I_NAI_G2x2y_Gx3z;
  abcd[827] = 2.0E0*I_NAI_I2xy3z_Gx3z_a-2*I_NAI_G2xyz_Gx3z;
  abcd[828] = 2.0E0*I_NAI_I2x4z_Gx3z_a-3*I_NAI_G2x2z_Gx3z;
  abcd[829] = 2.0E0*I_NAI_Ix4yz_Gx3z_a;
  abcd[830] = 2.0E0*I_NAI_Ix3y2z_Gx3z_a-1*I_NAI_Gx3y_Gx3z;
  abcd[831] = 2.0E0*I_NAI_Ix2y3z_Gx3z_a-2*I_NAI_Gx2yz_Gx3z;
  abcd[832] = 2.0E0*I_NAI_Ixy4z_Gx3z_a-3*I_NAI_Gxy2z_Gx3z;
  abcd[833] = 2.0E0*I_NAI_Ix5z_Gx3z_a-4*I_NAI_Gx3z_Gx3z;
  abcd[834] = 2.0E0*I_NAI_I5yz_Gx3z_a;
  abcd[835] = 2.0E0*I_NAI_I4y2z_Gx3z_a-1*I_NAI_G4y_Gx3z;
  abcd[836] = 2.0E0*I_NAI_I3y3z_Gx3z_a-2*I_NAI_G3yz_Gx3z;
  abcd[837] = 2.0E0*I_NAI_I2y4z_Gx3z_a-3*I_NAI_G2y2z_Gx3z;
  abcd[838] = 2.0E0*I_NAI_Iy5z_Gx3z_a-4*I_NAI_Gy3z_Gx3z;
  abcd[839] = 2.0E0*I_NAI_I6z_Gx3z_a-5*I_NAI_G4z_Gx3z;
  abcd[840] = 2.0E0*I_NAI_I5xz_G4y_a;
  abcd[841] = 2.0E0*I_NAI_I4xyz_G4y_a;
  abcd[842] = 2.0E0*I_NAI_I4x2z_G4y_a-1*I_NAI_G4x_G4y;
  abcd[843] = 2.0E0*I_NAI_I3x2yz_G4y_a;
  abcd[844] = 2.0E0*I_NAI_I3xy2z_G4y_a-1*I_NAI_G3xy_G4y;
  abcd[845] = 2.0E0*I_NAI_I3x3z_G4y_a-2*I_NAI_G3xz_G4y;
  abcd[846] = 2.0E0*I_NAI_I2x3yz_G4y_a;
  abcd[847] = 2.0E0*I_NAI_I2x2y2z_G4y_a-1*I_NAI_G2x2y_G4y;
  abcd[848] = 2.0E0*I_NAI_I2xy3z_G4y_a-2*I_NAI_G2xyz_G4y;
  abcd[849] = 2.0E0*I_NAI_I2x4z_G4y_a-3*I_NAI_G2x2z_G4y;
  abcd[850] = 2.0E0*I_NAI_Ix4yz_G4y_a;
  abcd[851] = 2.0E0*I_NAI_Ix3y2z_G4y_a-1*I_NAI_Gx3y_G4y;
  abcd[852] = 2.0E0*I_NAI_Ix2y3z_G4y_a-2*I_NAI_Gx2yz_G4y;
  abcd[853] = 2.0E0*I_NAI_Ixy4z_G4y_a-3*I_NAI_Gxy2z_G4y;
  abcd[854] = 2.0E0*I_NAI_Ix5z_G4y_a-4*I_NAI_Gx3z_G4y;
  abcd[855] = 2.0E0*I_NAI_I5yz_G4y_a;
  abcd[856] = 2.0E0*I_NAI_I4y2z_G4y_a-1*I_NAI_G4y_G4y;
  abcd[857] = 2.0E0*I_NAI_I3y3z_G4y_a-2*I_NAI_G3yz_G4y;
  abcd[858] = 2.0E0*I_NAI_I2y4z_G4y_a-3*I_NAI_G2y2z_G4y;
  abcd[859] = 2.0E0*I_NAI_Iy5z_G4y_a-4*I_NAI_Gy3z_G4y;
  abcd[860] = 2.0E0*I_NAI_I6z_G4y_a-5*I_NAI_G4z_G4y;
  abcd[861] = 2.0E0*I_NAI_I5xz_G3yz_a;
  abcd[862] = 2.0E0*I_NAI_I4xyz_G3yz_a;
  abcd[863] = 2.0E0*I_NAI_I4x2z_G3yz_a-1*I_NAI_G4x_G3yz;
  abcd[864] = 2.0E0*I_NAI_I3x2yz_G3yz_a;
  abcd[865] = 2.0E0*I_NAI_I3xy2z_G3yz_a-1*I_NAI_G3xy_G3yz;
  abcd[866] = 2.0E0*I_NAI_I3x3z_G3yz_a-2*I_NAI_G3xz_G3yz;
  abcd[867] = 2.0E0*I_NAI_I2x3yz_G3yz_a;
  abcd[868] = 2.0E0*I_NAI_I2x2y2z_G3yz_a-1*I_NAI_G2x2y_G3yz;
  abcd[869] = 2.0E0*I_NAI_I2xy3z_G3yz_a-2*I_NAI_G2xyz_G3yz;
  abcd[870] = 2.0E0*I_NAI_I2x4z_G3yz_a-3*I_NAI_G2x2z_G3yz;
  abcd[871] = 2.0E0*I_NAI_Ix4yz_G3yz_a;
  abcd[872] = 2.0E0*I_NAI_Ix3y2z_G3yz_a-1*I_NAI_Gx3y_G3yz;
  abcd[873] = 2.0E0*I_NAI_Ix2y3z_G3yz_a-2*I_NAI_Gx2yz_G3yz;
  abcd[874] = 2.0E0*I_NAI_Ixy4z_G3yz_a-3*I_NAI_Gxy2z_G3yz;
  abcd[875] = 2.0E0*I_NAI_Ix5z_G3yz_a-4*I_NAI_Gx3z_G3yz;
  abcd[876] = 2.0E0*I_NAI_I5yz_G3yz_a;
  abcd[877] = 2.0E0*I_NAI_I4y2z_G3yz_a-1*I_NAI_G4y_G3yz;
  abcd[878] = 2.0E0*I_NAI_I3y3z_G3yz_a-2*I_NAI_G3yz_G3yz;
  abcd[879] = 2.0E0*I_NAI_I2y4z_G3yz_a-3*I_NAI_G2y2z_G3yz;
  abcd[880] = 2.0E0*I_NAI_Iy5z_G3yz_a-4*I_NAI_Gy3z_G3yz;
  abcd[881] = 2.0E0*I_NAI_I6z_G3yz_a-5*I_NAI_G4z_G3yz;
  abcd[882] = 2.0E0*I_NAI_I5xz_G2y2z_a;
  abcd[883] = 2.0E0*I_NAI_I4xyz_G2y2z_a;
  abcd[884] = 2.0E0*I_NAI_I4x2z_G2y2z_a-1*I_NAI_G4x_G2y2z;
  abcd[885] = 2.0E0*I_NAI_I3x2yz_G2y2z_a;
  abcd[886] = 2.0E0*I_NAI_I3xy2z_G2y2z_a-1*I_NAI_G3xy_G2y2z;
  abcd[887] = 2.0E0*I_NAI_I3x3z_G2y2z_a-2*I_NAI_G3xz_G2y2z;
  abcd[888] = 2.0E0*I_NAI_I2x3yz_G2y2z_a;
  abcd[889] = 2.0E0*I_NAI_I2x2y2z_G2y2z_a-1*I_NAI_G2x2y_G2y2z;
  abcd[890] = 2.0E0*I_NAI_I2xy3z_G2y2z_a-2*I_NAI_G2xyz_G2y2z;
  abcd[891] = 2.0E0*I_NAI_I2x4z_G2y2z_a-3*I_NAI_G2x2z_G2y2z;
  abcd[892] = 2.0E0*I_NAI_Ix4yz_G2y2z_a;
  abcd[893] = 2.0E0*I_NAI_Ix3y2z_G2y2z_a-1*I_NAI_Gx3y_G2y2z;
  abcd[894] = 2.0E0*I_NAI_Ix2y3z_G2y2z_a-2*I_NAI_Gx2yz_G2y2z;
  abcd[895] = 2.0E0*I_NAI_Ixy4z_G2y2z_a-3*I_NAI_Gxy2z_G2y2z;
  abcd[896] = 2.0E0*I_NAI_Ix5z_G2y2z_a-4*I_NAI_Gx3z_G2y2z;
  abcd[897] = 2.0E0*I_NAI_I5yz_G2y2z_a;
  abcd[898] = 2.0E0*I_NAI_I4y2z_G2y2z_a-1*I_NAI_G4y_G2y2z;
  abcd[899] = 2.0E0*I_NAI_I3y3z_G2y2z_a-2*I_NAI_G3yz_G2y2z;
  abcd[900] = 2.0E0*I_NAI_I2y4z_G2y2z_a-3*I_NAI_G2y2z_G2y2z;
  abcd[901] = 2.0E0*I_NAI_Iy5z_G2y2z_a-4*I_NAI_Gy3z_G2y2z;
  abcd[902] = 2.0E0*I_NAI_I6z_G2y2z_a-5*I_NAI_G4z_G2y2z;
  abcd[903] = 2.0E0*I_NAI_I5xz_Gy3z_a;
  abcd[904] = 2.0E0*I_NAI_I4xyz_Gy3z_a;
  abcd[905] = 2.0E0*I_NAI_I4x2z_Gy3z_a-1*I_NAI_G4x_Gy3z;
  abcd[906] = 2.0E0*I_NAI_I3x2yz_Gy3z_a;
  abcd[907] = 2.0E0*I_NAI_I3xy2z_Gy3z_a-1*I_NAI_G3xy_Gy3z;
  abcd[908] = 2.0E0*I_NAI_I3x3z_Gy3z_a-2*I_NAI_G3xz_Gy3z;
  abcd[909] = 2.0E0*I_NAI_I2x3yz_Gy3z_a;
  abcd[910] = 2.0E0*I_NAI_I2x2y2z_Gy3z_a-1*I_NAI_G2x2y_Gy3z;
  abcd[911] = 2.0E0*I_NAI_I2xy3z_Gy3z_a-2*I_NAI_G2xyz_Gy3z;
  abcd[912] = 2.0E0*I_NAI_I2x4z_Gy3z_a-3*I_NAI_G2x2z_Gy3z;
  abcd[913] = 2.0E0*I_NAI_Ix4yz_Gy3z_a;
  abcd[914] = 2.0E0*I_NAI_Ix3y2z_Gy3z_a-1*I_NAI_Gx3y_Gy3z;
  abcd[915] = 2.0E0*I_NAI_Ix2y3z_Gy3z_a-2*I_NAI_Gx2yz_Gy3z;
  abcd[916] = 2.0E0*I_NAI_Ixy4z_Gy3z_a-3*I_NAI_Gxy2z_Gy3z;
  abcd[917] = 2.0E0*I_NAI_Ix5z_Gy3z_a-4*I_NAI_Gx3z_Gy3z;
  abcd[918] = 2.0E0*I_NAI_I5yz_Gy3z_a;
  abcd[919] = 2.0E0*I_NAI_I4y2z_Gy3z_a-1*I_NAI_G4y_Gy3z;
  abcd[920] = 2.0E0*I_NAI_I3y3z_Gy3z_a-2*I_NAI_G3yz_Gy3z;
  abcd[921] = 2.0E0*I_NAI_I2y4z_Gy3z_a-3*I_NAI_G2y2z_Gy3z;
  abcd[922] = 2.0E0*I_NAI_Iy5z_Gy3z_a-4*I_NAI_Gy3z_Gy3z;
  abcd[923] = 2.0E0*I_NAI_I6z_Gy3z_a-5*I_NAI_G4z_Gy3z;
  abcd[924] = 2.0E0*I_NAI_I5xz_G4z_a;
  abcd[925] = 2.0E0*I_NAI_I4xyz_G4z_a;
  abcd[926] = 2.0E0*I_NAI_I4x2z_G4z_a-1*I_NAI_G4x_G4z;
  abcd[927] = 2.0E0*I_NAI_I3x2yz_G4z_a;
  abcd[928] = 2.0E0*I_NAI_I3xy2z_G4z_a-1*I_NAI_G3xy_G4z;
  abcd[929] = 2.0E0*I_NAI_I3x3z_G4z_a-2*I_NAI_G3xz_G4z;
  abcd[930] = 2.0E0*I_NAI_I2x3yz_G4z_a;
  abcd[931] = 2.0E0*I_NAI_I2x2y2z_G4z_a-1*I_NAI_G2x2y_G4z;
  abcd[932] = 2.0E0*I_NAI_I2xy3z_G4z_a-2*I_NAI_G2xyz_G4z;
  abcd[933] = 2.0E0*I_NAI_I2x4z_G4z_a-3*I_NAI_G2x2z_G4z;
  abcd[934] = 2.0E0*I_NAI_Ix4yz_G4z_a;
  abcd[935] = 2.0E0*I_NAI_Ix3y2z_G4z_a-1*I_NAI_Gx3y_G4z;
  abcd[936] = 2.0E0*I_NAI_Ix2y3z_G4z_a-2*I_NAI_Gx2yz_G4z;
  abcd[937] = 2.0E0*I_NAI_Ixy4z_G4z_a-3*I_NAI_Gxy2z_G4z;
  abcd[938] = 2.0E0*I_NAI_Ix5z_G4z_a-4*I_NAI_Gx3z_G4z;
  abcd[939] = 2.0E0*I_NAI_I5yz_G4z_a;
  abcd[940] = 2.0E0*I_NAI_I4y2z_G4z_a-1*I_NAI_G4y_G4z;
  abcd[941] = 2.0E0*I_NAI_I3y3z_G4z_a-2*I_NAI_G3yz_G4z;
  abcd[942] = 2.0E0*I_NAI_I2y4z_G4z_a-3*I_NAI_G2y2z_G4z;
  abcd[943] = 2.0E0*I_NAI_Iy5z_G4z_a-4*I_NAI_Gy3z_G4z;
  abcd[944] = 2.0E0*I_NAI_I6z_G4z_a-5*I_NAI_G4z_G4z;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_H_b
   * RHS shell quartet name: SQ_NAI_H_F
   ************************************************************/
  abcd[945] = 2.0E0*I_NAI_H5x_H5x_b-4*I_NAI_H5x_F3x;
  abcd[946] = 2.0E0*I_NAI_H4xy_H5x_b-4*I_NAI_H4xy_F3x;
  abcd[947] = 2.0E0*I_NAI_H4xz_H5x_b-4*I_NAI_H4xz_F3x;
  abcd[948] = 2.0E0*I_NAI_H3x2y_H5x_b-4*I_NAI_H3x2y_F3x;
  abcd[949] = 2.0E0*I_NAI_H3xyz_H5x_b-4*I_NAI_H3xyz_F3x;
  abcd[950] = 2.0E0*I_NAI_H3x2z_H5x_b-4*I_NAI_H3x2z_F3x;
  abcd[951] = 2.0E0*I_NAI_H2x3y_H5x_b-4*I_NAI_H2x3y_F3x;
  abcd[952] = 2.0E0*I_NAI_H2x2yz_H5x_b-4*I_NAI_H2x2yz_F3x;
  abcd[953] = 2.0E0*I_NAI_H2xy2z_H5x_b-4*I_NAI_H2xy2z_F3x;
  abcd[954] = 2.0E0*I_NAI_H2x3z_H5x_b-4*I_NAI_H2x3z_F3x;
  abcd[955] = 2.0E0*I_NAI_Hx4y_H5x_b-4*I_NAI_Hx4y_F3x;
  abcd[956] = 2.0E0*I_NAI_Hx3yz_H5x_b-4*I_NAI_Hx3yz_F3x;
  abcd[957] = 2.0E0*I_NAI_Hx2y2z_H5x_b-4*I_NAI_Hx2y2z_F3x;
  abcd[958] = 2.0E0*I_NAI_Hxy3z_H5x_b-4*I_NAI_Hxy3z_F3x;
  abcd[959] = 2.0E0*I_NAI_Hx4z_H5x_b-4*I_NAI_Hx4z_F3x;
  abcd[960] = 2.0E0*I_NAI_H5y_H5x_b-4*I_NAI_H5y_F3x;
  abcd[961] = 2.0E0*I_NAI_H4yz_H5x_b-4*I_NAI_H4yz_F3x;
  abcd[962] = 2.0E0*I_NAI_H3y2z_H5x_b-4*I_NAI_H3y2z_F3x;
  abcd[963] = 2.0E0*I_NAI_H2y3z_H5x_b-4*I_NAI_H2y3z_F3x;
  abcd[964] = 2.0E0*I_NAI_Hy4z_H5x_b-4*I_NAI_Hy4z_F3x;
  abcd[965] = 2.0E0*I_NAI_H5z_H5x_b-4*I_NAI_H5z_F3x;
  abcd[966] = 2.0E0*I_NAI_H5x_H4xy_b-3*I_NAI_H5x_F2xy;
  abcd[967] = 2.0E0*I_NAI_H4xy_H4xy_b-3*I_NAI_H4xy_F2xy;
  abcd[968] = 2.0E0*I_NAI_H4xz_H4xy_b-3*I_NAI_H4xz_F2xy;
  abcd[969] = 2.0E0*I_NAI_H3x2y_H4xy_b-3*I_NAI_H3x2y_F2xy;
  abcd[970] = 2.0E0*I_NAI_H3xyz_H4xy_b-3*I_NAI_H3xyz_F2xy;
  abcd[971] = 2.0E0*I_NAI_H3x2z_H4xy_b-3*I_NAI_H3x2z_F2xy;
  abcd[972] = 2.0E0*I_NAI_H2x3y_H4xy_b-3*I_NAI_H2x3y_F2xy;
  abcd[973] = 2.0E0*I_NAI_H2x2yz_H4xy_b-3*I_NAI_H2x2yz_F2xy;
  abcd[974] = 2.0E0*I_NAI_H2xy2z_H4xy_b-3*I_NAI_H2xy2z_F2xy;
  abcd[975] = 2.0E0*I_NAI_H2x3z_H4xy_b-3*I_NAI_H2x3z_F2xy;
  abcd[976] = 2.0E0*I_NAI_Hx4y_H4xy_b-3*I_NAI_Hx4y_F2xy;
  abcd[977] = 2.0E0*I_NAI_Hx3yz_H4xy_b-3*I_NAI_Hx3yz_F2xy;
  abcd[978] = 2.0E0*I_NAI_Hx2y2z_H4xy_b-3*I_NAI_Hx2y2z_F2xy;
  abcd[979] = 2.0E0*I_NAI_Hxy3z_H4xy_b-3*I_NAI_Hxy3z_F2xy;
  abcd[980] = 2.0E0*I_NAI_Hx4z_H4xy_b-3*I_NAI_Hx4z_F2xy;
  abcd[981] = 2.0E0*I_NAI_H5y_H4xy_b-3*I_NAI_H5y_F2xy;
  abcd[982] = 2.0E0*I_NAI_H4yz_H4xy_b-3*I_NAI_H4yz_F2xy;
  abcd[983] = 2.0E0*I_NAI_H3y2z_H4xy_b-3*I_NAI_H3y2z_F2xy;
  abcd[984] = 2.0E0*I_NAI_H2y3z_H4xy_b-3*I_NAI_H2y3z_F2xy;
  abcd[985] = 2.0E0*I_NAI_Hy4z_H4xy_b-3*I_NAI_Hy4z_F2xy;
  abcd[986] = 2.0E0*I_NAI_H5z_H4xy_b-3*I_NAI_H5z_F2xy;
  abcd[987] = 2.0E0*I_NAI_H5x_H4xz_b-3*I_NAI_H5x_F2xz;
  abcd[988] = 2.0E0*I_NAI_H4xy_H4xz_b-3*I_NAI_H4xy_F2xz;
  abcd[989] = 2.0E0*I_NAI_H4xz_H4xz_b-3*I_NAI_H4xz_F2xz;
  abcd[990] = 2.0E0*I_NAI_H3x2y_H4xz_b-3*I_NAI_H3x2y_F2xz;
  abcd[991] = 2.0E0*I_NAI_H3xyz_H4xz_b-3*I_NAI_H3xyz_F2xz;
  abcd[992] = 2.0E0*I_NAI_H3x2z_H4xz_b-3*I_NAI_H3x2z_F2xz;
  abcd[993] = 2.0E0*I_NAI_H2x3y_H4xz_b-3*I_NAI_H2x3y_F2xz;
  abcd[994] = 2.0E0*I_NAI_H2x2yz_H4xz_b-3*I_NAI_H2x2yz_F2xz;
  abcd[995] = 2.0E0*I_NAI_H2xy2z_H4xz_b-3*I_NAI_H2xy2z_F2xz;
  abcd[996] = 2.0E0*I_NAI_H2x3z_H4xz_b-3*I_NAI_H2x3z_F2xz;
  abcd[997] = 2.0E0*I_NAI_Hx4y_H4xz_b-3*I_NAI_Hx4y_F2xz;
  abcd[998] = 2.0E0*I_NAI_Hx3yz_H4xz_b-3*I_NAI_Hx3yz_F2xz;
  abcd[999] = 2.0E0*I_NAI_Hx2y2z_H4xz_b-3*I_NAI_Hx2y2z_F2xz;
  abcd[1000] = 2.0E0*I_NAI_Hxy3z_H4xz_b-3*I_NAI_Hxy3z_F2xz;
  abcd[1001] = 2.0E0*I_NAI_Hx4z_H4xz_b-3*I_NAI_Hx4z_F2xz;
  abcd[1002] = 2.0E0*I_NAI_H5y_H4xz_b-3*I_NAI_H5y_F2xz;
  abcd[1003] = 2.0E0*I_NAI_H4yz_H4xz_b-3*I_NAI_H4yz_F2xz;
  abcd[1004] = 2.0E0*I_NAI_H3y2z_H4xz_b-3*I_NAI_H3y2z_F2xz;
  abcd[1005] = 2.0E0*I_NAI_H2y3z_H4xz_b-3*I_NAI_H2y3z_F2xz;
  abcd[1006] = 2.0E0*I_NAI_Hy4z_H4xz_b-3*I_NAI_Hy4z_F2xz;
  abcd[1007] = 2.0E0*I_NAI_H5z_H4xz_b-3*I_NAI_H5z_F2xz;
  abcd[1008] = 2.0E0*I_NAI_H5x_H3x2y_b-2*I_NAI_H5x_Fx2y;
  abcd[1009] = 2.0E0*I_NAI_H4xy_H3x2y_b-2*I_NAI_H4xy_Fx2y;
  abcd[1010] = 2.0E0*I_NAI_H4xz_H3x2y_b-2*I_NAI_H4xz_Fx2y;
  abcd[1011] = 2.0E0*I_NAI_H3x2y_H3x2y_b-2*I_NAI_H3x2y_Fx2y;
  abcd[1012] = 2.0E0*I_NAI_H3xyz_H3x2y_b-2*I_NAI_H3xyz_Fx2y;
  abcd[1013] = 2.0E0*I_NAI_H3x2z_H3x2y_b-2*I_NAI_H3x2z_Fx2y;
  abcd[1014] = 2.0E0*I_NAI_H2x3y_H3x2y_b-2*I_NAI_H2x3y_Fx2y;
  abcd[1015] = 2.0E0*I_NAI_H2x2yz_H3x2y_b-2*I_NAI_H2x2yz_Fx2y;
  abcd[1016] = 2.0E0*I_NAI_H2xy2z_H3x2y_b-2*I_NAI_H2xy2z_Fx2y;
  abcd[1017] = 2.0E0*I_NAI_H2x3z_H3x2y_b-2*I_NAI_H2x3z_Fx2y;
  abcd[1018] = 2.0E0*I_NAI_Hx4y_H3x2y_b-2*I_NAI_Hx4y_Fx2y;
  abcd[1019] = 2.0E0*I_NAI_Hx3yz_H3x2y_b-2*I_NAI_Hx3yz_Fx2y;
  abcd[1020] = 2.0E0*I_NAI_Hx2y2z_H3x2y_b-2*I_NAI_Hx2y2z_Fx2y;
  abcd[1021] = 2.0E0*I_NAI_Hxy3z_H3x2y_b-2*I_NAI_Hxy3z_Fx2y;
  abcd[1022] = 2.0E0*I_NAI_Hx4z_H3x2y_b-2*I_NAI_Hx4z_Fx2y;
  abcd[1023] = 2.0E0*I_NAI_H5y_H3x2y_b-2*I_NAI_H5y_Fx2y;
  abcd[1024] = 2.0E0*I_NAI_H4yz_H3x2y_b-2*I_NAI_H4yz_Fx2y;
  abcd[1025] = 2.0E0*I_NAI_H3y2z_H3x2y_b-2*I_NAI_H3y2z_Fx2y;
  abcd[1026] = 2.0E0*I_NAI_H2y3z_H3x2y_b-2*I_NAI_H2y3z_Fx2y;
  abcd[1027] = 2.0E0*I_NAI_Hy4z_H3x2y_b-2*I_NAI_Hy4z_Fx2y;
  abcd[1028] = 2.0E0*I_NAI_H5z_H3x2y_b-2*I_NAI_H5z_Fx2y;
  abcd[1029] = 2.0E0*I_NAI_H5x_H3xyz_b-2*I_NAI_H5x_Fxyz;
  abcd[1030] = 2.0E0*I_NAI_H4xy_H3xyz_b-2*I_NAI_H4xy_Fxyz;
  abcd[1031] = 2.0E0*I_NAI_H4xz_H3xyz_b-2*I_NAI_H4xz_Fxyz;
  abcd[1032] = 2.0E0*I_NAI_H3x2y_H3xyz_b-2*I_NAI_H3x2y_Fxyz;
  abcd[1033] = 2.0E0*I_NAI_H3xyz_H3xyz_b-2*I_NAI_H3xyz_Fxyz;
  abcd[1034] = 2.0E0*I_NAI_H3x2z_H3xyz_b-2*I_NAI_H3x2z_Fxyz;
  abcd[1035] = 2.0E0*I_NAI_H2x3y_H3xyz_b-2*I_NAI_H2x3y_Fxyz;
  abcd[1036] = 2.0E0*I_NAI_H2x2yz_H3xyz_b-2*I_NAI_H2x2yz_Fxyz;
  abcd[1037] = 2.0E0*I_NAI_H2xy2z_H3xyz_b-2*I_NAI_H2xy2z_Fxyz;
  abcd[1038] = 2.0E0*I_NAI_H2x3z_H3xyz_b-2*I_NAI_H2x3z_Fxyz;
  abcd[1039] = 2.0E0*I_NAI_Hx4y_H3xyz_b-2*I_NAI_Hx4y_Fxyz;
  abcd[1040] = 2.0E0*I_NAI_Hx3yz_H3xyz_b-2*I_NAI_Hx3yz_Fxyz;
  abcd[1041] = 2.0E0*I_NAI_Hx2y2z_H3xyz_b-2*I_NAI_Hx2y2z_Fxyz;
  abcd[1042] = 2.0E0*I_NAI_Hxy3z_H3xyz_b-2*I_NAI_Hxy3z_Fxyz;
  abcd[1043] = 2.0E0*I_NAI_Hx4z_H3xyz_b-2*I_NAI_Hx4z_Fxyz;
  abcd[1044] = 2.0E0*I_NAI_H5y_H3xyz_b-2*I_NAI_H5y_Fxyz;
  abcd[1045] = 2.0E0*I_NAI_H4yz_H3xyz_b-2*I_NAI_H4yz_Fxyz;
  abcd[1046] = 2.0E0*I_NAI_H3y2z_H3xyz_b-2*I_NAI_H3y2z_Fxyz;
  abcd[1047] = 2.0E0*I_NAI_H2y3z_H3xyz_b-2*I_NAI_H2y3z_Fxyz;
  abcd[1048] = 2.0E0*I_NAI_Hy4z_H3xyz_b-2*I_NAI_Hy4z_Fxyz;
  abcd[1049] = 2.0E0*I_NAI_H5z_H3xyz_b-2*I_NAI_H5z_Fxyz;
  abcd[1050] = 2.0E0*I_NAI_H5x_H3x2z_b-2*I_NAI_H5x_Fx2z;
  abcd[1051] = 2.0E0*I_NAI_H4xy_H3x2z_b-2*I_NAI_H4xy_Fx2z;
  abcd[1052] = 2.0E0*I_NAI_H4xz_H3x2z_b-2*I_NAI_H4xz_Fx2z;
  abcd[1053] = 2.0E0*I_NAI_H3x2y_H3x2z_b-2*I_NAI_H3x2y_Fx2z;
  abcd[1054] = 2.0E0*I_NAI_H3xyz_H3x2z_b-2*I_NAI_H3xyz_Fx2z;
  abcd[1055] = 2.0E0*I_NAI_H3x2z_H3x2z_b-2*I_NAI_H3x2z_Fx2z;
  abcd[1056] = 2.0E0*I_NAI_H2x3y_H3x2z_b-2*I_NAI_H2x3y_Fx2z;
  abcd[1057] = 2.0E0*I_NAI_H2x2yz_H3x2z_b-2*I_NAI_H2x2yz_Fx2z;
  abcd[1058] = 2.0E0*I_NAI_H2xy2z_H3x2z_b-2*I_NAI_H2xy2z_Fx2z;
  abcd[1059] = 2.0E0*I_NAI_H2x3z_H3x2z_b-2*I_NAI_H2x3z_Fx2z;
  abcd[1060] = 2.0E0*I_NAI_Hx4y_H3x2z_b-2*I_NAI_Hx4y_Fx2z;
  abcd[1061] = 2.0E0*I_NAI_Hx3yz_H3x2z_b-2*I_NAI_Hx3yz_Fx2z;
  abcd[1062] = 2.0E0*I_NAI_Hx2y2z_H3x2z_b-2*I_NAI_Hx2y2z_Fx2z;
  abcd[1063] = 2.0E0*I_NAI_Hxy3z_H3x2z_b-2*I_NAI_Hxy3z_Fx2z;
  abcd[1064] = 2.0E0*I_NAI_Hx4z_H3x2z_b-2*I_NAI_Hx4z_Fx2z;
  abcd[1065] = 2.0E0*I_NAI_H5y_H3x2z_b-2*I_NAI_H5y_Fx2z;
  abcd[1066] = 2.0E0*I_NAI_H4yz_H3x2z_b-2*I_NAI_H4yz_Fx2z;
  abcd[1067] = 2.0E0*I_NAI_H3y2z_H3x2z_b-2*I_NAI_H3y2z_Fx2z;
  abcd[1068] = 2.0E0*I_NAI_H2y3z_H3x2z_b-2*I_NAI_H2y3z_Fx2z;
  abcd[1069] = 2.0E0*I_NAI_Hy4z_H3x2z_b-2*I_NAI_Hy4z_Fx2z;
  abcd[1070] = 2.0E0*I_NAI_H5z_H3x2z_b-2*I_NAI_H5z_Fx2z;
  abcd[1071] = 2.0E0*I_NAI_H5x_H2x3y_b-1*I_NAI_H5x_F3y;
  abcd[1072] = 2.0E0*I_NAI_H4xy_H2x3y_b-1*I_NAI_H4xy_F3y;
  abcd[1073] = 2.0E0*I_NAI_H4xz_H2x3y_b-1*I_NAI_H4xz_F3y;
  abcd[1074] = 2.0E0*I_NAI_H3x2y_H2x3y_b-1*I_NAI_H3x2y_F3y;
  abcd[1075] = 2.0E0*I_NAI_H3xyz_H2x3y_b-1*I_NAI_H3xyz_F3y;
  abcd[1076] = 2.0E0*I_NAI_H3x2z_H2x3y_b-1*I_NAI_H3x2z_F3y;
  abcd[1077] = 2.0E0*I_NAI_H2x3y_H2x3y_b-1*I_NAI_H2x3y_F3y;
  abcd[1078] = 2.0E0*I_NAI_H2x2yz_H2x3y_b-1*I_NAI_H2x2yz_F3y;
  abcd[1079] = 2.0E0*I_NAI_H2xy2z_H2x3y_b-1*I_NAI_H2xy2z_F3y;
  abcd[1080] = 2.0E0*I_NAI_H2x3z_H2x3y_b-1*I_NAI_H2x3z_F3y;
  abcd[1081] = 2.0E0*I_NAI_Hx4y_H2x3y_b-1*I_NAI_Hx4y_F3y;
  abcd[1082] = 2.0E0*I_NAI_Hx3yz_H2x3y_b-1*I_NAI_Hx3yz_F3y;
  abcd[1083] = 2.0E0*I_NAI_Hx2y2z_H2x3y_b-1*I_NAI_Hx2y2z_F3y;
  abcd[1084] = 2.0E0*I_NAI_Hxy3z_H2x3y_b-1*I_NAI_Hxy3z_F3y;
  abcd[1085] = 2.0E0*I_NAI_Hx4z_H2x3y_b-1*I_NAI_Hx4z_F3y;
  abcd[1086] = 2.0E0*I_NAI_H5y_H2x3y_b-1*I_NAI_H5y_F3y;
  abcd[1087] = 2.0E0*I_NAI_H4yz_H2x3y_b-1*I_NAI_H4yz_F3y;
  abcd[1088] = 2.0E0*I_NAI_H3y2z_H2x3y_b-1*I_NAI_H3y2z_F3y;
  abcd[1089] = 2.0E0*I_NAI_H2y3z_H2x3y_b-1*I_NAI_H2y3z_F3y;
  abcd[1090] = 2.0E0*I_NAI_Hy4z_H2x3y_b-1*I_NAI_Hy4z_F3y;
  abcd[1091] = 2.0E0*I_NAI_H5z_H2x3y_b-1*I_NAI_H5z_F3y;
  abcd[1092] = 2.0E0*I_NAI_H5x_H2x2yz_b-1*I_NAI_H5x_F2yz;
  abcd[1093] = 2.0E0*I_NAI_H4xy_H2x2yz_b-1*I_NAI_H4xy_F2yz;
  abcd[1094] = 2.0E0*I_NAI_H4xz_H2x2yz_b-1*I_NAI_H4xz_F2yz;
  abcd[1095] = 2.0E0*I_NAI_H3x2y_H2x2yz_b-1*I_NAI_H3x2y_F2yz;
  abcd[1096] = 2.0E0*I_NAI_H3xyz_H2x2yz_b-1*I_NAI_H3xyz_F2yz;
  abcd[1097] = 2.0E0*I_NAI_H3x2z_H2x2yz_b-1*I_NAI_H3x2z_F2yz;
  abcd[1098] = 2.0E0*I_NAI_H2x3y_H2x2yz_b-1*I_NAI_H2x3y_F2yz;
  abcd[1099] = 2.0E0*I_NAI_H2x2yz_H2x2yz_b-1*I_NAI_H2x2yz_F2yz;
  abcd[1100] = 2.0E0*I_NAI_H2xy2z_H2x2yz_b-1*I_NAI_H2xy2z_F2yz;
  abcd[1101] = 2.0E0*I_NAI_H2x3z_H2x2yz_b-1*I_NAI_H2x3z_F2yz;
  abcd[1102] = 2.0E0*I_NAI_Hx4y_H2x2yz_b-1*I_NAI_Hx4y_F2yz;
  abcd[1103] = 2.0E0*I_NAI_Hx3yz_H2x2yz_b-1*I_NAI_Hx3yz_F2yz;
  abcd[1104] = 2.0E0*I_NAI_Hx2y2z_H2x2yz_b-1*I_NAI_Hx2y2z_F2yz;
  abcd[1105] = 2.0E0*I_NAI_Hxy3z_H2x2yz_b-1*I_NAI_Hxy3z_F2yz;
  abcd[1106] = 2.0E0*I_NAI_Hx4z_H2x2yz_b-1*I_NAI_Hx4z_F2yz;
  abcd[1107] = 2.0E0*I_NAI_H5y_H2x2yz_b-1*I_NAI_H5y_F2yz;
  abcd[1108] = 2.0E0*I_NAI_H4yz_H2x2yz_b-1*I_NAI_H4yz_F2yz;
  abcd[1109] = 2.0E0*I_NAI_H3y2z_H2x2yz_b-1*I_NAI_H3y2z_F2yz;
  abcd[1110] = 2.0E0*I_NAI_H2y3z_H2x2yz_b-1*I_NAI_H2y3z_F2yz;
  abcd[1111] = 2.0E0*I_NAI_Hy4z_H2x2yz_b-1*I_NAI_Hy4z_F2yz;
  abcd[1112] = 2.0E0*I_NAI_H5z_H2x2yz_b-1*I_NAI_H5z_F2yz;
  abcd[1113] = 2.0E0*I_NAI_H5x_H2xy2z_b-1*I_NAI_H5x_Fy2z;
  abcd[1114] = 2.0E0*I_NAI_H4xy_H2xy2z_b-1*I_NAI_H4xy_Fy2z;
  abcd[1115] = 2.0E0*I_NAI_H4xz_H2xy2z_b-1*I_NAI_H4xz_Fy2z;
  abcd[1116] = 2.0E0*I_NAI_H3x2y_H2xy2z_b-1*I_NAI_H3x2y_Fy2z;
  abcd[1117] = 2.0E0*I_NAI_H3xyz_H2xy2z_b-1*I_NAI_H3xyz_Fy2z;
  abcd[1118] = 2.0E0*I_NAI_H3x2z_H2xy2z_b-1*I_NAI_H3x2z_Fy2z;
  abcd[1119] = 2.0E0*I_NAI_H2x3y_H2xy2z_b-1*I_NAI_H2x3y_Fy2z;
  abcd[1120] = 2.0E0*I_NAI_H2x2yz_H2xy2z_b-1*I_NAI_H2x2yz_Fy2z;
  abcd[1121] = 2.0E0*I_NAI_H2xy2z_H2xy2z_b-1*I_NAI_H2xy2z_Fy2z;
  abcd[1122] = 2.0E0*I_NAI_H2x3z_H2xy2z_b-1*I_NAI_H2x3z_Fy2z;
  abcd[1123] = 2.0E0*I_NAI_Hx4y_H2xy2z_b-1*I_NAI_Hx4y_Fy2z;
  abcd[1124] = 2.0E0*I_NAI_Hx3yz_H2xy2z_b-1*I_NAI_Hx3yz_Fy2z;
  abcd[1125] = 2.0E0*I_NAI_Hx2y2z_H2xy2z_b-1*I_NAI_Hx2y2z_Fy2z;
  abcd[1126] = 2.0E0*I_NAI_Hxy3z_H2xy2z_b-1*I_NAI_Hxy3z_Fy2z;
  abcd[1127] = 2.0E0*I_NAI_Hx4z_H2xy2z_b-1*I_NAI_Hx4z_Fy2z;
  abcd[1128] = 2.0E0*I_NAI_H5y_H2xy2z_b-1*I_NAI_H5y_Fy2z;
  abcd[1129] = 2.0E0*I_NAI_H4yz_H2xy2z_b-1*I_NAI_H4yz_Fy2z;
  abcd[1130] = 2.0E0*I_NAI_H3y2z_H2xy2z_b-1*I_NAI_H3y2z_Fy2z;
  abcd[1131] = 2.0E0*I_NAI_H2y3z_H2xy2z_b-1*I_NAI_H2y3z_Fy2z;
  abcd[1132] = 2.0E0*I_NAI_Hy4z_H2xy2z_b-1*I_NAI_Hy4z_Fy2z;
  abcd[1133] = 2.0E0*I_NAI_H5z_H2xy2z_b-1*I_NAI_H5z_Fy2z;
  abcd[1134] = 2.0E0*I_NAI_H5x_H2x3z_b-1*I_NAI_H5x_F3z;
  abcd[1135] = 2.0E0*I_NAI_H4xy_H2x3z_b-1*I_NAI_H4xy_F3z;
  abcd[1136] = 2.0E0*I_NAI_H4xz_H2x3z_b-1*I_NAI_H4xz_F3z;
  abcd[1137] = 2.0E0*I_NAI_H3x2y_H2x3z_b-1*I_NAI_H3x2y_F3z;
  abcd[1138] = 2.0E0*I_NAI_H3xyz_H2x3z_b-1*I_NAI_H3xyz_F3z;
  abcd[1139] = 2.0E0*I_NAI_H3x2z_H2x3z_b-1*I_NAI_H3x2z_F3z;
  abcd[1140] = 2.0E0*I_NAI_H2x3y_H2x3z_b-1*I_NAI_H2x3y_F3z;
  abcd[1141] = 2.0E0*I_NAI_H2x2yz_H2x3z_b-1*I_NAI_H2x2yz_F3z;
  abcd[1142] = 2.0E0*I_NAI_H2xy2z_H2x3z_b-1*I_NAI_H2xy2z_F3z;
  abcd[1143] = 2.0E0*I_NAI_H2x3z_H2x3z_b-1*I_NAI_H2x3z_F3z;
  abcd[1144] = 2.0E0*I_NAI_Hx4y_H2x3z_b-1*I_NAI_Hx4y_F3z;
  abcd[1145] = 2.0E0*I_NAI_Hx3yz_H2x3z_b-1*I_NAI_Hx3yz_F3z;
  abcd[1146] = 2.0E0*I_NAI_Hx2y2z_H2x3z_b-1*I_NAI_Hx2y2z_F3z;
  abcd[1147] = 2.0E0*I_NAI_Hxy3z_H2x3z_b-1*I_NAI_Hxy3z_F3z;
  abcd[1148] = 2.0E0*I_NAI_Hx4z_H2x3z_b-1*I_NAI_Hx4z_F3z;
  abcd[1149] = 2.0E0*I_NAI_H5y_H2x3z_b-1*I_NAI_H5y_F3z;
  abcd[1150] = 2.0E0*I_NAI_H4yz_H2x3z_b-1*I_NAI_H4yz_F3z;
  abcd[1151] = 2.0E0*I_NAI_H3y2z_H2x3z_b-1*I_NAI_H3y2z_F3z;
  abcd[1152] = 2.0E0*I_NAI_H2y3z_H2x3z_b-1*I_NAI_H2y3z_F3z;
  abcd[1153] = 2.0E0*I_NAI_Hy4z_H2x3z_b-1*I_NAI_Hy4z_F3z;
  abcd[1154] = 2.0E0*I_NAI_H5z_H2x3z_b-1*I_NAI_H5z_F3z;
  abcd[1155] = 2.0E0*I_NAI_H5x_Hx4y_b;
  abcd[1156] = 2.0E0*I_NAI_H4xy_Hx4y_b;
  abcd[1157] = 2.0E0*I_NAI_H4xz_Hx4y_b;
  abcd[1158] = 2.0E0*I_NAI_H3x2y_Hx4y_b;
  abcd[1159] = 2.0E0*I_NAI_H3xyz_Hx4y_b;
  abcd[1160] = 2.0E0*I_NAI_H3x2z_Hx4y_b;
  abcd[1161] = 2.0E0*I_NAI_H2x3y_Hx4y_b;
  abcd[1162] = 2.0E0*I_NAI_H2x2yz_Hx4y_b;
  abcd[1163] = 2.0E0*I_NAI_H2xy2z_Hx4y_b;
  abcd[1164] = 2.0E0*I_NAI_H2x3z_Hx4y_b;
  abcd[1165] = 2.0E0*I_NAI_Hx4y_Hx4y_b;
  abcd[1166] = 2.0E0*I_NAI_Hx3yz_Hx4y_b;
  abcd[1167] = 2.0E0*I_NAI_Hx2y2z_Hx4y_b;
  abcd[1168] = 2.0E0*I_NAI_Hxy3z_Hx4y_b;
  abcd[1169] = 2.0E0*I_NAI_Hx4z_Hx4y_b;
  abcd[1170] = 2.0E0*I_NAI_H5y_Hx4y_b;
  abcd[1171] = 2.0E0*I_NAI_H4yz_Hx4y_b;
  abcd[1172] = 2.0E0*I_NAI_H3y2z_Hx4y_b;
  abcd[1173] = 2.0E0*I_NAI_H2y3z_Hx4y_b;
  abcd[1174] = 2.0E0*I_NAI_Hy4z_Hx4y_b;
  abcd[1175] = 2.0E0*I_NAI_H5z_Hx4y_b;
  abcd[1176] = 2.0E0*I_NAI_H5x_Hx3yz_b;
  abcd[1177] = 2.0E0*I_NAI_H4xy_Hx3yz_b;
  abcd[1178] = 2.0E0*I_NAI_H4xz_Hx3yz_b;
  abcd[1179] = 2.0E0*I_NAI_H3x2y_Hx3yz_b;
  abcd[1180] = 2.0E0*I_NAI_H3xyz_Hx3yz_b;
  abcd[1181] = 2.0E0*I_NAI_H3x2z_Hx3yz_b;
  abcd[1182] = 2.0E0*I_NAI_H2x3y_Hx3yz_b;
  abcd[1183] = 2.0E0*I_NAI_H2x2yz_Hx3yz_b;
  abcd[1184] = 2.0E0*I_NAI_H2xy2z_Hx3yz_b;
  abcd[1185] = 2.0E0*I_NAI_H2x3z_Hx3yz_b;
  abcd[1186] = 2.0E0*I_NAI_Hx4y_Hx3yz_b;
  abcd[1187] = 2.0E0*I_NAI_Hx3yz_Hx3yz_b;
  abcd[1188] = 2.0E0*I_NAI_Hx2y2z_Hx3yz_b;
  abcd[1189] = 2.0E0*I_NAI_Hxy3z_Hx3yz_b;
  abcd[1190] = 2.0E0*I_NAI_Hx4z_Hx3yz_b;
  abcd[1191] = 2.0E0*I_NAI_H5y_Hx3yz_b;
  abcd[1192] = 2.0E0*I_NAI_H4yz_Hx3yz_b;
  abcd[1193] = 2.0E0*I_NAI_H3y2z_Hx3yz_b;
  abcd[1194] = 2.0E0*I_NAI_H2y3z_Hx3yz_b;
  abcd[1195] = 2.0E0*I_NAI_Hy4z_Hx3yz_b;
  abcd[1196] = 2.0E0*I_NAI_H5z_Hx3yz_b;
  abcd[1197] = 2.0E0*I_NAI_H5x_Hx2y2z_b;
  abcd[1198] = 2.0E0*I_NAI_H4xy_Hx2y2z_b;
  abcd[1199] = 2.0E0*I_NAI_H4xz_Hx2y2z_b;
  abcd[1200] = 2.0E0*I_NAI_H3x2y_Hx2y2z_b;
  abcd[1201] = 2.0E0*I_NAI_H3xyz_Hx2y2z_b;
  abcd[1202] = 2.0E0*I_NAI_H3x2z_Hx2y2z_b;
  abcd[1203] = 2.0E0*I_NAI_H2x3y_Hx2y2z_b;
  abcd[1204] = 2.0E0*I_NAI_H2x2yz_Hx2y2z_b;
  abcd[1205] = 2.0E0*I_NAI_H2xy2z_Hx2y2z_b;
  abcd[1206] = 2.0E0*I_NAI_H2x3z_Hx2y2z_b;
  abcd[1207] = 2.0E0*I_NAI_Hx4y_Hx2y2z_b;
  abcd[1208] = 2.0E0*I_NAI_Hx3yz_Hx2y2z_b;
  abcd[1209] = 2.0E0*I_NAI_Hx2y2z_Hx2y2z_b;
  abcd[1210] = 2.0E0*I_NAI_Hxy3z_Hx2y2z_b;
  abcd[1211] = 2.0E0*I_NAI_Hx4z_Hx2y2z_b;
  abcd[1212] = 2.0E0*I_NAI_H5y_Hx2y2z_b;
  abcd[1213] = 2.0E0*I_NAI_H4yz_Hx2y2z_b;
  abcd[1214] = 2.0E0*I_NAI_H3y2z_Hx2y2z_b;
  abcd[1215] = 2.0E0*I_NAI_H2y3z_Hx2y2z_b;
  abcd[1216] = 2.0E0*I_NAI_Hy4z_Hx2y2z_b;
  abcd[1217] = 2.0E0*I_NAI_H5z_Hx2y2z_b;
  abcd[1218] = 2.0E0*I_NAI_H5x_Hxy3z_b;
  abcd[1219] = 2.0E0*I_NAI_H4xy_Hxy3z_b;
  abcd[1220] = 2.0E0*I_NAI_H4xz_Hxy3z_b;
  abcd[1221] = 2.0E0*I_NAI_H3x2y_Hxy3z_b;
  abcd[1222] = 2.0E0*I_NAI_H3xyz_Hxy3z_b;
  abcd[1223] = 2.0E0*I_NAI_H3x2z_Hxy3z_b;
  abcd[1224] = 2.0E0*I_NAI_H2x3y_Hxy3z_b;
  abcd[1225] = 2.0E0*I_NAI_H2x2yz_Hxy3z_b;
  abcd[1226] = 2.0E0*I_NAI_H2xy2z_Hxy3z_b;
  abcd[1227] = 2.0E0*I_NAI_H2x3z_Hxy3z_b;
  abcd[1228] = 2.0E0*I_NAI_Hx4y_Hxy3z_b;
  abcd[1229] = 2.0E0*I_NAI_Hx3yz_Hxy3z_b;
  abcd[1230] = 2.0E0*I_NAI_Hx2y2z_Hxy3z_b;
  abcd[1231] = 2.0E0*I_NAI_Hxy3z_Hxy3z_b;
  abcd[1232] = 2.0E0*I_NAI_Hx4z_Hxy3z_b;
  abcd[1233] = 2.0E0*I_NAI_H5y_Hxy3z_b;
  abcd[1234] = 2.0E0*I_NAI_H4yz_Hxy3z_b;
  abcd[1235] = 2.0E0*I_NAI_H3y2z_Hxy3z_b;
  abcd[1236] = 2.0E0*I_NAI_H2y3z_Hxy3z_b;
  abcd[1237] = 2.0E0*I_NAI_Hy4z_Hxy3z_b;
  abcd[1238] = 2.0E0*I_NAI_H5z_Hxy3z_b;
  abcd[1239] = 2.0E0*I_NAI_H5x_Hx4z_b;
  abcd[1240] = 2.0E0*I_NAI_H4xy_Hx4z_b;
  abcd[1241] = 2.0E0*I_NAI_H4xz_Hx4z_b;
  abcd[1242] = 2.0E0*I_NAI_H3x2y_Hx4z_b;
  abcd[1243] = 2.0E0*I_NAI_H3xyz_Hx4z_b;
  abcd[1244] = 2.0E0*I_NAI_H3x2z_Hx4z_b;
  abcd[1245] = 2.0E0*I_NAI_H2x3y_Hx4z_b;
  abcd[1246] = 2.0E0*I_NAI_H2x2yz_Hx4z_b;
  abcd[1247] = 2.0E0*I_NAI_H2xy2z_Hx4z_b;
  abcd[1248] = 2.0E0*I_NAI_H2x3z_Hx4z_b;
  abcd[1249] = 2.0E0*I_NAI_Hx4y_Hx4z_b;
  abcd[1250] = 2.0E0*I_NAI_Hx3yz_Hx4z_b;
  abcd[1251] = 2.0E0*I_NAI_Hx2y2z_Hx4z_b;
  abcd[1252] = 2.0E0*I_NAI_Hxy3z_Hx4z_b;
  abcd[1253] = 2.0E0*I_NAI_Hx4z_Hx4z_b;
  abcd[1254] = 2.0E0*I_NAI_H5y_Hx4z_b;
  abcd[1255] = 2.0E0*I_NAI_H4yz_Hx4z_b;
  abcd[1256] = 2.0E0*I_NAI_H3y2z_Hx4z_b;
  abcd[1257] = 2.0E0*I_NAI_H2y3z_Hx4z_b;
  abcd[1258] = 2.0E0*I_NAI_Hy4z_Hx4z_b;
  abcd[1259] = 2.0E0*I_NAI_H5z_Hx4z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_H_b
   * RHS shell quartet name: SQ_NAI_H_F
   ************************************************************/
  abcd[1260] = 2.0E0*I_NAI_H5x_H4xy_b;
  abcd[1261] = 2.0E0*I_NAI_H4xy_H4xy_b;
  abcd[1262] = 2.0E0*I_NAI_H4xz_H4xy_b;
  abcd[1263] = 2.0E0*I_NAI_H3x2y_H4xy_b;
  abcd[1264] = 2.0E0*I_NAI_H3xyz_H4xy_b;
  abcd[1265] = 2.0E0*I_NAI_H3x2z_H4xy_b;
  abcd[1266] = 2.0E0*I_NAI_H2x3y_H4xy_b;
  abcd[1267] = 2.0E0*I_NAI_H2x2yz_H4xy_b;
  abcd[1268] = 2.0E0*I_NAI_H2xy2z_H4xy_b;
  abcd[1269] = 2.0E0*I_NAI_H2x3z_H4xy_b;
  abcd[1270] = 2.0E0*I_NAI_Hx4y_H4xy_b;
  abcd[1271] = 2.0E0*I_NAI_Hx3yz_H4xy_b;
  abcd[1272] = 2.0E0*I_NAI_Hx2y2z_H4xy_b;
  abcd[1273] = 2.0E0*I_NAI_Hxy3z_H4xy_b;
  abcd[1274] = 2.0E0*I_NAI_Hx4z_H4xy_b;
  abcd[1275] = 2.0E0*I_NAI_H5y_H4xy_b;
  abcd[1276] = 2.0E0*I_NAI_H4yz_H4xy_b;
  abcd[1277] = 2.0E0*I_NAI_H3y2z_H4xy_b;
  abcd[1278] = 2.0E0*I_NAI_H2y3z_H4xy_b;
  abcd[1279] = 2.0E0*I_NAI_Hy4z_H4xy_b;
  abcd[1280] = 2.0E0*I_NAI_H5z_H4xy_b;
  abcd[1281] = 2.0E0*I_NAI_H5x_H3x2y_b-1*I_NAI_H5x_F3x;
  abcd[1282] = 2.0E0*I_NAI_H4xy_H3x2y_b-1*I_NAI_H4xy_F3x;
  abcd[1283] = 2.0E0*I_NAI_H4xz_H3x2y_b-1*I_NAI_H4xz_F3x;
  abcd[1284] = 2.0E0*I_NAI_H3x2y_H3x2y_b-1*I_NAI_H3x2y_F3x;
  abcd[1285] = 2.0E0*I_NAI_H3xyz_H3x2y_b-1*I_NAI_H3xyz_F3x;
  abcd[1286] = 2.0E0*I_NAI_H3x2z_H3x2y_b-1*I_NAI_H3x2z_F3x;
  abcd[1287] = 2.0E0*I_NAI_H2x3y_H3x2y_b-1*I_NAI_H2x3y_F3x;
  abcd[1288] = 2.0E0*I_NAI_H2x2yz_H3x2y_b-1*I_NAI_H2x2yz_F3x;
  abcd[1289] = 2.0E0*I_NAI_H2xy2z_H3x2y_b-1*I_NAI_H2xy2z_F3x;
  abcd[1290] = 2.0E0*I_NAI_H2x3z_H3x2y_b-1*I_NAI_H2x3z_F3x;
  abcd[1291] = 2.0E0*I_NAI_Hx4y_H3x2y_b-1*I_NAI_Hx4y_F3x;
  abcd[1292] = 2.0E0*I_NAI_Hx3yz_H3x2y_b-1*I_NAI_Hx3yz_F3x;
  abcd[1293] = 2.0E0*I_NAI_Hx2y2z_H3x2y_b-1*I_NAI_Hx2y2z_F3x;
  abcd[1294] = 2.0E0*I_NAI_Hxy3z_H3x2y_b-1*I_NAI_Hxy3z_F3x;
  abcd[1295] = 2.0E0*I_NAI_Hx4z_H3x2y_b-1*I_NAI_Hx4z_F3x;
  abcd[1296] = 2.0E0*I_NAI_H5y_H3x2y_b-1*I_NAI_H5y_F3x;
  abcd[1297] = 2.0E0*I_NAI_H4yz_H3x2y_b-1*I_NAI_H4yz_F3x;
  abcd[1298] = 2.0E0*I_NAI_H3y2z_H3x2y_b-1*I_NAI_H3y2z_F3x;
  abcd[1299] = 2.0E0*I_NAI_H2y3z_H3x2y_b-1*I_NAI_H2y3z_F3x;
  abcd[1300] = 2.0E0*I_NAI_Hy4z_H3x2y_b-1*I_NAI_Hy4z_F3x;
  abcd[1301] = 2.0E0*I_NAI_H5z_H3x2y_b-1*I_NAI_H5z_F3x;
  abcd[1302] = 2.0E0*I_NAI_H5x_H3xyz_b;
  abcd[1303] = 2.0E0*I_NAI_H4xy_H3xyz_b;
  abcd[1304] = 2.0E0*I_NAI_H4xz_H3xyz_b;
  abcd[1305] = 2.0E0*I_NAI_H3x2y_H3xyz_b;
  abcd[1306] = 2.0E0*I_NAI_H3xyz_H3xyz_b;
  abcd[1307] = 2.0E0*I_NAI_H3x2z_H3xyz_b;
  abcd[1308] = 2.0E0*I_NAI_H2x3y_H3xyz_b;
  abcd[1309] = 2.0E0*I_NAI_H2x2yz_H3xyz_b;
  abcd[1310] = 2.0E0*I_NAI_H2xy2z_H3xyz_b;
  abcd[1311] = 2.0E0*I_NAI_H2x3z_H3xyz_b;
  abcd[1312] = 2.0E0*I_NAI_Hx4y_H3xyz_b;
  abcd[1313] = 2.0E0*I_NAI_Hx3yz_H3xyz_b;
  abcd[1314] = 2.0E0*I_NAI_Hx2y2z_H3xyz_b;
  abcd[1315] = 2.0E0*I_NAI_Hxy3z_H3xyz_b;
  abcd[1316] = 2.0E0*I_NAI_Hx4z_H3xyz_b;
  abcd[1317] = 2.0E0*I_NAI_H5y_H3xyz_b;
  abcd[1318] = 2.0E0*I_NAI_H4yz_H3xyz_b;
  abcd[1319] = 2.0E0*I_NAI_H3y2z_H3xyz_b;
  abcd[1320] = 2.0E0*I_NAI_H2y3z_H3xyz_b;
  abcd[1321] = 2.0E0*I_NAI_Hy4z_H3xyz_b;
  abcd[1322] = 2.0E0*I_NAI_H5z_H3xyz_b;
  abcd[1323] = 2.0E0*I_NAI_H5x_H2x3y_b-2*I_NAI_H5x_F2xy;
  abcd[1324] = 2.0E0*I_NAI_H4xy_H2x3y_b-2*I_NAI_H4xy_F2xy;
  abcd[1325] = 2.0E0*I_NAI_H4xz_H2x3y_b-2*I_NAI_H4xz_F2xy;
  abcd[1326] = 2.0E0*I_NAI_H3x2y_H2x3y_b-2*I_NAI_H3x2y_F2xy;
  abcd[1327] = 2.0E0*I_NAI_H3xyz_H2x3y_b-2*I_NAI_H3xyz_F2xy;
  abcd[1328] = 2.0E0*I_NAI_H3x2z_H2x3y_b-2*I_NAI_H3x2z_F2xy;
  abcd[1329] = 2.0E0*I_NAI_H2x3y_H2x3y_b-2*I_NAI_H2x3y_F2xy;
  abcd[1330] = 2.0E0*I_NAI_H2x2yz_H2x3y_b-2*I_NAI_H2x2yz_F2xy;
  abcd[1331] = 2.0E0*I_NAI_H2xy2z_H2x3y_b-2*I_NAI_H2xy2z_F2xy;
  abcd[1332] = 2.0E0*I_NAI_H2x3z_H2x3y_b-2*I_NAI_H2x3z_F2xy;
  abcd[1333] = 2.0E0*I_NAI_Hx4y_H2x3y_b-2*I_NAI_Hx4y_F2xy;
  abcd[1334] = 2.0E0*I_NAI_Hx3yz_H2x3y_b-2*I_NAI_Hx3yz_F2xy;
  abcd[1335] = 2.0E0*I_NAI_Hx2y2z_H2x3y_b-2*I_NAI_Hx2y2z_F2xy;
  abcd[1336] = 2.0E0*I_NAI_Hxy3z_H2x3y_b-2*I_NAI_Hxy3z_F2xy;
  abcd[1337] = 2.0E0*I_NAI_Hx4z_H2x3y_b-2*I_NAI_Hx4z_F2xy;
  abcd[1338] = 2.0E0*I_NAI_H5y_H2x3y_b-2*I_NAI_H5y_F2xy;
  abcd[1339] = 2.0E0*I_NAI_H4yz_H2x3y_b-2*I_NAI_H4yz_F2xy;
  abcd[1340] = 2.0E0*I_NAI_H3y2z_H2x3y_b-2*I_NAI_H3y2z_F2xy;
  abcd[1341] = 2.0E0*I_NAI_H2y3z_H2x3y_b-2*I_NAI_H2y3z_F2xy;
  abcd[1342] = 2.0E0*I_NAI_Hy4z_H2x3y_b-2*I_NAI_Hy4z_F2xy;
  abcd[1343] = 2.0E0*I_NAI_H5z_H2x3y_b-2*I_NAI_H5z_F2xy;
  abcd[1344] = 2.0E0*I_NAI_H5x_H2x2yz_b-1*I_NAI_H5x_F2xz;
  abcd[1345] = 2.0E0*I_NAI_H4xy_H2x2yz_b-1*I_NAI_H4xy_F2xz;
  abcd[1346] = 2.0E0*I_NAI_H4xz_H2x2yz_b-1*I_NAI_H4xz_F2xz;
  abcd[1347] = 2.0E0*I_NAI_H3x2y_H2x2yz_b-1*I_NAI_H3x2y_F2xz;
  abcd[1348] = 2.0E0*I_NAI_H3xyz_H2x2yz_b-1*I_NAI_H3xyz_F2xz;
  abcd[1349] = 2.0E0*I_NAI_H3x2z_H2x2yz_b-1*I_NAI_H3x2z_F2xz;
  abcd[1350] = 2.0E0*I_NAI_H2x3y_H2x2yz_b-1*I_NAI_H2x3y_F2xz;
  abcd[1351] = 2.0E0*I_NAI_H2x2yz_H2x2yz_b-1*I_NAI_H2x2yz_F2xz;
  abcd[1352] = 2.0E0*I_NAI_H2xy2z_H2x2yz_b-1*I_NAI_H2xy2z_F2xz;
  abcd[1353] = 2.0E0*I_NAI_H2x3z_H2x2yz_b-1*I_NAI_H2x3z_F2xz;
  abcd[1354] = 2.0E0*I_NAI_Hx4y_H2x2yz_b-1*I_NAI_Hx4y_F2xz;
  abcd[1355] = 2.0E0*I_NAI_Hx3yz_H2x2yz_b-1*I_NAI_Hx3yz_F2xz;
  abcd[1356] = 2.0E0*I_NAI_Hx2y2z_H2x2yz_b-1*I_NAI_Hx2y2z_F2xz;
  abcd[1357] = 2.0E0*I_NAI_Hxy3z_H2x2yz_b-1*I_NAI_Hxy3z_F2xz;
  abcd[1358] = 2.0E0*I_NAI_Hx4z_H2x2yz_b-1*I_NAI_Hx4z_F2xz;
  abcd[1359] = 2.0E0*I_NAI_H5y_H2x2yz_b-1*I_NAI_H5y_F2xz;
  abcd[1360] = 2.0E0*I_NAI_H4yz_H2x2yz_b-1*I_NAI_H4yz_F2xz;
  abcd[1361] = 2.0E0*I_NAI_H3y2z_H2x2yz_b-1*I_NAI_H3y2z_F2xz;
  abcd[1362] = 2.0E0*I_NAI_H2y3z_H2x2yz_b-1*I_NAI_H2y3z_F2xz;
  abcd[1363] = 2.0E0*I_NAI_Hy4z_H2x2yz_b-1*I_NAI_Hy4z_F2xz;
  abcd[1364] = 2.0E0*I_NAI_H5z_H2x2yz_b-1*I_NAI_H5z_F2xz;
  abcd[1365] = 2.0E0*I_NAI_H5x_H2xy2z_b;
  abcd[1366] = 2.0E0*I_NAI_H4xy_H2xy2z_b;
  abcd[1367] = 2.0E0*I_NAI_H4xz_H2xy2z_b;
  abcd[1368] = 2.0E0*I_NAI_H3x2y_H2xy2z_b;
  abcd[1369] = 2.0E0*I_NAI_H3xyz_H2xy2z_b;
  abcd[1370] = 2.0E0*I_NAI_H3x2z_H2xy2z_b;
  abcd[1371] = 2.0E0*I_NAI_H2x3y_H2xy2z_b;
  abcd[1372] = 2.0E0*I_NAI_H2x2yz_H2xy2z_b;
  abcd[1373] = 2.0E0*I_NAI_H2xy2z_H2xy2z_b;
  abcd[1374] = 2.0E0*I_NAI_H2x3z_H2xy2z_b;
  abcd[1375] = 2.0E0*I_NAI_Hx4y_H2xy2z_b;
  abcd[1376] = 2.0E0*I_NAI_Hx3yz_H2xy2z_b;
  abcd[1377] = 2.0E0*I_NAI_Hx2y2z_H2xy2z_b;
  abcd[1378] = 2.0E0*I_NAI_Hxy3z_H2xy2z_b;
  abcd[1379] = 2.0E0*I_NAI_Hx4z_H2xy2z_b;
  abcd[1380] = 2.0E0*I_NAI_H5y_H2xy2z_b;
  abcd[1381] = 2.0E0*I_NAI_H4yz_H2xy2z_b;
  abcd[1382] = 2.0E0*I_NAI_H3y2z_H2xy2z_b;
  abcd[1383] = 2.0E0*I_NAI_H2y3z_H2xy2z_b;
  abcd[1384] = 2.0E0*I_NAI_Hy4z_H2xy2z_b;
  abcd[1385] = 2.0E0*I_NAI_H5z_H2xy2z_b;
  abcd[1386] = 2.0E0*I_NAI_H5x_Hx4y_b-3*I_NAI_H5x_Fx2y;
  abcd[1387] = 2.0E0*I_NAI_H4xy_Hx4y_b-3*I_NAI_H4xy_Fx2y;
  abcd[1388] = 2.0E0*I_NAI_H4xz_Hx4y_b-3*I_NAI_H4xz_Fx2y;
  abcd[1389] = 2.0E0*I_NAI_H3x2y_Hx4y_b-3*I_NAI_H3x2y_Fx2y;
  abcd[1390] = 2.0E0*I_NAI_H3xyz_Hx4y_b-3*I_NAI_H3xyz_Fx2y;
  abcd[1391] = 2.0E0*I_NAI_H3x2z_Hx4y_b-3*I_NAI_H3x2z_Fx2y;
  abcd[1392] = 2.0E0*I_NAI_H2x3y_Hx4y_b-3*I_NAI_H2x3y_Fx2y;
  abcd[1393] = 2.0E0*I_NAI_H2x2yz_Hx4y_b-3*I_NAI_H2x2yz_Fx2y;
  abcd[1394] = 2.0E0*I_NAI_H2xy2z_Hx4y_b-3*I_NAI_H2xy2z_Fx2y;
  abcd[1395] = 2.0E0*I_NAI_H2x3z_Hx4y_b-3*I_NAI_H2x3z_Fx2y;
  abcd[1396] = 2.0E0*I_NAI_Hx4y_Hx4y_b-3*I_NAI_Hx4y_Fx2y;
  abcd[1397] = 2.0E0*I_NAI_Hx3yz_Hx4y_b-3*I_NAI_Hx3yz_Fx2y;
  abcd[1398] = 2.0E0*I_NAI_Hx2y2z_Hx4y_b-3*I_NAI_Hx2y2z_Fx2y;
  abcd[1399] = 2.0E0*I_NAI_Hxy3z_Hx4y_b-3*I_NAI_Hxy3z_Fx2y;
  abcd[1400] = 2.0E0*I_NAI_Hx4z_Hx4y_b-3*I_NAI_Hx4z_Fx2y;
  abcd[1401] = 2.0E0*I_NAI_H5y_Hx4y_b-3*I_NAI_H5y_Fx2y;
  abcd[1402] = 2.0E0*I_NAI_H4yz_Hx4y_b-3*I_NAI_H4yz_Fx2y;
  abcd[1403] = 2.0E0*I_NAI_H3y2z_Hx4y_b-3*I_NAI_H3y2z_Fx2y;
  abcd[1404] = 2.0E0*I_NAI_H2y3z_Hx4y_b-3*I_NAI_H2y3z_Fx2y;
  abcd[1405] = 2.0E0*I_NAI_Hy4z_Hx4y_b-3*I_NAI_Hy4z_Fx2y;
  abcd[1406] = 2.0E0*I_NAI_H5z_Hx4y_b-3*I_NAI_H5z_Fx2y;
  abcd[1407] = 2.0E0*I_NAI_H5x_Hx3yz_b-2*I_NAI_H5x_Fxyz;
  abcd[1408] = 2.0E0*I_NAI_H4xy_Hx3yz_b-2*I_NAI_H4xy_Fxyz;
  abcd[1409] = 2.0E0*I_NAI_H4xz_Hx3yz_b-2*I_NAI_H4xz_Fxyz;
  abcd[1410] = 2.0E0*I_NAI_H3x2y_Hx3yz_b-2*I_NAI_H3x2y_Fxyz;
  abcd[1411] = 2.0E0*I_NAI_H3xyz_Hx3yz_b-2*I_NAI_H3xyz_Fxyz;
  abcd[1412] = 2.0E0*I_NAI_H3x2z_Hx3yz_b-2*I_NAI_H3x2z_Fxyz;
  abcd[1413] = 2.0E0*I_NAI_H2x3y_Hx3yz_b-2*I_NAI_H2x3y_Fxyz;
  abcd[1414] = 2.0E0*I_NAI_H2x2yz_Hx3yz_b-2*I_NAI_H2x2yz_Fxyz;
  abcd[1415] = 2.0E0*I_NAI_H2xy2z_Hx3yz_b-2*I_NAI_H2xy2z_Fxyz;
  abcd[1416] = 2.0E0*I_NAI_H2x3z_Hx3yz_b-2*I_NAI_H2x3z_Fxyz;
  abcd[1417] = 2.0E0*I_NAI_Hx4y_Hx3yz_b-2*I_NAI_Hx4y_Fxyz;
  abcd[1418] = 2.0E0*I_NAI_Hx3yz_Hx3yz_b-2*I_NAI_Hx3yz_Fxyz;
  abcd[1419] = 2.0E0*I_NAI_Hx2y2z_Hx3yz_b-2*I_NAI_Hx2y2z_Fxyz;
  abcd[1420] = 2.0E0*I_NAI_Hxy3z_Hx3yz_b-2*I_NAI_Hxy3z_Fxyz;
  abcd[1421] = 2.0E0*I_NAI_Hx4z_Hx3yz_b-2*I_NAI_Hx4z_Fxyz;
  abcd[1422] = 2.0E0*I_NAI_H5y_Hx3yz_b-2*I_NAI_H5y_Fxyz;
  abcd[1423] = 2.0E0*I_NAI_H4yz_Hx3yz_b-2*I_NAI_H4yz_Fxyz;
  abcd[1424] = 2.0E0*I_NAI_H3y2z_Hx3yz_b-2*I_NAI_H3y2z_Fxyz;
  abcd[1425] = 2.0E0*I_NAI_H2y3z_Hx3yz_b-2*I_NAI_H2y3z_Fxyz;
  abcd[1426] = 2.0E0*I_NAI_Hy4z_Hx3yz_b-2*I_NAI_Hy4z_Fxyz;
  abcd[1427] = 2.0E0*I_NAI_H5z_Hx3yz_b-2*I_NAI_H5z_Fxyz;
  abcd[1428] = 2.0E0*I_NAI_H5x_Hx2y2z_b-1*I_NAI_H5x_Fx2z;
  abcd[1429] = 2.0E0*I_NAI_H4xy_Hx2y2z_b-1*I_NAI_H4xy_Fx2z;
  abcd[1430] = 2.0E0*I_NAI_H4xz_Hx2y2z_b-1*I_NAI_H4xz_Fx2z;
  abcd[1431] = 2.0E0*I_NAI_H3x2y_Hx2y2z_b-1*I_NAI_H3x2y_Fx2z;
  abcd[1432] = 2.0E0*I_NAI_H3xyz_Hx2y2z_b-1*I_NAI_H3xyz_Fx2z;
  abcd[1433] = 2.0E0*I_NAI_H3x2z_Hx2y2z_b-1*I_NAI_H3x2z_Fx2z;
  abcd[1434] = 2.0E0*I_NAI_H2x3y_Hx2y2z_b-1*I_NAI_H2x3y_Fx2z;
  abcd[1435] = 2.0E0*I_NAI_H2x2yz_Hx2y2z_b-1*I_NAI_H2x2yz_Fx2z;
  abcd[1436] = 2.0E0*I_NAI_H2xy2z_Hx2y2z_b-1*I_NAI_H2xy2z_Fx2z;
  abcd[1437] = 2.0E0*I_NAI_H2x3z_Hx2y2z_b-1*I_NAI_H2x3z_Fx2z;
  abcd[1438] = 2.0E0*I_NAI_Hx4y_Hx2y2z_b-1*I_NAI_Hx4y_Fx2z;
  abcd[1439] = 2.0E0*I_NAI_Hx3yz_Hx2y2z_b-1*I_NAI_Hx3yz_Fx2z;
  abcd[1440] = 2.0E0*I_NAI_Hx2y2z_Hx2y2z_b-1*I_NAI_Hx2y2z_Fx2z;
  abcd[1441] = 2.0E0*I_NAI_Hxy3z_Hx2y2z_b-1*I_NAI_Hxy3z_Fx2z;
  abcd[1442] = 2.0E0*I_NAI_Hx4z_Hx2y2z_b-1*I_NAI_Hx4z_Fx2z;
  abcd[1443] = 2.0E0*I_NAI_H5y_Hx2y2z_b-1*I_NAI_H5y_Fx2z;
  abcd[1444] = 2.0E0*I_NAI_H4yz_Hx2y2z_b-1*I_NAI_H4yz_Fx2z;
  abcd[1445] = 2.0E0*I_NAI_H3y2z_Hx2y2z_b-1*I_NAI_H3y2z_Fx2z;
  abcd[1446] = 2.0E0*I_NAI_H2y3z_Hx2y2z_b-1*I_NAI_H2y3z_Fx2z;
  abcd[1447] = 2.0E0*I_NAI_Hy4z_Hx2y2z_b-1*I_NAI_Hy4z_Fx2z;
  abcd[1448] = 2.0E0*I_NAI_H5z_Hx2y2z_b-1*I_NAI_H5z_Fx2z;
  abcd[1449] = 2.0E0*I_NAI_H5x_Hxy3z_b;
  abcd[1450] = 2.0E0*I_NAI_H4xy_Hxy3z_b;
  abcd[1451] = 2.0E0*I_NAI_H4xz_Hxy3z_b;
  abcd[1452] = 2.0E0*I_NAI_H3x2y_Hxy3z_b;
  abcd[1453] = 2.0E0*I_NAI_H3xyz_Hxy3z_b;
  abcd[1454] = 2.0E0*I_NAI_H3x2z_Hxy3z_b;
  abcd[1455] = 2.0E0*I_NAI_H2x3y_Hxy3z_b;
  abcd[1456] = 2.0E0*I_NAI_H2x2yz_Hxy3z_b;
  abcd[1457] = 2.0E0*I_NAI_H2xy2z_Hxy3z_b;
  abcd[1458] = 2.0E0*I_NAI_H2x3z_Hxy3z_b;
  abcd[1459] = 2.0E0*I_NAI_Hx4y_Hxy3z_b;
  abcd[1460] = 2.0E0*I_NAI_Hx3yz_Hxy3z_b;
  abcd[1461] = 2.0E0*I_NAI_Hx2y2z_Hxy3z_b;
  abcd[1462] = 2.0E0*I_NAI_Hxy3z_Hxy3z_b;
  abcd[1463] = 2.0E0*I_NAI_Hx4z_Hxy3z_b;
  abcd[1464] = 2.0E0*I_NAI_H5y_Hxy3z_b;
  abcd[1465] = 2.0E0*I_NAI_H4yz_Hxy3z_b;
  abcd[1466] = 2.0E0*I_NAI_H3y2z_Hxy3z_b;
  abcd[1467] = 2.0E0*I_NAI_H2y3z_Hxy3z_b;
  abcd[1468] = 2.0E0*I_NAI_Hy4z_Hxy3z_b;
  abcd[1469] = 2.0E0*I_NAI_H5z_Hxy3z_b;
  abcd[1470] = 2.0E0*I_NAI_H5x_H5y_b-4*I_NAI_H5x_F3y;
  abcd[1471] = 2.0E0*I_NAI_H4xy_H5y_b-4*I_NAI_H4xy_F3y;
  abcd[1472] = 2.0E0*I_NAI_H4xz_H5y_b-4*I_NAI_H4xz_F3y;
  abcd[1473] = 2.0E0*I_NAI_H3x2y_H5y_b-4*I_NAI_H3x2y_F3y;
  abcd[1474] = 2.0E0*I_NAI_H3xyz_H5y_b-4*I_NAI_H3xyz_F3y;
  abcd[1475] = 2.0E0*I_NAI_H3x2z_H5y_b-4*I_NAI_H3x2z_F3y;
  abcd[1476] = 2.0E0*I_NAI_H2x3y_H5y_b-4*I_NAI_H2x3y_F3y;
  abcd[1477] = 2.0E0*I_NAI_H2x2yz_H5y_b-4*I_NAI_H2x2yz_F3y;
  abcd[1478] = 2.0E0*I_NAI_H2xy2z_H5y_b-4*I_NAI_H2xy2z_F3y;
  abcd[1479] = 2.0E0*I_NAI_H2x3z_H5y_b-4*I_NAI_H2x3z_F3y;
  abcd[1480] = 2.0E0*I_NAI_Hx4y_H5y_b-4*I_NAI_Hx4y_F3y;
  abcd[1481] = 2.0E0*I_NAI_Hx3yz_H5y_b-4*I_NAI_Hx3yz_F3y;
  abcd[1482] = 2.0E0*I_NAI_Hx2y2z_H5y_b-4*I_NAI_Hx2y2z_F3y;
  abcd[1483] = 2.0E0*I_NAI_Hxy3z_H5y_b-4*I_NAI_Hxy3z_F3y;
  abcd[1484] = 2.0E0*I_NAI_Hx4z_H5y_b-4*I_NAI_Hx4z_F3y;
  abcd[1485] = 2.0E0*I_NAI_H5y_H5y_b-4*I_NAI_H5y_F3y;
  abcd[1486] = 2.0E0*I_NAI_H4yz_H5y_b-4*I_NAI_H4yz_F3y;
  abcd[1487] = 2.0E0*I_NAI_H3y2z_H5y_b-4*I_NAI_H3y2z_F3y;
  abcd[1488] = 2.0E0*I_NAI_H2y3z_H5y_b-4*I_NAI_H2y3z_F3y;
  abcd[1489] = 2.0E0*I_NAI_Hy4z_H5y_b-4*I_NAI_Hy4z_F3y;
  abcd[1490] = 2.0E0*I_NAI_H5z_H5y_b-4*I_NAI_H5z_F3y;
  abcd[1491] = 2.0E0*I_NAI_H5x_H4yz_b-3*I_NAI_H5x_F2yz;
  abcd[1492] = 2.0E0*I_NAI_H4xy_H4yz_b-3*I_NAI_H4xy_F2yz;
  abcd[1493] = 2.0E0*I_NAI_H4xz_H4yz_b-3*I_NAI_H4xz_F2yz;
  abcd[1494] = 2.0E0*I_NAI_H3x2y_H4yz_b-3*I_NAI_H3x2y_F2yz;
  abcd[1495] = 2.0E0*I_NAI_H3xyz_H4yz_b-3*I_NAI_H3xyz_F2yz;
  abcd[1496] = 2.0E0*I_NAI_H3x2z_H4yz_b-3*I_NAI_H3x2z_F2yz;
  abcd[1497] = 2.0E0*I_NAI_H2x3y_H4yz_b-3*I_NAI_H2x3y_F2yz;
  abcd[1498] = 2.0E0*I_NAI_H2x2yz_H4yz_b-3*I_NAI_H2x2yz_F2yz;
  abcd[1499] = 2.0E0*I_NAI_H2xy2z_H4yz_b-3*I_NAI_H2xy2z_F2yz;
  abcd[1500] = 2.0E0*I_NAI_H2x3z_H4yz_b-3*I_NAI_H2x3z_F2yz;
  abcd[1501] = 2.0E0*I_NAI_Hx4y_H4yz_b-3*I_NAI_Hx4y_F2yz;
  abcd[1502] = 2.0E0*I_NAI_Hx3yz_H4yz_b-3*I_NAI_Hx3yz_F2yz;
  abcd[1503] = 2.0E0*I_NAI_Hx2y2z_H4yz_b-3*I_NAI_Hx2y2z_F2yz;
  abcd[1504] = 2.0E0*I_NAI_Hxy3z_H4yz_b-3*I_NAI_Hxy3z_F2yz;
  abcd[1505] = 2.0E0*I_NAI_Hx4z_H4yz_b-3*I_NAI_Hx4z_F2yz;
  abcd[1506] = 2.0E0*I_NAI_H5y_H4yz_b-3*I_NAI_H5y_F2yz;
  abcd[1507] = 2.0E0*I_NAI_H4yz_H4yz_b-3*I_NAI_H4yz_F2yz;
  abcd[1508] = 2.0E0*I_NAI_H3y2z_H4yz_b-3*I_NAI_H3y2z_F2yz;
  abcd[1509] = 2.0E0*I_NAI_H2y3z_H4yz_b-3*I_NAI_H2y3z_F2yz;
  abcd[1510] = 2.0E0*I_NAI_Hy4z_H4yz_b-3*I_NAI_Hy4z_F2yz;
  abcd[1511] = 2.0E0*I_NAI_H5z_H4yz_b-3*I_NAI_H5z_F2yz;
  abcd[1512] = 2.0E0*I_NAI_H5x_H3y2z_b-2*I_NAI_H5x_Fy2z;
  abcd[1513] = 2.0E0*I_NAI_H4xy_H3y2z_b-2*I_NAI_H4xy_Fy2z;
  abcd[1514] = 2.0E0*I_NAI_H4xz_H3y2z_b-2*I_NAI_H4xz_Fy2z;
  abcd[1515] = 2.0E0*I_NAI_H3x2y_H3y2z_b-2*I_NAI_H3x2y_Fy2z;
  abcd[1516] = 2.0E0*I_NAI_H3xyz_H3y2z_b-2*I_NAI_H3xyz_Fy2z;
  abcd[1517] = 2.0E0*I_NAI_H3x2z_H3y2z_b-2*I_NAI_H3x2z_Fy2z;
  abcd[1518] = 2.0E0*I_NAI_H2x3y_H3y2z_b-2*I_NAI_H2x3y_Fy2z;
  abcd[1519] = 2.0E0*I_NAI_H2x2yz_H3y2z_b-2*I_NAI_H2x2yz_Fy2z;
  abcd[1520] = 2.0E0*I_NAI_H2xy2z_H3y2z_b-2*I_NAI_H2xy2z_Fy2z;
  abcd[1521] = 2.0E0*I_NAI_H2x3z_H3y2z_b-2*I_NAI_H2x3z_Fy2z;
  abcd[1522] = 2.0E0*I_NAI_Hx4y_H3y2z_b-2*I_NAI_Hx4y_Fy2z;
  abcd[1523] = 2.0E0*I_NAI_Hx3yz_H3y2z_b-2*I_NAI_Hx3yz_Fy2z;
  abcd[1524] = 2.0E0*I_NAI_Hx2y2z_H3y2z_b-2*I_NAI_Hx2y2z_Fy2z;
  abcd[1525] = 2.0E0*I_NAI_Hxy3z_H3y2z_b-2*I_NAI_Hxy3z_Fy2z;
  abcd[1526] = 2.0E0*I_NAI_Hx4z_H3y2z_b-2*I_NAI_Hx4z_Fy2z;
  abcd[1527] = 2.0E0*I_NAI_H5y_H3y2z_b-2*I_NAI_H5y_Fy2z;
  abcd[1528] = 2.0E0*I_NAI_H4yz_H3y2z_b-2*I_NAI_H4yz_Fy2z;
  abcd[1529] = 2.0E0*I_NAI_H3y2z_H3y2z_b-2*I_NAI_H3y2z_Fy2z;
  abcd[1530] = 2.0E0*I_NAI_H2y3z_H3y2z_b-2*I_NAI_H2y3z_Fy2z;
  abcd[1531] = 2.0E0*I_NAI_Hy4z_H3y2z_b-2*I_NAI_Hy4z_Fy2z;
  abcd[1532] = 2.0E0*I_NAI_H5z_H3y2z_b-2*I_NAI_H5z_Fy2z;
  abcd[1533] = 2.0E0*I_NAI_H5x_H2y3z_b-1*I_NAI_H5x_F3z;
  abcd[1534] = 2.0E0*I_NAI_H4xy_H2y3z_b-1*I_NAI_H4xy_F3z;
  abcd[1535] = 2.0E0*I_NAI_H4xz_H2y3z_b-1*I_NAI_H4xz_F3z;
  abcd[1536] = 2.0E0*I_NAI_H3x2y_H2y3z_b-1*I_NAI_H3x2y_F3z;
  abcd[1537] = 2.0E0*I_NAI_H3xyz_H2y3z_b-1*I_NAI_H3xyz_F3z;
  abcd[1538] = 2.0E0*I_NAI_H3x2z_H2y3z_b-1*I_NAI_H3x2z_F3z;
  abcd[1539] = 2.0E0*I_NAI_H2x3y_H2y3z_b-1*I_NAI_H2x3y_F3z;
  abcd[1540] = 2.0E0*I_NAI_H2x2yz_H2y3z_b-1*I_NAI_H2x2yz_F3z;
  abcd[1541] = 2.0E0*I_NAI_H2xy2z_H2y3z_b-1*I_NAI_H2xy2z_F3z;
  abcd[1542] = 2.0E0*I_NAI_H2x3z_H2y3z_b-1*I_NAI_H2x3z_F3z;
  abcd[1543] = 2.0E0*I_NAI_Hx4y_H2y3z_b-1*I_NAI_Hx4y_F3z;
  abcd[1544] = 2.0E0*I_NAI_Hx3yz_H2y3z_b-1*I_NAI_Hx3yz_F3z;
  abcd[1545] = 2.0E0*I_NAI_Hx2y2z_H2y3z_b-1*I_NAI_Hx2y2z_F3z;
  abcd[1546] = 2.0E0*I_NAI_Hxy3z_H2y3z_b-1*I_NAI_Hxy3z_F3z;
  abcd[1547] = 2.0E0*I_NAI_Hx4z_H2y3z_b-1*I_NAI_Hx4z_F3z;
  abcd[1548] = 2.0E0*I_NAI_H5y_H2y3z_b-1*I_NAI_H5y_F3z;
  abcd[1549] = 2.0E0*I_NAI_H4yz_H2y3z_b-1*I_NAI_H4yz_F3z;
  abcd[1550] = 2.0E0*I_NAI_H3y2z_H2y3z_b-1*I_NAI_H3y2z_F3z;
  abcd[1551] = 2.0E0*I_NAI_H2y3z_H2y3z_b-1*I_NAI_H2y3z_F3z;
  abcd[1552] = 2.0E0*I_NAI_Hy4z_H2y3z_b-1*I_NAI_Hy4z_F3z;
  abcd[1553] = 2.0E0*I_NAI_H5z_H2y3z_b-1*I_NAI_H5z_F3z;
  abcd[1554] = 2.0E0*I_NAI_H5x_Hy4z_b;
  abcd[1555] = 2.0E0*I_NAI_H4xy_Hy4z_b;
  abcd[1556] = 2.0E0*I_NAI_H4xz_Hy4z_b;
  abcd[1557] = 2.0E0*I_NAI_H3x2y_Hy4z_b;
  abcd[1558] = 2.0E0*I_NAI_H3xyz_Hy4z_b;
  abcd[1559] = 2.0E0*I_NAI_H3x2z_Hy4z_b;
  abcd[1560] = 2.0E0*I_NAI_H2x3y_Hy4z_b;
  abcd[1561] = 2.0E0*I_NAI_H2x2yz_Hy4z_b;
  abcd[1562] = 2.0E0*I_NAI_H2xy2z_Hy4z_b;
  abcd[1563] = 2.0E0*I_NAI_H2x3z_Hy4z_b;
  abcd[1564] = 2.0E0*I_NAI_Hx4y_Hy4z_b;
  abcd[1565] = 2.0E0*I_NAI_Hx3yz_Hy4z_b;
  abcd[1566] = 2.0E0*I_NAI_Hx2y2z_Hy4z_b;
  abcd[1567] = 2.0E0*I_NAI_Hxy3z_Hy4z_b;
  abcd[1568] = 2.0E0*I_NAI_Hx4z_Hy4z_b;
  abcd[1569] = 2.0E0*I_NAI_H5y_Hy4z_b;
  abcd[1570] = 2.0E0*I_NAI_H4yz_Hy4z_b;
  abcd[1571] = 2.0E0*I_NAI_H3y2z_Hy4z_b;
  abcd[1572] = 2.0E0*I_NAI_H2y3z_Hy4z_b;
  abcd[1573] = 2.0E0*I_NAI_Hy4z_Hy4z_b;
  abcd[1574] = 2.0E0*I_NAI_H5z_Hy4z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_G_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_H_b
   * RHS shell quartet name: SQ_NAI_H_F
   ************************************************************/
  abcd[1575] = 2.0E0*I_NAI_H5x_H4xz_b;
  abcd[1576] = 2.0E0*I_NAI_H4xy_H4xz_b;
  abcd[1577] = 2.0E0*I_NAI_H4xz_H4xz_b;
  abcd[1578] = 2.0E0*I_NAI_H3x2y_H4xz_b;
  abcd[1579] = 2.0E0*I_NAI_H3xyz_H4xz_b;
  abcd[1580] = 2.0E0*I_NAI_H3x2z_H4xz_b;
  abcd[1581] = 2.0E0*I_NAI_H2x3y_H4xz_b;
  abcd[1582] = 2.0E0*I_NAI_H2x2yz_H4xz_b;
  abcd[1583] = 2.0E0*I_NAI_H2xy2z_H4xz_b;
  abcd[1584] = 2.0E0*I_NAI_H2x3z_H4xz_b;
  abcd[1585] = 2.0E0*I_NAI_Hx4y_H4xz_b;
  abcd[1586] = 2.0E0*I_NAI_Hx3yz_H4xz_b;
  abcd[1587] = 2.0E0*I_NAI_Hx2y2z_H4xz_b;
  abcd[1588] = 2.0E0*I_NAI_Hxy3z_H4xz_b;
  abcd[1589] = 2.0E0*I_NAI_Hx4z_H4xz_b;
  abcd[1590] = 2.0E0*I_NAI_H5y_H4xz_b;
  abcd[1591] = 2.0E0*I_NAI_H4yz_H4xz_b;
  abcd[1592] = 2.0E0*I_NAI_H3y2z_H4xz_b;
  abcd[1593] = 2.0E0*I_NAI_H2y3z_H4xz_b;
  abcd[1594] = 2.0E0*I_NAI_Hy4z_H4xz_b;
  abcd[1595] = 2.0E0*I_NAI_H5z_H4xz_b;
  abcd[1596] = 2.0E0*I_NAI_H5x_H3xyz_b;
  abcd[1597] = 2.0E0*I_NAI_H4xy_H3xyz_b;
  abcd[1598] = 2.0E0*I_NAI_H4xz_H3xyz_b;
  abcd[1599] = 2.0E0*I_NAI_H3x2y_H3xyz_b;
  abcd[1600] = 2.0E0*I_NAI_H3xyz_H3xyz_b;
  abcd[1601] = 2.0E0*I_NAI_H3x2z_H3xyz_b;
  abcd[1602] = 2.0E0*I_NAI_H2x3y_H3xyz_b;
  abcd[1603] = 2.0E0*I_NAI_H2x2yz_H3xyz_b;
  abcd[1604] = 2.0E0*I_NAI_H2xy2z_H3xyz_b;
  abcd[1605] = 2.0E0*I_NAI_H2x3z_H3xyz_b;
  abcd[1606] = 2.0E0*I_NAI_Hx4y_H3xyz_b;
  abcd[1607] = 2.0E0*I_NAI_Hx3yz_H3xyz_b;
  abcd[1608] = 2.0E0*I_NAI_Hx2y2z_H3xyz_b;
  abcd[1609] = 2.0E0*I_NAI_Hxy3z_H3xyz_b;
  abcd[1610] = 2.0E0*I_NAI_Hx4z_H3xyz_b;
  abcd[1611] = 2.0E0*I_NAI_H5y_H3xyz_b;
  abcd[1612] = 2.0E0*I_NAI_H4yz_H3xyz_b;
  abcd[1613] = 2.0E0*I_NAI_H3y2z_H3xyz_b;
  abcd[1614] = 2.0E0*I_NAI_H2y3z_H3xyz_b;
  abcd[1615] = 2.0E0*I_NAI_Hy4z_H3xyz_b;
  abcd[1616] = 2.0E0*I_NAI_H5z_H3xyz_b;
  abcd[1617] = 2.0E0*I_NAI_H5x_H3x2z_b-1*I_NAI_H5x_F3x;
  abcd[1618] = 2.0E0*I_NAI_H4xy_H3x2z_b-1*I_NAI_H4xy_F3x;
  abcd[1619] = 2.0E0*I_NAI_H4xz_H3x2z_b-1*I_NAI_H4xz_F3x;
  abcd[1620] = 2.0E0*I_NAI_H3x2y_H3x2z_b-1*I_NAI_H3x2y_F3x;
  abcd[1621] = 2.0E0*I_NAI_H3xyz_H3x2z_b-1*I_NAI_H3xyz_F3x;
  abcd[1622] = 2.0E0*I_NAI_H3x2z_H3x2z_b-1*I_NAI_H3x2z_F3x;
  abcd[1623] = 2.0E0*I_NAI_H2x3y_H3x2z_b-1*I_NAI_H2x3y_F3x;
  abcd[1624] = 2.0E0*I_NAI_H2x2yz_H3x2z_b-1*I_NAI_H2x2yz_F3x;
  abcd[1625] = 2.0E0*I_NAI_H2xy2z_H3x2z_b-1*I_NAI_H2xy2z_F3x;
  abcd[1626] = 2.0E0*I_NAI_H2x3z_H3x2z_b-1*I_NAI_H2x3z_F3x;
  abcd[1627] = 2.0E0*I_NAI_Hx4y_H3x2z_b-1*I_NAI_Hx4y_F3x;
  abcd[1628] = 2.0E0*I_NAI_Hx3yz_H3x2z_b-1*I_NAI_Hx3yz_F3x;
  abcd[1629] = 2.0E0*I_NAI_Hx2y2z_H3x2z_b-1*I_NAI_Hx2y2z_F3x;
  abcd[1630] = 2.0E0*I_NAI_Hxy3z_H3x2z_b-1*I_NAI_Hxy3z_F3x;
  abcd[1631] = 2.0E0*I_NAI_Hx4z_H3x2z_b-1*I_NAI_Hx4z_F3x;
  abcd[1632] = 2.0E0*I_NAI_H5y_H3x2z_b-1*I_NAI_H5y_F3x;
  abcd[1633] = 2.0E0*I_NAI_H4yz_H3x2z_b-1*I_NAI_H4yz_F3x;
  abcd[1634] = 2.0E0*I_NAI_H3y2z_H3x2z_b-1*I_NAI_H3y2z_F3x;
  abcd[1635] = 2.0E0*I_NAI_H2y3z_H3x2z_b-1*I_NAI_H2y3z_F3x;
  abcd[1636] = 2.0E0*I_NAI_Hy4z_H3x2z_b-1*I_NAI_Hy4z_F3x;
  abcd[1637] = 2.0E0*I_NAI_H5z_H3x2z_b-1*I_NAI_H5z_F3x;
  abcd[1638] = 2.0E0*I_NAI_H5x_H2x2yz_b;
  abcd[1639] = 2.0E0*I_NAI_H4xy_H2x2yz_b;
  abcd[1640] = 2.0E0*I_NAI_H4xz_H2x2yz_b;
  abcd[1641] = 2.0E0*I_NAI_H3x2y_H2x2yz_b;
  abcd[1642] = 2.0E0*I_NAI_H3xyz_H2x2yz_b;
  abcd[1643] = 2.0E0*I_NAI_H3x2z_H2x2yz_b;
  abcd[1644] = 2.0E0*I_NAI_H2x3y_H2x2yz_b;
  abcd[1645] = 2.0E0*I_NAI_H2x2yz_H2x2yz_b;
  abcd[1646] = 2.0E0*I_NAI_H2xy2z_H2x2yz_b;
  abcd[1647] = 2.0E0*I_NAI_H2x3z_H2x2yz_b;
  abcd[1648] = 2.0E0*I_NAI_Hx4y_H2x2yz_b;
  abcd[1649] = 2.0E0*I_NAI_Hx3yz_H2x2yz_b;
  abcd[1650] = 2.0E0*I_NAI_Hx2y2z_H2x2yz_b;
  abcd[1651] = 2.0E0*I_NAI_Hxy3z_H2x2yz_b;
  abcd[1652] = 2.0E0*I_NAI_Hx4z_H2x2yz_b;
  abcd[1653] = 2.0E0*I_NAI_H5y_H2x2yz_b;
  abcd[1654] = 2.0E0*I_NAI_H4yz_H2x2yz_b;
  abcd[1655] = 2.0E0*I_NAI_H3y2z_H2x2yz_b;
  abcd[1656] = 2.0E0*I_NAI_H2y3z_H2x2yz_b;
  abcd[1657] = 2.0E0*I_NAI_Hy4z_H2x2yz_b;
  abcd[1658] = 2.0E0*I_NAI_H5z_H2x2yz_b;
  abcd[1659] = 2.0E0*I_NAI_H5x_H2xy2z_b-1*I_NAI_H5x_F2xy;
  abcd[1660] = 2.0E0*I_NAI_H4xy_H2xy2z_b-1*I_NAI_H4xy_F2xy;
  abcd[1661] = 2.0E0*I_NAI_H4xz_H2xy2z_b-1*I_NAI_H4xz_F2xy;
  abcd[1662] = 2.0E0*I_NAI_H3x2y_H2xy2z_b-1*I_NAI_H3x2y_F2xy;
  abcd[1663] = 2.0E0*I_NAI_H3xyz_H2xy2z_b-1*I_NAI_H3xyz_F2xy;
  abcd[1664] = 2.0E0*I_NAI_H3x2z_H2xy2z_b-1*I_NAI_H3x2z_F2xy;
  abcd[1665] = 2.0E0*I_NAI_H2x3y_H2xy2z_b-1*I_NAI_H2x3y_F2xy;
  abcd[1666] = 2.0E0*I_NAI_H2x2yz_H2xy2z_b-1*I_NAI_H2x2yz_F2xy;
  abcd[1667] = 2.0E0*I_NAI_H2xy2z_H2xy2z_b-1*I_NAI_H2xy2z_F2xy;
  abcd[1668] = 2.0E0*I_NAI_H2x3z_H2xy2z_b-1*I_NAI_H2x3z_F2xy;
  abcd[1669] = 2.0E0*I_NAI_Hx4y_H2xy2z_b-1*I_NAI_Hx4y_F2xy;
  abcd[1670] = 2.0E0*I_NAI_Hx3yz_H2xy2z_b-1*I_NAI_Hx3yz_F2xy;
  abcd[1671] = 2.0E0*I_NAI_Hx2y2z_H2xy2z_b-1*I_NAI_Hx2y2z_F2xy;
  abcd[1672] = 2.0E0*I_NAI_Hxy3z_H2xy2z_b-1*I_NAI_Hxy3z_F2xy;
  abcd[1673] = 2.0E0*I_NAI_Hx4z_H2xy2z_b-1*I_NAI_Hx4z_F2xy;
  abcd[1674] = 2.0E0*I_NAI_H5y_H2xy2z_b-1*I_NAI_H5y_F2xy;
  abcd[1675] = 2.0E0*I_NAI_H4yz_H2xy2z_b-1*I_NAI_H4yz_F2xy;
  abcd[1676] = 2.0E0*I_NAI_H3y2z_H2xy2z_b-1*I_NAI_H3y2z_F2xy;
  abcd[1677] = 2.0E0*I_NAI_H2y3z_H2xy2z_b-1*I_NAI_H2y3z_F2xy;
  abcd[1678] = 2.0E0*I_NAI_Hy4z_H2xy2z_b-1*I_NAI_Hy4z_F2xy;
  abcd[1679] = 2.0E0*I_NAI_H5z_H2xy2z_b-1*I_NAI_H5z_F2xy;
  abcd[1680] = 2.0E0*I_NAI_H5x_H2x3z_b-2*I_NAI_H5x_F2xz;
  abcd[1681] = 2.0E0*I_NAI_H4xy_H2x3z_b-2*I_NAI_H4xy_F2xz;
  abcd[1682] = 2.0E0*I_NAI_H4xz_H2x3z_b-2*I_NAI_H4xz_F2xz;
  abcd[1683] = 2.0E0*I_NAI_H3x2y_H2x3z_b-2*I_NAI_H3x2y_F2xz;
  abcd[1684] = 2.0E0*I_NAI_H3xyz_H2x3z_b-2*I_NAI_H3xyz_F2xz;
  abcd[1685] = 2.0E0*I_NAI_H3x2z_H2x3z_b-2*I_NAI_H3x2z_F2xz;
  abcd[1686] = 2.0E0*I_NAI_H2x3y_H2x3z_b-2*I_NAI_H2x3y_F2xz;
  abcd[1687] = 2.0E0*I_NAI_H2x2yz_H2x3z_b-2*I_NAI_H2x2yz_F2xz;
  abcd[1688] = 2.0E0*I_NAI_H2xy2z_H2x3z_b-2*I_NAI_H2xy2z_F2xz;
  abcd[1689] = 2.0E0*I_NAI_H2x3z_H2x3z_b-2*I_NAI_H2x3z_F2xz;
  abcd[1690] = 2.0E0*I_NAI_Hx4y_H2x3z_b-2*I_NAI_Hx4y_F2xz;
  abcd[1691] = 2.0E0*I_NAI_Hx3yz_H2x3z_b-2*I_NAI_Hx3yz_F2xz;
  abcd[1692] = 2.0E0*I_NAI_Hx2y2z_H2x3z_b-2*I_NAI_Hx2y2z_F2xz;
  abcd[1693] = 2.0E0*I_NAI_Hxy3z_H2x3z_b-2*I_NAI_Hxy3z_F2xz;
  abcd[1694] = 2.0E0*I_NAI_Hx4z_H2x3z_b-2*I_NAI_Hx4z_F2xz;
  abcd[1695] = 2.0E0*I_NAI_H5y_H2x3z_b-2*I_NAI_H5y_F2xz;
  abcd[1696] = 2.0E0*I_NAI_H4yz_H2x3z_b-2*I_NAI_H4yz_F2xz;
  abcd[1697] = 2.0E0*I_NAI_H3y2z_H2x3z_b-2*I_NAI_H3y2z_F2xz;
  abcd[1698] = 2.0E0*I_NAI_H2y3z_H2x3z_b-2*I_NAI_H2y3z_F2xz;
  abcd[1699] = 2.0E0*I_NAI_Hy4z_H2x3z_b-2*I_NAI_Hy4z_F2xz;
  abcd[1700] = 2.0E0*I_NAI_H5z_H2x3z_b-2*I_NAI_H5z_F2xz;
  abcd[1701] = 2.0E0*I_NAI_H5x_Hx3yz_b;
  abcd[1702] = 2.0E0*I_NAI_H4xy_Hx3yz_b;
  abcd[1703] = 2.0E0*I_NAI_H4xz_Hx3yz_b;
  abcd[1704] = 2.0E0*I_NAI_H3x2y_Hx3yz_b;
  abcd[1705] = 2.0E0*I_NAI_H3xyz_Hx3yz_b;
  abcd[1706] = 2.0E0*I_NAI_H3x2z_Hx3yz_b;
  abcd[1707] = 2.0E0*I_NAI_H2x3y_Hx3yz_b;
  abcd[1708] = 2.0E0*I_NAI_H2x2yz_Hx3yz_b;
  abcd[1709] = 2.0E0*I_NAI_H2xy2z_Hx3yz_b;
  abcd[1710] = 2.0E0*I_NAI_H2x3z_Hx3yz_b;
  abcd[1711] = 2.0E0*I_NAI_Hx4y_Hx3yz_b;
  abcd[1712] = 2.0E0*I_NAI_Hx3yz_Hx3yz_b;
  abcd[1713] = 2.0E0*I_NAI_Hx2y2z_Hx3yz_b;
  abcd[1714] = 2.0E0*I_NAI_Hxy3z_Hx3yz_b;
  abcd[1715] = 2.0E0*I_NAI_Hx4z_Hx3yz_b;
  abcd[1716] = 2.0E0*I_NAI_H5y_Hx3yz_b;
  abcd[1717] = 2.0E0*I_NAI_H4yz_Hx3yz_b;
  abcd[1718] = 2.0E0*I_NAI_H3y2z_Hx3yz_b;
  abcd[1719] = 2.0E0*I_NAI_H2y3z_Hx3yz_b;
  abcd[1720] = 2.0E0*I_NAI_Hy4z_Hx3yz_b;
  abcd[1721] = 2.0E0*I_NAI_H5z_Hx3yz_b;
  abcd[1722] = 2.0E0*I_NAI_H5x_Hx2y2z_b-1*I_NAI_H5x_Fx2y;
  abcd[1723] = 2.0E0*I_NAI_H4xy_Hx2y2z_b-1*I_NAI_H4xy_Fx2y;
  abcd[1724] = 2.0E0*I_NAI_H4xz_Hx2y2z_b-1*I_NAI_H4xz_Fx2y;
  abcd[1725] = 2.0E0*I_NAI_H3x2y_Hx2y2z_b-1*I_NAI_H3x2y_Fx2y;
  abcd[1726] = 2.0E0*I_NAI_H3xyz_Hx2y2z_b-1*I_NAI_H3xyz_Fx2y;
  abcd[1727] = 2.0E0*I_NAI_H3x2z_Hx2y2z_b-1*I_NAI_H3x2z_Fx2y;
  abcd[1728] = 2.0E0*I_NAI_H2x3y_Hx2y2z_b-1*I_NAI_H2x3y_Fx2y;
  abcd[1729] = 2.0E0*I_NAI_H2x2yz_Hx2y2z_b-1*I_NAI_H2x2yz_Fx2y;
  abcd[1730] = 2.0E0*I_NAI_H2xy2z_Hx2y2z_b-1*I_NAI_H2xy2z_Fx2y;
  abcd[1731] = 2.0E0*I_NAI_H2x3z_Hx2y2z_b-1*I_NAI_H2x3z_Fx2y;
  abcd[1732] = 2.0E0*I_NAI_Hx4y_Hx2y2z_b-1*I_NAI_Hx4y_Fx2y;
  abcd[1733] = 2.0E0*I_NAI_Hx3yz_Hx2y2z_b-1*I_NAI_Hx3yz_Fx2y;
  abcd[1734] = 2.0E0*I_NAI_Hx2y2z_Hx2y2z_b-1*I_NAI_Hx2y2z_Fx2y;
  abcd[1735] = 2.0E0*I_NAI_Hxy3z_Hx2y2z_b-1*I_NAI_Hxy3z_Fx2y;
  abcd[1736] = 2.0E0*I_NAI_Hx4z_Hx2y2z_b-1*I_NAI_Hx4z_Fx2y;
  abcd[1737] = 2.0E0*I_NAI_H5y_Hx2y2z_b-1*I_NAI_H5y_Fx2y;
  abcd[1738] = 2.0E0*I_NAI_H4yz_Hx2y2z_b-1*I_NAI_H4yz_Fx2y;
  abcd[1739] = 2.0E0*I_NAI_H3y2z_Hx2y2z_b-1*I_NAI_H3y2z_Fx2y;
  abcd[1740] = 2.0E0*I_NAI_H2y3z_Hx2y2z_b-1*I_NAI_H2y3z_Fx2y;
  abcd[1741] = 2.0E0*I_NAI_Hy4z_Hx2y2z_b-1*I_NAI_Hy4z_Fx2y;
  abcd[1742] = 2.0E0*I_NAI_H5z_Hx2y2z_b-1*I_NAI_H5z_Fx2y;
  abcd[1743] = 2.0E0*I_NAI_H5x_Hxy3z_b-2*I_NAI_H5x_Fxyz;
  abcd[1744] = 2.0E0*I_NAI_H4xy_Hxy3z_b-2*I_NAI_H4xy_Fxyz;
  abcd[1745] = 2.0E0*I_NAI_H4xz_Hxy3z_b-2*I_NAI_H4xz_Fxyz;
  abcd[1746] = 2.0E0*I_NAI_H3x2y_Hxy3z_b-2*I_NAI_H3x2y_Fxyz;
  abcd[1747] = 2.0E0*I_NAI_H3xyz_Hxy3z_b-2*I_NAI_H3xyz_Fxyz;
  abcd[1748] = 2.0E0*I_NAI_H3x2z_Hxy3z_b-2*I_NAI_H3x2z_Fxyz;
  abcd[1749] = 2.0E0*I_NAI_H2x3y_Hxy3z_b-2*I_NAI_H2x3y_Fxyz;
  abcd[1750] = 2.0E0*I_NAI_H2x2yz_Hxy3z_b-2*I_NAI_H2x2yz_Fxyz;
  abcd[1751] = 2.0E0*I_NAI_H2xy2z_Hxy3z_b-2*I_NAI_H2xy2z_Fxyz;
  abcd[1752] = 2.0E0*I_NAI_H2x3z_Hxy3z_b-2*I_NAI_H2x3z_Fxyz;
  abcd[1753] = 2.0E0*I_NAI_Hx4y_Hxy3z_b-2*I_NAI_Hx4y_Fxyz;
  abcd[1754] = 2.0E0*I_NAI_Hx3yz_Hxy3z_b-2*I_NAI_Hx3yz_Fxyz;
  abcd[1755] = 2.0E0*I_NAI_Hx2y2z_Hxy3z_b-2*I_NAI_Hx2y2z_Fxyz;
  abcd[1756] = 2.0E0*I_NAI_Hxy3z_Hxy3z_b-2*I_NAI_Hxy3z_Fxyz;
  abcd[1757] = 2.0E0*I_NAI_Hx4z_Hxy3z_b-2*I_NAI_Hx4z_Fxyz;
  abcd[1758] = 2.0E0*I_NAI_H5y_Hxy3z_b-2*I_NAI_H5y_Fxyz;
  abcd[1759] = 2.0E0*I_NAI_H4yz_Hxy3z_b-2*I_NAI_H4yz_Fxyz;
  abcd[1760] = 2.0E0*I_NAI_H3y2z_Hxy3z_b-2*I_NAI_H3y2z_Fxyz;
  abcd[1761] = 2.0E0*I_NAI_H2y3z_Hxy3z_b-2*I_NAI_H2y3z_Fxyz;
  abcd[1762] = 2.0E0*I_NAI_Hy4z_Hxy3z_b-2*I_NAI_Hy4z_Fxyz;
  abcd[1763] = 2.0E0*I_NAI_H5z_Hxy3z_b-2*I_NAI_H5z_Fxyz;
  abcd[1764] = 2.0E0*I_NAI_H5x_Hx4z_b-3*I_NAI_H5x_Fx2z;
  abcd[1765] = 2.0E0*I_NAI_H4xy_Hx4z_b-3*I_NAI_H4xy_Fx2z;
  abcd[1766] = 2.0E0*I_NAI_H4xz_Hx4z_b-3*I_NAI_H4xz_Fx2z;
  abcd[1767] = 2.0E0*I_NAI_H3x2y_Hx4z_b-3*I_NAI_H3x2y_Fx2z;
  abcd[1768] = 2.0E0*I_NAI_H3xyz_Hx4z_b-3*I_NAI_H3xyz_Fx2z;
  abcd[1769] = 2.0E0*I_NAI_H3x2z_Hx4z_b-3*I_NAI_H3x2z_Fx2z;
  abcd[1770] = 2.0E0*I_NAI_H2x3y_Hx4z_b-3*I_NAI_H2x3y_Fx2z;
  abcd[1771] = 2.0E0*I_NAI_H2x2yz_Hx4z_b-3*I_NAI_H2x2yz_Fx2z;
  abcd[1772] = 2.0E0*I_NAI_H2xy2z_Hx4z_b-3*I_NAI_H2xy2z_Fx2z;
  abcd[1773] = 2.0E0*I_NAI_H2x3z_Hx4z_b-3*I_NAI_H2x3z_Fx2z;
  abcd[1774] = 2.0E0*I_NAI_Hx4y_Hx4z_b-3*I_NAI_Hx4y_Fx2z;
  abcd[1775] = 2.0E0*I_NAI_Hx3yz_Hx4z_b-3*I_NAI_Hx3yz_Fx2z;
  abcd[1776] = 2.0E0*I_NAI_Hx2y2z_Hx4z_b-3*I_NAI_Hx2y2z_Fx2z;
  abcd[1777] = 2.0E0*I_NAI_Hxy3z_Hx4z_b-3*I_NAI_Hxy3z_Fx2z;
  abcd[1778] = 2.0E0*I_NAI_Hx4z_Hx4z_b-3*I_NAI_Hx4z_Fx2z;
  abcd[1779] = 2.0E0*I_NAI_H5y_Hx4z_b-3*I_NAI_H5y_Fx2z;
  abcd[1780] = 2.0E0*I_NAI_H4yz_Hx4z_b-3*I_NAI_H4yz_Fx2z;
  abcd[1781] = 2.0E0*I_NAI_H3y2z_Hx4z_b-3*I_NAI_H3y2z_Fx2z;
  abcd[1782] = 2.0E0*I_NAI_H2y3z_Hx4z_b-3*I_NAI_H2y3z_Fx2z;
  abcd[1783] = 2.0E0*I_NAI_Hy4z_Hx4z_b-3*I_NAI_Hy4z_Fx2z;
  abcd[1784] = 2.0E0*I_NAI_H5z_Hx4z_b-3*I_NAI_H5z_Fx2z;
  abcd[1785] = 2.0E0*I_NAI_H5x_H4yz_b;
  abcd[1786] = 2.0E0*I_NAI_H4xy_H4yz_b;
  abcd[1787] = 2.0E0*I_NAI_H4xz_H4yz_b;
  abcd[1788] = 2.0E0*I_NAI_H3x2y_H4yz_b;
  abcd[1789] = 2.0E0*I_NAI_H3xyz_H4yz_b;
  abcd[1790] = 2.0E0*I_NAI_H3x2z_H4yz_b;
  abcd[1791] = 2.0E0*I_NAI_H2x3y_H4yz_b;
  abcd[1792] = 2.0E0*I_NAI_H2x2yz_H4yz_b;
  abcd[1793] = 2.0E0*I_NAI_H2xy2z_H4yz_b;
  abcd[1794] = 2.0E0*I_NAI_H2x3z_H4yz_b;
  abcd[1795] = 2.0E0*I_NAI_Hx4y_H4yz_b;
  abcd[1796] = 2.0E0*I_NAI_Hx3yz_H4yz_b;
  abcd[1797] = 2.0E0*I_NAI_Hx2y2z_H4yz_b;
  abcd[1798] = 2.0E0*I_NAI_Hxy3z_H4yz_b;
  abcd[1799] = 2.0E0*I_NAI_Hx4z_H4yz_b;
  abcd[1800] = 2.0E0*I_NAI_H5y_H4yz_b;
  abcd[1801] = 2.0E0*I_NAI_H4yz_H4yz_b;
  abcd[1802] = 2.0E0*I_NAI_H3y2z_H4yz_b;
  abcd[1803] = 2.0E0*I_NAI_H2y3z_H4yz_b;
  abcd[1804] = 2.0E0*I_NAI_Hy4z_H4yz_b;
  abcd[1805] = 2.0E0*I_NAI_H5z_H4yz_b;
  abcd[1806] = 2.0E0*I_NAI_H5x_H3y2z_b-1*I_NAI_H5x_F3y;
  abcd[1807] = 2.0E0*I_NAI_H4xy_H3y2z_b-1*I_NAI_H4xy_F3y;
  abcd[1808] = 2.0E0*I_NAI_H4xz_H3y2z_b-1*I_NAI_H4xz_F3y;
  abcd[1809] = 2.0E0*I_NAI_H3x2y_H3y2z_b-1*I_NAI_H3x2y_F3y;
  abcd[1810] = 2.0E0*I_NAI_H3xyz_H3y2z_b-1*I_NAI_H3xyz_F3y;
  abcd[1811] = 2.0E0*I_NAI_H3x2z_H3y2z_b-1*I_NAI_H3x2z_F3y;
  abcd[1812] = 2.0E0*I_NAI_H2x3y_H3y2z_b-1*I_NAI_H2x3y_F3y;
  abcd[1813] = 2.0E0*I_NAI_H2x2yz_H3y2z_b-1*I_NAI_H2x2yz_F3y;
  abcd[1814] = 2.0E0*I_NAI_H2xy2z_H3y2z_b-1*I_NAI_H2xy2z_F3y;
  abcd[1815] = 2.0E0*I_NAI_H2x3z_H3y2z_b-1*I_NAI_H2x3z_F3y;
  abcd[1816] = 2.0E0*I_NAI_Hx4y_H3y2z_b-1*I_NAI_Hx4y_F3y;
  abcd[1817] = 2.0E0*I_NAI_Hx3yz_H3y2z_b-1*I_NAI_Hx3yz_F3y;
  abcd[1818] = 2.0E0*I_NAI_Hx2y2z_H3y2z_b-1*I_NAI_Hx2y2z_F3y;
  abcd[1819] = 2.0E0*I_NAI_Hxy3z_H3y2z_b-1*I_NAI_Hxy3z_F3y;
  abcd[1820] = 2.0E0*I_NAI_Hx4z_H3y2z_b-1*I_NAI_Hx4z_F3y;
  abcd[1821] = 2.0E0*I_NAI_H5y_H3y2z_b-1*I_NAI_H5y_F3y;
  abcd[1822] = 2.0E0*I_NAI_H4yz_H3y2z_b-1*I_NAI_H4yz_F3y;
  abcd[1823] = 2.0E0*I_NAI_H3y2z_H3y2z_b-1*I_NAI_H3y2z_F3y;
  abcd[1824] = 2.0E0*I_NAI_H2y3z_H3y2z_b-1*I_NAI_H2y3z_F3y;
  abcd[1825] = 2.0E0*I_NAI_Hy4z_H3y2z_b-1*I_NAI_Hy4z_F3y;
  abcd[1826] = 2.0E0*I_NAI_H5z_H3y2z_b-1*I_NAI_H5z_F3y;
  abcd[1827] = 2.0E0*I_NAI_H5x_H2y3z_b-2*I_NAI_H5x_F2yz;
  abcd[1828] = 2.0E0*I_NAI_H4xy_H2y3z_b-2*I_NAI_H4xy_F2yz;
  abcd[1829] = 2.0E0*I_NAI_H4xz_H2y3z_b-2*I_NAI_H4xz_F2yz;
  abcd[1830] = 2.0E0*I_NAI_H3x2y_H2y3z_b-2*I_NAI_H3x2y_F2yz;
  abcd[1831] = 2.0E0*I_NAI_H3xyz_H2y3z_b-2*I_NAI_H3xyz_F2yz;
  abcd[1832] = 2.0E0*I_NAI_H3x2z_H2y3z_b-2*I_NAI_H3x2z_F2yz;
  abcd[1833] = 2.0E0*I_NAI_H2x3y_H2y3z_b-2*I_NAI_H2x3y_F2yz;
  abcd[1834] = 2.0E0*I_NAI_H2x2yz_H2y3z_b-2*I_NAI_H2x2yz_F2yz;
  abcd[1835] = 2.0E0*I_NAI_H2xy2z_H2y3z_b-2*I_NAI_H2xy2z_F2yz;
  abcd[1836] = 2.0E0*I_NAI_H2x3z_H2y3z_b-2*I_NAI_H2x3z_F2yz;
  abcd[1837] = 2.0E0*I_NAI_Hx4y_H2y3z_b-2*I_NAI_Hx4y_F2yz;
  abcd[1838] = 2.0E0*I_NAI_Hx3yz_H2y3z_b-2*I_NAI_Hx3yz_F2yz;
  abcd[1839] = 2.0E0*I_NAI_Hx2y2z_H2y3z_b-2*I_NAI_Hx2y2z_F2yz;
  abcd[1840] = 2.0E0*I_NAI_Hxy3z_H2y3z_b-2*I_NAI_Hxy3z_F2yz;
  abcd[1841] = 2.0E0*I_NAI_Hx4z_H2y3z_b-2*I_NAI_Hx4z_F2yz;
  abcd[1842] = 2.0E0*I_NAI_H5y_H2y3z_b-2*I_NAI_H5y_F2yz;
  abcd[1843] = 2.0E0*I_NAI_H4yz_H2y3z_b-2*I_NAI_H4yz_F2yz;
  abcd[1844] = 2.0E0*I_NAI_H3y2z_H2y3z_b-2*I_NAI_H3y2z_F2yz;
  abcd[1845] = 2.0E0*I_NAI_H2y3z_H2y3z_b-2*I_NAI_H2y3z_F2yz;
  abcd[1846] = 2.0E0*I_NAI_Hy4z_H2y3z_b-2*I_NAI_Hy4z_F2yz;
  abcd[1847] = 2.0E0*I_NAI_H5z_H2y3z_b-2*I_NAI_H5z_F2yz;
  abcd[1848] = 2.0E0*I_NAI_H5x_Hy4z_b-3*I_NAI_H5x_Fy2z;
  abcd[1849] = 2.0E0*I_NAI_H4xy_Hy4z_b-3*I_NAI_H4xy_Fy2z;
  abcd[1850] = 2.0E0*I_NAI_H4xz_Hy4z_b-3*I_NAI_H4xz_Fy2z;
  abcd[1851] = 2.0E0*I_NAI_H3x2y_Hy4z_b-3*I_NAI_H3x2y_Fy2z;
  abcd[1852] = 2.0E0*I_NAI_H3xyz_Hy4z_b-3*I_NAI_H3xyz_Fy2z;
  abcd[1853] = 2.0E0*I_NAI_H3x2z_Hy4z_b-3*I_NAI_H3x2z_Fy2z;
  abcd[1854] = 2.0E0*I_NAI_H2x3y_Hy4z_b-3*I_NAI_H2x3y_Fy2z;
  abcd[1855] = 2.0E0*I_NAI_H2x2yz_Hy4z_b-3*I_NAI_H2x2yz_Fy2z;
  abcd[1856] = 2.0E0*I_NAI_H2xy2z_Hy4z_b-3*I_NAI_H2xy2z_Fy2z;
  abcd[1857] = 2.0E0*I_NAI_H2x3z_Hy4z_b-3*I_NAI_H2x3z_Fy2z;
  abcd[1858] = 2.0E0*I_NAI_Hx4y_Hy4z_b-3*I_NAI_Hx4y_Fy2z;
  abcd[1859] = 2.0E0*I_NAI_Hx3yz_Hy4z_b-3*I_NAI_Hx3yz_Fy2z;
  abcd[1860] = 2.0E0*I_NAI_Hx2y2z_Hy4z_b-3*I_NAI_Hx2y2z_Fy2z;
  abcd[1861] = 2.0E0*I_NAI_Hxy3z_Hy4z_b-3*I_NAI_Hxy3z_Fy2z;
  abcd[1862] = 2.0E0*I_NAI_Hx4z_Hy4z_b-3*I_NAI_Hx4z_Fy2z;
  abcd[1863] = 2.0E0*I_NAI_H5y_Hy4z_b-3*I_NAI_H5y_Fy2z;
  abcd[1864] = 2.0E0*I_NAI_H4yz_Hy4z_b-3*I_NAI_H4yz_Fy2z;
  abcd[1865] = 2.0E0*I_NAI_H3y2z_Hy4z_b-3*I_NAI_H3y2z_Fy2z;
  abcd[1866] = 2.0E0*I_NAI_H2y3z_Hy4z_b-3*I_NAI_H2y3z_Fy2z;
  abcd[1867] = 2.0E0*I_NAI_Hy4z_Hy4z_b-3*I_NAI_Hy4z_Fy2z;
  abcd[1868] = 2.0E0*I_NAI_H5z_Hy4z_b-3*I_NAI_H5z_Fy2z;
  abcd[1869] = 2.0E0*I_NAI_H5x_H5z_b-4*I_NAI_H5x_F3z;
  abcd[1870] = 2.0E0*I_NAI_H4xy_H5z_b-4*I_NAI_H4xy_F3z;
  abcd[1871] = 2.0E0*I_NAI_H4xz_H5z_b-4*I_NAI_H4xz_F3z;
  abcd[1872] = 2.0E0*I_NAI_H3x2y_H5z_b-4*I_NAI_H3x2y_F3z;
  abcd[1873] = 2.0E0*I_NAI_H3xyz_H5z_b-4*I_NAI_H3xyz_F3z;
  abcd[1874] = 2.0E0*I_NAI_H3x2z_H5z_b-4*I_NAI_H3x2z_F3z;
  abcd[1875] = 2.0E0*I_NAI_H2x3y_H5z_b-4*I_NAI_H2x3y_F3z;
  abcd[1876] = 2.0E0*I_NAI_H2x2yz_H5z_b-4*I_NAI_H2x2yz_F3z;
  abcd[1877] = 2.0E0*I_NAI_H2xy2z_H5z_b-4*I_NAI_H2xy2z_F3z;
  abcd[1878] = 2.0E0*I_NAI_H2x3z_H5z_b-4*I_NAI_H2x3z_F3z;
  abcd[1879] = 2.0E0*I_NAI_Hx4y_H5z_b-4*I_NAI_Hx4y_F3z;
  abcd[1880] = 2.0E0*I_NAI_Hx3yz_H5z_b-4*I_NAI_Hx3yz_F3z;
  abcd[1881] = 2.0E0*I_NAI_Hx2y2z_H5z_b-4*I_NAI_Hx2y2z_F3z;
  abcd[1882] = 2.0E0*I_NAI_Hxy3z_H5z_b-4*I_NAI_Hxy3z_F3z;
  abcd[1883] = 2.0E0*I_NAI_Hx4z_H5z_b-4*I_NAI_Hx4z_F3z;
  abcd[1884] = 2.0E0*I_NAI_H5y_H5z_b-4*I_NAI_H5y_F3z;
  abcd[1885] = 2.0E0*I_NAI_H4yz_H5z_b-4*I_NAI_H4yz_F3z;
  abcd[1886] = 2.0E0*I_NAI_H3y2z_H5z_b-4*I_NAI_H3y2z_F3z;
  abcd[1887] = 2.0E0*I_NAI_H2y3z_H5z_b-4*I_NAI_H2y3z_F3z;
  abcd[1888] = 2.0E0*I_NAI_Hy4z_H5z_b-4*I_NAI_Hy4z_F3z;
  abcd[1889] = 2.0E0*I_NAI_H5z_H5z_b-4*I_NAI_H5z_F3z;
}
