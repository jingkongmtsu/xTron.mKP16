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

void hgp_os_nai_h_sp_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_I6x_S_C5_a = 0.0E0;
  Double I_NAI_I5xy_S_C5_a = 0.0E0;
  Double I_NAI_I5xz_S_C5_a = 0.0E0;
  Double I_NAI_I4x2y_S_C5_a = 0.0E0;
  Double I_NAI_I4xyz_S_C5_a = 0.0E0;
  Double I_NAI_I4x2z_S_C5_a = 0.0E0;
  Double I_NAI_I3x3y_S_C5_a = 0.0E0;
  Double I_NAI_I3x2yz_S_C5_a = 0.0E0;
  Double I_NAI_I3xy2z_S_C5_a = 0.0E0;
  Double I_NAI_I3x3z_S_C5_a = 0.0E0;
  Double I_NAI_I2x4y_S_C5_a = 0.0E0;
  Double I_NAI_I2x3yz_S_C5_a = 0.0E0;
  Double I_NAI_I2x2y2z_S_C5_a = 0.0E0;
  Double I_NAI_I2xy3z_S_C5_a = 0.0E0;
  Double I_NAI_I2x4z_S_C5_a = 0.0E0;
  Double I_NAI_Ix5y_S_C5_a = 0.0E0;
  Double I_NAI_Ix4yz_S_C5_a = 0.0E0;
  Double I_NAI_Ix3y2z_S_C5_a = 0.0E0;
  Double I_NAI_Ix2y3z_S_C5_a = 0.0E0;
  Double I_NAI_Ixy4z_S_C5_a = 0.0E0;
  Double I_NAI_Ix5z_S_C5_a = 0.0E0;
  Double I_NAI_I6y_S_C5_a = 0.0E0;
  Double I_NAI_I5yz_S_C5_a = 0.0E0;
  Double I_NAI_I4y2z_S_C5_a = 0.0E0;
  Double I_NAI_I3y3z_S_C5_a = 0.0E0;
  Double I_NAI_I2y4z_S_C5_a = 0.0E0;
  Double I_NAI_Iy5z_S_C5_a = 0.0E0;
  Double I_NAI_I6z_S_C5_a = 0.0E0;
  Double I_NAI_G4x_S_C5 = 0.0E0;
  Double I_NAI_G3xy_S_C5 = 0.0E0;
  Double I_NAI_G3xz_S_C5 = 0.0E0;
  Double I_NAI_G2x2y_S_C5 = 0.0E0;
  Double I_NAI_G2xyz_S_C5 = 0.0E0;
  Double I_NAI_G2x2z_S_C5 = 0.0E0;
  Double I_NAI_Gx3y_S_C5 = 0.0E0;
  Double I_NAI_Gx2yz_S_C5 = 0.0E0;
  Double I_NAI_Gxy2z_S_C5 = 0.0E0;
  Double I_NAI_Gx3z_S_C5 = 0.0E0;
  Double I_NAI_G4y_S_C5 = 0.0E0;
  Double I_NAI_G3yz_S_C5 = 0.0E0;
  Double I_NAI_G2y2z_S_C5 = 0.0E0;
  Double I_NAI_Gy3z_S_C5 = 0.0E0;
  Double I_NAI_G4z_S_C5 = 0.0E0;
  Double I_NAI_H5x_S_C1005 = 0.0E0;
  Double I_NAI_H4xy_S_C1005 = 0.0E0;
  Double I_NAI_H4xz_S_C1005 = 0.0E0;
  Double I_NAI_H3x2y_S_C1005 = 0.0E0;
  Double I_NAI_H3xyz_S_C1005 = 0.0E0;
  Double I_NAI_H3x2z_S_C1005 = 0.0E0;
  Double I_NAI_H2x3y_S_C1005 = 0.0E0;
  Double I_NAI_H2x2yz_S_C1005 = 0.0E0;
  Double I_NAI_H2xy2z_S_C1005 = 0.0E0;
  Double I_NAI_H2x3z_S_C1005 = 0.0E0;
  Double I_NAI_Hx4y_S_C1005 = 0.0E0;
  Double I_NAI_Hx3yz_S_C1005 = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1005 = 0.0E0;
  Double I_NAI_Hxy3z_S_C1005 = 0.0E0;
  Double I_NAI_Hx4z_S_C1005 = 0.0E0;
  Double I_NAI_H5y_S_C1005 = 0.0E0;
  Double I_NAI_H4yz_S_C1005 = 0.0E0;
  Double I_NAI_H3y2z_S_C1005 = 0.0E0;
  Double I_NAI_H2y3z_S_C1005 = 0.0E0;
  Double I_NAI_Hy4z_S_C1005 = 0.0E0;
  Double I_NAI_H5z_S_C1005 = 0.0E0;
  Double I_NAI_K7x_S_C1005_a = 0.0E0;
  Double I_NAI_K6xy_S_C1005_a = 0.0E0;
  Double I_NAI_K6xz_S_C1005_a = 0.0E0;
  Double I_NAI_K5x2y_S_C1005_a = 0.0E0;
  Double I_NAI_K5xyz_S_C1005_a = 0.0E0;
  Double I_NAI_K5x2z_S_C1005_a = 0.0E0;
  Double I_NAI_K4x3y_S_C1005_a = 0.0E0;
  Double I_NAI_K4x2yz_S_C1005_a = 0.0E0;
  Double I_NAI_K4xy2z_S_C1005_a = 0.0E0;
  Double I_NAI_K4x3z_S_C1005_a = 0.0E0;
  Double I_NAI_K3x4y_S_C1005_a = 0.0E0;
  Double I_NAI_K3x3yz_S_C1005_a = 0.0E0;
  Double I_NAI_K3x2y2z_S_C1005_a = 0.0E0;
  Double I_NAI_K3xy3z_S_C1005_a = 0.0E0;
  Double I_NAI_K3x4z_S_C1005_a = 0.0E0;
  Double I_NAI_K2x5y_S_C1005_a = 0.0E0;
  Double I_NAI_K2x4yz_S_C1005_a = 0.0E0;
  Double I_NAI_K2x3y2z_S_C1005_a = 0.0E0;
  Double I_NAI_K2x2y3z_S_C1005_a = 0.0E0;
  Double I_NAI_K2xy4z_S_C1005_a = 0.0E0;
  Double I_NAI_K2x5z_S_C1005_a = 0.0E0;
  Double I_NAI_Kx6y_S_C1005_a = 0.0E0;
  Double I_NAI_Kx5yz_S_C1005_a = 0.0E0;
  Double I_NAI_Kx4y2z_S_C1005_a = 0.0E0;
  Double I_NAI_Kx3y3z_S_C1005_a = 0.0E0;
  Double I_NAI_Kx2y4z_S_C1005_a = 0.0E0;
  Double I_NAI_Kxy5z_S_C1005_a = 0.0E0;
  Double I_NAI_Kx6z_S_C1005_a = 0.0E0;
  Double I_NAI_K7y_S_C1005_a = 0.0E0;
  Double I_NAI_K6yz_S_C1005_a = 0.0E0;
  Double I_NAI_K5y2z_S_C1005_a = 0.0E0;
  Double I_NAI_K4y3z_S_C1005_a = 0.0E0;
  Double I_NAI_K3y4z_S_C1005_a = 0.0E0;
  Double I_NAI_K2y5z_S_C1005_a = 0.0E0;
  Double I_NAI_Ky6z_S_C1005_a = 0.0E0;
  Double I_NAI_K7z_S_C1005_a = 0.0E0;
  Double I_NAI_I6x_S_C1005_a = 0.0E0;
  Double I_NAI_I5xy_S_C1005_a = 0.0E0;
  Double I_NAI_I5xz_S_C1005_a = 0.0E0;
  Double I_NAI_I4x2y_S_C1005_a = 0.0E0;
  Double I_NAI_I4xyz_S_C1005_a = 0.0E0;
  Double I_NAI_I4x2z_S_C1005_a = 0.0E0;
  Double I_NAI_I3x3y_S_C1005_a = 0.0E0;
  Double I_NAI_I3x2yz_S_C1005_a = 0.0E0;
  Double I_NAI_I3xy2z_S_C1005_a = 0.0E0;
  Double I_NAI_I3x3z_S_C1005_a = 0.0E0;
  Double I_NAI_I2x4y_S_C1005_a = 0.0E0;
  Double I_NAI_I2x3yz_S_C1005_a = 0.0E0;
  Double I_NAI_I2x2y2z_S_C1005_a = 0.0E0;
  Double I_NAI_I2xy3z_S_C1005_a = 0.0E0;
  Double I_NAI_I2x4z_S_C1005_a = 0.0E0;
  Double I_NAI_Ix5y_S_C1005_a = 0.0E0;
  Double I_NAI_Ix4yz_S_C1005_a = 0.0E0;
  Double I_NAI_Ix3y2z_S_C1005_a = 0.0E0;
  Double I_NAI_Ix2y3z_S_C1005_a = 0.0E0;
  Double I_NAI_Ixy4z_S_C1005_a = 0.0E0;
  Double I_NAI_Ix5z_S_C1005_a = 0.0E0;
  Double I_NAI_I6y_S_C1005_a = 0.0E0;
  Double I_NAI_I5yz_S_C1005_a = 0.0E0;
  Double I_NAI_I4y2z_S_C1005_a = 0.0E0;
  Double I_NAI_I3y3z_S_C1005_a = 0.0E0;
  Double I_NAI_I2y4z_S_C1005_a = 0.0E0;
  Double I_NAI_Iy5z_S_C1005_a = 0.0E0;
  Double I_NAI_I6z_S_C1005_a = 0.0E0;
  Double I_NAI_G4x_S_C1005 = 0.0E0;
  Double I_NAI_G3xy_S_C1005 = 0.0E0;
  Double I_NAI_G3xz_S_C1005 = 0.0E0;
  Double I_NAI_G2x2y_S_C1005 = 0.0E0;
  Double I_NAI_G2xyz_S_C1005 = 0.0E0;
  Double I_NAI_G2x2z_S_C1005 = 0.0E0;
  Double I_NAI_Gx3y_S_C1005 = 0.0E0;
  Double I_NAI_Gx2yz_S_C1005 = 0.0E0;
  Double I_NAI_Gxy2z_S_C1005 = 0.0E0;
  Double I_NAI_Gx3z_S_C1005 = 0.0E0;
  Double I_NAI_G4y_S_C1005 = 0.0E0;
  Double I_NAI_G3yz_S_C1005 = 0.0E0;
  Double I_NAI_G2y2z_S_C1005 = 0.0E0;
  Double I_NAI_Gy3z_S_C1005 = 0.0E0;
  Double I_NAI_G4z_S_C1005 = 0.0E0;
  Double I_NAI_I6x_S_C5_b = 0.0E0;
  Double I_NAI_I5xy_S_C5_b = 0.0E0;
  Double I_NAI_I5xz_S_C5_b = 0.0E0;
  Double I_NAI_I4x2y_S_C5_b = 0.0E0;
  Double I_NAI_I4xyz_S_C5_b = 0.0E0;
  Double I_NAI_I4x2z_S_C5_b = 0.0E0;
  Double I_NAI_I3x3y_S_C5_b = 0.0E0;
  Double I_NAI_I3x2yz_S_C5_b = 0.0E0;
  Double I_NAI_I3xy2z_S_C5_b = 0.0E0;
  Double I_NAI_I3x3z_S_C5_b = 0.0E0;
  Double I_NAI_I2x4y_S_C5_b = 0.0E0;
  Double I_NAI_I2x3yz_S_C5_b = 0.0E0;
  Double I_NAI_I2x2y2z_S_C5_b = 0.0E0;
  Double I_NAI_I2xy3z_S_C5_b = 0.0E0;
  Double I_NAI_I2x4z_S_C5_b = 0.0E0;
  Double I_NAI_Ix5y_S_C5_b = 0.0E0;
  Double I_NAI_Ix4yz_S_C5_b = 0.0E0;
  Double I_NAI_Ix3y2z_S_C5_b = 0.0E0;
  Double I_NAI_Ix2y3z_S_C5_b = 0.0E0;
  Double I_NAI_Ixy4z_S_C5_b = 0.0E0;
  Double I_NAI_Ix5z_S_C5_b = 0.0E0;
  Double I_NAI_I6y_S_C5_b = 0.0E0;
  Double I_NAI_I5yz_S_C5_b = 0.0E0;
  Double I_NAI_I4y2z_S_C5_b = 0.0E0;
  Double I_NAI_I3y3z_S_C5_b = 0.0E0;
  Double I_NAI_I2y4z_S_C5_b = 0.0E0;
  Double I_NAI_Iy5z_S_C5_b = 0.0E0;
  Double I_NAI_I6z_S_C5_b = 0.0E0;
  Double I_NAI_H5x_S_C5_b = 0.0E0;
  Double I_NAI_H4xy_S_C5_b = 0.0E0;
  Double I_NAI_H4xz_S_C5_b = 0.0E0;
  Double I_NAI_H3x2y_S_C5_b = 0.0E0;
  Double I_NAI_H3xyz_S_C5_b = 0.0E0;
  Double I_NAI_H3x2z_S_C5_b = 0.0E0;
  Double I_NAI_H2x3y_S_C5_b = 0.0E0;
  Double I_NAI_H2x2yz_S_C5_b = 0.0E0;
  Double I_NAI_H2xy2z_S_C5_b = 0.0E0;
  Double I_NAI_H2x3z_S_C5_b = 0.0E0;
  Double I_NAI_Hx4y_S_C5_b = 0.0E0;
  Double I_NAI_Hx3yz_S_C5_b = 0.0E0;
  Double I_NAI_Hx2y2z_S_C5_b = 0.0E0;
  Double I_NAI_Hxy3z_S_C5_b = 0.0E0;
  Double I_NAI_Hx4z_S_C5_b = 0.0E0;
  Double I_NAI_H5y_S_C5_b = 0.0E0;
  Double I_NAI_H4yz_S_C5_b = 0.0E0;
  Double I_NAI_H3y2z_S_C5_b = 0.0E0;
  Double I_NAI_H2y3z_S_C5_b = 0.0E0;
  Double I_NAI_Hy4z_S_C5_b = 0.0E0;
  Double I_NAI_H5z_S_C5_b = 0.0E0;
  Double I_NAI_K7x_S_C1005_b = 0.0E0;
  Double I_NAI_K6xy_S_C1005_b = 0.0E0;
  Double I_NAI_K6xz_S_C1005_b = 0.0E0;
  Double I_NAI_K5x2y_S_C1005_b = 0.0E0;
  Double I_NAI_K5xyz_S_C1005_b = 0.0E0;
  Double I_NAI_K5x2z_S_C1005_b = 0.0E0;
  Double I_NAI_K4x3y_S_C1005_b = 0.0E0;
  Double I_NAI_K4x2yz_S_C1005_b = 0.0E0;
  Double I_NAI_K4xy2z_S_C1005_b = 0.0E0;
  Double I_NAI_K4x3z_S_C1005_b = 0.0E0;
  Double I_NAI_K3x4y_S_C1005_b = 0.0E0;
  Double I_NAI_K3x3yz_S_C1005_b = 0.0E0;
  Double I_NAI_K3x2y2z_S_C1005_b = 0.0E0;
  Double I_NAI_K3xy3z_S_C1005_b = 0.0E0;
  Double I_NAI_K3x4z_S_C1005_b = 0.0E0;
  Double I_NAI_K2x5y_S_C1005_b = 0.0E0;
  Double I_NAI_K2x4yz_S_C1005_b = 0.0E0;
  Double I_NAI_K2x3y2z_S_C1005_b = 0.0E0;
  Double I_NAI_K2x2y3z_S_C1005_b = 0.0E0;
  Double I_NAI_K2xy4z_S_C1005_b = 0.0E0;
  Double I_NAI_K2x5z_S_C1005_b = 0.0E0;
  Double I_NAI_Kx6y_S_C1005_b = 0.0E0;
  Double I_NAI_Kx5yz_S_C1005_b = 0.0E0;
  Double I_NAI_Kx4y2z_S_C1005_b = 0.0E0;
  Double I_NAI_Kx3y3z_S_C1005_b = 0.0E0;
  Double I_NAI_Kx2y4z_S_C1005_b = 0.0E0;
  Double I_NAI_Kxy5z_S_C1005_b = 0.0E0;
  Double I_NAI_Kx6z_S_C1005_b = 0.0E0;
  Double I_NAI_K7y_S_C1005_b = 0.0E0;
  Double I_NAI_K6yz_S_C1005_b = 0.0E0;
  Double I_NAI_K5y2z_S_C1005_b = 0.0E0;
  Double I_NAI_K4y3z_S_C1005_b = 0.0E0;
  Double I_NAI_K3y4z_S_C1005_b = 0.0E0;
  Double I_NAI_K2y5z_S_C1005_b = 0.0E0;
  Double I_NAI_Ky6z_S_C1005_b = 0.0E0;
  Double I_NAI_K7z_S_C1005_b = 0.0E0;
  Double I_NAI_I6x_S_C1005_b = 0.0E0;
  Double I_NAI_I5xy_S_C1005_b = 0.0E0;
  Double I_NAI_I5xz_S_C1005_b = 0.0E0;
  Double I_NAI_I4x2y_S_C1005_b = 0.0E0;
  Double I_NAI_I4xyz_S_C1005_b = 0.0E0;
  Double I_NAI_I4x2z_S_C1005_b = 0.0E0;
  Double I_NAI_I3x3y_S_C1005_b = 0.0E0;
  Double I_NAI_I3x2yz_S_C1005_b = 0.0E0;
  Double I_NAI_I3xy2z_S_C1005_b = 0.0E0;
  Double I_NAI_I3x3z_S_C1005_b = 0.0E0;
  Double I_NAI_I2x4y_S_C1005_b = 0.0E0;
  Double I_NAI_I2x3yz_S_C1005_b = 0.0E0;
  Double I_NAI_I2x2y2z_S_C1005_b = 0.0E0;
  Double I_NAI_I2xy3z_S_C1005_b = 0.0E0;
  Double I_NAI_I2x4z_S_C1005_b = 0.0E0;
  Double I_NAI_Ix5y_S_C1005_b = 0.0E0;
  Double I_NAI_Ix4yz_S_C1005_b = 0.0E0;
  Double I_NAI_Ix3y2z_S_C1005_b = 0.0E0;
  Double I_NAI_Ix2y3z_S_C1005_b = 0.0E0;
  Double I_NAI_Ixy4z_S_C1005_b = 0.0E0;
  Double I_NAI_Ix5z_S_C1005_b = 0.0E0;
  Double I_NAI_I6y_S_C1005_b = 0.0E0;
  Double I_NAI_I5yz_S_C1005_b = 0.0E0;
  Double I_NAI_I4y2z_S_C1005_b = 0.0E0;
  Double I_NAI_I3y3z_S_C1005_b = 0.0E0;
  Double I_NAI_I2y4z_S_C1005_b = 0.0E0;
  Double I_NAI_Iy5z_S_C1005_b = 0.0E0;
  Double I_NAI_I6z_S_C1005_b = 0.0E0;
  Double I_NAI_H5x_S_C1005_b = 0.0E0;
  Double I_NAI_H4xy_S_C1005_b = 0.0E0;
  Double I_NAI_H4xz_S_C1005_b = 0.0E0;
  Double I_NAI_H3x2y_S_C1005_b = 0.0E0;
  Double I_NAI_H3xyz_S_C1005_b = 0.0E0;
  Double I_NAI_H3x2z_S_C1005_b = 0.0E0;
  Double I_NAI_H2x3y_S_C1005_b = 0.0E0;
  Double I_NAI_H2x2yz_S_C1005_b = 0.0E0;
  Double I_NAI_H2xy2z_S_C1005_b = 0.0E0;
  Double I_NAI_H2x3z_S_C1005_b = 0.0E0;
  Double I_NAI_Hx4y_S_C1005_b = 0.0E0;
  Double I_NAI_Hx3yz_S_C1005_b = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1005_b = 0.0E0;
  Double I_NAI_Hxy3z_S_C1005_b = 0.0E0;
  Double I_NAI_Hx4z_S_C1005_b = 0.0E0;
  Double I_NAI_H5y_S_C1005_b = 0.0E0;
  Double I_NAI_H4yz_S_C1005_b = 0.0E0;
  Double I_NAI_H3y2z_S_C1005_b = 0.0E0;
  Double I_NAI_H2y3z_S_C1005_b = 0.0E0;
  Double I_NAI_Hy4z_S_C1005_b = 0.0E0;
  Double I_NAI_H5z_S_C1005_b = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
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
      Double prefactor = -charge*fbra;

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
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER49;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER47*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = ONEOVER15*I_NAI_S_S_M7_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M7_vrr  = f*I_NAI_S_S_M7_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);
        I_NAI_S_S_M5_vrr = static_cast<Double>(I_NAI_S_S_M5_vrr_d);
        I_NAI_S_S_M6_vrr = static_cast<Double>(I_NAI_S_S_M6_vrr_d);
        I_NAI_S_S_M7_vrr = static_cast<Double>(I_NAI_S_S_M7_vrr_d);

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

#endif

      }


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
       * shell quartet name: SQ_NAI_I_S_C5_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C5_a_coefs = ic2*alpha;
      I_NAI_I6x_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C5_a += SQ_NAI_I_S_C5_a_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C5
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C5_coefs = ic2;
      I_NAI_G4x_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C5 += SQ_NAI_G_S_C5_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1005
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1005_coefs = ic2_1;
      I_NAI_H5x_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1005 += SQ_NAI_H_S_C1005_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_C1005_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_C1005_a_coefs = ic2_1*alpha;
      I_NAI_K7x_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_C1005_a += SQ_NAI_K_S_C1005_a_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C1005_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C1005_a_coefs = ic2_1*alpha;
      I_NAI_I6x_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C1005_a += SQ_NAI_I_S_C1005_a_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1005
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1005_coefs = ic2_1;
      I_NAI_G4x_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1005 += SQ_NAI_G_S_C1005_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C5_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C5_b_coefs = ic2*beta;
      I_NAI_I6x_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C5_b += SQ_NAI_I_S_C5_b_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C5_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C5_b_coefs = ic2*beta;
      I_NAI_H5x_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C5_b += SQ_NAI_H_S_C5_b_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_C1005_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_C1005_b_coefs = ic2_1*beta;
      I_NAI_K7x_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_C1005_b += SQ_NAI_K_S_C1005_b_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C1005_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C1005_b_coefs = ic2_1*beta;
      I_NAI_I6x_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C1005_b += SQ_NAI_I_S_C1005_b_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1005_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1005_b_coefs = ic2_1*beta;
      I_NAI_H5x_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1005_b += SQ_NAI_H_S_C1005_b_coefs*I_NAI_H5z_S_vrr;
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
   * shell quartet name: SQ_NAI_G_P_C1005
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C1005
   * RHS shell quartet name: SQ_NAI_G_S_C1005
   ************************************************************/
  Double I_NAI_G4x_Px_C1005 = I_NAI_H5x_S_C1005+ABX*I_NAI_G4x_S_C1005;
  Double I_NAI_G3xy_Px_C1005 = I_NAI_H4xy_S_C1005+ABX*I_NAI_G3xy_S_C1005;
  Double I_NAI_G3xz_Px_C1005 = I_NAI_H4xz_S_C1005+ABX*I_NAI_G3xz_S_C1005;
  Double I_NAI_G2x2y_Px_C1005 = I_NAI_H3x2y_S_C1005+ABX*I_NAI_G2x2y_S_C1005;
  Double I_NAI_G2xyz_Px_C1005 = I_NAI_H3xyz_S_C1005+ABX*I_NAI_G2xyz_S_C1005;
  Double I_NAI_G2x2z_Px_C1005 = I_NAI_H3x2z_S_C1005+ABX*I_NAI_G2x2z_S_C1005;
  Double I_NAI_Gx3y_Px_C1005 = I_NAI_H2x3y_S_C1005+ABX*I_NAI_Gx3y_S_C1005;
  Double I_NAI_Gx2yz_Px_C1005 = I_NAI_H2x2yz_S_C1005+ABX*I_NAI_Gx2yz_S_C1005;
  Double I_NAI_Gxy2z_Px_C1005 = I_NAI_H2xy2z_S_C1005+ABX*I_NAI_Gxy2z_S_C1005;
  Double I_NAI_Gx3z_Px_C1005 = I_NAI_H2x3z_S_C1005+ABX*I_NAI_Gx3z_S_C1005;
  Double I_NAI_G4y_Px_C1005 = I_NAI_Hx4y_S_C1005+ABX*I_NAI_G4y_S_C1005;
  Double I_NAI_G3yz_Px_C1005 = I_NAI_Hx3yz_S_C1005+ABX*I_NAI_G3yz_S_C1005;
  Double I_NAI_G2y2z_Px_C1005 = I_NAI_Hx2y2z_S_C1005+ABX*I_NAI_G2y2z_S_C1005;
  Double I_NAI_Gy3z_Px_C1005 = I_NAI_Hxy3z_S_C1005+ABX*I_NAI_Gy3z_S_C1005;
  Double I_NAI_G4z_Px_C1005 = I_NAI_Hx4z_S_C1005+ABX*I_NAI_G4z_S_C1005;
  Double I_NAI_G4x_Py_C1005 = I_NAI_H4xy_S_C1005+ABY*I_NAI_G4x_S_C1005;
  Double I_NAI_G3xy_Py_C1005 = I_NAI_H3x2y_S_C1005+ABY*I_NAI_G3xy_S_C1005;
  Double I_NAI_G3xz_Py_C1005 = I_NAI_H3xyz_S_C1005+ABY*I_NAI_G3xz_S_C1005;
  Double I_NAI_G2x2y_Py_C1005 = I_NAI_H2x3y_S_C1005+ABY*I_NAI_G2x2y_S_C1005;
  Double I_NAI_G2xyz_Py_C1005 = I_NAI_H2x2yz_S_C1005+ABY*I_NAI_G2xyz_S_C1005;
  Double I_NAI_G2x2z_Py_C1005 = I_NAI_H2xy2z_S_C1005+ABY*I_NAI_G2x2z_S_C1005;
  Double I_NAI_Gx3y_Py_C1005 = I_NAI_Hx4y_S_C1005+ABY*I_NAI_Gx3y_S_C1005;
  Double I_NAI_Gx2yz_Py_C1005 = I_NAI_Hx3yz_S_C1005+ABY*I_NAI_Gx2yz_S_C1005;
  Double I_NAI_Gxy2z_Py_C1005 = I_NAI_Hx2y2z_S_C1005+ABY*I_NAI_Gxy2z_S_C1005;
  Double I_NAI_Gx3z_Py_C1005 = I_NAI_Hxy3z_S_C1005+ABY*I_NAI_Gx3z_S_C1005;
  Double I_NAI_G4y_Py_C1005 = I_NAI_H5y_S_C1005+ABY*I_NAI_G4y_S_C1005;
  Double I_NAI_G3yz_Py_C1005 = I_NAI_H4yz_S_C1005+ABY*I_NAI_G3yz_S_C1005;
  Double I_NAI_G2y2z_Py_C1005 = I_NAI_H3y2z_S_C1005+ABY*I_NAI_G2y2z_S_C1005;
  Double I_NAI_Gy3z_Py_C1005 = I_NAI_H2y3z_S_C1005+ABY*I_NAI_Gy3z_S_C1005;
  Double I_NAI_G4z_Py_C1005 = I_NAI_Hy4z_S_C1005+ABY*I_NAI_G4z_S_C1005;
  Double I_NAI_G4x_Pz_C1005 = I_NAI_H4xz_S_C1005+ABZ*I_NAI_G4x_S_C1005;
  Double I_NAI_G3xy_Pz_C1005 = I_NAI_H3xyz_S_C1005+ABZ*I_NAI_G3xy_S_C1005;
  Double I_NAI_G3xz_Pz_C1005 = I_NAI_H3x2z_S_C1005+ABZ*I_NAI_G3xz_S_C1005;
  Double I_NAI_G2x2y_Pz_C1005 = I_NAI_H2x2yz_S_C1005+ABZ*I_NAI_G2x2y_S_C1005;
  Double I_NAI_G2xyz_Pz_C1005 = I_NAI_H2xy2z_S_C1005+ABZ*I_NAI_G2xyz_S_C1005;
  Double I_NAI_G2x2z_Pz_C1005 = I_NAI_H2x3z_S_C1005+ABZ*I_NAI_G2x2z_S_C1005;
  Double I_NAI_Gx3y_Pz_C1005 = I_NAI_Hx3yz_S_C1005+ABZ*I_NAI_Gx3y_S_C1005;
  Double I_NAI_Gx2yz_Pz_C1005 = I_NAI_Hx2y2z_S_C1005+ABZ*I_NAI_Gx2yz_S_C1005;
  Double I_NAI_Gxy2z_Pz_C1005 = I_NAI_Hxy3z_S_C1005+ABZ*I_NAI_Gxy2z_S_C1005;
  Double I_NAI_Gx3z_Pz_C1005 = I_NAI_Hx4z_S_C1005+ABZ*I_NAI_Gx3z_S_C1005;
  Double I_NAI_G4y_Pz_C1005 = I_NAI_H4yz_S_C1005+ABZ*I_NAI_G4y_S_C1005;
  Double I_NAI_G3yz_Pz_C1005 = I_NAI_H3y2z_S_C1005+ABZ*I_NAI_G3yz_S_C1005;
  Double I_NAI_G2y2z_Pz_C1005 = I_NAI_H2y3z_S_C1005+ABZ*I_NAI_G2y2z_S_C1005;
  Double I_NAI_Gy3z_Pz_C1005 = I_NAI_Hy4z_S_C1005+ABZ*I_NAI_Gy3z_S_C1005;
  Double I_NAI_G4z_Pz_C1005 = I_NAI_H5z_S_C1005+ABZ*I_NAI_G4z_S_C1005;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_C1005_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_C1005_a
   * RHS shell quartet name: SQ_NAI_I_S_C1005_a
   ************************************************************/
  Double I_NAI_I6x_Px_C1005_a = I_NAI_K7x_S_C1005_a+ABX*I_NAI_I6x_S_C1005_a;
  Double I_NAI_I5xy_Px_C1005_a = I_NAI_K6xy_S_C1005_a+ABX*I_NAI_I5xy_S_C1005_a;
  Double I_NAI_I5xz_Px_C1005_a = I_NAI_K6xz_S_C1005_a+ABX*I_NAI_I5xz_S_C1005_a;
  Double I_NAI_I4x2y_Px_C1005_a = I_NAI_K5x2y_S_C1005_a+ABX*I_NAI_I4x2y_S_C1005_a;
  Double I_NAI_I4xyz_Px_C1005_a = I_NAI_K5xyz_S_C1005_a+ABX*I_NAI_I4xyz_S_C1005_a;
  Double I_NAI_I4x2z_Px_C1005_a = I_NAI_K5x2z_S_C1005_a+ABX*I_NAI_I4x2z_S_C1005_a;
  Double I_NAI_I3x3y_Px_C1005_a = I_NAI_K4x3y_S_C1005_a+ABX*I_NAI_I3x3y_S_C1005_a;
  Double I_NAI_I3x2yz_Px_C1005_a = I_NAI_K4x2yz_S_C1005_a+ABX*I_NAI_I3x2yz_S_C1005_a;
  Double I_NAI_I3xy2z_Px_C1005_a = I_NAI_K4xy2z_S_C1005_a+ABX*I_NAI_I3xy2z_S_C1005_a;
  Double I_NAI_I3x3z_Px_C1005_a = I_NAI_K4x3z_S_C1005_a+ABX*I_NAI_I3x3z_S_C1005_a;
  Double I_NAI_I2x4y_Px_C1005_a = I_NAI_K3x4y_S_C1005_a+ABX*I_NAI_I2x4y_S_C1005_a;
  Double I_NAI_I2x3yz_Px_C1005_a = I_NAI_K3x3yz_S_C1005_a+ABX*I_NAI_I2x3yz_S_C1005_a;
  Double I_NAI_I2x2y2z_Px_C1005_a = I_NAI_K3x2y2z_S_C1005_a+ABX*I_NAI_I2x2y2z_S_C1005_a;
  Double I_NAI_I2xy3z_Px_C1005_a = I_NAI_K3xy3z_S_C1005_a+ABX*I_NAI_I2xy3z_S_C1005_a;
  Double I_NAI_I2x4z_Px_C1005_a = I_NAI_K3x4z_S_C1005_a+ABX*I_NAI_I2x4z_S_C1005_a;
  Double I_NAI_Ix5y_Px_C1005_a = I_NAI_K2x5y_S_C1005_a+ABX*I_NAI_Ix5y_S_C1005_a;
  Double I_NAI_Ix4yz_Px_C1005_a = I_NAI_K2x4yz_S_C1005_a+ABX*I_NAI_Ix4yz_S_C1005_a;
  Double I_NAI_Ix3y2z_Px_C1005_a = I_NAI_K2x3y2z_S_C1005_a+ABX*I_NAI_Ix3y2z_S_C1005_a;
  Double I_NAI_Ix2y3z_Px_C1005_a = I_NAI_K2x2y3z_S_C1005_a+ABX*I_NAI_Ix2y3z_S_C1005_a;
  Double I_NAI_Ixy4z_Px_C1005_a = I_NAI_K2xy4z_S_C1005_a+ABX*I_NAI_Ixy4z_S_C1005_a;
  Double I_NAI_Ix5z_Px_C1005_a = I_NAI_K2x5z_S_C1005_a+ABX*I_NAI_Ix5z_S_C1005_a;
  Double I_NAI_I6y_Px_C1005_a = I_NAI_Kx6y_S_C1005_a+ABX*I_NAI_I6y_S_C1005_a;
  Double I_NAI_I5yz_Px_C1005_a = I_NAI_Kx5yz_S_C1005_a+ABX*I_NAI_I5yz_S_C1005_a;
  Double I_NAI_I4y2z_Px_C1005_a = I_NAI_Kx4y2z_S_C1005_a+ABX*I_NAI_I4y2z_S_C1005_a;
  Double I_NAI_I3y3z_Px_C1005_a = I_NAI_Kx3y3z_S_C1005_a+ABX*I_NAI_I3y3z_S_C1005_a;
  Double I_NAI_I2y4z_Px_C1005_a = I_NAI_Kx2y4z_S_C1005_a+ABX*I_NAI_I2y4z_S_C1005_a;
  Double I_NAI_Iy5z_Px_C1005_a = I_NAI_Kxy5z_S_C1005_a+ABX*I_NAI_Iy5z_S_C1005_a;
  Double I_NAI_I6z_Px_C1005_a = I_NAI_Kx6z_S_C1005_a+ABX*I_NAI_I6z_S_C1005_a;
  Double I_NAI_I6x_Py_C1005_a = I_NAI_K6xy_S_C1005_a+ABY*I_NAI_I6x_S_C1005_a;
  Double I_NAI_I5xy_Py_C1005_a = I_NAI_K5x2y_S_C1005_a+ABY*I_NAI_I5xy_S_C1005_a;
  Double I_NAI_I5xz_Py_C1005_a = I_NAI_K5xyz_S_C1005_a+ABY*I_NAI_I5xz_S_C1005_a;
  Double I_NAI_I4x2y_Py_C1005_a = I_NAI_K4x3y_S_C1005_a+ABY*I_NAI_I4x2y_S_C1005_a;
  Double I_NAI_I4xyz_Py_C1005_a = I_NAI_K4x2yz_S_C1005_a+ABY*I_NAI_I4xyz_S_C1005_a;
  Double I_NAI_I4x2z_Py_C1005_a = I_NAI_K4xy2z_S_C1005_a+ABY*I_NAI_I4x2z_S_C1005_a;
  Double I_NAI_I3x3y_Py_C1005_a = I_NAI_K3x4y_S_C1005_a+ABY*I_NAI_I3x3y_S_C1005_a;
  Double I_NAI_I3x2yz_Py_C1005_a = I_NAI_K3x3yz_S_C1005_a+ABY*I_NAI_I3x2yz_S_C1005_a;
  Double I_NAI_I3xy2z_Py_C1005_a = I_NAI_K3x2y2z_S_C1005_a+ABY*I_NAI_I3xy2z_S_C1005_a;
  Double I_NAI_I3x3z_Py_C1005_a = I_NAI_K3xy3z_S_C1005_a+ABY*I_NAI_I3x3z_S_C1005_a;
  Double I_NAI_I2x4y_Py_C1005_a = I_NAI_K2x5y_S_C1005_a+ABY*I_NAI_I2x4y_S_C1005_a;
  Double I_NAI_I2x3yz_Py_C1005_a = I_NAI_K2x4yz_S_C1005_a+ABY*I_NAI_I2x3yz_S_C1005_a;
  Double I_NAI_I2x2y2z_Py_C1005_a = I_NAI_K2x3y2z_S_C1005_a+ABY*I_NAI_I2x2y2z_S_C1005_a;
  Double I_NAI_I2xy3z_Py_C1005_a = I_NAI_K2x2y3z_S_C1005_a+ABY*I_NAI_I2xy3z_S_C1005_a;
  Double I_NAI_I2x4z_Py_C1005_a = I_NAI_K2xy4z_S_C1005_a+ABY*I_NAI_I2x4z_S_C1005_a;
  Double I_NAI_Ix5y_Py_C1005_a = I_NAI_Kx6y_S_C1005_a+ABY*I_NAI_Ix5y_S_C1005_a;
  Double I_NAI_Ix4yz_Py_C1005_a = I_NAI_Kx5yz_S_C1005_a+ABY*I_NAI_Ix4yz_S_C1005_a;
  Double I_NAI_Ix3y2z_Py_C1005_a = I_NAI_Kx4y2z_S_C1005_a+ABY*I_NAI_Ix3y2z_S_C1005_a;
  Double I_NAI_Ix2y3z_Py_C1005_a = I_NAI_Kx3y3z_S_C1005_a+ABY*I_NAI_Ix2y3z_S_C1005_a;
  Double I_NAI_Ixy4z_Py_C1005_a = I_NAI_Kx2y4z_S_C1005_a+ABY*I_NAI_Ixy4z_S_C1005_a;
  Double I_NAI_Ix5z_Py_C1005_a = I_NAI_Kxy5z_S_C1005_a+ABY*I_NAI_Ix5z_S_C1005_a;
  Double I_NAI_I6y_Py_C1005_a = I_NAI_K7y_S_C1005_a+ABY*I_NAI_I6y_S_C1005_a;
  Double I_NAI_I5yz_Py_C1005_a = I_NAI_K6yz_S_C1005_a+ABY*I_NAI_I5yz_S_C1005_a;
  Double I_NAI_I4y2z_Py_C1005_a = I_NAI_K5y2z_S_C1005_a+ABY*I_NAI_I4y2z_S_C1005_a;
  Double I_NAI_I3y3z_Py_C1005_a = I_NAI_K4y3z_S_C1005_a+ABY*I_NAI_I3y3z_S_C1005_a;
  Double I_NAI_I2y4z_Py_C1005_a = I_NAI_K3y4z_S_C1005_a+ABY*I_NAI_I2y4z_S_C1005_a;
  Double I_NAI_Iy5z_Py_C1005_a = I_NAI_K2y5z_S_C1005_a+ABY*I_NAI_Iy5z_S_C1005_a;
  Double I_NAI_I6z_Py_C1005_a = I_NAI_Ky6z_S_C1005_a+ABY*I_NAI_I6z_S_C1005_a;
  Double I_NAI_I6x_Pz_C1005_a = I_NAI_K6xz_S_C1005_a+ABZ*I_NAI_I6x_S_C1005_a;
  Double I_NAI_I5xy_Pz_C1005_a = I_NAI_K5xyz_S_C1005_a+ABZ*I_NAI_I5xy_S_C1005_a;
  Double I_NAI_I5xz_Pz_C1005_a = I_NAI_K5x2z_S_C1005_a+ABZ*I_NAI_I5xz_S_C1005_a;
  Double I_NAI_I4x2y_Pz_C1005_a = I_NAI_K4x2yz_S_C1005_a+ABZ*I_NAI_I4x2y_S_C1005_a;
  Double I_NAI_I4xyz_Pz_C1005_a = I_NAI_K4xy2z_S_C1005_a+ABZ*I_NAI_I4xyz_S_C1005_a;
  Double I_NAI_I4x2z_Pz_C1005_a = I_NAI_K4x3z_S_C1005_a+ABZ*I_NAI_I4x2z_S_C1005_a;
  Double I_NAI_I3x3y_Pz_C1005_a = I_NAI_K3x3yz_S_C1005_a+ABZ*I_NAI_I3x3y_S_C1005_a;
  Double I_NAI_I3x2yz_Pz_C1005_a = I_NAI_K3x2y2z_S_C1005_a+ABZ*I_NAI_I3x2yz_S_C1005_a;
  Double I_NAI_I3xy2z_Pz_C1005_a = I_NAI_K3xy3z_S_C1005_a+ABZ*I_NAI_I3xy2z_S_C1005_a;
  Double I_NAI_I3x3z_Pz_C1005_a = I_NAI_K3x4z_S_C1005_a+ABZ*I_NAI_I3x3z_S_C1005_a;
  Double I_NAI_I2x4y_Pz_C1005_a = I_NAI_K2x4yz_S_C1005_a+ABZ*I_NAI_I2x4y_S_C1005_a;
  Double I_NAI_I2x3yz_Pz_C1005_a = I_NAI_K2x3y2z_S_C1005_a+ABZ*I_NAI_I2x3yz_S_C1005_a;
  Double I_NAI_I2x2y2z_Pz_C1005_a = I_NAI_K2x2y3z_S_C1005_a+ABZ*I_NAI_I2x2y2z_S_C1005_a;
  Double I_NAI_I2xy3z_Pz_C1005_a = I_NAI_K2xy4z_S_C1005_a+ABZ*I_NAI_I2xy3z_S_C1005_a;
  Double I_NAI_I2x4z_Pz_C1005_a = I_NAI_K2x5z_S_C1005_a+ABZ*I_NAI_I2x4z_S_C1005_a;
  Double I_NAI_Ix5y_Pz_C1005_a = I_NAI_Kx5yz_S_C1005_a+ABZ*I_NAI_Ix5y_S_C1005_a;
  Double I_NAI_Ix4yz_Pz_C1005_a = I_NAI_Kx4y2z_S_C1005_a+ABZ*I_NAI_Ix4yz_S_C1005_a;
  Double I_NAI_Ix3y2z_Pz_C1005_a = I_NAI_Kx3y3z_S_C1005_a+ABZ*I_NAI_Ix3y2z_S_C1005_a;
  Double I_NAI_Ix2y3z_Pz_C1005_a = I_NAI_Kx2y4z_S_C1005_a+ABZ*I_NAI_Ix2y3z_S_C1005_a;
  Double I_NAI_Ixy4z_Pz_C1005_a = I_NAI_Kxy5z_S_C1005_a+ABZ*I_NAI_Ixy4z_S_C1005_a;
  Double I_NAI_Ix5z_Pz_C1005_a = I_NAI_Kx6z_S_C1005_a+ABZ*I_NAI_Ix5z_S_C1005_a;
  Double I_NAI_I6y_Pz_C1005_a = I_NAI_K6yz_S_C1005_a+ABZ*I_NAI_I6y_S_C1005_a;
  Double I_NAI_I5yz_Pz_C1005_a = I_NAI_K5y2z_S_C1005_a+ABZ*I_NAI_I5yz_S_C1005_a;
  Double I_NAI_I4y2z_Pz_C1005_a = I_NAI_K4y3z_S_C1005_a+ABZ*I_NAI_I4y2z_S_C1005_a;
  Double I_NAI_I3y3z_Pz_C1005_a = I_NAI_K3y4z_S_C1005_a+ABZ*I_NAI_I3y3z_S_C1005_a;
  Double I_NAI_I2y4z_Pz_C1005_a = I_NAI_K2y5z_S_C1005_a+ABZ*I_NAI_I2y4z_S_C1005_a;
  Double I_NAI_Iy5z_Pz_C1005_a = I_NAI_Ky6z_S_C1005_a+ABZ*I_NAI_Iy5z_S_C1005_a;
  Double I_NAI_I6z_Pz_C1005_a = I_NAI_K7z_S_C1005_a+ABZ*I_NAI_I6z_S_C1005_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C5_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C5_b
   * RHS shell quartet name: SQ_NAI_H_S_C5_b
   ************************************************************/
  Double I_NAI_H5x_Px_C5_b = I_NAI_I6x_S_C5_b+ABX*I_NAI_H5x_S_C5_b;
  Double I_NAI_H4xy_Px_C5_b = I_NAI_I5xy_S_C5_b+ABX*I_NAI_H4xy_S_C5_b;
  Double I_NAI_H4xz_Px_C5_b = I_NAI_I5xz_S_C5_b+ABX*I_NAI_H4xz_S_C5_b;
  Double I_NAI_H3x2y_Px_C5_b = I_NAI_I4x2y_S_C5_b+ABX*I_NAI_H3x2y_S_C5_b;
  Double I_NAI_H3xyz_Px_C5_b = I_NAI_I4xyz_S_C5_b+ABX*I_NAI_H3xyz_S_C5_b;
  Double I_NAI_H3x2z_Px_C5_b = I_NAI_I4x2z_S_C5_b+ABX*I_NAI_H3x2z_S_C5_b;
  Double I_NAI_H2x3y_Px_C5_b = I_NAI_I3x3y_S_C5_b+ABX*I_NAI_H2x3y_S_C5_b;
  Double I_NAI_H2x2yz_Px_C5_b = I_NAI_I3x2yz_S_C5_b+ABX*I_NAI_H2x2yz_S_C5_b;
  Double I_NAI_H2xy2z_Px_C5_b = I_NAI_I3xy2z_S_C5_b+ABX*I_NAI_H2xy2z_S_C5_b;
  Double I_NAI_H2x3z_Px_C5_b = I_NAI_I3x3z_S_C5_b+ABX*I_NAI_H2x3z_S_C5_b;
  Double I_NAI_Hx4y_Px_C5_b = I_NAI_I2x4y_S_C5_b+ABX*I_NAI_Hx4y_S_C5_b;
  Double I_NAI_Hx3yz_Px_C5_b = I_NAI_I2x3yz_S_C5_b+ABX*I_NAI_Hx3yz_S_C5_b;
  Double I_NAI_Hx2y2z_Px_C5_b = I_NAI_I2x2y2z_S_C5_b+ABX*I_NAI_Hx2y2z_S_C5_b;
  Double I_NAI_Hxy3z_Px_C5_b = I_NAI_I2xy3z_S_C5_b+ABX*I_NAI_Hxy3z_S_C5_b;
  Double I_NAI_Hx4z_Px_C5_b = I_NAI_I2x4z_S_C5_b+ABX*I_NAI_Hx4z_S_C5_b;
  Double I_NAI_H5y_Px_C5_b = I_NAI_Ix5y_S_C5_b+ABX*I_NAI_H5y_S_C5_b;
  Double I_NAI_H4yz_Px_C5_b = I_NAI_Ix4yz_S_C5_b+ABX*I_NAI_H4yz_S_C5_b;
  Double I_NAI_H3y2z_Px_C5_b = I_NAI_Ix3y2z_S_C5_b+ABX*I_NAI_H3y2z_S_C5_b;
  Double I_NAI_H2y3z_Px_C5_b = I_NAI_Ix2y3z_S_C5_b+ABX*I_NAI_H2y3z_S_C5_b;
  Double I_NAI_Hy4z_Px_C5_b = I_NAI_Ixy4z_S_C5_b+ABX*I_NAI_Hy4z_S_C5_b;
  Double I_NAI_H5z_Px_C5_b = I_NAI_Ix5z_S_C5_b+ABX*I_NAI_H5z_S_C5_b;
  Double I_NAI_H5x_Py_C5_b = I_NAI_I5xy_S_C5_b+ABY*I_NAI_H5x_S_C5_b;
  Double I_NAI_H4xy_Py_C5_b = I_NAI_I4x2y_S_C5_b+ABY*I_NAI_H4xy_S_C5_b;
  Double I_NAI_H4xz_Py_C5_b = I_NAI_I4xyz_S_C5_b+ABY*I_NAI_H4xz_S_C5_b;
  Double I_NAI_H3x2y_Py_C5_b = I_NAI_I3x3y_S_C5_b+ABY*I_NAI_H3x2y_S_C5_b;
  Double I_NAI_H3xyz_Py_C5_b = I_NAI_I3x2yz_S_C5_b+ABY*I_NAI_H3xyz_S_C5_b;
  Double I_NAI_H3x2z_Py_C5_b = I_NAI_I3xy2z_S_C5_b+ABY*I_NAI_H3x2z_S_C5_b;
  Double I_NAI_H2x3y_Py_C5_b = I_NAI_I2x4y_S_C5_b+ABY*I_NAI_H2x3y_S_C5_b;
  Double I_NAI_H2x2yz_Py_C5_b = I_NAI_I2x3yz_S_C5_b+ABY*I_NAI_H2x2yz_S_C5_b;
  Double I_NAI_H2xy2z_Py_C5_b = I_NAI_I2x2y2z_S_C5_b+ABY*I_NAI_H2xy2z_S_C5_b;
  Double I_NAI_H2x3z_Py_C5_b = I_NAI_I2xy3z_S_C5_b+ABY*I_NAI_H2x3z_S_C5_b;
  Double I_NAI_Hx4y_Py_C5_b = I_NAI_Ix5y_S_C5_b+ABY*I_NAI_Hx4y_S_C5_b;
  Double I_NAI_Hx3yz_Py_C5_b = I_NAI_Ix4yz_S_C5_b+ABY*I_NAI_Hx3yz_S_C5_b;
  Double I_NAI_Hx2y2z_Py_C5_b = I_NAI_Ix3y2z_S_C5_b+ABY*I_NAI_Hx2y2z_S_C5_b;
  Double I_NAI_Hxy3z_Py_C5_b = I_NAI_Ix2y3z_S_C5_b+ABY*I_NAI_Hxy3z_S_C5_b;
  Double I_NAI_Hx4z_Py_C5_b = I_NAI_Ixy4z_S_C5_b+ABY*I_NAI_Hx4z_S_C5_b;
  Double I_NAI_H5y_Py_C5_b = I_NAI_I6y_S_C5_b+ABY*I_NAI_H5y_S_C5_b;
  Double I_NAI_H4yz_Py_C5_b = I_NAI_I5yz_S_C5_b+ABY*I_NAI_H4yz_S_C5_b;
  Double I_NAI_H3y2z_Py_C5_b = I_NAI_I4y2z_S_C5_b+ABY*I_NAI_H3y2z_S_C5_b;
  Double I_NAI_H2y3z_Py_C5_b = I_NAI_I3y3z_S_C5_b+ABY*I_NAI_H2y3z_S_C5_b;
  Double I_NAI_Hy4z_Py_C5_b = I_NAI_I2y4z_S_C5_b+ABY*I_NAI_Hy4z_S_C5_b;
  Double I_NAI_H5z_Py_C5_b = I_NAI_Iy5z_S_C5_b+ABY*I_NAI_H5z_S_C5_b;
  Double I_NAI_H5x_Pz_C5_b = I_NAI_I5xz_S_C5_b+ABZ*I_NAI_H5x_S_C5_b;
  Double I_NAI_H4xy_Pz_C5_b = I_NAI_I4xyz_S_C5_b+ABZ*I_NAI_H4xy_S_C5_b;
  Double I_NAI_H4xz_Pz_C5_b = I_NAI_I4x2z_S_C5_b+ABZ*I_NAI_H4xz_S_C5_b;
  Double I_NAI_H3x2y_Pz_C5_b = I_NAI_I3x2yz_S_C5_b+ABZ*I_NAI_H3x2y_S_C5_b;
  Double I_NAI_H3xyz_Pz_C5_b = I_NAI_I3xy2z_S_C5_b+ABZ*I_NAI_H3xyz_S_C5_b;
  Double I_NAI_H3x2z_Pz_C5_b = I_NAI_I3x3z_S_C5_b+ABZ*I_NAI_H3x2z_S_C5_b;
  Double I_NAI_H2x3y_Pz_C5_b = I_NAI_I2x3yz_S_C5_b+ABZ*I_NAI_H2x3y_S_C5_b;
  Double I_NAI_H2x2yz_Pz_C5_b = I_NAI_I2x2y2z_S_C5_b+ABZ*I_NAI_H2x2yz_S_C5_b;
  Double I_NAI_H2xy2z_Pz_C5_b = I_NAI_I2xy3z_S_C5_b+ABZ*I_NAI_H2xy2z_S_C5_b;
  Double I_NAI_H2x3z_Pz_C5_b = I_NAI_I2x4z_S_C5_b+ABZ*I_NAI_H2x3z_S_C5_b;
  Double I_NAI_Hx4y_Pz_C5_b = I_NAI_Ix4yz_S_C5_b+ABZ*I_NAI_Hx4y_S_C5_b;
  Double I_NAI_Hx3yz_Pz_C5_b = I_NAI_Ix3y2z_S_C5_b+ABZ*I_NAI_Hx3yz_S_C5_b;
  Double I_NAI_Hx2y2z_Pz_C5_b = I_NAI_Ix2y3z_S_C5_b+ABZ*I_NAI_Hx2y2z_S_C5_b;
  Double I_NAI_Hxy3z_Pz_C5_b = I_NAI_Ixy4z_S_C5_b+ABZ*I_NAI_Hxy3z_S_C5_b;
  Double I_NAI_Hx4z_Pz_C5_b = I_NAI_Ix5z_S_C5_b+ABZ*I_NAI_Hx4z_S_C5_b;
  Double I_NAI_H5y_Pz_C5_b = I_NAI_I5yz_S_C5_b+ABZ*I_NAI_H5y_S_C5_b;
  Double I_NAI_H4yz_Pz_C5_b = I_NAI_I4y2z_S_C5_b+ABZ*I_NAI_H4yz_S_C5_b;
  Double I_NAI_H3y2z_Pz_C5_b = I_NAI_I3y3z_S_C5_b+ABZ*I_NAI_H3y2z_S_C5_b;
  Double I_NAI_H2y3z_Pz_C5_b = I_NAI_I2y4z_S_C5_b+ABZ*I_NAI_H2y3z_S_C5_b;
  Double I_NAI_Hy4z_Pz_C5_b = I_NAI_Iy5z_S_C5_b+ABZ*I_NAI_Hy4z_S_C5_b;
  Double I_NAI_H5z_Pz_C5_b = I_NAI_I6z_S_C5_b+ABZ*I_NAI_H5z_S_C5_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C1005_b
   * RHS shell quartet name: SQ_NAI_H_S_C1005_b
   ************************************************************/
  Double I_NAI_H5x_Px_C1005_b = I_NAI_I6x_S_C1005_b+ABX*I_NAI_H5x_S_C1005_b;
  Double I_NAI_H4xy_Px_C1005_b = I_NAI_I5xy_S_C1005_b+ABX*I_NAI_H4xy_S_C1005_b;
  Double I_NAI_H4xz_Px_C1005_b = I_NAI_I5xz_S_C1005_b+ABX*I_NAI_H4xz_S_C1005_b;
  Double I_NAI_H3x2y_Px_C1005_b = I_NAI_I4x2y_S_C1005_b+ABX*I_NAI_H3x2y_S_C1005_b;
  Double I_NAI_H3xyz_Px_C1005_b = I_NAI_I4xyz_S_C1005_b+ABX*I_NAI_H3xyz_S_C1005_b;
  Double I_NAI_H3x2z_Px_C1005_b = I_NAI_I4x2z_S_C1005_b+ABX*I_NAI_H3x2z_S_C1005_b;
  Double I_NAI_H2x3y_Px_C1005_b = I_NAI_I3x3y_S_C1005_b+ABX*I_NAI_H2x3y_S_C1005_b;
  Double I_NAI_H2x2yz_Px_C1005_b = I_NAI_I3x2yz_S_C1005_b+ABX*I_NAI_H2x2yz_S_C1005_b;
  Double I_NAI_H2xy2z_Px_C1005_b = I_NAI_I3xy2z_S_C1005_b+ABX*I_NAI_H2xy2z_S_C1005_b;
  Double I_NAI_H2x3z_Px_C1005_b = I_NAI_I3x3z_S_C1005_b+ABX*I_NAI_H2x3z_S_C1005_b;
  Double I_NAI_Hx4y_Px_C1005_b = I_NAI_I2x4y_S_C1005_b+ABX*I_NAI_Hx4y_S_C1005_b;
  Double I_NAI_Hx3yz_Px_C1005_b = I_NAI_I2x3yz_S_C1005_b+ABX*I_NAI_Hx3yz_S_C1005_b;
  Double I_NAI_Hx2y2z_Px_C1005_b = I_NAI_I2x2y2z_S_C1005_b+ABX*I_NAI_Hx2y2z_S_C1005_b;
  Double I_NAI_Hxy3z_Px_C1005_b = I_NAI_I2xy3z_S_C1005_b+ABX*I_NAI_Hxy3z_S_C1005_b;
  Double I_NAI_Hx4z_Px_C1005_b = I_NAI_I2x4z_S_C1005_b+ABX*I_NAI_Hx4z_S_C1005_b;
  Double I_NAI_H5y_Px_C1005_b = I_NAI_Ix5y_S_C1005_b+ABX*I_NAI_H5y_S_C1005_b;
  Double I_NAI_H4yz_Px_C1005_b = I_NAI_Ix4yz_S_C1005_b+ABX*I_NAI_H4yz_S_C1005_b;
  Double I_NAI_H3y2z_Px_C1005_b = I_NAI_Ix3y2z_S_C1005_b+ABX*I_NAI_H3y2z_S_C1005_b;
  Double I_NAI_H2y3z_Px_C1005_b = I_NAI_Ix2y3z_S_C1005_b+ABX*I_NAI_H2y3z_S_C1005_b;
  Double I_NAI_Hy4z_Px_C1005_b = I_NAI_Ixy4z_S_C1005_b+ABX*I_NAI_Hy4z_S_C1005_b;
  Double I_NAI_H5z_Px_C1005_b = I_NAI_Ix5z_S_C1005_b+ABX*I_NAI_H5z_S_C1005_b;
  Double I_NAI_H5x_Py_C1005_b = I_NAI_I5xy_S_C1005_b+ABY*I_NAI_H5x_S_C1005_b;
  Double I_NAI_H4xy_Py_C1005_b = I_NAI_I4x2y_S_C1005_b+ABY*I_NAI_H4xy_S_C1005_b;
  Double I_NAI_H4xz_Py_C1005_b = I_NAI_I4xyz_S_C1005_b+ABY*I_NAI_H4xz_S_C1005_b;
  Double I_NAI_H3x2y_Py_C1005_b = I_NAI_I3x3y_S_C1005_b+ABY*I_NAI_H3x2y_S_C1005_b;
  Double I_NAI_H3xyz_Py_C1005_b = I_NAI_I3x2yz_S_C1005_b+ABY*I_NAI_H3xyz_S_C1005_b;
  Double I_NAI_H3x2z_Py_C1005_b = I_NAI_I3xy2z_S_C1005_b+ABY*I_NAI_H3x2z_S_C1005_b;
  Double I_NAI_H2x3y_Py_C1005_b = I_NAI_I2x4y_S_C1005_b+ABY*I_NAI_H2x3y_S_C1005_b;
  Double I_NAI_H2x2yz_Py_C1005_b = I_NAI_I2x3yz_S_C1005_b+ABY*I_NAI_H2x2yz_S_C1005_b;
  Double I_NAI_H2xy2z_Py_C1005_b = I_NAI_I2x2y2z_S_C1005_b+ABY*I_NAI_H2xy2z_S_C1005_b;
  Double I_NAI_H2x3z_Py_C1005_b = I_NAI_I2xy3z_S_C1005_b+ABY*I_NAI_H2x3z_S_C1005_b;
  Double I_NAI_Hx4y_Py_C1005_b = I_NAI_Ix5y_S_C1005_b+ABY*I_NAI_Hx4y_S_C1005_b;
  Double I_NAI_Hx3yz_Py_C1005_b = I_NAI_Ix4yz_S_C1005_b+ABY*I_NAI_Hx3yz_S_C1005_b;
  Double I_NAI_Hx2y2z_Py_C1005_b = I_NAI_Ix3y2z_S_C1005_b+ABY*I_NAI_Hx2y2z_S_C1005_b;
  Double I_NAI_Hxy3z_Py_C1005_b = I_NAI_Ix2y3z_S_C1005_b+ABY*I_NAI_Hxy3z_S_C1005_b;
  Double I_NAI_Hx4z_Py_C1005_b = I_NAI_Ixy4z_S_C1005_b+ABY*I_NAI_Hx4z_S_C1005_b;
  Double I_NAI_H5y_Py_C1005_b = I_NAI_I6y_S_C1005_b+ABY*I_NAI_H5y_S_C1005_b;
  Double I_NAI_H4yz_Py_C1005_b = I_NAI_I5yz_S_C1005_b+ABY*I_NAI_H4yz_S_C1005_b;
  Double I_NAI_H3y2z_Py_C1005_b = I_NAI_I4y2z_S_C1005_b+ABY*I_NAI_H3y2z_S_C1005_b;
  Double I_NAI_H2y3z_Py_C1005_b = I_NAI_I3y3z_S_C1005_b+ABY*I_NAI_H2y3z_S_C1005_b;
  Double I_NAI_Hy4z_Py_C1005_b = I_NAI_I2y4z_S_C1005_b+ABY*I_NAI_Hy4z_S_C1005_b;
  Double I_NAI_H5z_Py_C1005_b = I_NAI_Iy5z_S_C1005_b+ABY*I_NAI_H5z_S_C1005_b;
  Double I_NAI_H5x_Pz_C1005_b = I_NAI_I5xz_S_C1005_b+ABZ*I_NAI_H5x_S_C1005_b;
  Double I_NAI_H4xy_Pz_C1005_b = I_NAI_I4xyz_S_C1005_b+ABZ*I_NAI_H4xy_S_C1005_b;
  Double I_NAI_H4xz_Pz_C1005_b = I_NAI_I4x2z_S_C1005_b+ABZ*I_NAI_H4xz_S_C1005_b;
  Double I_NAI_H3x2y_Pz_C1005_b = I_NAI_I3x2yz_S_C1005_b+ABZ*I_NAI_H3x2y_S_C1005_b;
  Double I_NAI_H3xyz_Pz_C1005_b = I_NAI_I3xy2z_S_C1005_b+ABZ*I_NAI_H3xyz_S_C1005_b;
  Double I_NAI_H3x2z_Pz_C1005_b = I_NAI_I3x3z_S_C1005_b+ABZ*I_NAI_H3x2z_S_C1005_b;
  Double I_NAI_H2x3y_Pz_C1005_b = I_NAI_I2x3yz_S_C1005_b+ABZ*I_NAI_H2x3y_S_C1005_b;
  Double I_NAI_H2x2yz_Pz_C1005_b = I_NAI_I2x2y2z_S_C1005_b+ABZ*I_NAI_H2x2yz_S_C1005_b;
  Double I_NAI_H2xy2z_Pz_C1005_b = I_NAI_I2xy3z_S_C1005_b+ABZ*I_NAI_H2xy2z_S_C1005_b;
  Double I_NAI_H2x3z_Pz_C1005_b = I_NAI_I2x4z_S_C1005_b+ABZ*I_NAI_H2x3z_S_C1005_b;
  Double I_NAI_Hx4y_Pz_C1005_b = I_NAI_Ix4yz_S_C1005_b+ABZ*I_NAI_Hx4y_S_C1005_b;
  Double I_NAI_Hx3yz_Pz_C1005_b = I_NAI_Ix3y2z_S_C1005_b+ABZ*I_NAI_Hx3yz_S_C1005_b;
  Double I_NAI_Hx2y2z_Pz_C1005_b = I_NAI_Ix2y3z_S_C1005_b+ABZ*I_NAI_Hx2y2z_S_C1005_b;
  Double I_NAI_Hxy3z_Pz_C1005_b = I_NAI_Ixy4z_S_C1005_b+ABZ*I_NAI_Hxy3z_S_C1005_b;
  Double I_NAI_Hx4z_Pz_C1005_b = I_NAI_Ix5z_S_C1005_b+ABZ*I_NAI_Hx4z_S_C1005_b;
  Double I_NAI_H5y_Pz_C1005_b = I_NAI_I5yz_S_C1005_b+ABZ*I_NAI_H5y_S_C1005_b;
  Double I_NAI_H4yz_Pz_C1005_b = I_NAI_I4y2z_S_C1005_b+ABZ*I_NAI_H4yz_S_C1005_b;
  Double I_NAI_H3y2z_Pz_C1005_b = I_NAI_I3y3z_S_C1005_b+ABZ*I_NAI_H3y2z_S_C1005_b;
  Double I_NAI_H2y3z_Pz_C1005_b = I_NAI_I2y4z_S_C1005_b+ABZ*I_NAI_H2y3z_S_C1005_b;
  Double I_NAI_Hy4z_Pz_C1005_b = I_NAI_Iy5z_S_C1005_b+ABZ*I_NAI_Hy4z_S_C1005_b;
  Double I_NAI_H5z_Pz_C1005_b = I_NAI_I6z_S_C1005_b+ABZ*I_NAI_H5z_S_C1005_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_C1005_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 8 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_C1005_b
   * RHS shell quartet name: SQ_NAI_I_S_C1005_b
   ************************************************************/
  Double I_NAI_I6x_Px_C1005_b = I_NAI_K7x_S_C1005_b+ABX*I_NAI_I6x_S_C1005_b;
  Double I_NAI_I5xy_Px_C1005_b = I_NAI_K6xy_S_C1005_b+ABX*I_NAI_I5xy_S_C1005_b;
  Double I_NAI_I5xz_Px_C1005_b = I_NAI_K6xz_S_C1005_b+ABX*I_NAI_I5xz_S_C1005_b;
  Double I_NAI_I4x2y_Px_C1005_b = I_NAI_K5x2y_S_C1005_b+ABX*I_NAI_I4x2y_S_C1005_b;
  Double I_NAI_I4xyz_Px_C1005_b = I_NAI_K5xyz_S_C1005_b+ABX*I_NAI_I4xyz_S_C1005_b;
  Double I_NAI_I4x2z_Px_C1005_b = I_NAI_K5x2z_S_C1005_b+ABX*I_NAI_I4x2z_S_C1005_b;
  Double I_NAI_I3x3y_Px_C1005_b = I_NAI_K4x3y_S_C1005_b+ABX*I_NAI_I3x3y_S_C1005_b;
  Double I_NAI_I3x2yz_Px_C1005_b = I_NAI_K4x2yz_S_C1005_b+ABX*I_NAI_I3x2yz_S_C1005_b;
  Double I_NAI_I3xy2z_Px_C1005_b = I_NAI_K4xy2z_S_C1005_b+ABX*I_NAI_I3xy2z_S_C1005_b;
  Double I_NAI_I3x3z_Px_C1005_b = I_NAI_K4x3z_S_C1005_b+ABX*I_NAI_I3x3z_S_C1005_b;
  Double I_NAI_I2x4y_Px_C1005_b = I_NAI_K3x4y_S_C1005_b+ABX*I_NAI_I2x4y_S_C1005_b;
  Double I_NAI_I2x3yz_Px_C1005_b = I_NAI_K3x3yz_S_C1005_b+ABX*I_NAI_I2x3yz_S_C1005_b;
  Double I_NAI_I2x2y2z_Px_C1005_b = I_NAI_K3x2y2z_S_C1005_b+ABX*I_NAI_I2x2y2z_S_C1005_b;
  Double I_NAI_I2xy3z_Px_C1005_b = I_NAI_K3xy3z_S_C1005_b+ABX*I_NAI_I2xy3z_S_C1005_b;
  Double I_NAI_I2x4z_Px_C1005_b = I_NAI_K3x4z_S_C1005_b+ABX*I_NAI_I2x4z_S_C1005_b;
  Double I_NAI_Ix5y_Px_C1005_b = I_NAI_K2x5y_S_C1005_b+ABX*I_NAI_Ix5y_S_C1005_b;
  Double I_NAI_Ix4yz_Px_C1005_b = I_NAI_K2x4yz_S_C1005_b+ABX*I_NAI_Ix4yz_S_C1005_b;
  Double I_NAI_Ix3y2z_Px_C1005_b = I_NAI_K2x3y2z_S_C1005_b+ABX*I_NAI_Ix3y2z_S_C1005_b;
  Double I_NAI_Ix2y3z_Px_C1005_b = I_NAI_K2x2y3z_S_C1005_b+ABX*I_NAI_Ix2y3z_S_C1005_b;
  Double I_NAI_Ixy4z_Px_C1005_b = I_NAI_K2xy4z_S_C1005_b+ABX*I_NAI_Ixy4z_S_C1005_b;
  Double I_NAI_Ix5z_Px_C1005_b = I_NAI_K2x5z_S_C1005_b+ABX*I_NAI_Ix5z_S_C1005_b;
  Double I_NAI_I6y_Px_C1005_b = I_NAI_Kx6y_S_C1005_b+ABX*I_NAI_I6y_S_C1005_b;
  Double I_NAI_I5yz_Px_C1005_b = I_NAI_Kx5yz_S_C1005_b+ABX*I_NAI_I5yz_S_C1005_b;
  Double I_NAI_I4y2z_Px_C1005_b = I_NAI_Kx4y2z_S_C1005_b+ABX*I_NAI_I4y2z_S_C1005_b;
  Double I_NAI_I3y3z_Px_C1005_b = I_NAI_Kx3y3z_S_C1005_b+ABX*I_NAI_I3y3z_S_C1005_b;
  Double I_NAI_I2y4z_Px_C1005_b = I_NAI_Kx2y4z_S_C1005_b+ABX*I_NAI_I2y4z_S_C1005_b;
  Double I_NAI_Iy5z_Px_C1005_b = I_NAI_Kxy5z_S_C1005_b+ABX*I_NAI_Iy5z_S_C1005_b;
  Double I_NAI_I6z_Px_C1005_b = I_NAI_Kx6z_S_C1005_b+ABX*I_NAI_I6z_S_C1005_b;
  Double I_NAI_I5xy_Py_C1005_b = I_NAI_K5x2y_S_C1005_b+ABY*I_NAI_I5xy_S_C1005_b;
  Double I_NAI_I5xz_Py_C1005_b = I_NAI_K5xyz_S_C1005_b+ABY*I_NAI_I5xz_S_C1005_b;
  Double I_NAI_I4x2y_Py_C1005_b = I_NAI_K4x3y_S_C1005_b+ABY*I_NAI_I4x2y_S_C1005_b;
  Double I_NAI_I4xyz_Py_C1005_b = I_NAI_K4x2yz_S_C1005_b+ABY*I_NAI_I4xyz_S_C1005_b;
  Double I_NAI_I4x2z_Py_C1005_b = I_NAI_K4xy2z_S_C1005_b+ABY*I_NAI_I4x2z_S_C1005_b;
  Double I_NAI_I3x3y_Py_C1005_b = I_NAI_K3x4y_S_C1005_b+ABY*I_NAI_I3x3y_S_C1005_b;
  Double I_NAI_I3x2yz_Py_C1005_b = I_NAI_K3x3yz_S_C1005_b+ABY*I_NAI_I3x2yz_S_C1005_b;
  Double I_NAI_I3xy2z_Py_C1005_b = I_NAI_K3x2y2z_S_C1005_b+ABY*I_NAI_I3xy2z_S_C1005_b;
  Double I_NAI_I3x3z_Py_C1005_b = I_NAI_K3xy3z_S_C1005_b+ABY*I_NAI_I3x3z_S_C1005_b;
  Double I_NAI_I2x4y_Py_C1005_b = I_NAI_K2x5y_S_C1005_b+ABY*I_NAI_I2x4y_S_C1005_b;
  Double I_NAI_I2x3yz_Py_C1005_b = I_NAI_K2x4yz_S_C1005_b+ABY*I_NAI_I2x3yz_S_C1005_b;
  Double I_NAI_I2x2y2z_Py_C1005_b = I_NAI_K2x3y2z_S_C1005_b+ABY*I_NAI_I2x2y2z_S_C1005_b;
  Double I_NAI_I2xy3z_Py_C1005_b = I_NAI_K2x2y3z_S_C1005_b+ABY*I_NAI_I2xy3z_S_C1005_b;
  Double I_NAI_I2x4z_Py_C1005_b = I_NAI_K2xy4z_S_C1005_b+ABY*I_NAI_I2x4z_S_C1005_b;
  Double I_NAI_Ix5y_Py_C1005_b = I_NAI_Kx6y_S_C1005_b+ABY*I_NAI_Ix5y_S_C1005_b;
  Double I_NAI_Ix4yz_Py_C1005_b = I_NAI_Kx5yz_S_C1005_b+ABY*I_NAI_Ix4yz_S_C1005_b;
  Double I_NAI_Ix3y2z_Py_C1005_b = I_NAI_Kx4y2z_S_C1005_b+ABY*I_NAI_Ix3y2z_S_C1005_b;
  Double I_NAI_Ix2y3z_Py_C1005_b = I_NAI_Kx3y3z_S_C1005_b+ABY*I_NAI_Ix2y3z_S_C1005_b;
  Double I_NAI_Ixy4z_Py_C1005_b = I_NAI_Kx2y4z_S_C1005_b+ABY*I_NAI_Ixy4z_S_C1005_b;
  Double I_NAI_Ix5z_Py_C1005_b = I_NAI_Kxy5z_S_C1005_b+ABY*I_NAI_Ix5z_S_C1005_b;
  Double I_NAI_I6y_Py_C1005_b = I_NAI_K7y_S_C1005_b+ABY*I_NAI_I6y_S_C1005_b;
  Double I_NAI_I5yz_Py_C1005_b = I_NAI_K6yz_S_C1005_b+ABY*I_NAI_I5yz_S_C1005_b;
  Double I_NAI_I4y2z_Py_C1005_b = I_NAI_K5y2z_S_C1005_b+ABY*I_NAI_I4y2z_S_C1005_b;
  Double I_NAI_I3y3z_Py_C1005_b = I_NAI_K4y3z_S_C1005_b+ABY*I_NAI_I3y3z_S_C1005_b;
  Double I_NAI_I2y4z_Py_C1005_b = I_NAI_K3y4z_S_C1005_b+ABY*I_NAI_I2y4z_S_C1005_b;
  Double I_NAI_Iy5z_Py_C1005_b = I_NAI_K2y5z_S_C1005_b+ABY*I_NAI_Iy5z_S_C1005_b;
  Double I_NAI_I6z_Py_C1005_b = I_NAI_Ky6z_S_C1005_b+ABY*I_NAI_I6z_S_C1005_b;
  Double I_NAI_I5xz_Pz_C1005_b = I_NAI_K5x2z_S_C1005_b+ABZ*I_NAI_I5xz_S_C1005_b;
  Double I_NAI_I4xyz_Pz_C1005_b = I_NAI_K4xy2z_S_C1005_b+ABZ*I_NAI_I4xyz_S_C1005_b;
  Double I_NAI_I4x2z_Pz_C1005_b = I_NAI_K4x3z_S_C1005_b+ABZ*I_NAI_I4x2z_S_C1005_b;
  Double I_NAI_I3x2yz_Pz_C1005_b = I_NAI_K3x2y2z_S_C1005_b+ABZ*I_NAI_I3x2yz_S_C1005_b;
  Double I_NAI_I3xy2z_Pz_C1005_b = I_NAI_K3xy3z_S_C1005_b+ABZ*I_NAI_I3xy2z_S_C1005_b;
  Double I_NAI_I3x3z_Pz_C1005_b = I_NAI_K3x4z_S_C1005_b+ABZ*I_NAI_I3x3z_S_C1005_b;
  Double I_NAI_I2x3yz_Pz_C1005_b = I_NAI_K2x3y2z_S_C1005_b+ABZ*I_NAI_I2x3yz_S_C1005_b;
  Double I_NAI_I2x2y2z_Pz_C1005_b = I_NAI_K2x2y3z_S_C1005_b+ABZ*I_NAI_I2x2y2z_S_C1005_b;
  Double I_NAI_I2xy3z_Pz_C1005_b = I_NAI_K2xy4z_S_C1005_b+ABZ*I_NAI_I2xy3z_S_C1005_b;
  Double I_NAI_I2x4z_Pz_C1005_b = I_NAI_K2x5z_S_C1005_b+ABZ*I_NAI_I2x4z_S_C1005_b;
  Double I_NAI_Ix4yz_Pz_C1005_b = I_NAI_Kx4y2z_S_C1005_b+ABZ*I_NAI_Ix4yz_S_C1005_b;
  Double I_NAI_Ix3y2z_Pz_C1005_b = I_NAI_Kx3y3z_S_C1005_b+ABZ*I_NAI_Ix3y2z_S_C1005_b;
  Double I_NAI_Ix2y3z_Pz_C1005_b = I_NAI_Kx2y4z_S_C1005_b+ABZ*I_NAI_Ix2y3z_S_C1005_b;
  Double I_NAI_Ixy4z_Pz_C1005_b = I_NAI_Kxy5z_S_C1005_b+ABZ*I_NAI_Ixy4z_S_C1005_b;
  Double I_NAI_Ix5z_Pz_C1005_b = I_NAI_Kx6z_S_C1005_b+ABZ*I_NAI_Ix5z_S_C1005_b;
  Double I_NAI_I5yz_Pz_C1005_b = I_NAI_K5y2z_S_C1005_b+ABZ*I_NAI_I5yz_S_C1005_b;
  Double I_NAI_I4y2z_Pz_C1005_b = I_NAI_K4y3z_S_C1005_b+ABZ*I_NAI_I4y2z_S_C1005_b;
  Double I_NAI_I3y3z_Pz_C1005_b = I_NAI_K3y4z_S_C1005_b+ABZ*I_NAI_I3y3z_S_C1005_b;
  Double I_NAI_I2y4z_Pz_C1005_b = I_NAI_K2y5z_S_C1005_b+ABZ*I_NAI_I2y4z_S_C1005_b;
  Double I_NAI_Iy5z_Pz_C1005_b = I_NAI_Ky6z_S_C1005_b+ABZ*I_NAI_Iy5z_S_C1005_b;
  Double I_NAI_I6z_Pz_C1005_b = I_NAI_K7z_S_C1005_b+ABZ*I_NAI_I6z_S_C1005_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_C1005_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_C1005_b
   * RHS shell quartet name: SQ_NAI_H_P_C1005_b
   ************************************************************/
  Double I_NAI_H5x_D2x_C1005_b = I_NAI_I6x_Px_C1005_b+ABX*I_NAI_H5x_Px_C1005_b;
  Double I_NAI_H4xy_D2x_C1005_b = I_NAI_I5xy_Px_C1005_b+ABX*I_NAI_H4xy_Px_C1005_b;
  Double I_NAI_H4xz_D2x_C1005_b = I_NAI_I5xz_Px_C1005_b+ABX*I_NAI_H4xz_Px_C1005_b;
  Double I_NAI_H3x2y_D2x_C1005_b = I_NAI_I4x2y_Px_C1005_b+ABX*I_NAI_H3x2y_Px_C1005_b;
  Double I_NAI_H3xyz_D2x_C1005_b = I_NAI_I4xyz_Px_C1005_b+ABX*I_NAI_H3xyz_Px_C1005_b;
  Double I_NAI_H3x2z_D2x_C1005_b = I_NAI_I4x2z_Px_C1005_b+ABX*I_NAI_H3x2z_Px_C1005_b;
  Double I_NAI_H2x3y_D2x_C1005_b = I_NAI_I3x3y_Px_C1005_b+ABX*I_NAI_H2x3y_Px_C1005_b;
  Double I_NAI_H2x2yz_D2x_C1005_b = I_NAI_I3x2yz_Px_C1005_b+ABX*I_NAI_H2x2yz_Px_C1005_b;
  Double I_NAI_H2xy2z_D2x_C1005_b = I_NAI_I3xy2z_Px_C1005_b+ABX*I_NAI_H2xy2z_Px_C1005_b;
  Double I_NAI_H2x3z_D2x_C1005_b = I_NAI_I3x3z_Px_C1005_b+ABX*I_NAI_H2x3z_Px_C1005_b;
  Double I_NAI_Hx4y_D2x_C1005_b = I_NAI_I2x4y_Px_C1005_b+ABX*I_NAI_Hx4y_Px_C1005_b;
  Double I_NAI_Hx3yz_D2x_C1005_b = I_NAI_I2x3yz_Px_C1005_b+ABX*I_NAI_Hx3yz_Px_C1005_b;
  Double I_NAI_Hx2y2z_D2x_C1005_b = I_NAI_I2x2y2z_Px_C1005_b+ABX*I_NAI_Hx2y2z_Px_C1005_b;
  Double I_NAI_Hxy3z_D2x_C1005_b = I_NAI_I2xy3z_Px_C1005_b+ABX*I_NAI_Hxy3z_Px_C1005_b;
  Double I_NAI_Hx4z_D2x_C1005_b = I_NAI_I2x4z_Px_C1005_b+ABX*I_NAI_Hx4z_Px_C1005_b;
  Double I_NAI_H5y_D2x_C1005_b = I_NAI_Ix5y_Px_C1005_b+ABX*I_NAI_H5y_Px_C1005_b;
  Double I_NAI_H4yz_D2x_C1005_b = I_NAI_Ix4yz_Px_C1005_b+ABX*I_NAI_H4yz_Px_C1005_b;
  Double I_NAI_H3y2z_D2x_C1005_b = I_NAI_Ix3y2z_Px_C1005_b+ABX*I_NAI_H3y2z_Px_C1005_b;
  Double I_NAI_H2y3z_D2x_C1005_b = I_NAI_Ix2y3z_Px_C1005_b+ABX*I_NAI_H2y3z_Px_C1005_b;
  Double I_NAI_Hy4z_D2x_C1005_b = I_NAI_Ixy4z_Px_C1005_b+ABX*I_NAI_Hy4z_Px_C1005_b;
  Double I_NAI_H5z_D2x_C1005_b = I_NAI_Ix5z_Px_C1005_b+ABX*I_NAI_H5z_Px_C1005_b;
  Double I_NAI_H5x_Dxy_C1005_b = I_NAI_I5xy_Px_C1005_b+ABY*I_NAI_H5x_Px_C1005_b;
  Double I_NAI_H4xy_Dxy_C1005_b = I_NAI_I4x2y_Px_C1005_b+ABY*I_NAI_H4xy_Px_C1005_b;
  Double I_NAI_H4xz_Dxy_C1005_b = I_NAI_I4xyz_Px_C1005_b+ABY*I_NAI_H4xz_Px_C1005_b;
  Double I_NAI_H3x2y_Dxy_C1005_b = I_NAI_I3x3y_Px_C1005_b+ABY*I_NAI_H3x2y_Px_C1005_b;
  Double I_NAI_H3xyz_Dxy_C1005_b = I_NAI_I3x2yz_Px_C1005_b+ABY*I_NAI_H3xyz_Px_C1005_b;
  Double I_NAI_H3x2z_Dxy_C1005_b = I_NAI_I3xy2z_Px_C1005_b+ABY*I_NAI_H3x2z_Px_C1005_b;
  Double I_NAI_H2x3y_Dxy_C1005_b = I_NAI_I2x4y_Px_C1005_b+ABY*I_NAI_H2x3y_Px_C1005_b;
  Double I_NAI_H2x2yz_Dxy_C1005_b = I_NAI_I2x3yz_Px_C1005_b+ABY*I_NAI_H2x2yz_Px_C1005_b;
  Double I_NAI_H2xy2z_Dxy_C1005_b = I_NAI_I2x2y2z_Px_C1005_b+ABY*I_NAI_H2xy2z_Px_C1005_b;
  Double I_NAI_H2x3z_Dxy_C1005_b = I_NAI_I2xy3z_Px_C1005_b+ABY*I_NAI_H2x3z_Px_C1005_b;
  Double I_NAI_Hx4y_Dxy_C1005_b = I_NAI_Ix5y_Px_C1005_b+ABY*I_NAI_Hx4y_Px_C1005_b;
  Double I_NAI_Hx3yz_Dxy_C1005_b = I_NAI_Ix4yz_Px_C1005_b+ABY*I_NAI_Hx3yz_Px_C1005_b;
  Double I_NAI_Hx2y2z_Dxy_C1005_b = I_NAI_Ix3y2z_Px_C1005_b+ABY*I_NAI_Hx2y2z_Px_C1005_b;
  Double I_NAI_Hxy3z_Dxy_C1005_b = I_NAI_Ix2y3z_Px_C1005_b+ABY*I_NAI_Hxy3z_Px_C1005_b;
  Double I_NAI_Hx4z_Dxy_C1005_b = I_NAI_Ixy4z_Px_C1005_b+ABY*I_NAI_Hx4z_Px_C1005_b;
  Double I_NAI_H5y_Dxy_C1005_b = I_NAI_I6y_Px_C1005_b+ABY*I_NAI_H5y_Px_C1005_b;
  Double I_NAI_H4yz_Dxy_C1005_b = I_NAI_I5yz_Px_C1005_b+ABY*I_NAI_H4yz_Px_C1005_b;
  Double I_NAI_H3y2z_Dxy_C1005_b = I_NAI_I4y2z_Px_C1005_b+ABY*I_NAI_H3y2z_Px_C1005_b;
  Double I_NAI_H2y3z_Dxy_C1005_b = I_NAI_I3y3z_Px_C1005_b+ABY*I_NAI_H2y3z_Px_C1005_b;
  Double I_NAI_Hy4z_Dxy_C1005_b = I_NAI_I2y4z_Px_C1005_b+ABY*I_NAI_Hy4z_Px_C1005_b;
  Double I_NAI_H5z_Dxy_C1005_b = I_NAI_Iy5z_Px_C1005_b+ABY*I_NAI_H5z_Px_C1005_b;
  Double I_NAI_H5x_Dxz_C1005_b = I_NAI_I5xz_Px_C1005_b+ABZ*I_NAI_H5x_Px_C1005_b;
  Double I_NAI_H4xy_Dxz_C1005_b = I_NAI_I4xyz_Px_C1005_b+ABZ*I_NAI_H4xy_Px_C1005_b;
  Double I_NAI_H4xz_Dxz_C1005_b = I_NAI_I4x2z_Px_C1005_b+ABZ*I_NAI_H4xz_Px_C1005_b;
  Double I_NAI_H3x2y_Dxz_C1005_b = I_NAI_I3x2yz_Px_C1005_b+ABZ*I_NAI_H3x2y_Px_C1005_b;
  Double I_NAI_H3xyz_Dxz_C1005_b = I_NAI_I3xy2z_Px_C1005_b+ABZ*I_NAI_H3xyz_Px_C1005_b;
  Double I_NAI_H3x2z_Dxz_C1005_b = I_NAI_I3x3z_Px_C1005_b+ABZ*I_NAI_H3x2z_Px_C1005_b;
  Double I_NAI_H2x3y_Dxz_C1005_b = I_NAI_I2x3yz_Px_C1005_b+ABZ*I_NAI_H2x3y_Px_C1005_b;
  Double I_NAI_H2x2yz_Dxz_C1005_b = I_NAI_I2x2y2z_Px_C1005_b+ABZ*I_NAI_H2x2yz_Px_C1005_b;
  Double I_NAI_H2xy2z_Dxz_C1005_b = I_NAI_I2xy3z_Px_C1005_b+ABZ*I_NAI_H2xy2z_Px_C1005_b;
  Double I_NAI_H2x3z_Dxz_C1005_b = I_NAI_I2x4z_Px_C1005_b+ABZ*I_NAI_H2x3z_Px_C1005_b;
  Double I_NAI_Hx4y_Dxz_C1005_b = I_NAI_Ix4yz_Px_C1005_b+ABZ*I_NAI_Hx4y_Px_C1005_b;
  Double I_NAI_Hx3yz_Dxz_C1005_b = I_NAI_Ix3y2z_Px_C1005_b+ABZ*I_NAI_Hx3yz_Px_C1005_b;
  Double I_NAI_Hx2y2z_Dxz_C1005_b = I_NAI_Ix2y3z_Px_C1005_b+ABZ*I_NAI_Hx2y2z_Px_C1005_b;
  Double I_NAI_Hxy3z_Dxz_C1005_b = I_NAI_Ixy4z_Px_C1005_b+ABZ*I_NAI_Hxy3z_Px_C1005_b;
  Double I_NAI_Hx4z_Dxz_C1005_b = I_NAI_Ix5z_Px_C1005_b+ABZ*I_NAI_Hx4z_Px_C1005_b;
  Double I_NAI_H5y_Dxz_C1005_b = I_NAI_I5yz_Px_C1005_b+ABZ*I_NAI_H5y_Px_C1005_b;
  Double I_NAI_H4yz_Dxz_C1005_b = I_NAI_I4y2z_Px_C1005_b+ABZ*I_NAI_H4yz_Px_C1005_b;
  Double I_NAI_H3y2z_Dxz_C1005_b = I_NAI_I3y3z_Px_C1005_b+ABZ*I_NAI_H3y2z_Px_C1005_b;
  Double I_NAI_H2y3z_Dxz_C1005_b = I_NAI_I2y4z_Px_C1005_b+ABZ*I_NAI_H2y3z_Px_C1005_b;
  Double I_NAI_Hy4z_Dxz_C1005_b = I_NAI_Iy5z_Px_C1005_b+ABZ*I_NAI_Hy4z_Px_C1005_b;
  Double I_NAI_H5z_Dxz_C1005_b = I_NAI_I6z_Px_C1005_b+ABZ*I_NAI_H5z_Px_C1005_b;
  Double I_NAI_H5x_D2y_C1005_b = I_NAI_I5xy_Py_C1005_b+ABY*I_NAI_H5x_Py_C1005_b;
  Double I_NAI_H4xy_D2y_C1005_b = I_NAI_I4x2y_Py_C1005_b+ABY*I_NAI_H4xy_Py_C1005_b;
  Double I_NAI_H4xz_D2y_C1005_b = I_NAI_I4xyz_Py_C1005_b+ABY*I_NAI_H4xz_Py_C1005_b;
  Double I_NAI_H3x2y_D2y_C1005_b = I_NAI_I3x3y_Py_C1005_b+ABY*I_NAI_H3x2y_Py_C1005_b;
  Double I_NAI_H3xyz_D2y_C1005_b = I_NAI_I3x2yz_Py_C1005_b+ABY*I_NAI_H3xyz_Py_C1005_b;
  Double I_NAI_H3x2z_D2y_C1005_b = I_NAI_I3xy2z_Py_C1005_b+ABY*I_NAI_H3x2z_Py_C1005_b;
  Double I_NAI_H2x3y_D2y_C1005_b = I_NAI_I2x4y_Py_C1005_b+ABY*I_NAI_H2x3y_Py_C1005_b;
  Double I_NAI_H2x2yz_D2y_C1005_b = I_NAI_I2x3yz_Py_C1005_b+ABY*I_NAI_H2x2yz_Py_C1005_b;
  Double I_NAI_H2xy2z_D2y_C1005_b = I_NAI_I2x2y2z_Py_C1005_b+ABY*I_NAI_H2xy2z_Py_C1005_b;
  Double I_NAI_H2x3z_D2y_C1005_b = I_NAI_I2xy3z_Py_C1005_b+ABY*I_NAI_H2x3z_Py_C1005_b;
  Double I_NAI_Hx4y_D2y_C1005_b = I_NAI_Ix5y_Py_C1005_b+ABY*I_NAI_Hx4y_Py_C1005_b;
  Double I_NAI_Hx3yz_D2y_C1005_b = I_NAI_Ix4yz_Py_C1005_b+ABY*I_NAI_Hx3yz_Py_C1005_b;
  Double I_NAI_Hx2y2z_D2y_C1005_b = I_NAI_Ix3y2z_Py_C1005_b+ABY*I_NAI_Hx2y2z_Py_C1005_b;
  Double I_NAI_Hxy3z_D2y_C1005_b = I_NAI_Ix2y3z_Py_C1005_b+ABY*I_NAI_Hxy3z_Py_C1005_b;
  Double I_NAI_Hx4z_D2y_C1005_b = I_NAI_Ixy4z_Py_C1005_b+ABY*I_NAI_Hx4z_Py_C1005_b;
  Double I_NAI_H5y_D2y_C1005_b = I_NAI_I6y_Py_C1005_b+ABY*I_NAI_H5y_Py_C1005_b;
  Double I_NAI_H4yz_D2y_C1005_b = I_NAI_I5yz_Py_C1005_b+ABY*I_NAI_H4yz_Py_C1005_b;
  Double I_NAI_H3y2z_D2y_C1005_b = I_NAI_I4y2z_Py_C1005_b+ABY*I_NAI_H3y2z_Py_C1005_b;
  Double I_NAI_H2y3z_D2y_C1005_b = I_NAI_I3y3z_Py_C1005_b+ABY*I_NAI_H2y3z_Py_C1005_b;
  Double I_NAI_Hy4z_D2y_C1005_b = I_NAI_I2y4z_Py_C1005_b+ABY*I_NAI_Hy4z_Py_C1005_b;
  Double I_NAI_H5z_D2y_C1005_b = I_NAI_Iy5z_Py_C1005_b+ABY*I_NAI_H5z_Py_C1005_b;
  Double I_NAI_H5x_Dyz_C1005_b = I_NAI_I5xz_Py_C1005_b+ABZ*I_NAI_H5x_Py_C1005_b;
  Double I_NAI_H4xy_Dyz_C1005_b = I_NAI_I4xyz_Py_C1005_b+ABZ*I_NAI_H4xy_Py_C1005_b;
  Double I_NAI_H4xz_Dyz_C1005_b = I_NAI_I4x2z_Py_C1005_b+ABZ*I_NAI_H4xz_Py_C1005_b;
  Double I_NAI_H3x2y_Dyz_C1005_b = I_NAI_I3x2yz_Py_C1005_b+ABZ*I_NAI_H3x2y_Py_C1005_b;
  Double I_NAI_H3xyz_Dyz_C1005_b = I_NAI_I3xy2z_Py_C1005_b+ABZ*I_NAI_H3xyz_Py_C1005_b;
  Double I_NAI_H3x2z_Dyz_C1005_b = I_NAI_I3x3z_Py_C1005_b+ABZ*I_NAI_H3x2z_Py_C1005_b;
  Double I_NAI_H2x3y_Dyz_C1005_b = I_NAI_I2x3yz_Py_C1005_b+ABZ*I_NAI_H2x3y_Py_C1005_b;
  Double I_NAI_H2x2yz_Dyz_C1005_b = I_NAI_I2x2y2z_Py_C1005_b+ABZ*I_NAI_H2x2yz_Py_C1005_b;
  Double I_NAI_H2xy2z_Dyz_C1005_b = I_NAI_I2xy3z_Py_C1005_b+ABZ*I_NAI_H2xy2z_Py_C1005_b;
  Double I_NAI_H2x3z_Dyz_C1005_b = I_NAI_I2x4z_Py_C1005_b+ABZ*I_NAI_H2x3z_Py_C1005_b;
  Double I_NAI_Hx4y_Dyz_C1005_b = I_NAI_Ix4yz_Py_C1005_b+ABZ*I_NAI_Hx4y_Py_C1005_b;
  Double I_NAI_Hx3yz_Dyz_C1005_b = I_NAI_Ix3y2z_Py_C1005_b+ABZ*I_NAI_Hx3yz_Py_C1005_b;
  Double I_NAI_Hx2y2z_Dyz_C1005_b = I_NAI_Ix2y3z_Py_C1005_b+ABZ*I_NAI_Hx2y2z_Py_C1005_b;
  Double I_NAI_Hxy3z_Dyz_C1005_b = I_NAI_Ixy4z_Py_C1005_b+ABZ*I_NAI_Hxy3z_Py_C1005_b;
  Double I_NAI_Hx4z_Dyz_C1005_b = I_NAI_Ix5z_Py_C1005_b+ABZ*I_NAI_Hx4z_Py_C1005_b;
  Double I_NAI_H5y_Dyz_C1005_b = I_NAI_I5yz_Py_C1005_b+ABZ*I_NAI_H5y_Py_C1005_b;
  Double I_NAI_H4yz_Dyz_C1005_b = I_NAI_I4y2z_Py_C1005_b+ABZ*I_NAI_H4yz_Py_C1005_b;
  Double I_NAI_H3y2z_Dyz_C1005_b = I_NAI_I3y3z_Py_C1005_b+ABZ*I_NAI_H3y2z_Py_C1005_b;
  Double I_NAI_H2y3z_Dyz_C1005_b = I_NAI_I2y4z_Py_C1005_b+ABZ*I_NAI_H2y3z_Py_C1005_b;
  Double I_NAI_Hy4z_Dyz_C1005_b = I_NAI_Iy5z_Py_C1005_b+ABZ*I_NAI_Hy4z_Py_C1005_b;
  Double I_NAI_H5z_Dyz_C1005_b = I_NAI_I6z_Py_C1005_b+ABZ*I_NAI_H5z_Py_C1005_b;
  Double I_NAI_H5x_D2z_C1005_b = I_NAI_I5xz_Pz_C1005_b+ABZ*I_NAI_H5x_Pz_C1005_b;
  Double I_NAI_H4xy_D2z_C1005_b = I_NAI_I4xyz_Pz_C1005_b+ABZ*I_NAI_H4xy_Pz_C1005_b;
  Double I_NAI_H4xz_D2z_C1005_b = I_NAI_I4x2z_Pz_C1005_b+ABZ*I_NAI_H4xz_Pz_C1005_b;
  Double I_NAI_H3x2y_D2z_C1005_b = I_NAI_I3x2yz_Pz_C1005_b+ABZ*I_NAI_H3x2y_Pz_C1005_b;
  Double I_NAI_H3xyz_D2z_C1005_b = I_NAI_I3xy2z_Pz_C1005_b+ABZ*I_NAI_H3xyz_Pz_C1005_b;
  Double I_NAI_H3x2z_D2z_C1005_b = I_NAI_I3x3z_Pz_C1005_b+ABZ*I_NAI_H3x2z_Pz_C1005_b;
  Double I_NAI_H2x3y_D2z_C1005_b = I_NAI_I2x3yz_Pz_C1005_b+ABZ*I_NAI_H2x3y_Pz_C1005_b;
  Double I_NAI_H2x2yz_D2z_C1005_b = I_NAI_I2x2y2z_Pz_C1005_b+ABZ*I_NAI_H2x2yz_Pz_C1005_b;
  Double I_NAI_H2xy2z_D2z_C1005_b = I_NAI_I2xy3z_Pz_C1005_b+ABZ*I_NAI_H2xy2z_Pz_C1005_b;
  Double I_NAI_H2x3z_D2z_C1005_b = I_NAI_I2x4z_Pz_C1005_b+ABZ*I_NAI_H2x3z_Pz_C1005_b;
  Double I_NAI_Hx4y_D2z_C1005_b = I_NAI_Ix4yz_Pz_C1005_b+ABZ*I_NAI_Hx4y_Pz_C1005_b;
  Double I_NAI_Hx3yz_D2z_C1005_b = I_NAI_Ix3y2z_Pz_C1005_b+ABZ*I_NAI_Hx3yz_Pz_C1005_b;
  Double I_NAI_Hx2y2z_D2z_C1005_b = I_NAI_Ix2y3z_Pz_C1005_b+ABZ*I_NAI_Hx2y2z_Pz_C1005_b;
  Double I_NAI_Hxy3z_D2z_C1005_b = I_NAI_Ixy4z_Pz_C1005_b+ABZ*I_NAI_Hxy3z_Pz_C1005_b;
  Double I_NAI_Hx4z_D2z_C1005_b = I_NAI_Ix5z_Pz_C1005_b+ABZ*I_NAI_Hx4z_Pz_C1005_b;
  Double I_NAI_H5y_D2z_C1005_b = I_NAI_I5yz_Pz_C1005_b+ABZ*I_NAI_H5y_Pz_C1005_b;
  Double I_NAI_H4yz_D2z_C1005_b = I_NAI_I4y2z_Pz_C1005_b+ABZ*I_NAI_H4yz_Pz_C1005_b;
  Double I_NAI_H3y2z_D2z_C1005_b = I_NAI_I3y3z_Pz_C1005_b+ABZ*I_NAI_H3y2z_Pz_C1005_b;
  Double I_NAI_H2y3z_D2z_C1005_b = I_NAI_I2y4z_Pz_C1005_b+ABZ*I_NAI_H2y3z_Pz_C1005_b;
  Double I_NAI_Hy4z_D2z_C1005_b = I_NAI_Iy5z_Pz_C1005_b+ABZ*I_NAI_Hy4z_Pz_C1005_b;
  Double I_NAI_H5z_D2z_C1005_b = I_NAI_I6z_Pz_C1005_b+ABZ*I_NAI_H5z_Pz_C1005_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C5_a
   * RHS shell quartet name: SQ_NAI_G_S_C5
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_I6x_S_C5_a-5*I_NAI_G4x_S_C5;
  abcd[1] = 2.0E0*I_NAI_I5xy_S_C5_a-4*I_NAI_G3xy_S_C5;
  abcd[2] = 2.0E0*I_NAI_I5xz_S_C5_a-4*I_NAI_G3xz_S_C5;
  abcd[3] = 2.0E0*I_NAI_I4x2y_S_C5_a-3*I_NAI_G2x2y_S_C5;
  abcd[4] = 2.0E0*I_NAI_I4xyz_S_C5_a-3*I_NAI_G2xyz_S_C5;
  abcd[5] = 2.0E0*I_NAI_I4x2z_S_C5_a-3*I_NAI_G2x2z_S_C5;
  abcd[6] = 2.0E0*I_NAI_I3x3y_S_C5_a-2*I_NAI_Gx3y_S_C5;
  abcd[7] = 2.0E0*I_NAI_I3x2yz_S_C5_a-2*I_NAI_Gx2yz_S_C5;
  abcd[8] = 2.0E0*I_NAI_I3xy2z_S_C5_a-2*I_NAI_Gxy2z_S_C5;
  abcd[9] = 2.0E0*I_NAI_I3x3z_S_C5_a-2*I_NAI_Gx3z_S_C5;
  abcd[10] = 2.0E0*I_NAI_I2x4y_S_C5_a-1*I_NAI_G4y_S_C5;
  abcd[11] = 2.0E0*I_NAI_I2x3yz_S_C5_a-1*I_NAI_G3yz_S_C5;
  abcd[12] = 2.0E0*I_NAI_I2x2y2z_S_C5_a-1*I_NAI_G2y2z_S_C5;
  abcd[13] = 2.0E0*I_NAI_I2xy3z_S_C5_a-1*I_NAI_Gy3z_S_C5;
  abcd[14] = 2.0E0*I_NAI_I2x4z_S_C5_a-1*I_NAI_G4z_S_C5;
  abcd[15] = 2.0E0*I_NAI_Ix5y_S_C5_a;
  abcd[16] = 2.0E0*I_NAI_Ix4yz_S_C5_a;
  abcd[17] = 2.0E0*I_NAI_Ix3y2z_S_C5_a;
  abcd[18] = 2.0E0*I_NAI_Ix2y3z_S_C5_a;
  abcd[19] = 2.0E0*I_NAI_Ixy4z_S_C5_a;
  abcd[20] = 2.0E0*I_NAI_Ix5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_C1005_a
   * RHS shell quartet name: SQ_NAI_G_P_C1005
   ************************************************************/
  abcd[21] = 2.0E0*I_NAI_I6x_Px_C1005_a-5*I_NAI_G4x_Px_C1005;
  abcd[22] = 2.0E0*I_NAI_I5xy_Px_C1005_a-4*I_NAI_G3xy_Px_C1005;
  abcd[23] = 2.0E0*I_NAI_I5xz_Px_C1005_a-4*I_NAI_G3xz_Px_C1005;
  abcd[24] = 2.0E0*I_NAI_I4x2y_Px_C1005_a-3*I_NAI_G2x2y_Px_C1005;
  abcd[25] = 2.0E0*I_NAI_I4xyz_Px_C1005_a-3*I_NAI_G2xyz_Px_C1005;
  abcd[26] = 2.0E0*I_NAI_I4x2z_Px_C1005_a-3*I_NAI_G2x2z_Px_C1005;
  abcd[27] = 2.0E0*I_NAI_I3x3y_Px_C1005_a-2*I_NAI_Gx3y_Px_C1005;
  abcd[28] = 2.0E0*I_NAI_I3x2yz_Px_C1005_a-2*I_NAI_Gx2yz_Px_C1005;
  abcd[29] = 2.0E0*I_NAI_I3xy2z_Px_C1005_a-2*I_NAI_Gxy2z_Px_C1005;
  abcd[30] = 2.0E0*I_NAI_I3x3z_Px_C1005_a-2*I_NAI_Gx3z_Px_C1005;
  abcd[31] = 2.0E0*I_NAI_I2x4y_Px_C1005_a-1*I_NAI_G4y_Px_C1005;
  abcd[32] = 2.0E0*I_NAI_I2x3yz_Px_C1005_a-1*I_NAI_G3yz_Px_C1005;
  abcd[33] = 2.0E0*I_NAI_I2x2y2z_Px_C1005_a-1*I_NAI_G2y2z_Px_C1005;
  abcd[34] = 2.0E0*I_NAI_I2xy3z_Px_C1005_a-1*I_NAI_Gy3z_Px_C1005;
  abcd[35] = 2.0E0*I_NAI_I2x4z_Px_C1005_a-1*I_NAI_G4z_Px_C1005;
  abcd[36] = 2.0E0*I_NAI_Ix5y_Px_C1005_a;
  abcd[37] = 2.0E0*I_NAI_Ix4yz_Px_C1005_a;
  abcd[38] = 2.0E0*I_NAI_Ix3y2z_Px_C1005_a;
  abcd[39] = 2.0E0*I_NAI_Ix2y3z_Px_C1005_a;
  abcd[40] = 2.0E0*I_NAI_Ixy4z_Px_C1005_a;
  abcd[41] = 2.0E0*I_NAI_Ix5z_Px_C1005_a;
  abcd[42] = 2.0E0*I_NAI_I6x_Py_C1005_a-5*I_NAI_G4x_Py_C1005;
  abcd[43] = 2.0E0*I_NAI_I5xy_Py_C1005_a-4*I_NAI_G3xy_Py_C1005;
  abcd[44] = 2.0E0*I_NAI_I5xz_Py_C1005_a-4*I_NAI_G3xz_Py_C1005;
  abcd[45] = 2.0E0*I_NAI_I4x2y_Py_C1005_a-3*I_NAI_G2x2y_Py_C1005;
  abcd[46] = 2.0E0*I_NAI_I4xyz_Py_C1005_a-3*I_NAI_G2xyz_Py_C1005;
  abcd[47] = 2.0E0*I_NAI_I4x2z_Py_C1005_a-3*I_NAI_G2x2z_Py_C1005;
  abcd[48] = 2.0E0*I_NAI_I3x3y_Py_C1005_a-2*I_NAI_Gx3y_Py_C1005;
  abcd[49] = 2.0E0*I_NAI_I3x2yz_Py_C1005_a-2*I_NAI_Gx2yz_Py_C1005;
  abcd[50] = 2.0E0*I_NAI_I3xy2z_Py_C1005_a-2*I_NAI_Gxy2z_Py_C1005;
  abcd[51] = 2.0E0*I_NAI_I3x3z_Py_C1005_a-2*I_NAI_Gx3z_Py_C1005;
  abcd[52] = 2.0E0*I_NAI_I2x4y_Py_C1005_a-1*I_NAI_G4y_Py_C1005;
  abcd[53] = 2.0E0*I_NAI_I2x3yz_Py_C1005_a-1*I_NAI_G3yz_Py_C1005;
  abcd[54] = 2.0E0*I_NAI_I2x2y2z_Py_C1005_a-1*I_NAI_G2y2z_Py_C1005;
  abcd[55] = 2.0E0*I_NAI_I2xy3z_Py_C1005_a-1*I_NAI_Gy3z_Py_C1005;
  abcd[56] = 2.0E0*I_NAI_I2x4z_Py_C1005_a-1*I_NAI_G4z_Py_C1005;
  abcd[57] = 2.0E0*I_NAI_Ix5y_Py_C1005_a;
  abcd[58] = 2.0E0*I_NAI_Ix4yz_Py_C1005_a;
  abcd[59] = 2.0E0*I_NAI_Ix3y2z_Py_C1005_a;
  abcd[60] = 2.0E0*I_NAI_Ix2y3z_Py_C1005_a;
  abcd[61] = 2.0E0*I_NAI_Ixy4z_Py_C1005_a;
  abcd[62] = 2.0E0*I_NAI_Ix5z_Py_C1005_a;
  abcd[63] = 2.0E0*I_NAI_I6x_Pz_C1005_a-5*I_NAI_G4x_Pz_C1005;
  abcd[64] = 2.0E0*I_NAI_I5xy_Pz_C1005_a-4*I_NAI_G3xy_Pz_C1005;
  abcd[65] = 2.0E0*I_NAI_I5xz_Pz_C1005_a-4*I_NAI_G3xz_Pz_C1005;
  abcd[66] = 2.0E0*I_NAI_I4x2y_Pz_C1005_a-3*I_NAI_G2x2y_Pz_C1005;
  abcd[67] = 2.0E0*I_NAI_I4xyz_Pz_C1005_a-3*I_NAI_G2xyz_Pz_C1005;
  abcd[68] = 2.0E0*I_NAI_I4x2z_Pz_C1005_a-3*I_NAI_G2x2z_Pz_C1005;
  abcd[69] = 2.0E0*I_NAI_I3x3y_Pz_C1005_a-2*I_NAI_Gx3y_Pz_C1005;
  abcd[70] = 2.0E0*I_NAI_I3x2yz_Pz_C1005_a-2*I_NAI_Gx2yz_Pz_C1005;
  abcd[71] = 2.0E0*I_NAI_I3xy2z_Pz_C1005_a-2*I_NAI_Gxy2z_Pz_C1005;
  abcd[72] = 2.0E0*I_NAI_I3x3z_Pz_C1005_a-2*I_NAI_Gx3z_Pz_C1005;
  abcd[73] = 2.0E0*I_NAI_I2x4y_Pz_C1005_a-1*I_NAI_G4y_Pz_C1005;
  abcd[74] = 2.0E0*I_NAI_I2x3yz_Pz_C1005_a-1*I_NAI_G3yz_Pz_C1005;
  abcd[75] = 2.0E0*I_NAI_I2x2y2z_Pz_C1005_a-1*I_NAI_G2y2z_Pz_C1005;
  abcd[76] = 2.0E0*I_NAI_I2xy3z_Pz_C1005_a-1*I_NAI_Gy3z_Pz_C1005;
  abcd[77] = 2.0E0*I_NAI_I2x4z_Pz_C1005_a-1*I_NAI_G4z_Pz_C1005;
  abcd[78] = 2.0E0*I_NAI_Ix5y_Pz_C1005_a;
  abcd[79] = 2.0E0*I_NAI_Ix4yz_Pz_C1005_a;
  abcd[80] = 2.0E0*I_NAI_Ix3y2z_Pz_C1005_a;
  abcd[81] = 2.0E0*I_NAI_Ix2y3z_Pz_C1005_a;
  abcd[82] = 2.0E0*I_NAI_Ixy4z_Pz_C1005_a;
  abcd[83] = 2.0E0*I_NAI_Ix5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C5_a
   * RHS shell quartet name: SQ_NAI_G_S_C5
   ************************************************************/
  abcd[84] = 2.0E0*I_NAI_I5xy_S_C5_a;
  abcd[85] = 2.0E0*I_NAI_I4x2y_S_C5_a-1*I_NAI_G4x_S_C5;
  abcd[86] = 2.0E0*I_NAI_I4xyz_S_C5_a;
  abcd[87] = 2.0E0*I_NAI_I3x3y_S_C5_a-2*I_NAI_G3xy_S_C5;
  abcd[88] = 2.0E0*I_NAI_I3x2yz_S_C5_a-1*I_NAI_G3xz_S_C5;
  abcd[89] = 2.0E0*I_NAI_I3xy2z_S_C5_a;
  abcd[90] = 2.0E0*I_NAI_I2x4y_S_C5_a-3*I_NAI_G2x2y_S_C5;
  abcd[91] = 2.0E0*I_NAI_I2x3yz_S_C5_a-2*I_NAI_G2xyz_S_C5;
  abcd[92] = 2.0E0*I_NAI_I2x2y2z_S_C5_a-1*I_NAI_G2x2z_S_C5;
  abcd[93] = 2.0E0*I_NAI_I2xy3z_S_C5_a;
  abcd[94] = 2.0E0*I_NAI_Ix5y_S_C5_a-4*I_NAI_Gx3y_S_C5;
  abcd[95] = 2.0E0*I_NAI_Ix4yz_S_C5_a-3*I_NAI_Gx2yz_S_C5;
  abcd[96] = 2.0E0*I_NAI_Ix3y2z_S_C5_a-2*I_NAI_Gxy2z_S_C5;
  abcd[97] = 2.0E0*I_NAI_Ix2y3z_S_C5_a-1*I_NAI_Gx3z_S_C5;
  abcd[98] = 2.0E0*I_NAI_Ixy4z_S_C5_a;
  abcd[99] = 2.0E0*I_NAI_I6y_S_C5_a-5*I_NAI_G4y_S_C5;
  abcd[100] = 2.0E0*I_NAI_I5yz_S_C5_a-4*I_NAI_G3yz_S_C5;
  abcd[101] = 2.0E0*I_NAI_I4y2z_S_C5_a-3*I_NAI_G2y2z_S_C5;
  abcd[102] = 2.0E0*I_NAI_I3y3z_S_C5_a-2*I_NAI_Gy3z_S_C5;
  abcd[103] = 2.0E0*I_NAI_I2y4z_S_C5_a-1*I_NAI_G4z_S_C5;
  abcd[104] = 2.0E0*I_NAI_Iy5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_C1005_a
   * RHS shell quartet name: SQ_NAI_G_P_C1005
   ************************************************************/
  abcd[105] = 2.0E0*I_NAI_I5xy_Px_C1005_a;
  abcd[106] = 2.0E0*I_NAI_I4x2y_Px_C1005_a-1*I_NAI_G4x_Px_C1005;
  abcd[107] = 2.0E0*I_NAI_I4xyz_Px_C1005_a;
  abcd[108] = 2.0E0*I_NAI_I3x3y_Px_C1005_a-2*I_NAI_G3xy_Px_C1005;
  abcd[109] = 2.0E0*I_NAI_I3x2yz_Px_C1005_a-1*I_NAI_G3xz_Px_C1005;
  abcd[110] = 2.0E0*I_NAI_I3xy2z_Px_C1005_a;
  abcd[111] = 2.0E0*I_NAI_I2x4y_Px_C1005_a-3*I_NAI_G2x2y_Px_C1005;
  abcd[112] = 2.0E0*I_NAI_I2x3yz_Px_C1005_a-2*I_NAI_G2xyz_Px_C1005;
  abcd[113] = 2.0E0*I_NAI_I2x2y2z_Px_C1005_a-1*I_NAI_G2x2z_Px_C1005;
  abcd[114] = 2.0E0*I_NAI_I2xy3z_Px_C1005_a;
  abcd[115] = 2.0E0*I_NAI_Ix5y_Px_C1005_a-4*I_NAI_Gx3y_Px_C1005;
  abcd[116] = 2.0E0*I_NAI_Ix4yz_Px_C1005_a-3*I_NAI_Gx2yz_Px_C1005;
  abcd[117] = 2.0E0*I_NAI_Ix3y2z_Px_C1005_a-2*I_NAI_Gxy2z_Px_C1005;
  abcd[118] = 2.0E0*I_NAI_Ix2y3z_Px_C1005_a-1*I_NAI_Gx3z_Px_C1005;
  abcd[119] = 2.0E0*I_NAI_Ixy4z_Px_C1005_a;
  abcd[120] = 2.0E0*I_NAI_I6y_Px_C1005_a-5*I_NAI_G4y_Px_C1005;
  abcd[121] = 2.0E0*I_NAI_I5yz_Px_C1005_a-4*I_NAI_G3yz_Px_C1005;
  abcd[122] = 2.0E0*I_NAI_I4y2z_Px_C1005_a-3*I_NAI_G2y2z_Px_C1005;
  abcd[123] = 2.0E0*I_NAI_I3y3z_Px_C1005_a-2*I_NAI_Gy3z_Px_C1005;
  abcd[124] = 2.0E0*I_NAI_I2y4z_Px_C1005_a-1*I_NAI_G4z_Px_C1005;
  abcd[125] = 2.0E0*I_NAI_Iy5z_Px_C1005_a;
  abcd[126] = 2.0E0*I_NAI_I5xy_Py_C1005_a;
  abcd[127] = 2.0E0*I_NAI_I4x2y_Py_C1005_a-1*I_NAI_G4x_Py_C1005;
  abcd[128] = 2.0E0*I_NAI_I4xyz_Py_C1005_a;
  abcd[129] = 2.0E0*I_NAI_I3x3y_Py_C1005_a-2*I_NAI_G3xy_Py_C1005;
  abcd[130] = 2.0E0*I_NAI_I3x2yz_Py_C1005_a-1*I_NAI_G3xz_Py_C1005;
  abcd[131] = 2.0E0*I_NAI_I3xy2z_Py_C1005_a;
  abcd[132] = 2.0E0*I_NAI_I2x4y_Py_C1005_a-3*I_NAI_G2x2y_Py_C1005;
  abcd[133] = 2.0E0*I_NAI_I2x3yz_Py_C1005_a-2*I_NAI_G2xyz_Py_C1005;
  abcd[134] = 2.0E0*I_NAI_I2x2y2z_Py_C1005_a-1*I_NAI_G2x2z_Py_C1005;
  abcd[135] = 2.0E0*I_NAI_I2xy3z_Py_C1005_a;
  abcd[136] = 2.0E0*I_NAI_Ix5y_Py_C1005_a-4*I_NAI_Gx3y_Py_C1005;
  abcd[137] = 2.0E0*I_NAI_Ix4yz_Py_C1005_a-3*I_NAI_Gx2yz_Py_C1005;
  abcd[138] = 2.0E0*I_NAI_Ix3y2z_Py_C1005_a-2*I_NAI_Gxy2z_Py_C1005;
  abcd[139] = 2.0E0*I_NAI_Ix2y3z_Py_C1005_a-1*I_NAI_Gx3z_Py_C1005;
  abcd[140] = 2.0E0*I_NAI_Ixy4z_Py_C1005_a;
  abcd[141] = 2.0E0*I_NAI_I6y_Py_C1005_a-5*I_NAI_G4y_Py_C1005;
  abcd[142] = 2.0E0*I_NAI_I5yz_Py_C1005_a-4*I_NAI_G3yz_Py_C1005;
  abcd[143] = 2.0E0*I_NAI_I4y2z_Py_C1005_a-3*I_NAI_G2y2z_Py_C1005;
  abcd[144] = 2.0E0*I_NAI_I3y3z_Py_C1005_a-2*I_NAI_Gy3z_Py_C1005;
  abcd[145] = 2.0E0*I_NAI_I2y4z_Py_C1005_a-1*I_NAI_G4z_Py_C1005;
  abcd[146] = 2.0E0*I_NAI_Iy5z_Py_C1005_a;
  abcd[147] = 2.0E0*I_NAI_I5xy_Pz_C1005_a;
  abcd[148] = 2.0E0*I_NAI_I4x2y_Pz_C1005_a-1*I_NAI_G4x_Pz_C1005;
  abcd[149] = 2.0E0*I_NAI_I4xyz_Pz_C1005_a;
  abcd[150] = 2.0E0*I_NAI_I3x3y_Pz_C1005_a-2*I_NAI_G3xy_Pz_C1005;
  abcd[151] = 2.0E0*I_NAI_I3x2yz_Pz_C1005_a-1*I_NAI_G3xz_Pz_C1005;
  abcd[152] = 2.0E0*I_NAI_I3xy2z_Pz_C1005_a;
  abcd[153] = 2.0E0*I_NAI_I2x4y_Pz_C1005_a-3*I_NAI_G2x2y_Pz_C1005;
  abcd[154] = 2.0E0*I_NAI_I2x3yz_Pz_C1005_a-2*I_NAI_G2xyz_Pz_C1005;
  abcd[155] = 2.0E0*I_NAI_I2x2y2z_Pz_C1005_a-1*I_NAI_G2x2z_Pz_C1005;
  abcd[156] = 2.0E0*I_NAI_I2xy3z_Pz_C1005_a;
  abcd[157] = 2.0E0*I_NAI_Ix5y_Pz_C1005_a-4*I_NAI_Gx3y_Pz_C1005;
  abcd[158] = 2.0E0*I_NAI_Ix4yz_Pz_C1005_a-3*I_NAI_Gx2yz_Pz_C1005;
  abcd[159] = 2.0E0*I_NAI_Ix3y2z_Pz_C1005_a-2*I_NAI_Gxy2z_Pz_C1005;
  abcd[160] = 2.0E0*I_NAI_Ix2y3z_Pz_C1005_a-1*I_NAI_Gx3z_Pz_C1005;
  abcd[161] = 2.0E0*I_NAI_Ixy4z_Pz_C1005_a;
  abcd[162] = 2.0E0*I_NAI_I6y_Pz_C1005_a-5*I_NAI_G4y_Pz_C1005;
  abcd[163] = 2.0E0*I_NAI_I5yz_Pz_C1005_a-4*I_NAI_G3yz_Pz_C1005;
  abcd[164] = 2.0E0*I_NAI_I4y2z_Pz_C1005_a-3*I_NAI_G2y2z_Pz_C1005;
  abcd[165] = 2.0E0*I_NAI_I3y3z_Pz_C1005_a-2*I_NAI_Gy3z_Pz_C1005;
  abcd[166] = 2.0E0*I_NAI_I2y4z_Pz_C1005_a-1*I_NAI_G4z_Pz_C1005;
  abcd[167] = 2.0E0*I_NAI_Iy5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C5_a
   * RHS shell quartet name: SQ_NAI_G_S_C5
   ************************************************************/
  abcd[168] = 2.0E0*I_NAI_I5xz_S_C5_a;
  abcd[169] = 2.0E0*I_NAI_I4xyz_S_C5_a;
  abcd[170] = 2.0E0*I_NAI_I4x2z_S_C5_a-1*I_NAI_G4x_S_C5;
  abcd[171] = 2.0E0*I_NAI_I3x2yz_S_C5_a;
  abcd[172] = 2.0E0*I_NAI_I3xy2z_S_C5_a-1*I_NAI_G3xy_S_C5;
  abcd[173] = 2.0E0*I_NAI_I3x3z_S_C5_a-2*I_NAI_G3xz_S_C5;
  abcd[174] = 2.0E0*I_NAI_I2x3yz_S_C5_a;
  abcd[175] = 2.0E0*I_NAI_I2x2y2z_S_C5_a-1*I_NAI_G2x2y_S_C5;
  abcd[176] = 2.0E0*I_NAI_I2xy3z_S_C5_a-2*I_NAI_G2xyz_S_C5;
  abcd[177] = 2.0E0*I_NAI_I2x4z_S_C5_a-3*I_NAI_G2x2z_S_C5;
  abcd[178] = 2.0E0*I_NAI_Ix4yz_S_C5_a;
  abcd[179] = 2.0E0*I_NAI_Ix3y2z_S_C5_a-1*I_NAI_Gx3y_S_C5;
  abcd[180] = 2.0E0*I_NAI_Ix2y3z_S_C5_a-2*I_NAI_Gx2yz_S_C5;
  abcd[181] = 2.0E0*I_NAI_Ixy4z_S_C5_a-3*I_NAI_Gxy2z_S_C5;
  abcd[182] = 2.0E0*I_NAI_Ix5z_S_C5_a-4*I_NAI_Gx3z_S_C5;
  abcd[183] = 2.0E0*I_NAI_I5yz_S_C5_a;
  abcd[184] = 2.0E0*I_NAI_I4y2z_S_C5_a-1*I_NAI_G4y_S_C5;
  abcd[185] = 2.0E0*I_NAI_I3y3z_S_C5_a-2*I_NAI_G3yz_S_C5;
  abcd[186] = 2.0E0*I_NAI_I2y4z_S_C5_a-3*I_NAI_G2y2z_S_C5;
  abcd[187] = 2.0E0*I_NAI_Iy5z_S_C5_a-4*I_NAI_Gy3z_S_C5;
  abcd[188] = 2.0E0*I_NAI_I6z_S_C5_a-5*I_NAI_G4z_S_C5;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_C1005_a
   * RHS shell quartet name: SQ_NAI_G_P_C1005
   ************************************************************/
  abcd[189] = 2.0E0*I_NAI_I5xz_Px_C1005_a;
  abcd[190] = 2.0E0*I_NAI_I4xyz_Px_C1005_a;
  abcd[191] = 2.0E0*I_NAI_I4x2z_Px_C1005_a-1*I_NAI_G4x_Px_C1005;
  abcd[192] = 2.0E0*I_NAI_I3x2yz_Px_C1005_a;
  abcd[193] = 2.0E0*I_NAI_I3xy2z_Px_C1005_a-1*I_NAI_G3xy_Px_C1005;
  abcd[194] = 2.0E0*I_NAI_I3x3z_Px_C1005_a-2*I_NAI_G3xz_Px_C1005;
  abcd[195] = 2.0E0*I_NAI_I2x3yz_Px_C1005_a;
  abcd[196] = 2.0E0*I_NAI_I2x2y2z_Px_C1005_a-1*I_NAI_G2x2y_Px_C1005;
  abcd[197] = 2.0E0*I_NAI_I2xy3z_Px_C1005_a-2*I_NAI_G2xyz_Px_C1005;
  abcd[198] = 2.0E0*I_NAI_I2x4z_Px_C1005_a-3*I_NAI_G2x2z_Px_C1005;
  abcd[199] = 2.0E0*I_NAI_Ix4yz_Px_C1005_a;
  abcd[200] = 2.0E0*I_NAI_Ix3y2z_Px_C1005_a-1*I_NAI_Gx3y_Px_C1005;
  abcd[201] = 2.0E0*I_NAI_Ix2y3z_Px_C1005_a-2*I_NAI_Gx2yz_Px_C1005;
  abcd[202] = 2.0E0*I_NAI_Ixy4z_Px_C1005_a-3*I_NAI_Gxy2z_Px_C1005;
  abcd[203] = 2.0E0*I_NAI_Ix5z_Px_C1005_a-4*I_NAI_Gx3z_Px_C1005;
  abcd[204] = 2.0E0*I_NAI_I5yz_Px_C1005_a;
  abcd[205] = 2.0E0*I_NAI_I4y2z_Px_C1005_a-1*I_NAI_G4y_Px_C1005;
  abcd[206] = 2.0E0*I_NAI_I3y3z_Px_C1005_a-2*I_NAI_G3yz_Px_C1005;
  abcd[207] = 2.0E0*I_NAI_I2y4z_Px_C1005_a-3*I_NAI_G2y2z_Px_C1005;
  abcd[208] = 2.0E0*I_NAI_Iy5z_Px_C1005_a-4*I_NAI_Gy3z_Px_C1005;
  abcd[209] = 2.0E0*I_NAI_I6z_Px_C1005_a-5*I_NAI_G4z_Px_C1005;
  abcd[210] = 2.0E0*I_NAI_I5xz_Py_C1005_a;
  abcd[211] = 2.0E0*I_NAI_I4xyz_Py_C1005_a;
  abcd[212] = 2.0E0*I_NAI_I4x2z_Py_C1005_a-1*I_NAI_G4x_Py_C1005;
  abcd[213] = 2.0E0*I_NAI_I3x2yz_Py_C1005_a;
  abcd[214] = 2.0E0*I_NAI_I3xy2z_Py_C1005_a-1*I_NAI_G3xy_Py_C1005;
  abcd[215] = 2.0E0*I_NAI_I3x3z_Py_C1005_a-2*I_NAI_G3xz_Py_C1005;
  abcd[216] = 2.0E0*I_NAI_I2x3yz_Py_C1005_a;
  abcd[217] = 2.0E0*I_NAI_I2x2y2z_Py_C1005_a-1*I_NAI_G2x2y_Py_C1005;
  abcd[218] = 2.0E0*I_NAI_I2xy3z_Py_C1005_a-2*I_NAI_G2xyz_Py_C1005;
  abcd[219] = 2.0E0*I_NAI_I2x4z_Py_C1005_a-3*I_NAI_G2x2z_Py_C1005;
  abcd[220] = 2.0E0*I_NAI_Ix4yz_Py_C1005_a;
  abcd[221] = 2.0E0*I_NAI_Ix3y2z_Py_C1005_a-1*I_NAI_Gx3y_Py_C1005;
  abcd[222] = 2.0E0*I_NAI_Ix2y3z_Py_C1005_a-2*I_NAI_Gx2yz_Py_C1005;
  abcd[223] = 2.0E0*I_NAI_Ixy4z_Py_C1005_a-3*I_NAI_Gxy2z_Py_C1005;
  abcd[224] = 2.0E0*I_NAI_Ix5z_Py_C1005_a-4*I_NAI_Gx3z_Py_C1005;
  abcd[225] = 2.0E0*I_NAI_I5yz_Py_C1005_a;
  abcd[226] = 2.0E0*I_NAI_I4y2z_Py_C1005_a-1*I_NAI_G4y_Py_C1005;
  abcd[227] = 2.0E0*I_NAI_I3y3z_Py_C1005_a-2*I_NAI_G3yz_Py_C1005;
  abcd[228] = 2.0E0*I_NAI_I2y4z_Py_C1005_a-3*I_NAI_G2y2z_Py_C1005;
  abcd[229] = 2.0E0*I_NAI_Iy5z_Py_C1005_a-4*I_NAI_Gy3z_Py_C1005;
  abcd[230] = 2.0E0*I_NAI_I6z_Py_C1005_a-5*I_NAI_G4z_Py_C1005;
  abcd[231] = 2.0E0*I_NAI_I5xz_Pz_C1005_a;
  abcd[232] = 2.0E0*I_NAI_I4xyz_Pz_C1005_a;
  abcd[233] = 2.0E0*I_NAI_I4x2z_Pz_C1005_a-1*I_NAI_G4x_Pz_C1005;
  abcd[234] = 2.0E0*I_NAI_I3x2yz_Pz_C1005_a;
  abcd[235] = 2.0E0*I_NAI_I3xy2z_Pz_C1005_a-1*I_NAI_G3xy_Pz_C1005;
  abcd[236] = 2.0E0*I_NAI_I3x3z_Pz_C1005_a-2*I_NAI_G3xz_Pz_C1005;
  abcd[237] = 2.0E0*I_NAI_I2x3yz_Pz_C1005_a;
  abcd[238] = 2.0E0*I_NAI_I2x2y2z_Pz_C1005_a-1*I_NAI_G2x2y_Pz_C1005;
  abcd[239] = 2.0E0*I_NAI_I2xy3z_Pz_C1005_a-2*I_NAI_G2xyz_Pz_C1005;
  abcd[240] = 2.0E0*I_NAI_I2x4z_Pz_C1005_a-3*I_NAI_G2x2z_Pz_C1005;
  abcd[241] = 2.0E0*I_NAI_Ix4yz_Pz_C1005_a;
  abcd[242] = 2.0E0*I_NAI_Ix3y2z_Pz_C1005_a-1*I_NAI_Gx3y_Pz_C1005;
  abcd[243] = 2.0E0*I_NAI_Ix2y3z_Pz_C1005_a-2*I_NAI_Gx2yz_Pz_C1005;
  abcd[244] = 2.0E0*I_NAI_Ixy4z_Pz_C1005_a-3*I_NAI_Gxy2z_Pz_C1005;
  abcd[245] = 2.0E0*I_NAI_Ix5z_Pz_C1005_a-4*I_NAI_Gx3z_Pz_C1005;
  abcd[246] = 2.0E0*I_NAI_I5yz_Pz_C1005_a;
  abcd[247] = 2.0E0*I_NAI_I4y2z_Pz_C1005_a-1*I_NAI_G4y_Pz_C1005;
  abcd[248] = 2.0E0*I_NAI_I3y3z_Pz_C1005_a-2*I_NAI_G3yz_Pz_C1005;
  abcd[249] = 2.0E0*I_NAI_I2y4z_Pz_C1005_a-3*I_NAI_G2y2z_Pz_C1005;
  abcd[250] = 2.0E0*I_NAI_Iy5z_Pz_C1005_a-4*I_NAI_Gy3z_Pz_C1005;
  abcd[251] = 2.0E0*I_NAI_I6z_Pz_C1005_a-5*I_NAI_G4z_Pz_C1005;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C5_b
   ************************************************************/
  abcd[252] = 2.0E0*I_NAI_H5x_Px_C5_b;
  abcd[253] = 2.0E0*I_NAI_H4xy_Px_C5_b;
  abcd[254] = 2.0E0*I_NAI_H4xz_Px_C5_b;
  abcd[255] = 2.0E0*I_NAI_H3x2y_Px_C5_b;
  abcd[256] = 2.0E0*I_NAI_H3xyz_Px_C5_b;
  abcd[257] = 2.0E0*I_NAI_H3x2z_Px_C5_b;
  abcd[258] = 2.0E0*I_NAI_H2x3y_Px_C5_b;
  abcd[259] = 2.0E0*I_NAI_H2x2yz_Px_C5_b;
  abcd[260] = 2.0E0*I_NAI_H2xy2z_Px_C5_b;
  abcd[261] = 2.0E0*I_NAI_H2x3z_Px_C5_b;
  abcd[262] = 2.0E0*I_NAI_Hx4y_Px_C5_b;
  abcd[263] = 2.0E0*I_NAI_Hx3yz_Px_C5_b;
  abcd[264] = 2.0E0*I_NAI_Hx2y2z_Px_C5_b;
  abcd[265] = 2.0E0*I_NAI_Hxy3z_Px_C5_b;
  abcd[266] = 2.0E0*I_NAI_Hx4z_Px_C5_b;
  abcd[267] = 2.0E0*I_NAI_H5y_Px_C5_b;
  abcd[268] = 2.0E0*I_NAI_H4yz_Px_C5_b;
  abcd[269] = 2.0E0*I_NAI_H3y2z_Px_C5_b;
  abcd[270] = 2.0E0*I_NAI_H2y3z_Px_C5_b;
  abcd[271] = 2.0E0*I_NAI_Hy4z_Px_C5_b;
  abcd[272] = 2.0E0*I_NAI_H5z_Px_C5_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_C1005_b
   * RHS shell quartet name: SQ_NAI_H_S_C1005
   ************************************************************/
  abcd[273] = 2.0E0*I_NAI_H5x_D2x_C1005_b-1*I_NAI_H5x_S_C1005;
  abcd[274] = 2.0E0*I_NAI_H4xy_D2x_C1005_b-1*I_NAI_H4xy_S_C1005;
  abcd[275] = 2.0E0*I_NAI_H4xz_D2x_C1005_b-1*I_NAI_H4xz_S_C1005;
  abcd[276] = 2.0E0*I_NAI_H3x2y_D2x_C1005_b-1*I_NAI_H3x2y_S_C1005;
  abcd[277] = 2.0E0*I_NAI_H3xyz_D2x_C1005_b-1*I_NAI_H3xyz_S_C1005;
  abcd[278] = 2.0E0*I_NAI_H3x2z_D2x_C1005_b-1*I_NAI_H3x2z_S_C1005;
  abcd[279] = 2.0E0*I_NAI_H2x3y_D2x_C1005_b-1*I_NAI_H2x3y_S_C1005;
  abcd[280] = 2.0E0*I_NAI_H2x2yz_D2x_C1005_b-1*I_NAI_H2x2yz_S_C1005;
  abcd[281] = 2.0E0*I_NAI_H2xy2z_D2x_C1005_b-1*I_NAI_H2xy2z_S_C1005;
  abcd[282] = 2.0E0*I_NAI_H2x3z_D2x_C1005_b-1*I_NAI_H2x3z_S_C1005;
  abcd[283] = 2.0E0*I_NAI_Hx4y_D2x_C1005_b-1*I_NAI_Hx4y_S_C1005;
  abcd[284] = 2.0E0*I_NAI_Hx3yz_D2x_C1005_b-1*I_NAI_Hx3yz_S_C1005;
  abcd[285] = 2.0E0*I_NAI_Hx2y2z_D2x_C1005_b-1*I_NAI_Hx2y2z_S_C1005;
  abcd[286] = 2.0E0*I_NAI_Hxy3z_D2x_C1005_b-1*I_NAI_Hxy3z_S_C1005;
  abcd[287] = 2.0E0*I_NAI_Hx4z_D2x_C1005_b-1*I_NAI_Hx4z_S_C1005;
  abcd[288] = 2.0E0*I_NAI_H5y_D2x_C1005_b-1*I_NAI_H5y_S_C1005;
  abcd[289] = 2.0E0*I_NAI_H4yz_D2x_C1005_b-1*I_NAI_H4yz_S_C1005;
  abcd[290] = 2.0E0*I_NAI_H3y2z_D2x_C1005_b-1*I_NAI_H3y2z_S_C1005;
  abcd[291] = 2.0E0*I_NAI_H2y3z_D2x_C1005_b-1*I_NAI_H2y3z_S_C1005;
  abcd[292] = 2.0E0*I_NAI_Hy4z_D2x_C1005_b-1*I_NAI_Hy4z_S_C1005;
  abcd[293] = 2.0E0*I_NAI_H5z_D2x_C1005_b-1*I_NAI_H5z_S_C1005;
  abcd[294] = 2.0E0*I_NAI_H5x_Dxy_C1005_b;
  abcd[295] = 2.0E0*I_NAI_H4xy_Dxy_C1005_b;
  abcd[296] = 2.0E0*I_NAI_H4xz_Dxy_C1005_b;
  abcd[297] = 2.0E0*I_NAI_H3x2y_Dxy_C1005_b;
  abcd[298] = 2.0E0*I_NAI_H3xyz_Dxy_C1005_b;
  abcd[299] = 2.0E0*I_NAI_H3x2z_Dxy_C1005_b;
  abcd[300] = 2.0E0*I_NAI_H2x3y_Dxy_C1005_b;
  abcd[301] = 2.0E0*I_NAI_H2x2yz_Dxy_C1005_b;
  abcd[302] = 2.0E0*I_NAI_H2xy2z_Dxy_C1005_b;
  abcd[303] = 2.0E0*I_NAI_H2x3z_Dxy_C1005_b;
  abcd[304] = 2.0E0*I_NAI_Hx4y_Dxy_C1005_b;
  abcd[305] = 2.0E0*I_NAI_Hx3yz_Dxy_C1005_b;
  abcd[306] = 2.0E0*I_NAI_Hx2y2z_Dxy_C1005_b;
  abcd[307] = 2.0E0*I_NAI_Hxy3z_Dxy_C1005_b;
  abcd[308] = 2.0E0*I_NAI_Hx4z_Dxy_C1005_b;
  abcd[309] = 2.0E0*I_NAI_H5y_Dxy_C1005_b;
  abcd[310] = 2.0E0*I_NAI_H4yz_Dxy_C1005_b;
  abcd[311] = 2.0E0*I_NAI_H3y2z_Dxy_C1005_b;
  abcd[312] = 2.0E0*I_NAI_H2y3z_Dxy_C1005_b;
  abcd[313] = 2.0E0*I_NAI_Hy4z_Dxy_C1005_b;
  abcd[314] = 2.0E0*I_NAI_H5z_Dxy_C1005_b;
  abcd[315] = 2.0E0*I_NAI_H5x_Dxz_C1005_b;
  abcd[316] = 2.0E0*I_NAI_H4xy_Dxz_C1005_b;
  abcd[317] = 2.0E0*I_NAI_H4xz_Dxz_C1005_b;
  abcd[318] = 2.0E0*I_NAI_H3x2y_Dxz_C1005_b;
  abcd[319] = 2.0E0*I_NAI_H3xyz_Dxz_C1005_b;
  abcd[320] = 2.0E0*I_NAI_H3x2z_Dxz_C1005_b;
  abcd[321] = 2.0E0*I_NAI_H2x3y_Dxz_C1005_b;
  abcd[322] = 2.0E0*I_NAI_H2x2yz_Dxz_C1005_b;
  abcd[323] = 2.0E0*I_NAI_H2xy2z_Dxz_C1005_b;
  abcd[324] = 2.0E0*I_NAI_H2x3z_Dxz_C1005_b;
  abcd[325] = 2.0E0*I_NAI_Hx4y_Dxz_C1005_b;
  abcd[326] = 2.0E0*I_NAI_Hx3yz_Dxz_C1005_b;
  abcd[327] = 2.0E0*I_NAI_Hx2y2z_Dxz_C1005_b;
  abcd[328] = 2.0E0*I_NAI_Hxy3z_Dxz_C1005_b;
  abcd[329] = 2.0E0*I_NAI_Hx4z_Dxz_C1005_b;
  abcd[330] = 2.0E0*I_NAI_H5y_Dxz_C1005_b;
  abcd[331] = 2.0E0*I_NAI_H4yz_Dxz_C1005_b;
  abcd[332] = 2.0E0*I_NAI_H3y2z_Dxz_C1005_b;
  abcd[333] = 2.0E0*I_NAI_H2y3z_Dxz_C1005_b;
  abcd[334] = 2.0E0*I_NAI_Hy4z_Dxz_C1005_b;
  abcd[335] = 2.0E0*I_NAI_H5z_Dxz_C1005_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C5_b
   ************************************************************/
  abcd[336] = 2.0E0*I_NAI_H5x_Py_C5_b;
  abcd[337] = 2.0E0*I_NAI_H4xy_Py_C5_b;
  abcd[338] = 2.0E0*I_NAI_H4xz_Py_C5_b;
  abcd[339] = 2.0E0*I_NAI_H3x2y_Py_C5_b;
  abcd[340] = 2.0E0*I_NAI_H3xyz_Py_C5_b;
  abcd[341] = 2.0E0*I_NAI_H3x2z_Py_C5_b;
  abcd[342] = 2.0E0*I_NAI_H2x3y_Py_C5_b;
  abcd[343] = 2.0E0*I_NAI_H2x2yz_Py_C5_b;
  abcd[344] = 2.0E0*I_NAI_H2xy2z_Py_C5_b;
  abcd[345] = 2.0E0*I_NAI_H2x3z_Py_C5_b;
  abcd[346] = 2.0E0*I_NAI_Hx4y_Py_C5_b;
  abcd[347] = 2.0E0*I_NAI_Hx3yz_Py_C5_b;
  abcd[348] = 2.0E0*I_NAI_Hx2y2z_Py_C5_b;
  abcd[349] = 2.0E0*I_NAI_Hxy3z_Py_C5_b;
  abcd[350] = 2.0E0*I_NAI_Hx4z_Py_C5_b;
  abcd[351] = 2.0E0*I_NAI_H5y_Py_C5_b;
  abcd[352] = 2.0E0*I_NAI_H4yz_Py_C5_b;
  abcd[353] = 2.0E0*I_NAI_H3y2z_Py_C5_b;
  abcd[354] = 2.0E0*I_NAI_H2y3z_Py_C5_b;
  abcd[355] = 2.0E0*I_NAI_Hy4z_Py_C5_b;
  abcd[356] = 2.0E0*I_NAI_H5z_Py_C5_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_C1005_b
   * RHS shell quartet name: SQ_NAI_H_S_C1005
   ************************************************************/
  abcd[357] = 2.0E0*I_NAI_H5x_Dxy_C1005_b;
  abcd[358] = 2.0E0*I_NAI_H4xy_Dxy_C1005_b;
  abcd[359] = 2.0E0*I_NAI_H4xz_Dxy_C1005_b;
  abcd[360] = 2.0E0*I_NAI_H3x2y_Dxy_C1005_b;
  abcd[361] = 2.0E0*I_NAI_H3xyz_Dxy_C1005_b;
  abcd[362] = 2.0E0*I_NAI_H3x2z_Dxy_C1005_b;
  abcd[363] = 2.0E0*I_NAI_H2x3y_Dxy_C1005_b;
  abcd[364] = 2.0E0*I_NAI_H2x2yz_Dxy_C1005_b;
  abcd[365] = 2.0E0*I_NAI_H2xy2z_Dxy_C1005_b;
  abcd[366] = 2.0E0*I_NAI_H2x3z_Dxy_C1005_b;
  abcd[367] = 2.0E0*I_NAI_Hx4y_Dxy_C1005_b;
  abcd[368] = 2.0E0*I_NAI_Hx3yz_Dxy_C1005_b;
  abcd[369] = 2.0E0*I_NAI_Hx2y2z_Dxy_C1005_b;
  abcd[370] = 2.0E0*I_NAI_Hxy3z_Dxy_C1005_b;
  abcd[371] = 2.0E0*I_NAI_Hx4z_Dxy_C1005_b;
  abcd[372] = 2.0E0*I_NAI_H5y_Dxy_C1005_b;
  abcd[373] = 2.0E0*I_NAI_H4yz_Dxy_C1005_b;
  abcd[374] = 2.0E0*I_NAI_H3y2z_Dxy_C1005_b;
  abcd[375] = 2.0E0*I_NAI_H2y3z_Dxy_C1005_b;
  abcd[376] = 2.0E0*I_NAI_Hy4z_Dxy_C1005_b;
  abcd[377] = 2.0E0*I_NAI_H5z_Dxy_C1005_b;
  abcd[378] = 2.0E0*I_NAI_H5x_D2y_C1005_b-1*I_NAI_H5x_S_C1005;
  abcd[379] = 2.0E0*I_NAI_H4xy_D2y_C1005_b-1*I_NAI_H4xy_S_C1005;
  abcd[380] = 2.0E0*I_NAI_H4xz_D2y_C1005_b-1*I_NAI_H4xz_S_C1005;
  abcd[381] = 2.0E0*I_NAI_H3x2y_D2y_C1005_b-1*I_NAI_H3x2y_S_C1005;
  abcd[382] = 2.0E0*I_NAI_H3xyz_D2y_C1005_b-1*I_NAI_H3xyz_S_C1005;
  abcd[383] = 2.0E0*I_NAI_H3x2z_D2y_C1005_b-1*I_NAI_H3x2z_S_C1005;
  abcd[384] = 2.0E0*I_NAI_H2x3y_D2y_C1005_b-1*I_NAI_H2x3y_S_C1005;
  abcd[385] = 2.0E0*I_NAI_H2x2yz_D2y_C1005_b-1*I_NAI_H2x2yz_S_C1005;
  abcd[386] = 2.0E0*I_NAI_H2xy2z_D2y_C1005_b-1*I_NAI_H2xy2z_S_C1005;
  abcd[387] = 2.0E0*I_NAI_H2x3z_D2y_C1005_b-1*I_NAI_H2x3z_S_C1005;
  abcd[388] = 2.0E0*I_NAI_Hx4y_D2y_C1005_b-1*I_NAI_Hx4y_S_C1005;
  abcd[389] = 2.0E0*I_NAI_Hx3yz_D2y_C1005_b-1*I_NAI_Hx3yz_S_C1005;
  abcd[390] = 2.0E0*I_NAI_Hx2y2z_D2y_C1005_b-1*I_NAI_Hx2y2z_S_C1005;
  abcd[391] = 2.0E0*I_NAI_Hxy3z_D2y_C1005_b-1*I_NAI_Hxy3z_S_C1005;
  abcd[392] = 2.0E0*I_NAI_Hx4z_D2y_C1005_b-1*I_NAI_Hx4z_S_C1005;
  abcd[393] = 2.0E0*I_NAI_H5y_D2y_C1005_b-1*I_NAI_H5y_S_C1005;
  abcd[394] = 2.0E0*I_NAI_H4yz_D2y_C1005_b-1*I_NAI_H4yz_S_C1005;
  abcd[395] = 2.0E0*I_NAI_H3y2z_D2y_C1005_b-1*I_NAI_H3y2z_S_C1005;
  abcd[396] = 2.0E0*I_NAI_H2y3z_D2y_C1005_b-1*I_NAI_H2y3z_S_C1005;
  abcd[397] = 2.0E0*I_NAI_Hy4z_D2y_C1005_b-1*I_NAI_Hy4z_S_C1005;
  abcd[398] = 2.0E0*I_NAI_H5z_D2y_C1005_b-1*I_NAI_H5z_S_C1005;
  abcd[399] = 2.0E0*I_NAI_H5x_Dyz_C1005_b;
  abcd[400] = 2.0E0*I_NAI_H4xy_Dyz_C1005_b;
  abcd[401] = 2.0E0*I_NAI_H4xz_Dyz_C1005_b;
  abcd[402] = 2.0E0*I_NAI_H3x2y_Dyz_C1005_b;
  abcd[403] = 2.0E0*I_NAI_H3xyz_Dyz_C1005_b;
  abcd[404] = 2.0E0*I_NAI_H3x2z_Dyz_C1005_b;
  abcd[405] = 2.0E0*I_NAI_H2x3y_Dyz_C1005_b;
  abcd[406] = 2.0E0*I_NAI_H2x2yz_Dyz_C1005_b;
  abcd[407] = 2.0E0*I_NAI_H2xy2z_Dyz_C1005_b;
  abcd[408] = 2.0E0*I_NAI_H2x3z_Dyz_C1005_b;
  abcd[409] = 2.0E0*I_NAI_Hx4y_Dyz_C1005_b;
  abcd[410] = 2.0E0*I_NAI_Hx3yz_Dyz_C1005_b;
  abcd[411] = 2.0E0*I_NAI_Hx2y2z_Dyz_C1005_b;
  abcd[412] = 2.0E0*I_NAI_Hxy3z_Dyz_C1005_b;
  abcd[413] = 2.0E0*I_NAI_Hx4z_Dyz_C1005_b;
  abcd[414] = 2.0E0*I_NAI_H5y_Dyz_C1005_b;
  abcd[415] = 2.0E0*I_NAI_H4yz_Dyz_C1005_b;
  abcd[416] = 2.0E0*I_NAI_H3y2z_Dyz_C1005_b;
  abcd[417] = 2.0E0*I_NAI_H2y3z_Dyz_C1005_b;
  abcd[418] = 2.0E0*I_NAI_Hy4z_Dyz_C1005_b;
  abcd[419] = 2.0E0*I_NAI_H5z_Dyz_C1005_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_S_C5_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C5_b
   ************************************************************/
  abcd[420] = 2.0E0*I_NAI_H5x_Pz_C5_b;
  abcd[421] = 2.0E0*I_NAI_H4xy_Pz_C5_b;
  abcd[422] = 2.0E0*I_NAI_H4xz_Pz_C5_b;
  abcd[423] = 2.0E0*I_NAI_H3x2y_Pz_C5_b;
  abcd[424] = 2.0E0*I_NAI_H3xyz_Pz_C5_b;
  abcd[425] = 2.0E0*I_NAI_H3x2z_Pz_C5_b;
  abcd[426] = 2.0E0*I_NAI_H2x3y_Pz_C5_b;
  abcd[427] = 2.0E0*I_NAI_H2x2yz_Pz_C5_b;
  abcd[428] = 2.0E0*I_NAI_H2xy2z_Pz_C5_b;
  abcd[429] = 2.0E0*I_NAI_H2x3z_Pz_C5_b;
  abcd[430] = 2.0E0*I_NAI_Hx4y_Pz_C5_b;
  abcd[431] = 2.0E0*I_NAI_Hx3yz_Pz_C5_b;
  abcd[432] = 2.0E0*I_NAI_Hx2y2z_Pz_C5_b;
  abcd[433] = 2.0E0*I_NAI_Hxy3z_Pz_C5_b;
  abcd[434] = 2.0E0*I_NAI_Hx4z_Pz_C5_b;
  abcd[435] = 2.0E0*I_NAI_H5y_Pz_C5_b;
  abcd[436] = 2.0E0*I_NAI_H4yz_Pz_C5_b;
  abcd[437] = 2.0E0*I_NAI_H3y2z_Pz_C5_b;
  abcd[438] = 2.0E0*I_NAI_H2y3z_Pz_C5_b;
  abcd[439] = 2.0E0*I_NAI_Hy4z_Pz_C5_b;
  abcd[440] = 2.0E0*I_NAI_H5z_Pz_C5_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1005_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_C1005_b
   * RHS shell quartet name: SQ_NAI_H_S_C1005
   ************************************************************/
  abcd[441] = 2.0E0*I_NAI_H5x_Dxz_C1005_b;
  abcd[442] = 2.0E0*I_NAI_H4xy_Dxz_C1005_b;
  abcd[443] = 2.0E0*I_NAI_H4xz_Dxz_C1005_b;
  abcd[444] = 2.0E0*I_NAI_H3x2y_Dxz_C1005_b;
  abcd[445] = 2.0E0*I_NAI_H3xyz_Dxz_C1005_b;
  abcd[446] = 2.0E0*I_NAI_H3x2z_Dxz_C1005_b;
  abcd[447] = 2.0E0*I_NAI_H2x3y_Dxz_C1005_b;
  abcd[448] = 2.0E0*I_NAI_H2x2yz_Dxz_C1005_b;
  abcd[449] = 2.0E0*I_NAI_H2xy2z_Dxz_C1005_b;
  abcd[450] = 2.0E0*I_NAI_H2x3z_Dxz_C1005_b;
  abcd[451] = 2.0E0*I_NAI_Hx4y_Dxz_C1005_b;
  abcd[452] = 2.0E0*I_NAI_Hx3yz_Dxz_C1005_b;
  abcd[453] = 2.0E0*I_NAI_Hx2y2z_Dxz_C1005_b;
  abcd[454] = 2.0E0*I_NAI_Hxy3z_Dxz_C1005_b;
  abcd[455] = 2.0E0*I_NAI_Hx4z_Dxz_C1005_b;
  abcd[456] = 2.0E0*I_NAI_H5y_Dxz_C1005_b;
  abcd[457] = 2.0E0*I_NAI_H4yz_Dxz_C1005_b;
  abcd[458] = 2.0E0*I_NAI_H3y2z_Dxz_C1005_b;
  abcd[459] = 2.0E0*I_NAI_H2y3z_Dxz_C1005_b;
  abcd[460] = 2.0E0*I_NAI_Hy4z_Dxz_C1005_b;
  abcd[461] = 2.0E0*I_NAI_H5z_Dxz_C1005_b;
  abcd[462] = 2.0E0*I_NAI_H5x_Dyz_C1005_b;
  abcd[463] = 2.0E0*I_NAI_H4xy_Dyz_C1005_b;
  abcd[464] = 2.0E0*I_NAI_H4xz_Dyz_C1005_b;
  abcd[465] = 2.0E0*I_NAI_H3x2y_Dyz_C1005_b;
  abcd[466] = 2.0E0*I_NAI_H3xyz_Dyz_C1005_b;
  abcd[467] = 2.0E0*I_NAI_H3x2z_Dyz_C1005_b;
  abcd[468] = 2.0E0*I_NAI_H2x3y_Dyz_C1005_b;
  abcd[469] = 2.0E0*I_NAI_H2x2yz_Dyz_C1005_b;
  abcd[470] = 2.0E0*I_NAI_H2xy2z_Dyz_C1005_b;
  abcd[471] = 2.0E0*I_NAI_H2x3z_Dyz_C1005_b;
  abcd[472] = 2.0E0*I_NAI_Hx4y_Dyz_C1005_b;
  abcd[473] = 2.0E0*I_NAI_Hx3yz_Dyz_C1005_b;
  abcd[474] = 2.0E0*I_NAI_Hx2y2z_Dyz_C1005_b;
  abcd[475] = 2.0E0*I_NAI_Hxy3z_Dyz_C1005_b;
  abcd[476] = 2.0E0*I_NAI_Hx4z_Dyz_C1005_b;
  abcd[477] = 2.0E0*I_NAI_H5y_Dyz_C1005_b;
  abcd[478] = 2.0E0*I_NAI_H4yz_Dyz_C1005_b;
  abcd[479] = 2.0E0*I_NAI_H3y2z_Dyz_C1005_b;
  abcd[480] = 2.0E0*I_NAI_H2y3z_Dyz_C1005_b;
  abcd[481] = 2.0E0*I_NAI_Hy4z_Dyz_C1005_b;
  abcd[482] = 2.0E0*I_NAI_H5z_Dyz_C1005_b;
  abcd[483] = 2.0E0*I_NAI_H5x_D2z_C1005_b-1*I_NAI_H5x_S_C1005;
  abcd[484] = 2.0E0*I_NAI_H4xy_D2z_C1005_b-1*I_NAI_H4xy_S_C1005;
  abcd[485] = 2.0E0*I_NAI_H4xz_D2z_C1005_b-1*I_NAI_H4xz_S_C1005;
  abcd[486] = 2.0E0*I_NAI_H3x2y_D2z_C1005_b-1*I_NAI_H3x2y_S_C1005;
  abcd[487] = 2.0E0*I_NAI_H3xyz_D2z_C1005_b-1*I_NAI_H3xyz_S_C1005;
  abcd[488] = 2.0E0*I_NAI_H3x2z_D2z_C1005_b-1*I_NAI_H3x2z_S_C1005;
  abcd[489] = 2.0E0*I_NAI_H2x3y_D2z_C1005_b-1*I_NAI_H2x3y_S_C1005;
  abcd[490] = 2.0E0*I_NAI_H2x2yz_D2z_C1005_b-1*I_NAI_H2x2yz_S_C1005;
  abcd[491] = 2.0E0*I_NAI_H2xy2z_D2z_C1005_b-1*I_NAI_H2xy2z_S_C1005;
  abcd[492] = 2.0E0*I_NAI_H2x3z_D2z_C1005_b-1*I_NAI_H2x3z_S_C1005;
  abcd[493] = 2.0E0*I_NAI_Hx4y_D2z_C1005_b-1*I_NAI_Hx4y_S_C1005;
  abcd[494] = 2.0E0*I_NAI_Hx3yz_D2z_C1005_b-1*I_NAI_Hx3yz_S_C1005;
  abcd[495] = 2.0E0*I_NAI_Hx2y2z_D2z_C1005_b-1*I_NAI_Hx2y2z_S_C1005;
  abcd[496] = 2.0E0*I_NAI_Hxy3z_D2z_C1005_b-1*I_NAI_Hxy3z_S_C1005;
  abcd[497] = 2.0E0*I_NAI_Hx4z_D2z_C1005_b-1*I_NAI_Hx4z_S_C1005;
  abcd[498] = 2.0E0*I_NAI_H5y_D2z_C1005_b-1*I_NAI_H5y_S_C1005;
  abcd[499] = 2.0E0*I_NAI_H4yz_D2z_C1005_b-1*I_NAI_H4yz_S_C1005;
  abcd[500] = 2.0E0*I_NAI_H3y2z_D2z_C1005_b-1*I_NAI_H3y2z_S_C1005;
  abcd[501] = 2.0E0*I_NAI_H2y3z_D2z_C1005_b-1*I_NAI_H2y3z_S_C1005;
  abcd[502] = 2.0E0*I_NAI_Hy4z_D2z_C1005_b-1*I_NAI_Hy4z_S_C1005;
  abcd[503] = 2.0E0*I_NAI_H5z_D2z_C1005_b-1*I_NAI_H5z_S_C1005;
}
