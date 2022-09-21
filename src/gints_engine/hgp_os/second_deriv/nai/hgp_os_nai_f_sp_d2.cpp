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

void hgp_os_nai_f_sp_d2(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_H5x_S_C3_aa = 0.0E0;
  Double I_NAI_H4xy_S_C3_aa = 0.0E0;
  Double I_NAI_H4xz_S_C3_aa = 0.0E0;
  Double I_NAI_H3x2y_S_C3_aa = 0.0E0;
  Double I_NAI_H3xyz_S_C3_aa = 0.0E0;
  Double I_NAI_H3x2z_S_C3_aa = 0.0E0;
  Double I_NAI_H2x3y_S_C3_aa = 0.0E0;
  Double I_NAI_H2x2yz_S_C3_aa = 0.0E0;
  Double I_NAI_H2xy2z_S_C3_aa = 0.0E0;
  Double I_NAI_H2x3z_S_C3_aa = 0.0E0;
  Double I_NAI_Hx4y_S_C3_aa = 0.0E0;
  Double I_NAI_Hx3yz_S_C3_aa = 0.0E0;
  Double I_NAI_Hx2y2z_S_C3_aa = 0.0E0;
  Double I_NAI_Hxy3z_S_C3_aa = 0.0E0;
  Double I_NAI_Hx4z_S_C3_aa = 0.0E0;
  Double I_NAI_H5y_S_C3_aa = 0.0E0;
  Double I_NAI_H4yz_S_C3_aa = 0.0E0;
  Double I_NAI_H3y2z_S_C3_aa = 0.0E0;
  Double I_NAI_H2y3z_S_C3_aa = 0.0E0;
  Double I_NAI_Hy4z_S_C3_aa = 0.0E0;
  Double I_NAI_H5z_S_C3_aa = 0.0E0;
  Double I_NAI_F3x_S_C3_a = 0.0E0;
  Double I_NAI_F2xy_S_C3_a = 0.0E0;
  Double I_NAI_F2xz_S_C3_a = 0.0E0;
  Double I_NAI_Fx2y_S_C3_a = 0.0E0;
  Double I_NAI_Fxyz_S_C3_a = 0.0E0;
  Double I_NAI_Fx2z_S_C3_a = 0.0E0;
  Double I_NAI_F3y_S_C3_a = 0.0E0;
  Double I_NAI_F2yz_S_C3_a = 0.0E0;
  Double I_NAI_Fy2z_S_C3_a = 0.0E0;
  Double I_NAI_F3z_S_C3_a = 0.0E0;
  Double I_NAI_Px_S_C3 = 0.0E0;
  Double I_NAI_Py_S_C3 = 0.0E0;
  Double I_NAI_Pz_S_C3 = 0.0E0;
  Double I_NAI_G4x_S_C1003_a = 0.0E0;
  Double I_NAI_G3xy_S_C1003_a = 0.0E0;
  Double I_NAI_G3xz_S_C1003_a = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_a = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_a = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_a = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_a = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_a = 0.0E0;
  Double I_NAI_G4y_S_C1003_a = 0.0E0;
  Double I_NAI_G3yz_S_C1003_a = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_a = 0.0E0;
  Double I_NAI_G4z_S_C1003_a = 0.0E0;
  Double I_NAI_D2x_S_C1003 = 0.0E0;
  Double I_NAI_Dxy_S_C1003 = 0.0E0;
  Double I_NAI_Dxz_S_C1003 = 0.0E0;
  Double I_NAI_D2y_S_C1003 = 0.0E0;
  Double I_NAI_Dyz_S_C1003 = 0.0E0;
  Double I_NAI_D2z_S_C1003 = 0.0E0;
  Double I_NAI_F3x_S_C3_b = 0.0E0;
  Double I_NAI_F2xy_S_C3_b = 0.0E0;
  Double I_NAI_F2xz_S_C3_b = 0.0E0;
  Double I_NAI_Fx2y_S_C3_b = 0.0E0;
  Double I_NAI_Fxyz_S_C3_b = 0.0E0;
  Double I_NAI_Fx2z_S_C3_b = 0.0E0;
  Double I_NAI_F3y_S_C3_b = 0.0E0;
  Double I_NAI_F2yz_S_C3_b = 0.0E0;
  Double I_NAI_Fy2z_S_C3_b = 0.0E0;
  Double I_NAI_F3z_S_C3_b = 0.0E0;
  Double I_NAI_I6x_S_C1003_aa = 0.0E0;
  Double I_NAI_I5xy_S_C1003_aa = 0.0E0;
  Double I_NAI_I5xz_S_C1003_aa = 0.0E0;
  Double I_NAI_I4x2y_S_C1003_aa = 0.0E0;
  Double I_NAI_I4xyz_S_C1003_aa = 0.0E0;
  Double I_NAI_I4x2z_S_C1003_aa = 0.0E0;
  Double I_NAI_I3x3y_S_C1003_aa = 0.0E0;
  Double I_NAI_I3x2yz_S_C1003_aa = 0.0E0;
  Double I_NAI_I3xy2z_S_C1003_aa = 0.0E0;
  Double I_NAI_I3x3z_S_C1003_aa = 0.0E0;
  Double I_NAI_I2x4y_S_C1003_aa = 0.0E0;
  Double I_NAI_I2x3yz_S_C1003_aa = 0.0E0;
  Double I_NAI_I2x2y2z_S_C1003_aa = 0.0E0;
  Double I_NAI_I2xy3z_S_C1003_aa = 0.0E0;
  Double I_NAI_I2x4z_S_C1003_aa = 0.0E0;
  Double I_NAI_Ix5y_S_C1003_aa = 0.0E0;
  Double I_NAI_Ix4yz_S_C1003_aa = 0.0E0;
  Double I_NAI_Ix3y2z_S_C1003_aa = 0.0E0;
  Double I_NAI_Ix2y3z_S_C1003_aa = 0.0E0;
  Double I_NAI_Ixy4z_S_C1003_aa = 0.0E0;
  Double I_NAI_Ix5z_S_C1003_aa = 0.0E0;
  Double I_NAI_I6y_S_C1003_aa = 0.0E0;
  Double I_NAI_I5yz_S_C1003_aa = 0.0E0;
  Double I_NAI_I4y2z_S_C1003_aa = 0.0E0;
  Double I_NAI_I3y3z_S_C1003_aa = 0.0E0;
  Double I_NAI_I2y4z_S_C1003_aa = 0.0E0;
  Double I_NAI_Iy5z_S_C1003_aa = 0.0E0;
  Double I_NAI_I6z_S_C1003_aa = 0.0E0;
  Double I_NAI_H5x_S_C1003_aa = 0.0E0;
  Double I_NAI_H4xy_S_C1003_aa = 0.0E0;
  Double I_NAI_H4xz_S_C1003_aa = 0.0E0;
  Double I_NAI_H3x2y_S_C1003_aa = 0.0E0;
  Double I_NAI_H3xyz_S_C1003_aa = 0.0E0;
  Double I_NAI_H3x2z_S_C1003_aa = 0.0E0;
  Double I_NAI_H2x3y_S_C1003_aa = 0.0E0;
  Double I_NAI_H2x2yz_S_C1003_aa = 0.0E0;
  Double I_NAI_H2xy2z_S_C1003_aa = 0.0E0;
  Double I_NAI_H2x3z_S_C1003_aa = 0.0E0;
  Double I_NAI_Hx4y_S_C1003_aa = 0.0E0;
  Double I_NAI_Hx3yz_S_C1003_aa = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1003_aa = 0.0E0;
  Double I_NAI_Hxy3z_S_C1003_aa = 0.0E0;
  Double I_NAI_Hx4z_S_C1003_aa = 0.0E0;
  Double I_NAI_H5y_S_C1003_aa = 0.0E0;
  Double I_NAI_H4yz_S_C1003_aa = 0.0E0;
  Double I_NAI_H3y2z_S_C1003_aa = 0.0E0;
  Double I_NAI_H2y3z_S_C1003_aa = 0.0E0;
  Double I_NAI_Hy4z_S_C1003_aa = 0.0E0;
  Double I_NAI_H5z_S_C1003_aa = 0.0E0;
  Double I_NAI_F3x_S_C1003_a = 0.0E0;
  Double I_NAI_F2xy_S_C1003_a = 0.0E0;
  Double I_NAI_F2xz_S_C1003_a = 0.0E0;
  Double I_NAI_Fx2y_S_C1003_a = 0.0E0;
  Double I_NAI_Fxyz_S_C1003_a = 0.0E0;
  Double I_NAI_Fx2z_S_C1003_a = 0.0E0;
  Double I_NAI_F3y_S_C1003_a = 0.0E0;
  Double I_NAI_F2yz_S_C1003_a = 0.0E0;
  Double I_NAI_Fy2z_S_C1003_a = 0.0E0;
  Double I_NAI_F3z_S_C1003_a = 0.0E0;
  Double I_NAI_Px_S_C1003 = 0.0E0;
  Double I_NAI_Py_S_C1003 = 0.0E0;
  Double I_NAI_Pz_S_C1003 = 0.0E0;
  Double I_NAI_H5x_S_C3_ab = 0.0E0;
  Double I_NAI_H4xy_S_C3_ab = 0.0E0;
  Double I_NAI_H4xz_S_C3_ab = 0.0E0;
  Double I_NAI_H3x2y_S_C3_ab = 0.0E0;
  Double I_NAI_H3xyz_S_C3_ab = 0.0E0;
  Double I_NAI_H3x2z_S_C3_ab = 0.0E0;
  Double I_NAI_H2x3y_S_C3_ab = 0.0E0;
  Double I_NAI_H2x2yz_S_C3_ab = 0.0E0;
  Double I_NAI_H2xy2z_S_C3_ab = 0.0E0;
  Double I_NAI_H2x3z_S_C3_ab = 0.0E0;
  Double I_NAI_Hx4y_S_C3_ab = 0.0E0;
  Double I_NAI_Hx3yz_S_C3_ab = 0.0E0;
  Double I_NAI_Hx2y2z_S_C3_ab = 0.0E0;
  Double I_NAI_Hxy3z_S_C3_ab = 0.0E0;
  Double I_NAI_Hx4z_S_C3_ab = 0.0E0;
  Double I_NAI_H5y_S_C3_ab = 0.0E0;
  Double I_NAI_H4yz_S_C3_ab = 0.0E0;
  Double I_NAI_H3y2z_S_C3_ab = 0.0E0;
  Double I_NAI_H2y3z_S_C3_ab = 0.0E0;
  Double I_NAI_Hy4z_S_C3_ab = 0.0E0;
  Double I_NAI_H5z_S_C3_ab = 0.0E0;
  Double I_NAI_G4x_S_C3_ab = 0.0E0;
  Double I_NAI_G3xy_S_C3_ab = 0.0E0;
  Double I_NAI_G3xz_S_C3_ab = 0.0E0;
  Double I_NAI_G2x2y_S_C3_ab = 0.0E0;
  Double I_NAI_G2xyz_S_C3_ab = 0.0E0;
  Double I_NAI_G2x2z_S_C3_ab = 0.0E0;
  Double I_NAI_Gx3y_S_C3_ab = 0.0E0;
  Double I_NAI_Gx2yz_S_C3_ab = 0.0E0;
  Double I_NAI_Gxy2z_S_C3_ab = 0.0E0;
  Double I_NAI_Gx3z_S_C3_ab = 0.0E0;
  Double I_NAI_G4y_S_C3_ab = 0.0E0;
  Double I_NAI_G3yz_S_C3_ab = 0.0E0;
  Double I_NAI_G2y2z_S_C3_ab = 0.0E0;
  Double I_NAI_Gy3z_S_C3_ab = 0.0E0;
  Double I_NAI_G4z_S_C3_ab = 0.0E0;
  Double I_NAI_D2x_S_C3_b = 0.0E0;
  Double I_NAI_Dxy_S_C3_b = 0.0E0;
  Double I_NAI_Dxz_S_C3_b = 0.0E0;
  Double I_NAI_D2y_S_C3_b = 0.0E0;
  Double I_NAI_Dyz_S_C3_b = 0.0E0;
  Double I_NAI_D2z_S_C3_b = 0.0E0;
  Double I_NAI_G4x_S_C1003_b = 0.0E0;
  Double I_NAI_G3xy_S_C1003_b = 0.0E0;
  Double I_NAI_G3xz_S_C1003_b = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_b = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_b = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_b = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_b = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_b = 0.0E0;
  Double I_NAI_G4y_S_C1003_b = 0.0E0;
  Double I_NAI_G3yz_S_C1003_b = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_b = 0.0E0;
  Double I_NAI_G4z_S_C1003_b = 0.0E0;
  Double I_NAI_F3x_S_C1003_b = 0.0E0;
  Double I_NAI_F2xy_S_C1003_b = 0.0E0;
  Double I_NAI_F2xz_S_C1003_b = 0.0E0;
  Double I_NAI_Fx2y_S_C1003_b = 0.0E0;
  Double I_NAI_Fxyz_S_C1003_b = 0.0E0;
  Double I_NAI_Fx2z_S_C1003_b = 0.0E0;
  Double I_NAI_F3y_S_C1003_b = 0.0E0;
  Double I_NAI_F2yz_S_C1003_b = 0.0E0;
  Double I_NAI_Fy2z_S_C1003_b = 0.0E0;
  Double I_NAI_F3z_S_C1003_b = 0.0E0;
  Double I_NAI_I6x_S_C1003_ab = 0.0E0;
  Double I_NAI_I5xy_S_C1003_ab = 0.0E0;
  Double I_NAI_I5xz_S_C1003_ab = 0.0E0;
  Double I_NAI_I4x2y_S_C1003_ab = 0.0E0;
  Double I_NAI_I4xyz_S_C1003_ab = 0.0E0;
  Double I_NAI_I4x2z_S_C1003_ab = 0.0E0;
  Double I_NAI_I3x3y_S_C1003_ab = 0.0E0;
  Double I_NAI_I3x2yz_S_C1003_ab = 0.0E0;
  Double I_NAI_I3xy2z_S_C1003_ab = 0.0E0;
  Double I_NAI_I3x3z_S_C1003_ab = 0.0E0;
  Double I_NAI_I2x4y_S_C1003_ab = 0.0E0;
  Double I_NAI_I2x3yz_S_C1003_ab = 0.0E0;
  Double I_NAI_I2x2y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_I2xy3z_S_C1003_ab = 0.0E0;
  Double I_NAI_I2x4z_S_C1003_ab = 0.0E0;
  Double I_NAI_Ix5y_S_C1003_ab = 0.0E0;
  Double I_NAI_Ix4yz_S_C1003_ab = 0.0E0;
  Double I_NAI_Ix3y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_Ix2y3z_S_C1003_ab = 0.0E0;
  Double I_NAI_Ixy4z_S_C1003_ab = 0.0E0;
  Double I_NAI_Ix5z_S_C1003_ab = 0.0E0;
  Double I_NAI_I6y_S_C1003_ab = 0.0E0;
  Double I_NAI_I5yz_S_C1003_ab = 0.0E0;
  Double I_NAI_I4y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_I3y3z_S_C1003_ab = 0.0E0;
  Double I_NAI_I2y4z_S_C1003_ab = 0.0E0;
  Double I_NAI_Iy5z_S_C1003_ab = 0.0E0;
  Double I_NAI_I6z_S_C1003_ab = 0.0E0;
  Double I_NAI_H5x_S_C1003_ab = 0.0E0;
  Double I_NAI_H4xy_S_C1003_ab = 0.0E0;
  Double I_NAI_H4xz_S_C1003_ab = 0.0E0;
  Double I_NAI_H3x2y_S_C1003_ab = 0.0E0;
  Double I_NAI_H3xyz_S_C1003_ab = 0.0E0;
  Double I_NAI_H3x2z_S_C1003_ab = 0.0E0;
  Double I_NAI_H2x3y_S_C1003_ab = 0.0E0;
  Double I_NAI_H2x2yz_S_C1003_ab = 0.0E0;
  Double I_NAI_H2xy2z_S_C1003_ab = 0.0E0;
  Double I_NAI_H2x3z_S_C1003_ab = 0.0E0;
  Double I_NAI_Hx4y_S_C1003_ab = 0.0E0;
  Double I_NAI_Hx3yz_S_C1003_ab = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_Hxy3z_S_C1003_ab = 0.0E0;
  Double I_NAI_Hx4z_S_C1003_ab = 0.0E0;
  Double I_NAI_H5y_S_C1003_ab = 0.0E0;
  Double I_NAI_H4yz_S_C1003_ab = 0.0E0;
  Double I_NAI_H3y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_H2y3z_S_C1003_ab = 0.0E0;
  Double I_NAI_Hy4z_S_C1003_ab = 0.0E0;
  Double I_NAI_H5z_S_C1003_ab = 0.0E0;
  Double I_NAI_G4x_S_C1003_ab = 0.0E0;
  Double I_NAI_G3xy_S_C1003_ab = 0.0E0;
  Double I_NAI_G3xz_S_C1003_ab = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_ab = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_ab = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_ab = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_ab = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_ab = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_ab = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_ab = 0.0E0;
  Double I_NAI_G4y_S_C1003_ab = 0.0E0;
  Double I_NAI_G3yz_S_C1003_ab = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_ab = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_ab = 0.0E0;
  Double I_NAI_G4z_S_C1003_ab = 0.0E0;
  Double I_NAI_D2x_S_C1003_b = 0.0E0;
  Double I_NAI_Dxy_S_C1003_b = 0.0E0;
  Double I_NAI_Dxz_S_C1003_b = 0.0E0;
  Double I_NAI_D2y_S_C1003_b = 0.0E0;
  Double I_NAI_Dyz_S_C1003_b = 0.0E0;
  Double I_NAI_D2z_S_C1003_b = 0.0E0;
  Double I_NAI_H5x_S_C3_bb = 0.0E0;
  Double I_NAI_H4xy_S_C3_bb = 0.0E0;
  Double I_NAI_H4xz_S_C3_bb = 0.0E0;
  Double I_NAI_H3x2y_S_C3_bb = 0.0E0;
  Double I_NAI_H3xyz_S_C3_bb = 0.0E0;
  Double I_NAI_H3x2z_S_C3_bb = 0.0E0;
  Double I_NAI_H2x3y_S_C3_bb = 0.0E0;
  Double I_NAI_H2x2yz_S_C3_bb = 0.0E0;
  Double I_NAI_H2xy2z_S_C3_bb = 0.0E0;
  Double I_NAI_H2x3z_S_C3_bb = 0.0E0;
  Double I_NAI_Hx4y_S_C3_bb = 0.0E0;
  Double I_NAI_Hx3yz_S_C3_bb = 0.0E0;
  Double I_NAI_Hx2y2z_S_C3_bb = 0.0E0;
  Double I_NAI_Hxy3z_S_C3_bb = 0.0E0;
  Double I_NAI_Hx4z_S_C3_bb = 0.0E0;
  Double I_NAI_H5y_S_C3_bb = 0.0E0;
  Double I_NAI_H4yz_S_C3_bb = 0.0E0;
  Double I_NAI_H3y2z_S_C3_bb = 0.0E0;
  Double I_NAI_H2y3z_S_C3_bb = 0.0E0;
  Double I_NAI_Hy4z_S_C3_bb = 0.0E0;
  Double I_NAI_H5z_S_C3_bb = 0.0E0;
  Double I_NAI_G4x_S_C3_bb = 0.0E0;
  Double I_NAI_G3xy_S_C3_bb = 0.0E0;
  Double I_NAI_G3xz_S_C3_bb = 0.0E0;
  Double I_NAI_G2x2y_S_C3_bb = 0.0E0;
  Double I_NAI_G2xyz_S_C3_bb = 0.0E0;
  Double I_NAI_G2x2z_S_C3_bb = 0.0E0;
  Double I_NAI_Gx3y_S_C3_bb = 0.0E0;
  Double I_NAI_Gx2yz_S_C3_bb = 0.0E0;
  Double I_NAI_Gxy2z_S_C3_bb = 0.0E0;
  Double I_NAI_Gx3z_S_C3_bb = 0.0E0;
  Double I_NAI_G4y_S_C3_bb = 0.0E0;
  Double I_NAI_G3yz_S_C3_bb = 0.0E0;
  Double I_NAI_G2y2z_S_C3_bb = 0.0E0;
  Double I_NAI_Gy3z_S_C3_bb = 0.0E0;
  Double I_NAI_G4z_S_C3_bb = 0.0E0;
  Double I_NAI_F3x_S_C3_bb = 0.0E0;
  Double I_NAI_F2xy_S_C3_bb = 0.0E0;
  Double I_NAI_F2xz_S_C3_bb = 0.0E0;
  Double I_NAI_Fx2y_S_C3_bb = 0.0E0;
  Double I_NAI_Fxyz_S_C3_bb = 0.0E0;
  Double I_NAI_Fx2z_S_C3_bb = 0.0E0;
  Double I_NAI_F3y_S_C3_bb = 0.0E0;
  Double I_NAI_F2yz_S_C3_bb = 0.0E0;
  Double I_NAI_Fy2z_S_C3_bb = 0.0E0;
  Double I_NAI_F3z_S_C3_bb = 0.0E0;
  Double I_NAI_I6x_S_C1003_bb = 0.0E0;
  Double I_NAI_I5xy_S_C1003_bb = 0.0E0;
  Double I_NAI_I5xz_S_C1003_bb = 0.0E0;
  Double I_NAI_I4x2y_S_C1003_bb = 0.0E0;
  Double I_NAI_I4xyz_S_C1003_bb = 0.0E0;
  Double I_NAI_I4x2z_S_C1003_bb = 0.0E0;
  Double I_NAI_I3x3y_S_C1003_bb = 0.0E0;
  Double I_NAI_I3x2yz_S_C1003_bb = 0.0E0;
  Double I_NAI_I3xy2z_S_C1003_bb = 0.0E0;
  Double I_NAI_I3x3z_S_C1003_bb = 0.0E0;
  Double I_NAI_I2x4y_S_C1003_bb = 0.0E0;
  Double I_NAI_I2x3yz_S_C1003_bb = 0.0E0;
  Double I_NAI_I2x2y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_I2xy3z_S_C1003_bb = 0.0E0;
  Double I_NAI_I2x4z_S_C1003_bb = 0.0E0;
  Double I_NAI_Ix5y_S_C1003_bb = 0.0E0;
  Double I_NAI_Ix4yz_S_C1003_bb = 0.0E0;
  Double I_NAI_Ix3y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_Ix2y3z_S_C1003_bb = 0.0E0;
  Double I_NAI_Ixy4z_S_C1003_bb = 0.0E0;
  Double I_NAI_Ix5z_S_C1003_bb = 0.0E0;
  Double I_NAI_I6y_S_C1003_bb = 0.0E0;
  Double I_NAI_I5yz_S_C1003_bb = 0.0E0;
  Double I_NAI_I4y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_I3y3z_S_C1003_bb = 0.0E0;
  Double I_NAI_I2y4z_S_C1003_bb = 0.0E0;
  Double I_NAI_Iy5z_S_C1003_bb = 0.0E0;
  Double I_NAI_I6z_S_C1003_bb = 0.0E0;
  Double I_NAI_H5x_S_C1003_bb = 0.0E0;
  Double I_NAI_H4xy_S_C1003_bb = 0.0E0;
  Double I_NAI_H4xz_S_C1003_bb = 0.0E0;
  Double I_NAI_H3x2y_S_C1003_bb = 0.0E0;
  Double I_NAI_H3xyz_S_C1003_bb = 0.0E0;
  Double I_NAI_H3x2z_S_C1003_bb = 0.0E0;
  Double I_NAI_H2x3y_S_C1003_bb = 0.0E0;
  Double I_NAI_H2x2yz_S_C1003_bb = 0.0E0;
  Double I_NAI_H2xy2z_S_C1003_bb = 0.0E0;
  Double I_NAI_H2x3z_S_C1003_bb = 0.0E0;
  Double I_NAI_Hx4y_S_C1003_bb = 0.0E0;
  Double I_NAI_Hx3yz_S_C1003_bb = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_Hxy3z_S_C1003_bb = 0.0E0;
  Double I_NAI_Hx4z_S_C1003_bb = 0.0E0;
  Double I_NAI_H5y_S_C1003_bb = 0.0E0;
  Double I_NAI_H4yz_S_C1003_bb = 0.0E0;
  Double I_NAI_H3y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_H2y3z_S_C1003_bb = 0.0E0;
  Double I_NAI_Hy4z_S_C1003_bb = 0.0E0;
  Double I_NAI_H5z_S_C1003_bb = 0.0E0;
  Double I_NAI_G4x_S_C1003_bb = 0.0E0;
  Double I_NAI_G3xy_S_C1003_bb = 0.0E0;
  Double I_NAI_G3xz_S_C1003_bb = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_bb = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_bb = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_bb = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_bb = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_bb = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_bb = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_bb = 0.0E0;
  Double I_NAI_G4y_S_C1003_bb = 0.0E0;
  Double I_NAI_G3yz_S_C1003_bb = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_bb = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_bb = 0.0E0;
  Double I_NAI_G4z_S_C1003_bb = 0.0E0;
  Double I_NAI_F3x_S_C1003_bb = 0.0E0;
  Double I_NAI_F2xy_S_C1003_bb = 0.0E0;
  Double I_NAI_F2xz_S_C1003_bb = 0.0E0;
  Double I_NAI_Fx2y_S_C1003_bb = 0.0E0;
  Double I_NAI_Fxyz_S_C1003_bb = 0.0E0;
  Double I_NAI_Fx2z_S_C1003_bb = 0.0E0;
  Double I_NAI_F3y_S_C1003_bb = 0.0E0;
  Double I_NAI_F2yz_S_C1003_bb = 0.0E0;
  Double I_NAI_Fy2z_S_C1003_bb = 0.0E0;
  Double I_NAI_F3z_S_C1003_bb = 0.0E0;

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

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
      double I_NAI_S_S_M5_vrr_d  = 0.0E0;
      double I_NAI_S_S_M6_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER47;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = 1.0E0+u2*ONEOVER15*I_NAI_S_S_M6_vrr;
        I_NAI_S_S_M6_vrr = ONEOVER13*I_NAI_S_S_M6_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M6_vrr  = f*I_NAI_S_S_M6_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);
        I_NAI_S_S_M5_vrr = static_cast<Double>(I_NAI_S_S_M5_vrr_d);
        I_NAI_S_S_M6_vrr = static_cast<Double>(I_NAI_S_S_M6_vrr_d);

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

#endif

      }


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
       * shell quartet name: SQ_NAI_H_S_C3_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C3_aa_coefs = ic2*alpha*alpha;
      I_NAI_H5x_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C3_aa += SQ_NAI_H_S_C3_aa_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C3_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C3_a_coefs = ic2*alpha;
      I_NAI_F3x_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C3_a += SQ_NAI_F_S_C3_a_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C3
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C3_coefs = ic2;
      I_NAI_Px_S_C3 += SQ_NAI_P_S_C3_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C3 += SQ_NAI_P_S_C3_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C3 += SQ_NAI_P_S_C3_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_a_coefs = ic2_1*alpha;
      I_NAI_G4x_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1003_coefs = ic2_1;
      I_NAI_D2x_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C3_b_coefs = ic2*beta;
      I_NAI_F3x_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C1003_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C1003_aa_coefs = ic2_1*alpha*alpha;
      I_NAI_I6x_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C1003_aa += SQ_NAI_I_S_C1003_aa_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1003_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1003_aa_coefs = ic2_1*alpha*alpha;
      I_NAI_H5x_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1003_aa += SQ_NAI_H_S_C1003_aa_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1003_a_coefs = ic2_1*alpha;
      I_NAI_F3x_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1003_a += SQ_NAI_F_S_C1003_a_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1003_coefs = ic2_1;
      I_NAI_Px_S_C1003 += SQ_NAI_P_S_C1003_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1003 += SQ_NAI_P_S_C1003_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1003 += SQ_NAI_P_S_C1003_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C3_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C3_ab_coefs = ic2*alpha*beta;
      I_NAI_H5x_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C3_ab += SQ_NAI_H_S_C3_ab_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C3_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C3_ab_coefs = ic2*alpha*beta;
      I_NAI_G4x_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C3_ab += SQ_NAI_G_S_C3_ab_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C3_b_coefs = ic2*beta;
      I_NAI_D2x_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C3_b += SQ_NAI_D_S_C3_b_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_G4x_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_F3x_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C1003_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C1003_ab_coefs = ic2_1*alpha*beta;
      I_NAI_I6x_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C1003_ab += SQ_NAI_I_S_C1003_ab_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1003_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1003_ab_coefs = ic2_1*alpha*beta;
      I_NAI_H5x_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1003_ab += SQ_NAI_H_S_C1003_ab_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_ab_coefs = ic2_1*alpha*beta;
      I_NAI_G4x_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_ab += SQ_NAI_G_S_C1003_ab_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_D2x_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1003_b += SQ_NAI_D_S_C1003_b_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C3_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C3_bb_coefs = ic2*beta*beta;
      I_NAI_H5x_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C3_bb += SQ_NAI_H_S_C3_bb_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C3_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C3_bb_coefs = ic2*beta*beta;
      I_NAI_G4x_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C3_bb += SQ_NAI_G_S_C3_bb_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C3_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C3_bb_coefs = ic2*beta*beta;
      I_NAI_F3x_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C3_bb += SQ_NAI_F_S_C3_bb_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_C1003_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_C1003_bb_coefs = ic2_1*beta*beta;
      I_NAI_I6x_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_C1003_bb += SQ_NAI_I_S_C1003_bb_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1003_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1003_bb_coefs = ic2_1*beta*beta;
      I_NAI_H5x_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1003_bb += SQ_NAI_H_S_C1003_bb_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_bb_coefs = ic2_1*beta*beta;
      I_NAI_G4x_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_bb += SQ_NAI_G_S_C1003_bb_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1003_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1003_bb_coefs = ic2_1*beta*beta;
      I_NAI_F3x_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1003_bb += SQ_NAI_F_S_C1003_bb_coefs*I_NAI_F3z_S_vrr;
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
   * shell quartet name: SQ_NAI_P_P_C1003
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   * RHS shell quartet name: SQ_NAI_P_S_C1003
   ************************************************************/
  Double I_NAI_Px_Px_C1003 = I_NAI_D2x_S_C1003+ABX*I_NAI_Px_S_C1003;
  Double I_NAI_Py_Px_C1003 = I_NAI_Dxy_S_C1003+ABX*I_NAI_Py_S_C1003;
  Double I_NAI_Pz_Px_C1003 = I_NAI_Dxz_S_C1003+ABX*I_NAI_Pz_S_C1003;
  Double I_NAI_Px_Py_C1003 = I_NAI_Dxy_S_C1003+ABY*I_NAI_Px_S_C1003;
  Double I_NAI_Py_Py_C1003 = I_NAI_D2y_S_C1003+ABY*I_NAI_Py_S_C1003;
  Double I_NAI_Pz_Py_C1003 = I_NAI_Dyz_S_C1003+ABY*I_NAI_Pz_S_C1003;
  Double I_NAI_Px_Pz_C1003 = I_NAI_Dxz_S_C1003+ABZ*I_NAI_Px_S_C1003;
  Double I_NAI_Py_Pz_C1003 = I_NAI_Dyz_S_C1003+ABZ*I_NAI_Py_S_C1003;
  Double I_NAI_Pz_Pz_C1003 = I_NAI_D2z_S_C1003+ABZ*I_NAI_Pz_S_C1003;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_F_S_C1003_a
   ************************************************************/
  Double I_NAI_F3x_Px_C1003_a = I_NAI_G4x_S_C1003_a+ABX*I_NAI_F3x_S_C1003_a;
  Double I_NAI_F2xy_Px_C1003_a = I_NAI_G3xy_S_C1003_a+ABX*I_NAI_F2xy_S_C1003_a;
  Double I_NAI_F2xz_Px_C1003_a = I_NAI_G3xz_S_C1003_a+ABX*I_NAI_F2xz_S_C1003_a;
  Double I_NAI_Fx2y_Px_C1003_a = I_NAI_G2x2y_S_C1003_a+ABX*I_NAI_Fx2y_S_C1003_a;
  Double I_NAI_Fxyz_Px_C1003_a = I_NAI_G2xyz_S_C1003_a+ABX*I_NAI_Fxyz_S_C1003_a;
  Double I_NAI_Fx2z_Px_C1003_a = I_NAI_G2x2z_S_C1003_a+ABX*I_NAI_Fx2z_S_C1003_a;
  Double I_NAI_F3y_Px_C1003_a = I_NAI_Gx3y_S_C1003_a+ABX*I_NAI_F3y_S_C1003_a;
  Double I_NAI_F2yz_Px_C1003_a = I_NAI_Gx2yz_S_C1003_a+ABX*I_NAI_F2yz_S_C1003_a;
  Double I_NAI_Fy2z_Px_C1003_a = I_NAI_Gxy2z_S_C1003_a+ABX*I_NAI_Fy2z_S_C1003_a;
  Double I_NAI_F3z_Px_C1003_a = I_NAI_Gx3z_S_C1003_a+ABX*I_NAI_F3z_S_C1003_a;
  Double I_NAI_F3x_Py_C1003_a = I_NAI_G3xy_S_C1003_a+ABY*I_NAI_F3x_S_C1003_a;
  Double I_NAI_F2xy_Py_C1003_a = I_NAI_G2x2y_S_C1003_a+ABY*I_NAI_F2xy_S_C1003_a;
  Double I_NAI_F2xz_Py_C1003_a = I_NAI_G2xyz_S_C1003_a+ABY*I_NAI_F2xz_S_C1003_a;
  Double I_NAI_Fx2y_Py_C1003_a = I_NAI_Gx3y_S_C1003_a+ABY*I_NAI_Fx2y_S_C1003_a;
  Double I_NAI_Fxyz_Py_C1003_a = I_NAI_Gx2yz_S_C1003_a+ABY*I_NAI_Fxyz_S_C1003_a;
  Double I_NAI_Fx2z_Py_C1003_a = I_NAI_Gxy2z_S_C1003_a+ABY*I_NAI_Fx2z_S_C1003_a;
  Double I_NAI_F3y_Py_C1003_a = I_NAI_G4y_S_C1003_a+ABY*I_NAI_F3y_S_C1003_a;
  Double I_NAI_F2yz_Py_C1003_a = I_NAI_G3yz_S_C1003_a+ABY*I_NAI_F2yz_S_C1003_a;
  Double I_NAI_Fy2z_Py_C1003_a = I_NAI_G2y2z_S_C1003_a+ABY*I_NAI_Fy2z_S_C1003_a;
  Double I_NAI_F3z_Py_C1003_a = I_NAI_Gy3z_S_C1003_a+ABY*I_NAI_F3z_S_C1003_a;
  Double I_NAI_F3x_Pz_C1003_a = I_NAI_G3xz_S_C1003_a+ABZ*I_NAI_F3x_S_C1003_a;
  Double I_NAI_F2xy_Pz_C1003_a = I_NAI_G2xyz_S_C1003_a+ABZ*I_NAI_F2xy_S_C1003_a;
  Double I_NAI_F2xz_Pz_C1003_a = I_NAI_G2x2z_S_C1003_a+ABZ*I_NAI_F2xz_S_C1003_a;
  Double I_NAI_Fx2y_Pz_C1003_a = I_NAI_Gx2yz_S_C1003_a+ABZ*I_NAI_Fx2y_S_C1003_a;
  Double I_NAI_Fxyz_Pz_C1003_a = I_NAI_Gxy2z_S_C1003_a+ABZ*I_NAI_Fxyz_S_C1003_a;
  Double I_NAI_Fx2z_Pz_C1003_a = I_NAI_Gx3z_S_C1003_a+ABZ*I_NAI_Fx2z_S_C1003_a;
  Double I_NAI_F3y_Pz_C1003_a = I_NAI_G3yz_S_C1003_a+ABZ*I_NAI_F3y_S_C1003_a;
  Double I_NAI_F2yz_Pz_C1003_a = I_NAI_G2y2z_S_C1003_a+ABZ*I_NAI_F2yz_S_C1003_a;
  Double I_NAI_Fy2z_Pz_C1003_a = I_NAI_Gy3z_S_C1003_a+ABZ*I_NAI_Fy2z_S_C1003_a;
  Double I_NAI_F3z_Pz_C1003_a = I_NAI_G4z_S_C1003_a+ABZ*I_NAI_F3z_S_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C3_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   * RHS shell quartet name: SQ_NAI_D_S_C3_b
   ************************************************************/
  Double I_NAI_D2x_Px_C3_b = I_NAI_F3x_S_C3_b+ABX*I_NAI_D2x_S_C3_b;
  Double I_NAI_Dxy_Px_C3_b = I_NAI_F2xy_S_C3_b+ABX*I_NAI_Dxy_S_C3_b;
  Double I_NAI_Dxz_Px_C3_b = I_NAI_F2xz_S_C3_b+ABX*I_NAI_Dxz_S_C3_b;
  Double I_NAI_D2y_Px_C3_b = I_NAI_Fx2y_S_C3_b+ABX*I_NAI_D2y_S_C3_b;
  Double I_NAI_Dyz_Px_C3_b = I_NAI_Fxyz_S_C3_b+ABX*I_NAI_Dyz_S_C3_b;
  Double I_NAI_D2z_Px_C3_b = I_NAI_Fx2z_S_C3_b+ABX*I_NAI_D2z_S_C3_b;
  Double I_NAI_D2x_Py_C3_b = I_NAI_F2xy_S_C3_b+ABY*I_NAI_D2x_S_C3_b;
  Double I_NAI_Dxy_Py_C3_b = I_NAI_Fx2y_S_C3_b+ABY*I_NAI_Dxy_S_C3_b;
  Double I_NAI_Dxz_Py_C3_b = I_NAI_Fxyz_S_C3_b+ABY*I_NAI_Dxz_S_C3_b;
  Double I_NAI_D2y_Py_C3_b = I_NAI_F3y_S_C3_b+ABY*I_NAI_D2y_S_C3_b;
  Double I_NAI_Dyz_Py_C3_b = I_NAI_F2yz_S_C3_b+ABY*I_NAI_Dyz_S_C3_b;
  Double I_NAI_D2z_Py_C3_b = I_NAI_Fy2z_S_C3_b+ABY*I_NAI_D2z_S_C3_b;
  Double I_NAI_D2x_Pz_C3_b = I_NAI_F2xz_S_C3_b+ABZ*I_NAI_D2x_S_C3_b;
  Double I_NAI_Dxy_Pz_C3_b = I_NAI_Fxyz_S_C3_b+ABZ*I_NAI_Dxy_S_C3_b;
  Double I_NAI_Dxz_Pz_C3_b = I_NAI_Fx2z_S_C3_b+ABZ*I_NAI_Dxz_S_C3_b;
  Double I_NAI_D2y_Pz_C3_b = I_NAI_F2yz_S_C3_b+ABZ*I_NAI_D2y_S_C3_b;
  Double I_NAI_Dyz_Pz_C3_b = I_NAI_Fy2z_S_C3_b+ABZ*I_NAI_Dyz_S_C3_b;
  Double I_NAI_D2z_Pz_C3_b = I_NAI_F3z_S_C3_b+ABZ*I_NAI_D2z_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003_b
   ************************************************************/
  Double I_NAI_D2x_Px_C1003_b = I_NAI_F3x_S_C1003_b+ABX*I_NAI_D2x_S_C1003_b;
  Double I_NAI_Dxy_Px_C1003_b = I_NAI_F2xy_S_C1003_b+ABX*I_NAI_Dxy_S_C1003_b;
  Double I_NAI_Dxz_Px_C1003_b = I_NAI_F2xz_S_C1003_b+ABX*I_NAI_Dxz_S_C1003_b;
  Double I_NAI_D2y_Px_C1003_b = I_NAI_Fx2y_S_C1003_b+ABX*I_NAI_D2y_S_C1003_b;
  Double I_NAI_Dyz_Px_C1003_b = I_NAI_Fxyz_S_C1003_b+ABX*I_NAI_Dyz_S_C1003_b;
  Double I_NAI_D2z_Px_C1003_b = I_NAI_Fx2z_S_C1003_b+ABX*I_NAI_D2z_S_C1003_b;
  Double I_NAI_D2x_Py_C1003_b = I_NAI_F2xy_S_C1003_b+ABY*I_NAI_D2x_S_C1003_b;
  Double I_NAI_Dxy_Py_C1003_b = I_NAI_Fx2y_S_C1003_b+ABY*I_NAI_Dxy_S_C1003_b;
  Double I_NAI_Dxz_Py_C1003_b = I_NAI_Fxyz_S_C1003_b+ABY*I_NAI_Dxz_S_C1003_b;
  Double I_NAI_D2y_Py_C1003_b = I_NAI_F3y_S_C1003_b+ABY*I_NAI_D2y_S_C1003_b;
  Double I_NAI_Dyz_Py_C1003_b = I_NAI_F2yz_S_C1003_b+ABY*I_NAI_Dyz_S_C1003_b;
  Double I_NAI_D2z_Py_C1003_b = I_NAI_Fy2z_S_C1003_b+ABY*I_NAI_D2z_S_C1003_b;
  Double I_NAI_D2x_Pz_C1003_b = I_NAI_F2xz_S_C1003_b+ABZ*I_NAI_D2x_S_C1003_b;
  Double I_NAI_Dxy_Pz_C1003_b = I_NAI_Fxyz_S_C1003_b+ABZ*I_NAI_Dxy_S_C1003_b;
  Double I_NAI_Dxz_Pz_C1003_b = I_NAI_Fx2z_S_C1003_b+ABZ*I_NAI_Dxz_S_C1003_b;
  Double I_NAI_D2y_Pz_C1003_b = I_NAI_F2yz_S_C1003_b+ABZ*I_NAI_D2y_S_C1003_b;
  Double I_NAI_Dyz_Pz_C1003_b = I_NAI_Fy2z_S_C1003_b+ABZ*I_NAI_Dyz_S_C1003_b;
  Double I_NAI_D2z_Pz_C1003_b = I_NAI_F3z_S_C1003_b+ABZ*I_NAI_D2z_S_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1003_b
   * RHS shell quartet name: SQ_NAI_F_S_C1003_b
   ************************************************************/
  Double I_NAI_F3x_Px_C1003_b = I_NAI_G4x_S_C1003_b+ABX*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Px_C1003_b = I_NAI_G3xy_S_C1003_b+ABX*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Px_C1003_b = I_NAI_G3xz_S_C1003_b+ABX*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Px_C1003_b = I_NAI_G2x2y_S_C1003_b+ABX*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Px_C1003_b = I_NAI_G2xyz_S_C1003_b+ABX*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Px_C1003_b = I_NAI_G2x2z_S_C1003_b+ABX*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Px_C1003_b = I_NAI_Gx3y_S_C1003_b+ABX*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Px_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABX*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Px_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABX*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Px_C1003_b = I_NAI_Gx3z_S_C1003_b+ABX*I_NAI_F3z_S_C1003_b;
  Double I_NAI_F3x_Py_C1003_b = I_NAI_G3xy_S_C1003_b+ABY*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Py_C1003_b = I_NAI_G2x2y_S_C1003_b+ABY*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Py_C1003_b = I_NAI_G2xyz_S_C1003_b+ABY*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Py_C1003_b = I_NAI_Gx3y_S_C1003_b+ABY*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Py_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABY*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Py_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABY*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Py_C1003_b = I_NAI_G4y_S_C1003_b+ABY*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Py_C1003_b = I_NAI_G3yz_S_C1003_b+ABY*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Py_C1003_b = I_NAI_G2y2z_S_C1003_b+ABY*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Py_C1003_b = I_NAI_Gy3z_S_C1003_b+ABY*I_NAI_F3z_S_C1003_b;
  Double I_NAI_F3x_Pz_C1003_b = I_NAI_G3xz_S_C1003_b+ABZ*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Pz_C1003_b = I_NAI_G2xyz_S_C1003_b+ABZ*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Pz_C1003_b = I_NAI_G2x2z_S_C1003_b+ABZ*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Pz_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABZ*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Pz_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABZ*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Pz_C1003_b = I_NAI_Gx3z_S_C1003_b+ABZ*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Pz_C1003_b = I_NAI_G3yz_S_C1003_b+ABZ*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Pz_C1003_b = I_NAI_G2y2z_S_C1003_b+ABZ*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Pz_C1003_b = I_NAI_Gy3z_S_C1003_b+ABZ*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Pz_C1003_b = I_NAI_G4z_S_C1003_b+ABZ*I_NAI_F3z_S_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_D_P_C1003_b
   ************************************************************/
  Double I_NAI_D2x_D2x_C1003_b = I_NAI_F3x_Px_C1003_b+ABX*I_NAI_D2x_Px_C1003_b;
  Double I_NAI_Dxy_D2x_C1003_b = I_NAI_F2xy_Px_C1003_b+ABX*I_NAI_Dxy_Px_C1003_b;
  Double I_NAI_Dxz_D2x_C1003_b = I_NAI_F2xz_Px_C1003_b+ABX*I_NAI_Dxz_Px_C1003_b;
  Double I_NAI_D2y_D2x_C1003_b = I_NAI_Fx2y_Px_C1003_b+ABX*I_NAI_D2y_Px_C1003_b;
  Double I_NAI_Dyz_D2x_C1003_b = I_NAI_Fxyz_Px_C1003_b+ABX*I_NAI_Dyz_Px_C1003_b;
  Double I_NAI_D2z_D2x_C1003_b = I_NAI_Fx2z_Px_C1003_b+ABX*I_NAI_D2z_Px_C1003_b;
  Double I_NAI_D2x_Dxy_C1003_b = I_NAI_F2xy_Px_C1003_b+ABY*I_NAI_D2x_Px_C1003_b;
  Double I_NAI_Dxy_Dxy_C1003_b = I_NAI_Fx2y_Px_C1003_b+ABY*I_NAI_Dxy_Px_C1003_b;
  Double I_NAI_Dxz_Dxy_C1003_b = I_NAI_Fxyz_Px_C1003_b+ABY*I_NAI_Dxz_Px_C1003_b;
  Double I_NAI_D2y_Dxy_C1003_b = I_NAI_F3y_Px_C1003_b+ABY*I_NAI_D2y_Px_C1003_b;
  Double I_NAI_Dyz_Dxy_C1003_b = I_NAI_F2yz_Px_C1003_b+ABY*I_NAI_Dyz_Px_C1003_b;
  Double I_NAI_D2z_Dxy_C1003_b = I_NAI_Fy2z_Px_C1003_b+ABY*I_NAI_D2z_Px_C1003_b;
  Double I_NAI_D2x_Dxz_C1003_b = I_NAI_F2xz_Px_C1003_b+ABZ*I_NAI_D2x_Px_C1003_b;
  Double I_NAI_Dxy_Dxz_C1003_b = I_NAI_Fxyz_Px_C1003_b+ABZ*I_NAI_Dxy_Px_C1003_b;
  Double I_NAI_Dxz_Dxz_C1003_b = I_NAI_Fx2z_Px_C1003_b+ABZ*I_NAI_Dxz_Px_C1003_b;
  Double I_NAI_D2y_Dxz_C1003_b = I_NAI_F2yz_Px_C1003_b+ABZ*I_NAI_D2y_Px_C1003_b;
  Double I_NAI_Dyz_Dxz_C1003_b = I_NAI_Fy2z_Px_C1003_b+ABZ*I_NAI_Dyz_Px_C1003_b;
  Double I_NAI_D2z_Dxz_C1003_b = I_NAI_F3z_Px_C1003_b+ABZ*I_NAI_D2z_Px_C1003_b;
  Double I_NAI_D2x_D2y_C1003_b = I_NAI_F2xy_Py_C1003_b+ABY*I_NAI_D2x_Py_C1003_b;
  Double I_NAI_Dxy_D2y_C1003_b = I_NAI_Fx2y_Py_C1003_b+ABY*I_NAI_Dxy_Py_C1003_b;
  Double I_NAI_Dxz_D2y_C1003_b = I_NAI_Fxyz_Py_C1003_b+ABY*I_NAI_Dxz_Py_C1003_b;
  Double I_NAI_D2y_D2y_C1003_b = I_NAI_F3y_Py_C1003_b+ABY*I_NAI_D2y_Py_C1003_b;
  Double I_NAI_Dyz_D2y_C1003_b = I_NAI_F2yz_Py_C1003_b+ABY*I_NAI_Dyz_Py_C1003_b;
  Double I_NAI_D2z_D2y_C1003_b = I_NAI_Fy2z_Py_C1003_b+ABY*I_NAI_D2z_Py_C1003_b;
  Double I_NAI_D2x_Dyz_C1003_b = I_NAI_F2xz_Py_C1003_b+ABZ*I_NAI_D2x_Py_C1003_b;
  Double I_NAI_Dxy_Dyz_C1003_b = I_NAI_Fxyz_Py_C1003_b+ABZ*I_NAI_Dxy_Py_C1003_b;
  Double I_NAI_Dxz_Dyz_C1003_b = I_NAI_Fx2z_Py_C1003_b+ABZ*I_NAI_Dxz_Py_C1003_b;
  Double I_NAI_D2y_Dyz_C1003_b = I_NAI_F2yz_Py_C1003_b+ABZ*I_NAI_D2y_Py_C1003_b;
  Double I_NAI_Dyz_Dyz_C1003_b = I_NAI_Fy2z_Py_C1003_b+ABZ*I_NAI_Dyz_Py_C1003_b;
  Double I_NAI_D2z_Dyz_C1003_b = I_NAI_F3z_Py_C1003_b+ABZ*I_NAI_D2z_Py_C1003_b;
  Double I_NAI_D2x_D2z_C1003_b = I_NAI_F2xz_Pz_C1003_b+ABZ*I_NAI_D2x_Pz_C1003_b;
  Double I_NAI_Dxy_D2z_C1003_b = I_NAI_Fxyz_Pz_C1003_b+ABZ*I_NAI_Dxy_Pz_C1003_b;
  Double I_NAI_Dxz_D2z_C1003_b = I_NAI_Fx2z_Pz_C1003_b+ABZ*I_NAI_Dxz_Pz_C1003_b;
  Double I_NAI_D2y_D2z_C1003_b = I_NAI_F2yz_Pz_C1003_b+ABZ*I_NAI_D2y_Pz_C1003_b;
  Double I_NAI_Dyz_D2z_C1003_b = I_NAI_Fy2z_Pz_C1003_b+ABZ*I_NAI_Dyz_Pz_C1003_b;
  Double I_NAI_D2z_D2z_C1003_b = I_NAI_F3z_Pz_C1003_b+ABZ*I_NAI_D2z_Pz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1003_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C1003_aa
   * RHS shell quartet name: SQ_NAI_H_S_C1003_aa
   ************************************************************/
  Double I_NAI_H5x_Px_C1003_aa = I_NAI_I6x_S_C1003_aa+ABX*I_NAI_H5x_S_C1003_aa;
  Double I_NAI_H4xy_Px_C1003_aa = I_NAI_I5xy_S_C1003_aa+ABX*I_NAI_H4xy_S_C1003_aa;
  Double I_NAI_H4xz_Px_C1003_aa = I_NAI_I5xz_S_C1003_aa+ABX*I_NAI_H4xz_S_C1003_aa;
  Double I_NAI_H3x2y_Px_C1003_aa = I_NAI_I4x2y_S_C1003_aa+ABX*I_NAI_H3x2y_S_C1003_aa;
  Double I_NAI_H3xyz_Px_C1003_aa = I_NAI_I4xyz_S_C1003_aa+ABX*I_NAI_H3xyz_S_C1003_aa;
  Double I_NAI_H3x2z_Px_C1003_aa = I_NAI_I4x2z_S_C1003_aa+ABX*I_NAI_H3x2z_S_C1003_aa;
  Double I_NAI_H2x3y_Px_C1003_aa = I_NAI_I3x3y_S_C1003_aa+ABX*I_NAI_H2x3y_S_C1003_aa;
  Double I_NAI_H2x2yz_Px_C1003_aa = I_NAI_I3x2yz_S_C1003_aa+ABX*I_NAI_H2x2yz_S_C1003_aa;
  Double I_NAI_H2xy2z_Px_C1003_aa = I_NAI_I3xy2z_S_C1003_aa+ABX*I_NAI_H2xy2z_S_C1003_aa;
  Double I_NAI_H2x3z_Px_C1003_aa = I_NAI_I3x3z_S_C1003_aa+ABX*I_NAI_H2x3z_S_C1003_aa;
  Double I_NAI_Hx4y_Px_C1003_aa = I_NAI_I2x4y_S_C1003_aa+ABX*I_NAI_Hx4y_S_C1003_aa;
  Double I_NAI_Hx3yz_Px_C1003_aa = I_NAI_I2x3yz_S_C1003_aa+ABX*I_NAI_Hx3yz_S_C1003_aa;
  Double I_NAI_Hx2y2z_Px_C1003_aa = I_NAI_I2x2y2z_S_C1003_aa+ABX*I_NAI_Hx2y2z_S_C1003_aa;
  Double I_NAI_Hxy3z_Px_C1003_aa = I_NAI_I2xy3z_S_C1003_aa+ABX*I_NAI_Hxy3z_S_C1003_aa;
  Double I_NAI_Hx4z_Px_C1003_aa = I_NAI_I2x4z_S_C1003_aa+ABX*I_NAI_Hx4z_S_C1003_aa;
  Double I_NAI_H5y_Px_C1003_aa = I_NAI_Ix5y_S_C1003_aa+ABX*I_NAI_H5y_S_C1003_aa;
  Double I_NAI_H4yz_Px_C1003_aa = I_NAI_Ix4yz_S_C1003_aa+ABX*I_NAI_H4yz_S_C1003_aa;
  Double I_NAI_H3y2z_Px_C1003_aa = I_NAI_Ix3y2z_S_C1003_aa+ABX*I_NAI_H3y2z_S_C1003_aa;
  Double I_NAI_H2y3z_Px_C1003_aa = I_NAI_Ix2y3z_S_C1003_aa+ABX*I_NAI_H2y3z_S_C1003_aa;
  Double I_NAI_Hy4z_Px_C1003_aa = I_NAI_Ixy4z_S_C1003_aa+ABX*I_NAI_Hy4z_S_C1003_aa;
  Double I_NAI_H5z_Px_C1003_aa = I_NAI_Ix5z_S_C1003_aa+ABX*I_NAI_H5z_S_C1003_aa;
  Double I_NAI_H5x_Py_C1003_aa = I_NAI_I5xy_S_C1003_aa+ABY*I_NAI_H5x_S_C1003_aa;
  Double I_NAI_H4xy_Py_C1003_aa = I_NAI_I4x2y_S_C1003_aa+ABY*I_NAI_H4xy_S_C1003_aa;
  Double I_NAI_H4xz_Py_C1003_aa = I_NAI_I4xyz_S_C1003_aa+ABY*I_NAI_H4xz_S_C1003_aa;
  Double I_NAI_H3x2y_Py_C1003_aa = I_NAI_I3x3y_S_C1003_aa+ABY*I_NAI_H3x2y_S_C1003_aa;
  Double I_NAI_H3xyz_Py_C1003_aa = I_NAI_I3x2yz_S_C1003_aa+ABY*I_NAI_H3xyz_S_C1003_aa;
  Double I_NAI_H3x2z_Py_C1003_aa = I_NAI_I3xy2z_S_C1003_aa+ABY*I_NAI_H3x2z_S_C1003_aa;
  Double I_NAI_H2x3y_Py_C1003_aa = I_NAI_I2x4y_S_C1003_aa+ABY*I_NAI_H2x3y_S_C1003_aa;
  Double I_NAI_H2x2yz_Py_C1003_aa = I_NAI_I2x3yz_S_C1003_aa+ABY*I_NAI_H2x2yz_S_C1003_aa;
  Double I_NAI_H2xy2z_Py_C1003_aa = I_NAI_I2x2y2z_S_C1003_aa+ABY*I_NAI_H2xy2z_S_C1003_aa;
  Double I_NAI_H2x3z_Py_C1003_aa = I_NAI_I2xy3z_S_C1003_aa+ABY*I_NAI_H2x3z_S_C1003_aa;
  Double I_NAI_Hx4y_Py_C1003_aa = I_NAI_Ix5y_S_C1003_aa+ABY*I_NAI_Hx4y_S_C1003_aa;
  Double I_NAI_Hx3yz_Py_C1003_aa = I_NAI_Ix4yz_S_C1003_aa+ABY*I_NAI_Hx3yz_S_C1003_aa;
  Double I_NAI_Hx2y2z_Py_C1003_aa = I_NAI_Ix3y2z_S_C1003_aa+ABY*I_NAI_Hx2y2z_S_C1003_aa;
  Double I_NAI_Hxy3z_Py_C1003_aa = I_NAI_Ix2y3z_S_C1003_aa+ABY*I_NAI_Hxy3z_S_C1003_aa;
  Double I_NAI_Hx4z_Py_C1003_aa = I_NAI_Ixy4z_S_C1003_aa+ABY*I_NAI_Hx4z_S_C1003_aa;
  Double I_NAI_H5y_Py_C1003_aa = I_NAI_I6y_S_C1003_aa+ABY*I_NAI_H5y_S_C1003_aa;
  Double I_NAI_H4yz_Py_C1003_aa = I_NAI_I5yz_S_C1003_aa+ABY*I_NAI_H4yz_S_C1003_aa;
  Double I_NAI_H3y2z_Py_C1003_aa = I_NAI_I4y2z_S_C1003_aa+ABY*I_NAI_H3y2z_S_C1003_aa;
  Double I_NAI_H2y3z_Py_C1003_aa = I_NAI_I3y3z_S_C1003_aa+ABY*I_NAI_H2y3z_S_C1003_aa;
  Double I_NAI_Hy4z_Py_C1003_aa = I_NAI_I2y4z_S_C1003_aa+ABY*I_NAI_Hy4z_S_C1003_aa;
  Double I_NAI_H5z_Py_C1003_aa = I_NAI_Iy5z_S_C1003_aa+ABY*I_NAI_H5z_S_C1003_aa;
  Double I_NAI_H5x_Pz_C1003_aa = I_NAI_I5xz_S_C1003_aa+ABZ*I_NAI_H5x_S_C1003_aa;
  Double I_NAI_H4xy_Pz_C1003_aa = I_NAI_I4xyz_S_C1003_aa+ABZ*I_NAI_H4xy_S_C1003_aa;
  Double I_NAI_H4xz_Pz_C1003_aa = I_NAI_I4x2z_S_C1003_aa+ABZ*I_NAI_H4xz_S_C1003_aa;
  Double I_NAI_H3x2y_Pz_C1003_aa = I_NAI_I3x2yz_S_C1003_aa+ABZ*I_NAI_H3x2y_S_C1003_aa;
  Double I_NAI_H3xyz_Pz_C1003_aa = I_NAI_I3xy2z_S_C1003_aa+ABZ*I_NAI_H3xyz_S_C1003_aa;
  Double I_NAI_H3x2z_Pz_C1003_aa = I_NAI_I3x3z_S_C1003_aa+ABZ*I_NAI_H3x2z_S_C1003_aa;
  Double I_NAI_H2x3y_Pz_C1003_aa = I_NAI_I2x3yz_S_C1003_aa+ABZ*I_NAI_H2x3y_S_C1003_aa;
  Double I_NAI_H2x2yz_Pz_C1003_aa = I_NAI_I2x2y2z_S_C1003_aa+ABZ*I_NAI_H2x2yz_S_C1003_aa;
  Double I_NAI_H2xy2z_Pz_C1003_aa = I_NAI_I2xy3z_S_C1003_aa+ABZ*I_NAI_H2xy2z_S_C1003_aa;
  Double I_NAI_H2x3z_Pz_C1003_aa = I_NAI_I2x4z_S_C1003_aa+ABZ*I_NAI_H2x3z_S_C1003_aa;
  Double I_NAI_Hx4y_Pz_C1003_aa = I_NAI_Ix4yz_S_C1003_aa+ABZ*I_NAI_Hx4y_S_C1003_aa;
  Double I_NAI_Hx3yz_Pz_C1003_aa = I_NAI_Ix3y2z_S_C1003_aa+ABZ*I_NAI_Hx3yz_S_C1003_aa;
  Double I_NAI_Hx2y2z_Pz_C1003_aa = I_NAI_Ix2y3z_S_C1003_aa+ABZ*I_NAI_Hx2y2z_S_C1003_aa;
  Double I_NAI_Hxy3z_Pz_C1003_aa = I_NAI_Ixy4z_S_C1003_aa+ABZ*I_NAI_Hxy3z_S_C1003_aa;
  Double I_NAI_Hx4z_Pz_C1003_aa = I_NAI_Ix5z_S_C1003_aa+ABZ*I_NAI_Hx4z_S_C1003_aa;
  Double I_NAI_H5y_Pz_C1003_aa = I_NAI_I5yz_S_C1003_aa+ABZ*I_NAI_H5y_S_C1003_aa;
  Double I_NAI_H4yz_Pz_C1003_aa = I_NAI_I4y2z_S_C1003_aa+ABZ*I_NAI_H4yz_S_C1003_aa;
  Double I_NAI_H3y2z_Pz_C1003_aa = I_NAI_I3y3z_S_C1003_aa+ABZ*I_NAI_H3y2z_S_C1003_aa;
  Double I_NAI_H2y3z_Pz_C1003_aa = I_NAI_I2y4z_S_C1003_aa+ABZ*I_NAI_H2y3z_S_C1003_aa;
  Double I_NAI_Hy4z_Pz_C1003_aa = I_NAI_Iy5z_S_C1003_aa+ABZ*I_NAI_Hy4z_S_C1003_aa;
  Double I_NAI_H5z_Pz_C1003_aa = I_NAI_I6z_S_C1003_aa+ABZ*I_NAI_H5z_S_C1003_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C3_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_ab
   * RHS shell quartet name: SQ_NAI_G_S_C3_ab
   ************************************************************/
  Double I_NAI_G4x_Px_C3_ab = I_NAI_H5x_S_C3_ab+ABX*I_NAI_G4x_S_C3_ab;
  Double I_NAI_G3xy_Px_C3_ab = I_NAI_H4xy_S_C3_ab+ABX*I_NAI_G3xy_S_C3_ab;
  Double I_NAI_G3xz_Px_C3_ab = I_NAI_H4xz_S_C3_ab+ABX*I_NAI_G3xz_S_C3_ab;
  Double I_NAI_G2x2y_Px_C3_ab = I_NAI_H3x2y_S_C3_ab+ABX*I_NAI_G2x2y_S_C3_ab;
  Double I_NAI_G2xyz_Px_C3_ab = I_NAI_H3xyz_S_C3_ab+ABX*I_NAI_G2xyz_S_C3_ab;
  Double I_NAI_G2x2z_Px_C3_ab = I_NAI_H3x2z_S_C3_ab+ABX*I_NAI_G2x2z_S_C3_ab;
  Double I_NAI_Gx3y_Px_C3_ab = I_NAI_H2x3y_S_C3_ab+ABX*I_NAI_Gx3y_S_C3_ab;
  Double I_NAI_Gx2yz_Px_C3_ab = I_NAI_H2x2yz_S_C3_ab+ABX*I_NAI_Gx2yz_S_C3_ab;
  Double I_NAI_Gxy2z_Px_C3_ab = I_NAI_H2xy2z_S_C3_ab+ABX*I_NAI_Gxy2z_S_C3_ab;
  Double I_NAI_Gx3z_Px_C3_ab = I_NAI_H2x3z_S_C3_ab+ABX*I_NAI_Gx3z_S_C3_ab;
  Double I_NAI_G4y_Px_C3_ab = I_NAI_Hx4y_S_C3_ab+ABX*I_NAI_G4y_S_C3_ab;
  Double I_NAI_G3yz_Px_C3_ab = I_NAI_Hx3yz_S_C3_ab+ABX*I_NAI_G3yz_S_C3_ab;
  Double I_NAI_G2y2z_Px_C3_ab = I_NAI_Hx2y2z_S_C3_ab+ABX*I_NAI_G2y2z_S_C3_ab;
  Double I_NAI_Gy3z_Px_C3_ab = I_NAI_Hxy3z_S_C3_ab+ABX*I_NAI_Gy3z_S_C3_ab;
  Double I_NAI_G4z_Px_C3_ab = I_NAI_Hx4z_S_C3_ab+ABX*I_NAI_G4z_S_C3_ab;
  Double I_NAI_G4x_Py_C3_ab = I_NAI_H4xy_S_C3_ab+ABY*I_NAI_G4x_S_C3_ab;
  Double I_NAI_G3xy_Py_C3_ab = I_NAI_H3x2y_S_C3_ab+ABY*I_NAI_G3xy_S_C3_ab;
  Double I_NAI_G3xz_Py_C3_ab = I_NAI_H3xyz_S_C3_ab+ABY*I_NAI_G3xz_S_C3_ab;
  Double I_NAI_G2x2y_Py_C3_ab = I_NAI_H2x3y_S_C3_ab+ABY*I_NAI_G2x2y_S_C3_ab;
  Double I_NAI_G2xyz_Py_C3_ab = I_NAI_H2x2yz_S_C3_ab+ABY*I_NAI_G2xyz_S_C3_ab;
  Double I_NAI_G2x2z_Py_C3_ab = I_NAI_H2xy2z_S_C3_ab+ABY*I_NAI_G2x2z_S_C3_ab;
  Double I_NAI_Gx3y_Py_C3_ab = I_NAI_Hx4y_S_C3_ab+ABY*I_NAI_Gx3y_S_C3_ab;
  Double I_NAI_Gx2yz_Py_C3_ab = I_NAI_Hx3yz_S_C3_ab+ABY*I_NAI_Gx2yz_S_C3_ab;
  Double I_NAI_Gxy2z_Py_C3_ab = I_NAI_Hx2y2z_S_C3_ab+ABY*I_NAI_Gxy2z_S_C3_ab;
  Double I_NAI_Gx3z_Py_C3_ab = I_NAI_Hxy3z_S_C3_ab+ABY*I_NAI_Gx3z_S_C3_ab;
  Double I_NAI_G4y_Py_C3_ab = I_NAI_H5y_S_C3_ab+ABY*I_NAI_G4y_S_C3_ab;
  Double I_NAI_G3yz_Py_C3_ab = I_NAI_H4yz_S_C3_ab+ABY*I_NAI_G3yz_S_C3_ab;
  Double I_NAI_G2y2z_Py_C3_ab = I_NAI_H3y2z_S_C3_ab+ABY*I_NAI_G2y2z_S_C3_ab;
  Double I_NAI_Gy3z_Py_C3_ab = I_NAI_H2y3z_S_C3_ab+ABY*I_NAI_Gy3z_S_C3_ab;
  Double I_NAI_G4z_Py_C3_ab = I_NAI_Hy4z_S_C3_ab+ABY*I_NAI_G4z_S_C3_ab;
  Double I_NAI_G4x_Pz_C3_ab = I_NAI_H4xz_S_C3_ab+ABZ*I_NAI_G4x_S_C3_ab;
  Double I_NAI_G3xy_Pz_C3_ab = I_NAI_H3xyz_S_C3_ab+ABZ*I_NAI_G3xy_S_C3_ab;
  Double I_NAI_G3xz_Pz_C3_ab = I_NAI_H3x2z_S_C3_ab+ABZ*I_NAI_G3xz_S_C3_ab;
  Double I_NAI_G2x2y_Pz_C3_ab = I_NAI_H2x2yz_S_C3_ab+ABZ*I_NAI_G2x2y_S_C3_ab;
  Double I_NAI_G2xyz_Pz_C3_ab = I_NAI_H2xy2z_S_C3_ab+ABZ*I_NAI_G2xyz_S_C3_ab;
  Double I_NAI_G2x2z_Pz_C3_ab = I_NAI_H2x3z_S_C3_ab+ABZ*I_NAI_G2x2z_S_C3_ab;
  Double I_NAI_Gx3y_Pz_C3_ab = I_NAI_Hx3yz_S_C3_ab+ABZ*I_NAI_Gx3y_S_C3_ab;
  Double I_NAI_Gx2yz_Pz_C3_ab = I_NAI_Hx2y2z_S_C3_ab+ABZ*I_NAI_Gx2yz_S_C3_ab;
  Double I_NAI_Gxy2z_Pz_C3_ab = I_NAI_Hxy3z_S_C3_ab+ABZ*I_NAI_Gxy2z_S_C3_ab;
  Double I_NAI_Gx3z_Pz_C3_ab = I_NAI_Hx4z_S_C3_ab+ABZ*I_NAI_Gx3z_S_C3_ab;
  Double I_NAI_G4y_Pz_C3_ab = I_NAI_H4yz_S_C3_ab+ABZ*I_NAI_G4y_S_C3_ab;
  Double I_NAI_G3yz_Pz_C3_ab = I_NAI_H3y2z_S_C3_ab+ABZ*I_NAI_G3yz_S_C3_ab;
  Double I_NAI_G2y2z_Pz_C3_ab = I_NAI_H2y3z_S_C3_ab+ABZ*I_NAI_G2y2z_S_C3_ab;
  Double I_NAI_Gy3z_Pz_C3_ab = I_NAI_Hy4z_S_C3_ab+ABZ*I_NAI_Gy3z_S_C3_ab;
  Double I_NAI_G4z_Pz_C3_ab = I_NAI_H5z_S_C3_ab+ABZ*I_NAI_G4z_S_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C1003_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_ab
   ************************************************************/
  Double I_NAI_G4x_Px_C1003_ab = I_NAI_H5x_S_C1003_ab+ABX*I_NAI_G4x_S_C1003_ab;
  Double I_NAI_G3xy_Px_C1003_ab = I_NAI_H4xy_S_C1003_ab+ABX*I_NAI_G3xy_S_C1003_ab;
  Double I_NAI_G3xz_Px_C1003_ab = I_NAI_H4xz_S_C1003_ab+ABX*I_NAI_G3xz_S_C1003_ab;
  Double I_NAI_G2x2y_Px_C1003_ab = I_NAI_H3x2y_S_C1003_ab+ABX*I_NAI_G2x2y_S_C1003_ab;
  Double I_NAI_G2xyz_Px_C1003_ab = I_NAI_H3xyz_S_C1003_ab+ABX*I_NAI_G2xyz_S_C1003_ab;
  Double I_NAI_G2x2z_Px_C1003_ab = I_NAI_H3x2z_S_C1003_ab+ABX*I_NAI_G2x2z_S_C1003_ab;
  Double I_NAI_Gx3y_Px_C1003_ab = I_NAI_H2x3y_S_C1003_ab+ABX*I_NAI_Gx3y_S_C1003_ab;
  Double I_NAI_Gx2yz_Px_C1003_ab = I_NAI_H2x2yz_S_C1003_ab+ABX*I_NAI_Gx2yz_S_C1003_ab;
  Double I_NAI_Gxy2z_Px_C1003_ab = I_NAI_H2xy2z_S_C1003_ab+ABX*I_NAI_Gxy2z_S_C1003_ab;
  Double I_NAI_Gx3z_Px_C1003_ab = I_NAI_H2x3z_S_C1003_ab+ABX*I_NAI_Gx3z_S_C1003_ab;
  Double I_NAI_G4y_Px_C1003_ab = I_NAI_Hx4y_S_C1003_ab+ABX*I_NAI_G4y_S_C1003_ab;
  Double I_NAI_G3yz_Px_C1003_ab = I_NAI_Hx3yz_S_C1003_ab+ABX*I_NAI_G3yz_S_C1003_ab;
  Double I_NAI_G2y2z_Px_C1003_ab = I_NAI_Hx2y2z_S_C1003_ab+ABX*I_NAI_G2y2z_S_C1003_ab;
  Double I_NAI_Gy3z_Px_C1003_ab = I_NAI_Hxy3z_S_C1003_ab+ABX*I_NAI_Gy3z_S_C1003_ab;
  Double I_NAI_G4z_Px_C1003_ab = I_NAI_Hx4z_S_C1003_ab+ABX*I_NAI_G4z_S_C1003_ab;
  Double I_NAI_G4x_Py_C1003_ab = I_NAI_H4xy_S_C1003_ab+ABY*I_NAI_G4x_S_C1003_ab;
  Double I_NAI_G3xy_Py_C1003_ab = I_NAI_H3x2y_S_C1003_ab+ABY*I_NAI_G3xy_S_C1003_ab;
  Double I_NAI_G3xz_Py_C1003_ab = I_NAI_H3xyz_S_C1003_ab+ABY*I_NAI_G3xz_S_C1003_ab;
  Double I_NAI_G2x2y_Py_C1003_ab = I_NAI_H2x3y_S_C1003_ab+ABY*I_NAI_G2x2y_S_C1003_ab;
  Double I_NAI_G2xyz_Py_C1003_ab = I_NAI_H2x2yz_S_C1003_ab+ABY*I_NAI_G2xyz_S_C1003_ab;
  Double I_NAI_G2x2z_Py_C1003_ab = I_NAI_H2xy2z_S_C1003_ab+ABY*I_NAI_G2x2z_S_C1003_ab;
  Double I_NAI_Gx3y_Py_C1003_ab = I_NAI_Hx4y_S_C1003_ab+ABY*I_NAI_Gx3y_S_C1003_ab;
  Double I_NAI_Gx2yz_Py_C1003_ab = I_NAI_Hx3yz_S_C1003_ab+ABY*I_NAI_Gx2yz_S_C1003_ab;
  Double I_NAI_Gxy2z_Py_C1003_ab = I_NAI_Hx2y2z_S_C1003_ab+ABY*I_NAI_Gxy2z_S_C1003_ab;
  Double I_NAI_Gx3z_Py_C1003_ab = I_NAI_Hxy3z_S_C1003_ab+ABY*I_NAI_Gx3z_S_C1003_ab;
  Double I_NAI_G4y_Py_C1003_ab = I_NAI_H5y_S_C1003_ab+ABY*I_NAI_G4y_S_C1003_ab;
  Double I_NAI_G3yz_Py_C1003_ab = I_NAI_H4yz_S_C1003_ab+ABY*I_NAI_G3yz_S_C1003_ab;
  Double I_NAI_G2y2z_Py_C1003_ab = I_NAI_H3y2z_S_C1003_ab+ABY*I_NAI_G2y2z_S_C1003_ab;
  Double I_NAI_Gy3z_Py_C1003_ab = I_NAI_H2y3z_S_C1003_ab+ABY*I_NAI_Gy3z_S_C1003_ab;
  Double I_NAI_G4z_Py_C1003_ab = I_NAI_Hy4z_S_C1003_ab+ABY*I_NAI_G4z_S_C1003_ab;
  Double I_NAI_G4x_Pz_C1003_ab = I_NAI_H4xz_S_C1003_ab+ABZ*I_NAI_G4x_S_C1003_ab;
  Double I_NAI_G3xy_Pz_C1003_ab = I_NAI_H3xyz_S_C1003_ab+ABZ*I_NAI_G3xy_S_C1003_ab;
  Double I_NAI_G3xz_Pz_C1003_ab = I_NAI_H3x2z_S_C1003_ab+ABZ*I_NAI_G3xz_S_C1003_ab;
  Double I_NAI_G2x2y_Pz_C1003_ab = I_NAI_H2x2yz_S_C1003_ab+ABZ*I_NAI_G2x2y_S_C1003_ab;
  Double I_NAI_G2xyz_Pz_C1003_ab = I_NAI_H2xy2z_S_C1003_ab+ABZ*I_NAI_G2xyz_S_C1003_ab;
  Double I_NAI_G2x2z_Pz_C1003_ab = I_NAI_H2x3z_S_C1003_ab+ABZ*I_NAI_G2x2z_S_C1003_ab;
  Double I_NAI_Gx3y_Pz_C1003_ab = I_NAI_Hx3yz_S_C1003_ab+ABZ*I_NAI_Gx3y_S_C1003_ab;
  Double I_NAI_Gx2yz_Pz_C1003_ab = I_NAI_Hx2y2z_S_C1003_ab+ABZ*I_NAI_Gx2yz_S_C1003_ab;
  Double I_NAI_Gxy2z_Pz_C1003_ab = I_NAI_Hxy3z_S_C1003_ab+ABZ*I_NAI_Gxy2z_S_C1003_ab;
  Double I_NAI_Gx3z_Pz_C1003_ab = I_NAI_Hx4z_S_C1003_ab+ABZ*I_NAI_Gx3z_S_C1003_ab;
  Double I_NAI_G4y_Pz_C1003_ab = I_NAI_H4yz_S_C1003_ab+ABZ*I_NAI_G4y_S_C1003_ab;
  Double I_NAI_G3yz_Pz_C1003_ab = I_NAI_H3y2z_S_C1003_ab+ABZ*I_NAI_G3yz_S_C1003_ab;
  Double I_NAI_G2y2z_Pz_C1003_ab = I_NAI_H2y3z_S_C1003_ab+ABZ*I_NAI_G2y2z_S_C1003_ab;
  Double I_NAI_Gy3z_Pz_C1003_ab = I_NAI_Hy4z_S_C1003_ab+ABZ*I_NAI_Gy3z_S_C1003_ab;
  Double I_NAI_G4z_Pz_C1003_ab = I_NAI_H5z_S_C1003_ab+ABZ*I_NAI_G4z_S_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1003_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 7 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C1003_ab
   * RHS shell quartet name: SQ_NAI_H_S_C1003_ab
   ************************************************************/
  Double I_NAI_H5x_Px_C1003_ab = I_NAI_I6x_S_C1003_ab+ABX*I_NAI_H5x_S_C1003_ab;
  Double I_NAI_H4xy_Px_C1003_ab = I_NAI_I5xy_S_C1003_ab+ABX*I_NAI_H4xy_S_C1003_ab;
  Double I_NAI_H4xz_Px_C1003_ab = I_NAI_I5xz_S_C1003_ab+ABX*I_NAI_H4xz_S_C1003_ab;
  Double I_NAI_H3x2y_Px_C1003_ab = I_NAI_I4x2y_S_C1003_ab+ABX*I_NAI_H3x2y_S_C1003_ab;
  Double I_NAI_H3xyz_Px_C1003_ab = I_NAI_I4xyz_S_C1003_ab+ABX*I_NAI_H3xyz_S_C1003_ab;
  Double I_NAI_H3x2z_Px_C1003_ab = I_NAI_I4x2z_S_C1003_ab+ABX*I_NAI_H3x2z_S_C1003_ab;
  Double I_NAI_H2x3y_Px_C1003_ab = I_NAI_I3x3y_S_C1003_ab+ABX*I_NAI_H2x3y_S_C1003_ab;
  Double I_NAI_H2x2yz_Px_C1003_ab = I_NAI_I3x2yz_S_C1003_ab+ABX*I_NAI_H2x2yz_S_C1003_ab;
  Double I_NAI_H2xy2z_Px_C1003_ab = I_NAI_I3xy2z_S_C1003_ab+ABX*I_NAI_H2xy2z_S_C1003_ab;
  Double I_NAI_H2x3z_Px_C1003_ab = I_NAI_I3x3z_S_C1003_ab+ABX*I_NAI_H2x3z_S_C1003_ab;
  Double I_NAI_Hx4y_Px_C1003_ab = I_NAI_I2x4y_S_C1003_ab+ABX*I_NAI_Hx4y_S_C1003_ab;
  Double I_NAI_Hx3yz_Px_C1003_ab = I_NAI_I2x3yz_S_C1003_ab+ABX*I_NAI_Hx3yz_S_C1003_ab;
  Double I_NAI_Hx2y2z_Px_C1003_ab = I_NAI_I2x2y2z_S_C1003_ab+ABX*I_NAI_Hx2y2z_S_C1003_ab;
  Double I_NAI_Hxy3z_Px_C1003_ab = I_NAI_I2xy3z_S_C1003_ab+ABX*I_NAI_Hxy3z_S_C1003_ab;
  Double I_NAI_Hx4z_Px_C1003_ab = I_NAI_I2x4z_S_C1003_ab+ABX*I_NAI_Hx4z_S_C1003_ab;
  Double I_NAI_H5y_Px_C1003_ab = I_NAI_Ix5y_S_C1003_ab+ABX*I_NAI_H5y_S_C1003_ab;
  Double I_NAI_H4yz_Px_C1003_ab = I_NAI_Ix4yz_S_C1003_ab+ABX*I_NAI_H4yz_S_C1003_ab;
  Double I_NAI_H3y2z_Px_C1003_ab = I_NAI_Ix3y2z_S_C1003_ab+ABX*I_NAI_H3y2z_S_C1003_ab;
  Double I_NAI_H2y3z_Px_C1003_ab = I_NAI_Ix2y3z_S_C1003_ab+ABX*I_NAI_H2y3z_S_C1003_ab;
  Double I_NAI_Hy4z_Px_C1003_ab = I_NAI_Ixy4z_S_C1003_ab+ABX*I_NAI_Hy4z_S_C1003_ab;
  Double I_NAI_H5z_Px_C1003_ab = I_NAI_Ix5z_S_C1003_ab+ABX*I_NAI_H5z_S_C1003_ab;
  Double I_NAI_H4xy_Py_C1003_ab = I_NAI_I4x2y_S_C1003_ab+ABY*I_NAI_H4xy_S_C1003_ab;
  Double I_NAI_H4xz_Py_C1003_ab = I_NAI_I4xyz_S_C1003_ab+ABY*I_NAI_H4xz_S_C1003_ab;
  Double I_NAI_H3x2y_Py_C1003_ab = I_NAI_I3x3y_S_C1003_ab+ABY*I_NAI_H3x2y_S_C1003_ab;
  Double I_NAI_H3xyz_Py_C1003_ab = I_NAI_I3x2yz_S_C1003_ab+ABY*I_NAI_H3xyz_S_C1003_ab;
  Double I_NAI_H3x2z_Py_C1003_ab = I_NAI_I3xy2z_S_C1003_ab+ABY*I_NAI_H3x2z_S_C1003_ab;
  Double I_NAI_H2x3y_Py_C1003_ab = I_NAI_I2x4y_S_C1003_ab+ABY*I_NAI_H2x3y_S_C1003_ab;
  Double I_NAI_H2x2yz_Py_C1003_ab = I_NAI_I2x3yz_S_C1003_ab+ABY*I_NAI_H2x2yz_S_C1003_ab;
  Double I_NAI_H2xy2z_Py_C1003_ab = I_NAI_I2x2y2z_S_C1003_ab+ABY*I_NAI_H2xy2z_S_C1003_ab;
  Double I_NAI_H2x3z_Py_C1003_ab = I_NAI_I2xy3z_S_C1003_ab+ABY*I_NAI_H2x3z_S_C1003_ab;
  Double I_NAI_Hx4y_Py_C1003_ab = I_NAI_Ix5y_S_C1003_ab+ABY*I_NAI_Hx4y_S_C1003_ab;
  Double I_NAI_Hx3yz_Py_C1003_ab = I_NAI_Ix4yz_S_C1003_ab+ABY*I_NAI_Hx3yz_S_C1003_ab;
  Double I_NAI_Hx2y2z_Py_C1003_ab = I_NAI_Ix3y2z_S_C1003_ab+ABY*I_NAI_Hx2y2z_S_C1003_ab;
  Double I_NAI_Hxy3z_Py_C1003_ab = I_NAI_Ix2y3z_S_C1003_ab+ABY*I_NAI_Hxy3z_S_C1003_ab;
  Double I_NAI_Hx4z_Py_C1003_ab = I_NAI_Ixy4z_S_C1003_ab+ABY*I_NAI_Hx4z_S_C1003_ab;
  Double I_NAI_H5y_Py_C1003_ab = I_NAI_I6y_S_C1003_ab+ABY*I_NAI_H5y_S_C1003_ab;
  Double I_NAI_H4yz_Py_C1003_ab = I_NAI_I5yz_S_C1003_ab+ABY*I_NAI_H4yz_S_C1003_ab;
  Double I_NAI_H3y2z_Py_C1003_ab = I_NAI_I4y2z_S_C1003_ab+ABY*I_NAI_H3y2z_S_C1003_ab;
  Double I_NAI_H2y3z_Py_C1003_ab = I_NAI_I3y3z_S_C1003_ab+ABY*I_NAI_H2y3z_S_C1003_ab;
  Double I_NAI_Hy4z_Py_C1003_ab = I_NAI_I2y4z_S_C1003_ab+ABY*I_NAI_Hy4z_S_C1003_ab;
  Double I_NAI_H5z_Py_C1003_ab = I_NAI_Iy5z_S_C1003_ab+ABY*I_NAI_H5z_S_C1003_ab;
  Double I_NAI_H4xz_Pz_C1003_ab = I_NAI_I4x2z_S_C1003_ab+ABZ*I_NAI_H4xz_S_C1003_ab;
  Double I_NAI_H3xyz_Pz_C1003_ab = I_NAI_I3xy2z_S_C1003_ab+ABZ*I_NAI_H3xyz_S_C1003_ab;
  Double I_NAI_H3x2z_Pz_C1003_ab = I_NAI_I3x3z_S_C1003_ab+ABZ*I_NAI_H3x2z_S_C1003_ab;
  Double I_NAI_H2x2yz_Pz_C1003_ab = I_NAI_I2x2y2z_S_C1003_ab+ABZ*I_NAI_H2x2yz_S_C1003_ab;
  Double I_NAI_H2xy2z_Pz_C1003_ab = I_NAI_I2xy3z_S_C1003_ab+ABZ*I_NAI_H2xy2z_S_C1003_ab;
  Double I_NAI_H2x3z_Pz_C1003_ab = I_NAI_I2x4z_S_C1003_ab+ABZ*I_NAI_H2x3z_S_C1003_ab;
  Double I_NAI_Hx3yz_Pz_C1003_ab = I_NAI_Ix3y2z_S_C1003_ab+ABZ*I_NAI_Hx3yz_S_C1003_ab;
  Double I_NAI_Hx2y2z_Pz_C1003_ab = I_NAI_Ix2y3z_S_C1003_ab+ABZ*I_NAI_Hx2y2z_S_C1003_ab;
  Double I_NAI_Hxy3z_Pz_C1003_ab = I_NAI_Ixy4z_S_C1003_ab+ABZ*I_NAI_Hxy3z_S_C1003_ab;
  Double I_NAI_Hx4z_Pz_C1003_ab = I_NAI_Ix5z_S_C1003_ab+ABZ*I_NAI_Hx4z_S_C1003_ab;
  Double I_NAI_H4yz_Pz_C1003_ab = I_NAI_I4y2z_S_C1003_ab+ABZ*I_NAI_H4yz_S_C1003_ab;
  Double I_NAI_H3y2z_Pz_C1003_ab = I_NAI_I3y3z_S_C1003_ab+ABZ*I_NAI_H3y2z_S_C1003_ab;
  Double I_NAI_H2y3z_Pz_C1003_ab = I_NAI_I2y4z_S_C1003_ab+ABZ*I_NAI_H2y3z_S_C1003_ab;
  Double I_NAI_Hy4z_Pz_C1003_ab = I_NAI_Iy5z_S_C1003_ab+ABZ*I_NAI_Hy4z_S_C1003_ab;
  Double I_NAI_H5z_Pz_C1003_ab = I_NAI_I6z_S_C1003_ab+ABZ*I_NAI_H5z_S_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_C1003_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_P_C1003_ab
   ************************************************************/
  Double I_NAI_G4x_D2x_C1003_ab = I_NAI_H5x_Px_C1003_ab+ABX*I_NAI_G4x_Px_C1003_ab;
  Double I_NAI_G3xy_D2x_C1003_ab = I_NAI_H4xy_Px_C1003_ab+ABX*I_NAI_G3xy_Px_C1003_ab;
  Double I_NAI_G3xz_D2x_C1003_ab = I_NAI_H4xz_Px_C1003_ab+ABX*I_NAI_G3xz_Px_C1003_ab;
  Double I_NAI_G2x2y_D2x_C1003_ab = I_NAI_H3x2y_Px_C1003_ab+ABX*I_NAI_G2x2y_Px_C1003_ab;
  Double I_NAI_G2xyz_D2x_C1003_ab = I_NAI_H3xyz_Px_C1003_ab+ABX*I_NAI_G2xyz_Px_C1003_ab;
  Double I_NAI_G2x2z_D2x_C1003_ab = I_NAI_H3x2z_Px_C1003_ab+ABX*I_NAI_G2x2z_Px_C1003_ab;
  Double I_NAI_Gx3y_D2x_C1003_ab = I_NAI_H2x3y_Px_C1003_ab+ABX*I_NAI_Gx3y_Px_C1003_ab;
  Double I_NAI_Gx2yz_D2x_C1003_ab = I_NAI_H2x2yz_Px_C1003_ab+ABX*I_NAI_Gx2yz_Px_C1003_ab;
  Double I_NAI_Gxy2z_D2x_C1003_ab = I_NAI_H2xy2z_Px_C1003_ab+ABX*I_NAI_Gxy2z_Px_C1003_ab;
  Double I_NAI_Gx3z_D2x_C1003_ab = I_NAI_H2x3z_Px_C1003_ab+ABX*I_NAI_Gx3z_Px_C1003_ab;
  Double I_NAI_G4y_D2x_C1003_ab = I_NAI_Hx4y_Px_C1003_ab+ABX*I_NAI_G4y_Px_C1003_ab;
  Double I_NAI_G3yz_D2x_C1003_ab = I_NAI_Hx3yz_Px_C1003_ab+ABX*I_NAI_G3yz_Px_C1003_ab;
  Double I_NAI_G2y2z_D2x_C1003_ab = I_NAI_Hx2y2z_Px_C1003_ab+ABX*I_NAI_G2y2z_Px_C1003_ab;
  Double I_NAI_Gy3z_D2x_C1003_ab = I_NAI_Hxy3z_Px_C1003_ab+ABX*I_NAI_Gy3z_Px_C1003_ab;
  Double I_NAI_G4z_D2x_C1003_ab = I_NAI_Hx4z_Px_C1003_ab+ABX*I_NAI_G4z_Px_C1003_ab;
  Double I_NAI_G4x_Dxy_C1003_ab = I_NAI_H4xy_Px_C1003_ab+ABY*I_NAI_G4x_Px_C1003_ab;
  Double I_NAI_G3xy_Dxy_C1003_ab = I_NAI_H3x2y_Px_C1003_ab+ABY*I_NAI_G3xy_Px_C1003_ab;
  Double I_NAI_G3xz_Dxy_C1003_ab = I_NAI_H3xyz_Px_C1003_ab+ABY*I_NAI_G3xz_Px_C1003_ab;
  Double I_NAI_G2x2y_Dxy_C1003_ab = I_NAI_H2x3y_Px_C1003_ab+ABY*I_NAI_G2x2y_Px_C1003_ab;
  Double I_NAI_G2xyz_Dxy_C1003_ab = I_NAI_H2x2yz_Px_C1003_ab+ABY*I_NAI_G2xyz_Px_C1003_ab;
  Double I_NAI_G2x2z_Dxy_C1003_ab = I_NAI_H2xy2z_Px_C1003_ab+ABY*I_NAI_G2x2z_Px_C1003_ab;
  Double I_NAI_Gx3y_Dxy_C1003_ab = I_NAI_Hx4y_Px_C1003_ab+ABY*I_NAI_Gx3y_Px_C1003_ab;
  Double I_NAI_Gx2yz_Dxy_C1003_ab = I_NAI_Hx3yz_Px_C1003_ab+ABY*I_NAI_Gx2yz_Px_C1003_ab;
  Double I_NAI_Gxy2z_Dxy_C1003_ab = I_NAI_Hx2y2z_Px_C1003_ab+ABY*I_NAI_Gxy2z_Px_C1003_ab;
  Double I_NAI_Gx3z_Dxy_C1003_ab = I_NAI_Hxy3z_Px_C1003_ab+ABY*I_NAI_Gx3z_Px_C1003_ab;
  Double I_NAI_G4y_Dxy_C1003_ab = I_NAI_H5y_Px_C1003_ab+ABY*I_NAI_G4y_Px_C1003_ab;
  Double I_NAI_G3yz_Dxy_C1003_ab = I_NAI_H4yz_Px_C1003_ab+ABY*I_NAI_G3yz_Px_C1003_ab;
  Double I_NAI_G2y2z_Dxy_C1003_ab = I_NAI_H3y2z_Px_C1003_ab+ABY*I_NAI_G2y2z_Px_C1003_ab;
  Double I_NAI_Gy3z_Dxy_C1003_ab = I_NAI_H2y3z_Px_C1003_ab+ABY*I_NAI_Gy3z_Px_C1003_ab;
  Double I_NAI_G4z_Dxy_C1003_ab = I_NAI_Hy4z_Px_C1003_ab+ABY*I_NAI_G4z_Px_C1003_ab;
  Double I_NAI_G4x_Dxz_C1003_ab = I_NAI_H4xz_Px_C1003_ab+ABZ*I_NAI_G4x_Px_C1003_ab;
  Double I_NAI_G3xy_Dxz_C1003_ab = I_NAI_H3xyz_Px_C1003_ab+ABZ*I_NAI_G3xy_Px_C1003_ab;
  Double I_NAI_G3xz_Dxz_C1003_ab = I_NAI_H3x2z_Px_C1003_ab+ABZ*I_NAI_G3xz_Px_C1003_ab;
  Double I_NAI_G2x2y_Dxz_C1003_ab = I_NAI_H2x2yz_Px_C1003_ab+ABZ*I_NAI_G2x2y_Px_C1003_ab;
  Double I_NAI_G2xyz_Dxz_C1003_ab = I_NAI_H2xy2z_Px_C1003_ab+ABZ*I_NAI_G2xyz_Px_C1003_ab;
  Double I_NAI_G2x2z_Dxz_C1003_ab = I_NAI_H2x3z_Px_C1003_ab+ABZ*I_NAI_G2x2z_Px_C1003_ab;
  Double I_NAI_Gx3y_Dxz_C1003_ab = I_NAI_Hx3yz_Px_C1003_ab+ABZ*I_NAI_Gx3y_Px_C1003_ab;
  Double I_NAI_Gx2yz_Dxz_C1003_ab = I_NAI_Hx2y2z_Px_C1003_ab+ABZ*I_NAI_Gx2yz_Px_C1003_ab;
  Double I_NAI_Gxy2z_Dxz_C1003_ab = I_NAI_Hxy3z_Px_C1003_ab+ABZ*I_NAI_Gxy2z_Px_C1003_ab;
  Double I_NAI_Gx3z_Dxz_C1003_ab = I_NAI_Hx4z_Px_C1003_ab+ABZ*I_NAI_Gx3z_Px_C1003_ab;
  Double I_NAI_G4y_Dxz_C1003_ab = I_NAI_H4yz_Px_C1003_ab+ABZ*I_NAI_G4y_Px_C1003_ab;
  Double I_NAI_G3yz_Dxz_C1003_ab = I_NAI_H3y2z_Px_C1003_ab+ABZ*I_NAI_G3yz_Px_C1003_ab;
  Double I_NAI_G2y2z_Dxz_C1003_ab = I_NAI_H2y3z_Px_C1003_ab+ABZ*I_NAI_G2y2z_Px_C1003_ab;
  Double I_NAI_Gy3z_Dxz_C1003_ab = I_NAI_Hy4z_Px_C1003_ab+ABZ*I_NAI_Gy3z_Px_C1003_ab;
  Double I_NAI_G4z_Dxz_C1003_ab = I_NAI_H5z_Px_C1003_ab+ABZ*I_NAI_G4z_Px_C1003_ab;
  Double I_NAI_G4x_D2y_C1003_ab = I_NAI_H4xy_Py_C1003_ab+ABY*I_NAI_G4x_Py_C1003_ab;
  Double I_NAI_G3xy_D2y_C1003_ab = I_NAI_H3x2y_Py_C1003_ab+ABY*I_NAI_G3xy_Py_C1003_ab;
  Double I_NAI_G3xz_D2y_C1003_ab = I_NAI_H3xyz_Py_C1003_ab+ABY*I_NAI_G3xz_Py_C1003_ab;
  Double I_NAI_G2x2y_D2y_C1003_ab = I_NAI_H2x3y_Py_C1003_ab+ABY*I_NAI_G2x2y_Py_C1003_ab;
  Double I_NAI_G2xyz_D2y_C1003_ab = I_NAI_H2x2yz_Py_C1003_ab+ABY*I_NAI_G2xyz_Py_C1003_ab;
  Double I_NAI_G2x2z_D2y_C1003_ab = I_NAI_H2xy2z_Py_C1003_ab+ABY*I_NAI_G2x2z_Py_C1003_ab;
  Double I_NAI_Gx3y_D2y_C1003_ab = I_NAI_Hx4y_Py_C1003_ab+ABY*I_NAI_Gx3y_Py_C1003_ab;
  Double I_NAI_Gx2yz_D2y_C1003_ab = I_NAI_Hx3yz_Py_C1003_ab+ABY*I_NAI_Gx2yz_Py_C1003_ab;
  Double I_NAI_Gxy2z_D2y_C1003_ab = I_NAI_Hx2y2z_Py_C1003_ab+ABY*I_NAI_Gxy2z_Py_C1003_ab;
  Double I_NAI_Gx3z_D2y_C1003_ab = I_NAI_Hxy3z_Py_C1003_ab+ABY*I_NAI_Gx3z_Py_C1003_ab;
  Double I_NAI_G4y_D2y_C1003_ab = I_NAI_H5y_Py_C1003_ab+ABY*I_NAI_G4y_Py_C1003_ab;
  Double I_NAI_G3yz_D2y_C1003_ab = I_NAI_H4yz_Py_C1003_ab+ABY*I_NAI_G3yz_Py_C1003_ab;
  Double I_NAI_G2y2z_D2y_C1003_ab = I_NAI_H3y2z_Py_C1003_ab+ABY*I_NAI_G2y2z_Py_C1003_ab;
  Double I_NAI_Gy3z_D2y_C1003_ab = I_NAI_H2y3z_Py_C1003_ab+ABY*I_NAI_Gy3z_Py_C1003_ab;
  Double I_NAI_G4z_D2y_C1003_ab = I_NAI_Hy4z_Py_C1003_ab+ABY*I_NAI_G4z_Py_C1003_ab;
  Double I_NAI_G4x_Dyz_C1003_ab = I_NAI_H4xz_Py_C1003_ab+ABZ*I_NAI_G4x_Py_C1003_ab;
  Double I_NAI_G3xy_Dyz_C1003_ab = I_NAI_H3xyz_Py_C1003_ab+ABZ*I_NAI_G3xy_Py_C1003_ab;
  Double I_NAI_G3xz_Dyz_C1003_ab = I_NAI_H3x2z_Py_C1003_ab+ABZ*I_NAI_G3xz_Py_C1003_ab;
  Double I_NAI_G2x2y_Dyz_C1003_ab = I_NAI_H2x2yz_Py_C1003_ab+ABZ*I_NAI_G2x2y_Py_C1003_ab;
  Double I_NAI_G2xyz_Dyz_C1003_ab = I_NAI_H2xy2z_Py_C1003_ab+ABZ*I_NAI_G2xyz_Py_C1003_ab;
  Double I_NAI_G2x2z_Dyz_C1003_ab = I_NAI_H2x3z_Py_C1003_ab+ABZ*I_NAI_G2x2z_Py_C1003_ab;
  Double I_NAI_Gx3y_Dyz_C1003_ab = I_NAI_Hx3yz_Py_C1003_ab+ABZ*I_NAI_Gx3y_Py_C1003_ab;
  Double I_NAI_Gx2yz_Dyz_C1003_ab = I_NAI_Hx2y2z_Py_C1003_ab+ABZ*I_NAI_Gx2yz_Py_C1003_ab;
  Double I_NAI_Gxy2z_Dyz_C1003_ab = I_NAI_Hxy3z_Py_C1003_ab+ABZ*I_NAI_Gxy2z_Py_C1003_ab;
  Double I_NAI_Gx3z_Dyz_C1003_ab = I_NAI_Hx4z_Py_C1003_ab+ABZ*I_NAI_Gx3z_Py_C1003_ab;
  Double I_NAI_G4y_Dyz_C1003_ab = I_NAI_H4yz_Py_C1003_ab+ABZ*I_NAI_G4y_Py_C1003_ab;
  Double I_NAI_G3yz_Dyz_C1003_ab = I_NAI_H3y2z_Py_C1003_ab+ABZ*I_NAI_G3yz_Py_C1003_ab;
  Double I_NAI_G2y2z_Dyz_C1003_ab = I_NAI_H2y3z_Py_C1003_ab+ABZ*I_NAI_G2y2z_Py_C1003_ab;
  Double I_NAI_Gy3z_Dyz_C1003_ab = I_NAI_Hy4z_Py_C1003_ab+ABZ*I_NAI_Gy3z_Py_C1003_ab;
  Double I_NAI_G4z_Dyz_C1003_ab = I_NAI_H5z_Py_C1003_ab+ABZ*I_NAI_G4z_Py_C1003_ab;
  Double I_NAI_G4x_D2z_C1003_ab = I_NAI_H4xz_Pz_C1003_ab+ABZ*I_NAI_G4x_Pz_C1003_ab;
  Double I_NAI_G3xy_D2z_C1003_ab = I_NAI_H3xyz_Pz_C1003_ab+ABZ*I_NAI_G3xy_Pz_C1003_ab;
  Double I_NAI_G3xz_D2z_C1003_ab = I_NAI_H3x2z_Pz_C1003_ab+ABZ*I_NAI_G3xz_Pz_C1003_ab;
  Double I_NAI_G2x2y_D2z_C1003_ab = I_NAI_H2x2yz_Pz_C1003_ab+ABZ*I_NAI_G2x2y_Pz_C1003_ab;
  Double I_NAI_G2xyz_D2z_C1003_ab = I_NAI_H2xy2z_Pz_C1003_ab+ABZ*I_NAI_G2xyz_Pz_C1003_ab;
  Double I_NAI_G2x2z_D2z_C1003_ab = I_NAI_H2x3z_Pz_C1003_ab+ABZ*I_NAI_G2x2z_Pz_C1003_ab;
  Double I_NAI_Gx3y_D2z_C1003_ab = I_NAI_Hx3yz_Pz_C1003_ab+ABZ*I_NAI_Gx3y_Pz_C1003_ab;
  Double I_NAI_Gx2yz_D2z_C1003_ab = I_NAI_Hx2y2z_Pz_C1003_ab+ABZ*I_NAI_Gx2yz_Pz_C1003_ab;
  Double I_NAI_Gxy2z_D2z_C1003_ab = I_NAI_Hxy3z_Pz_C1003_ab+ABZ*I_NAI_Gxy2z_Pz_C1003_ab;
  Double I_NAI_Gx3z_D2z_C1003_ab = I_NAI_Hx4z_Pz_C1003_ab+ABZ*I_NAI_Gx3z_Pz_C1003_ab;
  Double I_NAI_G4y_D2z_C1003_ab = I_NAI_H4yz_Pz_C1003_ab+ABZ*I_NAI_G4y_Pz_C1003_ab;
  Double I_NAI_G3yz_D2z_C1003_ab = I_NAI_H3y2z_Pz_C1003_ab+ABZ*I_NAI_G3yz_Pz_C1003_ab;
  Double I_NAI_G2y2z_D2z_C1003_ab = I_NAI_H2y3z_Pz_C1003_ab+ABZ*I_NAI_G2y2z_Pz_C1003_ab;
  Double I_NAI_Gy3z_D2z_C1003_ab = I_NAI_Hy4z_Pz_C1003_ab+ABZ*I_NAI_Gy3z_Pz_C1003_ab;
  Double I_NAI_G4z_D2z_C1003_ab = I_NAI_H5z_Pz_C1003_ab+ABZ*I_NAI_G4z_Pz_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C3_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_bb
   ************************************************************/
  Double I_NAI_F3x_Px_C3_bb = I_NAI_G4x_S_C3_bb+ABX*I_NAI_F3x_S_C3_bb;
  Double I_NAI_F2xy_Px_C3_bb = I_NAI_G3xy_S_C3_bb+ABX*I_NAI_F2xy_S_C3_bb;
  Double I_NAI_F2xz_Px_C3_bb = I_NAI_G3xz_S_C3_bb+ABX*I_NAI_F2xz_S_C3_bb;
  Double I_NAI_Fx2y_Px_C3_bb = I_NAI_G2x2y_S_C3_bb+ABX*I_NAI_Fx2y_S_C3_bb;
  Double I_NAI_Fxyz_Px_C3_bb = I_NAI_G2xyz_S_C3_bb+ABX*I_NAI_Fxyz_S_C3_bb;
  Double I_NAI_Fx2z_Px_C3_bb = I_NAI_G2x2z_S_C3_bb+ABX*I_NAI_Fx2z_S_C3_bb;
  Double I_NAI_F3y_Px_C3_bb = I_NAI_Gx3y_S_C3_bb+ABX*I_NAI_F3y_S_C3_bb;
  Double I_NAI_F2yz_Px_C3_bb = I_NAI_Gx2yz_S_C3_bb+ABX*I_NAI_F2yz_S_C3_bb;
  Double I_NAI_Fy2z_Px_C3_bb = I_NAI_Gxy2z_S_C3_bb+ABX*I_NAI_Fy2z_S_C3_bb;
  Double I_NAI_F3z_Px_C3_bb = I_NAI_Gx3z_S_C3_bb+ABX*I_NAI_F3z_S_C3_bb;
  Double I_NAI_F3x_Py_C3_bb = I_NAI_G3xy_S_C3_bb+ABY*I_NAI_F3x_S_C3_bb;
  Double I_NAI_F2xy_Py_C3_bb = I_NAI_G2x2y_S_C3_bb+ABY*I_NAI_F2xy_S_C3_bb;
  Double I_NAI_F2xz_Py_C3_bb = I_NAI_G2xyz_S_C3_bb+ABY*I_NAI_F2xz_S_C3_bb;
  Double I_NAI_Fx2y_Py_C3_bb = I_NAI_Gx3y_S_C3_bb+ABY*I_NAI_Fx2y_S_C3_bb;
  Double I_NAI_Fxyz_Py_C3_bb = I_NAI_Gx2yz_S_C3_bb+ABY*I_NAI_Fxyz_S_C3_bb;
  Double I_NAI_Fx2z_Py_C3_bb = I_NAI_Gxy2z_S_C3_bb+ABY*I_NAI_Fx2z_S_C3_bb;
  Double I_NAI_F3y_Py_C3_bb = I_NAI_G4y_S_C3_bb+ABY*I_NAI_F3y_S_C3_bb;
  Double I_NAI_F2yz_Py_C3_bb = I_NAI_G3yz_S_C3_bb+ABY*I_NAI_F2yz_S_C3_bb;
  Double I_NAI_Fy2z_Py_C3_bb = I_NAI_G2y2z_S_C3_bb+ABY*I_NAI_Fy2z_S_C3_bb;
  Double I_NAI_F3z_Py_C3_bb = I_NAI_Gy3z_S_C3_bb+ABY*I_NAI_F3z_S_C3_bb;
  Double I_NAI_F3x_Pz_C3_bb = I_NAI_G3xz_S_C3_bb+ABZ*I_NAI_F3x_S_C3_bb;
  Double I_NAI_F2xy_Pz_C3_bb = I_NAI_G2xyz_S_C3_bb+ABZ*I_NAI_F2xy_S_C3_bb;
  Double I_NAI_F2xz_Pz_C3_bb = I_NAI_G2x2z_S_C3_bb+ABZ*I_NAI_F2xz_S_C3_bb;
  Double I_NAI_Fx2y_Pz_C3_bb = I_NAI_Gx2yz_S_C3_bb+ABZ*I_NAI_Fx2y_S_C3_bb;
  Double I_NAI_Fxyz_Pz_C3_bb = I_NAI_Gxy2z_S_C3_bb+ABZ*I_NAI_Fxyz_S_C3_bb;
  Double I_NAI_Fx2z_Pz_C3_bb = I_NAI_Gx3z_S_C3_bb+ABZ*I_NAI_Fx2z_S_C3_bb;
  Double I_NAI_F3y_Pz_C3_bb = I_NAI_G3yz_S_C3_bb+ABZ*I_NAI_F3y_S_C3_bb;
  Double I_NAI_F2yz_Pz_C3_bb = I_NAI_G2y2z_S_C3_bb+ABZ*I_NAI_F2yz_S_C3_bb;
  Double I_NAI_Fy2z_Pz_C3_bb = I_NAI_Gy3z_S_C3_bb+ABZ*I_NAI_Fy2z_S_C3_bb;
  Double I_NAI_F3z_Pz_C3_bb = I_NAI_G4z_S_C3_bb+ABZ*I_NAI_F3z_S_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C3_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_bb
   * RHS shell quartet name: SQ_NAI_G_S_C3_bb
   ************************************************************/
  Double I_NAI_G4x_Px_C3_bb = I_NAI_H5x_S_C3_bb+ABX*I_NAI_G4x_S_C3_bb;
  Double I_NAI_G3xy_Px_C3_bb = I_NAI_H4xy_S_C3_bb+ABX*I_NAI_G3xy_S_C3_bb;
  Double I_NAI_G3xz_Px_C3_bb = I_NAI_H4xz_S_C3_bb+ABX*I_NAI_G3xz_S_C3_bb;
  Double I_NAI_G2x2y_Px_C3_bb = I_NAI_H3x2y_S_C3_bb+ABX*I_NAI_G2x2y_S_C3_bb;
  Double I_NAI_G2xyz_Px_C3_bb = I_NAI_H3xyz_S_C3_bb+ABX*I_NAI_G2xyz_S_C3_bb;
  Double I_NAI_G2x2z_Px_C3_bb = I_NAI_H3x2z_S_C3_bb+ABX*I_NAI_G2x2z_S_C3_bb;
  Double I_NAI_Gx3y_Px_C3_bb = I_NAI_H2x3y_S_C3_bb+ABX*I_NAI_Gx3y_S_C3_bb;
  Double I_NAI_Gx2yz_Px_C3_bb = I_NAI_H2x2yz_S_C3_bb+ABX*I_NAI_Gx2yz_S_C3_bb;
  Double I_NAI_Gxy2z_Px_C3_bb = I_NAI_H2xy2z_S_C3_bb+ABX*I_NAI_Gxy2z_S_C3_bb;
  Double I_NAI_Gx3z_Px_C3_bb = I_NAI_H2x3z_S_C3_bb+ABX*I_NAI_Gx3z_S_C3_bb;
  Double I_NAI_G4y_Px_C3_bb = I_NAI_Hx4y_S_C3_bb+ABX*I_NAI_G4y_S_C3_bb;
  Double I_NAI_G3yz_Px_C3_bb = I_NAI_Hx3yz_S_C3_bb+ABX*I_NAI_G3yz_S_C3_bb;
  Double I_NAI_G2y2z_Px_C3_bb = I_NAI_Hx2y2z_S_C3_bb+ABX*I_NAI_G2y2z_S_C3_bb;
  Double I_NAI_Gy3z_Px_C3_bb = I_NAI_Hxy3z_S_C3_bb+ABX*I_NAI_Gy3z_S_C3_bb;
  Double I_NAI_G4z_Px_C3_bb = I_NAI_Hx4z_S_C3_bb+ABX*I_NAI_G4z_S_C3_bb;
  Double I_NAI_G3xy_Py_C3_bb = I_NAI_H3x2y_S_C3_bb+ABY*I_NAI_G3xy_S_C3_bb;
  Double I_NAI_G3xz_Py_C3_bb = I_NAI_H3xyz_S_C3_bb+ABY*I_NAI_G3xz_S_C3_bb;
  Double I_NAI_G2x2y_Py_C3_bb = I_NAI_H2x3y_S_C3_bb+ABY*I_NAI_G2x2y_S_C3_bb;
  Double I_NAI_G2xyz_Py_C3_bb = I_NAI_H2x2yz_S_C3_bb+ABY*I_NAI_G2xyz_S_C3_bb;
  Double I_NAI_G2x2z_Py_C3_bb = I_NAI_H2xy2z_S_C3_bb+ABY*I_NAI_G2x2z_S_C3_bb;
  Double I_NAI_Gx3y_Py_C3_bb = I_NAI_Hx4y_S_C3_bb+ABY*I_NAI_Gx3y_S_C3_bb;
  Double I_NAI_Gx2yz_Py_C3_bb = I_NAI_Hx3yz_S_C3_bb+ABY*I_NAI_Gx2yz_S_C3_bb;
  Double I_NAI_Gxy2z_Py_C3_bb = I_NAI_Hx2y2z_S_C3_bb+ABY*I_NAI_Gxy2z_S_C3_bb;
  Double I_NAI_Gx3z_Py_C3_bb = I_NAI_Hxy3z_S_C3_bb+ABY*I_NAI_Gx3z_S_C3_bb;
  Double I_NAI_G4y_Py_C3_bb = I_NAI_H5y_S_C3_bb+ABY*I_NAI_G4y_S_C3_bb;
  Double I_NAI_G3yz_Py_C3_bb = I_NAI_H4yz_S_C3_bb+ABY*I_NAI_G3yz_S_C3_bb;
  Double I_NAI_G2y2z_Py_C3_bb = I_NAI_H3y2z_S_C3_bb+ABY*I_NAI_G2y2z_S_C3_bb;
  Double I_NAI_Gy3z_Py_C3_bb = I_NAI_H2y3z_S_C3_bb+ABY*I_NAI_Gy3z_S_C3_bb;
  Double I_NAI_G4z_Py_C3_bb = I_NAI_Hy4z_S_C3_bb+ABY*I_NAI_G4z_S_C3_bb;
  Double I_NAI_G3xz_Pz_C3_bb = I_NAI_H3x2z_S_C3_bb+ABZ*I_NAI_G3xz_S_C3_bb;
  Double I_NAI_G2xyz_Pz_C3_bb = I_NAI_H2xy2z_S_C3_bb+ABZ*I_NAI_G2xyz_S_C3_bb;
  Double I_NAI_G2x2z_Pz_C3_bb = I_NAI_H2x3z_S_C3_bb+ABZ*I_NAI_G2x2z_S_C3_bb;
  Double I_NAI_Gx2yz_Pz_C3_bb = I_NAI_Hx2y2z_S_C3_bb+ABZ*I_NAI_Gx2yz_S_C3_bb;
  Double I_NAI_Gxy2z_Pz_C3_bb = I_NAI_Hxy3z_S_C3_bb+ABZ*I_NAI_Gxy2z_S_C3_bb;
  Double I_NAI_Gx3z_Pz_C3_bb = I_NAI_Hx4z_S_C3_bb+ABZ*I_NAI_Gx3z_S_C3_bb;
  Double I_NAI_G3yz_Pz_C3_bb = I_NAI_H3y2z_S_C3_bb+ABZ*I_NAI_G3yz_S_C3_bb;
  Double I_NAI_G2y2z_Pz_C3_bb = I_NAI_H2y3z_S_C3_bb+ABZ*I_NAI_G2y2z_S_C3_bb;
  Double I_NAI_Gy3z_Pz_C3_bb = I_NAI_Hy4z_S_C3_bb+ABZ*I_NAI_Gy3z_S_C3_bb;
  Double I_NAI_G4z_Pz_C3_bb = I_NAI_H5z_S_C3_bb+ABZ*I_NAI_G4z_S_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_C3_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_bb
   * RHS shell quartet name: SQ_NAI_F_P_C3_bb
   ************************************************************/
  Double I_NAI_F3x_D2x_C3_bb = I_NAI_G4x_Px_C3_bb+ABX*I_NAI_F3x_Px_C3_bb;
  Double I_NAI_F2xy_D2x_C3_bb = I_NAI_G3xy_Px_C3_bb+ABX*I_NAI_F2xy_Px_C3_bb;
  Double I_NAI_F2xz_D2x_C3_bb = I_NAI_G3xz_Px_C3_bb+ABX*I_NAI_F2xz_Px_C3_bb;
  Double I_NAI_Fx2y_D2x_C3_bb = I_NAI_G2x2y_Px_C3_bb+ABX*I_NAI_Fx2y_Px_C3_bb;
  Double I_NAI_Fxyz_D2x_C3_bb = I_NAI_G2xyz_Px_C3_bb+ABX*I_NAI_Fxyz_Px_C3_bb;
  Double I_NAI_Fx2z_D2x_C3_bb = I_NAI_G2x2z_Px_C3_bb+ABX*I_NAI_Fx2z_Px_C3_bb;
  Double I_NAI_F3y_D2x_C3_bb = I_NAI_Gx3y_Px_C3_bb+ABX*I_NAI_F3y_Px_C3_bb;
  Double I_NAI_F2yz_D2x_C3_bb = I_NAI_Gx2yz_Px_C3_bb+ABX*I_NAI_F2yz_Px_C3_bb;
  Double I_NAI_Fy2z_D2x_C3_bb = I_NAI_Gxy2z_Px_C3_bb+ABX*I_NAI_Fy2z_Px_C3_bb;
  Double I_NAI_F3z_D2x_C3_bb = I_NAI_Gx3z_Px_C3_bb+ABX*I_NAI_F3z_Px_C3_bb;
  Double I_NAI_F3x_Dxy_C3_bb = I_NAI_G3xy_Px_C3_bb+ABY*I_NAI_F3x_Px_C3_bb;
  Double I_NAI_F2xy_Dxy_C3_bb = I_NAI_G2x2y_Px_C3_bb+ABY*I_NAI_F2xy_Px_C3_bb;
  Double I_NAI_F2xz_Dxy_C3_bb = I_NAI_G2xyz_Px_C3_bb+ABY*I_NAI_F2xz_Px_C3_bb;
  Double I_NAI_Fx2y_Dxy_C3_bb = I_NAI_Gx3y_Px_C3_bb+ABY*I_NAI_Fx2y_Px_C3_bb;
  Double I_NAI_Fxyz_Dxy_C3_bb = I_NAI_Gx2yz_Px_C3_bb+ABY*I_NAI_Fxyz_Px_C3_bb;
  Double I_NAI_Fx2z_Dxy_C3_bb = I_NAI_Gxy2z_Px_C3_bb+ABY*I_NAI_Fx2z_Px_C3_bb;
  Double I_NAI_F3y_Dxy_C3_bb = I_NAI_G4y_Px_C3_bb+ABY*I_NAI_F3y_Px_C3_bb;
  Double I_NAI_F2yz_Dxy_C3_bb = I_NAI_G3yz_Px_C3_bb+ABY*I_NAI_F2yz_Px_C3_bb;
  Double I_NAI_Fy2z_Dxy_C3_bb = I_NAI_G2y2z_Px_C3_bb+ABY*I_NAI_Fy2z_Px_C3_bb;
  Double I_NAI_F3z_Dxy_C3_bb = I_NAI_Gy3z_Px_C3_bb+ABY*I_NAI_F3z_Px_C3_bb;
  Double I_NAI_F3x_Dxz_C3_bb = I_NAI_G3xz_Px_C3_bb+ABZ*I_NAI_F3x_Px_C3_bb;
  Double I_NAI_F2xy_Dxz_C3_bb = I_NAI_G2xyz_Px_C3_bb+ABZ*I_NAI_F2xy_Px_C3_bb;
  Double I_NAI_F2xz_Dxz_C3_bb = I_NAI_G2x2z_Px_C3_bb+ABZ*I_NAI_F2xz_Px_C3_bb;
  Double I_NAI_Fx2y_Dxz_C3_bb = I_NAI_Gx2yz_Px_C3_bb+ABZ*I_NAI_Fx2y_Px_C3_bb;
  Double I_NAI_Fxyz_Dxz_C3_bb = I_NAI_Gxy2z_Px_C3_bb+ABZ*I_NAI_Fxyz_Px_C3_bb;
  Double I_NAI_Fx2z_Dxz_C3_bb = I_NAI_Gx3z_Px_C3_bb+ABZ*I_NAI_Fx2z_Px_C3_bb;
  Double I_NAI_F3y_Dxz_C3_bb = I_NAI_G3yz_Px_C3_bb+ABZ*I_NAI_F3y_Px_C3_bb;
  Double I_NAI_F2yz_Dxz_C3_bb = I_NAI_G2y2z_Px_C3_bb+ABZ*I_NAI_F2yz_Px_C3_bb;
  Double I_NAI_Fy2z_Dxz_C3_bb = I_NAI_Gy3z_Px_C3_bb+ABZ*I_NAI_Fy2z_Px_C3_bb;
  Double I_NAI_F3z_Dxz_C3_bb = I_NAI_G4z_Px_C3_bb+ABZ*I_NAI_F3z_Px_C3_bb;
  Double I_NAI_F3x_D2y_C3_bb = I_NAI_G3xy_Py_C3_bb+ABY*I_NAI_F3x_Py_C3_bb;
  Double I_NAI_F2xy_D2y_C3_bb = I_NAI_G2x2y_Py_C3_bb+ABY*I_NAI_F2xy_Py_C3_bb;
  Double I_NAI_F2xz_D2y_C3_bb = I_NAI_G2xyz_Py_C3_bb+ABY*I_NAI_F2xz_Py_C3_bb;
  Double I_NAI_Fx2y_D2y_C3_bb = I_NAI_Gx3y_Py_C3_bb+ABY*I_NAI_Fx2y_Py_C3_bb;
  Double I_NAI_Fxyz_D2y_C3_bb = I_NAI_Gx2yz_Py_C3_bb+ABY*I_NAI_Fxyz_Py_C3_bb;
  Double I_NAI_Fx2z_D2y_C3_bb = I_NAI_Gxy2z_Py_C3_bb+ABY*I_NAI_Fx2z_Py_C3_bb;
  Double I_NAI_F3y_D2y_C3_bb = I_NAI_G4y_Py_C3_bb+ABY*I_NAI_F3y_Py_C3_bb;
  Double I_NAI_F2yz_D2y_C3_bb = I_NAI_G3yz_Py_C3_bb+ABY*I_NAI_F2yz_Py_C3_bb;
  Double I_NAI_Fy2z_D2y_C3_bb = I_NAI_G2y2z_Py_C3_bb+ABY*I_NAI_Fy2z_Py_C3_bb;
  Double I_NAI_F3z_D2y_C3_bb = I_NAI_Gy3z_Py_C3_bb+ABY*I_NAI_F3z_Py_C3_bb;
  Double I_NAI_F3x_Dyz_C3_bb = I_NAI_G3xz_Py_C3_bb+ABZ*I_NAI_F3x_Py_C3_bb;
  Double I_NAI_F2xy_Dyz_C3_bb = I_NAI_G2xyz_Py_C3_bb+ABZ*I_NAI_F2xy_Py_C3_bb;
  Double I_NAI_F2xz_Dyz_C3_bb = I_NAI_G2x2z_Py_C3_bb+ABZ*I_NAI_F2xz_Py_C3_bb;
  Double I_NAI_Fx2y_Dyz_C3_bb = I_NAI_Gx2yz_Py_C3_bb+ABZ*I_NAI_Fx2y_Py_C3_bb;
  Double I_NAI_Fxyz_Dyz_C3_bb = I_NAI_Gxy2z_Py_C3_bb+ABZ*I_NAI_Fxyz_Py_C3_bb;
  Double I_NAI_Fx2z_Dyz_C3_bb = I_NAI_Gx3z_Py_C3_bb+ABZ*I_NAI_Fx2z_Py_C3_bb;
  Double I_NAI_F3y_Dyz_C3_bb = I_NAI_G3yz_Py_C3_bb+ABZ*I_NAI_F3y_Py_C3_bb;
  Double I_NAI_F2yz_Dyz_C3_bb = I_NAI_G2y2z_Py_C3_bb+ABZ*I_NAI_F2yz_Py_C3_bb;
  Double I_NAI_Fy2z_Dyz_C3_bb = I_NAI_Gy3z_Py_C3_bb+ABZ*I_NAI_Fy2z_Py_C3_bb;
  Double I_NAI_F3z_Dyz_C3_bb = I_NAI_G4z_Py_C3_bb+ABZ*I_NAI_F3z_Py_C3_bb;
  Double I_NAI_F3x_D2z_C3_bb = I_NAI_G3xz_Pz_C3_bb+ABZ*I_NAI_F3x_Pz_C3_bb;
  Double I_NAI_F2xy_D2z_C3_bb = I_NAI_G2xyz_Pz_C3_bb+ABZ*I_NAI_F2xy_Pz_C3_bb;
  Double I_NAI_F2xz_D2z_C3_bb = I_NAI_G2x2z_Pz_C3_bb+ABZ*I_NAI_F2xz_Pz_C3_bb;
  Double I_NAI_Fx2y_D2z_C3_bb = I_NAI_Gx2yz_Pz_C3_bb+ABZ*I_NAI_Fx2y_Pz_C3_bb;
  Double I_NAI_Fxyz_D2z_C3_bb = I_NAI_Gxy2z_Pz_C3_bb+ABZ*I_NAI_Fxyz_Pz_C3_bb;
  Double I_NAI_Fx2z_D2z_C3_bb = I_NAI_Gx3z_Pz_C3_bb+ABZ*I_NAI_Fx2z_Pz_C3_bb;
  Double I_NAI_F3y_D2z_C3_bb = I_NAI_G3yz_Pz_C3_bb+ABZ*I_NAI_F3y_Pz_C3_bb;
  Double I_NAI_F2yz_D2z_C3_bb = I_NAI_G2y2z_Pz_C3_bb+ABZ*I_NAI_F2yz_Pz_C3_bb;
  Double I_NAI_Fy2z_D2z_C3_bb = I_NAI_Gy3z_Pz_C3_bb+ABZ*I_NAI_Fy2z_Pz_C3_bb;
  Double I_NAI_F3z_D2z_C3_bb = I_NAI_G4z_Pz_C3_bb+ABZ*I_NAI_F3z_Pz_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_S_C1003_bb
   ************************************************************/
  Double I_NAI_F3x_Px_C1003_bb = I_NAI_G4x_S_C1003_bb+ABX*I_NAI_F3x_S_C1003_bb;
  Double I_NAI_F2xy_Px_C1003_bb = I_NAI_G3xy_S_C1003_bb+ABX*I_NAI_F2xy_S_C1003_bb;
  Double I_NAI_F2xz_Px_C1003_bb = I_NAI_G3xz_S_C1003_bb+ABX*I_NAI_F2xz_S_C1003_bb;
  Double I_NAI_Fx2y_Px_C1003_bb = I_NAI_G2x2y_S_C1003_bb+ABX*I_NAI_Fx2y_S_C1003_bb;
  Double I_NAI_Fxyz_Px_C1003_bb = I_NAI_G2xyz_S_C1003_bb+ABX*I_NAI_Fxyz_S_C1003_bb;
  Double I_NAI_Fx2z_Px_C1003_bb = I_NAI_G2x2z_S_C1003_bb+ABX*I_NAI_Fx2z_S_C1003_bb;
  Double I_NAI_F3y_Px_C1003_bb = I_NAI_Gx3y_S_C1003_bb+ABX*I_NAI_F3y_S_C1003_bb;
  Double I_NAI_F2yz_Px_C1003_bb = I_NAI_Gx2yz_S_C1003_bb+ABX*I_NAI_F2yz_S_C1003_bb;
  Double I_NAI_Fy2z_Px_C1003_bb = I_NAI_Gxy2z_S_C1003_bb+ABX*I_NAI_Fy2z_S_C1003_bb;
  Double I_NAI_F3z_Px_C1003_bb = I_NAI_Gx3z_S_C1003_bb+ABX*I_NAI_F3z_S_C1003_bb;
  Double I_NAI_F3x_Py_C1003_bb = I_NAI_G3xy_S_C1003_bb+ABY*I_NAI_F3x_S_C1003_bb;
  Double I_NAI_F2xy_Py_C1003_bb = I_NAI_G2x2y_S_C1003_bb+ABY*I_NAI_F2xy_S_C1003_bb;
  Double I_NAI_F2xz_Py_C1003_bb = I_NAI_G2xyz_S_C1003_bb+ABY*I_NAI_F2xz_S_C1003_bb;
  Double I_NAI_Fx2y_Py_C1003_bb = I_NAI_Gx3y_S_C1003_bb+ABY*I_NAI_Fx2y_S_C1003_bb;
  Double I_NAI_Fxyz_Py_C1003_bb = I_NAI_Gx2yz_S_C1003_bb+ABY*I_NAI_Fxyz_S_C1003_bb;
  Double I_NAI_Fx2z_Py_C1003_bb = I_NAI_Gxy2z_S_C1003_bb+ABY*I_NAI_Fx2z_S_C1003_bb;
  Double I_NAI_F3y_Py_C1003_bb = I_NAI_G4y_S_C1003_bb+ABY*I_NAI_F3y_S_C1003_bb;
  Double I_NAI_F2yz_Py_C1003_bb = I_NAI_G3yz_S_C1003_bb+ABY*I_NAI_F2yz_S_C1003_bb;
  Double I_NAI_Fy2z_Py_C1003_bb = I_NAI_G2y2z_S_C1003_bb+ABY*I_NAI_Fy2z_S_C1003_bb;
  Double I_NAI_F3z_Py_C1003_bb = I_NAI_Gy3z_S_C1003_bb+ABY*I_NAI_F3z_S_C1003_bb;
  Double I_NAI_F3x_Pz_C1003_bb = I_NAI_G3xz_S_C1003_bb+ABZ*I_NAI_F3x_S_C1003_bb;
  Double I_NAI_F2xy_Pz_C1003_bb = I_NAI_G2xyz_S_C1003_bb+ABZ*I_NAI_F2xy_S_C1003_bb;
  Double I_NAI_F2xz_Pz_C1003_bb = I_NAI_G2x2z_S_C1003_bb+ABZ*I_NAI_F2xz_S_C1003_bb;
  Double I_NAI_Fx2y_Pz_C1003_bb = I_NAI_Gx2yz_S_C1003_bb+ABZ*I_NAI_Fx2y_S_C1003_bb;
  Double I_NAI_Fxyz_Pz_C1003_bb = I_NAI_Gxy2z_S_C1003_bb+ABZ*I_NAI_Fxyz_S_C1003_bb;
  Double I_NAI_Fx2z_Pz_C1003_bb = I_NAI_Gx3z_S_C1003_bb+ABZ*I_NAI_Fx2z_S_C1003_bb;
  Double I_NAI_F3y_Pz_C1003_bb = I_NAI_G3yz_S_C1003_bb+ABZ*I_NAI_F3y_S_C1003_bb;
  Double I_NAI_F2yz_Pz_C1003_bb = I_NAI_G2y2z_S_C1003_bb+ABZ*I_NAI_F2yz_S_C1003_bb;
  Double I_NAI_Fy2z_Pz_C1003_bb = I_NAI_Gy3z_S_C1003_bb+ABZ*I_NAI_Fy2z_S_C1003_bb;
  Double I_NAI_F3z_Pz_C1003_bb = I_NAI_G4z_S_C1003_bb+ABZ*I_NAI_F3z_S_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C1003_bb
   * RHS shell quartet name: SQ_NAI_G_S_C1003_bb
   ************************************************************/
  Double I_NAI_G4x_Px_C1003_bb = I_NAI_H5x_S_C1003_bb+ABX*I_NAI_G4x_S_C1003_bb;
  Double I_NAI_G3xy_Px_C1003_bb = I_NAI_H4xy_S_C1003_bb+ABX*I_NAI_G3xy_S_C1003_bb;
  Double I_NAI_G3xz_Px_C1003_bb = I_NAI_H4xz_S_C1003_bb+ABX*I_NAI_G3xz_S_C1003_bb;
  Double I_NAI_G2x2y_Px_C1003_bb = I_NAI_H3x2y_S_C1003_bb+ABX*I_NAI_G2x2y_S_C1003_bb;
  Double I_NAI_G2xyz_Px_C1003_bb = I_NAI_H3xyz_S_C1003_bb+ABX*I_NAI_G2xyz_S_C1003_bb;
  Double I_NAI_G2x2z_Px_C1003_bb = I_NAI_H3x2z_S_C1003_bb+ABX*I_NAI_G2x2z_S_C1003_bb;
  Double I_NAI_Gx3y_Px_C1003_bb = I_NAI_H2x3y_S_C1003_bb+ABX*I_NAI_Gx3y_S_C1003_bb;
  Double I_NAI_Gx2yz_Px_C1003_bb = I_NAI_H2x2yz_S_C1003_bb+ABX*I_NAI_Gx2yz_S_C1003_bb;
  Double I_NAI_Gxy2z_Px_C1003_bb = I_NAI_H2xy2z_S_C1003_bb+ABX*I_NAI_Gxy2z_S_C1003_bb;
  Double I_NAI_Gx3z_Px_C1003_bb = I_NAI_H2x3z_S_C1003_bb+ABX*I_NAI_Gx3z_S_C1003_bb;
  Double I_NAI_G4y_Px_C1003_bb = I_NAI_Hx4y_S_C1003_bb+ABX*I_NAI_G4y_S_C1003_bb;
  Double I_NAI_G3yz_Px_C1003_bb = I_NAI_Hx3yz_S_C1003_bb+ABX*I_NAI_G3yz_S_C1003_bb;
  Double I_NAI_G2y2z_Px_C1003_bb = I_NAI_Hx2y2z_S_C1003_bb+ABX*I_NAI_G2y2z_S_C1003_bb;
  Double I_NAI_Gy3z_Px_C1003_bb = I_NAI_Hxy3z_S_C1003_bb+ABX*I_NAI_Gy3z_S_C1003_bb;
  Double I_NAI_G4z_Px_C1003_bb = I_NAI_Hx4z_S_C1003_bb+ABX*I_NAI_G4z_S_C1003_bb;
  Double I_NAI_G4x_Py_C1003_bb = I_NAI_H4xy_S_C1003_bb+ABY*I_NAI_G4x_S_C1003_bb;
  Double I_NAI_G3xy_Py_C1003_bb = I_NAI_H3x2y_S_C1003_bb+ABY*I_NAI_G3xy_S_C1003_bb;
  Double I_NAI_G3xz_Py_C1003_bb = I_NAI_H3xyz_S_C1003_bb+ABY*I_NAI_G3xz_S_C1003_bb;
  Double I_NAI_G2x2y_Py_C1003_bb = I_NAI_H2x3y_S_C1003_bb+ABY*I_NAI_G2x2y_S_C1003_bb;
  Double I_NAI_G2xyz_Py_C1003_bb = I_NAI_H2x2yz_S_C1003_bb+ABY*I_NAI_G2xyz_S_C1003_bb;
  Double I_NAI_G2x2z_Py_C1003_bb = I_NAI_H2xy2z_S_C1003_bb+ABY*I_NAI_G2x2z_S_C1003_bb;
  Double I_NAI_Gx3y_Py_C1003_bb = I_NAI_Hx4y_S_C1003_bb+ABY*I_NAI_Gx3y_S_C1003_bb;
  Double I_NAI_Gx2yz_Py_C1003_bb = I_NAI_Hx3yz_S_C1003_bb+ABY*I_NAI_Gx2yz_S_C1003_bb;
  Double I_NAI_Gxy2z_Py_C1003_bb = I_NAI_Hx2y2z_S_C1003_bb+ABY*I_NAI_Gxy2z_S_C1003_bb;
  Double I_NAI_Gx3z_Py_C1003_bb = I_NAI_Hxy3z_S_C1003_bb+ABY*I_NAI_Gx3z_S_C1003_bb;
  Double I_NAI_G4y_Py_C1003_bb = I_NAI_H5y_S_C1003_bb+ABY*I_NAI_G4y_S_C1003_bb;
  Double I_NAI_G3yz_Py_C1003_bb = I_NAI_H4yz_S_C1003_bb+ABY*I_NAI_G3yz_S_C1003_bb;
  Double I_NAI_G2y2z_Py_C1003_bb = I_NAI_H3y2z_S_C1003_bb+ABY*I_NAI_G2y2z_S_C1003_bb;
  Double I_NAI_Gy3z_Py_C1003_bb = I_NAI_H2y3z_S_C1003_bb+ABY*I_NAI_Gy3z_S_C1003_bb;
  Double I_NAI_G4z_Py_C1003_bb = I_NAI_Hy4z_S_C1003_bb+ABY*I_NAI_G4z_S_C1003_bb;
  Double I_NAI_G4x_Pz_C1003_bb = I_NAI_H4xz_S_C1003_bb+ABZ*I_NAI_G4x_S_C1003_bb;
  Double I_NAI_G3xy_Pz_C1003_bb = I_NAI_H3xyz_S_C1003_bb+ABZ*I_NAI_G3xy_S_C1003_bb;
  Double I_NAI_G3xz_Pz_C1003_bb = I_NAI_H3x2z_S_C1003_bb+ABZ*I_NAI_G3xz_S_C1003_bb;
  Double I_NAI_G2x2y_Pz_C1003_bb = I_NAI_H2x2yz_S_C1003_bb+ABZ*I_NAI_G2x2y_S_C1003_bb;
  Double I_NAI_G2xyz_Pz_C1003_bb = I_NAI_H2xy2z_S_C1003_bb+ABZ*I_NAI_G2xyz_S_C1003_bb;
  Double I_NAI_G2x2z_Pz_C1003_bb = I_NAI_H2x3z_S_C1003_bb+ABZ*I_NAI_G2x2z_S_C1003_bb;
  Double I_NAI_Gx3y_Pz_C1003_bb = I_NAI_Hx3yz_S_C1003_bb+ABZ*I_NAI_Gx3y_S_C1003_bb;
  Double I_NAI_Gx2yz_Pz_C1003_bb = I_NAI_Hx2y2z_S_C1003_bb+ABZ*I_NAI_Gx2yz_S_C1003_bb;
  Double I_NAI_Gxy2z_Pz_C1003_bb = I_NAI_Hxy3z_S_C1003_bb+ABZ*I_NAI_Gxy2z_S_C1003_bb;
  Double I_NAI_Gx3z_Pz_C1003_bb = I_NAI_Hx4z_S_C1003_bb+ABZ*I_NAI_Gx3z_S_C1003_bb;
  Double I_NAI_G4y_Pz_C1003_bb = I_NAI_H4yz_S_C1003_bb+ABZ*I_NAI_G4y_S_C1003_bb;
  Double I_NAI_G3yz_Pz_C1003_bb = I_NAI_H3y2z_S_C1003_bb+ABZ*I_NAI_G3yz_S_C1003_bb;
  Double I_NAI_G2y2z_Pz_C1003_bb = I_NAI_H2y3z_S_C1003_bb+ABZ*I_NAI_G2y2z_S_C1003_bb;
  Double I_NAI_Gy3z_Pz_C1003_bb = I_NAI_Hy4z_S_C1003_bb+ABZ*I_NAI_Gy3z_S_C1003_bb;
  Double I_NAI_G4z_Pz_C1003_bb = I_NAI_H5z_S_C1003_bb+ABZ*I_NAI_G4z_S_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_bb
   ************************************************************/
  Double I_NAI_F3x_D2x_C1003_bb = I_NAI_G4x_Px_C1003_bb+ABX*I_NAI_F3x_Px_C1003_bb;
  Double I_NAI_F2xy_D2x_C1003_bb = I_NAI_G3xy_Px_C1003_bb+ABX*I_NAI_F2xy_Px_C1003_bb;
  Double I_NAI_F2xz_D2x_C1003_bb = I_NAI_G3xz_Px_C1003_bb+ABX*I_NAI_F2xz_Px_C1003_bb;
  Double I_NAI_Fx2y_D2x_C1003_bb = I_NAI_G2x2y_Px_C1003_bb+ABX*I_NAI_Fx2y_Px_C1003_bb;
  Double I_NAI_Fxyz_D2x_C1003_bb = I_NAI_G2xyz_Px_C1003_bb+ABX*I_NAI_Fxyz_Px_C1003_bb;
  Double I_NAI_Fx2z_D2x_C1003_bb = I_NAI_G2x2z_Px_C1003_bb+ABX*I_NAI_Fx2z_Px_C1003_bb;
  Double I_NAI_F3y_D2x_C1003_bb = I_NAI_Gx3y_Px_C1003_bb+ABX*I_NAI_F3y_Px_C1003_bb;
  Double I_NAI_F2yz_D2x_C1003_bb = I_NAI_Gx2yz_Px_C1003_bb+ABX*I_NAI_F2yz_Px_C1003_bb;
  Double I_NAI_Fy2z_D2x_C1003_bb = I_NAI_Gxy2z_Px_C1003_bb+ABX*I_NAI_Fy2z_Px_C1003_bb;
  Double I_NAI_F3z_D2x_C1003_bb = I_NAI_Gx3z_Px_C1003_bb+ABX*I_NAI_F3z_Px_C1003_bb;
  Double I_NAI_F3x_Dxy_C1003_bb = I_NAI_G3xy_Px_C1003_bb+ABY*I_NAI_F3x_Px_C1003_bb;
  Double I_NAI_F2xy_Dxy_C1003_bb = I_NAI_G2x2y_Px_C1003_bb+ABY*I_NAI_F2xy_Px_C1003_bb;
  Double I_NAI_F2xz_Dxy_C1003_bb = I_NAI_G2xyz_Px_C1003_bb+ABY*I_NAI_F2xz_Px_C1003_bb;
  Double I_NAI_Fx2y_Dxy_C1003_bb = I_NAI_Gx3y_Px_C1003_bb+ABY*I_NAI_Fx2y_Px_C1003_bb;
  Double I_NAI_Fxyz_Dxy_C1003_bb = I_NAI_Gx2yz_Px_C1003_bb+ABY*I_NAI_Fxyz_Px_C1003_bb;
  Double I_NAI_Fx2z_Dxy_C1003_bb = I_NAI_Gxy2z_Px_C1003_bb+ABY*I_NAI_Fx2z_Px_C1003_bb;
  Double I_NAI_F3y_Dxy_C1003_bb = I_NAI_G4y_Px_C1003_bb+ABY*I_NAI_F3y_Px_C1003_bb;
  Double I_NAI_F2yz_Dxy_C1003_bb = I_NAI_G3yz_Px_C1003_bb+ABY*I_NAI_F2yz_Px_C1003_bb;
  Double I_NAI_Fy2z_Dxy_C1003_bb = I_NAI_G2y2z_Px_C1003_bb+ABY*I_NAI_Fy2z_Px_C1003_bb;
  Double I_NAI_F3z_Dxy_C1003_bb = I_NAI_Gy3z_Px_C1003_bb+ABY*I_NAI_F3z_Px_C1003_bb;
  Double I_NAI_F3x_D2y_C1003_bb = I_NAI_G3xy_Py_C1003_bb+ABY*I_NAI_F3x_Py_C1003_bb;
  Double I_NAI_F2xy_D2y_C1003_bb = I_NAI_G2x2y_Py_C1003_bb+ABY*I_NAI_F2xy_Py_C1003_bb;
  Double I_NAI_F2xz_D2y_C1003_bb = I_NAI_G2xyz_Py_C1003_bb+ABY*I_NAI_F2xz_Py_C1003_bb;
  Double I_NAI_Fx2y_D2y_C1003_bb = I_NAI_Gx3y_Py_C1003_bb+ABY*I_NAI_Fx2y_Py_C1003_bb;
  Double I_NAI_Fxyz_D2y_C1003_bb = I_NAI_Gx2yz_Py_C1003_bb+ABY*I_NAI_Fxyz_Py_C1003_bb;
  Double I_NAI_Fx2z_D2y_C1003_bb = I_NAI_Gxy2z_Py_C1003_bb+ABY*I_NAI_Fx2z_Py_C1003_bb;
  Double I_NAI_F3y_D2y_C1003_bb = I_NAI_G4y_Py_C1003_bb+ABY*I_NAI_F3y_Py_C1003_bb;
  Double I_NAI_F2yz_D2y_C1003_bb = I_NAI_G3yz_Py_C1003_bb+ABY*I_NAI_F2yz_Py_C1003_bb;
  Double I_NAI_Fy2z_D2y_C1003_bb = I_NAI_G2y2z_Py_C1003_bb+ABY*I_NAI_Fy2z_Py_C1003_bb;
  Double I_NAI_F3z_D2y_C1003_bb = I_NAI_Gy3z_Py_C1003_bb+ABY*I_NAI_F3z_Py_C1003_bb;
  Double I_NAI_F3x_D2z_C1003_bb = I_NAI_G3xz_Pz_C1003_bb+ABZ*I_NAI_F3x_Pz_C1003_bb;
  Double I_NAI_F2xy_D2z_C1003_bb = I_NAI_G2xyz_Pz_C1003_bb+ABZ*I_NAI_F2xy_Pz_C1003_bb;
  Double I_NAI_F2xz_D2z_C1003_bb = I_NAI_G2x2z_Pz_C1003_bb+ABZ*I_NAI_F2xz_Pz_C1003_bb;
  Double I_NAI_Fx2y_D2z_C1003_bb = I_NAI_Gx2yz_Pz_C1003_bb+ABZ*I_NAI_Fx2y_Pz_C1003_bb;
  Double I_NAI_Fxyz_D2z_C1003_bb = I_NAI_Gxy2z_Pz_C1003_bb+ABZ*I_NAI_Fxyz_Pz_C1003_bb;
  Double I_NAI_Fx2z_D2z_C1003_bb = I_NAI_Gx3z_Pz_C1003_bb+ABZ*I_NAI_Fx2z_Pz_C1003_bb;
  Double I_NAI_F3y_D2z_C1003_bb = I_NAI_G3yz_Pz_C1003_bb+ABZ*I_NAI_F3y_Pz_C1003_bb;
  Double I_NAI_F2yz_D2z_C1003_bb = I_NAI_G2y2z_Pz_C1003_bb+ABZ*I_NAI_F2yz_Pz_C1003_bb;
  Double I_NAI_Fy2z_D2z_C1003_bb = I_NAI_Gy3z_Pz_C1003_bb+ABZ*I_NAI_Fy2z_Pz_C1003_bb;
  Double I_NAI_F3z_D2z_C1003_bb = I_NAI_G4z_Pz_C1003_bb+ABZ*I_NAI_F3z_Pz_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 14 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_C1003_bb
   * RHS shell quartet name: SQ_NAI_H_S_C1003_bb
   ************************************************************/
  Double I_NAI_H5x_Px_C1003_bb = I_NAI_I6x_S_C1003_bb+ABX*I_NAI_H5x_S_C1003_bb;
  Double I_NAI_H4xy_Px_C1003_bb = I_NAI_I5xy_S_C1003_bb+ABX*I_NAI_H4xy_S_C1003_bb;
  Double I_NAI_H4xz_Px_C1003_bb = I_NAI_I5xz_S_C1003_bb+ABX*I_NAI_H4xz_S_C1003_bb;
  Double I_NAI_H3x2y_Px_C1003_bb = I_NAI_I4x2y_S_C1003_bb+ABX*I_NAI_H3x2y_S_C1003_bb;
  Double I_NAI_H3xyz_Px_C1003_bb = I_NAI_I4xyz_S_C1003_bb+ABX*I_NAI_H3xyz_S_C1003_bb;
  Double I_NAI_H3x2z_Px_C1003_bb = I_NAI_I4x2z_S_C1003_bb+ABX*I_NAI_H3x2z_S_C1003_bb;
  Double I_NAI_H2x3y_Px_C1003_bb = I_NAI_I3x3y_S_C1003_bb+ABX*I_NAI_H2x3y_S_C1003_bb;
  Double I_NAI_H2x2yz_Px_C1003_bb = I_NAI_I3x2yz_S_C1003_bb+ABX*I_NAI_H2x2yz_S_C1003_bb;
  Double I_NAI_H2xy2z_Px_C1003_bb = I_NAI_I3xy2z_S_C1003_bb+ABX*I_NAI_H2xy2z_S_C1003_bb;
  Double I_NAI_H2x3z_Px_C1003_bb = I_NAI_I3x3z_S_C1003_bb+ABX*I_NAI_H2x3z_S_C1003_bb;
  Double I_NAI_Hx4y_Px_C1003_bb = I_NAI_I2x4y_S_C1003_bb+ABX*I_NAI_Hx4y_S_C1003_bb;
  Double I_NAI_Hx3yz_Px_C1003_bb = I_NAI_I2x3yz_S_C1003_bb+ABX*I_NAI_Hx3yz_S_C1003_bb;
  Double I_NAI_Hx2y2z_Px_C1003_bb = I_NAI_I2x2y2z_S_C1003_bb+ABX*I_NAI_Hx2y2z_S_C1003_bb;
  Double I_NAI_Hxy3z_Px_C1003_bb = I_NAI_I2xy3z_S_C1003_bb+ABX*I_NAI_Hxy3z_S_C1003_bb;
  Double I_NAI_Hx4z_Px_C1003_bb = I_NAI_I2x4z_S_C1003_bb+ABX*I_NAI_Hx4z_S_C1003_bb;
  Double I_NAI_H4yz_Px_C1003_bb = I_NAI_Ix4yz_S_C1003_bb+ABX*I_NAI_H4yz_S_C1003_bb;
  Double I_NAI_H3y2z_Px_C1003_bb = I_NAI_Ix3y2z_S_C1003_bb+ABX*I_NAI_H3y2z_S_C1003_bb;
  Double I_NAI_H2y3z_Px_C1003_bb = I_NAI_Ix2y3z_S_C1003_bb+ABX*I_NAI_H2y3z_S_C1003_bb;
  Double I_NAI_Hy4z_Px_C1003_bb = I_NAI_Ixy4z_S_C1003_bb+ABX*I_NAI_Hy4z_S_C1003_bb;
  Double I_NAI_H4xy_Py_C1003_bb = I_NAI_I4x2y_S_C1003_bb+ABY*I_NAI_H4xy_S_C1003_bb;
  Double I_NAI_H3x2y_Py_C1003_bb = I_NAI_I3x3y_S_C1003_bb+ABY*I_NAI_H3x2y_S_C1003_bb;
  Double I_NAI_H3xyz_Py_C1003_bb = I_NAI_I3x2yz_S_C1003_bb+ABY*I_NAI_H3xyz_S_C1003_bb;
  Double I_NAI_H2x3y_Py_C1003_bb = I_NAI_I2x4y_S_C1003_bb+ABY*I_NAI_H2x3y_S_C1003_bb;
  Double I_NAI_H2x2yz_Py_C1003_bb = I_NAI_I2x3yz_S_C1003_bb+ABY*I_NAI_H2x2yz_S_C1003_bb;
  Double I_NAI_H2xy2z_Py_C1003_bb = I_NAI_I2x2y2z_S_C1003_bb+ABY*I_NAI_H2xy2z_S_C1003_bb;
  Double I_NAI_Hx4y_Py_C1003_bb = I_NAI_Ix5y_S_C1003_bb+ABY*I_NAI_Hx4y_S_C1003_bb;
  Double I_NAI_Hx3yz_Py_C1003_bb = I_NAI_Ix4yz_S_C1003_bb+ABY*I_NAI_Hx3yz_S_C1003_bb;
  Double I_NAI_Hx2y2z_Py_C1003_bb = I_NAI_Ix3y2z_S_C1003_bb+ABY*I_NAI_Hx2y2z_S_C1003_bb;
  Double I_NAI_Hxy3z_Py_C1003_bb = I_NAI_Ix2y3z_S_C1003_bb+ABY*I_NAI_Hxy3z_S_C1003_bb;
  Double I_NAI_H5y_Py_C1003_bb = I_NAI_I6y_S_C1003_bb+ABY*I_NAI_H5y_S_C1003_bb;
  Double I_NAI_H4yz_Py_C1003_bb = I_NAI_I5yz_S_C1003_bb+ABY*I_NAI_H4yz_S_C1003_bb;
  Double I_NAI_H3y2z_Py_C1003_bb = I_NAI_I4y2z_S_C1003_bb+ABY*I_NAI_H3y2z_S_C1003_bb;
  Double I_NAI_H2y3z_Py_C1003_bb = I_NAI_I3y3z_S_C1003_bb+ABY*I_NAI_H2y3z_S_C1003_bb;
  Double I_NAI_Hy4z_Py_C1003_bb = I_NAI_I2y4z_S_C1003_bb+ABY*I_NAI_Hy4z_S_C1003_bb;
  Double I_NAI_H4xz_Pz_C1003_bb = I_NAI_I4x2z_S_C1003_bb+ABZ*I_NAI_H4xz_S_C1003_bb;
  Double I_NAI_H3xyz_Pz_C1003_bb = I_NAI_I3xy2z_S_C1003_bb+ABZ*I_NAI_H3xyz_S_C1003_bb;
  Double I_NAI_H3x2z_Pz_C1003_bb = I_NAI_I3x3z_S_C1003_bb+ABZ*I_NAI_H3x2z_S_C1003_bb;
  Double I_NAI_H2x2yz_Pz_C1003_bb = I_NAI_I2x2y2z_S_C1003_bb+ABZ*I_NAI_H2x2yz_S_C1003_bb;
  Double I_NAI_H2xy2z_Pz_C1003_bb = I_NAI_I2xy3z_S_C1003_bb+ABZ*I_NAI_H2xy2z_S_C1003_bb;
  Double I_NAI_H2x3z_Pz_C1003_bb = I_NAI_I2x4z_S_C1003_bb+ABZ*I_NAI_H2x3z_S_C1003_bb;
  Double I_NAI_Hx3yz_Pz_C1003_bb = I_NAI_Ix3y2z_S_C1003_bb+ABZ*I_NAI_Hx3yz_S_C1003_bb;
  Double I_NAI_Hx2y2z_Pz_C1003_bb = I_NAI_Ix2y3z_S_C1003_bb+ABZ*I_NAI_Hx2y2z_S_C1003_bb;
  Double I_NAI_Hxy3z_Pz_C1003_bb = I_NAI_Ixy4z_S_C1003_bb+ABZ*I_NAI_Hxy3z_S_C1003_bb;
  Double I_NAI_Hx4z_Pz_C1003_bb = I_NAI_Ix5z_S_C1003_bb+ABZ*I_NAI_Hx4z_S_C1003_bb;
  Double I_NAI_H4yz_Pz_C1003_bb = I_NAI_I4y2z_S_C1003_bb+ABZ*I_NAI_H4yz_S_C1003_bb;
  Double I_NAI_H3y2z_Pz_C1003_bb = I_NAI_I3y3z_S_C1003_bb+ABZ*I_NAI_H3y2z_S_C1003_bb;
  Double I_NAI_H2y3z_Pz_C1003_bb = I_NAI_I2y4z_S_C1003_bb+ABZ*I_NAI_H2y3z_S_C1003_bb;
  Double I_NAI_Hy4z_Pz_C1003_bb = I_NAI_Iy5z_S_C1003_bb+ABZ*I_NAI_Hy4z_S_C1003_bb;
  Double I_NAI_H5z_Pz_C1003_bb = I_NAI_I6z_S_C1003_bb+ABZ*I_NAI_H5z_S_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 35 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_bb
   * RHS shell quartet name: SQ_NAI_G_P_C1003_bb
   ************************************************************/
  Double I_NAI_G4x_D2x_C1003_bb = I_NAI_H5x_Px_C1003_bb+ABX*I_NAI_G4x_Px_C1003_bb;
  Double I_NAI_G3xy_D2x_C1003_bb = I_NAI_H4xy_Px_C1003_bb+ABX*I_NAI_G3xy_Px_C1003_bb;
  Double I_NAI_G3xz_D2x_C1003_bb = I_NAI_H4xz_Px_C1003_bb+ABX*I_NAI_G3xz_Px_C1003_bb;
  Double I_NAI_G2x2y_D2x_C1003_bb = I_NAI_H3x2y_Px_C1003_bb+ABX*I_NAI_G2x2y_Px_C1003_bb;
  Double I_NAI_G2xyz_D2x_C1003_bb = I_NAI_H3xyz_Px_C1003_bb+ABX*I_NAI_G2xyz_Px_C1003_bb;
  Double I_NAI_G2x2z_D2x_C1003_bb = I_NAI_H3x2z_Px_C1003_bb+ABX*I_NAI_G2x2z_Px_C1003_bb;
  Double I_NAI_Gx3y_D2x_C1003_bb = I_NAI_H2x3y_Px_C1003_bb+ABX*I_NAI_Gx3y_Px_C1003_bb;
  Double I_NAI_Gx2yz_D2x_C1003_bb = I_NAI_H2x2yz_Px_C1003_bb+ABX*I_NAI_Gx2yz_Px_C1003_bb;
  Double I_NAI_Gxy2z_D2x_C1003_bb = I_NAI_H2xy2z_Px_C1003_bb+ABX*I_NAI_Gxy2z_Px_C1003_bb;
  Double I_NAI_Gx3z_D2x_C1003_bb = I_NAI_H2x3z_Px_C1003_bb+ABX*I_NAI_Gx3z_Px_C1003_bb;
  Double I_NAI_G4y_D2x_C1003_bb = I_NAI_Hx4y_Px_C1003_bb+ABX*I_NAI_G4y_Px_C1003_bb;
  Double I_NAI_G3yz_D2x_C1003_bb = I_NAI_Hx3yz_Px_C1003_bb+ABX*I_NAI_G3yz_Px_C1003_bb;
  Double I_NAI_G2y2z_D2x_C1003_bb = I_NAI_Hx2y2z_Px_C1003_bb+ABX*I_NAI_G2y2z_Px_C1003_bb;
  Double I_NAI_Gy3z_D2x_C1003_bb = I_NAI_Hxy3z_Px_C1003_bb+ABX*I_NAI_Gy3z_Px_C1003_bb;
  Double I_NAI_G4z_D2x_C1003_bb = I_NAI_Hx4z_Px_C1003_bb+ABX*I_NAI_G4z_Px_C1003_bb;
  Double I_NAI_G3xz_Dxy_C1003_bb = I_NAI_H3xyz_Px_C1003_bb+ABY*I_NAI_G3xz_Px_C1003_bb;
  Double I_NAI_G2xyz_Dxy_C1003_bb = I_NAI_H2x2yz_Px_C1003_bb+ABY*I_NAI_G2xyz_Px_C1003_bb;
  Double I_NAI_G2x2z_Dxy_C1003_bb = I_NAI_H2xy2z_Px_C1003_bb+ABY*I_NAI_G2x2z_Px_C1003_bb;
  Double I_NAI_Gx2yz_Dxy_C1003_bb = I_NAI_Hx3yz_Px_C1003_bb+ABY*I_NAI_Gx2yz_Px_C1003_bb;
  Double I_NAI_Gxy2z_Dxy_C1003_bb = I_NAI_Hx2y2z_Px_C1003_bb+ABY*I_NAI_Gxy2z_Px_C1003_bb;
  Double I_NAI_Gx3z_Dxy_C1003_bb = I_NAI_Hxy3z_Px_C1003_bb+ABY*I_NAI_Gx3z_Px_C1003_bb;
  Double I_NAI_G3yz_Dxy_C1003_bb = I_NAI_H4yz_Px_C1003_bb+ABY*I_NAI_G3yz_Px_C1003_bb;
  Double I_NAI_G2y2z_Dxy_C1003_bb = I_NAI_H3y2z_Px_C1003_bb+ABY*I_NAI_G2y2z_Px_C1003_bb;
  Double I_NAI_Gy3z_Dxy_C1003_bb = I_NAI_H2y3z_Px_C1003_bb+ABY*I_NAI_Gy3z_Px_C1003_bb;
  Double I_NAI_G4z_Dxy_C1003_bb = I_NAI_Hy4z_Px_C1003_bb+ABY*I_NAI_G4z_Px_C1003_bb;
  Double I_NAI_G4x_D2y_C1003_bb = I_NAI_H4xy_Py_C1003_bb+ABY*I_NAI_G4x_Py_C1003_bb;
  Double I_NAI_G3xy_D2y_C1003_bb = I_NAI_H3x2y_Py_C1003_bb+ABY*I_NAI_G3xy_Py_C1003_bb;
  Double I_NAI_G3xz_D2y_C1003_bb = I_NAI_H3xyz_Py_C1003_bb+ABY*I_NAI_G3xz_Py_C1003_bb;
  Double I_NAI_G2x2y_D2y_C1003_bb = I_NAI_H2x3y_Py_C1003_bb+ABY*I_NAI_G2x2y_Py_C1003_bb;
  Double I_NAI_G2xyz_D2y_C1003_bb = I_NAI_H2x2yz_Py_C1003_bb+ABY*I_NAI_G2xyz_Py_C1003_bb;
  Double I_NAI_G2x2z_D2y_C1003_bb = I_NAI_H2xy2z_Py_C1003_bb+ABY*I_NAI_G2x2z_Py_C1003_bb;
  Double I_NAI_Gx3y_D2y_C1003_bb = I_NAI_Hx4y_Py_C1003_bb+ABY*I_NAI_Gx3y_Py_C1003_bb;
  Double I_NAI_Gx2yz_D2y_C1003_bb = I_NAI_Hx3yz_Py_C1003_bb+ABY*I_NAI_Gx2yz_Py_C1003_bb;
  Double I_NAI_Gxy2z_D2y_C1003_bb = I_NAI_Hx2y2z_Py_C1003_bb+ABY*I_NAI_Gxy2z_Py_C1003_bb;
  Double I_NAI_Gx3z_D2y_C1003_bb = I_NAI_Hxy3z_Py_C1003_bb+ABY*I_NAI_Gx3z_Py_C1003_bb;
  Double I_NAI_G4y_D2y_C1003_bb = I_NAI_H5y_Py_C1003_bb+ABY*I_NAI_G4y_Py_C1003_bb;
  Double I_NAI_G3yz_D2y_C1003_bb = I_NAI_H4yz_Py_C1003_bb+ABY*I_NAI_G3yz_Py_C1003_bb;
  Double I_NAI_G2y2z_D2y_C1003_bb = I_NAI_H3y2z_Py_C1003_bb+ABY*I_NAI_G2y2z_Py_C1003_bb;
  Double I_NAI_Gy3z_D2y_C1003_bb = I_NAI_H2y3z_Py_C1003_bb+ABY*I_NAI_Gy3z_Py_C1003_bb;
  Double I_NAI_G4z_D2y_C1003_bb = I_NAI_Hy4z_Py_C1003_bb+ABY*I_NAI_G4z_Py_C1003_bb;
  Double I_NAI_G4x_D2z_C1003_bb = I_NAI_H4xz_Pz_C1003_bb+ABZ*I_NAI_G4x_Pz_C1003_bb;
  Double I_NAI_G3xy_D2z_C1003_bb = I_NAI_H3xyz_Pz_C1003_bb+ABZ*I_NAI_G3xy_Pz_C1003_bb;
  Double I_NAI_G3xz_D2z_C1003_bb = I_NAI_H3x2z_Pz_C1003_bb+ABZ*I_NAI_G3xz_Pz_C1003_bb;
  Double I_NAI_G2x2y_D2z_C1003_bb = I_NAI_H2x2yz_Pz_C1003_bb+ABZ*I_NAI_G2x2y_Pz_C1003_bb;
  Double I_NAI_G2xyz_D2z_C1003_bb = I_NAI_H2xy2z_Pz_C1003_bb+ABZ*I_NAI_G2xyz_Pz_C1003_bb;
  Double I_NAI_G2x2z_D2z_C1003_bb = I_NAI_H2x3z_Pz_C1003_bb+ABZ*I_NAI_G2x2z_Pz_C1003_bb;
  Double I_NAI_Gx3y_D2z_C1003_bb = I_NAI_Hx3yz_Pz_C1003_bb+ABZ*I_NAI_Gx3y_Pz_C1003_bb;
  Double I_NAI_Gx2yz_D2z_C1003_bb = I_NAI_Hx2y2z_Pz_C1003_bb+ABZ*I_NAI_Gx2yz_Pz_C1003_bb;
  Double I_NAI_Gxy2z_D2z_C1003_bb = I_NAI_Hxy3z_Pz_C1003_bb+ABZ*I_NAI_Gxy2z_Pz_C1003_bb;
  Double I_NAI_Gx3z_D2z_C1003_bb = I_NAI_Hx4z_Pz_C1003_bb+ABZ*I_NAI_Gx3z_Pz_C1003_bb;
  Double I_NAI_G4y_D2z_C1003_bb = I_NAI_H4yz_Pz_C1003_bb+ABZ*I_NAI_G4y_Pz_C1003_bb;
  Double I_NAI_G3yz_D2z_C1003_bb = I_NAI_H3y2z_Pz_C1003_bb+ABZ*I_NAI_G3yz_Pz_C1003_bb;
  Double I_NAI_G2y2z_D2z_C1003_bb = I_NAI_H2y3z_Pz_C1003_bb+ABZ*I_NAI_G2y2z_Pz_C1003_bb;
  Double I_NAI_Gy3z_D2z_C1003_bb = I_NAI_Hy4z_Pz_C1003_bb+ABZ*I_NAI_Gy3z_Pz_C1003_bb;
  Double I_NAI_G4z_D2z_C1003_bb = I_NAI_H5z_Pz_C1003_bb+ABZ*I_NAI_G4z_Pz_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_C1003_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_D_C1003_bb
   ************************************************************/
  Double I_NAI_F3x_F3x_C1003_bb = I_NAI_G4x_D2x_C1003_bb+ABX*I_NAI_F3x_D2x_C1003_bb;
  Double I_NAI_F2xy_F3x_C1003_bb = I_NAI_G3xy_D2x_C1003_bb+ABX*I_NAI_F2xy_D2x_C1003_bb;
  Double I_NAI_F2xz_F3x_C1003_bb = I_NAI_G3xz_D2x_C1003_bb+ABX*I_NAI_F2xz_D2x_C1003_bb;
  Double I_NAI_Fx2y_F3x_C1003_bb = I_NAI_G2x2y_D2x_C1003_bb+ABX*I_NAI_Fx2y_D2x_C1003_bb;
  Double I_NAI_Fxyz_F3x_C1003_bb = I_NAI_G2xyz_D2x_C1003_bb+ABX*I_NAI_Fxyz_D2x_C1003_bb;
  Double I_NAI_Fx2z_F3x_C1003_bb = I_NAI_G2x2z_D2x_C1003_bb+ABX*I_NAI_Fx2z_D2x_C1003_bb;
  Double I_NAI_F3y_F3x_C1003_bb = I_NAI_Gx3y_D2x_C1003_bb+ABX*I_NAI_F3y_D2x_C1003_bb;
  Double I_NAI_F2yz_F3x_C1003_bb = I_NAI_Gx2yz_D2x_C1003_bb+ABX*I_NAI_F2yz_D2x_C1003_bb;
  Double I_NAI_Fy2z_F3x_C1003_bb = I_NAI_Gxy2z_D2x_C1003_bb+ABX*I_NAI_Fy2z_D2x_C1003_bb;
  Double I_NAI_F3z_F3x_C1003_bb = I_NAI_Gx3z_D2x_C1003_bb+ABX*I_NAI_F3z_D2x_C1003_bb;
  Double I_NAI_F3x_F2xy_C1003_bb = I_NAI_G3xy_D2x_C1003_bb+ABY*I_NAI_F3x_D2x_C1003_bb;
  Double I_NAI_F2xy_F2xy_C1003_bb = I_NAI_G2x2y_D2x_C1003_bb+ABY*I_NAI_F2xy_D2x_C1003_bb;
  Double I_NAI_F2xz_F2xy_C1003_bb = I_NAI_G2xyz_D2x_C1003_bb+ABY*I_NAI_F2xz_D2x_C1003_bb;
  Double I_NAI_Fx2y_F2xy_C1003_bb = I_NAI_Gx3y_D2x_C1003_bb+ABY*I_NAI_Fx2y_D2x_C1003_bb;
  Double I_NAI_Fxyz_F2xy_C1003_bb = I_NAI_Gx2yz_D2x_C1003_bb+ABY*I_NAI_Fxyz_D2x_C1003_bb;
  Double I_NAI_Fx2z_F2xy_C1003_bb = I_NAI_Gxy2z_D2x_C1003_bb+ABY*I_NAI_Fx2z_D2x_C1003_bb;
  Double I_NAI_F3y_F2xy_C1003_bb = I_NAI_G4y_D2x_C1003_bb+ABY*I_NAI_F3y_D2x_C1003_bb;
  Double I_NAI_F2yz_F2xy_C1003_bb = I_NAI_G3yz_D2x_C1003_bb+ABY*I_NAI_F2yz_D2x_C1003_bb;
  Double I_NAI_Fy2z_F2xy_C1003_bb = I_NAI_G2y2z_D2x_C1003_bb+ABY*I_NAI_Fy2z_D2x_C1003_bb;
  Double I_NAI_F3z_F2xy_C1003_bb = I_NAI_Gy3z_D2x_C1003_bb+ABY*I_NAI_F3z_D2x_C1003_bb;
  Double I_NAI_F3x_F2xz_C1003_bb = I_NAI_G3xz_D2x_C1003_bb+ABZ*I_NAI_F3x_D2x_C1003_bb;
  Double I_NAI_F2xy_F2xz_C1003_bb = I_NAI_G2xyz_D2x_C1003_bb+ABZ*I_NAI_F2xy_D2x_C1003_bb;
  Double I_NAI_F2xz_F2xz_C1003_bb = I_NAI_G2x2z_D2x_C1003_bb+ABZ*I_NAI_F2xz_D2x_C1003_bb;
  Double I_NAI_Fx2y_F2xz_C1003_bb = I_NAI_Gx2yz_D2x_C1003_bb+ABZ*I_NAI_Fx2y_D2x_C1003_bb;
  Double I_NAI_Fxyz_F2xz_C1003_bb = I_NAI_Gxy2z_D2x_C1003_bb+ABZ*I_NAI_Fxyz_D2x_C1003_bb;
  Double I_NAI_Fx2z_F2xz_C1003_bb = I_NAI_Gx3z_D2x_C1003_bb+ABZ*I_NAI_Fx2z_D2x_C1003_bb;
  Double I_NAI_F3y_F2xz_C1003_bb = I_NAI_G3yz_D2x_C1003_bb+ABZ*I_NAI_F3y_D2x_C1003_bb;
  Double I_NAI_F2yz_F2xz_C1003_bb = I_NAI_G2y2z_D2x_C1003_bb+ABZ*I_NAI_F2yz_D2x_C1003_bb;
  Double I_NAI_Fy2z_F2xz_C1003_bb = I_NAI_Gy3z_D2x_C1003_bb+ABZ*I_NAI_Fy2z_D2x_C1003_bb;
  Double I_NAI_F3z_F2xz_C1003_bb = I_NAI_G4z_D2x_C1003_bb+ABZ*I_NAI_F3z_D2x_C1003_bb;
  Double I_NAI_F3x_Fx2y_C1003_bb = I_NAI_G4x_D2y_C1003_bb+ABX*I_NAI_F3x_D2y_C1003_bb;
  Double I_NAI_F2xy_Fx2y_C1003_bb = I_NAI_G3xy_D2y_C1003_bb+ABX*I_NAI_F2xy_D2y_C1003_bb;
  Double I_NAI_F2xz_Fx2y_C1003_bb = I_NAI_G3xz_D2y_C1003_bb+ABX*I_NAI_F2xz_D2y_C1003_bb;
  Double I_NAI_Fx2y_Fx2y_C1003_bb = I_NAI_G2x2y_D2y_C1003_bb+ABX*I_NAI_Fx2y_D2y_C1003_bb;
  Double I_NAI_Fxyz_Fx2y_C1003_bb = I_NAI_G2xyz_D2y_C1003_bb+ABX*I_NAI_Fxyz_D2y_C1003_bb;
  Double I_NAI_Fx2z_Fx2y_C1003_bb = I_NAI_G2x2z_D2y_C1003_bb+ABX*I_NAI_Fx2z_D2y_C1003_bb;
  Double I_NAI_F3y_Fx2y_C1003_bb = I_NAI_Gx3y_D2y_C1003_bb+ABX*I_NAI_F3y_D2y_C1003_bb;
  Double I_NAI_F2yz_Fx2y_C1003_bb = I_NAI_Gx2yz_D2y_C1003_bb+ABX*I_NAI_F2yz_D2y_C1003_bb;
  Double I_NAI_Fy2z_Fx2y_C1003_bb = I_NAI_Gxy2z_D2y_C1003_bb+ABX*I_NAI_Fy2z_D2y_C1003_bb;
  Double I_NAI_F3z_Fx2y_C1003_bb = I_NAI_Gx3z_D2y_C1003_bb+ABX*I_NAI_F3z_D2y_C1003_bb;
  Double I_NAI_F3x_Fxyz_C1003_bb = I_NAI_G3xz_Dxy_C1003_bb+ABZ*I_NAI_F3x_Dxy_C1003_bb;
  Double I_NAI_F2xy_Fxyz_C1003_bb = I_NAI_G2xyz_Dxy_C1003_bb+ABZ*I_NAI_F2xy_Dxy_C1003_bb;
  Double I_NAI_F2xz_Fxyz_C1003_bb = I_NAI_G2x2z_Dxy_C1003_bb+ABZ*I_NAI_F2xz_Dxy_C1003_bb;
  Double I_NAI_Fx2y_Fxyz_C1003_bb = I_NAI_Gx2yz_Dxy_C1003_bb+ABZ*I_NAI_Fx2y_Dxy_C1003_bb;
  Double I_NAI_Fxyz_Fxyz_C1003_bb = I_NAI_Gxy2z_Dxy_C1003_bb+ABZ*I_NAI_Fxyz_Dxy_C1003_bb;
  Double I_NAI_Fx2z_Fxyz_C1003_bb = I_NAI_Gx3z_Dxy_C1003_bb+ABZ*I_NAI_Fx2z_Dxy_C1003_bb;
  Double I_NAI_F3y_Fxyz_C1003_bb = I_NAI_G3yz_Dxy_C1003_bb+ABZ*I_NAI_F3y_Dxy_C1003_bb;
  Double I_NAI_F2yz_Fxyz_C1003_bb = I_NAI_G2y2z_Dxy_C1003_bb+ABZ*I_NAI_F2yz_Dxy_C1003_bb;
  Double I_NAI_Fy2z_Fxyz_C1003_bb = I_NAI_Gy3z_Dxy_C1003_bb+ABZ*I_NAI_Fy2z_Dxy_C1003_bb;
  Double I_NAI_F3z_Fxyz_C1003_bb = I_NAI_G4z_Dxy_C1003_bb+ABZ*I_NAI_F3z_Dxy_C1003_bb;
  Double I_NAI_F3x_Fx2z_C1003_bb = I_NAI_G4x_D2z_C1003_bb+ABX*I_NAI_F3x_D2z_C1003_bb;
  Double I_NAI_F2xy_Fx2z_C1003_bb = I_NAI_G3xy_D2z_C1003_bb+ABX*I_NAI_F2xy_D2z_C1003_bb;
  Double I_NAI_F2xz_Fx2z_C1003_bb = I_NAI_G3xz_D2z_C1003_bb+ABX*I_NAI_F2xz_D2z_C1003_bb;
  Double I_NAI_Fx2y_Fx2z_C1003_bb = I_NAI_G2x2y_D2z_C1003_bb+ABX*I_NAI_Fx2y_D2z_C1003_bb;
  Double I_NAI_Fxyz_Fx2z_C1003_bb = I_NAI_G2xyz_D2z_C1003_bb+ABX*I_NAI_Fxyz_D2z_C1003_bb;
  Double I_NAI_Fx2z_Fx2z_C1003_bb = I_NAI_G2x2z_D2z_C1003_bb+ABX*I_NAI_Fx2z_D2z_C1003_bb;
  Double I_NAI_F3y_Fx2z_C1003_bb = I_NAI_Gx3y_D2z_C1003_bb+ABX*I_NAI_F3y_D2z_C1003_bb;
  Double I_NAI_F2yz_Fx2z_C1003_bb = I_NAI_Gx2yz_D2z_C1003_bb+ABX*I_NAI_F2yz_D2z_C1003_bb;
  Double I_NAI_Fy2z_Fx2z_C1003_bb = I_NAI_Gxy2z_D2z_C1003_bb+ABX*I_NAI_Fy2z_D2z_C1003_bb;
  Double I_NAI_F3z_Fx2z_C1003_bb = I_NAI_Gx3z_D2z_C1003_bb+ABX*I_NAI_F3z_D2z_C1003_bb;
  Double I_NAI_F3x_F3y_C1003_bb = I_NAI_G3xy_D2y_C1003_bb+ABY*I_NAI_F3x_D2y_C1003_bb;
  Double I_NAI_F2xy_F3y_C1003_bb = I_NAI_G2x2y_D2y_C1003_bb+ABY*I_NAI_F2xy_D2y_C1003_bb;
  Double I_NAI_F2xz_F3y_C1003_bb = I_NAI_G2xyz_D2y_C1003_bb+ABY*I_NAI_F2xz_D2y_C1003_bb;
  Double I_NAI_Fx2y_F3y_C1003_bb = I_NAI_Gx3y_D2y_C1003_bb+ABY*I_NAI_Fx2y_D2y_C1003_bb;
  Double I_NAI_Fxyz_F3y_C1003_bb = I_NAI_Gx2yz_D2y_C1003_bb+ABY*I_NAI_Fxyz_D2y_C1003_bb;
  Double I_NAI_Fx2z_F3y_C1003_bb = I_NAI_Gxy2z_D2y_C1003_bb+ABY*I_NAI_Fx2z_D2y_C1003_bb;
  Double I_NAI_F3y_F3y_C1003_bb = I_NAI_G4y_D2y_C1003_bb+ABY*I_NAI_F3y_D2y_C1003_bb;
  Double I_NAI_F2yz_F3y_C1003_bb = I_NAI_G3yz_D2y_C1003_bb+ABY*I_NAI_F2yz_D2y_C1003_bb;
  Double I_NAI_Fy2z_F3y_C1003_bb = I_NAI_G2y2z_D2y_C1003_bb+ABY*I_NAI_Fy2z_D2y_C1003_bb;
  Double I_NAI_F3z_F3y_C1003_bb = I_NAI_Gy3z_D2y_C1003_bb+ABY*I_NAI_F3z_D2y_C1003_bb;
  Double I_NAI_F3x_F2yz_C1003_bb = I_NAI_G3xz_D2y_C1003_bb+ABZ*I_NAI_F3x_D2y_C1003_bb;
  Double I_NAI_F2xy_F2yz_C1003_bb = I_NAI_G2xyz_D2y_C1003_bb+ABZ*I_NAI_F2xy_D2y_C1003_bb;
  Double I_NAI_F2xz_F2yz_C1003_bb = I_NAI_G2x2z_D2y_C1003_bb+ABZ*I_NAI_F2xz_D2y_C1003_bb;
  Double I_NAI_Fx2y_F2yz_C1003_bb = I_NAI_Gx2yz_D2y_C1003_bb+ABZ*I_NAI_Fx2y_D2y_C1003_bb;
  Double I_NAI_Fxyz_F2yz_C1003_bb = I_NAI_Gxy2z_D2y_C1003_bb+ABZ*I_NAI_Fxyz_D2y_C1003_bb;
  Double I_NAI_Fx2z_F2yz_C1003_bb = I_NAI_Gx3z_D2y_C1003_bb+ABZ*I_NAI_Fx2z_D2y_C1003_bb;
  Double I_NAI_F3y_F2yz_C1003_bb = I_NAI_G3yz_D2y_C1003_bb+ABZ*I_NAI_F3y_D2y_C1003_bb;
  Double I_NAI_F2yz_F2yz_C1003_bb = I_NAI_G2y2z_D2y_C1003_bb+ABZ*I_NAI_F2yz_D2y_C1003_bb;
  Double I_NAI_Fy2z_F2yz_C1003_bb = I_NAI_Gy3z_D2y_C1003_bb+ABZ*I_NAI_Fy2z_D2y_C1003_bb;
  Double I_NAI_F3z_F2yz_C1003_bb = I_NAI_G4z_D2y_C1003_bb+ABZ*I_NAI_F3z_D2y_C1003_bb;
  Double I_NAI_F3x_Fy2z_C1003_bb = I_NAI_G3xy_D2z_C1003_bb+ABY*I_NAI_F3x_D2z_C1003_bb;
  Double I_NAI_F2xy_Fy2z_C1003_bb = I_NAI_G2x2y_D2z_C1003_bb+ABY*I_NAI_F2xy_D2z_C1003_bb;
  Double I_NAI_F2xz_Fy2z_C1003_bb = I_NAI_G2xyz_D2z_C1003_bb+ABY*I_NAI_F2xz_D2z_C1003_bb;
  Double I_NAI_Fx2y_Fy2z_C1003_bb = I_NAI_Gx3y_D2z_C1003_bb+ABY*I_NAI_Fx2y_D2z_C1003_bb;
  Double I_NAI_Fxyz_Fy2z_C1003_bb = I_NAI_Gx2yz_D2z_C1003_bb+ABY*I_NAI_Fxyz_D2z_C1003_bb;
  Double I_NAI_Fx2z_Fy2z_C1003_bb = I_NAI_Gxy2z_D2z_C1003_bb+ABY*I_NAI_Fx2z_D2z_C1003_bb;
  Double I_NAI_F3y_Fy2z_C1003_bb = I_NAI_G4y_D2z_C1003_bb+ABY*I_NAI_F3y_D2z_C1003_bb;
  Double I_NAI_F2yz_Fy2z_C1003_bb = I_NAI_G3yz_D2z_C1003_bb+ABY*I_NAI_F2yz_D2z_C1003_bb;
  Double I_NAI_Fy2z_Fy2z_C1003_bb = I_NAI_G2y2z_D2z_C1003_bb+ABY*I_NAI_Fy2z_D2z_C1003_bb;
  Double I_NAI_F3z_Fy2z_C1003_bb = I_NAI_Gy3z_D2z_C1003_bb+ABY*I_NAI_F3z_D2z_C1003_bb;
  Double I_NAI_F3x_F3z_C1003_bb = I_NAI_G3xz_D2z_C1003_bb+ABZ*I_NAI_F3x_D2z_C1003_bb;
  Double I_NAI_F2xy_F3z_C1003_bb = I_NAI_G2xyz_D2z_C1003_bb+ABZ*I_NAI_F2xy_D2z_C1003_bb;
  Double I_NAI_F2xz_F3z_C1003_bb = I_NAI_G2x2z_D2z_C1003_bb+ABZ*I_NAI_F2xz_D2z_C1003_bb;
  Double I_NAI_Fx2y_F3z_C1003_bb = I_NAI_Gx2yz_D2z_C1003_bb+ABZ*I_NAI_Fx2y_D2z_C1003_bb;
  Double I_NAI_Fxyz_F3z_C1003_bb = I_NAI_Gxy2z_D2z_C1003_bb+ABZ*I_NAI_Fxyz_D2z_C1003_bb;
  Double I_NAI_Fx2z_F3z_C1003_bb = I_NAI_Gx3z_D2z_C1003_bb+ABZ*I_NAI_Fx2z_D2z_C1003_bb;
  Double I_NAI_F3y_F3z_C1003_bb = I_NAI_G3yz_D2z_C1003_bb+ABZ*I_NAI_F3y_D2z_C1003_bb;
  Double I_NAI_F2yz_F3z_C1003_bb = I_NAI_G2y2z_D2z_C1003_bb+ABZ*I_NAI_F2yz_D2z_C1003_bb;
  Double I_NAI_Fy2z_F3z_C1003_bb = I_NAI_Gy3z_D2z_C1003_bb+ABZ*I_NAI_Fy2z_D2z_C1003_bb;
  Double I_NAI_F3z_F3z_C1003_bb = I_NAI_G4z_D2z_C1003_bb+ABZ*I_NAI_F3z_D2z_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[0] = 4.0E0*I_NAI_H5x_S_C3_aa-2.0E0*3*I_NAI_F3x_S_C3_a-2.0E0*4*I_NAI_F3x_S_C3_a+3*2*I_NAI_Px_S_C3;
  abcd[1] = 4.0E0*I_NAI_H4xy_S_C3_aa-2.0E0*2*I_NAI_F2xy_S_C3_a-2.0E0*3*I_NAI_F2xy_S_C3_a+2*1*I_NAI_Py_S_C3;
  abcd[2] = 4.0E0*I_NAI_H4xz_S_C3_aa-2.0E0*2*I_NAI_F2xz_S_C3_a-2.0E0*3*I_NAI_F2xz_S_C3_a+2*1*I_NAI_Pz_S_C3;
  abcd[3] = 4.0E0*I_NAI_H3x2y_S_C3_aa-2.0E0*1*I_NAI_Fx2y_S_C3_a-2.0E0*2*I_NAI_Fx2y_S_C3_a;
  abcd[4] = 4.0E0*I_NAI_H3xyz_S_C3_aa-2.0E0*1*I_NAI_Fxyz_S_C3_a-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[5] = 4.0E0*I_NAI_H3x2z_S_C3_aa-2.0E0*1*I_NAI_Fx2z_S_C3_a-2.0E0*2*I_NAI_Fx2z_S_C3_a;
  abcd[6] = 4.0E0*I_NAI_H2x3y_S_C3_aa-2.0E0*1*I_NAI_F3y_S_C3_a;
  abcd[7] = 4.0E0*I_NAI_H2x2yz_S_C3_aa-2.0E0*1*I_NAI_F2yz_S_C3_a;
  abcd[8] = 4.0E0*I_NAI_H2xy2z_S_C3_aa-2.0E0*1*I_NAI_Fy2z_S_C3_a;
  abcd[9] = 4.0E0*I_NAI_H2x3z_S_C3_aa-2.0E0*1*I_NAI_F3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[10] = 4.0E0*I_NAI_H5x_Px_C1003_aa-2.0E0*3*I_NAI_F3x_Px_C1003_a-2.0E0*4*I_NAI_F3x_Px_C1003_a+3*2*I_NAI_Px_Px_C1003;
  abcd[11] = 4.0E0*I_NAI_H4xy_Px_C1003_aa-2.0E0*2*I_NAI_F2xy_Px_C1003_a-2.0E0*3*I_NAI_F2xy_Px_C1003_a+2*1*I_NAI_Py_Px_C1003;
  abcd[12] = 4.0E0*I_NAI_H4xz_Px_C1003_aa-2.0E0*2*I_NAI_F2xz_Px_C1003_a-2.0E0*3*I_NAI_F2xz_Px_C1003_a+2*1*I_NAI_Pz_Px_C1003;
  abcd[13] = 4.0E0*I_NAI_H3x2y_Px_C1003_aa-2.0E0*1*I_NAI_Fx2y_Px_C1003_a-2.0E0*2*I_NAI_Fx2y_Px_C1003_a;
  abcd[14] = 4.0E0*I_NAI_H3xyz_Px_C1003_aa-2.0E0*1*I_NAI_Fxyz_Px_C1003_a-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[15] = 4.0E0*I_NAI_H3x2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2z_Px_C1003_a-2.0E0*2*I_NAI_Fx2z_Px_C1003_a;
  abcd[16] = 4.0E0*I_NAI_H2x3y_Px_C1003_aa-2.0E0*1*I_NAI_F3y_Px_C1003_a;
  abcd[17] = 4.0E0*I_NAI_H2x2yz_Px_C1003_aa-2.0E0*1*I_NAI_F2yz_Px_C1003_a;
  abcd[18] = 4.0E0*I_NAI_H2xy2z_Px_C1003_aa-2.0E0*1*I_NAI_Fy2z_Px_C1003_a;
  abcd[19] = 4.0E0*I_NAI_H2x3z_Px_C1003_aa-2.0E0*1*I_NAI_F3z_Px_C1003_a;
  abcd[20] = 4.0E0*I_NAI_H5x_Py_C1003_aa-2.0E0*3*I_NAI_F3x_Py_C1003_a-2.0E0*4*I_NAI_F3x_Py_C1003_a+3*2*I_NAI_Px_Py_C1003;
  abcd[21] = 4.0E0*I_NAI_H4xy_Py_C1003_aa-2.0E0*2*I_NAI_F2xy_Py_C1003_a-2.0E0*3*I_NAI_F2xy_Py_C1003_a+2*1*I_NAI_Py_Py_C1003;
  abcd[22] = 4.0E0*I_NAI_H4xz_Py_C1003_aa-2.0E0*2*I_NAI_F2xz_Py_C1003_a-2.0E0*3*I_NAI_F2xz_Py_C1003_a+2*1*I_NAI_Pz_Py_C1003;
  abcd[23] = 4.0E0*I_NAI_H3x2y_Py_C1003_aa-2.0E0*1*I_NAI_Fx2y_Py_C1003_a-2.0E0*2*I_NAI_Fx2y_Py_C1003_a;
  abcd[24] = 4.0E0*I_NAI_H3xyz_Py_C1003_aa-2.0E0*1*I_NAI_Fxyz_Py_C1003_a-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[25] = 4.0E0*I_NAI_H3x2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2z_Py_C1003_a-2.0E0*2*I_NAI_Fx2z_Py_C1003_a;
  abcd[26] = 4.0E0*I_NAI_H2x3y_Py_C1003_aa-2.0E0*1*I_NAI_F3y_Py_C1003_a;
  abcd[27] = 4.0E0*I_NAI_H2x2yz_Py_C1003_aa-2.0E0*1*I_NAI_F2yz_Py_C1003_a;
  abcd[28] = 4.0E0*I_NAI_H2xy2z_Py_C1003_aa-2.0E0*1*I_NAI_Fy2z_Py_C1003_a;
  abcd[29] = 4.0E0*I_NAI_H2x3z_Py_C1003_aa-2.0E0*1*I_NAI_F3z_Py_C1003_a;
  abcd[30] = 4.0E0*I_NAI_H5x_Pz_C1003_aa-2.0E0*3*I_NAI_F3x_Pz_C1003_a-2.0E0*4*I_NAI_F3x_Pz_C1003_a+3*2*I_NAI_Px_Pz_C1003;
  abcd[31] = 4.0E0*I_NAI_H4xy_Pz_C1003_aa-2.0E0*2*I_NAI_F2xy_Pz_C1003_a-2.0E0*3*I_NAI_F2xy_Pz_C1003_a+2*1*I_NAI_Py_Pz_C1003;
  abcd[32] = 4.0E0*I_NAI_H4xz_Pz_C1003_aa-2.0E0*2*I_NAI_F2xz_Pz_C1003_a-2.0E0*3*I_NAI_F2xz_Pz_C1003_a+2*1*I_NAI_Pz_Pz_C1003;
  abcd[33] = 4.0E0*I_NAI_H3x2y_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2y_Pz_C1003_a-2.0E0*2*I_NAI_Fx2y_Pz_C1003_a;
  abcd[34] = 4.0E0*I_NAI_H3xyz_Pz_C1003_aa-2.0E0*1*I_NAI_Fxyz_Pz_C1003_a-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[35] = 4.0E0*I_NAI_H3x2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2z_Pz_C1003_a-2.0E0*2*I_NAI_Fx2z_Pz_C1003_a;
  abcd[36] = 4.0E0*I_NAI_H2x3y_Pz_C1003_aa-2.0E0*1*I_NAI_F3y_Pz_C1003_a;
  abcd[37] = 4.0E0*I_NAI_H2x2yz_Pz_C1003_aa-2.0E0*1*I_NAI_F2yz_Pz_C1003_a;
  abcd[38] = 4.0E0*I_NAI_H2xy2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fy2z_Pz_C1003_a;
  abcd[39] = 4.0E0*I_NAI_H2x3z_Pz_C1003_aa-2.0E0*1*I_NAI_F3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[40] = 4.0E0*I_NAI_H4xy_S_C3_aa-2.0E0*3*I_NAI_F2xy_S_C3_a;
  abcd[41] = 4.0E0*I_NAI_H3x2y_S_C3_aa-2.0E0*1*I_NAI_F3x_S_C3_a-2.0E0*2*I_NAI_Fx2y_S_C3_a+2*1*I_NAI_Px_S_C3;
  abcd[42] = 4.0E0*I_NAI_H3xyz_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[43] = 4.0E0*I_NAI_H2x3y_S_C3_aa-2.0E0*2*I_NAI_F2xy_S_C3_a-2.0E0*1*I_NAI_F3y_S_C3_a+2*I_NAI_Py_S_C3;
  abcd[44] = 4.0E0*I_NAI_H2x2yz_S_C3_aa-2.0E0*1*I_NAI_F2xz_S_C3_a-2.0E0*1*I_NAI_F2yz_S_C3_a+1*I_NAI_Pz_S_C3;
  abcd[45] = 4.0E0*I_NAI_H2xy2z_S_C3_aa-2.0E0*1*I_NAI_Fy2z_S_C3_a;
  abcd[46] = 4.0E0*I_NAI_Hx4y_S_C3_aa-2.0E0*3*I_NAI_Fx2y_S_C3_a;
  abcd[47] = 4.0E0*I_NAI_Hx3yz_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[48] = 4.0E0*I_NAI_Hx2y2z_S_C3_aa-2.0E0*1*I_NAI_Fx2z_S_C3_a;
  abcd[49] = 4.0E0*I_NAI_Hxy3z_S_C3_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[50] = 4.0E0*I_NAI_H4xy_Px_C1003_aa-2.0E0*3*I_NAI_F2xy_Px_C1003_a;
  abcd[51] = 4.0E0*I_NAI_H3x2y_Px_C1003_aa-2.0E0*1*I_NAI_F3x_Px_C1003_a-2.0E0*2*I_NAI_Fx2y_Px_C1003_a+2*1*I_NAI_Px_Px_C1003;
  abcd[52] = 4.0E0*I_NAI_H3xyz_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[53] = 4.0E0*I_NAI_H2x3y_Px_C1003_aa-2.0E0*2*I_NAI_F2xy_Px_C1003_a-2.0E0*1*I_NAI_F3y_Px_C1003_a+2*I_NAI_Py_Px_C1003;
  abcd[54] = 4.0E0*I_NAI_H2x2yz_Px_C1003_aa-2.0E0*1*I_NAI_F2xz_Px_C1003_a-2.0E0*1*I_NAI_F2yz_Px_C1003_a+1*I_NAI_Pz_Px_C1003;
  abcd[55] = 4.0E0*I_NAI_H2xy2z_Px_C1003_aa-2.0E0*1*I_NAI_Fy2z_Px_C1003_a;
  abcd[56] = 4.0E0*I_NAI_Hx4y_Px_C1003_aa-2.0E0*3*I_NAI_Fx2y_Px_C1003_a;
  abcd[57] = 4.0E0*I_NAI_Hx3yz_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[58] = 4.0E0*I_NAI_Hx2y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2z_Px_C1003_a;
  abcd[59] = 4.0E0*I_NAI_Hxy3z_Px_C1003_aa;
  abcd[60] = 4.0E0*I_NAI_H4xy_Py_C1003_aa-2.0E0*3*I_NAI_F2xy_Py_C1003_a;
  abcd[61] = 4.0E0*I_NAI_H3x2y_Py_C1003_aa-2.0E0*1*I_NAI_F3x_Py_C1003_a-2.0E0*2*I_NAI_Fx2y_Py_C1003_a+2*1*I_NAI_Px_Py_C1003;
  abcd[62] = 4.0E0*I_NAI_H3xyz_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[63] = 4.0E0*I_NAI_H2x3y_Py_C1003_aa-2.0E0*2*I_NAI_F2xy_Py_C1003_a-2.0E0*1*I_NAI_F3y_Py_C1003_a+2*I_NAI_Py_Py_C1003;
  abcd[64] = 4.0E0*I_NAI_H2x2yz_Py_C1003_aa-2.0E0*1*I_NAI_F2xz_Py_C1003_a-2.0E0*1*I_NAI_F2yz_Py_C1003_a+1*I_NAI_Pz_Py_C1003;
  abcd[65] = 4.0E0*I_NAI_H2xy2z_Py_C1003_aa-2.0E0*1*I_NAI_Fy2z_Py_C1003_a;
  abcd[66] = 4.0E0*I_NAI_Hx4y_Py_C1003_aa-2.0E0*3*I_NAI_Fx2y_Py_C1003_a;
  abcd[67] = 4.0E0*I_NAI_Hx3yz_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[68] = 4.0E0*I_NAI_Hx2y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2z_Py_C1003_a;
  abcd[69] = 4.0E0*I_NAI_Hxy3z_Py_C1003_aa;
  abcd[70] = 4.0E0*I_NAI_H4xy_Pz_C1003_aa-2.0E0*3*I_NAI_F2xy_Pz_C1003_a;
  abcd[71] = 4.0E0*I_NAI_H3x2y_Pz_C1003_aa-2.0E0*1*I_NAI_F3x_Pz_C1003_a-2.0E0*2*I_NAI_Fx2y_Pz_C1003_a+2*1*I_NAI_Px_Pz_C1003;
  abcd[72] = 4.0E0*I_NAI_H3xyz_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[73] = 4.0E0*I_NAI_H2x3y_Pz_C1003_aa-2.0E0*2*I_NAI_F2xy_Pz_C1003_a-2.0E0*1*I_NAI_F3y_Pz_C1003_a+2*I_NAI_Py_Pz_C1003;
  abcd[74] = 4.0E0*I_NAI_H2x2yz_Pz_C1003_aa-2.0E0*1*I_NAI_F2xz_Pz_C1003_a-2.0E0*1*I_NAI_F2yz_Pz_C1003_a+1*I_NAI_Pz_Pz_C1003;
  abcd[75] = 4.0E0*I_NAI_H2xy2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fy2z_Pz_C1003_a;
  abcd[76] = 4.0E0*I_NAI_Hx4y_Pz_C1003_aa-2.0E0*3*I_NAI_Fx2y_Pz_C1003_a;
  abcd[77] = 4.0E0*I_NAI_Hx3yz_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[78] = 4.0E0*I_NAI_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2z_Pz_C1003_a;
  abcd[79] = 4.0E0*I_NAI_Hxy3z_Pz_C1003_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[80] = 4.0E0*I_NAI_H4xz_S_C3_aa-2.0E0*3*I_NAI_F2xz_S_C3_a;
  abcd[81] = 4.0E0*I_NAI_H3xyz_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[82] = 4.0E0*I_NAI_H3x2z_S_C3_aa-2.0E0*1*I_NAI_F3x_S_C3_a-2.0E0*2*I_NAI_Fx2z_S_C3_a+2*1*I_NAI_Px_S_C3;
  abcd[83] = 4.0E0*I_NAI_H2x2yz_S_C3_aa-2.0E0*1*I_NAI_F2yz_S_C3_a;
  abcd[84] = 4.0E0*I_NAI_H2xy2z_S_C3_aa-2.0E0*1*I_NAI_F2xy_S_C3_a-2.0E0*1*I_NAI_Fy2z_S_C3_a+1*I_NAI_Py_S_C3;
  abcd[85] = 4.0E0*I_NAI_H2x3z_S_C3_aa-2.0E0*2*I_NAI_F2xz_S_C3_a-2.0E0*1*I_NAI_F3z_S_C3_a+2*I_NAI_Pz_S_C3;
  abcd[86] = 4.0E0*I_NAI_Hx3yz_S_C3_aa;
  abcd[87] = 4.0E0*I_NAI_Hx2y2z_S_C3_aa-2.0E0*1*I_NAI_Fx2y_S_C3_a;
  abcd[88] = 4.0E0*I_NAI_Hxy3z_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[89] = 4.0E0*I_NAI_Hx4z_S_C3_aa-2.0E0*3*I_NAI_Fx2z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[90] = 4.0E0*I_NAI_H4xz_Px_C1003_aa-2.0E0*3*I_NAI_F2xz_Px_C1003_a;
  abcd[91] = 4.0E0*I_NAI_H3xyz_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[92] = 4.0E0*I_NAI_H3x2z_Px_C1003_aa-2.0E0*1*I_NAI_F3x_Px_C1003_a-2.0E0*2*I_NAI_Fx2z_Px_C1003_a+2*1*I_NAI_Px_Px_C1003;
  abcd[93] = 4.0E0*I_NAI_H2x2yz_Px_C1003_aa-2.0E0*1*I_NAI_F2yz_Px_C1003_a;
  abcd[94] = 4.0E0*I_NAI_H2xy2z_Px_C1003_aa-2.0E0*1*I_NAI_F2xy_Px_C1003_a-2.0E0*1*I_NAI_Fy2z_Px_C1003_a+1*I_NAI_Py_Px_C1003;
  abcd[95] = 4.0E0*I_NAI_H2x3z_Px_C1003_aa-2.0E0*2*I_NAI_F2xz_Px_C1003_a-2.0E0*1*I_NAI_F3z_Px_C1003_a+2*I_NAI_Pz_Px_C1003;
  abcd[96] = 4.0E0*I_NAI_Hx3yz_Px_C1003_aa;
  abcd[97] = 4.0E0*I_NAI_Hx2y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2y_Px_C1003_a;
  abcd[98] = 4.0E0*I_NAI_Hxy3z_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[99] = 4.0E0*I_NAI_Hx4z_Px_C1003_aa-2.0E0*3*I_NAI_Fx2z_Px_C1003_a;
  abcd[100] = 4.0E0*I_NAI_H4xz_Py_C1003_aa-2.0E0*3*I_NAI_F2xz_Py_C1003_a;
  abcd[101] = 4.0E0*I_NAI_H3xyz_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[102] = 4.0E0*I_NAI_H3x2z_Py_C1003_aa-2.0E0*1*I_NAI_F3x_Py_C1003_a-2.0E0*2*I_NAI_Fx2z_Py_C1003_a+2*1*I_NAI_Px_Py_C1003;
  abcd[103] = 4.0E0*I_NAI_H2x2yz_Py_C1003_aa-2.0E0*1*I_NAI_F2yz_Py_C1003_a;
  abcd[104] = 4.0E0*I_NAI_H2xy2z_Py_C1003_aa-2.0E0*1*I_NAI_F2xy_Py_C1003_a-2.0E0*1*I_NAI_Fy2z_Py_C1003_a+1*I_NAI_Py_Py_C1003;
  abcd[105] = 4.0E0*I_NAI_H2x3z_Py_C1003_aa-2.0E0*2*I_NAI_F2xz_Py_C1003_a-2.0E0*1*I_NAI_F3z_Py_C1003_a+2*I_NAI_Pz_Py_C1003;
  abcd[106] = 4.0E0*I_NAI_Hx3yz_Py_C1003_aa;
  abcd[107] = 4.0E0*I_NAI_Hx2y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2y_Py_C1003_a;
  abcd[108] = 4.0E0*I_NAI_Hxy3z_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[109] = 4.0E0*I_NAI_Hx4z_Py_C1003_aa-2.0E0*3*I_NAI_Fx2z_Py_C1003_a;
  abcd[110] = 4.0E0*I_NAI_H4xz_Pz_C1003_aa-2.0E0*3*I_NAI_F2xz_Pz_C1003_a;
  abcd[111] = 4.0E0*I_NAI_H3xyz_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[112] = 4.0E0*I_NAI_H3x2z_Pz_C1003_aa-2.0E0*1*I_NAI_F3x_Pz_C1003_a-2.0E0*2*I_NAI_Fx2z_Pz_C1003_a+2*1*I_NAI_Px_Pz_C1003;
  abcd[113] = 4.0E0*I_NAI_H2x2yz_Pz_C1003_aa-2.0E0*1*I_NAI_F2yz_Pz_C1003_a;
  abcd[114] = 4.0E0*I_NAI_H2xy2z_Pz_C1003_aa-2.0E0*1*I_NAI_F2xy_Pz_C1003_a-2.0E0*1*I_NAI_Fy2z_Pz_C1003_a+1*I_NAI_Py_Pz_C1003;
  abcd[115] = 4.0E0*I_NAI_H2x3z_Pz_C1003_aa-2.0E0*2*I_NAI_F2xz_Pz_C1003_a-2.0E0*1*I_NAI_F3z_Pz_C1003_a+2*I_NAI_Pz_Pz_C1003;
  abcd[116] = 4.0E0*I_NAI_Hx3yz_Pz_C1003_aa;
  abcd[117] = 4.0E0*I_NAI_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2y_Pz_C1003_a;
  abcd[118] = 4.0E0*I_NAI_Hxy3z_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[119] = 4.0E0*I_NAI_Hx4z_Pz_C1003_aa-2.0E0*3*I_NAI_Fx2z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[120] = 4.0E0*I_NAI_H3x2y_S_C3_aa-2.0E0*1*I_NAI_F3x_S_C3_a;
  abcd[121] = 4.0E0*I_NAI_H2x3y_S_C3_aa-2.0E0*1*I_NAI_F2xy_S_C3_a-2.0E0*2*I_NAI_F2xy_S_C3_a;
  abcd[122] = 4.0E0*I_NAI_H2x2yz_S_C3_aa-2.0E0*1*I_NAI_F2xz_S_C3_a;
  abcd[123] = 4.0E0*I_NAI_Hx4y_S_C3_aa-2.0E0*2*I_NAI_Fx2y_S_C3_a-2.0E0*3*I_NAI_Fx2y_S_C3_a+2*1*I_NAI_Px_S_C3;
  abcd[124] = 4.0E0*I_NAI_Hx3yz_S_C3_aa-2.0E0*1*I_NAI_Fxyz_S_C3_a-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[125] = 4.0E0*I_NAI_Hx2y2z_S_C3_aa-2.0E0*1*I_NAI_Fx2z_S_C3_a;
  abcd[126] = 4.0E0*I_NAI_H5y_S_C3_aa-2.0E0*3*I_NAI_F3y_S_C3_a-2.0E0*4*I_NAI_F3y_S_C3_a+3*2*I_NAI_Py_S_C3;
  abcd[127] = 4.0E0*I_NAI_H4yz_S_C3_aa-2.0E0*2*I_NAI_F2yz_S_C3_a-2.0E0*3*I_NAI_F2yz_S_C3_a+2*1*I_NAI_Pz_S_C3;
  abcd[128] = 4.0E0*I_NAI_H3y2z_S_C3_aa-2.0E0*1*I_NAI_Fy2z_S_C3_a-2.0E0*2*I_NAI_Fy2z_S_C3_a;
  abcd[129] = 4.0E0*I_NAI_H2y3z_S_C3_aa-2.0E0*1*I_NAI_F3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[130] = 4.0E0*I_NAI_H3x2y_Px_C1003_aa-2.0E0*1*I_NAI_F3x_Px_C1003_a;
  abcd[131] = 4.0E0*I_NAI_H2x3y_Px_C1003_aa-2.0E0*1*I_NAI_F2xy_Px_C1003_a-2.0E0*2*I_NAI_F2xy_Px_C1003_a;
  abcd[132] = 4.0E0*I_NAI_H2x2yz_Px_C1003_aa-2.0E0*1*I_NAI_F2xz_Px_C1003_a;
  abcd[133] = 4.0E0*I_NAI_Hx4y_Px_C1003_aa-2.0E0*2*I_NAI_Fx2y_Px_C1003_a-2.0E0*3*I_NAI_Fx2y_Px_C1003_a+2*1*I_NAI_Px_Px_C1003;
  abcd[134] = 4.0E0*I_NAI_Hx3yz_Px_C1003_aa-2.0E0*1*I_NAI_Fxyz_Px_C1003_a-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[135] = 4.0E0*I_NAI_Hx2y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2z_Px_C1003_a;
  abcd[136] = 4.0E0*I_NAI_H5y_Px_C1003_aa-2.0E0*3*I_NAI_F3y_Px_C1003_a-2.0E0*4*I_NAI_F3y_Px_C1003_a+3*2*I_NAI_Py_Px_C1003;
  abcd[137] = 4.0E0*I_NAI_H4yz_Px_C1003_aa-2.0E0*2*I_NAI_F2yz_Px_C1003_a-2.0E0*3*I_NAI_F2yz_Px_C1003_a+2*1*I_NAI_Pz_Px_C1003;
  abcd[138] = 4.0E0*I_NAI_H3y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fy2z_Px_C1003_a-2.0E0*2*I_NAI_Fy2z_Px_C1003_a;
  abcd[139] = 4.0E0*I_NAI_H2y3z_Px_C1003_aa-2.0E0*1*I_NAI_F3z_Px_C1003_a;
  abcd[140] = 4.0E0*I_NAI_H3x2y_Py_C1003_aa-2.0E0*1*I_NAI_F3x_Py_C1003_a;
  abcd[141] = 4.0E0*I_NAI_H2x3y_Py_C1003_aa-2.0E0*1*I_NAI_F2xy_Py_C1003_a-2.0E0*2*I_NAI_F2xy_Py_C1003_a;
  abcd[142] = 4.0E0*I_NAI_H2x2yz_Py_C1003_aa-2.0E0*1*I_NAI_F2xz_Py_C1003_a;
  abcd[143] = 4.0E0*I_NAI_Hx4y_Py_C1003_aa-2.0E0*2*I_NAI_Fx2y_Py_C1003_a-2.0E0*3*I_NAI_Fx2y_Py_C1003_a+2*1*I_NAI_Px_Py_C1003;
  abcd[144] = 4.0E0*I_NAI_Hx3yz_Py_C1003_aa-2.0E0*1*I_NAI_Fxyz_Py_C1003_a-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[145] = 4.0E0*I_NAI_Hx2y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2z_Py_C1003_a;
  abcd[146] = 4.0E0*I_NAI_H5y_Py_C1003_aa-2.0E0*3*I_NAI_F3y_Py_C1003_a-2.0E0*4*I_NAI_F3y_Py_C1003_a+3*2*I_NAI_Py_Py_C1003;
  abcd[147] = 4.0E0*I_NAI_H4yz_Py_C1003_aa-2.0E0*2*I_NAI_F2yz_Py_C1003_a-2.0E0*3*I_NAI_F2yz_Py_C1003_a+2*1*I_NAI_Pz_Py_C1003;
  abcd[148] = 4.0E0*I_NAI_H3y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fy2z_Py_C1003_a-2.0E0*2*I_NAI_Fy2z_Py_C1003_a;
  abcd[149] = 4.0E0*I_NAI_H2y3z_Py_C1003_aa-2.0E0*1*I_NAI_F3z_Py_C1003_a;
  abcd[150] = 4.0E0*I_NAI_H3x2y_Pz_C1003_aa-2.0E0*1*I_NAI_F3x_Pz_C1003_a;
  abcd[151] = 4.0E0*I_NAI_H2x3y_Pz_C1003_aa-2.0E0*1*I_NAI_F2xy_Pz_C1003_a-2.0E0*2*I_NAI_F2xy_Pz_C1003_a;
  abcd[152] = 4.0E0*I_NAI_H2x2yz_Pz_C1003_aa-2.0E0*1*I_NAI_F2xz_Pz_C1003_a;
  abcd[153] = 4.0E0*I_NAI_Hx4y_Pz_C1003_aa-2.0E0*2*I_NAI_Fx2y_Pz_C1003_a-2.0E0*3*I_NAI_Fx2y_Pz_C1003_a+2*1*I_NAI_Px_Pz_C1003;
  abcd[154] = 4.0E0*I_NAI_Hx3yz_Pz_C1003_aa-2.0E0*1*I_NAI_Fxyz_Pz_C1003_a-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[155] = 4.0E0*I_NAI_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2z_Pz_C1003_a;
  abcd[156] = 4.0E0*I_NAI_H5y_Pz_C1003_aa-2.0E0*3*I_NAI_F3y_Pz_C1003_a-2.0E0*4*I_NAI_F3y_Pz_C1003_a+3*2*I_NAI_Py_Pz_C1003;
  abcd[157] = 4.0E0*I_NAI_H4yz_Pz_C1003_aa-2.0E0*2*I_NAI_F2yz_Pz_C1003_a-2.0E0*3*I_NAI_F2yz_Pz_C1003_a+2*1*I_NAI_Pz_Pz_C1003;
  abcd[158] = 4.0E0*I_NAI_H3y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fy2z_Pz_C1003_a-2.0E0*2*I_NAI_Fy2z_Pz_C1003_a;
  abcd[159] = 4.0E0*I_NAI_H2y3z_Pz_C1003_aa-2.0E0*1*I_NAI_F3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[160] = 4.0E0*I_NAI_H3xyz_S_C3_aa;
  abcd[161] = 4.0E0*I_NAI_H2x2yz_S_C3_aa-2.0E0*1*I_NAI_F2xz_S_C3_a;
  abcd[162] = 4.0E0*I_NAI_H2xy2z_S_C3_aa-2.0E0*1*I_NAI_F2xy_S_C3_a;
  abcd[163] = 4.0E0*I_NAI_Hx3yz_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[164] = 4.0E0*I_NAI_Hx2y2z_S_C3_aa-2.0E0*1*I_NAI_Fx2y_S_C3_a-2.0E0*1*I_NAI_Fx2z_S_C3_a+1*I_NAI_Px_S_C3;
  abcd[165] = 4.0E0*I_NAI_Hxy3z_S_C3_aa-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[166] = 4.0E0*I_NAI_H4yz_S_C3_aa-2.0E0*3*I_NAI_F2yz_S_C3_a;
  abcd[167] = 4.0E0*I_NAI_H3y2z_S_C3_aa-2.0E0*1*I_NAI_F3y_S_C3_a-2.0E0*2*I_NAI_Fy2z_S_C3_a+2*1*I_NAI_Py_S_C3;
  abcd[168] = 4.0E0*I_NAI_H2y3z_S_C3_aa-2.0E0*2*I_NAI_F2yz_S_C3_a-2.0E0*1*I_NAI_F3z_S_C3_a+2*I_NAI_Pz_S_C3;
  abcd[169] = 4.0E0*I_NAI_Hy4z_S_C3_aa-2.0E0*3*I_NAI_Fy2z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[170] = 4.0E0*I_NAI_H3xyz_Px_C1003_aa;
  abcd[171] = 4.0E0*I_NAI_H2x2yz_Px_C1003_aa-2.0E0*1*I_NAI_F2xz_Px_C1003_a;
  abcd[172] = 4.0E0*I_NAI_H2xy2z_Px_C1003_aa-2.0E0*1*I_NAI_F2xy_Px_C1003_a;
  abcd[173] = 4.0E0*I_NAI_Hx3yz_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[174] = 4.0E0*I_NAI_Hx2y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2y_Px_C1003_a-2.0E0*1*I_NAI_Fx2z_Px_C1003_a+1*I_NAI_Px_Px_C1003;
  abcd[175] = 4.0E0*I_NAI_Hxy3z_Px_C1003_aa-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[176] = 4.0E0*I_NAI_H4yz_Px_C1003_aa-2.0E0*3*I_NAI_F2yz_Px_C1003_a;
  abcd[177] = 4.0E0*I_NAI_H3y2z_Px_C1003_aa-2.0E0*1*I_NAI_F3y_Px_C1003_a-2.0E0*2*I_NAI_Fy2z_Px_C1003_a+2*1*I_NAI_Py_Px_C1003;
  abcd[178] = 4.0E0*I_NAI_H2y3z_Px_C1003_aa-2.0E0*2*I_NAI_F2yz_Px_C1003_a-2.0E0*1*I_NAI_F3z_Px_C1003_a+2*I_NAI_Pz_Px_C1003;
  abcd[179] = 4.0E0*I_NAI_Hy4z_Px_C1003_aa-2.0E0*3*I_NAI_Fy2z_Px_C1003_a;
  abcd[180] = 4.0E0*I_NAI_H3xyz_Py_C1003_aa;
  abcd[181] = 4.0E0*I_NAI_H2x2yz_Py_C1003_aa-2.0E0*1*I_NAI_F2xz_Py_C1003_a;
  abcd[182] = 4.0E0*I_NAI_H2xy2z_Py_C1003_aa-2.0E0*1*I_NAI_F2xy_Py_C1003_a;
  abcd[183] = 4.0E0*I_NAI_Hx3yz_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[184] = 4.0E0*I_NAI_Hx2y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2y_Py_C1003_a-2.0E0*1*I_NAI_Fx2z_Py_C1003_a+1*I_NAI_Px_Py_C1003;
  abcd[185] = 4.0E0*I_NAI_Hxy3z_Py_C1003_aa-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[186] = 4.0E0*I_NAI_H4yz_Py_C1003_aa-2.0E0*3*I_NAI_F2yz_Py_C1003_a;
  abcd[187] = 4.0E0*I_NAI_H3y2z_Py_C1003_aa-2.0E0*1*I_NAI_F3y_Py_C1003_a-2.0E0*2*I_NAI_Fy2z_Py_C1003_a+2*1*I_NAI_Py_Py_C1003;
  abcd[188] = 4.0E0*I_NAI_H2y3z_Py_C1003_aa-2.0E0*2*I_NAI_F2yz_Py_C1003_a-2.0E0*1*I_NAI_F3z_Py_C1003_a+2*I_NAI_Pz_Py_C1003;
  abcd[189] = 4.0E0*I_NAI_Hy4z_Py_C1003_aa-2.0E0*3*I_NAI_Fy2z_Py_C1003_a;
  abcd[190] = 4.0E0*I_NAI_H3xyz_Pz_C1003_aa;
  abcd[191] = 4.0E0*I_NAI_H2x2yz_Pz_C1003_aa-2.0E0*1*I_NAI_F2xz_Pz_C1003_a;
  abcd[192] = 4.0E0*I_NAI_H2xy2z_Pz_C1003_aa-2.0E0*1*I_NAI_F2xy_Pz_C1003_a;
  abcd[193] = 4.0E0*I_NAI_Hx3yz_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[194] = 4.0E0*I_NAI_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2y_Pz_C1003_a-2.0E0*1*I_NAI_Fx2z_Pz_C1003_a+1*I_NAI_Px_Pz_C1003;
  abcd[195] = 4.0E0*I_NAI_Hxy3z_Pz_C1003_aa-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[196] = 4.0E0*I_NAI_H4yz_Pz_C1003_aa-2.0E0*3*I_NAI_F2yz_Pz_C1003_a;
  abcd[197] = 4.0E0*I_NAI_H3y2z_Pz_C1003_aa-2.0E0*1*I_NAI_F3y_Pz_C1003_a-2.0E0*2*I_NAI_Fy2z_Pz_C1003_a+2*1*I_NAI_Py_Pz_C1003;
  abcd[198] = 4.0E0*I_NAI_H2y3z_Pz_C1003_aa-2.0E0*2*I_NAI_F2yz_Pz_C1003_a-2.0E0*1*I_NAI_F3z_Pz_C1003_a+2*I_NAI_Pz_Pz_C1003;
  abcd[199] = 4.0E0*I_NAI_Hy4z_Pz_C1003_aa-2.0E0*3*I_NAI_Fy2z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C3_aa
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_F_S_C3_a
   * RHS shell quartet name: SQ_NAI_P_S_C3
   ************************************************************/
  abcd[200] = 4.0E0*I_NAI_H3x2z_S_C3_aa-2.0E0*1*I_NAI_F3x_S_C3_a;
  abcd[201] = 4.0E0*I_NAI_H2xy2z_S_C3_aa-2.0E0*1*I_NAI_F2xy_S_C3_a;
  abcd[202] = 4.0E0*I_NAI_H2x3z_S_C3_aa-2.0E0*1*I_NAI_F2xz_S_C3_a-2.0E0*2*I_NAI_F2xz_S_C3_a;
  abcd[203] = 4.0E0*I_NAI_Hx2y2z_S_C3_aa-2.0E0*1*I_NAI_Fx2y_S_C3_a;
  abcd[204] = 4.0E0*I_NAI_Hxy3z_S_C3_aa-2.0E0*1*I_NAI_Fxyz_S_C3_a-2.0E0*2*I_NAI_Fxyz_S_C3_a;
  abcd[205] = 4.0E0*I_NAI_Hx4z_S_C3_aa-2.0E0*2*I_NAI_Fx2z_S_C3_a-2.0E0*3*I_NAI_Fx2z_S_C3_a+2*1*I_NAI_Px_S_C3;
  abcd[206] = 4.0E0*I_NAI_H3y2z_S_C3_aa-2.0E0*1*I_NAI_F3y_S_C3_a;
  abcd[207] = 4.0E0*I_NAI_H2y3z_S_C3_aa-2.0E0*1*I_NAI_F2yz_S_C3_a-2.0E0*2*I_NAI_F2yz_S_C3_a;
  abcd[208] = 4.0E0*I_NAI_Hy4z_S_C3_aa-2.0E0*2*I_NAI_Fy2z_S_C3_a-2.0E0*3*I_NAI_Fy2z_S_C3_a+2*1*I_NAI_Py_S_C3;
  abcd[209] = 4.0E0*I_NAI_H5z_S_C3_aa-2.0E0*3*I_NAI_F3z_S_C3_a-2.0E0*4*I_NAI_F3z_S_C3_a+3*2*I_NAI_Pz_S_C3;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_C1003_aa
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_F_P_C1003_a
   * RHS shell quartet name: SQ_NAI_P_P_C1003
   ************************************************************/
  abcd[210] = 4.0E0*I_NAI_H3x2z_Px_C1003_aa-2.0E0*1*I_NAI_F3x_Px_C1003_a;
  abcd[211] = 4.0E0*I_NAI_H2xy2z_Px_C1003_aa-2.0E0*1*I_NAI_F2xy_Px_C1003_a;
  abcd[212] = 4.0E0*I_NAI_H2x3z_Px_C1003_aa-2.0E0*1*I_NAI_F2xz_Px_C1003_a-2.0E0*2*I_NAI_F2xz_Px_C1003_a;
  abcd[213] = 4.0E0*I_NAI_Hx2y2z_Px_C1003_aa-2.0E0*1*I_NAI_Fx2y_Px_C1003_a;
  abcd[214] = 4.0E0*I_NAI_Hxy3z_Px_C1003_aa-2.0E0*1*I_NAI_Fxyz_Px_C1003_a-2.0E0*2*I_NAI_Fxyz_Px_C1003_a;
  abcd[215] = 4.0E0*I_NAI_Hx4z_Px_C1003_aa-2.0E0*2*I_NAI_Fx2z_Px_C1003_a-2.0E0*3*I_NAI_Fx2z_Px_C1003_a+2*1*I_NAI_Px_Px_C1003;
  abcd[216] = 4.0E0*I_NAI_H3y2z_Px_C1003_aa-2.0E0*1*I_NAI_F3y_Px_C1003_a;
  abcd[217] = 4.0E0*I_NAI_H2y3z_Px_C1003_aa-2.0E0*1*I_NAI_F2yz_Px_C1003_a-2.0E0*2*I_NAI_F2yz_Px_C1003_a;
  abcd[218] = 4.0E0*I_NAI_Hy4z_Px_C1003_aa-2.0E0*2*I_NAI_Fy2z_Px_C1003_a-2.0E0*3*I_NAI_Fy2z_Px_C1003_a+2*1*I_NAI_Py_Px_C1003;
  abcd[219] = 4.0E0*I_NAI_H5z_Px_C1003_aa-2.0E0*3*I_NAI_F3z_Px_C1003_a-2.0E0*4*I_NAI_F3z_Px_C1003_a+3*2*I_NAI_Pz_Px_C1003;
  abcd[220] = 4.0E0*I_NAI_H3x2z_Py_C1003_aa-2.0E0*1*I_NAI_F3x_Py_C1003_a;
  abcd[221] = 4.0E0*I_NAI_H2xy2z_Py_C1003_aa-2.0E0*1*I_NAI_F2xy_Py_C1003_a;
  abcd[222] = 4.0E0*I_NAI_H2x3z_Py_C1003_aa-2.0E0*1*I_NAI_F2xz_Py_C1003_a-2.0E0*2*I_NAI_F2xz_Py_C1003_a;
  abcd[223] = 4.0E0*I_NAI_Hx2y2z_Py_C1003_aa-2.0E0*1*I_NAI_Fx2y_Py_C1003_a;
  abcd[224] = 4.0E0*I_NAI_Hxy3z_Py_C1003_aa-2.0E0*1*I_NAI_Fxyz_Py_C1003_a-2.0E0*2*I_NAI_Fxyz_Py_C1003_a;
  abcd[225] = 4.0E0*I_NAI_Hx4z_Py_C1003_aa-2.0E0*2*I_NAI_Fx2z_Py_C1003_a-2.0E0*3*I_NAI_Fx2z_Py_C1003_a+2*1*I_NAI_Px_Py_C1003;
  abcd[226] = 4.0E0*I_NAI_H3y2z_Py_C1003_aa-2.0E0*1*I_NAI_F3y_Py_C1003_a;
  abcd[227] = 4.0E0*I_NAI_H2y3z_Py_C1003_aa-2.0E0*1*I_NAI_F2yz_Py_C1003_a-2.0E0*2*I_NAI_F2yz_Py_C1003_a;
  abcd[228] = 4.0E0*I_NAI_Hy4z_Py_C1003_aa-2.0E0*2*I_NAI_Fy2z_Py_C1003_a-2.0E0*3*I_NAI_Fy2z_Py_C1003_a+2*1*I_NAI_Py_Py_C1003;
  abcd[229] = 4.0E0*I_NAI_H5z_Py_C1003_aa-2.0E0*3*I_NAI_F3z_Py_C1003_a-2.0E0*4*I_NAI_F3z_Py_C1003_a+3*2*I_NAI_Pz_Py_C1003;
  abcd[230] = 4.0E0*I_NAI_H3x2z_Pz_C1003_aa-2.0E0*1*I_NAI_F3x_Pz_C1003_a;
  abcd[231] = 4.0E0*I_NAI_H2xy2z_Pz_C1003_aa-2.0E0*1*I_NAI_F2xy_Pz_C1003_a;
  abcd[232] = 4.0E0*I_NAI_H2x3z_Pz_C1003_aa-2.0E0*1*I_NAI_F2xz_Pz_C1003_a-2.0E0*2*I_NAI_F2xz_Pz_C1003_a;
  abcd[233] = 4.0E0*I_NAI_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_NAI_Fx2y_Pz_C1003_a;
  abcd[234] = 4.0E0*I_NAI_Hxy3z_Pz_C1003_aa-2.0E0*1*I_NAI_Fxyz_Pz_C1003_a-2.0E0*2*I_NAI_Fxyz_Pz_C1003_a;
  abcd[235] = 4.0E0*I_NAI_Hx4z_Pz_C1003_aa-2.0E0*2*I_NAI_Fx2z_Pz_C1003_a-2.0E0*3*I_NAI_Fx2z_Pz_C1003_a+2*1*I_NAI_Px_Pz_C1003;
  abcd[236] = 4.0E0*I_NAI_H3y2z_Pz_C1003_aa-2.0E0*1*I_NAI_F3y_Pz_C1003_a;
  abcd[237] = 4.0E0*I_NAI_H2y3z_Pz_C1003_aa-2.0E0*1*I_NAI_F2yz_Pz_C1003_a-2.0E0*2*I_NAI_F2yz_Pz_C1003_a;
  abcd[238] = 4.0E0*I_NAI_Hy4z_Pz_C1003_aa-2.0E0*2*I_NAI_Fy2z_Pz_C1003_a-2.0E0*3*I_NAI_Fy2z_Pz_C1003_a+2*1*I_NAI_Py_Pz_C1003;
  abcd[239] = 4.0E0*I_NAI_H5z_Pz_C1003_aa-2.0E0*3*I_NAI_F3z_Pz_C1003_a-2.0E0*4*I_NAI_F3z_Pz_C1003_a+3*2*I_NAI_Pz_Pz_C1003;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[240] = 4.0E0*I_NAI_G4x_Px_C3_ab-2.0E0*3*I_NAI_D2x_Px_C3_b;
  abcd[241] = 4.0E0*I_NAI_G3xy_Px_C3_ab-2.0E0*2*I_NAI_Dxy_Px_C3_b;
  abcd[242] = 4.0E0*I_NAI_G3xz_Px_C3_ab-2.0E0*2*I_NAI_Dxz_Px_C3_b;
  abcd[243] = 4.0E0*I_NAI_G2x2y_Px_C3_ab-2.0E0*1*I_NAI_D2y_Px_C3_b;
  abcd[244] = 4.0E0*I_NAI_G2xyz_Px_C3_ab-2.0E0*1*I_NAI_Dyz_Px_C3_b;
  abcd[245] = 4.0E0*I_NAI_G2x2z_Px_C3_ab-2.0E0*1*I_NAI_D2z_Px_C3_b;
  abcd[246] = 4.0E0*I_NAI_Gx3y_Px_C3_ab;
  abcd[247] = 4.0E0*I_NAI_Gx2yz_Px_C3_ab;
  abcd[248] = 4.0E0*I_NAI_Gxy2z_Px_C3_ab;
  abcd[249] = 4.0E0*I_NAI_Gx3z_Px_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[250] = 4.0E0*I_NAI_G4x_D2x_C1003_ab-2.0E0*1*I_NAI_G4x_S_C1003_a-2.0E0*3*I_NAI_D2x_D2x_C1003_b+3*1*I_NAI_D2x_S_C1003;
  abcd[251] = 4.0E0*I_NAI_G3xy_D2x_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2x_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[252] = 4.0E0*I_NAI_G3xz_D2x_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2x_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[253] = 4.0E0*I_NAI_G2x2y_D2x_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2y_D2x_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[254] = 4.0E0*I_NAI_G2xyz_D2x_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a-2.0E0*1*I_NAI_Dyz_D2x_C1003_b+1*I_NAI_Dyz_S_C1003;
  abcd[255] = 4.0E0*I_NAI_G2x2z_D2x_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2x_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[256] = 4.0E0*I_NAI_Gx3y_D2x_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a;
  abcd[257] = 4.0E0*I_NAI_Gx2yz_D2x_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[258] = 4.0E0*I_NAI_Gxy2z_D2x_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[259] = 4.0E0*I_NAI_Gx3z_D2x_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a;
  abcd[260] = 4.0E0*I_NAI_G4x_Dxy_C1003_ab-2.0E0*3*I_NAI_D2x_Dxy_C1003_b;
  abcd[261] = 4.0E0*I_NAI_G3xy_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxy_C1003_b;
  abcd[262] = 4.0E0*I_NAI_G3xz_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxy_C1003_b;
  abcd[263] = 4.0E0*I_NAI_G2x2y_Dxy_C1003_ab-2.0E0*1*I_NAI_D2y_Dxy_C1003_b;
  abcd[264] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab-2.0E0*1*I_NAI_Dyz_Dxy_C1003_b;
  abcd[265] = 4.0E0*I_NAI_G2x2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2z_Dxy_C1003_b;
  abcd[266] = 4.0E0*I_NAI_Gx3y_Dxy_C1003_ab;
  abcd[267] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab;
  abcd[268] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab;
  abcd[269] = 4.0E0*I_NAI_Gx3z_Dxy_C1003_ab;
  abcd[270] = 4.0E0*I_NAI_G4x_Dxz_C1003_ab-2.0E0*3*I_NAI_D2x_Dxz_C1003_b;
  abcd[271] = 4.0E0*I_NAI_G3xy_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxz_C1003_b;
  abcd[272] = 4.0E0*I_NAI_G3xz_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxz_C1003_b;
  abcd[273] = 4.0E0*I_NAI_G2x2y_Dxz_C1003_ab-2.0E0*1*I_NAI_D2y_Dxz_C1003_b;
  abcd[274] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab-2.0E0*1*I_NAI_Dyz_Dxz_C1003_b;
  abcd[275] = 4.0E0*I_NAI_G2x2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2z_Dxz_C1003_b;
  abcd[276] = 4.0E0*I_NAI_Gx3y_Dxz_C1003_ab;
  abcd[277] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab;
  abcd[278] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab;
  abcd[279] = 4.0E0*I_NAI_Gx3z_Dxz_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[280] = 4.0E0*I_NAI_G4x_Py_C3_ab-2.0E0*3*I_NAI_D2x_Py_C3_b;
  abcd[281] = 4.0E0*I_NAI_G3xy_Py_C3_ab-2.0E0*2*I_NAI_Dxy_Py_C3_b;
  abcd[282] = 4.0E0*I_NAI_G3xz_Py_C3_ab-2.0E0*2*I_NAI_Dxz_Py_C3_b;
  abcd[283] = 4.0E0*I_NAI_G2x2y_Py_C3_ab-2.0E0*1*I_NAI_D2y_Py_C3_b;
  abcd[284] = 4.0E0*I_NAI_G2xyz_Py_C3_ab-2.0E0*1*I_NAI_Dyz_Py_C3_b;
  abcd[285] = 4.0E0*I_NAI_G2x2z_Py_C3_ab-2.0E0*1*I_NAI_D2z_Py_C3_b;
  abcd[286] = 4.0E0*I_NAI_Gx3y_Py_C3_ab;
  abcd[287] = 4.0E0*I_NAI_Gx2yz_Py_C3_ab;
  abcd[288] = 4.0E0*I_NAI_Gxy2z_Py_C3_ab;
  abcd[289] = 4.0E0*I_NAI_Gx3z_Py_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[290] = 4.0E0*I_NAI_G4x_Dxy_C1003_ab-2.0E0*3*I_NAI_D2x_Dxy_C1003_b;
  abcd[291] = 4.0E0*I_NAI_G3xy_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxy_C1003_b;
  abcd[292] = 4.0E0*I_NAI_G3xz_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxy_C1003_b;
  abcd[293] = 4.0E0*I_NAI_G2x2y_Dxy_C1003_ab-2.0E0*1*I_NAI_D2y_Dxy_C1003_b;
  abcd[294] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab-2.0E0*1*I_NAI_Dyz_Dxy_C1003_b;
  abcd[295] = 4.0E0*I_NAI_G2x2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2z_Dxy_C1003_b;
  abcd[296] = 4.0E0*I_NAI_Gx3y_Dxy_C1003_ab;
  abcd[297] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab;
  abcd[298] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab;
  abcd[299] = 4.0E0*I_NAI_Gx3z_Dxy_C1003_ab;
  abcd[300] = 4.0E0*I_NAI_G4x_D2y_C1003_ab-2.0E0*1*I_NAI_G4x_S_C1003_a-2.0E0*3*I_NAI_D2x_D2y_C1003_b+3*1*I_NAI_D2x_S_C1003;
  abcd[301] = 4.0E0*I_NAI_G3xy_D2y_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2y_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[302] = 4.0E0*I_NAI_G3xz_D2y_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2y_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[303] = 4.0E0*I_NAI_G2x2y_D2y_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2y_D2y_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[304] = 4.0E0*I_NAI_G2xyz_D2y_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a-2.0E0*1*I_NAI_Dyz_D2y_C1003_b+1*I_NAI_Dyz_S_C1003;
  abcd[305] = 4.0E0*I_NAI_G2x2z_D2y_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2y_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[306] = 4.0E0*I_NAI_Gx3y_D2y_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a;
  abcd[307] = 4.0E0*I_NAI_Gx2yz_D2y_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[308] = 4.0E0*I_NAI_Gxy2z_D2y_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[309] = 4.0E0*I_NAI_Gx3z_D2y_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a;
  abcd[310] = 4.0E0*I_NAI_G4x_Dyz_C1003_ab-2.0E0*3*I_NAI_D2x_Dyz_C1003_b;
  abcd[311] = 4.0E0*I_NAI_G3xy_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dyz_C1003_b;
  abcd[312] = 4.0E0*I_NAI_G3xz_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dyz_C1003_b;
  abcd[313] = 4.0E0*I_NAI_G2x2y_Dyz_C1003_ab-2.0E0*1*I_NAI_D2y_Dyz_C1003_b;
  abcd[314] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab-2.0E0*1*I_NAI_Dyz_Dyz_C1003_b;
  abcd[315] = 4.0E0*I_NAI_G2x2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2z_Dyz_C1003_b;
  abcd[316] = 4.0E0*I_NAI_Gx3y_Dyz_C1003_ab;
  abcd[317] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab;
  abcd[318] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab;
  abcd[319] = 4.0E0*I_NAI_Gx3z_Dyz_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[320] = 4.0E0*I_NAI_G4x_Pz_C3_ab-2.0E0*3*I_NAI_D2x_Pz_C3_b;
  abcd[321] = 4.0E0*I_NAI_G3xy_Pz_C3_ab-2.0E0*2*I_NAI_Dxy_Pz_C3_b;
  abcd[322] = 4.0E0*I_NAI_G3xz_Pz_C3_ab-2.0E0*2*I_NAI_Dxz_Pz_C3_b;
  abcd[323] = 4.0E0*I_NAI_G2x2y_Pz_C3_ab-2.0E0*1*I_NAI_D2y_Pz_C3_b;
  abcd[324] = 4.0E0*I_NAI_G2xyz_Pz_C3_ab-2.0E0*1*I_NAI_Dyz_Pz_C3_b;
  abcd[325] = 4.0E0*I_NAI_G2x2z_Pz_C3_ab-2.0E0*1*I_NAI_D2z_Pz_C3_b;
  abcd[326] = 4.0E0*I_NAI_Gx3y_Pz_C3_ab;
  abcd[327] = 4.0E0*I_NAI_Gx2yz_Pz_C3_ab;
  abcd[328] = 4.0E0*I_NAI_Gxy2z_Pz_C3_ab;
  abcd[329] = 4.0E0*I_NAI_Gx3z_Pz_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[330] = 4.0E0*I_NAI_G4x_Dxz_C1003_ab-2.0E0*3*I_NAI_D2x_Dxz_C1003_b;
  abcd[331] = 4.0E0*I_NAI_G3xy_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxz_C1003_b;
  abcd[332] = 4.0E0*I_NAI_G3xz_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxz_C1003_b;
  abcd[333] = 4.0E0*I_NAI_G2x2y_Dxz_C1003_ab-2.0E0*1*I_NAI_D2y_Dxz_C1003_b;
  abcd[334] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab-2.0E0*1*I_NAI_Dyz_Dxz_C1003_b;
  abcd[335] = 4.0E0*I_NAI_G2x2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2z_Dxz_C1003_b;
  abcd[336] = 4.0E0*I_NAI_Gx3y_Dxz_C1003_ab;
  abcd[337] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab;
  abcd[338] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab;
  abcd[339] = 4.0E0*I_NAI_Gx3z_Dxz_C1003_ab;
  abcd[340] = 4.0E0*I_NAI_G4x_Dyz_C1003_ab-2.0E0*3*I_NAI_D2x_Dyz_C1003_b;
  abcd[341] = 4.0E0*I_NAI_G3xy_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dyz_C1003_b;
  abcd[342] = 4.0E0*I_NAI_G3xz_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dyz_C1003_b;
  abcd[343] = 4.0E0*I_NAI_G2x2y_Dyz_C1003_ab-2.0E0*1*I_NAI_D2y_Dyz_C1003_b;
  abcd[344] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab-2.0E0*1*I_NAI_Dyz_Dyz_C1003_b;
  abcd[345] = 4.0E0*I_NAI_G2x2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2z_Dyz_C1003_b;
  abcd[346] = 4.0E0*I_NAI_Gx3y_Dyz_C1003_ab;
  abcd[347] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab;
  abcd[348] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab;
  abcd[349] = 4.0E0*I_NAI_Gx3z_Dyz_C1003_ab;
  abcd[350] = 4.0E0*I_NAI_G4x_D2z_C1003_ab-2.0E0*1*I_NAI_G4x_S_C1003_a-2.0E0*3*I_NAI_D2x_D2z_C1003_b+3*1*I_NAI_D2x_S_C1003;
  abcd[351] = 4.0E0*I_NAI_G3xy_D2z_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2z_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[352] = 4.0E0*I_NAI_G3xz_D2z_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2z_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[353] = 4.0E0*I_NAI_G2x2y_D2z_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2y_D2z_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[354] = 4.0E0*I_NAI_G2xyz_D2z_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a-2.0E0*1*I_NAI_Dyz_D2z_C1003_b+1*I_NAI_Dyz_S_C1003;
  abcd[355] = 4.0E0*I_NAI_G2x2z_D2z_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2z_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[356] = 4.0E0*I_NAI_Gx3y_D2z_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a;
  abcd[357] = 4.0E0*I_NAI_Gx2yz_D2z_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[358] = 4.0E0*I_NAI_Gxy2z_D2z_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[359] = 4.0E0*I_NAI_Gx3z_D2z_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[360] = 4.0E0*I_NAI_G3xy_Px_C3_ab;
  abcd[361] = 4.0E0*I_NAI_G2x2y_Px_C3_ab-2.0E0*1*I_NAI_D2x_Px_C3_b;
  abcd[362] = 4.0E0*I_NAI_G2xyz_Px_C3_ab;
  abcd[363] = 4.0E0*I_NAI_Gx3y_Px_C3_ab-2.0E0*2*I_NAI_Dxy_Px_C3_b;
  abcd[364] = 4.0E0*I_NAI_Gx2yz_Px_C3_ab-2.0E0*1*I_NAI_Dxz_Px_C3_b;
  abcd[365] = 4.0E0*I_NAI_Gxy2z_Px_C3_ab;
  abcd[366] = 4.0E0*I_NAI_G4y_Px_C3_ab-2.0E0*3*I_NAI_D2y_Px_C3_b;
  abcd[367] = 4.0E0*I_NAI_G3yz_Px_C3_ab-2.0E0*2*I_NAI_Dyz_Px_C3_b;
  abcd[368] = 4.0E0*I_NAI_G2y2z_Px_C3_ab-2.0E0*1*I_NAI_D2z_Px_C3_b;
  abcd[369] = 4.0E0*I_NAI_Gy3z_Px_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[370] = 4.0E0*I_NAI_G3xy_D2x_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a;
  abcd[371] = 4.0E0*I_NAI_G2x2y_D2x_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2x_D2x_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[372] = 4.0E0*I_NAI_G2xyz_D2x_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[373] = 4.0E0*I_NAI_Gx3y_D2x_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2x_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[374] = 4.0E0*I_NAI_Gx2yz_D2x_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a-2.0E0*1*I_NAI_Dxz_D2x_C1003_b+1*I_NAI_Dxz_S_C1003;
  abcd[375] = 4.0E0*I_NAI_Gxy2z_D2x_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[376] = 4.0E0*I_NAI_G4y_D2x_C1003_ab-2.0E0*1*I_NAI_G4y_S_C1003_a-2.0E0*3*I_NAI_D2y_D2x_C1003_b+3*1*I_NAI_D2y_S_C1003;
  abcd[377] = 4.0E0*I_NAI_G3yz_D2x_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2x_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[378] = 4.0E0*I_NAI_G2y2z_D2x_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2x_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[379] = 4.0E0*I_NAI_Gy3z_D2x_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a;
  abcd[380] = 4.0E0*I_NAI_G3xy_Dxy_C1003_ab;
  abcd[381] = 4.0E0*I_NAI_G2x2y_Dxy_C1003_ab-2.0E0*1*I_NAI_D2x_Dxy_C1003_b;
  abcd[382] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab;
  abcd[383] = 4.0E0*I_NAI_Gx3y_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxy_C1003_b;
  abcd[384] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab-2.0E0*1*I_NAI_Dxz_Dxy_C1003_b;
  abcd[385] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab;
  abcd[386] = 4.0E0*I_NAI_G4y_Dxy_C1003_ab-2.0E0*3*I_NAI_D2y_Dxy_C1003_b;
  abcd[387] = 4.0E0*I_NAI_G3yz_Dxy_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxy_C1003_b;
  abcd[388] = 4.0E0*I_NAI_G2y2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2z_Dxy_C1003_b;
  abcd[389] = 4.0E0*I_NAI_Gy3z_Dxy_C1003_ab;
  abcd[390] = 4.0E0*I_NAI_G3xy_Dxz_C1003_ab;
  abcd[391] = 4.0E0*I_NAI_G2x2y_Dxz_C1003_ab-2.0E0*1*I_NAI_D2x_Dxz_C1003_b;
  abcd[392] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab;
  abcd[393] = 4.0E0*I_NAI_Gx3y_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxz_C1003_b;
  abcd[394] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab-2.0E0*1*I_NAI_Dxz_Dxz_C1003_b;
  abcd[395] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab;
  abcd[396] = 4.0E0*I_NAI_G4y_Dxz_C1003_ab-2.0E0*3*I_NAI_D2y_Dxz_C1003_b;
  abcd[397] = 4.0E0*I_NAI_G3yz_Dxz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxz_C1003_b;
  abcd[398] = 4.0E0*I_NAI_G2y2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2z_Dxz_C1003_b;
  abcd[399] = 4.0E0*I_NAI_Gy3z_Dxz_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[400] = 4.0E0*I_NAI_G3xy_Py_C3_ab;
  abcd[401] = 4.0E0*I_NAI_G2x2y_Py_C3_ab-2.0E0*1*I_NAI_D2x_Py_C3_b;
  abcd[402] = 4.0E0*I_NAI_G2xyz_Py_C3_ab;
  abcd[403] = 4.0E0*I_NAI_Gx3y_Py_C3_ab-2.0E0*2*I_NAI_Dxy_Py_C3_b;
  abcd[404] = 4.0E0*I_NAI_Gx2yz_Py_C3_ab-2.0E0*1*I_NAI_Dxz_Py_C3_b;
  abcd[405] = 4.0E0*I_NAI_Gxy2z_Py_C3_ab;
  abcd[406] = 4.0E0*I_NAI_G4y_Py_C3_ab-2.0E0*3*I_NAI_D2y_Py_C3_b;
  abcd[407] = 4.0E0*I_NAI_G3yz_Py_C3_ab-2.0E0*2*I_NAI_Dyz_Py_C3_b;
  abcd[408] = 4.0E0*I_NAI_G2y2z_Py_C3_ab-2.0E0*1*I_NAI_D2z_Py_C3_b;
  abcd[409] = 4.0E0*I_NAI_Gy3z_Py_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[410] = 4.0E0*I_NAI_G3xy_Dxy_C1003_ab;
  abcd[411] = 4.0E0*I_NAI_G2x2y_Dxy_C1003_ab-2.0E0*1*I_NAI_D2x_Dxy_C1003_b;
  abcd[412] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab;
  abcd[413] = 4.0E0*I_NAI_Gx3y_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxy_C1003_b;
  abcd[414] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab-2.0E0*1*I_NAI_Dxz_Dxy_C1003_b;
  abcd[415] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab;
  abcd[416] = 4.0E0*I_NAI_G4y_Dxy_C1003_ab-2.0E0*3*I_NAI_D2y_Dxy_C1003_b;
  abcd[417] = 4.0E0*I_NAI_G3yz_Dxy_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxy_C1003_b;
  abcd[418] = 4.0E0*I_NAI_G2y2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2z_Dxy_C1003_b;
  abcd[419] = 4.0E0*I_NAI_Gy3z_Dxy_C1003_ab;
  abcd[420] = 4.0E0*I_NAI_G3xy_D2y_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a;
  abcd[421] = 4.0E0*I_NAI_G2x2y_D2y_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2x_D2y_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[422] = 4.0E0*I_NAI_G2xyz_D2y_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[423] = 4.0E0*I_NAI_Gx3y_D2y_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2y_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[424] = 4.0E0*I_NAI_Gx2yz_D2y_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a-2.0E0*1*I_NAI_Dxz_D2y_C1003_b+1*I_NAI_Dxz_S_C1003;
  abcd[425] = 4.0E0*I_NAI_Gxy2z_D2y_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[426] = 4.0E0*I_NAI_G4y_D2y_C1003_ab-2.0E0*1*I_NAI_G4y_S_C1003_a-2.0E0*3*I_NAI_D2y_D2y_C1003_b+3*1*I_NAI_D2y_S_C1003;
  abcd[427] = 4.0E0*I_NAI_G3yz_D2y_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2y_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[428] = 4.0E0*I_NAI_G2y2z_D2y_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2y_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[429] = 4.0E0*I_NAI_Gy3z_D2y_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a;
  abcd[430] = 4.0E0*I_NAI_G3xy_Dyz_C1003_ab;
  abcd[431] = 4.0E0*I_NAI_G2x2y_Dyz_C1003_ab-2.0E0*1*I_NAI_D2x_Dyz_C1003_b;
  abcd[432] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab;
  abcd[433] = 4.0E0*I_NAI_Gx3y_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dyz_C1003_b;
  abcd[434] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab-2.0E0*1*I_NAI_Dxz_Dyz_C1003_b;
  abcd[435] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab;
  abcd[436] = 4.0E0*I_NAI_G4y_Dyz_C1003_ab-2.0E0*3*I_NAI_D2y_Dyz_C1003_b;
  abcd[437] = 4.0E0*I_NAI_G3yz_Dyz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dyz_C1003_b;
  abcd[438] = 4.0E0*I_NAI_G2y2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2z_Dyz_C1003_b;
  abcd[439] = 4.0E0*I_NAI_Gy3z_Dyz_C1003_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[440] = 4.0E0*I_NAI_G3xy_Pz_C3_ab;
  abcd[441] = 4.0E0*I_NAI_G2x2y_Pz_C3_ab-2.0E0*1*I_NAI_D2x_Pz_C3_b;
  abcd[442] = 4.0E0*I_NAI_G2xyz_Pz_C3_ab;
  abcd[443] = 4.0E0*I_NAI_Gx3y_Pz_C3_ab-2.0E0*2*I_NAI_Dxy_Pz_C3_b;
  abcd[444] = 4.0E0*I_NAI_Gx2yz_Pz_C3_ab-2.0E0*1*I_NAI_Dxz_Pz_C3_b;
  abcd[445] = 4.0E0*I_NAI_Gxy2z_Pz_C3_ab;
  abcd[446] = 4.0E0*I_NAI_G4y_Pz_C3_ab-2.0E0*3*I_NAI_D2y_Pz_C3_b;
  abcd[447] = 4.0E0*I_NAI_G3yz_Pz_C3_ab-2.0E0*2*I_NAI_Dyz_Pz_C3_b;
  abcd[448] = 4.0E0*I_NAI_G2y2z_Pz_C3_ab-2.0E0*1*I_NAI_D2z_Pz_C3_b;
  abcd[449] = 4.0E0*I_NAI_Gy3z_Pz_C3_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[450] = 4.0E0*I_NAI_G3xy_Dxz_C1003_ab;
  abcd[451] = 4.0E0*I_NAI_G2x2y_Dxz_C1003_ab-2.0E0*1*I_NAI_D2x_Dxz_C1003_b;
  abcd[452] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab;
  abcd[453] = 4.0E0*I_NAI_Gx3y_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dxz_C1003_b;
  abcd[454] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab-2.0E0*1*I_NAI_Dxz_Dxz_C1003_b;
  abcd[455] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab;
  abcd[456] = 4.0E0*I_NAI_G4y_Dxz_C1003_ab-2.0E0*3*I_NAI_D2y_Dxz_C1003_b;
  abcd[457] = 4.0E0*I_NAI_G3yz_Dxz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxz_C1003_b;
  abcd[458] = 4.0E0*I_NAI_G2y2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2z_Dxz_C1003_b;
  abcd[459] = 4.0E0*I_NAI_Gy3z_Dxz_C1003_ab;
  abcd[460] = 4.0E0*I_NAI_G3xy_Dyz_C1003_ab;
  abcd[461] = 4.0E0*I_NAI_G2x2y_Dyz_C1003_ab-2.0E0*1*I_NAI_D2x_Dyz_C1003_b;
  abcd[462] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab;
  abcd[463] = 4.0E0*I_NAI_Gx3y_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxy_Dyz_C1003_b;
  abcd[464] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab-2.0E0*1*I_NAI_Dxz_Dyz_C1003_b;
  abcd[465] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab;
  abcd[466] = 4.0E0*I_NAI_G4y_Dyz_C1003_ab-2.0E0*3*I_NAI_D2y_Dyz_C1003_b;
  abcd[467] = 4.0E0*I_NAI_G3yz_Dyz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dyz_C1003_b;
  abcd[468] = 4.0E0*I_NAI_G2y2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2z_Dyz_C1003_b;
  abcd[469] = 4.0E0*I_NAI_Gy3z_Dyz_C1003_ab;
  abcd[470] = 4.0E0*I_NAI_G3xy_D2z_C1003_ab-2.0E0*1*I_NAI_G3xy_S_C1003_a;
  abcd[471] = 4.0E0*I_NAI_G2x2y_D2z_C1003_ab-2.0E0*1*I_NAI_G2x2y_S_C1003_a-2.0E0*1*I_NAI_D2x_D2z_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[472] = 4.0E0*I_NAI_G2xyz_D2z_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[473] = 4.0E0*I_NAI_Gx3y_D2z_C1003_ab-2.0E0*1*I_NAI_Gx3y_S_C1003_a-2.0E0*2*I_NAI_Dxy_D2z_C1003_b+2*1*I_NAI_Dxy_S_C1003;
  abcd[474] = 4.0E0*I_NAI_Gx2yz_D2z_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a-2.0E0*1*I_NAI_Dxz_D2z_C1003_b+1*I_NAI_Dxz_S_C1003;
  abcd[475] = 4.0E0*I_NAI_Gxy2z_D2z_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a;
  abcd[476] = 4.0E0*I_NAI_G4y_D2z_C1003_ab-2.0E0*1*I_NAI_G4y_S_C1003_a-2.0E0*3*I_NAI_D2y_D2z_C1003_b+3*1*I_NAI_D2y_S_C1003;
  abcd[477] = 4.0E0*I_NAI_G3yz_D2z_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2z_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[478] = 4.0E0*I_NAI_G2y2z_D2z_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2z_D2z_C1003_b+1*I_NAI_D2z_S_C1003;
  abcd[479] = 4.0E0*I_NAI_Gy3z_D2z_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[480] = 4.0E0*I_NAI_G3xz_Px_C3_ab;
  abcd[481] = 4.0E0*I_NAI_G2xyz_Px_C3_ab;
  abcd[482] = 4.0E0*I_NAI_G2x2z_Px_C3_ab-2.0E0*1*I_NAI_D2x_Px_C3_b;
  abcd[483] = 4.0E0*I_NAI_Gx2yz_Px_C3_ab;
  abcd[484] = 4.0E0*I_NAI_Gxy2z_Px_C3_ab-2.0E0*1*I_NAI_Dxy_Px_C3_b;
  abcd[485] = 4.0E0*I_NAI_Gx3z_Px_C3_ab-2.0E0*2*I_NAI_Dxz_Px_C3_b;
  abcd[486] = 4.0E0*I_NAI_G3yz_Px_C3_ab;
  abcd[487] = 4.0E0*I_NAI_G2y2z_Px_C3_ab-2.0E0*1*I_NAI_D2y_Px_C3_b;
  abcd[488] = 4.0E0*I_NAI_Gy3z_Px_C3_ab-2.0E0*2*I_NAI_Dyz_Px_C3_b;
  abcd[489] = 4.0E0*I_NAI_G4z_Px_C3_ab-2.0E0*3*I_NAI_D2z_Px_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[490] = 4.0E0*I_NAI_G3xz_D2x_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a;
  abcd[491] = 4.0E0*I_NAI_G2xyz_D2x_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[492] = 4.0E0*I_NAI_G2x2z_D2x_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2x_D2x_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[493] = 4.0E0*I_NAI_Gx2yz_D2x_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[494] = 4.0E0*I_NAI_Gxy2z_D2x_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a-2.0E0*1*I_NAI_Dxy_D2x_C1003_b+1*I_NAI_Dxy_S_C1003;
  abcd[495] = 4.0E0*I_NAI_Gx3z_D2x_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2x_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[496] = 4.0E0*I_NAI_G3yz_D2x_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a;
  abcd[497] = 4.0E0*I_NAI_G2y2z_D2x_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2y_D2x_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[498] = 4.0E0*I_NAI_Gy3z_D2x_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2x_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[499] = 4.0E0*I_NAI_G4z_D2x_C1003_ab-2.0E0*1*I_NAI_G4z_S_C1003_a-2.0E0*3*I_NAI_D2z_D2x_C1003_b+3*1*I_NAI_D2z_S_C1003;
  abcd[500] = 4.0E0*I_NAI_G3xz_Dxy_C1003_ab;
  abcd[501] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab;
  abcd[502] = 4.0E0*I_NAI_G2x2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2x_Dxy_C1003_b;
  abcd[503] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab;
  abcd[504] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab-2.0E0*1*I_NAI_Dxy_Dxy_C1003_b;
  abcd[505] = 4.0E0*I_NAI_Gx3z_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxy_C1003_b;
  abcd[506] = 4.0E0*I_NAI_G3yz_Dxy_C1003_ab;
  abcd[507] = 4.0E0*I_NAI_G2y2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2y_Dxy_C1003_b;
  abcd[508] = 4.0E0*I_NAI_Gy3z_Dxy_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxy_C1003_b;
  abcd[509] = 4.0E0*I_NAI_G4z_Dxy_C1003_ab-2.0E0*3*I_NAI_D2z_Dxy_C1003_b;
  abcd[510] = 4.0E0*I_NAI_G3xz_Dxz_C1003_ab;
  abcd[511] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab;
  abcd[512] = 4.0E0*I_NAI_G2x2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2x_Dxz_C1003_b;
  abcd[513] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab;
  abcd[514] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab-2.0E0*1*I_NAI_Dxy_Dxz_C1003_b;
  abcd[515] = 4.0E0*I_NAI_Gx3z_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxz_C1003_b;
  abcd[516] = 4.0E0*I_NAI_G3yz_Dxz_C1003_ab;
  abcd[517] = 4.0E0*I_NAI_G2y2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2y_Dxz_C1003_b;
  abcd[518] = 4.0E0*I_NAI_Gy3z_Dxz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxz_C1003_b;
  abcd[519] = 4.0E0*I_NAI_G4z_Dxz_C1003_ab-2.0E0*3*I_NAI_D2z_Dxz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[520] = 4.0E0*I_NAI_G3xz_Py_C3_ab;
  abcd[521] = 4.0E0*I_NAI_G2xyz_Py_C3_ab;
  abcd[522] = 4.0E0*I_NAI_G2x2z_Py_C3_ab-2.0E0*1*I_NAI_D2x_Py_C3_b;
  abcd[523] = 4.0E0*I_NAI_Gx2yz_Py_C3_ab;
  abcd[524] = 4.0E0*I_NAI_Gxy2z_Py_C3_ab-2.0E0*1*I_NAI_Dxy_Py_C3_b;
  abcd[525] = 4.0E0*I_NAI_Gx3z_Py_C3_ab-2.0E0*2*I_NAI_Dxz_Py_C3_b;
  abcd[526] = 4.0E0*I_NAI_G3yz_Py_C3_ab;
  abcd[527] = 4.0E0*I_NAI_G2y2z_Py_C3_ab-2.0E0*1*I_NAI_D2y_Py_C3_b;
  abcd[528] = 4.0E0*I_NAI_Gy3z_Py_C3_ab-2.0E0*2*I_NAI_Dyz_Py_C3_b;
  abcd[529] = 4.0E0*I_NAI_G4z_Py_C3_ab-2.0E0*3*I_NAI_D2z_Py_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[530] = 4.0E0*I_NAI_G3xz_Dxy_C1003_ab;
  abcd[531] = 4.0E0*I_NAI_G2xyz_Dxy_C1003_ab;
  abcd[532] = 4.0E0*I_NAI_G2x2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2x_Dxy_C1003_b;
  abcd[533] = 4.0E0*I_NAI_Gx2yz_Dxy_C1003_ab;
  abcd[534] = 4.0E0*I_NAI_Gxy2z_Dxy_C1003_ab-2.0E0*1*I_NAI_Dxy_Dxy_C1003_b;
  abcd[535] = 4.0E0*I_NAI_Gx3z_Dxy_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxy_C1003_b;
  abcd[536] = 4.0E0*I_NAI_G3yz_Dxy_C1003_ab;
  abcd[537] = 4.0E0*I_NAI_G2y2z_Dxy_C1003_ab-2.0E0*1*I_NAI_D2y_Dxy_C1003_b;
  abcd[538] = 4.0E0*I_NAI_Gy3z_Dxy_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxy_C1003_b;
  abcd[539] = 4.0E0*I_NAI_G4z_Dxy_C1003_ab-2.0E0*3*I_NAI_D2z_Dxy_C1003_b;
  abcd[540] = 4.0E0*I_NAI_G3xz_D2y_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a;
  abcd[541] = 4.0E0*I_NAI_G2xyz_D2y_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[542] = 4.0E0*I_NAI_G2x2z_D2y_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2x_D2y_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[543] = 4.0E0*I_NAI_Gx2yz_D2y_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[544] = 4.0E0*I_NAI_Gxy2z_D2y_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a-2.0E0*1*I_NAI_Dxy_D2y_C1003_b+1*I_NAI_Dxy_S_C1003;
  abcd[545] = 4.0E0*I_NAI_Gx3z_D2y_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2y_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[546] = 4.0E0*I_NAI_G3yz_D2y_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a;
  abcd[547] = 4.0E0*I_NAI_G2y2z_D2y_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2y_D2y_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[548] = 4.0E0*I_NAI_Gy3z_D2y_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2y_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[549] = 4.0E0*I_NAI_G4z_D2y_C1003_ab-2.0E0*1*I_NAI_G4z_S_C1003_a-2.0E0*3*I_NAI_D2z_D2y_C1003_b+3*1*I_NAI_D2z_S_C1003;
  abcd[550] = 4.0E0*I_NAI_G3xz_Dyz_C1003_ab;
  abcd[551] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab;
  abcd[552] = 4.0E0*I_NAI_G2x2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2x_Dyz_C1003_b;
  abcd[553] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab;
  abcd[554] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab-2.0E0*1*I_NAI_Dxy_Dyz_C1003_b;
  abcd[555] = 4.0E0*I_NAI_Gx3z_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dyz_C1003_b;
  abcd[556] = 4.0E0*I_NAI_G3yz_Dyz_C1003_ab;
  abcd[557] = 4.0E0*I_NAI_G2y2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2y_Dyz_C1003_b;
  abcd[558] = 4.0E0*I_NAI_Gy3z_Dyz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dyz_C1003_b;
  abcd[559] = 4.0E0*I_NAI_G4z_Dyz_C1003_ab-2.0E0*3*I_NAI_D2z_Dyz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C3_ab
   * RHS shell quartet name: SQ_NAI_D_P_C3_b
   ************************************************************/
  abcd[560] = 4.0E0*I_NAI_G3xz_Pz_C3_ab;
  abcd[561] = 4.0E0*I_NAI_G2xyz_Pz_C3_ab;
  abcd[562] = 4.0E0*I_NAI_G2x2z_Pz_C3_ab-2.0E0*1*I_NAI_D2x_Pz_C3_b;
  abcd[563] = 4.0E0*I_NAI_Gx2yz_Pz_C3_ab;
  abcd[564] = 4.0E0*I_NAI_Gxy2z_Pz_C3_ab-2.0E0*1*I_NAI_Dxy_Pz_C3_b;
  abcd[565] = 4.0E0*I_NAI_Gx3z_Pz_C3_ab-2.0E0*2*I_NAI_Dxz_Pz_C3_b;
  abcd[566] = 4.0E0*I_NAI_G3yz_Pz_C3_ab;
  abcd[567] = 4.0E0*I_NAI_G2y2z_Pz_C3_ab-2.0E0*1*I_NAI_D2y_Pz_C3_b;
  abcd[568] = 4.0E0*I_NAI_Gy3z_Pz_C3_ab-2.0E0*2*I_NAI_Dyz_Pz_C3_b;
  abcd[569] = 4.0E0*I_NAI_G4z_Pz_C3_ab-2.0E0*3*I_NAI_D2z_Pz_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_C1003_ab
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   * RHS shell quartet name: SQ_NAI_D_D_C1003_b
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  abcd[570] = 4.0E0*I_NAI_G3xz_Dxz_C1003_ab;
  abcd[571] = 4.0E0*I_NAI_G2xyz_Dxz_C1003_ab;
  abcd[572] = 4.0E0*I_NAI_G2x2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2x_Dxz_C1003_b;
  abcd[573] = 4.0E0*I_NAI_Gx2yz_Dxz_C1003_ab;
  abcd[574] = 4.0E0*I_NAI_Gxy2z_Dxz_C1003_ab-2.0E0*1*I_NAI_Dxy_Dxz_C1003_b;
  abcd[575] = 4.0E0*I_NAI_Gx3z_Dxz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dxz_C1003_b;
  abcd[576] = 4.0E0*I_NAI_G3yz_Dxz_C1003_ab;
  abcd[577] = 4.0E0*I_NAI_G2y2z_Dxz_C1003_ab-2.0E0*1*I_NAI_D2y_Dxz_C1003_b;
  abcd[578] = 4.0E0*I_NAI_Gy3z_Dxz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dxz_C1003_b;
  abcd[579] = 4.0E0*I_NAI_G4z_Dxz_C1003_ab-2.0E0*3*I_NAI_D2z_Dxz_C1003_b;
  abcd[580] = 4.0E0*I_NAI_G3xz_Dyz_C1003_ab;
  abcd[581] = 4.0E0*I_NAI_G2xyz_Dyz_C1003_ab;
  abcd[582] = 4.0E0*I_NAI_G2x2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2x_Dyz_C1003_b;
  abcd[583] = 4.0E0*I_NAI_Gx2yz_Dyz_C1003_ab;
  abcd[584] = 4.0E0*I_NAI_Gxy2z_Dyz_C1003_ab-2.0E0*1*I_NAI_Dxy_Dyz_C1003_b;
  abcd[585] = 4.0E0*I_NAI_Gx3z_Dyz_C1003_ab-2.0E0*2*I_NAI_Dxz_Dyz_C1003_b;
  abcd[586] = 4.0E0*I_NAI_G3yz_Dyz_C1003_ab;
  abcd[587] = 4.0E0*I_NAI_G2y2z_Dyz_C1003_ab-2.0E0*1*I_NAI_D2y_Dyz_C1003_b;
  abcd[588] = 4.0E0*I_NAI_Gy3z_Dyz_C1003_ab-2.0E0*2*I_NAI_Dyz_Dyz_C1003_b;
  abcd[589] = 4.0E0*I_NAI_G4z_Dyz_C1003_ab-2.0E0*3*I_NAI_D2z_Dyz_C1003_b;
  abcd[590] = 4.0E0*I_NAI_G3xz_D2z_C1003_ab-2.0E0*1*I_NAI_G3xz_S_C1003_a;
  abcd[591] = 4.0E0*I_NAI_G2xyz_D2z_C1003_ab-2.0E0*1*I_NAI_G2xyz_S_C1003_a;
  abcd[592] = 4.0E0*I_NAI_G2x2z_D2z_C1003_ab-2.0E0*1*I_NAI_G2x2z_S_C1003_a-2.0E0*1*I_NAI_D2x_D2z_C1003_b+1*I_NAI_D2x_S_C1003;
  abcd[593] = 4.0E0*I_NAI_Gx2yz_D2z_C1003_ab-2.0E0*1*I_NAI_Gx2yz_S_C1003_a;
  abcd[594] = 4.0E0*I_NAI_Gxy2z_D2z_C1003_ab-2.0E0*1*I_NAI_Gxy2z_S_C1003_a-2.0E0*1*I_NAI_Dxy_D2z_C1003_b+1*I_NAI_Dxy_S_C1003;
  abcd[595] = 4.0E0*I_NAI_Gx3z_D2z_C1003_ab-2.0E0*1*I_NAI_Gx3z_S_C1003_a-2.0E0*2*I_NAI_Dxz_D2z_C1003_b+2*1*I_NAI_Dxz_S_C1003;
  abcd[596] = 4.0E0*I_NAI_G3yz_D2z_C1003_ab-2.0E0*1*I_NAI_G3yz_S_C1003_a;
  abcd[597] = 4.0E0*I_NAI_G2y2z_D2z_C1003_ab-2.0E0*1*I_NAI_G2y2z_S_C1003_a-2.0E0*1*I_NAI_D2y_D2z_C1003_b+1*I_NAI_D2y_S_C1003;
  abcd[598] = 4.0E0*I_NAI_Gy3z_D2z_C1003_ab-2.0E0*1*I_NAI_Gy3z_S_C1003_a-2.0E0*2*I_NAI_Dyz_D2z_C1003_b+2*1*I_NAI_Dyz_S_C1003;
  abcd[599] = 4.0E0*I_NAI_G4z_D2z_C1003_ab-2.0E0*1*I_NAI_G4z_S_C1003_a-2.0E0*3*I_NAI_D2z_D2z_C1003_b+3*1*I_NAI_D2z_S_C1003;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[600] = 4.0E0*I_NAI_F3x_D2x_C3_bb-2.0E0*1*I_NAI_F3x_S_C3_b;
  abcd[601] = 4.0E0*I_NAI_F2xy_D2x_C3_bb-2.0E0*1*I_NAI_F2xy_S_C3_b;
  abcd[602] = 4.0E0*I_NAI_F2xz_D2x_C3_bb-2.0E0*1*I_NAI_F2xz_S_C3_b;
  abcd[603] = 4.0E0*I_NAI_Fx2y_D2x_C3_bb-2.0E0*1*I_NAI_Fx2y_S_C3_b;
  abcd[604] = 4.0E0*I_NAI_Fxyz_D2x_C3_bb-2.0E0*1*I_NAI_Fxyz_S_C3_b;
  abcd[605] = 4.0E0*I_NAI_Fx2z_D2x_C3_bb-2.0E0*1*I_NAI_Fx2z_S_C3_b;
  abcd[606] = 4.0E0*I_NAI_F3y_D2x_C3_bb-2.0E0*1*I_NAI_F3y_S_C3_b;
  abcd[607] = 4.0E0*I_NAI_F2yz_D2x_C3_bb-2.0E0*1*I_NAI_F2yz_S_C3_b;
  abcd[608] = 4.0E0*I_NAI_Fy2z_D2x_C3_bb-2.0E0*1*I_NAI_Fy2z_S_C3_b;
  abcd[609] = 4.0E0*I_NAI_F3z_D2x_C3_bb-2.0E0*1*I_NAI_F3z_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[610] = 4.0E0*I_NAI_F3x_F3x_C1003_bb-2.0E0*1*I_NAI_F3x_Px_C1003_b-2.0E0*2*I_NAI_F3x_Px_C1003_b;
  abcd[611] = 4.0E0*I_NAI_F2xy_F3x_C1003_bb-2.0E0*1*I_NAI_F2xy_Px_C1003_b-2.0E0*2*I_NAI_F2xy_Px_C1003_b;
  abcd[612] = 4.0E0*I_NAI_F2xz_F3x_C1003_bb-2.0E0*1*I_NAI_F2xz_Px_C1003_b-2.0E0*2*I_NAI_F2xz_Px_C1003_b;
  abcd[613] = 4.0E0*I_NAI_Fx2y_F3x_C1003_bb-2.0E0*1*I_NAI_Fx2y_Px_C1003_b-2.0E0*2*I_NAI_Fx2y_Px_C1003_b;
  abcd[614] = 4.0E0*I_NAI_Fxyz_F3x_C1003_bb-2.0E0*1*I_NAI_Fxyz_Px_C1003_b-2.0E0*2*I_NAI_Fxyz_Px_C1003_b;
  abcd[615] = 4.0E0*I_NAI_Fx2z_F3x_C1003_bb-2.0E0*1*I_NAI_Fx2z_Px_C1003_b-2.0E0*2*I_NAI_Fx2z_Px_C1003_b;
  abcd[616] = 4.0E0*I_NAI_F3y_F3x_C1003_bb-2.0E0*1*I_NAI_F3y_Px_C1003_b-2.0E0*2*I_NAI_F3y_Px_C1003_b;
  abcd[617] = 4.0E0*I_NAI_F2yz_F3x_C1003_bb-2.0E0*1*I_NAI_F2yz_Px_C1003_b-2.0E0*2*I_NAI_F2yz_Px_C1003_b;
  abcd[618] = 4.0E0*I_NAI_Fy2z_F3x_C1003_bb-2.0E0*1*I_NAI_Fy2z_Px_C1003_b-2.0E0*2*I_NAI_Fy2z_Px_C1003_b;
  abcd[619] = 4.0E0*I_NAI_F3z_F3x_C1003_bb-2.0E0*1*I_NAI_F3z_Px_C1003_b-2.0E0*2*I_NAI_F3z_Px_C1003_b;
  abcd[620] = 4.0E0*I_NAI_F3x_F2xy_C1003_bb-2.0E0*1*I_NAI_F3x_Py_C1003_b;
  abcd[621] = 4.0E0*I_NAI_F2xy_F2xy_C1003_bb-2.0E0*1*I_NAI_F2xy_Py_C1003_b;
  abcd[622] = 4.0E0*I_NAI_F2xz_F2xy_C1003_bb-2.0E0*1*I_NAI_F2xz_Py_C1003_b;
  abcd[623] = 4.0E0*I_NAI_Fx2y_F2xy_C1003_bb-2.0E0*1*I_NAI_Fx2y_Py_C1003_b;
  abcd[624] = 4.0E0*I_NAI_Fxyz_F2xy_C1003_bb-2.0E0*1*I_NAI_Fxyz_Py_C1003_b;
  abcd[625] = 4.0E0*I_NAI_Fx2z_F2xy_C1003_bb-2.0E0*1*I_NAI_Fx2z_Py_C1003_b;
  abcd[626] = 4.0E0*I_NAI_F3y_F2xy_C1003_bb-2.0E0*1*I_NAI_F3y_Py_C1003_b;
  abcd[627] = 4.0E0*I_NAI_F2yz_F2xy_C1003_bb-2.0E0*1*I_NAI_F2yz_Py_C1003_b;
  abcd[628] = 4.0E0*I_NAI_Fy2z_F2xy_C1003_bb-2.0E0*1*I_NAI_Fy2z_Py_C1003_b;
  abcd[629] = 4.0E0*I_NAI_F3z_F2xy_C1003_bb-2.0E0*1*I_NAI_F3z_Py_C1003_b;
  abcd[630] = 4.0E0*I_NAI_F3x_F2xz_C1003_bb-2.0E0*1*I_NAI_F3x_Pz_C1003_b;
  abcd[631] = 4.0E0*I_NAI_F2xy_F2xz_C1003_bb-2.0E0*1*I_NAI_F2xy_Pz_C1003_b;
  abcd[632] = 4.0E0*I_NAI_F2xz_F2xz_C1003_bb-2.0E0*1*I_NAI_F2xz_Pz_C1003_b;
  abcd[633] = 4.0E0*I_NAI_Fx2y_F2xz_C1003_bb-2.0E0*1*I_NAI_Fx2y_Pz_C1003_b;
  abcd[634] = 4.0E0*I_NAI_Fxyz_F2xz_C1003_bb-2.0E0*1*I_NAI_Fxyz_Pz_C1003_b;
  abcd[635] = 4.0E0*I_NAI_Fx2z_F2xz_C1003_bb-2.0E0*1*I_NAI_Fx2z_Pz_C1003_b;
  abcd[636] = 4.0E0*I_NAI_F3y_F2xz_C1003_bb-2.0E0*1*I_NAI_F3y_Pz_C1003_b;
  abcd[637] = 4.0E0*I_NAI_F2yz_F2xz_C1003_bb-2.0E0*1*I_NAI_F2yz_Pz_C1003_b;
  abcd[638] = 4.0E0*I_NAI_Fy2z_F2xz_C1003_bb-2.0E0*1*I_NAI_Fy2z_Pz_C1003_b;
  abcd[639] = 4.0E0*I_NAI_F3z_F2xz_C1003_bb-2.0E0*1*I_NAI_F3z_Pz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[640] = 4.0E0*I_NAI_F3x_Dxy_C3_bb;
  abcd[641] = 4.0E0*I_NAI_F2xy_Dxy_C3_bb;
  abcd[642] = 4.0E0*I_NAI_F2xz_Dxy_C3_bb;
  abcd[643] = 4.0E0*I_NAI_Fx2y_Dxy_C3_bb;
  abcd[644] = 4.0E0*I_NAI_Fxyz_Dxy_C3_bb;
  abcd[645] = 4.0E0*I_NAI_Fx2z_Dxy_C3_bb;
  abcd[646] = 4.0E0*I_NAI_F3y_Dxy_C3_bb;
  abcd[647] = 4.0E0*I_NAI_F2yz_Dxy_C3_bb;
  abcd[648] = 4.0E0*I_NAI_Fy2z_Dxy_C3_bb;
  abcd[649] = 4.0E0*I_NAI_F3z_Dxy_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[650] = 4.0E0*I_NAI_F3x_F2xy_C1003_bb-2.0E0*1*I_NAI_F3x_Py_C1003_b;
  abcd[651] = 4.0E0*I_NAI_F2xy_F2xy_C1003_bb-2.0E0*1*I_NAI_F2xy_Py_C1003_b;
  abcd[652] = 4.0E0*I_NAI_F2xz_F2xy_C1003_bb-2.0E0*1*I_NAI_F2xz_Py_C1003_b;
  abcd[653] = 4.0E0*I_NAI_Fx2y_F2xy_C1003_bb-2.0E0*1*I_NAI_Fx2y_Py_C1003_b;
  abcd[654] = 4.0E0*I_NAI_Fxyz_F2xy_C1003_bb-2.0E0*1*I_NAI_Fxyz_Py_C1003_b;
  abcd[655] = 4.0E0*I_NAI_Fx2z_F2xy_C1003_bb-2.0E0*1*I_NAI_Fx2z_Py_C1003_b;
  abcd[656] = 4.0E0*I_NAI_F3y_F2xy_C1003_bb-2.0E0*1*I_NAI_F3y_Py_C1003_b;
  abcd[657] = 4.0E0*I_NAI_F2yz_F2xy_C1003_bb-2.0E0*1*I_NAI_F2yz_Py_C1003_b;
  abcd[658] = 4.0E0*I_NAI_Fy2z_F2xy_C1003_bb-2.0E0*1*I_NAI_Fy2z_Py_C1003_b;
  abcd[659] = 4.0E0*I_NAI_F3z_F2xy_C1003_bb-2.0E0*1*I_NAI_F3z_Py_C1003_b;
  abcd[660] = 4.0E0*I_NAI_F3x_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3x_Px_C1003_b;
  abcd[661] = 4.0E0*I_NAI_F2xy_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2xy_Px_C1003_b;
  abcd[662] = 4.0E0*I_NAI_F2xz_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2xz_Px_C1003_b;
  abcd[663] = 4.0E0*I_NAI_Fx2y_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fx2y_Px_C1003_b;
  abcd[664] = 4.0E0*I_NAI_Fxyz_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fxyz_Px_C1003_b;
  abcd[665] = 4.0E0*I_NAI_Fx2z_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fx2z_Px_C1003_b;
  abcd[666] = 4.0E0*I_NAI_F3y_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3y_Px_C1003_b;
  abcd[667] = 4.0E0*I_NAI_F2yz_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2yz_Px_C1003_b;
  abcd[668] = 4.0E0*I_NAI_Fy2z_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fy2z_Px_C1003_b;
  abcd[669] = 4.0E0*I_NAI_F3z_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3z_Px_C1003_b;
  abcd[670] = 4.0E0*I_NAI_F3x_Fxyz_C1003_bb;
  abcd[671] = 4.0E0*I_NAI_F2xy_Fxyz_C1003_bb;
  abcd[672] = 4.0E0*I_NAI_F2xz_Fxyz_C1003_bb;
  abcd[673] = 4.0E0*I_NAI_Fx2y_Fxyz_C1003_bb;
  abcd[674] = 4.0E0*I_NAI_Fxyz_Fxyz_C1003_bb;
  abcd[675] = 4.0E0*I_NAI_Fx2z_Fxyz_C1003_bb;
  abcd[676] = 4.0E0*I_NAI_F3y_Fxyz_C1003_bb;
  abcd[677] = 4.0E0*I_NAI_F2yz_Fxyz_C1003_bb;
  abcd[678] = 4.0E0*I_NAI_Fy2z_Fxyz_C1003_bb;
  abcd[679] = 4.0E0*I_NAI_F3z_Fxyz_C1003_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[680] = 4.0E0*I_NAI_F3x_Dxz_C3_bb;
  abcd[681] = 4.0E0*I_NAI_F2xy_Dxz_C3_bb;
  abcd[682] = 4.0E0*I_NAI_F2xz_Dxz_C3_bb;
  abcd[683] = 4.0E0*I_NAI_Fx2y_Dxz_C3_bb;
  abcd[684] = 4.0E0*I_NAI_Fxyz_Dxz_C3_bb;
  abcd[685] = 4.0E0*I_NAI_Fx2z_Dxz_C3_bb;
  abcd[686] = 4.0E0*I_NAI_F3y_Dxz_C3_bb;
  abcd[687] = 4.0E0*I_NAI_F2yz_Dxz_C3_bb;
  abcd[688] = 4.0E0*I_NAI_Fy2z_Dxz_C3_bb;
  abcd[689] = 4.0E0*I_NAI_F3z_Dxz_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[690] = 4.0E0*I_NAI_F3x_F2xz_C1003_bb-2.0E0*1*I_NAI_F3x_Pz_C1003_b;
  abcd[691] = 4.0E0*I_NAI_F2xy_F2xz_C1003_bb-2.0E0*1*I_NAI_F2xy_Pz_C1003_b;
  abcd[692] = 4.0E0*I_NAI_F2xz_F2xz_C1003_bb-2.0E0*1*I_NAI_F2xz_Pz_C1003_b;
  abcd[693] = 4.0E0*I_NAI_Fx2y_F2xz_C1003_bb-2.0E0*1*I_NAI_Fx2y_Pz_C1003_b;
  abcd[694] = 4.0E0*I_NAI_Fxyz_F2xz_C1003_bb-2.0E0*1*I_NAI_Fxyz_Pz_C1003_b;
  abcd[695] = 4.0E0*I_NAI_Fx2z_F2xz_C1003_bb-2.0E0*1*I_NAI_Fx2z_Pz_C1003_b;
  abcd[696] = 4.0E0*I_NAI_F3y_F2xz_C1003_bb-2.0E0*1*I_NAI_F3y_Pz_C1003_b;
  abcd[697] = 4.0E0*I_NAI_F2yz_F2xz_C1003_bb-2.0E0*1*I_NAI_F2yz_Pz_C1003_b;
  abcd[698] = 4.0E0*I_NAI_Fy2z_F2xz_C1003_bb-2.0E0*1*I_NAI_Fy2z_Pz_C1003_b;
  abcd[699] = 4.0E0*I_NAI_F3z_F2xz_C1003_bb-2.0E0*1*I_NAI_F3z_Pz_C1003_b;
  abcd[700] = 4.0E0*I_NAI_F3x_Fxyz_C1003_bb;
  abcd[701] = 4.0E0*I_NAI_F2xy_Fxyz_C1003_bb;
  abcd[702] = 4.0E0*I_NAI_F2xz_Fxyz_C1003_bb;
  abcd[703] = 4.0E0*I_NAI_Fx2y_Fxyz_C1003_bb;
  abcd[704] = 4.0E0*I_NAI_Fxyz_Fxyz_C1003_bb;
  abcd[705] = 4.0E0*I_NAI_Fx2z_Fxyz_C1003_bb;
  abcd[706] = 4.0E0*I_NAI_F3y_Fxyz_C1003_bb;
  abcd[707] = 4.0E0*I_NAI_F2yz_Fxyz_C1003_bb;
  abcd[708] = 4.0E0*I_NAI_Fy2z_Fxyz_C1003_bb;
  abcd[709] = 4.0E0*I_NAI_F3z_Fxyz_C1003_bb;
  abcd[710] = 4.0E0*I_NAI_F3x_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3x_Px_C1003_b;
  abcd[711] = 4.0E0*I_NAI_F2xy_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2xy_Px_C1003_b;
  abcd[712] = 4.0E0*I_NAI_F2xz_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2xz_Px_C1003_b;
  abcd[713] = 4.0E0*I_NAI_Fx2y_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fx2y_Px_C1003_b;
  abcd[714] = 4.0E0*I_NAI_Fxyz_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fxyz_Px_C1003_b;
  abcd[715] = 4.0E0*I_NAI_Fx2z_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fx2z_Px_C1003_b;
  abcd[716] = 4.0E0*I_NAI_F3y_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3y_Px_C1003_b;
  abcd[717] = 4.0E0*I_NAI_F2yz_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2yz_Px_C1003_b;
  abcd[718] = 4.0E0*I_NAI_Fy2z_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fy2z_Px_C1003_b;
  abcd[719] = 4.0E0*I_NAI_F3z_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3z_Px_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[720] = 4.0E0*I_NAI_F3x_D2y_C3_bb-2.0E0*1*I_NAI_F3x_S_C3_b;
  abcd[721] = 4.0E0*I_NAI_F2xy_D2y_C3_bb-2.0E0*1*I_NAI_F2xy_S_C3_b;
  abcd[722] = 4.0E0*I_NAI_F2xz_D2y_C3_bb-2.0E0*1*I_NAI_F2xz_S_C3_b;
  abcd[723] = 4.0E0*I_NAI_Fx2y_D2y_C3_bb-2.0E0*1*I_NAI_Fx2y_S_C3_b;
  abcd[724] = 4.0E0*I_NAI_Fxyz_D2y_C3_bb-2.0E0*1*I_NAI_Fxyz_S_C3_b;
  abcd[725] = 4.0E0*I_NAI_Fx2z_D2y_C3_bb-2.0E0*1*I_NAI_Fx2z_S_C3_b;
  abcd[726] = 4.0E0*I_NAI_F3y_D2y_C3_bb-2.0E0*1*I_NAI_F3y_S_C3_b;
  abcd[727] = 4.0E0*I_NAI_F2yz_D2y_C3_bb-2.0E0*1*I_NAI_F2yz_S_C3_b;
  abcd[728] = 4.0E0*I_NAI_Fy2z_D2y_C3_bb-2.0E0*1*I_NAI_Fy2z_S_C3_b;
  abcd[729] = 4.0E0*I_NAI_F3z_D2y_C3_bb-2.0E0*1*I_NAI_F3z_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[730] = 4.0E0*I_NAI_F3x_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3x_Px_C1003_b;
  abcd[731] = 4.0E0*I_NAI_F2xy_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2xy_Px_C1003_b;
  abcd[732] = 4.0E0*I_NAI_F2xz_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2xz_Px_C1003_b;
  abcd[733] = 4.0E0*I_NAI_Fx2y_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fx2y_Px_C1003_b;
  abcd[734] = 4.0E0*I_NAI_Fxyz_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fxyz_Px_C1003_b;
  abcd[735] = 4.0E0*I_NAI_Fx2z_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fx2z_Px_C1003_b;
  abcd[736] = 4.0E0*I_NAI_F3y_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3y_Px_C1003_b;
  abcd[737] = 4.0E0*I_NAI_F2yz_Fx2y_C1003_bb-2.0E0*1*I_NAI_F2yz_Px_C1003_b;
  abcd[738] = 4.0E0*I_NAI_Fy2z_Fx2y_C1003_bb-2.0E0*1*I_NAI_Fy2z_Px_C1003_b;
  abcd[739] = 4.0E0*I_NAI_F3z_Fx2y_C1003_bb-2.0E0*1*I_NAI_F3z_Px_C1003_b;
  abcd[740] = 4.0E0*I_NAI_F3x_F3y_C1003_bb-2.0E0*1*I_NAI_F3x_Py_C1003_b-2.0E0*2*I_NAI_F3x_Py_C1003_b;
  abcd[741] = 4.0E0*I_NAI_F2xy_F3y_C1003_bb-2.0E0*1*I_NAI_F2xy_Py_C1003_b-2.0E0*2*I_NAI_F2xy_Py_C1003_b;
  abcd[742] = 4.0E0*I_NAI_F2xz_F3y_C1003_bb-2.0E0*1*I_NAI_F2xz_Py_C1003_b-2.0E0*2*I_NAI_F2xz_Py_C1003_b;
  abcd[743] = 4.0E0*I_NAI_Fx2y_F3y_C1003_bb-2.0E0*1*I_NAI_Fx2y_Py_C1003_b-2.0E0*2*I_NAI_Fx2y_Py_C1003_b;
  abcd[744] = 4.0E0*I_NAI_Fxyz_F3y_C1003_bb-2.0E0*1*I_NAI_Fxyz_Py_C1003_b-2.0E0*2*I_NAI_Fxyz_Py_C1003_b;
  abcd[745] = 4.0E0*I_NAI_Fx2z_F3y_C1003_bb-2.0E0*1*I_NAI_Fx2z_Py_C1003_b-2.0E0*2*I_NAI_Fx2z_Py_C1003_b;
  abcd[746] = 4.0E0*I_NAI_F3y_F3y_C1003_bb-2.0E0*1*I_NAI_F3y_Py_C1003_b-2.0E0*2*I_NAI_F3y_Py_C1003_b;
  abcd[747] = 4.0E0*I_NAI_F2yz_F3y_C1003_bb-2.0E0*1*I_NAI_F2yz_Py_C1003_b-2.0E0*2*I_NAI_F2yz_Py_C1003_b;
  abcd[748] = 4.0E0*I_NAI_Fy2z_F3y_C1003_bb-2.0E0*1*I_NAI_Fy2z_Py_C1003_b-2.0E0*2*I_NAI_Fy2z_Py_C1003_b;
  abcd[749] = 4.0E0*I_NAI_F3z_F3y_C1003_bb-2.0E0*1*I_NAI_F3z_Py_C1003_b-2.0E0*2*I_NAI_F3z_Py_C1003_b;
  abcd[750] = 4.0E0*I_NAI_F3x_F2yz_C1003_bb-2.0E0*1*I_NAI_F3x_Pz_C1003_b;
  abcd[751] = 4.0E0*I_NAI_F2xy_F2yz_C1003_bb-2.0E0*1*I_NAI_F2xy_Pz_C1003_b;
  abcd[752] = 4.0E0*I_NAI_F2xz_F2yz_C1003_bb-2.0E0*1*I_NAI_F2xz_Pz_C1003_b;
  abcd[753] = 4.0E0*I_NAI_Fx2y_F2yz_C1003_bb-2.0E0*1*I_NAI_Fx2y_Pz_C1003_b;
  abcd[754] = 4.0E0*I_NAI_Fxyz_F2yz_C1003_bb-2.0E0*1*I_NAI_Fxyz_Pz_C1003_b;
  abcd[755] = 4.0E0*I_NAI_Fx2z_F2yz_C1003_bb-2.0E0*1*I_NAI_Fx2z_Pz_C1003_b;
  abcd[756] = 4.0E0*I_NAI_F3y_F2yz_C1003_bb-2.0E0*1*I_NAI_F3y_Pz_C1003_b;
  abcd[757] = 4.0E0*I_NAI_F2yz_F2yz_C1003_bb-2.0E0*1*I_NAI_F2yz_Pz_C1003_b;
  abcd[758] = 4.0E0*I_NAI_Fy2z_F2yz_C1003_bb-2.0E0*1*I_NAI_Fy2z_Pz_C1003_b;
  abcd[759] = 4.0E0*I_NAI_F3z_F2yz_C1003_bb-2.0E0*1*I_NAI_F3z_Pz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[760] = 4.0E0*I_NAI_F3x_Dyz_C3_bb;
  abcd[761] = 4.0E0*I_NAI_F2xy_Dyz_C3_bb;
  abcd[762] = 4.0E0*I_NAI_F2xz_Dyz_C3_bb;
  abcd[763] = 4.0E0*I_NAI_Fx2y_Dyz_C3_bb;
  abcd[764] = 4.0E0*I_NAI_Fxyz_Dyz_C3_bb;
  abcd[765] = 4.0E0*I_NAI_Fx2z_Dyz_C3_bb;
  abcd[766] = 4.0E0*I_NAI_F3y_Dyz_C3_bb;
  abcd[767] = 4.0E0*I_NAI_F2yz_Dyz_C3_bb;
  abcd[768] = 4.0E0*I_NAI_Fy2z_Dyz_C3_bb;
  abcd[769] = 4.0E0*I_NAI_F3z_Dyz_C3_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[770] = 4.0E0*I_NAI_F3x_Fxyz_C1003_bb;
  abcd[771] = 4.0E0*I_NAI_F2xy_Fxyz_C1003_bb;
  abcd[772] = 4.0E0*I_NAI_F2xz_Fxyz_C1003_bb;
  abcd[773] = 4.0E0*I_NAI_Fx2y_Fxyz_C1003_bb;
  abcd[774] = 4.0E0*I_NAI_Fxyz_Fxyz_C1003_bb;
  abcd[775] = 4.0E0*I_NAI_Fx2z_Fxyz_C1003_bb;
  abcd[776] = 4.0E0*I_NAI_F3y_Fxyz_C1003_bb;
  abcd[777] = 4.0E0*I_NAI_F2yz_Fxyz_C1003_bb;
  abcd[778] = 4.0E0*I_NAI_Fy2z_Fxyz_C1003_bb;
  abcd[779] = 4.0E0*I_NAI_F3z_Fxyz_C1003_bb;
  abcd[780] = 4.0E0*I_NAI_F3x_F2yz_C1003_bb-2.0E0*1*I_NAI_F3x_Pz_C1003_b;
  abcd[781] = 4.0E0*I_NAI_F2xy_F2yz_C1003_bb-2.0E0*1*I_NAI_F2xy_Pz_C1003_b;
  abcd[782] = 4.0E0*I_NAI_F2xz_F2yz_C1003_bb-2.0E0*1*I_NAI_F2xz_Pz_C1003_b;
  abcd[783] = 4.0E0*I_NAI_Fx2y_F2yz_C1003_bb-2.0E0*1*I_NAI_Fx2y_Pz_C1003_b;
  abcd[784] = 4.0E0*I_NAI_Fxyz_F2yz_C1003_bb-2.0E0*1*I_NAI_Fxyz_Pz_C1003_b;
  abcd[785] = 4.0E0*I_NAI_Fx2z_F2yz_C1003_bb-2.0E0*1*I_NAI_Fx2z_Pz_C1003_b;
  abcd[786] = 4.0E0*I_NAI_F3y_F2yz_C1003_bb-2.0E0*1*I_NAI_F3y_Pz_C1003_b;
  abcd[787] = 4.0E0*I_NAI_F2yz_F2yz_C1003_bb-2.0E0*1*I_NAI_F2yz_Pz_C1003_b;
  abcd[788] = 4.0E0*I_NAI_Fy2z_F2yz_C1003_bb-2.0E0*1*I_NAI_Fy2z_Pz_C1003_b;
  abcd[789] = 4.0E0*I_NAI_F3z_F2yz_C1003_bb-2.0E0*1*I_NAI_F3z_Pz_C1003_b;
  abcd[790] = 4.0E0*I_NAI_F3x_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3x_Py_C1003_b;
  abcd[791] = 4.0E0*I_NAI_F2xy_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2xy_Py_C1003_b;
  abcd[792] = 4.0E0*I_NAI_F2xz_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2xz_Py_C1003_b;
  abcd[793] = 4.0E0*I_NAI_Fx2y_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fx2y_Py_C1003_b;
  abcd[794] = 4.0E0*I_NAI_Fxyz_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fxyz_Py_C1003_b;
  abcd[795] = 4.0E0*I_NAI_Fx2z_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fx2z_Py_C1003_b;
  abcd[796] = 4.0E0*I_NAI_F3y_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3y_Py_C1003_b;
  abcd[797] = 4.0E0*I_NAI_F2yz_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2yz_Py_C1003_b;
  abcd[798] = 4.0E0*I_NAI_Fy2z_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fy2z_Py_C1003_b;
  abcd[799] = 4.0E0*I_NAI_F3z_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3z_Py_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C3_bb
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  abcd[800] = 4.0E0*I_NAI_F3x_D2z_C3_bb-2.0E0*1*I_NAI_F3x_S_C3_b;
  abcd[801] = 4.0E0*I_NAI_F2xy_D2z_C3_bb-2.0E0*1*I_NAI_F2xy_S_C3_b;
  abcd[802] = 4.0E0*I_NAI_F2xz_D2z_C3_bb-2.0E0*1*I_NAI_F2xz_S_C3_b;
  abcd[803] = 4.0E0*I_NAI_Fx2y_D2z_C3_bb-2.0E0*1*I_NAI_Fx2y_S_C3_b;
  abcd[804] = 4.0E0*I_NAI_Fxyz_D2z_C3_bb-2.0E0*1*I_NAI_Fxyz_S_C3_b;
  abcd[805] = 4.0E0*I_NAI_Fx2z_D2z_C3_bb-2.0E0*1*I_NAI_Fx2z_S_C3_b;
  abcd[806] = 4.0E0*I_NAI_F3y_D2z_C3_bb-2.0E0*1*I_NAI_F3y_S_C3_b;
  abcd[807] = 4.0E0*I_NAI_F2yz_D2z_C3_bb-2.0E0*1*I_NAI_F2yz_S_C3_b;
  abcd[808] = 4.0E0*I_NAI_Fy2z_D2z_C3_bb-2.0E0*1*I_NAI_Fy2z_S_C3_b;
  abcd[809] = 4.0E0*I_NAI_F3z_D2z_C3_bb-2.0E0*1*I_NAI_F3z_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_C1003_bb
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  abcd[810] = 4.0E0*I_NAI_F3x_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3x_Px_C1003_b;
  abcd[811] = 4.0E0*I_NAI_F2xy_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2xy_Px_C1003_b;
  abcd[812] = 4.0E0*I_NAI_F2xz_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2xz_Px_C1003_b;
  abcd[813] = 4.0E0*I_NAI_Fx2y_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fx2y_Px_C1003_b;
  abcd[814] = 4.0E0*I_NAI_Fxyz_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fxyz_Px_C1003_b;
  abcd[815] = 4.0E0*I_NAI_Fx2z_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fx2z_Px_C1003_b;
  abcd[816] = 4.0E0*I_NAI_F3y_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3y_Px_C1003_b;
  abcd[817] = 4.0E0*I_NAI_F2yz_Fx2z_C1003_bb-2.0E0*1*I_NAI_F2yz_Px_C1003_b;
  abcd[818] = 4.0E0*I_NAI_Fy2z_Fx2z_C1003_bb-2.0E0*1*I_NAI_Fy2z_Px_C1003_b;
  abcd[819] = 4.0E0*I_NAI_F3z_Fx2z_C1003_bb-2.0E0*1*I_NAI_F3z_Px_C1003_b;
  abcd[820] = 4.0E0*I_NAI_F3x_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3x_Py_C1003_b;
  abcd[821] = 4.0E0*I_NAI_F2xy_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2xy_Py_C1003_b;
  abcd[822] = 4.0E0*I_NAI_F2xz_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2xz_Py_C1003_b;
  abcd[823] = 4.0E0*I_NAI_Fx2y_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fx2y_Py_C1003_b;
  abcd[824] = 4.0E0*I_NAI_Fxyz_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fxyz_Py_C1003_b;
  abcd[825] = 4.0E0*I_NAI_Fx2z_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fx2z_Py_C1003_b;
  abcd[826] = 4.0E0*I_NAI_F3y_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3y_Py_C1003_b;
  abcd[827] = 4.0E0*I_NAI_F2yz_Fy2z_C1003_bb-2.0E0*1*I_NAI_F2yz_Py_C1003_b;
  abcd[828] = 4.0E0*I_NAI_Fy2z_Fy2z_C1003_bb-2.0E0*1*I_NAI_Fy2z_Py_C1003_b;
  abcd[829] = 4.0E0*I_NAI_F3z_Fy2z_C1003_bb-2.0E0*1*I_NAI_F3z_Py_C1003_b;
  abcd[830] = 4.0E0*I_NAI_F3x_F3z_C1003_bb-2.0E0*1*I_NAI_F3x_Pz_C1003_b-2.0E0*2*I_NAI_F3x_Pz_C1003_b;
  abcd[831] = 4.0E0*I_NAI_F2xy_F3z_C1003_bb-2.0E0*1*I_NAI_F2xy_Pz_C1003_b-2.0E0*2*I_NAI_F2xy_Pz_C1003_b;
  abcd[832] = 4.0E0*I_NAI_F2xz_F3z_C1003_bb-2.0E0*1*I_NAI_F2xz_Pz_C1003_b-2.0E0*2*I_NAI_F2xz_Pz_C1003_b;
  abcd[833] = 4.0E0*I_NAI_Fx2y_F3z_C1003_bb-2.0E0*1*I_NAI_Fx2y_Pz_C1003_b-2.0E0*2*I_NAI_Fx2y_Pz_C1003_b;
  abcd[834] = 4.0E0*I_NAI_Fxyz_F3z_C1003_bb-2.0E0*1*I_NAI_Fxyz_Pz_C1003_b-2.0E0*2*I_NAI_Fxyz_Pz_C1003_b;
  abcd[835] = 4.0E0*I_NAI_Fx2z_F3z_C1003_bb-2.0E0*1*I_NAI_Fx2z_Pz_C1003_b-2.0E0*2*I_NAI_Fx2z_Pz_C1003_b;
  abcd[836] = 4.0E0*I_NAI_F3y_F3z_C1003_bb-2.0E0*1*I_NAI_F3y_Pz_C1003_b-2.0E0*2*I_NAI_F3y_Pz_C1003_b;
  abcd[837] = 4.0E0*I_NAI_F2yz_F3z_C1003_bb-2.0E0*1*I_NAI_F2yz_Pz_C1003_b-2.0E0*2*I_NAI_F2yz_Pz_C1003_b;
  abcd[838] = 4.0E0*I_NAI_Fy2z_F3z_C1003_bb-2.0E0*1*I_NAI_Fy2z_Pz_C1003_b-2.0E0*2*I_NAI_Fy2z_Pz_C1003_b;
  abcd[839] = 4.0E0*I_NAI_F3z_F3z_C1003_bb-2.0E0*1*I_NAI_F3z_Pz_C1003_b-2.0E0*2*I_NAI_F3z_Pz_C1003_b;
}
