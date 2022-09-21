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

void hgp_os_nai_d_d_d2(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_S_D2x = 0.0E0;
  Double I_NAI_S_Dxy = 0.0E0;
  Double I_NAI_S_Dxz = 0.0E0;
  Double I_NAI_S_D2y = 0.0E0;
  Double I_NAI_S_Dyz = 0.0E0;
  Double I_NAI_S_D2z = 0.0E0;
  Double I_NAI_D2x_S = 0.0E0;
  Double I_NAI_Dxy_S = 0.0E0;
  Double I_NAI_Dxz_S = 0.0E0;
  Double I_NAI_D2y_S = 0.0E0;
  Double I_NAI_Dyz_S = 0.0E0;
  Double I_NAI_D2z_S = 0.0E0;
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
  Double I_NAI_F3x_S_a = 0.0E0;
  Double I_NAI_F2xy_S_a = 0.0E0;
  Double I_NAI_F2xz_S_a = 0.0E0;
  Double I_NAI_Fx2y_S_a = 0.0E0;
  Double I_NAI_Fxyz_S_a = 0.0E0;
  Double I_NAI_Fx2z_S_a = 0.0E0;
  Double I_NAI_F3y_S_a = 0.0E0;
  Double I_NAI_F2yz_S_a = 0.0E0;
  Double I_NAI_Fy2z_S_a = 0.0E0;
  Double I_NAI_F3z_S_a = 0.0E0;
  Double I_NAI_Px_S = 0.0E0;
  Double I_NAI_Py_S = 0.0E0;
  Double I_NAI_Pz_S = 0.0E0;
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
  Double I_NAI_H5x_S_aa = 0.0E0;
  Double I_NAI_H4xy_S_aa = 0.0E0;
  Double I_NAI_H4xz_S_aa = 0.0E0;
  Double I_NAI_H3x2y_S_aa = 0.0E0;
  Double I_NAI_H3xyz_S_aa = 0.0E0;
  Double I_NAI_H3x2z_S_aa = 0.0E0;
  Double I_NAI_H2x3y_S_aa = 0.0E0;
  Double I_NAI_H2x2yz_S_aa = 0.0E0;
  Double I_NAI_H2xy2z_S_aa = 0.0E0;
  Double I_NAI_H2x3z_S_aa = 0.0E0;
  Double I_NAI_Hx4y_S_aa = 0.0E0;
  Double I_NAI_Hx3yz_S_aa = 0.0E0;
  Double I_NAI_Hx2y2z_S_aa = 0.0E0;
  Double I_NAI_Hxy3z_S_aa = 0.0E0;
  Double I_NAI_Hx4z_S_aa = 0.0E0;
  Double I_NAI_H5y_S_aa = 0.0E0;
  Double I_NAI_H4yz_S_aa = 0.0E0;
  Double I_NAI_H3y2z_S_aa = 0.0E0;
  Double I_NAI_H2y3z_S_aa = 0.0E0;
  Double I_NAI_Hy4z_S_aa = 0.0E0;
  Double I_NAI_H5z_S_aa = 0.0E0;
  Double I_NAI_G4x_S_aa = 0.0E0;
  Double I_NAI_G3xy_S_aa = 0.0E0;
  Double I_NAI_G3xz_S_aa = 0.0E0;
  Double I_NAI_G2x2y_S_aa = 0.0E0;
  Double I_NAI_G2xyz_S_aa = 0.0E0;
  Double I_NAI_G2x2z_S_aa = 0.0E0;
  Double I_NAI_Gx3y_S_aa = 0.0E0;
  Double I_NAI_Gx2yz_S_aa = 0.0E0;
  Double I_NAI_Gxy2z_S_aa = 0.0E0;
  Double I_NAI_Gx3z_S_aa = 0.0E0;
  Double I_NAI_G4y_S_aa = 0.0E0;
  Double I_NAI_G3yz_S_aa = 0.0E0;
  Double I_NAI_G2y2z_S_aa = 0.0E0;
  Double I_NAI_Gy3z_S_aa = 0.0E0;
  Double I_NAI_G4z_S_aa = 0.0E0;
  Double I_NAI_D2x_S_a = 0.0E0;
  Double I_NAI_Dxy_S_a = 0.0E0;
  Double I_NAI_Dxz_S_a = 0.0E0;
  Double I_NAI_D2y_S_a = 0.0E0;
  Double I_NAI_Dyz_S_a = 0.0E0;
  Double I_NAI_D2z_S_a = 0.0E0;
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
  Double I_NAI_D2x_S_b = 0.0E0;
  Double I_NAI_Dxy_S_b = 0.0E0;
  Double I_NAI_Dxz_S_b = 0.0E0;
  Double I_NAI_D2y_S_b = 0.0E0;
  Double I_NAI_Dyz_S_b = 0.0E0;
  Double I_NAI_D2z_S_b = 0.0E0;
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
  Double I_NAI_G4x_S_ab = 0.0E0;
  Double I_NAI_G3xy_S_ab = 0.0E0;
  Double I_NAI_G3xz_S_ab = 0.0E0;
  Double I_NAI_G2x2y_S_ab = 0.0E0;
  Double I_NAI_G2xyz_S_ab = 0.0E0;
  Double I_NAI_G2x2z_S_ab = 0.0E0;
  Double I_NAI_Gx3y_S_ab = 0.0E0;
  Double I_NAI_Gx2yz_S_ab = 0.0E0;
  Double I_NAI_Gxy2z_S_ab = 0.0E0;
  Double I_NAI_Gx3z_S_ab = 0.0E0;
  Double I_NAI_G4y_S_ab = 0.0E0;
  Double I_NAI_G3yz_S_ab = 0.0E0;
  Double I_NAI_G2y2z_S_ab = 0.0E0;
  Double I_NAI_Gy3z_S_ab = 0.0E0;
  Double I_NAI_G4z_S_ab = 0.0E0;
  Double I_NAI_F3x_S_ab = 0.0E0;
  Double I_NAI_F2xy_S_ab = 0.0E0;
  Double I_NAI_F2xz_S_ab = 0.0E0;
  Double I_NAI_Fx2y_S_ab = 0.0E0;
  Double I_NAI_Fxyz_S_ab = 0.0E0;
  Double I_NAI_Fx2z_S_ab = 0.0E0;
  Double I_NAI_F3y_S_ab = 0.0E0;
  Double I_NAI_F2yz_S_ab = 0.0E0;
  Double I_NAI_Fy2z_S_ab = 0.0E0;
  Double I_NAI_F3z_S_ab = 0.0E0;
  Double I_NAI_Px_S_b = 0.0E0;
  Double I_NAI_Py_S_b = 0.0E0;
  Double I_NAI_Pz_S_b = 0.0E0;
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
  Double I_NAI_F3x_S_bb = 0.0E0;
  Double I_NAI_F2xy_S_bb = 0.0E0;
  Double I_NAI_F2xz_S_bb = 0.0E0;
  Double I_NAI_Fx2y_S_bb = 0.0E0;
  Double I_NAI_Fxyz_S_bb = 0.0E0;
  Double I_NAI_Fx2z_S_bb = 0.0E0;
  Double I_NAI_F3y_S_bb = 0.0E0;
  Double I_NAI_F2yz_S_bb = 0.0E0;
  Double I_NAI_Fy2z_S_bb = 0.0E0;
  Double I_NAI_F3z_S_bb = 0.0E0;
  Double I_NAI_D2x_S_bb = 0.0E0;
  Double I_NAI_Dxy_S_bb = 0.0E0;
  Double I_NAI_Dxz_S_bb = 0.0E0;
  Double I_NAI_D2y_S_bb = 0.0E0;
  Double I_NAI_Dyz_S_bb = 0.0E0;
  Double I_NAI_D2z_S_bb = 0.0E0;

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
    Double PBX   = PX - B[0];
    Double PBY   = PY - B[1];
    Double PBZ   = PZ - B[2];
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
       * shell quartet name: SQ_NAI_S_P_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_S_Px_M1_vrr = PBX*I_NAI_S_S_M1_vrr-PNX*I_NAI_S_S_M2_vrr;
      Double I_NAI_S_Py_M1_vrr = PBY*I_NAI_S_S_M1_vrr-PNY*I_NAI_S_S_M2_vrr;
      Double I_NAI_S_Pz_M1_vrr = PBZ*I_NAI_S_S_M1_vrr-PNZ*I_NAI_S_S_M2_vrr;

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
       * shell quartet name: SQ_NAI_S_P
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_S_Px_vrr = PBX*I_NAI_S_S_vrr-PNX*I_NAI_S_S_M1_vrr;
      Double I_NAI_S_Py_vrr = PBY*I_NAI_S_S_vrr-PNY*I_NAI_S_S_M1_vrr;
      Double I_NAI_S_Pz_vrr = PBZ*I_NAI_S_S_vrr-PNZ*I_NAI_S_S_M1_vrr;

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
       * shell quartet name: SQ_NAI_S_D
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_P
       * RHS shell quartet name: SQ_NAI_S_P_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_S_D2x_vrr = PBX*I_NAI_S_Px_vrr-PNX*I_NAI_S_Px_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_S_Dxy_vrr = PBY*I_NAI_S_Px_vrr-PNY*I_NAI_S_Px_M1_vrr;
      Double I_NAI_S_Dxz_vrr = PBZ*I_NAI_S_Px_vrr-PNZ*I_NAI_S_Px_M1_vrr;
      Double I_NAI_S_D2y_vrr = PBY*I_NAI_S_Py_vrr-PNY*I_NAI_S_Py_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_S_Dyz_vrr = PBZ*I_NAI_S_Py_vrr-PNZ*I_NAI_S_Py_M1_vrr;
      Double I_NAI_S_D2z_vrr = PBZ*I_NAI_S_Pz_vrr-PNZ*I_NAI_S_Pz_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

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
       * shell quartet name: SQ_NAI_S_D
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_S_D2x += I_NAI_S_D2x_vrr;
      I_NAI_S_Dxy += I_NAI_S_Dxy_vrr;
      I_NAI_S_Dxz += I_NAI_S_Dxz_vrr;
      I_NAI_S_D2y += I_NAI_S_D2y_vrr;
      I_NAI_S_Dyz += I_NAI_S_Dyz_vrr;
      I_NAI_S_D2z += I_NAI_S_D2z_vrr;

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
       * shell quartet name: SQ_NAI_F_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_a_coefs = alpha;
      I_NAI_F3x_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_a += SQ_NAI_F_S_a_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_a += SQ_NAI_F_S_a_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_a += SQ_NAI_F_S_a_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_a += SQ_NAI_F_S_a_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_a += SQ_NAI_F_S_a_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_Px_S += I_NAI_Px_S_vrr;
      I_NAI_Py_S += I_NAI_Py_S_vrr;
      I_NAI_Pz_S += I_NAI_Pz_S_vrr;

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
       * shell quartet name: SQ_NAI_H_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_aa_coefs = alpha*alpha;
      I_NAI_H5x_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_aa += SQ_NAI_H_S_aa_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_aa_coefs = alpha*alpha;
      I_NAI_G4x_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_aa += SQ_NAI_G_S_aa_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_a_coefs = alpha;
      I_NAI_D2x_S_a += SQ_NAI_D_S_a_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_a += SQ_NAI_D_S_a_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_a += SQ_NAI_D_S_a_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_a += SQ_NAI_D_S_a_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_a += SQ_NAI_D_S_a_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_a += SQ_NAI_D_S_a_coefs*I_NAI_D2z_S_vrr;

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
       * shell quartet name: SQ_NAI_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_b_coefs = beta;
      I_NAI_D2x_S_b += SQ_NAI_D_S_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_b += SQ_NAI_D_S_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_b += SQ_NAI_D_S_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_b += SQ_NAI_D_S_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_b += SQ_NAI_D_S_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_b += SQ_NAI_D_S_b_coefs*I_NAI_D2z_S_vrr;

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
       * shell quartet name: SQ_NAI_G_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_ab_coefs = alpha*beta;
      I_NAI_G4x_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_ab += SQ_NAI_G_S_ab_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_ab_coefs = alpha*beta;
      I_NAI_F3x_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_ab += SQ_NAI_F_S_ab_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_b_coefs = beta;
      I_NAI_Px_S_b += SQ_NAI_P_S_b_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_b += SQ_NAI_P_S_b_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_b += SQ_NAI_P_S_b_coefs*I_NAI_Pz_S_vrr;

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

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_bb_coefs = beta*beta;
      I_NAI_F3x_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_bb += SQ_NAI_F_S_bb_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_bb_coefs = beta*beta;
      I_NAI_D2x_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_bb += SQ_NAI_D_S_bb_coefs*I_NAI_D2z_S_vrr;
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
   * shell quartet name: SQ_NAI_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S
   * RHS shell quartet name: SQ_NAI_P_S
   ************************************************************/
  Double I_NAI_Px_Px = I_NAI_D2x_S+ABX*I_NAI_Px_S;
  Double I_NAI_Py_Px = I_NAI_Dxy_S+ABX*I_NAI_Py_S;
  Double I_NAI_Pz_Px = I_NAI_Dxz_S+ABX*I_NAI_Pz_S;
  Double I_NAI_Px_Py = I_NAI_Dxy_S+ABY*I_NAI_Px_S;
  Double I_NAI_Py_Py = I_NAI_D2y_S+ABY*I_NAI_Py_S;
  Double I_NAI_Pz_Py = I_NAI_Dyz_S+ABY*I_NAI_Pz_S;
  Double I_NAI_Px_Pz = I_NAI_Dxz_S+ABZ*I_NAI_Px_S;
  Double I_NAI_Py_Pz = I_NAI_Dyz_S+ABZ*I_NAI_Py_S;
  Double I_NAI_Pz_Pz = I_NAI_D2z_S+ABZ*I_NAI_Pz_S;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_a
   * RHS shell quartet name: SQ_NAI_D_S_a
   ************************************************************/
  Double I_NAI_D2x_Px_a = I_NAI_F3x_S_a+ABX*I_NAI_D2x_S_a;
  Double I_NAI_Dxy_Px_a = I_NAI_F2xy_S_a+ABX*I_NAI_Dxy_S_a;
  Double I_NAI_Dxz_Px_a = I_NAI_F2xz_S_a+ABX*I_NAI_Dxz_S_a;
  Double I_NAI_D2y_Px_a = I_NAI_Fx2y_S_a+ABX*I_NAI_D2y_S_a;
  Double I_NAI_Dyz_Px_a = I_NAI_Fxyz_S_a+ABX*I_NAI_Dyz_S_a;
  Double I_NAI_D2z_Px_a = I_NAI_Fx2z_S_a+ABX*I_NAI_D2z_S_a;
  Double I_NAI_D2x_Py_a = I_NAI_F2xy_S_a+ABY*I_NAI_D2x_S_a;
  Double I_NAI_Dxy_Py_a = I_NAI_Fx2y_S_a+ABY*I_NAI_Dxy_S_a;
  Double I_NAI_Dxz_Py_a = I_NAI_Fxyz_S_a+ABY*I_NAI_Dxz_S_a;
  Double I_NAI_D2y_Py_a = I_NAI_F3y_S_a+ABY*I_NAI_D2y_S_a;
  Double I_NAI_Dyz_Py_a = I_NAI_F2yz_S_a+ABY*I_NAI_Dyz_S_a;
  Double I_NAI_D2z_Py_a = I_NAI_Fy2z_S_a+ABY*I_NAI_D2z_S_a;
  Double I_NAI_D2x_Pz_a = I_NAI_F2xz_S_a+ABZ*I_NAI_D2x_S_a;
  Double I_NAI_Dxy_Pz_a = I_NAI_Fxyz_S_a+ABZ*I_NAI_Dxy_S_a;
  Double I_NAI_Dxz_Pz_a = I_NAI_Fx2z_S_a+ABZ*I_NAI_Dxz_S_a;
  Double I_NAI_D2y_Pz_a = I_NAI_F2yz_S_a+ABZ*I_NAI_D2y_S_a;
  Double I_NAI_Dyz_Pz_a = I_NAI_Fy2z_S_a+ABZ*I_NAI_Dyz_S_a;
  Double I_NAI_D2z_Pz_a = I_NAI_F3z_S_a+ABZ*I_NAI_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_a
   * RHS shell quartet name: SQ_NAI_F_S_a
   ************************************************************/
  Double I_NAI_F3x_Px_a = I_NAI_G4x_S_a+ABX*I_NAI_F3x_S_a;
  Double I_NAI_F2xy_Px_a = I_NAI_G3xy_S_a+ABX*I_NAI_F2xy_S_a;
  Double I_NAI_F2xz_Px_a = I_NAI_G3xz_S_a+ABX*I_NAI_F2xz_S_a;
  Double I_NAI_Fx2y_Px_a = I_NAI_G2x2y_S_a+ABX*I_NAI_Fx2y_S_a;
  Double I_NAI_Fxyz_Px_a = I_NAI_G2xyz_S_a+ABX*I_NAI_Fxyz_S_a;
  Double I_NAI_Fx2z_Px_a = I_NAI_G2x2z_S_a+ABX*I_NAI_Fx2z_S_a;
  Double I_NAI_F3y_Px_a = I_NAI_Gx3y_S_a+ABX*I_NAI_F3y_S_a;
  Double I_NAI_F2yz_Px_a = I_NAI_Gx2yz_S_a+ABX*I_NAI_F2yz_S_a;
  Double I_NAI_Fy2z_Px_a = I_NAI_Gxy2z_S_a+ABX*I_NAI_Fy2z_S_a;
  Double I_NAI_F3z_Px_a = I_NAI_Gx3z_S_a+ABX*I_NAI_F3z_S_a;
  Double I_NAI_F3x_Py_a = I_NAI_G3xy_S_a+ABY*I_NAI_F3x_S_a;
  Double I_NAI_F2xy_Py_a = I_NAI_G2x2y_S_a+ABY*I_NAI_F2xy_S_a;
  Double I_NAI_F2xz_Py_a = I_NAI_G2xyz_S_a+ABY*I_NAI_F2xz_S_a;
  Double I_NAI_Fx2y_Py_a = I_NAI_Gx3y_S_a+ABY*I_NAI_Fx2y_S_a;
  Double I_NAI_Fxyz_Py_a = I_NAI_Gx2yz_S_a+ABY*I_NAI_Fxyz_S_a;
  Double I_NAI_Fx2z_Py_a = I_NAI_Gxy2z_S_a+ABY*I_NAI_Fx2z_S_a;
  Double I_NAI_F3y_Py_a = I_NAI_G4y_S_a+ABY*I_NAI_F3y_S_a;
  Double I_NAI_F2yz_Py_a = I_NAI_G3yz_S_a+ABY*I_NAI_F2yz_S_a;
  Double I_NAI_Fy2z_Py_a = I_NAI_G2y2z_S_a+ABY*I_NAI_Fy2z_S_a;
  Double I_NAI_F3z_Py_a = I_NAI_Gy3z_S_a+ABY*I_NAI_F3z_S_a;
  Double I_NAI_F3x_Pz_a = I_NAI_G3xz_S_a+ABZ*I_NAI_F3x_S_a;
  Double I_NAI_F2xy_Pz_a = I_NAI_G2xyz_S_a+ABZ*I_NAI_F2xy_S_a;
  Double I_NAI_F2xz_Pz_a = I_NAI_G2x2z_S_a+ABZ*I_NAI_F2xz_S_a;
  Double I_NAI_Fx2y_Pz_a = I_NAI_Gx2yz_S_a+ABZ*I_NAI_Fx2y_S_a;
  Double I_NAI_Fxyz_Pz_a = I_NAI_Gxy2z_S_a+ABZ*I_NAI_Fxyz_S_a;
  Double I_NAI_Fx2z_Pz_a = I_NAI_Gx3z_S_a+ABZ*I_NAI_Fx2z_S_a;
  Double I_NAI_F3y_Pz_a = I_NAI_G3yz_S_a+ABZ*I_NAI_F3y_S_a;
  Double I_NAI_F2yz_Pz_a = I_NAI_G2y2z_S_a+ABZ*I_NAI_F2yz_S_a;
  Double I_NAI_Fy2z_Pz_a = I_NAI_Gy3z_S_a+ABZ*I_NAI_Fy2z_S_a;
  Double I_NAI_F3z_Pz_a = I_NAI_G4z_S_a+ABZ*I_NAI_F3z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_D_P_a
   ************************************************************/
  Double I_NAI_D2x_D2x_a = I_NAI_F3x_Px_a+ABX*I_NAI_D2x_Px_a;
  Double I_NAI_Dxy_D2x_a = I_NAI_F2xy_Px_a+ABX*I_NAI_Dxy_Px_a;
  Double I_NAI_Dxz_D2x_a = I_NAI_F2xz_Px_a+ABX*I_NAI_Dxz_Px_a;
  Double I_NAI_D2y_D2x_a = I_NAI_Fx2y_Px_a+ABX*I_NAI_D2y_Px_a;
  Double I_NAI_Dyz_D2x_a = I_NAI_Fxyz_Px_a+ABX*I_NAI_Dyz_Px_a;
  Double I_NAI_D2z_D2x_a = I_NAI_Fx2z_Px_a+ABX*I_NAI_D2z_Px_a;
  Double I_NAI_D2x_Dxy_a = I_NAI_F2xy_Px_a+ABY*I_NAI_D2x_Px_a;
  Double I_NAI_Dxy_Dxy_a = I_NAI_Fx2y_Px_a+ABY*I_NAI_Dxy_Px_a;
  Double I_NAI_Dxz_Dxy_a = I_NAI_Fxyz_Px_a+ABY*I_NAI_Dxz_Px_a;
  Double I_NAI_D2y_Dxy_a = I_NAI_F3y_Px_a+ABY*I_NAI_D2y_Px_a;
  Double I_NAI_Dyz_Dxy_a = I_NAI_F2yz_Px_a+ABY*I_NAI_Dyz_Px_a;
  Double I_NAI_D2z_Dxy_a = I_NAI_Fy2z_Px_a+ABY*I_NAI_D2z_Px_a;
  Double I_NAI_D2x_Dxz_a = I_NAI_F2xz_Px_a+ABZ*I_NAI_D2x_Px_a;
  Double I_NAI_Dxy_Dxz_a = I_NAI_Fxyz_Px_a+ABZ*I_NAI_Dxy_Px_a;
  Double I_NAI_Dxz_Dxz_a = I_NAI_Fx2z_Px_a+ABZ*I_NAI_Dxz_Px_a;
  Double I_NAI_D2y_Dxz_a = I_NAI_F2yz_Px_a+ABZ*I_NAI_D2y_Px_a;
  Double I_NAI_Dyz_Dxz_a = I_NAI_Fy2z_Px_a+ABZ*I_NAI_Dyz_Px_a;
  Double I_NAI_D2z_Dxz_a = I_NAI_F3z_Px_a+ABZ*I_NAI_D2z_Px_a;
  Double I_NAI_D2x_D2y_a = I_NAI_F2xy_Py_a+ABY*I_NAI_D2x_Py_a;
  Double I_NAI_Dxy_D2y_a = I_NAI_Fx2y_Py_a+ABY*I_NAI_Dxy_Py_a;
  Double I_NAI_Dxz_D2y_a = I_NAI_Fxyz_Py_a+ABY*I_NAI_Dxz_Py_a;
  Double I_NAI_D2y_D2y_a = I_NAI_F3y_Py_a+ABY*I_NAI_D2y_Py_a;
  Double I_NAI_Dyz_D2y_a = I_NAI_F2yz_Py_a+ABY*I_NAI_Dyz_Py_a;
  Double I_NAI_D2z_D2y_a = I_NAI_Fy2z_Py_a+ABY*I_NAI_D2z_Py_a;
  Double I_NAI_D2x_Dyz_a = I_NAI_F2xz_Py_a+ABZ*I_NAI_D2x_Py_a;
  Double I_NAI_Dxy_Dyz_a = I_NAI_Fxyz_Py_a+ABZ*I_NAI_Dxy_Py_a;
  Double I_NAI_Dxz_Dyz_a = I_NAI_Fx2z_Py_a+ABZ*I_NAI_Dxz_Py_a;
  Double I_NAI_D2y_Dyz_a = I_NAI_F2yz_Py_a+ABZ*I_NAI_D2y_Py_a;
  Double I_NAI_Dyz_Dyz_a = I_NAI_Fy2z_Py_a+ABZ*I_NAI_Dyz_Py_a;
  Double I_NAI_D2z_Dyz_a = I_NAI_F3z_Py_a+ABZ*I_NAI_D2z_Py_a;
  Double I_NAI_D2x_D2z_a = I_NAI_F2xz_Pz_a+ABZ*I_NAI_D2x_Pz_a;
  Double I_NAI_Dxy_D2z_a = I_NAI_Fxyz_Pz_a+ABZ*I_NAI_Dxy_Pz_a;
  Double I_NAI_Dxz_D2z_a = I_NAI_Fx2z_Pz_a+ABZ*I_NAI_Dxz_Pz_a;
  Double I_NAI_D2y_D2z_a = I_NAI_F2yz_Pz_a+ABZ*I_NAI_D2y_Pz_a;
  Double I_NAI_Dyz_D2z_a = I_NAI_Fy2z_Pz_a+ABZ*I_NAI_Dyz_Pz_a;
  Double I_NAI_D2z_D2z_a = I_NAI_F3z_Pz_a+ABZ*I_NAI_D2z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_b
   * RHS shell quartet name: SQ_NAI_P_S_b
   ************************************************************/
  Double I_NAI_Px_Px_b = I_NAI_D2x_S_b+ABX*I_NAI_Px_S_b;
  Double I_NAI_Py_Px_b = I_NAI_Dxy_S_b+ABX*I_NAI_Py_S_b;
  Double I_NAI_Pz_Px_b = I_NAI_Dxz_S_b+ABX*I_NAI_Pz_S_b;
  Double I_NAI_Px_Py_b = I_NAI_Dxy_S_b+ABY*I_NAI_Px_S_b;
  Double I_NAI_Py_Py_b = I_NAI_D2y_S_b+ABY*I_NAI_Py_S_b;
  Double I_NAI_Pz_Py_b = I_NAI_Dyz_S_b+ABY*I_NAI_Pz_S_b;
  Double I_NAI_Px_Pz_b = I_NAI_Dxz_S_b+ABZ*I_NAI_Px_S_b;
  Double I_NAI_Py_Pz_b = I_NAI_Dyz_S_b+ABZ*I_NAI_Py_S_b;
  Double I_NAI_Pz_Pz_b = I_NAI_D2z_S_b+ABZ*I_NAI_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_b
   * RHS shell quartet name: SQ_NAI_D_S_b
   ************************************************************/
  Double I_NAI_D2x_Px_b = I_NAI_F3x_S_b+ABX*I_NAI_D2x_S_b;
  Double I_NAI_Dxy_Px_b = I_NAI_F2xy_S_b+ABX*I_NAI_Dxy_S_b;
  Double I_NAI_Dxz_Px_b = I_NAI_F2xz_S_b+ABX*I_NAI_Dxz_S_b;
  Double I_NAI_D2y_Px_b = I_NAI_Fx2y_S_b+ABX*I_NAI_D2y_S_b;
  Double I_NAI_Dyz_Px_b = I_NAI_Fxyz_S_b+ABX*I_NAI_Dyz_S_b;
  Double I_NAI_D2z_Px_b = I_NAI_Fx2z_S_b+ABX*I_NAI_D2z_S_b;
  Double I_NAI_D2x_Py_b = I_NAI_F2xy_S_b+ABY*I_NAI_D2x_S_b;
  Double I_NAI_Dxy_Py_b = I_NAI_Fx2y_S_b+ABY*I_NAI_Dxy_S_b;
  Double I_NAI_Dxz_Py_b = I_NAI_Fxyz_S_b+ABY*I_NAI_Dxz_S_b;
  Double I_NAI_D2y_Py_b = I_NAI_F3y_S_b+ABY*I_NAI_D2y_S_b;
  Double I_NAI_Dyz_Py_b = I_NAI_F2yz_S_b+ABY*I_NAI_Dyz_S_b;
  Double I_NAI_D2z_Py_b = I_NAI_Fy2z_S_b+ABY*I_NAI_D2z_S_b;
  Double I_NAI_D2x_Pz_b = I_NAI_F2xz_S_b+ABZ*I_NAI_D2x_S_b;
  Double I_NAI_Dxy_Pz_b = I_NAI_Fxyz_S_b+ABZ*I_NAI_Dxy_S_b;
  Double I_NAI_Dxz_Pz_b = I_NAI_Fx2z_S_b+ABZ*I_NAI_Dxz_S_b;
  Double I_NAI_D2y_Pz_b = I_NAI_F2yz_S_b+ABZ*I_NAI_D2y_S_b;
  Double I_NAI_Dyz_Pz_b = I_NAI_Fy2z_S_b+ABZ*I_NAI_Dyz_S_b;
  Double I_NAI_D2z_Pz_b = I_NAI_F3z_S_b+ABZ*I_NAI_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_b
   * RHS shell quartet name: SQ_NAI_P_P_b
   ************************************************************/
  Double I_NAI_Px_D2x_b = I_NAI_D2x_Px_b+ABX*I_NAI_Px_Px_b;
  Double I_NAI_Py_D2x_b = I_NAI_Dxy_Px_b+ABX*I_NAI_Py_Px_b;
  Double I_NAI_Pz_D2x_b = I_NAI_Dxz_Px_b+ABX*I_NAI_Pz_Px_b;
  Double I_NAI_Px_Dxy_b = I_NAI_Dxy_Px_b+ABY*I_NAI_Px_Px_b;
  Double I_NAI_Py_Dxy_b = I_NAI_D2y_Px_b+ABY*I_NAI_Py_Px_b;
  Double I_NAI_Pz_Dxy_b = I_NAI_Dyz_Px_b+ABY*I_NAI_Pz_Px_b;
  Double I_NAI_Px_D2y_b = I_NAI_Dxy_Py_b+ABY*I_NAI_Px_Py_b;
  Double I_NAI_Py_D2y_b = I_NAI_D2y_Py_b+ABY*I_NAI_Py_Py_b;
  Double I_NAI_Pz_D2y_b = I_NAI_Dyz_Py_b+ABY*I_NAI_Pz_Py_b;
  Double I_NAI_Px_D2z_b = I_NAI_Dxz_Pz_b+ABZ*I_NAI_Px_Pz_b;
  Double I_NAI_Py_D2z_b = I_NAI_Dyz_Pz_b+ABZ*I_NAI_Py_Pz_b;
  Double I_NAI_Pz_D2z_b = I_NAI_D2z_Pz_b+ABZ*I_NAI_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
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
  Double I_NAI_F2xy_Py_b = I_NAI_G2x2y_S_b+ABY*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Py_b = I_NAI_G2xyz_S_b+ABY*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Py_b = I_NAI_Gx3y_S_b+ABY*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Py_b = I_NAI_Gx2yz_S_b+ABY*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Py_b = I_NAI_Gxy2z_S_b+ABY*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Py_b = I_NAI_G4y_S_b+ABY*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Py_b = I_NAI_G3yz_S_b+ABY*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Py_b = I_NAI_G2y2z_S_b+ABY*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Py_b = I_NAI_Gy3z_S_b+ABY*I_NAI_F3z_S_b;
  Double I_NAI_F2xz_Pz_b = I_NAI_G2x2z_S_b+ABZ*I_NAI_F2xz_S_b;
  Double I_NAI_Fxyz_Pz_b = I_NAI_Gxy2z_S_b+ABZ*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Pz_b = I_NAI_Gx3z_S_b+ABZ*I_NAI_Fx2z_S_b;
  Double I_NAI_F2yz_Pz_b = I_NAI_G2y2z_S_b+ABZ*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Pz_b = I_NAI_Gy3z_S_b+ABZ*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Pz_b = I_NAI_G4z_S_b+ABZ*I_NAI_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_b
   * RHS shell quartet name: SQ_NAI_D_P_b
   ************************************************************/
  Double I_NAI_D2x_D2x_b = I_NAI_F3x_Px_b+ABX*I_NAI_D2x_Px_b;
  Double I_NAI_Dxy_D2x_b = I_NAI_F2xy_Px_b+ABX*I_NAI_Dxy_Px_b;
  Double I_NAI_Dxz_D2x_b = I_NAI_F2xz_Px_b+ABX*I_NAI_Dxz_Px_b;
  Double I_NAI_D2y_D2x_b = I_NAI_Fx2y_Px_b+ABX*I_NAI_D2y_Px_b;
  Double I_NAI_Dyz_D2x_b = I_NAI_Fxyz_Px_b+ABX*I_NAI_Dyz_Px_b;
  Double I_NAI_D2z_D2x_b = I_NAI_Fx2z_Px_b+ABX*I_NAI_D2z_Px_b;
  Double I_NAI_D2x_Dxy_b = I_NAI_F2xy_Px_b+ABY*I_NAI_D2x_Px_b;
  Double I_NAI_Dxy_Dxy_b = I_NAI_Fx2y_Px_b+ABY*I_NAI_Dxy_Px_b;
  Double I_NAI_Dxz_Dxy_b = I_NAI_Fxyz_Px_b+ABY*I_NAI_Dxz_Px_b;
  Double I_NAI_D2y_Dxy_b = I_NAI_F3y_Px_b+ABY*I_NAI_D2y_Px_b;
  Double I_NAI_Dyz_Dxy_b = I_NAI_F2yz_Px_b+ABY*I_NAI_Dyz_Px_b;
  Double I_NAI_D2z_Dxy_b = I_NAI_Fy2z_Px_b+ABY*I_NAI_D2z_Px_b;
  Double I_NAI_D2x_Dxz_b = I_NAI_F2xz_Px_b+ABZ*I_NAI_D2x_Px_b;
  Double I_NAI_Dxy_Dxz_b = I_NAI_Fxyz_Px_b+ABZ*I_NAI_Dxy_Px_b;
  Double I_NAI_Dxz_Dxz_b = I_NAI_Fx2z_Px_b+ABZ*I_NAI_Dxz_Px_b;
  Double I_NAI_D2y_Dxz_b = I_NAI_F2yz_Px_b+ABZ*I_NAI_D2y_Px_b;
  Double I_NAI_Dyz_Dxz_b = I_NAI_Fy2z_Px_b+ABZ*I_NAI_Dyz_Px_b;
  Double I_NAI_D2z_Dxz_b = I_NAI_F3z_Px_b+ABZ*I_NAI_D2z_Px_b;
  Double I_NAI_D2x_D2y_b = I_NAI_F2xy_Py_b+ABY*I_NAI_D2x_Py_b;
  Double I_NAI_Dxy_D2y_b = I_NAI_Fx2y_Py_b+ABY*I_NAI_Dxy_Py_b;
  Double I_NAI_Dxz_D2y_b = I_NAI_Fxyz_Py_b+ABY*I_NAI_Dxz_Py_b;
  Double I_NAI_D2y_D2y_b = I_NAI_F3y_Py_b+ABY*I_NAI_D2y_Py_b;
  Double I_NAI_Dyz_D2y_b = I_NAI_F2yz_Py_b+ABY*I_NAI_Dyz_Py_b;
  Double I_NAI_D2z_D2y_b = I_NAI_Fy2z_Py_b+ABY*I_NAI_D2z_Py_b;
  Double I_NAI_D2x_Dyz_b = I_NAI_F2xz_Py_b+ABZ*I_NAI_D2x_Py_b;
  Double I_NAI_Dxy_Dyz_b = I_NAI_Fxyz_Py_b+ABZ*I_NAI_Dxy_Py_b;
  Double I_NAI_Dxz_Dyz_b = I_NAI_Fx2z_Py_b+ABZ*I_NAI_Dxz_Py_b;
  Double I_NAI_D2y_Dyz_b = I_NAI_F2yz_Py_b+ABZ*I_NAI_D2y_Py_b;
  Double I_NAI_Dyz_Dyz_b = I_NAI_Fy2z_Py_b+ABZ*I_NAI_Dyz_Py_b;
  Double I_NAI_D2z_Dyz_b = I_NAI_F3z_Py_b+ABZ*I_NAI_D2z_Py_b;
  Double I_NAI_D2x_D2z_b = I_NAI_F2xz_Pz_b+ABZ*I_NAI_D2x_Pz_b;
  Double I_NAI_Dxy_D2z_b = I_NAI_Fxyz_Pz_b+ABZ*I_NAI_Dxy_Pz_b;
  Double I_NAI_Dxz_D2z_b = I_NAI_Fx2z_Pz_b+ABZ*I_NAI_Dxz_Pz_b;
  Double I_NAI_D2y_D2z_b = I_NAI_F2yz_Pz_b+ABZ*I_NAI_D2y_Pz_b;
  Double I_NAI_Dyz_D2z_b = I_NAI_Fy2z_Pz_b+ABZ*I_NAI_Dyz_Pz_b;
  Double I_NAI_D2z_D2z_b = I_NAI_F3z_Pz_b+ABZ*I_NAI_D2z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_P_D_b
   ************************************************************/
  Double I_NAI_Px_F3x_b = I_NAI_D2x_D2x_b+ABX*I_NAI_Px_D2x_b;
  Double I_NAI_Py_F3x_b = I_NAI_Dxy_D2x_b+ABX*I_NAI_Py_D2x_b;
  Double I_NAI_Pz_F3x_b = I_NAI_Dxz_D2x_b+ABX*I_NAI_Pz_D2x_b;
  Double I_NAI_Px_F2xy_b = I_NAI_Dxy_D2x_b+ABY*I_NAI_Px_D2x_b;
  Double I_NAI_Py_F2xy_b = I_NAI_D2y_D2x_b+ABY*I_NAI_Py_D2x_b;
  Double I_NAI_Pz_F2xy_b = I_NAI_Dyz_D2x_b+ABY*I_NAI_Pz_D2x_b;
  Double I_NAI_Px_F2xz_b = I_NAI_Dxz_D2x_b+ABZ*I_NAI_Px_D2x_b;
  Double I_NAI_Py_F2xz_b = I_NAI_Dyz_D2x_b+ABZ*I_NAI_Py_D2x_b;
  Double I_NAI_Pz_F2xz_b = I_NAI_D2z_D2x_b+ABZ*I_NAI_Pz_D2x_b;
  Double I_NAI_Px_Fx2y_b = I_NAI_D2x_D2y_b+ABX*I_NAI_Px_D2y_b;
  Double I_NAI_Py_Fx2y_b = I_NAI_Dxy_D2y_b+ABX*I_NAI_Py_D2y_b;
  Double I_NAI_Pz_Fx2y_b = I_NAI_Dxz_D2y_b+ABX*I_NAI_Pz_D2y_b;
  Double I_NAI_Px_Fxyz_b = I_NAI_Dxz_Dxy_b+ABZ*I_NAI_Px_Dxy_b;
  Double I_NAI_Py_Fxyz_b = I_NAI_Dyz_Dxy_b+ABZ*I_NAI_Py_Dxy_b;
  Double I_NAI_Pz_Fxyz_b = I_NAI_D2z_Dxy_b+ABZ*I_NAI_Pz_Dxy_b;
  Double I_NAI_Px_Fx2z_b = I_NAI_D2x_D2z_b+ABX*I_NAI_Px_D2z_b;
  Double I_NAI_Py_Fx2z_b = I_NAI_Dxy_D2z_b+ABX*I_NAI_Py_D2z_b;
  Double I_NAI_Pz_Fx2z_b = I_NAI_Dxz_D2z_b+ABX*I_NAI_Pz_D2z_b;
  Double I_NAI_Px_F3y_b = I_NAI_Dxy_D2y_b+ABY*I_NAI_Px_D2y_b;
  Double I_NAI_Py_F3y_b = I_NAI_D2y_D2y_b+ABY*I_NAI_Py_D2y_b;
  Double I_NAI_Pz_F3y_b = I_NAI_Dyz_D2y_b+ABY*I_NAI_Pz_D2y_b;
  Double I_NAI_Px_F2yz_b = I_NAI_Dxz_D2y_b+ABZ*I_NAI_Px_D2y_b;
  Double I_NAI_Py_F2yz_b = I_NAI_Dyz_D2y_b+ABZ*I_NAI_Py_D2y_b;
  Double I_NAI_Pz_F2yz_b = I_NAI_D2z_D2y_b+ABZ*I_NAI_Pz_D2y_b;
  Double I_NAI_Px_Fy2z_b = I_NAI_Dxy_D2z_b+ABY*I_NAI_Px_D2z_b;
  Double I_NAI_Py_Fy2z_b = I_NAI_D2y_D2z_b+ABY*I_NAI_Py_D2z_b;
  Double I_NAI_Pz_Fy2z_b = I_NAI_Dyz_D2z_b+ABY*I_NAI_Pz_D2z_b;
  Double I_NAI_Px_F3z_b = I_NAI_Dxz_D2z_b+ABZ*I_NAI_Px_D2z_b;
  Double I_NAI_Py_F3z_b = I_NAI_Dyz_D2z_b+ABZ*I_NAI_Py_D2z_b;
  Double I_NAI_Pz_F3z_b = I_NAI_D2z_D2z_b+ABZ*I_NAI_Pz_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_aa
   * RHS shell quartet name: SQ_NAI_G_S_aa
   ************************************************************/
  Double I_NAI_G4x_Px_aa = I_NAI_H5x_S_aa+ABX*I_NAI_G4x_S_aa;
  Double I_NAI_G3xy_Px_aa = I_NAI_H4xy_S_aa+ABX*I_NAI_G3xy_S_aa;
  Double I_NAI_G3xz_Px_aa = I_NAI_H4xz_S_aa+ABX*I_NAI_G3xz_S_aa;
  Double I_NAI_G2x2y_Px_aa = I_NAI_H3x2y_S_aa+ABX*I_NAI_G2x2y_S_aa;
  Double I_NAI_G2xyz_Px_aa = I_NAI_H3xyz_S_aa+ABX*I_NAI_G2xyz_S_aa;
  Double I_NAI_G2x2z_Px_aa = I_NAI_H3x2z_S_aa+ABX*I_NAI_G2x2z_S_aa;
  Double I_NAI_Gx3y_Px_aa = I_NAI_H2x3y_S_aa+ABX*I_NAI_Gx3y_S_aa;
  Double I_NAI_Gx2yz_Px_aa = I_NAI_H2x2yz_S_aa+ABX*I_NAI_Gx2yz_S_aa;
  Double I_NAI_Gxy2z_Px_aa = I_NAI_H2xy2z_S_aa+ABX*I_NAI_Gxy2z_S_aa;
  Double I_NAI_Gx3z_Px_aa = I_NAI_H2x3z_S_aa+ABX*I_NAI_Gx3z_S_aa;
  Double I_NAI_G4y_Px_aa = I_NAI_Hx4y_S_aa+ABX*I_NAI_G4y_S_aa;
  Double I_NAI_G3yz_Px_aa = I_NAI_Hx3yz_S_aa+ABX*I_NAI_G3yz_S_aa;
  Double I_NAI_G2y2z_Px_aa = I_NAI_Hx2y2z_S_aa+ABX*I_NAI_G2y2z_S_aa;
  Double I_NAI_Gy3z_Px_aa = I_NAI_Hxy3z_S_aa+ABX*I_NAI_Gy3z_S_aa;
  Double I_NAI_G4z_Px_aa = I_NAI_Hx4z_S_aa+ABX*I_NAI_G4z_S_aa;
  Double I_NAI_G4x_Py_aa = I_NAI_H4xy_S_aa+ABY*I_NAI_G4x_S_aa;
  Double I_NAI_G3xy_Py_aa = I_NAI_H3x2y_S_aa+ABY*I_NAI_G3xy_S_aa;
  Double I_NAI_G3xz_Py_aa = I_NAI_H3xyz_S_aa+ABY*I_NAI_G3xz_S_aa;
  Double I_NAI_G2x2y_Py_aa = I_NAI_H2x3y_S_aa+ABY*I_NAI_G2x2y_S_aa;
  Double I_NAI_G2xyz_Py_aa = I_NAI_H2x2yz_S_aa+ABY*I_NAI_G2xyz_S_aa;
  Double I_NAI_G2x2z_Py_aa = I_NAI_H2xy2z_S_aa+ABY*I_NAI_G2x2z_S_aa;
  Double I_NAI_Gx3y_Py_aa = I_NAI_Hx4y_S_aa+ABY*I_NAI_Gx3y_S_aa;
  Double I_NAI_Gx2yz_Py_aa = I_NAI_Hx3yz_S_aa+ABY*I_NAI_Gx2yz_S_aa;
  Double I_NAI_Gxy2z_Py_aa = I_NAI_Hx2y2z_S_aa+ABY*I_NAI_Gxy2z_S_aa;
  Double I_NAI_Gx3z_Py_aa = I_NAI_Hxy3z_S_aa+ABY*I_NAI_Gx3z_S_aa;
  Double I_NAI_G4y_Py_aa = I_NAI_H5y_S_aa+ABY*I_NAI_G4y_S_aa;
  Double I_NAI_G3yz_Py_aa = I_NAI_H4yz_S_aa+ABY*I_NAI_G3yz_S_aa;
  Double I_NAI_G2y2z_Py_aa = I_NAI_H3y2z_S_aa+ABY*I_NAI_G2y2z_S_aa;
  Double I_NAI_Gy3z_Py_aa = I_NAI_H2y3z_S_aa+ABY*I_NAI_Gy3z_S_aa;
  Double I_NAI_G4z_Py_aa = I_NAI_Hy4z_S_aa+ABY*I_NAI_G4z_S_aa;
  Double I_NAI_G4x_Pz_aa = I_NAI_H4xz_S_aa+ABZ*I_NAI_G4x_S_aa;
  Double I_NAI_G3xy_Pz_aa = I_NAI_H3xyz_S_aa+ABZ*I_NAI_G3xy_S_aa;
  Double I_NAI_G3xz_Pz_aa = I_NAI_H3x2z_S_aa+ABZ*I_NAI_G3xz_S_aa;
  Double I_NAI_G2x2y_Pz_aa = I_NAI_H2x2yz_S_aa+ABZ*I_NAI_G2x2y_S_aa;
  Double I_NAI_G2xyz_Pz_aa = I_NAI_H2xy2z_S_aa+ABZ*I_NAI_G2xyz_S_aa;
  Double I_NAI_G2x2z_Pz_aa = I_NAI_H2x3z_S_aa+ABZ*I_NAI_G2x2z_S_aa;
  Double I_NAI_Gx3y_Pz_aa = I_NAI_Hx3yz_S_aa+ABZ*I_NAI_Gx3y_S_aa;
  Double I_NAI_Gx2yz_Pz_aa = I_NAI_Hx2y2z_S_aa+ABZ*I_NAI_Gx2yz_S_aa;
  Double I_NAI_Gxy2z_Pz_aa = I_NAI_Hxy3z_S_aa+ABZ*I_NAI_Gxy2z_S_aa;
  Double I_NAI_Gx3z_Pz_aa = I_NAI_Hx4z_S_aa+ABZ*I_NAI_Gx3z_S_aa;
  Double I_NAI_G4y_Pz_aa = I_NAI_H4yz_S_aa+ABZ*I_NAI_G4y_S_aa;
  Double I_NAI_G3yz_Pz_aa = I_NAI_H3y2z_S_aa+ABZ*I_NAI_G3yz_S_aa;
  Double I_NAI_G2y2z_Pz_aa = I_NAI_H2y3z_S_aa+ABZ*I_NAI_G2y2z_S_aa;
  Double I_NAI_Gy3z_Pz_aa = I_NAI_Hy4z_S_aa+ABZ*I_NAI_Gy3z_S_aa;
  Double I_NAI_G4z_Pz_aa = I_NAI_H5z_S_aa+ABZ*I_NAI_G4z_S_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 7 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_aa
   * RHS shell quartet name: SQ_NAI_H_S_aa
   ************************************************************/
  Double I_NAI_H5x_Px_aa = I_NAI_I6x_S_aa+ABX*I_NAI_H5x_S_aa;
  Double I_NAI_H4xy_Px_aa = I_NAI_I5xy_S_aa+ABX*I_NAI_H4xy_S_aa;
  Double I_NAI_H4xz_Px_aa = I_NAI_I5xz_S_aa+ABX*I_NAI_H4xz_S_aa;
  Double I_NAI_H3x2y_Px_aa = I_NAI_I4x2y_S_aa+ABX*I_NAI_H3x2y_S_aa;
  Double I_NAI_H3xyz_Px_aa = I_NAI_I4xyz_S_aa+ABX*I_NAI_H3xyz_S_aa;
  Double I_NAI_H3x2z_Px_aa = I_NAI_I4x2z_S_aa+ABX*I_NAI_H3x2z_S_aa;
  Double I_NAI_H2x3y_Px_aa = I_NAI_I3x3y_S_aa+ABX*I_NAI_H2x3y_S_aa;
  Double I_NAI_H2x2yz_Px_aa = I_NAI_I3x2yz_S_aa+ABX*I_NAI_H2x2yz_S_aa;
  Double I_NAI_H2xy2z_Px_aa = I_NAI_I3xy2z_S_aa+ABX*I_NAI_H2xy2z_S_aa;
  Double I_NAI_H2x3z_Px_aa = I_NAI_I3x3z_S_aa+ABX*I_NAI_H2x3z_S_aa;
  Double I_NAI_Hx4y_Px_aa = I_NAI_I2x4y_S_aa+ABX*I_NAI_Hx4y_S_aa;
  Double I_NAI_Hx3yz_Px_aa = I_NAI_I2x3yz_S_aa+ABX*I_NAI_Hx3yz_S_aa;
  Double I_NAI_Hx2y2z_Px_aa = I_NAI_I2x2y2z_S_aa+ABX*I_NAI_Hx2y2z_S_aa;
  Double I_NAI_Hxy3z_Px_aa = I_NAI_I2xy3z_S_aa+ABX*I_NAI_Hxy3z_S_aa;
  Double I_NAI_Hx4z_Px_aa = I_NAI_I2x4z_S_aa+ABX*I_NAI_Hx4z_S_aa;
  Double I_NAI_H5y_Px_aa = I_NAI_Ix5y_S_aa+ABX*I_NAI_H5y_S_aa;
  Double I_NAI_H4yz_Px_aa = I_NAI_Ix4yz_S_aa+ABX*I_NAI_H4yz_S_aa;
  Double I_NAI_H3y2z_Px_aa = I_NAI_Ix3y2z_S_aa+ABX*I_NAI_H3y2z_S_aa;
  Double I_NAI_H2y3z_Px_aa = I_NAI_Ix2y3z_S_aa+ABX*I_NAI_H2y3z_S_aa;
  Double I_NAI_Hy4z_Px_aa = I_NAI_Ixy4z_S_aa+ABX*I_NAI_Hy4z_S_aa;
  Double I_NAI_H5z_Px_aa = I_NAI_Ix5z_S_aa+ABX*I_NAI_H5z_S_aa;
  Double I_NAI_H4xy_Py_aa = I_NAI_I4x2y_S_aa+ABY*I_NAI_H4xy_S_aa;
  Double I_NAI_H4xz_Py_aa = I_NAI_I4xyz_S_aa+ABY*I_NAI_H4xz_S_aa;
  Double I_NAI_H3x2y_Py_aa = I_NAI_I3x3y_S_aa+ABY*I_NAI_H3x2y_S_aa;
  Double I_NAI_H3xyz_Py_aa = I_NAI_I3x2yz_S_aa+ABY*I_NAI_H3xyz_S_aa;
  Double I_NAI_H3x2z_Py_aa = I_NAI_I3xy2z_S_aa+ABY*I_NAI_H3x2z_S_aa;
  Double I_NAI_H2x3y_Py_aa = I_NAI_I2x4y_S_aa+ABY*I_NAI_H2x3y_S_aa;
  Double I_NAI_H2x2yz_Py_aa = I_NAI_I2x3yz_S_aa+ABY*I_NAI_H2x2yz_S_aa;
  Double I_NAI_H2xy2z_Py_aa = I_NAI_I2x2y2z_S_aa+ABY*I_NAI_H2xy2z_S_aa;
  Double I_NAI_H2x3z_Py_aa = I_NAI_I2xy3z_S_aa+ABY*I_NAI_H2x3z_S_aa;
  Double I_NAI_Hx4y_Py_aa = I_NAI_Ix5y_S_aa+ABY*I_NAI_Hx4y_S_aa;
  Double I_NAI_Hx3yz_Py_aa = I_NAI_Ix4yz_S_aa+ABY*I_NAI_Hx3yz_S_aa;
  Double I_NAI_Hx2y2z_Py_aa = I_NAI_Ix3y2z_S_aa+ABY*I_NAI_Hx2y2z_S_aa;
  Double I_NAI_Hxy3z_Py_aa = I_NAI_Ix2y3z_S_aa+ABY*I_NAI_Hxy3z_S_aa;
  Double I_NAI_Hx4z_Py_aa = I_NAI_Ixy4z_S_aa+ABY*I_NAI_Hx4z_S_aa;
  Double I_NAI_H5y_Py_aa = I_NAI_I6y_S_aa+ABY*I_NAI_H5y_S_aa;
  Double I_NAI_H4yz_Py_aa = I_NAI_I5yz_S_aa+ABY*I_NAI_H4yz_S_aa;
  Double I_NAI_H3y2z_Py_aa = I_NAI_I4y2z_S_aa+ABY*I_NAI_H3y2z_S_aa;
  Double I_NAI_H2y3z_Py_aa = I_NAI_I3y3z_S_aa+ABY*I_NAI_H2y3z_S_aa;
  Double I_NAI_Hy4z_Py_aa = I_NAI_I2y4z_S_aa+ABY*I_NAI_Hy4z_S_aa;
  Double I_NAI_H5z_Py_aa = I_NAI_Iy5z_S_aa+ABY*I_NAI_H5z_S_aa;
  Double I_NAI_H4xz_Pz_aa = I_NAI_I4x2z_S_aa+ABZ*I_NAI_H4xz_S_aa;
  Double I_NAI_H3xyz_Pz_aa = I_NAI_I3xy2z_S_aa+ABZ*I_NAI_H3xyz_S_aa;
  Double I_NAI_H3x2z_Pz_aa = I_NAI_I3x3z_S_aa+ABZ*I_NAI_H3x2z_S_aa;
  Double I_NAI_H2x2yz_Pz_aa = I_NAI_I2x2y2z_S_aa+ABZ*I_NAI_H2x2yz_S_aa;
  Double I_NAI_H2xy2z_Pz_aa = I_NAI_I2xy3z_S_aa+ABZ*I_NAI_H2xy2z_S_aa;
  Double I_NAI_H2x3z_Pz_aa = I_NAI_I2x4z_S_aa+ABZ*I_NAI_H2x3z_S_aa;
  Double I_NAI_Hx3yz_Pz_aa = I_NAI_Ix3y2z_S_aa+ABZ*I_NAI_Hx3yz_S_aa;
  Double I_NAI_Hx2y2z_Pz_aa = I_NAI_Ix2y3z_S_aa+ABZ*I_NAI_Hx2y2z_S_aa;
  Double I_NAI_Hxy3z_Pz_aa = I_NAI_Ixy4z_S_aa+ABZ*I_NAI_Hxy3z_S_aa;
  Double I_NAI_Hx4z_Pz_aa = I_NAI_Ix5z_S_aa+ABZ*I_NAI_Hx4z_S_aa;
  Double I_NAI_H4yz_Pz_aa = I_NAI_I4y2z_S_aa+ABZ*I_NAI_H4yz_S_aa;
  Double I_NAI_H3y2z_Pz_aa = I_NAI_I3y3z_S_aa+ABZ*I_NAI_H3y2z_S_aa;
  Double I_NAI_H2y3z_Pz_aa = I_NAI_I2y4z_S_aa+ABZ*I_NAI_H2y3z_S_aa;
  Double I_NAI_Hy4z_Pz_aa = I_NAI_Iy5z_S_aa+ABZ*I_NAI_Hy4z_S_aa;
  Double I_NAI_H5z_Pz_aa = I_NAI_I6z_S_aa+ABZ*I_NAI_H5z_S_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_aa
   * RHS shell quartet name: SQ_NAI_G_P_aa
   ************************************************************/
  Double I_NAI_G4x_D2x_aa = I_NAI_H5x_Px_aa+ABX*I_NAI_G4x_Px_aa;
  Double I_NAI_G3xy_D2x_aa = I_NAI_H4xy_Px_aa+ABX*I_NAI_G3xy_Px_aa;
  Double I_NAI_G3xz_D2x_aa = I_NAI_H4xz_Px_aa+ABX*I_NAI_G3xz_Px_aa;
  Double I_NAI_G2x2y_D2x_aa = I_NAI_H3x2y_Px_aa+ABX*I_NAI_G2x2y_Px_aa;
  Double I_NAI_G2xyz_D2x_aa = I_NAI_H3xyz_Px_aa+ABX*I_NAI_G2xyz_Px_aa;
  Double I_NAI_G2x2z_D2x_aa = I_NAI_H3x2z_Px_aa+ABX*I_NAI_G2x2z_Px_aa;
  Double I_NAI_Gx3y_D2x_aa = I_NAI_H2x3y_Px_aa+ABX*I_NAI_Gx3y_Px_aa;
  Double I_NAI_Gx2yz_D2x_aa = I_NAI_H2x2yz_Px_aa+ABX*I_NAI_Gx2yz_Px_aa;
  Double I_NAI_Gxy2z_D2x_aa = I_NAI_H2xy2z_Px_aa+ABX*I_NAI_Gxy2z_Px_aa;
  Double I_NAI_Gx3z_D2x_aa = I_NAI_H2x3z_Px_aa+ABX*I_NAI_Gx3z_Px_aa;
  Double I_NAI_G4y_D2x_aa = I_NAI_Hx4y_Px_aa+ABX*I_NAI_G4y_Px_aa;
  Double I_NAI_G3yz_D2x_aa = I_NAI_Hx3yz_Px_aa+ABX*I_NAI_G3yz_Px_aa;
  Double I_NAI_G2y2z_D2x_aa = I_NAI_Hx2y2z_Px_aa+ABX*I_NAI_G2y2z_Px_aa;
  Double I_NAI_Gy3z_D2x_aa = I_NAI_Hxy3z_Px_aa+ABX*I_NAI_Gy3z_Px_aa;
  Double I_NAI_G4z_D2x_aa = I_NAI_Hx4z_Px_aa+ABX*I_NAI_G4z_Px_aa;
  Double I_NAI_G4x_Dxy_aa = I_NAI_H4xy_Px_aa+ABY*I_NAI_G4x_Px_aa;
  Double I_NAI_G3xy_Dxy_aa = I_NAI_H3x2y_Px_aa+ABY*I_NAI_G3xy_Px_aa;
  Double I_NAI_G3xz_Dxy_aa = I_NAI_H3xyz_Px_aa+ABY*I_NAI_G3xz_Px_aa;
  Double I_NAI_G2x2y_Dxy_aa = I_NAI_H2x3y_Px_aa+ABY*I_NAI_G2x2y_Px_aa;
  Double I_NAI_G2xyz_Dxy_aa = I_NAI_H2x2yz_Px_aa+ABY*I_NAI_G2xyz_Px_aa;
  Double I_NAI_G2x2z_Dxy_aa = I_NAI_H2xy2z_Px_aa+ABY*I_NAI_G2x2z_Px_aa;
  Double I_NAI_Gx3y_Dxy_aa = I_NAI_Hx4y_Px_aa+ABY*I_NAI_Gx3y_Px_aa;
  Double I_NAI_Gx2yz_Dxy_aa = I_NAI_Hx3yz_Px_aa+ABY*I_NAI_Gx2yz_Px_aa;
  Double I_NAI_Gxy2z_Dxy_aa = I_NAI_Hx2y2z_Px_aa+ABY*I_NAI_Gxy2z_Px_aa;
  Double I_NAI_Gx3z_Dxy_aa = I_NAI_Hxy3z_Px_aa+ABY*I_NAI_Gx3z_Px_aa;
  Double I_NAI_G4y_Dxy_aa = I_NAI_H5y_Px_aa+ABY*I_NAI_G4y_Px_aa;
  Double I_NAI_G3yz_Dxy_aa = I_NAI_H4yz_Px_aa+ABY*I_NAI_G3yz_Px_aa;
  Double I_NAI_G2y2z_Dxy_aa = I_NAI_H3y2z_Px_aa+ABY*I_NAI_G2y2z_Px_aa;
  Double I_NAI_Gy3z_Dxy_aa = I_NAI_H2y3z_Px_aa+ABY*I_NAI_Gy3z_Px_aa;
  Double I_NAI_G4z_Dxy_aa = I_NAI_Hy4z_Px_aa+ABY*I_NAI_G4z_Px_aa;
  Double I_NAI_G4x_Dxz_aa = I_NAI_H4xz_Px_aa+ABZ*I_NAI_G4x_Px_aa;
  Double I_NAI_G3xy_Dxz_aa = I_NAI_H3xyz_Px_aa+ABZ*I_NAI_G3xy_Px_aa;
  Double I_NAI_G3xz_Dxz_aa = I_NAI_H3x2z_Px_aa+ABZ*I_NAI_G3xz_Px_aa;
  Double I_NAI_G2x2y_Dxz_aa = I_NAI_H2x2yz_Px_aa+ABZ*I_NAI_G2x2y_Px_aa;
  Double I_NAI_G2xyz_Dxz_aa = I_NAI_H2xy2z_Px_aa+ABZ*I_NAI_G2xyz_Px_aa;
  Double I_NAI_G2x2z_Dxz_aa = I_NAI_H2x3z_Px_aa+ABZ*I_NAI_G2x2z_Px_aa;
  Double I_NAI_Gx3y_Dxz_aa = I_NAI_Hx3yz_Px_aa+ABZ*I_NAI_Gx3y_Px_aa;
  Double I_NAI_Gx2yz_Dxz_aa = I_NAI_Hx2y2z_Px_aa+ABZ*I_NAI_Gx2yz_Px_aa;
  Double I_NAI_Gxy2z_Dxz_aa = I_NAI_Hxy3z_Px_aa+ABZ*I_NAI_Gxy2z_Px_aa;
  Double I_NAI_Gx3z_Dxz_aa = I_NAI_Hx4z_Px_aa+ABZ*I_NAI_Gx3z_Px_aa;
  Double I_NAI_G4y_Dxz_aa = I_NAI_H4yz_Px_aa+ABZ*I_NAI_G4y_Px_aa;
  Double I_NAI_G3yz_Dxz_aa = I_NAI_H3y2z_Px_aa+ABZ*I_NAI_G3yz_Px_aa;
  Double I_NAI_G2y2z_Dxz_aa = I_NAI_H2y3z_Px_aa+ABZ*I_NAI_G2y2z_Px_aa;
  Double I_NAI_Gy3z_Dxz_aa = I_NAI_Hy4z_Px_aa+ABZ*I_NAI_Gy3z_Px_aa;
  Double I_NAI_G4z_Dxz_aa = I_NAI_H5z_Px_aa+ABZ*I_NAI_G4z_Px_aa;
  Double I_NAI_G4x_D2y_aa = I_NAI_H4xy_Py_aa+ABY*I_NAI_G4x_Py_aa;
  Double I_NAI_G3xy_D2y_aa = I_NAI_H3x2y_Py_aa+ABY*I_NAI_G3xy_Py_aa;
  Double I_NAI_G3xz_D2y_aa = I_NAI_H3xyz_Py_aa+ABY*I_NAI_G3xz_Py_aa;
  Double I_NAI_G2x2y_D2y_aa = I_NAI_H2x3y_Py_aa+ABY*I_NAI_G2x2y_Py_aa;
  Double I_NAI_G2xyz_D2y_aa = I_NAI_H2x2yz_Py_aa+ABY*I_NAI_G2xyz_Py_aa;
  Double I_NAI_G2x2z_D2y_aa = I_NAI_H2xy2z_Py_aa+ABY*I_NAI_G2x2z_Py_aa;
  Double I_NAI_Gx3y_D2y_aa = I_NAI_Hx4y_Py_aa+ABY*I_NAI_Gx3y_Py_aa;
  Double I_NAI_Gx2yz_D2y_aa = I_NAI_Hx3yz_Py_aa+ABY*I_NAI_Gx2yz_Py_aa;
  Double I_NAI_Gxy2z_D2y_aa = I_NAI_Hx2y2z_Py_aa+ABY*I_NAI_Gxy2z_Py_aa;
  Double I_NAI_Gx3z_D2y_aa = I_NAI_Hxy3z_Py_aa+ABY*I_NAI_Gx3z_Py_aa;
  Double I_NAI_G4y_D2y_aa = I_NAI_H5y_Py_aa+ABY*I_NAI_G4y_Py_aa;
  Double I_NAI_G3yz_D2y_aa = I_NAI_H4yz_Py_aa+ABY*I_NAI_G3yz_Py_aa;
  Double I_NAI_G2y2z_D2y_aa = I_NAI_H3y2z_Py_aa+ABY*I_NAI_G2y2z_Py_aa;
  Double I_NAI_Gy3z_D2y_aa = I_NAI_H2y3z_Py_aa+ABY*I_NAI_Gy3z_Py_aa;
  Double I_NAI_G4z_D2y_aa = I_NAI_Hy4z_Py_aa+ABY*I_NAI_G4z_Py_aa;
  Double I_NAI_G4x_Dyz_aa = I_NAI_H4xz_Py_aa+ABZ*I_NAI_G4x_Py_aa;
  Double I_NAI_G3xy_Dyz_aa = I_NAI_H3xyz_Py_aa+ABZ*I_NAI_G3xy_Py_aa;
  Double I_NAI_G3xz_Dyz_aa = I_NAI_H3x2z_Py_aa+ABZ*I_NAI_G3xz_Py_aa;
  Double I_NAI_G2x2y_Dyz_aa = I_NAI_H2x2yz_Py_aa+ABZ*I_NAI_G2x2y_Py_aa;
  Double I_NAI_G2xyz_Dyz_aa = I_NAI_H2xy2z_Py_aa+ABZ*I_NAI_G2xyz_Py_aa;
  Double I_NAI_G2x2z_Dyz_aa = I_NAI_H2x3z_Py_aa+ABZ*I_NAI_G2x2z_Py_aa;
  Double I_NAI_Gx3y_Dyz_aa = I_NAI_Hx3yz_Py_aa+ABZ*I_NAI_Gx3y_Py_aa;
  Double I_NAI_Gx2yz_Dyz_aa = I_NAI_Hx2y2z_Py_aa+ABZ*I_NAI_Gx2yz_Py_aa;
  Double I_NAI_Gxy2z_Dyz_aa = I_NAI_Hxy3z_Py_aa+ABZ*I_NAI_Gxy2z_Py_aa;
  Double I_NAI_Gx3z_Dyz_aa = I_NAI_Hx4z_Py_aa+ABZ*I_NAI_Gx3z_Py_aa;
  Double I_NAI_G4y_Dyz_aa = I_NAI_H4yz_Py_aa+ABZ*I_NAI_G4y_Py_aa;
  Double I_NAI_G3yz_Dyz_aa = I_NAI_H3y2z_Py_aa+ABZ*I_NAI_G3yz_Py_aa;
  Double I_NAI_G2y2z_Dyz_aa = I_NAI_H2y3z_Py_aa+ABZ*I_NAI_G2y2z_Py_aa;
  Double I_NAI_Gy3z_Dyz_aa = I_NAI_Hy4z_Py_aa+ABZ*I_NAI_Gy3z_Py_aa;
  Double I_NAI_G4z_Dyz_aa = I_NAI_H5z_Py_aa+ABZ*I_NAI_G4z_Py_aa;
  Double I_NAI_G4x_D2z_aa = I_NAI_H4xz_Pz_aa+ABZ*I_NAI_G4x_Pz_aa;
  Double I_NAI_G3xy_D2z_aa = I_NAI_H3xyz_Pz_aa+ABZ*I_NAI_G3xy_Pz_aa;
  Double I_NAI_G3xz_D2z_aa = I_NAI_H3x2z_Pz_aa+ABZ*I_NAI_G3xz_Pz_aa;
  Double I_NAI_G2x2y_D2z_aa = I_NAI_H2x2yz_Pz_aa+ABZ*I_NAI_G2x2y_Pz_aa;
  Double I_NAI_G2xyz_D2z_aa = I_NAI_H2xy2z_Pz_aa+ABZ*I_NAI_G2xyz_Pz_aa;
  Double I_NAI_G2x2z_D2z_aa = I_NAI_H2x3z_Pz_aa+ABZ*I_NAI_G2x2z_Pz_aa;
  Double I_NAI_Gx3y_D2z_aa = I_NAI_Hx3yz_Pz_aa+ABZ*I_NAI_Gx3y_Pz_aa;
  Double I_NAI_Gx2yz_D2z_aa = I_NAI_Hx2y2z_Pz_aa+ABZ*I_NAI_Gx2yz_Pz_aa;
  Double I_NAI_Gxy2z_D2z_aa = I_NAI_Hxy3z_Pz_aa+ABZ*I_NAI_Gxy2z_Pz_aa;
  Double I_NAI_Gx3z_D2z_aa = I_NAI_Hx4z_Pz_aa+ABZ*I_NAI_Gx3z_Pz_aa;
  Double I_NAI_G4y_D2z_aa = I_NAI_H4yz_Pz_aa+ABZ*I_NAI_G4y_Pz_aa;
  Double I_NAI_G3yz_D2z_aa = I_NAI_H3y2z_Pz_aa+ABZ*I_NAI_G3yz_Pz_aa;
  Double I_NAI_G2y2z_D2z_aa = I_NAI_H2y3z_Pz_aa+ABZ*I_NAI_G2y2z_Pz_aa;
  Double I_NAI_Gy3z_D2z_aa = I_NAI_Hy4z_Pz_aa+ABZ*I_NAI_Gy3z_Pz_aa;
  Double I_NAI_G4z_D2z_aa = I_NAI_H5z_Pz_aa+ABZ*I_NAI_G4z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_ab
   * RHS shell quartet name: SQ_NAI_F_S_ab
   ************************************************************/
  Double I_NAI_F3x_Px_ab = I_NAI_G4x_S_ab+ABX*I_NAI_F3x_S_ab;
  Double I_NAI_F2xy_Px_ab = I_NAI_G3xy_S_ab+ABX*I_NAI_F2xy_S_ab;
  Double I_NAI_F2xz_Px_ab = I_NAI_G3xz_S_ab+ABX*I_NAI_F2xz_S_ab;
  Double I_NAI_Fx2y_Px_ab = I_NAI_G2x2y_S_ab+ABX*I_NAI_Fx2y_S_ab;
  Double I_NAI_Fxyz_Px_ab = I_NAI_G2xyz_S_ab+ABX*I_NAI_Fxyz_S_ab;
  Double I_NAI_Fx2z_Px_ab = I_NAI_G2x2z_S_ab+ABX*I_NAI_Fx2z_S_ab;
  Double I_NAI_F3y_Px_ab = I_NAI_Gx3y_S_ab+ABX*I_NAI_F3y_S_ab;
  Double I_NAI_F2yz_Px_ab = I_NAI_Gx2yz_S_ab+ABX*I_NAI_F2yz_S_ab;
  Double I_NAI_Fy2z_Px_ab = I_NAI_Gxy2z_S_ab+ABX*I_NAI_Fy2z_S_ab;
  Double I_NAI_F3z_Px_ab = I_NAI_Gx3z_S_ab+ABX*I_NAI_F3z_S_ab;
  Double I_NAI_F3x_Py_ab = I_NAI_G3xy_S_ab+ABY*I_NAI_F3x_S_ab;
  Double I_NAI_F2xy_Py_ab = I_NAI_G2x2y_S_ab+ABY*I_NAI_F2xy_S_ab;
  Double I_NAI_F2xz_Py_ab = I_NAI_G2xyz_S_ab+ABY*I_NAI_F2xz_S_ab;
  Double I_NAI_Fx2y_Py_ab = I_NAI_Gx3y_S_ab+ABY*I_NAI_Fx2y_S_ab;
  Double I_NAI_Fxyz_Py_ab = I_NAI_Gx2yz_S_ab+ABY*I_NAI_Fxyz_S_ab;
  Double I_NAI_Fx2z_Py_ab = I_NAI_Gxy2z_S_ab+ABY*I_NAI_Fx2z_S_ab;
  Double I_NAI_F3y_Py_ab = I_NAI_G4y_S_ab+ABY*I_NAI_F3y_S_ab;
  Double I_NAI_F2yz_Py_ab = I_NAI_G3yz_S_ab+ABY*I_NAI_F2yz_S_ab;
  Double I_NAI_Fy2z_Py_ab = I_NAI_G2y2z_S_ab+ABY*I_NAI_Fy2z_S_ab;
  Double I_NAI_F3z_Py_ab = I_NAI_Gy3z_S_ab+ABY*I_NAI_F3z_S_ab;
  Double I_NAI_F3x_Pz_ab = I_NAI_G3xz_S_ab+ABZ*I_NAI_F3x_S_ab;
  Double I_NAI_F2xy_Pz_ab = I_NAI_G2xyz_S_ab+ABZ*I_NAI_F2xy_S_ab;
  Double I_NAI_F2xz_Pz_ab = I_NAI_G2x2z_S_ab+ABZ*I_NAI_F2xz_S_ab;
  Double I_NAI_Fx2y_Pz_ab = I_NAI_Gx2yz_S_ab+ABZ*I_NAI_Fx2y_S_ab;
  Double I_NAI_Fxyz_Pz_ab = I_NAI_Gxy2z_S_ab+ABZ*I_NAI_Fxyz_S_ab;
  Double I_NAI_Fx2z_Pz_ab = I_NAI_Gx3z_S_ab+ABZ*I_NAI_Fx2z_S_ab;
  Double I_NAI_F3y_Pz_ab = I_NAI_G3yz_S_ab+ABZ*I_NAI_F3y_S_ab;
  Double I_NAI_F2yz_Pz_ab = I_NAI_G2y2z_S_ab+ABZ*I_NAI_F2yz_S_ab;
  Double I_NAI_Fy2z_Pz_ab = I_NAI_Gy3z_S_ab+ABZ*I_NAI_Fy2z_S_ab;
  Double I_NAI_F3z_Pz_ab = I_NAI_G4z_S_ab+ABZ*I_NAI_F3z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_ab
   * RHS shell quartet name: SQ_NAI_G_S_ab
   ************************************************************/
  Double I_NAI_G4x_Px_ab = I_NAI_H5x_S_ab+ABX*I_NAI_G4x_S_ab;
  Double I_NAI_G3xy_Px_ab = I_NAI_H4xy_S_ab+ABX*I_NAI_G3xy_S_ab;
  Double I_NAI_G3xz_Px_ab = I_NAI_H4xz_S_ab+ABX*I_NAI_G3xz_S_ab;
  Double I_NAI_G2x2y_Px_ab = I_NAI_H3x2y_S_ab+ABX*I_NAI_G2x2y_S_ab;
  Double I_NAI_G2xyz_Px_ab = I_NAI_H3xyz_S_ab+ABX*I_NAI_G2xyz_S_ab;
  Double I_NAI_G2x2z_Px_ab = I_NAI_H3x2z_S_ab+ABX*I_NAI_G2x2z_S_ab;
  Double I_NAI_Gx3y_Px_ab = I_NAI_H2x3y_S_ab+ABX*I_NAI_Gx3y_S_ab;
  Double I_NAI_Gx2yz_Px_ab = I_NAI_H2x2yz_S_ab+ABX*I_NAI_Gx2yz_S_ab;
  Double I_NAI_Gxy2z_Px_ab = I_NAI_H2xy2z_S_ab+ABX*I_NAI_Gxy2z_S_ab;
  Double I_NAI_Gx3z_Px_ab = I_NAI_H2x3z_S_ab+ABX*I_NAI_Gx3z_S_ab;
  Double I_NAI_G4y_Px_ab = I_NAI_Hx4y_S_ab+ABX*I_NAI_G4y_S_ab;
  Double I_NAI_G3yz_Px_ab = I_NAI_Hx3yz_S_ab+ABX*I_NAI_G3yz_S_ab;
  Double I_NAI_G2y2z_Px_ab = I_NAI_Hx2y2z_S_ab+ABX*I_NAI_G2y2z_S_ab;
  Double I_NAI_Gy3z_Px_ab = I_NAI_Hxy3z_S_ab+ABX*I_NAI_Gy3z_S_ab;
  Double I_NAI_G4z_Px_ab = I_NAI_Hx4z_S_ab+ABX*I_NAI_G4z_S_ab;
  Double I_NAI_G4x_Py_ab = I_NAI_H4xy_S_ab+ABY*I_NAI_G4x_S_ab;
  Double I_NAI_G3xy_Py_ab = I_NAI_H3x2y_S_ab+ABY*I_NAI_G3xy_S_ab;
  Double I_NAI_G3xz_Py_ab = I_NAI_H3xyz_S_ab+ABY*I_NAI_G3xz_S_ab;
  Double I_NAI_G2x2y_Py_ab = I_NAI_H2x3y_S_ab+ABY*I_NAI_G2x2y_S_ab;
  Double I_NAI_G2xyz_Py_ab = I_NAI_H2x2yz_S_ab+ABY*I_NAI_G2xyz_S_ab;
  Double I_NAI_G2x2z_Py_ab = I_NAI_H2xy2z_S_ab+ABY*I_NAI_G2x2z_S_ab;
  Double I_NAI_Gx3y_Py_ab = I_NAI_Hx4y_S_ab+ABY*I_NAI_Gx3y_S_ab;
  Double I_NAI_Gx2yz_Py_ab = I_NAI_Hx3yz_S_ab+ABY*I_NAI_Gx2yz_S_ab;
  Double I_NAI_Gxy2z_Py_ab = I_NAI_Hx2y2z_S_ab+ABY*I_NAI_Gxy2z_S_ab;
  Double I_NAI_Gx3z_Py_ab = I_NAI_Hxy3z_S_ab+ABY*I_NAI_Gx3z_S_ab;
  Double I_NAI_G4y_Py_ab = I_NAI_H5y_S_ab+ABY*I_NAI_G4y_S_ab;
  Double I_NAI_G3yz_Py_ab = I_NAI_H4yz_S_ab+ABY*I_NAI_G3yz_S_ab;
  Double I_NAI_G2y2z_Py_ab = I_NAI_H3y2z_S_ab+ABY*I_NAI_G2y2z_S_ab;
  Double I_NAI_Gy3z_Py_ab = I_NAI_H2y3z_S_ab+ABY*I_NAI_Gy3z_S_ab;
  Double I_NAI_G4z_Py_ab = I_NAI_Hy4z_S_ab+ABY*I_NAI_G4z_S_ab;
  Double I_NAI_G4x_Pz_ab = I_NAI_H4xz_S_ab+ABZ*I_NAI_G4x_S_ab;
  Double I_NAI_G3xy_Pz_ab = I_NAI_H3xyz_S_ab+ABZ*I_NAI_G3xy_S_ab;
  Double I_NAI_G3xz_Pz_ab = I_NAI_H3x2z_S_ab+ABZ*I_NAI_G3xz_S_ab;
  Double I_NAI_G2x2y_Pz_ab = I_NAI_H2x2yz_S_ab+ABZ*I_NAI_G2x2y_S_ab;
  Double I_NAI_G2xyz_Pz_ab = I_NAI_H2xy2z_S_ab+ABZ*I_NAI_G2xyz_S_ab;
  Double I_NAI_G2x2z_Pz_ab = I_NAI_H2x3z_S_ab+ABZ*I_NAI_G2x2z_S_ab;
  Double I_NAI_Gx3y_Pz_ab = I_NAI_Hx3yz_S_ab+ABZ*I_NAI_Gx3y_S_ab;
  Double I_NAI_Gx2yz_Pz_ab = I_NAI_Hx2y2z_S_ab+ABZ*I_NAI_Gx2yz_S_ab;
  Double I_NAI_Gxy2z_Pz_ab = I_NAI_Hxy3z_S_ab+ABZ*I_NAI_Gxy2z_S_ab;
  Double I_NAI_Gx3z_Pz_ab = I_NAI_Hx4z_S_ab+ABZ*I_NAI_Gx3z_S_ab;
  Double I_NAI_G4y_Pz_ab = I_NAI_H4yz_S_ab+ABZ*I_NAI_G4y_S_ab;
  Double I_NAI_G3yz_Pz_ab = I_NAI_H3y2z_S_ab+ABZ*I_NAI_G3yz_S_ab;
  Double I_NAI_G2y2z_Pz_ab = I_NAI_H2y3z_S_ab+ABZ*I_NAI_G2y2z_S_ab;
  Double I_NAI_Gy3z_Pz_ab = I_NAI_Hy4z_S_ab+ABZ*I_NAI_Gy3z_S_ab;
  Double I_NAI_G4z_Pz_ab = I_NAI_H5z_S_ab+ABZ*I_NAI_G4z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_ab
   * RHS shell quartet name: SQ_NAI_F_P_ab
   ************************************************************/
  Double I_NAI_F3x_D2x_ab = I_NAI_G4x_Px_ab+ABX*I_NAI_F3x_Px_ab;
  Double I_NAI_F2xy_D2x_ab = I_NAI_G3xy_Px_ab+ABX*I_NAI_F2xy_Px_ab;
  Double I_NAI_F2xz_D2x_ab = I_NAI_G3xz_Px_ab+ABX*I_NAI_F2xz_Px_ab;
  Double I_NAI_Fx2y_D2x_ab = I_NAI_G2x2y_Px_ab+ABX*I_NAI_Fx2y_Px_ab;
  Double I_NAI_Fxyz_D2x_ab = I_NAI_G2xyz_Px_ab+ABX*I_NAI_Fxyz_Px_ab;
  Double I_NAI_Fx2z_D2x_ab = I_NAI_G2x2z_Px_ab+ABX*I_NAI_Fx2z_Px_ab;
  Double I_NAI_F3y_D2x_ab = I_NAI_Gx3y_Px_ab+ABX*I_NAI_F3y_Px_ab;
  Double I_NAI_F2yz_D2x_ab = I_NAI_Gx2yz_Px_ab+ABX*I_NAI_F2yz_Px_ab;
  Double I_NAI_Fy2z_D2x_ab = I_NAI_Gxy2z_Px_ab+ABX*I_NAI_Fy2z_Px_ab;
  Double I_NAI_F3z_D2x_ab = I_NAI_Gx3z_Px_ab+ABX*I_NAI_F3z_Px_ab;
  Double I_NAI_F3x_Dxy_ab = I_NAI_G3xy_Px_ab+ABY*I_NAI_F3x_Px_ab;
  Double I_NAI_F2xy_Dxy_ab = I_NAI_G2x2y_Px_ab+ABY*I_NAI_F2xy_Px_ab;
  Double I_NAI_F2xz_Dxy_ab = I_NAI_G2xyz_Px_ab+ABY*I_NAI_F2xz_Px_ab;
  Double I_NAI_Fx2y_Dxy_ab = I_NAI_Gx3y_Px_ab+ABY*I_NAI_Fx2y_Px_ab;
  Double I_NAI_Fxyz_Dxy_ab = I_NAI_Gx2yz_Px_ab+ABY*I_NAI_Fxyz_Px_ab;
  Double I_NAI_Fx2z_Dxy_ab = I_NAI_Gxy2z_Px_ab+ABY*I_NAI_Fx2z_Px_ab;
  Double I_NAI_F3y_Dxy_ab = I_NAI_G4y_Px_ab+ABY*I_NAI_F3y_Px_ab;
  Double I_NAI_F2yz_Dxy_ab = I_NAI_G3yz_Px_ab+ABY*I_NAI_F2yz_Px_ab;
  Double I_NAI_Fy2z_Dxy_ab = I_NAI_G2y2z_Px_ab+ABY*I_NAI_Fy2z_Px_ab;
  Double I_NAI_F3z_Dxy_ab = I_NAI_Gy3z_Px_ab+ABY*I_NAI_F3z_Px_ab;
  Double I_NAI_F3x_D2y_ab = I_NAI_G3xy_Py_ab+ABY*I_NAI_F3x_Py_ab;
  Double I_NAI_F2xy_D2y_ab = I_NAI_G2x2y_Py_ab+ABY*I_NAI_F2xy_Py_ab;
  Double I_NAI_F2xz_D2y_ab = I_NAI_G2xyz_Py_ab+ABY*I_NAI_F2xz_Py_ab;
  Double I_NAI_Fx2y_D2y_ab = I_NAI_Gx3y_Py_ab+ABY*I_NAI_Fx2y_Py_ab;
  Double I_NAI_Fxyz_D2y_ab = I_NAI_Gx2yz_Py_ab+ABY*I_NAI_Fxyz_Py_ab;
  Double I_NAI_Fx2z_D2y_ab = I_NAI_Gxy2z_Py_ab+ABY*I_NAI_Fx2z_Py_ab;
  Double I_NAI_F3y_D2y_ab = I_NAI_G4y_Py_ab+ABY*I_NAI_F3y_Py_ab;
  Double I_NAI_F2yz_D2y_ab = I_NAI_G3yz_Py_ab+ABY*I_NAI_F2yz_Py_ab;
  Double I_NAI_Fy2z_D2y_ab = I_NAI_G2y2z_Py_ab+ABY*I_NAI_Fy2z_Py_ab;
  Double I_NAI_F3z_D2y_ab = I_NAI_Gy3z_Py_ab+ABY*I_NAI_F3z_Py_ab;
  Double I_NAI_F3x_D2z_ab = I_NAI_G3xz_Pz_ab+ABZ*I_NAI_F3x_Pz_ab;
  Double I_NAI_F2xy_D2z_ab = I_NAI_G2xyz_Pz_ab+ABZ*I_NAI_F2xy_Pz_ab;
  Double I_NAI_F2xz_D2z_ab = I_NAI_G2x2z_Pz_ab+ABZ*I_NAI_F2xz_Pz_ab;
  Double I_NAI_Fx2y_D2z_ab = I_NAI_Gx2yz_Pz_ab+ABZ*I_NAI_Fx2y_Pz_ab;
  Double I_NAI_Fxyz_D2z_ab = I_NAI_Gxy2z_Pz_ab+ABZ*I_NAI_Fxyz_Pz_ab;
  Double I_NAI_Fx2z_D2z_ab = I_NAI_Gx3z_Pz_ab+ABZ*I_NAI_Fx2z_Pz_ab;
  Double I_NAI_F3y_D2z_ab = I_NAI_G3yz_Pz_ab+ABZ*I_NAI_F3y_Pz_ab;
  Double I_NAI_F2yz_D2z_ab = I_NAI_G2y2z_Pz_ab+ABZ*I_NAI_F2yz_Pz_ab;
  Double I_NAI_Fy2z_D2z_ab = I_NAI_Gy3z_Pz_ab+ABZ*I_NAI_Fy2z_Pz_ab;
  Double I_NAI_F3z_D2z_ab = I_NAI_G4z_Pz_ab+ABZ*I_NAI_F3z_Pz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 14 integrals are omitted 
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
  Double I_NAI_H4yz_Px_ab = I_NAI_Ix4yz_S_ab+ABX*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Px_ab = I_NAI_Ix3y2z_S_ab+ABX*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Px_ab = I_NAI_Ix2y3z_S_ab+ABX*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Px_ab = I_NAI_Ixy4z_S_ab+ABX*I_NAI_Hy4z_S_ab;
  Double I_NAI_H4xy_Py_ab = I_NAI_I4x2y_S_ab+ABY*I_NAI_H4xy_S_ab;
  Double I_NAI_H3x2y_Py_ab = I_NAI_I3x3y_S_ab+ABY*I_NAI_H3x2y_S_ab;
  Double I_NAI_H3xyz_Py_ab = I_NAI_I3x2yz_S_ab+ABY*I_NAI_H3xyz_S_ab;
  Double I_NAI_H2x3y_Py_ab = I_NAI_I2x4y_S_ab+ABY*I_NAI_H2x3y_S_ab;
  Double I_NAI_H2x2yz_Py_ab = I_NAI_I2x3yz_S_ab+ABY*I_NAI_H2x2yz_S_ab;
  Double I_NAI_H2xy2z_Py_ab = I_NAI_I2x2y2z_S_ab+ABY*I_NAI_H2xy2z_S_ab;
  Double I_NAI_Hx4y_Py_ab = I_NAI_Ix5y_S_ab+ABY*I_NAI_Hx4y_S_ab;
  Double I_NAI_Hx3yz_Py_ab = I_NAI_Ix4yz_S_ab+ABY*I_NAI_Hx3yz_S_ab;
  Double I_NAI_Hx2y2z_Py_ab = I_NAI_Ix3y2z_S_ab+ABY*I_NAI_Hx2y2z_S_ab;
  Double I_NAI_Hxy3z_Py_ab = I_NAI_Ix2y3z_S_ab+ABY*I_NAI_Hxy3z_S_ab;
  Double I_NAI_H5y_Py_ab = I_NAI_I6y_S_ab+ABY*I_NAI_H5y_S_ab;
  Double I_NAI_H4yz_Py_ab = I_NAI_I5yz_S_ab+ABY*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Py_ab = I_NAI_I4y2z_S_ab+ABY*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Py_ab = I_NAI_I3y3z_S_ab+ABY*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Py_ab = I_NAI_I2y4z_S_ab+ABY*I_NAI_Hy4z_S_ab;
  Double I_NAI_H4xz_Pz_ab = I_NAI_I4x2z_S_ab+ABZ*I_NAI_H4xz_S_ab;
  Double I_NAI_H3xyz_Pz_ab = I_NAI_I3xy2z_S_ab+ABZ*I_NAI_H3xyz_S_ab;
  Double I_NAI_H3x2z_Pz_ab = I_NAI_I3x3z_S_ab+ABZ*I_NAI_H3x2z_S_ab;
  Double I_NAI_H2x2yz_Pz_ab = I_NAI_I2x2y2z_S_ab+ABZ*I_NAI_H2x2yz_S_ab;
  Double I_NAI_H2xy2z_Pz_ab = I_NAI_I2xy3z_S_ab+ABZ*I_NAI_H2xy2z_S_ab;
  Double I_NAI_H2x3z_Pz_ab = I_NAI_I2x4z_S_ab+ABZ*I_NAI_H2x3z_S_ab;
  Double I_NAI_Hx3yz_Pz_ab = I_NAI_Ix3y2z_S_ab+ABZ*I_NAI_Hx3yz_S_ab;
  Double I_NAI_Hx2y2z_Pz_ab = I_NAI_Ix2y3z_S_ab+ABZ*I_NAI_Hx2y2z_S_ab;
  Double I_NAI_Hxy3z_Pz_ab = I_NAI_Ixy4z_S_ab+ABZ*I_NAI_Hxy3z_S_ab;
  Double I_NAI_Hx4z_Pz_ab = I_NAI_Ix5z_S_ab+ABZ*I_NAI_Hx4z_S_ab;
  Double I_NAI_H4yz_Pz_ab = I_NAI_I4y2z_S_ab+ABZ*I_NAI_H4yz_S_ab;
  Double I_NAI_H3y2z_Pz_ab = I_NAI_I3y3z_S_ab+ABZ*I_NAI_H3y2z_S_ab;
  Double I_NAI_H2y3z_Pz_ab = I_NAI_I2y4z_S_ab+ABZ*I_NAI_H2y3z_S_ab;
  Double I_NAI_Hy4z_Pz_ab = I_NAI_Iy5z_S_ab+ABZ*I_NAI_Hy4z_S_ab;
  Double I_NAI_H5z_Pz_ab = I_NAI_I6z_S_ab+ABZ*I_NAI_H5z_S_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 35 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_ab
   * RHS shell quartet name: SQ_NAI_G_P_ab
   ************************************************************/
  Double I_NAI_G4x_D2x_ab = I_NAI_H5x_Px_ab+ABX*I_NAI_G4x_Px_ab;
  Double I_NAI_G3xy_D2x_ab = I_NAI_H4xy_Px_ab+ABX*I_NAI_G3xy_Px_ab;
  Double I_NAI_G3xz_D2x_ab = I_NAI_H4xz_Px_ab+ABX*I_NAI_G3xz_Px_ab;
  Double I_NAI_G2x2y_D2x_ab = I_NAI_H3x2y_Px_ab+ABX*I_NAI_G2x2y_Px_ab;
  Double I_NAI_G2xyz_D2x_ab = I_NAI_H3xyz_Px_ab+ABX*I_NAI_G2xyz_Px_ab;
  Double I_NAI_G2x2z_D2x_ab = I_NAI_H3x2z_Px_ab+ABX*I_NAI_G2x2z_Px_ab;
  Double I_NAI_Gx3y_D2x_ab = I_NAI_H2x3y_Px_ab+ABX*I_NAI_Gx3y_Px_ab;
  Double I_NAI_Gx2yz_D2x_ab = I_NAI_H2x2yz_Px_ab+ABX*I_NAI_Gx2yz_Px_ab;
  Double I_NAI_Gxy2z_D2x_ab = I_NAI_H2xy2z_Px_ab+ABX*I_NAI_Gxy2z_Px_ab;
  Double I_NAI_Gx3z_D2x_ab = I_NAI_H2x3z_Px_ab+ABX*I_NAI_Gx3z_Px_ab;
  Double I_NAI_G4y_D2x_ab = I_NAI_Hx4y_Px_ab+ABX*I_NAI_G4y_Px_ab;
  Double I_NAI_G3yz_D2x_ab = I_NAI_Hx3yz_Px_ab+ABX*I_NAI_G3yz_Px_ab;
  Double I_NAI_G2y2z_D2x_ab = I_NAI_Hx2y2z_Px_ab+ABX*I_NAI_G2y2z_Px_ab;
  Double I_NAI_Gy3z_D2x_ab = I_NAI_Hxy3z_Px_ab+ABX*I_NAI_Gy3z_Px_ab;
  Double I_NAI_G4z_D2x_ab = I_NAI_Hx4z_Px_ab+ABX*I_NAI_G4z_Px_ab;
  Double I_NAI_G3xz_Dxy_ab = I_NAI_H3xyz_Px_ab+ABY*I_NAI_G3xz_Px_ab;
  Double I_NAI_G2xyz_Dxy_ab = I_NAI_H2x2yz_Px_ab+ABY*I_NAI_G2xyz_Px_ab;
  Double I_NAI_G2x2z_Dxy_ab = I_NAI_H2xy2z_Px_ab+ABY*I_NAI_G2x2z_Px_ab;
  Double I_NAI_Gx2yz_Dxy_ab = I_NAI_Hx3yz_Px_ab+ABY*I_NAI_Gx2yz_Px_ab;
  Double I_NAI_Gxy2z_Dxy_ab = I_NAI_Hx2y2z_Px_ab+ABY*I_NAI_Gxy2z_Px_ab;
  Double I_NAI_Gx3z_Dxy_ab = I_NAI_Hxy3z_Px_ab+ABY*I_NAI_Gx3z_Px_ab;
  Double I_NAI_G3yz_Dxy_ab = I_NAI_H4yz_Px_ab+ABY*I_NAI_G3yz_Px_ab;
  Double I_NAI_G2y2z_Dxy_ab = I_NAI_H3y2z_Px_ab+ABY*I_NAI_G2y2z_Px_ab;
  Double I_NAI_Gy3z_Dxy_ab = I_NAI_H2y3z_Px_ab+ABY*I_NAI_Gy3z_Px_ab;
  Double I_NAI_G4z_Dxy_ab = I_NAI_Hy4z_Px_ab+ABY*I_NAI_G4z_Px_ab;
  Double I_NAI_G4x_D2y_ab = I_NAI_H4xy_Py_ab+ABY*I_NAI_G4x_Py_ab;
  Double I_NAI_G3xy_D2y_ab = I_NAI_H3x2y_Py_ab+ABY*I_NAI_G3xy_Py_ab;
  Double I_NAI_G3xz_D2y_ab = I_NAI_H3xyz_Py_ab+ABY*I_NAI_G3xz_Py_ab;
  Double I_NAI_G2x2y_D2y_ab = I_NAI_H2x3y_Py_ab+ABY*I_NAI_G2x2y_Py_ab;
  Double I_NAI_G2xyz_D2y_ab = I_NAI_H2x2yz_Py_ab+ABY*I_NAI_G2xyz_Py_ab;
  Double I_NAI_G2x2z_D2y_ab = I_NAI_H2xy2z_Py_ab+ABY*I_NAI_G2x2z_Py_ab;
  Double I_NAI_Gx3y_D2y_ab = I_NAI_Hx4y_Py_ab+ABY*I_NAI_Gx3y_Py_ab;
  Double I_NAI_Gx2yz_D2y_ab = I_NAI_Hx3yz_Py_ab+ABY*I_NAI_Gx2yz_Py_ab;
  Double I_NAI_Gxy2z_D2y_ab = I_NAI_Hx2y2z_Py_ab+ABY*I_NAI_Gxy2z_Py_ab;
  Double I_NAI_Gx3z_D2y_ab = I_NAI_Hxy3z_Py_ab+ABY*I_NAI_Gx3z_Py_ab;
  Double I_NAI_G4y_D2y_ab = I_NAI_H5y_Py_ab+ABY*I_NAI_G4y_Py_ab;
  Double I_NAI_G3yz_D2y_ab = I_NAI_H4yz_Py_ab+ABY*I_NAI_G3yz_Py_ab;
  Double I_NAI_G2y2z_D2y_ab = I_NAI_H3y2z_Py_ab+ABY*I_NAI_G2y2z_Py_ab;
  Double I_NAI_Gy3z_D2y_ab = I_NAI_H2y3z_Py_ab+ABY*I_NAI_Gy3z_Py_ab;
  Double I_NAI_G4z_D2y_ab = I_NAI_Hy4z_Py_ab+ABY*I_NAI_G4z_Py_ab;
  Double I_NAI_G4x_D2z_ab = I_NAI_H4xz_Pz_ab+ABZ*I_NAI_G4x_Pz_ab;
  Double I_NAI_G3xy_D2z_ab = I_NAI_H3xyz_Pz_ab+ABZ*I_NAI_G3xy_Pz_ab;
  Double I_NAI_G3xz_D2z_ab = I_NAI_H3x2z_Pz_ab+ABZ*I_NAI_G3xz_Pz_ab;
  Double I_NAI_G2x2y_D2z_ab = I_NAI_H2x2yz_Pz_ab+ABZ*I_NAI_G2x2y_Pz_ab;
  Double I_NAI_G2xyz_D2z_ab = I_NAI_H2xy2z_Pz_ab+ABZ*I_NAI_G2xyz_Pz_ab;
  Double I_NAI_G2x2z_D2z_ab = I_NAI_H2x3z_Pz_ab+ABZ*I_NAI_G2x2z_Pz_ab;
  Double I_NAI_Gx3y_D2z_ab = I_NAI_Hx3yz_Pz_ab+ABZ*I_NAI_Gx3y_Pz_ab;
  Double I_NAI_Gx2yz_D2z_ab = I_NAI_Hx2y2z_Pz_ab+ABZ*I_NAI_Gx2yz_Pz_ab;
  Double I_NAI_Gxy2z_D2z_ab = I_NAI_Hxy3z_Pz_ab+ABZ*I_NAI_Gxy2z_Pz_ab;
  Double I_NAI_Gx3z_D2z_ab = I_NAI_Hx4z_Pz_ab+ABZ*I_NAI_Gx3z_Pz_ab;
  Double I_NAI_G4y_D2z_ab = I_NAI_H4yz_Pz_ab+ABZ*I_NAI_G4y_Pz_ab;
  Double I_NAI_G3yz_D2z_ab = I_NAI_H3y2z_Pz_ab+ABZ*I_NAI_G3yz_Pz_ab;
  Double I_NAI_G2y2z_D2z_ab = I_NAI_H2y3z_Pz_ab+ABZ*I_NAI_G2y2z_Pz_ab;
  Double I_NAI_Gy3z_D2z_ab = I_NAI_Hy4z_Pz_ab+ABZ*I_NAI_Gy3z_Pz_ab;
  Double I_NAI_G4z_D2z_ab = I_NAI_H5z_Pz_ab+ABZ*I_NAI_G4z_Pz_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_ab
   * RHS shell quartet name: SQ_NAI_F_D_ab
   ************************************************************/
  Double I_NAI_F3x_F3x_ab = I_NAI_G4x_D2x_ab+ABX*I_NAI_F3x_D2x_ab;
  Double I_NAI_F2xy_F3x_ab = I_NAI_G3xy_D2x_ab+ABX*I_NAI_F2xy_D2x_ab;
  Double I_NAI_F2xz_F3x_ab = I_NAI_G3xz_D2x_ab+ABX*I_NAI_F2xz_D2x_ab;
  Double I_NAI_Fx2y_F3x_ab = I_NAI_G2x2y_D2x_ab+ABX*I_NAI_Fx2y_D2x_ab;
  Double I_NAI_Fxyz_F3x_ab = I_NAI_G2xyz_D2x_ab+ABX*I_NAI_Fxyz_D2x_ab;
  Double I_NAI_Fx2z_F3x_ab = I_NAI_G2x2z_D2x_ab+ABX*I_NAI_Fx2z_D2x_ab;
  Double I_NAI_F3y_F3x_ab = I_NAI_Gx3y_D2x_ab+ABX*I_NAI_F3y_D2x_ab;
  Double I_NAI_F2yz_F3x_ab = I_NAI_Gx2yz_D2x_ab+ABX*I_NAI_F2yz_D2x_ab;
  Double I_NAI_Fy2z_F3x_ab = I_NAI_Gxy2z_D2x_ab+ABX*I_NAI_Fy2z_D2x_ab;
  Double I_NAI_F3z_F3x_ab = I_NAI_Gx3z_D2x_ab+ABX*I_NAI_F3z_D2x_ab;
  Double I_NAI_F3x_F2xy_ab = I_NAI_G3xy_D2x_ab+ABY*I_NAI_F3x_D2x_ab;
  Double I_NAI_F2xy_F2xy_ab = I_NAI_G2x2y_D2x_ab+ABY*I_NAI_F2xy_D2x_ab;
  Double I_NAI_F2xz_F2xy_ab = I_NAI_G2xyz_D2x_ab+ABY*I_NAI_F2xz_D2x_ab;
  Double I_NAI_Fx2y_F2xy_ab = I_NAI_Gx3y_D2x_ab+ABY*I_NAI_Fx2y_D2x_ab;
  Double I_NAI_Fxyz_F2xy_ab = I_NAI_Gx2yz_D2x_ab+ABY*I_NAI_Fxyz_D2x_ab;
  Double I_NAI_Fx2z_F2xy_ab = I_NAI_Gxy2z_D2x_ab+ABY*I_NAI_Fx2z_D2x_ab;
  Double I_NAI_F3y_F2xy_ab = I_NAI_G4y_D2x_ab+ABY*I_NAI_F3y_D2x_ab;
  Double I_NAI_F2yz_F2xy_ab = I_NAI_G3yz_D2x_ab+ABY*I_NAI_F2yz_D2x_ab;
  Double I_NAI_Fy2z_F2xy_ab = I_NAI_G2y2z_D2x_ab+ABY*I_NAI_Fy2z_D2x_ab;
  Double I_NAI_F3z_F2xy_ab = I_NAI_Gy3z_D2x_ab+ABY*I_NAI_F3z_D2x_ab;
  Double I_NAI_F3x_F2xz_ab = I_NAI_G3xz_D2x_ab+ABZ*I_NAI_F3x_D2x_ab;
  Double I_NAI_F2xy_F2xz_ab = I_NAI_G2xyz_D2x_ab+ABZ*I_NAI_F2xy_D2x_ab;
  Double I_NAI_F2xz_F2xz_ab = I_NAI_G2x2z_D2x_ab+ABZ*I_NAI_F2xz_D2x_ab;
  Double I_NAI_Fx2y_F2xz_ab = I_NAI_Gx2yz_D2x_ab+ABZ*I_NAI_Fx2y_D2x_ab;
  Double I_NAI_Fxyz_F2xz_ab = I_NAI_Gxy2z_D2x_ab+ABZ*I_NAI_Fxyz_D2x_ab;
  Double I_NAI_Fx2z_F2xz_ab = I_NAI_Gx3z_D2x_ab+ABZ*I_NAI_Fx2z_D2x_ab;
  Double I_NAI_F3y_F2xz_ab = I_NAI_G3yz_D2x_ab+ABZ*I_NAI_F3y_D2x_ab;
  Double I_NAI_F2yz_F2xz_ab = I_NAI_G2y2z_D2x_ab+ABZ*I_NAI_F2yz_D2x_ab;
  Double I_NAI_Fy2z_F2xz_ab = I_NAI_Gy3z_D2x_ab+ABZ*I_NAI_Fy2z_D2x_ab;
  Double I_NAI_F3z_F2xz_ab = I_NAI_G4z_D2x_ab+ABZ*I_NAI_F3z_D2x_ab;
  Double I_NAI_F3x_Fx2y_ab = I_NAI_G4x_D2y_ab+ABX*I_NAI_F3x_D2y_ab;
  Double I_NAI_F2xy_Fx2y_ab = I_NAI_G3xy_D2y_ab+ABX*I_NAI_F2xy_D2y_ab;
  Double I_NAI_F2xz_Fx2y_ab = I_NAI_G3xz_D2y_ab+ABX*I_NAI_F2xz_D2y_ab;
  Double I_NAI_Fx2y_Fx2y_ab = I_NAI_G2x2y_D2y_ab+ABX*I_NAI_Fx2y_D2y_ab;
  Double I_NAI_Fxyz_Fx2y_ab = I_NAI_G2xyz_D2y_ab+ABX*I_NAI_Fxyz_D2y_ab;
  Double I_NAI_Fx2z_Fx2y_ab = I_NAI_G2x2z_D2y_ab+ABX*I_NAI_Fx2z_D2y_ab;
  Double I_NAI_F3y_Fx2y_ab = I_NAI_Gx3y_D2y_ab+ABX*I_NAI_F3y_D2y_ab;
  Double I_NAI_F2yz_Fx2y_ab = I_NAI_Gx2yz_D2y_ab+ABX*I_NAI_F2yz_D2y_ab;
  Double I_NAI_Fy2z_Fx2y_ab = I_NAI_Gxy2z_D2y_ab+ABX*I_NAI_Fy2z_D2y_ab;
  Double I_NAI_F3z_Fx2y_ab = I_NAI_Gx3z_D2y_ab+ABX*I_NAI_F3z_D2y_ab;
  Double I_NAI_F3x_Fxyz_ab = I_NAI_G3xz_Dxy_ab+ABZ*I_NAI_F3x_Dxy_ab;
  Double I_NAI_F2xy_Fxyz_ab = I_NAI_G2xyz_Dxy_ab+ABZ*I_NAI_F2xy_Dxy_ab;
  Double I_NAI_F2xz_Fxyz_ab = I_NAI_G2x2z_Dxy_ab+ABZ*I_NAI_F2xz_Dxy_ab;
  Double I_NAI_Fx2y_Fxyz_ab = I_NAI_Gx2yz_Dxy_ab+ABZ*I_NAI_Fx2y_Dxy_ab;
  Double I_NAI_Fxyz_Fxyz_ab = I_NAI_Gxy2z_Dxy_ab+ABZ*I_NAI_Fxyz_Dxy_ab;
  Double I_NAI_Fx2z_Fxyz_ab = I_NAI_Gx3z_Dxy_ab+ABZ*I_NAI_Fx2z_Dxy_ab;
  Double I_NAI_F3y_Fxyz_ab = I_NAI_G3yz_Dxy_ab+ABZ*I_NAI_F3y_Dxy_ab;
  Double I_NAI_F2yz_Fxyz_ab = I_NAI_G2y2z_Dxy_ab+ABZ*I_NAI_F2yz_Dxy_ab;
  Double I_NAI_Fy2z_Fxyz_ab = I_NAI_Gy3z_Dxy_ab+ABZ*I_NAI_Fy2z_Dxy_ab;
  Double I_NAI_F3z_Fxyz_ab = I_NAI_G4z_Dxy_ab+ABZ*I_NAI_F3z_Dxy_ab;
  Double I_NAI_F3x_Fx2z_ab = I_NAI_G4x_D2z_ab+ABX*I_NAI_F3x_D2z_ab;
  Double I_NAI_F2xy_Fx2z_ab = I_NAI_G3xy_D2z_ab+ABX*I_NAI_F2xy_D2z_ab;
  Double I_NAI_F2xz_Fx2z_ab = I_NAI_G3xz_D2z_ab+ABX*I_NAI_F2xz_D2z_ab;
  Double I_NAI_Fx2y_Fx2z_ab = I_NAI_G2x2y_D2z_ab+ABX*I_NAI_Fx2y_D2z_ab;
  Double I_NAI_Fxyz_Fx2z_ab = I_NAI_G2xyz_D2z_ab+ABX*I_NAI_Fxyz_D2z_ab;
  Double I_NAI_Fx2z_Fx2z_ab = I_NAI_G2x2z_D2z_ab+ABX*I_NAI_Fx2z_D2z_ab;
  Double I_NAI_F3y_Fx2z_ab = I_NAI_Gx3y_D2z_ab+ABX*I_NAI_F3y_D2z_ab;
  Double I_NAI_F2yz_Fx2z_ab = I_NAI_Gx2yz_D2z_ab+ABX*I_NAI_F2yz_D2z_ab;
  Double I_NAI_Fy2z_Fx2z_ab = I_NAI_Gxy2z_D2z_ab+ABX*I_NAI_Fy2z_D2z_ab;
  Double I_NAI_F3z_Fx2z_ab = I_NAI_Gx3z_D2z_ab+ABX*I_NAI_F3z_D2z_ab;
  Double I_NAI_F3x_F3y_ab = I_NAI_G3xy_D2y_ab+ABY*I_NAI_F3x_D2y_ab;
  Double I_NAI_F2xy_F3y_ab = I_NAI_G2x2y_D2y_ab+ABY*I_NAI_F2xy_D2y_ab;
  Double I_NAI_F2xz_F3y_ab = I_NAI_G2xyz_D2y_ab+ABY*I_NAI_F2xz_D2y_ab;
  Double I_NAI_Fx2y_F3y_ab = I_NAI_Gx3y_D2y_ab+ABY*I_NAI_Fx2y_D2y_ab;
  Double I_NAI_Fxyz_F3y_ab = I_NAI_Gx2yz_D2y_ab+ABY*I_NAI_Fxyz_D2y_ab;
  Double I_NAI_Fx2z_F3y_ab = I_NAI_Gxy2z_D2y_ab+ABY*I_NAI_Fx2z_D2y_ab;
  Double I_NAI_F3y_F3y_ab = I_NAI_G4y_D2y_ab+ABY*I_NAI_F3y_D2y_ab;
  Double I_NAI_F2yz_F3y_ab = I_NAI_G3yz_D2y_ab+ABY*I_NAI_F2yz_D2y_ab;
  Double I_NAI_Fy2z_F3y_ab = I_NAI_G2y2z_D2y_ab+ABY*I_NAI_Fy2z_D2y_ab;
  Double I_NAI_F3z_F3y_ab = I_NAI_Gy3z_D2y_ab+ABY*I_NAI_F3z_D2y_ab;
  Double I_NAI_F3x_F2yz_ab = I_NAI_G3xz_D2y_ab+ABZ*I_NAI_F3x_D2y_ab;
  Double I_NAI_F2xy_F2yz_ab = I_NAI_G2xyz_D2y_ab+ABZ*I_NAI_F2xy_D2y_ab;
  Double I_NAI_F2xz_F2yz_ab = I_NAI_G2x2z_D2y_ab+ABZ*I_NAI_F2xz_D2y_ab;
  Double I_NAI_Fx2y_F2yz_ab = I_NAI_Gx2yz_D2y_ab+ABZ*I_NAI_Fx2y_D2y_ab;
  Double I_NAI_Fxyz_F2yz_ab = I_NAI_Gxy2z_D2y_ab+ABZ*I_NAI_Fxyz_D2y_ab;
  Double I_NAI_Fx2z_F2yz_ab = I_NAI_Gx3z_D2y_ab+ABZ*I_NAI_Fx2z_D2y_ab;
  Double I_NAI_F3y_F2yz_ab = I_NAI_G3yz_D2y_ab+ABZ*I_NAI_F3y_D2y_ab;
  Double I_NAI_F2yz_F2yz_ab = I_NAI_G2y2z_D2y_ab+ABZ*I_NAI_F2yz_D2y_ab;
  Double I_NAI_Fy2z_F2yz_ab = I_NAI_Gy3z_D2y_ab+ABZ*I_NAI_Fy2z_D2y_ab;
  Double I_NAI_F3z_F2yz_ab = I_NAI_G4z_D2y_ab+ABZ*I_NAI_F3z_D2y_ab;
  Double I_NAI_F3x_Fy2z_ab = I_NAI_G3xy_D2z_ab+ABY*I_NAI_F3x_D2z_ab;
  Double I_NAI_F2xy_Fy2z_ab = I_NAI_G2x2y_D2z_ab+ABY*I_NAI_F2xy_D2z_ab;
  Double I_NAI_F2xz_Fy2z_ab = I_NAI_G2xyz_D2z_ab+ABY*I_NAI_F2xz_D2z_ab;
  Double I_NAI_Fx2y_Fy2z_ab = I_NAI_Gx3y_D2z_ab+ABY*I_NAI_Fx2y_D2z_ab;
  Double I_NAI_Fxyz_Fy2z_ab = I_NAI_Gx2yz_D2z_ab+ABY*I_NAI_Fxyz_D2z_ab;
  Double I_NAI_Fx2z_Fy2z_ab = I_NAI_Gxy2z_D2z_ab+ABY*I_NAI_Fx2z_D2z_ab;
  Double I_NAI_F3y_Fy2z_ab = I_NAI_G4y_D2z_ab+ABY*I_NAI_F3y_D2z_ab;
  Double I_NAI_F2yz_Fy2z_ab = I_NAI_G3yz_D2z_ab+ABY*I_NAI_F2yz_D2z_ab;
  Double I_NAI_Fy2z_Fy2z_ab = I_NAI_G2y2z_D2z_ab+ABY*I_NAI_Fy2z_D2z_ab;
  Double I_NAI_F3z_Fy2z_ab = I_NAI_Gy3z_D2z_ab+ABY*I_NAI_F3z_D2z_ab;
  Double I_NAI_F3x_F3z_ab = I_NAI_G3xz_D2z_ab+ABZ*I_NAI_F3x_D2z_ab;
  Double I_NAI_F2xy_F3z_ab = I_NAI_G2xyz_D2z_ab+ABZ*I_NAI_F2xy_D2z_ab;
  Double I_NAI_F2xz_F3z_ab = I_NAI_G2x2z_D2z_ab+ABZ*I_NAI_F2xz_D2z_ab;
  Double I_NAI_Fx2y_F3z_ab = I_NAI_Gx2yz_D2z_ab+ABZ*I_NAI_Fx2y_D2z_ab;
  Double I_NAI_Fxyz_F3z_ab = I_NAI_Gxy2z_D2z_ab+ABZ*I_NAI_Fxyz_D2z_ab;
  Double I_NAI_Fx2z_F3z_ab = I_NAI_Gx3z_D2z_ab+ABZ*I_NAI_Fx2z_D2z_ab;
  Double I_NAI_F3y_F3z_ab = I_NAI_G3yz_D2z_ab+ABZ*I_NAI_F3y_D2z_ab;
  Double I_NAI_F2yz_F3z_ab = I_NAI_G2y2z_D2z_ab+ABZ*I_NAI_F2yz_D2z_ab;
  Double I_NAI_Fy2z_F3z_ab = I_NAI_Gy3z_D2z_ab+ABZ*I_NAI_Fy2z_D2z_ab;
  Double I_NAI_F3z_F3z_ab = I_NAI_G4z_D2z_ab+ABZ*I_NAI_F3z_D2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_bb
   * RHS shell quartet name: SQ_NAI_D_S_bb
   ************************************************************/
  Double I_NAI_D2x_Px_bb = I_NAI_F3x_S_bb+ABX*I_NAI_D2x_S_bb;
  Double I_NAI_Dxy_Px_bb = I_NAI_F2xy_S_bb+ABX*I_NAI_Dxy_S_bb;
  Double I_NAI_Dxz_Px_bb = I_NAI_F2xz_S_bb+ABX*I_NAI_Dxz_S_bb;
  Double I_NAI_D2y_Px_bb = I_NAI_Fx2y_S_bb+ABX*I_NAI_D2y_S_bb;
  Double I_NAI_Dyz_Px_bb = I_NAI_Fxyz_S_bb+ABX*I_NAI_Dyz_S_bb;
  Double I_NAI_D2z_Px_bb = I_NAI_Fx2z_S_bb+ABX*I_NAI_D2z_S_bb;
  Double I_NAI_D2x_Py_bb = I_NAI_F2xy_S_bb+ABY*I_NAI_D2x_S_bb;
  Double I_NAI_Dxy_Py_bb = I_NAI_Fx2y_S_bb+ABY*I_NAI_Dxy_S_bb;
  Double I_NAI_Dxz_Py_bb = I_NAI_Fxyz_S_bb+ABY*I_NAI_Dxz_S_bb;
  Double I_NAI_D2y_Py_bb = I_NAI_F3y_S_bb+ABY*I_NAI_D2y_S_bb;
  Double I_NAI_Dyz_Py_bb = I_NAI_F2yz_S_bb+ABY*I_NAI_Dyz_S_bb;
  Double I_NAI_D2z_Py_bb = I_NAI_Fy2z_S_bb+ABY*I_NAI_D2z_S_bb;
  Double I_NAI_D2x_Pz_bb = I_NAI_F2xz_S_bb+ABZ*I_NAI_D2x_S_bb;
  Double I_NAI_Dxy_Pz_bb = I_NAI_Fxyz_S_bb+ABZ*I_NAI_Dxy_S_bb;
  Double I_NAI_Dxz_Pz_bb = I_NAI_Fx2z_S_bb+ABZ*I_NAI_Dxz_S_bb;
  Double I_NAI_D2y_Pz_bb = I_NAI_F2yz_S_bb+ABZ*I_NAI_D2y_S_bb;
  Double I_NAI_Dyz_Pz_bb = I_NAI_Fy2z_S_bb+ABZ*I_NAI_Dyz_S_bb;
  Double I_NAI_D2z_Pz_bb = I_NAI_F3z_S_bb+ABZ*I_NAI_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_bb
   * RHS shell quartet name: SQ_NAI_F_S_bb
   ************************************************************/
  Double I_NAI_F3x_Px_bb = I_NAI_G4x_S_bb+ABX*I_NAI_F3x_S_bb;
  Double I_NAI_F2xy_Px_bb = I_NAI_G3xy_S_bb+ABX*I_NAI_F2xy_S_bb;
  Double I_NAI_F2xz_Px_bb = I_NAI_G3xz_S_bb+ABX*I_NAI_F2xz_S_bb;
  Double I_NAI_Fx2y_Px_bb = I_NAI_G2x2y_S_bb+ABX*I_NAI_Fx2y_S_bb;
  Double I_NAI_Fxyz_Px_bb = I_NAI_G2xyz_S_bb+ABX*I_NAI_Fxyz_S_bb;
  Double I_NAI_Fx2z_Px_bb = I_NAI_G2x2z_S_bb+ABX*I_NAI_Fx2z_S_bb;
  Double I_NAI_F3y_Px_bb = I_NAI_Gx3y_S_bb+ABX*I_NAI_F3y_S_bb;
  Double I_NAI_F2yz_Px_bb = I_NAI_Gx2yz_S_bb+ABX*I_NAI_F2yz_S_bb;
  Double I_NAI_Fy2z_Px_bb = I_NAI_Gxy2z_S_bb+ABX*I_NAI_Fy2z_S_bb;
  Double I_NAI_F3z_Px_bb = I_NAI_Gx3z_S_bb+ABX*I_NAI_F3z_S_bb;
  Double I_NAI_F3x_Py_bb = I_NAI_G3xy_S_bb+ABY*I_NAI_F3x_S_bb;
  Double I_NAI_F2xy_Py_bb = I_NAI_G2x2y_S_bb+ABY*I_NAI_F2xy_S_bb;
  Double I_NAI_F2xz_Py_bb = I_NAI_G2xyz_S_bb+ABY*I_NAI_F2xz_S_bb;
  Double I_NAI_Fx2y_Py_bb = I_NAI_Gx3y_S_bb+ABY*I_NAI_Fx2y_S_bb;
  Double I_NAI_Fxyz_Py_bb = I_NAI_Gx2yz_S_bb+ABY*I_NAI_Fxyz_S_bb;
  Double I_NAI_Fx2z_Py_bb = I_NAI_Gxy2z_S_bb+ABY*I_NAI_Fx2z_S_bb;
  Double I_NAI_F3y_Py_bb = I_NAI_G4y_S_bb+ABY*I_NAI_F3y_S_bb;
  Double I_NAI_F2yz_Py_bb = I_NAI_G3yz_S_bb+ABY*I_NAI_F2yz_S_bb;
  Double I_NAI_Fy2z_Py_bb = I_NAI_G2y2z_S_bb+ABY*I_NAI_Fy2z_S_bb;
  Double I_NAI_F3z_Py_bb = I_NAI_Gy3z_S_bb+ABY*I_NAI_F3z_S_bb;
  Double I_NAI_F3x_Pz_bb = I_NAI_G3xz_S_bb+ABZ*I_NAI_F3x_S_bb;
  Double I_NAI_F2xy_Pz_bb = I_NAI_G2xyz_S_bb+ABZ*I_NAI_F2xy_S_bb;
  Double I_NAI_F2xz_Pz_bb = I_NAI_G2x2z_S_bb+ABZ*I_NAI_F2xz_S_bb;
  Double I_NAI_Fx2y_Pz_bb = I_NAI_Gx2yz_S_bb+ABZ*I_NAI_Fx2y_S_bb;
  Double I_NAI_Fxyz_Pz_bb = I_NAI_Gxy2z_S_bb+ABZ*I_NAI_Fxyz_S_bb;
  Double I_NAI_Fx2z_Pz_bb = I_NAI_Gx3z_S_bb+ABZ*I_NAI_Fx2z_S_bb;
  Double I_NAI_F3y_Pz_bb = I_NAI_G3yz_S_bb+ABZ*I_NAI_F3y_S_bb;
  Double I_NAI_F2yz_Pz_bb = I_NAI_G2y2z_S_bb+ABZ*I_NAI_F2yz_S_bb;
  Double I_NAI_Fy2z_Pz_bb = I_NAI_Gy3z_S_bb+ABZ*I_NAI_Fy2z_S_bb;
  Double I_NAI_F3z_Pz_bb = I_NAI_G4z_S_bb+ABZ*I_NAI_F3z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_bb
   * RHS shell quartet name: SQ_NAI_D_P_bb
   ************************************************************/
  Double I_NAI_D2x_D2x_bb = I_NAI_F3x_Px_bb+ABX*I_NAI_D2x_Px_bb;
  Double I_NAI_Dxy_D2x_bb = I_NAI_F2xy_Px_bb+ABX*I_NAI_Dxy_Px_bb;
  Double I_NAI_Dxz_D2x_bb = I_NAI_F2xz_Px_bb+ABX*I_NAI_Dxz_Px_bb;
  Double I_NAI_D2y_D2x_bb = I_NAI_Fx2y_Px_bb+ABX*I_NAI_D2y_Px_bb;
  Double I_NAI_Dyz_D2x_bb = I_NAI_Fxyz_Px_bb+ABX*I_NAI_Dyz_Px_bb;
  Double I_NAI_D2z_D2x_bb = I_NAI_Fx2z_Px_bb+ABX*I_NAI_D2z_Px_bb;
  Double I_NAI_D2x_D2y_bb = I_NAI_F2xy_Py_bb+ABY*I_NAI_D2x_Py_bb;
  Double I_NAI_Dxy_D2y_bb = I_NAI_Fx2y_Py_bb+ABY*I_NAI_Dxy_Py_bb;
  Double I_NAI_Dxz_D2y_bb = I_NAI_Fxyz_Py_bb+ABY*I_NAI_Dxz_Py_bb;
  Double I_NAI_D2y_D2y_bb = I_NAI_F3y_Py_bb+ABY*I_NAI_D2y_Py_bb;
  Double I_NAI_Dyz_D2y_bb = I_NAI_F2yz_Py_bb+ABY*I_NAI_Dyz_Py_bb;
  Double I_NAI_D2z_D2y_bb = I_NAI_Fy2z_Py_bb+ABY*I_NAI_D2z_Py_bb;
  Double I_NAI_D2x_D2z_bb = I_NAI_F2xz_Pz_bb+ABZ*I_NAI_D2x_Pz_bb;
  Double I_NAI_Dxy_D2z_bb = I_NAI_Fxyz_Pz_bb+ABZ*I_NAI_Dxy_Pz_bb;
  Double I_NAI_Dxz_D2z_bb = I_NAI_Fx2z_Pz_bb+ABZ*I_NAI_Dxz_Pz_bb;
  Double I_NAI_D2y_D2z_bb = I_NAI_F2yz_Pz_bb+ABZ*I_NAI_D2y_Pz_bb;
  Double I_NAI_Dyz_D2z_bb = I_NAI_Fy2z_Pz_bb+ABZ*I_NAI_Dyz_Pz_bb;
  Double I_NAI_D2z_D2z_bb = I_NAI_F3z_Pz_bb+ABZ*I_NAI_D2z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
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
  Double I_NAI_G3xy_Pz_bb = I_NAI_H3xyz_S_bb+ABZ*I_NAI_G3xy_S_bb;
  Double I_NAI_G3xz_Pz_bb = I_NAI_H3x2z_S_bb+ABZ*I_NAI_G3xz_S_bb;
  Double I_NAI_G2x2y_Pz_bb = I_NAI_H2x2yz_S_bb+ABZ*I_NAI_G2x2y_S_bb;
  Double I_NAI_G2xyz_Pz_bb = I_NAI_H2xy2z_S_bb+ABZ*I_NAI_G2xyz_S_bb;
  Double I_NAI_G2x2z_Pz_bb = I_NAI_H2x3z_S_bb+ABZ*I_NAI_G2x2z_S_bb;
  Double I_NAI_Gx3y_Pz_bb = I_NAI_Hx3yz_S_bb+ABZ*I_NAI_Gx3y_S_bb;
  Double I_NAI_Gx2yz_Pz_bb = I_NAI_Hx2y2z_S_bb+ABZ*I_NAI_Gx2yz_S_bb;
  Double I_NAI_Gxy2z_Pz_bb = I_NAI_Hxy3z_S_bb+ABZ*I_NAI_Gxy2z_S_bb;
  Double I_NAI_Gx3z_Pz_bb = I_NAI_Hx4z_S_bb+ABZ*I_NAI_Gx3z_S_bb;
  Double I_NAI_G3yz_Pz_bb = I_NAI_H3y2z_S_bb+ABZ*I_NAI_G3yz_S_bb;
  Double I_NAI_G2y2z_Pz_bb = I_NAI_H2y3z_S_bb+ABZ*I_NAI_G2y2z_S_bb;
  Double I_NAI_Gy3z_Pz_bb = I_NAI_Hy4z_S_bb+ABZ*I_NAI_Gy3z_S_bb;
  Double I_NAI_G4z_Pz_bb = I_NAI_H5z_S_bb+ABZ*I_NAI_G4z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_bb
   * RHS shell quartet name: SQ_NAI_F_P_bb
   ************************************************************/
  Double I_NAI_F3x_D2x_bb = I_NAI_G4x_Px_bb+ABX*I_NAI_F3x_Px_bb;
  Double I_NAI_F2xy_D2x_bb = I_NAI_G3xy_Px_bb+ABX*I_NAI_F2xy_Px_bb;
  Double I_NAI_F2xz_D2x_bb = I_NAI_G3xz_Px_bb+ABX*I_NAI_F2xz_Px_bb;
  Double I_NAI_Fx2y_D2x_bb = I_NAI_G2x2y_Px_bb+ABX*I_NAI_Fx2y_Px_bb;
  Double I_NAI_Fxyz_D2x_bb = I_NAI_G2xyz_Px_bb+ABX*I_NAI_Fxyz_Px_bb;
  Double I_NAI_Fx2z_D2x_bb = I_NAI_G2x2z_Px_bb+ABX*I_NAI_Fx2z_Px_bb;
  Double I_NAI_F3y_D2x_bb = I_NAI_Gx3y_Px_bb+ABX*I_NAI_F3y_Px_bb;
  Double I_NAI_F2yz_D2x_bb = I_NAI_Gx2yz_Px_bb+ABX*I_NAI_F2yz_Px_bb;
  Double I_NAI_Fy2z_D2x_bb = I_NAI_Gxy2z_Px_bb+ABX*I_NAI_Fy2z_Px_bb;
  Double I_NAI_F3z_D2x_bb = I_NAI_Gx3z_Px_bb+ABX*I_NAI_F3z_Px_bb;
  Double I_NAI_F3x_D2y_bb = I_NAI_G3xy_Py_bb+ABY*I_NAI_F3x_Py_bb;
  Double I_NAI_F2xy_D2y_bb = I_NAI_G2x2y_Py_bb+ABY*I_NAI_F2xy_Py_bb;
  Double I_NAI_F2xz_D2y_bb = I_NAI_G2xyz_Py_bb+ABY*I_NAI_F2xz_Py_bb;
  Double I_NAI_Fx2y_D2y_bb = I_NAI_Gx3y_Py_bb+ABY*I_NAI_Fx2y_Py_bb;
  Double I_NAI_Fxyz_D2y_bb = I_NAI_Gx2yz_Py_bb+ABY*I_NAI_Fxyz_Py_bb;
  Double I_NAI_Fx2z_D2y_bb = I_NAI_Gxy2z_Py_bb+ABY*I_NAI_Fx2z_Py_bb;
  Double I_NAI_F3y_D2y_bb = I_NAI_G4y_Py_bb+ABY*I_NAI_F3y_Py_bb;
  Double I_NAI_F2yz_D2y_bb = I_NAI_G3yz_Py_bb+ABY*I_NAI_F2yz_Py_bb;
  Double I_NAI_Fy2z_D2y_bb = I_NAI_G2y2z_Py_bb+ABY*I_NAI_Fy2z_Py_bb;
  Double I_NAI_F3z_D2y_bb = I_NAI_Gy3z_Py_bb+ABY*I_NAI_F3z_Py_bb;
  Double I_NAI_F3x_D2z_bb = I_NAI_G3xz_Pz_bb+ABZ*I_NAI_F3x_Pz_bb;
  Double I_NAI_F2xy_D2z_bb = I_NAI_G2xyz_Pz_bb+ABZ*I_NAI_F2xy_Pz_bb;
  Double I_NAI_F2xz_D2z_bb = I_NAI_G2x2z_Pz_bb+ABZ*I_NAI_F2xz_Pz_bb;
  Double I_NAI_Fx2y_D2z_bb = I_NAI_Gx2yz_Pz_bb+ABZ*I_NAI_Fx2y_Pz_bb;
  Double I_NAI_Fxyz_D2z_bb = I_NAI_Gxy2z_Pz_bb+ABZ*I_NAI_Fxyz_Pz_bb;
  Double I_NAI_Fx2z_D2z_bb = I_NAI_Gx3z_Pz_bb+ABZ*I_NAI_Fx2z_Pz_bb;
  Double I_NAI_F3y_D2z_bb = I_NAI_G3yz_Pz_bb+ABZ*I_NAI_F3y_Pz_bb;
  Double I_NAI_F2yz_D2z_bb = I_NAI_G2y2z_Pz_bb+ABZ*I_NAI_F2yz_Pz_bb;
  Double I_NAI_Fy2z_D2z_bb = I_NAI_Gy3z_Pz_bb+ABZ*I_NAI_Fy2z_Pz_bb;
  Double I_NAI_F3z_D2z_bb = I_NAI_G4z_Pz_bb+ABZ*I_NAI_F3z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_F_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_bb
   * RHS shell quartet name: SQ_NAI_D_D_bb
   ************************************************************/
  Double I_NAI_D2x_F3x_bb = I_NAI_F3x_D2x_bb+ABX*I_NAI_D2x_D2x_bb;
  Double I_NAI_Dxy_F3x_bb = I_NAI_F2xy_D2x_bb+ABX*I_NAI_Dxy_D2x_bb;
  Double I_NAI_Dxz_F3x_bb = I_NAI_F2xz_D2x_bb+ABX*I_NAI_Dxz_D2x_bb;
  Double I_NAI_D2y_F3x_bb = I_NAI_Fx2y_D2x_bb+ABX*I_NAI_D2y_D2x_bb;
  Double I_NAI_Dyz_F3x_bb = I_NAI_Fxyz_D2x_bb+ABX*I_NAI_Dyz_D2x_bb;
  Double I_NAI_D2z_F3x_bb = I_NAI_Fx2z_D2x_bb+ABX*I_NAI_D2z_D2x_bb;
  Double I_NAI_D2x_F2xy_bb = I_NAI_F2xy_D2x_bb+ABY*I_NAI_D2x_D2x_bb;
  Double I_NAI_Dxy_F2xy_bb = I_NAI_Fx2y_D2x_bb+ABY*I_NAI_Dxy_D2x_bb;
  Double I_NAI_Dxz_F2xy_bb = I_NAI_Fxyz_D2x_bb+ABY*I_NAI_Dxz_D2x_bb;
  Double I_NAI_D2y_F2xy_bb = I_NAI_F3y_D2x_bb+ABY*I_NAI_D2y_D2x_bb;
  Double I_NAI_Dyz_F2xy_bb = I_NAI_F2yz_D2x_bb+ABY*I_NAI_Dyz_D2x_bb;
  Double I_NAI_D2z_F2xy_bb = I_NAI_Fy2z_D2x_bb+ABY*I_NAI_D2z_D2x_bb;
  Double I_NAI_D2x_F2xz_bb = I_NAI_F2xz_D2x_bb+ABZ*I_NAI_D2x_D2x_bb;
  Double I_NAI_Dxy_F2xz_bb = I_NAI_Fxyz_D2x_bb+ABZ*I_NAI_Dxy_D2x_bb;
  Double I_NAI_Dxz_F2xz_bb = I_NAI_Fx2z_D2x_bb+ABZ*I_NAI_Dxz_D2x_bb;
  Double I_NAI_D2y_F2xz_bb = I_NAI_F2yz_D2x_bb+ABZ*I_NAI_D2y_D2x_bb;
  Double I_NAI_Dyz_F2xz_bb = I_NAI_Fy2z_D2x_bb+ABZ*I_NAI_Dyz_D2x_bb;
  Double I_NAI_D2z_F2xz_bb = I_NAI_F3z_D2x_bb+ABZ*I_NAI_D2z_D2x_bb;
  Double I_NAI_D2x_Fx2y_bb = I_NAI_F3x_D2y_bb+ABX*I_NAI_D2x_D2y_bb;
  Double I_NAI_Dxy_Fx2y_bb = I_NAI_F2xy_D2y_bb+ABX*I_NAI_Dxy_D2y_bb;
  Double I_NAI_Dxz_Fx2y_bb = I_NAI_F2xz_D2y_bb+ABX*I_NAI_Dxz_D2y_bb;
  Double I_NAI_D2y_Fx2y_bb = I_NAI_Fx2y_D2y_bb+ABX*I_NAI_D2y_D2y_bb;
  Double I_NAI_Dyz_Fx2y_bb = I_NAI_Fxyz_D2y_bb+ABX*I_NAI_Dyz_D2y_bb;
  Double I_NAI_D2z_Fx2y_bb = I_NAI_Fx2z_D2y_bb+ABX*I_NAI_D2z_D2y_bb;
  Double I_NAI_D2x_Fx2z_bb = I_NAI_F3x_D2z_bb+ABX*I_NAI_D2x_D2z_bb;
  Double I_NAI_Dxy_Fx2z_bb = I_NAI_F2xy_D2z_bb+ABX*I_NAI_Dxy_D2z_bb;
  Double I_NAI_Dxz_Fx2z_bb = I_NAI_F2xz_D2z_bb+ABX*I_NAI_Dxz_D2z_bb;
  Double I_NAI_D2y_Fx2z_bb = I_NAI_Fx2y_D2z_bb+ABX*I_NAI_D2y_D2z_bb;
  Double I_NAI_Dyz_Fx2z_bb = I_NAI_Fxyz_D2z_bb+ABX*I_NAI_Dyz_D2z_bb;
  Double I_NAI_D2z_Fx2z_bb = I_NAI_Fx2z_D2z_bb+ABX*I_NAI_D2z_D2z_bb;
  Double I_NAI_D2x_F3y_bb = I_NAI_F2xy_D2y_bb+ABY*I_NAI_D2x_D2y_bb;
  Double I_NAI_Dxy_F3y_bb = I_NAI_Fx2y_D2y_bb+ABY*I_NAI_Dxy_D2y_bb;
  Double I_NAI_Dxz_F3y_bb = I_NAI_Fxyz_D2y_bb+ABY*I_NAI_Dxz_D2y_bb;
  Double I_NAI_D2y_F3y_bb = I_NAI_F3y_D2y_bb+ABY*I_NAI_D2y_D2y_bb;
  Double I_NAI_Dyz_F3y_bb = I_NAI_F2yz_D2y_bb+ABY*I_NAI_Dyz_D2y_bb;
  Double I_NAI_D2z_F3y_bb = I_NAI_Fy2z_D2y_bb+ABY*I_NAI_D2z_D2y_bb;
  Double I_NAI_D2x_F2yz_bb = I_NAI_F2xz_D2y_bb+ABZ*I_NAI_D2x_D2y_bb;
  Double I_NAI_Dxy_F2yz_bb = I_NAI_Fxyz_D2y_bb+ABZ*I_NAI_Dxy_D2y_bb;
  Double I_NAI_Dxz_F2yz_bb = I_NAI_Fx2z_D2y_bb+ABZ*I_NAI_Dxz_D2y_bb;
  Double I_NAI_D2y_F2yz_bb = I_NAI_F2yz_D2y_bb+ABZ*I_NAI_D2y_D2y_bb;
  Double I_NAI_Dyz_F2yz_bb = I_NAI_Fy2z_D2y_bb+ABZ*I_NAI_Dyz_D2y_bb;
  Double I_NAI_D2z_F2yz_bb = I_NAI_F3z_D2y_bb+ABZ*I_NAI_D2z_D2y_bb;
  Double I_NAI_D2x_F3z_bb = I_NAI_F2xz_D2z_bb+ABZ*I_NAI_D2x_D2z_bb;
  Double I_NAI_Dxy_F3z_bb = I_NAI_Fxyz_D2z_bb+ABZ*I_NAI_Dxy_D2z_bb;
  Double I_NAI_Dxz_F3z_bb = I_NAI_Fx2z_D2z_bb+ABZ*I_NAI_Dxz_D2z_bb;
  Double I_NAI_D2y_F3z_bb = I_NAI_F2yz_D2z_bb+ABZ*I_NAI_D2y_D2z_bb;
  Double I_NAI_Dyz_F3z_bb = I_NAI_Fy2z_D2z_bb+ABZ*I_NAI_Dyz_D2z_bb;
  Double I_NAI_D2z_F3z_bb = I_NAI_F3z_D2z_bb+ABZ*I_NAI_D2z_D2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 21 integrals are omitted 
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
  Double I_NAI_H3x2y_Py_bb = I_NAI_I3x3y_S_bb+ABY*I_NAI_H3x2y_S_bb;
  Double I_NAI_H3xyz_Py_bb = I_NAI_I3x2yz_S_bb+ABY*I_NAI_H3xyz_S_bb;
  Double I_NAI_H2x3y_Py_bb = I_NAI_I2x4y_S_bb+ABY*I_NAI_H2x3y_S_bb;
  Double I_NAI_H2x2yz_Py_bb = I_NAI_I2x3yz_S_bb+ABY*I_NAI_H2x2yz_S_bb;
  Double I_NAI_H2xy2z_Py_bb = I_NAI_I2x2y2z_S_bb+ABY*I_NAI_H2xy2z_S_bb;
  Double I_NAI_Hx4y_Py_bb = I_NAI_Ix5y_S_bb+ABY*I_NAI_Hx4y_S_bb;
  Double I_NAI_Hx3yz_Py_bb = I_NAI_Ix4yz_S_bb+ABY*I_NAI_Hx3yz_S_bb;
  Double I_NAI_Hx2y2z_Py_bb = I_NAI_Ix3y2z_S_bb+ABY*I_NAI_Hx2y2z_S_bb;
  Double I_NAI_Hxy3z_Py_bb = I_NAI_Ix2y3z_S_bb+ABY*I_NAI_Hxy3z_S_bb;
  Double I_NAI_H5y_Py_bb = I_NAI_I6y_S_bb+ABY*I_NAI_H5y_S_bb;
  Double I_NAI_H4yz_Py_bb = I_NAI_I5yz_S_bb+ABY*I_NAI_H4yz_S_bb;
  Double I_NAI_H3y2z_Py_bb = I_NAI_I4y2z_S_bb+ABY*I_NAI_H3y2z_S_bb;
  Double I_NAI_H2y3z_Py_bb = I_NAI_I3y3z_S_bb+ABY*I_NAI_H2y3z_S_bb;
  Double I_NAI_Hy4z_Py_bb = I_NAI_I2y4z_S_bb+ABY*I_NAI_Hy4z_S_bb;
  Double I_NAI_H3xyz_Pz_bb = I_NAI_I3xy2z_S_bb+ABZ*I_NAI_H3xyz_S_bb;
  Double I_NAI_H3x2z_Pz_bb = I_NAI_I3x3z_S_bb+ABZ*I_NAI_H3x2z_S_bb;
  Double I_NAI_H2x2yz_Pz_bb = I_NAI_I2x2y2z_S_bb+ABZ*I_NAI_H2x2yz_S_bb;
  Double I_NAI_H2xy2z_Pz_bb = I_NAI_I2xy3z_S_bb+ABZ*I_NAI_H2xy2z_S_bb;
  Double I_NAI_H2x3z_Pz_bb = I_NAI_I2x4z_S_bb+ABZ*I_NAI_H2x3z_S_bb;
  Double I_NAI_Hx3yz_Pz_bb = I_NAI_Ix3y2z_S_bb+ABZ*I_NAI_Hx3yz_S_bb;
  Double I_NAI_Hx2y2z_Pz_bb = I_NAI_Ix2y3z_S_bb+ABZ*I_NAI_Hx2y2z_S_bb;
  Double I_NAI_Hxy3z_Pz_bb = I_NAI_Ixy4z_S_bb+ABZ*I_NAI_Hxy3z_S_bb;
  Double I_NAI_Hx4z_Pz_bb = I_NAI_Ix5z_S_bb+ABZ*I_NAI_Hx4z_S_bb;
  Double I_NAI_H3y2z_Pz_bb = I_NAI_I3y3z_S_bb+ABZ*I_NAI_H3y2z_S_bb;
  Double I_NAI_H2y3z_Pz_bb = I_NAI_I2y4z_S_bb+ABZ*I_NAI_H2y3z_S_bb;
  Double I_NAI_Hy4z_Pz_bb = I_NAI_Iy5z_S_bb+ABZ*I_NAI_Hy4z_S_bb;
  Double I_NAI_H5z_Pz_bb = I_NAI_I6z_S_bb+ABZ*I_NAI_H5z_S_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 48 integrals are omitted 
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
  Double I_NAI_G3xy_D2z_bb = I_NAI_H3xyz_Pz_bb+ABZ*I_NAI_G3xy_Pz_bb;
  Double I_NAI_G3xz_D2z_bb = I_NAI_H3x2z_Pz_bb+ABZ*I_NAI_G3xz_Pz_bb;
  Double I_NAI_G2x2y_D2z_bb = I_NAI_H2x2yz_Pz_bb+ABZ*I_NAI_G2x2y_Pz_bb;
  Double I_NAI_G2xyz_D2z_bb = I_NAI_H2xy2z_Pz_bb+ABZ*I_NAI_G2xyz_Pz_bb;
  Double I_NAI_G2x2z_D2z_bb = I_NAI_H2x3z_Pz_bb+ABZ*I_NAI_G2x2z_Pz_bb;
  Double I_NAI_Gx3y_D2z_bb = I_NAI_Hx3yz_Pz_bb+ABZ*I_NAI_Gx3y_Pz_bb;
  Double I_NAI_Gx2yz_D2z_bb = I_NAI_Hx2y2z_Pz_bb+ABZ*I_NAI_Gx2yz_Pz_bb;
  Double I_NAI_Gxy2z_D2z_bb = I_NAI_Hxy3z_Pz_bb+ABZ*I_NAI_Gxy2z_Pz_bb;
  Double I_NAI_Gx3z_D2z_bb = I_NAI_Hx4z_Pz_bb+ABZ*I_NAI_Gx3z_Pz_bb;
  Double I_NAI_G3yz_D2z_bb = I_NAI_H3y2z_Pz_bb+ABZ*I_NAI_G3yz_Pz_bb;
  Double I_NAI_G2y2z_D2z_bb = I_NAI_H2y3z_Pz_bb+ABZ*I_NAI_G2y2z_Pz_bb;
  Double I_NAI_Gy3z_D2z_bb = I_NAI_Hy4z_Pz_bb+ABZ*I_NAI_Gy3z_Pz_bb;
  Double I_NAI_G4z_D2z_bb = I_NAI_H5z_Pz_bb+ABZ*I_NAI_G4z_Pz_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 37 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_bb
   * RHS shell quartet name: SQ_NAI_F_D_bb
   ************************************************************/
  Double I_NAI_F3x_F3x_bb = I_NAI_G4x_D2x_bb+ABX*I_NAI_F3x_D2x_bb;
  Double I_NAI_F2xy_F3x_bb = I_NAI_G3xy_D2x_bb+ABX*I_NAI_F2xy_D2x_bb;
  Double I_NAI_F2xz_F3x_bb = I_NAI_G3xz_D2x_bb+ABX*I_NAI_F2xz_D2x_bb;
  Double I_NAI_Fx2y_F3x_bb = I_NAI_G2x2y_D2x_bb+ABX*I_NAI_Fx2y_D2x_bb;
  Double I_NAI_Fxyz_F3x_bb = I_NAI_G2xyz_D2x_bb+ABX*I_NAI_Fxyz_D2x_bb;
  Double I_NAI_Fx2z_F3x_bb = I_NAI_G2x2z_D2x_bb+ABX*I_NAI_Fx2z_D2x_bb;
  Double I_NAI_F3y_F3x_bb = I_NAI_Gx3y_D2x_bb+ABX*I_NAI_F3y_D2x_bb;
  Double I_NAI_F2yz_F3x_bb = I_NAI_Gx2yz_D2x_bb+ABX*I_NAI_F2yz_D2x_bb;
  Double I_NAI_Fy2z_F3x_bb = I_NAI_Gxy2z_D2x_bb+ABX*I_NAI_Fy2z_D2x_bb;
  Double I_NAI_F3z_F3x_bb = I_NAI_Gx3z_D2x_bb+ABX*I_NAI_F3z_D2x_bb;
  Double I_NAI_F2xy_F2xy_bb = I_NAI_G2x2y_D2x_bb+ABY*I_NAI_F2xy_D2x_bb;
  Double I_NAI_F2xz_F2xy_bb = I_NAI_G2xyz_D2x_bb+ABY*I_NAI_F2xz_D2x_bb;
  Double I_NAI_Fx2y_F2xy_bb = I_NAI_Gx3y_D2x_bb+ABY*I_NAI_Fx2y_D2x_bb;
  Double I_NAI_Fxyz_F2xy_bb = I_NAI_Gx2yz_D2x_bb+ABY*I_NAI_Fxyz_D2x_bb;
  Double I_NAI_Fx2z_F2xy_bb = I_NAI_Gxy2z_D2x_bb+ABY*I_NAI_Fx2z_D2x_bb;
  Double I_NAI_F3y_F2xy_bb = I_NAI_G4y_D2x_bb+ABY*I_NAI_F3y_D2x_bb;
  Double I_NAI_F2yz_F2xy_bb = I_NAI_G3yz_D2x_bb+ABY*I_NAI_F2yz_D2x_bb;
  Double I_NAI_Fy2z_F2xy_bb = I_NAI_G2y2z_D2x_bb+ABY*I_NAI_Fy2z_D2x_bb;
  Double I_NAI_F3z_F2xy_bb = I_NAI_Gy3z_D2x_bb+ABY*I_NAI_F3z_D2x_bb;
  Double I_NAI_F2xz_F2xz_bb = I_NAI_G2x2z_D2x_bb+ABZ*I_NAI_F2xz_D2x_bb;
  Double I_NAI_Fxyz_F2xz_bb = I_NAI_Gxy2z_D2x_bb+ABZ*I_NAI_Fxyz_D2x_bb;
  Double I_NAI_Fx2z_F2xz_bb = I_NAI_Gx3z_D2x_bb+ABZ*I_NAI_Fx2z_D2x_bb;
  Double I_NAI_F2yz_F2xz_bb = I_NAI_G2y2z_D2x_bb+ABZ*I_NAI_F2yz_D2x_bb;
  Double I_NAI_Fy2z_F2xz_bb = I_NAI_Gy3z_D2x_bb+ABZ*I_NAI_Fy2z_D2x_bb;
  Double I_NAI_F3z_F2xz_bb = I_NAI_G4z_D2x_bb+ABZ*I_NAI_F3z_D2x_bb;
  Double I_NAI_F2xz_Fx2y_bb = I_NAI_G3xz_D2y_bb+ABX*I_NAI_F2xz_D2y_bb;
  Double I_NAI_Fxyz_Fx2y_bb = I_NAI_G2xyz_D2y_bb+ABX*I_NAI_Fxyz_D2y_bb;
  Double I_NAI_Fx2z_Fx2y_bb = I_NAI_G2x2z_D2y_bb+ABX*I_NAI_Fx2z_D2y_bb;
  Double I_NAI_F2yz_Fx2y_bb = I_NAI_Gx2yz_D2y_bb+ABX*I_NAI_F2yz_D2y_bb;
  Double I_NAI_Fy2z_Fx2y_bb = I_NAI_Gxy2z_D2y_bb+ABX*I_NAI_Fy2z_D2y_bb;
  Double I_NAI_F3z_Fx2y_bb = I_NAI_Gx3z_D2y_bb+ABX*I_NAI_F3z_D2y_bb;
  Double I_NAI_F2xy_Fx2z_bb = I_NAI_G3xy_D2z_bb+ABX*I_NAI_F2xy_D2z_bb;
  Double I_NAI_Fx2y_Fx2z_bb = I_NAI_G2x2y_D2z_bb+ABX*I_NAI_Fx2y_D2z_bb;
  Double I_NAI_Fxyz_Fx2z_bb = I_NAI_G2xyz_D2z_bb+ABX*I_NAI_Fxyz_D2z_bb;
  Double I_NAI_F3y_Fx2z_bb = I_NAI_Gx3y_D2z_bb+ABX*I_NAI_F3y_D2z_bb;
  Double I_NAI_F2yz_Fx2z_bb = I_NAI_Gx2yz_D2z_bb+ABX*I_NAI_F2yz_D2z_bb;
  Double I_NAI_Fy2z_Fx2z_bb = I_NAI_Gxy2z_D2z_bb+ABX*I_NAI_Fy2z_D2z_bb;
  Double I_NAI_F3x_F3y_bb = I_NAI_G3xy_D2y_bb+ABY*I_NAI_F3x_D2y_bb;
  Double I_NAI_F2xy_F3y_bb = I_NAI_G2x2y_D2y_bb+ABY*I_NAI_F2xy_D2y_bb;
  Double I_NAI_F2xz_F3y_bb = I_NAI_G2xyz_D2y_bb+ABY*I_NAI_F2xz_D2y_bb;
  Double I_NAI_Fx2y_F3y_bb = I_NAI_Gx3y_D2y_bb+ABY*I_NAI_Fx2y_D2y_bb;
  Double I_NAI_Fxyz_F3y_bb = I_NAI_Gx2yz_D2y_bb+ABY*I_NAI_Fxyz_D2y_bb;
  Double I_NAI_Fx2z_F3y_bb = I_NAI_Gxy2z_D2y_bb+ABY*I_NAI_Fx2z_D2y_bb;
  Double I_NAI_F3y_F3y_bb = I_NAI_G4y_D2y_bb+ABY*I_NAI_F3y_D2y_bb;
  Double I_NAI_F2yz_F3y_bb = I_NAI_G3yz_D2y_bb+ABY*I_NAI_F2yz_D2y_bb;
  Double I_NAI_Fy2z_F3y_bb = I_NAI_G2y2z_D2y_bb+ABY*I_NAI_Fy2z_D2y_bb;
  Double I_NAI_F3z_F3y_bb = I_NAI_Gy3z_D2y_bb+ABY*I_NAI_F3z_D2y_bb;
  Double I_NAI_F2xz_F2yz_bb = I_NAI_G2x2z_D2y_bb+ABZ*I_NAI_F2xz_D2y_bb;
  Double I_NAI_Fxyz_F2yz_bb = I_NAI_Gxy2z_D2y_bb+ABZ*I_NAI_Fxyz_D2y_bb;
  Double I_NAI_Fx2z_F2yz_bb = I_NAI_Gx3z_D2y_bb+ABZ*I_NAI_Fx2z_D2y_bb;
  Double I_NAI_F2yz_F2yz_bb = I_NAI_G2y2z_D2y_bb+ABZ*I_NAI_F2yz_D2y_bb;
  Double I_NAI_Fy2z_F2yz_bb = I_NAI_Gy3z_D2y_bb+ABZ*I_NAI_Fy2z_D2y_bb;
  Double I_NAI_F3z_F2yz_bb = I_NAI_G4z_D2y_bb+ABZ*I_NAI_F3z_D2y_bb;
  Double I_NAI_F3x_F3z_bb = I_NAI_G3xz_D2z_bb+ABZ*I_NAI_F3x_D2z_bb;
  Double I_NAI_F2xy_F3z_bb = I_NAI_G2xyz_D2z_bb+ABZ*I_NAI_F2xy_D2z_bb;
  Double I_NAI_F2xz_F3z_bb = I_NAI_G2x2z_D2z_bb+ABZ*I_NAI_F2xz_D2z_bb;
  Double I_NAI_Fx2y_F3z_bb = I_NAI_Gx2yz_D2z_bb+ABZ*I_NAI_Fx2y_D2z_bb;
  Double I_NAI_Fxyz_F3z_bb = I_NAI_Gxy2z_D2z_bb+ABZ*I_NAI_Fxyz_D2z_bb;
  Double I_NAI_Fx2z_F3z_bb = I_NAI_Gx3z_D2z_bb+ABZ*I_NAI_Fx2z_D2z_bb;
  Double I_NAI_F3y_F3z_bb = I_NAI_G3yz_D2z_bb+ABZ*I_NAI_F3y_D2z_bb;
  Double I_NAI_F2yz_F3z_bb = I_NAI_G2y2z_D2z_bb+ABZ*I_NAI_F2yz_D2z_bb;
  Double I_NAI_Fy2z_F3z_bb = I_NAI_Gy3z_D2z_bb+ABZ*I_NAI_Fy2z_D2z_bb;
  Double I_NAI_F3z_F3z_bb = I_NAI_G4z_D2z_bb+ABZ*I_NAI_F3z_D2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_G_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_bb
   * RHS shell quartet name: SQ_NAI_D_F_bb
   ************************************************************/
  Double I_NAI_D2x_G4x_bb = I_NAI_F3x_F3x_bb+ABX*I_NAI_D2x_F3x_bb;
  Double I_NAI_Dxy_G4x_bb = I_NAI_F2xy_F3x_bb+ABX*I_NAI_Dxy_F3x_bb;
  Double I_NAI_Dxz_G4x_bb = I_NAI_F2xz_F3x_bb+ABX*I_NAI_Dxz_F3x_bb;
  Double I_NAI_D2y_G4x_bb = I_NAI_Fx2y_F3x_bb+ABX*I_NAI_D2y_F3x_bb;
  Double I_NAI_Dyz_G4x_bb = I_NAI_Fxyz_F3x_bb+ABX*I_NAI_Dyz_F3x_bb;
  Double I_NAI_D2z_G4x_bb = I_NAI_Fx2z_F3x_bb+ABX*I_NAI_D2z_F3x_bb;
  Double I_NAI_D2x_G3xy_bb = I_NAI_F2xy_F3x_bb+ABY*I_NAI_D2x_F3x_bb;
  Double I_NAI_Dxy_G3xy_bb = I_NAI_Fx2y_F3x_bb+ABY*I_NAI_Dxy_F3x_bb;
  Double I_NAI_Dxz_G3xy_bb = I_NAI_Fxyz_F3x_bb+ABY*I_NAI_Dxz_F3x_bb;
  Double I_NAI_D2y_G3xy_bb = I_NAI_F3y_F3x_bb+ABY*I_NAI_D2y_F3x_bb;
  Double I_NAI_Dyz_G3xy_bb = I_NAI_F2yz_F3x_bb+ABY*I_NAI_Dyz_F3x_bb;
  Double I_NAI_D2z_G3xy_bb = I_NAI_Fy2z_F3x_bb+ABY*I_NAI_D2z_F3x_bb;
  Double I_NAI_D2x_G3xz_bb = I_NAI_F2xz_F3x_bb+ABZ*I_NAI_D2x_F3x_bb;
  Double I_NAI_Dxy_G3xz_bb = I_NAI_Fxyz_F3x_bb+ABZ*I_NAI_Dxy_F3x_bb;
  Double I_NAI_Dxz_G3xz_bb = I_NAI_Fx2z_F3x_bb+ABZ*I_NAI_Dxz_F3x_bb;
  Double I_NAI_D2y_G3xz_bb = I_NAI_F2yz_F3x_bb+ABZ*I_NAI_D2y_F3x_bb;
  Double I_NAI_Dyz_G3xz_bb = I_NAI_Fy2z_F3x_bb+ABZ*I_NAI_Dyz_F3x_bb;
  Double I_NAI_D2z_G3xz_bb = I_NAI_F3z_F3x_bb+ABZ*I_NAI_D2z_F3x_bb;
  Double I_NAI_D2x_G2x2y_bb = I_NAI_F2xy_F2xy_bb+ABY*I_NAI_D2x_F2xy_bb;
  Double I_NAI_Dxy_G2x2y_bb = I_NAI_Fx2y_F2xy_bb+ABY*I_NAI_Dxy_F2xy_bb;
  Double I_NAI_Dxz_G2x2y_bb = I_NAI_Fxyz_F2xy_bb+ABY*I_NAI_Dxz_F2xy_bb;
  Double I_NAI_D2y_G2x2y_bb = I_NAI_F3y_F2xy_bb+ABY*I_NAI_D2y_F2xy_bb;
  Double I_NAI_Dyz_G2x2y_bb = I_NAI_F2yz_F2xy_bb+ABY*I_NAI_Dyz_F2xy_bb;
  Double I_NAI_D2z_G2x2y_bb = I_NAI_Fy2z_F2xy_bb+ABY*I_NAI_D2z_F2xy_bb;
  Double I_NAI_D2x_G2xyz_bb = I_NAI_F2xz_F2xy_bb+ABZ*I_NAI_D2x_F2xy_bb;
  Double I_NAI_Dxy_G2xyz_bb = I_NAI_Fxyz_F2xy_bb+ABZ*I_NAI_Dxy_F2xy_bb;
  Double I_NAI_Dxz_G2xyz_bb = I_NAI_Fx2z_F2xy_bb+ABZ*I_NAI_Dxz_F2xy_bb;
  Double I_NAI_D2y_G2xyz_bb = I_NAI_F2yz_F2xy_bb+ABZ*I_NAI_D2y_F2xy_bb;
  Double I_NAI_Dyz_G2xyz_bb = I_NAI_Fy2z_F2xy_bb+ABZ*I_NAI_Dyz_F2xy_bb;
  Double I_NAI_D2z_G2xyz_bb = I_NAI_F3z_F2xy_bb+ABZ*I_NAI_D2z_F2xy_bb;
  Double I_NAI_D2x_G2x2z_bb = I_NAI_F2xz_F2xz_bb+ABZ*I_NAI_D2x_F2xz_bb;
  Double I_NAI_Dxy_G2x2z_bb = I_NAI_Fxyz_F2xz_bb+ABZ*I_NAI_Dxy_F2xz_bb;
  Double I_NAI_Dxz_G2x2z_bb = I_NAI_Fx2z_F2xz_bb+ABZ*I_NAI_Dxz_F2xz_bb;
  Double I_NAI_D2y_G2x2z_bb = I_NAI_F2yz_F2xz_bb+ABZ*I_NAI_D2y_F2xz_bb;
  Double I_NAI_Dyz_G2x2z_bb = I_NAI_Fy2z_F2xz_bb+ABZ*I_NAI_Dyz_F2xz_bb;
  Double I_NAI_D2z_G2x2z_bb = I_NAI_F3z_F2xz_bb+ABZ*I_NAI_D2z_F2xz_bb;
  Double I_NAI_D2x_Gx3y_bb = I_NAI_F3x_F3y_bb+ABX*I_NAI_D2x_F3y_bb;
  Double I_NAI_Dxy_Gx3y_bb = I_NAI_F2xy_F3y_bb+ABX*I_NAI_Dxy_F3y_bb;
  Double I_NAI_Dxz_Gx3y_bb = I_NAI_F2xz_F3y_bb+ABX*I_NAI_Dxz_F3y_bb;
  Double I_NAI_D2y_Gx3y_bb = I_NAI_Fx2y_F3y_bb+ABX*I_NAI_D2y_F3y_bb;
  Double I_NAI_Dyz_Gx3y_bb = I_NAI_Fxyz_F3y_bb+ABX*I_NAI_Dyz_F3y_bb;
  Double I_NAI_D2z_Gx3y_bb = I_NAI_Fx2z_F3y_bb+ABX*I_NAI_D2z_F3y_bb;
  Double I_NAI_D2x_Gx2yz_bb = I_NAI_F2xz_Fx2y_bb+ABZ*I_NAI_D2x_Fx2y_bb;
  Double I_NAI_Dxy_Gx2yz_bb = I_NAI_Fxyz_Fx2y_bb+ABZ*I_NAI_Dxy_Fx2y_bb;
  Double I_NAI_Dxz_Gx2yz_bb = I_NAI_Fx2z_Fx2y_bb+ABZ*I_NAI_Dxz_Fx2y_bb;
  Double I_NAI_D2y_Gx2yz_bb = I_NAI_F2yz_Fx2y_bb+ABZ*I_NAI_D2y_Fx2y_bb;
  Double I_NAI_Dyz_Gx2yz_bb = I_NAI_Fy2z_Fx2y_bb+ABZ*I_NAI_Dyz_Fx2y_bb;
  Double I_NAI_D2z_Gx2yz_bb = I_NAI_F3z_Fx2y_bb+ABZ*I_NAI_D2z_Fx2y_bb;
  Double I_NAI_D2x_Gxy2z_bb = I_NAI_F2xy_Fx2z_bb+ABY*I_NAI_D2x_Fx2z_bb;
  Double I_NAI_Dxy_Gxy2z_bb = I_NAI_Fx2y_Fx2z_bb+ABY*I_NAI_Dxy_Fx2z_bb;
  Double I_NAI_Dxz_Gxy2z_bb = I_NAI_Fxyz_Fx2z_bb+ABY*I_NAI_Dxz_Fx2z_bb;
  Double I_NAI_D2y_Gxy2z_bb = I_NAI_F3y_Fx2z_bb+ABY*I_NAI_D2y_Fx2z_bb;
  Double I_NAI_Dyz_Gxy2z_bb = I_NAI_F2yz_Fx2z_bb+ABY*I_NAI_Dyz_Fx2z_bb;
  Double I_NAI_D2z_Gxy2z_bb = I_NAI_Fy2z_Fx2z_bb+ABY*I_NAI_D2z_Fx2z_bb;
  Double I_NAI_D2x_Gx3z_bb = I_NAI_F3x_F3z_bb+ABX*I_NAI_D2x_F3z_bb;
  Double I_NAI_Dxy_Gx3z_bb = I_NAI_F2xy_F3z_bb+ABX*I_NAI_Dxy_F3z_bb;
  Double I_NAI_Dxz_Gx3z_bb = I_NAI_F2xz_F3z_bb+ABX*I_NAI_Dxz_F3z_bb;
  Double I_NAI_D2y_Gx3z_bb = I_NAI_Fx2y_F3z_bb+ABX*I_NAI_D2y_F3z_bb;
  Double I_NAI_Dyz_Gx3z_bb = I_NAI_Fxyz_F3z_bb+ABX*I_NAI_Dyz_F3z_bb;
  Double I_NAI_D2z_Gx3z_bb = I_NAI_Fx2z_F3z_bb+ABX*I_NAI_D2z_F3z_bb;
  Double I_NAI_D2x_G4y_bb = I_NAI_F2xy_F3y_bb+ABY*I_NAI_D2x_F3y_bb;
  Double I_NAI_Dxy_G4y_bb = I_NAI_Fx2y_F3y_bb+ABY*I_NAI_Dxy_F3y_bb;
  Double I_NAI_Dxz_G4y_bb = I_NAI_Fxyz_F3y_bb+ABY*I_NAI_Dxz_F3y_bb;
  Double I_NAI_D2y_G4y_bb = I_NAI_F3y_F3y_bb+ABY*I_NAI_D2y_F3y_bb;
  Double I_NAI_Dyz_G4y_bb = I_NAI_F2yz_F3y_bb+ABY*I_NAI_Dyz_F3y_bb;
  Double I_NAI_D2z_G4y_bb = I_NAI_Fy2z_F3y_bb+ABY*I_NAI_D2z_F3y_bb;
  Double I_NAI_D2x_G3yz_bb = I_NAI_F2xz_F3y_bb+ABZ*I_NAI_D2x_F3y_bb;
  Double I_NAI_Dxy_G3yz_bb = I_NAI_Fxyz_F3y_bb+ABZ*I_NAI_Dxy_F3y_bb;
  Double I_NAI_Dxz_G3yz_bb = I_NAI_Fx2z_F3y_bb+ABZ*I_NAI_Dxz_F3y_bb;
  Double I_NAI_D2y_G3yz_bb = I_NAI_F2yz_F3y_bb+ABZ*I_NAI_D2y_F3y_bb;
  Double I_NAI_Dyz_G3yz_bb = I_NAI_Fy2z_F3y_bb+ABZ*I_NAI_Dyz_F3y_bb;
  Double I_NAI_D2z_G3yz_bb = I_NAI_F3z_F3y_bb+ABZ*I_NAI_D2z_F3y_bb;
  Double I_NAI_D2x_G2y2z_bb = I_NAI_F2xz_F2yz_bb+ABZ*I_NAI_D2x_F2yz_bb;
  Double I_NAI_Dxy_G2y2z_bb = I_NAI_Fxyz_F2yz_bb+ABZ*I_NAI_Dxy_F2yz_bb;
  Double I_NAI_Dxz_G2y2z_bb = I_NAI_Fx2z_F2yz_bb+ABZ*I_NAI_Dxz_F2yz_bb;
  Double I_NAI_D2y_G2y2z_bb = I_NAI_F2yz_F2yz_bb+ABZ*I_NAI_D2y_F2yz_bb;
  Double I_NAI_Dyz_G2y2z_bb = I_NAI_Fy2z_F2yz_bb+ABZ*I_NAI_Dyz_F2yz_bb;
  Double I_NAI_D2z_G2y2z_bb = I_NAI_F3z_F2yz_bb+ABZ*I_NAI_D2z_F2yz_bb;
  Double I_NAI_D2x_Gy3z_bb = I_NAI_F2xy_F3z_bb+ABY*I_NAI_D2x_F3z_bb;
  Double I_NAI_Dxy_Gy3z_bb = I_NAI_Fx2y_F3z_bb+ABY*I_NAI_Dxy_F3z_bb;
  Double I_NAI_Dxz_Gy3z_bb = I_NAI_Fxyz_F3z_bb+ABY*I_NAI_Dxz_F3z_bb;
  Double I_NAI_D2y_Gy3z_bb = I_NAI_F3y_F3z_bb+ABY*I_NAI_D2y_F3z_bb;
  Double I_NAI_Dyz_Gy3z_bb = I_NAI_F2yz_F3z_bb+ABY*I_NAI_Dyz_F3z_bb;
  Double I_NAI_D2z_Gy3z_bb = I_NAI_Fy2z_F3z_bb+ABY*I_NAI_D2z_F3z_bb;
  Double I_NAI_D2x_G4z_bb = I_NAI_F2xz_F3z_bb+ABZ*I_NAI_D2x_F3z_bb;
  Double I_NAI_Dxy_G4z_bb = I_NAI_Fxyz_F3z_bb+ABZ*I_NAI_Dxy_F3z_bb;
  Double I_NAI_Dxz_G4z_bb = I_NAI_Fx2z_F3z_bb+ABZ*I_NAI_Dxz_F3z_bb;
  Double I_NAI_D2y_G4z_bb = I_NAI_F2yz_F3z_bb+ABZ*I_NAI_D2y_F3z_bb;
  Double I_NAI_Dyz_G4z_bb = I_NAI_Fy2z_F3z_bb+ABZ*I_NAI_Dyz_F3z_bb;
  Double I_NAI_D2z_G4z_bb = I_NAI_F3z_F3z_bb+ABZ*I_NAI_D2z_F3z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[0] = 4.0E0*I_NAI_G4x_D2x_aa-2.0E0*2*I_NAI_D2x_D2x_a-2.0E0*3*I_NAI_D2x_D2x_a+2*1*I_NAI_S_D2x;
  abcd[1] = 4.0E0*I_NAI_G3xy_D2x_aa-2.0E0*1*I_NAI_Dxy_D2x_a-2.0E0*2*I_NAI_Dxy_D2x_a;
  abcd[2] = 4.0E0*I_NAI_G3xz_D2x_aa-2.0E0*1*I_NAI_Dxz_D2x_a-2.0E0*2*I_NAI_Dxz_D2x_a;
  abcd[3] = 4.0E0*I_NAI_G2x2y_D2x_aa-2.0E0*1*I_NAI_D2y_D2x_a;
  abcd[4] = 4.0E0*I_NAI_G2xyz_D2x_aa-2.0E0*1*I_NAI_Dyz_D2x_a;
  abcd[5] = 4.0E0*I_NAI_G2x2z_D2x_aa-2.0E0*1*I_NAI_D2z_D2x_a;
  abcd[6] = 4.0E0*I_NAI_G4x_Dxy_aa-2.0E0*2*I_NAI_D2x_Dxy_a-2.0E0*3*I_NAI_D2x_Dxy_a+2*1*I_NAI_S_Dxy;
  abcd[7] = 4.0E0*I_NAI_G3xy_Dxy_aa-2.0E0*1*I_NAI_Dxy_Dxy_a-2.0E0*2*I_NAI_Dxy_Dxy_a;
  abcd[8] = 4.0E0*I_NAI_G3xz_Dxy_aa-2.0E0*1*I_NAI_Dxz_Dxy_a-2.0E0*2*I_NAI_Dxz_Dxy_a;
  abcd[9] = 4.0E0*I_NAI_G2x2y_Dxy_aa-2.0E0*1*I_NAI_D2y_Dxy_a;
  abcd[10] = 4.0E0*I_NAI_G2xyz_Dxy_aa-2.0E0*1*I_NAI_Dyz_Dxy_a;
  abcd[11] = 4.0E0*I_NAI_G2x2z_Dxy_aa-2.0E0*1*I_NAI_D2z_Dxy_a;
  abcd[12] = 4.0E0*I_NAI_G4x_Dxz_aa-2.0E0*2*I_NAI_D2x_Dxz_a-2.0E0*3*I_NAI_D2x_Dxz_a+2*1*I_NAI_S_Dxz;
  abcd[13] = 4.0E0*I_NAI_G3xy_Dxz_aa-2.0E0*1*I_NAI_Dxy_Dxz_a-2.0E0*2*I_NAI_Dxy_Dxz_a;
  abcd[14] = 4.0E0*I_NAI_G3xz_Dxz_aa-2.0E0*1*I_NAI_Dxz_Dxz_a-2.0E0*2*I_NAI_Dxz_Dxz_a;
  abcd[15] = 4.0E0*I_NAI_G2x2y_Dxz_aa-2.0E0*1*I_NAI_D2y_Dxz_a;
  abcd[16] = 4.0E0*I_NAI_G2xyz_Dxz_aa-2.0E0*1*I_NAI_Dyz_Dxz_a;
  abcd[17] = 4.0E0*I_NAI_G2x2z_Dxz_aa-2.0E0*1*I_NAI_D2z_Dxz_a;
  abcd[18] = 4.0E0*I_NAI_G4x_D2y_aa-2.0E0*2*I_NAI_D2x_D2y_a-2.0E0*3*I_NAI_D2x_D2y_a+2*1*I_NAI_S_D2y;
  abcd[19] = 4.0E0*I_NAI_G3xy_D2y_aa-2.0E0*1*I_NAI_Dxy_D2y_a-2.0E0*2*I_NAI_Dxy_D2y_a;
  abcd[20] = 4.0E0*I_NAI_G3xz_D2y_aa-2.0E0*1*I_NAI_Dxz_D2y_a-2.0E0*2*I_NAI_Dxz_D2y_a;
  abcd[21] = 4.0E0*I_NAI_G2x2y_D2y_aa-2.0E0*1*I_NAI_D2y_D2y_a;
  abcd[22] = 4.0E0*I_NAI_G2xyz_D2y_aa-2.0E0*1*I_NAI_Dyz_D2y_a;
  abcd[23] = 4.0E0*I_NAI_G2x2z_D2y_aa-2.0E0*1*I_NAI_D2z_D2y_a;
  abcd[24] = 4.0E0*I_NAI_G4x_Dyz_aa-2.0E0*2*I_NAI_D2x_Dyz_a-2.0E0*3*I_NAI_D2x_Dyz_a+2*1*I_NAI_S_Dyz;
  abcd[25] = 4.0E0*I_NAI_G3xy_Dyz_aa-2.0E0*1*I_NAI_Dxy_Dyz_a-2.0E0*2*I_NAI_Dxy_Dyz_a;
  abcd[26] = 4.0E0*I_NAI_G3xz_Dyz_aa-2.0E0*1*I_NAI_Dxz_Dyz_a-2.0E0*2*I_NAI_Dxz_Dyz_a;
  abcd[27] = 4.0E0*I_NAI_G2x2y_Dyz_aa-2.0E0*1*I_NAI_D2y_Dyz_a;
  abcd[28] = 4.0E0*I_NAI_G2xyz_Dyz_aa-2.0E0*1*I_NAI_Dyz_Dyz_a;
  abcd[29] = 4.0E0*I_NAI_G2x2z_Dyz_aa-2.0E0*1*I_NAI_D2z_Dyz_a;
  abcd[30] = 4.0E0*I_NAI_G4x_D2z_aa-2.0E0*2*I_NAI_D2x_D2z_a-2.0E0*3*I_NAI_D2x_D2z_a+2*1*I_NAI_S_D2z;
  abcd[31] = 4.0E0*I_NAI_G3xy_D2z_aa-2.0E0*1*I_NAI_Dxy_D2z_a-2.0E0*2*I_NAI_Dxy_D2z_a;
  abcd[32] = 4.0E0*I_NAI_G3xz_D2z_aa-2.0E0*1*I_NAI_Dxz_D2z_a-2.0E0*2*I_NAI_Dxz_D2z_a;
  abcd[33] = 4.0E0*I_NAI_G2x2y_D2z_aa-2.0E0*1*I_NAI_D2y_D2z_a;
  abcd[34] = 4.0E0*I_NAI_G2xyz_D2z_aa-2.0E0*1*I_NAI_Dyz_D2z_a;
  abcd[35] = 4.0E0*I_NAI_G2x2z_D2z_aa-2.0E0*1*I_NAI_D2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[36] = 4.0E0*I_NAI_G3xy_D2x_aa-2.0E0*2*I_NAI_Dxy_D2x_a;
  abcd[37] = 4.0E0*I_NAI_G2x2y_D2x_aa-2.0E0*1*I_NAI_D2x_D2x_a-2.0E0*1*I_NAI_D2y_D2x_a+1*I_NAI_S_D2x;
  abcd[38] = 4.0E0*I_NAI_G2xyz_D2x_aa-2.0E0*1*I_NAI_Dyz_D2x_a;
  abcd[39] = 4.0E0*I_NAI_Gx3y_D2x_aa-2.0E0*2*I_NAI_Dxy_D2x_a;
  abcd[40] = 4.0E0*I_NAI_Gx2yz_D2x_aa-2.0E0*1*I_NAI_Dxz_D2x_a;
  abcd[41] = 4.0E0*I_NAI_Gxy2z_D2x_aa;
  abcd[42] = 4.0E0*I_NAI_G3xy_Dxy_aa-2.0E0*2*I_NAI_Dxy_Dxy_a;
  abcd[43] = 4.0E0*I_NAI_G2x2y_Dxy_aa-2.0E0*1*I_NAI_D2x_Dxy_a-2.0E0*1*I_NAI_D2y_Dxy_a+1*I_NAI_S_Dxy;
  abcd[44] = 4.0E0*I_NAI_G2xyz_Dxy_aa-2.0E0*1*I_NAI_Dyz_Dxy_a;
  abcd[45] = 4.0E0*I_NAI_Gx3y_Dxy_aa-2.0E0*2*I_NAI_Dxy_Dxy_a;
  abcd[46] = 4.0E0*I_NAI_Gx2yz_Dxy_aa-2.0E0*1*I_NAI_Dxz_Dxy_a;
  abcd[47] = 4.0E0*I_NAI_Gxy2z_Dxy_aa;
  abcd[48] = 4.0E0*I_NAI_G3xy_Dxz_aa-2.0E0*2*I_NAI_Dxy_Dxz_a;
  abcd[49] = 4.0E0*I_NAI_G2x2y_Dxz_aa-2.0E0*1*I_NAI_D2x_Dxz_a-2.0E0*1*I_NAI_D2y_Dxz_a+1*I_NAI_S_Dxz;
  abcd[50] = 4.0E0*I_NAI_G2xyz_Dxz_aa-2.0E0*1*I_NAI_Dyz_Dxz_a;
  abcd[51] = 4.0E0*I_NAI_Gx3y_Dxz_aa-2.0E0*2*I_NAI_Dxy_Dxz_a;
  abcd[52] = 4.0E0*I_NAI_Gx2yz_Dxz_aa-2.0E0*1*I_NAI_Dxz_Dxz_a;
  abcd[53] = 4.0E0*I_NAI_Gxy2z_Dxz_aa;
  abcd[54] = 4.0E0*I_NAI_G3xy_D2y_aa-2.0E0*2*I_NAI_Dxy_D2y_a;
  abcd[55] = 4.0E0*I_NAI_G2x2y_D2y_aa-2.0E0*1*I_NAI_D2x_D2y_a-2.0E0*1*I_NAI_D2y_D2y_a+1*I_NAI_S_D2y;
  abcd[56] = 4.0E0*I_NAI_G2xyz_D2y_aa-2.0E0*1*I_NAI_Dyz_D2y_a;
  abcd[57] = 4.0E0*I_NAI_Gx3y_D2y_aa-2.0E0*2*I_NAI_Dxy_D2y_a;
  abcd[58] = 4.0E0*I_NAI_Gx2yz_D2y_aa-2.0E0*1*I_NAI_Dxz_D2y_a;
  abcd[59] = 4.0E0*I_NAI_Gxy2z_D2y_aa;
  abcd[60] = 4.0E0*I_NAI_G3xy_Dyz_aa-2.0E0*2*I_NAI_Dxy_Dyz_a;
  abcd[61] = 4.0E0*I_NAI_G2x2y_Dyz_aa-2.0E0*1*I_NAI_D2x_Dyz_a-2.0E0*1*I_NAI_D2y_Dyz_a+1*I_NAI_S_Dyz;
  abcd[62] = 4.0E0*I_NAI_G2xyz_Dyz_aa-2.0E0*1*I_NAI_Dyz_Dyz_a;
  abcd[63] = 4.0E0*I_NAI_Gx3y_Dyz_aa-2.0E0*2*I_NAI_Dxy_Dyz_a;
  abcd[64] = 4.0E0*I_NAI_Gx2yz_Dyz_aa-2.0E0*1*I_NAI_Dxz_Dyz_a;
  abcd[65] = 4.0E0*I_NAI_Gxy2z_Dyz_aa;
  abcd[66] = 4.0E0*I_NAI_G3xy_D2z_aa-2.0E0*2*I_NAI_Dxy_D2z_a;
  abcd[67] = 4.0E0*I_NAI_G2x2y_D2z_aa-2.0E0*1*I_NAI_D2x_D2z_a-2.0E0*1*I_NAI_D2y_D2z_a+1*I_NAI_S_D2z;
  abcd[68] = 4.0E0*I_NAI_G2xyz_D2z_aa-2.0E0*1*I_NAI_Dyz_D2z_a;
  abcd[69] = 4.0E0*I_NAI_Gx3y_D2z_aa-2.0E0*2*I_NAI_Dxy_D2z_a;
  abcd[70] = 4.0E0*I_NAI_Gx2yz_D2z_aa-2.0E0*1*I_NAI_Dxz_D2z_a;
  abcd[71] = 4.0E0*I_NAI_Gxy2z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[72] = 4.0E0*I_NAI_G3xz_D2x_aa-2.0E0*2*I_NAI_Dxz_D2x_a;
  abcd[73] = 4.0E0*I_NAI_G2xyz_D2x_aa-2.0E0*1*I_NAI_Dyz_D2x_a;
  abcd[74] = 4.0E0*I_NAI_G2x2z_D2x_aa-2.0E0*1*I_NAI_D2x_D2x_a-2.0E0*1*I_NAI_D2z_D2x_a+1*I_NAI_S_D2x;
  abcd[75] = 4.0E0*I_NAI_Gx2yz_D2x_aa;
  abcd[76] = 4.0E0*I_NAI_Gxy2z_D2x_aa-2.0E0*1*I_NAI_Dxy_D2x_a;
  abcd[77] = 4.0E0*I_NAI_Gx3z_D2x_aa-2.0E0*2*I_NAI_Dxz_D2x_a;
  abcd[78] = 4.0E0*I_NAI_G3xz_Dxy_aa-2.0E0*2*I_NAI_Dxz_Dxy_a;
  abcd[79] = 4.0E0*I_NAI_G2xyz_Dxy_aa-2.0E0*1*I_NAI_Dyz_Dxy_a;
  abcd[80] = 4.0E0*I_NAI_G2x2z_Dxy_aa-2.0E0*1*I_NAI_D2x_Dxy_a-2.0E0*1*I_NAI_D2z_Dxy_a+1*I_NAI_S_Dxy;
  abcd[81] = 4.0E0*I_NAI_Gx2yz_Dxy_aa;
  abcd[82] = 4.0E0*I_NAI_Gxy2z_Dxy_aa-2.0E0*1*I_NAI_Dxy_Dxy_a;
  abcd[83] = 4.0E0*I_NAI_Gx3z_Dxy_aa-2.0E0*2*I_NAI_Dxz_Dxy_a;
  abcd[84] = 4.0E0*I_NAI_G3xz_Dxz_aa-2.0E0*2*I_NAI_Dxz_Dxz_a;
  abcd[85] = 4.0E0*I_NAI_G2xyz_Dxz_aa-2.0E0*1*I_NAI_Dyz_Dxz_a;
  abcd[86] = 4.0E0*I_NAI_G2x2z_Dxz_aa-2.0E0*1*I_NAI_D2x_Dxz_a-2.0E0*1*I_NAI_D2z_Dxz_a+1*I_NAI_S_Dxz;
  abcd[87] = 4.0E0*I_NAI_Gx2yz_Dxz_aa;
  abcd[88] = 4.0E0*I_NAI_Gxy2z_Dxz_aa-2.0E0*1*I_NAI_Dxy_Dxz_a;
  abcd[89] = 4.0E0*I_NAI_Gx3z_Dxz_aa-2.0E0*2*I_NAI_Dxz_Dxz_a;
  abcd[90] = 4.0E0*I_NAI_G3xz_D2y_aa-2.0E0*2*I_NAI_Dxz_D2y_a;
  abcd[91] = 4.0E0*I_NAI_G2xyz_D2y_aa-2.0E0*1*I_NAI_Dyz_D2y_a;
  abcd[92] = 4.0E0*I_NAI_G2x2z_D2y_aa-2.0E0*1*I_NAI_D2x_D2y_a-2.0E0*1*I_NAI_D2z_D2y_a+1*I_NAI_S_D2y;
  abcd[93] = 4.0E0*I_NAI_Gx2yz_D2y_aa;
  abcd[94] = 4.0E0*I_NAI_Gxy2z_D2y_aa-2.0E0*1*I_NAI_Dxy_D2y_a;
  abcd[95] = 4.0E0*I_NAI_Gx3z_D2y_aa-2.0E0*2*I_NAI_Dxz_D2y_a;
  abcd[96] = 4.0E0*I_NAI_G3xz_Dyz_aa-2.0E0*2*I_NAI_Dxz_Dyz_a;
  abcd[97] = 4.0E0*I_NAI_G2xyz_Dyz_aa-2.0E0*1*I_NAI_Dyz_Dyz_a;
  abcd[98] = 4.0E0*I_NAI_G2x2z_Dyz_aa-2.0E0*1*I_NAI_D2x_Dyz_a-2.0E0*1*I_NAI_D2z_Dyz_a+1*I_NAI_S_Dyz;
  abcd[99] = 4.0E0*I_NAI_Gx2yz_Dyz_aa;
  abcd[100] = 4.0E0*I_NAI_Gxy2z_Dyz_aa-2.0E0*1*I_NAI_Dxy_Dyz_a;
  abcd[101] = 4.0E0*I_NAI_Gx3z_Dyz_aa-2.0E0*2*I_NAI_Dxz_Dyz_a;
  abcd[102] = 4.0E0*I_NAI_G3xz_D2z_aa-2.0E0*2*I_NAI_Dxz_D2z_a;
  abcd[103] = 4.0E0*I_NAI_G2xyz_D2z_aa-2.0E0*1*I_NAI_Dyz_D2z_a;
  abcd[104] = 4.0E0*I_NAI_G2x2z_D2z_aa-2.0E0*1*I_NAI_D2x_D2z_a-2.0E0*1*I_NAI_D2z_D2z_a+1*I_NAI_S_D2z;
  abcd[105] = 4.0E0*I_NAI_Gx2yz_D2z_aa;
  abcd[106] = 4.0E0*I_NAI_Gxy2z_D2z_aa-2.0E0*1*I_NAI_Dxy_D2z_a;
  abcd[107] = 4.0E0*I_NAI_Gx3z_D2z_aa-2.0E0*2*I_NAI_Dxz_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[108] = 4.0E0*I_NAI_G2x2y_D2x_aa-2.0E0*1*I_NAI_D2x_D2x_a;
  abcd[109] = 4.0E0*I_NAI_Gx3y_D2x_aa-2.0E0*1*I_NAI_Dxy_D2x_a-2.0E0*2*I_NAI_Dxy_D2x_a;
  abcd[110] = 4.0E0*I_NAI_Gx2yz_D2x_aa-2.0E0*1*I_NAI_Dxz_D2x_a;
  abcd[111] = 4.0E0*I_NAI_G4y_D2x_aa-2.0E0*2*I_NAI_D2y_D2x_a-2.0E0*3*I_NAI_D2y_D2x_a+2*1*I_NAI_S_D2x;
  abcd[112] = 4.0E0*I_NAI_G3yz_D2x_aa-2.0E0*1*I_NAI_Dyz_D2x_a-2.0E0*2*I_NAI_Dyz_D2x_a;
  abcd[113] = 4.0E0*I_NAI_G2y2z_D2x_aa-2.0E0*1*I_NAI_D2z_D2x_a;
  abcd[114] = 4.0E0*I_NAI_G2x2y_Dxy_aa-2.0E0*1*I_NAI_D2x_Dxy_a;
  abcd[115] = 4.0E0*I_NAI_Gx3y_Dxy_aa-2.0E0*1*I_NAI_Dxy_Dxy_a-2.0E0*2*I_NAI_Dxy_Dxy_a;
  abcd[116] = 4.0E0*I_NAI_Gx2yz_Dxy_aa-2.0E0*1*I_NAI_Dxz_Dxy_a;
  abcd[117] = 4.0E0*I_NAI_G4y_Dxy_aa-2.0E0*2*I_NAI_D2y_Dxy_a-2.0E0*3*I_NAI_D2y_Dxy_a+2*1*I_NAI_S_Dxy;
  abcd[118] = 4.0E0*I_NAI_G3yz_Dxy_aa-2.0E0*1*I_NAI_Dyz_Dxy_a-2.0E0*2*I_NAI_Dyz_Dxy_a;
  abcd[119] = 4.0E0*I_NAI_G2y2z_Dxy_aa-2.0E0*1*I_NAI_D2z_Dxy_a;
  abcd[120] = 4.0E0*I_NAI_G2x2y_Dxz_aa-2.0E0*1*I_NAI_D2x_Dxz_a;
  abcd[121] = 4.0E0*I_NAI_Gx3y_Dxz_aa-2.0E0*1*I_NAI_Dxy_Dxz_a-2.0E0*2*I_NAI_Dxy_Dxz_a;
  abcd[122] = 4.0E0*I_NAI_Gx2yz_Dxz_aa-2.0E0*1*I_NAI_Dxz_Dxz_a;
  abcd[123] = 4.0E0*I_NAI_G4y_Dxz_aa-2.0E0*2*I_NAI_D2y_Dxz_a-2.0E0*3*I_NAI_D2y_Dxz_a+2*1*I_NAI_S_Dxz;
  abcd[124] = 4.0E0*I_NAI_G3yz_Dxz_aa-2.0E0*1*I_NAI_Dyz_Dxz_a-2.0E0*2*I_NAI_Dyz_Dxz_a;
  abcd[125] = 4.0E0*I_NAI_G2y2z_Dxz_aa-2.0E0*1*I_NAI_D2z_Dxz_a;
  abcd[126] = 4.0E0*I_NAI_G2x2y_D2y_aa-2.0E0*1*I_NAI_D2x_D2y_a;
  abcd[127] = 4.0E0*I_NAI_Gx3y_D2y_aa-2.0E0*1*I_NAI_Dxy_D2y_a-2.0E0*2*I_NAI_Dxy_D2y_a;
  abcd[128] = 4.0E0*I_NAI_Gx2yz_D2y_aa-2.0E0*1*I_NAI_Dxz_D2y_a;
  abcd[129] = 4.0E0*I_NAI_G4y_D2y_aa-2.0E0*2*I_NAI_D2y_D2y_a-2.0E0*3*I_NAI_D2y_D2y_a+2*1*I_NAI_S_D2y;
  abcd[130] = 4.0E0*I_NAI_G3yz_D2y_aa-2.0E0*1*I_NAI_Dyz_D2y_a-2.0E0*2*I_NAI_Dyz_D2y_a;
  abcd[131] = 4.0E0*I_NAI_G2y2z_D2y_aa-2.0E0*1*I_NAI_D2z_D2y_a;
  abcd[132] = 4.0E0*I_NAI_G2x2y_Dyz_aa-2.0E0*1*I_NAI_D2x_Dyz_a;
  abcd[133] = 4.0E0*I_NAI_Gx3y_Dyz_aa-2.0E0*1*I_NAI_Dxy_Dyz_a-2.0E0*2*I_NAI_Dxy_Dyz_a;
  abcd[134] = 4.0E0*I_NAI_Gx2yz_Dyz_aa-2.0E0*1*I_NAI_Dxz_Dyz_a;
  abcd[135] = 4.0E0*I_NAI_G4y_Dyz_aa-2.0E0*2*I_NAI_D2y_Dyz_a-2.0E0*3*I_NAI_D2y_Dyz_a+2*1*I_NAI_S_Dyz;
  abcd[136] = 4.0E0*I_NAI_G3yz_Dyz_aa-2.0E0*1*I_NAI_Dyz_Dyz_a-2.0E0*2*I_NAI_Dyz_Dyz_a;
  abcd[137] = 4.0E0*I_NAI_G2y2z_Dyz_aa-2.0E0*1*I_NAI_D2z_Dyz_a;
  abcd[138] = 4.0E0*I_NAI_G2x2y_D2z_aa-2.0E0*1*I_NAI_D2x_D2z_a;
  abcd[139] = 4.0E0*I_NAI_Gx3y_D2z_aa-2.0E0*1*I_NAI_Dxy_D2z_a-2.0E0*2*I_NAI_Dxy_D2z_a;
  abcd[140] = 4.0E0*I_NAI_Gx2yz_D2z_aa-2.0E0*1*I_NAI_Dxz_D2z_a;
  abcd[141] = 4.0E0*I_NAI_G4y_D2z_aa-2.0E0*2*I_NAI_D2y_D2z_a-2.0E0*3*I_NAI_D2y_D2z_a+2*1*I_NAI_S_D2z;
  abcd[142] = 4.0E0*I_NAI_G3yz_D2z_aa-2.0E0*1*I_NAI_Dyz_D2z_a-2.0E0*2*I_NAI_Dyz_D2z_a;
  abcd[143] = 4.0E0*I_NAI_G2y2z_D2z_aa-2.0E0*1*I_NAI_D2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[144] = 4.0E0*I_NAI_G2xyz_D2x_aa;
  abcd[145] = 4.0E0*I_NAI_Gx2yz_D2x_aa-2.0E0*1*I_NAI_Dxz_D2x_a;
  abcd[146] = 4.0E0*I_NAI_Gxy2z_D2x_aa-2.0E0*1*I_NAI_Dxy_D2x_a;
  abcd[147] = 4.0E0*I_NAI_G3yz_D2x_aa-2.0E0*2*I_NAI_Dyz_D2x_a;
  abcd[148] = 4.0E0*I_NAI_G2y2z_D2x_aa-2.0E0*1*I_NAI_D2y_D2x_a-2.0E0*1*I_NAI_D2z_D2x_a+1*I_NAI_S_D2x;
  abcd[149] = 4.0E0*I_NAI_Gy3z_D2x_aa-2.0E0*2*I_NAI_Dyz_D2x_a;
  abcd[150] = 4.0E0*I_NAI_G2xyz_Dxy_aa;
  abcd[151] = 4.0E0*I_NAI_Gx2yz_Dxy_aa-2.0E0*1*I_NAI_Dxz_Dxy_a;
  abcd[152] = 4.0E0*I_NAI_Gxy2z_Dxy_aa-2.0E0*1*I_NAI_Dxy_Dxy_a;
  abcd[153] = 4.0E0*I_NAI_G3yz_Dxy_aa-2.0E0*2*I_NAI_Dyz_Dxy_a;
  abcd[154] = 4.0E0*I_NAI_G2y2z_Dxy_aa-2.0E0*1*I_NAI_D2y_Dxy_a-2.0E0*1*I_NAI_D2z_Dxy_a+1*I_NAI_S_Dxy;
  abcd[155] = 4.0E0*I_NAI_Gy3z_Dxy_aa-2.0E0*2*I_NAI_Dyz_Dxy_a;
  abcd[156] = 4.0E0*I_NAI_G2xyz_Dxz_aa;
  abcd[157] = 4.0E0*I_NAI_Gx2yz_Dxz_aa-2.0E0*1*I_NAI_Dxz_Dxz_a;
  abcd[158] = 4.0E0*I_NAI_Gxy2z_Dxz_aa-2.0E0*1*I_NAI_Dxy_Dxz_a;
  abcd[159] = 4.0E0*I_NAI_G3yz_Dxz_aa-2.0E0*2*I_NAI_Dyz_Dxz_a;
  abcd[160] = 4.0E0*I_NAI_G2y2z_Dxz_aa-2.0E0*1*I_NAI_D2y_Dxz_a-2.0E0*1*I_NAI_D2z_Dxz_a+1*I_NAI_S_Dxz;
  abcd[161] = 4.0E0*I_NAI_Gy3z_Dxz_aa-2.0E0*2*I_NAI_Dyz_Dxz_a;
  abcd[162] = 4.0E0*I_NAI_G2xyz_D2y_aa;
  abcd[163] = 4.0E0*I_NAI_Gx2yz_D2y_aa-2.0E0*1*I_NAI_Dxz_D2y_a;
  abcd[164] = 4.0E0*I_NAI_Gxy2z_D2y_aa-2.0E0*1*I_NAI_Dxy_D2y_a;
  abcd[165] = 4.0E0*I_NAI_G3yz_D2y_aa-2.0E0*2*I_NAI_Dyz_D2y_a;
  abcd[166] = 4.0E0*I_NAI_G2y2z_D2y_aa-2.0E0*1*I_NAI_D2y_D2y_a-2.0E0*1*I_NAI_D2z_D2y_a+1*I_NAI_S_D2y;
  abcd[167] = 4.0E0*I_NAI_Gy3z_D2y_aa-2.0E0*2*I_NAI_Dyz_D2y_a;
  abcd[168] = 4.0E0*I_NAI_G2xyz_Dyz_aa;
  abcd[169] = 4.0E0*I_NAI_Gx2yz_Dyz_aa-2.0E0*1*I_NAI_Dxz_Dyz_a;
  abcd[170] = 4.0E0*I_NAI_Gxy2z_Dyz_aa-2.0E0*1*I_NAI_Dxy_Dyz_a;
  abcd[171] = 4.0E0*I_NAI_G3yz_Dyz_aa-2.0E0*2*I_NAI_Dyz_Dyz_a;
  abcd[172] = 4.0E0*I_NAI_G2y2z_Dyz_aa-2.0E0*1*I_NAI_D2y_Dyz_a-2.0E0*1*I_NAI_D2z_Dyz_a+1*I_NAI_S_Dyz;
  abcd[173] = 4.0E0*I_NAI_Gy3z_Dyz_aa-2.0E0*2*I_NAI_Dyz_Dyz_a;
  abcd[174] = 4.0E0*I_NAI_G2xyz_D2z_aa;
  abcd[175] = 4.0E0*I_NAI_Gx2yz_D2z_aa-2.0E0*1*I_NAI_Dxz_D2z_a;
  abcd[176] = 4.0E0*I_NAI_Gxy2z_D2z_aa-2.0E0*1*I_NAI_Dxy_D2z_a;
  abcd[177] = 4.0E0*I_NAI_G3yz_D2z_aa-2.0E0*2*I_NAI_Dyz_D2z_a;
  abcd[178] = 4.0E0*I_NAI_G2y2z_D2z_aa-2.0E0*1*I_NAI_D2y_D2z_a-2.0E0*1*I_NAI_D2z_D2z_a+1*I_NAI_S_D2z;
  abcd[179] = 4.0E0*I_NAI_Gy3z_D2z_aa-2.0E0*2*I_NAI_Dyz_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_aa
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_D_D_a
   * RHS shell quartet name: SQ_NAI_S_D
   ************************************************************/
  abcd[180] = 4.0E0*I_NAI_G2x2z_D2x_aa-2.0E0*1*I_NAI_D2x_D2x_a;
  abcd[181] = 4.0E0*I_NAI_Gxy2z_D2x_aa-2.0E0*1*I_NAI_Dxy_D2x_a;
  abcd[182] = 4.0E0*I_NAI_Gx3z_D2x_aa-2.0E0*1*I_NAI_Dxz_D2x_a-2.0E0*2*I_NAI_Dxz_D2x_a;
  abcd[183] = 4.0E0*I_NAI_G2y2z_D2x_aa-2.0E0*1*I_NAI_D2y_D2x_a;
  abcd[184] = 4.0E0*I_NAI_Gy3z_D2x_aa-2.0E0*1*I_NAI_Dyz_D2x_a-2.0E0*2*I_NAI_Dyz_D2x_a;
  abcd[185] = 4.0E0*I_NAI_G4z_D2x_aa-2.0E0*2*I_NAI_D2z_D2x_a-2.0E0*3*I_NAI_D2z_D2x_a+2*1*I_NAI_S_D2x;
  abcd[186] = 4.0E0*I_NAI_G2x2z_Dxy_aa-2.0E0*1*I_NAI_D2x_Dxy_a;
  abcd[187] = 4.0E0*I_NAI_Gxy2z_Dxy_aa-2.0E0*1*I_NAI_Dxy_Dxy_a;
  abcd[188] = 4.0E0*I_NAI_Gx3z_Dxy_aa-2.0E0*1*I_NAI_Dxz_Dxy_a-2.0E0*2*I_NAI_Dxz_Dxy_a;
  abcd[189] = 4.0E0*I_NAI_G2y2z_Dxy_aa-2.0E0*1*I_NAI_D2y_Dxy_a;
  abcd[190] = 4.0E0*I_NAI_Gy3z_Dxy_aa-2.0E0*1*I_NAI_Dyz_Dxy_a-2.0E0*2*I_NAI_Dyz_Dxy_a;
  abcd[191] = 4.0E0*I_NAI_G4z_Dxy_aa-2.0E0*2*I_NAI_D2z_Dxy_a-2.0E0*3*I_NAI_D2z_Dxy_a+2*1*I_NAI_S_Dxy;
  abcd[192] = 4.0E0*I_NAI_G2x2z_Dxz_aa-2.0E0*1*I_NAI_D2x_Dxz_a;
  abcd[193] = 4.0E0*I_NAI_Gxy2z_Dxz_aa-2.0E0*1*I_NAI_Dxy_Dxz_a;
  abcd[194] = 4.0E0*I_NAI_Gx3z_Dxz_aa-2.0E0*1*I_NAI_Dxz_Dxz_a-2.0E0*2*I_NAI_Dxz_Dxz_a;
  abcd[195] = 4.0E0*I_NAI_G2y2z_Dxz_aa-2.0E0*1*I_NAI_D2y_Dxz_a;
  abcd[196] = 4.0E0*I_NAI_Gy3z_Dxz_aa-2.0E0*1*I_NAI_Dyz_Dxz_a-2.0E0*2*I_NAI_Dyz_Dxz_a;
  abcd[197] = 4.0E0*I_NAI_G4z_Dxz_aa-2.0E0*2*I_NAI_D2z_Dxz_a-2.0E0*3*I_NAI_D2z_Dxz_a+2*1*I_NAI_S_Dxz;
  abcd[198] = 4.0E0*I_NAI_G2x2z_D2y_aa-2.0E0*1*I_NAI_D2x_D2y_a;
  abcd[199] = 4.0E0*I_NAI_Gxy2z_D2y_aa-2.0E0*1*I_NAI_Dxy_D2y_a;
  abcd[200] = 4.0E0*I_NAI_Gx3z_D2y_aa-2.0E0*1*I_NAI_Dxz_D2y_a-2.0E0*2*I_NAI_Dxz_D2y_a;
  abcd[201] = 4.0E0*I_NAI_G2y2z_D2y_aa-2.0E0*1*I_NAI_D2y_D2y_a;
  abcd[202] = 4.0E0*I_NAI_Gy3z_D2y_aa-2.0E0*1*I_NAI_Dyz_D2y_a-2.0E0*2*I_NAI_Dyz_D2y_a;
  abcd[203] = 4.0E0*I_NAI_G4z_D2y_aa-2.0E0*2*I_NAI_D2z_D2y_a-2.0E0*3*I_NAI_D2z_D2y_a+2*1*I_NAI_S_D2y;
  abcd[204] = 4.0E0*I_NAI_G2x2z_Dyz_aa-2.0E0*1*I_NAI_D2x_Dyz_a;
  abcd[205] = 4.0E0*I_NAI_Gxy2z_Dyz_aa-2.0E0*1*I_NAI_Dxy_Dyz_a;
  abcd[206] = 4.0E0*I_NAI_Gx3z_Dyz_aa-2.0E0*1*I_NAI_Dxz_Dyz_a-2.0E0*2*I_NAI_Dxz_Dyz_a;
  abcd[207] = 4.0E0*I_NAI_G2y2z_Dyz_aa-2.0E0*1*I_NAI_D2y_Dyz_a;
  abcd[208] = 4.0E0*I_NAI_Gy3z_Dyz_aa-2.0E0*1*I_NAI_Dyz_Dyz_a-2.0E0*2*I_NAI_Dyz_Dyz_a;
  abcd[209] = 4.0E0*I_NAI_G4z_Dyz_aa-2.0E0*2*I_NAI_D2z_Dyz_a-2.0E0*3*I_NAI_D2z_Dyz_a+2*1*I_NAI_S_Dyz;
  abcd[210] = 4.0E0*I_NAI_G2x2z_D2z_aa-2.0E0*1*I_NAI_D2x_D2z_a;
  abcd[211] = 4.0E0*I_NAI_Gxy2z_D2z_aa-2.0E0*1*I_NAI_Dxy_D2z_a;
  abcd[212] = 4.0E0*I_NAI_Gx3z_D2z_aa-2.0E0*1*I_NAI_Dxz_D2z_a-2.0E0*2*I_NAI_Dxz_D2z_a;
  abcd[213] = 4.0E0*I_NAI_G2y2z_D2z_aa-2.0E0*1*I_NAI_D2y_D2z_a;
  abcd[214] = 4.0E0*I_NAI_Gy3z_D2z_aa-2.0E0*1*I_NAI_Dyz_D2z_a-2.0E0*2*I_NAI_Dyz_D2z_a;
  abcd[215] = 4.0E0*I_NAI_G4z_D2z_aa-2.0E0*2*I_NAI_D2z_D2z_a-2.0E0*3*I_NAI_D2z_D2z_a+2*1*I_NAI_S_D2z;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[216] = 4.0E0*I_NAI_F3x_F3x_ab-2.0E0*2*I_NAI_F3x_Px_a-2.0E0*2*I_NAI_Px_F3x_b+2*2*I_NAI_Px_Px;
  abcd[217] = 4.0E0*I_NAI_F2xy_F3x_ab-2.0E0*2*I_NAI_F2xy_Px_a-2.0E0*1*I_NAI_Py_F3x_b+2*I_NAI_Py_Px;
  abcd[218] = 4.0E0*I_NAI_F2xz_F3x_ab-2.0E0*2*I_NAI_F2xz_Px_a-2.0E0*1*I_NAI_Pz_F3x_b+2*I_NAI_Pz_Px;
  abcd[219] = 4.0E0*I_NAI_Fx2y_F3x_ab-2.0E0*2*I_NAI_Fx2y_Px_a;
  abcd[220] = 4.0E0*I_NAI_Fxyz_F3x_ab-2.0E0*2*I_NAI_Fxyz_Px_a;
  abcd[221] = 4.0E0*I_NAI_Fx2z_F3x_ab-2.0E0*2*I_NAI_Fx2z_Px_a;
  abcd[222] = 4.0E0*I_NAI_F3x_F2xy_ab-2.0E0*1*I_NAI_F3x_Py_a-2.0E0*2*I_NAI_Px_F2xy_b+2*1*I_NAI_Px_Py;
  abcd[223] = 4.0E0*I_NAI_F2xy_F2xy_ab-2.0E0*1*I_NAI_F2xy_Py_a-2.0E0*1*I_NAI_Py_F2xy_b+1*I_NAI_Py_Py;
  abcd[224] = 4.0E0*I_NAI_F2xz_F2xy_ab-2.0E0*1*I_NAI_F2xz_Py_a-2.0E0*1*I_NAI_Pz_F2xy_b+1*I_NAI_Pz_Py;
  abcd[225] = 4.0E0*I_NAI_Fx2y_F2xy_ab-2.0E0*1*I_NAI_Fx2y_Py_a;
  abcd[226] = 4.0E0*I_NAI_Fxyz_F2xy_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[227] = 4.0E0*I_NAI_Fx2z_F2xy_ab-2.0E0*1*I_NAI_Fx2z_Py_a;
  abcd[228] = 4.0E0*I_NAI_F3x_F2xz_ab-2.0E0*1*I_NAI_F3x_Pz_a-2.0E0*2*I_NAI_Px_F2xz_b+2*1*I_NAI_Px_Pz;
  abcd[229] = 4.0E0*I_NAI_F2xy_F2xz_ab-2.0E0*1*I_NAI_F2xy_Pz_a-2.0E0*1*I_NAI_Py_F2xz_b+1*I_NAI_Py_Pz;
  abcd[230] = 4.0E0*I_NAI_F2xz_F2xz_ab-2.0E0*1*I_NAI_F2xz_Pz_a-2.0E0*1*I_NAI_Pz_F2xz_b+1*I_NAI_Pz_Pz;
  abcd[231] = 4.0E0*I_NAI_Fx2y_F2xz_ab-2.0E0*1*I_NAI_Fx2y_Pz_a;
  abcd[232] = 4.0E0*I_NAI_Fxyz_F2xz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[233] = 4.0E0*I_NAI_Fx2z_F2xz_ab-2.0E0*1*I_NAI_Fx2z_Pz_a;
  abcd[234] = 4.0E0*I_NAI_F3x_Fx2y_ab-2.0E0*2*I_NAI_Px_Fx2y_b;
  abcd[235] = 4.0E0*I_NAI_F2xy_Fx2y_ab-2.0E0*1*I_NAI_Py_Fx2y_b;
  abcd[236] = 4.0E0*I_NAI_F2xz_Fx2y_ab-2.0E0*1*I_NAI_Pz_Fx2y_b;
  abcd[237] = 4.0E0*I_NAI_Fx2y_Fx2y_ab;
  abcd[238] = 4.0E0*I_NAI_Fxyz_Fx2y_ab;
  abcd[239] = 4.0E0*I_NAI_Fx2z_Fx2y_ab;
  abcd[240] = 4.0E0*I_NAI_F3x_Fxyz_ab-2.0E0*2*I_NAI_Px_Fxyz_b;
  abcd[241] = 4.0E0*I_NAI_F2xy_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[242] = 4.0E0*I_NAI_F2xz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[243] = 4.0E0*I_NAI_Fx2y_Fxyz_ab;
  abcd[244] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[245] = 4.0E0*I_NAI_Fx2z_Fxyz_ab;
  abcd[246] = 4.0E0*I_NAI_F3x_Fx2z_ab-2.0E0*2*I_NAI_Px_Fx2z_b;
  abcd[247] = 4.0E0*I_NAI_F2xy_Fx2z_ab-2.0E0*1*I_NAI_Py_Fx2z_b;
  abcd[248] = 4.0E0*I_NAI_F2xz_Fx2z_ab-2.0E0*1*I_NAI_Pz_Fx2z_b;
  abcd[249] = 4.0E0*I_NAI_Fx2y_Fx2z_ab;
  abcd[250] = 4.0E0*I_NAI_Fxyz_Fx2z_ab;
  abcd[251] = 4.0E0*I_NAI_Fx2z_Fx2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[252] = 4.0E0*I_NAI_F3x_F2xy_ab-2.0E0*2*I_NAI_Px_F2xy_b;
  abcd[253] = 4.0E0*I_NAI_F2xy_F2xy_ab-2.0E0*1*I_NAI_Py_F2xy_b;
  abcd[254] = 4.0E0*I_NAI_F2xz_F2xy_ab-2.0E0*1*I_NAI_Pz_F2xy_b;
  abcd[255] = 4.0E0*I_NAI_Fx2y_F2xy_ab;
  abcd[256] = 4.0E0*I_NAI_Fxyz_F2xy_ab;
  abcd[257] = 4.0E0*I_NAI_Fx2z_F2xy_ab;
  abcd[258] = 4.0E0*I_NAI_F3x_Fx2y_ab-2.0E0*1*I_NAI_F3x_Px_a-2.0E0*2*I_NAI_Px_Fx2y_b+2*1*I_NAI_Px_Px;
  abcd[259] = 4.0E0*I_NAI_F2xy_Fx2y_ab-2.0E0*1*I_NAI_F2xy_Px_a-2.0E0*1*I_NAI_Py_Fx2y_b+1*I_NAI_Py_Px;
  abcd[260] = 4.0E0*I_NAI_F2xz_Fx2y_ab-2.0E0*1*I_NAI_F2xz_Px_a-2.0E0*1*I_NAI_Pz_Fx2y_b+1*I_NAI_Pz_Px;
  abcd[261] = 4.0E0*I_NAI_Fx2y_Fx2y_ab-2.0E0*1*I_NAI_Fx2y_Px_a;
  abcd[262] = 4.0E0*I_NAI_Fxyz_Fx2y_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[263] = 4.0E0*I_NAI_Fx2z_Fx2y_ab-2.0E0*1*I_NAI_Fx2z_Px_a;
  abcd[264] = 4.0E0*I_NAI_F3x_Fxyz_ab-2.0E0*2*I_NAI_Px_Fxyz_b;
  abcd[265] = 4.0E0*I_NAI_F2xy_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[266] = 4.0E0*I_NAI_F2xz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[267] = 4.0E0*I_NAI_Fx2y_Fxyz_ab;
  abcd[268] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[269] = 4.0E0*I_NAI_Fx2z_Fxyz_ab;
  abcd[270] = 4.0E0*I_NAI_F3x_F3y_ab-2.0E0*2*I_NAI_F3x_Py_a-2.0E0*2*I_NAI_Px_F3y_b+2*2*I_NAI_Px_Py;
  abcd[271] = 4.0E0*I_NAI_F2xy_F3y_ab-2.0E0*2*I_NAI_F2xy_Py_a-2.0E0*1*I_NAI_Py_F3y_b+2*I_NAI_Py_Py;
  abcd[272] = 4.0E0*I_NAI_F2xz_F3y_ab-2.0E0*2*I_NAI_F2xz_Py_a-2.0E0*1*I_NAI_Pz_F3y_b+2*I_NAI_Pz_Py;
  abcd[273] = 4.0E0*I_NAI_Fx2y_F3y_ab-2.0E0*2*I_NAI_Fx2y_Py_a;
  abcd[274] = 4.0E0*I_NAI_Fxyz_F3y_ab-2.0E0*2*I_NAI_Fxyz_Py_a;
  abcd[275] = 4.0E0*I_NAI_Fx2z_F3y_ab-2.0E0*2*I_NAI_Fx2z_Py_a;
  abcd[276] = 4.0E0*I_NAI_F3x_F2yz_ab-2.0E0*1*I_NAI_F3x_Pz_a-2.0E0*2*I_NAI_Px_F2yz_b+2*1*I_NAI_Px_Pz;
  abcd[277] = 4.0E0*I_NAI_F2xy_F2yz_ab-2.0E0*1*I_NAI_F2xy_Pz_a-2.0E0*1*I_NAI_Py_F2yz_b+1*I_NAI_Py_Pz;
  abcd[278] = 4.0E0*I_NAI_F2xz_F2yz_ab-2.0E0*1*I_NAI_F2xz_Pz_a-2.0E0*1*I_NAI_Pz_F2yz_b+1*I_NAI_Pz_Pz;
  abcd[279] = 4.0E0*I_NAI_Fx2y_F2yz_ab-2.0E0*1*I_NAI_Fx2y_Pz_a;
  abcd[280] = 4.0E0*I_NAI_Fxyz_F2yz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[281] = 4.0E0*I_NAI_Fx2z_F2yz_ab-2.0E0*1*I_NAI_Fx2z_Pz_a;
  abcd[282] = 4.0E0*I_NAI_F3x_Fy2z_ab-2.0E0*2*I_NAI_Px_Fy2z_b;
  abcd[283] = 4.0E0*I_NAI_F2xy_Fy2z_ab-2.0E0*1*I_NAI_Py_Fy2z_b;
  abcd[284] = 4.0E0*I_NAI_F2xz_Fy2z_ab-2.0E0*1*I_NAI_Pz_Fy2z_b;
  abcd[285] = 4.0E0*I_NAI_Fx2y_Fy2z_ab;
  abcd[286] = 4.0E0*I_NAI_Fxyz_Fy2z_ab;
  abcd[287] = 4.0E0*I_NAI_Fx2z_Fy2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[288] = 4.0E0*I_NAI_F3x_F2xz_ab-2.0E0*2*I_NAI_Px_F2xz_b;
  abcd[289] = 4.0E0*I_NAI_F2xy_F2xz_ab-2.0E0*1*I_NAI_Py_F2xz_b;
  abcd[290] = 4.0E0*I_NAI_F2xz_F2xz_ab-2.0E0*1*I_NAI_Pz_F2xz_b;
  abcd[291] = 4.0E0*I_NAI_Fx2y_F2xz_ab;
  abcd[292] = 4.0E0*I_NAI_Fxyz_F2xz_ab;
  abcd[293] = 4.0E0*I_NAI_Fx2z_F2xz_ab;
  abcd[294] = 4.0E0*I_NAI_F3x_Fxyz_ab-2.0E0*2*I_NAI_Px_Fxyz_b;
  abcd[295] = 4.0E0*I_NAI_F2xy_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[296] = 4.0E0*I_NAI_F2xz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[297] = 4.0E0*I_NAI_Fx2y_Fxyz_ab;
  abcd[298] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[299] = 4.0E0*I_NAI_Fx2z_Fxyz_ab;
  abcd[300] = 4.0E0*I_NAI_F3x_Fx2z_ab-2.0E0*1*I_NAI_F3x_Px_a-2.0E0*2*I_NAI_Px_Fx2z_b+2*1*I_NAI_Px_Px;
  abcd[301] = 4.0E0*I_NAI_F2xy_Fx2z_ab-2.0E0*1*I_NAI_F2xy_Px_a-2.0E0*1*I_NAI_Py_Fx2z_b+1*I_NAI_Py_Px;
  abcd[302] = 4.0E0*I_NAI_F2xz_Fx2z_ab-2.0E0*1*I_NAI_F2xz_Px_a-2.0E0*1*I_NAI_Pz_Fx2z_b+1*I_NAI_Pz_Px;
  abcd[303] = 4.0E0*I_NAI_Fx2y_Fx2z_ab-2.0E0*1*I_NAI_Fx2y_Px_a;
  abcd[304] = 4.0E0*I_NAI_Fxyz_Fx2z_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[305] = 4.0E0*I_NAI_Fx2z_Fx2z_ab-2.0E0*1*I_NAI_Fx2z_Px_a;
  abcd[306] = 4.0E0*I_NAI_F3x_F2yz_ab-2.0E0*2*I_NAI_Px_F2yz_b;
  abcd[307] = 4.0E0*I_NAI_F2xy_F2yz_ab-2.0E0*1*I_NAI_Py_F2yz_b;
  abcd[308] = 4.0E0*I_NAI_F2xz_F2yz_ab-2.0E0*1*I_NAI_Pz_F2yz_b;
  abcd[309] = 4.0E0*I_NAI_Fx2y_F2yz_ab;
  abcd[310] = 4.0E0*I_NAI_Fxyz_F2yz_ab;
  abcd[311] = 4.0E0*I_NAI_Fx2z_F2yz_ab;
  abcd[312] = 4.0E0*I_NAI_F3x_Fy2z_ab-2.0E0*1*I_NAI_F3x_Py_a-2.0E0*2*I_NAI_Px_Fy2z_b+2*1*I_NAI_Px_Py;
  abcd[313] = 4.0E0*I_NAI_F2xy_Fy2z_ab-2.0E0*1*I_NAI_F2xy_Py_a-2.0E0*1*I_NAI_Py_Fy2z_b+1*I_NAI_Py_Py;
  abcd[314] = 4.0E0*I_NAI_F2xz_Fy2z_ab-2.0E0*1*I_NAI_F2xz_Py_a-2.0E0*1*I_NAI_Pz_Fy2z_b+1*I_NAI_Pz_Py;
  abcd[315] = 4.0E0*I_NAI_Fx2y_Fy2z_ab-2.0E0*1*I_NAI_Fx2y_Py_a;
  abcd[316] = 4.0E0*I_NAI_Fxyz_Fy2z_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[317] = 4.0E0*I_NAI_Fx2z_Fy2z_ab-2.0E0*1*I_NAI_Fx2z_Py_a;
  abcd[318] = 4.0E0*I_NAI_F3x_F3z_ab-2.0E0*2*I_NAI_F3x_Pz_a-2.0E0*2*I_NAI_Px_F3z_b+2*2*I_NAI_Px_Pz;
  abcd[319] = 4.0E0*I_NAI_F2xy_F3z_ab-2.0E0*2*I_NAI_F2xy_Pz_a-2.0E0*1*I_NAI_Py_F3z_b+2*I_NAI_Py_Pz;
  abcd[320] = 4.0E0*I_NAI_F2xz_F3z_ab-2.0E0*2*I_NAI_F2xz_Pz_a-2.0E0*1*I_NAI_Pz_F3z_b+2*I_NAI_Pz_Pz;
  abcd[321] = 4.0E0*I_NAI_Fx2y_F3z_ab-2.0E0*2*I_NAI_Fx2y_Pz_a;
  abcd[322] = 4.0E0*I_NAI_Fxyz_F3z_ab-2.0E0*2*I_NAI_Fxyz_Pz_a;
  abcd[323] = 4.0E0*I_NAI_Fx2z_F3z_ab-2.0E0*2*I_NAI_Fx2z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[324] = 4.0E0*I_NAI_F2xy_F3x_ab-2.0E0*2*I_NAI_F2xy_Px_a;
  abcd[325] = 4.0E0*I_NAI_Fx2y_F3x_ab-2.0E0*2*I_NAI_Fx2y_Px_a-2.0E0*1*I_NAI_Px_F3x_b+2*I_NAI_Px_Px;
  abcd[326] = 4.0E0*I_NAI_Fxyz_F3x_ab-2.0E0*2*I_NAI_Fxyz_Px_a;
  abcd[327] = 4.0E0*I_NAI_F3y_F3x_ab-2.0E0*2*I_NAI_F3y_Px_a-2.0E0*2*I_NAI_Py_F3x_b+2*2*I_NAI_Py_Px;
  abcd[328] = 4.0E0*I_NAI_F2yz_F3x_ab-2.0E0*2*I_NAI_F2yz_Px_a-2.0E0*1*I_NAI_Pz_F3x_b+2*I_NAI_Pz_Px;
  abcd[329] = 4.0E0*I_NAI_Fy2z_F3x_ab-2.0E0*2*I_NAI_Fy2z_Px_a;
  abcd[330] = 4.0E0*I_NAI_F2xy_F2xy_ab-2.0E0*1*I_NAI_F2xy_Py_a;
  abcd[331] = 4.0E0*I_NAI_Fx2y_F2xy_ab-2.0E0*1*I_NAI_Fx2y_Py_a-2.0E0*1*I_NAI_Px_F2xy_b+1*I_NAI_Px_Py;
  abcd[332] = 4.0E0*I_NAI_Fxyz_F2xy_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[333] = 4.0E0*I_NAI_F3y_F2xy_ab-2.0E0*1*I_NAI_F3y_Py_a-2.0E0*2*I_NAI_Py_F2xy_b+2*1*I_NAI_Py_Py;
  abcd[334] = 4.0E0*I_NAI_F2yz_F2xy_ab-2.0E0*1*I_NAI_F2yz_Py_a-2.0E0*1*I_NAI_Pz_F2xy_b+1*I_NAI_Pz_Py;
  abcd[335] = 4.0E0*I_NAI_Fy2z_F2xy_ab-2.0E0*1*I_NAI_Fy2z_Py_a;
  abcd[336] = 4.0E0*I_NAI_F2xy_F2xz_ab-2.0E0*1*I_NAI_F2xy_Pz_a;
  abcd[337] = 4.0E0*I_NAI_Fx2y_F2xz_ab-2.0E0*1*I_NAI_Fx2y_Pz_a-2.0E0*1*I_NAI_Px_F2xz_b+1*I_NAI_Px_Pz;
  abcd[338] = 4.0E0*I_NAI_Fxyz_F2xz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[339] = 4.0E0*I_NAI_F3y_F2xz_ab-2.0E0*1*I_NAI_F3y_Pz_a-2.0E0*2*I_NAI_Py_F2xz_b+2*1*I_NAI_Py_Pz;
  abcd[340] = 4.0E0*I_NAI_F2yz_F2xz_ab-2.0E0*1*I_NAI_F2yz_Pz_a-2.0E0*1*I_NAI_Pz_F2xz_b+1*I_NAI_Pz_Pz;
  abcd[341] = 4.0E0*I_NAI_Fy2z_F2xz_ab-2.0E0*1*I_NAI_Fy2z_Pz_a;
  abcd[342] = 4.0E0*I_NAI_F2xy_Fx2y_ab;
  abcd[343] = 4.0E0*I_NAI_Fx2y_Fx2y_ab-2.0E0*1*I_NAI_Px_Fx2y_b;
  abcd[344] = 4.0E0*I_NAI_Fxyz_Fx2y_ab;
  abcd[345] = 4.0E0*I_NAI_F3y_Fx2y_ab-2.0E0*2*I_NAI_Py_Fx2y_b;
  abcd[346] = 4.0E0*I_NAI_F2yz_Fx2y_ab-2.0E0*1*I_NAI_Pz_Fx2y_b;
  abcd[347] = 4.0E0*I_NAI_Fy2z_Fx2y_ab;
  abcd[348] = 4.0E0*I_NAI_F2xy_Fxyz_ab;
  abcd[349] = 4.0E0*I_NAI_Fx2y_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[350] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[351] = 4.0E0*I_NAI_F3y_Fxyz_ab-2.0E0*2*I_NAI_Py_Fxyz_b;
  abcd[352] = 4.0E0*I_NAI_F2yz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[353] = 4.0E0*I_NAI_Fy2z_Fxyz_ab;
  abcd[354] = 4.0E0*I_NAI_F2xy_Fx2z_ab;
  abcd[355] = 4.0E0*I_NAI_Fx2y_Fx2z_ab-2.0E0*1*I_NAI_Px_Fx2z_b;
  abcd[356] = 4.0E0*I_NAI_Fxyz_Fx2z_ab;
  abcd[357] = 4.0E0*I_NAI_F3y_Fx2z_ab-2.0E0*2*I_NAI_Py_Fx2z_b;
  abcd[358] = 4.0E0*I_NAI_F2yz_Fx2z_ab-2.0E0*1*I_NAI_Pz_Fx2z_b;
  abcd[359] = 4.0E0*I_NAI_Fy2z_Fx2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[360] = 4.0E0*I_NAI_F2xy_F2xy_ab;
  abcd[361] = 4.0E0*I_NAI_Fx2y_F2xy_ab-2.0E0*1*I_NAI_Px_F2xy_b;
  abcd[362] = 4.0E0*I_NAI_Fxyz_F2xy_ab;
  abcd[363] = 4.0E0*I_NAI_F3y_F2xy_ab-2.0E0*2*I_NAI_Py_F2xy_b;
  abcd[364] = 4.0E0*I_NAI_F2yz_F2xy_ab-2.0E0*1*I_NAI_Pz_F2xy_b;
  abcd[365] = 4.0E0*I_NAI_Fy2z_F2xy_ab;
  abcd[366] = 4.0E0*I_NAI_F2xy_Fx2y_ab-2.0E0*1*I_NAI_F2xy_Px_a;
  abcd[367] = 4.0E0*I_NAI_Fx2y_Fx2y_ab-2.0E0*1*I_NAI_Fx2y_Px_a-2.0E0*1*I_NAI_Px_Fx2y_b+1*I_NAI_Px_Px;
  abcd[368] = 4.0E0*I_NAI_Fxyz_Fx2y_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[369] = 4.0E0*I_NAI_F3y_Fx2y_ab-2.0E0*1*I_NAI_F3y_Px_a-2.0E0*2*I_NAI_Py_Fx2y_b+2*1*I_NAI_Py_Px;
  abcd[370] = 4.0E0*I_NAI_F2yz_Fx2y_ab-2.0E0*1*I_NAI_F2yz_Px_a-2.0E0*1*I_NAI_Pz_Fx2y_b+1*I_NAI_Pz_Px;
  abcd[371] = 4.0E0*I_NAI_Fy2z_Fx2y_ab-2.0E0*1*I_NAI_Fy2z_Px_a;
  abcd[372] = 4.0E0*I_NAI_F2xy_Fxyz_ab;
  abcd[373] = 4.0E0*I_NAI_Fx2y_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[374] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[375] = 4.0E0*I_NAI_F3y_Fxyz_ab-2.0E0*2*I_NAI_Py_Fxyz_b;
  abcd[376] = 4.0E0*I_NAI_F2yz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[377] = 4.0E0*I_NAI_Fy2z_Fxyz_ab;
  abcd[378] = 4.0E0*I_NAI_F2xy_F3y_ab-2.0E0*2*I_NAI_F2xy_Py_a;
  abcd[379] = 4.0E0*I_NAI_Fx2y_F3y_ab-2.0E0*2*I_NAI_Fx2y_Py_a-2.0E0*1*I_NAI_Px_F3y_b+2*I_NAI_Px_Py;
  abcd[380] = 4.0E0*I_NAI_Fxyz_F3y_ab-2.0E0*2*I_NAI_Fxyz_Py_a;
  abcd[381] = 4.0E0*I_NAI_F3y_F3y_ab-2.0E0*2*I_NAI_F3y_Py_a-2.0E0*2*I_NAI_Py_F3y_b+2*2*I_NAI_Py_Py;
  abcd[382] = 4.0E0*I_NAI_F2yz_F3y_ab-2.0E0*2*I_NAI_F2yz_Py_a-2.0E0*1*I_NAI_Pz_F3y_b+2*I_NAI_Pz_Py;
  abcd[383] = 4.0E0*I_NAI_Fy2z_F3y_ab-2.0E0*2*I_NAI_Fy2z_Py_a;
  abcd[384] = 4.0E0*I_NAI_F2xy_F2yz_ab-2.0E0*1*I_NAI_F2xy_Pz_a;
  abcd[385] = 4.0E0*I_NAI_Fx2y_F2yz_ab-2.0E0*1*I_NAI_Fx2y_Pz_a-2.0E0*1*I_NAI_Px_F2yz_b+1*I_NAI_Px_Pz;
  abcd[386] = 4.0E0*I_NAI_Fxyz_F2yz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[387] = 4.0E0*I_NAI_F3y_F2yz_ab-2.0E0*1*I_NAI_F3y_Pz_a-2.0E0*2*I_NAI_Py_F2yz_b+2*1*I_NAI_Py_Pz;
  abcd[388] = 4.0E0*I_NAI_F2yz_F2yz_ab-2.0E0*1*I_NAI_F2yz_Pz_a-2.0E0*1*I_NAI_Pz_F2yz_b+1*I_NAI_Pz_Pz;
  abcd[389] = 4.0E0*I_NAI_Fy2z_F2yz_ab-2.0E0*1*I_NAI_Fy2z_Pz_a;
  abcd[390] = 4.0E0*I_NAI_F2xy_Fy2z_ab;
  abcd[391] = 4.0E0*I_NAI_Fx2y_Fy2z_ab-2.0E0*1*I_NAI_Px_Fy2z_b;
  abcd[392] = 4.0E0*I_NAI_Fxyz_Fy2z_ab;
  abcd[393] = 4.0E0*I_NAI_F3y_Fy2z_ab-2.0E0*2*I_NAI_Py_Fy2z_b;
  abcd[394] = 4.0E0*I_NAI_F2yz_Fy2z_ab-2.0E0*1*I_NAI_Pz_Fy2z_b;
  abcd[395] = 4.0E0*I_NAI_Fy2z_Fy2z_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[396] = 4.0E0*I_NAI_F2xy_F2xz_ab;
  abcd[397] = 4.0E0*I_NAI_Fx2y_F2xz_ab-2.0E0*1*I_NAI_Px_F2xz_b;
  abcd[398] = 4.0E0*I_NAI_Fxyz_F2xz_ab;
  abcd[399] = 4.0E0*I_NAI_F3y_F2xz_ab-2.0E0*2*I_NAI_Py_F2xz_b;
  abcd[400] = 4.0E0*I_NAI_F2yz_F2xz_ab-2.0E0*1*I_NAI_Pz_F2xz_b;
  abcd[401] = 4.0E0*I_NAI_Fy2z_F2xz_ab;
  abcd[402] = 4.0E0*I_NAI_F2xy_Fxyz_ab;
  abcd[403] = 4.0E0*I_NAI_Fx2y_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[404] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[405] = 4.0E0*I_NAI_F3y_Fxyz_ab-2.0E0*2*I_NAI_Py_Fxyz_b;
  abcd[406] = 4.0E0*I_NAI_F2yz_Fxyz_ab-2.0E0*1*I_NAI_Pz_Fxyz_b;
  abcd[407] = 4.0E0*I_NAI_Fy2z_Fxyz_ab;
  abcd[408] = 4.0E0*I_NAI_F2xy_Fx2z_ab-2.0E0*1*I_NAI_F2xy_Px_a;
  abcd[409] = 4.0E0*I_NAI_Fx2y_Fx2z_ab-2.0E0*1*I_NAI_Fx2y_Px_a-2.0E0*1*I_NAI_Px_Fx2z_b+1*I_NAI_Px_Px;
  abcd[410] = 4.0E0*I_NAI_Fxyz_Fx2z_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[411] = 4.0E0*I_NAI_F3y_Fx2z_ab-2.0E0*1*I_NAI_F3y_Px_a-2.0E0*2*I_NAI_Py_Fx2z_b+2*1*I_NAI_Py_Px;
  abcd[412] = 4.0E0*I_NAI_F2yz_Fx2z_ab-2.0E0*1*I_NAI_F2yz_Px_a-2.0E0*1*I_NAI_Pz_Fx2z_b+1*I_NAI_Pz_Px;
  abcd[413] = 4.0E0*I_NAI_Fy2z_Fx2z_ab-2.0E0*1*I_NAI_Fy2z_Px_a;
  abcd[414] = 4.0E0*I_NAI_F2xy_F2yz_ab;
  abcd[415] = 4.0E0*I_NAI_Fx2y_F2yz_ab-2.0E0*1*I_NAI_Px_F2yz_b;
  abcd[416] = 4.0E0*I_NAI_Fxyz_F2yz_ab;
  abcd[417] = 4.0E0*I_NAI_F3y_F2yz_ab-2.0E0*2*I_NAI_Py_F2yz_b;
  abcd[418] = 4.0E0*I_NAI_F2yz_F2yz_ab-2.0E0*1*I_NAI_Pz_F2yz_b;
  abcd[419] = 4.0E0*I_NAI_Fy2z_F2yz_ab;
  abcd[420] = 4.0E0*I_NAI_F2xy_Fy2z_ab-2.0E0*1*I_NAI_F2xy_Py_a;
  abcd[421] = 4.0E0*I_NAI_Fx2y_Fy2z_ab-2.0E0*1*I_NAI_Fx2y_Py_a-2.0E0*1*I_NAI_Px_Fy2z_b+1*I_NAI_Px_Py;
  abcd[422] = 4.0E0*I_NAI_Fxyz_Fy2z_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[423] = 4.0E0*I_NAI_F3y_Fy2z_ab-2.0E0*1*I_NAI_F3y_Py_a-2.0E0*2*I_NAI_Py_Fy2z_b+2*1*I_NAI_Py_Py;
  abcd[424] = 4.0E0*I_NAI_F2yz_Fy2z_ab-2.0E0*1*I_NAI_F2yz_Py_a-2.0E0*1*I_NAI_Pz_Fy2z_b+1*I_NAI_Pz_Py;
  abcd[425] = 4.0E0*I_NAI_Fy2z_Fy2z_ab-2.0E0*1*I_NAI_Fy2z_Py_a;
  abcd[426] = 4.0E0*I_NAI_F2xy_F3z_ab-2.0E0*2*I_NAI_F2xy_Pz_a;
  abcd[427] = 4.0E0*I_NAI_Fx2y_F3z_ab-2.0E0*2*I_NAI_Fx2y_Pz_a-2.0E0*1*I_NAI_Px_F3z_b+2*I_NAI_Px_Pz;
  abcd[428] = 4.0E0*I_NAI_Fxyz_F3z_ab-2.0E0*2*I_NAI_Fxyz_Pz_a;
  abcd[429] = 4.0E0*I_NAI_F3y_F3z_ab-2.0E0*2*I_NAI_F3y_Pz_a-2.0E0*2*I_NAI_Py_F3z_b+2*2*I_NAI_Py_Pz;
  abcd[430] = 4.0E0*I_NAI_F2yz_F3z_ab-2.0E0*2*I_NAI_F2yz_Pz_a-2.0E0*1*I_NAI_Pz_F3z_b+2*I_NAI_Pz_Pz;
  abcd[431] = 4.0E0*I_NAI_Fy2z_F3z_ab-2.0E0*2*I_NAI_Fy2z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[432] = 4.0E0*I_NAI_F2xz_F3x_ab-2.0E0*2*I_NAI_F2xz_Px_a;
  abcd[433] = 4.0E0*I_NAI_Fxyz_F3x_ab-2.0E0*2*I_NAI_Fxyz_Px_a;
  abcd[434] = 4.0E0*I_NAI_Fx2z_F3x_ab-2.0E0*2*I_NAI_Fx2z_Px_a-2.0E0*1*I_NAI_Px_F3x_b+2*I_NAI_Px_Px;
  abcd[435] = 4.0E0*I_NAI_F2yz_F3x_ab-2.0E0*2*I_NAI_F2yz_Px_a;
  abcd[436] = 4.0E0*I_NAI_Fy2z_F3x_ab-2.0E0*2*I_NAI_Fy2z_Px_a-2.0E0*1*I_NAI_Py_F3x_b+2*I_NAI_Py_Px;
  abcd[437] = 4.0E0*I_NAI_F3z_F3x_ab-2.0E0*2*I_NAI_F3z_Px_a-2.0E0*2*I_NAI_Pz_F3x_b+2*2*I_NAI_Pz_Px;
  abcd[438] = 4.0E0*I_NAI_F2xz_F2xy_ab-2.0E0*1*I_NAI_F2xz_Py_a;
  abcd[439] = 4.0E0*I_NAI_Fxyz_F2xy_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[440] = 4.0E0*I_NAI_Fx2z_F2xy_ab-2.0E0*1*I_NAI_Fx2z_Py_a-2.0E0*1*I_NAI_Px_F2xy_b+1*I_NAI_Px_Py;
  abcd[441] = 4.0E0*I_NAI_F2yz_F2xy_ab-2.0E0*1*I_NAI_F2yz_Py_a;
  abcd[442] = 4.0E0*I_NAI_Fy2z_F2xy_ab-2.0E0*1*I_NAI_Fy2z_Py_a-2.0E0*1*I_NAI_Py_F2xy_b+1*I_NAI_Py_Py;
  abcd[443] = 4.0E0*I_NAI_F3z_F2xy_ab-2.0E0*1*I_NAI_F3z_Py_a-2.0E0*2*I_NAI_Pz_F2xy_b+2*1*I_NAI_Pz_Py;
  abcd[444] = 4.0E0*I_NAI_F2xz_F2xz_ab-2.0E0*1*I_NAI_F2xz_Pz_a;
  abcd[445] = 4.0E0*I_NAI_Fxyz_F2xz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[446] = 4.0E0*I_NAI_Fx2z_F2xz_ab-2.0E0*1*I_NAI_Fx2z_Pz_a-2.0E0*1*I_NAI_Px_F2xz_b+1*I_NAI_Px_Pz;
  abcd[447] = 4.0E0*I_NAI_F2yz_F2xz_ab-2.0E0*1*I_NAI_F2yz_Pz_a;
  abcd[448] = 4.0E0*I_NAI_Fy2z_F2xz_ab-2.0E0*1*I_NAI_Fy2z_Pz_a-2.0E0*1*I_NAI_Py_F2xz_b+1*I_NAI_Py_Pz;
  abcd[449] = 4.0E0*I_NAI_F3z_F2xz_ab-2.0E0*1*I_NAI_F3z_Pz_a-2.0E0*2*I_NAI_Pz_F2xz_b+2*1*I_NAI_Pz_Pz;
  abcd[450] = 4.0E0*I_NAI_F2xz_Fx2y_ab;
  abcd[451] = 4.0E0*I_NAI_Fxyz_Fx2y_ab;
  abcd[452] = 4.0E0*I_NAI_Fx2z_Fx2y_ab-2.0E0*1*I_NAI_Px_Fx2y_b;
  abcd[453] = 4.0E0*I_NAI_F2yz_Fx2y_ab;
  abcd[454] = 4.0E0*I_NAI_Fy2z_Fx2y_ab-2.0E0*1*I_NAI_Py_Fx2y_b;
  abcd[455] = 4.0E0*I_NAI_F3z_Fx2y_ab-2.0E0*2*I_NAI_Pz_Fx2y_b;
  abcd[456] = 4.0E0*I_NAI_F2xz_Fxyz_ab;
  abcd[457] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[458] = 4.0E0*I_NAI_Fx2z_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[459] = 4.0E0*I_NAI_F2yz_Fxyz_ab;
  abcd[460] = 4.0E0*I_NAI_Fy2z_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[461] = 4.0E0*I_NAI_F3z_Fxyz_ab-2.0E0*2*I_NAI_Pz_Fxyz_b;
  abcd[462] = 4.0E0*I_NAI_F2xz_Fx2z_ab;
  abcd[463] = 4.0E0*I_NAI_Fxyz_Fx2z_ab;
  abcd[464] = 4.0E0*I_NAI_Fx2z_Fx2z_ab-2.0E0*1*I_NAI_Px_Fx2z_b;
  abcd[465] = 4.0E0*I_NAI_F2yz_Fx2z_ab;
  abcd[466] = 4.0E0*I_NAI_Fy2z_Fx2z_ab-2.0E0*1*I_NAI_Py_Fx2z_b;
  abcd[467] = 4.0E0*I_NAI_F3z_Fx2z_ab-2.0E0*2*I_NAI_Pz_Fx2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[468] = 4.0E0*I_NAI_F2xz_F2xy_ab;
  abcd[469] = 4.0E0*I_NAI_Fxyz_F2xy_ab;
  abcd[470] = 4.0E0*I_NAI_Fx2z_F2xy_ab-2.0E0*1*I_NAI_Px_F2xy_b;
  abcd[471] = 4.0E0*I_NAI_F2yz_F2xy_ab;
  abcd[472] = 4.0E0*I_NAI_Fy2z_F2xy_ab-2.0E0*1*I_NAI_Py_F2xy_b;
  abcd[473] = 4.0E0*I_NAI_F3z_F2xy_ab-2.0E0*2*I_NAI_Pz_F2xy_b;
  abcd[474] = 4.0E0*I_NAI_F2xz_Fx2y_ab-2.0E0*1*I_NAI_F2xz_Px_a;
  abcd[475] = 4.0E0*I_NAI_Fxyz_Fx2y_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[476] = 4.0E0*I_NAI_Fx2z_Fx2y_ab-2.0E0*1*I_NAI_Fx2z_Px_a-2.0E0*1*I_NAI_Px_Fx2y_b+1*I_NAI_Px_Px;
  abcd[477] = 4.0E0*I_NAI_F2yz_Fx2y_ab-2.0E0*1*I_NAI_F2yz_Px_a;
  abcd[478] = 4.0E0*I_NAI_Fy2z_Fx2y_ab-2.0E0*1*I_NAI_Fy2z_Px_a-2.0E0*1*I_NAI_Py_Fx2y_b+1*I_NAI_Py_Px;
  abcd[479] = 4.0E0*I_NAI_F3z_Fx2y_ab-2.0E0*1*I_NAI_F3z_Px_a-2.0E0*2*I_NAI_Pz_Fx2y_b+2*1*I_NAI_Pz_Px;
  abcd[480] = 4.0E0*I_NAI_F2xz_Fxyz_ab;
  abcd[481] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[482] = 4.0E0*I_NAI_Fx2z_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[483] = 4.0E0*I_NAI_F2yz_Fxyz_ab;
  abcd[484] = 4.0E0*I_NAI_Fy2z_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[485] = 4.0E0*I_NAI_F3z_Fxyz_ab-2.0E0*2*I_NAI_Pz_Fxyz_b;
  abcd[486] = 4.0E0*I_NAI_F2xz_F3y_ab-2.0E0*2*I_NAI_F2xz_Py_a;
  abcd[487] = 4.0E0*I_NAI_Fxyz_F3y_ab-2.0E0*2*I_NAI_Fxyz_Py_a;
  abcd[488] = 4.0E0*I_NAI_Fx2z_F3y_ab-2.0E0*2*I_NAI_Fx2z_Py_a-2.0E0*1*I_NAI_Px_F3y_b+2*I_NAI_Px_Py;
  abcd[489] = 4.0E0*I_NAI_F2yz_F3y_ab-2.0E0*2*I_NAI_F2yz_Py_a;
  abcd[490] = 4.0E0*I_NAI_Fy2z_F3y_ab-2.0E0*2*I_NAI_Fy2z_Py_a-2.0E0*1*I_NAI_Py_F3y_b+2*I_NAI_Py_Py;
  abcd[491] = 4.0E0*I_NAI_F3z_F3y_ab-2.0E0*2*I_NAI_F3z_Py_a-2.0E0*2*I_NAI_Pz_F3y_b+2*2*I_NAI_Pz_Py;
  abcd[492] = 4.0E0*I_NAI_F2xz_F2yz_ab-2.0E0*1*I_NAI_F2xz_Pz_a;
  abcd[493] = 4.0E0*I_NAI_Fxyz_F2yz_ab-2.0E0*1*I_NAI_Fxyz_Pz_a;
  abcd[494] = 4.0E0*I_NAI_Fx2z_F2yz_ab-2.0E0*1*I_NAI_Fx2z_Pz_a-2.0E0*1*I_NAI_Px_F2yz_b+1*I_NAI_Px_Pz;
  abcd[495] = 4.0E0*I_NAI_F2yz_F2yz_ab-2.0E0*1*I_NAI_F2yz_Pz_a;
  abcd[496] = 4.0E0*I_NAI_Fy2z_F2yz_ab-2.0E0*1*I_NAI_Fy2z_Pz_a-2.0E0*1*I_NAI_Py_F2yz_b+1*I_NAI_Py_Pz;
  abcd[497] = 4.0E0*I_NAI_F3z_F2yz_ab-2.0E0*1*I_NAI_F3z_Pz_a-2.0E0*2*I_NAI_Pz_F2yz_b+2*1*I_NAI_Pz_Pz;
  abcd[498] = 4.0E0*I_NAI_F2xz_Fy2z_ab;
  abcd[499] = 4.0E0*I_NAI_Fxyz_Fy2z_ab;
  abcd[500] = 4.0E0*I_NAI_Fx2z_Fy2z_ab-2.0E0*1*I_NAI_Px_Fy2z_b;
  abcd[501] = 4.0E0*I_NAI_F2yz_Fy2z_ab;
  abcd[502] = 4.0E0*I_NAI_Fy2z_Fy2z_ab-2.0E0*1*I_NAI_Py_Fy2z_b;
  abcd[503] = 4.0E0*I_NAI_F3z_Fy2z_ab-2.0E0*2*I_NAI_Pz_Fy2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_F_ab
   * RHS shell quartet name: SQ_NAI_F_P_a
   * RHS shell quartet name: SQ_NAI_P_F_b
   * RHS shell quartet name: SQ_NAI_P_P
   ************************************************************/
  abcd[504] = 4.0E0*I_NAI_F2xz_F2xz_ab;
  abcd[505] = 4.0E0*I_NAI_Fxyz_F2xz_ab;
  abcd[506] = 4.0E0*I_NAI_Fx2z_F2xz_ab-2.0E0*1*I_NAI_Px_F2xz_b;
  abcd[507] = 4.0E0*I_NAI_F2yz_F2xz_ab;
  abcd[508] = 4.0E0*I_NAI_Fy2z_F2xz_ab-2.0E0*1*I_NAI_Py_F2xz_b;
  abcd[509] = 4.0E0*I_NAI_F3z_F2xz_ab-2.0E0*2*I_NAI_Pz_F2xz_b;
  abcd[510] = 4.0E0*I_NAI_F2xz_Fxyz_ab;
  abcd[511] = 4.0E0*I_NAI_Fxyz_Fxyz_ab;
  abcd[512] = 4.0E0*I_NAI_Fx2z_Fxyz_ab-2.0E0*1*I_NAI_Px_Fxyz_b;
  abcd[513] = 4.0E0*I_NAI_F2yz_Fxyz_ab;
  abcd[514] = 4.0E0*I_NAI_Fy2z_Fxyz_ab-2.0E0*1*I_NAI_Py_Fxyz_b;
  abcd[515] = 4.0E0*I_NAI_F3z_Fxyz_ab-2.0E0*2*I_NAI_Pz_Fxyz_b;
  abcd[516] = 4.0E0*I_NAI_F2xz_Fx2z_ab-2.0E0*1*I_NAI_F2xz_Px_a;
  abcd[517] = 4.0E0*I_NAI_Fxyz_Fx2z_ab-2.0E0*1*I_NAI_Fxyz_Px_a;
  abcd[518] = 4.0E0*I_NAI_Fx2z_Fx2z_ab-2.0E0*1*I_NAI_Fx2z_Px_a-2.0E0*1*I_NAI_Px_Fx2z_b+1*I_NAI_Px_Px;
  abcd[519] = 4.0E0*I_NAI_F2yz_Fx2z_ab-2.0E0*1*I_NAI_F2yz_Px_a;
  abcd[520] = 4.0E0*I_NAI_Fy2z_Fx2z_ab-2.0E0*1*I_NAI_Fy2z_Px_a-2.0E0*1*I_NAI_Py_Fx2z_b+1*I_NAI_Py_Px;
  abcd[521] = 4.0E0*I_NAI_F3z_Fx2z_ab-2.0E0*1*I_NAI_F3z_Px_a-2.0E0*2*I_NAI_Pz_Fx2z_b+2*1*I_NAI_Pz_Px;
  abcd[522] = 4.0E0*I_NAI_F2xz_F2yz_ab;
  abcd[523] = 4.0E0*I_NAI_Fxyz_F2yz_ab;
  abcd[524] = 4.0E0*I_NAI_Fx2z_F2yz_ab-2.0E0*1*I_NAI_Px_F2yz_b;
  abcd[525] = 4.0E0*I_NAI_F2yz_F2yz_ab;
  abcd[526] = 4.0E0*I_NAI_Fy2z_F2yz_ab-2.0E0*1*I_NAI_Py_F2yz_b;
  abcd[527] = 4.0E0*I_NAI_F3z_F2yz_ab-2.0E0*2*I_NAI_Pz_F2yz_b;
  abcd[528] = 4.0E0*I_NAI_F2xz_Fy2z_ab-2.0E0*1*I_NAI_F2xz_Py_a;
  abcd[529] = 4.0E0*I_NAI_Fxyz_Fy2z_ab-2.0E0*1*I_NAI_Fxyz_Py_a;
  abcd[530] = 4.0E0*I_NAI_Fx2z_Fy2z_ab-2.0E0*1*I_NAI_Fx2z_Py_a-2.0E0*1*I_NAI_Px_Fy2z_b+1*I_NAI_Px_Py;
  abcd[531] = 4.0E0*I_NAI_F2yz_Fy2z_ab-2.0E0*1*I_NAI_F2yz_Py_a;
  abcd[532] = 4.0E0*I_NAI_Fy2z_Fy2z_ab-2.0E0*1*I_NAI_Fy2z_Py_a-2.0E0*1*I_NAI_Py_Fy2z_b+1*I_NAI_Py_Py;
  abcd[533] = 4.0E0*I_NAI_F3z_Fy2z_ab-2.0E0*1*I_NAI_F3z_Py_a-2.0E0*2*I_NAI_Pz_Fy2z_b+2*1*I_NAI_Pz_Py;
  abcd[534] = 4.0E0*I_NAI_F2xz_F3z_ab-2.0E0*2*I_NAI_F2xz_Pz_a;
  abcd[535] = 4.0E0*I_NAI_Fxyz_F3z_ab-2.0E0*2*I_NAI_Fxyz_Pz_a;
  abcd[536] = 4.0E0*I_NAI_Fx2z_F3z_ab-2.0E0*2*I_NAI_Fx2z_Pz_a-2.0E0*1*I_NAI_Px_F3z_b+2*I_NAI_Px_Pz;
  abcd[537] = 4.0E0*I_NAI_F2yz_F3z_ab-2.0E0*2*I_NAI_F2yz_Pz_a;
  abcd[538] = 4.0E0*I_NAI_Fy2z_F3z_ab-2.0E0*2*I_NAI_Fy2z_Pz_a-2.0E0*1*I_NAI_Py_F3z_b+2*I_NAI_Py_Pz;
  abcd[539] = 4.0E0*I_NAI_F3z_F3z_ab-2.0E0*2*I_NAI_F3z_Pz_a-2.0E0*2*I_NAI_Pz_F3z_b+2*2*I_NAI_Pz_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[540] = 4.0E0*I_NAI_D2x_G4x_bb-2.0E0*2*I_NAI_D2x_D2x_b-2.0E0*3*I_NAI_D2x_D2x_b+2*1*I_NAI_D2x_S;
  abcd[541] = 4.0E0*I_NAI_Dxy_G4x_bb-2.0E0*2*I_NAI_Dxy_D2x_b-2.0E0*3*I_NAI_Dxy_D2x_b+2*1*I_NAI_Dxy_S;
  abcd[542] = 4.0E0*I_NAI_Dxz_G4x_bb-2.0E0*2*I_NAI_Dxz_D2x_b-2.0E0*3*I_NAI_Dxz_D2x_b+2*1*I_NAI_Dxz_S;
  abcd[543] = 4.0E0*I_NAI_D2y_G4x_bb-2.0E0*2*I_NAI_D2y_D2x_b-2.0E0*3*I_NAI_D2y_D2x_b+2*1*I_NAI_D2y_S;
  abcd[544] = 4.0E0*I_NAI_Dyz_G4x_bb-2.0E0*2*I_NAI_Dyz_D2x_b-2.0E0*3*I_NAI_Dyz_D2x_b+2*1*I_NAI_Dyz_S;
  abcd[545] = 4.0E0*I_NAI_D2z_G4x_bb-2.0E0*2*I_NAI_D2z_D2x_b-2.0E0*3*I_NAI_D2z_D2x_b+2*1*I_NAI_D2z_S;
  abcd[546] = 4.0E0*I_NAI_D2x_G3xy_bb-2.0E0*1*I_NAI_D2x_Dxy_b-2.0E0*2*I_NAI_D2x_Dxy_b;
  abcd[547] = 4.0E0*I_NAI_Dxy_G3xy_bb-2.0E0*1*I_NAI_Dxy_Dxy_b-2.0E0*2*I_NAI_Dxy_Dxy_b;
  abcd[548] = 4.0E0*I_NAI_Dxz_G3xy_bb-2.0E0*1*I_NAI_Dxz_Dxy_b-2.0E0*2*I_NAI_Dxz_Dxy_b;
  abcd[549] = 4.0E0*I_NAI_D2y_G3xy_bb-2.0E0*1*I_NAI_D2y_Dxy_b-2.0E0*2*I_NAI_D2y_Dxy_b;
  abcd[550] = 4.0E0*I_NAI_Dyz_G3xy_bb-2.0E0*1*I_NAI_Dyz_Dxy_b-2.0E0*2*I_NAI_Dyz_Dxy_b;
  abcd[551] = 4.0E0*I_NAI_D2z_G3xy_bb-2.0E0*1*I_NAI_D2z_Dxy_b-2.0E0*2*I_NAI_D2z_Dxy_b;
  abcd[552] = 4.0E0*I_NAI_D2x_G3xz_bb-2.0E0*1*I_NAI_D2x_Dxz_b-2.0E0*2*I_NAI_D2x_Dxz_b;
  abcd[553] = 4.0E0*I_NAI_Dxy_G3xz_bb-2.0E0*1*I_NAI_Dxy_Dxz_b-2.0E0*2*I_NAI_Dxy_Dxz_b;
  abcd[554] = 4.0E0*I_NAI_Dxz_G3xz_bb-2.0E0*1*I_NAI_Dxz_Dxz_b-2.0E0*2*I_NAI_Dxz_Dxz_b;
  abcd[555] = 4.0E0*I_NAI_D2y_G3xz_bb-2.0E0*1*I_NAI_D2y_Dxz_b-2.0E0*2*I_NAI_D2y_Dxz_b;
  abcd[556] = 4.0E0*I_NAI_Dyz_G3xz_bb-2.0E0*1*I_NAI_Dyz_Dxz_b-2.0E0*2*I_NAI_Dyz_Dxz_b;
  abcd[557] = 4.0E0*I_NAI_D2z_G3xz_bb-2.0E0*1*I_NAI_D2z_Dxz_b-2.0E0*2*I_NAI_D2z_Dxz_b;
  abcd[558] = 4.0E0*I_NAI_D2x_G2x2y_bb-2.0E0*1*I_NAI_D2x_D2y_b;
  abcd[559] = 4.0E0*I_NAI_Dxy_G2x2y_bb-2.0E0*1*I_NAI_Dxy_D2y_b;
  abcd[560] = 4.0E0*I_NAI_Dxz_G2x2y_bb-2.0E0*1*I_NAI_Dxz_D2y_b;
  abcd[561] = 4.0E0*I_NAI_D2y_G2x2y_bb-2.0E0*1*I_NAI_D2y_D2y_b;
  abcd[562] = 4.0E0*I_NAI_Dyz_G2x2y_bb-2.0E0*1*I_NAI_Dyz_D2y_b;
  abcd[563] = 4.0E0*I_NAI_D2z_G2x2y_bb-2.0E0*1*I_NAI_D2z_D2y_b;
  abcd[564] = 4.0E0*I_NAI_D2x_G2xyz_bb-2.0E0*1*I_NAI_D2x_Dyz_b;
  abcd[565] = 4.0E0*I_NAI_Dxy_G2xyz_bb-2.0E0*1*I_NAI_Dxy_Dyz_b;
  abcd[566] = 4.0E0*I_NAI_Dxz_G2xyz_bb-2.0E0*1*I_NAI_Dxz_Dyz_b;
  abcd[567] = 4.0E0*I_NAI_D2y_G2xyz_bb-2.0E0*1*I_NAI_D2y_Dyz_b;
  abcd[568] = 4.0E0*I_NAI_Dyz_G2xyz_bb-2.0E0*1*I_NAI_Dyz_Dyz_b;
  abcd[569] = 4.0E0*I_NAI_D2z_G2xyz_bb-2.0E0*1*I_NAI_D2z_Dyz_b;
  abcd[570] = 4.0E0*I_NAI_D2x_G2x2z_bb-2.0E0*1*I_NAI_D2x_D2z_b;
  abcd[571] = 4.0E0*I_NAI_Dxy_G2x2z_bb-2.0E0*1*I_NAI_Dxy_D2z_b;
  abcd[572] = 4.0E0*I_NAI_Dxz_G2x2z_bb-2.0E0*1*I_NAI_Dxz_D2z_b;
  abcd[573] = 4.0E0*I_NAI_D2y_G2x2z_bb-2.0E0*1*I_NAI_D2y_D2z_b;
  abcd[574] = 4.0E0*I_NAI_Dyz_G2x2z_bb-2.0E0*1*I_NAI_Dyz_D2z_b;
  abcd[575] = 4.0E0*I_NAI_D2z_G2x2z_bb-2.0E0*1*I_NAI_D2z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[576] = 4.0E0*I_NAI_D2x_G3xy_bb-2.0E0*2*I_NAI_D2x_Dxy_b;
  abcd[577] = 4.0E0*I_NAI_Dxy_G3xy_bb-2.0E0*2*I_NAI_Dxy_Dxy_b;
  abcd[578] = 4.0E0*I_NAI_Dxz_G3xy_bb-2.0E0*2*I_NAI_Dxz_Dxy_b;
  abcd[579] = 4.0E0*I_NAI_D2y_G3xy_bb-2.0E0*2*I_NAI_D2y_Dxy_b;
  abcd[580] = 4.0E0*I_NAI_Dyz_G3xy_bb-2.0E0*2*I_NAI_Dyz_Dxy_b;
  abcd[581] = 4.0E0*I_NAI_D2z_G3xy_bb-2.0E0*2*I_NAI_D2z_Dxy_b;
  abcd[582] = 4.0E0*I_NAI_D2x_G2x2y_bb-2.0E0*1*I_NAI_D2x_D2x_b-2.0E0*1*I_NAI_D2x_D2y_b+1*I_NAI_D2x_S;
  abcd[583] = 4.0E0*I_NAI_Dxy_G2x2y_bb-2.0E0*1*I_NAI_Dxy_D2x_b-2.0E0*1*I_NAI_Dxy_D2y_b+1*I_NAI_Dxy_S;
  abcd[584] = 4.0E0*I_NAI_Dxz_G2x2y_bb-2.0E0*1*I_NAI_Dxz_D2x_b-2.0E0*1*I_NAI_Dxz_D2y_b+1*I_NAI_Dxz_S;
  abcd[585] = 4.0E0*I_NAI_D2y_G2x2y_bb-2.0E0*1*I_NAI_D2y_D2x_b-2.0E0*1*I_NAI_D2y_D2y_b+1*I_NAI_D2y_S;
  abcd[586] = 4.0E0*I_NAI_Dyz_G2x2y_bb-2.0E0*1*I_NAI_Dyz_D2x_b-2.0E0*1*I_NAI_Dyz_D2y_b+1*I_NAI_Dyz_S;
  abcd[587] = 4.0E0*I_NAI_D2z_G2x2y_bb-2.0E0*1*I_NAI_D2z_D2x_b-2.0E0*1*I_NAI_D2z_D2y_b+1*I_NAI_D2z_S;
  abcd[588] = 4.0E0*I_NAI_D2x_G2xyz_bb-2.0E0*1*I_NAI_D2x_Dyz_b;
  abcd[589] = 4.0E0*I_NAI_Dxy_G2xyz_bb-2.0E0*1*I_NAI_Dxy_Dyz_b;
  abcd[590] = 4.0E0*I_NAI_Dxz_G2xyz_bb-2.0E0*1*I_NAI_Dxz_Dyz_b;
  abcd[591] = 4.0E0*I_NAI_D2y_G2xyz_bb-2.0E0*1*I_NAI_D2y_Dyz_b;
  abcd[592] = 4.0E0*I_NAI_Dyz_G2xyz_bb-2.0E0*1*I_NAI_Dyz_Dyz_b;
  abcd[593] = 4.0E0*I_NAI_D2z_G2xyz_bb-2.0E0*1*I_NAI_D2z_Dyz_b;
  abcd[594] = 4.0E0*I_NAI_D2x_Gx3y_bb-2.0E0*2*I_NAI_D2x_Dxy_b;
  abcd[595] = 4.0E0*I_NAI_Dxy_Gx3y_bb-2.0E0*2*I_NAI_Dxy_Dxy_b;
  abcd[596] = 4.0E0*I_NAI_Dxz_Gx3y_bb-2.0E0*2*I_NAI_Dxz_Dxy_b;
  abcd[597] = 4.0E0*I_NAI_D2y_Gx3y_bb-2.0E0*2*I_NAI_D2y_Dxy_b;
  abcd[598] = 4.0E0*I_NAI_Dyz_Gx3y_bb-2.0E0*2*I_NAI_Dyz_Dxy_b;
  abcd[599] = 4.0E0*I_NAI_D2z_Gx3y_bb-2.0E0*2*I_NAI_D2z_Dxy_b;
  abcd[600] = 4.0E0*I_NAI_D2x_Gx2yz_bb-2.0E0*1*I_NAI_D2x_Dxz_b;
  abcd[601] = 4.0E0*I_NAI_Dxy_Gx2yz_bb-2.0E0*1*I_NAI_Dxy_Dxz_b;
  abcd[602] = 4.0E0*I_NAI_Dxz_Gx2yz_bb-2.0E0*1*I_NAI_Dxz_Dxz_b;
  abcd[603] = 4.0E0*I_NAI_D2y_Gx2yz_bb-2.0E0*1*I_NAI_D2y_Dxz_b;
  abcd[604] = 4.0E0*I_NAI_Dyz_Gx2yz_bb-2.0E0*1*I_NAI_Dyz_Dxz_b;
  abcd[605] = 4.0E0*I_NAI_D2z_Gx2yz_bb-2.0E0*1*I_NAI_D2z_Dxz_b;
  abcd[606] = 4.0E0*I_NAI_D2x_Gxy2z_bb;
  abcd[607] = 4.0E0*I_NAI_Dxy_Gxy2z_bb;
  abcd[608] = 4.0E0*I_NAI_Dxz_Gxy2z_bb;
  abcd[609] = 4.0E0*I_NAI_D2y_Gxy2z_bb;
  abcd[610] = 4.0E0*I_NAI_Dyz_Gxy2z_bb;
  abcd[611] = 4.0E0*I_NAI_D2z_Gxy2z_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[612] = 4.0E0*I_NAI_D2x_G3xz_bb-2.0E0*2*I_NAI_D2x_Dxz_b;
  abcd[613] = 4.0E0*I_NAI_Dxy_G3xz_bb-2.0E0*2*I_NAI_Dxy_Dxz_b;
  abcd[614] = 4.0E0*I_NAI_Dxz_G3xz_bb-2.0E0*2*I_NAI_Dxz_Dxz_b;
  abcd[615] = 4.0E0*I_NAI_D2y_G3xz_bb-2.0E0*2*I_NAI_D2y_Dxz_b;
  abcd[616] = 4.0E0*I_NAI_Dyz_G3xz_bb-2.0E0*2*I_NAI_Dyz_Dxz_b;
  abcd[617] = 4.0E0*I_NAI_D2z_G3xz_bb-2.0E0*2*I_NAI_D2z_Dxz_b;
  abcd[618] = 4.0E0*I_NAI_D2x_G2xyz_bb-2.0E0*1*I_NAI_D2x_Dyz_b;
  abcd[619] = 4.0E0*I_NAI_Dxy_G2xyz_bb-2.0E0*1*I_NAI_Dxy_Dyz_b;
  abcd[620] = 4.0E0*I_NAI_Dxz_G2xyz_bb-2.0E0*1*I_NAI_Dxz_Dyz_b;
  abcd[621] = 4.0E0*I_NAI_D2y_G2xyz_bb-2.0E0*1*I_NAI_D2y_Dyz_b;
  abcd[622] = 4.0E0*I_NAI_Dyz_G2xyz_bb-2.0E0*1*I_NAI_Dyz_Dyz_b;
  abcd[623] = 4.0E0*I_NAI_D2z_G2xyz_bb-2.0E0*1*I_NAI_D2z_Dyz_b;
  abcd[624] = 4.0E0*I_NAI_D2x_G2x2z_bb-2.0E0*1*I_NAI_D2x_D2x_b-2.0E0*1*I_NAI_D2x_D2z_b+1*I_NAI_D2x_S;
  abcd[625] = 4.0E0*I_NAI_Dxy_G2x2z_bb-2.0E0*1*I_NAI_Dxy_D2x_b-2.0E0*1*I_NAI_Dxy_D2z_b+1*I_NAI_Dxy_S;
  abcd[626] = 4.0E0*I_NAI_Dxz_G2x2z_bb-2.0E0*1*I_NAI_Dxz_D2x_b-2.0E0*1*I_NAI_Dxz_D2z_b+1*I_NAI_Dxz_S;
  abcd[627] = 4.0E0*I_NAI_D2y_G2x2z_bb-2.0E0*1*I_NAI_D2y_D2x_b-2.0E0*1*I_NAI_D2y_D2z_b+1*I_NAI_D2y_S;
  abcd[628] = 4.0E0*I_NAI_Dyz_G2x2z_bb-2.0E0*1*I_NAI_Dyz_D2x_b-2.0E0*1*I_NAI_Dyz_D2z_b+1*I_NAI_Dyz_S;
  abcd[629] = 4.0E0*I_NAI_D2z_G2x2z_bb-2.0E0*1*I_NAI_D2z_D2x_b-2.0E0*1*I_NAI_D2z_D2z_b+1*I_NAI_D2z_S;
  abcd[630] = 4.0E0*I_NAI_D2x_Gx2yz_bb;
  abcd[631] = 4.0E0*I_NAI_Dxy_Gx2yz_bb;
  abcd[632] = 4.0E0*I_NAI_Dxz_Gx2yz_bb;
  abcd[633] = 4.0E0*I_NAI_D2y_Gx2yz_bb;
  abcd[634] = 4.0E0*I_NAI_Dyz_Gx2yz_bb;
  abcd[635] = 4.0E0*I_NAI_D2z_Gx2yz_bb;
  abcd[636] = 4.0E0*I_NAI_D2x_Gxy2z_bb-2.0E0*1*I_NAI_D2x_Dxy_b;
  abcd[637] = 4.0E0*I_NAI_Dxy_Gxy2z_bb-2.0E0*1*I_NAI_Dxy_Dxy_b;
  abcd[638] = 4.0E0*I_NAI_Dxz_Gxy2z_bb-2.0E0*1*I_NAI_Dxz_Dxy_b;
  abcd[639] = 4.0E0*I_NAI_D2y_Gxy2z_bb-2.0E0*1*I_NAI_D2y_Dxy_b;
  abcd[640] = 4.0E0*I_NAI_Dyz_Gxy2z_bb-2.0E0*1*I_NAI_Dyz_Dxy_b;
  abcd[641] = 4.0E0*I_NAI_D2z_Gxy2z_bb-2.0E0*1*I_NAI_D2z_Dxy_b;
  abcd[642] = 4.0E0*I_NAI_D2x_Gx3z_bb-2.0E0*2*I_NAI_D2x_Dxz_b;
  abcd[643] = 4.0E0*I_NAI_Dxy_Gx3z_bb-2.0E0*2*I_NAI_Dxy_Dxz_b;
  abcd[644] = 4.0E0*I_NAI_Dxz_Gx3z_bb-2.0E0*2*I_NAI_Dxz_Dxz_b;
  abcd[645] = 4.0E0*I_NAI_D2y_Gx3z_bb-2.0E0*2*I_NAI_D2y_Dxz_b;
  abcd[646] = 4.0E0*I_NAI_Dyz_Gx3z_bb-2.0E0*2*I_NAI_Dyz_Dxz_b;
  abcd[647] = 4.0E0*I_NAI_D2z_Gx3z_bb-2.0E0*2*I_NAI_D2z_Dxz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[648] = 4.0E0*I_NAI_D2x_G2x2y_bb-2.0E0*1*I_NAI_D2x_D2x_b;
  abcd[649] = 4.0E0*I_NAI_Dxy_G2x2y_bb-2.0E0*1*I_NAI_Dxy_D2x_b;
  abcd[650] = 4.0E0*I_NAI_Dxz_G2x2y_bb-2.0E0*1*I_NAI_Dxz_D2x_b;
  abcd[651] = 4.0E0*I_NAI_D2y_G2x2y_bb-2.0E0*1*I_NAI_D2y_D2x_b;
  abcd[652] = 4.0E0*I_NAI_Dyz_G2x2y_bb-2.0E0*1*I_NAI_Dyz_D2x_b;
  abcd[653] = 4.0E0*I_NAI_D2z_G2x2y_bb-2.0E0*1*I_NAI_D2z_D2x_b;
  abcd[654] = 4.0E0*I_NAI_D2x_Gx3y_bb-2.0E0*1*I_NAI_D2x_Dxy_b-2.0E0*2*I_NAI_D2x_Dxy_b;
  abcd[655] = 4.0E0*I_NAI_Dxy_Gx3y_bb-2.0E0*1*I_NAI_Dxy_Dxy_b-2.0E0*2*I_NAI_Dxy_Dxy_b;
  abcd[656] = 4.0E0*I_NAI_Dxz_Gx3y_bb-2.0E0*1*I_NAI_Dxz_Dxy_b-2.0E0*2*I_NAI_Dxz_Dxy_b;
  abcd[657] = 4.0E0*I_NAI_D2y_Gx3y_bb-2.0E0*1*I_NAI_D2y_Dxy_b-2.0E0*2*I_NAI_D2y_Dxy_b;
  abcd[658] = 4.0E0*I_NAI_Dyz_Gx3y_bb-2.0E0*1*I_NAI_Dyz_Dxy_b-2.0E0*2*I_NAI_Dyz_Dxy_b;
  abcd[659] = 4.0E0*I_NAI_D2z_Gx3y_bb-2.0E0*1*I_NAI_D2z_Dxy_b-2.0E0*2*I_NAI_D2z_Dxy_b;
  abcd[660] = 4.0E0*I_NAI_D2x_Gx2yz_bb-2.0E0*1*I_NAI_D2x_Dxz_b;
  abcd[661] = 4.0E0*I_NAI_Dxy_Gx2yz_bb-2.0E0*1*I_NAI_Dxy_Dxz_b;
  abcd[662] = 4.0E0*I_NAI_Dxz_Gx2yz_bb-2.0E0*1*I_NAI_Dxz_Dxz_b;
  abcd[663] = 4.0E0*I_NAI_D2y_Gx2yz_bb-2.0E0*1*I_NAI_D2y_Dxz_b;
  abcd[664] = 4.0E0*I_NAI_Dyz_Gx2yz_bb-2.0E0*1*I_NAI_Dyz_Dxz_b;
  abcd[665] = 4.0E0*I_NAI_D2z_Gx2yz_bb-2.0E0*1*I_NAI_D2z_Dxz_b;
  abcd[666] = 4.0E0*I_NAI_D2x_G4y_bb-2.0E0*2*I_NAI_D2x_D2y_b-2.0E0*3*I_NAI_D2x_D2y_b+2*1*I_NAI_D2x_S;
  abcd[667] = 4.0E0*I_NAI_Dxy_G4y_bb-2.0E0*2*I_NAI_Dxy_D2y_b-2.0E0*3*I_NAI_Dxy_D2y_b+2*1*I_NAI_Dxy_S;
  abcd[668] = 4.0E0*I_NAI_Dxz_G4y_bb-2.0E0*2*I_NAI_Dxz_D2y_b-2.0E0*3*I_NAI_Dxz_D2y_b+2*1*I_NAI_Dxz_S;
  abcd[669] = 4.0E0*I_NAI_D2y_G4y_bb-2.0E0*2*I_NAI_D2y_D2y_b-2.0E0*3*I_NAI_D2y_D2y_b+2*1*I_NAI_D2y_S;
  abcd[670] = 4.0E0*I_NAI_Dyz_G4y_bb-2.0E0*2*I_NAI_Dyz_D2y_b-2.0E0*3*I_NAI_Dyz_D2y_b+2*1*I_NAI_Dyz_S;
  abcd[671] = 4.0E0*I_NAI_D2z_G4y_bb-2.0E0*2*I_NAI_D2z_D2y_b-2.0E0*3*I_NAI_D2z_D2y_b+2*1*I_NAI_D2z_S;
  abcd[672] = 4.0E0*I_NAI_D2x_G3yz_bb-2.0E0*1*I_NAI_D2x_Dyz_b-2.0E0*2*I_NAI_D2x_Dyz_b;
  abcd[673] = 4.0E0*I_NAI_Dxy_G3yz_bb-2.0E0*1*I_NAI_Dxy_Dyz_b-2.0E0*2*I_NAI_Dxy_Dyz_b;
  abcd[674] = 4.0E0*I_NAI_Dxz_G3yz_bb-2.0E0*1*I_NAI_Dxz_Dyz_b-2.0E0*2*I_NAI_Dxz_Dyz_b;
  abcd[675] = 4.0E0*I_NAI_D2y_G3yz_bb-2.0E0*1*I_NAI_D2y_Dyz_b-2.0E0*2*I_NAI_D2y_Dyz_b;
  abcd[676] = 4.0E0*I_NAI_Dyz_G3yz_bb-2.0E0*1*I_NAI_Dyz_Dyz_b-2.0E0*2*I_NAI_Dyz_Dyz_b;
  abcd[677] = 4.0E0*I_NAI_D2z_G3yz_bb-2.0E0*1*I_NAI_D2z_Dyz_b-2.0E0*2*I_NAI_D2z_Dyz_b;
  abcd[678] = 4.0E0*I_NAI_D2x_G2y2z_bb-2.0E0*1*I_NAI_D2x_D2z_b;
  abcd[679] = 4.0E0*I_NAI_Dxy_G2y2z_bb-2.0E0*1*I_NAI_Dxy_D2z_b;
  abcd[680] = 4.0E0*I_NAI_Dxz_G2y2z_bb-2.0E0*1*I_NAI_Dxz_D2z_b;
  abcd[681] = 4.0E0*I_NAI_D2y_G2y2z_bb-2.0E0*1*I_NAI_D2y_D2z_b;
  abcd[682] = 4.0E0*I_NAI_Dyz_G2y2z_bb-2.0E0*1*I_NAI_Dyz_D2z_b;
  abcd[683] = 4.0E0*I_NAI_D2z_G2y2z_bb-2.0E0*1*I_NAI_D2z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[684] = 4.0E0*I_NAI_D2x_G2xyz_bb;
  abcd[685] = 4.0E0*I_NAI_Dxy_G2xyz_bb;
  abcd[686] = 4.0E0*I_NAI_Dxz_G2xyz_bb;
  abcd[687] = 4.0E0*I_NAI_D2y_G2xyz_bb;
  abcd[688] = 4.0E0*I_NAI_Dyz_G2xyz_bb;
  abcd[689] = 4.0E0*I_NAI_D2z_G2xyz_bb;
  abcd[690] = 4.0E0*I_NAI_D2x_Gx2yz_bb-2.0E0*1*I_NAI_D2x_Dxz_b;
  abcd[691] = 4.0E0*I_NAI_Dxy_Gx2yz_bb-2.0E0*1*I_NAI_Dxy_Dxz_b;
  abcd[692] = 4.0E0*I_NAI_Dxz_Gx2yz_bb-2.0E0*1*I_NAI_Dxz_Dxz_b;
  abcd[693] = 4.0E0*I_NAI_D2y_Gx2yz_bb-2.0E0*1*I_NAI_D2y_Dxz_b;
  abcd[694] = 4.0E0*I_NAI_Dyz_Gx2yz_bb-2.0E0*1*I_NAI_Dyz_Dxz_b;
  abcd[695] = 4.0E0*I_NAI_D2z_Gx2yz_bb-2.0E0*1*I_NAI_D2z_Dxz_b;
  abcd[696] = 4.0E0*I_NAI_D2x_Gxy2z_bb-2.0E0*1*I_NAI_D2x_Dxy_b;
  abcd[697] = 4.0E0*I_NAI_Dxy_Gxy2z_bb-2.0E0*1*I_NAI_Dxy_Dxy_b;
  abcd[698] = 4.0E0*I_NAI_Dxz_Gxy2z_bb-2.0E0*1*I_NAI_Dxz_Dxy_b;
  abcd[699] = 4.0E0*I_NAI_D2y_Gxy2z_bb-2.0E0*1*I_NAI_D2y_Dxy_b;
  abcd[700] = 4.0E0*I_NAI_Dyz_Gxy2z_bb-2.0E0*1*I_NAI_Dyz_Dxy_b;
  abcd[701] = 4.0E0*I_NAI_D2z_Gxy2z_bb-2.0E0*1*I_NAI_D2z_Dxy_b;
  abcd[702] = 4.0E0*I_NAI_D2x_G3yz_bb-2.0E0*2*I_NAI_D2x_Dyz_b;
  abcd[703] = 4.0E0*I_NAI_Dxy_G3yz_bb-2.0E0*2*I_NAI_Dxy_Dyz_b;
  abcd[704] = 4.0E0*I_NAI_Dxz_G3yz_bb-2.0E0*2*I_NAI_Dxz_Dyz_b;
  abcd[705] = 4.0E0*I_NAI_D2y_G3yz_bb-2.0E0*2*I_NAI_D2y_Dyz_b;
  abcd[706] = 4.0E0*I_NAI_Dyz_G3yz_bb-2.0E0*2*I_NAI_Dyz_Dyz_b;
  abcd[707] = 4.0E0*I_NAI_D2z_G3yz_bb-2.0E0*2*I_NAI_D2z_Dyz_b;
  abcd[708] = 4.0E0*I_NAI_D2x_G2y2z_bb-2.0E0*1*I_NAI_D2x_D2y_b-2.0E0*1*I_NAI_D2x_D2z_b+1*I_NAI_D2x_S;
  abcd[709] = 4.0E0*I_NAI_Dxy_G2y2z_bb-2.0E0*1*I_NAI_Dxy_D2y_b-2.0E0*1*I_NAI_Dxy_D2z_b+1*I_NAI_Dxy_S;
  abcd[710] = 4.0E0*I_NAI_Dxz_G2y2z_bb-2.0E0*1*I_NAI_Dxz_D2y_b-2.0E0*1*I_NAI_Dxz_D2z_b+1*I_NAI_Dxz_S;
  abcd[711] = 4.0E0*I_NAI_D2y_G2y2z_bb-2.0E0*1*I_NAI_D2y_D2y_b-2.0E0*1*I_NAI_D2y_D2z_b+1*I_NAI_D2y_S;
  abcd[712] = 4.0E0*I_NAI_Dyz_G2y2z_bb-2.0E0*1*I_NAI_Dyz_D2y_b-2.0E0*1*I_NAI_Dyz_D2z_b+1*I_NAI_Dyz_S;
  abcd[713] = 4.0E0*I_NAI_D2z_G2y2z_bb-2.0E0*1*I_NAI_D2z_D2y_b-2.0E0*1*I_NAI_D2z_D2z_b+1*I_NAI_D2z_S;
  abcd[714] = 4.0E0*I_NAI_D2x_Gy3z_bb-2.0E0*2*I_NAI_D2x_Dyz_b;
  abcd[715] = 4.0E0*I_NAI_Dxy_Gy3z_bb-2.0E0*2*I_NAI_Dxy_Dyz_b;
  abcd[716] = 4.0E0*I_NAI_Dxz_Gy3z_bb-2.0E0*2*I_NAI_Dxz_Dyz_b;
  abcd[717] = 4.0E0*I_NAI_D2y_Gy3z_bb-2.0E0*2*I_NAI_D2y_Dyz_b;
  abcd[718] = 4.0E0*I_NAI_Dyz_Gy3z_bb-2.0E0*2*I_NAI_Dyz_Dyz_b;
  abcd[719] = 4.0E0*I_NAI_D2z_Gy3z_bb-2.0E0*2*I_NAI_D2z_Dyz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_G_bb
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_D_b
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  abcd[720] = 4.0E0*I_NAI_D2x_G2x2z_bb-2.0E0*1*I_NAI_D2x_D2x_b;
  abcd[721] = 4.0E0*I_NAI_Dxy_G2x2z_bb-2.0E0*1*I_NAI_Dxy_D2x_b;
  abcd[722] = 4.0E0*I_NAI_Dxz_G2x2z_bb-2.0E0*1*I_NAI_Dxz_D2x_b;
  abcd[723] = 4.0E0*I_NAI_D2y_G2x2z_bb-2.0E0*1*I_NAI_D2y_D2x_b;
  abcd[724] = 4.0E0*I_NAI_Dyz_G2x2z_bb-2.0E0*1*I_NAI_Dyz_D2x_b;
  abcd[725] = 4.0E0*I_NAI_D2z_G2x2z_bb-2.0E0*1*I_NAI_D2z_D2x_b;
  abcd[726] = 4.0E0*I_NAI_D2x_Gxy2z_bb-2.0E0*1*I_NAI_D2x_Dxy_b;
  abcd[727] = 4.0E0*I_NAI_Dxy_Gxy2z_bb-2.0E0*1*I_NAI_Dxy_Dxy_b;
  abcd[728] = 4.0E0*I_NAI_Dxz_Gxy2z_bb-2.0E0*1*I_NAI_Dxz_Dxy_b;
  abcd[729] = 4.0E0*I_NAI_D2y_Gxy2z_bb-2.0E0*1*I_NAI_D2y_Dxy_b;
  abcd[730] = 4.0E0*I_NAI_Dyz_Gxy2z_bb-2.0E0*1*I_NAI_Dyz_Dxy_b;
  abcd[731] = 4.0E0*I_NAI_D2z_Gxy2z_bb-2.0E0*1*I_NAI_D2z_Dxy_b;
  abcd[732] = 4.0E0*I_NAI_D2x_Gx3z_bb-2.0E0*1*I_NAI_D2x_Dxz_b-2.0E0*2*I_NAI_D2x_Dxz_b;
  abcd[733] = 4.0E0*I_NAI_Dxy_Gx3z_bb-2.0E0*1*I_NAI_Dxy_Dxz_b-2.0E0*2*I_NAI_Dxy_Dxz_b;
  abcd[734] = 4.0E0*I_NAI_Dxz_Gx3z_bb-2.0E0*1*I_NAI_Dxz_Dxz_b-2.0E0*2*I_NAI_Dxz_Dxz_b;
  abcd[735] = 4.0E0*I_NAI_D2y_Gx3z_bb-2.0E0*1*I_NAI_D2y_Dxz_b-2.0E0*2*I_NAI_D2y_Dxz_b;
  abcd[736] = 4.0E0*I_NAI_Dyz_Gx3z_bb-2.0E0*1*I_NAI_Dyz_Dxz_b-2.0E0*2*I_NAI_Dyz_Dxz_b;
  abcd[737] = 4.0E0*I_NAI_D2z_Gx3z_bb-2.0E0*1*I_NAI_D2z_Dxz_b-2.0E0*2*I_NAI_D2z_Dxz_b;
  abcd[738] = 4.0E0*I_NAI_D2x_G2y2z_bb-2.0E0*1*I_NAI_D2x_D2y_b;
  abcd[739] = 4.0E0*I_NAI_Dxy_G2y2z_bb-2.0E0*1*I_NAI_Dxy_D2y_b;
  abcd[740] = 4.0E0*I_NAI_Dxz_G2y2z_bb-2.0E0*1*I_NAI_Dxz_D2y_b;
  abcd[741] = 4.0E0*I_NAI_D2y_G2y2z_bb-2.0E0*1*I_NAI_D2y_D2y_b;
  abcd[742] = 4.0E0*I_NAI_Dyz_G2y2z_bb-2.0E0*1*I_NAI_Dyz_D2y_b;
  abcd[743] = 4.0E0*I_NAI_D2z_G2y2z_bb-2.0E0*1*I_NAI_D2z_D2y_b;
  abcd[744] = 4.0E0*I_NAI_D2x_Gy3z_bb-2.0E0*1*I_NAI_D2x_Dyz_b-2.0E0*2*I_NAI_D2x_Dyz_b;
  abcd[745] = 4.0E0*I_NAI_Dxy_Gy3z_bb-2.0E0*1*I_NAI_Dxy_Dyz_b-2.0E0*2*I_NAI_Dxy_Dyz_b;
  abcd[746] = 4.0E0*I_NAI_Dxz_Gy3z_bb-2.0E0*1*I_NAI_Dxz_Dyz_b-2.0E0*2*I_NAI_Dxz_Dyz_b;
  abcd[747] = 4.0E0*I_NAI_D2y_Gy3z_bb-2.0E0*1*I_NAI_D2y_Dyz_b-2.0E0*2*I_NAI_D2y_Dyz_b;
  abcd[748] = 4.0E0*I_NAI_Dyz_Gy3z_bb-2.0E0*1*I_NAI_Dyz_Dyz_b-2.0E0*2*I_NAI_Dyz_Dyz_b;
  abcd[749] = 4.0E0*I_NAI_D2z_Gy3z_bb-2.0E0*1*I_NAI_D2z_Dyz_b-2.0E0*2*I_NAI_D2z_Dyz_b;
  abcd[750] = 4.0E0*I_NAI_D2x_G4z_bb-2.0E0*2*I_NAI_D2x_D2z_b-2.0E0*3*I_NAI_D2x_D2z_b+2*1*I_NAI_D2x_S;
  abcd[751] = 4.0E0*I_NAI_Dxy_G4z_bb-2.0E0*2*I_NAI_Dxy_D2z_b-2.0E0*3*I_NAI_Dxy_D2z_b+2*1*I_NAI_Dxy_S;
  abcd[752] = 4.0E0*I_NAI_Dxz_G4z_bb-2.0E0*2*I_NAI_Dxz_D2z_b-2.0E0*3*I_NAI_Dxz_D2z_b+2*1*I_NAI_Dxz_S;
  abcd[753] = 4.0E0*I_NAI_D2y_G4z_bb-2.0E0*2*I_NAI_D2y_D2z_b-2.0E0*3*I_NAI_D2y_D2z_b+2*1*I_NAI_D2y_S;
  abcd[754] = 4.0E0*I_NAI_Dyz_G4z_bb-2.0E0*2*I_NAI_Dyz_D2z_b-2.0E0*3*I_NAI_Dyz_D2z_b+2*1*I_NAI_Dyz_S;
  abcd[755] = 4.0E0*I_NAI_D2z_G4z_bb-2.0E0*2*I_NAI_D2z_D2z_b-2.0E0*3*I_NAI_D2z_D2z_b+2*1*I_NAI_D2z_S;
}
