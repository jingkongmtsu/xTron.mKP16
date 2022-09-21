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

void hgp_os_nai_sp_sp_d2(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_D2x_S_C0_aa = 0.0E0;
  Double I_NAI_Dxy_S_C0_aa = 0.0E0;
  Double I_NAI_Dxz_S_C0_aa = 0.0E0;
  Double I_NAI_D2y_S_C0_aa = 0.0E0;
  Double I_NAI_Dyz_S_C0_aa = 0.0E0;
  Double I_NAI_D2z_S_C0_aa = 0.0E0;
  Double I_NAI_S_S_C0_a = 0.0E0;
  Double I_NAI_F3x_S_C1_aa = 0.0E0;
  Double I_NAI_F2xy_S_C1_aa = 0.0E0;
  Double I_NAI_F2xz_S_C1_aa = 0.0E0;
  Double I_NAI_Fx2y_S_C1_aa = 0.0E0;
  Double I_NAI_Fxyz_S_C1_aa = 0.0E0;
  Double I_NAI_Fx2z_S_C1_aa = 0.0E0;
  Double I_NAI_F3y_S_C1_aa = 0.0E0;
  Double I_NAI_F2yz_S_C1_aa = 0.0E0;
  Double I_NAI_Fy2z_S_C1_aa = 0.0E0;
  Double I_NAI_F3z_S_C1_aa = 0.0E0;
  Double I_NAI_Px_S_C1_a = 0.0E0;
  Double I_NAI_Py_S_C1_a = 0.0E0;
  Double I_NAI_Pz_S_C1_a = 0.0E0;
  Double I_NAI_S_Px_C1000_a = 0.0E0;
  Double I_NAI_S_Py_C1000_a = 0.0E0;
  Double I_NAI_S_Pz_C1000_a = 0.0E0;
  Double I_NAI_S_Px_C1_b = 0.0E0;
  Double I_NAI_S_Py_C1_b = 0.0E0;
  Double I_NAI_S_Pz_C1_b = 0.0E0;
  Double I_NAI_Px_S_C1000_a = 0.0E0;
  Double I_NAI_Py_S_C1000_a = 0.0E0;
  Double I_NAI_Pz_S_C1000_a = 0.0E0;
  Double I_NAI_D2x_S_C1001_a = 0.0E0;
  Double I_NAI_Dxy_S_C1001_a = 0.0E0;
  Double I_NAI_Dxz_S_C1001_a = 0.0E0;
  Double I_NAI_D2y_S_C1001_a = 0.0E0;
  Double I_NAI_Dyz_S_C1001_a = 0.0E0;
  Double I_NAI_D2z_S_C1001_a = 0.0E0;
  Double I_NAI_S_D2x_C1001_b = 0.0E0;
  Double I_NAI_S_Dxy_C1001_b = 0.0E0;
  Double I_NAI_S_Dxz_C1001_b = 0.0E0;
  Double I_NAI_S_D2y_C1001_b = 0.0E0;
  Double I_NAI_S_Dyz_C1001_b = 0.0E0;
  Double I_NAI_S_D2z_C1001_b = 0.0E0;
  Double I_NAI_S_S_C1001 = 0.0E0;
  Double I_NAI_S_D2x_C0_bb = 0.0E0;
  Double I_NAI_S_Dxy_C0_bb = 0.0E0;
  Double I_NAI_S_Dxz_C0_bb = 0.0E0;
  Double I_NAI_S_D2y_C0_bb = 0.0E0;
  Double I_NAI_S_Dyz_C0_bb = 0.0E0;
  Double I_NAI_S_D2z_C0_bb = 0.0E0;
  Double I_NAI_S_S_C0_b = 0.0E0;
  Double I_NAI_Px_S_C1_b = 0.0E0;
  Double I_NAI_Py_S_C1_b = 0.0E0;
  Double I_NAI_Pz_S_C1_b = 0.0E0;
  Double I_NAI_S_F3x_C1000_bb = 0.0E0;
  Double I_NAI_S_F2xy_C1000_bb = 0.0E0;
  Double I_NAI_S_F2xz_C1000_bb = 0.0E0;
  Double I_NAI_S_Fx2y_C1000_bb = 0.0E0;
  Double I_NAI_S_Fxyz_C1000_bb = 0.0E0;
  Double I_NAI_S_Fx2z_C1000_bb = 0.0E0;
  Double I_NAI_S_F3y_C1000_bb = 0.0E0;
  Double I_NAI_S_F2yz_C1000_bb = 0.0E0;
  Double I_NAI_S_Fy2z_C1000_bb = 0.0E0;
  Double I_NAI_S_F3z_C1000_bb = 0.0E0;
  Double I_NAI_S_Px_C1000_b = 0.0E0;
  Double I_NAI_S_Py_C1000_b = 0.0E0;
  Double I_NAI_S_Pz_C1000_b = 0.0E0;
  Double I_NAI_F3x_S_C1000_aa = 0.0E0;
  Double I_NAI_F2xy_S_C1000_aa = 0.0E0;
  Double I_NAI_F2xz_S_C1000_aa = 0.0E0;
  Double I_NAI_Fx2y_S_C1000_aa = 0.0E0;
  Double I_NAI_Fxyz_S_C1000_aa = 0.0E0;
  Double I_NAI_Fx2z_S_C1000_aa = 0.0E0;
  Double I_NAI_F3y_S_C1000_aa = 0.0E0;
  Double I_NAI_F2yz_S_C1000_aa = 0.0E0;
  Double I_NAI_Fy2z_S_C1000_aa = 0.0E0;
  Double I_NAI_F3z_S_C1000_aa = 0.0E0;
  Double I_NAI_D2x_S_C1000_aa = 0.0E0;
  Double I_NAI_Dxy_S_C1000_aa = 0.0E0;
  Double I_NAI_Dxz_S_C1000_aa = 0.0E0;
  Double I_NAI_D2y_S_C1000_aa = 0.0E0;
  Double I_NAI_Dyz_S_C1000_aa = 0.0E0;
  Double I_NAI_D2z_S_C1000_aa = 0.0E0;
  Double I_NAI_G4x_S_C1001_aa = 0.0E0;
  Double I_NAI_G3xy_S_C1001_aa = 0.0E0;
  Double I_NAI_G3xz_S_C1001_aa = 0.0E0;
  Double I_NAI_G2x2y_S_C1001_aa = 0.0E0;
  Double I_NAI_G2xyz_S_C1001_aa = 0.0E0;
  Double I_NAI_G2x2z_S_C1001_aa = 0.0E0;
  Double I_NAI_Gx3y_S_C1001_aa = 0.0E0;
  Double I_NAI_Gx2yz_S_C1001_aa = 0.0E0;
  Double I_NAI_Gxy2z_S_C1001_aa = 0.0E0;
  Double I_NAI_Gx3z_S_C1001_aa = 0.0E0;
  Double I_NAI_G4y_S_C1001_aa = 0.0E0;
  Double I_NAI_G3yz_S_C1001_aa = 0.0E0;
  Double I_NAI_G2y2z_S_C1001_aa = 0.0E0;
  Double I_NAI_Gy3z_S_C1001_aa = 0.0E0;
  Double I_NAI_G4z_S_C1001_aa = 0.0E0;
  Double I_NAI_F3x_S_C1001_aa = 0.0E0;
  Double I_NAI_F2xy_S_C1001_aa = 0.0E0;
  Double I_NAI_F2xz_S_C1001_aa = 0.0E0;
  Double I_NAI_Fx2y_S_C1001_aa = 0.0E0;
  Double I_NAI_Fxyz_S_C1001_aa = 0.0E0;
  Double I_NAI_Fx2z_S_C1001_aa = 0.0E0;
  Double I_NAI_F3y_S_C1001_aa = 0.0E0;
  Double I_NAI_F2yz_S_C1001_aa = 0.0E0;
  Double I_NAI_Fy2z_S_C1001_aa = 0.0E0;
  Double I_NAI_F3z_S_C1001_aa = 0.0E0;
  Double I_NAI_Px_S_C1001_a = 0.0E0;
  Double I_NAI_Py_S_C1001_a = 0.0E0;
  Double I_NAI_Pz_S_C1001_a = 0.0E0;
  Double I_NAI_D2x_S_C0_ab = 0.0E0;
  Double I_NAI_Dxy_S_C0_ab = 0.0E0;
  Double I_NAI_Dxz_S_C0_ab = 0.0E0;
  Double I_NAI_D2y_S_C0_ab = 0.0E0;
  Double I_NAI_Dyz_S_C0_ab = 0.0E0;
  Double I_NAI_D2z_S_C0_ab = 0.0E0;
  Double I_NAI_Px_S_C0_ab = 0.0E0;
  Double I_NAI_Py_S_C0_ab = 0.0E0;
  Double I_NAI_Pz_S_C0_ab = 0.0E0;
  Double I_NAI_F3x_S_C1_ab = 0.0E0;
  Double I_NAI_F2xy_S_C1_ab = 0.0E0;
  Double I_NAI_F2xz_S_C1_ab = 0.0E0;
  Double I_NAI_Fx2y_S_C1_ab = 0.0E0;
  Double I_NAI_Fxyz_S_C1_ab = 0.0E0;
  Double I_NAI_Fx2z_S_C1_ab = 0.0E0;
  Double I_NAI_F3y_S_C1_ab = 0.0E0;
  Double I_NAI_F2yz_S_C1_ab = 0.0E0;
  Double I_NAI_Fy2z_S_C1_ab = 0.0E0;
  Double I_NAI_F3z_S_C1_ab = 0.0E0;
  Double I_NAI_D2x_S_C1_ab = 0.0E0;
  Double I_NAI_Dxy_S_C1_ab = 0.0E0;
  Double I_NAI_Dxz_S_C1_ab = 0.0E0;
  Double I_NAI_D2y_S_C1_ab = 0.0E0;
  Double I_NAI_Dyz_S_C1_ab = 0.0E0;
  Double I_NAI_D2z_S_C1_ab = 0.0E0;
  Double I_NAI_D2x_S_C1001_b = 0.0E0;
  Double I_NAI_Dxy_S_C1001_b = 0.0E0;
  Double I_NAI_Dxz_S_C1001_b = 0.0E0;
  Double I_NAI_D2y_S_C1001_b = 0.0E0;
  Double I_NAI_Dyz_S_C1001_b = 0.0E0;
  Double I_NAI_D2z_S_C1001_b = 0.0E0;
  Double I_NAI_Px_S_C1001_b = 0.0E0;
  Double I_NAI_Py_S_C1001_b = 0.0E0;
  Double I_NAI_Pz_S_C1001_b = 0.0E0;
  Double I_NAI_F3x_S_C1000_ab = 0.0E0;
  Double I_NAI_F2xy_S_C1000_ab = 0.0E0;
  Double I_NAI_F2xz_S_C1000_ab = 0.0E0;
  Double I_NAI_Fx2y_S_C1000_ab = 0.0E0;
  Double I_NAI_Fxyz_S_C1000_ab = 0.0E0;
  Double I_NAI_Fx2z_S_C1000_ab = 0.0E0;
  Double I_NAI_F3y_S_C1000_ab = 0.0E0;
  Double I_NAI_F2yz_S_C1000_ab = 0.0E0;
  Double I_NAI_Fy2z_S_C1000_ab = 0.0E0;
  Double I_NAI_F3z_S_C1000_ab = 0.0E0;
  Double I_NAI_D2x_S_C1000_ab = 0.0E0;
  Double I_NAI_Dxy_S_C1000_ab = 0.0E0;
  Double I_NAI_Dxz_S_C1000_ab = 0.0E0;
  Double I_NAI_D2y_S_C1000_ab = 0.0E0;
  Double I_NAI_Dyz_S_C1000_ab = 0.0E0;
  Double I_NAI_D2z_S_C1000_ab = 0.0E0;
  Double I_NAI_Px_S_C1000_ab = 0.0E0;
  Double I_NAI_Py_S_C1000_ab = 0.0E0;
  Double I_NAI_Pz_S_C1000_ab = 0.0E0;
  Double I_NAI_G4x_S_C1001_ab = 0.0E0;
  Double I_NAI_G3xy_S_C1001_ab = 0.0E0;
  Double I_NAI_G3xz_S_C1001_ab = 0.0E0;
  Double I_NAI_G2x2y_S_C1001_ab = 0.0E0;
  Double I_NAI_G2xyz_S_C1001_ab = 0.0E0;
  Double I_NAI_G2x2z_S_C1001_ab = 0.0E0;
  Double I_NAI_Gx3y_S_C1001_ab = 0.0E0;
  Double I_NAI_Gx2yz_S_C1001_ab = 0.0E0;
  Double I_NAI_Gxy2z_S_C1001_ab = 0.0E0;
  Double I_NAI_Gx3z_S_C1001_ab = 0.0E0;
  Double I_NAI_G4y_S_C1001_ab = 0.0E0;
  Double I_NAI_G3yz_S_C1001_ab = 0.0E0;
  Double I_NAI_G2y2z_S_C1001_ab = 0.0E0;
  Double I_NAI_Gy3z_S_C1001_ab = 0.0E0;
  Double I_NAI_G4z_S_C1001_ab = 0.0E0;
  Double I_NAI_F3x_S_C1001_ab = 0.0E0;
  Double I_NAI_F2xy_S_C1001_ab = 0.0E0;
  Double I_NAI_F2xz_S_C1001_ab = 0.0E0;
  Double I_NAI_Fx2y_S_C1001_ab = 0.0E0;
  Double I_NAI_Fxyz_S_C1001_ab = 0.0E0;
  Double I_NAI_Fx2z_S_C1001_ab = 0.0E0;
  Double I_NAI_F3y_S_C1001_ab = 0.0E0;
  Double I_NAI_F2yz_S_C1001_ab = 0.0E0;
  Double I_NAI_Fy2z_S_C1001_ab = 0.0E0;
  Double I_NAI_F3z_S_C1001_ab = 0.0E0;
  Double I_NAI_D2x_S_C1001_ab = 0.0E0;
  Double I_NAI_Dxy_S_C1001_ab = 0.0E0;
  Double I_NAI_Dxz_S_C1001_ab = 0.0E0;
  Double I_NAI_D2y_S_C1001_ab = 0.0E0;
  Double I_NAI_Dyz_S_C1001_ab = 0.0E0;
  Double I_NAI_D2z_S_C1001_ab = 0.0E0;
  Double I_NAI_F3x_S_C1_bb = 0.0E0;
  Double I_NAI_F2xy_S_C1_bb = 0.0E0;
  Double I_NAI_F2xz_S_C1_bb = 0.0E0;
  Double I_NAI_Fx2y_S_C1_bb = 0.0E0;
  Double I_NAI_Fxyz_S_C1_bb = 0.0E0;
  Double I_NAI_Fx2z_S_C1_bb = 0.0E0;
  Double I_NAI_F3y_S_C1_bb = 0.0E0;
  Double I_NAI_F2yz_S_C1_bb = 0.0E0;
  Double I_NAI_Fy2z_S_C1_bb = 0.0E0;
  Double I_NAI_F3z_S_C1_bb = 0.0E0;
  Double I_NAI_D2x_S_C1_bb = 0.0E0;
  Double I_NAI_Dxy_S_C1_bb = 0.0E0;
  Double I_NAI_Dxz_S_C1_bb = 0.0E0;
  Double I_NAI_D2y_S_C1_bb = 0.0E0;
  Double I_NAI_Dyz_S_C1_bb = 0.0E0;
  Double I_NAI_D2z_S_C1_bb = 0.0E0;
  Double I_NAI_Px_S_C1_bb = 0.0E0;
  Double I_NAI_Py_S_C1_bb = 0.0E0;
  Double I_NAI_Pz_S_C1_bb = 0.0E0;
  Double I_NAI_G4x_S_C1001_bb = 0.0E0;
  Double I_NAI_G3xy_S_C1001_bb = 0.0E0;
  Double I_NAI_G3xz_S_C1001_bb = 0.0E0;
  Double I_NAI_G2x2y_S_C1001_bb = 0.0E0;
  Double I_NAI_G2xyz_S_C1001_bb = 0.0E0;
  Double I_NAI_G2x2z_S_C1001_bb = 0.0E0;
  Double I_NAI_Gx3y_S_C1001_bb = 0.0E0;
  Double I_NAI_Gx2yz_S_C1001_bb = 0.0E0;
  Double I_NAI_Gxy2z_S_C1001_bb = 0.0E0;
  Double I_NAI_Gx3z_S_C1001_bb = 0.0E0;
  Double I_NAI_G4y_S_C1001_bb = 0.0E0;
  Double I_NAI_G3yz_S_C1001_bb = 0.0E0;
  Double I_NAI_G2y2z_S_C1001_bb = 0.0E0;
  Double I_NAI_Gy3z_S_C1001_bb = 0.0E0;
  Double I_NAI_G4z_S_C1001_bb = 0.0E0;
  Double I_NAI_F3x_S_C1001_bb = 0.0E0;
  Double I_NAI_F2xy_S_C1001_bb = 0.0E0;
  Double I_NAI_F2xz_S_C1001_bb = 0.0E0;
  Double I_NAI_Fx2y_S_C1001_bb = 0.0E0;
  Double I_NAI_Fxyz_S_C1001_bb = 0.0E0;
  Double I_NAI_Fx2z_S_C1001_bb = 0.0E0;
  Double I_NAI_F3y_S_C1001_bb = 0.0E0;
  Double I_NAI_F2yz_S_C1001_bb = 0.0E0;
  Double I_NAI_Fy2z_S_C1001_bb = 0.0E0;
  Double I_NAI_F3z_S_C1001_bb = 0.0E0;
  Double I_NAI_D2x_S_C1001_bb = 0.0E0;
  Double I_NAI_Dxy_S_C1001_bb = 0.0E0;
  Double I_NAI_Dxz_S_C1001_bb = 0.0E0;
  Double I_NAI_D2y_S_C1001_bb = 0.0E0;
  Double I_NAI_Dyz_S_C1001_bb = 0.0E0;
  Double I_NAI_D2z_S_C1001_bb = 0.0E0;
  Double I_NAI_Px_S_C1001_bb = 0.0E0;
  Double I_NAI_Py_S_C1001_bb = 0.0E0;
  Double I_NAI_Pz_S_C1001_bb = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
    Double ic2_2 = icoe[ip2+2*inp2];
    Double ic2_3 = icoe[ip2+3*inp2];
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

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER43;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER15*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER13*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER11*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = ONEOVER9*I_NAI_S_S_M4_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M4_vrr  = f*I_NAI_S_S_M4_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);

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

#endif

      }


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
       * shell quartet name: SQ_NAI_S_P_M2
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_S_Px_M2_vrr = PBX*I_NAI_S_S_M2_vrr-PNX*I_NAI_S_S_M3_vrr;
      Double I_NAI_S_Py_M2_vrr = PBY*I_NAI_S_S_M2_vrr-PNY*I_NAI_S_S_M3_vrr;
      Double I_NAI_S_Pz_M2_vrr = PBZ*I_NAI_S_S_M2_vrr-PNZ*I_NAI_S_S_M3_vrr;

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
       * shell quartet name: SQ_NAI_S_D_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_P_M1
       * RHS shell quartet name: SQ_NAI_S_P_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_S_D2x_M1_vrr = PBX*I_NAI_S_Px_M1_vrr-PNX*I_NAI_S_Px_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_S_Dxy_M1_vrr = PBY*I_NAI_S_Px_M1_vrr-PNY*I_NAI_S_Px_M2_vrr;
      Double I_NAI_S_D2y_M1_vrr = PBY*I_NAI_S_Py_M1_vrr-PNY*I_NAI_S_Py_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_S_D2z_M1_vrr = PBZ*I_NAI_S_Pz_M1_vrr-PNZ*I_NAI_S_Pz_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;

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
       * shell quartet name: SQ_NAI_S_F
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_D
       * RHS shell quartet name: SQ_NAI_S_D_M1
       * RHS shell quartet name: SQ_NAI_S_P
       * RHS shell quartet name: SQ_NAI_S_P_M1
       ************************************************************/
      Double I_NAI_S_F3x_vrr = PBX*I_NAI_S_D2x_vrr-PNX*I_NAI_S_D2x_M1_vrr+2*oned2z*I_NAI_S_Px_vrr-2*oned2z*I_NAI_S_Px_M1_vrr;
      Double I_NAI_S_F2xy_vrr = PBY*I_NAI_S_D2x_vrr-PNY*I_NAI_S_D2x_M1_vrr;
      Double I_NAI_S_F2xz_vrr = PBZ*I_NAI_S_D2x_vrr-PNZ*I_NAI_S_D2x_M1_vrr;
      Double I_NAI_S_Fx2y_vrr = PBX*I_NAI_S_D2y_vrr-PNX*I_NAI_S_D2y_M1_vrr;
      Double I_NAI_S_Fxyz_vrr = PBZ*I_NAI_S_Dxy_vrr-PNZ*I_NAI_S_Dxy_M1_vrr;
      Double I_NAI_S_Fx2z_vrr = PBX*I_NAI_S_D2z_vrr-PNX*I_NAI_S_D2z_M1_vrr;
      Double I_NAI_S_F3y_vrr = PBY*I_NAI_S_D2y_vrr-PNY*I_NAI_S_D2y_M1_vrr+2*oned2z*I_NAI_S_Py_vrr-2*oned2z*I_NAI_S_Py_M1_vrr;
      Double I_NAI_S_F2yz_vrr = PBZ*I_NAI_S_D2y_vrr-PNZ*I_NAI_S_D2y_M1_vrr;
      Double I_NAI_S_Fy2z_vrr = PBY*I_NAI_S_D2z_vrr-PNY*I_NAI_S_D2z_M1_vrr;
      Double I_NAI_S_F3z_vrr = PBZ*I_NAI_S_D2z_vrr-PNZ*I_NAI_S_D2z_M1_vrr+2*oned2z*I_NAI_S_Pz_vrr-2*oned2z*I_NAI_S_Pz_M1_vrr;

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
       * shell quartet name: SQ_NAI_D_S_C0_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C0_aa_coefs = ic2*alpha*alpha;
      I_NAI_D2x_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C0_aa += SQ_NAI_D_S_C0_aa_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_S_C0_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_S_C0_a_coefs = ic2*alpha;
      I_NAI_S_S_C0_a += SQ_NAI_S_S_C0_a_coefs*I_NAI_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1_aa_coefs = ic2_1*alpha*alpha;
      I_NAI_F3x_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1_aa += SQ_NAI_F_S_C1_aa_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1_a_coefs = ic2_1*alpha;
      I_NAI_Px_S_C1_a += SQ_NAI_P_S_C1_a_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1_a += SQ_NAI_P_S_C1_a_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1_a += SQ_NAI_P_S_C1_a_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_P_C1000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_P_C1000_a_coefs = ic2_2*alpha;
      I_NAI_S_Px_C1000_a += SQ_NAI_S_P_C1000_a_coefs*I_NAI_S_Px_vrr;
      I_NAI_S_Py_C1000_a += SQ_NAI_S_P_C1000_a_coefs*I_NAI_S_Py_vrr;
      I_NAI_S_Pz_C1000_a += SQ_NAI_S_P_C1000_a_coefs*I_NAI_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_P_C1_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_P_C1_b_coefs = ic2_1*beta;
      I_NAI_S_Px_C1_b += SQ_NAI_S_P_C1_b_coefs*I_NAI_S_Px_vrr;
      I_NAI_S_Py_C1_b += SQ_NAI_S_P_C1_b_coefs*I_NAI_S_Py_vrr;
      I_NAI_S_Pz_C1_b += SQ_NAI_S_P_C1_b_coefs*I_NAI_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1000_a_coefs = ic2_2*alpha;
      I_NAI_Px_S_C1000_a += SQ_NAI_P_S_C1000_a_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1000_a += SQ_NAI_P_S_C1000_a_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1000_a += SQ_NAI_P_S_C1000_a_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1001_a_coefs = ic2_3*alpha;
      I_NAI_D2x_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1001_a += SQ_NAI_D_S_C1001_a_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_D_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_D_C1001_b_coefs = ic2_3*beta;
      I_NAI_S_D2x_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_D2x_vrr;
      I_NAI_S_Dxy_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_Dxy_vrr;
      I_NAI_S_Dxz_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_Dxz_vrr;
      I_NAI_S_D2y_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_D2y_vrr;
      I_NAI_S_Dyz_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_Dyz_vrr;
      I_NAI_S_D2z_C1001_b += SQ_NAI_S_D_C1001_b_coefs*I_NAI_S_D2z_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_S_C1001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_S_C1001_coefs = ic2_3;
      I_NAI_S_S_C1001 += SQ_NAI_S_S_C1001_coefs*I_NAI_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_D_C0_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_D_C0_bb_coefs = ic2*beta*beta;
      I_NAI_S_D2x_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_D2x_vrr;
      I_NAI_S_Dxy_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_Dxy_vrr;
      I_NAI_S_Dxz_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_Dxz_vrr;
      I_NAI_S_D2y_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_D2y_vrr;
      I_NAI_S_Dyz_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_Dyz_vrr;
      I_NAI_S_D2z_C0_bb += SQ_NAI_S_D_C0_bb_coefs*I_NAI_S_D2z_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_S_C0_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_S_C0_b_coefs = ic2*beta;
      I_NAI_S_S_C0_b += SQ_NAI_S_S_C0_b_coefs*I_NAI_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1_b_coefs = ic2_1*beta;
      I_NAI_Px_S_C1_b += SQ_NAI_P_S_C1_b_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1_b += SQ_NAI_P_S_C1_b_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1_b += SQ_NAI_P_S_C1_b_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_F_C1000_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_F_C1000_bb_coefs = ic2_2*beta*beta;
      I_NAI_S_F3x_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F3x_vrr;
      I_NAI_S_F2xy_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F2xy_vrr;
      I_NAI_S_F2xz_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F2xz_vrr;
      I_NAI_S_Fx2y_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_Fx2y_vrr;
      I_NAI_S_Fxyz_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_Fxyz_vrr;
      I_NAI_S_Fx2z_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_Fx2z_vrr;
      I_NAI_S_F3y_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F3y_vrr;
      I_NAI_S_F2yz_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F2yz_vrr;
      I_NAI_S_Fy2z_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_Fy2z_vrr;
      I_NAI_S_F3z_C1000_bb += SQ_NAI_S_F_C1000_bb_coefs*I_NAI_S_F3z_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_S_P_C1000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_S_P_C1000_b_coefs = ic2_2*beta;
      I_NAI_S_Px_C1000_b += SQ_NAI_S_P_C1000_b_coefs*I_NAI_S_Px_vrr;
      I_NAI_S_Py_C1000_b += SQ_NAI_S_P_C1000_b_coefs*I_NAI_S_Py_vrr;
      I_NAI_S_Pz_C1000_b += SQ_NAI_S_P_C1000_b_coefs*I_NAI_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1000_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1000_aa_coefs = ic2_2*alpha*alpha;
      I_NAI_F3x_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1000_aa += SQ_NAI_F_S_C1000_aa_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1000_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1000_aa_coefs = ic2_2*alpha*alpha;
      I_NAI_D2x_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1000_aa += SQ_NAI_D_S_C1000_aa_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1001_aa_coefs = ic2_3*alpha*alpha;
      I_NAI_G4x_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1001_aa += SQ_NAI_G_S_C1001_aa_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1001_aa_coefs = ic2_3*alpha*alpha;
      I_NAI_F3x_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1001_aa += SQ_NAI_F_S_C1001_aa_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1001_a_coefs = ic2_3*alpha;
      I_NAI_Px_S_C1001_a += SQ_NAI_P_S_C1001_a_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1001_a += SQ_NAI_P_S_C1001_a_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1001_a += SQ_NAI_P_S_C1001_a_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C0_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C0_ab_coefs = ic2*alpha*beta;
      I_NAI_D2x_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C0_ab += SQ_NAI_D_S_C0_ab_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C0_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C0_ab_coefs = ic2*alpha*beta;
      I_NAI_Px_S_C0_ab += SQ_NAI_P_S_C0_ab_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C0_ab += SQ_NAI_P_S_C0_ab_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C0_ab += SQ_NAI_P_S_C0_ab_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1_ab_coefs = ic2_1*alpha*beta;
      I_NAI_F3x_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1_ab += SQ_NAI_F_S_C1_ab_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1_ab_coefs = ic2_1*alpha*beta;
      I_NAI_D2x_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1_ab += SQ_NAI_D_S_C1_ab_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1001_b_coefs = ic2_3*beta;
      I_NAI_D2x_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1001_b += SQ_NAI_D_S_C1001_b_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1001_b_coefs = ic2_3*beta;
      I_NAI_Px_S_C1001_b += SQ_NAI_P_S_C1001_b_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1001_b += SQ_NAI_P_S_C1001_b_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1001_b += SQ_NAI_P_S_C1001_b_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1000_ab_coefs = ic2_2*alpha*beta;
      I_NAI_F3x_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1000_ab += SQ_NAI_F_S_C1000_ab_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1000_ab_coefs = ic2_2*alpha*beta;
      I_NAI_D2x_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1000_ab += SQ_NAI_D_S_C1000_ab_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1000_ab_coefs = ic2_2*alpha*beta;
      I_NAI_Px_S_C1000_ab += SQ_NAI_P_S_C1000_ab_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1000_ab += SQ_NAI_P_S_C1000_ab_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1000_ab += SQ_NAI_P_S_C1000_ab_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1001_ab_coefs = ic2_3*alpha*beta;
      I_NAI_G4x_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1001_ab += SQ_NAI_G_S_C1001_ab_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1001_ab_coefs = ic2_3*alpha*beta;
      I_NAI_F3x_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1001_ab += SQ_NAI_F_S_C1001_ab_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1001_ab_coefs = ic2_3*alpha*beta;
      I_NAI_D2x_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1001_ab += SQ_NAI_D_S_C1001_ab_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1_bb_coefs = ic2_1*beta*beta;
      I_NAI_F3x_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1_bb += SQ_NAI_F_S_C1_bb_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1_bb_coefs = ic2_1*beta*beta;
      I_NAI_D2x_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1_bb += SQ_NAI_D_S_C1_bb_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1_bb_coefs = ic2_1*beta*beta;
      I_NAI_Px_S_C1_bb += SQ_NAI_P_S_C1_bb_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1_bb += SQ_NAI_P_S_C1_bb_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1_bb += SQ_NAI_P_S_C1_bb_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1001_bb_coefs = ic2_3*beta*beta;
      I_NAI_G4x_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1001_bb += SQ_NAI_G_S_C1001_bb_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1001_bb_coefs = ic2_3*beta*beta;
      I_NAI_F3x_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1001_bb += SQ_NAI_F_S_C1001_bb_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1001_bb_coefs = ic2_3*beta*beta;
      I_NAI_D2x_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1001_bb += SQ_NAI_D_S_C1001_bb_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1001_bb_coefs = ic2_3*beta*beta;
      I_NAI_Px_S_C1001_bb += SQ_NAI_P_S_C1001_bb_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1001_bb += SQ_NAI_P_S_C1001_bb_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1001_bb += SQ_NAI_P_S_C1001_bb_coefs*I_NAI_Pz_S_vrr;
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
   * shell quartet name: SQ_NAI_P_P_C1001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_P_S_C1001_a
   ************************************************************/
  Double I_NAI_Px_Px_C1001_a = I_NAI_D2x_S_C1001_a+ABX*I_NAI_Px_S_C1001_a;
  Double I_NAI_Py_Px_C1001_a = I_NAI_Dxy_S_C1001_a+ABX*I_NAI_Py_S_C1001_a;
  Double I_NAI_Pz_Px_C1001_a = I_NAI_Dxz_S_C1001_a+ABX*I_NAI_Pz_S_C1001_a;
  Double I_NAI_Px_Py_C1001_a = I_NAI_Dxy_S_C1001_a+ABY*I_NAI_Px_S_C1001_a;
  Double I_NAI_Py_Py_C1001_a = I_NAI_D2y_S_C1001_a+ABY*I_NAI_Py_S_C1001_a;
  Double I_NAI_Pz_Py_C1001_a = I_NAI_Dyz_S_C1001_a+ABY*I_NAI_Pz_S_C1001_a;
  Double I_NAI_Px_Pz_C1001_a = I_NAI_Dxz_S_C1001_a+ABZ*I_NAI_Px_S_C1001_a;
  Double I_NAI_Py_Pz_C1001_a = I_NAI_Dyz_S_C1001_a+ABZ*I_NAI_Py_S_C1001_a;
  Double I_NAI_Pz_Pz_C1001_a = I_NAI_D2z_S_C1001_a+ABZ*I_NAI_Pz_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1001_b
   * RHS shell quartet name: SQ_NAI_P_S_C1001_b
   ************************************************************/
  Double I_NAI_Px_Px_C1001_b = I_NAI_D2x_S_C1001_b+ABX*I_NAI_Px_S_C1001_b;
  Double I_NAI_Py_Px_C1001_b = I_NAI_Dxy_S_C1001_b+ABX*I_NAI_Py_S_C1001_b;
  Double I_NAI_Pz_Px_C1001_b = I_NAI_Dxz_S_C1001_b+ABX*I_NAI_Pz_S_C1001_b;
  Double I_NAI_Px_Py_C1001_b = I_NAI_Dxy_S_C1001_b+ABY*I_NAI_Px_S_C1001_b;
  Double I_NAI_Py_Py_C1001_b = I_NAI_D2y_S_C1001_b+ABY*I_NAI_Py_S_C1001_b;
  Double I_NAI_Pz_Py_C1001_b = I_NAI_Dyz_S_C1001_b+ABY*I_NAI_Pz_S_C1001_b;
  Double I_NAI_Px_Pz_C1001_b = I_NAI_Dxz_S_C1001_b+ABZ*I_NAI_Px_S_C1001_b;
  Double I_NAI_Py_Pz_C1001_b = I_NAI_Dyz_S_C1001_b+ABZ*I_NAI_Py_S_C1001_b;
  Double I_NAI_Pz_Pz_C1001_b = I_NAI_D2z_S_C1001_b+ABZ*I_NAI_Pz_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1000_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1000_aa
   * RHS shell quartet name: SQ_NAI_D_S_C1000_aa
   ************************************************************/
  Double I_NAI_D2x_Px_C1000_aa = I_NAI_F3x_S_C1000_aa+ABX*I_NAI_D2x_S_C1000_aa;
  Double I_NAI_Dxy_Px_C1000_aa = I_NAI_F2xy_S_C1000_aa+ABX*I_NAI_Dxy_S_C1000_aa;
  Double I_NAI_Dxz_Px_C1000_aa = I_NAI_F2xz_S_C1000_aa+ABX*I_NAI_Dxz_S_C1000_aa;
  Double I_NAI_D2y_Px_C1000_aa = I_NAI_Fx2y_S_C1000_aa+ABX*I_NAI_D2y_S_C1000_aa;
  Double I_NAI_Dyz_Px_C1000_aa = I_NAI_Fxyz_S_C1000_aa+ABX*I_NAI_Dyz_S_C1000_aa;
  Double I_NAI_D2z_Px_C1000_aa = I_NAI_Fx2z_S_C1000_aa+ABX*I_NAI_D2z_S_C1000_aa;
  Double I_NAI_D2x_Py_C1000_aa = I_NAI_F2xy_S_C1000_aa+ABY*I_NAI_D2x_S_C1000_aa;
  Double I_NAI_Dxy_Py_C1000_aa = I_NAI_Fx2y_S_C1000_aa+ABY*I_NAI_Dxy_S_C1000_aa;
  Double I_NAI_Dxz_Py_C1000_aa = I_NAI_Fxyz_S_C1000_aa+ABY*I_NAI_Dxz_S_C1000_aa;
  Double I_NAI_D2y_Py_C1000_aa = I_NAI_F3y_S_C1000_aa+ABY*I_NAI_D2y_S_C1000_aa;
  Double I_NAI_Dyz_Py_C1000_aa = I_NAI_F2yz_S_C1000_aa+ABY*I_NAI_Dyz_S_C1000_aa;
  Double I_NAI_D2z_Py_C1000_aa = I_NAI_Fy2z_S_C1000_aa+ABY*I_NAI_D2z_S_C1000_aa;
  Double I_NAI_D2x_Pz_C1000_aa = I_NAI_F2xz_S_C1000_aa+ABZ*I_NAI_D2x_S_C1000_aa;
  Double I_NAI_Dxy_Pz_C1000_aa = I_NAI_Fxyz_S_C1000_aa+ABZ*I_NAI_Dxy_S_C1000_aa;
  Double I_NAI_Dxz_Pz_C1000_aa = I_NAI_Fx2z_S_C1000_aa+ABZ*I_NAI_Dxz_S_C1000_aa;
  Double I_NAI_D2y_Pz_C1000_aa = I_NAI_F2yz_S_C1000_aa+ABZ*I_NAI_D2y_S_C1000_aa;
  Double I_NAI_Dyz_Pz_C1000_aa = I_NAI_Fy2z_S_C1000_aa+ABZ*I_NAI_Dyz_S_C1000_aa;
  Double I_NAI_D2z_Pz_C1000_aa = I_NAI_F3z_S_C1000_aa+ABZ*I_NAI_D2z_S_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1001_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1001_aa
   * RHS shell quartet name: SQ_NAI_F_S_C1001_aa
   ************************************************************/
  Double I_NAI_F3x_Px_C1001_aa = I_NAI_G4x_S_C1001_aa+ABX*I_NAI_F3x_S_C1001_aa;
  Double I_NAI_F2xy_Px_C1001_aa = I_NAI_G3xy_S_C1001_aa+ABX*I_NAI_F2xy_S_C1001_aa;
  Double I_NAI_F2xz_Px_C1001_aa = I_NAI_G3xz_S_C1001_aa+ABX*I_NAI_F2xz_S_C1001_aa;
  Double I_NAI_Fx2y_Px_C1001_aa = I_NAI_G2x2y_S_C1001_aa+ABX*I_NAI_Fx2y_S_C1001_aa;
  Double I_NAI_Fxyz_Px_C1001_aa = I_NAI_G2xyz_S_C1001_aa+ABX*I_NAI_Fxyz_S_C1001_aa;
  Double I_NAI_Fx2z_Px_C1001_aa = I_NAI_G2x2z_S_C1001_aa+ABX*I_NAI_Fx2z_S_C1001_aa;
  Double I_NAI_F3y_Px_C1001_aa = I_NAI_Gx3y_S_C1001_aa+ABX*I_NAI_F3y_S_C1001_aa;
  Double I_NAI_F2yz_Px_C1001_aa = I_NAI_Gx2yz_S_C1001_aa+ABX*I_NAI_F2yz_S_C1001_aa;
  Double I_NAI_Fy2z_Px_C1001_aa = I_NAI_Gxy2z_S_C1001_aa+ABX*I_NAI_Fy2z_S_C1001_aa;
  Double I_NAI_F3z_Px_C1001_aa = I_NAI_Gx3z_S_C1001_aa+ABX*I_NAI_F3z_S_C1001_aa;
  Double I_NAI_F3x_Py_C1001_aa = I_NAI_G3xy_S_C1001_aa+ABY*I_NAI_F3x_S_C1001_aa;
  Double I_NAI_F2xy_Py_C1001_aa = I_NAI_G2x2y_S_C1001_aa+ABY*I_NAI_F2xy_S_C1001_aa;
  Double I_NAI_F2xz_Py_C1001_aa = I_NAI_G2xyz_S_C1001_aa+ABY*I_NAI_F2xz_S_C1001_aa;
  Double I_NAI_Fx2y_Py_C1001_aa = I_NAI_Gx3y_S_C1001_aa+ABY*I_NAI_Fx2y_S_C1001_aa;
  Double I_NAI_Fxyz_Py_C1001_aa = I_NAI_Gx2yz_S_C1001_aa+ABY*I_NAI_Fxyz_S_C1001_aa;
  Double I_NAI_Fx2z_Py_C1001_aa = I_NAI_Gxy2z_S_C1001_aa+ABY*I_NAI_Fx2z_S_C1001_aa;
  Double I_NAI_F3y_Py_C1001_aa = I_NAI_G4y_S_C1001_aa+ABY*I_NAI_F3y_S_C1001_aa;
  Double I_NAI_F2yz_Py_C1001_aa = I_NAI_G3yz_S_C1001_aa+ABY*I_NAI_F2yz_S_C1001_aa;
  Double I_NAI_Fy2z_Py_C1001_aa = I_NAI_G2y2z_S_C1001_aa+ABY*I_NAI_Fy2z_S_C1001_aa;
  Double I_NAI_F3z_Py_C1001_aa = I_NAI_Gy3z_S_C1001_aa+ABY*I_NAI_F3z_S_C1001_aa;
  Double I_NAI_F3x_Pz_C1001_aa = I_NAI_G3xz_S_C1001_aa+ABZ*I_NAI_F3x_S_C1001_aa;
  Double I_NAI_F2xy_Pz_C1001_aa = I_NAI_G2xyz_S_C1001_aa+ABZ*I_NAI_F2xy_S_C1001_aa;
  Double I_NAI_F2xz_Pz_C1001_aa = I_NAI_G2x2z_S_C1001_aa+ABZ*I_NAI_F2xz_S_C1001_aa;
  Double I_NAI_Fx2y_Pz_C1001_aa = I_NAI_Gx2yz_S_C1001_aa+ABZ*I_NAI_Fx2y_S_C1001_aa;
  Double I_NAI_Fxyz_Pz_C1001_aa = I_NAI_Gxy2z_S_C1001_aa+ABZ*I_NAI_Fxyz_S_C1001_aa;
  Double I_NAI_Fx2z_Pz_C1001_aa = I_NAI_Gx3z_S_C1001_aa+ABZ*I_NAI_Fx2z_S_C1001_aa;
  Double I_NAI_F3y_Pz_C1001_aa = I_NAI_G3yz_S_C1001_aa+ABZ*I_NAI_F3y_S_C1001_aa;
  Double I_NAI_F2yz_Pz_C1001_aa = I_NAI_G2y2z_S_C1001_aa+ABZ*I_NAI_F2yz_S_C1001_aa;
  Double I_NAI_Fy2z_Pz_C1001_aa = I_NAI_Gy3z_S_C1001_aa+ABZ*I_NAI_Fy2z_S_C1001_aa;
  Double I_NAI_F3z_Pz_C1001_aa = I_NAI_G4z_S_C1001_aa+ABZ*I_NAI_F3z_S_C1001_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C0_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_ab
   * RHS shell quartet name: SQ_NAI_P_S_C0_ab
   ************************************************************/
  Double I_NAI_Px_Px_C0_ab = I_NAI_D2x_S_C0_ab+ABX*I_NAI_Px_S_C0_ab;
  Double I_NAI_Py_Px_C0_ab = I_NAI_Dxy_S_C0_ab+ABX*I_NAI_Py_S_C0_ab;
  Double I_NAI_Pz_Px_C0_ab = I_NAI_Dxz_S_C0_ab+ABX*I_NAI_Pz_S_C0_ab;
  Double I_NAI_Px_Py_C0_ab = I_NAI_Dxy_S_C0_ab+ABY*I_NAI_Px_S_C0_ab;
  Double I_NAI_Py_Py_C0_ab = I_NAI_D2y_S_C0_ab+ABY*I_NAI_Py_S_C0_ab;
  Double I_NAI_Pz_Py_C0_ab = I_NAI_Dyz_S_C0_ab+ABY*I_NAI_Pz_S_C0_ab;
  Double I_NAI_Px_Pz_C0_ab = I_NAI_Dxz_S_C0_ab+ABZ*I_NAI_Px_S_C0_ab;
  Double I_NAI_Py_Pz_C0_ab = I_NAI_Dyz_S_C0_ab+ABZ*I_NAI_Py_S_C0_ab;
  Double I_NAI_Pz_Pz_C0_ab = I_NAI_D2z_S_C0_ab+ABZ*I_NAI_Pz_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1_ab
   ************************************************************/
  Double I_NAI_D2x_Px_C1_ab = I_NAI_F3x_S_C1_ab+ABX*I_NAI_D2x_S_C1_ab;
  Double I_NAI_Dxy_Px_C1_ab = I_NAI_F2xy_S_C1_ab+ABX*I_NAI_Dxy_S_C1_ab;
  Double I_NAI_Dxz_Px_C1_ab = I_NAI_F2xz_S_C1_ab+ABX*I_NAI_Dxz_S_C1_ab;
  Double I_NAI_D2y_Px_C1_ab = I_NAI_Fx2y_S_C1_ab+ABX*I_NAI_D2y_S_C1_ab;
  Double I_NAI_Dyz_Px_C1_ab = I_NAI_Fxyz_S_C1_ab+ABX*I_NAI_Dyz_S_C1_ab;
  Double I_NAI_D2z_Px_C1_ab = I_NAI_Fx2z_S_C1_ab+ABX*I_NAI_D2z_S_C1_ab;
  Double I_NAI_D2x_Py_C1_ab = I_NAI_F2xy_S_C1_ab+ABY*I_NAI_D2x_S_C1_ab;
  Double I_NAI_Dxy_Py_C1_ab = I_NAI_Fx2y_S_C1_ab+ABY*I_NAI_Dxy_S_C1_ab;
  Double I_NAI_Dxz_Py_C1_ab = I_NAI_Fxyz_S_C1_ab+ABY*I_NAI_Dxz_S_C1_ab;
  Double I_NAI_D2y_Py_C1_ab = I_NAI_F3y_S_C1_ab+ABY*I_NAI_D2y_S_C1_ab;
  Double I_NAI_Dyz_Py_C1_ab = I_NAI_F2yz_S_C1_ab+ABY*I_NAI_Dyz_S_C1_ab;
  Double I_NAI_D2z_Py_C1_ab = I_NAI_Fy2z_S_C1_ab+ABY*I_NAI_D2z_S_C1_ab;
  Double I_NAI_D2x_Pz_C1_ab = I_NAI_F2xz_S_C1_ab+ABZ*I_NAI_D2x_S_C1_ab;
  Double I_NAI_Dxy_Pz_C1_ab = I_NAI_Fxyz_S_C1_ab+ABZ*I_NAI_Dxy_S_C1_ab;
  Double I_NAI_Dxz_Pz_C1_ab = I_NAI_Fx2z_S_C1_ab+ABZ*I_NAI_Dxz_S_C1_ab;
  Double I_NAI_D2y_Pz_C1_ab = I_NAI_F2yz_S_C1_ab+ABZ*I_NAI_D2y_S_C1_ab;
  Double I_NAI_Dyz_Pz_C1_ab = I_NAI_Fy2z_S_C1_ab+ABZ*I_NAI_Dyz_S_C1_ab;
  Double I_NAI_D2z_Pz_C1_ab = I_NAI_F3z_S_C1_ab+ABZ*I_NAI_D2z_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_ab
   ************************************************************/
  Double I_NAI_Px_Px_C1000_ab = I_NAI_D2x_S_C1000_ab+ABX*I_NAI_Px_S_C1000_ab;
  Double I_NAI_Py_Px_C1000_ab = I_NAI_Dxy_S_C1000_ab+ABX*I_NAI_Py_S_C1000_ab;
  Double I_NAI_Pz_Px_C1000_ab = I_NAI_Dxz_S_C1000_ab+ABX*I_NAI_Pz_S_C1000_ab;
  Double I_NAI_Px_Py_C1000_ab = I_NAI_Dxy_S_C1000_ab+ABY*I_NAI_Px_S_C1000_ab;
  Double I_NAI_Py_Py_C1000_ab = I_NAI_D2y_S_C1000_ab+ABY*I_NAI_Py_S_C1000_ab;
  Double I_NAI_Pz_Py_C1000_ab = I_NAI_Dyz_S_C1000_ab+ABY*I_NAI_Pz_S_C1000_ab;
  Double I_NAI_Px_Pz_C1000_ab = I_NAI_Dxz_S_C1000_ab+ABZ*I_NAI_Px_S_C1000_ab;
  Double I_NAI_Py_Pz_C1000_ab = I_NAI_Dyz_S_C1000_ab+ABZ*I_NAI_Py_S_C1000_ab;
  Double I_NAI_Pz_Pz_C1000_ab = I_NAI_D2z_S_C1000_ab+ABZ*I_NAI_Pz_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 4 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1000_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1000_ab
   ************************************************************/
  Double I_NAI_D2x_Px_C1000_ab = I_NAI_F3x_S_C1000_ab+ABX*I_NAI_D2x_S_C1000_ab;
  Double I_NAI_Dxy_Px_C1000_ab = I_NAI_F2xy_S_C1000_ab+ABX*I_NAI_Dxy_S_C1000_ab;
  Double I_NAI_Dxz_Px_C1000_ab = I_NAI_F2xz_S_C1000_ab+ABX*I_NAI_Dxz_S_C1000_ab;
  Double I_NAI_D2y_Px_C1000_ab = I_NAI_Fx2y_S_C1000_ab+ABX*I_NAI_D2y_S_C1000_ab;
  Double I_NAI_Dyz_Px_C1000_ab = I_NAI_Fxyz_S_C1000_ab+ABX*I_NAI_Dyz_S_C1000_ab;
  Double I_NAI_D2z_Px_C1000_ab = I_NAI_Fx2z_S_C1000_ab+ABX*I_NAI_D2z_S_C1000_ab;
  Double I_NAI_Dxy_Py_C1000_ab = I_NAI_Fx2y_S_C1000_ab+ABY*I_NAI_Dxy_S_C1000_ab;
  Double I_NAI_Dxz_Py_C1000_ab = I_NAI_Fxyz_S_C1000_ab+ABY*I_NAI_Dxz_S_C1000_ab;
  Double I_NAI_D2y_Py_C1000_ab = I_NAI_F3y_S_C1000_ab+ABY*I_NAI_D2y_S_C1000_ab;
  Double I_NAI_Dyz_Py_C1000_ab = I_NAI_F2yz_S_C1000_ab+ABY*I_NAI_Dyz_S_C1000_ab;
  Double I_NAI_D2z_Py_C1000_ab = I_NAI_Fy2z_S_C1000_ab+ABY*I_NAI_D2z_S_C1000_ab;
  Double I_NAI_Dxz_Pz_C1000_ab = I_NAI_Fx2z_S_C1000_ab+ABZ*I_NAI_Dxz_S_C1000_ab;
  Double I_NAI_Dyz_Pz_C1000_ab = I_NAI_Fy2z_S_C1000_ab+ABZ*I_NAI_Dyz_S_C1000_ab;
  Double I_NAI_D2z_Pz_C1000_ab = I_NAI_F3z_S_C1000_ab+ABZ*I_NAI_D2z_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_D_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_P_C1000_ab
   ************************************************************/
  Double I_NAI_Px_D2x_C1000_ab = I_NAI_D2x_Px_C1000_ab+ABX*I_NAI_Px_Px_C1000_ab;
  Double I_NAI_Py_D2x_C1000_ab = I_NAI_Dxy_Px_C1000_ab+ABX*I_NAI_Py_Px_C1000_ab;
  Double I_NAI_Pz_D2x_C1000_ab = I_NAI_Dxz_Px_C1000_ab+ABX*I_NAI_Pz_Px_C1000_ab;
  Double I_NAI_Px_Dxy_C1000_ab = I_NAI_Dxy_Px_C1000_ab+ABY*I_NAI_Px_Px_C1000_ab;
  Double I_NAI_Py_Dxy_C1000_ab = I_NAI_D2y_Px_C1000_ab+ABY*I_NAI_Py_Px_C1000_ab;
  Double I_NAI_Pz_Dxy_C1000_ab = I_NAI_Dyz_Px_C1000_ab+ABY*I_NAI_Pz_Px_C1000_ab;
  Double I_NAI_Px_Dxz_C1000_ab = I_NAI_Dxz_Px_C1000_ab+ABZ*I_NAI_Px_Px_C1000_ab;
  Double I_NAI_Py_Dxz_C1000_ab = I_NAI_Dyz_Px_C1000_ab+ABZ*I_NAI_Py_Px_C1000_ab;
  Double I_NAI_Pz_Dxz_C1000_ab = I_NAI_D2z_Px_C1000_ab+ABZ*I_NAI_Pz_Px_C1000_ab;
  Double I_NAI_Px_D2y_C1000_ab = I_NAI_Dxy_Py_C1000_ab+ABY*I_NAI_Px_Py_C1000_ab;
  Double I_NAI_Py_D2y_C1000_ab = I_NAI_D2y_Py_C1000_ab+ABY*I_NAI_Py_Py_C1000_ab;
  Double I_NAI_Pz_D2y_C1000_ab = I_NAI_Dyz_Py_C1000_ab+ABY*I_NAI_Pz_Py_C1000_ab;
  Double I_NAI_Px_Dyz_C1000_ab = I_NAI_Dxz_Py_C1000_ab+ABZ*I_NAI_Px_Py_C1000_ab;
  Double I_NAI_Py_Dyz_C1000_ab = I_NAI_Dyz_Py_C1000_ab+ABZ*I_NAI_Py_Py_C1000_ab;
  Double I_NAI_Pz_Dyz_C1000_ab = I_NAI_D2z_Py_C1000_ab+ABZ*I_NAI_Pz_Py_C1000_ab;
  Double I_NAI_Px_D2z_C1000_ab = I_NAI_Dxz_Pz_C1000_ab+ABZ*I_NAI_Px_Pz_C1000_ab;
  Double I_NAI_Py_D2z_C1000_ab = I_NAI_Dyz_Pz_C1000_ab+ABZ*I_NAI_Py_Pz_C1000_ab;
  Double I_NAI_Pz_D2z_C1000_ab = I_NAI_D2z_Pz_C1000_ab+ABZ*I_NAI_Pz_Pz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_ab
   ************************************************************/
  Double I_NAI_D2x_Px_C1001_ab = I_NAI_F3x_S_C1001_ab+ABX*I_NAI_D2x_S_C1001_ab;
  Double I_NAI_Dxy_Px_C1001_ab = I_NAI_F2xy_S_C1001_ab+ABX*I_NAI_Dxy_S_C1001_ab;
  Double I_NAI_Dxz_Px_C1001_ab = I_NAI_F2xz_S_C1001_ab+ABX*I_NAI_Dxz_S_C1001_ab;
  Double I_NAI_D2y_Px_C1001_ab = I_NAI_Fx2y_S_C1001_ab+ABX*I_NAI_D2y_S_C1001_ab;
  Double I_NAI_Dyz_Px_C1001_ab = I_NAI_Fxyz_S_C1001_ab+ABX*I_NAI_Dyz_S_C1001_ab;
  Double I_NAI_D2z_Px_C1001_ab = I_NAI_Fx2z_S_C1001_ab+ABX*I_NAI_D2z_S_C1001_ab;
  Double I_NAI_D2x_Py_C1001_ab = I_NAI_F2xy_S_C1001_ab+ABY*I_NAI_D2x_S_C1001_ab;
  Double I_NAI_Dxy_Py_C1001_ab = I_NAI_Fx2y_S_C1001_ab+ABY*I_NAI_Dxy_S_C1001_ab;
  Double I_NAI_Dxz_Py_C1001_ab = I_NAI_Fxyz_S_C1001_ab+ABY*I_NAI_Dxz_S_C1001_ab;
  Double I_NAI_D2y_Py_C1001_ab = I_NAI_F3y_S_C1001_ab+ABY*I_NAI_D2y_S_C1001_ab;
  Double I_NAI_Dyz_Py_C1001_ab = I_NAI_F2yz_S_C1001_ab+ABY*I_NAI_Dyz_S_C1001_ab;
  Double I_NAI_D2z_Py_C1001_ab = I_NAI_Fy2z_S_C1001_ab+ABY*I_NAI_D2z_S_C1001_ab;
  Double I_NAI_D2x_Pz_C1001_ab = I_NAI_F2xz_S_C1001_ab+ABZ*I_NAI_D2x_S_C1001_ab;
  Double I_NAI_Dxy_Pz_C1001_ab = I_NAI_Fxyz_S_C1001_ab+ABZ*I_NAI_Dxy_S_C1001_ab;
  Double I_NAI_Dxz_Pz_C1001_ab = I_NAI_Fx2z_S_C1001_ab+ABZ*I_NAI_Dxz_S_C1001_ab;
  Double I_NAI_D2y_Pz_C1001_ab = I_NAI_F2yz_S_C1001_ab+ABZ*I_NAI_D2y_S_C1001_ab;
  Double I_NAI_Dyz_Pz_C1001_ab = I_NAI_Fy2z_S_C1001_ab+ABZ*I_NAI_Dyz_S_C1001_ab;
  Double I_NAI_D2z_Pz_C1001_ab = I_NAI_F3z_S_C1001_ab+ABZ*I_NAI_D2z_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1001_ab
   * RHS shell quartet name: SQ_NAI_F_S_C1001_ab
   ************************************************************/
  Double I_NAI_F3x_Px_C1001_ab = I_NAI_G4x_S_C1001_ab+ABX*I_NAI_F3x_S_C1001_ab;
  Double I_NAI_F2xy_Px_C1001_ab = I_NAI_G3xy_S_C1001_ab+ABX*I_NAI_F2xy_S_C1001_ab;
  Double I_NAI_F2xz_Px_C1001_ab = I_NAI_G3xz_S_C1001_ab+ABX*I_NAI_F2xz_S_C1001_ab;
  Double I_NAI_Fx2y_Px_C1001_ab = I_NAI_G2x2y_S_C1001_ab+ABX*I_NAI_Fx2y_S_C1001_ab;
  Double I_NAI_Fxyz_Px_C1001_ab = I_NAI_G2xyz_S_C1001_ab+ABX*I_NAI_Fxyz_S_C1001_ab;
  Double I_NAI_Fx2z_Px_C1001_ab = I_NAI_G2x2z_S_C1001_ab+ABX*I_NAI_Fx2z_S_C1001_ab;
  Double I_NAI_F3y_Px_C1001_ab = I_NAI_Gx3y_S_C1001_ab+ABX*I_NAI_F3y_S_C1001_ab;
  Double I_NAI_F2yz_Px_C1001_ab = I_NAI_Gx2yz_S_C1001_ab+ABX*I_NAI_F2yz_S_C1001_ab;
  Double I_NAI_Fy2z_Px_C1001_ab = I_NAI_Gxy2z_S_C1001_ab+ABX*I_NAI_Fy2z_S_C1001_ab;
  Double I_NAI_F3z_Px_C1001_ab = I_NAI_Gx3z_S_C1001_ab+ABX*I_NAI_F3z_S_C1001_ab;
  Double I_NAI_F2xy_Py_C1001_ab = I_NAI_G2x2y_S_C1001_ab+ABY*I_NAI_F2xy_S_C1001_ab;
  Double I_NAI_F2xz_Py_C1001_ab = I_NAI_G2xyz_S_C1001_ab+ABY*I_NAI_F2xz_S_C1001_ab;
  Double I_NAI_Fx2y_Py_C1001_ab = I_NAI_Gx3y_S_C1001_ab+ABY*I_NAI_Fx2y_S_C1001_ab;
  Double I_NAI_Fxyz_Py_C1001_ab = I_NAI_Gx2yz_S_C1001_ab+ABY*I_NAI_Fxyz_S_C1001_ab;
  Double I_NAI_Fx2z_Py_C1001_ab = I_NAI_Gxy2z_S_C1001_ab+ABY*I_NAI_Fx2z_S_C1001_ab;
  Double I_NAI_F3y_Py_C1001_ab = I_NAI_G4y_S_C1001_ab+ABY*I_NAI_F3y_S_C1001_ab;
  Double I_NAI_F2yz_Py_C1001_ab = I_NAI_G3yz_S_C1001_ab+ABY*I_NAI_F2yz_S_C1001_ab;
  Double I_NAI_Fy2z_Py_C1001_ab = I_NAI_G2y2z_S_C1001_ab+ABY*I_NAI_Fy2z_S_C1001_ab;
  Double I_NAI_F3z_Py_C1001_ab = I_NAI_Gy3z_S_C1001_ab+ABY*I_NAI_F3z_S_C1001_ab;
  Double I_NAI_F2xz_Pz_C1001_ab = I_NAI_G2x2z_S_C1001_ab+ABZ*I_NAI_F2xz_S_C1001_ab;
  Double I_NAI_Fxyz_Pz_C1001_ab = I_NAI_Gxy2z_S_C1001_ab+ABZ*I_NAI_Fxyz_S_C1001_ab;
  Double I_NAI_Fx2z_Pz_C1001_ab = I_NAI_Gx3z_S_C1001_ab+ABZ*I_NAI_Fx2z_S_C1001_ab;
  Double I_NAI_F2yz_Pz_C1001_ab = I_NAI_G2y2z_S_C1001_ab+ABZ*I_NAI_F2yz_S_C1001_ab;
  Double I_NAI_Fy2z_Pz_C1001_ab = I_NAI_Gy3z_S_C1001_ab+ABZ*I_NAI_Fy2z_S_C1001_ab;
  Double I_NAI_F3z_Pz_C1001_ab = I_NAI_G4z_S_C1001_ab+ABZ*I_NAI_F3z_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_P_C1001_ab
   ************************************************************/
  Double I_NAI_D2x_D2x_C1001_ab = I_NAI_F3x_Px_C1001_ab+ABX*I_NAI_D2x_Px_C1001_ab;
  Double I_NAI_Dxy_D2x_C1001_ab = I_NAI_F2xy_Px_C1001_ab+ABX*I_NAI_Dxy_Px_C1001_ab;
  Double I_NAI_Dxz_D2x_C1001_ab = I_NAI_F2xz_Px_C1001_ab+ABX*I_NAI_Dxz_Px_C1001_ab;
  Double I_NAI_D2y_D2x_C1001_ab = I_NAI_Fx2y_Px_C1001_ab+ABX*I_NAI_D2y_Px_C1001_ab;
  Double I_NAI_Dyz_D2x_C1001_ab = I_NAI_Fxyz_Px_C1001_ab+ABX*I_NAI_Dyz_Px_C1001_ab;
  Double I_NAI_D2z_D2x_C1001_ab = I_NAI_Fx2z_Px_C1001_ab+ABX*I_NAI_D2z_Px_C1001_ab;
  Double I_NAI_D2x_Dxy_C1001_ab = I_NAI_F2xy_Px_C1001_ab+ABY*I_NAI_D2x_Px_C1001_ab;
  Double I_NAI_Dxy_Dxy_C1001_ab = I_NAI_Fx2y_Px_C1001_ab+ABY*I_NAI_Dxy_Px_C1001_ab;
  Double I_NAI_Dxz_Dxy_C1001_ab = I_NAI_Fxyz_Px_C1001_ab+ABY*I_NAI_Dxz_Px_C1001_ab;
  Double I_NAI_D2y_Dxy_C1001_ab = I_NAI_F3y_Px_C1001_ab+ABY*I_NAI_D2y_Px_C1001_ab;
  Double I_NAI_Dyz_Dxy_C1001_ab = I_NAI_F2yz_Px_C1001_ab+ABY*I_NAI_Dyz_Px_C1001_ab;
  Double I_NAI_D2z_Dxy_C1001_ab = I_NAI_Fy2z_Px_C1001_ab+ABY*I_NAI_D2z_Px_C1001_ab;
  Double I_NAI_D2x_Dxz_C1001_ab = I_NAI_F2xz_Px_C1001_ab+ABZ*I_NAI_D2x_Px_C1001_ab;
  Double I_NAI_Dxy_Dxz_C1001_ab = I_NAI_Fxyz_Px_C1001_ab+ABZ*I_NAI_Dxy_Px_C1001_ab;
  Double I_NAI_Dxz_Dxz_C1001_ab = I_NAI_Fx2z_Px_C1001_ab+ABZ*I_NAI_Dxz_Px_C1001_ab;
  Double I_NAI_D2y_Dxz_C1001_ab = I_NAI_F2yz_Px_C1001_ab+ABZ*I_NAI_D2y_Px_C1001_ab;
  Double I_NAI_Dyz_Dxz_C1001_ab = I_NAI_Fy2z_Px_C1001_ab+ABZ*I_NAI_Dyz_Px_C1001_ab;
  Double I_NAI_D2z_Dxz_C1001_ab = I_NAI_F3z_Px_C1001_ab+ABZ*I_NAI_D2z_Px_C1001_ab;
  Double I_NAI_D2x_D2y_C1001_ab = I_NAI_F2xy_Py_C1001_ab+ABY*I_NAI_D2x_Py_C1001_ab;
  Double I_NAI_Dxy_D2y_C1001_ab = I_NAI_Fx2y_Py_C1001_ab+ABY*I_NAI_Dxy_Py_C1001_ab;
  Double I_NAI_Dxz_D2y_C1001_ab = I_NAI_Fxyz_Py_C1001_ab+ABY*I_NAI_Dxz_Py_C1001_ab;
  Double I_NAI_D2y_D2y_C1001_ab = I_NAI_F3y_Py_C1001_ab+ABY*I_NAI_D2y_Py_C1001_ab;
  Double I_NAI_Dyz_D2y_C1001_ab = I_NAI_F2yz_Py_C1001_ab+ABY*I_NAI_Dyz_Py_C1001_ab;
  Double I_NAI_D2z_D2y_C1001_ab = I_NAI_Fy2z_Py_C1001_ab+ABY*I_NAI_D2z_Py_C1001_ab;
  Double I_NAI_D2x_Dyz_C1001_ab = I_NAI_F2xz_Py_C1001_ab+ABZ*I_NAI_D2x_Py_C1001_ab;
  Double I_NAI_Dxy_Dyz_C1001_ab = I_NAI_Fxyz_Py_C1001_ab+ABZ*I_NAI_Dxy_Py_C1001_ab;
  Double I_NAI_Dxz_Dyz_C1001_ab = I_NAI_Fx2z_Py_C1001_ab+ABZ*I_NAI_Dxz_Py_C1001_ab;
  Double I_NAI_D2y_Dyz_C1001_ab = I_NAI_F2yz_Py_C1001_ab+ABZ*I_NAI_D2y_Py_C1001_ab;
  Double I_NAI_Dyz_Dyz_C1001_ab = I_NAI_Fy2z_Py_C1001_ab+ABZ*I_NAI_Dyz_Py_C1001_ab;
  Double I_NAI_D2z_Dyz_C1001_ab = I_NAI_F3z_Py_C1001_ab+ABZ*I_NAI_D2z_Py_C1001_ab;
  Double I_NAI_D2x_D2z_C1001_ab = I_NAI_F2xz_Pz_C1001_ab+ABZ*I_NAI_D2x_Pz_C1001_ab;
  Double I_NAI_Dxy_D2z_C1001_ab = I_NAI_Fxyz_Pz_C1001_ab+ABZ*I_NAI_Dxy_Pz_C1001_ab;
  Double I_NAI_Dxz_D2z_C1001_ab = I_NAI_Fx2z_Pz_C1001_ab+ABZ*I_NAI_Dxz_Pz_C1001_ab;
  Double I_NAI_D2y_D2z_C1001_ab = I_NAI_F2yz_Pz_C1001_ab+ABZ*I_NAI_D2y_Pz_C1001_ab;
  Double I_NAI_Dyz_D2z_C1001_ab = I_NAI_Fy2z_Pz_C1001_ab+ABZ*I_NAI_Dyz_Pz_C1001_ab;
  Double I_NAI_D2z_D2z_C1001_ab = I_NAI_F3z_Pz_C1001_ab+ABZ*I_NAI_D2z_Pz_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_bb
   ************************************************************/
  Double I_NAI_Px_Px_C1_bb = I_NAI_D2x_S_C1_bb+ABX*I_NAI_Px_S_C1_bb;
  Double I_NAI_Py_Px_C1_bb = I_NAI_Dxy_S_C1_bb+ABX*I_NAI_Py_S_C1_bb;
  Double I_NAI_Pz_Px_C1_bb = I_NAI_Dxz_S_C1_bb+ABX*I_NAI_Pz_S_C1_bb;
  Double I_NAI_Px_Py_C1_bb = I_NAI_Dxy_S_C1_bb+ABY*I_NAI_Px_S_C1_bb;
  Double I_NAI_Py_Py_C1_bb = I_NAI_D2y_S_C1_bb+ABY*I_NAI_Py_S_C1_bb;
  Double I_NAI_Pz_Py_C1_bb = I_NAI_Dyz_S_C1_bb+ABY*I_NAI_Pz_S_C1_bb;
  Double I_NAI_Px_Pz_C1_bb = I_NAI_Dxz_S_C1_bb+ABZ*I_NAI_Px_S_C1_bb;
  Double I_NAI_Py_Pz_C1_bb = I_NAI_Dyz_S_C1_bb+ABZ*I_NAI_Py_S_C1_bb;
  Double I_NAI_Pz_Pz_C1_bb = I_NAI_D2z_S_C1_bb+ABZ*I_NAI_Pz_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 4 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_bb
   * RHS shell quartet name: SQ_NAI_D_S_C1_bb
   ************************************************************/
  Double I_NAI_D2x_Px_C1_bb = I_NAI_F3x_S_C1_bb+ABX*I_NAI_D2x_S_C1_bb;
  Double I_NAI_Dxy_Px_C1_bb = I_NAI_F2xy_S_C1_bb+ABX*I_NAI_Dxy_S_C1_bb;
  Double I_NAI_Dxz_Px_C1_bb = I_NAI_F2xz_S_C1_bb+ABX*I_NAI_Dxz_S_C1_bb;
  Double I_NAI_D2y_Px_C1_bb = I_NAI_Fx2y_S_C1_bb+ABX*I_NAI_D2y_S_C1_bb;
  Double I_NAI_Dyz_Px_C1_bb = I_NAI_Fxyz_S_C1_bb+ABX*I_NAI_Dyz_S_C1_bb;
  Double I_NAI_D2z_Px_C1_bb = I_NAI_Fx2z_S_C1_bb+ABX*I_NAI_D2z_S_C1_bb;
  Double I_NAI_Dxy_Py_C1_bb = I_NAI_Fx2y_S_C1_bb+ABY*I_NAI_Dxy_S_C1_bb;
  Double I_NAI_Dxz_Py_C1_bb = I_NAI_Fxyz_S_C1_bb+ABY*I_NAI_Dxz_S_C1_bb;
  Double I_NAI_D2y_Py_C1_bb = I_NAI_F3y_S_C1_bb+ABY*I_NAI_D2y_S_C1_bb;
  Double I_NAI_Dyz_Py_C1_bb = I_NAI_F2yz_S_C1_bb+ABY*I_NAI_Dyz_S_C1_bb;
  Double I_NAI_D2z_Py_C1_bb = I_NAI_Fy2z_S_C1_bb+ABY*I_NAI_D2z_S_C1_bb;
  Double I_NAI_Dxz_Pz_C1_bb = I_NAI_Fx2z_S_C1_bb+ABZ*I_NAI_Dxz_S_C1_bb;
  Double I_NAI_Dyz_Pz_C1_bb = I_NAI_Fy2z_S_C1_bb+ABZ*I_NAI_Dyz_S_C1_bb;
  Double I_NAI_D2z_Pz_C1_bb = I_NAI_F3z_S_C1_bb+ABZ*I_NAI_D2z_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_D_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1_bb
   ************************************************************/
  Double I_NAI_Px_D2x_C1_bb = I_NAI_D2x_Px_C1_bb+ABX*I_NAI_Px_Px_C1_bb;
  Double I_NAI_Py_D2x_C1_bb = I_NAI_Dxy_Px_C1_bb+ABX*I_NAI_Py_Px_C1_bb;
  Double I_NAI_Pz_D2x_C1_bb = I_NAI_Dxz_Px_C1_bb+ABX*I_NAI_Pz_Px_C1_bb;
  Double I_NAI_Px_Dxy_C1_bb = I_NAI_Dxy_Px_C1_bb+ABY*I_NAI_Px_Px_C1_bb;
  Double I_NAI_Py_Dxy_C1_bb = I_NAI_D2y_Px_C1_bb+ABY*I_NAI_Py_Px_C1_bb;
  Double I_NAI_Pz_Dxy_C1_bb = I_NAI_Dyz_Px_C1_bb+ABY*I_NAI_Pz_Px_C1_bb;
  Double I_NAI_Px_Dxz_C1_bb = I_NAI_Dxz_Px_C1_bb+ABZ*I_NAI_Px_Px_C1_bb;
  Double I_NAI_Py_Dxz_C1_bb = I_NAI_Dyz_Px_C1_bb+ABZ*I_NAI_Py_Px_C1_bb;
  Double I_NAI_Pz_Dxz_C1_bb = I_NAI_D2z_Px_C1_bb+ABZ*I_NAI_Pz_Px_C1_bb;
  Double I_NAI_Px_D2y_C1_bb = I_NAI_Dxy_Py_C1_bb+ABY*I_NAI_Px_Py_C1_bb;
  Double I_NAI_Py_D2y_C1_bb = I_NAI_D2y_Py_C1_bb+ABY*I_NAI_Py_Py_C1_bb;
  Double I_NAI_Pz_D2y_C1_bb = I_NAI_Dyz_Py_C1_bb+ABY*I_NAI_Pz_Py_C1_bb;
  Double I_NAI_Px_Dyz_C1_bb = I_NAI_Dxz_Py_C1_bb+ABZ*I_NAI_Px_Py_C1_bb;
  Double I_NAI_Py_Dyz_C1_bb = I_NAI_Dyz_Py_C1_bb+ABZ*I_NAI_Py_Py_C1_bb;
  Double I_NAI_Pz_Dyz_C1_bb = I_NAI_D2z_Py_C1_bb+ABZ*I_NAI_Pz_Py_C1_bb;
  Double I_NAI_Px_D2z_C1_bb = I_NAI_Dxz_Pz_C1_bb+ABZ*I_NAI_Px_Pz_C1_bb;
  Double I_NAI_Py_D2z_C1_bb = I_NAI_Dyz_Pz_C1_bb+ABZ*I_NAI_Py_Pz_C1_bb;
  Double I_NAI_Pz_D2z_C1_bb = I_NAI_D2z_Pz_C1_bb+ABZ*I_NAI_Pz_Pz_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1001_bb
   ************************************************************/
  Double I_NAI_Px_Px_C1001_bb = I_NAI_D2x_S_C1001_bb+ABX*I_NAI_Px_S_C1001_bb;
  Double I_NAI_Py_Px_C1001_bb = I_NAI_Dxy_S_C1001_bb+ABX*I_NAI_Py_S_C1001_bb;
  Double I_NAI_Pz_Px_C1001_bb = I_NAI_Dxz_S_C1001_bb+ABX*I_NAI_Pz_S_C1001_bb;
  Double I_NAI_Px_Py_C1001_bb = I_NAI_Dxy_S_C1001_bb+ABY*I_NAI_Px_S_C1001_bb;
  Double I_NAI_Py_Py_C1001_bb = I_NAI_D2y_S_C1001_bb+ABY*I_NAI_Py_S_C1001_bb;
  Double I_NAI_Pz_Py_C1001_bb = I_NAI_Dyz_S_C1001_bb+ABY*I_NAI_Pz_S_C1001_bb;
  Double I_NAI_Px_Pz_C1001_bb = I_NAI_Dxz_S_C1001_bb+ABZ*I_NAI_Px_S_C1001_bb;
  Double I_NAI_Py_Pz_C1001_bb = I_NAI_Dyz_S_C1001_bb+ABZ*I_NAI_Py_S_C1001_bb;
  Double I_NAI_Pz_Pz_C1001_bb = I_NAI_D2z_S_C1001_bb+ABZ*I_NAI_Pz_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1001_bb
   * RHS shell quartet name: SQ_NAI_D_S_C1001_bb
   ************************************************************/
  Double I_NAI_D2x_Px_C1001_bb = I_NAI_F3x_S_C1001_bb+ABX*I_NAI_D2x_S_C1001_bb;
  Double I_NAI_Dxy_Px_C1001_bb = I_NAI_F2xy_S_C1001_bb+ABX*I_NAI_Dxy_S_C1001_bb;
  Double I_NAI_Dxz_Px_C1001_bb = I_NAI_F2xz_S_C1001_bb+ABX*I_NAI_Dxz_S_C1001_bb;
  Double I_NAI_D2y_Px_C1001_bb = I_NAI_Fx2y_S_C1001_bb+ABX*I_NAI_D2y_S_C1001_bb;
  Double I_NAI_Dyz_Px_C1001_bb = I_NAI_Fxyz_S_C1001_bb+ABX*I_NAI_Dyz_S_C1001_bb;
  Double I_NAI_D2z_Px_C1001_bb = I_NAI_Fx2z_S_C1001_bb+ABX*I_NAI_D2z_S_C1001_bb;
  Double I_NAI_D2x_Py_C1001_bb = I_NAI_F2xy_S_C1001_bb+ABY*I_NAI_D2x_S_C1001_bb;
  Double I_NAI_Dxy_Py_C1001_bb = I_NAI_Fx2y_S_C1001_bb+ABY*I_NAI_Dxy_S_C1001_bb;
  Double I_NAI_Dxz_Py_C1001_bb = I_NAI_Fxyz_S_C1001_bb+ABY*I_NAI_Dxz_S_C1001_bb;
  Double I_NAI_D2y_Py_C1001_bb = I_NAI_F3y_S_C1001_bb+ABY*I_NAI_D2y_S_C1001_bb;
  Double I_NAI_Dyz_Py_C1001_bb = I_NAI_F2yz_S_C1001_bb+ABY*I_NAI_Dyz_S_C1001_bb;
  Double I_NAI_D2z_Py_C1001_bb = I_NAI_Fy2z_S_C1001_bb+ABY*I_NAI_D2z_S_C1001_bb;
  Double I_NAI_D2x_Pz_C1001_bb = I_NAI_F2xz_S_C1001_bb+ABZ*I_NAI_D2x_S_C1001_bb;
  Double I_NAI_Dxy_Pz_C1001_bb = I_NAI_Fxyz_S_C1001_bb+ABZ*I_NAI_Dxy_S_C1001_bb;
  Double I_NAI_Dxz_Pz_C1001_bb = I_NAI_Fx2z_S_C1001_bb+ABZ*I_NAI_Dxz_S_C1001_bb;
  Double I_NAI_D2y_Pz_C1001_bb = I_NAI_F2yz_S_C1001_bb+ABZ*I_NAI_D2y_S_C1001_bb;
  Double I_NAI_Dyz_Pz_C1001_bb = I_NAI_Fy2z_S_C1001_bb+ABZ*I_NAI_Dyz_S_C1001_bb;
  Double I_NAI_D2z_Pz_C1001_bb = I_NAI_F3z_S_C1001_bb+ABZ*I_NAI_D2z_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_D_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_bb
   ************************************************************/
  Double I_NAI_Px_D2x_C1001_bb = I_NAI_D2x_Px_C1001_bb+ABX*I_NAI_Px_Px_C1001_bb;
  Double I_NAI_Py_D2x_C1001_bb = I_NAI_Dxy_Px_C1001_bb+ABX*I_NAI_Py_Px_C1001_bb;
  Double I_NAI_Pz_D2x_C1001_bb = I_NAI_Dxz_Px_C1001_bb+ABX*I_NAI_Pz_Px_C1001_bb;
  Double I_NAI_Px_Dxy_C1001_bb = I_NAI_Dxy_Px_C1001_bb+ABY*I_NAI_Px_Px_C1001_bb;
  Double I_NAI_Py_Dxy_C1001_bb = I_NAI_D2y_Px_C1001_bb+ABY*I_NAI_Py_Px_C1001_bb;
  Double I_NAI_Pz_Dxy_C1001_bb = I_NAI_Dyz_Px_C1001_bb+ABY*I_NAI_Pz_Px_C1001_bb;
  Double I_NAI_Px_D2y_C1001_bb = I_NAI_Dxy_Py_C1001_bb+ABY*I_NAI_Px_Py_C1001_bb;
  Double I_NAI_Py_D2y_C1001_bb = I_NAI_D2y_Py_C1001_bb+ABY*I_NAI_Py_Py_C1001_bb;
  Double I_NAI_Pz_D2y_C1001_bb = I_NAI_Dyz_Py_C1001_bb+ABY*I_NAI_Pz_Py_C1001_bb;
  Double I_NAI_Px_D2z_C1001_bb = I_NAI_Dxz_Pz_C1001_bb+ABZ*I_NAI_Px_Pz_C1001_bb;
  Double I_NAI_Py_D2z_C1001_bb = I_NAI_Dyz_Pz_C1001_bb+ABZ*I_NAI_Py_Pz_C1001_bb;
  Double I_NAI_Pz_D2z_C1001_bb = I_NAI_D2z_Pz_C1001_bb+ABZ*I_NAI_Pz_Pz_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 10 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1001_bb
   * RHS shell quartet name: SQ_NAI_F_S_C1001_bb
   ************************************************************/
  Double I_NAI_F3x_Px_C1001_bb = I_NAI_G4x_S_C1001_bb+ABX*I_NAI_F3x_S_C1001_bb;
  Double I_NAI_F2xy_Px_C1001_bb = I_NAI_G3xy_S_C1001_bb+ABX*I_NAI_F2xy_S_C1001_bb;
  Double I_NAI_F2xz_Px_C1001_bb = I_NAI_G3xz_S_C1001_bb+ABX*I_NAI_F2xz_S_C1001_bb;
  Double I_NAI_Fx2y_Px_C1001_bb = I_NAI_G2x2y_S_C1001_bb+ABX*I_NAI_Fx2y_S_C1001_bb;
  Double I_NAI_Fxyz_Px_C1001_bb = I_NAI_G2xyz_S_C1001_bb+ABX*I_NAI_Fxyz_S_C1001_bb;
  Double I_NAI_Fx2z_Px_C1001_bb = I_NAI_G2x2z_S_C1001_bb+ABX*I_NAI_Fx2z_S_C1001_bb;
  Double I_NAI_F2yz_Px_C1001_bb = I_NAI_Gx2yz_S_C1001_bb+ABX*I_NAI_F2yz_S_C1001_bb;
  Double I_NAI_Fy2z_Px_C1001_bb = I_NAI_Gxy2z_S_C1001_bb+ABX*I_NAI_Fy2z_S_C1001_bb;
  Double I_NAI_F2xy_Py_C1001_bb = I_NAI_G2x2y_S_C1001_bb+ABY*I_NAI_F2xy_S_C1001_bb;
  Double I_NAI_Fx2y_Py_C1001_bb = I_NAI_Gx3y_S_C1001_bb+ABY*I_NAI_Fx2y_S_C1001_bb;
  Double I_NAI_Fxyz_Py_C1001_bb = I_NAI_Gx2yz_S_C1001_bb+ABY*I_NAI_Fxyz_S_C1001_bb;
  Double I_NAI_F3y_Py_C1001_bb = I_NAI_G4y_S_C1001_bb+ABY*I_NAI_F3y_S_C1001_bb;
  Double I_NAI_F2yz_Py_C1001_bb = I_NAI_G3yz_S_C1001_bb+ABY*I_NAI_F2yz_S_C1001_bb;
  Double I_NAI_Fy2z_Py_C1001_bb = I_NAI_G2y2z_S_C1001_bb+ABY*I_NAI_Fy2z_S_C1001_bb;
  Double I_NAI_F2xz_Pz_C1001_bb = I_NAI_G2x2z_S_C1001_bb+ABZ*I_NAI_F2xz_S_C1001_bb;
  Double I_NAI_Fxyz_Pz_C1001_bb = I_NAI_Gxy2z_S_C1001_bb+ABZ*I_NAI_Fxyz_S_C1001_bb;
  Double I_NAI_Fx2z_Pz_C1001_bb = I_NAI_Gx3z_S_C1001_bb+ABZ*I_NAI_Fx2z_S_C1001_bb;
  Double I_NAI_F2yz_Pz_C1001_bb = I_NAI_G2y2z_S_C1001_bb+ABZ*I_NAI_F2yz_S_C1001_bb;
  Double I_NAI_Fy2z_Pz_C1001_bb = I_NAI_Gy3z_S_C1001_bb+ABZ*I_NAI_Fy2z_S_C1001_bb;
  Double I_NAI_F3z_Pz_C1001_bb = I_NAI_G4z_S_C1001_bb+ABZ*I_NAI_F3z_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_bb
   * RHS shell quartet name: SQ_NAI_D_P_C1001_bb
   ************************************************************/
  Double I_NAI_D2x_D2x_C1001_bb = I_NAI_F3x_Px_C1001_bb+ABX*I_NAI_D2x_Px_C1001_bb;
  Double I_NAI_Dxy_D2x_C1001_bb = I_NAI_F2xy_Px_C1001_bb+ABX*I_NAI_Dxy_Px_C1001_bb;
  Double I_NAI_Dxz_D2x_C1001_bb = I_NAI_F2xz_Px_C1001_bb+ABX*I_NAI_Dxz_Px_C1001_bb;
  Double I_NAI_D2y_D2x_C1001_bb = I_NAI_Fx2y_Px_C1001_bb+ABX*I_NAI_D2y_Px_C1001_bb;
  Double I_NAI_Dyz_D2x_C1001_bb = I_NAI_Fxyz_Px_C1001_bb+ABX*I_NAI_Dyz_Px_C1001_bb;
  Double I_NAI_D2z_D2x_C1001_bb = I_NAI_Fx2z_Px_C1001_bb+ABX*I_NAI_D2z_Px_C1001_bb;
  Double I_NAI_Dxz_Dxy_C1001_bb = I_NAI_Fxyz_Px_C1001_bb+ABY*I_NAI_Dxz_Px_C1001_bb;
  Double I_NAI_Dyz_Dxy_C1001_bb = I_NAI_F2yz_Px_C1001_bb+ABY*I_NAI_Dyz_Px_C1001_bb;
  Double I_NAI_D2z_Dxy_C1001_bb = I_NAI_Fy2z_Px_C1001_bb+ABY*I_NAI_D2z_Px_C1001_bb;
  Double I_NAI_D2x_D2y_C1001_bb = I_NAI_F2xy_Py_C1001_bb+ABY*I_NAI_D2x_Py_C1001_bb;
  Double I_NAI_Dxy_D2y_C1001_bb = I_NAI_Fx2y_Py_C1001_bb+ABY*I_NAI_Dxy_Py_C1001_bb;
  Double I_NAI_Dxz_D2y_C1001_bb = I_NAI_Fxyz_Py_C1001_bb+ABY*I_NAI_Dxz_Py_C1001_bb;
  Double I_NAI_D2y_D2y_C1001_bb = I_NAI_F3y_Py_C1001_bb+ABY*I_NAI_D2y_Py_C1001_bb;
  Double I_NAI_Dyz_D2y_C1001_bb = I_NAI_F2yz_Py_C1001_bb+ABY*I_NAI_Dyz_Py_C1001_bb;
  Double I_NAI_D2z_D2y_C1001_bb = I_NAI_Fy2z_Py_C1001_bb+ABY*I_NAI_D2z_Py_C1001_bb;
  Double I_NAI_D2x_D2z_C1001_bb = I_NAI_F2xz_Pz_C1001_bb+ABZ*I_NAI_D2x_Pz_C1001_bb;
  Double I_NAI_Dxy_D2z_C1001_bb = I_NAI_Fxyz_Pz_C1001_bb+ABZ*I_NAI_Dxy_Pz_C1001_bb;
  Double I_NAI_Dxz_D2z_C1001_bb = I_NAI_Fx2z_Pz_C1001_bb+ABZ*I_NAI_Dxz_Pz_C1001_bb;
  Double I_NAI_D2y_D2z_C1001_bb = I_NAI_F2yz_Pz_C1001_bb+ABZ*I_NAI_D2y_Pz_C1001_bb;
  Double I_NAI_Dyz_D2z_C1001_bb = I_NAI_Fy2z_Pz_C1001_bb+ABZ*I_NAI_Dyz_Pz_C1001_bb;
  Double I_NAI_D2z_D2z_C1001_bb = I_NAI_F3z_Pz_C1001_bb+ABZ*I_NAI_D2z_Pz_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_F_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_D_C1001_bb
   ************************************************************/
  Double I_NAI_Px_F3x_C1001_bb = I_NAI_D2x_D2x_C1001_bb+ABX*I_NAI_Px_D2x_C1001_bb;
  Double I_NAI_Py_F3x_C1001_bb = I_NAI_Dxy_D2x_C1001_bb+ABX*I_NAI_Py_D2x_C1001_bb;
  Double I_NAI_Pz_F3x_C1001_bb = I_NAI_Dxz_D2x_C1001_bb+ABX*I_NAI_Pz_D2x_C1001_bb;
  Double I_NAI_Px_F2xy_C1001_bb = I_NAI_Dxy_D2x_C1001_bb+ABY*I_NAI_Px_D2x_C1001_bb;
  Double I_NAI_Py_F2xy_C1001_bb = I_NAI_D2y_D2x_C1001_bb+ABY*I_NAI_Py_D2x_C1001_bb;
  Double I_NAI_Pz_F2xy_C1001_bb = I_NAI_Dyz_D2x_C1001_bb+ABY*I_NAI_Pz_D2x_C1001_bb;
  Double I_NAI_Px_F2xz_C1001_bb = I_NAI_Dxz_D2x_C1001_bb+ABZ*I_NAI_Px_D2x_C1001_bb;
  Double I_NAI_Py_F2xz_C1001_bb = I_NAI_Dyz_D2x_C1001_bb+ABZ*I_NAI_Py_D2x_C1001_bb;
  Double I_NAI_Pz_F2xz_C1001_bb = I_NAI_D2z_D2x_C1001_bb+ABZ*I_NAI_Pz_D2x_C1001_bb;
  Double I_NAI_Px_Fx2y_C1001_bb = I_NAI_D2x_D2y_C1001_bb+ABX*I_NAI_Px_D2y_C1001_bb;
  Double I_NAI_Py_Fx2y_C1001_bb = I_NAI_Dxy_D2y_C1001_bb+ABX*I_NAI_Py_D2y_C1001_bb;
  Double I_NAI_Pz_Fx2y_C1001_bb = I_NAI_Dxz_D2y_C1001_bb+ABX*I_NAI_Pz_D2y_C1001_bb;
  Double I_NAI_Px_Fxyz_C1001_bb = I_NAI_Dxz_Dxy_C1001_bb+ABZ*I_NAI_Px_Dxy_C1001_bb;
  Double I_NAI_Py_Fxyz_C1001_bb = I_NAI_Dyz_Dxy_C1001_bb+ABZ*I_NAI_Py_Dxy_C1001_bb;
  Double I_NAI_Pz_Fxyz_C1001_bb = I_NAI_D2z_Dxy_C1001_bb+ABZ*I_NAI_Pz_Dxy_C1001_bb;
  Double I_NAI_Px_Fx2z_C1001_bb = I_NAI_D2x_D2z_C1001_bb+ABX*I_NAI_Px_D2z_C1001_bb;
  Double I_NAI_Py_Fx2z_C1001_bb = I_NAI_Dxy_D2z_C1001_bb+ABX*I_NAI_Py_D2z_C1001_bb;
  Double I_NAI_Pz_Fx2z_C1001_bb = I_NAI_Dxz_D2z_C1001_bb+ABX*I_NAI_Pz_D2z_C1001_bb;
  Double I_NAI_Px_F3y_C1001_bb = I_NAI_Dxy_D2y_C1001_bb+ABY*I_NAI_Px_D2y_C1001_bb;
  Double I_NAI_Py_F3y_C1001_bb = I_NAI_D2y_D2y_C1001_bb+ABY*I_NAI_Py_D2y_C1001_bb;
  Double I_NAI_Pz_F3y_C1001_bb = I_NAI_Dyz_D2y_C1001_bb+ABY*I_NAI_Pz_D2y_C1001_bb;
  Double I_NAI_Px_F2yz_C1001_bb = I_NAI_Dxz_D2y_C1001_bb+ABZ*I_NAI_Px_D2y_C1001_bb;
  Double I_NAI_Py_F2yz_C1001_bb = I_NAI_Dyz_D2y_C1001_bb+ABZ*I_NAI_Py_D2y_C1001_bb;
  Double I_NAI_Pz_F2yz_C1001_bb = I_NAI_D2z_D2y_C1001_bb+ABZ*I_NAI_Pz_D2y_C1001_bb;
  Double I_NAI_Px_Fy2z_C1001_bb = I_NAI_Dxy_D2z_C1001_bb+ABY*I_NAI_Px_D2z_C1001_bb;
  Double I_NAI_Py_Fy2z_C1001_bb = I_NAI_D2y_D2z_C1001_bb+ABY*I_NAI_Py_D2z_C1001_bb;
  Double I_NAI_Pz_Fy2z_C1001_bb = I_NAI_Dyz_D2z_C1001_bb+ABY*I_NAI_Pz_D2z_C1001_bb;
  Double I_NAI_Px_F3z_C1001_bb = I_NAI_Dxz_D2z_C1001_bb+ABZ*I_NAI_Px_D2z_C1001_bb;
  Double I_NAI_Py_F3z_C1001_bb = I_NAI_Dyz_D2z_C1001_bb+ABZ*I_NAI_Py_D2z_C1001_bb;
  Double I_NAI_Pz_F3z_C1001_bb = I_NAI_D2z_D2z_C1001_bb+ABZ*I_NAI_Pz_D2z_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[0] = 4.0E0*I_NAI_D2x_S_C0_aa-2.0E0*1*I_NAI_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[1] = 4.0E0*I_NAI_F3x_S_C1_aa-2.0E0*1*I_NAI_Px_S_C1_a-2.0E0*2*I_NAI_Px_S_C1_a;
  abcd[2] = 4.0E0*I_NAI_F2xy_S_C1_aa-2.0E0*1*I_NAI_Py_S_C1_a;
  abcd[3] = 4.0E0*I_NAI_F2xz_S_C1_aa-2.0E0*1*I_NAI_Pz_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[4] = 4.0E0*I_NAI_D2x_Px_C1000_aa-2.0E0*1*I_NAI_S_Px_C1000_a;
  abcd[8] = 4.0E0*I_NAI_D2x_Py_C1000_aa-2.0E0*1*I_NAI_S_Py_C1000_a;
  abcd[12] = 4.0E0*I_NAI_D2x_Pz_C1000_aa-2.0E0*1*I_NAI_S_Pz_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[5] = 4.0E0*I_NAI_F3x_Px_C1001_aa-2.0E0*1*I_NAI_Px_Px_C1001_a-2.0E0*2*I_NAI_Px_Px_C1001_a;
  abcd[6] = 4.0E0*I_NAI_F2xy_Px_C1001_aa-2.0E0*1*I_NAI_Py_Px_C1001_a;
  abcd[7] = 4.0E0*I_NAI_F2xz_Px_C1001_aa-2.0E0*1*I_NAI_Pz_Px_C1001_a;
  abcd[9] = 4.0E0*I_NAI_F3x_Py_C1001_aa-2.0E0*1*I_NAI_Px_Py_C1001_a-2.0E0*2*I_NAI_Px_Py_C1001_a;
  abcd[10] = 4.0E0*I_NAI_F2xy_Py_C1001_aa-2.0E0*1*I_NAI_Py_Py_C1001_a;
  abcd[11] = 4.0E0*I_NAI_F2xz_Py_C1001_aa-2.0E0*1*I_NAI_Pz_Py_C1001_a;
  abcd[13] = 4.0E0*I_NAI_F3x_Pz_C1001_aa-2.0E0*1*I_NAI_Px_Pz_C1001_a-2.0E0*2*I_NAI_Px_Pz_C1001_a;
  abcd[14] = 4.0E0*I_NAI_F2xy_Pz_C1001_aa-2.0E0*1*I_NAI_Py_Pz_C1001_a;
  abcd[15] = 4.0E0*I_NAI_F2xz_Pz_C1001_aa-2.0E0*1*I_NAI_Pz_Pz_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[16] = 4.0E0*I_NAI_Dxy_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[17] = 4.0E0*I_NAI_F2xy_S_C1_aa-2.0E0*1*I_NAI_Py_S_C1_a;
  abcd[18] = 4.0E0*I_NAI_Fx2y_S_C1_aa-2.0E0*1*I_NAI_Px_S_C1_a;
  abcd[19] = 4.0E0*I_NAI_Fxyz_S_C1_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[20] = 4.0E0*I_NAI_Dxy_Px_C1000_aa;
  abcd[24] = 4.0E0*I_NAI_Dxy_Py_C1000_aa;
  abcd[28] = 4.0E0*I_NAI_Dxy_Pz_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[21] = 4.0E0*I_NAI_F2xy_Px_C1001_aa-2.0E0*1*I_NAI_Py_Px_C1001_a;
  abcd[22] = 4.0E0*I_NAI_Fx2y_Px_C1001_aa-2.0E0*1*I_NAI_Px_Px_C1001_a;
  abcd[23] = 4.0E0*I_NAI_Fxyz_Px_C1001_aa;
  abcd[25] = 4.0E0*I_NAI_F2xy_Py_C1001_aa-2.0E0*1*I_NAI_Py_Py_C1001_a;
  abcd[26] = 4.0E0*I_NAI_Fx2y_Py_C1001_aa-2.0E0*1*I_NAI_Px_Py_C1001_a;
  abcd[27] = 4.0E0*I_NAI_Fxyz_Py_C1001_aa;
  abcd[29] = 4.0E0*I_NAI_F2xy_Pz_C1001_aa-2.0E0*1*I_NAI_Py_Pz_C1001_a;
  abcd[30] = 4.0E0*I_NAI_Fx2y_Pz_C1001_aa-2.0E0*1*I_NAI_Px_Pz_C1001_a;
  abcd[31] = 4.0E0*I_NAI_Fxyz_Pz_C1001_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[32] = 4.0E0*I_NAI_Dxz_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[33] = 4.0E0*I_NAI_F2xz_S_C1_aa-2.0E0*1*I_NAI_Pz_S_C1_a;
  abcd[34] = 4.0E0*I_NAI_Fxyz_S_C1_aa;
  abcd[35] = 4.0E0*I_NAI_Fx2z_S_C1_aa-2.0E0*1*I_NAI_Px_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[36] = 4.0E0*I_NAI_Dxz_Px_C1000_aa;
  abcd[40] = 4.0E0*I_NAI_Dxz_Py_C1000_aa;
  abcd[44] = 4.0E0*I_NAI_Dxz_Pz_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[37] = 4.0E0*I_NAI_F2xz_Px_C1001_aa-2.0E0*1*I_NAI_Pz_Px_C1001_a;
  abcd[38] = 4.0E0*I_NAI_Fxyz_Px_C1001_aa;
  abcd[39] = 4.0E0*I_NAI_Fx2z_Px_C1001_aa-2.0E0*1*I_NAI_Px_Px_C1001_a;
  abcd[41] = 4.0E0*I_NAI_F2xz_Py_C1001_aa-2.0E0*1*I_NAI_Pz_Py_C1001_a;
  abcd[42] = 4.0E0*I_NAI_Fxyz_Py_C1001_aa;
  abcd[43] = 4.0E0*I_NAI_Fx2z_Py_C1001_aa-2.0E0*1*I_NAI_Px_Py_C1001_a;
  abcd[45] = 4.0E0*I_NAI_F2xz_Pz_C1001_aa-2.0E0*1*I_NAI_Pz_Pz_C1001_a;
  abcd[46] = 4.0E0*I_NAI_Fxyz_Pz_C1001_aa;
  abcd[47] = 4.0E0*I_NAI_Fx2z_Pz_C1001_aa-2.0E0*1*I_NAI_Px_Pz_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[48] = 4.0E0*I_NAI_D2y_S_C0_aa-2.0E0*1*I_NAI_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[49] = 4.0E0*I_NAI_Fx2y_S_C1_aa-2.0E0*1*I_NAI_Px_S_C1_a;
  abcd[50] = 4.0E0*I_NAI_F3y_S_C1_aa-2.0E0*1*I_NAI_Py_S_C1_a-2.0E0*2*I_NAI_Py_S_C1_a;
  abcd[51] = 4.0E0*I_NAI_F2yz_S_C1_aa-2.0E0*1*I_NAI_Pz_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[52] = 4.0E0*I_NAI_D2y_Px_C1000_aa-2.0E0*1*I_NAI_S_Px_C1000_a;
  abcd[56] = 4.0E0*I_NAI_D2y_Py_C1000_aa-2.0E0*1*I_NAI_S_Py_C1000_a;
  abcd[60] = 4.0E0*I_NAI_D2y_Pz_C1000_aa-2.0E0*1*I_NAI_S_Pz_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[53] = 4.0E0*I_NAI_Fx2y_Px_C1001_aa-2.0E0*1*I_NAI_Px_Px_C1001_a;
  abcd[54] = 4.0E0*I_NAI_F3y_Px_C1001_aa-2.0E0*1*I_NAI_Py_Px_C1001_a-2.0E0*2*I_NAI_Py_Px_C1001_a;
  abcd[55] = 4.0E0*I_NAI_F2yz_Px_C1001_aa-2.0E0*1*I_NAI_Pz_Px_C1001_a;
  abcd[57] = 4.0E0*I_NAI_Fx2y_Py_C1001_aa-2.0E0*1*I_NAI_Px_Py_C1001_a;
  abcd[58] = 4.0E0*I_NAI_F3y_Py_C1001_aa-2.0E0*1*I_NAI_Py_Py_C1001_a-2.0E0*2*I_NAI_Py_Py_C1001_a;
  abcd[59] = 4.0E0*I_NAI_F2yz_Py_C1001_aa-2.0E0*1*I_NAI_Pz_Py_C1001_a;
  abcd[61] = 4.0E0*I_NAI_Fx2y_Pz_C1001_aa-2.0E0*1*I_NAI_Px_Pz_C1001_a;
  abcd[62] = 4.0E0*I_NAI_F3y_Pz_C1001_aa-2.0E0*1*I_NAI_Py_Pz_C1001_a-2.0E0*2*I_NAI_Py_Pz_C1001_a;
  abcd[63] = 4.0E0*I_NAI_F2yz_Pz_C1001_aa-2.0E0*1*I_NAI_Pz_Pz_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[64] = 4.0E0*I_NAI_Dyz_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[65] = 4.0E0*I_NAI_Fxyz_S_C1_aa;
  abcd[66] = 4.0E0*I_NAI_F2yz_S_C1_aa-2.0E0*1*I_NAI_Pz_S_C1_a;
  abcd[67] = 4.0E0*I_NAI_Fy2z_S_C1_aa-2.0E0*1*I_NAI_Py_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[68] = 4.0E0*I_NAI_Dyz_Px_C1000_aa;
  abcd[72] = 4.0E0*I_NAI_Dyz_Py_C1000_aa;
  abcd[76] = 4.0E0*I_NAI_Dyz_Pz_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[69] = 4.0E0*I_NAI_Fxyz_Px_C1001_aa;
  abcd[70] = 4.0E0*I_NAI_F2yz_Px_C1001_aa-2.0E0*1*I_NAI_Pz_Px_C1001_a;
  abcd[71] = 4.0E0*I_NAI_Fy2z_Px_C1001_aa-2.0E0*1*I_NAI_Py_Px_C1001_a;
  abcd[73] = 4.0E0*I_NAI_Fxyz_Py_C1001_aa;
  abcd[74] = 4.0E0*I_NAI_F2yz_Py_C1001_aa-2.0E0*1*I_NAI_Pz_Py_C1001_a;
  abcd[75] = 4.0E0*I_NAI_Fy2z_Py_C1001_aa-2.0E0*1*I_NAI_Py_Py_C1001_a;
  abcd[77] = 4.0E0*I_NAI_Fxyz_Pz_C1001_aa;
  abcd[78] = 4.0E0*I_NAI_F2yz_Pz_C1001_aa-2.0E0*1*I_NAI_Pz_Pz_C1001_a;
  abcd[79] = 4.0E0*I_NAI_Fy2z_Pz_C1001_aa-2.0E0*1*I_NAI_Py_Pz_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C0_aa
   * RHS shell quartet name: SQ_NAI_S_S_C0_a
   ************************************************************/
  abcd[80] = 4.0E0*I_NAI_D2z_S_C0_aa-2.0E0*1*I_NAI_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1_aa
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   * RHS shell quartet name: SQ_NAI_P_S_C1_a
   ************************************************************/
  abcd[81] = 4.0E0*I_NAI_Fx2z_S_C1_aa-2.0E0*1*I_NAI_Px_S_C1_a;
  abcd[82] = 4.0E0*I_NAI_Fy2z_S_C1_aa-2.0E0*1*I_NAI_Py_S_C1_a;
  abcd[83] = 4.0E0*I_NAI_F3z_S_C1_aa-2.0E0*1*I_NAI_Pz_S_C1_a-2.0E0*2*I_NAI_Pz_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1000_aa
   * RHS shell quartet name: SQ_NAI_S_P_C1000_a
   ************************************************************/
  abcd[84] = 4.0E0*I_NAI_D2z_Px_C1000_aa-2.0E0*1*I_NAI_S_Px_C1000_a;
  abcd[88] = 4.0E0*I_NAI_D2z_Py_C1000_aa-2.0E0*1*I_NAI_S_Py_C1000_a;
  abcd[92] = 4.0E0*I_NAI_D2z_Pz_C1000_aa-2.0E0*1*I_NAI_S_Pz_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1001_aa
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   * RHS shell quartet name: SQ_NAI_P_P_C1001_a
   ************************************************************/
  abcd[85] = 4.0E0*I_NAI_Fx2z_Px_C1001_aa-2.0E0*1*I_NAI_Px_Px_C1001_a;
  abcd[86] = 4.0E0*I_NAI_Fy2z_Px_C1001_aa-2.0E0*1*I_NAI_Py_Px_C1001_a;
  abcd[87] = 4.0E0*I_NAI_F3z_Px_C1001_aa-2.0E0*1*I_NAI_Pz_Px_C1001_a-2.0E0*2*I_NAI_Pz_Px_C1001_a;
  abcd[89] = 4.0E0*I_NAI_Fx2z_Py_C1001_aa-2.0E0*1*I_NAI_Px_Py_C1001_a;
  abcd[90] = 4.0E0*I_NAI_Fy2z_Py_C1001_aa-2.0E0*1*I_NAI_Py_Py_C1001_a;
  abcd[91] = 4.0E0*I_NAI_F3z_Py_C1001_aa-2.0E0*1*I_NAI_Pz_Py_C1001_a-2.0E0*2*I_NAI_Pz_Py_C1001_a;
  abcd[93] = 4.0E0*I_NAI_Fx2z_Pz_C1001_aa-2.0E0*1*I_NAI_Px_Pz_C1001_a;
  abcd[94] = 4.0E0*I_NAI_Fy2z_Pz_C1001_aa-2.0E0*1*I_NAI_Py_Pz_C1001_a;
  abcd[95] = 4.0E0*I_NAI_F3z_Pz_C1001_aa-2.0E0*1*I_NAI_Pz_Pz_C1001_a-2.0E0*2*I_NAI_Pz_Pz_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[96] = 4.0E0*I_NAI_Px_Px_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[97] = 4.0E0*I_NAI_D2x_Px_C1_ab-2.0E0*1*I_NAI_S_Px_C1_b;
  abcd[98] = 4.0E0*I_NAI_Dxy_Px_C1_ab;
  abcd[99] = 4.0E0*I_NAI_Dxz_Px_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[100] = 4.0E0*I_NAI_Px_D2x_C1000_ab-2.0E0*1*I_NAI_Px_S_C1000_a;
  abcd[104] = 4.0E0*I_NAI_Px_Dxy_C1000_ab;
  abcd[108] = 4.0E0*I_NAI_Px_Dxz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[101] = 4.0E0*I_NAI_D2x_D2x_C1001_ab-2.0E0*1*I_NAI_D2x_S_C1001_a-2.0E0*1*I_NAI_S_D2x_C1001_b+1*I_NAI_S_S_C1001;
  abcd[102] = 4.0E0*I_NAI_Dxy_D2x_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[103] = 4.0E0*I_NAI_Dxz_D2x_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;
  abcd[105] = 4.0E0*I_NAI_D2x_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[106] = 4.0E0*I_NAI_Dxy_Dxy_C1001_ab;
  abcd[107] = 4.0E0*I_NAI_Dxz_Dxy_C1001_ab;
  abcd[109] = 4.0E0*I_NAI_D2x_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;
  abcd[110] = 4.0E0*I_NAI_Dxy_Dxz_C1001_ab;
  abcd[111] = 4.0E0*I_NAI_Dxz_Dxz_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[112] = 4.0E0*I_NAI_Px_Py_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[113] = 4.0E0*I_NAI_D2x_Py_C1_ab-2.0E0*1*I_NAI_S_Py_C1_b;
  abcd[114] = 4.0E0*I_NAI_Dxy_Py_C1_ab;
  abcd[115] = 4.0E0*I_NAI_Dxz_Py_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[116] = 4.0E0*I_NAI_Px_Dxy_C1000_ab;
  abcd[120] = 4.0E0*I_NAI_Px_D2y_C1000_ab-2.0E0*1*I_NAI_Px_S_C1000_a;
  abcd[124] = 4.0E0*I_NAI_Px_Dyz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[117] = 4.0E0*I_NAI_D2x_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[118] = 4.0E0*I_NAI_Dxy_Dxy_C1001_ab;
  abcd[119] = 4.0E0*I_NAI_Dxz_Dxy_C1001_ab;
  abcd[121] = 4.0E0*I_NAI_D2x_D2y_C1001_ab-2.0E0*1*I_NAI_D2x_S_C1001_a-2.0E0*1*I_NAI_S_D2y_C1001_b+1*I_NAI_S_S_C1001;
  abcd[122] = 4.0E0*I_NAI_Dxy_D2y_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[123] = 4.0E0*I_NAI_Dxz_D2y_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;
  abcd[125] = 4.0E0*I_NAI_D2x_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;
  abcd[126] = 4.0E0*I_NAI_Dxy_Dyz_C1001_ab;
  abcd[127] = 4.0E0*I_NAI_Dxz_Dyz_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[128] = 4.0E0*I_NAI_Px_Pz_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[129] = 4.0E0*I_NAI_D2x_Pz_C1_ab-2.0E0*1*I_NAI_S_Pz_C1_b;
  abcd[130] = 4.0E0*I_NAI_Dxy_Pz_C1_ab;
  abcd[131] = 4.0E0*I_NAI_Dxz_Pz_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[132] = 4.0E0*I_NAI_Px_Dxz_C1000_ab;
  abcd[136] = 4.0E0*I_NAI_Px_Dyz_C1000_ab;
  abcd[140] = 4.0E0*I_NAI_Px_D2z_C1000_ab-2.0E0*1*I_NAI_Px_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[133] = 4.0E0*I_NAI_D2x_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;
  abcd[134] = 4.0E0*I_NAI_Dxy_Dxz_C1001_ab;
  abcd[135] = 4.0E0*I_NAI_Dxz_Dxz_C1001_ab;
  abcd[137] = 4.0E0*I_NAI_D2x_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;
  abcd[138] = 4.0E0*I_NAI_Dxy_Dyz_C1001_ab;
  abcd[139] = 4.0E0*I_NAI_Dxz_Dyz_C1001_ab;
  abcd[141] = 4.0E0*I_NAI_D2x_D2z_C1001_ab-2.0E0*1*I_NAI_D2x_S_C1001_a-2.0E0*1*I_NAI_S_D2z_C1001_b+1*I_NAI_S_S_C1001;
  abcd[142] = 4.0E0*I_NAI_Dxy_D2z_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[143] = 4.0E0*I_NAI_Dxz_D2z_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[144] = 4.0E0*I_NAI_Py_Px_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[145] = 4.0E0*I_NAI_Dxy_Px_C1_ab;
  abcd[146] = 4.0E0*I_NAI_D2y_Px_C1_ab-2.0E0*1*I_NAI_S_Px_C1_b;
  abcd[147] = 4.0E0*I_NAI_Dyz_Px_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[148] = 4.0E0*I_NAI_Py_D2x_C1000_ab-2.0E0*1*I_NAI_Py_S_C1000_a;
  abcd[152] = 4.0E0*I_NAI_Py_Dxy_C1000_ab;
  abcd[156] = 4.0E0*I_NAI_Py_Dxz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[149] = 4.0E0*I_NAI_Dxy_D2x_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[150] = 4.0E0*I_NAI_D2y_D2x_C1001_ab-2.0E0*1*I_NAI_D2y_S_C1001_a-2.0E0*1*I_NAI_S_D2x_C1001_b+1*I_NAI_S_S_C1001;
  abcd[151] = 4.0E0*I_NAI_Dyz_D2x_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;
  abcd[153] = 4.0E0*I_NAI_Dxy_Dxy_C1001_ab;
  abcd[154] = 4.0E0*I_NAI_D2y_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[155] = 4.0E0*I_NAI_Dyz_Dxy_C1001_ab;
  abcd[157] = 4.0E0*I_NAI_Dxy_Dxz_C1001_ab;
  abcd[158] = 4.0E0*I_NAI_D2y_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;
  abcd[159] = 4.0E0*I_NAI_Dyz_Dxz_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[160] = 4.0E0*I_NAI_Py_Py_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[161] = 4.0E0*I_NAI_Dxy_Py_C1_ab;
  abcd[162] = 4.0E0*I_NAI_D2y_Py_C1_ab-2.0E0*1*I_NAI_S_Py_C1_b;
  abcd[163] = 4.0E0*I_NAI_Dyz_Py_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[164] = 4.0E0*I_NAI_Py_Dxy_C1000_ab;
  abcd[168] = 4.0E0*I_NAI_Py_D2y_C1000_ab-2.0E0*1*I_NAI_Py_S_C1000_a;
  abcd[172] = 4.0E0*I_NAI_Py_Dyz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[165] = 4.0E0*I_NAI_Dxy_Dxy_C1001_ab;
  abcd[166] = 4.0E0*I_NAI_D2y_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[167] = 4.0E0*I_NAI_Dyz_Dxy_C1001_ab;
  abcd[169] = 4.0E0*I_NAI_Dxy_D2y_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[170] = 4.0E0*I_NAI_D2y_D2y_C1001_ab-2.0E0*1*I_NAI_D2y_S_C1001_a-2.0E0*1*I_NAI_S_D2y_C1001_b+1*I_NAI_S_S_C1001;
  abcd[171] = 4.0E0*I_NAI_Dyz_D2y_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;
  abcd[173] = 4.0E0*I_NAI_Dxy_Dyz_C1001_ab;
  abcd[174] = 4.0E0*I_NAI_D2y_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;
  abcd[175] = 4.0E0*I_NAI_Dyz_Dyz_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[176] = 4.0E0*I_NAI_Py_Pz_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[177] = 4.0E0*I_NAI_Dxy_Pz_C1_ab;
  abcd[178] = 4.0E0*I_NAI_D2y_Pz_C1_ab-2.0E0*1*I_NAI_S_Pz_C1_b;
  abcd[179] = 4.0E0*I_NAI_Dyz_Pz_C1_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[180] = 4.0E0*I_NAI_Py_Dxz_C1000_ab;
  abcd[184] = 4.0E0*I_NAI_Py_Dyz_C1000_ab;
  abcd[188] = 4.0E0*I_NAI_Py_D2z_C1000_ab-2.0E0*1*I_NAI_Py_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[181] = 4.0E0*I_NAI_Dxy_Dxz_C1001_ab;
  abcd[182] = 4.0E0*I_NAI_D2y_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;
  abcd[183] = 4.0E0*I_NAI_Dyz_Dxz_C1001_ab;
  abcd[185] = 4.0E0*I_NAI_Dxy_Dyz_C1001_ab;
  abcd[186] = 4.0E0*I_NAI_D2y_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;
  abcd[187] = 4.0E0*I_NAI_Dyz_Dyz_C1001_ab;
  abcd[189] = 4.0E0*I_NAI_Dxy_D2z_C1001_ab-2.0E0*1*I_NAI_Dxy_S_C1001_a;
  abcd[190] = 4.0E0*I_NAI_D2y_D2z_C1001_ab-2.0E0*1*I_NAI_D2y_S_C1001_a-2.0E0*1*I_NAI_S_D2z_C1001_b+1*I_NAI_S_S_C1001;
  abcd[191] = 4.0E0*I_NAI_Dyz_D2z_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[192] = 4.0E0*I_NAI_Pz_Px_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[193] = 4.0E0*I_NAI_Dxz_Px_C1_ab;
  abcd[194] = 4.0E0*I_NAI_Dyz_Px_C1_ab;
  abcd[195] = 4.0E0*I_NAI_D2z_Px_C1_ab-2.0E0*1*I_NAI_S_Px_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[196] = 4.0E0*I_NAI_Pz_D2x_C1000_ab-2.0E0*1*I_NAI_Pz_S_C1000_a;
  abcd[200] = 4.0E0*I_NAI_Pz_Dxy_C1000_ab;
  abcd[204] = 4.0E0*I_NAI_Pz_Dxz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[197] = 4.0E0*I_NAI_Dxz_D2x_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;
  abcd[198] = 4.0E0*I_NAI_Dyz_D2x_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;
  abcd[199] = 4.0E0*I_NAI_D2z_D2x_C1001_ab-2.0E0*1*I_NAI_D2z_S_C1001_a-2.0E0*1*I_NAI_S_D2x_C1001_b+1*I_NAI_S_S_C1001;
  abcd[201] = 4.0E0*I_NAI_Dxz_Dxy_C1001_ab;
  abcd[202] = 4.0E0*I_NAI_Dyz_Dxy_C1001_ab;
  abcd[203] = 4.0E0*I_NAI_D2z_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[205] = 4.0E0*I_NAI_Dxz_Dxz_C1001_ab;
  abcd[206] = 4.0E0*I_NAI_Dyz_Dxz_C1001_ab;
  abcd[207] = 4.0E0*I_NAI_D2z_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[208] = 4.0E0*I_NAI_Pz_Py_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[209] = 4.0E0*I_NAI_Dxz_Py_C1_ab;
  abcd[210] = 4.0E0*I_NAI_Dyz_Py_C1_ab;
  abcd[211] = 4.0E0*I_NAI_D2z_Py_C1_ab-2.0E0*1*I_NAI_S_Py_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[212] = 4.0E0*I_NAI_Pz_Dxy_C1000_ab;
  abcd[216] = 4.0E0*I_NAI_Pz_D2y_C1000_ab-2.0E0*1*I_NAI_Pz_S_C1000_a;
  abcd[220] = 4.0E0*I_NAI_Pz_Dyz_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[213] = 4.0E0*I_NAI_Dxz_Dxy_C1001_ab;
  abcd[214] = 4.0E0*I_NAI_Dyz_Dxy_C1001_ab;
  abcd[215] = 4.0E0*I_NAI_D2z_Dxy_C1001_ab-2.0E0*1*I_NAI_S_Dxy_C1001_b;
  abcd[217] = 4.0E0*I_NAI_Dxz_D2y_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;
  abcd[218] = 4.0E0*I_NAI_Dyz_D2y_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;
  abcd[219] = 4.0E0*I_NAI_D2z_D2y_C1001_ab-2.0E0*1*I_NAI_D2z_S_C1001_a-2.0E0*1*I_NAI_S_D2y_C1001_b+1*I_NAI_S_S_C1001;
  abcd[221] = 4.0E0*I_NAI_Dxz_Dyz_C1001_ab;
  abcd[222] = 4.0E0*I_NAI_Dyz_Dyz_C1001_ab;
  abcd[223] = 4.0E0*I_NAI_D2z_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_P_C0_ab
   ************************************************************/
  abcd[224] = 4.0E0*I_NAI_Pz_Pz_C0_ab;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C1_ab
   * RHS shell quartet name: SQ_NAI_S_P_C1_b
   ************************************************************/
  abcd[225] = 4.0E0*I_NAI_Dxz_Pz_C1_ab;
  abcd[226] = 4.0E0*I_NAI_Dyz_Pz_C1_ab;
  abcd[227] = 4.0E0*I_NAI_D2z_Pz_C1_ab-2.0E0*1*I_NAI_S_Pz_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1000_ab
   * RHS shell quartet name: SQ_NAI_P_S_C1000_a
   ************************************************************/
  abcd[228] = 4.0E0*I_NAI_Pz_Dxz_C1000_ab;
  abcd[232] = 4.0E0*I_NAI_Pz_Dyz_C1000_ab;
  abcd[236] = 4.0E0*I_NAI_Pz_D2z_C1000_ab-2.0E0*1*I_NAI_Pz_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1001_ab
   * RHS shell quartet name: SQ_NAI_D_S_C1001_a
   * RHS shell quartet name: SQ_NAI_S_D_C1001_b
   * RHS shell quartet name: SQ_NAI_S_S_C1001
   ************************************************************/
  abcd[229] = 4.0E0*I_NAI_Dxz_Dxz_C1001_ab;
  abcd[230] = 4.0E0*I_NAI_Dyz_Dxz_C1001_ab;
  abcd[231] = 4.0E0*I_NAI_D2z_Dxz_C1001_ab-2.0E0*1*I_NAI_S_Dxz_C1001_b;
  abcd[233] = 4.0E0*I_NAI_Dxz_Dyz_C1001_ab;
  abcd[234] = 4.0E0*I_NAI_Dyz_Dyz_C1001_ab;
  abcd[235] = 4.0E0*I_NAI_D2z_Dyz_C1001_ab-2.0E0*1*I_NAI_S_Dyz_C1001_b;
  abcd[237] = 4.0E0*I_NAI_Dxz_D2z_C1001_ab-2.0E0*1*I_NAI_Dxz_S_C1001_a;
  abcd[238] = 4.0E0*I_NAI_Dyz_D2z_C1001_ab-2.0E0*1*I_NAI_Dyz_S_C1001_a;
  abcd[239] = 4.0E0*I_NAI_D2z_D2z_C1001_ab-2.0E0*1*I_NAI_D2z_S_C1001_a-2.0E0*1*I_NAI_S_D2z_C1001_b+1*I_NAI_S_S_C1001;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[240] = 4.0E0*I_NAI_S_D2x_C0_bb-2.0E0*1*I_NAI_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[241] = 4.0E0*I_NAI_Px_D2x_C1_bb-2.0E0*1*I_NAI_Px_S_C1_b;
  abcd[242] = 4.0E0*I_NAI_Py_D2x_C1_bb-2.0E0*1*I_NAI_Py_S_C1_b;
  abcd[243] = 4.0E0*I_NAI_Pz_D2x_C1_bb-2.0E0*1*I_NAI_Pz_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[244] = 4.0E0*I_NAI_S_F3x_C1000_bb-2.0E0*1*I_NAI_S_Px_C1000_b-2.0E0*2*I_NAI_S_Px_C1000_b;
  abcd[248] = 4.0E0*I_NAI_S_F2xy_C1000_bb-2.0E0*1*I_NAI_S_Py_C1000_b;
  abcd[252] = 4.0E0*I_NAI_S_F2xz_C1000_bb-2.0E0*1*I_NAI_S_Pz_C1000_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[245] = 4.0E0*I_NAI_Px_F3x_C1001_bb-2.0E0*1*I_NAI_Px_Px_C1001_b-2.0E0*2*I_NAI_Px_Px_C1001_b;
  abcd[246] = 4.0E0*I_NAI_Py_F3x_C1001_bb-2.0E0*1*I_NAI_Py_Px_C1001_b-2.0E0*2*I_NAI_Py_Px_C1001_b;
  abcd[247] = 4.0E0*I_NAI_Pz_F3x_C1001_bb-2.0E0*1*I_NAI_Pz_Px_C1001_b-2.0E0*2*I_NAI_Pz_Px_C1001_b;
  abcd[249] = 4.0E0*I_NAI_Px_F2xy_C1001_bb-2.0E0*1*I_NAI_Px_Py_C1001_b;
  abcd[250] = 4.0E0*I_NAI_Py_F2xy_C1001_bb-2.0E0*1*I_NAI_Py_Py_C1001_b;
  abcd[251] = 4.0E0*I_NAI_Pz_F2xy_C1001_bb-2.0E0*1*I_NAI_Pz_Py_C1001_b;
  abcd[253] = 4.0E0*I_NAI_Px_F2xz_C1001_bb-2.0E0*1*I_NAI_Px_Pz_C1001_b;
  abcd[254] = 4.0E0*I_NAI_Py_F2xz_C1001_bb-2.0E0*1*I_NAI_Py_Pz_C1001_b;
  abcd[255] = 4.0E0*I_NAI_Pz_F2xz_C1001_bb-2.0E0*1*I_NAI_Pz_Pz_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[256] = 4.0E0*I_NAI_S_Dxy_C0_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[257] = 4.0E0*I_NAI_Px_Dxy_C1_bb;
  abcd[258] = 4.0E0*I_NAI_Py_Dxy_C1_bb;
  abcd[259] = 4.0E0*I_NAI_Pz_Dxy_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[260] = 4.0E0*I_NAI_S_F2xy_C1000_bb-2.0E0*1*I_NAI_S_Py_C1000_b;
  abcd[264] = 4.0E0*I_NAI_S_Fx2y_C1000_bb-2.0E0*1*I_NAI_S_Px_C1000_b;
  abcd[268] = 4.0E0*I_NAI_S_Fxyz_C1000_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[261] = 4.0E0*I_NAI_Px_F2xy_C1001_bb-2.0E0*1*I_NAI_Px_Py_C1001_b;
  abcd[262] = 4.0E0*I_NAI_Py_F2xy_C1001_bb-2.0E0*1*I_NAI_Py_Py_C1001_b;
  abcd[263] = 4.0E0*I_NAI_Pz_F2xy_C1001_bb-2.0E0*1*I_NAI_Pz_Py_C1001_b;
  abcd[265] = 4.0E0*I_NAI_Px_Fx2y_C1001_bb-2.0E0*1*I_NAI_Px_Px_C1001_b;
  abcd[266] = 4.0E0*I_NAI_Py_Fx2y_C1001_bb-2.0E0*1*I_NAI_Py_Px_C1001_b;
  abcd[267] = 4.0E0*I_NAI_Pz_Fx2y_C1001_bb-2.0E0*1*I_NAI_Pz_Px_C1001_b;
  abcd[269] = 4.0E0*I_NAI_Px_Fxyz_C1001_bb;
  abcd[270] = 4.0E0*I_NAI_Py_Fxyz_C1001_bb;
  abcd[271] = 4.0E0*I_NAI_Pz_Fxyz_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[272] = 4.0E0*I_NAI_S_Dxz_C0_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[273] = 4.0E0*I_NAI_Px_Dxz_C1_bb;
  abcd[274] = 4.0E0*I_NAI_Py_Dxz_C1_bb;
  abcd[275] = 4.0E0*I_NAI_Pz_Dxz_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[276] = 4.0E0*I_NAI_S_F2xz_C1000_bb-2.0E0*1*I_NAI_S_Pz_C1000_b;
  abcd[280] = 4.0E0*I_NAI_S_Fxyz_C1000_bb;
  abcd[284] = 4.0E0*I_NAI_S_Fx2z_C1000_bb-2.0E0*1*I_NAI_S_Px_C1000_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[277] = 4.0E0*I_NAI_Px_F2xz_C1001_bb-2.0E0*1*I_NAI_Px_Pz_C1001_b;
  abcd[278] = 4.0E0*I_NAI_Py_F2xz_C1001_bb-2.0E0*1*I_NAI_Py_Pz_C1001_b;
  abcd[279] = 4.0E0*I_NAI_Pz_F2xz_C1001_bb-2.0E0*1*I_NAI_Pz_Pz_C1001_b;
  abcd[281] = 4.0E0*I_NAI_Px_Fxyz_C1001_bb;
  abcd[282] = 4.0E0*I_NAI_Py_Fxyz_C1001_bb;
  abcd[283] = 4.0E0*I_NAI_Pz_Fxyz_C1001_bb;
  abcd[285] = 4.0E0*I_NAI_Px_Fx2z_C1001_bb-2.0E0*1*I_NAI_Px_Px_C1001_b;
  abcd[286] = 4.0E0*I_NAI_Py_Fx2z_C1001_bb-2.0E0*1*I_NAI_Py_Px_C1001_b;
  abcd[287] = 4.0E0*I_NAI_Pz_Fx2z_C1001_bb-2.0E0*1*I_NAI_Pz_Px_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[288] = 4.0E0*I_NAI_S_D2y_C0_bb-2.0E0*1*I_NAI_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[289] = 4.0E0*I_NAI_Px_D2y_C1_bb-2.0E0*1*I_NAI_Px_S_C1_b;
  abcd[290] = 4.0E0*I_NAI_Py_D2y_C1_bb-2.0E0*1*I_NAI_Py_S_C1_b;
  abcd[291] = 4.0E0*I_NAI_Pz_D2y_C1_bb-2.0E0*1*I_NAI_Pz_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[292] = 4.0E0*I_NAI_S_Fx2y_C1000_bb-2.0E0*1*I_NAI_S_Px_C1000_b;
  abcd[296] = 4.0E0*I_NAI_S_F3y_C1000_bb-2.0E0*1*I_NAI_S_Py_C1000_b-2.0E0*2*I_NAI_S_Py_C1000_b;
  abcd[300] = 4.0E0*I_NAI_S_F2yz_C1000_bb-2.0E0*1*I_NAI_S_Pz_C1000_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[293] = 4.0E0*I_NAI_Px_Fx2y_C1001_bb-2.0E0*1*I_NAI_Px_Px_C1001_b;
  abcd[294] = 4.0E0*I_NAI_Py_Fx2y_C1001_bb-2.0E0*1*I_NAI_Py_Px_C1001_b;
  abcd[295] = 4.0E0*I_NAI_Pz_Fx2y_C1001_bb-2.0E0*1*I_NAI_Pz_Px_C1001_b;
  abcd[297] = 4.0E0*I_NAI_Px_F3y_C1001_bb-2.0E0*1*I_NAI_Px_Py_C1001_b-2.0E0*2*I_NAI_Px_Py_C1001_b;
  abcd[298] = 4.0E0*I_NAI_Py_F3y_C1001_bb-2.0E0*1*I_NAI_Py_Py_C1001_b-2.0E0*2*I_NAI_Py_Py_C1001_b;
  abcd[299] = 4.0E0*I_NAI_Pz_F3y_C1001_bb-2.0E0*1*I_NAI_Pz_Py_C1001_b-2.0E0*2*I_NAI_Pz_Py_C1001_b;
  abcd[301] = 4.0E0*I_NAI_Px_F2yz_C1001_bb-2.0E0*1*I_NAI_Px_Pz_C1001_b;
  abcd[302] = 4.0E0*I_NAI_Py_F2yz_C1001_bb-2.0E0*1*I_NAI_Py_Pz_C1001_b;
  abcd[303] = 4.0E0*I_NAI_Pz_F2yz_C1001_bb-2.0E0*1*I_NAI_Pz_Pz_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[304] = 4.0E0*I_NAI_S_Dyz_C0_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[305] = 4.0E0*I_NAI_Px_Dyz_C1_bb;
  abcd[306] = 4.0E0*I_NAI_Py_Dyz_C1_bb;
  abcd[307] = 4.0E0*I_NAI_Pz_Dyz_C1_bb;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[308] = 4.0E0*I_NAI_S_Fxyz_C1000_bb;
  abcd[312] = 4.0E0*I_NAI_S_F2yz_C1000_bb-2.0E0*1*I_NAI_S_Pz_C1000_b;
  abcd[316] = 4.0E0*I_NAI_S_Fy2z_C1000_bb-2.0E0*1*I_NAI_S_Py_C1000_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[309] = 4.0E0*I_NAI_Px_Fxyz_C1001_bb;
  abcd[310] = 4.0E0*I_NAI_Py_Fxyz_C1001_bb;
  abcd[311] = 4.0E0*I_NAI_Pz_Fxyz_C1001_bb;
  abcd[313] = 4.0E0*I_NAI_Px_F2yz_C1001_bb-2.0E0*1*I_NAI_Px_Pz_C1001_b;
  abcd[314] = 4.0E0*I_NAI_Py_F2yz_C1001_bb-2.0E0*1*I_NAI_Py_Pz_C1001_b;
  abcd[315] = 4.0E0*I_NAI_Pz_F2yz_C1001_bb-2.0E0*1*I_NAI_Pz_Pz_C1001_b;
  abcd[317] = 4.0E0*I_NAI_Px_Fy2z_C1001_bb-2.0E0*1*I_NAI_Px_Py_C1001_b;
  abcd[318] = 4.0E0*I_NAI_Py_Fy2z_C1001_bb-2.0E0*1*I_NAI_Py_Py_C1001_b;
  abcd[319] = 4.0E0*I_NAI_Pz_Fy2z_C1001_bb-2.0E0*1*I_NAI_Pz_Py_C1001_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_S_C0_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_D_C0_bb
   * RHS shell quartet name: SQ_NAI_S_S_C0_b
   ************************************************************/
  abcd[320] = 4.0E0*I_NAI_S_D2z_C0_bb-2.0E0*1*I_NAI_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_S_C1_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_D_C1_bb
   * RHS shell quartet name: SQ_NAI_P_S_C1_b
   ************************************************************/
  abcd[321] = 4.0E0*I_NAI_Px_D2z_C1_bb-2.0E0*1*I_NAI_Px_S_C1_b;
  abcd[322] = 4.0E0*I_NAI_Py_D2z_C1_bb-2.0E0*1*I_NAI_Py_S_C1_b;
  abcd[323] = 4.0E0*I_NAI_Pz_D2z_C1_bb-2.0E0*1*I_NAI_Pz_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_NAI_S_P_C1000_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_S_F_C1000_bb
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   * RHS shell quartet name: SQ_NAI_S_P_C1000_b
   ************************************************************/
  abcd[324] = 4.0E0*I_NAI_S_Fx2z_C1000_bb-2.0E0*1*I_NAI_S_Px_C1000_b;
  abcd[328] = 4.0E0*I_NAI_S_Fy2z_C1000_bb-2.0E0*1*I_NAI_S_Py_C1000_b;
  abcd[332] = 4.0E0*I_NAI_S_F3z_C1000_bb-2.0E0*1*I_NAI_S_Pz_C1000_b-2.0E0*2*I_NAI_S_Pz_C1000_b;

  /************************************************************
   * shell quartet name: SQ_NAI_P_P_C1001_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_P_F_C1001_bb
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   * RHS shell quartet name: SQ_NAI_P_P_C1001_b
   ************************************************************/
  abcd[325] = 4.0E0*I_NAI_Px_Fx2z_C1001_bb-2.0E0*1*I_NAI_Px_Px_C1001_b;
  abcd[326] = 4.0E0*I_NAI_Py_Fx2z_C1001_bb-2.0E0*1*I_NAI_Py_Px_C1001_b;
  abcd[327] = 4.0E0*I_NAI_Pz_Fx2z_C1001_bb-2.0E0*1*I_NAI_Pz_Px_C1001_b;
  abcd[329] = 4.0E0*I_NAI_Px_Fy2z_C1001_bb-2.0E0*1*I_NAI_Px_Py_C1001_b;
  abcd[330] = 4.0E0*I_NAI_Py_Fy2z_C1001_bb-2.0E0*1*I_NAI_Py_Py_C1001_b;
  abcd[331] = 4.0E0*I_NAI_Pz_Fy2z_C1001_bb-2.0E0*1*I_NAI_Pz_Py_C1001_b;
  abcd[333] = 4.0E0*I_NAI_Px_F3z_C1001_bb-2.0E0*1*I_NAI_Px_Pz_C1001_b-2.0E0*2*I_NAI_Px_Pz_C1001_b;
  abcd[334] = 4.0E0*I_NAI_Py_F3z_C1001_bb-2.0E0*1*I_NAI_Py_Pz_C1001_b-2.0E0*2*I_NAI_Py_Pz_C1001_b;
  abcd[335] = 4.0E0*I_NAI_Pz_F3z_C1001_bb-2.0E0*1*I_NAI_Pz_Pz_C1001_b-2.0E0*2*I_NAI_Pz_Pz_C1001_b;
}
