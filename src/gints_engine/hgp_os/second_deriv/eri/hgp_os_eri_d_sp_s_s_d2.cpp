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
// BRA1 as redundant position, total RHS integrals evaluated as: 54068
// BRA2 as redundant position, total RHS integrals evaluated as: 55248
// KET1 as redundant position, total RHS integrals evaluated as: 34226
// KET2 as redundant position, total RHS integrals evaluated as: 34226
// the redundant position is: KET2
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
// BRA1  KET1
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
// BRA2  KET1
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
// KET1  KET1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_eri_d_sp_s_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_G4x_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G4y_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_G4z_S_S_S_C2_aa = 0.0E0;
  Double I_ERI_D2x_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C2_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C2_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C2_a = 0.0E0;
  Double I_ERI_S_S_S_S_C2 = 0.0E0;
  Double I_ERI_S_Px_S_S_C1002 = 0.0E0;
  Double I_ERI_S_Py_S_S_C1002 = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1002 = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Px_S_S_S_C1002 = 0.0E0;
  Double I_ERI_Py_S_S_S_C1002 = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1002 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C2_ac = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C2_ac = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C2_ac = 0.0E0;
  Double I_ERI_Px_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_D2x_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C2_cc = 0.0E0;
  Double I_ERI_D2x_S_S_S_C2_c = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C2_c = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C2_c = 0.0E0;
  Double I_ERI_D2y_S_S_S_C2_c = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C2_c = 0.0E0;
  Double I_ERI_D2z_S_S_S_C2_c = 0.0E0;
  Double I_ERI_H5x_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H4xy_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H4xz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H5y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H4yz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_H5z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1002_aa = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G4x_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G4y_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_G4z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F3x_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_C2_ab = 0.0E0;
  Double I_ERI_Px_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1002_ac = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1002_ac = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1002_ac = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C2_bc = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C2_bc = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C2_bc = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1002_cc = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1002_cc = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1002_c = 0.0E0;
  Double I_ERI_H5x_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H4xy_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H4xz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H5y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H4yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_H5z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1002_ab = 0.0E0;
  Double I_ERI_Px_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G4x_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G4y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G4z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F3x_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F3y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_F3z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_D2x_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_D2y_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_D2z_S_S_S_C2_bb = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1002_bc = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1002_bc = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1002_bc = 0.0E0;
  Double I_ERI_H5x_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H4xy_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H4xz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H5y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H4yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_H5z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1002_bb = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1002_bb = 0.0E0;

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
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
    for(UInt jp2=0; jp2<jnp2; jp2++) {
      Double onede = jexp[jp2];
      Double eta   = 1.0E0/onede;
      Double ediff = jexpdiff[jp2];
      Double gamma = 0.5E0*(eta+ediff);
      Double delta = 0.5E0*(eta-ediff);
      Double jc2   = jcoe[jp2];
      Double fket  = jfac[jp2];
      Double pref      = fbra*fket;
      Double prefactor = pref;

      // 
      // here below the code is performing significance test for integrals on
      // primitive integrals. Here we use the overlap integrals to roughly 
      // estimate the order of the result integrals
      // the threshold value should be for primitive function quartet, we compare
      // the value against machine precision for significance test
      // 
      Double I_ERI_S_S_S_S_vrr_IntegralTest = pref;
      if (fabs(ic2*jc2)>1.0E0) {
        I_ERI_S_S_S_S_vrr_IntegralTest = prefactor;
      }

      // test the integrals with the pMax, which is the maximum value
      // of the corresponding density matrix block(or it may be maximum
      // value pair of the corresponding density matrix block)
      if(fabs(I_ERI_S_S_S_S_vrr_IntegralTest*pMax)<THRESHOLD_MATH) continue;
      isSignificant = true;


      UInt offsetQ  = 3*jp2;
      Double QX    = Q[offsetQ  ];
      Double QY    = Q[offsetQ+1];
      Double QZ    = Q[offsetQ+2];
      Double rho   = 1.0E0/(onedz+onede);
      Double sqrho = sqrt(rho);
      Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);
      Double u     = rho*PQ2;
      if (withErfR12) u = PQ2/(1.0E0/(omega*omega)+1.0E0/rho);
      Double squ   = sqrt(u);
      Double QCX   = QX - C[0];
      Double QCY   = QY - C[1];
      Double QCZ   = QZ - C[2];
      Double WX    = rho*(PX*onede + QX*onedz);
      Double WY    = rho*(PY*onede + QY*onedz);
      Double WZ    = rho*(PZ*onede + QZ*onedz);
      Double oned2k= 0.5E0*rho*onede*onedz;
      Double WPX   = WX - PX;
      Double WPY   = WY - PY;
      Double WPZ   = WZ - PZ;
      Double rhod2zsq = rho*oned2z*onedz;
      Double WQX   = WX - QX;
      Double WQY   = WY - QY;
      Double WQZ   = WZ - QZ;
      Double oned2e= 0.5E0*onede;
      Double rhod2esq= rho*oned2e*onede;


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

      Double I_ERI_S_S_S_S_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M1_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M2_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M3_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M4_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M5_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M5_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER45;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER17*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER15*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = 1.0E0+u2*ONEOVER13*I_ERI_S_S_S_S_M5_vrr;
        I_ERI_S_S_S_S_M5_vrr = ONEOVER11*I_ERI_S_S_S_S_M5_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M5_vrr  = f*I_ERI_S_S_S_S_M5_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M4_vrr  = ONEOVER9*(u2*I_ERI_S_S_S_S_M5_vrr+f);
        I_ERI_S_S_S_S_M3_vrr  = ONEOVER7*(u2*I_ERI_S_S_S_S_M4_vrr+f);
        I_ERI_S_S_S_S_M2_vrr  = ONEOVER5*(u2*I_ERI_S_S_S_S_M3_vrr+f);
        I_ERI_S_S_S_S_M1_vrr  = ONEOVER3*(u2*I_ERI_S_S_S_S_M2_vrr+f);
        I_ERI_S_S_S_S_vrr  = ONEOVER1*(u2*I_ERI_S_S_S_S_M1_vrr+f);

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
          I_ERI_S_S_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ERI_S_S_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ERI_S_S_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ERI_S_S_S_S_vrr_d-f);
        I_ERI_S_S_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ERI_S_S_S_S_M1_vrr_d-f);
        I_ERI_S_S_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ERI_S_S_S_S_M2_vrr_d-f);
        I_ERI_S_S_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ERI_S_S_S_S_M3_vrr_d-f);
        I_ERI_S_S_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ERI_S_S_S_S_M4_vrr_d-f);

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);
        I_ERI_S_S_S_S_M5_vrr = static_cast<Double>(I_ERI_S_S_S_S_M5_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ERI_S_S_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ERI_S_S_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M1_vrr = oneO2u*(1.0E0*I_ERI_S_S_S_S_vrr-f);
        I_ERI_S_S_S_S_M2_vrr = oneO2u*(3.0E0*I_ERI_S_S_S_S_M1_vrr-f);
        I_ERI_S_S_S_S_M3_vrr = oneO2u*(5.0E0*I_ERI_S_S_S_S_M2_vrr-f);
        I_ERI_S_S_S_S_M4_vrr = oneO2u*(7.0E0*I_ERI_S_S_S_S_M3_vrr-f);
        I_ERI_S_S_S_S_M5_vrr = oneO2u*(9.0E0*I_ERI_S_S_S_S_M4_vrr-f);

#endif

      }


      // now scale the bottom integral if oper in erf(r12)/r12 form
      if (withErfR12) {
        Double erfPref0   = 1.0E0+rho/(omega*omega);
        Double erfPref1   = 1.0E0/erfPref0;
        Double erfp       = sqrt(erfPref1);
        Double erfp2      = erfp*erfp;
        Double erfPref_1  = erfp;
        I_ERI_S_S_S_S_vrr = I_ERI_S_S_S_S_vrr*erfPref_1;
        Double erfPref_3 = erfPref_1*erfp2;
        Double erfPref_5 = erfPref_3*erfp2;
        Double erfPref_7 = erfPref_5*erfp2;
        Double erfPref_9 = erfPref_7*erfp2;
        Double erfPref_11 = erfPref_9*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_Px_S_S_S_M4_vrr = PAX*I_ERI_S_S_S_S_M4_vrr+WPX*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Py_S_S_S_M4_vrr = PAY*I_ERI_S_S_S_S_M4_vrr+WPY*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Pz_S_S_S_M4_vrr = PAZ*I_ERI_S_S_S_S_M4_vrr+WPZ*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_Px_S_S_S_M3_vrr = PAX*I_ERI_S_S_S_S_M3_vrr+WPX*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Py_S_S_S_M3_vrr = PAY*I_ERI_S_S_S_S_M3_vrr+WPY*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Pz_S_S_S_M3_vrr = PAZ*I_ERI_S_S_S_S_M3_vrr+WPZ*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M3_vrr = PAX*I_ERI_Px_S_S_S_M3_vrr+WPX*I_ERI_Px_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_S_S_M3_vrr = PAY*I_ERI_Px_S_S_S_M3_vrr+WPY*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_S_S_M3_vrr = PAY*I_ERI_Py_S_S_S_M3_vrr+WPY*I_ERI_Py_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_S_S_M3_vrr = PAZ*I_ERI_Pz_S_S_S_M3_vrr+WPZ*I_ERI_Pz_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_S_Px_S_M2_vrr = QCX*I_ERI_S_S_S_S_M2_vrr+WQX*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Py_S_M2_vrr = QCY*I_ERI_S_S_S_S_M2_vrr+WQY*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Pz_S_M2_vrr = QCZ*I_ERI_S_S_S_S_M2_vrr+WQZ*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_Px_S_S_S_M2_vrr = PAX*I_ERI_S_S_S_S_M2_vrr+WPX*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Py_S_S_S_M2_vrr = PAY*I_ERI_S_S_S_S_M2_vrr+WPY*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Pz_S_S_S_M2_vrr = PAZ*I_ERI_S_S_S_S_M2_vrr+WPZ*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M2_vrr = PAX*I_ERI_Px_S_S_S_M2_vrr+WPX*I_ERI_Px_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_S_S_M2_vrr = PAY*I_ERI_Px_S_S_S_M2_vrr+WPY*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_S_S_M2_vrr = PAZ*I_ERI_Px_S_S_S_M2_vrr+WPZ*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_S_S_M2_vrr = PAY*I_ERI_Py_S_S_S_M2_vrr+WPY*I_ERI_Py_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_S_S_M2_vrr = PAZ*I_ERI_Py_S_S_S_M2_vrr+WPZ*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_S_S_M2_vrr = PAZ*I_ERI_Pz_S_S_S_M2_vrr+WPZ*I_ERI_Pz_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M2_vrr = PAX*I_ERI_D2x_S_S_S_M2_vrr+WPX*I_ERI_D2x_S_S_S_M3_vrr+2*oned2z*I_ERI_Px_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_S_S_M2_vrr = PAY*I_ERI_D2x_S_S_S_M2_vrr+WPY*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_S_S_M2_vrr = PAZ*I_ERI_D2x_S_S_S_M2_vrr+WPZ*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_S_S_M2_vrr = PAX*I_ERI_D2y_S_S_S_M2_vrr+WPX*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_S_S_M2_vrr = PAZ*I_ERI_Dxy_S_S_S_M2_vrr+WPZ*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_S_S_M2_vrr = PAX*I_ERI_D2z_S_S_S_M2_vrr+WPX*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_S_S_M2_vrr = PAY*I_ERI_D2y_S_S_S_M2_vrr+WPY*I_ERI_D2y_S_S_S_M3_vrr+2*oned2z*I_ERI_Py_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_S_S_M2_vrr = PAZ*I_ERI_D2y_S_S_S_M2_vrr+WPZ*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_S_S_M2_vrr = PAY*I_ERI_D2z_S_S_S_M2_vrr+WPY*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_S_S_M2_vrr = PAZ*I_ERI_D2z_S_S_S_M2_vrr+WPZ*I_ERI_D2z_S_S_S_M3_vrr+2*oned2z*I_ERI_Pz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_Px_S_M1_vrr = QCX*I_ERI_S_S_S_S_M1_vrr+WQX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Py_S_M1_vrr = QCY*I_ERI_S_S_S_S_M1_vrr+WQY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Pz_S_M1_vrr = QCZ*I_ERI_S_S_S_S_M1_vrr+WQZ*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_Px_S_S_S_M1_vrr = PAX*I_ERI_S_S_S_S_M1_vrr+WPX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Py_S_S_S_M1_vrr = PAY*I_ERI_S_S_S_S_M1_vrr+WPY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_S_S_M1_vrr = PAZ*I_ERI_S_S_S_S_M1_vrr+WPZ*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M1_vrr = PAX*I_ERI_S_S_Px_S_M1_vrr+WPX*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Py_S_Px_S_M1_vrr = PAY*I_ERI_S_S_Px_S_M1_vrr+WPY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Px_S_M1_vrr = PAZ*I_ERI_S_S_Px_S_M1_vrr+WPZ*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_Py_S_M1_vrr = PAX*I_ERI_S_S_Py_S_M1_vrr+WPX*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Py_S_M1_vrr = PAY*I_ERI_S_S_Py_S_M1_vrr+WPY*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_Py_S_M1_vrr = PAZ*I_ERI_S_S_Py_S_M1_vrr+WPZ*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_Pz_S_M1_vrr = PAX*I_ERI_S_S_Pz_S_M1_vrr+WPX*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Py_S_Pz_S_M1_vrr = PAY*I_ERI_S_S_Pz_S_M1_vrr+WPY*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Pz_S_Pz_S_M1_vrr = PAZ*I_ERI_S_S_Pz_S_M1_vrr+WPZ*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M1_vrr = PAX*I_ERI_Px_S_S_S_M1_vrr+WPX*I_ERI_Px_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_S_S_M1_vrr = PAY*I_ERI_Px_S_S_S_M1_vrr+WPY*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_S_S_M1_vrr = PAZ*I_ERI_Px_S_S_S_M1_vrr+WPZ*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_S_S_M1_vrr = PAY*I_ERI_Py_S_S_S_M1_vrr+WPY*I_ERI_Py_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_S_S_M1_vrr = PAZ*I_ERI_Py_S_S_S_M1_vrr+WPZ*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_S_S_M1_vrr = PAZ*I_ERI_Pz_S_S_S_M1_vrr+WPZ*I_ERI_Pz_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M1_vrr = QCX*I_ERI_D2x_S_S_S_M1_vrr+WQX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Px_S_M1_vrr = QCX*I_ERI_Dxy_S_S_S_M1_vrr+WQX*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Px_S_M1_vrr = QCX*I_ERI_Dxz_S_S_S_M1_vrr+WQX*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Px_S_M1_vrr = QCX*I_ERI_D2y_S_S_S_M1_vrr+WQX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Px_S_M1_vrr = QCX*I_ERI_Dyz_S_S_S_M1_vrr+WQX*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Px_S_M1_vrr = QCX*I_ERI_D2z_S_S_S_M1_vrr+WQX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Py_S_M1_vrr = QCY*I_ERI_D2x_S_S_S_M1_vrr+WQY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Py_S_M1_vrr = QCY*I_ERI_Dxy_S_S_S_M1_vrr+WQY*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Py_S_M1_vrr = QCY*I_ERI_Dxz_S_S_S_M1_vrr+WQY*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Py_S_M1_vrr = QCY*I_ERI_D2y_S_S_S_M1_vrr+WQY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Py_S_M1_vrr = QCY*I_ERI_Dyz_S_S_S_M1_vrr+WQY*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Py_S_M1_vrr = QCY*I_ERI_D2z_S_S_S_M1_vrr+WQY*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Pz_S_M1_vrr = QCZ*I_ERI_D2x_S_S_S_M1_vrr+WQZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Pz_S_M1_vrr = QCZ*I_ERI_Dxy_S_S_S_M1_vrr+WQZ*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Pz_S_M1_vrr = QCZ*I_ERI_Dxz_S_S_S_M1_vrr+WQZ*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Pz_S_M1_vrr = QCZ*I_ERI_D2y_S_S_S_M1_vrr+WQZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Pz_S_M1_vrr = QCZ*I_ERI_Dyz_S_S_S_M1_vrr+WQZ*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Pz_S_M1_vrr = QCZ*I_ERI_D2z_S_S_S_M1_vrr+WQZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Pz_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M1_vrr = PAX*I_ERI_D2x_S_S_S_M1_vrr+WPX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2z*I_ERI_Px_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_S_S_M1_vrr = PAY*I_ERI_D2x_S_S_S_M1_vrr+WPY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_S_S_M1_vrr = PAZ*I_ERI_D2x_S_S_S_M1_vrr+WPZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_S_S_M1_vrr = PAX*I_ERI_D2y_S_S_S_M1_vrr+WPX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fxyz_S_S_S_M1_vrr = PAZ*I_ERI_Dxy_S_S_S_M1_vrr+WPZ*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_S_S_M1_vrr = PAX*I_ERI_D2z_S_S_S_M1_vrr+WPX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_S_S_M1_vrr = PAY*I_ERI_D2y_S_S_S_M1_vrr+WPY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2z*I_ERI_Py_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_S_S_M1_vrr = PAZ*I_ERI_D2y_S_S_S_M1_vrr+WPZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_S_S_M1_vrr = PAY*I_ERI_D2z_S_S_S_M1_vrr+WPY*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_S_S_M1_vrr = PAZ*I_ERI_D2z_S_S_S_M1_vrr+WPZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2z*I_ERI_Pz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M1_vrr = QCX*I_ERI_F3x_S_S_S_M1_vrr+WQX*I_ERI_F3x_S_S_S_M2_vrr+3*oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_Px_S_M1_vrr = QCX*I_ERI_F2xy_S_S_S_M1_vrr+WQX*I_ERI_F2xy_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_Px_S_M1_vrr = QCX*I_ERI_F2xz_S_S_S_M1_vrr+WQX*I_ERI_F2xz_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_Px_S_M1_vrr = QCX*I_ERI_Fx2y_S_S_S_M1_vrr+WQX*I_ERI_Fx2y_S_S_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fxyz_S_Px_S_M1_vrr = QCX*I_ERI_Fxyz_S_S_S_M1_vrr+WQX*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_Px_S_M1_vrr = QCX*I_ERI_Fx2z_S_S_S_M1_vrr+WQX*I_ERI_Fx2z_S_S_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_Px_S_M1_vrr = QCX*I_ERI_F3y_S_S_S_M1_vrr+WQX*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_Px_S_M1_vrr = QCX*I_ERI_F2yz_S_S_S_M1_vrr+WQX*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_Px_S_M1_vrr = QCX*I_ERI_Fy2z_S_S_S_M1_vrr+WQX*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_Px_S_M1_vrr = QCX*I_ERI_F3z_S_S_S_M1_vrr+WQX*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_F3x_S_Py_S_M1_vrr = QCY*I_ERI_F3x_S_S_S_M1_vrr+WQY*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_Py_S_M1_vrr = QCY*I_ERI_F2xy_S_S_S_M1_vrr+WQY*I_ERI_F2xy_S_S_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_Py_S_M1_vrr = QCY*I_ERI_F2xz_S_S_S_M1_vrr+WQY*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_Py_S_M1_vrr = QCY*I_ERI_Fx2y_S_S_S_M1_vrr+WQY*I_ERI_Fx2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Fxyz_S_Py_S_M1_vrr = QCY*I_ERI_Fxyz_S_S_S_M1_vrr+WQY*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_Py_S_M1_vrr = QCY*I_ERI_Fx2z_S_S_S_M1_vrr+WQY*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_Py_S_M1_vrr = QCY*I_ERI_F3y_S_S_S_M1_vrr+WQY*I_ERI_F3y_S_S_S_M2_vrr+3*oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_Py_S_M1_vrr = QCY*I_ERI_F2yz_S_S_S_M1_vrr+WQY*I_ERI_F2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_Py_S_M1_vrr = QCY*I_ERI_Fy2z_S_S_S_M1_vrr+WQY*I_ERI_Fy2z_S_S_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_Py_S_M1_vrr = QCY*I_ERI_F3z_S_S_S_M1_vrr+WQY*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_F3x_S_Pz_S_M1_vrr = QCZ*I_ERI_F3x_S_S_S_M1_vrr+WQZ*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_Pz_S_M1_vrr = QCZ*I_ERI_F2xy_S_S_S_M1_vrr+WQZ*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_Pz_S_M1_vrr = QCZ*I_ERI_F2xz_S_S_S_M1_vrr+WQZ*I_ERI_F2xz_S_S_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M1_vrr = QCZ*I_ERI_Fx2y_S_S_S_M1_vrr+WQZ*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Fxyz_S_Pz_S_M1_vrr = QCZ*I_ERI_Fxyz_S_S_S_M1_vrr+WQZ*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Fx2z_S_S_S_M1_vrr+WQZ*I_ERI_Fx2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_Pz_S_M1_vrr = QCZ*I_ERI_F3y_S_S_S_M1_vrr+WQZ*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_Pz_S_M1_vrr = QCZ*I_ERI_F2yz_S_S_S_M1_vrr+WQZ*I_ERI_F2yz_S_S_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Fy2z_S_S_S_M1_vrr+WQZ*I_ERI_Fy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_Pz_S_M1_vrr = QCZ*I_ERI_F3z_S_S_S_M1_vrr+WQZ*I_ERI_F3z_S_S_S_M2_vrr+3*oned2k*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M1_vrr = PAX*I_ERI_F3x_S_S_S_M1_vrr+WPX*I_ERI_F3x_S_S_S_M2_vrr+3*oned2z*I_ERI_D2x_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_S_S_M1_vrr = PAY*I_ERI_F3x_S_S_S_M1_vrr+WPY*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_S_S_M1_vrr = PAZ*I_ERI_F3x_S_S_S_M1_vrr+WPZ*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_S_S_M1_vrr = PAY*I_ERI_F2xy_S_S_S_M1_vrr+WPY*I_ERI_F2xy_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_S_S_M1_vrr = PAZ*I_ERI_F2xz_S_S_S_M1_vrr+WPZ*I_ERI_F2xz_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_S_S_M1_vrr = PAX*I_ERI_F3y_S_S_S_M1_vrr+WPX*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_S_S_M1_vrr = PAX*I_ERI_F3z_S_S_S_M1_vrr+WPX*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_S_S_M1_vrr = PAY*I_ERI_F3y_S_S_S_M1_vrr+WPY*I_ERI_F3y_S_S_S_M2_vrr+3*oned2z*I_ERI_D2y_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_S_S_M1_vrr = PAZ*I_ERI_F3y_S_S_S_M1_vrr+WPZ*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_S_S_M1_vrr = PAZ*I_ERI_F2yz_S_S_S_M1_vrr+WPZ*I_ERI_F2yz_S_S_S_M2_vrr+oned2z*I_ERI_D2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_S_S_M1_vrr = PAY*I_ERI_F3z_S_S_S_M1_vrr+WPY*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_S_S_M1_vrr = PAZ*I_ERI_F3z_S_S_S_M1_vrr+WPZ*I_ERI_F3z_S_S_S_M2_vrr+3*oned2z*I_ERI_D2z_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_Px_S_vrr = QCX*I_ERI_S_S_S_S_vrr+WQX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Py_S_vrr = QCY*I_ERI_S_S_S_S_vrr+WQY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Pz_S_vrr = QCZ*I_ERI_S_S_S_S_vrr+WQZ*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_Px_S_S_vrr = PBX*I_ERI_S_S_S_S_vrr+WPX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Py_S_S_vrr = PBY*I_ERI_S_S_S_S_vrr+WPY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Pz_S_S_vrr = PBZ*I_ERI_S_S_S_S_vrr+WPZ*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_Px_S_S_S_vrr = PAX*I_ERI_S_S_S_S_vrr+WPX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Py_S_S_S_vrr = PAY*I_ERI_S_S_S_S_vrr+WPY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_S_S_vrr = PAZ*I_ERI_S_S_S_S_vrr+WPZ*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_Px_S_Px_S_vrr = PAX*I_ERI_S_S_Px_S_vrr+WPX*I_ERI_S_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Py_S_Px_S_vrr = PAY*I_ERI_S_S_Px_S_vrr+WPY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Px_S_vrr = PAZ*I_ERI_S_S_Px_S_vrr+WPZ*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_Py_S_vrr = PAX*I_ERI_S_S_Py_S_vrr+WPX*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_Py_S_vrr = PAY*I_ERI_S_S_Py_S_vrr+WPY*I_ERI_S_S_Py_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_Py_S_vrr = PAZ*I_ERI_S_S_Py_S_vrr+WPZ*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_Pz_S_vrr = PAX*I_ERI_S_S_Pz_S_vrr+WPX*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Py_S_Pz_S_vrr = PAY*I_ERI_S_S_Pz_S_vrr+WPY*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Pz_S_Pz_S_vrr = PAZ*I_ERI_S_S_Pz_S_vrr+WPZ*I_ERI_S_S_Pz_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_S_vrr = PAY*I_ERI_Px_S_S_S_vrr+WPY*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_S_S_vrr = PAZ*I_ERI_Px_S_S_S_vrr+WPZ*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_S_S_vrr = PAZ*I_ERI_Py_S_S_S_vrr+WPZ*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_S_vrr = PAZ*I_ERI_Pz_S_S_S_vrr+WPZ*I_ERI_Pz_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_vrr = QCX*I_ERI_D2x_S_S_S_vrr+WQX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Px_S_vrr = QCX*I_ERI_Dxy_S_S_S_vrr+WQX*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Px_S_vrr = QCX*I_ERI_Dxz_S_S_S_vrr+WQX*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Px_S_vrr = QCX*I_ERI_D2y_S_S_S_vrr+WQX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Px_S_vrr = QCX*I_ERI_Dyz_S_S_S_vrr+WQX*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Px_S_vrr = QCX*I_ERI_D2z_S_S_S_vrr+WQX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Py_S_vrr = QCY*I_ERI_D2x_S_S_S_vrr+WQY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Py_S_vrr = QCY*I_ERI_Dxy_S_S_S_vrr+WQY*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Py_S_vrr = QCY*I_ERI_Dxz_S_S_S_vrr+WQY*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Py_S_vrr = QCY*I_ERI_D2y_S_S_S_vrr+WQY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Py_S_vrr = QCY*I_ERI_Dyz_S_S_S_vrr+WQY*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Py_S_vrr = QCY*I_ERI_D2z_S_S_S_vrr+WQY*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Pz_S_vrr = QCZ*I_ERI_D2x_S_S_S_vrr+WQZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Pz_S_vrr = QCZ*I_ERI_Dxy_S_S_S_vrr+WQZ*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Pz_S_vrr = QCZ*I_ERI_Dxz_S_S_S_vrr+WQZ*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Pz_S_vrr = QCZ*I_ERI_D2y_S_S_S_vrr+WQZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Pz_S_vrr = QCZ*I_ERI_Dyz_S_S_S_vrr+WQZ*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Pz_S_vrr = QCZ*I_ERI_D2z_S_S_S_vrr+WQZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_S_S_vrr = PAX*I_ERI_D2x_S_S_S_vrr+WPX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2z*I_ERI_Px_S_S_S_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_S_vrr = PAY*I_ERI_D2x_S_S_S_vrr+WPY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_S_vrr = PAZ*I_ERI_D2x_S_S_S_vrr+WPZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_S_vrr = PAX*I_ERI_D2y_S_S_S_vrr+WPX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_S_S_vrr = PAZ*I_ERI_Dxy_S_S_S_vrr+WPZ*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_S_vrr = PAX*I_ERI_D2z_S_S_S_vrr+WPX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_S_vrr = PAY*I_ERI_D2y_S_S_S_vrr+WPY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2z*I_ERI_Py_S_S_S_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_S_vrr = PAZ*I_ERI_D2y_S_S_S_vrr+WPZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_S_S_vrr = PAY*I_ERI_D2z_S_S_S_vrr+WPY*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_S_vrr = PAZ*I_ERI_D2z_S_S_S_vrr+WPZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2z*I_ERI_Pz_S_S_S_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_vrr = QCX*I_ERI_D2x_S_Px_S_vrr+WQX*I_ERI_D2x_S_Px_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_D2x_S_vrr = QCX*I_ERI_Dxy_S_Px_S_vrr+WQX*I_ERI_Dxy_S_Px_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_D2x_S_vrr = QCX*I_ERI_Dxz_S_Px_S_vrr+WQX*I_ERI_Dxz_S_Px_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_D2x_S_vrr = QCX*I_ERI_D2y_S_Px_S_vrr+WQX*I_ERI_D2y_S_Px_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_D2x_S_vrr = QCX*I_ERI_Dyz_S_Px_S_vrr+WQX*I_ERI_Dyz_S_Px_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_D2x_S_vrr = QCX*I_ERI_D2z_S_Px_S_vrr+WQX*I_ERI_D2z_S_Px_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dxy_S_vrr = QCY*I_ERI_D2x_S_Px_S_vrr+WQY*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxy_S_vrr = QCY*I_ERI_Dxy_S_Px_S_vrr+WQY*I_ERI_Dxy_S_Px_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxy_S_vrr = QCY*I_ERI_Dxz_S_Px_S_vrr+WQY*I_ERI_Dxz_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxy_S_vrr = QCY*I_ERI_D2y_S_Px_S_vrr+WQY*I_ERI_D2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxy_S_vrr = QCY*I_ERI_Dyz_S_Px_S_vrr+WQY*I_ERI_Dyz_S_Px_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxy_S_vrr = QCY*I_ERI_D2z_S_Px_S_vrr+WQY*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_Dxz_S_vrr = QCZ*I_ERI_D2x_S_Px_S_vrr+WQZ*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxz_S_vrr = QCZ*I_ERI_Dxy_S_Px_S_vrr+WQZ*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxz_S_vrr = QCZ*I_ERI_Dxz_S_Px_S_vrr+WQZ*I_ERI_Dxz_S_Px_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxz_S_vrr = QCZ*I_ERI_D2y_S_Px_S_vrr+WQZ*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxz_S_vrr = QCZ*I_ERI_Dyz_S_Px_S_vrr+WQZ*I_ERI_Dyz_S_Px_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxz_S_vrr = QCZ*I_ERI_D2z_S_Px_S_vrr+WQZ*I_ERI_D2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_D2y_S_vrr = QCY*I_ERI_D2x_S_Py_S_vrr+WQY*I_ERI_D2x_S_Py_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_D2y_S_vrr = QCY*I_ERI_Dxy_S_Py_S_vrr+WQY*I_ERI_Dxy_S_Py_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_D2y_S_vrr = QCY*I_ERI_Dxz_S_Py_S_vrr+WQY*I_ERI_Dxz_S_Py_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = QCY*I_ERI_D2y_S_Py_S_vrr+WQY*I_ERI_D2y_S_Py_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_D2y_S_vrr = QCY*I_ERI_Dyz_S_Py_S_vrr+WQY*I_ERI_Dyz_S_Py_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = QCY*I_ERI_D2z_S_Py_S_vrr+WQY*I_ERI_D2z_S_Py_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dyz_S_vrr = QCZ*I_ERI_D2x_S_Py_S_vrr+WQZ*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_Dyz_S_vrr = QCZ*I_ERI_Dxy_S_Py_S_vrr+WQZ*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_Dyz_S_vrr = QCZ*I_ERI_Dxz_S_Py_S_vrr+WQZ*I_ERI_Dxz_S_Py_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Dyz_S_vrr = QCZ*I_ERI_D2y_S_Py_S_vrr+WQZ*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_Dyz_S_vrr = QCZ*I_ERI_Dyz_S_Py_S_vrr+WQZ*I_ERI_Dyz_S_Py_S_M1_vrr+oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Dyz_S_vrr = QCZ*I_ERI_D2z_S_Py_S_vrr+WQZ*I_ERI_D2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = QCZ*I_ERI_D2x_S_Pz_S_vrr+WQZ*I_ERI_D2x_S_Pz_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_D2z_S_vrr = QCZ*I_ERI_Dxy_S_Pz_S_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M1_vrr+oned2e*I_ERI_Dxy_S_S_S_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_D2z_S_vrr = QCZ*I_ERI_Dxz_S_Pz_S_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M1_vrr+oned2e*I_ERI_Dxz_S_S_S_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = QCZ*I_ERI_D2y_S_Pz_S_vrr+WQZ*I_ERI_D2y_S_Pz_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_D2z_S_vrr = QCZ*I_ERI_Dyz_S_Pz_S_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M1_vrr+oned2e*I_ERI_Dyz_S_S_S_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = QCZ*I_ERI_D2z_S_Pz_S_vrr+WQZ*I_ERI_D2z_S_Pz_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_vrr = PAX*I_ERI_D2x_S_Px_S_vrr+WPX*I_ERI_D2x_S_Px_S_M1_vrr+2*oned2z*I_ERI_Px_S_Px_S_vrr-2*rhod2zsq*I_ERI_Px_S_Px_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_Px_S_vrr = PAY*I_ERI_D2x_S_Px_S_vrr+WPY*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Px_S_vrr = PAZ*I_ERI_D2x_S_Px_S_vrr+WPZ*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Px_S_vrr = PAX*I_ERI_D2y_S_Px_S_vrr+WPX*I_ERI_D2y_S_Px_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_Px_S_vrr = PAZ*I_ERI_Dxy_S_Px_S_vrr+WPZ*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Px_S_vrr = PAX*I_ERI_D2z_S_Px_S_vrr+WPX*I_ERI_D2z_S_Px_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_Px_S_vrr = PAY*I_ERI_D2y_S_Px_S_vrr+WPY*I_ERI_D2y_S_Px_S_M1_vrr+2*oned2z*I_ERI_Py_S_Px_S_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Px_S_vrr = PAZ*I_ERI_D2y_S_Px_S_vrr+WPZ*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Px_S_vrr = PAY*I_ERI_D2z_S_Px_S_vrr+WPY*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Px_S_vrr = PAZ*I_ERI_D2z_S_Px_S_vrr+WPZ*I_ERI_D2z_S_Px_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Px_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_Py_S_vrr = PAX*I_ERI_D2x_S_Py_S_vrr+WPX*I_ERI_D2x_S_Py_S_M1_vrr+2*oned2z*I_ERI_Px_S_Py_S_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_Py_S_vrr = PAY*I_ERI_D2x_S_Py_S_vrr+WPY*I_ERI_D2x_S_Py_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_Py_S_vrr = PAZ*I_ERI_D2x_S_Py_S_vrr+WPZ*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_Py_S_vrr = PAX*I_ERI_D2y_S_Py_S_vrr+WPX*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_Py_S_vrr = PAZ*I_ERI_Dxy_S_Py_S_vrr+WPZ*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_Py_S_vrr = PAX*I_ERI_D2z_S_Py_S_vrr+WPX*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_Py_S_vrr = PAY*I_ERI_D2y_S_Py_S_vrr+WPY*I_ERI_D2y_S_Py_S_M1_vrr+2*oned2z*I_ERI_Py_S_Py_S_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_Py_S_vrr = PAZ*I_ERI_D2y_S_Py_S_vrr+WPZ*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_Py_S_vrr = PAY*I_ERI_D2z_S_Py_S_vrr+WPY*I_ERI_D2z_S_Py_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_Py_S_vrr = PAZ*I_ERI_D2z_S_Py_S_vrr+WPZ*I_ERI_D2z_S_Py_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Py_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_Pz_S_vrr = PAX*I_ERI_D2x_S_Pz_S_vrr+WPX*I_ERI_D2x_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_F2xy_S_Pz_S_vrr = PAY*I_ERI_D2x_S_Pz_S_vrr+WPY*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xz_S_Pz_S_vrr = PAZ*I_ERI_D2x_S_Pz_S_vrr+WPZ*I_ERI_D2x_S_Pz_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_Pz_S_vrr = PAX*I_ERI_D2y_S_Pz_S_vrr+WPX*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Pz_S_vrr = PAZ*I_ERI_Dxy_S_Pz_S_vrr+WPZ*I_ERI_Dxy_S_Pz_S_M1_vrr+oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_Pz_S_vrr = PAX*I_ERI_D2z_S_Pz_S_vrr+WPX*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3y_S_Pz_S_vrr = PAY*I_ERI_D2y_S_Pz_S_vrr+WPY*I_ERI_D2y_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_F2yz_S_Pz_S_vrr = PAZ*I_ERI_D2y_S_Pz_S_vrr+WPZ*I_ERI_D2y_S_Pz_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_Pz_S_vrr = PAY*I_ERI_D2z_S_Pz_S_vrr+WPY*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3z_S_Pz_S_vrr = PAZ*I_ERI_D2z_S_Pz_S_vrr+WPZ*I_ERI_D2z_S_Pz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_S_S_vrr = PAX*I_ERI_F3x_S_S_S_vrr+WPX*I_ERI_F3x_S_S_S_M1_vrr+3*oned2z*I_ERI_D2x_S_S_S_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_S_S_vrr = PAY*I_ERI_F3x_S_S_S_vrr+WPY*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_S_S_vrr = PAZ*I_ERI_F3x_S_S_S_vrr+WPZ*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_S_S_vrr = PAY*I_ERI_F2xy_S_S_S_vrr+WPY*I_ERI_F2xy_S_S_S_M1_vrr+oned2z*I_ERI_D2x_S_S_S_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_S_S_vrr = PAZ*I_ERI_F2xy_S_S_S_vrr+WPZ*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_S_S_vrr = PAZ*I_ERI_F2xz_S_S_S_vrr+WPZ*I_ERI_F2xz_S_S_S_M1_vrr+oned2z*I_ERI_D2x_S_S_S_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_S_S_vrr = PAX*I_ERI_F3y_S_S_S_vrr+WPX*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_S_S_vrr = PAZ*I_ERI_Fx2y_S_S_S_vrr+WPZ*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_S_S_vrr = PAY*I_ERI_Fx2z_S_S_S_vrr+WPY*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_S_S_vrr = PAX*I_ERI_F3z_S_S_S_vrr+WPX*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_S_S_vrr = PAY*I_ERI_F3y_S_S_S_vrr+WPY*I_ERI_F3y_S_S_S_M1_vrr+3*oned2z*I_ERI_D2y_S_S_S_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_S_S_vrr = PAZ*I_ERI_F3y_S_S_S_vrr+WPZ*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_S_S_vrr = PAZ*I_ERI_F2yz_S_S_S_vrr+WPZ*I_ERI_F2yz_S_S_S_M1_vrr+oned2z*I_ERI_D2y_S_S_S_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_S_S_vrr = PAY*I_ERI_F3z_S_S_S_vrr+WPY*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_S_S_vrr = PAZ*I_ERI_F3z_S_S_S_vrr+WPZ*I_ERI_F3z_S_S_S_M1_vrr+3*oned2z*I_ERI_D2z_S_S_S_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_vrr = QCX*I_ERI_F3x_S_Px_S_vrr+WQX*I_ERI_F3x_S_Px_S_M1_vrr+oned2e*I_ERI_F3x_S_S_S_vrr-rhod2esq*I_ERI_F3x_S_S_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xy_S_D2x_S_vrr = QCX*I_ERI_F2xy_S_Px_S_vrr+WQX*I_ERI_F2xy_S_Px_S_M1_vrr+oned2e*I_ERI_F2xy_S_S_S_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_D2x_S_vrr = QCX*I_ERI_F2xz_S_Px_S_vrr+WQX*I_ERI_F2xz_S_Px_S_M1_vrr+oned2e*I_ERI_F2xz_S_S_S_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2x_S_vrr = QCX*I_ERI_Fx2y_S_Px_S_vrr+WQX*I_ERI_Fx2y_S_Px_S_M1_vrr+oned2e*I_ERI_Fx2y_S_S_S_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2x_S_vrr = QCX*I_ERI_Fxyz_S_Px_S_vrr+WQX*I_ERI_Fxyz_S_Px_S_M1_vrr+oned2e*I_ERI_Fxyz_S_S_S_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dyz_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2x_S_vrr = QCX*I_ERI_Fx2z_S_Px_S_vrr+WQX*I_ERI_Fx2z_S_Px_S_M1_vrr+oned2e*I_ERI_Fx2z_S_S_S_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3y_S_D2x_S_vrr = QCX*I_ERI_F3y_S_Px_S_vrr+WQX*I_ERI_F3y_S_Px_S_M1_vrr+oned2e*I_ERI_F3y_S_S_S_vrr-rhod2esq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_D2x_S_vrr = QCX*I_ERI_F2yz_S_Px_S_vrr+WQX*I_ERI_F2yz_S_Px_S_M1_vrr+oned2e*I_ERI_F2yz_S_S_S_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2x_S_vrr = QCX*I_ERI_Fy2z_S_Px_S_vrr+WQX*I_ERI_Fy2z_S_Px_S_M1_vrr+oned2e*I_ERI_Fy2z_S_S_S_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_D2x_S_vrr = QCX*I_ERI_F3z_S_Px_S_vrr+WQX*I_ERI_F3z_S_Px_S_M1_vrr+oned2e*I_ERI_F3z_S_S_S_vrr-rhod2esq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_Dxy_S_vrr = QCY*I_ERI_F3x_S_Px_S_vrr+WQY*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxy_S_vrr = QCY*I_ERI_F2xy_S_Px_S_vrr+WQY*I_ERI_F2xy_S_Px_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxy_S_vrr = QCY*I_ERI_F2xz_S_Px_S_vrr+WQY*I_ERI_F2xz_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_vrr = QCY*I_ERI_Fx2y_S_Px_S_vrr+WQY*I_ERI_Fx2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxy_S_vrr = QCY*I_ERI_Fxyz_S_Px_S_vrr+WQY*I_ERI_Fxyz_S_Px_S_M1_vrr+oned2k*I_ERI_Dxz_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_vrr = QCY*I_ERI_Fx2z_S_Px_S_vrr+WQY*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_F3y_S_Dxy_S_vrr = QCY*I_ERI_F3y_S_Px_S_vrr+WQY*I_ERI_F3y_S_Px_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxy_S_vrr = QCY*I_ERI_F2yz_S_Px_S_vrr+WQY*I_ERI_F2yz_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxy_S_vrr = QCY*I_ERI_Fy2z_S_Px_S_vrr+WQY*I_ERI_Fy2z_S_Px_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Dxy_S_vrr = QCY*I_ERI_F3z_S_Px_S_vrr+WQY*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_Dxz_S_vrr = QCZ*I_ERI_F3x_S_Px_S_vrr+WQZ*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxz_S_vrr = QCZ*I_ERI_F2xy_S_Px_S_vrr+WQZ*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxz_S_vrr = QCZ*I_ERI_F2xz_S_Px_S_vrr+WQZ*I_ERI_F2xz_S_Px_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_vrr = QCZ*I_ERI_Fx2y_S_Px_S_vrr+WQZ*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxz_S_vrr = QCZ*I_ERI_Fxyz_S_Px_S_vrr+WQZ*I_ERI_Fxyz_S_Px_S_M1_vrr+oned2k*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_vrr = QCZ*I_ERI_Fx2z_S_Px_S_vrr+WQZ*I_ERI_Fx2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Px_S_M1_vrr;
      Double I_ERI_F3y_S_Dxz_S_vrr = QCZ*I_ERI_F3y_S_Px_S_vrr+WQZ*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxz_S_vrr = QCZ*I_ERI_F2yz_S_Px_S_vrr+WQZ*I_ERI_F2yz_S_Px_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxz_S_vrr = QCZ*I_ERI_Fy2z_S_Px_S_vrr+WQZ*I_ERI_Fy2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Dxz_S_vrr = QCZ*I_ERI_F3z_S_Px_S_vrr+WQZ*I_ERI_F3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_D2y_S_vrr = QCY*I_ERI_F3x_S_Py_S_vrr+WQY*I_ERI_F3x_S_Py_S_M1_vrr+oned2e*I_ERI_F3x_S_S_S_vrr-rhod2esq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_D2y_S_vrr = QCY*I_ERI_F2xy_S_Py_S_vrr+WQY*I_ERI_F2xy_S_Py_S_M1_vrr+oned2e*I_ERI_F2xy_S_S_S_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_F2xz_S_D2y_S_vrr = QCY*I_ERI_F2xz_S_Py_S_vrr+WQY*I_ERI_F2xz_S_Py_S_M1_vrr+oned2e*I_ERI_F2xz_S_S_S_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2y_S_vrr = QCY*I_ERI_Fx2y_S_Py_S_vrr+WQY*I_ERI_Fx2y_S_Py_S_M1_vrr+oned2e*I_ERI_Fx2y_S_S_S_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2y_S_vrr = QCY*I_ERI_Fxyz_S_Py_S_vrr+WQY*I_ERI_Fxyz_S_Py_S_M1_vrr+oned2e*I_ERI_Fxyz_S_S_S_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxz_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2y_S_vrr = QCY*I_ERI_Fx2z_S_Py_S_vrr+WQY*I_ERI_Fx2z_S_Py_S_M1_vrr+oned2e*I_ERI_Fx2z_S_S_S_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_D2y_S_vrr = QCY*I_ERI_F3y_S_Py_S_vrr+WQY*I_ERI_F3y_S_Py_S_M1_vrr+oned2e*I_ERI_F3y_S_S_S_vrr-rhod2esq*I_ERI_F3y_S_S_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_F2yz_S_D2y_S_vrr = QCY*I_ERI_F2yz_S_Py_S_vrr+WQY*I_ERI_F2yz_S_Py_S_M1_vrr+oned2e*I_ERI_F2yz_S_S_S_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2y_S_vrr = QCY*I_ERI_Fy2z_S_Py_S_vrr+WQY*I_ERI_Fy2z_S_Py_S_M1_vrr+oned2e*I_ERI_Fy2z_S_S_S_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3z_S_D2y_S_vrr = QCY*I_ERI_F3z_S_Py_S_vrr+WQY*I_ERI_F3z_S_Py_S_M1_vrr+oned2e*I_ERI_F3z_S_S_S_vrr-rhod2esq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_Dyz_S_vrr = QCZ*I_ERI_F3x_S_Py_S_vrr+WQZ*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_Dyz_S_vrr = QCZ*I_ERI_F2xy_S_Py_S_vrr+WQZ*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_F2xz_S_Dyz_S_vrr = QCZ*I_ERI_F2xz_S_Py_S_vrr+WQZ*I_ERI_F2xz_S_Py_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_vrr = QCZ*I_ERI_Fx2y_S_Py_S_vrr+WQZ*I_ERI_Fx2y_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dyz_S_vrr = QCZ*I_ERI_Fxyz_S_Py_S_vrr+WQZ*I_ERI_Fxyz_S_Py_S_M1_vrr+oned2k*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_vrr = QCZ*I_ERI_Fx2z_S_Py_S_vrr+WQZ*I_ERI_Fx2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_Dyz_S_vrr = QCZ*I_ERI_F3y_S_Py_S_vrr+WQZ*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_F2yz_S_Dyz_S_vrr = QCZ*I_ERI_F2yz_S_Py_S_vrr+WQZ*I_ERI_F2yz_S_Py_S_M1_vrr+oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dyz_S_vrr = QCZ*I_ERI_Fy2z_S_Py_S_vrr+WQZ*I_ERI_Fy2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Py_S_M1_vrr;
      Double I_ERI_F3z_S_Dyz_S_vrr = QCZ*I_ERI_F3z_S_Py_S_vrr+WQZ*I_ERI_F3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_D2z_S_vrr = QCZ*I_ERI_F3x_S_Pz_S_vrr+WQZ*I_ERI_F3x_S_Pz_S_M1_vrr+oned2e*I_ERI_F3x_S_S_S_vrr-rhod2esq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_D2z_S_vrr = QCZ*I_ERI_F2xy_S_Pz_S_vrr+WQZ*I_ERI_F2xy_S_Pz_S_M1_vrr+oned2e*I_ERI_F2xy_S_S_S_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_D2z_S_vrr = QCZ*I_ERI_F2xz_S_Pz_S_vrr+WQZ*I_ERI_F2xz_S_Pz_S_M1_vrr+oned2e*I_ERI_F2xz_S_S_S_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2z_S_vrr = QCZ*I_ERI_Fx2y_S_Pz_S_vrr+WQZ*I_ERI_Fx2y_S_Pz_S_M1_vrr+oned2e*I_ERI_Fx2y_S_S_S_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2z_S_vrr = QCZ*I_ERI_Fxyz_S_Pz_S_vrr+WQZ*I_ERI_Fxyz_S_Pz_S_M1_vrr+oned2e*I_ERI_Fxyz_S_S_S_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxy_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2z_S_vrr = QCZ*I_ERI_Fx2z_S_Pz_S_vrr+WQZ*I_ERI_Fx2z_S_Pz_S_M1_vrr+oned2e*I_ERI_Fx2z_S_S_S_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Pz_S_M1_vrr;
      Double I_ERI_F3y_S_D2z_S_vrr = QCZ*I_ERI_F3y_S_Pz_S_vrr+WQZ*I_ERI_F3y_S_Pz_S_M1_vrr+oned2e*I_ERI_F3y_S_S_S_vrr-rhod2esq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_D2z_S_vrr = QCZ*I_ERI_F2yz_S_Pz_S_vrr+WQZ*I_ERI_F2yz_S_Pz_S_M1_vrr+oned2e*I_ERI_F2yz_S_S_S_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2z_S_vrr = QCZ*I_ERI_Fy2z_S_Pz_S_vrr+WQZ*I_ERI_Fy2z_S_Pz_S_M1_vrr+oned2e*I_ERI_Fy2z_S_S_S_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Pz_S_M1_vrr;
      Double I_ERI_F3z_S_D2z_S_vrr = QCZ*I_ERI_F3z_S_Pz_S_vrr+WQZ*I_ERI_F3z_S_Pz_S_M1_vrr+oned2e*I_ERI_F3z_S_S_S_vrr-rhod2esq*I_ERI_F3z_S_S_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_vrr = PAX*I_ERI_F3x_S_Px_S_vrr+WPX*I_ERI_F3x_S_Px_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Px_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Px_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Px_S_vrr = PAY*I_ERI_F3x_S_Px_S_vrr+WPY*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G3xz_S_Px_S_vrr = PAZ*I_ERI_F3x_S_Px_S_vrr+WPZ*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G2x2y_S_Px_S_vrr = PAY*I_ERI_F2xy_S_Px_S_vrr+WPY*I_ERI_F2xy_S_Px_S_M1_vrr+oned2z*I_ERI_D2x_S_Px_S_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_G2xyz_S_Px_S_vrr = PAZ*I_ERI_F2xy_S_Px_S_vrr+WPZ*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G2x2z_S_Px_S_vrr = PAZ*I_ERI_F2xz_S_Px_S_vrr+WPZ*I_ERI_F2xz_S_Px_S_M1_vrr+oned2z*I_ERI_D2x_S_Px_S_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Gx3y_S_Px_S_vrr = PAX*I_ERI_F3y_S_Px_S_vrr+WPX*I_ERI_F3y_S_Px_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Px_S_vrr = PAZ*I_ERI_Fx2y_S_Px_S_vrr+WPZ*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Px_S_vrr = PAY*I_ERI_Fx2z_S_Px_S_vrr+WPY*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3z_S_Px_S_vrr = PAX*I_ERI_F3z_S_Px_S_vrr+WPX*I_ERI_F3z_S_Px_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Px_S_vrr = PAY*I_ERI_F3y_S_Px_S_vrr+WPY*I_ERI_F3y_S_Px_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Px_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_Px_S_vrr = PAZ*I_ERI_F3y_S_Px_S_vrr+WPZ*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_Px_S_vrr = PAZ*I_ERI_F2yz_S_Px_S_vrr+WPZ*I_ERI_F2yz_S_Px_S_M1_vrr+oned2z*I_ERI_D2y_S_Px_S_vrr-rhod2zsq*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_Px_S_vrr = PAY*I_ERI_F3z_S_Px_S_vrr+WPY*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_Px_S_vrr = PAZ*I_ERI_F3z_S_Px_S_vrr+WPZ*I_ERI_F3z_S_Px_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Px_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_Py_S_vrr = PAX*I_ERI_F3x_S_Py_S_vrr+WPX*I_ERI_F3x_S_Py_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Py_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_G3xy_S_Py_S_vrr = PAY*I_ERI_F3x_S_Py_S_vrr+WPY*I_ERI_F3x_S_Py_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Py_S_vrr = PAZ*I_ERI_F3x_S_Py_S_vrr+WPZ*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G2x2y_S_Py_S_vrr = PAY*I_ERI_F2xy_S_Py_S_vrr+WPY*I_ERI_F2xy_S_Py_S_M1_vrr+oned2z*I_ERI_D2x_S_Py_S_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M1_vrr+oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Py_S_vrr = PAZ*I_ERI_F2xy_S_Py_S_vrr+WPZ*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_G2x2z_S_Py_S_vrr = PAZ*I_ERI_F2xz_S_Py_S_vrr+WPZ*I_ERI_F2xz_S_Py_S_M1_vrr+oned2z*I_ERI_D2x_S_Py_S_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Gx3y_S_Py_S_vrr = PAX*I_ERI_F3y_S_Py_S_vrr+WPX*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Py_S_vrr = PAZ*I_ERI_Fx2y_S_Py_S_vrr+WPZ*I_ERI_Fx2y_S_Py_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Py_S_vrr = PAY*I_ERI_Fx2z_S_Py_S_vrr+WPY*I_ERI_Fx2z_S_Py_S_M1_vrr+oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Py_S_vrr = PAX*I_ERI_F3z_S_Py_S_vrr+WPX*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4y_S_Py_S_vrr = PAY*I_ERI_F3y_S_Py_S_vrr+WPY*I_ERI_F3y_S_Py_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Py_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Py_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Py_S_vrr = PAZ*I_ERI_F3y_S_Py_S_vrr+WPZ*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_G2y2z_S_Py_S_vrr = PAZ*I_ERI_F2yz_S_Py_S_vrr+WPZ*I_ERI_F2yz_S_Py_S_M1_vrr+oned2z*I_ERI_D2y_S_Py_S_vrr-rhod2zsq*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Gy3z_S_Py_S_vrr = PAY*I_ERI_F3z_S_Py_S_vrr+WPY*I_ERI_F3z_S_Py_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Py_S_vrr = PAZ*I_ERI_F3z_S_Py_S_vrr+WPZ*I_ERI_F3z_S_Py_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Py_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_G4x_S_Pz_S_vrr = PAX*I_ERI_F3x_S_Pz_S_vrr+WPX*I_ERI_F3x_S_Pz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Pz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xy_S_Pz_S_vrr = PAY*I_ERI_F3x_S_Pz_S_vrr+WPY*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xz_S_Pz_S_vrr = PAZ*I_ERI_F3x_S_Pz_S_vrr+WPZ*I_ERI_F3x_S_Pz_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Pz_S_vrr = PAY*I_ERI_F2xy_S_Pz_S_vrr+WPY*I_ERI_F2xy_S_Pz_S_M1_vrr+oned2z*I_ERI_D2x_S_Pz_S_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Pz_S_vrr = PAZ*I_ERI_F2xy_S_Pz_S_vrr+WPZ*I_ERI_F2xy_S_Pz_S_M1_vrr+oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Pz_S_vrr = PAZ*I_ERI_F2xz_S_Pz_S_vrr+WPZ*I_ERI_F2xz_S_Pz_S_M1_vrr+oned2z*I_ERI_D2x_S_Pz_S_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M1_vrr+oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Pz_S_vrr = PAX*I_ERI_F3y_S_Pz_S_vrr+WPX*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_vrr = PAZ*I_ERI_Fx2y_S_Pz_S_vrr+WPZ*I_ERI_Fx2y_S_Pz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_vrr = PAY*I_ERI_Fx2z_S_Pz_S_vrr+WPY*I_ERI_Fx2z_S_Pz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Pz_S_vrr = PAX*I_ERI_F3z_S_Pz_S_vrr+WPX*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_G4y_S_Pz_S_vrr = PAY*I_ERI_F3y_S_Pz_S_vrr+WPY*I_ERI_F3y_S_Pz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Pz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_G3yz_S_Pz_S_vrr = PAZ*I_ERI_F3y_S_Pz_S_vrr+WPZ*I_ERI_F3y_S_Pz_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Pz_S_vrr = PAZ*I_ERI_F2yz_S_Pz_S_vrr+WPZ*I_ERI_F2yz_S_Pz_S_M1_vrr+oned2z*I_ERI_D2y_S_Pz_S_vrr-rhod2zsq*I_ERI_D2y_S_Pz_S_M1_vrr+oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Pz_S_vrr = PAY*I_ERI_F3z_S_Pz_S_vrr+WPY*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_G4z_S_Pz_S_vrr = PAZ*I_ERI_F3z_S_Pz_S_vrr+WPZ*I_ERI_F3z_S_Pz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Pz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Pz_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_S_S_vrr = PAX*I_ERI_G4x_S_S_S_vrr+WPX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2z*I_ERI_F3x_S_S_S_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_S_S_vrr = PAY*I_ERI_G4x_S_S_S_vrr+WPY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_S_S_vrr = PAZ*I_ERI_G4x_S_S_S_vrr+WPZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_S_S_vrr = PAY*I_ERI_G3xy_S_S_S_vrr+WPY*I_ERI_G3xy_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_S_S_vrr = PAZ*I_ERI_G3xy_S_S_S_vrr+WPZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_S_S_vrr = PAZ*I_ERI_G3xz_S_S_S_vrr+WPZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_S_S_vrr = PAX*I_ERI_Gx3y_S_S_S_vrr+WPX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_S_S_vrr = PAZ*I_ERI_G2x2y_S_S_S_vrr+WPZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_S_S_vrr = PAY*I_ERI_G2x2z_S_S_S_vrr+WPY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_S_S_vrr = PAX*I_ERI_Gx3z_S_S_S_vrr+WPX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_S_S_vrr = PAX*I_ERI_G4y_S_S_S_vrr+WPX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_S_S_vrr = PAZ*I_ERI_Gx3y_S_S_S_vrr+WPZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_S_S_vrr = PAX*I_ERI_G2y2z_S_S_S_vrr+WPX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_S_S_vrr = PAY*I_ERI_Gx3z_S_S_S_vrr+WPY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_S_S_vrr = PAX*I_ERI_G4z_S_S_S_vrr+WPX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_S_S_vrr = PAY*I_ERI_G4y_S_S_S_vrr+WPY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2z*I_ERI_F3y_S_S_S_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_S_S_vrr = PAZ*I_ERI_G4y_S_S_S_vrr+WPZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_S_S_vrr = PAZ*I_ERI_G3yz_S_S_S_vrr+WPZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_S_S_vrr = PAY*I_ERI_Gy3z_S_S_S_vrr+WPY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_S_S_vrr = PAY*I_ERI_G4z_S_S_S_vrr+WPY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_S_S_vrr = PAZ*I_ERI_G4z_S_S_S_vrr+WPZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2z*I_ERI_F3z_S_S_S_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C2_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C2_aa_coefs = ic2*jc2*alpha*alpha;
      I_ERI_G4x_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C2_aa += SQ_ERI_G_S_S_S_C2_aa_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C2_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C2_a_coefs = ic2*jc2*alpha;
      I_ERI_D2x_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C2_a += SQ_ERI_D_S_S_S_C2_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C2
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C2_coefs = ic2*jc2;
      I_ERI_S_S_S_S_C2 += SQ_ERI_S_S_S_S_C2_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1002_coefs = ic2_1*jc2;
      I_ERI_S_Px_S_S_C1002 += SQ_ERI_S_P_S_S_C1002_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1002 += SQ_ERI_S_P_S_S_C1002_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1002 += SQ_ERI_S_P_S_S_C1002_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1002_a_coefs = ic2_1*jc2*alpha;
      I_ERI_F3x_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1002_a += SQ_ERI_F_S_S_S_C1002_a_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1002_coefs = ic2_1*jc2;
      I_ERI_Px_S_S_S_C1002 += SQ_ERI_P_S_S_S_C1002_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1002 += SQ_ERI_P_S_S_S_C1002_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1002 += SQ_ERI_P_S_S_S_C1002_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C2_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C2_ac_coefs = ic2*jc2*alpha*gamma;
      I_ERI_F3x_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C2_ac += SQ_ERI_F_S_P_S_C2_ac_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C2_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C2_c_coefs = ic2*jc2*gamma;
      I_ERI_Px_S_Px_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C2_c += SQ_ERI_P_S_P_S_C2_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C2_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C2_b_coefs = ic2*jc2*beta;
      I_ERI_D2x_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C2_b += SQ_ERI_D_S_S_S_C2_b_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1002_c_coefs = ic2_1*jc2*gamma;
      I_ERI_D2x_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1002_c += SQ_ERI_D_S_P_S_C1002_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C2_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C2_cc_coefs = ic2*jc2*gamma*gamma;
      I_ERI_D2x_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C2_cc += SQ_ERI_D_S_D_S_C2_cc_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C2_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C2_c_coefs = ic2*jc2*gamma;
      I_ERI_D2x_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C2_c += SQ_ERI_D_S_S_S_C2_c_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_C1002_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_C1002_aa_coefs = ic2_1*jc2*alpha*alpha;
      I_ERI_H5x_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_C1002_aa += SQ_ERI_H_S_S_S_C1002_aa_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1002_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1002_aa_coefs = ic2_1*jc2*alpha*alpha;
      I_ERI_G4x_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1002_aa += SQ_ERI_G_S_S_S_C1002_aa_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1002_a_coefs = ic2_1*jc2*alpha;
      I_ERI_D2x_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1002_a += SQ_ERI_D_S_S_S_C1002_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C2_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C2_ab_coefs = ic2*jc2*alpha*beta;
      I_ERI_G4x_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C2_ab += SQ_ERI_G_S_S_S_C2_ab_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C2_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C2_ab_coefs = ic2*jc2*alpha*beta;
      I_ERI_F3x_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C2_ab += SQ_ERI_F_S_S_S_C2_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C2_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C2_b_coefs = ic2*jc2*beta;
      I_ERI_Px_S_S_S_C2_b += SQ_ERI_P_S_S_S_C2_b_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C2_b += SQ_ERI_P_S_S_S_C2_b_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C2_b += SQ_ERI_P_S_S_S_C2_b_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1002_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1002_ac_coefs = ic2_1*jc2*alpha*gamma;
      I_ERI_G4x_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1002_ac += SQ_ERI_G_S_P_S_C1002_ac_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1002_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1002_ac_coefs = ic2_1*jc2*alpha*gamma;
      I_ERI_F3x_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1002_ac += SQ_ERI_F_S_P_S_C1002_ac_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1002_c_coefs = ic2_1*jc2*gamma;
      I_ERI_Px_S_Px_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1002_c += SQ_ERI_P_S_P_S_C1002_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1002_b_coefs = ic2_1*jc2*beta;
      I_ERI_F3x_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1002_b += SQ_ERI_F_S_S_S_C1002_b_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1002_b_coefs = ic2_1*jc2*beta;
      I_ERI_D2x_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1002_b += SQ_ERI_D_S_S_S_C1002_b_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C2_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C2_bc_coefs = ic2*jc2*beta*gamma;
      I_ERI_F3x_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C2_bc += SQ_ERI_F_S_P_S_C2_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C2_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C2_bc_coefs = ic2*jc2*beta*gamma;
      I_ERI_D2x_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C2_bc += SQ_ERI_D_S_P_S_C2_bc_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1002_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1002_cc_coefs = ic2_1*jc2*gamma*gamma;
      I_ERI_F3x_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1002_cc += SQ_ERI_F_S_D_S_C1002_cc_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1002_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1002_cc_coefs = ic2_1*jc2*gamma*gamma;
      I_ERI_D2x_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1002_cc += SQ_ERI_D_S_D_S_C1002_cc_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1002_c_coefs = ic2_1*jc2*gamma;
      I_ERI_F3x_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1002_c += SQ_ERI_F_S_S_S_C1002_c_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1002_c_coefs = ic2_1*jc2*gamma;
      I_ERI_D2x_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1002_c += SQ_ERI_D_S_S_S_C1002_c_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_C1002_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_C1002_ab_coefs = ic2_1*jc2*alpha*beta;
      I_ERI_H5x_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_C1002_ab += SQ_ERI_H_S_S_S_C1002_ab_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1002_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1002_ab_coefs = ic2_1*jc2*alpha*beta;
      I_ERI_G4x_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1002_ab += SQ_ERI_G_S_S_S_C1002_ab_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1002_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1002_ab_coefs = ic2_1*jc2*alpha*beta;
      I_ERI_F3x_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1002_ab += SQ_ERI_F_S_S_S_C1002_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1002_b_coefs = ic2_1*jc2*beta;
      I_ERI_Px_S_S_S_C1002_b += SQ_ERI_P_S_S_S_C1002_b_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1002_b += SQ_ERI_P_S_S_S_C1002_b_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1002_b += SQ_ERI_P_S_S_S_C1002_b_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C2_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C2_bb_coefs = ic2*jc2*beta*beta;
      I_ERI_G4x_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C2_bb += SQ_ERI_G_S_S_S_C2_bb_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C2_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C2_bb_coefs = ic2*jc2*beta*beta;
      I_ERI_F3x_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C2_bb += SQ_ERI_F_S_S_S_C2_bb_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C2_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C2_bb_coefs = ic2*jc2*beta*beta;
      I_ERI_D2x_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C2_bb += SQ_ERI_D_S_S_S_C2_bb_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1002_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1002_bc_coefs = ic2_1*jc2*beta*gamma;
      I_ERI_G4x_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1002_bc += SQ_ERI_G_S_P_S_C1002_bc_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1002_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1002_bc_coefs = ic2_1*jc2*beta*gamma;
      I_ERI_F3x_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1002_bc += SQ_ERI_F_S_P_S_C1002_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1002_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1002_bc_coefs = ic2_1*jc2*beta*gamma;
      I_ERI_D2x_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1002_bc += SQ_ERI_D_S_P_S_C1002_bc_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_C1002_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_C1002_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_H5x_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_C1002_bb += SQ_ERI_H_S_S_S_C1002_bb_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1002_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1002_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_G4x_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1002_bb += SQ_ERI_G_S_S_S_C1002_bb_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1002_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1002_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_F3x_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1002_bb += SQ_ERI_F_S_S_S_C1002_bb_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1002_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1002_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_D2x_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1002_bb += SQ_ERI_D_S_S_S_C1002_bb_coefs*I_ERI_D2z_S_S_S_vrr;
    }
  }

  /************************************************************
   * let's see the significance test result. if VRR result is
   * insignificant, there's no need to do following codes
   ************************************************************/
  if (! isSignificant) return;

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
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002_a
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1002_a = I_ERI_F3x_S_S_S_C1002_a+ABX*I_ERI_D2x_S_S_S_C1002_a;
  Double I_ERI_Dxy_Px_S_S_C1002_a = I_ERI_F2xy_S_S_S_C1002_a+ABX*I_ERI_Dxy_S_S_S_C1002_a;
  Double I_ERI_Dxz_Px_S_S_C1002_a = I_ERI_F2xz_S_S_S_C1002_a+ABX*I_ERI_Dxz_S_S_S_C1002_a;
  Double I_ERI_D2y_Px_S_S_C1002_a = I_ERI_Fx2y_S_S_S_C1002_a+ABX*I_ERI_D2y_S_S_S_C1002_a;
  Double I_ERI_Dyz_Px_S_S_C1002_a = I_ERI_Fxyz_S_S_S_C1002_a+ABX*I_ERI_Dyz_S_S_S_C1002_a;
  Double I_ERI_D2z_Px_S_S_C1002_a = I_ERI_Fx2z_S_S_S_C1002_a+ABX*I_ERI_D2z_S_S_S_C1002_a;
  Double I_ERI_D2x_Py_S_S_C1002_a = I_ERI_F2xy_S_S_S_C1002_a+ABY*I_ERI_D2x_S_S_S_C1002_a;
  Double I_ERI_Dxy_Py_S_S_C1002_a = I_ERI_Fx2y_S_S_S_C1002_a+ABY*I_ERI_Dxy_S_S_S_C1002_a;
  Double I_ERI_Dxz_Py_S_S_C1002_a = I_ERI_Fxyz_S_S_S_C1002_a+ABY*I_ERI_Dxz_S_S_S_C1002_a;
  Double I_ERI_D2y_Py_S_S_C1002_a = I_ERI_F3y_S_S_S_C1002_a+ABY*I_ERI_D2y_S_S_S_C1002_a;
  Double I_ERI_Dyz_Py_S_S_C1002_a = I_ERI_F2yz_S_S_S_C1002_a+ABY*I_ERI_Dyz_S_S_S_C1002_a;
  Double I_ERI_D2z_Py_S_S_C1002_a = I_ERI_Fy2z_S_S_S_C1002_a+ABY*I_ERI_D2z_S_S_S_C1002_a;
  Double I_ERI_D2x_Pz_S_S_C1002_a = I_ERI_F2xz_S_S_S_C1002_a+ABZ*I_ERI_D2x_S_S_S_C1002_a;
  Double I_ERI_Dxy_Pz_S_S_C1002_a = I_ERI_Fxyz_S_S_S_C1002_a+ABZ*I_ERI_Dxy_S_S_S_C1002_a;
  Double I_ERI_Dxz_Pz_S_S_C1002_a = I_ERI_Fx2z_S_S_S_C1002_a+ABZ*I_ERI_Dxz_S_S_S_C1002_a;
  Double I_ERI_D2y_Pz_S_S_C1002_a = I_ERI_F2yz_S_S_S_C1002_a+ABZ*I_ERI_D2y_S_S_S_C1002_a;
  Double I_ERI_Dyz_Pz_S_S_C1002_a = I_ERI_Fy2z_S_S_S_C1002_a+ABZ*I_ERI_Dyz_S_S_S_C1002_a;
  Double I_ERI_D2z_Pz_S_S_C1002_a = I_ERI_F3z_S_S_S_C1002_a+ABZ*I_ERI_D2z_S_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C2_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C2_b
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C2_b = I_ERI_D2x_S_S_S_C2_b+ABX*I_ERI_Px_S_S_S_C2_b;
  Double I_ERI_Py_Px_S_S_C2_b = I_ERI_Dxy_S_S_S_C2_b+ABX*I_ERI_Py_S_S_S_C2_b;
  Double I_ERI_Pz_Px_S_S_C2_b = I_ERI_Dxz_S_S_S_C2_b+ABX*I_ERI_Pz_S_S_S_C2_b;
  Double I_ERI_Px_Py_S_S_C2_b = I_ERI_Dxy_S_S_S_C2_b+ABY*I_ERI_Px_S_S_S_C2_b;
  Double I_ERI_Py_Py_S_S_C2_b = I_ERI_D2y_S_S_S_C2_b+ABY*I_ERI_Py_S_S_S_C2_b;
  Double I_ERI_Pz_Py_S_S_C2_b = I_ERI_Dyz_S_S_S_C2_b+ABY*I_ERI_Pz_S_S_S_C2_b;
  Double I_ERI_Px_Pz_S_S_C2_b = I_ERI_Dxz_S_S_S_C2_b+ABZ*I_ERI_Px_S_S_S_C2_b;
  Double I_ERI_Py_Pz_S_S_C2_b = I_ERI_Dyz_S_S_S_C2_b+ABZ*I_ERI_Py_S_S_S_C2_b;
  Double I_ERI_Pz_Pz_S_S_C2_b = I_ERI_D2z_S_S_S_C2_b+ABZ*I_ERI_Pz_S_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002_b
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1002_b = I_ERI_D2x_S_S_S_C1002_b+ABX*I_ERI_Px_S_S_S_C1002_b;
  Double I_ERI_Py_Px_S_S_C1002_b = I_ERI_Dxy_S_S_S_C1002_b+ABX*I_ERI_Py_S_S_S_C1002_b;
  Double I_ERI_Pz_Px_S_S_C1002_b = I_ERI_Dxz_S_S_S_C1002_b+ABX*I_ERI_Pz_S_S_S_C1002_b;
  Double I_ERI_Px_Py_S_S_C1002_b = I_ERI_Dxy_S_S_S_C1002_b+ABY*I_ERI_Px_S_S_S_C1002_b;
  Double I_ERI_Py_Py_S_S_C1002_b = I_ERI_D2y_S_S_S_C1002_b+ABY*I_ERI_Py_S_S_S_C1002_b;
  Double I_ERI_Pz_Py_S_S_C1002_b = I_ERI_Dyz_S_S_S_C1002_b+ABY*I_ERI_Pz_S_S_S_C1002_b;
  Double I_ERI_Px_Pz_S_S_C1002_b = I_ERI_Dxz_S_S_S_C1002_b+ABZ*I_ERI_Px_S_S_S_C1002_b;
  Double I_ERI_Py_Pz_S_S_C1002_b = I_ERI_Dyz_S_S_S_C1002_b+ABZ*I_ERI_Py_S_S_S_C1002_b;
  Double I_ERI_Pz_Pz_S_S_C1002_b = I_ERI_D2z_S_S_S_C1002_b+ABZ*I_ERI_Pz_S_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002_b
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1002_b = I_ERI_F3x_S_S_S_C1002_b+ABX*I_ERI_D2x_S_S_S_C1002_b;
  Double I_ERI_Dxy_Px_S_S_C1002_b = I_ERI_F2xy_S_S_S_C1002_b+ABX*I_ERI_Dxy_S_S_S_C1002_b;
  Double I_ERI_Dxz_Px_S_S_C1002_b = I_ERI_F2xz_S_S_S_C1002_b+ABX*I_ERI_Dxz_S_S_S_C1002_b;
  Double I_ERI_D2y_Px_S_S_C1002_b = I_ERI_Fx2y_S_S_S_C1002_b+ABX*I_ERI_D2y_S_S_S_C1002_b;
  Double I_ERI_Dyz_Px_S_S_C1002_b = I_ERI_Fxyz_S_S_S_C1002_b+ABX*I_ERI_Dyz_S_S_S_C1002_b;
  Double I_ERI_D2z_Px_S_S_C1002_b = I_ERI_Fx2z_S_S_S_C1002_b+ABX*I_ERI_D2z_S_S_S_C1002_b;
  Double I_ERI_D2x_Py_S_S_C1002_b = I_ERI_F2xy_S_S_S_C1002_b+ABY*I_ERI_D2x_S_S_S_C1002_b;
  Double I_ERI_Dxy_Py_S_S_C1002_b = I_ERI_Fx2y_S_S_S_C1002_b+ABY*I_ERI_Dxy_S_S_S_C1002_b;
  Double I_ERI_Dxz_Py_S_S_C1002_b = I_ERI_Fxyz_S_S_S_C1002_b+ABY*I_ERI_Dxz_S_S_S_C1002_b;
  Double I_ERI_D2y_Py_S_S_C1002_b = I_ERI_F3y_S_S_S_C1002_b+ABY*I_ERI_D2y_S_S_S_C1002_b;
  Double I_ERI_Dyz_Py_S_S_C1002_b = I_ERI_F2yz_S_S_S_C1002_b+ABY*I_ERI_Dyz_S_S_S_C1002_b;
  Double I_ERI_D2z_Py_S_S_C1002_b = I_ERI_Fy2z_S_S_S_C1002_b+ABY*I_ERI_D2z_S_S_S_C1002_b;
  Double I_ERI_D2x_Pz_S_S_C1002_b = I_ERI_F2xz_S_S_S_C1002_b+ABZ*I_ERI_D2x_S_S_S_C1002_b;
  Double I_ERI_Dxy_Pz_S_S_C1002_b = I_ERI_Fxyz_S_S_S_C1002_b+ABZ*I_ERI_Dxy_S_S_S_C1002_b;
  Double I_ERI_Dxz_Pz_S_S_C1002_b = I_ERI_Fx2z_S_S_S_C1002_b+ABZ*I_ERI_Dxz_S_S_S_C1002_b;
  Double I_ERI_D2y_Pz_S_S_C1002_b = I_ERI_F2yz_S_S_S_C1002_b+ABZ*I_ERI_D2y_S_S_S_C1002_b;
  Double I_ERI_Dyz_Pz_S_S_C1002_b = I_ERI_Fy2z_S_S_S_C1002_b+ABZ*I_ERI_Dyz_S_S_S_C1002_b;
  Double I_ERI_D2z_Pz_S_S_C1002_b = I_ERI_F3z_S_S_S_C1002_b+ABZ*I_ERI_D2z_S_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1002_b
   ************************************************************/
  Double I_ERI_Px_D2x_S_S_C1002_b = I_ERI_D2x_Px_S_S_C1002_b+ABX*I_ERI_Px_Px_S_S_C1002_b;
  Double I_ERI_Py_D2x_S_S_C1002_b = I_ERI_Dxy_Px_S_S_C1002_b+ABX*I_ERI_Py_Px_S_S_C1002_b;
  Double I_ERI_Pz_D2x_S_S_C1002_b = I_ERI_Dxz_Px_S_S_C1002_b+ABX*I_ERI_Pz_Px_S_S_C1002_b;
  Double I_ERI_Px_Dxy_S_S_C1002_b = I_ERI_Dxy_Px_S_S_C1002_b+ABY*I_ERI_Px_Px_S_S_C1002_b;
  Double I_ERI_Py_Dxy_S_S_C1002_b = I_ERI_D2y_Px_S_S_C1002_b+ABY*I_ERI_Py_Px_S_S_C1002_b;
  Double I_ERI_Pz_Dxy_S_S_C1002_b = I_ERI_Dyz_Px_S_S_C1002_b+ABY*I_ERI_Pz_Px_S_S_C1002_b;
  Double I_ERI_Px_Dxz_S_S_C1002_b = I_ERI_Dxz_Px_S_S_C1002_b+ABZ*I_ERI_Px_Px_S_S_C1002_b;
  Double I_ERI_Py_Dxz_S_S_C1002_b = I_ERI_Dyz_Px_S_S_C1002_b+ABZ*I_ERI_Py_Px_S_S_C1002_b;
  Double I_ERI_Pz_Dxz_S_S_C1002_b = I_ERI_D2z_Px_S_S_C1002_b+ABZ*I_ERI_Pz_Px_S_S_C1002_b;
  Double I_ERI_Px_D2y_S_S_C1002_b = I_ERI_Dxy_Py_S_S_C1002_b+ABY*I_ERI_Px_Py_S_S_C1002_b;
  Double I_ERI_Py_D2y_S_S_C1002_b = I_ERI_D2y_Py_S_S_C1002_b+ABY*I_ERI_Py_Py_S_S_C1002_b;
  Double I_ERI_Pz_D2y_S_S_C1002_b = I_ERI_Dyz_Py_S_S_C1002_b+ABY*I_ERI_Pz_Py_S_S_C1002_b;
  Double I_ERI_Px_Dyz_S_S_C1002_b = I_ERI_Dxz_Py_S_S_C1002_b+ABZ*I_ERI_Px_Py_S_S_C1002_b;
  Double I_ERI_Py_Dyz_S_S_C1002_b = I_ERI_Dyz_Py_S_S_C1002_b+ABZ*I_ERI_Py_Py_S_S_C1002_b;
  Double I_ERI_Pz_Dyz_S_S_C1002_b = I_ERI_D2z_Py_S_S_C1002_b+ABZ*I_ERI_Pz_Py_S_S_C1002_b;
  Double I_ERI_Px_D2z_S_S_C1002_b = I_ERI_Dxz_Pz_S_S_C1002_b+ABZ*I_ERI_Px_Pz_S_S_C1002_b;
  Double I_ERI_Py_D2z_S_S_C1002_b = I_ERI_Dyz_Pz_S_S_C1002_b+ABZ*I_ERI_Py_Pz_S_S_C1002_b;
  Double I_ERI_Pz_D2z_S_S_C1002_b = I_ERI_D2z_Pz_S_S_C1002_b+ABZ*I_ERI_Pz_Pz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1002_c
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1002_c = I_ERI_D2x_S_Px_S_C1002_c+ABX*I_ERI_Px_S_Px_S_C1002_c;
  Double I_ERI_Py_Px_Px_S_C1002_c = I_ERI_Dxy_S_Px_S_C1002_c+ABX*I_ERI_Py_S_Px_S_C1002_c;
  Double I_ERI_Pz_Px_Px_S_C1002_c = I_ERI_Dxz_S_Px_S_C1002_c+ABX*I_ERI_Pz_S_Px_S_C1002_c;
  Double I_ERI_Px_Py_Px_S_C1002_c = I_ERI_Dxy_S_Px_S_C1002_c+ABY*I_ERI_Px_S_Px_S_C1002_c;
  Double I_ERI_Py_Py_Px_S_C1002_c = I_ERI_D2y_S_Px_S_C1002_c+ABY*I_ERI_Py_S_Px_S_C1002_c;
  Double I_ERI_Pz_Py_Px_S_C1002_c = I_ERI_Dyz_S_Px_S_C1002_c+ABY*I_ERI_Pz_S_Px_S_C1002_c;
  Double I_ERI_Px_Pz_Px_S_C1002_c = I_ERI_Dxz_S_Px_S_C1002_c+ABZ*I_ERI_Px_S_Px_S_C1002_c;
  Double I_ERI_Py_Pz_Px_S_C1002_c = I_ERI_Dyz_S_Px_S_C1002_c+ABZ*I_ERI_Py_S_Px_S_C1002_c;
  Double I_ERI_Pz_Pz_Px_S_C1002_c = I_ERI_D2z_S_Px_S_C1002_c+ABZ*I_ERI_Pz_S_Px_S_C1002_c;
  Double I_ERI_Px_Px_Py_S_C1002_c = I_ERI_D2x_S_Py_S_C1002_c+ABX*I_ERI_Px_S_Py_S_C1002_c;
  Double I_ERI_Py_Px_Py_S_C1002_c = I_ERI_Dxy_S_Py_S_C1002_c+ABX*I_ERI_Py_S_Py_S_C1002_c;
  Double I_ERI_Pz_Px_Py_S_C1002_c = I_ERI_Dxz_S_Py_S_C1002_c+ABX*I_ERI_Pz_S_Py_S_C1002_c;
  Double I_ERI_Px_Py_Py_S_C1002_c = I_ERI_Dxy_S_Py_S_C1002_c+ABY*I_ERI_Px_S_Py_S_C1002_c;
  Double I_ERI_Py_Py_Py_S_C1002_c = I_ERI_D2y_S_Py_S_C1002_c+ABY*I_ERI_Py_S_Py_S_C1002_c;
  Double I_ERI_Pz_Py_Py_S_C1002_c = I_ERI_Dyz_S_Py_S_C1002_c+ABY*I_ERI_Pz_S_Py_S_C1002_c;
  Double I_ERI_Px_Pz_Py_S_C1002_c = I_ERI_Dxz_S_Py_S_C1002_c+ABZ*I_ERI_Px_S_Py_S_C1002_c;
  Double I_ERI_Py_Pz_Py_S_C1002_c = I_ERI_Dyz_S_Py_S_C1002_c+ABZ*I_ERI_Py_S_Py_S_C1002_c;
  Double I_ERI_Pz_Pz_Py_S_C1002_c = I_ERI_D2z_S_Py_S_C1002_c+ABZ*I_ERI_Pz_S_Py_S_C1002_c;
  Double I_ERI_Px_Px_Pz_S_C1002_c = I_ERI_D2x_S_Pz_S_C1002_c+ABX*I_ERI_Px_S_Pz_S_C1002_c;
  Double I_ERI_Py_Px_Pz_S_C1002_c = I_ERI_Dxy_S_Pz_S_C1002_c+ABX*I_ERI_Py_S_Pz_S_C1002_c;
  Double I_ERI_Pz_Px_Pz_S_C1002_c = I_ERI_Dxz_S_Pz_S_C1002_c+ABX*I_ERI_Pz_S_Pz_S_C1002_c;
  Double I_ERI_Px_Py_Pz_S_C1002_c = I_ERI_Dxy_S_Pz_S_C1002_c+ABY*I_ERI_Px_S_Pz_S_C1002_c;
  Double I_ERI_Py_Py_Pz_S_C1002_c = I_ERI_D2y_S_Pz_S_C1002_c+ABY*I_ERI_Py_S_Pz_S_C1002_c;
  Double I_ERI_Pz_Py_Pz_S_C1002_c = I_ERI_Dyz_S_Pz_S_C1002_c+ABY*I_ERI_Pz_S_Pz_S_C1002_c;
  Double I_ERI_Px_Pz_Pz_S_C1002_c = I_ERI_Dxz_S_Pz_S_C1002_c+ABZ*I_ERI_Px_S_Pz_S_C1002_c;
  Double I_ERI_Py_Pz_Pz_S_C1002_c = I_ERI_Dyz_S_Pz_S_C1002_c+ABZ*I_ERI_Py_S_Pz_S_C1002_c;
  Double I_ERI_Pz_Pz_Pz_S_C1002_c = I_ERI_D2z_S_Pz_S_C1002_c+ABZ*I_ERI_Pz_S_Pz_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002_c
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1002_c = I_ERI_F3x_S_S_S_C1002_c+ABX*I_ERI_D2x_S_S_S_C1002_c;
  Double I_ERI_Dxy_Px_S_S_C1002_c = I_ERI_F2xy_S_S_S_C1002_c+ABX*I_ERI_Dxy_S_S_S_C1002_c;
  Double I_ERI_Dxz_Px_S_S_C1002_c = I_ERI_F2xz_S_S_S_C1002_c+ABX*I_ERI_Dxz_S_S_S_C1002_c;
  Double I_ERI_D2y_Px_S_S_C1002_c = I_ERI_Fx2y_S_S_S_C1002_c+ABX*I_ERI_D2y_S_S_S_C1002_c;
  Double I_ERI_Dyz_Px_S_S_C1002_c = I_ERI_Fxyz_S_S_S_C1002_c+ABX*I_ERI_Dyz_S_S_S_C1002_c;
  Double I_ERI_D2z_Px_S_S_C1002_c = I_ERI_Fx2z_S_S_S_C1002_c+ABX*I_ERI_D2z_S_S_S_C1002_c;
  Double I_ERI_D2x_Py_S_S_C1002_c = I_ERI_F2xy_S_S_S_C1002_c+ABY*I_ERI_D2x_S_S_S_C1002_c;
  Double I_ERI_Dxy_Py_S_S_C1002_c = I_ERI_Fx2y_S_S_S_C1002_c+ABY*I_ERI_Dxy_S_S_S_C1002_c;
  Double I_ERI_Dxz_Py_S_S_C1002_c = I_ERI_Fxyz_S_S_S_C1002_c+ABY*I_ERI_Dxz_S_S_S_C1002_c;
  Double I_ERI_D2y_Py_S_S_C1002_c = I_ERI_F3y_S_S_S_C1002_c+ABY*I_ERI_D2y_S_S_S_C1002_c;
  Double I_ERI_Dyz_Py_S_S_C1002_c = I_ERI_F2yz_S_S_S_C1002_c+ABY*I_ERI_Dyz_S_S_S_C1002_c;
  Double I_ERI_D2z_Py_S_S_C1002_c = I_ERI_Fy2z_S_S_S_C1002_c+ABY*I_ERI_D2z_S_S_S_C1002_c;
  Double I_ERI_D2x_Pz_S_S_C1002_c = I_ERI_F2xz_S_S_S_C1002_c+ABZ*I_ERI_D2x_S_S_S_C1002_c;
  Double I_ERI_Dxy_Pz_S_S_C1002_c = I_ERI_Fxyz_S_S_S_C1002_c+ABZ*I_ERI_Dxy_S_S_S_C1002_c;
  Double I_ERI_Dxz_Pz_S_S_C1002_c = I_ERI_Fx2z_S_S_S_C1002_c+ABZ*I_ERI_Dxz_S_S_S_C1002_c;
  Double I_ERI_D2y_Pz_S_S_C1002_c = I_ERI_F2yz_S_S_S_C1002_c+ABZ*I_ERI_D2y_S_S_S_C1002_c;
  Double I_ERI_Dyz_Pz_S_S_C1002_c = I_ERI_Fy2z_S_S_S_C1002_c+ABZ*I_ERI_Dyz_S_S_S_C1002_c;
  Double I_ERI_D2z_Pz_S_S_C1002_c = I_ERI_F3z_S_S_S_C1002_c+ABZ*I_ERI_D2z_S_S_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_aa
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_C1002_aa = I_ERI_H5x_S_S_S_C1002_aa+ABX*I_ERI_G4x_S_S_S_C1002_aa;
  Double I_ERI_G3xy_Px_S_S_C1002_aa = I_ERI_H4xy_S_S_S_C1002_aa+ABX*I_ERI_G3xy_S_S_S_C1002_aa;
  Double I_ERI_G3xz_Px_S_S_C1002_aa = I_ERI_H4xz_S_S_S_C1002_aa+ABX*I_ERI_G3xz_S_S_S_C1002_aa;
  Double I_ERI_G2x2y_Px_S_S_C1002_aa = I_ERI_H3x2y_S_S_S_C1002_aa+ABX*I_ERI_G2x2y_S_S_S_C1002_aa;
  Double I_ERI_G2xyz_Px_S_S_C1002_aa = I_ERI_H3xyz_S_S_S_C1002_aa+ABX*I_ERI_G2xyz_S_S_S_C1002_aa;
  Double I_ERI_G2x2z_Px_S_S_C1002_aa = I_ERI_H3x2z_S_S_S_C1002_aa+ABX*I_ERI_G2x2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3y_Px_S_S_C1002_aa = I_ERI_H2x3y_S_S_S_C1002_aa+ABX*I_ERI_Gx3y_S_S_S_C1002_aa;
  Double I_ERI_Gx2yz_Px_S_S_C1002_aa = I_ERI_H2x2yz_S_S_S_C1002_aa+ABX*I_ERI_Gx2yz_S_S_S_C1002_aa;
  Double I_ERI_Gxy2z_Px_S_S_C1002_aa = I_ERI_H2xy2z_S_S_S_C1002_aa+ABX*I_ERI_Gxy2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3z_Px_S_S_C1002_aa = I_ERI_H2x3z_S_S_S_C1002_aa+ABX*I_ERI_Gx3z_S_S_S_C1002_aa;
  Double I_ERI_G4y_Px_S_S_C1002_aa = I_ERI_Hx4y_S_S_S_C1002_aa+ABX*I_ERI_G4y_S_S_S_C1002_aa;
  Double I_ERI_G3yz_Px_S_S_C1002_aa = I_ERI_Hx3yz_S_S_S_C1002_aa+ABX*I_ERI_G3yz_S_S_S_C1002_aa;
  Double I_ERI_G2y2z_Px_S_S_C1002_aa = I_ERI_Hx2y2z_S_S_S_C1002_aa+ABX*I_ERI_G2y2z_S_S_S_C1002_aa;
  Double I_ERI_Gy3z_Px_S_S_C1002_aa = I_ERI_Hxy3z_S_S_S_C1002_aa+ABX*I_ERI_Gy3z_S_S_S_C1002_aa;
  Double I_ERI_G4z_Px_S_S_C1002_aa = I_ERI_Hx4z_S_S_S_C1002_aa+ABX*I_ERI_G4z_S_S_S_C1002_aa;
  Double I_ERI_G4x_Py_S_S_C1002_aa = I_ERI_H4xy_S_S_S_C1002_aa+ABY*I_ERI_G4x_S_S_S_C1002_aa;
  Double I_ERI_G3xy_Py_S_S_C1002_aa = I_ERI_H3x2y_S_S_S_C1002_aa+ABY*I_ERI_G3xy_S_S_S_C1002_aa;
  Double I_ERI_G3xz_Py_S_S_C1002_aa = I_ERI_H3xyz_S_S_S_C1002_aa+ABY*I_ERI_G3xz_S_S_S_C1002_aa;
  Double I_ERI_G2x2y_Py_S_S_C1002_aa = I_ERI_H2x3y_S_S_S_C1002_aa+ABY*I_ERI_G2x2y_S_S_S_C1002_aa;
  Double I_ERI_G2xyz_Py_S_S_C1002_aa = I_ERI_H2x2yz_S_S_S_C1002_aa+ABY*I_ERI_G2xyz_S_S_S_C1002_aa;
  Double I_ERI_G2x2z_Py_S_S_C1002_aa = I_ERI_H2xy2z_S_S_S_C1002_aa+ABY*I_ERI_G2x2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3y_Py_S_S_C1002_aa = I_ERI_Hx4y_S_S_S_C1002_aa+ABY*I_ERI_Gx3y_S_S_S_C1002_aa;
  Double I_ERI_Gx2yz_Py_S_S_C1002_aa = I_ERI_Hx3yz_S_S_S_C1002_aa+ABY*I_ERI_Gx2yz_S_S_S_C1002_aa;
  Double I_ERI_Gxy2z_Py_S_S_C1002_aa = I_ERI_Hx2y2z_S_S_S_C1002_aa+ABY*I_ERI_Gxy2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3z_Py_S_S_C1002_aa = I_ERI_Hxy3z_S_S_S_C1002_aa+ABY*I_ERI_Gx3z_S_S_S_C1002_aa;
  Double I_ERI_G4y_Py_S_S_C1002_aa = I_ERI_H5y_S_S_S_C1002_aa+ABY*I_ERI_G4y_S_S_S_C1002_aa;
  Double I_ERI_G3yz_Py_S_S_C1002_aa = I_ERI_H4yz_S_S_S_C1002_aa+ABY*I_ERI_G3yz_S_S_S_C1002_aa;
  Double I_ERI_G2y2z_Py_S_S_C1002_aa = I_ERI_H3y2z_S_S_S_C1002_aa+ABY*I_ERI_G2y2z_S_S_S_C1002_aa;
  Double I_ERI_Gy3z_Py_S_S_C1002_aa = I_ERI_H2y3z_S_S_S_C1002_aa+ABY*I_ERI_Gy3z_S_S_S_C1002_aa;
  Double I_ERI_G4z_Py_S_S_C1002_aa = I_ERI_Hy4z_S_S_S_C1002_aa+ABY*I_ERI_G4z_S_S_S_C1002_aa;
  Double I_ERI_G4x_Pz_S_S_C1002_aa = I_ERI_H4xz_S_S_S_C1002_aa+ABZ*I_ERI_G4x_S_S_S_C1002_aa;
  Double I_ERI_G3xy_Pz_S_S_C1002_aa = I_ERI_H3xyz_S_S_S_C1002_aa+ABZ*I_ERI_G3xy_S_S_S_C1002_aa;
  Double I_ERI_G3xz_Pz_S_S_C1002_aa = I_ERI_H3x2z_S_S_S_C1002_aa+ABZ*I_ERI_G3xz_S_S_S_C1002_aa;
  Double I_ERI_G2x2y_Pz_S_S_C1002_aa = I_ERI_H2x2yz_S_S_S_C1002_aa+ABZ*I_ERI_G2x2y_S_S_S_C1002_aa;
  Double I_ERI_G2xyz_Pz_S_S_C1002_aa = I_ERI_H2xy2z_S_S_S_C1002_aa+ABZ*I_ERI_G2xyz_S_S_S_C1002_aa;
  Double I_ERI_G2x2z_Pz_S_S_C1002_aa = I_ERI_H2x3z_S_S_S_C1002_aa+ABZ*I_ERI_G2x2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3y_Pz_S_S_C1002_aa = I_ERI_Hx3yz_S_S_S_C1002_aa+ABZ*I_ERI_Gx3y_S_S_S_C1002_aa;
  Double I_ERI_Gx2yz_Pz_S_S_C1002_aa = I_ERI_Hx2y2z_S_S_S_C1002_aa+ABZ*I_ERI_Gx2yz_S_S_S_C1002_aa;
  Double I_ERI_Gxy2z_Pz_S_S_C1002_aa = I_ERI_Hxy3z_S_S_S_C1002_aa+ABZ*I_ERI_Gxy2z_S_S_S_C1002_aa;
  Double I_ERI_Gx3z_Pz_S_S_C1002_aa = I_ERI_Hx4z_S_S_S_C1002_aa+ABZ*I_ERI_Gx3z_S_S_S_C1002_aa;
  Double I_ERI_G4y_Pz_S_S_C1002_aa = I_ERI_H4yz_S_S_S_C1002_aa+ABZ*I_ERI_G4y_S_S_S_C1002_aa;
  Double I_ERI_G3yz_Pz_S_S_C1002_aa = I_ERI_H3y2z_S_S_S_C1002_aa+ABZ*I_ERI_G3yz_S_S_S_C1002_aa;
  Double I_ERI_G2y2z_Pz_S_S_C1002_aa = I_ERI_H2y3z_S_S_S_C1002_aa+ABZ*I_ERI_G2y2z_S_S_S_C1002_aa;
  Double I_ERI_Gy3z_Pz_S_S_C1002_aa = I_ERI_Hy4z_S_S_S_C1002_aa+ABZ*I_ERI_Gy3z_S_S_S_C1002_aa;
  Double I_ERI_G4z_Pz_S_S_C1002_aa = I_ERI_H5z_S_S_S_C1002_aa+ABZ*I_ERI_G4z_S_S_S_C1002_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_ab
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C2_ab = I_ERI_G4x_S_S_S_C2_ab+ABX*I_ERI_F3x_S_S_S_C2_ab;
  Double I_ERI_F2xy_Px_S_S_C2_ab = I_ERI_G3xy_S_S_S_C2_ab+ABX*I_ERI_F2xy_S_S_S_C2_ab;
  Double I_ERI_F2xz_Px_S_S_C2_ab = I_ERI_G3xz_S_S_S_C2_ab+ABX*I_ERI_F2xz_S_S_S_C2_ab;
  Double I_ERI_Fx2y_Px_S_S_C2_ab = I_ERI_G2x2y_S_S_S_C2_ab+ABX*I_ERI_Fx2y_S_S_S_C2_ab;
  Double I_ERI_Fxyz_Px_S_S_C2_ab = I_ERI_G2xyz_S_S_S_C2_ab+ABX*I_ERI_Fxyz_S_S_S_C2_ab;
  Double I_ERI_Fx2z_Px_S_S_C2_ab = I_ERI_G2x2z_S_S_S_C2_ab+ABX*I_ERI_Fx2z_S_S_S_C2_ab;
  Double I_ERI_F3y_Px_S_S_C2_ab = I_ERI_Gx3y_S_S_S_C2_ab+ABX*I_ERI_F3y_S_S_S_C2_ab;
  Double I_ERI_F2yz_Px_S_S_C2_ab = I_ERI_Gx2yz_S_S_S_C2_ab+ABX*I_ERI_F2yz_S_S_S_C2_ab;
  Double I_ERI_Fy2z_Px_S_S_C2_ab = I_ERI_Gxy2z_S_S_S_C2_ab+ABX*I_ERI_Fy2z_S_S_S_C2_ab;
  Double I_ERI_F3z_Px_S_S_C2_ab = I_ERI_Gx3z_S_S_S_C2_ab+ABX*I_ERI_F3z_S_S_S_C2_ab;
  Double I_ERI_F3x_Py_S_S_C2_ab = I_ERI_G3xy_S_S_S_C2_ab+ABY*I_ERI_F3x_S_S_S_C2_ab;
  Double I_ERI_F2xy_Py_S_S_C2_ab = I_ERI_G2x2y_S_S_S_C2_ab+ABY*I_ERI_F2xy_S_S_S_C2_ab;
  Double I_ERI_F2xz_Py_S_S_C2_ab = I_ERI_G2xyz_S_S_S_C2_ab+ABY*I_ERI_F2xz_S_S_S_C2_ab;
  Double I_ERI_Fx2y_Py_S_S_C2_ab = I_ERI_Gx3y_S_S_S_C2_ab+ABY*I_ERI_Fx2y_S_S_S_C2_ab;
  Double I_ERI_Fxyz_Py_S_S_C2_ab = I_ERI_Gx2yz_S_S_S_C2_ab+ABY*I_ERI_Fxyz_S_S_S_C2_ab;
  Double I_ERI_Fx2z_Py_S_S_C2_ab = I_ERI_Gxy2z_S_S_S_C2_ab+ABY*I_ERI_Fx2z_S_S_S_C2_ab;
  Double I_ERI_F3y_Py_S_S_C2_ab = I_ERI_G4y_S_S_S_C2_ab+ABY*I_ERI_F3y_S_S_S_C2_ab;
  Double I_ERI_F2yz_Py_S_S_C2_ab = I_ERI_G3yz_S_S_S_C2_ab+ABY*I_ERI_F2yz_S_S_S_C2_ab;
  Double I_ERI_Fy2z_Py_S_S_C2_ab = I_ERI_G2y2z_S_S_S_C2_ab+ABY*I_ERI_Fy2z_S_S_S_C2_ab;
  Double I_ERI_F3z_Py_S_S_C2_ab = I_ERI_Gy3z_S_S_S_C2_ab+ABY*I_ERI_F3z_S_S_S_C2_ab;
  Double I_ERI_F3x_Pz_S_S_C2_ab = I_ERI_G3xz_S_S_S_C2_ab+ABZ*I_ERI_F3x_S_S_S_C2_ab;
  Double I_ERI_F2xy_Pz_S_S_C2_ab = I_ERI_G2xyz_S_S_S_C2_ab+ABZ*I_ERI_F2xy_S_S_S_C2_ab;
  Double I_ERI_F2xz_Pz_S_S_C2_ab = I_ERI_G2x2z_S_S_S_C2_ab+ABZ*I_ERI_F2xz_S_S_S_C2_ab;
  Double I_ERI_Fx2y_Pz_S_S_C2_ab = I_ERI_Gx2yz_S_S_S_C2_ab+ABZ*I_ERI_Fx2y_S_S_S_C2_ab;
  Double I_ERI_Fxyz_Pz_S_S_C2_ab = I_ERI_Gxy2z_S_S_S_C2_ab+ABZ*I_ERI_Fxyz_S_S_S_C2_ab;
  Double I_ERI_Fx2z_Pz_S_S_C2_ab = I_ERI_Gx3z_S_S_S_C2_ab+ABZ*I_ERI_Fx2z_S_S_S_C2_ab;
  Double I_ERI_F3y_Pz_S_S_C2_ab = I_ERI_G3yz_S_S_S_C2_ab+ABZ*I_ERI_F3y_S_S_S_C2_ab;
  Double I_ERI_F2yz_Pz_S_S_C2_ab = I_ERI_G2y2z_S_S_S_C2_ab+ABZ*I_ERI_F2yz_S_S_S_C2_ab;
  Double I_ERI_Fy2z_Pz_S_S_C2_ab = I_ERI_Gy3z_S_S_S_C2_ab+ABZ*I_ERI_Fy2z_S_S_S_C2_ab;
  Double I_ERI_F3z_Pz_S_S_C2_ab = I_ERI_G4z_S_S_S_C2_ab+ABZ*I_ERI_F3z_S_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1002_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_ab
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1002_ab = I_ERI_G4x_S_S_S_C1002_ab+ABX*I_ERI_F3x_S_S_S_C1002_ab;
  Double I_ERI_F2xy_Px_S_S_C1002_ab = I_ERI_G3xy_S_S_S_C1002_ab+ABX*I_ERI_F2xy_S_S_S_C1002_ab;
  Double I_ERI_F2xz_Px_S_S_C1002_ab = I_ERI_G3xz_S_S_S_C1002_ab+ABX*I_ERI_F2xz_S_S_S_C1002_ab;
  Double I_ERI_Fx2y_Px_S_S_C1002_ab = I_ERI_G2x2y_S_S_S_C1002_ab+ABX*I_ERI_Fx2y_S_S_S_C1002_ab;
  Double I_ERI_Fxyz_Px_S_S_C1002_ab = I_ERI_G2xyz_S_S_S_C1002_ab+ABX*I_ERI_Fxyz_S_S_S_C1002_ab;
  Double I_ERI_Fx2z_Px_S_S_C1002_ab = I_ERI_G2x2z_S_S_S_C1002_ab+ABX*I_ERI_Fx2z_S_S_S_C1002_ab;
  Double I_ERI_F3y_Px_S_S_C1002_ab = I_ERI_Gx3y_S_S_S_C1002_ab+ABX*I_ERI_F3y_S_S_S_C1002_ab;
  Double I_ERI_F2yz_Px_S_S_C1002_ab = I_ERI_Gx2yz_S_S_S_C1002_ab+ABX*I_ERI_F2yz_S_S_S_C1002_ab;
  Double I_ERI_Fy2z_Px_S_S_C1002_ab = I_ERI_Gxy2z_S_S_S_C1002_ab+ABX*I_ERI_Fy2z_S_S_S_C1002_ab;
  Double I_ERI_F3z_Px_S_S_C1002_ab = I_ERI_Gx3z_S_S_S_C1002_ab+ABX*I_ERI_F3z_S_S_S_C1002_ab;
  Double I_ERI_F3x_Py_S_S_C1002_ab = I_ERI_G3xy_S_S_S_C1002_ab+ABY*I_ERI_F3x_S_S_S_C1002_ab;
  Double I_ERI_F2xy_Py_S_S_C1002_ab = I_ERI_G2x2y_S_S_S_C1002_ab+ABY*I_ERI_F2xy_S_S_S_C1002_ab;
  Double I_ERI_F2xz_Py_S_S_C1002_ab = I_ERI_G2xyz_S_S_S_C1002_ab+ABY*I_ERI_F2xz_S_S_S_C1002_ab;
  Double I_ERI_Fx2y_Py_S_S_C1002_ab = I_ERI_Gx3y_S_S_S_C1002_ab+ABY*I_ERI_Fx2y_S_S_S_C1002_ab;
  Double I_ERI_Fxyz_Py_S_S_C1002_ab = I_ERI_Gx2yz_S_S_S_C1002_ab+ABY*I_ERI_Fxyz_S_S_S_C1002_ab;
  Double I_ERI_Fx2z_Py_S_S_C1002_ab = I_ERI_Gxy2z_S_S_S_C1002_ab+ABY*I_ERI_Fx2z_S_S_S_C1002_ab;
  Double I_ERI_F3y_Py_S_S_C1002_ab = I_ERI_G4y_S_S_S_C1002_ab+ABY*I_ERI_F3y_S_S_S_C1002_ab;
  Double I_ERI_F2yz_Py_S_S_C1002_ab = I_ERI_G3yz_S_S_S_C1002_ab+ABY*I_ERI_F2yz_S_S_S_C1002_ab;
  Double I_ERI_Fy2z_Py_S_S_C1002_ab = I_ERI_G2y2z_S_S_S_C1002_ab+ABY*I_ERI_Fy2z_S_S_S_C1002_ab;
  Double I_ERI_F3z_Py_S_S_C1002_ab = I_ERI_Gy3z_S_S_S_C1002_ab+ABY*I_ERI_F3z_S_S_S_C1002_ab;
  Double I_ERI_F3x_Pz_S_S_C1002_ab = I_ERI_G3xz_S_S_S_C1002_ab+ABZ*I_ERI_F3x_S_S_S_C1002_ab;
  Double I_ERI_F2xy_Pz_S_S_C1002_ab = I_ERI_G2xyz_S_S_S_C1002_ab+ABZ*I_ERI_F2xy_S_S_S_C1002_ab;
  Double I_ERI_F2xz_Pz_S_S_C1002_ab = I_ERI_G2x2z_S_S_S_C1002_ab+ABZ*I_ERI_F2xz_S_S_S_C1002_ab;
  Double I_ERI_Fx2y_Pz_S_S_C1002_ab = I_ERI_Gx2yz_S_S_S_C1002_ab+ABZ*I_ERI_Fx2y_S_S_S_C1002_ab;
  Double I_ERI_Fxyz_Pz_S_S_C1002_ab = I_ERI_Gxy2z_S_S_S_C1002_ab+ABZ*I_ERI_Fxyz_S_S_S_C1002_ab;
  Double I_ERI_Fx2z_Pz_S_S_C1002_ab = I_ERI_Gx3z_S_S_S_C1002_ab+ABZ*I_ERI_Fx2z_S_S_S_C1002_ab;
  Double I_ERI_F3y_Pz_S_S_C1002_ab = I_ERI_G3yz_S_S_S_C1002_ab+ABZ*I_ERI_F3y_S_S_S_C1002_ab;
  Double I_ERI_F2yz_Pz_S_S_C1002_ab = I_ERI_G2y2z_S_S_S_C1002_ab+ABZ*I_ERI_F2yz_S_S_S_C1002_ab;
  Double I_ERI_Fy2z_Pz_S_S_C1002_ab = I_ERI_Gy3z_S_S_S_C1002_ab+ABZ*I_ERI_Fy2z_S_S_S_C1002_ab;
  Double I_ERI_F3z_Pz_S_S_C1002_ab = I_ERI_G4z_S_S_S_C1002_ab+ABZ*I_ERI_F3z_S_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_C1002_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_ab
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_C1002_ab = I_ERI_H5x_S_S_S_C1002_ab+ABX*I_ERI_G4x_S_S_S_C1002_ab;
  Double I_ERI_G3xy_Px_S_S_C1002_ab = I_ERI_H4xy_S_S_S_C1002_ab+ABX*I_ERI_G3xy_S_S_S_C1002_ab;
  Double I_ERI_G3xz_Px_S_S_C1002_ab = I_ERI_H4xz_S_S_S_C1002_ab+ABX*I_ERI_G3xz_S_S_S_C1002_ab;
  Double I_ERI_G2x2y_Px_S_S_C1002_ab = I_ERI_H3x2y_S_S_S_C1002_ab+ABX*I_ERI_G2x2y_S_S_S_C1002_ab;
  Double I_ERI_G2xyz_Px_S_S_C1002_ab = I_ERI_H3xyz_S_S_S_C1002_ab+ABX*I_ERI_G2xyz_S_S_S_C1002_ab;
  Double I_ERI_G2x2z_Px_S_S_C1002_ab = I_ERI_H3x2z_S_S_S_C1002_ab+ABX*I_ERI_G2x2z_S_S_S_C1002_ab;
  Double I_ERI_Gx3y_Px_S_S_C1002_ab = I_ERI_H2x3y_S_S_S_C1002_ab+ABX*I_ERI_Gx3y_S_S_S_C1002_ab;
  Double I_ERI_Gx2yz_Px_S_S_C1002_ab = I_ERI_H2x2yz_S_S_S_C1002_ab+ABX*I_ERI_Gx2yz_S_S_S_C1002_ab;
  Double I_ERI_Gxy2z_Px_S_S_C1002_ab = I_ERI_H2xy2z_S_S_S_C1002_ab+ABX*I_ERI_Gxy2z_S_S_S_C1002_ab;
  Double I_ERI_Gx3z_Px_S_S_C1002_ab = I_ERI_H2x3z_S_S_S_C1002_ab+ABX*I_ERI_Gx3z_S_S_S_C1002_ab;
  Double I_ERI_G4y_Px_S_S_C1002_ab = I_ERI_Hx4y_S_S_S_C1002_ab+ABX*I_ERI_G4y_S_S_S_C1002_ab;
  Double I_ERI_G3yz_Px_S_S_C1002_ab = I_ERI_Hx3yz_S_S_S_C1002_ab+ABX*I_ERI_G3yz_S_S_S_C1002_ab;
  Double I_ERI_G2y2z_Px_S_S_C1002_ab = I_ERI_Hx2y2z_S_S_S_C1002_ab+ABX*I_ERI_G2y2z_S_S_S_C1002_ab;
  Double I_ERI_Gy3z_Px_S_S_C1002_ab = I_ERI_Hxy3z_S_S_S_C1002_ab+ABX*I_ERI_Gy3z_S_S_S_C1002_ab;
  Double I_ERI_G4z_Px_S_S_C1002_ab = I_ERI_Hx4z_S_S_S_C1002_ab+ABX*I_ERI_G4z_S_S_S_C1002_ab;
  Double I_ERI_G3xy_Py_S_S_C1002_ab = I_ERI_H3x2y_S_S_S_C1002_ab+ABY*I_ERI_G3xy_S_S_S_C1002_ab;
  Double I_ERI_G3xz_Py_S_S_C1002_ab = I_ERI_H3xyz_S_S_S_C1002_ab+ABY*I_ERI_G3xz_S_S_S_C1002_ab;
  Double I_ERI_G2x2y_Py_S_S_C1002_ab = I_ERI_H2x3y_S_S_S_C1002_ab+ABY*I_ERI_G2x2y_S_S_S_C1002_ab;
  Double I_ERI_G2xyz_Py_S_S_C1002_ab = I_ERI_H2x2yz_S_S_S_C1002_ab+ABY*I_ERI_G2xyz_S_S_S_C1002_ab;
  Double I_ERI_G2x2z_Py_S_S_C1002_ab = I_ERI_H2xy2z_S_S_S_C1002_ab+ABY*I_ERI_G2x2z_S_S_S_C1002_ab;
  Double I_ERI_Gx3y_Py_S_S_C1002_ab = I_ERI_Hx4y_S_S_S_C1002_ab+ABY*I_ERI_Gx3y_S_S_S_C1002_ab;
  Double I_ERI_Gx2yz_Py_S_S_C1002_ab = I_ERI_Hx3yz_S_S_S_C1002_ab+ABY*I_ERI_Gx2yz_S_S_S_C1002_ab;
  Double I_ERI_Gxy2z_Py_S_S_C1002_ab = I_ERI_Hx2y2z_S_S_S_C1002_ab+ABY*I_ERI_Gxy2z_S_S_S_C1002_ab;
  Double I_ERI_Gx3z_Py_S_S_C1002_ab = I_ERI_Hxy3z_S_S_S_C1002_ab+ABY*I_ERI_Gx3z_S_S_S_C1002_ab;
  Double I_ERI_G4y_Py_S_S_C1002_ab = I_ERI_H5y_S_S_S_C1002_ab+ABY*I_ERI_G4y_S_S_S_C1002_ab;
  Double I_ERI_G3yz_Py_S_S_C1002_ab = I_ERI_H4yz_S_S_S_C1002_ab+ABY*I_ERI_G3yz_S_S_S_C1002_ab;
  Double I_ERI_G2y2z_Py_S_S_C1002_ab = I_ERI_H3y2z_S_S_S_C1002_ab+ABY*I_ERI_G2y2z_S_S_S_C1002_ab;
  Double I_ERI_Gy3z_Py_S_S_C1002_ab = I_ERI_H2y3z_S_S_S_C1002_ab+ABY*I_ERI_Gy3z_S_S_S_C1002_ab;
  Double I_ERI_G4z_Py_S_S_C1002_ab = I_ERI_Hy4z_S_S_S_C1002_ab+ABY*I_ERI_G4z_S_S_S_C1002_ab;
  Double I_ERI_G3xz_Pz_S_S_C1002_ab = I_ERI_H3x2z_S_S_S_C1002_ab+ABZ*I_ERI_G3xz_S_S_S_C1002_ab;
  Double I_ERI_G2xyz_Pz_S_S_C1002_ab = I_ERI_H2xy2z_S_S_S_C1002_ab+ABZ*I_ERI_G2xyz_S_S_S_C1002_ab;
  Double I_ERI_G2x2z_Pz_S_S_C1002_ab = I_ERI_H2x3z_S_S_S_C1002_ab+ABZ*I_ERI_G2x2z_S_S_S_C1002_ab;
  Double I_ERI_Gx2yz_Pz_S_S_C1002_ab = I_ERI_Hx2y2z_S_S_S_C1002_ab+ABZ*I_ERI_Gx2yz_S_S_S_C1002_ab;
  Double I_ERI_Gxy2z_Pz_S_S_C1002_ab = I_ERI_Hxy3z_S_S_S_C1002_ab+ABZ*I_ERI_Gxy2z_S_S_S_C1002_ab;
  Double I_ERI_Gx3z_Pz_S_S_C1002_ab = I_ERI_Hx4z_S_S_S_C1002_ab+ABZ*I_ERI_Gx3z_S_S_S_C1002_ab;
  Double I_ERI_G3yz_Pz_S_S_C1002_ab = I_ERI_H3y2z_S_S_S_C1002_ab+ABZ*I_ERI_G3yz_S_S_S_C1002_ab;
  Double I_ERI_G2y2z_Pz_S_S_C1002_ab = I_ERI_H2y3z_S_S_S_C1002_ab+ABZ*I_ERI_G2y2z_S_S_S_C1002_ab;
  Double I_ERI_Gy3z_Pz_S_S_C1002_ab = I_ERI_Hy4z_S_S_S_C1002_ab+ABZ*I_ERI_Gy3z_S_S_S_C1002_ab;
  Double I_ERI_G4z_Pz_S_S_C1002_ab = I_ERI_H5z_S_S_S_C1002_ab+ABZ*I_ERI_G4z_S_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_ab
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_C1002_ab = I_ERI_G4x_Px_S_S_C1002_ab+ABX*I_ERI_F3x_Px_S_S_C1002_ab;
  Double I_ERI_F2xy_D2x_S_S_C1002_ab = I_ERI_G3xy_Px_S_S_C1002_ab+ABX*I_ERI_F2xy_Px_S_S_C1002_ab;
  Double I_ERI_F2xz_D2x_S_S_C1002_ab = I_ERI_G3xz_Px_S_S_C1002_ab+ABX*I_ERI_F2xz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2y_D2x_S_S_C1002_ab = I_ERI_G2x2y_Px_S_S_C1002_ab+ABX*I_ERI_Fx2y_Px_S_S_C1002_ab;
  Double I_ERI_Fxyz_D2x_S_S_C1002_ab = I_ERI_G2xyz_Px_S_S_C1002_ab+ABX*I_ERI_Fxyz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2z_D2x_S_S_C1002_ab = I_ERI_G2x2z_Px_S_S_C1002_ab+ABX*I_ERI_Fx2z_Px_S_S_C1002_ab;
  Double I_ERI_F3y_D2x_S_S_C1002_ab = I_ERI_Gx3y_Px_S_S_C1002_ab+ABX*I_ERI_F3y_Px_S_S_C1002_ab;
  Double I_ERI_F2yz_D2x_S_S_C1002_ab = I_ERI_Gx2yz_Px_S_S_C1002_ab+ABX*I_ERI_F2yz_Px_S_S_C1002_ab;
  Double I_ERI_Fy2z_D2x_S_S_C1002_ab = I_ERI_Gxy2z_Px_S_S_C1002_ab+ABX*I_ERI_Fy2z_Px_S_S_C1002_ab;
  Double I_ERI_F3z_D2x_S_S_C1002_ab = I_ERI_Gx3z_Px_S_S_C1002_ab+ABX*I_ERI_F3z_Px_S_S_C1002_ab;
  Double I_ERI_F3x_Dxy_S_S_C1002_ab = I_ERI_G3xy_Px_S_S_C1002_ab+ABY*I_ERI_F3x_Px_S_S_C1002_ab;
  Double I_ERI_F2xy_Dxy_S_S_C1002_ab = I_ERI_G2x2y_Px_S_S_C1002_ab+ABY*I_ERI_F2xy_Px_S_S_C1002_ab;
  Double I_ERI_F2xz_Dxy_S_S_C1002_ab = I_ERI_G2xyz_Px_S_S_C1002_ab+ABY*I_ERI_F2xz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2y_Dxy_S_S_C1002_ab = I_ERI_Gx3y_Px_S_S_C1002_ab+ABY*I_ERI_Fx2y_Px_S_S_C1002_ab;
  Double I_ERI_Fxyz_Dxy_S_S_C1002_ab = I_ERI_Gx2yz_Px_S_S_C1002_ab+ABY*I_ERI_Fxyz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2z_Dxy_S_S_C1002_ab = I_ERI_Gxy2z_Px_S_S_C1002_ab+ABY*I_ERI_Fx2z_Px_S_S_C1002_ab;
  Double I_ERI_F3y_Dxy_S_S_C1002_ab = I_ERI_G4y_Px_S_S_C1002_ab+ABY*I_ERI_F3y_Px_S_S_C1002_ab;
  Double I_ERI_F2yz_Dxy_S_S_C1002_ab = I_ERI_G3yz_Px_S_S_C1002_ab+ABY*I_ERI_F2yz_Px_S_S_C1002_ab;
  Double I_ERI_Fy2z_Dxy_S_S_C1002_ab = I_ERI_G2y2z_Px_S_S_C1002_ab+ABY*I_ERI_Fy2z_Px_S_S_C1002_ab;
  Double I_ERI_F3z_Dxy_S_S_C1002_ab = I_ERI_Gy3z_Px_S_S_C1002_ab+ABY*I_ERI_F3z_Px_S_S_C1002_ab;
  Double I_ERI_F3x_Dxz_S_S_C1002_ab = I_ERI_G3xz_Px_S_S_C1002_ab+ABZ*I_ERI_F3x_Px_S_S_C1002_ab;
  Double I_ERI_F2xy_Dxz_S_S_C1002_ab = I_ERI_G2xyz_Px_S_S_C1002_ab+ABZ*I_ERI_F2xy_Px_S_S_C1002_ab;
  Double I_ERI_F2xz_Dxz_S_S_C1002_ab = I_ERI_G2x2z_Px_S_S_C1002_ab+ABZ*I_ERI_F2xz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2y_Dxz_S_S_C1002_ab = I_ERI_Gx2yz_Px_S_S_C1002_ab+ABZ*I_ERI_Fx2y_Px_S_S_C1002_ab;
  Double I_ERI_Fxyz_Dxz_S_S_C1002_ab = I_ERI_Gxy2z_Px_S_S_C1002_ab+ABZ*I_ERI_Fxyz_Px_S_S_C1002_ab;
  Double I_ERI_Fx2z_Dxz_S_S_C1002_ab = I_ERI_Gx3z_Px_S_S_C1002_ab+ABZ*I_ERI_Fx2z_Px_S_S_C1002_ab;
  Double I_ERI_F3y_Dxz_S_S_C1002_ab = I_ERI_G3yz_Px_S_S_C1002_ab+ABZ*I_ERI_F3y_Px_S_S_C1002_ab;
  Double I_ERI_F2yz_Dxz_S_S_C1002_ab = I_ERI_G2y2z_Px_S_S_C1002_ab+ABZ*I_ERI_F2yz_Px_S_S_C1002_ab;
  Double I_ERI_Fy2z_Dxz_S_S_C1002_ab = I_ERI_Gy3z_Px_S_S_C1002_ab+ABZ*I_ERI_Fy2z_Px_S_S_C1002_ab;
  Double I_ERI_F3z_Dxz_S_S_C1002_ab = I_ERI_G4z_Px_S_S_C1002_ab+ABZ*I_ERI_F3z_Px_S_S_C1002_ab;
  Double I_ERI_F3x_D2y_S_S_C1002_ab = I_ERI_G3xy_Py_S_S_C1002_ab+ABY*I_ERI_F3x_Py_S_S_C1002_ab;
  Double I_ERI_F2xy_D2y_S_S_C1002_ab = I_ERI_G2x2y_Py_S_S_C1002_ab+ABY*I_ERI_F2xy_Py_S_S_C1002_ab;
  Double I_ERI_F2xz_D2y_S_S_C1002_ab = I_ERI_G2xyz_Py_S_S_C1002_ab+ABY*I_ERI_F2xz_Py_S_S_C1002_ab;
  Double I_ERI_Fx2y_D2y_S_S_C1002_ab = I_ERI_Gx3y_Py_S_S_C1002_ab+ABY*I_ERI_Fx2y_Py_S_S_C1002_ab;
  Double I_ERI_Fxyz_D2y_S_S_C1002_ab = I_ERI_Gx2yz_Py_S_S_C1002_ab+ABY*I_ERI_Fxyz_Py_S_S_C1002_ab;
  Double I_ERI_Fx2z_D2y_S_S_C1002_ab = I_ERI_Gxy2z_Py_S_S_C1002_ab+ABY*I_ERI_Fx2z_Py_S_S_C1002_ab;
  Double I_ERI_F3y_D2y_S_S_C1002_ab = I_ERI_G4y_Py_S_S_C1002_ab+ABY*I_ERI_F3y_Py_S_S_C1002_ab;
  Double I_ERI_F2yz_D2y_S_S_C1002_ab = I_ERI_G3yz_Py_S_S_C1002_ab+ABY*I_ERI_F2yz_Py_S_S_C1002_ab;
  Double I_ERI_Fy2z_D2y_S_S_C1002_ab = I_ERI_G2y2z_Py_S_S_C1002_ab+ABY*I_ERI_Fy2z_Py_S_S_C1002_ab;
  Double I_ERI_F3z_D2y_S_S_C1002_ab = I_ERI_Gy3z_Py_S_S_C1002_ab+ABY*I_ERI_F3z_Py_S_S_C1002_ab;
  Double I_ERI_F3x_Dyz_S_S_C1002_ab = I_ERI_G3xz_Py_S_S_C1002_ab+ABZ*I_ERI_F3x_Py_S_S_C1002_ab;
  Double I_ERI_F2xy_Dyz_S_S_C1002_ab = I_ERI_G2xyz_Py_S_S_C1002_ab+ABZ*I_ERI_F2xy_Py_S_S_C1002_ab;
  Double I_ERI_F2xz_Dyz_S_S_C1002_ab = I_ERI_G2x2z_Py_S_S_C1002_ab+ABZ*I_ERI_F2xz_Py_S_S_C1002_ab;
  Double I_ERI_Fx2y_Dyz_S_S_C1002_ab = I_ERI_Gx2yz_Py_S_S_C1002_ab+ABZ*I_ERI_Fx2y_Py_S_S_C1002_ab;
  Double I_ERI_Fxyz_Dyz_S_S_C1002_ab = I_ERI_Gxy2z_Py_S_S_C1002_ab+ABZ*I_ERI_Fxyz_Py_S_S_C1002_ab;
  Double I_ERI_Fx2z_Dyz_S_S_C1002_ab = I_ERI_Gx3z_Py_S_S_C1002_ab+ABZ*I_ERI_Fx2z_Py_S_S_C1002_ab;
  Double I_ERI_F3y_Dyz_S_S_C1002_ab = I_ERI_G3yz_Py_S_S_C1002_ab+ABZ*I_ERI_F3y_Py_S_S_C1002_ab;
  Double I_ERI_F2yz_Dyz_S_S_C1002_ab = I_ERI_G2y2z_Py_S_S_C1002_ab+ABZ*I_ERI_F2yz_Py_S_S_C1002_ab;
  Double I_ERI_Fy2z_Dyz_S_S_C1002_ab = I_ERI_Gy3z_Py_S_S_C1002_ab+ABZ*I_ERI_Fy2z_Py_S_S_C1002_ab;
  Double I_ERI_F3z_Dyz_S_S_C1002_ab = I_ERI_G4z_Py_S_S_C1002_ab+ABZ*I_ERI_F3z_Py_S_S_C1002_ab;
  Double I_ERI_F3x_D2z_S_S_C1002_ab = I_ERI_G3xz_Pz_S_S_C1002_ab+ABZ*I_ERI_F3x_Pz_S_S_C1002_ab;
  Double I_ERI_F2xy_D2z_S_S_C1002_ab = I_ERI_G2xyz_Pz_S_S_C1002_ab+ABZ*I_ERI_F2xy_Pz_S_S_C1002_ab;
  Double I_ERI_F2xz_D2z_S_S_C1002_ab = I_ERI_G2x2z_Pz_S_S_C1002_ab+ABZ*I_ERI_F2xz_Pz_S_S_C1002_ab;
  Double I_ERI_Fx2y_D2z_S_S_C1002_ab = I_ERI_Gx2yz_Pz_S_S_C1002_ab+ABZ*I_ERI_Fx2y_Pz_S_S_C1002_ab;
  Double I_ERI_Fxyz_D2z_S_S_C1002_ab = I_ERI_Gxy2z_Pz_S_S_C1002_ab+ABZ*I_ERI_Fxyz_Pz_S_S_C1002_ab;
  Double I_ERI_Fx2z_D2z_S_S_C1002_ab = I_ERI_Gx3z_Pz_S_S_C1002_ab+ABZ*I_ERI_Fx2z_Pz_S_S_C1002_ab;
  Double I_ERI_F3y_D2z_S_S_C1002_ab = I_ERI_G3yz_Pz_S_S_C1002_ab+ABZ*I_ERI_F3y_Pz_S_S_C1002_ab;
  Double I_ERI_F2yz_D2z_S_S_C1002_ab = I_ERI_G2y2z_Pz_S_S_C1002_ab+ABZ*I_ERI_F2yz_Pz_S_S_C1002_ab;
  Double I_ERI_Fy2z_D2z_S_S_C1002_ab = I_ERI_Gy3z_Pz_S_S_C1002_ab+ABZ*I_ERI_Fy2z_Pz_S_S_C1002_ab;
  Double I_ERI_F3z_D2z_S_S_C1002_ab = I_ERI_G4z_Pz_S_S_C1002_ab+ABZ*I_ERI_F3z_Pz_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1002_ac
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1002_ac = I_ERI_G4x_S_Px_S_C1002_ac+ABX*I_ERI_F3x_S_Px_S_C1002_ac;
  Double I_ERI_F2xy_Px_Px_S_C1002_ac = I_ERI_G3xy_S_Px_S_C1002_ac+ABX*I_ERI_F2xy_S_Px_S_C1002_ac;
  Double I_ERI_F2xz_Px_Px_S_C1002_ac = I_ERI_G3xz_S_Px_S_C1002_ac+ABX*I_ERI_F2xz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2y_Px_Px_S_C1002_ac = I_ERI_G2x2y_S_Px_S_C1002_ac+ABX*I_ERI_Fx2y_S_Px_S_C1002_ac;
  Double I_ERI_Fxyz_Px_Px_S_C1002_ac = I_ERI_G2xyz_S_Px_S_C1002_ac+ABX*I_ERI_Fxyz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2z_Px_Px_S_C1002_ac = I_ERI_G2x2z_S_Px_S_C1002_ac+ABX*I_ERI_Fx2z_S_Px_S_C1002_ac;
  Double I_ERI_F3y_Px_Px_S_C1002_ac = I_ERI_Gx3y_S_Px_S_C1002_ac+ABX*I_ERI_F3y_S_Px_S_C1002_ac;
  Double I_ERI_F2yz_Px_Px_S_C1002_ac = I_ERI_Gx2yz_S_Px_S_C1002_ac+ABX*I_ERI_F2yz_S_Px_S_C1002_ac;
  Double I_ERI_Fy2z_Px_Px_S_C1002_ac = I_ERI_Gxy2z_S_Px_S_C1002_ac+ABX*I_ERI_Fy2z_S_Px_S_C1002_ac;
  Double I_ERI_F3z_Px_Px_S_C1002_ac = I_ERI_Gx3z_S_Px_S_C1002_ac+ABX*I_ERI_F3z_S_Px_S_C1002_ac;
  Double I_ERI_F3x_Py_Px_S_C1002_ac = I_ERI_G3xy_S_Px_S_C1002_ac+ABY*I_ERI_F3x_S_Px_S_C1002_ac;
  Double I_ERI_F2xy_Py_Px_S_C1002_ac = I_ERI_G2x2y_S_Px_S_C1002_ac+ABY*I_ERI_F2xy_S_Px_S_C1002_ac;
  Double I_ERI_F2xz_Py_Px_S_C1002_ac = I_ERI_G2xyz_S_Px_S_C1002_ac+ABY*I_ERI_F2xz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2y_Py_Px_S_C1002_ac = I_ERI_Gx3y_S_Px_S_C1002_ac+ABY*I_ERI_Fx2y_S_Px_S_C1002_ac;
  Double I_ERI_Fxyz_Py_Px_S_C1002_ac = I_ERI_Gx2yz_S_Px_S_C1002_ac+ABY*I_ERI_Fxyz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2z_Py_Px_S_C1002_ac = I_ERI_Gxy2z_S_Px_S_C1002_ac+ABY*I_ERI_Fx2z_S_Px_S_C1002_ac;
  Double I_ERI_F3y_Py_Px_S_C1002_ac = I_ERI_G4y_S_Px_S_C1002_ac+ABY*I_ERI_F3y_S_Px_S_C1002_ac;
  Double I_ERI_F2yz_Py_Px_S_C1002_ac = I_ERI_G3yz_S_Px_S_C1002_ac+ABY*I_ERI_F2yz_S_Px_S_C1002_ac;
  Double I_ERI_Fy2z_Py_Px_S_C1002_ac = I_ERI_G2y2z_S_Px_S_C1002_ac+ABY*I_ERI_Fy2z_S_Px_S_C1002_ac;
  Double I_ERI_F3z_Py_Px_S_C1002_ac = I_ERI_Gy3z_S_Px_S_C1002_ac+ABY*I_ERI_F3z_S_Px_S_C1002_ac;
  Double I_ERI_F3x_Pz_Px_S_C1002_ac = I_ERI_G3xz_S_Px_S_C1002_ac+ABZ*I_ERI_F3x_S_Px_S_C1002_ac;
  Double I_ERI_F2xy_Pz_Px_S_C1002_ac = I_ERI_G2xyz_S_Px_S_C1002_ac+ABZ*I_ERI_F2xy_S_Px_S_C1002_ac;
  Double I_ERI_F2xz_Pz_Px_S_C1002_ac = I_ERI_G2x2z_S_Px_S_C1002_ac+ABZ*I_ERI_F2xz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2y_Pz_Px_S_C1002_ac = I_ERI_Gx2yz_S_Px_S_C1002_ac+ABZ*I_ERI_Fx2y_S_Px_S_C1002_ac;
  Double I_ERI_Fxyz_Pz_Px_S_C1002_ac = I_ERI_Gxy2z_S_Px_S_C1002_ac+ABZ*I_ERI_Fxyz_S_Px_S_C1002_ac;
  Double I_ERI_Fx2z_Pz_Px_S_C1002_ac = I_ERI_Gx3z_S_Px_S_C1002_ac+ABZ*I_ERI_Fx2z_S_Px_S_C1002_ac;
  Double I_ERI_F3y_Pz_Px_S_C1002_ac = I_ERI_G3yz_S_Px_S_C1002_ac+ABZ*I_ERI_F3y_S_Px_S_C1002_ac;
  Double I_ERI_F2yz_Pz_Px_S_C1002_ac = I_ERI_G2y2z_S_Px_S_C1002_ac+ABZ*I_ERI_F2yz_S_Px_S_C1002_ac;
  Double I_ERI_Fy2z_Pz_Px_S_C1002_ac = I_ERI_Gy3z_S_Px_S_C1002_ac+ABZ*I_ERI_Fy2z_S_Px_S_C1002_ac;
  Double I_ERI_F3z_Pz_Px_S_C1002_ac = I_ERI_G4z_S_Px_S_C1002_ac+ABZ*I_ERI_F3z_S_Px_S_C1002_ac;
  Double I_ERI_F3x_Px_Py_S_C1002_ac = I_ERI_G4x_S_Py_S_C1002_ac+ABX*I_ERI_F3x_S_Py_S_C1002_ac;
  Double I_ERI_F2xy_Px_Py_S_C1002_ac = I_ERI_G3xy_S_Py_S_C1002_ac+ABX*I_ERI_F2xy_S_Py_S_C1002_ac;
  Double I_ERI_F2xz_Px_Py_S_C1002_ac = I_ERI_G3xz_S_Py_S_C1002_ac+ABX*I_ERI_F2xz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2y_Px_Py_S_C1002_ac = I_ERI_G2x2y_S_Py_S_C1002_ac+ABX*I_ERI_Fx2y_S_Py_S_C1002_ac;
  Double I_ERI_Fxyz_Px_Py_S_C1002_ac = I_ERI_G2xyz_S_Py_S_C1002_ac+ABX*I_ERI_Fxyz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2z_Px_Py_S_C1002_ac = I_ERI_G2x2z_S_Py_S_C1002_ac+ABX*I_ERI_Fx2z_S_Py_S_C1002_ac;
  Double I_ERI_F3y_Px_Py_S_C1002_ac = I_ERI_Gx3y_S_Py_S_C1002_ac+ABX*I_ERI_F3y_S_Py_S_C1002_ac;
  Double I_ERI_F2yz_Px_Py_S_C1002_ac = I_ERI_Gx2yz_S_Py_S_C1002_ac+ABX*I_ERI_F2yz_S_Py_S_C1002_ac;
  Double I_ERI_Fy2z_Px_Py_S_C1002_ac = I_ERI_Gxy2z_S_Py_S_C1002_ac+ABX*I_ERI_Fy2z_S_Py_S_C1002_ac;
  Double I_ERI_F3z_Px_Py_S_C1002_ac = I_ERI_Gx3z_S_Py_S_C1002_ac+ABX*I_ERI_F3z_S_Py_S_C1002_ac;
  Double I_ERI_F3x_Py_Py_S_C1002_ac = I_ERI_G3xy_S_Py_S_C1002_ac+ABY*I_ERI_F3x_S_Py_S_C1002_ac;
  Double I_ERI_F2xy_Py_Py_S_C1002_ac = I_ERI_G2x2y_S_Py_S_C1002_ac+ABY*I_ERI_F2xy_S_Py_S_C1002_ac;
  Double I_ERI_F2xz_Py_Py_S_C1002_ac = I_ERI_G2xyz_S_Py_S_C1002_ac+ABY*I_ERI_F2xz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2y_Py_Py_S_C1002_ac = I_ERI_Gx3y_S_Py_S_C1002_ac+ABY*I_ERI_Fx2y_S_Py_S_C1002_ac;
  Double I_ERI_Fxyz_Py_Py_S_C1002_ac = I_ERI_Gx2yz_S_Py_S_C1002_ac+ABY*I_ERI_Fxyz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2z_Py_Py_S_C1002_ac = I_ERI_Gxy2z_S_Py_S_C1002_ac+ABY*I_ERI_Fx2z_S_Py_S_C1002_ac;
  Double I_ERI_F3y_Py_Py_S_C1002_ac = I_ERI_G4y_S_Py_S_C1002_ac+ABY*I_ERI_F3y_S_Py_S_C1002_ac;
  Double I_ERI_F2yz_Py_Py_S_C1002_ac = I_ERI_G3yz_S_Py_S_C1002_ac+ABY*I_ERI_F2yz_S_Py_S_C1002_ac;
  Double I_ERI_Fy2z_Py_Py_S_C1002_ac = I_ERI_G2y2z_S_Py_S_C1002_ac+ABY*I_ERI_Fy2z_S_Py_S_C1002_ac;
  Double I_ERI_F3z_Py_Py_S_C1002_ac = I_ERI_Gy3z_S_Py_S_C1002_ac+ABY*I_ERI_F3z_S_Py_S_C1002_ac;
  Double I_ERI_F3x_Pz_Py_S_C1002_ac = I_ERI_G3xz_S_Py_S_C1002_ac+ABZ*I_ERI_F3x_S_Py_S_C1002_ac;
  Double I_ERI_F2xy_Pz_Py_S_C1002_ac = I_ERI_G2xyz_S_Py_S_C1002_ac+ABZ*I_ERI_F2xy_S_Py_S_C1002_ac;
  Double I_ERI_F2xz_Pz_Py_S_C1002_ac = I_ERI_G2x2z_S_Py_S_C1002_ac+ABZ*I_ERI_F2xz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2y_Pz_Py_S_C1002_ac = I_ERI_Gx2yz_S_Py_S_C1002_ac+ABZ*I_ERI_Fx2y_S_Py_S_C1002_ac;
  Double I_ERI_Fxyz_Pz_Py_S_C1002_ac = I_ERI_Gxy2z_S_Py_S_C1002_ac+ABZ*I_ERI_Fxyz_S_Py_S_C1002_ac;
  Double I_ERI_Fx2z_Pz_Py_S_C1002_ac = I_ERI_Gx3z_S_Py_S_C1002_ac+ABZ*I_ERI_Fx2z_S_Py_S_C1002_ac;
  Double I_ERI_F3y_Pz_Py_S_C1002_ac = I_ERI_G3yz_S_Py_S_C1002_ac+ABZ*I_ERI_F3y_S_Py_S_C1002_ac;
  Double I_ERI_F2yz_Pz_Py_S_C1002_ac = I_ERI_G2y2z_S_Py_S_C1002_ac+ABZ*I_ERI_F2yz_S_Py_S_C1002_ac;
  Double I_ERI_Fy2z_Pz_Py_S_C1002_ac = I_ERI_Gy3z_S_Py_S_C1002_ac+ABZ*I_ERI_Fy2z_S_Py_S_C1002_ac;
  Double I_ERI_F3z_Pz_Py_S_C1002_ac = I_ERI_G4z_S_Py_S_C1002_ac+ABZ*I_ERI_F3z_S_Py_S_C1002_ac;
  Double I_ERI_F3x_Px_Pz_S_C1002_ac = I_ERI_G4x_S_Pz_S_C1002_ac+ABX*I_ERI_F3x_S_Pz_S_C1002_ac;
  Double I_ERI_F2xy_Px_Pz_S_C1002_ac = I_ERI_G3xy_S_Pz_S_C1002_ac+ABX*I_ERI_F2xy_S_Pz_S_C1002_ac;
  Double I_ERI_F2xz_Px_Pz_S_C1002_ac = I_ERI_G3xz_S_Pz_S_C1002_ac+ABX*I_ERI_F2xz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2y_Px_Pz_S_C1002_ac = I_ERI_G2x2y_S_Pz_S_C1002_ac+ABX*I_ERI_Fx2y_S_Pz_S_C1002_ac;
  Double I_ERI_Fxyz_Px_Pz_S_C1002_ac = I_ERI_G2xyz_S_Pz_S_C1002_ac+ABX*I_ERI_Fxyz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2z_Px_Pz_S_C1002_ac = I_ERI_G2x2z_S_Pz_S_C1002_ac+ABX*I_ERI_Fx2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3y_Px_Pz_S_C1002_ac = I_ERI_Gx3y_S_Pz_S_C1002_ac+ABX*I_ERI_F3y_S_Pz_S_C1002_ac;
  Double I_ERI_F2yz_Px_Pz_S_C1002_ac = I_ERI_Gx2yz_S_Pz_S_C1002_ac+ABX*I_ERI_F2yz_S_Pz_S_C1002_ac;
  Double I_ERI_Fy2z_Px_Pz_S_C1002_ac = I_ERI_Gxy2z_S_Pz_S_C1002_ac+ABX*I_ERI_Fy2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3z_Px_Pz_S_C1002_ac = I_ERI_Gx3z_S_Pz_S_C1002_ac+ABX*I_ERI_F3z_S_Pz_S_C1002_ac;
  Double I_ERI_F3x_Py_Pz_S_C1002_ac = I_ERI_G3xy_S_Pz_S_C1002_ac+ABY*I_ERI_F3x_S_Pz_S_C1002_ac;
  Double I_ERI_F2xy_Py_Pz_S_C1002_ac = I_ERI_G2x2y_S_Pz_S_C1002_ac+ABY*I_ERI_F2xy_S_Pz_S_C1002_ac;
  Double I_ERI_F2xz_Py_Pz_S_C1002_ac = I_ERI_G2xyz_S_Pz_S_C1002_ac+ABY*I_ERI_F2xz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2y_Py_Pz_S_C1002_ac = I_ERI_Gx3y_S_Pz_S_C1002_ac+ABY*I_ERI_Fx2y_S_Pz_S_C1002_ac;
  Double I_ERI_Fxyz_Py_Pz_S_C1002_ac = I_ERI_Gx2yz_S_Pz_S_C1002_ac+ABY*I_ERI_Fxyz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2z_Py_Pz_S_C1002_ac = I_ERI_Gxy2z_S_Pz_S_C1002_ac+ABY*I_ERI_Fx2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3y_Py_Pz_S_C1002_ac = I_ERI_G4y_S_Pz_S_C1002_ac+ABY*I_ERI_F3y_S_Pz_S_C1002_ac;
  Double I_ERI_F2yz_Py_Pz_S_C1002_ac = I_ERI_G3yz_S_Pz_S_C1002_ac+ABY*I_ERI_F2yz_S_Pz_S_C1002_ac;
  Double I_ERI_Fy2z_Py_Pz_S_C1002_ac = I_ERI_G2y2z_S_Pz_S_C1002_ac+ABY*I_ERI_Fy2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3z_Py_Pz_S_C1002_ac = I_ERI_Gy3z_S_Pz_S_C1002_ac+ABY*I_ERI_F3z_S_Pz_S_C1002_ac;
  Double I_ERI_F3x_Pz_Pz_S_C1002_ac = I_ERI_G3xz_S_Pz_S_C1002_ac+ABZ*I_ERI_F3x_S_Pz_S_C1002_ac;
  Double I_ERI_F2xy_Pz_Pz_S_C1002_ac = I_ERI_G2xyz_S_Pz_S_C1002_ac+ABZ*I_ERI_F2xy_S_Pz_S_C1002_ac;
  Double I_ERI_F2xz_Pz_Pz_S_C1002_ac = I_ERI_G2x2z_S_Pz_S_C1002_ac+ABZ*I_ERI_F2xz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2y_Pz_Pz_S_C1002_ac = I_ERI_Gx2yz_S_Pz_S_C1002_ac+ABZ*I_ERI_Fx2y_S_Pz_S_C1002_ac;
  Double I_ERI_Fxyz_Pz_Pz_S_C1002_ac = I_ERI_Gxy2z_S_Pz_S_C1002_ac+ABZ*I_ERI_Fxyz_S_Pz_S_C1002_ac;
  Double I_ERI_Fx2z_Pz_Pz_S_C1002_ac = I_ERI_Gx3z_S_Pz_S_C1002_ac+ABZ*I_ERI_Fx2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3y_Pz_Pz_S_C1002_ac = I_ERI_G3yz_S_Pz_S_C1002_ac+ABZ*I_ERI_F3y_S_Pz_S_C1002_ac;
  Double I_ERI_F2yz_Pz_Pz_S_C1002_ac = I_ERI_G2y2z_S_Pz_S_C1002_ac+ABZ*I_ERI_F2yz_S_Pz_S_C1002_ac;
  Double I_ERI_Fy2z_Pz_Pz_S_C1002_ac = I_ERI_Gy3z_S_Pz_S_C1002_ac+ABZ*I_ERI_Fy2z_S_Pz_S_C1002_ac;
  Double I_ERI_F3z_Pz_Pz_S_C1002_ac = I_ERI_G4z_S_Pz_S_C1002_ac+ABZ*I_ERI_F3z_S_Pz_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C2_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_bb
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C2_bb = I_ERI_F3x_S_S_S_C2_bb+ABX*I_ERI_D2x_S_S_S_C2_bb;
  Double I_ERI_Dxy_Px_S_S_C2_bb = I_ERI_F2xy_S_S_S_C2_bb+ABX*I_ERI_Dxy_S_S_S_C2_bb;
  Double I_ERI_Dxz_Px_S_S_C2_bb = I_ERI_F2xz_S_S_S_C2_bb+ABX*I_ERI_Dxz_S_S_S_C2_bb;
  Double I_ERI_D2y_Px_S_S_C2_bb = I_ERI_Fx2y_S_S_S_C2_bb+ABX*I_ERI_D2y_S_S_S_C2_bb;
  Double I_ERI_Dyz_Px_S_S_C2_bb = I_ERI_Fxyz_S_S_S_C2_bb+ABX*I_ERI_Dyz_S_S_S_C2_bb;
  Double I_ERI_D2z_Px_S_S_C2_bb = I_ERI_Fx2z_S_S_S_C2_bb+ABX*I_ERI_D2z_S_S_S_C2_bb;
  Double I_ERI_D2x_Py_S_S_C2_bb = I_ERI_F2xy_S_S_S_C2_bb+ABY*I_ERI_D2x_S_S_S_C2_bb;
  Double I_ERI_Dxy_Py_S_S_C2_bb = I_ERI_Fx2y_S_S_S_C2_bb+ABY*I_ERI_Dxy_S_S_S_C2_bb;
  Double I_ERI_Dxz_Py_S_S_C2_bb = I_ERI_Fxyz_S_S_S_C2_bb+ABY*I_ERI_Dxz_S_S_S_C2_bb;
  Double I_ERI_D2y_Py_S_S_C2_bb = I_ERI_F3y_S_S_S_C2_bb+ABY*I_ERI_D2y_S_S_S_C2_bb;
  Double I_ERI_Dyz_Py_S_S_C2_bb = I_ERI_F2yz_S_S_S_C2_bb+ABY*I_ERI_Dyz_S_S_S_C2_bb;
  Double I_ERI_D2z_Py_S_S_C2_bb = I_ERI_Fy2z_S_S_S_C2_bb+ABY*I_ERI_D2z_S_S_S_C2_bb;
  Double I_ERI_D2x_Pz_S_S_C2_bb = I_ERI_F2xz_S_S_S_C2_bb+ABZ*I_ERI_D2x_S_S_S_C2_bb;
  Double I_ERI_Dxy_Pz_S_S_C2_bb = I_ERI_Fxyz_S_S_S_C2_bb+ABZ*I_ERI_Dxy_S_S_S_C2_bb;
  Double I_ERI_Dxz_Pz_S_S_C2_bb = I_ERI_Fx2z_S_S_S_C2_bb+ABZ*I_ERI_Dxz_S_S_S_C2_bb;
  Double I_ERI_D2y_Pz_S_S_C2_bb = I_ERI_F2yz_S_S_S_C2_bb+ABZ*I_ERI_D2y_S_S_S_C2_bb;
  Double I_ERI_Dyz_Pz_S_S_C2_bb = I_ERI_Fy2z_S_S_S_C2_bb+ABZ*I_ERI_Dyz_S_S_S_C2_bb;
  Double I_ERI_D2z_Pz_S_S_C2_bb = I_ERI_F3z_S_S_S_C2_bb+ABZ*I_ERI_D2z_S_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C2_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_bb
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C2_bb = I_ERI_G4x_S_S_S_C2_bb+ABX*I_ERI_F3x_S_S_S_C2_bb;
  Double I_ERI_F2xy_Px_S_S_C2_bb = I_ERI_G3xy_S_S_S_C2_bb+ABX*I_ERI_F2xy_S_S_S_C2_bb;
  Double I_ERI_F2xz_Px_S_S_C2_bb = I_ERI_G3xz_S_S_S_C2_bb+ABX*I_ERI_F2xz_S_S_S_C2_bb;
  Double I_ERI_Fx2y_Px_S_S_C2_bb = I_ERI_G2x2y_S_S_S_C2_bb+ABX*I_ERI_Fx2y_S_S_S_C2_bb;
  Double I_ERI_Fxyz_Px_S_S_C2_bb = I_ERI_G2xyz_S_S_S_C2_bb+ABX*I_ERI_Fxyz_S_S_S_C2_bb;
  Double I_ERI_Fx2z_Px_S_S_C2_bb = I_ERI_G2x2z_S_S_S_C2_bb+ABX*I_ERI_Fx2z_S_S_S_C2_bb;
  Double I_ERI_F3y_Px_S_S_C2_bb = I_ERI_Gx3y_S_S_S_C2_bb+ABX*I_ERI_F3y_S_S_S_C2_bb;
  Double I_ERI_F2yz_Px_S_S_C2_bb = I_ERI_Gx2yz_S_S_S_C2_bb+ABX*I_ERI_F2yz_S_S_S_C2_bb;
  Double I_ERI_Fy2z_Px_S_S_C2_bb = I_ERI_Gxy2z_S_S_S_C2_bb+ABX*I_ERI_Fy2z_S_S_S_C2_bb;
  Double I_ERI_F3z_Px_S_S_C2_bb = I_ERI_Gx3z_S_S_S_C2_bb+ABX*I_ERI_F3z_S_S_S_C2_bb;
  Double I_ERI_F2xy_Py_S_S_C2_bb = I_ERI_G2x2y_S_S_S_C2_bb+ABY*I_ERI_F2xy_S_S_S_C2_bb;
  Double I_ERI_F2xz_Py_S_S_C2_bb = I_ERI_G2xyz_S_S_S_C2_bb+ABY*I_ERI_F2xz_S_S_S_C2_bb;
  Double I_ERI_Fx2y_Py_S_S_C2_bb = I_ERI_Gx3y_S_S_S_C2_bb+ABY*I_ERI_Fx2y_S_S_S_C2_bb;
  Double I_ERI_Fxyz_Py_S_S_C2_bb = I_ERI_Gx2yz_S_S_S_C2_bb+ABY*I_ERI_Fxyz_S_S_S_C2_bb;
  Double I_ERI_Fx2z_Py_S_S_C2_bb = I_ERI_Gxy2z_S_S_S_C2_bb+ABY*I_ERI_Fx2z_S_S_S_C2_bb;
  Double I_ERI_F3y_Py_S_S_C2_bb = I_ERI_G4y_S_S_S_C2_bb+ABY*I_ERI_F3y_S_S_S_C2_bb;
  Double I_ERI_F2yz_Py_S_S_C2_bb = I_ERI_G3yz_S_S_S_C2_bb+ABY*I_ERI_F2yz_S_S_S_C2_bb;
  Double I_ERI_Fy2z_Py_S_S_C2_bb = I_ERI_G2y2z_S_S_S_C2_bb+ABY*I_ERI_Fy2z_S_S_S_C2_bb;
  Double I_ERI_F3z_Py_S_S_C2_bb = I_ERI_Gy3z_S_S_S_C2_bb+ABY*I_ERI_F3z_S_S_S_C2_bb;
  Double I_ERI_F2xz_Pz_S_S_C2_bb = I_ERI_G2x2z_S_S_S_C2_bb+ABZ*I_ERI_F2xz_S_S_S_C2_bb;
  Double I_ERI_Fxyz_Pz_S_S_C2_bb = I_ERI_Gxy2z_S_S_S_C2_bb+ABZ*I_ERI_Fxyz_S_S_S_C2_bb;
  Double I_ERI_Fx2z_Pz_S_S_C2_bb = I_ERI_Gx3z_S_S_S_C2_bb+ABZ*I_ERI_Fx2z_S_S_S_C2_bb;
  Double I_ERI_F2yz_Pz_S_S_C2_bb = I_ERI_G2y2z_S_S_S_C2_bb+ABZ*I_ERI_F2yz_S_S_S_C2_bb;
  Double I_ERI_Fy2z_Pz_S_S_C2_bb = I_ERI_Gy3z_S_S_S_C2_bb+ABZ*I_ERI_Fy2z_S_S_S_C2_bb;
  Double I_ERI_F3z_Pz_S_S_C2_bb = I_ERI_G4z_S_S_S_C2_bb+ABZ*I_ERI_F3z_S_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C2_bb
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_C2_bb = I_ERI_F3x_Px_S_S_C2_bb+ABX*I_ERI_D2x_Px_S_S_C2_bb;
  Double I_ERI_Dxy_D2x_S_S_C2_bb = I_ERI_F2xy_Px_S_S_C2_bb+ABX*I_ERI_Dxy_Px_S_S_C2_bb;
  Double I_ERI_Dxz_D2x_S_S_C2_bb = I_ERI_F2xz_Px_S_S_C2_bb+ABX*I_ERI_Dxz_Px_S_S_C2_bb;
  Double I_ERI_D2y_D2x_S_S_C2_bb = I_ERI_Fx2y_Px_S_S_C2_bb+ABX*I_ERI_D2y_Px_S_S_C2_bb;
  Double I_ERI_Dyz_D2x_S_S_C2_bb = I_ERI_Fxyz_Px_S_S_C2_bb+ABX*I_ERI_Dyz_Px_S_S_C2_bb;
  Double I_ERI_D2z_D2x_S_S_C2_bb = I_ERI_Fx2z_Px_S_S_C2_bb+ABX*I_ERI_D2z_Px_S_S_C2_bb;
  Double I_ERI_D2x_Dxy_S_S_C2_bb = I_ERI_F2xy_Px_S_S_C2_bb+ABY*I_ERI_D2x_Px_S_S_C2_bb;
  Double I_ERI_Dxy_Dxy_S_S_C2_bb = I_ERI_Fx2y_Px_S_S_C2_bb+ABY*I_ERI_Dxy_Px_S_S_C2_bb;
  Double I_ERI_Dxz_Dxy_S_S_C2_bb = I_ERI_Fxyz_Px_S_S_C2_bb+ABY*I_ERI_Dxz_Px_S_S_C2_bb;
  Double I_ERI_D2y_Dxy_S_S_C2_bb = I_ERI_F3y_Px_S_S_C2_bb+ABY*I_ERI_D2y_Px_S_S_C2_bb;
  Double I_ERI_Dyz_Dxy_S_S_C2_bb = I_ERI_F2yz_Px_S_S_C2_bb+ABY*I_ERI_Dyz_Px_S_S_C2_bb;
  Double I_ERI_D2z_Dxy_S_S_C2_bb = I_ERI_Fy2z_Px_S_S_C2_bb+ABY*I_ERI_D2z_Px_S_S_C2_bb;
  Double I_ERI_D2x_Dxz_S_S_C2_bb = I_ERI_F2xz_Px_S_S_C2_bb+ABZ*I_ERI_D2x_Px_S_S_C2_bb;
  Double I_ERI_Dxy_Dxz_S_S_C2_bb = I_ERI_Fxyz_Px_S_S_C2_bb+ABZ*I_ERI_Dxy_Px_S_S_C2_bb;
  Double I_ERI_Dxz_Dxz_S_S_C2_bb = I_ERI_Fx2z_Px_S_S_C2_bb+ABZ*I_ERI_Dxz_Px_S_S_C2_bb;
  Double I_ERI_D2y_Dxz_S_S_C2_bb = I_ERI_F2yz_Px_S_S_C2_bb+ABZ*I_ERI_D2y_Px_S_S_C2_bb;
  Double I_ERI_Dyz_Dxz_S_S_C2_bb = I_ERI_Fy2z_Px_S_S_C2_bb+ABZ*I_ERI_Dyz_Px_S_S_C2_bb;
  Double I_ERI_D2z_Dxz_S_S_C2_bb = I_ERI_F3z_Px_S_S_C2_bb+ABZ*I_ERI_D2z_Px_S_S_C2_bb;
  Double I_ERI_D2x_D2y_S_S_C2_bb = I_ERI_F2xy_Py_S_S_C2_bb+ABY*I_ERI_D2x_Py_S_S_C2_bb;
  Double I_ERI_Dxy_D2y_S_S_C2_bb = I_ERI_Fx2y_Py_S_S_C2_bb+ABY*I_ERI_Dxy_Py_S_S_C2_bb;
  Double I_ERI_Dxz_D2y_S_S_C2_bb = I_ERI_Fxyz_Py_S_S_C2_bb+ABY*I_ERI_Dxz_Py_S_S_C2_bb;
  Double I_ERI_D2y_D2y_S_S_C2_bb = I_ERI_F3y_Py_S_S_C2_bb+ABY*I_ERI_D2y_Py_S_S_C2_bb;
  Double I_ERI_Dyz_D2y_S_S_C2_bb = I_ERI_F2yz_Py_S_S_C2_bb+ABY*I_ERI_Dyz_Py_S_S_C2_bb;
  Double I_ERI_D2z_D2y_S_S_C2_bb = I_ERI_Fy2z_Py_S_S_C2_bb+ABY*I_ERI_D2z_Py_S_S_C2_bb;
  Double I_ERI_D2x_Dyz_S_S_C2_bb = I_ERI_F2xz_Py_S_S_C2_bb+ABZ*I_ERI_D2x_Py_S_S_C2_bb;
  Double I_ERI_Dxy_Dyz_S_S_C2_bb = I_ERI_Fxyz_Py_S_S_C2_bb+ABZ*I_ERI_Dxy_Py_S_S_C2_bb;
  Double I_ERI_Dxz_Dyz_S_S_C2_bb = I_ERI_Fx2z_Py_S_S_C2_bb+ABZ*I_ERI_Dxz_Py_S_S_C2_bb;
  Double I_ERI_D2y_Dyz_S_S_C2_bb = I_ERI_F2yz_Py_S_S_C2_bb+ABZ*I_ERI_D2y_Py_S_S_C2_bb;
  Double I_ERI_Dyz_Dyz_S_S_C2_bb = I_ERI_Fy2z_Py_S_S_C2_bb+ABZ*I_ERI_Dyz_Py_S_S_C2_bb;
  Double I_ERI_D2z_Dyz_S_S_C2_bb = I_ERI_F3z_Py_S_S_C2_bb+ABZ*I_ERI_D2z_Py_S_S_C2_bb;
  Double I_ERI_D2x_D2z_S_S_C2_bb = I_ERI_F2xz_Pz_S_S_C2_bb+ABZ*I_ERI_D2x_Pz_S_S_C2_bb;
  Double I_ERI_Dxy_D2z_S_S_C2_bb = I_ERI_Fxyz_Pz_S_S_C2_bb+ABZ*I_ERI_Dxy_Pz_S_S_C2_bb;
  Double I_ERI_Dxz_D2z_S_S_C2_bb = I_ERI_Fx2z_Pz_S_S_C2_bb+ABZ*I_ERI_Dxz_Pz_S_S_C2_bb;
  Double I_ERI_D2y_D2z_S_S_C2_bb = I_ERI_F2yz_Pz_S_S_C2_bb+ABZ*I_ERI_D2y_Pz_S_S_C2_bb;
  Double I_ERI_Dyz_D2z_S_S_C2_bb = I_ERI_Fy2z_Pz_S_S_C2_bb+ABZ*I_ERI_Dyz_Pz_S_S_C2_bb;
  Double I_ERI_D2z_D2z_S_S_C2_bb = I_ERI_F3z_Pz_S_S_C2_bb+ABZ*I_ERI_D2z_Pz_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1002_bb = I_ERI_F3x_S_S_S_C1002_bb+ABX*I_ERI_D2x_S_S_S_C1002_bb;
  Double I_ERI_Dxy_Px_S_S_C1002_bb = I_ERI_F2xy_S_S_S_C1002_bb+ABX*I_ERI_Dxy_S_S_S_C1002_bb;
  Double I_ERI_Dxz_Px_S_S_C1002_bb = I_ERI_F2xz_S_S_S_C1002_bb+ABX*I_ERI_Dxz_S_S_S_C1002_bb;
  Double I_ERI_D2y_Px_S_S_C1002_bb = I_ERI_Fx2y_S_S_S_C1002_bb+ABX*I_ERI_D2y_S_S_S_C1002_bb;
  Double I_ERI_Dyz_Px_S_S_C1002_bb = I_ERI_Fxyz_S_S_S_C1002_bb+ABX*I_ERI_Dyz_S_S_S_C1002_bb;
  Double I_ERI_D2z_Px_S_S_C1002_bb = I_ERI_Fx2z_S_S_S_C1002_bb+ABX*I_ERI_D2z_S_S_S_C1002_bb;
  Double I_ERI_D2x_Py_S_S_C1002_bb = I_ERI_F2xy_S_S_S_C1002_bb+ABY*I_ERI_D2x_S_S_S_C1002_bb;
  Double I_ERI_Dxy_Py_S_S_C1002_bb = I_ERI_Fx2y_S_S_S_C1002_bb+ABY*I_ERI_Dxy_S_S_S_C1002_bb;
  Double I_ERI_Dxz_Py_S_S_C1002_bb = I_ERI_Fxyz_S_S_S_C1002_bb+ABY*I_ERI_Dxz_S_S_S_C1002_bb;
  Double I_ERI_D2y_Py_S_S_C1002_bb = I_ERI_F3y_S_S_S_C1002_bb+ABY*I_ERI_D2y_S_S_S_C1002_bb;
  Double I_ERI_Dyz_Py_S_S_C1002_bb = I_ERI_F2yz_S_S_S_C1002_bb+ABY*I_ERI_Dyz_S_S_S_C1002_bb;
  Double I_ERI_D2z_Py_S_S_C1002_bb = I_ERI_Fy2z_S_S_S_C1002_bb+ABY*I_ERI_D2z_S_S_S_C1002_bb;
  Double I_ERI_D2x_Pz_S_S_C1002_bb = I_ERI_F2xz_S_S_S_C1002_bb+ABZ*I_ERI_D2x_S_S_S_C1002_bb;
  Double I_ERI_Dxy_Pz_S_S_C1002_bb = I_ERI_Fxyz_S_S_S_C1002_bb+ABZ*I_ERI_Dxy_S_S_S_C1002_bb;
  Double I_ERI_Dxz_Pz_S_S_C1002_bb = I_ERI_Fx2z_S_S_S_C1002_bb+ABZ*I_ERI_Dxz_S_S_S_C1002_bb;
  Double I_ERI_D2y_Pz_S_S_C1002_bb = I_ERI_F2yz_S_S_S_C1002_bb+ABZ*I_ERI_D2y_S_S_S_C1002_bb;
  Double I_ERI_Dyz_Pz_S_S_C1002_bb = I_ERI_Fy2z_S_S_S_C1002_bb+ABZ*I_ERI_Dyz_S_S_S_C1002_bb;
  Double I_ERI_D2z_Pz_S_S_C1002_bb = I_ERI_F3z_S_S_S_C1002_bb+ABZ*I_ERI_D2z_S_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1002_bb = I_ERI_G4x_S_S_S_C1002_bb+ABX*I_ERI_F3x_S_S_S_C1002_bb;
  Double I_ERI_F2xy_Px_S_S_C1002_bb = I_ERI_G3xy_S_S_S_C1002_bb+ABX*I_ERI_F2xy_S_S_S_C1002_bb;
  Double I_ERI_F2xz_Px_S_S_C1002_bb = I_ERI_G3xz_S_S_S_C1002_bb+ABX*I_ERI_F2xz_S_S_S_C1002_bb;
  Double I_ERI_Fx2y_Px_S_S_C1002_bb = I_ERI_G2x2y_S_S_S_C1002_bb+ABX*I_ERI_Fx2y_S_S_S_C1002_bb;
  Double I_ERI_Fxyz_Px_S_S_C1002_bb = I_ERI_G2xyz_S_S_S_C1002_bb+ABX*I_ERI_Fxyz_S_S_S_C1002_bb;
  Double I_ERI_Fx2z_Px_S_S_C1002_bb = I_ERI_G2x2z_S_S_S_C1002_bb+ABX*I_ERI_Fx2z_S_S_S_C1002_bb;
  Double I_ERI_F3y_Px_S_S_C1002_bb = I_ERI_Gx3y_S_S_S_C1002_bb+ABX*I_ERI_F3y_S_S_S_C1002_bb;
  Double I_ERI_F2yz_Px_S_S_C1002_bb = I_ERI_Gx2yz_S_S_S_C1002_bb+ABX*I_ERI_F2yz_S_S_S_C1002_bb;
  Double I_ERI_Fy2z_Px_S_S_C1002_bb = I_ERI_Gxy2z_S_S_S_C1002_bb+ABX*I_ERI_Fy2z_S_S_S_C1002_bb;
  Double I_ERI_F3z_Px_S_S_C1002_bb = I_ERI_Gx3z_S_S_S_C1002_bb+ABX*I_ERI_F3z_S_S_S_C1002_bb;
  Double I_ERI_F3x_Py_S_S_C1002_bb = I_ERI_G3xy_S_S_S_C1002_bb+ABY*I_ERI_F3x_S_S_S_C1002_bb;
  Double I_ERI_F2xy_Py_S_S_C1002_bb = I_ERI_G2x2y_S_S_S_C1002_bb+ABY*I_ERI_F2xy_S_S_S_C1002_bb;
  Double I_ERI_F2xz_Py_S_S_C1002_bb = I_ERI_G2xyz_S_S_S_C1002_bb+ABY*I_ERI_F2xz_S_S_S_C1002_bb;
  Double I_ERI_Fx2y_Py_S_S_C1002_bb = I_ERI_Gx3y_S_S_S_C1002_bb+ABY*I_ERI_Fx2y_S_S_S_C1002_bb;
  Double I_ERI_Fxyz_Py_S_S_C1002_bb = I_ERI_Gx2yz_S_S_S_C1002_bb+ABY*I_ERI_Fxyz_S_S_S_C1002_bb;
  Double I_ERI_Fx2z_Py_S_S_C1002_bb = I_ERI_Gxy2z_S_S_S_C1002_bb+ABY*I_ERI_Fx2z_S_S_S_C1002_bb;
  Double I_ERI_F3y_Py_S_S_C1002_bb = I_ERI_G4y_S_S_S_C1002_bb+ABY*I_ERI_F3y_S_S_S_C1002_bb;
  Double I_ERI_F2yz_Py_S_S_C1002_bb = I_ERI_G3yz_S_S_S_C1002_bb+ABY*I_ERI_F2yz_S_S_S_C1002_bb;
  Double I_ERI_Fy2z_Py_S_S_C1002_bb = I_ERI_G2y2z_S_S_S_C1002_bb+ABY*I_ERI_Fy2z_S_S_S_C1002_bb;
  Double I_ERI_F3z_Py_S_S_C1002_bb = I_ERI_Gy3z_S_S_S_C1002_bb+ABY*I_ERI_F3z_S_S_S_C1002_bb;
  Double I_ERI_F3x_Pz_S_S_C1002_bb = I_ERI_G3xz_S_S_S_C1002_bb+ABZ*I_ERI_F3x_S_S_S_C1002_bb;
  Double I_ERI_F2xy_Pz_S_S_C1002_bb = I_ERI_G2xyz_S_S_S_C1002_bb+ABZ*I_ERI_F2xy_S_S_S_C1002_bb;
  Double I_ERI_F2xz_Pz_S_S_C1002_bb = I_ERI_G2x2z_S_S_S_C1002_bb+ABZ*I_ERI_F2xz_S_S_S_C1002_bb;
  Double I_ERI_Fx2y_Pz_S_S_C1002_bb = I_ERI_Gx2yz_S_S_S_C1002_bb+ABZ*I_ERI_Fx2y_S_S_S_C1002_bb;
  Double I_ERI_Fxyz_Pz_S_S_C1002_bb = I_ERI_Gxy2z_S_S_S_C1002_bb+ABZ*I_ERI_Fxyz_S_S_S_C1002_bb;
  Double I_ERI_Fx2z_Pz_S_S_C1002_bb = I_ERI_Gx3z_S_S_S_C1002_bb+ABZ*I_ERI_Fx2z_S_S_S_C1002_bb;
  Double I_ERI_F3y_Pz_S_S_C1002_bb = I_ERI_G3yz_S_S_S_C1002_bb+ABZ*I_ERI_F3y_S_S_S_C1002_bb;
  Double I_ERI_F2yz_Pz_S_S_C1002_bb = I_ERI_G2y2z_S_S_S_C1002_bb+ABZ*I_ERI_F2yz_S_S_S_C1002_bb;
  Double I_ERI_Fy2z_Pz_S_S_C1002_bb = I_ERI_Gy3z_S_S_S_C1002_bb+ABZ*I_ERI_Fy2z_S_S_S_C1002_bb;
  Double I_ERI_F3z_Pz_S_S_C1002_bb = I_ERI_G4z_S_S_S_C1002_bb+ABZ*I_ERI_F3z_S_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_C1002_bb = I_ERI_F3x_Px_S_S_C1002_bb+ABX*I_ERI_D2x_Px_S_S_C1002_bb;
  Double I_ERI_Dxy_D2x_S_S_C1002_bb = I_ERI_F2xy_Px_S_S_C1002_bb+ABX*I_ERI_Dxy_Px_S_S_C1002_bb;
  Double I_ERI_Dxz_D2x_S_S_C1002_bb = I_ERI_F2xz_Px_S_S_C1002_bb+ABX*I_ERI_Dxz_Px_S_S_C1002_bb;
  Double I_ERI_D2y_D2x_S_S_C1002_bb = I_ERI_Fx2y_Px_S_S_C1002_bb+ABX*I_ERI_D2y_Px_S_S_C1002_bb;
  Double I_ERI_Dyz_D2x_S_S_C1002_bb = I_ERI_Fxyz_Px_S_S_C1002_bb+ABX*I_ERI_Dyz_Px_S_S_C1002_bb;
  Double I_ERI_D2z_D2x_S_S_C1002_bb = I_ERI_Fx2z_Px_S_S_C1002_bb+ABX*I_ERI_D2z_Px_S_S_C1002_bb;
  Double I_ERI_D2x_Dxy_S_S_C1002_bb = I_ERI_F2xy_Px_S_S_C1002_bb+ABY*I_ERI_D2x_Px_S_S_C1002_bb;
  Double I_ERI_Dxy_Dxy_S_S_C1002_bb = I_ERI_Fx2y_Px_S_S_C1002_bb+ABY*I_ERI_Dxy_Px_S_S_C1002_bb;
  Double I_ERI_Dxz_Dxy_S_S_C1002_bb = I_ERI_Fxyz_Px_S_S_C1002_bb+ABY*I_ERI_Dxz_Px_S_S_C1002_bb;
  Double I_ERI_D2y_Dxy_S_S_C1002_bb = I_ERI_F3y_Px_S_S_C1002_bb+ABY*I_ERI_D2y_Px_S_S_C1002_bb;
  Double I_ERI_Dyz_Dxy_S_S_C1002_bb = I_ERI_F2yz_Px_S_S_C1002_bb+ABY*I_ERI_Dyz_Px_S_S_C1002_bb;
  Double I_ERI_D2z_Dxy_S_S_C1002_bb = I_ERI_Fy2z_Px_S_S_C1002_bb+ABY*I_ERI_D2z_Px_S_S_C1002_bb;
  Double I_ERI_D2x_D2y_S_S_C1002_bb = I_ERI_F2xy_Py_S_S_C1002_bb+ABY*I_ERI_D2x_Py_S_S_C1002_bb;
  Double I_ERI_Dxy_D2y_S_S_C1002_bb = I_ERI_Fx2y_Py_S_S_C1002_bb+ABY*I_ERI_Dxy_Py_S_S_C1002_bb;
  Double I_ERI_Dxz_D2y_S_S_C1002_bb = I_ERI_Fxyz_Py_S_S_C1002_bb+ABY*I_ERI_Dxz_Py_S_S_C1002_bb;
  Double I_ERI_D2y_D2y_S_S_C1002_bb = I_ERI_F3y_Py_S_S_C1002_bb+ABY*I_ERI_D2y_Py_S_S_C1002_bb;
  Double I_ERI_Dyz_D2y_S_S_C1002_bb = I_ERI_F2yz_Py_S_S_C1002_bb+ABY*I_ERI_Dyz_Py_S_S_C1002_bb;
  Double I_ERI_D2z_D2y_S_S_C1002_bb = I_ERI_Fy2z_Py_S_S_C1002_bb+ABY*I_ERI_D2z_Py_S_S_C1002_bb;
  Double I_ERI_D2x_D2z_S_S_C1002_bb = I_ERI_F2xz_Pz_S_S_C1002_bb+ABZ*I_ERI_D2x_Pz_S_S_C1002_bb;
  Double I_ERI_Dxy_D2z_S_S_C1002_bb = I_ERI_Fxyz_Pz_S_S_C1002_bb+ABZ*I_ERI_Dxy_Pz_S_S_C1002_bb;
  Double I_ERI_Dxz_D2z_S_S_C1002_bb = I_ERI_Fx2z_Pz_S_S_C1002_bb+ABZ*I_ERI_Dxz_Pz_S_S_C1002_bb;
  Double I_ERI_D2y_D2z_S_S_C1002_bb = I_ERI_F2yz_Pz_S_S_C1002_bb+ABZ*I_ERI_D2y_Pz_S_S_C1002_bb;
  Double I_ERI_Dyz_D2z_S_S_C1002_bb = I_ERI_Fy2z_Pz_S_S_C1002_bb+ABZ*I_ERI_Dyz_Pz_S_S_C1002_bb;
  Double I_ERI_D2z_D2z_S_S_C1002_bb = I_ERI_F3z_Pz_S_S_C1002_bb+ABZ*I_ERI_D2z_Pz_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_C1002_bb = I_ERI_H5x_S_S_S_C1002_bb+ABX*I_ERI_G4x_S_S_S_C1002_bb;
  Double I_ERI_G3xy_Px_S_S_C1002_bb = I_ERI_H4xy_S_S_S_C1002_bb+ABX*I_ERI_G3xy_S_S_S_C1002_bb;
  Double I_ERI_G3xz_Px_S_S_C1002_bb = I_ERI_H4xz_S_S_S_C1002_bb+ABX*I_ERI_G3xz_S_S_S_C1002_bb;
  Double I_ERI_G2x2y_Px_S_S_C1002_bb = I_ERI_H3x2y_S_S_S_C1002_bb+ABX*I_ERI_G2x2y_S_S_S_C1002_bb;
  Double I_ERI_G2xyz_Px_S_S_C1002_bb = I_ERI_H3xyz_S_S_S_C1002_bb+ABX*I_ERI_G2xyz_S_S_S_C1002_bb;
  Double I_ERI_G2x2z_Px_S_S_C1002_bb = I_ERI_H3x2z_S_S_S_C1002_bb+ABX*I_ERI_G2x2z_S_S_S_C1002_bb;
  Double I_ERI_Gx3y_Px_S_S_C1002_bb = I_ERI_H2x3y_S_S_S_C1002_bb+ABX*I_ERI_Gx3y_S_S_S_C1002_bb;
  Double I_ERI_Gx2yz_Px_S_S_C1002_bb = I_ERI_H2x2yz_S_S_S_C1002_bb+ABX*I_ERI_Gx2yz_S_S_S_C1002_bb;
  Double I_ERI_Gxy2z_Px_S_S_C1002_bb = I_ERI_H2xy2z_S_S_S_C1002_bb+ABX*I_ERI_Gxy2z_S_S_S_C1002_bb;
  Double I_ERI_Gx3z_Px_S_S_C1002_bb = I_ERI_H2x3z_S_S_S_C1002_bb+ABX*I_ERI_Gx3z_S_S_S_C1002_bb;
  Double I_ERI_G3yz_Px_S_S_C1002_bb = I_ERI_Hx3yz_S_S_S_C1002_bb+ABX*I_ERI_G3yz_S_S_S_C1002_bb;
  Double I_ERI_G2y2z_Px_S_S_C1002_bb = I_ERI_Hx2y2z_S_S_S_C1002_bb+ABX*I_ERI_G2y2z_S_S_S_C1002_bb;
  Double I_ERI_Gy3z_Px_S_S_C1002_bb = I_ERI_Hxy3z_S_S_S_C1002_bb+ABX*I_ERI_Gy3z_S_S_S_C1002_bb;
  Double I_ERI_G3xy_Py_S_S_C1002_bb = I_ERI_H3x2y_S_S_S_C1002_bb+ABY*I_ERI_G3xy_S_S_S_C1002_bb;
  Double I_ERI_G2x2y_Py_S_S_C1002_bb = I_ERI_H2x3y_S_S_S_C1002_bb+ABY*I_ERI_G2x2y_S_S_S_C1002_bb;
  Double I_ERI_G2xyz_Py_S_S_C1002_bb = I_ERI_H2x2yz_S_S_S_C1002_bb+ABY*I_ERI_G2xyz_S_S_S_C1002_bb;
  Double I_ERI_Gx3y_Py_S_S_C1002_bb = I_ERI_Hx4y_S_S_S_C1002_bb+ABY*I_ERI_Gx3y_S_S_S_C1002_bb;
  Double I_ERI_Gx2yz_Py_S_S_C1002_bb = I_ERI_Hx3yz_S_S_S_C1002_bb+ABY*I_ERI_Gx2yz_S_S_S_C1002_bb;
  Double I_ERI_Gxy2z_Py_S_S_C1002_bb = I_ERI_Hx2y2z_S_S_S_C1002_bb+ABY*I_ERI_Gxy2z_S_S_S_C1002_bb;
  Double I_ERI_G4y_Py_S_S_C1002_bb = I_ERI_H5y_S_S_S_C1002_bb+ABY*I_ERI_G4y_S_S_S_C1002_bb;
  Double I_ERI_G3yz_Py_S_S_C1002_bb = I_ERI_H4yz_S_S_S_C1002_bb+ABY*I_ERI_G3yz_S_S_S_C1002_bb;
  Double I_ERI_G2y2z_Py_S_S_C1002_bb = I_ERI_H3y2z_S_S_S_C1002_bb+ABY*I_ERI_G2y2z_S_S_S_C1002_bb;
  Double I_ERI_Gy3z_Py_S_S_C1002_bb = I_ERI_H2y3z_S_S_S_C1002_bb+ABY*I_ERI_Gy3z_S_S_S_C1002_bb;
  Double I_ERI_G3xz_Pz_S_S_C1002_bb = I_ERI_H3x2z_S_S_S_C1002_bb+ABZ*I_ERI_G3xz_S_S_S_C1002_bb;
  Double I_ERI_G2xyz_Pz_S_S_C1002_bb = I_ERI_H2xy2z_S_S_S_C1002_bb+ABZ*I_ERI_G2xyz_S_S_S_C1002_bb;
  Double I_ERI_G2x2z_Pz_S_S_C1002_bb = I_ERI_H2x3z_S_S_S_C1002_bb+ABZ*I_ERI_G2x2z_S_S_S_C1002_bb;
  Double I_ERI_Gx2yz_Pz_S_S_C1002_bb = I_ERI_Hx2y2z_S_S_S_C1002_bb+ABZ*I_ERI_Gx2yz_S_S_S_C1002_bb;
  Double I_ERI_Gxy2z_Pz_S_S_C1002_bb = I_ERI_Hxy3z_S_S_S_C1002_bb+ABZ*I_ERI_Gxy2z_S_S_S_C1002_bb;
  Double I_ERI_Gx3z_Pz_S_S_C1002_bb = I_ERI_Hx4z_S_S_S_C1002_bb+ABZ*I_ERI_Gx3z_S_S_S_C1002_bb;
  Double I_ERI_G3yz_Pz_S_S_C1002_bb = I_ERI_H3y2z_S_S_S_C1002_bb+ABZ*I_ERI_G3yz_S_S_S_C1002_bb;
  Double I_ERI_G2y2z_Pz_S_S_C1002_bb = I_ERI_H2y3z_S_S_S_C1002_bb+ABZ*I_ERI_G2y2z_S_S_S_C1002_bb;
  Double I_ERI_Gy3z_Pz_S_S_C1002_bb = I_ERI_Hy4z_S_S_S_C1002_bb+ABZ*I_ERI_Gy3z_S_S_S_C1002_bb;
  Double I_ERI_G4z_Pz_S_S_C1002_bb = I_ERI_H5z_S_S_S_C1002_bb+ABZ*I_ERI_G4z_S_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 24 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_C1002_bb = I_ERI_G4x_Px_S_S_C1002_bb+ABX*I_ERI_F3x_Px_S_S_C1002_bb;
  Double I_ERI_F2xy_D2x_S_S_C1002_bb = I_ERI_G3xy_Px_S_S_C1002_bb+ABX*I_ERI_F2xy_Px_S_S_C1002_bb;
  Double I_ERI_F2xz_D2x_S_S_C1002_bb = I_ERI_G3xz_Px_S_S_C1002_bb+ABX*I_ERI_F2xz_Px_S_S_C1002_bb;
  Double I_ERI_Fx2y_D2x_S_S_C1002_bb = I_ERI_G2x2y_Px_S_S_C1002_bb+ABX*I_ERI_Fx2y_Px_S_S_C1002_bb;
  Double I_ERI_Fxyz_D2x_S_S_C1002_bb = I_ERI_G2xyz_Px_S_S_C1002_bb+ABX*I_ERI_Fxyz_Px_S_S_C1002_bb;
  Double I_ERI_Fx2z_D2x_S_S_C1002_bb = I_ERI_G2x2z_Px_S_S_C1002_bb+ABX*I_ERI_Fx2z_Px_S_S_C1002_bb;
  Double I_ERI_F3y_D2x_S_S_C1002_bb = I_ERI_Gx3y_Px_S_S_C1002_bb+ABX*I_ERI_F3y_Px_S_S_C1002_bb;
  Double I_ERI_F2yz_D2x_S_S_C1002_bb = I_ERI_Gx2yz_Px_S_S_C1002_bb+ABX*I_ERI_F2yz_Px_S_S_C1002_bb;
  Double I_ERI_Fy2z_D2x_S_S_C1002_bb = I_ERI_Gxy2z_Px_S_S_C1002_bb+ABX*I_ERI_Fy2z_Px_S_S_C1002_bb;
  Double I_ERI_F3z_D2x_S_S_C1002_bb = I_ERI_Gx3z_Px_S_S_C1002_bb+ABX*I_ERI_F3z_Px_S_S_C1002_bb;
  Double I_ERI_F2xz_Dxy_S_S_C1002_bb = I_ERI_G2xyz_Px_S_S_C1002_bb+ABY*I_ERI_F2xz_Px_S_S_C1002_bb;
  Double I_ERI_Fxyz_Dxy_S_S_C1002_bb = I_ERI_Gx2yz_Px_S_S_C1002_bb+ABY*I_ERI_Fxyz_Px_S_S_C1002_bb;
  Double I_ERI_Fx2z_Dxy_S_S_C1002_bb = I_ERI_Gxy2z_Px_S_S_C1002_bb+ABY*I_ERI_Fx2z_Px_S_S_C1002_bb;
  Double I_ERI_F2yz_Dxy_S_S_C1002_bb = I_ERI_G3yz_Px_S_S_C1002_bb+ABY*I_ERI_F2yz_Px_S_S_C1002_bb;
  Double I_ERI_Fy2z_Dxy_S_S_C1002_bb = I_ERI_G2y2z_Px_S_S_C1002_bb+ABY*I_ERI_Fy2z_Px_S_S_C1002_bb;
  Double I_ERI_F3z_Dxy_S_S_C1002_bb = I_ERI_Gy3z_Px_S_S_C1002_bb+ABY*I_ERI_F3z_Px_S_S_C1002_bb;
  Double I_ERI_F3x_D2y_S_S_C1002_bb = I_ERI_G3xy_Py_S_S_C1002_bb+ABY*I_ERI_F3x_Py_S_S_C1002_bb;
  Double I_ERI_F2xy_D2y_S_S_C1002_bb = I_ERI_G2x2y_Py_S_S_C1002_bb+ABY*I_ERI_F2xy_Py_S_S_C1002_bb;
  Double I_ERI_F2xz_D2y_S_S_C1002_bb = I_ERI_G2xyz_Py_S_S_C1002_bb+ABY*I_ERI_F2xz_Py_S_S_C1002_bb;
  Double I_ERI_Fx2y_D2y_S_S_C1002_bb = I_ERI_Gx3y_Py_S_S_C1002_bb+ABY*I_ERI_Fx2y_Py_S_S_C1002_bb;
  Double I_ERI_Fxyz_D2y_S_S_C1002_bb = I_ERI_Gx2yz_Py_S_S_C1002_bb+ABY*I_ERI_Fxyz_Py_S_S_C1002_bb;
  Double I_ERI_Fx2z_D2y_S_S_C1002_bb = I_ERI_Gxy2z_Py_S_S_C1002_bb+ABY*I_ERI_Fx2z_Py_S_S_C1002_bb;
  Double I_ERI_F3y_D2y_S_S_C1002_bb = I_ERI_G4y_Py_S_S_C1002_bb+ABY*I_ERI_F3y_Py_S_S_C1002_bb;
  Double I_ERI_F2yz_D2y_S_S_C1002_bb = I_ERI_G3yz_Py_S_S_C1002_bb+ABY*I_ERI_F2yz_Py_S_S_C1002_bb;
  Double I_ERI_Fy2z_D2y_S_S_C1002_bb = I_ERI_G2y2z_Py_S_S_C1002_bb+ABY*I_ERI_Fy2z_Py_S_S_C1002_bb;
  Double I_ERI_F3z_D2y_S_S_C1002_bb = I_ERI_Gy3z_Py_S_S_C1002_bb+ABY*I_ERI_F3z_Py_S_S_C1002_bb;
  Double I_ERI_F3x_D2z_S_S_C1002_bb = I_ERI_G3xz_Pz_S_S_C1002_bb+ABZ*I_ERI_F3x_Pz_S_S_C1002_bb;
  Double I_ERI_F2xy_D2z_S_S_C1002_bb = I_ERI_G2xyz_Pz_S_S_C1002_bb+ABZ*I_ERI_F2xy_Pz_S_S_C1002_bb;
  Double I_ERI_F2xz_D2z_S_S_C1002_bb = I_ERI_G2x2z_Pz_S_S_C1002_bb+ABZ*I_ERI_F2xz_Pz_S_S_C1002_bb;
  Double I_ERI_Fx2y_D2z_S_S_C1002_bb = I_ERI_Gx2yz_Pz_S_S_C1002_bb+ABZ*I_ERI_Fx2y_Pz_S_S_C1002_bb;
  Double I_ERI_Fxyz_D2z_S_S_C1002_bb = I_ERI_Gxy2z_Pz_S_S_C1002_bb+ABZ*I_ERI_Fxyz_Pz_S_S_C1002_bb;
  Double I_ERI_Fx2z_D2z_S_S_C1002_bb = I_ERI_Gx3z_Pz_S_S_C1002_bb+ABZ*I_ERI_Fx2z_Pz_S_S_C1002_bb;
  Double I_ERI_F3y_D2z_S_S_C1002_bb = I_ERI_G3yz_Pz_S_S_C1002_bb+ABZ*I_ERI_F3y_Pz_S_S_C1002_bb;
  Double I_ERI_F2yz_D2z_S_S_C1002_bb = I_ERI_G2y2z_Pz_S_S_C1002_bb+ABZ*I_ERI_F2yz_Pz_S_S_C1002_bb;
  Double I_ERI_Fy2z_D2z_S_S_C1002_bb = I_ERI_Gy3z_Pz_S_S_C1002_bb+ABZ*I_ERI_Fy2z_Pz_S_S_C1002_bb;
  Double I_ERI_F3z_D2z_S_S_C1002_bb = I_ERI_G4z_Pz_S_S_C1002_bb+ABZ*I_ERI_F3z_Pz_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1002_bb
   ************************************************************/
  Double I_ERI_D2x_F3x_S_S_C1002_bb = I_ERI_F3x_D2x_S_S_C1002_bb+ABX*I_ERI_D2x_D2x_S_S_C1002_bb;
  Double I_ERI_Dxy_F3x_S_S_C1002_bb = I_ERI_F2xy_D2x_S_S_C1002_bb+ABX*I_ERI_Dxy_D2x_S_S_C1002_bb;
  Double I_ERI_Dxz_F3x_S_S_C1002_bb = I_ERI_F2xz_D2x_S_S_C1002_bb+ABX*I_ERI_Dxz_D2x_S_S_C1002_bb;
  Double I_ERI_D2y_F3x_S_S_C1002_bb = I_ERI_Fx2y_D2x_S_S_C1002_bb+ABX*I_ERI_D2y_D2x_S_S_C1002_bb;
  Double I_ERI_Dyz_F3x_S_S_C1002_bb = I_ERI_Fxyz_D2x_S_S_C1002_bb+ABX*I_ERI_Dyz_D2x_S_S_C1002_bb;
  Double I_ERI_D2z_F3x_S_S_C1002_bb = I_ERI_Fx2z_D2x_S_S_C1002_bb+ABX*I_ERI_D2z_D2x_S_S_C1002_bb;
  Double I_ERI_D2x_F2xy_S_S_C1002_bb = I_ERI_F2xy_D2x_S_S_C1002_bb+ABY*I_ERI_D2x_D2x_S_S_C1002_bb;
  Double I_ERI_Dxy_F2xy_S_S_C1002_bb = I_ERI_Fx2y_D2x_S_S_C1002_bb+ABY*I_ERI_Dxy_D2x_S_S_C1002_bb;
  Double I_ERI_Dxz_F2xy_S_S_C1002_bb = I_ERI_Fxyz_D2x_S_S_C1002_bb+ABY*I_ERI_Dxz_D2x_S_S_C1002_bb;
  Double I_ERI_D2y_F2xy_S_S_C1002_bb = I_ERI_F3y_D2x_S_S_C1002_bb+ABY*I_ERI_D2y_D2x_S_S_C1002_bb;
  Double I_ERI_Dyz_F2xy_S_S_C1002_bb = I_ERI_F2yz_D2x_S_S_C1002_bb+ABY*I_ERI_Dyz_D2x_S_S_C1002_bb;
  Double I_ERI_D2z_F2xy_S_S_C1002_bb = I_ERI_Fy2z_D2x_S_S_C1002_bb+ABY*I_ERI_D2z_D2x_S_S_C1002_bb;
  Double I_ERI_D2x_F2xz_S_S_C1002_bb = I_ERI_F2xz_D2x_S_S_C1002_bb+ABZ*I_ERI_D2x_D2x_S_S_C1002_bb;
  Double I_ERI_Dxy_F2xz_S_S_C1002_bb = I_ERI_Fxyz_D2x_S_S_C1002_bb+ABZ*I_ERI_Dxy_D2x_S_S_C1002_bb;
  Double I_ERI_Dxz_F2xz_S_S_C1002_bb = I_ERI_Fx2z_D2x_S_S_C1002_bb+ABZ*I_ERI_Dxz_D2x_S_S_C1002_bb;
  Double I_ERI_D2y_F2xz_S_S_C1002_bb = I_ERI_F2yz_D2x_S_S_C1002_bb+ABZ*I_ERI_D2y_D2x_S_S_C1002_bb;
  Double I_ERI_Dyz_F2xz_S_S_C1002_bb = I_ERI_Fy2z_D2x_S_S_C1002_bb+ABZ*I_ERI_Dyz_D2x_S_S_C1002_bb;
  Double I_ERI_D2z_F2xz_S_S_C1002_bb = I_ERI_F3z_D2x_S_S_C1002_bb+ABZ*I_ERI_D2z_D2x_S_S_C1002_bb;
  Double I_ERI_D2x_Fx2y_S_S_C1002_bb = I_ERI_F3x_D2y_S_S_C1002_bb+ABX*I_ERI_D2x_D2y_S_S_C1002_bb;
  Double I_ERI_Dxy_Fx2y_S_S_C1002_bb = I_ERI_F2xy_D2y_S_S_C1002_bb+ABX*I_ERI_Dxy_D2y_S_S_C1002_bb;
  Double I_ERI_Dxz_Fx2y_S_S_C1002_bb = I_ERI_F2xz_D2y_S_S_C1002_bb+ABX*I_ERI_Dxz_D2y_S_S_C1002_bb;
  Double I_ERI_D2y_Fx2y_S_S_C1002_bb = I_ERI_Fx2y_D2y_S_S_C1002_bb+ABX*I_ERI_D2y_D2y_S_S_C1002_bb;
  Double I_ERI_Dyz_Fx2y_S_S_C1002_bb = I_ERI_Fxyz_D2y_S_S_C1002_bb+ABX*I_ERI_Dyz_D2y_S_S_C1002_bb;
  Double I_ERI_D2z_Fx2y_S_S_C1002_bb = I_ERI_Fx2z_D2y_S_S_C1002_bb+ABX*I_ERI_D2z_D2y_S_S_C1002_bb;
  Double I_ERI_D2x_Fxyz_S_S_C1002_bb = I_ERI_F2xz_Dxy_S_S_C1002_bb+ABZ*I_ERI_D2x_Dxy_S_S_C1002_bb;
  Double I_ERI_Dxy_Fxyz_S_S_C1002_bb = I_ERI_Fxyz_Dxy_S_S_C1002_bb+ABZ*I_ERI_Dxy_Dxy_S_S_C1002_bb;
  Double I_ERI_Dxz_Fxyz_S_S_C1002_bb = I_ERI_Fx2z_Dxy_S_S_C1002_bb+ABZ*I_ERI_Dxz_Dxy_S_S_C1002_bb;
  Double I_ERI_D2y_Fxyz_S_S_C1002_bb = I_ERI_F2yz_Dxy_S_S_C1002_bb+ABZ*I_ERI_D2y_Dxy_S_S_C1002_bb;
  Double I_ERI_Dyz_Fxyz_S_S_C1002_bb = I_ERI_Fy2z_Dxy_S_S_C1002_bb+ABZ*I_ERI_Dyz_Dxy_S_S_C1002_bb;
  Double I_ERI_D2z_Fxyz_S_S_C1002_bb = I_ERI_F3z_Dxy_S_S_C1002_bb+ABZ*I_ERI_D2z_Dxy_S_S_C1002_bb;
  Double I_ERI_D2x_Fx2z_S_S_C1002_bb = I_ERI_F3x_D2z_S_S_C1002_bb+ABX*I_ERI_D2x_D2z_S_S_C1002_bb;
  Double I_ERI_Dxy_Fx2z_S_S_C1002_bb = I_ERI_F2xy_D2z_S_S_C1002_bb+ABX*I_ERI_Dxy_D2z_S_S_C1002_bb;
  Double I_ERI_Dxz_Fx2z_S_S_C1002_bb = I_ERI_F2xz_D2z_S_S_C1002_bb+ABX*I_ERI_Dxz_D2z_S_S_C1002_bb;
  Double I_ERI_D2y_Fx2z_S_S_C1002_bb = I_ERI_Fx2y_D2z_S_S_C1002_bb+ABX*I_ERI_D2y_D2z_S_S_C1002_bb;
  Double I_ERI_Dyz_Fx2z_S_S_C1002_bb = I_ERI_Fxyz_D2z_S_S_C1002_bb+ABX*I_ERI_Dyz_D2z_S_S_C1002_bb;
  Double I_ERI_D2z_Fx2z_S_S_C1002_bb = I_ERI_Fx2z_D2z_S_S_C1002_bb+ABX*I_ERI_D2z_D2z_S_S_C1002_bb;
  Double I_ERI_D2x_F3y_S_S_C1002_bb = I_ERI_F2xy_D2y_S_S_C1002_bb+ABY*I_ERI_D2x_D2y_S_S_C1002_bb;
  Double I_ERI_Dxy_F3y_S_S_C1002_bb = I_ERI_Fx2y_D2y_S_S_C1002_bb+ABY*I_ERI_Dxy_D2y_S_S_C1002_bb;
  Double I_ERI_Dxz_F3y_S_S_C1002_bb = I_ERI_Fxyz_D2y_S_S_C1002_bb+ABY*I_ERI_Dxz_D2y_S_S_C1002_bb;
  Double I_ERI_D2y_F3y_S_S_C1002_bb = I_ERI_F3y_D2y_S_S_C1002_bb+ABY*I_ERI_D2y_D2y_S_S_C1002_bb;
  Double I_ERI_Dyz_F3y_S_S_C1002_bb = I_ERI_F2yz_D2y_S_S_C1002_bb+ABY*I_ERI_Dyz_D2y_S_S_C1002_bb;
  Double I_ERI_D2z_F3y_S_S_C1002_bb = I_ERI_Fy2z_D2y_S_S_C1002_bb+ABY*I_ERI_D2z_D2y_S_S_C1002_bb;
  Double I_ERI_D2x_F2yz_S_S_C1002_bb = I_ERI_F2xz_D2y_S_S_C1002_bb+ABZ*I_ERI_D2x_D2y_S_S_C1002_bb;
  Double I_ERI_Dxy_F2yz_S_S_C1002_bb = I_ERI_Fxyz_D2y_S_S_C1002_bb+ABZ*I_ERI_Dxy_D2y_S_S_C1002_bb;
  Double I_ERI_Dxz_F2yz_S_S_C1002_bb = I_ERI_Fx2z_D2y_S_S_C1002_bb+ABZ*I_ERI_Dxz_D2y_S_S_C1002_bb;
  Double I_ERI_D2y_F2yz_S_S_C1002_bb = I_ERI_F2yz_D2y_S_S_C1002_bb+ABZ*I_ERI_D2y_D2y_S_S_C1002_bb;
  Double I_ERI_Dyz_F2yz_S_S_C1002_bb = I_ERI_Fy2z_D2y_S_S_C1002_bb+ABZ*I_ERI_Dyz_D2y_S_S_C1002_bb;
  Double I_ERI_D2z_F2yz_S_S_C1002_bb = I_ERI_F3z_D2y_S_S_C1002_bb+ABZ*I_ERI_D2z_D2y_S_S_C1002_bb;
  Double I_ERI_D2x_Fy2z_S_S_C1002_bb = I_ERI_F2xy_D2z_S_S_C1002_bb+ABY*I_ERI_D2x_D2z_S_S_C1002_bb;
  Double I_ERI_Dxy_Fy2z_S_S_C1002_bb = I_ERI_Fx2y_D2z_S_S_C1002_bb+ABY*I_ERI_Dxy_D2z_S_S_C1002_bb;
  Double I_ERI_Dxz_Fy2z_S_S_C1002_bb = I_ERI_Fxyz_D2z_S_S_C1002_bb+ABY*I_ERI_Dxz_D2z_S_S_C1002_bb;
  Double I_ERI_D2y_Fy2z_S_S_C1002_bb = I_ERI_F3y_D2z_S_S_C1002_bb+ABY*I_ERI_D2y_D2z_S_S_C1002_bb;
  Double I_ERI_Dyz_Fy2z_S_S_C1002_bb = I_ERI_F2yz_D2z_S_S_C1002_bb+ABY*I_ERI_Dyz_D2z_S_S_C1002_bb;
  Double I_ERI_D2z_Fy2z_S_S_C1002_bb = I_ERI_Fy2z_D2z_S_S_C1002_bb+ABY*I_ERI_D2z_D2z_S_S_C1002_bb;
  Double I_ERI_D2x_F3z_S_S_C1002_bb = I_ERI_F2xz_D2z_S_S_C1002_bb+ABZ*I_ERI_D2x_D2z_S_S_C1002_bb;
  Double I_ERI_Dxy_F3z_S_S_C1002_bb = I_ERI_Fxyz_D2z_S_S_C1002_bb+ABZ*I_ERI_Dxy_D2z_S_S_C1002_bb;
  Double I_ERI_Dxz_F3z_S_S_C1002_bb = I_ERI_Fx2z_D2z_S_S_C1002_bb+ABZ*I_ERI_Dxz_D2z_S_S_C1002_bb;
  Double I_ERI_D2y_F3z_S_S_C1002_bb = I_ERI_F2yz_D2z_S_S_C1002_bb+ABZ*I_ERI_D2y_D2z_S_S_C1002_bb;
  Double I_ERI_Dyz_F3z_S_S_C1002_bb = I_ERI_Fy2z_D2z_S_S_C1002_bb+ABZ*I_ERI_Dyz_D2z_S_S_C1002_bb;
  Double I_ERI_D2z_F3z_S_S_C1002_bb = I_ERI_F3z_D2z_S_S_C1002_bb+ABZ*I_ERI_D2z_D2z_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C2_bc
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C2_bc = I_ERI_F3x_S_Px_S_C2_bc+ABX*I_ERI_D2x_S_Px_S_C2_bc;
  Double I_ERI_Dxy_Px_Px_S_C2_bc = I_ERI_F2xy_S_Px_S_C2_bc+ABX*I_ERI_Dxy_S_Px_S_C2_bc;
  Double I_ERI_Dxz_Px_Px_S_C2_bc = I_ERI_F2xz_S_Px_S_C2_bc+ABX*I_ERI_Dxz_S_Px_S_C2_bc;
  Double I_ERI_D2y_Px_Px_S_C2_bc = I_ERI_Fx2y_S_Px_S_C2_bc+ABX*I_ERI_D2y_S_Px_S_C2_bc;
  Double I_ERI_Dyz_Px_Px_S_C2_bc = I_ERI_Fxyz_S_Px_S_C2_bc+ABX*I_ERI_Dyz_S_Px_S_C2_bc;
  Double I_ERI_D2z_Px_Px_S_C2_bc = I_ERI_Fx2z_S_Px_S_C2_bc+ABX*I_ERI_D2z_S_Px_S_C2_bc;
  Double I_ERI_D2x_Py_Px_S_C2_bc = I_ERI_F2xy_S_Px_S_C2_bc+ABY*I_ERI_D2x_S_Px_S_C2_bc;
  Double I_ERI_Dxy_Py_Px_S_C2_bc = I_ERI_Fx2y_S_Px_S_C2_bc+ABY*I_ERI_Dxy_S_Px_S_C2_bc;
  Double I_ERI_Dxz_Py_Px_S_C2_bc = I_ERI_Fxyz_S_Px_S_C2_bc+ABY*I_ERI_Dxz_S_Px_S_C2_bc;
  Double I_ERI_D2y_Py_Px_S_C2_bc = I_ERI_F3y_S_Px_S_C2_bc+ABY*I_ERI_D2y_S_Px_S_C2_bc;
  Double I_ERI_Dyz_Py_Px_S_C2_bc = I_ERI_F2yz_S_Px_S_C2_bc+ABY*I_ERI_Dyz_S_Px_S_C2_bc;
  Double I_ERI_D2z_Py_Px_S_C2_bc = I_ERI_Fy2z_S_Px_S_C2_bc+ABY*I_ERI_D2z_S_Px_S_C2_bc;
  Double I_ERI_D2x_Pz_Px_S_C2_bc = I_ERI_F2xz_S_Px_S_C2_bc+ABZ*I_ERI_D2x_S_Px_S_C2_bc;
  Double I_ERI_Dxy_Pz_Px_S_C2_bc = I_ERI_Fxyz_S_Px_S_C2_bc+ABZ*I_ERI_Dxy_S_Px_S_C2_bc;
  Double I_ERI_Dxz_Pz_Px_S_C2_bc = I_ERI_Fx2z_S_Px_S_C2_bc+ABZ*I_ERI_Dxz_S_Px_S_C2_bc;
  Double I_ERI_D2y_Pz_Px_S_C2_bc = I_ERI_F2yz_S_Px_S_C2_bc+ABZ*I_ERI_D2y_S_Px_S_C2_bc;
  Double I_ERI_Dyz_Pz_Px_S_C2_bc = I_ERI_Fy2z_S_Px_S_C2_bc+ABZ*I_ERI_Dyz_S_Px_S_C2_bc;
  Double I_ERI_D2z_Pz_Px_S_C2_bc = I_ERI_F3z_S_Px_S_C2_bc+ABZ*I_ERI_D2z_S_Px_S_C2_bc;
  Double I_ERI_D2x_Px_Py_S_C2_bc = I_ERI_F3x_S_Py_S_C2_bc+ABX*I_ERI_D2x_S_Py_S_C2_bc;
  Double I_ERI_Dxy_Px_Py_S_C2_bc = I_ERI_F2xy_S_Py_S_C2_bc+ABX*I_ERI_Dxy_S_Py_S_C2_bc;
  Double I_ERI_Dxz_Px_Py_S_C2_bc = I_ERI_F2xz_S_Py_S_C2_bc+ABX*I_ERI_Dxz_S_Py_S_C2_bc;
  Double I_ERI_D2y_Px_Py_S_C2_bc = I_ERI_Fx2y_S_Py_S_C2_bc+ABX*I_ERI_D2y_S_Py_S_C2_bc;
  Double I_ERI_Dyz_Px_Py_S_C2_bc = I_ERI_Fxyz_S_Py_S_C2_bc+ABX*I_ERI_Dyz_S_Py_S_C2_bc;
  Double I_ERI_D2z_Px_Py_S_C2_bc = I_ERI_Fx2z_S_Py_S_C2_bc+ABX*I_ERI_D2z_S_Py_S_C2_bc;
  Double I_ERI_D2x_Py_Py_S_C2_bc = I_ERI_F2xy_S_Py_S_C2_bc+ABY*I_ERI_D2x_S_Py_S_C2_bc;
  Double I_ERI_Dxy_Py_Py_S_C2_bc = I_ERI_Fx2y_S_Py_S_C2_bc+ABY*I_ERI_Dxy_S_Py_S_C2_bc;
  Double I_ERI_Dxz_Py_Py_S_C2_bc = I_ERI_Fxyz_S_Py_S_C2_bc+ABY*I_ERI_Dxz_S_Py_S_C2_bc;
  Double I_ERI_D2y_Py_Py_S_C2_bc = I_ERI_F3y_S_Py_S_C2_bc+ABY*I_ERI_D2y_S_Py_S_C2_bc;
  Double I_ERI_Dyz_Py_Py_S_C2_bc = I_ERI_F2yz_S_Py_S_C2_bc+ABY*I_ERI_Dyz_S_Py_S_C2_bc;
  Double I_ERI_D2z_Py_Py_S_C2_bc = I_ERI_Fy2z_S_Py_S_C2_bc+ABY*I_ERI_D2z_S_Py_S_C2_bc;
  Double I_ERI_D2x_Pz_Py_S_C2_bc = I_ERI_F2xz_S_Py_S_C2_bc+ABZ*I_ERI_D2x_S_Py_S_C2_bc;
  Double I_ERI_Dxy_Pz_Py_S_C2_bc = I_ERI_Fxyz_S_Py_S_C2_bc+ABZ*I_ERI_Dxy_S_Py_S_C2_bc;
  Double I_ERI_Dxz_Pz_Py_S_C2_bc = I_ERI_Fx2z_S_Py_S_C2_bc+ABZ*I_ERI_Dxz_S_Py_S_C2_bc;
  Double I_ERI_D2y_Pz_Py_S_C2_bc = I_ERI_F2yz_S_Py_S_C2_bc+ABZ*I_ERI_D2y_S_Py_S_C2_bc;
  Double I_ERI_Dyz_Pz_Py_S_C2_bc = I_ERI_Fy2z_S_Py_S_C2_bc+ABZ*I_ERI_Dyz_S_Py_S_C2_bc;
  Double I_ERI_D2z_Pz_Py_S_C2_bc = I_ERI_F3z_S_Py_S_C2_bc+ABZ*I_ERI_D2z_S_Py_S_C2_bc;
  Double I_ERI_D2x_Px_Pz_S_C2_bc = I_ERI_F3x_S_Pz_S_C2_bc+ABX*I_ERI_D2x_S_Pz_S_C2_bc;
  Double I_ERI_Dxy_Px_Pz_S_C2_bc = I_ERI_F2xy_S_Pz_S_C2_bc+ABX*I_ERI_Dxy_S_Pz_S_C2_bc;
  Double I_ERI_Dxz_Px_Pz_S_C2_bc = I_ERI_F2xz_S_Pz_S_C2_bc+ABX*I_ERI_Dxz_S_Pz_S_C2_bc;
  Double I_ERI_D2y_Px_Pz_S_C2_bc = I_ERI_Fx2y_S_Pz_S_C2_bc+ABX*I_ERI_D2y_S_Pz_S_C2_bc;
  Double I_ERI_Dyz_Px_Pz_S_C2_bc = I_ERI_Fxyz_S_Pz_S_C2_bc+ABX*I_ERI_Dyz_S_Pz_S_C2_bc;
  Double I_ERI_D2z_Px_Pz_S_C2_bc = I_ERI_Fx2z_S_Pz_S_C2_bc+ABX*I_ERI_D2z_S_Pz_S_C2_bc;
  Double I_ERI_D2x_Py_Pz_S_C2_bc = I_ERI_F2xy_S_Pz_S_C2_bc+ABY*I_ERI_D2x_S_Pz_S_C2_bc;
  Double I_ERI_Dxy_Py_Pz_S_C2_bc = I_ERI_Fx2y_S_Pz_S_C2_bc+ABY*I_ERI_Dxy_S_Pz_S_C2_bc;
  Double I_ERI_Dxz_Py_Pz_S_C2_bc = I_ERI_Fxyz_S_Pz_S_C2_bc+ABY*I_ERI_Dxz_S_Pz_S_C2_bc;
  Double I_ERI_D2y_Py_Pz_S_C2_bc = I_ERI_F3y_S_Pz_S_C2_bc+ABY*I_ERI_D2y_S_Pz_S_C2_bc;
  Double I_ERI_Dyz_Py_Pz_S_C2_bc = I_ERI_F2yz_S_Pz_S_C2_bc+ABY*I_ERI_Dyz_S_Pz_S_C2_bc;
  Double I_ERI_D2z_Py_Pz_S_C2_bc = I_ERI_Fy2z_S_Pz_S_C2_bc+ABY*I_ERI_D2z_S_Pz_S_C2_bc;
  Double I_ERI_D2x_Pz_Pz_S_C2_bc = I_ERI_F2xz_S_Pz_S_C2_bc+ABZ*I_ERI_D2x_S_Pz_S_C2_bc;
  Double I_ERI_Dxy_Pz_Pz_S_C2_bc = I_ERI_Fxyz_S_Pz_S_C2_bc+ABZ*I_ERI_Dxy_S_Pz_S_C2_bc;
  Double I_ERI_Dxz_Pz_Pz_S_C2_bc = I_ERI_Fx2z_S_Pz_S_C2_bc+ABZ*I_ERI_Dxz_S_Pz_S_C2_bc;
  Double I_ERI_D2y_Pz_Pz_S_C2_bc = I_ERI_F2yz_S_Pz_S_C2_bc+ABZ*I_ERI_D2y_S_Pz_S_C2_bc;
  Double I_ERI_Dyz_Pz_Pz_S_C2_bc = I_ERI_Fy2z_S_Pz_S_C2_bc+ABZ*I_ERI_Dyz_S_Pz_S_C2_bc;
  Double I_ERI_D2z_Pz_Pz_S_C2_bc = I_ERI_F3z_S_Pz_S_C2_bc+ABZ*I_ERI_D2z_S_Pz_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1002_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_bc
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1002_bc = I_ERI_F3x_S_Px_S_C1002_bc+ABX*I_ERI_D2x_S_Px_S_C1002_bc;
  Double I_ERI_Dxy_Px_Px_S_C1002_bc = I_ERI_F2xy_S_Px_S_C1002_bc+ABX*I_ERI_Dxy_S_Px_S_C1002_bc;
  Double I_ERI_Dxz_Px_Px_S_C1002_bc = I_ERI_F2xz_S_Px_S_C1002_bc+ABX*I_ERI_Dxz_S_Px_S_C1002_bc;
  Double I_ERI_D2y_Px_Px_S_C1002_bc = I_ERI_Fx2y_S_Px_S_C1002_bc+ABX*I_ERI_D2y_S_Px_S_C1002_bc;
  Double I_ERI_Dyz_Px_Px_S_C1002_bc = I_ERI_Fxyz_S_Px_S_C1002_bc+ABX*I_ERI_Dyz_S_Px_S_C1002_bc;
  Double I_ERI_D2z_Px_Px_S_C1002_bc = I_ERI_Fx2z_S_Px_S_C1002_bc+ABX*I_ERI_D2z_S_Px_S_C1002_bc;
  Double I_ERI_D2x_Py_Px_S_C1002_bc = I_ERI_F2xy_S_Px_S_C1002_bc+ABY*I_ERI_D2x_S_Px_S_C1002_bc;
  Double I_ERI_Dxy_Py_Px_S_C1002_bc = I_ERI_Fx2y_S_Px_S_C1002_bc+ABY*I_ERI_Dxy_S_Px_S_C1002_bc;
  Double I_ERI_Dxz_Py_Px_S_C1002_bc = I_ERI_Fxyz_S_Px_S_C1002_bc+ABY*I_ERI_Dxz_S_Px_S_C1002_bc;
  Double I_ERI_D2y_Py_Px_S_C1002_bc = I_ERI_F3y_S_Px_S_C1002_bc+ABY*I_ERI_D2y_S_Px_S_C1002_bc;
  Double I_ERI_Dyz_Py_Px_S_C1002_bc = I_ERI_F2yz_S_Px_S_C1002_bc+ABY*I_ERI_Dyz_S_Px_S_C1002_bc;
  Double I_ERI_D2z_Py_Px_S_C1002_bc = I_ERI_Fy2z_S_Px_S_C1002_bc+ABY*I_ERI_D2z_S_Px_S_C1002_bc;
  Double I_ERI_D2x_Pz_Px_S_C1002_bc = I_ERI_F2xz_S_Px_S_C1002_bc+ABZ*I_ERI_D2x_S_Px_S_C1002_bc;
  Double I_ERI_Dxy_Pz_Px_S_C1002_bc = I_ERI_Fxyz_S_Px_S_C1002_bc+ABZ*I_ERI_Dxy_S_Px_S_C1002_bc;
  Double I_ERI_Dxz_Pz_Px_S_C1002_bc = I_ERI_Fx2z_S_Px_S_C1002_bc+ABZ*I_ERI_Dxz_S_Px_S_C1002_bc;
  Double I_ERI_D2y_Pz_Px_S_C1002_bc = I_ERI_F2yz_S_Px_S_C1002_bc+ABZ*I_ERI_D2y_S_Px_S_C1002_bc;
  Double I_ERI_Dyz_Pz_Px_S_C1002_bc = I_ERI_Fy2z_S_Px_S_C1002_bc+ABZ*I_ERI_Dyz_S_Px_S_C1002_bc;
  Double I_ERI_D2z_Pz_Px_S_C1002_bc = I_ERI_F3z_S_Px_S_C1002_bc+ABZ*I_ERI_D2z_S_Px_S_C1002_bc;
  Double I_ERI_D2x_Px_Py_S_C1002_bc = I_ERI_F3x_S_Py_S_C1002_bc+ABX*I_ERI_D2x_S_Py_S_C1002_bc;
  Double I_ERI_Dxy_Px_Py_S_C1002_bc = I_ERI_F2xy_S_Py_S_C1002_bc+ABX*I_ERI_Dxy_S_Py_S_C1002_bc;
  Double I_ERI_Dxz_Px_Py_S_C1002_bc = I_ERI_F2xz_S_Py_S_C1002_bc+ABX*I_ERI_Dxz_S_Py_S_C1002_bc;
  Double I_ERI_D2y_Px_Py_S_C1002_bc = I_ERI_Fx2y_S_Py_S_C1002_bc+ABX*I_ERI_D2y_S_Py_S_C1002_bc;
  Double I_ERI_Dyz_Px_Py_S_C1002_bc = I_ERI_Fxyz_S_Py_S_C1002_bc+ABX*I_ERI_Dyz_S_Py_S_C1002_bc;
  Double I_ERI_D2z_Px_Py_S_C1002_bc = I_ERI_Fx2z_S_Py_S_C1002_bc+ABX*I_ERI_D2z_S_Py_S_C1002_bc;
  Double I_ERI_D2x_Py_Py_S_C1002_bc = I_ERI_F2xy_S_Py_S_C1002_bc+ABY*I_ERI_D2x_S_Py_S_C1002_bc;
  Double I_ERI_Dxy_Py_Py_S_C1002_bc = I_ERI_Fx2y_S_Py_S_C1002_bc+ABY*I_ERI_Dxy_S_Py_S_C1002_bc;
  Double I_ERI_Dxz_Py_Py_S_C1002_bc = I_ERI_Fxyz_S_Py_S_C1002_bc+ABY*I_ERI_Dxz_S_Py_S_C1002_bc;
  Double I_ERI_D2y_Py_Py_S_C1002_bc = I_ERI_F3y_S_Py_S_C1002_bc+ABY*I_ERI_D2y_S_Py_S_C1002_bc;
  Double I_ERI_Dyz_Py_Py_S_C1002_bc = I_ERI_F2yz_S_Py_S_C1002_bc+ABY*I_ERI_Dyz_S_Py_S_C1002_bc;
  Double I_ERI_D2z_Py_Py_S_C1002_bc = I_ERI_Fy2z_S_Py_S_C1002_bc+ABY*I_ERI_D2z_S_Py_S_C1002_bc;
  Double I_ERI_D2x_Pz_Py_S_C1002_bc = I_ERI_F2xz_S_Py_S_C1002_bc+ABZ*I_ERI_D2x_S_Py_S_C1002_bc;
  Double I_ERI_Dxy_Pz_Py_S_C1002_bc = I_ERI_Fxyz_S_Py_S_C1002_bc+ABZ*I_ERI_Dxy_S_Py_S_C1002_bc;
  Double I_ERI_Dxz_Pz_Py_S_C1002_bc = I_ERI_Fx2z_S_Py_S_C1002_bc+ABZ*I_ERI_Dxz_S_Py_S_C1002_bc;
  Double I_ERI_D2y_Pz_Py_S_C1002_bc = I_ERI_F2yz_S_Py_S_C1002_bc+ABZ*I_ERI_D2y_S_Py_S_C1002_bc;
  Double I_ERI_Dyz_Pz_Py_S_C1002_bc = I_ERI_Fy2z_S_Py_S_C1002_bc+ABZ*I_ERI_Dyz_S_Py_S_C1002_bc;
  Double I_ERI_D2z_Pz_Py_S_C1002_bc = I_ERI_F3z_S_Py_S_C1002_bc+ABZ*I_ERI_D2z_S_Py_S_C1002_bc;
  Double I_ERI_D2x_Px_Pz_S_C1002_bc = I_ERI_F3x_S_Pz_S_C1002_bc+ABX*I_ERI_D2x_S_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Px_Pz_S_C1002_bc = I_ERI_F2xy_S_Pz_S_C1002_bc+ABX*I_ERI_Dxy_S_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Px_Pz_S_C1002_bc = I_ERI_F2xz_S_Pz_S_C1002_bc+ABX*I_ERI_Dxz_S_Pz_S_C1002_bc;
  Double I_ERI_D2y_Px_Pz_S_C1002_bc = I_ERI_Fx2y_S_Pz_S_C1002_bc+ABX*I_ERI_D2y_S_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Px_Pz_S_C1002_bc = I_ERI_Fxyz_S_Pz_S_C1002_bc+ABX*I_ERI_Dyz_S_Pz_S_C1002_bc;
  Double I_ERI_D2z_Px_Pz_S_C1002_bc = I_ERI_Fx2z_S_Pz_S_C1002_bc+ABX*I_ERI_D2z_S_Pz_S_C1002_bc;
  Double I_ERI_D2x_Py_Pz_S_C1002_bc = I_ERI_F2xy_S_Pz_S_C1002_bc+ABY*I_ERI_D2x_S_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Py_Pz_S_C1002_bc = I_ERI_Fx2y_S_Pz_S_C1002_bc+ABY*I_ERI_Dxy_S_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Py_Pz_S_C1002_bc = I_ERI_Fxyz_S_Pz_S_C1002_bc+ABY*I_ERI_Dxz_S_Pz_S_C1002_bc;
  Double I_ERI_D2y_Py_Pz_S_C1002_bc = I_ERI_F3y_S_Pz_S_C1002_bc+ABY*I_ERI_D2y_S_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Py_Pz_S_C1002_bc = I_ERI_F2yz_S_Pz_S_C1002_bc+ABY*I_ERI_Dyz_S_Pz_S_C1002_bc;
  Double I_ERI_D2z_Py_Pz_S_C1002_bc = I_ERI_Fy2z_S_Pz_S_C1002_bc+ABY*I_ERI_D2z_S_Pz_S_C1002_bc;
  Double I_ERI_D2x_Pz_Pz_S_C1002_bc = I_ERI_F2xz_S_Pz_S_C1002_bc+ABZ*I_ERI_D2x_S_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Pz_Pz_S_C1002_bc = I_ERI_Fxyz_S_Pz_S_C1002_bc+ABZ*I_ERI_Dxy_S_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Pz_Pz_S_C1002_bc = I_ERI_Fx2z_S_Pz_S_C1002_bc+ABZ*I_ERI_Dxz_S_Pz_S_C1002_bc;
  Double I_ERI_D2y_Pz_Pz_S_C1002_bc = I_ERI_F2yz_S_Pz_S_C1002_bc+ABZ*I_ERI_D2y_S_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Pz_Pz_S_C1002_bc = I_ERI_Fy2z_S_Pz_S_C1002_bc+ABZ*I_ERI_Dyz_S_Pz_S_C1002_bc;
  Double I_ERI_D2z_Pz_Pz_S_C1002_bc = I_ERI_F3z_S_Pz_S_C1002_bc+ABZ*I_ERI_D2z_S_Pz_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_C1002_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1002_bc
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1002_bc = I_ERI_G4x_S_Px_S_C1002_bc+ABX*I_ERI_F3x_S_Px_S_C1002_bc;
  Double I_ERI_F2xy_Px_Px_S_C1002_bc = I_ERI_G3xy_S_Px_S_C1002_bc+ABX*I_ERI_F2xy_S_Px_S_C1002_bc;
  Double I_ERI_F2xz_Px_Px_S_C1002_bc = I_ERI_G3xz_S_Px_S_C1002_bc+ABX*I_ERI_F2xz_S_Px_S_C1002_bc;
  Double I_ERI_Fx2y_Px_Px_S_C1002_bc = I_ERI_G2x2y_S_Px_S_C1002_bc+ABX*I_ERI_Fx2y_S_Px_S_C1002_bc;
  Double I_ERI_Fxyz_Px_Px_S_C1002_bc = I_ERI_G2xyz_S_Px_S_C1002_bc+ABX*I_ERI_Fxyz_S_Px_S_C1002_bc;
  Double I_ERI_Fx2z_Px_Px_S_C1002_bc = I_ERI_G2x2z_S_Px_S_C1002_bc+ABX*I_ERI_Fx2z_S_Px_S_C1002_bc;
  Double I_ERI_F3y_Px_Px_S_C1002_bc = I_ERI_Gx3y_S_Px_S_C1002_bc+ABX*I_ERI_F3y_S_Px_S_C1002_bc;
  Double I_ERI_F2yz_Px_Px_S_C1002_bc = I_ERI_Gx2yz_S_Px_S_C1002_bc+ABX*I_ERI_F2yz_S_Px_S_C1002_bc;
  Double I_ERI_Fy2z_Px_Px_S_C1002_bc = I_ERI_Gxy2z_S_Px_S_C1002_bc+ABX*I_ERI_Fy2z_S_Px_S_C1002_bc;
  Double I_ERI_F3z_Px_Px_S_C1002_bc = I_ERI_Gx3z_S_Px_S_C1002_bc+ABX*I_ERI_F3z_S_Px_S_C1002_bc;
  Double I_ERI_F2xy_Py_Px_S_C1002_bc = I_ERI_G2x2y_S_Px_S_C1002_bc+ABY*I_ERI_F2xy_S_Px_S_C1002_bc;
  Double I_ERI_F2xz_Py_Px_S_C1002_bc = I_ERI_G2xyz_S_Px_S_C1002_bc+ABY*I_ERI_F2xz_S_Px_S_C1002_bc;
  Double I_ERI_Fx2y_Py_Px_S_C1002_bc = I_ERI_Gx3y_S_Px_S_C1002_bc+ABY*I_ERI_Fx2y_S_Px_S_C1002_bc;
  Double I_ERI_Fxyz_Py_Px_S_C1002_bc = I_ERI_Gx2yz_S_Px_S_C1002_bc+ABY*I_ERI_Fxyz_S_Px_S_C1002_bc;
  Double I_ERI_Fx2z_Py_Px_S_C1002_bc = I_ERI_Gxy2z_S_Px_S_C1002_bc+ABY*I_ERI_Fx2z_S_Px_S_C1002_bc;
  Double I_ERI_F3y_Py_Px_S_C1002_bc = I_ERI_G4y_S_Px_S_C1002_bc+ABY*I_ERI_F3y_S_Px_S_C1002_bc;
  Double I_ERI_F2yz_Py_Px_S_C1002_bc = I_ERI_G3yz_S_Px_S_C1002_bc+ABY*I_ERI_F2yz_S_Px_S_C1002_bc;
  Double I_ERI_Fy2z_Py_Px_S_C1002_bc = I_ERI_G2y2z_S_Px_S_C1002_bc+ABY*I_ERI_Fy2z_S_Px_S_C1002_bc;
  Double I_ERI_F3z_Py_Px_S_C1002_bc = I_ERI_Gy3z_S_Px_S_C1002_bc+ABY*I_ERI_F3z_S_Px_S_C1002_bc;
  Double I_ERI_F2xz_Pz_Px_S_C1002_bc = I_ERI_G2x2z_S_Px_S_C1002_bc+ABZ*I_ERI_F2xz_S_Px_S_C1002_bc;
  Double I_ERI_Fxyz_Pz_Px_S_C1002_bc = I_ERI_Gxy2z_S_Px_S_C1002_bc+ABZ*I_ERI_Fxyz_S_Px_S_C1002_bc;
  Double I_ERI_Fx2z_Pz_Px_S_C1002_bc = I_ERI_Gx3z_S_Px_S_C1002_bc+ABZ*I_ERI_Fx2z_S_Px_S_C1002_bc;
  Double I_ERI_F2yz_Pz_Px_S_C1002_bc = I_ERI_G2y2z_S_Px_S_C1002_bc+ABZ*I_ERI_F2yz_S_Px_S_C1002_bc;
  Double I_ERI_Fy2z_Pz_Px_S_C1002_bc = I_ERI_Gy3z_S_Px_S_C1002_bc+ABZ*I_ERI_Fy2z_S_Px_S_C1002_bc;
  Double I_ERI_F3z_Pz_Px_S_C1002_bc = I_ERI_G4z_S_Px_S_C1002_bc+ABZ*I_ERI_F3z_S_Px_S_C1002_bc;
  Double I_ERI_F3x_Px_Py_S_C1002_bc = I_ERI_G4x_S_Py_S_C1002_bc+ABX*I_ERI_F3x_S_Py_S_C1002_bc;
  Double I_ERI_F2xy_Px_Py_S_C1002_bc = I_ERI_G3xy_S_Py_S_C1002_bc+ABX*I_ERI_F2xy_S_Py_S_C1002_bc;
  Double I_ERI_F2xz_Px_Py_S_C1002_bc = I_ERI_G3xz_S_Py_S_C1002_bc+ABX*I_ERI_F2xz_S_Py_S_C1002_bc;
  Double I_ERI_Fx2y_Px_Py_S_C1002_bc = I_ERI_G2x2y_S_Py_S_C1002_bc+ABX*I_ERI_Fx2y_S_Py_S_C1002_bc;
  Double I_ERI_Fxyz_Px_Py_S_C1002_bc = I_ERI_G2xyz_S_Py_S_C1002_bc+ABX*I_ERI_Fxyz_S_Py_S_C1002_bc;
  Double I_ERI_Fx2z_Px_Py_S_C1002_bc = I_ERI_G2x2z_S_Py_S_C1002_bc+ABX*I_ERI_Fx2z_S_Py_S_C1002_bc;
  Double I_ERI_F3y_Px_Py_S_C1002_bc = I_ERI_Gx3y_S_Py_S_C1002_bc+ABX*I_ERI_F3y_S_Py_S_C1002_bc;
  Double I_ERI_F2yz_Px_Py_S_C1002_bc = I_ERI_Gx2yz_S_Py_S_C1002_bc+ABX*I_ERI_F2yz_S_Py_S_C1002_bc;
  Double I_ERI_Fy2z_Px_Py_S_C1002_bc = I_ERI_Gxy2z_S_Py_S_C1002_bc+ABX*I_ERI_Fy2z_S_Py_S_C1002_bc;
  Double I_ERI_F3z_Px_Py_S_C1002_bc = I_ERI_Gx3z_S_Py_S_C1002_bc+ABX*I_ERI_F3z_S_Py_S_C1002_bc;
  Double I_ERI_F2xy_Py_Py_S_C1002_bc = I_ERI_G2x2y_S_Py_S_C1002_bc+ABY*I_ERI_F2xy_S_Py_S_C1002_bc;
  Double I_ERI_F2xz_Py_Py_S_C1002_bc = I_ERI_G2xyz_S_Py_S_C1002_bc+ABY*I_ERI_F2xz_S_Py_S_C1002_bc;
  Double I_ERI_Fx2y_Py_Py_S_C1002_bc = I_ERI_Gx3y_S_Py_S_C1002_bc+ABY*I_ERI_Fx2y_S_Py_S_C1002_bc;
  Double I_ERI_Fxyz_Py_Py_S_C1002_bc = I_ERI_Gx2yz_S_Py_S_C1002_bc+ABY*I_ERI_Fxyz_S_Py_S_C1002_bc;
  Double I_ERI_Fx2z_Py_Py_S_C1002_bc = I_ERI_Gxy2z_S_Py_S_C1002_bc+ABY*I_ERI_Fx2z_S_Py_S_C1002_bc;
  Double I_ERI_F3y_Py_Py_S_C1002_bc = I_ERI_G4y_S_Py_S_C1002_bc+ABY*I_ERI_F3y_S_Py_S_C1002_bc;
  Double I_ERI_F2yz_Py_Py_S_C1002_bc = I_ERI_G3yz_S_Py_S_C1002_bc+ABY*I_ERI_F2yz_S_Py_S_C1002_bc;
  Double I_ERI_Fy2z_Py_Py_S_C1002_bc = I_ERI_G2y2z_S_Py_S_C1002_bc+ABY*I_ERI_Fy2z_S_Py_S_C1002_bc;
  Double I_ERI_F3z_Py_Py_S_C1002_bc = I_ERI_Gy3z_S_Py_S_C1002_bc+ABY*I_ERI_F3z_S_Py_S_C1002_bc;
  Double I_ERI_F2xz_Pz_Py_S_C1002_bc = I_ERI_G2x2z_S_Py_S_C1002_bc+ABZ*I_ERI_F2xz_S_Py_S_C1002_bc;
  Double I_ERI_Fxyz_Pz_Py_S_C1002_bc = I_ERI_Gxy2z_S_Py_S_C1002_bc+ABZ*I_ERI_Fxyz_S_Py_S_C1002_bc;
  Double I_ERI_Fx2z_Pz_Py_S_C1002_bc = I_ERI_Gx3z_S_Py_S_C1002_bc+ABZ*I_ERI_Fx2z_S_Py_S_C1002_bc;
  Double I_ERI_F2yz_Pz_Py_S_C1002_bc = I_ERI_G2y2z_S_Py_S_C1002_bc+ABZ*I_ERI_F2yz_S_Py_S_C1002_bc;
  Double I_ERI_Fy2z_Pz_Py_S_C1002_bc = I_ERI_Gy3z_S_Py_S_C1002_bc+ABZ*I_ERI_Fy2z_S_Py_S_C1002_bc;
  Double I_ERI_F3z_Pz_Py_S_C1002_bc = I_ERI_G4z_S_Py_S_C1002_bc+ABZ*I_ERI_F3z_S_Py_S_C1002_bc;
  Double I_ERI_F3x_Px_Pz_S_C1002_bc = I_ERI_G4x_S_Pz_S_C1002_bc+ABX*I_ERI_F3x_S_Pz_S_C1002_bc;
  Double I_ERI_F2xy_Px_Pz_S_C1002_bc = I_ERI_G3xy_S_Pz_S_C1002_bc+ABX*I_ERI_F2xy_S_Pz_S_C1002_bc;
  Double I_ERI_F2xz_Px_Pz_S_C1002_bc = I_ERI_G3xz_S_Pz_S_C1002_bc+ABX*I_ERI_F2xz_S_Pz_S_C1002_bc;
  Double I_ERI_Fx2y_Px_Pz_S_C1002_bc = I_ERI_G2x2y_S_Pz_S_C1002_bc+ABX*I_ERI_Fx2y_S_Pz_S_C1002_bc;
  Double I_ERI_Fxyz_Px_Pz_S_C1002_bc = I_ERI_G2xyz_S_Pz_S_C1002_bc+ABX*I_ERI_Fxyz_S_Pz_S_C1002_bc;
  Double I_ERI_Fx2z_Px_Pz_S_C1002_bc = I_ERI_G2x2z_S_Pz_S_C1002_bc+ABX*I_ERI_Fx2z_S_Pz_S_C1002_bc;
  Double I_ERI_F3y_Px_Pz_S_C1002_bc = I_ERI_Gx3y_S_Pz_S_C1002_bc+ABX*I_ERI_F3y_S_Pz_S_C1002_bc;
  Double I_ERI_F2yz_Px_Pz_S_C1002_bc = I_ERI_Gx2yz_S_Pz_S_C1002_bc+ABX*I_ERI_F2yz_S_Pz_S_C1002_bc;
  Double I_ERI_Fy2z_Px_Pz_S_C1002_bc = I_ERI_Gxy2z_S_Pz_S_C1002_bc+ABX*I_ERI_Fy2z_S_Pz_S_C1002_bc;
  Double I_ERI_F3z_Px_Pz_S_C1002_bc = I_ERI_Gx3z_S_Pz_S_C1002_bc+ABX*I_ERI_F3z_S_Pz_S_C1002_bc;
  Double I_ERI_F2xy_Py_Pz_S_C1002_bc = I_ERI_G2x2y_S_Pz_S_C1002_bc+ABY*I_ERI_F2xy_S_Pz_S_C1002_bc;
  Double I_ERI_F2xz_Py_Pz_S_C1002_bc = I_ERI_G2xyz_S_Pz_S_C1002_bc+ABY*I_ERI_F2xz_S_Pz_S_C1002_bc;
  Double I_ERI_Fx2y_Py_Pz_S_C1002_bc = I_ERI_Gx3y_S_Pz_S_C1002_bc+ABY*I_ERI_Fx2y_S_Pz_S_C1002_bc;
  Double I_ERI_Fxyz_Py_Pz_S_C1002_bc = I_ERI_Gx2yz_S_Pz_S_C1002_bc+ABY*I_ERI_Fxyz_S_Pz_S_C1002_bc;
  Double I_ERI_Fx2z_Py_Pz_S_C1002_bc = I_ERI_Gxy2z_S_Pz_S_C1002_bc+ABY*I_ERI_Fx2z_S_Pz_S_C1002_bc;
  Double I_ERI_F3y_Py_Pz_S_C1002_bc = I_ERI_G4y_S_Pz_S_C1002_bc+ABY*I_ERI_F3y_S_Pz_S_C1002_bc;
  Double I_ERI_F2yz_Py_Pz_S_C1002_bc = I_ERI_G3yz_S_Pz_S_C1002_bc+ABY*I_ERI_F2yz_S_Pz_S_C1002_bc;
  Double I_ERI_Fy2z_Py_Pz_S_C1002_bc = I_ERI_G2y2z_S_Pz_S_C1002_bc+ABY*I_ERI_Fy2z_S_Pz_S_C1002_bc;
  Double I_ERI_F3z_Py_Pz_S_C1002_bc = I_ERI_Gy3z_S_Pz_S_C1002_bc+ABY*I_ERI_F3z_S_Pz_S_C1002_bc;
  Double I_ERI_F2xz_Pz_Pz_S_C1002_bc = I_ERI_G2x2z_S_Pz_S_C1002_bc+ABZ*I_ERI_F2xz_S_Pz_S_C1002_bc;
  Double I_ERI_Fxyz_Pz_Pz_S_C1002_bc = I_ERI_Gxy2z_S_Pz_S_C1002_bc+ABZ*I_ERI_Fxyz_S_Pz_S_C1002_bc;
  Double I_ERI_Fx2z_Pz_Pz_S_C1002_bc = I_ERI_Gx3z_S_Pz_S_C1002_bc+ABZ*I_ERI_Fx2z_S_Pz_S_C1002_bc;
  Double I_ERI_F2yz_Pz_Pz_S_C1002_bc = I_ERI_G2y2z_S_Pz_S_C1002_bc+ABZ*I_ERI_F2yz_S_Pz_S_C1002_bc;
  Double I_ERI_Fy2z_Pz_Pz_S_C1002_bc = I_ERI_Gy3z_S_Pz_S_C1002_bc+ABZ*I_ERI_Fy2z_S_Pz_S_C1002_bc;
  Double I_ERI_F3z_Pz_Pz_S_C1002_bc = I_ERI_G4z_S_Pz_S_C1002_bc+ABZ*I_ERI_F3z_S_Pz_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1002_bc
   ************************************************************/
  Double I_ERI_D2x_D2x_Px_S_C1002_bc = I_ERI_F3x_Px_Px_S_C1002_bc+ABX*I_ERI_D2x_Px_Px_S_C1002_bc;
  Double I_ERI_Dxy_D2x_Px_S_C1002_bc = I_ERI_F2xy_Px_Px_S_C1002_bc+ABX*I_ERI_Dxy_Px_Px_S_C1002_bc;
  Double I_ERI_Dxz_D2x_Px_S_C1002_bc = I_ERI_F2xz_Px_Px_S_C1002_bc+ABX*I_ERI_Dxz_Px_Px_S_C1002_bc;
  Double I_ERI_D2y_D2x_Px_S_C1002_bc = I_ERI_Fx2y_Px_Px_S_C1002_bc+ABX*I_ERI_D2y_Px_Px_S_C1002_bc;
  Double I_ERI_Dyz_D2x_Px_S_C1002_bc = I_ERI_Fxyz_Px_Px_S_C1002_bc+ABX*I_ERI_Dyz_Px_Px_S_C1002_bc;
  Double I_ERI_D2z_D2x_Px_S_C1002_bc = I_ERI_Fx2z_Px_Px_S_C1002_bc+ABX*I_ERI_D2z_Px_Px_S_C1002_bc;
  Double I_ERI_D2x_Dxy_Px_S_C1002_bc = I_ERI_F2xy_Px_Px_S_C1002_bc+ABY*I_ERI_D2x_Px_Px_S_C1002_bc;
  Double I_ERI_Dxy_Dxy_Px_S_C1002_bc = I_ERI_Fx2y_Px_Px_S_C1002_bc+ABY*I_ERI_Dxy_Px_Px_S_C1002_bc;
  Double I_ERI_Dxz_Dxy_Px_S_C1002_bc = I_ERI_Fxyz_Px_Px_S_C1002_bc+ABY*I_ERI_Dxz_Px_Px_S_C1002_bc;
  Double I_ERI_D2y_Dxy_Px_S_C1002_bc = I_ERI_F3y_Px_Px_S_C1002_bc+ABY*I_ERI_D2y_Px_Px_S_C1002_bc;
  Double I_ERI_Dyz_Dxy_Px_S_C1002_bc = I_ERI_F2yz_Px_Px_S_C1002_bc+ABY*I_ERI_Dyz_Px_Px_S_C1002_bc;
  Double I_ERI_D2z_Dxy_Px_S_C1002_bc = I_ERI_Fy2z_Px_Px_S_C1002_bc+ABY*I_ERI_D2z_Px_Px_S_C1002_bc;
  Double I_ERI_D2x_Dxz_Px_S_C1002_bc = I_ERI_F2xz_Px_Px_S_C1002_bc+ABZ*I_ERI_D2x_Px_Px_S_C1002_bc;
  Double I_ERI_Dxy_Dxz_Px_S_C1002_bc = I_ERI_Fxyz_Px_Px_S_C1002_bc+ABZ*I_ERI_Dxy_Px_Px_S_C1002_bc;
  Double I_ERI_Dxz_Dxz_Px_S_C1002_bc = I_ERI_Fx2z_Px_Px_S_C1002_bc+ABZ*I_ERI_Dxz_Px_Px_S_C1002_bc;
  Double I_ERI_D2y_Dxz_Px_S_C1002_bc = I_ERI_F2yz_Px_Px_S_C1002_bc+ABZ*I_ERI_D2y_Px_Px_S_C1002_bc;
  Double I_ERI_Dyz_Dxz_Px_S_C1002_bc = I_ERI_Fy2z_Px_Px_S_C1002_bc+ABZ*I_ERI_Dyz_Px_Px_S_C1002_bc;
  Double I_ERI_D2z_Dxz_Px_S_C1002_bc = I_ERI_F3z_Px_Px_S_C1002_bc+ABZ*I_ERI_D2z_Px_Px_S_C1002_bc;
  Double I_ERI_D2x_D2y_Px_S_C1002_bc = I_ERI_F2xy_Py_Px_S_C1002_bc+ABY*I_ERI_D2x_Py_Px_S_C1002_bc;
  Double I_ERI_Dxy_D2y_Px_S_C1002_bc = I_ERI_Fx2y_Py_Px_S_C1002_bc+ABY*I_ERI_Dxy_Py_Px_S_C1002_bc;
  Double I_ERI_Dxz_D2y_Px_S_C1002_bc = I_ERI_Fxyz_Py_Px_S_C1002_bc+ABY*I_ERI_Dxz_Py_Px_S_C1002_bc;
  Double I_ERI_D2y_D2y_Px_S_C1002_bc = I_ERI_F3y_Py_Px_S_C1002_bc+ABY*I_ERI_D2y_Py_Px_S_C1002_bc;
  Double I_ERI_Dyz_D2y_Px_S_C1002_bc = I_ERI_F2yz_Py_Px_S_C1002_bc+ABY*I_ERI_Dyz_Py_Px_S_C1002_bc;
  Double I_ERI_D2z_D2y_Px_S_C1002_bc = I_ERI_Fy2z_Py_Px_S_C1002_bc+ABY*I_ERI_D2z_Py_Px_S_C1002_bc;
  Double I_ERI_D2x_Dyz_Px_S_C1002_bc = I_ERI_F2xz_Py_Px_S_C1002_bc+ABZ*I_ERI_D2x_Py_Px_S_C1002_bc;
  Double I_ERI_Dxy_Dyz_Px_S_C1002_bc = I_ERI_Fxyz_Py_Px_S_C1002_bc+ABZ*I_ERI_Dxy_Py_Px_S_C1002_bc;
  Double I_ERI_Dxz_Dyz_Px_S_C1002_bc = I_ERI_Fx2z_Py_Px_S_C1002_bc+ABZ*I_ERI_Dxz_Py_Px_S_C1002_bc;
  Double I_ERI_D2y_Dyz_Px_S_C1002_bc = I_ERI_F2yz_Py_Px_S_C1002_bc+ABZ*I_ERI_D2y_Py_Px_S_C1002_bc;
  Double I_ERI_Dyz_Dyz_Px_S_C1002_bc = I_ERI_Fy2z_Py_Px_S_C1002_bc+ABZ*I_ERI_Dyz_Py_Px_S_C1002_bc;
  Double I_ERI_D2z_Dyz_Px_S_C1002_bc = I_ERI_F3z_Py_Px_S_C1002_bc+ABZ*I_ERI_D2z_Py_Px_S_C1002_bc;
  Double I_ERI_D2x_D2z_Px_S_C1002_bc = I_ERI_F2xz_Pz_Px_S_C1002_bc+ABZ*I_ERI_D2x_Pz_Px_S_C1002_bc;
  Double I_ERI_Dxy_D2z_Px_S_C1002_bc = I_ERI_Fxyz_Pz_Px_S_C1002_bc+ABZ*I_ERI_Dxy_Pz_Px_S_C1002_bc;
  Double I_ERI_Dxz_D2z_Px_S_C1002_bc = I_ERI_Fx2z_Pz_Px_S_C1002_bc+ABZ*I_ERI_Dxz_Pz_Px_S_C1002_bc;
  Double I_ERI_D2y_D2z_Px_S_C1002_bc = I_ERI_F2yz_Pz_Px_S_C1002_bc+ABZ*I_ERI_D2y_Pz_Px_S_C1002_bc;
  Double I_ERI_Dyz_D2z_Px_S_C1002_bc = I_ERI_Fy2z_Pz_Px_S_C1002_bc+ABZ*I_ERI_Dyz_Pz_Px_S_C1002_bc;
  Double I_ERI_D2z_D2z_Px_S_C1002_bc = I_ERI_F3z_Pz_Px_S_C1002_bc+ABZ*I_ERI_D2z_Pz_Px_S_C1002_bc;
  Double I_ERI_D2x_D2x_Py_S_C1002_bc = I_ERI_F3x_Px_Py_S_C1002_bc+ABX*I_ERI_D2x_Px_Py_S_C1002_bc;
  Double I_ERI_Dxy_D2x_Py_S_C1002_bc = I_ERI_F2xy_Px_Py_S_C1002_bc+ABX*I_ERI_Dxy_Px_Py_S_C1002_bc;
  Double I_ERI_Dxz_D2x_Py_S_C1002_bc = I_ERI_F2xz_Px_Py_S_C1002_bc+ABX*I_ERI_Dxz_Px_Py_S_C1002_bc;
  Double I_ERI_D2y_D2x_Py_S_C1002_bc = I_ERI_Fx2y_Px_Py_S_C1002_bc+ABX*I_ERI_D2y_Px_Py_S_C1002_bc;
  Double I_ERI_Dyz_D2x_Py_S_C1002_bc = I_ERI_Fxyz_Px_Py_S_C1002_bc+ABX*I_ERI_Dyz_Px_Py_S_C1002_bc;
  Double I_ERI_D2z_D2x_Py_S_C1002_bc = I_ERI_Fx2z_Px_Py_S_C1002_bc+ABX*I_ERI_D2z_Px_Py_S_C1002_bc;
  Double I_ERI_D2x_Dxy_Py_S_C1002_bc = I_ERI_F2xy_Px_Py_S_C1002_bc+ABY*I_ERI_D2x_Px_Py_S_C1002_bc;
  Double I_ERI_Dxy_Dxy_Py_S_C1002_bc = I_ERI_Fx2y_Px_Py_S_C1002_bc+ABY*I_ERI_Dxy_Px_Py_S_C1002_bc;
  Double I_ERI_Dxz_Dxy_Py_S_C1002_bc = I_ERI_Fxyz_Px_Py_S_C1002_bc+ABY*I_ERI_Dxz_Px_Py_S_C1002_bc;
  Double I_ERI_D2y_Dxy_Py_S_C1002_bc = I_ERI_F3y_Px_Py_S_C1002_bc+ABY*I_ERI_D2y_Px_Py_S_C1002_bc;
  Double I_ERI_Dyz_Dxy_Py_S_C1002_bc = I_ERI_F2yz_Px_Py_S_C1002_bc+ABY*I_ERI_Dyz_Px_Py_S_C1002_bc;
  Double I_ERI_D2z_Dxy_Py_S_C1002_bc = I_ERI_Fy2z_Px_Py_S_C1002_bc+ABY*I_ERI_D2z_Px_Py_S_C1002_bc;
  Double I_ERI_D2x_Dxz_Py_S_C1002_bc = I_ERI_F2xz_Px_Py_S_C1002_bc+ABZ*I_ERI_D2x_Px_Py_S_C1002_bc;
  Double I_ERI_Dxy_Dxz_Py_S_C1002_bc = I_ERI_Fxyz_Px_Py_S_C1002_bc+ABZ*I_ERI_Dxy_Px_Py_S_C1002_bc;
  Double I_ERI_Dxz_Dxz_Py_S_C1002_bc = I_ERI_Fx2z_Px_Py_S_C1002_bc+ABZ*I_ERI_Dxz_Px_Py_S_C1002_bc;
  Double I_ERI_D2y_Dxz_Py_S_C1002_bc = I_ERI_F2yz_Px_Py_S_C1002_bc+ABZ*I_ERI_D2y_Px_Py_S_C1002_bc;
  Double I_ERI_Dyz_Dxz_Py_S_C1002_bc = I_ERI_Fy2z_Px_Py_S_C1002_bc+ABZ*I_ERI_Dyz_Px_Py_S_C1002_bc;
  Double I_ERI_D2z_Dxz_Py_S_C1002_bc = I_ERI_F3z_Px_Py_S_C1002_bc+ABZ*I_ERI_D2z_Px_Py_S_C1002_bc;
  Double I_ERI_D2x_D2y_Py_S_C1002_bc = I_ERI_F2xy_Py_Py_S_C1002_bc+ABY*I_ERI_D2x_Py_Py_S_C1002_bc;
  Double I_ERI_Dxy_D2y_Py_S_C1002_bc = I_ERI_Fx2y_Py_Py_S_C1002_bc+ABY*I_ERI_Dxy_Py_Py_S_C1002_bc;
  Double I_ERI_Dxz_D2y_Py_S_C1002_bc = I_ERI_Fxyz_Py_Py_S_C1002_bc+ABY*I_ERI_Dxz_Py_Py_S_C1002_bc;
  Double I_ERI_D2y_D2y_Py_S_C1002_bc = I_ERI_F3y_Py_Py_S_C1002_bc+ABY*I_ERI_D2y_Py_Py_S_C1002_bc;
  Double I_ERI_Dyz_D2y_Py_S_C1002_bc = I_ERI_F2yz_Py_Py_S_C1002_bc+ABY*I_ERI_Dyz_Py_Py_S_C1002_bc;
  Double I_ERI_D2z_D2y_Py_S_C1002_bc = I_ERI_Fy2z_Py_Py_S_C1002_bc+ABY*I_ERI_D2z_Py_Py_S_C1002_bc;
  Double I_ERI_D2x_Dyz_Py_S_C1002_bc = I_ERI_F2xz_Py_Py_S_C1002_bc+ABZ*I_ERI_D2x_Py_Py_S_C1002_bc;
  Double I_ERI_Dxy_Dyz_Py_S_C1002_bc = I_ERI_Fxyz_Py_Py_S_C1002_bc+ABZ*I_ERI_Dxy_Py_Py_S_C1002_bc;
  Double I_ERI_Dxz_Dyz_Py_S_C1002_bc = I_ERI_Fx2z_Py_Py_S_C1002_bc+ABZ*I_ERI_Dxz_Py_Py_S_C1002_bc;
  Double I_ERI_D2y_Dyz_Py_S_C1002_bc = I_ERI_F2yz_Py_Py_S_C1002_bc+ABZ*I_ERI_D2y_Py_Py_S_C1002_bc;
  Double I_ERI_Dyz_Dyz_Py_S_C1002_bc = I_ERI_Fy2z_Py_Py_S_C1002_bc+ABZ*I_ERI_Dyz_Py_Py_S_C1002_bc;
  Double I_ERI_D2z_Dyz_Py_S_C1002_bc = I_ERI_F3z_Py_Py_S_C1002_bc+ABZ*I_ERI_D2z_Py_Py_S_C1002_bc;
  Double I_ERI_D2x_D2z_Py_S_C1002_bc = I_ERI_F2xz_Pz_Py_S_C1002_bc+ABZ*I_ERI_D2x_Pz_Py_S_C1002_bc;
  Double I_ERI_Dxy_D2z_Py_S_C1002_bc = I_ERI_Fxyz_Pz_Py_S_C1002_bc+ABZ*I_ERI_Dxy_Pz_Py_S_C1002_bc;
  Double I_ERI_Dxz_D2z_Py_S_C1002_bc = I_ERI_Fx2z_Pz_Py_S_C1002_bc+ABZ*I_ERI_Dxz_Pz_Py_S_C1002_bc;
  Double I_ERI_D2y_D2z_Py_S_C1002_bc = I_ERI_F2yz_Pz_Py_S_C1002_bc+ABZ*I_ERI_D2y_Pz_Py_S_C1002_bc;
  Double I_ERI_Dyz_D2z_Py_S_C1002_bc = I_ERI_Fy2z_Pz_Py_S_C1002_bc+ABZ*I_ERI_Dyz_Pz_Py_S_C1002_bc;
  Double I_ERI_D2z_D2z_Py_S_C1002_bc = I_ERI_F3z_Pz_Py_S_C1002_bc+ABZ*I_ERI_D2z_Pz_Py_S_C1002_bc;
  Double I_ERI_D2x_D2x_Pz_S_C1002_bc = I_ERI_F3x_Px_Pz_S_C1002_bc+ABX*I_ERI_D2x_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxy_D2x_Pz_S_C1002_bc = I_ERI_F2xy_Px_Pz_S_C1002_bc+ABX*I_ERI_Dxy_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxz_D2x_Pz_S_C1002_bc = I_ERI_F2xz_Px_Pz_S_C1002_bc+ABX*I_ERI_Dxz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2y_D2x_Pz_S_C1002_bc = I_ERI_Fx2y_Px_Pz_S_C1002_bc+ABX*I_ERI_D2y_Px_Pz_S_C1002_bc;
  Double I_ERI_Dyz_D2x_Pz_S_C1002_bc = I_ERI_Fxyz_Px_Pz_S_C1002_bc+ABX*I_ERI_Dyz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2z_D2x_Pz_S_C1002_bc = I_ERI_Fx2z_Px_Pz_S_C1002_bc+ABX*I_ERI_D2z_Px_Pz_S_C1002_bc;
  Double I_ERI_D2x_Dxy_Pz_S_C1002_bc = I_ERI_F2xy_Px_Pz_S_C1002_bc+ABY*I_ERI_D2x_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Dxy_Pz_S_C1002_bc = I_ERI_Fx2y_Px_Pz_S_C1002_bc+ABY*I_ERI_Dxy_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Dxy_Pz_S_C1002_bc = I_ERI_Fxyz_Px_Pz_S_C1002_bc+ABY*I_ERI_Dxz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2y_Dxy_Pz_S_C1002_bc = I_ERI_F3y_Px_Pz_S_C1002_bc+ABY*I_ERI_D2y_Px_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Dxy_Pz_S_C1002_bc = I_ERI_F2yz_Px_Pz_S_C1002_bc+ABY*I_ERI_Dyz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2z_Dxy_Pz_S_C1002_bc = I_ERI_Fy2z_Px_Pz_S_C1002_bc+ABY*I_ERI_D2z_Px_Pz_S_C1002_bc;
  Double I_ERI_D2x_Dxz_Pz_S_C1002_bc = I_ERI_F2xz_Px_Pz_S_C1002_bc+ABZ*I_ERI_D2x_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Dxz_Pz_S_C1002_bc = I_ERI_Fxyz_Px_Pz_S_C1002_bc+ABZ*I_ERI_Dxy_Px_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Dxz_Pz_S_C1002_bc = I_ERI_Fx2z_Px_Pz_S_C1002_bc+ABZ*I_ERI_Dxz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2y_Dxz_Pz_S_C1002_bc = I_ERI_F2yz_Px_Pz_S_C1002_bc+ABZ*I_ERI_D2y_Px_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Dxz_Pz_S_C1002_bc = I_ERI_Fy2z_Px_Pz_S_C1002_bc+ABZ*I_ERI_Dyz_Px_Pz_S_C1002_bc;
  Double I_ERI_D2z_Dxz_Pz_S_C1002_bc = I_ERI_F3z_Px_Pz_S_C1002_bc+ABZ*I_ERI_D2z_Px_Pz_S_C1002_bc;
  Double I_ERI_D2x_D2y_Pz_S_C1002_bc = I_ERI_F2xy_Py_Pz_S_C1002_bc+ABY*I_ERI_D2x_Py_Pz_S_C1002_bc;
  Double I_ERI_Dxy_D2y_Pz_S_C1002_bc = I_ERI_Fx2y_Py_Pz_S_C1002_bc+ABY*I_ERI_Dxy_Py_Pz_S_C1002_bc;
  Double I_ERI_Dxz_D2y_Pz_S_C1002_bc = I_ERI_Fxyz_Py_Pz_S_C1002_bc+ABY*I_ERI_Dxz_Py_Pz_S_C1002_bc;
  Double I_ERI_D2y_D2y_Pz_S_C1002_bc = I_ERI_F3y_Py_Pz_S_C1002_bc+ABY*I_ERI_D2y_Py_Pz_S_C1002_bc;
  Double I_ERI_Dyz_D2y_Pz_S_C1002_bc = I_ERI_F2yz_Py_Pz_S_C1002_bc+ABY*I_ERI_Dyz_Py_Pz_S_C1002_bc;
  Double I_ERI_D2z_D2y_Pz_S_C1002_bc = I_ERI_Fy2z_Py_Pz_S_C1002_bc+ABY*I_ERI_D2z_Py_Pz_S_C1002_bc;
  Double I_ERI_D2x_Dyz_Pz_S_C1002_bc = I_ERI_F2xz_Py_Pz_S_C1002_bc+ABZ*I_ERI_D2x_Py_Pz_S_C1002_bc;
  Double I_ERI_Dxy_Dyz_Pz_S_C1002_bc = I_ERI_Fxyz_Py_Pz_S_C1002_bc+ABZ*I_ERI_Dxy_Py_Pz_S_C1002_bc;
  Double I_ERI_Dxz_Dyz_Pz_S_C1002_bc = I_ERI_Fx2z_Py_Pz_S_C1002_bc+ABZ*I_ERI_Dxz_Py_Pz_S_C1002_bc;
  Double I_ERI_D2y_Dyz_Pz_S_C1002_bc = I_ERI_F2yz_Py_Pz_S_C1002_bc+ABZ*I_ERI_D2y_Py_Pz_S_C1002_bc;
  Double I_ERI_Dyz_Dyz_Pz_S_C1002_bc = I_ERI_Fy2z_Py_Pz_S_C1002_bc+ABZ*I_ERI_Dyz_Py_Pz_S_C1002_bc;
  Double I_ERI_D2z_Dyz_Pz_S_C1002_bc = I_ERI_F3z_Py_Pz_S_C1002_bc+ABZ*I_ERI_D2z_Py_Pz_S_C1002_bc;
  Double I_ERI_D2x_D2z_Pz_S_C1002_bc = I_ERI_F2xz_Pz_Pz_S_C1002_bc+ABZ*I_ERI_D2x_Pz_Pz_S_C1002_bc;
  Double I_ERI_Dxy_D2z_Pz_S_C1002_bc = I_ERI_Fxyz_Pz_Pz_S_C1002_bc+ABZ*I_ERI_Dxy_Pz_Pz_S_C1002_bc;
  Double I_ERI_Dxz_D2z_Pz_S_C1002_bc = I_ERI_Fx2z_Pz_Pz_S_C1002_bc+ABZ*I_ERI_Dxz_Pz_Pz_S_C1002_bc;
  Double I_ERI_D2y_D2z_Pz_S_C1002_bc = I_ERI_F2yz_Pz_Pz_S_C1002_bc+ABZ*I_ERI_D2y_Pz_Pz_S_C1002_bc;
  Double I_ERI_Dyz_D2z_Pz_S_C1002_bc = I_ERI_Fy2z_Pz_Pz_S_C1002_bc+ABZ*I_ERI_Dyz_Pz_Pz_S_C1002_bc;
  Double I_ERI_D2z_D2z_Pz_S_C1002_bc = I_ERI_F3z_Pz_Pz_S_C1002_bc+ABZ*I_ERI_D2z_Pz_Pz_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1002_cc
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S_C1002_cc = I_ERI_F3x_S_D2x_S_C1002_cc+ABX*I_ERI_D2x_S_D2x_S_C1002_cc;
  Double I_ERI_Dxy_Px_D2x_S_C1002_cc = I_ERI_F2xy_S_D2x_S_C1002_cc+ABX*I_ERI_Dxy_S_D2x_S_C1002_cc;
  Double I_ERI_Dxz_Px_D2x_S_C1002_cc = I_ERI_F2xz_S_D2x_S_C1002_cc+ABX*I_ERI_Dxz_S_D2x_S_C1002_cc;
  Double I_ERI_D2y_Px_D2x_S_C1002_cc = I_ERI_Fx2y_S_D2x_S_C1002_cc+ABX*I_ERI_D2y_S_D2x_S_C1002_cc;
  Double I_ERI_Dyz_Px_D2x_S_C1002_cc = I_ERI_Fxyz_S_D2x_S_C1002_cc+ABX*I_ERI_Dyz_S_D2x_S_C1002_cc;
  Double I_ERI_D2z_Px_D2x_S_C1002_cc = I_ERI_Fx2z_S_D2x_S_C1002_cc+ABX*I_ERI_D2z_S_D2x_S_C1002_cc;
  Double I_ERI_D2x_Py_D2x_S_C1002_cc = I_ERI_F2xy_S_D2x_S_C1002_cc+ABY*I_ERI_D2x_S_D2x_S_C1002_cc;
  Double I_ERI_Dxy_Py_D2x_S_C1002_cc = I_ERI_Fx2y_S_D2x_S_C1002_cc+ABY*I_ERI_Dxy_S_D2x_S_C1002_cc;
  Double I_ERI_Dxz_Py_D2x_S_C1002_cc = I_ERI_Fxyz_S_D2x_S_C1002_cc+ABY*I_ERI_Dxz_S_D2x_S_C1002_cc;
  Double I_ERI_D2y_Py_D2x_S_C1002_cc = I_ERI_F3y_S_D2x_S_C1002_cc+ABY*I_ERI_D2y_S_D2x_S_C1002_cc;
  Double I_ERI_Dyz_Py_D2x_S_C1002_cc = I_ERI_F2yz_S_D2x_S_C1002_cc+ABY*I_ERI_Dyz_S_D2x_S_C1002_cc;
  Double I_ERI_D2z_Py_D2x_S_C1002_cc = I_ERI_Fy2z_S_D2x_S_C1002_cc+ABY*I_ERI_D2z_S_D2x_S_C1002_cc;
  Double I_ERI_D2x_Pz_D2x_S_C1002_cc = I_ERI_F2xz_S_D2x_S_C1002_cc+ABZ*I_ERI_D2x_S_D2x_S_C1002_cc;
  Double I_ERI_Dxy_Pz_D2x_S_C1002_cc = I_ERI_Fxyz_S_D2x_S_C1002_cc+ABZ*I_ERI_Dxy_S_D2x_S_C1002_cc;
  Double I_ERI_Dxz_Pz_D2x_S_C1002_cc = I_ERI_Fx2z_S_D2x_S_C1002_cc+ABZ*I_ERI_Dxz_S_D2x_S_C1002_cc;
  Double I_ERI_D2y_Pz_D2x_S_C1002_cc = I_ERI_F2yz_S_D2x_S_C1002_cc+ABZ*I_ERI_D2y_S_D2x_S_C1002_cc;
  Double I_ERI_Dyz_Pz_D2x_S_C1002_cc = I_ERI_Fy2z_S_D2x_S_C1002_cc+ABZ*I_ERI_Dyz_S_D2x_S_C1002_cc;
  Double I_ERI_D2z_Pz_D2x_S_C1002_cc = I_ERI_F3z_S_D2x_S_C1002_cc+ABZ*I_ERI_D2z_S_D2x_S_C1002_cc;
  Double I_ERI_D2x_Px_Dxy_S_C1002_cc = I_ERI_F3x_S_Dxy_S_C1002_cc+ABX*I_ERI_D2x_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxy_Px_Dxy_S_C1002_cc = I_ERI_F2xy_S_Dxy_S_C1002_cc+ABX*I_ERI_Dxy_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxz_Px_Dxy_S_C1002_cc = I_ERI_F2xz_S_Dxy_S_C1002_cc+ABX*I_ERI_Dxz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2y_Px_Dxy_S_C1002_cc = I_ERI_Fx2y_S_Dxy_S_C1002_cc+ABX*I_ERI_D2y_S_Dxy_S_C1002_cc;
  Double I_ERI_Dyz_Px_Dxy_S_C1002_cc = I_ERI_Fxyz_S_Dxy_S_C1002_cc+ABX*I_ERI_Dyz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2z_Px_Dxy_S_C1002_cc = I_ERI_Fx2z_S_Dxy_S_C1002_cc+ABX*I_ERI_D2z_S_Dxy_S_C1002_cc;
  Double I_ERI_D2x_Py_Dxy_S_C1002_cc = I_ERI_F2xy_S_Dxy_S_C1002_cc+ABY*I_ERI_D2x_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxy_Py_Dxy_S_C1002_cc = I_ERI_Fx2y_S_Dxy_S_C1002_cc+ABY*I_ERI_Dxy_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxz_Py_Dxy_S_C1002_cc = I_ERI_Fxyz_S_Dxy_S_C1002_cc+ABY*I_ERI_Dxz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2y_Py_Dxy_S_C1002_cc = I_ERI_F3y_S_Dxy_S_C1002_cc+ABY*I_ERI_D2y_S_Dxy_S_C1002_cc;
  Double I_ERI_Dyz_Py_Dxy_S_C1002_cc = I_ERI_F2yz_S_Dxy_S_C1002_cc+ABY*I_ERI_Dyz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2z_Py_Dxy_S_C1002_cc = I_ERI_Fy2z_S_Dxy_S_C1002_cc+ABY*I_ERI_D2z_S_Dxy_S_C1002_cc;
  Double I_ERI_D2x_Pz_Dxy_S_C1002_cc = I_ERI_F2xz_S_Dxy_S_C1002_cc+ABZ*I_ERI_D2x_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxy_Pz_Dxy_S_C1002_cc = I_ERI_Fxyz_S_Dxy_S_C1002_cc+ABZ*I_ERI_Dxy_S_Dxy_S_C1002_cc;
  Double I_ERI_Dxz_Pz_Dxy_S_C1002_cc = I_ERI_Fx2z_S_Dxy_S_C1002_cc+ABZ*I_ERI_Dxz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2y_Pz_Dxy_S_C1002_cc = I_ERI_F2yz_S_Dxy_S_C1002_cc+ABZ*I_ERI_D2y_S_Dxy_S_C1002_cc;
  Double I_ERI_Dyz_Pz_Dxy_S_C1002_cc = I_ERI_Fy2z_S_Dxy_S_C1002_cc+ABZ*I_ERI_Dyz_S_Dxy_S_C1002_cc;
  Double I_ERI_D2z_Pz_Dxy_S_C1002_cc = I_ERI_F3z_S_Dxy_S_C1002_cc+ABZ*I_ERI_D2z_S_Dxy_S_C1002_cc;
  Double I_ERI_D2x_Px_Dxz_S_C1002_cc = I_ERI_F3x_S_Dxz_S_C1002_cc+ABX*I_ERI_D2x_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxy_Px_Dxz_S_C1002_cc = I_ERI_F2xy_S_Dxz_S_C1002_cc+ABX*I_ERI_Dxy_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxz_Px_Dxz_S_C1002_cc = I_ERI_F2xz_S_Dxz_S_C1002_cc+ABX*I_ERI_Dxz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2y_Px_Dxz_S_C1002_cc = I_ERI_Fx2y_S_Dxz_S_C1002_cc+ABX*I_ERI_D2y_S_Dxz_S_C1002_cc;
  Double I_ERI_Dyz_Px_Dxz_S_C1002_cc = I_ERI_Fxyz_S_Dxz_S_C1002_cc+ABX*I_ERI_Dyz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2z_Px_Dxz_S_C1002_cc = I_ERI_Fx2z_S_Dxz_S_C1002_cc+ABX*I_ERI_D2z_S_Dxz_S_C1002_cc;
  Double I_ERI_D2x_Py_Dxz_S_C1002_cc = I_ERI_F2xy_S_Dxz_S_C1002_cc+ABY*I_ERI_D2x_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxy_Py_Dxz_S_C1002_cc = I_ERI_Fx2y_S_Dxz_S_C1002_cc+ABY*I_ERI_Dxy_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxz_Py_Dxz_S_C1002_cc = I_ERI_Fxyz_S_Dxz_S_C1002_cc+ABY*I_ERI_Dxz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2y_Py_Dxz_S_C1002_cc = I_ERI_F3y_S_Dxz_S_C1002_cc+ABY*I_ERI_D2y_S_Dxz_S_C1002_cc;
  Double I_ERI_Dyz_Py_Dxz_S_C1002_cc = I_ERI_F2yz_S_Dxz_S_C1002_cc+ABY*I_ERI_Dyz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2z_Py_Dxz_S_C1002_cc = I_ERI_Fy2z_S_Dxz_S_C1002_cc+ABY*I_ERI_D2z_S_Dxz_S_C1002_cc;
  Double I_ERI_D2x_Pz_Dxz_S_C1002_cc = I_ERI_F2xz_S_Dxz_S_C1002_cc+ABZ*I_ERI_D2x_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxy_Pz_Dxz_S_C1002_cc = I_ERI_Fxyz_S_Dxz_S_C1002_cc+ABZ*I_ERI_Dxy_S_Dxz_S_C1002_cc;
  Double I_ERI_Dxz_Pz_Dxz_S_C1002_cc = I_ERI_Fx2z_S_Dxz_S_C1002_cc+ABZ*I_ERI_Dxz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2y_Pz_Dxz_S_C1002_cc = I_ERI_F2yz_S_Dxz_S_C1002_cc+ABZ*I_ERI_D2y_S_Dxz_S_C1002_cc;
  Double I_ERI_Dyz_Pz_Dxz_S_C1002_cc = I_ERI_Fy2z_S_Dxz_S_C1002_cc+ABZ*I_ERI_Dyz_S_Dxz_S_C1002_cc;
  Double I_ERI_D2z_Pz_Dxz_S_C1002_cc = I_ERI_F3z_S_Dxz_S_C1002_cc+ABZ*I_ERI_D2z_S_Dxz_S_C1002_cc;
  Double I_ERI_D2x_Px_D2y_S_C1002_cc = I_ERI_F3x_S_D2y_S_C1002_cc+ABX*I_ERI_D2x_S_D2y_S_C1002_cc;
  Double I_ERI_Dxy_Px_D2y_S_C1002_cc = I_ERI_F2xy_S_D2y_S_C1002_cc+ABX*I_ERI_Dxy_S_D2y_S_C1002_cc;
  Double I_ERI_Dxz_Px_D2y_S_C1002_cc = I_ERI_F2xz_S_D2y_S_C1002_cc+ABX*I_ERI_Dxz_S_D2y_S_C1002_cc;
  Double I_ERI_D2y_Px_D2y_S_C1002_cc = I_ERI_Fx2y_S_D2y_S_C1002_cc+ABX*I_ERI_D2y_S_D2y_S_C1002_cc;
  Double I_ERI_Dyz_Px_D2y_S_C1002_cc = I_ERI_Fxyz_S_D2y_S_C1002_cc+ABX*I_ERI_Dyz_S_D2y_S_C1002_cc;
  Double I_ERI_D2z_Px_D2y_S_C1002_cc = I_ERI_Fx2z_S_D2y_S_C1002_cc+ABX*I_ERI_D2z_S_D2y_S_C1002_cc;
  Double I_ERI_D2x_Py_D2y_S_C1002_cc = I_ERI_F2xy_S_D2y_S_C1002_cc+ABY*I_ERI_D2x_S_D2y_S_C1002_cc;
  Double I_ERI_Dxy_Py_D2y_S_C1002_cc = I_ERI_Fx2y_S_D2y_S_C1002_cc+ABY*I_ERI_Dxy_S_D2y_S_C1002_cc;
  Double I_ERI_Dxz_Py_D2y_S_C1002_cc = I_ERI_Fxyz_S_D2y_S_C1002_cc+ABY*I_ERI_Dxz_S_D2y_S_C1002_cc;
  Double I_ERI_D2y_Py_D2y_S_C1002_cc = I_ERI_F3y_S_D2y_S_C1002_cc+ABY*I_ERI_D2y_S_D2y_S_C1002_cc;
  Double I_ERI_Dyz_Py_D2y_S_C1002_cc = I_ERI_F2yz_S_D2y_S_C1002_cc+ABY*I_ERI_Dyz_S_D2y_S_C1002_cc;
  Double I_ERI_D2z_Py_D2y_S_C1002_cc = I_ERI_Fy2z_S_D2y_S_C1002_cc+ABY*I_ERI_D2z_S_D2y_S_C1002_cc;
  Double I_ERI_D2x_Pz_D2y_S_C1002_cc = I_ERI_F2xz_S_D2y_S_C1002_cc+ABZ*I_ERI_D2x_S_D2y_S_C1002_cc;
  Double I_ERI_Dxy_Pz_D2y_S_C1002_cc = I_ERI_Fxyz_S_D2y_S_C1002_cc+ABZ*I_ERI_Dxy_S_D2y_S_C1002_cc;
  Double I_ERI_Dxz_Pz_D2y_S_C1002_cc = I_ERI_Fx2z_S_D2y_S_C1002_cc+ABZ*I_ERI_Dxz_S_D2y_S_C1002_cc;
  Double I_ERI_D2y_Pz_D2y_S_C1002_cc = I_ERI_F2yz_S_D2y_S_C1002_cc+ABZ*I_ERI_D2y_S_D2y_S_C1002_cc;
  Double I_ERI_Dyz_Pz_D2y_S_C1002_cc = I_ERI_Fy2z_S_D2y_S_C1002_cc+ABZ*I_ERI_Dyz_S_D2y_S_C1002_cc;
  Double I_ERI_D2z_Pz_D2y_S_C1002_cc = I_ERI_F3z_S_D2y_S_C1002_cc+ABZ*I_ERI_D2z_S_D2y_S_C1002_cc;
  Double I_ERI_D2x_Px_Dyz_S_C1002_cc = I_ERI_F3x_S_Dyz_S_C1002_cc+ABX*I_ERI_D2x_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxy_Px_Dyz_S_C1002_cc = I_ERI_F2xy_S_Dyz_S_C1002_cc+ABX*I_ERI_Dxy_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxz_Px_Dyz_S_C1002_cc = I_ERI_F2xz_S_Dyz_S_C1002_cc+ABX*I_ERI_Dxz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2y_Px_Dyz_S_C1002_cc = I_ERI_Fx2y_S_Dyz_S_C1002_cc+ABX*I_ERI_D2y_S_Dyz_S_C1002_cc;
  Double I_ERI_Dyz_Px_Dyz_S_C1002_cc = I_ERI_Fxyz_S_Dyz_S_C1002_cc+ABX*I_ERI_Dyz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2z_Px_Dyz_S_C1002_cc = I_ERI_Fx2z_S_Dyz_S_C1002_cc+ABX*I_ERI_D2z_S_Dyz_S_C1002_cc;
  Double I_ERI_D2x_Py_Dyz_S_C1002_cc = I_ERI_F2xy_S_Dyz_S_C1002_cc+ABY*I_ERI_D2x_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxy_Py_Dyz_S_C1002_cc = I_ERI_Fx2y_S_Dyz_S_C1002_cc+ABY*I_ERI_Dxy_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxz_Py_Dyz_S_C1002_cc = I_ERI_Fxyz_S_Dyz_S_C1002_cc+ABY*I_ERI_Dxz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2y_Py_Dyz_S_C1002_cc = I_ERI_F3y_S_Dyz_S_C1002_cc+ABY*I_ERI_D2y_S_Dyz_S_C1002_cc;
  Double I_ERI_Dyz_Py_Dyz_S_C1002_cc = I_ERI_F2yz_S_Dyz_S_C1002_cc+ABY*I_ERI_Dyz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2z_Py_Dyz_S_C1002_cc = I_ERI_Fy2z_S_Dyz_S_C1002_cc+ABY*I_ERI_D2z_S_Dyz_S_C1002_cc;
  Double I_ERI_D2x_Pz_Dyz_S_C1002_cc = I_ERI_F2xz_S_Dyz_S_C1002_cc+ABZ*I_ERI_D2x_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxy_Pz_Dyz_S_C1002_cc = I_ERI_Fxyz_S_Dyz_S_C1002_cc+ABZ*I_ERI_Dxy_S_Dyz_S_C1002_cc;
  Double I_ERI_Dxz_Pz_Dyz_S_C1002_cc = I_ERI_Fx2z_S_Dyz_S_C1002_cc+ABZ*I_ERI_Dxz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2y_Pz_Dyz_S_C1002_cc = I_ERI_F2yz_S_Dyz_S_C1002_cc+ABZ*I_ERI_D2y_S_Dyz_S_C1002_cc;
  Double I_ERI_Dyz_Pz_Dyz_S_C1002_cc = I_ERI_Fy2z_S_Dyz_S_C1002_cc+ABZ*I_ERI_Dyz_S_Dyz_S_C1002_cc;
  Double I_ERI_D2z_Pz_Dyz_S_C1002_cc = I_ERI_F3z_S_Dyz_S_C1002_cc+ABZ*I_ERI_D2z_S_Dyz_S_C1002_cc;
  Double I_ERI_D2x_Px_D2z_S_C1002_cc = I_ERI_F3x_S_D2z_S_C1002_cc+ABX*I_ERI_D2x_S_D2z_S_C1002_cc;
  Double I_ERI_Dxy_Px_D2z_S_C1002_cc = I_ERI_F2xy_S_D2z_S_C1002_cc+ABX*I_ERI_Dxy_S_D2z_S_C1002_cc;
  Double I_ERI_Dxz_Px_D2z_S_C1002_cc = I_ERI_F2xz_S_D2z_S_C1002_cc+ABX*I_ERI_Dxz_S_D2z_S_C1002_cc;
  Double I_ERI_D2y_Px_D2z_S_C1002_cc = I_ERI_Fx2y_S_D2z_S_C1002_cc+ABX*I_ERI_D2y_S_D2z_S_C1002_cc;
  Double I_ERI_Dyz_Px_D2z_S_C1002_cc = I_ERI_Fxyz_S_D2z_S_C1002_cc+ABX*I_ERI_Dyz_S_D2z_S_C1002_cc;
  Double I_ERI_D2z_Px_D2z_S_C1002_cc = I_ERI_Fx2z_S_D2z_S_C1002_cc+ABX*I_ERI_D2z_S_D2z_S_C1002_cc;
  Double I_ERI_D2x_Py_D2z_S_C1002_cc = I_ERI_F2xy_S_D2z_S_C1002_cc+ABY*I_ERI_D2x_S_D2z_S_C1002_cc;
  Double I_ERI_Dxy_Py_D2z_S_C1002_cc = I_ERI_Fx2y_S_D2z_S_C1002_cc+ABY*I_ERI_Dxy_S_D2z_S_C1002_cc;
  Double I_ERI_Dxz_Py_D2z_S_C1002_cc = I_ERI_Fxyz_S_D2z_S_C1002_cc+ABY*I_ERI_Dxz_S_D2z_S_C1002_cc;
  Double I_ERI_D2y_Py_D2z_S_C1002_cc = I_ERI_F3y_S_D2z_S_C1002_cc+ABY*I_ERI_D2y_S_D2z_S_C1002_cc;
  Double I_ERI_Dyz_Py_D2z_S_C1002_cc = I_ERI_F2yz_S_D2z_S_C1002_cc+ABY*I_ERI_Dyz_S_D2z_S_C1002_cc;
  Double I_ERI_D2z_Py_D2z_S_C1002_cc = I_ERI_Fy2z_S_D2z_S_C1002_cc+ABY*I_ERI_D2z_S_D2z_S_C1002_cc;
  Double I_ERI_D2x_Pz_D2z_S_C1002_cc = I_ERI_F2xz_S_D2z_S_C1002_cc+ABZ*I_ERI_D2x_S_D2z_S_C1002_cc;
  Double I_ERI_Dxy_Pz_D2z_S_C1002_cc = I_ERI_Fxyz_S_D2z_S_C1002_cc+ABZ*I_ERI_Dxy_S_D2z_S_C1002_cc;
  Double I_ERI_Dxz_Pz_D2z_S_C1002_cc = I_ERI_Fx2z_S_D2z_S_C1002_cc+ABZ*I_ERI_Dxz_S_D2z_S_C1002_cc;
  Double I_ERI_D2y_Pz_D2z_S_C1002_cc = I_ERI_F2yz_S_D2z_S_C1002_cc+ABZ*I_ERI_D2y_S_D2z_S_C1002_cc;
  Double I_ERI_Dyz_Pz_D2z_S_C1002_cc = I_ERI_Fy2z_S_D2z_S_C1002_cc+ABZ*I_ERI_Dyz_S_D2z_S_C1002_cc;
  Double I_ERI_D2z_Pz_D2z_S_C1002_cc = I_ERI_F3z_S_D2z_S_C1002_cc+ABZ*I_ERI_D2z_S_D2z_S_C1002_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_G4x_S_S_S_C2_aa-2.0E0*2*I_ERI_D2x_S_S_S_C2_a-2.0E0*3*I_ERI_D2x_S_S_S_C2_a+2*1*I_ERI_S_S_S_S_C2;
  abcd[1] = 4.0E0*I_ERI_G3xy_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxy_S_S_S_C2_a-2.0E0*2*I_ERI_Dxy_S_S_S_C2_a;
  abcd[2] = 4.0E0*I_ERI_G3xz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxz_S_S_S_C2_a-2.0E0*2*I_ERI_Dxz_S_S_S_C2_a;
  abcd[3] = 4.0E0*I_ERI_G2x2y_S_S_S_C2_aa-2.0E0*1*I_ERI_D2y_S_S_S_C2_a;
  abcd[4] = 4.0E0*I_ERI_G2xyz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dyz_S_S_S_C2_a;
  abcd[5] = 4.0E0*I_ERI_G2x2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2z_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[6] = 4.0E0*I_ERI_G4x_Px_S_S_C1002_aa-2.0E0*2*I_ERI_D2x_Px_S_S_C1002_a-2.0E0*3*I_ERI_D2x_Px_S_S_C1002_a+2*1*I_ERI_S_Px_S_S_C1002;
  abcd[7] = 4.0E0*I_ERI_G3xy_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[8] = 4.0E0*I_ERI_G3xz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[9] = 4.0E0*I_ERI_G2x2y_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_a;
  abcd[10] = 4.0E0*I_ERI_G2xyz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[11] = 4.0E0*I_ERI_G2x2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_a;
  abcd[12] = 4.0E0*I_ERI_G4x_Py_S_S_C1002_aa-2.0E0*2*I_ERI_D2x_Py_S_S_C1002_a-2.0E0*3*I_ERI_D2x_Py_S_S_C1002_a+2*1*I_ERI_S_Py_S_S_C1002;
  abcd[13] = 4.0E0*I_ERI_G3xy_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[14] = 4.0E0*I_ERI_G3xz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[15] = 4.0E0*I_ERI_G2x2y_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_a;
  abcd[16] = 4.0E0*I_ERI_G2xyz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[17] = 4.0E0*I_ERI_G2x2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_a;
  abcd[18] = 4.0E0*I_ERI_G4x_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_D2x_Pz_S_S_C1002_a-2.0E0*3*I_ERI_D2x_Pz_S_S_C1002_a+2*1*I_ERI_S_Pz_S_S_C1002;
  abcd[19] = 4.0E0*I_ERI_G3xy_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[20] = 4.0E0*I_ERI_G3xz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[21] = 4.0E0*I_ERI_G2x2y_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_a;
  abcd[22] = 4.0E0*I_ERI_G2xyz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[23] = 4.0E0*I_ERI_G2x2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[24] = 4.0E0*I_ERI_G3xy_S_S_S_C2_aa-2.0E0*2*I_ERI_Dxy_S_S_S_C2_a;
  abcd[25] = 4.0E0*I_ERI_G2x2y_S_S_S_C2_aa-2.0E0*1*I_ERI_D2x_S_S_S_C2_a-2.0E0*1*I_ERI_D2y_S_S_S_C2_a+1*I_ERI_S_S_S_S_C2;
  abcd[26] = 4.0E0*I_ERI_G2xyz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dyz_S_S_S_C2_a;
  abcd[27] = 4.0E0*I_ERI_Gx3y_S_S_S_C2_aa-2.0E0*2*I_ERI_Dxy_S_S_S_C2_a;
  abcd[28] = 4.0E0*I_ERI_Gx2yz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxz_S_S_S_C2_a;
  abcd[29] = 4.0E0*I_ERI_Gxy2z_S_S_S_C2_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[30] = 4.0E0*I_ERI_G3xy_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[31] = 4.0E0*I_ERI_G2x2y_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_a-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_a+1*I_ERI_S_Px_S_S_C1002;
  abcd[32] = 4.0E0*I_ERI_G2xyz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[33] = 4.0E0*I_ERI_Gx3y_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[34] = 4.0E0*I_ERI_Gx2yz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[35] = 4.0E0*I_ERI_Gxy2z_Px_S_S_C1002_aa;
  abcd[36] = 4.0E0*I_ERI_G3xy_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[37] = 4.0E0*I_ERI_G2x2y_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_a-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_a+1*I_ERI_S_Py_S_S_C1002;
  abcd[38] = 4.0E0*I_ERI_G2xyz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[39] = 4.0E0*I_ERI_Gx3y_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[40] = 4.0E0*I_ERI_Gx2yz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[41] = 4.0E0*I_ERI_Gxy2z_Py_S_S_C1002_aa;
  abcd[42] = 4.0E0*I_ERI_G3xy_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[43] = 4.0E0*I_ERI_G2x2y_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_a-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_a+1*I_ERI_S_Pz_S_S_C1002;
  abcd[44] = 4.0E0*I_ERI_G2xyz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[45] = 4.0E0*I_ERI_Gx3y_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[46] = 4.0E0*I_ERI_Gx2yz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[47] = 4.0E0*I_ERI_Gxy2z_Pz_S_S_C1002_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[48] = 4.0E0*I_ERI_G3xz_S_S_S_C2_aa-2.0E0*2*I_ERI_Dxz_S_S_S_C2_a;
  abcd[49] = 4.0E0*I_ERI_G2xyz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dyz_S_S_S_C2_a;
  abcd[50] = 4.0E0*I_ERI_G2x2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2x_S_S_S_C2_a-2.0E0*1*I_ERI_D2z_S_S_S_C2_a+1*I_ERI_S_S_S_S_C2;
  abcd[51] = 4.0E0*I_ERI_Gx2yz_S_S_S_C2_aa;
  abcd[52] = 4.0E0*I_ERI_Gxy2z_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxy_S_S_S_C2_a;
  abcd[53] = 4.0E0*I_ERI_Gx3z_S_S_S_C2_aa-2.0E0*2*I_ERI_Dxz_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[54] = 4.0E0*I_ERI_G3xz_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[55] = 4.0E0*I_ERI_G2xyz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[56] = 4.0E0*I_ERI_G2x2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_a+1*I_ERI_S_Px_S_S_C1002;
  abcd[57] = 4.0E0*I_ERI_Gx2yz_Px_S_S_C1002_aa;
  abcd[58] = 4.0E0*I_ERI_Gxy2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[59] = 4.0E0*I_ERI_Gx3z_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[60] = 4.0E0*I_ERI_G3xz_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[61] = 4.0E0*I_ERI_G2xyz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[62] = 4.0E0*I_ERI_G2x2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_a+1*I_ERI_S_Py_S_S_C1002;
  abcd[63] = 4.0E0*I_ERI_Gx2yz_Py_S_S_C1002_aa;
  abcd[64] = 4.0E0*I_ERI_Gxy2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[65] = 4.0E0*I_ERI_Gx3z_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[66] = 4.0E0*I_ERI_G3xz_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[67] = 4.0E0*I_ERI_G2xyz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[68] = 4.0E0*I_ERI_G2x2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_a+1*I_ERI_S_Pz_S_S_C1002;
  abcd[69] = 4.0E0*I_ERI_Gx2yz_Pz_S_S_C1002_aa;
  abcd[70] = 4.0E0*I_ERI_Gxy2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[71] = 4.0E0*I_ERI_Gx3z_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dxz_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[72] = 4.0E0*I_ERI_G2x2y_S_S_S_C2_aa-2.0E0*1*I_ERI_D2x_S_S_S_C2_a;
  abcd[73] = 4.0E0*I_ERI_Gx3y_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxy_S_S_S_C2_a-2.0E0*2*I_ERI_Dxy_S_S_S_C2_a;
  abcd[74] = 4.0E0*I_ERI_Gx2yz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxz_S_S_S_C2_a;
  abcd[75] = 4.0E0*I_ERI_G4y_S_S_S_C2_aa-2.0E0*2*I_ERI_D2y_S_S_S_C2_a-2.0E0*3*I_ERI_D2y_S_S_S_C2_a+2*1*I_ERI_S_S_S_S_C2;
  abcd[76] = 4.0E0*I_ERI_G3yz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dyz_S_S_S_C2_a-2.0E0*2*I_ERI_Dyz_S_S_S_C2_a;
  abcd[77] = 4.0E0*I_ERI_G2y2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2z_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[78] = 4.0E0*I_ERI_G2x2y_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_a;
  abcd[79] = 4.0E0*I_ERI_Gx3y_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[80] = 4.0E0*I_ERI_Gx2yz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[81] = 4.0E0*I_ERI_G4y_Px_S_S_C1002_aa-2.0E0*2*I_ERI_D2y_Px_S_S_C1002_a-2.0E0*3*I_ERI_D2y_Px_S_S_C1002_a+2*1*I_ERI_S_Px_S_S_C1002;
  abcd[82] = 4.0E0*I_ERI_G3yz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[83] = 4.0E0*I_ERI_G2y2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_a;
  abcd[84] = 4.0E0*I_ERI_G2x2y_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_a;
  abcd[85] = 4.0E0*I_ERI_Gx3y_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[86] = 4.0E0*I_ERI_Gx2yz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[87] = 4.0E0*I_ERI_G4y_Py_S_S_C1002_aa-2.0E0*2*I_ERI_D2y_Py_S_S_C1002_a-2.0E0*3*I_ERI_D2y_Py_S_S_C1002_a+2*1*I_ERI_S_Py_S_S_C1002;
  abcd[88] = 4.0E0*I_ERI_G3yz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[89] = 4.0E0*I_ERI_G2y2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_a;
  abcd[90] = 4.0E0*I_ERI_G2x2y_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_a;
  abcd[91] = 4.0E0*I_ERI_Gx3y_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[92] = 4.0E0*I_ERI_Gx2yz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[93] = 4.0E0*I_ERI_G4y_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_D2y_Pz_S_S_C1002_a-2.0E0*3*I_ERI_D2y_Pz_S_S_C1002_a+2*1*I_ERI_S_Pz_S_S_C1002;
  abcd[94] = 4.0E0*I_ERI_G3yz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[95] = 4.0E0*I_ERI_G2y2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[96] = 4.0E0*I_ERI_G2xyz_S_S_S_C2_aa;
  abcd[97] = 4.0E0*I_ERI_Gx2yz_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxz_S_S_S_C2_a;
  abcd[98] = 4.0E0*I_ERI_Gxy2z_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxy_S_S_S_C2_a;
  abcd[99] = 4.0E0*I_ERI_G3yz_S_S_S_C2_aa-2.0E0*2*I_ERI_Dyz_S_S_S_C2_a;
  abcd[100] = 4.0E0*I_ERI_G2y2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2y_S_S_S_C2_a-2.0E0*1*I_ERI_D2z_S_S_S_C2_a+1*I_ERI_S_S_S_S_C2;
  abcd[101] = 4.0E0*I_ERI_Gy3z_S_S_S_C2_aa-2.0E0*2*I_ERI_Dyz_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[102] = 4.0E0*I_ERI_G2xyz_Px_S_S_C1002_aa;
  abcd[103] = 4.0E0*I_ERI_Gx2yz_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[104] = 4.0E0*I_ERI_Gxy2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[105] = 4.0E0*I_ERI_G3yz_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[106] = 4.0E0*I_ERI_G2y2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_a+1*I_ERI_S_Px_S_S_C1002;
  abcd[107] = 4.0E0*I_ERI_Gy3z_Px_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[108] = 4.0E0*I_ERI_G2xyz_Py_S_S_C1002_aa;
  abcd[109] = 4.0E0*I_ERI_Gx2yz_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[110] = 4.0E0*I_ERI_Gxy2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[111] = 4.0E0*I_ERI_G3yz_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[112] = 4.0E0*I_ERI_G2y2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_a+1*I_ERI_S_Py_S_S_C1002;
  abcd[113] = 4.0E0*I_ERI_Gy3z_Py_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[114] = 4.0E0*I_ERI_G2xyz_Pz_S_S_C1002_aa;
  abcd[115] = 4.0E0*I_ERI_Gx2yz_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[116] = 4.0E0*I_ERI_Gxy2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[117] = 4.0E0*I_ERI_G3yz_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[118] = 4.0E0*I_ERI_G2y2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_a-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_a+1*I_ERI_S_Pz_S_S_C1002;
  abcd[119] = 4.0E0*I_ERI_Gy3z_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_Dyz_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C2_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C2
   ************************************************************/
  abcd[120] = 4.0E0*I_ERI_G2x2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2x_S_S_S_C2_a;
  abcd[121] = 4.0E0*I_ERI_Gxy2z_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxy_S_S_S_C2_a;
  abcd[122] = 4.0E0*I_ERI_Gx3z_S_S_S_C2_aa-2.0E0*1*I_ERI_Dxz_S_S_S_C2_a-2.0E0*2*I_ERI_Dxz_S_S_S_C2_a;
  abcd[123] = 4.0E0*I_ERI_G2y2z_S_S_S_C2_aa-2.0E0*1*I_ERI_D2y_S_S_S_C2_a;
  abcd[124] = 4.0E0*I_ERI_Gy3z_S_S_S_C2_aa-2.0E0*1*I_ERI_Dyz_S_S_S_C2_a-2.0E0*2*I_ERI_Dyz_S_S_S_C2_a;
  abcd[125] = 4.0E0*I_ERI_G4z_S_S_S_C2_aa-2.0E0*2*I_ERI_D2z_S_S_S_C2_a-2.0E0*3*I_ERI_D2z_S_S_S_C2_a+2*1*I_ERI_S_S_S_S_C2;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_C1002_aa
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1002
   ************************************************************/
  abcd[126] = 4.0E0*I_ERI_G2x2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_a;
  abcd[127] = 4.0E0*I_ERI_Gxy2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_a;
  abcd[128] = 4.0E0*I_ERI_Gx3z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Px_S_S_C1002_a;
  abcd[129] = 4.0E0*I_ERI_G2y2z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_a;
  abcd[130] = 4.0E0*I_ERI_Gy3z_Px_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Px_S_S_C1002_a;
  abcd[131] = 4.0E0*I_ERI_G4z_Px_S_S_C1002_aa-2.0E0*2*I_ERI_D2z_Px_S_S_C1002_a-2.0E0*3*I_ERI_D2z_Px_S_S_C1002_a+2*1*I_ERI_S_Px_S_S_C1002;
  abcd[132] = 4.0E0*I_ERI_G2x2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_a;
  abcd[133] = 4.0E0*I_ERI_Gxy2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_a;
  abcd[134] = 4.0E0*I_ERI_Gx3z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Py_S_S_C1002_a;
  abcd[135] = 4.0E0*I_ERI_G2y2z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_a;
  abcd[136] = 4.0E0*I_ERI_Gy3z_Py_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Py_S_S_C1002_a;
  abcd[137] = 4.0E0*I_ERI_G4z_Py_S_S_C1002_aa-2.0E0*2*I_ERI_D2z_Py_S_S_C1002_a-2.0E0*3*I_ERI_D2z_Py_S_S_C1002_a+2*1*I_ERI_S_Py_S_S_C1002;
  abcd[138] = 4.0E0*I_ERI_G2x2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_a;
  abcd[139] = 4.0E0*I_ERI_Gxy2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_a;
  abcd[140] = 4.0E0*I_ERI_Gx3z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dxz_Pz_S_S_C1002_a;
  abcd[141] = 4.0E0*I_ERI_G2y2z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_a;
  abcd[142] = 4.0E0*I_ERI_Gy3z_Pz_S_S_C1002_aa-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_a-2.0E0*2*I_ERI_Dyz_Pz_S_S_C1002_a;
  abcd[143] = 4.0E0*I_ERI_G4z_Pz_S_S_C1002_aa-2.0E0*2*I_ERI_D2z_Pz_S_S_C1002_a-2.0E0*3*I_ERI_D2z_Pz_S_S_C1002_a+2*1*I_ERI_S_Pz_S_S_C1002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[144] = 4.0E0*I_ERI_F3x_Px_S_S_C2_ab-2.0E0*2*I_ERI_Px_Px_S_S_C2_b;
  abcd[145] = 4.0E0*I_ERI_F2xy_Px_S_S_C2_ab-2.0E0*1*I_ERI_Py_Px_S_S_C2_b;
  abcd[146] = 4.0E0*I_ERI_F2xz_Px_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Px_S_S_C2_b;
  abcd[147] = 4.0E0*I_ERI_Fx2y_Px_S_S_C2_ab;
  abcd[148] = 4.0E0*I_ERI_Fxyz_Px_S_S_C2_ab;
  abcd[149] = 4.0E0*I_ERI_Fx2z_Px_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[150] = 4.0E0*I_ERI_F3x_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F3x_S_S_S_C1002_a-2.0E0*2*I_ERI_Px_D2x_S_S_C1002_b+2*1*I_ERI_Px_S_S_S_C1002;
  abcd[151] = 4.0E0*I_ERI_F2xy_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2x_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[152] = 4.0E0*I_ERI_F2xz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2x_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[153] = 4.0E0*I_ERI_Fx2y_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a;
  abcd[154] = 4.0E0*I_ERI_Fxyz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[155] = 4.0E0*I_ERI_Fx2z_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a;
  abcd[156] = 4.0E0*I_ERI_F3x_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[157] = 4.0E0*I_ERI_F2xy_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[158] = 4.0E0*I_ERI_F2xz_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[159] = 4.0E0*I_ERI_Fx2y_Dxy_S_S_C1002_ab;
  abcd[160] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[161] = 4.0E0*I_ERI_Fx2z_Dxy_S_S_C1002_ab;
  abcd[162] = 4.0E0*I_ERI_F3x_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[163] = 4.0E0*I_ERI_F2xy_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[164] = 4.0E0*I_ERI_F2xz_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxz_S_S_C1002_b;
  abcd[165] = 4.0E0*I_ERI_Fx2y_Dxz_S_S_C1002_ab;
  abcd[166] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[167] = 4.0E0*I_ERI_Fx2z_Dxz_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[168] = 4.0E0*I_ERI_F3x_Py_S_S_C2_ab-2.0E0*2*I_ERI_Px_Py_S_S_C2_b;
  abcd[169] = 4.0E0*I_ERI_F2xy_Py_S_S_C2_ab-2.0E0*1*I_ERI_Py_Py_S_S_C2_b;
  abcd[170] = 4.0E0*I_ERI_F2xz_Py_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Py_S_S_C2_b;
  abcd[171] = 4.0E0*I_ERI_Fx2y_Py_S_S_C2_ab;
  abcd[172] = 4.0E0*I_ERI_Fxyz_Py_S_S_C2_ab;
  abcd[173] = 4.0E0*I_ERI_Fx2z_Py_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[174] = 4.0E0*I_ERI_F3x_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[175] = 4.0E0*I_ERI_F2xy_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[176] = 4.0E0*I_ERI_F2xz_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[177] = 4.0E0*I_ERI_Fx2y_Dxy_S_S_C1002_ab;
  abcd[178] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[179] = 4.0E0*I_ERI_Fx2z_Dxy_S_S_C1002_ab;
  abcd[180] = 4.0E0*I_ERI_F3x_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F3x_S_S_S_C1002_a-2.0E0*2*I_ERI_Px_D2y_S_S_C1002_b+2*1*I_ERI_Px_S_S_S_C1002;
  abcd[181] = 4.0E0*I_ERI_F2xy_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2y_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[182] = 4.0E0*I_ERI_F2xz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2y_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[183] = 4.0E0*I_ERI_Fx2y_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a;
  abcd[184] = 4.0E0*I_ERI_Fxyz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[185] = 4.0E0*I_ERI_Fx2z_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a;
  abcd[186] = 4.0E0*I_ERI_F3x_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[187] = 4.0E0*I_ERI_F2xy_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[188] = 4.0E0*I_ERI_F2xz_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dyz_S_S_C1002_b;
  abcd[189] = 4.0E0*I_ERI_Fx2y_Dyz_S_S_C1002_ab;
  abcd[190] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[191] = 4.0E0*I_ERI_Fx2z_Dyz_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[192] = 4.0E0*I_ERI_F3x_Pz_S_S_C2_ab-2.0E0*2*I_ERI_Px_Pz_S_S_C2_b;
  abcd[193] = 4.0E0*I_ERI_F2xy_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Py_Pz_S_S_C2_b;
  abcd[194] = 4.0E0*I_ERI_F2xz_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Pz_S_S_C2_b;
  abcd[195] = 4.0E0*I_ERI_Fx2y_Pz_S_S_C2_ab;
  abcd[196] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C2_ab;
  abcd[197] = 4.0E0*I_ERI_Fx2z_Pz_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[198] = 4.0E0*I_ERI_F3x_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[199] = 4.0E0*I_ERI_F2xy_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[200] = 4.0E0*I_ERI_F2xz_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxz_S_S_C1002_b;
  abcd[201] = 4.0E0*I_ERI_Fx2y_Dxz_S_S_C1002_ab;
  abcd[202] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[203] = 4.0E0*I_ERI_Fx2z_Dxz_S_S_C1002_ab;
  abcd[204] = 4.0E0*I_ERI_F3x_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[205] = 4.0E0*I_ERI_F2xy_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[206] = 4.0E0*I_ERI_F2xz_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dyz_S_S_C1002_b;
  abcd[207] = 4.0E0*I_ERI_Fx2y_Dyz_S_S_C1002_ab;
  abcd[208] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[209] = 4.0E0*I_ERI_Fx2z_Dyz_S_S_C1002_ab;
  abcd[210] = 4.0E0*I_ERI_F3x_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F3x_S_S_S_C1002_a-2.0E0*2*I_ERI_Px_D2z_S_S_C1002_b+2*1*I_ERI_Px_S_S_S_C1002;
  abcd[211] = 4.0E0*I_ERI_F2xy_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2z_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[212] = 4.0E0*I_ERI_F2xz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2z_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[213] = 4.0E0*I_ERI_Fx2y_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a;
  abcd[214] = 4.0E0*I_ERI_Fxyz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[215] = 4.0E0*I_ERI_Fx2z_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[216] = 4.0E0*I_ERI_F2xy_Px_S_S_C2_ab;
  abcd[217] = 4.0E0*I_ERI_Fx2y_Px_S_S_C2_ab-2.0E0*1*I_ERI_Px_Px_S_S_C2_b;
  abcd[218] = 4.0E0*I_ERI_Fxyz_Px_S_S_C2_ab;
  abcd[219] = 4.0E0*I_ERI_F3y_Px_S_S_C2_ab-2.0E0*2*I_ERI_Py_Px_S_S_C2_b;
  abcd[220] = 4.0E0*I_ERI_F2yz_Px_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Px_S_S_C2_b;
  abcd[221] = 4.0E0*I_ERI_Fy2z_Px_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[222] = 4.0E0*I_ERI_F2xy_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a;
  abcd[223] = 4.0E0*I_ERI_Fx2y_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2x_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[224] = 4.0E0*I_ERI_Fxyz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[225] = 4.0E0*I_ERI_F3y_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F3y_S_S_S_C1002_a-2.0E0*2*I_ERI_Py_D2x_S_S_C1002_b+2*1*I_ERI_Py_S_S_S_C1002;
  abcd[226] = 4.0E0*I_ERI_F2yz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2x_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[227] = 4.0E0*I_ERI_Fy2z_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a;
  abcd[228] = 4.0E0*I_ERI_F2xy_Dxy_S_S_C1002_ab;
  abcd[229] = 4.0E0*I_ERI_Fx2y_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[230] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[231] = 4.0E0*I_ERI_F3y_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[232] = 4.0E0*I_ERI_F2yz_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[233] = 4.0E0*I_ERI_Fy2z_Dxy_S_S_C1002_ab;
  abcd[234] = 4.0E0*I_ERI_F2xy_Dxz_S_S_C1002_ab;
  abcd[235] = 4.0E0*I_ERI_Fx2y_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[236] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[237] = 4.0E0*I_ERI_F3y_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[238] = 4.0E0*I_ERI_F2yz_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxz_S_S_C1002_b;
  abcd[239] = 4.0E0*I_ERI_Fy2z_Dxz_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_F2xy_Py_S_S_C2_ab;
  abcd[241] = 4.0E0*I_ERI_Fx2y_Py_S_S_C2_ab-2.0E0*1*I_ERI_Px_Py_S_S_C2_b;
  abcd[242] = 4.0E0*I_ERI_Fxyz_Py_S_S_C2_ab;
  abcd[243] = 4.0E0*I_ERI_F3y_Py_S_S_C2_ab-2.0E0*2*I_ERI_Py_Py_S_S_C2_b;
  abcd[244] = 4.0E0*I_ERI_F2yz_Py_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Py_S_S_C2_b;
  abcd[245] = 4.0E0*I_ERI_Fy2z_Py_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[246] = 4.0E0*I_ERI_F2xy_Dxy_S_S_C1002_ab;
  abcd[247] = 4.0E0*I_ERI_Fx2y_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[248] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[249] = 4.0E0*I_ERI_F3y_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[250] = 4.0E0*I_ERI_F2yz_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[251] = 4.0E0*I_ERI_Fy2z_Dxy_S_S_C1002_ab;
  abcd[252] = 4.0E0*I_ERI_F2xy_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a;
  abcd[253] = 4.0E0*I_ERI_Fx2y_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2y_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[254] = 4.0E0*I_ERI_Fxyz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[255] = 4.0E0*I_ERI_F3y_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F3y_S_S_S_C1002_a-2.0E0*2*I_ERI_Py_D2y_S_S_C1002_b+2*1*I_ERI_Py_S_S_S_C1002;
  abcd[256] = 4.0E0*I_ERI_F2yz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2y_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[257] = 4.0E0*I_ERI_Fy2z_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a;
  abcd[258] = 4.0E0*I_ERI_F2xy_Dyz_S_S_C1002_ab;
  abcd[259] = 4.0E0*I_ERI_Fx2y_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[260] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[261] = 4.0E0*I_ERI_F3y_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[262] = 4.0E0*I_ERI_F2yz_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dyz_S_S_C1002_b;
  abcd[263] = 4.0E0*I_ERI_Fy2z_Dyz_S_S_C1002_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[264] = 4.0E0*I_ERI_F2xy_Pz_S_S_C2_ab;
  abcd[265] = 4.0E0*I_ERI_Fx2y_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Px_Pz_S_S_C2_b;
  abcd[266] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C2_ab;
  abcd[267] = 4.0E0*I_ERI_F3y_Pz_S_S_C2_ab-2.0E0*2*I_ERI_Py_Pz_S_S_C2_b;
  abcd[268] = 4.0E0*I_ERI_F2yz_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Pz_Pz_S_S_C2_b;
  abcd[269] = 4.0E0*I_ERI_Fy2z_Pz_S_S_C2_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[270] = 4.0E0*I_ERI_F2xy_Dxz_S_S_C1002_ab;
  abcd[271] = 4.0E0*I_ERI_Fx2y_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[272] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[273] = 4.0E0*I_ERI_F3y_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[274] = 4.0E0*I_ERI_F2yz_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dxz_S_S_C1002_b;
  abcd[275] = 4.0E0*I_ERI_Fy2z_Dxz_S_S_C1002_ab;
  abcd[276] = 4.0E0*I_ERI_F2xy_Dyz_S_S_C1002_ab;
  abcd[277] = 4.0E0*I_ERI_Fx2y_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[278] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[279] = 4.0E0*I_ERI_F3y_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[280] = 4.0E0*I_ERI_F2yz_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Pz_Dyz_S_S_C1002_b;
  abcd[281] = 4.0E0*I_ERI_Fy2z_Dyz_S_S_C1002_ab;
  abcd[282] = 4.0E0*I_ERI_F2xy_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2xy_S_S_S_C1002_a;
  abcd[283] = 4.0E0*I_ERI_Fx2y_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2y_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2z_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[284] = 4.0E0*I_ERI_Fxyz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[285] = 4.0E0*I_ERI_F3y_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F3y_S_S_S_C1002_a-2.0E0*2*I_ERI_Py_D2z_S_S_C1002_b+2*1*I_ERI_Py_S_S_S_C1002;
  abcd[286] = 4.0E0*I_ERI_F2yz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a-2.0E0*1*I_ERI_Pz_D2z_S_S_C1002_b+1*I_ERI_Pz_S_S_S_C1002;
  abcd[287] = 4.0E0*I_ERI_Fy2z_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[288] = 4.0E0*I_ERI_F2xz_Px_S_S_C2_ab;
  abcd[289] = 4.0E0*I_ERI_Fxyz_Px_S_S_C2_ab;
  abcd[290] = 4.0E0*I_ERI_Fx2z_Px_S_S_C2_ab-2.0E0*1*I_ERI_Px_Px_S_S_C2_b;
  abcd[291] = 4.0E0*I_ERI_F2yz_Px_S_S_C2_ab;
  abcd[292] = 4.0E0*I_ERI_Fy2z_Px_S_S_C2_ab-2.0E0*1*I_ERI_Py_Px_S_S_C2_b;
  abcd[293] = 4.0E0*I_ERI_F3z_Px_S_S_C2_ab-2.0E0*2*I_ERI_Pz_Px_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[294] = 4.0E0*I_ERI_F2xz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a;
  abcd[295] = 4.0E0*I_ERI_Fxyz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[296] = 4.0E0*I_ERI_Fx2z_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2x_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[297] = 4.0E0*I_ERI_F2yz_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a;
  abcd[298] = 4.0E0*I_ERI_Fy2z_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2x_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[299] = 4.0E0*I_ERI_F3z_D2x_S_S_C1002_ab-2.0E0*1*I_ERI_F3z_S_S_S_C1002_a-2.0E0*2*I_ERI_Pz_D2x_S_S_C1002_b+2*1*I_ERI_Pz_S_S_S_C1002;
  abcd[300] = 4.0E0*I_ERI_F2xz_Dxy_S_S_C1002_ab;
  abcd[301] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[302] = 4.0E0*I_ERI_Fx2z_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[303] = 4.0E0*I_ERI_F2yz_Dxy_S_S_C1002_ab;
  abcd[304] = 4.0E0*I_ERI_Fy2z_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[305] = 4.0E0*I_ERI_F3z_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[306] = 4.0E0*I_ERI_F2xz_Dxz_S_S_C1002_ab;
  abcd[307] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[308] = 4.0E0*I_ERI_Fx2z_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[309] = 4.0E0*I_ERI_F2yz_Dxz_S_S_C1002_ab;
  abcd[310] = 4.0E0*I_ERI_Fy2z_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[311] = 4.0E0*I_ERI_F3z_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dxz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[312] = 4.0E0*I_ERI_F2xz_Py_S_S_C2_ab;
  abcd[313] = 4.0E0*I_ERI_Fxyz_Py_S_S_C2_ab;
  abcd[314] = 4.0E0*I_ERI_Fx2z_Py_S_S_C2_ab-2.0E0*1*I_ERI_Px_Py_S_S_C2_b;
  abcd[315] = 4.0E0*I_ERI_F2yz_Py_S_S_C2_ab;
  abcd[316] = 4.0E0*I_ERI_Fy2z_Py_S_S_C2_ab-2.0E0*1*I_ERI_Py_Py_S_S_C2_b;
  abcd[317] = 4.0E0*I_ERI_F3z_Py_S_S_C2_ab-2.0E0*2*I_ERI_Pz_Py_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[318] = 4.0E0*I_ERI_F2xz_Dxy_S_S_C1002_ab;
  abcd[319] = 4.0E0*I_ERI_Fxyz_Dxy_S_S_C1002_ab;
  abcd[320] = 4.0E0*I_ERI_Fx2z_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxy_S_S_C1002_b;
  abcd[321] = 4.0E0*I_ERI_F2yz_Dxy_S_S_C1002_ab;
  abcd[322] = 4.0E0*I_ERI_Fy2z_Dxy_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxy_S_S_C1002_b;
  abcd[323] = 4.0E0*I_ERI_F3z_Dxy_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dxy_S_S_C1002_b;
  abcd[324] = 4.0E0*I_ERI_F2xz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a;
  abcd[325] = 4.0E0*I_ERI_Fxyz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[326] = 4.0E0*I_ERI_Fx2z_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2y_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[327] = 4.0E0*I_ERI_F2yz_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a;
  abcd[328] = 4.0E0*I_ERI_Fy2z_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2y_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[329] = 4.0E0*I_ERI_F3z_D2y_S_S_C1002_ab-2.0E0*1*I_ERI_F3z_S_S_S_C1002_a-2.0E0*2*I_ERI_Pz_D2y_S_S_C1002_b+2*1*I_ERI_Pz_S_S_S_C1002;
  abcd[330] = 4.0E0*I_ERI_F2xz_Dyz_S_S_C1002_ab;
  abcd[331] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[332] = 4.0E0*I_ERI_Fx2z_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[333] = 4.0E0*I_ERI_F2yz_Dyz_S_S_C1002_ab;
  abcd[334] = 4.0E0*I_ERI_Fy2z_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[335] = 4.0E0*I_ERI_F3z_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dyz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C2_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C2_b
   ************************************************************/
  abcd[336] = 4.0E0*I_ERI_F2xz_Pz_S_S_C2_ab;
  abcd[337] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C2_ab;
  abcd[338] = 4.0E0*I_ERI_Fx2z_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Px_Pz_S_S_C2_b;
  abcd[339] = 4.0E0*I_ERI_F2yz_Pz_S_S_C2_ab;
  abcd[340] = 4.0E0*I_ERI_Fy2z_Pz_S_S_C2_ab-2.0E0*1*I_ERI_Py_Pz_S_S_C2_b;
  abcd[341] = 4.0E0*I_ERI_F3z_Pz_S_S_C2_ab-2.0E0*2*I_ERI_Pz_Pz_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_C1002_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  abcd[342] = 4.0E0*I_ERI_F2xz_Dxz_S_S_C1002_ab;
  abcd[343] = 4.0E0*I_ERI_Fxyz_Dxz_S_S_C1002_ab;
  abcd[344] = 4.0E0*I_ERI_Fx2z_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dxz_S_S_C1002_b;
  abcd[345] = 4.0E0*I_ERI_F2yz_Dxz_S_S_C1002_ab;
  abcd[346] = 4.0E0*I_ERI_Fy2z_Dxz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dxz_S_S_C1002_b;
  abcd[347] = 4.0E0*I_ERI_F3z_Dxz_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dxz_S_S_C1002_b;
  abcd[348] = 4.0E0*I_ERI_F2xz_Dyz_S_S_C1002_ab;
  abcd[349] = 4.0E0*I_ERI_Fxyz_Dyz_S_S_C1002_ab;
  abcd[350] = 4.0E0*I_ERI_Fx2z_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Px_Dyz_S_S_C1002_b;
  abcd[351] = 4.0E0*I_ERI_F2yz_Dyz_S_S_C1002_ab;
  abcd[352] = 4.0E0*I_ERI_Fy2z_Dyz_S_S_C1002_ab-2.0E0*1*I_ERI_Py_Dyz_S_S_C1002_b;
  abcd[353] = 4.0E0*I_ERI_F3z_Dyz_S_S_C1002_ab-2.0E0*2*I_ERI_Pz_Dyz_S_S_C1002_b;
  abcd[354] = 4.0E0*I_ERI_F2xz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2xz_S_S_S_C1002_a;
  abcd[355] = 4.0E0*I_ERI_Fxyz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fxyz_S_S_S_C1002_a;
  abcd[356] = 4.0E0*I_ERI_Fx2z_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fx2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Px_D2z_S_S_C1002_b+1*I_ERI_Px_S_S_S_C1002;
  abcd[357] = 4.0E0*I_ERI_F2yz_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F2yz_S_S_S_C1002_a;
  abcd[358] = 4.0E0*I_ERI_Fy2z_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_Fy2z_S_S_S_C1002_a-2.0E0*1*I_ERI_Py_D2z_S_S_C1002_b+1*I_ERI_Py_S_S_S_C1002;
  abcd[359] = 4.0E0*I_ERI_F3z_D2z_S_S_C1002_ab-2.0E0*1*I_ERI_F3z_S_S_S_C1002_a-2.0E0*2*I_ERI_Pz_D2z_S_S_C1002_b+2*1*I_ERI_Pz_S_S_S_C1002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[360] = 4.0E0*I_ERI_F3x_S_Px_S_C2_ac-2.0E0*2*I_ERI_Px_S_Px_S_C2_c;
  abcd[361] = 4.0E0*I_ERI_F2xy_S_Px_S_C2_ac-2.0E0*1*I_ERI_Py_S_Px_S_C2_c;
  abcd[362] = 4.0E0*I_ERI_F2xz_S_Px_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Px_S_C2_c;
  abcd[363] = 4.0E0*I_ERI_Fx2y_S_Px_S_C2_ac;
  abcd[364] = 4.0E0*I_ERI_Fxyz_S_Px_S_C2_ac;
  abcd[365] = 4.0E0*I_ERI_Fx2z_S_Px_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[366] = 4.0E0*I_ERI_F3x_Px_Px_S_C1002_ac-2.0E0*2*I_ERI_Px_Px_Px_S_C1002_c;
  abcd[367] = 4.0E0*I_ERI_F2xy_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Px_S_C1002_c;
  abcd[368] = 4.0E0*I_ERI_F2xz_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Px_S_C1002_c;
  abcd[369] = 4.0E0*I_ERI_Fx2y_Px_Px_S_C1002_ac;
  abcd[370] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1002_ac;
  abcd[371] = 4.0E0*I_ERI_Fx2z_Px_Px_S_C1002_ac;
  abcd[372] = 4.0E0*I_ERI_F3x_Py_Px_S_C1002_ac-2.0E0*2*I_ERI_Px_Py_Px_S_C1002_c;
  abcd[373] = 4.0E0*I_ERI_F2xy_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Px_S_C1002_c;
  abcd[374] = 4.0E0*I_ERI_F2xz_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Px_S_C1002_c;
  abcd[375] = 4.0E0*I_ERI_Fx2y_Py_Px_S_C1002_ac;
  abcd[376] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1002_ac;
  abcd[377] = 4.0E0*I_ERI_Fx2z_Py_Px_S_C1002_ac;
  abcd[378] = 4.0E0*I_ERI_F3x_Pz_Px_S_C1002_ac-2.0E0*2*I_ERI_Px_Pz_Px_S_C1002_c;
  abcd[379] = 4.0E0*I_ERI_F2xy_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Px_S_C1002_c;
  abcd[380] = 4.0E0*I_ERI_F2xz_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1002_c;
  abcd[381] = 4.0E0*I_ERI_Fx2y_Pz_Px_S_C1002_ac;
  abcd[382] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1002_ac;
  abcd[383] = 4.0E0*I_ERI_Fx2z_Pz_Px_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[384] = 4.0E0*I_ERI_F3x_S_Py_S_C2_ac-2.0E0*2*I_ERI_Px_S_Py_S_C2_c;
  abcd[385] = 4.0E0*I_ERI_F2xy_S_Py_S_C2_ac-2.0E0*1*I_ERI_Py_S_Py_S_C2_c;
  abcd[386] = 4.0E0*I_ERI_F2xz_S_Py_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Py_S_C2_c;
  abcd[387] = 4.0E0*I_ERI_Fx2y_S_Py_S_C2_ac;
  abcd[388] = 4.0E0*I_ERI_Fxyz_S_Py_S_C2_ac;
  abcd[389] = 4.0E0*I_ERI_Fx2z_S_Py_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[390] = 4.0E0*I_ERI_F3x_Px_Py_S_C1002_ac-2.0E0*2*I_ERI_Px_Px_Py_S_C1002_c;
  abcd[391] = 4.0E0*I_ERI_F2xy_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Py_S_C1002_c;
  abcd[392] = 4.0E0*I_ERI_F2xz_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Py_S_C1002_c;
  abcd[393] = 4.0E0*I_ERI_Fx2y_Px_Py_S_C1002_ac;
  abcd[394] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1002_ac;
  abcd[395] = 4.0E0*I_ERI_Fx2z_Px_Py_S_C1002_ac;
  abcd[396] = 4.0E0*I_ERI_F3x_Py_Py_S_C1002_ac-2.0E0*2*I_ERI_Px_Py_Py_S_C1002_c;
  abcd[397] = 4.0E0*I_ERI_F2xy_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Py_S_C1002_c;
  abcd[398] = 4.0E0*I_ERI_F2xz_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Py_S_C1002_c;
  abcd[399] = 4.0E0*I_ERI_Fx2y_Py_Py_S_C1002_ac;
  abcd[400] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1002_ac;
  abcd[401] = 4.0E0*I_ERI_Fx2z_Py_Py_S_C1002_ac;
  abcd[402] = 4.0E0*I_ERI_F3x_Pz_Py_S_C1002_ac-2.0E0*2*I_ERI_Px_Pz_Py_S_C1002_c;
  abcd[403] = 4.0E0*I_ERI_F2xy_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Py_S_C1002_c;
  abcd[404] = 4.0E0*I_ERI_F2xz_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1002_c;
  abcd[405] = 4.0E0*I_ERI_Fx2y_Pz_Py_S_C1002_ac;
  abcd[406] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1002_ac;
  abcd[407] = 4.0E0*I_ERI_Fx2z_Pz_Py_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[408] = 4.0E0*I_ERI_F3x_S_Pz_S_C2_ac-2.0E0*2*I_ERI_Px_S_Pz_S_C2_c;
  abcd[409] = 4.0E0*I_ERI_F2xy_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Py_S_Pz_S_C2_c;
  abcd[410] = 4.0E0*I_ERI_F2xz_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Pz_S_C2_c;
  abcd[411] = 4.0E0*I_ERI_Fx2y_S_Pz_S_C2_ac;
  abcd[412] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C2_ac;
  abcd[413] = 4.0E0*I_ERI_Fx2z_S_Pz_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[414] = 4.0E0*I_ERI_F3x_Px_Pz_S_C1002_ac-2.0E0*2*I_ERI_Px_Px_Pz_S_C1002_c;
  abcd[415] = 4.0E0*I_ERI_F2xy_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Pz_S_C1002_c;
  abcd[416] = 4.0E0*I_ERI_F2xz_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1002_c;
  abcd[417] = 4.0E0*I_ERI_Fx2y_Px_Pz_S_C1002_ac;
  abcd[418] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1002_ac;
  abcd[419] = 4.0E0*I_ERI_Fx2z_Px_Pz_S_C1002_ac;
  abcd[420] = 4.0E0*I_ERI_F3x_Py_Pz_S_C1002_ac-2.0E0*2*I_ERI_Px_Py_Pz_S_C1002_c;
  abcd[421] = 4.0E0*I_ERI_F2xy_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Pz_S_C1002_c;
  abcd[422] = 4.0E0*I_ERI_F2xz_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1002_c;
  abcd[423] = 4.0E0*I_ERI_Fx2y_Py_Pz_S_C1002_ac;
  abcd[424] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1002_ac;
  abcd[425] = 4.0E0*I_ERI_Fx2z_Py_Pz_S_C1002_ac;
  abcd[426] = 4.0E0*I_ERI_F3x_Pz_Pz_S_C1002_ac-2.0E0*2*I_ERI_Px_Pz_Pz_S_C1002_c;
  abcd[427] = 4.0E0*I_ERI_F2xy_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1002_c;
  abcd[428] = 4.0E0*I_ERI_F2xz_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1002_c;
  abcd[429] = 4.0E0*I_ERI_Fx2y_Pz_Pz_S_C1002_ac;
  abcd[430] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1002_ac;
  abcd[431] = 4.0E0*I_ERI_Fx2z_Pz_Pz_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[432] = 4.0E0*I_ERI_F2xy_S_Px_S_C2_ac;
  abcd[433] = 4.0E0*I_ERI_Fx2y_S_Px_S_C2_ac-2.0E0*1*I_ERI_Px_S_Px_S_C2_c;
  abcd[434] = 4.0E0*I_ERI_Fxyz_S_Px_S_C2_ac;
  abcd[435] = 4.0E0*I_ERI_F3y_S_Px_S_C2_ac-2.0E0*2*I_ERI_Py_S_Px_S_C2_c;
  abcd[436] = 4.0E0*I_ERI_F2yz_S_Px_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Px_S_C2_c;
  abcd[437] = 4.0E0*I_ERI_Fy2z_S_Px_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[438] = 4.0E0*I_ERI_F2xy_Px_Px_S_C1002_ac;
  abcd[439] = 4.0E0*I_ERI_Fx2y_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Px_S_C1002_c;
  abcd[440] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1002_ac;
  abcd[441] = 4.0E0*I_ERI_F3y_Px_Px_S_C1002_ac-2.0E0*2*I_ERI_Py_Px_Px_S_C1002_c;
  abcd[442] = 4.0E0*I_ERI_F2yz_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Px_S_C1002_c;
  abcd[443] = 4.0E0*I_ERI_Fy2z_Px_Px_S_C1002_ac;
  abcd[444] = 4.0E0*I_ERI_F2xy_Py_Px_S_C1002_ac;
  abcd[445] = 4.0E0*I_ERI_Fx2y_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Px_S_C1002_c;
  abcd[446] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1002_ac;
  abcd[447] = 4.0E0*I_ERI_F3y_Py_Px_S_C1002_ac-2.0E0*2*I_ERI_Py_Py_Px_S_C1002_c;
  abcd[448] = 4.0E0*I_ERI_F2yz_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Px_S_C1002_c;
  abcd[449] = 4.0E0*I_ERI_Fy2z_Py_Px_S_C1002_ac;
  abcd[450] = 4.0E0*I_ERI_F2xy_Pz_Px_S_C1002_ac;
  abcd[451] = 4.0E0*I_ERI_Fx2y_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Px_S_C1002_c;
  abcd[452] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1002_ac;
  abcd[453] = 4.0E0*I_ERI_F3y_Pz_Px_S_C1002_ac-2.0E0*2*I_ERI_Py_Pz_Px_S_C1002_c;
  abcd[454] = 4.0E0*I_ERI_F2yz_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1002_c;
  abcd[455] = 4.0E0*I_ERI_Fy2z_Pz_Px_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[456] = 4.0E0*I_ERI_F2xy_S_Py_S_C2_ac;
  abcd[457] = 4.0E0*I_ERI_Fx2y_S_Py_S_C2_ac-2.0E0*1*I_ERI_Px_S_Py_S_C2_c;
  abcd[458] = 4.0E0*I_ERI_Fxyz_S_Py_S_C2_ac;
  abcd[459] = 4.0E0*I_ERI_F3y_S_Py_S_C2_ac-2.0E0*2*I_ERI_Py_S_Py_S_C2_c;
  abcd[460] = 4.0E0*I_ERI_F2yz_S_Py_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Py_S_C2_c;
  abcd[461] = 4.0E0*I_ERI_Fy2z_S_Py_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[462] = 4.0E0*I_ERI_F2xy_Px_Py_S_C1002_ac;
  abcd[463] = 4.0E0*I_ERI_Fx2y_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Py_S_C1002_c;
  abcd[464] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1002_ac;
  abcd[465] = 4.0E0*I_ERI_F3y_Px_Py_S_C1002_ac-2.0E0*2*I_ERI_Py_Px_Py_S_C1002_c;
  abcd[466] = 4.0E0*I_ERI_F2yz_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Py_S_C1002_c;
  abcd[467] = 4.0E0*I_ERI_Fy2z_Px_Py_S_C1002_ac;
  abcd[468] = 4.0E0*I_ERI_F2xy_Py_Py_S_C1002_ac;
  abcd[469] = 4.0E0*I_ERI_Fx2y_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Py_S_C1002_c;
  abcd[470] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1002_ac;
  abcd[471] = 4.0E0*I_ERI_F3y_Py_Py_S_C1002_ac-2.0E0*2*I_ERI_Py_Py_Py_S_C1002_c;
  abcd[472] = 4.0E0*I_ERI_F2yz_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Py_S_C1002_c;
  abcd[473] = 4.0E0*I_ERI_Fy2z_Py_Py_S_C1002_ac;
  abcd[474] = 4.0E0*I_ERI_F2xy_Pz_Py_S_C1002_ac;
  abcd[475] = 4.0E0*I_ERI_Fx2y_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Py_S_C1002_c;
  abcd[476] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1002_ac;
  abcd[477] = 4.0E0*I_ERI_F3y_Pz_Py_S_C1002_ac-2.0E0*2*I_ERI_Py_Pz_Py_S_C1002_c;
  abcd[478] = 4.0E0*I_ERI_F2yz_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1002_c;
  abcd[479] = 4.0E0*I_ERI_Fy2z_Pz_Py_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[480] = 4.0E0*I_ERI_F2xy_S_Pz_S_C2_ac;
  abcd[481] = 4.0E0*I_ERI_Fx2y_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Px_S_Pz_S_C2_c;
  abcd[482] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C2_ac;
  abcd[483] = 4.0E0*I_ERI_F3y_S_Pz_S_C2_ac-2.0E0*2*I_ERI_Py_S_Pz_S_C2_c;
  abcd[484] = 4.0E0*I_ERI_F2yz_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Pz_S_Pz_S_C2_c;
  abcd[485] = 4.0E0*I_ERI_Fy2z_S_Pz_S_C2_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[486] = 4.0E0*I_ERI_F2xy_Px_Pz_S_C1002_ac;
  abcd[487] = 4.0E0*I_ERI_Fx2y_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Pz_S_C1002_c;
  abcd[488] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1002_ac;
  abcd[489] = 4.0E0*I_ERI_F3y_Px_Pz_S_C1002_ac-2.0E0*2*I_ERI_Py_Px_Pz_S_C1002_c;
  abcd[490] = 4.0E0*I_ERI_F2yz_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1002_c;
  abcd[491] = 4.0E0*I_ERI_Fy2z_Px_Pz_S_C1002_ac;
  abcd[492] = 4.0E0*I_ERI_F2xy_Py_Pz_S_C1002_ac;
  abcd[493] = 4.0E0*I_ERI_Fx2y_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Pz_S_C1002_c;
  abcd[494] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1002_ac;
  abcd[495] = 4.0E0*I_ERI_F3y_Py_Pz_S_C1002_ac-2.0E0*2*I_ERI_Py_Py_Pz_S_C1002_c;
  abcd[496] = 4.0E0*I_ERI_F2yz_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1002_c;
  abcd[497] = 4.0E0*I_ERI_Fy2z_Py_Pz_S_C1002_ac;
  abcd[498] = 4.0E0*I_ERI_F2xy_Pz_Pz_S_C1002_ac;
  abcd[499] = 4.0E0*I_ERI_Fx2y_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1002_c;
  abcd[500] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1002_ac;
  abcd[501] = 4.0E0*I_ERI_F3y_Pz_Pz_S_C1002_ac-2.0E0*2*I_ERI_Py_Pz_Pz_S_C1002_c;
  abcd[502] = 4.0E0*I_ERI_F2yz_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1002_c;
  abcd[503] = 4.0E0*I_ERI_Fy2z_Pz_Pz_S_C1002_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[504] = 4.0E0*I_ERI_F2xz_S_Px_S_C2_ac;
  abcd[505] = 4.0E0*I_ERI_Fxyz_S_Px_S_C2_ac;
  abcd[506] = 4.0E0*I_ERI_Fx2z_S_Px_S_C2_ac-2.0E0*1*I_ERI_Px_S_Px_S_C2_c;
  abcd[507] = 4.0E0*I_ERI_F2yz_S_Px_S_C2_ac;
  abcd[508] = 4.0E0*I_ERI_Fy2z_S_Px_S_C2_ac-2.0E0*1*I_ERI_Py_S_Px_S_C2_c;
  abcd[509] = 4.0E0*I_ERI_F3z_S_Px_S_C2_ac-2.0E0*2*I_ERI_Pz_S_Px_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[510] = 4.0E0*I_ERI_F2xz_Px_Px_S_C1002_ac;
  abcd[511] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1002_ac;
  abcd[512] = 4.0E0*I_ERI_Fx2z_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Px_S_C1002_c;
  abcd[513] = 4.0E0*I_ERI_F2yz_Px_Px_S_C1002_ac;
  abcd[514] = 4.0E0*I_ERI_Fy2z_Px_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Px_S_C1002_c;
  abcd[515] = 4.0E0*I_ERI_F3z_Px_Px_S_C1002_ac-2.0E0*2*I_ERI_Pz_Px_Px_S_C1002_c;
  abcd[516] = 4.0E0*I_ERI_F2xz_Py_Px_S_C1002_ac;
  abcd[517] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1002_ac;
  abcd[518] = 4.0E0*I_ERI_Fx2z_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Px_S_C1002_c;
  abcd[519] = 4.0E0*I_ERI_F2yz_Py_Px_S_C1002_ac;
  abcd[520] = 4.0E0*I_ERI_Fy2z_Py_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Px_S_C1002_c;
  abcd[521] = 4.0E0*I_ERI_F3z_Py_Px_S_C1002_ac-2.0E0*2*I_ERI_Pz_Py_Px_S_C1002_c;
  abcd[522] = 4.0E0*I_ERI_F2xz_Pz_Px_S_C1002_ac;
  abcd[523] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1002_ac;
  abcd[524] = 4.0E0*I_ERI_Fx2z_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Px_S_C1002_c;
  abcd[525] = 4.0E0*I_ERI_F2yz_Pz_Px_S_C1002_ac;
  abcd[526] = 4.0E0*I_ERI_Fy2z_Pz_Px_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Px_S_C1002_c;
  abcd[527] = 4.0E0*I_ERI_F3z_Pz_Px_S_C1002_ac-2.0E0*2*I_ERI_Pz_Pz_Px_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[528] = 4.0E0*I_ERI_F2xz_S_Py_S_C2_ac;
  abcd[529] = 4.0E0*I_ERI_Fxyz_S_Py_S_C2_ac;
  abcd[530] = 4.0E0*I_ERI_Fx2z_S_Py_S_C2_ac-2.0E0*1*I_ERI_Px_S_Py_S_C2_c;
  abcd[531] = 4.0E0*I_ERI_F2yz_S_Py_S_C2_ac;
  abcd[532] = 4.0E0*I_ERI_Fy2z_S_Py_S_C2_ac-2.0E0*1*I_ERI_Py_S_Py_S_C2_c;
  abcd[533] = 4.0E0*I_ERI_F3z_S_Py_S_C2_ac-2.0E0*2*I_ERI_Pz_S_Py_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[534] = 4.0E0*I_ERI_F2xz_Px_Py_S_C1002_ac;
  abcd[535] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1002_ac;
  abcd[536] = 4.0E0*I_ERI_Fx2z_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Py_S_C1002_c;
  abcd[537] = 4.0E0*I_ERI_F2yz_Px_Py_S_C1002_ac;
  abcd[538] = 4.0E0*I_ERI_Fy2z_Px_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Py_S_C1002_c;
  abcd[539] = 4.0E0*I_ERI_F3z_Px_Py_S_C1002_ac-2.0E0*2*I_ERI_Pz_Px_Py_S_C1002_c;
  abcd[540] = 4.0E0*I_ERI_F2xz_Py_Py_S_C1002_ac;
  abcd[541] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1002_ac;
  abcd[542] = 4.0E0*I_ERI_Fx2z_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Py_S_C1002_c;
  abcd[543] = 4.0E0*I_ERI_F2yz_Py_Py_S_C1002_ac;
  abcd[544] = 4.0E0*I_ERI_Fy2z_Py_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Py_S_C1002_c;
  abcd[545] = 4.0E0*I_ERI_F3z_Py_Py_S_C1002_ac-2.0E0*2*I_ERI_Pz_Py_Py_S_C1002_c;
  abcd[546] = 4.0E0*I_ERI_F2xz_Pz_Py_S_C1002_ac;
  abcd[547] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1002_ac;
  abcd[548] = 4.0E0*I_ERI_Fx2z_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Py_S_C1002_c;
  abcd[549] = 4.0E0*I_ERI_F2yz_Pz_Py_S_C1002_ac;
  abcd[550] = 4.0E0*I_ERI_Fy2z_Pz_Py_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Py_S_C1002_c;
  abcd[551] = 4.0E0*I_ERI_F3z_Pz_Py_S_C1002_ac-2.0E0*2*I_ERI_Pz_Pz_Py_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C2_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C2_c
   ************************************************************/
  abcd[552] = 4.0E0*I_ERI_F2xz_S_Pz_S_C2_ac;
  abcd[553] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C2_ac;
  abcd[554] = 4.0E0*I_ERI_Fx2z_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Px_S_Pz_S_C2_c;
  abcd[555] = 4.0E0*I_ERI_F2yz_S_Pz_S_C2_ac;
  abcd[556] = 4.0E0*I_ERI_Fy2z_S_Pz_S_C2_ac-2.0E0*1*I_ERI_Py_S_Pz_S_C2_c;
  abcd[557] = 4.0E0*I_ERI_F3z_S_Pz_S_C2_ac-2.0E0*2*I_ERI_Pz_S_Pz_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1002_ac
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1002_c
   ************************************************************/
  abcd[558] = 4.0E0*I_ERI_F2xz_Px_Pz_S_C1002_ac;
  abcd[559] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1002_ac;
  abcd[560] = 4.0E0*I_ERI_Fx2z_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Px_Pz_S_C1002_c;
  abcd[561] = 4.0E0*I_ERI_F2yz_Px_Pz_S_C1002_ac;
  abcd[562] = 4.0E0*I_ERI_Fy2z_Px_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Px_Pz_S_C1002_c;
  abcd[563] = 4.0E0*I_ERI_F3z_Px_Pz_S_C1002_ac-2.0E0*2*I_ERI_Pz_Px_Pz_S_C1002_c;
  abcd[564] = 4.0E0*I_ERI_F2xz_Py_Pz_S_C1002_ac;
  abcd[565] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1002_ac;
  abcd[566] = 4.0E0*I_ERI_Fx2z_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Py_Pz_S_C1002_c;
  abcd[567] = 4.0E0*I_ERI_F2yz_Py_Pz_S_C1002_ac;
  abcd[568] = 4.0E0*I_ERI_Fy2z_Py_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Py_Pz_S_C1002_c;
  abcd[569] = 4.0E0*I_ERI_F3z_Py_Pz_S_C1002_ac-2.0E0*2*I_ERI_Pz_Py_Pz_S_C1002_c;
  abcd[570] = 4.0E0*I_ERI_F2xz_Pz_Pz_S_C1002_ac;
  abcd[571] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1002_ac;
  abcd[572] = 4.0E0*I_ERI_Fx2z_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1002_c;
  abcd[573] = 4.0E0*I_ERI_F2yz_Pz_Pz_S_C1002_ac;
  abcd[574] = 4.0E0*I_ERI_Fy2z_Pz_Pz_S_C1002_ac-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1002_c;
  abcd[575] = 4.0E0*I_ERI_F3z_Pz_Pz_S_C1002_ac-2.0E0*2*I_ERI_Pz_Pz_Pz_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[576] = 4.0E0*I_ERI_D2x_D2x_S_S_C2_bb-2.0E0*1*I_ERI_D2x_S_S_S_C2_b;
  abcd[577] = 4.0E0*I_ERI_Dxy_D2x_S_S_C2_bb-2.0E0*1*I_ERI_Dxy_S_S_S_C2_b;
  abcd[578] = 4.0E0*I_ERI_Dxz_D2x_S_S_C2_bb-2.0E0*1*I_ERI_Dxz_S_S_S_C2_b;
  abcd[579] = 4.0E0*I_ERI_D2y_D2x_S_S_C2_bb-2.0E0*1*I_ERI_D2y_S_S_S_C2_b;
  abcd[580] = 4.0E0*I_ERI_Dyz_D2x_S_S_C2_bb-2.0E0*1*I_ERI_Dyz_S_S_S_C2_b;
  abcd[581] = 4.0E0*I_ERI_D2z_D2x_S_S_C2_bb-2.0E0*1*I_ERI_D2z_S_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[582] = 4.0E0*I_ERI_D2x_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_b-2.0E0*2*I_ERI_D2x_Px_S_S_C1002_b;
  abcd[583] = 4.0E0*I_ERI_Dxy_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_b-2.0E0*2*I_ERI_Dxy_Px_S_S_C1002_b;
  abcd[584] = 4.0E0*I_ERI_Dxz_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_b-2.0E0*2*I_ERI_Dxz_Px_S_S_C1002_b;
  abcd[585] = 4.0E0*I_ERI_D2y_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_b-2.0E0*2*I_ERI_D2y_Px_S_S_C1002_b;
  abcd[586] = 4.0E0*I_ERI_Dyz_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_b-2.0E0*2*I_ERI_Dyz_Px_S_S_C1002_b;
  abcd[587] = 4.0E0*I_ERI_D2z_F3x_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_b-2.0E0*2*I_ERI_D2z_Px_S_S_C1002_b;
  abcd[588] = 4.0E0*I_ERI_D2x_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_b;
  abcd[589] = 4.0E0*I_ERI_Dxy_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_b;
  abcd[590] = 4.0E0*I_ERI_Dxz_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_b;
  abcd[591] = 4.0E0*I_ERI_D2y_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_b;
  abcd[592] = 4.0E0*I_ERI_Dyz_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_b;
  abcd[593] = 4.0E0*I_ERI_D2z_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_b;
  abcd[594] = 4.0E0*I_ERI_D2x_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_b;
  abcd[595] = 4.0E0*I_ERI_Dxy_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_b;
  abcd[596] = 4.0E0*I_ERI_Dxz_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_b;
  abcd[597] = 4.0E0*I_ERI_D2y_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_b;
  abcd[598] = 4.0E0*I_ERI_Dyz_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_b;
  abcd[599] = 4.0E0*I_ERI_D2z_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[600] = 4.0E0*I_ERI_D2x_Dxy_S_S_C2_bb;
  abcd[601] = 4.0E0*I_ERI_Dxy_Dxy_S_S_C2_bb;
  abcd[602] = 4.0E0*I_ERI_Dxz_Dxy_S_S_C2_bb;
  abcd[603] = 4.0E0*I_ERI_D2y_Dxy_S_S_C2_bb;
  abcd[604] = 4.0E0*I_ERI_Dyz_Dxy_S_S_C2_bb;
  abcd[605] = 4.0E0*I_ERI_D2z_Dxy_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[606] = 4.0E0*I_ERI_D2x_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_b;
  abcd[607] = 4.0E0*I_ERI_Dxy_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_b;
  abcd[608] = 4.0E0*I_ERI_Dxz_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_b;
  abcd[609] = 4.0E0*I_ERI_D2y_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_b;
  abcd[610] = 4.0E0*I_ERI_Dyz_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_b;
  abcd[611] = 4.0E0*I_ERI_D2z_F2xy_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_b;
  abcd[612] = 4.0E0*I_ERI_D2x_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_b;
  abcd[613] = 4.0E0*I_ERI_Dxy_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_b;
  abcd[614] = 4.0E0*I_ERI_Dxz_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_b;
  abcd[615] = 4.0E0*I_ERI_D2y_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_b;
  abcd[616] = 4.0E0*I_ERI_Dyz_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_b;
  abcd[617] = 4.0E0*I_ERI_D2z_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_b;
  abcd[618] = 4.0E0*I_ERI_D2x_Fxyz_S_S_C1002_bb;
  abcd[619] = 4.0E0*I_ERI_Dxy_Fxyz_S_S_C1002_bb;
  abcd[620] = 4.0E0*I_ERI_Dxz_Fxyz_S_S_C1002_bb;
  abcd[621] = 4.0E0*I_ERI_D2y_Fxyz_S_S_C1002_bb;
  abcd[622] = 4.0E0*I_ERI_Dyz_Fxyz_S_S_C1002_bb;
  abcd[623] = 4.0E0*I_ERI_D2z_Fxyz_S_S_C1002_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[624] = 4.0E0*I_ERI_D2x_Dxz_S_S_C2_bb;
  abcd[625] = 4.0E0*I_ERI_Dxy_Dxz_S_S_C2_bb;
  abcd[626] = 4.0E0*I_ERI_Dxz_Dxz_S_S_C2_bb;
  abcd[627] = 4.0E0*I_ERI_D2y_Dxz_S_S_C2_bb;
  abcd[628] = 4.0E0*I_ERI_Dyz_Dxz_S_S_C2_bb;
  abcd[629] = 4.0E0*I_ERI_D2z_Dxz_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[630] = 4.0E0*I_ERI_D2x_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_b;
  abcd[631] = 4.0E0*I_ERI_Dxy_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_b;
  abcd[632] = 4.0E0*I_ERI_Dxz_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_b;
  abcd[633] = 4.0E0*I_ERI_D2y_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_b;
  abcd[634] = 4.0E0*I_ERI_Dyz_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_b;
  abcd[635] = 4.0E0*I_ERI_D2z_F2xz_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_b;
  abcd[636] = 4.0E0*I_ERI_D2x_Fxyz_S_S_C1002_bb;
  abcd[637] = 4.0E0*I_ERI_Dxy_Fxyz_S_S_C1002_bb;
  abcd[638] = 4.0E0*I_ERI_Dxz_Fxyz_S_S_C1002_bb;
  abcd[639] = 4.0E0*I_ERI_D2y_Fxyz_S_S_C1002_bb;
  abcd[640] = 4.0E0*I_ERI_Dyz_Fxyz_S_S_C1002_bb;
  abcd[641] = 4.0E0*I_ERI_D2z_Fxyz_S_S_C1002_bb;
  abcd[642] = 4.0E0*I_ERI_D2x_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_b;
  abcd[643] = 4.0E0*I_ERI_Dxy_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_b;
  abcd[644] = 4.0E0*I_ERI_Dxz_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_b;
  abcd[645] = 4.0E0*I_ERI_D2y_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_b;
  abcd[646] = 4.0E0*I_ERI_Dyz_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_b;
  abcd[647] = 4.0E0*I_ERI_D2z_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[648] = 4.0E0*I_ERI_D2x_D2y_S_S_C2_bb-2.0E0*1*I_ERI_D2x_S_S_S_C2_b;
  abcd[649] = 4.0E0*I_ERI_Dxy_D2y_S_S_C2_bb-2.0E0*1*I_ERI_Dxy_S_S_S_C2_b;
  abcd[650] = 4.0E0*I_ERI_Dxz_D2y_S_S_C2_bb-2.0E0*1*I_ERI_Dxz_S_S_S_C2_b;
  abcd[651] = 4.0E0*I_ERI_D2y_D2y_S_S_C2_bb-2.0E0*1*I_ERI_D2y_S_S_S_C2_b;
  abcd[652] = 4.0E0*I_ERI_Dyz_D2y_S_S_C2_bb-2.0E0*1*I_ERI_Dyz_S_S_S_C2_b;
  abcd[653] = 4.0E0*I_ERI_D2z_D2y_S_S_C2_bb-2.0E0*1*I_ERI_D2z_S_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[654] = 4.0E0*I_ERI_D2x_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_b;
  abcd[655] = 4.0E0*I_ERI_Dxy_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_b;
  abcd[656] = 4.0E0*I_ERI_Dxz_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_b;
  abcd[657] = 4.0E0*I_ERI_D2y_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_b;
  abcd[658] = 4.0E0*I_ERI_Dyz_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_b;
  abcd[659] = 4.0E0*I_ERI_D2z_Fx2y_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_b;
  abcd[660] = 4.0E0*I_ERI_D2x_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_b-2.0E0*2*I_ERI_D2x_Py_S_S_C1002_b;
  abcd[661] = 4.0E0*I_ERI_Dxy_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_b-2.0E0*2*I_ERI_Dxy_Py_S_S_C1002_b;
  abcd[662] = 4.0E0*I_ERI_Dxz_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_b-2.0E0*2*I_ERI_Dxz_Py_S_S_C1002_b;
  abcd[663] = 4.0E0*I_ERI_D2y_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_b-2.0E0*2*I_ERI_D2y_Py_S_S_C1002_b;
  abcd[664] = 4.0E0*I_ERI_Dyz_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_b-2.0E0*2*I_ERI_Dyz_Py_S_S_C1002_b;
  abcd[665] = 4.0E0*I_ERI_D2z_F3y_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_b-2.0E0*2*I_ERI_D2z_Py_S_S_C1002_b;
  abcd[666] = 4.0E0*I_ERI_D2x_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_b;
  abcd[667] = 4.0E0*I_ERI_Dxy_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_b;
  abcd[668] = 4.0E0*I_ERI_Dxz_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_b;
  abcd[669] = 4.0E0*I_ERI_D2y_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_b;
  abcd[670] = 4.0E0*I_ERI_Dyz_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_b;
  abcd[671] = 4.0E0*I_ERI_D2z_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[672] = 4.0E0*I_ERI_D2x_Dyz_S_S_C2_bb;
  abcd[673] = 4.0E0*I_ERI_Dxy_Dyz_S_S_C2_bb;
  abcd[674] = 4.0E0*I_ERI_Dxz_Dyz_S_S_C2_bb;
  abcd[675] = 4.0E0*I_ERI_D2y_Dyz_S_S_C2_bb;
  abcd[676] = 4.0E0*I_ERI_Dyz_Dyz_S_S_C2_bb;
  abcd[677] = 4.0E0*I_ERI_D2z_Dyz_S_S_C2_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[678] = 4.0E0*I_ERI_D2x_Fxyz_S_S_C1002_bb;
  abcd[679] = 4.0E0*I_ERI_Dxy_Fxyz_S_S_C1002_bb;
  abcd[680] = 4.0E0*I_ERI_Dxz_Fxyz_S_S_C1002_bb;
  abcd[681] = 4.0E0*I_ERI_D2y_Fxyz_S_S_C1002_bb;
  abcd[682] = 4.0E0*I_ERI_Dyz_Fxyz_S_S_C1002_bb;
  abcd[683] = 4.0E0*I_ERI_D2z_Fxyz_S_S_C1002_bb;
  abcd[684] = 4.0E0*I_ERI_D2x_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_b;
  abcd[685] = 4.0E0*I_ERI_Dxy_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_b;
  abcd[686] = 4.0E0*I_ERI_Dxz_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_b;
  abcd[687] = 4.0E0*I_ERI_D2y_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_b;
  abcd[688] = 4.0E0*I_ERI_Dyz_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_b;
  abcd[689] = 4.0E0*I_ERI_D2z_F2yz_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_b;
  abcd[690] = 4.0E0*I_ERI_D2x_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_b;
  abcd[691] = 4.0E0*I_ERI_Dxy_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_b;
  abcd[692] = 4.0E0*I_ERI_Dxz_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_b;
  abcd[693] = 4.0E0*I_ERI_D2y_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_b;
  abcd[694] = 4.0E0*I_ERI_Dyz_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_b;
  abcd[695] = 4.0E0*I_ERI_D2z_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C2_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  abcd[696] = 4.0E0*I_ERI_D2x_D2z_S_S_C2_bb-2.0E0*1*I_ERI_D2x_S_S_S_C2_b;
  abcd[697] = 4.0E0*I_ERI_Dxy_D2z_S_S_C2_bb-2.0E0*1*I_ERI_Dxy_S_S_S_C2_b;
  abcd[698] = 4.0E0*I_ERI_Dxz_D2z_S_S_C2_bb-2.0E0*1*I_ERI_Dxz_S_S_S_C2_b;
  abcd[699] = 4.0E0*I_ERI_D2y_D2z_S_S_C2_bb-2.0E0*1*I_ERI_D2y_S_S_S_C2_b;
  abcd[700] = 4.0E0*I_ERI_Dyz_D2z_S_S_C2_bb-2.0E0*1*I_ERI_Dyz_S_S_S_C2_b;
  abcd[701] = 4.0E0*I_ERI_D2z_D2z_S_S_C2_bb-2.0E0*1*I_ERI_D2z_S_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_S_S_C1002_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  abcd[702] = 4.0E0*I_ERI_D2x_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_b;
  abcd[703] = 4.0E0*I_ERI_Dxy_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_b;
  abcd[704] = 4.0E0*I_ERI_Dxz_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_b;
  abcd[705] = 4.0E0*I_ERI_D2y_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_b;
  abcd[706] = 4.0E0*I_ERI_Dyz_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_b;
  abcd[707] = 4.0E0*I_ERI_D2z_Fx2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_b;
  abcd[708] = 4.0E0*I_ERI_D2x_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_b;
  abcd[709] = 4.0E0*I_ERI_Dxy_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_b;
  abcd[710] = 4.0E0*I_ERI_Dxz_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_b;
  abcd[711] = 4.0E0*I_ERI_D2y_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_b;
  abcd[712] = 4.0E0*I_ERI_Dyz_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_b;
  abcd[713] = 4.0E0*I_ERI_D2z_Fy2z_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_b;
  abcd[714] = 4.0E0*I_ERI_D2x_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_b-2.0E0*2*I_ERI_D2x_Pz_S_S_C1002_b;
  abcd[715] = 4.0E0*I_ERI_Dxy_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_b-2.0E0*2*I_ERI_Dxy_Pz_S_S_C1002_b;
  abcd[716] = 4.0E0*I_ERI_Dxz_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_b-2.0E0*2*I_ERI_Dxz_Pz_S_S_C1002_b;
  abcd[717] = 4.0E0*I_ERI_D2y_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_b-2.0E0*2*I_ERI_D2y_Pz_S_S_C1002_b;
  abcd[718] = 4.0E0*I_ERI_Dyz_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_b-2.0E0*2*I_ERI_Dyz_Pz_S_S_C1002_b;
  abcd[719] = 4.0E0*I_ERI_D2z_F3z_S_S_C1002_bb-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_b-2.0E0*2*I_ERI_D2z_Pz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[720] = 4.0E0*I_ERI_D2x_Px_Px_S_C2_bc;
  abcd[721] = 4.0E0*I_ERI_Dxy_Px_Px_S_C2_bc;
  abcd[722] = 4.0E0*I_ERI_Dxz_Px_Px_S_C2_bc;
  abcd[723] = 4.0E0*I_ERI_D2y_Px_Px_S_C2_bc;
  abcd[724] = 4.0E0*I_ERI_Dyz_Px_Px_S_C2_bc;
  abcd[725] = 4.0E0*I_ERI_D2z_Px_Px_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[726] = 4.0E0*I_ERI_D2x_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Px_S_C1002_c;
  abcd[727] = 4.0E0*I_ERI_Dxy_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Px_S_C1002_c;
  abcd[728] = 4.0E0*I_ERI_Dxz_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Px_S_C1002_c;
  abcd[729] = 4.0E0*I_ERI_D2y_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Px_S_C1002_c;
  abcd[730] = 4.0E0*I_ERI_Dyz_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Px_S_C1002_c;
  abcd[731] = 4.0E0*I_ERI_D2z_D2x_Px_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Px_S_C1002_c;
  abcd[732] = 4.0E0*I_ERI_D2x_Dxy_Px_S_C1002_bc;
  abcd[733] = 4.0E0*I_ERI_Dxy_Dxy_Px_S_C1002_bc;
  abcd[734] = 4.0E0*I_ERI_Dxz_Dxy_Px_S_C1002_bc;
  abcd[735] = 4.0E0*I_ERI_D2y_Dxy_Px_S_C1002_bc;
  abcd[736] = 4.0E0*I_ERI_Dyz_Dxy_Px_S_C1002_bc;
  abcd[737] = 4.0E0*I_ERI_D2z_Dxy_Px_S_C1002_bc;
  abcd[738] = 4.0E0*I_ERI_D2x_Dxz_Px_S_C1002_bc;
  abcd[739] = 4.0E0*I_ERI_Dxy_Dxz_Px_S_C1002_bc;
  abcd[740] = 4.0E0*I_ERI_Dxz_Dxz_Px_S_C1002_bc;
  abcd[741] = 4.0E0*I_ERI_D2y_Dxz_Px_S_C1002_bc;
  abcd[742] = 4.0E0*I_ERI_Dyz_Dxz_Px_S_C1002_bc;
  abcd[743] = 4.0E0*I_ERI_D2z_Dxz_Px_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[744] = 4.0E0*I_ERI_D2x_Px_Py_S_C2_bc;
  abcd[745] = 4.0E0*I_ERI_Dxy_Px_Py_S_C2_bc;
  abcd[746] = 4.0E0*I_ERI_Dxz_Px_Py_S_C2_bc;
  abcd[747] = 4.0E0*I_ERI_D2y_Px_Py_S_C2_bc;
  abcd[748] = 4.0E0*I_ERI_Dyz_Px_Py_S_C2_bc;
  abcd[749] = 4.0E0*I_ERI_D2z_Px_Py_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[750] = 4.0E0*I_ERI_D2x_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Py_S_C1002_c;
  abcd[751] = 4.0E0*I_ERI_Dxy_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Py_S_C1002_c;
  abcd[752] = 4.0E0*I_ERI_Dxz_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Py_S_C1002_c;
  abcd[753] = 4.0E0*I_ERI_D2y_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Py_S_C1002_c;
  abcd[754] = 4.0E0*I_ERI_Dyz_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Py_S_C1002_c;
  abcd[755] = 4.0E0*I_ERI_D2z_D2x_Py_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Py_S_C1002_c;
  abcd[756] = 4.0E0*I_ERI_D2x_Dxy_Py_S_C1002_bc;
  abcd[757] = 4.0E0*I_ERI_Dxy_Dxy_Py_S_C1002_bc;
  abcd[758] = 4.0E0*I_ERI_Dxz_Dxy_Py_S_C1002_bc;
  abcd[759] = 4.0E0*I_ERI_D2y_Dxy_Py_S_C1002_bc;
  abcd[760] = 4.0E0*I_ERI_Dyz_Dxy_Py_S_C1002_bc;
  abcd[761] = 4.0E0*I_ERI_D2z_Dxy_Py_S_C1002_bc;
  abcd[762] = 4.0E0*I_ERI_D2x_Dxz_Py_S_C1002_bc;
  abcd[763] = 4.0E0*I_ERI_Dxy_Dxz_Py_S_C1002_bc;
  abcd[764] = 4.0E0*I_ERI_Dxz_Dxz_Py_S_C1002_bc;
  abcd[765] = 4.0E0*I_ERI_D2y_Dxz_Py_S_C1002_bc;
  abcd[766] = 4.0E0*I_ERI_Dyz_Dxz_Py_S_C1002_bc;
  abcd[767] = 4.0E0*I_ERI_D2z_Dxz_Py_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[768] = 4.0E0*I_ERI_D2x_Px_Pz_S_C2_bc;
  abcd[769] = 4.0E0*I_ERI_Dxy_Px_Pz_S_C2_bc;
  abcd[770] = 4.0E0*I_ERI_Dxz_Px_Pz_S_C2_bc;
  abcd[771] = 4.0E0*I_ERI_D2y_Px_Pz_S_C2_bc;
  abcd[772] = 4.0E0*I_ERI_Dyz_Px_Pz_S_C2_bc;
  abcd[773] = 4.0E0*I_ERI_D2z_Px_Pz_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[774] = 4.0E0*I_ERI_D2x_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Pz_S_C1002_c;
  abcd[775] = 4.0E0*I_ERI_Dxy_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Pz_S_C1002_c;
  abcd[776] = 4.0E0*I_ERI_Dxz_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Pz_S_C1002_c;
  abcd[777] = 4.0E0*I_ERI_D2y_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Pz_S_C1002_c;
  abcd[778] = 4.0E0*I_ERI_Dyz_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Pz_S_C1002_c;
  abcd[779] = 4.0E0*I_ERI_D2z_D2x_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Pz_S_C1002_c;
  abcd[780] = 4.0E0*I_ERI_D2x_Dxy_Pz_S_C1002_bc;
  abcd[781] = 4.0E0*I_ERI_Dxy_Dxy_Pz_S_C1002_bc;
  abcd[782] = 4.0E0*I_ERI_Dxz_Dxy_Pz_S_C1002_bc;
  abcd[783] = 4.0E0*I_ERI_D2y_Dxy_Pz_S_C1002_bc;
  abcd[784] = 4.0E0*I_ERI_Dyz_Dxy_Pz_S_C1002_bc;
  abcd[785] = 4.0E0*I_ERI_D2z_Dxy_Pz_S_C1002_bc;
  abcd[786] = 4.0E0*I_ERI_D2x_Dxz_Pz_S_C1002_bc;
  abcd[787] = 4.0E0*I_ERI_Dxy_Dxz_Pz_S_C1002_bc;
  abcd[788] = 4.0E0*I_ERI_Dxz_Dxz_Pz_S_C1002_bc;
  abcd[789] = 4.0E0*I_ERI_D2y_Dxz_Pz_S_C1002_bc;
  abcd[790] = 4.0E0*I_ERI_Dyz_Dxz_Pz_S_C1002_bc;
  abcd[791] = 4.0E0*I_ERI_D2z_Dxz_Pz_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[792] = 4.0E0*I_ERI_D2x_Py_Px_S_C2_bc;
  abcd[793] = 4.0E0*I_ERI_Dxy_Py_Px_S_C2_bc;
  abcd[794] = 4.0E0*I_ERI_Dxz_Py_Px_S_C2_bc;
  abcd[795] = 4.0E0*I_ERI_D2y_Py_Px_S_C2_bc;
  abcd[796] = 4.0E0*I_ERI_Dyz_Py_Px_S_C2_bc;
  abcd[797] = 4.0E0*I_ERI_D2z_Py_Px_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[798] = 4.0E0*I_ERI_D2x_Dxy_Px_S_C1002_bc;
  abcd[799] = 4.0E0*I_ERI_Dxy_Dxy_Px_S_C1002_bc;
  abcd[800] = 4.0E0*I_ERI_Dxz_Dxy_Px_S_C1002_bc;
  abcd[801] = 4.0E0*I_ERI_D2y_Dxy_Px_S_C1002_bc;
  abcd[802] = 4.0E0*I_ERI_Dyz_Dxy_Px_S_C1002_bc;
  abcd[803] = 4.0E0*I_ERI_D2z_Dxy_Px_S_C1002_bc;
  abcd[804] = 4.0E0*I_ERI_D2x_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Px_S_C1002_c;
  abcd[805] = 4.0E0*I_ERI_Dxy_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Px_S_C1002_c;
  abcd[806] = 4.0E0*I_ERI_Dxz_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Px_S_C1002_c;
  abcd[807] = 4.0E0*I_ERI_D2y_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Px_S_C1002_c;
  abcd[808] = 4.0E0*I_ERI_Dyz_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Px_S_C1002_c;
  abcd[809] = 4.0E0*I_ERI_D2z_D2y_Px_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Px_S_C1002_c;
  abcd[810] = 4.0E0*I_ERI_D2x_Dyz_Px_S_C1002_bc;
  abcd[811] = 4.0E0*I_ERI_Dxy_Dyz_Px_S_C1002_bc;
  abcd[812] = 4.0E0*I_ERI_Dxz_Dyz_Px_S_C1002_bc;
  abcd[813] = 4.0E0*I_ERI_D2y_Dyz_Px_S_C1002_bc;
  abcd[814] = 4.0E0*I_ERI_Dyz_Dyz_Px_S_C1002_bc;
  abcd[815] = 4.0E0*I_ERI_D2z_Dyz_Px_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[816] = 4.0E0*I_ERI_D2x_Py_Py_S_C2_bc;
  abcd[817] = 4.0E0*I_ERI_Dxy_Py_Py_S_C2_bc;
  abcd[818] = 4.0E0*I_ERI_Dxz_Py_Py_S_C2_bc;
  abcd[819] = 4.0E0*I_ERI_D2y_Py_Py_S_C2_bc;
  abcd[820] = 4.0E0*I_ERI_Dyz_Py_Py_S_C2_bc;
  abcd[821] = 4.0E0*I_ERI_D2z_Py_Py_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[822] = 4.0E0*I_ERI_D2x_Dxy_Py_S_C1002_bc;
  abcd[823] = 4.0E0*I_ERI_Dxy_Dxy_Py_S_C1002_bc;
  abcd[824] = 4.0E0*I_ERI_Dxz_Dxy_Py_S_C1002_bc;
  abcd[825] = 4.0E0*I_ERI_D2y_Dxy_Py_S_C1002_bc;
  abcd[826] = 4.0E0*I_ERI_Dyz_Dxy_Py_S_C1002_bc;
  abcd[827] = 4.0E0*I_ERI_D2z_Dxy_Py_S_C1002_bc;
  abcd[828] = 4.0E0*I_ERI_D2x_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Py_S_C1002_c;
  abcd[829] = 4.0E0*I_ERI_Dxy_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Py_S_C1002_c;
  abcd[830] = 4.0E0*I_ERI_Dxz_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Py_S_C1002_c;
  abcd[831] = 4.0E0*I_ERI_D2y_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Py_S_C1002_c;
  abcd[832] = 4.0E0*I_ERI_Dyz_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Py_S_C1002_c;
  abcd[833] = 4.0E0*I_ERI_D2z_D2y_Py_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Py_S_C1002_c;
  abcd[834] = 4.0E0*I_ERI_D2x_Dyz_Py_S_C1002_bc;
  abcd[835] = 4.0E0*I_ERI_Dxy_Dyz_Py_S_C1002_bc;
  abcd[836] = 4.0E0*I_ERI_Dxz_Dyz_Py_S_C1002_bc;
  abcd[837] = 4.0E0*I_ERI_D2y_Dyz_Py_S_C1002_bc;
  abcd[838] = 4.0E0*I_ERI_Dyz_Dyz_Py_S_C1002_bc;
  abcd[839] = 4.0E0*I_ERI_D2z_Dyz_Py_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[840] = 4.0E0*I_ERI_D2x_Py_Pz_S_C2_bc;
  abcd[841] = 4.0E0*I_ERI_Dxy_Py_Pz_S_C2_bc;
  abcd[842] = 4.0E0*I_ERI_Dxz_Py_Pz_S_C2_bc;
  abcd[843] = 4.0E0*I_ERI_D2y_Py_Pz_S_C2_bc;
  abcd[844] = 4.0E0*I_ERI_Dyz_Py_Pz_S_C2_bc;
  abcd[845] = 4.0E0*I_ERI_D2z_Py_Pz_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[846] = 4.0E0*I_ERI_D2x_Dxy_Pz_S_C1002_bc;
  abcd[847] = 4.0E0*I_ERI_Dxy_Dxy_Pz_S_C1002_bc;
  abcd[848] = 4.0E0*I_ERI_Dxz_Dxy_Pz_S_C1002_bc;
  abcd[849] = 4.0E0*I_ERI_D2y_Dxy_Pz_S_C1002_bc;
  abcd[850] = 4.0E0*I_ERI_Dyz_Dxy_Pz_S_C1002_bc;
  abcd[851] = 4.0E0*I_ERI_D2z_Dxy_Pz_S_C1002_bc;
  abcd[852] = 4.0E0*I_ERI_D2x_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Pz_S_C1002_c;
  abcd[853] = 4.0E0*I_ERI_Dxy_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Pz_S_C1002_c;
  abcd[854] = 4.0E0*I_ERI_Dxz_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Pz_S_C1002_c;
  abcd[855] = 4.0E0*I_ERI_D2y_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Pz_S_C1002_c;
  abcd[856] = 4.0E0*I_ERI_Dyz_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Pz_S_C1002_c;
  abcd[857] = 4.0E0*I_ERI_D2z_D2y_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Pz_S_C1002_c;
  abcd[858] = 4.0E0*I_ERI_D2x_Dyz_Pz_S_C1002_bc;
  abcd[859] = 4.0E0*I_ERI_Dxy_Dyz_Pz_S_C1002_bc;
  abcd[860] = 4.0E0*I_ERI_Dxz_Dyz_Pz_S_C1002_bc;
  abcd[861] = 4.0E0*I_ERI_D2y_Dyz_Pz_S_C1002_bc;
  abcd[862] = 4.0E0*I_ERI_Dyz_Dyz_Pz_S_C1002_bc;
  abcd[863] = 4.0E0*I_ERI_D2z_Dyz_Pz_S_C1002_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[864] = 4.0E0*I_ERI_D2x_Pz_Px_S_C2_bc;
  abcd[865] = 4.0E0*I_ERI_Dxy_Pz_Px_S_C2_bc;
  abcd[866] = 4.0E0*I_ERI_Dxz_Pz_Px_S_C2_bc;
  abcd[867] = 4.0E0*I_ERI_D2y_Pz_Px_S_C2_bc;
  abcd[868] = 4.0E0*I_ERI_Dyz_Pz_Px_S_C2_bc;
  abcd[869] = 4.0E0*I_ERI_D2z_Pz_Px_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[870] = 4.0E0*I_ERI_D2x_Dxz_Px_S_C1002_bc;
  abcd[871] = 4.0E0*I_ERI_Dxy_Dxz_Px_S_C1002_bc;
  abcd[872] = 4.0E0*I_ERI_Dxz_Dxz_Px_S_C1002_bc;
  abcd[873] = 4.0E0*I_ERI_D2y_Dxz_Px_S_C1002_bc;
  abcd[874] = 4.0E0*I_ERI_Dyz_Dxz_Px_S_C1002_bc;
  abcd[875] = 4.0E0*I_ERI_D2z_Dxz_Px_S_C1002_bc;
  abcd[876] = 4.0E0*I_ERI_D2x_Dyz_Px_S_C1002_bc;
  abcd[877] = 4.0E0*I_ERI_Dxy_Dyz_Px_S_C1002_bc;
  abcd[878] = 4.0E0*I_ERI_Dxz_Dyz_Px_S_C1002_bc;
  abcd[879] = 4.0E0*I_ERI_D2y_Dyz_Px_S_C1002_bc;
  abcd[880] = 4.0E0*I_ERI_Dyz_Dyz_Px_S_C1002_bc;
  abcd[881] = 4.0E0*I_ERI_D2z_Dyz_Px_S_C1002_bc;
  abcd[882] = 4.0E0*I_ERI_D2x_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Px_S_C1002_c;
  abcd[883] = 4.0E0*I_ERI_Dxy_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Px_S_C1002_c;
  abcd[884] = 4.0E0*I_ERI_Dxz_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Px_S_C1002_c;
  abcd[885] = 4.0E0*I_ERI_D2y_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Px_S_C1002_c;
  abcd[886] = 4.0E0*I_ERI_Dyz_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Px_S_C1002_c;
  abcd[887] = 4.0E0*I_ERI_D2z_D2z_Px_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Px_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[888] = 4.0E0*I_ERI_D2x_Pz_Py_S_C2_bc;
  abcd[889] = 4.0E0*I_ERI_Dxy_Pz_Py_S_C2_bc;
  abcd[890] = 4.0E0*I_ERI_Dxz_Pz_Py_S_C2_bc;
  abcd[891] = 4.0E0*I_ERI_D2y_Pz_Py_S_C2_bc;
  abcd[892] = 4.0E0*I_ERI_Dyz_Pz_Py_S_C2_bc;
  abcd[893] = 4.0E0*I_ERI_D2z_Pz_Py_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[894] = 4.0E0*I_ERI_D2x_Dxz_Py_S_C1002_bc;
  abcd[895] = 4.0E0*I_ERI_Dxy_Dxz_Py_S_C1002_bc;
  abcd[896] = 4.0E0*I_ERI_Dxz_Dxz_Py_S_C1002_bc;
  abcd[897] = 4.0E0*I_ERI_D2y_Dxz_Py_S_C1002_bc;
  abcd[898] = 4.0E0*I_ERI_Dyz_Dxz_Py_S_C1002_bc;
  abcd[899] = 4.0E0*I_ERI_D2z_Dxz_Py_S_C1002_bc;
  abcd[900] = 4.0E0*I_ERI_D2x_Dyz_Py_S_C1002_bc;
  abcd[901] = 4.0E0*I_ERI_Dxy_Dyz_Py_S_C1002_bc;
  abcd[902] = 4.0E0*I_ERI_Dxz_Dyz_Py_S_C1002_bc;
  abcd[903] = 4.0E0*I_ERI_D2y_Dyz_Py_S_C1002_bc;
  abcd[904] = 4.0E0*I_ERI_Dyz_Dyz_Py_S_C1002_bc;
  abcd[905] = 4.0E0*I_ERI_D2z_Dyz_Py_S_C1002_bc;
  abcd[906] = 4.0E0*I_ERI_D2x_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Py_S_C1002_c;
  abcd[907] = 4.0E0*I_ERI_Dxy_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Py_S_C1002_c;
  abcd[908] = 4.0E0*I_ERI_Dxz_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Py_S_C1002_c;
  abcd[909] = 4.0E0*I_ERI_D2y_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Py_S_C1002_c;
  abcd[910] = 4.0E0*I_ERI_Dyz_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Py_S_C1002_c;
  abcd[911] = 4.0E0*I_ERI_D2z_D2z_Py_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Py_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C2_bc
   ************************************************************/
  abcd[912] = 4.0E0*I_ERI_D2x_Pz_Pz_S_C2_bc;
  abcd[913] = 4.0E0*I_ERI_Dxy_Pz_Pz_S_C2_bc;
  abcd[914] = 4.0E0*I_ERI_Dxz_Pz_Pz_S_C2_bc;
  abcd[915] = 4.0E0*I_ERI_D2y_Pz_Pz_S_C2_bc;
  abcd[916] = 4.0E0*I_ERI_Dyz_Pz_Pz_S_C2_bc;
  abcd[917] = 4.0E0*I_ERI_D2z_Pz_Pz_S_C2_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1002_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  abcd[918] = 4.0E0*I_ERI_D2x_Dxz_Pz_S_C1002_bc;
  abcd[919] = 4.0E0*I_ERI_Dxy_Dxz_Pz_S_C1002_bc;
  abcd[920] = 4.0E0*I_ERI_Dxz_Dxz_Pz_S_C1002_bc;
  abcd[921] = 4.0E0*I_ERI_D2y_Dxz_Pz_S_C1002_bc;
  abcd[922] = 4.0E0*I_ERI_Dyz_Dxz_Pz_S_C1002_bc;
  abcd[923] = 4.0E0*I_ERI_D2z_Dxz_Pz_S_C1002_bc;
  abcd[924] = 4.0E0*I_ERI_D2x_Dyz_Pz_S_C1002_bc;
  abcd[925] = 4.0E0*I_ERI_Dxy_Dyz_Pz_S_C1002_bc;
  abcd[926] = 4.0E0*I_ERI_Dxz_Dyz_Pz_S_C1002_bc;
  abcd[927] = 4.0E0*I_ERI_D2y_Dyz_Pz_S_C1002_bc;
  abcd[928] = 4.0E0*I_ERI_Dyz_Dyz_Pz_S_C1002_bc;
  abcd[929] = 4.0E0*I_ERI_D2z_Dyz_Pz_S_C1002_bc;
  abcd[930] = 4.0E0*I_ERI_D2x_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2x_S_Pz_S_C1002_c;
  abcd[931] = 4.0E0*I_ERI_Dxy_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxy_S_Pz_S_C1002_c;
  abcd[932] = 4.0E0*I_ERI_Dxz_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dxz_S_Pz_S_C1002_c;
  abcd[933] = 4.0E0*I_ERI_D2y_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2y_S_Pz_S_C1002_c;
  abcd[934] = 4.0E0*I_ERI_Dyz_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_Dyz_S_Pz_S_C1002_c;
  abcd[935] = 4.0E0*I_ERI_D2z_D2z_Pz_S_C1002_bc-2.0E0*1*I_ERI_D2z_S_Pz_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[936] = 4.0E0*I_ERI_D2x_S_D2x_S_C2_cc-2.0E0*1*I_ERI_D2x_S_S_S_C2_c;
  abcd[937] = 4.0E0*I_ERI_Dxy_S_D2x_S_C2_cc-2.0E0*1*I_ERI_Dxy_S_S_S_C2_c;
  abcd[938] = 4.0E0*I_ERI_Dxz_S_D2x_S_C2_cc-2.0E0*1*I_ERI_Dxz_S_S_S_C2_c;
  abcd[939] = 4.0E0*I_ERI_D2y_S_D2x_S_C2_cc-2.0E0*1*I_ERI_D2y_S_S_S_C2_c;
  abcd[940] = 4.0E0*I_ERI_Dyz_S_D2x_S_C2_cc-2.0E0*1*I_ERI_Dyz_S_S_S_C2_c;
  abcd[941] = 4.0E0*I_ERI_D2z_S_D2x_S_C2_cc-2.0E0*1*I_ERI_D2z_S_S_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[942] = 4.0E0*I_ERI_D2x_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_c;
  abcd[943] = 4.0E0*I_ERI_Dxy_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_c;
  abcd[944] = 4.0E0*I_ERI_Dxz_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_c;
  abcd[945] = 4.0E0*I_ERI_D2y_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_c;
  abcd[946] = 4.0E0*I_ERI_Dyz_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_c;
  abcd[947] = 4.0E0*I_ERI_D2z_Px_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_c;
  abcd[948] = 4.0E0*I_ERI_D2x_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_c;
  abcd[949] = 4.0E0*I_ERI_Dxy_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_c;
  abcd[950] = 4.0E0*I_ERI_Dxz_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_c;
  abcd[951] = 4.0E0*I_ERI_D2y_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_c;
  abcd[952] = 4.0E0*I_ERI_Dyz_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_c;
  abcd[953] = 4.0E0*I_ERI_D2z_Py_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_c;
  abcd[954] = 4.0E0*I_ERI_D2x_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_c;
  abcd[955] = 4.0E0*I_ERI_Dxy_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_c;
  abcd[956] = 4.0E0*I_ERI_Dxz_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_c;
  abcd[957] = 4.0E0*I_ERI_D2y_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_c;
  abcd[958] = 4.0E0*I_ERI_Dyz_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_c;
  abcd[959] = 4.0E0*I_ERI_D2z_Pz_D2x_S_C1002_cc-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[960] = 4.0E0*I_ERI_D2x_S_Dxy_S_C2_cc;
  abcd[961] = 4.0E0*I_ERI_Dxy_S_Dxy_S_C2_cc;
  abcd[962] = 4.0E0*I_ERI_Dxz_S_Dxy_S_C2_cc;
  abcd[963] = 4.0E0*I_ERI_D2y_S_Dxy_S_C2_cc;
  abcd[964] = 4.0E0*I_ERI_Dyz_S_Dxy_S_C2_cc;
  abcd[965] = 4.0E0*I_ERI_D2z_S_Dxy_S_C2_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[966] = 4.0E0*I_ERI_D2x_Px_Dxy_S_C1002_cc;
  abcd[967] = 4.0E0*I_ERI_Dxy_Px_Dxy_S_C1002_cc;
  abcd[968] = 4.0E0*I_ERI_Dxz_Px_Dxy_S_C1002_cc;
  abcd[969] = 4.0E0*I_ERI_D2y_Px_Dxy_S_C1002_cc;
  abcd[970] = 4.0E0*I_ERI_Dyz_Px_Dxy_S_C1002_cc;
  abcd[971] = 4.0E0*I_ERI_D2z_Px_Dxy_S_C1002_cc;
  abcd[972] = 4.0E0*I_ERI_D2x_Py_Dxy_S_C1002_cc;
  abcd[973] = 4.0E0*I_ERI_Dxy_Py_Dxy_S_C1002_cc;
  abcd[974] = 4.0E0*I_ERI_Dxz_Py_Dxy_S_C1002_cc;
  abcd[975] = 4.0E0*I_ERI_D2y_Py_Dxy_S_C1002_cc;
  abcd[976] = 4.0E0*I_ERI_Dyz_Py_Dxy_S_C1002_cc;
  abcd[977] = 4.0E0*I_ERI_D2z_Py_Dxy_S_C1002_cc;
  abcd[978] = 4.0E0*I_ERI_D2x_Pz_Dxy_S_C1002_cc;
  abcd[979] = 4.0E0*I_ERI_Dxy_Pz_Dxy_S_C1002_cc;
  abcd[980] = 4.0E0*I_ERI_Dxz_Pz_Dxy_S_C1002_cc;
  abcd[981] = 4.0E0*I_ERI_D2y_Pz_Dxy_S_C1002_cc;
  abcd[982] = 4.0E0*I_ERI_Dyz_Pz_Dxy_S_C1002_cc;
  abcd[983] = 4.0E0*I_ERI_D2z_Pz_Dxy_S_C1002_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[984] = 4.0E0*I_ERI_D2x_S_Dxz_S_C2_cc;
  abcd[985] = 4.0E0*I_ERI_Dxy_S_Dxz_S_C2_cc;
  abcd[986] = 4.0E0*I_ERI_Dxz_S_Dxz_S_C2_cc;
  abcd[987] = 4.0E0*I_ERI_D2y_S_Dxz_S_C2_cc;
  abcd[988] = 4.0E0*I_ERI_Dyz_S_Dxz_S_C2_cc;
  abcd[989] = 4.0E0*I_ERI_D2z_S_Dxz_S_C2_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[990] = 4.0E0*I_ERI_D2x_Px_Dxz_S_C1002_cc;
  abcd[991] = 4.0E0*I_ERI_Dxy_Px_Dxz_S_C1002_cc;
  abcd[992] = 4.0E0*I_ERI_Dxz_Px_Dxz_S_C1002_cc;
  abcd[993] = 4.0E0*I_ERI_D2y_Px_Dxz_S_C1002_cc;
  abcd[994] = 4.0E0*I_ERI_Dyz_Px_Dxz_S_C1002_cc;
  abcd[995] = 4.0E0*I_ERI_D2z_Px_Dxz_S_C1002_cc;
  abcd[996] = 4.0E0*I_ERI_D2x_Py_Dxz_S_C1002_cc;
  abcd[997] = 4.0E0*I_ERI_Dxy_Py_Dxz_S_C1002_cc;
  abcd[998] = 4.0E0*I_ERI_Dxz_Py_Dxz_S_C1002_cc;
  abcd[999] = 4.0E0*I_ERI_D2y_Py_Dxz_S_C1002_cc;
  abcd[1000] = 4.0E0*I_ERI_Dyz_Py_Dxz_S_C1002_cc;
  abcd[1001] = 4.0E0*I_ERI_D2z_Py_Dxz_S_C1002_cc;
  abcd[1002] = 4.0E0*I_ERI_D2x_Pz_Dxz_S_C1002_cc;
  abcd[1003] = 4.0E0*I_ERI_Dxy_Pz_Dxz_S_C1002_cc;
  abcd[1004] = 4.0E0*I_ERI_Dxz_Pz_Dxz_S_C1002_cc;
  abcd[1005] = 4.0E0*I_ERI_D2y_Pz_Dxz_S_C1002_cc;
  abcd[1006] = 4.0E0*I_ERI_Dyz_Pz_Dxz_S_C1002_cc;
  abcd[1007] = 4.0E0*I_ERI_D2z_Pz_Dxz_S_C1002_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[1008] = 4.0E0*I_ERI_D2x_S_D2y_S_C2_cc-2.0E0*1*I_ERI_D2x_S_S_S_C2_c;
  abcd[1009] = 4.0E0*I_ERI_Dxy_S_D2y_S_C2_cc-2.0E0*1*I_ERI_Dxy_S_S_S_C2_c;
  abcd[1010] = 4.0E0*I_ERI_Dxz_S_D2y_S_C2_cc-2.0E0*1*I_ERI_Dxz_S_S_S_C2_c;
  abcd[1011] = 4.0E0*I_ERI_D2y_S_D2y_S_C2_cc-2.0E0*1*I_ERI_D2y_S_S_S_C2_c;
  abcd[1012] = 4.0E0*I_ERI_Dyz_S_D2y_S_C2_cc-2.0E0*1*I_ERI_Dyz_S_S_S_C2_c;
  abcd[1013] = 4.0E0*I_ERI_D2z_S_D2y_S_C2_cc-2.0E0*1*I_ERI_D2z_S_S_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[1014] = 4.0E0*I_ERI_D2x_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_c;
  abcd[1015] = 4.0E0*I_ERI_Dxy_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_c;
  abcd[1016] = 4.0E0*I_ERI_Dxz_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_c;
  abcd[1017] = 4.0E0*I_ERI_D2y_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_c;
  abcd[1018] = 4.0E0*I_ERI_Dyz_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_c;
  abcd[1019] = 4.0E0*I_ERI_D2z_Px_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_c;
  abcd[1020] = 4.0E0*I_ERI_D2x_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_c;
  abcd[1021] = 4.0E0*I_ERI_Dxy_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_c;
  abcd[1022] = 4.0E0*I_ERI_Dxz_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_c;
  abcd[1023] = 4.0E0*I_ERI_D2y_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_c;
  abcd[1024] = 4.0E0*I_ERI_Dyz_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_c;
  abcd[1025] = 4.0E0*I_ERI_D2z_Py_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_c;
  abcd[1026] = 4.0E0*I_ERI_D2x_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_c;
  abcd[1027] = 4.0E0*I_ERI_Dxy_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_c;
  abcd[1028] = 4.0E0*I_ERI_Dxz_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_c;
  abcd[1029] = 4.0E0*I_ERI_D2y_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_c;
  abcd[1030] = 4.0E0*I_ERI_Dyz_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_c;
  abcd[1031] = 4.0E0*I_ERI_D2z_Pz_D2y_S_C1002_cc-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[1032] = 4.0E0*I_ERI_D2x_S_Dyz_S_C2_cc;
  abcd[1033] = 4.0E0*I_ERI_Dxy_S_Dyz_S_C2_cc;
  abcd[1034] = 4.0E0*I_ERI_Dxz_S_Dyz_S_C2_cc;
  abcd[1035] = 4.0E0*I_ERI_D2y_S_Dyz_S_C2_cc;
  abcd[1036] = 4.0E0*I_ERI_Dyz_S_Dyz_S_C2_cc;
  abcd[1037] = 4.0E0*I_ERI_D2z_S_Dyz_S_C2_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[1038] = 4.0E0*I_ERI_D2x_Px_Dyz_S_C1002_cc;
  abcd[1039] = 4.0E0*I_ERI_Dxy_Px_Dyz_S_C1002_cc;
  abcd[1040] = 4.0E0*I_ERI_Dxz_Px_Dyz_S_C1002_cc;
  abcd[1041] = 4.0E0*I_ERI_D2y_Px_Dyz_S_C1002_cc;
  abcd[1042] = 4.0E0*I_ERI_Dyz_Px_Dyz_S_C1002_cc;
  abcd[1043] = 4.0E0*I_ERI_D2z_Px_Dyz_S_C1002_cc;
  abcd[1044] = 4.0E0*I_ERI_D2x_Py_Dyz_S_C1002_cc;
  abcd[1045] = 4.0E0*I_ERI_Dxy_Py_Dyz_S_C1002_cc;
  abcd[1046] = 4.0E0*I_ERI_Dxz_Py_Dyz_S_C1002_cc;
  abcd[1047] = 4.0E0*I_ERI_D2y_Py_Dyz_S_C1002_cc;
  abcd[1048] = 4.0E0*I_ERI_Dyz_Py_Dyz_S_C1002_cc;
  abcd[1049] = 4.0E0*I_ERI_D2z_Py_Dyz_S_C1002_cc;
  abcd[1050] = 4.0E0*I_ERI_D2x_Pz_Dyz_S_C1002_cc;
  abcd[1051] = 4.0E0*I_ERI_Dxy_Pz_Dyz_S_C1002_cc;
  abcd[1052] = 4.0E0*I_ERI_Dxz_Pz_Dyz_S_C1002_cc;
  abcd[1053] = 4.0E0*I_ERI_D2y_Pz_Dyz_S_C1002_cc;
  abcd[1054] = 4.0E0*I_ERI_Dyz_Pz_Dyz_S_C1002_cc;
  abcd[1055] = 4.0E0*I_ERI_D2z_Pz_Dyz_S_C1002_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_c
   ************************************************************/
  abcd[1056] = 4.0E0*I_ERI_D2x_S_D2z_S_C2_cc-2.0E0*1*I_ERI_D2x_S_S_S_C2_c;
  abcd[1057] = 4.0E0*I_ERI_Dxy_S_D2z_S_C2_cc-2.0E0*1*I_ERI_Dxy_S_S_S_C2_c;
  abcd[1058] = 4.0E0*I_ERI_Dxz_S_D2z_S_C2_cc-2.0E0*1*I_ERI_Dxz_S_S_S_C2_c;
  abcd[1059] = 4.0E0*I_ERI_D2y_S_D2z_S_C2_cc-2.0E0*1*I_ERI_D2y_S_S_S_C2_c;
  abcd[1060] = 4.0E0*I_ERI_Dyz_S_D2z_S_C2_cc-2.0E0*1*I_ERI_Dyz_S_S_S_C2_c;
  abcd[1061] = 4.0E0*I_ERI_D2z_S_D2z_S_C2_cc-2.0E0*1*I_ERI_D2z_S_S_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1002_cc
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_c
   ************************************************************/
  abcd[1062] = 4.0E0*I_ERI_D2x_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2x_Px_S_S_C1002_c;
  abcd[1063] = 4.0E0*I_ERI_Dxy_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Px_S_S_C1002_c;
  abcd[1064] = 4.0E0*I_ERI_Dxz_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Px_S_S_C1002_c;
  abcd[1065] = 4.0E0*I_ERI_D2y_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2y_Px_S_S_C1002_c;
  abcd[1066] = 4.0E0*I_ERI_Dyz_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Px_S_S_C1002_c;
  abcd[1067] = 4.0E0*I_ERI_D2z_Px_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2z_Px_S_S_C1002_c;
  abcd[1068] = 4.0E0*I_ERI_D2x_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2x_Py_S_S_C1002_c;
  abcd[1069] = 4.0E0*I_ERI_Dxy_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Py_S_S_C1002_c;
  abcd[1070] = 4.0E0*I_ERI_Dxz_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Py_S_S_C1002_c;
  abcd[1071] = 4.0E0*I_ERI_D2y_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2y_Py_S_S_C1002_c;
  abcd[1072] = 4.0E0*I_ERI_Dyz_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Py_S_S_C1002_c;
  abcd[1073] = 4.0E0*I_ERI_D2z_Py_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2z_Py_S_S_C1002_c;
  abcd[1074] = 4.0E0*I_ERI_D2x_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2x_Pz_S_S_C1002_c;
  abcd[1075] = 4.0E0*I_ERI_Dxy_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1002_c;
  abcd[1076] = 4.0E0*I_ERI_Dxz_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1002_c;
  abcd[1077] = 4.0E0*I_ERI_D2y_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2y_Pz_S_S_C1002_c;
  abcd[1078] = 4.0E0*I_ERI_Dyz_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1002_c;
  abcd[1079] = 4.0E0*I_ERI_D2z_Pz_D2z_S_C1002_cc-2.0E0*1*I_ERI_D2z_Pz_S_S_C1002_c;
}
