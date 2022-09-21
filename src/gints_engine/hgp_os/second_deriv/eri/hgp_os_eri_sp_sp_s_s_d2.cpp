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
// BRA1 as redundant position, total RHS integrals evaluated as: 101418
// BRA2 as redundant position, total RHS integrals evaluated as: 101418
// KET1 as redundant position, total RHS integrals evaluated as: 72308
// KET2 as redundant position, total RHS integrals evaluated as: 72308
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

void hgp_os_eri_sp_sp_s_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_D2x_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_D2y_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_D2z_S_S_S_C0_aa = 0.0E0;
  Double I_ERI_S_S_S_S_C0_a = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1_aa = 0.0E0;
  Double I_ERI_Px_S_S_S_C1_a = 0.0E0;
  Double I_ERI_Py_S_S_S_C1_a = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1_a = 0.0E0;
  Double I_ERI_S_Px_S_S_C1000_a = 0.0E0;
  Double I_ERI_S_Py_S_S_C1000_a = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1000_a = 0.0E0;
  Double I_ERI_S_Px_S_S_C1_b = 0.0E0;
  Double I_ERI_S_Py_S_S_C1_b = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1_b = 0.0E0;
  Double I_ERI_Px_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Py_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_S_D2x_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_Dxy_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_Dxz_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_D2y_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_Dyz_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_D2z_S_S_C1001_b = 0.0E0;
  Double I_ERI_S_S_S_S_C1001 = 0.0E0;
  Double I_ERI_Px_S_Px_S_C0_ac = 0.0E0;
  Double I_ERI_Py_S_Px_S_C0_ac = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C0_ac = 0.0E0;
  Double I_ERI_Px_S_Py_S_C0_ac = 0.0E0;
  Double I_ERI_Py_S_Py_S_C0_ac = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C0_ac = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C0_ac = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C0_ac = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C0_ac = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1_ac = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1_ac = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1_ac = 0.0E0;
  Double I_ERI_S_S_Px_S_C1_c = 0.0E0;
  Double I_ERI_S_S_Py_S_C1_c = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1_c = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001_c = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001_c = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001_c = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001_c = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001_c = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001_c = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_S_D2x_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_Dxy_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_Dxz_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_D2y_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_Dyz_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_D2z_S_S_C0_bb = 0.0E0;
  Double I_ERI_S_S_S_S_C0_b = 0.0E0;
  Double I_ERI_Px_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1_b = 0.0E0;
  Double I_ERI_S_F3x_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_F2xy_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_F2xz_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_Fx2y_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_Fxyz_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_Fx2z_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_F3y_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_F2yz_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_Fy2z_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_F3z_S_S_C1000_bb = 0.0E0;
  Double I_ERI_S_Px_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Py_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Px_Px_S_C0_bc = 0.0E0;
  Double I_ERI_S_Py_Px_S_C0_bc = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C0_bc = 0.0E0;
  Double I_ERI_S_Px_Py_S_C0_bc = 0.0E0;
  Double I_ERI_S_Py_Py_S_C0_bc = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C0_bc = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C0_bc = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C0_bc = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C0_bc = 0.0E0;
  Double I_ERI_S_D2x_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxy_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxz_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2y_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dyz_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2z_Px_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2x_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxy_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxz_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2y_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dyz_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2z_Py_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2x_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxy_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dxz_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2y_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_Dyz_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_D2z_Pz_S_C1000_bc = 0.0E0;
  Double I_ERI_S_S_Px_S_C1000_c = 0.0E0;
  Double I_ERI_S_S_Py_S_C1000_c = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1000_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_S_S_D2x_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_D2y_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_D2z_S_C0_cc = 0.0E0;
  Double I_ERI_S_S_S_S_C0_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1_cc = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1_cc = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1_cc = 0.0E0;
  Double I_ERI_Px_S_S_S_C1_c = 0.0E0;
  Double I_ERI_Py_S_S_S_C1_c = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1_c = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1000_cc = 0.0E0;
  Double I_ERI_S_Px_S_S_C1000_c = 0.0E0;
  Double I_ERI_S_Py_S_S_C1000_c = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1000_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1000_aa = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001_aa = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_D2y_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_D2z_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Px_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Py_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_Pz_S_S_S_C0_ab = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1_ab = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000_ac = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000_ac = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000_ac = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001_ac = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001_ac = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001_ac = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1_bc = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1_bc = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1_bc = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001_cc = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001_cc = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Px_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Py_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1000_ab = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_ab = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Px_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Py_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1_bb = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001_bc = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001_bc = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001_bc = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001_bb = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001_bb = 0.0E0;

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
    Double ic2_2 = icoe[ip2+2*inp2];
    Double ic2_3 = icoe[ip2+3*inp2];
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

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER43;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER17*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER15*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER13*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = 1.0E0+u2*ONEOVER11*I_ERI_S_S_S_S_M4_vrr;
        I_ERI_S_S_S_S_M4_vrr = ONEOVER9*I_ERI_S_S_S_S_M4_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M4_vrr  = f*I_ERI_S_S_S_S_M4_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);

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
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
      }

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
       * shell quartet name: SQ_ERI_S_P_S_S_M2
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_Px_S_S_M2_vrr = PBX*I_ERI_S_S_S_S_M2_vrr+WPX*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_Py_S_S_M2_vrr = PBY*I_ERI_S_S_S_S_M2_vrr+WPY*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_Pz_S_S_M2_vrr = PBZ*I_ERI_S_S_S_S_M2_vrr+WPZ*I_ERI_S_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_S_P_S_S_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_Px_S_S_M1_vrr = PBX*I_ERI_S_S_S_S_M1_vrr+WPX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Py_S_S_M1_vrr = PBY*I_ERI_S_S_S_S_M1_vrr+WPY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Pz_S_S_M1_vrr = PBZ*I_ERI_S_S_S_S_M1_vrr+WPZ*I_ERI_S_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_P_P_S_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_Px_Px_S_M1_vrr = PBX*I_ERI_S_S_Px_S_M1_vrr+WPX*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Py_Px_S_M1_vrr = PBY*I_ERI_S_S_Px_S_M1_vrr+WPY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_Pz_Px_S_M1_vrr = PBZ*I_ERI_S_S_Px_S_M1_vrr+WPZ*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_Px_Py_S_M1_vrr = PBX*I_ERI_S_S_Py_S_M1_vrr+WPX*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_Py_Py_S_M1_vrr = PBY*I_ERI_S_S_Py_S_M1_vrr+WPY*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Pz_Py_S_M1_vrr = PBZ*I_ERI_S_S_Py_S_M1_vrr+WPZ*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_Px_Pz_S_M1_vrr = PBX*I_ERI_S_S_Pz_S_M1_vrr+WPX*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_S_Py_Pz_S_M1_vrr = PBY*I_ERI_S_S_Pz_S_M1_vrr+WPY*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_S_Pz_Pz_S_M1_vrr = PBZ*I_ERI_S_S_Pz_S_M1_vrr+WPZ*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_S_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_D_S_S_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_D2x_S_S_M1_vrr = PBX*I_ERI_S_Px_S_S_M1_vrr+WPX*I_ERI_S_Px_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Dxy_S_S_M1_vrr = PBY*I_ERI_S_Px_S_S_M1_vrr+WPY*I_ERI_S_Px_S_S_M2_vrr;
      Double I_ERI_S_D2y_S_S_M1_vrr = PBY*I_ERI_S_Py_S_S_M1_vrr+WPY*I_ERI_S_Py_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_D2z_S_S_M1_vrr = PBZ*I_ERI_S_Pz_S_S_M1_vrr+WPZ*I_ERI_S_Pz_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M1_vrr = PAX*I_ERI_D2x_S_S_S_M1_vrr+WPX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2z*I_ERI_Px_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_S_S_M1_vrr = PAY*I_ERI_D2x_S_S_S_M1_vrr+WPY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_S_S_M1_vrr = PAZ*I_ERI_D2x_S_S_S_M1_vrr+WPZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_S_S_M1_vrr = PAX*I_ERI_D2y_S_S_S_M1_vrr+WPX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_S_S_M1_vrr = PAX*I_ERI_D2z_S_S_S_M1_vrr+WPX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_S_S_M1_vrr = PAY*I_ERI_D2y_S_S_S_M1_vrr+WPY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2z*I_ERI_Py_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_S_S_M1_vrr = PAZ*I_ERI_D2y_S_S_S_M1_vrr+WPZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_S_S_M1_vrr = PAZ*I_ERI_D2z_S_S_S_M1_vrr+WPZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2z*I_ERI_Pz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_P_P_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_Px_Px_S_vrr = PBX*I_ERI_S_S_Px_S_vrr+WPX*I_ERI_S_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Py_Px_S_vrr = PBY*I_ERI_S_S_Px_S_vrr+WPY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Pz_Px_S_vrr = PBZ*I_ERI_S_S_Px_S_vrr+WPZ*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Px_Py_S_vrr = PBX*I_ERI_S_S_Py_S_vrr+WPX*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Py_Py_S_vrr = PBY*I_ERI_S_S_Py_S_vrr+WPY*I_ERI_S_S_Py_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Pz_Py_S_vrr = PBZ*I_ERI_S_S_Py_S_vrr+WPZ*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Px_Pz_S_vrr = PBX*I_ERI_S_S_Pz_S_vrr+WPX*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Py_Pz_S_vrr = PBY*I_ERI_S_S_Pz_S_vrr+WPY*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Pz_Pz_S_vrr = PBZ*I_ERI_S_S_Pz_S_vrr+WPZ*I_ERI_S_S_Pz_S_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dxy_S_vrr = QCY*I_ERI_S_S_Px_S_vrr+WQY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_Dxz_S_vrr = QCZ*I_ERI_S_S_Px_S_vrr+WQZ*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dyz_S_vrr = QCZ*I_ERI_S_S_Py_S_vrr+WQZ*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_S_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_S_S
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_D2x_S_S_vrr = PBX*I_ERI_S_Px_S_S_vrr+WPX*I_ERI_S_Px_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Dxy_S_S_vrr = PBY*I_ERI_S_Px_S_S_vrr+WPY*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Dxz_S_S_vrr = PBZ*I_ERI_S_Px_S_S_vrr+WPZ*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_D2y_S_S_vrr = PBY*I_ERI_S_Py_S_S_vrr+WPY*I_ERI_S_Py_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Dyz_S_S_vrr = PBZ*I_ERI_S_Py_S_S_vrr+WPZ*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_D2z_S_S_vrr = PBZ*I_ERI_S_Pz_S_S_vrr+WPZ*I_ERI_S_Pz_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_P_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_P_S
       * RHS shell quartet name: SQ_ERI_S_P_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_Px_D2x_S_vrr = QCX*I_ERI_S_Px_Px_S_vrr+WQX*I_ERI_S_Px_Px_S_M1_vrr+oned2e*I_ERI_S_Px_S_S_vrr-rhod2esq*I_ERI_S_Px_S_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Py_D2x_S_vrr = QCX*I_ERI_S_Py_Px_S_vrr+WQX*I_ERI_S_Py_Px_S_M1_vrr+oned2e*I_ERI_S_Py_S_S_vrr-rhod2esq*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_Pz_D2x_S_vrr = QCX*I_ERI_S_Pz_Px_S_vrr+WQX*I_ERI_S_Pz_Px_S_M1_vrr+oned2e*I_ERI_S_Pz_S_S_vrr-rhod2esq*I_ERI_S_Pz_S_S_M1_vrr;
      Double I_ERI_S_Px_Dxy_S_vrr = QCY*I_ERI_S_Px_Px_S_vrr+WQY*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_Py_Dxy_S_vrr = QCY*I_ERI_S_Py_Px_S_vrr+WQY*I_ERI_S_Py_Px_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Pz_Dxy_S_vrr = QCY*I_ERI_S_Pz_Px_S_vrr+WQY*I_ERI_S_Pz_Px_S_M1_vrr;
      Double I_ERI_S_Px_Dxz_S_vrr = QCZ*I_ERI_S_Px_Px_S_vrr+WQZ*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_Py_Dxz_S_vrr = QCZ*I_ERI_S_Py_Px_S_vrr+WQZ*I_ERI_S_Py_Px_S_M1_vrr;
      Double I_ERI_S_Pz_Dxz_S_vrr = QCZ*I_ERI_S_Pz_Px_S_vrr+WQZ*I_ERI_S_Pz_Px_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Px_D2y_S_vrr = QCY*I_ERI_S_Px_Py_S_vrr+WQY*I_ERI_S_Px_Py_S_M1_vrr+oned2e*I_ERI_S_Px_S_S_vrr-rhod2esq*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Py_D2y_S_vrr = QCY*I_ERI_S_Py_Py_S_vrr+WQY*I_ERI_S_Py_Py_S_M1_vrr+oned2e*I_ERI_S_Py_S_S_vrr-rhod2esq*I_ERI_S_Py_S_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Pz_D2y_S_vrr = QCY*I_ERI_S_Pz_Py_S_vrr+WQY*I_ERI_S_Pz_Py_S_M1_vrr+oned2e*I_ERI_S_Pz_S_S_vrr-rhod2esq*I_ERI_S_Pz_S_S_M1_vrr;
      Double I_ERI_S_Px_Dyz_S_vrr = QCZ*I_ERI_S_Px_Py_S_vrr+WQZ*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_Py_Dyz_S_vrr = QCZ*I_ERI_S_Py_Py_S_vrr+WQZ*I_ERI_S_Py_Py_S_M1_vrr;
      Double I_ERI_S_Pz_Dyz_S_vrr = QCZ*I_ERI_S_Pz_Py_S_vrr+WQZ*I_ERI_S_Pz_Py_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Px_D2z_S_vrr = QCZ*I_ERI_S_Px_Pz_S_vrr+WQZ*I_ERI_S_Px_Pz_S_M1_vrr+oned2e*I_ERI_S_Px_S_S_vrr-rhod2esq*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Py_D2z_S_vrr = QCZ*I_ERI_S_Py_Pz_S_vrr+WQZ*I_ERI_S_Py_Pz_S_M1_vrr+oned2e*I_ERI_S_Py_S_S_vrr-rhod2esq*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_Pz_D2z_S_vrr = QCZ*I_ERI_S_Pz_Pz_S_vrr+WQZ*I_ERI_S_Pz_Pz_S_M1_vrr+oned2e*I_ERI_S_Pz_S_S_vrr-rhod2esq*I_ERI_S_Pz_S_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_P_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_P_S
       * RHS shell quartet name: SQ_ERI_S_P_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       ************************************************************/
      Double I_ERI_S_D2x_Px_S_vrr = PBX*I_ERI_S_Px_Px_S_vrr+WPX*I_ERI_S_Px_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr+oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Dxy_Px_S_vrr = PBY*I_ERI_S_Px_Px_S_vrr+WPY*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_Dxz_Px_S_vrr = PBZ*I_ERI_S_Px_Px_S_vrr+WPZ*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_D2y_Px_S_vrr = PBY*I_ERI_S_Py_Px_S_vrr+WPY*I_ERI_S_Py_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Dyz_Px_S_vrr = PBZ*I_ERI_S_Py_Px_S_vrr+WPZ*I_ERI_S_Py_Px_S_M1_vrr;
      Double I_ERI_S_D2z_Px_S_vrr = PBZ*I_ERI_S_Pz_Px_S_vrr+WPZ*I_ERI_S_Pz_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_D2x_Py_S_vrr = PBX*I_ERI_S_Px_Py_S_vrr+WPX*I_ERI_S_Px_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Dxy_Py_S_vrr = PBY*I_ERI_S_Px_Py_S_vrr+WPY*I_ERI_S_Px_Py_S_M1_vrr+oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Dxz_Py_S_vrr = PBZ*I_ERI_S_Px_Py_S_vrr+WPZ*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_D2y_Py_S_vrr = PBY*I_ERI_S_Py_Py_S_vrr+WPY*I_ERI_S_Py_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr+oned2k*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_Dyz_Py_S_vrr = PBZ*I_ERI_S_Py_Py_S_vrr+WPZ*I_ERI_S_Py_Py_S_M1_vrr;
      Double I_ERI_S_D2z_Py_S_vrr = PBZ*I_ERI_S_Pz_Py_S_vrr+WPZ*I_ERI_S_Pz_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_D2x_Pz_S_vrr = PBX*I_ERI_S_Px_Pz_S_vrr+WPX*I_ERI_S_Px_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Dxy_Pz_S_vrr = PBY*I_ERI_S_Px_Pz_S_vrr+WPY*I_ERI_S_Px_Pz_S_M1_vrr;
      Double I_ERI_S_Dxz_Pz_S_vrr = PBZ*I_ERI_S_Px_Pz_S_vrr+WPZ*I_ERI_S_Px_Pz_S_M1_vrr+oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_D2y_Pz_S_vrr = PBY*I_ERI_S_Py_Pz_S_vrr+WPY*I_ERI_S_Py_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Dyz_Pz_S_vrr = PBZ*I_ERI_S_Py_Pz_S_vrr+WPZ*I_ERI_S_Py_Pz_S_M1_vrr+oned2k*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_D2z_Pz_S_vrr = PBZ*I_ERI_S_Pz_Pz_S_vrr+WPZ*I_ERI_S_Pz_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr+oned2k*I_ERI_S_Pz_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_vrr = QCX*I_ERI_Px_S_Px_S_vrr+WQX*I_ERI_Px_S_Px_S_M1_vrr+oned2e*I_ERI_Px_S_S_S_vrr-rhod2esq*I_ERI_Px_S_S_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_D2x_S_vrr = QCX*I_ERI_Py_S_Px_S_vrr+WQX*I_ERI_Py_S_Px_S_M1_vrr+oned2e*I_ERI_Py_S_S_S_vrr-rhod2esq*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_D2x_S_vrr = QCX*I_ERI_Pz_S_Px_S_vrr+WQX*I_ERI_Pz_S_Px_S_M1_vrr+oned2e*I_ERI_Pz_S_S_S_vrr-rhod2esq*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_Px_S_Dxy_S_vrr = QCY*I_ERI_Px_S_Px_S_vrr+WQY*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_Dxy_S_vrr = QCY*I_ERI_Py_S_Px_S_vrr+WQY*I_ERI_Py_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Dxy_S_vrr = QCY*I_ERI_Pz_S_Px_S_vrr+WQY*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_Dxz_S_vrr = QCZ*I_ERI_Px_S_Px_S_vrr+WQZ*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_Dxz_S_vrr = QCZ*I_ERI_Py_S_Px_S_vrr+WQZ*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Dxz_S_vrr = QCZ*I_ERI_Pz_S_Px_S_vrr+WQZ*I_ERI_Pz_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_D2y_S_vrr = QCY*I_ERI_Px_S_Py_S_vrr+WQY*I_ERI_Px_S_Py_S_M1_vrr+oned2e*I_ERI_Px_S_S_S_vrr-rhod2esq*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Py_S_D2y_S_vrr = QCY*I_ERI_Py_S_Py_S_vrr+WQY*I_ERI_Py_S_Py_S_M1_vrr+oned2e*I_ERI_Py_S_S_S_vrr-rhod2esq*I_ERI_Py_S_S_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Pz_S_D2y_S_vrr = QCY*I_ERI_Pz_S_Py_S_vrr+WQY*I_ERI_Pz_S_Py_S_M1_vrr+oned2e*I_ERI_Pz_S_S_S_vrr-rhod2esq*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_Px_S_Dyz_S_vrr = QCZ*I_ERI_Px_S_Py_S_vrr+WQZ*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_Dyz_S_vrr = QCZ*I_ERI_Py_S_Py_S_vrr+WQZ*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Pz_S_Dyz_S_vrr = QCZ*I_ERI_Pz_S_Py_S_vrr+WQZ*I_ERI_Pz_S_Py_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_D2z_S_vrr = QCZ*I_ERI_Px_S_Pz_S_vrr+WQZ*I_ERI_Px_S_Pz_S_M1_vrr+oned2e*I_ERI_Px_S_S_S_vrr-rhod2esq*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Py_S_D2z_S_vrr = QCZ*I_ERI_Py_S_Pz_S_vrr+WQZ*I_ERI_Py_S_Pz_S_M1_vrr+oned2e*I_ERI_Py_S_S_S_vrr-rhod2esq*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_D2z_S_vrr = QCZ*I_ERI_Pz_S_Pz_S_vrr+WQZ*I_ERI_Pz_S_Pz_S_M1_vrr+oned2e*I_ERI_Pz_S_S_S_vrr-rhod2esq*I_ERI_Pz_S_S_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_vrr = PAX*I_ERI_Px_S_Px_S_vrr+WPX*I_ERI_Px_S_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_Px_S_vrr = PAY*I_ERI_Px_S_Px_S_vrr+WPY*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Px_S_vrr = PAZ*I_ERI_Px_S_Px_S_vrr+WPZ*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Px_S_vrr = PAY*I_ERI_Py_S_Px_S_vrr+WPY*I_ERI_Py_S_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Px_S_vrr = PAZ*I_ERI_Py_S_Px_S_vrr+WPZ*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Px_S_vrr = PAZ*I_ERI_Pz_S_Px_S_vrr+WPZ*I_ERI_Pz_S_Px_S_M1_vrr+oned2z*I_ERI_S_S_Px_S_vrr-rhod2zsq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_Py_S_vrr = PAX*I_ERI_Px_S_Py_S_vrr+WPX*I_ERI_Px_S_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_Py_S_vrr = PAY*I_ERI_Px_S_Py_S_vrr+WPY*I_ERI_Px_S_Py_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_Py_S_vrr = PAZ*I_ERI_Px_S_Py_S_vrr+WPZ*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Py_S_vrr = PAY*I_ERI_Py_S_Py_S_vrr+WPY*I_ERI_Py_S_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_Py_S_vrr = PAZ*I_ERI_Py_S_Py_S_vrr+WPZ*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Py_S_vrr = PAZ*I_ERI_Pz_S_Py_S_vrr+WPZ*I_ERI_Pz_S_Py_S_M1_vrr+oned2z*I_ERI_S_S_Py_S_vrr-rhod2zsq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_Pz_S_vrr = PAX*I_ERI_Px_S_Pz_S_vrr+WPX*I_ERI_Px_S_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Dxy_S_Pz_S_vrr = PAY*I_ERI_Px_S_Pz_S_vrr+WPY*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Dxz_S_Pz_S_vrr = PAZ*I_ERI_Px_S_Pz_S_vrr+WPZ*I_ERI_Px_S_Pz_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Pz_S_vrr = PAY*I_ERI_Py_S_Pz_S_vrr+WPY*I_ERI_Py_S_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Dyz_S_Pz_S_vrr = PAZ*I_ERI_Py_S_Pz_S_vrr+WPZ*I_ERI_Py_S_Pz_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Pz_S_vrr = PAZ*I_ERI_Pz_S_Pz_S_vrr+WPZ*I_ERI_Pz_S_Pz_S_M1_vrr+oned2z*I_ERI_S_S_Pz_S_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_F_S_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_D_S_S
       * RHS shell quartet name: SQ_ERI_S_D_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       ************************************************************/
      Double I_ERI_S_F3x_S_S_vrr = PBX*I_ERI_S_D2x_S_S_vrr+WPX*I_ERI_S_D2x_S_S_M1_vrr+2*oned2z*I_ERI_S_Px_S_S_vrr-2*rhod2zsq*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_F2xy_S_S_vrr = PBY*I_ERI_S_D2x_S_S_vrr+WPY*I_ERI_S_D2x_S_S_M1_vrr;
      Double I_ERI_S_F2xz_S_S_vrr = PBZ*I_ERI_S_D2x_S_S_vrr+WPZ*I_ERI_S_D2x_S_S_M1_vrr;
      Double I_ERI_S_Fx2y_S_S_vrr = PBX*I_ERI_S_D2y_S_S_vrr+WPX*I_ERI_S_D2y_S_S_M1_vrr;
      Double I_ERI_S_Fxyz_S_S_vrr = PBZ*I_ERI_S_Dxy_S_S_vrr+WPZ*I_ERI_S_Dxy_S_S_M1_vrr;
      Double I_ERI_S_Fx2z_S_S_vrr = PBX*I_ERI_S_D2z_S_S_vrr+WPX*I_ERI_S_D2z_S_S_M1_vrr;
      Double I_ERI_S_F3y_S_S_vrr = PBY*I_ERI_S_D2y_S_S_vrr+WPY*I_ERI_S_D2y_S_S_M1_vrr+2*oned2z*I_ERI_S_Py_S_S_vrr-2*rhod2zsq*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_F2yz_S_S_vrr = PBZ*I_ERI_S_D2y_S_S_vrr+WPZ*I_ERI_S_D2y_S_S_M1_vrr;
      Double I_ERI_S_Fy2z_S_S_vrr = PBY*I_ERI_S_D2z_S_S_vrr+WPY*I_ERI_S_D2z_S_S_M1_vrr;
      Double I_ERI_S_F3z_S_S_vrr = PBZ*I_ERI_S_D2z_S_S_vrr+WPZ*I_ERI_S_D2z_S_S_M1_vrr+2*oned2z*I_ERI_S_Pz_S_S_vrr-2*rhod2zsq*I_ERI_S_Pz_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_C0_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C0_aa_coefs = ic2*jc2*alpha*alpha;
      I_ERI_D2x_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C0_aa += SQ_ERI_D_S_S_S_C0_aa_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C0_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C0_a_coefs = ic2*jc2*alpha;
      I_ERI_S_S_S_S_C0_a += SQ_ERI_S_S_S_S_C0_a_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1_aa_coefs = ic2_1*jc2*alpha*alpha;
      I_ERI_F3x_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1_aa += SQ_ERI_F_S_S_S_C1_aa_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1_a_coefs = ic2_1*jc2*alpha;
      I_ERI_Px_S_S_S_C1_a += SQ_ERI_P_S_S_S_C1_a_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1_a += SQ_ERI_P_S_S_S_C1_a_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1_a += SQ_ERI_P_S_S_S_C1_a_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1000_a_coefs = ic2_2*jc2*alpha;
      I_ERI_S_Px_S_S_C1000_a += SQ_ERI_S_P_S_S_C1000_a_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1000_a += SQ_ERI_S_P_S_S_C1000_a_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1000_a += SQ_ERI_S_P_S_S_C1000_a_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1_b_coefs = ic2_1*jc2*beta;
      I_ERI_S_Px_S_S_C1_b += SQ_ERI_S_P_S_S_C1_b_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1_b += SQ_ERI_S_P_S_S_C1_b_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1_b += SQ_ERI_S_P_S_S_C1_b_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1000_a_coefs = ic2_2*jc2*alpha;
      I_ERI_Px_S_S_S_C1000_a += SQ_ERI_P_S_S_S_C1000_a_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1000_a += SQ_ERI_P_S_S_S_C1000_a_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1000_a += SQ_ERI_P_S_S_S_C1000_a_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001_a_coefs = ic2_3*jc2*alpha;
      I_ERI_D2x_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001_a += SQ_ERI_D_S_S_S_C1001_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_S_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_S_S_C1001_b_coefs = ic2_3*jc2*beta;
      I_ERI_S_D2x_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_D2x_S_S_vrr;
      I_ERI_S_Dxy_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_Dxy_S_S_vrr;
      I_ERI_S_Dxz_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_Dxz_S_S_vrr;
      I_ERI_S_D2y_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_D2y_S_S_vrr;
      I_ERI_S_Dyz_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_Dyz_S_S_vrr;
      I_ERI_S_D2z_S_S_C1001_b += SQ_ERI_S_D_S_S_C1001_b_coefs*I_ERI_S_D2z_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C1001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C1001_coefs = ic2_3*jc2;
      I_ERI_S_S_S_S_C1001 += SQ_ERI_S_S_S_S_C1001_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C0_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C0_ac_coefs = ic2*jc2*alpha*gamma;
      I_ERI_Px_S_Px_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C0_ac += SQ_ERI_P_S_P_S_C0_ac_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1_ac_coefs = ic2_1*jc2*alpha*gamma;
      I_ERI_D2x_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1_ac += SQ_ERI_D_S_P_S_C1_ac_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1_c_coefs = ic2_1*jc2*gamma;
      I_ERI_S_S_Px_S_C1_c += SQ_ERI_S_S_P_S_C1_c_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1_c += SQ_ERI_S_S_P_S_C1_c_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1_c += SQ_ERI_S_S_P_S_C1_c_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001_c_coefs = ic2_3*jc2*gamma;
      I_ERI_S_Px_Px_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001_c += SQ_ERI_S_P_P_S_C1001_c_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_S_S_C0_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_S_S_C0_bb_coefs = ic2*jc2*beta*beta;
      I_ERI_S_D2x_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_D2x_S_S_vrr;
      I_ERI_S_Dxy_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_Dxy_S_S_vrr;
      I_ERI_S_Dxz_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_Dxz_S_S_vrr;
      I_ERI_S_D2y_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_D2y_S_S_vrr;
      I_ERI_S_Dyz_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_Dyz_S_S_vrr;
      I_ERI_S_D2z_S_S_C0_bb += SQ_ERI_S_D_S_S_C0_bb_coefs*I_ERI_S_D2z_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C0_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C0_b_coefs = ic2*jc2*beta;
      I_ERI_S_S_S_S_C0_b += SQ_ERI_S_S_S_S_C0_b_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1_b_coefs = ic2_1*jc2*beta;
      I_ERI_Px_S_S_S_C1_b += SQ_ERI_P_S_S_S_C1_b_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1_b += SQ_ERI_P_S_S_S_C1_b_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1_b += SQ_ERI_P_S_S_S_C1_b_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_F_S_S_C1000_bb_coefs = ic2_2*jc2*beta*beta;
      I_ERI_S_F3x_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F3x_S_S_vrr;
      I_ERI_S_F2xy_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F2xy_S_S_vrr;
      I_ERI_S_F2xz_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F2xz_S_S_vrr;
      I_ERI_S_Fx2y_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_Fx2y_S_S_vrr;
      I_ERI_S_Fxyz_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_Fxyz_S_S_vrr;
      I_ERI_S_Fx2z_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_Fx2z_S_S_vrr;
      I_ERI_S_F3y_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F3y_S_S_vrr;
      I_ERI_S_F2yz_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F2yz_S_S_vrr;
      I_ERI_S_Fy2z_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_Fy2z_S_S_vrr;
      I_ERI_S_F3z_S_S_C1000_bb += SQ_ERI_S_F_S_S_C1000_bb_coefs*I_ERI_S_F3z_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1000_b_coefs = ic2_2*jc2*beta;
      I_ERI_S_Px_S_S_C1000_b += SQ_ERI_S_P_S_S_C1000_b_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1000_b += SQ_ERI_S_P_S_S_C1000_b_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1000_b += SQ_ERI_S_P_S_S_C1000_b_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C0_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C0_bc_coefs = ic2*jc2*beta*gamma;
      I_ERI_S_Px_Px_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C0_bc += SQ_ERI_S_P_P_S_C0_bc_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_P_S_C1000_bc_coefs = ic2_2*jc2*beta*gamma;
      I_ERI_S_D2x_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2x_Px_S_vrr;
      I_ERI_S_Dxy_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxy_Px_S_vrr;
      I_ERI_S_Dxz_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxz_Px_S_vrr;
      I_ERI_S_D2y_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2y_Px_S_vrr;
      I_ERI_S_Dyz_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dyz_Px_S_vrr;
      I_ERI_S_D2z_Px_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2z_Px_S_vrr;
      I_ERI_S_D2x_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2x_Py_S_vrr;
      I_ERI_S_Dxy_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxy_Py_S_vrr;
      I_ERI_S_Dxz_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxz_Py_S_vrr;
      I_ERI_S_D2y_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2y_Py_S_vrr;
      I_ERI_S_Dyz_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dyz_Py_S_vrr;
      I_ERI_S_D2z_Py_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2z_Py_S_vrr;
      I_ERI_S_D2x_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2x_Pz_S_vrr;
      I_ERI_S_Dxy_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxy_Pz_S_vrr;
      I_ERI_S_Dxz_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dxz_Pz_S_vrr;
      I_ERI_S_D2y_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2y_Pz_S_vrr;
      I_ERI_S_Dyz_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_Dyz_Pz_S_vrr;
      I_ERI_S_D2z_Pz_S_C1000_bc += SQ_ERI_S_D_P_S_C1000_bc_coefs*I_ERI_S_D2z_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1000_c_coefs = ic2_2*jc2*gamma;
      I_ERI_S_S_Px_S_C1000_c += SQ_ERI_S_S_P_S_C1000_c_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000_c += SQ_ERI_S_S_P_S_C1000_c_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000_c += SQ_ERI_S_S_P_S_C1000_c_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001_c_coefs = ic2_3*jc2*gamma;
      I_ERI_Px_S_Px_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001_c += SQ_ERI_P_S_P_S_C1001_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C0_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C0_cc_coefs = ic2*jc2*gamma*gamma;
      I_ERI_S_S_D2x_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C0_cc += SQ_ERI_S_S_D_S_C0_cc_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C0_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C0_c_coefs = ic2*jc2*gamma;
      I_ERI_S_S_S_S_C0_c += SQ_ERI_S_S_S_S_C0_c_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1_cc_coefs = ic2_1*jc2*gamma*gamma;
      I_ERI_Px_S_D2x_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1_cc += SQ_ERI_P_S_D_S_C1_cc_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1_c_coefs = ic2_1*jc2*gamma;
      I_ERI_Px_S_S_S_C1_c += SQ_ERI_P_S_S_S_C1_c_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1_c += SQ_ERI_P_S_S_S_C1_c_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1_c += SQ_ERI_P_S_S_S_C1_c_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1000_cc_coefs = ic2_2*jc2*gamma*gamma;
      I_ERI_S_Px_D2x_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1000_cc += SQ_ERI_S_P_D_S_C1000_cc_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1000_c_coefs = ic2_2*jc2*gamma;
      I_ERI_S_Px_S_S_C1000_c += SQ_ERI_S_P_S_S_C1000_c_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1000_c += SQ_ERI_S_P_S_S_C1000_c_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1000_c += SQ_ERI_S_P_S_S_C1000_c_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1000_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1000_aa_coefs = ic2_2*jc2*alpha*alpha;
      I_ERI_F3x_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1000_aa += SQ_ERI_F_S_S_S_C1000_aa_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1000_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1000_aa_coefs = ic2_2*jc2*alpha*alpha;
      I_ERI_D2x_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1000_aa += SQ_ERI_D_S_S_S_C1000_aa_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1001_aa_coefs = ic2_3*jc2*alpha*alpha;
      I_ERI_G4x_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1001_aa += SQ_ERI_G_S_S_S_C1001_aa_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001_aa_coefs = ic2_3*jc2*alpha*alpha;
      I_ERI_F3x_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001_aa += SQ_ERI_F_S_S_S_C1001_aa_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001_a_coefs = ic2_3*jc2*alpha;
      I_ERI_Px_S_S_S_C1001_a += SQ_ERI_P_S_S_S_C1001_a_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001_a += SQ_ERI_P_S_S_S_C1001_a_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001_a += SQ_ERI_P_S_S_S_C1001_a_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C0_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C0_ab_coefs = ic2*jc2*alpha*beta;
      I_ERI_D2x_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C0_ab += SQ_ERI_D_S_S_S_C0_ab_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C0_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C0_ab_coefs = ic2*jc2*alpha*beta;
      I_ERI_Px_S_S_S_C0_ab += SQ_ERI_P_S_S_S_C0_ab_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C0_ab += SQ_ERI_P_S_S_S_C0_ab_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C0_ab += SQ_ERI_P_S_S_S_C0_ab_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1_ab_coefs = ic2_1*jc2*alpha*beta;
      I_ERI_F3x_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1_ab += SQ_ERI_F_S_S_S_C1_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1_ab_coefs = ic2_1*jc2*alpha*beta;
      I_ERI_D2x_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1_ab += SQ_ERI_D_S_S_S_C1_ab_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000_ac_coefs = ic2_2*jc2*alpha*gamma;
      I_ERI_D2x_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000_ac += SQ_ERI_D_S_P_S_C1000_ac_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000_ac_coefs = ic2_2*jc2*alpha*gamma;
      I_ERI_Px_S_Px_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000_ac += SQ_ERI_P_S_P_S_C1000_ac_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001_ac_coefs = ic2_3*jc2*alpha*gamma;
      I_ERI_F3x_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001_ac += SQ_ERI_F_S_P_S_C1001_ac_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001_ac_coefs = ic2_3*jc2*alpha*gamma;
      I_ERI_D2x_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001_ac += SQ_ERI_D_S_P_S_C1001_ac_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001_b_coefs = ic2_3*jc2*beta;
      I_ERI_D2x_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001_b += SQ_ERI_D_S_S_S_C1001_b_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001_b_coefs = ic2_3*jc2*beta;
      I_ERI_Px_S_S_S_C1001_b += SQ_ERI_P_S_S_S_C1001_b_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001_b += SQ_ERI_P_S_S_S_C1001_b_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001_b += SQ_ERI_P_S_S_S_C1001_b_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1_bc_coefs = ic2_1*jc2*beta*gamma;
      I_ERI_D2x_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1_bc += SQ_ERI_D_S_P_S_C1_bc_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1_bc_coefs = ic2_1*jc2*beta*gamma;
      I_ERI_Px_S_Px_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1_bc += SQ_ERI_P_S_P_S_C1_bc_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001_cc_coefs = ic2_3*jc2*gamma*gamma;
      I_ERI_D2x_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001_cc += SQ_ERI_D_S_D_S_C1001_cc_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001_cc_coefs = ic2_3*jc2*gamma*gamma;
      I_ERI_Px_S_D2x_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001_cc += SQ_ERI_P_S_D_S_C1001_cc_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001_c_coefs = ic2_3*jc2*gamma;
      I_ERI_D2x_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001_c += SQ_ERI_D_S_S_S_C1001_c_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001_c_coefs = ic2_3*jc2*gamma;
      I_ERI_Px_S_S_S_C1001_c += SQ_ERI_P_S_S_S_C1001_c_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001_c += SQ_ERI_P_S_S_S_C1001_c_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001_c += SQ_ERI_P_S_S_S_C1001_c_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1000_ab_coefs = ic2_2*jc2*alpha*beta;
      I_ERI_F3x_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1000_ab += SQ_ERI_F_S_S_S_C1000_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1000_ab_coefs = ic2_2*jc2*alpha*beta;
      I_ERI_D2x_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1000_ab += SQ_ERI_D_S_S_S_C1000_ab_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1000_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1000_ab_coefs = ic2_2*jc2*alpha*beta;
      I_ERI_Px_S_S_S_C1000_ab += SQ_ERI_P_S_S_S_C1000_ab_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1000_ab += SQ_ERI_P_S_S_S_C1000_ab_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1000_ab += SQ_ERI_P_S_S_S_C1000_ab_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1001_ab_coefs = ic2_3*jc2*alpha*beta;
      I_ERI_G4x_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1001_ab += SQ_ERI_G_S_S_S_C1001_ab_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001_ab_coefs = ic2_3*jc2*alpha*beta;
      I_ERI_F3x_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001_ab += SQ_ERI_F_S_S_S_C1001_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001_ab_coefs = ic2_3*jc2*alpha*beta;
      I_ERI_D2x_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001_ab += SQ_ERI_D_S_S_S_C1001_ab_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_F3x_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1_bb += SQ_ERI_F_S_S_S_C1_bb_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_D2x_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1_bb += SQ_ERI_D_S_S_S_C1_bb_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1_bb_coefs = ic2_1*jc2*beta*beta;
      I_ERI_Px_S_S_S_C1_bb += SQ_ERI_P_S_S_S_C1_bb_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1_bb += SQ_ERI_P_S_S_S_C1_bb_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1_bb += SQ_ERI_P_S_S_S_C1_bb_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001_bc_coefs = ic2_3*jc2*beta*gamma;
      I_ERI_F3x_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001_bc += SQ_ERI_F_S_P_S_C1001_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001_bc_coefs = ic2_3*jc2*beta*gamma;
      I_ERI_D2x_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001_bc += SQ_ERI_D_S_P_S_C1001_bc_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001_bc_coefs = ic2_3*jc2*beta*gamma;
      I_ERI_Px_S_Px_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001_bc += SQ_ERI_P_S_P_S_C1001_bc_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1001_bb_coefs = ic2_3*jc2*beta*beta;
      I_ERI_G4x_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1001_bb += SQ_ERI_G_S_S_S_C1001_bb_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001_bb_coefs = ic2_3*jc2*beta*beta;
      I_ERI_F3x_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001_bb += SQ_ERI_F_S_S_S_C1001_bb_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001_bb_coefs = ic2_3*jc2*beta*beta;
      I_ERI_D2x_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001_bb += SQ_ERI_D_S_S_S_C1001_bb_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001_bb_coefs = ic2_3*jc2*beta*beta;
      I_ERI_Px_S_S_S_C1001_bb += SQ_ERI_P_S_S_S_C1001_bb_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001_bb += SQ_ERI_P_S_S_S_C1001_bb_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001_bb += SQ_ERI_P_S_S_S_C1001_bb_coefs*I_ERI_Pz_S_S_S_vrr;
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
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001_a
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1001_a = I_ERI_D2x_S_S_S_C1001_a+ABX*I_ERI_Px_S_S_S_C1001_a;
  Double I_ERI_Py_Px_S_S_C1001_a = I_ERI_Dxy_S_S_S_C1001_a+ABX*I_ERI_Py_S_S_S_C1001_a;
  Double I_ERI_Pz_Px_S_S_C1001_a = I_ERI_Dxz_S_S_S_C1001_a+ABX*I_ERI_Pz_S_S_S_C1001_a;
  Double I_ERI_Px_Py_S_S_C1001_a = I_ERI_Dxy_S_S_S_C1001_a+ABY*I_ERI_Px_S_S_S_C1001_a;
  Double I_ERI_Py_Py_S_S_C1001_a = I_ERI_D2y_S_S_S_C1001_a+ABY*I_ERI_Py_S_S_S_C1001_a;
  Double I_ERI_Pz_Py_S_S_C1001_a = I_ERI_Dyz_S_S_S_C1001_a+ABY*I_ERI_Pz_S_S_S_C1001_a;
  Double I_ERI_Px_Pz_S_S_C1001_a = I_ERI_Dxz_S_S_S_C1001_a+ABZ*I_ERI_Px_S_S_S_C1001_a;
  Double I_ERI_Py_Pz_S_S_C1001_a = I_ERI_Dyz_S_S_S_C1001_a+ABZ*I_ERI_Py_S_S_S_C1001_a;
  Double I_ERI_Pz_Pz_S_S_C1001_a = I_ERI_D2z_S_S_S_C1001_a+ABZ*I_ERI_Pz_S_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001_b
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1001_b = I_ERI_D2x_S_S_S_C1001_b+ABX*I_ERI_Px_S_S_S_C1001_b;
  Double I_ERI_Py_Px_S_S_C1001_b = I_ERI_Dxy_S_S_S_C1001_b+ABX*I_ERI_Py_S_S_S_C1001_b;
  Double I_ERI_Pz_Px_S_S_C1001_b = I_ERI_Dxz_S_S_S_C1001_b+ABX*I_ERI_Pz_S_S_S_C1001_b;
  Double I_ERI_Px_Py_S_S_C1001_b = I_ERI_Dxy_S_S_S_C1001_b+ABY*I_ERI_Px_S_S_S_C1001_b;
  Double I_ERI_Py_Py_S_S_C1001_b = I_ERI_D2y_S_S_S_C1001_b+ABY*I_ERI_Py_S_S_S_C1001_b;
  Double I_ERI_Pz_Py_S_S_C1001_b = I_ERI_Dyz_S_S_S_C1001_b+ABY*I_ERI_Pz_S_S_S_C1001_b;
  Double I_ERI_Px_Pz_S_S_C1001_b = I_ERI_Dxz_S_S_S_C1001_b+ABZ*I_ERI_Px_S_S_S_C1001_b;
  Double I_ERI_Py_Pz_S_S_C1001_b = I_ERI_Dyz_S_S_S_C1001_b+ABZ*I_ERI_Py_S_S_S_C1001_b;
  Double I_ERI_Pz_Pz_S_S_C1001_b = I_ERI_D2z_S_S_S_C1001_b+ABZ*I_ERI_Pz_S_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001_c
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1001_c = I_ERI_D2x_S_S_S_C1001_c+ABX*I_ERI_Px_S_S_S_C1001_c;
  Double I_ERI_Py_Px_S_S_C1001_c = I_ERI_Dxy_S_S_S_C1001_c+ABX*I_ERI_Py_S_S_S_C1001_c;
  Double I_ERI_Pz_Px_S_S_C1001_c = I_ERI_Dxz_S_S_S_C1001_c+ABX*I_ERI_Pz_S_S_S_C1001_c;
  Double I_ERI_Px_Py_S_S_C1001_c = I_ERI_Dxy_S_S_S_C1001_c+ABY*I_ERI_Px_S_S_S_C1001_c;
  Double I_ERI_Py_Py_S_S_C1001_c = I_ERI_D2y_S_S_S_C1001_c+ABY*I_ERI_Py_S_S_S_C1001_c;
  Double I_ERI_Pz_Py_S_S_C1001_c = I_ERI_Dyz_S_S_S_C1001_c+ABY*I_ERI_Pz_S_S_S_C1001_c;
  Double I_ERI_Px_Pz_S_S_C1001_c = I_ERI_Dxz_S_S_S_C1001_c+ABZ*I_ERI_Px_S_S_S_C1001_c;
  Double I_ERI_Py_Pz_S_S_C1001_c = I_ERI_Dyz_S_S_S_C1001_c+ABZ*I_ERI_Py_S_S_S_C1001_c;
  Double I_ERI_Pz_Pz_S_S_C1001_c = I_ERI_D2z_S_S_S_C1001_c+ABZ*I_ERI_Pz_S_S_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000_aa
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1000_aa = I_ERI_F3x_S_S_S_C1000_aa+ABX*I_ERI_D2x_S_S_S_C1000_aa;
  Double I_ERI_Dxy_Px_S_S_C1000_aa = I_ERI_F2xy_S_S_S_C1000_aa+ABX*I_ERI_Dxy_S_S_S_C1000_aa;
  Double I_ERI_Dxz_Px_S_S_C1000_aa = I_ERI_F2xz_S_S_S_C1000_aa+ABX*I_ERI_Dxz_S_S_S_C1000_aa;
  Double I_ERI_D2y_Px_S_S_C1000_aa = I_ERI_Fx2y_S_S_S_C1000_aa+ABX*I_ERI_D2y_S_S_S_C1000_aa;
  Double I_ERI_Dyz_Px_S_S_C1000_aa = I_ERI_Fxyz_S_S_S_C1000_aa+ABX*I_ERI_Dyz_S_S_S_C1000_aa;
  Double I_ERI_D2z_Px_S_S_C1000_aa = I_ERI_Fx2z_S_S_S_C1000_aa+ABX*I_ERI_D2z_S_S_S_C1000_aa;
  Double I_ERI_D2x_Py_S_S_C1000_aa = I_ERI_F2xy_S_S_S_C1000_aa+ABY*I_ERI_D2x_S_S_S_C1000_aa;
  Double I_ERI_Dxy_Py_S_S_C1000_aa = I_ERI_Fx2y_S_S_S_C1000_aa+ABY*I_ERI_Dxy_S_S_S_C1000_aa;
  Double I_ERI_Dxz_Py_S_S_C1000_aa = I_ERI_Fxyz_S_S_S_C1000_aa+ABY*I_ERI_Dxz_S_S_S_C1000_aa;
  Double I_ERI_D2y_Py_S_S_C1000_aa = I_ERI_F3y_S_S_S_C1000_aa+ABY*I_ERI_D2y_S_S_S_C1000_aa;
  Double I_ERI_Dyz_Py_S_S_C1000_aa = I_ERI_F2yz_S_S_S_C1000_aa+ABY*I_ERI_Dyz_S_S_S_C1000_aa;
  Double I_ERI_D2z_Py_S_S_C1000_aa = I_ERI_Fy2z_S_S_S_C1000_aa+ABY*I_ERI_D2z_S_S_S_C1000_aa;
  Double I_ERI_D2x_Pz_S_S_C1000_aa = I_ERI_F2xz_S_S_S_C1000_aa+ABZ*I_ERI_D2x_S_S_S_C1000_aa;
  Double I_ERI_Dxy_Pz_S_S_C1000_aa = I_ERI_Fxyz_S_S_S_C1000_aa+ABZ*I_ERI_Dxy_S_S_S_C1000_aa;
  Double I_ERI_Dxz_Pz_S_S_C1000_aa = I_ERI_Fx2z_S_S_S_C1000_aa+ABZ*I_ERI_Dxz_S_S_S_C1000_aa;
  Double I_ERI_D2y_Pz_S_S_C1000_aa = I_ERI_F2yz_S_S_S_C1000_aa+ABZ*I_ERI_D2y_S_S_S_C1000_aa;
  Double I_ERI_Dyz_Pz_S_S_C1000_aa = I_ERI_Fy2z_S_S_S_C1000_aa+ABZ*I_ERI_Dyz_S_S_S_C1000_aa;
  Double I_ERI_D2z_Pz_S_S_C1000_aa = I_ERI_F3z_S_S_S_C1000_aa+ABZ*I_ERI_D2z_S_S_S_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_aa
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1001_aa = I_ERI_G4x_S_S_S_C1001_aa+ABX*I_ERI_F3x_S_S_S_C1001_aa;
  Double I_ERI_F2xy_Px_S_S_C1001_aa = I_ERI_G3xy_S_S_S_C1001_aa+ABX*I_ERI_F2xy_S_S_S_C1001_aa;
  Double I_ERI_F2xz_Px_S_S_C1001_aa = I_ERI_G3xz_S_S_S_C1001_aa+ABX*I_ERI_F2xz_S_S_S_C1001_aa;
  Double I_ERI_Fx2y_Px_S_S_C1001_aa = I_ERI_G2x2y_S_S_S_C1001_aa+ABX*I_ERI_Fx2y_S_S_S_C1001_aa;
  Double I_ERI_Fxyz_Px_S_S_C1001_aa = I_ERI_G2xyz_S_S_S_C1001_aa+ABX*I_ERI_Fxyz_S_S_S_C1001_aa;
  Double I_ERI_Fx2z_Px_S_S_C1001_aa = I_ERI_G2x2z_S_S_S_C1001_aa+ABX*I_ERI_Fx2z_S_S_S_C1001_aa;
  Double I_ERI_F3y_Px_S_S_C1001_aa = I_ERI_Gx3y_S_S_S_C1001_aa+ABX*I_ERI_F3y_S_S_S_C1001_aa;
  Double I_ERI_F2yz_Px_S_S_C1001_aa = I_ERI_Gx2yz_S_S_S_C1001_aa+ABX*I_ERI_F2yz_S_S_S_C1001_aa;
  Double I_ERI_Fy2z_Px_S_S_C1001_aa = I_ERI_Gxy2z_S_S_S_C1001_aa+ABX*I_ERI_Fy2z_S_S_S_C1001_aa;
  Double I_ERI_F3z_Px_S_S_C1001_aa = I_ERI_Gx3z_S_S_S_C1001_aa+ABX*I_ERI_F3z_S_S_S_C1001_aa;
  Double I_ERI_F3x_Py_S_S_C1001_aa = I_ERI_G3xy_S_S_S_C1001_aa+ABY*I_ERI_F3x_S_S_S_C1001_aa;
  Double I_ERI_F2xy_Py_S_S_C1001_aa = I_ERI_G2x2y_S_S_S_C1001_aa+ABY*I_ERI_F2xy_S_S_S_C1001_aa;
  Double I_ERI_F2xz_Py_S_S_C1001_aa = I_ERI_G2xyz_S_S_S_C1001_aa+ABY*I_ERI_F2xz_S_S_S_C1001_aa;
  Double I_ERI_Fx2y_Py_S_S_C1001_aa = I_ERI_Gx3y_S_S_S_C1001_aa+ABY*I_ERI_Fx2y_S_S_S_C1001_aa;
  Double I_ERI_Fxyz_Py_S_S_C1001_aa = I_ERI_Gx2yz_S_S_S_C1001_aa+ABY*I_ERI_Fxyz_S_S_S_C1001_aa;
  Double I_ERI_Fx2z_Py_S_S_C1001_aa = I_ERI_Gxy2z_S_S_S_C1001_aa+ABY*I_ERI_Fx2z_S_S_S_C1001_aa;
  Double I_ERI_F3y_Py_S_S_C1001_aa = I_ERI_G4y_S_S_S_C1001_aa+ABY*I_ERI_F3y_S_S_S_C1001_aa;
  Double I_ERI_F2yz_Py_S_S_C1001_aa = I_ERI_G3yz_S_S_S_C1001_aa+ABY*I_ERI_F2yz_S_S_S_C1001_aa;
  Double I_ERI_Fy2z_Py_S_S_C1001_aa = I_ERI_G2y2z_S_S_S_C1001_aa+ABY*I_ERI_Fy2z_S_S_S_C1001_aa;
  Double I_ERI_F3z_Py_S_S_C1001_aa = I_ERI_Gy3z_S_S_S_C1001_aa+ABY*I_ERI_F3z_S_S_S_C1001_aa;
  Double I_ERI_F3x_Pz_S_S_C1001_aa = I_ERI_G3xz_S_S_S_C1001_aa+ABZ*I_ERI_F3x_S_S_S_C1001_aa;
  Double I_ERI_F2xy_Pz_S_S_C1001_aa = I_ERI_G2xyz_S_S_S_C1001_aa+ABZ*I_ERI_F2xy_S_S_S_C1001_aa;
  Double I_ERI_F2xz_Pz_S_S_C1001_aa = I_ERI_G2x2z_S_S_S_C1001_aa+ABZ*I_ERI_F2xz_S_S_S_C1001_aa;
  Double I_ERI_Fx2y_Pz_S_S_C1001_aa = I_ERI_Gx2yz_S_S_S_C1001_aa+ABZ*I_ERI_Fx2y_S_S_S_C1001_aa;
  Double I_ERI_Fxyz_Pz_S_S_C1001_aa = I_ERI_Gxy2z_S_S_S_C1001_aa+ABZ*I_ERI_Fxyz_S_S_S_C1001_aa;
  Double I_ERI_Fx2z_Pz_S_S_C1001_aa = I_ERI_Gx3z_S_S_S_C1001_aa+ABZ*I_ERI_Fx2z_S_S_S_C1001_aa;
  Double I_ERI_F3y_Pz_S_S_C1001_aa = I_ERI_G3yz_S_S_S_C1001_aa+ABZ*I_ERI_F3y_S_S_S_C1001_aa;
  Double I_ERI_F2yz_Pz_S_S_C1001_aa = I_ERI_G2y2z_S_S_S_C1001_aa+ABZ*I_ERI_F2yz_S_S_S_C1001_aa;
  Double I_ERI_Fy2z_Pz_S_S_C1001_aa = I_ERI_Gy3z_S_S_S_C1001_aa+ABZ*I_ERI_Fy2z_S_S_S_C1001_aa;
  Double I_ERI_F3z_Pz_S_S_C1001_aa = I_ERI_G4z_S_S_S_C1001_aa+ABZ*I_ERI_F3z_S_S_S_C1001_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C0_ab
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C0_ab = I_ERI_D2x_S_S_S_C0_ab+ABX*I_ERI_Px_S_S_S_C0_ab;
  Double I_ERI_Py_Px_S_S_C0_ab = I_ERI_Dxy_S_S_S_C0_ab+ABX*I_ERI_Py_S_S_S_C0_ab;
  Double I_ERI_Pz_Px_S_S_C0_ab = I_ERI_Dxz_S_S_S_C0_ab+ABX*I_ERI_Pz_S_S_S_C0_ab;
  Double I_ERI_Px_Py_S_S_C0_ab = I_ERI_Dxy_S_S_S_C0_ab+ABY*I_ERI_Px_S_S_S_C0_ab;
  Double I_ERI_Py_Py_S_S_C0_ab = I_ERI_D2y_S_S_S_C0_ab+ABY*I_ERI_Py_S_S_S_C0_ab;
  Double I_ERI_Pz_Py_S_S_C0_ab = I_ERI_Dyz_S_S_S_C0_ab+ABY*I_ERI_Pz_S_S_S_C0_ab;
  Double I_ERI_Px_Pz_S_S_C0_ab = I_ERI_Dxz_S_S_S_C0_ab+ABZ*I_ERI_Px_S_S_S_C0_ab;
  Double I_ERI_Py_Pz_S_S_C0_ab = I_ERI_Dyz_S_S_S_C0_ab+ABZ*I_ERI_Py_S_S_S_C0_ab;
  Double I_ERI_Pz_Pz_S_S_C0_ab = I_ERI_D2z_S_S_S_C0_ab+ABZ*I_ERI_Pz_S_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_ab
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1_ab = I_ERI_F3x_S_S_S_C1_ab+ABX*I_ERI_D2x_S_S_S_C1_ab;
  Double I_ERI_Dxy_Px_S_S_C1_ab = I_ERI_F2xy_S_S_S_C1_ab+ABX*I_ERI_Dxy_S_S_S_C1_ab;
  Double I_ERI_Dxz_Px_S_S_C1_ab = I_ERI_F2xz_S_S_S_C1_ab+ABX*I_ERI_Dxz_S_S_S_C1_ab;
  Double I_ERI_D2y_Px_S_S_C1_ab = I_ERI_Fx2y_S_S_S_C1_ab+ABX*I_ERI_D2y_S_S_S_C1_ab;
  Double I_ERI_Dyz_Px_S_S_C1_ab = I_ERI_Fxyz_S_S_S_C1_ab+ABX*I_ERI_Dyz_S_S_S_C1_ab;
  Double I_ERI_D2z_Px_S_S_C1_ab = I_ERI_Fx2z_S_S_S_C1_ab+ABX*I_ERI_D2z_S_S_S_C1_ab;
  Double I_ERI_D2x_Py_S_S_C1_ab = I_ERI_F2xy_S_S_S_C1_ab+ABY*I_ERI_D2x_S_S_S_C1_ab;
  Double I_ERI_Dxy_Py_S_S_C1_ab = I_ERI_Fx2y_S_S_S_C1_ab+ABY*I_ERI_Dxy_S_S_S_C1_ab;
  Double I_ERI_Dxz_Py_S_S_C1_ab = I_ERI_Fxyz_S_S_S_C1_ab+ABY*I_ERI_Dxz_S_S_S_C1_ab;
  Double I_ERI_D2y_Py_S_S_C1_ab = I_ERI_F3y_S_S_S_C1_ab+ABY*I_ERI_D2y_S_S_S_C1_ab;
  Double I_ERI_Dyz_Py_S_S_C1_ab = I_ERI_F2yz_S_S_S_C1_ab+ABY*I_ERI_Dyz_S_S_S_C1_ab;
  Double I_ERI_D2z_Py_S_S_C1_ab = I_ERI_Fy2z_S_S_S_C1_ab+ABY*I_ERI_D2z_S_S_S_C1_ab;
  Double I_ERI_D2x_Pz_S_S_C1_ab = I_ERI_F2xz_S_S_S_C1_ab+ABZ*I_ERI_D2x_S_S_S_C1_ab;
  Double I_ERI_Dxy_Pz_S_S_C1_ab = I_ERI_Fxyz_S_S_S_C1_ab+ABZ*I_ERI_Dxy_S_S_S_C1_ab;
  Double I_ERI_Dxz_Pz_S_S_C1_ab = I_ERI_Fx2z_S_S_S_C1_ab+ABZ*I_ERI_Dxz_S_S_S_C1_ab;
  Double I_ERI_D2y_Pz_S_S_C1_ab = I_ERI_F2yz_S_S_S_C1_ab+ABZ*I_ERI_D2y_S_S_S_C1_ab;
  Double I_ERI_Dyz_Pz_S_S_C1_ab = I_ERI_Fy2z_S_S_S_C1_ab+ABZ*I_ERI_Dyz_S_S_S_C1_ab;
  Double I_ERI_D2z_Pz_S_S_C1_ab = I_ERI_F3z_S_S_S_C1_ab+ABZ*I_ERI_D2z_S_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_ab
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1000_ab = I_ERI_D2x_S_S_S_C1000_ab+ABX*I_ERI_Px_S_S_S_C1000_ab;
  Double I_ERI_Py_Px_S_S_C1000_ab = I_ERI_Dxy_S_S_S_C1000_ab+ABX*I_ERI_Py_S_S_S_C1000_ab;
  Double I_ERI_Pz_Px_S_S_C1000_ab = I_ERI_Dxz_S_S_S_C1000_ab+ABX*I_ERI_Pz_S_S_S_C1000_ab;
  Double I_ERI_Px_Py_S_S_C1000_ab = I_ERI_Dxy_S_S_S_C1000_ab+ABY*I_ERI_Px_S_S_S_C1000_ab;
  Double I_ERI_Py_Py_S_S_C1000_ab = I_ERI_D2y_S_S_S_C1000_ab+ABY*I_ERI_Py_S_S_S_C1000_ab;
  Double I_ERI_Pz_Py_S_S_C1000_ab = I_ERI_Dyz_S_S_S_C1000_ab+ABY*I_ERI_Pz_S_S_S_C1000_ab;
  Double I_ERI_Px_Pz_S_S_C1000_ab = I_ERI_Dxz_S_S_S_C1000_ab+ABZ*I_ERI_Px_S_S_S_C1000_ab;
  Double I_ERI_Py_Pz_S_S_C1000_ab = I_ERI_Dyz_S_S_S_C1000_ab+ABZ*I_ERI_Py_S_S_S_C1000_ab;
  Double I_ERI_Pz_Pz_S_S_C1000_ab = I_ERI_D2z_S_S_S_C1000_ab+ABZ*I_ERI_Pz_S_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 4 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000_ab
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1000_ab = I_ERI_F3x_S_S_S_C1000_ab+ABX*I_ERI_D2x_S_S_S_C1000_ab;
  Double I_ERI_Dxy_Px_S_S_C1000_ab = I_ERI_F2xy_S_S_S_C1000_ab+ABX*I_ERI_Dxy_S_S_S_C1000_ab;
  Double I_ERI_Dxz_Px_S_S_C1000_ab = I_ERI_F2xz_S_S_S_C1000_ab+ABX*I_ERI_Dxz_S_S_S_C1000_ab;
  Double I_ERI_D2y_Px_S_S_C1000_ab = I_ERI_Fx2y_S_S_S_C1000_ab+ABX*I_ERI_D2y_S_S_S_C1000_ab;
  Double I_ERI_Dyz_Px_S_S_C1000_ab = I_ERI_Fxyz_S_S_S_C1000_ab+ABX*I_ERI_Dyz_S_S_S_C1000_ab;
  Double I_ERI_D2z_Px_S_S_C1000_ab = I_ERI_Fx2z_S_S_S_C1000_ab+ABX*I_ERI_D2z_S_S_S_C1000_ab;
  Double I_ERI_Dxy_Py_S_S_C1000_ab = I_ERI_Fx2y_S_S_S_C1000_ab+ABY*I_ERI_Dxy_S_S_S_C1000_ab;
  Double I_ERI_Dxz_Py_S_S_C1000_ab = I_ERI_Fxyz_S_S_S_C1000_ab+ABY*I_ERI_Dxz_S_S_S_C1000_ab;
  Double I_ERI_D2y_Py_S_S_C1000_ab = I_ERI_F3y_S_S_S_C1000_ab+ABY*I_ERI_D2y_S_S_S_C1000_ab;
  Double I_ERI_Dyz_Py_S_S_C1000_ab = I_ERI_F2yz_S_S_S_C1000_ab+ABY*I_ERI_Dyz_S_S_S_C1000_ab;
  Double I_ERI_D2z_Py_S_S_C1000_ab = I_ERI_Fy2z_S_S_S_C1000_ab+ABY*I_ERI_D2z_S_S_S_C1000_ab;
  Double I_ERI_Dxz_Pz_S_S_C1000_ab = I_ERI_Fx2z_S_S_S_C1000_ab+ABZ*I_ERI_Dxz_S_S_S_C1000_ab;
  Double I_ERI_Dyz_Pz_S_S_C1000_ab = I_ERI_Fy2z_S_S_S_C1000_ab+ABZ*I_ERI_Dyz_S_S_S_C1000_ab;
  Double I_ERI_D2z_Pz_S_S_C1000_ab = I_ERI_F3z_S_S_S_C1000_ab+ABZ*I_ERI_D2z_S_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1000_ab
   ************************************************************/
  Double I_ERI_Px_D2x_S_S_C1000_ab = I_ERI_D2x_Px_S_S_C1000_ab+ABX*I_ERI_Px_Px_S_S_C1000_ab;
  Double I_ERI_Py_D2x_S_S_C1000_ab = I_ERI_Dxy_Px_S_S_C1000_ab+ABX*I_ERI_Py_Px_S_S_C1000_ab;
  Double I_ERI_Pz_D2x_S_S_C1000_ab = I_ERI_Dxz_Px_S_S_C1000_ab+ABX*I_ERI_Pz_Px_S_S_C1000_ab;
  Double I_ERI_Px_Dxy_S_S_C1000_ab = I_ERI_Dxy_Px_S_S_C1000_ab+ABY*I_ERI_Px_Px_S_S_C1000_ab;
  Double I_ERI_Py_Dxy_S_S_C1000_ab = I_ERI_D2y_Px_S_S_C1000_ab+ABY*I_ERI_Py_Px_S_S_C1000_ab;
  Double I_ERI_Pz_Dxy_S_S_C1000_ab = I_ERI_Dyz_Px_S_S_C1000_ab+ABY*I_ERI_Pz_Px_S_S_C1000_ab;
  Double I_ERI_Px_Dxz_S_S_C1000_ab = I_ERI_Dxz_Px_S_S_C1000_ab+ABZ*I_ERI_Px_Px_S_S_C1000_ab;
  Double I_ERI_Py_Dxz_S_S_C1000_ab = I_ERI_Dyz_Px_S_S_C1000_ab+ABZ*I_ERI_Py_Px_S_S_C1000_ab;
  Double I_ERI_Pz_Dxz_S_S_C1000_ab = I_ERI_D2z_Px_S_S_C1000_ab+ABZ*I_ERI_Pz_Px_S_S_C1000_ab;
  Double I_ERI_Px_D2y_S_S_C1000_ab = I_ERI_Dxy_Py_S_S_C1000_ab+ABY*I_ERI_Px_Py_S_S_C1000_ab;
  Double I_ERI_Py_D2y_S_S_C1000_ab = I_ERI_D2y_Py_S_S_C1000_ab+ABY*I_ERI_Py_Py_S_S_C1000_ab;
  Double I_ERI_Pz_D2y_S_S_C1000_ab = I_ERI_Dyz_Py_S_S_C1000_ab+ABY*I_ERI_Pz_Py_S_S_C1000_ab;
  Double I_ERI_Px_Dyz_S_S_C1000_ab = I_ERI_Dxz_Py_S_S_C1000_ab+ABZ*I_ERI_Px_Py_S_S_C1000_ab;
  Double I_ERI_Py_Dyz_S_S_C1000_ab = I_ERI_Dyz_Py_S_S_C1000_ab+ABZ*I_ERI_Py_Py_S_S_C1000_ab;
  Double I_ERI_Pz_Dyz_S_S_C1000_ab = I_ERI_D2z_Py_S_S_C1000_ab+ABZ*I_ERI_Pz_Py_S_S_C1000_ab;
  Double I_ERI_Px_D2z_S_S_C1000_ab = I_ERI_Dxz_Pz_S_S_C1000_ab+ABZ*I_ERI_Px_Pz_S_S_C1000_ab;
  Double I_ERI_Py_D2z_S_S_C1000_ab = I_ERI_Dyz_Pz_S_S_C1000_ab+ABZ*I_ERI_Py_Pz_S_S_C1000_ab;
  Double I_ERI_Pz_D2z_S_S_C1000_ab = I_ERI_D2z_Pz_S_S_C1000_ab+ABZ*I_ERI_Pz_Pz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_ab
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001_ab = I_ERI_F3x_S_S_S_C1001_ab+ABX*I_ERI_D2x_S_S_S_C1001_ab;
  Double I_ERI_Dxy_Px_S_S_C1001_ab = I_ERI_F2xy_S_S_S_C1001_ab+ABX*I_ERI_Dxy_S_S_S_C1001_ab;
  Double I_ERI_Dxz_Px_S_S_C1001_ab = I_ERI_F2xz_S_S_S_C1001_ab+ABX*I_ERI_Dxz_S_S_S_C1001_ab;
  Double I_ERI_D2y_Px_S_S_C1001_ab = I_ERI_Fx2y_S_S_S_C1001_ab+ABX*I_ERI_D2y_S_S_S_C1001_ab;
  Double I_ERI_Dyz_Px_S_S_C1001_ab = I_ERI_Fxyz_S_S_S_C1001_ab+ABX*I_ERI_Dyz_S_S_S_C1001_ab;
  Double I_ERI_D2z_Px_S_S_C1001_ab = I_ERI_Fx2z_S_S_S_C1001_ab+ABX*I_ERI_D2z_S_S_S_C1001_ab;
  Double I_ERI_D2x_Py_S_S_C1001_ab = I_ERI_F2xy_S_S_S_C1001_ab+ABY*I_ERI_D2x_S_S_S_C1001_ab;
  Double I_ERI_Dxy_Py_S_S_C1001_ab = I_ERI_Fx2y_S_S_S_C1001_ab+ABY*I_ERI_Dxy_S_S_S_C1001_ab;
  Double I_ERI_Dxz_Py_S_S_C1001_ab = I_ERI_Fxyz_S_S_S_C1001_ab+ABY*I_ERI_Dxz_S_S_S_C1001_ab;
  Double I_ERI_D2y_Py_S_S_C1001_ab = I_ERI_F3y_S_S_S_C1001_ab+ABY*I_ERI_D2y_S_S_S_C1001_ab;
  Double I_ERI_Dyz_Py_S_S_C1001_ab = I_ERI_F2yz_S_S_S_C1001_ab+ABY*I_ERI_Dyz_S_S_S_C1001_ab;
  Double I_ERI_D2z_Py_S_S_C1001_ab = I_ERI_Fy2z_S_S_S_C1001_ab+ABY*I_ERI_D2z_S_S_S_C1001_ab;
  Double I_ERI_D2x_Pz_S_S_C1001_ab = I_ERI_F2xz_S_S_S_C1001_ab+ABZ*I_ERI_D2x_S_S_S_C1001_ab;
  Double I_ERI_Dxy_Pz_S_S_C1001_ab = I_ERI_Fxyz_S_S_S_C1001_ab+ABZ*I_ERI_Dxy_S_S_S_C1001_ab;
  Double I_ERI_Dxz_Pz_S_S_C1001_ab = I_ERI_Fx2z_S_S_S_C1001_ab+ABZ*I_ERI_Dxz_S_S_S_C1001_ab;
  Double I_ERI_D2y_Pz_S_S_C1001_ab = I_ERI_F2yz_S_S_S_C1001_ab+ABZ*I_ERI_D2y_S_S_S_C1001_ab;
  Double I_ERI_Dyz_Pz_S_S_C1001_ab = I_ERI_Fy2z_S_S_S_C1001_ab+ABZ*I_ERI_Dyz_S_S_S_C1001_ab;
  Double I_ERI_D2z_Pz_S_S_C1001_ab = I_ERI_F3z_S_S_S_C1001_ab+ABZ*I_ERI_D2z_S_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_ab
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1001_ab = I_ERI_G4x_S_S_S_C1001_ab+ABX*I_ERI_F3x_S_S_S_C1001_ab;
  Double I_ERI_F2xy_Px_S_S_C1001_ab = I_ERI_G3xy_S_S_S_C1001_ab+ABX*I_ERI_F2xy_S_S_S_C1001_ab;
  Double I_ERI_F2xz_Px_S_S_C1001_ab = I_ERI_G3xz_S_S_S_C1001_ab+ABX*I_ERI_F2xz_S_S_S_C1001_ab;
  Double I_ERI_Fx2y_Px_S_S_C1001_ab = I_ERI_G2x2y_S_S_S_C1001_ab+ABX*I_ERI_Fx2y_S_S_S_C1001_ab;
  Double I_ERI_Fxyz_Px_S_S_C1001_ab = I_ERI_G2xyz_S_S_S_C1001_ab+ABX*I_ERI_Fxyz_S_S_S_C1001_ab;
  Double I_ERI_Fx2z_Px_S_S_C1001_ab = I_ERI_G2x2z_S_S_S_C1001_ab+ABX*I_ERI_Fx2z_S_S_S_C1001_ab;
  Double I_ERI_F3y_Px_S_S_C1001_ab = I_ERI_Gx3y_S_S_S_C1001_ab+ABX*I_ERI_F3y_S_S_S_C1001_ab;
  Double I_ERI_F2yz_Px_S_S_C1001_ab = I_ERI_Gx2yz_S_S_S_C1001_ab+ABX*I_ERI_F2yz_S_S_S_C1001_ab;
  Double I_ERI_Fy2z_Px_S_S_C1001_ab = I_ERI_Gxy2z_S_S_S_C1001_ab+ABX*I_ERI_Fy2z_S_S_S_C1001_ab;
  Double I_ERI_F3z_Px_S_S_C1001_ab = I_ERI_Gx3z_S_S_S_C1001_ab+ABX*I_ERI_F3z_S_S_S_C1001_ab;
  Double I_ERI_F2xy_Py_S_S_C1001_ab = I_ERI_G2x2y_S_S_S_C1001_ab+ABY*I_ERI_F2xy_S_S_S_C1001_ab;
  Double I_ERI_F2xz_Py_S_S_C1001_ab = I_ERI_G2xyz_S_S_S_C1001_ab+ABY*I_ERI_F2xz_S_S_S_C1001_ab;
  Double I_ERI_Fx2y_Py_S_S_C1001_ab = I_ERI_Gx3y_S_S_S_C1001_ab+ABY*I_ERI_Fx2y_S_S_S_C1001_ab;
  Double I_ERI_Fxyz_Py_S_S_C1001_ab = I_ERI_Gx2yz_S_S_S_C1001_ab+ABY*I_ERI_Fxyz_S_S_S_C1001_ab;
  Double I_ERI_Fx2z_Py_S_S_C1001_ab = I_ERI_Gxy2z_S_S_S_C1001_ab+ABY*I_ERI_Fx2z_S_S_S_C1001_ab;
  Double I_ERI_F3y_Py_S_S_C1001_ab = I_ERI_G4y_S_S_S_C1001_ab+ABY*I_ERI_F3y_S_S_S_C1001_ab;
  Double I_ERI_F2yz_Py_S_S_C1001_ab = I_ERI_G3yz_S_S_S_C1001_ab+ABY*I_ERI_F2yz_S_S_S_C1001_ab;
  Double I_ERI_Fy2z_Py_S_S_C1001_ab = I_ERI_G2y2z_S_S_S_C1001_ab+ABY*I_ERI_Fy2z_S_S_S_C1001_ab;
  Double I_ERI_F3z_Py_S_S_C1001_ab = I_ERI_Gy3z_S_S_S_C1001_ab+ABY*I_ERI_F3z_S_S_S_C1001_ab;
  Double I_ERI_F2xz_Pz_S_S_C1001_ab = I_ERI_G2x2z_S_S_S_C1001_ab+ABZ*I_ERI_F2xz_S_S_S_C1001_ab;
  Double I_ERI_Fxyz_Pz_S_S_C1001_ab = I_ERI_Gxy2z_S_S_S_C1001_ab+ABZ*I_ERI_Fxyz_S_S_S_C1001_ab;
  Double I_ERI_Fx2z_Pz_S_S_C1001_ab = I_ERI_Gx3z_S_S_S_C1001_ab+ABZ*I_ERI_Fx2z_S_S_S_C1001_ab;
  Double I_ERI_F2yz_Pz_S_S_C1001_ab = I_ERI_G2y2z_S_S_S_C1001_ab+ABZ*I_ERI_F2yz_S_S_S_C1001_ab;
  Double I_ERI_Fy2z_Pz_S_S_C1001_ab = I_ERI_Gy3z_S_S_S_C1001_ab+ABZ*I_ERI_Fy2z_S_S_S_C1001_ab;
  Double I_ERI_F3z_Pz_S_S_C1001_ab = I_ERI_G4z_S_S_S_C1001_ab+ABZ*I_ERI_F3z_S_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_ab
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_C1001_ab = I_ERI_F3x_Px_S_S_C1001_ab+ABX*I_ERI_D2x_Px_S_S_C1001_ab;
  Double I_ERI_Dxy_D2x_S_S_C1001_ab = I_ERI_F2xy_Px_S_S_C1001_ab+ABX*I_ERI_Dxy_Px_S_S_C1001_ab;
  Double I_ERI_Dxz_D2x_S_S_C1001_ab = I_ERI_F2xz_Px_S_S_C1001_ab+ABX*I_ERI_Dxz_Px_S_S_C1001_ab;
  Double I_ERI_D2y_D2x_S_S_C1001_ab = I_ERI_Fx2y_Px_S_S_C1001_ab+ABX*I_ERI_D2y_Px_S_S_C1001_ab;
  Double I_ERI_Dyz_D2x_S_S_C1001_ab = I_ERI_Fxyz_Px_S_S_C1001_ab+ABX*I_ERI_Dyz_Px_S_S_C1001_ab;
  Double I_ERI_D2z_D2x_S_S_C1001_ab = I_ERI_Fx2z_Px_S_S_C1001_ab+ABX*I_ERI_D2z_Px_S_S_C1001_ab;
  Double I_ERI_D2x_Dxy_S_S_C1001_ab = I_ERI_F2xy_Px_S_S_C1001_ab+ABY*I_ERI_D2x_Px_S_S_C1001_ab;
  Double I_ERI_Dxy_Dxy_S_S_C1001_ab = I_ERI_Fx2y_Px_S_S_C1001_ab+ABY*I_ERI_Dxy_Px_S_S_C1001_ab;
  Double I_ERI_Dxz_Dxy_S_S_C1001_ab = I_ERI_Fxyz_Px_S_S_C1001_ab+ABY*I_ERI_Dxz_Px_S_S_C1001_ab;
  Double I_ERI_D2y_Dxy_S_S_C1001_ab = I_ERI_F3y_Px_S_S_C1001_ab+ABY*I_ERI_D2y_Px_S_S_C1001_ab;
  Double I_ERI_Dyz_Dxy_S_S_C1001_ab = I_ERI_F2yz_Px_S_S_C1001_ab+ABY*I_ERI_Dyz_Px_S_S_C1001_ab;
  Double I_ERI_D2z_Dxy_S_S_C1001_ab = I_ERI_Fy2z_Px_S_S_C1001_ab+ABY*I_ERI_D2z_Px_S_S_C1001_ab;
  Double I_ERI_D2x_Dxz_S_S_C1001_ab = I_ERI_F2xz_Px_S_S_C1001_ab+ABZ*I_ERI_D2x_Px_S_S_C1001_ab;
  Double I_ERI_Dxy_Dxz_S_S_C1001_ab = I_ERI_Fxyz_Px_S_S_C1001_ab+ABZ*I_ERI_Dxy_Px_S_S_C1001_ab;
  Double I_ERI_Dxz_Dxz_S_S_C1001_ab = I_ERI_Fx2z_Px_S_S_C1001_ab+ABZ*I_ERI_Dxz_Px_S_S_C1001_ab;
  Double I_ERI_D2y_Dxz_S_S_C1001_ab = I_ERI_F2yz_Px_S_S_C1001_ab+ABZ*I_ERI_D2y_Px_S_S_C1001_ab;
  Double I_ERI_Dyz_Dxz_S_S_C1001_ab = I_ERI_Fy2z_Px_S_S_C1001_ab+ABZ*I_ERI_Dyz_Px_S_S_C1001_ab;
  Double I_ERI_D2z_Dxz_S_S_C1001_ab = I_ERI_F3z_Px_S_S_C1001_ab+ABZ*I_ERI_D2z_Px_S_S_C1001_ab;
  Double I_ERI_D2x_D2y_S_S_C1001_ab = I_ERI_F2xy_Py_S_S_C1001_ab+ABY*I_ERI_D2x_Py_S_S_C1001_ab;
  Double I_ERI_Dxy_D2y_S_S_C1001_ab = I_ERI_Fx2y_Py_S_S_C1001_ab+ABY*I_ERI_Dxy_Py_S_S_C1001_ab;
  Double I_ERI_Dxz_D2y_S_S_C1001_ab = I_ERI_Fxyz_Py_S_S_C1001_ab+ABY*I_ERI_Dxz_Py_S_S_C1001_ab;
  Double I_ERI_D2y_D2y_S_S_C1001_ab = I_ERI_F3y_Py_S_S_C1001_ab+ABY*I_ERI_D2y_Py_S_S_C1001_ab;
  Double I_ERI_Dyz_D2y_S_S_C1001_ab = I_ERI_F2yz_Py_S_S_C1001_ab+ABY*I_ERI_Dyz_Py_S_S_C1001_ab;
  Double I_ERI_D2z_D2y_S_S_C1001_ab = I_ERI_Fy2z_Py_S_S_C1001_ab+ABY*I_ERI_D2z_Py_S_S_C1001_ab;
  Double I_ERI_D2x_Dyz_S_S_C1001_ab = I_ERI_F2xz_Py_S_S_C1001_ab+ABZ*I_ERI_D2x_Py_S_S_C1001_ab;
  Double I_ERI_Dxy_Dyz_S_S_C1001_ab = I_ERI_Fxyz_Py_S_S_C1001_ab+ABZ*I_ERI_Dxy_Py_S_S_C1001_ab;
  Double I_ERI_Dxz_Dyz_S_S_C1001_ab = I_ERI_Fx2z_Py_S_S_C1001_ab+ABZ*I_ERI_Dxz_Py_S_S_C1001_ab;
  Double I_ERI_D2y_Dyz_S_S_C1001_ab = I_ERI_F2yz_Py_S_S_C1001_ab+ABZ*I_ERI_D2y_Py_S_S_C1001_ab;
  Double I_ERI_Dyz_Dyz_S_S_C1001_ab = I_ERI_Fy2z_Py_S_S_C1001_ab+ABZ*I_ERI_Dyz_Py_S_S_C1001_ab;
  Double I_ERI_D2z_Dyz_S_S_C1001_ab = I_ERI_F3z_Py_S_S_C1001_ab+ABZ*I_ERI_D2z_Py_S_S_C1001_ab;
  Double I_ERI_D2x_D2z_S_S_C1001_ab = I_ERI_F2xz_Pz_S_S_C1001_ab+ABZ*I_ERI_D2x_Pz_S_S_C1001_ab;
  Double I_ERI_Dxy_D2z_S_S_C1001_ab = I_ERI_Fxyz_Pz_S_S_C1001_ab+ABZ*I_ERI_Dxy_Pz_S_S_C1001_ab;
  Double I_ERI_Dxz_D2z_S_S_C1001_ab = I_ERI_Fx2z_Pz_S_S_C1001_ab+ABZ*I_ERI_Dxz_Pz_S_S_C1001_ab;
  Double I_ERI_D2y_D2z_S_S_C1001_ab = I_ERI_F2yz_Pz_S_S_C1001_ab+ABZ*I_ERI_D2y_Pz_S_S_C1001_ab;
  Double I_ERI_Dyz_D2z_S_S_C1001_ab = I_ERI_Fy2z_Pz_S_S_C1001_ab+ABZ*I_ERI_Dyz_Pz_S_S_C1001_ab;
  Double I_ERI_D2z_D2z_S_S_C1001_ab = I_ERI_F3z_Pz_S_S_C1001_ab+ABZ*I_ERI_D2z_Pz_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000_ac
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1000_ac = I_ERI_D2x_S_Px_S_C1000_ac+ABX*I_ERI_Px_S_Px_S_C1000_ac;
  Double I_ERI_Py_Px_Px_S_C1000_ac = I_ERI_Dxy_S_Px_S_C1000_ac+ABX*I_ERI_Py_S_Px_S_C1000_ac;
  Double I_ERI_Pz_Px_Px_S_C1000_ac = I_ERI_Dxz_S_Px_S_C1000_ac+ABX*I_ERI_Pz_S_Px_S_C1000_ac;
  Double I_ERI_Px_Py_Px_S_C1000_ac = I_ERI_Dxy_S_Px_S_C1000_ac+ABY*I_ERI_Px_S_Px_S_C1000_ac;
  Double I_ERI_Py_Py_Px_S_C1000_ac = I_ERI_D2y_S_Px_S_C1000_ac+ABY*I_ERI_Py_S_Px_S_C1000_ac;
  Double I_ERI_Pz_Py_Px_S_C1000_ac = I_ERI_Dyz_S_Px_S_C1000_ac+ABY*I_ERI_Pz_S_Px_S_C1000_ac;
  Double I_ERI_Px_Pz_Px_S_C1000_ac = I_ERI_Dxz_S_Px_S_C1000_ac+ABZ*I_ERI_Px_S_Px_S_C1000_ac;
  Double I_ERI_Py_Pz_Px_S_C1000_ac = I_ERI_Dyz_S_Px_S_C1000_ac+ABZ*I_ERI_Py_S_Px_S_C1000_ac;
  Double I_ERI_Pz_Pz_Px_S_C1000_ac = I_ERI_D2z_S_Px_S_C1000_ac+ABZ*I_ERI_Pz_S_Px_S_C1000_ac;
  Double I_ERI_Px_Px_Py_S_C1000_ac = I_ERI_D2x_S_Py_S_C1000_ac+ABX*I_ERI_Px_S_Py_S_C1000_ac;
  Double I_ERI_Py_Px_Py_S_C1000_ac = I_ERI_Dxy_S_Py_S_C1000_ac+ABX*I_ERI_Py_S_Py_S_C1000_ac;
  Double I_ERI_Pz_Px_Py_S_C1000_ac = I_ERI_Dxz_S_Py_S_C1000_ac+ABX*I_ERI_Pz_S_Py_S_C1000_ac;
  Double I_ERI_Px_Py_Py_S_C1000_ac = I_ERI_Dxy_S_Py_S_C1000_ac+ABY*I_ERI_Px_S_Py_S_C1000_ac;
  Double I_ERI_Py_Py_Py_S_C1000_ac = I_ERI_D2y_S_Py_S_C1000_ac+ABY*I_ERI_Py_S_Py_S_C1000_ac;
  Double I_ERI_Pz_Py_Py_S_C1000_ac = I_ERI_Dyz_S_Py_S_C1000_ac+ABY*I_ERI_Pz_S_Py_S_C1000_ac;
  Double I_ERI_Px_Pz_Py_S_C1000_ac = I_ERI_Dxz_S_Py_S_C1000_ac+ABZ*I_ERI_Px_S_Py_S_C1000_ac;
  Double I_ERI_Py_Pz_Py_S_C1000_ac = I_ERI_Dyz_S_Py_S_C1000_ac+ABZ*I_ERI_Py_S_Py_S_C1000_ac;
  Double I_ERI_Pz_Pz_Py_S_C1000_ac = I_ERI_D2z_S_Py_S_C1000_ac+ABZ*I_ERI_Pz_S_Py_S_C1000_ac;
  Double I_ERI_Px_Px_Pz_S_C1000_ac = I_ERI_D2x_S_Pz_S_C1000_ac+ABX*I_ERI_Px_S_Pz_S_C1000_ac;
  Double I_ERI_Py_Px_Pz_S_C1000_ac = I_ERI_Dxy_S_Pz_S_C1000_ac+ABX*I_ERI_Py_S_Pz_S_C1000_ac;
  Double I_ERI_Pz_Px_Pz_S_C1000_ac = I_ERI_Dxz_S_Pz_S_C1000_ac+ABX*I_ERI_Pz_S_Pz_S_C1000_ac;
  Double I_ERI_Px_Py_Pz_S_C1000_ac = I_ERI_Dxy_S_Pz_S_C1000_ac+ABY*I_ERI_Px_S_Pz_S_C1000_ac;
  Double I_ERI_Py_Py_Pz_S_C1000_ac = I_ERI_D2y_S_Pz_S_C1000_ac+ABY*I_ERI_Py_S_Pz_S_C1000_ac;
  Double I_ERI_Pz_Py_Pz_S_C1000_ac = I_ERI_Dyz_S_Pz_S_C1000_ac+ABY*I_ERI_Pz_S_Pz_S_C1000_ac;
  Double I_ERI_Px_Pz_Pz_S_C1000_ac = I_ERI_Dxz_S_Pz_S_C1000_ac+ABZ*I_ERI_Px_S_Pz_S_C1000_ac;
  Double I_ERI_Py_Pz_Pz_S_C1000_ac = I_ERI_Dyz_S_Pz_S_C1000_ac+ABZ*I_ERI_Py_S_Pz_S_C1000_ac;
  Double I_ERI_Pz_Pz_Pz_S_C1000_ac = I_ERI_D2z_S_Pz_S_C1000_ac+ABZ*I_ERI_Pz_S_Pz_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001_ac
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1001_ac = I_ERI_F3x_S_Px_S_C1001_ac+ABX*I_ERI_D2x_S_Px_S_C1001_ac;
  Double I_ERI_Dxy_Px_Px_S_C1001_ac = I_ERI_F2xy_S_Px_S_C1001_ac+ABX*I_ERI_Dxy_S_Px_S_C1001_ac;
  Double I_ERI_Dxz_Px_Px_S_C1001_ac = I_ERI_F2xz_S_Px_S_C1001_ac+ABX*I_ERI_Dxz_S_Px_S_C1001_ac;
  Double I_ERI_D2y_Px_Px_S_C1001_ac = I_ERI_Fx2y_S_Px_S_C1001_ac+ABX*I_ERI_D2y_S_Px_S_C1001_ac;
  Double I_ERI_Dyz_Px_Px_S_C1001_ac = I_ERI_Fxyz_S_Px_S_C1001_ac+ABX*I_ERI_Dyz_S_Px_S_C1001_ac;
  Double I_ERI_D2z_Px_Px_S_C1001_ac = I_ERI_Fx2z_S_Px_S_C1001_ac+ABX*I_ERI_D2z_S_Px_S_C1001_ac;
  Double I_ERI_D2x_Py_Px_S_C1001_ac = I_ERI_F2xy_S_Px_S_C1001_ac+ABY*I_ERI_D2x_S_Px_S_C1001_ac;
  Double I_ERI_Dxy_Py_Px_S_C1001_ac = I_ERI_Fx2y_S_Px_S_C1001_ac+ABY*I_ERI_Dxy_S_Px_S_C1001_ac;
  Double I_ERI_Dxz_Py_Px_S_C1001_ac = I_ERI_Fxyz_S_Px_S_C1001_ac+ABY*I_ERI_Dxz_S_Px_S_C1001_ac;
  Double I_ERI_D2y_Py_Px_S_C1001_ac = I_ERI_F3y_S_Px_S_C1001_ac+ABY*I_ERI_D2y_S_Px_S_C1001_ac;
  Double I_ERI_Dyz_Py_Px_S_C1001_ac = I_ERI_F2yz_S_Px_S_C1001_ac+ABY*I_ERI_Dyz_S_Px_S_C1001_ac;
  Double I_ERI_D2z_Py_Px_S_C1001_ac = I_ERI_Fy2z_S_Px_S_C1001_ac+ABY*I_ERI_D2z_S_Px_S_C1001_ac;
  Double I_ERI_D2x_Pz_Px_S_C1001_ac = I_ERI_F2xz_S_Px_S_C1001_ac+ABZ*I_ERI_D2x_S_Px_S_C1001_ac;
  Double I_ERI_Dxy_Pz_Px_S_C1001_ac = I_ERI_Fxyz_S_Px_S_C1001_ac+ABZ*I_ERI_Dxy_S_Px_S_C1001_ac;
  Double I_ERI_Dxz_Pz_Px_S_C1001_ac = I_ERI_Fx2z_S_Px_S_C1001_ac+ABZ*I_ERI_Dxz_S_Px_S_C1001_ac;
  Double I_ERI_D2y_Pz_Px_S_C1001_ac = I_ERI_F2yz_S_Px_S_C1001_ac+ABZ*I_ERI_D2y_S_Px_S_C1001_ac;
  Double I_ERI_Dyz_Pz_Px_S_C1001_ac = I_ERI_Fy2z_S_Px_S_C1001_ac+ABZ*I_ERI_Dyz_S_Px_S_C1001_ac;
  Double I_ERI_D2z_Pz_Px_S_C1001_ac = I_ERI_F3z_S_Px_S_C1001_ac+ABZ*I_ERI_D2z_S_Px_S_C1001_ac;
  Double I_ERI_D2x_Px_Py_S_C1001_ac = I_ERI_F3x_S_Py_S_C1001_ac+ABX*I_ERI_D2x_S_Py_S_C1001_ac;
  Double I_ERI_Dxy_Px_Py_S_C1001_ac = I_ERI_F2xy_S_Py_S_C1001_ac+ABX*I_ERI_Dxy_S_Py_S_C1001_ac;
  Double I_ERI_Dxz_Px_Py_S_C1001_ac = I_ERI_F2xz_S_Py_S_C1001_ac+ABX*I_ERI_Dxz_S_Py_S_C1001_ac;
  Double I_ERI_D2y_Px_Py_S_C1001_ac = I_ERI_Fx2y_S_Py_S_C1001_ac+ABX*I_ERI_D2y_S_Py_S_C1001_ac;
  Double I_ERI_Dyz_Px_Py_S_C1001_ac = I_ERI_Fxyz_S_Py_S_C1001_ac+ABX*I_ERI_Dyz_S_Py_S_C1001_ac;
  Double I_ERI_D2z_Px_Py_S_C1001_ac = I_ERI_Fx2z_S_Py_S_C1001_ac+ABX*I_ERI_D2z_S_Py_S_C1001_ac;
  Double I_ERI_D2x_Py_Py_S_C1001_ac = I_ERI_F2xy_S_Py_S_C1001_ac+ABY*I_ERI_D2x_S_Py_S_C1001_ac;
  Double I_ERI_Dxy_Py_Py_S_C1001_ac = I_ERI_Fx2y_S_Py_S_C1001_ac+ABY*I_ERI_Dxy_S_Py_S_C1001_ac;
  Double I_ERI_Dxz_Py_Py_S_C1001_ac = I_ERI_Fxyz_S_Py_S_C1001_ac+ABY*I_ERI_Dxz_S_Py_S_C1001_ac;
  Double I_ERI_D2y_Py_Py_S_C1001_ac = I_ERI_F3y_S_Py_S_C1001_ac+ABY*I_ERI_D2y_S_Py_S_C1001_ac;
  Double I_ERI_Dyz_Py_Py_S_C1001_ac = I_ERI_F2yz_S_Py_S_C1001_ac+ABY*I_ERI_Dyz_S_Py_S_C1001_ac;
  Double I_ERI_D2z_Py_Py_S_C1001_ac = I_ERI_Fy2z_S_Py_S_C1001_ac+ABY*I_ERI_D2z_S_Py_S_C1001_ac;
  Double I_ERI_D2x_Pz_Py_S_C1001_ac = I_ERI_F2xz_S_Py_S_C1001_ac+ABZ*I_ERI_D2x_S_Py_S_C1001_ac;
  Double I_ERI_Dxy_Pz_Py_S_C1001_ac = I_ERI_Fxyz_S_Py_S_C1001_ac+ABZ*I_ERI_Dxy_S_Py_S_C1001_ac;
  Double I_ERI_Dxz_Pz_Py_S_C1001_ac = I_ERI_Fx2z_S_Py_S_C1001_ac+ABZ*I_ERI_Dxz_S_Py_S_C1001_ac;
  Double I_ERI_D2y_Pz_Py_S_C1001_ac = I_ERI_F2yz_S_Py_S_C1001_ac+ABZ*I_ERI_D2y_S_Py_S_C1001_ac;
  Double I_ERI_Dyz_Pz_Py_S_C1001_ac = I_ERI_Fy2z_S_Py_S_C1001_ac+ABZ*I_ERI_Dyz_S_Py_S_C1001_ac;
  Double I_ERI_D2z_Pz_Py_S_C1001_ac = I_ERI_F3z_S_Py_S_C1001_ac+ABZ*I_ERI_D2z_S_Py_S_C1001_ac;
  Double I_ERI_D2x_Px_Pz_S_C1001_ac = I_ERI_F3x_S_Pz_S_C1001_ac+ABX*I_ERI_D2x_S_Pz_S_C1001_ac;
  Double I_ERI_Dxy_Px_Pz_S_C1001_ac = I_ERI_F2xy_S_Pz_S_C1001_ac+ABX*I_ERI_Dxy_S_Pz_S_C1001_ac;
  Double I_ERI_Dxz_Px_Pz_S_C1001_ac = I_ERI_F2xz_S_Pz_S_C1001_ac+ABX*I_ERI_Dxz_S_Pz_S_C1001_ac;
  Double I_ERI_D2y_Px_Pz_S_C1001_ac = I_ERI_Fx2y_S_Pz_S_C1001_ac+ABX*I_ERI_D2y_S_Pz_S_C1001_ac;
  Double I_ERI_Dyz_Px_Pz_S_C1001_ac = I_ERI_Fxyz_S_Pz_S_C1001_ac+ABX*I_ERI_Dyz_S_Pz_S_C1001_ac;
  Double I_ERI_D2z_Px_Pz_S_C1001_ac = I_ERI_Fx2z_S_Pz_S_C1001_ac+ABX*I_ERI_D2z_S_Pz_S_C1001_ac;
  Double I_ERI_D2x_Py_Pz_S_C1001_ac = I_ERI_F2xy_S_Pz_S_C1001_ac+ABY*I_ERI_D2x_S_Pz_S_C1001_ac;
  Double I_ERI_Dxy_Py_Pz_S_C1001_ac = I_ERI_Fx2y_S_Pz_S_C1001_ac+ABY*I_ERI_Dxy_S_Pz_S_C1001_ac;
  Double I_ERI_Dxz_Py_Pz_S_C1001_ac = I_ERI_Fxyz_S_Pz_S_C1001_ac+ABY*I_ERI_Dxz_S_Pz_S_C1001_ac;
  Double I_ERI_D2y_Py_Pz_S_C1001_ac = I_ERI_F3y_S_Pz_S_C1001_ac+ABY*I_ERI_D2y_S_Pz_S_C1001_ac;
  Double I_ERI_Dyz_Py_Pz_S_C1001_ac = I_ERI_F2yz_S_Pz_S_C1001_ac+ABY*I_ERI_Dyz_S_Pz_S_C1001_ac;
  Double I_ERI_D2z_Py_Pz_S_C1001_ac = I_ERI_Fy2z_S_Pz_S_C1001_ac+ABY*I_ERI_D2z_S_Pz_S_C1001_ac;
  Double I_ERI_D2x_Pz_Pz_S_C1001_ac = I_ERI_F2xz_S_Pz_S_C1001_ac+ABZ*I_ERI_D2x_S_Pz_S_C1001_ac;
  Double I_ERI_Dxy_Pz_Pz_S_C1001_ac = I_ERI_Fxyz_S_Pz_S_C1001_ac+ABZ*I_ERI_Dxy_S_Pz_S_C1001_ac;
  Double I_ERI_Dxz_Pz_Pz_S_C1001_ac = I_ERI_Fx2z_S_Pz_S_C1001_ac+ABZ*I_ERI_Dxz_S_Pz_S_C1001_ac;
  Double I_ERI_D2y_Pz_Pz_S_C1001_ac = I_ERI_F2yz_S_Pz_S_C1001_ac+ABZ*I_ERI_D2y_S_Pz_S_C1001_ac;
  Double I_ERI_Dyz_Pz_Pz_S_C1001_ac = I_ERI_Fy2z_S_Pz_S_C1001_ac+ABZ*I_ERI_Dyz_S_Pz_S_C1001_ac;
  Double I_ERI_D2z_Pz_Pz_S_C1001_ac = I_ERI_F3z_S_Pz_S_C1001_ac+ABZ*I_ERI_D2z_S_Pz_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_bb
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1_bb = I_ERI_D2x_S_S_S_C1_bb+ABX*I_ERI_Px_S_S_S_C1_bb;
  Double I_ERI_Py_Px_S_S_C1_bb = I_ERI_Dxy_S_S_S_C1_bb+ABX*I_ERI_Py_S_S_S_C1_bb;
  Double I_ERI_Pz_Px_S_S_C1_bb = I_ERI_Dxz_S_S_S_C1_bb+ABX*I_ERI_Pz_S_S_S_C1_bb;
  Double I_ERI_Px_Py_S_S_C1_bb = I_ERI_Dxy_S_S_S_C1_bb+ABY*I_ERI_Px_S_S_S_C1_bb;
  Double I_ERI_Py_Py_S_S_C1_bb = I_ERI_D2y_S_S_S_C1_bb+ABY*I_ERI_Py_S_S_S_C1_bb;
  Double I_ERI_Pz_Py_S_S_C1_bb = I_ERI_Dyz_S_S_S_C1_bb+ABY*I_ERI_Pz_S_S_S_C1_bb;
  Double I_ERI_Px_Pz_S_S_C1_bb = I_ERI_Dxz_S_S_S_C1_bb+ABZ*I_ERI_Px_S_S_S_C1_bb;
  Double I_ERI_Py_Pz_S_S_C1_bb = I_ERI_Dyz_S_S_S_C1_bb+ABZ*I_ERI_Py_S_S_S_C1_bb;
  Double I_ERI_Pz_Pz_S_S_C1_bb = I_ERI_D2z_S_S_S_C1_bb+ABZ*I_ERI_Pz_S_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 4 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_bb
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1_bb = I_ERI_F3x_S_S_S_C1_bb+ABX*I_ERI_D2x_S_S_S_C1_bb;
  Double I_ERI_Dxy_Px_S_S_C1_bb = I_ERI_F2xy_S_S_S_C1_bb+ABX*I_ERI_Dxy_S_S_S_C1_bb;
  Double I_ERI_Dxz_Px_S_S_C1_bb = I_ERI_F2xz_S_S_S_C1_bb+ABX*I_ERI_Dxz_S_S_S_C1_bb;
  Double I_ERI_D2y_Px_S_S_C1_bb = I_ERI_Fx2y_S_S_S_C1_bb+ABX*I_ERI_D2y_S_S_S_C1_bb;
  Double I_ERI_Dyz_Px_S_S_C1_bb = I_ERI_Fxyz_S_S_S_C1_bb+ABX*I_ERI_Dyz_S_S_S_C1_bb;
  Double I_ERI_D2z_Px_S_S_C1_bb = I_ERI_Fx2z_S_S_S_C1_bb+ABX*I_ERI_D2z_S_S_S_C1_bb;
  Double I_ERI_Dxy_Py_S_S_C1_bb = I_ERI_Fx2y_S_S_S_C1_bb+ABY*I_ERI_Dxy_S_S_S_C1_bb;
  Double I_ERI_Dxz_Py_S_S_C1_bb = I_ERI_Fxyz_S_S_S_C1_bb+ABY*I_ERI_Dxz_S_S_S_C1_bb;
  Double I_ERI_D2y_Py_S_S_C1_bb = I_ERI_F3y_S_S_S_C1_bb+ABY*I_ERI_D2y_S_S_S_C1_bb;
  Double I_ERI_Dyz_Py_S_S_C1_bb = I_ERI_F2yz_S_S_S_C1_bb+ABY*I_ERI_Dyz_S_S_S_C1_bb;
  Double I_ERI_D2z_Py_S_S_C1_bb = I_ERI_Fy2z_S_S_S_C1_bb+ABY*I_ERI_D2z_S_S_S_C1_bb;
  Double I_ERI_Dxz_Pz_S_S_C1_bb = I_ERI_Fx2z_S_S_S_C1_bb+ABZ*I_ERI_Dxz_S_S_S_C1_bb;
  Double I_ERI_Dyz_Pz_S_S_C1_bb = I_ERI_Fy2z_S_S_S_C1_bb+ABZ*I_ERI_Dyz_S_S_S_C1_bb;
  Double I_ERI_D2z_Pz_S_S_C1_bb = I_ERI_F3z_S_S_S_C1_bb+ABZ*I_ERI_D2z_S_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1_bb
   ************************************************************/
  Double I_ERI_Px_D2x_S_S_C1_bb = I_ERI_D2x_Px_S_S_C1_bb+ABX*I_ERI_Px_Px_S_S_C1_bb;
  Double I_ERI_Py_D2x_S_S_C1_bb = I_ERI_Dxy_Px_S_S_C1_bb+ABX*I_ERI_Py_Px_S_S_C1_bb;
  Double I_ERI_Pz_D2x_S_S_C1_bb = I_ERI_Dxz_Px_S_S_C1_bb+ABX*I_ERI_Pz_Px_S_S_C1_bb;
  Double I_ERI_Px_Dxy_S_S_C1_bb = I_ERI_Dxy_Px_S_S_C1_bb+ABY*I_ERI_Px_Px_S_S_C1_bb;
  Double I_ERI_Py_Dxy_S_S_C1_bb = I_ERI_D2y_Px_S_S_C1_bb+ABY*I_ERI_Py_Px_S_S_C1_bb;
  Double I_ERI_Pz_Dxy_S_S_C1_bb = I_ERI_Dyz_Px_S_S_C1_bb+ABY*I_ERI_Pz_Px_S_S_C1_bb;
  Double I_ERI_Px_Dxz_S_S_C1_bb = I_ERI_Dxz_Px_S_S_C1_bb+ABZ*I_ERI_Px_Px_S_S_C1_bb;
  Double I_ERI_Py_Dxz_S_S_C1_bb = I_ERI_Dyz_Px_S_S_C1_bb+ABZ*I_ERI_Py_Px_S_S_C1_bb;
  Double I_ERI_Pz_Dxz_S_S_C1_bb = I_ERI_D2z_Px_S_S_C1_bb+ABZ*I_ERI_Pz_Px_S_S_C1_bb;
  Double I_ERI_Px_D2y_S_S_C1_bb = I_ERI_Dxy_Py_S_S_C1_bb+ABY*I_ERI_Px_Py_S_S_C1_bb;
  Double I_ERI_Py_D2y_S_S_C1_bb = I_ERI_D2y_Py_S_S_C1_bb+ABY*I_ERI_Py_Py_S_S_C1_bb;
  Double I_ERI_Pz_D2y_S_S_C1_bb = I_ERI_Dyz_Py_S_S_C1_bb+ABY*I_ERI_Pz_Py_S_S_C1_bb;
  Double I_ERI_Px_Dyz_S_S_C1_bb = I_ERI_Dxz_Py_S_S_C1_bb+ABZ*I_ERI_Px_Py_S_S_C1_bb;
  Double I_ERI_Py_Dyz_S_S_C1_bb = I_ERI_Dyz_Py_S_S_C1_bb+ABZ*I_ERI_Py_Py_S_S_C1_bb;
  Double I_ERI_Pz_Dyz_S_S_C1_bb = I_ERI_D2z_Py_S_S_C1_bb+ABZ*I_ERI_Pz_Py_S_S_C1_bb;
  Double I_ERI_Px_D2z_S_S_C1_bb = I_ERI_Dxz_Pz_S_S_C1_bb+ABZ*I_ERI_Px_Pz_S_S_C1_bb;
  Double I_ERI_Py_D2z_S_S_C1_bb = I_ERI_Dyz_Pz_S_S_C1_bb+ABZ*I_ERI_Py_Pz_S_S_C1_bb;
  Double I_ERI_Pz_D2z_S_S_C1_bb = I_ERI_D2z_Pz_S_S_C1_bb+ABZ*I_ERI_Pz_Pz_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1001_bb = I_ERI_D2x_S_S_S_C1001_bb+ABX*I_ERI_Px_S_S_S_C1001_bb;
  Double I_ERI_Py_Px_S_S_C1001_bb = I_ERI_Dxy_S_S_S_C1001_bb+ABX*I_ERI_Py_S_S_S_C1001_bb;
  Double I_ERI_Pz_Px_S_S_C1001_bb = I_ERI_Dxz_S_S_S_C1001_bb+ABX*I_ERI_Pz_S_S_S_C1001_bb;
  Double I_ERI_Px_Py_S_S_C1001_bb = I_ERI_Dxy_S_S_S_C1001_bb+ABY*I_ERI_Px_S_S_S_C1001_bb;
  Double I_ERI_Py_Py_S_S_C1001_bb = I_ERI_D2y_S_S_S_C1001_bb+ABY*I_ERI_Py_S_S_S_C1001_bb;
  Double I_ERI_Pz_Py_S_S_C1001_bb = I_ERI_Dyz_S_S_S_C1001_bb+ABY*I_ERI_Pz_S_S_S_C1001_bb;
  Double I_ERI_Px_Pz_S_S_C1001_bb = I_ERI_Dxz_S_S_S_C1001_bb+ABZ*I_ERI_Px_S_S_S_C1001_bb;
  Double I_ERI_Py_Pz_S_S_C1001_bb = I_ERI_Dyz_S_S_S_C1001_bb+ABZ*I_ERI_Py_S_S_S_C1001_bb;
  Double I_ERI_Pz_Pz_S_S_C1001_bb = I_ERI_D2z_S_S_S_C1001_bb+ABZ*I_ERI_Pz_S_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001_bb = I_ERI_F3x_S_S_S_C1001_bb+ABX*I_ERI_D2x_S_S_S_C1001_bb;
  Double I_ERI_Dxy_Px_S_S_C1001_bb = I_ERI_F2xy_S_S_S_C1001_bb+ABX*I_ERI_Dxy_S_S_S_C1001_bb;
  Double I_ERI_Dxz_Px_S_S_C1001_bb = I_ERI_F2xz_S_S_S_C1001_bb+ABX*I_ERI_Dxz_S_S_S_C1001_bb;
  Double I_ERI_D2y_Px_S_S_C1001_bb = I_ERI_Fx2y_S_S_S_C1001_bb+ABX*I_ERI_D2y_S_S_S_C1001_bb;
  Double I_ERI_Dyz_Px_S_S_C1001_bb = I_ERI_Fxyz_S_S_S_C1001_bb+ABX*I_ERI_Dyz_S_S_S_C1001_bb;
  Double I_ERI_D2z_Px_S_S_C1001_bb = I_ERI_Fx2z_S_S_S_C1001_bb+ABX*I_ERI_D2z_S_S_S_C1001_bb;
  Double I_ERI_D2x_Py_S_S_C1001_bb = I_ERI_F2xy_S_S_S_C1001_bb+ABY*I_ERI_D2x_S_S_S_C1001_bb;
  Double I_ERI_Dxy_Py_S_S_C1001_bb = I_ERI_Fx2y_S_S_S_C1001_bb+ABY*I_ERI_Dxy_S_S_S_C1001_bb;
  Double I_ERI_Dxz_Py_S_S_C1001_bb = I_ERI_Fxyz_S_S_S_C1001_bb+ABY*I_ERI_Dxz_S_S_S_C1001_bb;
  Double I_ERI_D2y_Py_S_S_C1001_bb = I_ERI_F3y_S_S_S_C1001_bb+ABY*I_ERI_D2y_S_S_S_C1001_bb;
  Double I_ERI_Dyz_Py_S_S_C1001_bb = I_ERI_F2yz_S_S_S_C1001_bb+ABY*I_ERI_Dyz_S_S_S_C1001_bb;
  Double I_ERI_D2z_Py_S_S_C1001_bb = I_ERI_Fy2z_S_S_S_C1001_bb+ABY*I_ERI_D2z_S_S_S_C1001_bb;
  Double I_ERI_D2x_Pz_S_S_C1001_bb = I_ERI_F2xz_S_S_S_C1001_bb+ABZ*I_ERI_D2x_S_S_S_C1001_bb;
  Double I_ERI_Dxy_Pz_S_S_C1001_bb = I_ERI_Fxyz_S_S_S_C1001_bb+ABZ*I_ERI_Dxy_S_S_S_C1001_bb;
  Double I_ERI_Dxz_Pz_S_S_C1001_bb = I_ERI_Fx2z_S_S_S_C1001_bb+ABZ*I_ERI_Dxz_S_S_S_C1001_bb;
  Double I_ERI_D2y_Pz_S_S_C1001_bb = I_ERI_F2yz_S_S_S_C1001_bb+ABZ*I_ERI_D2y_S_S_S_C1001_bb;
  Double I_ERI_Dyz_Pz_S_S_C1001_bb = I_ERI_Fy2z_S_S_S_C1001_bb+ABZ*I_ERI_Dyz_S_S_S_C1001_bb;
  Double I_ERI_D2z_Pz_S_S_C1001_bb = I_ERI_F3z_S_S_S_C1001_bb+ABZ*I_ERI_D2z_S_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_Px_D2x_S_S_C1001_bb = I_ERI_D2x_Px_S_S_C1001_bb+ABX*I_ERI_Px_Px_S_S_C1001_bb;
  Double I_ERI_Py_D2x_S_S_C1001_bb = I_ERI_Dxy_Px_S_S_C1001_bb+ABX*I_ERI_Py_Px_S_S_C1001_bb;
  Double I_ERI_Pz_D2x_S_S_C1001_bb = I_ERI_Dxz_Px_S_S_C1001_bb+ABX*I_ERI_Pz_Px_S_S_C1001_bb;
  Double I_ERI_Px_Dxy_S_S_C1001_bb = I_ERI_Dxy_Px_S_S_C1001_bb+ABY*I_ERI_Px_Px_S_S_C1001_bb;
  Double I_ERI_Py_Dxy_S_S_C1001_bb = I_ERI_D2y_Px_S_S_C1001_bb+ABY*I_ERI_Py_Px_S_S_C1001_bb;
  Double I_ERI_Pz_Dxy_S_S_C1001_bb = I_ERI_Dyz_Px_S_S_C1001_bb+ABY*I_ERI_Pz_Px_S_S_C1001_bb;
  Double I_ERI_Px_D2y_S_S_C1001_bb = I_ERI_Dxy_Py_S_S_C1001_bb+ABY*I_ERI_Px_Py_S_S_C1001_bb;
  Double I_ERI_Py_D2y_S_S_C1001_bb = I_ERI_D2y_Py_S_S_C1001_bb+ABY*I_ERI_Py_Py_S_S_C1001_bb;
  Double I_ERI_Pz_D2y_S_S_C1001_bb = I_ERI_Dyz_Py_S_S_C1001_bb+ABY*I_ERI_Pz_Py_S_S_C1001_bb;
  Double I_ERI_Px_D2z_S_S_C1001_bb = I_ERI_Dxz_Pz_S_S_C1001_bb+ABZ*I_ERI_Px_Pz_S_S_C1001_bb;
  Double I_ERI_Py_D2z_S_S_C1001_bb = I_ERI_Dyz_Pz_S_S_C1001_bb+ABZ*I_ERI_Py_Pz_S_S_C1001_bb;
  Double I_ERI_Pz_D2z_S_S_C1001_bb = I_ERI_D2z_Pz_S_S_C1001_bb+ABZ*I_ERI_Pz_Pz_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 10 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1001_bb = I_ERI_G4x_S_S_S_C1001_bb+ABX*I_ERI_F3x_S_S_S_C1001_bb;
  Double I_ERI_F2xy_Px_S_S_C1001_bb = I_ERI_G3xy_S_S_S_C1001_bb+ABX*I_ERI_F2xy_S_S_S_C1001_bb;
  Double I_ERI_F2xz_Px_S_S_C1001_bb = I_ERI_G3xz_S_S_S_C1001_bb+ABX*I_ERI_F2xz_S_S_S_C1001_bb;
  Double I_ERI_Fx2y_Px_S_S_C1001_bb = I_ERI_G2x2y_S_S_S_C1001_bb+ABX*I_ERI_Fx2y_S_S_S_C1001_bb;
  Double I_ERI_Fxyz_Px_S_S_C1001_bb = I_ERI_G2xyz_S_S_S_C1001_bb+ABX*I_ERI_Fxyz_S_S_S_C1001_bb;
  Double I_ERI_Fx2z_Px_S_S_C1001_bb = I_ERI_G2x2z_S_S_S_C1001_bb+ABX*I_ERI_Fx2z_S_S_S_C1001_bb;
  Double I_ERI_F2yz_Px_S_S_C1001_bb = I_ERI_Gx2yz_S_S_S_C1001_bb+ABX*I_ERI_F2yz_S_S_S_C1001_bb;
  Double I_ERI_Fy2z_Px_S_S_C1001_bb = I_ERI_Gxy2z_S_S_S_C1001_bb+ABX*I_ERI_Fy2z_S_S_S_C1001_bb;
  Double I_ERI_F2xy_Py_S_S_C1001_bb = I_ERI_G2x2y_S_S_S_C1001_bb+ABY*I_ERI_F2xy_S_S_S_C1001_bb;
  Double I_ERI_Fx2y_Py_S_S_C1001_bb = I_ERI_Gx3y_S_S_S_C1001_bb+ABY*I_ERI_Fx2y_S_S_S_C1001_bb;
  Double I_ERI_Fxyz_Py_S_S_C1001_bb = I_ERI_Gx2yz_S_S_S_C1001_bb+ABY*I_ERI_Fxyz_S_S_S_C1001_bb;
  Double I_ERI_F3y_Py_S_S_C1001_bb = I_ERI_G4y_S_S_S_C1001_bb+ABY*I_ERI_F3y_S_S_S_C1001_bb;
  Double I_ERI_F2yz_Py_S_S_C1001_bb = I_ERI_G3yz_S_S_S_C1001_bb+ABY*I_ERI_F2yz_S_S_S_C1001_bb;
  Double I_ERI_Fy2z_Py_S_S_C1001_bb = I_ERI_G2y2z_S_S_S_C1001_bb+ABY*I_ERI_Fy2z_S_S_S_C1001_bb;
  Double I_ERI_F2xz_Pz_S_S_C1001_bb = I_ERI_G2x2z_S_S_S_C1001_bb+ABZ*I_ERI_F2xz_S_S_S_C1001_bb;
  Double I_ERI_Fxyz_Pz_S_S_C1001_bb = I_ERI_Gxy2z_S_S_S_C1001_bb+ABZ*I_ERI_Fxyz_S_S_S_C1001_bb;
  Double I_ERI_Fx2z_Pz_S_S_C1001_bb = I_ERI_Gx3z_S_S_S_C1001_bb+ABZ*I_ERI_Fx2z_S_S_S_C1001_bb;
  Double I_ERI_F2yz_Pz_S_S_C1001_bb = I_ERI_G2y2z_S_S_S_C1001_bb+ABZ*I_ERI_F2yz_S_S_S_C1001_bb;
  Double I_ERI_Fy2z_Pz_S_S_C1001_bb = I_ERI_Gy3z_S_S_S_C1001_bb+ABZ*I_ERI_Fy2z_S_S_S_C1001_bb;
  Double I_ERI_F3z_Pz_S_S_C1001_bb = I_ERI_G4z_S_S_S_C1001_bb+ABZ*I_ERI_F3z_S_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_C1001_bb = I_ERI_F3x_Px_S_S_C1001_bb+ABX*I_ERI_D2x_Px_S_S_C1001_bb;
  Double I_ERI_Dxy_D2x_S_S_C1001_bb = I_ERI_F2xy_Px_S_S_C1001_bb+ABX*I_ERI_Dxy_Px_S_S_C1001_bb;
  Double I_ERI_Dxz_D2x_S_S_C1001_bb = I_ERI_F2xz_Px_S_S_C1001_bb+ABX*I_ERI_Dxz_Px_S_S_C1001_bb;
  Double I_ERI_D2y_D2x_S_S_C1001_bb = I_ERI_Fx2y_Px_S_S_C1001_bb+ABX*I_ERI_D2y_Px_S_S_C1001_bb;
  Double I_ERI_Dyz_D2x_S_S_C1001_bb = I_ERI_Fxyz_Px_S_S_C1001_bb+ABX*I_ERI_Dyz_Px_S_S_C1001_bb;
  Double I_ERI_D2z_D2x_S_S_C1001_bb = I_ERI_Fx2z_Px_S_S_C1001_bb+ABX*I_ERI_D2z_Px_S_S_C1001_bb;
  Double I_ERI_Dxz_Dxy_S_S_C1001_bb = I_ERI_Fxyz_Px_S_S_C1001_bb+ABY*I_ERI_Dxz_Px_S_S_C1001_bb;
  Double I_ERI_Dyz_Dxy_S_S_C1001_bb = I_ERI_F2yz_Px_S_S_C1001_bb+ABY*I_ERI_Dyz_Px_S_S_C1001_bb;
  Double I_ERI_D2z_Dxy_S_S_C1001_bb = I_ERI_Fy2z_Px_S_S_C1001_bb+ABY*I_ERI_D2z_Px_S_S_C1001_bb;
  Double I_ERI_D2x_D2y_S_S_C1001_bb = I_ERI_F2xy_Py_S_S_C1001_bb+ABY*I_ERI_D2x_Py_S_S_C1001_bb;
  Double I_ERI_Dxy_D2y_S_S_C1001_bb = I_ERI_Fx2y_Py_S_S_C1001_bb+ABY*I_ERI_Dxy_Py_S_S_C1001_bb;
  Double I_ERI_Dxz_D2y_S_S_C1001_bb = I_ERI_Fxyz_Py_S_S_C1001_bb+ABY*I_ERI_Dxz_Py_S_S_C1001_bb;
  Double I_ERI_D2y_D2y_S_S_C1001_bb = I_ERI_F3y_Py_S_S_C1001_bb+ABY*I_ERI_D2y_Py_S_S_C1001_bb;
  Double I_ERI_Dyz_D2y_S_S_C1001_bb = I_ERI_F2yz_Py_S_S_C1001_bb+ABY*I_ERI_Dyz_Py_S_S_C1001_bb;
  Double I_ERI_D2z_D2y_S_S_C1001_bb = I_ERI_Fy2z_Py_S_S_C1001_bb+ABY*I_ERI_D2z_Py_S_S_C1001_bb;
  Double I_ERI_D2x_D2z_S_S_C1001_bb = I_ERI_F2xz_Pz_S_S_C1001_bb+ABZ*I_ERI_D2x_Pz_S_S_C1001_bb;
  Double I_ERI_Dxy_D2z_S_S_C1001_bb = I_ERI_Fxyz_Pz_S_S_C1001_bb+ABZ*I_ERI_Dxy_Pz_S_S_C1001_bb;
  Double I_ERI_Dxz_D2z_S_S_C1001_bb = I_ERI_Fx2z_Pz_S_S_C1001_bb+ABZ*I_ERI_Dxz_Pz_S_S_C1001_bb;
  Double I_ERI_D2y_D2z_S_S_C1001_bb = I_ERI_F2yz_Pz_S_S_C1001_bb+ABZ*I_ERI_D2y_Pz_S_S_C1001_bb;
  Double I_ERI_Dyz_D2z_S_S_C1001_bb = I_ERI_Fy2z_Pz_S_S_C1001_bb+ABZ*I_ERI_Dyz_Pz_S_S_C1001_bb;
  Double I_ERI_D2z_D2z_S_S_C1001_bb = I_ERI_F3z_Pz_S_S_C1001_bb+ABZ*I_ERI_D2z_Pz_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1001_bb
   ************************************************************/
  Double I_ERI_Px_F3x_S_S_C1001_bb = I_ERI_D2x_D2x_S_S_C1001_bb+ABX*I_ERI_Px_D2x_S_S_C1001_bb;
  Double I_ERI_Py_F3x_S_S_C1001_bb = I_ERI_Dxy_D2x_S_S_C1001_bb+ABX*I_ERI_Py_D2x_S_S_C1001_bb;
  Double I_ERI_Pz_F3x_S_S_C1001_bb = I_ERI_Dxz_D2x_S_S_C1001_bb+ABX*I_ERI_Pz_D2x_S_S_C1001_bb;
  Double I_ERI_Px_F2xy_S_S_C1001_bb = I_ERI_Dxy_D2x_S_S_C1001_bb+ABY*I_ERI_Px_D2x_S_S_C1001_bb;
  Double I_ERI_Py_F2xy_S_S_C1001_bb = I_ERI_D2y_D2x_S_S_C1001_bb+ABY*I_ERI_Py_D2x_S_S_C1001_bb;
  Double I_ERI_Pz_F2xy_S_S_C1001_bb = I_ERI_Dyz_D2x_S_S_C1001_bb+ABY*I_ERI_Pz_D2x_S_S_C1001_bb;
  Double I_ERI_Px_F2xz_S_S_C1001_bb = I_ERI_Dxz_D2x_S_S_C1001_bb+ABZ*I_ERI_Px_D2x_S_S_C1001_bb;
  Double I_ERI_Py_F2xz_S_S_C1001_bb = I_ERI_Dyz_D2x_S_S_C1001_bb+ABZ*I_ERI_Py_D2x_S_S_C1001_bb;
  Double I_ERI_Pz_F2xz_S_S_C1001_bb = I_ERI_D2z_D2x_S_S_C1001_bb+ABZ*I_ERI_Pz_D2x_S_S_C1001_bb;
  Double I_ERI_Px_Fx2y_S_S_C1001_bb = I_ERI_D2x_D2y_S_S_C1001_bb+ABX*I_ERI_Px_D2y_S_S_C1001_bb;
  Double I_ERI_Py_Fx2y_S_S_C1001_bb = I_ERI_Dxy_D2y_S_S_C1001_bb+ABX*I_ERI_Py_D2y_S_S_C1001_bb;
  Double I_ERI_Pz_Fx2y_S_S_C1001_bb = I_ERI_Dxz_D2y_S_S_C1001_bb+ABX*I_ERI_Pz_D2y_S_S_C1001_bb;
  Double I_ERI_Px_Fxyz_S_S_C1001_bb = I_ERI_Dxz_Dxy_S_S_C1001_bb+ABZ*I_ERI_Px_Dxy_S_S_C1001_bb;
  Double I_ERI_Py_Fxyz_S_S_C1001_bb = I_ERI_Dyz_Dxy_S_S_C1001_bb+ABZ*I_ERI_Py_Dxy_S_S_C1001_bb;
  Double I_ERI_Pz_Fxyz_S_S_C1001_bb = I_ERI_D2z_Dxy_S_S_C1001_bb+ABZ*I_ERI_Pz_Dxy_S_S_C1001_bb;
  Double I_ERI_Px_Fx2z_S_S_C1001_bb = I_ERI_D2x_D2z_S_S_C1001_bb+ABX*I_ERI_Px_D2z_S_S_C1001_bb;
  Double I_ERI_Py_Fx2z_S_S_C1001_bb = I_ERI_Dxy_D2z_S_S_C1001_bb+ABX*I_ERI_Py_D2z_S_S_C1001_bb;
  Double I_ERI_Pz_Fx2z_S_S_C1001_bb = I_ERI_Dxz_D2z_S_S_C1001_bb+ABX*I_ERI_Pz_D2z_S_S_C1001_bb;
  Double I_ERI_Px_F3y_S_S_C1001_bb = I_ERI_Dxy_D2y_S_S_C1001_bb+ABY*I_ERI_Px_D2y_S_S_C1001_bb;
  Double I_ERI_Py_F3y_S_S_C1001_bb = I_ERI_D2y_D2y_S_S_C1001_bb+ABY*I_ERI_Py_D2y_S_S_C1001_bb;
  Double I_ERI_Pz_F3y_S_S_C1001_bb = I_ERI_Dyz_D2y_S_S_C1001_bb+ABY*I_ERI_Pz_D2y_S_S_C1001_bb;
  Double I_ERI_Px_F2yz_S_S_C1001_bb = I_ERI_Dxz_D2y_S_S_C1001_bb+ABZ*I_ERI_Px_D2y_S_S_C1001_bb;
  Double I_ERI_Py_F2yz_S_S_C1001_bb = I_ERI_Dyz_D2y_S_S_C1001_bb+ABZ*I_ERI_Py_D2y_S_S_C1001_bb;
  Double I_ERI_Pz_F2yz_S_S_C1001_bb = I_ERI_D2z_D2y_S_S_C1001_bb+ABZ*I_ERI_Pz_D2y_S_S_C1001_bb;
  Double I_ERI_Px_Fy2z_S_S_C1001_bb = I_ERI_Dxy_D2z_S_S_C1001_bb+ABY*I_ERI_Px_D2z_S_S_C1001_bb;
  Double I_ERI_Py_Fy2z_S_S_C1001_bb = I_ERI_D2y_D2z_S_S_C1001_bb+ABY*I_ERI_Py_D2z_S_S_C1001_bb;
  Double I_ERI_Pz_Fy2z_S_S_C1001_bb = I_ERI_Dyz_D2z_S_S_C1001_bb+ABY*I_ERI_Pz_D2z_S_S_C1001_bb;
  Double I_ERI_Px_F3z_S_S_C1001_bb = I_ERI_Dxz_D2z_S_S_C1001_bb+ABZ*I_ERI_Px_D2z_S_S_C1001_bb;
  Double I_ERI_Py_F3z_S_S_C1001_bb = I_ERI_Dyz_D2z_S_S_C1001_bb+ABZ*I_ERI_Py_D2z_S_S_C1001_bb;
  Double I_ERI_Pz_F3z_S_S_C1001_bb = I_ERI_D2z_D2z_S_S_C1001_bb+ABZ*I_ERI_Pz_D2z_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1_bc
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1_bc = I_ERI_D2x_S_Px_S_C1_bc+ABX*I_ERI_Px_S_Px_S_C1_bc;
  Double I_ERI_Py_Px_Px_S_C1_bc = I_ERI_Dxy_S_Px_S_C1_bc+ABX*I_ERI_Py_S_Px_S_C1_bc;
  Double I_ERI_Pz_Px_Px_S_C1_bc = I_ERI_Dxz_S_Px_S_C1_bc+ABX*I_ERI_Pz_S_Px_S_C1_bc;
  Double I_ERI_Px_Py_Px_S_C1_bc = I_ERI_Dxy_S_Px_S_C1_bc+ABY*I_ERI_Px_S_Px_S_C1_bc;
  Double I_ERI_Py_Py_Px_S_C1_bc = I_ERI_D2y_S_Px_S_C1_bc+ABY*I_ERI_Py_S_Px_S_C1_bc;
  Double I_ERI_Pz_Py_Px_S_C1_bc = I_ERI_Dyz_S_Px_S_C1_bc+ABY*I_ERI_Pz_S_Px_S_C1_bc;
  Double I_ERI_Px_Pz_Px_S_C1_bc = I_ERI_Dxz_S_Px_S_C1_bc+ABZ*I_ERI_Px_S_Px_S_C1_bc;
  Double I_ERI_Py_Pz_Px_S_C1_bc = I_ERI_Dyz_S_Px_S_C1_bc+ABZ*I_ERI_Py_S_Px_S_C1_bc;
  Double I_ERI_Pz_Pz_Px_S_C1_bc = I_ERI_D2z_S_Px_S_C1_bc+ABZ*I_ERI_Pz_S_Px_S_C1_bc;
  Double I_ERI_Px_Px_Py_S_C1_bc = I_ERI_D2x_S_Py_S_C1_bc+ABX*I_ERI_Px_S_Py_S_C1_bc;
  Double I_ERI_Py_Px_Py_S_C1_bc = I_ERI_Dxy_S_Py_S_C1_bc+ABX*I_ERI_Py_S_Py_S_C1_bc;
  Double I_ERI_Pz_Px_Py_S_C1_bc = I_ERI_Dxz_S_Py_S_C1_bc+ABX*I_ERI_Pz_S_Py_S_C1_bc;
  Double I_ERI_Px_Py_Py_S_C1_bc = I_ERI_Dxy_S_Py_S_C1_bc+ABY*I_ERI_Px_S_Py_S_C1_bc;
  Double I_ERI_Py_Py_Py_S_C1_bc = I_ERI_D2y_S_Py_S_C1_bc+ABY*I_ERI_Py_S_Py_S_C1_bc;
  Double I_ERI_Pz_Py_Py_S_C1_bc = I_ERI_Dyz_S_Py_S_C1_bc+ABY*I_ERI_Pz_S_Py_S_C1_bc;
  Double I_ERI_Px_Pz_Py_S_C1_bc = I_ERI_Dxz_S_Py_S_C1_bc+ABZ*I_ERI_Px_S_Py_S_C1_bc;
  Double I_ERI_Py_Pz_Py_S_C1_bc = I_ERI_Dyz_S_Py_S_C1_bc+ABZ*I_ERI_Py_S_Py_S_C1_bc;
  Double I_ERI_Pz_Pz_Py_S_C1_bc = I_ERI_D2z_S_Py_S_C1_bc+ABZ*I_ERI_Pz_S_Py_S_C1_bc;
  Double I_ERI_Px_Px_Pz_S_C1_bc = I_ERI_D2x_S_Pz_S_C1_bc+ABX*I_ERI_Px_S_Pz_S_C1_bc;
  Double I_ERI_Py_Px_Pz_S_C1_bc = I_ERI_Dxy_S_Pz_S_C1_bc+ABX*I_ERI_Py_S_Pz_S_C1_bc;
  Double I_ERI_Pz_Px_Pz_S_C1_bc = I_ERI_Dxz_S_Pz_S_C1_bc+ABX*I_ERI_Pz_S_Pz_S_C1_bc;
  Double I_ERI_Px_Py_Pz_S_C1_bc = I_ERI_Dxy_S_Pz_S_C1_bc+ABY*I_ERI_Px_S_Pz_S_C1_bc;
  Double I_ERI_Py_Py_Pz_S_C1_bc = I_ERI_D2y_S_Pz_S_C1_bc+ABY*I_ERI_Py_S_Pz_S_C1_bc;
  Double I_ERI_Pz_Py_Pz_S_C1_bc = I_ERI_Dyz_S_Pz_S_C1_bc+ABY*I_ERI_Pz_S_Pz_S_C1_bc;
  Double I_ERI_Px_Pz_Pz_S_C1_bc = I_ERI_Dxz_S_Pz_S_C1_bc+ABZ*I_ERI_Px_S_Pz_S_C1_bc;
  Double I_ERI_Py_Pz_Pz_S_C1_bc = I_ERI_Dyz_S_Pz_S_C1_bc+ABZ*I_ERI_Py_S_Pz_S_C1_bc;
  Double I_ERI_Pz_Pz_Pz_S_C1_bc = I_ERI_D2z_S_Pz_S_C1_bc+ABZ*I_ERI_Pz_S_Pz_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_bc
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001_bc = I_ERI_D2x_S_Px_S_C1001_bc+ABX*I_ERI_Px_S_Px_S_C1001_bc;
  Double I_ERI_Py_Px_Px_S_C1001_bc = I_ERI_Dxy_S_Px_S_C1001_bc+ABX*I_ERI_Py_S_Px_S_C1001_bc;
  Double I_ERI_Pz_Px_Px_S_C1001_bc = I_ERI_Dxz_S_Px_S_C1001_bc+ABX*I_ERI_Pz_S_Px_S_C1001_bc;
  Double I_ERI_Px_Py_Px_S_C1001_bc = I_ERI_Dxy_S_Px_S_C1001_bc+ABY*I_ERI_Px_S_Px_S_C1001_bc;
  Double I_ERI_Py_Py_Px_S_C1001_bc = I_ERI_D2y_S_Px_S_C1001_bc+ABY*I_ERI_Py_S_Px_S_C1001_bc;
  Double I_ERI_Pz_Py_Px_S_C1001_bc = I_ERI_Dyz_S_Px_S_C1001_bc+ABY*I_ERI_Pz_S_Px_S_C1001_bc;
  Double I_ERI_Px_Pz_Px_S_C1001_bc = I_ERI_Dxz_S_Px_S_C1001_bc+ABZ*I_ERI_Px_S_Px_S_C1001_bc;
  Double I_ERI_Py_Pz_Px_S_C1001_bc = I_ERI_Dyz_S_Px_S_C1001_bc+ABZ*I_ERI_Py_S_Px_S_C1001_bc;
  Double I_ERI_Pz_Pz_Px_S_C1001_bc = I_ERI_D2z_S_Px_S_C1001_bc+ABZ*I_ERI_Pz_S_Px_S_C1001_bc;
  Double I_ERI_Px_Px_Py_S_C1001_bc = I_ERI_D2x_S_Py_S_C1001_bc+ABX*I_ERI_Px_S_Py_S_C1001_bc;
  Double I_ERI_Py_Px_Py_S_C1001_bc = I_ERI_Dxy_S_Py_S_C1001_bc+ABX*I_ERI_Py_S_Py_S_C1001_bc;
  Double I_ERI_Pz_Px_Py_S_C1001_bc = I_ERI_Dxz_S_Py_S_C1001_bc+ABX*I_ERI_Pz_S_Py_S_C1001_bc;
  Double I_ERI_Px_Py_Py_S_C1001_bc = I_ERI_Dxy_S_Py_S_C1001_bc+ABY*I_ERI_Px_S_Py_S_C1001_bc;
  Double I_ERI_Py_Py_Py_S_C1001_bc = I_ERI_D2y_S_Py_S_C1001_bc+ABY*I_ERI_Py_S_Py_S_C1001_bc;
  Double I_ERI_Pz_Py_Py_S_C1001_bc = I_ERI_Dyz_S_Py_S_C1001_bc+ABY*I_ERI_Pz_S_Py_S_C1001_bc;
  Double I_ERI_Px_Pz_Py_S_C1001_bc = I_ERI_Dxz_S_Py_S_C1001_bc+ABZ*I_ERI_Px_S_Py_S_C1001_bc;
  Double I_ERI_Py_Pz_Py_S_C1001_bc = I_ERI_Dyz_S_Py_S_C1001_bc+ABZ*I_ERI_Py_S_Py_S_C1001_bc;
  Double I_ERI_Pz_Pz_Py_S_C1001_bc = I_ERI_D2z_S_Py_S_C1001_bc+ABZ*I_ERI_Pz_S_Py_S_C1001_bc;
  Double I_ERI_Px_Px_Pz_S_C1001_bc = I_ERI_D2x_S_Pz_S_C1001_bc+ABX*I_ERI_Px_S_Pz_S_C1001_bc;
  Double I_ERI_Py_Px_Pz_S_C1001_bc = I_ERI_Dxy_S_Pz_S_C1001_bc+ABX*I_ERI_Py_S_Pz_S_C1001_bc;
  Double I_ERI_Pz_Px_Pz_S_C1001_bc = I_ERI_Dxz_S_Pz_S_C1001_bc+ABX*I_ERI_Pz_S_Pz_S_C1001_bc;
  Double I_ERI_Px_Py_Pz_S_C1001_bc = I_ERI_Dxy_S_Pz_S_C1001_bc+ABY*I_ERI_Px_S_Pz_S_C1001_bc;
  Double I_ERI_Py_Py_Pz_S_C1001_bc = I_ERI_D2y_S_Pz_S_C1001_bc+ABY*I_ERI_Py_S_Pz_S_C1001_bc;
  Double I_ERI_Pz_Py_Pz_S_C1001_bc = I_ERI_Dyz_S_Pz_S_C1001_bc+ABY*I_ERI_Pz_S_Pz_S_C1001_bc;
  Double I_ERI_Px_Pz_Pz_S_C1001_bc = I_ERI_Dxz_S_Pz_S_C1001_bc+ABZ*I_ERI_Px_S_Pz_S_C1001_bc;
  Double I_ERI_Py_Pz_Pz_S_C1001_bc = I_ERI_Dyz_S_Pz_S_C1001_bc+ABZ*I_ERI_Py_S_Pz_S_C1001_bc;
  Double I_ERI_Pz_Pz_Pz_S_C1001_bc = I_ERI_D2z_S_Pz_S_C1001_bc+ABZ*I_ERI_Pz_S_Pz_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001_bc
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1001_bc = I_ERI_F3x_S_Px_S_C1001_bc+ABX*I_ERI_D2x_S_Px_S_C1001_bc;
  Double I_ERI_Dxy_Px_Px_S_C1001_bc = I_ERI_F2xy_S_Px_S_C1001_bc+ABX*I_ERI_Dxy_S_Px_S_C1001_bc;
  Double I_ERI_Dxz_Px_Px_S_C1001_bc = I_ERI_F2xz_S_Px_S_C1001_bc+ABX*I_ERI_Dxz_S_Px_S_C1001_bc;
  Double I_ERI_D2y_Px_Px_S_C1001_bc = I_ERI_Fx2y_S_Px_S_C1001_bc+ABX*I_ERI_D2y_S_Px_S_C1001_bc;
  Double I_ERI_Dyz_Px_Px_S_C1001_bc = I_ERI_Fxyz_S_Px_S_C1001_bc+ABX*I_ERI_Dyz_S_Px_S_C1001_bc;
  Double I_ERI_D2z_Px_Px_S_C1001_bc = I_ERI_Fx2z_S_Px_S_C1001_bc+ABX*I_ERI_D2z_S_Px_S_C1001_bc;
  Double I_ERI_Dxy_Py_Px_S_C1001_bc = I_ERI_Fx2y_S_Px_S_C1001_bc+ABY*I_ERI_Dxy_S_Px_S_C1001_bc;
  Double I_ERI_Dxz_Py_Px_S_C1001_bc = I_ERI_Fxyz_S_Px_S_C1001_bc+ABY*I_ERI_Dxz_S_Px_S_C1001_bc;
  Double I_ERI_D2y_Py_Px_S_C1001_bc = I_ERI_F3y_S_Px_S_C1001_bc+ABY*I_ERI_D2y_S_Px_S_C1001_bc;
  Double I_ERI_Dyz_Py_Px_S_C1001_bc = I_ERI_F2yz_S_Px_S_C1001_bc+ABY*I_ERI_Dyz_S_Px_S_C1001_bc;
  Double I_ERI_D2z_Py_Px_S_C1001_bc = I_ERI_Fy2z_S_Px_S_C1001_bc+ABY*I_ERI_D2z_S_Px_S_C1001_bc;
  Double I_ERI_Dxz_Pz_Px_S_C1001_bc = I_ERI_Fx2z_S_Px_S_C1001_bc+ABZ*I_ERI_Dxz_S_Px_S_C1001_bc;
  Double I_ERI_Dyz_Pz_Px_S_C1001_bc = I_ERI_Fy2z_S_Px_S_C1001_bc+ABZ*I_ERI_Dyz_S_Px_S_C1001_bc;
  Double I_ERI_D2z_Pz_Px_S_C1001_bc = I_ERI_F3z_S_Px_S_C1001_bc+ABZ*I_ERI_D2z_S_Px_S_C1001_bc;
  Double I_ERI_D2x_Px_Py_S_C1001_bc = I_ERI_F3x_S_Py_S_C1001_bc+ABX*I_ERI_D2x_S_Py_S_C1001_bc;
  Double I_ERI_Dxy_Px_Py_S_C1001_bc = I_ERI_F2xy_S_Py_S_C1001_bc+ABX*I_ERI_Dxy_S_Py_S_C1001_bc;
  Double I_ERI_Dxz_Px_Py_S_C1001_bc = I_ERI_F2xz_S_Py_S_C1001_bc+ABX*I_ERI_Dxz_S_Py_S_C1001_bc;
  Double I_ERI_D2y_Px_Py_S_C1001_bc = I_ERI_Fx2y_S_Py_S_C1001_bc+ABX*I_ERI_D2y_S_Py_S_C1001_bc;
  Double I_ERI_Dyz_Px_Py_S_C1001_bc = I_ERI_Fxyz_S_Py_S_C1001_bc+ABX*I_ERI_Dyz_S_Py_S_C1001_bc;
  Double I_ERI_D2z_Px_Py_S_C1001_bc = I_ERI_Fx2z_S_Py_S_C1001_bc+ABX*I_ERI_D2z_S_Py_S_C1001_bc;
  Double I_ERI_Dxy_Py_Py_S_C1001_bc = I_ERI_Fx2y_S_Py_S_C1001_bc+ABY*I_ERI_Dxy_S_Py_S_C1001_bc;
  Double I_ERI_Dxz_Py_Py_S_C1001_bc = I_ERI_Fxyz_S_Py_S_C1001_bc+ABY*I_ERI_Dxz_S_Py_S_C1001_bc;
  Double I_ERI_D2y_Py_Py_S_C1001_bc = I_ERI_F3y_S_Py_S_C1001_bc+ABY*I_ERI_D2y_S_Py_S_C1001_bc;
  Double I_ERI_Dyz_Py_Py_S_C1001_bc = I_ERI_F2yz_S_Py_S_C1001_bc+ABY*I_ERI_Dyz_S_Py_S_C1001_bc;
  Double I_ERI_D2z_Py_Py_S_C1001_bc = I_ERI_Fy2z_S_Py_S_C1001_bc+ABY*I_ERI_D2z_S_Py_S_C1001_bc;
  Double I_ERI_Dxz_Pz_Py_S_C1001_bc = I_ERI_Fx2z_S_Py_S_C1001_bc+ABZ*I_ERI_Dxz_S_Py_S_C1001_bc;
  Double I_ERI_Dyz_Pz_Py_S_C1001_bc = I_ERI_Fy2z_S_Py_S_C1001_bc+ABZ*I_ERI_Dyz_S_Py_S_C1001_bc;
  Double I_ERI_D2z_Pz_Py_S_C1001_bc = I_ERI_F3z_S_Py_S_C1001_bc+ABZ*I_ERI_D2z_S_Py_S_C1001_bc;
  Double I_ERI_D2x_Px_Pz_S_C1001_bc = I_ERI_F3x_S_Pz_S_C1001_bc+ABX*I_ERI_D2x_S_Pz_S_C1001_bc;
  Double I_ERI_Dxy_Px_Pz_S_C1001_bc = I_ERI_F2xy_S_Pz_S_C1001_bc+ABX*I_ERI_Dxy_S_Pz_S_C1001_bc;
  Double I_ERI_Dxz_Px_Pz_S_C1001_bc = I_ERI_F2xz_S_Pz_S_C1001_bc+ABX*I_ERI_Dxz_S_Pz_S_C1001_bc;
  Double I_ERI_D2y_Px_Pz_S_C1001_bc = I_ERI_Fx2y_S_Pz_S_C1001_bc+ABX*I_ERI_D2y_S_Pz_S_C1001_bc;
  Double I_ERI_Dyz_Px_Pz_S_C1001_bc = I_ERI_Fxyz_S_Pz_S_C1001_bc+ABX*I_ERI_Dyz_S_Pz_S_C1001_bc;
  Double I_ERI_D2z_Px_Pz_S_C1001_bc = I_ERI_Fx2z_S_Pz_S_C1001_bc+ABX*I_ERI_D2z_S_Pz_S_C1001_bc;
  Double I_ERI_Dxy_Py_Pz_S_C1001_bc = I_ERI_Fx2y_S_Pz_S_C1001_bc+ABY*I_ERI_Dxy_S_Pz_S_C1001_bc;
  Double I_ERI_Dxz_Py_Pz_S_C1001_bc = I_ERI_Fxyz_S_Pz_S_C1001_bc+ABY*I_ERI_Dxz_S_Pz_S_C1001_bc;
  Double I_ERI_D2y_Py_Pz_S_C1001_bc = I_ERI_F3y_S_Pz_S_C1001_bc+ABY*I_ERI_D2y_S_Pz_S_C1001_bc;
  Double I_ERI_Dyz_Py_Pz_S_C1001_bc = I_ERI_F2yz_S_Pz_S_C1001_bc+ABY*I_ERI_Dyz_S_Pz_S_C1001_bc;
  Double I_ERI_D2z_Py_Pz_S_C1001_bc = I_ERI_Fy2z_S_Pz_S_C1001_bc+ABY*I_ERI_D2z_S_Pz_S_C1001_bc;
  Double I_ERI_Dxz_Pz_Pz_S_C1001_bc = I_ERI_Fx2z_S_Pz_S_C1001_bc+ABZ*I_ERI_Dxz_S_Pz_S_C1001_bc;
  Double I_ERI_Dyz_Pz_Pz_S_C1001_bc = I_ERI_Fy2z_S_Pz_S_C1001_bc+ABZ*I_ERI_Dyz_S_Pz_S_C1001_bc;
  Double I_ERI_D2z_Pz_Pz_S_C1001_bc = I_ERI_F3z_S_Pz_S_C1001_bc+ABZ*I_ERI_D2z_S_Pz_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001_bc
   ************************************************************/
  Double I_ERI_Px_D2x_Px_S_C1001_bc = I_ERI_D2x_Px_Px_S_C1001_bc+ABX*I_ERI_Px_Px_Px_S_C1001_bc;
  Double I_ERI_Py_D2x_Px_S_C1001_bc = I_ERI_Dxy_Px_Px_S_C1001_bc+ABX*I_ERI_Py_Px_Px_S_C1001_bc;
  Double I_ERI_Pz_D2x_Px_S_C1001_bc = I_ERI_Dxz_Px_Px_S_C1001_bc+ABX*I_ERI_Pz_Px_Px_S_C1001_bc;
  Double I_ERI_Px_Dxy_Px_S_C1001_bc = I_ERI_Dxy_Px_Px_S_C1001_bc+ABY*I_ERI_Px_Px_Px_S_C1001_bc;
  Double I_ERI_Py_Dxy_Px_S_C1001_bc = I_ERI_D2y_Px_Px_S_C1001_bc+ABY*I_ERI_Py_Px_Px_S_C1001_bc;
  Double I_ERI_Pz_Dxy_Px_S_C1001_bc = I_ERI_Dyz_Px_Px_S_C1001_bc+ABY*I_ERI_Pz_Px_Px_S_C1001_bc;
  Double I_ERI_Px_Dxz_Px_S_C1001_bc = I_ERI_Dxz_Px_Px_S_C1001_bc+ABZ*I_ERI_Px_Px_Px_S_C1001_bc;
  Double I_ERI_Py_Dxz_Px_S_C1001_bc = I_ERI_Dyz_Px_Px_S_C1001_bc+ABZ*I_ERI_Py_Px_Px_S_C1001_bc;
  Double I_ERI_Pz_Dxz_Px_S_C1001_bc = I_ERI_D2z_Px_Px_S_C1001_bc+ABZ*I_ERI_Pz_Px_Px_S_C1001_bc;
  Double I_ERI_Px_D2y_Px_S_C1001_bc = I_ERI_Dxy_Py_Px_S_C1001_bc+ABY*I_ERI_Px_Py_Px_S_C1001_bc;
  Double I_ERI_Py_D2y_Px_S_C1001_bc = I_ERI_D2y_Py_Px_S_C1001_bc+ABY*I_ERI_Py_Py_Px_S_C1001_bc;
  Double I_ERI_Pz_D2y_Px_S_C1001_bc = I_ERI_Dyz_Py_Px_S_C1001_bc+ABY*I_ERI_Pz_Py_Px_S_C1001_bc;
  Double I_ERI_Px_Dyz_Px_S_C1001_bc = I_ERI_Dxz_Py_Px_S_C1001_bc+ABZ*I_ERI_Px_Py_Px_S_C1001_bc;
  Double I_ERI_Py_Dyz_Px_S_C1001_bc = I_ERI_Dyz_Py_Px_S_C1001_bc+ABZ*I_ERI_Py_Py_Px_S_C1001_bc;
  Double I_ERI_Pz_Dyz_Px_S_C1001_bc = I_ERI_D2z_Py_Px_S_C1001_bc+ABZ*I_ERI_Pz_Py_Px_S_C1001_bc;
  Double I_ERI_Px_D2z_Px_S_C1001_bc = I_ERI_Dxz_Pz_Px_S_C1001_bc+ABZ*I_ERI_Px_Pz_Px_S_C1001_bc;
  Double I_ERI_Py_D2z_Px_S_C1001_bc = I_ERI_Dyz_Pz_Px_S_C1001_bc+ABZ*I_ERI_Py_Pz_Px_S_C1001_bc;
  Double I_ERI_Pz_D2z_Px_S_C1001_bc = I_ERI_D2z_Pz_Px_S_C1001_bc+ABZ*I_ERI_Pz_Pz_Px_S_C1001_bc;
  Double I_ERI_Px_D2x_Py_S_C1001_bc = I_ERI_D2x_Px_Py_S_C1001_bc+ABX*I_ERI_Px_Px_Py_S_C1001_bc;
  Double I_ERI_Py_D2x_Py_S_C1001_bc = I_ERI_Dxy_Px_Py_S_C1001_bc+ABX*I_ERI_Py_Px_Py_S_C1001_bc;
  Double I_ERI_Pz_D2x_Py_S_C1001_bc = I_ERI_Dxz_Px_Py_S_C1001_bc+ABX*I_ERI_Pz_Px_Py_S_C1001_bc;
  Double I_ERI_Px_Dxy_Py_S_C1001_bc = I_ERI_Dxy_Px_Py_S_C1001_bc+ABY*I_ERI_Px_Px_Py_S_C1001_bc;
  Double I_ERI_Py_Dxy_Py_S_C1001_bc = I_ERI_D2y_Px_Py_S_C1001_bc+ABY*I_ERI_Py_Px_Py_S_C1001_bc;
  Double I_ERI_Pz_Dxy_Py_S_C1001_bc = I_ERI_Dyz_Px_Py_S_C1001_bc+ABY*I_ERI_Pz_Px_Py_S_C1001_bc;
  Double I_ERI_Px_Dxz_Py_S_C1001_bc = I_ERI_Dxz_Px_Py_S_C1001_bc+ABZ*I_ERI_Px_Px_Py_S_C1001_bc;
  Double I_ERI_Py_Dxz_Py_S_C1001_bc = I_ERI_Dyz_Px_Py_S_C1001_bc+ABZ*I_ERI_Py_Px_Py_S_C1001_bc;
  Double I_ERI_Pz_Dxz_Py_S_C1001_bc = I_ERI_D2z_Px_Py_S_C1001_bc+ABZ*I_ERI_Pz_Px_Py_S_C1001_bc;
  Double I_ERI_Px_D2y_Py_S_C1001_bc = I_ERI_Dxy_Py_Py_S_C1001_bc+ABY*I_ERI_Px_Py_Py_S_C1001_bc;
  Double I_ERI_Py_D2y_Py_S_C1001_bc = I_ERI_D2y_Py_Py_S_C1001_bc+ABY*I_ERI_Py_Py_Py_S_C1001_bc;
  Double I_ERI_Pz_D2y_Py_S_C1001_bc = I_ERI_Dyz_Py_Py_S_C1001_bc+ABY*I_ERI_Pz_Py_Py_S_C1001_bc;
  Double I_ERI_Px_Dyz_Py_S_C1001_bc = I_ERI_Dxz_Py_Py_S_C1001_bc+ABZ*I_ERI_Px_Py_Py_S_C1001_bc;
  Double I_ERI_Py_Dyz_Py_S_C1001_bc = I_ERI_Dyz_Py_Py_S_C1001_bc+ABZ*I_ERI_Py_Py_Py_S_C1001_bc;
  Double I_ERI_Pz_Dyz_Py_S_C1001_bc = I_ERI_D2z_Py_Py_S_C1001_bc+ABZ*I_ERI_Pz_Py_Py_S_C1001_bc;
  Double I_ERI_Px_D2z_Py_S_C1001_bc = I_ERI_Dxz_Pz_Py_S_C1001_bc+ABZ*I_ERI_Px_Pz_Py_S_C1001_bc;
  Double I_ERI_Py_D2z_Py_S_C1001_bc = I_ERI_Dyz_Pz_Py_S_C1001_bc+ABZ*I_ERI_Py_Pz_Py_S_C1001_bc;
  Double I_ERI_Pz_D2z_Py_S_C1001_bc = I_ERI_D2z_Pz_Py_S_C1001_bc+ABZ*I_ERI_Pz_Pz_Py_S_C1001_bc;
  Double I_ERI_Px_D2x_Pz_S_C1001_bc = I_ERI_D2x_Px_Pz_S_C1001_bc+ABX*I_ERI_Px_Px_Pz_S_C1001_bc;
  Double I_ERI_Py_D2x_Pz_S_C1001_bc = I_ERI_Dxy_Px_Pz_S_C1001_bc+ABX*I_ERI_Py_Px_Pz_S_C1001_bc;
  Double I_ERI_Pz_D2x_Pz_S_C1001_bc = I_ERI_Dxz_Px_Pz_S_C1001_bc+ABX*I_ERI_Pz_Px_Pz_S_C1001_bc;
  Double I_ERI_Px_Dxy_Pz_S_C1001_bc = I_ERI_Dxy_Px_Pz_S_C1001_bc+ABY*I_ERI_Px_Px_Pz_S_C1001_bc;
  Double I_ERI_Py_Dxy_Pz_S_C1001_bc = I_ERI_D2y_Px_Pz_S_C1001_bc+ABY*I_ERI_Py_Px_Pz_S_C1001_bc;
  Double I_ERI_Pz_Dxy_Pz_S_C1001_bc = I_ERI_Dyz_Px_Pz_S_C1001_bc+ABY*I_ERI_Pz_Px_Pz_S_C1001_bc;
  Double I_ERI_Px_Dxz_Pz_S_C1001_bc = I_ERI_Dxz_Px_Pz_S_C1001_bc+ABZ*I_ERI_Px_Px_Pz_S_C1001_bc;
  Double I_ERI_Py_Dxz_Pz_S_C1001_bc = I_ERI_Dyz_Px_Pz_S_C1001_bc+ABZ*I_ERI_Py_Px_Pz_S_C1001_bc;
  Double I_ERI_Pz_Dxz_Pz_S_C1001_bc = I_ERI_D2z_Px_Pz_S_C1001_bc+ABZ*I_ERI_Pz_Px_Pz_S_C1001_bc;
  Double I_ERI_Px_D2y_Pz_S_C1001_bc = I_ERI_Dxy_Py_Pz_S_C1001_bc+ABY*I_ERI_Px_Py_Pz_S_C1001_bc;
  Double I_ERI_Py_D2y_Pz_S_C1001_bc = I_ERI_D2y_Py_Pz_S_C1001_bc+ABY*I_ERI_Py_Py_Pz_S_C1001_bc;
  Double I_ERI_Pz_D2y_Pz_S_C1001_bc = I_ERI_Dyz_Py_Pz_S_C1001_bc+ABY*I_ERI_Pz_Py_Pz_S_C1001_bc;
  Double I_ERI_Px_Dyz_Pz_S_C1001_bc = I_ERI_Dxz_Py_Pz_S_C1001_bc+ABZ*I_ERI_Px_Py_Pz_S_C1001_bc;
  Double I_ERI_Py_Dyz_Pz_S_C1001_bc = I_ERI_Dyz_Py_Pz_S_C1001_bc+ABZ*I_ERI_Py_Py_Pz_S_C1001_bc;
  Double I_ERI_Pz_Dyz_Pz_S_C1001_bc = I_ERI_D2z_Py_Pz_S_C1001_bc+ABZ*I_ERI_Pz_Py_Pz_S_C1001_bc;
  Double I_ERI_Px_D2z_Pz_S_C1001_bc = I_ERI_Dxz_Pz_Pz_S_C1001_bc+ABZ*I_ERI_Px_Pz_Pz_S_C1001_bc;
  Double I_ERI_Py_D2z_Pz_S_C1001_bc = I_ERI_Dyz_Pz_Pz_S_C1001_bc+ABZ*I_ERI_Py_Pz_Pz_S_C1001_bc;
  Double I_ERI_Pz_D2z_Pz_S_C1001_bc = I_ERI_D2z_Pz_Pz_S_C1001_bc+ABZ*I_ERI_Pz_Pz_Pz_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001_cc
   ************************************************************/
  Double I_ERI_Px_Px_D2x_S_C1001_cc = I_ERI_D2x_S_D2x_S_C1001_cc+ABX*I_ERI_Px_S_D2x_S_C1001_cc;
  Double I_ERI_Py_Px_D2x_S_C1001_cc = I_ERI_Dxy_S_D2x_S_C1001_cc+ABX*I_ERI_Py_S_D2x_S_C1001_cc;
  Double I_ERI_Pz_Px_D2x_S_C1001_cc = I_ERI_Dxz_S_D2x_S_C1001_cc+ABX*I_ERI_Pz_S_D2x_S_C1001_cc;
  Double I_ERI_Px_Py_D2x_S_C1001_cc = I_ERI_Dxy_S_D2x_S_C1001_cc+ABY*I_ERI_Px_S_D2x_S_C1001_cc;
  Double I_ERI_Py_Py_D2x_S_C1001_cc = I_ERI_D2y_S_D2x_S_C1001_cc+ABY*I_ERI_Py_S_D2x_S_C1001_cc;
  Double I_ERI_Pz_Py_D2x_S_C1001_cc = I_ERI_Dyz_S_D2x_S_C1001_cc+ABY*I_ERI_Pz_S_D2x_S_C1001_cc;
  Double I_ERI_Px_Pz_D2x_S_C1001_cc = I_ERI_Dxz_S_D2x_S_C1001_cc+ABZ*I_ERI_Px_S_D2x_S_C1001_cc;
  Double I_ERI_Py_Pz_D2x_S_C1001_cc = I_ERI_Dyz_S_D2x_S_C1001_cc+ABZ*I_ERI_Py_S_D2x_S_C1001_cc;
  Double I_ERI_Pz_Pz_D2x_S_C1001_cc = I_ERI_D2z_S_D2x_S_C1001_cc+ABZ*I_ERI_Pz_S_D2x_S_C1001_cc;
  Double I_ERI_Px_Px_Dxy_S_C1001_cc = I_ERI_D2x_S_Dxy_S_C1001_cc+ABX*I_ERI_Px_S_Dxy_S_C1001_cc;
  Double I_ERI_Py_Px_Dxy_S_C1001_cc = I_ERI_Dxy_S_Dxy_S_C1001_cc+ABX*I_ERI_Py_S_Dxy_S_C1001_cc;
  Double I_ERI_Pz_Px_Dxy_S_C1001_cc = I_ERI_Dxz_S_Dxy_S_C1001_cc+ABX*I_ERI_Pz_S_Dxy_S_C1001_cc;
  Double I_ERI_Px_Py_Dxy_S_C1001_cc = I_ERI_Dxy_S_Dxy_S_C1001_cc+ABY*I_ERI_Px_S_Dxy_S_C1001_cc;
  Double I_ERI_Py_Py_Dxy_S_C1001_cc = I_ERI_D2y_S_Dxy_S_C1001_cc+ABY*I_ERI_Py_S_Dxy_S_C1001_cc;
  Double I_ERI_Pz_Py_Dxy_S_C1001_cc = I_ERI_Dyz_S_Dxy_S_C1001_cc+ABY*I_ERI_Pz_S_Dxy_S_C1001_cc;
  Double I_ERI_Px_Pz_Dxy_S_C1001_cc = I_ERI_Dxz_S_Dxy_S_C1001_cc+ABZ*I_ERI_Px_S_Dxy_S_C1001_cc;
  Double I_ERI_Py_Pz_Dxy_S_C1001_cc = I_ERI_Dyz_S_Dxy_S_C1001_cc+ABZ*I_ERI_Py_S_Dxy_S_C1001_cc;
  Double I_ERI_Pz_Pz_Dxy_S_C1001_cc = I_ERI_D2z_S_Dxy_S_C1001_cc+ABZ*I_ERI_Pz_S_Dxy_S_C1001_cc;
  Double I_ERI_Px_Px_Dxz_S_C1001_cc = I_ERI_D2x_S_Dxz_S_C1001_cc+ABX*I_ERI_Px_S_Dxz_S_C1001_cc;
  Double I_ERI_Py_Px_Dxz_S_C1001_cc = I_ERI_Dxy_S_Dxz_S_C1001_cc+ABX*I_ERI_Py_S_Dxz_S_C1001_cc;
  Double I_ERI_Pz_Px_Dxz_S_C1001_cc = I_ERI_Dxz_S_Dxz_S_C1001_cc+ABX*I_ERI_Pz_S_Dxz_S_C1001_cc;
  Double I_ERI_Px_Py_Dxz_S_C1001_cc = I_ERI_Dxy_S_Dxz_S_C1001_cc+ABY*I_ERI_Px_S_Dxz_S_C1001_cc;
  Double I_ERI_Py_Py_Dxz_S_C1001_cc = I_ERI_D2y_S_Dxz_S_C1001_cc+ABY*I_ERI_Py_S_Dxz_S_C1001_cc;
  Double I_ERI_Pz_Py_Dxz_S_C1001_cc = I_ERI_Dyz_S_Dxz_S_C1001_cc+ABY*I_ERI_Pz_S_Dxz_S_C1001_cc;
  Double I_ERI_Px_Pz_Dxz_S_C1001_cc = I_ERI_Dxz_S_Dxz_S_C1001_cc+ABZ*I_ERI_Px_S_Dxz_S_C1001_cc;
  Double I_ERI_Py_Pz_Dxz_S_C1001_cc = I_ERI_Dyz_S_Dxz_S_C1001_cc+ABZ*I_ERI_Py_S_Dxz_S_C1001_cc;
  Double I_ERI_Pz_Pz_Dxz_S_C1001_cc = I_ERI_D2z_S_Dxz_S_C1001_cc+ABZ*I_ERI_Pz_S_Dxz_S_C1001_cc;
  Double I_ERI_Px_Px_D2y_S_C1001_cc = I_ERI_D2x_S_D2y_S_C1001_cc+ABX*I_ERI_Px_S_D2y_S_C1001_cc;
  Double I_ERI_Py_Px_D2y_S_C1001_cc = I_ERI_Dxy_S_D2y_S_C1001_cc+ABX*I_ERI_Py_S_D2y_S_C1001_cc;
  Double I_ERI_Pz_Px_D2y_S_C1001_cc = I_ERI_Dxz_S_D2y_S_C1001_cc+ABX*I_ERI_Pz_S_D2y_S_C1001_cc;
  Double I_ERI_Px_Py_D2y_S_C1001_cc = I_ERI_Dxy_S_D2y_S_C1001_cc+ABY*I_ERI_Px_S_D2y_S_C1001_cc;
  Double I_ERI_Py_Py_D2y_S_C1001_cc = I_ERI_D2y_S_D2y_S_C1001_cc+ABY*I_ERI_Py_S_D2y_S_C1001_cc;
  Double I_ERI_Pz_Py_D2y_S_C1001_cc = I_ERI_Dyz_S_D2y_S_C1001_cc+ABY*I_ERI_Pz_S_D2y_S_C1001_cc;
  Double I_ERI_Px_Pz_D2y_S_C1001_cc = I_ERI_Dxz_S_D2y_S_C1001_cc+ABZ*I_ERI_Px_S_D2y_S_C1001_cc;
  Double I_ERI_Py_Pz_D2y_S_C1001_cc = I_ERI_Dyz_S_D2y_S_C1001_cc+ABZ*I_ERI_Py_S_D2y_S_C1001_cc;
  Double I_ERI_Pz_Pz_D2y_S_C1001_cc = I_ERI_D2z_S_D2y_S_C1001_cc+ABZ*I_ERI_Pz_S_D2y_S_C1001_cc;
  Double I_ERI_Px_Px_Dyz_S_C1001_cc = I_ERI_D2x_S_Dyz_S_C1001_cc+ABX*I_ERI_Px_S_Dyz_S_C1001_cc;
  Double I_ERI_Py_Px_Dyz_S_C1001_cc = I_ERI_Dxy_S_Dyz_S_C1001_cc+ABX*I_ERI_Py_S_Dyz_S_C1001_cc;
  Double I_ERI_Pz_Px_Dyz_S_C1001_cc = I_ERI_Dxz_S_Dyz_S_C1001_cc+ABX*I_ERI_Pz_S_Dyz_S_C1001_cc;
  Double I_ERI_Px_Py_Dyz_S_C1001_cc = I_ERI_Dxy_S_Dyz_S_C1001_cc+ABY*I_ERI_Px_S_Dyz_S_C1001_cc;
  Double I_ERI_Py_Py_Dyz_S_C1001_cc = I_ERI_D2y_S_Dyz_S_C1001_cc+ABY*I_ERI_Py_S_Dyz_S_C1001_cc;
  Double I_ERI_Pz_Py_Dyz_S_C1001_cc = I_ERI_Dyz_S_Dyz_S_C1001_cc+ABY*I_ERI_Pz_S_Dyz_S_C1001_cc;
  Double I_ERI_Px_Pz_Dyz_S_C1001_cc = I_ERI_Dxz_S_Dyz_S_C1001_cc+ABZ*I_ERI_Px_S_Dyz_S_C1001_cc;
  Double I_ERI_Py_Pz_Dyz_S_C1001_cc = I_ERI_Dyz_S_Dyz_S_C1001_cc+ABZ*I_ERI_Py_S_Dyz_S_C1001_cc;
  Double I_ERI_Pz_Pz_Dyz_S_C1001_cc = I_ERI_D2z_S_Dyz_S_C1001_cc+ABZ*I_ERI_Pz_S_Dyz_S_C1001_cc;
  Double I_ERI_Px_Px_D2z_S_C1001_cc = I_ERI_D2x_S_D2z_S_C1001_cc+ABX*I_ERI_Px_S_D2z_S_C1001_cc;
  Double I_ERI_Py_Px_D2z_S_C1001_cc = I_ERI_Dxy_S_D2z_S_C1001_cc+ABX*I_ERI_Py_S_D2z_S_C1001_cc;
  Double I_ERI_Pz_Px_D2z_S_C1001_cc = I_ERI_Dxz_S_D2z_S_C1001_cc+ABX*I_ERI_Pz_S_D2z_S_C1001_cc;
  Double I_ERI_Px_Py_D2z_S_C1001_cc = I_ERI_Dxy_S_D2z_S_C1001_cc+ABY*I_ERI_Px_S_D2z_S_C1001_cc;
  Double I_ERI_Py_Py_D2z_S_C1001_cc = I_ERI_D2y_S_D2z_S_C1001_cc+ABY*I_ERI_Py_S_D2z_S_C1001_cc;
  Double I_ERI_Pz_Py_D2z_S_C1001_cc = I_ERI_Dyz_S_D2z_S_C1001_cc+ABY*I_ERI_Pz_S_D2z_S_C1001_cc;
  Double I_ERI_Px_Pz_D2z_S_C1001_cc = I_ERI_Dxz_S_D2z_S_C1001_cc+ABZ*I_ERI_Px_S_D2z_S_C1001_cc;
  Double I_ERI_Py_Pz_D2z_S_C1001_cc = I_ERI_Dyz_S_D2z_S_C1001_cc+ABZ*I_ERI_Py_S_D2z_S_C1001_cc;
  Double I_ERI_Pz_Pz_D2z_S_C1001_cc = I_ERI_D2z_S_D2z_S_C1001_cc+ABZ*I_ERI_Pz_S_D2z_S_C1001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_D2x_S_S_S_C0_aa-2.0E0*1*I_ERI_S_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[1] = 4.0E0*I_ERI_F3x_S_S_S_C1_aa-2.0E0*1*I_ERI_Px_S_S_S_C1_a-2.0E0*2*I_ERI_Px_S_S_S_C1_a;
  abcd[2] = 4.0E0*I_ERI_F2xy_S_S_S_C1_aa-2.0E0*1*I_ERI_Py_S_S_S_C1_a;
  abcd[3] = 4.0E0*I_ERI_F2xz_S_S_S_C1_aa-2.0E0*1*I_ERI_Pz_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[4] = 4.0E0*I_ERI_D2x_Px_S_S_C1000_aa-2.0E0*1*I_ERI_S_Px_S_S_C1000_a;
  abcd[8] = 4.0E0*I_ERI_D2x_Py_S_S_C1000_aa-2.0E0*1*I_ERI_S_Py_S_S_C1000_a;
  abcd[12] = 4.0E0*I_ERI_D2x_Pz_S_S_C1000_aa-2.0E0*1*I_ERI_S_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[5] = 4.0E0*I_ERI_F3x_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Px_S_S_C1001_a-2.0E0*2*I_ERI_Px_Px_S_S_C1001_a;
  abcd[6] = 4.0E0*I_ERI_F2xy_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Px_S_S_C1001_a;
  abcd[7] = 4.0E0*I_ERI_F2xz_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_a;
  abcd[9] = 4.0E0*I_ERI_F3x_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Py_S_S_C1001_a-2.0E0*2*I_ERI_Px_Py_S_S_C1001_a;
  abcd[10] = 4.0E0*I_ERI_F2xy_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Py_S_S_C1001_a;
  abcd[11] = 4.0E0*I_ERI_F2xz_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_a;
  abcd[13] = 4.0E0*I_ERI_F3x_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_a-2.0E0*2*I_ERI_Px_Pz_S_S_C1001_a;
  abcd[14] = 4.0E0*I_ERI_F2xy_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_a;
  abcd[15] = 4.0E0*I_ERI_F2xz_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[16] = 4.0E0*I_ERI_Dxy_S_S_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[17] = 4.0E0*I_ERI_F2xy_S_S_S_C1_aa-2.0E0*1*I_ERI_Py_S_S_S_C1_a;
  abcd[18] = 4.0E0*I_ERI_Fx2y_S_S_S_C1_aa-2.0E0*1*I_ERI_Px_S_S_S_C1_a;
  abcd[19] = 4.0E0*I_ERI_Fxyz_S_S_S_C1_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[20] = 4.0E0*I_ERI_Dxy_Px_S_S_C1000_aa;
  abcd[24] = 4.0E0*I_ERI_Dxy_Py_S_S_C1000_aa;
  abcd[28] = 4.0E0*I_ERI_Dxy_Pz_S_S_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[21] = 4.0E0*I_ERI_F2xy_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Px_S_S_C1001_a;
  abcd[22] = 4.0E0*I_ERI_Fx2y_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Px_S_S_C1001_a;
  abcd[23] = 4.0E0*I_ERI_Fxyz_Px_S_S_C1001_aa;
  abcd[25] = 4.0E0*I_ERI_F2xy_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Py_S_S_C1001_a;
  abcd[26] = 4.0E0*I_ERI_Fx2y_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Py_S_S_C1001_a;
  abcd[27] = 4.0E0*I_ERI_Fxyz_Py_S_S_C1001_aa;
  abcd[29] = 4.0E0*I_ERI_F2xy_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_a;
  abcd[30] = 4.0E0*I_ERI_Fx2y_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_a;
  abcd[31] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C1001_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[32] = 4.0E0*I_ERI_Dxz_S_S_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[33] = 4.0E0*I_ERI_F2xz_S_S_S_C1_aa-2.0E0*1*I_ERI_Pz_S_S_S_C1_a;
  abcd[34] = 4.0E0*I_ERI_Fxyz_S_S_S_C1_aa;
  abcd[35] = 4.0E0*I_ERI_Fx2z_S_S_S_C1_aa-2.0E0*1*I_ERI_Px_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[36] = 4.0E0*I_ERI_Dxz_Px_S_S_C1000_aa;
  abcd[40] = 4.0E0*I_ERI_Dxz_Py_S_S_C1000_aa;
  abcd[44] = 4.0E0*I_ERI_Dxz_Pz_S_S_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[37] = 4.0E0*I_ERI_F2xz_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_a;
  abcd[38] = 4.0E0*I_ERI_Fxyz_Px_S_S_C1001_aa;
  abcd[39] = 4.0E0*I_ERI_Fx2z_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Px_S_S_C1001_a;
  abcd[41] = 4.0E0*I_ERI_F2xz_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_a;
  abcd[42] = 4.0E0*I_ERI_Fxyz_Py_S_S_C1001_aa;
  abcd[43] = 4.0E0*I_ERI_Fx2z_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Py_S_S_C1001_a;
  abcd[45] = 4.0E0*I_ERI_F2xz_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_a;
  abcd[46] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C1001_aa;
  abcd[47] = 4.0E0*I_ERI_Fx2z_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[48] = 4.0E0*I_ERI_D2y_S_S_S_C0_aa-2.0E0*1*I_ERI_S_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[49] = 4.0E0*I_ERI_Fx2y_S_S_S_C1_aa-2.0E0*1*I_ERI_Px_S_S_S_C1_a;
  abcd[50] = 4.0E0*I_ERI_F3y_S_S_S_C1_aa-2.0E0*1*I_ERI_Py_S_S_S_C1_a-2.0E0*2*I_ERI_Py_S_S_S_C1_a;
  abcd[51] = 4.0E0*I_ERI_F2yz_S_S_S_C1_aa-2.0E0*1*I_ERI_Pz_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[52] = 4.0E0*I_ERI_D2y_Px_S_S_C1000_aa-2.0E0*1*I_ERI_S_Px_S_S_C1000_a;
  abcd[56] = 4.0E0*I_ERI_D2y_Py_S_S_C1000_aa-2.0E0*1*I_ERI_S_Py_S_S_C1000_a;
  abcd[60] = 4.0E0*I_ERI_D2y_Pz_S_S_C1000_aa-2.0E0*1*I_ERI_S_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[53] = 4.0E0*I_ERI_Fx2y_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Px_S_S_C1001_a;
  abcd[54] = 4.0E0*I_ERI_F3y_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Px_S_S_C1001_a-2.0E0*2*I_ERI_Py_Px_S_S_C1001_a;
  abcd[55] = 4.0E0*I_ERI_F2yz_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_a;
  abcd[57] = 4.0E0*I_ERI_Fx2y_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Py_S_S_C1001_a;
  abcd[58] = 4.0E0*I_ERI_F3y_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Py_S_S_C1001_a-2.0E0*2*I_ERI_Py_Py_S_S_C1001_a;
  abcd[59] = 4.0E0*I_ERI_F2yz_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_a;
  abcd[61] = 4.0E0*I_ERI_Fx2y_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_a;
  abcd[62] = 4.0E0*I_ERI_F3y_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_a-2.0E0*2*I_ERI_Py_Pz_S_S_C1001_a;
  abcd[63] = 4.0E0*I_ERI_F2yz_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[64] = 4.0E0*I_ERI_Dyz_S_S_S_C0_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[65] = 4.0E0*I_ERI_Fxyz_S_S_S_C1_aa;
  abcd[66] = 4.0E0*I_ERI_F2yz_S_S_S_C1_aa-2.0E0*1*I_ERI_Pz_S_S_S_C1_a;
  abcd[67] = 4.0E0*I_ERI_Fy2z_S_S_S_C1_aa-2.0E0*1*I_ERI_Py_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[68] = 4.0E0*I_ERI_Dyz_Px_S_S_C1000_aa;
  abcd[72] = 4.0E0*I_ERI_Dyz_Py_S_S_C1000_aa;
  abcd[76] = 4.0E0*I_ERI_Dyz_Pz_S_S_C1000_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[69] = 4.0E0*I_ERI_Fxyz_Px_S_S_C1001_aa;
  abcd[70] = 4.0E0*I_ERI_F2yz_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_a;
  abcd[71] = 4.0E0*I_ERI_Fy2z_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Px_S_S_C1001_a;
  abcd[73] = 4.0E0*I_ERI_Fxyz_Py_S_S_C1001_aa;
  abcd[74] = 4.0E0*I_ERI_F2yz_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_a;
  abcd[75] = 4.0E0*I_ERI_Fy2z_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Py_S_S_C1001_a;
  abcd[77] = 4.0E0*I_ERI_Fxyz_Pz_S_S_C1001_aa;
  abcd[78] = 4.0E0*I_ERI_F2yz_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_a;
  abcd[79] = 4.0E0*I_ERI_Fy2z_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C0_aa
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_a
   ************************************************************/
  abcd[80] = 4.0E0*I_ERI_D2z_S_S_S_C0_aa-2.0E0*1*I_ERI_S_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1_aa
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_a
   ************************************************************/
  abcd[81] = 4.0E0*I_ERI_Fx2z_S_S_S_C1_aa-2.0E0*1*I_ERI_Px_S_S_S_C1_a;
  abcd[82] = 4.0E0*I_ERI_Fy2z_S_S_S_C1_aa-2.0E0*1*I_ERI_Py_S_S_S_C1_a;
  abcd[83] = 4.0E0*I_ERI_F3z_S_S_S_C1_aa-2.0E0*1*I_ERI_Pz_S_S_S_C1_a-2.0E0*2*I_ERI_Pz_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1000_aa
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_a
   ************************************************************/
  abcd[84] = 4.0E0*I_ERI_D2z_Px_S_S_C1000_aa-2.0E0*1*I_ERI_S_Px_S_S_C1000_a;
  abcd[88] = 4.0E0*I_ERI_D2z_Py_S_S_C1000_aa-2.0E0*1*I_ERI_S_Py_S_S_C1000_a;
  abcd[92] = 4.0E0*I_ERI_D2z_Pz_S_S_C1000_aa-2.0E0*1*I_ERI_S_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1001_aa
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_a
   ************************************************************/
  abcd[85] = 4.0E0*I_ERI_Fx2z_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Px_S_S_C1001_a;
  abcd[86] = 4.0E0*I_ERI_Fy2z_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Px_S_S_C1001_a;
  abcd[87] = 4.0E0*I_ERI_F3z_Px_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_a-2.0E0*2*I_ERI_Pz_Px_S_S_C1001_a;
  abcd[89] = 4.0E0*I_ERI_Fx2z_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Py_S_S_C1001_a;
  abcd[90] = 4.0E0*I_ERI_Fy2z_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Py_S_S_C1001_a;
  abcd[91] = 4.0E0*I_ERI_F3z_Py_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_a-2.0E0*2*I_ERI_Pz_Py_S_S_C1001_a;
  abcd[93] = 4.0E0*I_ERI_Fx2z_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_a;
  abcd[94] = 4.0E0*I_ERI_Fy2z_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_a;
  abcd[95] = 4.0E0*I_ERI_F3z_Pz_S_S_C1001_aa-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_a-2.0E0*2*I_ERI_Pz_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[96] = 4.0E0*I_ERI_Px_Px_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[97] = 4.0E0*I_ERI_D2x_Px_S_S_C1_ab-2.0E0*1*I_ERI_S_Px_S_S_C1_b;
  abcd[98] = 4.0E0*I_ERI_Dxy_Px_S_S_C1_ab;
  abcd[99] = 4.0E0*I_ERI_Dxz_Px_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[100] = 4.0E0*I_ERI_Px_D2x_S_S_C1000_ab-2.0E0*1*I_ERI_Px_S_S_S_C1000_a;
  abcd[104] = 4.0E0*I_ERI_Px_Dxy_S_S_C1000_ab;
  abcd[108] = 4.0E0*I_ERI_Px_Dxz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[101] = 4.0E0*I_ERI_D2x_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_D2x_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2x_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[102] = 4.0E0*I_ERI_Dxy_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[103] = 4.0E0*I_ERI_Dxz_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;
  abcd[105] = 4.0E0*I_ERI_D2x_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[106] = 4.0E0*I_ERI_Dxy_Dxy_S_S_C1001_ab;
  abcd[107] = 4.0E0*I_ERI_Dxz_Dxy_S_S_C1001_ab;
  abcd[109] = 4.0E0*I_ERI_D2x_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;
  abcd[110] = 4.0E0*I_ERI_Dxy_Dxz_S_S_C1001_ab;
  abcd[111] = 4.0E0*I_ERI_Dxz_Dxz_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[112] = 4.0E0*I_ERI_Px_Py_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[113] = 4.0E0*I_ERI_D2x_Py_S_S_C1_ab-2.0E0*1*I_ERI_S_Py_S_S_C1_b;
  abcd[114] = 4.0E0*I_ERI_Dxy_Py_S_S_C1_ab;
  abcd[115] = 4.0E0*I_ERI_Dxz_Py_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[116] = 4.0E0*I_ERI_Px_Dxy_S_S_C1000_ab;
  abcd[120] = 4.0E0*I_ERI_Px_D2y_S_S_C1000_ab-2.0E0*1*I_ERI_Px_S_S_S_C1000_a;
  abcd[124] = 4.0E0*I_ERI_Px_Dyz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[117] = 4.0E0*I_ERI_D2x_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[118] = 4.0E0*I_ERI_Dxy_Dxy_S_S_C1001_ab;
  abcd[119] = 4.0E0*I_ERI_Dxz_Dxy_S_S_C1001_ab;
  abcd[121] = 4.0E0*I_ERI_D2x_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_D2x_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2y_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[122] = 4.0E0*I_ERI_Dxy_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[123] = 4.0E0*I_ERI_Dxz_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;
  abcd[125] = 4.0E0*I_ERI_D2x_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;
  abcd[126] = 4.0E0*I_ERI_Dxy_Dyz_S_S_C1001_ab;
  abcd[127] = 4.0E0*I_ERI_Dxz_Dyz_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[128] = 4.0E0*I_ERI_Px_Pz_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[129] = 4.0E0*I_ERI_D2x_Pz_S_S_C1_ab-2.0E0*1*I_ERI_S_Pz_S_S_C1_b;
  abcd[130] = 4.0E0*I_ERI_Dxy_Pz_S_S_C1_ab;
  abcd[131] = 4.0E0*I_ERI_Dxz_Pz_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[132] = 4.0E0*I_ERI_Px_Dxz_S_S_C1000_ab;
  abcd[136] = 4.0E0*I_ERI_Px_Dyz_S_S_C1000_ab;
  abcd[140] = 4.0E0*I_ERI_Px_D2z_S_S_C1000_ab-2.0E0*1*I_ERI_Px_S_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[133] = 4.0E0*I_ERI_D2x_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;
  abcd[134] = 4.0E0*I_ERI_Dxy_Dxz_S_S_C1001_ab;
  abcd[135] = 4.0E0*I_ERI_Dxz_Dxz_S_S_C1001_ab;
  abcd[137] = 4.0E0*I_ERI_D2x_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;
  abcd[138] = 4.0E0*I_ERI_Dxy_Dyz_S_S_C1001_ab;
  abcd[139] = 4.0E0*I_ERI_Dxz_Dyz_S_S_C1001_ab;
  abcd[141] = 4.0E0*I_ERI_D2x_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_D2x_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2z_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[142] = 4.0E0*I_ERI_Dxy_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[143] = 4.0E0*I_ERI_Dxz_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[144] = 4.0E0*I_ERI_Py_Px_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[145] = 4.0E0*I_ERI_Dxy_Px_S_S_C1_ab;
  abcd[146] = 4.0E0*I_ERI_D2y_Px_S_S_C1_ab-2.0E0*1*I_ERI_S_Px_S_S_C1_b;
  abcd[147] = 4.0E0*I_ERI_Dyz_Px_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[148] = 4.0E0*I_ERI_Py_D2x_S_S_C1000_ab-2.0E0*1*I_ERI_Py_S_S_S_C1000_a;
  abcd[152] = 4.0E0*I_ERI_Py_Dxy_S_S_C1000_ab;
  abcd[156] = 4.0E0*I_ERI_Py_Dxz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[149] = 4.0E0*I_ERI_Dxy_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[150] = 4.0E0*I_ERI_D2y_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_D2y_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2x_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[151] = 4.0E0*I_ERI_Dyz_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;
  abcd[153] = 4.0E0*I_ERI_Dxy_Dxy_S_S_C1001_ab;
  abcd[154] = 4.0E0*I_ERI_D2y_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[155] = 4.0E0*I_ERI_Dyz_Dxy_S_S_C1001_ab;
  abcd[157] = 4.0E0*I_ERI_Dxy_Dxz_S_S_C1001_ab;
  abcd[158] = 4.0E0*I_ERI_D2y_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;
  abcd[159] = 4.0E0*I_ERI_Dyz_Dxz_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[160] = 4.0E0*I_ERI_Py_Py_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[161] = 4.0E0*I_ERI_Dxy_Py_S_S_C1_ab;
  abcd[162] = 4.0E0*I_ERI_D2y_Py_S_S_C1_ab-2.0E0*1*I_ERI_S_Py_S_S_C1_b;
  abcd[163] = 4.0E0*I_ERI_Dyz_Py_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[164] = 4.0E0*I_ERI_Py_Dxy_S_S_C1000_ab;
  abcd[168] = 4.0E0*I_ERI_Py_D2y_S_S_C1000_ab-2.0E0*1*I_ERI_Py_S_S_S_C1000_a;
  abcd[172] = 4.0E0*I_ERI_Py_Dyz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[165] = 4.0E0*I_ERI_Dxy_Dxy_S_S_C1001_ab;
  abcd[166] = 4.0E0*I_ERI_D2y_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[167] = 4.0E0*I_ERI_Dyz_Dxy_S_S_C1001_ab;
  abcd[169] = 4.0E0*I_ERI_Dxy_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[170] = 4.0E0*I_ERI_D2y_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_D2y_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2y_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[171] = 4.0E0*I_ERI_Dyz_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;
  abcd[173] = 4.0E0*I_ERI_Dxy_Dyz_S_S_C1001_ab;
  abcd[174] = 4.0E0*I_ERI_D2y_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;
  abcd[175] = 4.0E0*I_ERI_Dyz_Dyz_S_S_C1001_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[176] = 4.0E0*I_ERI_Py_Pz_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[177] = 4.0E0*I_ERI_Dxy_Pz_S_S_C1_ab;
  abcd[178] = 4.0E0*I_ERI_D2y_Pz_S_S_C1_ab-2.0E0*1*I_ERI_S_Pz_S_S_C1_b;
  abcd[179] = 4.0E0*I_ERI_Dyz_Pz_S_S_C1_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_Py_Dxz_S_S_C1000_ab;
  abcd[184] = 4.0E0*I_ERI_Py_Dyz_S_S_C1000_ab;
  abcd[188] = 4.0E0*I_ERI_Py_D2z_S_S_C1000_ab-2.0E0*1*I_ERI_Py_S_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[181] = 4.0E0*I_ERI_Dxy_Dxz_S_S_C1001_ab;
  abcd[182] = 4.0E0*I_ERI_D2y_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;
  abcd[183] = 4.0E0*I_ERI_Dyz_Dxz_S_S_C1001_ab;
  abcd[185] = 4.0E0*I_ERI_Dxy_Dyz_S_S_C1001_ab;
  abcd[186] = 4.0E0*I_ERI_D2y_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;
  abcd[187] = 4.0E0*I_ERI_Dyz_Dyz_S_S_C1001_ab;
  abcd[189] = 4.0E0*I_ERI_Dxy_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dxy_S_S_S_C1001_a;
  abcd[190] = 4.0E0*I_ERI_D2y_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_D2y_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2z_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[191] = 4.0E0*I_ERI_Dyz_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[192] = 4.0E0*I_ERI_Pz_Px_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[193] = 4.0E0*I_ERI_Dxz_Px_S_S_C1_ab;
  abcd[194] = 4.0E0*I_ERI_Dyz_Px_S_S_C1_ab;
  abcd[195] = 4.0E0*I_ERI_D2z_Px_S_S_C1_ab-2.0E0*1*I_ERI_S_Px_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[196] = 4.0E0*I_ERI_Pz_D2x_S_S_C1000_ab-2.0E0*1*I_ERI_Pz_S_S_S_C1000_a;
  abcd[200] = 4.0E0*I_ERI_Pz_Dxy_S_S_C1000_ab;
  abcd[204] = 4.0E0*I_ERI_Pz_Dxz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[197] = 4.0E0*I_ERI_Dxz_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;
  abcd[198] = 4.0E0*I_ERI_Dyz_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;
  abcd[199] = 4.0E0*I_ERI_D2z_D2x_S_S_C1001_ab-2.0E0*1*I_ERI_D2z_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2x_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[201] = 4.0E0*I_ERI_Dxz_Dxy_S_S_C1001_ab;
  abcd[202] = 4.0E0*I_ERI_Dyz_Dxy_S_S_C1001_ab;
  abcd[203] = 4.0E0*I_ERI_D2z_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[205] = 4.0E0*I_ERI_Dxz_Dxz_S_S_C1001_ab;
  abcd[206] = 4.0E0*I_ERI_Dyz_Dxz_S_S_C1001_ab;
  abcd[207] = 4.0E0*I_ERI_D2z_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[208] = 4.0E0*I_ERI_Pz_Py_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[209] = 4.0E0*I_ERI_Dxz_Py_S_S_C1_ab;
  abcd[210] = 4.0E0*I_ERI_Dyz_Py_S_S_C1_ab;
  abcd[211] = 4.0E0*I_ERI_D2z_Py_S_S_C1_ab-2.0E0*1*I_ERI_S_Py_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[212] = 4.0E0*I_ERI_Pz_Dxy_S_S_C1000_ab;
  abcd[216] = 4.0E0*I_ERI_Pz_D2y_S_S_C1000_ab-2.0E0*1*I_ERI_Pz_S_S_S_C1000_a;
  abcd[220] = 4.0E0*I_ERI_Pz_Dyz_S_S_C1000_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[213] = 4.0E0*I_ERI_Dxz_Dxy_S_S_C1001_ab;
  abcd[214] = 4.0E0*I_ERI_Dyz_Dxy_S_S_C1001_ab;
  abcd[215] = 4.0E0*I_ERI_D2z_Dxy_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxy_S_S_C1001_b;
  abcd[217] = 4.0E0*I_ERI_Dxz_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;
  abcd[218] = 4.0E0*I_ERI_Dyz_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;
  abcd[219] = 4.0E0*I_ERI_D2z_D2y_S_S_C1001_ab-2.0E0*1*I_ERI_D2z_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2y_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;
  abcd[221] = 4.0E0*I_ERI_Dxz_Dyz_S_S_C1001_ab;
  abcd[222] = 4.0E0*I_ERI_Dyz_Dyz_S_S_C1001_ab;
  abcd[223] = 4.0E0*I_ERI_D2z_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C0_ab
   ************************************************************/
  abcd[224] = 4.0E0*I_ERI_Pz_Pz_S_S_C0_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1_ab
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1_b
   ************************************************************/
  abcd[225] = 4.0E0*I_ERI_Dxz_Pz_S_S_C1_ab;
  abcd[226] = 4.0E0*I_ERI_Dyz_Pz_S_S_C1_ab;
  abcd[227] = 4.0E0*I_ERI_D2z_Pz_S_S_C1_ab-2.0E0*1*I_ERI_S_Pz_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1000_ab
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  abcd[228] = 4.0E0*I_ERI_Pz_Dxz_S_S_C1000_ab;
  abcd[232] = 4.0E0*I_ERI_Pz_Dyz_S_S_C1000_ab;
  abcd[236] = 4.0E0*I_ERI_Pz_D2z_S_S_C1000_ab-2.0E0*1*I_ERI_Pz_S_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1001_ab
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1001
   ************************************************************/
  abcd[229] = 4.0E0*I_ERI_Dxz_Dxz_S_S_C1001_ab;
  abcd[230] = 4.0E0*I_ERI_Dyz_Dxz_S_S_C1001_ab;
  abcd[231] = 4.0E0*I_ERI_D2z_Dxz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dxz_S_S_C1001_b;
  abcd[233] = 4.0E0*I_ERI_Dxz_Dyz_S_S_C1001_ab;
  abcd[234] = 4.0E0*I_ERI_Dyz_Dyz_S_S_C1001_ab;
  abcd[235] = 4.0E0*I_ERI_D2z_Dyz_S_S_C1001_ab-2.0E0*1*I_ERI_S_Dyz_S_S_C1001_b;
  abcd[237] = 4.0E0*I_ERI_Dxz_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dxz_S_S_S_C1001_a;
  abcd[238] = 4.0E0*I_ERI_Dyz_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_Dyz_S_S_S_C1001_a;
  abcd[239] = 4.0E0*I_ERI_D2z_D2z_S_S_C1001_ab-2.0E0*1*I_ERI_D2z_S_S_S_C1001_a-2.0E0*1*I_ERI_S_D2z_S_S_C1001_b+1*I_ERI_S_S_S_S_C1001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_Px_S_Px_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[241] = 4.0E0*I_ERI_D2x_S_Px_S_C1_ac-2.0E0*1*I_ERI_S_S_Px_S_C1_c;
  abcd[242] = 4.0E0*I_ERI_Dxy_S_Px_S_C1_ac;
  abcd[243] = 4.0E0*I_ERI_Dxz_S_Px_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[244] = 4.0E0*I_ERI_Px_Px_Px_S_C1000_ac;
  abcd[248] = 4.0E0*I_ERI_Px_Py_Px_S_C1000_ac;
  abcd[252] = 4.0E0*I_ERI_Px_Pz_Px_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[245] = 4.0E0*I_ERI_D2x_Px_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Px_S_C1001_c;
  abcd[246] = 4.0E0*I_ERI_Dxy_Px_Px_S_C1001_ac;
  abcd[247] = 4.0E0*I_ERI_Dxz_Px_Px_S_C1001_ac;
  abcd[249] = 4.0E0*I_ERI_D2x_Py_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Px_S_C1001_c;
  abcd[250] = 4.0E0*I_ERI_Dxy_Py_Px_S_C1001_ac;
  abcd[251] = 4.0E0*I_ERI_Dxz_Py_Px_S_C1001_ac;
  abcd[253] = 4.0E0*I_ERI_D2x_Pz_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Px_S_C1001_c;
  abcd[254] = 4.0E0*I_ERI_Dxy_Pz_Px_S_C1001_ac;
  abcd[255] = 4.0E0*I_ERI_Dxz_Pz_Px_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[256] = 4.0E0*I_ERI_Px_S_Py_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[257] = 4.0E0*I_ERI_D2x_S_Py_S_C1_ac-2.0E0*1*I_ERI_S_S_Py_S_C1_c;
  abcd[258] = 4.0E0*I_ERI_Dxy_S_Py_S_C1_ac;
  abcd[259] = 4.0E0*I_ERI_Dxz_S_Py_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[260] = 4.0E0*I_ERI_Px_Px_Py_S_C1000_ac;
  abcd[264] = 4.0E0*I_ERI_Px_Py_Py_S_C1000_ac;
  abcd[268] = 4.0E0*I_ERI_Px_Pz_Py_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[261] = 4.0E0*I_ERI_D2x_Px_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Py_S_C1001_c;
  abcd[262] = 4.0E0*I_ERI_Dxy_Px_Py_S_C1001_ac;
  abcd[263] = 4.0E0*I_ERI_Dxz_Px_Py_S_C1001_ac;
  abcd[265] = 4.0E0*I_ERI_D2x_Py_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Py_S_C1001_c;
  abcd[266] = 4.0E0*I_ERI_Dxy_Py_Py_S_C1001_ac;
  abcd[267] = 4.0E0*I_ERI_Dxz_Py_Py_S_C1001_ac;
  abcd[269] = 4.0E0*I_ERI_D2x_Pz_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Py_S_C1001_c;
  abcd[270] = 4.0E0*I_ERI_Dxy_Pz_Py_S_C1001_ac;
  abcd[271] = 4.0E0*I_ERI_Dxz_Pz_Py_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[272] = 4.0E0*I_ERI_Px_S_Pz_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[273] = 4.0E0*I_ERI_D2x_S_Pz_S_C1_ac-2.0E0*1*I_ERI_S_S_Pz_S_C1_c;
  abcd[274] = 4.0E0*I_ERI_Dxy_S_Pz_S_C1_ac;
  abcd[275] = 4.0E0*I_ERI_Dxz_S_Pz_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[276] = 4.0E0*I_ERI_Px_Px_Pz_S_C1000_ac;
  abcd[280] = 4.0E0*I_ERI_Px_Py_Pz_S_C1000_ac;
  abcd[284] = 4.0E0*I_ERI_Px_Pz_Pz_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[277] = 4.0E0*I_ERI_D2x_Px_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Pz_S_C1001_c;
  abcd[278] = 4.0E0*I_ERI_Dxy_Px_Pz_S_C1001_ac;
  abcd[279] = 4.0E0*I_ERI_Dxz_Px_Pz_S_C1001_ac;
  abcd[281] = 4.0E0*I_ERI_D2x_Py_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Pz_S_C1001_c;
  abcd[282] = 4.0E0*I_ERI_Dxy_Py_Pz_S_C1001_ac;
  abcd[283] = 4.0E0*I_ERI_Dxz_Py_Pz_S_C1001_ac;
  abcd[285] = 4.0E0*I_ERI_D2x_Pz_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Pz_S_C1001_c;
  abcd[286] = 4.0E0*I_ERI_Dxy_Pz_Pz_S_C1001_ac;
  abcd[287] = 4.0E0*I_ERI_Dxz_Pz_Pz_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[288] = 4.0E0*I_ERI_Py_S_Px_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[289] = 4.0E0*I_ERI_Dxy_S_Px_S_C1_ac;
  abcd[290] = 4.0E0*I_ERI_D2y_S_Px_S_C1_ac-2.0E0*1*I_ERI_S_S_Px_S_C1_c;
  abcd[291] = 4.0E0*I_ERI_Dyz_S_Px_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[292] = 4.0E0*I_ERI_Py_Px_Px_S_C1000_ac;
  abcd[296] = 4.0E0*I_ERI_Py_Py_Px_S_C1000_ac;
  abcd[300] = 4.0E0*I_ERI_Py_Pz_Px_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[293] = 4.0E0*I_ERI_Dxy_Px_Px_S_C1001_ac;
  abcd[294] = 4.0E0*I_ERI_D2y_Px_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Px_S_C1001_c;
  abcd[295] = 4.0E0*I_ERI_Dyz_Px_Px_S_C1001_ac;
  abcd[297] = 4.0E0*I_ERI_Dxy_Py_Px_S_C1001_ac;
  abcd[298] = 4.0E0*I_ERI_D2y_Py_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Px_S_C1001_c;
  abcd[299] = 4.0E0*I_ERI_Dyz_Py_Px_S_C1001_ac;
  abcd[301] = 4.0E0*I_ERI_Dxy_Pz_Px_S_C1001_ac;
  abcd[302] = 4.0E0*I_ERI_D2y_Pz_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Px_S_C1001_c;
  abcd[303] = 4.0E0*I_ERI_Dyz_Pz_Px_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[304] = 4.0E0*I_ERI_Py_S_Py_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[305] = 4.0E0*I_ERI_Dxy_S_Py_S_C1_ac;
  abcd[306] = 4.0E0*I_ERI_D2y_S_Py_S_C1_ac-2.0E0*1*I_ERI_S_S_Py_S_C1_c;
  abcd[307] = 4.0E0*I_ERI_Dyz_S_Py_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[308] = 4.0E0*I_ERI_Py_Px_Py_S_C1000_ac;
  abcd[312] = 4.0E0*I_ERI_Py_Py_Py_S_C1000_ac;
  abcd[316] = 4.0E0*I_ERI_Py_Pz_Py_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[309] = 4.0E0*I_ERI_Dxy_Px_Py_S_C1001_ac;
  abcd[310] = 4.0E0*I_ERI_D2y_Px_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Py_S_C1001_c;
  abcd[311] = 4.0E0*I_ERI_Dyz_Px_Py_S_C1001_ac;
  abcd[313] = 4.0E0*I_ERI_Dxy_Py_Py_S_C1001_ac;
  abcd[314] = 4.0E0*I_ERI_D2y_Py_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Py_S_C1001_c;
  abcd[315] = 4.0E0*I_ERI_Dyz_Py_Py_S_C1001_ac;
  abcd[317] = 4.0E0*I_ERI_Dxy_Pz_Py_S_C1001_ac;
  abcd[318] = 4.0E0*I_ERI_D2y_Pz_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Py_S_C1001_c;
  abcd[319] = 4.0E0*I_ERI_Dyz_Pz_Py_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[320] = 4.0E0*I_ERI_Py_S_Pz_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[321] = 4.0E0*I_ERI_Dxy_S_Pz_S_C1_ac;
  abcd[322] = 4.0E0*I_ERI_D2y_S_Pz_S_C1_ac-2.0E0*1*I_ERI_S_S_Pz_S_C1_c;
  abcd[323] = 4.0E0*I_ERI_Dyz_S_Pz_S_C1_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[324] = 4.0E0*I_ERI_Py_Px_Pz_S_C1000_ac;
  abcd[328] = 4.0E0*I_ERI_Py_Py_Pz_S_C1000_ac;
  abcd[332] = 4.0E0*I_ERI_Py_Pz_Pz_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[325] = 4.0E0*I_ERI_Dxy_Px_Pz_S_C1001_ac;
  abcd[326] = 4.0E0*I_ERI_D2y_Px_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Pz_S_C1001_c;
  abcd[327] = 4.0E0*I_ERI_Dyz_Px_Pz_S_C1001_ac;
  abcd[329] = 4.0E0*I_ERI_Dxy_Py_Pz_S_C1001_ac;
  abcd[330] = 4.0E0*I_ERI_D2y_Py_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Pz_S_C1001_c;
  abcd[331] = 4.0E0*I_ERI_Dyz_Py_Pz_S_C1001_ac;
  abcd[333] = 4.0E0*I_ERI_Dxy_Pz_Pz_S_C1001_ac;
  abcd[334] = 4.0E0*I_ERI_D2y_Pz_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Pz_S_C1001_c;
  abcd[335] = 4.0E0*I_ERI_Dyz_Pz_Pz_S_C1001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[336] = 4.0E0*I_ERI_Pz_S_Px_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[337] = 4.0E0*I_ERI_Dxz_S_Px_S_C1_ac;
  abcd[338] = 4.0E0*I_ERI_Dyz_S_Px_S_C1_ac;
  abcd[339] = 4.0E0*I_ERI_D2z_S_Px_S_C1_ac-2.0E0*1*I_ERI_S_S_Px_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[340] = 4.0E0*I_ERI_Pz_Px_Px_S_C1000_ac;
  abcd[344] = 4.0E0*I_ERI_Pz_Py_Px_S_C1000_ac;
  abcd[348] = 4.0E0*I_ERI_Pz_Pz_Px_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[341] = 4.0E0*I_ERI_Dxz_Px_Px_S_C1001_ac;
  abcd[342] = 4.0E0*I_ERI_Dyz_Px_Px_S_C1001_ac;
  abcd[343] = 4.0E0*I_ERI_D2z_Px_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Px_S_C1001_c;
  abcd[345] = 4.0E0*I_ERI_Dxz_Py_Px_S_C1001_ac;
  abcd[346] = 4.0E0*I_ERI_Dyz_Py_Px_S_C1001_ac;
  abcd[347] = 4.0E0*I_ERI_D2z_Py_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Px_S_C1001_c;
  abcd[349] = 4.0E0*I_ERI_Dxz_Pz_Px_S_C1001_ac;
  abcd[350] = 4.0E0*I_ERI_Dyz_Pz_Px_S_C1001_ac;
  abcd[351] = 4.0E0*I_ERI_D2z_Pz_Px_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Px_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[352] = 4.0E0*I_ERI_Pz_S_Py_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[353] = 4.0E0*I_ERI_Dxz_S_Py_S_C1_ac;
  abcd[354] = 4.0E0*I_ERI_Dyz_S_Py_S_C1_ac;
  abcd[355] = 4.0E0*I_ERI_D2z_S_Py_S_C1_ac-2.0E0*1*I_ERI_S_S_Py_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[356] = 4.0E0*I_ERI_Pz_Px_Py_S_C1000_ac;
  abcd[360] = 4.0E0*I_ERI_Pz_Py_Py_S_C1000_ac;
  abcd[364] = 4.0E0*I_ERI_Pz_Pz_Py_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[357] = 4.0E0*I_ERI_Dxz_Px_Py_S_C1001_ac;
  abcd[358] = 4.0E0*I_ERI_Dyz_Px_Py_S_C1001_ac;
  abcd[359] = 4.0E0*I_ERI_D2z_Px_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Py_S_C1001_c;
  abcd[361] = 4.0E0*I_ERI_Dxz_Py_Py_S_C1001_ac;
  abcd[362] = 4.0E0*I_ERI_Dyz_Py_Py_S_C1001_ac;
  abcd[363] = 4.0E0*I_ERI_D2z_Py_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Py_S_C1001_c;
  abcd[365] = 4.0E0*I_ERI_Dxz_Pz_Py_S_C1001_ac;
  abcd[366] = 4.0E0*I_ERI_Dyz_Pz_Py_S_C1001_ac;
  abcd[367] = 4.0E0*I_ERI_D2z_Pz_Py_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Py_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C0_ac
   ************************************************************/
  abcd[368] = 4.0E0*I_ERI_Pz_S_Pz_S_C0_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1_ac
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1_c
   ************************************************************/
  abcd[369] = 4.0E0*I_ERI_Dxz_S_Pz_S_C1_ac;
  abcd[370] = 4.0E0*I_ERI_Dyz_S_Pz_S_C1_ac;
  abcd[371] = 4.0E0*I_ERI_D2z_S_Pz_S_C1_ac-2.0E0*1*I_ERI_S_S_Pz_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000_ac
   ************************************************************/
  abcd[372] = 4.0E0*I_ERI_Pz_Px_Pz_S_C1000_ac;
  abcd[376] = 4.0E0*I_ERI_Pz_Py_Pz_S_C1000_ac;
  abcd[380] = 4.0E0*I_ERI_Pz_Pz_Pz_S_C1000_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001_ac
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001_c
   ************************************************************/
  abcd[373] = 4.0E0*I_ERI_Dxz_Px_Pz_S_C1001_ac;
  abcd[374] = 4.0E0*I_ERI_Dyz_Px_Pz_S_C1001_ac;
  abcd[375] = 4.0E0*I_ERI_D2z_Px_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Px_Pz_S_C1001_c;
  abcd[377] = 4.0E0*I_ERI_Dxz_Py_Pz_S_C1001_ac;
  abcd[378] = 4.0E0*I_ERI_Dyz_Py_Pz_S_C1001_ac;
  abcd[379] = 4.0E0*I_ERI_D2z_Py_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Py_Pz_S_C1001_c;
  abcd[381] = 4.0E0*I_ERI_Dxz_Pz_Pz_S_C1001_ac;
  abcd[382] = 4.0E0*I_ERI_Dyz_Pz_Pz_S_C1001_ac;
  abcd[383] = 4.0E0*I_ERI_D2z_Pz_Pz_S_C1001_ac-2.0E0*1*I_ERI_S_Pz_Pz_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[384] = 4.0E0*I_ERI_S_D2x_S_S_C0_bb-2.0E0*1*I_ERI_S_S_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[385] = 4.0E0*I_ERI_Px_D2x_S_S_C1_bb-2.0E0*1*I_ERI_Px_S_S_S_C1_b;
  abcd[386] = 4.0E0*I_ERI_Py_D2x_S_S_C1_bb-2.0E0*1*I_ERI_Py_S_S_S_C1_b;
  abcd[387] = 4.0E0*I_ERI_Pz_D2x_S_S_C1_bb-2.0E0*1*I_ERI_Pz_S_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[388] = 4.0E0*I_ERI_S_F3x_S_S_C1000_bb-2.0E0*1*I_ERI_S_Px_S_S_C1000_b-2.0E0*2*I_ERI_S_Px_S_S_C1000_b;
  abcd[392] = 4.0E0*I_ERI_S_F2xy_S_S_C1000_bb-2.0E0*1*I_ERI_S_Py_S_S_C1000_b;
  abcd[396] = 4.0E0*I_ERI_S_F2xz_S_S_C1000_bb-2.0E0*1*I_ERI_S_Pz_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[389] = 4.0E0*I_ERI_Px_F3x_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Px_S_S_C1001_b-2.0E0*2*I_ERI_Px_Px_S_S_C1001_b;
  abcd[390] = 4.0E0*I_ERI_Py_F3x_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Px_S_S_C1001_b-2.0E0*2*I_ERI_Py_Px_S_S_C1001_b;
  abcd[391] = 4.0E0*I_ERI_Pz_F3x_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_b-2.0E0*2*I_ERI_Pz_Px_S_S_C1001_b;
  abcd[393] = 4.0E0*I_ERI_Px_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Py_S_S_C1001_b;
  abcd[394] = 4.0E0*I_ERI_Py_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Py_S_S_C1001_b;
  abcd[395] = 4.0E0*I_ERI_Pz_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_b;
  abcd[397] = 4.0E0*I_ERI_Px_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_b;
  abcd[398] = 4.0E0*I_ERI_Py_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_b;
  abcd[399] = 4.0E0*I_ERI_Pz_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[400] = 4.0E0*I_ERI_S_Dxy_S_S_C0_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[401] = 4.0E0*I_ERI_Px_Dxy_S_S_C1_bb;
  abcd[402] = 4.0E0*I_ERI_Py_Dxy_S_S_C1_bb;
  abcd[403] = 4.0E0*I_ERI_Pz_Dxy_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[404] = 4.0E0*I_ERI_S_F2xy_S_S_C1000_bb-2.0E0*1*I_ERI_S_Py_S_S_C1000_b;
  abcd[408] = 4.0E0*I_ERI_S_Fx2y_S_S_C1000_bb-2.0E0*1*I_ERI_S_Px_S_S_C1000_b;
  abcd[412] = 4.0E0*I_ERI_S_Fxyz_S_S_C1000_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[405] = 4.0E0*I_ERI_Px_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Py_S_S_C1001_b;
  abcd[406] = 4.0E0*I_ERI_Py_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Py_S_S_C1001_b;
  abcd[407] = 4.0E0*I_ERI_Pz_F2xy_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_b;
  abcd[409] = 4.0E0*I_ERI_Px_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Px_S_S_C1001_b;
  abcd[410] = 4.0E0*I_ERI_Py_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Px_S_S_C1001_b;
  abcd[411] = 4.0E0*I_ERI_Pz_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_b;
  abcd[413] = 4.0E0*I_ERI_Px_Fxyz_S_S_C1001_bb;
  abcd[414] = 4.0E0*I_ERI_Py_Fxyz_S_S_C1001_bb;
  abcd[415] = 4.0E0*I_ERI_Pz_Fxyz_S_S_C1001_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[416] = 4.0E0*I_ERI_S_Dxz_S_S_C0_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[417] = 4.0E0*I_ERI_Px_Dxz_S_S_C1_bb;
  abcd[418] = 4.0E0*I_ERI_Py_Dxz_S_S_C1_bb;
  abcd[419] = 4.0E0*I_ERI_Pz_Dxz_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[420] = 4.0E0*I_ERI_S_F2xz_S_S_C1000_bb-2.0E0*1*I_ERI_S_Pz_S_S_C1000_b;
  abcd[424] = 4.0E0*I_ERI_S_Fxyz_S_S_C1000_bb;
  abcd[428] = 4.0E0*I_ERI_S_Fx2z_S_S_C1000_bb-2.0E0*1*I_ERI_S_Px_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[421] = 4.0E0*I_ERI_Px_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_b;
  abcd[422] = 4.0E0*I_ERI_Py_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_b;
  abcd[423] = 4.0E0*I_ERI_Pz_F2xz_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_b;
  abcd[425] = 4.0E0*I_ERI_Px_Fxyz_S_S_C1001_bb;
  abcd[426] = 4.0E0*I_ERI_Py_Fxyz_S_S_C1001_bb;
  abcd[427] = 4.0E0*I_ERI_Pz_Fxyz_S_S_C1001_bb;
  abcd[429] = 4.0E0*I_ERI_Px_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Px_S_S_C1001_b;
  abcd[430] = 4.0E0*I_ERI_Py_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Px_S_S_C1001_b;
  abcd[431] = 4.0E0*I_ERI_Pz_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[432] = 4.0E0*I_ERI_S_D2y_S_S_C0_bb-2.0E0*1*I_ERI_S_S_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[433] = 4.0E0*I_ERI_Px_D2y_S_S_C1_bb-2.0E0*1*I_ERI_Px_S_S_S_C1_b;
  abcd[434] = 4.0E0*I_ERI_Py_D2y_S_S_C1_bb-2.0E0*1*I_ERI_Py_S_S_S_C1_b;
  abcd[435] = 4.0E0*I_ERI_Pz_D2y_S_S_C1_bb-2.0E0*1*I_ERI_Pz_S_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[436] = 4.0E0*I_ERI_S_Fx2y_S_S_C1000_bb-2.0E0*1*I_ERI_S_Px_S_S_C1000_b;
  abcd[440] = 4.0E0*I_ERI_S_F3y_S_S_C1000_bb-2.0E0*1*I_ERI_S_Py_S_S_C1000_b-2.0E0*2*I_ERI_S_Py_S_S_C1000_b;
  abcd[444] = 4.0E0*I_ERI_S_F2yz_S_S_C1000_bb-2.0E0*1*I_ERI_S_Pz_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[437] = 4.0E0*I_ERI_Px_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Px_S_S_C1001_b;
  abcd[438] = 4.0E0*I_ERI_Py_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Px_S_S_C1001_b;
  abcd[439] = 4.0E0*I_ERI_Pz_Fx2y_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_b;
  abcd[441] = 4.0E0*I_ERI_Px_F3y_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Py_S_S_C1001_b-2.0E0*2*I_ERI_Px_Py_S_S_C1001_b;
  abcd[442] = 4.0E0*I_ERI_Py_F3y_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Py_S_S_C1001_b-2.0E0*2*I_ERI_Py_Py_S_S_C1001_b;
  abcd[443] = 4.0E0*I_ERI_Pz_F3y_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_b-2.0E0*2*I_ERI_Pz_Py_S_S_C1001_b;
  abcd[445] = 4.0E0*I_ERI_Px_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_b;
  abcd[446] = 4.0E0*I_ERI_Py_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_b;
  abcd[447] = 4.0E0*I_ERI_Pz_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[448] = 4.0E0*I_ERI_S_Dyz_S_S_C0_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[449] = 4.0E0*I_ERI_Px_Dyz_S_S_C1_bb;
  abcd[450] = 4.0E0*I_ERI_Py_Dyz_S_S_C1_bb;
  abcd[451] = 4.0E0*I_ERI_Pz_Dyz_S_S_C1_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[452] = 4.0E0*I_ERI_S_Fxyz_S_S_C1000_bb;
  abcd[456] = 4.0E0*I_ERI_S_F2yz_S_S_C1000_bb-2.0E0*1*I_ERI_S_Pz_S_S_C1000_b;
  abcd[460] = 4.0E0*I_ERI_S_Fy2z_S_S_C1000_bb-2.0E0*1*I_ERI_S_Py_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[453] = 4.0E0*I_ERI_Px_Fxyz_S_S_C1001_bb;
  abcd[454] = 4.0E0*I_ERI_Py_Fxyz_S_S_C1001_bb;
  abcd[455] = 4.0E0*I_ERI_Pz_Fxyz_S_S_C1001_bb;
  abcd[457] = 4.0E0*I_ERI_Px_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_b;
  abcd[458] = 4.0E0*I_ERI_Py_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_b;
  abcd[459] = 4.0E0*I_ERI_Pz_F2yz_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_b;
  abcd[461] = 4.0E0*I_ERI_Px_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Py_S_S_C1001_b;
  abcd[462] = 4.0E0*I_ERI_Py_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Py_S_S_C1001_b;
  abcd[463] = 4.0E0*I_ERI_Pz_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C0_bb
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_b
   ************************************************************/
  abcd[464] = 4.0E0*I_ERI_S_D2z_S_S_C0_bb-2.0E0*1*I_ERI_S_S_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1_bb
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  abcd[465] = 4.0E0*I_ERI_Px_D2z_S_S_C1_bb-2.0E0*1*I_ERI_Px_S_S_S_C1_b;
  abcd[466] = 4.0E0*I_ERI_Py_D2z_S_S_C1_bb-2.0E0*1*I_ERI_Py_S_S_S_C1_b;
  abcd[467] = 4.0E0*I_ERI_Pz_D2z_S_S_C1_bb-2.0E0*1*I_ERI_Pz_S_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_F_S_S_C1000_bb
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_b
   ************************************************************/
  abcd[468] = 4.0E0*I_ERI_S_Fx2z_S_S_C1000_bb-2.0E0*1*I_ERI_S_Px_S_S_C1000_b;
  abcd[472] = 4.0E0*I_ERI_S_Fy2z_S_S_C1000_bb-2.0E0*1*I_ERI_S_Py_S_S_C1000_b;
  abcd[476] = 4.0E0*I_ERI_S_F3z_S_S_C1000_bb-2.0E0*1*I_ERI_S_Pz_S_S_C1000_b-2.0E0*2*I_ERI_S_Pz_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_F_S_S_C1001_bb
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  abcd[469] = 4.0E0*I_ERI_Px_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Px_S_S_C1001_b;
  abcd[470] = 4.0E0*I_ERI_Py_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Px_S_S_C1001_b;
  abcd[471] = 4.0E0*I_ERI_Pz_Fx2z_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_b;
  abcd[473] = 4.0E0*I_ERI_Px_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Py_S_S_C1001_b;
  abcd[474] = 4.0E0*I_ERI_Py_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Py_S_S_C1001_b;
  abcd[475] = 4.0E0*I_ERI_Pz_Fy2z_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_b;
  abcd[477] = 4.0E0*I_ERI_Px_F3z_S_S_C1001_bb-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_b-2.0E0*2*I_ERI_Px_Pz_S_S_C1001_b;
  abcd[478] = 4.0E0*I_ERI_Py_F3z_S_S_C1001_bb-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_b-2.0E0*2*I_ERI_Py_Pz_S_S_C1001_b;
  abcd[479] = 4.0E0*I_ERI_Pz_F3z_S_S_C1001_bb-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_b-2.0E0*2*I_ERI_Pz_Pz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[480] = 4.0E0*I_ERI_S_Px_Px_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[481] = 4.0E0*I_ERI_Px_Px_Px_S_C1_bc;
  abcd[482] = 4.0E0*I_ERI_Py_Px_Px_S_C1_bc;
  abcd[483] = 4.0E0*I_ERI_Pz_Px_Px_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[484] = 4.0E0*I_ERI_S_D2x_Px_S_C1000_bc-2.0E0*1*I_ERI_S_S_Px_S_C1000_c;
  abcd[488] = 4.0E0*I_ERI_S_Dxy_Px_S_C1000_bc;
  abcd[492] = 4.0E0*I_ERI_S_Dxz_Px_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[485] = 4.0E0*I_ERI_Px_D2x_Px_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Px_S_C1001_c;
  abcd[486] = 4.0E0*I_ERI_Py_D2x_Px_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Px_S_C1001_c;
  abcd[487] = 4.0E0*I_ERI_Pz_D2x_Px_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Px_S_C1001_c;
  abcd[489] = 4.0E0*I_ERI_Px_Dxy_Px_S_C1001_bc;
  abcd[490] = 4.0E0*I_ERI_Py_Dxy_Px_S_C1001_bc;
  abcd[491] = 4.0E0*I_ERI_Pz_Dxy_Px_S_C1001_bc;
  abcd[493] = 4.0E0*I_ERI_Px_Dxz_Px_S_C1001_bc;
  abcd[494] = 4.0E0*I_ERI_Py_Dxz_Px_S_C1001_bc;
  abcd[495] = 4.0E0*I_ERI_Pz_Dxz_Px_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[496] = 4.0E0*I_ERI_S_Px_Py_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[497] = 4.0E0*I_ERI_Px_Px_Py_S_C1_bc;
  abcd[498] = 4.0E0*I_ERI_Py_Px_Py_S_C1_bc;
  abcd[499] = 4.0E0*I_ERI_Pz_Px_Py_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[500] = 4.0E0*I_ERI_S_D2x_Py_S_C1000_bc-2.0E0*1*I_ERI_S_S_Py_S_C1000_c;
  abcd[504] = 4.0E0*I_ERI_S_Dxy_Py_S_C1000_bc;
  abcd[508] = 4.0E0*I_ERI_S_Dxz_Py_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[501] = 4.0E0*I_ERI_Px_D2x_Py_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Py_S_C1001_c;
  abcd[502] = 4.0E0*I_ERI_Py_D2x_Py_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Py_S_C1001_c;
  abcd[503] = 4.0E0*I_ERI_Pz_D2x_Py_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Py_S_C1001_c;
  abcd[505] = 4.0E0*I_ERI_Px_Dxy_Py_S_C1001_bc;
  abcd[506] = 4.0E0*I_ERI_Py_Dxy_Py_S_C1001_bc;
  abcd[507] = 4.0E0*I_ERI_Pz_Dxy_Py_S_C1001_bc;
  abcd[509] = 4.0E0*I_ERI_Px_Dxz_Py_S_C1001_bc;
  abcd[510] = 4.0E0*I_ERI_Py_Dxz_Py_S_C1001_bc;
  abcd[511] = 4.0E0*I_ERI_Pz_Dxz_Py_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[512] = 4.0E0*I_ERI_S_Px_Pz_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[513] = 4.0E0*I_ERI_Px_Px_Pz_S_C1_bc;
  abcd[514] = 4.0E0*I_ERI_Py_Px_Pz_S_C1_bc;
  abcd[515] = 4.0E0*I_ERI_Pz_Px_Pz_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[516] = 4.0E0*I_ERI_S_D2x_Pz_S_C1000_bc-2.0E0*1*I_ERI_S_S_Pz_S_C1000_c;
  abcd[520] = 4.0E0*I_ERI_S_Dxy_Pz_S_C1000_bc;
  abcd[524] = 4.0E0*I_ERI_S_Dxz_Pz_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[517] = 4.0E0*I_ERI_Px_D2x_Pz_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Pz_S_C1001_c;
  abcd[518] = 4.0E0*I_ERI_Py_D2x_Pz_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Pz_S_C1001_c;
  abcd[519] = 4.0E0*I_ERI_Pz_D2x_Pz_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1001_c;
  abcd[521] = 4.0E0*I_ERI_Px_Dxy_Pz_S_C1001_bc;
  abcd[522] = 4.0E0*I_ERI_Py_Dxy_Pz_S_C1001_bc;
  abcd[523] = 4.0E0*I_ERI_Pz_Dxy_Pz_S_C1001_bc;
  abcd[525] = 4.0E0*I_ERI_Px_Dxz_Pz_S_C1001_bc;
  abcd[526] = 4.0E0*I_ERI_Py_Dxz_Pz_S_C1001_bc;
  abcd[527] = 4.0E0*I_ERI_Pz_Dxz_Pz_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[528] = 4.0E0*I_ERI_S_Py_Px_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[529] = 4.0E0*I_ERI_Px_Py_Px_S_C1_bc;
  abcd[530] = 4.0E0*I_ERI_Py_Py_Px_S_C1_bc;
  abcd[531] = 4.0E0*I_ERI_Pz_Py_Px_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[532] = 4.0E0*I_ERI_S_Dxy_Px_S_C1000_bc;
  abcd[536] = 4.0E0*I_ERI_S_D2y_Px_S_C1000_bc-2.0E0*1*I_ERI_S_S_Px_S_C1000_c;
  abcd[540] = 4.0E0*I_ERI_S_Dyz_Px_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[533] = 4.0E0*I_ERI_Px_Dxy_Px_S_C1001_bc;
  abcd[534] = 4.0E0*I_ERI_Py_Dxy_Px_S_C1001_bc;
  abcd[535] = 4.0E0*I_ERI_Pz_Dxy_Px_S_C1001_bc;
  abcd[537] = 4.0E0*I_ERI_Px_D2y_Px_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Px_S_C1001_c;
  abcd[538] = 4.0E0*I_ERI_Py_D2y_Px_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Px_S_C1001_c;
  abcd[539] = 4.0E0*I_ERI_Pz_D2y_Px_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Px_S_C1001_c;
  abcd[541] = 4.0E0*I_ERI_Px_Dyz_Px_S_C1001_bc;
  abcd[542] = 4.0E0*I_ERI_Py_Dyz_Px_S_C1001_bc;
  abcd[543] = 4.0E0*I_ERI_Pz_Dyz_Px_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[544] = 4.0E0*I_ERI_S_Py_Py_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[545] = 4.0E0*I_ERI_Px_Py_Py_S_C1_bc;
  abcd[546] = 4.0E0*I_ERI_Py_Py_Py_S_C1_bc;
  abcd[547] = 4.0E0*I_ERI_Pz_Py_Py_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[548] = 4.0E0*I_ERI_S_Dxy_Py_S_C1000_bc;
  abcd[552] = 4.0E0*I_ERI_S_D2y_Py_S_C1000_bc-2.0E0*1*I_ERI_S_S_Py_S_C1000_c;
  abcd[556] = 4.0E0*I_ERI_S_Dyz_Py_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[549] = 4.0E0*I_ERI_Px_Dxy_Py_S_C1001_bc;
  abcd[550] = 4.0E0*I_ERI_Py_Dxy_Py_S_C1001_bc;
  abcd[551] = 4.0E0*I_ERI_Pz_Dxy_Py_S_C1001_bc;
  abcd[553] = 4.0E0*I_ERI_Px_D2y_Py_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Py_S_C1001_c;
  abcd[554] = 4.0E0*I_ERI_Py_D2y_Py_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Py_S_C1001_c;
  abcd[555] = 4.0E0*I_ERI_Pz_D2y_Py_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Py_S_C1001_c;
  abcd[557] = 4.0E0*I_ERI_Px_Dyz_Py_S_C1001_bc;
  abcd[558] = 4.0E0*I_ERI_Py_Dyz_Py_S_C1001_bc;
  abcd[559] = 4.0E0*I_ERI_Pz_Dyz_Py_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[560] = 4.0E0*I_ERI_S_Py_Pz_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[561] = 4.0E0*I_ERI_Px_Py_Pz_S_C1_bc;
  abcd[562] = 4.0E0*I_ERI_Py_Py_Pz_S_C1_bc;
  abcd[563] = 4.0E0*I_ERI_Pz_Py_Pz_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[564] = 4.0E0*I_ERI_S_Dxy_Pz_S_C1000_bc;
  abcd[568] = 4.0E0*I_ERI_S_D2y_Pz_S_C1000_bc-2.0E0*1*I_ERI_S_S_Pz_S_C1000_c;
  abcd[572] = 4.0E0*I_ERI_S_Dyz_Pz_S_C1000_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[565] = 4.0E0*I_ERI_Px_Dxy_Pz_S_C1001_bc;
  abcd[566] = 4.0E0*I_ERI_Py_Dxy_Pz_S_C1001_bc;
  abcd[567] = 4.0E0*I_ERI_Pz_Dxy_Pz_S_C1001_bc;
  abcd[569] = 4.0E0*I_ERI_Px_D2y_Pz_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Pz_S_C1001_c;
  abcd[570] = 4.0E0*I_ERI_Py_D2y_Pz_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Pz_S_C1001_c;
  abcd[571] = 4.0E0*I_ERI_Pz_D2y_Pz_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1001_c;
  abcd[573] = 4.0E0*I_ERI_Px_Dyz_Pz_S_C1001_bc;
  abcd[574] = 4.0E0*I_ERI_Py_Dyz_Pz_S_C1001_bc;
  abcd[575] = 4.0E0*I_ERI_Pz_Dyz_Pz_S_C1001_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[576] = 4.0E0*I_ERI_S_Pz_Px_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[577] = 4.0E0*I_ERI_Px_Pz_Px_S_C1_bc;
  abcd[578] = 4.0E0*I_ERI_Py_Pz_Px_S_C1_bc;
  abcd[579] = 4.0E0*I_ERI_Pz_Pz_Px_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[580] = 4.0E0*I_ERI_S_Dxz_Px_S_C1000_bc;
  abcd[584] = 4.0E0*I_ERI_S_Dyz_Px_S_C1000_bc;
  abcd[588] = 4.0E0*I_ERI_S_D2z_Px_S_C1000_bc-2.0E0*1*I_ERI_S_S_Px_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[581] = 4.0E0*I_ERI_Px_Dxz_Px_S_C1001_bc;
  abcd[582] = 4.0E0*I_ERI_Py_Dxz_Px_S_C1001_bc;
  abcd[583] = 4.0E0*I_ERI_Pz_Dxz_Px_S_C1001_bc;
  abcd[585] = 4.0E0*I_ERI_Px_Dyz_Px_S_C1001_bc;
  abcd[586] = 4.0E0*I_ERI_Py_Dyz_Px_S_C1001_bc;
  abcd[587] = 4.0E0*I_ERI_Pz_Dyz_Px_S_C1001_bc;
  abcd[589] = 4.0E0*I_ERI_Px_D2z_Px_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Px_S_C1001_c;
  abcd[590] = 4.0E0*I_ERI_Py_D2z_Px_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Px_S_C1001_c;
  abcd[591] = 4.0E0*I_ERI_Pz_D2z_Px_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Px_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[592] = 4.0E0*I_ERI_S_Pz_Py_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[593] = 4.0E0*I_ERI_Px_Pz_Py_S_C1_bc;
  abcd[594] = 4.0E0*I_ERI_Py_Pz_Py_S_C1_bc;
  abcd[595] = 4.0E0*I_ERI_Pz_Pz_Py_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[596] = 4.0E0*I_ERI_S_Dxz_Py_S_C1000_bc;
  abcd[600] = 4.0E0*I_ERI_S_Dyz_Py_S_C1000_bc;
  abcd[604] = 4.0E0*I_ERI_S_D2z_Py_S_C1000_bc-2.0E0*1*I_ERI_S_S_Py_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[597] = 4.0E0*I_ERI_Px_Dxz_Py_S_C1001_bc;
  abcd[598] = 4.0E0*I_ERI_Py_Dxz_Py_S_C1001_bc;
  abcd[599] = 4.0E0*I_ERI_Pz_Dxz_Py_S_C1001_bc;
  abcd[601] = 4.0E0*I_ERI_Px_Dyz_Py_S_C1001_bc;
  abcd[602] = 4.0E0*I_ERI_Py_Dyz_Py_S_C1001_bc;
  abcd[603] = 4.0E0*I_ERI_Pz_Dyz_Py_S_C1001_bc;
  abcd[605] = 4.0E0*I_ERI_Px_D2z_Py_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Py_S_C1001_c;
  abcd[606] = 4.0E0*I_ERI_Py_D2z_Py_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Py_S_C1001_c;
  abcd[607] = 4.0E0*I_ERI_Pz_D2z_Py_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Py_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C0_bc
   ************************************************************/
  abcd[608] = 4.0E0*I_ERI_S_Pz_Pz_S_C0_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1_bc
   ************************************************************/
  abcd[609] = 4.0E0*I_ERI_Px_Pz_Pz_S_C1_bc;
  abcd[610] = 4.0E0*I_ERI_Py_Pz_Pz_S_C1_bc;
  abcd[611] = 4.0E0*I_ERI_Pz_Pz_Pz_S_C1_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1000_bc
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000_c
   ************************************************************/
  abcd[612] = 4.0E0*I_ERI_S_Dxz_Pz_S_C1000_bc;
  abcd[616] = 4.0E0*I_ERI_S_Dyz_Pz_S_C1000_bc;
  abcd[620] = 4.0E0*I_ERI_S_D2z_Pz_S_C1000_bc-2.0E0*1*I_ERI_S_S_Pz_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001_bc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  abcd[613] = 4.0E0*I_ERI_Px_Dxz_Pz_S_C1001_bc;
  abcd[614] = 4.0E0*I_ERI_Py_Dxz_Pz_S_C1001_bc;
  abcd[615] = 4.0E0*I_ERI_Pz_Dxz_Pz_S_C1001_bc;
  abcd[617] = 4.0E0*I_ERI_Px_Dyz_Pz_S_C1001_bc;
  abcd[618] = 4.0E0*I_ERI_Py_Dyz_Pz_S_C1001_bc;
  abcd[619] = 4.0E0*I_ERI_Pz_Dyz_Pz_S_C1001_bc;
  abcd[621] = 4.0E0*I_ERI_Px_D2z_Pz_S_C1001_bc-2.0E0*1*I_ERI_Px_S_Pz_S_C1001_c;
  abcd[622] = 4.0E0*I_ERI_Py_D2z_Pz_S_C1001_bc-2.0E0*1*I_ERI_Py_S_Pz_S_C1001_c;
  abcd[623] = 4.0E0*I_ERI_Pz_D2z_Pz_S_C1001_bc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[624] = 4.0E0*I_ERI_S_S_D2x_S_C0_cc-2.0E0*1*I_ERI_S_S_S_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[625] = 4.0E0*I_ERI_Px_S_D2x_S_C1_cc-2.0E0*1*I_ERI_Px_S_S_S_C1_c;
  abcd[626] = 4.0E0*I_ERI_Py_S_D2x_S_C1_cc-2.0E0*1*I_ERI_Py_S_S_S_C1_c;
  abcd[627] = 4.0E0*I_ERI_Pz_S_D2x_S_C1_cc-2.0E0*1*I_ERI_Pz_S_S_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[628] = 4.0E0*I_ERI_S_Px_D2x_S_C1000_cc-2.0E0*1*I_ERI_S_Px_S_S_C1000_c;
  abcd[632] = 4.0E0*I_ERI_S_Py_D2x_S_C1000_cc-2.0E0*1*I_ERI_S_Py_S_S_C1000_c;
  abcd[636] = 4.0E0*I_ERI_S_Pz_D2x_S_C1000_cc-2.0E0*1*I_ERI_S_Pz_S_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[629] = 4.0E0*I_ERI_Px_Px_D2x_S_C1001_cc-2.0E0*1*I_ERI_Px_Px_S_S_C1001_c;
  abcd[630] = 4.0E0*I_ERI_Py_Px_D2x_S_C1001_cc-2.0E0*1*I_ERI_Py_Px_S_S_C1001_c;
  abcd[631] = 4.0E0*I_ERI_Pz_Px_D2x_S_C1001_cc-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_c;
  abcd[633] = 4.0E0*I_ERI_Px_Py_D2x_S_C1001_cc-2.0E0*1*I_ERI_Px_Py_S_S_C1001_c;
  abcd[634] = 4.0E0*I_ERI_Py_Py_D2x_S_C1001_cc-2.0E0*1*I_ERI_Py_Py_S_S_C1001_c;
  abcd[635] = 4.0E0*I_ERI_Pz_Py_D2x_S_C1001_cc-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_c;
  abcd[637] = 4.0E0*I_ERI_Px_Pz_D2x_S_C1001_cc-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_c;
  abcd[638] = 4.0E0*I_ERI_Py_Pz_D2x_S_C1001_cc-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_c;
  abcd[639] = 4.0E0*I_ERI_Pz_Pz_D2x_S_C1001_cc-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[640] = 4.0E0*I_ERI_S_S_Dxy_S_C0_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[641] = 4.0E0*I_ERI_Px_S_Dxy_S_C1_cc;
  abcd[642] = 4.0E0*I_ERI_Py_S_Dxy_S_C1_cc;
  abcd[643] = 4.0E0*I_ERI_Pz_S_Dxy_S_C1_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[644] = 4.0E0*I_ERI_S_Px_Dxy_S_C1000_cc;
  abcd[648] = 4.0E0*I_ERI_S_Py_Dxy_S_C1000_cc;
  abcd[652] = 4.0E0*I_ERI_S_Pz_Dxy_S_C1000_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[645] = 4.0E0*I_ERI_Px_Px_Dxy_S_C1001_cc;
  abcd[646] = 4.0E0*I_ERI_Py_Px_Dxy_S_C1001_cc;
  abcd[647] = 4.0E0*I_ERI_Pz_Px_Dxy_S_C1001_cc;
  abcd[649] = 4.0E0*I_ERI_Px_Py_Dxy_S_C1001_cc;
  abcd[650] = 4.0E0*I_ERI_Py_Py_Dxy_S_C1001_cc;
  abcd[651] = 4.0E0*I_ERI_Pz_Py_Dxy_S_C1001_cc;
  abcd[653] = 4.0E0*I_ERI_Px_Pz_Dxy_S_C1001_cc;
  abcd[654] = 4.0E0*I_ERI_Py_Pz_Dxy_S_C1001_cc;
  abcd[655] = 4.0E0*I_ERI_Pz_Pz_Dxy_S_C1001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[656] = 4.0E0*I_ERI_S_S_Dxz_S_C0_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[657] = 4.0E0*I_ERI_Px_S_Dxz_S_C1_cc;
  abcd[658] = 4.0E0*I_ERI_Py_S_Dxz_S_C1_cc;
  abcd[659] = 4.0E0*I_ERI_Pz_S_Dxz_S_C1_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[660] = 4.0E0*I_ERI_S_Px_Dxz_S_C1000_cc;
  abcd[664] = 4.0E0*I_ERI_S_Py_Dxz_S_C1000_cc;
  abcd[668] = 4.0E0*I_ERI_S_Pz_Dxz_S_C1000_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[661] = 4.0E0*I_ERI_Px_Px_Dxz_S_C1001_cc;
  abcd[662] = 4.0E0*I_ERI_Py_Px_Dxz_S_C1001_cc;
  abcd[663] = 4.0E0*I_ERI_Pz_Px_Dxz_S_C1001_cc;
  abcd[665] = 4.0E0*I_ERI_Px_Py_Dxz_S_C1001_cc;
  abcd[666] = 4.0E0*I_ERI_Py_Py_Dxz_S_C1001_cc;
  abcd[667] = 4.0E0*I_ERI_Pz_Py_Dxz_S_C1001_cc;
  abcd[669] = 4.0E0*I_ERI_Px_Pz_Dxz_S_C1001_cc;
  abcd[670] = 4.0E0*I_ERI_Py_Pz_Dxz_S_C1001_cc;
  abcd[671] = 4.0E0*I_ERI_Pz_Pz_Dxz_S_C1001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[672] = 4.0E0*I_ERI_S_S_D2y_S_C0_cc-2.0E0*1*I_ERI_S_S_S_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[673] = 4.0E0*I_ERI_Px_S_D2y_S_C1_cc-2.0E0*1*I_ERI_Px_S_S_S_C1_c;
  abcd[674] = 4.0E0*I_ERI_Py_S_D2y_S_C1_cc-2.0E0*1*I_ERI_Py_S_S_S_C1_c;
  abcd[675] = 4.0E0*I_ERI_Pz_S_D2y_S_C1_cc-2.0E0*1*I_ERI_Pz_S_S_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[676] = 4.0E0*I_ERI_S_Px_D2y_S_C1000_cc-2.0E0*1*I_ERI_S_Px_S_S_C1000_c;
  abcd[680] = 4.0E0*I_ERI_S_Py_D2y_S_C1000_cc-2.0E0*1*I_ERI_S_Py_S_S_C1000_c;
  abcd[684] = 4.0E0*I_ERI_S_Pz_D2y_S_C1000_cc-2.0E0*1*I_ERI_S_Pz_S_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[677] = 4.0E0*I_ERI_Px_Px_D2y_S_C1001_cc-2.0E0*1*I_ERI_Px_Px_S_S_C1001_c;
  abcd[678] = 4.0E0*I_ERI_Py_Px_D2y_S_C1001_cc-2.0E0*1*I_ERI_Py_Px_S_S_C1001_c;
  abcd[679] = 4.0E0*I_ERI_Pz_Px_D2y_S_C1001_cc-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_c;
  abcd[681] = 4.0E0*I_ERI_Px_Py_D2y_S_C1001_cc-2.0E0*1*I_ERI_Px_Py_S_S_C1001_c;
  abcd[682] = 4.0E0*I_ERI_Py_Py_D2y_S_C1001_cc-2.0E0*1*I_ERI_Py_Py_S_S_C1001_c;
  abcd[683] = 4.0E0*I_ERI_Pz_Py_D2y_S_C1001_cc-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_c;
  abcd[685] = 4.0E0*I_ERI_Px_Pz_D2y_S_C1001_cc-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_c;
  abcd[686] = 4.0E0*I_ERI_Py_Pz_D2y_S_C1001_cc-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_c;
  abcd[687] = 4.0E0*I_ERI_Pz_Pz_D2y_S_C1001_cc-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[688] = 4.0E0*I_ERI_S_S_Dyz_S_C0_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[689] = 4.0E0*I_ERI_Px_S_Dyz_S_C1_cc;
  abcd[690] = 4.0E0*I_ERI_Py_S_Dyz_S_C1_cc;
  abcd[691] = 4.0E0*I_ERI_Pz_S_Dyz_S_C1_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[692] = 4.0E0*I_ERI_S_Px_Dyz_S_C1000_cc;
  abcd[696] = 4.0E0*I_ERI_S_Py_Dyz_S_C1000_cc;
  abcd[700] = 4.0E0*I_ERI_S_Pz_Dyz_S_C1000_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[693] = 4.0E0*I_ERI_Px_Px_Dyz_S_C1001_cc;
  abcd[694] = 4.0E0*I_ERI_Py_Px_Dyz_S_C1001_cc;
  abcd[695] = 4.0E0*I_ERI_Pz_Px_Dyz_S_C1001_cc;
  abcd[697] = 4.0E0*I_ERI_Px_Py_Dyz_S_C1001_cc;
  abcd[698] = 4.0E0*I_ERI_Py_Py_Dyz_S_C1001_cc;
  abcd[699] = 4.0E0*I_ERI_Pz_Py_Dyz_S_C1001_cc;
  abcd[701] = 4.0E0*I_ERI_Px_Pz_Dyz_S_C1001_cc;
  abcd[702] = 4.0E0*I_ERI_Py_Pz_Dyz_S_C1001_cc;
  abcd[703] = 4.0E0*I_ERI_Pz_Pz_Dyz_S_C1001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C0_cc
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C0_c
   ************************************************************/
  abcd[704] = 4.0E0*I_ERI_S_S_D2z_S_C0_cc-2.0E0*1*I_ERI_S_S_S_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1_cc
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_c
   ************************************************************/
  abcd[705] = 4.0E0*I_ERI_Px_S_D2z_S_C1_cc-2.0E0*1*I_ERI_Px_S_S_S_C1_c;
  abcd[706] = 4.0E0*I_ERI_Py_S_D2z_S_C1_cc-2.0E0*1*I_ERI_Py_S_S_S_C1_c;
  abcd[707] = 4.0E0*I_ERI_Pz_S_D2z_S_C1_cc-2.0E0*1*I_ERI_Pz_S_S_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000_cc
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1000_c
   ************************************************************/
  abcd[708] = 4.0E0*I_ERI_S_Px_D2z_S_C1000_cc-2.0E0*1*I_ERI_S_Px_S_S_C1000_c;
  abcd[712] = 4.0E0*I_ERI_S_Py_D2z_S_C1000_cc-2.0E0*1*I_ERI_S_Py_S_S_C1000_c;
  abcd[716] = 4.0E0*I_ERI_S_Pz_D2z_S_C1000_cc-2.0E0*1*I_ERI_S_Pz_S_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001_cc
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_c
   ************************************************************/
  abcd[709] = 4.0E0*I_ERI_Px_Px_D2z_S_C1001_cc-2.0E0*1*I_ERI_Px_Px_S_S_C1001_c;
  abcd[710] = 4.0E0*I_ERI_Py_Px_D2z_S_C1001_cc-2.0E0*1*I_ERI_Py_Px_S_S_C1001_c;
  abcd[711] = 4.0E0*I_ERI_Pz_Px_D2z_S_C1001_cc-2.0E0*1*I_ERI_Pz_Px_S_S_C1001_c;
  abcd[713] = 4.0E0*I_ERI_Px_Py_D2z_S_C1001_cc-2.0E0*1*I_ERI_Px_Py_S_S_C1001_c;
  abcd[714] = 4.0E0*I_ERI_Py_Py_D2z_S_C1001_cc-2.0E0*1*I_ERI_Py_Py_S_S_C1001_c;
  abcd[715] = 4.0E0*I_ERI_Pz_Py_D2z_S_C1001_cc-2.0E0*1*I_ERI_Pz_Py_S_S_C1001_c;
  abcd[717] = 4.0E0*I_ERI_Px_Pz_D2z_S_C1001_cc-2.0E0*1*I_ERI_Px_Pz_S_S_C1001_c;
  abcd[718] = 4.0E0*I_ERI_Py_Pz_D2z_S_C1001_cc-2.0E0*1*I_ERI_Py_Pz_S_S_C1001_c;
  abcd[719] = 4.0E0*I_ERI_Pz_Pz_D2z_S_C1001_cc-2.0E0*1*I_ERI_Pz_Pz_S_S_C1001_c;
}
