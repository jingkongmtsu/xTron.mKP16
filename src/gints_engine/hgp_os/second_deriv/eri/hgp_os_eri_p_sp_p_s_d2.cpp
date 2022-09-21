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
// BRA1 as redundant position, total RHS integrals evaluated as: 156018
// BRA2 as redundant position, total RHS integrals evaluated as: 143688
// KET1 as redundant position, total RHS integrals evaluated as: 162408
// KET2 as redundant position, total RHS integrals evaluated as: 143992
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
// BRA1  KET2
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
// KET1  KET2
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
// KET2  KET2
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_eri_p_sp_p_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000001_aa = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1000001_ac = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1000001_a = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1000001_c = 0.0E0;
  Double I_ERI_S_S_S_S_C1000001 = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_S_Px_S_S_C1001001 = 0.0E0;
  Double I_ERI_S_Py_S_S_C1001001 = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1000001_cc = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001_aa = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_ac = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001001_a = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1001001_cc = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000001_ad = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000001_cd = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000001_cd = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_ad = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001_cd = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_dd = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_dd = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_dd = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_dd = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_dd = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_dd = 0.0E0;

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
      Double QDX   = QX - D[0];
      Double QDY   = QY - D[1];
      Double QDZ   = QZ - D[2];
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
       * shell quartet name: SQ_ERI_S_S_P_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_S_S_Px_S_M3_vrr = QCX*I_ERI_S_S_S_S_M3_vrr+WQX*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Py_S_M3_vrr = QCY*I_ERI_S_S_S_S_M3_vrr+WQY*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Pz_S_M3_vrr = QCZ*I_ERI_S_S_S_S_M3_vrr+WQZ*I_ERI_S_S_S_S_M4_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M3_vrr = PAX*I_ERI_Px_S_S_S_M3_vrr+WPX*I_ERI_Px_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_S_S_M3_vrr = PAY*I_ERI_Px_S_S_S_M3_vrr+WPY*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_Dxz_S_S_S_M3_vrr = PAZ*I_ERI_Px_S_S_S_M3_vrr+WPZ*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_S_S_M3_vrr = PAY*I_ERI_Py_S_S_S_M3_vrr+WPY*I_ERI_Py_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Dyz_S_S_S_M3_vrr = PAZ*I_ERI_Py_S_S_S_M3_vrr+WPZ*I_ERI_Py_S_S_S_M4_vrr;
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
       * shell quartet name: SQ_ERI_P_S_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M2_vrr = PAX*I_ERI_S_S_Px_S_M2_vrr+WPX*I_ERI_S_S_Px_S_M3_vrr+oned2k*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Py_S_Px_S_M2_vrr = PAY*I_ERI_S_S_Px_S_M2_vrr+WPY*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Pz_S_Px_S_M2_vrr = PAZ*I_ERI_S_S_Px_S_M2_vrr+WPZ*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Px_S_Py_S_M2_vrr = PAX*I_ERI_S_S_Py_S_M2_vrr+WPX*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Py_S_Py_S_M2_vrr = PAY*I_ERI_S_S_Py_S_M2_vrr+WPY*I_ERI_S_S_Py_S_M3_vrr+oned2k*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Pz_S_Py_S_M2_vrr = PAZ*I_ERI_S_S_Py_S_M2_vrr+WPZ*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Px_S_Pz_S_M2_vrr = PAX*I_ERI_S_S_Pz_S_M2_vrr+WPX*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Py_S_Pz_S_M2_vrr = PAY*I_ERI_S_S_Pz_S_M2_vrr+WPY*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Pz_S_Pz_S_M2_vrr = PAZ*I_ERI_S_S_Pz_S_M2_vrr+WPZ*I_ERI_S_S_Pz_S_M3_vrr+oned2k*I_ERI_S_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M2_vrr = QCX*I_ERI_D2x_S_S_S_M2_vrr+WQX*I_ERI_D2x_S_S_S_M3_vrr+2*oned2k*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_Px_S_M2_vrr = QCX*I_ERI_Dxy_S_S_S_M2_vrr+WQX*I_ERI_Dxy_S_S_S_M3_vrr+oned2k*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_Px_S_M2_vrr = QCX*I_ERI_Dxz_S_S_S_M2_vrr+WQX*I_ERI_Dxz_S_S_S_M3_vrr+oned2k*I_ERI_Pz_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_Px_S_M2_vrr = QCX*I_ERI_D2y_S_S_S_M2_vrr+WQX*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_Px_S_M2_vrr = QCX*I_ERI_Dyz_S_S_S_M2_vrr+WQX*I_ERI_Dyz_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_Px_S_M2_vrr = QCX*I_ERI_D2z_S_S_S_M2_vrr+WQX*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_D2x_S_Py_S_M2_vrr = QCY*I_ERI_D2x_S_S_S_M2_vrr+WQY*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_Py_S_M2_vrr = QCY*I_ERI_Dxy_S_S_S_M2_vrr+WQY*I_ERI_Dxy_S_S_S_M3_vrr+oned2k*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_Py_S_M2_vrr = QCY*I_ERI_Dxz_S_S_S_M2_vrr+WQY*I_ERI_Dxz_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_Py_S_M2_vrr = QCY*I_ERI_D2y_S_S_S_M2_vrr+WQY*I_ERI_D2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_Py_S_M2_vrr = QCY*I_ERI_Dyz_S_S_S_M2_vrr+WQY*I_ERI_Dyz_S_S_S_M3_vrr+oned2k*I_ERI_Pz_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_Py_S_M2_vrr = QCY*I_ERI_D2z_S_S_S_M2_vrr+WQY*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_D2x_S_Pz_S_M2_vrr = QCZ*I_ERI_D2x_S_S_S_M2_vrr+WQZ*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_Pz_S_M2_vrr = QCZ*I_ERI_Dxy_S_S_S_M2_vrr+WQZ*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_Pz_S_M2_vrr = QCZ*I_ERI_Dxz_S_S_S_M2_vrr+WQZ*I_ERI_Dxz_S_S_S_M3_vrr+oned2k*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_Pz_S_M2_vrr = QCZ*I_ERI_D2y_S_S_S_M2_vrr+WQZ*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_Pz_S_M2_vrr = QCZ*I_ERI_Dyz_S_S_S_M2_vrr+WQZ*I_ERI_Dyz_S_S_S_M3_vrr+oned2k*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_Pz_S_M2_vrr = QCZ*I_ERI_D2z_S_S_S_M2_vrr+WQZ*I_ERI_D2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Pz_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_P_M1
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_S_Px_M1_vrr = QDX*I_ERI_S_S_S_S_M1_vrr+WQX*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_S_Py_M1_vrr = QDY*I_ERI_S_S_S_S_M1_vrr+WQY*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_S_Pz_M1_vrr = QDZ*I_ERI_S_S_S_S_M1_vrr+WQZ*I_ERI_S_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M1_vrr = QCX*I_ERI_S_S_Px_S_M1_vrr+WQX*I_ERI_S_S_Px_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dxy_S_M1_vrr = QCY*I_ERI_S_S_Px_S_M1_vrr+WQY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_Dxz_S_M1_vrr = QCZ*I_ERI_S_S_Px_S_M1_vrr+WQZ*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_D2y_S_M1_vrr = QCY*I_ERI_S_S_Py_S_M1_vrr+WQY*I_ERI_S_S_Py_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dyz_S_M1_vrr = QCZ*I_ERI_S_S_Py_S_M1_vrr+WQZ*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_S_D2z_S_M1_vrr = QCZ*I_ERI_S_S_Pz_S_M1_vrr+WQZ*I_ERI_S_S_Pz_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_P_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M1_vrr = QCX*I_ERI_Px_S_Px_S_M1_vrr+WQX*I_ERI_Px_S_Px_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_D2x_S_M1_vrr = QCX*I_ERI_Py_S_Px_S_M1_vrr+WQX*I_ERI_Py_S_Px_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_D2x_S_M1_vrr = QCX*I_ERI_Pz_S_Px_S_M1_vrr+WQX*I_ERI_Pz_S_Px_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_Px_S_Dxy_S_M1_vrr = QCY*I_ERI_Px_S_Px_S_M1_vrr+WQY*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_Dxy_S_M1_vrr = QCY*I_ERI_Py_S_Px_S_M1_vrr+WQY*I_ERI_Py_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxy_S_M1_vrr = QCY*I_ERI_Pz_S_Px_S_M1_vrr+WQY*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_Dxz_S_M1_vrr = QCZ*I_ERI_Px_S_Px_S_M1_vrr+WQZ*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_Dxz_S_M1_vrr = QCZ*I_ERI_Py_S_Px_S_M1_vrr+WQZ*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxz_S_M1_vrr = QCZ*I_ERI_Pz_S_Px_S_M1_vrr+WQZ*I_ERI_Pz_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_D2y_S_M1_vrr = QCY*I_ERI_Px_S_Py_S_M1_vrr+WQY*I_ERI_Px_S_Py_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Py_S_D2y_S_M1_vrr = QCY*I_ERI_Py_S_Py_S_M1_vrr+WQY*I_ERI_Py_S_Py_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_D2y_S_M1_vrr = QCY*I_ERI_Pz_S_Py_S_M1_vrr+WQY*I_ERI_Pz_S_Py_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_Px_S_Dyz_S_M1_vrr = QCZ*I_ERI_Px_S_Py_S_M1_vrr+WQZ*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Dyz_S_M1_vrr = QCZ*I_ERI_Py_S_Py_S_M1_vrr+WQZ*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_Dyz_S_M1_vrr = QCZ*I_ERI_Pz_S_Py_S_M1_vrr+WQZ*I_ERI_Pz_S_Py_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_D2z_S_M1_vrr = QCZ*I_ERI_Px_S_Pz_S_M1_vrr+WQZ*I_ERI_Px_S_Pz_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Py_S_D2z_S_M1_vrr = QCZ*I_ERI_Py_S_Pz_S_M1_vrr+WQZ*I_ERI_Py_S_Pz_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_D2z_S_M1_vrr = QCZ*I_ERI_Pz_S_Pz_S_M1_vrr+WQZ*I_ERI_Pz_S_Pz_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M1_vrr = PAX*I_ERI_Px_S_Px_S_M1_vrr+WPX*I_ERI_Px_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Px_S_M1_vrr = PAY*I_ERI_Px_S_Px_S_M1_vrr+WPY*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxz_S_Px_S_M1_vrr = PAZ*I_ERI_Px_S_Px_S_M1_vrr+WPZ*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Px_S_M1_vrr = PAY*I_ERI_Py_S_Px_S_M1_vrr+WPY*I_ERI_Py_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Dyz_S_Px_S_M1_vrr = PAZ*I_ERI_Py_S_Px_S_M1_vrr+WPZ*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Px_S_M1_vrr = PAZ*I_ERI_Pz_S_Px_S_M1_vrr+WPZ*I_ERI_Pz_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_Py_S_M1_vrr = PAX*I_ERI_Px_S_Py_S_M1_vrr+WPX*I_ERI_Px_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Py_S_M1_vrr = PAY*I_ERI_Px_S_Py_S_M1_vrr+WPY*I_ERI_Px_S_Py_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_Py_S_M1_vrr = PAZ*I_ERI_Px_S_Py_S_M1_vrr+WPZ*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_D2y_S_Py_S_M1_vrr = PAY*I_ERI_Py_S_Py_S_M1_vrr+WPY*I_ERI_Py_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_Py_S_M1_vrr = PAZ*I_ERI_Py_S_Py_S_M1_vrr+WPZ*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_Py_S_M1_vrr = PAZ*I_ERI_Pz_S_Py_S_M1_vrr+WPZ*I_ERI_Pz_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_Pz_S_M1_vrr = PAX*I_ERI_Px_S_Pz_S_M1_vrr+WPX*I_ERI_Px_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_Pz_S_M1_vrr = PAY*I_ERI_Px_S_Pz_S_M1_vrr+WPY*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_Dxz_S_Pz_S_M1_vrr = PAZ*I_ERI_Px_S_Pz_S_M1_vrr+WPZ*I_ERI_Px_S_Pz_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Pz_S_M1_vrr = PAY*I_ERI_Py_S_Pz_S_M1_vrr+WPY*I_ERI_Py_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Dyz_S_Pz_S_M1_vrr = PAZ*I_ERI_Py_S_Pz_S_M1_vrr+WPZ*I_ERI_Py_S_Pz_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Pz_S_M1_vrr = PAZ*I_ERI_Pz_S_Pz_S_M1_vrr+WPZ*I_ERI_Pz_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M1_vrr = QCX*I_ERI_D2x_S_Px_S_M1_vrr+WQX*I_ERI_D2x_S_Px_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr+2*oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxy_S_D2x_S_M1_vrr = QCX*I_ERI_Dxy_S_Px_S_M1_vrr+WQX*I_ERI_Dxy_S_Px_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_Dxz_S_D2x_S_M1_vrr = QCX*I_ERI_Dxz_S_Px_S_M1_vrr+WQX*I_ERI_Dxz_S_Px_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_D2x_S_M1_vrr = QCX*I_ERI_D2y_S_Px_S_M1_vrr+WQX*I_ERI_D2y_S_Px_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_D2x_S_M1_vrr = QCX*I_ERI_Dyz_S_Px_S_M1_vrr+WQX*I_ERI_Dyz_S_Px_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_D2x_S_M1_vrr = QCX*I_ERI_D2z_S_Px_S_M1_vrr+WQX*I_ERI_D2z_S_Px_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Dxy_S_M1_vrr = QCY*I_ERI_D2x_S_Px_S_M1_vrr+WQY*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxy_S_M1_vrr = QCY*I_ERI_Dxy_S_Px_S_M1_vrr+WQY*I_ERI_Dxy_S_Px_S_M2_vrr+oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxz_S_Dxy_S_M1_vrr = QCY*I_ERI_Dxz_S_Px_S_M1_vrr+WQY*I_ERI_Dxz_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Dxy_S_M1_vrr = QCY*I_ERI_D2y_S_Px_S_M1_vrr+WQY*I_ERI_D2y_S_Px_S_M2_vrr+2*oned2k*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_Dyz_S_Dxy_S_M1_vrr = QCY*I_ERI_Dyz_S_Px_S_M1_vrr+WQY*I_ERI_Dyz_S_Px_S_M2_vrr+oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Dxy_S_M1_vrr = QCY*I_ERI_D2z_S_Px_S_M1_vrr+WQY*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2x_S_Px_S_M1_vrr+WQZ*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxz_S_M1_vrr = QCZ*I_ERI_Dxy_S_Px_S_M1_vrr+WQZ*I_ERI_Dxy_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2y_S_Px_S_M1_vrr+WQZ*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2z_S_Px_S_M1_vrr+WQZ*I_ERI_D2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_D2y_S_M1_vrr = QCY*I_ERI_D2x_S_Py_S_M1_vrr+WQY*I_ERI_D2x_S_Py_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = QCY*I_ERI_Dxy_S_Py_S_M1_vrr+WQY*I_ERI_Dxy_S_Py_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = QCY*I_ERI_Dxz_S_Py_S_M1_vrr+WQY*I_ERI_Dxz_S_Py_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = QCY*I_ERI_D2y_S_Py_S_M1_vrr+WQY*I_ERI_D2y_S_Py_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = QCY*I_ERI_Dyz_S_Py_S_M1_vrr+WQY*I_ERI_Dyz_S_Py_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = QCY*I_ERI_D2z_S_Py_S_M1_vrr+WQY*I_ERI_D2z_S_Py_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Dyz_S_M1_vrr = QCZ*I_ERI_D2x_S_Py_S_M1_vrr+WQZ*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Dyz_S_M1_vrr = QCZ*I_ERI_Dxy_S_Py_S_M1_vrr+WQZ*I_ERI_Dxy_S_Py_S_M2_vrr;
      Double I_ERI_D2y_S_Dyz_S_M1_vrr = QCZ*I_ERI_D2y_S_Py_S_M1_vrr+WQZ*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_Dyz_S_M1_vrr = QCZ*I_ERI_D2z_S_Py_S_M1_vrr+WQZ*I_ERI_D2z_S_Py_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_D2z_S_M1_vrr = QCZ*I_ERI_D2x_S_Pz_S_M1_vrr+WQZ*I_ERI_D2x_S_Pz_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxy_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_D2z_S_M1_vrr = QCZ*I_ERI_D2y_S_Pz_S_M1_vrr+WQZ*I_ERI_D2y_S_Pz_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dyz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_D2z_S_M1_vrr = QCZ*I_ERI_D2z_S_Pz_S_M1_vrr+WQZ*I_ERI_D2z_S_Pz_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M1_vrr = PAX*I_ERI_D2x_S_Px_S_M1_vrr+WPX*I_ERI_D2x_S_Px_S_M2_vrr+2*oned2z*I_ERI_Px_S_Px_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Px_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_Px_S_M1_vrr = PAY*I_ERI_D2x_S_Px_S_M1_vrr+WPY*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_Px_S_M1_vrr = PAZ*I_ERI_D2x_S_Px_S_M1_vrr+WPZ*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_Px_S_M1_vrr = PAX*I_ERI_D2y_S_Px_S_M1_vrr+WPX*I_ERI_D2y_S_Px_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_Px_S_M1_vrr = PAX*I_ERI_D2z_S_Px_S_M1_vrr+WPX*I_ERI_D2z_S_Px_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_Px_S_M1_vrr = PAY*I_ERI_D2y_S_Px_S_M1_vrr+WPY*I_ERI_D2y_S_Px_S_M2_vrr+2*oned2z*I_ERI_Py_S_Px_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Px_S_M1_vrr = PAZ*I_ERI_D2y_S_Px_S_M1_vrr+WPZ*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Px_S_M1_vrr = PAZ*I_ERI_D2z_S_Px_S_M1_vrr+WPZ*I_ERI_D2z_S_Px_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Px_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_Py_S_M1_vrr = PAX*I_ERI_D2x_S_Py_S_M1_vrr+WPX*I_ERI_D2x_S_Py_S_M2_vrr+2*oned2z*I_ERI_Px_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_Py_S_M1_vrr = PAY*I_ERI_D2x_S_Py_S_M1_vrr+WPY*I_ERI_D2x_S_Py_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_Py_S_M1_vrr = PAZ*I_ERI_D2x_S_Py_S_M1_vrr+WPZ*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_Py_S_M1_vrr = PAX*I_ERI_D2y_S_Py_S_M1_vrr+WPX*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_Py_S_M1_vrr = PAX*I_ERI_D2z_S_Py_S_M1_vrr+WPX*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_Py_S_M1_vrr = PAY*I_ERI_D2y_S_Py_S_M1_vrr+WPY*I_ERI_D2y_S_Py_S_M2_vrr+2*oned2z*I_ERI_Py_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_Py_S_M1_vrr = PAZ*I_ERI_D2y_S_Py_S_M1_vrr+WPZ*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_Py_S_M1_vrr = PAZ*I_ERI_D2z_S_Py_S_M1_vrr+WPZ*I_ERI_D2z_S_Py_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_F3x_S_Pz_S_M1_vrr = PAX*I_ERI_D2x_S_Pz_S_M1_vrr+WPX*I_ERI_D2x_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_F2xy_S_Pz_S_M1_vrr = PAY*I_ERI_D2x_S_Pz_S_M1_vrr+WPY*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xz_S_Pz_S_M1_vrr = PAZ*I_ERI_D2x_S_Pz_S_M1_vrr+WPZ*I_ERI_D2x_S_Pz_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M1_vrr = PAX*I_ERI_D2y_S_Pz_S_M1_vrr+WPX*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M1_vrr = PAX*I_ERI_D2z_S_Pz_S_M1_vrr+WPX*I_ERI_D2z_S_Pz_S_M2_vrr;
      Double I_ERI_F3y_S_Pz_S_M1_vrr = PAY*I_ERI_D2y_S_Pz_S_M1_vrr+WPY*I_ERI_D2y_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_F2yz_S_Pz_S_M1_vrr = PAZ*I_ERI_D2y_S_Pz_S_M1_vrr+WPZ*I_ERI_D2y_S_Pz_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_Pz_S_M1_vrr = PAZ*I_ERI_D2z_S_Pz_S_M1_vrr+WPZ*I_ERI_D2z_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_P
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_S_Px_vrr = QDX*I_ERI_S_S_S_S_vrr+WQX*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_S_Py_vrr = QDY*I_ERI_S_S_S_S_vrr+WQY*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_S_Pz_vrr = QDZ*I_ERI_S_S_S_S_vrr+WQZ*I_ERI_S_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_P
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_P
       * RHS shell quartet name: SQ_ERI_S_S_S_P_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_Px_S_S_Px_vrr = PAX*I_ERI_S_S_S_Px_vrr+WPX*I_ERI_S_S_S_Px_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Py_S_S_Px_vrr = PAY*I_ERI_S_S_S_Px_vrr+WPY*I_ERI_S_S_S_Px_M1_vrr;
      Double I_ERI_Pz_S_S_Px_vrr = PAZ*I_ERI_S_S_S_Px_vrr+WPZ*I_ERI_S_S_S_Px_M1_vrr;
      Double I_ERI_Px_S_S_Py_vrr = PAX*I_ERI_S_S_S_Py_vrr+WPX*I_ERI_S_S_S_Py_M1_vrr;
      Double I_ERI_Py_S_S_Py_vrr = PAY*I_ERI_S_S_S_Py_vrr+WPY*I_ERI_S_S_S_Py_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Pz_S_S_Py_vrr = PAZ*I_ERI_S_S_S_Py_vrr+WPZ*I_ERI_S_S_S_Py_M1_vrr;
      Double I_ERI_Px_S_S_Pz_vrr = PAX*I_ERI_S_S_S_Pz_vrr+WPX*I_ERI_S_S_S_Pz_M1_vrr;
      Double I_ERI_Py_S_S_Pz_vrr = PAY*I_ERI_S_S_S_Pz_vrr+WPY*I_ERI_S_S_S_Pz_M1_vrr;
      Double I_ERI_Pz_S_S_Pz_vrr = PAZ*I_ERI_S_S_S_Pz_vrr+WPZ*I_ERI_S_S_S_Pz_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;

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
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_Px_D2x_S_vrr = PBX*I_ERI_S_S_D2x_S_vrr+WPX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Py_D2x_S_vrr = PBY*I_ERI_S_S_D2x_S_vrr+WPY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Pz_D2x_S_vrr = PBZ*I_ERI_S_S_D2x_S_vrr+WPZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Px_Dxy_S_vrr = PBX*I_ERI_S_S_Dxy_S_vrr+WPX*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Py_Dxy_S_vrr = PBY*I_ERI_S_S_Dxy_S_vrr+WPY*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Pz_Dxy_S_vrr = PBZ*I_ERI_S_S_Dxy_S_vrr+WPZ*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_S_Px_Dxz_S_vrr = PBX*I_ERI_S_S_Dxz_S_vrr+WPX*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Py_Dxz_S_vrr = PBY*I_ERI_S_S_Dxz_S_vrr+WPY*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_S_Pz_Dxz_S_vrr = PBZ*I_ERI_S_S_Dxz_S_vrr+WPZ*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_Px_D2y_S_vrr = PBX*I_ERI_S_S_D2y_S_vrr+WPX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Py_D2y_S_vrr = PBY*I_ERI_S_S_D2y_S_vrr+WPY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Pz_D2y_S_vrr = PBZ*I_ERI_S_S_D2y_S_vrr+WPZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Px_Dyz_S_vrr = PBX*I_ERI_S_S_Dyz_S_vrr+WPX*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_S_Py_Dyz_S_vrr = PBY*I_ERI_S_S_Dyz_S_vrr+WPY*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_S_Pz_Dyz_S_vrr = PBZ*I_ERI_S_S_Dyz_S_vrr+WPZ*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_Px_D2z_S_vrr = PBX*I_ERI_S_S_D2z_S_vrr+WPX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Py_D2z_S_vrr = PBY*I_ERI_S_S_D2z_S_vrr+WPY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Pz_D2z_S_vrr = PBZ*I_ERI_S_S_D2z_S_vrr+WPZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_vrr = PAX*I_ERI_S_S_D2x_S_vrr+WPX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_D2x_S_vrr = PAY*I_ERI_S_S_D2x_S_vrr+WPY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_D2x_S_vrr = PAZ*I_ERI_S_S_D2x_S_vrr+WPZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_Dxy_S_vrr = PAX*I_ERI_S_S_Dxy_S_vrr+WPX*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_Dxy_S_vrr = PAY*I_ERI_S_S_Dxy_S_vrr+WPY*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Dxy_S_vrr = PAZ*I_ERI_S_S_Dxy_S_vrr+WPZ*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Px_S_Dxz_S_vrr = PAX*I_ERI_S_S_Dxz_S_vrr+WPX*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Py_S_Dxz_S_vrr = PAY*I_ERI_S_S_Dxz_S_vrr+WPY*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Pz_S_Dxz_S_vrr = PAZ*I_ERI_S_S_Dxz_S_vrr+WPZ*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_D2y_S_vrr = PAX*I_ERI_S_S_D2y_S_vrr+WPX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_D2y_S_vrr = PAY*I_ERI_S_S_D2y_S_vrr+WPY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Pz_S_D2y_S_vrr = PAZ*I_ERI_S_S_D2y_S_vrr+WPZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Dyz_S_vrr = PAX*I_ERI_S_S_Dyz_S_vrr+WPX*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Py_S_Dyz_S_vrr = PAY*I_ERI_S_S_Dyz_S_vrr+WPY*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Pz_S_Dyz_S_vrr = PAZ*I_ERI_S_S_Dyz_S_vrr+WPZ*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_D2z_S_vrr = PAX*I_ERI_S_S_D2z_S_vrr+WPX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_D2z_S_vrr = PAY*I_ERI_S_S_D2z_S_vrr+WPY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_D2z_S_vrr = PAZ*I_ERI_S_S_D2z_S_vrr+WPZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_Px_vrr = QDX*I_ERI_D2x_S_S_S_vrr+WQX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_Px_vrr = QDX*I_ERI_Dxy_S_S_S_vrr+WQX*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_S_Px_vrr = QDX*I_ERI_Dxz_S_S_S_vrr+WQX*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_Px_vrr = QDX*I_ERI_D2y_S_S_S_vrr+WQX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_S_Px_vrr = QDX*I_ERI_Dyz_S_S_S_vrr+WQX*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_Px_vrr = QDX*I_ERI_D2z_S_S_S_vrr+WQX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_S_Py_vrr = QDY*I_ERI_D2x_S_S_S_vrr+WQY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_Py_vrr = QDY*I_ERI_Dxy_S_S_S_vrr+WQY*I_ERI_Dxy_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_S_Py_vrr = QDY*I_ERI_Dxz_S_S_S_vrr+WQY*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_Py_vrr = QDY*I_ERI_D2y_S_S_S_vrr+WQY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_S_Py_vrr = QDY*I_ERI_Dyz_S_S_S_vrr+WQY*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Pz_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_Py_vrr = QDY*I_ERI_D2z_S_S_S_vrr+WQY*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_S_Pz_vrr = QDZ*I_ERI_D2x_S_S_S_vrr+WQZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_Pz_vrr = QDZ*I_ERI_Dxy_S_S_S_vrr+WQZ*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Dxz_S_S_Pz_vrr = QDZ*I_ERI_Dxz_S_S_S_vrr+WQZ*I_ERI_Dxz_S_S_S_M1_vrr+oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_Pz_vrr = QDZ*I_ERI_D2y_S_S_S_vrr+WQZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Dyz_S_S_Pz_vrr = QDZ*I_ERI_Dyz_S_S_S_vrr+WQZ*I_ERI_Dyz_S_S_S_M1_vrr+oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_Pz_vrr = QDZ*I_ERI_D2z_S_S_S_vrr+WQZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_S_S_M1_vrr;

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
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_vrr = PAX*I_ERI_Px_S_D2x_S_vrr+WPX*I_ERI_Px_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_D2x_S_vrr = PAY*I_ERI_Px_S_D2x_S_vrr+WPY*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_D2x_S_vrr = PAZ*I_ERI_Px_S_D2x_S_vrr+WPZ*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_D2x_S_vrr = PAY*I_ERI_Py_S_D2x_S_vrr+WPY*I_ERI_Py_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_D2x_S_vrr = PAZ*I_ERI_Py_S_D2x_S_vrr+WPZ*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_D2x_S_vrr = PAZ*I_ERI_Pz_S_D2x_S_vrr+WPZ*I_ERI_Pz_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_Dxy_S_vrr = PAX*I_ERI_Px_S_Dxy_S_vrr+WPX*I_ERI_Px_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxy_S_vrr = PAY*I_ERI_Px_S_Dxy_S_vrr+WPY*I_ERI_Px_S_Dxy_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxy_S_vrr = PAZ*I_ERI_Px_S_Dxy_S_vrr+WPZ*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_D2y_S_Dxy_S_vrr = PAY*I_ERI_Py_S_Dxy_S_vrr+WPY*I_ERI_Py_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxy_S_vrr = PAZ*I_ERI_Py_S_Dxy_S_vrr+WPZ*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_D2z_S_Dxy_S_vrr = PAZ*I_ERI_Pz_S_Dxy_S_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_D2x_S_Dxz_S_vrr = PAX*I_ERI_Px_S_Dxz_S_vrr+WPX*I_ERI_Px_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxz_S_vrr = PAY*I_ERI_Px_S_Dxz_S_vrr+WPY*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxz_S_vrr = PAZ*I_ERI_Px_S_Dxz_S_vrr+WPZ*I_ERI_Px_S_Dxz_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxz_S_vrr = PAY*I_ERI_Py_S_Dxz_S_vrr+WPY*I_ERI_Py_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxz_S_vrr = PAZ*I_ERI_Py_S_Dxz_S_vrr+WPZ*I_ERI_Py_S_Dxz_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxz_S_vrr = PAZ*I_ERI_Pz_S_Dxz_S_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_D2y_S_vrr = PAX*I_ERI_Px_S_D2y_S_vrr+WPX*I_ERI_Px_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_D2y_S_vrr = PAY*I_ERI_Px_S_D2y_S_vrr+WPY*I_ERI_Px_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_D2y_S_vrr = PAZ*I_ERI_Px_S_D2y_S_vrr+WPZ*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = PAY*I_ERI_Py_S_D2y_S_vrr+WPY*I_ERI_Py_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_D2y_S_vrr = PAZ*I_ERI_Py_S_D2y_S_vrr+WPZ*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = PAZ*I_ERI_Pz_S_D2y_S_vrr+WPZ*I_ERI_Pz_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Dyz_S_vrr = PAX*I_ERI_Px_S_Dyz_S_vrr+WPX*I_ERI_Px_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxy_S_Dyz_S_vrr = PAY*I_ERI_Px_S_Dyz_S_vrr+WPY*I_ERI_Px_S_Dyz_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Dxz_S_Dyz_S_vrr = PAZ*I_ERI_Px_S_Dyz_S_vrr+WPZ*I_ERI_Px_S_Dyz_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Dyz_S_vrr = PAY*I_ERI_Py_S_Dyz_S_vrr+WPY*I_ERI_Py_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_Dyz_S_Dyz_S_vrr = PAZ*I_ERI_Py_S_Dyz_S_vrr+WPZ*I_ERI_Py_S_Dyz_S_M1_vrr+oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Dyz_S_vrr = PAZ*I_ERI_Pz_S_Dyz_S_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = PAX*I_ERI_Px_S_D2z_S_vrr+WPX*I_ERI_Px_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_D2z_S_vrr = PAY*I_ERI_Px_S_D2z_S_vrr+WPY*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_D2z_S_vrr = PAZ*I_ERI_Px_S_D2z_S_vrr+WPZ*I_ERI_Px_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = PAY*I_ERI_Py_S_D2z_S_vrr+WPY*I_ERI_Py_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_D2z_S_vrr = PAZ*I_ERI_Py_S_D2z_S_vrr+WPZ*I_ERI_Py_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = PAZ*I_ERI_Pz_S_D2z_S_vrr+WPZ*I_ERI_Pz_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_vrr = QCX*I_ERI_Px_S_D2x_S_vrr+WQX*I_ERI_Px_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Px_S_Px_S_vrr-2*rhod2esq*I_ERI_Px_S_Px_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Py_S_F3x_S_vrr = QCX*I_ERI_Py_S_D2x_S_vrr+WQX*I_ERI_Py_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Py_S_Px_S_vrr-2*rhod2esq*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_F3x_S_vrr = QCX*I_ERI_Pz_S_D2x_S_vrr+WQX*I_ERI_Pz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Pz_S_Px_S_vrr-2*rhod2esq*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_F2xy_S_vrr = QCY*I_ERI_Px_S_D2x_S_vrr+WQY*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Py_S_F2xy_S_vrr = QCY*I_ERI_Py_S_D2x_S_vrr+WQY*I_ERI_Py_S_D2x_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_F2xy_S_vrr = QCY*I_ERI_Pz_S_D2x_S_vrr+WQY*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_F2xz_S_vrr = QCZ*I_ERI_Px_S_D2x_S_vrr+WQZ*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Py_S_F2xz_S_vrr = QCZ*I_ERI_Py_S_D2x_S_vrr+WQZ*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_F2xz_S_vrr = QCZ*I_ERI_Pz_S_D2x_S_vrr+WQZ*I_ERI_Pz_S_D2x_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_Fx2y_S_vrr = QCX*I_ERI_Px_S_D2y_S_vrr+WQX*I_ERI_Px_S_D2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_Fx2y_S_vrr = QCX*I_ERI_Py_S_D2y_S_vrr+WQX*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2y_S_vrr = QCX*I_ERI_Pz_S_D2y_S_vrr+WQX*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Fxyz_S_vrr = QCZ*I_ERI_Px_S_Dxy_S_vrr+WQZ*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_Py_S_Fxyz_S_vrr = QCZ*I_ERI_Py_S_Dxy_S_vrr+WQZ*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_Pz_S_Fxyz_S_vrr = QCZ*I_ERI_Pz_S_Dxy_S_vrr+WQZ*I_ERI_Pz_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Px_S_Fx2z_S_vrr = QCX*I_ERI_Px_S_D2z_S_vrr+WQX*I_ERI_Px_S_D2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_Fx2z_S_vrr = QCX*I_ERI_Py_S_D2z_S_vrr+WQX*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2z_S_vrr = QCX*I_ERI_Pz_S_D2z_S_vrr+WQX*I_ERI_Pz_S_D2z_S_M1_vrr;
      Double I_ERI_Px_S_F3y_S_vrr = QCY*I_ERI_Px_S_D2y_S_vrr+WQY*I_ERI_Px_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Px_S_Py_S_vrr-2*rhod2esq*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_F3y_S_vrr = QCY*I_ERI_Py_S_D2y_S_vrr+WQY*I_ERI_Py_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Py_S_Py_S_vrr-2*rhod2esq*I_ERI_Py_S_Py_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Pz_S_F3y_S_vrr = QCY*I_ERI_Pz_S_D2y_S_vrr+WQY*I_ERI_Pz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Pz_S_Py_S_vrr-2*rhod2esq*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_F2yz_S_vrr = QCZ*I_ERI_Px_S_D2y_S_vrr+WQZ*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_F2yz_S_vrr = QCZ*I_ERI_Py_S_D2y_S_vrr+WQZ*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Pz_S_F2yz_S_vrr = QCZ*I_ERI_Pz_S_D2y_S_vrr+WQZ*I_ERI_Pz_S_D2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Fy2z_S_vrr = QCY*I_ERI_Px_S_D2z_S_vrr+WQY*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_Fy2z_S_vrr = QCY*I_ERI_Py_S_D2z_S_vrr+WQY*I_ERI_Py_S_D2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fy2z_S_vrr = QCY*I_ERI_Pz_S_D2z_S_vrr+WQY*I_ERI_Pz_S_D2z_S_M1_vrr;
      Double I_ERI_Px_S_F3z_S_vrr = QCZ*I_ERI_Px_S_D2z_S_vrr+WQZ*I_ERI_Px_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Px_S_Pz_S_vrr-2*rhod2esq*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Py_S_F3z_S_vrr = QCZ*I_ERI_Py_S_D2z_S_vrr+WQZ*I_ERI_Py_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Py_S_Pz_S_vrr-2*rhod2esq*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_Pz_S_F3z_S_vrr = QCZ*I_ERI_Pz_S_D2z_S_vrr+WQZ*I_ERI_Pz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Pz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Pz_S_Pz_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_vrr = QCX*I_ERI_D2x_S_D2x_S_vrr+WQX*I_ERI_D2x_S_D2x_S_M1_vrr+2*oned2e*I_ERI_D2x_S_Px_S_vrr-2*rhod2esq*I_ERI_D2x_S_Px_S_M1_vrr+2*oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_F3x_S_vrr = QCX*I_ERI_Dxy_S_D2x_S_vrr+WQX*I_ERI_Dxy_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Dxy_S_Px_S_vrr-2*rhod2esq*I_ERI_Dxy_S_Px_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_F3x_S_vrr = QCX*I_ERI_Dxz_S_D2x_S_vrr+WQX*I_ERI_Dxz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Dxz_S_Px_S_vrr-2*rhod2esq*I_ERI_Dxz_S_Px_S_M1_vrr+oned2k*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_F3x_S_vrr = QCX*I_ERI_D2y_S_D2x_S_vrr+WQX*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_D2y_S_Px_S_vrr-2*rhod2esq*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_F3x_S_vrr = QCX*I_ERI_Dyz_S_D2x_S_vrr+WQX*I_ERI_Dyz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Dyz_S_Px_S_vrr-2*rhod2esq*I_ERI_Dyz_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_F3x_S_vrr = QCX*I_ERI_D2z_S_D2x_S_vrr+WQX*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_D2z_S_Px_S_vrr-2*rhod2esq*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_F2xy_S_vrr = QCY*I_ERI_D2x_S_D2x_S_vrr+WQY*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xy_S_vrr = QCY*I_ERI_Dxy_S_D2x_S_vrr+WQY*I_ERI_Dxy_S_D2x_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xy_S_vrr = QCY*I_ERI_Dxz_S_D2x_S_vrr+WQY*I_ERI_Dxz_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_F2xy_S_vrr = QCY*I_ERI_D2y_S_D2x_S_vrr+WQY*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xy_S_vrr = QCY*I_ERI_Dyz_S_D2x_S_vrr+WQY*I_ERI_Dyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_F2xy_S_vrr = QCY*I_ERI_D2z_S_D2x_S_vrr+WQY*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_F2xz_S_vrr = QCZ*I_ERI_D2x_S_D2x_S_vrr+WQZ*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xz_S_vrr = QCZ*I_ERI_Dxy_S_D2x_S_vrr+WQZ*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xz_S_vrr = QCZ*I_ERI_Dxz_S_D2x_S_vrr+WQZ*I_ERI_Dxz_S_D2x_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_F2xz_S_vrr = QCZ*I_ERI_D2y_S_D2x_S_vrr+WQZ*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xz_S_vrr = QCZ*I_ERI_Dyz_S_D2x_S_vrr+WQZ*I_ERI_Dyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_F2xz_S_vrr = QCZ*I_ERI_D2z_S_D2x_S_vrr+WQZ*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2y_S_vrr = QCX*I_ERI_D2x_S_D2y_S_vrr+WQX*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_vrr = QCX*I_ERI_Dxy_S_D2y_S_vrr+WQX*I_ERI_Dxy_S_D2y_S_M1_vrr+oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2y_S_vrr = QCX*I_ERI_Dxz_S_D2y_S_vrr+WQX*I_ERI_Dxz_S_D2y_S_M1_vrr+oned2k*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2y_S_vrr = QCX*I_ERI_D2y_S_D2y_S_vrr+WQX*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2y_S_vrr = QCX*I_ERI_Dyz_S_D2y_S_vrr+WQX*I_ERI_Dyz_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2y_S_vrr = QCX*I_ERI_D2z_S_D2y_S_vrr+WQX*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fxyz_S_vrr = QCZ*I_ERI_D2x_S_Dxy_S_vrr+WQZ*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxy_S_Fxyz_S_vrr = QCZ*I_ERI_Dxy_S_Dxy_S_vrr+WQZ*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxz_S_Fxyz_S_vrr = QCZ*I_ERI_Dxz_S_Dxy_S_vrr+WQZ*I_ERI_Dxz_S_Dxy_S_M1_vrr+oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_D2y_S_Fxyz_S_vrr = QCZ*I_ERI_D2y_S_Dxy_S_vrr+WQZ*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Dyz_S_Fxyz_S_vrr = QCZ*I_ERI_Dyz_S_Dxy_S_vrr+WQZ*I_ERI_Dyz_S_Dxy_S_M1_vrr+oned2k*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_D2z_S_Fxyz_S_vrr = QCZ*I_ERI_D2z_S_Dxy_S_vrr+WQZ*I_ERI_D2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Dxy_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2z_S_vrr = QCX*I_ERI_D2x_S_D2z_S_vrr+WQX*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_vrr = QCX*I_ERI_Dxy_S_D2z_S_vrr+WQX*I_ERI_Dxy_S_D2z_S_M1_vrr+oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2z_S_vrr = QCX*I_ERI_Dxz_S_D2z_S_vrr+WQX*I_ERI_Dxz_S_D2z_S_M1_vrr+oned2k*I_ERI_Pz_S_D2z_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2z_S_vrr = QCX*I_ERI_D2y_S_D2z_S_vrr+WQX*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_vrr = QCX*I_ERI_Dyz_S_D2z_S_vrr+WQX*I_ERI_Dyz_S_D2z_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2z_S_vrr = QCX*I_ERI_D2z_S_D2z_S_vrr+WQX*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_D2x_S_F3y_S_vrr = QCY*I_ERI_D2x_S_D2y_S_vrr+WQY*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2e*I_ERI_D2x_S_Py_S_vrr-2*rhod2esq*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_F3y_S_vrr = QCY*I_ERI_Dxy_S_D2y_S_vrr+WQY*I_ERI_Dxy_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Dxy_S_Py_S_vrr-2*rhod2esq*I_ERI_Dxy_S_Py_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_F3y_S_vrr = QCY*I_ERI_Dxz_S_D2y_S_vrr+WQY*I_ERI_Dxz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Dxz_S_Py_S_vrr-2*rhod2esq*I_ERI_Dxz_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_F3y_S_vrr = QCY*I_ERI_D2y_S_D2y_S_vrr+WQY*I_ERI_D2y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_D2y_S_Py_S_vrr-2*rhod2esq*I_ERI_D2y_S_Py_S_M1_vrr+2*oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_F3y_S_vrr = QCY*I_ERI_Dyz_S_D2y_S_vrr+WQY*I_ERI_Dyz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Dyz_S_Py_S_vrr-2*rhod2esq*I_ERI_Dyz_S_Py_S_M1_vrr+oned2k*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_F3y_S_vrr = QCY*I_ERI_D2z_S_D2y_S_vrr+WQY*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_D2z_S_Py_S_vrr-2*rhod2esq*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_F2yz_S_vrr = QCZ*I_ERI_D2x_S_D2y_S_vrr+WQZ*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_F2yz_S_vrr = QCZ*I_ERI_Dxy_S_D2y_S_vrr+WQZ*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_F2yz_S_vrr = QCZ*I_ERI_Dxz_S_D2y_S_vrr+WQZ*I_ERI_Dxz_S_D2y_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_F2yz_S_vrr = QCZ*I_ERI_D2y_S_D2y_S_vrr+WQZ*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_F2yz_S_vrr = QCZ*I_ERI_Dyz_S_D2y_S_vrr+WQZ*I_ERI_Dyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_F2yz_S_vrr = QCZ*I_ERI_D2z_S_D2y_S_vrr+WQZ*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fy2z_S_vrr = QCY*I_ERI_D2x_S_D2z_S_vrr+WQY*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fy2z_S_vrr = QCY*I_ERI_Dxy_S_D2z_S_vrr+WQY*I_ERI_Dxy_S_D2z_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fy2z_S_vrr = QCY*I_ERI_Dxz_S_D2z_S_vrr+WQY*I_ERI_Dxz_S_D2z_S_M1_vrr;
      Double I_ERI_D2y_S_Fy2z_S_vrr = QCY*I_ERI_D2y_S_D2z_S_vrr+WQY*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fy2z_S_vrr = QCY*I_ERI_Dyz_S_D2z_S_vrr+WQY*I_ERI_Dyz_S_D2z_S_M1_vrr+oned2k*I_ERI_Pz_S_D2z_S_M1_vrr;
      Double I_ERI_D2z_S_Fy2z_S_vrr = QCY*I_ERI_D2z_S_D2z_S_vrr+WQY*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_D2x_S_F3z_S_vrr = QCZ*I_ERI_D2x_S_D2z_S_vrr+WQZ*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2e*I_ERI_D2x_S_Pz_S_vrr-2*rhod2esq*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_Dxy_S_F3z_S_vrr = QCZ*I_ERI_Dxy_S_D2z_S_vrr+WQZ*I_ERI_Dxy_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Dxy_S_Pz_S_vrr-2*rhod2esq*I_ERI_Dxy_S_Pz_S_M1_vrr;
      Double I_ERI_Dxz_S_F3z_S_vrr = QCZ*I_ERI_Dxz_S_D2z_S_vrr+WQZ*I_ERI_Dxz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Dxz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Dxz_S_Pz_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_D2y_S_F3z_S_vrr = QCZ*I_ERI_D2y_S_D2z_S_vrr+WQZ*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_D2y_S_Pz_S_vrr-2*rhod2esq*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Dyz_S_F3z_S_vrr = QCZ*I_ERI_Dyz_S_D2z_S_vrr+WQZ*I_ERI_Dyz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Dyz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Dyz_S_Pz_S_M1_vrr+oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_D2z_S_F3z_S_vrr = QCZ*I_ERI_D2z_S_D2z_S_vrr+WQZ*I_ERI_D2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_D2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_D2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_Pz_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_vrr = PAX*I_ERI_D2x_S_D2x_S_vrr+WPX*I_ERI_D2x_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xy_S_D2x_S_vrr = PAY*I_ERI_D2x_S_D2x_S_vrr+WPY*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_D2x_S_vrr = PAZ*I_ERI_D2x_S_D2x_S_vrr+WPZ*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2x_S_vrr = PAX*I_ERI_D2y_S_D2x_S_vrr+WPX*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2x_S_vrr = PAZ*I_ERI_Dxy_S_D2x_S_vrr+WPZ*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2x_S_vrr = PAX*I_ERI_D2z_S_D2x_S_vrr+WPX*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3y_S_D2x_S_vrr = PAY*I_ERI_D2y_S_D2x_S_vrr+WPY*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_D2x_S_vrr = PAZ*I_ERI_D2y_S_D2x_S_vrr+WPZ*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2x_S_vrr = PAY*I_ERI_D2z_S_D2x_S_vrr+WPY*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_D2x_S_vrr = PAZ*I_ERI_D2z_S_D2x_S_vrr+WPZ*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_Dxy_S_vrr = PAX*I_ERI_D2x_S_Dxy_S_vrr+WPX*I_ERI_D2x_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxy_S_vrr = PAY*I_ERI_D2x_S_Dxy_S_vrr+WPY*I_ERI_D2x_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxy_S_vrr = PAZ*I_ERI_D2x_S_Dxy_S_vrr+WPZ*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_vrr = PAX*I_ERI_D2y_S_Dxy_S_vrr+WPX*I_ERI_D2y_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxy_S_vrr = PAZ*I_ERI_Dxy_S_Dxy_S_vrr+WPZ*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_vrr = PAX*I_ERI_D2z_S_Dxy_S_vrr+WPX*I_ERI_D2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_Dxy_S_vrr = PAY*I_ERI_D2y_S_Dxy_S_vrr+WPY*I_ERI_D2y_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxy_S_vrr = PAZ*I_ERI_D2y_S_Dxy_S_vrr+WPZ*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxy_S_vrr = PAY*I_ERI_D2z_S_Dxy_S_vrr+WPY*I_ERI_D2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Dxy_S_vrr = PAZ*I_ERI_D2z_S_Dxy_S_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M1_vrr;
      Double I_ERI_F3x_S_Dxz_S_vrr = PAX*I_ERI_D2x_S_Dxz_S_vrr+WPX*I_ERI_D2x_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxz_S_vrr = PAY*I_ERI_D2x_S_Dxz_S_vrr+WPY*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxz_S_vrr = PAZ*I_ERI_D2x_S_Dxz_S_vrr+WPZ*I_ERI_D2x_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_vrr = PAX*I_ERI_D2y_S_Dxz_S_vrr+WPX*I_ERI_D2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxz_S_vrr = PAZ*I_ERI_Dxy_S_Dxz_S_vrr+WPZ*I_ERI_Dxy_S_Dxz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_vrr = PAX*I_ERI_D2z_S_Dxz_S_vrr+WPX*I_ERI_D2z_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3y_S_Dxz_S_vrr = PAY*I_ERI_D2y_S_Dxz_S_vrr+WPY*I_ERI_D2y_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxz_S_vrr = PAZ*I_ERI_D2y_S_Dxz_S_vrr+WPZ*I_ERI_D2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxz_S_vrr = PAY*I_ERI_D2z_S_Dxz_S_vrr+WPY*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3z_S_Dxz_S_vrr = PAZ*I_ERI_D2z_S_Dxz_S_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_D2y_S_vrr = PAX*I_ERI_D2x_S_D2y_S_vrr+WPX*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_D2y_S_vrr = PAY*I_ERI_D2x_S_D2y_S_vrr+WPY*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_F2xz_S_D2y_S_vrr = PAZ*I_ERI_D2x_S_D2y_S_vrr+WPZ*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2y_S_vrr = PAX*I_ERI_D2y_S_D2y_S_vrr+WPX*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2y_S_vrr = PAZ*I_ERI_Dxy_S_D2y_S_vrr+WPZ*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2y_S_vrr = PAX*I_ERI_D2z_S_D2y_S_vrr+WPX*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_D2y_S_vrr = PAY*I_ERI_D2y_S_D2y_S_vrr+WPY*I_ERI_D2y_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_F2yz_S_D2y_S_vrr = PAZ*I_ERI_D2y_S_D2y_S_vrr+WPZ*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2y_S_vrr = PAY*I_ERI_D2z_S_D2y_S_vrr+WPY*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3z_S_D2y_S_vrr = PAZ*I_ERI_D2z_S_D2y_S_vrr+WPZ*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Dyz_S_vrr = PAX*I_ERI_D2x_S_Dyz_S_vrr+WPX*I_ERI_D2x_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xy_S_Dyz_S_vrr = PAY*I_ERI_D2x_S_Dyz_S_vrr+WPY*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xz_S_Dyz_S_vrr = PAZ*I_ERI_D2x_S_Dyz_S_vrr+WPZ*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_vrr = PAX*I_ERI_D2y_S_Dyz_S_vrr+WPX*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dyz_S_vrr = PAZ*I_ERI_Dxy_S_Dyz_S_vrr+WPZ*I_ERI_Dxy_S_Dyz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_vrr = PAX*I_ERI_D2z_S_Dyz_S_vrr+WPX*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3y_S_Dyz_S_vrr = PAY*I_ERI_D2y_S_Dyz_S_vrr+WPY*I_ERI_D2y_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_F2yz_S_Dyz_S_vrr = PAZ*I_ERI_D2y_S_Dyz_S_vrr+WPZ*I_ERI_D2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dyz_S_vrr = PAY*I_ERI_D2z_S_Dyz_S_vrr+WPY*I_ERI_D2z_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3z_S_Dyz_S_vrr = PAZ*I_ERI_D2z_S_Dyz_S_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_D2z_S_vrr = PAX*I_ERI_D2x_S_D2z_S_vrr+WPX*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_D2z_S_vrr = PAY*I_ERI_D2x_S_D2z_S_vrr+WPY*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_D2z_S_vrr = PAZ*I_ERI_D2x_S_D2z_S_vrr+WPZ*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2z_S_vrr = PAX*I_ERI_D2y_S_D2z_S_vrr+WPX*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2z_S_vrr = PAZ*I_ERI_Dxy_S_D2z_S_vrr+WPZ*I_ERI_Dxy_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2z_S_vrr = PAX*I_ERI_D2z_S_D2z_S_vrr+WPX*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_D2z_S_vrr = PAY*I_ERI_D2y_S_D2z_S_vrr+WPY*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_D2z_S_vrr = PAZ*I_ERI_D2y_S_D2z_S_vrr+WPZ*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2z_S_vrr = PAY*I_ERI_D2z_S_D2z_S_vrr+WPY*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_D2z_S_vrr = PAZ*I_ERI_D2z_S_D2z_S_vrr+WPZ*I_ERI_D2z_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000001_aa_coefs = ic2*jc2*alpha*alpha;
      I_ERI_F3x_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000001_aa += SQ_ERI_F_S_P_S_C1000001_aa_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_a_coefs = ic2*jc2*alpha;
      I_ERI_Px_S_Px_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001_a += SQ_ERI_P_S_P_S_C1000001_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1000001_ac_coefs = ic2*jc2*alpha*gamma;
      I_ERI_D2x_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1000001_ac += SQ_ERI_D_S_D_S_C1000001_ac_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1000001_a_coefs = ic2*jc2*alpha;
      I_ERI_D2x_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1000001_a += SQ_ERI_D_S_S_S_C1000001_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000001_c_coefs = ic2*jc2*gamma;
      I_ERI_S_S_D2x_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1000001_c += SQ_ERI_S_S_D_S_C1000001_c_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C1000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C1000001_coefs = ic2*jc2;
      I_ERI_S_S_S_S_C1000001 += SQ_ERI_S_S_S_S_C1000001_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001001_c_coefs = ic2_1*jc2*gamma;
      I_ERI_S_Px_D2x_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001001_c += SQ_ERI_S_P_D_S_C1001001_c_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1001001_coefs = ic2_1*jc2;
      I_ERI_S_Px_S_S_C1001001 += SQ_ERI_S_P_S_S_C1001001_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1001001 += SQ_ERI_S_P_S_S_C1001001_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1001001 += SQ_ERI_S_P_S_S_C1001001_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1000001_cc_coefs = ic2*jc2*gamma*gamma;
      I_ERI_Px_S_F3x_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1000001_cc += SQ_ERI_P_S_F_S_C1000001_cc_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_c_coefs = ic2*jc2*gamma;
      I_ERI_Px_S_Px_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001_c += SQ_ERI_P_S_P_S_C1000001_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000001_d_coefs = ic2*jc2*delta;
      I_ERI_Px_S_S_Px_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000001_d += SQ_ERI_P_S_S_P_C1000001_d_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_d_coefs = ic2*jc2*delta;
      I_ERI_Px_S_Px_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001_d += SQ_ERI_P_S_P_S_C1000001_d_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1001001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1001001_aa_coefs = ic2_1*jc2*alpha*alpha;
      I_ERI_G4x_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1001001_aa += SQ_ERI_G_S_P_S_C1001001_aa_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_aa_coefs = ic2_1*jc2*alpha*alpha;
      I_ERI_F3x_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001_aa += SQ_ERI_F_S_P_S_C1001001_aa_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_a_coefs = ic2_1*jc2*alpha;
      I_ERI_D2x_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_a += SQ_ERI_D_S_P_S_C1001001_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_a_coefs = ic2_1*jc2*alpha;
      I_ERI_Px_S_Px_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_a += SQ_ERI_P_S_P_S_C1001001_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001001_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001001_ac_coefs = ic2_1*jc2*alpha*gamma;
      I_ERI_F3x_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001001_ac += SQ_ERI_F_S_D_S_C1001001_ac_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_ac_coefs = ic2_1*jc2*alpha*gamma;
      I_ERI_D2x_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_ac += SQ_ERI_D_S_D_S_C1001001_ac_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001001_a_coefs = ic2_1*jc2*alpha;
      I_ERI_F3x_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001001_a += SQ_ERI_F_S_S_S_C1001001_a_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001001_a_coefs = ic2_1*jc2*alpha;
      I_ERI_D2x_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001001_a += SQ_ERI_D_S_S_S_C1001001_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_C1001001_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_C1001001_cc_coefs = ic2_1*jc2*gamma*gamma;
      I_ERI_D2x_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_C1001001_cc += SQ_ERI_D_S_F_S_C1001001_cc_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1001001_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1001001_cc_coefs = ic2_1*jc2*gamma*gamma;
      I_ERI_Px_S_F3x_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1001001_cc += SQ_ERI_P_S_F_S_C1001001_cc_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_c_coefs = ic2_1*jc2*gamma;
      I_ERI_D2x_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_c += SQ_ERI_D_S_P_S_C1001001_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_c_coefs = ic2_1*jc2*gamma;
      I_ERI_Px_S_Px_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_c += SQ_ERI_P_S_P_S_C1001001_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_D2x_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1001001_d += SQ_ERI_D_S_S_P_C1001001_d_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_Px_S_S_Px_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1001001_d += SQ_ERI_P_S_S_P_C1001001_d_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_D2x_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_d += SQ_ERI_D_S_P_S_C1001001_d_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_Px_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1000001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1000001_ad_coefs = ic2*jc2*alpha*delta;
      I_ERI_D2x_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1000001_ad += SQ_ERI_D_S_D_S_C1000001_ad_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001_ad_coefs = ic2*jc2*alpha*delta;
      I_ERI_D2x_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000001_ad += SQ_ERI_D_S_P_S_C1000001_ad_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000001_d_coefs = ic2*jc2*delta;
      I_ERI_S_S_D2x_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1000001_d += SQ_ERI_S_S_D_S_C1000001_d_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1000001_d_coefs = ic2*jc2*delta;
      I_ERI_S_S_Px_S_C1000001_d += SQ_ERI_S_S_P_S_C1000001_d_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000001_d += SQ_ERI_S_S_P_S_C1000001_d_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000001_d += SQ_ERI_S_S_P_S_C1000001_d_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_S_Px_D2x_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001001_d += SQ_ERI_S_P_D_S_C1001001_d_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001001_d_coefs = ic2_1*jc2*delta;
      I_ERI_S_Px_Px_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001001_d += SQ_ERI_S_P_P_S_C1001001_d_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1000001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1000001_cd_coefs = ic2*jc2*gamma*delta;
      I_ERI_Px_S_F3x_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1000001_cd += SQ_ERI_P_S_F_S_C1000001_cd_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_cd_coefs = ic2*jc2*gamma*delta;
      I_ERI_Px_S_D2x_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000001_cd += SQ_ERI_P_S_D_S_C1000001_cd_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001001_ad_coefs = ic2_1*jc2*alpha*delta;
      I_ERI_F3x_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001001_ad += SQ_ERI_F_S_D_S_C1001001_ad_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_ad_coefs = ic2_1*jc2*alpha*delta;
      I_ERI_F3x_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001_ad += SQ_ERI_F_S_P_S_C1001001_ad_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_ad_coefs = ic2_1*jc2*alpha*delta;
      I_ERI_D2x_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_ad += SQ_ERI_D_S_D_S_C1001001_ad_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_ad
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_ad_coefs = ic2_1*jc2*alpha*delta;
      I_ERI_D2x_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_ad += SQ_ERI_D_S_P_S_C1001001_ad_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_C1001001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_C1001001_cd_coefs = ic2_1*jc2*gamma*delta;
      I_ERI_D2x_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_C1001001_cd += SQ_ERI_D_S_F_S_C1001001_cd_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_cd_coefs = ic2_1*jc2*gamma*delta;
      I_ERI_D2x_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_cd += SQ_ERI_D_S_D_S_C1001001_cd_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1001001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1001001_cd_coefs = ic2_1*jc2*gamma*delta;
      I_ERI_Px_S_F3x_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1001001_cd += SQ_ERI_P_S_F_S_C1001001_cd_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001_cd_coefs = ic2_1*jc2*gamma*delta;
      I_ERI_Px_S_D2x_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001_cd += SQ_ERI_P_S_D_S_C1001001_cd_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1000001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1000001_dd_coefs = ic2*jc2*delta*delta;
      I_ERI_Px_S_F3x_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1000001_dd += SQ_ERI_P_S_F_S_C1000001_dd_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_dd_coefs = ic2*jc2*delta*delta;
      I_ERI_Px_S_D2x_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000001_dd += SQ_ERI_P_S_D_S_C1000001_dd_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_dd_coefs = ic2*jc2*delta*delta;
      I_ERI_Px_S_Px_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001_dd += SQ_ERI_P_S_P_S_C1000001_dd_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_D2x_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_C1001001_dd += SQ_ERI_D_S_F_S_C1001001_dd_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_D2x_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_dd += SQ_ERI_D_S_D_S_C1001001_dd_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_D2x_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_dd += SQ_ERI_D_S_P_S_C1001001_dd_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_Px_S_F3x_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1001001_dd += SQ_ERI_P_S_F_S_C1001001_dd_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_Px_S_D2x_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001_dd += SQ_ERI_P_S_D_S_C1001001_dd_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_dd_coefs = ic2_1*jc2*delta*delta;
      I_ERI_Px_S_Px_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_dd += SQ_ERI_P_S_P_S_C1001001_dd_coefs*I_ERI_Pz_S_Pz_S_vrr;
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
  Double CDX = C[0] - D[0];
  Double CDY = C[1] - D[1];
  Double CDZ = C[2] - D[2];

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_d
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001_d
   ************************************************************/
  Double I_ERI_S_S_Px_Px_C1000001_d = I_ERI_S_S_D2x_S_C1000001_d+CDX*I_ERI_S_S_Px_S_C1000001_d;
  Double I_ERI_S_S_Py_Px_C1000001_d = I_ERI_S_S_Dxy_S_C1000001_d+CDX*I_ERI_S_S_Py_S_C1000001_d;
  Double I_ERI_S_S_Pz_Px_C1000001_d = I_ERI_S_S_Dxz_S_C1000001_d+CDX*I_ERI_S_S_Pz_S_C1000001_d;
  Double I_ERI_S_S_Px_Py_C1000001_d = I_ERI_S_S_Dxy_S_C1000001_d+CDY*I_ERI_S_S_Px_S_C1000001_d;
  Double I_ERI_S_S_Py_Py_C1000001_d = I_ERI_S_S_D2y_S_C1000001_d+CDY*I_ERI_S_S_Py_S_C1000001_d;
  Double I_ERI_S_S_Pz_Py_C1000001_d = I_ERI_S_S_Dyz_S_C1000001_d+CDY*I_ERI_S_S_Pz_S_C1000001_d;
  Double I_ERI_S_S_Px_Pz_C1000001_d = I_ERI_S_S_Dxz_S_C1000001_d+CDZ*I_ERI_S_S_Px_S_C1000001_d;
  Double I_ERI_S_S_Py_Pz_C1000001_d = I_ERI_S_S_Dyz_S_C1000001_d+CDZ*I_ERI_S_S_Py_S_C1000001_d;
  Double I_ERI_S_S_Pz_Pz_C1000001_d = I_ERI_S_S_D2z_S_C1000001_d+CDZ*I_ERI_S_S_Pz_S_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_d
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001_d
   ************************************************************/
  Double I_ERI_S_Px_Px_Px_C1001001_d = I_ERI_S_Px_D2x_S_C1001001_d+CDX*I_ERI_S_Px_Px_S_C1001001_d;
  Double I_ERI_S_Py_Px_Px_C1001001_d = I_ERI_S_Py_D2x_S_C1001001_d+CDX*I_ERI_S_Py_Px_S_C1001001_d;
  Double I_ERI_S_Pz_Px_Px_C1001001_d = I_ERI_S_Pz_D2x_S_C1001001_d+CDX*I_ERI_S_Pz_Px_S_C1001001_d;
  Double I_ERI_S_Px_Py_Px_C1001001_d = I_ERI_S_Px_Dxy_S_C1001001_d+CDX*I_ERI_S_Px_Py_S_C1001001_d;
  Double I_ERI_S_Py_Py_Px_C1001001_d = I_ERI_S_Py_Dxy_S_C1001001_d+CDX*I_ERI_S_Py_Py_S_C1001001_d;
  Double I_ERI_S_Pz_Py_Px_C1001001_d = I_ERI_S_Pz_Dxy_S_C1001001_d+CDX*I_ERI_S_Pz_Py_S_C1001001_d;
  Double I_ERI_S_Px_Pz_Px_C1001001_d = I_ERI_S_Px_Dxz_S_C1001001_d+CDX*I_ERI_S_Px_Pz_S_C1001001_d;
  Double I_ERI_S_Py_Pz_Px_C1001001_d = I_ERI_S_Py_Dxz_S_C1001001_d+CDX*I_ERI_S_Py_Pz_S_C1001001_d;
  Double I_ERI_S_Pz_Pz_Px_C1001001_d = I_ERI_S_Pz_Dxz_S_C1001001_d+CDX*I_ERI_S_Pz_Pz_S_C1001001_d;
  Double I_ERI_S_Px_Px_Py_C1001001_d = I_ERI_S_Px_Dxy_S_C1001001_d+CDY*I_ERI_S_Px_Px_S_C1001001_d;
  Double I_ERI_S_Py_Px_Py_C1001001_d = I_ERI_S_Py_Dxy_S_C1001001_d+CDY*I_ERI_S_Py_Px_S_C1001001_d;
  Double I_ERI_S_Pz_Px_Py_C1001001_d = I_ERI_S_Pz_Dxy_S_C1001001_d+CDY*I_ERI_S_Pz_Px_S_C1001001_d;
  Double I_ERI_S_Px_Py_Py_C1001001_d = I_ERI_S_Px_D2y_S_C1001001_d+CDY*I_ERI_S_Px_Py_S_C1001001_d;
  Double I_ERI_S_Py_Py_Py_C1001001_d = I_ERI_S_Py_D2y_S_C1001001_d+CDY*I_ERI_S_Py_Py_S_C1001001_d;
  Double I_ERI_S_Pz_Py_Py_C1001001_d = I_ERI_S_Pz_D2y_S_C1001001_d+CDY*I_ERI_S_Pz_Py_S_C1001001_d;
  Double I_ERI_S_Px_Pz_Py_C1001001_d = I_ERI_S_Px_Dyz_S_C1001001_d+CDY*I_ERI_S_Px_Pz_S_C1001001_d;
  Double I_ERI_S_Py_Pz_Py_C1001001_d = I_ERI_S_Py_Dyz_S_C1001001_d+CDY*I_ERI_S_Py_Pz_S_C1001001_d;
  Double I_ERI_S_Pz_Pz_Py_C1001001_d = I_ERI_S_Pz_Dyz_S_C1001001_d+CDY*I_ERI_S_Pz_Pz_S_C1001001_d;
  Double I_ERI_S_Px_Px_Pz_C1001001_d = I_ERI_S_Px_Dxz_S_C1001001_d+CDZ*I_ERI_S_Px_Px_S_C1001001_d;
  Double I_ERI_S_Py_Px_Pz_C1001001_d = I_ERI_S_Py_Dxz_S_C1001001_d+CDZ*I_ERI_S_Py_Px_S_C1001001_d;
  Double I_ERI_S_Pz_Px_Pz_C1001001_d = I_ERI_S_Pz_Dxz_S_C1001001_d+CDZ*I_ERI_S_Pz_Px_S_C1001001_d;
  Double I_ERI_S_Px_Py_Pz_C1001001_d = I_ERI_S_Px_Dyz_S_C1001001_d+CDZ*I_ERI_S_Px_Py_S_C1001001_d;
  Double I_ERI_S_Py_Py_Pz_C1001001_d = I_ERI_S_Py_Dyz_S_C1001001_d+CDZ*I_ERI_S_Py_Py_S_C1001001_d;
  Double I_ERI_S_Pz_Py_Pz_C1001001_d = I_ERI_S_Pz_Dyz_S_C1001001_d+CDZ*I_ERI_S_Pz_Py_S_C1001001_d;
  Double I_ERI_S_Px_Pz_Pz_C1001001_d = I_ERI_S_Px_D2z_S_C1001001_d+CDZ*I_ERI_S_Px_Pz_S_C1001001_d;
  Double I_ERI_S_Py_Pz_Pz_C1001001_d = I_ERI_S_Py_D2z_S_C1001001_d+CDZ*I_ERI_S_Py_Pz_S_C1001001_d;
  Double I_ERI_S_Pz_Pz_Pz_C1001001_d = I_ERI_S_Pz_D2z_S_C1001001_d+CDZ*I_ERI_S_Pz_Pz_S_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ad
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_ad
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1000001_ad = I_ERI_D2x_S_D2x_S_C1000001_ad+CDX*I_ERI_D2x_S_Px_S_C1000001_ad;
  Double I_ERI_Dxy_S_Px_Px_C1000001_ad = I_ERI_Dxy_S_D2x_S_C1000001_ad+CDX*I_ERI_Dxy_S_Px_S_C1000001_ad;
  Double I_ERI_Dxz_S_Px_Px_C1000001_ad = I_ERI_Dxz_S_D2x_S_C1000001_ad+CDX*I_ERI_Dxz_S_Px_S_C1000001_ad;
  Double I_ERI_D2y_S_Px_Px_C1000001_ad = I_ERI_D2y_S_D2x_S_C1000001_ad+CDX*I_ERI_D2y_S_Px_S_C1000001_ad;
  Double I_ERI_Dyz_S_Px_Px_C1000001_ad = I_ERI_Dyz_S_D2x_S_C1000001_ad+CDX*I_ERI_Dyz_S_Px_S_C1000001_ad;
  Double I_ERI_D2z_S_Px_Px_C1000001_ad = I_ERI_D2z_S_D2x_S_C1000001_ad+CDX*I_ERI_D2z_S_Px_S_C1000001_ad;
  Double I_ERI_D2x_S_Py_Px_C1000001_ad = I_ERI_D2x_S_Dxy_S_C1000001_ad+CDX*I_ERI_D2x_S_Py_S_C1000001_ad;
  Double I_ERI_Dxy_S_Py_Px_C1000001_ad = I_ERI_Dxy_S_Dxy_S_C1000001_ad+CDX*I_ERI_Dxy_S_Py_S_C1000001_ad;
  Double I_ERI_Dxz_S_Py_Px_C1000001_ad = I_ERI_Dxz_S_Dxy_S_C1000001_ad+CDX*I_ERI_Dxz_S_Py_S_C1000001_ad;
  Double I_ERI_D2y_S_Py_Px_C1000001_ad = I_ERI_D2y_S_Dxy_S_C1000001_ad+CDX*I_ERI_D2y_S_Py_S_C1000001_ad;
  Double I_ERI_Dyz_S_Py_Px_C1000001_ad = I_ERI_Dyz_S_Dxy_S_C1000001_ad+CDX*I_ERI_Dyz_S_Py_S_C1000001_ad;
  Double I_ERI_D2z_S_Py_Px_C1000001_ad = I_ERI_D2z_S_Dxy_S_C1000001_ad+CDX*I_ERI_D2z_S_Py_S_C1000001_ad;
  Double I_ERI_D2x_S_Pz_Px_C1000001_ad = I_ERI_D2x_S_Dxz_S_C1000001_ad+CDX*I_ERI_D2x_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxy_S_Pz_Px_C1000001_ad = I_ERI_Dxy_S_Dxz_S_C1000001_ad+CDX*I_ERI_Dxy_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxz_S_Pz_Px_C1000001_ad = I_ERI_Dxz_S_Dxz_S_C1000001_ad+CDX*I_ERI_Dxz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2y_S_Pz_Px_C1000001_ad = I_ERI_D2y_S_Dxz_S_C1000001_ad+CDX*I_ERI_D2y_S_Pz_S_C1000001_ad;
  Double I_ERI_Dyz_S_Pz_Px_C1000001_ad = I_ERI_Dyz_S_Dxz_S_C1000001_ad+CDX*I_ERI_Dyz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2z_S_Pz_Px_C1000001_ad = I_ERI_D2z_S_Dxz_S_C1000001_ad+CDX*I_ERI_D2z_S_Pz_S_C1000001_ad;
  Double I_ERI_D2x_S_Px_Py_C1000001_ad = I_ERI_D2x_S_Dxy_S_C1000001_ad+CDY*I_ERI_D2x_S_Px_S_C1000001_ad;
  Double I_ERI_Dxy_S_Px_Py_C1000001_ad = I_ERI_Dxy_S_Dxy_S_C1000001_ad+CDY*I_ERI_Dxy_S_Px_S_C1000001_ad;
  Double I_ERI_Dxz_S_Px_Py_C1000001_ad = I_ERI_Dxz_S_Dxy_S_C1000001_ad+CDY*I_ERI_Dxz_S_Px_S_C1000001_ad;
  Double I_ERI_D2y_S_Px_Py_C1000001_ad = I_ERI_D2y_S_Dxy_S_C1000001_ad+CDY*I_ERI_D2y_S_Px_S_C1000001_ad;
  Double I_ERI_Dyz_S_Px_Py_C1000001_ad = I_ERI_Dyz_S_Dxy_S_C1000001_ad+CDY*I_ERI_Dyz_S_Px_S_C1000001_ad;
  Double I_ERI_D2z_S_Px_Py_C1000001_ad = I_ERI_D2z_S_Dxy_S_C1000001_ad+CDY*I_ERI_D2z_S_Px_S_C1000001_ad;
  Double I_ERI_D2x_S_Py_Py_C1000001_ad = I_ERI_D2x_S_D2y_S_C1000001_ad+CDY*I_ERI_D2x_S_Py_S_C1000001_ad;
  Double I_ERI_Dxy_S_Py_Py_C1000001_ad = I_ERI_Dxy_S_D2y_S_C1000001_ad+CDY*I_ERI_Dxy_S_Py_S_C1000001_ad;
  Double I_ERI_Dxz_S_Py_Py_C1000001_ad = I_ERI_Dxz_S_D2y_S_C1000001_ad+CDY*I_ERI_Dxz_S_Py_S_C1000001_ad;
  Double I_ERI_D2y_S_Py_Py_C1000001_ad = I_ERI_D2y_S_D2y_S_C1000001_ad+CDY*I_ERI_D2y_S_Py_S_C1000001_ad;
  Double I_ERI_Dyz_S_Py_Py_C1000001_ad = I_ERI_Dyz_S_D2y_S_C1000001_ad+CDY*I_ERI_Dyz_S_Py_S_C1000001_ad;
  Double I_ERI_D2z_S_Py_Py_C1000001_ad = I_ERI_D2z_S_D2y_S_C1000001_ad+CDY*I_ERI_D2z_S_Py_S_C1000001_ad;
  Double I_ERI_D2x_S_Pz_Py_C1000001_ad = I_ERI_D2x_S_Dyz_S_C1000001_ad+CDY*I_ERI_D2x_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxy_S_Pz_Py_C1000001_ad = I_ERI_Dxy_S_Dyz_S_C1000001_ad+CDY*I_ERI_Dxy_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxz_S_Pz_Py_C1000001_ad = I_ERI_Dxz_S_Dyz_S_C1000001_ad+CDY*I_ERI_Dxz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2y_S_Pz_Py_C1000001_ad = I_ERI_D2y_S_Dyz_S_C1000001_ad+CDY*I_ERI_D2y_S_Pz_S_C1000001_ad;
  Double I_ERI_Dyz_S_Pz_Py_C1000001_ad = I_ERI_Dyz_S_Dyz_S_C1000001_ad+CDY*I_ERI_Dyz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2z_S_Pz_Py_C1000001_ad = I_ERI_D2z_S_Dyz_S_C1000001_ad+CDY*I_ERI_D2z_S_Pz_S_C1000001_ad;
  Double I_ERI_D2x_S_Px_Pz_C1000001_ad = I_ERI_D2x_S_Dxz_S_C1000001_ad+CDZ*I_ERI_D2x_S_Px_S_C1000001_ad;
  Double I_ERI_Dxy_S_Px_Pz_C1000001_ad = I_ERI_Dxy_S_Dxz_S_C1000001_ad+CDZ*I_ERI_Dxy_S_Px_S_C1000001_ad;
  Double I_ERI_Dxz_S_Px_Pz_C1000001_ad = I_ERI_Dxz_S_Dxz_S_C1000001_ad+CDZ*I_ERI_Dxz_S_Px_S_C1000001_ad;
  Double I_ERI_D2y_S_Px_Pz_C1000001_ad = I_ERI_D2y_S_Dxz_S_C1000001_ad+CDZ*I_ERI_D2y_S_Px_S_C1000001_ad;
  Double I_ERI_Dyz_S_Px_Pz_C1000001_ad = I_ERI_Dyz_S_Dxz_S_C1000001_ad+CDZ*I_ERI_Dyz_S_Px_S_C1000001_ad;
  Double I_ERI_D2z_S_Px_Pz_C1000001_ad = I_ERI_D2z_S_Dxz_S_C1000001_ad+CDZ*I_ERI_D2z_S_Px_S_C1000001_ad;
  Double I_ERI_D2x_S_Py_Pz_C1000001_ad = I_ERI_D2x_S_Dyz_S_C1000001_ad+CDZ*I_ERI_D2x_S_Py_S_C1000001_ad;
  Double I_ERI_Dxy_S_Py_Pz_C1000001_ad = I_ERI_Dxy_S_Dyz_S_C1000001_ad+CDZ*I_ERI_Dxy_S_Py_S_C1000001_ad;
  Double I_ERI_Dxz_S_Py_Pz_C1000001_ad = I_ERI_Dxz_S_Dyz_S_C1000001_ad+CDZ*I_ERI_Dxz_S_Py_S_C1000001_ad;
  Double I_ERI_D2y_S_Py_Pz_C1000001_ad = I_ERI_D2y_S_Dyz_S_C1000001_ad+CDZ*I_ERI_D2y_S_Py_S_C1000001_ad;
  Double I_ERI_Dyz_S_Py_Pz_C1000001_ad = I_ERI_Dyz_S_Dyz_S_C1000001_ad+CDZ*I_ERI_Dyz_S_Py_S_C1000001_ad;
  Double I_ERI_D2z_S_Py_Pz_C1000001_ad = I_ERI_D2z_S_Dyz_S_C1000001_ad+CDZ*I_ERI_D2z_S_Py_S_C1000001_ad;
  Double I_ERI_D2x_S_Pz_Pz_C1000001_ad = I_ERI_D2x_S_D2z_S_C1000001_ad+CDZ*I_ERI_D2x_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxy_S_Pz_Pz_C1000001_ad = I_ERI_Dxy_S_D2z_S_C1000001_ad+CDZ*I_ERI_Dxy_S_Pz_S_C1000001_ad;
  Double I_ERI_Dxz_S_Pz_Pz_C1000001_ad = I_ERI_Dxz_S_D2z_S_C1000001_ad+CDZ*I_ERI_Dxz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2y_S_Pz_Pz_C1000001_ad = I_ERI_D2y_S_D2z_S_C1000001_ad+CDZ*I_ERI_D2y_S_Pz_S_C1000001_ad;
  Double I_ERI_Dyz_S_Pz_Pz_C1000001_ad = I_ERI_Dyz_S_D2z_S_C1000001_ad+CDZ*I_ERI_Dyz_S_Pz_S_C1000001_ad;
  Double I_ERI_D2z_S_Pz_Pz_C1000001_ad = I_ERI_D2z_S_D2z_S_C1000001_ad+CDZ*I_ERI_D2z_S_Pz_S_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001_ad
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_ad
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_ad
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001_ad = I_ERI_D2x_S_D2x_S_C1001001_ad+CDX*I_ERI_D2x_S_Px_S_C1001001_ad;
  Double I_ERI_Dxy_S_Px_Px_C1001001_ad = I_ERI_Dxy_S_D2x_S_C1001001_ad+CDX*I_ERI_Dxy_S_Px_S_C1001001_ad;
  Double I_ERI_Dxz_S_Px_Px_C1001001_ad = I_ERI_Dxz_S_D2x_S_C1001001_ad+CDX*I_ERI_Dxz_S_Px_S_C1001001_ad;
  Double I_ERI_D2y_S_Px_Px_C1001001_ad = I_ERI_D2y_S_D2x_S_C1001001_ad+CDX*I_ERI_D2y_S_Px_S_C1001001_ad;
  Double I_ERI_Dyz_S_Px_Px_C1001001_ad = I_ERI_Dyz_S_D2x_S_C1001001_ad+CDX*I_ERI_Dyz_S_Px_S_C1001001_ad;
  Double I_ERI_D2z_S_Px_Px_C1001001_ad = I_ERI_D2z_S_D2x_S_C1001001_ad+CDX*I_ERI_D2z_S_Px_S_C1001001_ad;
  Double I_ERI_D2x_S_Py_Px_C1001001_ad = I_ERI_D2x_S_Dxy_S_C1001001_ad+CDX*I_ERI_D2x_S_Py_S_C1001001_ad;
  Double I_ERI_Dxy_S_Py_Px_C1001001_ad = I_ERI_Dxy_S_Dxy_S_C1001001_ad+CDX*I_ERI_Dxy_S_Py_S_C1001001_ad;
  Double I_ERI_Dxz_S_Py_Px_C1001001_ad = I_ERI_Dxz_S_Dxy_S_C1001001_ad+CDX*I_ERI_Dxz_S_Py_S_C1001001_ad;
  Double I_ERI_D2y_S_Py_Px_C1001001_ad = I_ERI_D2y_S_Dxy_S_C1001001_ad+CDX*I_ERI_D2y_S_Py_S_C1001001_ad;
  Double I_ERI_Dyz_S_Py_Px_C1001001_ad = I_ERI_Dyz_S_Dxy_S_C1001001_ad+CDX*I_ERI_Dyz_S_Py_S_C1001001_ad;
  Double I_ERI_D2z_S_Py_Px_C1001001_ad = I_ERI_D2z_S_Dxy_S_C1001001_ad+CDX*I_ERI_D2z_S_Py_S_C1001001_ad;
  Double I_ERI_D2x_S_Pz_Px_C1001001_ad = I_ERI_D2x_S_Dxz_S_C1001001_ad+CDX*I_ERI_D2x_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxy_S_Pz_Px_C1001001_ad = I_ERI_Dxy_S_Dxz_S_C1001001_ad+CDX*I_ERI_Dxy_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxz_S_Pz_Px_C1001001_ad = I_ERI_Dxz_S_Dxz_S_C1001001_ad+CDX*I_ERI_Dxz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2y_S_Pz_Px_C1001001_ad = I_ERI_D2y_S_Dxz_S_C1001001_ad+CDX*I_ERI_D2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Dyz_S_Pz_Px_C1001001_ad = I_ERI_Dyz_S_Dxz_S_C1001001_ad+CDX*I_ERI_Dyz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2z_S_Pz_Px_C1001001_ad = I_ERI_D2z_S_Dxz_S_C1001001_ad+CDX*I_ERI_D2z_S_Pz_S_C1001001_ad;
  Double I_ERI_D2x_S_Px_Py_C1001001_ad = I_ERI_D2x_S_Dxy_S_C1001001_ad+CDY*I_ERI_D2x_S_Px_S_C1001001_ad;
  Double I_ERI_Dxy_S_Px_Py_C1001001_ad = I_ERI_Dxy_S_Dxy_S_C1001001_ad+CDY*I_ERI_Dxy_S_Px_S_C1001001_ad;
  Double I_ERI_Dxz_S_Px_Py_C1001001_ad = I_ERI_Dxz_S_Dxy_S_C1001001_ad+CDY*I_ERI_Dxz_S_Px_S_C1001001_ad;
  Double I_ERI_D2y_S_Px_Py_C1001001_ad = I_ERI_D2y_S_Dxy_S_C1001001_ad+CDY*I_ERI_D2y_S_Px_S_C1001001_ad;
  Double I_ERI_Dyz_S_Px_Py_C1001001_ad = I_ERI_Dyz_S_Dxy_S_C1001001_ad+CDY*I_ERI_Dyz_S_Px_S_C1001001_ad;
  Double I_ERI_D2z_S_Px_Py_C1001001_ad = I_ERI_D2z_S_Dxy_S_C1001001_ad+CDY*I_ERI_D2z_S_Px_S_C1001001_ad;
  Double I_ERI_D2x_S_Py_Py_C1001001_ad = I_ERI_D2x_S_D2y_S_C1001001_ad+CDY*I_ERI_D2x_S_Py_S_C1001001_ad;
  Double I_ERI_Dxy_S_Py_Py_C1001001_ad = I_ERI_Dxy_S_D2y_S_C1001001_ad+CDY*I_ERI_Dxy_S_Py_S_C1001001_ad;
  Double I_ERI_Dxz_S_Py_Py_C1001001_ad = I_ERI_Dxz_S_D2y_S_C1001001_ad+CDY*I_ERI_Dxz_S_Py_S_C1001001_ad;
  Double I_ERI_D2y_S_Py_Py_C1001001_ad = I_ERI_D2y_S_D2y_S_C1001001_ad+CDY*I_ERI_D2y_S_Py_S_C1001001_ad;
  Double I_ERI_Dyz_S_Py_Py_C1001001_ad = I_ERI_Dyz_S_D2y_S_C1001001_ad+CDY*I_ERI_Dyz_S_Py_S_C1001001_ad;
  Double I_ERI_D2z_S_Py_Py_C1001001_ad = I_ERI_D2z_S_D2y_S_C1001001_ad+CDY*I_ERI_D2z_S_Py_S_C1001001_ad;
  Double I_ERI_D2x_S_Pz_Py_C1001001_ad = I_ERI_D2x_S_Dyz_S_C1001001_ad+CDY*I_ERI_D2x_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxy_S_Pz_Py_C1001001_ad = I_ERI_Dxy_S_Dyz_S_C1001001_ad+CDY*I_ERI_Dxy_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxz_S_Pz_Py_C1001001_ad = I_ERI_Dxz_S_Dyz_S_C1001001_ad+CDY*I_ERI_Dxz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2y_S_Pz_Py_C1001001_ad = I_ERI_D2y_S_Dyz_S_C1001001_ad+CDY*I_ERI_D2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Dyz_S_Pz_Py_C1001001_ad = I_ERI_Dyz_S_Dyz_S_C1001001_ad+CDY*I_ERI_Dyz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2z_S_Pz_Py_C1001001_ad = I_ERI_D2z_S_Dyz_S_C1001001_ad+CDY*I_ERI_D2z_S_Pz_S_C1001001_ad;
  Double I_ERI_D2x_S_Px_Pz_C1001001_ad = I_ERI_D2x_S_Dxz_S_C1001001_ad+CDZ*I_ERI_D2x_S_Px_S_C1001001_ad;
  Double I_ERI_Dxy_S_Px_Pz_C1001001_ad = I_ERI_Dxy_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Dxy_S_Px_S_C1001001_ad;
  Double I_ERI_Dxz_S_Px_Pz_C1001001_ad = I_ERI_Dxz_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Dxz_S_Px_S_C1001001_ad;
  Double I_ERI_D2y_S_Px_Pz_C1001001_ad = I_ERI_D2y_S_Dxz_S_C1001001_ad+CDZ*I_ERI_D2y_S_Px_S_C1001001_ad;
  Double I_ERI_Dyz_S_Px_Pz_C1001001_ad = I_ERI_Dyz_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Dyz_S_Px_S_C1001001_ad;
  Double I_ERI_D2z_S_Px_Pz_C1001001_ad = I_ERI_D2z_S_Dxz_S_C1001001_ad+CDZ*I_ERI_D2z_S_Px_S_C1001001_ad;
  Double I_ERI_D2x_S_Py_Pz_C1001001_ad = I_ERI_D2x_S_Dyz_S_C1001001_ad+CDZ*I_ERI_D2x_S_Py_S_C1001001_ad;
  Double I_ERI_Dxy_S_Py_Pz_C1001001_ad = I_ERI_Dxy_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Dxy_S_Py_S_C1001001_ad;
  Double I_ERI_Dxz_S_Py_Pz_C1001001_ad = I_ERI_Dxz_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Dxz_S_Py_S_C1001001_ad;
  Double I_ERI_D2y_S_Py_Pz_C1001001_ad = I_ERI_D2y_S_Dyz_S_C1001001_ad+CDZ*I_ERI_D2y_S_Py_S_C1001001_ad;
  Double I_ERI_Dyz_S_Py_Pz_C1001001_ad = I_ERI_Dyz_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Dyz_S_Py_S_C1001001_ad;
  Double I_ERI_D2z_S_Py_Pz_C1001001_ad = I_ERI_D2z_S_Dyz_S_C1001001_ad+CDZ*I_ERI_D2z_S_Py_S_C1001001_ad;
  Double I_ERI_D2x_S_Pz_Pz_C1001001_ad = I_ERI_D2x_S_D2z_S_C1001001_ad+CDZ*I_ERI_D2x_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001_ad = I_ERI_Dxy_S_D2z_S_C1001001_ad+CDZ*I_ERI_Dxy_S_Pz_S_C1001001_ad;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001_ad = I_ERI_Dxz_S_D2z_S_C1001001_ad+CDZ*I_ERI_Dxz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2y_S_Pz_Pz_C1001001_ad = I_ERI_D2y_S_D2z_S_C1001001_ad+CDZ*I_ERI_D2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001_ad = I_ERI_Dyz_S_D2z_S_C1001001_ad+CDZ*I_ERI_Dyz_S_Pz_S_C1001001_ad;
  Double I_ERI_D2z_S_Pz_Pz_C1001001_ad = I_ERI_D2z_S_D2z_S_C1001001_ad+CDZ*I_ERI_D2z_S_Pz_S_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001001_ad
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001001_ad
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001_ad
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1001001_ad = I_ERI_F3x_S_D2x_S_C1001001_ad+CDX*I_ERI_F3x_S_Px_S_C1001001_ad;
  Double I_ERI_F2xy_S_Px_Px_C1001001_ad = I_ERI_F2xy_S_D2x_S_C1001001_ad+CDX*I_ERI_F2xy_S_Px_S_C1001001_ad;
  Double I_ERI_F2xz_S_Px_Px_C1001001_ad = I_ERI_F2xz_S_D2x_S_C1001001_ad+CDX*I_ERI_F2xz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Px_Px_C1001001_ad = I_ERI_Fx2y_S_D2x_S_C1001001_ad+CDX*I_ERI_Fx2y_S_Px_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Px_Px_C1001001_ad = I_ERI_Fxyz_S_D2x_S_C1001001_ad+CDX*I_ERI_Fxyz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Px_Px_C1001001_ad = I_ERI_Fx2z_S_D2x_S_C1001001_ad+CDX*I_ERI_Fx2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3y_S_Px_Px_C1001001_ad = I_ERI_F3y_S_D2x_S_C1001001_ad+CDX*I_ERI_F3y_S_Px_S_C1001001_ad;
  Double I_ERI_F2yz_S_Px_Px_C1001001_ad = I_ERI_F2yz_S_D2x_S_C1001001_ad+CDX*I_ERI_F2yz_S_Px_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Px_Px_C1001001_ad = I_ERI_Fy2z_S_D2x_S_C1001001_ad+CDX*I_ERI_Fy2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3z_S_Px_Px_C1001001_ad = I_ERI_F3z_S_D2x_S_C1001001_ad+CDX*I_ERI_F3z_S_Px_S_C1001001_ad;
  Double I_ERI_F3x_S_Py_Px_C1001001_ad = I_ERI_F3x_S_Dxy_S_C1001001_ad+CDX*I_ERI_F3x_S_Py_S_C1001001_ad;
  Double I_ERI_F2xy_S_Py_Px_C1001001_ad = I_ERI_F2xy_S_Dxy_S_C1001001_ad+CDX*I_ERI_F2xy_S_Py_S_C1001001_ad;
  Double I_ERI_F2xz_S_Py_Px_C1001001_ad = I_ERI_F2xz_S_Dxy_S_C1001001_ad+CDX*I_ERI_F2xz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Py_Px_C1001001_ad = I_ERI_Fx2y_S_Dxy_S_C1001001_ad+CDX*I_ERI_Fx2y_S_Py_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Py_Px_C1001001_ad = I_ERI_Fxyz_S_Dxy_S_C1001001_ad+CDX*I_ERI_Fxyz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Py_Px_C1001001_ad = I_ERI_Fx2z_S_Dxy_S_C1001001_ad+CDX*I_ERI_Fx2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3y_S_Py_Px_C1001001_ad = I_ERI_F3y_S_Dxy_S_C1001001_ad+CDX*I_ERI_F3y_S_Py_S_C1001001_ad;
  Double I_ERI_F2yz_S_Py_Px_C1001001_ad = I_ERI_F2yz_S_Dxy_S_C1001001_ad+CDX*I_ERI_F2yz_S_Py_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Py_Px_C1001001_ad = I_ERI_Fy2z_S_Dxy_S_C1001001_ad+CDX*I_ERI_Fy2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3z_S_Py_Px_C1001001_ad = I_ERI_F3z_S_Dxy_S_C1001001_ad+CDX*I_ERI_F3z_S_Py_S_C1001001_ad;
  Double I_ERI_F3x_S_Pz_Px_C1001001_ad = I_ERI_F3x_S_Dxz_S_C1001001_ad+CDX*I_ERI_F3x_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xy_S_Pz_Px_C1001001_ad = I_ERI_F2xy_S_Dxz_S_C1001001_ad+CDX*I_ERI_F2xy_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xz_S_Pz_Px_C1001001_ad = I_ERI_F2xz_S_Dxz_S_C1001001_ad+CDX*I_ERI_F2xz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Pz_Px_C1001001_ad = I_ERI_Fx2y_S_Dxz_S_C1001001_ad+CDX*I_ERI_Fx2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Pz_Px_C1001001_ad = I_ERI_Fxyz_S_Dxz_S_C1001001_ad+CDX*I_ERI_Fxyz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Pz_Px_C1001001_ad = I_ERI_Fx2z_S_Dxz_S_C1001001_ad+CDX*I_ERI_Fx2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3y_S_Pz_Px_C1001001_ad = I_ERI_F3y_S_Dxz_S_C1001001_ad+CDX*I_ERI_F3y_S_Pz_S_C1001001_ad;
  Double I_ERI_F2yz_S_Pz_Px_C1001001_ad = I_ERI_F2yz_S_Dxz_S_C1001001_ad+CDX*I_ERI_F2yz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Pz_Px_C1001001_ad = I_ERI_Fy2z_S_Dxz_S_C1001001_ad+CDX*I_ERI_Fy2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3z_S_Pz_Px_C1001001_ad = I_ERI_F3z_S_Dxz_S_C1001001_ad+CDX*I_ERI_F3z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3x_S_Px_Py_C1001001_ad = I_ERI_F3x_S_Dxy_S_C1001001_ad+CDY*I_ERI_F3x_S_Px_S_C1001001_ad;
  Double I_ERI_F2xy_S_Px_Py_C1001001_ad = I_ERI_F2xy_S_Dxy_S_C1001001_ad+CDY*I_ERI_F2xy_S_Px_S_C1001001_ad;
  Double I_ERI_F2xz_S_Px_Py_C1001001_ad = I_ERI_F2xz_S_Dxy_S_C1001001_ad+CDY*I_ERI_F2xz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Px_Py_C1001001_ad = I_ERI_Fx2y_S_Dxy_S_C1001001_ad+CDY*I_ERI_Fx2y_S_Px_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Px_Py_C1001001_ad = I_ERI_Fxyz_S_Dxy_S_C1001001_ad+CDY*I_ERI_Fxyz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Px_Py_C1001001_ad = I_ERI_Fx2z_S_Dxy_S_C1001001_ad+CDY*I_ERI_Fx2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3y_S_Px_Py_C1001001_ad = I_ERI_F3y_S_Dxy_S_C1001001_ad+CDY*I_ERI_F3y_S_Px_S_C1001001_ad;
  Double I_ERI_F2yz_S_Px_Py_C1001001_ad = I_ERI_F2yz_S_Dxy_S_C1001001_ad+CDY*I_ERI_F2yz_S_Px_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Px_Py_C1001001_ad = I_ERI_Fy2z_S_Dxy_S_C1001001_ad+CDY*I_ERI_Fy2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3z_S_Px_Py_C1001001_ad = I_ERI_F3z_S_Dxy_S_C1001001_ad+CDY*I_ERI_F3z_S_Px_S_C1001001_ad;
  Double I_ERI_F3x_S_Py_Py_C1001001_ad = I_ERI_F3x_S_D2y_S_C1001001_ad+CDY*I_ERI_F3x_S_Py_S_C1001001_ad;
  Double I_ERI_F2xy_S_Py_Py_C1001001_ad = I_ERI_F2xy_S_D2y_S_C1001001_ad+CDY*I_ERI_F2xy_S_Py_S_C1001001_ad;
  Double I_ERI_F2xz_S_Py_Py_C1001001_ad = I_ERI_F2xz_S_D2y_S_C1001001_ad+CDY*I_ERI_F2xz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Py_Py_C1001001_ad = I_ERI_Fx2y_S_D2y_S_C1001001_ad+CDY*I_ERI_Fx2y_S_Py_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Py_Py_C1001001_ad = I_ERI_Fxyz_S_D2y_S_C1001001_ad+CDY*I_ERI_Fxyz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Py_Py_C1001001_ad = I_ERI_Fx2z_S_D2y_S_C1001001_ad+CDY*I_ERI_Fx2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3y_S_Py_Py_C1001001_ad = I_ERI_F3y_S_D2y_S_C1001001_ad+CDY*I_ERI_F3y_S_Py_S_C1001001_ad;
  Double I_ERI_F2yz_S_Py_Py_C1001001_ad = I_ERI_F2yz_S_D2y_S_C1001001_ad+CDY*I_ERI_F2yz_S_Py_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Py_Py_C1001001_ad = I_ERI_Fy2z_S_D2y_S_C1001001_ad+CDY*I_ERI_Fy2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3z_S_Py_Py_C1001001_ad = I_ERI_F3z_S_D2y_S_C1001001_ad+CDY*I_ERI_F3z_S_Py_S_C1001001_ad;
  Double I_ERI_F3x_S_Pz_Py_C1001001_ad = I_ERI_F3x_S_Dyz_S_C1001001_ad+CDY*I_ERI_F3x_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xy_S_Pz_Py_C1001001_ad = I_ERI_F2xy_S_Dyz_S_C1001001_ad+CDY*I_ERI_F2xy_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xz_S_Pz_Py_C1001001_ad = I_ERI_F2xz_S_Dyz_S_C1001001_ad+CDY*I_ERI_F2xz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Pz_Py_C1001001_ad = I_ERI_Fx2y_S_Dyz_S_C1001001_ad+CDY*I_ERI_Fx2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Pz_Py_C1001001_ad = I_ERI_Fxyz_S_Dyz_S_C1001001_ad+CDY*I_ERI_Fxyz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Pz_Py_C1001001_ad = I_ERI_Fx2z_S_Dyz_S_C1001001_ad+CDY*I_ERI_Fx2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3y_S_Pz_Py_C1001001_ad = I_ERI_F3y_S_Dyz_S_C1001001_ad+CDY*I_ERI_F3y_S_Pz_S_C1001001_ad;
  Double I_ERI_F2yz_S_Pz_Py_C1001001_ad = I_ERI_F2yz_S_Dyz_S_C1001001_ad+CDY*I_ERI_F2yz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Pz_Py_C1001001_ad = I_ERI_Fy2z_S_Dyz_S_C1001001_ad+CDY*I_ERI_Fy2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3z_S_Pz_Py_C1001001_ad = I_ERI_F3z_S_Dyz_S_C1001001_ad+CDY*I_ERI_F3z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3x_S_Px_Pz_C1001001_ad = I_ERI_F3x_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F3x_S_Px_S_C1001001_ad;
  Double I_ERI_F2xy_S_Px_Pz_C1001001_ad = I_ERI_F2xy_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F2xy_S_Px_S_C1001001_ad;
  Double I_ERI_F2xz_S_Px_Pz_C1001001_ad = I_ERI_F2xz_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F2xz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Px_Pz_C1001001_ad = I_ERI_Fx2y_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Fx2y_S_Px_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Px_Pz_C1001001_ad = I_ERI_Fxyz_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Fxyz_S_Px_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Px_Pz_C1001001_ad = I_ERI_Fx2z_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Fx2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3y_S_Px_Pz_C1001001_ad = I_ERI_F3y_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F3y_S_Px_S_C1001001_ad;
  Double I_ERI_F2yz_S_Px_Pz_C1001001_ad = I_ERI_F2yz_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F2yz_S_Px_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Px_Pz_C1001001_ad = I_ERI_Fy2z_S_Dxz_S_C1001001_ad+CDZ*I_ERI_Fy2z_S_Px_S_C1001001_ad;
  Double I_ERI_F3z_S_Px_Pz_C1001001_ad = I_ERI_F3z_S_Dxz_S_C1001001_ad+CDZ*I_ERI_F3z_S_Px_S_C1001001_ad;
  Double I_ERI_F3x_S_Py_Pz_C1001001_ad = I_ERI_F3x_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F3x_S_Py_S_C1001001_ad;
  Double I_ERI_F2xy_S_Py_Pz_C1001001_ad = I_ERI_F2xy_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F2xy_S_Py_S_C1001001_ad;
  Double I_ERI_F2xz_S_Py_Pz_C1001001_ad = I_ERI_F2xz_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F2xz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Py_Pz_C1001001_ad = I_ERI_Fx2y_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Fx2y_S_Py_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Py_Pz_C1001001_ad = I_ERI_Fxyz_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Fxyz_S_Py_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Py_Pz_C1001001_ad = I_ERI_Fx2z_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Fx2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3y_S_Py_Pz_C1001001_ad = I_ERI_F3y_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F3y_S_Py_S_C1001001_ad;
  Double I_ERI_F2yz_S_Py_Pz_C1001001_ad = I_ERI_F2yz_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F2yz_S_Py_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Py_Pz_C1001001_ad = I_ERI_Fy2z_S_Dyz_S_C1001001_ad+CDZ*I_ERI_Fy2z_S_Py_S_C1001001_ad;
  Double I_ERI_F3z_S_Py_Pz_C1001001_ad = I_ERI_F3z_S_Dyz_S_C1001001_ad+CDZ*I_ERI_F3z_S_Py_S_C1001001_ad;
  Double I_ERI_F3x_S_Pz_Pz_C1001001_ad = I_ERI_F3x_S_D2z_S_C1001001_ad+CDZ*I_ERI_F3x_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xy_S_Pz_Pz_C1001001_ad = I_ERI_F2xy_S_D2z_S_C1001001_ad+CDZ*I_ERI_F2xy_S_Pz_S_C1001001_ad;
  Double I_ERI_F2xz_S_Pz_Pz_C1001001_ad = I_ERI_F2xz_S_D2z_S_C1001001_ad+CDZ*I_ERI_F2xz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2y_S_Pz_Pz_C1001001_ad = I_ERI_Fx2y_S_D2z_S_C1001001_ad+CDZ*I_ERI_Fx2y_S_Pz_S_C1001001_ad;
  Double I_ERI_Fxyz_S_Pz_Pz_C1001001_ad = I_ERI_Fxyz_S_D2z_S_C1001001_ad+CDZ*I_ERI_Fxyz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fx2z_S_Pz_Pz_C1001001_ad = I_ERI_Fx2z_S_D2z_S_C1001001_ad+CDZ*I_ERI_Fx2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3y_S_Pz_Pz_C1001001_ad = I_ERI_F3y_S_D2z_S_C1001001_ad+CDZ*I_ERI_F3y_S_Pz_S_C1001001_ad;
  Double I_ERI_F2yz_S_Pz_Pz_C1001001_ad = I_ERI_F2yz_S_D2z_S_C1001001_ad+CDZ*I_ERI_F2yz_S_Pz_S_C1001001_ad;
  Double I_ERI_Fy2z_S_Pz_Pz_C1001001_ad = I_ERI_Fy2z_S_D2z_S_C1001001_ad+CDZ*I_ERI_Fy2z_S_Pz_S_C1001001_ad;
  Double I_ERI_F3z_S_Pz_Pz_C1001001_ad = I_ERI_F3z_S_D2z_S_C1001001_ad+CDZ*I_ERI_F3z_S_Pz_S_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_cd
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1000001_cd = I_ERI_Px_S_F3x_S_C1000001_cd+CDX*I_ERI_Px_S_D2x_S_C1000001_cd;
  Double I_ERI_Py_S_D2x_Px_C1000001_cd = I_ERI_Py_S_F3x_S_C1000001_cd+CDX*I_ERI_Py_S_D2x_S_C1000001_cd;
  Double I_ERI_Pz_S_D2x_Px_C1000001_cd = I_ERI_Pz_S_F3x_S_C1000001_cd+CDX*I_ERI_Pz_S_D2x_S_C1000001_cd;
  Double I_ERI_Px_S_Dxy_Px_C1000001_cd = I_ERI_Px_S_F2xy_S_C1000001_cd+CDX*I_ERI_Px_S_Dxy_S_C1000001_cd;
  Double I_ERI_Py_S_Dxy_Px_C1000001_cd = I_ERI_Py_S_F2xy_S_C1000001_cd+CDX*I_ERI_Py_S_Dxy_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxy_Px_C1000001_cd = I_ERI_Pz_S_F2xy_S_C1000001_cd+CDX*I_ERI_Pz_S_Dxy_S_C1000001_cd;
  Double I_ERI_Px_S_Dxz_Px_C1000001_cd = I_ERI_Px_S_F2xz_S_C1000001_cd+CDX*I_ERI_Px_S_Dxz_S_C1000001_cd;
  Double I_ERI_Py_S_Dxz_Px_C1000001_cd = I_ERI_Py_S_F2xz_S_C1000001_cd+CDX*I_ERI_Py_S_Dxz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxz_Px_C1000001_cd = I_ERI_Pz_S_F2xz_S_C1000001_cd+CDX*I_ERI_Pz_S_Dxz_S_C1000001_cd;
  Double I_ERI_Px_S_D2y_Px_C1000001_cd = I_ERI_Px_S_Fx2y_S_C1000001_cd+CDX*I_ERI_Px_S_D2y_S_C1000001_cd;
  Double I_ERI_Py_S_D2y_Px_C1000001_cd = I_ERI_Py_S_Fx2y_S_C1000001_cd+CDX*I_ERI_Py_S_D2y_S_C1000001_cd;
  Double I_ERI_Pz_S_D2y_Px_C1000001_cd = I_ERI_Pz_S_Fx2y_S_C1000001_cd+CDX*I_ERI_Pz_S_D2y_S_C1000001_cd;
  Double I_ERI_Px_S_Dyz_Px_C1000001_cd = I_ERI_Px_S_Fxyz_S_C1000001_cd+CDX*I_ERI_Px_S_Dyz_S_C1000001_cd;
  Double I_ERI_Py_S_Dyz_Px_C1000001_cd = I_ERI_Py_S_Fxyz_S_C1000001_cd+CDX*I_ERI_Py_S_Dyz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dyz_Px_C1000001_cd = I_ERI_Pz_S_Fxyz_S_C1000001_cd+CDX*I_ERI_Pz_S_Dyz_S_C1000001_cd;
  Double I_ERI_Px_S_D2z_Px_C1000001_cd = I_ERI_Px_S_Fx2z_S_C1000001_cd+CDX*I_ERI_Px_S_D2z_S_C1000001_cd;
  Double I_ERI_Py_S_D2z_Px_C1000001_cd = I_ERI_Py_S_Fx2z_S_C1000001_cd+CDX*I_ERI_Py_S_D2z_S_C1000001_cd;
  Double I_ERI_Pz_S_D2z_Px_C1000001_cd = I_ERI_Pz_S_Fx2z_S_C1000001_cd+CDX*I_ERI_Pz_S_D2z_S_C1000001_cd;
  Double I_ERI_Px_S_D2x_Py_C1000001_cd = I_ERI_Px_S_F2xy_S_C1000001_cd+CDY*I_ERI_Px_S_D2x_S_C1000001_cd;
  Double I_ERI_Py_S_D2x_Py_C1000001_cd = I_ERI_Py_S_F2xy_S_C1000001_cd+CDY*I_ERI_Py_S_D2x_S_C1000001_cd;
  Double I_ERI_Pz_S_D2x_Py_C1000001_cd = I_ERI_Pz_S_F2xy_S_C1000001_cd+CDY*I_ERI_Pz_S_D2x_S_C1000001_cd;
  Double I_ERI_Px_S_Dxy_Py_C1000001_cd = I_ERI_Px_S_Fx2y_S_C1000001_cd+CDY*I_ERI_Px_S_Dxy_S_C1000001_cd;
  Double I_ERI_Py_S_Dxy_Py_C1000001_cd = I_ERI_Py_S_Fx2y_S_C1000001_cd+CDY*I_ERI_Py_S_Dxy_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxy_Py_C1000001_cd = I_ERI_Pz_S_Fx2y_S_C1000001_cd+CDY*I_ERI_Pz_S_Dxy_S_C1000001_cd;
  Double I_ERI_Px_S_Dxz_Py_C1000001_cd = I_ERI_Px_S_Fxyz_S_C1000001_cd+CDY*I_ERI_Px_S_Dxz_S_C1000001_cd;
  Double I_ERI_Py_S_Dxz_Py_C1000001_cd = I_ERI_Py_S_Fxyz_S_C1000001_cd+CDY*I_ERI_Py_S_Dxz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxz_Py_C1000001_cd = I_ERI_Pz_S_Fxyz_S_C1000001_cd+CDY*I_ERI_Pz_S_Dxz_S_C1000001_cd;
  Double I_ERI_Px_S_D2y_Py_C1000001_cd = I_ERI_Px_S_F3y_S_C1000001_cd+CDY*I_ERI_Px_S_D2y_S_C1000001_cd;
  Double I_ERI_Py_S_D2y_Py_C1000001_cd = I_ERI_Py_S_F3y_S_C1000001_cd+CDY*I_ERI_Py_S_D2y_S_C1000001_cd;
  Double I_ERI_Pz_S_D2y_Py_C1000001_cd = I_ERI_Pz_S_F3y_S_C1000001_cd+CDY*I_ERI_Pz_S_D2y_S_C1000001_cd;
  Double I_ERI_Px_S_Dyz_Py_C1000001_cd = I_ERI_Px_S_F2yz_S_C1000001_cd+CDY*I_ERI_Px_S_Dyz_S_C1000001_cd;
  Double I_ERI_Py_S_Dyz_Py_C1000001_cd = I_ERI_Py_S_F2yz_S_C1000001_cd+CDY*I_ERI_Py_S_Dyz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dyz_Py_C1000001_cd = I_ERI_Pz_S_F2yz_S_C1000001_cd+CDY*I_ERI_Pz_S_Dyz_S_C1000001_cd;
  Double I_ERI_Px_S_D2z_Py_C1000001_cd = I_ERI_Px_S_Fy2z_S_C1000001_cd+CDY*I_ERI_Px_S_D2z_S_C1000001_cd;
  Double I_ERI_Py_S_D2z_Py_C1000001_cd = I_ERI_Py_S_Fy2z_S_C1000001_cd+CDY*I_ERI_Py_S_D2z_S_C1000001_cd;
  Double I_ERI_Pz_S_D2z_Py_C1000001_cd = I_ERI_Pz_S_Fy2z_S_C1000001_cd+CDY*I_ERI_Pz_S_D2z_S_C1000001_cd;
  Double I_ERI_Px_S_D2x_Pz_C1000001_cd = I_ERI_Px_S_F2xz_S_C1000001_cd+CDZ*I_ERI_Px_S_D2x_S_C1000001_cd;
  Double I_ERI_Py_S_D2x_Pz_C1000001_cd = I_ERI_Py_S_F2xz_S_C1000001_cd+CDZ*I_ERI_Py_S_D2x_S_C1000001_cd;
  Double I_ERI_Pz_S_D2x_Pz_C1000001_cd = I_ERI_Pz_S_F2xz_S_C1000001_cd+CDZ*I_ERI_Pz_S_D2x_S_C1000001_cd;
  Double I_ERI_Px_S_Dxy_Pz_C1000001_cd = I_ERI_Px_S_Fxyz_S_C1000001_cd+CDZ*I_ERI_Px_S_Dxy_S_C1000001_cd;
  Double I_ERI_Py_S_Dxy_Pz_C1000001_cd = I_ERI_Py_S_Fxyz_S_C1000001_cd+CDZ*I_ERI_Py_S_Dxy_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxy_Pz_C1000001_cd = I_ERI_Pz_S_Fxyz_S_C1000001_cd+CDZ*I_ERI_Pz_S_Dxy_S_C1000001_cd;
  Double I_ERI_Px_S_Dxz_Pz_C1000001_cd = I_ERI_Px_S_Fx2z_S_C1000001_cd+CDZ*I_ERI_Px_S_Dxz_S_C1000001_cd;
  Double I_ERI_Py_S_Dxz_Pz_C1000001_cd = I_ERI_Py_S_Fx2z_S_C1000001_cd+CDZ*I_ERI_Py_S_Dxz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dxz_Pz_C1000001_cd = I_ERI_Pz_S_Fx2z_S_C1000001_cd+CDZ*I_ERI_Pz_S_Dxz_S_C1000001_cd;
  Double I_ERI_Px_S_D2y_Pz_C1000001_cd = I_ERI_Px_S_F2yz_S_C1000001_cd+CDZ*I_ERI_Px_S_D2y_S_C1000001_cd;
  Double I_ERI_Py_S_D2y_Pz_C1000001_cd = I_ERI_Py_S_F2yz_S_C1000001_cd+CDZ*I_ERI_Py_S_D2y_S_C1000001_cd;
  Double I_ERI_Pz_S_D2y_Pz_C1000001_cd = I_ERI_Pz_S_F2yz_S_C1000001_cd+CDZ*I_ERI_Pz_S_D2y_S_C1000001_cd;
  Double I_ERI_Px_S_Dyz_Pz_C1000001_cd = I_ERI_Px_S_Fy2z_S_C1000001_cd+CDZ*I_ERI_Px_S_Dyz_S_C1000001_cd;
  Double I_ERI_Py_S_Dyz_Pz_C1000001_cd = I_ERI_Py_S_Fy2z_S_C1000001_cd+CDZ*I_ERI_Py_S_Dyz_S_C1000001_cd;
  Double I_ERI_Pz_S_Dyz_Pz_C1000001_cd = I_ERI_Pz_S_Fy2z_S_C1000001_cd+CDZ*I_ERI_Pz_S_Dyz_S_C1000001_cd;
  Double I_ERI_Px_S_D2z_Pz_C1000001_cd = I_ERI_Px_S_F3z_S_C1000001_cd+CDZ*I_ERI_Px_S_D2z_S_C1000001_cd;
  Double I_ERI_Py_S_D2z_Pz_C1000001_cd = I_ERI_Py_S_F3z_S_C1000001_cd+CDZ*I_ERI_Py_S_D2z_S_C1000001_cd;
  Double I_ERI_Pz_S_D2z_Pz_C1000001_cd = I_ERI_Pz_S_F3z_S_C1000001_cd+CDZ*I_ERI_Pz_S_D2z_S_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1001001_cd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001_cd
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1001001_cd = I_ERI_Px_S_F3x_S_C1001001_cd+CDX*I_ERI_Px_S_D2x_S_C1001001_cd;
  Double I_ERI_Py_S_D2x_Px_C1001001_cd = I_ERI_Py_S_F3x_S_C1001001_cd+CDX*I_ERI_Py_S_D2x_S_C1001001_cd;
  Double I_ERI_Pz_S_D2x_Px_C1001001_cd = I_ERI_Pz_S_F3x_S_C1001001_cd+CDX*I_ERI_Pz_S_D2x_S_C1001001_cd;
  Double I_ERI_Px_S_Dxy_Px_C1001001_cd = I_ERI_Px_S_F2xy_S_C1001001_cd+CDX*I_ERI_Px_S_Dxy_S_C1001001_cd;
  Double I_ERI_Py_S_Dxy_Px_C1001001_cd = I_ERI_Py_S_F2xy_S_C1001001_cd+CDX*I_ERI_Py_S_Dxy_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxy_Px_C1001001_cd = I_ERI_Pz_S_F2xy_S_C1001001_cd+CDX*I_ERI_Pz_S_Dxy_S_C1001001_cd;
  Double I_ERI_Px_S_Dxz_Px_C1001001_cd = I_ERI_Px_S_F2xz_S_C1001001_cd+CDX*I_ERI_Px_S_Dxz_S_C1001001_cd;
  Double I_ERI_Py_S_Dxz_Px_C1001001_cd = I_ERI_Py_S_F2xz_S_C1001001_cd+CDX*I_ERI_Py_S_Dxz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxz_Px_C1001001_cd = I_ERI_Pz_S_F2xz_S_C1001001_cd+CDX*I_ERI_Pz_S_Dxz_S_C1001001_cd;
  Double I_ERI_Px_S_D2y_Px_C1001001_cd = I_ERI_Px_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Px_S_D2y_S_C1001001_cd;
  Double I_ERI_Py_S_D2y_Px_C1001001_cd = I_ERI_Py_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Py_S_D2y_S_C1001001_cd;
  Double I_ERI_Pz_S_D2y_Px_C1001001_cd = I_ERI_Pz_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Pz_S_D2y_S_C1001001_cd;
  Double I_ERI_Px_S_Dyz_Px_C1001001_cd = I_ERI_Px_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Px_S_Dyz_S_C1001001_cd;
  Double I_ERI_Py_S_Dyz_Px_C1001001_cd = I_ERI_Py_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Py_S_Dyz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dyz_Px_C1001001_cd = I_ERI_Pz_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Pz_S_Dyz_S_C1001001_cd;
  Double I_ERI_Px_S_D2z_Px_C1001001_cd = I_ERI_Px_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Px_S_D2z_S_C1001001_cd;
  Double I_ERI_Py_S_D2z_Px_C1001001_cd = I_ERI_Py_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Py_S_D2z_S_C1001001_cd;
  Double I_ERI_Pz_S_D2z_Px_C1001001_cd = I_ERI_Pz_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Pz_S_D2z_S_C1001001_cd;
  Double I_ERI_Px_S_D2x_Py_C1001001_cd = I_ERI_Px_S_F2xy_S_C1001001_cd+CDY*I_ERI_Px_S_D2x_S_C1001001_cd;
  Double I_ERI_Py_S_D2x_Py_C1001001_cd = I_ERI_Py_S_F2xy_S_C1001001_cd+CDY*I_ERI_Py_S_D2x_S_C1001001_cd;
  Double I_ERI_Pz_S_D2x_Py_C1001001_cd = I_ERI_Pz_S_F2xy_S_C1001001_cd+CDY*I_ERI_Pz_S_D2x_S_C1001001_cd;
  Double I_ERI_Px_S_Dxy_Py_C1001001_cd = I_ERI_Px_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Px_S_Dxy_S_C1001001_cd;
  Double I_ERI_Py_S_Dxy_Py_C1001001_cd = I_ERI_Py_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Py_S_Dxy_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxy_Py_C1001001_cd = I_ERI_Pz_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Pz_S_Dxy_S_C1001001_cd;
  Double I_ERI_Px_S_Dxz_Py_C1001001_cd = I_ERI_Px_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Px_S_Dxz_S_C1001001_cd;
  Double I_ERI_Py_S_Dxz_Py_C1001001_cd = I_ERI_Py_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Py_S_Dxz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxz_Py_C1001001_cd = I_ERI_Pz_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Pz_S_Dxz_S_C1001001_cd;
  Double I_ERI_Px_S_D2y_Py_C1001001_cd = I_ERI_Px_S_F3y_S_C1001001_cd+CDY*I_ERI_Px_S_D2y_S_C1001001_cd;
  Double I_ERI_Py_S_D2y_Py_C1001001_cd = I_ERI_Py_S_F3y_S_C1001001_cd+CDY*I_ERI_Py_S_D2y_S_C1001001_cd;
  Double I_ERI_Pz_S_D2y_Py_C1001001_cd = I_ERI_Pz_S_F3y_S_C1001001_cd+CDY*I_ERI_Pz_S_D2y_S_C1001001_cd;
  Double I_ERI_Px_S_Dyz_Py_C1001001_cd = I_ERI_Px_S_F2yz_S_C1001001_cd+CDY*I_ERI_Px_S_Dyz_S_C1001001_cd;
  Double I_ERI_Py_S_Dyz_Py_C1001001_cd = I_ERI_Py_S_F2yz_S_C1001001_cd+CDY*I_ERI_Py_S_Dyz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dyz_Py_C1001001_cd = I_ERI_Pz_S_F2yz_S_C1001001_cd+CDY*I_ERI_Pz_S_Dyz_S_C1001001_cd;
  Double I_ERI_Px_S_D2z_Py_C1001001_cd = I_ERI_Px_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Px_S_D2z_S_C1001001_cd;
  Double I_ERI_Py_S_D2z_Py_C1001001_cd = I_ERI_Py_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Py_S_D2z_S_C1001001_cd;
  Double I_ERI_Pz_S_D2z_Py_C1001001_cd = I_ERI_Pz_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Pz_S_D2z_S_C1001001_cd;
  Double I_ERI_Px_S_D2x_Pz_C1001001_cd = I_ERI_Px_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Px_S_D2x_S_C1001001_cd;
  Double I_ERI_Py_S_D2x_Pz_C1001001_cd = I_ERI_Py_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Py_S_D2x_S_C1001001_cd;
  Double I_ERI_Pz_S_D2x_Pz_C1001001_cd = I_ERI_Pz_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Pz_S_D2x_S_C1001001_cd;
  Double I_ERI_Px_S_Dxy_Pz_C1001001_cd = I_ERI_Px_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Px_S_Dxy_S_C1001001_cd;
  Double I_ERI_Py_S_Dxy_Pz_C1001001_cd = I_ERI_Py_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Py_S_Dxy_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxy_Pz_C1001001_cd = I_ERI_Pz_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Pz_S_Dxy_S_C1001001_cd;
  Double I_ERI_Px_S_Dxz_Pz_C1001001_cd = I_ERI_Px_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Px_S_Dxz_S_C1001001_cd;
  Double I_ERI_Py_S_Dxz_Pz_C1001001_cd = I_ERI_Py_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Py_S_Dxz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dxz_Pz_C1001001_cd = I_ERI_Pz_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Pz_S_Dxz_S_C1001001_cd;
  Double I_ERI_Px_S_D2y_Pz_C1001001_cd = I_ERI_Px_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Px_S_D2y_S_C1001001_cd;
  Double I_ERI_Py_S_D2y_Pz_C1001001_cd = I_ERI_Py_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Py_S_D2y_S_C1001001_cd;
  Double I_ERI_Pz_S_D2y_Pz_C1001001_cd = I_ERI_Pz_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Pz_S_D2y_S_C1001001_cd;
  Double I_ERI_Px_S_Dyz_Pz_C1001001_cd = I_ERI_Px_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Px_S_Dyz_S_C1001001_cd;
  Double I_ERI_Py_S_Dyz_Pz_C1001001_cd = I_ERI_Py_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Py_S_Dyz_S_C1001001_cd;
  Double I_ERI_Pz_S_Dyz_Pz_C1001001_cd = I_ERI_Pz_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Pz_S_Dyz_S_C1001001_cd;
  Double I_ERI_Px_S_D2z_Pz_C1001001_cd = I_ERI_Px_S_F3z_S_C1001001_cd+CDZ*I_ERI_Px_S_D2z_S_C1001001_cd;
  Double I_ERI_Py_S_D2z_Pz_C1001001_cd = I_ERI_Py_S_F3z_S_C1001001_cd+CDZ*I_ERI_Py_S_D2z_S_C1001001_cd;
  Double I_ERI_Pz_S_D2z_Pz_C1001001_cd = I_ERI_Pz_S_F3z_S_C1001001_cd+CDZ*I_ERI_Pz_S_D2z_S_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_D_P_C1001001_cd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_C1001001_cd
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_cd
   ************************************************************/
  Double I_ERI_D2x_S_D2x_Px_C1001001_cd = I_ERI_D2x_S_F3x_S_C1001001_cd+CDX*I_ERI_D2x_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2x_Px_C1001001_cd = I_ERI_Dxy_S_F3x_S_C1001001_cd+CDX*I_ERI_Dxy_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2x_Px_C1001001_cd = I_ERI_Dxz_S_F3x_S_C1001001_cd+CDX*I_ERI_Dxz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2y_S_D2x_Px_C1001001_cd = I_ERI_D2y_S_F3x_S_C1001001_cd+CDX*I_ERI_D2y_S_D2x_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2x_Px_C1001001_cd = I_ERI_Dyz_S_F3x_S_C1001001_cd+CDX*I_ERI_Dyz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2z_S_D2x_Px_C1001001_cd = I_ERI_D2z_S_F3x_S_C1001001_cd+CDX*I_ERI_D2z_S_D2x_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxy_Px_C1001001_cd = I_ERI_D2x_S_F2xy_S_C1001001_cd+CDX*I_ERI_D2x_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxy_Px_C1001001_cd = I_ERI_Dxy_S_F2xy_S_C1001001_cd+CDX*I_ERI_Dxy_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxy_Px_C1001001_cd = I_ERI_Dxz_S_F2xy_S_C1001001_cd+CDX*I_ERI_Dxz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxy_Px_C1001001_cd = I_ERI_D2y_S_F2xy_S_C1001001_cd+CDX*I_ERI_D2y_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxy_Px_C1001001_cd = I_ERI_Dyz_S_F2xy_S_C1001001_cd+CDX*I_ERI_Dyz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxy_Px_C1001001_cd = I_ERI_D2z_S_F2xy_S_C1001001_cd+CDX*I_ERI_D2z_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxz_Px_C1001001_cd = I_ERI_D2x_S_F2xz_S_C1001001_cd+CDX*I_ERI_D2x_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxz_Px_C1001001_cd = I_ERI_Dxy_S_F2xz_S_C1001001_cd+CDX*I_ERI_Dxy_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxz_Px_C1001001_cd = I_ERI_Dxz_S_F2xz_S_C1001001_cd+CDX*I_ERI_Dxz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxz_Px_C1001001_cd = I_ERI_D2y_S_F2xz_S_C1001001_cd+CDX*I_ERI_D2y_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxz_Px_C1001001_cd = I_ERI_Dyz_S_F2xz_S_C1001001_cd+CDX*I_ERI_Dyz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxz_Px_C1001001_cd = I_ERI_D2z_S_F2xz_S_C1001001_cd+CDX*I_ERI_D2z_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2y_Px_C1001001_cd = I_ERI_D2x_S_Fx2y_S_C1001001_cd+CDX*I_ERI_D2x_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2y_Px_C1001001_cd = I_ERI_Dxy_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Dxy_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2y_Px_C1001001_cd = I_ERI_Dxz_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Dxz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2y_S_D2y_Px_C1001001_cd = I_ERI_D2y_S_Fx2y_S_C1001001_cd+CDX*I_ERI_D2y_S_D2y_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2y_Px_C1001001_cd = I_ERI_Dyz_S_Fx2y_S_C1001001_cd+CDX*I_ERI_Dyz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2z_S_D2y_Px_C1001001_cd = I_ERI_D2z_S_Fx2y_S_C1001001_cd+CDX*I_ERI_D2z_S_D2y_S_C1001001_cd;
  Double I_ERI_D2x_S_Dyz_Px_C1001001_cd = I_ERI_D2x_S_Fxyz_S_C1001001_cd+CDX*I_ERI_D2x_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dyz_Px_C1001001_cd = I_ERI_Dxy_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Dxy_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dyz_Px_C1001001_cd = I_ERI_Dxz_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Dxz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dyz_Px_C1001001_cd = I_ERI_D2y_S_Fxyz_S_C1001001_cd+CDX*I_ERI_D2y_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dyz_Px_C1001001_cd = I_ERI_Dyz_S_Fxyz_S_C1001001_cd+CDX*I_ERI_Dyz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dyz_Px_C1001001_cd = I_ERI_D2z_S_Fxyz_S_C1001001_cd+CDX*I_ERI_D2z_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2z_Px_C1001001_cd = I_ERI_D2x_S_Fx2z_S_C1001001_cd+CDX*I_ERI_D2x_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2z_Px_C1001001_cd = I_ERI_Dxy_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Dxy_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2z_Px_C1001001_cd = I_ERI_Dxz_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Dxz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2y_S_D2z_Px_C1001001_cd = I_ERI_D2y_S_Fx2z_S_C1001001_cd+CDX*I_ERI_D2y_S_D2z_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2z_Px_C1001001_cd = I_ERI_Dyz_S_Fx2z_S_C1001001_cd+CDX*I_ERI_Dyz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2z_S_D2z_Px_C1001001_cd = I_ERI_D2z_S_Fx2z_S_C1001001_cd+CDX*I_ERI_D2z_S_D2z_S_C1001001_cd;
  Double I_ERI_D2x_S_D2x_Py_C1001001_cd = I_ERI_D2x_S_F2xy_S_C1001001_cd+CDY*I_ERI_D2x_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2x_Py_C1001001_cd = I_ERI_Dxy_S_F2xy_S_C1001001_cd+CDY*I_ERI_Dxy_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2x_Py_C1001001_cd = I_ERI_Dxz_S_F2xy_S_C1001001_cd+CDY*I_ERI_Dxz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2y_S_D2x_Py_C1001001_cd = I_ERI_D2y_S_F2xy_S_C1001001_cd+CDY*I_ERI_D2y_S_D2x_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2x_Py_C1001001_cd = I_ERI_Dyz_S_F2xy_S_C1001001_cd+CDY*I_ERI_Dyz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2z_S_D2x_Py_C1001001_cd = I_ERI_D2z_S_F2xy_S_C1001001_cd+CDY*I_ERI_D2z_S_D2x_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxy_Py_C1001001_cd = I_ERI_D2x_S_Fx2y_S_C1001001_cd+CDY*I_ERI_D2x_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxy_Py_C1001001_cd = I_ERI_Dxy_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Dxy_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxy_Py_C1001001_cd = I_ERI_Dxz_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Dxz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxy_Py_C1001001_cd = I_ERI_D2y_S_Fx2y_S_C1001001_cd+CDY*I_ERI_D2y_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxy_Py_C1001001_cd = I_ERI_Dyz_S_Fx2y_S_C1001001_cd+CDY*I_ERI_Dyz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxy_Py_C1001001_cd = I_ERI_D2z_S_Fx2y_S_C1001001_cd+CDY*I_ERI_D2z_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxz_Py_C1001001_cd = I_ERI_D2x_S_Fxyz_S_C1001001_cd+CDY*I_ERI_D2x_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxz_Py_C1001001_cd = I_ERI_Dxy_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Dxy_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxz_Py_C1001001_cd = I_ERI_Dxz_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Dxz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxz_Py_C1001001_cd = I_ERI_D2y_S_Fxyz_S_C1001001_cd+CDY*I_ERI_D2y_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxz_Py_C1001001_cd = I_ERI_Dyz_S_Fxyz_S_C1001001_cd+CDY*I_ERI_Dyz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxz_Py_C1001001_cd = I_ERI_D2z_S_Fxyz_S_C1001001_cd+CDY*I_ERI_D2z_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2y_Py_C1001001_cd = I_ERI_D2x_S_F3y_S_C1001001_cd+CDY*I_ERI_D2x_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2y_Py_C1001001_cd = I_ERI_Dxy_S_F3y_S_C1001001_cd+CDY*I_ERI_Dxy_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2y_Py_C1001001_cd = I_ERI_Dxz_S_F3y_S_C1001001_cd+CDY*I_ERI_Dxz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2y_S_D2y_Py_C1001001_cd = I_ERI_D2y_S_F3y_S_C1001001_cd+CDY*I_ERI_D2y_S_D2y_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2y_Py_C1001001_cd = I_ERI_Dyz_S_F3y_S_C1001001_cd+CDY*I_ERI_Dyz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2z_S_D2y_Py_C1001001_cd = I_ERI_D2z_S_F3y_S_C1001001_cd+CDY*I_ERI_D2z_S_D2y_S_C1001001_cd;
  Double I_ERI_D2x_S_Dyz_Py_C1001001_cd = I_ERI_D2x_S_F2yz_S_C1001001_cd+CDY*I_ERI_D2x_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dyz_Py_C1001001_cd = I_ERI_Dxy_S_F2yz_S_C1001001_cd+CDY*I_ERI_Dxy_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dyz_Py_C1001001_cd = I_ERI_Dxz_S_F2yz_S_C1001001_cd+CDY*I_ERI_Dxz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dyz_Py_C1001001_cd = I_ERI_D2y_S_F2yz_S_C1001001_cd+CDY*I_ERI_D2y_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dyz_Py_C1001001_cd = I_ERI_Dyz_S_F2yz_S_C1001001_cd+CDY*I_ERI_Dyz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dyz_Py_C1001001_cd = I_ERI_D2z_S_F2yz_S_C1001001_cd+CDY*I_ERI_D2z_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2z_Py_C1001001_cd = I_ERI_D2x_S_Fy2z_S_C1001001_cd+CDY*I_ERI_D2x_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2z_Py_C1001001_cd = I_ERI_Dxy_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Dxy_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2z_Py_C1001001_cd = I_ERI_Dxz_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Dxz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2y_S_D2z_Py_C1001001_cd = I_ERI_D2y_S_Fy2z_S_C1001001_cd+CDY*I_ERI_D2y_S_D2z_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2z_Py_C1001001_cd = I_ERI_Dyz_S_Fy2z_S_C1001001_cd+CDY*I_ERI_Dyz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2z_S_D2z_Py_C1001001_cd = I_ERI_D2z_S_Fy2z_S_C1001001_cd+CDY*I_ERI_D2z_S_D2z_S_C1001001_cd;
  Double I_ERI_D2x_S_D2x_Pz_C1001001_cd = I_ERI_D2x_S_F2xz_S_C1001001_cd+CDZ*I_ERI_D2x_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2x_Pz_C1001001_cd = I_ERI_Dxy_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Dxy_S_D2x_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2x_Pz_C1001001_cd = I_ERI_Dxz_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Dxz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2y_S_D2x_Pz_C1001001_cd = I_ERI_D2y_S_F2xz_S_C1001001_cd+CDZ*I_ERI_D2y_S_D2x_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2x_Pz_C1001001_cd = I_ERI_Dyz_S_F2xz_S_C1001001_cd+CDZ*I_ERI_Dyz_S_D2x_S_C1001001_cd;
  Double I_ERI_D2z_S_D2x_Pz_C1001001_cd = I_ERI_D2z_S_F2xz_S_C1001001_cd+CDZ*I_ERI_D2z_S_D2x_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxy_Pz_C1001001_cd = I_ERI_D2x_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_D2x_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxy_Pz_C1001001_cd = I_ERI_Dxy_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Dxy_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxy_Pz_C1001001_cd = I_ERI_Dxz_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Dxz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxy_Pz_C1001001_cd = I_ERI_D2y_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_D2y_S_Dxy_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxy_Pz_C1001001_cd = I_ERI_Dyz_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_Dyz_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxy_Pz_C1001001_cd = I_ERI_D2z_S_Fxyz_S_C1001001_cd+CDZ*I_ERI_D2z_S_Dxy_S_C1001001_cd;
  Double I_ERI_D2x_S_Dxz_Pz_C1001001_cd = I_ERI_D2x_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_D2x_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dxz_Pz_C1001001_cd = I_ERI_Dxy_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Dxy_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dxz_Pz_C1001001_cd = I_ERI_Dxz_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Dxz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dxz_Pz_C1001001_cd = I_ERI_D2y_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_D2y_S_Dxz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dxz_Pz_C1001001_cd = I_ERI_Dyz_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_Dyz_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dxz_Pz_C1001001_cd = I_ERI_D2z_S_Fx2z_S_C1001001_cd+CDZ*I_ERI_D2z_S_Dxz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2y_Pz_C1001001_cd = I_ERI_D2x_S_F2yz_S_C1001001_cd+CDZ*I_ERI_D2x_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2y_Pz_C1001001_cd = I_ERI_Dxy_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Dxy_S_D2y_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2y_Pz_C1001001_cd = I_ERI_Dxz_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Dxz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2y_S_D2y_Pz_C1001001_cd = I_ERI_D2y_S_F2yz_S_C1001001_cd+CDZ*I_ERI_D2y_S_D2y_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2y_Pz_C1001001_cd = I_ERI_Dyz_S_F2yz_S_C1001001_cd+CDZ*I_ERI_Dyz_S_D2y_S_C1001001_cd;
  Double I_ERI_D2z_S_D2y_Pz_C1001001_cd = I_ERI_D2z_S_F2yz_S_C1001001_cd+CDZ*I_ERI_D2z_S_D2y_S_C1001001_cd;
  Double I_ERI_D2x_S_Dyz_Pz_C1001001_cd = I_ERI_D2x_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_D2x_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxy_S_Dyz_Pz_C1001001_cd = I_ERI_Dxy_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Dxy_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dxz_S_Dyz_Pz_C1001001_cd = I_ERI_Dxz_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Dxz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2y_S_Dyz_Pz_C1001001_cd = I_ERI_D2y_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_D2y_S_Dyz_S_C1001001_cd;
  Double I_ERI_Dyz_S_Dyz_Pz_C1001001_cd = I_ERI_Dyz_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_Dyz_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2z_S_Dyz_Pz_C1001001_cd = I_ERI_D2z_S_Fy2z_S_C1001001_cd+CDZ*I_ERI_D2z_S_Dyz_S_C1001001_cd;
  Double I_ERI_D2x_S_D2z_Pz_C1001001_cd = I_ERI_D2x_S_F3z_S_C1001001_cd+CDZ*I_ERI_D2x_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxy_S_D2z_Pz_C1001001_cd = I_ERI_Dxy_S_F3z_S_C1001001_cd+CDZ*I_ERI_Dxy_S_D2z_S_C1001001_cd;
  Double I_ERI_Dxz_S_D2z_Pz_C1001001_cd = I_ERI_Dxz_S_F3z_S_C1001001_cd+CDZ*I_ERI_Dxz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2y_S_D2z_Pz_C1001001_cd = I_ERI_D2y_S_F3z_S_C1001001_cd+CDZ*I_ERI_D2y_S_D2z_S_C1001001_cd;
  Double I_ERI_Dyz_S_D2z_Pz_C1001001_cd = I_ERI_Dyz_S_F3z_S_C1001001_cd+CDZ*I_ERI_Dyz_S_D2z_S_C1001001_cd;
  Double I_ERI_D2z_S_D2z_Pz_C1001001_cd = I_ERI_D2z_S_F3z_S_C1001001_cd+CDZ*I_ERI_D2z_S_D2z_S_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1000001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_dd
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1000001_dd = I_ERI_Px_S_D2x_S_C1000001_dd+CDX*I_ERI_Px_S_Px_S_C1000001_dd;
  Double I_ERI_Py_S_Px_Px_C1000001_dd = I_ERI_Py_S_D2x_S_C1000001_dd+CDX*I_ERI_Py_S_Px_S_C1000001_dd;
  Double I_ERI_Pz_S_Px_Px_C1000001_dd = I_ERI_Pz_S_D2x_S_C1000001_dd+CDX*I_ERI_Pz_S_Px_S_C1000001_dd;
  Double I_ERI_Px_S_Py_Px_C1000001_dd = I_ERI_Px_S_Dxy_S_C1000001_dd+CDX*I_ERI_Px_S_Py_S_C1000001_dd;
  Double I_ERI_Py_S_Py_Px_C1000001_dd = I_ERI_Py_S_Dxy_S_C1000001_dd+CDX*I_ERI_Py_S_Py_S_C1000001_dd;
  Double I_ERI_Pz_S_Py_Px_C1000001_dd = I_ERI_Pz_S_Dxy_S_C1000001_dd+CDX*I_ERI_Pz_S_Py_S_C1000001_dd;
  Double I_ERI_Px_S_Pz_Px_C1000001_dd = I_ERI_Px_S_Dxz_S_C1000001_dd+CDX*I_ERI_Px_S_Pz_S_C1000001_dd;
  Double I_ERI_Py_S_Pz_Px_C1000001_dd = I_ERI_Py_S_Dxz_S_C1000001_dd+CDX*I_ERI_Py_S_Pz_S_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Px_C1000001_dd = I_ERI_Pz_S_Dxz_S_C1000001_dd+CDX*I_ERI_Pz_S_Pz_S_C1000001_dd;
  Double I_ERI_Px_S_Px_Py_C1000001_dd = I_ERI_Px_S_Dxy_S_C1000001_dd+CDY*I_ERI_Px_S_Px_S_C1000001_dd;
  Double I_ERI_Py_S_Px_Py_C1000001_dd = I_ERI_Py_S_Dxy_S_C1000001_dd+CDY*I_ERI_Py_S_Px_S_C1000001_dd;
  Double I_ERI_Pz_S_Px_Py_C1000001_dd = I_ERI_Pz_S_Dxy_S_C1000001_dd+CDY*I_ERI_Pz_S_Px_S_C1000001_dd;
  Double I_ERI_Px_S_Py_Py_C1000001_dd = I_ERI_Px_S_D2y_S_C1000001_dd+CDY*I_ERI_Px_S_Py_S_C1000001_dd;
  Double I_ERI_Py_S_Py_Py_C1000001_dd = I_ERI_Py_S_D2y_S_C1000001_dd+CDY*I_ERI_Py_S_Py_S_C1000001_dd;
  Double I_ERI_Pz_S_Py_Py_C1000001_dd = I_ERI_Pz_S_D2y_S_C1000001_dd+CDY*I_ERI_Pz_S_Py_S_C1000001_dd;
  Double I_ERI_Px_S_Pz_Py_C1000001_dd = I_ERI_Px_S_Dyz_S_C1000001_dd+CDY*I_ERI_Px_S_Pz_S_C1000001_dd;
  Double I_ERI_Py_S_Pz_Py_C1000001_dd = I_ERI_Py_S_Dyz_S_C1000001_dd+CDY*I_ERI_Py_S_Pz_S_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Py_C1000001_dd = I_ERI_Pz_S_Dyz_S_C1000001_dd+CDY*I_ERI_Pz_S_Pz_S_C1000001_dd;
  Double I_ERI_Px_S_Px_Pz_C1000001_dd = I_ERI_Px_S_Dxz_S_C1000001_dd+CDZ*I_ERI_Px_S_Px_S_C1000001_dd;
  Double I_ERI_Py_S_Px_Pz_C1000001_dd = I_ERI_Py_S_Dxz_S_C1000001_dd+CDZ*I_ERI_Py_S_Px_S_C1000001_dd;
  Double I_ERI_Pz_S_Px_Pz_C1000001_dd = I_ERI_Pz_S_Dxz_S_C1000001_dd+CDZ*I_ERI_Pz_S_Px_S_C1000001_dd;
  Double I_ERI_Px_S_Py_Pz_C1000001_dd = I_ERI_Px_S_Dyz_S_C1000001_dd+CDZ*I_ERI_Px_S_Py_S_C1000001_dd;
  Double I_ERI_Py_S_Py_Pz_C1000001_dd = I_ERI_Py_S_Dyz_S_C1000001_dd+CDZ*I_ERI_Py_S_Py_S_C1000001_dd;
  Double I_ERI_Pz_S_Py_Pz_C1000001_dd = I_ERI_Pz_S_Dyz_S_C1000001_dd+CDZ*I_ERI_Pz_S_Py_S_C1000001_dd;
  Double I_ERI_Px_S_Pz_Pz_C1000001_dd = I_ERI_Px_S_D2z_S_C1000001_dd+CDZ*I_ERI_Px_S_Pz_S_C1000001_dd;
  Double I_ERI_Py_S_Pz_Pz_C1000001_dd = I_ERI_Py_S_D2z_S_C1000001_dd+CDZ*I_ERI_Py_S_Pz_S_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Pz_C1000001_dd = I_ERI_Pz_S_D2z_S_C1000001_dd+CDZ*I_ERI_Pz_S_Pz_S_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1000001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_dd
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1000001_dd = I_ERI_Px_S_F3x_S_C1000001_dd+CDX*I_ERI_Px_S_D2x_S_C1000001_dd;
  Double I_ERI_Py_S_D2x_Px_C1000001_dd = I_ERI_Py_S_F3x_S_C1000001_dd+CDX*I_ERI_Py_S_D2x_S_C1000001_dd;
  Double I_ERI_Pz_S_D2x_Px_C1000001_dd = I_ERI_Pz_S_F3x_S_C1000001_dd+CDX*I_ERI_Pz_S_D2x_S_C1000001_dd;
  Double I_ERI_Px_S_Dxy_Px_C1000001_dd = I_ERI_Px_S_F2xy_S_C1000001_dd+CDX*I_ERI_Px_S_Dxy_S_C1000001_dd;
  Double I_ERI_Py_S_Dxy_Px_C1000001_dd = I_ERI_Py_S_F2xy_S_C1000001_dd+CDX*I_ERI_Py_S_Dxy_S_C1000001_dd;
  Double I_ERI_Pz_S_Dxy_Px_C1000001_dd = I_ERI_Pz_S_F2xy_S_C1000001_dd+CDX*I_ERI_Pz_S_Dxy_S_C1000001_dd;
  Double I_ERI_Px_S_Dxz_Px_C1000001_dd = I_ERI_Px_S_F2xz_S_C1000001_dd+CDX*I_ERI_Px_S_Dxz_S_C1000001_dd;
  Double I_ERI_Py_S_Dxz_Px_C1000001_dd = I_ERI_Py_S_F2xz_S_C1000001_dd+CDX*I_ERI_Py_S_Dxz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dxz_Px_C1000001_dd = I_ERI_Pz_S_F2xz_S_C1000001_dd+CDX*I_ERI_Pz_S_Dxz_S_C1000001_dd;
  Double I_ERI_Px_S_D2y_Px_C1000001_dd = I_ERI_Px_S_Fx2y_S_C1000001_dd+CDX*I_ERI_Px_S_D2y_S_C1000001_dd;
  Double I_ERI_Py_S_D2y_Px_C1000001_dd = I_ERI_Py_S_Fx2y_S_C1000001_dd+CDX*I_ERI_Py_S_D2y_S_C1000001_dd;
  Double I_ERI_Pz_S_D2y_Px_C1000001_dd = I_ERI_Pz_S_Fx2y_S_C1000001_dd+CDX*I_ERI_Pz_S_D2y_S_C1000001_dd;
  Double I_ERI_Px_S_Dyz_Px_C1000001_dd = I_ERI_Px_S_Fxyz_S_C1000001_dd+CDX*I_ERI_Px_S_Dyz_S_C1000001_dd;
  Double I_ERI_Py_S_Dyz_Px_C1000001_dd = I_ERI_Py_S_Fxyz_S_C1000001_dd+CDX*I_ERI_Py_S_Dyz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dyz_Px_C1000001_dd = I_ERI_Pz_S_Fxyz_S_C1000001_dd+CDX*I_ERI_Pz_S_Dyz_S_C1000001_dd;
  Double I_ERI_Px_S_D2z_Px_C1000001_dd = I_ERI_Px_S_Fx2z_S_C1000001_dd+CDX*I_ERI_Px_S_D2z_S_C1000001_dd;
  Double I_ERI_Py_S_D2z_Px_C1000001_dd = I_ERI_Py_S_Fx2z_S_C1000001_dd+CDX*I_ERI_Py_S_D2z_S_C1000001_dd;
  Double I_ERI_Pz_S_D2z_Px_C1000001_dd = I_ERI_Pz_S_Fx2z_S_C1000001_dd+CDX*I_ERI_Pz_S_D2z_S_C1000001_dd;
  Double I_ERI_Px_S_Dxy_Py_C1000001_dd = I_ERI_Px_S_Fx2y_S_C1000001_dd+CDY*I_ERI_Px_S_Dxy_S_C1000001_dd;
  Double I_ERI_Py_S_Dxy_Py_C1000001_dd = I_ERI_Py_S_Fx2y_S_C1000001_dd+CDY*I_ERI_Py_S_Dxy_S_C1000001_dd;
  Double I_ERI_Pz_S_Dxy_Py_C1000001_dd = I_ERI_Pz_S_Fx2y_S_C1000001_dd+CDY*I_ERI_Pz_S_Dxy_S_C1000001_dd;
  Double I_ERI_Px_S_Dxz_Py_C1000001_dd = I_ERI_Px_S_Fxyz_S_C1000001_dd+CDY*I_ERI_Px_S_Dxz_S_C1000001_dd;
  Double I_ERI_Py_S_Dxz_Py_C1000001_dd = I_ERI_Py_S_Fxyz_S_C1000001_dd+CDY*I_ERI_Py_S_Dxz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dxz_Py_C1000001_dd = I_ERI_Pz_S_Fxyz_S_C1000001_dd+CDY*I_ERI_Pz_S_Dxz_S_C1000001_dd;
  Double I_ERI_Px_S_D2y_Py_C1000001_dd = I_ERI_Px_S_F3y_S_C1000001_dd+CDY*I_ERI_Px_S_D2y_S_C1000001_dd;
  Double I_ERI_Py_S_D2y_Py_C1000001_dd = I_ERI_Py_S_F3y_S_C1000001_dd+CDY*I_ERI_Py_S_D2y_S_C1000001_dd;
  Double I_ERI_Pz_S_D2y_Py_C1000001_dd = I_ERI_Pz_S_F3y_S_C1000001_dd+CDY*I_ERI_Pz_S_D2y_S_C1000001_dd;
  Double I_ERI_Px_S_Dyz_Py_C1000001_dd = I_ERI_Px_S_F2yz_S_C1000001_dd+CDY*I_ERI_Px_S_Dyz_S_C1000001_dd;
  Double I_ERI_Py_S_Dyz_Py_C1000001_dd = I_ERI_Py_S_F2yz_S_C1000001_dd+CDY*I_ERI_Py_S_Dyz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dyz_Py_C1000001_dd = I_ERI_Pz_S_F2yz_S_C1000001_dd+CDY*I_ERI_Pz_S_Dyz_S_C1000001_dd;
  Double I_ERI_Px_S_D2z_Py_C1000001_dd = I_ERI_Px_S_Fy2z_S_C1000001_dd+CDY*I_ERI_Px_S_D2z_S_C1000001_dd;
  Double I_ERI_Py_S_D2z_Py_C1000001_dd = I_ERI_Py_S_Fy2z_S_C1000001_dd+CDY*I_ERI_Py_S_D2z_S_C1000001_dd;
  Double I_ERI_Pz_S_D2z_Py_C1000001_dd = I_ERI_Pz_S_Fy2z_S_C1000001_dd+CDY*I_ERI_Pz_S_D2z_S_C1000001_dd;
  Double I_ERI_Px_S_Dxz_Pz_C1000001_dd = I_ERI_Px_S_Fx2z_S_C1000001_dd+CDZ*I_ERI_Px_S_Dxz_S_C1000001_dd;
  Double I_ERI_Py_S_Dxz_Pz_C1000001_dd = I_ERI_Py_S_Fx2z_S_C1000001_dd+CDZ*I_ERI_Py_S_Dxz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dxz_Pz_C1000001_dd = I_ERI_Pz_S_Fx2z_S_C1000001_dd+CDZ*I_ERI_Pz_S_Dxz_S_C1000001_dd;
  Double I_ERI_Px_S_Dyz_Pz_C1000001_dd = I_ERI_Px_S_Fy2z_S_C1000001_dd+CDZ*I_ERI_Px_S_Dyz_S_C1000001_dd;
  Double I_ERI_Py_S_Dyz_Pz_C1000001_dd = I_ERI_Py_S_Fy2z_S_C1000001_dd+CDZ*I_ERI_Py_S_Dyz_S_C1000001_dd;
  Double I_ERI_Pz_S_Dyz_Pz_C1000001_dd = I_ERI_Pz_S_Fy2z_S_C1000001_dd+CDZ*I_ERI_Pz_S_Dyz_S_C1000001_dd;
  Double I_ERI_Px_S_D2z_Pz_C1000001_dd = I_ERI_Px_S_F3z_S_C1000001_dd+CDZ*I_ERI_Px_S_D2z_S_C1000001_dd;
  Double I_ERI_Py_S_D2z_Pz_C1000001_dd = I_ERI_Py_S_F3z_S_C1000001_dd+CDZ*I_ERI_Py_S_D2z_S_C1000001_dd;
  Double I_ERI_Pz_S_D2z_Pz_C1000001_dd = I_ERI_Pz_S_F3z_S_C1000001_dd+CDZ*I_ERI_Pz_S_D2z_S_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000001_dd
   ************************************************************/
  Double I_ERI_Px_S_Px_D2x_C1000001_dd = I_ERI_Px_S_D2x_Px_C1000001_dd+CDX*I_ERI_Px_S_Px_Px_C1000001_dd;
  Double I_ERI_Py_S_Px_D2x_C1000001_dd = I_ERI_Py_S_D2x_Px_C1000001_dd+CDX*I_ERI_Py_S_Px_Px_C1000001_dd;
  Double I_ERI_Pz_S_Px_D2x_C1000001_dd = I_ERI_Pz_S_D2x_Px_C1000001_dd+CDX*I_ERI_Pz_S_Px_Px_C1000001_dd;
  Double I_ERI_Px_S_Py_D2x_C1000001_dd = I_ERI_Px_S_Dxy_Px_C1000001_dd+CDX*I_ERI_Px_S_Py_Px_C1000001_dd;
  Double I_ERI_Py_S_Py_D2x_C1000001_dd = I_ERI_Py_S_Dxy_Px_C1000001_dd+CDX*I_ERI_Py_S_Py_Px_C1000001_dd;
  Double I_ERI_Pz_S_Py_D2x_C1000001_dd = I_ERI_Pz_S_Dxy_Px_C1000001_dd+CDX*I_ERI_Pz_S_Py_Px_C1000001_dd;
  Double I_ERI_Px_S_Pz_D2x_C1000001_dd = I_ERI_Px_S_Dxz_Px_C1000001_dd+CDX*I_ERI_Px_S_Pz_Px_C1000001_dd;
  Double I_ERI_Py_S_Pz_D2x_C1000001_dd = I_ERI_Py_S_Dxz_Px_C1000001_dd+CDX*I_ERI_Py_S_Pz_Px_C1000001_dd;
  Double I_ERI_Pz_S_Pz_D2x_C1000001_dd = I_ERI_Pz_S_Dxz_Px_C1000001_dd+CDX*I_ERI_Pz_S_Pz_Px_C1000001_dd;
  Double I_ERI_Px_S_Px_Dxy_C1000001_dd = I_ERI_Px_S_Dxy_Px_C1000001_dd+CDY*I_ERI_Px_S_Px_Px_C1000001_dd;
  Double I_ERI_Py_S_Px_Dxy_C1000001_dd = I_ERI_Py_S_Dxy_Px_C1000001_dd+CDY*I_ERI_Py_S_Px_Px_C1000001_dd;
  Double I_ERI_Pz_S_Px_Dxy_C1000001_dd = I_ERI_Pz_S_Dxy_Px_C1000001_dd+CDY*I_ERI_Pz_S_Px_Px_C1000001_dd;
  Double I_ERI_Px_S_Py_Dxy_C1000001_dd = I_ERI_Px_S_D2y_Px_C1000001_dd+CDY*I_ERI_Px_S_Py_Px_C1000001_dd;
  Double I_ERI_Py_S_Py_Dxy_C1000001_dd = I_ERI_Py_S_D2y_Px_C1000001_dd+CDY*I_ERI_Py_S_Py_Px_C1000001_dd;
  Double I_ERI_Pz_S_Py_Dxy_C1000001_dd = I_ERI_Pz_S_D2y_Px_C1000001_dd+CDY*I_ERI_Pz_S_Py_Px_C1000001_dd;
  Double I_ERI_Px_S_Pz_Dxy_C1000001_dd = I_ERI_Px_S_Dyz_Px_C1000001_dd+CDY*I_ERI_Px_S_Pz_Px_C1000001_dd;
  Double I_ERI_Py_S_Pz_Dxy_C1000001_dd = I_ERI_Py_S_Dyz_Px_C1000001_dd+CDY*I_ERI_Py_S_Pz_Px_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Dxy_C1000001_dd = I_ERI_Pz_S_Dyz_Px_C1000001_dd+CDY*I_ERI_Pz_S_Pz_Px_C1000001_dd;
  Double I_ERI_Px_S_Px_Dxz_C1000001_dd = I_ERI_Px_S_Dxz_Px_C1000001_dd+CDZ*I_ERI_Px_S_Px_Px_C1000001_dd;
  Double I_ERI_Py_S_Px_Dxz_C1000001_dd = I_ERI_Py_S_Dxz_Px_C1000001_dd+CDZ*I_ERI_Py_S_Px_Px_C1000001_dd;
  Double I_ERI_Pz_S_Px_Dxz_C1000001_dd = I_ERI_Pz_S_Dxz_Px_C1000001_dd+CDZ*I_ERI_Pz_S_Px_Px_C1000001_dd;
  Double I_ERI_Px_S_Py_Dxz_C1000001_dd = I_ERI_Px_S_Dyz_Px_C1000001_dd+CDZ*I_ERI_Px_S_Py_Px_C1000001_dd;
  Double I_ERI_Py_S_Py_Dxz_C1000001_dd = I_ERI_Py_S_Dyz_Px_C1000001_dd+CDZ*I_ERI_Py_S_Py_Px_C1000001_dd;
  Double I_ERI_Pz_S_Py_Dxz_C1000001_dd = I_ERI_Pz_S_Dyz_Px_C1000001_dd+CDZ*I_ERI_Pz_S_Py_Px_C1000001_dd;
  Double I_ERI_Px_S_Pz_Dxz_C1000001_dd = I_ERI_Px_S_D2z_Px_C1000001_dd+CDZ*I_ERI_Px_S_Pz_Px_C1000001_dd;
  Double I_ERI_Py_S_Pz_Dxz_C1000001_dd = I_ERI_Py_S_D2z_Px_C1000001_dd+CDZ*I_ERI_Py_S_Pz_Px_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Dxz_C1000001_dd = I_ERI_Pz_S_D2z_Px_C1000001_dd+CDZ*I_ERI_Pz_S_Pz_Px_C1000001_dd;
  Double I_ERI_Px_S_Px_D2y_C1000001_dd = I_ERI_Px_S_Dxy_Py_C1000001_dd+CDY*I_ERI_Px_S_Px_Py_C1000001_dd;
  Double I_ERI_Py_S_Px_D2y_C1000001_dd = I_ERI_Py_S_Dxy_Py_C1000001_dd+CDY*I_ERI_Py_S_Px_Py_C1000001_dd;
  Double I_ERI_Pz_S_Px_D2y_C1000001_dd = I_ERI_Pz_S_Dxy_Py_C1000001_dd+CDY*I_ERI_Pz_S_Px_Py_C1000001_dd;
  Double I_ERI_Px_S_Py_D2y_C1000001_dd = I_ERI_Px_S_D2y_Py_C1000001_dd+CDY*I_ERI_Px_S_Py_Py_C1000001_dd;
  Double I_ERI_Py_S_Py_D2y_C1000001_dd = I_ERI_Py_S_D2y_Py_C1000001_dd+CDY*I_ERI_Py_S_Py_Py_C1000001_dd;
  Double I_ERI_Pz_S_Py_D2y_C1000001_dd = I_ERI_Pz_S_D2y_Py_C1000001_dd+CDY*I_ERI_Pz_S_Py_Py_C1000001_dd;
  Double I_ERI_Px_S_Pz_D2y_C1000001_dd = I_ERI_Px_S_Dyz_Py_C1000001_dd+CDY*I_ERI_Px_S_Pz_Py_C1000001_dd;
  Double I_ERI_Py_S_Pz_D2y_C1000001_dd = I_ERI_Py_S_Dyz_Py_C1000001_dd+CDY*I_ERI_Py_S_Pz_Py_C1000001_dd;
  Double I_ERI_Pz_S_Pz_D2y_C1000001_dd = I_ERI_Pz_S_Dyz_Py_C1000001_dd+CDY*I_ERI_Pz_S_Pz_Py_C1000001_dd;
  Double I_ERI_Px_S_Px_Dyz_C1000001_dd = I_ERI_Px_S_Dxz_Py_C1000001_dd+CDZ*I_ERI_Px_S_Px_Py_C1000001_dd;
  Double I_ERI_Py_S_Px_Dyz_C1000001_dd = I_ERI_Py_S_Dxz_Py_C1000001_dd+CDZ*I_ERI_Py_S_Px_Py_C1000001_dd;
  Double I_ERI_Pz_S_Px_Dyz_C1000001_dd = I_ERI_Pz_S_Dxz_Py_C1000001_dd+CDZ*I_ERI_Pz_S_Px_Py_C1000001_dd;
  Double I_ERI_Px_S_Py_Dyz_C1000001_dd = I_ERI_Px_S_Dyz_Py_C1000001_dd+CDZ*I_ERI_Px_S_Py_Py_C1000001_dd;
  Double I_ERI_Py_S_Py_Dyz_C1000001_dd = I_ERI_Py_S_Dyz_Py_C1000001_dd+CDZ*I_ERI_Py_S_Py_Py_C1000001_dd;
  Double I_ERI_Pz_S_Py_Dyz_C1000001_dd = I_ERI_Pz_S_Dyz_Py_C1000001_dd+CDZ*I_ERI_Pz_S_Py_Py_C1000001_dd;
  Double I_ERI_Px_S_Pz_Dyz_C1000001_dd = I_ERI_Px_S_D2z_Py_C1000001_dd+CDZ*I_ERI_Px_S_Pz_Py_C1000001_dd;
  Double I_ERI_Py_S_Pz_Dyz_C1000001_dd = I_ERI_Py_S_D2z_Py_C1000001_dd+CDZ*I_ERI_Py_S_Pz_Py_C1000001_dd;
  Double I_ERI_Pz_S_Pz_Dyz_C1000001_dd = I_ERI_Pz_S_D2z_Py_C1000001_dd+CDZ*I_ERI_Pz_S_Pz_Py_C1000001_dd;
  Double I_ERI_Px_S_Px_D2z_C1000001_dd = I_ERI_Px_S_Dxz_Pz_C1000001_dd+CDZ*I_ERI_Px_S_Px_Pz_C1000001_dd;
  Double I_ERI_Py_S_Px_D2z_C1000001_dd = I_ERI_Py_S_Dxz_Pz_C1000001_dd+CDZ*I_ERI_Py_S_Px_Pz_C1000001_dd;
  Double I_ERI_Pz_S_Px_D2z_C1000001_dd = I_ERI_Pz_S_Dxz_Pz_C1000001_dd+CDZ*I_ERI_Pz_S_Px_Pz_C1000001_dd;
  Double I_ERI_Px_S_Py_D2z_C1000001_dd = I_ERI_Px_S_Dyz_Pz_C1000001_dd+CDZ*I_ERI_Px_S_Py_Pz_C1000001_dd;
  Double I_ERI_Py_S_Py_D2z_C1000001_dd = I_ERI_Py_S_Dyz_Pz_C1000001_dd+CDZ*I_ERI_Py_S_Py_Pz_C1000001_dd;
  Double I_ERI_Pz_S_Py_D2z_C1000001_dd = I_ERI_Pz_S_Dyz_Pz_C1000001_dd+CDZ*I_ERI_Pz_S_Py_Pz_C1000001_dd;
  Double I_ERI_Px_S_Pz_D2z_C1000001_dd = I_ERI_Px_S_D2z_Pz_C1000001_dd+CDZ*I_ERI_Px_S_Pz_Pz_C1000001_dd;
  Double I_ERI_Py_S_Pz_D2z_C1000001_dd = I_ERI_Py_S_D2z_Pz_C1000001_dd+CDZ*I_ERI_Py_S_Pz_Pz_C1000001_dd;
  Double I_ERI_Pz_S_Pz_D2z_C1000001_dd = I_ERI_Pz_S_D2z_Pz_C1000001_dd+CDZ*I_ERI_Pz_S_Pz_Pz_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_dd
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001001_dd = I_ERI_Px_S_D2x_S_C1001001_dd+CDX*I_ERI_Px_S_Px_S_C1001001_dd;
  Double I_ERI_Py_S_Px_Px_C1001001_dd = I_ERI_Py_S_D2x_S_C1001001_dd+CDX*I_ERI_Py_S_Px_S_C1001001_dd;
  Double I_ERI_Pz_S_Px_Px_C1001001_dd = I_ERI_Pz_S_D2x_S_C1001001_dd+CDX*I_ERI_Pz_S_Px_S_C1001001_dd;
  Double I_ERI_Px_S_Py_Px_C1001001_dd = I_ERI_Px_S_Dxy_S_C1001001_dd+CDX*I_ERI_Px_S_Py_S_C1001001_dd;
  Double I_ERI_Py_S_Py_Px_C1001001_dd = I_ERI_Py_S_Dxy_S_C1001001_dd+CDX*I_ERI_Py_S_Py_S_C1001001_dd;
  Double I_ERI_Pz_S_Py_Px_C1001001_dd = I_ERI_Pz_S_Dxy_S_C1001001_dd+CDX*I_ERI_Pz_S_Py_S_C1001001_dd;
  Double I_ERI_Px_S_Pz_Px_C1001001_dd = I_ERI_Px_S_Dxz_S_C1001001_dd+CDX*I_ERI_Px_S_Pz_S_C1001001_dd;
  Double I_ERI_Py_S_Pz_Px_C1001001_dd = I_ERI_Py_S_Dxz_S_C1001001_dd+CDX*I_ERI_Py_S_Pz_S_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Px_C1001001_dd = I_ERI_Pz_S_Dxz_S_C1001001_dd+CDX*I_ERI_Pz_S_Pz_S_C1001001_dd;
  Double I_ERI_Px_S_Px_Py_C1001001_dd = I_ERI_Px_S_Dxy_S_C1001001_dd+CDY*I_ERI_Px_S_Px_S_C1001001_dd;
  Double I_ERI_Py_S_Px_Py_C1001001_dd = I_ERI_Py_S_Dxy_S_C1001001_dd+CDY*I_ERI_Py_S_Px_S_C1001001_dd;
  Double I_ERI_Pz_S_Px_Py_C1001001_dd = I_ERI_Pz_S_Dxy_S_C1001001_dd+CDY*I_ERI_Pz_S_Px_S_C1001001_dd;
  Double I_ERI_Px_S_Py_Py_C1001001_dd = I_ERI_Px_S_D2y_S_C1001001_dd+CDY*I_ERI_Px_S_Py_S_C1001001_dd;
  Double I_ERI_Py_S_Py_Py_C1001001_dd = I_ERI_Py_S_D2y_S_C1001001_dd+CDY*I_ERI_Py_S_Py_S_C1001001_dd;
  Double I_ERI_Pz_S_Py_Py_C1001001_dd = I_ERI_Pz_S_D2y_S_C1001001_dd+CDY*I_ERI_Pz_S_Py_S_C1001001_dd;
  Double I_ERI_Px_S_Pz_Py_C1001001_dd = I_ERI_Px_S_Dyz_S_C1001001_dd+CDY*I_ERI_Px_S_Pz_S_C1001001_dd;
  Double I_ERI_Py_S_Pz_Py_C1001001_dd = I_ERI_Py_S_Dyz_S_C1001001_dd+CDY*I_ERI_Py_S_Pz_S_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Py_C1001001_dd = I_ERI_Pz_S_Dyz_S_C1001001_dd+CDY*I_ERI_Pz_S_Pz_S_C1001001_dd;
  Double I_ERI_Px_S_Px_Pz_C1001001_dd = I_ERI_Px_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Px_S_Px_S_C1001001_dd;
  Double I_ERI_Py_S_Px_Pz_C1001001_dd = I_ERI_Py_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Py_S_Px_S_C1001001_dd;
  Double I_ERI_Pz_S_Px_Pz_C1001001_dd = I_ERI_Pz_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Pz_S_Px_S_C1001001_dd;
  Double I_ERI_Px_S_Py_Pz_C1001001_dd = I_ERI_Px_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Px_S_Py_S_C1001001_dd;
  Double I_ERI_Py_S_Py_Pz_C1001001_dd = I_ERI_Py_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Py_S_Py_S_C1001001_dd;
  Double I_ERI_Pz_S_Py_Pz_C1001001_dd = I_ERI_Pz_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Pz_S_Py_S_C1001001_dd;
  Double I_ERI_Px_S_Pz_Pz_C1001001_dd = I_ERI_Px_S_D2z_S_C1001001_dd+CDZ*I_ERI_Px_S_Pz_S_C1001001_dd;
  Double I_ERI_Py_S_Pz_Pz_C1001001_dd = I_ERI_Py_S_D2z_S_C1001001_dd+CDZ*I_ERI_Py_S_Pz_S_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Pz_C1001001_dd = I_ERI_Pz_S_D2z_S_C1001001_dd+CDZ*I_ERI_Pz_S_Pz_S_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_dd
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_dd
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001_dd = I_ERI_D2x_S_D2x_S_C1001001_dd+CDX*I_ERI_D2x_S_Px_S_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Px_C1001001_dd = I_ERI_Dxy_S_D2x_S_C1001001_dd+CDX*I_ERI_Dxy_S_Px_S_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Px_C1001001_dd = I_ERI_Dxz_S_D2x_S_C1001001_dd+CDX*I_ERI_Dxz_S_Px_S_C1001001_dd;
  Double I_ERI_D2y_S_Px_Px_C1001001_dd = I_ERI_D2y_S_D2x_S_C1001001_dd+CDX*I_ERI_D2y_S_Px_S_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Px_C1001001_dd = I_ERI_Dyz_S_D2x_S_C1001001_dd+CDX*I_ERI_Dyz_S_Px_S_C1001001_dd;
  Double I_ERI_D2z_S_Px_Px_C1001001_dd = I_ERI_D2z_S_D2x_S_C1001001_dd+CDX*I_ERI_D2z_S_Px_S_C1001001_dd;
  Double I_ERI_D2x_S_Py_Px_C1001001_dd = I_ERI_D2x_S_Dxy_S_C1001001_dd+CDX*I_ERI_D2x_S_Py_S_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Px_C1001001_dd = I_ERI_Dxy_S_Dxy_S_C1001001_dd+CDX*I_ERI_Dxy_S_Py_S_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Px_C1001001_dd = I_ERI_Dxz_S_Dxy_S_C1001001_dd+CDX*I_ERI_Dxz_S_Py_S_C1001001_dd;
  Double I_ERI_D2y_S_Py_Px_C1001001_dd = I_ERI_D2y_S_Dxy_S_C1001001_dd+CDX*I_ERI_D2y_S_Py_S_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Px_C1001001_dd = I_ERI_Dyz_S_Dxy_S_C1001001_dd+CDX*I_ERI_Dyz_S_Py_S_C1001001_dd;
  Double I_ERI_D2z_S_Py_Px_C1001001_dd = I_ERI_D2z_S_Dxy_S_C1001001_dd+CDX*I_ERI_D2z_S_Py_S_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Px_C1001001_dd = I_ERI_D2x_S_Dxz_S_C1001001_dd+CDX*I_ERI_D2x_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Px_C1001001_dd = I_ERI_Dxy_S_Dxz_S_C1001001_dd+CDX*I_ERI_Dxy_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Px_C1001001_dd = I_ERI_Dxz_S_Dxz_S_C1001001_dd+CDX*I_ERI_Dxz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Px_C1001001_dd = I_ERI_D2y_S_Dxz_S_C1001001_dd+CDX*I_ERI_D2y_S_Pz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Px_C1001001_dd = I_ERI_Dyz_S_Dxz_S_C1001001_dd+CDX*I_ERI_Dyz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Px_C1001001_dd = I_ERI_D2z_S_Dxz_S_C1001001_dd+CDX*I_ERI_D2z_S_Pz_S_C1001001_dd;
  Double I_ERI_D2x_S_Px_Py_C1001001_dd = I_ERI_D2x_S_Dxy_S_C1001001_dd+CDY*I_ERI_D2x_S_Px_S_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Py_C1001001_dd = I_ERI_Dxy_S_Dxy_S_C1001001_dd+CDY*I_ERI_Dxy_S_Px_S_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Py_C1001001_dd = I_ERI_Dxz_S_Dxy_S_C1001001_dd+CDY*I_ERI_Dxz_S_Px_S_C1001001_dd;
  Double I_ERI_D2y_S_Px_Py_C1001001_dd = I_ERI_D2y_S_Dxy_S_C1001001_dd+CDY*I_ERI_D2y_S_Px_S_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Py_C1001001_dd = I_ERI_Dyz_S_Dxy_S_C1001001_dd+CDY*I_ERI_Dyz_S_Px_S_C1001001_dd;
  Double I_ERI_D2z_S_Px_Py_C1001001_dd = I_ERI_D2z_S_Dxy_S_C1001001_dd+CDY*I_ERI_D2z_S_Px_S_C1001001_dd;
  Double I_ERI_D2x_S_Py_Py_C1001001_dd = I_ERI_D2x_S_D2y_S_C1001001_dd+CDY*I_ERI_D2x_S_Py_S_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Py_C1001001_dd = I_ERI_Dxy_S_D2y_S_C1001001_dd+CDY*I_ERI_Dxy_S_Py_S_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Py_C1001001_dd = I_ERI_Dxz_S_D2y_S_C1001001_dd+CDY*I_ERI_Dxz_S_Py_S_C1001001_dd;
  Double I_ERI_D2y_S_Py_Py_C1001001_dd = I_ERI_D2y_S_D2y_S_C1001001_dd+CDY*I_ERI_D2y_S_Py_S_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Py_C1001001_dd = I_ERI_Dyz_S_D2y_S_C1001001_dd+CDY*I_ERI_Dyz_S_Py_S_C1001001_dd;
  Double I_ERI_D2z_S_Py_Py_C1001001_dd = I_ERI_D2z_S_D2y_S_C1001001_dd+CDY*I_ERI_D2z_S_Py_S_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Py_C1001001_dd = I_ERI_D2x_S_Dyz_S_C1001001_dd+CDY*I_ERI_D2x_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Py_C1001001_dd = I_ERI_Dxy_S_Dyz_S_C1001001_dd+CDY*I_ERI_Dxy_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Py_C1001001_dd = I_ERI_Dxz_S_Dyz_S_C1001001_dd+CDY*I_ERI_Dxz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Py_C1001001_dd = I_ERI_D2y_S_Dyz_S_C1001001_dd+CDY*I_ERI_D2y_S_Pz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Py_C1001001_dd = I_ERI_Dyz_S_Dyz_S_C1001001_dd+CDY*I_ERI_Dyz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Py_C1001001_dd = I_ERI_D2z_S_Dyz_S_C1001001_dd+CDY*I_ERI_D2z_S_Pz_S_C1001001_dd;
  Double I_ERI_D2x_S_Px_Pz_C1001001_dd = I_ERI_D2x_S_Dxz_S_C1001001_dd+CDZ*I_ERI_D2x_S_Px_S_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Pz_C1001001_dd = I_ERI_Dxy_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Dxy_S_Px_S_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Pz_C1001001_dd = I_ERI_Dxz_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Dxz_S_Px_S_C1001001_dd;
  Double I_ERI_D2y_S_Px_Pz_C1001001_dd = I_ERI_D2y_S_Dxz_S_C1001001_dd+CDZ*I_ERI_D2y_S_Px_S_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Pz_C1001001_dd = I_ERI_Dyz_S_Dxz_S_C1001001_dd+CDZ*I_ERI_Dyz_S_Px_S_C1001001_dd;
  Double I_ERI_D2z_S_Px_Pz_C1001001_dd = I_ERI_D2z_S_Dxz_S_C1001001_dd+CDZ*I_ERI_D2z_S_Px_S_C1001001_dd;
  Double I_ERI_D2x_S_Py_Pz_C1001001_dd = I_ERI_D2x_S_Dyz_S_C1001001_dd+CDZ*I_ERI_D2x_S_Py_S_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Pz_C1001001_dd = I_ERI_Dxy_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Dxy_S_Py_S_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Pz_C1001001_dd = I_ERI_Dxz_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Dxz_S_Py_S_C1001001_dd;
  Double I_ERI_D2y_S_Py_Pz_C1001001_dd = I_ERI_D2y_S_Dyz_S_C1001001_dd+CDZ*I_ERI_D2y_S_Py_S_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Pz_C1001001_dd = I_ERI_Dyz_S_Dyz_S_C1001001_dd+CDZ*I_ERI_Dyz_S_Py_S_C1001001_dd;
  Double I_ERI_D2z_S_Py_Pz_C1001001_dd = I_ERI_D2z_S_Dyz_S_C1001001_dd+CDZ*I_ERI_D2z_S_Py_S_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Pz_C1001001_dd = I_ERI_D2x_S_D2z_S_C1001001_dd+CDZ*I_ERI_D2x_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001_dd = I_ERI_Dxy_S_D2z_S_C1001001_dd+CDZ*I_ERI_Dxy_S_Pz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001_dd = I_ERI_Dxz_S_D2z_S_C1001001_dd+CDZ*I_ERI_Dxz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Pz_C1001001_dd = I_ERI_D2y_S_D2z_S_C1001001_dd+CDZ*I_ERI_D2y_S_Pz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001_dd = I_ERI_Dyz_S_D2z_S_C1001001_dd+CDZ*I_ERI_Dyz_S_Pz_S_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Pz_C1001001_dd = I_ERI_D2z_S_D2z_S_C1001001_dd+CDZ*I_ERI_D2z_S_Pz_S_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001_dd
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1001001_dd = I_ERI_Px_S_F3x_S_C1001001_dd+CDX*I_ERI_Px_S_D2x_S_C1001001_dd;
  Double I_ERI_Py_S_D2x_Px_C1001001_dd = I_ERI_Py_S_F3x_S_C1001001_dd+CDX*I_ERI_Py_S_D2x_S_C1001001_dd;
  Double I_ERI_Pz_S_D2x_Px_C1001001_dd = I_ERI_Pz_S_F3x_S_C1001001_dd+CDX*I_ERI_Pz_S_D2x_S_C1001001_dd;
  Double I_ERI_Px_S_Dxy_Px_C1001001_dd = I_ERI_Px_S_F2xy_S_C1001001_dd+CDX*I_ERI_Px_S_Dxy_S_C1001001_dd;
  Double I_ERI_Py_S_Dxy_Px_C1001001_dd = I_ERI_Py_S_F2xy_S_C1001001_dd+CDX*I_ERI_Py_S_Dxy_S_C1001001_dd;
  Double I_ERI_Pz_S_Dxy_Px_C1001001_dd = I_ERI_Pz_S_F2xy_S_C1001001_dd+CDX*I_ERI_Pz_S_Dxy_S_C1001001_dd;
  Double I_ERI_Px_S_Dxz_Px_C1001001_dd = I_ERI_Px_S_F2xz_S_C1001001_dd+CDX*I_ERI_Px_S_Dxz_S_C1001001_dd;
  Double I_ERI_Py_S_Dxz_Px_C1001001_dd = I_ERI_Py_S_F2xz_S_C1001001_dd+CDX*I_ERI_Py_S_Dxz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dxz_Px_C1001001_dd = I_ERI_Pz_S_F2xz_S_C1001001_dd+CDX*I_ERI_Pz_S_Dxz_S_C1001001_dd;
  Double I_ERI_Px_S_D2y_Px_C1001001_dd = I_ERI_Px_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Px_S_D2y_S_C1001001_dd;
  Double I_ERI_Py_S_D2y_Px_C1001001_dd = I_ERI_Py_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Py_S_D2y_S_C1001001_dd;
  Double I_ERI_Pz_S_D2y_Px_C1001001_dd = I_ERI_Pz_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Pz_S_D2y_S_C1001001_dd;
  Double I_ERI_Px_S_Dyz_Px_C1001001_dd = I_ERI_Px_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Px_S_Dyz_S_C1001001_dd;
  Double I_ERI_Py_S_Dyz_Px_C1001001_dd = I_ERI_Py_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Py_S_Dyz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dyz_Px_C1001001_dd = I_ERI_Pz_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Pz_S_Dyz_S_C1001001_dd;
  Double I_ERI_Px_S_D2z_Px_C1001001_dd = I_ERI_Px_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Px_S_D2z_S_C1001001_dd;
  Double I_ERI_Py_S_D2z_Px_C1001001_dd = I_ERI_Py_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Py_S_D2z_S_C1001001_dd;
  Double I_ERI_Pz_S_D2z_Px_C1001001_dd = I_ERI_Pz_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Pz_S_D2z_S_C1001001_dd;
  Double I_ERI_Px_S_Dxy_Py_C1001001_dd = I_ERI_Px_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Px_S_Dxy_S_C1001001_dd;
  Double I_ERI_Py_S_Dxy_Py_C1001001_dd = I_ERI_Py_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Py_S_Dxy_S_C1001001_dd;
  Double I_ERI_Pz_S_Dxy_Py_C1001001_dd = I_ERI_Pz_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Pz_S_Dxy_S_C1001001_dd;
  Double I_ERI_Px_S_Dxz_Py_C1001001_dd = I_ERI_Px_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Px_S_Dxz_S_C1001001_dd;
  Double I_ERI_Py_S_Dxz_Py_C1001001_dd = I_ERI_Py_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Py_S_Dxz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dxz_Py_C1001001_dd = I_ERI_Pz_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Pz_S_Dxz_S_C1001001_dd;
  Double I_ERI_Px_S_D2y_Py_C1001001_dd = I_ERI_Px_S_F3y_S_C1001001_dd+CDY*I_ERI_Px_S_D2y_S_C1001001_dd;
  Double I_ERI_Py_S_D2y_Py_C1001001_dd = I_ERI_Py_S_F3y_S_C1001001_dd+CDY*I_ERI_Py_S_D2y_S_C1001001_dd;
  Double I_ERI_Pz_S_D2y_Py_C1001001_dd = I_ERI_Pz_S_F3y_S_C1001001_dd+CDY*I_ERI_Pz_S_D2y_S_C1001001_dd;
  Double I_ERI_Px_S_Dyz_Py_C1001001_dd = I_ERI_Px_S_F2yz_S_C1001001_dd+CDY*I_ERI_Px_S_Dyz_S_C1001001_dd;
  Double I_ERI_Py_S_Dyz_Py_C1001001_dd = I_ERI_Py_S_F2yz_S_C1001001_dd+CDY*I_ERI_Py_S_Dyz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dyz_Py_C1001001_dd = I_ERI_Pz_S_F2yz_S_C1001001_dd+CDY*I_ERI_Pz_S_Dyz_S_C1001001_dd;
  Double I_ERI_Px_S_D2z_Py_C1001001_dd = I_ERI_Px_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Px_S_D2z_S_C1001001_dd;
  Double I_ERI_Py_S_D2z_Py_C1001001_dd = I_ERI_Py_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Py_S_D2z_S_C1001001_dd;
  Double I_ERI_Pz_S_D2z_Py_C1001001_dd = I_ERI_Pz_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Pz_S_D2z_S_C1001001_dd;
  Double I_ERI_Px_S_Dxz_Pz_C1001001_dd = I_ERI_Px_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Px_S_Dxz_S_C1001001_dd;
  Double I_ERI_Py_S_Dxz_Pz_C1001001_dd = I_ERI_Py_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Py_S_Dxz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dxz_Pz_C1001001_dd = I_ERI_Pz_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Pz_S_Dxz_S_C1001001_dd;
  Double I_ERI_Px_S_Dyz_Pz_C1001001_dd = I_ERI_Px_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Px_S_Dyz_S_C1001001_dd;
  Double I_ERI_Py_S_Dyz_Pz_C1001001_dd = I_ERI_Py_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Py_S_Dyz_S_C1001001_dd;
  Double I_ERI_Pz_S_Dyz_Pz_C1001001_dd = I_ERI_Pz_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Pz_S_Dyz_S_C1001001_dd;
  Double I_ERI_Px_S_D2z_Pz_C1001001_dd = I_ERI_Px_S_F3z_S_C1001001_dd+CDZ*I_ERI_Px_S_D2z_S_C1001001_dd;
  Double I_ERI_Py_S_D2z_Pz_C1001001_dd = I_ERI_Py_S_F3z_S_C1001001_dd+CDZ*I_ERI_Py_S_D2z_S_C1001001_dd;
  Double I_ERI_Pz_S_D2z_Pz_C1001001_dd = I_ERI_Pz_S_F3z_S_C1001001_dd+CDZ*I_ERI_Pz_S_D2z_S_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_D_P_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 24 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_C1001001_dd
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_dd
   ************************************************************/
  Double I_ERI_D2x_S_D2x_Px_C1001001_dd = I_ERI_D2x_S_F3x_S_C1001001_dd+CDX*I_ERI_D2x_S_D2x_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2x_Px_C1001001_dd = I_ERI_Dxy_S_F3x_S_C1001001_dd+CDX*I_ERI_Dxy_S_D2x_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2x_Px_C1001001_dd = I_ERI_Dxz_S_F3x_S_C1001001_dd+CDX*I_ERI_Dxz_S_D2x_S_C1001001_dd;
  Double I_ERI_D2y_S_D2x_Px_C1001001_dd = I_ERI_D2y_S_F3x_S_C1001001_dd+CDX*I_ERI_D2y_S_D2x_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2x_Px_C1001001_dd = I_ERI_Dyz_S_F3x_S_C1001001_dd+CDX*I_ERI_Dyz_S_D2x_S_C1001001_dd;
  Double I_ERI_D2z_S_D2x_Px_C1001001_dd = I_ERI_D2z_S_F3x_S_C1001001_dd+CDX*I_ERI_D2z_S_D2x_S_C1001001_dd;
  Double I_ERI_D2x_S_Dxy_Px_C1001001_dd = I_ERI_D2x_S_F2xy_S_C1001001_dd+CDX*I_ERI_D2x_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dxy_Px_C1001001_dd = I_ERI_Dxy_S_F2xy_S_C1001001_dd+CDX*I_ERI_Dxy_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dxy_Px_C1001001_dd = I_ERI_Dxz_S_F2xy_S_C1001001_dd+CDX*I_ERI_Dxz_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2y_S_Dxy_Px_C1001001_dd = I_ERI_D2y_S_F2xy_S_C1001001_dd+CDX*I_ERI_D2y_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dxy_Px_C1001001_dd = I_ERI_Dyz_S_F2xy_S_C1001001_dd+CDX*I_ERI_Dyz_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2z_S_Dxy_Px_C1001001_dd = I_ERI_D2z_S_F2xy_S_C1001001_dd+CDX*I_ERI_D2z_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2x_S_Dxz_Px_C1001001_dd = I_ERI_D2x_S_F2xz_S_C1001001_dd+CDX*I_ERI_D2x_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dxz_Px_C1001001_dd = I_ERI_Dxy_S_F2xz_S_C1001001_dd+CDX*I_ERI_Dxy_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dxz_Px_C1001001_dd = I_ERI_Dxz_S_F2xz_S_C1001001_dd+CDX*I_ERI_Dxz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dxz_Px_C1001001_dd = I_ERI_D2y_S_F2xz_S_C1001001_dd+CDX*I_ERI_D2y_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dxz_Px_C1001001_dd = I_ERI_Dyz_S_F2xz_S_C1001001_dd+CDX*I_ERI_Dyz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dxz_Px_C1001001_dd = I_ERI_D2z_S_F2xz_S_C1001001_dd+CDX*I_ERI_D2z_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2x_S_D2y_Px_C1001001_dd = I_ERI_D2x_S_Fx2y_S_C1001001_dd+CDX*I_ERI_D2x_S_D2y_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2y_Px_C1001001_dd = I_ERI_Dxy_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Dxy_S_D2y_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2y_Px_C1001001_dd = I_ERI_Dxz_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Dxz_S_D2y_S_C1001001_dd;
  Double I_ERI_D2y_S_D2y_Px_C1001001_dd = I_ERI_D2y_S_Fx2y_S_C1001001_dd+CDX*I_ERI_D2y_S_D2y_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2y_Px_C1001001_dd = I_ERI_Dyz_S_Fx2y_S_C1001001_dd+CDX*I_ERI_Dyz_S_D2y_S_C1001001_dd;
  Double I_ERI_D2z_S_D2y_Px_C1001001_dd = I_ERI_D2z_S_Fx2y_S_C1001001_dd+CDX*I_ERI_D2z_S_D2y_S_C1001001_dd;
  Double I_ERI_D2x_S_Dyz_Px_C1001001_dd = I_ERI_D2x_S_Fxyz_S_C1001001_dd+CDX*I_ERI_D2x_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dyz_Px_C1001001_dd = I_ERI_Dxy_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Dxy_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dyz_Px_C1001001_dd = I_ERI_Dxz_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Dxz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dyz_Px_C1001001_dd = I_ERI_D2y_S_Fxyz_S_C1001001_dd+CDX*I_ERI_D2y_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dyz_Px_C1001001_dd = I_ERI_Dyz_S_Fxyz_S_C1001001_dd+CDX*I_ERI_Dyz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dyz_Px_C1001001_dd = I_ERI_D2z_S_Fxyz_S_C1001001_dd+CDX*I_ERI_D2z_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2x_S_D2z_Px_C1001001_dd = I_ERI_D2x_S_Fx2z_S_C1001001_dd+CDX*I_ERI_D2x_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2z_Px_C1001001_dd = I_ERI_Dxy_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Dxy_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2z_Px_C1001001_dd = I_ERI_Dxz_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Dxz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2y_S_D2z_Px_C1001001_dd = I_ERI_D2y_S_Fx2z_S_C1001001_dd+CDX*I_ERI_D2y_S_D2z_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2z_Px_C1001001_dd = I_ERI_Dyz_S_Fx2z_S_C1001001_dd+CDX*I_ERI_Dyz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2z_S_D2z_Px_C1001001_dd = I_ERI_D2z_S_Fx2z_S_C1001001_dd+CDX*I_ERI_D2z_S_D2z_S_C1001001_dd;
  Double I_ERI_D2x_S_Dxy_Py_C1001001_dd = I_ERI_D2x_S_Fx2y_S_C1001001_dd+CDY*I_ERI_D2x_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dxy_Py_C1001001_dd = I_ERI_Dxy_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Dxy_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dxy_Py_C1001001_dd = I_ERI_Dxz_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Dxz_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2y_S_Dxy_Py_C1001001_dd = I_ERI_D2y_S_Fx2y_S_C1001001_dd+CDY*I_ERI_D2y_S_Dxy_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dxy_Py_C1001001_dd = I_ERI_Dyz_S_Fx2y_S_C1001001_dd+CDY*I_ERI_Dyz_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2z_S_Dxy_Py_C1001001_dd = I_ERI_D2z_S_Fx2y_S_C1001001_dd+CDY*I_ERI_D2z_S_Dxy_S_C1001001_dd;
  Double I_ERI_D2x_S_Dxz_Py_C1001001_dd = I_ERI_D2x_S_Fxyz_S_C1001001_dd+CDY*I_ERI_D2x_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dxz_Py_C1001001_dd = I_ERI_Dxy_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Dxy_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dxz_Py_C1001001_dd = I_ERI_Dxz_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Dxz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dxz_Py_C1001001_dd = I_ERI_D2y_S_Fxyz_S_C1001001_dd+CDY*I_ERI_D2y_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dxz_Py_C1001001_dd = I_ERI_Dyz_S_Fxyz_S_C1001001_dd+CDY*I_ERI_Dyz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dxz_Py_C1001001_dd = I_ERI_D2z_S_Fxyz_S_C1001001_dd+CDY*I_ERI_D2z_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2x_S_D2y_Py_C1001001_dd = I_ERI_D2x_S_F3y_S_C1001001_dd+CDY*I_ERI_D2x_S_D2y_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2y_Py_C1001001_dd = I_ERI_Dxy_S_F3y_S_C1001001_dd+CDY*I_ERI_Dxy_S_D2y_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2y_Py_C1001001_dd = I_ERI_Dxz_S_F3y_S_C1001001_dd+CDY*I_ERI_Dxz_S_D2y_S_C1001001_dd;
  Double I_ERI_D2y_S_D2y_Py_C1001001_dd = I_ERI_D2y_S_F3y_S_C1001001_dd+CDY*I_ERI_D2y_S_D2y_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2y_Py_C1001001_dd = I_ERI_Dyz_S_F3y_S_C1001001_dd+CDY*I_ERI_Dyz_S_D2y_S_C1001001_dd;
  Double I_ERI_D2z_S_D2y_Py_C1001001_dd = I_ERI_D2z_S_F3y_S_C1001001_dd+CDY*I_ERI_D2z_S_D2y_S_C1001001_dd;
  Double I_ERI_D2x_S_Dyz_Py_C1001001_dd = I_ERI_D2x_S_F2yz_S_C1001001_dd+CDY*I_ERI_D2x_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dyz_Py_C1001001_dd = I_ERI_Dxy_S_F2yz_S_C1001001_dd+CDY*I_ERI_Dxy_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dyz_Py_C1001001_dd = I_ERI_Dxz_S_F2yz_S_C1001001_dd+CDY*I_ERI_Dxz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dyz_Py_C1001001_dd = I_ERI_D2y_S_F2yz_S_C1001001_dd+CDY*I_ERI_D2y_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dyz_Py_C1001001_dd = I_ERI_Dyz_S_F2yz_S_C1001001_dd+CDY*I_ERI_Dyz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dyz_Py_C1001001_dd = I_ERI_D2z_S_F2yz_S_C1001001_dd+CDY*I_ERI_D2z_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2x_S_D2z_Py_C1001001_dd = I_ERI_D2x_S_Fy2z_S_C1001001_dd+CDY*I_ERI_D2x_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2z_Py_C1001001_dd = I_ERI_Dxy_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Dxy_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2z_Py_C1001001_dd = I_ERI_Dxz_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Dxz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2y_S_D2z_Py_C1001001_dd = I_ERI_D2y_S_Fy2z_S_C1001001_dd+CDY*I_ERI_D2y_S_D2z_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2z_Py_C1001001_dd = I_ERI_Dyz_S_Fy2z_S_C1001001_dd+CDY*I_ERI_Dyz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2z_S_D2z_Py_C1001001_dd = I_ERI_D2z_S_Fy2z_S_C1001001_dd+CDY*I_ERI_D2z_S_D2z_S_C1001001_dd;
  Double I_ERI_D2x_S_Dxz_Pz_C1001001_dd = I_ERI_D2x_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_D2x_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dxz_Pz_C1001001_dd = I_ERI_Dxy_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Dxy_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dxz_Pz_C1001001_dd = I_ERI_Dxz_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Dxz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dxz_Pz_C1001001_dd = I_ERI_D2y_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_D2y_S_Dxz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dxz_Pz_C1001001_dd = I_ERI_Dyz_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_Dyz_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dxz_Pz_C1001001_dd = I_ERI_D2z_S_Fx2z_S_C1001001_dd+CDZ*I_ERI_D2z_S_Dxz_S_C1001001_dd;
  Double I_ERI_D2x_S_Dyz_Pz_C1001001_dd = I_ERI_D2x_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_D2x_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxy_S_Dyz_Pz_C1001001_dd = I_ERI_Dxy_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Dxy_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dxz_S_Dyz_Pz_C1001001_dd = I_ERI_Dxz_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Dxz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2y_S_Dyz_Pz_C1001001_dd = I_ERI_D2y_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_D2y_S_Dyz_S_C1001001_dd;
  Double I_ERI_Dyz_S_Dyz_Pz_C1001001_dd = I_ERI_Dyz_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_Dyz_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2z_S_Dyz_Pz_C1001001_dd = I_ERI_D2z_S_Fy2z_S_C1001001_dd+CDZ*I_ERI_D2z_S_Dyz_S_C1001001_dd;
  Double I_ERI_D2x_S_D2z_Pz_C1001001_dd = I_ERI_D2x_S_F3z_S_C1001001_dd+CDZ*I_ERI_D2x_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxy_S_D2z_Pz_C1001001_dd = I_ERI_Dxy_S_F3z_S_C1001001_dd+CDZ*I_ERI_Dxy_S_D2z_S_C1001001_dd;
  Double I_ERI_Dxz_S_D2z_Pz_C1001001_dd = I_ERI_Dxz_S_F3z_S_C1001001_dd+CDZ*I_ERI_Dxz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2y_S_D2z_Pz_C1001001_dd = I_ERI_D2y_S_F3z_S_C1001001_dd+CDZ*I_ERI_D2y_S_D2z_S_C1001001_dd;
  Double I_ERI_Dyz_S_D2z_Pz_C1001001_dd = I_ERI_Dyz_S_F3z_S_C1001001_dd+CDZ*I_ERI_Dyz_S_D2z_S_C1001001_dd;
  Double I_ERI_D2z_S_D2z_Pz_C1001001_dd = I_ERI_D2z_S_F3z_S_C1001001_dd+CDZ*I_ERI_D2z_S_D2z_S_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_D_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001_dd
   ************************************************************/
  Double I_ERI_Px_S_Px_D2x_C1001001_dd = I_ERI_Px_S_D2x_Px_C1001001_dd+CDX*I_ERI_Px_S_Px_Px_C1001001_dd;
  Double I_ERI_Py_S_Px_D2x_C1001001_dd = I_ERI_Py_S_D2x_Px_C1001001_dd+CDX*I_ERI_Py_S_Px_Px_C1001001_dd;
  Double I_ERI_Pz_S_Px_D2x_C1001001_dd = I_ERI_Pz_S_D2x_Px_C1001001_dd+CDX*I_ERI_Pz_S_Px_Px_C1001001_dd;
  Double I_ERI_Px_S_Py_D2x_C1001001_dd = I_ERI_Px_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Px_S_Py_Px_C1001001_dd;
  Double I_ERI_Py_S_Py_D2x_C1001001_dd = I_ERI_Py_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Py_S_Py_Px_C1001001_dd;
  Double I_ERI_Pz_S_Py_D2x_C1001001_dd = I_ERI_Pz_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Pz_S_Py_Px_C1001001_dd;
  Double I_ERI_Px_S_Pz_D2x_C1001001_dd = I_ERI_Px_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Px_S_Pz_Px_C1001001_dd;
  Double I_ERI_Py_S_Pz_D2x_C1001001_dd = I_ERI_Py_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Py_S_Pz_Px_C1001001_dd;
  Double I_ERI_Pz_S_Pz_D2x_C1001001_dd = I_ERI_Pz_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Pz_S_Pz_Px_C1001001_dd;
  Double I_ERI_Px_S_Px_Dxy_C1001001_dd = I_ERI_Px_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Px_S_Px_Px_C1001001_dd;
  Double I_ERI_Py_S_Px_Dxy_C1001001_dd = I_ERI_Py_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Py_S_Px_Px_C1001001_dd;
  Double I_ERI_Pz_S_Px_Dxy_C1001001_dd = I_ERI_Pz_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Pz_S_Px_Px_C1001001_dd;
  Double I_ERI_Px_S_Py_Dxy_C1001001_dd = I_ERI_Px_S_D2y_Px_C1001001_dd+CDY*I_ERI_Px_S_Py_Px_C1001001_dd;
  Double I_ERI_Py_S_Py_Dxy_C1001001_dd = I_ERI_Py_S_D2y_Px_C1001001_dd+CDY*I_ERI_Py_S_Py_Px_C1001001_dd;
  Double I_ERI_Pz_S_Py_Dxy_C1001001_dd = I_ERI_Pz_S_D2y_Px_C1001001_dd+CDY*I_ERI_Pz_S_Py_Px_C1001001_dd;
  Double I_ERI_Px_S_Pz_Dxy_C1001001_dd = I_ERI_Px_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Px_S_Pz_Px_C1001001_dd;
  Double I_ERI_Py_S_Pz_Dxy_C1001001_dd = I_ERI_Py_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Py_S_Pz_Px_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Dxy_C1001001_dd = I_ERI_Pz_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Pz_S_Pz_Px_C1001001_dd;
  Double I_ERI_Px_S_Px_Dxz_C1001001_dd = I_ERI_Px_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Px_S_Px_Px_C1001001_dd;
  Double I_ERI_Py_S_Px_Dxz_C1001001_dd = I_ERI_Py_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Py_S_Px_Px_C1001001_dd;
  Double I_ERI_Pz_S_Px_Dxz_C1001001_dd = I_ERI_Pz_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Pz_S_Px_Px_C1001001_dd;
  Double I_ERI_Px_S_Py_Dxz_C1001001_dd = I_ERI_Px_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Px_S_Py_Px_C1001001_dd;
  Double I_ERI_Py_S_Py_Dxz_C1001001_dd = I_ERI_Py_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Py_S_Py_Px_C1001001_dd;
  Double I_ERI_Pz_S_Py_Dxz_C1001001_dd = I_ERI_Pz_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Pz_S_Py_Px_C1001001_dd;
  Double I_ERI_Px_S_Pz_Dxz_C1001001_dd = I_ERI_Px_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Px_S_Pz_Px_C1001001_dd;
  Double I_ERI_Py_S_Pz_Dxz_C1001001_dd = I_ERI_Py_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Py_S_Pz_Px_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Dxz_C1001001_dd = I_ERI_Pz_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Pz_S_Pz_Px_C1001001_dd;
  Double I_ERI_Px_S_Px_D2y_C1001001_dd = I_ERI_Px_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Px_S_Px_Py_C1001001_dd;
  Double I_ERI_Py_S_Px_D2y_C1001001_dd = I_ERI_Py_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Py_S_Px_Py_C1001001_dd;
  Double I_ERI_Pz_S_Px_D2y_C1001001_dd = I_ERI_Pz_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Pz_S_Px_Py_C1001001_dd;
  Double I_ERI_Px_S_Py_D2y_C1001001_dd = I_ERI_Px_S_D2y_Py_C1001001_dd+CDY*I_ERI_Px_S_Py_Py_C1001001_dd;
  Double I_ERI_Py_S_Py_D2y_C1001001_dd = I_ERI_Py_S_D2y_Py_C1001001_dd+CDY*I_ERI_Py_S_Py_Py_C1001001_dd;
  Double I_ERI_Pz_S_Py_D2y_C1001001_dd = I_ERI_Pz_S_D2y_Py_C1001001_dd+CDY*I_ERI_Pz_S_Py_Py_C1001001_dd;
  Double I_ERI_Px_S_Pz_D2y_C1001001_dd = I_ERI_Px_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Px_S_Pz_Py_C1001001_dd;
  Double I_ERI_Py_S_Pz_D2y_C1001001_dd = I_ERI_Py_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Py_S_Pz_Py_C1001001_dd;
  Double I_ERI_Pz_S_Pz_D2y_C1001001_dd = I_ERI_Pz_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Pz_S_Pz_Py_C1001001_dd;
  Double I_ERI_Px_S_Px_Dyz_C1001001_dd = I_ERI_Px_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Px_S_Px_Py_C1001001_dd;
  Double I_ERI_Py_S_Px_Dyz_C1001001_dd = I_ERI_Py_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Py_S_Px_Py_C1001001_dd;
  Double I_ERI_Pz_S_Px_Dyz_C1001001_dd = I_ERI_Pz_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Pz_S_Px_Py_C1001001_dd;
  Double I_ERI_Px_S_Py_Dyz_C1001001_dd = I_ERI_Px_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Px_S_Py_Py_C1001001_dd;
  Double I_ERI_Py_S_Py_Dyz_C1001001_dd = I_ERI_Py_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Py_S_Py_Py_C1001001_dd;
  Double I_ERI_Pz_S_Py_Dyz_C1001001_dd = I_ERI_Pz_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Pz_S_Py_Py_C1001001_dd;
  Double I_ERI_Px_S_Pz_Dyz_C1001001_dd = I_ERI_Px_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Px_S_Pz_Py_C1001001_dd;
  Double I_ERI_Py_S_Pz_Dyz_C1001001_dd = I_ERI_Py_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Py_S_Pz_Py_C1001001_dd;
  Double I_ERI_Pz_S_Pz_Dyz_C1001001_dd = I_ERI_Pz_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Pz_S_Pz_Py_C1001001_dd;
  Double I_ERI_Px_S_Px_D2z_C1001001_dd = I_ERI_Px_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Px_S_Px_Pz_C1001001_dd;
  Double I_ERI_Py_S_Px_D2z_C1001001_dd = I_ERI_Py_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Py_S_Px_Pz_C1001001_dd;
  Double I_ERI_Pz_S_Px_D2z_C1001001_dd = I_ERI_Pz_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Pz_S_Px_Pz_C1001001_dd;
  Double I_ERI_Px_S_Py_D2z_C1001001_dd = I_ERI_Px_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Px_S_Py_Pz_C1001001_dd;
  Double I_ERI_Py_S_Py_D2z_C1001001_dd = I_ERI_Py_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Py_S_Py_Pz_C1001001_dd;
  Double I_ERI_Pz_S_Py_D2z_C1001001_dd = I_ERI_Pz_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Pz_S_Py_Pz_C1001001_dd;
  Double I_ERI_Px_S_Pz_D2z_C1001001_dd = I_ERI_Px_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Px_S_Pz_Pz_C1001001_dd;
  Double I_ERI_Py_S_Pz_D2z_C1001001_dd = I_ERI_Py_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Py_S_Pz_Pz_C1001001_dd;
  Double I_ERI_Pz_S_Pz_D2z_C1001001_dd = I_ERI_Pz_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Pz_S_Pz_Pz_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_D_C1001001_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_P_C1001001_dd
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001_dd
   ************************************************************/
  Double I_ERI_D2x_S_Px_D2x_C1001001_dd = I_ERI_D2x_S_D2x_Px_C1001001_dd+CDX*I_ERI_D2x_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Px_D2x_C1001001_dd = I_ERI_Dxy_S_D2x_Px_C1001001_dd+CDX*I_ERI_Dxy_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Px_D2x_C1001001_dd = I_ERI_Dxz_S_D2x_Px_C1001001_dd+CDX*I_ERI_Dxz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2y_S_Px_D2x_C1001001_dd = I_ERI_D2y_S_D2x_Px_C1001001_dd+CDX*I_ERI_D2y_S_Px_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Px_D2x_C1001001_dd = I_ERI_Dyz_S_D2x_Px_C1001001_dd+CDX*I_ERI_Dyz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2z_S_Px_D2x_C1001001_dd = I_ERI_D2z_S_D2x_Px_C1001001_dd+CDX*I_ERI_D2z_S_Px_Px_C1001001_dd;
  Double I_ERI_D2x_S_Py_D2x_C1001001_dd = I_ERI_D2x_S_Dxy_Px_C1001001_dd+CDX*I_ERI_D2x_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Py_D2x_C1001001_dd = I_ERI_Dxy_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Dxy_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Py_D2x_C1001001_dd = I_ERI_Dxz_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Dxz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2y_S_Py_D2x_C1001001_dd = I_ERI_D2y_S_Dxy_Px_C1001001_dd+CDX*I_ERI_D2y_S_Py_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Py_D2x_C1001001_dd = I_ERI_Dyz_S_Dxy_Px_C1001001_dd+CDX*I_ERI_Dyz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2z_S_Py_D2x_C1001001_dd = I_ERI_D2z_S_Dxy_Px_C1001001_dd+CDX*I_ERI_D2z_S_Py_Px_C1001001_dd;
  Double I_ERI_D2x_S_Pz_D2x_C1001001_dd = I_ERI_D2x_S_Dxz_Px_C1001001_dd+CDX*I_ERI_D2x_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_D2x_C1001001_dd = I_ERI_Dxy_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Dxy_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_D2x_C1001001_dd = I_ERI_Dxz_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Dxz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2y_S_Pz_D2x_C1001001_dd = I_ERI_D2y_S_Dxz_Px_C1001001_dd+CDX*I_ERI_D2y_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_D2x_C1001001_dd = I_ERI_Dyz_S_Dxz_Px_C1001001_dd+CDX*I_ERI_Dyz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2z_S_Pz_D2x_C1001001_dd = I_ERI_D2z_S_Dxz_Px_C1001001_dd+CDX*I_ERI_D2z_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2x_S_Px_Dxy_C1001001_dd = I_ERI_D2x_S_Dxy_Px_C1001001_dd+CDY*I_ERI_D2x_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Dxy_C1001001_dd = I_ERI_Dxy_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Dxy_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Dxy_C1001001_dd = I_ERI_Dxz_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Dxz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2y_S_Px_Dxy_C1001001_dd = I_ERI_D2y_S_Dxy_Px_C1001001_dd+CDY*I_ERI_D2y_S_Px_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Dxy_C1001001_dd = I_ERI_Dyz_S_Dxy_Px_C1001001_dd+CDY*I_ERI_Dyz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2z_S_Px_Dxy_C1001001_dd = I_ERI_D2z_S_Dxy_Px_C1001001_dd+CDY*I_ERI_D2z_S_Px_Px_C1001001_dd;
  Double I_ERI_D2x_S_Py_Dxy_C1001001_dd = I_ERI_D2x_S_D2y_Px_C1001001_dd+CDY*I_ERI_D2x_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Dxy_C1001001_dd = I_ERI_Dxy_S_D2y_Px_C1001001_dd+CDY*I_ERI_Dxy_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Dxy_C1001001_dd = I_ERI_Dxz_S_D2y_Px_C1001001_dd+CDY*I_ERI_Dxz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2y_S_Py_Dxy_C1001001_dd = I_ERI_D2y_S_D2y_Px_C1001001_dd+CDY*I_ERI_D2y_S_Py_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Dxy_C1001001_dd = I_ERI_Dyz_S_D2y_Px_C1001001_dd+CDY*I_ERI_Dyz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2z_S_Py_Dxy_C1001001_dd = I_ERI_D2z_S_D2y_Px_C1001001_dd+CDY*I_ERI_D2z_S_Py_Px_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Dxy_C1001001_dd = I_ERI_D2x_S_Dyz_Px_C1001001_dd+CDY*I_ERI_D2x_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Dxy_C1001001_dd = I_ERI_Dxy_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Dxy_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Dxy_C1001001_dd = I_ERI_Dxz_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Dxz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Dxy_C1001001_dd = I_ERI_D2y_S_Dyz_Px_C1001001_dd+CDY*I_ERI_D2y_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Dxy_C1001001_dd = I_ERI_Dyz_S_Dyz_Px_C1001001_dd+CDY*I_ERI_Dyz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Dxy_C1001001_dd = I_ERI_D2z_S_Dyz_Px_C1001001_dd+CDY*I_ERI_D2z_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2x_S_Px_Dxz_C1001001_dd = I_ERI_D2x_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_D2x_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Dxz_C1001001_dd = I_ERI_Dxy_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Dxy_S_Px_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Dxz_C1001001_dd = I_ERI_Dxz_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Dxz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2y_S_Px_Dxz_C1001001_dd = I_ERI_D2y_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_D2y_S_Px_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Dxz_C1001001_dd = I_ERI_Dyz_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_Dyz_S_Px_Px_C1001001_dd;
  Double I_ERI_D2z_S_Px_Dxz_C1001001_dd = I_ERI_D2z_S_Dxz_Px_C1001001_dd+CDZ*I_ERI_D2z_S_Px_Px_C1001001_dd;
  Double I_ERI_D2x_S_Py_Dxz_C1001001_dd = I_ERI_D2x_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_D2x_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Dxz_C1001001_dd = I_ERI_Dxy_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Dxy_S_Py_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Dxz_C1001001_dd = I_ERI_Dxz_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Dxz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2y_S_Py_Dxz_C1001001_dd = I_ERI_D2y_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_D2y_S_Py_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Dxz_C1001001_dd = I_ERI_Dyz_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_Dyz_S_Py_Px_C1001001_dd;
  Double I_ERI_D2z_S_Py_Dxz_C1001001_dd = I_ERI_D2z_S_Dyz_Px_C1001001_dd+CDZ*I_ERI_D2z_S_Py_Px_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Dxz_C1001001_dd = I_ERI_D2x_S_D2z_Px_C1001001_dd+CDZ*I_ERI_D2x_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Dxz_C1001001_dd = I_ERI_Dxy_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Dxy_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Dxz_C1001001_dd = I_ERI_Dxz_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Dxz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Dxz_C1001001_dd = I_ERI_D2y_S_D2z_Px_C1001001_dd+CDZ*I_ERI_D2y_S_Pz_Px_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Dxz_C1001001_dd = I_ERI_Dyz_S_D2z_Px_C1001001_dd+CDZ*I_ERI_Dyz_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Dxz_C1001001_dd = I_ERI_D2z_S_D2z_Px_C1001001_dd+CDZ*I_ERI_D2z_S_Pz_Px_C1001001_dd;
  Double I_ERI_D2x_S_Px_D2y_C1001001_dd = I_ERI_D2x_S_Dxy_Py_C1001001_dd+CDY*I_ERI_D2x_S_Px_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Px_D2y_C1001001_dd = I_ERI_Dxy_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Dxy_S_Px_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Px_D2y_C1001001_dd = I_ERI_Dxz_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Dxz_S_Px_Py_C1001001_dd;
  Double I_ERI_D2y_S_Px_D2y_C1001001_dd = I_ERI_D2y_S_Dxy_Py_C1001001_dd+CDY*I_ERI_D2y_S_Px_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Px_D2y_C1001001_dd = I_ERI_Dyz_S_Dxy_Py_C1001001_dd+CDY*I_ERI_Dyz_S_Px_Py_C1001001_dd;
  Double I_ERI_D2z_S_Px_D2y_C1001001_dd = I_ERI_D2z_S_Dxy_Py_C1001001_dd+CDY*I_ERI_D2z_S_Px_Py_C1001001_dd;
  Double I_ERI_D2x_S_Py_D2y_C1001001_dd = I_ERI_D2x_S_D2y_Py_C1001001_dd+CDY*I_ERI_D2x_S_Py_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Py_D2y_C1001001_dd = I_ERI_Dxy_S_D2y_Py_C1001001_dd+CDY*I_ERI_Dxy_S_Py_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Py_D2y_C1001001_dd = I_ERI_Dxz_S_D2y_Py_C1001001_dd+CDY*I_ERI_Dxz_S_Py_Py_C1001001_dd;
  Double I_ERI_D2y_S_Py_D2y_C1001001_dd = I_ERI_D2y_S_D2y_Py_C1001001_dd+CDY*I_ERI_D2y_S_Py_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Py_D2y_C1001001_dd = I_ERI_Dyz_S_D2y_Py_C1001001_dd+CDY*I_ERI_Dyz_S_Py_Py_C1001001_dd;
  Double I_ERI_D2z_S_Py_D2y_C1001001_dd = I_ERI_D2z_S_D2y_Py_C1001001_dd+CDY*I_ERI_D2z_S_Py_Py_C1001001_dd;
  Double I_ERI_D2x_S_Pz_D2y_C1001001_dd = I_ERI_D2x_S_Dyz_Py_C1001001_dd+CDY*I_ERI_D2x_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_D2y_C1001001_dd = I_ERI_Dxy_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Dxy_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_D2y_C1001001_dd = I_ERI_Dxz_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Dxz_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2y_S_Pz_D2y_C1001001_dd = I_ERI_D2y_S_Dyz_Py_C1001001_dd+CDY*I_ERI_D2y_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_D2y_C1001001_dd = I_ERI_Dyz_S_Dyz_Py_C1001001_dd+CDY*I_ERI_Dyz_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2z_S_Pz_D2y_C1001001_dd = I_ERI_D2z_S_Dyz_Py_C1001001_dd+CDY*I_ERI_D2z_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2x_S_Px_Dyz_C1001001_dd = I_ERI_D2x_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_D2x_S_Px_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Px_Dyz_C1001001_dd = I_ERI_Dxy_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Dxy_S_Px_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Px_Dyz_C1001001_dd = I_ERI_Dxz_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Dxz_S_Px_Py_C1001001_dd;
  Double I_ERI_D2y_S_Px_Dyz_C1001001_dd = I_ERI_D2y_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_D2y_S_Px_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Px_Dyz_C1001001_dd = I_ERI_Dyz_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_Dyz_S_Px_Py_C1001001_dd;
  Double I_ERI_D2z_S_Px_Dyz_C1001001_dd = I_ERI_D2z_S_Dxz_Py_C1001001_dd+CDZ*I_ERI_D2z_S_Px_Py_C1001001_dd;
  Double I_ERI_D2x_S_Py_Dyz_C1001001_dd = I_ERI_D2x_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_D2x_S_Py_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Py_Dyz_C1001001_dd = I_ERI_Dxy_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Dxy_S_Py_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Py_Dyz_C1001001_dd = I_ERI_Dxz_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Dxz_S_Py_Py_C1001001_dd;
  Double I_ERI_D2y_S_Py_Dyz_C1001001_dd = I_ERI_D2y_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_D2y_S_Py_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Py_Dyz_C1001001_dd = I_ERI_Dyz_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_Dyz_S_Py_Py_C1001001_dd;
  Double I_ERI_D2z_S_Py_Dyz_C1001001_dd = I_ERI_D2z_S_Dyz_Py_C1001001_dd+CDZ*I_ERI_D2z_S_Py_Py_C1001001_dd;
  Double I_ERI_D2x_S_Pz_Dyz_C1001001_dd = I_ERI_D2x_S_D2z_Py_C1001001_dd+CDZ*I_ERI_D2x_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_Dyz_C1001001_dd = I_ERI_Dxy_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Dxy_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_Dyz_C1001001_dd = I_ERI_Dxz_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Dxz_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2y_S_Pz_Dyz_C1001001_dd = I_ERI_D2y_S_D2z_Py_C1001001_dd+CDZ*I_ERI_D2y_S_Pz_Py_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_Dyz_C1001001_dd = I_ERI_Dyz_S_D2z_Py_C1001001_dd+CDZ*I_ERI_Dyz_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2z_S_Pz_Dyz_C1001001_dd = I_ERI_D2z_S_D2z_Py_C1001001_dd+CDZ*I_ERI_D2z_S_Pz_Py_C1001001_dd;
  Double I_ERI_D2x_S_Px_D2z_C1001001_dd = I_ERI_D2x_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_D2x_S_Px_Pz_C1001001_dd;
  Double I_ERI_Dxy_S_Px_D2z_C1001001_dd = I_ERI_Dxy_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Dxy_S_Px_Pz_C1001001_dd;
  Double I_ERI_Dxz_S_Px_D2z_C1001001_dd = I_ERI_Dxz_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Dxz_S_Px_Pz_C1001001_dd;
  Double I_ERI_D2y_S_Px_D2z_C1001001_dd = I_ERI_D2y_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_D2y_S_Px_Pz_C1001001_dd;
  Double I_ERI_Dyz_S_Px_D2z_C1001001_dd = I_ERI_Dyz_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_Dyz_S_Px_Pz_C1001001_dd;
  Double I_ERI_D2z_S_Px_D2z_C1001001_dd = I_ERI_D2z_S_Dxz_Pz_C1001001_dd+CDZ*I_ERI_D2z_S_Px_Pz_C1001001_dd;
  Double I_ERI_D2x_S_Py_D2z_C1001001_dd = I_ERI_D2x_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_D2x_S_Py_Pz_C1001001_dd;
  Double I_ERI_Dxy_S_Py_D2z_C1001001_dd = I_ERI_Dxy_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Dxy_S_Py_Pz_C1001001_dd;
  Double I_ERI_Dxz_S_Py_D2z_C1001001_dd = I_ERI_Dxz_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Dxz_S_Py_Pz_C1001001_dd;
  Double I_ERI_D2y_S_Py_D2z_C1001001_dd = I_ERI_D2y_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_D2y_S_Py_Pz_C1001001_dd;
  Double I_ERI_Dyz_S_Py_D2z_C1001001_dd = I_ERI_Dyz_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_Dyz_S_Py_Pz_C1001001_dd;
  Double I_ERI_D2z_S_Py_D2z_C1001001_dd = I_ERI_D2z_S_Dyz_Pz_C1001001_dd+CDZ*I_ERI_D2z_S_Py_Pz_C1001001_dd;
  Double I_ERI_D2x_S_Pz_D2z_C1001001_dd = I_ERI_D2x_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_D2x_S_Pz_Pz_C1001001_dd;
  Double I_ERI_Dxy_S_Pz_D2z_C1001001_dd = I_ERI_Dxy_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Dxy_S_Pz_Pz_C1001001_dd;
  Double I_ERI_Dxz_S_Pz_D2z_C1001001_dd = I_ERI_Dxz_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Dxz_S_Pz_Pz_C1001001_dd;
  Double I_ERI_D2y_S_Pz_D2z_C1001001_dd = I_ERI_D2y_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_D2y_S_Pz_Pz_C1001001_dd;
  Double I_ERI_Dyz_S_Pz_D2z_C1001001_dd = I_ERI_Dyz_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_Dyz_S_Pz_Pz_C1001001_dd;
  Double I_ERI_D2z_S_Pz_D2z_C1001001_dd = I_ERI_D2z_S_D2z_Pz_C1001001_dd+CDZ*I_ERI_D2z_S_Pz_Pz_C1001001_dd;

  /************************************************************
   * declare the HRR2 result shell quartets in array form
   ************************************************************/

  /************************************************************
   * initilize the HRR steps : build the AB/CD variables
   ************************************************************/
  Double ABX = A[0] - B[0];
  Double ABY = A[1] - B[1];
  Double ABZ = A[2] - B[2];

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_a
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001001_a = I_ERI_D2x_S_Px_S_C1001001_a+ABX*I_ERI_Px_S_Px_S_C1001001_a;
  Double I_ERI_Py_Px_Px_S_C1001001_a = I_ERI_Dxy_S_Px_S_C1001001_a+ABX*I_ERI_Py_S_Px_S_C1001001_a;
  Double I_ERI_Pz_Px_Px_S_C1001001_a = I_ERI_Dxz_S_Px_S_C1001001_a+ABX*I_ERI_Pz_S_Px_S_C1001001_a;
  Double I_ERI_Px_Py_Px_S_C1001001_a = I_ERI_Dxy_S_Px_S_C1001001_a+ABY*I_ERI_Px_S_Px_S_C1001001_a;
  Double I_ERI_Py_Py_Px_S_C1001001_a = I_ERI_D2y_S_Px_S_C1001001_a+ABY*I_ERI_Py_S_Px_S_C1001001_a;
  Double I_ERI_Pz_Py_Px_S_C1001001_a = I_ERI_Dyz_S_Px_S_C1001001_a+ABY*I_ERI_Pz_S_Px_S_C1001001_a;
  Double I_ERI_Px_Pz_Px_S_C1001001_a = I_ERI_Dxz_S_Px_S_C1001001_a+ABZ*I_ERI_Px_S_Px_S_C1001001_a;
  Double I_ERI_Py_Pz_Px_S_C1001001_a = I_ERI_Dyz_S_Px_S_C1001001_a+ABZ*I_ERI_Py_S_Px_S_C1001001_a;
  Double I_ERI_Pz_Pz_Px_S_C1001001_a = I_ERI_D2z_S_Px_S_C1001001_a+ABZ*I_ERI_Pz_S_Px_S_C1001001_a;
  Double I_ERI_Px_Px_Py_S_C1001001_a = I_ERI_D2x_S_Py_S_C1001001_a+ABX*I_ERI_Px_S_Py_S_C1001001_a;
  Double I_ERI_Py_Px_Py_S_C1001001_a = I_ERI_Dxy_S_Py_S_C1001001_a+ABX*I_ERI_Py_S_Py_S_C1001001_a;
  Double I_ERI_Pz_Px_Py_S_C1001001_a = I_ERI_Dxz_S_Py_S_C1001001_a+ABX*I_ERI_Pz_S_Py_S_C1001001_a;
  Double I_ERI_Px_Py_Py_S_C1001001_a = I_ERI_Dxy_S_Py_S_C1001001_a+ABY*I_ERI_Px_S_Py_S_C1001001_a;
  Double I_ERI_Py_Py_Py_S_C1001001_a = I_ERI_D2y_S_Py_S_C1001001_a+ABY*I_ERI_Py_S_Py_S_C1001001_a;
  Double I_ERI_Pz_Py_Py_S_C1001001_a = I_ERI_Dyz_S_Py_S_C1001001_a+ABY*I_ERI_Pz_S_Py_S_C1001001_a;
  Double I_ERI_Px_Pz_Py_S_C1001001_a = I_ERI_Dxz_S_Py_S_C1001001_a+ABZ*I_ERI_Px_S_Py_S_C1001001_a;
  Double I_ERI_Py_Pz_Py_S_C1001001_a = I_ERI_Dyz_S_Py_S_C1001001_a+ABZ*I_ERI_Py_S_Py_S_C1001001_a;
  Double I_ERI_Pz_Pz_Py_S_C1001001_a = I_ERI_D2z_S_Py_S_C1001001_a+ABZ*I_ERI_Pz_S_Py_S_C1001001_a;
  Double I_ERI_Px_Px_Pz_S_C1001001_a = I_ERI_D2x_S_Pz_S_C1001001_a+ABX*I_ERI_Px_S_Pz_S_C1001001_a;
  Double I_ERI_Py_Px_Pz_S_C1001001_a = I_ERI_Dxy_S_Pz_S_C1001001_a+ABX*I_ERI_Py_S_Pz_S_C1001001_a;
  Double I_ERI_Pz_Px_Pz_S_C1001001_a = I_ERI_Dxz_S_Pz_S_C1001001_a+ABX*I_ERI_Pz_S_Pz_S_C1001001_a;
  Double I_ERI_Px_Py_Pz_S_C1001001_a = I_ERI_Dxy_S_Pz_S_C1001001_a+ABY*I_ERI_Px_S_Pz_S_C1001001_a;
  Double I_ERI_Py_Py_Pz_S_C1001001_a = I_ERI_D2y_S_Pz_S_C1001001_a+ABY*I_ERI_Py_S_Pz_S_C1001001_a;
  Double I_ERI_Pz_Py_Pz_S_C1001001_a = I_ERI_Dyz_S_Pz_S_C1001001_a+ABY*I_ERI_Pz_S_Pz_S_C1001001_a;
  Double I_ERI_Px_Pz_Pz_S_C1001001_a = I_ERI_Dxz_S_Pz_S_C1001001_a+ABZ*I_ERI_Px_S_Pz_S_C1001001_a;
  Double I_ERI_Py_Pz_Pz_S_C1001001_a = I_ERI_Dyz_S_Pz_S_C1001001_a+ABZ*I_ERI_Py_S_Pz_S_C1001001_a;
  Double I_ERI_Pz_Pz_Pz_S_C1001001_a = I_ERI_D2z_S_Pz_S_C1001001_a+ABZ*I_ERI_Pz_S_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001001_a
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001001_a = I_ERI_F3x_S_S_S_C1001001_a+ABX*I_ERI_D2x_S_S_S_C1001001_a;
  Double I_ERI_Dxy_Px_S_S_C1001001_a = I_ERI_F2xy_S_S_S_C1001001_a+ABX*I_ERI_Dxy_S_S_S_C1001001_a;
  Double I_ERI_Dxz_Px_S_S_C1001001_a = I_ERI_F2xz_S_S_S_C1001001_a+ABX*I_ERI_Dxz_S_S_S_C1001001_a;
  Double I_ERI_D2y_Px_S_S_C1001001_a = I_ERI_Fx2y_S_S_S_C1001001_a+ABX*I_ERI_D2y_S_S_S_C1001001_a;
  Double I_ERI_Dyz_Px_S_S_C1001001_a = I_ERI_Fxyz_S_S_S_C1001001_a+ABX*I_ERI_Dyz_S_S_S_C1001001_a;
  Double I_ERI_D2z_Px_S_S_C1001001_a = I_ERI_Fx2z_S_S_S_C1001001_a+ABX*I_ERI_D2z_S_S_S_C1001001_a;
  Double I_ERI_D2x_Py_S_S_C1001001_a = I_ERI_F2xy_S_S_S_C1001001_a+ABY*I_ERI_D2x_S_S_S_C1001001_a;
  Double I_ERI_Dxy_Py_S_S_C1001001_a = I_ERI_Fx2y_S_S_S_C1001001_a+ABY*I_ERI_Dxy_S_S_S_C1001001_a;
  Double I_ERI_Dxz_Py_S_S_C1001001_a = I_ERI_Fxyz_S_S_S_C1001001_a+ABY*I_ERI_Dxz_S_S_S_C1001001_a;
  Double I_ERI_D2y_Py_S_S_C1001001_a = I_ERI_F3y_S_S_S_C1001001_a+ABY*I_ERI_D2y_S_S_S_C1001001_a;
  Double I_ERI_Dyz_Py_S_S_C1001001_a = I_ERI_F2yz_S_S_S_C1001001_a+ABY*I_ERI_Dyz_S_S_S_C1001001_a;
  Double I_ERI_D2z_Py_S_S_C1001001_a = I_ERI_Fy2z_S_S_S_C1001001_a+ABY*I_ERI_D2z_S_S_S_C1001001_a;
  Double I_ERI_D2x_Pz_S_S_C1001001_a = I_ERI_F2xz_S_S_S_C1001001_a+ABZ*I_ERI_D2x_S_S_S_C1001001_a;
  Double I_ERI_Dxy_Pz_S_S_C1001001_a = I_ERI_Fxyz_S_S_S_C1001001_a+ABZ*I_ERI_Dxy_S_S_S_C1001001_a;
  Double I_ERI_Dxz_Pz_S_S_C1001001_a = I_ERI_Fx2z_S_S_S_C1001001_a+ABZ*I_ERI_Dxz_S_S_S_C1001001_a;
  Double I_ERI_D2y_Pz_S_S_C1001001_a = I_ERI_F2yz_S_S_S_C1001001_a+ABZ*I_ERI_D2y_S_S_S_C1001001_a;
  Double I_ERI_Dyz_Pz_S_S_C1001001_a = I_ERI_Fy2z_S_S_S_C1001001_a+ABZ*I_ERI_Dyz_S_S_S_C1001001_a;
  Double I_ERI_D2z_Pz_S_S_C1001001_a = I_ERI_F3z_S_S_S_C1001001_a+ABZ*I_ERI_D2z_S_S_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_c
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001001_c = I_ERI_D2x_S_Px_S_C1001001_c+ABX*I_ERI_Px_S_Px_S_C1001001_c;
  Double I_ERI_Py_Px_Px_S_C1001001_c = I_ERI_Dxy_S_Px_S_C1001001_c+ABX*I_ERI_Py_S_Px_S_C1001001_c;
  Double I_ERI_Pz_Px_Px_S_C1001001_c = I_ERI_Dxz_S_Px_S_C1001001_c+ABX*I_ERI_Pz_S_Px_S_C1001001_c;
  Double I_ERI_Px_Py_Px_S_C1001001_c = I_ERI_Dxy_S_Px_S_C1001001_c+ABY*I_ERI_Px_S_Px_S_C1001001_c;
  Double I_ERI_Py_Py_Px_S_C1001001_c = I_ERI_D2y_S_Px_S_C1001001_c+ABY*I_ERI_Py_S_Px_S_C1001001_c;
  Double I_ERI_Pz_Py_Px_S_C1001001_c = I_ERI_Dyz_S_Px_S_C1001001_c+ABY*I_ERI_Pz_S_Px_S_C1001001_c;
  Double I_ERI_Px_Pz_Px_S_C1001001_c = I_ERI_Dxz_S_Px_S_C1001001_c+ABZ*I_ERI_Px_S_Px_S_C1001001_c;
  Double I_ERI_Py_Pz_Px_S_C1001001_c = I_ERI_Dyz_S_Px_S_C1001001_c+ABZ*I_ERI_Py_S_Px_S_C1001001_c;
  Double I_ERI_Pz_Pz_Px_S_C1001001_c = I_ERI_D2z_S_Px_S_C1001001_c+ABZ*I_ERI_Pz_S_Px_S_C1001001_c;
  Double I_ERI_Px_Px_Py_S_C1001001_c = I_ERI_D2x_S_Py_S_C1001001_c+ABX*I_ERI_Px_S_Py_S_C1001001_c;
  Double I_ERI_Py_Px_Py_S_C1001001_c = I_ERI_Dxy_S_Py_S_C1001001_c+ABX*I_ERI_Py_S_Py_S_C1001001_c;
  Double I_ERI_Pz_Px_Py_S_C1001001_c = I_ERI_Dxz_S_Py_S_C1001001_c+ABX*I_ERI_Pz_S_Py_S_C1001001_c;
  Double I_ERI_Px_Py_Py_S_C1001001_c = I_ERI_Dxy_S_Py_S_C1001001_c+ABY*I_ERI_Px_S_Py_S_C1001001_c;
  Double I_ERI_Py_Py_Py_S_C1001001_c = I_ERI_D2y_S_Py_S_C1001001_c+ABY*I_ERI_Py_S_Py_S_C1001001_c;
  Double I_ERI_Pz_Py_Py_S_C1001001_c = I_ERI_Dyz_S_Py_S_C1001001_c+ABY*I_ERI_Pz_S_Py_S_C1001001_c;
  Double I_ERI_Px_Pz_Py_S_C1001001_c = I_ERI_Dxz_S_Py_S_C1001001_c+ABZ*I_ERI_Px_S_Py_S_C1001001_c;
  Double I_ERI_Py_Pz_Py_S_C1001001_c = I_ERI_Dyz_S_Py_S_C1001001_c+ABZ*I_ERI_Py_S_Py_S_C1001001_c;
  Double I_ERI_Pz_Pz_Py_S_C1001001_c = I_ERI_D2z_S_Py_S_C1001001_c+ABZ*I_ERI_Pz_S_Py_S_C1001001_c;
  Double I_ERI_Px_Px_Pz_S_C1001001_c = I_ERI_D2x_S_Pz_S_C1001001_c+ABX*I_ERI_Px_S_Pz_S_C1001001_c;
  Double I_ERI_Py_Px_Pz_S_C1001001_c = I_ERI_Dxy_S_Pz_S_C1001001_c+ABX*I_ERI_Py_S_Pz_S_C1001001_c;
  Double I_ERI_Pz_Px_Pz_S_C1001001_c = I_ERI_Dxz_S_Pz_S_C1001001_c+ABX*I_ERI_Pz_S_Pz_S_C1001001_c;
  Double I_ERI_Px_Py_Pz_S_C1001001_c = I_ERI_Dxy_S_Pz_S_C1001001_c+ABY*I_ERI_Px_S_Pz_S_C1001001_c;
  Double I_ERI_Py_Py_Pz_S_C1001001_c = I_ERI_D2y_S_Pz_S_C1001001_c+ABY*I_ERI_Py_S_Pz_S_C1001001_c;
  Double I_ERI_Pz_Py_Pz_S_C1001001_c = I_ERI_Dyz_S_Pz_S_C1001001_c+ABY*I_ERI_Pz_S_Pz_S_C1001001_c;
  Double I_ERI_Px_Pz_Pz_S_C1001001_c = I_ERI_Dxz_S_Pz_S_C1001001_c+ABZ*I_ERI_Px_S_Pz_S_C1001001_c;
  Double I_ERI_Py_Pz_Pz_S_C1001001_c = I_ERI_Dyz_S_Pz_S_C1001001_c+ABZ*I_ERI_Py_S_Pz_S_C1001001_c;
  Double I_ERI_Pz_Pz_Pz_S_C1001001_c = I_ERI_D2z_S_Pz_S_C1001001_c+ABZ*I_ERI_Pz_S_Pz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_d
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_d
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001001_d = I_ERI_D2x_S_Px_S_C1001001_d+ABX*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_Px_Px_S_C1001001_d = I_ERI_Dxy_S_Px_S_C1001001_d+ABX*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_Px_Px_S_C1001001_d = I_ERI_Dxz_S_Px_S_C1001001_d+ABX*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_Py_Px_S_C1001001_d = I_ERI_Dxy_S_Px_S_C1001001_d+ABY*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_Py_Px_S_C1001001_d = I_ERI_D2y_S_Px_S_C1001001_d+ABY*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_Py_Px_S_C1001001_d = I_ERI_Dyz_S_Px_S_C1001001_d+ABY*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_Pz_Px_S_C1001001_d = I_ERI_Dxz_S_Px_S_C1001001_d+ABZ*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_Pz_Px_S_C1001001_d = I_ERI_Dyz_S_Px_S_C1001001_d+ABZ*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_Pz_Px_S_C1001001_d = I_ERI_D2z_S_Px_S_C1001001_d+ABZ*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_Px_Py_S_C1001001_d = I_ERI_D2x_S_Py_S_C1001001_d+ABX*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_Px_Py_S_C1001001_d = I_ERI_Dxy_S_Py_S_C1001001_d+ABX*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_Px_Py_S_C1001001_d = I_ERI_Dxz_S_Py_S_C1001001_d+ABX*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_Py_Py_S_C1001001_d = I_ERI_Dxy_S_Py_S_C1001001_d+ABY*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_Py_Py_S_C1001001_d = I_ERI_D2y_S_Py_S_C1001001_d+ABY*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_Py_Py_S_C1001001_d = I_ERI_Dyz_S_Py_S_C1001001_d+ABY*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_Pz_Py_S_C1001001_d = I_ERI_Dxz_S_Py_S_C1001001_d+ABZ*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_Pz_Py_S_C1001001_d = I_ERI_Dyz_S_Py_S_C1001001_d+ABZ*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_Pz_Py_S_C1001001_d = I_ERI_D2z_S_Py_S_C1001001_d+ABZ*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_Px_Pz_S_C1001001_d = I_ERI_D2x_S_Pz_S_C1001001_d+ABX*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_Px_Pz_S_C1001001_d = I_ERI_Dxy_S_Pz_S_C1001001_d+ABX*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_Px_Pz_S_C1001001_d = I_ERI_Dxz_S_Pz_S_C1001001_d+ABX*I_ERI_Pz_S_Pz_S_C1001001_d;
  Double I_ERI_Px_Py_Pz_S_C1001001_d = I_ERI_Dxy_S_Pz_S_C1001001_d+ABY*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_Py_Pz_S_C1001001_d = I_ERI_D2y_S_Pz_S_C1001001_d+ABY*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_Py_Pz_S_C1001001_d = I_ERI_Dyz_S_Pz_S_C1001001_d+ABY*I_ERI_Pz_S_Pz_S_C1001001_d;
  Double I_ERI_Px_Pz_Pz_S_C1001001_d = I_ERI_Dxz_S_Pz_S_C1001001_d+ABZ*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_Pz_Pz_S_C1001001_d = I_ERI_Dyz_S_Pz_S_C1001001_d+ABZ*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_Pz_Pz_S_C1001001_d = I_ERI_D2z_S_Pz_S_C1001001_d+ABZ*I_ERI_Pz_S_Pz_S_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1001001_d
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1001001_d
   ************************************************************/
  Double I_ERI_Px_Px_S_Px_C1001001_d = I_ERI_D2x_S_S_Px_C1001001_d+ABX*I_ERI_Px_S_S_Px_C1001001_d;
  Double I_ERI_Py_Px_S_Px_C1001001_d = I_ERI_Dxy_S_S_Px_C1001001_d+ABX*I_ERI_Py_S_S_Px_C1001001_d;
  Double I_ERI_Pz_Px_S_Px_C1001001_d = I_ERI_Dxz_S_S_Px_C1001001_d+ABX*I_ERI_Pz_S_S_Px_C1001001_d;
  Double I_ERI_Px_Py_S_Px_C1001001_d = I_ERI_Dxy_S_S_Px_C1001001_d+ABY*I_ERI_Px_S_S_Px_C1001001_d;
  Double I_ERI_Py_Py_S_Px_C1001001_d = I_ERI_D2y_S_S_Px_C1001001_d+ABY*I_ERI_Py_S_S_Px_C1001001_d;
  Double I_ERI_Pz_Py_S_Px_C1001001_d = I_ERI_Dyz_S_S_Px_C1001001_d+ABY*I_ERI_Pz_S_S_Px_C1001001_d;
  Double I_ERI_Px_Pz_S_Px_C1001001_d = I_ERI_Dxz_S_S_Px_C1001001_d+ABZ*I_ERI_Px_S_S_Px_C1001001_d;
  Double I_ERI_Py_Pz_S_Px_C1001001_d = I_ERI_Dyz_S_S_Px_C1001001_d+ABZ*I_ERI_Py_S_S_Px_C1001001_d;
  Double I_ERI_Pz_Pz_S_Px_C1001001_d = I_ERI_D2z_S_S_Px_C1001001_d+ABZ*I_ERI_Pz_S_S_Px_C1001001_d;
  Double I_ERI_Px_Px_S_Py_C1001001_d = I_ERI_D2x_S_S_Py_C1001001_d+ABX*I_ERI_Px_S_S_Py_C1001001_d;
  Double I_ERI_Py_Px_S_Py_C1001001_d = I_ERI_Dxy_S_S_Py_C1001001_d+ABX*I_ERI_Py_S_S_Py_C1001001_d;
  Double I_ERI_Pz_Px_S_Py_C1001001_d = I_ERI_Dxz_S_S_Py_C1001001_d+ABX*I_ERI_Pz_S_S_Py_C1001001_d;
  Double I_ERI_Px_Py_S_Py_C1001001_d = I_ERI_Dxy_S_S_Py_C1001001_d+ABY*I_ERI_Px_S_S_Py_C1001001_d;
  Double I_ERI_Py_Py_S_Py_C1001001_d = I_ERI_D2y_S_S_Py_C1001001_d+ABY*I_ERI_Py_S_S_Py_C1001001_d;
  Double I_ERI_Pz_Py_S_Py_C1001001_d = I_ERI_Dyz_S_S_Py_C1001001_d+ABY*I_ERI_Pz_S_S_Py_C1001001_d;
  Double I_ERI_Px_Pz_S_Py_C1001001_d = I_ERI_Dxz_S_S_Py_C1001001_d+ABZ*I_ERI_Px_S_S_Py_C1001001_d;
  Double I_ERI_Py_Pz_S_Py_C1001001_d = I_ERI_Dyz_S_S_Py_C1001001_d+ABZ*I_ERI_Py_S_S_Py_C1001001_d;
  Double I_ERI_Pz_Pz_S_Py_C1001001_d = I_ERI_D2z_S_S_Py_C1001001_d+ABZ*I_ERI_Pz_S_S_Py_C1001001_d;
  Double I_ERI_Px_Px_S_Pz_C1001001_d = I_ERI_D2x_S_S_Pz_C1001001_d+ABX*I_ERI_Px_S_S_Pz_C1001001_d;
  Double I_ERI_Py_Px_S_Pz_C1001001_d = I_ERI_Dxy_S_S_Pz_C1001001_d+ABX*I_ERI_Py_S_S_Pz_C1001001_d;
  Double I_ERI_Pz_Px_S_Pz_C1001001_d = I_ERI_Dxz_S_S_Pz_C1001001_d+ABX*I_ERI_Pz_S_S_Pz_C1001001_d;
  Double I_ERI_Px_Py_S_Pz_C1001001_d = I_ERI_Dxy_S_S_Pz_C1001001_d+ABY*I_ERI_Px_S_S_Pz_C1001001_d;
  Double I_ERI_Py_Py_S_Pz_C1001001_d = I_ERI_D2y_S_S_Pz_C1001001_d+ABY*I_ERI_Py_S_S_Pz_C1001001_d;
  Double I_ERI_Pz_Py_S_Pz_C1001001_d = I_ERI_Dyz_S_S_Pz_C1001001_d+ABY*I_ERI_Pz_S_S_Pz_C1001001_d;
  Double I_ERI_Px_Pz_S_Pz_C1001001_d = I_ERI_Dxz_S_S_Pz_C1001001_d+ABZ*I_ERI_Px_S_S_Pz_C1001001_d;
  Double I_ERI_Py_Pz_S_Pz_C1001001_d = I_ERI_Dyz_S_S_Pz_C1001001_d+ABZ*I_ERI_Py_S_S_Pz_C1001001_d;
  Double I_ERI_Pz_Pz_S_Pz_C1001001_d = I_ERI_D2z_S_S_Pz_C1001001_d+ABZ*I_ERI_Pz_S_S_Pz_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001_aa
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1001001_aa = I_ERI_G4x_S_Px_S_C1001001_aa+ABX*I_ERI_F3x_S_Px_S_C1001001_aa;
  Double I_ERI_F2xy_Px_Px_S_C1001001_aa = I_ERI_G3xy_S_Px_S_C1001001_aa+ABX*I_ERI_F2xy_S_Px_S_C1001001_aa;
  Double I_ERI_F2xz_Px_Px_S_C1001001_aa = I_ERI_G3xz_S_Px_S_C1001001_aa+ABX*I_ERI_F2xz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2y_Px_Px_S_C1001001_aa = I_ERI_G2x2y_S_Px_S_C1001001_aa+ABX*I_ERI_Fx2y_S_Px_S_C1001001_aa;
  Double I_ERI_Fxyz_Px_Px_S_C1001001_aa = I_ERI_G2xyz_S_Px_S_C1001001_aa+ABX*I_ERI_Fxyz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2z_Px_Px_S_C1001001_aa = I_ERI_G2x2z_S_Px_S_C1001001_aa+ABX*I_ERI_Fx2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3y_Px_Px_S_C1001001_aa = I_ERI_Gx3y_S_Px_S_C1001001_aa+ABX*I_ERI_F3y_S_Px_S_C1001001_aa;
  Double I_ERI_F2yz_Px_Px_S_C1001001_aa = I_ERI_Gx2yz_S_Px_S_C1001001_aa+ABX*I_ERI_F2yz_S_Px_S_C1001001_aa;
  Double I_ERI_Fy2z_Px_Px_S_C1001001_aa = I_ERI_Gxy2z_S_Px_S_C1001001_aa+ABX*I_ERI_Fy2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3z_Px_Px_S_C1001001_aa = I_ERI_Gx3z_S_Px_S_C1001001_aa+ABX*I_ERI_F3z_S_Px_S_C1001001_aa;
  Double I_ERI_F3x_Py_Px_S_C1001001_aa = I_ERI_G3xy_S_Px_S_C1001001_aa+ABY*I_ERI_F3x_S_Px_S_C1001001_aa;
  Double I_ERI_F2xy_Py_Px_S_C1001001_aa = I_ERI_G2x2y_S_Px_S_C1001001_aa+ABY*I_ERI_F2xy_S_Px_S_C1001001_aa;
  Double I_ERI_F2xz_Py_Px_S_C1001001_aa = I_ERI_G2xyz_S_Px_S_C1001001_aa+ABY*I_ERI_F2xz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2y_Py_Px_S_C1001001_aa = I_ERI_Gx3y_S_Px_S_C1001001_aa+ABY*I_ERI_Fx2y_S_Px_S_C1001001_aa;
  Double I_ERI_Fxyz_Py_Px_S_C1001001_aa = I_ERI_Gx2yz_S_Px_S_C1001001_aa+ABY*I_ERI_Fxyz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2z_Py_Px_S_C1001001_aa = I_ERI_Gxy2z_S_Px_S_C1001001_aa+ABY*I_ERI_Fx2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3y_Py_Px_S_C1001001_aa = I_ERI_G4y_S_Px_S_C1001001_aa+ABY*I_ERI_F3y_S_Px_S_C1001001_aa;
  Double I_ERI_F2yz_Py_Px_S_C1001001_aa = I_ERI_G3yz_S_Px_S_C1001001_aa+ABY*I_ERI_F2yz_S_Px_S_C1001001_aa;
  Double I_ERI_Fy2z_Py_Px_S_C1001001_aa = I_ERI_G2y2z_S_Px_S_C1001001_aa+ABY*I_ERI_Fy2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3z_Py_Px_S_C1001001_aa = I_ERI_Gy3z_S_Px_S_C1001001_aa+ABY*I_ERI_F3z_S_Px_S_C1001001_aa;
  Double I_ERI_F3x_Pz_Px_S_C1001001_aa = I_ERI_G3xz_S_Px_S_C1001001_aa+ABZ*I_ERI_F3x_S_Px_S_C1001001_aa;
  Double I_ERI_F2xy_Pz_Px_S_C1001001_aa = I_ERI_G2xyz_S_Px_S_C1001001_aa+ABZ*I_ERI_F2xy_S_Px_S_C1001001_aa;
  Double I_ERI_F2xz_Pz_Px_S_C1001001_aa = I_ERI_G2x2z_S_Px_S_C1001001_aa+ABZ*I_ERI_F2xz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2y_Pz_Px_S_C1001001_aa = I_ERI_Gx2yz_S_Px_S_C1001001_aa+ABZ*I_ERI_Fx2y_S_Px_S_C1001001_aa;
  Double I_ERI_Fxyz_Pz_Px_S_C1001001_aa = I_ERI_Gxy2z_S_Px_S_C1001001_aa+ABZ*I_ERI_Fxyz_S_Px_S_C1001001_aa;
  Double I_ERI_Fx2z_Pz_Px_S_C1001001_aa = I_ERI_Gx3z_S_Px_S_C1001001_aa+ABZ*I_ERI_Fx2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3y_Pz_Px_S_C1001001_aa = I_ERI_G3yz_S_Px_S_C1001001_aa+ABZ*I_ERI_F3y_S_Px_S_C1001001_aa;
  Double I_ERI_F2yz_Pz_Px_S_C1001001_aa = I_ERI_G2y2z_S_Px_S_C1001001_aa+ABZ*I_ERI_F2yz_S_Px_S_C1001001_aa;
  Double I_ERI_Fy2z_Pz_Px_S_C1001001_aa = I_ERI_Gy3z_S_Px_S_C1001001_aa+ABZ*I_ERI_Fy2z_S_Px_S_C1001001_aa;
  Double I_ERI_F3z_Pz_Px_S_C1001001_aa = I_ERI_G4z_S_Px_S_C1001001_aa+ABZ*I_ERI_F3z_S_Px_S_C1001001_aa;
  Double I_ERI_F3x_Px_Py_S_C1001001_aa = I_ERI_G4x_S_Py_S_C1001001_aa+ABX*I_ERI_F3x_S_Py_S_C1001001_aa;
  Double I_ERI_F2xy_Px_Py_S_C1001001_aa = I_ERI_G3xy_S_Py_S_C1001001_aa+ABX*I_ERI_F2xy_S_Py_S_C1001001_aa;
  Double I_ERI_F2xz_Px_Py_S_C1001001_aa = I_ERI_G3xz_S_Py_S_C1001001_aa+ABX*I_ERI_F2xz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2y_Px_Py_S_C1001001_aa = I_ERI_G2x2y_S_Py_S_C1001001_aa+ABX*I_ERI_Fx2y_S_Py_S_C1001001_aa;
  Double I_ERI_Fxyz_Px_Py_S_C1001001_aa = I_ERI_G2xyz_S_Py_S_C1001001_aa+ABX*I_ERI_Fxyz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2z_Px_Py_S_C1001001_aa = I_ERI_G2x2z_S_Py_S_C1001001_aa+ABX*I_ERI_Fx2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3y_Px_Py_S_C1001001_aa = I_ERI_Gx3y_S_Py_S_C1001001_aa+ABX*I_ERI_F3y_S_Py_S_C1001001_aa;
  Double I_ERI_F2yz_Px_Py_S_C1001001_aa = I_ERI_Gx2yz_S_Py_S_C1001001_aa+ABX*I_ERI_F2yz_S_Py_S_C1001001_aa;
  Double I_ERI_Fy2z_Px_Py_S_C1001001_aa = I_ERI_Gxy2z_S_Py_S_C1001001_aa+ABX*I_ERI_Fy2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3z_Px_Py_S_C1001001_aa = I_ERI_Gx3z_S_Py_S_C1001001_aa+ABX*I_ERI_F3z_S_Py_S_C1001001_aa;
  Double I_ERI_F3x_Py_Py_S_C1001001_aa = I_ERI_G3xy_S_Py_S_C1001001_aa+ABY*I_ERI_F3x_S_Py_S_C1001001_aa;
  Double I_ERI_F2xy_Py_Py_S_C1001001_aa = I_ERI_G2x2y_S_Py_S_C1001001_aa+ABY*I_ERI_F2xy_S_Py_S_C1001001_aa;
  Double I_ERI_F2xz_Py_Py_S_C1001001_aa = I_ERI_G2xyz_S_Py_S_C1001001_aa+ABY*I_ERI_F2xz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2y_Py_Py_S_C1001001_aa = I_ERI_Gx3y_S_Py_S_C1001001_aa+ABY*I_ERI_Fx2y_S_Py_S_C1001001_aa;
  Double I_ERI_Fxyz_Py_Py_S_C1001001_aa = I_ERI_Gx2yz_S_Py_S_C1001001_aa+ABY*I_ERI_Fxyz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2z_Py_Py_S_C1001001_aa = I_ERI_Gxy2z_S_Py_S_C1001001_aa+ABY*I_ERI_Fx2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3y_Py_Py_S_C1001001_aa = I_ERI_G4y_S_Py_S_C1001001_aa+ABY*I_ERI_F3y_S_Py_S_C1001001_aa;
  Double I_ERI_F2yz_Py_Py_S_C1001001_aa = I_ERI_G3yz_S_Py_S_C1001001_aa+ABY*I_ERI_F2yz_S_Py_S_C1001001_aa;
  Double I_ERI_Fy2z_Py_Py_S_C1001001_aa = I_ERI_G2y2z_S_Py_S_C1001001_aa+ABY*I_ERI_Fy2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3z_Py_Py_S_C1001001_aa = I_ERI_Gy3z_S_Py_S_C1001001_aa+ABY*I_ERI_F3z_S_Py_S_C1001001_aa;
  Double I_ERI_F3x_Pz_Py_S_C1001001_aa = I_ERI_G3xz_S_Py_S_C1001001_aa+ABZ*I_ERI_F3x_S_Py_S_C1001001_aa;
  Double I_ERI_F2xy_Pz_Py_S_C1001001_aa = I_ERI_G2xyz_S_Py_S_C1001001_aa+ABZ*I_ERI_F2xy_S_Py_S_C1001001_aa;
  Double I_ERI_F2xz_Pz_Py_S_C1001001_aa = I_ERI_G2x2z_S_Py_S_C1001001_aa+ABZ*I_ERI_F2xz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2y_Pz_Py_S_C1001001_aa = I_ERI_Gx2yz_S_Py_S_C1001001_aa+ABZ*I_ERI_Fx2y_S_Py_S_C1001001_aa;
  Double I_ERI_Fxyz_Pz_Py_S_C1001001_aa = I_ERI_Gxy2z_S_Py_S_C1001001_aa+ABZ*I_ERI_Fxyz_S_Py_S_C1001001_aa;
  Double I_ERI_Fx2z_Pz_Py_S_C1001001_aa = I_ERI_Gx3z_S_Py_S_C1001001_aa+ABZ*I_ERI_Fx2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3y_Pz_Py_S_C1001001_aa = I_ERI_G3yz_S_Py_S_C1001001_aa+ABZ*I_ERI_F3y_S_Py_S_C1001001_aa;
  Double I_ERI_F2yz_Pz_Py_S_C1001001_aa = I_ERI_G2y2z_S_Py_S_C1001001_aa+ABZ*I_ERI_F2yz_S_Py_S_C1001001_aa;
  Double I_ERI_Fy2z_Pz_Py_S_C1001001_aa = I_ERI_Gy3z_S_Py_S_C1001001_aa+ABZ*I_ERI_Fy2z_S_Py_S_C1001001_aa;
  Double I_ERI_F3z_Pz_Py_S_C1001001_aa = I_ERI_G4z_S_Py_S_C1001001_aa+ABZ*I_ERI_F3z_S_Py_S_C1001001_aa;
  Double I_ERI_F3x_Px_Pz_S_C1001001_aa = I_ERI_G4x_S_Pz_S_C1001001_aa+ABX*I_ERI_F3x_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xy_Px_Pz_S_C1001001_aa = I_ERI_G3xy_S_Pz_S_C1001001_aa+ABX*I_ERI_F2xy_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xz_Px_Pz_S_C1001001_aa = I_ERI_G3xz_S_Pz_S_C1001001_aa+ABX*I_ERI_F2xz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2y_Px_Pz_S_C1001001_aa = I_ERI_G2x2y_S_Pz_S_C1001001_aa+ABX*I_ERI_Fx2y_S_Pz_S_C1001001_aa;
  Double I_ERI_Fxyz_Px_Pz_S_C1001001_aa = I_ERI_G2xyz_S_Pz_S_C1001001_aa+ABX*I_ERI_Fxyz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2z_Px_Pz_S_C1001001_aa = I_ERI_G2x2z_S_Pz_S_C1001001_aa+ABX*I_ERI_Fx2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3y_Px_Pz_S_C1001001_aa = I_ERI_Gx3y_S_Pz_S_C1001001_aa+ABX*I_ERI_F3y_S_Pz_S_C1001001_aa;
  Double I_ERI_F2yz_Px_Pz_S_C1001001_aa = I_ERI_Gx2yz_S_Pz_S_C1001001_aa+ABX*I_ERI_F2yz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fy2z_Px_Pz_S_C1001001_aa = I_ERI_Gxy2z_S_Pz_S_C1001001_aa+ABX*I_ERI_Fy2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3z_Px_Pz_S_C1001001_aa = I_ERI_Gx3z_S_Pz_S_C1001001_aa+ABX*I_ERI_F3z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3x_Py_Pz_S_C1001001_aa = I_ERI_G3xy_S_Pz_S_C1001001_aa+ABY*I_ERI_F3x_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xy_Py_Pz_S_C1001001_aa = I_ERI_G2x2y_S_Pz_S_C1001001_aa+ABY*I_ERI_F2xy_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xz_Py_Pz_S_C1001001_aa = I_ERI_G2xyz_S_Pz_S_C1001001_aa+ABY*I_ERI_F2xz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2y_Py_Pz_S_C1001001_aa = I_ERI_Gx3y_S_Pz_S_C1001001_aa+ABY*I_ERI_Fx2y_S_Pz_S_C1001001_aa;
  Double I_ERI_Fxyz_Py_Pz_S_C1001001_aa = I_ERI_Gx2yz_S_Pz_S_C1001001_aa+ABY*I_ERI_Fxyz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2z_Py_Pz_S_C1001001_aa = I_ERI_Gxy2z_S_Pz_S_C1001001_aa+ABY*I_ERI_Fx2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3y_Py_Pz_S_C1001001_aa = I_ERI_G4y_S_Pz_S_C1001001_aa+ABY*I_ERI_F3y_S_Pz_S_C1001001_aa;
  Double I_ERI_F2yz_Py_Pz_S_C1001001_aa = I_ERI_G3yz_S_Pz_S_C1001001_aa+ABY*I_ERI_F2yz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fy2z_Py_Pz_S_C1001001_aa = I_ERI_G2y2z_S_Pz_S_C1001001_aa+ABY*I_ERI_Fy2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3z_Py_Pz_S_C1001001_aa = I_ERI_Gy3z_S_Pz_S_C1001001_aa+ABY*I_ERI_F3z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3x_Pz_Pz_S_C1001001_aa = I_ERI_G3xz_S_Pz_S_C1001001_aa+ABZ*I_ERI_F3x_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xy_Pz_Pz_S_C1001001_aa = I_ERI_G2xyz_S_Pz_S_C1001001_aa+ABZ*I_ERI_F2xy_S_Pz_S_C1001001_aa;
  Double I_ERI_F2xz_Pz_Pz_S_C1001001_aa = I_ERI_G2x2z_S_Pz_S_C1001001_aa+ABZ*I_ERI_F2xz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2y_Pz_Pz_S_C1001001_aa = I_ERI_Gx2yz_S_Pz_S_C1001001_aa+ABZ*I_ERI_Fx2y_S_Pz_S_C1001001_aa;
  Double I_ERI_Fxyz_Pz_Pz_S_C1001001_aa = I_ERI_Gxy2z_S_Pz_S_C1001001_aa+ABZ*I_ERI_Fxyz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fx2z_Pz_Pz_S_C1001001_aa = I_ERI_Gx3z_S_Pz_S_C1001001_aa+ABZ*I_ERI_Fx2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3y_Pz_Pz_S_C1001001_aa = I_ERI_G3yz_S_Pz_S_C1001001_aa+ABZ*I_ERI_F3y_S_Pz_S_C1001001_aa;
  Double I_ERI_F2yz_Pz_Pz_S_C1001001_aa = I_ERI_G2y2z_S_Pz_S_C1001001_aa+ABZ*I_ERI_F2yz_S_Pz_S_C1001001_aa;
  Double I_ERI_Fy2z_Pz_Pz_S_C1001001_aa = I_ERI_Gy3z_S_Pz_S_C1001001_aa+ABZ*I_ERI_Fy2z_S_Pz_S_C1001001_aa;
  Double I_ERI_F3z_Pz_Pz_S_C1001001_aa = I_ERI_G4z_S_Pz_S_C1001001_aa+ABZ*I_ERI_F3z_S_Pz_S_C1001001_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_ac
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S_C1001001_ac = I_ERI_F3x_S_D2x_S_C1001001_ac+ABX*I_ERI_D2x_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxy_Px_D2x_S_C1001001_ac = I_ERI_F2xy_S_D2x_S_C1001001_ac+ABX*I_ERI_Dxy_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxz_Px_D2x_S_C1001001_ac = I_ERI_F2xz_S_D2x_S_C1001001_ac+ABX*I_ERI_Dxz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2y_Px_D2x_S_C1001001_ac = I_ERI_Fx2y_S_D2x_S_C1001001_ac+ABX*I_ERI_D2y_S_D2x_S_C1001001_ac;
  Double I_ERI_Dyz_Px_D2x_S_C1001001_ac = I_ERI_Fxyz_S_D2x_S_C1001001_ac+ABX*I_ERI_Dyz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2z_Px_D2x_S_C1001001_ac = I_ERI_Fx2z_S_D2x_S_C1001001_ac+ABX*I_ERI_D2z_S_D2x_S_C1001001_ac;
  Double I_ERI_D2x_Py_D2x_S_C1001001_ac = I_ERI_F2xy_S_D2x_S_C1001001_ac+ABY*I_ERI_D2x_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxy_Py_D2x_S_C1001001_ac = I_ERI_Fx2y_S_D2x_S_C1001001_ac+ABY*I_ERI_Dxy_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxz_Py_D2x_S_C1001001_ac = I_ERI_Fxyz_S_D2x_S_C1001001_ac+ABY*I_ERI_Dxz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2y_Py_D2x_S_C1001001_ac = I_ERI_F3y_S_D2x_S_C1001001_ac+ABY*I_ERI_D2y_S_D2x_S_C1001001_ac;
  Double I_ERI_Dyz_Py_D2x_S_C1001001_ac = I_ERI_F2yz_S_D2x_S_C1001001_ac+ABY*I_ERI_Dyz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2z_Py_D2x_S_C1001001_ac = I_ERI_Fy2z_S_D2x_S_C1001001_ac+ABY*I_ERI_D2z_S_D2x_S_C1001001_ac;
  Double I_ERI_D2x_Pz_D2x_S_C1001001_ac = I_ERI_F2xz_S_D2x_S_C1001001_ac+ABZ*I_ERI_D2x_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_D2x_S_C1001001_ac = I_ERI_Fxyz_S_D2x_S_C1001001_ac+ABZ*I_ERI_Dxy_S_D2x_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_D2x_S_C1001001_ac = I_ERI_Fx2z_S_D2x_S_C1001001_ac+ABZ*I_ERI_Dxz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2y_Pz_D2x_S_C1001001_ac = I_ERI_F2yz_S_D2x_S_C1001001_ac+ABZ*I_ERI_D2y_S_D2x_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_D2x_S_C1001001_ac = I_ERI_Fy2z_S_D2x_S_C1001001_ac+ABZ*I_ERI_Dyz_S_D2x_S_C1001001_ac;
  Double I_ERI_D2z_Pz_D2x_S_C1001001_ac = I_ERI_F3z_S_D2x_S_C1001001_ac+ABZ*I_ERI_D2z_S_D2x_S_C1001001_ac;
  Double I_ERI_D2x_Px_Dxy_S_C1001001_ac = I_ERI_F3x_S_Dxy_S_C1001001_ac+ABX*I_ERI_D2x_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxy_Px_Dxy_S_C1001001_ac = I_ERI_F2xy_S_Dxy_S_C1001001_ac+ABX*I_ERI_Dxy_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxz_Px_Dxy_S_C1001001_ac = I_ERI_F2xz_S_Dxy_S_C1001001_ac+ABX*I_ERI_Dxz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2y_Px_Dxy_S_C1001001_ac = I_ERI_Fx2y_S_Dxy_S_C1001001_ac+ABX*I_ERI_D2y_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dyz_Px_Dxy_S_C1001001_ac = I_ERI_Fxyz_S_Dxy_S_C1001001_ac+ABX*I_ERI_Dyz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2z_Px_Dxy_S_C1001001_ac = I_ERI_Fx2z_S_Dxy_S_C1001001_ac+ABX*I_ERI_D2z_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2x_Py_Dxy_S_C1001001_ac = I_ERI_F2xy_S_Dxy_S_C1001001_ac+ABY*I_ERI_D2x_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxy_Py_Dxy_S_C1001001_ac = I_ERI_Fx2y_S_Dxy_S_C1001001_ac+ABY*I_ERI_Dxy_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxz_Py_Dxy_S_C1001001_ac = I_ERI_Fxyz_S_Dxy_S_C1001001_ac+ABY*I_ERI_Dxz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2y_Py_Dxy_S_C1001001_ac = I_ERI_F3y_S_Dxy_S_C1001001_ac+ABY*I_ERI_D2y_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dyz_Py_Dxy_S_C1001001_ac = I_ERI_F2yz_S_Dxy_S_C1001001_ac+ABY*I_ERI_Dyz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2z_Py_Dxy_S_C1001001_ac = I_ERI_Fy2z_S_Dxy_S_C1001001_ac+ABY*I_ERI_D2z_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2x_Pz_Dxy_S_C1001001_ac = I_ERI_F2xz_S_Dxy_S_C1001001_ac+ABZ*I_ERI_D2x_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_Dxy_S_C1001001_ac = I_ERI_Fxyz_S_Dxy_S_C1001001_ac+ABZ*I_ERI_Dxy_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_Dxy_S_C1001001_ac = I_ERI_Fx2z_S_Dxy_S_C1001001_ac+ABZ*I_ERI_Dxz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2y_Pz_Dxy_S_C1001001_ac = I_ERI_F2yz_S_Dxy_S_C1001001_ac+ABZ*I_ERI_D2y_S_Dxy_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_Dxy_S_C1001001_ac = I_ERI_Fy2z_S_Dxy_S_C1001001_ac+ABZ*I_ERI_Dyz_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2z_Pz_Dxy_S_C1001001_ac = I_ERI_F3z_S_Dxy_S_C1001001_ac+ABZ*I_ERI_D2z_S_Dxy_S_C1001001_ac;
  Double I_ERI_D2x_Px_Dxz_S_C1001001_ac = I_ERI_F3x_S_Dxz_S_C1001001_ac+ABX*I_ERI_D2x_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxy_Px_Dxz_S_C1001001_ac = I_ERI_F2xy_S_Dxz_S_C1001001_ac+ABX*I_ERI_Dxy_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxz_Px_Dxz_S_C1001001_ac = I_ERI_F2xz_S_Dxz_S_C1001001_ac+ABX*I_ERI_Dxz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2y_Px_Dxz_S_C1001001_ac = I_ERI_Fx2y_S_Dxz_S_C1001001_ac+ABX*I_ERI_D2y_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dyz_Px_Dxz_S_C1001001_ac = I_ERI_Fxyz_S_Dxz_S_C1001001_ac+ABX*I_ERI_Dyz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2z_Px_Dxz_S_C1001001_ac = I_ERI_Fx2z_S_Dxz_S_C1001001_ac+ABX*I_ERI_D2z_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2x_Py_Dxz_S_C1001001_ac = I_ERI_F2xy_S_Dxz_S_C1001001_ac+ABY*I_ERI_D2x_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxy_Py_Dxz_S_C1001001_ac = I_ERI_Fx2y_S_Dxz_S_C1001001_ac+ABY*I_ERI_Dxy_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxz_Py_Dxz_S_C1001001_ac = I_ERI_Fxyz_S_Dxz_S_C1001001_ac+ABY*I_ERI_Dxz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2y_Py_Dxz_S_C1001001_ac = I_ERI_F3y_S_Dxz_S_C1001001_ac+ABY*I_ERI_D2y_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dyz_Py_Dxz_S_C1001001_ac = I_ERI_F2yz_S_Dxz_S_C1001001_ac+ABY*I_ERI_Dyz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2z_Py_Dxz_S_C1001001_ac = I_ERI_Fy2z_S_Dxz_S_C1001001_ac+ABY*I_ERI_D2z_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2x_Pz_Dxz_S_C1001001_ac = I_ERI_F2xz_S_Dxz_S_C1001001_ac+ABZ*I_ERI_D2x_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_Dxz_S_C1001001_ac = I_ERI_Fxyz_S_Dxz_S_C1001001_ac+ABZ*I_ERI_Dxy_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_Dxz_S_C1001001_ac = I_ERI_Fx2z_S_Dxz_S_C1001001_ac+ABZ*I_ERI_Dxz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2y_Pz_Dxz_S_C1001001_ac = I_ERI_F2yz_S_Dxz_S_C1001001_ac+ABZ*I_ERI_D2y_S_Dxz_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_Dxz_S_C1001001_ac = I_ERI_Fy2z_S_Dxz_S_C1001001_ac+ABZ*I_ERI_Dyz_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2z_Pz_Dxz_S_C1001001_ac = I_ERI_F3z_S_Dxz_S_C1001001_ac+ABZ*I_ERI_D2z_S_Dxz_S_C1001001_ac;
  Double I_ERI_D2x_Px_D2y_S_C1001001_ac = I_ERI_F3x_S_D2y_S_C1001001_ac+ABX*I_ERI_D2x_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxy_Px_D2y_S_C1001001_ac = I_ERI_F2xy_S_D2y_S_C1001001_ac+ABX*I_ERI_Dxy_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxz_Px_D2y_S_C1001001_ac = I_ERI_F2xz_S_D2y_S_C1001001_ac+ABX*I_ERI_Dxz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2y_Px_D2y_S_C1001001_ac = I_ERI_Fx2y_S_D2y_S_C1001001_ac+ABX*I_ERI_D2y_S_D2y_S_C1001001_ac;
  Double I_ERI_Dyz_Px_D2y_S_C1001001_ac = I_ERI_Fxyz_S_D2y_S_C1001001_ac+ABX*I_ERI_Dyz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2z_Px_D2y_S_C1001001_ac = I_ERI_Fx2z_S_D2y_S_C1001001_ac+ABX*I_ERI_D2z_S_D2y_S_C1001001_ac;
  Double I_ERI_D2x_Py_D2y_S_C1001001_ac = I_ERI_F2xy_S_D2y_S_C1001001_ac+ABY*I_ERI_D2x_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxy_Py_D2y_S_C1001001_ac = I_ERI_Fx2y_S_D2y_S_C1001001_ac+ABY*I_ERI_Dxy_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxz_Py_D2y_S_C1001001_ac = I_ERI_Fxyz_S_D2y_S_C1001001_ac+ABY*I_ERI_Dxz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2y_Py_D2y_S_C1001001_ac = I_ERI_F3y_S_D2y_S_C1001001_ac+ABY*I_ERI_D2y_S_D2y_S_C1001001_ac;
  Double I_ERI_Dyz_Py_D2y_S_C1001001_ac = I_ERI_F2yz_S_D2y_S_C1001001_ac+ABY*I_ERI_Dyz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2z_Py_D2y_S_C1001001_ac = I_ERI_Fy2z_S_D2y_S_C1001001_ac+ABY*I_ERI_D2z_S_D2y_S_C1001001_ac;
  Double I_ERI_D2x_Pz_D2y_S_C1001001_ac = I_ERI_F2xz_S_D2y_S_C1001001_ac+ABZ*I_ERI_D2x_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_D2y_S_C1001001_ac = I_ERI_Fxyz_S_D2y_S_C1001001_ac+ABZ*I_ERI_Dxy_S_D2y_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_D2y_S_C1001001_ac = I_ERI_Fx2z_S_D2y_S_C1001001_ac+ABZ*I_ERI_Dxz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2y_Pz_D2y_S_C1001001_ac = I_ERI_F2yz_S_D2y_S_C1001001_ac+ABZ*I_ERI_D2y_S_D2y_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_D2y_S_C1001001_ac = I_ERI_Fy2z_S_D2y_S_C1001001_ac+ABZ*I_ERI_Dyz_S_D2y_S_C1001001_ac;
  Double I_ERI_D2z_Pz_D2y_S_C1001001_ac = I_ERI_F3z_S_D2y_S_C1001001_ac+ABZ*I_ERI_D2z_S_D2y_S_C1001001_ac;
  Double I_ERI_D2x_Px_Dyz_S_C1001001_ac = I_ERI_F3x_S_Dyz_S_C1001001_ac+ABX*I_ERI_D2x_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxy_Px_Dyz_S_C1001001_ac = I_ERI_F2xy_S_Dyz_S_C1001001_ac+ABX*I_ERI_Dxy_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxz_Px_Dyz_S_C1001001_ac = I_ERI_F2xz_S_Dyz_S_C1001001_ac+ABX*I_ERI_Dxz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2y_Px_Dyz_S_C1001001_ac = I_ERI_Fx2y_S_Dyz_S_C1001001_ac+ABX*I_ERI_D2y_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dyz_Px_Dyz_S_C1001001_ac = I_ERI_Fxyz_S_Dyz_S_C1001001_ac+ABX*I_ERI_Dyz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2z_Px_Dyz_S_C1001001_ac = I_ERI_Fx2z_S_Dyz_S_C1001001_ac+ABX*I_ERI_D2z_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2x_Py_Dyz_S_C1001001_ac = I_ERI_F2xy_S_Dyz_S_C1001001_ac+ABY*I_ERI_D2x_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxy_Py_Dyz_S_C1001001_ac = I_ERI_Fx2y_S_Dyz_S_C1001001_ac+ABY*I_ERI_Dxy_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxz_Py_Dyz_S_C1001001_ac = I_ERI_Fxyz_S_Dyz_S_C1001001_ac+ABY*I_ERI_Dxz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2y_Py_Dyz_S_C1001001_ac = I_ERI_F3y_S_Dyz_S_C1001001_ac+ABY*I_ERI_D2y_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dyz_Py_Dyz_S_C1001001_ac = I_ERI_F2yz_S_Dyz_S_C1001001_ac+ABY*I_ERI_Dyz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2z_Py_Dyz_S_C1001001_ac = I_ERI_Fy2z_S_Dyz_S_C1001001_ac+ABY*I_ERI_D2z_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2x_Pz_Dyz_S_C1001001_ac = I_ERI_F2xz_S_Dyz_S_C1001001_ac+ABZ*I_ERI_D2x_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_Dyz_S_C1001001_ac = I_ERI_Fxyz_S_Dyz_S_C1001001_ac+ABZ*I_ERI_Dxy_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_Dyz_S_C1001001_ac = I_ERI_Fx2z_S_Dyz_S_C1001001_ac+ABZ*I_ERI_Dxz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2y_Pz_Dyz_S_C1001001_ac = I_ERI_F2yz_S_Dyz_S_C1001001_ac+ABZ*I_ERI_D2y_S_Dyz_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_Dyz_S_C1001001_ac = I_ERI_Fy2z_S_Dyz_S_C1001001_ac+ABZ*I_ERI_Dyz_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2z_Pz_Dyz_S_C1001001_ac = I_ERI_F3z_S_Dyz_S_C1001001_ac+ABZ*I_ERI_D2z_S_Dyz_S_C1001001_ac;
  Double I_ERI_D2x_Px_D2z_S_C1001001_ac = I_ERI_F3x_S_D2z_S_C1001001_ac+ABX*I_ERI_D2x_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxy_Px_D2z_S_C1001001_ac = I_ERI_F2xy_S_D2z_S_C1001001_ac+ABX*I_ERI_Dxy_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxz_Px_D2z_S_C1001001_ac = I_ERI_F2xz_S_D2z_S_C1001001_ac+ABX*I_ERI_Dxz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2y_Px_D2z_S_C1001001_ac = I_ERI_Fx2y_S_D2z_S_C1001001_ac+ABX*I_ERI_D2y_S_D2z_S_C1001001_ac;
  Double I_ERI_Dyz_Px_D2z_S_C1001001_ac = I_ERI_Fxyz_S_D2z_S_C1001001_ac+ABX*I_ERI_Dyz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2z_Px_D2z_S_C1001001_ac = I_ERI_Fx2z_S_D2z_S_C1001001_ac+ABX*I_ERI_D2z_S_D2z_S_C1001001_ac;
  Double I_ERI_D2x_Py_D2z_S_C1001001_ac = I_ERI_F2xy_S_D2z_S_C1001001_ac+ABY*I_ERI_D2x_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxy_Py_D2z_S_C1001001_ac = I_ERI_Fx2y_S_D2z_S_C1001001_ac+ABY*I_ERI_Dxy_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxz_Py_D2z_S_C1001001_ac = I_ERI_Fxyz_S_D2z_S_C1001001_ac+ABY*I_ERI_Dxz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2y_Py_D2z_S_C1001001_ac = I_ERI_F3y_S_D2z_S_C1001001_ac+ABY*I_ERI_D2y_S_D2z_S_C1001001_ac;
  Double I_ERI_Dyz_Py_D2z_S_C1001001_ac = I_ERI_F2yz_S_D2z_S_C1001001_ac+ABY*I_ERI_Dyz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2z_Py_D2z_S_C1001001_ac = I_ERI_Fy2z_S_D2z_S_C1001001_ac+ABY*I_ERI_D2z_S_D2z_S_C1001001_ac;
  Double I_ERI_D2x_Pz_D2z_S_C1001001_ac = I_ERI_F2xz_S_D2z_S_C1001001_ac+ABZ*I_ERI_D2x_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxy_Pz_D2z_S_C1001001_ac = I_ERI_Fxyz_S_D2z_S_C1001001_ac+ABZ*I_ERI_Dxy_S_D2z_S_C1001001_ac;
  Double I_ERI_Dxz_Pz_D2z_S_C1001001_ac = I_ERI_Fx2z_S_D2z_S_C1001001_ac+ABZ*I_ERI_Dxz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2y_Pz_D2z_S_C1001001_ac = I_ERI_F2yz_S_D2z_S_C1001001_ac+ABZ*I_ERI_D2y_S_D2z_S_C1001001_ac;
  Double I_ERI_Dyz_Pz_D2z_S_C1001001_ac = I_ERI_Fy2z_S_D2z_S_C1001001_ac+ABZ*I_ERI_Dyz_S_D2z_S_C1001001_ac;
  Double I_ERI_D2z_Pz_D2z_S_C1001001_ac = I_ERI_F3z_S_D2z_S_C1001001_ac+ABZ*I_ERI_D2z_S_D2z_S_C1001001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001_ad
   ************************************************************/
  Double I_ERI_D2x_Px_Px_Px_C1001001_ad = I_ERI_F3x_S_Px_Px_C1001001_ad+ABX*I_ERI_D2x_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxy_Px_Px_Px_C1001001_ad = I_ERI_F2xy_S_Px_Px_C1001001_ad+ABX*I_ERI_Dxy_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxz_Px_Px_Px_C1001001_ad = I_ERI_F2xz_S_Px_Px_C1001001_ad+ABX*I_ERI_Dxz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2y_Px_Px_Px_C1001001_ad = I_ERI_Fx2y_S_Px_Px_C1001001_ad+ABX*I_ERI_D2y_S_Px_Px_C1001001_ad;
  Double I_ERI_Dyz_Px_Px_Px_C1001001_ad = I_ERI_Fxyz_S_Px_Px_C1001001_ad+ABX*I_ERI_Dyz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2z_Px_Px_Px_C1001001_ad = I_ERI_Fx2z_S_Px_Px_C1001001_ad+ABX*I_ERI_D2z_S_Px_Px_C1001001_ad;
  Double I_ERI_D2x_Py_Px_Px_C1001001_ad = I_ERI_F2xy_S_Px_Px_C1001001_ad+ABY*I_ERI_D2x_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxy_Py_Px_Px_C1001001_ad = I_ERI_Fx2y_S_Px_Px_C1001001_ad+ABY*I_ERI_Dxy_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxz_Py_Px_Px_C1001001_ad = I_ERI_Fxyz_S_Px_Px_C1001001_ad+ABY*I_ERI_Dxz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2y_Py_Px_Px_C1001001_ad = I_ERI_F3y_S_Px_Px_C1001001_ad+ABY*I_ERI_D2y_S_Px_Px_C1001001_ad;
  Double I_ERI_Dyz_Py_Px_Px_C1001001_ad = I_ERI_F2yz_S_Px_Px_C1001001_ad+ABY*I_ERI_Dyz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2z_Py_Px_Px_C1001001_ad = I_ERI_Fy2z_S_Px_Px_C1001001_ad+ABY*I_ERI_D2z_S_Px_Px_C1001001_ad;
  Double I_ERI_D2x_Pz_Px_Px_C1001001_ad = I_ERI_F2xz_S_Px_Px_C1001001_ad+ABZ*I_ERI_D2x_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxy_Pz_Px_Px_C1001001_ad = I_ERI_Fxyz_S_Px_Px_C1001001_ad+ABZ*I_ERI_Dxy_S_Px_Px_C1001001_ad;
  Double I_ERI_Dxz_Pz_Px_Px_C1001001_ad = I_ERI_Fx2z_S_Px_Px_C1001001_ad+ABZ*I_ERI_Dxz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2y_Pz_Px_Px_C1001001_ad = I_ERI_F2yz_S_Px_Px_C1001001_ad+ABZ*I_ERI_D2y_S_Px_Px_C1001001_ad;
  Double I_ERI_Dyz_Pz_Px_Px_C1001001_ad = I_ERI_Fy2z_S_Px_Px_C1001001_ad+ABZ*I_ERI_Dyz_S_Px_Px_C1001001_ad;
  Double I_ERI_D2z_Pz_Px_Px_C1001001_ad = I_ERI_F3z_S_Px_Px_C1001001_ad+ABZ*I_ERI_D2z_S_Px_Px_C1001001_ad;
  Double I_ERI_D2x_Px_Py_Px_C1001001_ad = I_ERI_F3x_S_Py_Px_C1001001_ad+ABX*I_ERI_D2x_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxy_Px_Py_Px_C1001001_ad = I_ERI_F2xy_S_Py_Px_C1001001_ad+ABX*I_ERI_Dxy_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxz_Px_Py_Px_C1001001_ad = I_ERI_F2xz_S_Py_Px_C1001001_ad+ABX*I_ERI_Dxz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2y_Px_Py_Px_C1001001_ad = I_ERI_Fx2y_S_Py_Px_C1001001_ad+ABX*I_ERI_D2y_S_Py_Px_C1001001_ad;
  Double I_ERI_Dyz_Px_Py_Px_C1001001_ad = I_ERI_Fxyz_S_Py_Px_C1001001_ad+ABX*I_ERI_Dyz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2z_Px_Py_Px_C1001001_ad = I_ERI_Fx2z_S_Py_Px_C1001001_ad+ABX*I_ERI_D2z_S_Py_Px_C1001001_ad;
  Double I_ERI_D2x_Py_Py_Px_C1001001_ad = I_ERI_F2xy_S_Py_Px_C1001001_ad+ABY*I_ERI_D2x_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxy_Py_Py_Px_C1001001_ad = I_ERI_Fx2y_S_Py_Px_C1001001_ad+ABY*I_ERI_Dxy_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxz_Py_Py_Px_C1001001_ad = I_ERI_Fxyz_S_Py_Px_C1001001_ad+ABY*I_ERI_Dxz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2y_Py_Py_Px_C1001001_ad = I_ERI_F3y_S_Py_Px_C1001001_ad+ABY*I_ERI_D2y_S_Py_Px_C1001001_ad;
  Double I_ERI_Dyz_Py_Py_Px_C1001001_ad = I_ERI_F2yz_S_Py_Px_C1001001_ad+ABY*I_ERI_Dyz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2z_Py_Py_Px_C1001001_ad = I_ERI_Fy2z_S_Py_Px_C1001001_ad+ABY*I_ERI_D2z_S_Py_Px_C1001001_ad;
  Double I_ERI_D2x_Pz_Py_Px_C1001001_ad = I_ERI_F2xz_S_Py_Px_C1001001_ad+ABZ*I_ERI_D2x_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxy_Pz_Py_Px_C1001001_ad = I_ERI_Fxyz_S_Py_Px_C1001001_ad+ABZ*I_ERI_Dxy_S_Py_Px_C1001001_ad;
  Double I_ERI_Dxz_Pz_Py_Px_C1001001_ad = I_ERI_Fx2z_S_Py_Px_C1001001_ad+ABZ*I_ERI_Dxz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2y_Pz_Py_Px_C1001001_ad = I_ERI_F2yz_S_Py_Px_C1001001_ad+ABZ*I_ERI_D2y_S_Py_Px_C1001001_ad;
  Double I_ERI_Dyz_Pz_Py_Px_C1001001_ad = I_ERI_Fy2z_S_Py_Px_C1001001_ad+ABZ*I_ERI_Dyz_S_Py_Px_C1001001_ad;
  Double I_ERI_D2z_Pz_Py_Px_C1001001_ad = I_ERI_F3z_S_Py_Px_C1001001_ad+ABZ*I_ERI_D2z_S_Py_Px_C1001001_ad;
  Double I_ERI_D2x_Px_Pz_Px_C1001001_ad = I_ERI_F3x_S_Pz_Px_C1001001_ad+ABX*I_ERI_D2x_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxy_Px_Pz_Px_C1001001_ad = I_ERI_F2xy_S_Pz_Px_C1001001_ad+ABX*I_ERI_Dxy_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxz_Px_Pz_Px_C1001001_ad = I_ERI_F2xz_S_Pz_Px_C1001001_ad+ABX*I_ERI_Dxz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2y_Px_Pz_Px_C1001001_ad = I_ERI_Fx2y_S_Pz_Px_C1001001_ad+ABX*I_ERI_D2y_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dyz_Px_Pz_Px_C1001001_ad = I_ERI_Fxyz_S_Pz_Px_C1001001_ad+ABX*I_ERI_Dyz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2z_Px_Pz_Px_C1001001_ad = I_ERI_Fx2z_S_Pz_Px_C1001001_ad+ABX*I_ERI_D2z_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2x_Py_Pz_Px_C1001001_ad = I_ERI_F2xy_S_Pz_Px_C1001001_ad+ABY*I_ERI_D2x_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxy_Py_Pz_Px_C1001001_ad = I_ERI_Fx2y_S_Pz_Px_C1001001_ad+ABY*I_ERI_Dxy_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxz_Py_Pz_Px_C1001001_ad = I_ERI_Fxyz_S_Pz_Px_C1001001_ad+ABY*I_ERI_Dxz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2y_Py_Pz_Px_C1001001_ad = I_ERI_F3y_S_Pz_Px_C1001001_ad+ABY*I_ERI_D2y_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dyz_Py_Pz_Px_C1001001_ad = I_ERI_F2yz_S_Pz_Px_C1001001_ad+ABY*I_ERI_Dyz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2z_Py_Pz_Px_C1001001_ad = I_ERI_Fy2z_S_Pz_Px_C1001001_ad+ABY*I_ERI_D2z_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2x_Pz_Pz_Px_C1001001_ad = I_ERI_F2xz_S_Pz_Px_C1001001_ad+ABZ*I_ERI_D2x_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxy_Pz_Pz_Px_C1001001_ad = I_ERI_Fxyz_S_Pz_Px_C1001001_ad+ABZ*I_ERI_Dxy_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dxz_Pz_Pz_Px_C1001001_ad = I_ERI_Fx2z_S_Pz_Px_C1001001_ad+ABZ*I_ERI_Dxz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2y_Pz_Pz_Px_C1001001_ad = I_ERI_F2yz_S_Pz_Px_C1001001_ad+ABZ*I_ERI_D2y_S_Pz_Px_C1001001_ad;
  Double I_ERI_Dyz_Pz_Pz_Px_C1001001_ad = I_ERI_Fy2z_S_Pz_Px_C1001001_ad+ABZ*I_ERI_Dyz_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2z_Pz_Pz_Px_C1001001_ad = I_ERI_F3z_S_Pz_Px_C1001001_ad+ABZ*I_ERI_D2z_S_Pz_Px_C1001001_ad;
  Double I_ERI_D2x_Px_Px_Py_C1001001_ad = I_ERI_F3x_S_Px_Py_C1001001_ad+ABX*I_ERI_D2x_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxy_Px_Px_Py_C1001001_ad = I_ERI_F2xy_S_Px_Py_C1001001_ad+ABX*I_ERI_Dxy_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxz_Px_Px_Py_C1001001_ad = I_ERI_F2xz_S_Px_Py_C1001001_ad+ABX*I_ERI_Dxz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2y_Px_Px_Py_C1001001_ad = I_ERI_Fx2y_S_Px_Py_C1001001_ad+ABX*I_ERI_D2y_S_Px_Py_C1001001_ad;
  Double I_ERI_Dyz_Px_Px_Py_C1001001_ad = I_ERI_Fxyz_S_Px_Py_C1001001_ad+ABX*I_ERI_Dyz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2z_Px_Px_Py_C1001001_ad = I_ERI_Fx2z_S_Px_Py_C1001001_ad+ABX*I_ERI_D2z_S_Px_Py_C1001001_ad;
  Double I_ERI_D2x_Py_Px_Py_C1001001_ad = I_ERI_F2xy_S_Px_Py_C1001001_ad+ABY*I_ERI_D2x_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxy_Py_Px_Py_C1001001_ad = I_ERI_Fx2y_S_Px_Py_C1001001_ad+ABY*I_ERI_Dxy_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxz_Py_Px_Py_C1001001_ad = I_ERI_Fxyz_S_Px_Py_C1001001_ad+ABY*I_ERI_Dxz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2y_Py_Px_Py_C1001001_ad = I_ERI_F3y_S_Px_Py_C1001001_ad+ABY*I_ERI_D2y_S_Px_Py_C1001001_ad;
  Double I_ERI_Dyz_Py_Px_Py_C1001001_ad = I_ERI_F2yz_S_Px_Py_C1001001_ad+ABY*I_ERI_Dyz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2z_Py_Px_Py_C1001001_ad = I_ERI_Fy2z_S_Px_Py_C1001001_ad+ABY*I_ERI_D2z_S_Px_Py_C1001001_ad;
  Double I_ERI_D2x_Pz_Px_Py_C1001001_ad = I_ERI_F2xz_S_Px_Py_C1001001_ad+ABZ*I_ERI_D2x_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxy_Pz_Px_Py_C1001001_ad = I_ERI_Fxyz_S_Px_Py_C1001001_ad+ABZ*I_ERI_Dxy_S_Px_Py_C1001001_ad;
  Double I_ERI_Dxz_Pz_Px_Py_C1001001_ad = I_ERI_Fx2z_S_Px_Py_C1001001_ad+ABZ*I_ERI_Dxz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2y_Pz_Px_Py_C1001001_ad = I_ERI_F2yz_S_Px_Py_C1001001_ad+ABZ*I_ERI_D2y_S_Px_Py_C1001001_ad;
  Double I_ERI_Dyz_Pz_Px_Py_C1001001_ad = I_ERI_Fy2z_S_Px_Py_C1001001_ad+ABZ*I_ERI_Dyz_S_Px_Py_C1001001_ad;
  Double I_ERI_D2z_Pz_Px_Py_C1001001_ad = I_ERI_F3z_S_Px_Py_C1001001_ad+ABZ*I_ERI_D2z_S_Px_Py_C1001001_ad;
  Double I_ERI_D2x_Px_Py_Py_C1001001_ad = I_ERI_F3x_S_Py_Py_C1001001_ad+ABX*I_ERI_D2x_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxy_Px_Py_Py_C1001001_ad = I_ERI_F2xy_S_Py_Py_C1001001_ad+ABX*I_ERI_Dxy_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxz_Px_Py_Py_C1001001_ad = I_ERI_F2xz_S_Py_Py_C1001001_ad+ABX*I_ERI_Dxz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2y_Px_Py_Py_C1001001_ad = I_ERI_Fx2y_S_Py_Py_C1001001_ad+ABX*I_ERI_D2y_S_Py_Py_C1001001_ad;
  Double I_ERI_Dyz_Px_Py_Py_C1001001_ad = I_ERI_Fxyz_S_Py_Py_C1001001_ad+ABX*I_ERI_Dyz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2z_Px_Py_Py_C1001001_ad = I_ERI_Fx2z_S_Py_Py_C1001001_ad+ABX*I_ERI_D2z_S_Py_Py_C1001001_ad;
  Double I_ERI_D2x_Py_Py_Py_C1001001_ad = I_ERI_F2xy_S_Py_Py_C1001001_ad+ABY*I_ERI_D2x_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxy_Py_Py_Py_C1001001_ad = I_ERI_Fx2y_S_Py_Py_C1001001_ad+ABY*I_ERI_Dxy_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxz_Py_Py_Py_C1001001_ad = I_ERI_Fxyz_S_Py_Py_C1001001_ad+ABY*I_ERI_Dxz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2y_Py_Py_Py_C1001001_ad = I_ERI_F3y_S_Py_Py_C1001001_ad+ABY*I_ERI_D2y_S_Py_Py_C1001001_ad;
  Double I_ERI_Dyz_Py_Py_Py_C1001001_ad = I_ERI_F2yz_S_Py_Py_C1001001_ad+ABY*I_ERI_Dyz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2z_Py_Py_Py_C1001001_ad = I_ERI_Fy2z_S_Py_Py_C1001001_ad+ABY*I_ERI_D2z_S_Py_Py_C1001001_ad;
  Double I_ERI_D2x_Pz_Py_Py_C1001001_ad = I_ERI_F2xz_S_Py_Py_C1001001_ad+ABZ*I_ERI_D2x_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxy_Pz_Py_Py_C1001001_ad = I_ERI_Fxyz_S_Py_Py_C1001001_ad+ABZ*I_ERI_Dxy_S_Py_Py_C1001001_ad;
  Double I_ERI_Dxz_Pz_Py_Py_C1001001_ad = I_ERI_Fx2z_S_Py_Py_C1001001_ad+ABZ*I_ERI_Dxz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2y_Pz_Py_Py_C1001001_ad = I_ERI_F2yz_S_Py_Py_C1001001_ad+ABZ*I_ERI_D2y_S_Py_Py_C1001001_ad;
  Double I_ERI_Dyz_Pz_Py_Py_C1001001_ad = I_ERI_Fy2z_S_Py_Py_C1001001_ad+ABZ*I_ERI_Dyz_S_Py_Py_C1001001_ad;
  Double I_ERI_D2z_Pz_Py_Py_C1001001_ad = I_ERI_F3z_S_Py_Py_C1001001_ad+ABZ*I_ERI_D2z_S_Py_Py_C1001001_ad;
  Double I_ERI_D2x_Px_Pz_Py_C1001001_ad = I_ERI_F3x_S_Pz_Py_C1001001_ad+ABX*I_ERI_D2x_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxy_Px_Pz_Py_C1001001_ad = I_ERI_F2xy_S_Pz_Py_C1001001_ad+ABX*I_ERI_Dxy_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxz_Px_Pz_Py_C1001001_ad = I_ERI_F2xz_S_Pz_Py_C1001001_ad+ABX*I_ERI_Dxz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2y_Px_Pz_Py_C1001001_ad = I_ERI_Fx2y_S_Pz_Py_C1001001_ad+ABX*I_ERI_D2y_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dyz_Px_Pz_Py_C1001001_ad = I_ERI_Fxyz_S_Pz_Py_C1001001_ad+ABX*I_ERI_Dyz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2z_Px_Pz_Py_C1001001_ad = I_ERI_Fx2z_S_Pz_Py_C1001001_ad+ABX*I_ERI_D2z_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2x_Py_Pz_Py_C1001001_ad = I_ERI_F2xy_S_Pz_Py_C1001001_ad+ABY*I_ERI_D2x_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxy_Py_Pz_Py_C1001001_ad = I_ERI_Fx2y_S_Pz_Py_C1001001_ad+ABY*I_ERI_Dxy_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxz_Py_Pz_Py_C1001001_ad = I_ERI_Fxyz_S_Pz_Py_C1001001_ad+ABY*I_ERI_Dxz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2y_Py_Pz_Py_C1001001_ad = I_ERI_F3y_S_Pz_Py_C1001001_ad+ABY*I_ERI_D2y_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dyz_Py_Pz_Py_C1001001_ad = I_ERI_F2yz_S_Pz_Py_C1001001_ad+ABY*I_ERI_Dyz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2z_Py_Pz_Py_C1001001_ad = I_ERI_Fy2z_S_Pz_Py_C1001001_ad+ABY*I_ERI_D2z_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2x_Pz_Pz_Py_C1001001_ad = I_ERI_F2xz_S_Pz_Py_C1001001_ad+ABZ*I_ERI_D2x_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxy_Pz_Pz_Py_C1001001_ad = I_ERI_Fxyz_S_Pz_Py_C1001001_ad+ABZ*I_ERI_Dxy_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dxz_Pz_Pz_Py_C1001001_ad = I_ERI_Fx2z_S_Pz_Py_C1001001_ad+ABZ*I_ERI_Dxz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2y_Pz_Pz_Py_C1001001_ad = I_ERI_F2yz_S_Pz_Py_C1001001_ad+ABZ*I_ERI_D2y_S_Pz_Py_C1001001_ad;
  Double I_ERI_Dyz_Pz_Pz_Py_C1001001_ad = I_ERI_Fy2z_S_Pz_Py_C1001001_ad+ABZ*I_ERI_Dyz_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2z_Pz_Pz_Py_C1001001_ad = I_ERI_F3z_S_Pz_Py_C1001001_ad+ABZ*I_ERI_D2z_S_Pz_Py_C1001001_ad;
  Double I_ERI_D2x_Px_Px_Pz_C1001001_ad = I_ERI_F3x_S_Px_Pz_C1001001_ad+ABX*I_ERI_D2x_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxy_Px_Px_Pz_C1001001_ad = I_ERI_F2xy_S_Px_Pz_C1001001_ad+ABX*I_ERI_Dxy_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxz_Px_Px_Pz_C1001001_ad = I_ERI_F2xz_S_Px_Pz_C1001001_ad+ABX*I_ERI_Dxz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2y_Px_Px_Pz_C1001001_ad = I_ERI_Fx2y_S_Px_Pz_C1001001_ad+ABX*I_ERI_D2y_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dyz_Px_Px_Pz_C1001001_ad = I_ERI_Fxyz_S_Px_Pz_C1001001_ad+ABX*I_ERI_Dyz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2z_Px_Px_Pz_C1001001_ad = I_ERI_Fx2z_S_Px_Pz_C1001001_ad+ABX*I_ERI_D2z_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2x_Py_Px_Pz_C1001001_ad = I_ERI_F2xy_S_Px_Pz_C1001001_ad+ABY*I_ERI_D2x_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxy_Py_Px_Pz_C1001001_ad = I_ERI_Fx2y_S_Px_Pz_C1001001_ad+ABY*I_ERI_Dxy_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxz_Py_Px_Pz_C1001001_ad = I_ERI_Fxyz_S_Px_Pz_C1001001_ad+ABY*I_ERI_Dxz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2y_Py_Px_Pz_C1001001_ad = I_ERI_F3y_S_Px_Pz_C1001001_ad+ABY*I_ERI_D2y_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dyz_Py_Px_Pz_C1001001_ad = I_ERI_F2yz_S_Px_Pz_C1001001_ad+ABY*I_ERI_Dyz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2z_Py_Px_Pz_C1001001_ad = I_ERI_Fy2z_S_Px_Pz_C1001001_ad+ABY*I_ERI_D2z_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2x_Pz_Px_Pz_C1001001_ad = I_ERI_F2xz_S_Px_Pz_C1001001_ad+ABZ*I_ERI_D2x_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxy_Pz_Px_Pz_C1001001_ad = I_ERI_Fxyz_S_Px_Pz_C1001001_ad+ABZ*I_ERI_Dxy_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dxz_Pz_Px_Pz_C1001001_ad = I_ERI_Fx2z_S_Px_Pz_C1001001_ad+ABZ*I_ERI_Dxz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2y_Pz_Px_Pz_C1001001_ad = I_ERI_F2yz_S_Px_Pz_C1001001_ad+ABZ*I_ERI_D2y_S_Px_Pz_C1001001_ad;
  Double I_ERI_Dyz_Pz_Px_Pz_C1001001_ad = I_ERI_Fy2z_S_Px_Pz_C1001001_ad+ABZ*I_ERI_Dyz_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2z_Pz_Px_Pz_C1001001_ad = I_ERI_F3z_S_Px_Pz_C1001001_ad+ABZ*I_ERI_D2z_S_Px_Pz_C1001001_ad;
  Double I_ERI_D2x_Px_Py_Pz_C1001001_ad = I_ERI_F3x_S_Py_Pz_C1001001_ad+ABX*I_ERI_D2x_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxy_Px_Py_Pz_C1001001_ad = I_ERI_F2xy_S_Py_Pz_C1001001_ad+ABX*I_ERI_Dxy_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxz_Px_Py_Pz_C1001001_ad = I_ERI_F2xz_S_Py_Pz_C1001001_ad+ABX*I_ERI_Dxz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2y_Px_Py_Pz_C1001001_ad = I_ERI_Fx2y_S_Py_Pz_C1001001_ad+ABX*I_ERI_D2y_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dyz_Px_Py_Pz_C1001001_ad = I_ERI_Fxyz_S_Py_Pz_C1001001_ad+ABX*I_ERI_Dyz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2z_Px_Py_Pz_C1001001_ad = I_ERI_Fx2z_S_Py_Pz_C1001001_ad+ABX*I_ERI_D2z_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2x_Py_Py_Pz_C1001001_ad = I_ERI_F2xy_S_Py_Pz_C1001001_ad+ABY*I_ERI_D2x_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxy_Py_Py_Pz_C1001001_ad = I_ERI_Fx2y_S_Py_Pz_C1001001_ad+ABY*I_ERI_Dxy_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxz_Py_Py_Pz_C1001001_ad = I_ERI_Fxyz_S_Py_Pz_C1001001_ad+ABY*I_ERI_Dxz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2y_Py_Py_Pz_C1001001_ad = I_ERI_F3y_S_Py_Pz_C1001001_ad+ABY*I_ERI_D2y_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dyz_Py_Py_Pz_C1001001_ad = I_ERI_F2yz_S_Py_Pz_C1001001_ad+ABY*I_ERI_Dyz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2z_Py_Py_Pz_C1001001_ad = I_ERI_Fy2z_S_Py_Pz_C1001001_ad+ABY*I_ERI_D2z_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2x_Pz_Py_Pz_C1001001_ad = I_ERI_F2xz_S_Py_Pz_C1001001_ad+ABZ*I_ERI_D2x_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxy_Pz_Py_Pz_C1001001_ad = I_ERI_Fxyz_S_Py_Pz_C1001001_ad+ABZ*I_ERI_Dxy_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dxz_Pz_Py_Pz_C1001001_ad = I_ERI_Fx2z_S_Py_Pz_C1001001_ad+ABZ*I_ERI_Dxz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2y_Pz_Py_Pz_C1001001_ad = I_ERI_F2yz_S_Py_Pz_C1001001_ad+ABZ*I_ERI_D2y_S_Py_Pz_C1001001_ad;
  Double I_ERI_Dyz_Pz_Py_Pz_C1001001_ad = I_ERI_Fy2z_S_Py_Pz_C1001001_ad+ABZ*I_ERI_Dyz_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2z_Pz_Py_Pz_C1001001_ad = I_ERI_F3z_S_Py_Pz_C1001001_ad+ABZ*I_ERI_D2z_S_Py_Pz_C1001001_ad;
  Double I_ERI_D2x_Px_Pz_Pz_C1001001_ad = I_ERI_F3x_S_Pz_Pz_C1001001_ad+ABX*I_ERI_D2x_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxy_Px_Pz_Pz_C1001001_ad = I_ERI_F2xy_S_Pz_Pz_C1001001_ad+ABX*I_ERI_Dxy_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxz_Px_Pz_Pz_C1001001_ad = I_ERI_F2xz_S_Pz_Pz_C1001001_ad+ABX*I_ERI_Dxz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2y_Px_Pz_Pz_C1001001_ad = I_ERI_Fx2y_S_Pz_Pz_C1001001_ad+ABX*I_ERI_D2y_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dyz_Px_Pz_Pz_C1001001_ad = I_ERI_Fxyz_S_Pz_Pz_C1001001_ad+ABX*I_ERI_Dyz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2z_Px_Pz_Pz_C1001001_ad = I_ERI_Fx2z_S_Pz_Pz_C1001001_ad+ABX*I_ERI_D2z_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2x_Py_Pz_Pz_C1001001_ad = I_ERI_F2xy_S_Pz_Pz_C1001001_ad+ABY*I_ERI_D2x_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxy_Py_Pz_Pz_C1001001_ad = I_ERI_Fx2y_S_Pz_Pz_C1001001_ad+ABY*I_ERI_Dxy_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxz_Py_Pz_Pz_C1001001_ad = I_ERI_Fxyz_S_Pz_Pz_C1001001_ad+ABY*I_ERI_Dxz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2y_Py_Pz_Pz_C1001001_ad = I_ERI_F3y_S_Pz_Pz_C1001001_ad+ABY*I_ERI_D2y_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dyz_Py_Pz_Pz_C1001001_ad = I_ERI_F2yz_S_Pz_Pz_C1001001_ad+ABY*I_ERI_Dyz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2z_Py_Pz_Pz_C1001001_ad = I_ERI_Fy2z_S_Pz_Pz_C1001001_ad+ABY*I_ERI_D2z_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2x_Pz_Pz_Pz_C1001001_ad = I_ERI_F2xz_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_D2x_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxy_Pz_Pz_Pz_C1001001_ad = I_ERI_Fxyz_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_Dxy_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dxz_Pz_Pz_Pz_C1001001_ad = I_ERI_Fx2z_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_Dxz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2y_Pz_Pz_Pz_C1001001_ad = I_ERI_F2yz_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_D2y_S_Pz_Pz_C1001001_ad;
  Double I_ERI_Dyz_Pz_Pz_Pz_C1001001_ad = I_ERI_Fy2z_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_Dyz_S_Pz_Pz_C1001001_ad;
  Double I_ERI_D2z_Pz_Pz_Pz_C1001001_ad = I_ERI_F3z_S_Pz_Pz_C1001001_ad+ABZ*I_ERI_D2z_S_Pz_Pz_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1001001_cc
   ************************************************************/
  Double I_ERI_Px_Px_F3x_S_C1001001_cc = I_ERI_D2x_S_F3x_S_C1001001_cc+ABX*I_ERI_Px_S_F3x_S_C1001001_cc;
  Double I_ERI_Py_Px_F3x_S_C1001001_cc = I_ERI_Dxy_S_F3x_S_C1001001_cc+ABX*I_ERI_Py_S_F3x_S_C1001001_cc;
  Double I_ERI_Pz_Px_F3x_S_C1001001_cc = I_ERI_Dxz_S_F3x_S_C1001001_cc+ABX*I_ERI_Pz_S_F3x_S_C1001001_cc;
  Double I_ERI_Px_Py_F3x_S_C1001001_cc = I_ERI_Dxy_S_F3x_S_C1001001_cc+ABY*I_ERI_Px_S_F3x_S_C1001001_cc;
  Double I_ERI_Py_Py_F3x_S_C1001001_cc = I_ERI_D2y_S_F3x_S_C1001001_cc+ABY*I_ERI_Py_S_F3x_S_C1001001_cc;
  Double I_ERI_Pz_Py_F3x_S_C1001001_cc = I_ERI_Dyz_S_F3x_S_C1001001_cc+ABY*I_ERI_Pz_S_F3x_S_C1001001_cc;
  Double I_ERI_Px_Pz_F3x_S_C1001001_cc = I_ERI_Dxz_S_F3x_S_C1001001_cc+ABZ*I_ERI_Px_S_F3x_S_C1001001_cc;
  Double I_ERI_Py_Pz_F3x_S_C1001001_cc = I_ERI_Dyz_S_F3x_S_C1001001_cc+ABZ*I_ERI_Py_S_F3x_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F3x_S_C1001001_cc = I_ERI_D2z_S_F3x_S_C1001001_cc+ABZ*I_ERI_Pz_S_F3x_S_C1001001_cc;
  Double I_ERI_Px_Px_F2xy_S_C1001001_cc = I_ERI_D2x_S_F2xy_S_C1001001_cc+ABX*I_ERI_Px_S_F2xy_S_C1001001_cc;
  Double I_ERI_Py_Px_F2xy_S_C1001001_cc = I_ERI_Dxy_S_F2xy_S_C1001001_cc+ABX*I_ERI_Py_S_F2xy_S_C1001001_cc;
  Double I_ERI_Pz_Px_F2xy_S_C1001001_cc = I_ERI_Dxz_S_F2xy_S_C1001001_cc+ABX*I_ERI_Pz_S_F2xy_S_C1001001_cc;
  Double I_ERI_Px_Py_F2xy_S_C1001001_cc = I_ERI_Dxy_S_F2xy_S_C1001001_cc+ABY*I_ERI_Px_S_F2xy_S_C1001001_cc;
  Double I_ERI_Py_Py_F2xy_S_C1001001_cc = I_ERI_D2y_S_F2xy_S_C1001001_cc+ABY*I_ERI_Py_S_F2xy_S_C1001001_cc;
  Double I_ERI_Pz_Py_F2xy_S_C1001001_cc = I_ERI_Dyz_S_F2xy_S_C1001001_cc+ABY*I_ERI_Pz_S_F2xy_S_C1001001_cc;
  Double I_ERI_Px_Pz_F2xy_S_C1001001_cc = I_ERI_Dxz_S_F2xy_S_C1001001_cc+ABZ*I_ERI_Px_S_F2xy_S_C1001001_cc;
  Double I_ERI_Py_Pz_F2xy_S_C1001001_cc = I_ERI_Dyz_S_F2xy_S_C1001001_cc+ABZ*I_ERI_Py_S_F2xy_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F2xy_S_C1001001_cc = I_ERI_D2z_S_F2xy_S_C1001001_cc+ABZ*I_ERI_Pz_S_F2xy_S_C1001001_cc;
  Double I_ERI_Px_Px_F2xz_S_C1001001_cc = I_ERI_D2x_S_F2xz_S_C1001001_cc+ABX*I_ERI_Px_S_F2xz_S_C1001001_cc;
  Double I_ERI_Py_Px_F2xz_S_C1001001_cc = I_ERI_Dxy_S_F2xz_S_C1001001_cc+ABX*I_ERI_Py_S_F2xz_S_C1001001_cc;
  Double I_ERI_Pz_Px_F2xz_S_C1001001_cc = I_ERI_Dxz_S_F2xz_S_C1001001_cc+ABX*I_ERI_Pz_S_F2xz_S_C1001001_cc;
  Double I_ERI_Px_Py_F2xz_S_C1001001_cc = I_ERI_Dxy_S_F2xz_S_C1001001_cc+ABY*I_ERI_Px_S_F2xz_S_C1001001_cc;
  Double I_ERI_Py_Py_F2xz_S_C1001001_cc = I_ERI_D2y_S_F2xz_S_C1001001_cc+ABY*I_ERI_Py_S_F2xz_S_C1001001_cc;
  Double I_ERI_Pz_Py_F2xz_S_C1001001_cc = I_ERI_Dyz_S_F2xz_S_C1001001_cc+ABY*I_ERI_Pz_S_F2xz_S_C1001001_cc;
  Double I_ERI_Px_Pz_F2xz_S_C1001001_cc = I_ERI_Dxz_S_F2xz_S_C1001001_cc+ABZ*I_ERI_Px_S_F2xz_S_C1001001_cc;
  Double I_ERI_Py_Pz_F2xz_S_C1001001_cc = I_ERI_Dyz_S_F2xz_S_C1001001_cc+ABZ*I_ERI_Py_S_F2xz_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F2xz_S_C1001001_cc = I_ERI_D2z_S_F2xz_S_C1001001_cc+ABZ*I_ERI_Pz_S_F2xz_S_C1001001_cc;
  Double I_ERI_Px_Px_Fx2y_S_C1001001_cc = I_ERI_D2x_S_Fx2y_S_C1001001_cc+ABX*I_ERI_Px_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Py_Px_Fx2y_S_C1001001_cc = I_ERI_Dxy_S_Fx2y_S_C1001001_cc+ABX*I_ERI_Py_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Pz_Px_Fx2y_S_C1001001_cc = I_ERI_Dxz_S_Fx2y_S_C1001001_cc+ABX*I_ERI_Pz_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Px_Py_Fx2y_S_C1001001_cc = I_ERI_Dxy_S_Fx2y_S_C1001001_cc+ABY*I_ERI_Px_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Py_Py_Fx2y_S_C1001001_cc = I_ERI_D2y_S_Fx2y_S_C1001001_cc+ABY*I_ERI_Py_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Pz_Py_Fx2y_S_C1001001_cc = I_ERI_Dyz_S_Fx2y_S_C1001001_cc+ABY*I_ERI_Pz_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Px_Pz_Fx2y_S_C1001001_cc = I_ERI_Dxz_S_Fx2y_S_C1001001_cc+ABZ*I_ERI_Px_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Py_Pz_Fx2y_S_C1001001_cc = I_ERI_Dyz_S_Fx2y_S_C1001001_cc+ABZ*I_ERI_Py_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Pz_Pz_Fx2y_S_C1001001_cc = I_ERI_D2z_S_Fx2y_S_C1001001_cc+ABZ*I_ERI_Pz_S_Fx2y_S_C1001001_cc;
  Double I_ERI_Px_Px_Fxyz_S_C1001001_cc = I_ERI_D2x_S_Fxyz_S_C1001001_cc+ABX*I_ERI_Px_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Py_Px_Fxyz_S_C1001001_cc = I_ERI_Dxy_S_Fxyz_S_C1001001_cc+ABX*I_ERI_Py_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Pz_Px_Fxyz_S_C1001001_cc = I_ERI_Dxz_S_Fxyz_S_C1001001_cc+ABX*I_ERI_Pz_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Px_Py_Fxyz_S_C1001001_cc = I_ERI_Dxy_S_Fxyz_S_C1001001_cc+ABY*I_ERI_Px_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Py_Py_Fxyz_S_C1001001_cc = I_ERI_D2y_S_Fxyz_S_C1001001_cc+ABY*I_ERI_Py_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Pz_Py_Fxyz_S_C1001001_cc = I_ERI_Dyz_S_Fxyz_S_C1001001_cc+ABY*I_ERI_Pz_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Px_Pz_Fxyz_S_C1001001_cc = I_ERI_Dxz_S_Fxyz_S_C1001001_cc+ABZ*I_ERI_Px_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Py_Pz_Fxyz_S_C1001001_cc = I_ERI_Dyz_S_Fxyz_S_C1001001_cc+ABZ*I_ERI_Py_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Pz_Pz_Fxyz_S_C1001001_cc = I_ERI_D2z_S_Fxyz_S_C1001001_cc+ABZ*I_ERI_Pz_S_Fxyz_S_C1001001_cc;
  Double I_ERI_Px_Px_Fx2z_S_C1001001_cc = I_ERI_D2x_S_Fx2z_S_C1001001_cc+ABX*I_ERI_Px_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Py_Px_Fx2z_S_C1001001_cc = I_ERI_Dxy_S_Fx2z_S_C1001001_cc+ABX*I_ERI_Py_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Pz_Px_Fx2z_S_C1001001_cc = I_ERI_Dxz_S_Fx2z_S_C1001001_cc+ABX*I_ERI_Pz_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Px_Py_Fx2z_S_C1001001_cc = I_ERI_Dxy_S_Fx2z_S_C1001001_cc+ABY*I_ERI_Px_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Py_Py_Fx2z_S_C1001001_cc = I_ERI_D2y_S_Fx2z_S_C1001001_cc+ABY*I_ERI_Py_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Pz_Py_Fx2z_S_C1001001_cc = I_ERI_Dyz_S_Fx2z_S_C1001001_cc+ABY*I_ERI_Pz_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Px_Pz_Fx2z_S_C1001001_cc = I_ERI_Dxz_S_Fx2z_S_C1001001_cc+ABZ*I_ERI_Px_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Py_Pz_Fx2z_S_C1001001_cc = I_ERI_Dyz_S_Fx2z_S_C1001001_cc+ABZ*I_ERI_Py_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Pz_Pz_Fx2z_S_C1001001_cc = I_ERI_D2z_S_Fx2z_S_C1001001_cc+ABZ*I_ERI_Pz_S_Fx2z_S_C1001001_cc;
  Double I_ERI_Px_Px_F3y_S_C1001001_cc = I_ERI_D2x_S_F3y_S_C1001001_cc+ABX*I_ERI_Px_S_F3y_S_C1001001_cc;
  Double I_ERI_Py_Px_F3y_S_C1001001_cc = I_ERI_Dxy_S_F3y_S_C1001001_cc+ABX*I_ERI_Py_S_F3y_S_C1001001_cc;
  Double I_ERI_Pz_Px_F3y_S_C1001001_cc = I_ERI_Dxz_S_F3y_S_C1001001_cc+ABX*I_ERI_Pz_S_F3y_S_C1001001_cc;
  Double I_ERI_Px_Py_F3y_S_C1001001_cc = I_ERI_Dxy_S_F3y_S_C1001001_cc+ABY*I_ERI_Px_S_F3y_S_C1001001_cc;
  Double I_ERI_Py_Py_F3y_S_C1001001_cc = I_ERI_D2y_S_F3y_S_C1001001_cc+ABY*I_ERI_Py_S_F3y_S_C1001001_cc;
  Double I_ERI_Pz_Py_F3y_S_C1001001_cc = I_ERI_Dyz_S_F3y_S_C1001001_cc+ABY*I_ERI_Pz_S_F3y_S_C1001001_cc;
  Double I_ERI_Px_Pz_F3y_S_C1001001_cc = I_ERI_Dxz_S_F3y_S_C1001001_cc+ABZ*I_ERI_Px_S_F3y_S_C1001001_cc;
  Double I_ERI_Py_Pz_F3y_S_C1001001_cc = I_ERI_Dyz_S_F3y_S_C1001001_cc+ABZ*I_ERI_Py_S_F3y_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F3y_S_C1001001_cc = I_ERI_D2z_S_F3y_S_C1001001_cc+ABZ*I_ERI_Pz_S_F3y_S_C1001001_cc;
  Double I_ERI_Px_Px_F2yz_S_C1001001_cc = I_ERI_D2x_S_F2yz_S_C1001001_cc+ABX*I_ERI_Px_S_F2yz_S_C1001001_cc;
  Double I_ERI_Py_Px_F2yz_S_C1001001_cc = I_ERI_Dxy_S_F2yz_S_C1001001_cc+ABX*I_ERI_Py_S_F2yz_S_C1001001_cc;
  Double I_ERI_Pz_Px_F2yz_S_C1001001_cc = I_ERI_Dxz_S_F2yz_S_C1001001_cc+ABX*I_ERI_Pz_S_F2yz_S_C1001001_cc;
  Double I_ERI_Px_Py_F2yz_S_C1001001_cc = I_ERI_Dxy_S_F2yz_S_C1001001_cc+ABY*I_ERI_Px_S_F2yz_S_C1001001_cc;
  Double I_ERI_Py_Py_F2yz_S_C1001001_cc = I_ERI_D2y_S_F2yz_S_C1001001_cc+ABY*I_ERI_Py_S_F2yz_S_C1001001_cc;
  Double I_ERI_Pz_Py_F2yz_S_C1001001_cc = I_ERI_Dyz_S_F2yz_S_C1001001_cc+ABY*I_ERI_Pz_S_F2yz_S_C1001001_cc;
  Double I_ERI_Px_Pz_F2yz_S_C1001001_cc = I_ERI_Dxz_S_F2yz_S_C1001001_cc+ABZ*I_ERI_Px_S_F2yz_S_C1001001_cc;
  Double I_ERI_Py_Pz_F2yz_S_C1001001_cc = I_ERI_Dyz_S_F2yz_S_C1001001_cc+ABZ*I_ERI_Py_S_F2yz_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F2yz_S_C1001001_cc = I_ERI_D2z_S_F2yz_S_C1001001_cc+ABZ*I_ERI_Pz_S_F2yz_S_C1001001_cc;
  Double I_ERI_Px_Px_Fy2z_S_C1001001_cc = I_ERI_D2x_S_Fy2z_S_C1001001_cc+ABX*I_ERI_Px_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Py_Px_Fy2z_S_C1001001_cc = I_ERI_Dxy_S_Fy2z_S_C1001001_cc+ABX*I_ERI_Py_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Pz_Px_Fy2z_S_C1001001_cc = I_ERI_Dxz_S_Fy2z_S_C1001001_cc+ABX*I_ERI_Pz_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Px_Py_Fy2z_S_C1001001_cc = I_ERI_Dxy_S_Fy2z_S_C1001001_cc+ABY*I_ERI_Px_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Py_Py_Fy2z_S_C1001001_cc = I_ERI_D2y_S_Fy2z_S_C1001001_cc+ABY*I_ERI_Py_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Pz_Py_Fy2z_S_C1001001_cc = I_ERI_Dyz_S_Fy2z_S_C1001001_cc+ABY*I_ERI_Pz_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Px_Pz_Fy2z_S_C1001001_cc = I_ERI_Dxz_S_Fy2z_S_C1001001_cc+ABZ*I_ERI_Px_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Py_Pz_Fy2z_S_C1001001_cc = I_ERI_Dyz_S_Fy2z_S_C1001001_cc+ABZ*I_ERI_Py_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Pz_Pz_Fy2z_S_C1001001_cc = I_ERI_D2z_S_Fy2z_S_C1001001_cc+ABZ*I_ERI_Pz_S_Fy2z_S_C1001001_cc;
  Double I_ERI_Px_Px_F3z_S_C1001001_cc = I_ERI_D2x_S_F3z_S_C1001001_cc+ABX*I_ERI_Px_S_F3z_S_C1001001_cc;
  Double I_ERI_Py_Px_F3z_S_C1001001_cc = I_ERI_Dxy_S_F3z_S_C1001001_cc+ABX*I_ERI_Py_S_F3z_S_C1001001_cc;
  Double I_ERI_Pz_Px_F3z_S_C1001001_cc = I_ERI_Dxz_S_F3z_S_C1001001_cc+ABX*I_ERI_Pz_S_F3z_S_C1001001_cc;
  Double I_ERI_Px_Py_F3z_S_C1001001_cc = I_ERI_Dxy_S_F3z_S_C1001001_cc+ABY*I_ERI_Px_S_F3z_S_C1001001_cc;
  Double I_ERI_Py_Py_F3z_S_C1001001_cc = I_ERI_D2y_S_F3z_S_C1001001_cc+ABY*I_ERI_Py_S_F3z_S_C1001001_cc;
  Double I_ERI_Pz_Py_F3z_S_C1001001_cc = I_ERI_Dyz_S_F3z_S_C1001001_cc+ABY*I_ERI_Pz_S_F3z_S_C1001001_cc;
  Double I_ERI_Px_Pz_F3z_S_C1001001_cc = I_ERI_Dxz_S_F3z_S_C1001001_cc+ABZ*I_ERI_Px_S_F3z_S_C1001001_cc;
  Double I_ERI_Py_Pz_F3z_S_C1001001_cc = I_ERI_Dyz_S_F3z_S_C1001001_cc+ABZ*I_ERI_Py_S_F3z_S_C1001001_cc;
  Double I_ERI_Pz_Pz_F3z_S_C1001001_cc = I_ERI_D2z_S_F3z_S_C1001001_cc+ABZ*I_ERI_Pz_S_F3z_S_C1001001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001001_cd
   ************************************************************/
  Double I_ERI_Px_Px_D2x_Px_C1001001_cd = I_ERI_D2x_S_D2x_Px_C1001001_cd+ABX*I_ERI_Px_S_D2x_Px_C1001001_cd;
  Double I_ERI_Py_Px_D2x_Px_C1001001_cd = I_ERI_Dxy_S_D2x_Px_C1001001_cd+ABX*I_ERI_Py_S_D2x_Px_C1001001_cd;
  Double I_ERI_Pz_Px_D2x_Px_C1001001_cd = I_ERI_Dxz_S_D2x_Px_C1001001_cd+ABX*I_ERI_Pz_S_D2x_Px_C1001001_cd;
  Double I_ERI_Px_Py_D2x_Px_C1001001_cd = I_ERI_Dxy_S_D2x_Px_C1001001_cd+ABY*I_ERI_Px_S_D2x_Px_C1001001_cd;
  Double I_ERI_Py_Py_D2x_Px_C1001001_cd = I_ERI_D2y_S_D2x_Px_C1001001_cd+ABY*I_ERI_Py_S_D2x_Px_C1001001_cd;
  Double I_ERI_Pz_Py_D2x_Px_C1001001_cd = I_ERI_Dyz_S_D2x_Px_C1001001_cd+ABY*I_ERI_Pz_S_D2x_Px_C1001001_cd;
  Double I_ERI_Px_Pz_D2x_Px_C1001001_cd = I_ERI_Dxz_S_D2x_Px_C1001001_cd+ABZ*I_ERI_Px_S_D2x_Px_C1001001_cd;
  Double I_ERI_Py_Pz_D2x_Px_C1001001_cd = I_ERI_Dyz_S_D2x_Px_C1001001_cd+ABZ*I_ERI_Py_S_D2x_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_D2x_Px_C1001001_cd = I_ERI_D2z_S_D2x_Px_C1001001_cd+ABZ*I_ERI_Pz_S_D2x_Px_C1001001_cd;
  Double I_ERI_Px_Px_Dxy_Px_C1001001_cd = I_ERI_D2x_S_Dxy_Px_C1001001_cd+ABX*I_ERI_Px_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Py_Px_Dxy_Px_C1001001_cd = I_ERI_Dxy_S_Dxy_Px_C1001001_cd+ABX*I_ERI_Py_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Pz_Px_Dxy_Px_C1001001_cd = I_ERI_Dxz_S_Dxy_Px_C1001001_cd+ABX*I_ERI_Pz_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Px_Py_Dxy_Px_C1001001_cd = I_ERI_Dxy_S_Dxy_Px_C1001001_cd+ABY*I_ERI_Px_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Py_Py_Dxy_Px_C1001001_cd = I_ERI_D2y_S_Dxy_Px_C1001001_cd+ABY*I_ERI_Py_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Pz_Py_Dxy_Px_C1001001_cd = I_ERI_Dyz_S_Dxy_Px_C1001001_cd+ABY*I_ERI_Pz_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Px_Pz_Dxy_Px_C1001001_cd = I_ERI_Dxz_S_Dxy_Px_C1001001_cd+ABZ*I_ERI_Px_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Py_Pz_Dxy_Px_C1001001_cd = I_ERI_Dyz_S_Dxy_Px_C1001001_cd+ABZ*I_ERI_Py_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxy_Px_C1001001_cd = I_ERI_D2z_S_Dxy_Px_C1001001_cd+ABZ*I_ERI_Pz_S_Dxy_Px_C1001001_cd;
  Double I_ERI_Px_Px_Dxz_Px_C1001001_cd = I_ERI_D2x_S_Dxz_Px_C1001001_cd+ABX*I_ERI_Px_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Py_Px_Dxz_Px_C1001001_cd = I_ERI_Dxy_S_Dxz_Px_C1001001_cd+ABX*I_ERI_Py_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Pz_Px_Dxz_Px_C1001001_cd = I_ERI_Dxz_S_Dxz_Px_C1001001_cd+ABX*I_ERI_Pz_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Px_Py_Dxz_Px_C1001001_cd = I_ERI_Dxy_S_Dxz_Px_C1001001_cd+ABY*I_ERI_Px_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Py_Py_Dxz_Px_C1001001_cd = I_ERI_D2y_S_Dxz_Px_C1001001_cd+ABY*I_ERI_Py_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Pz_Py_Dxz_Px_C1001001_cd = I_ERI_Dyz_S_Dxz_Px_C1001001_cd+ABY*I_ERI_Pz_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Px_Pz_Dxz_Px_C1001001_cd = I_ERI_Dxz_S_Dxz_Px_C1001001_cd+ABZ*I_ERI_Px_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Py_Pz_Dxz_Px_C1001001_cd = I_ERI_Dyz_S_Dxz_Px_C1001001_cd+ABZ*I_ERI_Py_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxz_Px_C1001001_cd = I_ERI_D2z_S_Dxz_Px_C1001001_cd+ABZ*I_ERI_Pz_S_Dxz_Px_C1001001_cd;
  Double I_ERI_Px_Px_D2y_Px_C1001001_cd = I_ERI_D2x_S_D2y_Px_C1001001_cd+ABX*I_ERI_Px_S_D2y_Px_C1001001_cd;
  Double I_ERI_Py_Px_D2y_Px_C1001001_cd = I_ERI_Dxy_S_D2y_Px_C1001001_cd+ABX*I_ERI_Py_S_D2y_Px_C1001001_cd;
  Double I_ERI_Pz_Px_D2y_Px_C1001001_cd = I_ERI_Dxz_S_D2y_Px_C1001001_cd+ABX*I_ERI_Pz_S_D2y_Px_C1001001_cd;
  Double I_ERI_Px_Py_D2y_Px_C1001001_cd = I_ERI_Dxy_S_D2y_Px_C1001001_cd+ABY*I_ERI_Px_S_D2y_Px_C1001001_cd;
  Double I_ERI_Py_Py_D2y_Px_C1001001_cd = I_ERI_D2y_S_D2y_Px_C1001001_cd+ABY*I_ERI_Py_S_D2y_Px_C1001001_cd;
  Double I_ERI_Pz_Py_D2y_Px_C1001001_cd = I_ERI_Dyz_S_D2y_Px_C1001001_cd+ABY*I_ERI_Pz_S_D2y_Px_C1001001_cd;
  Double I_ERI_Px_Pz_D2y_Px_C1001001_cd = I_ERI_Dxz_S_D2y_Px_C1001001_cd+ABZ*I_ERI_Px_S_D2y_Px_C1001001_cd;
  Double I_ERI_Py_Pz_D2y_Px_C1001001_cd = I_ERI_Dyz_S_D2y_Px_C1001001_cd+ABZ*I_ERI_Py_S_D2y_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_D2y_Px_C1001001_cd = I_ERI_D2z_S_D2y_Px_C1001001_cd+ABZ*I_ERI_Pz_S_D2y_Px_C1001001_cd;
  Double I_ERI_Px_Px_Dyz_Px_C1001001_cd = I_ERI_D2x_S_Dyz_Px_C1001001_cd+ABX*I_ERI_Px_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Py_Px_Dyz_Px_C1001001_cd = I_ERI_Dxy_S_Dyz_Px_C1001001_cd+ABX*I_ERI_Py_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Pz_Px_Dyz_Px_C1001001_cd = I_ERI_Dxz_S_Dyz_Px_C1001001_cd+ABX*I_ERI_Pz_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Px_Py_Dyz_Px_C1001001_cd = I_ERI_Dxy_S_Dyz_Px_C1001001_cd+ABY*I_ERI_Px_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Py_Py_Dyz_Px_C1001001_cd = I_ERI_D2y_S_Dyz_Px_C1001001_cd+ABY*I_ERI_Py_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Pz_Py_Dyz_Px_C1001001_cd = I_ERI_Dyz_S_Dyz_Px_C1001001_cd+ABY*I_ERI_Pz_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Px_Pz_Dyz_Px_C1001001_cd = I_ERI_Dxz_S_Dyz_Px_C1001001_cd+ABZ*I_ERI_Px_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Py_Pz_Dyz_Px_C1001001_cd = I_ERI_Dyz_S_Dyz_Px_C1001001_cd+ABZ*I_ERI_Py_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_Dyz_Px_C1001001_cd = I_ERI_D2z_S_Dyz_Px_C1001001_cd+ABZ*I_ERI_Pz_S_Dyz_Px_C1001001_cd;
  Double I_ERI_Px_Px_D2z_Px_C1001001_cd = I_ERI_D2x_S_D2z_Px_C1001001_cd+ABX*I_ERI_Px_S_D2z_Px_C1001001_cd;
  Double I_ERI_Py_Px_D2z_Px_C1001001_cd = I_ERI_Dxy_S_D2z_Px_C1001001_cd+ABX*I_ERI_Py_S_D2z_Px_C1001001_cd;
  Double I_ERI_Pz_Px_D2z_Px_C1001001_cd = I_ERI_Dxz_S_D2z_Px_C1001001_cd+ABX*I_ERI_Pz_S_D2z_Px_C1001001_cd;
  Double I_ERI_Px_Py_D2z_Px_C1001001_cd = I_ERI_Dxy_S_D2z_Px_C1001001_cd+ABY*I_ERI_Px_S_D2z_Px_C1001001_cd;
  Double I_ERI_Py_Py_D2z_Px_C1001001_cd = I_ERI_D2y_S_D2z_Px_C1001001_cd+ABY*I_ERI_Py_S_D2z_Px_C1001001_cd;
  Double I_ERI_Pz_Py_D2z_Px_C1001001_cd = I_ERI_Dyz_S_D2z_Px_C1001001_cd+ABY*I_ERI_Pz_S_D2z_Px_C1001001_cd;
  Double I_ERI_Px_Pz_D2z_Px_C1001001_cd = I_ERI_Dxz_S_D2z_Px_C1001001_cd+ABZ*I_ERI_Px_S_D2z_Px_C1001001_cd;
  Double I_ERI_Py_Pz_D2z_Px_C1001001_cd = I_ERI_Dyz_S_D2z_Px_C1001001_cd+ABZ*I_ERI_Py_S_D2z_Px_C1001001_cd;
  Double I_ERI_Pz_Pz_D2z_Px_C1001001_cd = I_ERI_D2z_S_D2z_Px_C1001001_cd+ABZ*I_ERI_Pz_S_D2z_Px_C1001001_cd;
  Double I_ERI_Px_Px_D2x_Py_C1001001_cd = I_ERI_D2x_S_D2x_Py_C1001001_cd+ABX*I_ERI_Px_S_D2x_Py_C1001001_cd;
  Double I_ERI_Py_Px_D2x_Py_C1001001_cd = I_ERI_Dxy_S_D2x_Py_C1001001_cd+ABX*I_ERI_Py_S_D2x_Py_C1001001_cd;
  Double I_ERI_Pz_Px_D2x_Py_C1001001_cd = I_ERI_Dxz_S_D2x_Py_C1001001_cd+ABX*I_ERI_Pz_S_D2x_Py_C1001001_cd;
  Double I_ERI_Px_Py_D2x_Py_C1001001_cd = I_ERI_Dxy_S_D2x_Py_C1001001_cd+ABY*I_ERI_Px_S_D2x_Py_C1001001_cd;
  Double I_ERI_Py_Py_D2x_Py_C1001001_cd = I_ERI_D2y_S_D2x_Py_C1001001_cd+ABY*I_ERI_Py_S_D2x_Py_C1001001_cd;
  Double I_ERI_Pz_Py_D2x_Py_C1001001_cd = I_ERI_Dyz_S_D2x_Py_C1001001_cd+ABY*I_ERI_Pz_S_D2x_Py_C1001001_cd;
  Double I_ERI_Px_Pz_D2x_Py_C1001001_cd = I_ERI_Dxz_S_D2x_Py_C1001001_cd+ABZ*I_ERI_Px_S_D2x_Py_C1001001_cd;
  Double I_ERI_Py_Pz_D2x_Py_C1001001_cd = I_ERI_Dyz_S_D2x_Py_C1001001_cd+ABZ*I_ERI_Py_S_D2x_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_D2x_Py_C1001001_cd = I_ERI_D2z_S_D2x_Py_C1001001_cd+ABZ*I_ERI_Pz_S_D2x_Py_C1001001_cd;
  Double I_ERI_Px_Px_Dxy_Py_C1001001_cd = I_ERI_D2x_S_Dxy_Py_C1001001_cd+ABX*I_ERI_Px_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Py_Px_Dxy_Py_C1001001_cd = I_ERI_Dxy_S_Dxy_Py_C1001001_cd+ABX*I_ERI_Py_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Pz_Px_Dxy_Py_C1001001_cd = I_ERI_Dxz_S_Dxy_Py_C1001001_cd+ABX*I_ERI_Pz_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Px_Py_Dxy_Py_C1001001_cd = I_ERI_Dxy_S_Dxy_Py_C1001001_cd+ABY*I_ERI_Px_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Py_Py_Dxy_Py_C1001001_cd = I_ERI_D2y_S_Dxy_Py_C1001001_cd+ABY*I_ERI_Py_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Pz_Py_Dxy_Py_C1001001_cd = I_ERI_Dyz_S_Dxy_Py_C1001001_cd+ABY*I_ERI_Pz_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Px_Pz_Dxy_Py_C1001001_cd = I_ERI_Dxz_S_Dxy_Py_C1001001_cd+ABZ*I_ERI_Px_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Py_Pz_Dxy_Py_C1001001_cd = I_ERI_Dyz_S_Dxy_Py_C1001001_cd+ABZ*I_ERI_Py_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxy_Py_C1001001_cd = I_ERI_D2z_S_Dxy_Py_C1001001_cd+ABZ*I_ERI_Pz_S_Dxy_Py_C1001001_cd;
  Double I_ERI_Px_Px_Dxz_Py_C1001001_cd = I_ERI_D2x_S_Dxz_Py_C1001001_cd+ABX*I_ERI_Px_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Py_Px_Dxz_Py_C1001001_cd = I_ERI_Dxy_S_Dxz_Py_C1001001_cd+ABX*I_ERI_Py_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Pz_Px_Dxz_Py_C1001001_cd = I_ERI_Dxz_S_Dxz_Py_C1001001_cd+ABX*I_ERI_Pz_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Px_Py_Dxz_Py_C1001001_cd = I_ERI_Dxy_S_Dxz_Py_C1001001_cd+ABY*I_ERI_Px_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Py_Py_Dxz_Py_C1001001_cd = I_ERI_D2y_S_Dxz_Py_C1001001_cd+ABY*I_ERI_Py_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Pz_Py_Dxz_Py_C1001001_cd = I_ERI_Dyz_S_Dxz_Py_C1001001_cd+ABY*I_ERI_Pz_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Px_Pz_Dxz_Py_C1001001_cd = I_ERI_Dxz_S_Dxz_Py_C1001001_cd+ABZ*I_ERI_Px_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Py_Pz_Dxz_Py_C1001001_cd = I_ERI_Dyz_S_Dxz_Py_C1001001_cd+ABZ*I_ERI_Py_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxz_Py_C1001001_cd = I_ERI_D2z_S_Dxz_Py_C1001001_cd+ABZ*I_ERI_Pz_S_Dxz_Py_C1001001_cd;
  Double I_ERI_Px_Px_D2y_Py_C1001001_cd = I_ERI_D2x_S_D2y_Py_C1001001_cd+ABX*I_ERI_Px_S_D2y_Py_C1001001_cd;
  Double I_ERI_Py_Px_D2y_Py_C1001001_cd = I_ERI_Dxy_S_D2y_Py_C1001001_cd+ABX*I_ERI_Py_S_D2y_Py_C1001001_cd;
  Double I_ERI_Pz_Px_D2y_Py_C1001001_cd = I_ERI_Dxz_S_D2y_Py_C1001001_cd+ABX*I_ERI_Pz_S_D2y_Py_C1001001_cd;
  Double I_ERI_Px_Py_D2y_Py_C1001001_cd = I_ERI_Dxy_S_D2y_Py_C1001001_cd+ABY*I_ERI_Px_S_D2y_Py_C1001001_cd;
  Double I_ERI_Py_Py_D2y_Py_C1001001_cd = I_ERI_D2y_S_D2y_Py_C1001001_cd+ABY*I_ERI_Py_S_D2y_Py_C1001001_cd;
  Double I_ERI_Pz_Py_D2y_Py_C1001001_cd = I_ERI_Dyz_S_D2y_Py_C1001001_cd+ABY*I_ERI_Pz_S_D2y_Py_C1001001_cd;
  Double I_ERI_Px_Pz_D2y_Py_C1001001_cd = I_ERI_Dxz_S_D2y_Py_C1001001_cd+ABZ*I_ERI_Px_S_D2y_Py_C1001001_cd;
  Double I_ERI_Py_Pz_D2y_Py_C1001001_cd = I_ERI_Dyz_S_D2y_Py_C1001001_cd+ABZ*I_ERI_Py_S_D2y_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_D2y_Py_C1001001_cd = I_ERI_D2z_S_D2y_Py_C1001001_cd+ABZ*I_ERI_Pz_S_D2y_Py_C1001001_cd;
  Double I_ERI_Px_Px_Dyz_Py_C1001001_cd = I_ERI_D2x_S_Dyz_Py_C1001001_cd+ABX*I_ERI_Px_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Py_Px_Dyz_Py_C1001001_cd = I_ERI_Dxy_S_Dyz_Py_C1001001_cd+ABX*I_ERI_Py_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Pz_Px_Dyz_Py_C1001001_cd = I_ERI_Dxz_S_Dyz_Py_C1001001_cd+ABX*I_ERI_Pz_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Px_Py_Dyz_Py_C1001001_cd = I_ERI_Dxy_S_Dyz_Py_C1001001_cd+ABY*I_ERI_Px_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Py_Py_Dyz_Py_C1001001_cd = I_ERI_D2y_S_Dyz_Py_C1001001_cd+ABY*I_ERI_Py_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Pz_Py_Dyz_Py_C1001001_cd = I_ERI_Dyz_S_Dyz_Py_C1001001_cd+ABY*I_ERI_Pz_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Px_Pz_Dyz_Py_C1001001_cd = I_ERI_Dxz_S_Dyz_Py_C1001001_cd+ABZ*I_ERI_Px_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Py_Pz_Dyz_Py_C1001001_cd = I_ERI_Dyz_S_Dyz_Py_C1001001_cd+ABZ*I_ERI_Py_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_Dyz_Py_C1001001_cd = I_ERI_D2z_S_Dyz_Py_C1001001_cd+ABZ*I_ERI_Pz_S_Dyz_Py_C1001001_cd;
  Double I_ERI_Px_Px_D2z_Py_C1001001_cd = I_ERI_D2x_S_D2z_Py_C1001001_cd+ABX*I_ERI_Px_S_D2z_Py_C1001001_cd;
  Double I_ERI_Py_Px_D2z_Py_C1001001_cd = I_ERI_Dxy_S_D2z_Py_C1001001_cd+ABX*I_ERI_Py_S_D2z_Py_C1001001_cd;
  Double I_ERI_Pz_Px_D2z_Py_C1001001_cd = I_ERI_Dxz_S_D2z_Py_C1001001_cd+ABX*I_ERI_Pz_S_D2z_Py_C1001001_cd;
  Double I_ERI_Px_Py_D2z_Py_C1001001_cd = I_ERI_Dxy_S_D2z_Py_C1001001_cd+ABY*I_ERI_Px_S_D2z_Py_C1001001_cd;
  Double I_ERI_Py_Py_D2z_Py_C1001001_cd = I_ERI_D2y_S_D2z_Py_C1001001_cd+ABY*I_ERI_Py_S_D2z_Py_C1001001_cd;
  Double I_ERI_Pz_Py_D2z_Py_C1001001_cd = I_ERI_Dyz_S_D2z_Py_C1001001_cd+ABY*I_ERI_Pz_S_D2z_Py_C1001001_cd;
  Double I_ERI_Px_Pz_D2z_Py_C1001001_cd = I_ERI_Dxz_S_D2z_Py_C1001001_cd+ABZ*I_ERI_Px_S_D2z_Py_C1001001_cd;
  Double I_ERI_Py_Pz_D2z_Py_C1001001_cd = I_ERI_Dyz_S_D2z_Py_C1001001_cd+ABZ*I_ERI_Py_S_D2z_Py_C1001001_cd;
  Double I_ERI_Pz_Pz_D2z_Py_C1001001_cd = I_ERI_D2z_S_D2z_Py_C1001001_cd+ABZ*I_ERI_Pz_S_D2z_Py_C1001001_cd;
  Double I_ERI_Px_Px_D2x_Pz_C1001001_cd = I_ERI_D2x_S_D2x_Pz_C1001001_cd+ABX*I_ERI_Px_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Py_Px_D2x_Pz_C1001001_cd = I_ERI_Dxy_S_D2x_Pz_C1001001_cd+ABX*I_ERI_Py_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_D2x_Pz_C1001001_cd = I_ERI_Dxz_S_D2x_Pz_C1001001_cd+ABX*I_ERI_Pz_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Px_Py_D2x_Pz_C1001001_cd = I_ERI_Dxy_S_D2x_Pz_C1001001_cd+ABY*I_ERI_Px_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Py_Py_D2x_Pz_C1001001_cd = I_ERI_D2y_S_D2x_Pz_C1001001_cd+ABY*I_ERI_Py_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_D2x_Pz_C1001001_cd = I_ERI_Dyz_S_D2x_Pz_C1001001_cd+ABY*I_ERI_Pz_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_D2x_Pz_C1001001_cd = I_ERI_Dxz_S_D2x_Pz_C1001001_cd+ABZ*I_ERI_Px_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_D2x_Pz_C1001001_cd = I_ERI_Dyz_S_D2x_Pz_C1001001_cd+ABZ*I_ERI_Py_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_D2x_Pz_C1001001_cd = I_ERI_D2z_S_D2x_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_D2x_Pz_C1001001_cd;
  Double I_ERI_Px_Px_Dxy_Pz_C1001001_cd = I_ERI_D2x_S_Dxy_Pz_C1001001_cd+ABX*I_ERI_Px_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Py_Px_Dxy_Pz_C1001001_cd = I_ERI_Dxy_S_Dxy_Pz_C1001001_cd+ABX*I_ERI_Py_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_Dxy_Pz_C1001001_cd = I_ERI_Dxz_S_Dxy_Pz_C1001001_cd+ABX*I_ERI_Pz_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Px_Py_Dxy_Pz_C1001001_cd = I_ERI_Dxy_S_Dxy_Pz_C1001001_cd+ABY*I_ERI_Px_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Py_Py_Dxy_Pz_C1001001_cd = I_ERI_D2y_S_Dxy_Pz_C1001001_cd+ABY*I_ERI_Py_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_Dxy_Pz_C1001001_cd = I_ERI_Dyz_S_Dxy_Pz_C1001001_cd+ABY*I_ERI_Pz_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_Dxy_Pz_C1001001_cd = I_ERI_Dxz_S_Dxy_Pz_C1001001_cd+ABZ*I_ERI_Px_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_Dxy_Pz_C1001001_cd = I_ERI_Dyz_S_Dxy_Pz_C1001001_cd+ABZ*I_ERI_Py_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxy_Pz_C1001001_cd = I_ERI_D2z_S_Dxy_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_Dxy_Pz_C1001001_cd;
  Double I_ERI_Px_Px_Dxz_Pz_C1001001_cd = I_ERI_D2x_S_Dxz_Pz_C1001001_cd+ABX*I_ERI_Px_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Py_Px_Dxz_Pz_C1001001_cd = I_ERI_Dxy_S_Dxz_Pz_C1001001_cd+ABX*I_ERI_Py_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_Dxz_Pz_C1001001_cd = I_ERI_Dxz_S_Dxz_Pz_C1001001_cd+ABX*I_ERI_Pz_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Px_Py_Dxz_Pz_C1001001_cd = I_ERI_Dxy_S_Dxz_Pz_C1001001_cd+ABY*I_ERI_Px_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Py_Py_Dxz_Pz_C1001001_cd = I_ERI_D2y_S_Dxz_Pz_C1001001_cd+ABY*I_ERI_Py_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_Dxz_Pz_C1001001_cd = I_ERI_Dyz_S_Dxz_Pz_C1001001_cd+ABY*I_ERI_Pz_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_Dxz_Pz_C1001001_cd = I_ERI_Dxz_S_Dxz_Pz_C1001001_cd+ABZ*I_ERI_Px_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_Dxz_Pz_C1001001_cd = I_ERI_Dyz_S_Dxz_Pz_C1001001_cd+ABZ*I_ERI_Py_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_Dxz_Pz_C1001001_cd = I_ERI_D2z_S_Dxz_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_Dxz_Pz_C1001001_cd;
  Double I_ERI_Px_Px_D2y_Pz_C1001001_cd = I_ERI_D2x_S_D2y_Pz_C1001001_cd+ABX*I_ERI_Px_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Py_Px_D2y_Pz_C1001001_cd = I_ERI_Dxy_S_D2y_Pz_C1001001_cd+ABX*I_ERI_Py_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_D2y_Pz_C1001001_cd = I_ERI_Dxz_S_D2y_Pz_C1001001_cd+ABX*I_ERI_Pz_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Px_Py_D2y_Pz_C1001001_cd = I_ERI_Dxy_S_D2y_Pz_C1001001_cd+ABY*I_ERI_Px_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Py_Py_D2y_Pz_C1001001_cd = I_ERI_D2y_S_D2y_Pz_C1001001_cd+ABY*I_ERI_Py_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_D2y_Pz_C1001001_cd = I_ERI_Dyz_S_D2y_Pz_C1001001_cd+ABY*I_ERI_Pz_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_D2y_Pz_C1001001_cd = I_ERI_Dxz_S_D2y_Pz_C1001001_cd+ABZ*I_ERI_Px_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_D2y_Pz_C1001001_cd = I_ERI_Dyz_S_D2y_Pz_C1001001_cd+ABZ*I_ERI_Py_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_D2y_Pz_C1001001_cd = I_ERI_D2z_S_D2y_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_D2y_Pz_C1001001_cd;
  Double I_ERI_Px_Px_Dyz_Pz_C1001001_cd = I_ERI_D2x_S_Dyz_Pz_C1001001_cd+ABX*I_ERI_Px_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Py_Px_Dyz_Pz_C1001001_cd = I_ERI_Dxy_S_Dyz_Pz_C1001001_cd+ABX*I_ERI_Py_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_Dyz_Pz_C1001001_cd = I_ERI_Dxz_S_Dyz_Pz_C1001001_cd+ABX*I_ERI_Pz_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Px_Py_Dyz_Pz_C1001001_cd = I_ERI_Dxy_S_Dyz_Pz_C1001001_cd+ABY*I_ERI_Px_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Py_Py_Dyz_Pz_C1001001_cd = I_ERI_D2y_S_Dyz_Pz_C1001001_cd+ABY*I_ERI_Py_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_Dyz_Pz_C1001001_cd = I_ERI_Dyz_S_Dyz_Pz_C1001001_cd+ABY*I_ERI_Pz_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_Dyz_Pz_C1001001_cd = I_ERI_Dxz_S_Dyz_Pz_C1001001_cd+ABZ*I_ERI_Px_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_Dyz_Pz_C1001001_cd = I_ERI_Dyz_S_Dyz_Pz_C1001001_cd+ABZ*I_ERI_Py_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_Dyz_Pz_C1001001_cd = I_ERI_D2z_S_Dyz_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_Dyz_Pz_C1001001_cd;
  Double I_ERI_Px_Px_D2z_Pz_C1001001_cd = I_ERI_D2x_S_D2z_Pz_C1001001_cd+ABX*I_ERI_Px_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Py_Px_D2z_Pz_C1001001_cd = I_ERI_Dxy_S_D2z_Pz_C1001001_cd+ABX*I_ERI_Py_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Pz_Px_D2z_Pz_C1001001_cd = I_ERI_Dxz_S_D2z_Pz_C1001001_cd+ABX*I_ERI_Pz_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Px_Py_D2z_Pz_C1001001_cd = I_ERI_Dxy_S_D2z_Pz_C1001001_cd+ABY*I_ERI_Px_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Py_Py_D2z_Pz_C1001001_cd = I_ERI_D2y_S_D2z_Pz_C1001001_cd+ABY*I_ERI_Py_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Pz_Py_D2z_Pz_C1001001_cd = I_ERI_Dyz_S_D2z_Pz_C1001001_cd+ABY*I_ERI_Pz_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Px_Pz_D2z_Pz_C1001001_cd = I_ERI_Dxz_S_D2z_Pz_C1001001_cd+ABZ*I_ERI_Px_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Py_Pz_D2z_Pz_C1001001_cd = I_ERI_Dyz_S_D2z_Pz_C1001001_cd+ABZ*I_ERI_Py_S_D2z_Pz_C1001001_cd;
  Double I_ERI_Pz_Pz_D2z_Pz_C1001001_cd = I_ERI_D2z_S_D2z_Pz_C1001001_cd+ABZ*I_ERI_Pz_S_D2z_Pz_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1001001_dd
   ************************************************************/
  Double I_ERI_Px_Px_Px_D2x_C1001001_dd = I_ERI_D2x_S_Px_D2x_C1001001_dd+ABX*I_ERI_Px_S_Px_D2x_C1001001_dd;
  Double I_ERI_Py_Px_Px_D2x_C1001001_dd = I_ERI_Dxy_S_Px_D2x_C1001001_dd+ABX*I_ERI_Py_S_Px_D2x_C1001001_dd;
  Double I_ERI_Pz_Px_Px_D2x_C1001001_dd = I_ERI_Dxz_S_Px_D2x_C1001001_dd+ABX*I_ERI_Pz_S_Px_D2x_C1001001_dd;
  Double I_ERI_Px_Py_Px_D2x_C1001001_dd = I_ERI_Dxy_S_Px_D2x_C1001001_dd+ABY*I_ERI_Px_S_Px_D2x_C1001001_dd;
  Double I_ERI_Py_Py_Px_D2x_C1001001_dd = I_ERI_D2y_S_Px_D2x_C1001001_dd+ABY*I_ERI_Py_S_Px_D2x_C1001001_dd;
  Double I_ERI_Pz_Py_Px_D2x_C1001001_dd = I_ERI_Dyz_S_Px_D2x_C1001001_dd+ABY*I_ERI_Pz_S_Px_D2x_C1001001_dd;
  Double I_ERI_Px_Pz_Px_D2x_C1001001_dd = I_ERI_Dxz_S_Px_D2x_C1001001_dd+ABZ*I_ERI_Px_S_Px_D2x_C1001001_dd;
  Double I_ERI_Py_Pz_Px_D2x_C1001001_dd = I_ERI_Dyz_S_Px_D2x_C1001001_dd+ABZ*I_ERI_Py_S_Px_D2x_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_D2x_C1001001_dd = I_ERI_D2z_S_Px_D2x_C1001001_dd+ABZ*I_ERI_Pz_S_Px_D2x_C1001001_dd;
  Double I_ERI_Px_Px_Py_D2x_C1001001_dd = I_ERI_D2x_S_Py_D2x_C1001001_dd+ABX*I_ERI_Px_S_Py_D2x_C1001001_dd;
  Double I_ERI_Py_Px_Py_D2x_C1001001_dd = I_ERI_Dxy_S_Py_D2x_C1001001_dd+ABX*I_ERI_Py_S_Py_D2x_C1001001_dd;
  Double I_ERI_Pz_Px_Py_D2x_C1001001_dd = I_ERI_Dxz_S_Py_D2x_C1001001_dd+ABX*I_ERI_Pz_S_Py_D2x_C1001001_dd;
  Double I_ERI_Px_Py_Py_D2x_C1001001_dd = I_ERI_Dxy_S_Py_D2x_C1001001_dd+ABY*I_ERI_Px_S_Py_D2x_C1001001_dd;
  Double I_ERI_Py_Py_Py_D2x_C1001001_dd = I_ERI_D2y_S_Py_D2x_C1001001_dd+ABY*I_ERI_Py_S_Py_D2x_C1001001_dd;
  Double I_ERI_Pz_Py_Py_D2x_C1001001_dd = I_ERI_Dyz_S_Py_D2x_C1001001_dd+ABY*I_ERI_Pz_S_Py_D2x_C1001001_dd;
  Double I_ERI_Px_Pz_Py_D2x_C1001001_dd = I_ERI_Dxz_S_Py_D2x_C1001001_dd+ABZ*I_ERI_Px_S_Py_D2x_C1001001_dd;
  Double I_ERI_Py_Pz_Py_D2x_C1001001_dd = I_ERI_Dyz_S_Py_D2x_C1001001_dd+ABZ*I_ERI_Py_S_Py_D2x_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_D2x_C1001001_dd = I_ERI_D2z_S_Py_D2x_C1001001_dd+ABZ*I_ERI_Pz_S_Py_D2x_C1001001_dd;
  Double I_ERI_Px_Px_Pz_D2x_C1001001_dd = I_ERI_D2x_S_Pz_D2x_C1001001_dd+ABX*I_ERI_Px_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Py_Px_Pz_D2x_C1001001_dd = I_ERI_Dxy_S_Pz_D2x_C1001001_dd+ABX*I_ERI_Py_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_D2x_C1001001_dd = I_ERI_Dxz_S_Pz_D2x_C1001001_dd+ABX*I_ERI_Pz_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Px_Py_Pz_D2x_C1001001_dd = I_ERI_Dxy_S_Pz_D2x_C1001001_dd+ABY*I_ERI_Px_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Py_Py_Pz_D2x_C1001001_dd = I_ERI_D2y_S_Pz_D2x_C1001001_dd+ABY*I_ERI_Py_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_D2x_C1001001_dd = I_ERI_Dyz_S_Pz_D2x_C1001001_dd+ABY*I_ERI_Pz_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_D2x_C1001001_dd = I_ERI_Dxz_S_Pz_D2x_C1001001_dd+ABZ*I_ERI_Px_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_D2x_C1001001_dd = I_ERI_Dyz_S_Pz_D2x_C1001001_dd+ABZ*I_ERI_Py_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_D2x_C1001001_dd = I_ERI_D2z_S_Pz_D2x_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_D2x_C1001001_dd;
  Double I_ERI_Px_Px_Px_Dxy_C1001001_dd = I_ERI_D2x_S_Px_Dxy_C1001001_dd+ABX*I_ERI_Px_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Py_Px_Px_Dxy_C1001001_dd = I_ERI_Dxy_S_Px_Dxy_C1001001_dd+ABX*I_ERI_Py_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Pz_Px_Px_Dxy_C1001001_dd = I_ERI_Dxz_S_Px_Dxy_C1001001_dd+ABX*I_ERI_Pz_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Px_Py_Px_Dxy_C1001001_dd = I_ERI_Dxy_S_Px_Dxy_C1001001_dd+ABY*I_ERI_Px_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Py_Py_Px_Dxy_C1001001_dd = I_ERI_D2y_S_Px_Dxy_C1001001_dd+ABY*I_ERI_Py_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Pz_Py_Px_Dxy_C1001001_dd = I_ERI_Dyz_S_Px_Dxy_C1001001_dd+ABY*I_ERI_Pz_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Px_Pz_Px_Dxy_C1001001_dd = I_ERI_Dxz_S_Px_Dxy_C1001001_dd+ABZ*I_ERI_Px_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Py_Pz_Px_Dxy_C1001001_dd = I_ERI_Dyz_S_Px_Dxy_C1001001_dd+ABZ*I_ERI_Py_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_Dxy_C1001001_dd = I_ERI_D2z_S_Px_Dxy_C1001001_dd+ABZ*I_ERI_Pz_S_Px_Dxy_C1001001_dd;
  Double I_ERI_Px_Px_Py_Dxy_C1001001_dd = I_ERI_D2x_S_Py_Dxy_C1001001_dd+ABX*I_ERI_Px_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Py_Px_Py_Dxy_C1001001_dd = I_ERI_Dxy_S_Py_Dxy_C1001001_dd+ABX*I_ERI_Py_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Pz_Px_Py_Dxy_C1001001_dd = I_ERI_Dxz_S_Py_Dxy_C1001001_dd+ABX*I_ERI_Pz_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Px_Py_Py_Dxy_C1001001_dd = I_ERI_Dxy_S_Py_Dxy_C1001001_dd+ABY*I_ERI_Px_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Py_Py_Py_Dxy_C1001001_dd = I_ERI_D2y_S_Py_Dxy_C1001001_dd+ABY*I_ERI_Py_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Pz_Py_Py_Dxy_C1001001_dd = I_ERI_Dyz_S_Py_Dxy_C1001001_dd+ABY*I_ERI_Pz_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Px_Pz_Py_Dxy_C1001001_dd = I_ERI_Dxz_S_Py_Dxy_C1001001_dd+ABZ*I_ERI_Px_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Py_Pz_Py_Dxy_C1001001_dd = I_ERI_Dyz_S_Py_Dxy_C1001001_dd+ABZ*I_ERI_Py_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_Dxy_C1001001_dd = I_ERI_D2z_S_Py_Dxy_C1001001_dd+ABZ*I_ERI_Pz_S_Py_Dxy_C1001001_dd;
  Double I_ERI_Px_Px_Pz_Dxy_C1001001_dd = I_ERI_D2x_S_Pz_Dxy_C1001001_dd+ABX*I_ERI_Px_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Py_Px_Pz_Dxy_C1001001_dd = I_ERI_Dxy_S_Pz_Dxy_C1001001_dd+ABX*I_ERI_Py_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_Dxy_C1001001_dd = I_ERI_Dxz_S_Pz_Dxy_C1001001_dd+ABX*I_ERI_Pz_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Px_Py_Pz_Dxy_C1001001_dd = I_ERI_Dxy_S_Pz_Dxy_C1001001_dd+ABY*I_ERI_Px_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Py_Py_Pz_Dxy_C1001001_dd = I_ERI_D2y_S_Pz_Dxy_C1001001_dd+ABY*I_ERI_Py_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_Dxy_C1001001_dd = I_ERI_Dyz_S_Pz_Dxy_C1001001_dd+ABY*I_ERI_Pz_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_Dxy_C1001001_dd = I_ERI_Dxz_S_Pz_Dxy_C1001001_dd+ABZ*I_ERI_Px_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_Dxy_C1001001_dd = I_ERI_Dyz_S_Pz_Dxy_C1001001_dd+ABZ*I_ERI_Py_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_Dxy_C1001001_dd = I_ERI_D2z_S_Pz_Dxy_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_Dxy_C1001001_dd;
  Double I_ERI_Px_Px_Px_Dxz_C1001001_dd = I_ERI_D2x_S_Px_Dxz_C1001001_dd+ABX*I_ERI_Px_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Py_Px_Px_Dxz_C1001001_dd = I_ERI_Dxy_S_Px_Dxz_C1001001_dd+ABX*I_ERI_Py_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Pz_Px_Px_Dxz_C1001001_dd = I_ERI_Dxz_S_Px_Dxz_C1001001_dd+ABX*I_ERI_Pz_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Px_Py_Px_Dxz_C1001001_dd = I_ERI_Dxy_S_Px_Dxz_C1001001_dd+ABY*I_ERI_Px_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Py_Py_Px_Dxz_C1001001_dd = I_ERI_D2y_S_Px_Dxz_C1001001_dd+ABY*I_ERI_Py_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Pz_Py_Px_Dxz_C1001001_dd = I_ERI_Dyz_S_Px_Dxz_C1001001_dd+ABY*I_ERI_Pz_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Px_Pz_Px_Dxz_C1001001_dd = I_ERI_Dxz_S_Px_Dxz_C1001001_dd+ABZ*I_ERI_Px_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Py_Pz_Px_Dxz_C1001001_dd = I_ERI_Dyz_S_Px_Dxz_C1001001_dd+ABZ*I_ERI_Py_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_Dxz_C1001001_dd = I_ERI_D2z_S_Px_Dxz_C1001001_dd+ABZ*I_ERI_Pz_S_Px_Dxz_C1001001_dd;
  Double I_ERI_Px_Px_Py_Dxz_C1001001_dd = I_ERI_D2x_S_Py_Dxz_C1001001_dd+ABX*I_ERI_Px_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Py_Px_Py_Dxz_C1001001_dd = I_ERI_Dxy_S_Py_Dxz_C1001001_dd+ABX*I_ERI_Py_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Pz_Px_Py_Dxz_C1001001_dd = I_ERI_Dxz_S_Py_Dxz_C1001001_dd+ABX*I_ERI_Pz_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Px_Py_Py_Dxz_C1001001_dd = I_ERI_Dxy_S_Py_Dxz_C1001001_dd+ABY*I_ERI_Px_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Py_Py_Py_Dxz_C1001001_dd = I_ERI_D2y_S_Py_Dxz_C1001001_dd+ABY*I_ERI_Py_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Pz_Py_Py_Dxz_C1001001_dd = I_ERI_Dyz_S_Py_Dxz_C1001001_dd+ABY*I_ERI_Pz_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Px_Pz_Py_Dxz_C1001001_dd = I_ERI_Dxz_S_Py_Dxz_C1001001_dd+ABZ*I_ERI_Px_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Py_Pz_Py_Dxz_C1001001_dd = I_ERI_Dyz_S_Py_Dxz_C1001001_dd+ABZ*I_ERI_Py_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_Dxz_C1001001_dd = I_ERI_D2z_S_Py_Dxz_C1001001_dd+ABZ*I_ERI_Pz_S_Py_Dxz_C1001001_dd;
  Double I_ERI_Px_Px_Pz_Dxz_C1001001_dd = I_ERI_D2x_S_Pz_Dxz_C1001001_dd+ABX*I_ERI_Px_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Py_Px_Pz_Dxz_C1001001_dd = I_ERI_Dxy_S_Pz_Dxz_C1001001_dd+ABX*I_ERI_Py_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_Dxz_C1001001_dd = I_ERI_Dxz_S_Pz_Dxz_C1001001_dd+ABX*I_ERI_Pz_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Px_Py_Pz_Dxz_C1001001_dd = I_ERI_Dxy_S_Pz_Dxz_C1001001_dd+ABY*I_ERI_Px_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Py_Py_Pz_Dxz_C1001001_dd = I_ERI_D2y_S_Pz_Dxz_C1001001_dd+ABY*I_ERI_Py_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_Dxz_C1001001_dd = I_ERI_Dyz_S_Pz_Dxz_C1001001_dd+ABY*I_ERI_Pz_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_Dxz_C1001001_dd = I_ERI_Dxz_S_Pz_Dxz_C1001001_dd+ABZ*I_ERI_Px_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_Dxz_C1001001_dd = I_ERI_Dyz_S_Pz_Dxz_C1001001_dd+ABZ*I_ERI_Py_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_Dxz_C1001001_dd = I_ERI_D2z_S_Pz_Dxz_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_Dxz_C1001001_dd;
  Double I_ERI_Px_Px_Px_D2y_C1001001_dd = I_ERI_D2x_S_Px_D2y_C1001001_dd+ABX*I_ERI_Px_S_Px_D2y_C1001001_dd;
  Double I_ERI_Py_Px_Px_D2y_C1001001_dd = I_ERI_Dxy_S_Px_D2y_C1001001_dd+ABX*I_ERI_Py_S_Px_D2y_C1001001_dd;
  Double I_ERI_Pz_Px_Px_D2y_C1001001_dd = I_ERI_Dxz_S_Px_D2y_C1001001_dd+ABX*I_ERI_Pz_S_Px_D2y_C1001001_dd;
  Double I_ERI_Px_Py_Px_D2y_C1001001_dd = I_ERI_Dxy_S_Px_D2y_C1001001_dd+ABY*I_ERI_Px_S_Px_D2y_C1001001_dd;
  Double I_ERI_Py_Py_Px_D2y_C1001001_dd = I_ERI_D2y_S_Px_D2y_C1001001_dd+ABY*I_ERI_Py_S_Px_D2y_C1001001_dd;
  Double I_ERI_Pz_Py_Px_D2y_C1001001_dd = I_ERI_Dyz_S_Px_D2y_C1001001_dd+ABY*I_ERI_Pz_S_Px_D2y_C1001001_dd;
  Double I_ERI_Px_Pz_Px_D2y_C1001001_dd = I_ERI_Dxz_S_Px_D2y_C1001001_dd+ABZ*I_ERI_Px_S_Px_D2y_C1001001_dd;
  Double I_ERI_Py_Pz_Px_D2y_C1001001_dd = I_ERI_Dyz_S_Px_D2y_C1001001_dd+ABZ*I_ERI_Py_S_Px_D2y_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_D2y_C1001001_dd = I_ERI_D2z_S_Px_D2y_C1001001_dd+ABZ*I_ERI_Pz_S_Px_D2y_C1001001_dd;
  Double I_ERI_Px_Px_Py_D2y_C1001001_dd = I_ERI_D2x_S_Py_D2y_C1001001_dd+ABX*I_ERI_Px_S_Py_D2y_C1001001_dd;
  Double I_ERI_Py_Px_Py_D2y_C1001001_dd = I_ERI_Dxy_S_Py_D2y_C1001001_dd+ABX*I_ERI_Py_S_Py_D2y_C1001001_dd;
  Double I_ERI_Pz_Px_Py_D2y_C1001001_dd = I_ERI_Dxz_S_Py_D2y_C1001001_dd+ABX*I_ERI_Pz_S_Py_D2y_C1001001_dd;
  Double I_ERI_Px_Py_Py_D2y_C1001001_dd = I_ERI_Dxy_S_Py_D2y_C1001001_dd+ABY*I_ERI_Px_S_Py_D2y_C1001001_dd;
  Double I_ERI_Py_Py_Py_D2y_C1001001_dd = I_ERI_D2y_S_Py_D2y_C1001001_dd+ABY*I_ERI_Py_S_Py_D2y_C1001001_dd;
  Double I_ERI_Pz_Py_Py_D2y_C1001001_dd = I_ERI_Dyz_S_Py_D2y_C1001001_dd+ABY*I_ERI_Pz_S_Py_D2y_C1001001_dd;
  Double I_ERI_Px_Pz_Py_D2y_C1001001_dd = I_ERI_Dxz_S_Py_D2y_C1001001_dd+ABZ*I_ERI_Px_S_Py_D2y_C1001001_dd;
  Double I_ERI_Py_Pz_Py_D2y_C1001001_dd = I_ERI_Dyz_S_Py_D2y_C1001001_dd+ABZ*I_ERI_Py_S_Py_D2y_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_D2y_C1001001_dd = I_ERI_D2z_S_Py_D2y_C1001001_dd+ABZ*I_ERI_Pz_S_Py_D2y_C1001001_dd;
  Double I_ERI_Px_Px_Pz_D2y_C1001001_dd = I_ERI_D2x_S_Pz_D2y_C1001001_dd+ABX*I_ERI_Px_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Py_Px_Pz_D2y_C1001001_dd = I_ERI_Dxy_S_Pz_D2y_C1001001_dd+ABX*I_ERI_Py_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_D2y_C1001001_dd = I_ERI_Dxz_S_Pz_D2y_C1001001_dd+ABX*I_ERI_Pz_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Px_Py_Pz_D2y_C1001001_dd = I_ERI_Dxy_S_Pz_D2y_C1001001_dd+ABY*I_ERI_Px_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Py_Py_Pz_D2y_C1001001_dd = I_ERI_D2y_S_Pz_D2y_C1001001_dd+ABY*I_ERI_Py_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_D2y_C1001001_dd = I_ERI_Dyz_S_Pz_D2y_C1001001_dd+ABY*I_ERI_Pz_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_D2y_C1001001_dd = I_ERI_Dxz_S_Pz_D2y_C1001001_dd+ABZ*I_ERI_Px_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_D2y_C1001001_dd = I_ERI_Dyz_S_Pz_D2y_C1001001_dd+ABZ*I_ERI_Py_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_D2y_C1001001_dd = I_ERI_D2z_S_Pz_D2y_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_D2y_C1001001_dd;
  Double I_ERI_Px_Px_Px_Dyz_C1001001_dd = I_ERI_D2x_S_Px_Dyz_C1001001_dd+ABX*I_ERI_Px_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Py_Px_Px_Dyz_C1001001_dd = I_ERI_Dxy_S_Px_Dyz_C1001001_dd+ABX*I_ERI_Py_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Pz_Px_Px_Dyz_C1001001_dd = I_ERI_Dxz_S_Px_Dyz_C1001001_dd+ABX*I_ERI_Pz_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Px_Py_Px_Dyz_C1001001_dd = I_ERI_Dxy_S_Px_Dyz_C1001001_dd+ABY*I_ERI_Px_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Py_Py_Px_Dyz_C1001001_dd = I_ERI_D2y_S_Px_Dyz_C1001001_dd+ABY*I_ERI_Py_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Pz_Py_Px_Dyz_C1001001_dd = I_ERI_Dyz_S_Px_Dyz_C1001001_dd+ABY*I_ERI_Pz_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Px_Pz_Px_Dyz_C1001001_dd = I_ERI_Dxz_S_Px_Dyz_C1001001_dd+ABZ*I_ERI_Px_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Py_Pz_Px_Dyz_C1001001_dd = I_ERI_Dyz_S_Px_Dyz_C1001001_dd+ABZ*I_ERI_Py_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_Dyz_C1001001_dd = I_ERI_D2z_S_Px_Dyz_C1001001_dd+ABZ*I_ERI_Pz_S_Px_Dyz_C1001001_dd;
  Double I_ERI_Px_Px_Py_Dyz_C1001001_dd = I_ERI_D2x_S_Py_Dyz_C1001001_dd+ABX*I_ERI_Px_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Py_Px_Py_Dyz_C1001001_dd = I_ERI_Dxy_S_Py_Dyz_C1001001_dd+ABX*I_ERI_Py_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Pz_Px_Py_Dyz_C1001001_dd = I_ERI_Dxz_S_Py_Dyz_C1001001_dd+ABX*I_ERI_Pz_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Px_Py_Py_Dyz_C1001001_dd = I_ERI_Dxy_S_Py_Dyz_C1001001_dd+ABY*I_ERI_Px_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Py_Py_Py_Dyz_C1001001_dd = I_ERI_D2y_S_Py_Dyz_C1001001_dd+ABY*I_ERI_Py_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Pz_Py_Py_Dyz_C1001001_dd = I_ERI_Dyz_S_Py_Dyz_C1001001_dd+ABY*I_ERI_Pz_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Px_Pz_Py_Dyz_C1001001_dd = I_ERI_Dxz_S_Py_Dyz_C1001001_dd+ABZ*I_ERI_Px_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Py_Pz_Py_Dyz_C1001001_dd = I_ERI_Dyz_S_Py_Dyz_C1001001_dd+ABZ*I_ERI_Py_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_Dyz_C1001001_dd = I_ERI_D2z_S_Py_Dyz_C1001001_dd+ABZ*I_ERI_Pz_S_Py_Dyz_C1001001_dd;
  Double I_ERI_Px_Px_Pz_Dyz_C1001001_dd = I_ERI_D2x_S_Pz_Dyz_C1001001_dd+ABX*I_ERI_Px_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Py_Px_Pz_Dyz_C1001001_dd = I_ERI_Dxy_S_Pz_Dyz_C1001001_dd+ABX*I_ERI_Py_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_Dyz_C1001001_dd = I_ERI_Dxz_S_Pz_Dyz_C1001001_dd+ABX*I_ERI_Pz_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Px_Py_Pz_Dyz_C1001001_dd = I_ERI_Dxy_S_Pz_Dyz_C1001001_dd+ABY*I_ERI_Px_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Py_Py_Pz_Dyz_C1001001_dd = I_ERI_D2y_S_Pz_Dyz_C1001001_dd+ABY*I_ERI_Py_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_Dyz_C1001001_dd = I_ERI_Dyz_S_Pz_Dyz_C1001001_dd+ABY*I_ERI_Pz_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_Dyz_C1001001_dd = I_ERI_Dxz_S_Pz_Dyz_C1001001_dd+ABZ*I_ERI_Px_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_Dyz_C1001001_dd = I_ERI_Dyz_S_Pz_Dyz_C1001001_dd+ABZ*I_ERI_Py_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_Dyz_C1001001_dd = I_ERI_D2z_S_Pz_Dyz_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_Dyz_C1001001_dd;
  Double I_ERI_Px_Px_Px_D2z_C1001001_dd = I_ERI_D2x_S_Px_D2z_C1001001_dd+ABX*I_ERI_Px_S_Px_D2z_C1001001_dd;
  Double I_ERI_Py_Px_Px_D2z_C1001001_dd = I_ERI_Dxy_S_Px_D2z_C1001001_dd+ABX*I_ERI_Py_S_Px_D2z_C1001001_dd;
  Double I_ERI_Pz_Px_Px_D2z_C1001001_dd = I_ERI_Dxz_S_Px_D2z_C1001001_dd+ABX*I_ERI_Pz_S_Px_D2z_C1001001_dd;
  Double I_ERI_Px_Py_Px_D2z_C1001001_dd = I_ERI_Dxy_S_Px_D2z_C1001001_dd+ABY*I_ERI_Px_S_Px_D2z_C1001001_dd;
  Double I_ERI_Py_Py_Px_D2z_C1001001_dd = I_ERI_D2y_S_Px_D2z_C1001001_dd+ABY*I_ERI_Py_S_Px_D2z_C1001001_dd;
  Double I_ERI_Pz_Py_Px_D2z_C1001001_dd = I_ERI_Dyz_S_Px_D2z_C1001001_dd+ABY*I_ERI_Pz_S_Px_D2z_C1001001_dd;
  Double I_ERI_Px_Pz_Px_D2z_C1001001_dd = I_ERI_Dxz_S_Px_D2z_C1001001_dd+ABZ*I_ERI_Px_S_Px_D2z_C1001001_dd;
  Double I_ERI_Py_Pz_Px_D2z_C1001001_dd = I_ERI_Dyz_S_Px_D2z_C1001001_dd+ABZ*I_ERI_Py_S_Px_D2z_C1001001_dd;
  Double I_ERI_Pz_Pz_Px_D2z_C1001001_dd = I_ERI_D2z_S_Px_D2z_C1001001_dd+ABZ*I_ERI_Pz_S_Px_D2z_C1001001_dd;
  Double I_ERI_Px_Px_Py_D2z_C1001001_dd = I_ERI_D2x_S_Py_D2z_C1001001_dd+ABX*I_ERI_Px_S_Py_D2z_C1001001_dd;
  Double I_ERI_Py_Px_Py_D2z_C1001001_dd = I_ERI_Dxy_S_Py_D2z_C1001001_dd+ABX*I_ERI_Py_S_Py_D2z_C1001001_dd;
  Double I_ERI_Pz_Px_Py_D2z_C1001001_dd = I_ERI_Dxz_S_Py_D2z_C1001001_dd+ABX*I_ERI_Pz_S_Py_D2z_C1001001_dd;
  Double I_ERI_Px_Py_Py_D2z_C1001001_dd = I_ERI_Dxy_S_Py_D2z_C1001001_dd+ABY*I_ERI_Px_S_Py_D2z_C1001001_dd;
  Double I_ERI_Py_Py_Py_D2z_C1001001_dd = I_ERI_D2y_S_Py_D2z_C1001001_dd+ABY*I_ERI_Py_S_Py_D2z_C1001001_dd;
  Double I_ERI_Pz_Py_Py_D2z_C1001001_dd = I_ERI_Dyz_S_Py_D2z_C1001001_dd+ABY*I_ERI_Pz_S_Py_D2z_C1001001_dd;
  Double I_ERI_Px_Pz_Py_D2z_C1001001_dd = I_ERI_Dxz_S_Py_D2z_C1001001_dd+ABZ*I_ERI_Px_S_Py_D2z_C1001001_dd;
  Double I_ERI_Py_Pz_Py_D2z_C1001001_dd = I_ERI_Dyz_S_Py_D2z_C1001001_dd+ABZ*I_ERI_Py_S_Py_D2z_C1001001_dd;
  Double I_ERI_Pz_Pz_Py_D2z_C1001001_dd = I_ERI_D2z_S_Py_D2z_C1001001_dd+ABZ*I_ERI_Pz_S_Py_D2z_C1001001_dd;
  Double I_ERI_Px_Px_Pz_D2z_C1001001_dd = I_ERI_D2x_S_Pz_D2z_C1001001_dd+ABX*I_ERI_Px_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Py_Px_Pz_D2z_C1001001_dd = I_ERI_Dxy_S_Pz_D2z_C1001001_dd+ABX*I_ERI_Py_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Pz_Px_Pz_D2z_C1001001_dd = I_ERI_Dxz_S_Pz_D2z_C1001001_dd+ABX*I_ERI_Pz_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Px_Py_Pz_D2z_C1001001_dd = I_ERI_Dxy_S_Pz_D2z_C1001001_dd+ABY*I_ERI_Px_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Py_Py_Pz_D2z_C1001001_dd = I_ERI_D2y_S_Pz_D2z_C1001001_dd+ABY*I_ERI_Py_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Pz_Py_Pz_D2z_C1001001_dd = I_ERI_Dyz_S_Pz_D2z_C1001001_dd+ABY*I_ERI_Pz_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Px_Pz_Pz_D2z_C1001001_dd = I_ERI_Dxz_S_Pz_D2z_C1001001_dd+ABZ*I_ERI_Px_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Py_Pz_Pz_D2z_C1001001_dd = I_ERI_Dyz_S_Pz_D2z_C1001001_dd+ABZ*I_ERI_Py_S_Pz_D2z_C1001001_dd;
  Double I_ERI_Pz_Pz_Pz_D2z_C1001001_dd = I_ERI_D2z_S_Pz_D2z_C1001001_dd+ABZ*I_ERI_Pz_S_Pz_D2z_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_F3x_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_a-2.0E0*2*I_ERI_Px_S_Px_S_C1000001_a;
  abcd[1] = 4.0E0*I_ERI_F2xy_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_a;
  abcd[2] = 4.0E0*I_ERI_F2xz_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_a;
  abcd[12] = 4.0E0*I_ERI_F3x_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_a-2.0E0*2*I_ERI_Px_S_Py_S_C1000001_a;
  abcd[13] = 4.0E0*I_ERI_F2xy_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_a;
  abcd[14] = 4.0E0*I_ERI_F2xz_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_a;
  abcd[24] = 4.0E0*I_ERI_F3x_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_a-2.0E0*2*I_ERI_Px_S_Pz_S_C1000001_a;
  abcd[25] = 4.0E0*I_ERI_F2xy_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_a;
  abcd[26] = 4.0E0*I_ERI_F2xz_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[3] = 4.0E0*I_ERI_F3x_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_a-2.0E0*2*I_ERI_Px_Px_Px_S_C1001001_a;
  abcd[4] = 4.0E0*I_ERI_F2xy_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_a;
  abcd[5] = 4.0E0*I_ERI_F2xz_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_a;
  abcd[6] = 4.0E0*I_ERI_F3x_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_a-2.0E0*2*I_ERI_Px_Py_Px_S_C1001001_a;
  abcd[7] = 4.0E0*I_ERI_F2xy_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_a;
  abcd[8] = 4.0E0*I_ERI_F2xz_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_a;
  abcd[9] = 4.0E0*I_ERI_F3x_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_a-2.0E0*2*I_ERI_Px_Pz_Px_S_C1001001_a;
  abcd[10] = 4.0E0*I_ERI_F2xy_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_a;
  abcd[11] = 4.0E0*I_ERI_F2xz_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_a;
  abcd[15] = 4.0E0*I_ERI_F3x_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_a-2.0E0*2*I_ERI_Px_Px_Py_S_C1001001_a;
  abcd[16] = 4.0E0*I_ERI_F2xy_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_a;
  abcd[17] = 4.0E0*I_ERI_F2xz_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_a;
  abcd[18] = 4.0E0*I_ERI_F3x_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_a-2.0E0*2*I_ERI_Px_Py_Py_S_C1001001_a;
  abcd[19] = 4.0E0*I_ERI_F2xy_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_a;
  abcd[20] = 4.0E0*I_ERI_F2xz_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_a;
  abcd[21] = 4.0E0*I_ERI_F3x_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_a-2.0E0*2*I_ERI_Px_Pz_Py_S_C1001001_a;
  abcd[22] = 4.0E0*I_ERI_F2xy_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_a;
  abcd[23] = 4.0E0*I_ERI_F2xz_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_a;
  abcd[27] = 4.0E0*I_ERI_F3x_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_a-2.0E0*2*I_ERI_Px_Px_Pz_S_C1001001_a;
  abcd[28] = 4.0E0*I_ERI_F2xy_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_a;
  abcd[29] = 4.0E0*I_ERI_F2xz_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_a;
  abcd[30] = 4.0E0*I_ERI_F3x_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_a-2.0E0*2*I_ERI_Px_Py_Pz_S_C1001001_a;
  abcd[31] = 4.0E0*I_ERI_F2xy_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_a;
  abcd[32] = 4.0E0*I_ERI_F2xz_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_a;
  abcd[33] = 4.0E0*I_ERI_F3x_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_a-2.0E0*2*I_ERI_Px_Pz_Pz_S_C1001001_a;
  abcd[34] = 4.0E0*I_ERI_F2xy_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_a;
  abcd[35] = 4.0E0*I_ERI_F2xz_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[36] = 4.0E0*I_ERI_F2xy_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_a;
  abcd[37] = 4.0E0*I_ERI_Fx2y_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_a;
  abcd[38] = 4.0E0*I_ERI_Fxyz_S_Px_S_C1000001_aa;
  abcd[48] = 4.0E0*I_ERI_F2xy_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_a;
  abcd[49] = 4.0E0*I_ERI_Fx2y_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_a;
  abcd[50] = 4.0E0*I_ERI_Fxyz_S_Py_S_C1000001_aa;
  abcd[60] = 4.0E0*I_ERI_F2xy_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_a;
  abcd[61] = 4.0E0*I_ERI_Fx2y_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_a;
  abcd[62] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C1000001_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[39] = 4.0E0*I_ERI_F2xy_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_a;
  abcd[40] = 4.0E0*I_ERI_Fx2y_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_a;
  abcd[41] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1001001_aa;
  abcd[42] = 4.0E0*I_ERI_F2xy_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_a;
  abcd[43] = 4.0E0*I_ERI_Fx2y_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_a;
  abcd[44] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1001001_aa;
  abcd[45] = 4.0E0*I_ERI_F2xy_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_a;
  abcd[46] = 4.0E0*I_ERI_Fx2y_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_a;
  abcd[47] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1001001_aa;
  abcd[51] = 4.0E0*I_ERI_F2xy_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_a;
  abcd[52] = 4.0E0*I_ERI_Fx2y_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_a;
  abcd[53] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1001001_aa;
  abcd[54] = 4.0E0*I_ERI_F2xy_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_a;
  abcd[55] = 4.0E0*I_ERI_Fx2y_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_a;
  abcd[56] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1001001_aa;
  abcd[57] = 4.0E0*I_ERI_F2xy_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_a;
  abcd[58] = 4.0E0*I_ERI_Fx2y_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_a;
  abcd[59] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1001001_aa;
  abcd[63] = 4.0E0*I_ERI_F2xy_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_a;
  abcd[64] = 4.0E0*I_ERI_Fx2y_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_a;
  abcd[65] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1001001_aa;
  abcd[66] = 4.0E0*I_ERI_F2xy_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_a;
  abcd[67] = 4.0E0*I_ERI_Fx2y_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_a;
  abcd[68] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1001001_aa;
  abcd[69] = 4.0E0*I_ERI_F2xy_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_a;
  abcd[70] = 4.0E0*I_ERI_Fx2y_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_a;
  abcd[71] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001001_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[72] = 4.0E0*I_ERI_F2xz_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_a;
  abcd[73] = 4.0E0*I_ERI_Fxyz_S_Px_S_C1000001_aa;
  abcd[74] = 4.0E0*I_ERI_Fx2z_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_a;
  abcd[84] = 4.0E0*I_ERI_F2xz_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_a;
  abcd[85] = 4.0E0*I_ERI_Fxyz_S_Py_S_C1000001_aa;
  abcd[86] = 4.0E0*I_ERI_Fx2z_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_a;
  abcd[96] = 4.0E0*I_ERI_F2xz_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_a;
  abcd[97] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C1000001_aa;
  abcd[98] = 4.0E0*I_ERI_Fx2z_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[75] = 4.0E0*I_ERI_F2xz_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_a;
  abcd[76] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1001001_aa;
  abcd[77] = 4.0E0*I_ERI_Fx2z_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_a;
  abcd[78] = 4.0E0*I_ERI_F2xz_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_a;
  abcd[79] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1001001_aa;
  abcd[80] = 4.0E0*I_ERI_Fx2z_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_a;
  abcd[81] = 4.0E0*I_ERI_F2xz_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_a;
  abcd[82] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1001001_aa;
  abcd[83] = 4.0E0*I_ERI_Fx2z_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_a;
  abcd[87] = 4.0E0*I_ERI_F2xz_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_a;
  abcd[88] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1001001_aa;
  abcd[89] = 4.0E0*I_ERI_Fx2z_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_a;
  abcd[90] = 4.0E0*I_ERI_F2xz_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_a;
  abcd[91] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1001001_aa;
  abcd[92] = 4.0E0*I_ERI_Fx2z_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_a;
  abcd[93] = 4.0E0*I_ERI_F2xz_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_a;
  abcd[94] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1001001_aa;
  abcd[95] = 4.0E0*I_ERI_Fx2z_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_a;
  abcd[99] = 4.0E0*I_ERI_F2xz_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_a;
  abcd[100] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1001001_aa;
  abcd[101] = 4.0E0*I_ERI_Fx2z_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_a;
  abcd[102] = 4.0E0*I_ERI_F2xz_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_a;
  abcd[103] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1001001_aa;
  abcd[104] = 4.0E0*I_ERI_Fx2z_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_a;
  abcd[105] = 4.0E0*I_ERI_F2xz_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_a;
  abcd[106] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001001_aa;
  abcd[107] = 4.0E0*I_ERI_Fx2z_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[108] = 4.0E0*I_ERI_Fx2y_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_a;
  abcd[109] = 4.0E0*I_ERI_F3y_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_a-2.0E0*2*I_ERI_Py_S_Px_S_C1000001_a;
  abcd[110] = 4.0E0*I_ERI_F2yz_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_a;
  abcd[120] = 4.0E0*I_ERI_Fx2y_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_a;
  abcd[121] = 4.0E0*I_ERI_F3y_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_a-2.0E0*2*I_ERI_Py_S_Py_S_C1000001_a;
  abcd[122] = 4.0E0*I_ERI_F2yz_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_a;
  abcd[132] = 4.0E0*I_ERI_Fx2y_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_a;
  abcd[133] = 4.0E0*I_ERI_F3y_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_a-2.0E0*2*I_ERI_Py_S_Pz_S_C1000001_a;
  abcd[134] = 4.0E0*I_ERI_F2yz_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[111] = 4.0E0*I_ERI_Fx2y_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_a;
  abcd[112] = 4.0E0*I_ERI_F3y_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_a-2.0E0*2*I_ERI_Py_Px_Px_S_C1001001_a;
  abcd[113] = 4.0E0*I_ERI_F2yz_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_a;
  abcd[114] = 4.0E0*I_ERI_Fx2y_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_a;
  abcd[115] = 4.0E0*I_ERI_F3y_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_a-2.0E0*2*I_ERI_Py_Py_Px_S_C1001001_a;
  abcd[116] = 4.0E0*I_ERI_F2yz_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_a;
  abcd[117] = 4.0E0*I_ERI_Fx2y_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_a;
  abcd[118] = 4.0E0*I_ERI_F3y_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_a-2.0E0*2*I_ERI_Py_Pz_Px_S_C1001001_a;
  abcd[119] = 4.0E0*I_ERI_F2yz_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_a;
  abcd[123] = 4.0E0*I_ERI_Fx2y_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_a;
  abcd[124] = 4.0E0*I_ERI_F3y_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_a-2.0E0*2*I_ERI_Py_Px_Py_S_C1001001_a;
  abcd[125] = 4.0E0*I_ERI_F2yz_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_a;
  abcd[126] = 4.0E0*I_ERI_Fx2y_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_a;
  abcd[127] = 4.0E0*I_ERI_F3y_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_a-2.0E0*2*I_ERI_Py_Py_Py_S_C1001001_a;
  abcd[128] = 4.0E0*I_ERI_F2yz_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_a;
  abcd[129] = 4.0E0*I_ERI_Fx2y_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_a;
  abcd[130] = 4.0E0*I_ERI_F3y_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_a-2.0E0*2*I_ERI_Py_Pz_Py_S_C1001001_a;
  abcd[131] = 4.0E0*I_ERI_F2yz_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_a;
  abcd[135] = 4.0E0*I_ERI_Fx2y_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_a;
  abcd[136] = 4.0E0*I_ERI_F3y_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_a-2.0E0*2*I_ERI_Py_Px_Pz_S_C1001001_a;
  abcd[137] = 4.0E0*I_ERI_F2yz_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_a;
  abcd[138] = 4.0E0*I_ERI_Fx2y_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_a;
  abcd[139] = 4.0E0*I_ERI_F3y_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_a-2.0E0*2*I_ERI_Py_Py_Pz_S_C1001001_a;
  abcd[140] = 4.0E0*I_ERI_F2yz_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_a;
  abcd[141] = 4.0E0*I_ERI_Fx2y_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_a;
  abcd[142] = 4.0E0*I_ERI_F3y_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_a-2.0E0*2*I_ERI_Py_Pz_Pz_S_C1001001_a;
  abcd[143] = 4.0E0*I_ERI_F2yz_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[144] = 4.0E0*I_ERI_Fxyz_S_Px_S_C1000001_aa;
  abcd[145] = 4.0E0*I_ERI_F2yz_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_a;
  abcd[146] = 4.0E0*I_ERI_Fy2z_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_a;
  abcd[156] = 4.0E0*I_ERI_Fxyz_S_Py_S_C1000001_aa;
  abcd[157] = 4.0E0*I_ERI_F2yz_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_a;
  abcd[158] = 4.0E0*I_ERI_Fy2z_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_a;
  abcd[168] = 4.0E0*I_ERI_Fxyz_S_Pz_S_C1000001_aa;
  abcd[169] = 4.0E0*I_ERI_F2yz_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_a;
  abcd[170] = 4.0E0*I_ERI_Fy2z_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[147] = 4.0E0*I_ERI_Fxyz_Px_Px_S_C1001001_aa;
  abcd[148] = 4.0E0*I_ERI_F2yz_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_a;
  abcd[149] = 4.0E0*I_ERI_Fy2z_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_a;
  abcd[150] = 4.0E0*I_ERI_Fxyz_Py_Px_S_C1001001_aa;
  abcd[151] = 4.0E0*I_ERI_F2yz_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_a;
  abcd[152] = 4.0E0*I_ERI_Fy2z_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_a;
  abcd[153] = 4.0E0*I_ERI_Fxyz_Pz_Px_S_C1001001_aa;
  abcd[154] = 4.0E0*I_ERI_F2yz_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_a;
  abcd[155] = 4.0E0*I_ERI_Fy2z_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_a;
  abcd[159] = 4.0E0*I_ERI_Fxyz_Px_Py_S_C1001001_aa;
  abcd[160] = 4.0E0*I_ERI_F2yz_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_a;
  abcd[161] = 4.0E0*I_ERI_Fy2z_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_a;
  abcd[162] = 4.0E0*I_ERI_Fxyz_Py_Py_S_C1001001_aa;
  abcd[163] = 4.0E0*I_ERI_F2yz_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_a;
  abcd[164] = 4.0E0*I_ERI_Fy2z_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_a;
  abcd[165] = 4.0E0*I_ERI_Fxyz_Pz_Py_S_C1001001_aa;
  abcd[166] = 4.0E0*I_ERI_F2yz_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_a;
  abcd[167] = 4.0E0*I_ERI_Fy2z_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_a;
  abcd[171] = 4.0E0*I_ERI_Fxyz_Px_Pz_S_C1001001_aa;
  abcd[172] = 4.0E0*I_ERI_F2yz_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_a;
  abcd[173] = 4.0E0*I_ERI_Fy2z_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_a;
  abcd[174] = 4.0E0*I_ERI_Fxyz_Py_Pz_S_C1001001_aa;
  abcd[175] = 4.0E0*I_ERI_F2yz_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_a;
  abcd[176] = 4.0E0*I_ERI_Fy2z_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_a;
  abcd[177] = 4.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001001_aa;
  abcd[178] = 4.0E0*I_ERI_F2yz_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_a;
  abcd[179] = 4.0E0*I_ERI_Fy2z_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000001_aa
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_a
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_Fx2z_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_a;
  abcd[181] = 4.0E0*I_ERI_Fy2z_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_a;
  abcd[182] = 4.0E0*I_ERI_F3z_S_Px_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_a-2.0E0*2*I_ERI_Pz_S_Px_S_C1000001_a;
  abcd[192] = 4.0E0*I_ERI_Fx2z_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_a;
  abcd[193] = 4.0E0*I_ERI_Fy2z_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_a;
  abcd[194] = 4.0E0*I_ERI_F3z_S_Py_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_a-2.0E0*2*I_ERI_Pz_S_Py_S_C1000001_a;
  abcd[204] = 4.0E0*I_ERI_Fx2z_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_a;
  abcd[205] = 4.0E0*I_ERI_Fy2z_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_a;
  abcd[206] = 4.0E0*I_ERI_F3z_S_Pz_S_C1000001_aa-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_a-2.0E0*2*I_ERI_Pz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001001_aa
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_a
   ************************************************************/
  abcd[183] = 4.0E0*I_ERI_Fx2z_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_a;
  abcd[184] = 4.0E0*I_ERI_Fy2z_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_a;
  abcd[185] = 4.0E0*I_ERI_F3z_Px_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_a-2.0E0*2*I_ERI_Pz_Px_Px_S_C1001001_a;
  abcd[186] = 4.0E0*I_ERI_Fx2z_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_a;
  abcd[187] = 4.0E0*I_ERI_Fy2z_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_a;
  abcd[188] = 4.0E0*I_ERI_F3z_Py_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_a-2.0E0*2*I_ERI_Pz_Py_Px_S_C1001001_a;
  abcd[189] = 4.0E0*I_ERI_Fx2z_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_a;
  abcd[190] = 4.0E0*I_ERI_Fy2z_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_a;
  abcd[191] = 4.0E0*I_ERI_F3z_Pz_Px_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_a-2.0E0*2*I_ERI_Pz_Pz_Px_S_C1001001_a;
  abcd[195] = 4.0E0*I_ERI_Fx2z_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_a;
  abcd[196] = 4.0E0*I_ERI_Fy2z_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_a;
  abcd[197] = 4.0E0*I_ERI_F3z_Px_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_a-2.0E0*2*I_ERI_Pz_Px_Py_S_C1001001_a;
  abcd[198] = 4.0E0*I_ERI_Fx2z_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_a;
  abcd[199] = 4.0E0*I_ERI_Fy2z_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_a;
  abcd[200] = 4.0E0*I_ERI_F3z_Py_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_a-2.0E0*2*I_ERI_Pz_Py_Py_S_C1001001_a;
  abcd[201] = 4.0E0*I_ERI_Fx2z_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_a;
  abcd[202] = 4.0E0*I_ERI_Fy2z_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_a;
  abcd[203] = 4.0E0*I_ERI_F3z_Pz_Py_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_a-2.0E0*2*I_ERI_Pz_Pz_Py_S_C1001001_a;
  abcd[207] = 4.0E0*I_ERI_Fx2z_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_a;
  abcd[208] = 4.0E0*I_ERI_Fy2z_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_a;
  abcd[209] = 4.0E0*I_ERI_F3z_Px_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_a-2.0E0*2*I_ERI_Pz_Px_Pz_S_C1001001_a;
  abcd[210] = 4.0E0*I_ERI_Fx2z_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_a;
  abcd[211] = 4.0E0*I_ERI_Fy2z_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_a;
  abcd[212] = 4.0E0*I_ERI_F3z_Py_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_a-2.0E0*2*I_ERI_Pz_Py_Pz_S_C1001001_a;
  abcd[213] = 4.0E0*I_ERI_Fx2z_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_a;
  abcd[214] = 4.0E0*I_ERI_Fy2z_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_a;
  abcd[215] = 4.0E0*I_ERI_F3z_Pz_Pz_S_C1001001_aa-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_a-2.0E0*2*I_ERI_Pz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[216] = 4.0E0*I_ERI_D2x_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_D2x_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2x_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[217] = 4.0E0*I_ERI_Dxy_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[218] = 4.0E0*I_ERI_Dxz_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;
  abcd[228] = 4.0E0*I_ERI_D2x_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[229] = 4.0E0*I_ERI_Dxy_S_Dxy_S_C1000001_ac;
  abcd[230] = 4.0E0*I_ERI_Dxz_S_Dxy_S_C1000001_ac;
  abcd[240] = 4.0E0*I_ERI_D2x_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;
  abcd[241] = 4.0E0*I_ERI_Dxy_S_Dxz_S_C1000001_ac;
  abcd[242] = 4.0E0*I_ERI_Dxz_S_Dxz_S_C1000001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[219] = 4.0E0*I_ERI_D2x_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2x_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[220] = 4.0E0*I_ERI_Dxy_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[221] = 4.0E0*I_ERI_Dxz_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[222] = 4.0E0*I_ERI_D2x_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2x_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[223] = 4.0E0*I_ERI_Dxy_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[224] = 4.0E0*I_ERI_Dxz_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[225] = 4.0E0*I_ERI_D2x_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2x_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[226] = 4.0E0*I_ERI_Dxy_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[227] = 4.0E0*I_ERI_Dxz_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;
  abcd[231] = 4.0E0*I_ERI_D2x_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[232] = 4.0E0*I_ERI_Dxy_Px_Dxy_S_C1001001_ac;
  abcd[233] = 4.0E0*I_ERI_Dxz_Px_Dxy_S_C1001001_ac;
  abcd[234] = 4.0E0*I_ERI_D2x_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[235] = 4.0E0*I_ERI_Dxy_Py_Dxy_S_C1001001_ac;
  abcd[236] = 4.0E0*I_ERI_Dxz_Py_Dxy_S_C1001001_ac;
  abcd[237] = 4.0E0*I_ERI_D2x_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[238] = 4.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001001_ac;
  abcd[239] = 4.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001001_ac;
  abcd[243] = 4.0E0*I_ERI_D2x_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[244] = 4.0E0*I_ERI_Dxy_Px_Dxz_S_C1001001_ac;
  abcd[245] = 4.0E0*I_ERI_Dxz_Px_Dxz_S_C1001001_ac;
  abcd[246] = 4.0E0*I_ERI_D2x_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[247] = 4.0E0*I_ERI_Dxy_Py_Dxz_S_C1001001_ac;
  abcd[248] = 4.0E0*I_ERI_Dxz_Py_Dxz_S_C1001001_ac;
  abcd[249] = 4.0E0*I_ERI_D2x_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;
  abcd[250] = 4.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001001_ac;
  abcd[251] = 4.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[252] = 4.0E0*I_ERI_D2x_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[253] = 4.0E0*I_ERI_Dxy_S_Dxy_S_C1000001_ac;
  abcd[254] = 4.0E0*I_ERI_Dxz_S_Dxy_S_C1000001_ac;
  abcd[264] = 4.0E0*I_ERI_D2x_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_D2x_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2y_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[265] = 4.0E0*I_ERI_Dxy_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[266] = 4.0E0*I_ERI_Dxz_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;
  abcd[276] = 4.0E0*I_ERI_D2x_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;
  abcd[277] = 4.0E0*I_ERI_Dxy_S_Dyz_S_C1000001_ac;
  abcd[278] = 4.0E0*I_ERI_Dxz_S_Dyz_S_C1000001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[255] = 4.0E0*I_ERI_D2x_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[256] = 4.0E0*I_ERI_Dxy_Px_Dxy_S_C1001001_ac;
  abcd[257] = 4.0E0*I_ERI_Dxz_Px_Dxy_S_C1001001_ac;
  abcd[258] = 4.0E0*I_ERI_D2x_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[259] = 4.0E0*I_ERI_Dxy_Py_Dxy_S_C1001001_ac;
  abcd[260] = 4.0E0*I_ERI_Dxz_Py_Dxy_S_C1001001_ac;
  abcd[261] = 4.0E0*I_ERI_D2x_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[262] = 4.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001001_ac;
  abcd[263] = 4.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001001_ac;
  abcd[267] = 4.0E0*I_ERI_D2x_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2y_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[268] = 4.0E0*I_ERI_Dxy_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[269] = 4.0E0*I_ERI_Dxz_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[270] = 4.0E0*I_ERI_D2x_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2y_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[271] = 4.0E0*I_ERI_Dxy_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[272] = 4.0E0*I_ERI_Dxz_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[273] = 4.0E0*I_ERI_D2x_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2y_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[274] = 4.0E0*I_ERI_Dxy_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[275] = 4.0E0*I_ERI_Dxz_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;
  abcd[279] = 4.0E0*I_ERI_D2x_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[280] = 4.0E0*I_ERI_Dxy_Px_Dyz_S_C1001001_ac;
  abcd[281] = 4.0E0*I_ERI_Dxz_Px_Dyz_S_C1001001_ac;
  abcd[282] = 4.0E0*I_ERI_D2x_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[283] = 4.0E0*I_ERI_Dxy_Py_Dyz_S_C1001001_ac;
  abcd[284] = 4.0E0*I_ERI_Dxz_Py_Dyz_S_C1001001_ac;
  abcd[285] = 4.0E0*I_ERI_D2x_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;
  abcd[286] = 4.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001001_ac;
  abcd[287] = 4.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[288] = 4.0E0*I_ERI_D2x_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;
  abcd[289] = 4.0E0*I_ERI_Dxy_S_Dxz_S_C1000001_ac;
  abcd[290] = 4.0E0*I_ERI_Dxz_S_Dxz_S_C1000001_ac;
  abcd[300] = 4.0E0*I_ERI_D2x_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;
  abcd[301] = 4.0E0*I_ERI_Dxy_S_Dyz_S_C1000001_ac;
  abcd[302] = 4.0E0*I_ERI_Dxz_S_Dyz_S_C1000001_ac;
  abcd[312] = 4.0E0*I_ERI_D2x_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_D2x_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2z_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[313] = 4.0E0*I_ERI_Dxy_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[314] = 4.0E0*I_ERI_Dxz_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[291] = 4.0E0*I_ERI_D2x_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[292] = 4.0E0*I_ERI_Dxy_Px_Dxz_S_C1001001_ac;
  abcd[293] = 4.0E0*I_ERI_Dxz_Px_Dxz_S_C1001001_ac;
  abcd[294] = 4.0E0*I_ERI_D2x_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[295] = 4.0E0*I_ERI_Dxy_Py_Dxz_S_C1001001_ac;
  abcd[296] = 4.0E0*I_ERI_Dxz_Py_Dxz_S_C1001001_ac;
  abcd[297] = 4.0E0*I_ERI_D2x_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;
  abcd[298] = 4.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001001_ac;
  abcd[299] = 4.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001001_ac;
  abcd[303] = 4.0E0*I_ERI_D2x_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[304] = 4.0E0*I_ERI_Dxy_Px_Dyz_S_C1001001_ac;
  abcd[305] = 4.0E0*I_ERI_Dxz_Px_Dyz_S_C1001001_ac;
  abcd[306] = 4.0E0*I_ERI_D2x_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[307] = 4.0E0*I_ERI_Dxy_Py_Dyz_S_C1001001_ac;
  abcd[308] = 4.0E0*I_ERI_Dxz_Py_Dyz_S_C1001001_ac;
  abcd[309] = 4.0E0*I_ERI_D2x_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;
  abcd[310] = 4.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001001_ac;
  abcd[311] = 4.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001001_ac;
  abcd[315] = 4.0E0*I_ERI_D2x_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2z_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[316] = 4.0E0*I_ERI_Dxy_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[317] = 4.0E0*I_ERI_Dxz_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[318] = 4.0E0*I_ERI_D2x_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2z_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[319] = 4.0E0*I_ERI_Dxy_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[320] = 4.0E0*I_ERI_Dxz_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[321] = 4.0E0*I_ERI_D2x_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2x_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2z_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[322] = 4.0E0*I_ERI_Dxy_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[323] = 4.0E0*I_ERI_Dxz_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[324] = 4.0E0*I_ERI_Dxy_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[325] = 4.0E0*I_ERI_D2y_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_D2y_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2x_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[326] = 4.0E0*I_ERI_Dyz_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;
  abcd[336] = 4.0E0*I_ERI_Dxy_S_Dxy_S_C1000001_ac;
  abcd[337] = 4.0E0*I_ERI_D2y_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[338] = 4.0E0*I_ERI_Dyz_S_Dxy_S_C1000001_ac;
  abcd[348] = 4.0E0*I_ERI_Dxy_S_Dxz_S_C1000001_ac;
  abcd[349] = 4.0E0*I_ERI_D2y_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;
  abcd[350] = 4.0E0*I_ERI_Dyz_S_Dxz_S_C1000001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[327] = 4.0E0*I_ERI_Dxy_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[328] = 4.0E0*I_ERI_D2y_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2x_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[329] = 4.0E0*I_ERI_Dyz_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[330] = 4.0E0*I_ERI_Dxy_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[331] = 4.0E0*I_ERI_D2y_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2x_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[332] = 4.0E0*I_ERI_Dyz_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[333] = 4.0E0*I_ERI_Dxy_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[334] = 4.0E0*I_ERI_D2y_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2x_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[335] = 4.0E0*I_ERI_Dyz_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;
  abcd[339] = 4.0E0*I_ERI_Dxy_Px_Dxy_S_C1001001_ac;
  abcd[340] = 4.0E0*I_ERI_D2y_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[341] = 4.0E0*I_ERI_Dyz_Px_Dxy_S_C1001001_ac;
  abcd[342] = 4.0E0*I_ERI_Dxy_Py_Dxy_S_C1001001_ac;
  abcd[343] = 4.0E0*I_ERI_D2y_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[344] = 4.0E0*I_ERI_Dyz_Py_Dxy_S_C1001001_ac;
  abcd[345] = 4.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001001_ac;
  abcd[346] = 4.0E0*I_ERI_D2y_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[347] = 4.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001001_ac;
  abcd[351] = 4.0E0*I_ERI_Dxy_Px_Dxz_S_C1001001_ac;
  abcd[352] = 4.0E0*I_ERI_D2y_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[353] = 4.0E0*I_ERI_Dyz_Px_Dxz_S_C1001001_ac;
  abcd[354] = 4.0E0*I_ERI_Dxy_Py_Dxz_S_C1001001_ac;
  abcd[355] = 4.0E0*I_ERI_D2y_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[356] = 4.0E0*I_ERI_Dyz_Py_Dxz_S_C1001001_ac;
  abcd[357] = 4.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001001_ac;
  abcd[358] = 4.0E0*I_ERI_D2y_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;
  abcd[359] = 4.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[360] = 4.0E0*I_ERI_Dxy_S_Dxy_S_C1000001_ac;
  abcd[361] = 4.0E0*I_ERI_D2y_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[362] = 4.0E0*I_ERI_Dyz_S_Dxy_S_C1000001_ac;
  abcd[372] = 4.0E0*I_ERI_Dxy_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[373] = 4.0E0*I_ERI_D2y_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_D2y_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2y_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[374] = 4.0E0*I_ERI_Dyz_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;
  abcd[384] = 4.0E0*I_ERI_Dxy_S_Dyz_S_C1000001_ac;
  abcd[385] = 4.0E0*I_ERI_D2y_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;
  abcd[386] = 4.0E0*I_ERI_Dyz_S_Dyz_S_C1000001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[363] = 4.0E0*I_ERI_Dxy_Px_Dxy_S_C1001001_ac;
  abcd[364] = 4.0E0*I_ERI_D2y_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[365] = 4.0E0*I_ERI_Dyz_Px_Dxy_S_C1001001_ac;
  abcd[366] = 4.0E0*I_ERI_Dxy_Py_Dxy_S_C1001001_ac;
  abcd[367] = 4.0E0*I_ERI_D2y_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[368] = 4.0E0*I_ERI_Dyz_Py_Dxy_S_C1001001_ac;
  abcd[369] = 4.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001001_ac;
  abcd[370] = 4.0E0*I_ERI_D2y_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[371] = 4.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001001_ac;
  abcd[375] = 4.0E0*I_ERI_Dxy_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[376] = 4.0E0*I_ERI_D2y_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2y_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[377] = 4.0E0*I_ERI_Dyz_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[378] = 4.0E0*I_ERI_Dxy_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[379] = 4.0E0*I_ERI_D2y_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2y_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[380] = 4.0E0*I_ERI_Dyz_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[381] = 4.0E0*I_ERI_Dxy_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[382] = 4.0E0*I_ERI_D2y_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2y_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[383] = 4.0E0*I_ERI_Dyz_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;
  abcd[387] = 4.0E0*I_ERI_Dxy_Px_Dyz_S_C1001001_ac;
  abcd[388] = 4.0E0*I_ERI_D2y_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[389] = 4.0E0*I_ERI_Dyz_Px_Dyz_S_C1001001_ac;
  abcd[390] = 4.0E0*I_ERI_Dxy_Py_Dyz_S_C1001001_ac;
  abcd[391] = 4.0E0*I_ERI_D2y_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[392] = 4.0E0*I_ERI_Dyz_Py_Dyz_S_C1001001_ac;
  abcd[393] = 4.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001001_ac;
  abcd[394] = 4.0E0*I_ERI_D2y_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;
  abcd[395] = 4.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001001_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[396] = 4.0E0*I_ERI_Dxy_S_Dxz_S_C1000001_ac;
  abcd[397] = 4.0E0*I_ERI_D2y_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;
  abcd[398] = 4.0E0*I_ERI_Dyz_S_Dxz_S_C1000001_ac;
  abcd[408] = 4.0E0*I_ERI_Dxy_S_Dyz_S_C1000001_ac;
  abcd[409] = 4.0E0*I_ERI_D2y_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;
  abcd[410] = 4.0E0*I_ERI_Dyz_S_Dyz_S_C1000001_ac;
  abcd[420] = 4.0E0*I_ERI_Dxy_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dxy_S_S_S_C1000001_a;
  abcd[421] = 4.0E0*I_ERI_D2y_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_D2y_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2z_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[422] = 4.0E0*I_ERI_Dyz_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[399] = 4.0E0*I_ERI_Dxy_Px_Dxz_S_C1001001_ac;
  abcd[400] = 4.0E0*I_ERI_D2y_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[401] = 4.0E0*I_ERI_Dyz_Px_Dxz_S_C1001001_ac;
  abcd[402] = 4.0E0*I_ERI_Dxy_Py_Dxz_S_C1001001_ac;
  abcd[403] = 4.0E0*I_ERI_D2y_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[404] = 4.0E0*I_ERI_Dyz_Py_Dxz_S_C1001001_ac;
  abcd[405] = 4.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001001_ac;
  abcd[406] = 4.0E0*I_ERI_D2y_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;
  abcd[407] = 4.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001001_ac;
  abcd[411] = 4.0E0*I_ERI_Dxy_Px_Dyz_S_C1001001_ac;
  abcd[412] = 4.0E0*I_ERI_D2y_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[413] = 4.0E0*I_ERI_Dyz_Px_Dyz_S_C1001001_ac;
  abcd[414] = 4.0E0*I_ERI_Dxy_Py_Dyz_S_C1001001_ac;
  abcd[415] = 4.0E0*I_ERI_D2y_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[416] = 4.0E0*I_ERI_Dyz_Py_Dyz_S_C1001001_ac;
  abcd[417] = 4.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001001_ac;
  abcd[418] = 4.0E0*I_ERI_D2y_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;
  abcd[419] = 4.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001001_ac;
  abcd[423] = 4.0E0*I_ERI_Dxy_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Px_S_S_C1001001_a;
  abcd[424] = 4.0E0*I_ERI_D2y_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2z_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[425] = 4.0E0*I_ERI_Dyz_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[426] = 4.0E0*I_ERI_Dxy_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Py_S_S_C1001001_a;
  abcd[427] = 4.0E0*I_ERI_D2y_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2z_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[428] = 4.0E0*I_ERI_Dyz_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[429] = 4.0E0*I_ERI_Dxy_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxy_Pz_S_S_C1001001_a;
  abcd[430] = 4.0E0*I_ERI_D2y_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2y_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2z_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[431] = 4.0E0*I_ERI_Dyz_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[432] = 4.0E0*I_ERI_Dxz_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;
  abcd[433] = 4.0E0*I_ERI_Dyz_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;
  abcd[434] = 4.0E0*I_ERI_D2z_S_D2x_S_C1000001_ac-2.0E0*1*I_ERI_D2z_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2x_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[444] = 4.0E0*I_ERI_Dxz_S_Dxy_S_C1000001_ac;
  abcd[445] = 4.0E0*I_ERI_Dyz_S_Dxy_S_C1000001_ac;
  abcd[446] = 4.0E0*I_ERI_D2z_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[456] = 4.0E0*I_ERI_Dxz_S_Dxz_S_C1000001_ac;
  abcd[457] = 4.0E0*I_ERI_Dyz_S_Dxz_S_C1000001_ac;
  abcd[458] = 4.0E0*I_ERI_D2z_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[435] = 4.0E0*I_ERI_Dxz_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[436] = 4.0E0*I_ERI_Dyz_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[437] = 4.0E0*I_ERI_D2z_Px_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2x_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[438] = 4.0E0*I_ERI_Dxz_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[439] = 4.0E0*I_ERI_Dyz_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[440] = 4.0E0*I_ERI_D2z_Py_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2x_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[441] = 4.0E0*I_ERI_Dxz_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;
  abcd[442] = 4.0E0*I_ERI_Dyz_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;
  abcd[443] = 4.0E0*I_ERI_D2z_Pz_D2x_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2x_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[447] = 4.0E0*I_ERI_Dxz_Px_Dxy_S_C1001001_ac;
  abcd[448] = 4.0E0*I_ERI_Dyz_Px_Dxy_S_C1001001_ac;
  abcd[449] = 4.0E0*I_ERI_D2z_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[450] = 4.0E0*I_ERI_Dxz_Py_Dxy_S_C1001001_ac;
  abcd[451] = 4.0E0*I_ERI_Dyz_Py_Dxy_S_C1001001_ac;
  abcd[452] = 4.0E0*I_ERI_D2z_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[453] = 4.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001001_ac;
  abcd[454] = 4.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001001_ac;
  abcd[455] = 4.0E0*I_ERI_D2z_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[459] = 4.0E0*I_ERI_Dxz_Px_Dxz_S_C1001001_ac;
  abcd[460] = 4.0E0*I_ERI_Dyz_Px_Dxz_S_C1001001_ac;
  abcd[461] = 4.0E0*I_ERI_D2z_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[462] = 4.0E0*I_ERI_Dxz_Py_Dxz_S_C1001001_ac;
  abcd[463] = 4.0E0*I_ERI_Dyz_Py_Dxz_S_C1001001_ac;
  abcd[464] = 4.0E0*I_ERI_D2z_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[465] = 4.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001001_ac;
  abcd[466] = 4.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001001_ac;
  abcd[467] = 4.0E0*I_ERI_D2z_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[468] = 4.0E0*I_ERI_Dxz_S_Dxy_S_C1000001_ac;
  abcd[469] = 4.0E0*I_ERI_Dyz_S_Dxy_S_C1000001_ac;
  abcd[470] = 4.0E0*I_ERI_D2z_S_Dxy_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxy_S_C1000001_c;
  abcd[480] = 4.0E0*I_ERI_Dxz_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;
  abcd[481] = 4.0E0*I_ERI_Dyz_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;
  abcd[482] = 4.0E0*I_ERI_D2z_S_D2y_S_C1000001_ac-2.0E0*1*I_ERI_D2z_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2y_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;
  abcd[492] = 4.0E0*I_ERI_Dxz_S_Dyz_S_C1000001_ac;
  abcd[493] = 4.0E0*I_ERI_Dyz_S_Dyz_S_C1000001_ac;
  abcd[494] = 4.0E0*I_ERI_D2z_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[471] = 4.0E0*I_ERI_Dxz_Px_Dxy_S_C1001001_ac;
  abcd[472] = 4.0E0*I_ERI_Dyz_Px_Dxy_S_C1001001_ac;
  abcd[473] = 4.0E0*I_ERI_D2z_Px_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxy_S_C1001001_c;
  abcd[474] = 4.0E0*I_ERI_Dxz_Py_Dxy_S_C1001001_ac;
  abcd[475] = 4.0E0*I_ERI_Dyz_Py_Dxy_S_C1001001_ac;
  abcd[476] = 4.0E0*I_ERI_D2z_Py_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxy_S_C1001001_c;
  abcd[477] = 4.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001001_ac;
  abcd[478] = 4.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001001_ac;
  abcd[479] = 4.0E0*I_ERI_D2z_Pz_Dxy_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxy_S_C1001001_c;
  abcd[483] = 4.0E0*I_ERI_Dxz_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[484] = 4.0E0*I_ERI_Dyz_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[485] = 4.0E0*I_ERI_D2z_Px_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2y_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[486] = 4.0E0*I_ERI_Dxz_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[487] = 4.0E0*I_ERI_Dyz_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[488] = 4.0E0*I_ERI_D2z_Py_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2y_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[489] = 4.0E0*I_ERI_Dxz_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;
  abcd[490] = 4.0E0*I_ERI_Dyz_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;
  abcd[491] = 4.0E0*I_ERI_D2z_Pz_D2y_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2y_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;
  abcd[495] = 4.0E0*I_ERI_Dxz_Px_Dyz_S_C1001001_ac;
  abcd[496] = 4.0E0*I_ERI_Dyz_Px_Dyz_S_C1001001_ac;
  abcd[497] = 4.0E0*I_ERI_D2z_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[498] = 4.0E0*I_ERI_Dxz_Py_Dyz_S_C1001001_ac;
  abcd[499] = 4.0E0*I_ERI_Dyz_Py_Dyz_S_C1001001_ac;
  abcd[500] = 4.0E0*I_ERI_D2z_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[501] = 4.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001001_ac;
  abcd[502] = 4.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001001_ac;
  abcd[503] = 4.0E0*I_ERI_D2z_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001_ac
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000001
   ************************************************************/
  abcd[504] = 4.0E0*I_ERI_Dxz_S_Dxz_S_C1000001_ac;
  abcd[505] = 4.0E0*I_ERI_Dyz_S_Dxz_S_C1000001_ac;
  abcd[506] = 4.0E0*I_ERI_D2z_S_Dxz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dxz_S_C1000001_c;
  abcd[516] = 4.0E0*I_ERI_Dxz_S_Dyz_S_C1000001_ac;
  abcd[517] = 4.0E0*I_ERI_Dyz_S_Dyz_S_C1000001_ac;
  abcd[518] = 4.0E0*I_ERI_D2z_S_Dyz_S_C1000001_ac-2.0E0*1*I_ERI_S_S_Dyz_S_C1000001_c;
  abcd[528] = 4.0E0*I_ERI_Dxz_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dxz_S_S_S_C1000001_a;
  abcd[529] = 4.0E0*I_ERI_Dyz_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_Dyz_S_S_S_C1000001_a;
  abcd[530] = 4.0E0*I_ERI_D2z_S_D2z_S_C1000001_ac-2.0E0*1*I_ERI_D2z_S_S_S_C1000001_a-2.0E0*1*I_ERI_S_S_D2z_S_C1000001_c+1*I_ERI_S_S_S_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001001_ac
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001001
   ************************************************************/
  abcd[507] = 4.0E0*I_ERI_Dxz_Px_Dxz_S_C1001001_ac;
  abcd[508] = 4.0E0*I_ERI_Dyz_Px_Dxz_S_C1001001_ac;
  abcd[509] = 4.0E0*I_ERI_D2z_Px_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dxz_S_C1001001_c;
  abcd[510] = 4.0E0*I_ERI_Dxz_Py_Dxz_S_C1001001_ac;
  abcd[511] = 4.0E0*I_ERI_Dyz_Py_Dxz_S_C1001001_ac;
  abcd[512] = 4.0E0*I_ERI_D2z_Py_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dxz_S_C1001001_c;
  abcd[513] = 4.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001001_ac;
  abcd[514] = 4.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001001_ac;
  abcd[515] = 4.0E0*I_ERI_D2z_Pz_Dxz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dxz_S_C1001001_c;
  abcd[519] = 4.0E0*I_ERI_Dxz_Px_Dyz_S_C1001001_ac;
  abcd[520] = 4.0E0*I_ERI_Dyz_Px_Dyz_S_C1001001_ac;
  abcd[521] = 4.0E0*I_ERI_D2z_Px_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Px_Dyz_S_C1001001_c;
  abcd[522] = 4.0E0*I_ERI_Dxz_Py_Dyz_S_C1001001_ac;
  abcd[523] = 4.0E0*I_ERI_Dyz_Py_Dyz_S_C1001001_ac;
  abcd[524] = 4.0E0*I_ERI_D2z_Py_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Py_Dyz_S_C1001001_c;
  abcd[525] = 4.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001001_ac;
  abcd[526] = 4.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001001_ac;
  abcd[527] = 4.0E0*I_ERI_D2z_Pz_Dyz_S_C1001001_ac-2.0E0*1*I_ERI_S_Pz_Dyz_S_C1001001_c;
  abcd[531] = 4.0E0*I_ERI_Dxz_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Px_S_S_C1001001_a;
  abcd[532] = 4.0E0*I_ERI_Dyz_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Px_S_S_C1001001_a;
  abcd[533] = 4.0E0*I_ERI_D2z_Px_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Px_S_S_C1001001_a-2.0E0*1*I_ERI_S_Px_D2z_S_C1001001_c+1*I_ERI_S_Px_S_S_C1001001;
  abcd[534] = 4.0E0*I_ERI_Dxz_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Py_S_S_C1001001_a;
  abcd[535] = 4.0E0*I_ERI_Dyz_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Py_S_S_C1001001_a;
  abcd[536] = 4.0E0*I_ERI_D2z_Py_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Py_S_S_C1001001_a-2.0E0*1*I_ERI_S_Py_D2z_S_C1001001_c+1*I_ERI_S_Py_S_S_C1001001;
  abcd[537] = 4.0E0*I_ERI_Dxz_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dxz_Pz_S_S_C1001001_a;
  abcd[538] = 4.0E0*I_ERI_Dyz_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_Dyz_Pz_S_S_C1001001_a;
  abcd[539] = 4.0E0*I_ERI_D2z_Pz_D2z_S_C1001001_ac-2.0E0*1*I_ERI_D2z_Pz_S_S_C1001001_a-2.0E0*1*I_ERI_S_Pz_D2z_S_C1001001_c+1*I_ERI_S_Pz_S_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[540] = 4.0E0*I_ERI_D2x_S_Px_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Px_C1000001_d;
  abcd[541] = 4.0E0*I_ERI_Dxy_S_Px_Px_C1000001_ad;
  abcd[542] = 4.0E0*I_ERI_Dxz_S_Px_Px_C1000001_ad;
  abcd[552] = 4.0E0*I_ERI_D2x_S_Py_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Px_C1000001_d;
  abcd[553] = 4.0E0*I_ERI_Dxy_S_Py_Px_C1000001_ad;
  abcd[554] = 4.0E0*I_ERI_Dxz_S_Py_Px_C1000001_ad;
  abcd[564] = 4.0E0*I_ERI_D2x_S_Pz_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Px_C1000001_d;
  abcd[565] = 4.0E0*I_ERI_Dxy_S_Pz_Px_C1000001_ad;
  abcd[566] = 4.0E0*I_ERI_Dxz_S_Pz_Px_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[543] = 4.0E0*I_ERI_D2x_Px_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Px_C1001001_d;
  abcd[544] = 4.0E0*I_ERI_Dxy_Px_Px_Px_C1001001_ad;
  abcd[545] = 4.0E0*I_ERI_Dxz_Px_Px_Px_C1001001_ad;
  abcd[546] = 4.0E0*I_ERI_D2x_Py_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Px_C1001001_d;
  abcd[547] = 4.0E0*I_ERI_Dxy_Py_Px_Px_C1001001_ad;
  abcd[548] = 4.0E0*I_ERI_Dxz_Py_Px_Px_C1001001_ad;
  abcd[549] = 4.0E0*I_ERI_D2x_Pz_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Px_C1001001_d;
  abcd[550] = 4.0E0*I_ERI_Dxy_Pz_Px_Px_C1001001_ad;
  abcd[551] = 4.0E0*I_ERI_Dxz_Pz_Px_Px_C1001001_ad;
  abcd[555] = 4.0E0*I_ERI_D2x_Px_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Px_C1001001_d;
  abcd[556] = 4.0E0*I_ERI_Dxy_Px_Py_Px_C1001001_ad;
  abcd[557] = 4.0E0*I_ERI_Dxz_Px_Py_Px_C1001001_ad;
  abcd[558] = 4.0E0*I_ERI_D2x_Py_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Px_C1001001_d;
  abcd[559] = 4.0E0*I_ERI_Dxy_Py_Py_Px_C1001001_ad;
  abcd[560] = 4.0E0*I_ERI_Dxz_Py_Py_Px_C1001001_ad;
  abcd[561] = 4.0E0*I_ERI_D2x_Pz_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Px_C1001001_d;
  abcd[562] = 4.0E0*I_ERI_Dxy_Pz_Py_Px_C1001001_ad;
  abcd[563] = 4.0E0*I_ERI_Dxz_Pz_Py_Px_C1001001_ad;
  abcd[567] = 4.0E0*I_ERI_D2x_Px_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Px_C1001001_d;
  abcd[568] = 4.0E0*I_ERI_Dxy_Px_Pz_Px_C1001001_ad;
  abcd[569] = 4.0E0*I_ERI_Dxz_Px_Pz_Px_C1001001_ad;
  abcd[570] = 4.0E0*I_ERI_D2x_Py_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Px_C1001001_d;
  abcd[571] = 4.0E0*I_ERI_Dxy_Py_Pz_Px_C1001001_ad;
  abcd[572] = 4.0E0*I_ERI_Dxz_Py_Pz_Px_C1001001_ad;
  abcd[573] = 4.0E0*I_ERI_D2x_Pz_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Px_C1001001_d;
  abcd[574] = 4.0E0*I_ERI_Dxy_Pz_Pz_Px_C1001001_ad;
  abcd[575] = 4.0E0*I_ERI_Dxz_Pz_Pz_Px_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[576] = 4.0E0*I_ERI_D2x_S_Px_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Py_C1000001_d;
  abcd[577] = 4.0E0*I_ERI_Dxy_S_Px_Py_C1000001_ad;
  abcd[578] = 4.0E0*I_ERI_Dxz_S_Px_Py_C1000001_ad;
  abcd[588] = 4.0E0*I_ERI_D2x_S_Py_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Py_C1000001_d;
  abcd[589] = 4.0E0*I_ERI_Dxy_S_Py_Py_C1000001_ad;
  abcd[590] = 4.0E0*I_ERI_Dxz_S_Py_Py_C1000001_ad;
  abcd[600] = 4.0E0*I_ERI_D2x_S_Pz_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Py_C1000001_d;
  abcd[601] = 4.0E0*I_ERI_Dxy_S_Pz_Py_C1000001_ad;
  abcd[602] = 4.0E0*I_ERI_Dxz_S_Pz_Py_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[579] = 4.0E0*I_ERI_D2x_Px_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Py_C1001001_d;
  abcd[580] = 4.0E0*I_ERI_Dxy_Px_Px_Py_C1001001_ad;
  abcd[581] = 4.0E0*I_ERI_Dxz_Px_Px_Py_C1001001_ad;
  abcd[582] = 4.0E0*I_ERI_D2x_Py_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Py_C1001001_d;
  abcd[583] = 4.0E0*I_ERI_Dxy_Py_Px_Py_C1001001_ad;
  abcd[584] = 4.0E0*I_ERI_Dxz_Py_Px_Py_C1001001_ad;
  abcd[585] = 4.0E0*I_ERI_D2x_Pz_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Py_C1001001_d;
  abcd[586] = 4.0E0*I_ERI_Dxy_Pz_Px_Py_C1001001_ad;
  abcd[587] = 4.0E0*I_ERI_Dxz_Pz_Px_Py_C1001001_ad;
  abcd[591] = 4.0E0*I_ERI_D2x_Px_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Py_C1001001_d;
  abcd[592] = 4.0E0*I_ERI_Dxy_Px_Py_Py_C1001001_ad;
  abcd[593] = 4.0E0*I_ERI_Dxz_Px_Py_Py_C1001001_ad;
  abcd[594] = 4.0E0*I_ERI_D2x_Py_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Py_C1001001_d;
  abcd[595] = 4.0E0*I_ERI_Dxy_Py_Py_Py_C1001001_ad;
  abcd[596] = 4.0E0*I_ERI_Dxz_Py_Py_Py_C1001001_ad;
  abcd[597] = 4.0E0*I_ERI_D2x_Pz_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Py_C1001001_d;
  abcd[598] = 4.0E0*I_ERI_Dxy_Pz_Py_Py_C1001001_ad;
  abcd[599] = 4.0E0*I_ERI_Dxz_Pz_Py_Py_C1001001_ad;
  abcd[603] = 4.0E0*I_ERI_D2x_Px_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Py_C1001001_d;
  abcd[604] = 4.0E0*I_ERI_Dxy_Px_Pz_Py_C1001001_ad;
  abcd[605] = 4.0E0*I_ERI_Dxz_Px_Pz_Py_C1001001_ad;
  abcd[606] = 4.0E0*I_ERI_D2x_Py_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Py_C1001001_d;
  abcd[607] = 4.0E0*I_ERI_Dxy_Py_Pz_Py_C1001001_ad;
  abcd[608] = 4.0E0*I_ERI_Dxz_Py_Pz_Py_C1001001_ad;
  abcd[609] = 4.0E0*I_ERI_D2x_Pz_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Py_C1001001_d;
  abcd[610] = 4.0E0*I_ERI_Dxy_Pz_Pz_Py_C1001001_ad;
  abcd[611] = 4.0E0*I_ERI_Dxz_Pz_Pz_Py_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[612] = 4.0E0*I_ERI_D2x_S_Px_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Pz_C1000001_d;
  abcd[613] = 4.0E0*I_ERI_Dxy_S_Px_Pz_C1000001_ad;
  abcd[614] = 4.0E0*I_ERI_Dxz_S_Px_Pz_C1000001_ad;
  abcd[624] = 4.0E0*I_ERI_D2x_S_Py_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Pz_C1000001_d;
  abcd[625] = 4.0E0*I_ERI_Dxy_S_Py_Pz_C1000001_ad;
  abcd[626] = 4.0E0*I_ERI_Dxz_S_Py_Pz_C1000001_ad;
  abcd[636] = 4.0E0*I_ERI_D2x_S_Pz_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Pz_C1000001_d;
  abcd[637] = 4.0E0*I_ERI_Dxy_S_Pz_Pz_C1000001_ad;
  abcd[638] = 4.0E0*I_ERI_Dxz_S_Pz_Pz_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[615] = 4.0E0*I_ERI_D2x_Px_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Pz_C1001001_d;
  abcd[616] = 4.0E0*I_ERI_Dxy_Px_Px_Pz_C1001001_ad;
  abcd[617] = 4.0E0*I_ERI_Dxz_Px_Px_Pz_C1001001_ad;
  abcd[618] = 4.0E0*I_ERI_D2x_Py_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Pz_C1001001_d;
  abcd[619] = 4.0E0*I_ERI_Dxy_Py_Px_Pz_C1001001_ad;
  abcd[620] = 4.0E0*I_ERI_Dxz_Py_Px_Pz_C1001001_ad;
  abcd[621] = 4.0E0*I_ERI_D2x_Pz_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Pz_C1001001_d;
  abcd[622] = 4.0E0*I_ERI_Dxy_Pz_Px_Pz_C1001001_ad;
  abcd[623] = 4.0E0*I_ERI_Dxz_Pz_Px_Pz_C1001001_ad;
  abcd[627] = 4.0E0*I_ERI_D2x_Px_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Pz_C1001001_d;
  abcd[628] = 4.0E0*I_ERI_Dxy_Px_Py_Pz_C1001001_ad;
  abcd[629] = 4.0E0*I_ERI_Dxz_Px_Py_Pz_C1001001_ad;
  abcd[630] = 4.0E0*I_ERI_D2x_Py_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Pz_C1001001_d;
  abcd[631] = 4.0E0*I_ERI_Dxy_Py_Py_Pz_C1001001_ad;
  abcd[632] = 4.0E0*I_ERI_Dxz_Py_Py_Pz_C1001001_ad;
  abcd[633] = 4.0E0*I_ERI_D2x_Pz_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Pz_C1001001_d;
  abcd[634] = 4.0E0*I_ERI_Dxy_Pz_Py_Pz_C1001001_ad;
  abcd[635] = 4.0E0*I_ERI_Dxz_Pz_Py_Pz_C1001001_ad;
  abcd[639] = 4.0E0*I_ERI_D2x_Px_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Pz_C1001001_d;
  abcd[640] = 4.0E0*I_ERI_Dxy_Px_Pz_Pz_C1001001_ad;
  abcd[641] = 4.0E0*I_ERI_Dxz_Px_Pz_Pz_C1001001_ad;
  abcd[642] = 4.0E0*I_ERI_D2x_Py_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Pz_C1001001_d;
  abcd[643] = 4.0E0*I_ERI_Dxy_Py_Pz_Pz_C1001001_ad;
  abcd[644] = 4.0E0*I_ERI_Dxz_Py_Pz_Pz_C1001001_ad;
  abcd[645] = 4.0E0*I_ERI_D2x_Pz_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Pz_C1001001_d;
  abcd[646] = 4.0E0*I_ERI_Dxy_Pz_Pz_Pz_C1001001_ad;
  abcd[647] = 4.0E0*I_ERI_Dxz_Pz_Pz_Pz_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[648] = 4.0E0*I_ERI_Dxy_S_Px_Px_C1000001_ad;
  abcd[649] = 4.0E0*I_ERI_D2y_S_Px_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Px_C1000001_d;
  abcd[650] = 4.0E0*I_ERI_Dyz_S_Px_Px_C1000001_ad;
  abcd[660] = 4.0E0*I_ERI_Dxy_S_Py_Px_C1000001_ad;
  abcd[661] = 4.0E0*I_ERI_D2y_S_Py_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Px_C1000001_d;
  abcd[662] = 4.0E0*I_ERI_Dyz_S_Py_Px_C1000001_ad;
  abcd[672] = 4.0E0*I_ERI_Dxy_S_Pz_Px_C1000001_ad;
  abcd[673] = 4.0E0*I_ERI_D2y_S_Pz_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Px_C1000001_d;
  abcd[674] = 4.0E0*I_ERI_Dyz_S_Pz_Px_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[651] = 4.0E0*I_ERI_Dxy_Px_Px_Px_C1001001_ad;
  abcd[652] = 4.0E0*I_ERI_D2y_Px_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Px_C1001001_d;
  abcd[653] = 4.0E0*I_ERI_Dyz_Px_Px_Px_C1001001_ad;
  abcd[654] = 4.0E0*I_ERI_Dxy_Py_Px_Px_C1001001_ad;
  abcd[655] = 4.0E0*I_ERI_D2y_Py_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Px_C1001001_d;
  abcd[656] = 4.0E0*I_ERI_Dyz_Py_Px_Px_C1001001_ad;
  abcd[657] = 4.0E0*I_ERI_Dxy_Pz_Px_Px_C1001001_ad;
  abcd[658] = 4.0E0*I_ERI_D2y_Pz_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Px_C1001001_d;
  abcd[659] = 4.0E0*I_ERI_Dyz_Pz_Px_Px_C1001001_ad;
  abcd[663] = 4.0E0*I_ERI_Dxy_Px_Py_Px_C1001001_ad;
  abcd[664] = 4.0E0*I_ERI_D2y_Px_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Px_C1001001_d;
  abcd[665] = 4.0E0*I_ERI_Dyz_Px_Py_Px_C1001001_ad;
  abcd[666] = 4.0E0*I_ERI_Dxy_Py_Py_Px_C1001001_ad;
  abcd[667] = 4.0E0*I_ERI_D2y_Py_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Px_C1001001_d;
  abcd[668] = 4.0E0*I_ERI_Dyz_Py_Py_Px_C1001001_ad;
  abcd[669] = 4.0E0*I_ERI_Dxy_Pz_Py_Px_C1001001_ad;
  abcd[670] = 4.0E0*I_ERI_D2y_Pz_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Px_C1001001_d;
  abcd[671] = 4.0E0*I_ERI_Dyz_Pz_Py_Px_C1001001_ad;
  abcd[675] = 4.0E0*I_ERI_Dxy_Px_Pz_Px_C1001001_ad;
  abcd[676] = 4.0E0*I_ERI_D2y_Px_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Px_C1001001_d;
  abcd[677] = 4.0E0*I_ERI_Dyz_Px_Pz_Px_C1001001_ad;
  abcd[678] = 4.0E0*I_ERI_Dxy_Py_Pz_Px_C1001001_ad;
  abcd[679] = 4.0E0*I_ERI_D2y_Py_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Px_C1001001_d;
  abcd[680] = 4.0E0*I_ERI_Dyz_Py_Pz_Px_C1001001_ad;
  abcd[681] = 4.0E0*I_ERI_Dxy_Pz_Pz_Px_C1001001_ad;
  abcd[682] = 4.0E0*I_ERI_D2y_Pz_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Px_C1001001_d;
  abcd[683] = 4.0E0*I_ERI_Dyz_Pz_Pz_Px_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[684] = 4.0E0*I_ERI_Dxy_S_Px_Py_C1000001_ad;
  abcd[685] = 4.0E0*I_ERI_D2y_S_Px_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Py_C1000001_d;
  abcd[686] = 4.0E0*I_ERI_Dyz_S_Px_Py_C1000001_ad;
  abcd[696] = 4.0E0*I_ERI_Dxy_S_Py_Py_C1000001_ad;
  abcd[697] = 4.0E0*I_ERI_D2y_S_Py_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Py_C1000001_d;
  abcd[698] = 4.0E0*I_ERI_Dyz_S_Py_Py_C1000001_ad;
  abcd[708] = 4.0E0*I_ERI_Dxy_S_Pz_Py_C1000001_ad;
  abcd[709] = 4.0E0*I_ERI_D2y_S_Pz_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Py_C1000001_d;
  abcd[710] = 4.0E0*I_ERI_Dyz_S_Pz_Py_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[687] = 4.0E0*I_ERI_Dxy_Px_Px_Py_C1001001_ad;
  abcd[688] = 4.0E0*I_ERI_D2y_Px_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Py_C1001001_d;
  abcd[689] = 4.0E0*I_ERI_Dyz_Px_Px_Py_C1001001_ad;
  abcd[690] = 4.0E0*I_ERI_Dxy_Py_Px_Py_C1001001_ad;
  abcd[691] = 4.0E0*I_ERI_D2y_Py_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Py_C1001001_d;
  abcd[692] = 4.0E0*I_ERI_Dyz_Py_Px_Py_C1001001_ad;
  abcd[693] = 4.0E0*I_ERI_Dxy_Pz_Px_Py_C1001001_ad;
  abcd[694] = 4.0E0*I_ERI_D2y_Pz_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Py_C1001001_d;
  abcd[695] = 4.0E0*I_ERI_Dyz_Pz_Px_Py_C1001001_ad;
  abcd[699] = 4.0E0*I_ERI_Dxy_Px_Py_Py_C1001001_ad;
  abcd[700] = 4.0E0*I_ERI_D2y_Px_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Py_C1001001_d;
  abcd[701] = 4.0E0*I_ERI_Dyz_Px_Py_Py_C1001001_ad;
  abcd[702] = 4.0E0*I_ERI_Dxy_Py_Py_Py_C1001001_ad;
  abcd[703] = 4.0E0*I_ERI_D2y_Py_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Py_C1001001_d;
  abcd[704] = 4.0E0*I_ERI_Dyz_Py_Py_Py_C1001001_ad;
  abcd[705] = 4.0E0*I_ERI_Dxy_Pz_Py_Py_C1001001_ad;
  abcd[706] = 4.0E0*I_ERI_D2y_Pz_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Py_C1001001_d;
  abcd[707] = 4.0E0*I_ERI_Dyz_Pz_Py_Py_C1001001_ad;
  abcd[711] = 4.0E0*I_ERI_Dxy_Px_Pz_Py_C1001001_ad;
  abcd[712] = 4.0E0*I_ERI_D2y_Px_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Py_C1001001_d;
  abcd[713] = 4.0E0*I_ERI_Dyz_Px_Pz_Py_C1001001_ad;
  abcd[714] = 4.0E0*I_ERI_Dxy_Py_Pz_Py_C1001001_ad;
  abcd[715] = 4.0E0*I_ERI_D2y_Py_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Py_C1001001_d;
  abcd[716] = 4.0E0*I_ERI_Dyz_Py_Pz_Py_C1001001_ad;
  abcd[717] = 4.0E0*I_ERI_Dxy_Pz_Pz_Py_C1001001_ad;
  abcd[718] = 4.0E0*I_ERI_D2y_Pz_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Py_C1001001_d;
  abcd[719] = 4.0E0*I_ERI_Dyz_Pz_Pz_Py_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[720] = 4.0E0*I_ERI_Dxy_S_Px_Pz_C1000001_ad;
  abcd[721] = 4.0E0*I_ERI_D2y_S_Px_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Pz_C1000001_d;
  abcd[722] = 4.0E0*I_ERI_Dyz_S_Px_Pz_C1000001_ad;
  abcd[732] = 4.0E0*I_ERI_Dxy_S_Py_Pz_C1000001_ad;
  abcd[733] = 4.0E0*I_ERI_D2y_S_Py_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Pz_C1000001_d;
  abcd[734] = 4.0E0*I_ERI_Dyz_S_Py_Pz_C1000001_ad;
  abcd[744] = 4.0E0*I_ERI_Dxy_S_Pz_Pz_C1000001_ad;
  abcd[745] = 4.0E0*I_ERI_D2y_S_Pz_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Pz_C1000001_d;
  abcd[746] = 4.0E0*I_ERI_Dyz_S_Pz_Pz_C1000001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[723] = 4.0E0*I_ERI_Dxy_Px_Px_Pz_C1001001_ad;
  abcd[724] = 4.0E0*I_ERI_D2y_Px_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Pz_C1001001_d;
  abcd[725] = 4.0E0*I_ERI_Dyz_Px_Px_Pz_C1001001_ad;
  abcd[726] = 4.0E0*I_ERI_Dxy_Py_Px_Pz_C1001001_ad;
  abcd[727] = 4.0E0*I_ERI_D2y_Py_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Pz_C1001001_d;
  abcd[728] = 4.0E0*I_ERI_Dyz_Py_Px_Pz_C1001001_ad;
  abcd[729] = 4.0E0*I_ERI_Dxy_Pz_Px_Pz_C1001001_ad;
  abcd[730] = 4.0E0*I_ERI_D2y_Pz_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Pz_C1001001_d;
  abcd[731] = 4.0E0*I_ERI_Dyz_Pz_Px_Pz_C1001001_ad;
  abcd[735] = 4.0E0*I_ERI_Dxy_Px_Py_Pz_C1001001_ad;
  abcd[736] = 4.0E0*I_ERI_D2y_Px_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Pz_C1001001_d;
  abcd[737] = 4.0E0*I_ERI_Dyz_Px_Py_Pz_C1001001_ad;
  abcd[738] = 4.0E0*I_ERI_Dxy_Py_Py_Pz_C1001001_ad;
  abcd[739] = 4.0E0*I_ERI_D2y_Py_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Pz_C1001001_d;
  abcd[740] = 4.0E0*I_ERI_Dyz_Py_Py_Pz_C1001001_ad;
  abcd[741] = 4.0E0*I_ERI_Dxy_Pz_Py_Pz_C1001001_ad;
  abcd[742] = 4.0E0*I_ERI_D2y_Pz_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Pz_C1001001_d;
  abcd[743] = 4.0E0*I_ERI_Dyz_Pz_Py_Pz_C1001001_ad;
  abcd[747] = 4.0E0*I_ERI_Dxy_Px_Pz_Pz_C1001001_ad;
  abcd[748] = 4.0E0*I_ERI_D2y_Px_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Pz_C1001001_d;
  abcd[749] = 4.0E0*I_ERI_Dyz_Px_Pz_Pz_C1001001_ad;
  abcd[750] = 4.0E0*I_ERI_Dxy_Py_Pz_Pz_C1001001_ad;
  abcd[751] = 4.0E0*I_ERI_D2y_Py_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Pz_C1001001_d;
  abcd[752] = 4.0E0*I_ERI_Dyz_Py_Pz_Pz_C1001001_ad;
  abcd[753] = 4.0E0*I_ERI_Dxy_Pz_Pz_Pz_C1001001_ad;
  abcd[754] = 4.0E0*I_ERI_D2y_Pz_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Pz_C1001001_d;
  abcd[755] = 4.0E0*I_ERI_Dyz_Pz_Pz_Pz_C1001001_ad;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[756] = 4.0E0*I_ERI_Dxz_S_Px_Px_C1000001_ad;
  abcd[757] = 4.0E0*I_ERI_Dyz_S_Px_Px_C1000001_ad;
  abcd[758] = 4.0E0*I_ERI_D2z_S_Px_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Px_C1000001_d;
  abcd[768] = 4.0E0*I_ERI_Dxz_S_Py_Px_C1000001_ad;
  abcd[769] = 4.0E0*I_ERI_Dyz_S_Py_Px_C1000001_ad;
  abcd[770] = 4.0E0*I_ERI_D2z_S_Py_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Px_C1000001_d;
  abcd[780] = 4.0E0*I_ERI_Dxz_S_Pz_Px_C1000001_ad;
  abcd[781] = 4.0E0*I_ERI_Dyz_S_Pz_Px_C1000001_ad;
  abcd[782] = 4.0E0*I_ERI_D2z_S_Pz_Px_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Px_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[759] = 4.0E0*I_ERI_Dxz_Px_Px_Px_C1001001_ad;
  abcd[760] = 4.0E0*I_ERI_Dyz_Px_Px_Px_C1001001_ad;
  abcd[761] = 4.0E0*I_ERI_D2z_Px_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Px_C1001001_d;
  abcd[762] = 4.0E0*I_ERI_Dxz_Py_Px_Px_C1001001_ad;
  abcd[763] = 4.0E0*I_ERI_Dyz_Py_Px_Px_C1001001_ad;
  abcd[764] = 4.0E0*I_ERI_D2z_Py_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Px_C1001001_d;
  abcd[765] = 4.0E0*I_ERI_Dxz_Pz_Px_Px_C1001001_ad;
  abcd[766] = 4.0E0*I_ERI_Dyz_Pz_Px_Px_C1001001_ad;
  abcd[767] = 4.0E0*I_ERI_D2z_Pz_Px_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Px_C1001001_d;
  abcd[771] = 4.0E0*I_ERI_Dxz_Px_Py_Px_C1001001_ad;
  abcd[772] = 4.0E0*I_ERI_Dyz_Px_Py_Px_C1001001_ad;
  abcd[773] = 4.0E0*I_ERI_D2z_Px_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Px_C1001001_d;
  abcd[774] = 4.0E0*I_ERI_Dxz_Py_Py_Px_C1001001_ad;
  abcd[775] = 4.0E0*I_ERI_Dyz_Py_Py_Px_C1001001_ad;
  abcd[776] = 4.0E0*I_ERI_D2z_Py_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Px_C1001001_d;
  abcd[777] = 4.0E0*I_ERI_Dxz_Pz_Py_Px_C1001001_ad;
  abcd[778] = 4.0E0*I_ERI_Dyz_Pz_Py_Px_C1001001_ad;
  abcd[779] = 4.0E0*I_ERI_D2z_Pz_Py_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Px_C1001001_d;
  abcd[783] = 4.0E0*I_ERI_Dxz_Px_Pz_Px_C1001001_ad;
  abcd[784] = 4.0E0*I_ERI_Dyz_Px_Pz_Px_C1001001_ad;
  abcd[785] = 4.0E0*I_ERI_D2z_Px_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Px_C1001001_d;
  abcd[786] = 4.0E0*I_ERI_Dxz_Py_Pz_Px_C1001001_ad;
  abcd[787] = 4.0E0*I_ERI_Dyz_Py_Pz_Px_C1001001_ad;
  abcd[788] = 4.0E0*I_ERI_D2z_Py_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Px_C1001001_d;
  abcd[789] = 4.0E0*I_ERI_Dxz_Pz_Pz_Px_C1001001_ad;
  abcd[790] = 4.0E0*I_ERI_Dyz_Pz_Pz_Px_C1001001_ad;
  abcd[791] = 4.0E0*I_ERI_D2z_Pz_Pz_Px_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Px_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[792] = 4.0E0*I_ERI_Dxz_S_Px_Py_C1000001_ad;
  abcd[793] = 4.0E0*I_ERI_Dyz_S_Px_Py_C1000001_ad;
  abcd[794] = 4.0E0*I_ERI_D2z_S_Px_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Py_C1000001_d;
  abcd[804] = 4.0E0*I_ERI_Dxz_S_Py_Py_C1000001_ad;
  abcd[805] = 4.0E0*I_ERI_Dyz_S_Py_Py_C1000001_ad;
  abcd[806] = 4.0E0*I_ERI_D2z_S_Py_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Py_C1000001_d;
  abcd[816] = 4.0E0*I_ERI_Dxz_S_Pz_Py_C1000001_ad;
  abcd[817] = 4.0E0*I_ERI_Dyz_S_Pz_Py_C1000001_ad;
  abcd[818] = 4.0E0*I_ERI_D2z_S_Pz_Py_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Py_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[795] = 4.0E0*I_ERI_Dxz_Px_Px_Py_C1001001_ad;
  abcd[796] = 4.0E0*I_ERI_Dyz_Px_Px_Py_C1001001_ad;
  abcd[797] = 4.0E0*I_ERI_D2z_Px_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Py_C1001001_d;
  abcd[798] = 4.0E0*I_ERI_Dxz_Py_Px_Py_C1001001_ad;
  abcd[799] = 4.0E0*I_ERI_Dyz_Py_Px_Py_C1001001_ad;
  abcd[800] = 4.0E0*I_ERI_D2z_Py_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Py_C1001001_d;
  abcd[801] = 4.0E0*I_ERI_Dxz_Pz_Px_Py_C1001001_ad;
  abcd[802] = 4.0E0*I_ERI_Dyz_Pz_Px_Py_C1001001_ad;
  abcd[803] = 4.0E0*I_ERI_D2z_Pz_Px_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Py_C1001001_d;
  abcd[807] = 4.0E0*I_ERI_Dxz_Px_Py_Py_C1001001_ad;
  abcd[808] = 4.0E0*I_ERI_Dyz_Px_Py_Py_C1001001_ad;
  abcd[809] = 4.0E0*I_ERI_D2z_Px_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Py_C1001001_d;
  abcd[810] = 4.0E0*I_ERI_Dxz_Py_Py_Py_C1001001_ad;
  abcd[811] = 4.0E0*I_ERI_Dyz_Py_Py_Py_C1001001_ad;
  abcd[812] = 4.0E0*I_ERI_D2z_Py_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Py_C1001001_d;
  abcd[813] = 4.0E0*I_ERI_Dxz_Pz_Py_Py_C1001001_ad;
  abcd[814] = 4.0E0*I_ERI_Dyz_Pz_Py_Py_C1001001_ad;
  abcd[815] = 4.0E0*I_ERI_D2z_Pz_Py_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Py_C1001001_d;
  abcd[819] = 4.0E0*I_ERI_Dxz_Px_Pz_Py_C1001001_ad;
  abcd[820] = 4.0E0*I_ERI_Dyz_Px_Pz_Py_C1001001_ad;
  abcd[821] = 4.0E0*I_ERI_D2z_Px_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Py_C1001001_d;
  abcd[822] = 4.0E0*I_ERI_Dxz_Py_Pz_Py_C1001001_ad;
  abcd[823] = 4.0E0*I_ERI_Dyz_Py_Pz_Py_C1001001_ad;
  abcd[824] = 4.0E0*I_ERI_D2z_Py_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Py_C1001001_d;
  abcd[825] = 4.0E0*I_ERI_Dxz_Pz_Pz_Py_C1001001_ad;
  abcd[826] = 4.0E0*I_ERI_Dyz_Pz_Pz_Py_C1001001_ad;
  abcd[827] = 4.0E0*I_ERI_D2z_Pz_Pz_Py_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Py_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001_ad
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000001_d
   ************************************************************/
  abcd[828] = 4.0E0*I_ERI_Dxz_S_Px_Pz_C1000001_ad;
  abcd[829] = 4.0E0*I_ERI_Dyz_S_Px_Pz_C1000001_ad;
  abcd[830] = 4.0E0*I_ERI_D2z_S_Px_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Px_Pz_C1000001_d;
  abcd[840] = 4.0E0*I_ERI_Dxz_S_Py_Pz_C1000001_ad;
  abcd[841] = 4.0E0*I_ERI_Dyz_S_Py_Pz_C1000001_ad;
  abcd[842] = 4.0E0*I_ERI_D2z_S_Py_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Py_Pz_C1000001_d;
  abcd[852] = 4.0E0*I_ERI_Dxz_S_Pz_Pz_C1000001_ad;
  abcd[853] = 4.0E0*I_ERI_Dyz_S_Pz_Pz_C1000001_ad;
  abcd[854] = 4.0E0*I_ERI_D2z_S_Pz_Pz_C1000001_ad-2.0E0*1*I_ERI_S_S_Pz_Pz_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001_ad
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001_d
   ************************************************************/
  abcd[831] = 4.0E0*I_ERI_Dxz_Px_Px_Pz_C1001001_ad;
  abcd[832] = 4.0E0*I_ERI_Dyz_Px_Px_Pz_C1001001_ad;
  abcd[833] = 4.0E0*I_ERI_D2z_Px_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Px_Pz_C1001001_d;
  abcd[834] = 4.0E0*I_ERI_Dxz_Py_Px_Pz_C1001001_ad;
  abcd[835] = 4.0E0*I_ERI_Dyz_Py_Px_Pz_C1001001_ad;
  abcd[836] = 4.0E0*I_ERI_D2z_Py_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Px_Pz_C1001001_d;
  abcd[837] = 4.0E0*I_ERI_Dxz_Pz_Px_Pz_C1001001_ad;
  abcd[838] = 4.0E0*I_ERI_Dyz_Pz_Px_Pz_C1001001_ad;
  abcd[839] = 4.0E0*I_ERI_D2z_Pz_Px_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Px_Pz_C1001001_d;
  abcd[843] = 4.0E0*I_ERI_Dxz_Px_Py_Pz_C1001001_ad;
  abcd[844] = 4.0E0*I_ERI_Dyz_Px_Py_Pz_C1001001_ad;
  abcd[845] = 4.0E0*I_ERI_D2z_Px_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Py_Pz_C1001001_d;
  abcd[846] = 4.0E0*I_ERI_Dxz_Py_Py_Pz_C1001001_ad;
  abcd[847] = 4.0E0*I_ERI_Dyz_Py_Py_Pz_C1001001_ad;
  abcd[848] = 4.0E0*I_ERI_D2z_Py_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Py_Pz_C1001001_d;
  abcd[849] = 4.0E0*I_ERI_Dxz_Pz_Py_Pz_C1001001_ad;
  abcd[850] = 4.0E0*I_ERI_Dyz_Pz_Py_Pz_C1001001_ad;
  abcd[851] = 4.0E0*I_ERI_D2z_Pz_Py_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Py_Pz_C1001001_d;
  abcd[855] = 4.0E0*I_ERI_Dxz_Px_Pz_Pz_C1001001_ad;
  abcd[856] = 4.0E0*I_ERI_Dyz_Px_Pz_Pz_C1001001_ad;
  abcd[857] = 4.0E0*I_ERI_D2z_Px_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Px_Pz_Pz_C1001001_d;
  abcd[858] = 4.0E0*I_ERI_Dxz_Py_Pz_Pz_C1001001_ad;
  abcd[859] = 4.0E0*I_ERI_Dyz_Py_Pz_Pz_C1001001_ad;
  abcd[860] = 4.0E0*I_ERI_D2z_Py_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Py_Pz_Pz_C1001001_d;
  abcd[861] = 4.0E0*I_ERI_Dxz_Pz_Pz_Pz_C1001001_ad;
  abcd[862] = 4.0E0*I_ERI_Dyz_Pz_Pz_Pz_C1001001_ad;
  abcd[863] = 4.0E0*I_ERI_D2z_Pz_Pz_Pz_C1001001_ad-2.0E0*1*I_ERI_S_Pz_Pz_Pz_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[864] = 4.0E0*I_ERI_Px_S_F3x_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_c-2.0E0*2*I_ERI_Px_S_Px_S_C1000001_c;
  abcd[865] = 4.0E0*I_ERI_Py_S_F3x_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_c-2.0E0*2*I_ERI_Py_S_Px_S_C1000001_c;
  abcd[866] = 4.0E0*I_ERI_Pz_S_F3x_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_c-2.0E0*2*I_ERI_Pz_S_Px_S_C1000001_c;
  abcd[876] = 4.0E0*I_ERI_Px_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_c;
  abcd[877] = 4.0E0*I_ERI_Py_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_c;
  abcd[878] = 4.0E0*I_ERI_Pz_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_c;
  abcd[888] = 4.0E0*I_ERI_Px_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_c;
  abcd[889] = 4.0E0*I_ERI_Py_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_c;
  abcd[890] = 4.0E0*I_ERI_Pz_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[867] = 4.0E0*I_ERI_Px_Px_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_c-2.0E0*2*I_ERI_Px_Px_Px_S_C1001001_c;
  abcd[868] = 4.0E0*I_ERI_Py_Px_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_c-2.0E0*2*I_ERI_Py_Px_Px_S_C1001001_c;
  abcd[869] = 4.0E0*I_ERI_Pz_Px_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_c-2.0E0*2*I_ERI_Pz_Px_Px_S_C1001001_c;
  abcd[870] = 4.0E0*I_ERI_Px_Py_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_c-2.0E0*2*I_ERI_Px_Py_Px_S_C1001001_c;
  abcd[871] = 4.0E0*I_ERI_Py_Py_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_c-2.0E0*2*I_ERI_Py_Py_Px_S_C1001001_c;
  abcd[872] = 4.0E0*I_ERI_Pz_Py_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_c-2.0E0*2*I_ERI_Pz_Py_Px_S_C1001001_c;
  abcd[873] = 4.0E0*I_ERI_Px_Pz_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_c-2.0E0*2*I_ERI_Px_Pz_Px_S_C1001001_c;
  abcd[874] = 4.0E0*I_ERI_Py_Pz_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_c-2.0E0*2*I_ERI_Py_Pz_Px_S_C1001001_c;
  abcd[875] = 4.0E0*I_ERI_Pz_Pz_F3x_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_c-2.0E0*2*I_ERI_Pz_Pz_Px_S_C1001001_c;
  abcd[879] = 4.0E0*I_ERI_Px_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_c;
  abcd[880] = 4.0E0*I_ERI_Py_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_c;
  abcd[881] = 4.0E0*I_ERI_Pz_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_c;
  abcd[882] = 4.0E0*I_ERI_Px_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_c;
  abcd[883] = 4.0E0*I_ERI_Py_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_c;
  abcd[884] = 4.0E0*I_ERI_Pz_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_c;
  abcd[885] = 4.0E0*I_ERI_Px_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_c;
  abcd[886] = 4.0E0*I_ERI_Py_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_c;
  abcd[887] = 4.0E0*I_ERI_Pz_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_c;
  abcd[891] = 4.0E0*I_ERI_Px_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_c;
  abcd[892] = 4.0E0*I_ERI_Py_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_c;
  abcd[893] = 4.0E0*I_ERI_Pz_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_c;
  abcd[894] = 4.0E0*I_ERI_Px_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_c;
  abcd[895] = 4.0E0*I_ERI_Py_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_c;
  abcd[896] = 4.0E0*I_ERI_Pz_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_c;
  abcd[897] = 4.0E0*I_ERI_Px_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_c;
  abcd[898] = 4.0E0*I_ERI_Py_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_c;
  abcd[899] = 4.0E0*I_ERI_Pz_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[900] = 4.0E0*I_ERI_Px_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_c;
  abcd[901] = 4.0E0*I_ERI_Py_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_c;
  abcd[902] = 4.0E0*I_ERI_Pz_S_F2xy_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_c;
  abcd[912] = 4.0E0*I_ERI_Px_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_c;
  abcd[913] = 4.0E0*I_ERI_Py_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_c;
  abcd[914] = 4.0E0*I_ERI_Pz_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_c;
  abcd[924] = 4.0E0*I_ERI_Px_S_Fxyz_S_C1000001_cc;
  abcd[925] = 4.0E0*I_ERI_Py_S_Fxyz_S_C1000001_cc;
  abcd[926] = 4.0E0*I_ERI_Pz_S_Fxyz_S_C1000001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[903] = 4.0E0*I_ERI_Px_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_c;
  abcd[904] = 4.0E0*I_ERI_Py_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_c;
  abcd[905] = 4.0E0*I_ERI_Pz_Px_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_c;
  abcd[906] = 4.0E0*I_ERI_Px_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_c;
  abcd[907] = 4.0E0*I_ERI_Py_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_c;
  abcd[908] = 4.0E0*I_ERI_Pz_Py_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_c;
  abcd[909] = 4.0E0*I_ERI_Px_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_c;
  abcd[910] = 4.0E0*I_ERI_Py_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_c;
  abcd[911] = 4.0E0*I_ERI_Pz_Pz_F2xy_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_c;
  abcd[915] = 4.0E0*I_ERI_Px_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_c;
  abcd[916] = 4.0E0*I_ERI_Py_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_c;
  abcd[917] = 4.0E0*I_ERI_Pz_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_c;
  abcd[918] = 4.0E0*I_ERI_Px_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_c;
  abcd[919] = 4.0E0*I_ERI_Py_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_c;
  abcd[920] = 4.0E0*I_ERI_Pz_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_c;
  abcd[921] = 4.0E0*I_ERI_Px_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_c;
  abcd[922] = 4.0E0*I_ERI_Py_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_c;
  abcd[923] = 4.0E0*I_ERI_Pz_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_c;
  abcd[927] = 4.0E0*I_ERI_Px_Px_Fxyz_S_C1001001_cc;
  abcd[928] = 4.0E0*I_ERI_Py_Px_Fxyz_S_C1001001_cc;
  abcd[929] = 4.0E0*I_ERI_Pz_Px_Fxyz_S_C1001001_cc;
  abcd[930] = 4.0E0*I_ERI_Px_Py_Fxyz_S_C1001001_cc;
  abcd[931] = 4.0E0*I_ERI_Py_Py_Fxyz_S_C1001001_cc;
  abcd[932] = 4.0E0*I_ERI_Pz_Py_Fxyz_S_C1001001_cc;
  abcd[933] = 4.0E0*I_ERI_Px_Pz_Fxyz_S_C1001001_cc;
  abcd[934] = 4.0E0*I_ERI_Py_Pz_Fxyz_S_C1001001_cc;
  abcd[935] = 4.0E0*I_ERI_Pz_Pz_Fxyz_S_C1001001_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[936] = 4.0E0*I_ERI_Px_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_c;
  abcd[937] = 4.0E0*I_ERI_Py_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_c;
  abcd[938] = 4.0E0*I_ERI_Pz_S_F2xz_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_c;
  abcd[948] = 4.0E0*I_ERI_Px_S_Fxyz_S_C1000001_cc;
  abcd[949] = 4.0E0*I_ERI_Py_S_Fxyz_S_C1000001_cc;
  abcd[950] = 4.0E0*I_ERI_Pz_S_Fxyz_S_C1000001_cc;
  abcd[960] = 4.0E0*I_ERI_Px_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_c;
  abcd[961] = 4.0E0*I_ERI_Py_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_c;
  abcd[962] = 4.0E0*I_ERI_Pz_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[939] = 4.0E0*I_ERI_Px_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_c;
  abcd[940] = 4.0E0*I_ERI_Py_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_c;
  abcd[941] = 4.0E0*I_ERI_Pz_Px_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_c;
  abcd[942] = 4.0E0*I_ERI_Px_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_c;
  abcd[943] = 4.0E0*I_ERI_Py_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_c;
  abcd[944] = 4.0E0*I_ERI_Pz_Py_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_c;
  abcd[945] = 4.0E0*I_ERI_Px_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_c;
  abcd[946] = 4.0E0*I_ERI_Py_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_c;
  abcd[947] = 4.0E0*I_ERI_Pz_Pz_F2xz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_c;
  abcd[951] = 4.0E0*I_ERI_Px_Px_Fxyz_S_C1001001_cc;
  abcd[952] = 4.0E0*I_ERI_Py_Px_Fxyz_S_C1001001_cc;
  abcd[953] = 4.0E0*I_ERI_Pz_Px_Fxyz_S_C1001001_cc;
  abcd[954] = 4.0E0*I_ERI_Px_Py_Fxyz_S_C1001001_cc;
  abcd[955] = 4.0E0*I_ERI_Py_Py_Fxyz_S_C1001001_cc;
  abcd[956] = 4.0E0*I_ERI_Pz_Py_Fxyz_S_C1001001_cc;
  abcd[957] = 4.0E0*I_ERI_Px_Pz_Fxyz_S_C1001001_cc;
  abcd[958] = 4.0E0*I_ERI_Py_Pz_Fxyz_S_C1001001_cc;
  abcd[959] = 4.0E0*I_ERI_Pz_Pz_Fxyz_S_C1001001_cc;
  abcd[963] = 4.0E0*I_ERI_Px_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_c;
  abcd[964] = 4.0E0*I_ERI_Py_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_c;
  abcd[965] = 4.0E0*I_ERI_Pz_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_c;
  abcd[966] = 4.0E0*I_ERI_Px_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_c;
  abcd[967] = 4.0E0*I_ERI_Py_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_c;
  abcd[968] = 4.0E0*I_ERI_Pz_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_c;
  abcd[969] = 4.0E0*I_ERI_Px_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_c;
  abcd[970] = 4.0E0*I_ERI_Py_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_c;
  abcd[971] = 4.0E0*I_ERI_Pz_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[972] = 4.0E0*I_ERI_Px_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_c;
  abcd[973] = 4.0E0*I_ERI_Py_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_c;
  abcd[974] = 4.0E0*I_ERI_Pz_S_Fx2y_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_c;
  abcd[984] = 4.0E0*I_ERI_Px_S_F3y_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_c-2.0E0*2*I_ERI_Px_S_Py_S_C1000001_c;
  abcd[985] = 4.0E0*I_ERI_Py_S_F3y_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_c-2.0E0*2*I_ERI_Py_S_Py_S_C1000001_c;
  abcd[986] = 4.0E0*I_ERI_Pz_S_F3y_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_c-2.0E0*2*I_ERI_Pz_S_Py_S_C1000001_c;
  abcd[996] = 4.0E0*I_ERI_Px_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_c;
  abcd[997] = 4.0E0*I_ERI_Py_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_c;
  abcd[998] = 4.0E0*I_ERI_Pz_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[975] = 4.0E0*I_ERI_Px_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_c;
  abcd[976] = 4.0E0*I_ERI_Py_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_c;
  abcd[977] = 4.0E0*I_ERI_Pz_Px_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_c;
  abcd[978] = 4.0E0*I_ERI_Px_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_c;
  abcd[979] = 4.0E0*I_ERI_Py_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_c;
  abcd[980] = 4.0E0*I_ERI_Pz_Py_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_c;
  abcd[981] = 4.0E0*I_ERI_Px_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_c;
  abcd[982] = 4.0E0*I_ERI_Py_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_c;
  abcd[983] = 4.0E0*I_ERI_Pz_Pz_Fx2y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_c;
  abcd[987] = 4.0E0*I_ERI_Px_Px_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_c-2.0E0*2*I_ERI_Px_Px_Py_S_C1001001_c;
  abcd[988] = 4.0E0*I_ERI_Py_Px_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_c-2.0E0*2*I_ERI_Py_Px_Py_S_C1001001_c;
  abcd[989] = 4.0E0*I_ERI_Pz_Px_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_c-2.0E0*2*I_ERI_Pz_Px_Py_S_C1001001_c;
  abcd[990] = 4.0E0*I_ERI_Px_Py_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_c-2.0E0*2*I_ERI_Px_Py_Py_S_C1001001_c;
  abcd[991] = 4.0E0*I_ERI_Py_Py_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_c-2.0E0*2*I_ERI_Py_Py_Py_S_C1001001_c;
  abcd[992] = 4.0E0*I_ERI_Pz_Py_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_c-2.0E0*2*I_ERI_Pz_Py_Py_S_C1001001_c;
  abcd[993] = 4.0E0*I_ERI_Px_Pz_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_c-2.0E0*2*I_ERI_Px_Pz_Py_S_C1001001_c;
  abcd[994] = 4.0E0*I_ERI_Py_Pz_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_c-2.0E0*2*I_ERI_Py_Pz_Py_S_C1001001_c;
  abcd[995] = 4.0E0*I_ERI_Pz_Pz_F3y_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_c-2.0E0*2*I_ERI_Pz_Pz_Py_S_C1001001_c;
  abcd[999] = 4.0E0*I_ERI_Px_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_c;
  abcd[1000] = 4.0E0*I_ERI_Py_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_c;
  abcd[1001] = 4.0E0*I_ERI_Pz_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_c;
  abcd[1002] = 4.0E0*I_ERI_Px_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_c;
  abcd[1003] = 4.0E0*I_ERI_Py_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_c;
  abcd[1004] = 4.0E0*I_ERI_Pz_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_c;
  abcd[1005] = 4.0E0*I_ERI_Px_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_c;
  abcd[1006] = 4.0E0*I_ERI_Py_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_c;
  abcd[1007] = 4.0E0*I_ERI_Pz_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[1008] = 4.0E0*I_ERI_Px_S_Fxyz_S_C1000001_cc;
  abcd[1009] = 4.0E0*I_ERI_Py_S_Fxyz_S_C1000001_cc;
  abcd[1010] = 4.0E0*I_ERI_Pz_S_Fxyz_S_C1000001_cc;
  abcd[1020] = 4.0E0*I_ERI_Px_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_c;
  abcd[1021] = 4.0E0*I_ERI_Py_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_c;
  abcd[1022] = 4.0E0*I_ERI_Pz_S_F2yz_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_c;
  abcd[1032] = 4.0E0*I_ERI_Px_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_c;
  abcd[1033] = 4.0E0*I_ERI_Py_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_c;
  abcd[1034] = 4.0E0*I_ERI_Pz_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[1011] = 4.0E0*I_ERI_Px_Px_Fxyz_S_C1001001_cc;
  abcd[1012] = 4.0E0*I_ERI_Py_Px_Fxyz_S_C1001001_cc;
  abcd[1013] = 4.0E0*I_ERI_Pz_Px_Fxyz_S_C1001001_cc;
  abcd[1014] = 4.0E0*I_ERI_Px_Py_Fxyz_S_C1001001_cc;
  abcd[1015] = 4.0E0*I_ERI_Py_Py_Fxyz_S_C1001001_cc;
  abcd[1016] = 4.0E0*I_ERI_Pz_Py_Fxyz_S_C1001001_cc;
  abcd[1017] = 4.0E0*I_ERI_Px_Pz_Fxyz_S_C1001001_cc;
  abcd[1018] = 4.0E0*I_ERI_Py_Pz_Fxyz_S_C1001001_cc;
  abcd[1019] = 4.0E0*I_ERI_Pz_Pz_Fxyz_S_C1001001_cc;
  abcd[1023] = 4.0E0*I_ERI_Px_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_c;
  abcd[1024] = 4.0E0*I_ERI_Py_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_c;
  abcd[1025] = 4.0E0*I_ERI_Pz_Px_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_c;
  abcd[1026] = 4.0E0*I_ERI_Px_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_c;
  abcd[1027] = 4.0E0*I_ERI_Py_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_c;
  abcd[1028] = 4.0E0*I_ERI_Pz_Py_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_c;
  abcd[1029] = 4.0E0*I_ERI_Px_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_c;
  abcd[1030] = 4.0E0*I_ERI_Py_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_c;
  abcd[1031] = 4.0E0*I_ERI_Pz_Pz_F2yz_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_c;
  abcd[1035] = 4.0E0*I_ERI_Px_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_c;
  abcd[1036] = 4.0E0*I_ERI_Py_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_c;
  abcd[1037] = 4.0E0*I_ERI_Pz_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_c;
  abcd[1038] = 4.0E0*I_ERI_Px_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_c;
  abcd[1039] = 4.0E0*I_ERI_Py_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_c;
  abcd[1040] = 4.0E0*I_ERI_Pz_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_c;
  abcd[1041] = 4.0E0*I_ERI_Px_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_c;
  abcd[1042] = 4.0E0*I_ERI_Py_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_c;
  abcd[1043] = 4.0E0*I_ERI_Pz_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1000001_cc
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_c
   ************************************************************/
  abcd[1044] = 4.0E0*I_ERI_Px_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_c;
  abcd[1045] = 4.0E0*I_ERI_Py_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_c;
  abcd[1046] = 4.0E0*I_ERI_Pz_S_Fx2z_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_c;
  abcd[1056] = 4.0E0*I_ERI_Px_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_c;
  abcd[1057] = 4.0E0*I_ERI_Py_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_c;
  abcd[1058] = 4.0E0*I_ERI_Pz_S_Fy2z_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_c;
  abcd[1068] = 4.0E0*I_ERI_Px_S_F3z_S_C1000001_cc-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_c-2.0E0*2*I_ERI_Px_S_Pz_S_C1000001_c;
  abcd[1069] = 4.0E0*I_ERI_Py_S_F3z_S_C1000001_cc-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_c-2.0E0*2*I_ERI_Py_S_Pz_S_C1000001_c;
  abcd[1070] = 4.0E0*I_ERI_Pz_S_F3z_S_C1000001_cc-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_c-2.0E0*2*I_ERI_Pz_S_Pz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_F_S_C1001001_cc
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_c
   ************************************************************/
  abcd[1047] = 4.0E0*I_ERI_Px_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_c;
  abcd[1048] = 4.0E0*I_ERI_Py_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_c;
  abcd[1049] = 4.0E0*I_ERI_Pz_Px_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_c;
  abcd[1050] = 4.0E0*I_ERI_Px_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_c;
  abcd[1051] = 4.0E0*I_ERI_Py_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_c;
  abcd[1052] = 4.0E0*I_ERI_Pz_Py_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_c;
  abcd[1053] = 4.0E0*I_ERI_Px_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_c;
  abcd[1054] = 4.0E0*I_ERI_Py_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_c;
  abcd[1055] = 4.0E0*I_ERI_Pz_Pz_Fx2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_c;
  abcd[1059] = 4.0E0*I_ERI_Px_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_c;
  abcd[1060] = 4.0E0*I_ERI_Py_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_c;
  abcd[1061] = 4.0E0*I_ERI_Pz_Px_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_c;
  abcd[1062] = 4.0E0*I_ERI_Px_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_c;
  abcd[1063] = 4.0E0*I_ERI_Py_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_c;
  abcd[1064] = 4.0E0*I_ERI_Pz_Py_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_c;
  abcd[1065] = 4.0E0*I_ERI_Px_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_c;
  abcd[1066] = 4.0E0*I_ERI_Py_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_c;
  abcd[1067] = 4.0E0*I_ERI_Pz_Pz_Fy2z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_c;
  abcd[1071] = 4.0E0*I_ERI_Px_Px_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_c-2.0E0*2*I_ERI_Px_Px_Pz_S_C1001001_c;
  abcd[1072] = 4.0E0*I_ERI_Py_Px_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_c-2.0E0*2*I_ERI_Py_Px_Pz_S_C1001001_c;
  abcd[1073] = 4.0E0*I_ERI_Pz_Px_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_c-2.0E0*2*I_ERI_Pz_Px_Pz_S_C1001001_c;
  abcd[1074] = 4.0E0*I_ERI_Px_Py_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_c-2.0E0*2*I_ERI_Px_Py_Pz_S_C1001001_c;
  abcd[1075] = 4.0E0*I_ERI_Py_Py_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_c-2.0E0*2*I_ERI_Py_Py_Pz_S_C1001001_c;
  abcd[1076] = 4.0E0*I_ERI_Pz_Py_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_c-2.0E0*2*I_ERI_Pz_Py_Pz_S_C1001001_c;
  abcd[1077] = 4.0E0*I_ERI_Px_Pz_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_c-2.0E0*2*I_ERI_Px_Pz_Pz_S_C1001001_c;
  abcd[1078] = 4.0E0*I_ERI_Py_Pz_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_c-2.0E0*2*I_ERI_Py_Pz_Pz_S_C1001001_c;
  abcd[1079] = 4.0E0*I_ERI_Pz_Pz_F3z_S_C1001001_cc-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_c-2.0E0*2*I_ERI_Pz_Pz_Pz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1080] = 4.0E0*I_ERI_Px_S_D2x_Px_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Px_C1000001_d;
  abcd[1081] = 4.0E0*I_ERI_Py_S_D2x_Px_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Px_C1000001_d;
  abcd[1082] = 4.0E0*I_ERI_Pz_S_D2x_Px_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Px_C1000001_d;
  abcd[1092] = 4.0E0*I_ERI_Px_S_Dxy_Px_C1000001_cd;
  abcd[1093] = 4.0E0*I_ERI_Py_S_Dxy_Px_C1000001_cd;
  abcd[1094] = 4.0E0*I_ERI_Pz_S_Dxy_Px_C1000001_cd;
  abcd[1104] = 4.0E0*I_ERI_Px_S_Dxz_Px_C1000001_cd;
  abcd[1105] = 4.0E0*I_ERI_Py_S_Dxz_Px_C1000001_cd;
  abcd[1106] = 4.0E0*I_ERI_Pz_S_Dxz_Px_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1083] = 4.0E0*I_ERI_Px_Px_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Px_C1001001_d;
  abcd[1084] = 4.0E0*I_ERI_Py_Px_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Px_C1001001_d;
  abcd[1085] = 4.0E0*I_ERI_Pz_Px_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Px_C1001001_d;
  abcd[1086] = 4.0E0*I_ERI_Px_Py_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Px_C1001001_d;
  abcd[1087] = 4.0E0*I_ERI_Py_Py_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Px_C1001001_d;
  abcd[1088] = 4.0E0*I_ERI_Pz_Py_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Px_C1001001_d;
  abcd[1089] = 4.0E0*I_ERI_Px_Pz_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Px_C1001001_d;
  abcd[1090] = 4.0E0*I_ERI_Py_Pz_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Px_C1001001_d;
  abcd[1091] = 4.0E0*I_ERI_Pz_Pz_D2x_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Px_C1001001_d;
  abcd[1095] = 4.0E0*I_ERI_Px_Px_Dxy_Px_C1001001_cd;
  abcd[1096] = 4.0E0*I_ERI_Py_Px_Dxy_Px_C1001001_cd;
  abcd[1097] = 4.0E0*I_ERI_Pz_Px_Dxy_Px_C1001001_cd;
  abcd[1098] = 4.0E0*I_ERI_Px_Py_Dxy_Px_C1001001_cd;
  abcd[1099] = 4.0E0*I_ERI_Py_Py_Dxy_Px_C1001001_cd;
  abcd[1100] = 4.0E0*I_ERI_Pz_Py_Dxy_Px_C1001001_cd;
  abcd[1101] = 4.0E0*I_ERI_Px_Pz_Dxy_Px_C1001001_cd;
  abcd[1102] = 4.0E0*I_ERI_Py_Pz_Dxy_Px_C1001001_cd;
  abcd[1103] = 4.0E0*I_ERI_Pz_Pz_Dxy_Px_C1001001_cd;
  abcd[1107] = 4.0E0*I_ERI_Px_Px_Dxz_Px_C1001001_cd;
  abcd[1108] = 4.0E0*I_ERI_Py_Px_Dxz_Px_C1001001_cd;
  abcd[1109] = 4.0E0*I_ERI_Pz_Px_Dxz_Px_C1001001_cd;
  abcd[1110] = 4.0E0*I_ERI_Px_Py_Dxz_Px_C1001001_cd;
  abcd[1111] = 4.0E0*I_ERI_Py_Py_Dxz_Px_C1001001_cd;
  abcd[1112] = 4.0E0*I_ERI_Pz_Py_Dxz_Px_C1001001_cd;
  abcd[1113] = 4.0E0*I_ERI_Px_Pz_Dxz_Px_C1001001_cd;
  abcd[1114] = 4.0E0*I_ERI_Py_Pz_Dxz_Px_C1001001_cd;
  abcd[1115] = 4.0E0*I_ERI_Pz_Pz_Dxz_Px_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1116] = 4.0E0*I_ERI_Px_S_D2x_Py_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Py_C1000001_d;
  abcd[1117] = 4.0E0*I_ERI_Py_S_D2x_Py_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Py_C1000001_d;
  abcd[1118] = 4.0E0*I_ERI_Pz_S_D2x_Py_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Py_C1000001_d;
  abcd[1128] = 4.0E0*I_ERI_Px_S_Dxy_Py_C1000001_cd;
  abcd[1129] = 4.0E0*I_ERI_Py_S_Dxy_Py_C1000001_cd;
  abcd[1130] = 4.0E0*I_ERI_Pz_S_Dxy_Py_C1000001_cd;
  abcd[1140] = 4.0E0*I_ERI_Px_S_Dxz_Py_C1000001_cd;
  abcd[1141] = 4.0E0*I_ERI_Py_S_Dxz_Py_C1000001_cd;
  abcd[1142] = 4.0E0*I_ERI_Pz_S_Dxz_Py_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1119] = 4.0E0*I_ERI_Px_Px_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Py_C1001001_d;
  abcd[1120] = 4.0E0*I_ERI_Py_Px_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Py_C1001001_d;
  abcd[1121] = 4.0E0*I_ERI_Pz_Px_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Py_C1001001_d;
  abcd[1122] = 4.0E0*I_ERI_Px_Py_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Py_C1001001_d;
  abcd[1123] = 4.0E0*I_ERI_Py_Py_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Py_C1001001_d;
  abcd[1124] = 4.0E0*I_ERI_Pz_Py_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Py_C1001001_d;
  abcd[1125] = 4.0E0*I_ERI_Px_Pz_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Py_C1001001_d;
  abcd[1126] = 4.0E0*I_ERI_Py_Pz_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Py_C1001001_d;
  abcd[1127] = 4.0E0*I_ERI_Pz_Pz_D2x_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Py_C1001001_d;
  abcd[1131] = 4.0E0*I_ERI_Px_Px_Dxy_Py_C1001001_cd;
  abcd[1132] = 4.0E0*I_ERI_Py_Px_Dxy_Py_C1001001_cd;
  abcd[1133] = 4.0E0*I_ERI_Pz_Px_Dxy_Py_C1001001_cd;
  abcd[1134] = 4.0E0*I_ERI_Px_Py_Dxy_Py_C1001001_cd;
  abcd[1135] = 4.0E0*I_ERI_Py_Py_Dxy_Py_C1001001_cd;
  abcd[1136] = 4.0E0*I_ERI_Pz_Py_Dxy_Py_C1001001_cd;
  abcd[1137] = 4.0E0*I_ERI_Px_Pz_Dxy_Py_C1001001_cd;
  abcd[1138] = 4.0E0*I_ERI_Py_Pz_Dxy_Py_C1001001_cd;
  abcd[1139] = 4.0E0*I_ERI_Pz_Pz_Dxy_Py_C1001001_cd;
  abcd[1143] = 4.0E0*I_ERI_Px_Px_Dxz_Py_C1001001_cd;
  abcd[1144] = 4.0E0*I_ERI_Py_Px_Dxz_Py_C1001001_cd;
  abcd[1145] = 4.0E0*I_ERI_Pz_Px_Dxz_Py_C1001001_cd;
  abcd[1146] = 4.0E0*I_ERI_Px_Py_Dxz_Py_C1001001_cd;
  abcd[1147] = 4.0E0*I_ERI_Py_Py_Dxz_Py_C1001001_cd;
  abcd[1148] = 4.0E0*I_ERI_Pz_Py_Dxz_Py_C1001001_cd;
  abcd[1149] = 4.0E0*I_ERI_Px_Pz_Dxz_Py_C1001001_cd;
  abcd[1150] = 4.0E0*I_ERI_Py_Pz_Dxz_Py_C1001001_cd;
  abcd[1151] = 4.0E0*I_ERI_Pz_Pz_Dxz_Py_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1152] = 4.0E0*I_ERI_Px_S_D2x_Pz_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Pz_C1000001_d;
  abcd[1153] = 4.0E0*I_ERI_Py_S_D2x_Pz_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Pz_C1000001_d;
  abcd[1154] = 4.0E0*I_ERI_Pz_S_D2x_Pz_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Pz_C1000001_d;
  abcd[1164] = 4.0E0*I_ERI_Px_S_Dxy_Pz_C1000001_cd;
  abcd[1165] = 4.0E0*I_ERI_Py_S_Dxy_Pz_C1000001_cd;
  abcd[1166] = 4.0E0*I_ERI_Pz_S_Dxy_Pz_C1000001_cd;
  abcd[1176] = 4.0E0*I_ERI_Px_S_Dxz_Pz_C1000001_cd;
  abcd[1177] = 4.0E0*I_ERI_Py_S_Dxz_Pz_C1000001_cd;
  abcd[1178] = 4.0E0*I_ERI_Pz_S_Dxz_Pz_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1155] = 4.0E0*I_ERI_Px_Px_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Pz_C1001001_d;
  abcd[1156] = 4.0E0*I_ERI_Py_Px_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Pz_C1001001_d;
  abcd[1157] = 4.0E0*I_ERI_Pz_Px_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Pz_C1001001_d;
  abcd[1158] = 4.0E0*I_ERI_Px_Py_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Pz_C1001001_d;
  abcd[1159] = 4.0E0*I_ERI_Py_Py_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Pz_C1001001_d;
  abcd[1160] = 4.0E0*I_ERI_Pz_Py_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Pz_C1001001_d;
  abcd[1161] = 4.0E0*I_ERI_Px_Pz_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Pz_C1001001_d;
  abcd[1162] = 4.0E0*I_ERI_Py_Pz_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Pz_C1001001_d;
  abcd[1163] = 4.0E0*I_ERI_Pz_Pz_D2x_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Pz_C1001001_d;
  abcd[1167] = 4.0E0*I_ERI_Px_Px_Dxy_Pz_C1001001_cd;
  abcd[1168] = 4.0E0*I_ERI_Py_Px_Dxy_Pz_C1001001_cd;
  abcd[1169] = 4.0E0*I_ERI_Pz_Px_Dxy_Pz_C1001001_cd;
  abcd[1170] = 4.0E0*I_ERI_Px_Py_Dxy_Pz_C1001001_cd;
  abcd[1171] = 4.0E0*I_ERI_Py_Py_Dxy_Pz_C1001001_cd;
  abcd[1172] = 4.0E0*I_ERI_Pz_Py_Dxy_Pz_C1001001_cd;
  abcd[1173] = 4.0E0*I_ERI_Px_Pz_Dxy_Pz_C1001001_cd;
  abcd[1174] = 4.0E0*I_ERI_Py_Pz_Dxy_Pz_C1001001_cd;
  abcd[1175] = 4.0E0*I_ERI_Pz_Pz_Dxy_Pz_C1001001_cd;
  abcd[1179] = 4.0E0*I_ERI_Px_Px_Dxz_Pz_C1001001_cd;
  abcd[1180] = 4.0E0*I_ERI_Py_Px_Dxz_Pz_C1001001_cd;
  abcd[1181] = 4.0E0*I_ERI_Pz_Px_Dxz_Pz_C1001001_cd;
  abcd[1182] = 4.0E0*I_ERI_Px_Py_Dxz_Pz_C1001001_cd;
  abcd[1183] = 4.0E0*I_ERI_Py_Py_Dxz_Pz_C1001001_cd;
  abcd[1184] = 4.0E0*I_ERI_Pz_Py_Dxz_Pz_C1001001_cd;
  abcd[1185] = 4.0E0*I_ERI_Px_Pz_Dxz_Pz_C1001001_cd;
  abcd[1186] = 4.0E0*I_ERI_Py_Pz_Dxz_Pz_C1001001_cd;
  abcd[1187] = 4.0E0*I_ERI_Pz_Pz_Dxz_Pz_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1188] = 4.0E0*I_ERI_Px_S_Dxy_Px_C1000001_cd;
  abcd[1189] = 4.0E0*I_ERI_Py_S_Dxy_Px_C1000001_cd;
  abcd[1190] = 4.0E0*I_ERI_Pz_S_Dxy_Px_C1000001_cd;
  abcd[1200] = 4.0E0*I_ERI_Px_S_D2y_Px_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Px_C1000001_d;
  abcd[1201] = 4.0E0*I_ERI_Py_S_D2y_Px_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Px_C1000001_d;
  abcd[1202] = 4.0E0*I_ERI_Pz_S_D2y_Px_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Px_C1000001_d;
  abcd[1212] = 4.0E0*I_ERI_Px_S_Dyz_Px_C1000001_cd;
  abcd[1213] = 4.0E0*I_ERI_Py_S_Dyz_Px_C1000001_cd;
  abcd[1214] = 4.0E0*I_ERI_Pz_S_Dyz_Px_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1191] = 4.0E0*I_ERI_Px_Px_Dxy_Px_C1001001_cd;
  abcd[1192] = 4.0E0*I_ERI_Py_Px_Dxy_Px_C1001001_cd;
  abcd[1193] = 4.0E0*I_ERI_Pz_Px_Dxy_Px_C1001001_cd;
  abcd[1194] = 4.0E0*I_ERI_Px_Py_Dxy_Px_C1001001_cd;
  abcd[1195] = 4.0E0*I_ERI_Py_Py_Dxy_Px_C1001001_cd;
  abcd[1196] = 4.0E0*I_ERI_Pz_Py_Dxy_Px_C1001001_cd;
  abcd[1197] = 4.0E0*I_ERI_Px_Pz_Dxy_Px_C1001001_cd;
  abcd[1198] = 4.0E0*I_ERI_Py_Pz_Dxy_Px_C1001001_cd;
  abcd[1199] = 4.0E0*I_ERI_Pz_Pz_Dxy_Px_C1001001_cd;
  abcd[1203] = 4.0E0*I_ERI_Px_Px_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Px_C1001001_d;
  abcd[1204] = 4.0E0*I_ERI_Py_Px_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Px_C1001001_d;
  abcd[1205] = 4.0E0*I_ERI_Pz_Px_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Px_C1001001_d;
  abcd[1206] = 4.0E0*I_ERI_Px_Py_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Px_C1001001_d;
  abcd[1207] = 4.0E0*I_ERI_Py_Py_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Px_C1001001_d;
  abcd[1208] = 4.0E0*I_ERI_Pz_Py_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Px_C1001001_d;
  abcd[1209] = 4.0E0*I_ERI_Px_Pz_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Px_C1001001_d;
  abcd[1210] = 4.0E0*I_ERI_Py_Pz_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Px_C1001001_d;
  abcd[1211] = 4.0E0*I_ERI_Pz_Pz_D2y_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Px_C1001001_d;
  abcd[1215] = 4.0E0*I_ERI_Px_Px_Dyz_Px_C1001001_cd;
  abcd[1216] = 4.0E0*I_ERI_Py_Px_Dyz_Px_C1001001_cd;
  abcd[1217] = 4.0E0*I_ERI_Pz_Px_Dyz_Px_C1001001_cd;
  abcd[1218] = 4.0E0*I_ERI_Px_Py_Dyz_Px_C1001001_cd;
  abcd[1219] = 4.0E0*I_ERI_Py_Py_Dyz_Px_C1001001_cd;
  abcd[1220] = 4.0E0*I_ERI_Pz_Py_Dyz_Px_C1001001_cd;
  abcd[1221] = 4.0E0*I_ERI_Px_Pz_Dyz_Px_C1001001_cd;
  abcd[1222] = 4.0E0*I_ERI_Py_Pz_Dyz_Px_C1001001_cd;
  abcd[1223] = 4.0E0*I_ERI_Pz_Pz_Dyz_Px_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1224] = 4.0E0*I_ERI_Px_S_Dxy_Py_C1000001_cd;
  abcd[1225] = 4.0E0*I_ERI_Py_S_Dxy_Py_C1000001_cd;
  abcd[1226] = 4.0E0*I_ERI_Pz_S_Dxy_Py_C1000001_cd;
  abcd[1236] = 4.0E0*I_ERI_Px_S_D2y_Py_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Py_C1000001_d;
  abcd[1237] = 4.0E0*I_ERI_Py_S_D2y_Py_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Py_C1000001_d;
  abcd[1238] = 4.0E0*I_ERI_Pz_S_D2y_Py_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Py_C1000001_d;
  abcd[1248] = 4.0E0*I_ERI_Px_S_Dyz_Py_C1000001_cd;
  abcd[1249] = 4.0E0*I_ERI_Py_S_Dyz_Py_C1000001_cd;
  abcd[1250] = 4.0E0*I_ERI_Pz_S_Dyz_Py_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1227] = 4.0E0*I_ERI_Px_Px_Dxy_Py_C1001001_cd;
  abcd[1228] = 4.0E0*I_ERI_Py_Px_Dxy_Py_C1001001_cd;
  abcd[1229] = 4.0E0*I_ERI_Pz_Px_Dxy_Py_C1001001_cd;
  abcd[1230] = 4.0E0*I_ERI_Px_Py_Dxy_Py_C1001001_cd;
  abcd[1231] = 4.0E0*I_ERI_Py_Py_Dxy_Py_C1001001_cd;
  abcd[1232] = 4.0E0*I_ERI_Pz_Py_Dxy_Py_C1001001_cd;
  abcd[1233] = 4.0E0*I_ERI_Px_Pz_Dxy_Py_C1001001_cd;
  abcd[1234] = 4.0E0*I_ERI_Py_Pz_Dxy_Py_C1001001_cd;
  abcd[1235] = 4.0E0*I_ERI_Pz_Pz_Dxy_Py_C1001001_cd;
  abcd[1239] = 4.0E0*I_ERI_Px_Px_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Py_C1001001_d;
  abcd[1240] = 4.0E0*I_ERI_Py_Px_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Py_C1001001_d;
  abcd[1241] = 4.0E0*I_ERI_Pz_Px_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Py_C1001001_d;
  abcd[1242] = 4.0E0*I_ERI_Px_Py_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Py_C1001001_d;
  abcd[1243] = 4.0E0*I_ERI_Py_Py_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Py_C1001001_d;
  abcd[1244] = 4.0E0*I_ERI_Pz_Py_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Py_C1001001_d;
  abcd[1245] = 4.0E0*I_ERI_Px_Pz_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Py_C1001001_d;
  abcd[1246] = 4.0E0*I_ERI_Py_Pz_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Py_C1001001_d;
  abcd[1247] = 4.0E0*I_ERI_Pz_Pz_D2y_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Py_C1001001_d;
  abcd[1251] = 4.0E0*I_ERI_Px_Px_Dyz_Py_C1001001_cd;
  abcd[1252] = 4.0E0*I_ERI_Py_Px_Dyz_Py_C1001001_cd;
  abcd[1253] = 4.0E0*I_ERI_Pz_Px_Dyz_Py_C1001001_cd;
  abcd[1254] = 4.0E0*I_ERI_Px_Py_Dyz_Py_C1001001_cd;
  abcd[1255] = 4.0E0*I_ERI_Py_Py_Dyz_Py_C1001001_cd;
  abcd[1256] = 4.0E0*I_ERI_Pz_Py_Dyz_Py_C1001001_cd;
  abcd[1257] = 4.0E0*I_ERI_Px_Pz_Dyz_Py_C1001001_cd;
  abcd[1258] = 4.0E0*I_ERI_Py_Pz_Dyz_Py_C1001001_cd;
  abcd[1259] = 4.0E0*I_ERI_Pz_Pz_Dyz_Py_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1260] = 4.0E0*I_ERI_Px_S_Dxy_Pz_C1000001_cd;
  abcd[1261] = 4.0E0*I_ERI_Py_S_Dxy_Pz_C1000001_cd;
  abcd[1262] = 4.0E0*I_ERI_Pz_S_Dxy_Pz_C1000001_cd;
  abcd[1272] = 4.0E0*I_ERI_Px_S_D2y_Pz_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Pz_C1000001_d;
  abcd[1273] = 4.0E0*I_ERI_Py_S_D2y_Pz_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Pz_C1000001_d;
  abcd[1274] = 4.0E0*I_ERI_Pz_S_D2y_Pz_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Pz_C1000001_d;
  abcd[1284] = 4.0E0*I_ERI_Px_S_Dyz_Pz_C1000001_cd;
  abcd[1285] = 4.0E0*I_ERI_Py_S_Dyz_Pz_C1000001_cd;
  abcd[1286] = 4.0E0*I_ERI_Pz_S_Dyz_Pz_C1000001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1263] = 4.0E0*I_ERI_Px_Px_Dxy_Pz_C1001001_cd;
  abcd[1264] = 4.0E0*I_ERI_Py_Px_Dxy_Pz_C1001001_cd;
  abcd[1265] = 4.0E0*I_ERI_Pz_Px_Dxy_Pz_C1001001_cd;
  abcd[1266] = 4.0E0*I_ERI_Px_Py_Dxy_Pz_C1001001_cd;
  abcd[1267] = 4.0E0*I_ERI_Py_Py_Dxy_Pz_C1001001_cd;
  abcd[1268] = 4.0E0*I_ERI_Pz_Py_Dxy_Pz_C1001001_cd;
  abcd[1269] = 4.0E0*I_ERI_Px_Pz_Dxy_Pz_C1001001_cd;
  abcd[1270] = 4.0E0*I_ERI_Py_Pz_Dxy_Pz_C1001001_cd;
  abcd[1271] = 4.0E0*I_ERI_Pz_Pz_Dxy_Pz_C1001001_cd;
  abcd[1275] = 4.0E0*I_ERI_Px_Px_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Pz_C1001001_d;
  abcd[1276] = 4.0E0*I_ERI_Py_Px_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Pz_C1001001_d;
  abcd[1277] = 4.0E0*I_ERI_Pz_Px_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Pz_C1001001_d;
  abcd[1278] = 4.0E0*I_ERI_Px_Py_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Pz_C1001001_d;
  abcd[1279] = 4.0E0*I_ERI_Py_Py_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Pz_C1001001_d;
  abcd[1280] = 4.0E0*I_ERI_Pz_Py_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Pz_C1001001_d;
  abcd[1281] = 4.0E0*I_ERI_Px_Pz_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Pz_C1001001_d;
  abcd[1282] = 4.0E0*I_ERI_Py_Pz_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Pz_C1001001_d;
  abcd[1283] = 4.0E0*I_ERI_Pz_Pz_D2y_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Pz_C1001001_d;
  abcd[1287] = 4.0E0*I_ERI_Px_Px_Dyz_Pz_C1001001_cd;
  abcd[1288] = 4.0E0*I_ERI_Py_Px_Dyz_Pz_C1001001_cd;
  abcd[1289] = 4.0E0*I_ERI_Pz_Px_Dyz_Pz_C1001001_cd;
  abcd[1290] = 4.0E0*I_ERI_Px_Py_Dyz_Pz_C1001001_cd;
  abcd[1291] = 4.0E0*I_ERI_Py_Py_Dyz_Pz_C1001001_cd;
  abcd[1292] = 4.0E0*I_ERI_Pz_Py_Dyz_Pz_C1001001_cd;
  abcd[1293] = 4.0E0*I_ERI_Px_Pz_Dyz_Pz_C1001001_cd;
  abcd[1294] = 4.0E0*I_ERI_Py_Pz_Dyz_Pz_C1001001_cd;
  abcd[1295] = 4.0E0*I_ERI_Pz_Pz_Dyz_Pz_C1001001_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1296] = 4.0E0*I_ERI_Px_S_Dxz_Px_C1000001_cd;
  abcd[1297] = 4.0E0*I_ERI_Py_S_Dxz_Px_C1000001_cd;
  abcd[1298] = 4.0E0*I_ERI_Pz_S_Dxz_Px_C1000001_cd;
  abcd[1308] = 4.0E0*I_ERI_Px_S_Dyz_Px_C1000001_cd;
  abcd[1309] = 4.0E0*I_ERI_Py_S_Dyz_Px_C1000001_cd;
  abcd[1310] = 4.0E0*I_ERI_Pz_S_Dyz_Px_C1000001_cd;
  abcd[1320] = 4.0E0*I_ERI_Px_S_D2z_Px_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Px_C1000001_d;
  abcd[1321] = 4.0E0*I_ERI_Py_S_D2z_Px_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Px_C1000001_d;
  abcd[1322] = 4.0E0*I_ERI_Pz_S_D2z_Px_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Px_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1299] = 4.0E0*I_ERI_Px_Px_Dxz_Px_C1001001_cd;
  abcd[1300] = 4.0E0*I_ERI_Py_Px_Dxz_Px_C1001001_cd;
  abcd[1301] = 4.0E0*I_ERI_Pz_Px_Dxz_Px_C1001001_cd;
  abcd[1302] = 4.0E0*I_ERI_Px_Py_Dxz_Px_C1001001_cd;
  abcd[1303] = 4.0E0*I_ERI_Py_Py_Dxz_Px_C1001001_cd;
  abcd[1304] = 4.0E0*I_ERI_Pz_Py_Dxz_Px_C1001001_cd;
  abcd[1305] = 4.0E0*I_ERI_Px_Pz_Dxz_Px_C1001001_cd;
  abcd[1306] = 4.0E0*I_ERI_Py_Pz_Dxz_Px_C1001001_cd;
  abcd[1307] = 4.0E0*I_ERI_Pz_Pz_Dxz_Px_C1001001_cd;
  abcd[1311] = 4.0E0*I_ERI_Px_Px_Dyz_Px_C1001001_cd;
  abcd[1312] = 4.0E0*I_ERI_Py_Px_Dyz_Px_C1001001_cd;
  abcd[1313] = 4.0E0*I_ERI_Pz_Px_Dyz_Px_C1001001_cd;
  abcd[1314] = 4.0E0*I_ERI_Px_Py_Dyz_Px_C1001001_cd;
  abcd[1315] = 4.0E0*I_ERI_Py_Py_Dyz_Px_C1001001_cd;
  abcd[1316] = 4.0E0*I_ERI_Pz_Py_Dyz_Px_C1001001_cd;
  abcd[1317] = 4.0E0*I_ERI_Px_Pz_Dyz_Px_C1001001_cd;
  abcd[1318] = 4.0E0*I_ERI_Py_Pz_Dyz_Px_C1001001_cd;
  abcd[1319] = 4.0E0*I_ERI_Pz_Pz_Dyz_Px_C1001001_cd;
  abcd[1323] = 4.0E0*I_ERI_Px_Px_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Px_C1001001_d;
  abcd[1324] = 4.0E0*I_ERI_Py_Px_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Px_C1001001_d;
  abcd[1325] = 4.0E0*I_ERI_Pz_Px_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Px_C1001001_d;
  abcd[1326] = 4.0E0*I_ERI_Px_Py_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Px_C1001001_d;
  abcd[1327] = 4.0E0*I_ERI_Py_Py_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Px_C1001001_d;
  abcd[1328] = 4.0E0*I_ERI_Pz_Py_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Px_C1001001_d;
  abcd[1329] = 4.0E0*I_ERI_Px_Pz_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Px_C1001001_d;
  abcd[1330] = 4.0E0*I_ERI_Py_Pz_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Px_C1001001_d;
  abcd[1331] = 4.0E0*I_ERI_Pz_Pz_D2z_Px_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Px_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1332] = 4.0E0*I_ERI_Px_S_Dxz_Py_C1000001_cd;
  abcd[1333] = 4.0E0*I_ERI_Py_S_Dxz_Py_C1000001_cd;
  abcd[1334] = 4.0E0*I_ERI_Pz_S_Dxz_Py_C1000001_cd;
  abcd[1344] = 4.0E0*I_ERI_Px_S_Dyz_Py_C1000001_cd;
  abcd[1345] = 4.0E0*I_ERI_Py_S_Dyz_Py_C1000001_cd;
  abcd[1346] = 4.0E0*I_ERI_Pz_S_Dyz_Py_C1000001_cd;
  abcd[1356] = 4.0E0*I_ERI_Px_S_D2z_Py_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Py_C1000001_d;
  abcd[1357] = 4.0E0*I_ERI_Py_S_D2z_Py_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Py_C1000001_d;
  abcd[1358] = 4.0E0*I_ERI_Pz_S_D2z_Py_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Py_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1335] = 4.0E0*I_ERI_Px_Px_Dxz_Py_C1001001_cd;
  abcd[1336] = 4.0E0*I_ERI_Py_Px_Dxz_Py_C1001001_cd;
  abcd[1337] = 4.0E0*I_ERI_Pz_Px_Dxz_Py_C1001001_cd;
  abcd[1338] = 4.0E0*I_ERI_Px_Py_Dxz_Py_C1001001_cd;
  abcd[1339] = 4.0E0*I_ERI_Py_Py_Dxz_Py_C1001001_cd;
  abcd[1340] = 4.0E0*I_ERI_Pz_Py_Dxz_Py_C1001001_cd;
  abcd[1341] = 4.0E0*I_ERI_Px_Pz_Dxz_Py_C1001001_cd;
  abcd[1342] = 4.0E0*I_ERI_Py_Pz_Dxz_Py_C1001001_cd;
  abcd[1343] = 4.0E0*I_ERI_Pz_Pz_Dxz_Py_C1001001_cd;
  abcd[1347] = 4.0E0*I_ERI_Px_Px_Dyz_Py_C1001001_cd;
  abcd[1348] = 4.0E0*I_ERI_Py_Px_Dyz_Py_C1001001_cd;
  abcd[1349] = 4.0E0*I_ERI_Pz_Px_Dyz_Py_C1001001_cd;
  abcd[1350] = 4.0E0*I_ERI_Px_Py_Dyz_Py_C1001001_cd;
  abcd[1351] = 4.0E0*I_ERI_Py_Py_Dyz_Py_C1001001_cd;
  abcd[1352] = 4.0E0*I_ERI_Pz_Py_Dyz_Py_C1001001_cd;
  abcd[1353] = 4.0E0*I_ERI_Px_Pz_Dyz_Py_C1001001_cd;
  abcd[1354] = 4.0E0*I_ERI_Py_Pz_Dyz_Py_C1001001_cd;
  abcd[1355] = 4.0E0*I_ERI_Pz_Pz_Dyz_Py_C1001001_cd;
  abcd[1359] = 4.0E0*I_ERI_Px_Px_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Py_C1001001_d;
  abcd[1360] = 4.0E0*I_ERI_Py_Px_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Py_C1001001_d;
  abcd[1361] = 4.0E0*I_ERI_Pz_Px_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Py_C1001001_d;
  abcd[1362] = 4.0E0*I_ERI_Px_Py_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Py_C1001001_d;
  abcd[1363] = 4.0E0*I_ERI_Py_Py_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Py_C1001001_d;
  abcd[1364] = 4.0E0*I_ERI_Pz_Py_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Py_C1001001_d;
  abcd[1365] = 4.0E0*I_ERI_Px_Pz_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Py_C1001001_d;
  abcd[1366] = 4.0E0*I_ERI_Py_Pz_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Py_C1001001_d;
  abcd[1367] = 4.0E0*I_ERI_Pz_Pz_D2z_Py_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Py_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1000001_cd
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001_d
   ************************************************************/
  abcd[1368] = 4.0E0*I_ERI_Px_S_Dxz_Pz_C1000001_cd;
  abcd[1369] = 4.0E0*I_ERI_Py_S_Dxz_Pz_C1000001_cd;
  abcd[1370] = 4.0E0*I_ERI_Pz_S_Dxz_Pz_C1000001_cd;
  abcd[1380] = 4.0E0*I_ERI_Px_S_Dyz_Pz_C1000001_cd;
  abcd[1381] = 4.0E0*I_ERI_Py_S_Dyz_Pz_C1000001_cd;
  abcd[1382] = 4.0E0*I_ERI_Pz_S_Dyz_Pz_C1000001_cd;
  abcd[1392] = 4.0E0*I_ERI_Px_S_D2z_Pz_C1000001_cd-2.0E0*1*I_ERI_Px_S_S_Pz_C1000001_d;
  abcd[1393] = 4.0E0*I_ERI_Py_S_D2z_Pz_C1000001_cd-2.0E0*1*I_ERI_Py_S_S_Pz_C1000001_d;
  abcd[1394] = 4.0E0*I_ERI_Pz_S_D2z_Pz_C1000001_cd-2.0E0*1*I_ERI_Pz_S_S_Pz_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001_cd
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001_d
   ************************************************************/
  abcd[1371] = 4.0E0*I_ERI_Px_Px_Dxz_Pz_C1001001_cd;
  abcd[1372] = 4.0E0*I_ERI_Py_Px_Dxz_Pz_C1001001_cd;
  abcd[1373] = 4.0E0*I_ERI_Pz_Px_Dxz_Pz_C1001001_cd;
  abcd[1374] = 4.0E0*I_ERI_Px_Py_Dxz_Pz_C1001001_cd;
  abcd[1375] = 4.0E0*I_ERI_Py_Py_Dxz_Pz_C1001001_cd;
  abcd[1376] = 4.0E0*I_ERI_Pz_Py_Dxz_Pz_C1001001_cd;
  abcd[1377] = 4.0E0*I_ERI_Px_Pz_Dxz_Pz_C1001001_cd;
  abcd[1378] = 4.0E0*I_ERI_Py_Pz_Dxz_Pz_C1001001_cd;
  abcd[1379] = 4.0E0*I_ERI_Pz_Pz_Dxz_Pz_C1001001_cd;
  abcd[1383] = 4.0E0*I_ERI_Px_Px_Dyz_Pz_C1001001_cd;
  abcd[1384] = 4.0E0*I_ERI_Py_Px_Dyz_Pz_C1001001_cd;
  abcd[1385] = 4.0E0*I_ERI_Pz_Px_Dyz_Pz_C1001001_cd;
  abcd[1386] = 4.0E0*I_ERI_Px_Py_Dyz_Pz_C1001001_cd;
  abcd[1387] = 4.0E0*I_ERI_Py_Py_Dyz_Pz_C1001001_cd;
  abcd[1388] = 4.0E0*I_ERI_Pz_Py_Dyz_Pz_C1001001_cd;
  abcd[1389] = 4.0E0*I_ERI_Px_Pz_Dyz_Pz_C1001001_cd;
  abcd[1390] = 4.0E0*I_ERI_Py_Pz_Dyz_Pz_C1001001_cd;
  abcd[1391] = 4.0E0*I_ERI_Pz_Pz_Dyz_Pz_C1001001_cd;
  abcd[1395] = 4.0E0*I_ERI_Px_Px_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Px_S_Pz_C1001001_d;
  abcd[1396] = 4.0E0*I_ERI_Py_Px_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Px_S_Pz_C1001001_d;
  abcd[1397] = 4.0E0*I_ERI_Pz_Px_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Px_S_Pz_C1001001_d;
  abcd[1398] = 4.0E0*I_ERI_Px_Py_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Py_S_Pz_C1001001_d;
  abcd[1399] = 4.0E0*I_ERI_Py_Py_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Py_S_Pz_C1001001_d;
  abcd[1400] = 4.0E0*I_ERI_Pz_Py_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Py_S_Pz_C1001001_d;
  abcd[1401] = 4.0E0*I_ERI_Px_Pz_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Px_Pz_S_Pz_C1001001_d;
  abcd[1402] = 4.0E0*I_ERI_Py_Pz_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Py_Pz_S_Pz_C1001001_d;
  abcd[1403] = 4.0E0*I_ERI_Pz_Pz_D2z_Pz_C1001001_cd-2.0E0*1*I_ERI_Pz_Pz_S_Pz_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1404] = 4.0E0*I_ERI_Px_S_Px_D2x_C1000001_dd-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_d;
  abcd[1405] = 4.0E0*I_ERI_Py_S_Px_D2x_C1000001_dd-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_d;
  abcd[1406] = 4.0E0*I_ERI_Pz_S_Px_D2x_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_d;
  abcd[1416] = 4.0E0*I_ERI_Px_S_Py_D2x_C1000001_dd-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_d;
  abcd[1417] = 4.0E0*I_ERI_Py_S_Py_D2x_C1000001_dd-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_d;
  abcd[1418] = 4.0E0*I_ERI_Pz_S_Py_D2x_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_d;
  abcd[1428] = 4.0E0*I_ERI_Px_S_Pz_D2x_C1000001_dd-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_d;
  abcd[1429] = 4.0E0*I_ERI_Py_S_Pz_D2x_C1000001_dd-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_d;
  abcd[1430] = 4.0E0*I_ERI_Pz_S_Pz_D2x_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1407] = 4.0E0*I_ERI_Px_Px_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_d;
  abcd[1408] = 4.0E0*I_ERI_Py_Px_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_d;
  abcd[1409] = 4.0E0*I_ERI_Pz_Px_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_d;
  abcd[1410] = 4.0E0*I_ERI_Px_Py_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_d;
  abcd[1411] = 4.0E0*I_ERI_Py_Py_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_d;
  abcd[1412] = 4.0E0*I_ERI_Pz_Py_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_d;
  abcd[1413] = 4.0E0*I_ERI_Px_Pz_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_d;
  abcd[1414] = 4.0E0*I_ERI_Py_Pz_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_d;
  abcd[1415] = 4.0E0*I_ERI_Pz_Pz_Px_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_d;
  abcd[1419] = 4.0E0*I_ERI_Px_Px_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_d;
  abcd[1420] = 4.0E0*I_ERI_Py_Px_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_d;
  abcd[1421] = 4.0E0*I_ERI_Pz_Px_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_d;
  abcd[1422] = 4.0E0*I_ERI_Px_Py_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_d;
  abcd[1423] = 4.0E0*I_ERI_Py_Py_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_d;
  abcd[1424] = 4.0E0*I_ERI_Pz_Py_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_d;
  abcd[1425] = 4.0E0*I_ERI_Px_Pz_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_d;
  abcd[1426] = 4.0E0*I_ERI_Py_Pz_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_d;
  abcd[1427] = 4.0E0*I_ERI_Pz_Pz_Py_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_d;
  abcd[1431] = 4.0E0*I_ERI_Px_Px_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_d;
  abcd[1432] = 4.0E0*I_ERI_Py_Px_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_d;
  abcd[1433] = 4.0E0*I_ERI_Pz_Px_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_d;
  abcd[1434] = 4.0E0*I_ERI_Px_Py_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_d;
  abcd[1435] = 4.0E0*I_ERI_Py_Py_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_d;
  abcd[1436] = 4.0E0*I_ERI_Pz_Py_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_d;
  abcd[1437] = 4.0E0*I_ERI_Px_Pz_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_d;
  abcd[1438] = 4.0E0*I_ERI_Py_Pz_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_d;
  abcd[1439] = 4.0E0*I_ERI_Pz_Pz_Pz_D2x_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1440] = 4.0E0*I_ERI_Px_S_Px_Dxy_C1000001_dd;
  abcd[1441] = 4.0E0*I_ERI_Py_S_Px_Dxy_C1000001_dd;
  abcd[1442] = 4.0E0*I_ERI_Pz_S_Px_Dxy_C1000001_dd;
  abcd[1452] = 4.0E0*I_ERI_Px_S_Py_Dxy_C1000001_dd;
  abcd[1453] = 4.0E0*I_ERI_Py_S_Py_Dxy_C1000001_dd;
  abcd[1454] = 4.0E0*I_ERI_Pz_S_Py_Dxy_C1000001_dd;
  abcd[1464] = 4.0E0*I_ERI_Px_S_Pz_Dxy_C1000001_dd;
  abcd[1465] = 4.0E0*I_ERI_Py_S_Pz_Dxy_C1000001_dd;
  abcd[1466] = 4.0E0*I_ERI_Pz_S_Pz_Dxy_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1443] = 4.0E0*I_ERI_Px_Px_Px_Dxy_C1001001_dd;
  abcd[1444] = 4.0E0*I_ERI_Py_Px_Px_Dxy_C1001001_dd;
  abcd[1445] = 4.0E0*I_ERI_Pz_Px_Px_Dxy_C1001001_dd;
  abcd[1446] = 4.0E0*I_ERI_Px_Py_Px_Dxy_C1001001_dd;
  abcd[1447] = 4.0E0*I_ERI_Py_Py_Px_Dxy_C1001001_dd;
  abcd[1448] = 4.0E0*I_ERI_Pz_Py_Px_Dxy_C1001001_dd;
  abcd[1449] = 4.0E0*I_ERI_Px_Pz_Px_Dxy_C1001001_dd;
  abcd[1450] = 4.0E0*I_ERI_Py_Pz_Px_Dxy_C1001001_dd;
  abcd[1451] = 4.0E0*I_ERI_Pz_Pz_Px_Dxy_C1001001_dd;
  abcd[1455] = 4.0E0*I_ERI_Px_Px_Py_Dxy_C1001001_dd;
  abcd[1456] = 4.0E0*I_ERI_Py_Px_Py_Dxy_C1001001_dd;
  abcd[1457] = 4.0E0*I_ERI_Pz_Px_Py_Dxy_C1001001_dd;
  abcd[1458] = 4.0E0*I_ERI_Px_Py_Py_Dxy_C1001001_dd;
  abcd[1459] = 4.0E0*I_ERI_Py_Py_Py_Dxy_C1001001_dd;
  abcd[1460] = 4.0E0*I_ERI_Pz_Py_Py_Dxy_C1001001_dd;
  abcd[1461] = 4.0E0*I_ERI_Px_Pz_Py_Dxy_C1001001_dd;
  abcd[1462] = 4.0E0*I_ERI_Py_Pz_Py_Dxy_C1001001_dd;
  abcd[1463] = 4.0E0*I_ERI_Pz_Pz_Py_Dxy_C1001001_dd;
  abcd[1467] = 4.0E0*I_ERI_Px_Px_Pz_Dxy_C1001001_dd;
  abcd[1468] = 4.0E0*I_ERI_Py_Px_Pz_Dxy_C1001001_dd;
  abcd[1469] = 4.0E0*I_ERI_Pz_Px_Pz_Dxy_C1001001_dd;
  abcd[1470] = 4.0E0*I_ERI_Px_Py_Pz_Dxy_C1001001_dd;
  abcd[1471] = 4.0E0*I_ERI_Py_Py_Pz_Dxy_C1001001_dd;
  abcd[1472] = 4.0E0*I_ERI_Pz_Py_Pz_Dxy_C1001001_dd;
  abcd[1473] = 4.0E0*I_ERI_Px_Pz_Pz_Dxy_C1001001_dd;
  abcd[1474] = 4.0E0*I_ERI_Py_Pz_Pz_Dxy_C1001001_dd;
  abcd[1475] = 4.0E0*I_ERI_Pz_Pz_Pz_Dxy_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1476] = 4.0E0*I_ERI_Px_S_Px_Dxz_C1000001_dd;
  abcd[1477] = 4.0E0*I_ERI_Py_S_Px_Dxz_C1000001_dd;
  abcd[1478] = 4.0E0*I_ERI_Pz_S_Px_Dxz_C1000001_dd;
  abcd[1488] = 4.0E0*I_ERI_Px_S_Py_Dxz_C1000001_dd;
  abcd[1489] = 4.0E0*I_ERI_Py_S_Py_Dxz_C1000001_dd;
  abcd[1490] = 4.0E0*I_ERI_Pz_S_Py_Dxz_C1000001_dd;
  abcd[1500] = 4.0E0*I_ERI_Px_S_Pz_Dxz_C1000001_dd;
  abcd[1501] = 4.0E0*I_ERI_Py_S_Pz_Dxz_C1000001_dd;
  abcd[1502] = 4.0E0*I_ERI_Pz_S_Pz_Dxz_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1479] = 4.0E0*I_ERI_Px_Px_Px_Dxz_C1001001_dd;
  abcd[1480] = 4.0E0*I_ERI_Py_Px_Px_Dxz_C1001001_dd;
  abcd[1481] = 4.0E0*I_ERI_Pz_Px_Px_Dxz_C1001001_dd;
  abcd[1482] = 4.0E0*I_ERI_Px_Py_Px_Dxz_C1001001_dd;
  abcd[1483] = 4.0E0*I_ERI_Py_Py_Px_Dxz_C1001001_dd;
  abcd[1484] = 4.0E0*I_ERI_Pz_Py_Px_Dxz_C1001001_dd;
  abcd[1485] = 4.0E0*I_ERI_Px_Pz_Px_Dxz_C1001001_dd;
  abcd[1486] = 4.0E0*I_ERI_Py_Pz_Px_Dxz_C1001001_dd;
  abcd[1487] = 4.0E0*I_ERI_Pz_Pz_Px_Dxz_C1001001_dd;
  abcd[1491] = 4.0E0*I_ERI_Px_Px_Py_Dxz_C1001001_dd;
  abcd[1492] = 4.0E0*I_ERI_Py_Px_Py_Dxz_C1001001_dd;
  abcd[1493] = 4.0E0*I_ERI_Pz_Px_Py_Dxz_C1001001_dd;
  abcd[1494] = 4.0E0*I_ERI_Px_Py_Py_Dxz_C1001001_dd;
  abcd[1495] = 4.0E0*I_ERI_Py_Py_Py_Dxz_C1001001_dd;
  abcd[1496] = 4.0E0*I_ERI_Pz_Py_Py_Dxz_C1001001_dd;
  abcd[1497] = 4.0E0*I_ERI_Px_Pz_Py_Dxz_C1001001_dd;
  abcd[1498] = 4.0E0*I_ERI_Py_Pz_Py_Dxz_C1001001_dd;
  abcd[1499] = 4.0E0*I_ERI_Pz_Pz_Py_Dxz_C1001001_dd;
  abcd[1503] = 4.0E0*I_ERI_Px_Px_Pz_Dxz_C1001001_dd;
  abcd[1504] = 4.0E0*I_ERI_Py_Px_Pz_Dxz_C1001001_dd;
  abcd[1505] = 4.0E0*I_ERI_Pz_Px_Pz_Dxz_C1001001_dd;
  abcd[1506] = 4.0E0*I_ERI_Px_Py_Pz_Dxz_C1001001_dd;
  abcd[1507] = 4.0E0*I_ERI_Py_Py_Pz_Dxz_C1001001_dd;
  abcd[1508] = 4.0E0*I_ERI_Pz_Py_Pz_Dxz_C1001001_dd;
  abcd[1509] = 4.0E0*I_ERI_Px_Pz_Pz_Dxz_C1001001_dd;
  abcd[1510] = 4.0E0*I_ERI_Py_Pz_Pz_Dxz_C1001001_dd;
  abcd[1511] = 4.0E0*I_ERI_Pz_Pz_Pz_Dxz_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1512] = 4.0E0*I_ERI_Px_S_Px_D2y_C1000001_dd-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_d;
  abcd[1513] = 4.0E0*I_ERI_Py_S_Px_D2y_C1000001_dd-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_d;
  abcd[1514] = 4.0E0*I_ERI_Pz_S_Px_D2y_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_d;
  abcd[1524] = 4.0E0*I_ERI_Px_S_Py_D2y_C1000001_dd-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_d;
  abcd[1525] = 4.0E0*I_ERI_Py_S_Py_D2y_C1000001_dd-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_d;
  abcd[1526] = 4.0E0*I_ERI_Pz_S_Py_D2y_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_d;
  abcd[1536] = 4.0E0*I_ERI_Px_S_Pz_D2y_C1000001_dd-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_d;
  abcd[1537] = 4.0E0*I_ERI_Py_S_Pz_D2y_C1000001_dd-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_d;
  abcd[1538] = 4.0E0*I_ERI_Pz_S_Pz_D2y_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1515] = 4.0E0*I_ERI_Px_Px_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_d;
  abcd[1516] = 4.0E0*I_ERI_Py_Px_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_d;
  abcd[1517] = 4.0E0*I_ERI_Pz_Px_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_d;
  abcd[1518] = 4.0E0*I_ERI_Px_Py_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_d;
  abcd[1519] = 4.0E0*I_ERI_Py_Py_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_d;
  abcd[1520] = 4.0E0*I_ERI_Pz_Py_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_d;
  abcd[1521] = 4.0E0*I_ERI_Px_Pz_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_d;
  abcd[1522] = 4.0E0*I_ERI_Py_Pz_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_d;
  abcd[1523] = 4.0E0*I_ERI_Pz_Pz_Px_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_d;
  abcd[1527] = 4.0E0*I_ERI_Px_Px_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_d;
  abcd[1528] = 4.0E0*I_ERI_Py_Px_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_d;
  abcd[1529] = 4.0E0*I_ERI_Pz_Px_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_d;
  abcd[1530] = 4.0E0*I_ERI_Px_Py_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_d;
  abcd[1531] = 4.0E0*I_ERI_Py_Py_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_d;
  abcd[1532] = 4.0E0*I_ERI_Pz_Py_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_d;
  abcd[1533] = 4.0E0*I_ERI_Px_Pz_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_d;
  abcd[1534] = 4.0E0*I_ERI_Py_Pz_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_d;
  abcd[1535] = 4.0E0*I_ERI_Pz_Pz_Py_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_d;
  abcd[1539] = 4.0E0*I_ERI_Px_Px_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_d;
  abcd[1540] = 4.0E0*I_ERI_Py_Px_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_d;
  abcd[1541] = 4.0E0*I_ERI_Pz_Px_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_d;
  abcd[1542] = 4.0E0*I_ERI_Px_Py_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_d;
  abcd[1543] = 4.0E0*I_ERI_Py_Py_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_d;
  abcd[1544] = 4.0E0*I_ERI_Pz_Py_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_d;
  abcd[1545] = 4.0E0*I_ERI_Px_Pz_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_d;
  abcd[1546] = 4.0E0*I_ERI_Py_Pz_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_d;
  abcd[1547] = 4.0E0*I_ERI_Pz_Pz_Pz_D2y_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1548] = 4.0E0*I_ERI_Px_S_Px_Dyz_C1000001_dd;
  abcd[1549] = 4.0E0*I_ERI_Py_S_Px_Dyz_C1000001_dd;
  abcd[1550] = 4.0E0*I_ERI_Pz_S_Px_Dyz_C1000001_dd;
  abcd[1560] = 4.0E0*I_ERI_Px_S_Py_Dyz_C1000001_dd;
  abcd[1561] = 4.0E0*I_ERI_Py_S_Py_Dyz_C1000001_dd;
  abcd[1562] = 4.0E0*I_ERI_Pz_S_Py_Dyz_C1000001_dd;
  abcd[1572] = 4.0E0*I_ERI_Px_S_Pz_Dyz_C1000001_dd;
  abcd[1573] = 4.0E0*I_ERI_Py_S_Pz_Dyz_C1000001_dd;
  abcd[1574] = 4.0E0*I_ERI_Pz_S_Pz_Dyz_C1000001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1551] = 4.0E0*I_ERI_Px_Px_Px_Dyz_C1001001_dd;
  abcd[1552] = 4.0E0*I_ERI_Py_Px_Px_Dyz_C1001001_dd;
  abcd[1553] = 4.0E0*I_ERI_Pz_Px_Px_Dyz_C1001001_dd;
  abcd[1554] = 4.0E0*I_ERI_Px_Py_Px_Dyz_C1001001_dd;
  abcd[1555] = 4.0E0*I_ERI_Py_Py_Px_Dyz_C1001001_dd;
  abcd[1556] = 4.0E0*I_ERI_Pz_Py_Px_Dyz_C1001001_dd;
  abcd[1557] = 4.0E0*I_ERI_Px_Pz_Px_Dyz_C1001001_dd;
  abcd[1558] = 4.0E0*I_ERI_Py_Pz_Px_Dyz_C1001001_dd;
  abcd[1559] = 4.0E0*I_ERI_Pz_Pz_Px_Dyz_C1001001_dd;
  abcd[1563] = 4.0E0*I_ERI_Px_Px_Py_Dyz_C1001001_dd;
  abcd[1564] = 4.0E0*I_ERI_Py_Px_Py_Dyz_C1001001_dd;
  abcd[1565] = 4.0E0*I_ERI_Pz_Px_Py_Dyz_C1001001_dd;
  abcd[1566] = 4.0E0*I_ERI_Px_Py_Py_Dyz_C1001001_dd;
  abcd[1567] = 4.0E0*I_ERI_Py_Py_Py_Dyz_C1001001_dd;
  abcd[1568] = 4.0E0*I_ERI_Pz_Py_Py_Dyz_C1001001_dd;
  abcd[1569] = 4.0E0*I_ERI_Px_Pz_Py_Dyz_C1001001_dd;
  abcd[1570] = 4.0E0*I_ERI_Py_Pz_Py_Dyz_C1001001_dd;
  abcd[1571] = 4.0E0*I_ERI_Pz_Pz_Py_Dyz_C1001001_dd;
  abcd[1575] = 4.0E0*I_ERI_Px_Px_Pz_Dyz_C1001001_dd;
  abcd[1576] = 4.0E0*I_ERI_Py_Px_Pz_Dyz_C1001001_dd;
  abcd[1577] = 4.0E0*I_ERI_Pz_Px_Pz_Dyz_C1001001_dd;
  abcd[1578] = 4.0E0*I_ERI_Px_Py_Pz_Dyz_C1001001_dd;
  abcd[1579] = 4.0E0*I_ERI_Py_Py_Pz_Dyz_C1001001_dd;
  abcd[1580] = 4.0E0*I_ERI_Pz_Py_Pz_Dyz_C1001001_dd;
  abcd[1581] = 4.0E0*I_ERI_Px_Pz_Pz_Dyz_C1001001_dd;
  abcd[1582] = 4.0E0*I_ERI_Py_Pz_Pz_Dyz_C1001001_dd;
  abcd[1583] = 4.0E0*I_ERI_Pz_Pz_Pz_Dyz_C1001001_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_D_C1000001_dd
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  abcd[1584] = 4.0E0*I_ERI_Px_S_Px_D2z_C1000001_dd-2.0E0*1*I_ERI_Px_S_Px_S_C1000001_d;
  abcd[1585] = 4.0E0*I_ERI_Py_S_Px_D2z_C1000001_dd-2.0E0*1*I_ERI_Py_S_Px_S_C1000001_d;
  abcd[1586] = 4.0E0*I_ERI_Pz_S_Px_D2z_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Px_S_C1000001_d;
  abcd[1596] = 4.0E0*I_ERI_Px_S_Py_D2z_C1000001_dd-2.0E0*1*I_ERI_Px_S_Py_S_C1000001_d;
  abcd[1597] = 4.0E0*I_ERI_Py_S_Py_D2z_C1000001_dd-2.0E0*1*I_ERI_Py_S_Py_S_C1000001_d;
  abcd[1598] = 4.0E0*I_ERI_Pz_S_Py_D2z_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Py_S_C1000001_d;
  abcd[1608] = 4.0E0*I_ERI_Px_S_Pz_D2z_C1000001_dd-2.0E0*1*I_ERI_Px_S_Pz_S_C1000001_d;
  abcd[1609] = 4.0E0*I_ERI_Py_S_Pz_D2z_C1000001_dd-2.0E0*1*I_ERI_Py_S_Pz_S_C1000001_d;
  abcd[1610] = 4.0E0*I_ERI_Pz_S_Pz_D2z_C1000001_dd-2.0E0*1*I_ERI_Pz_S_Pz_S_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_D_C1001001_dd
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_d
   ************************************************************/
  abcd[1587] = 4.0E0*I_ERI_Px_Px_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Px_S_C1001001_d;
  abcd[1588] = 4.0E0*I_ERI_Py_Px_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Px_S_C1001001_d;
  abcd[1589] = 4.0E0*I_ERI_Pz_Px_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Px_S_C1001001_d;
  abcd[1590] = 4.0E0*I_ERI_Px_Py_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Px_S_C1001001_d;
  abcd[1591] = 4.0E0*I_ERI_Py_Py_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Px_S_C1001001_d;
  abcd[1592] = 4.0E0*I_ERI_Pz_Py_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Px_S_C1001001_d;
  abcd[1593] = 4.0E0*I_ERI_Px_Pz_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Px_S_C1001001_d;
  abcd[1594] = 4.0E0*I_ERI_Py_Pz_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Px_S_C1001001_d;
  abcd[1595] = 4.0E0*I_ERI_Pz_Pz_Px_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Px_S_C1001001_d;
  abcd[1599] = 4.0E0*I_ERI_Px_Px_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Py_S_C1001001_d;
  abcd[1600] = 4.0E0*I_ERI_Py_Px_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Py_S_C1001001_d;
  abcd[1601] = 4.0E0*I_ERI_Pz_Px_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Py_S_C1001001_d;
  abcd[1602] = 4.0E0*I_ERI_Px_Py_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Py_S_C1001001_d;
  abcd[1603] = 4.0E0*I_ERI_Py_Py_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Py_S_C1001001_d;
  abcd[1604] = 4.0E0*I_ERI_Pz_Py_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Py_S_C1001001_d;
  abcd[1605] = 4.0E0*I_ERI_Px_Pz_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Py_S_C1001001_d;
  abcd[1606] = 4.0E0*I_ERI_Py_Pz_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Py_S_C1001001_d;
  abcd[1607] = 4.0E0*I_ERI_Pz_Pz_Py_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Py_S_C1001001_d;
  abcd[1611] = 4.0E0*I_ERI_Px_Px_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Px_Pz_S_C1001001_d;
  abcd[1612] = 4.0E0*I_ERI_Py_Px_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Px_Pz_S_C1001001_d;
  abcd[1613] = 4.0E0*I_ERI_Pz_Px_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Px_Pz_S_C1001001_d;
  abcd[1614] = 4.0E0*I_ERI_Px_Py_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Py_Pz_S_C1001001_d;
  abcd[1615] = 4.0E0*I_ERI_Py_Py_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Py_Pz_S_C1001001_d;
  abcd[1616] = 4.0E0*I_ERI_Pz_Py_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Py_Pz_S_C1001001_d;
  abcd[1617] = 4.0E0*I_ERI_Px_Pz_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Px_Pz_Pz_S_C1001001_d;
  abcd[1618] = 4.0E0*I_ERI_Py_Pz_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Py_Pz_Pz_S_C1001001_d;
  abcd[1619] = 4.0E0*I_ERI_Pz_Pz_Pz_D2z_C1001001_dd-2.0E0*1*I_ERI_Pz_Pz_Pz_S_C1001001_d;
}
