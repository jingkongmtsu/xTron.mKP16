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
// BRA1 as redundant position, total RHS integrals evaluated as: 68137
// BRA2 as redundant position, total RHS integrals evaluated as: 74377
// KET1 as redundant position, total RHS integrals evaluated as: 79722
// KET2 as redundant position, total RHS integrals evaluated as: 72285
// the redundant position is: BRA1
//

//
// @@@@ derivative position-direction information
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
// BRA2  KET2
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

void hgp_os_eri_f_s_p_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_Px_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_b = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Px_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3y_S_Px_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3x_S_Py_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3y_S_Py_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_F3x_S_S_Px_d = 0.0E0;
  Double I_ERI_F2xy_S_S_Px_d = 0.0E0;
  Double I_ERI_F2xz_S_S_Px_d = 0.0E0;
  Double I_ERI_Fx2y_S_S_Px_d = 0.0E0;
  Double I_ERI_Fxyz_S_S_Px_d = 0.0E0;
  Double I_ERI_Fx2z_S_S_Px_d = 0.0E0;
  Double I_ERI_F3y_S_S_Px_d = 0.0E0;
  Double I_ERI_F2yz_S_S_Px_d = 0.0E0;
  Double I_ERI_Fy2z_S_S_Px_d = 0.0E0;
  Double I_ERI_F3z_S_S_Px_d = 0.0E0;
  Double I_ERI_F3x_S_S_Py_d = 0.0E0;
  Double I_ERI_F2xy_S_S_Py_d = 0.0E0;
  Double I_ERI_F2xz_S_S_Py_d = 0.0E0;
  Double I_ERI_Fx2y_S_S_Py_d = 0.0E0;
  Double I_ERI_Fxyz_S_S_Py_d = 0.0E0;
  Double I_ERI_Fx2z_S_S_Py_d = 0.0E0;
  Double I_ERI_F3y_S_S_Py_d = 0.0E0;
  Double I_ERI_F2yz_S_S_Py_d = 0.0E0;
  Double I_ERI_Fy2z_S_S_Py_d = 0.0E0;
  Double I_ERI_F3z_S_S_Py_d = 0.0E0;
  Double I_ERI_F3x_S_S_Pz_d = 0.0E0;
  Double I_ERI_F2xy_S_S_Pz_d = 0.0E0;
  Double I_ERI_F2xz_S_S_Pz_d = 0.0E0;
  Double I_ERI_Fx2y_S_S_Pz_d = 0.0E0;
  Double I_ERI_Fxyz_S_S_Pz_d = 0.0E0;
  Double I_ERI_Fx2z_S_S_Pz_d = 0.0E0;
  Double I_ERI_F3y_S_S_Pz_d = 0.0E0;
  Double I_ERI_F2yz_S_S_Pz_d = 0.0E0;
  Double I_ERI_Fy2z_S_S_Pz_d = 0.0E0;
  Double I_ERI_F3z_S_S_Pz_d = 0.0E0;
  Double I_ERI_F3x_S_Px_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_d = 0.0E0;
  Double I_ERI_F3y_S_Px_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_d = 0.0E0;
  Double I_ERI_F3z_S_Px_S_d = 0.0E0;
  Double I_ERI_F3x_S_Py_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_d = 0.0E0;
  Double I_ERI_F3y_S_Py_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_d = 0.0E0;
  Double I_ERI_F3z_S_Py_S_d = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_d = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_d = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_d = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4x_S_S_S_b = 0.0E0;
  Double I_ERI_G3xy_S_S_S_b = 0.0E0;
  Double I_ERI_G3xz_S_S_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_b = 0.0E0;
  Double I_ERI_G4y_S_S_S_b = 0.0E0;
  Double I_ERI_G3yz_S_S_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_b = 0.0E0;
  Double I_ERI_G4z_S_S_S_b = 0.0E0;
  Double I_ERI_F3x_S_S_S_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_b = 0.0E0;
  Double I_ERI_H5x_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_cd = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_cd = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_cd = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Px_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Px_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Px_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Py_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Py_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Py_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Px_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Px_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Px_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Py_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Py_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Py_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_dd = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_dd = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Px_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Px_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Px_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Py_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Py_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Py_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_dd = 0.0E0;

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
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt jp2=0; jp2<jnp2; jp2++) {
      Double onede = jexp[jp2];
      Double eta   = 1.0E0/onede;
      Double ediff = jexpdiff[jp2];
      Double gamma = 0.5E0*(eta+ediff);
      Double delta = 0.5E0*(eta-ediff);
      Double jc2   = jcoe[jp2];
      Double fket  = jfac[jp2];
      Double pref      = fbra*fket;
      Double prefactor = ic2*jc2*pref;

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
      Double I_ERI_S_S_S_S_M6_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M5_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M6_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER47;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER17*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = 1.0E0+u2*ONEOVER15*I_ERI_S_S_S_S_M6_vrr;
        I_ERI_S_S_S_S_M6_vrr = ONEOVER13*I_ERI_S_S_S_S_M6_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M6_vrr  = f*I_ERI_S_S_S_S_M6_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M5_vrr  = ONEOVER11*(u2*I_ERI_S_S_S_S_M6_vrr+f);
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
        I_ERI_S_S_S_S_M6_vrr_d = oneO2u*(11.0E0*I_ERI_S_S_S_S_M5_vrr_d-f);

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);
        I_ERI_S_S_S_S_M5_vrr = static_cast<Double>(I_ERI_S_S_S_S_M5_vrr_d);
        I_ERI_S_S_S_S_M6_vrr = static_cast<Double>(I_ERI_S_S_S_S_M6_vrr_d);

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
        I_ERI_S_S_S_S_M6_vrr = oneO2u*(11.0E0*I_ERI_S_S_S_S_M5_vrr-f);

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
        Double erfPref_13 = erfPref_11*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_Px_S_S_S_M5_vrr = PAX*I_ERI_S_S_S_S_M5_vrr+WPX*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Py_S_S_S_M5_vrr = PAY*I_ERI_S_S_S_S_M5_vrr+WPY*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Pz_S_S_S_M5_vrr = PAZ*I_ERI_S_S_S_S_M5_vrr+WPZ*I_ERI_S_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M4_vrr = PAX*I_ERI_Px_S_S_S_M4_vrr+WPX*I_ERI_Px_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Dxy_S_S_S_M4_vrr = PAY*I_ERI_Px_S_S_S_M4_vrr+WPY*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_D2y_S_S_S_M4_vrr = PAY*I_ERI_Py_S_S_S_M4_vrr+WPY*I_ERI_Py_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_D2z_S_S_S_M4_vrr = PAZ*I_ERI_Pz_S_S_S_M4_vrr+WPZ*I_ERI_Pz_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M3_vrr = PAX*I_ERI_D2x_S_S_S_M3_vrr+WPX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2z*I_ERI_Px_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_S_S_M3_vrr = PAY*I_ERI_D2x_S_S_S_M3_vrr+WPY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_S_S_M3_vrr = PAZ*I_ERI_D2x_S_S_S_M3_vrr+WPZ*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_S_S_M3_vrr = PAX*I_ERI_D2y_S_S_S_M3_vrr+WPX*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fxyz_S_S_S_M3_vrr = PAZ*I_ERI_Dxy_S_S_S_M3_vrr+WPZ*I_ERI_Dxy_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_S_S_M3_vrr = PAX*I_ERI_D2z_S_S_S_M3_vrr+WPX*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_S_S_M3_vrr = PAY*I_ERI_D2y_S_S_S_M3_vrr+WPY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2z*I_ERI_Py_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_S_S_M3_vrr = PAZ*I_ERI_D2y_S_S_S_M3_vrr+WPZ*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fy2z_S_S_S_M3_vrr = PAY*I_ERI_D2z_S_S_S_M3_vrr+WPY*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_S_S_M3_vrr = PAZ*I_ERI_D2z_S_S_S_M3_vrr+WPZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2z*I_ERI_Pz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_F_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M2_vrr = QCX*I_ERI_F3x_S_S_S_M2_vrr+WQX*I_ERI_F3x_S_S_S_M3_vrr+3*oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_Px_S_M2_vrr = QCX*I_ERI_F2xy_S_S_S_M2_vrr+WQX*I_ERI_F2xy_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_Px_S_M2_vrr = QCX*I_ERI_F2xz_S_S_S_M2_vrr+WQX*I_ERI_F2xz_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_Px_S_M2_vrr = QCX*I_ERI_Fx2y_S_S_S_M2_vrr+WQX*I_ERI_Fx2y_S_S_S_M3_vrr+oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_Px_S_M2_vrr = QCX*I_ERI_Fxyz_S_S_S_M2_vrr+WQX*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dyz_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_Px_S_M2_vrr = QCX*I_ERI_Fx2z_S_S_S_M2_vrr+WQX*I_ERI_Fx2z_S_S_S_M3_vrr+oned2k*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_Px_S_M2_vrr = QCX*I_ERI_F3y_S_S_S_M2_vrr+WQX*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_Px_S_M2_vrr = QCX*I_ERI_F2yz_S_S_S_M2_vrr+WQX*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_Px_S_M2_vrr = QCX*I_ERI_Fy2z_S_S_S_M2_vrr+WQX*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_Px_S_M2_vrr = QCX*I_ERI_F3z_S_S_S_M2_vrr+WQX*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_F3x_S_Py_S_M2_vrr = QCY*I_ERI_F3x_S_S_S_M2_vrr+WQY*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_Py_S_M2_vrr = QCY*I_ERI_F2xy_S_S_S_M2_vrr+WQY*I_ERI_F2xy_S_S_S_M3_vrr+oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_Py_S_M2_vrr = QCY*I_ERI_F2xz_S_S_S_M2_vrr+WQY*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_Py_S_M2_vrr = QCY*I_ERI_Fx2y_S_S_S_M2_vrr+WQY*I_ERI_Fx2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_Py_S_M2_vrr = QCY*I_ERI_Fxyz_S_S_S_M2_vrr+WQY*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dxz_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_Py_S_M2_vrr = QCY*I_ERI_Fx2z_S_S_S_M2_vrr+WQY*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_Py_S_M2_vrr = QCY*I_ERI_F3y_S_S_S_M2_vrr+WQY*I_ERI_F3y_S_S_S_M3_vrr+3*oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_Py_S_M2_vrr = QCY*I_ERI_F2yz_S_S_S_M2_vrr+WQY*I_ERI_F2yz_S_S_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_Py_S_M2_vrr = QCY*I_ERI_Fy2z_S_S_S_M2_vrr+WQY*I_ERI_Fy2z_S_S_S_M3_vrr+oned2k*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_Py_S_M2_vrr = QCY*I_ERI_F3z_S_S_S_M2_vrr+WQY*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_F3x_S_Pz_S_M2_vrr = QCZ*I_ERI_F3x_S_S_S_M2_vrr+WQZ*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_Pz_S_M2_vrr = QCZ*I_ERI_F2xy_S_S_S_M2_vrr+WQZ*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_Pz_S_M2_vrr = QCZ*I_ERI_F2xz_S_S_S_M2_vrr+WQZ*I_ERI_F2xz_S_S_S_M3_vrr+oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M2_vrr = QCZ*I_ERI_Fx2y_S_S_S_M2_vrr+WQZ*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_Pz_S_M2_vrr = QCZ*I_ERI_Fxyz_S_S_S_M2_vrr+WQZ*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M2_vrr = QCZ*I_ERI_Fx2z_S_S_S_M2_vrr+WQZ*I_ERI_Fx2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_Pz_S_M2_vrr = QCZ*I_ERI_F3y_S_S_S_M2_vrr+WQZ*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_Pz_S_M2_vrr = QCZ*I_ERI_F2yz_S_S_S_M2_vrr+WQZ*I_ERI_F2yz_S_S_S_M3_vrr+oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_Pz_S_M2_vrr = QCZ*I_ERI_Fy2z_S_S_S_M2_vrr+WQZ*I_ERI_Fy2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_Pz_S_M2_vrr = QCZ*I_ERI_F3z_S_S_S_M2_vrr+WQZ*I_ERI_F3z_S_S_S_M3_vrr+3*oned2k*I_ERI_D2z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 12 integrals are omitted 
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
      Double I_ERI_D2x_S_D2y_S_M1_vrr = QCY*I_ERI_D2x_S_Py_S_M1_vrr+WQY*I_ERI_D2x_S_Py_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = QCY*I_ERI_Dxy_S_Py_S_M1_vrr+WQY*I_ERI_Dxy_S_Py_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = QCY*I_ERI_Dxz_S_Py_S_M1_vrr+WQY*I_ERI_Dxz_S_Py_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = QCY*I_ERI_D2y_S_Py_S_M1_vrr+WQY*I_ERI_D2y_S_Py_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = QCY*I_ERI_Dyz_S_Py_S_M1_vrr+WQY*I_ERI_Dyz_S_Py_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = QCY*I_ERI_D2z_S_Py_S_M1_vrr+WQY*I_ERI_D2z_S_Py_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_D2z_S_M1_vrr = QCZ*I_ERI_D2x_S_Pz_S_M1_vrr+WQZ*I_ERI_D2x_S_Pz_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxy_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_D2z_S_M1_vrr = QCZ*I_ERI_D2y_S_Pz_S_M1_vrr+WQZ*I_ERI_D2y_S_Pz_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dyz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_D2z_S_M1_vrr = QCZ*I_ERI_D2z_S_Pz_S_M1_vrr+WQZ*I_ERI_D2z_S_Pz_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M2_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M1_vrr = PAX*I_ERI_F3x_S_S_S_M1_vrr+WPX*I_ERI_F3x_S_S_S_M2_vrr+3*oned2z*I_ERI_D2x_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_S_S_M1_vrr = PAY*I_ERI_F3x_S_S_S_M1_vrr+WPY*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_S_S_M1_vrr = PAZ*I_ERI_F3x_S_S_S_M1_vrr+WPZ*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_S_S_M1_vrr = PAY*I_ERI_F2xy_S_S_S_M1_vrr+WPY*I_ERI_F2xy_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_S_S_M1_vrr = PAZ*I_ERI_F2xy_S_S_S_M1_vrr+WPZ*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_S_S_M1_vrr = PAZ*I_ERI_F2xz_S_S_S_M1_vrr+WPZ*I_ERI_F2xz_S_S_S_M2_vrr+oned2z*I_ERI_D2x_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_S_S_M1_vrr = PAX*I_ERI_F3y_S_S_S_M1_vrr+WPX*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_S_S_M1_vrr = PAZ*I_ERI_Fx2y_S_S_S_M1_vrr+WPZ*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_S_S_M1_vrr = PAY*I_ERI_Fx2z_S_S_S_M1_vrr+WPY*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_S_S_M1_vrr = PAX*I_ERI_F3z_S_S_S_M1_vrr+WPX*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_S_S_M1_vrr = PAY*I_ERI_F3y_S_S_S_M1_vrr+WPY*I_ERI_F3y_S_S_S_M2_vrr+3*oned2z*I_ERI_D2y_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_S_S_M1_vrr = PAZ*I_ERI_F3y_S_S_S_M1_vrr+WPZ*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_S_S_M1_vrr = PAZ*I_ERI_F2yz_S_S_S_M1_vrr+WPZ*I_ERI_F2yz_S_S_S_M2_vrr+oned2z*I_ERI_D2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_S_S_M1_vrr = PAY*I_ERI_F3z_S_S_S_M1_vrr+WPY*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_S_S_M1_vrr = PAZ*I_ERI_F3z_S_S_S_M1_vrr+WPZ*I_ERI_F3z_S_S_S_M2_vrr+3*oned2z*I_ERI_D2z_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 20 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M1_vrr = QCX*I_ERI_F3x_S_Px_S_M1_vrr+WQX*I_ERI_F3x_S_Px_S_M2_vrr+oned2e*I_ERI_F3x_S_S_S_M1_vrr-rhod2esq*I_ERI_F3x_S_S_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xy_S_D2x_S_M1_vrr = QCX*I_ERI_F2xy_S_Px_S_M1_vrr+WQX*I_ERI_F2xy_S_Px_S_M2_vrr+oned2e*I_ERI_F2xy_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_D2x_S_M1_vrr = QCX*I_ERI_F2xz_S_Px_S_M1_vrr+WQX*I_ERI_F2xz_S_Px_S_M2_vrr+oned2e*I_ERI_F2xz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2x_S_M1_vrr = QCX*I_ERI_Fx2y_S_Px_S_M1_vrr+WQX*I_ERI_Fx2y_S_Px_S_M2_vrr+oned2e*I_ERI_Fx2y_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2x_S_M1_vrr = QCX*I_ERI_Fxyz_S_Px_S_M1_vrr+WQX*I_ERI_Fxyz_S_Px_S_M2_vrr+oned2e*I_ERI_Fxyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dyz_S_Px_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2x_S_M1_vrr = QCX*I_ERI_Fx2z_S_Px_S_M1_vrr+WQX*I_ERI_Fx2z_S_Px_S_M2_vrr+oned2e*I_ERI_Fx2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M2_vrr+oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3y_S_D2x_S_M1_vrr = QCX*I_ERI_F3y_S_Px_S_M1_vrr+WQX*I_ERI_F3y_S_Px_S_M2_vrr+oned2e*I_ERI_F3y_S_S_S_M1_vrr-rhod2esq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_D2x_S_M1_vrr = QCX*I_ERI_F2yz_S_Px_S_M1_vrr+WQX*I_ERI_F2yz_S_Px_S_M2_vrr+oned2e*I_ERI_F2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2x_S_M1_vrr = QCX*I_ERI_Fy2z_S_Px_S_M1_vrr+WQX*I_ERI_Fy2z_S_Px_S_M2_vrr+oned2e*I_ERI_Fy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_D2x_S_M1_vrr = QCX*I_ERI_F3z_S_Px_S_M1_vrr+WQX*I_ERI_F3z_S_Px_S_M2_vrr+oned2e*I_ERI_F3z_S_S_S_M1_vrr-rhod2esq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_F3x_S_Dxy_S_M1_vrr = QCY*I_ERI_F3x_S_Px_S_M1_vrr+WQY*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M1_vrr = QCY*I_ERI_F2xy_S_Px_S_M1_vrr+WQY*I_ERI_F2xy_S_Px_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxy_S_M1_vrr = QCY*I_ERI_F2xz_S_Px_S_M1_vrr+WQY*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_M1_vrr = QCY*I_ERI_Fx2y_S_Px_S_M1_vrr+WQY*I_ERI_Fx2y_S_Px_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M2_vrr;
      Double I_ERI_Fxyz_S_Dxy_S_M1_vrr = QCY*I_ERI_Fxyz_S_Px_S_M1_vrr+WQY*I_ERI_Fxyz_S_Px_S_M2_vrr+oned2k*I_ERI_Dxz_S_Px_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_M1_vrr = QCY*I_ERI_Fx2z_S_Px_S_M1_vrr+WQY*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_F3y_S_Dxy_S_M1_vrr = QCY*I_ERI_F3y_S_Px_S_M1_vrr+WQY*I_ERI_F3y_S_Px_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxy_S_M1_vrr = QCY*I_ERI_F2yz_S_Px_S_M1_vrr+WQY*I_ERI_F2yz_S_Px_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_Px_S_M2_vrr;
      Double I_ERI_Fy2z_S_Dxy_S_M1_vrr = QCY*I_ERI_Fy2z_S_Px_S_M1_vrr+WQY*I_ERI_Fy2z_S_Px_S_M2_vrr+oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Dxy_S_M1_vrr = QCY*I_ERI_F3z_S_Px_S_M1_vrr+WQY*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_D2y_S_M1_vrr = QCY*I_ERI_F3x_S_Py_S_M1_vrr+WQY*I_ERI_F3x_S_Py_S_M2_vrr+oned2e*I_ERI_F3x_S_S_S_M1_vrr-rhod2esq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_D2y_S_M1_vrr = QCY*I_ERI_F2xy_S_Py_S_M1_vrr+WQY*I_ERI_F2xy_S_Py_S_M2_vrr+oned2e*I_ERI_F2xy_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_F2xz_S_D2y_S_M1_vrr = QCY*I_ERI_F2xz_S_Py_S_M1_vrr+WQY*I_ERI_F2xz_S_Py_S_M2_vrr+oned2e*I_ERI_F2xz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2y_S_M1_vrr = QCY*I_ERI_Fx2y_S_Py_S_M1_vrr+WQY*I_ERI_Fx2y_S_Py_S_M2_vrr+oned2e*I_ERI_Fx2y_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_Py_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2y_S_M1_vrr = QCY*I_ERI_Fxyz_S_Py_S_M1_vrr+WQY*I_ERI_Fxyz_S_Py_S_M2_vrr+oned2e*I_ERI_Fxyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dxz_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2y_S_M1_vrr = QCY*I_ERI_Fx2z_S_Py_S_M1_vrr+WQY*I_ERI_Fx2z_S_Py_S_M2_vrr+oned2e*I_ERI_Fx2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_D2y_S_M1_vrr = QCY*I_ERI_F3y_S_Py_S_M1_vrr+WQY*I_ERI_F3y_S_Py_S_M2_vrr+oned2e*I_ERI_F3y_S_S_S_M1_vrr-rhod2esq*I_ERI_F3y_S_S_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F2yz_S_D2y_S_M1_vrr = QCY*I_ERI_F2yz_S_Py_S_M1_vrr+WQY*I_ERI_F2yz_S_Py_S_M2_vrr+oned2e*I_ERI_F2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_Py_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2y_S_M1_vrr = QCY*I_ERI_Fy2z_S_Py_S_M1_vrr+WQY*I_ERI_Fy2z_S_Py_S_M2_vrr+oned2e*I_ERI_Fy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M2_vrr+oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_D2y_S_M1_vrr = QCY*I_ERI_F3z_S_Py_S_M1_vrr+WQY*I_ERI_F3z_S_Py_S_M2_vrr+oned2e*I_ERI_F3z_S_S_S_M1_vrr-rhod2esq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_F3x_S_D2z_S_M1_vrr = QCZ*I_ERI_F3x_S_Pz_S_M1_vrr+WQZ*I_ERI_F3x_S_Pz_S_M2_vrr+oned2e*I_ERI_F3x_S_S_S_M1_vrr-rhod2esq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_D2z_S_M1_vrr = QCZ*I_ERI_F2xy_S_Pz_S_M1_vrr+WQZ*I_ERI_F2xy_S_Pz_S_M2_vrr+oned2e*I_ERI_F2xy_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_D2z_S_M1_vrr = QCZ*I_ERI_F2xz_S_Pz_S_M1_vrr+WQZ*I_ERI_F2xz_S_Pz_S_M2_vrr+oned2e*I_ERI_F2xz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M2_vrr+oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2z_S_M1_vrr = QCZ*I_ERI_Fx2y_S_Pz_S_M1_vrr+WQZ*I_ERI_Fx2y_S_Pz_S_M2_vrr+oned2e*I_ERI_Fx2y_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2z_S_M1_vrr = QCZ*I_ERI_Fxyz_S_Pz_S_M1_vrr+WQZ*I_ERI_Fxyz_S_Pz_S_M2_vrr+oned2e*I_ERI_Fxyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M2_vrr+oned2k*I_ERI_Dxy_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2z_S_M1_vrr = QCZ*I_ERI_Fx2z_S_Pz_S_M1_vrr+WQZ*I_ERI_Fx2z_S_Pz_S_M2_vrr+oned2e*I_ERI_Fx2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_Pz_S_M2_vrr;
      Double I_ERI_F3y_S_D2z_S_M1_vrr = QCZ*I_ERI_F3y_S_Pz_S_M1_vrr+WQZ*I_ERI_F3y_S_Pz_S_M2_vrr+oned2e*I_ERI_F3y_S_S_S_M1_vrr-rhod2esq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_D2z_S_M1_vrr = QCZ*I_ERI_F2yz_S_Pz_S_M1_vrr+WQZ*I_ERI_F2yz_S_Pz_S_M2_vrr+oned2e*I_ERI_F2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M2_vrr+oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2z_S_M1_vrr = QCZ*I_ERI_Fy2z_S_Pz_S_M1_vrr+WQZ*I_ERI_Fy2z_S_Pz_S_M2_vrr+oned2e*I_ERI_Fy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_Pz_S_M2_vrr;
      Double I_ERI_F3z_S_D2z_S_M1_vrr = QCZ*I_ERI_F3z_S_Pz_S_M1_vrr+WQZ*I_ERI_F3z_S_Pz_S_M2_vrr+oned2e*I_ERI_F3z_S_S_S_M1_vrr-rhod2esq*I_ERI_F3z_S_S_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M1_vrr = PAX*I_ERI_F3x_S_Px_S_M1_vrr+WPX*I_ERI_F3x_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Px_S_M1_vrr = PAY*I_ERI_F3x_S_Px_S_M1_vrr+WPY*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Px_S_M1_vrr = PAZ*I_ERI_F3x_S_Px_S_M1_vrr+WPZ*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Px_S_M1_vrr = PAY*I_ERI_F2xy_S_Px_S_M1_vrr+WPY*I_ERI_F2xy_S_Px_S_M2_vrr+oned2z*I_ERI_D2x_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_Px_S_M1_vrr = PAZ*I_ERI_F2xy_S_Px_S_M1_vrr+WPZ*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_Px_S_M1_vrr = PAZ*I_ERI_F2xz_S_Px_S_M1_vrr+WPZ*I_ERI_F2xz_S_Px_S_M2_vrr+oned2z*I_ERI_D2x_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Px_S_M1_vrr = PAX*I_ERI_F3y_S_Px_S_M1_vrr+WPX*I_ERI_F3y_S_Px_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Px_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M1_vrr = PAY*I_ERI_Fx2z_S_Px_S_M1_vrr+WPY*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Px_S_M1_vrr = PAX*I_ERI_F3z_S_Px_S_M1_vrr+WPX*I_ERI_F3z_S_Px_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Px_S_M1_vrr = PAY*I_ERI_F3y_S_Px_S_M1_vrr+WPY*I_ERI_F3y_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Px_S_M1_vrr = PAZ*I_ERI_F3y_S_Px_S_M1_vrr+WPZ*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Px_S_M1_vrr = PAZ*I_ERI_F2yz_S_Px_S_M1_vrr+WPZ*I_ERI_F2yz_S_Px_S_M2_vrr+oned2z*I_ERI_D2y_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Px_S_M1_vrr = PAY*I_ERI_F3z_S_Px_S_M1_vrr+WPY*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Px_S_M1_vrr = PAZ*I_ERI_F3z_S_Px_S_M1_vrr+WPZ*I_ERI_F3z_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_Py_S_M1_vrr = PAX*I_ERI_F3x_S_Py_S_M1_vrr+WPX*I_ERI_F3x_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_G3xy_S_Py_S_M1_vrr = PAY*I_ERI_F3x_S_Py_S_M1_vrr+WPY*I_ERI_F3x_S_Py_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Py_S_M1_vrr = PAZ*I_ERI_F3x_S_Py_S_M1_vrr+WPZ*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G2x2y_S_Py_S_M1_vrr = PAY*I_ERI_F2xy_S_Py_S_M1_vrr+WPY*I_ERI_F2xy_S_Py_S_M2_vrr+oned2z*I_ERI_D2x_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Py_S_M1_vrr = PAZ*I_ERI_F2xy_S_Py_S_M1_vrr+WPZ*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_G2x2z_S_Py_S_M1_vrr = PAZ*I_ERI_F2xz_S_Py_S_M1_vrr+WPZ*I_ERI_F2xz_S_Py_S_M2_vrr+oned2z*I_ERI_D2x_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Gx3y_S_Py_S_M1_vrr = PAX*I_ERI_F3y_S_Py_S_M1_vrr+WPX*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Py_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Py_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M1_vrr = PAY*I_ERI_Fx2z_S_Py_S_M1_vrr+WPY*I_ERI_Fx2z_S_Py_S_M2_vrr+oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Py_S_M1_vrr = PAX*I_ERI_F3z_S_Py_S_M1_vrr+WPX*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4y_S_Py_S_M1_vrr = PAY*I_ERI_F3y_S_Py_S_M1_vrr+WPY*I_ERI_F3y_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Py_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Py_S_M1_vrr = PAZ*I_ERI_F3y_S_Py_S_M1_vrr+WPZ*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_Py_S_M1_vrr = PAZ*I_ERI_F2yz_S_Py_S_M1_vrr+WPZ*I_ERI_F2yz_S_Py_S_M2_vrr+oned2z*I_ERI_D2y_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_Py_S_M1_vrr = PAY*I_ERI_F3z_S_Py_S_M1_vrr+WPY*I_ERI_F3z_S_Py_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Py_S_M1_vrr = PAZ*I_ERI_F3z_S_Py_S_M1_vrr+WPZ*I_ERI_F3z_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_G4x_S_Pz_S_M1_vrr = PAX*I_ERI_F3x_S_Pz_S_M1_vrr+WPX*I_ERI_F3x_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xy_S_Pz_S_M1_vrr = PAY*I_ERI_F3x_S_Pz_S_M1_vrr+WPY*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xz_S_Pz_S_M1_vrr = PAZ*I_ERI_F3x_S_Pz_S_M1_vrr+WPZ*I_ERI_F3x_S_Pz_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M1_vrr = PAY*I_ERI_F2xy_S_Pz_S_M1_vrr+WPY*I_ERI_F2xy_S_Pz_S_M2_vrr+oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M1_vrr = PAZ*I_ERI_F2xy_S_Pz_S_M1_vrr+WPZ*I_ERI_F2xy_S_Pz_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Pz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Pz_S_M2_vrr+oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr+oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M1_vrr = PAX*I_ERI_F3y_S_Pz_S_M1_vrr+WPX*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Pz_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Pz_S_M2_vrr+oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M1_vrr = PAY*I_ERI_Fx2z_S_Pz_S_M1_vrr+WPY*I_ERI_Fx2z_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M1_vrr = PAX*I_ERI_F3z_S_Pz_S_M1_vrr+WPX*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4y_S_Pz_S_M1_vrr = PAY*I_ERI_F3y_S_Pz_S_M1_vrr+WPY*I_ERI_F3y_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_G3yz_S_Pz_S_M1_vrr = PAZ*I_ERI_F3y_S_Pz_S_M1_vrr+WPZ*I_ERI_F3y_S_Pz_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Pz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Pz_S_M2_vrr+oned2z*I_ERI_D2y_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Pz_S_M2_vrr+oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M1_vrr = PAY*I_ERI_F3z_S_Pz_S_M1_vrr+WPY*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4z_S_Pz_S_M1_vrr = PAZ*I_ERI_F3z_S_Pz_S_M1_vrr+WPZ*I_ERI_F3z_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Pz_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_S_vrr = PAY*I_ERI_Px_S_S_S_vrr+WPY*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_S_vrr = PAZ*I_ERI_Pz_S_S_S_vrr+WPZ*I_ERI_Pz_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_P
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_S_Px_vrr = QDX*I_ERI_F3x_S_S_S_vrr+WQX*I_ERI_F3x_S_S_S_M1_vrr+3*oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_Px_vrr = QDX*I_ERI_F2xy_S_S_S_vrr+WQX*I_ERI_F2xy_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_Px_vrr = QDX*I_ERI_F2xz_S_S_S_vrr+WQX*I_ERI_F2xz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_Px_vrr = QDX*I_ERI_Fx2y_S_S_S_vrr+WQX*I_ERI_Fx2y_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_S_Px_vrr = QDX*I_ERI_Fxyz_S_S_S_vrr+WQX*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_Px_vrr = QDX*I_ERI_Fx2z_S_S_S_vrr+WQX*I_ERI_Fx2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_Px_vrr = QDX*I_ERI_F3y_S_S_S_vrr+WQX*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_Px_vrr = QDX*I_ERI_F2yz_S_S_S_vrr+WQX*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_S_Px_vrr = QDX*I_ERI_Fy2z_S_S_S_vrr+WQX*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_Px_vrr = QDX*I_ERI_F3z_S_S_S_vrr+WQX*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_S_Py_vrr = QDY*I_ERI_F3x_S_S_S_vrr+WQY*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_Py_vrr = QDY*I_ERI_F2xy_S_S_S_vrr+WQY*I_ERI_F2xy_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_Py_vrr = QDY*I_ERI_F2xz_S_S_S_vrr+WQY*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_Py_vrr = QDY*I_ERI_Fx2y_S_S_S_vrr+WQY*I_ERI_Fx2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_S_Py_vrr = QDY*I_ERI_Fxyz_S_S_S_vrr+WQY*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_Py_vrr = QDY*I_ERI_Fx2z_S_S_S_vrr+WQY*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_Py_vrr = QDY*I_ERI_F3y_S_S_S_vrr+WQY*I_ERI_F3y_S_S_S_M1_vrr+3*oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_Py_vrr = QDY*I_ERI_F2yz_S_S_S_vrr+WQY*I_ERI_F2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_S_Py_vrr = QDY*I_ERI_Fy2z_S_S_S_vrr+WQY*I_ERI_Fy2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_Py_vrr = QDY*I_ERI_F3z_S_S_S_vrr+WQY*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_S_Pz_vrr = QDZ*I_ERI_F3x_S_S_S_vrr+WQZ*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_Pz_vrr = QDZ*I_ERI_F2xy_S_S_S_vrr+WQZ*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_Pz_vrr = QDZ*I_ERI_F2xz_S_S_S_vrr+WQZ*I_ERI_F2xz_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_Pz_vrr = QDZ*I_ERI_Fx2y_S_S_S_vrr+WQZ*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_S_Pz_vrr = QDZ*I_ERI_Fxyz_S_S_S_vrr+WQZ*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_Pz_vrr = QDZ*I_ERI_Fx2z_S_S_S_vrr+WQZ*I_ERI_Fx2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_Pz_vrr = QDZ*I_ERI_F3y_S_S_S_vrr+WQZ*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_Pz_vrr = QDZ*I_ERI_F2yz_S_S_S_vrr+WQZ*I_ERI_F2yz_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_S_Pz_vrr = QDZ*I_ERI_Fy2z_S_S_S_vrr+WQZ*I_ERI_Fy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_Pz_vrr = QDZ*I_ERI_F3z_S_S_S_vrr+WQZ*I_ERI_F3z_S_S_S_M1_vrr+3*oned2k*I_ERI_D2z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_vrr = QCX*I_ERI_F3x_S_S_S_vrr+WQX*I_ERI_F3x_S_S_S_M1_vrr+3*oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_Px_S_vrr = QCX*I_ERI_F2xy_S_S_S_vrr+WQX*I_ERI_F2xy_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_Px_S_vrr = QCX*I_ERI_F2xz_S_S_S_vrr+WQX*I_ERI_F2xz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_Px_S_vrr = QCX*I_ERI_Fx2y_S_S_S_vrr+WQX*I_ERI_Fx2y_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_Px_S_vrr = QCX*I_ERI_Fxyz_S_S_S_vrr+WQX*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_Px_S_vrr = QCX*I_ERI_Fx2z_S_S_S_vrr+WQX*I_ERI_Fx2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_Px_S_vrr = QCX*I_ERI_F3y_S_S_S_vrr+WQX*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_Px_S_vrr = QCX*I_ERI_F2yz_S_S_S_vrr+WQX*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_Px_S_vrr = QCX*I_ERI_Fy2z_S_S_S_vrr+WQX*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_Px_S_vrr = QCX*I_ERI_F3z_S_S_S_vrr+WQX*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_Py_S_vrr = QCY*I_ERI_F3x_S_S_S_vrr+WQY*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_Py_S_vrr = QCY*I_ERI_F2xy_S_S_S_vrr+WQY*I_ERI_F2xy_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_Py_S_vrr = QCY*I_ERI_F2xz_S_S_S_vrr+WQY*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_Py_S_vrr = QCY*I_ERI_Fx2y_S_S_S_vrr+WQY*I_ERI_Fx2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_Py_S_vrr = QCY*I_ERI_Fxyz_S_S_S_vrr+WQY*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_Py_S_vrr = QCY*I_ERI_Fx2z_S_S_S_vrr+WQY*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_Py_S_vrr = QCY*I_ERI_F3y_S_S_S_vrr+WQY*I_ERI_F3y_S_S_S_M1_vrr+3*oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_Py_S_vrr = QCY*I_ERI_F2yz_S_S_S_vrr+WQY*I_ERI_F2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_Py_S_vrr = QCY*I_ERI_Fy2z_S_S_S_vrr+WQY*I_ERI_Fy2z_S_S_S_M1_vrr+oned2k*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_Py_S_vrr = QCY*I_ERI_F3z_S_S_S_vrr+WQY*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_F3x_S_Pz_S_vrr = QCZ*I_ERI_F3x_S_S_S_vrr+WQZ*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_Pz_S_vrr = QCZ*I_ERI_F2xy_S_S_S_vrr+WQZ*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_Pz_S_vrr = QCZ*I_ERI_F2xz_S_S_S_vrr+WQZ*I_ERI_F2xz_S_S_S_M1_vrr+oned2k*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_Pz_S_vrr = QCZ*I_ERI_Fx2y_S_S_S_vrr+WQZ*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Fxyz_S_Pz_S_vrr = QCZ*I_ERI_Fxyz_S_S_S_vrr+WQZ*I_ERI_Fxyz_S_S_S_M1_vrr+oned2k*I_ERI_Dxy_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_Pz_S_vrr = QCZ*I_ERI_Fx2z_S_S_S_vrr+WQZ*I_ERI_Fx2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_Pz_S_vrr = QCZ*I_ERI_F3y_S_S_S_vrr+WQZ*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_Pz_S_vrr = QCZ*I_ERI_F2yz_S_S_S_vrr+WQZ*I_ERI_F2yz_S_S_S_M1_vrr+oned2k*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fy2z_S_Pz_S_vrr = QCZ*I_ERI_Fy2z_S_S_S_vrr+WQZ*I_ERI_Fy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_Pz_S_vrr = QCZ*I_ERI_F3z_S_S_S_vrr+WQZ*I_ERI_F3z_S_S_S_M1_vrr+3*oned2k*I_ERI_D2z_S_S_S_M1_vrr;

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
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_vrr = QCX*I_ERI_G4x_S_S_S_vrr+WQX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Px_S_vrr = QCX*I_ERI_G3xy_S_S_S_vrr+WQX*I_ERI_G3xy_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Px_S_vrr = QCX*I_ERI_G3xz_S_S_S_vrr+WQX*I_ERI_G3xz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Px_S_vrr = QCX*I_ERI_G2x2y_S_S_S_vrr+WQX*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Px_S_vrr = QCX*I_ERI_G2xyz_S_S_S_vrr+WQX*I_ERI_G2xyz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Px_S_vrr = QCX*I_ERI_G2x2z_S_S_S_vrr+WQX*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Px_S_vrr = QCX*I_ERI_Gx3y_S_S_S_vrr+WQX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Px_S_vrr = QCX*I_ERI_Gx2yz_S_S_S_vrr+WQX*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Px_S_vrr = QCX*I_ERI_Gxy2z_S_S_S_vrr+WQX*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Px_S_vrr = QCX*I_ERI_Gx3z_S_S_S_vrr+WQX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Px_S_vrr = QCX*I_ERI_G4y_S_S_S_vrr+WQX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Px_S_vrr = QCX*I_ERI_G3yz_S_S_S_vrr+WQX*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Px_S_vrr = QCX*I_ERI_G2y2z_S_S_S_vrr+WQX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Px_S_vrr = QCX*I_ERI_Gy3z_S_S_S_vrr+WQX*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Px_S_vrr = QCX*I_ERI_G4z_S_S_S_vrr+WQX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Py_S_vrr = QCY*I_ERI_G4x_S_S_S_vrr+WQY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Py_S_vrr = QCY*I_ERI_G3xy_S_S_S_vrr+WQY*I_ERI_G3xy_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Py_S_vrr = QCY*I_ERI_G3xz_S_S_S_vrr+WQY*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Py_S_vrr = QCY*I_ERI_G2x2y_S_S_S_vrr+WQY*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Py_S_vrr = QCY*I_ERI_G2xyz_S_S_S_vrr+WQY*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Py_S_vrr = QCY*I_ERI_G2x2z_S_S_S_vrr+WQY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Py_S_vrr = QCY*I_ERI_Gx3y_S_S_S_vrr+WQY*I_ERI_Gx3y_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Py_S_vrr = QCY*I_ERI_Gx2yz_S_S_S_vrr+WQY*I_ERI_Gx2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Py_S_vrr = QCY*I_ERI_Gxy2z_S_S_S_vrr+WQY*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Py_S_vrr = QCY*I_ERI_Gx3z_S_S_S_vrr+WQY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Py_S_vrr = QCY*I_ERI_G4y_S_S_S_vrr+WQY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Py_S_vrr = QCY*I_ERI_G3yz_S_S_S_vrr+WQY*I_ERI_G3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Py_S_vrr = QCY*I_ERI_G2y2z_S_S_S_vrr+WQY*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Py_S_vrr = QCY*I_ERI_Gy3z_S_S_S_vrr+WQY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Py_S_vrr = QCY*I_ERI_G4z_S_S_S_vrr+WQY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Pz_S_vrr = QCZ*I_ERI_G4x_S_S_S_vrr+WQZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Pz_S_vrr = QCZ*I_ERI_G3xy_S_S_S_vrr+WQZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Pz_S_vrr = QCZ*I_ERI_G3xz_S_S_S_vrr+WQZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Pz_S_vrr = QCZ*I_ERI_G2x2y_S_S_S_vrr+WQZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Pz_S_vrr = QCZ*I_ERI_G2xyz_S_S_S_vrr+WQZ*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Pz_S_vrr = QCZ*I_ERI_G2x2z_S_S_S_vrr+WQZ*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Pz_S_vrr = QCZ*I_ERI_Gx3y_S_S_S_vrr+WQZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_vrr = QCZ*I_ERI_Gx2yz_S_S_S_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_vrr = QCZ*I_ERI_Gxy2z_S_S_S_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Pz_S_vrr = QCZ*I_ERI_Gx3z_S_S_S_vrr+WQZ*I_ERI_Gx3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Pz_S_vrr = QCZ*I_ERI_G4y_S_S_S_vrr+WQZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Pz_S_vrr = QCZ*I_ERI_G3yz_S_S_S_vrr+WQZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Pz_S_vrr = QCZ*I_ERI_G2y2z_S_S_S_vrr+WQZ*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Pz_S_vrr = QCZ*I_ERI_Gy3z_S_S_S_vrr+WQZ*I_ERI_Gy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Pz_S_vrr = QCZ*I_ERI_G4z_S_S_S_vrr+WQZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2k*I_ERI_F3z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_D_S
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_vrr = QCX*I_ERI_F3x_S_D2x_S_vrr+WQX*I_ERI_F3x_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F3x_S_Px_S_vrr-2*rhod2esq*I_ERI_F3x_S_Px_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_F3x_S_vrr = QCX*I_ERI_F2xy_S_D2x_S_vrr+WQX*I_ERI_F2xy_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F2xy_S_Px_S_vrr-2*rhod2esq*I_ERI_F2xy_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_F3x_S_vrr = QCX*I_ERI_F2xz_S_D2x_S_vrr+WQX*I_ERI_F2xz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F2xz_S_Px_S_vrr-2*rhod2esq*I_ERI_F2xz_S_Px_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3x_S_vrr = QCX*I_ERI_Fx2y_S_D2x_S_vrr+WQX*I_ERI_Fx2y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Fx2y_S_Px_S_vrr-2*rhod2esq*I_ERI_Fx2y_S_Px_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3x_S_vrr = QCX*I_ERI_Fxyz_S_D2x_S_vrr+WQX*I_ERI_Fxyz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Fxyz_S_Px_S_vrr-2*rhod2esq*I_ERI_Fxyz_S_Px_S_M1_vrr+oned2k*I_ERI_Dyz_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3x_S_vrr = QCX*I_ERI_Fx2z_S_D2x_S_vrr+WQX*I_ERI_Fx2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Fx2z_S_Px_S_vrr-2*rhod2esq*I_ERI_Fx2z_S_Px_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_F3x_S_vrr = QCX*I_ERI_F3y_S_D2x_S_vrr+WQX*I_ERI_F3y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F3y_S_Px_S_vrr-2*rhod2esq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_F3x_S_vrr = QCX*I_ERI_F2yz_S_D2x_S_vrr+WQX*I_ERI_F2yz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F2yz_S_Px_S_vrr-2*rhod2esq*I_ERI_F2yz_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3x_S_vrr = QCX*I_ERI_Fy2z_S_D2x_S_vrr+WQX*I_ERI_Fy2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Fy2z_S_Px_S_vrr-2*rhod2esq*I_ERI_Fy2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_F3x_S_vrr = QCX*I_ERI_F3z_S_D2x_S_vrr+WQX*I_ERI_F3z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_F3z_S_Px_S_vrr-2*rhod2esq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_F2xy_S_vrr = QCY*I_ERI_F3x_S_D2x_S_vrr+WQY*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xy_S_vrr = QCY*I_ERI_F2xy_S_D2x_S_vrr+WQY*I_ERI_F2xy_S_D2x_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xy_S_vrr = QCY*I_ERI_F2xz_S_D2x_S_vrr+WQY*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_vrr = QCY*I_ERI_Fx2y_S_D2x_S_vrr+WQY*I_ERI_Fx2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xy_S_vrr = QCY*I_ERI_Fxyz_S_D2x_S_vrr+WQY*I_ERI_Fxyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Dxz_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_vrr = QCY*I_ERI_Fx2z_S_D2x_S_vrr+WQY*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_F2xy_S_vrr = QCY*I_ERI_F3y_S_D2x_S_vrr+WQY*I_ERI_F3y_S_D2x_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xy_S_vrr = QCY*I_ERI_F2yz_S_D2x_S_vrr+WQY*I_ERI_F2yz_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xy_S_vrr = QCY*I_ERI_Fy2z_S_D2x_S_vrr+WQY*I_ERI_Fy2z_S_D2x_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_F2xy_S_vrr = QCY*I_ERI_F3z_S_D2x_S_vrr+WQY*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_F2xz_S_vrr = QCZ*I_ERI_F3x_S_D2x_S_vrr+WQZ*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xz_S_vrr = QCZ*I_ERI_F2xy_S_D2x_S_vrr+WQZ*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xz_S_vrr = QCZ*I_ERI_F2xz_S_D2x_S_vrr+WQZ*I_ERI_F2xz_S_D2x_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_vrr = QCZ*I_ERI_Fx2y_S_D2x_S_vrr+WQZ*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xz_S_vrr = QCZ*I_ERI_Fxyz_S_D2x_S_vrr+WQZ*I_ERI_Fxyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_vrr = QCZ*I_ERI_Fx2z_S_D2x_S_vrr+WQZ*I_ERI_Fx2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_F2xz_S_vrr = QCZ*I_ERI_F3y_S_D2x_S_vrr+WQZ*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xz_S_vrr = QCZ*I_ERI_F2yz_S_D2x_S_vrr+WQZ*I_ERI_F2yz_S_D2x_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xz_S_vrr = QCZ*I_ERI_Fy2z_S_D2x_S_vrr+WQZ*I_ERI_Fy2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_F2xz_S_vrr = QCZ*I_ERI_F3z_S_D2x_S_vrr+WQZ*I_ERI_F3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2y_S_vrr = QCX*I_ERI_F3x_S_D2y_S_vrr+WQX*I_ERI_F3x_S_D2y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_vrr = QCX*I_ERI_F2xy_S_D2y_S_vrr+WQX*I_ERI_F2xy_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_vrr = QCX*I_ERI_F2xz_S_D2y_S_vrr+WQX*I_ERI_F2xz_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_vrr = QCX*I_ERI_Fx2y_S_D2y_S_vrr+WQX*I_ERI_Fx2y_S_D2y_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2y_S_vrr = QCX*I_ERI_Fxyz_S_D2y_S_vrr+WQX*I_ERI_Fxyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Dyz_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_vrr = QCX*I_ERI_Fx2z_S_D2y_S_vrr+WQX*I_ERI_Fx2z_S_D2y_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2y_S_vrr = QCX*I_ERI_F3y_S_D2y_S_vrr+WQX*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_vrr = QCX*I_ERI_F2yz_S_D2y_S_vrr+WQX*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2y_S_vrr = QCX*I_ERI_Fy2z_S_D2y_S_vrr+WQX*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2y_S_vrr = QCX*I_ERI_F3z_S_D2y_S_vrr+WQX*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fxyz_S_vrr = QCZ*I_ERI_F3x_S_Dxy_S_vrr+WQZ*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_vrr = QCZ*I_ERI_F2xy_S_Dxy_S_vrr+WQZ*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_vrr = QCZ*I_ERI_F2xz_S_Dxy_S_vrr+WQZ*I_ERI_F2xz_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_vrr = QCZ*I_ERI_Fx2y_S_Dxy_S_vrr+WQZ*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fxyz_S_vrr = QCZ*I_ERI_Fxyz_S_Dxy_S_vrr+WQZ*I_ERI_Fxyz_S_Dxy_S_M1_vrr+oned2k*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_vrr = QCZ*I_ERI_Fx2z_S_Dxy_S_vrr+WQZ*I_ERI_Fx2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Dxy_S_M1_vrr;
      Double I_ERI_F3y_S_Fxyz_S_vrr = QCZ*I_ERI_F3y_S_Dxy_S_vrr+WQZ*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_vrr = QCZ*I_ERI_F2yz_S_Dxy_S_vrr+WQZ*I_ERI_F2yz_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fxyz_S_vrr = QCZ*I_ERI_Fy2z_S_Dxy_S_vrr+WQZ*I_ERI_Fy2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Dxy_S_M1_vrr;
      Double I_ERI_F3z_S_Fxyz_S_vrr = QCZ*I_ERI_F3z_S_Dxy_S_vrr+WQZ*I_ERI_F3z_S_Dxy_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2z_S_vrr = QCX*I_ERI_F3x_S_D2z_S_vrr+WQX*I_ERI_F3x_S_D2z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_vrr = QCX*I_ERI_F2xy_S_D2z_S_vrr+WQX*I_ERI_F2xy_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_vrr = QCX*I_ERI_F2xz_S_D2z_S_vrr+WQX*I_ERI_F2xz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_vrr = QCX*I_ERI_Fx2y_S_D2z_S_vrr+WQX*I_ERI_Fx2y_S_D2z_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2z_S_vrr = QCX*I_ERI_Fxyz_S_D2z_S_vrr+WQX*I_ERI_Fxyz_S_D2z_S_M1_vrr+oned2k*I_ERI_Dyz_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_vrr = QCX*I_ERI_Fx2z_S_D2z_S_vrr+WQX*I_ERI_Fx2z_S_D2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2z_S_vrr = QCX*I_ERI_F3y_S_D2z_S_vrr+WQX*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_vrr = QCX*I_ERI_F2yz_S_D2z_S_vrr+WQX*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2z_S_vrr = QCX*I_ERI_Fy2z_S_D2z_S_vrr+WQX*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2z_S_vrr = QCX*I_ERI_F3z_S_D2z_S_vrr+WQX*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_F3x_S_F3y_S_vrr = QCY*I_ERI_F3x_S_D2y_S_vrr+WQY*I_ERI_F3x_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F3x_S_Py_S_vrr-2*rhod2esq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_F3y_S_vrr = QCY*I_ERI_F2xy_S_D2y_S_vrr+WQY*I_ERI_F2xy_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F2xy_S_Py_S_vrr-2*rhod2esq*I_ERI_F2xy_S_Py_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_F3y_S_vrr = QCY*I_ERI_F2xz_S_D2y_S_vrr+WQY*I_ERI_F2xz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F2xz_S_Py_S_vrr-2*rhod2esq*I_ERI_F2xz_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3y_S_vrr = QCY*I_ERI_Fx2y_S_D2y_S_vrr+WQY*I_ERI_Fx2y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Fx2y_S_Py_S_vrr-2*rhod2esq*I_ERI_Fx2y_S_Py_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3y_S_vrr = QCY*I_ERI_Fxyz_S_D2y_S_vrr+WQY*I_ERI_Fxyz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Fxyz_S_Py_S_vrr-2*rhod2esq*I_ERI_Fxyz_S_Py_S_M1_vrr+oned2k*I_ERI_Dxz_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3y_S_vrr = QCY*I_ERI_Fx2z_S_D2y_S_vrr+WQY*I_ERI_Fx2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Fx2z_S_Py_S_vrr-2*rhod2esq*I_ERI_Fx2z_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_F3y_S_vrr = QCY*I_ERI_F3y_S_D2y_S_vrr+WQY*I_ERI_F3y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F3y_S_Py_S_vrr-2*rhod2esq*I_ERI_F3y_S_Py_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_F3y_S_vrr = QCY*I_ERI_F2yz_S_D2y_S_vrr+WQY*I_ERI_F2yz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F2yz_S_Py_S_vrr-2*rhod2esq*I_ERI_F2yz_S_Py_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3y_S_vrr = QCY*I_ERI_Fy2z_S_D2y_S_vrr+WQY*I_ERI_Fy2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Fy2z_S_Py_S_vrr-2*rhod2esq*I_ERI_Fy2z_S_Py_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3z_S_F3y_S_vrr = QCY*I_ERI_F3z_S_D2y_S_vrr+WQY*I_ERI_F3z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_F3z_S_Py_S_vrr-2*rhod2esq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_F2yz_S_vrr = QCZ*I_ERI_F3x_S_D2y_S_vrr+WQZ*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_F2yz_S_vrr = QCZ*I_ERI_F2xy_S_D2y_S_vrr+WQZ*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_F2yz_S_vrr = QCZ*I_ERI_F2xz_S_D2y_S_vrr+WQZ*I_ERI_F2xz_S_D2y_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_vrr = QCZ*I_ERI_Fx2y_S_D2y_S_vrr+WQZ*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2yz_S_vrr = QCZ*I_ERI_Fxyz_S_D2y_S_vrr+WQZ*I_ERI_Fxyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_vrr = QCZ*I_ERI_Fx2z_S_D2y_S_vrr+WQZ*I_ERI_Fx2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_F2yz_S_vrr = QCZ*I_ERI_F3y_S_D2y_S_vrr+WQZ*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_F2yz_S_vrr = QCZ*I_ERI_F2yz_S_D2y_S_vrr+WQZ*I_ERI_F2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2yz_S_vrr = QCZ*I_ERI_Fy2z_S_D2y_S_vrr+WQZ*I_ERI_Fy2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2y_S_M1_vrr;
      Double I_ERI_F3z_S_F2yz_S_vrr = QCZ*I_ERI_F3z_S_D2y_S_vrr+WQZ*I_ERI_F3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fy2z_S_vrr = QCY*I_ERI_F3x_S_D2z_S_vrr+WQY*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_vrr = QCY*I_ERI_F2xy_S_D2z_S_vrr+WQY*I_ERI_F2xy_S_D2z_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_vrr = QCY*I_ERI_F2xz_S_D2z_S_vrr+WQY*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_vrr = QCY*I_ERI_Fx2y_S_D2z_S_vrr+WQY*I_ERI_Fx2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fy2z_S_vrr = QCY*I_ERI_Fxyz_S_D2z_S_vrr+WQY*I_ERI_Fxyz_S_D2z_S_M1_vrr+oned2k*I_ERI_Dxz_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_vrr = QCY*I_ERI_Fx2z_S_D2z_S_vrr+WQY*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fy2z_S_vrr = QCY*I_ERI_F3y_S_D2z_S_vrr+WQY*I_ERI_F3y_S_D2z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_vrr = QCY*I_ERI_F2yz_S_D2z_S_vrr+WQY*I_ERI_F2yz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fy2z_S_vrr = QCY*I_ERI_Fy2z_S_D2z_S_vrr+WQY*I_ERI_Fy2z_S_D2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fy2z_S_vrr = QCY*I_ERI_F3z_S_D2z_S_vrr+WQY*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_F3x_S_F3z_S_vrr = QCZ*I_ERI_F3x_S_D2z_S_vrr+WQZ*I_ERI_F3x_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F3x_S_Pz_S_vrr-2*rhod2esq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xy_S_F3z_S_vrr = QCZ*I_ERI_F2xy_S_D2z_S_vrr+WQZ*I_ERI_F2xy_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F2xy_S_Pz_S_vrr-2*rhod2esq*I_ERI_F2xy_S_Pz_S_M1_vrr;
      Double I_ERI_F2xz_S_F3z_S_vrr = QCZ*I_ERI_F2xz_S_D2z_S_vrr+WQZ*I_ERI_F2xz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F2xz_S_Pz_S_vrr-2*rhod2esq*I_ERI_F2xz_S_Pz_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3z_S_vrr = QCZ*I_ERI_Fx2y_S_D2z_S_vrr+WQZ*I_ERI_Fx2y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Fx2y_S_Pz_S_vrr-2*rhod2esq*I_ERI_Fx2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3z_S_vrr = QCZ*I_ERI_Fxyz_S_D2z_S_vrr+WQZ*I_ERI_Fxyz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Fxyz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Fxyz_S_Pz_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3z_S_vrr = QCZ*I_ERI_Fx2z_S_D2z_S_vrr+WQZ*I_ERI_Fx2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Fx2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Fx2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_F3z_S_vrr = QCZ*I_ERI_F3y_S_D2z_S_vrr+WQZ*I_ERI_F3y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F3y_S_Pz_S_vrr-2*rhod2esq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_F2yz_S_F3z_S_vrr = QCZ*I_ERI_F2yz_S_D2z_S_vrr+WQZ*I_ERI_F2yz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F2yz_S_Pz_S_vrr-2*rhod2esq*I_ERI_F2yz_S_Pz_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3z_S_vrr = QCZ*I_ERI_Fy2z_S_D2z_S_vrr+WQZ*I_ERI_Fy2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Fy2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Fy2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_F3z_S_vrr = QCZ*I_ERI_F3z_S_D2z_S_vrr+WQZ*I_ERI_F3z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_F3z_S_Pz_S_vrr-2*rhod2esq*I_ERI_F3z_S_Pz_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_vrr = QCX*I_ERI_G4x_S_Px_S_vrr+WQX*I_ERI_G4x_S_Px_S_M1_vrr+oned2e*I_ERI_G4x_S_S_S_vrr-rhod2esq*I_ERI_G4x_S_S_S_M1_vrr+4*oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G3xy_S_D2x_S_vrr = QCX*I_ERI_G3xy_S_Px_S_vrr+WQX*I_ERI_G3xy_S_Px_S_M1_vrr+oned2e*I_ERI_G3xy_S_S_S_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G3xz_S_D2x_S_vrr = QCX*I_ERI_G3xz_S_Px_S_vrr+WQX*I_ERI_G3xz_S_Px_S_M1_vrr+oned2e*I_ERI_G3xz_S_S_S_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_Px_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2x_S_vrr = QCX*I_ERI_G2x2y_S_Px_S_vrr+WQX*I_ERI_G2x2y_S_Px_S_M1_vrr+oned2e*I_ERI_G2x2y_S_S_S_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2x_S_vrr = QCX*I_ERI_G2xyz_S_Px_S_vrr+WQX*I_ERI_G2xyz_S_Px_S_M1_vrr+oned2e*I_ERI_G2xyz_S_S_S_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2x_S_vrr = QCX*I_ERI_G2x2z_S_Px_S_vrr+WQX*I_ERI_G2x2z_S_Px_S_M1_vrr+oned2e*I_ERI_G2x2z_S_S_S_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2x_S_vrr = QCX*I_ERI_Gx3y_S_Px_S_vrr+WQX*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2e*I_ERI_Gx3y_S_S_S_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2x_S_vrr = QCX*I_ERI_Gx2yz_S_Px_S_vrr+WQX*I_ERI_Gx2yz_S_Px_S_M1_vrr+oned2e*I_ERI_Gx2yz_S_S_S_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_F2yz_S_Px_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2x_S_vrr = QCX*I_ERI_Gxy2z_S_Px_S_vrr+WQX*I_ERI_Gxy2z_S_Px_S_M1_vrr+oned2e*I_ERI_Gxy2z_S_S_S_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fy2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2x_S_vrr = QCX*I_ERI_Gx3z_S_Px_S_vrr+WQX*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2e*I_ERI_Gx3z_S_S_S_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4y_S_D2x_S_vrr = QCX*I_ERI_G4y_S_Px_S_vrr+WQX*I_ERI_G4y_S_Px_S_M1_vrr+oned2e*I_ERI_G4y_S_S_S_vrr-rhod2esq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_D2x_S_vrr = QCX*I_ERI_G3yz_S_Px_S_vrr+WQX*I_ERI_G3yz_S_Px_S_M1_vrr+oned2e*I_ERI_G3yz_S_S_S_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2x_S_vrr = QCX*I_ERI_G2y2z_S_Px_S_vrr+WQX*I_ERI_G2y2z_S_Px_S_M1_vrr+oned2e*I_ERI_G2y2z_S_S_S_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2x_S_vrr = QCX*I_ERI_Gy3z_S_Px_S_vrr+WQX*I_ERI_Gy3z_S_Px_S_M1_vrr+oned2e*I_ERI_Gy3z_S_S_S_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_D2x_S_vrr = QCX*I_ERI_G4z_S_Px_S_vrr+WQX*I_ERI_G4z_S_Px_S_M1_vrr+oned2e*I_ERI_G4z_S_S_S_vrr-rhod2esq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Dxy_S_vrr = QCY*I_ERI_G4x_S_Px_S_vrr+WQY*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_G3xy_S_Dxy_S_vrr = QCY*I_ERI_G3xy_S_Px_S_vrr+WQY*I_ERI_G3xy_S_Px_S_M1_vrr+oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G3xz_S_Dxy_S_vrr = QCY*I_ERI_G3xz_S_Px_S_vrr+WQY*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_vrr = QCY*I_ERI_G2x2y_S_Px_S_vrr+WQY*I_ERI_G2x2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dxy_S_vrr = QCY*I_ERI_G2xyz_S_Px_S_vrr+WQY*I_ERI_G2xyz_S_Px_S_M1_vrr+oned2k*I_ERI_F2xz_S_Px_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dxy_S_vrr = QCY*I_ERI_G2x2z_S_Px_S_vrr+WQY*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_vrr = QCY*I_ERI_Gx3y_S_Px_S_vrr+WQY*I_ERI_Gx3y_S_Px_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dxy_S_vrr = QCY*I_ERI_Gx2yz_S_Px_S_vrr+WQY*I_ERI_Gx2yz_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dxy_S_vrr = QCY*I_ERI_Gxy2z_S_Px_S_vrr+WQY*I_ERI_Gxy2z_S_Px_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_vrr = QCY*I_ERI_Gx3z_S_Px_S_vrr+WQY*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_G4y_S_Dxy_S_vrr = QCY*I_ERI_G4y_S_Px_S_vrr+WQY*I_ERI_G4y_S_Px_S_M1_vrr+4*oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_Dxy_S_vrr = QCY*I_ERI_G3yz_S_Px_S_vrr+WQY*I_ERI_G3yz_S_Px_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dxy_S_vrr = QCY*I_ERI_G2y2z_S_Px_S_vrr+WQY*I_ERI_G2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_vrr = QCY*I_ERI_Gy3z_S_Px_S_vrr+WQY*I_ERI_Gy3z_S_Px_S_M1_vrr+oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_Dxy_S_vrr = QCY*I_ERI_G4z_S_Px_S_vrr+WQY*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_Dxz_S_vrr = QCZ*I_ERI_G4x_S_Px_S_vrr+WQZ*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_G3xy_S_Dxz_S_vrr = QCZ*I_ERI_G3xy_S_Px_S_vrr+WQZ*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_G3xz_S_Dxz_S_vrr = QCZ*I_ERI_G3xz_S_Px_S_vrr+WQZ*I_ERI_G3xz_S_Px_S_M1_vrr+oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dxz_S_vrr = QCZ*I_ERI_G2x2y_S_Px_S_vrr+WQZ*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dxz_S_vrr = QCZ*I_ERI_G2xyz_S_Px_S_vrr+WQZ*I_ERI_G2xyz_S_Px_S_M1_vrr+oned2k*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dxz_S_vrr = QCZ*I_ERI_G2x2z_S_Px_S_vrr+WQZ*I_ERI_G2x2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Px_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dxz_S_vrr = QCZ*I_ERI_Gx3y_S_Px_S_vrr+WQZ*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dxz_S_vrr = QCZ*I_ERI_Gx2yz_S_Px_S_vrr+WQZ*I_ERI_Gx2yz_S_Px_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dxz_S_vrr = QCZ*I_ERI_Gxy2z_S_Px_S_vrr+WQZ*I_ERI_Gxy2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dxz_S_vrr = QCZ*I_ERI_Gx3z_S_Px_S_vrr+WQZ*I_ERI_Gx3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_G4y_S_Dxz_S_vrr = QCZ*I_ERI_G4y_S_Px_S_vrr+WQZ*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_Dxz_S_vrr = QCZ*I_ERI_G3yz_S_Px_S_vrr+WQZ*I_ERI_G3yz_S_Px_S_M1_vrr+oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dxz_S_vrr = QCZ*I_ERI_G2y2z_S_Px_S_vrr+WQZ*I_ERI_G2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dxz_S_vrr = QCZ*I_ERI_Gy3z_S_Px_S_vrr+WQZ*I_ERI_Gy3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_Dxz_S_vrr = QCZ*I_ERI_G4z_S_Px_S_vrr+WQZ*I_ERI_G4z_S_Px_S_M1_vrr+4*oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_D2y_S_vrr = QCY*I_ERI_G4x_S_Py_S_vrr+WQY*I_ERI_G4x_S_Py_S_M1_vrr+oned2e*I_ERI_G4x_S_S_S_vrr-rhod2esq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_D2y_S_vrr = QCY*I_ERI_G3xy_S_Py_S_vrr+WQY*I_ERI_G3xy_S_Py_S_M1_vrr+oned2e*I_ERI_G3xy_S_S_S_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G3xz_S_D2y_S_vrr = QCY*I_ERI_G3xz_S_Py_S_vrr+WQY*I_ERI_G3xz_S_Py_S_M1_vrr+oned2e*I_ERI_G3xz_S_S_S_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2y_S_vrr = QCY*I_ERI_G2x2y_S_Py_S_vrr+WQY*I_ERI_G2x2y_S_Py_S_M1_vrr+oned2e*I_ERI_G2x2y_S_S_S_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2y_S_vrr = QCY*I_ERI_G2xyz_S_Py_S_vrr+WQY*I_ERI_G2xyz_S_Py_S_M1_vrr+oned2e*I_ERI_G2xyz_S_S_S_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xz_S_Py_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2y_S_vrr = QCY*I_ERI_G2x2z_S_Py_S_vrr+WQY*I_ERI_G2x2z_S_Py_S_M1_vrr+oned2e*I_ERI_G2x2z_S_S_S_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2y_S_vrr = QCY*I_ERI_Gx3y_S_Py_S_vrr+WQY*I_ERI_Gx3y_S_Py_S_M1_vrr+oned2e*I_ERI_Gx3y_S_S_S_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Py_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2y_S_vrr = QCY*I_ERI_Gx2yz_S_Py_S_vrr+WQY*I_ERI_Gx2yz_S_Py_S_M1_vrr+oned2e*I_ERI_Gx2yz_S_S_S_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Py_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2y_S_vrr = QCY*I_ERI_Gxy2z_S_Py_S_vrr+WQY*I_ERI_Gxy2z_S_Py_S_M1_vrr+oned2e*I_ERI_Gxy2z_S_S_S_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Py_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2y_S_vrr = QCY*I_ERI_Gx3z_S_Py_S_vrr+WQY*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2e*I_ERI_Gx3z_S_S_S_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_D2y_S_vrr = QCY*I_ERI_G4y_S_Py_S_vrr+WQY*I_ERI_G4y_S_Py_S_M1_vrr+oned2e*I_ERI_G4y_S_S_S_vrr-rhod2esq*I_ERI_G4y_S_S_S_M1_vrr+4*oned2k*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_G3yz_S_D2y_S_vrr = QCY*I_ERI_G3yz_S_Py_S_vrr+WQY*I_ERI_G3yz_S_Py_S_M1_vrr+oned2e*I_ERI_G3yz_S_S_S_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Py_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2y_S_vrr = QCY*I_ERI_G2y2z_S_Py_S_vrr+WQY*I_ERI_G2y2z_S_Py_S_M1_vrr+oned2e*I_ERI_G2y2z_S_S_S_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_Py_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2y_S_vrr = QCY*I_ERI_Gy3z_S_Py_S_vrr+WQY*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2e*I_ERI_Gy3z_S_S_S_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4z_S_D2y_S_vrr = QCY*I_ERI_G4z_S_Py_S_vrr+WQY*I_ERI_G4z_S_Py_S_M1_vrr+oned2e*I_ERI_G4z_S_S_S_vrr-rhod2esq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Dyz_S_vrr = QCZ*I_ERI_G4x_S_Py_S_vrr+WQZ*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_G3xy_S_Dyz_S_vrr = QCZ*I_ERI_G3xy_S_Py_S_vrr+WQZ*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_G3xz_S_Dyz_S_vrr = QCZ*I_ERI_G3xz_S_Py_S_vrr+WQZ*I_ERI_G3xz_S_Py_S_M1_vrr+oned2k*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dyz_S_vrr = QCZ*I_ERI_G2x2y_S_Py_S_vrr+WQZ*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dyz_S_vrr = QCZ*I_ERI_G2xyz_S_Py_S_vrr+WQZ*I_ERI_G2xyz_S_Py_S_M1_vrr+oned2k*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dyz_S_vrr = QCZ*I_ERI_G2x2z_S_Py_S_vrr+WQZ*I_ERI_G2x2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Py_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dyz_S_vrr = QCZ*I_ERI_Gx3y_S_Py_S_vrr+WQZ*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dyz_S_vrr = QCZ*I_ERI_Gx2yz_S_Py_S_vrr+WQZ*I_ERI_Gx2yz_S_Py_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Py_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dyz_S_vrr = QCZ*I_ERI_Gxy2z_S_Py_S_vrr+WQZ*I_ERI_Gxy2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Py_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dyz_S_vrr = QCZ*I_ERI_Gx3z_S_Py_S_vrr+WQZ*I_ERI_Gx3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Py_S_M1_vrr;
      Double I_ERI_G4y_S_Dyz_S_vrr = QCZ*I_ERI_G4y_S_Py_S_vrr+WQZ*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_G3yz_S_Dyz_S_vrr = QCZ*I_ERI_G3yz_S_Py_S_vrr+WQZ*I_ERI_G3yz_S_Py_S_M1_vrr+oned2k*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dyz_S_vrr = QCZ*I_ERI_G2y2z_S_Py_S_vrr+WQZ*I_ERI_G2y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Py_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dyz_S_vrr = QCZ*I_ERI_Gy3z_S_Py_S_vrr+WQZ*I_ERI_Gy3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_Py_S_M1_vrr;
      Double I_ERI_G4z_S_Dyz_S_vrr = QCZ*I_ERI_G4z_S_Py_S_vrr+WQZ*I_ERI_G4z_S_Py_S_M1_vrr+4*oned2k*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4x_S_D2z_S_vrr = QCZ*I_ERI_G4x_S_Pz_S_vrr+WQZ*I_ERI_G4x_S_Pz_S_M1_vrr+oned2e*I_ERI_G4x_S_S_S_vrr-rhod2esq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_D2z_S_vrr = QCZ*I_ERI_G3xy_S_Pz_S_vrr+WQZ*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2e*I_ERI_G3xy_S_S_S_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_D2z_S_vrr = QCZ*I_ERI_G3xz_S_Pz_S_vrr+WQZ*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2e*I_ERI_G3xz_S_S_S_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2z_S_vrr = QCZ*I_ERI_G2x2y_S_Pz_S_vrr+WQZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+oned2e*I_ERI_G2x2y_S_S_S_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2z_S_vrr = QCZ*I_ERI_G2xyz_S_Pz_S_vrr+WQZ*I_ERI_G2xyz_S_Pz_S_M1_vrr+oned2e*I_ERI_G2xyz_S_S_S_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2z_S_vrr = QCZ*I_ERI_G2x2z_S_Pz_S_vrr+WQZ*I_ERI_G2x2z_S_Pz_S_M1_vrr+oned2e*I_ERI_G2x2z_S_S_S_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Pz_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2z_S_vrr = QCZ*I_ERI_Gx3y_S_Pz_S_vrr+WQZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2e*I_ERI_Gx3y_S_S_S_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2z_S_vrr = QCZ*I_ERI_Gx2yz_S_Pz_S_vrr+WQZ*I_ERI_Gx2yz_S_Pz_S_M1_vrr+oned2e*I_ERI_Gx2yz_S_S_S_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Pz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2z_S_vrr = QCZ*I_ERI_Gxy2z_S_Pz_S_vrr+WQZ*I_ERI_Gxy2z_S_Pz_S_M1_vrr+oned2e*I_ERI_Gxy2z_S_S_S_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Pz_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2z_S_vrr = QCZ*I_ERI_Gx3z_S_Pz_S_vrr+WQZ*I_ERI_Gx3z_S_Pz_S_M1_vrr+oned2e*I_ERI_Gx3z_S_S_S_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Pz_S_M1_vrr;
      Double I_ERI_G4y_S_D2z_S_vrr = QCZ*I_ERI_G4y_S_Pz_S_vrr+WQZ*I_ERI_G4y_S_Pz_S_M1_vrr+oned2e*I_ERI_G4y_S_S_S_vrr-rhod2esq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_D2z_S_vrr = QCZ*I_ERI_G3yz_S_Pz_S_vrr+WQZ*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2e*I_ERI_G3yz_S_S_S_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2z_S_vrr = QCZ*I_ERI_G2y2z_S_Pz_S_vrr+WQZ*I_ERI_G2y2z_S_Pz_S_M1_vrr+oned2e*I_ERI_G2y2z_S_S_S_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Pz_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2z_S_vrr = QCZ*I_ERI_Gy3z_S_Pz_S_vrr+WQZ*I_ERI_Gy3z_S_Pz_S_M1_vrr+oned2e*I_ERI_Gy3z_S_S_S_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_Pz_S_M1_vrr;
      Double I_ERI_G4z_S_D2z_S_vrr = QCZ*I_ERI_G4z_S_Pz_S_vrr+WQZ*I_ERI_G4z_S_Pz_S_M1_vrr+oned2e*I_ERI_G4z_S_S_S_vrr-rhod2esq*I_ERI_G4z_S_S_S_M1_vrr+4*oned2k*I_ERI_F3z_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_vrr = PAX*I_ERI_G4x_S_Px_S_vrr+WPX*I_ERI_G4x_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Px_S_vrr = PAY*I_ERI_G4x_S_Px_S_vrr+WPY*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Px_S_vrr = PAZ*I_ERI_G4x_S_Px_S_vrr+WPZ*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Px_S_vrr = PAY*I_ERI_G3xy_S_Px_S_vrr+WPY*I_ERI_G3xy_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Px_S_vrr = PAZ*I_ERI_G3xy_S_Px_S_vrr+WPZ*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Px_S_vrr = PAZ*I_ERI_G3xz_S_Px_S_vrr+WPZ*I_ERI_G3xz_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Px_S_vrr = PAX*I_ERI_Gx3y_S_Px_S_vrr+WPX*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Px_S_vrr = PAZ*I_ERI_G2x2y_S_Px_S_vrr+WPZ*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Px_S_vrr = PAY*I_ERI_G2x2z_S_Px_S_vrr+WPY*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Px_S_vrr = PAX*I_ERI_Gx3z_S_Px_S_vrr+WPX*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Px_S_vrr = PAX*I_ERI_G4y_S_Px_S_vrr+WPX*I_ERI_G4y_S_Px_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Px_S_vrr = PAZ*I_ERI_Gx3y_S_Px_S_vrr+WPZ*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_vrr = PAX*I_ERI_G2y2z_S_Px_S_vrr+WPX*I_ERI_G2y2z_S_Px_S_M1_vrr+oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Px_S_vrr = PAY*I_ERI_Gx3z_S_Px_S_vrr+WPY*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Px_S_vrr = PAX*I_ERI_G4z_S_Px_S_vrr+WPX*I_ERI_G4z_S_Px_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Px_S_vrr = PAY*I_ERI_G4y_S_Px_S_vrr+WPY*I_ERI_G4y_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Px_S_vrr = PAZ*I_ERI_G4y_S_Px_S_vrr+WPZ*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Px_S_vrr = PAZ*I_ERI_G3yz_S_Px_S_vrr+WPZ*I_ERI_G3yz_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Px_S_vrr = PAY*I_ERI_Gy3z_S_Px_S_vrr+WPY*I_ERI_Gy3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Px_S_vrr = PAY*I_ERI_G4z_S_Px_S_vrr+WPY*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Px_S_vrr = PAZ*I_ERI_G4z_S_Px_S_vrr+WPZ*I_ERI_G4z_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_Py_S_vrr = PAX*I_ERI_G4x_S_Py_S_vrr+WPX*I_ERI_G4x_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Py_S_vrr = PAY*I_ERI_G4x_S_Py_S_vrr+WPY*I_ERI_G4x_S_Py_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Py_S_vrr = PAZ*I_ERI_G4x_S_Py_S_vrr+WPZ*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Py_S_vrr = PAY*I_ERI_G3xy_S_Py_S_vrr+WPY*I_ERI_G3xy_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Py_S_vrr = PAZ*I_ERI_G3xy_S_Py_S_vrr+WPZ*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Py_S_vrr = PAZ*I_ERI_G3xz_S_Py_S_vrr+WPZ*I_ERI_G3xz_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Py_S_vrr = PAX*I_ERI_Gx3y_S_Py_S_vrr+WPX*I_ERI_Gx3y_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Py_S_vrr = PAZ*I_ERI_G2x2y_S_Py_S_vrr+WPZ*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Py_S_vrr = PAY*I_ERI_G2x2z_S_Py_S_vrr+WPY*I_ERI_G2x2z_S_Py_S_M1_vrr+oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Py_S_vrr = PAX*I_ERI_Gx3z_S_Py_S_vrr+WPX*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Py_S_vrr = PAX*I_ERI_G4y_S_Py_S_vrr+WPX*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Py_S_vrr = PAZ*I_ERI_Gx3y_S_Py_S_vrr+WPZ*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_vrr = PAX*I_ERI_G2y2z_S_Py_S_vrr+WPX*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Py_S_vrr = PAY*I_ERI_Gx3z_S_Py_S_vrr+WPY*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Py_S_vrr = PAX*I_ERI_G4z_S_Py_S_vrr+WPX*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Py_S_vrr = PAY*I_ERI_G4y_S_Py_S_vrr+WPY*I_ERI_G4y_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Py_S_vrr = PAZ*I_ERI_G4y_S_Py_S_vrr+WPZ*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Py_S_vrr = PAZ*I_ERI_G3yz_S_Py_S_vrr+WPZ*I_ERI_G3yz_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Py_S_vrr = PAY*I_ERI_Gy3z_S_Py_S_vrr+WPY*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Py_S_vrr = PAY*I_ERI_G4z_S_Py_S_vrr+WPY*I_ERI_G4z_S_Py_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Py_S_vrr = PAZ*I_ERI_G4z_S_Py_S_vrr+WPZ*I_ERI_G4z_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_Pz_S_vrr = PAX*I_ERI_G4x_S_Pz_S_vrr+WPX*I_ERI_G4x_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xy_S_Pz_S_vrr = PAY*I_ERI_G4x_S_Pz_S_vrr+WPY*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xz_S_Pz_S_vrr = PAZ*I_ERI_G4x_S_Pz_S_vrr+WPZ*I_ERI_G4x_S_Pz_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Pz_S_vrr = PAY*I_ERI_G3xy_S_Pz_S_vrr+WPY*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Pz_S_vrr = PAZ*I_ERI_G3xy_S_Pz_S_vrr+WPZ*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Pz_S_vrr = PAZ*I_ERI_G3xz_S_Pz_S_vrr+WPZ*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Pz_S_vrr = PAX*I_ERI_Gx3y_S_Pz_S_vrr+WPX*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_vrr = PAZ*I_ERI_G2x2y_S_Pz_S_vrr+WPZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_vrr = PAY*I_ERI_G2x2z_S_Pz_S_vrr+WPY*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Pz_S_vrr = PAX*I_ERI_Gx3z_S_Pz_S_vrr+WPX*I_ERI_Gx3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Pz_S_vrr = PAX*I_ERI_G4y_S_Pz_S_vrr+WPX*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_vrr = PAZ*I_ERI_Gx3y_S_Pz_S_vrr+WPZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_vrr = PAX*I_ERI_G2y2z_S_Pz_S_vrr+WPX*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_vrr = PAY*I_ERI_Gx3z_S_Pz_S_vrr+WPY*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Pz_S_vrr = PAX*I_ERI_G4z_S_Pz_S_vrr+WPX*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_Pz_S_vrr = PAY*I_ERI_G4y_S_Pz_S_vrr+WPY*I_ERI_G4y_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H4yz_S_Pz_S_vrr = PAZ*I_ERI_G4y_S_Pz_S_vrr+WPZ*I_ERI_G4y_S_Pz_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Pz_S_vrr = PAZ*I_ERI_G3yz_S_Pz_S_vrr+WPZ*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr+oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Pz_S_vrr = PAY*I_ERI_Gy3z_S_Pz_S_vrr+WPY*I_ERI_Gy3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Pz_S_vrr = PAY*I_ERI_G4z_S_Pz_S_vrr+WPY*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_Pz_S_vrr = PAZ*I_ERI_G4z_S_Pz_S_vrr+WPZ*I_ERI_G4z_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_b_coefs = beta;
      I_ERI_F3x_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_b += SQ_ERI_F_S_P_S_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_cc_coefs = gamma*gamma;
      I_ERI_F3x_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_c_coefs = gamma;
      I_ERI_F3x_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_P_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_P_d_coefs = delta;
      I_ERI_F3x_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3x_S_S_Px_vrr;
      I_ERI_F2xy_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xy_S_S_Px_vrr;
      I_ERI_F2xz_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xz_S_S_Px_vrr;
      I_ERI_Fx2y_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2y_S_S_Px_vrr;
      I_ERI_Fxyz_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fxyz_S_S_Px_vrr;
      I_ERI_Fx2z_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2z_S_S_Px_vrr;
      I_ERI_F3y_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3y_S_S_Px_vrr;
      I_ERI_F2yz_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2yz_S_S_Px_vrr;
      I_ERI_Fy2z_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fy2z_S_S_Px_vrr;
      I_ERI_F3z_S_S_Px_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3z_S_S_Px_vrr;
      I_ERI_F3x_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3x_S_S_Py_vrr;
      I_ERI_F2xy_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xy_S_S_Py_vrr;
      I_ERI_F2xz_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xz_S_S_Py_vrr;
      I_ERI_Fx2y_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2y_S_S_Py_vrr;
      I_ERI_Fxyz_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fxyz_S_S_Py_vrr;
      I_ERI_Fx2z_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2z_S_S_Py_vrr;
      I_ERI_F3y_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3y_S_S_Py_vrr;
      I_ERI_F2yz_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2yz_S_S_Py_vrr;
      I_ERI_Fy2z_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fy2z_S_S_Py_vrr;
      I_ERI_F3z_S_S_Py_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3z_S_S_Py_vrr;
      I_ERI_F3x_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3x_S_S_Pz_vrr;
      I_ERI_F2xy_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xy_S_S_Pz_vrr;
      I_ERI_F2xz_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2xz_S_S_Pz_vrr;
      I_ERI_Fx2y_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2y_S_S_Pz_vrr;
      I_ERI_Fxyz_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fxyz_S_S_Pz_vrr;
      I_ERI_Fx2z_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fx2z_S_S_Pz_vrr;
      I_ERI_F3y_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3y_S_S_Pz_vrr;
      I_ERI_F2yz_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F2yz_S_S_Pz_vrr;
      I_ERI_Fy2z_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_Fy2z_S_S_Pz_vrr;
      I_ERI_F3z_S_S_Pz_d += SQ_ERI_F_S_S_P_d_coefs*I_ERI_F3z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_d_coefs = delta;
      I_ERI_F3x_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_d += SQ_ERI_F_S_P_S_d_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_bc_coefs = beta*gamma;
      I_ERI_G4x_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_bc_coefs = beta*gamma;
      I_ERI_F3x_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_b_coefs = beta;
      I_ERI_G4x_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_b_coefs = beta;
      I_ERI_F3x_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_bb_coefs = beta*beta;
      I_ERI_H5x_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_bb_coefs = beta*beta;
      I_ERI_G4x_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_bb_coefs = beta*beta;
      I_ERI_F3x_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_cd_coefs = gamma*delta;
      I_ERI_F3x_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_cd += SQ_ERI_F_S_F_S_cd_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_cd_coefs = gamma*delta;
      I_ERI_F3x_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_cd += SQ_ERI_F_S_D_S_cd_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_bd_coefs = beta*delta;
      I_ERI_G4x_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_bd += SQ_ERI_G_S_D_S_bd_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_bd_coefs = beta*delta;
      I_ERI_G4x_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_bd += SQ_ERI_G_S_P_S_bd_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_bd_coefs = beta*delta;
      I_ERI_F3x_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_bd += SQ_ERI_F_S_D_S_bd_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_bd_coefs = beta*delta;
      I_ERI_F3x_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_bd += SQ_ERI_F_S_P_S_bd_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_dd_coefs = delta*delta;
      I_ERI_F3x_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_dd += SQ_ERI_F_S_F_S_dd_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_dd_coefs = delta*delta;
      I_ERI_F3x_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_dd += SQ_ERI_F_S_D_S_dd_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_dd_coefs = delta*delta;
      I_ERI_F3x_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_dd += SQ_ERI_F_S_P_S_dd_coefs*I_ERI_F3z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_F_S_P_P_bd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bd
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_bd = I_ERI_F3x_S_D2x_S_bd+CDX*I_ERI_F3x_S_Px_S_bd;
  Double I_ERI_F2xy_S_Px_Px_bd = I_ERI_F2xy_S_D2x_S_bd+CDX*I_ERI_F2xy_S_Px_S_bd;
  Double I_ERI_F2xz_S_Px_Px_bd = I_ERI_F2xz_S_D2x_S_bd+CDX*I_ERI_F2xz_S_Px_S_bd;
  Double I_ERI_Fx2y_S_Px_Px_bd = I_ERI_Fx2y_S_D2x_S_bd+CDX*I_ERI_Fx2y_S_Px_S_bd;
  Double I_ERI_Fxyz_S_Px_Px_bd = I_ERI_Fxyz_S_D2x_S_bd+CDX*I_ERI_Fxyz_S_Px_S_bd;
  Double I_ERI_Fx2z_S_Px_Px_bd = I_ERI_Fx2z_S_D2x_S_bd+CDX*I_ERI_Fx2z_S_Px_S_bd;
  Double I_ERI_F3y_S_Px_Px_bd = I_ERI_F3y_S_D2x_S_bd+CDX*I_ERI_F3y_S_Px_S_bd;
  Double I_ERI_F2yz_S_Px_Px_bd = I_ERI_F2yz_S_D2x_S_bd+CDX*I_ERI_F2yz_S_Px_S_bd;
  Double I_ERI_Fy2z_S_Px_Px_bd = I_ERI_Fy2z_S_D2x_S_bd+CDX*I_ERI_Fy2z_S_Px_S_bd;
  Double I_ERI_F3z_S_Px_Px_bd = I_ERI_F3z_S_D2x_S_bd+CDX*I_ERI_F3z_S_Px_S_bd;
  Double I_ERI_F3x_S_Py_Px_bd = I_ERI_F3x_S_Dxy_S_bd+CDX*I_ERI_F3x_S_Py_S_bd;
  Double I_ERI_F2xy_S_Py_Px_bd = I_ERI_F2xy_S_Dxy_S_bd+CDX*I_ERI_F2xy_S_Py_S_bd;
  Double I_ERI_F2xz_S_Py_Px_bd = I_ERI_F2xz_S_Dxy_S_bd+CDX*I_ERI_F2xz_S_Py_S_bd;
  Double I_ERI_Fx2y_S_Py_Px_bd = I_ERI_Fx2y_S_Dxy_S_bd+CDX*I_ERI_Fx2y_S_Py_S_bd;
  Double I_ERI_Fxyz_S_Py_Px_bd = I_ERI_Fxyz_S_Dxy_S_bd+CDX*I_ERI_Fxyz_S_Py_S_bd;
  Double I_ERI_Fx2z_S_Py_Px_bd = I_ERI_Fx2z_S_Dxy_S_bd+CDX*I_ERI_Fx2z_S_Py_S_bd;
  Double I_ERI_F3y_S_Py_Px_bd = I_ERI_F3y_S_Dxy_S_bd+CDX*I_ERI_F3y_S_Py_S_bd;
  Double I_ERI_F2yz_S_Py_Px_bd = I_ERI_F2yz_S_Dxy_S_bd+CDX*I_ERI_F2yz_S_Py_S_bd;
  Double I_ERI_Fy2z_S_Py_Px_bd = I_ERI_Fy2z_S_Dxy_S_bd+CDX*I_ERI_Fy2z_S_Py_S_bd;
  Double I_ERI_F3z_S_Py_Px_bd = I_ERI_F3z_S_Dxy_S_bd+CDX*I_ERI_F3z_S_Py_S_bd;
  Double I_ERI_F3x_S_Pz_Px_bd = I_ERI_F3x_S_Dxz_S_bd+CDX*I_ERI_F3x_S_Pz_S_bd;
  Double I_ERI_F2xy_S_Pz_Px_bd = I_ERI_F2xy_S_Dxz_S_bd+CDX*I_ERI_F2xy_S_Pz_S_bd;
  Double I_ERI_F2xz_S_Pz_Px_bd = I_ERI_F2xz_S_Dxz_S_bd+CDX*I_ERI_F2xz_S_Pz_S_bd;
  Double I_ERI_Fx2y_S_Pz_Px_bd = I_ERI_Fx2y_S_Dxz_S_bd+CDX*I_ERI_Fx2y_S_Pz_S_bd;
  Double I_ERI_Fxyz_S_Pz_Px_bd = I_ERI_Fxyz_S_Dxz_S_bd+CDX*I_ERI_Fxyz_S_Pz_S_bd;
  Double I_ERI_Fx2z_S_Pz_Px_bd = I_ERI_Fx2z_S_Dxz_S_bd+CDX*I_ERI_Fx2z_S_Pz_S_bd;
  Double I_ERI_F3y_S_Pz_Px_bd = I_ERI_F3y_S_Dxz_S_bd+CDX*I_ERI_F3y_S_Pz_S_bd;
  Double I_ERI_F2yz_S_Pz_Px_bd = I_ERI_F2yz_S_Dxz_S_bd+CDX*I_ERI_F2yz_S_Pz_S_bd;
  Double I_ERI_Fy2z_S_Pz_Px_bd = I_ERI_Fy2z_S_Dxz_S_bd+CDX*I_ERI_Fy2z_S_Pz_S_bd;
  Double I_ERI_F3z_S_Pz_Px_bd = I_ERI_F3z_S_Dxz_S_bd+CDX*I_ERI_F3z_S_Pz_S_bd;
  Double I_ERI_F3x_S_Px_Py_bd = I_ERI_F3x_S_Dxy_S_bd+CDY*I_ERI_F3x_S_Px_S_bd;
  Double I_ERI_F2xy_S_Px_Py_bd = I_ERI_F2xy_S_Dxy_S_bd+CDY*I_ERI_F2xy_S_Px_S_bd;
  Double I_ERI_F2xz_S_Px_Py_bd = I_ERI_F2xz_S_Dxy_S_bd+CDY*I_ERI_F2xz_S_Px_S_bd;
  Double I_ERI_Fx2y_S_Px_Py_bd = I_ERI_Fx2y_S_Dxy_S_bd+CDY*I_ERI_Fx2y_S_Px_S_bd;
  Double I_ERI_Fxyz_S_Px_Py_bd = I_ERI_Fxyz_S_Dxy_S_bd+CDY*I_ERI_Fxyz_S_Px_S_bd;
  Double I_ERI_Fx2z_S_Px_Py_bd = I_ERI_Fx2z_S_Dxy_S_bd+CDY*I_ERI_Fx2z_S_Px_S_bd;
  Double I_ERI_F3y_S_Px_Py_bd = I_ERI_F3y_S_Dxy_S_bd+CDY*I_ERI_F3y_S_Px_S_bd;
  Double I_ERI_F2yz_S_Px_Py_bd = I_ERI_F2yz_S_Dxy_S_bd+CDY*I_ERI_F2yz_S_Px_S_bd;
  Double I_ERI_Fy2z_S_Px_Py_bd = I_ERI_Fy2z_S_Dxy_S_bd+CDY*I_ERI_Fy2z_S_Px_S_bd;
  Double I_ERI_F3z_S_Px_Py_bd = I_ERI_F3z_S_Dxy_S_bd+CDY*I_ERI_F3z_S_Px_S_bd;
  Double I_ERI_F3x_S_Py_Py_bd = I_ERI_F3x_S_D2y_S_bd+CDY*I_ERI_F3x_S_Py_S_bd;
  Double I_ERI_F2xy_S_Py_Py_bd = I_ERI_F2xy_S_D2y_S_bd+CDY*I_ERI_F2xy_S_Py_S_bd;
  Double I_ERI_F2xz_S_Py_Py_bd = I_ERI_F2xz_S_D2y_S_bd+CDY*I_ERI_F2xz_S_Py_S_bd;
  Double I_ERI_Fx2y_S_Py_Py_bd = I_ERI_Fx2y_S_D2y_S_bd+CDY*I_ERI_Fx2y_S_Py_S_bd;
  Double I_ERI_Fxyz_S_Py_Py_bd = I_ERI_Fxyz_S_D2y_S_bd+CDY*I_ERI_Fxyz_S_Py_S_bd;
  Double I_ERI_Fx2z_S_Py_Py_bd = I_ERI_Fx2z_S_D2y_S_bd+CDY*I_ERI_Fx2z_S_Py_S_bd;
  Double I_ERI_F3y_S_Py_Py_bd = I_ERI_F3y_S_D2y_S_bd+CDY*I_ERI_F3y_S_Py_S_bd;
  Double I_ERI_F2yz_S_Py_Py_bd = I_ERI_F2yz_S_D2y_S_bd+CDY*I_ERI_F2yz_S_Py_S_bd;
  Double I_ERI_Fy2z_S_Py_Py_bd = I_ERI_Fy2z_S_D2y_S_bd+CDY*I_ERI_Fy2z_S_Py_S_bd;
  Double I_ERI_F3z_S_Py_Py_bd = I_ERI_F3z_S_D2y_S_bd+CDY*I_ERI_F3z_S_Py_S_bd;
  Double I_ERI_F3x_S_Pz_Py_bd = I_ERI_F3x_S_Dyz_S_bd+CDY*I_ERI_F3x_S_Pz_S_bd;
  Double I_ERI_F2xy_S_Pz_Py_bd = I_ERI_F2xy_S_Dyz_S_bd+CDY*I_ERI_F2xy_S_Pz_S_bd;
  Double I_ERI_F2xz_S_Pz_Py_bd = I_ERI_F2xz_S_Dyz_S_bd+CDY*I_ERI_F2xz_S_Pz_S_bd;
  Double I_ERI_Fx2y_S_Pz_Py_bd = I_ERI_Fx2y_S_Dyz_S_bd+CDY*I_ERI_Fx2y_S_Pz_S_bd;
  Double I_ERI_Fxyz_S_Pz_Py_bd = I_ERI_Fxyz_S_Dyz_S_bd+CDY*I_ERI_Fxyz_S_Pz_S_bd;
  Double I_ERI_Fx2z_S_Pz_Py_bd = I_ERI_Fx2z_S_Dyz_S_bd+CDY*I_ERI_Fx2z_S_Pz_S_bd;
  Double I_ERI_F3y_S_Pz_Py_bd = I_ERI_F3y_S_Dyz_S_bd+CDY*I_ERI_F3y_S_Pz_S_bd;
  Double I_ERI_F2yz_S_Pz_Py_bd = I_ERI_F2yz_S_Dyz_S_bd+CDY*I_ERI_F2yz_S_Pz_S_bd;
  Double I_ERI_Fy2z_S_Pz_Py_bd = I_ERI_Fy2z_S_Dyz_S_bd+CDY*I_ERI_Fy2z_S_Pz_S_bd;
  Double I_ERI_F3z_S_Pz_Py_bd = I_ERI_F3z_S_Dyz_S_bd+CDY*I_ERI_F3z_S_Pz_S_bd;
  Double I_ERI_F3x_S_Px_Pz_bd = I_ERI_F3x_S_Dxz_S_bd+CDZ*I_ERI_F3x_S_Px_S_bd;
  Double I_ERI_F2xy_S_Px_Pz_bd = I_ERI_F2xy_S_Dxz_S_bd+CDZ*I_ERI_F2xy_S_Px_S_bd;
  Double I_ERI_F2xz_S_Px_Pz_bd = I_ERI_F2xz_S_Dxz_S_bd+CDZ*I_ERI_F2xz_S_Px_S_bd;
  Double I_ERI_Fx2y_S_Px_Pz_bd = I_ERI_Fx2y_S_Dxz_S_bd+CDZ*I_ERI_Fx2y_S_Px_S_bd;
  Double I_ERI_Fxyz_S_Px_Pz_bd = I_ERI_Fxyz_S_Dxz_S_bd+CDZ*I_ERI_Fxyz_S_Px_S_bd;
  Double I_ERI_Fx2z_S_Px_Pz_bd = I_ERI_Fx2z_S_Dxz_S_bd+CDZ*I_ERI_Fx2z_S_Px_S_bd;
  Double I_ERI_F3y_S_Px_Pz_bd = I_ERI_F3y_S_Dxz_S_bd+CDZ*I_ERI_F3y_S_Px_S_bd;
  Double I_ERI_F2yz_S_Px_Pz_bd = I_ERI_F2yz_S_Dxz_S_bd+CDZ*I_ERI_F2yz_S_Px_S_bd;
  Double I_ERI_Fy2z_S_Px_Pz_bd = I_ERI_Fy2z_S_Dxz_S_bd+CDZ*I_ERI_Fy2z_S_Px_S_bd;
  Double I_ERI_F3z_S_Px_Pz_bd = I_ERI_F3z_S_Dxz_S_bd+CDZ*I_ERI_F3z_S_Px_S_bd;
  Double I_ERI_F3x_S_Py_Pz_bd = I_ERI_F3x_S_Dyz_S_bd+CDZ*I_ERI_F3x_S_Py_S_bd;
  Double I_ERI_F2xy_S_Py_Pz_bd = I_ERI_F2xy_S_Dyz_S_bd+CDZ*I_ERI_F2xy_S_Py_S_bd;
  Double I_ERI_F2xz_S_Py_Pz_bd = I_ERI_F2xz_S_Dyz_S_bd+CDZ*I_ERI_F2xz_S_Py_S_bd;
  Double I_ERI_Fx2y_S_Py_Pz_bd = I_ERI_Fx2y_S_Dyz_S_bd+CDZ*I_ERI_Fx2y_S_Py_S_bd;
  Double I_ERI_Fxyz_S_Py_Pz_bd = I_ERI_Fxyz_S_Dyz_S_bd+CDZ*I_ERI_Fxyz_S_Py_S_bd;
  Double I_ERI_Fx2z_S_Py_Pz_bd = I_ERI_Fx2z_S_Dyz_S_bd+CDZ*I_ERI_Fx2z_S_Py_S_bd;
  Double I_ERI_F3y_S_Py_Pz_bd = I_ERI_F3y_S_Dyz_S_bd+CDZ*I_ERI_F3y_S_Py_S_bd;
  Double I_ERI_F2yz_S_Py_Pz_bd = I_ERI_F2yz_S_Dyz_S_bd+CDZ*I_ERI_F2yz_S_Py_S_bd;
  Double I_ERI_Fy2z_S_Py_Pz_bd = I_ERI_Fy2z_S_Dyz_S_bd+CDZ*I_ERI_Fy2z_S_Py_S_bd;
  Double I_ERI_F3z_S_Py_Pz_bd = I_ERI_F3z_S_Dyz_S_bd+CDZ*I_ERI_F3z_S_Py_S_bd;
  Double I_ERI_F3x_S_Pz_Pz_bd = I_ERI_F3x_S_D2z_S_bd+CDZ*I_ERI_F3x_S_Pz_S_bd;
  Double I_ERI_F2xy_S_Pz_Pz_bd = I_ERI_F2xy_S_D2z_S_bd+CDZ*I_ERI_F2xy_S_Pz_S_bd;
  Double I_ERI_F2xz_S_Pz_Pz_bd = I_ERI_F2xz_S_D2z_S_bd+CDZ*I_ERI_F2xz_S_Pz_S_bd;
  Double I_ERI_Fx2y_S_Pz_Pz_bd = I_ERI_Fx2y_S_D2z_S_bd+CDZ*I_ERI_Fx2y_S_Pz_S_bd;
  Double I_ERI_Fxyz_S_Pz_Pz_bd = I_ERI_Fxyz_S_D2z_S_bd+CDZ*I_ERI_Fxyz_S_Pz_S_bd;
  Double I_ERI_Fx2z_S_Pz_Pz_bd = I_ERI_Fx2z_S_D2z_S_bd+CDZ*I_ERI_Fx2z_S_Pz_S_bd;
  Double I_ERI_F3y_S_Pz_Pz_bd = I_ERI_F3y_S_D2z_S_bd+CDZ*I_ERI_F3y_S_Pz_S_bd;
  Double I_ERI_F2yz_S_Pz_Pz_bd = I_ERI_F2yz_S_D2z_S_bd+CDZ*I_ERI_F2yz_S_Pz_S_bd;
  Double I_ERI_Fy2z_S_Pz_Pz_bd = I_ERI_Fy2z_S_D2z_S_bd+CDZ*I_ERI_Fy2z_S_Pz_S_bd;
  Double I_ERI_F3z_S_Pz_Pz_bd = I_ERI_F3z_S_D2z_S_bd+CDZ*I_ERI_F3z_S_Pz_S_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_P_bd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bd
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bd
   ************************************************************/
  Double I_ERI_G4x_S_Px_Px_bd = I_ERI_G4x_S_D2x_S_bd+CDX*I_ERI_G4x_S_Px_S_bd;
  Double I_ERI_G3xy_S_Px_Px_bd = I_ERI_G3xy_S_D2x_S_bd+CDX*I_ERI_G3xy_S_Px_S_bd;
  Double I_ERI_G3xz_S_Px_Px_bd = I_ERI_G3xz_S_D2x_S_bd+CDX*I_ERI_G3xz_S_Px_S_bd;
  Double I_ERI_G2x2y_S_Px_Px_bd = I_ERI_G2x2y_S_D2x_S_bd+CDX*I_ERI_G2x2y_S_Px_S_bd;
  Double I_ERI_G2xyz_S_Px_Px_bd = I_ERI_G2xyz_S_D2x_S_bd+CDX*I_ERI_G2xyz_S_Px_S_bd;
  Double I_ERI_G2x2z_S_Px_Px_bd = I_ERI_G2x2z_S_D2x_S_bd+CDX*I_ERI_G2x2z_S_Px_S_bd;
  Double I_ERI_Gx3y_S_Px_Px_bd = I_ERI_Gx3y_S_D2x_S_bd+CDX*I_ERI_Gx3y_S_Px_S_bd;
  Double I_ERI_Gx2yz_S_Px_Px_bd = I_ERI_Gx2yz_S_D2x_S_bd+CDX*I_ERI_Gx2yz_S_Px_S_bd;
  Double I_ERI_Gxy2z_S_Px_Px_bd = I_ERI_Gxy2z_S_D2x_S_bd+CDX*I_ERI_Gxy2z_S_Px_S_bd;
  Double I_ERI_Gx3z_S_Px_Px_bd = I_ERI_Gx3z_S_D2x_S_bd+CDX*I_ERI_Gx3z_S_Px_S_bd;
  Double I_ERI_G4y_S_Px_Px_bd = I_ERI_G4y_S_D2x_S_bd+CDX*I_ERI_G4y_S_Px_S_bd;
  Double I_ERI_G3yz_S_Px_Px_bd = I_ERI_G3yz_S_D2x_S_bd+CDX*I_ERI_G3yz_S_Px_S_bd;
  Double I_ERI_G2y2z_S_Px_Px_bd = I_ERI_G2y2z_S_D2x_S_bd+CDX*I_ERI_G2y2z_S_Px_S_bd;
  Double I_ERI_Gy3z_S_Px_Px_bd = I_ERI_Gy3z_S_D2x_S_bd+CDX*I_ERI_Gy3z_S_Px_S_bd;
  Double I_ERI_G4z_S_Px_Px_bd = I_ERI_G4z_S_D2x_S_bd+CDX*I_ERI_G4z_S_Px_S_bd;
  Double I_ERI_G4x_S_Py_Px_bd = I_ERI_G4x_S_Dxy_S_bd+CDX*I_ERI_G4x_S_Py_S_bd;
  Double I_ERI_G3xy_S_Py_Px_bd = I_ERI_G3xy_S_Dxy_S_bd+CDX*I_ERI_G3xy_S_Py_S_bd;
  Double I_ERI_G3xz_S_Py_Px_bd = I_ERI_G3xz_S_Dxy_S_bd+CDX*I_ERI_G3xz_S_Py_S_bd;
  Double I_ERI_G2x2y_S_Py_Px_bd = I_ERI_G2x2y_S_Dxy_S_bd+CDX*I_ERI_G2x2y_S_Py_S_bd;
  Double I_ERI_G2xyz_S_Py_Px_bd = I_ERI_G2xyz_S_Dxy_S_bd+CDX*I_ERI_G2xyz_S_Py_S_bd;
  Double I_ERI_G2x2z_S_Py_Px_bd = I_ERI_G2x2z_S_Dxy_S_bd+CDX*I_ERI_G2x2z_S_Py_S_bd;
  Double I_ERI_Gx3y_S_Py_Px_bd = I_ERI_Gx3y_S_Dxy_S_bd+CDX*I_ERI_Gx3y_S_Py_S_bd;
  Double I_ERI_Gx2yz_S_Py_Px_bd = I_ERI_Gx2yz_S_Dxy_S_bd+CDX*I_ERI_Gx2yz_S_Py_S_bd;
  Double I_ERI_Gxy2z_S_Py_Px_bd = I_ERI_Gxy2z_S_Dxy_S_bd+CDX*I_ERI_Gxy2z_S_Py_S_bd;
  Double I_ERI_Gx3z_S_Py_Px_bd = I_ERI_Gx3z_S_Dxy_S_bd+CDX*I_ERI_Gx3z_S_Py_S_bd;
  Double I_ERI_G4y_S_Py_Px_bd = I_ERI_G4y_S_Dxy_S_bd+CDX*I_ERI_G4y_S_Py_S_bd;
  Double I_ERI_G3yz_S_Py_Px_bd = I_ERI_G3yz_S_Dxy_S_bd+CDX*I_ERI_G3yz_S_Py_S_bd;
  Double I_ERI_G2y2z_S_Py_Px_bd = I_ERI_G2y2z_S_Dxy_S_bd+CDX*I_ERI_G2y2z_S_Py_S_bd;
  Double I_ERI_Gy3z_S_Py_Px_bd = I_ERI_Gy3z_S_Dxy_S_bd+CDX*I_ERI_Gy3z_S_Py_S_bd;
  Double I_ERI_G4z_S_Py_Px_bd = I_ERI_G4z_S_Dxy_S_bd+CDX*I_ERI_G4z_S_Py_S_bd;
  Double I_ERI_G4x_S_Pz_Px_bd = I_ERI_G4x_S_Dxz_S_bd+CDX*I_ERI_G4x_S_Pz_S_bd;
  Double I_ERI_G3xy_S_Pz_Px_bd = I_ERI_G3xy_S_Dxz_S_bd+CDX*I_ERI_G3xy_S_Pz_S_bd;
  Double I_ERI_G3xz_S_Pz_Px_bd = I_ERI_G3xz_S_Dxz_S_bd+CDX*I_ERI_G3xz_S_Pz_S_bd;
  Double I_ERI_G2x2y_S_Pz_Px_bd = I_ERI_G2x2y_S_Dxz_S_bd+CDX*I_ERI_G2x2y_S_Pz_S_bd;
  Double I_ERI_G2xyz_S_Pz_Px_bd = I_ERI_G2xyz_S_Dxz_S_bd+CDX*I_ERI_G2xyz_S_Pz_S_bd;
  Double I_ERI_G2x2z_S_Pz_Px_bd = I_ERI_G2x2z_S_Dxz_S_bd+CDX*I_ERI_G2x2z_S_Pz_S_bd;
  Double I_ERI_Gx3y_S_Pz_Px_bd = I_ERI_Gx3y_S_Dxz_S_bd+CDX*I_ERI_Gx3y_S_Pz_S_bd;
  Double I_ERI_Gx2yz_S_Pz_Px_bd = I_ERI_Gx2yz_S_Dxz_S_bd+CDX*I_ERI_Gx2yz_S_Pz_S_bd;
  Double I_ERI_Gxy2z_S_Pz_Px_bd = I_ERI_Gxy2z_S_Dxz_S_bd+CDX*I_ERI_Gxy2z_S_Pz_S_bd;
  Double I_ERI_Gx3z_S_Pz_Px_bd = I_ERI_Gx3z_S_Dxz_S_bd+CDX*I_ERI_Gx3z_S_Pz_S_bd;
  Double I_ERI_G4y_S_Pz_Px_bd = I_ERI_G4y_S_Dxz_S_bd+CDX*I_ERI_G4y_S_Pz_S_bd;
  Double I_ERI_G3yz_S_Pz_Px_bd = I_ERI_G3yz_S_Dxz_S_bd+CDX*I_ERI_G3yz_S_Pz_S_bd;
  Double I_ERI_G2y2z_S_Pz_Px_bd = I_ERI_G2y2z_S_Dxz_S_bd+CDX*I_ERI_G2y2z_S_Pz_S_bd;
  Double I_ERI_Gy3z_S_Pz_Px_bd = I_ERI_Gy3z_S_Dxz_S_bd+CDX*I_ERI_Gy3z_S_Pz_S_bd;
  Double I_ERI_G4z_S_Pz_Px_bd = I_ERI_G4z_S_Dxz_S_bd+CDX*I_ERI_G4z_S_Pz_S_bd;
  Double I_ERI_G4x_S_Px_Py_bd = I_ERI_G4x_S_Dxy_S_bd+CDY*I_ERI_G4x_S_Px_S_bd;
  Double I_ERI_G3xy_S_Px_Py_bd = I_ERI_G3xy_S_Dxy_S_bd+CDY*I_ERI_G3xy_S_Px_S_bd;
  Double I_ERI_G3xz_S_Px_Py_bd = I_ERI_G3xz_S_Dxy_S_bd+CDY*I_ERI_G3xz_S_Px_S_bd;
  Double I_ERI_G2x2y_S_Px_Py_bd = I_ERI_G2x2y_S_Dxy_S_bd+CDY*I_ERI_G2x2y_S_Px_S_bd;
  Double I_ERI_G2xyz_S_Px_Py_bd = I_ERI_G2xyz_S_Dxy_S_bd+CDY*I_ERI_G2xyz_S_Px_S_bd;
  Double I_ERI_G2x2z_S_Px_Py_bd = I_ERI_G2x2z_S_Dxy_S_bd+CDY*I_ERI_G2x2z_S_Px_S_bd;
  Double I_ERI_Gx3y_S_Px_Py_bd = I_ERI_Gx3y_S_Dxy_S_bd+CDY*I_ERI_Gx3y_S_Px_S_bd;
  Double I_ERI_Gx2yz_S_Px_Py_bd = I_ERI_Gx2yz_S_Dxy_S_bd+CDY*I_ERI_Gx2yz_S_Px_S_bd;
  Double I_ERI_Gxy2z_S_Px_Py_bd = I_ERI_Gxy2z_S_Dxy_S_bd+CDY*I_ERI_Gxy2z_S_Px_S_bd;
  Double I_ERI_Gx3z_S_Px_Py_bd = I_ERI_Gx3z_S_Dxy_S_bd+CDY*I_ERI_Gx3z_S_Px_S_bd;
  Double I_ERI_G4y_S_Px_Py_bd = I_ERI_G4y_S_Dxy_S_bd+CDY*I_ERI_G4y_S_Px_S_bd;
  Double I_ERI_G3yz_S_Px_Py_bd = I_ERI_G3yz_S_Dxy_S_bd+CDY*I_ERI_G3yz_S_Px_S_bd;
  Double I_ERI_G2y2z_S_Px_Py_bd = I_ERI_G2y2z_S_Dxy_S_bd+CDY*I_ERI_G2y2z_S_Px_S_bd;
  Double I_ERI_Gy3z_S_Px_Py_bd = I_ERI_Gy3z_S_Dxy_S_bd+CDY*I_ERI_Gy3z_S_Px_S_bd;
  Double I_ERI_G4z_S_Px_Py_bd = I_ERI_G4z_S_Dxy_S_bd+CDY*I_ERI_G4z_S_Px_S_bd;
  Double I_ERI_G4x_S_Py_Py_bd = I_ERI_G4x_S_D2y_S_bd+CDY*I_ERI_G4x_S_Py_S_bd;
  Double I_ERI_G3xy_S_Py_Py_bd = I_ERI_G3xy_S_D2y_S_bd+CDY*I_ERI_G3xy_S_Py_S_bd;
  Double I_ERI_G3xz_S_Py_Py_bd = I_ERI_G3xz_S_D2y_S_bd+CDY*I_ERI_G3xz_S_Py_S_bd;
  Double I_ERI_G2x2y_S_Py_Py_bd = I_ERI_G2x2y_S_D2y_S_bd+CDY*I_ERI_G2x2y_S_Py_S_bd;
  Double I_ERI_G2xyz_S_Py_Py_bd = I_ERI_G2xyz_S_D2y_S_bd+CDY*I_ERI_G2xyz_S_Py_S_bd;
  Double I_ERI_G2x2z_S_Py_Py_bd = I_ERI_G2x2z_S_D2y_S_bd+CDY*I_ERI_G2x2z_S_Py_S_bd;
  Double I_ERI_Gx3y_S_Py_Py_bd = I_ERI_Gx3y_S_D2y_S_bd+CDY*I_ERI_Gx3y_S_Py_S_bd;
  Double I_ERI_Gx2yz_S_Py_Py_bd = I_ERI_Gx2yz_S_D2y_S_bd+CDY*I_ERI_Gx2yz_S_Py_S_bd;
  Double I_ERI_Gxy2z_S_Py_Py_bd = I_ERI_Gxy2z_S_D2y_S_bd+CDY*I_ERI_Gxy2z_S_Py_S_bd;
  Double I_ERI_Gx3z_S_Py_Py_bd = I_ERI_Gx3z_S_D2y_S_bd+CDY*I_ERI_Gx3z_S_Py_S_bd;
  Double I_ERI_G4y_S_Py_Py_bd = I_ERI_G4y_S_D2y_S_bd+CDY*I_ERI_G4y_S_Py_S_bd;
  Double I_ERI_G3yz_S_Py_Py_bd = I_ERI_G3yz_S_D2y_S_bd+CDY*I_ERI_G3yz_S_Py_S_bd;
  Double I_ERI_G2y2z_S_Py_Py_bd = I_ERI_G2y2z_S_D2y_S_bd+CDY*I_ERI_G2y2z_S_Py_S_bd;
  Double I_ERI_Gy3z_S_Py_Py_bd = I_ERI_Gy3z_S_D2y_S_bd+CDY*I_ERI_Gy3z_S_Py_S_bd;
  Double I_ERI_G4z_S_Py_Py_bd = I_ERI_G4z_S_D2y_S_bd+CDY*I_ERI_G4z_S_Py_S_bd;
  Double I_ERI_G4x_S_Pz_Py_bd = I_ERI_G4x_S_Dyz_S_bd+CDY*I_ERI_G4x_S_Pz_S_bd;
  Double I_ERI_G3xy_S_Pz_Py_bd = I_ERI_G3xy_S_Dyz_S_bd+CDY*I_ERI_G3xy_S_Pz_S_bd;
  Double I_ERI_G3xz_S_Pz_Py_bd = I_ERI_G3xz_S_Dyz_S_bd+CDY*I_ERI_G3xz_S_Pz_S_bd;
  Double I_ERI_G2x2y_S_Pz_Py_bd = I_ERI_G2x2y_S_Dyz_S_bd+CDY*I_ERI_G2x2y_S_Pz_S_bd;
  Double I_ERI_G2xyz_S_Pz_Py_bd = I_ERI_G2xyz_S_Dyz_S_bd+CDY*I_ERI_G2xyz_S_Pz_S_bd;
  Double I_ERI_G2x2z_S_Pz_Py_bd = I_ERI_G2x2z_S_Dyz_S_bd+CDY*I_ERI_G2x2z_S_Pz_S_bd;
  Double I_ERI_Gx3y_S_Pz_Py_bd = I_ERI_Gx3y_S_Dyz_S_bd+CDY*I_ERI_Gx3y_S_Pz_S_bd;
  Double I_ERI_Gx2yz_S_Pz_Py_bd = I_ERI_Gx2yz_S_Dyz_S_bd+CDY*I_ERI_Gx2yz_S_Pz_S_bd;
  Double I_ERI_Gxy2z_S_Pz_Py_bd = I_ERI_Gxy2z_S_Dyz_S_bd+CDY*I_ERI_Gxy2z_S_Pz_S_bd;
  Double I_ERI_Gx3z_S_Pz_Py_bd = I_ERI_Gx3z_S_Dyz_S_bd+CDY*I_ERI_Gx3z_S_Pz_S_bd;
  Double I_ERI_G4y_S_Pz_Py_bd = I_ERI_G4y_S_Dyz_S_bd+CDY*I_ERI_G4y_S_Pz_S_bd;
  Double I_ERI_G3yz_S_Pz_Py_bd = I_ERI_G3yz_S_Dyz_S_bd+CDY*I_ERI_G3yz_S_Pz_S_bd;
  Double I_ERI_G2y2z_S_Pz_Py_bd = I_ERI_G2y2z_S_Dyz_S_bd+CDY*I_ERI_G2y2z_S_Pz_S_bd;
  Double I_ERI_Gy3z_S_Pz_Py_bd = I_ERI_Gy3z_S_Dyz_S_bd+CDY*I_ERI_Gy3z_S_Pz_S_bd;
  Double I_ERI_G4z_S_Pz_Py_bd = I_ERI_G4z_S_Dyz_S_bd+CDY*I_ERI_G4z_S_Pz_S_bd;
  Double I_ERI_G4x_S_Px_Pz_bd = I_ERI_G4x_S_Dxz_S_bd+CDZ*I_ERI_G4x_S_Px_S_bd;
  Double I_ERI_G3xy_S_Px_Pz_bd = I_ERI_G3xy_S_Dxz_S_bd+CDZ*I_ERI_G3xy_S_Px_S_bd;
  Double I_ERI_G3xz_S_Px_Pz_bd = I_ERI_G3xz_S_Dxz_S_bd+CDZ*I_ERI_G3xz_S_Px_S_bd;
  Double I_ERI_G2x2y_S_Px_Pz_bd = I_ERI_G2x2y_S_Dxz_S_bd+CDZ*I_ERI_G2x2y_S_Px_S_bd;
  Double I_ERI_G2xyz_S_Px_Pz_bd = I_ERI_G2xyz_S_Dxz_S_bd+CDZ*I_ERI_G2xyz_S_Px_S_bd;
  Double I_ERI_G2x2z_S_Px_Pz_bd = I_ERI_G2x2z_S_Dxz_S_bd+CDZ*I_ERI_G2x2z_S_Px_S_bd;
  Double I_ERI_Gx3y_S_Px_Pz_bd = I_ERI_Gx3y_S_Dxz_S_bd+CDZ*I_ERI_Gx3y_S_Px_S_bd;
  Double I_ERI_Gx2yz_S_Px_Pz_bd = I_ERI_Gx2yz_S_Dxz_S_bd+CDZ*I_ERI_Gx2yz_S_Px_S_bd;
  Double I_ERI_Gxy2z_S_Px_Pz_bd = I_ERI_Gxy2z_S_Dxz_S_bd+CDZ*I_ERI_Gxy2z_S_Px_S_bd;
  Double I_ERI_Gx3z_S_Px_Pz_bd = I_ERI_Gx3z_S_Dxz_S_bd+CDZ*I_ERI_Gx3z_S_Px_S_bd;
  Double I_ERI_G4y_S_Px_Pz_bd = I_ERI_G4y_S_Dxz_S_bd+CDZ*I_ERI_G4y_S_Px_S_bd;
  Double I_ERI_G3yz_S_Px_Pz_bd = I_ERI_G3yz_S_Dxz_S_bd+CDZ*I_ERI_G3yz_S_Px_S_bd;
  Double I_ERI_G2y2z_S_Px_Pz_bd = I_ERI_G2y2z_S_Dxz_S_bd+CDZ*I_ERI_G2y2z_S_Px_S_bd;
  Double I_ERI_Gy3z_S_Px_Pz_bd = I_ERI_Gy3z_S_Dxz_S_bd+CDZ*I_ERI_Gy3z_S_Px_S_bd;
  Double I_ERI_G4z_S_Px_Pz_bd = I_ERI_G4z_S_Dxz_S_bd+CDZ*I_ERI_G4z_S_Px_S_bd;
  Double I_ERI_G4x_S_Py_Pz_bd = I_ERI_G4x_S_Dyz_S_bd+CDZ*I_ERI_G4x_S_Py_S_bd;
  Double I_ERI_G3xy_S_Py_Pz_bd = I_ERI_G3xy_S_Dyz_S_bd+CDZ*I_ERI_G3xy_S_Py_S_bd;
  Double I_ERI_G3xz_S_Py_Pz_bd = I_ERI_G3xz_S_Dyz_S_bd+CDZ*I_ERI_G3xz_S_Py_S_bd;
  Double I_ERI_G2x2y_S_Py_Pz_bd = I_ERI_G2x2y_S_Dyz_S_bd+CDZ*I_ERI_G2x2y_S_Py_S_bd;
  Double I_ERI_G2xyz_S_Py_Pz_bd = I_ERI_G2xyz_S_Dyz_S_bd+CDZ*I_ERI_G2xyz_S_Py_S_bd;
  Double I_ERI_G2x2z_S_Py_Pz_bd = I_ERI_G2x2z_S_Dyz_S_bd+CDZ*I_ERI_G2x2z_S_Py_S_bd;
  Double I_ERI_Gx3y_S_Py_Pz_bd = I_ERI_Gx3y_S_Dyz_S_bd+CDZ*I_ERI_Gx3y_S_Py_S_bd;
  Double I_ERI_Gx2yz_S_Py_Pz_bd = I_ERI_Gx2yz_S_Dyz_S_bd+CDZ*I_ERI_Gx2yz_S_Py_S_bd;
  Double I_ERI_Gxy2z_S_Py_Pz_bd = I_ERI_Gxy2z_S_Dyz_S_bd+CDZ*I_ERI_Gxy2z_S_Py_S_bd;
  Double I_ERI_Gx3z_S_Py_Pz_bd = I_ERI_Gx3z_S_Dyz_S_bd+CDZ*I_ERI_Gx3z_S_Py_S_bd;
  Double I_ERI_G4y_S_Py_Pz_bd = I_ERI_G4y_S_Dyz_S_bd+CDZ*I_ERI_G4y_S_Py_S_bd;
  Double I_ERI_G3yz_S_Py_Pz_bd = I_ERI_G3yz_S_Dyz_S_bd+CDZ*I_ERI_G3yz_S_Py_S_bd;
  Double I_ERI_G2y2z_S_Py_Pz_bd = I_ERI_G2y2z_S_Dyz_S_bd+CDZ*I_ERI_G2y2z_S_Py_S_bd;
  Double I_ERI_Gy3z_S_Py_Pz_bd = I_ERI_Gy3z_S_Dyz_S_bd+CDZ*I_ERI_Gy3z_S_Py_S_bd;
  Double I_ERI_G4z_S_Py_Pz_bd = I_ERI_G4z_S_Dyz_S_bd+CDZ*I_ERI_G4z_S_Py_S_bd;
  Double I_ERI_G4x_S_Pz_Pz_bd = I_ERI_G4x_S_D2z_S_bd+CDZ*I_ERI_G4x_S_Pz_S_bd;
  Double I_ERI_G3xy_S_Pz_Pz_bd = I_ERI_G3xy_S_D2z_S_bd+CDZ*I_ERI_G3xy_S_Pz_S_bd;
  Double I_ERI_G3xz_S_Pz_Pz_bd = I_ERI_G3xz_S_D2z_S_bd+CDZ*I_ERI_G3xz_S_Pz_S_bd;
  Double I_ERI_G2x2y_S_Pz_Pz_bd = I_ERI_G2x2y_S_D2z_S_bd+CDZ*I_ERI_G2x2y_S_Pz_S_bd;
  Double I_ERI_G2xyz_S_Pz_Pz_bd = I_ERI_G2xyz_S_D2z_S_bd+CDZ*I_ERI_G2xyz_S_Pz_S_bd;
  Double I_ERI_G2x2z_S_Pz_Pz_bd = I_ERI_G2x2z_S_D2z_S_bd+CDZ*I_ERI_G2x2z_S_Pz_S_bd;
  Double I_ERI_Gx3y_S_Pz_Pz_bd = I_ERI_Gx3y_S_D2z_S_bd+CDZ*I_ERI_Gx3y_S_Pz_S_bd;
  Double I_ERI_Gx2yz_S_Pz_Pz_bd = I_ERI_Gx2yz_S_D2z_S_bd+CDZ*I_ERI_Gx2yz_S_Pz_S_bd;
  Double I_ERI_Gxy2z_S_Pz_Pz_bd = I_ERI_Gxy2z_S_D2z_S_bd+CDZ*I_ERI_Gxy2z_S_Pz_S_bd;
  Double I_ERI_Gx3z_S_Pz_Pz_bd = I_ERI_Gx3z_S_D2z_S_bd+CDZ*I_ERI_Gx3z_S_Pz_S_bd;
  Double I_ERI_G4y_S_Pz_Pz_bd = I_ERI_G4y_S_D2z_S_bd+CDZ*I_ERI_G4y_S_Pz_S_bd;
  Double I_ERI_G3yz_S_Pz_Pz_bd = I_ERI_G3yz_S_D2z_S_bd+CDZ*I_ERI_G3yz_S_Pz_S_bd;
  Double I_ERI_G2y2z_S_Pz_Pz_bd = I_ERI_G2y2z_S_D2z_S_bd+CDZ*I_ERI_G2y2z_S_Pz_S_bd;
  Double I_ERI_Gy3z_S_Pz_Pz_bd = I_ERI_Gy3z_S_D2z_S_bd+CDZ*I_ERI_Gy3z_S_Pz_S_bd;
  Double I_ERI_G4z_S_Pz_Pz_bd = I_ERI_G4z_S_D2z_S_bd+CDZ*I_ERI_G4z_S_Pz_S_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_cd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_cd
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_cd = I_ERI_F3x_S_F3x_S_cd+CDX*I_ERI_F3x_S_D2x_S_cd;
  Double I_ERI_F2xy_S_D2x_Px_cd = I_ERI_F2xy_S_F3x_S_cd+CDX*I_ERI_F2xy_S_D2x_S_cd;
  Double I_ERI_F2xz_S_D2x_Px_cd = I_ERI_F2xz_S_F3x_S_cd+CDX*I_ERI_F2xz_S_D2x_S_cd;
  Double I_ERI_Fx2y_S_D2x_Px_cd = I_ERI_Fx2y_S_F3x_S_cd+CDX*I_ERI_Fx2y_S_D2x_S_cd;
  Double I_ERI_Fxyz_S_D2x_Px_cd = I_ERI_Fxyz_S_F3x_S_cd+CDX*I_ERI_Fxyz_S_D2x_S_cd;
  Double I_ERI_Fx2z_S_D2x_Px_cd = I_ERI_Fx2z_S_F3x_S_cd+CDX*I_ERI_Fx2z_S_D2x_S_cd;
  Double I_ERI_F3y_S_D2x_Px_cd = I_ERI_F3y_S_F3x_S_cd+CDX*I_ERI_F3y_S_D2x_S_cd;
  Double I_ERI_F2yz_S_D2x_Px_cd = I_ERI_F2yz_S_F3x_S_cd+CDX*I_ERI_F2yz_S_D2x_S_cd;
  Double I_ERI_Fy2z_S_D2x_Px_cd = I_ERI_Fy2z_S_F3x_S_cd+CDX*I_ERI_Fy2z_S_D2x_S_cd;
  Double I_ERI_F3z_S_D2x_Px_cd = I_ERI_F3z_S_F3x_S_cd+CDX*I_ERI_F3z_S_D2x_S_cd;
  Double I_ERI_F3x_S_Dxy_Px_cd = I_ERI_F3x_S_F2xy_S_cd+CDX*I_ERI_F3x_S_Dxy_S_cd;
  Double I_ERI_F2xy_S_Dxy_Px_cd = I_ERI_F2xy_S_F2xy_S_cd+CDX*I_ERI_F2xy_S_Dxy_S_cd;
  Double I_ERI_F2xz_S_Dxy_Px_cd = I_ERI_F2xz_S_F2xy_S_cd+CDX*I_ERI_F2xz_S_Dxy_S_cd;
  Double I_ERI_Fx2y_S_Dxy_Px_cd = I_ERI_Fx2y_S_F2xy_S_cd+CDX*I_ERI_Fx2y_S_Dxy_S_cd;
  Double I_ERI_Fxyz_S_Dxy_Px_cd = I_ERI_Fxyz_S_F2xy_S_cd+CDX*I_ERI_Fxyz_S_Dxy_S_cd;
  Double I_ERI_Fx2z_S_Dxy_Px_cd = I_ERI_Fx2z_S_F2xy_S_cd+CDX*I_ERI_Fx2z_S_Dxy_S_cd;
  Double I_ERI_F3y_S_Dxy_Px_cd = I_ERI_F3y_S_F2xy_S_cd+CDX*I_ERI_F3y_S_Dxy_S_cd;
  Double I_ERI_F2yz_S_Dxy_Px_cd = I_ERI_F2yz_S_F2xy_S_cd+CDX*I_ERI_F2yz_S_Dxy_S_cd;
  Double I_ERI_Fy2z_S_Dxy_Px_cd = I_ERI_Fy2z_S_F2xy_S_cd+CDX*I_ERI_Fy2z_S_Dxy_S_cd;
  Double I_ERI_F3z_S_Dxy_Px_cd = I_ERI_F3z_S_F2xy_S_cd+CDX*I_ERI_F3z_S_Dxy_S_cd;
  Double I_ERI_F3x_S_Dxz_Px_cd = I_ERI_F3x_S_F2xz_S_cd+CDX*I_ERI_F3x_S_Dxz_S_cd;
  Double I_ERI_F2xy_S_Dxz_Px_cd = I_ERI_F2xy_S_F2xz_S_cd+CDX*I_ERI_F2xy_S_Dxz_S_cd;
  Double I_ERI_F2xz_S_Dxz_Px_cd = I_ERI_F2xz_S_F2xz_S_cd+CDX*I_ERI_F2xz_S_Dxz_S_cd;
  Double I_ERI_Fx2y_S_Dxz_Px_cd = I_ERI_Fx2y_S_F2xz_S_cd+CDX*I_ERI_Fx2y_S_Dxz_S_cd;
  Double I_ERI_Fxyz_S_Dxz_Px_cd = I_ERI_Fxyz_S_F2xz_S_cd+CDX*I_ERI_Fxyz_S_Dxz_S_cd;
  Double I_ERI_Fx2z_S_Dxz_Px_cd = I_ERI_Fx2z_S_F2xz_S_cd+CDX*I_ERI_Fx2z_S_Dxz_S_cd;
  Double I_ERI_F3y_S_Dxz_Px_cd = I_ERI_F3y_S_F2xz_S_cd+CDX*I_ERI_F3y_S_Dxz_S_cd;
  Double I_ERI_F2yz_S_Dxz_Px_cd = I_ERI_F2yz_S_F2xz_S_cd+CDX*I_ERI_F2yz_S_Dxz_S_cd;
  Double I_ERI_Fy2z_S_Dxz_Px_cd = I_ERI_Fy2z_S_F2xz_S_cd+CDX*I_ERI_Fy2z_S_Dxz_S_cd;
  Double I_ERI_F3z_S_Dxz_Px_cd = I_ERI_F3z_S_F2xz_S_cd+CDX*I_ERI_F3z_S_Dxz_S_cd;
  Double I_ERI_F3x_S_D2y_Px_cd = I_ERI_F3x_S_Fx2y_S_cd+CDX*I_ERI_F3x_S_D2y_S_cd;
  Double I_ERI_F2xy_S_D2y_Px_cd = I_ERI_F2xy_S_Fx2y_S_cd+CDX*I_ERI_F2xy_S_D2y_S_cd;
  Double I_ERI_F2xz_S_D2y_Px_cd = I_ERI_F2xz_S_Fx2y_S_cd+CDX*I_ERI_F2xz_S_D2y_S_cd;
  Double I_ERI_Fx2y_S_D2y_Px_cd = I_ERI_Fx2y_S_Fx2y_S_cd+CDX*I_ERI_Fx2y_S_D2y_S_cd;
  Double I_ERI_Fxyz_S_D2y_Px_cd = I_ERI_Fxyz_S_Fx2y_S_cd+CDX*I_ERI_Fxyz_S_D2y_S_cd;
  Double I_ERI_Fx2z_S_D2y_Px_cd = I_ERI_Fx2z_S_Fx2y_S_cd+CDX*I_ERI_Fx2z_S_D2y_S_cd;
  Double I_ERI_F3y_S_D2y_Px_cd = I_ERI_F3y_S_Fx2y_S_cd+CDX*I_ERI_F3y_S_D2y_S_cd;
  Double I_ERI_F2yz_S_D2y_Px_cd = I_ERI_F2yz_S_Fx2y_S_cd+CDX*I_ERI_F2yz_S_D2y_S_cd;
  Double I_ERI_Fy2z_S_D2y_Px_cd = I_ERI_Fy2z_S_Fx2y_S_cd+CDX*I_ERI_Fy2z_S_D2y_S_cd;
  Double I_ERI_F3z_S_D2y_Px_cd = I_ERI_F3z_S_Fx2y_S_cd+CDX*I_ERI_F3z_S_D2y_S_cd;
  Double I_ERI_F3x_S_Dyz_Px_cd = I_ERI_F3x_S_Fxyz_S_cd+CDX*I_ERI_F3x_S_Dyz_S_cd;
  Double I_ERI_F2xy_S_Dyz_Px_cd = I_ERI_F2xy_S_Fxyz_S_cd+CDX*I_ERI_F2xy_S_Dyz_S_cd;
  Double I_ERI_F2xz_S_Dyz_Px_cd = I_ERI_F2xz_S_Fxyz_S_cd+CDX*I_ERI_F2xz_S_Dyz_S_cd;
  Double I_ERI_Fx2y_S_Dyz_Px_cd = I_ERI_Fx2y_S_Fxyz_S_cd+CDX*I_ERI_Fx2y_S_Dyz_S_cd;
  Double I_ERI_Fxyz_S_Dyz_Px_cd = I_ERI_Fxyz_S_Fxyz_S_cd+CDX*I_ERI_Fxyz_S_Dyz_S_cd;
  Double I_ERI_Fx2z_S_Dyz_Px_cd = I_ERI_Fx2z_S_Fxyz_S_cd+CDX*I_ERI_Fx2z_S_Dyz_S_cd;
  Double I_ERI_F3y_S_Dyz_Px_cd = I_ERI_F3y_S_Fxyz_S_cd+CDX*I_ERI_F3y_S_Dyz_S_cd;
  Double I_ERI_F2yz_S_Dyz_Px_cd = I_ERI_F2yz_S_Fxyz_S_cd+CDX*I_ERI_F2yz_S_Dyz_S_cd;
  Double I_ERI_Fy2z_S_Dyz_Px_cd = I_ERI_Fy2z_S_Fxyz_S_cd+CDX*I_ERI_Fy2z_S_Dyz_S_cd;
  Double I_ERI_F3z_S_Dyz_Px_cd = I_ERI_F3z_S_Fxyz_S_cd+CDX*I_ERI_F3z_S_Dyz_S_cd;
  Double I_ERI_F3x_S_D2z_Px_cd = I_ERI_F3x_S_Fx2z_S_cd+CDX*I_ERI_F3x_S_D2z_S_cd;
  Double I_ERI_F2xy_S_D2z_Px_cd = I_ERI_F2xy_S_Fx2z_S_cd+CDX*I_ERI_F2xy_S_D2z_S_cd;
  Double I_ERI_F2xz_S_D2z_Px_cd = I_ERI_F2xz_S_Fx2z_S_cd+CDX*I_ERI_F2xz_S_D2z_S_cd;
  Double I_ERI_Fx2y_S_D2z_Px_cd = I_ERI_Fx2y_S_Fx2z_S_cd+CDX*I_ERI_Fx2y_S_D2z_S_cd;
  Double I_ERI_Fxyz_S_D2z_Px_cd = I_ERI_Fxyz_S_Fx2z_S_cd+CDX*I_ERI_Fxyz_S_D2z_S_cd;
  Double I_ERI_Fx2z_S_D2z_Px_cd = I_ERI_Fx2z_S_Fx2z_S_cd+CDX*I_ERI_Fx2z_S_D2z_S_cd;
  Double I_ERI_F3y_S_D2z_Px_cd = I_ERI_F3y_S_Fx2z_S_cd+CDX*I_ERI_F3y_S_D2z_S_cd;
  Double I_ERI_F2yz_S_D2z_Px_cd = I_ERI_F2yz_S_Fx2z_S_cd+CDX*I_ERI_F2yz_S_D2z_S_cd;
  Double I_ERI_Fy2z_S_D2z_Px_cd = I_ERI_Fy2z_S_Fx2z_S_cd+CDX*I_ERI_Fy2z_S_D2z_S_cd;
  Double I_ERI_F3z_S_D2z_Px_cd = I_ERI_F3z_S_Fx2z_S_cd+CDX*I_ERI_F3z_S_D2z_S_cd;
  Double I_ERI_F3x_S_D2x_Py_cd = I_ERI_F3x_S_F2xy_S_cd+CDY*I_ERI_F3x_S_D2x_S_cd;
  Double I_ERI_F2xy_S_D2x_Py_cd = I_ERI_F2xy_S_F2xy_S_cd+CDY*I_ERI_F2xy_S_D2x_S_cd;
  Double I_ERI_F2xz_S_D2x_Py_cd = I_ERI_F2xz_S_F2xy_S_cd+CDY*I_ERI_F2xz_S_D2x_S_cd;
  Double I_ERI_Fx2y_S_D2x_Py_cd = I_ERI_Fx2y_S_F2xy_S_cd+CDY*I_ERI_Fx2y_S_D2x_S_cd;
  Double I_ERI_Fxyz_S_D2x_Py_cd = I_ERI_Fxyz_S_F2xy_S_cd+CDY*I_ERI_Fxyz_S_D2x_S_cd;
  Double I_ERI_Fx2z_S_D2x_Py_cd = I_ERI_Fx2z_S_F2xy_S_cd+CDY*I_ERI_Fx2z_S_D2x_S_cd;
  Double I_ERI_F3y_S_D2x_Py_cd = I_ERI_F3y_S_F2xy_S_cd+CDY*I_ERI_F3y_S_D2x_S_cd;
  Double I_ERI_F2yz_S_D2x_Py_cd = I_ERI_F2yz_S_F2xy_S_cd+CDY*I_ERI_F2yz_S_D2x_S_cd;
  Double I_ERI_Fy2z_S_D2x_Py_cd = I_ERI_Fy2z_S_F2xy_S_cd+CDY*I_ERI_Fy2z_S_D2x_S_cd;
  Double I_ERI_F3z_S_D2x_Py_cd = I_ERI_F3z_S_F2xy_S_cd+CDY*I_ERI_F3z_S_D2x_S_cd;
  Double I_ERI_F3x_S_Dxy_Py_cd = I_ERI_F3x_S_Fx2y_S_cd+CDY*I_ERI_F3x_S_Dxy_S_cd;
  Double I_ERI_F2xy_S_Dxy_Py_cd = I_ERI_F2xy_S_Fx2y_S_cd+CDY*I_ERI_F2xy_S_Dxy_S_cd;
  Double I_ERI_F2xz_S_Dxy_Py_cd = I_ERI_F2xz_S_Fx2y_S_cd+CDY*I_ERI_F2xz_S_Dxy_S_cd;
  Double I_ERI_Fx2y_S_Dxy_Py_cd = I_ERI_Fx2y_S_Fx2y_S_cd+CDY*I_ERI_Fx2y_S_Dxy_S_cd;
  Double I_ERI_Fxyz_S_Dxy_Py_cd = I_ERI_Fxyz_S_Fx2y_S_cd+CDY*I_ERI_Fxyz_S_Dxy_S_cd;
  Double I_ERI_Fx2z_S_Dxy_Py_cd = I_ERI_Fx2z_S_Fx2y_S_cd+CDY*I_ERI_Fx2z_S_Dxy_S_cd;
  Double I_ERI_F3y_S_Dxy_Py_cd = I_ERI_F3y_S_Fx2y_S_cd+CDY*I_ERI_F3y_S_Dxy_S_cd;
  Double I_ERI_F2yz_S_Dxy_Py_cd = I_ERI_F2yz_S_Fx2y_S_cd+CDY*I_ERI_F2yz_S_Dxy_S_cd;
  Double I_ERI_Fy2z_S_Dxy_Py_cd = I_ERI_Fy2z_S_Fx2y_S_cd+CDY*I_ERI_Fy2z_S_Dxy_S_cd;
  Double I_ERI_F3z_S_Dxy_Py_cd = I_ERI_F3z_S_Fx2y_S_cd+CDY*I_ERI_F3z_S_Dxy_S_cd;
  Double I_ERI_F3x_S_Dxz_Py_cd = I_ERI_F3x_S_Fxyz_S_cd+CDY*I_ERI_F3x_S_Dxz_S_cd;
  Double I_ERI_F2xy_S_Dxz_Py_cd = I_ERI_F2xy_S_Fxyz_S_cd+CDY*I_ERI_F2xy_S_Dxz_S_cd;
  Double I_ERI_F2xz_S_Dxz_Py_cd = I_ERI_F2xz_S_Fxyz_S_cd+CDY*I_ERI_F2xz_S_Dxz_S_cd;
  Double I_ERI_Fx2y_S_Dxz_Py_cd = I_ERI_Fx2y_S_Fxyz_S_cd+CDY*I_ERI_Fx2y_S_Dxz_S_cd;
  Double I_ERI_Fxyz_S_Dxz_Py_cd = I_ERI_Fxyz_S_Fxyz_S_cd+CDY*I_ERI_Fxyz_S_Dxz_S_cd;
  Double I_ERI_Fx2z_S_Dxz_Py_cd = I_ERI_Fx2z_S_Fxyz_S_cd+CDY*I_ERI_Fx2z_S_Dxz_S_cd;
  Double I_ERI_F3y_S_Dxz_Py_cd = I_ERI_F3y_S_Fxyz_S_cd+CDY*I_ERI_F3y_S_Dxz_S_cd;
  Double I_ERI_F2yz_S_Dxz_Py_cd = I_ERI_F2yz_S_Fxyz_S_cd+CDY*I_ERI_F2yz_S_Dxz_S_cd;
  Double I_ERI_Fy2z_S_Dxz_Py_cd = I_ERI_Fy2z_S_Fxyz_S_cd+CDY*I_ERI_Fy2z_S_Dxz_S_cd;
  Double I_ERI_F3z_S_Dxz_Py_cd = I_ERI_F3z_S_Fxyz_S_cd+CDY*I_ERI_F3z_S_Dxz_S_cd;
  Double I_ERI_F3x_S_D2y_Py_cd = I_ERI_F3x_S_F3y_S_cd+CDY*I_ERI_F3x_S_D2y_S_cd;
  Double I_ERI_F2xy_S_D2y_Py_cd = I_ERI_F2xy_S_F3y_S_cd+CDY*I_ERI_F2xy_S_D2y_S_cd;
  Double I_ERI_F2xz_S_D2y_Py_cd = I_ERI_F2xz_S_F3y_S_cd+CDY*I_ERI_F2xz_S_D2y_S_cd;
  Double I_ERI_Fx2y_S_D2y_Py_cd = I_ERI_Fx2y_S_F3y_S_cd+CDY*I_ERI_Fx2y_S_D2y_S_cd;
  Double I_ERI_Fxyz_S_D2y_Py_cd = I_ERI_Fxyz_S_F3y_S_cd+CDY*I_ERI_Fxyz_S_D2y_S_cd;
  Double I_ERI_Fx2z_S_D2y_Py_cd = I_ERI_Fx2z_S_F3y_S_cd+CDY*I_ERI_Fx2z_S_D2y_S_cd;
  Double I_ERI_F3y_S_D2y_Py_cd = I_ERI_F3y_S_F3y_S_cd+CDY*I_ERI_F3y_S_D2y_S_cd;
  Double I_ERI_F2yz_S_D2y_Py_cd = I_ERI_F2yz_S_F3y_S_cd+CDY*I_ERI_F2yz_S_D2y_S_cd;
  Double I_ERI_Fy2z_S_D2y_Py_cd = I_ERI_Fy2z_S_F3y_S_cd+CDY*I_ERI_Fy2z_S_D2y_S_cd;
  Double I_ERI_F3z_S_D2y_Py_cd = I_ERI_F3z_S_F3y_S_cd+CDY*I_ERI_F3z_S_D2y_S_cd;
  Double I_ERI_F3x_S_Dyz_Py_cd = I_ERI_F3x_S_F2yz_S_cd+CDY*I_ERI_F3x_S_Dyz_S_cd;
  Double I_ERI_F2xy_S_Dyz_Py_cd = I_ERI_F2xy_S_F2yz_S_cd+CDY*I_ERI_F2xy_S_Dyz_S_cd;
  Double I_ERI_F2xz_S_Dyz_Py_cd = I_ERI_F2xz_S_F2yz_S_cd+CDY*I_ERI_F2xz_S_Dyz_S_cd;
  Double I_ERI_Fx2y_S_Dyz_Py_cd = I_ERI_Fx2y_S_F2yz_S_cd+CDY*I_ERI_Fx2y_S_Dyz_S_cd;
  Double I_ERI_Fxyz_S_Dyz_Py_cd = I_ERI_Fxyz_S_F2yz_S_cd+CDY*I_ERI_Fxyz_S_Dyz_S_cd;
  Double I_ERI_Fx2z_S_Dyz_Py_cd = I_ERI_Fx2z_S_F2yz_S_cd+CDY*I_ERI_Fx2z_S_Dyz_S_cd;
  Double I_ERI_F3y_S_Dyz_Py_cd = I_ERI_F3y_S_F2yz_S_cd+CDY*I_ERI_F3y_S_Dyz_S_cd;
  Double I_ERI_F2yz_S_Dyz_Py_cd = I_ERI_F2yz_S_F2yz_S_cd+CDY*I_ERI_F2yz_S_Dyz_S_cd;
  Double I_ERI_Fy2z_S_Dyz_Py_cd = I_ERI_Fy2z_S_F2yz_S_cd+CDY*I_ERI_Fy2z_S_Dyz_S_cd;
  Double I_ERI_F3z_S_Dyz_Py_cd = I_ERI_F3z_S_F2yz_S_cd+CDY*I_ERI_F3z_S_Dyz_S_cd;
  Double I_ERI_F3x_S_D2z_Py_cd = I_ERI_F3x_S_Fy2z_S_cd+CDY*I_ERI_F3x_S_D2z_S_cd;
  Double I_ERI_F2xy_S_D2z_Py_cd = I_ERI_F2xy_S_Fy2z_S_cd+CDY*I_ERI_F2xy_S_D2z_S_cd;
  Double I_ERI_F2xz_S_D2z_Py_cd = I_ERI_F2xz_S_Fy2z_S_cd+CDY*I_ERI_F2xz_S_D2z_S_cd;
  Double I_ERI_Fx2y_S_D2z_Py_cd = I_ERI_Fx2y_S_Fy2z_S_cd+CDY*I_ERI_Fx2y_S_D2z_S_cd;
  Double I_ERI_Fxyz_S_D2z_Py_cd = I_ERI_Fxyz_S_Fy2z_S_cd+CDY*I_ERI_Fxyz_S_D2z_S_cd;
  Double I_ERI_Fx2z_S_D2z_Py_cd = I_ERI_Fx2z_S_Fy2z_S_cd+CDY*I_ERI_Fx2z_S_D2z_S_cd;
  Double I_ERI_F3y_S_D2z_Py_cd = I_ERI_F3y_S_Fy2z_S_cd+CDY*I_ERI_F3y_S_D2z_S_cd;
  Double I_ERI_F2yz_S_D2z_Py_cd = I_ERI_F2yz_S_Fy2z_S_cd+CDY*I_ERI_F2yz_S_D2z_S_cd;
  Double I_ERI_Fy2z_S_D2z_Py_cd = I_ERI_Fy2z_S_Fy2z_S_cd+CDY*I_ERI_Fy2z_S_D2z_S_cd;
  Double I_ERI_F3z_S_D2z_Py_cd = I_ERI_F3z_S_Fy2z_S_cd+CDY*I_ERI_F3z_S_D2z_S_cd;
  Double I_ERI_F3x_S_D2x_Pz_cd = I_ERI_F3x_S_F2xz_S_cd+CDZ*I_ERI_F3x_S_D2x_S_cd;
  Double I_ERI_F2xy_S_D2x_Pz_cd = I_ERI_F2xy_S_F2xz_S_cd+CDZ*I_ERI_F2xy_S_D2x_S_cd;
  Double I_ERI_F2xz_S_D2x_Pz_cd = I_ERI_F2xz_S_F2xz_S_cd+CDZ*I_ERI_F2xz_S_D2x_S_cd;
  Double I_ERI_Fx2y_S_D2x_Pz_cd = I_ERI_Fx2y_S_F2xz_S_cd+CDZ*I_ERI_Fx2y_S_D2x_S_cd;
  Double I_ERI_Fxyz_S_D2x_Pz_cd = I_ERI_Fxyz_S_F2xz_S_cd+CDZ*I_ERI_Fxyz_S_D2x_S_cd;
  Double I_ERI_Fx2z_S_D2x_Pz_cd = I_ERI_Fx2z_S_F2xz_S_cd+CDZ*I_ERI_Fx2z_S_D2x_S_cd;
  Double I_ERI_F3y_S_D2x_Pz_cd = I_ERI_F3y_S_F2xz_S_cd+CDZ*I_ERI_F3y_S_D2x_S_cd;
  Double I_ERI_F2yz_S_D2x_Pz_cd = I_ERI_F2yz_S_F2xz_S_cd+CDZ*I_ERI_F2yz_S_D2x_S_cd;
  Double I_ERI_Fy2z_S_D2x_Pz_cd = I_ERI_Fy2z_S_F2xz_S_cd+CDZ*I_ERI_Fy2z_S_D2x_S_cd;
  Double I_ERI_F3z_S_D2x_Pz_cd = I_ERI_F3z_S_F2xz_S_cd+CDZ*I_ERI_F3z_S_D2x_S_cd;
  Double I_ERI_F3x_S_Dxy_Pz_cd = I_ERI_F3x_S_Fxyz_S_cd+CDZ*I_ERI_F3x_S_Dxy_S_cd;
  Double I_ERI_F2xy_S_Dxy_Pz_cd = I_ERI_F2xy_S_Fxyz_S_cd+CDZ*I_ERI_F2xy_S_Dxy_S_cd;
  Double I_ERI_F2xz_S_Dxy_Pz_cd = I_ERI_F2xz_S_Fxyz_S_cd+CDZ*I_ERI_F2xz_S_Dxy_S_cd;
  Double I_ERI_Fx2y_S_Dxy_Pz_cd = I_ERI_Fx2y_S_Fxyz_S_cd+CDZ*I_ERI_Fx2y_S_Dxy_S_cd;
  Double I_ERI_Fxyz_S_Dxy_Pz_cd = I_ERI_Fxyz_S_Fxyz_S_cd+CDZ*I_ERI_Fxyz_S_Dxy_S_cd;
  Double I_ERI_Fx2z_S_Dxy_Pz_cd = I_ERI_Fx2z_S_Fxyz_S_cd+CDZ*I_ERI_Fx2z_S_Dxy_S_cd;
  Double I_ERI_F3y_S_Dxy_Pz_cd = I_ERI_F3y_S_Fxyz_S_cd+CDZ*I_ERI_F3y_S_Dxy_S_cd;
  Double I_ERI_F2yz_S_Dxy_Pz_cd = I_ERI_F2yz_S_Fxyz_S_cd+CDZ*I_ERI_F2yz_S_Dxy_S_cd;
  Double I_ERI_Fy2z_S_Dxy_Pz_cd = I_ERI_Fy2z_S_Fxyz_S_cd+CDZ*I_ERI_Fy2z_S_Dxy_S_cd;
  Double I_ERI_F3z_S_Dxy_Pz_cd = I_ERI_F3z_S_Fxyz_S_cd+CDZ*I_ERI_F3z_S_Dxy_S_cd;
  Double I_ERI_F3x_S_Dxz_Pz_cd = I_ERI_F3x_S_Fx2z_S_cd+CDZ*I_ERI_F3x_S_Dxz_S_cd;
  Double I_ERI_F2xy_S_Dxz_Pz_cd = I_ERI_F2xy_S_Fx2z_S_cd+CDZ*I_ERI_F2xy_S_Dxz_S_cd;
  Double I_ERI_F2xz_S_Dxz_Pz_cd = I_ERI_F2xz_S_Fx2z_S_cd+CDZ*I_ERI_F2xz_S_Dxz_S_cd;
  Double I_ERI_Fx2y_S_Dxz_Pz_cd = I_ERI_Fx2y_S_Fx2z_S_cd+CDZ*I_ERI_Fx2y_S_Dxz_S_cd;
  Double I_ERI_Fxyz_S_Dxz_Pz_cd = I_ERI_Fxyz_S_Fx2z_S_cd+CDZ*I_ERI_Fxyz_S_Dxz_S_cd;
  Double I_ERI_Fx2z_S_Dxz_Pz_cd = I_ERI_Fx2z_S_Fx2z_S_cd+CDZ*I_ERI_Fx2z_S_Dxz_S_cd;
  Double I_ERI_F3y_S_Dxz_Pz_cd = I_ERI_F3y_S_Fx2z_S_cd+CDZ*I_ERI_F3y_S_Dxz_S_cd;
  Double I_ERI_F2yz_S_Dxz_Pz_cd = I_ERI_F2yz_S_Fx2z_S_cd+CDZ*I_ERI_F2yz_S_Dxz_S_cd;
  Double I_ERI_Fy2z_S_Dxz_Pz_cd = I_ERI_Fy2z_S_Fx2z_S_cd+CDZ*I_ERI_Fy2z_S_Dxz_S_cd;
  Double I_ERI_F3z_S_Dxz_Pz_cd = I_ERI_F3z_S_Fx2z_S_cd+CDZ*I_ERI_F3z_S_Dxz_S_cd;
  Double I_ERI_F3x_S_D2y_Pz_cd = I_ERI_F3x_S_F2yz_S_cd+CDZ*I_ERI_F3x_S_D2y_S_cd;
  Double I_ERI_F2xy_S_D2y_Pz_cd = I_ERI_F2xy_S_F2yz_S_cd+CDZ*I_ERI_F2xy_S_D2y_S_cd;
  Double I_ERI_F2xz_S_D2y_Pz_cd = I_ERI_F2xz_S_F2yz_S_cd+CDZ*I_ERI_F2xz_S_D2y_S_cd;
  Double I_ERI_Fx2y_S_D2y_Pz_cd = I_ERI_Fx2y_S_F2yz_S_cd+CDZ*I_ERI_Fx2y_S_D2y_S_cd;
  Double I_ERI_Fxyz_S_D2y_Pz_cd = I_ERI_Fxyz_S_F2yz_S_cd+CDZ*I_ERI_Fxyz_S_D2y_S_cd;
  Double I_ERI_Fx2z_S_D2y_Pz_cd = I_ERI_Fx2z_S_F2yz_S_cd+CDZ*I_ERI_Fx2z_S_D2y_S_cd;
  Double I_ERI_F3y_S_D2y_Pz_cd = I_ERI_F3y_S_F2yz_S_cd+CDZ*I_ERI_F3y_S_D2y_S_cd;
  Double I_ERI_F2yz_S_D2y_Pz_cd = I_ERI_F2yz_S_F2yz_S_cd+CDZ*I_ERI_F2yz_S_D2y_S_cd;
  Double I_ERI_Fy2z_S_D2y_Pz_cd = I_ERI_Fy2z_S_F2yz_S_cd+CDZ*I_ERI_Fy2z_S_D2y_S_cd;
  Double I_ERI_F3z_S_D2y_Pz_cd = I_ERI_F3z_S_F2yz_S_cd+CDZ*I_ERI_F3z_S_D2y_S_cd;
  Double I_ERI_F3x_S_Dyz_Pz_cd = I_ERI_F3x_S_Fy2z_S_cd+CDZ*I_ERI_F3x_S_Dyz_S_cd;
  Double I_ERI_F2xy_S_Dyz_Pz_cd = I_ERI_F2xy_S_Fy2z_S_cd+CDZ*I_ERI_F2xy_S_Dyz_S_cd;
  Double I_ERI_F2xz_S_Dyz_Pz_cd = I_ERI_F2xz_S_Fy2z_S_cd+CDZ*I_ERI_F2xz_S_Dyz_S_cd;
  Double I_ERI_Fx2y_S_Dyz_Pz_cd = I_ERI_Fx2y_S_Fy2z_S_cd+CDZ*I_ERI_Fx2y_S_Dyz_S_cd;
  Double I_ERI_Fxyz_S_Dyz_Pz_cd = I_ERI_Fxyz_S_Fy2z_S_cd+CDZ*I_ERI_Fxyz_S_Dyz_S_cd;
  Double I_ERI_Fx2z_S_Dyz_Pz_cd = I_ERI_Fx2z_S_Fy2z_S_cd+CDZ*I_ERI_Fx2z_S_Dyz_S_cd;
  Double I_ERI_F3y_S_Dyz_Pz_cd = I_ERI_F3y_S_Fy2z_S_cd+CDZ*I_ERI_F3y_S_Dyz_S_cd;
  Double I_ERI_F2yz_S_Dyz_Pz_cd = I_ERI_F2yz_S_Fy2z_S_cd+CDZ*I_ERI_F2yz_S_Dyz_S_cd;
  Double I_ERI_Fy2z_S_Dyz_Pz_cd = I_ERI_Fy2z_S_Fy2z_S_cd+CDZ*I_ERI_Fy2z_S_Dyz_S_cd;
  Double I_ERI_F3z_S_Dyz_Pz_cd = I_ERI_F3z_S_Fy2z_S_cd+CDZ*I_ERI_F3z_S_Dyz_S_cd;
  Double I_ERI_F3x_S_D2z_Pz_cd = I_ERI_F3x_S_F3z_S_cd+CDZ*I_ERI_F3x_S_D2z_S_cd;
  Double I_ERI_F2xy_S_D2z_Pz_cd = I_ERI_F2xy_S_F3z_S_cd+CDZ*I_ERI_F2xy_S_D2z_S_cd;
  Double I_ERI_F2xz_S_D2z_Pz_cd = I_ERI_F2xz_S_F3z_S_cd+CDZ*I_ERI_F2xz_S_D2z_S_cd;
  Double I_ERI_Fx2y_S_D2z_Pz_cd = I_ERI_Fx2y_S_F3z_S_cd+CDZ*I_ERI_Fx2y_S_D2z_S_cd;
  Double I_ERI_Fxyz_S_D2z_Pz_cd = I_ERI_Fxyz_S_F3z_S_cd+CDZ*I_ERI_Fxyz_S_D2z_S_cd;
  Double I_ERI_Fx2z_S_D2z_Pz_cd = I_ERI_Fx2z_S_F3z_S_cd+CDZ*I_ERI_Fx2z_S_D2z_S_cd;
  Double I_ERI_F3y_S_D2z_Pz_cd = I_ERI_F3y_S_F3z_S_cd+CDZ*I_ERI_F3y_S_D2z_S_cd;
  Double I_ERI_F2yz_S_D2z_Pz_cd = I_ERI_F2yz_S_F3z_S_cd+CDZ*I_ERI_F2yz_S_D2z_S_cd;
  Double I_ERI_Fy2z_S_D2z_Pz_cd = I_ERI_Fy2z_S_F3z_S_cd+CDZ*I_ERI_Fy2z_S_D2z_S_cd;
  Double I_ERI_F3z_S_D2z_Pz_cd = I_ERI_F3z_S_F3z_S_cd+CDZ*I_ERI_F3z_S_D2z_S_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_dd
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_dd = I_ERI_F3x_S_D2x_S_dd+CDX*I_ERI_F3x_S_Px_S_dd;
  Double I_ERI_F2xy_S_Px_Px_dd = I_ERI_F2xy_S_D2x_S_dd+CDX*I_ERI_F2xy_S_Px_S_dd;
  Double I_ERI_F2xz_S_Px_Px_dd = I_ERI_F2xz_S_D2x_S_dd+CDX*I_ERI_F2xz_S_Px_S_dd;
  Double I_ERI_Fx2y_S_Px_Px_dd = I_ERI_Fx2y_S_D2x_S_dd+CDX*I_ERI_Fx2y_S_Px_S_dd;
  Double I_ERI_Fxyz_S_Px_Px_dd = I_ERI_Fxyz_S_D2x_S_dd+CDX*I_ERI_Fxyz_S_Px_S_dd;
  Double I_ERI_Fx2z_S_Px_Px_dd = I_ERI_Fx2z_S_D2x_S_dd+CDX*I_ERI_Fx2z_S_Px_S_dd;
  Double I_ERI_F3y_S_Px_Px_dd = I_ERI_F3y_S_D2x_S_dd+CDX*I_ERI_F3y_S_Px_S_dd;
  Double I_ERI_F2yz_S_Px_Px_dd = I_ERI_F2yz_S_D2x_S_dd+CDX*I_ERI_F2yz_S_Px_S_dd;
  Double I_ERI_Fy2z_S_Px_Px_dd = I_ERI_Fy2z_S_D2x_S_dd+CDX*I_ERI_Fy2z_S_Px_S_dd;
  Double I_ERI_F3z_S_Px_Px_dd = I_ERI_F3z_S_D2x_S_dd+CDX*I_ERI_F3z_S_Px_S_dd;
  Double I_ERI_F3x_S_Py_Px_dd = I_ERI_F3x_S_Dxy_S_dd+CDX*I_ERI_F3x_S_Py_S_dd;
  Double I_ERI_F2xy_S_Py_Px_dd = I_ERI_F2xy_S_Dxy_S_dd+CDX*I_ERI_F2xy_S_Py_S_dd;
  Double I_ERI_F2xz_S_Py_Px_dd = I_ERI_F2xz_S_Dxy_S_dd+CDX*I_ERI_F2xz_S_Py_S_dd;
  Double I_ERI_Fx2y_S_Py_Px_dd = I_ERI_Fx2y_S_Dxy_S_dd+CDX*I_ERI_Fx2y_S_Py_S_dd;
  Double I_ERI_Fxyz_S_Py_Px_dd = I_ERI_Fxyz_S_Dxy_S_dd+CDX*I_ERI_Fxyz_S_Py_S_dd;
  Double I_ERI_Fx2z_S_Py_Px_dd = I_ERI_Fx2z_S_Dxy_S_dd+CDX*I_ERI_Fx2z_S_Py_S_dd;
  Double I_ERI_F3y_S_Py_Px_dd = I_ERI_F3y_S_Dxy_S_dd+CDX*I_ERI_F3y_S_Py_S_dd;
  Double I_ERI_F2yz_S_Py_Px_dd = I_ERI_F2yz_S_Dxy_S_dd+CDX*I_ERI_F2yz_S_Py_S_dd;
  Double I_ERI_Fy2z_S_Py_Px_dd = I_ERI_Fy2z_S_Dxy_S_dd+CDX*I_ERI_Fy2z_S_Py_S_dd;
  Double I_ERI_F3z_S_Py_Px_dd = I_ERI_F3z_S_Dxy_S_dd+CDX*I_ERI_F3z_S_Py_S_dd;
  Double I_ERI_F3x_S_Pz_Px_dd = I_ERI_F3x_S_Dxz_S_dd+CDX*I_ERI_F3x_S_Pz_S_dd;
  Double I_ERI_F2xy_S_Pz_Px_dd = I_ERI_F2xy_S_Dxz_S_dd+CDX*I_ERI_F2xy_S_Pz_S_dd;
  Double I_ERI_F2xz_S_Pz_Px_dd = I_ERI_F2xz_S_Dxz_S_dd+CDX*I_ERI_F2xz_S_Pz_S_dd;
  Double I_ERI_Fx2y_S_Pz_Px_dd = I_ERI_Fx2y_S_Dxz_S_dd+CDX*I_ERI_Fx2y_S_Pz_S_dd;
  Double I_ERI_Fxyz_S_Pz_Px_dd = I_ERI_Fxyz_S_Dxz_S_dd+CDX*I_ERI_Fxyz_S_Pz_S_dd;
  Double I_ERI_Fx2z_S_Pz_Px_dd = I_ERI_Fx2z_S_Dxz_S_dd+CDX*I_ERI_Fx2z_S_Pz_S_dd;
  Double I_ERI_F3y_S_Pz_Px_dd = I_ERI_F3y_S_Dxz_S_dd+CDX*I_ERI_F3y_S_Pz_S_dd;
  Double I_ERI_F2yz_S_Pz_Px_dd = I_ERI_F2yz_S_Dxz_S_dd+CDX*I_ERI_F2yz_S_Pz_S_dd;
  Double I_ERI_Fy2z_S_Pz_Px_dd = I_ERI_Fy2z_S_Dxz_S_dd+CDX*I_ERI_Fy2z_S_Pz_S_dd;
  Double I_ERI_F3z_S_Pz_Px_dd = I_ERI_F3z_S_Dxz_S_dd+CDX*I_ERI_F3z_S_Pz_S_dd;
  Double I_ERI_F3x_S_Px_Py_dd = I_ERI_F3x_S_Dxy_S_dd+CDY*I_ERI_F3x_S_Px_S_dd;
  Double I_ERI_F2xy_S_Px_Py_dd = I_ERI_F2xy_S_Dxy_S_dd+CDY*I_ERI_F2xy_S_Px_S_dd;
  Double I_ERI_F2xz_S_Px_Py_dd = I_ERI_F2xz_S_Dxy_S_dd+CDY*I_ERI_F2xz_S_Px_S_dd;
  Double I_ERI_Fx2y_S_Px_Py_dd = I_ERI_Fx2y_S_Dxy_S_dd+CDY*I_ERI_Fx2y_S_Px_S_dd;
  Double I_ERI_Fxyz_S_Px_Py_dd = I_ERI_Fxyz_S_Dxy_S_dd+CDY*I_ERI_Fxyz_S_Px_S_dd;
  Double I_ERI_Fx2z_S_Px_Py_dd = I_ERI_Fx2z_S_Dxy_S_dd+CDY*I_ERI_Fx2z_S_Px_S_dd;
  Double I_ERI_F3y_S_Px_Py_dd = I_ERI_F3y_S_Dxy_S_dd+CDY*I_ERI_F3y_S_Px_S_dd;
  Double I_ERI_F2yz_S_Px_Py_dd = I_ERI_F2yz_S_Dxy_S_dd+CDY*I_ERI_F2yz_S_Px_S_dd;
  Double I_ERI_Fy2z_S_Px_Py_dd = I_ERI_Fy2z_S_Dxy_S_dd+CDY*I_ERI_Fy2z_S_Px_S_dd;
  Double I_ERI_F3z_S_Px_Py_dd = I_ERI_F3z_S_Dxy_S_dd+CDY*I_ERI_F3z_S_Px_S_dd;
  Double I_ERI_F3x_S_Py_Py_dd = I_ERI_F3x_S_D2y_S_dd+CDY*I_ERI_F3x_S_Py_S_dd;
  Double I_ERI_F2xy_S_Py_Py_dd = I_ERI_F2xy_S_D2y_S_dd+CDY*I_ERI_F2xy_S_Py_S_dd;
  Double I_ERI_F2xz_S_Py_Py_dd = I_ERI_F2xz_S_D2y_S_dd+CDY*I_ERI_F2xz_S_Py_S_dd;
  Double I_ERI_Fx2y_S_Py_Py_dd = I_ERI_Fx2y_S_D2y_S_dd+CDY*I_ERI_Fx2y_S_Py_S_dd;
  Double I_ERI_Fxyz_S_Py_Py_dd = I_ERI_Fxyz_S_D2y_S_dd+CDY*I_ERI_Fxyz_S_Py_S_dd;
  Double I_ERI_Fx2z_S_Py_Py_dd = I_ERI_Fx2z_S_D2y_S_dd+CDY*I_ERI_Fx2z_S_Py_S_dd;
  Double I_ERI_F3y_S_Py_Py_dd = I_ERI_F3y_S_D2y_S_dd+CDY*I_ERI_F3y_S_Py_S_dd;
  Double I_ERI_F2yz_S_Py_Py_dd = I_ERI_F2yz_S_D2y_S_dd+CDY*I_ERI_F2yz_S_Py_S_dd;
  Double I_ERI_Fy2z_S_Py_Py_dd = I_ERI_Fy2z_S_D2y_S_dd+CDY*I_ERI_Fy2z_S_Py_S_dd;
  Double I_ERI_F3z_S_Py_Py_dd = I_ERI_F3z_S_D2y_S_dd+CDY*I_ERI_F3z_S_Py_S_dd;
  Double I_ERI_F3x_S_Pz_Py_dd = I_ERI_F3x_S_Dyz_S_dd+CDY*I_ERI_F3x_S_Pz_S_dd;
  Double I_ERI_F2xy_S_Pz_Py_dd = I_ERI_F2xy_S_Dyz_S_dd+CDY*I_ERI_F2xy_S_Pz_S_dd;
  Double I_ERI_F2xz_S_Pz_Py_dd = I_ERI_F2xz_S_Dyz_S_dd+CDY*I_ERI_F2xz_S_Pz_S_dd;
  Double I_ERI_Fx2y_S_Pz_Py_dd = I_ERI_Fx2y_S_Dyz_S_dd+CDY*I_ERI_Fx2y_S_Pz_S_dd;
  Double I_ERI_Fxyz_S_Pz_Py_dd = I_ERI_Fxyz_S_Dyz_S_dd+CDY*I_ERI_Fxyz_S_Pz_S_dd;
  Double I_ERI_Fx2z_S_Pz_Py_dd = I_ERI_Fx2z_S_Dyz_S_dd+CDY*I_ERI_Fx2z_S_Pz_S_dd;
  Double I_ERI_F3y_S_Pz_Py_dd = I_ERI_F3y_S_Dyz_S_dd+CDY*I_ERI_F3y_S_Pz_S_dd;
  Double I_ERI_F2yz_S_Pz_Py_dd = I_ERI_F2yz_S_Dyz_S_dd+CDY*I_ERI_F2yz_S_Pz_S_dd;
  Double I_ERI_Fy2z_S_Pz_Py_dd = I_ERI_Fy2z_S_Dyz_S_dd+CDY*I_ERI_Fy2z_S_Pz_S_dd;
  Double I_ERI_F3z_S_Pz_Py_dd = I_ERI_F3z_S_Dyz_S_dd+CDY*I_ERI_F3z_S_Pz_S_dd;
  Double I_ERI_F3x_S_Px_Pz_dd = I_ERI_F3x_S_Dxz_S_dd+CDZ*I_ERI_F3x_S_Px_S_dd;
  Double I_ERI_F2xy_S_Px_Pz_dd = I_ERI_F2xy_S_Dxz_S_dd+CDZ*I_ERI_F2xy_S_Px_S_dd;
  Double I_ERI_F2xz_S_Px_Pz_dd = I_ERI_F2xz_S_Dxz_S_dd+CDZ*I_ERI_F2xz_S_Px_S_dd;
  Double I_ERI_Fx2y_S_Px_Pz_dd = I_ERI_Fx2y_S_Dxz_S_dd+CDZ*I_ERI_Fx2y_S_Px_S_dd;
  Double I_ERI_Fxyz_S_Px_Pz_dd = I_ERI_Fxyz_S_Dxz_S_dd+CDZ*I_ERI_Fxyz_S_Px_S_dd;
  Double I_ERI_Fx2z_S_Px_Pz_dd = I_ERI_Fx2z_S_Dxz_S_dd+CDZ*I_ERI_Fx2z_S_Px_S_dd;
  Double I_ERI_F3y_S_Px_Pz_dd = I_ERI_F3y_S_Dxz_S_dd+CDZ*I_ERI_F3y_S_Px_S_dd;
  Double I_ERI_F2yz_S_Px_Pz_dd = I_ERI_F2yz_S_Dxz_S_dd+CDZ*I_ERI_F2yz_S_Px_S_dd;
  Double I_ERI_Fy2z_S_Px_Pz_dd = I_ERI_Fy2z_S_Dxz_S_dd+CDZ*I_ERI_Fy2z_S_Px_S_dd;
  Double I_ERI_F3z_S_Px_Pz_dd = I_ERI_F3z_S_Dxz_S_dd+CDZ*I_ERI_F3z_S_Px_S_dd;
  Double I_ERI_F3x_S_Py_Pz_dd = I_ERI_F3x_S_Dyz_S_dd+CDZ*I_ERI_F3x_S_Py_S_dd;
  Double I_ERI_F2xy_S_Py_Pz_dd = I_ERI_F2xy_S_Dyz_S_dd+CDZ*I_ERI_F2xy_S_Py_S_dd;
  Double I_ERI_F2xz_S_Py_Pz_dd = I_ERI_F2xz_S_Dyz_S_dd+CDZ*I_ERI_F2xz_S_Py_S_dd;
  Double I_ERI_Fx2y_S_Py_Pz_dd = I_ERI_Fx2y_S_Dyz_S_dd+CDZ*I_ERI_Fx2y_S_Py_S_dd;
  Double I_ERI_Fxyz_S_Py_Pz_dd = I_ERI_Fxyz_S_Dyz_S_dd+CDZ*I_ERI_Fxyz_S_Py_S_dd;
  Double I_ERI_Fx2z_S_Py_Pz_dd = I_ERI_Fx2z_S_Dyz_S_dd+CDZ*I_ERI_Fx2z_S_Py_S_dd;
  Double I_ERI_F3y_S_Py_Pz_dd = I_ERI_F3y_S_Dyz_S_dd+CDZ*I_ERI_F3y_S_Py_S_dd;
  Double I_ERI_F2yz_S_Py_Pz_dd = I_ERI_F2yz_S_Dyz_S_dd+CDZ*I_ERI_F2yz_S_Py_S_dd;
  Double I_ERI_Fy2z_S_Py_Pz_dd = I_ERI_Fy2z_S_Dyz_S_dd+CDZ*I_ERI_Fy2z_S_Py_S_dd;
  Double I_ERI_F3z_S_Py_Pz_dd = I_ERI_F3z_S_Dyz_S_dd+CDZ*I_ERI_F3z_S_Py_S_dd;
  Double I_ERI_F3x_S_Pz_Pz_dd = I_ERI_F3x_S_D2z_S_dd+CDZ*I_ERI_F3x_S_Pz_S_dd;
  Double I_ERI_F2xy_S_Pz_Pz_dd = I_ERI_F2xy_S_D2z_S_dd+CDZ*I_ERI_F2xy_S_Pz_S_dd;
  Double I_ERI_F2xz_S_Pz_Pz_dd = I_ERI_F2xz_S_D2z_S_dd+CDZ*I_ERI_F2xz_S_Pz_S_dd;
  Double I_ERI_Fx2y_S_Pz_Pz_dd = I_ERI_Fx2y_S_D2z_S_dd+CDZ*I_ERI_Fx2y_S_Pz_S_dd;
  Double I_ERI_Fxyz_S_Pz_Pz_dd = I_ERI_Fxyz_S_D2z_S_dd+CDZ*I_ERI_Fxyz_S_Pz_S_dd;
  Double I_ERI_Fx2z_S_Pz_Pz_dd = I_ERI_Fx2z_S_D2z_S_dd+CDZ*I_ERI_Fx2z_S_Pz_S_dd;
  Double I_ERI_F3y_S_Pz_Pz_dd = I_ERI_F3y_S_D2z_S_dd+CDZ*I_ERI_F3y_S_Pz_S_dd;
  Double I_ERI_F2yz_S_Pz_Pz_dd = I_ERI_F2yz_S_D2z_S_dd+CDZ*I_ERI_F2yz_S_Pz_S_dd;
  Double I_ERI_Fy2z_S_Pz_Pz_dd = I_ERI_Fy2z_S_D2z_S_dd+CDZ*I_ERI_Fy2z_S_Pz_S_dd;
  Double I_ERI_F3z_S_Pz_Pz_dd = I_ERI_F3z_S_D2z_S_dd+CDZ*I_ERI_F3z_S_Pz_S_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 40 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_dd
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_dd = I_ERI_F3x_S_F3x_S_dd+CDX*I_ERI_F3x_S_D2x_S_dd;
  Double I_ERI_F2xy_S_D2x_Px_dd = I_ERI_F2xy_S_F3x_S_dd+CDX*I_ERI_F2xy_S_D2x_S_dd;
  Double I_ERI_F2xz_S_D2x_Px_dd = I_ERI_F2xz_S_F3x_S_dd+CDX*I_ERI_F2xz_S_D2x_S_dd;
  Double I_ERI_Fx2y_S_D2x_Px_dd = I_ERI_Fx2y_S_F3x_S_dd+CDX*I_ERI_Fx2y_S_D2x_S_dd;
  Double I_ERI_Fxyz_S_D2x_Px_dd = I_ERI_Fxyz_S_F3x_S_dd+CDX*I_ERI_Fxyz_S_D2x_S_dd;
  Double I_ERI_Fx2z_S_D2x_Px_dd = I_ERI_Fx2z_S_F3x_S_dd+CDX*I_ERI_Fx2z_S_D2x_S_dd;
  Double I_ERI_F3y_S_D2x_Px_dd = I_ERI_F3y_S_F3x_S_dd+CDX*I_ERI_F3y_S_D2x_S_dd;
  Double I_ERI_F2yz_S_D2x_Px_dd = I_ERI_F2yz_S_F3x_S_dd+CDX*I_ERI_F2yz_S_D2x_S_dd;
  Double I_ERI_Fy2z_S_D2x_Px_dd = I_ERI_Fy2z_S_F3x_S_dd+CDX*I_ERI_Fy2z_S_D2x_S_dd;
  Double I_ERI_F3z_S_D2x_Px_dd = I_ERI_F3z_S_F3x_S_dd+CDX*I_ERI_F3z_S_D2x_S_dd;
  Double I_ERI_F3x_S_Dxy_Px_dd = I_ERI_F3x_S_F2xy_S_dd+CDX*I_ERI_F3x_S_Dxy_S_dd;
  Double I_ERI_F2xy_S_Dxy_Px_dd = I_ERI_F2xy_S_F2xy_S_dd+CDX*I_ERI_F2xy_S_Dxy_S_dd;
  Double I_ERI_F2xz_S_Dxy_Px_dd = I_ERI_F2xz_S_F2xy_S_dd+CDX*I_ERI_F2xz_S_Dxy_S_dd;
  Double I_ERI_Fx2y_S_Dxy_Px_dd = I_ERI_Fx2y_S_F2xy_S_dd+CDX*I_ERI_Fx2y_S_Dxy_S_dd;
  Double I_ERI_Fxyz_S_Dxy_Px_dd = I_ERI_Fxyz_S_F2xy_S_dd+CDX*I_ERI_Fxyz_S_Dxy_S_dd;
  Double I_ERI_Fx2z_S_Dxy_Px_dd = I_ERI_Fx2z_S_F2xy_S_dd+CDX*I_ERI_Fx2z_S_Dxy_S_dd;
  Double I_ERI_F3y_S_Dxy_Px_dd = I_ERI_F3y_S_F2xy_S_dd+CDX*I_ERI_F3y_S_Dxy_S_dd;
  Double I_ERI_F2yz_S_Dxy_Px_dd = I_ERI_F2yz_S_F2xy_S_dd+CDX*I_ERI_F2yz_S_Dxy_S_dd;
  Double I_ERI_Fy2z_S_Dxy_Px_dd = I_ERI_Fy2z_S_F2xy_S_dd+CDX*I_ERI_Fy2z_S_Dxy_S_dd;
  Double I_ERI_F3z_S_Dxy_Px_dd = I_ERI_F3z_S_F2xy_S_dd+CDX*I_ERI_F3z_S_Dxy_S_dd;
  Double I_ERI_F3x_S_Dxz_Px_dd = I_ERI_F3x_S_F2xz_S_dd+CDX*I_ERI_F3x_S_Dxz_S_dd;
  Double I_ERI_F2xy_S_Dxz_Px_dd = I_ERI_F2xy_S_F2xz_S_dd+CDX*I_ERI_F2xy_S_Dxz_S_dd;
  Double I_ERI_F2xz_S_Dxz_Px_dd = I_ERI_F2xz_S_F2xz_S_dd+CDX*I_ERI_F2xz_S_Dxz_S_dd;
  Double I_ERI_Fx2y_S_Dxz_Px_dd = I_ERI_Fx2y_S_F2xz_S_dd+CDX*I_ERI_Fx2y_S_Dxz_S_dd;
  Double I_ERI_Fxyz_S_Dxz_Px_dd = I_ERI_Fxyz_S_F2xz_S_dd+CDX*I_ERI_Fxyz_S_Dxz_S_dd;
  Double I_ERI_Fx2z_S_Dxz_Px_dd = I_ERI_Fx2z_S_F2xz_S_dd+CDX*I_ERI_Fx2z_S_Dxz_S_dd;
  Double I_ERI_F3y_S_Dxz_Px_dd = I_ERI_F3y_S_F2xz_S_dd+CDX*I_ERI_F3y_S_Dxz_S_dd;
  Double I_ERI_F2yz_S_Dxz_Px_dd = I_ERI_F2yz_S_F2xz_S_dd+CDX*I_ERI_F2yz_S_Dxz_S_dd;
  Double I_ERI_Fy2z_S_Dxz_Px_dd = I_ERI_Fy2z_S_F2xz_S_dd+CDX*I_ERI_Fy2z_S_Dxz_S_dd;
  Double I_ERI_F3z_S_Dxz_Px_dd = I_ERI_F3z_S_F2xz_S_dd+CDX*I_ERI_F3z_S_Dxz_S_dd;
  Double I_ERI_F3x_S_D2y_Px_dd = I_ERI_F3x_S_Fx2y_S_dd+CDX*I_ERI_F3x_S_D2y_S_dd;
  Double I_ERI_F2xy_S_D2y_Px_dd = I_ERI_F2xy_S_Fx2y_S_dd+CDX*I_ERI_F2xy_S_D2y_S_dd;
  Double I_ERI_F2xz_S_D2y_Px_dd = I_ERI_F2xz_S_Fx2y_S_dd+CDX*I_ERI_F2xz_S_D2y_S_dd;
  Double I_ERI_Fx2y_S_D2y_Px_dd = I_ERI_Fx2y_S_Fx2y_S_dd+CDX*I_ERI_Fx2y_S_D2y_S_dd;
  Double I_ERI_Fxyz_S_D2y_Px_dd = I_ERI_Fxyz_S_Fx2y_S_dd+CDX*I_ERI_Fxyz_S_D2y_S_dd;
  Double I_ERI_Fx2z_S_D2y_Px_dd = I_ERI_Fx2z_S_Fx2y_S_dd+CDX*I_ERI_Fx2z_S_D2y_S_dd;
  Double I_ERI_F3y_S_D2y_Px_dd = I_ERI_F3y_S_Fx2y_S_dd+CDX*I_ERI_F3y_S_D2y_S_dd;
  Double I_ERI_F2yz_S_D2y_Px_dd = I_ERI_F2yz_S_Fx2y_S_dd+CDX*I_ERI_F2yz_S_D2y_S_dd;
  Double I_ERI_Fy2z_S_D2y_Px_dd = I_ERI_Fy2z_S_Fx2y_S_dd+CDX*I_ERI_Fy2z_S_D2y_S_dd;
  Double I_ERI_F3z_S_D2y_Px_dd = I_ERI_F3z_S_Fx2y_S_dd+CDX*I_ERI_F3z_S_D2y_S_dd;
  Double I_ERI_F3x_S_Dyz_Px_dd = I_ERI_F3x_S_Fxyz_S_dd+CDX*I_ERI_F3x_S_Dyz_S_dd;
  Double I_ERI_F2xy_S_Dyz_Px_dd = I_ERI_F2xy_S_Fxyz_S_dd+CDX*I_ERI_F2xy_S_Dyz_S_dd;
  Double I_ERI_F2xz_S_Dyz_Px_dd = I_ERI_F2xz_S_Fxyz_S_dd+CDX*I_ERI_F2xz_S_Dyz_S_dd;
  Double I_ERI_Fx2y_S_Dyz_Px_dd = I_ERI_Fx2y_S_Fxyz_S_dd+CDX*I_ERI_Fx2y_S_Dyz_S_dd;
  Double I_ERI_Fxyz_S_Dyz_Px_dd = I_ERI_Fxyz_S_Fxyz_S_dd+CDX*I_ERI_Fxyz_S_Dyz_S_dd;
  Double I_ERI_Fx2z_S_Dyz_Px_dd = I_ERI_Fx2z_S_Fxyz_S_dd+CDX*I_ERI_Fx2z_S_Dyz_S_dd;
  Double I_ERI_F3y_S_Dyz_Px_dd = I_ERI_F3y_S_Fxyz_S_dd+CDX*I_ERI_F3y_S_Dyz_S_dd;
  Double I_ERI_F2yz_S_Dyz_Px_dd = I_ERI_F2yz_S_Fxyz_S_dd+CDX*I_ERI_F2yz_S_Dyz_S_dd;
  Double I_ERI_Fy2z_S_Dyz_Px_dd = I_ERI_Fy2z_S_Fxyz_S_dd+CDX*I_ERI_Fy2z_S_Dyz_S_dd;
  Double I_ERI_F3z_S_Dyz_Px_dd = I_ERI_F3z_S_Fxyz_S_dd+CDX*I_ERI_F3z_S_Dyz_S_dd;
  Double I_ERI_F3x_S_D2z_Px_dd = I_ERI_F3x_S_Fx2z_S_dd+CDX*I_ERI_F3x_S_D2z_S_dd;
  Double I_ERI_F2xy_S_D2z_Px_dd = I_ERI_F2xy_S_Fx2z_S_dd+CDX*I_ERI_F2xy_S_D2z_S_dd;
  Double I_ERI_F2xz_S_D2z_Px_dd = I_ERI_F2xz_S_Fx2z_S_dd+CDX*I_ERI_F2xz_S_D2z_S_dd;
  Double I_ERI_Fx2y_S_D2z_Px_dd = I_ERI_Fx2y_S_Fx2z_S_dd+CDX*I_ERI_Fx2y_S_D2z_S_dd;
  Double I_ERI_Fxyz_S_D2z_Px_dd = I_ERI_Fxyz_S_Fx2z_S_dd+CDX*I_ERI_Fxyz_S_D2z_S_dd;
  Double I_ERI_Fx2z_S_D2z_Px_dd = I_ERI_Fx2z_S_Fx2z_S_dd+CDX*I_ERI_Fx2z_S_D2z_S_dd;
  Double I_ERI_F3y_S_D2z_Px_dd = I_ERI_F3y_S_Fx2z_S_dd+CDX*I_ERI_F3y_S_D2z_S_dd;
  Double I_ERI_F2yz_S_D2z_Px_dd = I_ERI_F2yz_S_Fx2z_S_dd+CDX*I_ERI_F2yz_S_D2z_S_dd;
  Double I_ERI_Fy2z_S_D2z_Px_dd = I_ERI_Fy2z_S_Fx2z_S_dd+CDX*I_ERI_Fy2z_S_D2z_S_dd;
  Double I_ERI_F3z_S_D2z_Px_dd = I_ERI_F3z_S_Fx2z_S_dd+CDX*I_ERI_F3z_S_D2z_S_dd;
  Double I_ERI_F3x_S_Dxy_Py_dd = I_ERI_F3x_S_Fx2y_S_dd+CDY*I_ERI_F3x_S_Dxy_S_dd;
  Double I_ERI_F2xy_S_Dxy_Py_dd = I_ERI_F2xy_S_Fx2y_S_dd+CDY*I_ERI_F2xy_S_Dxy_S_dd;
  Double I_ERI_F2xz_S_Dxy_Py_dd = I_ERI_F2xz_S_Fx2y_S_dd+CDY*I_ERI_F2xz_S_Dxy_S_dd;
  Double I_ERI_Fx2y_S_Dxy_Py_dd = I_ERI_Fx2y_S_Fx2y_S_dd+CDY*I_ERI_Fx2y_S_Dxy_S_dd;
  Double I_ERI_Fxyz_S_Dxy_Py_dd = I_ERI_Fxyz_S_Fx2y_S_dd+CDY*I_ERI_Fxyz_S_Dxy_S_dd;
  Double I_ERI_Fx2z_S_Dxy_Py_dd = I_ERI_Fx2z_S_Fx2y_S_dd+CDY*I_ERI_Fx2z_S_Dxy_S_dd;
  Double I_ERI_F3y_S_Dxy_Py_dd = I_ERI_F3y_S_Fx2y_S_dd+CDY*I_ERI_F3y_S_Dxy_S_dd;
  Double I_ERI_F2yz_S_Dxy_Py_dd = I_ERI_F2yz_S_Fx2y_S_dd+CDY*I_ERI_F2yz_S_Dxy_S_dd;
  Double I_ERI_Fy2z_S_Dxy_Py_dd = I_ERI_Fy2z_S_Fx2y_S_dd+CDY*I_ERI_Fy2z_S_Dxy_S_dd;
  Double I_ERI_F3z_S_Dxy_Py_dd = I_ERI_F3z_S_Fx2y_S_dd+CDY*I_ERI_F3z_S_Dxy_S_dd;
  Double I_ERI_F3x_S_Dxz_Py_dd = I_ERI_F3x_S_Fxyz_S_dd+CDY*I_ERI_F3x_S_Dxz_S_dd;
  Double I_ERI_F2xy_S_Dxz_Py_dd = I_ERI_F2xy_S_Fxyz_S_dd+CDY*I_ERI_F2xy_S_Dxz_S_dd;
  Double I_ERI_F2xz_S_Dxz_Py_dd = I_ERI_F2xz_S_Fxyz_S_dd+CDY*I_ERI_F2xz_S_Dxz_S_dd;
  Double I_ERI_Fx2y_S_Dxz_Py_dd = I_ERI_Fx2y_S_Fxyz_S_dd+CDY*I_ERI_Fx2y_S_Dxz_S_dd;
  Double I_ERI_Fxyz_S_Dxz_Py_dd = I_ERI_Fxyz_S_Fxyz_S_dd+CDY*I_ERI_Fxyz_S_Dxz_S_dd;
  Double I_ERI_Fx2z_S_Dxz_Py_dd = I_ERI_Fx2z_S_Fxyz_S_dd+CDY*I_ERI_Fx2z_S_Dxz_S_dd;
  Double I_ERI_F3y_S_Dxz_Py_dd = I_ERI_F3y_S_Fxyz_S_dd+CDY*I_ERI_F3y_S_Dxz_S_dd;
  Double I_ERI_F2yz_S_Dxz_Py_dd = I_ERI_F2yz_S_Fxyz_S_dd+CDY*I_ERI_F2yz_S_Dxz_S_dd;
  Double I_ERI_Fy2z_S_Dxz_Py_dd = I_ERI_Fy2z_S_Fxyz_S_dd+CDY*I_ERI_Fy2z_S_Dxz_S_dd;
  Double I_ERI_F3z_S_Dxz_Py_dd = I_ERI_F3z_S_Fxyz_S_dd+CDY*I_ERI_F3z_S_Dxz_S_dd;
  Double I_ERI_F3x_S_D2y_Py_dd = I_ERI_F3x_S_F3y_S_dd+CDY*I_ERI_F3x_S_D2y_S_dd;
  Double I_ERI_F2xy_S_D2y_Py_dd = I_ERI_F2xy_S_F3y_S_dd+CDY*I_ERI_F2xy_S_D2y_S_dd;
  Double I_ERI_F2xz_S_D2y_Py_dd = I_ERI_F2xz_S_F3y_S_dd+CDY*I_ERI_F2xz_S_D2y_S_dd;
  Double I_ERI_Fx2y_S_D2y_Py_dd = I_ERI_Fx2y_S_F3y_S_dd+CDY*I_ERI_Fx2y_S_D2y_S_dd;
  Double I_ERI_Fxyz_S_D2y_Py_dd = I_ERI_Fxyz_S_F3y_S_dd+CDY*I_ERI_Fxyz_S_D2y_S_dd;
  Double I_ERI_Fx2z_S_D2y_Py_dd = I_ERI_Fx2z_S_F3y_S_dd+CDY*I_ERI_Fx2z_S_D2y_S_dd;
  Double I_ERI_F3y_S_D2y_Py_dd = I_ERI_F3y_S_F3y_S_dd+CDY*I_ERI_F3y_S_D2y_S_dd;
  Double I_ERI_F2yz_S_D2y_Py_dd = I_ERI_F2yz_S_F3y_S_dd+CDY*I_ERI_F2yz_S_D2y_S_dd;
  Double I_ERI_Fy2z_S_D2y_Py_dd = I_ERI_Fy2z_S_F3y_S_dd+CDY*I_ERI_Fy2z_S_D2y_S_dd;
  Double I_ERI_F3z_S_D2y_Py_dd = I_ERI_F3z_S_F3y_S_dd+CDY*I_ERI_F3z_S_D2y_S_dd;
  Double I_ERI_F3x_S_Dyz_Py_dd = I_ERI_F3x_S_F2yz_S_dd+CDY*I_ERI_F3x_S_Dyz_S_dd;
  Double I_ERI_F2xy_S_Dyz_Py_dd = I_ERI_F2xy_S_F2yz_S_dd+CDY*I_ERI_F2xy_S_Dyz_S_dd;
  Double I_ERI_F2xz_S_Dyz_Py_dd = I_ERI_F2xz_S_F2yz_S_dd+CDY*I_ERI_F2xz_S_Dyz_S_dd;
  Double I_ERI_Fx2y_S_Dyz_Py_dd = I_ERI_Fx2y_S_F2yz_S_dd+CDY*I_ERI_Fx2y_S_Dyz_S_dd;
  Double I_ERI_Fxyz_S_Dyz_Py_dd = I_ERI_Fxyz_S_F2yz_S_dd+CDY*I_ERI_Fxyz_S_Dyz_S_dd;
  Double I_ERI_Fx2z_S_Dyz_Py_dd = I_ERI_Fx2z_S_F2yz_S_dd+CDY*I_ERI_Fx2z_S_Dyz_S_dd;
  Double I_ERI_F3y_S_Dyz_Py_dd = I_ERI_F3y_S_F2yz_S_dd+CDY*I_ERI_F3y_S_Dyz_S_dd;
  Double I_ERI_F2yz_S_Dyz_Py_dd = I_ERI_F2yz_S_F2yz_S_dd+CDY*I_ERI_F2yz_S_Dyz_S_dd;
  Double I_ERI_Fy2z_S_Dyz_Py_dd = I_ERI_Fy2z_S_F2yz_S_dd+CDY*I_ERI_Fy2z_S_Dyz_S_dd;
  Double I_ERI_F3z_S_Dyz_Py_dd = I_ERI_F3z_S_F2yz_S_dd+CDY*I_ERI_F3z_S_Dyz_S_dd;
  Double I_ERI_F3x_S_D2z_Py_dd = I_ERI_F3x_S_Fy2z_S_dd+CDY*I_ERI_F3x_S_D2z_S_dd;
  Double I_ERI_F2xy_S_D2z_Py_dd = I_ERI_F2xy_S_Fy2z_S_dd+CDY*I_ERI_F2xy_S_D2z_S_dd;
  Double I_ERI_F2xz_S_D2z_Py_dd = I_ERI_F2xz_S_Fy2z_S_dd+CDY*I_ERI_F2xz_S_D2z_S_dd;
  Double I_ERI_Fx2y_S_D2z_Py_dd = I_ERI_Fx2y_S_Fy2z_S_dd+CDY*I_ERI_Fx2y_S_D2z_S_dd;
  Double I_ERI_Fxyz_S_D2z_Py_dd = I_ERI_Fxyz_S_Fy2z_S_dd+CDY*I_ERI_Fxyz_S_D2z_S_dd;
  Double I_ERI_Fx2z_S_D2z_Py_dd = I_ERI_Fx2z_S_Fy2z_S_dd+CDY*I_ERI_Fx2z_S_D2z_S_dd;
  Double I_ERI_F3y_S_D2z_Py_dd = I_ERI_F3y_S_Fy2z_S_dd+CDY*I_ERI_F3y_S_D2z_S_dd;
  Double I_ERI_F2yz_S_D2z_Py_dd = I_ERI_F2yz_S_Fy2z_S_dd+CDY*I_ERI_F2yz_S_D2z_S_dd;
  Double I_ERI_Fy2z_S_D2z_Py_dd = I_ERI_Fy2z_S_Fy2z_S_dd+CDY*I_ERI_Fy2z_S_D2z_S_dd;
  Double I_ERI_F3z_S_D2z_Py_dd = I_ERI_F3z_S_Fy2z_S_dd+CDY*I_ERI_F3z_S_D2z_S_dd;
  Double I_ERI_F3x_S_Dxz_Pz_dd = I_ERI_F3x_S_Fx2z_S_dd+CDZ*I_ERI_F3x_S_Dxz_S_dd;
  Double I_ERI_F2xy_S_Dxz_Pz_dd = I_ERI_F2xy_S_Fx2z_S_dd+CDZ*I_ERI_F2xy_S_Dxz_S_dd;
  Double I_ERI_F2xz_S_Dxz_Pz_dd = I_ERI_F2xz_S_Fx2z_S_dd+CDZ*I_ERI_F2xz_S_Dxz_S_dd;
  Double I_ERI_Fx2y_S_Dxz_Pz_dd = I_ERI_Fx2y_S_Fx2z_S_dd+CDZ*I_ERI_Fx2y_S_Dxz_S_dd;
  Double I_ERI_Fxyz_S_Dxz_Pz_dd = I_ERI_Fxyz_S_Fx2z_S_dd+CDZ*I_ERI_Fxyz_S_Dxz_S_dd;
  Double I_ERI_Fx2z_S_Dxz_Pz_dd = I_ERI_Fx2z_S_Fx2z_S_dd+CDZ*I_ERI_Fx2z_S_Dxz_S_dd;
  Double I_ERI_F3y_S_Dxz_Pz_dd = I_ERI_F3y_S_Fx2z_S_dd+CDZ*I_ERI_F3y_S_Dxz_S_dd;
  Double I_ERI_F2yz_S_Dxz_Pz_dd = I_ERI_F2yz_S_Fx2z_S_dd+CDZ*I_ERI_F2yz_S_Dxz_S_dd;
  Double I_ERI_Fy2z_S_Dxz_Pz_dd = I_ERI_Fy2z_S_Fx2z_S_dd+CDZ*I_ERI_Fy2z_S_Dxz_S_dd;
  Double I_ERI_F3z_S_Dxz_Pz_dd = I_ERI_F3z_S_Fx2z_S_dd+CDZ*I_ERI_F3z_S_Dxz_S_dd;
  Double I_ERI_F3x_S_Dyz_Pz_dd = I_ERI_F3x_S_Fy2z_S_dd+CDZ*I_ERI_F3x_S_Dyz_S_dd;
  Double I_ERI_F2xy_S_Dyz_Pz_dd = I_ERI_F2xy_S_Fy2z_S_dd+CDZ*I_ERI_F2xy_S_Dyz_S_dd;
  Double I_ERI_F2xz_S_Dyz_Pz_dd = I_ERI_F2xz_S_Fy2z_S_dd+CDZ*I_ERI_F2xz_S_Dyz_S_dd;
  Double I_ERI_Fx2y_S_Dyz_Pz_dd = I_ERI_Fx2y_S_Fy2z_S_dd+CDZ*I_ERI_Fx2y_S_Dyz_S_dd;
  Double I_ERI_Fxyz_S_Dyz_Pz_dd = I_ERI_Fxyz_S_Fy2z_S_dd+CDZ*I_ERI_Fxyz_S_Dyz_S_dd;
  Double I_ERI_Fx2z_S_Dyz_Pz_dd = I_ERI_Fx2z_S_Fy2z_S_dd+CDZ*I_ERI_Fx2z_S_Dyz_S_dd;
  Double I_ERI_F3y_S_Dyz_Pz_dd = I_ERI_F3y_S_Fy2z_S_dd+CDZ*I_ERI_F3y_S_Dyz_S_dd;
  Double I_ERI_F2yz_S_Dyz_Pz_dd = I_ERI_F2yz_S_Fy2z_S_dd+CDZ*I_ERI_F2yz_S_Dyz_S_dd;
  Double I_ERI_Fy2z_S_Dyz_Pz_dd = I_ERI_Fy2z_S_Fy2z_S_dd+CDZ*I_ERI_Fy2z_S_Dyz_S_dd;
  Double I_ERI_F3z_S_Dyz_Pz_dd = I_ERI_F3z_S_Fy2z_S_dd+CDZ*I_ERI_F3z_S_Dyz_S_dd;
  Double I_ERI_F3x_S_D2z_Pz_dd = I_ERI_F3x_S_F3z_S_dd+CDZ*I_ERI_F3x_S_D2z_S_dd;
  Double I_ERI_F2xy_S_D2z_Pz_dd = I_ERI_F2xy_S_F3z_S_dd+CDZ*I_ERI_F2xy_S_D2z_S_dd;
  Double I_ERI_F2xz_S_D2z_Pz_dd = I_ERI_F2xz_S_F3z_S_dd+CDZ*I_ERI_F2xz_S_D2z_S_dd;
  Double I_ERI_Fx2y_S_D2z_Pz_dd = I_ERI_Fx2y_S_F3z_S_dd+CDZ*I_ERI_Fx2y_S_D2z_S_dd;
  Double I_ERI_Fxyz_S_D2z_Pz_dd = I_ERI_Fxyz_S_F3z_S_dd+CDZ*I_ERI_Fxyz_S_D2z_S_dd;
  Double I_ERI_Fx2z_S_D2z_Pz_dd = I_ERI_Fx2z_S_F3z_S_dd+CDZ*I_ERI_Fx2z_S_D2z_S_dd;
  Double I_ERI_F3y_S_D2z_Pz_dd = I_ERI_F3y_S_F3z_S_dd+CDZ*I_ERI_F3y_S_D2z_S_dd;
  Double I_ERI_F2yz_S_D2z_Pz_dd = I_ERI_F2yz_S_F3z_S_dd+CDZ*I_ERI_F2yz_S_D2z_S_dd;
  Double I_ERI_Fy2z_S_D2z_Pz_dd = I_ERI_Fy2z_S_F3z_S_dd+CDZ*I_ERI_Fy2z_S_D2z_S_dd;
  Double I_ERI_F3z_S_D2z_Pz_dd = I_ERI_F3z_S_F3z_S_dd+CDZ*I_ERI_F3z_S_D2z_S_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_D_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_dd
   ************************************************************/
  Double I_ERI_F3x_S_Px_D2x_dd = I_ERI_F3x_S_D2x_Px_dd+CDX*I_ERI_F3x_S_Px_Px_dd;
  Double I_ERI_F2xy_S_Px_D2x_dd = I_ERI_F2xy_S_D2x_Px_dd+CDX*I_ERI_F2xy_S_Px_Px_dd;
  Double I_ERI_F2xz_S_Px_D2x_dd = I_ERI_F2xz_S_D2x_Px_dd+CDX*I_ERI_F2xz_S_Px_Px_dd;
  Double I_ERI_Fx2y_S_Px_D2x_dd = I_ERI_Fx2y_S_D2x_Px_dd+CDX*I_ERI_Fx2y_S_Px_Px_dd;
  Double I_ERI_Fxyz_S_Px_D2x_dd = I_ERI_Fxyz_S_D2x_Px_dd+CDX*I_ERI_Fxyz_S_Px_Px_dd;
  Double I_ERI_Fx2z_S_Px_D2x_dd = I_ERI_Fx2z_S_D2x_Px_dd+CDX*I_ERI_Fx2z_S_Px_Px_dd;
  Double I_ERI_F3y_S_Px_D2x_dd = I_ERI_F3y_S_D2x_Px_dd+CDX*I_ERI_F3y_S_Px_Px_dd;
  Double I_ERI_F2yz_S_Px_D2x_dd = I_ERI_F2yz_S_D2x_Px_dd+CDX*I_ERI_F2yz_S_Px_Px_dd;
  Double I_ERI_Fy2z_S_Px_D2x_dd = I_ERI_Fy2z_S_D2x_Px_dd+CDX*I_ERI_Fy2z_S_Px_Px_dd;
  Double I_ERI_F3z_S_Px_D2x_dd = I_ERI_F3z_S_D2x_Px_dd+CDX*I_ERI_F3z_S_Px_Px_dd;
  Double I_ERI_F3x_S_Py_D2x_dd = I_ERI_F3x_S_Dxy_Px_dd+CDX*I_ERI_F3x_S_Py_Px_dd;
  Double I_ERI_F2xy_S_Py_D2x_dd = I_ERI_F2xy_S_Dxy_Px_dd+CDX*I_ERI_F2xy_S_Py_Px_dd;
  Double I_ERI_F2xz_S_Py_D2x_dd = I_ERI_F2xz_S_Dxy_Px_dd+CDX*I_ERI_F2xz_S_Py_Px_dd;
  Double I_ERI_Fx2y_S_Py_D2x_dd = I_ERI_Fx2y_S_Dxy_Px_dd+CDX*I_ERI_Fx2y_S_Py_Px_dd;
  Double I_ERI_Fxyz_S_Py_D2x_dd = I_ERI_Fxyz_S_Dxy_Px_dd+CDX*I_ERI_Fxyz_S_Py_Px_dd;
  Double I_ERI_Fx2z_S_Py_D2x_dd = I_ERI_Fx2z_S_Dxy_Px_dd+CDX*I_ERI_Fx2z_S_Py_Px_dd;
  Double I_ERI_F3y_S_Py_D2x_dd = I_ERI_F3y_S_Dxy_Px_dd+CDX*I_ERI_F3y_S_Py_Px_dd;
  Double I_ERI_F2yz_S_Py_D2x_dd = I_ERI_F2yz_S_Dxy_Px_dd+CDX*I_ERI_F2yz_S_Py_Px_dd;
  Double I_ERI_Fy2z_S_Py_D2x_dd = I_ERI_Fy2z_S_Dxy_Px_dd+CDX*I_ERI_Fy2z_S_Py_Px_dd;
  Double I_ERI_F3z_S_Py_D2x_dd = I_ERI_F3z_S_Dxy_Px_dd+CDX*I_ERI_F3z_S_Py_Px_dd;
  Double I_ERI_F3x_S_Pz_D2x_dd = I_ERI_F3x_S_Dxz_Px_dd+CDX*I_ERI_F3x_S_Pz_Px_dd;
  Double I_ERI_F2xy_S_Pz_D2x_dd = I_ERI_F2xy_S_Dxz_Px_dd+CDX*I_ERI_F2xy_S_Pz_Px_dd;
  Double I_ERI_F2xz_S_Pz_D2x_dd = I_ERI_F2xz_S_Dxz_Px_dd+CDX*I_ERI_F2xz_S_Pz_Px_dd;
  Double I_ERI_Fx2y_S_Pz_D2x_dd = I_ERI_Fx2y_S_Dxz_Px_dd+CDX*I_ERI_Fx2y_S_Pz_Px_dd;
  Double I_ERI_Fxyz_S_Pz_D2x_dd = I_ERI_Fxyz_S_Dxz_Px_dd+CDX*I_ERI_Fxyz_S_Pz_Px_dd;
  Double I_ERI_Fx2z_S_Pz_D2x_dd = I_ERI_Fx2z_S_Dxz_Px_dd+CDX*I_ERI_Fx2z_S_Pz_Px_dd;
  Double I_ERI_F3y_S_Pz_D2x_dd = I_ERI_F3y_S_Dxz_Px_dd+CDX*I_ERI_F3y_S_Pz_Px_dd;
  Double I_ERI_F2yz_S_Pz_D2x_dd = I_ERI_F2yz_S_Dxz_Px_dd+CDX*I_ERI_F2yz_S_Pz_Px_dd;
  Double I_ERI_Fy2z_S_Pz_D2x_dd = I_ERI_Fy2z_S_Dxz_Px_dd+CDX*I_ERI_Fy2z_S_Pz_Px_dd;
  Double I_ERI_F3z_S_Pz_D2x_dd = I_ERI_F3z_S_Dxz_Px_dd+CDX*I_ERI_F3z_S_Pz_Px_dd;
  Double I_ERI_F3x_S_Px_Dxy_dd = I_ERI_F3x_S_Dxy_Px_dd+CDY*I_ERI_F3x_S_Px_Px_dd;
  Double I_ERI_F2xy_S_Px_Dxy_dd = I_ERI_F2xy_S_Dxy_Px_dd+CDY*I_ERI_F2xy_S_Px_Px_dd;
  Double I_ERI_F2xz_S_Px_Dxy_dd = I_ERI_F2xz_S_Dxy_Px_dd+CDY*I_ERI_F2xz_S_Px_Px_dd;
  Double I_ERI_Fx2y_S_Px_Dxy_dd = I_ERI_Fx2y_S_Dxy_Px_dd+CDY*I_ERI_Fx2y_S_Px_Px_dd;
  Double I_ERI_Fxyz_S_Px_Dxy_dd = I_ERI_Fxyz_S_Dxy_Px_dd+CDY*I_ERI_Fxyz_S_Px_Px_dd;
  Double I_ERI_Fx2z_S_Px_Dxy_dd = I_ERI_Fx2z_S_Dxy_Px_dd+CDY*I_ERI_Fx2z_S_Px_Px_dd;
  Double I_ERI_F3y_S_Px_Dxy_dd = I_ERI_F3y_S_Dxy_Px_dd+CDY*I_ERI_F3y_S_Px_Px_dd;
  Double I_ERI_F2yz_S_Px_Dxy_dd = I_ERI_F2yz_S_Dxy_Px_dd+CDY*I_ERI_F2yz_S_Px_Px_dd;
  Double I_ERI_Fy2z_S_Px_Dxy_dd = I_ERI_Fy2z_S_Dxy_Px_dd+CDY*I_ERI_Fy2z_S_Px_Px_dd;
  Double I_ERI_F3z_S_Px_Dxy_dd = I_ERI_F3z_S_Dxy_Px_dd+CDY*I_ERI_F3z_S_Px_Px_dd;
  Double I_ERI_F3x_S_Py_Dxy_dd = I_ERI_F3x_S_D2y_Px_dd+CDY*I_ERI_F3x_S_Py_Px_dd;
  Double I_ERI_F2xy_S_Py_Dxy_dd = I_ERI_F2xy_S_D2y_Px_dd+CDY*I_ERI_F2xy_S_Py_Px_dd;
  Double I_ERI_F2xz_S_Py_Dxy_dd = I_ERI_F2xz_S_D2y_Px_dd+CDY*I_ERI_F2xz_S_Py_Px_dd;
  Double I_ERI_Fx2y_S_Py_Dxy_dd = I_ERI_Fx2y_S_D2y_Px_dd+CDY*I_ERI_Fx2y_S_Py_Px_dd;
  Double I_ERI_Fxyz_S_Py_Dxy_dd = I_ERI_Fxyz_S_D2y_Px_dd+CDY*I_ERI_Fxyz_S_Py_Px_dd;
  Double I_ERI_Fx2z_S_Py_Dxy_dd = I_ERI_Fx2z_S_D2y_Px_dd+CDY*I_ERI_Fx2z_S_Py_Px_dd;
  Double I_ERI_F3y_S_Py_Dxy_dd = I_ERI_F3y_S_D2y_Px_dd+CDY*I_ERI_F3y_S_Py_Px_dd;
  Double I_ERI_F2yz_S_Py_Dxy_dd = I_ERI_F2yz_S_D2y_Px_dd+CDY*I_ERI_F2yz_S_Py_Px_dd;
  Double I_ERI_Fy2z_S_Py_Dxy_dd = I_ERI_Fy2z_S_D2y_Px_dd+CDY*I_ERI_Fy2z_S_Py_Px_dd;
  Double I_ERI_F3z_S_Py_Dxy_dd = I_ERI_F3z_S_D2y_Px_dd+CDY*I_ERI_F3z_S_Py_Px_dd;
  Double I_ERI_F3x_S_Pz_Dxy_dd = I_ERI_F3x_S_Dyz_Px_dd+CDY*I_ERI_F3x_S_Pz_Px_dd;
  Double I_ERI_F2xy_S_Pz_Dxy_dd = I_ERI_F2xy_S_Dyz_Px_dd+CDY*I_ERI_F2xy_S_Pz_Px_dd;
  Double I_ERI_F2xz_S_Pz_Dxy_dd = I_ERI_F2xz_S_Dyz_Px_dd+CDY*I_ERI_F2xz_S_Pz_Px_dd;
  Double I_ERI_Fx2y_S_Pz_Dxy_dd = I_ERI_Fx2y_S_Dyz_Px_dd+CDY*I_ERI_Fx2y_S_Pz_Px_dd;
  Double I_ERI_Fxyz_S_Pz_Dxy_dd = I_ERI_Fxyz_S_Dyz_Px_dd+CDY*I_ERI_Fxyz_S_Pz_Px_dd;
  Double I_ERI_Fx2z_S_Pz_Dxy_dd = I_ERI_Fx2z_S_Dyz_Px_dd+CDY*I_ERI_Fx2z_S_Pz_Px_dd;
  Double I_ERI_F3y_S_Pz_Dxy_dd = I_ERI_F3y_S_Dyz_Px_dd+CDY*I_ERI_F3y_S_Pz_Px_dd;
  Double I_ERI_F2yz_S_Pz_Dxy_dd = I_ERI_F2yz_S_Dyz_Px_dd+CDY*I_ERI_F2yz_S_Pz_Px_dd;
  Double I_ERI_Fy2z_S_Pz_Dxy_dd = I_ERI_Fy2z_S_Dyz_Px_dd+CDY*I_ERI_Fy2z_S_Pz_Px_dd;
  Double I_ERI_F3z_S_Pz_Dxy_dd = I_ERI_F3z_S_Dyz_Px_dd+CDY*I_ERI_F3z_S_Pz_Px_dd;
  Double I_ERI_F3x_S_Px_Dxz_dd = I_ERI_F3x_S_Dxz_Px_dd+CDZ*I_ERI_F3x_S_Px_Px_dd;
  Double I_ERI_F2xy_S_Px_Dxz_dd = I_ERI_F2xy_S_Dxz_Px_dd+CDZ*I_ERI_F2xy_S_Px_Px_dd;
  Double I_ERI_F2xz_S_Px_Dxz_dd = I_ERI_F2xz_S_Dxz_Px_dd+CDZ*I_ERI_F2xz_S_Px_Px_dd;
  Double I_ERI_Fx2y_S_Px_Dxz_dd = I_ERI_Fx2y_S_Dxz_Px_dd+CDZ*I_ERI_Fx2y_S_Px_Px_dd;
  Double I_ERI_Fxyz_S_Px_Dxz_dd = I_ERI_Fxyz_S_Dxz_Px_dd+CDZ*I_ERI_Fxyz_S_Px_Px_dd;
  Double I_ERI_Fx2z_S_Px_Dxz_dd = I_ERI_Fx2z_S_Dxz_Px_dd+CDZ*I_ERI_Fx2z_S_Px_Px_dd;
  Double I_ERI_F3y_S_Px_Dxz_dd = I_ERI_F3y_S_Dxz_Px_dd+CDZ*I_ERI_F3y_S_Px_Px_dd;
  Double I_ERI_F2yz_S_Px_Dxz_dd = I_ERI_F2yz_S_Dxz_Px_dd+CDZ*I_ERI_F2yz_S_Px_Px_dd;
  Double I_ERI_Fy2z_S_Px_Dxz_dd = I_ERI_Fy2z_S_Dxz_Px_dd+CDZ*I_ERI_Fy2z_S_Px_Px_dd;
  Double I_ERI_F3z_S_Px_Dxz_dd = I_ERI_F3z_S_Dxz_Px_dd+CDZ*I_ERI_F3z_S_Px_Px_dd;
  Double I_ERI_F3x_S_Py_Dxz_dd = I_ERI_F3x_S_Dyz_Px_dd+CDZ*I_ERI_F3x_S_Py_Px_dd;
  Double I_ERI_F2xy_S_Py_Dxz_dd = I_ERI_F2xy_S_Dyz_Px_dd+CDZ*I_ERI_F2xy_S_Py_Px_dd;
  Double I_ERI_F2xz_S_Py_Dxz_dd = I_ERI_F2xz_S_Dyz_Px_dd+CDZ*I_ERI_F2xz_S_Py_Px_dd;
  Double I_ERI_Fx2y_S_Py_Dxz_dd = I_ERI_Fx2y_S_Dyz_Px_dd+CDZ*I_ERI_Fx2y_S_Py_Px_dd;
  Double I_ERI_Fxyz_S_Py_Dxz_dd = I_ERI_Fxyz_S_Dyz_Px_dd+CDZ*I_ERI_Fxyz_S_Py_Px_dd;
  Double I_ERI_Fx2z_S_Py_Dxz_dd = I_ERI_Fx2z_S_Dyz_Px_dd+CDZ*I_ERI_Fx2z_S_Py_Px_dd;
  Double I_ERI_F3y_S_Py_Dxz_dd = I_ERI_F3y_S_Dyz_Px_dd+CDZ*I_ERI_F3y_S_Py_Px_dd;
  Double I_ERI_F2yz_S_Py_Dxz_dd = I_ERI_F2yz_S_Dyz_Px_dd+CDZ*I_ERI_F2yz_S_Py_Px_dd;
  Double I_ERI_Fy2z_S_Py_Dxz_dd = I_ERI_Fy2z_S_Dyz_Px_dd+CDZ*I_ERI_Fy2z_S_Py_Px_dd;
  Double I_ERI_F3z_S_Py_Dxz_dd = I_ERI_F3z_S_Dyz_Px_dd+CDZ*I_ERI_F3z_S_Py_Px_dd;
  Double I_ERI_F3x_S_Pz_Dxz_dd = I_ERI_F3x_S_D2z_Px_dd+CDZ*I_ERI_F3x_S_Pz_Px_dd;
  Double I_ERI_F2xy_S_Pz_Dxz_dd = I_ERI_F2xy_S_D2z_Px_dd+CDZ*I_ERI_F2xy_S_Pz_Px_dd;
  Double I_ERI_F2xz_S_Pz_Dxz_dd = I_ERI_F2xz_S_D2z_Px_dd+CDZ*I_ERI_F2xz_S_Pz_Px_dd;
  Double I_ERI_Fx2y_S_Pz_Dxz_dd = I_ERI_Fx2y_S_D2z_Px_dd+CDZ*I_ERI_Fx2y_S_Pz_Px_dd;
  Double I_ERI_Fxyz_S_Pz_Dxz_dd = I_ERI_Fxyz_S_D2z_Px_dd+CDZ*I_ERI_Fxyz_S_Pz_Px_dd;
  Double I_ERI_Fx2z_S_Pz_Dxz_dd = I_ERI_Fx2z_S_D2z_Px_dd+CDZ*I_ERI_Fx2z_S_Pz_Px_dd;
  Double I_ERI_F3y_S_Pz_Dxz_dd = I_ERI_F3y_S_D2z_Px_dd+CDZ*I_ERI_F3y_S_Pz_Px_dd;
  Double I_ERI_F2yz_S_Pz_Dxz_dd = I_ERI_F2yz_S_D2z_Px_dd+CDZ*I_ERI_F2yz_S_Pz_Px_dd;
  Double I_ERI_Fy2z_S_Pz_Dxz_dd = I_ERI_Fy2z_S_D2z_Px_dd+CDZ*I_ERI_Fy2z_S_Pz_Px_dd;
  Double I_ERI_F3z_S_Pz_Dxz_dd = I_ERI_F3z_S_D2z_Px_dd+CDZ*I_ERI_F3z_S_Pz_Px_dd;
  Double I_ERI_F3x_S_Px_D2y_dd = I_ERI_F3x_S_Dxy_Py_dd+CDY*I_ERI_F3x_S_Px_Py_dd;
  Double I_ERI_F2xy_S_Px_D2y_dd = I_ERI_F2xy_S_Dxy_Py_dd+CDY*I_ERI_F2xy_S_Px_Py_dd;
  Double I_ERI_F2xz_S_Px_D2y_dd = I_ERI_F2xz_S_Dxy_Py_dd+CDY*I_ERI_F2xz_S_Px_Py_dd;
  Double I_ERI_Fx2y_S_Px_D2y_dd = I_ERI_Fx2y_S_Dxy_Py_dd+CDY*I_ERI_Fx2y_S_Px_Py_dd;
  Double I_ERI_Fxyz_S_Px_D2y_dd = I_ERI_Fxyz_S_Dxy_Py_dd+CDY*I_ERI_Fxyz_S_Px_Py_dd;
  Double I_ERI_Fx2z_S_Px_D2y_dd = I_ERI_Fx2z_S_Dxy_Py_dd+CDY*I_ERI_Fx2z_S_Px_Py_dd;
  Double I_ERI_F3y_S_Px_D2y_dd = I_ERI_F3y_S_Dxy_Py_dd+CDY*I_ERI_F3y_S_Px_Py_dd;
  Double I_ERI_F2yz_S_Px_D2y_dd = I_ERI_F2yz_S_Dxy_Py_dd+CDY*I_ERI_F2yz_S_Px_Py_dd;
  Double I_ERI_Fy2z_S_Px_D2y_dd = I_ERI_Fy2z_S_Dxy_Py_dd+CDY*I_ERI_Fy2z_S_Px_Py_dd;
  Double I_ERI_F3z_S_Px_D2y_dd = I_ERI_F3z_S_Dxy_Py_dd+CDY*I_ERI_F3z_S_Px_Py_dd;
  Double I_ERI_F3x_S_Py_D2y_dd = I_ERI_F3x_S_D2y_Py_dd+CDY*I_ERI_F3x_S_Py_Py_dd;
  Double I_ERI_F2xy_S_Py_D2y_dd = I_ERI_F2xy_S_D2y_Py_dd+CDY*I_ERI_F2xy_S_Py_Py_dd;
  Double I_ERI_F2xz_S_Py_D2y_dd = I_ERI_F2xz_S_D2y_Py_dd+CDY*I_ERI_F2xz_S_Py_Py_dd;
  Double I_ERI_Fx2y_S_Py_D2y_dd = I_ERI_Fx2y_S_D2y_Py_dd+CDY*I_ERI_Fx2y_S_Py_Py_dd;
  Double I_ERI_Fxyz_S_Py_D2y_dd = I_ERI_Fxyz_S_D2y_Py_dd+CDY*I_ERI_Fxyz_S_Py_Py_dd;
  Double I_ERI_Fx2z_S_Py_D2y_dd = I_ERI_Fx2z_S_D2y_Py_dd+CDY*I_ERI_Fx2z_S_Py_Py_dd;
  Double I_ERI_F3y_S_Py_D2y_dd = I_ERI_F3y_S_D2y_Py_dd+CDY*I_ERI_F3y_S_Py_Py_dd;
  Double I_ERI_F2yz_S_Py_D2y_dd = I_ERI_F2yz_S_D2y_Py_dd+CDY*I_ERI_F2yz_S_Py_Py_dd;
  Double I_ERI_Fy2z_S_Py_D2y_dd = I_ERI_Fy2z_S_D2y_Py_dd+CDY*I_ERI_Fy2z_S_Py_Py_dd;
  Double I_ERI_F3z_S_Py_D2y_dd = I_ERI_F3z_S_D2y_Py_dd+CDY*I_ERI_F3z_S_Py_Py_dd;
  Double I_ERI_F3x_S_Pz_D2y_dd = I_ERI_F3x_S_Dyz_Py_dd+CDY*I_ERI_F3x_S_Pz_Py_dd;
  Double I_ERI_F2xy_S_Pz_D2y_dd = I_ERI_F2xy_S_Dyz_Py_dd+CDY*I_ERI_F2xy_S_Pz_Py_dd;
  Double I_ERI_F2xz_S_Pz_D2y_dd = I_ERI_F2xz_S_Dyz_Py_dd+CDY*I_ERI_F2xz_S_Pz_Py_dd;
  Double I_ERI_Fx2y_S_Pz_D2y_dd = I_ERI_Fx2y_S_Dyz_Py_dd+CDY*I_ERI_Fx2y_S_Pz_Py_dd;
  Double I_ERI_Fxyz_S_Pz_D2y_dd = I_ERI_Fxyz_S_Dyz_Py_dd+CDY*I_ERI_Fxyz_S_Pz_Py_dd;
  Double I_ERI_Fx2z_S_Pz_D2y_dd = I_ERI_Fx2z_S_Dyz_Py_dd+CDY*I_ERI_Fx2z_S_Pz_Py_dd;
  Double I_ERI_F3y_S_Pz_D2y_dd = I_ERI_F3y_S_Dyz_Py_dd+CDY*I_ERI_F3y_S_Pz_Py_dd;
  Double I_ERI_F2yz_S_Pz_D2y_dd = I_ERI_F2yz_S_Dyz_Py_dd+CDY*I_ERI_F2yz_S_Pz_Py_dd;
  Double I_ERI_Fy2z_S_Pz_D2y_dd = I_ERI_Fy2z_S_Dyz_Py_dd+CDY*I_ERI_Fy2z_S_Pz_Py_dd;
  Double I_ERI_F3z_S_Pz_D2y_dd = I_ERI_F3z_S_Dyz_Py_dd+CDY*I_ERI_F3z_S_Pz_Py_dd;
  Double I_ERI_F3x_S_Px_Dyz_dd = I_ERI_F3x_S_Dxz_Py_dd+CDZ*I_ERI_F3x_S_Px_Py_dd;
  Double I_ERI_F2xy_S_Px_Dyz_dd = I_ERI_F2xy_S_Dxz_Py_dd+CDZ*I_ERI_F2xy_S_Px_Py_dd;
  Double I_ERI_F2xz_S_Px_Dyz_dd = I_ERI_F2xz_S_Dxz_Py_dd+CDZ*I_ERI_F2xz_S_Px_Py_dd;
  Double I_ERI_Fx2y_S_Px_Dyz_dd = I_ERI_Fx2y_S_Dxz_Py_dd+CDZ*I_ERI_Fx2y_S_Px_Py_dd;
  Double I_ERI_Fxyz_S_Px_Dyz_dd = I_ERI_Fxyz_S_Dxz_Py_dd+CDZ*I_ERI_Fxyz_S_Px_Py_dd;
  Double I_ERI_Fx2z_S_Px_Dyz_dd = I_ERI_Fx2z_S_Dxz_Py_dd+CDZ*I_ERI_Fx2z_S_Px_Py_dd;
  Double I_ERI_F3y_S_Px_Dyz_dd = I_ERI_F3y_S_Dxz_Py_dd+CDZ*I_ERI_F3y_S_Px_Py_dd;
  Double I_ERI_F2yz_S_Px_Dyz_dd = I_ERI_F2yz_S_Dxz_Py_dd+CDZ*I_ERI_F2yz_S_Px_Py_dd;
  Double I_ERI_Fy2z_S_Px_Dyz_dd = I_ERI_Fy2z_S_Dxz_Py_dd+CDZ*I_ERI_Fy2z_S_Px_Py_dd;
  Double I_ERI_F3z_S_Px_Dyz_dd = I_ERI_F3z_S_Dxz_Py_dd+CDZ*I_ERI_F3z_S_Px_Py_dd;
  Double I_ERI_F3x_S_Py_Dyz_dd = I_ERI_F3x_S_Dyz_Py_dd+CDZ*I_ERI_F3x_S_Py_Py_dd;
  Double I_ERI_F2xy_S_Py_Dyz_dd = I_ERI_F2xy_S_Dyz_Py_dd+CDZ*I_ERI_F2xy_S_Py_Py_dd;
  Double I_ERI_F2xz_S_Py_Dyz_dd = I_ERI_F2xz_S_Dyz_Py_dd+CDZ*I_ERI_F2xz_S_Py_Py_dd;
  Double I_ERI_Fx2y_S_Py_Dyz_dd = I_ERI_Fx2y_S_Dyz_Py_dd+CDZ*I_ERI_Fx2y_S_Py_Py_dd;
  Double I_ERI_Fxyz_S_Py_Dyz_dd = I_ERI_Fxyz_S_Dyz_Py_dd+CDZ*I_ERI_Fxyz_S_Py_Py_dd;
  Double I_ERI_Fx2z_S_Py_Dyz_dd = I_ERI_Fx2z_S_Dyz_Py_dd+CDZ*I_ERI_Fx2z_S_Py_Py_dd;
  Double I_ERI_F3y_S_Py_Dyz_dd = I_ERI_F3y_S_Dyz_Py_dd+CDZ*I_ERI_F3y_S_Py_Py_dd;
  Double I_ERI_F2yz_S_Py_Dyz_dd = I_ERI_F2yz_S_Dyz_Py_dd+CDZ*I_ERI_F2yz_S_Py_Py_dd;
  Double I_ERI_Fy2z_S_Py_Dyz_dd = I_ERI_Fy2z_S_Dyz_Py_dd+CDZ*I_ERI_Fy2z_S_Py_Py_dd;
  Double I_ERI_F3z_S_Py_Dyz_dd = I_ERI_F3z_S_Dyz_Py_dd+CDZ*I_ERI_F3z_S_Py_Py_dd;
  Double I_ERI_F3x_S_Pz_Dyz_dd = I_ERI_F3x_S_D2z_Py_dd+CDZ*I_ERI_F3x_S_Pz_Py_dd;
  Double I_ERI_F2xy_S_Pz_Dyz_dd = I_ERI_F2xy_S_D2z_Py_dd+CDZ*I_ERI_F2xy_S_Pz_Py_dd;
  Double I_ERI_F2xz_S_Pz_Dyz_dd = I_ERI_F2xz_S_D2z_Py_dd+CDZ*I_ERI_F2xz_S_Pz_Py_dd;
  Double I_ERI_Fx2y_S_Pz_Dyz_dd = I_ERI_Fx2y_S_D2z_Py_dd+CDZ*I_ERI_Fx2y_S_Pz_Py_dd;
  Double I_ERI_Fxyz_S_Pz_Dyz_dd = I_ERI_Fxyz_S_D2z_Py_dd+CDZ*I_ERI_Fxyz_S_Pz_Py_dd;
  Double I_ERI_Fx2z_S_Pz_Dyz_dd = I_ERI_Fx2z_S_D2z_Py_dd+CDZ*I_ERI_Fx2z_S_Pz_Py_dd;
  Double I_ERI_F3y_S_Pz_Dyz_dd = I_ERI_F3y_S_D2z_Py_dd+CDZ*I_ERI_F3y_S_Pz_Py_dd;
  Double I_ERI_F2yz_S_Pz_Dyz_dd = I_ERI_F2yz_S_D2z_Py_dd+CDZ*I_ERI_F2yz_S_Pz_Py_dd;
  Double I_ERI_Fy2z_S_Pz_Dyz_dd = I_ERI_Fy2z_S_D2z_Py_dd+CDZ*I_ERI_Fy2z_S_Pz_Py_dd;
  Double I_ERI_F3z_S_Pz_Dyz_dd = I_ERI_F3z_S_D2z_Py_dd+CDZ*I_ERI_F3z_S_Pz_Py_dd;
  Double I_ERI_F3x_S_Px_D2z_dd = I_ERI_F3x_S_Dxz_Pz_dd+CDZ*I_ERI_F3x_S_Px_Pz_dd;
  Double I_ERI_F2xy_S_Px_D2z_dd = I_ERI_F2xy_S_Dxz_Pz_dd+CDZ*I_ERI_F2xy_S_Px_Pz_dd;
  Double I_ERI_F2xz_S_Px_D2z_dd = I_ERI_F2xz_S_Dxz_Pz_dd+CDZ*I_ERI_F2xz_S_Px_Pz_dd;
  Double I_ERI_Fx2y_S_Px_D2z_dd = I_ERI_Fx2y_S_Dxz_Pz_dd+CDZ*I_ERI_Fx2y_S_Px_Pz_dd;
  Double I_ERI_Fxyz_S_Px_D2z_dd = I_ERI_Fxyz_S_Dxz_Pz_dd+CDZ*I_ERI_Fxyz_S_Px_Pz_dd;
  Double I_ERI_Fx2z_S_Px_D2z_dd = I_ERI_Fx2z_S_Dxz_Pz_dd+CDZ*I_ERI_Fx2z_S_Px_Pz_dd;
  Double I_ERI_F3y_S_Px_D2z_dd = I_ERI_F3y_S_Dxz_Pz_dd+CDZ*I_ERI_F3y_S_Px_Pz_dd;
  Double I_ERI_F2yz_S_Px_D2z_dd = I_ERI_F2yz_S_Dxz_Pz_dd+CDZ*I_ERI_F2yz_S_Px_Pz_dd;
  Double I_ERI_Fy2z_S_Px_D2z_dd = I_ERI_Fy2z_S_Dxz_Pz_dd+CDZ*I_ERI_Fy2z_S_Px_Pz_dd;
  Double I_ERI_F3z_S_Px_D2z_dd = I_ERI_F3z_S_Dxz_Pz_dd+CDZ*I_ERI_F3z_S_Px_Pz_dd;
  Double I_ERI_F3x_S_Py_D2z_dd = I_ERI_F3x_S_Dyz_Pz_dd+CDZ*I_ERI_F3x_S_Py_Pz_dd;
  Double I_ERI_F2xy_S_Py_D2z_dd = I_ERI_F2xy_S_Dyz_Pz_dd+CDZ*I_ERI_F2xy_S_Py_Pz_dd;
  Double I_ERI_F2xz_S_Py_D2z_dd = I_ERI_F2xz_S_Dyz_Pz_dd+CDZ*I_ERI_F2xz_S_Py_Pz_dd;
  Double I_ERI_Fx2y_S_Py_D2z_dd = I_ERI_Fx2y_S_Dyz_Pz_dd+CDZ*I_ERI_Fx2y_S_Py_Pz_dd;
  Double I_ERI_Fxyz_S_Py_D2z_dd = I_ERI_Fxyz_S_Dyz_Pz_dd+CDZ*I_ERI_Fxyz_S_Py_Pz_dd;
  Double I_ERI_Fx2z_S_Py_D2z_dd = I_ERI_Fx2z_S_Dyz_Pz_dd+CDZ*I_ERI_Fx2z_S_Py_Pz_dd;
  Double I_ERI_F3y_S_Py_D2z_dd = I_ERI_F3y_S_Dyz_Pz_dd+CDZ*I_ERI_F3y_S_Py_Pz_dd;
  Double I_ERI_F2yz_S_Py_D2z_dd = I_ERI_F2yz_S_Dyz_Pz_dd+CDZ*I_ERI_F2yz_S_Py_Pz_dd;
  Double I_ERI_Fy2z_S_Py_D2z_dd = I_ERI_Fy2z_S_Dyz_Pz_dd+CDZ*I_ERI_Fy2z_S_Py_Pz_dd;
  Double I_ERI_F3z_S_Py_D2z_dd = I_ERI_F3z_S_Dyz_Pz_dd+CDZ*I_ERI_F3z_S_Py_Pz_dd;
  Double I_ERI_F3x_S_Pz_D2z_dd = I_ERI_F3x_S_D2z_Pz_dd+CDZ*I_ERI_F3x_S_Pz_Pz_dd;
  Double I_ERI_F2xy_S_Pz_D2z_dd = I_ERI_F2xy_S_D2z_Pz_dd+CDZ*I_ERI_F2xy_S_Pz_Pz_dd;
  Double I_ERI_F2xz_S_Pz_D2z_dd = I_ERI_F2xz_S_D2z_Pz_dd+CDZ*I_ERI_F2xz_S_Pz_Pz_dd;
  Double I_ERI_Fx2y_S_Pz_D2z_dd = I_ERI_Fx2y_S_D2z_Pz_dd+CDZ*I_ERI_Fx2y_S_Pz_Pz_dd;
  Double I_ERI_Fxyz_S_Pz_D2z_dd = I_ERI_Fxyz_S_D2z_Pz_dd+CDZ*I_ERI_Fxyz_S_Pz_Pz_dd;
  Double I_ERI_Fx2z_S_Pz_D2z_dd = I_ERI_Fx2z_S_D2z_Pz_dd+CDZ*I_ERI_Fx2z_S_Pz_Pz_dd;
  Double I_ERI_F3y_S_Pz_D2z_dd = I_ERI_F3y_S_D2z_Pz_dd+CDZ*I_ERI_F3y_S_Pz_Pz_dd;
  Double I_ERI_F2yz_S_Pz_D2z_dd = I_ERI_F2yz_S_D2z_Pz_dd+CDZ*I_ERI_F2yz_S_Pz_Pz_dd;
  Double I_ERI_Fy2z_S_Pz_D2z_dd = I_ERI_Fy2z_S_D2z_Pz_dd+CDZ*I_ERI_Fy2z_S_Pz_Pz_dd;
  Double I_ERI_F3z_S_Pz_D2z_dd = I_ERI_F3z_S_D2z_Pz_dd+CDZ*I_ERI_F3z_S_Pz_Pz_dd;

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
   * shell quartet name: SQ_ERI_F_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S_b
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_b = I_ERI_G4x_S_S_S_b+ABX*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Px_S_S_b = I_ERI_G3xy_S_S_S_b+ABX*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Px_S_S_b = I_ERI_G3xz_S_S_S_b+ABX*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Px_S_S_b = I_ERI_G2x2y_S_S_S_b+ABX*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Px_S_S_b = I_ERI_G2xyz_S_S_S_b+ABX*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Px_S_S_b = I_ERI_G2x2z_S_S_S_b+ABX*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Px_S_S_b = I_ERI_Gx3y_S_S_S_b+ABX*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Px_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABX*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Px_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABX*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Px_S_S_b = I_ERI_Gx3z_S_S_S_b+ABX*I_ERI_F3z_S_S_S_b;
  Double I_ERI_F3x_Py_S_S_b = I_ERI_G3xy_S_S_S_b+ABY*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Py_S_S_b = I_ERI_G2x2y_S_S_S_b+ABY*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Py_S_S_b = I_ERI_G2xyz_S_S_S_b+ABY*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Py_S_S_b = I_ERI_Gx3y_S_S_S_b+ABY*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Py_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABY*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Py_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABY*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Py_S_S_b = I_ERI_G4y_S_S_S_b+ABY*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Py_S_S_b = I_ERI_G3yz_S_S_S_b+ABY*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Py_S_S_b = I_ERI_G2y2z_S_S_S_b+ABY*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Py_S_S_b = I_ERI_Gy3z_S_S_S_b+ABY*I_ERI_F3z_S_S_S_b;
  Double I_ERI_F3x_Pz_S_S_b = I_ERI_G3xz_S_S_S_b+ABZ*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Pz_S_S_b = I_ERI_G2xyz_S_S_S_b+ABZ*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Pz_S_S_b = I_ERI_G2x2z_S_S_S_b+ABZ*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Pz_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABZ*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Pz_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABZ*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Pz_S_S_b = I_ERI_Gx3z_S_S_S_b+ABZ*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Pz_S_S_b = I_ERI_G3yz_S_S_S_b+ABZ*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Pz_S_S_b = I_ERI_G2y2z_S_S_S_b+ABZ*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Pz_S_S_b = I_ERI_Gy3z_S_S_S_b+ABZ*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Pz_S_S_b = I_ERI_G4z_S_S_S_b+ABZ*I_ERI_F3z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bb
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_bb = I_ERI_G4x_S_Px_S_bb+ABX*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Px_Px_S_bb = I_ERI_G3xy_S_Px_S_bb+ABX*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Px_Px_S_bb = I_ERI_G3xz_S_Px_S_bb+ABX*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Px_Px_S_bb = I_ERI_G2x2y_S_Px_S_bb+ABX*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Px_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABX*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Px_Px_S_bb = I_ERI_G2x2z_S_Px_S_bb+ABX*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Px_Px_S_bb = I_ERI_Gx3y_S_Px_S_bb+ABX*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Px_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABX*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Px_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABX*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Px_Px_S_bb = I_ERI_Gx3z_S_Px_S_bb+ABX*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Py_Px_S_bb = I_ERI_G3xy_S_Px_S_bb+ABY*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Py_Px_S_bb = I_ERI_G2x2y_S_Px_S_bb+ABY*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Py_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABY*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Py_Px_S_bb = I_ERI_Gx3y_S_Px_S_bb+ABY*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Py_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABY*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Py_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABY*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Py_Px_S_bb = I_ERI_G4y_S_Px_S_bb+ABY*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Py_Px_S_bb = I_ERI_G3yz_S_Px_S_bb+ABY*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Py_Px_S_bb = I_ERI_G2y2z_S_Px_S_bb+ABY*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Py_Px_S_bb = I_ERI_Gy3z_S_Px_S_bb+ABY*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Pz_Px_S_bb = I_ERI_G3xz_S_Px_S_bb+ABZ*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Pz_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABZ*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Pz_Px_S_bb = I_ERI_G2x2z_S_Px_S_bb+ABZ*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Pz_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABZ*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Pz_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABZ*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Pz_Px_S_bb = I_ERI_Gx3z_S_Px_S_bb+ABZ*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Pz_Px_S_bb = I_ERI_G3yz_S_Px_S_bb+ABZ*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Pz_Px_S_bb = I_ERI_G2y2z_S_Px_S_bb+ABZ*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Pz_Px_S_bb = I_ERI_Gy3z_S_Px_S_bb+ABZ*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Pz_Px_S_bb = I_ERI_G4z_S_Px_S_bb+ABZ*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Px_Py_S_bb = I_ERI_G4x_S_Py_S_bb+ABX*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Px_Py_S_bb = I_ERI_G3xy_S_Py_S_bb+ABX*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Px_Py_S_bb = I_ERI_G3xz_S_Py_S_bb+ABX*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Px_Py_S_bb = I_ERI_G2x2y_S_Py_S_bb+ABX*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Px_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABX*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Px_Py_S_bb = I_ERI_G2x2z_S_Py_S_bb+ABX*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Px_Py_S_bb = I_ERI_Gx3y_S_Py_S_bb+ABX*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Px_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABX*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Px_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABX*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Px_Py_S_bb = I_ERI_Gx3z_S_Py_S_bb+ABX*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Py_Py_S_bb = I_ERI_G3xy_S_Py_S_bb+ABY*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Py_Py_S_bb = I_ERI_G2x2y_S_Py_S_bb+ABY*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Py_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABY*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Py_Py_S_bb = I_ERI_Gx3y_S_Py_S_bb+ABY*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Py_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABY*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Py_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABY*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Py_Py_S_bb = I_ERI_G4y_S_Py_S_bb+ABY*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Py_Py_S_bb = I_ERI_G3yz_S_Py_S_bb+ABY*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Py_Py_S_bb = I_ERI_G2y2z_S_Py_S_bb+ABY*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Py_Py_S_bb = I_ERI_Gy3z_S_Py_S_bb+ABY*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Pz_Py_S_bb = I_ERI_G3xz_S_Py_S_bb+ABZ*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Pz_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABZ*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Pz_Py_S_bb = I_ERI_G2x2z_S_Py_S_bb+ABZ*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Pz_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABZ*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Pz_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABZ*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Pz_Py_S_bb = I_ERI_Gx3z_S_Py_S_bb+ABZ*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Pz_Py_S_bb = I_ERI_G3yz_S_Py_S_bb+ABZ*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Pz_Py_S_bb = I_ERI_G2y2z_S_Py_S_bb+ABZ*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Pz_Py_S_bb = I_ERI_Gy3z_S_Py_S_bb+ABZ*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Pz_Py_S_bb = I_ERI_G4z_S_Py_S_bb+ABZ*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Px_Pz_S_bb = I_ERI_G4x_S_Pz_S_bb+ABX*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Px_Pz_S_bb = I_ERI_G3xy_S_Pz_S_bb+ABX*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Px_Pz_S_bb = I_ERI_G3xz_S_Pz_S_bb+ABX*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Px_Pz_S_bb = I_ERI_G2x2y_S_Pz_S_bb+ABX*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Px_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABX*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Px_Pz_S_bb = I_ERI_G2x2z_S_Pz_S_bb+ABX*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Px_Pz_S_bb = I_ERI_Gx3y_S_Pz_S_bb+ABX*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Px_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABX*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Px_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABX*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Px_Pz_S_bb = I_ERI_Gx3z_S_Pz_S_bb+ABX*I_ERI_F3z_S_Pz_S_bb;
  Double I_ERI_F3x_Py_Pz_S_bb = I_ERI_G3xy_S_Pz_S_bb+ABY*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Py_Pz_S_bb = I_ERI_G2x2y_S_Pz_S_bb+ABY*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Py_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABY*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Py_Pz_S_bb = I_ERI_Gx3y_S_Pz_S_bb+ABY*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Py_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABY*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Py_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABY*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Py_Pz_S_bb = I_ERI_G4y_S_Pz_S_bb+ABY*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Py_Pz_S_bb = I_ERI_G3yz_S_Pz_S_bb+ABY*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Py_Pz_S_bb = I_ERI_G2y2z_S_Pz_S_bb+ABY*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Py_Pz_S_bb = I_ERI_Gy3z_S_Pz_S_bb+ABY*I_ERI_F3z_S_Pz_S_bb;
  Double I_ERI_F3x_Pz_Pz_S_bb = I_ERI_G3xz_S_Pz_S_bb+ABZ*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Pz_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABZ*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Pz_Pz_S_bb = I_ERI_G2x2z_S_Pz_S_bb+ABZ*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Pz_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABZ*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Pz_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABZ*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Pz_Pz_S_bb = I_ERI_Gx3z_S_Pz_S_bb+ABZ*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Pz_Pz_S_bb = I_ERI_G3yz_S_Pz_S_bb+ABZ*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Pz_Pz_S_bb = I_ERI_G2y2z_S_Pz_S_bb+ABZ*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Pz_Pz_S_bb = I_ERI_Gy3z_S_Pz_S_bb+ABZ*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Pz_Pz_S_bb = I_ERI_G4z_S_Pz_S_bb+ABZ*I_ERI_F3z_S_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_bb
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bb
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_bb = I_ERI_H5x_S_Px_S_bb+ABX*I_ERI_G4x_S_Px_S_bb;
  Double I_ERI_G3xy_Px_Px_S_bb = I_ERI_H4xy_S_Px_S_bb+ABX*I_ERI_G3xy_S_Px_S_bb;
  Double I_ERI_G3xz_Px_Px_S_bb = I_ERI_H4xz_S_Px_S_bb+ABX*I_ERI_G3xz_S_Px_S_bb;
  Double I_ERI_G2x2y_Px_Px_S_bb = I_ERI_H3x2y_S_Px_S_bb+ABX*I_ERI_G2x2y_S_Px_S_bb;
  Double I_ERI_G2xyz_Px_Px_S_bb = I_ERI_H3xyz_S_Px_S_bb+ABX*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_G2x2z_Px_Px_S_bb = I_ERI_H3x2z_S_Px_S_bb+ABX*I_ERI_G2x2z_S_Px_S_bb;
  Double I_ERI_Gx3y_Px_Px_S_bb = I_ERI_H2x3y_S_Px_S_bb+ABX*I_ERI_Gx3y_S_Px_S_bb;
  Double I_ERI_Gx2yz_Px_Px_S_bb = I_ERI_H2x2yz_S_Px_S_bb+ABX*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Px_Px_S_bb = I_ERI_H2xy2z_S_Px_S_bb+ABX*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_Gx3z_Px_Px_S_bb = I_ERI_H2x3z_S_Px_S_bb+ABX*I_ERI_Gx3z_S_Px_S_bb;
  Double I_ERI_G4y_Px_Px_S_bb = I_ERI_Hx4y_S_Px_S_bb+ABX*I_ERI_G4y_S_Px_S_bb;
  Double I_ERI_G3yz_Px_Px_S_bb = I_ERI_Hx3yz_S_Px_S_bb+ABX*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Px_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABX*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Px_Px_S_bb = I_ERI_Hxy3z_S_Px_S_bb+ABX*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G4z_Px_Px_S_bb = I_ERI_Hx4z_S_Px_S_bb+ABX*I_ERI_G4z_S_Px_S_bb;
  Double I_ERI_G3xy_Py_Px_S_bb = I_ERI_H3x2y_S_Px_S_bb+ABY*I_ERI_G3xy_S_Px_S_bb;
  Double I_ERI_G3xz_Py_Px_S_bb = I_ERI_H3xyz_S_Px_S_bb+ABY*I_ERI_G3xz_S_Px_S_bb;
  Double I_ERI_G2x2y_Py_Px_S_bb = I_ERI_H2x3y_S_Px_S_bb+ABY*I_ERI_G2x2y_S_Px_S_bb;
  Double I_ERI_G2xyz_Py_Px_S_bb = I_ERI_H2x2yz_S_Px_S_bb+ABY*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_G2x2z_Py_Px_S_bb = I_ERI_H2xy2z_S_Px_S_bb+ABY*I_ERI_G2x2z_S_Px_S_bb;
  Double I_ERI_Gx3y_Py_Px_S_bb = I_ERI_Hx4y_S_Px_S_bb+ABY*I_ERI_Gx3y_S_Px_S_bb;
  Double I_ERI_Gx2yz_Py_Px_S_bb = I_ERI_Hx3yz_S_Px_S_bb+ABY*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Py_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABY*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_Gx3z_Py_Px_S_bb = I_ERI_Hxy3z_S_Px_S_bb+ABY*I_ERI_Gx3z_S_Px_S_bb;
  Double I_ERI_G4y_Py_Px_S_bb = I_ERI_H5y_S_Px_S_bb+ABY*I_ERI_G4y_S_Px_S_bb;
  Double I_ERI_G3yz_Py_Px_S_bb = I_ERI_H4yz_S_Px_S_bb+ABY*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Py_Px_S_bb = I_ERI_H3y2z_S_Px_S_bb+ABY*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Py_Px_S_bb = I_ERI_H2y3z_S_Px_S_bb+ABY*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G4z_Py_Px_S_bb = I_ERI_Hy4z_S_Px_S_bb+ABY*I_ERI_G4z_S_Px_S_bb;
  Double I_ERI_G3xz_Pz_Px_S_bb = I_ERI_H3x2z_S_Px_S_bb+ABZ*I_ERI_G3xz_S_Px_S_bb;
  Double I_ERI_G2xyz_Pz_Px_S_bb = I_ERI_H2xy2z_S_Px_S_bb+ABZ*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_G2x2z_Pz_Px_S_bb = I_ERI_H2x3z_S_Px_S_bb+ABZ*I_ERI_G2x2z_S_Px_S_bb;
  Double I_ERI_Gx2yz_Pz_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABZ*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Pz_Px_S_bb = I_ERI_Hxy3z_S_Px_S_bb+ABZ*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_Gx3z_Pz_Px_S_bb = I_ERI_Hx4z_S_Px_S_bb+ABZ*I_ERI_Gx3z_S_Px_S_bb;
  Double I_ERI_G3yz_Pz_Px_S_bb = I_ERI_H3y2z_S_Px_S_bb+ABZ*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Pz_Px_S_bb = I_ERI_H2y3z_S_Px_S_bb+ABZ*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Pz_Px_S_bb = I_ERI_Hy4z_S_Px_S_bb+ABZ*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G4z_Pz_Px_S_bb = I_ERI_H5z_S_Px_S_bb+ABZ*I_ERI_G4z_S_Px_S_bb;
  Double I_ERI_G4x_Px_Py_S_bb = I_ERI_H5x_S_Py_S_bb+ABX*I_ERI_G4x_S_Py_S_bb;
  Double I_ERI_G3xy_Px_Py_S_bb = I_ERI_H4xy_S_Py_S_bb+ABX*I_ERI_G3xy_S_Py_S_bb;
  Double I_ERI_G3xz_Px_Py_S_bb = I_ERI_H4xz_S_Py_S_bb+ABX*I_ERI_G3xz_S_Py_S_bb;
  Double I_ERI_G2x2y_Px_Py_S_bb = I_ERI_H3x2y_S_Py_S_bb+ABX*I_ERI_G2x2y_S_Py_S_bb;
  Double I_ERI_G2xyz_Px_Py_S_bb = I_ERI_H3xyz_S_Py_S_bb+ABX*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_G2x2z_Px_Py_S_bb = I_ERI_H3x2z_S_Py_S_bb+ABX*I_ERI_G2x2z_S_Py_S_bb;
  Double I_ERI_Gx3y_Px_Py_S_bb = I_ERI_H2x3y_S_Py_S_bb+ABX*I_ERI_Gx3y_S_Py_S_bb;
  Double I_ERI_Gx2yz_Px_Py_S_bb = I_ERI_H2x2yz_S_Py_S_bb+ABX*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Px_Py_S_bb = I_ERI_H2xy2z_S_Py_S_bb+ABX*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_Gx3z_Px_Py_S_bb = I_ERI_H2x3z_S_Py_S_bb+ABX*I_ERI_Gx3z_S_Py_S_bb;
  Double I_ERI_G4y_Px_Py_S_bb = I_ERI_Hx4y_S_Py_S_bb+ABX*I_ERI_G4y_S_Py_S_bb;
  Double I_ERI_G3yz_Px_Py_S_bb = I_ERI_Hx3yz_S_Py_S_bb+ABX*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Px_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABX*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Px_Py_S_bb = I_ERI_Hxy3z_S_Py_S_bb+ABX*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G4z_Px_Py_S_bb = I_ERI_Hx4z_S_Py_S_bb+ABX*I_ERI_G4z_S_Py_S_bb;
  Double I_ERI_G3xy_Py_Py_S_bb = I_ERI_H3x2y_S_Py_S_bb+ABY*I_ERI_G3xy_S_Py_S_bb;
  Double I_ERI_G3xz_Py_Py_S_bb = I_ERI_H3xyz_S_Py_S_bb+ABY*I_ERI_G3xz_S_Py_S_bb;
  Double I_ERI_G2x2y_Py_Py_S_bb = I_ERI_H2x3y_S_Py_S_bb+ABY*I_ERI_G2x2y_S_Py_S_bb;
  Double I_ERI_G2xyz_Py_Py_S_bb = I_ERI_H2x2yz_S_Py_S_bb+ABY*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_G2x2z_Py_Py_S_bb = I_ERI_H2xy2z_S_Py_S_bb+ABY*I_ERI_G2x2z_S_Py_S_bb;
  Double I_ERI_Gx3y_Py_Py_S_bb = I_ERI_Hx4y_S_Py_S_bb+ABY*I_ERI_Gx3y_S_Py_S_bb;
  Double I_ERI_Gx2yz_Py_Py_S_bb = I_ERI_Hx3yz_S_Py_S_bb+ABY*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Py_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABY*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_Gx3z_Py_Py_S_bb = I_ERI_Hxy3z_S_Py_S_bb+ABY*I_ERI_Gx3z_S_Py_S_bb;
  Double I_ERI_G4y_Py_Py_S_bb = I_ERI_H5y_S_Py_S_bb+ABY*I_ERI_G4y_S_Py_S_bb;
  Double I_ERI_G3yz_Py_Py_S_bb = I_ERI_H4yz_S_Py_S_bb+ABY*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Py_Py_S_bb = I_ERI_H3y2z_S_Py_S_bb+ABY*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Py_Py_S_bb = I_ERI_H2y3z_S_Py_S_bb+ABY*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G4z_Py_Py_S_bb = I_ERI_Hy4z_S_Py_S_bb+ABY*I_ERI_G4z_S_Py_S_bb;
  Double I_ERI_G3xz_Pz_Py_S_bb = I_ERI_H3x2z_S_Py_S_bb+ABZ*I_ERI_G3xz_S_Py_S_bb;
  Double I_ERI_G2xyz_Pz_Py_S_bb = I_ERI_H2xy2z_S_Py_S_bb+ABZ*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_G2x2z_Pz_Py_S_bb = I_ERI_H2x3z_S_Py_S_bb+ABZ*I_ERI_G2x2z_S_Py_S_bb;
  Double I_ERI_Gx2yz_Pz_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABZ*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Pz_Py_S_bb = I_ERI_Hxy3z_S_Py_S_bb+ABZ*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_Gx3z_Pz_Py_S_bb = I_ERI_Hx4z_S_Py_S_bb+ABZ*I_ERI_Gx3z_S_Py_S_bb;
  Double I_ERI_G3yz_Pz_Py_S_bb = I_ERI_H3y2z_S_Py_S_bb+ABZ*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Pz_Py_S_bb = I_ERI_H2y3z_S_Py_S_bb+ABZ*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Pz_Py_S_bb = I_ERI_Hy4z_S_Py_S_bb+ABZ*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G4z_Pz_Py_S_bb = I_ERI_H5z_S_Py_S_bb+ABZ*I_ERI_G4z_S_Py_S_bb;
  Double I_ERI_G4x_Px_Pz_S_bb = I_ERI_H5x_S_Pz_S_bb+ABX*I_ERI_G4x_S_Pz_S_bb;
  Double I_ERI_G3xy_Px_Pz_S_bb = I_ERI_H4xy_S_Pz_S_bb+ABX*I_ERI_G3xy_S_Pz_S_bb;
  Double I_ERI_G3xz_Px_Pz_S_bb = I_ERI_H4xz_S_Pz_S_bb+ABX*I_ERI_G3xz_S_Pz_S_bb;
  Double I_ERI_G2x2y_Px_Pz_S_bb = I_ERI_H3x2y_S_Pz_S_bb+ABX*I_ERI_G2x2y_S_Pz_S_bb;
  Double I_ERI_G2xyz_Px_Pz_S_bb = I_ERI_H3xyz_S_Pz_S_bb+ABX*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_G2x2z_Px_Pz_S_bb = I_ERI_H3x2z_S_Pz_S_bb+ABX*I_ERI_G2x2z_S_Pz_S_bb;
  Double I_ERI_Gx3y_Px_Pz_S_bb = I_ERI_H2x3y_S_Pz_S_bb+ABX*I_ERI_Gx3y_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Px_Pz_S_bb = I_ERI_H2x2yz_S_Pz_S_bb+ABX*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Px_Pz_S_bb = I_ERI_H2xy2z_S_Pz_S_bb+ABX*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_Gx3z_Px_Pz_S_bb = I_ERI_H2x3z_S_Pz_S_bb+ABX*I_ERI_Gx3z_S_Pz_S_bb;
  Double I_ERI_G4y_Px_Pz_S_bb = I_ERI_Hx4y_S_Pz_S_bb+ABX*I_ERI_G4y_S_Pz_S_bb;
  Double I_ERI_G3yz_Px_Pz_S_bb = I_ERI_Hx3yz_S_Pz_S_bb+ABX*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Px_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABX*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Px_Pz_S_bb = I_ERI_Hxy3z_S_Pz_S_bb+ABX*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G4z_Px_Pz_S_bb = I_ERI_Hx4z_S_Pz_S_bb+ABX*I_ERI_G4z_S_Pz_S_bb;
  Double I_ERI_G3xy_Py_Pz_S_bb = I_ERI_H3x2y_S_Pz_S_bb+ABY*I_ERI_G3xy_S_Pz_S_bb;
  Double I_ERI_G3xz_Py_Pz_S_bb = I_ERI_H3xyz_S_Pz_S_bb+ABY*I_ERI_G3xz_S_Pz_S_bb;
  Double I_ERI_G2x2y_Py_Pz_S_bb = I_ERI_H2x3y_S_Pz_S_bb+ABY*I_ERI_G2x2y_S_Pz_S_bb;
  Double I_ERI_G2xyz_Py_Pz_S_bb = I_ERI_H2x2yz_S_Pz_S_bb+ABY*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_G2x2z_Py_Pz_S_bb = I_ERI_H2xy2z_S_Pz_S_bb+ABY*I_ERI_G2x2z_S_Pz_S_bb;
  Double I_ERI_Gx3y_Py_Pz_S_bb = I_ERI_Hx4y_S_Pz_S_bb+ABY*I_ERI_Gx3y_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Py_Pz_S_bb = I_ERI_Hx3yz_S_Pz_S_bb+ABY*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Py_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABY*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_Gx3z_Py_Pz_S_bb = I_ERI_Hxy3z_S_Pz_S_bb+ABY*I_ERI_Gx3z_S_Pz_S_bb;
  Double I_ERI_G4y_Py_Pz_S_bb = I_ERI_H5y_S_Pz_S_bb+ABY*I_ERI_G4y_S_Pz_S_bb;
  Double I_ERI_G3yz_Py_Pz_S_bb = I_ERI_H4yz_S_Pz_S_bb+ABY*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Py_Pz_S_bb = I_ERI_H3y2z_S_Pz_S_bb+ABY*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Py_Pz_S_bb = I_ERI_H2y3z_S_Pz_S_bb+ABY*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G4z_Py_Pz_S_bb = I_ERI_Hy4z_S_Pz_S_bb+ABY*I_ERI_G4z_S_Pz_S_bb;
  Double I_ERI_G3xz_Pz_Pz_S_bb = I_ERI_H3x2z_S_Pz_S_bb+ABZ*I_ERI_G3xz_S_Pz_S_bb;
  Double I_ERI_G2xyz_Pz_Pz_S_bb = I_ERI_H2xy2z_S_Pz_S_bb+ABZ*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_G2x2z_Pz_Pz_S_bb = I_ERI_H2x3z_S_Pz_S_bb+ABZ*I_ERI_G2x2z_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Pz_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABZ*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Pz_Pz_S_bb = I_ERI_Hxy3z_S_Pz_S_bb+ABZ*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_Gx3z_Pz_Pz_S_bb = I_ERI_Hx4z_S_Pz_S_bb+ABZ*I_ERI_Gx3z_S_Pz_S_bb;
  Double I_ERI_G3yz_Pz_Pz_S_bb = I_ERI_H3y2z_S_Pz_S_bb+ABZ*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Pz_Pz_S_bb = I_ERI_H2y3z_S_Pz_S_bb+ABZ*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Pz_Pz_S_bb = I_ERI_Hy4z_S_Pz_S_bb+ABZ*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G4z_Pz_Pz_S_bb = I_ERI_H5z_S_Pz_S_bb+ABZ*I_ERI_G4z_S_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_P_P_S_bb
   ************************************************************/
  Double I_ERI_F3x_D2x_Px_S_bb = I_ERI_G4x_Px_Px_S_bb+ABX*I_ERI_F3x_Px_Px_S_bb;
  Double I_ERI_F2xy_D2x_Px_S_bb = I_ERI_G3xy_Px_Px_S_bb+ABX*I_ERI_F2xy_Px_Px_S_bb;
  Double I_ERI_F2xz_D2x_Px_S_bb = I_ERI_G3xz_Px_Px_S_bb+ABX*I_ERI_F2xz_Px_Px_S_bb;
  Double I_ERI_Fx2y_D2x_Px_S_bb = I_ERI_G2x2y_Px_Px_S_bb+ABX*I_ERI_Fx2y_Px_Px_S_bb;
  Double I_ERI_Fxyz_D2x_Px_S_bb = I_ERI_G2xyz_Px_Px_S_bb+ABX*I_ERI_Fxyz_Px_Px_S_bb;
  Double I_ERI_Fx2z_D2x_Px_S_bb = I_ERI_G2x2z_Px_Px_S_bb+ABX*I_ERI_Fx2z_Px_Px_S_bb;
  Double I_ERI_F3y_D2x_Px_S_bb = I_ERI_Gx3y_Px_Px_S_bb+ABX*I_ERI_F3y_Px_Px_S_bb;
  Double I_ERI_F2yz_D2x_Px_S_bb = I_ERI_Gx2yz_Px_Px_S_bb+ABX*I_ERI_F2yz_Px_Px_S_bb;
  Double I_ERI_Fy2z_D2x_Px_S_bb = I_ERI_Gxy2z_Px_Px_S_bb+ABX*I_ERI_Fy2z_Px_Px_S_bb;
  Double I_ERI_F3z_D2x_Px_S_bb = I_ERI_Gx3z_Px_Px_S_bb+ABX*I_ERI_F3z_Px_Px_S_bb;
  Double I_ERI_F3x_Dxy_Px_S_bb = I_ERI_G3xy_Px_Px_S_bb+ABY*I_ERI_F3x_Px_Px_S_bb;
  Double I_ERI_F2xy_Dxy_Px_S_bb = I_ERI_G2x2y_Px_Px_S_bb+ABY*I_ERI_F2xy_Px_Px_S_bb;
  Double I_ERI_F2xz_Dxy_Px_S_bb = I_ERI_G2xyz_Px_Px_S_bb+ABY*I_ERI_F2xz_Px_Px_S_bb;
  Double I_ERI_Fx2y_Dxy_Px_S_bb = I_ERI_Gx3y_Px_Px_S_bb+ABY*I_ERI_Fx2y_Px_Px_S_bb;
  Double I_ERI_Fxyz_Dxy_Px_S_bb = I_ERI_Gx2yz_Px_Px_S_bb+ABY*I_ERI_Fxyz_Px_Px_S_bb;
  Double I_ERI_Fx2z_Dxy_Px_S_bb = I_ERI_Gxy2z_Px_Px_S_bb+ABY*I_ERI_Fx2z_Px_Px_S_bb;
  Double I_ERI_F3y_Dxy_Px_S_bb = I_ERI_G4y_Px_Px_S_bb+ABY*I_ERI_F3y_Px_Px_S_bb;
  Double I_ERI_F2yz_Dxy_Px_S_bb = I_ERI_G3yz_Px_Px_S_bb+ABY*I_ERI_F2yz_Px_Px_S_bb;
  Double I_ERI_Fy2z_Dxy_Px_S_bb = I_ERI_G2y2z_Px_Px_S_bb+ABY*I_ERI_Fy2z_Px_Px_S_bb;
  Double I_ERI_F3z_Dxy_Px_S_bb = I_ERI_Gy3z_Px_Px_S_bb+ABY*I_ERI_F3z_Px_Px_S_bb;
  Double I_ERI_F3x_Dxz_Px_S_bb = I_ERI_G3xz_Px_Px_S_bb+ABZ*I_ERI_F3x_Px_Px_S_bb;
  Double I_ERI_F2xy_Dxz_Px_S_bb = I_ERI_G2xyz_Px_Px_S_bb+ABZ*I_ERI_F2xy_Px_Px_S_bb;
  Double I_ERI_F2xz_Dxz_Px_S_bb = I_ERI_G2x2z_Px_Px_S_bb+ABZ*I_ERI_F2xz_Px_Px_S_bb;
  Double I_ERI_Fx2y_Dxz_Px_S_bb = I_ERI_Gx2yz_Px_Px_S_bb+ABZ*I_ERI_Fx2y_Px_Px_S_bb;
  Double I_ERI_Fxyz_Dxz_Px_S_bb = I_ERI_Gxy2z_Px_Px_S_bb+ABZ*I_ERI_Fxyz_Px_Px_S_bb;
  Double I_ERI_Fx2z_Dxz_Px_S_bb = I_ERI_Gx3z_Px_Px_S_bb+ABZ*I_ERI_Fx2z_Px_Px_S_bb;
  Double I_ERI_F3y_Dxz_Px_S_bb = I_ERI_G3yz_Px_Px_S_bb+ABZ*I_ERI_F3y_Px_Px_S_bb;
  Double I_ERI_F2yz_Dxz_Px_S_bb = I_ERI_G2y2z_Px_Px_S_bb+ABZ*I_ERI_F2yz_Px_Px_S_bb;
  Double I_ERI_Fy2z_Dxz_Px_S_bb = I_ERI_Gy3z_Px_Px_S_bb+ABZ*I_ERI_Fy2z_Px_Px_S_bb;
  Double I_ERI_F3z_Dxz_Px_S_bb = I_ERI_G4z_Px_Px_S_bb+ABZ*I_ERI_F3z_Px_Px_S_bb;
  Double I_ERI_F3x_D2y_Px_S_bb = I_ERI_G3xy_Py_Px_S_bb+ABY*I_ERI_F3x_Py_Px_S_bb;
  Double I_ERI_F2xy_D2y_Px_S_bb = I_ERI_G2x2y_Py_Px_S_bb+ABY*I_ERI_F2xy_Py_Px_S_bb;
  Double I_ERI_F2xz_D2y_Px_S_bb = I_ERI_G2xyz_Py_Px_S_bb+ABY*I_ERI_F2xz_Py_Px_S_bb;
  Double I_ERI_Fx2y_D2y_Px_S_bb = I_ERI_Gx3y_Py_Px_S_bb+ABY*I_ERI_Fx2y_Py_Px_S_bb;
  Double I_ERI_Fxyz_D2y_Px_S_bb = I_ERI_Gx2yz_Py_Px_S_bb+ABY*I_ERI_Fxyz_Py_Px_S_bb;
  Double I_ERI_Fx2z_D2y_Px_S_bb = I_ERI_Gxy2z_Py_Px_S_bb+ABY*I_ERI_Fx2z_Py_Px_S_bb;
  Double I_ERI_F3y_D2y_Px_S_bb = I_ERI_G4y_Py_Px_S_bb+ABY*I_ERI_F3y_Py_Px_S_bb;
  Double I_ERI_F2yz_D2y_Px_S_bb = I_ERI_G3yz_Py_Px_S_bb+ABY*I_ERI_F2yz_Py_Px_S_bb;
  Double I_ERI_Fy2z_D2y_Px_S_bb = I_ERI_G2y2z_Py_Px_S_bb+ABY*I_ERI_Fy2z_Py_Px_S_bb;
  Double I_ERI_F3z_D2y_Px_S_bb = I_ERI_Gy3z_Py_Px_S_bb+ABY*I_ERI_F3z_Py_Px_S_bb;
  Double I_ERI_F3x_Dyz_Px_S_bb = I_ERI_G3xz_Py_Px_S_bb+ABZ*I_ERI_F3x_Py_Px_S_bb;
  Double I_ERI_F2xy_Dyz_Px_S_bb = I_ERI_G2xyz_Py_Px_S_bb+ABZ*I_ERI_F2xy_Py_Px_S_bb;
  Double I_ERI_F2xz_Dyz_Px_S_bb = I_ERI_G2x2z_Py_Px_S_bb+ABZ*I_ERI_F2xz_Py_Px_S_bb;
  Double I_ERI_Fx2y_Dyz_Px_S_bb = I_ERI_Gx2yz_Py_Px_S_bb+ABZ*I_ERI_Fx2y_Py_Px_S_bb;
  Double I_ERI_Fxyz_Dyz_Px_S_bb = I_ERI_Gxy2z_Py_Px_S_bb+ABZ*I_ERI_Fxyz_Py_Px_S_bb;
  Double I_ERI_Fx2z_Dyz_Px_S_bb = I_ERI_Gx3z_Py_Px_S_bb+ABZ*I_ERI_Fx2z_Py_Px_S_bb;
  Double I_ERI_F3y_Dyz_Px_S_bb = I_ERI_G3yz_Py_Px_S_bb+ABZ*I_ERI_F3y_Py_Px_S_bb;
  Double I_ERI_F2yz_Dyz_Px_S_bb = I_ERI_G2y2z_Py_Px_S_bb+ABZ*I_ERI_F2yz_Py_Px_S_bb;
  Double I_ERI_Fy2z_Dyz_Px_S_bb = I_ERI_Gy3z_Py_Px_S_bb+ABZ*I_ERI_Fy2z_Py_Px_S_bb;
  Double I_ERI_F3z_Dyz_Px_S_bb = I_ERI_G4z_Py_Px_S_bb+ABZ*I_ERI_F3z_Py_Px_S_bb;
  Double I_ERI_F3x_D2z_Px_S_bb = I_ERI_G3xz_Pz_Px_S_bb+ABZ*I_ERI_F3x_Pz_Px_S_bb;
  Double I_ERI_F2xy_D2z_Px_S_bb = I_ERI_G2xyz_Pz_Px_S_bb+ABZ*I_ERI_F2xy_Pz_Px_S_bb;
  Double I_ERI_F2xz_D2z_Px_S_bb = I_ERI_G2x2z_Pz_Px_S_bb+ABZ*I_ERI_F2xz_Pz_Px_S_bb;
  Double I_ERI_Fx2y_D2z_Px_S_bb = I_ERI_Gx2yz_Pz_Px_S_bb+ABZ*I_ERI_Fx2y_Pz_Px_S_bb;
  Double I_ERI_Fxyz_D2z_Px_S_bb = I_ERI_Gxy2z_Pz_Px_S_bb+ABZ*I_ERI_Fxyz_Pz_Px_S_bb;
  Double I_ERI_Fx2z_D2z_Px_S_bb = I_ERI_Gx3z_Pz_Px_S_bb+ABZ*I_ERI_Fx2z_Pz_Px_S_bb;
  Double I_ERI_F3y_D2z_Px_S_bb = I_ERI_G3yz_Pz_Px_S_bb+ABZ*I_ERI_F3y_Pz_Px_S_bb;
  Double I_ERI_F2yz_D2z_Px_S_bb = I_ERI_G2y2z_Pz_Px_S_bb+ABZ*I_ERI_F2yz_Pz_Px_S_bb;
  Double I_ERI_Fy2z_D2z_Px_S_bb = I_ERI_Gy3z_Pz_Px_S_bb+ABZ*I_ERI_Fy2z_Pz_Px_S_bb;
  Double I_ERI_F3z_D2z_Px_S_bb = I_ERI_G4z_Pz_Px_S_bb+ABZ*I_ERI_F3z_Pz_Px_S_bb;
  Double I_ERI_F3x_D2x_Py_S_bb = I_ERI_G4x_Px_Py_S_bb+ABX*I_ERI_F3x_Px_Py_S_bb;
  Double I_ERI_F2xy_D2x_Py_S_bb = I_ERI_G3xy_Px_Py_S_bb+ABX*I_ERI_F2xy_Px_Py_S_bb;
  Double I_ERI_F2xz_D2x_Py_S_bb = I_ERI_G3xz_Px_Py_S_bb+ABX*I_ERI_F2xz_Px_Py_S_bb;
  Double I_ERI_Fx2y_D2x_Py_S_bb = I_ERI_G2x2y_Px_Py_S_bb+ABX*I_ERI_Fx2y_Px_Py_S_bb;
  Double I_ERI_Fxyz_D2x_Py_S_bb = I_ERI_G2xyz_Px_Py_S_bb+ABX*I_ERI_Fxyz_Px_Py_S_bb;
  Double I_ERI_Fx2z_D2x_Py_S_bb = I_ERI_G2x2z_Px_Py_S_bb+ABX*I_ERI_Fx2z_Px_Py_S_bb;
  Double I_ERI_F3y_D2x_Py_S_bb = I_ERI_Gx3y_Px_Py_S_bb+ABX*I_ERI_F3y_Px_Py_S_bb;
  Double I_ERI_F2yz_D2x_Py_S_bb = I_ERI_Gx2yz_Px_Py_S_bb+ABX*I_ERI_F2yz_Px_Py_S_bb;
  Double I_ERI_Fy2z_D2x_Py_S_bb = I_ERI_Gxy2z_Px_Py_S_bb+ABX*I_ERI_Fy2z_Px_Py_S_bb;
  Double I_ERI_F3z_D2x_Py_S_bb = I_ERI_Gx3z_Px_Py_S_bb+ABX*I_ERI_F3z_Px_Py_S_bb;
  Double I_ERI_F3x_Dxy_Py_S_bb = I_ERI_G3xy_Px_Py_S_bb+ABY*I_ERI_F3x_Px_Py_S_bb;
  Double I_ERI_F2xy_Dxy_Py_S_bb = I_ERI_G2x2y_Px_Py_S_bb+ABY*I_ERI_F2xy_Px_Py_S_bb;
  Double I_ERI_F2xz_Dxy_Py_S_bb = I_ERI_G2xyz_Px_Py_S_bb+ABY*I_ERI_F2xz_Px_Py_S_bb;
  Double I_ERI_Fx2y_Dxy_Py_S_bb = I_ERI_Gx3y_Px_Py_S_bb+ABY*I_ERI_Fx2y_Px_Py_S_bb;
  Double I_ERI_Fxyz_Dxy_Py_S_bb = I_ERI_Gx2yz_Px_Py_S_bb+ABY*I_ERI_Fxyz_Px_Py_S_bb;
  Double I_ERI_Fx2z_Dxy_Py_S_bb = I_ERI_Gxy2z_Px_Py_S_bb+ABY*I_ERI_Fx2z_Px_Py_S_bb;
  Double I_ERI_F3y_Dxy_Py_S_bb = I_ERI_G4y_Px_Py_S_bb+ABY*I_ERI_F3y_Px_Py_S_bb;
  Double I_ERI_F2yz_Dxy_Py_S_bb = I_ERI_G3yz_Px_Py_S_bb+ABY*I_ERI_F2yz_Px_Py_S_bb;
  Double I_ERI_Fy2z_Dxy_Py_S_bb = I_ERI_G2y2z_Px_Py_S_bb+ABY*I_ERI_Fy2z_Px_Py_S_bb;
  Double I_ERI_F3z_Dxy_Py_S_bb = I_ERI_Gy3z_Px_Py_S_bb+ABY*I_ERI_F3z_Px_Py_S_bb;
  Double I_ERI_F3x_Dxz_Py_S_bb = I_ERI_G3xz_Px_Py_S_bb+ABZ*I_ERI_F3x_Px_Py_S_bb;
  Double I_ERI_F2xy_Dxz_Py_S_bb = I_ERI_G2xyz_Px_Py_S_bb+ABZ*I_ERI_F2xy_Px_Py_S_bb;
  Double I_ERI_F2xz_Dxz_Py_S_bb = I_ERI_G2x2z_Px_Py_S_bb+ABZ*I_ERI_F2xz_Px_Py_S_bb;
  Double I_ERI_Fx2y_Dxz_Py_S_bb = I_ERI_Gx2yz_Px_Py_S_bb+ABZ*I_ERI_Fx2y_Px_Py_S_bb;
  Double I_ERI_Fxyz_Dxz_Py_S_bb = I_ERI_Gxy2z_Px_Py_S_bb+ABZ*I_ERI_Fxyz_Px_Py_S_bb;
  Double I_ERI_Fx2z_Dxz_Py_S_bb = I_ERI_Gx3z_Px_Py_S_bb+ABZ*I_ERI_Fx2z_Px_Py_S_bb;
  Double I_ERI_F3y_Dxz_Py_S_bb = I_ERI_G3yz_Px_Py_S_bb+ABZ*I_ERI_F3y_Px_Py_S_bb;
  Double I_ERI_F2yz_Dxz_Py_S_bb = I_ERI_G2y2z_Px_Py_S_bb+ABZ*I_ERI_F2yz_Px_Py_S_bb;
  Double I_ERI_Fy2z_Dxz_Py_S_bb = I_ERI_Gy3z_Px_Py_S_bb+ABZ*I_ERI_Fy2z_Px_Py_S_bb;
  Double I_ERI_F3z_Dxz_Py_S_bb = I_ERI_G4z_Px_Py_S_bb+ABZ*I_ERI_F3z_Px_Py_S_bb;
  Double I_ERI_F3x_D2y_Py_S_bb = I_ERI_G3xy_Py_Py_S_bb+ABY*I_ERI_F3x_Py_Py_S_bb;
  Double I_ERI_F2xy_D2y_Py_S_bb = I_ERI_G2x2y_Py_Py_S_bb+ABY*I_ERI_F2xy_Py_Py_S_bb;
  Double I_ERI_F2xz_D2y_Py_S_bb = I_ERI_G2xyz_Py_Py_S_bb+ABY*I_ERI_F2xz_Py_Py_S_bb;
  Double I_ERI_Fx2y_D2y_Py_S_bb = I_ERI_Gx3y_Py_Py_S_bb+ABY*I_ERI_Fx2y_Py_Py_S_bb;
  Double I_ERI_Fxyz_D2y_Py_S_bb = I_ERI_Gx2yz_Py_Py_S_bb+ABY*I_ERI_Fxyz_Py_Py_S_bb;
  Double I_ERI_Fx2z_D2y_Py_S_bb = I_ERI_Gxy2z_Py_Py_S_bb+ABY*I_ERI_Fx2z_Py_Py_S_bb;
  Double I_ERI_F3y_D2y_Py_S_bb = I_ERI_G4y_Py_Py_S_bb+ABY*I_ERI_F3y_Py_Py_S_bb;
  Double I_ERI_F2yz_D2y_Py_S_bb = I_ERI_G3yz_Py_Py_S_bb+ABY*I_ERI_F2yz_Py_Py_S_bb;
  Double I_ERI_Fy2z_D2y_Py_S_bb = I_ERI_G2y2z_Py_Py_S_bb+ABY*I_ERI_Fy2z_Py_Py_S_bb;
  Double I_ERI_F3z_D2y_Py_S_bb = I_ERI_Gy3z_Py_Py_S_bb+ABY*I_ERI_F3z_Py_Py_S_bb;
  Double I_ERI_F3x_Dyz_Py_S_bb = I_ERI_G3xz_Py_Py_S_bb+ABZ*I_ERI_F3x_Py_Py_S_bb;
  Double I_ERI_F2xy_Dyz_Py_S_bb = I_ERI_G2xyz_Py_Py_S_bb+ABZ*I_ERI_F2xy_Py_Py_S_bb;
  Double I_ERI_F2xz_Dyz_Py_S_bb = I_ERI_G2x2z_Py_Py_S_bb+ABZ*I_ERI_F2xz_Py_Py_S_bb;
  Double I_ERI_Fx2y_Dyz_Py_S_bb = I_ERI_Gx2yz_Py_Py_S_bb+ABZ*I_ERI_Fx2y_Py_Py_S_bb;
  Double I_ERI_Fxyz_Dyz_Py_S_bb = I_ERI_Gxy2z_Py_Py_S_bb+ABZ*I_ERI_Fxyz_Py_Py_S_bb;
  Double I_ERI_Fx2z_Dyz_Py_S_bb = I_ERI_Gx3z_Py_Py_S_bb+ABZ*I_ERI_Fx2z_Py_Py_S_bb;
  Double I_ERI_F3y_Dyz_Py_S_bb = I_ERI_G3yz_Py_Py_S_bb+ABZ*I_ERI_F3y_Py_Py_S_bb;
  Double I_ERI_F2yz_Dyz_Py_S_bb = I_ERI_G2y2z_Py_Py_S_bb+ABZ*I_ERI_F2yz_Py_Py_S_bb;
  Double I_ERI_Fy2z_Dyz_Py_S_bb = I_ERI_Gy3z_Py_Py_S_bb+ABZ*I_ERI_Fy2z_Py_Py_S_bb;
  Double I_ERI_F3z_Dyz_Py_S_bb = I_ERI_G4z_Py_Py_S_bb+ABZ*I_ERI_F3z_Py_Py_S_bb;
  Double I_ERI_F3x_D2z_Py_S_bb = I_ERI_G3xz_Pz_Py_S_bb+ABZ*I_ERI_F3x_Pz_Py_S_bb;
  Double I_ERI_F2xy_D2z_Py_S_bb = I_ERI_G2xyz_Pz_Py_S_bb+ABZ*I_ERI_F2xy_Pz_Py_S_bb;
  Double I_ERI_F2xz_D2z_Py_S_bb = I_ERI_G2x2z_Pz_Py_S_bb+ABZ*I_ERI_F2xz_Pz_Py_S_bb;
  Double I_ERI_Fx2y_D2z_Py_S_bb = I_ERI_Gx2yz_Pz_Py_S_bb+ABZ*I_ERI_Fx2y_Pz_Py_S_bb;
  Double I_ERI_Fxyz_D2z_Py_S_bb = I_ERI_Gxy2z_Pz_Py_S_bb+ABZ*I_ERI_Fxyz_Pz_Py_S_bb;
  Double I_ERI_Fx2z_D2z_Py_S_bb = I_ERI_Gx3z_Pz_Py_S_bb+ABZ*I_ERI_Fx2z_Pz_Py_S_bb;
  Double I_ERI_F3y_D2z_Py_S_bb = I_ERI_G3yz_Pz_Py_S_bb+ABZ*I_ERI_F3y_Pz_Py_S_bb;
  Double I_ERI_F2yz_D2z_Py_S_bb = I_ERI_G2y2z_Pz_Py_S_bb+ABZ*I_ERI_F2yz_Pz_Py_S_bb;
  Double I_ERI_Fy2z_D2z_Py_S_bb = I_ERI_Gy3z_Pz_Py_S_bb+ABZ*I_ERI_Fy2z_Pz_Py_S_bb;
  Double I_ERI_F3z_D2z_Py_S_bb = I_ERI_G4z_Pz_Py_S_bb+ABZ*I_ERI_F3z_Pz_Py_S_bb;
  Double I_ERI_F3x_D2x_Pz_S_bb = I_ERI_G4x_Px_Pz_S_bb+ABX*I_ERI_F3x_Px_Pz_S_bb;
  Double I_ERI_F2xy_D2x_Pz_S_bb = I_ERI_G3xy_Px_Pz_S_bb+ABX*I_ERI_F2xy_Px_Pz_S_bb;
  Double I_ERI_F2xz_D2x_Pz_S_bb = I_ERI_G3xz_Px_Pz_S_bb+ABX*I_ERI_F2xz_Px_Pz_S_bb;
  Double I_ERI_Fx2y_D2x_Pz_S_bb = I_ERI_G2x2y_Px_Pz_S_bb+ABX*I_ERI_Fx2y_Px_Pz_S_bb;
  Double I_ERI_Fxyz_D2x_Pz_S_bb = I_ERI_G2xyz_Px_Pz_S_bb+ABX*I_ERI_Fxyz_Px_Pz_S_bb;
  Double I_ERI_Fx2z_D2x_Pz_S_bb = I_ERI_G2x2z_Px_Pz_S_bb+ABX*I_ERI_Fx2z_Px_Pz_S_bb;
  Double I_ERI_F3y_D2x_Pz_S_bb = I_ERI_Gx3y_Px_Pz_S_bb+ABX*I_ERI_F3y_Px_Pz_S_bb;
  Double I_ERI_F2yz_D2x_Pz_S_bb = I_ERI_Gx2yz_Px_Pz_S_bb+ABX*I_ERI_F2yz_Px_Pz_S_bb;
  Double I_ERI_Fy2z_D2x_Pz_S_bb = I_ERI_Gxy2z_Px_Pz_S_bb+ABX*I_ERI_Fy2z_Px_Pz_S_bb;
  Double I_ERI_F3z_D2x_Pz_S_bb = I_ERI_Gx3z_Px_Pz_S_bb+ABX*I_ERI_F3z_Px_Pz_S_bb;
  Double I_ERI_F3x_Dxy_Pz_S_bb = I_ERI_G3xy_Px_Pz_S_bb+ABY*I_ERI_F3x_Px_Pz_S_bb;
  Double I_ERI_F2xy_Dxy_Pz_S_bb = I_ERI_G2x2y_Px_Pz_S_bb+ABY*I_ERI_F2xy_Px_Pz_S_bb;
  Double I_ERI_F2xz_Dxy_Pz_S_bb = I_ERI_G2xyz_Px_Pz_S_bb+ABY*I_ERI_F2xz_Px_Pz_S_bb;
  Double I_ERI_Fx2y_Dxy_Pz_S_bb = I_ERI_Gx3y_Px_Pz_S_bb+ABY*I_ERI_Fx2y_Px_Pz_S_bb;
  Double I_ERI_Fxyz_Dxy_Pz_S_bb = I_ERI_Gx2yz_Px_Pz_S_bb+ABY*I_ERI_Fxyz_Px_Pz_S_bb;
  Double I_ERI_Fx2z_Dxy_Pz_S_bb = I_ERI_Gxy2z_Px_Pz_S_bb+ABY*I_ERI_Fx2z_Px_Pz_S_bb;
  Double I_ERI_F3y_Dxy_Pz_S_bb = I_ERI_G4y_Px_Pz_S_bb+ABY*I_ERI_F3y_Px_Pz_S_bb;
  Double I_ERI_F2yz_Dxy_Pz_S_bb = I_ERI_G3yz_Px_Pz_S_bb+ABY*I_ERI_F2yz_Px_Pz_S_bb;
  Double I_ERI_Fy2z_Dxy_Pz_S_bb = I_ERI_G2y2z_Px_Pz_S_bb+ABY*I_ERI_Fy2z_Px_Pz_S_bb;
  Double I_ERI_F3z_Dxy_Pz_S_bb = I_ERI_Gy3z_Px_Pz_S_bb+ABY*I_ERI_F3z_Px_Pz_S_bb;
  Double I_ERI_F3x_Dxz_Pz_S_bb = I_ERI_G3xz_Px_Pz_S_bb+ABZ*I_ERI_F3x_Px_Pz_S_bb;
  Double I_ERI_F2xy_Dxz_Pz_S_bb = I_ERI_G2xyz_Px_Pz_S_bb+ABZ*I_ERI_F2xy_Px_Pz_S_bb;
  Double I_ERI_F2xz_Dxz_Pz_S_bb = I_ERI_G2x2z_Px_Pz_S_bb+ABZ*I_ERI_F2xz_Px_Pz_S_bb;
  Double I_ERI_Fx2y_Dxz_Pz_S_bb = I_ERI_Gx2yz_Px_Pz_S_bb+ABZ*I_ERI_Fx2y_Px_Pz_S_bb;
  Double I_ERI_Fxyz_Dxz_Pz_S_bb = I_ERI_Gxy2z_Px_Pz_S_bb+ABZ*I_ERI_Fxyz_Px_Pz_S_bb;
  Double I_ERI_Fx2z_Dxz_Pz_S_bb = I_ERI_Gx3z_Px_Pz_S_bb+ABZ*I_ERI_Fx2z_Px_Pz_S_bb;
  Double I_ERI_F3y_Dxz_Pz_S_bb = I_ERI_G3yz_Px_Pz_S_bb+ABZ*I_ERI_F3y_Px_Pz_S_bb;
  Double I_ERI_F2yz_Dxz_Pz_S_bb = I_ERI_G2y2z_Px_Pz_S_bb+ABZ*I_ERI_F2yz_Px_Pz_S_bb;
  Double I_ERI_Fy2z_Dxz_Pz_S_bb = I_ERI_Gy3z_Px_Pz_S_bb+ABZ*I_ERI_Fy2z_Px_Pz_S_bb;
  Double I_ERI_F3z_Dxz_Pz_S_bb = I_ERI_G4z_Px_Pz_S_bb+ABZ*I_ERI_F3z_Px_Pz_S_bb;
  Double I_ERI_F3x_D2y_Pz_S_bb = I_ERI_G3xy_Py_Pz_S_bb+ABY*I_ERI_F3x_Py_Pz_S_bb;
  Double I_ERI_F2xy_D2y_Pz_S_bb = I_ERI_G2x2y_Py_Pz_S_bb+ABY*I_ERI_F2xy_Py_Pz_S_bb;
  Double I_ERI_F2xz_D2y_Pz_S_bb = I_ERI_G2xyz_Py_Pz_S_bb+ABY*I_ERI_F2xz_Py_Pz_S_bb;
  Double I_ERI_Fx2y_D2y_Pz_S_bb = I_ERI_Gx3y_Py_Pz_S_bb+ABY*I_ERI_Fx2y_Py_Pz_S_bb;
  Double I_ERI_Fxyz_D2y_Pz_S_bb = I_ERI_Gx2yz_Py_Pz_S_bb+ABY*I_ERI_Fxyz_Py_Pz_S_bb;
  Double I_ERI_Fx2z_D2y_Pz_S_bb = I_ERI_Gxy2z_Py_Pz_S_bb+ABY*I_ERI_Fx2z_Py_Pz_S_bb;
  Double I_ERI_F3y_D2y_Pz_S_bb = I_ERI_G4y_Py_Pz_S_bb+ABY*I_ERI_F3y_Py_Pz_S_bb;
  Double I_ERI_F2yz_D2y_Pz_S_bb = I_ERI_G3yz_Py_Pz_S_bb+ABY*I_ERI_F2yz_Py_Pz_S_bb;
  Double I_ERI_Fy2z_D2y_Pz_S_bb = I_ERI_G2y2z_Py_Pz_S_bb+ABY*I_ERI_Fy2z_Py_Pz_S_bb;
  Double I_ERI_F3z_D2y_Pz_S_bb = I_ERI_Gy3z_Py_Pz_S_bb+ABY*I_ERI_F3z_Py_Pz_S_bb;
  Double I_ERI_F3x_Dyz_Pz_S_bb = I_ERI_G3xz_Py_Pz_S_bb+ABZ*I_ERI_F3x_Py_Pz_S_bb;
  Double I_ERI_F2xy_Dyz_Pz_S_bb = I_ERI_G2xyz_Py_Pz_S_bb+ABZ*I_ERI_F2xy_Py_Pz_S_bb;
  Double I_ERI_F2xz_Dyz_Pz_S_bb = I_ERI_G2x2z_Py_Pz_S_bb+ABZ*I_ERI_F2xz_Py_Pz_S_bb;
  Double I_ERI_Fx2y_Dyz_Pz_S_bb = I_ERI_Gx2yz_Py_Pz_S_bb+ABZ*I_ERI_Fx2y_Py_Pz_S_bb;
  Double I_ERI_Fxyz_Dyz_Pz_S_bb = I_ERI_Gxy2z_Py_Pz_S_bb+ABZ*I_ERI_Fxyz_Py_Pz_S_bb;
  Double I_ERI_Fx2z_Dyz_Pz_S_bb = I_ERI_Gx3z_Py_Pz_S_bb+ABZ*I_ERI_Fx2z_Py_Pz_S_bb;
  Double I_ERI_F3y_Dyz_Pz_S_bb = I_ERI_G3yz_Py_Pz_S_bb+ABZ*I_ERI_F3y_Py_Pz_S_bb;
  Double I_ERI_F2yz_Dyz_Pz_S_bb = I_ERI_G2y2z_Py_Pz_S_bb+ABZ*I_ERI_F2yz_Py_Pz_S_bb;
  Double I_ERI_Fy2z_Dyz_Pz_S_bb = I_ERI_Gy3z_Py_Pz_S_bb+ABZ*I_ERI_Fy2z_Py_Pz_S_bb;
  Double I_ERI_F3z_Dyz_Pz_S_bb = I_ERI_G4z_Py_Pz_S_bb+ABZ*I_ERI_F3z_Py_Pz_S_bb;
  Double I_ERI_F3x_D2z_Pz_S_bb = I_ERI_G3xz_Pz_Pz_S_bb+ABZ*I_ERI_F3x_Pz_Pz_S_bb;
  Double I_ERI_F2xy_D2z_Pz_S_bb = I_ERI_G2xyz_Pz_Pz_S_bb+ABZ*I_ERI_F2xy_Pz_Pz_S_bb;
  Double I_ERI_F2xz_D2z_Pz_S_bb = I_ERI_G2x2z_Pz_Pz_S_bb+ABZ*I_ERI_F2xz_Pz_Pz_S_bb;
  Double I_ERI_Fx2y_D2z_Pz_S_bb = I_ERI_Gx2yz_Pz_Pz_S_bb+ABZ*I_ERI_Fx2y_Pz_Pz_S_bb;
  Double I_ERI_Fxyz_D2z_Pz_S_bb = I_ERI_Gxy2z_Pz_Pz_S_bb+ABZ*I_ERI_Fxyz_Pz_Pz_S_bb;
  Double I_ERI_Fx2z_D2z_Pz_S_bb = I_ERI_Gx3z_Pz_Pz_S_bb+ABZ*I_ERI_Fx2z_Pz_Pz_S_bb;
  Double I_ERI_F3y_D2z_Pz_S_bb = I_ERI_G3yz_Pz_Pz_S_bb+ABZ*I_ERI_F3y_Pz_Pz_S_bb;
  Double I_ERI_F2yz_D2z_Pz_S_bb = I_ERI_G2y2z_Pz_Pz_S_bb+ABZ*I_ERI_F2yz_Pz_Pz_S_bb;
  Double I_ERI_Fy2z_D2z_Pz_S_bb = I_ERI_Gy3z_Pz_Pz_S_bb+ABZ*I_ERI_Fy2z_Pz_Pz_S_bb;
  Double I_ERI_F3z_D2z_Pz_S_bb = I_ERI_G4z_Pz_Pz_S_bb+ABZ*I_ERI_F3z_Pz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bc
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_bc = I_ERI_G4x_S_D2x_S_bc+ABX*I_ERI_F3x_S_D2x_S_bc;
  Double I_ERI_F2xy_Px_D2x_S_bc = I_ERI_G3xy_S_D2x_S_bc+ABX*I_ERI_F2xy_S_D2x_S_bc;
  Double I_ERI_F2xz_Px_D2x_S_bc = I_ERI_G3xz_S_D2x_S_bc+ABX*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fx2y_Px_D2x_S_bc = I_ERI_G2x2y_S_D2x_S_bc+ABX*I_ERI_Fx2y_S_D2x_S_bc;
  Double I_ERI_Fxyz_Px_D2x_S_bc = I_ERI_G2xyz_S_D2x_S_bc+ABX*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Px_D2x_S_bc = I_ERI_G2x2z_S_D2x_S_bc+ABX*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F3y_Px_D2x_S_bc = I_ERI_Gx3y_S_D2x_S_bc+ABX*I_ERI_F3y_S_D2x_S_bc;
  Double I_ERI_F2yz_Px_D2x_S_bc = I_ERI_Gx2yz_S_D2x_S_bc+ABX*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Px_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABX*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Px_D2x_S_bc = I_ERI_Gx3z_S_D2x_S_bc+ABX*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F3x_Py_D2x_S_bc = I_ERI_G3xy_S_D2x_S_bc+ABY*I_ERI_F3x_S_D2x_S_bc;
  Double I_ERI_F2xy_Py_D2x_S_bc = I_ERI_G2x2y_S_D2x_S_bc+ABY*I_ERI_F2xy_S_D2x_S_bc;
  Double I_ERI_F2xz_Py_D2x_S_bc = I_ERI_G2xyz_S_D2x_S_bc+ABY*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fx2y_Py_D2x_S_bc = I_ERI_Gx3y_S_D2x_S_bc+ABY*I_ERI_Fx2y_S_D2x_S_bc;
  Double I_ERI_Fxyz_Py_D2x_S_bc = I_ERI_Gx2yz_S_D2x_S_bc+ABY*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Py_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABY*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F3y_Py_D2x_S_bc = I_ERI_G4y_S_D2x_S_bc+ABY*I_ERI_F3y_S_D2x_S_bc;
  Double I_ERI_F2yz_Py_D2x_S_bc = I_ERI_G3yz_S_D2x_S_bc+ABY*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Py_D2x_S_bc = I_ERI_G2y2z_S_D2x_S_bc+ABY*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Py_D2x_S_bc = I_ERI_Gy3z_S_D2x_S_bc+ABY*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F3x_Pz_D2x_S_bc = I_ERI_G3xz_S_D2x_S_bc+ABZ*I_ERI_F3x_S_D2x_S_bc;
  Double I_ERI_F2xy_Pz_D2x_S_bc = I_ERI_G2xyz_S_D2x_S_bc+ABZ*I_ERI_F2xy_S_D2x_S_bc;
  Double I_ERI_F2xz_Pz_D2x_S_bc = I_ERI_G2x2z_S_D2x_S_bc+ABZ*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fx2y_Pz_D2x_S_bc = I_ERI_Gx2yz_S_D2x_S_bc+ABZ*I_ERI_Fx2y_S_D2x_S_bc;
  Double I_ERI_Fxyz_Pz_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABZ*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Pz_D2x_S_bc = I_ERI_Gx3z_S_D2x_S_bc+ABZ*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F3y_Pz_D2x_S_bc = I_ERI_G3yz_S_D2x_S_bc+ABZ*I_ERI_F3y_S_D2x_S_bc;
  Double I_ERI_F2yz_Pz_D2x_S_bc = I_ERI_G2y2z_S_D2x_S_bc+ABZ*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Pz_D2x_S_bc = I_ERI_Gy3z_S_D2x_S_bc+ABZ*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Pz_D2x_S_bc = I_ERI_G4z_S_D2x_S_bc+ABZ*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F3x_Px_Dxy_S_bc = I_ERI_G4x_S_Dxy_S_bc+ABX*I_ERI_F3x_S_Dxy_S_bc;
  Double I_ERI_F2xy_Px_Dxy_S_bc = I_ERI_G3xy_S_Dxy_S_bc+ABX*I_ERI_F2xy_S_Dxy_S_bc;
  Double I_ERI_F2xz_Px_Dxy_S_bc = I_ERI_G3xz_S_Dxy_S_bc+ABX*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fx2y_Px_Dxy_S_bc = I_ERI_G2x2y_S_Dxy_S_bc+ABX*I_ERI_Fx2y_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Px_Dxy_S_bc = I_ERI_G2xyz_S_Dxy_S_bc+ABX*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Px_Dxy_S_bc = I_ERI_G2x2z_S_Dxy_S_bc+ABX*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F3y_Px_Dxy_S_bc = I_ERI_Gx3y_S_Dxy_S_bc+ABX*I_ERI_F3y_S_Dxy_S_bc;
  Double I_ERI_F2yz_Px_Dxy_S_bc = I_ERI_Gx2yz_S_Dxy_S_bc+ABX*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Px_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABX*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Px_Dxy_S_bc = I_ERI_Gx3z_S_Dxy_S_bc+ABX*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F3x_Py_Dxy_S_bc = I_ERI_G3xy_S_Dxy_S_bc+ABY*I_ERI_F3x_S_Dxy_S_bc;
  Double I_ERI_F2xy_Py_Dxy_S_bc = I_ERI_G2x2y_S_Dxy_S_bc+ABY*I_ERI_F2xy_S_Dxy_S_bc;
  Double I_ERI_F2xz_Py_Dxy_S_bc = I_ERI_G2xyz_S_Dxy_S_bc+ABY*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fx2y_Py_Dxy_S_bc = I_ERI_Gx3y_S_Dxy_S_bc+ABY*I_ERI_Fx2y_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Py_Dxy_S_bc = I_ERI_Gx2yz_S_Dxy_S_bc+ABY*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Py_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABY*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F3y_Py_Dxy_S_bc = I_ERI_G4y_S_Dxy_S_bc+ABY*I_ERI_F3y_S_Dxy_S_bc;
  Double I_ERI_F2yz_Py_Dxy_S_bc = I_ERI_G3yz_S_Dxy_S_bc+ABY*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Py_Dxy_S_bc = I_ERI_G2y2z_S_Dxy_S_bc+ABY*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Py_Dxy_S_bc = I_ERI_Gy3z_S_Dxy_S_bc+ABY*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F3x_Pz_Dxy_S_bc = I_ERI_G3xz_S_Dxy_S_bc+ABZ*I_ERI_F3x_S_Dxy_S_bc;
  Double I_ERI_F2xy_Pz_Dxy_S_bc = I_ERI_G2xyz_S_Dxy_S_bc+ABZ*I_ERI_F2xy_S_Dxy_S_bc;
  Double I_ERI_F2xz_Pz_Dxy_S_bc = I_ERI_G2x2z_S_Dxy_S_bc+ABZ*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fx2y_Pz_Dxy_S_bc = I_ERI_Gx2yz_S_Dxy_S_bc+ABZ*I_ERI_Fx2y_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Pz_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABZ*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Pz_Dxy_S_bc = I_ERI_Gx3z_S_Dxy_S_bc+ABZ*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F3y_Pz_Dxy_S_bc = I_ERI_G3yz_S_Dxy_S_bc+ABZ*I_ERI_F3y_S_Dxy_S_bc;
  Double I_ERI_F2yz_Pz_Dxy_S_bc = I_ERI_G2y2z_S_Dxy_S_bc+ABZ*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Pz_Dxy_S_bc = I_ERI_Gy3z_S_Dxy_S_bc+ABZ*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Pz_Dxy_S_bc = I_ERI_G4z_S_Dxy_S_bc+ABZ*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F3x_Px_Dxz_S_bc = I_ERI_G4x_S_Dxz_S_bc+ABX*I_ERI_F3x_S_Dxz_S_bc;
  Double I_ERI_F2xy_Px_Dxz_S_bc = I_ERI_G3xy_S_Dxz_S_bc+ABX*I_ERI_F2xy_S_Dxz_S_bc;
  Double I_ERI_F2xz_Px_Dxz_S_bc = I_ERI_G3xz_S_Dxz_S_bc+ABX*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fx2y_Px_Dxz_S_bc = I_ERI_G2x2y_S_Dxz_S_bc+ABX*I_ERI_Fx2y_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Px_Dxz_S_bc = I_ERI_G2xyz_S_Dxz_S_bc+ABX*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Px_Dxz_S_bc = I_ERI_G2x2z_S_Dxz_S_bc+ABX*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F3y_Px_Dxz_S_bc = I_ERI_Gx3y_S_Dxz_S_bc+ABX*I_ERI_F3y_S_Dxz_S_bc;
  Double I_ERI_F2yz_Px_Dxz_S_bc = I_ERI_Gx2yz_S_Dxz_S_bc+ABX*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Px_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABX*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Px_Dxz_S_bc = I_ERI_Gx3z_S_Dxz_S_bc+ABX*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F3x_Py_Dxz_S_bc = I_ERI_G3xy_S_Dxz_S_bc+ABY*I_ERI_F3x_S_Dxz_S_bc;
  Double I_ERI_F2xy_Py_Dxz_S_bc = I_ERI_G2x2y_S_Dxz_S_bc+ABY*I_ERI_F2xy_S_Dxz_S_bc;
  Double I_ERI_F2xz_Py_Dxz_S_bc = I_ERI_G2xyz_S_Dxz_S_bc+ABY*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fx2y_Py_Dxz_S_bc = I_ERI_Gx3y_S_Dxz_S_bc+ABY*I_ERI_Fx2y_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Py_Dxz_S_bc = I_ERI_Gx2yz_S_Dxz_S_bc+ABY*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Py_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABY*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F3y_Py_Dxz_S_bc = I_ERI_G4y_S_Dxz_S_bc+ABY*I_ERI_F3y_S_Dxz_S_bc;
  Double I_ERI_F2yz_Py_Dxz_S_bc = I_ERI_G3yz_S_Dxz_S_bc+ABY*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Py_Dxz_S_bc = I_ERI_G2y2z_S_Dxz_S_bc+ABY*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Py_Dxz_S_bc = I_ERI_Gy3z_S_Dxz_S_bc+ABY*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F3x_Pz_Dxz_S_bc = I_ERI_G3xz_S_Dxz_S_bc+ABZ*I_ERI_F3x_S_Dxz_S_bc;
  Double I_ERI_F2xy_Pz_Dxz_S_bc = I_ERI_G2xyz_S_Dxz_S_bc+ABZ*I_ERI_F2xy_S_Dxz_S_bc;
  Double I_ERI_F2xz_Pz_Dxz_S_bc = I_ERI_G2x2z_S_Dxz_S_bc+ABZ*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fx2y_Pz_Dxz_S_bc = I_ERI_Gx2yz_S_Dxz_S_bc+ABZ*I_ERI_Fx2y_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Pz_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABZ*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Pz_Dxz_S_bc = I_ERI_Gx3z_S_Dxz_S_bc+ABZ*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F3y_Pz_Dxz_S_bc = I_ERI_G3yz_S_Dxz_S_bc+ABZ*I_ERI_F3y_S_Dxz_S_bc;
  Double I_ERI_F2yz_Pz_Dxz_S_bc = I_ERI_G2y2z_S_Dxz_S_bc+ABZ*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Pz_Dxz_S_bc = I_ERI_Gy3z_S_Dxz_S_bc+ABZ*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Pz_Dxz_S_bc = I_ERI_G4z_S_Dxz_S_bc+ABZ*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F3x_Px_D2y_S_bc = I_ERI_G4x_S_D2y_S_bc+ABX*I_ERI_F3x_S_D2y_S_bc;
  Double I_ERI_F2xy_Px_D2y_S_bc = I_ERI_G3xy_S_D2y_S_bc+ABX*I_ERI_F2xy_S_D2y_S_bc;
  Double I_ERI_F2xz_Px_D2y_S_bc = I_ERI_G3xz_S_D2y_S_bc+ABX*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fx2y_Px_D2y_S_bc = I_ERI_G2x2y_S_D2y_S_bc+ABX*I_ERI_Fx2y_S_D2y_S_bc;
  Double I_ERI_Fxyz_Px_D2y_S_bc = I_ERI_G2xyz_S_D2y_S_bc+ABX*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Px_D2y_S_bc = I_ERI_G2x2z_S_D2y_S_bc+ABX*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F3y_Px_D2y_S_bc = I_ERI_Gx3y_S_D2y_S_bc+ABX*I_ERI_F3y_S_D2y_S_bc;
  Double I_ERI_F2yz_Px_D2y_S_bc = I_ERI_Gx2yz_S_D2y_S_bc+ABX*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Px_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABX*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Px_D2y_S_bc = I_ERI_Gx3z_S_D2y_S_bc+ABX*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F3x_Py_D2y_S_bc = I_ERI_G3xy_S_D2y_S_bc+ABY*I_ERI_F3x_S_D2y_S_bc;
  Double I_ERI_F2xy_Py_D2y_S_bc = I_ERI_G2x2y_S_D2y_S_bc+ABY*I_ERI_F2xy_S_D2y_S_bc;
  Double I_ERI_F2xz_Py_D2y_S_bc = I_ERI_G2xyz_S_D2y_S_bc+ABY*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fx2y_Py_D2y_S_bc = I_ERI_Gx3y_S_D2y_S_bc+ABY*I_ERI_Fx2y_S_D2y_S_bc;
  Double I_ERI_Fxyz_Py_D2y_S_bc = I_ERI_Gx2yz_S_D2y_S_bc+ABY*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Py_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABY*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F3y_Py_D2y_S_bc = I_ERI_G4y_S_D2y_S_bc+ABY*I_ERI_F3y_S_D2y_S_bc;
  Double I_ERI_F2yz_Py_D2y_S_bc = I_ERI_G3yz_S_D2y_S_bc+ABY*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Py_D2y_S_bc = I_ERI_G2y2z_S_D2y_S_bc+ABY*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Py_D2y_S_bc = I_ERI_Gy3z_S_D2y_S_bc+ABY*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F3x_Pz_D2y_S_bc = I_ERI_G3xz_S_D2y_S_bc+ABZ*I_ERI_F3x_S_D2y_S_bc;
  Double I_ERI_F2xy_Pz_D2y_S_bc = I_ERI_G2xyz_S_D2y_S_bc+ABZ*I_ERI_F2xy_S_D2y_S_bc;
  Double I_ERI_F2xz_Pz_D2y_S_bc = I_ERI_G2x2z_S_D2y_S_bc+ABZ*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fx2y_Pz_D2y_S_bc = I_ERI_Gx2yz_S_D2y_S_bc+ABZ*I_ERI_Fx2y_S_D2y_S_bc;
  Double I_ERI_Fxyz_Pz_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABZ*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Pz_D2y_S_bc = I_ERI_Gx3z_S_D2y_S_bc+ABZ*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F3y_Pz_D2y_S_bc = I_ERI_G3yz_S_D2y_S_bc+ABZ*I_ERI_F3y_S_D2y_S_bc;
  Double I_ERI_F2yz_Pz_D2y_S_bc = I_ERI_G2y2z_S_D2y_S_bc+ABZ*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Pz_D2y_S_bc = I_ERI_Gy3z_S_D2y_S_bc+ABZ*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Pz_D2y_S_bc = I_ERI_G4z_S_D2y_S_bc+ABZ*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F3x_Px_Dyz_S_bc = I_ERI_G4x_S_Dyz_S_bc+ABX*I_ERI_F3x_S_Dyz_S_bc;
  Double I_ERI_F2xy_Px_Dyz_S_bc = I_ERI_G3xy_S_Dyz_S_bc+ABX*I_ERI_F2xy_S_Dyz_S_bc;
  Double I_ERI_F2xz_Px_Dyz_S_bc = I_ERI_G3xz_S_Dyz_S_bc+ABX*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fx2y_Px_Dyz_S_bc = I_ERI_G2x2y_S_Dyz_S_bc+ABX*I_ERI_Fx2y_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Px_Dyz_S_bc = I_ERI_G2xyz_S_Dyz_S_bc+ABX*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Px_Dyz_S_bc = I_ERI_G2x2z_S_Dyz_S_bc+ABX*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F3y_Px_Dyz_S_bc = I_ERI_Gx3y_S_Dyz_S_bc+ABX*I_ERI_F3y_S_Dyz_S_bc;
  Double I_ERI_F2yz_Px_Dyz_S_bc = I_ERI_Gx2yz_S_Dyz_S_bc+ABX*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Px_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABX*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Px_Dyz_S_bc = I_ERI_Gx3z_S_Dyz_S_bc+ABX*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F3x_Py_Dyz_S_bc = I_ERI_G3xy_S_Dyz_S_bc+ABY*I_ERI_F3x_S_Dyz_S_bc;
  Double I_ERI_F2xy_Py_Dyz_S_bc = I_ERI_G2x2y_S_Dyz_S_bc+ABY*I_ERI_F2xy_S_Dyz_S_bc;
  Double I_ERI_F2xz_Py_Dyz_S_bc = I_ERI_G2xyz_S_Dyz_S_bc+ABY*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fx2y_Py_Dyz_S_bc = I_ERI_Gx3y_S_Dyz_S_bc+ABY*I_ERI_Fx2y_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Py_Dyz_S_bc = I_ERI_Gx2yz_S_Dyz_S_bc+ABY*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Py_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABY*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F3y_Py_Dyz_S_bc = I_ERI_G4y_S_Dyz_S_bc+ABY*I_ERI_F3y_S_Dyz_S_bc;
  Double I_ERI_F2yz_Py_Dyz_S_bc = I_ERI_G3yz_S_Dyz_S_bc+ABY*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Py_Dyz_S_bc = I_ERI_G2y2z_S_Dyz_S_bc+ABY*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Py_Dyz_S_bc = I_ERI_Gy3z_S_Dyz_S_bc+ABY*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F3x_Pz_Dyz_S_bc = I_ERI_G3xz_S_Dyz_S_bc+ABZ*I_ERI_F3x_S_Dyz_S_bc;
  Double I_ERI_F2xy_Pz_Dyz_S_bc = I_ERI_G2xyz_S_Dyz_S_bc+ABZ*I_ERI_F2xy_S_Dyz_S_bc;
  Double I_ERI_F2xz_Pz_Dyz_S_bc = I_ERI_G2x2z_S_Dyz_S_bc+ABZ*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fx2y_Pz_Dyz_S_bc = I_ERI_Gx2yz_S_Dyz_S_bc+ABZ*I_ERI_Fx2y_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Pz_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABZ*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Pz_Dyz_S_bc = I_ERI_Gx3z_S_Dyz_S_bc+ABZ*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F3y_Pz_Dyz_S_bc = I_ERI_G3yz_S_Dyz_S_bc+ABZ*I_ERI_F3y_S_Dyz_S_bc;
  Double I_ERI_F2yz_Pz_Dyz_S_bc = I_ERI_G2y2z_S_Dyz_S_bc+ABZ*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Pz_Dyz_S_bc = I_ERI_Gy3z_S_Dyz_S_bc+ABZ*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Pz_Dyz_S_bc = I_ERI_G4z_S_Dyz_S_bc+ABZ*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F3x_Px_D2z_S_bc = I_ERI_G4x_S_D2z_S_bc+ABX*I_ERI_F3x_S_D2z_S_bc;
  Double I_ERI_F2xy_Px_D2z_S_bc = I_ERI_G3xy_S_D2z_S_bc+ABX*I_ERI_F2xy_S_D2z_S_bc;
  Double I_ERI_F2xz_Px_D2z_S_bc = I_ERI_G3xz_S_D2z_S_bc+ABX*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fx2y_Px_D2z_S_bc = I_ERI_G2x2y_S_D2z_S_bc+ABX*I_ERI_Fx2y_S_D2z_S_bc;
  Double I_ERI_Fxyz_Px_D2z_S_bc = I_ERI_G2xyz_S_D2z_S_bc+ABX*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Px_D2z_S_bc = I_ERI_G2x2z_S_D2z_S_bc+ABX*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F3y_Px_D2z_S_bc = I_ERI_Gx3y_S_D2z_S_bc+ABX*I_ERI_F3y_S_D2z_S_bc;
  Double I_ERI_F2yz_Px_D2z_S_bc = I_ERI_Gx2yz_S_D2z_S_bc+ABX*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Px_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABX*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Px_D2z_S_bc = I_ERI_Gx3z_S_D2z_S_bc+ABX*I_ERI_F3z_S_D2z_S_bc;
  Double I_ERI_F3x_Py_D2z_S_bc = I_ERI_G3xy_S_D2z_S_bc+ABY*I_ERI_F3x_S_D2z_S_bc;
  Double I_ERI_F2xy_Py_D2z_S_bc = I_ERI_G2x2y_S_D2z_S_bc+ABY*I_ERI_F2xy_S_D2z_S_bc;
  Double I_ERI_F2xz_Py_D2z_S_bc = I_ERI_G2xyz_S_D2z_S_bc+ABY*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fx2y_Py_D2z_S_bc = I_ERI_Gx3y_S_D2z_S_bc+ABY*I_ERI_Fx2y_S_D2z_S_bc;
  Double I_ERI_Fxyz_Py_D2z_S_bc = I_ERI_Gx2yz_S_D2z_S_bc+ABY*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Py_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABY*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F3y_Py_D2z_S_bc = I_ERI_G4y_S_D2z_S_bc+ABY*I_ERI_F3y_S_D2z_S_bc;
  Double I_ERI_F2yz_Py_D2z_S_bc = I_ERI_G3yz_S_D2z_S_bc+ABY*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Py_D2z_S_bc = I_ERI_G2y2z_S_D2z_S_bc+ABY*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Py_D2z_S_bc = I_ERI_Gy3z_S_D2z_S_bc+ABY*I_ERI_F3z_S_D2z_S_bc;
  Double I_ERI_F3x_Pz_D2z_S_bc = I_ERI_G3xz_S_D2z_S_bc+ABZ*I_ERI_F3x_S_D2z_S_bc;
  Double I_ERI_F2xy_Pz_D2z_S_bc = I_ERI_G2xyz_S_D2z_S_bc+ABZ*I_ERI_F2xy_S_D2z_S_bc;
  Double I_ERI_F2xz_Pz_D2z_S_bc = I_ERI_G2x2z_S_D2z_S_bc+ABZ*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fx2y_Pz_D2z_S_bc = I_ERI_Gx2yz_S_D2z_S_bc+ABZ*I_ERI_Fx2y_S_D2z_S_bc;
  Double I_ERI_Fxyz_Pz_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABZ*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Pz_D2z_S_bc = I_ERI_Gx3z_S_D2z_S_bc+ABZ*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F3y_Pz_D2z_S_bc = I_ERI_G3yz_S_D2z_S_bc+ABZ*I_ERI_F3y_S_D2z_S_bc;
  Double I_ERI_F2yz_Pz_D2z_S_bc = I_ERI_G2y2z_S_D2z_S_bc+ABZ*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Pz_D2z_S_bc = I_ERI_Gy3z_S_D2z_S_bc+ABZ*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Pz_D2z_S_bc = I_ERI_G4z_S_D2z_S_bc+ABZ*I_ERI_F3z_S_D2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_P_bd
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P_bd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_bd
   ************************************************************/
  Double I_ERI_F3x_Px_Px_Px_bd = I_ERI_G4x_S_Px_Px_bd+ABX*I_ERI_F3x_S_Px_Px_bd;
  Double I_ERI_F2xy_Px_Px_Px_bd = I_ERI_G3xy_S_Px_Px_bd+ABX*I_ERI_F2xy_S_Px_Px_bd;
  Double I_ERI_F2xz_Px_Px_Px_bd = I_ERI_G3xz_S_Px_Px_bd+ABX*I_ERI_F2xz_S_Px_Px_bd;
  Double I_ERI_Fx2y_Px_Px_Px_bd = I_ERI_G2x2y_S_Px_Px_bd+ABX*I_ERI_Fx2y_S_Px_Px_bd;
  Double I_ERI_Fxyz_Px_Px_Px_bd = I_ERI_G2xyz_S_Px_Px_bd+ABX*I_ERI_Fxyz_S_Px_Px_bd;
  Double I_ERI_Fx2z_Px_Px_Px_bd = I_ERI_G2x2z_S_Px_Px_bd+ABX*I_ERI_Fx2z_S_Px_Px_bd;
  Double I_ERI_F3y_Px_Px_Px_bd = I_ERI_Gx3y_S_Px_Px_bd+ABX*I_ERI_F3y_S_Px_Px_bd;
  Double I_ERI_F2yz_Px_Px_Px_bd = I_ERI_Gx2yz_S_Px_Px_bd+ABX*I_ERI_F2yz_S_Px_Px_bd;
  Double I_ERI_Fy2z_Px_Px_Px_bd = I_ERI_Gxy2z_S_Px_Px_bd+ABX*I_ERI_Fy2z_S_Px_Px_bd;
  Double I_ERI_F3z_Px_Px_Px_bd = I_ERI_Gx3z_S_Px_Px_bd+ABX*I_ERI_F3z_S_Px_Px_bd;
  Double I_ERI_F3x_Py_Px_Px_bd = I_ERI_G3xy_S_Px_Px_bd+ABY*I_ERI_F3x_S_Px_Px_bd;
  Double I_ERI_F2xy_Py_Px_Px_bd = I_ERI_G2x2y_S_Px_Px_bd+ABY*I_ERI_F2xy_S_Px_Px_bd;
  Double I_ERI_F2xz_Py_Px_Px_bd = I_ERI_G2xyz_S_Px_Px_bd+ABY*I_ERI_F2xz_S_Px_Px_bd;
  Double I_ERI_Fx2y_Py_Px_Px_bd = I_ERI_Gx3y_S_Px_Px_bd+ABY*I_ERI_Fx2y_S_Px_Px_bd;
  Double I_ERI_Fxyz_Py_Px_Px_bd = I_ERI_Gx2yz_S_Px_Px_bd+ABY*I_ERI_Fxyz_S_Px_Px_bd;
  Double I_ERI_Fx2z_Py_Px_Px_bd = I_ERI_Gxy2z_S_Px_Px_bd+ABY*I_ERI_Fx2z_S_Px_Px_bd;
  Double I_ERI_F3y_Py_Px_Px_bd = I_ERI_G4y_S_Px_Px_bd+ABY*I_ERI_F3y_S_Px_Px_bd;
  Double I_ERI_F2yz_Py_Px_Px_bd = I_ERI_G3yz_S_Px_Px_bd+ABY*I_ERI_F2yz_S_Px_Px_bd;
  Double I_ERI_Fy2z_Py_Px_Px_bd = I_ERI_G2y2z_S_Px_Px_bd+ABY*I_ERI_Fy2z_S_Px_Px_bd;
  Double I_ERI_F3z_Py_Px_Px_bd = I_ERI_Gy3z_S_Px_Px_bd+ABY*I_ERI_F3z_S_Px_Px_bd;
  Double I_ERI_F3x_Pz_Px_Px_bd = I_ERI_G3xz_S_Px_Px_bd+ABZ*I_ERI_F3x_S_Px_Px_bd;
  Double I_ERI_F2xy_Pz_Px_Px_bd = I_ERI_G2xyz_S_Px_Px_bd+ABZ*I_ERI_F2xy_S_Px_Px_bd;
  Double I_ERI_F2xz_Pz_Px_Px_bd = I_ERI_G2x2z_S_Px_Px_bd+ABZ*I_ERI_F2xz_S_Px_Px_bd;
  Double I_ERI_Fx2y_Pz_Px_Px_bd = I_ERI_Gx2yz_S_Px_Px_bd+ABZ*I_ERI_Fx2y_S_Px_Px_bd;
  Double I_ERI_Fxyz_Pz_Px_Px_bd = I_ERI_Gxy2z_S_Px_Px_bd+ABZ*I_ERI_Fxyz_S_Px_Px_bd;
  Double I_ERI_Fx2z_Pz_Px_Px_bd = I_ERI_Gx3z_S_Px_Px_bd+ABZ*I_ERI_Fx2z_S_Px_Px_bd;
  Double I_ERI_F3y_Pz_Px_Px_bd = I_ERI_G3yz_S_Px_Px_bd+ABZ*I_ERI_F3y_S_Px_Px_bd;
  Double I_ERI_F2yz_Pz_Px_Px_bd = I_ERI_G2y2z_S_Px_Px_bd+ABZ*I_ERI_F2yz_S_Px_Px_bd;
  Double I_ERI_Fy2z_Pz_Px_Px_bd = I_ERI_Gy3z_S_Px_Px_bd+ABZ*I_ERI_Fy2z_S_Px_Px_bd;
  Double I_ERI_F3z_Pz_Px_Px_bd = I_ERI_G4z_S_Px_Px_bd+ABZ*I_ERI_F3z_S_Px_Px_bd;
  Double I_ERI_F3x_Px_Py_Px_bd = I_ERI_G4x_S_Py_Px_bd+ABX*I_ERI_F3x_S_Py_Px_bd;
  Double I_ERI_F2xy_Px_Py_Px_bd = I_ERI_G3xy_S_Py_Px_bd+ABX*I_ERI_F2xy_S_Py_Px_bd;
  Double I_ERI_F2xz_Px_Py_Px_bd = I_ERI_G3xz_S_Py_Px_bd+ABX*I_ERI_F2xz_S_Py_Px_bd;
  Double I_ERI_Fx2y_Px_Py_Px_bd = I_ERI_G2x2y_S_Py_Px_bd+ABX*I_ERI_Fx2y_S_Py_Px_bd;
  Double I_ERI_Fxyz_Px_Py_Px_bd = I_ERI_G2xyz_S_Py_Px_bd+ABX*I_ERI_Fxyz_S_Py_Px_bd;
  Double I_ERI_Fx2z_Px_Py_Px_bd = I_ERI_G2x2z_S_Py_Px_bd+ABX*I_ERI_Fx2z_S_Py_Px_bd;
  Double I_ERI_F3y_Px_Py_Px_bd = I_ERI_Gx3y_S_Py_Px_bd+ABX*I_ERI_F3y_S_Py_Px_bd;
  Double I_ERI_F2yz_Px_Py_Px_bd = I_ERI_Gx2yz_S_Py_Px_bd+ABX*I_ERI_F2yz_S_Py_Px_bd;
  Double I_ERI_Fy2z_Px_Py_Px_bd = I_ERI_Gxy2z_S_Py_Px_bd+ABX*I_ERI_Fy2z_S_Py_Px_bd;
  Double I_ERI_F3z_Px_Py_Px_bd = I_ERI_Gx3z_S_Py_Px_bd+ABX*I_ERI_F3z_S_Py_Px_bd;
  Double I_ERI_F3x_Py_Py_Px_bd = I_ERI_G3xy_S_Py_Px_bd+ABY*I_ERI_F3x_S_Py_Px_bd;
  Double I_ERI_F2xy_Py_Py_Px_bd = I_ERI_G2x2y_S_Py_Px_bd+ABY*I_ERI_F2xy_S_Py_Px_bd;
  Double I_ERI_F2xz_Py_Py_Px_bd = I_ERI_G2xyz_S_Py_Px_bd+ABY*I_ERI_F2xz_S_Py_Px_bd;
  Double I_ERI_Fx2y_Py_Py_Px_bd = I_ERI_Gx3y_S_Py_Px_bd+ABY*I_ERI_Fx2y_S_Py_Px_bd;
  Double I_ERI_Fxyz_Py_Py_Px_bd = I_ERI_Gx2yz_S_Py_Px_bd+ABY*I_ERI_Fxyz_S_Py_Px_bd;
  Double I_ERI_Fx2z_Py_Py_Px_bd = I_ERI_Gxy2z_S_Py_Px_bd+ABY*I_ERI_Fx2z_S_Py_Px_bd;
  Double I_ERI_F3y_Py_Py_Px_bd = I_ERI_G4y_S_Py_Px_bd+ABY*I_ERI_F3y_S_Py_Px_bd;
  Double I_ERI_F2yz_Py_Py_Px_bd = I_ERI_G3yz_S_Py_Px_bd+ABY*I_ERI_F2yz_S_Py_Px_bd;
  Double I_ERI_Fy2z_Py_Py_Px_bd = I_ERI_G2y2z_S_Py_Px_bd+ABY*I_ERI_Fy2z_S_Py_Px_bd;
  Double I_ERI_F3z_Py_Py_Px_bd = I_ERI_Gy3z_S_Py_Px_bd+ABY*I_ERI_F3z_S_Py_Px_bd;
  Double I_ERI_F3x_Pz_Py_Px_bd = I_ERI_G3xz_S_Py_Px_bd+ABZ*I_ERI_F3x_S_Py_Px_bd;
  Double I_ERI_F2xy_Pz_Py_Px_bd = I_ERI_G2xyz_S_Py_Px_bd+ABZ*I_ERI_F2xy_S_Py_Px_bd;
  Double I_ERI_F2xz_Pz_Py_Px_bd = I_ERI_G2x2z_S_Py_Px_bd+ABZ*I_ERI_F2xz_S_Py_Px_bd;
  Double I_ERI_Fx2y_Pz_Py_Px_bd = I_ERI_Gx2yz_S_Py_Px_bd+ABZ*I_ERI_Fx2y_S_Py_Px_bd;
  Double I_ERI_Fxyz_Pz_Py_Px_bd = I_ERI_Gxy2z_S_Py_Px_bd+ABZ*I_ERI_Fxyz_S_Py_Px_bd;
  Double I_ERI_Fx2z_Pz_Py_Px_bd = I_ERI_Gx3z_S_Py_Px_bd+ABZ*I_ERI_Fx2z_S_Py_Px_bd;
  Double I_ERI_F3y_Pz_Py_Px_bd = I_ERI_G3yz_S_Py_Px_bd+ABZ*I_ERI_F3y_S_Py_Px_bd;
  Double I_ERI_F2yz_Pz_Py_Px_bd = I_ERI_G2y2z_S_Py_Px_bd+ABZ*I_ERI_F2yz_S_Py_Px_bd;
  Double I_ERI_Fy2z_Pz_Py_Px_bd = I_ERI_Gy3z_S_Py_Px_bd+ABZ*I_ERI_Fy2z_S_Py_Px_bd;
  Double I_ERI_F3z_Pz_Py_Px_bd = I_ERI_G4z_S_Py_Px_bd+ABZ*I_ERI_F3z_S_Py_Px_bd;
  Double I_ERI_F3x_Px_Pz_Px_bd = I_ERI_G4x_S_Pz_Px_bd+ABX*I_ERI_F3x_S_Pz_Px_bd;
  Double I_ERI_F2xy_Px_Pz_Px_bd = I_ERI_G3xy_S_Pz_Px_bd+ABX*I_ERI_F2xy_S_Pz_Px_bd;
  Double I_ERI_F2xz_Px_Pz_Px_bd = I_ERI_G3xz_S_Pz_Px_bd+ABX*I_ERI_F2xz_S_Pz_Px_bd;
  Double I_ERI_Fx2y_Px_Pz_Px_bd = I_ERI_G2x2y_S_Pz_Px_bd+ABX*I_ERI_Fx2y_S_Pz_Px_bd;
  Double I_ERI_Fxyz_Px_Pz_Px_bd = I_ERI_G2xyz_S_Pz_Px_bd+ABX*I_ERI_Fxyz_S_Pz_Px_bd;
  Double I_ERI_Fx2z_Px_Pz_Px_bd = I_ERI_G2x2z_S_Pz_Px_bd+ABX*I_ERI_Fx2z_S_Pz_Px_bd;
  Double I_ERI_F3y_Px_Pz_Px_bd = I_ERI_Gx3y_S_Pz_Px_bd+ABX*I_ERI_F3y_S_Pz_Px_bd;
  Double I_ERI_F2yz_Px_Pz_Px_bd = I_ERI_Gx2yz_S_Pz_Px_bd+ABX*I_ERI_F2yz_S_Pz_Px_bd;
  Double I_ERI_Fy2z_Px_Pz_Px_bd = I_ERI_Gxy2z_S_Pz_Px_bd+ABX*I_ERI_Fy2z_S_Pz_Px_bd;
  Double I_ERI_F3z_Px_Pz_Px_bd = I_ERI_Gx3z_S_Pz_Px_bd+ABX*I_ERI_F3z_S_Pz_Px_bd;
  Double I_ERI_F3x_Py_Pz_Px_bd = I_ERI_G3xy_S_Pz_Px_bd+ABY*I_ERI_F3x_S_Pz_Px_bd;
  Double I_ERI_F2xy_Py_Pz_Px_bd = I_ERI_G2x2y_S_Pz_Px_bd+ABY*I_ERI_F2xy_S_Pz_Px_bd;
  Double I_ERI_F2xz_Py_Pz_Px_bd = I_ERI_G2xyz_S_Pz_Px_bd+ABY*I_ERI_F2xz_S_Pz_Px_bd;
  Double I_ERI_Fx2y_Py_Pz_Px_bd = I_ERI_Gx3y_S_Pz_Px_bd+ABY*I_ERI_Fx2y_S_Pz_Px_bd;
  Double I_ERI_Fxyz_Py_Pz_Px_bd = I_ERI_Gx2yz_S_Pz_Px_bd+ABY*I_ERI_Fxyz_S_Pz_Px_bd;
  Double I_ERI_Fx2z_Py_Pz_Px_bd = I_ERI_Gxy2z_S_Pz_Px_bd+ABY*I_ERI_Fx2z_S_Pz_Px_bd;
  Double I_ERI_F3y_Py_Pz_Px_bd = I_ERI_G4y_S_Pz_Px_bd+ABY*I_ERI_F3y_S_Pz_Px_bd;
  Double I_ERI_F2yz_Py_Pz_Px_bd = I_ERI_G3yz_S_Pz_Px_bd+ABY*I_ERI_F2yz_S_Pz_Px_bd;
  Double I_ERI_Fy2z_Py_Pz_Px_bd = I_ERI_G2y2z_S_Pz_Px_bd+ABY*I_ERI_Fy2z_S_Pz_Px_bd;
  Double I_ERI_F3z_Py_Pz_Px_bd = I_ERI_Gy3z_S_Pz_Px_bd+ABY*I_ERI_F3z_S_Pz_Px_bd;
  Double I_ERI_F3x_Pz_Pz_Px_bd = I_ERI_G3xz_S_Pz_Px_bd+ABZ*I_ERI_F3x_S_Pz_Px_bd;
  Double I_ERI_F2xy_Pz_Pz_Px_bd = I_ERI_G2xyz_S_Pz_Px_bd+ABZ*I_ERI_F2xy_S_Pz_Px_bd;
  Double I_ERI_F2xz_Pz_Pz_Px_bd = I_ERI_G2x2z_S_Pz_Px_bd+ABZ*I_ERI_F2xz_S_Pz_Px_bd;
  Double I_ERI_Fx2y_Pz_Pz_Px_bd = I_ERI_Gx2yz_S_Pz_Px_bd+ABZ*I_ERI_Fx2y_S_Pz_Px_bd;
  Double I_ERI_Fxyz_Pz_Pz_Px_bd = I_ERI_Gxy2z_S_Pz_Px_bd+ABZ*I_ERI_Fxyz_S_Pz_Px_bd;
  Double I_ERI_Fx2z_Pz_Pz_Px_bd = I_ERI_Gx3z_S_Pz_Px_bd+ABZ*I_ERI_Fx2z_S_Pz_Px_bd;
  Double I_ERI_F3y_Pz_Pz_Px_bd = I_ERI_G3yz_S_Pz_Px_bd+ABZ*I_ERI_F3y_S_Pz_Px_bd;
  Double I_ERI_F2yz_Pz_Pz_Px_bd = I_ERI_G2y2z_S_Pz_Px_bd+ABZ*I_ERI_F2yz_S_Pz_Px_bd;
  Double I_ERI_Fy2z_Pz_Pz_Px_bd = I_ERI_Gy3z_S_Pz_Px_bd+ABZ*I_ERI_Fy2z_S_Pz_Px_bd;
  Double I_ERI_F3z_Pz_Pz_Px_bd = I_ERI_G4z_S_Pz_Px_bd+ABZ*I_ERI_F3z_S_Pz_Px_bd;
  Double I_ERI_F3x_Px_Px_Py_bd = I_ERI_G4x_S_Px_Py_bd+ABX*I_ERI_F3x_S_Px_Py_bd;
  Double I_ERI_F2xy_Px_Px_Py_bd = I_ERI_G3xy_S_Px_Py_bd+ABX*I_ERI_F2xy_S_Px_Py_bd;
  Double I_ERI_F2xz_Px_Px_Py_bd = I_ERI_G3xz_S_Px_Py_bd+ABX*I_ERI_F2xz_S_Px_Py_bd;
  Double I_ERI_Fx2y_Px_Px_Py_bd = I_ERI_G2x2y_S_Px_Py_bd+ABX*I_ERI_Fx2y_S_Px_Py_bd;
  Double I_ERI_Fxyz_Px_Px_Py_bd = I_ERI_G2xyz_S_Px_Py_bd+ABX*I_ERI_Fxyz_S_Px_Py_bd;
  Double I_ERI_Fx2z_Px_Px_Py_bd = I_ERI_G2x2z_S_Px_Py_bd+ABX*I_ERI_Fx2z_S_Px_Py_bd;
  Double I_ERI_F3y_Px_Px_Py_bd = I_ERI_Gx3y_S_Px_Py_bd+ABX*I_ERI_F3y_S_Px_Py_bd;
  Double I_ERI_F2yz_Px_Px_Py_bd = I_ERI_Gx2yz_S_Px_Py_bd+ABX*I_ERI_F2yz_S_Px_Py_bd;
  Double I_ERI_Fy2z_Px_Px_Py_bd = I_ERI_Gxy2z_S_Px_Py_bd+ABX*I_ERI_Fy2z_S_Px_Py_bd;
  Double I_ERI_F3z_Px_Px_Py_bd = I_ERI_Gx3z_S_Px_Py_bd+ABX*I_ERI_F3z_S_Px_Py_bd;
  Double I_ERI_F3x_Py_Px_Py_bd = I_ERI_G3xy_S_Px_Py_bd+ABY*I_ERI_F3x_S_Px_Py_bd;
  Double I_ERI_F2xy_Py_Px_Py_bd = I_ERI_G2x2y_S_Px_Py_bd+ABY*I_ERI_F2xy_S_Px_Py_bd;
  Double I_ERI_F2xz_Py_Px_Py_bd = I_ERI_G2xyz_S_Px_Py_bd+ABY*I_ERI_F2xz_S_Px_Py_bd;
  Double I_ERI_Fx2y_Py_Px_Py_bd = I_ERI_Gx3y_S_Px_Py_bd+ABY*I_ERI_Fx2y_S_Px_Py_bd;
  Double I_ERI_Fxyz_Py_Px_Py_bd = I_ERI_Gx2yz_S_Px_Py_bd+ABY*I_ERI_Fxyz_S_Px_Py_bd;
  Double I_ERI_Fx2z_Py_Px_Py_bd = I_ERI_Gxy2z_S_Px_Py_bd+ABY*I_ERI_Fx2z_S_Px_Py_bd;
  Double I_ERI_F3y_Py_Px_Py_bd = I_ERI_G4y_S_Px_Py_bd+ABY*I_ERI_F3y_S_Px_Py_bd;
  Double I_ERI_F2yz_Py_Px_Py_bd = I_ERI_G3yz_S_Px_Py_bd+ABY*I_ERI_F2yz_S_Px_Py_bd;
  Double I_ERI_Fy2z_Py_Px_Py_bd = I_ERI_G2y2z_S_Px_Py_bd+ABY*I_ERI_Fy2z_S_Px_Py_bd;
  Double I_ERI_F3z_Py_Px_Py_bd = I_ERI_Gy3z_S_Px_Py_bd+ABY*I_ERI_F3z_S_Px_Py_bd;
  Double I_ERI_F3x_Pz_Px_Py_bd = I_ERI_G3xz_S_Px_Py_bd+ABZ*I_ERI_F3x_S_Px_Py_bd;
  Double I_ERI_F2xy_Pz_Px_Py_bd = I_ERI_G2xyz_S_Px_Py_bd+ABZ*I_ERI_F2xy_S_Px_Py_bd;
  Double I_ERI_F2xz_Pz_Px_Py_bd = I_ERI_G2x2z_S_Px_Py_bd+ABZ*I_ERI_F2xz_S_Px_Py_bd;
  Double I_ERI_Fx2y_Pz_Px_Py_bd = I_ERI_Gx2yz_S_Px_Py_bd+ABZ*I_ERI_Fx2y_S_Px_Py_bd;
  Double I_ERI_Fxyz_Pz_Px_Py_bd = I_ERI_Gxy2z_S_Px_Py_bd+ABZ*I_ERI_Fxyz_S_Px_Py_bd;
  Double I_ERI_Fx2z_Pz_Px_Py_bd = I_ERI_Gx3z_S_Px_Py_bd+ABZ*I_ERI_Fx2z_S_Px_Py_bd;
  Double I_ERI_F3y_Pz_Px_Py_bd = I_ERI_G3yz_S_Px_Py_bd+ABZ*I_ERI_F3y_S_Px_Py_bd;
  Double I_ERI_F2yz_Pz_Px_Py_bd = I_ERI_G2y2z_S_Px_Py_bd+ABZ*I_ERI_F2yz_S_Px_Py_bd;
  Double I_ERI_Fy2z_Pz_Px_Py_bd = I_ERI_Gy3z_S_Px_Py_bd+ABZ*I_ERI_Fy2z_S_Px_Py_bd;
  Double I_ERI_F3z_Pz_Px_Py_bd = I_ERI_G4z_S_Px_Py_bd+ABZ*I_ERI_F3z_S_Px_Py_bd;
  Double I_ERI_F3x_Px_Py_Py_bd = I_ERI_G4x_S_Py_Py_bd+ABX*I_ERI_F3x_S_Py_Py_bd;
  Double I_ERI_F2xy_Px_Py_Py_bd = I_ERI_G3xy_S_Py_Py_bd+ABX*I_ERI_F2xy_S_Py_Py_bd;
  Double I_ERI_F2xz_Px_Py_Py_bd = I_ERI_G3xz_S_Py_Py_bd+ABX*I_ERI_F2xz_S_Py_Py_bd;
  Double I_ERI_Fx2y_Px_Py_Py_bd = I_ERI_G2x2y_S_Py_Py_bd+ABX*I_ERI_Fx2y_S_Py_Py_bd;
  Double I_ERI_Fxyz_Px_Py_Py_bd = I_ERI_G2xyz_S_Py_Py_bd+ABX*I_ERI_Fxyz_S_Py_Py_bd;
  Double I_ERI_Fx2z_Px_Py_Py_bd = I_ERI_G2x2z_S_Py_Py_bd+ABX*I_ERI_Fx2z_S_Py_Py_bd;
  Double I_ERI_F3y_Px_Py_Py_bd = I_ERI_Gx3y_S_Py_Py_bd+ABX*I_ERI_F3y_S_Py_Py_bd;
  Double I_ERI_F2yz_Px_Py_Py_bd = I_ERI_Gx2yz_S_Py_Py_bd+ABX*I_ERI_F2yz_S_Py_Py_bd;
  Double I_ERI_Fy2z_Px_Py_Py_bd = I_ERI_Gxy2z_S_Py_Py_bd+ABX*I_ERI_Fy2z_S_Py_Py_bd;
  Double I_ERI_F3z_Px_Py_Py_bd = I_ERI_Gx3z_S_Py_Py_bd+ABX*I_ERI_F3z_S_Py_Py_bd;
  Double I_ERI_F3x_Py_Py_Py_bd = I_ERI_G3xy_S_Py_Py_bd+ABY*I_ERI_F3x_S_Py_Py_bd;
  Double I_ERI_F2xy_Py_Py_Py_bd = I_ERI_G2x2y_S_Py_Py_bd+ABY*I_ERI_F2xy_S_Py_Py_bd;
  Double I_ERI_F2xz_Py_Py_Py_bd = I_ERI_G2xyz_S_Py_Py_bd+ABY*I_ERI_F2xz_S_Py_Py_bd;
  Double I_ERI_Fx2y_Py_Py_Py_bd = I_ERI_Gx3y_S_Py_Py_bd+ABY*I_ERI_Fx2y_S_Py_Py_bd;
  Double I_ERI_Fxyz_Py_Py_Py_bd = I_ERI_Gx2yz_S_Py_Py_bd+ABY*I_ERI_Fxyz_S_Py_Py_bd;
  Double I_ERI_Fx2z_Py_Py_Py_bd = I_ERI_Gxy2z_S_Py_Py_bd+ABY*I_ERI_Fx2z_S_Py_Py_bd;
  Double I_ERI_F3y_Py_Py_Py_bd = I_ERI_G4y_S_Py_Py_bd+ABY*I_ERI_F3y_S_Py_Py_bd;
  Double I_ERI_F2yz_Py_Py_Py_bd = I_ERI_G3yz_S_Py_Py_bd+ABY*I_ERI_F2yz_S_Py_Py_bd;
  Double I_ERI_Fy2z_Py_Py_Py_bd = I_ERI_G2y2z_S_Py_Py_bd+ABY*I_ERI_Fy2z_S_Py_Py_bd;
  Double I_ERI_F3z_Py_Py_Py_bd = I_ERI_Gy3z_S_Py_Py_bd+ABY*I_ERI_F3z_S_Py_Py_bd;
  Double I_ERI_F3x_Pz_Py_Py_bd = I_ERI_G3xz_S_Py_Py_bd+ABZ*I_ERI_F3x_S_Py_Py_bd;
  Double I_ERI_F2xy_Pz_Py_Py_bd = I_ERI_G2xyz_S_Py_Py_bd+ABZ*I_ERI_F2xy_S_Py_Py_bd;
  Double I_ERI_F2xz_Pz_Py_Py_bd = I_ERI_G2x2z_S_Py_Py_bd+ABZ*I_ERI_F2xz_S_Py_Py_bd;
  Double I_ERI_Fx2y_Pz_Py_Py_bd = I_ERI_Gx2yz_S_Py_Py_bd+ABZ*I_ERI_Fx2y_S_Py_Py_bd;
  Double I_ERI_Fxyz_Pz_Py_Py_bd = I_ERI_Gxy2z_S_Py_Py_bd+ABZ*I_ERI_Fxyz_S_Py_Py_bd;
  Double I_ERI_Fx2z_Pz_Py_Py_bd = I_ERI_Gx3z_S_Py_Py_bd+ABZ*I_ERI_Fx2z_S_Py_Py_bd;
  Double I_ERI_F3y_Pz_Py_Py_bd = I_ERI_G3yz_S_Py_Py_bd+ABZ*I_ERI_F3y_S_Py_Py_bd;
  Double I_ERI_F2yz_Pz_Py_Py_bd = I_ERI_G2y2z_S_Py_Py_bd+ABZ*I_ERI_F2yz_S_Py_Py_bd;
  Double I_ERI_Fy2z_Pz_Py_Py_bd = I_ERI_Gy3z_S_Py_Py_bd+ABZ*I_ERI_Fy2z_S_Py_Py_bd;
  Double I_ERI_F3z_Pz_Py_Py_bd = I_ERI_G4z_S_Py_Py_bd+ABZ*I_ERI_F3z_S_Py_Py_bd;
  Double I_ERI_F3x_Px_Pz_Py_bd = I_ERI_G4x_S_Pz_Py_bd+ABX*I_ERI_F3x_S_Pz_Py_bd;
  Double I_ERI_F2xy_Px_Pz_Py_bd = I_ERI_G3xy_S_Pz_Py_bd+ABX*I_ERI_F2xy_S_Pz_Py_bd;
  Double I_ERI_F2xz_Px_Pz_Py_bd = I_ERI_G3xz_S_Pz_Py_bd+ABX*I_ERI_F2xz_S_Pz_Py_bd;
  Double I_ERI_Fx2y_Px_Pz_Py_bd = I_ERI_G2x2y_S_Pz_Py_bd+ABX*I_ERI_Fx2y_S_Pz_Py_bd;
  Double I_ERI_Fxyz_Px_Pz_Py_bd = I_ERI_G2xyz_S_Pz_Py_bd+ABX*I_ERI_Fxyz_S_Pz_Py_bd;
  Double I_ERI_Fx2z_Px_Pz_Py_bd = I_ERI_G2x2z_S_Pz_Py_bd+ABX*I_ERI_Fx2z_S_Pz_Py_bd;
  Double I_ERI_F3y_Px_Pz_Py_bd = I_ERI_Gx3y_S_Pz_Py_bd+ABX*I_ERI_F3y_S_Pz_Py_bd;
  Double I_ERI_F2yz_Px_Pz_Py_bd = I_ERI_Gx2yz_S_Pz_Py_bd+ABX*I_ERI_F2yz_S_Pz_Py_bd;
  Double I_ERI_Fy2z_Px_Pz_Py_bd = I_ERI_Gxy2z_S_Pz_Py_bd+ABX*I_ERI_Fy2z_S_Pz_Py_bd;
  Double I_ERI_F3z_Px_Pz_Py_bd = I_ERI_Gx3z_S_Pz_Py_bd+ABX*I_ERI_F3z_S_Pz_Py_bd;
  Double I_ERI_F3x_Py_Pz_Py_bd = I_ERI_G3xy_S_Pz_Py_bd+ABY*I_ERI_F3x_S_Pz_Py_bd;
  Double I_ERI_F2xy_Py_Pz_Py_bd = I_ERI_G2x2y_S_Pz_Py_bd+ABY*I_ERI_F2xy_S_Pz_Py_bd;
  Double I_ERI_F2xz_Py_Pz_Py_bd = I_ERI_G2xyz_S_Pz_Py_bd+ABY*I_ERI_F2xz_S_Pz_Py_bd;
  Double I_ERI_Fx2y_Py_Pz_Py_bd = I_ERI_Gx3y_S_Pz_Py_bd+ABY*I_ERI_Fx2y_S_Pz_Py_bd;
  Double I_ERI_Fxyz_Py_Pz_Py_bd = I_ERI_Gx2yz_S_Pz_Py_bd+ABY*I_ERI_Fxyz_S_Pz_Py_bd;
  Double I_ERI_Fx2z_Py_Pz_Py_bd = I_ERI_Gxy2z_S_Pz_Py_bd+ABY*I_ERI_Fx2z_S_Pz_Py_bd;
  Double I_ERI_F3y_Py_Pz_Py_bd = I_ERI_G4y_S_Pz_Py_bd+ABY*I_ERI_F3y_S_Pz_Py_bd;
  Double I_ERI_F2yz_Py_Pz_Py_bd = I_ERI_G3yz_S_Pz_Py_bd+ABY*I_ERI_F2yz_S_Pz_Py_bd;
  Double I_ERI_Fy2z_Py_Pz_Py_bd = I_ERI_G2y2z_S_Pz_Py_bd+ABY*I_ERI_Fy2z_S_Pz_Py_bd;
  Double I_ERI_F3z_Py_Pz_Py_bd = I_ERI_Gy3z_S_Pz_Py_bd+ABY*I_ERI_F3z_S_Pz_Py_bd;
  Double I_ERI_F3x_Pz_Pz_Py_bd = I_ERI_G3xz_S_Pz_Py_bd+ABZ*I_ERI_F3x_S_Pz_Py_bd;
  Double I_ERI_F2xy_Pz_Pz_Py_bd = I_ERI_G2xyz_S_Pz_Py_bd+ABZ*I_ERI_F2xy_S_Pz_Py_bd;
  Double I_ERI_F2xz_Pz_Pz_Py_bd = I_ERI_G2x2z_S_Pz_Py_bd+ABZ*I_ERI_F2xz_S_Pz_Py_bd;
  Double I_ERI_Fx2y_Pz_Pz_Py_bd = I_ERI_Gx2yz_S_Pz_Py_bd+ABZ*I_ERI_Fx2y_S_Pz_Py_bd;
  Double I_ERI_Fxyz_Pz_Pz_Py_bd = I_ERI_Gxy2z_S_Pz_Py_bd+ABZ*I_ERI_Fxyz_S_Pz_Py_bd;
  Double I_ERI_Fx2z_Pz_Pz_Py_bd = I_ERI_Gx3z_S_Pz_Py_bd+ABZ*I_ERI_Fx2z_S_Pz_Py_bd;
  Double I_ERI_F3y_Pz_Pz_Py_bd = I_ERI_G3yz_S_Pz_Py_bd+ABZ*I_ERI_F3y_S_Pz_Py_bd;
  Double I_ERI_F2yz_Pz_Pz_Py_bd = I_ERI_G2y2z_S_Pz_Py_bd+ABZ*I_ERI_F2yz_S_Pz_Py_bd;
  Double I_ERI_Fy2z_Pz_Pz_Py_bd = I_ERI_Gy3z_S_Pz_Py_bd+ABZ*I_ERI_Fy2z_S_Pz_Py_bd;
  Double I_ERI_F3z_Pz_Pz_Py_bd = I_ERI_G4z_S_Pz_Py_bd+ABZ*I_ERI_F3z_S_Pz_Py_bd;
  Double I_ERI_F3x_Px_Px_Pz_bd = I_ERI_G4x_S_Px_Pz_bd+ABX*I_ERI_F3x_S_Px_Pz_bd;
  Double I_ERI_F2xy_Px_Px_Pz_bd = I_ERI_G3xy_S_Px_Pz_bd+ABX*I_ERI_F2xy_S_Px_Pz_bd;
  Double I_ERI_F2xz_Px_Px_Pz_bd = I_ERI_G3xz_S_Px_Pz_bd+ABX*I_ERI_F2xz_S_Px_Pz_bd;
  Double I_ERI_Fx2y_Px_Px_Pz_bd = I_ERI_G2x2y_S_Px_Pz_bd+ABX*I_ERI_Fx2y_S_Px_Pz_bd;
  Double I_ERI_Fxyz_Px_Px_Pz_bd = I_ERI_G2xyz_S_Px_Pz_bd+ABX*I_ERI_Fxyz_S_Px_Pz_bd;
  Double I_ERI_Fx2z_Px_Px_Pz_bd = I_ERI_G2x2z_S_Px_Pz_bd+ABX*I_ERI_Fx2z_S_Px_Pz_bd;
  Double I_ERI_F3y_Px_Px_Pz_bd = I_ERI_Gx3y_S_Px_Pz_bd+ABX*I_ERI_F3y_S_Px_Pz_bd;
  Double I_ERI_F2yz_Px_Px_Pz_bd = I_ERI_Gx2yz_S_Px_Pz_bd+ABX*I_ERI_F2yz_S_Px_Pz_bd;
  Double I_ERI_Fy2z_Px_Px_Pz_bd = I_ERI_Gxy2z_S_Px_Pz_bd+ABX*I_ERI_Fy2z_S_Px_Pz_bd;
  Double I_ERI_F3z_Px_Px_Pz_bd = I_ERI_Gx3z_S_Px_Pz_bd+ABX*I_ERI_F3z_S_Px_Pz_bd;
  Double I_ERI_F3x_Py_Px_Pz_bd = I_ERI_G3xy_S_Px_Pz_bd+ABY*I_ERI_F3x_S_Px_Pz_bd;
  Double I_ERI_F2xy_Py_Px_Pz_bd = I_ERI_G2x2y_S_Px_Pz_bd+ABY*I_ERI_F2xy_S_Px_Pz_bd;
  Double I_ERI_F2xz_Py_Px_Pz_bd = I_ERI_G2xyz_S_Px_Pz_bd+ABY*I_ERI_F2xz_S_Px_Pz_bd;
  Double I_ERI_Fx2y_Py_Px_Pz_bd = I_ERI_Gx3y_S_Px_Pz_bd+ABY*I_ERI_Fx2y_S_Px_Pz_bd;
  Double I_ERI_Fxyz_Py_Px_Pz_bd = I_ERI_Gx2yz_S_Px_Pz_bd+ABY*I_ERI_Fxyz_S_Px_Pz_bd;
  Double I_ERI_Fx2z_Py_Px_Pz_bd = I_ERI_Gxy2z_S_Px_Pz_bd+ABY*I_ERI_Fx2z_S_Px_Pz_bd;
  Double I_ERI_F3y_Py_Px_Pz_bd = I_ERI_G4y_S_Px_Pz_bd+ABY*I_ERI_F3y_S_Px_Pz_bd;
  Double I_ERI_F2yz_Py_Px_Pz_bd = I_ERI_G3yz_S_Px_Pz_bd+ABY*I_ERI_F2yz_S_Px_Pz_bd;
  Double I_ERI_Fy2z_Py_Px_Pz_bd = I_ERI_G2y2z_S_Px_Pz_bd+ABY*I_ERI_Fy2z_S_Px_Pz_bd;
  Double I_ERI_F3z_Py_Px_Pz_bd = I_ERI_Gy3z_S_Px_Pz_bd+ABY*I_ERI_F3z_S_Px_Pz_bd;
  Double I_ERI_F3x_Pz_Px_Pz_bd = I_ERI_G3xz_S_Px_Pz_bd+ABZ*I_ERI_F3x_S_Px_Pz_bd;
  Double I_ERI_F2xy_Pz_Px_Pz_bd = I_ERI_G2xyz_S_Px_Pz_bd+ABZ*I_ERI_F2xy_S_Px_Pz_bd;
  Double I_ERI_F2xz_Pz_Px_Pz_bd = I_ERI_G2x2z_S_Px_Pz_bd+ABZ*I_ERI_F2xz_S_Px_Pz_bd;
  Double I_ERI_Fx2y_Pz_Px_Pz_bd = I_ERI_Gx2yz_S_Px_Pz_bd+ABZ*I_ERI_Fx2y_S_Px_Pz_bd;
  Double I_ERI_Fxyz_Pz_Px_Pz_bd = I_ERI_Gxy2z_S_Px_Pz_bd+ABZ*I_ERI_Fxyz_S_Px_Pz_bd;
  Double I_ERI_Fx2z_Pz_Px_Pz_bd = I_ERI_Gx3z_S_Px_Pz_bd+ABZ*I_ERI_Fx2z_S_Px_Pz_bd;
  Double I_ERI_F3y_Pz_Px_Pz_bd = I_ERI_G3yz_S_Px_Pz_bd+ABZ*I_ERI_F3y_S_Px_Pz_bd;
  Double I_ERI_F2yz_Pz_Px_Pz_bd = I_ERI_G2y2z_S_Px_Pz_bd+ABZ*I_ERI_F2yz_S_Px_Pz_bd;
  Double I_ERI_Fy2z_Pz_Px_Pz_bd = I_ERI_Gy3z_S_Px_Pz_bd+ABZ*I_ERI_Fy2z_S_Px_Pz_bd;
  Double I_ERI_F3z_Pz_Px_Pz_bd = I_ERI_G4z_S_Px_Pz_bd+ABZ*I_ERI_F3z_S_Px_Pz_bd;
  Double I_ERI_F3x_Px_Py_Pz_bd = I_ERI_G4x_S_Py_Pz_bd+ABX*I_ERI_F3x_S_Py_Pz_bd;
  Double I_ERI_F2xy_Px_Py_Pz_bd = I_ERI_G3xy_S_Py_Pz_bd+ABX*I_ERI_F2xy_S_Py_Pz_bd;
  Double I_ERI_F2xz_Px_Py_Pz_bd = I_ERI_G3xz_S_Py_Pz_bd+ABX*I_ERI_F2xz_S_Py_Pz_bd;
  Double I_ERI_Fx2y_Px_Py_Pz_bd = I_ERI_G2x2y_S_Py_Pz_bd+ABX*I_ERI_Fx2y_S_Py_Pz_bd;
  Double I_ERI_Fxyz_Px_Py_Pz_bd = I_ERI_G2xyz_S_Py_Pz_bd+ABX*I_ERI_Fxyz_S_Py_Pz_bd;
  Double I_ERI_Fx2z_Px_Py_Pz_bd = I_ERI_G2x2z_S_Py_Pz_bd+ABX*I_ERI_Fx2z_S_Py_Pz_bd;
  Double I_ERI_F3y_Px_Py_Pz_bd = I_ERI_Gx3y_S_Py_Pz_bd+ABX*I_ERI_F3y_S_Py_Pz_bd;
  Double I_ERI_F2yz_Px_Py_Pz_bd = I_ERI_Gx2yz_S_Py_Pz_bd+ABX*I_ERI_F2yz_S_Py_Pz_bd;
  Double I_ERI_Fy2z_Px_Py_Pz_bd = I_ERI_Gxy2z_S_Py_Pz_bd+ABX*I_ERI_Fy2z_S_Py_Pz_bd;
  Double I_ERI_F3z_Px_Py_Pz_bd = I_ERI_Gx3z_S_Py_Pz_bd+ABX*I_ERI_F3z_S_Py_Pz_bd;
  Double I_ERI_F3x_Py_Py_Pz_bd = I_ERI_G3xy_S_Py_Pz_bd+ABY*I_ERI_F3x_S_Py_Pz_bd;
  Double I_ERI_F2xy_Py_Py_Pz_bd = I_ERI_G2x2y_S_Py_Pz_bd+ABY*I_ERI_F2xy_S_Py_Pz_bd;
  Double I_ERI_F2xz_Py_Py_Pz_bd = I_ERI_G2xyz_S_Py_Pz_bd+ABY*I_ERI_F2xz_S_Py_Pz_bd;
  Double I_ERI_Fx2y_Py_Py_Pz_bd = I_ERI_Gx3y_S_Py_Pz_bd+ABY*I_ERI_Fx2y_S_Py_Pz_bd;
  Double I_ERI_Fxyz_Py_Py_Pz_bd = I_ERI_Gx2yz_S_Py_Pz_bd+ABY*I_ERI_Fxyz_S_Py_Pz_bd;
  Double I_ERI_Fx2z_Py_Py_Pz_bd = I_ERI_Gxy2z_S_Py_Pz_bd+ABY*I_ERI_Fx2z_S_Py_Pz_bd;
  Double I_ERI_F3y_Py_Py_Pz_bd = I_ERI_G4y_S_Py_Pz_bd+ABY*I_ERI_F3y_S_Py_Pz_bd;
  Double I_ERI_F2yz_Py_Py_Pz_bd = I_ERI_G3yz_S_Py_Pz_bd+ABY*I_ERI_F2yz_S_Py_Pz_bd;
  Double I_ERI_Fy2z_Py_Py_Pz_bd = I_ERI_G2y2z_S_Py_Pz_bd+ABY*I_ERI_Fy2z_S_Py_Pz_bd;
  Double I_ERI_F3z_Py_Py_Pz_bd = I_ERI_Gy3z_S_Py_Pz_bd+ABY*I_ERI_F3z_S_Py_Pz_bd;
  Double I_ERI_F3x_Pz_Py_Pz_bd = I_ERI_G3xz_S_Py_Pz_bd+ABZ*I_ERI_F3x_S_Py_Pz_bd;
  Double I_ERI_F2xy_Pz_Py_Pz_bd = I_ERI_G2xyz_S_Py_Pz_bd+ABZ*I_ERI_F2xy_S_Py_Pz_bd;
  Double I_ERI_F2xz_Pz_Py_Pz_bd = I_ERI_G2x2z_S_Py_Pz_bd+ABZ*I_ERI_F2xz_S_Py_Pz_bd;
  Double I_ERI_Fx2y_Pz_Py_Pz_bd = I_ERI_Gx2yz_S_Py_Pz_bd+ABZ*I_ERI_Fx2y_S_Py_Pz_bd;
  Double I_ERI_Fxyz_Pz_Py_Pz_bd = I_ERI_Gxy2z_S_Py_Pz_bd+ABZ*I_ERI_Fxyz_S_Py_Pz_bd;
  Double I_ERI_Fx2z_Pz_Py_Pz_bd = I_ERI_Gx3z_S_Py_Pz_bd+ABZ*I_ERI_Fx2z_S_Py_Pz_bd;
  Double I_ERI_F3y_Pz_Py_Pz_bd = I_ERI_G3yz_S_Py_Pz_bd+ABZ*I_ERI_F3y_S_Py_Pz_bd;
  Double I_ERI_F2yz_Pz_Py_Pz_bd = I_ERI_G2y2z_S_Py_Pz_bd+ABZ*I_ERI_F2yz_S_Py_Pz_bd;
  Double I_ERI_Fy2z_Pz_Py_Pz_bd = I_ERI_Gy3z_S_Py_Pz_bd+ABZ*I_ERI_Fy2z_S_Py_Pz_bd;
  Double I_ERI_F3z_Pz_Py_Pz_bd = I_ERI_G4z_S_Py_Pz_bd+ABZ*I_ERI_F3z_S_Py_Pz_bd;
  Double I_ERI_F3x_Px_Pz_Pz_bd = I_ERI_G4x_S_Pz_Pz_bd+ABX*I_ERI_F3x_S_Pz_Pz_bd;
  Double I_ERI_F2xy_Px_Pz_Pz_bd = I_ERI_G3xy_S_Pz_Pz_bd+ABX*I_ERI_F2xy_S_Pz_Pz_bd;
  Double I_ERI_F2xz_Px_Pz_Pz_bd = I_ERI_G3xz_S_Pz_Pz_bd+ABX*I_ERI_F2xz_S_Pz_Pz_bd;
  Double I_ERI_Fx2y_Px_Pz_Pz_bd = I_ERI_G2x2y_S_Pz_Pz_bd+ABX*I_ERI_Fx2y_S_Pz_Pz_bd;
  Double I_ERI_Fxyz_Px_Pz_Pz_bd = I_ERI_G2xyz_S_Pz_Pz_bd+ABX*I_ERI_Fxyz_S_Pz_Pz_bd;
  Double I_ERI_Fx2z_Px_Pz_Pz_bd = I_ERI_G2x2z_S_Pz_Pz_bd+ABX*I_ERI_Fx2z_S_Pz_Pz_bd;
  Double I_ERI_F3y_Px_Pz_Pz_bd = I_ERI_Gx3y_S_Pz_Pz_bd+ABX*I_ERI_F3y_S_Pz_Pz_bd;
  Double I_ERI_F2yz_Px_Pz_Pz_bd = I_ERI_Gx2yz_S_Pz_Pz_bd+ABX*I_ERI_F2yz_S_Pz_Pz_bd;
  Double I_ERI_Fy2z_Px_Pz_Pz_bd = I_ERI_Gxy2z_S_Pz_Pz_bd+ABX*I_ERI_Fy2z_S_Pz_Pz_bd;
  Double I_ERI_F3z_Px_Pz_Pz_bd = I_ERI_Gx3z_S_Pz_Pz_bd+ABX*I_ERI_F3z_S_Pz_Pz_bd;
  Double I_ERI_F3x_Py_Pz_Pz_bd = I_ERI_G3xy_S_Pz_Pz_bd+ABY*I_ERI_F3x_S_Pz_Pz_bd;
  Double I_ERI_F2xy_Py_Pz_Pz_bd = I_ERI_G2x2y_S_Pz_Pz_bd+ABY*I_ERI_F2xy_S_Pz_Pz_bd;
  Double I_ERI_F2xz_Py_Pz_Pz_bd = I_ERI_G2xyz_S_Pz_Pz_bd+ABY*I_ERI_F2xz_S_Pz_Pz_bd;
  Double I_ERI_Fx2y_Py_Pz_Pz_bd = I_ERI_Gx3y_S_Pz_Pz_bd+ABY*I_ERI_Fx2y_S_Pz_Pz_bd;
  Double I_ERI_Fxyz_Py_Pz_Pz_bd = I_ERI_Gx2yz_S_Pz_Pz_bd+ABY*I_ERI_Fxyz_S_Pz_Pz_bd;
  Double I_ERI_Fx2z_Py_Pz_Pz_bd = I_ERI_Gxy2z_S_Pz_Pz_bd+ABY*I_ERI_Fx2z_S_Pz_Pz_bd;
  Double I_ERI_F3y_Py_Pz_Pz_bd = I_ERI_G4y_S_Pz_Pz_bd+ABY*I_ERI_F3y_S_Pz_Pz_bd;
  Double I_ERI_F2yz_Py_Pz_Pz_bd = I_ERI_G3yz_S_Pz_Pz_bd+ABY*I_ERI_F2yz_S_Pz_Pz_bd;
  Double I_ERI_Fy2z_Py_Pz_Pz_bd = I_ERI_G2y2z_S_Pz_Pz_bd+ABY*I_ERI_Fy2z_S_Pz_Pz_bd;
  Double I_ERI_F3z_Py_Pz_Pz_bd = I_ERI_Gy3z_S_Pz_Pz_bd+ABY*I_ERI_F3z_S_Pz_Pz_bd;
  Double I_ERI_F3x_Pz_Pz_Pz_bd = I_ERI_G3xz_S_Pz_Pz_bd+ABZ*I_ERI_F3x_S_Pz_Pz_bd;
  Double I_ERI_F2xy_Pz_Pz_Pz_bd = I_ERI_G2xyz_S_Pz_Pz_bd+ABZ*I_ERI_F2xy_S_Pz_Pz_bd;
  Double I_ERI_F2xz_Pz_Pz_Pz_bd = I_ERI_G2x2z_S_Pz_Pz_bd+ABZ*I_ERI_F2xz_S_Pz_Pz_bd;
  Double I_ERI_Fx2y_Pz_Pz_Pz_bd = I_ERI_Gx2yz_S_Pz_Pz_bd+ABZ*I_ERI_Fx2y_S_Pz_Pz_bd;
  Double I_ERI_Fxyz_Pz_Pz_Pz_bd = I_ERI_Gxy2z_S_Pz_Pz_bd+ABZ*I_ERI_Fxyz_S_Pz_Pz_bd;
  Double I_ERI_Fx2z_Pz_Pz_Pz_bd = I_ERI_Gx3z_S_Pz_Pz_bd+ABZ*I_ERI_Fx2z_S_Pz_Pz_bd;
  Double I_ERI_F3y_Pz_Pz_Pz_bd = I_ERI_G3yz_S_Pz_Pz_bd+ABZ*I_ERI_F3y_S_Pz_Pz_bd;
  Double I_ERI_F2yz_Pz_Pz_Pz_bd = I_ERI_G2y2z_S_Pz_Pz_bd+ABZ*I_ERI_F2yz_S_Pz_Pz_bd;
  Double I_ERI_Fy2z_Pz_Pz_Pz_bd = I_ERI_Gy3z_S_Pz_Pz_bd+ABZ*I_ERI_Fy2z_S_Pz_Pz_bd;
  Double I_ERI_F3z_Pz_Pz_Pz_bd = I_ERI_G4z_S_Pz_Pz_bd+ABZ*I_ERI_F3z_S_Pz_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_F3x_D2x_Px_S_bb-2.0E0*1*I_ERI_F3x_S_Px_S_b;
  abcd[1] = 4.0E0*I_ERI_F2xy_D2x_Px_S_bb-2.0E0*1*I_ERI_F2xy_S_Px_S_b;
  abcd[2] = 4.0E0*I_ERI_F2xz_D2x_Px_S_bb-2.0E0*1*I_ERI_F2xz_S_Px_S_b;
  abcd[3] = 4.0E0*I_ERI_Fx2y_D2x_Px_S_bb-2.0E0*1*I_ERI_Fx2y_S_Px_S_b;
  abcd[4] = 4.0E0*I_ERI_Fxyz_D2x_Px_S_bb-2.0E0*1*I_ERI_Fxyz_S_Px_S_b;
  abcd[5] = 4.0E0*I_ERI_Fx2z_D2x_Px_S_bb-2.0E0*1*I_ERI_Fx2z_S_Px_S_b;
  abcd[6] = 4.0E0*I_ERI_F3y_D2x_Px_S_bb-2.0E0*1*I_ERI_F3y_S_Px_S_b;
  abcd[7] = 4.0E0*I_ERI_F2yz_D2x_Px_S_bb-2.0E0*1*I_ERI_F2yz_S_Px_S_b;
  abcd[8] = 4.0E0*I_ERI_Fy2z_D2x_Px_S_bb-2.0E0*1*I_ERI_Fy2z_S_Px_S_b;
  abcd[9] = 4.0E0*I_ERI_F3z_D2x_Px_S_bb-2.0E0*1*I_ERI_F3z_S_Px_S_b;
  abcd[10] = 4.0E0*I_ERI_F3x_D2x_Py_S_bb-2.0E0*1*I_ERI_F3x_S_Py_S_b;
  abcd[11] = 4.0E0*I_ERI_F2xy_D2x_Py_S_bb-2.0E0*1*I_ERI_F2xy_S_Py_S_b;
  abcd[12] = 4.0E0*I_ERI_F2xz_D2x_Py_S_bb-2.0E0*1*I_ERI_F2xz_S_Py_S_b;
  abcd[13] = 4.0E0*I_ERI_Fx2y_D2x_Py_S_bb-2.0E0*1*I_ERI_Fx2y_S_Py_S_b;
  abcd[14] = 4.0E0*I_ERI_Fxyz_D2x_Py_S_bb-2.0E0*1*I_ERI_Fxyz_S_Py_S_b;
  abcd[15] = 4.0E0*I_ERI_Fx2z_D2x_Py_S_bb-2.0E0*1*I_ERI_Fx2z_S_Py_S_b;
  abcd[16] = 4.0E0*I_ERI_F3y_D2x_Py_S_bb-2.0E0*1*I_ERI_F3y_S_Py_S_b;
  abcd[17] = 4.0E0*I_ERI_F2yz_D2x_Py_S_bb-2.0E0*1*I_ERI_F2yz_S_Py_S_b;
  abcd[18] = 4.0E0*I_ERI_Fy2z_D2x_Py_S_bb-2.0E0*1*I_ERI_Fy2z_S_Py_S_b;
  abcd[19] = 4.0E0*I_ERI_F3z_D2x_Py_S_bb-2.0E0*1*I_ERI_F3z_S_Py_S_b;
  abcd[20] = 4.0E0*I_ERI_F3x_D2x_Pz_S_bb-2.0E0*1*I_ERI_F3x_S_Pz_S_b;
  abcd[21] = 4.0E0*I_ERI_F2xy_D2x_Pz_S_bb-2.0E0*1*I_ERI_F2xy_S_Pz_S_b;
  abcd[22] = 4.0E0*I_ERI_F2xz_D2x_Pz_S_bb-2.0E0*1*I_ERI_F2xz_S_Pz_S_b;
  abcd[23] = 4.0E0*I_ERI_Fx2y_D2x_Pz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Pz_S_b;
  abcd[24] = 4.0E0*I_ERI_Fxyz_D2x_Pz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Pz_S_b;
  abcd[25] = 4.0E0*I_ERI_Fx2z_D2x_Pz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Pz_S_b;
  abcd[26] = 4.0E0*I_ERI_F3y_D2x_Pz_S_bb-2.0E0*1*I_ERI_F3y_S_Pz_S_b;
  abcd[27] = 4.0E0*I_ERI_F2yz_D2x_Pz_S_bb-2.0E0*1*I_ERI_F2yz_S_Pz_S_b;
  abcd[28] = 4.0E0*I_ERI_Fy2z_D2x_Pz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Pz_S_b;
  abcd[29] = 4.0E0*I_ERI_F3z_D2x_Pz_S_bb-2.0E0*1*I_ERI_F3z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[30] = 4.0E0*I_ERI_F3x_Dxy_Px_S_bb;
  abcd[31] = 4.0E0*I_ERI_F2xy_Dxy_Px_S_bb;
  abcd[32] = 4.0E0*I_ERI_F2xz_Dxy_Px_S_bb;
  abcd[33] = 4.0E0*I_ERI_Fx2y_Dxy_Px_S_bb;
  abcd[34] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_bb;
  abcd[35] = 4.0E0*I_ERI_Fx2z_Dxy_Px_S_bb;
  abcd[36] = 4.0E0*I_ERI_F3y_Dxy_Px_S_bb;
  abcd[37] = 4.0E0*I_ERI_F2yz_Dxy_Px_S_bb;
  abcd[38] = 4.0E0*I_ERI_Fy2z_Dxy_Px_S_bb;
  abcd[39] = 4.0E0*I_ERI_F3z_Dxy_Px_S_bb;
  abcd[40] = 4.0E0*I_ERI_F3x_Dxy_Py_S_bb;
  abcd[41] = 4.0E0*I_ERI_F2xy_Dxy_Py_S_bb;
  abcd[42] = 4.0E0*I_ERI_F2xz_Dxy_Py_S_bb;
  abcd[43] = 4.0E0*I_ERI_Fx2y_Dxy_Py_S_bb;
  abcd[44] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_bb;
  abcd[45] = 4.0E0*I_ERI_Fx2z_Dxy_Py_S_bb;
  abcd[46] = 4.0E0*I_ERI_F3y_Dxy_Py_S_bb;
  abcd[47] = 4.0E0*I_ERI_F2yz_Dxy_Py_S_bb;
  abcd[48] = 4.0E0*I_ERI_Fy2z_Dxy_Py_S_bb;
  abcd[49] = 4.0E0*I_ERI_F3z_Dxy_Py_S_bb;
  abcd[50] = 4.0E0*I_ERI_F3x_Dxy_Pz_S_bb;
  abcd[51] = 4.0E0*I_ERI_F2xy_Dxy_Pz_S_bb;
  abcd[52] = 4.0E0*I_ERI_F2xz_Dxy_Pz_S_bb;
  abcd[53] = 4.0E0*I_ERI_Fx2y_Dxy_Pz_S_bb;
  abcd[54] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_bb;
  abcd[55] = 4.0E0*I_ERI_Fx2z_Dxy_Pz_S_bb;
  abcd[56] = 4.0E0*I_ERI_F3y_Dxy_Pz_S_bb;
  abcd[57] = 4.0E0*I_ERI_F2yz_Dxy_Pz_S_bb;
  abcd[58] = 4.0E0*I_ERI_Fy2z_Dxy_Pz_S_bb;
  abcd[59] = 4.0E0*I_ERI_F3z_Dxy_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[60] = 4.0E0*I_ERI_F3x_Dxz_Px_S_bb;
  abcd[61] = 4.0E0*I_ERI_F2xy_Dxz_Px_S_bb;
  abcd[62] = 4.0E0*I_ERI_F2xz_Dxz_Px_S_bb;
  abcd[63] = 4.0E0*I_ERI_Fx2y_Dxz_Px_S_bb;
  abcd[64] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_bb;
  abcd[65] = 4.0E0*I_ERI_Fx2z_Dxz_Px_S_bb;
  abcd[66] = 4.0E0*I_ERI_F3y_Dxz_Px_S_bb;
  abcd[67] = 4.0E0*I_ERI_F2yz_Dxz_Px_S_bb;
  abcd[68] = 4.0E0*I_ERI_Fy2z_Dxz_Px_S_bb;
  abcd[69] = 4.0E0*I_ERI_F3z_Dxz_Px_S_bb;
  abcd[70] = 4.0E0*I_ERI_F3x_Dxz_Py_S_bb;
  abcd[71] = 4.0E0*I_ERI_F2xy_Dxz_Py_S_bb;
  abcd[72] = 4.0E0*I_ERI_F2xz_Dxz_Py_S_bb;
  abcd[73] = 4.0E0*I_ERI_Fx2y_Dxz_Py_S_bb;
  abcd[74] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_bb;
  abcd[75] = 4.0E0*I_ERI_Fx2z_Dxz_Py_S_bb;
  abcd[76] = 4.0E0*I_ERI_F3y_Dxz_Py_S_bb;
  abcd[77] = 4.0E0*I_ERI_F2yz_Dxz_Py_S_bb;
  abcd[78] = 4.0E0*I_ERI_Fy2z_Dxz_Py_S_bb;
  abcd[79] = 4.0E0*I_ERI_F3z_Dxz_Py_S_bb;
  abcd[80] = 4.0E0*I_ERI_F3x_Dxz_Pz_S_bb;
  abcd[81] = 4.0E0*I_ERI_F2xy_Dxz_Pz_S_bb;
  abcd[82] = 4.0E0*I_ERI_F2xz_Dxz_Pz_S_bb;
  abcd[83] = 4.0E0*I_ERI_Fx2y_Dxz_Pz_S_bb;
  abcd[84] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_bb;
  abcd[85] = 4.0E0*I_ERI_Fx2z_Dxz_Pz_S_bb;
  abcd[86] = 4.0E0*I_ERI_F3y_Dxz_Pz_S_bb;
  abcd[87] = 4.0E0*I_ERI_F2yz_Dxz_Pz_S_bb;
  abcd[88] = 4.0E0*I_ERI_Fy2z_Dxz_Pz_S_bb;
  abcd[89] = 4.0E0*I_ERI_F3z_Dxz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[90] = 4.0E0*I_ERI_F3x_D2y_Px_S_bb-2.0E0*1*I_ERI_F3x_S_Px_S_b;
  abcd[91] = 4.0E0*I_ERI_F2xy_D2y_Px_S_bb-2.0E0*1*I_ERI_F2xy_S_Px_S_b;
  abcd[92] = 4.0E0*I_ERI_F2xz_D2y_Px_S_bb-2.0E0*1*I_ERI_F2xz_S_Px_S_b;
  abcd[93] = 4.0E0*I_ERI_Fx2y_D2y_Px_S_bb-2.0E0*1*I_ERI_Fx2y_S_Px_S_b;
  abcd[94] = 4.0E0*I_ERI_Fxyz_D2y_Px_S_bb-2.0E0*1*I_ERI_Fxyz_S_Px_S_b;
  abcd[95] = 4.0E0*I_ERI_Fx2z_D2y_Px_S_bb-2.0E0*1*I_ERI_Fx2z_S_Px_S_b;
  abcd[96] = 4.0E0*I_ERI_F3y_D2y_Px_S_bb-2.0E0*1*I_ERI_F3y_S_Px_S_b;
  abcd[97] = 4.0E0*I_ERI_F2yz_D2y_Px_S_bb-2.0E0*1*I_ERI_F2yz_S_Px_S_b;
  abcd[98] = 4.0E0*I_ERI_Fy2z_D2y_Px_S_bb-2.0E0*1*I_ERI_Fy2z_S_Px_S_b;
  abcd[99] = 4.0E0*I_ERI_F3z_D2y_Px_S_bb-2.0E0*1*I_ERI_F3z_S_Px_S_b;
  abcd[100] = 4.0E0*I_ERI_F3x_D2y_Py_S_bb-2.0E0*1*I_ERI_F3x_S_Py_S_b;
  abcd[101] = 4.0E0*I_ERI_F2xy_D2y_Py_S_bb-2.0E0*1*I_ERI_F2xy_S_Py_S_b;
  abcd[102] = 4.0E0*I_ERI_F2xz_D2y_Py_S_bb-2.0E0*1*I_ERI_F2xz_S_Py_S_b;
  abcd[103] = 4.0E0*I_ERI_Fx2y_D2y_Py_S_bb-2.0E0*1*I_ERI_Fx2y_S_Py_S_b;
  abcd[104] = 4.0E0*I_ERI_Fxyz_D2y_Py_S_bb-2.0E0*1*I_ERI_Fxyz_S_Py_S_b;
  abcd[105] = 4.0E0*I_ERI_Fx2z_D2y_Py_S_bb-2.0E0*1*I_ERI_Fx2z_S_Py_S_b;
  abcd[106] = 4.0E0*I_ERI_F3y_D2y_Py_S_bb-2.0E0*1*I_ERI_F3y_S_Py_S_b;
  abcd[107] = 4.0E0*I_ERI_F2yz_D2y_Py_S_bb-2.0E0*1*I_ERI_F2yz_S_Py_S_b;
  abcd[108] = 4.0E0*I_ERI_Fy2z_D2y_Py_S_bb-2.0E0*1*I_ERI_Fy2z_S_Py_S_b;
  abcd[109] = 4.0E0*I_ERI_F3z_D2y_Py_S_bb-2.0E0*1*I_ERI_F3z_S_Py_S_b;
  abcd[110] = 4.0E0*I_ERI_F3x_D2y_Pz_S_bb-2.0E0*1*I_ERI_F3x_S_Pz_S_b;
  abcd[111] = 4.0E0*I_ERI_F2xy_D2y_Pz_S_bb-2.0E0*1*I_ERI_F2xy_S_Pz_S_b;
  abcd[112] = 4.0E0*I_ERI_F2xz_D2y_Pz_S_bb-2.0E0*1*I_ERI_F2xz_S_Pz_S_b;
  abcd[113] = 4.0E0*I_ERI_Fx2y_D2y_Pz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Pz_S_b;
  abcd[114] = 4.0E0*I_ERI_Fxyz_D2y_Pz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Pz_S_b;
  abcd[115] = 4.0E0*I_ERI_Fx2z_D2y_Pz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Pz_S_b;
  abcd[116] = 4.0E0*I_ERI_F3y_D2y_Pz_S_bb-2.0E0*1*I_ERI_F3y_S_Pz_S_b;
  abcd[117] = 4.0E0*I_ERI_F2yz_D2y_Pz_S_bb-2.0E0*1*I_ERI_F2yz_S_Pz_S_b;
  abcd[118] = 4.0E0*I_ERI_Fy2z_D2y_Pz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Pz_S_b;
  abcd[119] = 4.0E0*I_ERI_F3z_D2y_Pz_S_bb-2.0E0*1*I_ERI_F3z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[120] = 4.0E0*I_ERI_F3x_Dyz_Px_S_bb;
  abcd[121] = 4.0E0*I_ERI_F2xy_Dyz_Px_S_bb;
  abcd[122] = 4.0E0*I_ERI_F2xz_Dyz_Px_S_bb;
  abcd[123] = 4.0E0*I_ERI_Fx2y_Dyz_Px_S_bb;
  abcd[124] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_bb;
  abcd[125] = 4.0E0*I_ERI_Fx2z_Dyz_Px_S_bb;
  abcd[126] = 4.0E0*I_ERI_F3y_Dyz_Px_S_bb;
  abcd[127] = 4.0E0*I_ERI_F2yz_Dyz_Px_S_bb;
  abcd[128] = 4.0E0*I_ERI_Fy2z_Dyz_Px_S_bb;
  abcd[129] = 4.0E0*I_ERI_F3z_Dyz_Px_S_bb;
  abcd[130] = 4.0E0*I_ERI_F3x_Dyz_Py_S_bb;
  abcd[131] = 4.0E0*I_ERI_F2xy_Dyz_Py_S_bb;
  abcd[132] = 4.0E0*I_ERI_F2xz_Dyz_Py_S_bb;
  abcd[133] = 4.0E0*I_ERI_Fx2y_Dyz_Py_S_bb;
  abcd[134] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_bb;
  abcd[135] = 4.0E0*I_ERI_Fx2z_Dyz_Py_S_bb;
  abcd[136] = 4.0E0*I_ERI_F3y_Dyz_Py_S_bb;
  abcd[137] = 4.0E0*I_ERI_F2yz_Dyz_Py_S_bb;
  abcd[138] = 4.0E0*I_ERI_Fy2z_Dyz_Py_S_bb;
  abcd[139] = 4.0E0*I_ERI_F3z_Dyz_Py_S_bb;
  abcd[140] = 4.0E0*I_ERI_F3x_Dyz_Pz_S_bb;
  abcd[141] = 4.0E0*I_ERI_F2xy_Dyz_Pz_S_bb;
  abcd[142] = 4.0E0*I_ERI_F2xz_Dyz_Pz_S_bb;
  abcd[143] = 4.0E0*I_ERI_Fx2y_Dyz_Pz_S_bb;
  abcd[144] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_bb;
  abcd[145] = 4.0E0*I_ERI_Fx2z_Dyz_Pz_S_bb;
  abcd[146] = 4.0E0*I_ERI_F3y_Dyz_Pz_S_bb;
  abcd[147] = 4.0E0*I_ERI_F2yz_Dyz_Pz_S_bb;
  abcd[148] = 4.0E0*I_ERI_Fy2z_Dyz_Pz_S_bb;
  abcd[149] = 4.0E0*I_ERI_F3z_Dyz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  abcd[150] = 4.0E0*I_ERI_F3x_D2z_Px_S_bb-2.0E0*1*I_ERI_F3x_S_Px_S_b;
  abcd[151] = 4.0E0*I_ERI_F2xy_D2z_Px_S_bb-2.0E0*1*I_ERI_F2xy_S_Px_S_b;
  abcd[152] = 4.0E0*I_ERI_F2xz_D2z_Px_S_bb-2.0E0*1*I_ERI_F2xz_S_Px_S_b;
  abcd[153] = 4.0E0*I_ERI_Fx2y_D2z_Px_S_bb-2.0E0*1*I_ERI_Fx2y_S_Px_S_b;
  abcd[154] = 4.0E0*I_ERI_Fxyz_D2z_Px_S_bb-2.0E0*1*I_ERI_Fxyz_S_Px_S_b;
  abcd[155] = 4.0E0*I_ERI_Fx2z_D2z_Px_S_bb-2.0E0*1*I_ERI_Fx2z_S_Px_S_b;
  abcd[156] = 4.0E0*I_ERI_F3y_D2z_Px_S_bb-2.0E0*1*I_ERI_F3y_S_Px_S_b;
  abcd[157] = 4.0E0*I_ERI_F2yz_D2z_Px_S_bb-2.0E0*1*I_ERI_F2yz_S_Px_S_b;
  abcd[158] = 4.0E0*I_ERI_Fy2z_D2z_Px_S_bb-2.0E0*1*I_ERI_Fy2z_S_Px_S_b;
  abcd[159] = 4.0E0*I_ERI_F3z_D2z_Px_S_bb-2.0E0*1*I_ERI_F3z_S_Px_S_b;
  abcd[160] = 4.0E0*I_ERI_F3x_D2z_Py_S_bb-2.0E0*1*I_ERI_F3x_S_Py_S_b;
  abcd[161] = 4.0E0*I_ERI_F2xy_D2z_Py_S_bb-2.0E0*1*I_ERI_F2xy_S_Py_S_b;
  abcd[162] = 4.0E0*I_ERI_F2xz_D2z_Py_S_bb-2.0E0*1*I_ERI_F2xz_S_Py_S_b;
  abcd[163] = 4.0E0*I_ERI_Fx2y_D2z_Py_S_bb-2.0E0*1*I_ERI_Fx2y_S_Py_S_b;
  abcd[164] = 4.0E0*I_ERI_Fxyz_D2z_Py_S_bb-2.0E0*1*I_ERI_Fxyz_S_Py_S_b;
  abcd[165] = 4.0E0*I_ERI_Fx2z_D2z_Py_S_bb-2.0E0*1*I_ERI_Fx2z_S_Py_S_b;
  abcd[166] = 4.0E0*I_ERI_F3y_D2z_Py_S_bb-2.0E0*1*I_ERI_F3y_S_Py_S_b;
  abcd[167] = 4.0E0*I_ERI_F2yz_D2z_Py_S_bb-2.0E0*1*I_ERI_F2yz_S_Py_S_b;
  abcd[168] = 4.0E0*I_ERI_Fy2z_D2z_Py_S_bb-2.0E0*1*I_ERI_Fy2z_S_Py_S_b;
  abcd[169] = 4.0E0*I_ERI_F3z_D2z_Py_S_bb-2.0E0*1*I_ERI_F3z_S_Py_S_b;
  abcd[170] = 4.0E0*I_ERI_F3x_D2z_Pz_S_bb-2.0E0*1*I_ERI_F3x_S_Pz_S_b;
  abcd[171] = 4.0E0*I_ERI_F2xy_D2z_Pz_S_bb-2.0E0*1*I_ERI_F2xy_S_Pz_S_b;
  abcd[172] = 4.0E0*I_ERI_F2xz_D2z_Pz_S_bb-2.0E0*1*I_ERI_F2xz_S_Pz_S_b;
  abcd[173] = 4.0E0*I_ERI_Fx2y_D2z_Pz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Pz_S_b;
  abcd[174] = 4.0E0*I_ERI_Fxyz_D2z_Pz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Pz_S_b;
  abcd[175] = 4.0E0*I_ERI_Fx2z_D2z_Pz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Pz_S_b;
  abcd[176] = 4.0E0*I_ERI_F3y_D2z_Pz_S_bb-2.0E0*1*I_ERI_F3y_S_Pz_S_b;
  abcd[177] = 4.0E0*I_ERI_F2yz_D2z_Pz_S_bb-2.0E0*1*I_ERI_F2yz_S_Pz_S_b;
  abcd[178] = 4.0E0*I_ERI_Fy2z_D2z_Pz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Pz_S_b;
  abcd[179] = 4.0E0*I_ERI_F3z_D2z_Pz_S_bb-2.0E0*1*I_ERI_F3z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_F3x_Px_D2x_S_bc-2.0E0*1*I_ERI_F3x_Px_S_S_b;
  abcd[181] = 4.0E0*I_ERI_F2xy_Px_D2x_S_bc-2.0E0*1*I_ERI_F2xy_Px_S_S_b;
  abcd[182] = 4.0E0*I_ERI_F2xz_Px_D2x_S_bc-2.0E0*1*I_ERI_F2xz_Px_S_S_b;
  abcd[183] = 4.0E0*I_ERI_Fx2y_Px_D2x_S_bc-2.0E0*1*I_ERI_Fx2y_Px_S_S_b;
  abcd[184] = 4.0E0*I_ERI_Fxyz_Px_D2x_S_bc-2.0E0*1*I_ERI_Fxyz_Px_S_S_b;
  abcd[185] = 4.0E0*I_ERI_Fx2z_Px_D2x_S_bc-2.0E0*1*I_ERI_Fx2z_Px_S_S_b;
  abcd[186] = 4.0E0*I_ERI_F3y_Px_D2x_S_bc-2.0E0*1*I_ERI_F3y_Px_S_S_b;
  abcd[187] = 4.0E0*I_ERI_F2yz_Px_D2x_S_bc-2.0E0*1*I_ERI_F2yz_Px_S_S_b;
  abcd[188] = 4.0E0*I_ERI_Fy2z_Px_D2x_S_bc-2.0E0*1*I_ERI_Fy2z_Px_S_S_b;
  abcd[189] = 4.0E0*I_ERI_F3z_Px_D2x_S_bc-2.0E0*1*I_ERI_F3z_Px_S_S_b;
  abcd[190] = 4.0E0*I_ERI_F3x_Px_Dxy_S_bc;
  abcd[191] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_bc;
  abcd[192] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_bc;
  abcd[193] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_bc;
  abcd[194] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_bc;
  abcd[195] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_bc;
  abcd[196] = 4.0E0*I_ERI_F3y_Px_Dxy_S_bc;
  abcd[197] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_bc;
  abcd[198] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_bc;
  abcd[199] = 4.0E0*I_ERI_F3z_Px_Dxy_S_bc;
  abcd[200] = 4.0E0*I_ERI_F3x_Px_Dxz_S_bc;
  abcd[201] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_bc;
  abcd[202] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_bc;
  abcd[203] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_bc;
  abcd[204] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_bc;
  abcd[205] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_bc;
  abcd[206] = 4.0E0*I_ERI_F3y_Px_Dxz_S_bc;
  abcd[207] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_bc;
  abcd[208] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_bc;
  abcd[209] = 4.0E0*I_ERI_F3z_Px_Dxz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[210] = 4.0E0*I_ERI_F3x_Px_Dxy_S_bc;
  abcd[211] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_bc;
  abcd[212] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_bc;
  abcd[213] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_bc;
  abcd[214] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_bc;
  abcd[215] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_bc;
  abcd[216] = 4.0E0*I_ERI_F3y_Px_Dxy_S_bc;
  abcd[217] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_bc;
  abcd[218] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_bc;
  abcd[219] = 4.0E0*I_ERI_F3z_Px_Dxy_S_bc;
  abcd[220] = 4.0E0*I_ERI_F3x_Px_D2y_S_bc-2.0E0*1*I_ERI_F3x_Px_S_S_b;
  abcd[221] = 4.0E0*I_ERI_F2xy_Px_D2y_S_bc-2.0E0*1*I_ERI_F2xy_Px_S_S_b;
  abcd[222] = 4.0E0*I_ERI_F2xz_Px_D2y_S_bc-2.0E0*1*I_ERI_F2xz_Px_S_S_b;
  abcd[223] = 4.0E0*I_ERI_Fx2y_Px_D2y_S_bc-2.0E0*1*I_ERI_Fx2y_Px_S_S_b;
  abcd[224] = 4.0E0*I_ERI_Fxyz_Px_D2y_S_bc-2.0E0*1*I_ERI_Fxyz_Px_S_S_b;
  abcd[225] = 4.0E0*I_ERI_Fx2z_Px_D2y_S_bc-2.0E0*1*I_ERI_Fx2z_Px_S_S_b;
  abcd[226] = 4.0E0*I_ERI_F3y_Px_D2y_S_bc-2.0E0*1*I_ERI_F3y_Px_S_S_b;
  abcd[227] = 4.0E0*I_ERI_F2yz_Px_D2y_S_bc-2.0E0*1*I_ERI_F2yz_Px_S_S_b;
  abcd[228] = 4.0E0*I_ERI_Fy2z_Px_D2y_S_bc-2.0E0*1*I_ERI_Fy2z_Px_S_S_b;
  abcd[229] = 4.0E0*I_ERI_F3z_Px_D2y_S_bc-2.0E0*1*I_ERI_F3z_Px_S_S_b;
  abcd[230] = 4.0E0*I_ERI_F3x_Px_Dyz_S_bc;
  abcd[231] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_bc;
  abcd[232] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_bc;
  abcd[233] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_bc;
  abcd[234] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_bc;
  abcd[235] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_bc;
  abcd[236] = 4.0E0*I_ERI_F3y_Px_Dyz_S_bc;
  abcd[237] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_bc;
  abcd[238] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_bc;
  abcd[239] = 4.0E0*I_ERI_F3z_Px_Dyz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_F3x_Px_Dxz_S_bc;
  abcd[241] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_bc;
  abcd[242] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_bc;
  abcd[243] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_bc;
  abcd[244] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_bc;
  abcd[245] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_bc;
  abcd[246] = 4.0E0*I_ERI_F3y_Px_Dxz_S_bc;
  abcd[247] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_bc;
  abcd[248] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_bc;
  abcd[249] = 4.0E0*I_ERI_F3z_Px_Dxz_S_bc;
  abcd[250] = 4.0E0*I_ERI_F3x_Px_Dyz_S_bc;
  abcd[251] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_bc;
  abcd[252] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_bc;
  abcd[253] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_bc;
  abcd[254] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_bc;
  abcd[255] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_bc;
  abcd[256] = 4.0E0*I_ERI_F3y_Px_Dyz_S_bc;
  abcd[257] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_bc;
  abcd[258] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_bc;
  abcd[259] = 4.0E0*I_ERI_F3z_Px_Dyz_S_bc;
  abcd[260] = 4.0E0*I_ERI_F3x_Px_D2z_S_bc-2.0E0*1*I_ERI_F3x_Px_S_S_b;
  abcd[261] = 4.0E0*I_ERI_F2xy_Px_D2z_S_bc-2.0E0*1*I_ERI_F2xy_Px_S_S_b;
  abcd[262] = 4.0E0*I_ERI_F2xz_Px_D2z_S_bc-2.0E0*1*I_ERI_F2xz_Px_S_S_b;
  abcd[263] = 4.0E0*I_ERI_Fx2y_Px_D2z_S_bc-2.0E0*1*I_ERI_Fx2y_Px_S_S_b;
  abcd[264] = 4.0E0*I_ERI_Fxyz_Px_D2z_S_bc-2.0E0*1*I_ERI_Fxyz_Px_S_S_b;
  abcd[265] = 4.0E0*I_ERI_Fx2z_Px_D2z_S_bc-2.0E0*1*I_ERI_Fx2z_Px_S_S_b;
  abcd[266] = 4.0E0*I_ERI_F3y_Px_D2z_S_bc-2.0E0*1*I_ERI_F3y_Px_S_S_b;
  abcd[267] = 4.0E0*I_ERI_F2yz_Px_D2z_S_bc-2.0E0*1*I_ERI_F2yz_Px_S_S_b;
  abcd[268] = 4.0E0*I_ERI_Fy2z_Px_D2z_S_bc-2.0E0*1*I_ERI_Fy2z_Px_S_S_b;
  abcd[269] = 4.0E0*I_ERI_F3z_Px_D2z_S_bc-2.0E0*1*I_ERI_F3z_Px_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[270] = 4.0E0*I_ERI_F3x_Py_D2x_S_bc-2.0E0*1*I_ERI_F3x_Py_S_S_b;
  abcd[271] = 4.0E0*I_ERI_F2xy_Py_D2x_S_bc-2.0E0*1*I_ERI_F2xy_Py_S_S_b;
  abcd[272] = 4.0E0*I_ERI_F2xz_Py_D2x_S_bc-2.0E0*1*I_ERI_F2xz_Py_S_S_b;
  abcd[273] = 4.0E0*I_ERI_Fx2y_Py_D2x_S_bc-2.0E0*1*I_ERI_Fx2y_Py_S_S_b;
  abcd[274] = 4.0E0*I_ERI_Fxyz_Py_D2x_S_bc-2.0E0*1*I_ERI_Fxyz_Py_S_S_b;
  abcd[275] = 4.0E0*I_ERI_Fx2z_Py_D2x_S_bc-2.0E0*1*I_ERI_Fx2z_Py_S_S_b;
  abcd[276] = 4.0E0*I_ERI_F3y_Py_D2x_S_bc-2.0E0*1*I_ERI_F3y_Py_S_S_b;
  abcd[277] = 4.0E0*I_ERI_F2yz_Py_D2x_S_bc-2.0E0*1*I_ERI_F2yz_Py_S_S_b;
  abcd[278] = 4.0E0*I_ERI_Fy2z_Py_D2x_S_bc-2.0E0*1*I_ERI_Fy2z_Py_S_S_b;
  abcd[279] = 4.0E0*I_ERI_F3z_Py_D2x_S_bc-2.0E0*1*I_ERI_F3z_Py_S_S_b;
  abcd[280] = 4.0E0*I_ERI_F3x_Py_Dxy_S_bc;
  abcd[281] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_bc;
  abcd[282] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_bc;
  abcd[283] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_bc;
  abcd[284] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_bc;
  abcd[285] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_bc;
  abcd[286] = 4.0E0*I_ERI_F3y_Py_Dxy_S_bc;
  abcd[287] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_bc;
  abcd[288] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_bc;
  abcd[289] = 4.0E0*I_ERI_F3z_Py_Dxy_S_bc;
  abcd[290] = 4.0E0*I_ERI_F3x_Py_Dxz_S_bc;
  abcd[291] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_bc;
  abcd[292] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_bc;
  abcd[293] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_bc;
  abcd[294] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_bc;
  abcd[295] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_bc;
  abcd[296] = 4.0E0*I_ERI_F3y_Py_Dxz_S_bc;
  abcd[297] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_bc;
  abcd[298] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_bc;
  abcd[299] = 4.0E0*I_ERI_F3z_Py_Dxz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[300] = 4.0E0*I_ERI_F3x_Py_Dxy_S_bc;
  abcd[301] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_bc;
  abcd[302] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_bc;
  abcd[303] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_bc;
  abcd[304] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_bc;
  abcd[305] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_bc;
  abcd[306] = 4.0E0*I_ERI_F3y_Py_Dxy_S_bc;
  abcd[307] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_bc;
  abcd[308] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_bc;
  abcd[309] = 4.0E0*I_ERI_F3z_Py_Dxy_S_bc;
  abcd[310] = 4.0E0*I_ERI_F3x_Py_D2y_S_bc-2.0E0*1*I_ERI_F3x_Py_S_S_b;
  abcd[311] = 4.0E0*I_ERI_F2xy_Py_D2y_S_bc-2.0E0*1*I_ERI_F2xy_Py_S_S_b;
  abcd[312] = 4.0E0*I_ERI_F2xz_Py_D2y_S_bc-2.0E0*1*I_ERI_F2xz_Py_S_S_b;
  abcd[313] = 4.0E0*I_ERI_Fx2y_Py_D2y_S_bc-2.0E0*1*I_ERI_Fx2y_Py_S_S_b;
  abcd[314] = 4.0E0*I_ERI_Fxyz_Py_D2y_S_bc-2.0E0*1*I_ERI_Fxyz_Py_S_S_b;
  abcd[315] = 4.0E0*I_ERI_Fx2z_Py_D2y_S_bc-2.0E0*1*I_ERI_Fx2z_Py_S_S_b;
  abcd[316] = 4.0E0*I_ERI_F3y_Py_D2y_S_bc-2.0E0*1*I_ERI_F3y_Py_S_S_b;
  abcd[317] = 4.0E0*I_ERI_F2yz_Py_D2y_S_bc-2.0E0*1*I_ERI_F2yz_Py_S_S_b;
  abcd[318] = 4.0E0*I_ERI_Fy2z_Py_D2y_S_bc-2.0E0*1*I_ERI_Fy2z_Py_S_S_b;
  abcd[319] = 4.0E0*I_ERI_F3z_Py_D2y_S_bc-2.0E0*1*I_ERI_F3z_Py_S_S_b;
  abcd[320] = 4.0E0*I_ERI_F3x_Py_Dyz_S_bc;
  abcd[321] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_bc;
  abcd[322] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_bc;
  abcd[323] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_bc;
  abcd[324] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_bc;
  abcd[325] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_bc;
  abcd[326] = 4.0E0*I_ERI_F3y_Py_Dyz_S_bc;
  abcd[327] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_bc;
  abcd[328] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_bc;
  abcd[329] = 4.0E0*I_ERI_F3z_Py_Dyz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[330] = 4.0E0*I_ERI_F3x_Py_Dxz_S_bc;
  abcd[331] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_bc;
  abcd[332] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_bc;
  abcd[333] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_bc;
  abcd[334] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_bc;
  abcd[335] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_bc;
  abcd[336] = 4.0E0*I_ERI_F3y_Py_Dxz_S_bc;
  abcd[337] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_bc;
  abcd[338] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_bc;
  abcd[339] = 4.0E0*I_ERI_F3z_Py_Dxz_S_bc;
  abcd[340] = 4.0E0*I_ERI_F3x_Py_Dyz_S_bc;
  abcd[341] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_bc;
  abcd[342] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_bc;
  abcd[343] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_bc;
  abcd[344] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_bc;
  abcd[345] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_bc;
  abcd[346] = 4.0E0*I_ERI_F3y_Py_Dyz_S_bc;
  abcd[347] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_bc;
  abcd[348] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_bc;
  abcd[349] = 4.0E0*I_ERI_F3z_Py_Dyz_S_bc;
  abcd[350] = 4.0E0*I_ERI_F3x_Py_D2z_S_bc-2.0E0*1*I_ERI_F3x_Py_S_S_b;
  abcd[351] = 4.0E0*I_ERI_F2xy_Py_D2z_S_bc-2.0E0*1*I_ERI_F2xy_Py_S_S_b;
  abcd[352] = 4.0E0*I_ERI_F2xz_Py_D2z_S_bc-2.0E0*1*I_ERI_F2xz_Py_S_S_b;
  abcd[353] = 4.0E0*I_ERI_Fx2y_Py_D2z_S_bc-2.0E0*1*I_ERI_Fx2y_Py_S_S_b;
  abcd[354] = 4.0E0*I_ERI_Fxyz_Py_D2z_S_bc-2.0E0*1*I_ERI_Fxyz_Py_S_S_b;
  abcd[355] = 4.0E0*I_ERI_Fx2z_Py_D2z_S_bc-2.0E0*1*I_ERI_Fx2z_Py_S_S_b;
  abcd[356] = 4.0E0*I_ERI_F3y_Py_D2z_S_bc-2.0E0*1*I_ERI_F3y_Py_S_S_b;
  abcd[357] = 4.0E0*I_ERI_F2yz_Py_D2z_S_bc-2.0E0*1*I_ERI_F2yz_Py_S_S_b;
  abcd[358] = 4.0E0*I_ERI_Fy2z_Py_D2z_S_bc-2.0E0*1*I_ERI_Fy2z_Py_S_S_b;
  abcd[359] = 4.0E0*I_ERI_F3z_Py_D2z_S_bc-2.0E0*1*I_ERI_F3z_Py_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[360] = 4.0E0*I_ERI_F3x_Pz_D2x_S_bc-2.0E0*1*I_ERI_F3x_Pz_S_S_b;
  abcd[361] = 4.0E0*I_ERI_F2xy_Pz_D2x_S_bc-2.0E0*1*I_ERI_F2xy_Pz_S_S_b;
  abcd[362] = 4.0E0*I_ERI_F2xz_Pz_D2x_S_bc-2.0E0*1*I_ERI_F2xz_Pz_S_S_b;
  abcd[363] = 4.0E0*I_ERI_Fx2y_Pz_D2x_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_S_S_b;
  abcd[364] = 4.0E0*I_ERI_Fxyz_Pz_D2x_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_S_S_b;
  abcd[365] = 4.0E0*I_ERI_Fx2z_Pz_D2x_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_S_S_b;
  abcd[366] = 4.0E0*I_ERI_F3y_Pz_D2x_S_bc-2.0E0*1*I_ERI_F3y_Pz_S_S_b;
  abcd[367] = 4.0E0*I_ERI_F2yz_Pz_D2x_S_bc-2.0E0*1*I_ERI_F2yz_Pz_S_S_b;
  abcd[368] = 4.0E0*I_ERI_Fy2z_Pz_D2x_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_S_S_b;
  abcd[369] = 4.0E0*I_ERI_F3z_Pz_D2x_S_bc-2.0E0*1*I_ERI_F3z_Pz_S_S_b;
  abcd[370] = 4.0E0*I_ERI_F3x_Pz_Dxy_S_bc;
  abcd[371] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_bc;
  abcd[372] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_bc;
  abcd[373] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_bc;
  abcd[374] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_bc;
  abcd[375] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_bc;
  abcd[376] = 4.0E0*I_ERI_F3y_Pz_Dxy_S_bc;
  abcd[377] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_bc;
  abcd[378] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_bc;
  abcd[379] = 4.0E0*I_ERI_F3z_Pz_Dxy_S_bc;
  abcd[380] = 4.0E0*I_ERI_F3x_Pz_Dxz_S_bc;
  abcd[381] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_bc;
  abcd[382] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_bc;
  abcd[383] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_bc;
  abcd[384] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_bc;
  abcd[385] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_bc;
  abcd[386] = 4.0E0*I_ERI_F3y_Pz_Dxz_S_bc;
  abcd[387] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_bc;
  abcd[388] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_bc;
  abcd[389] = 4.0E0*I_ERI_F3z_Pz_Dxz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[390] = 4.0E0*I_ERI_F3x_Pz_Dxy_S_bc;
  abcd[391] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_bc;
  abcd[392] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_bc;
  abcd[393] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_bc;
  abcd[394] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_bc;
  abcd[395] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_bc;
  abcd[396] = 4.0E0*I_ERI_F3y_Pz_Dxy_S_bc;
  abcd[397] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_bc;
  abcd[398] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_bc;
  abcd[399] = 4.0E0*I_ERI_F3z_Pz_Dxy_S_bc;
  abcd[400] = 4.0E0*I_ERI_F3x_Pz_D2y_S_bc-2.0E0*1*I_ERI_F3x_Pz_S_S_b;
  abcd[401] = 4.0E0*I_ERI_F2xy_Pz_D2y_S_bc-2.0E0*1*I_ERI_F2xy_Pz_S_S_b;
  abcd[402] = 4.0E0*I_ERI_F2xz_Pz_D2y_S_bc-2.0E0*1*I_ERI_F2xz_Pz_S_S_b;
  abcd[403] = 4.0E0*I_ERI_Fx2y_Pz_D2y_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_S_S_b;
  abcd[404] = 4.0E0*I_ERI_Fxyz_Pz_D2y_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_S_S_b;
  abcd[405] = 4.0E0*I_ERI_Fx2z_Pz_D2y_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_S_S_b;
  abcd[406] = 4.0E0*I_ERI_F3y_Pz_D2y_S_bc-2.0E0*1*I_ERI_F3y_Pz_S_S_b;
  abcd[407] = 4.0E0*I_ERI_F2yz_Pz_D2y_S_bc-2.0E0*1*I_ERI_F2yz_Pz_S_S_b;
  abcd[408] = 4.0E0*I_ERI_Fy2z_Pz_D2y_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_S_S_b;
  abcd[409] = 4.0E0*I_ERI_F3z_Pz_D2y_S_bc-2.0E0*1*I_ERI_F3z_Pz_S_S_b;
  abcd[410] = 4.0E0*I_ERI_F3x_Pz_Dyz_S_bc;
  abcd[411] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_bc;
  abcd[412] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_bc;
  abcd[413] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_bc;
  abcd[414] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_bc;
  abcd[415] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_bc;
  abcd[416] = 4.0E0*I_ERI_F3y_Pz_Dyz_S_bc;
  abcd[417] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_bc;
  abcd[418] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_bc;
  abcd[419] = 4.0E0*I_ERI_F3z_Pz_Dyz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  abcd[420] = 4.0E0*I_ERI_F3x_Pz_Dxz_S_bc;
  abcd[421] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_bc;
  abcd[422] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_bc;
  abcd[423] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_bc;
  abcd[424] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_bc;
  abcd[425] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_bc;
  abcd[426] = 4.0E0*I_ERI_F3y_Pz_Dxz_S_bc;
  abcd[427] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_bc;
  abcd[428] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_bc;
  abcd[429] = 4.0E0*I_ERI_F3z_Pz_Dxz_S_bc;
  abcd[430] = 4.0E0*I_ERI_F3x_Pz_Dyz_S_bc;
  abcd[431] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_bc;
  abcd[432] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_bc;
  abcd[433] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_bc;
  abcd[434] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_bc;
  abcd[435] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_bc;
  abcd[436] = 4.0E0*I_ERI_F3y_Pz_Dyz_S_bc;
  abcd[437] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_bc;
  abcd[438] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_bc;
  abcd[439] = 4.0E0*I_ERI_F3z_Pz_Dyz_S_bc;
  abcd[440] = 4.0E0*I_ERI_F3x_Pz_D2z_S_bc-2.0E0*1*I_ERI_F3x_Pz_S_S_b;
  abcd[441] = 4.0E0*I_ERI_F2xy_Pz_D2z_S_bc-2.0E0*1*I_ERI_F2xy_Pz_S_S_b;
  abcd[442] = 4.0E0*I_ERI_F2xz_Pz_D2z_S_bc-2.0E0*1*I_ERI_F2xz_Pz_S_S_b;
  abcd[443] = 4.0E0*I_ERI_Fx2y_Pz_D2z_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_S_S_b;
  abcd[444] = 4.0E0*I_ERI_Fxyz_Pz_D2z_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_S_S_b;
  abcd[445] = 4.0E0*I_ERI_Fx2z_Pz_D2z_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_S_S_b;
  abcd[446] = 4.0E0*I_ERI_F3y_Pz_D2z_S_bc-2.0E0*1*I_ERI_F3y_Pz_S_S_b;
  abcd[447] = 4.0E0*I_ERI_F2yz_Pz_D2z_S_bc-2.0E0*1*I_ERI_F2yz_Pz_S_S_b;
  abcd[448] = 4.0E0*I_ERI_Fy2z_Pz_D2z_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_S_S_b;
  abcd[449] = 4.0E0*I_ERI_F3z_Pz_D2z_S_bc-2.0E0*1*I_ERI_F3z_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[450] = 4.0E0*I_ERI_F3x_Px_Px_Px_bd;
  abcd[451] = 4.0E0*I_ERI_F2xy_Px_Px_Px_bd;
  abcd[452] = 4.0E0*I_ERI_F2xz_Px_Px_Px_bd;
  abcd[453] = 4.0E0*I_ERI_Fx2y_Px_Px_Px_bd;
  abcd[454] = 4.0E0*I_ERI_Fxyz_Px_Px_Px_bd;
  abcd[455] = 4.0E0*I_ERI_Fx2z_Px_Px_Px_bd;
  abcd[456] = 4.0E0*I_ERI_F3y_Px_Px_Px_bd;
  abcd[457] = 4.0E0*I_ERI_F2yz_Px_Px_Px_bd;
  abcd[458] = 4.0E0*I_ERI_Fy2z_Px_Px_Px_bd;
  abcd[459] = 4.0E0*I_ERI_F3z_Px_Px_Px_bd;
  abcd[460] = 4.0E0*I_ERI_F3x_Px_Py_Px_bd;
  abcd[461] = 4.0E0*I_ERI_F2xy_Px_Py_Px_bd;
  abcd[462] = 4.0E0*I_ERI_F2xz_Px_Py_Px_bd;
  abcd[463] = 4.0E0*I_ERI_Fx2y_Px_Py_Px_bd;
  abcd[464] = 4.0E0*I_ERI_Fxyz_Px_Py_Px_bd;
  abcd[465] = 4.0E0*I_ERI_Fx2z_Px_Py_Px_bd;
  abcd[466] = 4.0E0*I_ERI_F3y_Px_Py_Px_bd;
  abcd[467] = 4.0E0*I_ERI_F2yz_Px_Py_Px_bd;
  abcd[468] = 4.0E0*I_ERI_Fy2z_Px_Py_Px_bd;
  abcd[469] = 4.0E0*I_ERI_F3z_Px_Py_Px_bd;
  abcd[470] = 4.0E0*I_ERI_F3x_Px_Pz_Px_bd;
  abcd[471] = 4.0E0*I_ERI_F2xy_Px_Pz_Px_bd;
  abcd[472] = 4.0E0*I_ERI_F2xz_Px_Pz_Px_bd;
  abcd[473] = 4.0E0*I_ERI_Fx2y_Px_Pz_Px_bd;
  abcd[474] = 4.0E0*I_ERI_Fxyz_Px_Pz_Px_bd;
  abcd[475] = 4.0E0*I_ERI_Fx2z_Px_Pz_Px_bd;
  abcd[476] = 4.0E0*I_ERI_F3y_Px_Pz_Px_bd;
  abcd[477] = 4.0E0*I_ERI_F2yz_Px_Pz_Px_bd;
  abcd[478] = 4.0E0*I_ERI_Fy2z_Px_Pz_Px_bd;
  abcd[479] = 4.0E0*I_ERI_F3z_Px_Pz_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[480] = 4.0E0*I_ERI_F3x_Px_Px_Py_bd;
  abcd[481] = 4.0E0*I_ERI_F2xy_Px_Px_Py_bd;
  abcd[482] = 4.0E0*I_ERI_F2xz_Px_Px_Py_bd;
  abcd[483] = 4.0E0*I_ERI_Fx2y_Px_Px_Py_bd;
  abcd[484] = 4.0E0*I_ERI_Fxyz_Px_Px_Py_bd;
  abcd[485] = 4.0E0*I_ERI_Fx2z_Px_Px_Py_bd;
  abcd[486] = 4.0E0*I_ERI_F3y_Px_Px_Py_bd;
  abcd[487] = 4.0E0*I_ERI_F2yz_Px_Px_Py_bd;
  abcd[488] = 4.0E0*I_ERI_Fy2z_Px_Px_Py_bd;
  abcd[489] = 4.0E0*I_ERI_F3z_Px_Px_Py_bd;
  abcd[490] = 4.0E0*I_ERI_F3x_Px_Py_Py_bd;
  abcd[491] = 4.0E0*I_ERI_F2xy_Px_Py_Py_bd;
  abcd[492] = 4.0E0*I_ERI_F2xz_Px_Py_Py_bd;
  abcd[493] = 4.0E0*I_ERI_Fx2y_Px_Py_Py_bd;
  abcd[494] = 4.0E0*I_ERI_Fxyz_Px_Py_Py_bd;
  abcd[495] = 4.0E0*I_ERI_Fx2z_Px_Py_Py_bd;
  abcd[496] = 4.0E0*I_ERI_F3y_Px_Py_Py_bd;
  abcd[497] = 4.0E0*I_ERI_F2yz_Px_Py_Py_bd;
  abcd[498] = 4.0E0*I_ERI_Fy2z_Px_Py_Py_bd;
  abcd[499] = 4.0E0*I_ERI_F3z_Px_Py_Py_bd;
  abcd[500] = 4.0E0*I_ERI_F3x_Px_Pz_Py_bd;
  abcd[501] = 4.0E0*I_ERI_F2xy_Px_Pz_Py_bd;
  abcd[502] = 4.0E0*I_ERI_F2xz_Px_Pz_Py_bd;
  abcd[503] = 4.0E0*I_ERI_Fx2y_Px_Pz_Py_bd;
  abcd[504] = 4.0E0*I_ERI_Fxyz_Px_Pz_Py_bd;
  abcd[505] = 4.0E0*I_ERI_Fx2z_Px_Pz_Py_bd;
  abcd[506] = 4.0E0*I_ERI_F3y_Px_Pz_Py_bd;
  abcd[507] = 4.0E0*I_ERI_F2yz_Px_Pz_Py_bd;
  abcd[508] = 4.0E0*I_ERI_Fy2z_Px_Pz_Py_bd;
  abcd[509] = 4.0E0*I_ERI_F3z_Px_Pz_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[510] = 4.0E0*I_ERI_F3x_Px_Px_Pz_bd;
  abcd[511] = 4.0E0*I_ERI_F2xy_Px_Px_Pz_bd;
  abcd[512] = 4.0E0*I_ERI_F2xz_Px_Px_Pz_bd;
  abcd[513] = 4.0E0*I_ERI_Fx2y_Px_Px_Pz_bd;
  abcd[514] = 4.0E0*I_ERI_Fxyz_Px_Px_Pz_bd;
  abcd[515] = 4.0E0*I_ERI_Fx2z_Px_Px_Pz_bd;
  abcd[516] = 4.0E0*I_ERI_F3y_Px_Px_Pz_bd;
  abcd[517] = 4.0E0*I_ERI_F2yz_Px_Px_Pz_bd;
  abcd[518] = 4.0E0*I_ERI_Fy2z_Px_Px_Pz_bd;
  abcd[519] = 4.0E0*I_ERI_F3z_Px_Px_Pz_bd;
  abcd[520] = 4.0E0*I_ERI_F3x_Px_Py_Pz_bd;
  abcd[521] = 4.0E0*I_ERI_F2xy_Px_Py_Pz_bd;
  abcd[522] = 4.0E0*I_ERI_F2xz_Px_Py_Pz_bd;
  abcd[523] = 4.0E0*I_ERI_Fx2y_Px_Py_Pz_bd;
  abcd[524] = 4.0E0*I_ERI_Fxyz_Px_Py_Pz_bd;
  abcd[525] = 4.0E0*I_ERI_Fx2z_Px_Py_Pz_bd;
  abcd[526] = 4.0E0*I_ERI_F3y_Px_Py_Pz_bd;
  abcd[527] = 4.0E0*I_ERI_F2yz_Px_Py_Pz_bd;
  abcd[528] = 4.0E0*I_ERI_Fy2z_Px_Py_Pz_bd;
  abcd[529] = 4.0E0*I_ERI_F3z_Px_Py_Pz_bd;
  abcd[530] = 4.0E0*I_ERI_F3x_Px_Pz_Pz_bd;
  abcd[531] = 4.0E0*I_ERI_F2xy_Px_Pz_Pz_bd;
  abcd[532] = 4.0E0*I_ERI_F2xz_Px_Pz_Pz_bd;
  abcd[533] = 4.0E0*I_ERI_Fx2y_Px_Pz_Pz_bd;
  abcd[534] = 4.0E0*I_ERI_Fxyz_Px_Pz_Pz_bd;
  abcd[535] = 4.0E0*I_ERI_Fx2z_Px_Pz_Pz_bd;
  abcd[536] = 4.0E0*I_ERI_F3y_Px_Pz_Pz_bd;
  abcd[537] = 4.0E0*I_ERI_F2yz_Px_Pz_Pz_bd;
  abcd[538] = 4.0E0*I_ERI_Fy2z_Px_Pz_Pz_bd;
  abcd[539] = 4.0E0*I_ERI_F3z_Px_Pz_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[540] = 4.0E0*I_ERI_F3x_Py_Px_Px_bd;
  abcd[541] = 4.0E0*I_ERI_F2xy_Py_Px_Px_bd;
  abcd[542] = 4.0E0*I_ERI_F2xz_Py_Px_Px_bd;
  abcd[543] = 4.0E0*I_ERI_Fx2y_Py_Px_Px_bd;
  abcd[544] = 4.0E0*I_ERI_Fxyz_Py_Px_Px_bd;
  abcd[545] = 4.0E0*I_ERI_Fx2z_Py_Px_Px_bd;
  abcd[546] = 4.0E0*I_ERI_F3y_Py_Px_Px_bd;
  abcd[547] = 4.0E0*I_ERI_F2yz_Py_Px_Px_bd;
  abcd[548] = 4.0E0*I_ERI_Fy2z_Py_Px_Px_bd;
  abcd[549] = 4.0E0*I_ERI_F3z_Py_Px_Px_bd;
  abcd[550] = 4.0E0*I_ERI_F3x_Py_Py_Px_bd;
  abcd[551] = 4.0E0*I_ERI_F2xy_Py_Py_Px_bd;
  abcd[552] = 4.0E0*I_ERI_F2xz_Py_Py_Px_bd;
  abcd[553] = 4.0E0*I_ERI_Fx2y_Py_Py_Px_bd;
  abcd[554] = 4.0E0*I_ERI_Fxyz_Py_Py_Px_bd;
  abcd[555] = 4.0E0*I_ERI_Fx2z_Py_Py_Px_bd;
  abcd[556] = 4.0E0*I_ERI_F3y_Py_Py_Px_bd;
  abcd[557] = 4.0E0*I_ERI_F2yz_Py_Py_Px_bd;
  abcd[558] = 4.0E0*I_ERI_Fy2z_Py_Py_Px_bd;
  abcd[559] = 4.0E0*I_ERI_F3z_Py_Py_Px_bd;
  abcd[560] = 4.0E0*I_ERI_F3x_Py_Pz_Px_bd;
  abcd[561] = 4.0E0*I_ERI_F2xy_Py_Pz_Px_bd;
  abcd[562] = 4.0E0*I_ERI_F2xz_Py_Pz_Px_bd;
  abcd[563] = 4.0E0*I_ERI_Fx2y_Py_Pz_Px_bd;
  abcd[564] = 4.0E0*I_ERI_Fxyz_Py_Pz_Px_bd;
  abcd[565] = 4.0E0*I_ERI_Fx2z_Py_Pz_Px_bd;
  abcd[566] = 4.0E0*I_ERI_F3y_Py_Pz_Px_bd;
  abcd[567] = 4.0E0*I_ERI_F2yz_Py_Pz_Px_bd;
  abcd[568] = 4.0E0*I_ERI_Fy2z_Py_Pz_Px_bd;
  abcd[569] = 4.0E0*I_ERI_F3z_Py_Pz_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[570] = 4.0E0*I_ERI_F3x_Py_Px_Py_bd;
  abcd[571] = 4.0E0*I_ERI_F2xy_Py_Px_Py_bd;
  abcd[572] = 4.0E0*I_ERI_F2xz_Py_Px_Py_bd;
  abcd[573] = 4.0E0*I_ERI_Fx2y_Py_Px_Py_bd;
  abcd[574] = 4.0E0*I_ERI_Fxyz_Py_Px_Py_bd;
  abcd[575] = 4.0E0*I_ERI_Fx2z_Py_Px_Py_bd;
  abcd[576] = 4.0E0*I_ERI_F3y_Py_Px_Py_bd;
  abcd[577] = 4.0E0*I_ERI_F2yz_Py_Px_Py_bd;
  abcd[578] = 4.0E0*I_ERI_Fy2z_Py_Px_Py_bd;
  abcd[579] = 4.0E0*I_ERI_F3z_Py_Px_Py_bd;
  abcd[580] = 4.0E0*I_ERI_F3x_Py_Py_Py_bd;
  abcd[581] = 4.0E0*I_ERI_F2xy_Py_Py_Py_bd;
  abcd[582] = 4.0E0*I_ERI_F2xz_Py_Py_Py_bd;
  abcd[583] = 4.0E0*I_ERI_Fx2y_Py_Py_Py_bd;
  abcd[584] = 4.0E0*I_ERI_Fxyz_Py_Py_Py_bd;
  abcd[585] = 4.0E0*I_ERI_Fx2z_Py_Py_Py_bd;
  abcd[586] = 4.0E0*I_ERI_F3y_Py_Py_Py_bd;
  abcd[587] = 4.0E0*I_ERI_F2yz_Py_Py_Py_bd;
  abcd[588] = 4.0E0*I_ERI_Fy2z_Py_Py_Py_bd;
  abcd[589] = 4.0E0*I_ERI_F3z_Py_Py_Py_bd;
  abcd[590] = 4.0E0*I_ERI_F3x_Py_Pz_Py_bd;
  abcd[591] = 4.0E0*I_ERI_F2xy_Py_Pz_Py_bd;
  abcd[592] = 4.0E0*I_ERI_F2xz_Py_Pz_Py_bd;
  abcd[593] = 4.0E0*I_ERI_Fx2y_Py_Pz_Py_bd;
  abcd[594] = 4.0E0*I_ERI_Fxyz_Py_Pz_Py_bd;
  abcd[595] = 4.0E0*I_ERI_Fx2z_Py_Pz_Py_bd;
  abcd[596] = 4.0E0*I_ERI_F3y_Py_Pz_Py_bd;
  abcd[597] = 4.0E0*I_ERI_F2yz_Py_Pz_Py_bd;
  abcd[598] = 4.0E0*I_ERI_Fy2z_Py_Pz_Py_bd;
  abcd[599] = 4.0E0*I_ERI_F3z_Py_Pz_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dby_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[600] = 4.0E0*I_ERI_F3x_Py_Px_Pz_bd;
  abcd[601] = 4.0E0*I_ERI_F2xy_Py_Px_Pz_bd;
  abcd[602] = 4.0E0*I_ERI_F2xz_Py_Px_Pz_bd;
  abcd[603] = 4.0E0*I_ERI_Fx2y_Py_Px_Pz_bd;
  abcd[604] = 4.0E0*I_ERI_Fxyz_Py_Px_Pz_bd;
  abcd[605] = 4.0E0*I_ERI_Fx2z_Py_Px_Pz_bd;
  abcd[606] = 4.0E0*I_ERI_F3y_Py_Px_Pz_bd;
  abcd[607] = 4.0E0*I_ERI_F2yz_Py_Px_Pz_bd;
  abcd[608] = 4.0E0*I_ERI_Fy2z_Py_Px_Pz_bd;
  abcd[609] = 4.0E0*I_ERI_F3z_Py_Px_Pz_bd;
  abcd[610] = 4.0E0*I_ERI_F3x_Py_Py_Pz_bd;
  abcd[611] = 4.0E0*I_ERI_F2xy_Py_Py_Pz_bd;
  abcd[612] = 4.0E0*I_ERI_F2xz_Py_Py_Pz_bd;
  abcd[613] = 4.0E0*I_ERI_Fx2y_Py_Py_Pz_bd;
  abcd[614] = 4.0E0*I_ERI_Fxyz_Py_Py_Pz_bd;
  abcd[615] = 4.0E0*I_ERI_Fx2z_Py_Py_Pz_bd;
  abcd[616] = 4.0E0*I_ERI_F3y_Py_Py_Pz_bd;
  abcd[617] = 4.0E0*I_ERI_F2yz_Py_Py_Pz_bd;
  abcd[618] = 4.0E0*I_ERI_Fy2z_Py_Py_Pz_bd;
  abcd[619] = 4.0E0*I_ERI_F3z_Py_Py_Pz_bd;
  abcd[620] = 4.0E0*I_ERI_F3x_Py_Pz_Pz_bd;
  abcd[621] = 4.0E0*I_ERI_F2xy_Py_Pz_Pz_bd;
  abcd[622] = 4.0E0*I_ERI_F2xz_Py_Pz_Pz_bd;
  abcd[623] = 4.0E0*I_ERI_Fx2y_Py_Pz_Pz_bd;
  abcd[624] = 4.0E0*I_ERI_Fxyz_Py_Pz_Pz_bd;
  abcd[625] = 4.0E0*I_ERI_Fx2z_Py_Pz_Pz_bd;
  abcd[626] = 4.0E0*I_ERI_F3y_Py_Pz_Pz_bd;
  abcd[627] = 4.0E0*I_ERI_F2yz_Py_Pz_Pz_bd;
  abcd[628] = 4.0E0*I_ERI_Fy2z_Py_Pz_Pz_bd;
  abcd[629] = 4.0E0*I_ERI_F3z_Py_Pz_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[630] = 4.0E0*I_ERI_F3x_Pz_Px_Px_bd;
  abcd[631] = 4.0E0*I_ERI_F2xy_Pz_Px_Px_bd;
  abcd[632] = 4.0E0*I_ERI_F2xz_Pz_Px_Px_bd;
  abcd[633] = 4.0E0*I_ERI_Fx2y_Pz_Px_Px_bd;
  abcd[634] = 4.0E0*I_ERI_Fxyz_Pz_Px_Px_bd;
  abcd[635] = 4.0E0*I_ERI_Fx2z_Pz_Px_Px_bd;
  abcd[636] = 4.0E0*I_ERI_F3y_Pz_Px_Px_bd;
  abcd[637] = 4.0E0*I_ERI_F2yz_Pz_Px_Px_bd;
  abcd[638] = 4.0E0*I_ERI_Fy2z_Pz_Px_Px_bd;
  abcd[639] = 4.0E0*I_ERI_F3z_Pz_Px_Px_bd;
  abcd[640] = 4.0E0*I_ERI_F3x_Pz_Py_Px_bd;
  abcd[641] = 4.0E0*I_ERI_F2xy_Pz_Py_Px_bd;
  abcd[642] = 4.0E0*I_ERI_F2xz_Pz_Py_Px_bd;
  abcd[643] = 4.0E0*I_ERI_Fx2y_Pz_Py_Px_bd;
  abcd[644] = 4.0E0*I_ERI_Fxyz_Pz_Py_Px_bd;
  abcd[645] = 4.0E0*I_ERI_Fx2z_Pz_Py_Px_bd;
  abcd[646] = 4.0E0*I_ERI_F3y_Pz_Py_Px_bd;
  abcd[647] = 4.0E0*I_ERI_F2yz_Pz_Py_Px_bd;
  abcd[648] = 4.0E0*I_ERI_Fy2z_Pz_Py_Px_bd;
  abcd[649] = 4.0E0*I_ERI_F3z_Pz_Py_Px_bd;
  abcd[650] = 4.0E0*I_ERI_F3x_Pz_Pz_Px_bd;
  abcd[651] = 4.0E0*I_ERI_F2xy_Pz_Pz_Px_bd;
  abcd[652] = 4.0E0*I_ERI_F2xz_Pz_Pz_Px_bd;
  abcd[653] = 4.0E0*I_ERI_Fx2y_Pz_Pz_Px_bd;
  abcd[654] = 4.0E0*I_ERI_Fxyz_Pz_Pz_Px_bd;
  abcd[655] = 4.0E0*I_ERI_Fx2z_Pz_Pz_Px_bd;
  abcd[656] = 4.0E0*I_ERI_F3y_Pz_Pz_Px_bd;
  abcd[657] = 4.0E0*I_ERI_F2yz_Pz_Pz_Px_bd;
  abcd[658] = 4.0E0*I_ERI_Fy2z_Pz_Pz_Px_bd;
  abcd[659] = 4.0E0*I_ERI_F3z_Pz_Pz_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[660] = 4.0E0*I_ERI_F3x_Pz_Px_Py_bd;
  abcd[661] = 4.0E0*I_ERI_F2xy_Pz_Px_Py_bd;
  abcd[662] = 4.0E0*I_ERI_F2xz_Pz_Px_Py_bd;
  abcd[663] = 4.0E0*I_ERI_Fx2y_Pz_Px_Py_bd;
  abcd[664] = 4.0E0*I_ERI_Fxyz_Pz_Px_Py_bd;
  abcd[665] = 4.0E0*I_ERI_Fx2z_Pz_Px_Py_bd;
  abcd[666] = 4.0E0*I_ERI_F3y_Pz_Px_Py_bd;
  abcd[667] = 4.0E0*I_ERI_F2yz_Pz_Px_Py_bd;
  abcd[668] = 4.0E0*I_ERI_Fy2z_Pz_Px_Py_bd;
  abcd[669] = 4.0E0*I_ERI_F3z_Pz_Px_Py_bd;
  abcd[670] = 4.0E0*I_ERI_F3x_Pz_Py_Py_bd;
  abcd[671] = 4.0E0*I_ERI_F2xy_Pz_Py_Py_bd;
  abcd[672] = 4.0E0*I_ERI_F2xz_Pz_Py_Py_bd;
  abcd[673] = 4.0E0*I_ERI_Fx2y_Pz_Py_Py_bd;
  abcd[674] = 4.0E0*I_ERI_Fxyz_Pz_Py_Py_bd;
  abcd[675] = 4.0E0*I_ERI_Fx2z_Pz_Py_Py_bd;
  abcd[676] = 4.0E0*I_ERI_F3y_Pz_Py_Py_bd;
  abcd[677] = 4.0E0*I_ERI_F2yz_Pz_Py_Py_bd;
  abcd[678] = 4.0E0*I_ERI_Fy2z_Pz_Py_Py_bd;
  abcd[679] = 4.0E0*I_ERI_F3z_Pz_Py_Py_bd;
  abcd[680] = 4.0E0*I_ERI_F3x_Pz_Pz_Py_bd;
  abcd[681] = 4.0E0*I_ERI_F2xy_Pz_Pz_Py_bd;
  abcd[682] = 4.0E0*I_ERI_F2xz_Pz_Pz_Py_bd;
  abcd[683] = 4.0E0*I_ERI_Fx2y_Pz_Pz_Py_bd;
  abcd[684] = 4.0E0*I_ERI_Fxyz_Pz_Pz_Py_bd;
  abcd[685] = 4.0E0*I_ERI_Fx2z_Pz_Pz_Py_bd;
  abcd[686] = 4.0E0*I_ERI_F3y_Pz_Pz_Py_bd;
  abcd[687] = 4.0E0*I_ERI_F2yz_Pz_Pz_Py_bd;
  abcd[688] = 4.0E0*I_ERI_Fy2z_Pz_Pz_Py_bd;
  abcd[689] = 4.0E0*I_ERI_F3z_Pz_Pz_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dbz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_bd
   ************************************************************/
  abcd[690] = 4.0E0*I_ERI_F3x_Pz_Px_Pz_bd;
  abcd[691] = 4.0E0*I_ERI_F2xy_Pz_Px_Pz_bd;
  abcd[692] = 4.0E0*I_ERI_F2xz_Pz_Px_Pz_bd;
  abcd[693] = 4.0E0*I_ERI_Fx2y_Pz_Px_Pz_bd;
  abcd[694] = 4.0E0*I_ERI_Fxyz_Pz_Px_Pz_bd;
  abcd[695] = 4.0E0*I_ERI_Fx2z_Pz_Px_Pz_bd;
  abcd[696] = 4.0E0*I_ERI_F3y_Pz_Px_Pz_bd;
  abcd[697] = 4.0E0*I_ERI_F2yz_Pz_Px_Pz_bd;
  abcd[698] = 4.0E0*I_ERI_Fy2z_Pz_Px_Pz_bd;
  abcd[699] = 4.0E0*I_ERI_F3z_Pz_Px_Pz_bd;
  abcd[700] = 4.0E0*I_ERI_F3x_Pz_Py_Pz_bd;
  abcd[701] = 4.0E0*I_ERI_F2xy_Pz_Py_Pz_bd;
  abcd[702] = 4.0E0*I_ERI_F2xz_Pz_Py_Pz_bd;
  abcd[703] = 4.0E0*I_ERI_Fx2y_Pz_Py_Pz_bd;
  abcd[704] = 4.0E0*I_ERI_Fxyz_Pz_Py_Pz_bd;
  abcd[705] = 4.0E0*I_ERI_Fx2z_Pz_Py_Pz_bd;
  abcd[706] = 4.0E0*I_ERI_F3y_Pz_Py_Pz_bd;
  abcd[707] = 4.0E0*I_ERI_F2yz_Pz_Py_Pz_bd;
  abcd[708] = 4.0E0*I_ERI_Fy2z_Pz_Py_Pz_bd;
  abcd[709] = 4.0E0*I_ERI_F3z_Pz_Py_Pz_bd;
  abcd[710] = 4.0E0*I_ERI_F3x_Pz_Pz_Pz_bd;
  abcd[711] = 4.0E0*I_ERI_F2xy_Pz_Pz_Pz_bd;
  abcd[712] = 4.0E0*I_ERI_F2xz_Pz_Pz_Pz_bd;
  abcd[713] = 4.0E0*I_ERI_Fx2y_Pz_Pz_Pz_bd;
  abcd[714] = 4.0E0*I_ERI_Fxyz_Pz_Pz_Pz_bd;
  abcd[715] = 4.0E0*I_ERI_Fx2z_Pz_Pz_Pz_bd;
  abcd[716] = 4.0E0*I_ERI_F3y_Pz_Pz_Pz_bd;
  abcd[717] = 4.0E0*I_ERI_F2yz_Pz_Pz_Pz_bd;
  abcd[718] = 4.0E0*I_ERI_Fy2z_Pz_Pz_Pz_bd;
  abcd[719] = 4.0E0*I_ERI_F3z_Pz_Pz_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[720] = 4.0E0*I_ERI_F3x_S_F3x_S_cc-2.0E0*1*I_ERI_F3x_S_Px_S_c-2.0E0*2*I_ERI_F3x_S_Px_S_c;
  abcd[721] = 4.0E0*I_ERI_F2xy_S_F3x_S_cc-2.0E0*1*I_ERI_F2xy_S_Px_S_c-2.0E0*2*I_ERI_F2xy_S_Px_S_c;
  abcd[722] = 4.0E0*I_ERI_F2xz_S_F3x_S_cc-2.0E0*1*I_ERI_F2xz_S_Px_S_c-2.0E0*2*I_ERI_F2xz_S_Px_S_c;
  abcd[723] = 4.0E0*I_ERI_Fx2y_S_F3x_S_cc-2.0E0*1*I_ERI_Fx2y_S_Px_S_c-2.0E0*2*I_ERI_Fx2y_S_Px_S_c;
  abcd[724] = 4.0E0*I_ERI_Fxyz_S_F3x_S_cc-2.0E0*1*I_ERI_Fxyz_S_Px_S_c-2.0E0*2*I_ERI_Fxyz_S_Px_S_c;
  abcd[725] = 4.0E0*I_ERI_Fx2z_S_F3x_S_cc-2.0E0*1*I_ERI_Fx2z_S_Px_S_c-2.0E0*2*I_ERI_Fx2z_S_Px_S_c;
  abcd[726] = 4.0E0*I_ERI_F3y_S_F3x_S_cc-2.0E0*1*I_ERI_F3y_S_Px_S_c-2.0E0*2*I_ERI_F3y_S_Px_S_c;
  abcd[727] = 4.0E0*I_ERI_F2yz_S_F3x_S_cc-2.0E0*1*I_ERI_F2yz_S_Px_S_c-2.0E0*2*I_ERI_F2yz_S_Px_S_c;
  abcd[728] = 4.0E0*I_ERI_Fy2z_S_F3x_S_cc-2.0E0*1*I_ERI_Fy2z_S_Px_S_c-2.0E0*2*I_ERI_Fy2z_S_Px_S_c;
  abcd[729] = 4.0E0*I_ERI_F3z_S_F3x_S_cc-2.0E0*1*I_ERI_F3z_S_Px_S_c-2.0E0*2*I_ERI_F3z_S_Px_S_c;
  abcd[730] = 4.0E0*I_ERI_F3x_S_F2xy_S_cc-2.0E0*1*I_ERI_F3x_S_Py_S_c;
  abcd[731] = 4.0E0*I_ERI_F2xy_S_F2xy_S_cc-2.0E0*1*I_ERI_F2xy_S_Py_S_c;
  abcd[732] = 4.0E0*I_ERI_F2xz_S_F2xy_S_cc-2.0E0*1*I_ERI_F2xz_S_Py_S_c;
  abcd[733] = 4.0E0*I_ERI_Fx2y_S_F2xy_S_cc-2.0E0*1*I_ERI_Fx2y_S_Py_S_c;
  abcd[734] = 4.0E0*I_ERI_Fxyz_S_F2xy_S_cc-2.0E0*1*I_ERI_Fxyz_S_Py_S_c;
  abcd[735] = 4.0E0*I_ERI_Fx2z_S_F2xy_S_cc-2.0E0*1*I_ERI_Fx2z_S_Py_S_c;
  abcd[736] = 4.0E0*I_ERI_F3y_S_F2xy_S_cc-2.0E0*1*I_ERI_F3y_S_Py_S_c;
  abcd[737] = 4.0E0*I_ERI_F2yz_S_F2xy_S_cc-2.0E0*1*I_ERI_F2yz_S_Py_S_c;
  abcd[738] = 4.0E0*I_ERI_Fy2z_S_F2xy_S_cc-2.0E0*1*I_ERI_Fy2z_S_Py_S_c;
  abcd[739] = 4.0E0*I_ERI_F3z_S_F2xy_S_cc-2.0E0*1*I_ERI_F3z_S_Py_S_c;
  abcd[740] = 4.0E0*I_ERI_F3x_S_F2xz_S_cc-2.0E0*1*I_ERI_F3x_S_Pz_S_c;
  abcd[741] = 4.0E0*I_ERI_F2xy_S_F2xz_S_cc-2.0E0*1*I_ERI_F2xy_S_Pz_S_c;
  abcd[742] = 4.0E0*I_ERI_F2xz_S_F2xz_S_cc-2.0E0*1*I_ERI_F2xz_S_Pz_S_c;
  abcd[743] = 4.0E0*I_ERI_Fx2y_S_F2xz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Pz_S_c;
  abcd[744] = 4.0E0*I_ERI_Fxyz_S_F2xz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Pz_S_c;
  abcd[745] = 4.0E0*I_ERI_Fx2z_S_F2xz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Pz_S_c;
  abcd[746] = 4.0E0*I_ERI_F3y_S_F2xz_S_cc-2.0E0*1*I_ERI_F3y_S_Pz_S_c;
  abcd[747] = 4.0E0*I_ERI_F2yz_S_F2xz_S_cc-2.0E0*1*I_ERI_F2yz_S_Pz_S_c;
  abcd[748] = 4.0E0*I_ERI_Fy2z_S_F2xz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Pz_S_c;
  abcd[749] = 4.0E0*I_ERI_F3z_S_F2xz_S_cc-2.0E0*1*I_ERI_F3z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[750] = 4.0E0*I_ERI_F3x_S_F2xy_S_cc-2.0E0*1*I_ERI_F3x_S_Py_S_c;
  abcd[751] = 4.0E0*I_ERI_F2xy_S_F2xy_S_cc-2.0E0*1*I_ERI_F2xy_S_Py_S_c;
  abcd[752] = 4.0E0*I_ERI_F2xz_S_F2xy_S_cc-2.0E0*1*I_ERI_F2xz_S_Py_S_c;
  abcd[753] = 4.0E0*I_ERI_Fx2y_S_F2xy_S_cc-2.0E0*1*I_ERI_Fx2y_S_Py_S_c;
  abcd[754] = 4.0E0*I_ERI_Fxyz_S_F2xy_S_cc-2.0E0*1*I_ERI_Fxyz_S_Py_S_c;
  abcd[755] = 4.0E0*I_ERI_Fx2z_S_F2xy_S_cc-2.0E0*1*I_ERI_Fx2z_S_Py_S_c;
  abcd[756] = 4.0E0*I_ERI_F3y_S_F2xy_S_cc-2.0E0*1*I_ERI_F3y_S_Py_S_c;
  abcd[757] = 4.0E0*I_ERI_F2yz_S_F2xy_S_cc-2.0E0*1*I_ERI_F2yz_S_Py_S_c;
  abcd[758] = 4.0E0*I_ERI_Fy2z_S_F2xy_S_cc-2.0E0*1*I_ERI_Fy2z_S_Py_S_c;
  abcd[759] = 4.0E0*I_ERI_F3z_S_F2xy_S_cc-2.0E0*1*I_ERI_F3z_S_Py_S_c;
  abcd[760] = 4.0E0*I_ERI_F3x_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3x_S_Px_S_c;
  abcd[761] = 4.0E0*I_ERI_F2xy_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2xy_S_Px_S_c;
  abcd[762] = 4.0E0*I_ERI_F2xz_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2xz_S_Px_S_c;
  abcd[763] = 4.0E0*I_ERI_Fx2y_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fx2y_S_Px_S_c;
  abcd[764] = 4.0E0*I_ERI_Fxyz_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fxyz_S_Px_S_c;
  abcd[765] = 4.0E0*I_ERI_Fx2z_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fx2z_S_Px_S_c;
  abcd[766] = 4.0E0*I_ERI_F3y_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3y_S_Px_S_c;
  abcd[767] = 4.0E0*I_ERI_F2yz_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2yz_S_Px_S_c;
  abcd[768] = 4.0E0*I_ERI_Fy2z_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fy2z_S_Px_S_c;
  abcd[769] = 4.0E0*I_ERI_F3z_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3z_S_Px_S_c;
  abcd[770] = 4.0E0*I_ERI_F3x_S_Fxyz_S_cc;
  abcd[771] = 4.0E0*I_ERI_F2xy_S_Fxyz_S_cc;
  abcd[772] = 4.0E0*I_ERI_F2xz_S_Fxyz_S_cc;
  abcd[773] = 4.0E0*I_ERI_Fx2y_S_Fxyz_S_cc;
  abcd[774] = 4.0E0*I_ERI_Fxyz_S_Fxyz_S_cc;
  abcd[775] = 4.0E0*I_ERI_Fx2z_S_Fxyz_S_cc;
  abcd[776] = 4.0E0*I_ERI_F3y_S_Fxyz_S_cc;
  abcd[777] = 4.0E0*I_ERI_F2yz_S_Fxyz_S_cc;
  abcd[778] = 4.0E0*I_ERI_Fy2z_S_Fxyz_S_cc;
  abcd[779] = 4.0E0*I_ERI_F3z_S_Fxyz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[780] = 4.0E0*I_ERI_F3x_S_F2xz_S_cc-2.0E0*1*I_ERI_F3x_S_Pz_S_c;
  abcd[781] = 4.0E0*I_ERI_F2xy_S_F2xz_S_cc-2.0E0*1*I_ERI_F2xy_S_Pz_S_c;
  abcd[782] = 4.0E0*I_ERI_F2xz_S_F2xz_S_cc-2.0E0*1*I_ERI_F2xz_S_Pz_S_c;
  abcd[783] = 4.0E0*I_ERI_Fx2y_S_F2xz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Pz_S_c;
  abcd[784] = 4.0E0*I_ERI_Fxyz_S_F2xz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Pz_S_c;
  abcd[785] = 4.0E0*I_ERI_Fx2z_S_F2xz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Pz_S_c;
  abcd[786] = 4.0E0*I_ERI_F3y_S_F2xz_S_cc-2.0E0*1*I_ERI_F3y_S_Pz_S_c;
  abcd[787] = 4.0E0*I_ERI_F2yz_S_F2xz_S_cc-2.0E0*1*I_ERI_F2yz_S_Pz_S_c;
  abcd[788] = 4.0E0*I_ERI_Fy2z_S_F2xz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Pz_S_c;
  abcd[789] = 4.0E0*I_ERI_F3z_S_F2xz_S_cc-2.0E0*1*I_ERI_F3z_S_Pz_S_c;
  abcd[790] = 4.0E0*I_ERI_F3x_S_Fxyz_S_cc;
  abcd[791] = 4.0E0*I_ERI_F2xy_S_Fxyz_S_cc;
  abcd[792] = 4.0E0*I_ERI_F2xz_S_Fxyz_S_cc;
  abcd[793] = 4.0E0*I_ERI_Fx2y_S_Fxyz_S_cc;
  abcd[794] = 4.0E0*I_ERI_Fxyz_S_Fxyz_S_cc;
  abcd[795] = 4.0E0*I_ERI_Fx2z_S_Fxyz_S_cc;
  abcd[796] = 4.0E0*I_ERI_F3y_S_Fxyz_S_cc;
  abcd[797] = 4.0E0*I_ERI_F2yz_S_Fxyz_S_cc;
  abcd[798] = 4.0E0*I_ERI_Fy2z_S_Fxyz_S_cc;
  abcd[799] = 4.0E0*I_ERI_F3z_S_Fxyz_S_cc;
  abcd[800] = 4.0E0*I_ERI_F3x_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3x_S_Px_S_c;
  abcd[801] = 4.0E0*I_ERI_F2xy_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Px_S_c;
  abcd[802] = 4.0E0*I_ERI_F2xz_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Px_S_c;
  abcd[803] = 4.0E0*I_ERI_Fx2y_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Px_S_c;
  abcd[804] = 4.0E0*I_ERI_Fxyz_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Px_S_c;
  abcd[805] = 4.0E0*I_ERI_Fx2z_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Px_S_c;
  abcd[806] = 4.0E0*I_ERI_F3y_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3y_S_Px_S_c;
  abcd[807] = 4.0E0*I_ERI_F2yz_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Px_S_c;
  abcd[808] = 4.0E0*I_ERI_Fy2z_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Px_S_c;
  abcd[809] = 4.0E0*I_ERI_F3z_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3z_S_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[810] = 4.0E0*I_ERI_F3x_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3x_S_Px_S_c;
  abcd[811] = 4.0E0*I_ERI_F2xy_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2xy_S_Px_S_c;
  abcd[812] = 4.0E0*I_ERI_F2xz_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2xz_S_Px_S_c;
  abcd[813] = 4.0E0*I_ERI_Fx2y_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fx2y_S_Px_S_c;
  abcd[814] = 4.0E0*I_ERI_Fxyz_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fxyz_S_Px_S_c;
  abcd[815] = 4.0E0*I_ERI_Fx2z_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fx2z_S_Px_S_c;
  abcd[816] = 4.0E0*I_ERI_F3y_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3y_S_Px_S_c;
  abcd[817] = 4.0E0*I_ERI_F2yz_S_Fx2y_S_cc-2.0E0*1*I_ERI_F2yz_S_Px_S_c;
  abcd[818] = 4.0E0*I_ERI_Fy2z_S_Fx2y_S_cc-2.0E0*1*I_ERI_Fy2z_S_Px_S_c;
  abcd[819] = 4.0E0*I_ERI_F3z_S_Fx2y_S_cc-2.0E0*1*I_ERI_F3z_S_Px_S_c;
  abcd[820] = 4.0E0*I_ERI_F3x_S_F3y_S_cc-2.0E0*1*I_ERI_F3x_S_Py_S_c-2.0E0*2*I_ERI_F3x_S_Py_S_c;
  abcd[821] = 4.0E0*I_ERI_F2xy_S_F3y_S_cc-2.0E0*1*I_ERI_F2xy_S_Py_S_c-2.0E0*2*I_ERI_F2xy_S_Py_S_c;
  abcd[822] = 4.0E0*I_ERI_F2xz_S_F3y_S_cc-2.0E0*1*I_ERI_F2xz_S_Py_S_c-2.0E0*2*I_ERI_F2xz_S_Py_S_c;
  abcd[823] = 4.0E0*I_ERI_Fx2y_S_F3y_S_cc-2.0E0*1*I_ERI_Fx2y_S_Py_S_c-2.0E0*2*I_ERI_Fx2y_S_Py_S_c;
  abcd[824] = 4.0E0*I_ERI_Fxyz_S_F3y_S_cc-2.0E0*1*I_ERI_Fxyz_S_Py_S_c-2.0E0*2*I_ERI_Fxyz_S_Py_S_c;
  abcd[825] = 4.0E0*I_ERI_Fx2z_S_F3y_S_cc-2.0E0*1*I_ERI_Fx2z_S_Py_S_c-2.0E0*2*I_ERI_Fx2z_S_Py_S_c;
  abcd[826] = 4.0E0*I_ERI_F3y_S_F3y_S_cc-2.0E0*1*I_ERI_F3y_S_Py_S_c-2.0E0*2*I_ERI_F3y_S_Py_S_c;
  abcd[827] = 4.0E0*I_ERI_F2yz_S_F3y_S_cc-2.0E0*1*I_ERI_F2yz_S_Py_S_c-2.0E0*2*I_ERI_F2yz_S_Py_S_c;
  abcd[828] = 4.0E0*I_ERI_Fy2z_S_F3y_S_cc-2.0E0*1*I_ERI_Fy2z_S_Py_S_c-2.0E0*2*I_ERI_Fy2z_S_Py_S_c;
  abcd[829] = 4.0E0*I_ERI_F3z_S_F3y_S_cc-2.0E0*1*I_ERI_F3z_S_Py_S_c-2.0E0*2*I_ERI_F3z_S_Py_S_c;
  abcd[830] = 4.0E0*I_ERI_F3x_S_F2yz_S_cc-2.0E0*1*I_ERI_F3x_S_Pz_S_c;
  abcd[831] = 4.0E0*I_ERI_F2xy_S_F2yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Pz_S_c;
  abcd[832] = 4.0E0*I_ERI_F2xz_S_F2yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Pz_S_c;
  abcd[833] = 4.0E0*I_ERI_Fx2y_S_F2yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Pz_S_c;
  abcd[834] = 4.0E0*I_ERI_Fxyz_S_F2yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Pz_S_c;
  abcd[835] = 4.0E0*I_ERI_Fx2z_S_F2yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Pz_S_c;
  abcd[836] = 4.0E0*I_ERI_F3y_S_F2yz_S_cc-2.0E0*1*I_ERI_F3y_S_Pz_S_c;
  abcd[837] = 4.0E0*I_ERI_F2yz_S_F2yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Pz_S_c;
  abcd[838] = 4.0E0*I_ERI_Fy2z_S_F2yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Pz_S_c;
  abcd[839] = 4.0E0*I_ERI_F3z_S_F2yz_S_cc-2.0E0*1*I_ERI_F3z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[840] = 4.0E0*I_ERI_F3x_S_Fxyz_S_cc;
  abcd[841] = 4.0E0*I_ERI_F2xy_S_Fxyz_S_cc;
  abcd[842] = 4.0E0*I_ERI_F2xz_S_Fxyz_S_cc;
  abcd[843] = 4.0E0*I_ERI_Fx2y_S_Fxyz_S_cc;
  abcd[844] = 4.0E0*I_ERI_Fxyz_S_Fxyz_S_cc;
  abcd[845] = 4.0E0*I_ERI_Fx2z_S_Fxyz_S_cc;
  abcd[846] = 4.0E0*I_ERI_F3y_S_Fxyz_S_cc;
  abcd[847] = 4.0E0*I_ERI_F2yz_S_Fxyz_S_cc;
  abcd[848] = 4.0E0*I_ERI_Fy2z_S_Fxyz_S_cc;
  abcd[849] = 4.0E0*I_ERI_F3z_S_Fxyz_S_cc;
  abcd[850] = 4.0E0*I_ERI_F3x_S_F2yz_S_cc-2.0E0*1*I_ERI_F3x_S_Pz_S_c;
  abcd[851] = 4.0E0*I_ERI_F2xy_S_F2yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Pz_S_c;
  abcd[852] = 4.0E0*I_ERI_F2xz_S_F2yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Pz_S_c;
  abcd[853] = 4.0E0*I_ERI_Fx2y_S_F2yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Pz_S_c;
  abcd[854] = 4.0E0*I_ERI_Fxyz_S_F2yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Pz_S_c;
  abcd[855] = 4.0E0*I_ERI_Fx2z_S_F2yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Pz_S_c;
  abcd[856] = 4.0E0*I_ERI_F3y_S_F2yz_S_cc-2.0E0*1*I_ERI_F3y_S_Pz_S_c;
  abcd[857] = 4.0E0*I_ERI_F2yz_S_F2yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Pz_S_c;
  abcd[858] = 4.0E0*I_ERI_Fy2z_S_F2yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Pz_S_c;
  abcd[859] = 4.0E0*I_ERI_F3z_S_F2yz_S_cc-2.0E0*1*I_ERI_F3z_S_Pz_S_c;
  abcd[860] = 4.0E0*I_ERI_F3x_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3x_S_Py_S_c;
  abcd[861] = 4.0E0*I_ERI_F2xy_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Py_S_c;
  abcd[862] = 4.0E0*I_ERI_F2xz_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Py_S_c;
  abcd[863] = 4.0E0*I_ERI_Fx2y_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Py_S_c;
  abcd[864] = 4.0E0*I_ERI_Fxyz_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Py_S_c;
  abcd[865] = 4.0E0*I_ERI_Fx2z_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Py_S_c;
  abcd[866] = 4.0E0*I_ERI_F3y_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3y_S_Py_S_c;
  abcd[867] = 4.0E0*I_ERI_F2yz_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Py_S_c;
  abcd[868] = 4.0E0*I_ERI_Fy2z_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Py_S_c;
  abcd[869] = 4.0E0*I_ERI_F3z_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3z_S_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  abcd[870] = 4.0E0*I_ERI_F3x_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3x_S_Px_S_c;
  abcd[871] = 4.0E0*I_ERI_F2xy_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Px_S_c;
  abcd[872] = 4.0E0*I_ERI_F2xz_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Px_S_c;
  abcd[873] = 4.0E0*I_ERI_Fx2y_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Px_S_c;
  abcd[874] = 4.0E0*I_ERI_Fxyz_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Px_S_c;
  abcd[875] = 4.0E0*I_ERI_Fx2z_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Px_S_c;
  abcd[876] = 4.0E0*I_ERI_F3y_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3y_S_Px_S_c;
  abcd[877] = 4.0E0*I_ERI_F2yz_S_Fx2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Px_S_c;
  abcd[878] = 4.0E0*I_ERI_Fy2z_S_Fx2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Px_S_c;
  abcd[879] = 4.0E0*I_ERI_F3z_S_Fx2z_S_cc-2.0E0*1*I_ERI_F3z_S_Px_S_c;
  abcd[880] = 4.0E0*I_ERI_F3x_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3x_S_Py_S_c;
  abcd[881] = 4.0E0*I_ERI_F2xy_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Py_S_c;
  abcd[882] = 4.0E0*I_ERI_F2xz_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Py_S_c;
  abcd[883] = 4.0E0*I_ERI_Fx2y_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Py_S_c;
  abcd[884] = 4.0E0*I_ERI_Fxyz_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Py_S_c;
  abcd[885] = 4.0E0*I_ERI_Fx2z_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Py_S_c;
  abcd[886] = 4.0E0*I_ERI_F3y_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3y_S_Py_S_c;
  abcd[887] = 4.0E0*I_ERI_F2yz_S_Fy2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Py_S_c;
  abcd[888] = 4.0E0*I_ERI_Fy2z_S_Fy2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Py_S_c;
  abcd[889] = 4.0E0*I_ERI_F3z_S_Fy2z_S_cc-2.0E0*1*I_ERI_F3z_S_Py_S_c;
  abcd[890] = 4.0E0*I_ERI_F3x_S_F3z_S_cc-2.0E0*1*I_ERI_F3x_S_Pz_S_c-2.0E0*2*I_ERI_F3x_S_Pz_S_c;
  abcd[891] = 4.0E0*I_ERI_F2xy_S_F3z_S_cc-2.0E0*1*I_ERI_F2xy_S_Pz_S_c-2.0E0*2*I_ERI_F2xy_S_Pz_S_c;
  abcd[892] = 4.0E0*I_ERI_F2xz_S_F3z_S_cc-2.0E0*1*I_ERI_F2xz_S_Pz_S_c-2.0E0*2*I_ERI_F2xz_S_Pz_S_c;
  abcd[893] = 4.0E0*I_ERI_Fx2y_S_F3z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Pz_S_c-2.0E0*2*I_ERI_Fx2y_S_Pz_S_c;
  abcd[894] = 4.0E0*I_ERI_Fxyz_S_F3z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Pz_S_c-2.0E0*2*I_ERI_Fxyz_S_Pz_S_c;
  abcd[895] = 4.0E0*I_ERI_Fx2z_S_F3z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Pz_S_c-2.0E0*2*I_ERI_Fx2z_S_Pz_S_c;
  abcd[896] = 4.0E0*I_ERI_F3y_S_F3z_S_cc-2.0E0*1*I_ERI_F3y_S_Pz_S_c-2.0E0*2*I_ERI_F3y_S_Pz_S_c;
  abcd[897] = 4.0E0*I_ERI_F2yz_S_F3z_S_cc-2.0E0*1*I_ERI_F2yz_S_Pz_S_c-2.0E0*2*I_ERI_F2yz_S_Pz_S_c;
  abcd[898] = 4.0E0*I_ERI_Fy2z_S_F3z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Pz_S_c-2.0E0*2*I_ERI_Fy2z_S_Pz_S_c;
  abcd[899] = 4.0E0*I_ERI_F3z_S_F3z_S_cc-2.0E0*1*I_ERI_F3z_S_Pz_S_c-2.0E0*2*I_ERI_F3z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[900] = 4.0E0*I_ERI_F3x_S_D2x_Px_cd-2.0E0*1*I_ERI_F3x_S_S_Px_d;
  abcd[901] = 4.0E0*I_ERI_F2xy_S_D2x_Px_cd-2.0E0*1*I_ERI_F2xy_S_S_Px_d;
  abcd[902] = 4.0E0*I_ERI_F2xz_S_D2x_Px_cd-2.0E0*1*I_ERI_F2xz_S_S_Px_d;
  abcd[903] = 4.0E0*I_ERI_Fx2y_S_D2x_Px_cd-2.0E0*1*I_ERI_Fx2y_S_S_Px_d;
  abcd[904] = 4.0E0*I_ERI_Fxyz_S_D2x_Px_cd-2.0E0*1*I_ERI_Fxyz_S_S_Px_d;
  abcd[905] = 4.0E0*I_ERI_Fx2z_S_D2x_Px_cd-2.0E0*1*I_ERI_Fx2z_S_S_Px_d;
  abcd[906] = 4.0E0*I_ERI_F3y_S_D2x_Px_cd-2.0E0*1*I_ERI_F3y_S_S_Px_d;
  abcd[907] = 4.0E0*I_ERI_F2yz_S_D2x_Px_cd-2.0E0*1*I_ERI_F2yz_S_S_Px_d;
  abcd[908] = 4.0E0*I_ERI_Fy2z_S_D2x_Px_cd-2.0E0*1*I_ERI_Fy2z_S_S_Px_d;
  abcd[909] = 4.0E0*I_ERI_F3z_S_D2x_Px_cd-2.0E0*1*I_ERI_F3z_S_S_Px_d;
  abcd[910] = 4.0E0*I_ERI_F3x_S_Dxy_Px_cd;
  abcd[911] = 4.0E0*I_ERI_F2xy_S_Dxy_Px_cd;
  abcd[912] = 4.0E0*I_ERI_F2xz_S_Dxy_Px_cd;
  abcd[913] = 4.0E0*I_ERI_Fx2y_S_Dxy_Px_cd;
  abcd[914] = 4.0E0*I_ERI_Fxyz_S_Dxy_Px_cd;
  abcd[915] = 4.0E0*I_ERI_Fx2z_S_Dxy_Px_cd;
  abcd[916] = 4.0E0*I_ERI_F3y_S_Dxy_Px_cd;
  abcd[917] = 4.0E0*I_ERI_F2yz_S_Dxy_Px_cd;
  abcd[918] = 4.0E0*I_ERI_Fy2z_S_Dxy_Px_cd;
  abcd[919] = 4.0E0*I_ERI_F3z_S_Dxy_Px_cd;
  abcd[920] = 4.0E0*I_ERI_F3x_S_Dxz_Px_cd;
  abcd[921] = 4.0E0*I_ERI_F2xy_S_Dxz_Px_cd;
  abcd[922] = 4.0E0*I_ERI_F2xz_S_Dxz_Px_cd;
  abcd[923] = 4.0E0*I_ERI_Fx2y_S_Dxz_Px_cd;
  abcd[924] = 4.0E0*I_ERI_Fxyz_S_Dxz_Px_cd;
  abcd[925] = 4.0E0*I_ERI_Fx2z_S_Dxz_Px_cd;
  abcd[926] = 4.0E0*I_ERI_F3y_S_Dxz_Px_cd;
  abcd[927] = 4.0E0*I_ERI_F2yz_S_Dxz_Px_cd;
  abcd[928] = 4.0E0*I_ERI_Fy2z_S_Dxz_Px_cd;
  abcd[929] = 4.0E0*I_ERI_F3z_S_Dxz_Px_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[930] = 4.0E0*I_ERI_F3x_S_D2x_Py_cd-2.0E0*1*I_ERI_F3x_S_S_Py_d;
  abcd[931] = 4.0E0*I_ERI_F2xy_S_D2x_Py_cd-2.0E0*1*I_ERI_F2xy_S_S_Py_d;
  abcd[932] = 4.0E0*I_ERI_F2xz_S_D2x_Py_cd-2.0E0*1*I_ERI_F2xz_S_S_Py_d;
  abcd[933] = 4.0E0*I_ERI_Fx2y_S_D2x_Py_cd-2.0E0*1*I_ERI_Fx2y_S_S_Py_d;
  abcd[934] = 4.0E0*I_ERI_Fxyz_S_D2x_Py_cd-2.0E0*1*I_ERI_Fxyz_S_S_Py_d;
  abcd[935] = 4.0E0*I_ERI_Fx2z_S_D2x_Py_cd-2.0E0*1*I_ERI_Fx2z_S_S_Py_d;
  abcd[936] = 4.0E0*I_ERI_F3y_S_D2x_Py_cd-2.0E0*1*I_ERI_F3y_S_S_Py_d;
  abcd[937] = 4.0E0*I_ERI_F2yz_S_D2x_Py_cd-2.0E0*1*I_ERI_F2yz_S_S_Py_d;
  abcd[938] = 4.0E0*I_ERI_Fy2z_S_D2x_Py_cd-2.0E0*1*I_ERI_Fy2z_S_S_Py_d;
  abcd[939] = 4.0E0*I_ERI_F3z_S_D2x_Py_cd-2.0E0*1*I_ERI_F3z_S_S_Py_d;
  abcd[940] = 4.0E0*I_ERI_F3x_S_Dxy_Py_cd;
  abcd[941] = 4.0E0*I_ERI_F2xy_S_Dxy_Py_cd;
  abcd[942] = 4.0E0*I_ERI_F2xz_S_Dxy_Py_cd;
  abcd[943] = 4.0E0*I_ERI_Fx2y_S_Dxy_Py_cd;
  abcd[944] = 4.0E0*I_ERI_Fxyz_S_Dxy_Py_cd;
  abcd[945] = 4.0E0*I_ERI_Fx2z_S_Dxy_Py_cd;
  abcd[946] = 4.0E0*I_ERI_F3y_S_Dxy_Py_cd;
  abcd[947] = 4.0E0*I_ERI_F2yz_S_Dxy_Py_cd;
  abcd[948] = 4.0E0*I_ERI_Fy2z_S_Dxy_Py_cd;
  abcd[949] = 4.0E0*I_ERI_F3z_S_Dxy_Py_cd;
  abcd[950] = 4.0E0*I_ERI_F3x_S_Dxz_Py_cd;
  abcd[951] = 4.0E0*I_ERI_F2xy_S_Dxz_Py_cd;
  abcd[952] = 4.0E0*I_ERI_F2xz_S_Dxz_Py_cd;
  abcd[953] = 4.0E0*I_ERI_Fx2y_S_Dxz_Py_cd;
  abcd[954] = 4.0E0*I_ERI_Fxyz_S_Dxz_Py_cd;
  abcd[955] = 4.0E0*I_ERI_Fx2z_S_Dxz_Py_cd;
  abcd[956] = 4.0E0*I_ERI_F3y_S_Dxz_Py_cd;
  abcd[957] = 4.0E0*I_ERI_F2yz_S_Dxz_Py_cd;
  abcd[958] = 4.0E0*I_ERI_Fy2z_S_Dxz_Py_cd;
  abcd[959] = 4.0E0*I_ERI_F3z_S_Dxz_Py_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[960] = 4.0E0*I_ERI_F3x_S_D2x_Pz_cd-2.0E0*1*I_ERI_F3x_S_S_Pz_d;
  abcd[961] = 4.0E0*I_ERI_F2xy_S_D2x_Pz_cd-2.0E0*1*I_ERI_F2xy_S_S_Pz_d;
  abcd[962] = 4.0E0*I_ERI_F2xz_S_D2x_Pz_cd-2.0E0*1*I_ERI_F2xz_S_S_Pz_d;
  abcd[963] = 4.0E0*I_ERI_Fx2y_S_D2x_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_S_Pz_d;
  abcd[964] = 4.0E0*I_ERI_Fxyz_S_D2x_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_S_Pz_d;
  abcd[965] = 4.0E0*I_ERI_Fx2z_S_D2x_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_S_Pz_d;
  abcd[966] = 4.0E0*I_ERI_F3y_S_D2x_Pz_cd-2.0E0*1*I_ERI_F3y_S_S_Pz_d;
  abcd[967] = 4.0E0*I_ERI_F2yz_S_D2x_Pz_cd-2.0E0*1*I_ERI_F2yz_S_S_Pz_d;
  abcd[968] = 4.0E0*I_ERI_Fy2z_S_D2x_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_S_Pz_d;
  abcd[969] = 4.0E0*I_ERI_F3z_S_D2x_Pz_cd-2.0E0*1*I_ERI_F3z_S_S_Pz_d;
  abcd[970] = 4.0E0*I_ERI_F3x_S_Dxy_Pz_cd;
  abcd[971] = 4.0E0*I_ERI_F2xy_S_Dxy_Pz_cd;
  abcd[972] = 4.0E0*I_ERI_F2xz_S_Dxy_Pz_cd;
  abcd[973] = 4.0E0*I_ERI_Fx2y_S_Dxy_Pz_cd;
  abcd[974] = 4.0E0*I_ERI_Fxyz_S_Dxy_Pz_cd;
  abcd[975] = 4.0E0*I_ERI_Fx2z_S_Dxy_Pz_cd;
  abcd[976] = 4.0E0*I_ERI_F3y_S_Dxy_Pz_cd;
  abcd[977] = 4.0E0*I_ERI_F2yz_S_Dxy_Pz_cd;
  abcd[978] = 4.0E0*I_ERI_Fy2z_S_Dxy_Pz_cd;
  abcd[979] = 4.0E0*I_ERI_F3z_S_Dxy_Pz_cd;
  abcd[980] = 4.0E0*I_ERI_F3x_S_Dxz_Pz_cd;
  abcd[981] = 4.0E0*I_ERI_F2xy_S_Dxz_Pz_cd;
  abcd[982] = 4.0E0*I_ERI_F2xz_S_Dxz_Pz_cd;
  abcd[983] = 4.0E0*I_ERI_Fx2y_S_Dxz_Pz_cd;
  abcd[984] = 4.0E0*I_ERI_Fxyz_S_Dxz_Pz_cd;
  abcd[985] = 4.0E0*I_ERI_Fx2z_S_Dxz_Pz_cd;
  abcd[986] = 4.0E0*I_ERI_F3y_S_Dxz_Pz_cd;
  abcd[987] = 4.0E0*I_ERI_F2yz_S_Dxz_Pz_cd;
  abcd[988] = 4.0E0*I_ERI_Fy2z_S_Dxz_Pz_cd;
  abcd[989] = 4.0E0*I_ERI_F3z_S_Dxz_Pz_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcy_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[990] = 4.0E0*I_ERI_F3x_S_Dxy_Px_cd;
  abcd[991] = 4.0E0*I_ERI_F2xy_S_Dxy_Px_cd;
  abcd[992] = 4.0E0*I_ERI_F2xz_S_Dxy_Px_cd;
  abcd[993] = 4.0E0*I_ERI_Fx2y_S_Dxy_Px_cd;
  abcd[994] = 4.0E0*I_ERI_Fxyz_S_Dxy_Px_cd;
  abcd[995] = 4.0E0*I_ERI_Fx2z_S_Dxy_Px_cd;
  abcd[996] = 4.0E0*I_ERI_F3y_S_Dxy_Px_cd;
  abcd[997] = 4.0E0*I_ERI_F2yz_S_Dxy_Px_cd;
  abcd[998] = 4.0E0*I_ERI_Fy2z_S_Dxy_Px_cd;
  abcd[999] = 4.0E0*I_ERI_F3z_S_Dxy_Px_cd;
  abcd[1000] = 4.0E0*I_ERI_F3x_S_D2y_Px_cd-2.0E0*1*I_ERI_F3x_S_S_Px_d;
  abcd[1001] = 4.0E0*I_ERI_F2xy_S_D2y_Px_cd-2.0E0*1*I_ERI_F2xy_S_S_Px_d;
  abcd[1002] = 4.0E0*I_ERI_F2xz_S_D2y_Px_cd-2.0E0*1*I_ERI_F2xz_S_S_Px_d;
  abcd[1003] = 4.0E0*I_ERI_Fx2y_S_D2y_Px_cd-2.0E0*1*I_ERI_Fx2y_S_S_Px_d;
  abcd[1004] = 4.0E0*I_ERI_Fxyz_S_D2y_Px_cd-2.0E0*1*I_ERI_Fxyz_S_S_Px_d;
  abcd[1005] = 4.0E0*I_ERI_Fx2z_S_D2y_Px_cd-2.0E0*1*I_ERI_Fx2z_S_S_Px_d;
  abcd[1006] = 4.0E0*I_ERI_F3y_S_D2y_Px_cd-2.0E0*1*I_ERI_F3y_S_S_Px_d;
  abcd[1007] = 4.0E0*I_ERI_F2yz_S_D2y_Px_cd-2.0E0*1*I_ERI_F2yz_S_S_Px_d;
  abcd[1008] = 4.0E0*I_ERI_Fy2z_S_D2y_Px_cd-2.0E0*1*I_ERI_Fy2z_S_S_Px_d;
  abcd[1009] = 4.0E0*I_ERI_F3z_S_D2y_Px_cd-2.0E0*1*I_ERI_F3z_S_S_Px_d;
  abcd[1010] = 4.0E0*I_ERI_F3x_S_Dyz_Px_cd;
  abcd[1011] = 4.0E0*I_ERI_F2xy_S_Dyz_Px_cd;
  abcd[1012] = 4.0E0*I_ERI_F2xz_S_Dyz_Px_cd;
  abcd[1013] = 4.0E0*I_ERI_Fx2y_S_Dyz_Px_cd;
  abcd[1014] = 4.0E0*I_ERI_Fxyz_S_Dyz_Px_cd;
  abcd[1015] = 4.0E0*I_ERI_Fx2z_S_Dyz_Px_cd;
  abcd[1016] = 4.0E0*I_ERI_F3y_S_Dyz_Px_cd;
  abcd[1017] = 4.0E0*I_ERI_F2yz_S_Dyz_Px_cd;
  abcd[1018] = 4.0E0*I_ERI_Fy2z_S_Dyz_Px_cd;
  abcd[1019] = 4.0E0*I_ERI_F3z_S_Dyz_Px_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[1020] = 4.0E0*I_ERI_F3x_S_Dxy_Py_cd;
  abcd[1021] = 4.0E0*I_ERI_F2xy_S_Dxy_Py_cd;
  abcd[1022] = 4.0E0*I_ERI_F2xz_S_Dxy_Py_cd;
  abcd[1023] = 4.0E0*I_ERI_Fx2y_S_Dxy_Py_cd;
  abcd[1024] = 4.0E0*I_ERI_Fxyz_S_Dxy_Py_cd;
  abcd[1025] = 4.0E0*I_ERI_Fx2z_S_Dxy_Py_cd;
  abcd[1026] = 4.0E0*I_ERI_F3y_S_Dxy_Py_cd;
  abcd[1027] = 4.0E0*I_ERI_F2yz_S_Dxy_Py_cd;
  abcd[1028] = 4.0E0*I_ERI_Fy2z_S_Dxy_Py_cd;
  abcd[1029] = 4.0E0*I_ERI_F3z_S_Dxy_Py_cd;
  abcd[1030] = 4.0E0*I_ERI_F3x_S_D2y_Py_cd-2.0E0*1*I_ERI_F3x_S_S_Py_d;
  abcd[1031] = 4.0E0*I_ERI_F2xy_S_D2y_Py_cd-2.0E0*1*I_ERI_F2xy_S_S_Py_d;
  abcd[1032] = 4.0E0*I_ERI_F2xz_S_D2y_Py_cd-2.0E0*1*I_ERI_F2xz_S_S_Py_d;
  abcd[1033] = 4.0E0*I_ERI_Fx2y_S_D2y_Py_cd-2.0E0*1*I_ERI_Fx2y_S_S_Py_d;
  abcd[1034] = 4.0E0*I_ERI_Fxyz_S_D2y_Py_cd-2.0E0*1*I_ERI_Fxyz_S_S_Py_d;
  abcd[1035] = 4.0E0*I_ERI_Fx2z_S_D2y_Py_cd-2.0E0*1*I_ERI_Fx2z_S_S_Py_d;
  abcd[1036] = 4.0E0*I_ERI_F3y_S_D2y_Py_cd-2.0E0*1*I_ERI_F3y_S_S_Py_d;
  abcd[1037] = 4.0E0*I_ERI_F2yz_S_D2y_Py_cd-2.0E0*1*I_ERI_F2yz_S_S_Py_d;
  abcd[1038] = 4.0E0*I_ERI_Fy2z_S_D2y_Py_cd-2.0E0*1*I_ERI_Fy2z_S_S_Py_d;
  abcd[1039] = 4.0E0*I_ERI_F3z_S_D2y_Py_cd-2.0E0*1*I_ERI_F3z_S_S_Py_d;
  abcd[1040] = 4.0E0*I_ERI_F3x_S_Dyz_Py_cd;
  abcd[1041] = 4.0E0*I_ERI_F2xy_S_Dyz_Py_cd;
  abcd[1042] = 4.0E0*I_ERI_F2xz_S_Dyz_Py_cd;
  abcd[1043] = 4.0E0*I_ERI_Fx2y_S_Dyz_Py_cd;
  abcd[1044] = 4.0E0*I_ERI_Fxyz_S_Dyz_Py_cd;
  abcd[1045] = 4.0E0*I_ERI_Fx2z_S_Dyz_Py_cd;
  abcd[1046] = 4.0E0*I_ERI_F3y_S_Dyz_Py_cd;
  abcd[1047] = 4.0E0*I_ERI_F2yz_S_Dyz_Py_cd;
  abcd[1048] = 4.0E0*I_ERI_Fy2z_S_Dyz_Py_cd;
  abcd[1049] = 4.0E0*I_ERI_F3z_S_Dyz_Py_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[1050] = 4.0E0*I_ERI_F3x_S_Dxy_Pz_cd;
  abcd[1051] = 4.0E0*I_ERI_F2xy_S_Dxy_Pz_cd;
  abcd[1052] = 4.0E0*I_ERI_F2xz_S_Dxy_Pz_cd;
  abcd[1053] = 4.0E0*I_ERI_Fx2y_S_Dxy_Pz_cd;
  abcd[1054] = 4.0E0*I_ERI_Fxyz_S_Dxy_Pz_cd;
  abcd[1055] = 4.0E0*I_ERI_Fx2z_S_Dxy_Pz_cd;
  abcd[1056] = 4.0E0*I_ERI_F3y_S_Dxy_Pz_cd;
  abcd[1057] = 4.0E0*I_ERI_F2yz_S_Dxy_Pz_cd;
  abcd[1058] = 4.0E0*I_ERI_Fy2z_S_Dxy_Pz_cd;
  abcd[1059] = 4.0E0*I_ERI_F3z_S_Dxy_Pz_cd;
  abcd[1060] = 4.0E0*I_ERI_F3x_S_D2y_Pz_cd-2.0E0*1*I_ERI_F3x_S_S_Pz_d;
  abcd[1061] = 4.0E0*I_ERI_F2xy_S_D2y_Pz_cd-2.0E0*1*I_ERI_F2xy_S_S_Pz_d;
  abcd[1062] = 4.0E0*I_ERI_F2xz_S_D2y_Pz_cd-2.0E0*1*I_ERI_F2xz_S_S_Pz_d;
  abcd[1063] = 4.0E0*I_ERI_Fx2y_S_D2y_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_S_Pz_d;
  abcd[1064] = 4.0E0*I_ERI_Fxyz_S_D2y_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_S_Pz_d;
  abcd[1065] = 4.0E0*I_ERI_Fx2z_S_D2y_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_S_Pz_d;
  abcd[1066] = 4.0E0*I_ERI_F3y_S_D2y_Pz_cd-2.0E0*1*I_ERI_F3y_S_S_Pz_d;
  abcd[1067] = 4.0E0*I_ERI_F2yz_S_D2y_Pz_cd-2.0E0*1*I_ERI_F2yz_S_S_Pz_d;
  abcd[1068] = 4.0E0*I_ERI_Fy2z_S_D2y_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_S_Pz_d;
  abcd[1069] = 4.0E0*I_ERI_F3z_S_D2y_Pz_cd-2.0E0*1*I_ERI_F3z_S_S_Pz_d;
  abcd[1070] = 4.0E0*I_ERI_F3x_S_Dyz_Pz_cd;
  abcd[1071] = 4.0E0*I_ERI_F2xy_S_Dyz_Pz_cd;
  abcd[1072] = 4.0E0*I_ERI_F2xz_S_Dyz_Pz_cd;
  abcd[1073] = 4.0E0*I_ERI_Fx2y_S_Dyz_Pz_cd;
  abcd[1074] = 4.0E0*I_ERI_Fxyz_S_Dyz_Pz_cd;
  abcd[1075] = 4.0E0*I_ERI_Fx2z_S_Dyz_Pz_cd;
  abcd[1076] = 4.0E0*I_ERI_F3y_S_Dyz_Pz_cd;
  abcd[1077] = 4.0E0*I_ERI_F2yz_S_Dyz_Pz_cd;
  abcd[1078] = 4.0E0*I_ERI_Fy2z_S_Dyz_Pz_cd;
  abcd[1079] = 4.0E0*I_ERI_F3z_S_Dyz_Pz_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[1080] = 4.0E0*I_ERI_F3x_S_Dxz_Px_cd;
  abcd[1081] = 4.0E0*I_ERI_F2xy_S_Dxz_Px_cd;
  abcd[1082] = 4.0E0*I_ERI_F2xz_S_Dxz_Px_cd;
  abcd[1083] = 4.0E0*I_ERI_Fx2y_S_Dxz_Px_cd;
  abcd[1084] = 4.0E0*I_ERI_Fxyz_S_Dxz_Px_cd;
  abcd[1085] = 4.0E0*I_ERI_Fx2z_S_Dxz_Px_cd;
  abcd[1086] = 4.0E0*I_ERI_F3y_S_Dxz_Px_cd;
  abcd[1087] = 4.0E0*I_ERI_F2yz_S_Dxz_Px_cd;
  abcd[1088] = 4.0E0*I_ERI_Fy2z_S_Dxz_Px_cd;
  abcd[1089] = 4.0E0*I_ERI_F3z_S_Dxz_Px_cd;
  abcd[1090] = 4.0E0*I_ERI_F3x_S_Dyz_Px_cd;
  abcd[1091] = 4.0E0*I_ERI_F2xy_S_Dyz_Px_cd;
  abcd[1092] = 4.0E0*I_ERI_F2xz_S_Dyz_Px_cd;
  abcd[1093] = 4.0E0*I_ERI_Fx2y_S_Dyz_Px_cd;
  abcd[1094] = 4.0E0*I_ERI_Fxyz_S_Dyz_Px_cd;
  abcd[1095] = 4.0E0*I_ERI_Fx2z_S_Dyz_Px_cd;
  abcd[1096] = 4.0E0*I_ERI_F3y_S_Dyz_Px_cd;
  abcd[1097] = 4.0E0*I_ERI_F2yz_S_Dyz_Px_cd;
  abcd[1098] = 4.0E0*I_ERI_Fy2z_S_Dyz_Px_cd;
  abcd[1099] = 4.0E0*I_ERI_F3z_S_Dyz_Px_cd;
  abcd[1100] = 4.0E0*I_ERI_F3x_S_D2z_Px_cd-2.0E0*1*I_ERI_F3x_S_S_Px_d;
  abcd[1101] = 4.0E0*I_ERI_F2xy_S_D2z_Px_cd-2.0E0*1*I_ERI_F2xy_S_S_Px_d;
  abcd[1102] = 4.0E0*I_ERI_F2xz_S_D2z_Px_cd-2.0E0*1*I_ERI_F2xz_S_S_Px_d;
  abcd[1103] = 4.0E0*I_ERI_Fx2y_S_D2z_Px_cd-2.0E0*1*I_ERI_Fx2y_S_S_Px_d;
  abcd[1104] = 4.0E0*I_ERI_Fxyz_S_D2z_Px_cd-2.0E0*1*I_ERI_Fxyz_S_S_Px_d;
  abcd[1105] = 4.0E0*I_ERI_Fx2z_S_D2z_Px_cd-2.0E0*1*I_ERI_Fx2z_S_S_Px_d;
  abcd[1106] = 4.0E0*I_ERI_F3y_S_D2z_Px_cd-2.0E0*1*I_ERI_F3y_S_S_Px_d;
  abcd[1107] = 4.0E0*I_ERI_F2yz_S_D2z_Px_cd-2.0E0*1*I_ERI_F2yz_S_S_Px_d;
  abcd[1108] = 4.0E0*I_ERI_Fy2z_S_D2z_Px_cd-2.0E0*1*I_ERI_Fy2z_S_S_Px_d;
  abcd[1109] = 4.0E0*I_ERI_F3z_S_D2z_Px_cd-2.0E0*1*I_ERI_F3z_S_S_Px_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[1110] = 4.0E0*I_ERI_F3x_S_Dxz_Py_cd;
  abcd[1111] = 4.0E0*I_ERI_F2xy_S_Dxz_Py_cd;
  abcd[1112] = 4.0E0*I_ERI_F2xz_S_Dxz_Py_cd;
  abcd[1113] = 4.0E0*I_ERI_Fx2y_S_Dxz_Py_cd;
  abcd[1114] = 4.0E0*I_ERI_Fxyz_S_Dxz_Py_cd;
  abcd[1115] = 4.0E0*I_ERI_Fx2z_S_Dxz_Py_cd;
  abcd[1116] = 4.0E0*I_ERI_F3y_S_Dxz_Py_cd;
  abcd[1117] = 4.0E0*I_ERI_F2yz_S_Dxz_Py_cd;
  abcd[1118] = 4.0E0*I_ERI_Fy2z_S_Dxz_Py_cd;
  abcd[1119] = 4.0E0*I_ERI_F3z_S_Dxz_Py_cd;
  abcd[1120] = 4.0E0*I_ERI_F3x_S_Dyz_Py_cd;
  abcd[1121] = 4.0E0*I_ERI_F2xy_S_Dyz_Py_cd;
  abcd[1122] = 4.0E0*I_ERI_F2xz_S_Dyz_Py_cd;
  abcd[1123] = 4.0E0*I_ERI_Fx2y_S_Dyz_Py_cd;
  abcd[1124] = 4.0E0*I_ERI_Fxyz_S_Dyz_Py_cd;
  abcd[1125] = 4.0E0*I_ERI_Fx2z_S_Dyz_Py_cd;
  abcd[1126] = 4.0E0*I_ERI_F3y_S_Dyz_Py_cd;
  abcd[1127] = 4.0E0*I_ERI_F2yz_S_Dyz_Py_cd;
  abcd[1128] = 4.0E0*I_ERI_Fy2z_S_Dyz_Py_cd;
  abcd[1129] = 4.0E0*I_ERI_F3z_S_Dyz_Py_cd;
  abcd[1130] = 4.0E0*I_ERI_F3x_S_D2z_Py_cd-2.0E0*1*I_ERI_F3x_S_S_Py_d;
  abcd[1131] = 4.0E0*I_ERI_F2xy_S_D2z_Py_cd-2.0E0*1*I_ERI_F2xy_S_S_Py_d;
  abcd[1132] = 4.0E0*I_ERI_F2xz_S_D2z_Py_cd-2.0E0*1*I_ERI_F2xz_S_S_Py_d;
  abcd[1133] = 4.0E0*I_ERI_Fx2y_S_D2z_Py_cd-2.0E0*1*I_ERI_Fx2y_S_S_Py_d;
  abcd[1134] = 4.0E0*I_ERI_Fxyz_S_D2z_Py_cd-2.0E0*1*I_ERI_Fxyz_S_S_Py_d;
  abcd[1135] = 4.0E0*I_ERI_Fx2z_S_D2z_Py_cd-2.0E0*1*I_ERI_Fx2z_S_S_Py_d;
  abcd[1136] = 4.0E0*I_ERI_F3y_S_D2z_Py_cd-2.0E0*1*I_ERI_F3y_S_S_Py_d;
  abcd[1137] = 4.0E0*I_ERI_F2yz_S_D2z_Py_cd-2.0E0*1*I_ERI_F2yz_S_S_Py_d;
  abcd[1138] = 4.0E0*I_ERI_Fy2z_S_D2z_Py_cd-2.0E0*1*I_ERI_Fy2z_S_S_Py_d;
  abcd[1139] = 4.0E0*I_ERI_F3z_S_D2z_Py_cd-2.0E0*1*I_ERI_F3z_S_S_Py_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_dcz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_S_P_d
   ************************************************************/
  abcd[1140] = 4.0E0*I_ERI_F3x_S_Dxz_Pz_cd;
  abcd[1141] = 4.0E0*I_ERI_F2xy_S_Dxz_Pz_cd;
  abcd[1142] = 4.0E0*I_ERI_F2xz_S_Dxz_Pz_cd;
  abcd[1143] = 4.0E0*I_ERI_Fx2y_S_Dxz_Pz_cd;
  abcd[1144] = 4.0E0*I_ERI_Fxyz_S_Dxz_Pz_cd;
  abcd[1145] = 4.0E0*I_ERI_Fx2z_S_Dxz_Pz_cd;
  abcd[1146] = 4.0E0*I_ERI_F3y_S_Dxz_Pz_cd;
  abcd[1147] = 4.0E0*I_ERI_F2yz_S_Dxz_Pz_cd;
  abcd[1148] = 4.0E0*I_ERI_Fy2z_S_Dxz_Pz_cd;
  abcd[1149] = 4.0E0*I_ERI_F3z_S_Dxz_Pz_cd;
  abcd[1150] = 4.0E0*I_ERI_F3x_S_Dyz_Pz_cd;
  abcd[1151] = 4.0E0*I_ERI_F2xy_S_Dyz_Pz_cd;
  abcd[1152] = 4.0E0*I_ERI_F2xz_S_Dyz_Pz_cd;
  abcd[1153] = 4.0E0*I_ERI_Fx2y_S_Dyz_Pz_cd;
  abcd[1154] = 4.0E0*I_ERI_Fxyz_S_Dyz_Pz_cd;
  abcd[1155] = 4.0E0*I_ERI_Fx2z_S_Dyz_Pz_cd;
  abcd[1156] = 4.0E0*I_ERI_F3y_S_Dyz_Pz_cd;
  abcd[1157] = 4.0E0*I_ERI_F2yz_S_Dyz_Pz_cd;
  abcd[1158] = 4.0E0*I_ERI_Fy2z_S_Dyz_Pz_cd;
  abcd[1159] = 4.0E0*I_ERI_F3z_S_Dyz_Pz_cd;
  abcd[1160] = 4.0E0*I_ERI_F3x_S_D2z_Pz_cd-2.0E0*1*I_ERI_F3x_S_S_Pz_d;
  abcd[1161] = 4.0E0*I_ERI_F2xy_S_D2z_Pz_cd-2.0E0*1*I_ERI_F2xy_S_S_Pz_d;
  abcd[1162] = 4.0E0*I_ERI_F2xz_S_D2z_Pz_cd-2.0E0*1*I_ERI_F2xz_S_S_Pz_d;
  abcd[1163] = 4.0E0*I_ERI_Fx2y_S_D2z_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_S_Pz_d;
  abcd[1164] = 4.0E0*I_ERI_Fxyz_S_D2z_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_S_Pz_d;
  abcd[1165] = 4.0E0*I_ERI_Fx2z_S_D2z_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_S_Pz_d;
  abcd[1166] = 4.0E0*I_ERI_F3y_S_D2z_Pz_cd-2.0E0*1*I_ERI_F3y_S_S_Pz_d;
  abcd[1167] = 4.0E0*I_ERI_F2yz_S_D2z_Pz_cd-2.0E0*1*I_ERI_F2yz_S_S_Pz_d;
  abcd[1168] = 4.0E0*I_ERI_Fy2z_S_D2z_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_S_Pz_d;
  abcd[1169] = 4.0E0*I_ERI_F3z_S_D2z_Pz_cd-2.0E0*1*I_ERI_F3z_S_S_Pz_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1170] = 4.0E0*I_ERI_F3x_S_Px_D2x_dd-2.0E0*1*I_ERI_F3x_S_Px_S_d;
  abcd[1171] = 4.0E0*I_ERI_F2xy_S_Px_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Px_S_d;
  abcd[1172] = 4.0E0*I_ERI_F2xz_S_Px_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Px_S_d;
  abcd[1173] = 4.0E0*I_ERI_Fx2y_S_Px_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Px_S_d;
  abcd[1174] = 4.0E0*I_ERI_Fxyz_S_Px_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Px_S_d;
  abcd[1175] = 4.0E0*I_ERI_Fx2z_S_Px_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Px_S_d;
  abcd[1176] = 4.0E0*I_ERI_F3y_S_Px_D2x_dd-2.0E0*1*I_ERI_F3y_S_Px_S_d;
  abcd[1177] = 4.0E0*I_ERI_F2yz_S_Px_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Px_S_d;
  abcd[1178] = 4.0E0*I_ERI_Fy2z_S_Px_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Px_S_d;
  abcd[1179] = 4.0E0*I_ERI_F3z_S_Px_D2x_dd-2.0E0*1*I_ERI_F3z_S_Px_S_d;
  abcd[1180] = 4.0E0*I_ERI_F3x_S_Py_D2x_dd-2.0E0*1*I_ERI_F3x_S_Py_S_d;
  abcd[1181] = 4.0E0*I_ERI_F2xy_S_Py_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Py_S_d;
  abcd[1182] = 4.0E0*I_ERI_F2xz_S_Py_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Py_S_d;
  abcd[1183] = 4.0E0*I_ERI_Fx2y_S_Py_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Py_S_d;
  abcd[1184] = 4.0E0*I_ERI_Fxyz_S_Py_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Py_S_d;
  abcd[1185] = 4.0E0*I_ERI_Fx2z_S_Py_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Py_S_d;
  abcd[1186] = 4.0E0*I_ERI_F3y_S_Py_D2x_dd-2.0E0*1*I_ERI_F3y_S_Py_S_d;
  abcd[1187] = 4.0E0*I_ERI_F2yz_S_Py_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Py_S_d;
  abcd[1188] = 4.0E0*I_ERI_Fy2z_S_Py_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Py_S_d;
  abcd[1189] = 4.0E0*I_ERI_F3z_S_Py_D2x_dd-2.0E0*1*I_ERI_F3z_S_Py_S_d;
  abcd[1190] = 4.0E0*I_ERI_F3x_S_Pz_D2x_dd-2.0E0*1*I_ERI_F3x_S_Pz_S_d;
  abcd[1191] = 4.0E0*I_ERI_F2xy_S_Pz_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Pz_S_d;
  abcd[1192] = 4.0E0*I_ERI_F2xz_S_Pz_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Pz_S_d;
  abcd[1193] = 4.0E0*I_ERI_Fx2y_S_Pz_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Pz_S_d;
  abcd[1194] = 4.0E0*I_ERI_Fxyz_S_Pz_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Pz_S_d;
  abcd[1195] = 4.0E0*I_ERI_Fx2z_S_Pz_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Pz_S_d;
  abcd[1196] = 4.0E0*I_ERI_F3y_S_Pz_D2x_dd-2.0E0*1*I_ERI_F3y_S_Pz_S_d;
  abcd[1197] = 4.0E0*I_ERI_F2yz_S_Pz_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Pz_S_d;
  abcd[1198] = 4.0E0*I_ERI_Fy2z_S_Pz_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Pz_S_d;
  abcd[1199] = 4.0E0*I_ERI_F3z_S_Pz_D2x_dd-2.0E0*1*I_ERI_F3z_S_Pz_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1200] = 4.0E0*I_ERI_F3x_S_Px_Dxy_dd;
  abcd[1201] = 4.0E0*I_ERI_F2xy_S_Px_Dxy_dd;
  abcd[1202] = 4.0E0*I_ERI_F2xz_S_Px_Dxy_dd;
  abcd[1203] = 4.0E0*I_ERI_Fx2y_S_Px_Dxy_dd;
  abcd[1204] = 4.0E0*I_ERI_Fxyz_S_Px_Dxy_dd;
  abcd[1205] = 4.0E0*I_ERI_Fx2z_S_Px_Dxy_dd;
  abcd[1206] = 4.0E0*I_ERI_F3y_S_Px_Dxy_dd;
  abcd[1207] = 4.0E0*I_ERI_F2yz_S_Px_Dxy_dd;
  abcd[1208] = 4.0E0*I_ERI_Fy2z_S_Px_Dxy_dd;
  abcd[1209] = 4.0E0*I_ERI_F3z_S_Px_Dxy_dd;
  abcd[1210] = 4.0E0*I_ERI_F3x_S_Py_Dxy_dd;
  abcd[1211] = 4.0E0*I_ERI_F2xy_S_Py_Dxy_dd;
  abcd[1212] = 4.0E0*I_ERI_F2xz_S_Py_Dxy_dd;
  abcd[1213] = 4.0E0*I_ERI_Fx2y_S_Py_Dxy_dd;
  abcd[1214] = 4.0E0*I_ERI_Fxyz_S_Py_Dxy_dd;
  abcd[1215] = 4.0E0*I_ERI_Fx2z_S_Py_Dxy_dd;
  abcd[1216] = 4.0E0*I_ERI_F3y_S_Py_Dxy_dd;
  abcd[1217] = 4.0E0*I_ERI_F2yz_S_Py_Dxy_dd;
  abcd[1218] = 4.0E0*I_ERI_Fy2z_S_Py_Dxy_dd;
  abcd[1219] = 4.0E0*I_ERI_F3z_S_Py_Dxy_dd;
  abcd[1220] = 4.0E0*I_ERI_F3x_S_Pz_Dxy_dd;
  abcd[1221] = 4.0E0*I_ERI_F2xy_S_Pz_Dxy_dd;
  abcd[1222] = 4.0E0*I_ERI_F2xz_S_Pz_Dxy_dd;
  abcd[1223] = 4.0E0*I_ERI_Fx2y_S_Pz_Dxy_dd;
  abcd[1224] = 4.0E0*I_ERI_Fxyz_S_Pz_Dxy_dd;
  abcd[1225] = 4.0E0*I_ERI_Fx2z_S_Pz_Dxy_dd;
  abcd[1226] = 4.0E0*I_ERI_F3y_S_Pz_Dxy_dd;
  abcd[1227] = 4.0E0*I_ERI_F2yz_S_Pz_Dxy_dd;
  abcd[1228] = 4.0E0*I_ERI_Fy2z_S_Pz_Dxy_dd;
  abcd[1229] = 4.0E0*I_ERI_F3z_S_Pz_Dxy_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1230] = 4.0E0*I_ERI_F3x_S_Px_Dxz_dd;
  abcd[1231] = 4.0E0*I_ERI_F2xy_S_Px_Dxz_dd;
  abcd[1232] = 4.0E0*I_ERI_F2xz_S_Px_Dxz_dd;
  abcd[1233] = 4.0E0*I_ERI_Fx2y_S_Px_Dxz_dd;
  abcd[1234] = 4.0E0*I_ERI_Fxyz_S_Px_Dxz_dd;
  abcd[1235] = 4.0E0*I_ERI_Fx2z_S_Px_Dxz_dd;
  abcd[1236] = 4.0E0*I_ERI_F3y_S_Px_Dxz_dd;
  abcd[1237] = 4.0E0*I_ERI_F2yz_S_Px_Dxz_dd;
  abcd[1238] = 4.0E0*I_ERI_Fy2z_S_Px_Dxz_dd;
  abcd[1239] = 4.0E0*I_ERI_F3z_S_Px_Dxz_dd;
  abcd[1240] = 4.0E0*I_ERI_F3x_S_Py_Dxz_dd;
  abcd[1241] = 4.0E0*I_ERI_F2xy_S_Py_Dxz_dd;
  abcd[1242] = 4.0E0*I_ERI_F2xz_S_Py_Dxz_dd;
  abcd[1243] = 4.0E0*I_ERI_Fx2y_S_Py_Dxz_dd;
  abcd[1244] = 4.0E0*I_ERI_Fxyz_S_Py_Dxz_dd;
  abcd[1245] = 4.0E0*I_ERI_Fx2z_S_Py_Dxz_dd;
  abcd[1246] = 4.0E0*I_ERI_F3y_S_Py_Dxz_dd;
  abcd[1247] = 4.0E0*I_ERI_F2yz_S_Py_Dxz_dd;
  abcd[1248] = 4.0E0*I_ERI_Fy2z_S_Py_Dxz_dd;
  abcd[1249] = 4.0E0*I_ERI_F3z_S_Py_Dxz_dd;
  abcd[1250] = 4.0E0*I_ERI_F3x_S_Pz_Dxz_dd;
  abcd[1251] = 4.0E0*I_ERI_F2xy_S_Pz_Dxz_dd;
  abcd[1252] = 4.0E0*I_ERI_F2xz_S_Pz_Dxz_dd;
  abcd[1253] = 4.0E0*I_ERI_Fx2y_S_Pz_Dxz_dd;
  abcd[1254] = 4.0E0*I_ERI_Fxyz_S_Pz_Dxz_dd;
  abcd[1255] = 4.0E0*I_ERI_Fx2z_S_Pz_Dxz_dd;
  abcd[1256] = 4.0E0*I_ERI_F3y_S_Pz_Dxz_dd;
  abcd[1257] = 4.0E0*I_ERI_F2yz_S_Pz_Dxz_dd;
  abcd[1258] = 4.0E0*I_ERI_Fy2z_S_Pz_Dxz_dd;
  abcd[1259] = 4.0E0*I_ERI_F3z_S_Pz_Dxz_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1260] = 4.0E0*I_ERI_F3x_S_Px_D2y_dd-2.0E0*1*I_ERI_F3x_S_Px_S_d;
  abcd[1261] = 4.0E0*I_ERI_F2xy_S_Px_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Px_S_d;
  abcd[1262] = 4.0E0*I_ERI_F2xz_S_Px_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Px_S_d;
  abcd[1263] = 4.0E0*I_ERI_Fx2y_S_Px_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Px_S_d;
  abcd[1264] = 4.0E0*I_ERI_Fxyz_S_Px_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Px_S_d;
  abcd[1265] = 4.0E0*I_ERI_Fx2z_S_Px_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Px_S_d;
  abcd[1266] = 4.0E0*I_ERI_F3y_S_Px_D2y_dd-2.0E0*1*I_ERI_F3y_S_Px_S_d;
  abcd[1267] = 4.0E0*I_ERI_F2yz_S_Px_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Px_S_d;
  abcd[1268] = 4.0E0*I_ERI_Fy2z_S_Px_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Px_S_d;
  abcd[1269] = 4.0E0*I_ERI_F3z_S_Px_D2y_dd-2.0E0*1*I_ERI_F3z_S_Px_S_d;
  abcd[1270] = 4.0E0*I_ERI_F3x_S_Py_D2y_dd-2.0E0*1*I_ERI_F3x_S_Py_S_d;
  abcd[1271] = 4.0E0*I_ERI_F2xy_S_Py_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Py_S_d;
  abcd[1272] = 4.0E0*I_ERI_F2xz_S_Py_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Py_S_d;
  abcd[1273] = 4.0E0*I_ERI_Fx2y_S_Py_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Py_S_d;
  abcd[1274] = 4.0E0*I_ERI_Fxyz_S_Py_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Py_S_d;
  abcd[1275] = 4.0E0*I_ERI_Fx2z_S_Py_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Py_S_d;
  abcd[1276] = 4.0E0*I_ERI_F3y_S_Py_D2y_dd-2.0E0*1*I_ERI_F3y_S_Py_S_d;
  abcd[1277] = 4.0E0*I_ERI_F2yz_S_Py_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Py_S_d;
  abcd[1278] = 4.0E0*I_ERI_Fy2z_S_Py_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Py_S_d;
  abcd[1279] = 4.0E0*I_ERI_F3z_S_Py_D2y_dd-2.0E0*1*I_ERI_F3z_S_Py_S_d;
  abcd[1280] = 4.0E0*I_ERI_F3x_S_Pz_D2y_dd-2.0E0*1*I_ERI_F3x_S_Pz_S_d;
  abcd[1281] = 4.0E0*I_ERI_F2xy_S_Pz_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Pz_S_d;
  abcd[1282] = 4.0E0*I_ERI_F2xz_S_Pz_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Pz_S_d;
  abcd[1283] = 4.0E0*I_ERI_Fx2y_S_Pz_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Pz_S_d;
  abcd[1284] = 4.0E0*I_ERI_Fxyz_S_Pz_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Pz_S_d;
  abcd[1285] = 4.0E0*I_ERI_Fx2z_S_Pz_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Pz_S_d;
  abcd[1286] = 4.0E0*I_ERI_F3y_S_Pz_D2y_dd-2.0E0*1*I_ERI_F3y_S_Pz_S_d;
  abcd[1287] = 4.0E0*I_ERI_F2yz_S_Pz_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Pz_S_d;
  abcd[1288] = 4.0E0*I_ERI_Fy2z_S_Pz_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Pz_S_d;
  abcd[1289] = 4.0E0*I_ERI_F3z_S_Pz_D2y_dd-2.0E0*1*I_ERI_F3z_S_Pz_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1290] = 4.0E0*I_ERI_F3x_S_Px_Dyz_dd;
  abcd[1291] = 4.0E0*I_ERI_F2xy_S_Px_Dyz_dd;
  abcd[1292] = 4.0E0*I_ERI_F2xz_S_Px_Dyz_dd;
  abcd[1293] = 4.0E0*I_ERI_Fx2y_S_Px_Dyz_dd;
  abcd[1294] = 4.0E0*I_ERI_Fxyz_S_Px_Dyz_dd;
  abcd[1295] = 4.0E0*I_ERI_Fx2z_S_Px_Dyz_dd;
  abcd[1296] = 4.0E0*I_ERI_F3y_S_Px_Dyz_dd;
  abcd[1297] = 4.0E0*I_ERI_F2yz_S_Px_Dyz_dd;
  abcd[1298] = 4.0E0*I_ERI_Fy2z_S_Px_Dyz_dd;
  abcd[1299] = 4.0E0*I_ERI_F3z_S_Px_Dyz_dd;
  abcd[1300] = 4.0E0*I_ERI_F3x_S_Py_Dyz_dd;
  abcd[1301] = 4.0E0*I_ERI_F2xy_S_Py_Dyz_dd;
  abcd[1302] = 4.0E0*I_ERI_F2xz_S_Py_Dyz_dd;
  abcd[1303] = 4.0E0*I_ERI_Fx2y_S_Py_Dyz_dd;
  abcd[1304] = 4.0E0*I_ERI_Fxyz_S_Py_Dyz_dd;
  abcd[1305] = 4.0E0*I_ERI_Fx2z_S_Py_Dyz_dd;
  abcd[1306] = 4.0E0*I_ERI_F3y_S_Py_Dyz_dd;
  abcd[1307] = 4.0E0*I_ERI_F2yz_S_Py_Dyz_dd;
  abcd[1308] = 4.0E0*I_ERI_Fy2z_S_Py_Dyz_dd;
  abcd[1309] = 4.0E0*I_ERI_F3z_S_Py_Dyz_dd;
  abcd[1310] = 4.0E0*I_ERI_F3x_S_Pz_Dyz_dd;
  abcd[1311] = 4.0E0*I_ERI_F2xy_S_Pz_Dyz_dd;
  abcd[1312] = 4.0E0*I_ERI_F2xz_S_Pz_Dyz_dd;
  abcd[1313] = 4.0E0*I_ERI_Fx2y_S_Pz_Dyz_dd;
  abcd[1314] = 4.0E0*I_ERI_Fxyz_S_Pz_Dyz_dd;
  abcd[1315] = 4.0E0*I_ERI_Fx2z_S_Pz_Dyz_dd;
  abcd[1316] = 4.0E0*I_ERI_F3y_S_Pz_Dyz_dd;
  abcd[1317] = 4.0E0*I_ERI_F2yz_S_Pz_Dyz_dd;
  abcd[1318] = 4.0E0*I_ERI_Fy2z_S_Pz_Dyz_dd;
  abcd[1319] = 4.0E0*I_ERI_F3z_S_Pz_Dyz_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_ddz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  abcd[1320] = 4.0E0*I_ERI_F3x_S_Px_D2z_dd-2.0E0*1*I_ERI_F3x_S_Px_S_d;
  abcd[1321] = 4.0E0*I_ERI_F2xy_S_Px_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Px_S_d;
  abcd[1322] = 4.0E0*I_ERI_F2xz_S_Px_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Px_S_d;
  abcd[1323] = 4.0E0*I_ERI_Fx2y_S_Px_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Px_S_d;
  abcd[1324] = 4.0E0*I_ERI_Fxyz_S_Px_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Px_S_d;
  abcd[1325] = 4.0E0*I_ERI_Fx2z_S_Px_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Px_S_d;
  abcd[1326] = 4.0E0*I_ERI_F3y_S_Px_D2z_dd-2.0E0*1*I_ERI_F3y_S_Px_S_d;
  abcd[1327] = 4.0E0*I_ERI_F2yz_S_Px_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Px_S_d;
  abcd[1328] = 4.0E0*I_ERI_Fy2z_S_Px_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Px_S_d;
  abcd[1329] = 4.0E0*I_ERI_F3z_S_Px_D2z_dd-2.0E0*1*I_ERI_F3z_S_Px_S_d;
  abcd[1330] = 4.0E0*I_ERI_F3x_S_Py_D2z_dd-2.0E0*1*I_ERI_F3x_S_Py_S_d;
  abcd[1331] = 4.0E0*I_ERI_F2xy_S_Py_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Py_S_d;
  abcd[1332] = 4.0E0*I_ERI_F2xz_S_Py_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Py_S_d;
  abcd[1333] = 4.0E0*I_ERI_Fx2y_S_Py_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Py_S_d;
  abcd[1334] = 4.0E0*I_ERI_Fxyz_S_Py_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Py_S_d;
  abcd[1335] = 4.0E0*I_ERI_Fx2z_S_Py_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Py_S_d;
  abcd[1336] = 4.0E0*I_ERI_F3y_S_Py_D2z_dd-2.0E0*1*I_ERI_F3y_S_Py_S_d;
  abcd[1337] = 4.0E0*I_ERI_F2yz_S_Py_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Py_S_d;
  abcd[1338] = 4.0E0*I_ERI_Fy2z_S_Py_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Py_S_d;
  abcd[1339] = 4.0E0*I_ERI_F3z_S_Py_D2z_dd-2.0E0*1*I_ERI_F3z_S_Py_S_d;
  abcd[1340] = 4.0E0*I_ERI_F3x_S_Pz_D2z_dd-2.0E0*1*I_ERI_F3x_S_Pz_S_d;
  abcd[1341] = 4.0E0*I_ERI_F2xy_S_Pz_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Pz_S_d;
  abcd[1342] = 4.0E0*I_ERI_F2xz_S_Pz_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Pz_S_d;
  abcd[1343] = 4.0E0*I_ERI_Fx2y_S_Pz_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Pz_S_d;
  abcd[1344] = 4.0E0*I_ERI_Fxyz_S_Pz_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Pz_S_d;
  abcd[1345] = 4.0E0*I_ERI_Fx2z_S_Pz_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Pz_S_d;
  abcd[1346] = 4.0E0*I_ERI_F3y_S_Pz_D2z_dd-2.0E0*1*I_ERI_F3y_S_Pz_S_d;
  abcd[1347] = 4.0E0*I_ERI_F2yz_S_Pz_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Pz_S_d;
  abcd[1348] = 4.0E0*I_ERI_Fy2z_S_Pz_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Pz_S_d;
  abcd[1349] = 4.0E0*I_ERI_F3z_S_Pz_D2z_dd-2.0E0*1*I_ERI_F3z_S_Pz_S_d;
}
