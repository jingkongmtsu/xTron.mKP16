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
// BRA1 as redundant position, total RHS integrals evaluated as: 49248
// BRA2 as redundant position, total RHS integrals evaluated as: 44586
// KET1 as redundant position, total RHS integrals evaluated as: 44751
// KET2 as redundant position, total RHS integrals evaluated as: 44535
// the redundant position is: KET2
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
// KET1
// X
// Y
// Z
// ####

void hgp_os_eri_p_sp_p_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_D2x_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000001_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000001_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000001_a = 0.0E0;
  Double I_ERI_S_S_Px_S_C1000001 = 0.0E0;
  Double I_ERI_S_S_Py_S_C1000001 = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1000001 = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001001 = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001001 = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001001 = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001001 = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001001 = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001001 = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000001_c = 0.0E0;
  Double I_ERI_Px_S_S_S_C1000001 = 0.0E0;
  Double I_ERI_Py_S_S_S_C1000001 = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1000001 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001_a = 0.0E0;
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
  Double I_ERI_D2x_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_b = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001_c = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001001 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_b = 0.0E0;

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
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001_a_coefs = ic2*jc2*alpha;
      I_ERI_D2x_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000001_a += SQ_ERI_D_S_P_S_C1000001_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1000001_coefs = ic2*jc2;
      I_ERI_S_S_Px_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001001_coefs = ic2_1*jc2;
      I_ERI_S_Px_Px_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001001 += SQ_ERI_S_P_P_S_C1001001_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_coefs = ic2_1*jc2;
      I_ERI_Px_S_Px_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001 += SQ_ERI_P_S_P_S_C1001001_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_c_coefs = ic2*jc2*gamma;
      I_ERI_Px_S_D2x_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000001_c += SQ_ERI_P_S_D_S_C1000001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1000001_coefs = ic2*jc2;
      I_ERI_Px_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_a_coefs = ic2_1*jc2*alpha;
      I_ERI_F3x_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001_a += SQ_ERI_F_S_P_S_C1001001_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001_b_coefs = ic2*jc2*beta;
      I_ERI_D2x_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000001_b += SQ_ERI_D_S_P_S_C1000001_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_b_coefs = ic2*jc2*beta;
      I_ERI_Px_S_Px_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001_b += SQ_ERI_P_S_P_S_C1000001_b_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_c_coefs = ic2_1*jc2*gamma;
      I_ERI_D2x_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_c += SQ_ERI_D_S_D_S_C1001001_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001_c_coefs = ic2_1*jc2*gamma;
      I_ERI_Px_S_D2x_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001_c += SQ_ERI_P_S_D_S_C1001001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001001_coefs = ic2_1*jc2;
      I_ERI_D2x_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001001 += SQ_ERI_D_S_S_S_C1001001_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001001_coefs = ic2_1*jc2;
      I_ERI_Px_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_b_coefs = ic2_1*jc2*beta;
      I_ERI_F3x_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001_b += SQ_ERI_F_S_P_S_C1001001_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_b_coefs = ic2_1*jc2*beta;
      I_ERI_D2x_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001_b += SQ_ERI_D_S_P_S_C1001001_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_b_coefs = ic2_1*jc2*beta;
      I_ERI_Px_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_P_P_S_S_C1001001
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001001
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001001
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1001001 = I_ERI_D2x_S_S_S_C1001001+ABX*I_ERI_Px_S_S_S_C1001001;
  Double I_ERI_Py_Px_S_S_C1001001 = I_ERI_Dxy_S_S_S_C1001001+ABX*I_ERI_Py_S_S_S_C1001001;
  Double I_ERI_Pz_Px_S_S_C1001001 = I_ERI_Dxz_S_S_S_C1001001+ABX*I_ERI_Pz_S_S_S_C1001001;
  Double I_ERI_Px_Py_S_S_C1001001 = I_ERI_Dxy_S_S_S_C1001001+ABY*I_ERI_Px_S_S_S_C1001001;
  Double I_ERI_Py_Py_S_S_C1001001 = I_ERI_D2y_S_S_S_C1001001+ABY*I_ERI_Py_S_S_S_C1001001;
  Double I_ERI_Pz_Py_S_S_C1001001 = I_ERI_Dyz_S_S_S_C1001001+ABY*I_ERI_Pz_S_S_S_C1001001;
  Double I_ERI_Px_Pz_S_S_C1001001 = I_ERI_Dxz_S_S_S_C1001001+ABZ*I_ERI_Px_S_S_S_C1001001;
  Double I_ERI_Py_Pz_S_S_C1001001 = I_ERI_Dyz_S_S_S_C1001001+ABZ*I_ERI_Py_S_S_S_C1001001;
  Double I_ERI_Pz_Pz_S_S_C1001001 = I_ERI_D2z_S_S_S_C1001001+ABZ*I_ERI_Pz_S_S_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_a
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1001001_a = I_ERI_F3x_S_Px_S_C1001001_a+ABX*I_ERI_D2x_S_Px_S_C1001001_a;
  Double I_ERI_Dxy_Px_Px_S_C1001001_a = I_ERI_F2xy_S_Px_S_C1001001_a+ABX*I_ERI_Dxy_S_Px_S_C1001001_a;
  Double I_ERI_Dxz_Px_Px_S_C1001001_a = I_ERI_F2xz_S_Px_S_C1001001_a+ABX*I_ERI_Dxz_S_Px_S_C1001001_a;
  Double I_ERI_D2y_Px_Px_S_C1001001_a = I_ERI_Fx2y_S_Px_S_C1001001_a+ABX*I_ERI_D2y_S_Px_S_C1001001_a;
  Double I_ERI_Dyz_Px_Px_S_C1001001_a = I_ERI_Fxyz_S_Px_S_C1001001_a+ABX*I_ERI_Dyz_S_Px_S_C1001001_a;
  Double I_ERI_D2z_Px_Px_S_C1001001_a = I_ERI_Fx2z_S_Px_S_C1001001_a+ABX*I_ERI_D2z_S_Px_S_C1001001_a;
  Double I_ERI_D2x_Py_Px_S_C1001001_a = I_ERI_F2xy_S_Px_S_C1001001_a+ABY*I_ERI_D2x_S_Px_S_C1001001_a;
  Double I_ERI_Dxy_Py_Px_S_C1001001_a = I_ERI_Fx2y_S_Px_S_C1001001_a+ABY*I_ERI_Dxy_S_Px_S_C1001001_a;
  Double I_ERI_Dxz_Py_Px_S_C1001001_a = I_ERI_Fxyz_S_Px_S_C1001001_a+ABY*I_ERI_Dxz_S_Px_S_C1001001_a;
  Double I_ERI_D2y_Py_Px_S_C1001001_a = I_ERI_F3y_S_Px_S_C1001001_a+ABY*I_ERI_D2y_S_Px_S_C1001001_a;
  Double I_ERI_Dyz_Py_Px_S_C1001001_a = I_ERI_F2yz_S_Px_S_C1001001_a+ABY*I_ERI_Dyz_S_Px_S_C1001001_a;
  Double I_ERI_D2z_Py_Px_S_C1001001_a = I_ERI_Fy2z_S_Px_S_C1001001_a+ABY*I_ERI_D2z_S_Px_S_C1001001_a;
  Double I_ERI_D2x_Pz_Px_S_C1001001_a = I_ERI_F2xz_S_Px_S_C1001001_a+ABZ*I_ERI_D2x_S_Px_S_C1001001_a;
  Double I_ERI_Dxy_Pz_Px_S_C1001001_a = I_ERI_Fxyz_S_Px_S_C1001001_a+ABZ*I_ERI_Dxy_S_Px_S_C1001001_a;
  Double I_ERI_Dxz_Pz_Px_S_C1001001_a = I_ERI_Fx2z_S_Px_S_C1001001_a+ABZ*I_ERI_Dxz_S_Px_S_C1001001_a;
  Double I_ERI_D2y_Pz_Px_S_C1001001_a = I_ERI_F2yz_S_Px_S_C1001001_a+ABZ*I_ERI_D2y_S_Px_S_C1001001_a;
  Double I_ERI_Dyz_Pz_Px_S_C1001001_a = I_ERI_Fy2z_S_Px_S_C1001001_a+ABZ*I_ERI_Dyz_S_Px_S_C1001001_a;
  Double I_ERI_D2z_Pz_Px_S_C1001001_a = I_ERI_F3z_S_Px_S_C1001001_a+ABZ*I_ERI_D2z_S_Px_S_C1001001_a;
  Double I_ERI_D2x_Px_Py_S_C1001001_a = I_ERI_F3x_S_Py_S_C1001001_a+ABX*I_ERI_D2x_S_Py_S_C1001001_a;
  Double I_ERI_Dxy_Px_Py_S_C1001001_a = I_ERI_F2xy_S_Py_S_C1001001_a+ABX*I_ERI_Dxy_S_Py_S_C1001001_a;
  Double I_ERI_Dxz_Px_Py_S_C1001001_a = I_ERI_F2xz_S_Py_S_C1001001_a+ABX*I_ERI_Dxz_S_Py_S_C1001001_a;
  Double I_ERI_D2y_Px_Py_S_C1001001_a = I_ERI_Fx2y_S_Py_S_C1001001_a+ABX*I_ERI_D2y_S_Py_S_C1001001_a;
  Double I_ERI_Dyz_Px_Py_S_C1001001_a = I_ERI_Fxyz_S_Py_S_C1001001_a+ABX*I_ERI_Dyz_S_Py_S_C1001001_a;
  Double I_ERI_D2z_Px_Py_S_C1001001_a = I_ERI_Fx2z_S_Py_S_C1001001_a+ABX*I_ERI_D2z_S_Py_S_C1001001_a;
  Double I_ERI_D2x_Py_Py_S_C1001001_a = I_ERI_F2xy_S_Py_S_C1001001_a+ABY*I_ERI_D2x_S_Py_S_C1001001_a;
  Double I_ERI_Dxy_Py_Py_S_C1001001_a = I_ERI_Fx2y_S_Py_S_C1001001_a+ABY*I_ERI_Dxy_S_Py_S_C1001001_a;
  Double I_ERI_Dxz_Py_Py_S_C1001001_a = I_ERI_Fxyz_S_Py_S_C1001001_a+ABY*I_ERI_Dxz_S_Py_S_C1001001_a;
  Double I_ERI_D2y_Py_Py_S_C1001001_a = I_ERI_F3y_S_Py_S_C1001001_a+ABY*I_ERI_D2y_S_Py_S_C1001001_a;
  Double I_ERI_Dyz_Py_Py_S_C1001001_a = I_ERI_F2yz_S_Py_S_C1001001_a+ABY*I_ERI_Dyz_S_Py_S_C1001001_a;
  Double I_ERI_D2z_Py_Py_S_C1001001_a = I_ERI_Fy2z_S_Py_S_C1001001_a+ABY*I_ERI_D2z_S_Py_S_C1001001_a;
  Double I_ERI_D2x_Pz_Py_S_C1001001_a = I_ERI_F2xz_S_Py_S_C1001001_a+ABZ*I_ERI_D2x_S_Py_S_C1001001_a;
  Double I_ERI_Dxy_Pz_Py_S_C1001001_a = I_ERI_Fxyz_S_Py_S_C1001001_a+ABZ*I_ERI_Dxy_S_Py_S_C1001001_a;
  Double I_ERI_Dxz_Pz_Py_S_C1001001_a = I_ERI_Fx2z_S_Py_S_C1001001_a+ABZ*I_ERI_Dxz_S_Py_S_C1001001_a;
  Double I_ERI_D2y_Pz_Py_S_C1001001_a = I_ERI_F2yz_S_Py_S_C1001001_a+ABZ*I_ERI_D2y_S_Py_S_C1001001_a;
  Double I_ERI_Dyz_Pz_Py_S_C1001001_a = I_ERI_Fy2z_S_Py_S_C1001001_a+ABZ*I_ERI_Dyz_S_Py_S_C1001001_a;
  Double I_ERI_D2z_Pz_Py_S_C1001001_a = I_ERI_F3z_S_Py_S_C1001001_a+ABZ*I_ERI_D2z_S_Py_S_C1001001_a;
  Double I_ERI_D2x_Px_Pz_S_C1001001_a = I_ERI_F3x_S_Pz_S_C1001001_a+ABX*I_ERI_D2x_S_Pz_S_C1001001_a;
  Double I_ERI_Dxy_Px_Pz_S_C1001001_a = I_ERI_F2xy_S_Pz_S_C1001001_a+ABX*I_ERI_Dxy_S_Pz_S_C1001001_a;
  Double I_ERI_Dxz_Px_Pz_S_C1001001_a = I_ERI_F2xz_S_Pz_S_C1001001_a+ABX*I_ERI_Dxz_S_Pz_S_C1001001_a;
  Double I_ERI_D2y_Px_Pz_S_C1001001_a = I_ERI_Fx2y_S_Pz_S_C1001001_a+ABX*I_ERI_D2y_S_Pz_S_C1001001_a;
  Double I_ERI_Dyz_Px_Pz_S_C1001001_a = I_ERI_Fxyz_S_Pz_S_C1001001_a+ABX*I_ERI_Dyz_S_Pz_S_C1001001_a;
  Double I_ERI_D2z_Px_Pz_S_C1001001_a = I_ERI_Fx2z_S_Pz_S_C1001001_a+ABX*I_ERI_D2z_S_Pz_S_C1001001_a;
  Double I_ERI_D2x_Py_Pz_S_C1001001_a = I_ERI_F2xy_S_Pz_S_C1001001_a+ABY*I_ERI_D2x_S_Pz_S_C1001001_a;
  Double I_ERI_Dxy_Py_Pz_S_C1001001_a = I_ERI_Fx2y_S_Pz_S_C1001001_a+ABY*I_ERI_Dxy_S_Pz_S_C1001001_a;
  Double I_ERI_Dxz_Py_Pz_S_C1001001_a = I_ERI_Fxyz_S_Pz_S_C1001001_a+ABY*I_ERI_Dxz_S_Pz_S_C1001001_a;
  Double I_ERI_D2y_Py_Pz_S_C1001001_a = I_ERI_F3y_S_Pz_S_C1001001_a+ABY*I_ERI_D2y_S_Pz_S_C1001001_a;
  Double I_ERI_Dyz_Py_Pz_S_C1001001_a = I_ERI_F2yz_S_Pz_S_C1001001_a+ABY*I_ERI_Dyz_S_Pz_S_C1001001_a;
  Double I_ERI_D2z_Py_Pz_S_C1001001_a = I_ERI_Fy2z_S_Pz_S_C1001001_a+ABY*I_ERI_D2z_S_Pz_S_C1001001_a;
  Double I_ERI_D2x_Pz_Pz_S_C1001001_a = I_ERI_F2xz_S_Pz_S_C1001001_a+ABZ*I_ERI_D2x_S_Pz_S_C1001001_a;
  Double I_ERI_Dxy_Pz_Pz_S_C1001001_a = I_ERI_Fxyz_S_Pz_S_C1001001_a+ABZ*I_ERI_Dxy_S_Pz_S_C1001001_a;
  Double I_ERI_Dxz_Pz_Pz_S_C1001001_a = I_ERI_Fx2z_S_Pz_S_C1001001_a+ABZ*I_ERI_Dxz_S_Pz_S_C1001001_a;
  Double I_ERI_D2y_Pz_Pz_S_C1001001_a = I_ERI_F2yz_S_Pz_S_C1001001_a+ABZ*I_ERI_D2y_S_Pz_S_C1001001_a;
  Double I_ERI_Dyz_Pz_Pz_S_C1001001_a = I_ERI_Fy2z_S_Pz_S_C1001001_a+ABZ*I_ERI_Dyz_S_Pz_S_C1001001_a;
  Double I_ERI_D2z_Pz_Pz_S_C1001001_a = I_ERI_F3z_S_Pz_S_C1001001_a+ABZ*I_ERI_D2z_S_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1000001_b = I_ERI_D2x_S_Px_S_C1000001_b+ABX*I_ERI_Px_S_Px_S_C1000001_b;
  Double I_ERI_Py_Px_Px_S_C1000001_b = I_ERI_Dxy_S_Px_S_C1000001_b+ABX*I_ERI_Py_S_Px_S_C1000001_b;
  Double I_ERI_Pz_Px_Px_S_C1000001_b = I_ERI_Dxz_S_Px_S_C1000001_b+ABX*I_ERI_Pz_S_Px_S_C1000001_b;
  Double I_ERI_Px_Py_Px_S_C1000001_b = I_ERI_Dxy_S_Px_S_C1000001_b+ABY*I_ERI_Px_S_Px_S_C1000001_b;
  Double I_ERI_Py_Py_Px_S_C1000001_b = I_ERI_D2y_S_Px_S_C1000001_b+ABY*I_ERI_Py_S_Px_S_C1000001_b;
  Double I_ERI_Pz_Py_Px_S_C1000001_b = I_ERI_Dyz_S_Px_S_C1000001_b+ABY*I_ERI_Pz_S_Px_S_C1000001_b;
  Double I_ERI_Px_Pz_Px_S_C1000001_b = I_ERI_Dxz_S_Px_S_C1000001_b+ABZ*I_ERI_Px_S_Px_S_C1000001_b;
  Double I_ERI_Py_Pz_Px_S_C1000001_b = I_ERI_Dyz_S_Px_S_C1000001_b+ABZ*I_ERI_Py_S_Px_S_C1000001_b;
  Double I_ERI_Pz_Pz_Px_S_C1000001_b = I_ERI_D2z_S_Px_S_C1000001_b+ABZ*I_ERI_Pz_S_Px_S_C1000001_b;
  Double I_ERI_Px_Px_Py_S_C1000001_b = I_ERI_D2x_S_Py_S_C1000001_b+ABX*I_ERI_Px_S_Py_S_C1000001_b;
  Double I_ERI_Py_Px_Py_S_C1000001_b = I_ERI_Dxy_S_Py_S_C1000001_b+ABX*I_ERI_Py_S_Py_S_C1000001_b;
  Double I_ERI_Pz_Px_Py_S_C1000001_b = I_ERI_Dxz_S_Py_S_C1000001_b+ABX*I_ERI_Pz_S_Py_S_C1000001_b;
  Double I_ERI_Px_Py_Py_S_C1000001_b = I_ERI_Dxy_S_Py_S_C1000001_b+ABY*I_ERI_Px_S_Py_S_C1000001_b;
  Double I_ERI_Py_Py_Py_S_C1000001_b = I_ERI_D2y_S_Py_S_C1000001_b+ABY*I_ERI_Py_S_Py_S_C1000001_b;
  Double I_ERI_Pz_Py_Py_S_C1000001_b = I_ERI_Dyz_S_Py_S_C1000001_b+ABY*I_ERI_Pz_S_Py_S_C1000001_b;
  Double I_ERI_Px_Pz_Py_S_C1000001_b = I_ERI_Dxz_S_Py_S_C1000001_b+ABZ*I_ERI_Px_S_Py_S_C1000001_b;
  Double I_ERI_Py_Pz_Py_S_C1000001_b = I_ERI_Dyz_S_Py_S_C1000001_b+ABZ*I_ERI_Py_S_Py_S_C1000001_b;
  Double I_ERI_Pz_Pz_Py_S_C1000001_b = I_ERI_D2z_S_Py_S_C1000001_b+ABZ*I_ERI_Pz_S_Py_S_C1000001_b;
  Double I_ERI_Px_Px_Pz_S_C1000001_b = I_ERI_D2x_S_Pz_S_C1000001_b+ABX*I_ERI_Px_S_Pz_S_C1000001_b;
  Double I_ERI_Py_Px_Pz_S_C1000001_b = I_ERI_Dxy_S_Pz_S_C1000001_b+ABX*I_ERI_Py_S_Pz_S_C1000001_b;
  Double I_ERI_Pz_Px_Pz_S_C1000001_b = I_ERI_Dxz_S_Pz_S_C1000001_b+ABX*I_ERI_Pz_S_Pz_S_C1000001_b;
  Double I_ERI_Px_Py_Pz_S_C1000001_b = I_ERI_Dxy_S_Pz_S_C1000001_b+ABY*I_ERI_Px_S_Pz_S_C1000001_b;
  Double I_ERI_Py_Py_Pz_S_C1000001_b = I_ERI_D2y_S_Pz_S_C1000001_b+ABY*I_ERI_Py_S_Pz_S_C1000001_b;
  Double I_ERI_Pz_Py_Pz_S_C1000001_b = I_ERI_Dyz_S_Pz_S_C1000001_b+ABY*I_ERI_Pz_S_Pz_S_C1000001_b;
  Double I_ERI_Px_Pz_Pz_S_C1000001_b = I_ERI_Dxz_S_Pz_S_C1000001_b+ABZ*I_ERI_Px_S_Pz_S_C1000001_b;
  Double I_ERI_Py_Pz_Pz_S_C1000001_b = I_ERI_Dyz_S_Pz_S_C1000001_b+ABZ*I_ERI_Py_S_Pz_S_C1000001_b;
  Double I_ERI_Pz_Pz_Pz_S_C1000001_b = I_ERI_D2z_S_Pz_S_C1000001_b+ABZ*I_ERI_Pz_S_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001001_b = I_ERI_D2x_S_Px_S_C1001001_b+ABX*I_ERI_Px_S_Px_S_C1001001_b;
  Double I_ERI_Py_Px_Px_S_C1001001_b = I_ERI_Dxy_S_Px_S_C1001001_b+ABX*I_ERI_Py_S_Px_S_C1001001_b;
  Double I_ERI_Pz_Px_Px_S_C1001001_b = I_ERI_Dxz_S_Px_S_C1001001_b+ABX*I_ERI_Pz_S_Px_S_C1001001_b;
  Double I_ERI_Px_Py_Px_S_C1001001_b = I_ERI_Dxy_S_Px_S_C1001001_b+ABY*I_ERI_Px_S_Px_S_C1001001_b;
  Double I_ERI_Py_Py_Px_S_C1001001_b = I_ERI_D2y_S_Px_S_C1001001_b+ABY*I_ERI_Py_S_Px_S_C1001001_b;
  Double I_ERI_Pz_Py_Px_S_C1001001_b = I_ERI_Dyz_S_Px_S_C1001001_b+ABY*I_ERI_Pz_S_Px_S_C1001001_b;
  Double I_ERI_Px_Pz_Px_S_C1001001_b = I_ERI_Dxz_S_Px_S_C1001001_b+ABZ*I_ERI_Px_S_Px_S_C1001001_b;
  Double I_ERI_Py_Pz_Px_S_C1001001_b = I_ERI_Dyz_S_Px_S_C1001001_b+ABZ*I_ERI_Py_S_Px_S_C1001001_b;
  Double I_ERI_Pz_Pz_Px_S_C1001001_b = I_ERI_D2z_S_Px_S_C1001001_b+ABZ*I_ERI_Pz_S_Px_S_C1001001_b;
  Double I_ERI_Px_Px_Py_S_C1001001_b = I_ERI_D2x_S_Py_S_C1001001_b+ABX*I_ERI_Px_S_Py_S_C1001001_b;
  Double I_ERI_Py_Px_Py_S_C1001001_b = I_ERI_Dxy_S_Py_S_C1001001_b+ABX*I_ERI_Py_S_Py_S_C1001001_b;
  Double I_ERI_Pz_Px_Py_S_C1001001_b = I_ERI_Dxz_S_Py_S_C1001001_b+ABX*I_ERI_Pz_S_Py_S_C1001001_b;
  Double I_ERI_Px_Py_Py_S_C1001001_b = I_ERI_Dxy_S_Py_S_C1001001_b+ABY*I_ERI_Px_S_Py_S_C1001001_b;
  Double I_ERI_Py_Py_Py_S_C1001001_b = I_ERI_D2y_S_Py_S_C1001001_b+ABY*I_ERI_Py_S_Py_S_C1001001_b;
  Double I_ERI_Pz_Py_Py_S_C1001001_b = I_ERI_Dyz_S_Py_S_C1001001_b+ABY*I_ERI_Pz_S_Py_S_C1001001_b;
  Double I_ERI_Px_Pz_Py_S_C1001001_b = I_ERI_Dxz_S_Py_S_C1001001_b+ABZ*I_ERI_Px_S_Py_S_C1001001_b;
  Double I_ERI_Py_Pz_Py_S_C1001001_b = I_ERI_Dyz_S_Py_S_C1001001_b+ABZ*I_ERI_Py_S_Py_S_C1001001_b;
  Double I_ERI_Pz_Pz_Py_S_C1001001_b = I_ERI_D2z_S_Py_S_C1001001_b+ABZ*I_ERI_Pz_S_Py_S_C1001001_b;
  Double I_ERI_Px_Px_Pz_S_C1001001_b = I_ERI_D2x_S_Pz_S_C1001001_b+ABX*I_ERI_Px_S_Pz_S_C1001001_b;
  Double I_ERI_Py_Px_Pz_S_C1001001_b = I_ERI_Dxy_S_Pz_S_C1001001_b+ABX*I_ERI_Py_S_Pz_S_C1001001_b;
  Double I_ERI_Pz_Px_Pz_S_C1001001_b = I_ERI_Dxz_S_Pz_S_C1001001_b+ABX*I_ERI_Pz_S_Pz_S_C1001001_b;
  Double I_ERI_Px_Py_Pz_S_C1001001_b = I_ERI_Dxy_S_Pz_S_C1001001_b+ABY*I_ERI_Px_S_Pz_S_C1001001_b;
  Double I_ERI_Py_Py_Pz_S_C1001001_b = I_ERI_D2y_S_Pz_S_C1001001_b+ABY*I_ERI_Py_S_Pz_S_C1001001_b;
  Double I_ERI_Pz_Py_Pz_S_C1001001_b = I_ERI_Dyz_S_Pz_S_C1001001_b+ABY*I_ERI_Pz_S_Pz_S_C1001001_b;
  Double I_ERI_Px_Pz_Pz_S_C1001001_b = I_ERI_Dxz_S_Pz_S_C1001001_b+ABZ*I_ERI_Px_S_Pz_S_C1001001_b;
  Double I_ERI_Py_Pz_Pz_S_C1001001_b = I_ERI_Dyz_S_Pz_S_C1001001_b+ABZ*I_ERI_Py_S_Pz_S_C1001001_b;
  Double I_ERI_Pz_Pz_Pz_S_C1001001_b = I_ERI_D2z_S_Pz_S_C1001001_b+ABZ*I_ERI_Pz_S_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1001001_b = I_ERI_F3x_S_Px_S_C1001001_b+ABX*I_ERI_D2x_S_Px_S_C1001001_b;
  Double I_ERI_Dxy_Px_Px_S_C1001001_b = I_ERI_F2xy_S_Px_S_C1001001_b+ABX*I_ERI_Dxy_S_Px_S_C1001001_b;
  Double I_ERI_Dxz_Px_Px_S_C1001001_b = I_ERI_F2xz_S_Px_S_C1001001_b+ABX*I_ERI_Dxz_S_Px_S_C1001001_b;
  Double I_ERI_D2y_Px_Px_S_C1001001_b = I_ERI_Fx2y_S_Px_S_C1001001_b+ABX*I_ERI_D2y_S_Px_S_C1001001_b;
  Double I_ERI_Dyz_Px_Px_S_C1001001_b = I_ERI_Fxyz_S_Px_S_C1001001_b+ABX*I_ERI_Dyz_S_Px_S_C1001001_b;
  Double I_ERI_D2z_Px_Px_S_C1001001_b = I_ERI_Fx2z_S_Px_S_C1001001_b+ABX*I_ERI_D2z_S_Px_S_C1001001_b;
  Double I_ERI_Dxy_Py_Px_S_C1001001_b = I_ERI_Fx2y_S_Px_S_C1001001_b+ABY*I_ERI_Dxy_S_Px_S_C1001001_b;
  Double I_ERI_Dxz_Py_Px_S_C1001001_b = I_ERI_Fxyz_S_Px_S_C1001001_b+ABY*I_ERI_Dxz_S_Px_S_C1001001_b;
  Double I_ERI_D2y_Py_Px_S_C1001001_b = I_ERI_F3y_S_Px_S_C1001001_b+ABY*I_ERI_D2y_S_Px_S_C1001001_b;
  Double I_ERI_Dyz_Py_Px_S_C1001001_b = I_ERI_F2yz_S_Px_S_C1001001_b+ABY*I_ERI_Dyz_S_Px_S_C1001001_b;
  Double I_ERI_D2z_Py_Px_S_C1001001_b = I_ERI_Fy2z_S_Px_S_C1001001_b+ABY*I_ERI_D2z_S_Px_S_C1001001_b;
  Double I_ERI_Dxz_Pz_Px_S_C1001001_b = I_ERI_Fx2z_S_Px_S_C1001001_b+ABZ*I_ERI_Dxz_S_Px_S_C1001001_b;
  Double I_ERI_Dyz_Pz_Px_S_C1001001_b = I_ERI_Fy2z_S_Px_S_C1001001_b+ABZ*I_ERI_Dyz_S_Px_S_C1001001_b;
  Double I_ERI_D2z_Pz_Px_S_C1001001_b = I_ERI_F3z_S_Px_S_C1001001_b+ABZ*I_ERI_D2z_S_Px_S_C1001001_b;
  Double I_ERI_D2x_Px_Py_S_C1001001_b = I_ERI_F3x_S_Py_S_C1001001_b+ABX*I_ERI_D2x_S_Py_S_C1001001_b;
  Double I_ERI_Dxy_Px_Py_S_C1001001_b = I_ERI_F2xy_S_Py_S_C1001001_b+ABX*I_ERI_Dxy_S_Py_S_C1001001_b;
  Double I_ERI_Dxz_Px_Py_S_C1001001_b = I_ERI_F2xz_S_Py_S_C1001001_b+ABX*I_ERI_Dxz_S_Py_S_C1001001_b;
  Double I_ERI_D2y_Px_Py_S_C1001001_b = I_ERI_Fx2y_S_Py_S_C1001001_b+ABX*I_ERI_D2y_S_Py_S_C1001001_b;
  Double I_ERI_Dyz_Px_Py_S_C1001001_b = I_ERI_Fxyz_S_Py_S_C1001001_b+ABX*I_ERI_Dyz_S_Py_S_C1001001_b;
  Double I_ERI_D2z_Px_Py_S_C1001001_b = I_ERI_Fx2z_S_Py_S_C1001001_b+ABX*I_ERI_D2z_S_Py_S_C1001001_b;
  Double I_ERI_Dxy_Py_Py_S_C1001001_b = I_ERI_Fx2y_S_Py_S_C1001001_b+ABY*I_ERI_Dxy_S_Py_S_C1001001_b;
  Double I_ERI_Dxz_Py_Py_S_C1001001_b = I_ERI_Fxyz_S_Py_S_C1001001_b+ABY*I_ERI_Dxz_S_Py_S_C1001001_b;
  Double I_ERI_D2y_Py_Py_S_C1001001_b = I_ERI_F3y_S_Py_S_C1001001_b+ABY*I_ERI_D2y_S_Py_S_C1001001_b;
  Double I_ERI_Dyz_Py_Py_S_C1001001_b = I_ERI_F2yz_S_Py_S_C1001001_b+ABY*I_ERI_Dyz_S_Py_S_C1001001_b;
  Double I_ERI_D2z_Py_Py_S_C1001001_b = I_ERI_Fy2z_S_Py_S_C1001001_b+ABY*I_ERI_D2z_S_Py_S_C1001001_b;
  Double I_ERI_Dxz_Pz_Py_S_C1001001_b = I_ERI_Fx2z_S_Py_S_C1001001_b+ABZ*I_ERI_Dxz_S_Py_S_C1001001_b;
  Double I_ERI_Dyz_Pz_Py_S_C1001001_b = I_ERI_Fy2z_S_Py_S_C1001001_b+ABZ*I_ERI_Dyz_S_Py_S_C1001001_b;
  Double I_ERI_D2z_Pz_Py_S_C1001001_b = I_ERI_F3z_S_Py_S_C1001001_b+ABZ*I_ERI_D2z_S_Py_S_C1001001_b;
  Double I_ERI_D2x_Px_Pz_S_C1001001_b = I_ERI_F3x_S_Pz_S_C1001001_b+ABX*I_ERI_D2x_S_Pz_S_C1001001_b;
  Double I_ERI_Dxy_Px_Pz_S_C1001001_b = I_ERI_F2xy_S_Pz_S_C1001001_b+ABX*I_ERI_Dxy_S_Pz_S_C1001001_b;
  Double I_ERI_Dxz_Px_Pz_S_C1001001_b = I_ERI_F2xz_S_Pz_S_C1001001_b+ABX*I_ERI_Dxz_S_Pz_S_C1001001_b;
  Double I_ERI_D2y_Px_Pz_S_C1001001_b = I_ERI_Fx2y_S_Pz_S_C1001001_b+ABX*I_ERI_D2y_S_Pz_S_C1001001_b;
  Double I_ERI_Dyz_Px_Pz_S_C1001001_b = I_ERI_Fxyz_S_Pz_S_C1001001_b+ABX*I_ERI_Dyz_S_Pz_S_C1001001_b;
  Double I_ERI_D2z_Px_Pz_S_C1001001_b = I_ERI_Fx2z_S_Pz_S_C1001001_b+ABX*I_ERI_D2z_S_Pz_S_C1001001_b;
  Double I_ERI_Dxy_Py_Pz_S_C1001001_b = I_ERI_Fx2y_S_Pz_S_C1001001_b+ABY*I_ERI_Dxy_S_Pz_S_C1001001_b;
  Double I_ERI_Dxz_Py_Pz_S_C1001001_b = I_ERI_Fxyz_S_Pz_S_C1001001_b+ABY*I_ERI_Dxz_S_Pz_S_C1001001_b;
  Double I_ERI_D2y_Py_Pz_S_C1001001_b = I_ERI_F3y_S_Pz_S_C1001001_b+ABY*I_ERI_D2y_S_Pz_S_C1001001_b;
  Double I_ERI_Dyz_Py_Pz_S_C1001001_b = I_ERI_F2yz_S_Pz_S_C1001001_b+ABY*I_ERI_Dyz_S_Pz_S_C1001001_b;
  Double I_ERI_D2z_Py_Pz_S_C1001001_b = I_ERI_Fy2z_S_Pz_S_C1001001_b+ABY*I_ERI_D2z_S_Pz_S_C1001001_b;
  Double I_ERI_Dxz_Pz_Pz_S_C1001001_b = I_ERI_Fx2z_S_Pz_S_C1001001_b+ABZ*I_ERI_Dxz_S_Pz_S_C1001001_b;
  Double I_ERI_Dyz_Pz_Pz_S_C1001001_b = I_ERI_Fy2z_S_Pz_S_C1001001_b+ABZ*I_ERI_Dyz_S_Pz_S_C1001001_b;
  Double I_ERI_D2z_Pz_Pz_S_C1001001_b = I_ERI_F3z_S_Pz_S_C1001001_b+ABZ*I_ERI_D2z_S_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001001_b
   ************************************************************/
  Double I_ERI_Px_D2x_Px_S_C1001001_b = I_ERI_D2x_Px_Px_S_C1001001_b+ABX*I_ERI_Px_Px_Px_S_C1001001_b;
  Double I_ERI_Py_D2x_Px_S_C1001001_b = I_ERI_Dxy_Px_Px_S_C1001001_b+ABX*I_ERI_Py_Px_Px_S_C1001001_b;
  Double I_ERI_Pz_D2x_Px_S_C1001001_b = I_ERI_Dxz_Px_Px_S_C1001001_b+ABX*I_ERI_Pz_Px_Px_S_C1001001_b;
  Double I_ERI_Px_Dxy_Px_S_C1001001_b = I_ERI_Dxy_Px_Px_S_C1001001_b+ABY*I_ERI_Px_Px_Px_S_C1001001_b;
  Double I_ERI_Py_Dxy_Px_S_C1001001_b = I_ERI_D2y_Px_Px_S_C1001001_b+ABY*I_ERI_Py_Px_Px_S_C1001001_b;
  Double I_ERI_Pz_Dxy_Px_S_C1001001_b = I_ERI_Dyz_Px_Px_S_C1001001_b+ABY*I_ERI_Pz_Px_Px_S_C1001001_b;
  Double I_ERI_Px_Dxz_Px_S_C1001001_b = I_ERI_Dxz_Px_Px_S_C1001001_b+ABZ*I_ERI_Px_Px_Px_S_C1001001_b;
  Double I_ERI_Py_Dxz_Px_S_C1001001_b = I_ERI_Dyz_Px_Px_S_C1001001_b+ABZ*I_ERI_Py_Px_Px_S_C1001001_b;
  Double I_ERI_Pz_Dxz_Px_S_C1001001_b = I_ERI_D2z_Px_Px_S_C1001001_b+ABZ*I_ERI_Pz_Px_Px_S_C1001001_b;
  Double I_ERI_Px_D2y_Px_S_C1001001_b = I_ERI_Dxy_Py_Px_S_C1001001_b+ABY*I_ERI_Px_Py_Px_S_C1001001_b;
  Double I_ERI_Py_D2y_Px_S_C1001001_b = I_ERI_D2y_Py_Px_S_C1001001_b+ABY*I_ERI_Py_Py_Px_S_C1001001_b;
  Double I_ERI_Pz_D2y_Px_S_C1001001_b = I_ERI_Dyz_Py_Px_S_C1001001_b+ABY*I_ERI_Pz_Py_Px_S_C1001001_b;
  Double I_ERI_Px_Dyz_Px_S_C1001001_b = I_ERI_Dxz_Py_Px_S_C1001001_b+ABZ*I_ERI_Px_Py_Px_S_C1001001_b;
  Double I_ERI_Py_Dyz_Px_S_C1001001_b = I_ERI_Dyz_Py_Px_S_C1001001_b+ABZ*I_ERI_Py_Py_Px_S_C1001001_b;
  Double I_ERI_Pz_Dyz_Px_S_C1001001_b = I_ERI_D2z_Py_Px_S_C1001001_b+ABZ*I_ERI_Pz_Py_Px_S_C1001001_b;
  Double I_ERI_Px_D2z_Px_S_C1001001_b = I_ERI_Dxz_Pz_Px_S_C1001001_b+ABZ*I_ERI_Px_Pz_Px_S_C1001001_b;
  Double I_ERI_Py_D2z_Px_S_C1001001_b = I_ERI_Dyz_Pz_Px_S_C1001001_b+ABZ*I_ERI_Py_Pz_Px_S_C1001001_b;
  Double I_ERI_Pz_D2z_Px_S_C1001001_b = I_ERI_D2z_Pz_Px_S_C1001001_b+ABZ*I_ERI_Pz_Pz_Px_S_C1001001_b;
  Double I_ERI_Px_D2x_Py_S_C1001001_b = I_ERI_D2x_Px_Py_S_C1001001_b+ABX*I_ERI_Px_Px_Py_S_C1001001_b;
  Double I_ERI_Py_D2x_Py_S_C1001001_b = I_ERI_Dxy_Px_Py_S_C1001001_b+ABX*I_ERI_Py_Px_Py_S_C1001001_b;
  Double I_ERI_Pz_D2x_Py_S_C1001001_b = I_ERI_Dxz_Px_Py_S_C1001001_b+ABX*I_ERI_Pz_Px_Py_S_C1001001_b;
  Double I_ERI_Px_Dxy_Py_S_C1001001_b = I_ERI_Dxy_Px_Py_S_C1001001_b+ABY*I_ERI_Px_Px_Py_S_C1001001_b;
  Double I_ERI_Py_Dxy_Py_S_C1001001_b = I_ERI_D2y_Px_Py_S_C1001001_b+ABY*I_ERI_Py_Px_Py_S_C1001001_b;
  Double I_ERI_Pz_Dxy_Py_S_C1001001_b = I_ERI_Dyz_Px_Py_S_C1001001_b+ABY*I_ERI_Pz_Px_Py_S_C1001001_b;
  Double I_ERI_Px_Dxz_Py_S_C1001001_b = I_ERI_Dxz_Px_Py_S_C1001001_b+ABZ*I_ERI_Px_Px_Py_S_C1001001_b;
  Double I_ERI_Py_Dxz_Py_S_C1001001_b = I_ERI_Dyz_Px_Py_S_C1001001_b+ABZ*I_ERI_Py_Px_Py_S_C1001001_b;
  Double I_ERI_Pz_Dxz_Py_S_C1001001_b = I_ERI_D2z_Px_Py_S_C1001001_b+ABZ*I_ERI_Pz_Px_Py_S_C1001001_b;
  Double I_ERI_Px_D2y_Py_S_C1001001_b = I_ERI_Dxy_Py_Py_S_C1001001_b+ABY*I_ERI_Px_Py_Py_S_C1001001_b;
  Double I_ERI_Py_D2y_Py_S_C1001001_b = I_ERI_D2y_Py_Py_S_C1001001_b+ABY*I_ERI_Py_Py_Py_S_C1001001_b;
  Double I_ERI_Pz_D2y_Py_S_C1001001_b = I_ERI_Dyz_Py_Py_S_C1001001_b+ABY*I_ERI_Pz_Py_Py_S_C1001001_b;
  Double I_ERI_Px_Dyz_Py_S_C1001001_b = I_ERI_Dxz_Py_Py_S_C1001001_b+ABZ*I_ERI_Px_Py_Py_S_C1001001_b;
  Double I_ERI_Py_Dyz_Py_S_C1001001_b = I_ERI_Dyz_Py_Py_S_C1001001_b+ABZ*I_ERI_Py_Py_Py_S_C1001001_b;
  Double I_ERI_Pz_Dyz_Py_S_C1001001_b = I_ERI_D2z_Py_Py_S_C1001001_b+ABZ*I_ERI_Pz_Py_Py_S_C1001001_b;
  Double I_ERI_Px_D2z_Py_S_C1001001_b = I_ERI_Dxz_Pz_Py_S_C1001001_b+ABZ*I_ERI_Px_Pz_Py_S_C1001001_b;
  Double I_ERI_Py_D2z_Py_S_C1001001_b = I_ERI_Dyz_Pz_Py_S_C1001001_b+ABZ*I_ERI_Py_Pz_Py_S_C1001001_b;
  Double I_ERI_Pz_D2z_Py_S_C1001001_b = I_ERI_D2z_Pz_Py_S_C1001001_b+ABZ*I_ERI_Pz_Pz_Py_S_C1001001_b;
  Double I_ERI_Px_D2x_Pz_S_C1001001_b = I_ERI_D2x_Px_Pz_S_C1001001_b+ABX*I_ERI_Px_Px_Pz_S_C1001001_b;
  Double I_ERI_Py_D2x_Pz_S_C1001001_b = I_ERI_Dxy_Px_Pz_S_C1001001_b+ABX*I_ERI_Py_Px_Pz_S_C1001001_b;
  Double I_ERI_Pz_D2x_Pz_S_C1001001_b = I_ERI_Dxz_Px_Pz_S_C1001001_b+ABX*I_ERI_Pz_Px_Pz_S_C1001001_b;
  Double I_ERI_Px_Dxy_Pz_S_C1001001_b = I_ERI_Dxy_Px_Pz_S_C1001001_b+ABY*I_ERI_Px_Px_Pz_S_C1001001_b;
  Double I_ERI_Py_Dxy_Pz_S_C1001001_b = I_ERI_D2y_Px_Pz_S_C1001001_b+ABY*I_ERI_Py_Px_Pz_S_C1001001_b;
  Double I_ERI_Pz_Dxy_Pz_S_C1001001_b = I_ERI_Dyz_Px_Pz_S_C1001001_b+ABY*I_ERI_Pz_Px_Pz_S_C1001001_b;
  Double I_ERI_Px_Dxz_Pz_S_C1001001_b = I_ERI_Dxz_Px_Pz_S_C1001001_b+ABZ*I_ERI_Px_Px_Pz_S_C1001001_b;
  Double I_ERI_Py_Dxz_Pz_S_C1001001_b = I_ERI_Dyz_Px_Pz_S_C1001001_b+ABZ*I_ERI_Py_Px_Pz_S_C1001001_b;
  Double I_ERI_Pz_Dxz_Pz_S_C1001001_b = I_ERI_D2z_Px_Pz_S_C1001001_b+ABZ*I_ERI_Pz_Px_Pz_S_C1001001_b;
  Double I_ERI_Px_D2y_Pz_S_C1001001_b = I_ERI_Dxy_Py_Pz_S_C1001001_b+ABY*I_ERI_Px_Py_Pz_S_C1001001_b;
  Double I_ERI_Py_D2y_Pz_S_C1001001_b = I_ERI_D2y_Py_Pz_S_C1001001_b+ABY*I_ERI_Py_Py_Pz_S_C1001001_b;
  Double I_ERI_Pz_D2y_Pz_S_C1001001_b = I_ERI_Dyz_Py_Pz_S_C1001001_b+ABY*I_ERI_Pz_Py_Pz_S_C1001001_b;
  Double I_ERI_Px_Dyz_Pz_S_C1001001_b = I_ERI_Dxz_Py_Pz_S_C1001001_b+ABZ*I_ERI_Px_Py_Pz_S_C1001001_b;
  Double I_ERI_Py_Dyz_Pz_S_C1001001_b = I_ERI_Dyz_Py_Pz_S_C1001001_b+ABZ*I_ERI_Py_Py_Pz_S_C1001001_b;
  Double I_ERI_Pz_Dyz_Pz_S_C1001001_b = I_ERI_D2z_Py_Pz_S_C1001001_b+ABZ*I_ERI_Pz_Py_Pz_S_C1001001_b;
  Double I_ERI_Px_D2z_Pz_S_C1001001_b = I_ERI_Dxz_Pz_Pz_S_C1001001_b+ABZ*I_ERI_Px_Pz_Pz_S_C1001001_b;
  Double I_ERI_Py_D2z_Pz_S_C1001001_b = I_ERI_Dyz_Pz_Pz_S_C1001001_b+ABZ*I_ERI_Py_Pz_Pz_S_C1001001_b;
  Double I_ERI_Pz_D2z_Pz_S_C1001001_b = I_ERI_D2z_Pz_Pz_S_C1001001_b+ABZ*I_ERI_Pz_Pz_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001_c
   ************************************************************/
  Double I_ERI_Px_Px_D2x_S_C1001001_c = I_ERI_D2x_S_D2x_S_C1001001_c+ABX*I_ERI_Px_S_D2x_S_C1001001_c;
  Double I_ERI_Py_Px_D2x_S_C1001001_c = I_ERI_Dxy_S_D2x_S_C1001001_c+ABX*I_ERI_Py_S_D2x_S_C1001001_c;
  Double I_ERI_Pz_Px_D2x_S_C1001001_c = I_ERI_Dxz_S_D2x_S_C1001001_c+ABX*I_ERI_Pz_S_D2x_S_C1001001_c;
  Double I_ERI_Px_Py_D2x_S_C1001001_c = I_ERI_Dxy_S_D2x_S_C1001001_c+ABY*I_ERI_Px_S_D2x_S_C1001001_c;
  Double I_ERI_Py_Py_D2x_S_C1001001_c = I_ERI_D2y_S_D2x_S_C1001001_c+ABY*I_ERI_Py_S_D2x_S_C1001001_c;
  Double I_ERI_Pz_Py_D2x_S_C1001001_c = I_ERI_Dyz_S_D2x_S_C1001001_c+ABY*I_ERI_Pz_S_D2x_S_C1001001_c;
  Double I_ERI_Px_Pz_D2x_S_C1001001_c = I_ERI_Dxz_S_D2x_S_C1001001_c+ABZ*I_ERI_Px_S_D2x_S_C1001001_c;
  Double I_ERI_Py_Pz_D2x_S_C1001001_c = I_ERI_Dyz_S_D2x_S_C1001001_c+ABZ*I_ERI_Py_S_D2x_S_C1001001_c;
  Double I_ERI_Pz_Pz_D2x_S_C1001001_c = I_ERI_D2z_S_D2x_S_C1001001_c+ABZ*I_ERI_Pz_S_D2x_S_C1001001_c;
  Double I_ERI_Px_Px_Dxy_S_C1001001_c = I_ERI_D2x_S_Dxy_S_C1001001_c+ABX*I_ERI_Px_S_Dxy_S_C1001001_c;
  Double I_ERI_Py_Px_Dxy_S_C1001001_c = I_ERI_Dxy_S_Dxy_S_C1001001_c+ABX*I_ERI_Py_S_Dxy_S_C1001001_c;
  Double I_ERI_Pz_Px_Dxy_S_C1001001_c = I_ERI_Dxz_S_Dxy_S_C1001001_c+ABX*I_ERI_Pz_S_Dxy_S_C1001001_c;
  Double I_ERI_Px_Py_Dxy_S_C1001001_c = I_ERI_Dxy_S_Dxy_S_C1001001_c+ABY*I_ERI_Px_S_Dxy_S_C1001001_c;
  Double I_ERI_Py_Py_Dxy_S_C1001001_c = I_ERI_D2y_S_Dxy_S_C1001001_c+ABY*I_ERI_Py_S_Dxy_S_C1001001_c;
  Double I_ERI_Pz_Py_Dxy_S_C1001001_c = I_ERI_Dyz_S_Dxy_S_C1001001_c+ABY*I_ERI_Pz_S_Dxy_S_C1001001_c;
  Double I_ERI_Px_Pz_Dxy_S_C1001001_c = I_ERI_Dxz_S_Dxy_S_C1001001_c+ABZ*I_ERI_Px_S_Dxy_S_C1001001_c;
  Double I_ERI_Py_Pz_Dxy_S_C1001001_c = I_ERI_Dyz_S_Dxy_S_C1001001_c+ABZ*I_ERI_Py_S_Dxy_S_C1001001_c;
  Double I_ERI_Pz_Pz_Dxy_S_C1001001_c = I_ERI_D2z_S_Dxy_S_C1001001_c+ABZ*I_ERI_Pz_S_Dxy_S_C1001001_c;
  Double I_ERI_Px_Px_Dxz_S_C1001001_c = I_ERI_D2x_S_Dxz_S_C1001001_c+ABX*I_ERI_Px_S_Dxz_S_C1001001_c;
  Double I_ERI_Py_Px_Dxz_S_C1001001_c = I_ERI_Dxy_S_Dxz_S_C1001001_c+ABX*I_ERI_Py_S_Dxz_S_C1001001_c;
  Double I_ERI_Pz_Px_Dxz_S_C1001001_c = I_ERI_Dxz_S_Dxz_S_C1001001_c+ABX*I_ERI_Pz_S_Dxz_S_C1001001_c;
  Double I_ERI_Px_Py_Dxz_S_C1001001_c = I_ERI_Dxy_S_Dxz_S_C1001001_c+ABY*I_ERI_Px_S_Dxz_S_C1001001_c;
  Double I_ERI_Py_Py_Dxz_S_C1001001_c = I_ERI_D2y_S_Dxz_S_C1001001_c+ABY*I_ERI_Py_S_Dxz_S_C1001001_c;
  Double I_ERI_Pz_Py_Dxz_S_C1001001_c = I_ERI_Dyz_S_Dxz_S_C1001001_c+ABY*I_ERI_Pz_S_Dxz_S_C1001001_c;
  Double I_ERI_Px_Pz_Dxz_S_C1001001_c = I_ERI_Dxz_S_Dxz_S_C1001001_c+ABZ*I_ERI_Px_S_Dxz_S_C1001001_c;
  Double I_ERI_Py_Pz_Dxz_S_C1001001_c = I_ERI_Dyz_S_Dxz_S_C1001001_c+ABZ*I_ERI_Py_S_Dxz_S_C1001001_c;
  Double I_ERI_Pz_Pz_Dxz_S_C1001001_c = I_ERI_D2z_S_Dxz_S_C1001001_c+ABZ*I_ERI_Pz_S_Dxz_S_C1001001_c;
  Double I_ERI_Px_Px_D2y_S_C1001001_c = I_ERI_D2x_S_D2y_S_C1001001_c+ABX*I_ERI_Px_S_D2y_S_C1001001_c;
  Double I_ERI_Py_Px_D2y_S_C1001001_c = I_ERI_Dxy_S_D2y_S_C1001001_c+ABX*I_ERI_Py_S_D2y_S_C1001001_c;
  Double I_ERI_Pz_Px_D2y_S_C1001001_c = I_ERI_Dxz_S_D2y_S_C1001001_c+ABX*I_ERI_Pz_S_D2y_S_C1001001_c;
  Double I_ERI_Px_Py_D2y_S_C1001001_c = I_ERI_Dxy_S_D2y_S_C1001001_c+ABY*I_ERI_Px_S_D2y_S_C1001001_c;
  Double I_ERI_Py_Py_D2y_S_C1001001_c = I_ERI_D2y_S_D2y_S_C1001001_c+ABY*I_ERI_Py_S_D2y_S_C1001001_c;
  Double I_ERI_Pz_Py_D2y_S_C1001001_c = I_ERI_Dyz_S_D2y_S_C1001001_c+ABY*I_ERI_Pz_S_D2y_S_C1001001_c;
  Double I_ERI_Px_Pz_D2y_S_C1001001_c = I_ERI_Dxz_S_D2y_S_C1001001_c+ABZ*I_ERI_Px_S_D2y_S_C1001001_c;
  Double I_ERI_Py_Pz_D2y_S_C1001001_c = I_ERI_Dyz_S_D2y_S_C1001001_c+ABZ*I_ERI_Py_S_D2y_S_C1001001_c;
  Double I_ERI_Pz_Pz_D2y_S_C1001001_c = I_ERI_D2z_S_D2y_S_C1001001_c+ABZ*I_ERI_Pz_S_D2y_S_C1001001_c;
  Double I_ERI_Px_Px_Dyz_S_C1001001_c = I_ERI_D2x_S_Dyz_S_C1001001_c+ABX*I_ERI_Px_S_Dyz_S_C1001001_c;
  Double I_ERI_Py_Px_Dyz_S_C1001001_c = I_ERI_Dxy_S_Dyz_S_C1001001_c+ABX*I_ERI_Py_S_Dyz_S_C1001001_c;
  Double I_ERI_Pz_Px_Dyz_S_C1001001_c = I_ERI_Dxz_S_Dyz_S_C1001001_c+ABX*I_ERI_Pz_S_Dyz_S_C1001001_c;
  Double I_ERI_Px_Py_Dyz_S_C1001001_c = I_ERI_Dxy_S_Dyz_S_C1001001_c+ABY*I_ERI_Px_S_Dyz_S_C1001001_c;
  Double I_ERI_Py_Py_Dyz_S_C1001001_c = I_ERI_D2y_S_Dyz_S_C1001001_c+ABY*I_ERI_Py_S_Dyz_S_C1001001_c;
  Double I_ERI_Pz_Py_Dyz_S_C1001001_c = I_ERI_Dyz_S_Dyz_S_C1001001_c+ABY*I_ERI_Pz_S_Dyz_S_C1001001_c;
  Double I_ERI_Px_Pz_Dyz_S_C1001001_c = I_ERI_Dxz_S_Dyz_S_C1001001_c+ABZ*I_ERI_Px_S_Dyz_S_C1001001_c;
  Double I_ERI_Py_Pz_Dyz_S_C1001001_c = I_ERI_Dyz_S_Dyz_S_C1001001_c+ABZ*I_ERI_Py_S_Dyz_S_C1001001_c;
  Double I_ERI_Pz_Pz_Dyz_S_C1001001_c = I_ERI_D2z_S_Dyz_S_C1001001_c+ABZ*I_ERI_Pz_S_Dyz_S_C1001001_c;
  Double I_ERI_Px_Px_D2z_S_C1001001_c = I_ERI_D2x_S_D2z_S_C1001001_c+ABX*I_ERI_Px_S_D2z_S_C1001001_c;
  Double I_ERI_Py_Px_D2z_S_C1001001_c = I_ERI_Dxy_S_D2z_S_C1001001_c+ABX*I_ERI_Py_S_D2z_S_C1001001_c;
  Double I_ERI_Pz_Px_D2z_S_C1001001_c = I_ERI_Dxz_S_D2z_S_C1001001_c+ABX*I_ERI_Pz_S_D2z_S_C1001001_c;
  Double I_ERI_Px_Py_D2z_S_C1001001_c = I_ERI_Dxy_S_D2z_S_C1001001_c+ABY*I_ERI_Px_S_D2z_S_C1001001_c;
  Double I_ERI_Py_Py_D2z_S_C1001001_c = I_ERI_D2y_S_D2z_S_C1001001_c+ABY*I_ERI_Py_S_D2z_S_C1001001_c;
  Double I_ERI_Pz_Py_D2z_S_C1001001_c = I_ERI_Dyz_S_D2z_S_C1001001_c+ABY*I_ERI_Pz_S_D2z_S_C1001001_c;
  Double I_ERI_Px_Pz_D2z_S_C1001001_c = I_ERI_Dxz_S_D2z_S_C1001001_c+ABZ*I_ERI_Px_S_D2z_S_C1001001_c;
  Double I_ERI_Py_Pz_D2z_S_C1001001_c = I_ERI_Dyz_S_D2z_S_C1001001_c+ABZ*I_ERI_Py_S_D2z_S_C1001001_c;
  Double I_ERI_Pz_Pz_D2z_S_C1001001_c = I_ERI_D2z_S_D2z_S_C1001001_c+ABZ*I_ERI_Pz_S_D2z_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_D2x_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[1] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[2] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[12] = 2.0E0*I_ERI_D2x_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[13] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[14] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[24] = 2.0E0*I_ERI_D2x_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[25] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[26] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[3] = 2.0E0*I_ERI_D2x_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[4] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[5] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[6] = 2.0E0*I_ERI_D2x_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[7] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[8] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[9] = 2.0E0*I_ERI_D2x_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[10] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[11] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[15] = 2.0E0*I_ERI_D2x_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[16] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[17] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[18] = 2.0E0*I_ERI_D2x_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[19] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[20] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[21] = 2.0E0*I_ERI_D2x_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[22] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[23] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[27] = 2.0E0*I_ERI_D2x_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[28] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[29] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[30] = 2.0E0*I_ERI_D2x_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[31] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[32] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[33] = 2.0E0*I_ERI_D2x_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[34] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[35] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[36] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[37] = 2.0E0*I_ERI_D2y_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[38] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[48] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[49] = 2.0E0*I_ERI_D2y_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[50] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[60] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[61] = 2.0E0*I_ERI_D2y_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[62] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[39] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[40] = 2.0E0*I_ERI_D2y_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[41] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[42] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[43] = 2.0E0*I_ERI_D2y_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[44] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[45] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[46] = 2.0E0*I_ERI_D2y_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[47] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[51] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[52] = 2.0E0*I_ERI_D2y_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[53] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[54] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[55] = 2.0E0*I_ERI_D2y_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[56] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[57] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[58] = 2.0E0*I_ERI_D2y_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[59] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[63] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[64] = 2.0E0*I_ERI_D2y_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[65] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[66] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[67] = 2.0E0*I_ERI_D2y_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[68] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[69] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[70] = 2.0E0*I_ERI_D2y_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[71] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[72] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[73] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[74] = 2.0E0*I_ERI_D2z_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[84] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[85] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[86] = 2.0E0*I_ERI_D2z_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[96] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;
  abcd[97] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;
  abcd[98] = 2.0E0*I_ERI_D2z_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[75] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[76] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[77] = 2.0E0*I_ERI_D2z_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[78] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[79] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[80] = 2.0E0*I_ERI_D2z_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[81] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[82] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[83] = 2.0E0*I_ERI_D2z_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[87] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[88] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[89] = 2.0E0*I_ERI_D2z_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[90] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[91] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[92] = 2.0E0*I_ERI_D2z_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[93] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[94] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[95] = 2.0E0*I_ERI_D2z_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[99] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[100] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[101] = 2.0E0*I_ERI_D2z_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[102] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[103] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[104] = 2.0E0*I_ERI_D2z_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[105] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;
  abcd[106] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;
  abcd[107] = 2.0E0*I_ERI_D2z_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[108] = 2.0E0*I_ERI_Px_Px_Px_S_C1000001_b;
  abcd[109] = 2.0E0*I_ERI_Py_Px_Px_S_C1000001_b;
  abcd[110] = 2.0E0*I_ERI_Pz_Px_Px_S_C1000001_b;
  abcd[120] = 2.0E0*I_ERI_Px_Px_Py_S_C1000001_b;
  abcd[121] = 2.0E0*I_ERI_Py_Px_Py_S_C1000001_b;
  abcd[122] = 2.0E0*I_ERI_Pz_Px_Py_S_C1000001_b;
  abcd[132] = 2.0E0*I_ERI_Px_Px_Pz_S_C1000001_b;
  abcd[133] = 2.0E0*I_ERI_Py_Px_Pz_S_C1000001_b;
  abcd[134] = 2.0E0*I_ERI_Pz_Px_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[111] = 2.0E0*I_ERI_Px_D2x_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[112] = 2.0E0*I_ERI_Py_D2x_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[113] = 2.0E0*I_ERI_Pz_D2x_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[114] = 2.0E0*I_ERI_Px_Dxy_Px_S_C1001001_b;
  abcd[115] = 2.0E0*I_ERI_Py_Dxy_Px_S_C1001001_b;
  abcd[116] = 2.0E0*I_ERI_Pz_Dxy_Px_S_C1001001_b;
  abcd[117] = 2.0E0*I_ERI_Px_Dxz_Px_S_C1001001_b;
  abcd[118] = 2.0E0*I_ERI_Py_Dxz_Px_S_C1001001_b;
  abcd[119] = 2.0E0*I_ERI_Pz_Dxz_Px_S_C1001001_b;
  abcd[123] = 2.0E0*I_ERI_Px_D2x_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[124] = 2.0E0*I_ERI_Py_D2x_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[125] = 2.0E0*I_ERI_Pz_D2x_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[126] = 2.0E0*I_ERI_Px_Dxy_Py_S_C1001001_b;
  abcd[127] = 2.0E0*I_ERI_Py_Dxy_Py_S_C1001001_b;
  abcd[128] = 2.0E0*I_ERI_Pz_Dxy_Py_S_C1001001_b;
  abcd[129] = 2.0E0*I_ERI_Px_Dxz_Py_S_C1001001_b;
  abcd[130] = 2.0E0*I_ERI_Py_Dxz_Py_S_C1001001_b;
  abcd[131] = 2.0E0*I_ERI_Pz_Dxz_Py_S_C1001001_b;
  abcd[135] = 2.0E0*I_ERI_Px_D2x_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[136] = 2.0E0*I_ERI_Py_D2x_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[137] = 2.0E0*I_ERI_Pz_D2x_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;
  abcd[138] = 2.0E0*I_ERI_Px_Dxy_Pz_S_C1001001_b;
  abcd[139] = 2.0E0*I_ERI_Py_Dxy_Pz_S_C1001001_b;
  abcd[140] = 2.0E0*I_ERI_Pz_Dxy_Pz_S_C1001001_b;
  abcd[141] = 2.0E0*I_ERI_Px_Dxz_Pz_S_C1001001_b;
  abcd[142] = 2.0E0*I_ERI_Py_Dxz_Pz_S_C1001001_b;
  abcd[143] = 2.0E0*I_ERI_Pz_Dxz_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[144] = 2.0E0*I_ERI_Px_Py_Px_S_C1000001_b;
  abcd[145] = 2.0E0*I_ERI_Py_Py_Px_S_C1000001_b;
  abcd[146] = 2.0E0*I_ERI_Pz_Py_Px_S_C1000001_b;
  abcd[156] = 2.0E0*I_ERI_Px_Py_Py_S_C1000001_b;
  abcd[157] = 2.0E0*I_ERI_Py_Py_Py_S_C1000001_b;
  abcd[158] = 2.0E0*I_ERI_Pz_Py_Py_S_C1000001_b;
  abcd[168] = 2.0E0*I_ERI_Px_Py_Pz_S_C1000001_b;
  abcd[169] = 2.0E0*I_ERI_Py_Py_Pz_S_C1000001_b;
  abcd[170] = 2.0E0*I_ERI_Pz_Py_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[147] = 2.0E0*I_ERI_Px_Dxy_Px_S_C1001001_b;
  abcd[148] = 2.0E0*I_ERI_Py_Dxy_Px_S_C1001001_b;
  abcd[149] = 2.0E0*I_ERI_Pz_Dxy_Px_S_C1001001_b;
  abcd[150] = 2.0E0*I_ERI_Px_D2y_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[151] = 2.0E0*I_ERI_Py_D2y_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[152] = 2.0E0*I_ERI_Pz_D2y_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[153] = 2.0E0*I_ERI_Px_Dyz_Px_S_C1001001_b;
  abcd[154] = 2.0E0*I_ERI_Py_Dyz_Px_S_C1001001_b;
  abcd[155] = 2.0E0*I_ERI_Pz_Dyz_Px_S_C1001001_b;
  abcd[159] = 2.0E0*I_ERI_Px_Dxy_Py_S_C1001001_b;
  abcd[160] = 2.0E0*I_ERI_Py_Dxy_Py_S_C1001001_b;
  abcd[161] = 2.0E0*I_ERI_Pz_Dxy_Py_S_C1001001_b;
  abcd[162] = 2.0E0*I_ERI_Px_D2y_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[163] = 2.0E0*I_ERI_Py_D2y_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[164] = 2.0E0*I_ERI_Pz_D2y_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[165] = 2.0E0*I_ERI_Px_Dyz_Py_S_C1001001_b;
  abcd[166] = 2.0E0*I_ERI_Py_Dyz_Py_S_C1001001_b;
  abcd[167] = 2.0E0*I_ERI_Pz_Dyz_Py_S_C1001001_b;
  abcd[171] = 2.0E0*I_ERI_Px_Dxy_Pz_S_C1001001_b;
  abcd[172] = 2.0E0*I_ERI_Py_Dxy_Pz_S_C1001001_b;
  abcd[173] = 2.0E0*I_ERI_Pz_Dxy_Pz_S_C1001001_b;
  abcd[174] = 2.0E0*I_ERI_Px_D2y_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[175] = 2.0E0*I_ERI_Py_D2y_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[176] = 2.0E0*I_ERI_Pz_D2y_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;
  abcd[177] = 2.0E0*I_ERI_Px_Dyz_Pz_S_C1001001_b;
  abcd[178] = 2.0E0*I_ERI_Py_Dyz_Pz_S_C1001001_b;
  abcd[179] = 2.0E0*I_ERI_Pz_Dyz_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[180] = 2.0E0*I_ERI_Px_Pz_Px_S_C1000001_b;
  abcd[181] = 2.0E0*I_ERI_Py_Pz_Px_S_C1000001_b;
  abcd[182] = 2.0E0*I_ERI_Pz_Pz_Px_S_C1000001_b;
  abcd[192] = 2.0E0*I_ERI_Px_Pz_Py_S_C1000001_b;
  abcd[193] = 2.0E0*I_ERI_Py_Pz_Py_S_C1000001_b;
  abcd[194] = 2.0E0*I_ERI_Pz_Pz_Py_S_C1000001_b;
  abcd[204] = 2.0E0*I_ERI_Px_Pz_Pz_S_C1000001_b;
  abcd[205] = 2.0E0*I_ERI_Py_Pz_Pz_S_C1000001_b;
  abcd[206] = 2.0E0*I_ERI_Pz_Pz_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[183] = 2.0E0*I_ERI_Px_Dxz_Px_S_C1001001_b;
  abcd[184] = 2.0E0*I_ERI_Py_Dxz_Px_S_C1001001_b;
  abcd[185] = 2.0E0*I_ERI_Pz_Dxz_Px_S_C1001001_b;
  abcd[186] = 2.0E0*I_ERI_Px_Dyz_Px_S_C1001001_b;
  abcd[187] = 2.0E0*I_ERI_Py_Dyz_Px_S_C1001001_b;
  abcd[188] = 2.0E0*I_ERI_Pz_Dyz_Px_S_C1001001_b;
  abcd[189] = 2.0E0*I_ERI_Px_D2z_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[190] = 2.0E0*I_ERI_Py_D2z_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[191] = 2.0E0*I_ERI_Pz_D2z_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[195] = 2.0E0*I_ERI_Px_Dxz_Py_S_C1001001_b;
  abcd[196] = 2.0E0*I_ERI_Py_Dxz_Py_S_C1001001_b;
  abcd[197] = 2.0E0*I_ERI_Pz_Dxz_Py_S_C1001001_b;
  abcd[198] = 2.0E0*I_ERI_Px_Dyz_Py_S_C1001001_b;
  abcd[199] = 2.0E0*I_ERI_Py_Dyz_Py_S_C1001001_b;
  abcd[200] = 2.0E0*I_ERI_Pz_Dyz_Py_S_C1001001_b;
  abcd[201] = 2.0E0*I_ERI_Px_D2z_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[202] = 2.0E0*I_ERI_Py_D2z_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[203] = 2.0E0*I_ERI_Pz_D2z_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[207] = 2.0E0*I_ERI_Px_Dxz_Pz_S_C1001001_b;
  abcd[208] = 2.0E0*I_ERI_Py_Dxz_Pz_S_C1001001_b;
  abcd[209] = 2.0E0*I_ERI_Pz_Dxz_Pz_S_C1001001_b;
  abcd[210] = 2.0E0*I_ERI_Px_Dyz_Pz_S_C1001001_b;
  abcd[211] = 2.0E0*I_ERI_Py_Dyz_Pz_S_C1001001_b;
  abcd[212] = 2.0E0*I_ERI_Pz_Dyz_Pz_S_C1001001_b;
  abcd[213] = 2.0E0*I_ERI_Px_D2z_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[214] = 2.0E0*I_ERI_Py_D2z_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[215] = 2.0E0*I_ERI_Pz_D2z_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[216] = 2.0E0*I_ERI_Px_S_D2x_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[217] = 2.0E0*I_ERI_Py_S_D2x_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[218] = 2.0E0*I_ERI_Pz_S_D2x_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[228] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[229] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[230] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[240] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[241] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[242] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[219] = 2.0E0*I_ERI_Px_Px_D2x_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[220] = 2.0E0*I_ERI_Py_Px_D2x_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[221] = 2.0E0*I_ERI_Pz_Px_D2x_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[222] = 2.0E0*I_ERI_Px_Py_D2x_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[223] = 2.0E0*I_ERI_Py_Py_D2x_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[224] = 2.0E0*I_ERI_Pz_Py_D2x_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[225] = 2.0E0*I_ERI_Px_Pz_D2x_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[226] = 2.0E0*I_ERI_Py_Pz_D2x_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[227] = 2.0E0*I_ERI_Pz_Pz_D2x_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[231] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[232] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[233] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[234] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[235] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[236] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[237] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[238] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[239] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[243] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[244] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[245] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[246] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[247] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[248] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[249] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[250] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[251] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[252] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[253] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[254] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[264] = 2.0E0*I_ERI_Px_S_D2y_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[265] = 2.0E0*I_ERI_Py_S_D2y_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[266] = 2.0E0*I_ERI_Pz_S_D2y_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[276] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[277] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[278] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[255] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[256] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[257] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[258] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[259] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[260] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[261] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[262] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[263] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[267] = 2.0E0*I_ERI_Px_Px_D2y_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[268] = 2.0E0*I_ERI_Py_Px_D2y_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[269] = 2.0E0*I_ERI_Pz_Px_D2y_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[270] = 2.0E0*I_ERI_Px_Py_D2y_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[271] = 2.0E0*I_ERI_Py_Py_D2y_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[272] = 2.0E0*I_ERI_Pz_Py_D2y_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[273] = 2.0E0*I_ERI_Px_Pz_D2y_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[274] = 2.0E0*I_ERI_Py_Pz_D2y_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[275] = 2.0E0*I_ERI_Pz_Pz_D2y_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[279] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[280] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[281] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[282] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[283] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[284] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[285] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[286] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[287] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[288] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[289] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[290] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;
  abcd[300] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[301] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[302] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;
  abcd[312] = 2.0E0*I_ERI_Px_S_D2z_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[313] = 2.0E0*I_ERI_Py_S_D2z_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[314] = 2.0E0*I_ERI_Pz_S_D2z_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[291] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[292] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[293] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[294] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[295] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[296] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[297] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[298] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[299] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;
  abcd[303] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[304] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[305] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[306] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[307] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[308] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[309] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[310] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[311] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;
  abcd[315] = 2.0E0*I_ERI_Px_Px_D2z_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[316] = 2.0E0*I_ERI_Py_Px_D2z_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[317] = 2.0E0*I_ERI_Pz_Px_D2z_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[318] = 2.0E0*I_ERI_Px_Py_D2z_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[319] = 2.0E0*I_ERI_Py_Py_D2z_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[320] = 2.0E0*I_ERI_Pz_Py_D2z_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[321] = 2.0E0*I_ERI_Px_Pz_D2z_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[322] = 2.0E0*I_ERI_Py_Pz_D2z_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[323] = 2.0E0*I_ERI_Pz_Pz_D2z_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
}
