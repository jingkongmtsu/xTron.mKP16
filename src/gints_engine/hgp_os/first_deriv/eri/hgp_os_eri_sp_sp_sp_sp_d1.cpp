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
// KET1
// X
// Y
// Z
// ####

void hgp_os_eri_sp_sp_sp_sp_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_Px_S_S_S_C0_a = 0.0E0;
  Double I_ERI_Py_S_S_S_C0_a = 0.0E0;
  Double I_ERI_Pz_S_S_S_C0_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1_a = 0.0E0;
  Double I_ERI_S_S_S_S_C1 = 0.0E0;
  Double I_ERI_S_Px_S_S_C1001 = 0.0E0;
  Double I_ERI_S_Py_S_S_C1001 = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1001 = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000000_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000000_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000000_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000000_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000000_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000000_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000000_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000000_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000000_a = 0.0E0;
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
  Double I_ERI_Px_S_S_Px_C1000000000_a = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000000000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000000000_a = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000000000_a = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000000000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000000000_a = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000000000_a = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000000000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000000000_a = 0.0E0;
  Double I_ERI_D2x_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000000001_a = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000000001_a = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000000001_a = 0.0E0;
  Double I_ERI_S_S_S_Px_C1000000001 = 0.0E0;
  Double I_ERI_S_S_S_Py_C1000000001 = 0.0E0;
  Double I_ERI_S_S_S_Pz_C1000000001 = 0.0E0;
  Double I_ERI_S_Px_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_S_Py_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_S_Pz_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_S_Px_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_S_Py_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_S_Pz_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_S_Px_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_S_Py_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_S_Pz_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_S_Px_S_S_C0_b = 0.0E0;
  Double I_ERI_S_Py_S_S_C0_b = 0.0E0;
  Double I_ERI_S_Pz_S_S_C0_b = 0.0E0;
  Double I_ERI_S_D2x_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Dxy_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Dxz_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_D2y_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_Dyz_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_D2z_S_S_C1000_b = 0.0E0;
  Double I_ERI_S_S_S_S_C1000 = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001 = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001 = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001 = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1000000_b = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1000000_b = 0.0E0;
  Double I_ERI_S_D2x_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2y_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2z_Px_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2x_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2y_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2z_Py_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2x_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2y_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_D2z_Pz_S_C1001000_b = 0.0E0;
  Double I_ERI_S_S_Px_S_C1001000 = 0.0E0;
  Double I_ERI_S_S_Py_S_C1001000 = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1001000 = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001 = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001 = 0.0E0;
  Double I_ERI_S_Px_S_Px_C1000000000_b = 0.0E0;
  Double I_ERI_S_Py_S_Px_C1000000000_b = 0.0E0;
  Double I_ERI_S_Pz_S_Px_C1000000000_b = 0.0E0;
  Double I_ERI_S_Px_S_Py_C1000000000_b = 0.0E0;
  Double I_ERI_S_Py_S_Py_C1000000000_b = 0.0E0;
  Double I_ERI_S_Pz_S_Py_C1000000000_b = 0.0E0;
  Double I_ERI_S_Px_S_Pz_C1000000000_b = 0.0E0;
  Double I_ERI_S_Py_S_Pz_C1000000000_b = 0.0E0;
  Double I_ERI_S_Pz_S_Pz_C1000000000_b = 0.0E0;
  Double I_ERI_S_D2x_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxy_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxz_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2y_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dyz_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2z_S_Px_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2x_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxy_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxz_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2y_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dyz_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2z_S_Py_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2x_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxy_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dxz_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2y_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_Dyz_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_D2z_S_Pz_C1000001000_b = 0.0E0;
  Double I_ERI_S_S_S_Px_C1000001000 = 0.0E0;
  Double I_ERI_S_S_S_Py_C1000001000 = 0.0E0;
  Double I_ERI_S_S_S_Pz_C1000001000 = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000001001 = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000001001 = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000001001 = 0.0E0;
  Double I_ERI_S_S_Px_S_C0_c = 0.0E0;
  Double I_ERI_S_S_Py_S_C0_c = 0.0E0;
  Double I_ERI_S_S_Pz_S_C0_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1_c = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1000_c = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1000_c = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1000_c = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1000_c = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1000_c = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1000_c = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1000_c = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1000_c = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1000_c = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1000000_c = 0.0E0;
  Double I_ERI_S_S_S_S_C1000000 = 0.0E0;
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
  Double I_ERI_S_Px_D2x_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001000_c = 0.0E0;
  Double I_ERI_S_Px_S_S_C1001000 = 0.0E0;
  Double I_ERI_S_Py_S_S_C1001000 = 0.0E0;
  Double I_ERI_S_Pz_S_S_C1001000 = 0.0E0;
  Double I_ERI_S_S_S_Px_C1001000000 = 0.0E0;
  Double I_ERI_S_S_S_Py_C1001000000 = 0.0E0;
  Double I_ERI_S_S_S_Pz_C1001000000 = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1001000001 = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1001000001 = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1001000001 = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1001000001 = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1001000001 = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1001000001 = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1001000001 = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1001000001 = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1001000001 = 0.0E0;
  Double I_ERI_S_Px_S_Px_C1001001000 = 0.0E0;
  Double I_ERI_S_Py_S_Px_C1001001000 = 0.0E0;
  Double I_ERI_S_Pz_S_Px_C1001001000 = 0.0E0;
  Double I_ERI_S_Px_S_Py_C1001001000 = 0.0E0;
  Double I_ERI_S_Py_S_Py_C1001001000 = 0.0E0;
  Double I_ERI_S_Pz_S_Py_C1001001000 = 0.0E0;
  Double I_ERI_S_Px_S_Pz_C1001001000 = 0.0E0;
  Double I_ERI_S_Py_S_Pz_C1001001000 = 0.0E0;
  Double I_ERI_S_Pz_S_Pz_C1001001000 = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Px_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Py_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1000_a = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001000_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001000_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001000_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001000_a = 0.0E0;
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
  Double I_ERI_D2x_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000001000_a = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000001000_a = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000001000_a = 0.0E0;
  Double I_ERI_F3x_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F2xy_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F2xz_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F3y_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F2yz_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F3z_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_F3x_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F2xy_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F2xz_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F3y_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F2yz_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F3z_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_F3x_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_F2xy_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_F2xz_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_F3y_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_F2yz_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_F3z_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_D2x_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000001001_a = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000001001_a = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000001001_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Px_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C1_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1_b = 0.0E0;
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
  Double I_ERI_D2x_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000000001_b = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000000001_b = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000000001_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001_c = 0.0E0;
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
  Double I_ERI_D2x_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1001001001 = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Px_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Py_S_S_S_C1001_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_C1001_b = 0.0E0;
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
  Double I_ERI_F3x_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F3y_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F3z_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_F3x_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F3y_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F3z_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_F3x_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_F3y_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_F3z_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_D2x_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Px_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Py_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Px_C1000001001_b = 0.0E0;
  Double I_ERI_Px_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Py_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Py_C1000001001_b = 0.0E0;
  Double I_ERI_Px_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Py_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Pz_S_S_Pz_C1000001001_b = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001000000_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001000000_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001000000_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001000001_a = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Px_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Py_S_C1001000001 = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1001000001 = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001000000_b = 0.0E0;
  Double I_ERI_S_D2x_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_D2x_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Dxy_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Dxz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_D2y_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Dyz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_D2z_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Px_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Py_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2x_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxy_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dxz_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2y_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_Dyz_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_D2z_Pz_S_C1001001000_b = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Px_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Py_S_C1001001000 = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1001001000 = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001001 = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001001 = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Px_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Py_S_C1000000000_c = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1000000000_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000000001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000000001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000000001_c = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1000001000_c = 0.0E0;
  Double I_ERI_S_S_F3x_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_F2xy_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_F2xz_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Fx2y_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Fxyz_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Fx2z_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_F3y_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_F2yz_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Fy2z_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_F3z_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_D2x_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1001000000_c = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1001000000_c = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001000001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001000001_c = 0.0E0;
  Double I_ERI_S_Px_F3x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F3x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F3x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_F2xy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F2xy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F2xy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_F2xz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F2xz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F2xz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Fx2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Fx2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Fx2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Fxyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Fxyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Fxyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Fx2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Fx2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Fx2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_F3y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F3y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F3y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_F2yz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F2yz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F2yz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Fy2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Fy2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Fy2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_F3z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_F3z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_F3z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001001000_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001000_a = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001000_a = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001001_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001000001_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001000001_b = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001001_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001001_c = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001001_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001001_c = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001001_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001001_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001001_b = 0.0E0;

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
      Double jc2_1 = jcoe[jp2+1*jnp2];
      Double jc2_2 = jcoe[jp2+2*jnp2];
      Double jc2_3 = jcoe[jp2+3*jnp2];
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
       * shell quartet name: SQ_ERI_S_S_P_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_Px_S_M4_vrr = QCX*I_ERI_S_S_S_S_M4_vrr+WQX*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Py_S_M4_vrr = QCY*I_ERI_S_S_S_S_M4_vrr+WQY*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Pz_S_M4_vrr = QCZ*I_ERI_S_S_S_S_M4_vrr+WQZ*I_ERI_S_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M3_vrr = QCX*I_ERI_S_S_Px_S_M3_vrr+WQX*I_ERI_S_S_Px_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dxy_S_M3_vrr = QCY*I_ERI_S_S_Px_S_M3_vrr+WQY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_Dxz_S_M3_vrr = QCZ*I_ERI_S_S_Px_S_M3_vrr+WQZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_D2y_S_M3_vrr = QCY*I_ERI_S_S_Py_S_M3_vrr+WQY*I_ERI_S_S_Py_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dyz_S_M3_vrr = QCZ*I_ERI_S_S_Py_S_M3_vrr+WQZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_D2z_S_M3_vrr = QCZ*I_ERI_S_S_Pz_S_M3_vrr+WQZ*I_ERI_S_S_Pz_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M2_vrr = QCX*I_ERI_S_S_Px_S_M2_vrr+WQX*I_ERI_S_S_Px_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dxy_S_M2_vrr = QCY*I_ERI_S_S_Px_S_M2_vrr+WQY*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_Dxz_S_M2_vrr = QCZ*I_ERI_S_S_Px_S_M2_vrr+WQZ*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_D2y_S_M2_vrr = QCY*I_ERI_S_S_Py_S_M2_vrr+WQY*I_ERI_S_S_Py_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dyz_S_M2_vrr = QCZ*I_ERI_S_S_Py_S_M2_vrr+WQZ*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_D2z_S_M2_vrr = QCZ*I_ERI_S_S_Pz_S_M2_vrr+WQZ*I_ERI_S_S_Pz_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M2_vrr = PAX*I_ERI_Px_S_S_S_M2_vrr+WPX*I_ERI_Px_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_S_S_M2_vrr = PAY*I_ERI_Px_S_S_S_M2_vrr+WPY*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_S_S_M2_vrr = PAY*I_ERI_Py_S_S_S_M2_vrr+WPY*I_ERI_Py_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_S_S_M2_vrr = PAZ*I_ERI_Pz_S_S_S_M2_vrr+WPZ*I_ERI_Pz_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M2_vrr = PAX*I_ERI_S_S_D2x_S_M2_vrr+WPX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Py_S_D2x_S_M2_vrr = PAY*I_ERI_S_S_D2x_S_M2_vrr+WPY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_D2x_S_M2_vrr = PAZ*I_ERI_S_S_D2x_S_M2_vrr+WPZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Dxy_S_M2_vrr = PAX*I_ERI_S_S_Dxy_S_M2_vrr+WPX*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Py_S_Dxy_S_M2_vrr = PAY*I_ERI_S_S_Dxy_S_M2_vrr+WPY*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Pz_S_Dxy_S_M2_vrr = PAZ*I_ERI_S_S_Dxy_S_M2_vrr+WPZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Dxz_S_M2_vrr = PAX*I_ERI_S_S_Dxz_S_M2_vrr+WPX*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Py_S_Dxz_S_M2_vrr = PAY*I_ERI_S_S_Dxz_S_M2_vrr+WPY*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Dxz_S_M2_vrr = PAZ*I_ERI_S_S_Dxz_S_M2_vrr+WPZ*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Px_S_D2y_S_M2_vrr = PAX*I_ERI_S_S_D2y_S_M2_vrr+WPX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_D2y_S_M2_vrr = PAY*I_ERI_S_S_D2y_S_M2_vrr+WPY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Pz_S_D2y_S_M2_vrr = PAZ*I_ERI_S_S_D2y_S_M2_vrr+WPZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Dyz_S_M2_vrr = PAX*I_ERI_S_S_Dyz_S_M2_vrr+WPX*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Dyz_S_M2_vrr = PAY*I_ERI_S_S_Dyz_S_M2_vrr+WPY*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Pz_S_Dyz_S_M2_vrr = PAZ*I_ERI_S_S_Dyz_S_M2_vrr+WPZ*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Px_S_D2z_S_M2_vrr = PAX*I_ERI_S_S_D2z_S_M2_vrr+WPX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_D2z_S_M2_vrr = PAY*I_ERI_S_S_D2z_S_M2_vrr+WPY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_D2z_S_M2_vrr = PAZ*I_ERI_S_S_D2z_S_M2_vrr+WPZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Pz_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_S_D_S_S_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_D2x_S_S_M1_vrr = PBX*I_ERI_S_Px_S_S_M1_vrr+WPX*I_ERI_S_Px_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Dxy_S_S_M1_vrr = PBY*I_ERI_S_Px_S_S_M1_vrr+WPY*I_ERI_S_Px_S_S_M2_vrr;
      Double I_ERI_S_Dxz_S_S_M1_vrr = PBZ*I_ERI_S_Px_S_S_M1_vrr+WPZ*I_ERI_S_Px_S_S_M2_vrr;
      Double I_ERI_S_D2y_S_S_M1_vrr = PBY*I_ERI_S_Py_S_S_M1_vrr+WPY*I_ERI_S_Py_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_Dyz_S_S_M1_vrr = PBZ*I_ERI_S_Py_S_S_M1_vrr+WPZ*I_ERI_S_Py_S_S_M2_vrr;
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
       * shell quartet name: SQ_ERI_S_P_D_S_M1
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_S_Px_D2x_S_M1_vrr = PBX*I_ERI_S_S_D2x_S_M1_vrr+WPX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_Py_D2x_S_M1_vrr = PBY*I_ERI_S_S_D2x_S_M1_vrr+WPY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_Pz_D2x_S_M1_vrr = PBZ*I_ERI_S_S_D2x_S_M1_vrr+WPZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_Px_Dxy_S_M1_vrr = PBX*I_ERI_S_S_Dxy_S_M1_vrr+WPX*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_Py_Dxy_S_M1_vrr = PBY*I_ERI_S_S_Dxy_S_M1_vrr+WPY*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_Pz_Dxy_S_M1_vrr = PBZ*I_ERI_S_S_Dxy_S_M1_vrr+WPZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_S_Px_Dxz_S_M1_vrr = PBX*I_ERI_S_S_Dxz_S_M1_vrr+WPX*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_S_Py_Dxz_S_M1_vrr = PBY*I_ERI_S_S_Dxz_S_M1_vrr+WPY*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_S_Pz_Dxz_S_M1_vrr = PBZ*I_ERI_S_S_Dxz_S_M1_vrr+WPZ*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_Px_D2y_S_M1_vrr = PBX*I_ERI_S_S_D2y_S_M1_vrr+WPX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_Py_D2y_S_M1_vrr = PBY*I_ERI_S_S_D2y_S_M1_vrr+WPY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_Pz_D2y_S_M1_vrr = PBZ*I_ERI_S_S_D2y_S_M1_vrr+WPZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_Px_Dyz_S_M1_vrr = PBX*I_ERI_S_S_Dyz_S_M1_vrr+WPX*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_S_Py_Dyz_S_M1_vrr = PBY*I_ERI_S_S_Dyz_S_M1_vrr+WPY*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_S_Pz_Dyz_S_M1_vrr = PBZ*I_ERI_S_S_Dyz_S_M1_vrr+WPZ*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_Px_D2z_S_M1_vrr = PBX*I_ERI_S_S_D2z_S_M1_vrr+WPX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_Py_D2z_S_M1_vrr = PBY*I_ERI_S_S_D2z_S_M1_vrr+WPY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_Pz_D2z_S_M1_vrr = PBZ*I_ERI_S_S_D2z_S_M1_vrr+WPZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M1_vrr = PAX*I_ERI_S_S_D2x_S_M1_vrr+WPX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_D2x_S_M1_vrr = PAY*I_ERI_S_S_D2x_S_M1_vrr+WPY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Pz_S_D2x_S_M1_vrr = PAZ*I_ERI_S_S_D2x_S_M1_vrr+WPZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Px_S_Dxy_S_M1_vrr = PAX*I_ERI_S_S_Dxy_S_M1_vrr+WPX*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Dxy_S_M1_vrr = PAY*I_ERI_S_S_Dxy_S_M1_vrr+WPY*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxy_S_M1_vrr = PAZ*I_ERI_S_S_Dxy_S_M1_vrr+WPZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Px_S_Dxz_S_M1_vrr = PAX*I_ERI_S_S_Dxz_S_M1_vrr+WPX*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Py_S_Dxz_S_M1_vrr = PAY*I_ERI_S_S_Dxz_S_M1_vrr+WPY*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Pz_S_Dxz_S_M1_vrr = PAZ*I_ERI_S_S_Dxz_S_M1_vrr+WPZ*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_D2y_S_M1_vrr = PAX*I_ERI_S_S_D2y_S_M1_vrr+WPX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Py_S_D2y_S_M1_vrr = PAY*I_ERI_S_S_D2y_S_M1_vrr+WPY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_D2y_S_M1_vrr = PAZ*I_ERI_S_S_D2y_S_M1_vrr+WPZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Px_S_Dyz_S_M1_vrr = PAX*I_ERI_S_S_Dyz_S_M1_vrr+WPX*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Py_S_Dyz_S_M1_vrr = PAY*I_ERI_S_S_Dyz_S_M1_vrr+WPY*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Pz_S_Dyz_S_M1_vrr = PAZ*I_ERI_S_S_Dyz_S_M1_vrr+WPZ*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_D2z_S_M1_vrr = PAX*I_ERI_S_S_D2z_S_M1_vrr+WPX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Py_S_D2z_S_M1_vrr = PAY*I_ERI_S_S_D2z_S_M1_vrr+WPY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Pz_S_D2z_S_M1_vrr = PAZ*I_ERI_S_S_D2z_S_M1_vrr+WPZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Pz_S_M2_vrr;

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
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M1_vrr = PAX*I_ERI_Px_S_D2x_S_M1_vrr+WPX*I_ERI_Px_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxy_S_D2x_S_M1_vrr = PAY*I_ERI_Px_S_D2x_S_M1_vrr+WPY*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxz_S_D2x_S_M1_vrr = PAZ*I_ERI_Px_S_D2x_S_M1_vrr+WPZ*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_D2x_S_M1_vrr = PAY*I_ERI_Py_S_D2x_S_M1_vrr+WPY*I_ERI_Py_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Dyz_S_D2x_S_M1_vrr = PAZ*I_ERI_Py_S_D2x_S_M1_vrr+WPZ*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_D2x_S_M1_vrr = PAZ*I_ERI_Pz_S_D2x_S_M1_vrr+WPZ*I_ERI_Pz_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Dxy_S_M1_vrr = PAX*I_ERI_Px_S_Dxy_S_M1_vrr+WPX*I_ERI_Px_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxy_S_M1_vrr = PAY*I_ERI_Px_S_Dxy_S_M1_vrr+WPY*I_ERI_Px_S_Dxy_S_M2_vrr+oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxz_S_Dxy_S_M1_vrr = PAZ*I_ERI_Px_S_Dxy_S_M1_vrr+WPZ*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_D2y_S_Dxy_S_M1_vrr = PAY*I_ERI_Py_S_Dxy_S_M1_vrr+WPY*I_ERI_Py_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_Dyz_S_Dxy_S_M1_vrr = PAZ*I_ERI_Py_S_Dxy_S_M1_vrr+WPZ*I_ERI_Py_S_Dxy_S_M2_vrr;
      Double I_ERI_D2z_S_Dxy_S_M1_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M1_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Dxz_S_M1_vrr = PAX*I_ERI_Px_S_Dxz_S_M1_vrr+WPX*I_ERI_Px_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxz_S_M1_vrr = PAY*I_ERI_Px_S_Dxz_S_M1_vrr+WPY*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_Dxz_S_M1_vrr = PAY*I_ERI_Py_S_Dxz_S_M1_vrr+WPY*I_ERI_Py_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Dxz_S_M1_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M1_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_D2y_S_M1_vrr = PAX*I_ERI_Px_S_D2y_S_M1_vrr+WPX*I_ERI_Px_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = PAY*I_ERI_Px_S_D2y_S_M1_vrr+WPY*I_ERI_Px_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = PAZ*I_ERI_Px_S_D2y_S_M1_vrr+WPZ*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = PAY*I_ERI_Py_S_D2y_S_M1_vrr+WPY*I_ERI_Py_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = PAZ*I_ERI_Py_S_D2y_S_M1_vrr+WPZ*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = PAZ*I_ERI_Pz_S_D2y_S_M1_vrr+WPZ*I_ERI_Pz_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Dyz_S_M1_vrr = PAX*I_ERI_Px_S_Dyz_S_M1_vrr+WPX*I_ERI_Px_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Dxy_S_Dyz_S_M1_vrr = PAY*I_ERI_Px_S_Dyz_S_M1_vrr+WPY*I_ERI_Px_S_Dyz_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_Dyz_S_M1_vrr = PAY*I_ERI_Py_S_Dyz_S_M1_vrr+WPY*I_ERI_Py_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_Dyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_D2z_S_M1_vrr = PAX*I_ERI_Px_S_D2z_S_M1_vrr+WPX*I_ERI_Px_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Dxy_S_D2z_S_M1_vrr = PAY*I_ERI_Px_S_D2z_S_M1_vrr+WPY*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_Dxz_S_D2z_S_M1_vrr = PAZ*I_ERI_Px_S_D2z_S_M1_vrr+WPZ*I_ERI_Px_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_D2z_S_M1_vrr = PAY*I_ERI_Py_S_D2z_S_M1_vrr+WPY*I_ERI_Py_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Dyz_S_D2z_S_M1_vrr = PAZ*I_ERI_Py_S_D2z_S_M1_vrr+WPZ*I_ERI_Py_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_D2z_S_M1_vrr = PAZ*I_ERI_Pz_S_D2z_S_M1_vrr+WPZ*I_ERI_Pz_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_P_S_P
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_P
       * RHS shell quartet name: SQ_ERI_S_S_S_P_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_Px_S_Px_vrr = PBX*I_ERI_S_S_S_Px_vrr+WPX*I_ERI_S_S_S_Px_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Py_S_Px_vrr = PBY*I_ERI_S_S_S_Px_vrr+WPY*I_ERI_S_S_S_Px_M1_vrr;
      Double I_ERI_S_Pz_S_Px_vrr = PBZ*I_ERI_S_S_S_Px_vrr+WPZ*I_ERI_S_S_S_Px_M1_vrr;
      Double I_ERI_S_Px_S_Py_vrr = PBX*I_ERI_S_S_S_Py_vrr+WPX*I_ERI_S_S_S_Py_M1_vrr;
      Double I_ERI_S_Py_S_Py_vrr = PBY*I_ERI_S_S_S_Py_vrr+WPY*I_ERI_S_S_S_Py_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_Pz_S_Py_vrr = PBZ*I_ERI_S_S_S_Py_vrr+WPZ*I_ERI_S_S_S_Py_M1_vrr;
      Double I_ERI_S_Px_S_Pz_vrr = PBX*I_ERI_S_S_S_Pz_vrr+WPX*I_ERI_S_S_S_Pz_M1_vrr;
      Double I_ERI_S_Py_S_Pz_vrr = PBY*I_ERI_S_S_S_Pz_vrr+WPY*I_ERI_S_S_S_Pz_M1_vrr;
      Double I_ERI_S_Pz_S_Pz_vrr = PBZ*I_ERI_S_S_S_Pz_vrr+WPZ*I_ERI_S_S_S_Pz_M1_vrr+oned2k*I_ERI_S_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_D_S_P
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_D_S_S
       * RHS shell quartet name: SQ_ERI_S_D_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_S_S_M1
       ************************************************************/
      Double I_ERI_S_D2x_S_Px_vrr = QDX*I_ERI_S_D2x_S_S_vrr+WQX*I_ERI_S_D2x_S_S_M1_vrr+2*oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Dxy_S_Px_vrr = QDX*I_ERI_S_Dxy_S_S_vrr+WQX*I_ERI_S_Dxy_S_S_M1_vrr+oned2k*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_Dxz_S_Px_vrr = QDX*I_ERI_S_Dxz_S_S_vrr+WQX*I_ERI_S_Dxz_S_S_M1_vrr+oned2k*I_ERI_S_Pz_S_S_M1_vrr;
      Double I_ERI_S_D2y_S_Px_vrr = QDX*I_ERI_S_D2y_S_S_vrr+WQX*I_ERI_S_D2y_S_S_M1_vrr;
      Double I_ERI_S_Dyz_S_Px_vrr = QDX*I_ERI_S_Dyz_S_S_vrr+WQX*I_ERI_S_Dyz_S_S_M1_vrr;
      Double I_ERI_S_D2z_S_Px_vrr = QDX*I_ERI_S_D2z_S_S_vrr+WQX*I_ERI_S_D2z_S_S_M1_vrr;
      Double I_ERI_S_D2x_S_Py_vrr = QDY*I_ERI_S_D2x_S_S_vrr+WQY*I_ERI_S_D2x_S_S_M1_vrr;
      Double I_ERI_S_Dxy_S_Py_vrr = QDY*I_ERI_S_Dxy_S_S_vrr+WQY*I_ERI_S_Dxy_S_S_M1_vrr+oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_Dxz_S_Py_vrr = QDY*I_ERI_S_Dxz_S_S_vrr+WQY*I_ERI_S_Dxz_S_S_M1_vrr;
      Double I_ERI_S_D2y_S_Py_vrr = QDY*I_ERI_S_D2y_S_S_vrr+WQY*I_ERI_S_D2y_S_S_M1_vrr+2*oned2k*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_Dyz_S_Py_vrr = QDY*I_ERI_S_Dyz_S_S_vrr+WQY*I_ERI_S_Dyz_S_S_M1_vrr+oned2k*I_ERI_S_Pz_S_S_M1_vrr;
      Double I_ERI_S_D2z_S_Py_vrr = QDY*I_ERI_S_D2z_S_S_vrr+WQY*I_ERI_S_D2z_S_S_M1_vrr;
      Double I_ERI_S_D2x_S_Pz_vrr = QDZ*I_ERI_S_D2x_S_S_vrr+WQZ*I_ERI_S_D2x_S_S_M1_vrr;
      Double I_ERI_S_Dxy_S_Pz_vrr = QDZ*I_ERI_S_Dxy_S_S_vrr+WQZ*I_ERI_S_Dxy_S_S_M1_vrr;
      Double I_ERI_S_Dxz_S_Pz_vrr = QDZ*I_ERI_S_Dxz_S_S_vrr+WQZ*I_ERI_S_Dxz_S_S_M1_vrr+oned2k*I_ERI_S_Px_S_S_M1_vrr;
      Double I_ERI_S_D2y_S_Pz_vrr = QDZ*I_ERI_S_D2y_S_S_vrr+WQZ*I_ERI_S_D2y_S_S_M1_vrr;
      Double I_ERI_S_Dyz_S_Pz_vrr = QDZ*I_ERI_S_Dyz_S_S_vrr+WQZ*I_ERI_S_Dyz_S_S_M1_vrr+oned2k*I_ERI_S_Py_S_S_M1_vrr;
      Double I_ERI_S_D2z_S_Pz_vrr = QDZ*I_ERI_S_D2z_S_S_vrr+WQZ*I_ERI_S_D2z_S_S_M1_vrr+2*oned2k*I_ERI_S_Pz_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_S_F3x_S_vrr = QCX*I_ERI_S_S_D2x_S_vrr+WQX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2e*I_ERI_S_S_Px_S_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_F2xy_S_vrr = QCY*I_ERI_S_S_D2x_S_vrr+WQY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_F2xz_S_vrr = QCZ*I_ERI_S_S_D2x_S_vrr+WQZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Fx2y_S_vrr = QCX*I_ERI_S_S_D2y_S_vrr+WQX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fxyz_S_vrr = QCZ*I_ERI_S_S_Dxy_S_vrr+WQZ*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_S_S_Fx2z_S_vrr = QCX*I_ERI_S_S_D2z_S_vrr+WQX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3y_S_vrr = QCY*I_ERI_S_S_D2y_S_vrr+WQY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2e*I_ERI_S_S_Py_S_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_F2yz_S_vrr = QCZ*I_ERI_S_S_D2y_S_vrr+WQZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fy2z_S_vrr = QCY*I_ERI_S_S_D2z_S_vrr+WQY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3z_S_vrr = QCZ*I_ERI_S_S_D2z_S_vrr+WQZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2e*I_ERI_S_S_Pz_S_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_D_D_S
       * expanding position: BRA2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_D_S
       * RHS shell quartet name: SQ_ERI_S_P_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_P_S_M1
       ************************************************************/
      Double I_ERI_S_D2x_D2x_S_vrr = PBX*I_ERI_S_Px_D2x_S_vrr+WPX*I_ERI_S_Px_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_Dxy_D2x_S_vrr = PBY*I_ERI_S_Px_D2x_S_vrr+WPY*I_ERI_S_Px_D2x_S_M1_vrr;
      Double I_ERI_S_Dxz_D2x_S_vrr = PBZ*I_ERI_S_Px_D2x_S_vrr+WPZ*I_ERI_S_Px_D2x_S_M1_vrr;
      Double I_ERI_S_D2y_D2x_S_vrr = PBY*I_ERI_S_Py_D2x_S_vrr+WPY*I_ERI_S_Py_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Dyz_D2x_S_vrr = PBZ*I_ERI_S_Py_D2x_S_vrr+WPZ*I_ERI_S_Py_D2x_S_M1_vrr;
      Double I_ERI_S_D2z_D2x_S_vrr = PBZ*I_ERI_S_Pz_D2x_S_vrr+WPZ*I_ERI_S_Pz_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_D2x_Dxy_S_vrr = PBX*I_ERI_S_Px_Dxy_S_vrr+WPX*I_ERI_S_Px_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_Dxy_Dxy_S_vrr = PBY*I_ERI_S_Px_Dxy_S_vrr+WPY*I_ERI_S_Px_Dxy_S_M1_vrr+oned2k*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_Dxz_Dxy_S_vrr = PBZ*I_ERI_S_Px_Dxy_S_vrr+WPZ*I_ERI_S_Px_Dxy_S_M1_vrr;
      Double I_ERI_S_D2y_Dxy_S_vrr = PBY*I_ERI_S_Py_Dxy_S_vrr+WPY*I_ERI_S_Py_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_Py_Px_S_M1_vrr;
      Double I_ERI_S_Dyz_Dxy_S_vrr = PBZ*I_ERI_S_Py_Dxy_S_vrr+WPZ*I_ERI_S_Py_Dxy_S_M1_vrr;
      Double I_ERI_S_D2z_Dxy_S_vrr = PBZ*I_ERI_S_Pz_Dxy_S_vrr+WPZ*I_ERI_S_Pz_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_S_D2x_Dxz_S_vrr = PBX*I_ERI_S_Px_Dxz_S_vrr+WPX*I_ERI_S_Px_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_Px_Pz_S_M1_vrr;
      Double I_ERI_S_Dxy_Dxz_S_vrr = PBY*I_ERI_S_Px_Dxz_S_vrr+WPY*I_ERI_S_Px_Dxz_S_M1_vrr;
      Double I_ERI_S_Dxz_Dxz_S_vrr = PBZ*I_ERI_S_Px_Dxz_S_vrr+WPZ*I_ERI_S_Px_Dxz_S_M1_vrr+oned2k*I_ERI_S_Px_Px_S_M1_vrr;
      Double I_ERI_S_D2y_Dxz_S_vrr = PBY*I_ERI_S_Py_Dxz_S_vrr+WPY*I_ERI_S_Py_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_S_Dyz_Dxz_S_vrr = PBZ*I_ERI_S_Py_Dxz_S_vrr+WPZ*I_ERI_S_Py_Dxz_S_M1_vrr+oned2k*I_ERI_S_Py_Px_S_M1_vrr;
      Double I_ERI_S_D2z_Dxz_S_vrr = PBZ*I_ERI_S_Pz_Dxz_S_vrr+WPZ*I_ERI_S_Pz_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_Pz_Px_S_M1_vrr;
      Double I_ERI_S_D2x_D2y_S_vrr = PBX*I_ERI_S_Px_D2y_S_vrr+WPX*I_ERI_S_Px_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Dxy_D2y_S_vrr = PBY*I_ERI_S_Px_D2y_S_vrr+WPY*I_ERI_S_Px_D2y_S_M1_vrr+2*oned2k*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_Dxz_D2y_S_vrr = PBZ*I_ERI_S_Px_D2y_S_vrr+WPZ*I_ERI_S_Px_D2y_S_M1_vrr;
      Double I_ERI_S_D2y_D2y_S_vrr = PBY*I_ERI_S_Py_D2y_S_vrr+WPY*I_ERI_S_Py_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_S_Py_Py_S_M1_vrr;
      Double I_ERI_S_Dyz_D2y_S_vrr = PBZ*I_ERI_S_Py_D2y_S_vrr+WPZ*I_ERI_S_Py_D2y_S_M1_vrr;
      Double I_ERI_S_D2z_D2y_S_vrr = PBZ*I_ERI_S_Pz_D2y_S_vrr+WPZ*I_ERI_S_Pz_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_D2x_Dyz_S_vrr = PBX*I_ERI_S_Px_Dyz_S_vrr+WPX*I_ERI_S_Px_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_S_Dxy_Dyz_S_vrr = PBY*I_ERI_S_Px_Dyz_S_vrr+WPY*I_ERI_S_Px_Dyz_S_M1_vrr+oned2k*I_ERI_S_Px_Pz_S_M1_vrr;
      Double I_ERI_S_Dxz_Dyz_S_vrr = PBZ*I_ERI_S_Px_Dyz_S_vrr+WPZ*I_ERI_S_Px_Dyz_S_M1_vrr+oned2k*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_D2y_Dyz_S_vrr = PBY*I_ERI_S_Py_Dyz_S_vrr+WPY*I_ERI_S_Py_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_Py_Pz_S_M1_vrr;
      Double I_ERI_S_Dyz_Dyz_S_vrr = PBZ*I_ERI_S_Py_Dyz_S_vrr+WPZ*I_ERI_S_Py_Dyz_S_M1_vrr+oned2k*I_ERI_S_Py_Py_S_M1_vrr;
      Double I_ERI_S_D2z_Dyz_S_vrr = PBZ*I_ERI_S_Pz_Dyz_S_vrr+WPZ*I_ERI_S_Pz_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_Pz_Py_S_M1_vrr;
      Double I_ERI_S_D2x_D2z_S_vrr = PBX*I_ERI_S_Px_D2z_S_vrr+WPX*I_ERI_S_Px_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Dxy_D2z_S_vrr = PBY*I_ERI_S_Px_D2z_S_vrr+WPY*I_ERI_S_Px_D2z_S_M1_vrr;
      Double I_ERI_S_Dxz_D2z_S_vrr = PBZ*I_ERI_S_Px_D2z_S_vrr+WPZ*I_ERI_S_Px_D2z_S_M1_vrr+2*oned2k*I_ERI_S_Px_Pz_S_M1_vrr;
      Double I_ERI_S_D2y_D2z_S_vrr = PBY*I_ERI_S_Py_D2z_S_vrr+WPY*I_ERI_S_Py_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Dyz_D2z_S_vrr = PBZ*I_ERI_S_Py_D2z_S_vrr+WPZ*I_ERI_S_Py_D2z_S_M1_vrr+2*oned2k*I_ERI_S_Py_Pz_S_M1_vrr;
      Double I_ERI_S_D2z_D2z_S_vrr = PBZ*I_ERI_S_Pz_D2z_S_vrr+WPZ*I_ERI_S_Pz_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_S_Pz_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_S_P_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_P_D_S
       * RHS shell quartet name: SQ_ERI_S_P_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_P_P_S
       * RHS shell quartet name: SQ_ERI_S_P_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_S_Px_F3x_S_vrr = QCX*I_ERI_S_Px_D2x_S_vrr+WQX*I_ERI_S_Px_D2x_S_M1_vrr+2*oned2e*I_ERI_S_Px_Px_S_vrr-2*rhod2esq*I_ERI_S_Px_Px_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Py_F3x_S_vrr = QCX*I_ERI_S_Py_D2x_S_vrr+WQX*I_ERI_S_Py_D2x_S_M1_vrr+2*oned2e*I_ERI_S_Py_Px_S_vrr-2*rhod2esq*I_ERI_S_Py_Px_S_M1_vrr;
      Double I_ERI_S_Pz_F3x_S_vrr = QCX*I_ERI_S_Pz_D2x_S_vrr+WQX*I_ERI_S_Pz_D2x_S_M1_vrr+2*oned2e*I_ERI_S_Pz_Px_S_vrr-2*rhod2esq*I_ERI_S_Pz_Px_S_M1_vrr;
      Double I_ERI_S_Px_F2xy_S_vrr = QCY*I_ERI_S_Px_D2x_S_vrr+WQY*I_ERI_S_Px_D2x_S_M1_vrr;
      Double I_ERI_S_Py_F2xy_S_vrr = QCY*I_ERI_S_Py_D2x_S_vrr+WQY*I_ERI_S_Py_D2x_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Pz_F2xy_S_vrr = QCY*I_ERI_S_Pz_D2x_S_vrr+WQY*I_ERI_S_Pz_D2x_S_M1_vrr;
      Double I_ERI_S_Px_F2xz_S_vrr = QCZ*I_ERI_S_Px_D2x_S_vrr+WQZ*I_ERI_S_Px_D2x_S_M1_vrr;
      Double I_ERI_S_Py_F2xz_S_vrr = QCZ*I_ERI_S_Py_D2x_S_vrr+WQZ*I_ERI_S_Py_D2x_S_M1_vrr;
      Double I_ERI_S_Pz_F2xz_S_vrr = QCZ*I_ERI_S_Pz_D2x_S_vrr+WQZ*I_ERI_S_Pz_D2x_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_Px_Fx2y_S_vrr = QCX*I_ERI_S_Px_D2y_S_vrr+WQX*I_ERI_S_Px_D2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Py_Fx2y_S_vrr = QCX*I_ERI_S_Py_D2y_S_vrr+WQX*I_ERI_S_Py_D2y_S_M1_vrr;
      Double I_ERI_S_Pz_Fx2y_S_vrr = QCX*I_ERI_S_Pz_D2y_S_vrr+WQX*I_ERI_S_Pz_D2y_S_M1_vrr;
      Double I_ERI_S_Px_Fxyz_S_vrr = QCZ*I_ERI_S_Px_Dxy_S_vrr+WQZ*I_ERI_S_Px_Dxy_S_M1_vrr;
      Double I_ERI_S_Py_Fxyz_S_vrr = QCZ*I_ERI_S_Py_Dxy_S_vrr+WQZ*I_ERI_S_Py_Dxy_S_M1_vrr;
      Double I_ERI_S_Pz_Fxyz_S_vrr = QCZ*I_ERI_S_Pz_Dxy_S_vrr+WQZ*I_ERI_S_Pz_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_S_Px_Fx2z_S_vrr = QCX*I_ERI_S_Px_D2z_S_vrr+WQX*I_ERI_S_Px_D2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Py_Fx2z_S_vrr = QCX*I_ERI_S_Py_D2z_S_vrr+WQX*I_ERI_S_Py_D2z_S_M1_vrr;
      Double I_ERI_S_Pz_Fx2z_S_vrr = QCX*I_ERI_S_Pz_D2z_S_vrr+WQX*I_ERI_S_Pz_D2z_S_M1_vrr;
      Double I_ERI_S_Px_F3y_S_vrr = QCY*I_ERI_S_Px_D2y_S_vrr+WQY*I_ERI_S_Px_D2y_S_M1_vrr+2*oned2e*I_ERI_S_Px_Py_S_vrr-2*rhod2esq*I_ERI_S_Px_Py_S_M1_vrr;
      Double I_ERI_S_Py_F3y_S_vrr = QCY*I_ERI_S_Py_D2y_S_vrr+WQY*I_ERI_S_Py_D2y_S_M1_vrr+2*oned2e*I_ERI_S_Py_Py_S_vrr-2*rhod2esq*I_ERI_S_Py_Py_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Pz_F3y_S_vrr = QCY*I_ERI_S_Pz_D2y_S_vrr+WQY*I_ERI_S_Pz_D2y_S_M1_vrr+2*oned2e*I_ERI_S_Pz_Py_S_vrr-2*rhod2esq*I_ERI_S_Pz_Py_S_M1_vrr;
      Double I_ERI_S_Px_F2yz_S_vrr = QCZ*I_ERI_S_Px_D2y_S_vrr+WQZ*I_ERI_S_Px_D2y_S_M1_vrr;
      Double I_ERI_S_Py_F2yz_S_vrr = QCZ*I_ERI_S_Py_D2y_S_vrr+WQZ*I_ERI_S_Py_D2y_S_M1_vrr;
      Double I_ERI_S_Pz_F2yz_S_vrr = QCZ*I_ERI_S_Pz_D2y_S_vrr+WQZ*I_ERI_S_Pz_D2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_Px_Fy2z_S_vrr = QCY*I_ERI_S_Px_D2z_S_vrr+WQY*I_ERI_S_Px_D2z_S_M1_vrr;
      Double I_ERI_S_Py_Fy2z_S_vrr = QCY*I_ERI_S_Py_D2z_S_vrr+WQY*I_ERI_S_Py_D2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_Pz_Fy2z_S_vrr = QCY*I_ERI_S_Pz_D2z_S_vrr+WQY*I_ERI_S_Pz_D2z_S_M1_vrr;
      Double I_ERI_S_Px_F3z_S_vrr = QCZ*I_ERI_S_Px_D2z_S_vrr+WQZ*I_ERI_S_Px_D2z_S_M1_vrr+2*oned2e*I_ERI_S_Px_Pz_S_vrr-2*rhod2esq*I_ERI_S_Px_Pz_S_M1_vrr;
      Double I_ERI_S_Py_F3z_S_vrr = QCZ*I_ERI_S_Py_D2z_S_vrr+WQZ*I_ERI_S_Py_D2z_S_M1_vrr+2*oned2e*I_ERI_S_Py_Pz_S_vrr-2*rhod2esq*I_ERI_S_Py_Pz_S_M1_vrr;
      Double I_ERI_S_Pz_F3z_S_vrr = QCZ*I_ERI_S_Pz_D2z_S_vrr+WQZ*I_ERI_S_Pz_D2z_S_M1_vrr+2*oned2e*I_ERI_S_Pz_Pz_S_vrr-2*rhod2esq*I_ERI_S_Pz_Pz_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S_C0_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C0_a_coefs = ic2*jc2*alpha;
      I_ERI_Px_S_S_S_C0_a += SQ_ERI_P_S_S_S_C0_a_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C0_a += SQ_ERI_P_S_S_S_C0_a_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C0_a += SQ_ERI_P_S_S_S_C0_a_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1_a_coefs = ic2_1*jc2*alpha;
      I_ERI_D2x_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1_a += SQ_ERI_D_S_S_S_C1_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C1
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C1_coefs = ic2_1*jc2;
      I_ERI_S_S_S_S_C1 += SQ_ERI_S_S_S_S_C1_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1001_coefs = ic2_3*jc2;
      I_ERI_S_Px_S_S_C1001 += SQ_ERI_S_P_S_S_C1001_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1001 += SQ_ERI_S_P_S_S_C1001_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1001 += SQ_ERI_S_P_S_S_C1001_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000000_a_coefs = ic2*jc2_1*alpha;
      I_ERI_Px_S_Px_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000000_a += SQ_ERI_P_S_P_S_C1000000_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001_a_coefs = ic2_1*jc2_1*alpha;
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
      Double SQ_ERI_S_S_P_S_C1000001_coefs = ic2_1*jc2_1;
      I_ERI_S_S_Px_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001001_coefs = ic2_3*jc2_1;
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
       * shell quartet name: SQ_ERI_P_S_S_P_C1000000000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000000000_a_coefs = ic2*jc2_2*alpha;
      I_ERI_Px_S_S_Px_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000000000_a += SQ_ERI_P_S_S_P_C1000000000_a_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1000000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000000001_a_coefs = ic2_1*jc2_2*alpha;
      I_ERI_D2x_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000000001_a += SQ_ERI_D_S_S_P_C1000000001_a_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_P_C1000000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_P_C1000000001_coefs = ic2_1*jc2_2;
      I_ERI_S_S_S_Px_C1000000001 += SQ_ERI_S_S_S_P_C1000000001_coefs*I_ERI_S_S_S_Px_vrr;
      I_ERI_S_S_S_Py_C1000000001 += SQ_ERI_S_S_S_P_C1000000001_coefs*I_ERI_S_S_S_Py_vrr;
      I_ERI_S_S_S_Pz_C1000000001 += SQ_ERI_S_S_S_P_C1000000001_coefs*I_ERI_S_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_P_C1000001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_P_C1000001001_coefs = ic2_3*jc2_2;
      I_ERI_S_Px_S_Px_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Px_S_Px_vrr;
      I_ERI_S_Py_S_Px_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Py_S_Px_vrr;
      I_ERI_S_Pz_S_Px_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Pz_S_Px_vrr;
      I_ERI_S_Px_S_Py_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Px_S_Py_vrr;
      I_ERI_S_Py_S_Py_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Py_S_Py_vrr;
      I_ERI_S_Pz_S_Py_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Pz_S_Py_vrr;
      I_ERI_S_Px_S_Pz_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Px_S_Pz_vrr;
      I_ERI_S_Py_S_Pz_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Py_S_Pz_vrr;
      I_ERI_S_Pz_S_Pz_C1000001001 += SQ_ERI_S_P_S_P_C1000001001_coefs*I_ERI_S_Pz_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C0_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C0_b_coefs = ic2*jc2*beta;
      I_ERI_S_Px_S_S_C0_b += SQ_ERI_S_P_S_S_C0_b_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C0_b += SQ_ERI_S_P_S_S_C0_b_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C0_b += SQ_ERI_S_P_S_S_C0_b_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_S_S_C1000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_S_S_C1000_b_coefs = ic2_2*jc2*beta;
      I_ERI_S_D2x_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_D2x_S_S_vrr;
      I_ERI_S_Dxy_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_Dxy_S_S_vrr;
      I_ERI_S_Dxz_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_Dxz_S_S_vrr;
      I_ERI_S_D2y_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_D2y_S_S_vrr;
      I_ERI_S_Dyz_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_Dyz_S_S_vrr;
      I_ERI_S_D2z_S_S_C1000_b += SQ_ERI_S_D_S_S_C1000_b_coefs*I_ERI_S_D2z_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C1000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C1000_coefs = ic2_2*jc2;
      I_ERI_S_S_S_S_C1000 += SQ_ERI_S_S_S_S_C1000_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C1001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C1001_coefs = ic2_3*jc2;
      I_ERI_Px_S_S_S_C1001 += SQ_ERI_P_S_S_S_C1001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001 += SQ_ERI_P_S_S_S_C1001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001 += SQ_ERI_P_S_S_S_C1001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1000000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1000000_b_coefs = ic2*jc2_1*beta;
      I_ERI_S_Px_Px_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1000000_b += SQ_ERI_S_P_P_S_C1000000_b_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_P_S_C1001000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_P_S_C1001000_b_coefs = ic2_2*jc2_1*beta;
      I_ERI_S_D2x_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2x_Px_S_vrr;
      I_ERI_S_Dxy_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxy_Px_S_vrr;
      I_ERI_S_Dxz_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxz_Px_S_vrr;
      I_ERI_S_D2y_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2y_Px_S_vrr;
      I_ERI_S_Dyz_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dyz_Px_S_vrr;
      I_ERI_S_D2z_Px_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2z_Px_S_vrr;
      I_ERI_S_D2x_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2x_Py_S_vrr;
      I_ERI_S_Dxy_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxy_Py_S_vrr;
      I_ERI_S_Dxz_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxz_Py_S_vrr;
      I_ERI_S_D2y_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2y_Py_S_vrr;
      I_ERI_S_Dyz_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dyz_Py_S_vrr;
      I_ERI_S_D2z_Py_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2z_Py_S_vrr;
      I_ERI_S_D2x_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2x_Pz_S_vrr;
      I_ERI_S_Dxy_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxy_Pz_S_vrr;
      I_ERI_S_Dxz_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dxz_Pz_S_vrr;
      I_ERI_S_D2y_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2y_Pz_S_vrr;
      I_ERI_S_Dyz_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_Dyz_Pz_S_vrr;
      I_ERI_S_D2z_Pz_S_C1001000_b += SQ_ERI_S_D_P_S_C1001000_b_coefs*I_ERI_S_D2z_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1001000_coefs = ic2_2*jc2_1;
      I_ERI_S_S_Px_S_C1001000 += SQ_ERI_S_S_P_S_C1001000_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1001000 += SQ_ERI_S_S_P_S_C1001000_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1001000 += SQ_ERI_S_S_P_S_C1001000_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_coefs = ic2_3*jc2_1;
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
       * shell quartet name: SQ_ERI_S_P_S_P_C1000000000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_P_C1000000000_b_coefs = ic2*jc2_2*beta;
      I_ERI_S_Px_S_Px_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Px_S_Px_vrr;
      I_ERI_S_Py_S_Px_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Py_S_Px_vrr;
      I_ERI_S_Pz_S_Px_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Pz_S_Px_vrr;
      I_ERI_S_Px_S_Py_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Px_S_Py_vrr;
      I_ERI_S_Py_S_Py_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Py_S_Py_vrr;
      I_ERI_S_Pz_S_Py_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Pz_S_Py_vrr;
      I_ERI_S_Px_S_Pz_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Px_S_Pz_vrr;
      I_ERI_S_Py_S_Pz_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Py_S_Pz_vrr;
      I_ERI_S_Pz_S_Pz_C1000000000_b += SQ_ERI_S_P_S_P_C1000000000_b_coefs*I_ERI_S_Pz_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_S_P_C1000001000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_S_P_C1000001000_b_coefs = ic2_2*jc2_2*beta;
      I_ERI_S_D2x_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2x_S_Px_vrr;
      I_ERI_S_Dxy_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxy_S_Px_vrr;
      I_ERI_S_Dxz_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxz_S_Px_vrr;
      I_ERI_S_D2y_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2y_S_Px_vrr;
      I_ERI_S_Dyz_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dyz_S_Px_vrr;
      I_ERI_S_D2z_S_Px_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2z_S_Px_vrr;
      I_ERI_S_D2x_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2x_S_Py_vrr;
      I_ERI_S_Dxy_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxy_S_Py_vrr;
      I_ERI_S_Dxz_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxz_S_Py_vrr;
      I_ERI_S_D2y_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2y_S_Py_vrr;
      I_ERI_S_Dyz_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dyz_S_Py_vrr;
      I_ERI_S_D2z_S_Py_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2z_S_Py_vrr;
      I_ERI_S_D2x_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2x_S_Pz_vrr;
      I_ERI_S_Dxy_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxy_S_Pz_vrr;
      I_ERI_S_Dxz_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dxz_S_Pz_vrr;
      I_ERI_S_D2y_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2y_S_Pz_vrr;
      I_ERI_S_Dyz_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_Dyz_S_Pz_vrr;
      I_ERI_S_D2z_S_Pz_C1000001000_b += SQ_ERI_S_D_S_P_C1000001000_b_coefs*I_ERI_S_D2z_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_P_C1000001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_P_C1000001000_coefs = ic2_2*jc2_2;
      I_ERI_S_S_S_Px_C1000001000 += SQ_ERI_S_S_S_P_C1000001000_coefs*I_ERI_S_S_S_Px_vrr;
      I_ERI_S_S_S_Py_C1000001000 += SQ_ERI_S_S_S_P_C1000001000_coefs*I_ERI_S_S_S_Py_vrr;
      I_ERI_S_S_S_Pz_C1000001000 += SQ_ERI_S_S_S_P_C1000001000_coefs*I_ERI_S_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1000001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000001001_coefs = ic2_3*jc2_2;
      I_ERI_Px_S_S_Px_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000001001 += SQ_ERI_P_S_S_P_C1000001001_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C0_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C0_c_coefs = ic2*jc2*gamma;
      I_ERI_S_S_Px_S_C0_c += SQ_ERI_S_S_P_S_C0_c_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C0_c += SQ_ERI_S_S_P_S_C0_c_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C0_c += SQ_ERI_S_S_P_S_C0_c_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1_c_coefs = ic2_1*jc2*gamma;
      I_ERI_Px_S_Px_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1_c += SQ_ERI_P_S_P_S_C1_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1000_c_coefs = ic2_2*jc2*gamma;
      I_ERI_S_Px_Px_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1000_c += SQ_ERI_S_P_P_S_C1000_c_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000000_c_coefs = ic2*jc2_1*gamma;
      I_ERI_S_S_D2x_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1000000_c += SQ_ERI_S_S_D_S_C1000000_c_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S_C1000000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_S_C1000000_coefs = ic2*jc2_1;
      I_ERI_S_S_S_S_C1000000 += SQ_ERI_S_S_S_S_C1000000_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_c_coefs = ic2_1*jc2_1*gamma;
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
      Double SQ_ERI_P_S_S_S_C1000001_coefs = ic2_1*jc2_1;
      I_ERI_Px_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001000_c_coefs = ic2_2*jc2_1*gamma;
      I_ERI_S_Px_D2x_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001000_c += SQ_ERI_S_P_D_S_C1001000_c_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_S_C1001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_S_C1001000_coefs = ic2_2*jc2_1;
      I_ERI_S_Px_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_P_C1001000000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_S_P_C1001000000_coefs = ic2*jc2_3;
      I_ERI_S_S_S_Px_C1001000000 += SQ_ERI_S_S_S_P_C1001000000_coefs*I_ERI_S_S_S_Px_vrr;
      I_ERI_S_S_S_Py_C1001000000 += SQ_ERI_S_S_S_P_C1001000000_coefs*I_ERI_S_S_S_Py_vrr;
      I_ERI_S_S_S_Pz_C1001000000 += SQ_ERI_S_S_S_P_C1001000000_coefs*I_ERI_S_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1001000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1001000001_coefs = ic2_1*jc2_3;
      I_ERI_Px_S_S_Px_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1001000001 += SQ_ERI_P_S_S_P_C1001000001_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_S_P_C1001001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_S_P_C1001001000_coefs = ic2_2*jc2_3;
      I_ERI_S_Px_S_Px_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Px_S_Px_vrr;
      I_ERI_S_Py_S_Px_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Py_S_Px_vrr;
      I_ERI_S_Pz_S_Px_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Pz_S_Px_vrr;
      I_ERI_S_Px_S_Py_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Px_S_Py_vrr;
      I_ERI_S_Py_S_Py_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Py_S_Py_vrr;
      I_ERI_S_Pz_S_Py_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Pz_S_Py_vrr;
      I_ERI_S_Px_S_Pz_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Px_S_Pz_vrr;
      I_ERI_S_Py_S_Pz_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Py_S_Pz_vrr;
      I_ERI_S_Pz_S_Pz_C1001001000 += SQ_ERI_S_P_S_P_C1001001000_coefs*I_ERI_S_Pz_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1000_a_coefs = ic2_2*jc2*alpha;
      I_ERI_D2x_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1000_a += SQ_ERI_D_S_S_S_C1000_a_coefs*I_ERI_D2z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_C1001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001_a_coefs = ic2_3*jc2*alpha;
      I_ERI_F3x_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001_a += SQ_ERI_F_S_S_S_C1001_a_coefs*I_ERI_F3z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_C1001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001000_a_coefs = ic2_2*jc2_1*alpha;
      I_ERI_D2x_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001000_a += SQ_ERI_D_S_P_S_C1001000_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001000_a_coefs = ic2_2*jc2_1*alpha;
      I_ERI_Px_S_Px_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001000_a += SQ_ERI_P_S_P_S_C1001000_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_a_coefs = ic2_3*jc2_1*alpha;
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
      Double SQ_ERI_D_S_P_S_C1001001_a_coefs = ic2_3*jc2_1*alpha;
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
       * shell quartet name: SQ_ERI_D_S_S_P_C1000001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000001000_a_coefs = ic2_2*jc2_2*alpha;
      I_ERI_D2x_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000001000_a += SQ_ERI_D_S_S_P_C1000001000_a_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1000001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000001000_a_coefs = ic2_2*jc2_2*alpha;
      I_ERI_Px_S_S_Px_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000001000_a += SQ_ERI_P_S_S_P_C1000001000_a_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_P_C1000001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_P_C1000001001_a_coefs = ic2_3*jc2_2*alpha;
      I_ERI_F3x_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3x_S_S_Px_vrr;
      I_ERI_F2xy_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xy_S_S_Px_vrr;
      I_ERI_F2xz_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xz_S_S_Px_vrr;
      I_ERI_Fx2y_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2y_S_S_Px_vrr;
      I_ERI_Fxyz_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fxyz_S_S_Px_vrr;
      I_ERI_Fx2z_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2z_S_S_Px_vrr;
      I_ERI_F3y_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3y_S_S_Px_vrr;
      I_ERI_F2yz_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2yz_S_S_Px_vrr;
      I_ERI_Fy2z_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fy2z_S_S_Px_vrr;
      I_ERI_F3z_S_S_Px_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3z_S_S_Px_vrr;
      I_ERI_F3x_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3x_S_S_Py_vrr;
      I_ERI_F2xy_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xy_S_S_Py_vrr;
      I_ERI_F2xz_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xz_S_S_Py_vrr;
      I_ERI_Fx2y_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2y_S_S_Py_vrr;
      I_ERI_Fxyz_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fxyz_S_S_Py_vrr;
      I_ERI_Fx2z_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2z_S_S_Py_vrr;
      I_ERI_F3y_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3y_S_S_Py_vrr;
      I_ERI_F2yz_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2yz_S_S_Py_vrr;
      I_ERI_Fy2z_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fy2z_S_S_Py_vrr;
      I_ERI_F3z_S_S_Py_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3z_S_S_Py_vrr;
      I_ERI_F3x_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3x_S_S_Pz_vrr;
      I_ERI_F2xy_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xy_S_S_Pz_vrr;
      I_ERI_F2xz_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2xz_S_S_Pz_vrr;
      I_ERI_Fx2y_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2y_S_S_Pz_vrr;
      I_ERI_Fxyz_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fxyz_S_S_Pz_vrr;
      I_ERI_Fx2z_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fx2z_S_S_Pz_vrr;
      I_ERI_F3y_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3y_S_S_Pz_vrr;
      I_ERI_F2yz_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F2yz_S_S_Pz_vrr;
      I_ERI_Fy2z_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_Fy2z_S_S_Pz_vrr;
      I_ERI_F3z_S_S_Pz_C1000001001_a += SQ_ERI_F_S_S_P_C1000001001_a_coefs*I_ERI_F3z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1000001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000001001_a_coefs = ic2_3*jc2_2*alpha;
      I_ERI_D2x_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000001001_a += SQ_ERI_D_S_S_P_C1000001001_a_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1_b_coefs = ic2_1*jc2*beta;
      I_ERI_D2x_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1_b += SQ_ERI_D_S_S_S_C1_b_coefs*I_ERI_D2z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001_b_coefs = ic2_1*jc2_1*beta;
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
      Double SQ_ERI_P_S_P_S_C1000001_b_coefs = ic2_1*jc2_1*beta;
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
       * shell quartet name: SQ_ERI_D_S_S_P_C1000000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000000001_b_coefs = ic2_1*jc2_2*beta;
      I_ERI_D2x_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000000001_b += SQ_ERI_D_S_S_P_C1000000001_b_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000000001_b_coefs = ic2_1*jc2_2*beta;
      I_ERI_Px_S_S_Px_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000000001_b += SQ_ERI_P_S_S_P_C1000000001_b_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001_c_coefs = ic2_3*jc2*gamma;
      I_ERI_D2x_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001_c += SQ_ERI_D_S_P_S_C1001_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_c_coefs = ic2_3*jc2_1*gamma;
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
      Double SQ_ERI_P_S_D_S_C1001001_c_coefs = ic2_3*jc2_1*gamma;
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
      Double SQ_ERI_D_S_S_S_C1001001_coefs = ic2_3*jc2_1;
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
      Double SQ_ERI_P_S_S_S_C1001001_coefs = ic2_3*jc2_1;
      I_ERI_Px_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_D2x_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1001001001 += SQ_ERI_D_S_S_P_C1001001001_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_Px_S_S_Px_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1001001001 += SQ_ERI_P_S_S_P_C1001001001_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001_b_coefs = ic2_3*jc2*beta;
      I_ERI_F3x_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001_b += SQ_ERI_F_S_S_S_C1001_b_coefs*I_ERI_F3z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001_b_coefs = ic2_3*jc2_1*beta;
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
      Double SQ_ERI_D_S_P_S_C1001001_b_coefs = ic2_3*jc2_1*beta;
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
      Double SQ_ERI_P_S_P_S_C1001001_b_coefs = ic2_3*jc2_1*beta;
      I_ERI_Px_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_b += SQ_ERI_P_S_P_S_C1001001_b_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_P_C1000001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_P_C1000001001_b_coefs = ic2_3*jc2_2*beta;
      I_ERI_F3x_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3x_S_S_Px_vrr;
      I_ERI_F2xy_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xy_S_S_Px_vrr;
      I_ERI_F2xz_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xz_S_S_Px_vrr;
      I_ERI_Fx2y_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2y_S_S_Px_vrr;
      I_ERI_Fxyz_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fxyz_S_S_Px_vrr;
      I_ERI_Fx2z_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2z_S_S_Px_vrr;
      I_ERI_F3y_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3y_S_S_Px_vrr;
      I_ERI_F2yz_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2yz_S_S_Px_vrr;
      I_ERI_Fy2z_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fy2z_S_S_Px_vrr;
      I_ERI_F3z_S_S_Px_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3z_S_S_Px_vrr;
      I_ERI_F3x_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3x_S_S_Py_vrr;
      I_ERI_F2xy_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xy_S_S_Py_vrr;
      I_ERI_F2xz_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xz_S_S_Py_vrr;
      I_ERI_Fx2y_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2y_S_S_Py_vrr;
      I_ERI_Fxyz_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fxyz_S_S_Py_vrr;
      I_ERI_Fx2z_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2z_S_S_Py_vrr;
      I_ERI_F3y_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3y_S_S_Py_vrr;
      I_ERI_F2yz_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2yz_S_S_Py_vrr;
      I_ERI_Fy2z_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fy2z_S_S_Py_vrr;
      I_ERI_F3z_S_S_Py_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3z_S_S_Py_vrr;
      I_ERI_F3x_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3x_S_S_Pz_vrr;
      I_ERI_F2xy_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xy_S_S_Pz_vrr;
      I_ERI_F2xz_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2xz_S_S_Pz_vrr;
      I_ERI_Fx2y_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2y_S_S_Pz_vrr;
      I_ERI_Fxyz_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fxyz_S_S_Pz_vrr;
      I_ERI_Fx2z_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fx2z_S_S_Pz_vrr;
      I_ERI_F3y_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3y_S_S_Pz_vrr;
      I_ERI_F2yz_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F2yz_S_S_Pz_vrr;
      I_ERI_Fy2z_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_Fy2z_S_S_Pz_vrr;
      I_ERI_F3z_S_S_Pz_C1000001001_b += SQ_ERI_F_S_S_P_C1000001001_b_coefs*I_ERI_F3z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1000001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000001001_b_coefs = ic2_3*jc2_2*beta;
      I_ERI_D2x_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000001001_b += SQ_ERI_D_S_S_P_C1000001001_b_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P_C1000001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_P_C1000001001_b_coefs = ic2_3*jc2_2*beta;
      I_ERI_Px_S_S_Px_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz_C1000001001_b += SQ_ERI_P_S_S_P_C1000001001_b_coefs*I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001000000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001000000_a_coefs = ic2*jc2_3*alpha;
      I_ERI_Px_S_D2x_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001000000_a += SQ_ERI_P_S_D_S_C1001000000_a_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001000000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001000000_a_coefs = ic2*jc2_3*alpha;
      I_ERI_Px_S_Px_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001000000_a += SQ_ERI_P_S_P_S_C1001000000_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001000001_a_coefs = ic2_1*jc2_3*alpha;
      I_ERI_D2x_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001000001_a += SQ_ERI_D_S_D_S_C1001000001_a_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001000001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001000001_a_coefs = ic2_1*jc2_3*alpha;
      I_ERI_D2x_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001000001_a += SQ_ERI_D_S_P_S_C1001000001_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1001000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1001000001_coefs = ic2_1*jc2_3;
      I_ERI_S_S_D2x_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1001000001 += SQ_ERI_S_S_D_S_C1001000001_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1001000001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1001000001_coefs = ic2_1*jc2_3;
      I_ERI_S_S_Px_S_C1001000001 += SQ_ERI_S_S_P_S_C1001000001_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1001000001 += SQ_ERI_S_S_P_S_C1001000001_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1001000001 += SQ_ERI_S_S_P_S_C1001000001_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_S_Px_D2x_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001001001 += SQ_ERI_S_P_D_S_C1001001001_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_S_Px_Px_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001001001 += SQ_ERI_S_P_P_S_C1001001001_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001000000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001000000_b_coefs = ic2*jc2_3*beta;
      I_ERI_S_Px_D2x_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001000000_b += SQ_ERI_S_P_D_S_C1001000000_b_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001000000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001000000_b_coefs = ic2*jc2_3*beta;
      I_ERI_S_Px_Px_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001000000_b += SQ_ERI_S_P_P_S_C1001000000_b_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_D_S_C1001001000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_D_S_C1001001000_b_coefs = ic2_2*jc2_3*beta;
      I_ERI_S_D2x_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_D2x_S_vrr;
      I_ERI_S_Dxy_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_D2x_S_vrr;
      I_ERI_S_Dxz_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_D2x_S_vrr;
      I_ERI_S_D2y_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_D2x_S_vrr;
      I_ERI_S_Dyz_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_D2x_S_vrr;
      I_ERI_S_D2z_D2x_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_D2x_S_vrr;
      I_ERI_S_D2x_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_Dxy_S_vrr;
      I_ERI_S_Dxy_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_Dxy_S_vrr;
      I_ERI_S_Dxz_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_Dxy_S_vrr;
      I_ERI_S_D2y_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_Dxy_S_vrr;
      I_ERI_S_Dyz_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_Dxy_S_vrr;
      I_ERI_S_D2z_Dxy_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_Dxy_S_vrr;
      I_ERI_S_D2x_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_Dxz_S_vrr;
      I_ERI_S_Dxy_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_Dxz_S_vrr;
      I_ERI_S_Dxz_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_Dxz_S_vrr;
      I_ERI_S_D2y_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_Dxz_S_vrr;
      I_ERI_S_Dyz_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_Dxz_S_vrr;
      I_ERI_S_D2z_Dxz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_Dxz_S_vrr;
      I_ERI_S_D2x_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_D2y_S_vrr;
      I_ERI_S_Dxy_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_D2y_S_vrr;
      I_ERI_S_Dxz_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_D2y_S_vrr;
      I_ERI_S_D2y_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_D2y_S_vrr;
      I_ERI_S_Dyz_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_D2y_S_vrr;
      I_ERI_S_D2z_D2y_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_D2y_S_vrr;
      I_ERI_S_D2x_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_Dyz_S_vrr;
      I_ERI_S_Dxy_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_Dyz_S_vrr;
      I_ERI_S_Dxz_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_Dyz_S_vrr;
      I_ERI_S_D2y_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_Dyz_S_vrr;
      I_ERI_S_Dyz_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_Dyz_S_vrr;
      I_ERI_S_D2z_Dyz_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_Dyz_S_vrr;
      I_ERI_S_D2x_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2x_D2z_S_vrr;
      I_ERI_S_Dxy_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxy_D2z_S_vrr;
      I_ERI_S_Dxz_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dxz_D2z_S_vrr;
      I_ERI_S_D2y_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2y_D2z_S_vrr;
      I_ERI_S_Dyz_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_Dyz_D2z_S_vrr;
      I_ERI_S_D2z_D2z_S_C1001001000_b += SQ_ERI_S_D_D_S_C1001001000_b_coefs*I_ERI_S_D2z_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_D_P_S_C1001001000_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_D_P_S_C1001001000_b_coefs = ic2_2*jc2_3*beta;
      I_ERI_S_D2x_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2x_Px_S_vrr;
      I_ERI_S_Dxy_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxy_Px_S_vrr;
      I_ERI_S_Dxz_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxz_Px_S_vrr;
      I_ERI_S_D2y_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2y_Px_S_vrr;
      I_ERI_S_Dyz_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dyz_Px_S_vrr;
      I_ERI_S_D2z_Px_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2z_Px_S_vrr;
      I_ERI_S_D2x_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2x_Py_S_vrr;
      I_ERI_S_Dxy_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxy_Py_S_vrr;
      I_ERI_S_Dxz_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxz_Py_S_vrr;
      I_ERI_S_D2y_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2y_Py_S_vrr;
      I_ERI_S_Dyz_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dyz_Py_S_vrr;
      I_ERI_S_D2z_Py_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2z_Py_S_vrr;
      I_ERI_S_D2x_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2x_Pz_S_vrr;
      I_ERI_S_Dxy_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxy_Pz_S_vrr;
      I_ERI_S_Dxz_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dxz_Pz_S_vrr;
      I_ERI_S_D2y_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2y_Pz_S_vrr;
      I_ERI_S_Dyz_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_Dyz_Pz_S_vrr;
      I_ERI_S_D2z_Pz_S_C1001001000_b += SQ_ERI_S_D_P_S_C1001001000_b_coefs*I_ERI_S_D2z_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1001001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1001001000_coefs = ic2_2*jc2_3;
      I_ERI_S_S_D2x_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1001001000 += SQ_ERI_S_S_D_S_C1001001000_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1001001000
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1001001000_coefs = ic2_2*jc2_3;
      I_ERI_S_S_Px_S_C1001001000 += SQ_ERI_S_S_P_S_C1001001000_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1001001000 += SQ_ERI_S_S_P_S_C1001001000_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1001001000 += SQ_ERI_S_S_P_S_C1001001000_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_Px_S_D2x_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001001 += SQ_ERI_P_S_D_S_C1001001001_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001001_coefs = ic2_3*jc2_3;
      I_ERI_Px_S_Px_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001001 += SQ_ERI_P_S_P_S_C1001001001_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1000000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000000000_c_coefs = ic2*jc2_2*gamma;
      I_ERI_S_S_D2x_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1000000000_c += SQ_ERI_S_S_D_S_C1000000000_c_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1000000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1000000000_c_coefs = ic2*jc2_2*gamma;
      I_ERI_S_S_Px_S_C1000000000_c += SQ_ERI_S_S_P_S_C1000000000_c_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000000000_c += SQ_ERI_S_S_P_S_C1000000000_c_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000000000_c += SQ_ERI_S_S_P_S_C1000000000_c_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000000001_c_coefs = ic2_1*jc2_2*gamma;
      I_ERI_Px_S_D2x_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000000001_c += SQ_ERI_P_S_D_S_C1000000001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000000001_c_coefs = ic2_1*jc2_2*gamma;
      I_ERI_Px_S_Px_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000000001_c += SQ_ERI_P_S_P_S_C1000000001_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1000001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1000001000_c_coefs = ic2_2*jc2_2*gamma;
      I_ERI_S_Px_D2x_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1000001000_c += SQ_ERI_S_P_D_S_C1000001000_c_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1000001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1000001000_c_coefs = ic2_2*jc2_2*gamma;
      I_ERI_S_Px_Px_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1000001000_c += SQ_ERI_S_P_P_S_C1000001000_c_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_C1001000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_F_S_C1001000000_c_coefs = ic2*jc2_3*gamma;
      I_ERI_S_S_F3x_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F3x_S_vrr;
      I_ERI_S_S_F2xy_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F2xy_S_vrr;
      I_ERI_S_S_F2xz_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F2xz_S_vrr;
      I_ERI_S_S_Fx2y_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_Fx2y_S_vrr;
      I_ERI_S_S_Fxyz_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_Fxyz_S_vrr;
      I_ERI_S_S_Fx2z_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_Fx2z_S_vrr;
      I_ERI_S_S_F3y_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F3y_S_vrr;
      I_ERI_S_S_F2yz_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F2yz_S_vrr;
      I_ERI_S_S_Fy2z_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_Fy2z_S_vrr;
      I_ERI_S_S_F3z_S_C1001000000_c += SQ_ERI_S_S_F_S_C1001000000_c_coefs*I_ERI_S_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1001000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1001000000_c_coefs = ic2*jc2_3*gamma;
      I_ERI_S_S_D2x_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1001000000_c += SQ_ERI_S_S_D_S_C1001000000_c_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1001000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1001000001_c_coefs = ic2_1*jc2_3*gamma;
      I_ERI_Px_S_F3x_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1001000001_c += SQ_ERI_P_S_F_S_C1001000001_c_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001000001_c_coefs = ic2_1*jc2_3*gamma;
      I_ERI_Px_S_D2x_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001000001_c += SQ_ERI_P_S_D_S_C1001000001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_F_S_C1001001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_F_S_C1001001000_c_coefs = ic2_2*jc2_3*gamma;
      I_ERI_S_Px_F3x_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F3x_S_vrr;
      I_ERI_S_Py_F3x_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F3x_S_vrr;
      I_ERI_S_Pz_F3x_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F3x_S_vrr;
      I_ERI_S_Px_F2xy_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F2xy_S_vrr;
      I_ERI_S_Py_F2xy_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F2xy_S_vrr;
      I_ERI_S_Pz_F2xy_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F2xy_S_vrr;
      I_ERI_S_Px_F2xz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F2xz_S_vrr;
      I_ERI_S_Py_F2xz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F2xz_S_vrr;
      I_ERI_S_Pz_F2xz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F2xz_S_vrr;
      I_ERI_S_Px_Fx2y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_Fx2y_S_vrr;
      I_ERI_S_Py_Fx2y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_Fx2y_S_vrr;
      I_ERI_S_Pz_Fx2y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_Fx2y_S_vrr;
      I_ERI_S_Px_Fxyz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_Fxyz_S_vrr;
      I_ERI_S_Py_Fxyz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_Fxyz_S_vrr;
      I_ERI_S_Pz_Fxyz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_Fxyz_S_vrr;
      I_ERI_S_Px_Fx2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_Fx2z_S_vrr;
      I_ERI_S_Py_Fx2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_Fx2z_S_vrr;
      I_ERI_S_Pz_Fx2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_Fx2z_S_vrr;
      I_ERI_S_Px_F3y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F3y_S_vrr;
      I_ERI_S_Py_F3y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F3y_S_vrr;
      I_ERI_S_Pz_F3y_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F3y_S_vrr;
      I_ERI_S_Px_F2yz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F2yz_S_vrr;
      I_ERI_S_Py_F2yz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F2yz_S_vrr;
      I_ERI_S_Pz_F2yz_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F2yz_S_vrr;
      I_ERI_S_Px_Fy2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_Fy2z_S_vrr;
      I_ERI_S_Py_Fy2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_Fy2z_S_vrr;
      I_ERI_S_Pz_Fy2z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_Fy2z_S_vrr;
      I_ERI_S_Px_F3z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Px_F3z_S_vrr;
      I_ERI_S_Py_F3z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Py_F3z_S_vrr;
      I_ERI_S_Pz_F3z_S_C1001001000_c += SQ_ERI_S_P_F_S_C1001001000_c_coefs*I_ERI_S_Pz_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001001000_c_coefs = ic2_2*jc2_3*gamma;
      I_ERI_S_Px_D2x_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001001000_c += SQ_ERI_S_P_D_S_C1001001000_c_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001000_a_coefs = ic2_2*jc2_3*alpha;
      I_ERI_D2x_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001000_a += SQ_ERI_D_S_D_S_C1001001000_a_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001000_a_coefs = ic2_2*jc2_3*alpha;
      I_ERI_D2x_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001000_a += SQ_ERI_D_S_P_S_C1001001000_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001000_a_coefs = ic2_2*jc2_3*alpha;
      I_ERI_Px_S_D2x_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001000_a += SQ_ERI_P_S_D_S_C1001001000_a_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001000_a_coefs = ic2_2*jc2_3*alpha;
      I_ERI_Px_S_Px_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001000_a += SQ_ERI_P_S_P_S_C1001001000_a_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001001001_a_coefs = ic2_3*jc2_3*alpha;
      I_ERI_F3x_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001001001_a += SQ_ERI_F_S_D_S_C1001001001_a_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001001_a_coefs = ic2_3*jc2_3*alpha;
      I_ERI_F3x_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001001_a += SQ_ERI_F_S_P_S_C1001001001_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001001_a_coefs = ic2_3*jc2_3*alpha;
      I_ERI_D2x_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001001_a += SQ_ERI_D_S_D_S_C1001001001_a_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001001_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001001_a_coefs = ic2_3*jc2_3*alpha;
      I_ERI_D2x_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001001_a += SQ_ERI_D_S_P_S_C1001001001_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001000001_b_coefs = ic2_1*jc2_3*beta;
      I_ERI_D2x_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001000001_b += SQ_ERI_D_S_D_S_C1001000001_b_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001000001_b_coefs = ic2_1*jc2_3*beta;
      I_ERI_D2x_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001000001_b += SQ_ERI_D_S_P_S_C1001000001_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001000001_b_coefs = ic2_1*jc2_3*beta;
      I_ERI_Px_S_D2x_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001000001_b += SQ_ERI_P_S_D_S_C1001000001_b_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001000001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001000001_b_coefs = ic2_1*jc2_3*beta;
      I_ERI_Px_S_Px_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001000001_b += SQ_ERI_P_S_P_S_C1001000001_b_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1000001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1000001001_c_coefs = ic2_3*jc2_2*gamma;
      I_ERI_D2x_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1000001001_c += SQ_ERI_D_S_D_S_C1000001001_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000001001_c_coefs = ic2_3*jc2_2*gamma;
      I_ERI_D2x_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000001001_c += SQ_ERI_D_S_P_S_C1000001001_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001001_c_coefs = ic2_3*jc2_2*gamma;
      I_ERI_Px_S_D2x_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000001001_c += SQ_ERI_P_S_D_S_C1000001001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001001_c_coefs = ic2_3*jc2_2*gamma;
      I_ERI_Px_S_Px_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000001001_c += SQ_ERI_P_S_P_S_C1000001001_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_C1001001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_C1001001001_c_coefs = ic2_3*jc2_3*gamma;
      I_ERI_D2x_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_C1001001001_c += SQ_ERI_D_S_F_S_C1001001001_c_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001001_c_coefs = ic2_3*jc2_3*gamma;
      I_ERI_D2x_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001001_c += SQ_ERI_D_S_D_S_C1001001001_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_C1001001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_C1001001001_c_coefs = ic2_3*jc2_3*gamma;
      I_ERI_Px_S_F3x_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_C1001001001_c += SQ_ERI_P_S_F_S_C1001001001_c_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001001_c_coefs = ic2_3*jc2_3*gamma;
      I_ERI_Px_S_D2x_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001001_c += SQ_ERI_P_S_D_S_C1001001001_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_F3x_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001001001_b += SQ_ERI_F_S_D_S_C1001001001_b_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_F3x_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001001001_b += SQ_ERI_F_S_P_S_C1001001001_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_D2x_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001001_b += SQ_ERI_D_S_D_S_C1001001001_b_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_D2x_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001001001_b += SQ_ERI_D_S_P_S_C1001001001_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_Px_S_D2x_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001001_b += SQ_ERI_P_S_D_S_C1001001001_b_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001001_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001001_b_coefs = ic2_3*jc2_3*beta;
      I_ERI_Px_S_Px_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001001_b += SQ_ERI_P_S_P_S_C1001001001_b_coefs*I_ERI_Pz_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000001
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1001000001
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1001000001
   ************************************************************/
  Double I_ERI_S_S_Px_Px_C1001000001 = I_ERI_S_S_D2x_S_C1001000001+CDX*I_ERI_S_S_Px_S_C1001000001;
  Double I_ERI_S_S_Py_Px_C1001000001 = I_ERI_S_S_Dxy_S_C1001000001+CDX*I_ERI_S_S_Py_S_C1001000001;
  Double I_ERI_S_S_Pz_Px_C1001000001 = I_ERI_S_S_Dxz_S_C1001000001+CDX*I_ERI_S_S_Pz_S_C1001000001;
  Double I_ERI_S_S_Px_Py_C1001000001 = I_ERI_S_S_Dxy_S_C1001000001+CDY*I_ERI_S_S_Px_S_C1001000001;
  Double I_ERI_S_S_Py_Py_C1001000001 = I_ERI_S_S_D2y_S_C1001000001+CDY*I_ERI_S_S_Py_S_C1001000001;
  Double I_ERI_S_S_Pz_Py_C1001000001 = I_ERI_S_S_Dyz_S_C1001000001+CDY*I_ERI_S_S_Pz_S_C1001000001;
  Double I_ERI_S_S_Px_Pz_C1001000001 = I_ERI_S_S_Dxz_S_C1001000001+CDZ*I_ERI_S_S_Px_S_C1001000001;
  Double I_ERI_S_S_Py_Pz_C1001000001 = I_ERI_S_S_Dyz_S_C1001000001+CDZ*I_ERI_S_S_Py_S_C1001000001;
  Double I_ERI_S_S_Pz_Pz_C1001000001 = I_ERI_S_S_D2z_S_C1001000001+CDZ*I_ERI_S_S_Pz_S_C1001000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001001000
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1001001000
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1001001000
   ************************************************************/
  Double I_ERI_S_S_Px_Px_C1001001000 = I_ERI_S_S_D2x_S_C1001001000+CDX*I_ERI_S_S_Px_S_C1001001000;
  Double I_ERI_S_S_Py_Px_C1001001000 = I_ERI_S_S_Dxy_S_C1001001000+CDX*I_ERI_S_S_Py_S_C1001001000;
  Double I_ERI_S_S_Pz_Px_C1001001000 = I_ERI_S_S_Dxz_S_C1001001000+CDX*I_ERI_S_S_Pz_S_C1001001000;
  Double I_ERI_S_S_Px_Py_C1001001000 = I_ERI_S_S_Dxy_S_C1001001000+CDY*I_ERI_S_S_Px_S_C1001001000;
  Double I_ERI_S_S_Py_Py_C1001001000 = I_ERI_S_S_D2y_S_C1001001000+CDY*I_ERI_S_S_Py_S_C1001001000;
  Double I_ERI_S_S_Pz_Py_C1001001000 = I_ERI_S_S_Dyz_S_C1001001000+CDY*I_ERI_S_S_Pz_S_C1001001000;
  Double I_ERI_S_S_Px_Pz_C1001001000 = I_ERI_S_S_Dxz_S_C1001001000+CDZ*I_ERI_S_S_Px_S_C1001001000;
  Double I_ERI_S_S_Py_Pz_C1001001000 = I_ERI_S_S_Dyz_S_C1001001000+CDZ*I_ERI_S_S_Py_S_C1001001000;
  Double I_ERI_S_S_Pz_Pz_C1001001000 = I_ERI_S_S_D2z_S_C1001001000+CDZ*I_ERI_S_S_Pz_S_C1001001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001001001
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001001
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001001
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001001001 = I_ERI_Px_S_D2x_S_C1001001001+CDX*I_ERI_Px_S_Px_S_C1001001001;
  Double I_ERI_Py_S_Px_Px_C1001001001 = I_ERI_Py_S_D2x_S_C1001001001+CDX*I_ERI_Py_S_Px_S_C1001001001;
  Double I_ERI_Pz_S_Px_Px_C1001001001 = I_ERI_Pz_S_D2x_S_C1001001001+CDX*I_ERI_Pz_S_Px_S_C1001001001;
  Double I_ERI_Px_S_Py_Px_C1001001001 = I_ERI_Px_S_Dxy_S_C1001001001+CDX*I_ERI_Px_S_Py_S_C1001001001;
  Double I_ERI_Py_S_Py_Px_C1001001001 = I_ERI_Py_S_Dxy_S_C1001001001+CDX*I_ERI_Py_S_Py_S_C1001001001;
  Double I_ERI_Pz_S_Py_Px_C1001001001 = I_ERI_Pz_S_Dxy_S_C1001001001+CDX*I_ERI_Pz_S_Py_S_C1001001001;
  Double I_ERI_Px_S_Pz_Px_C1001001001 = I_ERI_Px_S_Dxz_S_C1001001001+CDX*I_ERI_Px_S_Pz_S_C1001001001;
  Double I_ERI_Py_S_Pz_Px_C1001001001 = I_ERI_Py_S_Dxz_S_C1001001001+CDX*I_ERI_Py_S_Pz_S_C1001001001;
  Double I_ERI_Pz_S_Pz_Px_C1001001001 = I_ERI_Pz_S_Dxz_S_C1001001001+CDX*I_ERI_Pz_S_Pz_S_C1001001001;
  Double I_ERI_Px_S_Px_Py_C1001001001 = I_ERI_Px_S_Dxy_S_C1001001001+CDY*I_ERI_Px_S_Px_S_C1001001001;
  Double I_ERI_Py_S_Px_Py_C1001001001 = I_ERI_Py_S_Dxy_S_C1001001001+CDY*I_ERI_Py_S_Px_S_C1001001001;
  Double I_ERI_Pz_S_Px_Py_C1001001001 = I_ERI_Pz_S_Dxy_S_C1001001001+CDY*I_ERI_Pz_S_Px_S_C1001001001;
  Double I_ERI_Px_S_Py_Py_C1001001001 = I_ERI_Px_S_D2y_S_C1001001001+CDY*I_ERI_Px_S_Py_S_C1001001001;
  Double I_ERI_Py_S_Py_Py_C1001001001 = I_ERI_Py_S_D2y_S_C1001001001+CDY*I_ERI_Py_S_Py_S_C1001001001;
  Double I_ERI_Pz_S_Py_Py_C1001001001 = I_ERI_Pz_S_D2y_S_C1001001001+CDY*I_ERI_Pz_S_Py_S_C1001001001;
  Double I_ERI_Px_S_Pz_Py_C1001001001 = I_ERI_Px_S_Dyz_S_C1001001001+CDY*I_ERI_Px_S_Pz_S_C1001001001;
  Double I_ERI_Py_S_Pz_Py_C1001001001 = I_ERI_Py_S_Dyz_S_C1001001001+CDY*I_ERI_Py_S_Pz_S_C1001001001;
  Double I_ERI_Pz_S_Pz_Py_C1001001001 = I_ERI_Pz_S_Dyz_S_C1001001001+CDY*I_ERI_Pz_S_Pz_S_C1001001001;
  Double I_ERI_Px_S_Px_Pz_C1001001001 = I_ERI_Px_S_Dxz_S_C1001001001+CDZ*I_ERI_Px_S_Px_S_C1001001001;
  Double I_ERI_Py_S_Px_Pz_C1001001001 = I_ERI_Py_S_Dxz_S_C1001001001+CDZ*I_ERI_Py_S_Px_S_C1001001001;
  Double I_ERI_Pz_S_Px_Pz_C1001001001 = I_ERI_Pz_S_Dxz_S_C1001001001+CDZ*I_ERI_Pz_S_Px_S_C1001001001;
  Double I_ERI_Px_S_Py_Pz_C1001001001 = I_ERI_Px_S_Dyz_S_C1001001001+CDZ*I_ERI_Px_S_Py_S_C1001001001;
  Double I_ERI_Py_S_Py_Pz_C1001001001 = I_ERI_Py_S_Dyz_S_C1001001001+CDZ*I_ERI_Py_S_Py_S_C1001001001;
  Double I_ERI_Pz_S_Py_Pz_C1001001001 = I_ERI_Pz_S_Dyz_S_C1001001001+CDZ*I_ERI_Pz_S_Py_S_C1001001001;
  Double I_ERI_Px_S_Pz_Pz_C1001001001 = I_ERI_Px_S_D2z_S_C1001001001+CDZ*I_ERI_Px_S_Pz_S_C1001001001;
  Double I_ERI_Py_S_Pz_Pz_C1001001001 = I_ERI_Py_S_D2z_S_C1001001001+CDZ*I_ERI_Py_S_Pz_S_C1001001001;
  Double I_ERI_Pz_S_Pz_Pz_C1001001001 = I_ERI_Pz_S_D2z_S_C1001001001+CDZ*I_ERI_Pz_S_Pz_S_C1001001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001001
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001001
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001001
   ************************************************************/
  Double I_ERI_S_Px_Px_Px_C1001001001 = I_ERI_S_Px_D2x_S_C1001001001+CDX*I_ERI_S_Px_Px_S_C1001001001;
  Double I_ERI_S_Py_Px_Px_C1001001001 = I_ERI_S_Py_D2x_S_C1001001001+CDX*I_ERI_S_Py_Px_S_C1001001001;
  Double I_ERI_S_Pz_Px_Px_C1001001001 = I_ERI_S_Pz_D2x_S_C1001001001+CDX*I_ERI_S_Pz_Px_S_C1001001001;
  Double I_ERI_S_Px_Py_Px_C1001001001 = I_ERI_S_Px_Dxy_S_C1001001001+CDX*I_ERI_S_Px_Py_S_C1001001001;
  Double I_ERI_S_Py_Py_Px_C1001001001 = I_ERI_S_Py_Dxy_S_C1001001001+CDX*I_ERI_S_Py_Py_S_C1001001001;
  Double I_ERI_S_Pz_Py_Px_C1001001001 = I_ERI_S_Pz_Dxy_S_C1001001001+CDX*I_ERI_S_Pz_Py_S_C1001001001;
  Double I_ERI_S_Px_Pz_Px_C1001001001 = I_ERI_S_Px_Dxz_S_C1001001001+CDX*I_ERI_S_Px_Pz_S_C1001001001;
  Double I_ERI_S_Py_Pz_Px_C1001001001 = I_ERI_S_Py_Dxz_S_C1001001001+CDX*I_ERI_S_Py_Pz_S_C1001001001;
  Double I_ERI_S_Pz_Pz_Px_C1001001001 = I_ERI_S_Pz_Dxz_S_C1001001001+CDX*I_ERI_S_Pz_Pz_S_C1001001001;
  Double I_ERI_S_Px_Px_Py_C1001001001 = I_ERI_S_Px_Dxy_S_C1001001001+CDY*I_ERI_S_Px_Px_S_C1001001001;
  Double I_ERI_S_Py_Px_Py_C1001001001 = I_ERI_S_Py_Dxy_S_C1001001001+CDY*I_ERI_S_Py_Px_S_C1001001001;
  Double I_ERI_S_Pz_Px_Py_C1001001001 = I_ERI_S_Pz_Dxy_S_C1001001001+CDY*I_ERI_S_Pz_Px_S_C1001001001;
  Double I_ERI_S_Px_Py_Py_C1001001001 = I_ERI_S_Px_D2y_S_C1001001001+CDY*I_ERI_S_Px_Py_S_C1001001001;
  Double I_ERI_S_Py_Py_Py_C1001001001 = I_ERI_S_Py_D2y_S_C1001001001+CDY*I_ERI_S_Py_Py_S_C1001001001;
  Double I_ERI_S_Pz_Py_Py_C1001001001 = I_ERI_S_Pz_D2y_S_C1001001001+CDY*I_ERI_S_Pz_Py_S_C1001001001;
  Double I_ERI_S_Px_Pz_Py_C1001001001 = I_ERI_S_Px_Dyz_S_C1001001001+CDY*I_ERI_S_Px_Pz_S_C1001001001;
  Double I_ERI_S_Py_Pz_Py_C1001001001 = I_ERI_S_Py_Dyz_S_C1001001001+CDY*I_ERI_S_Py_Pz_S_C1001001001;
  Double I_ERI_S_Pz_Pz_Py_C1001001001 = I_ERI_S_Pz_Dyz_S_C1001001001+CDY*I_ERI_S_Pz_Pz_S_C1001001001;
  Double I_ERI_S_Px_Px_Pz_C1001001001 = I_ERI_S_Px_Dxz_S_C1001001001+CDZ*I_ERI_S_Px_Px_S_C1001001001;
  Double I_ERI_S_Py_Px_Pz_C1001001001 = I_ERI_S_Py_Dxz_S_C1001001001+CDZ*I_ERI_S_Py_Px_S_C1001001001;
  Double I_ERI_S_Pz_Px_Pz_C1001001001 = I_ERI_S_Pz_Dxz_S_C1001001001+CDZ*I_ERI_S_Pz_Px_S_C1001001001;
  Double I_ERI_S_Px_Py_Pz_C1001001001 = I_ERI_S_Px_Dyz_S_C1001001001+CDZ*I_ERI_S_Px_Py_S_C1001001001;
  Double I_ERI_S_Py_Py_Pz_C1001001001 = I_ERI_S_Py_Dyz_S_C1001001001+CDZ*I_ERI_S_Py_Py_S_C1001001001;
  Double I_ERI_S_Pz_Py_Pz_C1001001001 = I_ERI_S_Pz_Dyz_S_C1001001001+CDZ*I_ERI_S_Pz_Py_S_C1001001001;
  Double I_ERI_S_Px_Pz_Pz_C1001001001 = I_ERI_S_Px_D2z_S_C1001001001+CDZ*I_ERI_S_Px_Pz_S_C1001001001;
  Double I_ERI_S_Py_Pz_Pz_C1001001001 = I_ERI_S_Py_D2z_S_C1001001001+CDZ*I_ERI_S_Py_Pz_S_C1001001001;
  Double I_ERI_S_Pz_Pz_Pz_C1001001001 = I_ERI_S_Pz_D2z_S_C1001001001+CDZ*I_ERI_S_Pz_Pz_S_C1001001001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000000_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001000000_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001000000_a
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001000000_a = I_ERI_Px_S_D2x_S_C1001000000_a+CDX*I_ERI_Px_S_Px_S_C1001000000_a;
  Double I_ERI_Py_S_Px_Px_C1001000000_a = I_ERI_Py_S_D2x_S_C1001000000_a+CDX*I_ERI_Py_S_Px_S_C1001000000_a;
  Double I_ERI_Pz_S_Px_Px_C1001000000_a = I_ERI_Pz_S_D2x_S_C1001000000_a+CDX*I_ERI_Pz_S_Px_S_C1001000000_a;
  Double I_ERI_Px_S_Py_Px_C1001000000_a = I_ERI_Px_S_Dxy_S_C1001000000_a+CDX*I_ERI_Px_S_Py_S_C1001000000_a;
  Double I_ERI_Py_S_Py_Px_C1001000000_a = I_ERI_Py_S_Dxy_S_C1001000000_a+CDX*I_ERI_Py_S_Py_S_C1001000000_a;
  Double I_ERI_Pz_S_Py_Px_C1001000000_a = I_ERI_Pz_S_Dxy_S_C1001000000_a+CDX*I_ERI_Pz_S_Py_S_C1001000000_a;
  Double I_ERI_Px_S_Pz_Px_C1001000000_a = I_ERI_Px_S_Dxz_S_C1001000000_a+CDX*I_ERI_Px_S_Pz_S_C1001000000_a;
  Double I_ERI_Py_S_Pz_Px_C1001000000_a = I_ERI_Py_S_Dxz_S_C1001000000_a+CDX*I_ERI_Py_S_Pz_S_C1001000000_a;
  Double I_ERI_Pz_S_Pz_Px_C1001000000_a = I_ERI_Pz_S_Dxz_S_C1001000000_a+CDX*I_ERI_Pz_S_Pz_S_C1001000000_a;
  Double I_ERI_Px_S_Px_Py_C1001000000_a = I_ERI_Px_S_Dxy_S_C1001000000_a+CDY*I_ERI_Px_S_Px_S_C1001000000_a;
  Double I_ERI_Py_S_Px_Py_C1001000000_a = I_ERI_Py_S_Dxy_S_C1001000000_a+CDY*I_ERI_Py_S_Px_S_C1001000000_a;
  Double I_ERI_Pz_S_Px_Py_C1001000000_a = I_ERI_Pz_S_Dxy_S_C1001000000_a+CDY*I_ERI_Pz_S_Px_S_C1001000000_a;
  Double I_ERI_Px_S_Py_Py_C1001000000_a = I_ERI_Px_S_D2y_S_C1001000000_a+CDY*I_ERI_Px_S_Py_S_C1001000000_a;
  Double I_ERI_Py_S_Py_Py_C1001000000_a = I_ERI_Py_S_D2y_S_C1001000000_a+CDY*I_ERI_Py_S_Py_S_C1001000000_a;
  Double I_ERI_Pz_S_Py_Py_C1001000000_a = I_ERI_Pz_S_D2y_S_C1001000000_a+CDY*I_ERI_Pz_S_Py_S_C1001000000_a;
  Double I_ERI_Px_S_Pz_Py_C1001000000_a = I_ERI_Px_S_Dyz_S_C1001000000_a+CDY*I_ERI_Px_S_Pz_S_C1001000000_a;
  Double I_ERI_Py_S_Pz_Py_C1001000000_a = I_ERI_Py_S_Dyz_S_C1001000000_a+CDY*I_ERI_Py_S_Pz_S_C1001000000_a;
  Double I_ERI_Pz_S_Pz_Py_C1001000000_a = I_ERI_Pz_S_Dyz_S_C1001000000_a+CDY*I_ERI_Pz_S_Pz_S_C1001000000_a;
  Double I_ERI_Px_S_Px_Pz_C1001000000_a = I_ERI_Px_S_Dxz_S_C1001000000_a+CDZ*I_ERI_Px_S_Px_S_C1001000000_a;
  Double I_ERI_Py_S_Px_Pz_C1001000000_a = I_ERI_Py_S_Dxz_S_C1001000000_a+CDZ*I_ERI_Py_S_Px_S_C1001000000_a;
  Double I_ERI_Pz_S_Px_Pz_C1001000000_a = I_ERI_Pz_S_Dxz_S_C1001000000_a+CDZ*I_ERI_Pz_S_Px_S_C1001000000_a;
  Double I_ERI_Px_S_Py_Pz_C1001000000_a = I_ERI_Px_S_Dyz_S_C1001000000_a+CDZ*I_ERI_Px_S_Py_S_C1001000000_a;
  Double I_ERI_Py_S_Py_Pz_C1001000000_a = I_ERI_Py_S_Dyz_S_C1001000000_a+CDZ*I_ERI_Py_S_Py_S_C1001000000_a;
  Double I_ERI_Pz_S_Py_Pz_C1001000000_a = I_ERI_Pz_S_Dyz_S_C1001000000_a+CDZ*I_ERI_Pz_S_Py_S_C1001000000_a;
  Double I_ERI_Px_S_Pz_Pz_C1001000000_a = I_ERI_Px_S_D2z_S_C1001000000_a+CDZ*I_ERI_Px_S_Pz_S_C1001000000_a;
  Double I_ERI_Py_S_Pz_Pz_C1001000000_a = I_ERI_Py_S_D2z_S_C1001000000_a+CDZ*I_ERI_Py_S_Pz_S_C1001000000_a;
  Double I_ERI_Pz_S_Pz_Pz_C1001000000_a = I_ERI_Pz_S_D2z_S_C1001000000_a+CDZ*I_ERI_Pz_S_Pz_S_C1001000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001000001_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001000001_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001000001_a
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001000001_a = I_ERI_D2x_S_D2x_S_C1001000001_a+CDX*I_ERI_D2x_S_Px_S_C1001000001_a;
  Double I_ERI_Dxy_S_Px_Px_C1001000001_a = I_ERI_Dxy_S_D2x_S_C1001000001_a+CDX*I_ERI_Dxy_S_Px_S_C1001000001_a;
  Double I_ERI_Dxz_S_Px_Px_C1001000001_a = I_ERI_Dxz_S_D2x_S_C1001000001_a+CDX*I_ERI_Dxz_S_Px_S_C1001000001_a;
  Double I_ERI_D2y_S_Px_Px_C1001000001_a = I_ERI_D2y_S_D2x_S_C1001000001_a+CDX*I_ERI_D2y_S_Px_S_C1001000001_a;
  Double I_ERI_Dyz_S_Px_Px_C1001000001_a = I_ERI_Dyz_S_D2x_S_C1001000001_a+CDX*I_ERI_Dyz_S_Px_S_C1001000001_a;
  Double I_ERI_D2z_S_Px_Px_C1001000001_a = I_ERI_D2z_S_D2x_S_C1001000001_a+CDX*I_ERI_D2z_S_Px_S_C1001000001_a;
  Double I_ERI_D2x_S_Py_Px_C1001000001_a = I_ERI_D2x_S_Dxy_S_C1001000001_a+CDX*I_ERI_D2x_S_Py_S_C1001000001_a;
  Double I_ERI_Dxy_S_Py_Px_C1001000001_a = I_ERI_Dxy_S_Dxy_S_C1001000001_a+CDX*I_ERI_Dxy_S_Py_S_C1001000001_a;
  Double I_ERI_Dxz_S_Py_Px_C1001000001_a = I_ERI_Dxz_S_Dxy_S_C1001000001_a+CDX*I_ERI_Dxz_S_Py_S_C1001000001_a;
  Double I_ERI_D2y_S_Py_Px_C1001000001_a = I_ERI_D2y_S_Dxy_S_C1001000001_a+CDX*I_ERI_D2y_S_Py_S_C1001000001_a;
  Double I_ERI_Dyz_S_Py_Px_C1001000001_a = I_ERI_Dyz_S_Dxy_S_C1001000001_a+CDX*I_ERI_Dyz_S_Py_S_C1001000001_a;
  Double I_ERI_D2z_S_Py_Px_C1001000001_a = I_ERI_D2z_S_Dxy_S_C1001000001_a+CDX*I_ERI_D2z_S_Py_S_C1001000001_a;
  Double I_ERI_D2x_S_Pz_Px_C1001000001_a = I_ERI_D2x_S_Dxz_S_C1001000001_a+CDX*I_ERI_D2x_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxy_S_Pz_Px_C1001000001_a = I_ERI_Dxy_S_Dxz_S_C1001000001_a+CDX*I_ERI_Dxy_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxz_S_Pz_Px_C1001000001_a = I_ERI_Dxz_S_Dxz_S_C1001000001_a+CDX*I_ERI_Dxz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2y_S_Pz_Px_C1001000001_a = I_ERI_D2y_S_Dxz_S_C1001000001_a+CDX*I_ERI_D2y_S_Pz_S_C1001000001_a;
  Double I_ERI_Dyz_S_Pz_Px_C1001000001_a = I_ERI_Dyz_S_Dxz_S_C1001000001_a+CDX*I_ERI_Dyz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2z_S_Pz_Px_C1001000001_a = I_ERI_D2z_S_Dxz_S_C1001000001_a+CDX*I_ERI_D2z_S_Pz_S_C1001000001_a;
  Double I_ERI_D2x_S_Px_Py_C1001000001_a = I_ERI_D2x_S_Dxy_S_C1001000001_a+CDY*I_ERI_D2x_S_Px_S_C1001000001_a;
  Double I_ERI_Dxy_S_Px_Py_C1001000001_a = I_ERI_Dxy_S_Dxy_S_C1001000001_a+CDY*I_ERI_Dxy_S_Px_S_C1001000001_a;
  Double I_ERI_Dxz_S_Px_Py_C1001000001_a = I_ERI_Dxz_S_Dxy_S_C1001000001_a+CDY*I_ERI_Dxz_S_Px_S_C1001000001_a;
  Double I_ERI_D2y_S_Px_Py_C1001000001_a = I_ERI_D2y_S_Dxy_S_C1001000001_a+CDY*I_ERI_D2y_S_Px_S_C1001000001_a;
  Double I_ERI_Dyz_S_Px_Py_C1001000001_a = I_ERI_Dyz_S_Dxy_S_C1001000001_a+CDY*I_ERI_Dyz_S_Px_S_C1001000001_a;
  Double I_ERI_D2z_S_Px_Py_C1001000001_a = I_ERI_D2z_S_Dxy_S_C1001000001_a+CDY*I_ERI_D2z_S_Px_S_C1001000001_a;
  Double I_ERI_D2x_S_Py_Py_C1001000001_a = I_ERI_D2x_S_D2y_S_C1001000001_a+CDY*I_ERI_D2x_S_Py_S_C1001000001_a;
  Double I_ERI_Dxy_S_Py_Py_C1001000001_a = I_ERI_Dxy_S_D2y_S_C1001000001_a+CDY*I_ERI_Dxy_S_Py_S_C1001000001_a;
  Double I_ERI_Dxz_S_Py_Py_C1001000001_a = I_ERI_Dxz_S_D2y_S_C1001000001_a+CDY*I_ERI_Dxz_S_Py_S_C1001000001_a;
  Double I_ERI_D2y_S_Py_Py_C1001000001_a = I_ERI_D2y_S_D2y_S_C1001000001_a+CDY*I_ERI_D2y_S_Py_S_C1001000001_a;
  Double I_ERI_Dyz_S_Py_Py_C1001000001_a = I_ERI_Dyz_S_D2y_S_C1001000001_a+CDY*I_ERI_Dyz_S_Py_S_C1001000001_a;
  Double I_ERI_D2z_S_Py_Py_C1001000001_a = I_ERI_D2z_S_D2y_S_C1001000001_a+CDY*I_ERI_D2z_S_Py_S_C1001000001_a;
  Double I_ERI_D2x_S_Pz_Py_C1001000001_a = I_ERI_D2x_S_Dyz_S_C1001000001_a+CDY*I_ERI_D2x_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxy_S_Pz_Py_C1001000001_a = I_ERI_Dxy_S_Dyz_S_C1001000001_a+CDY*I_ERI_Dxy_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxz_S_Pz_Py_C1001000001_a = I_ERI_Dxz_S_Dyz_S_C1001000001_a+CDY*I_ERI_Dxz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2y_S_Pz_Py_C1001000001_a = I_ERI_D2y_S_Dyz_S_C1001000001_a+CDY*I_ERI_D2y_S_Pz_S_C1001000001_a;
  Double I_ERI_Dyz_S_Pz_Py_C1001000001_a = I_ERI_Dyz_S_Dyz_S_C1001000001_a+CDY*I_ERI_Dyz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2z_S_Pz_Py_C1001000001_a = I_ERI_D2z_S_Dyz_S_C1001000001_a+CDY*I_ERI_D2z_S_Pz_S_C1001000001_a;
  Double I_ERI_D2x_S_Px_Pz_C1001000001_a = I_ERI_D2x_S_Dxz_S_C1001000001_a+CDZ*I_ERI_D2x_S_Px_S_C1001000001_a;
  Double I_ERI_Dxy_S_Px_Pz_C1001000001_a = I_ERI_Dxy_S_Dxz_S_C1001000001_a+CDZ*I_ERI_Dxy_S_Px_S_C1001000001_a;
  Double I_ERI_Dxz_S_Px_Pz_C1001000001_a = I_ERI_Dxz_S_Dxz_S_C1001000001_a+CDZ*I_ERI_Dxz_S_Px_S_C1001000001_a;
  Double I_ERI_D2y_S_Px_Pz_C1001000001_a = I_ERI_D2y_S_Dxz_S_C1001000001_a+CDZ*I_ERI_D2y_S_Px_S_C1001000001_a;
  Double I_ERI_Dyz_S_Px_Pz_C1001000001_a = I_ERI_Dyz_S_Dxz_S_C1001000001_a+CDZ*I_ERI_Dyz_S_Px_S_C1001000001_a;
  Double I_ERI_D2z_S_Px_Pz_C1001000001_a = I_ERI_D2z_S_Dxz_S_C1001000001_a+CDZ*I_ERI_D2z_S_Px_S_C1001000001_a;
  Double I_ERI_D2x_S_Py_Pz_C1001000001_a = I_ERI_D2x_S_Dyz_S_C1001000001_a+CDZ*I_ERI_D2x_S_Py_S_C1001000001_a;
  Double I_ERI_Dxy_S_Py_Pz_C1001000001_a = I_ERI_Dxy_S_Dyz_S_C1001000001_a+CDZ*I_ERI_Dxy_S_Py_S_C1001000001_a;
  Double I_ERI_Dxz_S_Py_Pz_C1001000001_a = I_ERI_Dxz_S_Dyz_S_C1001000001_a+CDZ*I_ERI_Dxz_S_Py_S_C1001000001_a;
  Double I_ERI_D2y_S_Py_Pz_C1001000001_a = I_ERI_D2y_S_Dyz_S_C1001000001_a+CDZ*I_ERI_D2y_S_Py_S_C1001000001_a;
  Double I_ERI_Dyz_S_Py_Pz_C1001000001_a = I_ERI_Dyz_S_Dyz_S_C1001000001_a+CDZ*I_ERI_Dyz_S_Py_S_C1001000001_a;
  Double I_ERI_D2z_S_Py_Pz_C1001000001_a = I_ERI_D2z_S_Dyz_S_C1001000001_a+CDZ*I_ERI_D2z_S_Py_S_C1001000001_a;
  Double I_ERI_D2x_S_Pz_Pz_C1001000001_a = I_ERI_D2x_S_D2z_S_C1001000001_a+CDZ*I_ERI_D2x_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxy_S_Pz_Pz_C1001000001_a = I_ERI_Dxy_S_D2z_S_C1001000001_a+CDZ*I_ERI_Dxy_S_Pz_S_C1001000001_a;
  Double I_ERI_Dxz_S_Pz_Pz_C1001000001_a = I_ERI_Dxz_S_D2z_S_C1001000001_a+CDZ*I_ERI_Dxz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2y_S_Pz_Pz_C1001000001_a = I_ERI_D2y_S_D2z_S_C1001000001_a+CDZ*I_ERI_D2y_S_Pz_S_C1001000001_a;
  Double I_ERI_Dyz_S_Pz_Pz_C1001000001_a = I_ERI_Dyz_S_D2z_S_C1001000001_a+CDZ*I_ERI_Dyz_S_Pz_S_C1001000001_a;
  Double I_ERI_D2z_S_Pz_Pz_C1001000001_a = I_ERI_D2z_S_D2z_S_C1001000001_a+CDZ*I_ERI_D2z_S_Pz_S_C1001000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001001000_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001000_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001000_a
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001001000_a = I_ERI_Px_S_D2x_S_C1001001000_a+CDX*I_ERI_Px_S_Px_S_C1001001000_a;
  Double I_ERI_Py_S_Px_Px_C1001001000_a = I_ERI_Py_S_D2x_S_C1001001000_a+CDX*I_ERI_Py_S_Px_S_C1001001000_a;
  Double I_ERI_Pz_S_Px_Px_C1001001000_a = I_ERI_Pz_S_D2x_S_C1001001000_a+CDX*I_ERI_Pz_S_Px_S_C1001001000_a;
  Double I_ERI_Px_S_Py_Px_C1001001000_a = I_ERI_Px_S_Dxy_S_C1001001000_a+CDX*I_ERI_Px_S_Py_S_C1001001000_a;
  Double I_ERI_Py_S_Py_Px_C1001001000_a = I_ERI_Py_S_Dxy_S_C1001001000_a+CDX*I_ERI_Py_S_Py_S_C1001001000_a;
  Double I_ERI_Pz_S_Py_Px_C1001001000_a = I_ERI_Pz_S_Dxy_S_C1001001000_a+CDX*I_ERI_Pz_S_Py_S_C1001001000_a;
  Double I_ERI_Px_S_Pz_Px_C1001001000_a = I_ERI_Px_S_Dxz_S_C1001001000_a+CDX*I_ERI_Px_S_Pz_S_C1001001000_a;
  Double I_ERI_Py_S_Pz_Px_C1001001000_a = I_ERI_Py_S_Dxz_S_C1001001000_a+CDX*I_ERI_Py_S_Pz_S_C1001001000_a;
  Double I_ERI_Pz_S_Pz_Px_C1001001000_a = I_ERI_Pz_S_Dxz_S_C1001001000_a+CDX*I_ERI_Pz_S_Pz_S_C1001001000_a;
  Double I_ERI_Px_S_Px_Py_C1001001000_a = I_ERI_Px_S_Dxy_S_C1001001000_a+CDY*I_ERI_Px_S_Px_S_C1001001000_a;
  Double I_ERI_Py_S_Px_Py_C1001001000_a = I_ERI_Py_S_Dxy_S_C1001001000_a+CDY*I_ERI_Py_S_Px_S_C1001001000_a;
  Double I_ERI_Pz_S_Px_Py_C1001001000_a = I_ERI_Pz_S_Dxy_S_C1001001000_a+CDY*I_ERI_Pz_S_Px_S_C1001001000_a;
  Double I_ERI_Px_S_Py_Py_C1001001000_a = I_ERI_Px_S_D2y_S_C1001001000_a+CDY*I_ERI_Px_S_Py_S_C1001001000_a;
  Double I_ERI_Py_S_Py_Py_C1001001000_a = I_ERI_Py_S_D2y_S_C1001001000_a+CDY*I_ERI_Py_S_Py_S_C1001001000_a;
  Double I_ERI_Pz_S_Py_Py_C1001001000_a = I_ERI_Pz_S_D2y_S_C1001001000_a+CDY*I_ERI_Pz_S_Py_S_C1001001000_a;
  Double I_ERI_Px_S_Pz_Py_C1001001000_a = I_ERI_Px_S_Dyz_S_C1001001000_a+CDY*I_ERI_Px_S_Pz_S_C1001001000_a;
  Double I_ERI_Py_S_Pz_Py_C1001001000_a = I_ERI_Py_S_Dyz_S_C1001001000_a+CDY*I_ERI_Py_S_Pz_S_C1001001000_a;
  Double I_ERI_Pz_S_Pz_Py_C1001001000_a = I_ERI_Pz_S_Dyz_S_C1001001000_a+CDY*I_ERI_Pz_S_Pz_S_C1001001000_a;
  Double I_ERI_Px_S_Px_Pz_C1001001000_a = I_ERI_Px_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Px_S_Px_S_C1001001000_a;
  Double I_ERI_Py_S_Px_Pz_C1001001000_a = I_ERI_Py_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Py_S_Px_S_C1001001000_a;
  Double I_ERI_Pz_S_Px_Pz_C1001001000_a = I_ERI_Pz_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Pz_S_Px_S_C1001001000_a;
  Double I_ERI_Px_S_Py_Pz_C1001001000_a = I_ERI_Px_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Px_S_Py_S_C1001001000_a;
  Double I_ERI_Py_S_Py_Pz_C1001001000_a = I_ERI_Py_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Py_S_Py_S_C1001001000_a;
  Double I_ERI_Pz_S_Py_Pz_C1001001000_a = I_ERI_Pz_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Pz_S_Py_S_C1001001000_a;
  Double I_ERI_Px_S_Pz_Pz_C1001001000_a = I_ERI_Px_S_D2z_S_C1001001000_a+CDZ*I_ERI_Px_S_Pz_S_C1001001000_a;
  Double I_ERI_Py_S_Pz_Pz_C1001001000_a = I_ERI_Py_S_D2z_S_C1001001000_a+CDZ*I_ERI_Py_S_Pz_S_C1001001000_a;
  Double I_ERI_Pz_S_Pz_Pz_C1001001000_a = I_ERI_Pz_S_D2z_S_C1001001000_a+CDZ*I_ERI_Pz_S_Pz_S_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001000_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001000_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001000_a
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001000_a = I_ERI_D2x_S_D2x_S_C1001001000_a+CDX*I_ERI_D2x_S_Px_S_C1001001000_a;
  Double I_ERI_Dxy_S_Px_Px_C1001001000_a = I_ERI_Dxy_S_D2x_S_C1001001000_a+CDX*I_ERI_Dxy_S_Px_S_C1001001000_a;
  Double I_ERI_Dxz_S_Px_Px_C1001001000_a = I_ERI_Dxz_S_D2x_S_C1001001000_a+CDX*I_ERI_Dxz_S_Px_S_C1001001000_a;
  Double I_ERI_D2y_S_Px_Px_C1001001000_a = I_ERI_D2y_S_D2x_S_C1001001000_a+CDX*I_ERI_D2y_S_Px_S_C1001001000_a;
  Double I_ERI_Dyz_S_Px_Px_C1001001000_a = I_ERI_Dyz_S_D2x_S_C1001001000_a+CDX*I_ERI_Dyz_S_Px_S_C1001001000_a;
  Double I_ERI_D2z_S_Px_Px_C1001001000_a = I_ERI_D2z_S_D2x_S_C1001001000_a+CDX*I_ERI_D2z_S_Px_S_C1001001000_a;
  Double I_ERI_D2x_S_Py_Px_C1001001000_a = I_ERI_D2x_S_Dxy_S_C1001001000_a+CDX*I_ERI_D2x_S_Py_S_C1001001000_a;
  Double I_ERI_Dxy_S_Py_Px_C1001001000_a = I_ERI_Dxy_S_Dxy_S_C1001001000_a+CDX*I_ERI_Dxy_S_Py_S_C1001001000_a;
  Double I_ERI_Dxz_S_Py_Px_C1001001000_a = I_ERI_Dxz_S_Dxy_S_C1001001000_a+CDX*I_ERI_Dxz_S_Py_S_C1001001000_a;
  Double I_ERI_D2y_S_Py_Px_C1001001000_a = I_ERI_D2y_S_Dxy_S_C1001001000_a+CDX*I_ERI_D2y_S_Py_S_C1001001000_a;
  Double I_ERI_Dyz_S_Py_Px_C1001001000_a = I_ERI_Dyz_S_Dxy_S_C1001001000_a+CDX*I_ERI_Dyz_S_Py_S_C1001001000_a;
  Double I_ERI_D2z_S_Py_Px_C1001001000_a = I_ERI_D2z_S_Dxy_S_C1001001000_a+CDX*I_ERI_D2z_S_Py_S_C1001001000_a;
  Double I_ERI_D2x_S_Pz_Px_C1001001000_a = I_ERI_D2x_S_Dxz_S_C1001001000_a+CDX*I_ERI_D2x_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxy_S_Pz_Px_C1001001000_a = I_ERI_Dxy_S_Dxz_S_C1001001000_a+CDX*I_ERI_Dxy_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxz_S_Pz_Px_C1001001000_a = I_ERI_Dxz_S_Dxz_S_C1001001000_a+CDX*I_ERI_Dxz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2y_S_Pz_Px_C1001001000_a = I_ERI_D2y_S_Dxz_S_C1001001000_a+CDX*I_ERI_D2y_S_Pz_S_C1001001000_a;
  Double I_ERI_Dyz_S_Pz_Px_C1001001000_a = I_ERI_Dyz_S_Dxz_S_C1001001000_a+CDX*I_ERI_Dyz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2z_S_Pz_Px_C1001001000_a = I_ERI_D2z_S_Dxz_S_C1001001000_a+CDX*I_ERI_D2z_S_Pz_S_C1001001000_a;
  Double I_ERI_D2x_S_Px_Py_C1001001000_a = I_ERI_D2x_S_Dxy_S_C1001001000_a+CDY*I_ERI_D2x_S_Px_S_C1001001000_a;
  Double I_ERI_Dxy_S_Px_Py_C1001001000_a = I_ERI_Dxy_S_Dxy_S_C1001001000_a+CDY*I_ERI_Dxy_S_Px_S_C1001001000_a;
  Double I_ERI_Dxz_S_Px_Py_C1001001000_a = I_ERI_Dxz_S_Dxy_S_C1001001000_a+CDY*I_ERI_Dxz_S_Px_S_C1001001000_a;
  Double I_ERI_D2y_S_Px_Py_C1001001000_a = I_ERI_D2y_S_Dxy_S_C1001001000_a+CDY*I_ERI_D2y_S_Px_S_C1001001000_a;
  Double I_ERI_Dyz_S_Px_Py_C1001001000_a = I_ERI_Dyz_S_Dxy_S_C1001001000_a+CDY*I_ERI_Dyz_S_Px_S_C1001001000_a;
  Double I_ERI_D2z_S_Px_Py_C1001001000_a = I_ERI_D2z_S_Dxy_S_C1001001000_a+CDY*I_ERI_D2z_S_Px_S_C1001001000_a;
  Double I_ERI_D2x_S_Py_Py_C1001001000_a = I_ERI_D2x_S_D2y_S_C1001001000_a+CDY*I_ERI_D2x_S_Py_S_C1001001000_a;
  Double I_ERI_Dxy_S_Py_Py_C1001001000_a = I_ERI_Dxy_S_D2y_S_C1001001000_a+CDY*I_ERI_Dxy_S_Py_S_C1001001000_a;
  Double I_ERI_Dxz_S_Py_Py_C1001001000_a = I_ERI_Dxz_S_D2y_S_C1001001000_a+CDY*I_ERI_Dxz_S_Py_S_C1001001000_a;
  Double I_ERI_D2y_S_Py_Py_C1001001000_a = I_ERI_D2y_S_D2y_S_C1001001000_a+CDY*I_ERI_D2y_S_Py_S_C1001001000_a;
  Double I_ERI_Dyz_S_Py_Py_C1001001000_a = I_ERI_Dyz_S_D2y_S_C1001001000_a+CDY*I_ERI_Dyz_S_Py_S_C1001001000_a;
  Double I_ERI_D2z_S_Py_Py_C1001001000_a = I_ERI_D2z_S_D2y_S_C1001001000_a+CDY*I_ERI_D2z_S_Py_S_C1001001000_a;
  Double I_ERI_D2x_S_Pz_Py_C1001001000_a = I_ERI_D2x_S_Dyz_S_C1001001000_a+CDY*I_ERI_D2x_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxy_S_Pz_Py_C1001001000_a = I_ERI_Dxy_S_Dyz_S_C1001001000_a+CDY*I_ERI_Dxy_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxz_S_Pz_Py_C1001001000_a = I_ERI_Dxz_S_Dyz_S_C1001001000_a+CDY*I_ERI_Dxz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2y_S_Pz_Py_C1001001000_a = I_ERI_D2y_S_Dyz_S_C1001001000_a+CDY*I_ERI_D2y_S_Pz_S_C1001001000_a;
  Double I_ERI_Dyz_S_Pz_Py_C1001001000_a = I_ERI_Dyz_S_Dyz_S_C1001001000_a+CDY*I_ERI_Dyz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2z_S_Pz_Py_C1001001000_a = I_ERI_D2z_S_Dyz_S_C1001001000_a+CDY*I_ERI_D2z_S_Pz_S_C1001001000_a;
  Double I_ERI_D2x_S_Px_Pz_C1001001000_a = I_ERI_D2x_S_Dxz_S_C1001001000_a+CDZ*I_ERI_D2x_S_Px_S_C1001001000_a;
  Double I_ERI_Dxy_S_Px_Pz_C1001001000_a = I_ERI_Dxy_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Dxy_S_Px_S_C1001001000_a;
  Double I_ERI_Dxz_S_Px_Pz_C1001001000_a = I_ERI_Dxz_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Dxz_S_Px_S_C1001001000_a;
  Double I_ERI_D2y_S_Px_Pz_C1001001000_a = I_ERI_D2y_S_Dxz_S_C1001001000_a+CDZ*I_ERI_D2y_S_Px_S_C1001001000_a;
  Double I_ERI_Dyz_S_Px_Pz_C1001001000_a = I_ERI_Dyz_S_Dxz_S_C1001001000_a+CDZ*I_ERI_Dyz_S_Px_S_C1001001000_a;
  Double I_ERI_D2z_S_Px_Pz_C1001001000_a = I_ERI_D2z_S_Dxz_S_C1001001000_a+CDZ*I_ERI_D2z_S_Px_S_C1001001000_a;
  Double I_ERI_D2x_S_Py_Pz_C1001001000_a = I_ERI_D2x_S_Dyz_S_C1001001000_a+CDZ*I_ERI_D2x_S_Py_S_C1001001000_a;
  Double I_ERI_Dxy_S_Py_Pz_C1001001000_a = I_ERI_Dxy_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Dxy_S_Py_S_C1001001000_a;
  Double I_ERI_Dxz_S_Py_Pz_C1001001000_a = I_ERI_Dxz_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Dxz_S_Py_S_C1001001000_a;
  Double I_ERI_D2y_S_Py_Pz_C1001001000_a = I_ERI_D2y_S_Dyz_S_C1001001000_a+CDZ*I_ERI_D2y_S_Py_S_C1001001000_a;
  Double I_ERI_Dyz_S_Py_Pz_C1001001000_a = I_ERI_Dyz_S_Dyz_S_C1001001000_a+CDZ*I_ERI_Dyz_S_Py_S_C1001001000_a;
  Double I_ERI_D2z_S_Py_Pz_C1001001000_a = I_ERI_D2z_S_Dyz_S_C1001001000_a+CDZ*I_ERI_D2z_S_Py_S_C1001001000_a;
  Double I_ERI_D2x_S_Pz_Pz_C1001001000_a = I_ERI_D2x_S_D2z_S_C1001001000_a+CDZ*I_ERI_D2x_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001000_a = I_ERI_Dxy_S_D2z_S_C1001001000_a+CDZ*I_ERI_Dxy_S_Pz_S_C1001001000_a;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001000_a = I_ERI_Dxz_S_D2z_S_C1001001000_a+CDZ*I_ERI_Dxz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2y_S_Pz_Pz_C1001001000_a = I_ERI_D2y_S_D2z_S_C1001001000_a+CDZ*I_ERI_D2y_S_Pz_S_C1001001000_a;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001000_a = I_ERI_Dyz_S_D2z_S_C1001001000_a+CDZ*I_ERI_Dyz_S_Pz_S_C1001001000_a;
  Double I_ERI_D2z_S_Pz_Pz_C1001001000_a = I_ERI_D2z_S_D2z_S_C1001001000_a+CDZ*I_ERI_D2z_S_Pz_S_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001001_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001001_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001001_a
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001001_a = I_ERI_D2x_S_D2x_S_C1001001001_a+CDX*I_ERI_D2x_S_Px_S_C1001001001_a;
  Double I_ERI_Dxy_S_Px_Px_C1001001001_a = I_ERI_Dxy_S_D2x_S_C1001001001_a+CDX*I_ERI_Dxy_S_Px_S_C1001001001_a;
  Double I_ERI_Dxz_S_Px_Px_C1001001001_a = I_ERI_Dxz_S_D2x_S_C1001001001_a+CDX*I_ERI_Dxz_S_Px_S_C1001001001_a;
  Double I_ERI_D2y_S_Px_Px_C1001001001_a = I_ERI_D2y_S_D2x_S_C1001001001_a+CDX*I_ERI_D2y_S_Px_S_C1001001001_a;
  Double I_ERI_Dyz_S_Px_Px_C1001001001_a = I_ERI_Dyz_S_D2x_S_C1001001001_a+CDX*I_ERI_Dyz_S_Px_S_C1001001001_a;
  Double I_ERI_D2z_S_Px_Px_C1001001001_a = I_ERI_D2z_S_D2x_S_C1001001001_a+CDX*I_ERI_D2z_S_Px_S_C1001001001_a;
  Double I_ERI_D2x_S_Py_Px_C1001001001_a = I_ERI_D2x_S_Dxy_S_C1001001001_a+CDX*I_ERI_D2x_S_Py_S_C1001001001_a;
  Double I_ERI_Dxy_S_Py_Px_C1001001001_a = I_ERI_Dxy_S_Dxy_S_C1001001001_a+CDX*I_ERI_Dxy_S_Py_S_C1001001001_a;
  Double I_ERI_Dxz_S_Py_Px_C1001001001_a = I_ERI_Dxz_S_Dxy_S_C1001001001_a+CDX*I_ERI_Dxz_S_Py_S_C1001001001_a;
  Double I_ERI_D2y_S_Py_Px_C1001001001_a = I_ERI_D2y_S_Dxy_S_C1001001001_a+CDX*I_ERI_D2y_S_Py_S_C1001001001_a;
  Double I_ERI_Dyz_S_Py_Px_C1001001001_a = I_ERI_Dyz_S_Dxy_S_C1001001001_a+CDX*I_ERI_Dyz_S_Py_S_C1001001001_a;
  Double I_ERI_D2z_S_Py_Px_C1001001001_a = I_ERI_D2z_S_Dxy_S_C1001001001_a+CDX*I_ERI_D2z_S_Py_S_C1001001001_a;
  Double I_ERI_D2x_S_Pz_Px_C1001001001_a = I_ERI_D2x_S_Dxz_S_C1001001001_a+CDX*I_ERI_D2x_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxy_S_Pz_Px_C1001001001_a = I_ERI_Dxy_S_Dxz_S_C1001001001_a+CDX*I_ERI_Dxy_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxz_S_Pz_Px_C1001001001_a = I_ERI_Dxz_S_Dxz_S_C1001001001_a+CDX*I_ERI_Dxz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2y_S_Pz_Px_C1001001001_a = I_ERI_D2y_S_Dxz_S_C1001001001_a+CDX*I_ERI_D2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Dyz_S_Pz_Px_C1001001001_a = I_ERI_Dyz_S_Dxz_S_C1001001001_a+CDX*I_ERI_Dyz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2z_S_Pz_Px_C1001001001_a = I_ERI_D2z_S_Dxz_S_C1001001001_a+CDX*I_ERI_D2z_S_Pz_S_C1001001001_a;
  Double I_ERI_D2x_S_Px_Py_C1001001001_a = I_ERI_D2x_S_Dxy_S_C1001001001_a+CDY*I_ERI_D2x_S_Px_S_C1001001001_a;
  Double I_ERI_Dxy_S_Px_Py_C1001001001_a = I_ERI_Dxy_S_Dxy_S_C1001001001_a+CDY*I_ERI_Dxy_S_Px_S_C1001001001_a;
  Double I_ERI_Dxz_S_Px_Py_C1001001001_a = I_ERI_Dxz_S_Dxy_S_C1001001001_a+CDY*I_ERI_Dxz_S_Px_S_C1001001001_a;
  Double I_ERI_D2y_S_Px_Py_C1001001001_a = I_ERI_D2y_S_Dxy_S_C1001001001_a+CDY*I_ERI_D2y_S_Px_S_C1001001001_a;
  Double I_ERI_Dyz_S_Px_Py_C1001001001_a = I_ERI_Dyz_S_Dxy_S_C1001001001_a+CDY*I_ERI_Dyz_S_Px_S_C1001001001_a;
  Double I_ERI_D2z_S_Px_Py_C1001001001_a = I_ERI_D2z_S_Dxy_S_C1001001001_a+CDY*I_ERI_D2z_S_Px_S_C1001001001_a;
  Double I_ERI_D2x_S_Py_Py_C1001001001_a = I_ERI_D2x_S_D2y_S_C1001001001_a+CDY*I_ERI_D2x_S_Py_S_C1001001001_a;
  Double I_ERI_Dxy_S_Py_Py_C1001001001_a = I_ERI_Dxy_S_D2y_S_C1001001001_a+CDY*I_ERI_Dxy_S_Py_S_C1001001001_a;
  Double I_ERI_Dxz_S_Py_Py_C1001001001_a = I_ERI_Dxz_S_D2y_S_C1001001001_a+CDY*I_ERI_Dxz_S_Py_S_C1001001001_a;
  Double I_ERI_D2y_S_Py_Py_C1001001001_a = I_ERI_D2y_S_D2y_S_C1001001001_a+CDY*I_ERI_D2y_S_Py_S_C1001001001_a;
  Double I_ERI_Dyz_S_Py_Py_C1001001001_a = I_ERI_Dyz_S_D2y_S_C1001001001_a+CDY*I_ERI_Dyz_S_Py_S_C1001001001_a;
  Double I_ERI_D2z_S_Py_Py_C1001001001_a = I_ERI_D2z_S_D2y_S_C1001001001_a+CDY*I_ERI_D2z_S_Py_S_C1001001001_a;
  Double I_ERI_D2x_S_Pz_Py_C1001001001_a = I_ERI_D2x_S_Dyz_S_C1001001001_a+CDY*I_ERI_D2x_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxy_S_Pz_Py_C1001001001_a = I_ERI_Dxy_S_Dyz_S_C1001001001_a+CDY*I_ERI_Dxy_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxz_S_Pz_Py_C1001001001_a = I_ERI_Dxz_S_Dyz_S_C1001001001_a+CDY*I_ERI_Dxz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2y_S_Pz_Py_C1001001001_a = I_ERI_D2y_S_Dyz_S_C1001001001_a+CDY*I_ERI_D2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Dyz_S_Pz_Py_C1001001001_a = I_ERI_Dyz_S_Dyz_S_C1001001001_a+CDY*I_ERI_Dyz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2z_S_Pz_Py_C1001001001_a = I_ERI_D2z_S_Dyz_S_C1001001001_a+CDY*I_ERI_D2z_S_Pz_S_C1001001001_a;
  Double I_ERI_D2x_S_Px_Pz_C1001001001_a = I_ERI_D2x_S_Dxz_S_C1001001001_a+CDZ*I_ERI_D2x_S_Px_S_C1001001001_a;
  Double I_ERI_Dxy_S_Px_Pz_C1001001001_a = I_ERI_Dxy_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Dxy_S_Px_S_C1001001001_a;
  Double I_ERI_Dxz_S_Px_Pz_C1001001001_a = I_ERI_Dxz_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Dxz_S_Px_S_C1001001001_a;
  Double I_ERI_D2y_S_Px_Pz_C1001001001_a = I_ERI_D2y_S_Dxz_S_C1001001001_a+CDZ*I_ERI_D2y_S_Px_S_C1001001001_a;
  Double I_ERI_Dyz_S_Px_Pz_C1001001001_a = I_ERI_Dyz_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Dyz_S_Px_S_C1001001001_a;
  Double I_ERI_D2z_S_Px_Pz_C1001001001_a = I_ERI_D2z_S_Dxz_S_C1001001001_a+CDZ*I_ERI_D2z_S_Px_S_C1001001001_a;
  Double I_ERI_D2x_S_Py_Pz_C1001001001_a = I_ERI_D2x_S_Dyz_S_C1001001001_a+CDZ*I_ERI_D2x_S_Py_S_C1001001001_a;
  Double I_ERI_Dxy_S_Py_Pz_C1001001001_a = I_ERI_Dxy_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Dxy_S_Py_S_C1001001001_a;
  Double I_ERI_Dxz_S_Py_Pz_C1001001001_a = I_ERI_Dxz_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Dxz_S_Py_S_C1001001001_a;
  Double I_ERI_D2y_S_Py_Pz_C1001001001_a = I_ERI_D2y_S_Dyz_S_C1001001001_a+CDZ*I_ERI_D2y_S_Py_S_C1001001001_a;
  Double I_ERI_Dyz_S_Py_Pz_C1001001001_a = I_ERI_Dyz_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Dyz_S_Py_S_C1001001001_a;
  Double I_ERI_D2z_S_Py_Pz_C1001001001_a = I_ERI_D2z_S_Dyz_S_C1001001001_a+CDZ*I_ERI_D2z_S_Py_S_C1001001001_a;
  Double I_ERI_D2x_S_Pz_Pz_C1001001001_a = I_ERI_D2x_S_D2z_S_C1001001001_a+CDZ*I_ERI_D2x_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001001_a = I_ERI_Dxy_S_D2z_S_C1001001001_a+CDZ*I_ERI_Dxy_S_Pz_S_C1001001001_a;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001001_a = I_ERI_Dxz_S_D2z_S_C1001001001_a+CDZ*I_ERI_Dxz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2y_S_Pz_Pz_C1001001001_a = I_ERI_D2y_S_D2z_S_C1001001001_a+CDZ*I_ERI_D2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001001_a = I_ERI_Dyz_S_D2z_S_C1001001001_a+CDZ*I_ERI_Dyz_S_Pz_S_C1001001001_a;
  Double I_ERI_D2z_S_Pz_Pz_C1001001001_a = I_ERI_D2z_S_D2z_S_C1001001001_a+CDZ*I_ERI_D2z_S_Pz_S_C1001001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001001001_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001001001_a
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001001_a
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1001001001_a = I_ERI_F3x_S_D2x_S_C1001001001_a+CDX*I_ERI_F3x_S_Px_S_C1001001001_a;
  Double I_ERI_F2xy_S_Px_Px_C1001001001_a = I_ERI_F2xy_S_D2x_S_C1001001001_a+CDX*I_ERI_F2xy_S_Px_S_C1001001001_a;
  Double I_ERI_F2xz_S_Px_Px_C1001001001_a = I_ERI_F2xz_S_D2x_S_C1001001001_a+CDX*I_ERI_F2xz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Px_Px_C1001001001_a = I_ERI_Fx2y_S_D2x_S_C1001001001_a+CDX*I_ERI_Fx2y_S_Px_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Px_Px_C1001001001_a = I_ERI_Fxyz_S_D2x_S_C1001001001_a+CDX*I_ERI_Fxyz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Px_Px_C1001001001_a = I_ERI_Fx2z_S_D2x_S_C1001001001_a+CDX*I_ERI_Fx2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3y_S_Px_Px_C1001001001_a = I_ERI_F3y_S_D2x_S_C1001001001_a+CDX*I_ERI_F3y_S_Px_S_C1001001001_a;
  Double I_ERI_F2yz_S_Px_Px_C1001001001_a = I_ERI_F2yz_S_D2x_S_C1001001001_a+CDX*I_ERI_F2yz_S_Px_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Px_Px_C1001001001_a = I_ERI_Fy2z_S_D2x_S_C1001001001_a+CDX*I_ERI_Fy2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3z_S_Px_Px_C1001001001_a = I_ERI_F3z_S_D2x_S_C1001001001_a+CDX*I_ERI_F3z_S_Px_S_C1001001001_a;
  Double I_ERI_F3x_S_Py_Px_C1001001001_a = I_ERI_F3x_S_Dxy_S_C1001001001_a+CDX*I_ERI_F3x_S_Py_S_C1001001001_a;
  Double I_ERI_F2xy_S_Py_Px_C1001001001_a = I_ERI_F2xy_S_Dxy_S_C1001001001_a+CDX*I_ERI_F2xy_S_Py_S_C1001001001_a;
  Double I_ERI_F2xz_S_Py_Px_C1001001001_a = I_ERI_F2xz_S_Dxy_S_C1001001001_a+CDX*I_ERI_F2xz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Py_Px_C1001001001_a = I_ERI_Fx2y_S_Dxy_S_C1001001001_a+CDX*I_ERI_Fx2y_S_Py_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Py_Px_C1001001001_a = I_ERI_Fxyz_S_Dxy_S_C1001001001_a+CDX*I_ERI_Fxyz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Py_Px_C1001001001_a = I_ERI_Fx2z_S_Dxy_S_C1001001001_a+CDX*I_ERI_Fx2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3y_S_Py_Px_C1001001001_a = I_ERI_F3y_S_Dxy_S_C1001001001_a+CDX*I_ERI_F3y_S_Py_S_C1001001001_a;
  Double I_ERI_F2yz_S_Py_Px_C1001001001_a = I_ERI_F2yz_S_Dxy_S_C1001001001_a+CDX*I_ERI_F2yz_S_Py_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Py_Px_C1001001001_a = I_ERI_Fy2z_S_Dxy_S_C1001001001_a+CDX*I_ERI_Fy2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3z_S_Py_Px_C1001001001_a = I_ERI_F3z_S_Dxy_S_C1001001001_a+CDX*I_ERI_F3z_S_Py_S_C1001001001_a;
  Double I_ERI_F3x_S_Pz_Px_C1001001001_a = I_ERI_F3x_S_Dxz_S_C1001001001_a+CDX*I_ERI_F3x_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xy_S_Pz_Px_C1001001001_a = I_ERI_F2xy_S_Dxz_S_C1001001001_a+CDX*I_ERI_F2xy_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xz_S_Pz_Px_C1001001001_a = I_ERI_F2xz_S_Dxz_S_C1001001001_a+CDX*I_ERI_F2xz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Pz_Px_C1001001001_a = I_ERI_Fx2y_S_Dxz_S_C1001001001_a+CDX*I_ERI_Fx2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Pz_Px_C1001001001_a = I_ERI_Fxyz_S_Dxz_S_C1001001001_a+CDX*I_ERI_Fxyz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Pz_Px_C1001001001_a = I_ERI_Fx2z_S_Dxz_S_C1001001001_a+CDX*I_ERI_Fx2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3y_S_Pz_Px_C1001001001_a = I_ERI_F3y_S_Dxz_S_C1001001001_a+CDX*I_ERI_F3y_S_Pz_S_C1001001001_a;
  Double I_ERI_F2yz_S_Pz_Px_C1001001001_a = I_ERI_F2yz_S_Dxz_S_C1001001001_a+CDX*I_ERI_F2yz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Pz_Px_C1001001001_a = I_ERI_Fy2z_S_Dxz_S_C1001001001_a+CDX*I_ERI_Fy2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3z_S_Pz_Px_C1001001001_a = I_ERI_F3z_S_Dxz_S_C1001001001_a+CDX*I_ERI_F3z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3x_S_Px_Py_C1001001001_a = I_ERI_F3x_S_Dxy_S_C1001001001_a+CDY*I_ERI_F3x_S_Px_S_C1001001001_a;
  Double I_ERI_F2xy_S_Px_Py_C1001001001_a = I_ERI_F2xy_S_Dxy_S_C1001001001_a+CDY*I_ERI_F2xy_S_Px_S_C1001001001_a;
  Double I_ERI_F2xz_S_Px_Py_C1001001001_a = I_ERI_F2xz_S_Dxy_S_C1001001001_a+CDY*I_ERI_F2xz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Px_Py_C1001001001_a = I_ERI_Fx2y_S_Dxy_S_C1001001001_a+CDY*I_ERI_Fx2y_S_Px_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Px_Py_C1001001001_a = I_ERI_Fxyz_S_Dxy_S_C1001001001_a+CDY*I_ERI_Fxyz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Px_Py_C1001001001_a = I_ERI_Fx2z_S_Dxy_S_C1001001001_a+CDY*I_ERI_Fx2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3y_S_Px_Py_C1001001001_a = I_ERI_F3y_S_Dxy_S_C1001001001_a+CDY*I_ERI_F3y_S_Px_S_C1001001001_a;
  Double I_ERI_F2yz_S_Px_Py_C1001001001_a = I_ERI_F2yz_S_Dxy_S_C1001001001_a+CDY*I_ERI_F2yz_S_Px_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Px_Py_C1001001001_a = I_ERI_Fy2z_S_Dxy_S_C1001001001_a+CDY*I_ERI_Fy2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3z_S_Px_Py_C1001001001_a = I_ERI_F3z_S_Dxy_S_C1001001001_a+CDY*I_ERI_F3z_S_Px_S_C1001001001_a;
  Double I_ERI_F3x_S_Py_Py_C1001001001_a = I_ERI_F3x_S_D2y_S_C1001001001_a+CDY*I_ERI_F3x_S_Py_S_C1001001001_a;
  Double I_ERI_F2xy_S_Py_Py_C1001001001_a = I_ERI_F2xy_S_D2y_S_C1001001001_a+CDY*I_ERI_F2xy_S_Py_S_C1001001001_a;
  Double I_ERI_F2xz_S_Py_Py_C1001001001_a = I_ERI_F2xz_S_D2y_S_C1001001001_a+CDY*I_ERI_F2xz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Py_Py_C1001001001_a = I_ERI_Fx2y_S_D2y_S_C1001001001_a+CDY*I_ERI_Fx2y_S_Py_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Py_Py_C1001001001_a = I_ERI_Fxyz_S_D2y_S_C1001001001_a+CDY*I_ERI_Fxyz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Py_Py_C1001001001_a = I_ERI_Fx2z_S_D2y_S_C1001001001_a+CDY*I_ERI_Fx2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3y_S_Py_Py_C1001001001_a = I_ERI_F3y_S_D2y_S_C1001001001_a+CDY*I_ERI_F3y_S_Py_S_C1001001001_a;
  Double I_ERI_F2yz_S_Py_Py_C1001001001_a = I_ERI_F2yz_S_D2y_S_C1001001001_a+CDY*I_ERI_F2yz_S_Py_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Py_Py_C1001001001_a = I_ERI_Fy2z_S_D2y_S_C1001001001_a+CDY*I_ERI_Fy2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3z_S_Py_Py_C1001001001_a = I_ERI_F3z_S_D2y_S_C1001001001_a+CDY*I_ERI_F3z_S_Py_S_C1001001001_a;
  Double I_ERI_F3x_S_Pz_Py_C1001001001_a = I_ERI_F3x_S_Dyz_S_C1001001001_a+CDY*I_ERI_F3x_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xy_S_Pz_Py_C1001001001_a = I_ERI_F2xy_S_Dyz_S_C1001001001_a+CDY*I_ERI_F2xy_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xz_S_Pz_Py_C1001001001_a = I_ERI_F2xz_S_Dyz_S_C1001001001_a+CDY*I_ERI_F2xz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Pz_Py_C1001001001_a = I_ERI_Fx2y_S_Dyz_S_C1001001001_a+CDY*I_ERI_Fx2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Pz_Py_C1001001001_a = I_ERI_Fxyz_S_Dyz_S_C1001001001_a+CDY*I_ERI_Fxyz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Pz_Py_C1001001001_a = I_ERI_Fx2z_S_Dyz_S_C1001001001_a+CDY*I_ERI_Fx2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3y_S_Pz_Py_C1001001001_a = I_ERI_F3y_S_Dyz_S_C1001001001_a+CDY*I_ERI_F3y_S_Pz_S_C1001001001_a;
  Double I_ERI_F2yz_S_Pz_Py_C1001001001_a = I_ERI_F2yz_S_Dyz_S_C1001001001_a+CDY*I_ERI_F2yz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Pz_Py_C1001001001_a = I_ERI_Fy2z_S_Dyz_S_C1001001001_a+CDY*I_ERI_Fy2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3z_S_Pz_Py_C1001001001_a = I_ERI_F3z_S_Dyz_S_C1001001001_a+CDY*I_ERI_F3z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3x_S_Px_Pz_C1001001001_a = I_ERI_F3x_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F3x_S_Px_S_C1001001001_a;
  Double I_ERI_F2xy_S_Px_Pz_C1001001001_a = I_ERI_F2xy_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F2xy_S_Px_S_C1001001001_a;
  Double I_ERI_F2xz_S_Px_Pz_C1001001001_a = I_ERI_F2xz_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F2xz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Px_Pz_C1001001001_a = I_ERI_Fx2y_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Fx2y_S_Px_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Px_Pz_C1001001001_a = I_ERI_Fxyz_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Fxyz_S_Px_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Px_Pz_C1001001001_a = I_ERI_Fx2z_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Fx2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3y_S_Px_Pz_C1001001001_a = I_ERI_F3y_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F3y_S_Px_S_C1001001001_a;
  Double I_ERI_F2yz_S_Px_Pz_C1001001001_a = I_ERI_F2yz_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F2yz_S_Px_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Px_Pz_C1001001001_a = I_ERI_Fy2z_S_Dxz_S_C1001001001_a+CDZ*I_ERI_Fy2z_S_Px_S_C1001001001_a;
  Double I_ERI_F3z_S_Px_Pz_C1001001001_a = I_ERI_F3z_S_Dxz_S_C1001001001_a+CDZ*I_ERI_F3z_S_Px_S_C1001001001_a;
  Double I_ERI_F3x_S_Py_Pz_C1001001001_a = I_ERI_F3x_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F3x_S_Py_S_C1001001001_a;
  Double I_ERI_F2xy_S_Py_Pz_C1001001001_a = I_ERI_F2xy_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F2xy_S_Py_S_C1001001001_a;
  Double I_ERI_F2xz_S_Py_Pz_C1001001001_a = I_ERI_F2xz_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F2xz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Py_Pz_C1001001001_a = I_ERI_Fx2y_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Fx2y_S_Py_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Py_Pz_C1001001001_a = I_ERI_Fxyz_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Fxyz_S_Py_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Py_Pz_C1001001001_a = I_ERI_Fx2z_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Fx2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3y_S_Py_Pz_C1001001001_a = I_ERI_F3y_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F3y_S_Py_S_C1001001001_a;
  Double I_ERI_F2yz_S_Py_Pz_C1001001001_a = I_ERI_F2yz_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F2yz_S_Py_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Py_Pz_C1001001001_a = I_ERI_Fy2z_S_Dyz_S_C1001001001_a+CDZ*I_ERI_Fy2z_S_Py_S_C1001001001_a;
  Double I_ERI_F3z_S_Py_Pz_C1001001001_a = I_ERI_F3z_S_Dyz_S_C1001001001_a+CDZ*I_ERI_F3z_S_Py_S_C1001001001_a;
  Double I_ERI_F3x_S_Pz_Pz_C1001001001_a = I_ERI_F3x_S_D2z_S_C1001001001_a+CDZ*I_ERI_F3x_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xy_S_Pz_Pz_C1001001001_a = I_ERI_F2xy_S_D2z_S_C1001001001_a+CDZ*I_ERI_F2xy_S_Pz_S_C1001001001_a;
  Double I_ERI_F2xz_S_Pz_Pz_C1001001001_a = I_ERI_F2xz_S_D2z_S_C1001001001_a+CDZ*I_ERI_F2xz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2y_S_Pz_Pz_C1001001001_a = I_ERI_Fx2y_S_D2z_S_C1001001001_a+CDZ*I_ERI_Fx2y_S_Pz_S_C1001001001_a;
  Double I_ERI_Fxyz_S_Pz_Pz_C1001001001_a = I_ERI_Fxyz_S_D2z_S_C1001001001_a+CDZ*I_ERI_Fxyz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fx2z_S_Pz_Pz_C1001001001_a = I_ERI_Fx2z_S_D2z_S_C1001001001_a+CDZ*I_ERI_Fx2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3y_S_Pz_Pz_C1001001001_a = I_ERI_F3y_S_D2z_S_C1001001001_a+CDZ*I_ERI_F3y_S_Pz_S_C1001001001_a;
  Double I_ERI_F2yz_S_Pz_Pz_C1001001001_a = I_ERI_F2yz_S_D2z_S_C1001001001_a+CDZ*I_ERI_F2yz_S_Pz_S_C1001001001_a;
  Double I_ERI_Fy2z_S_Pz_Pz_C1001001001_a = I_ERI_Fy2z_S_D2z_S_C1001001001_a+CDZ*I_ERI_Fy2z_S_Pz_S_C1001001001_a;
  Double I_ERI_F3z_S_Pz_Pz_C1001001001_a = I_ERI_F3z_S_D2z_S_C1001001001_a+CDZ*I_ERI_F3z_S_Pz_S_C1001001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001000000_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000000_b
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001000000_b
   ************************************************************/
  Double I_ERI_S_Px_Px_Px_C1001000000_b = I_ERI_S_Px_D2x_S_C1001000000_b+CDX*I_ERI_S_Px_Px_S_C1001000000_b;
  Double I_ERI_S_Py_Px_Px_C1001000000_b = I_ERI_S_Py_D2x_S_C1001000000_b+CDX*I_ERI_S_Py_Px_S_C1001000000_b;
  Double I_ERI_S_Pz_Px_Px_C1001000000_b = I_ERI_S_Pz_D2x_S_C1001000000_b+CDX*I_ERI_S_Pz_Px_S_C1001000000_b;
  Double I_ERI_S_Px_Py_Px_C1001000000_b = I_ERI_S_Px_Dxy_S_C1001000000_b+CDX*I_ERI_S_Px_Py_S_C1001000000_b;
  Double I_ERI_S_Py_Py_Px_C1001000000_b = I_ERI_S_Py_Dxy_S_C1001000000_b+CDX*I_ERI_S_Py_Py_S_C1001000000_b;
  Double I_ERI_S_Pz_Py_Px_C1001000000_b = I_ERI_S_Pz_Dxy_S_C1001000000_b+CDX*I_ERI_S_Pz_Py_S_C1001000000_b;
  Double I_ERI_S_Px_Pz_Px_C1001000000_b = I_ERI_S_Px_Dxz_S_C1001000000_b+CDX*I_ERI_S_Px_Pz_S_C1001000000_b;
  Double I_ERI_S_Py_Pz_Px_C1001000000_b = I_ERI_S_Py_Dxz_S_C1001000000_b+CDX*I_ERI_S_Py_Pz_S_C1001000000_b;
  Double I_ERI_S_Pz_Pz_Px_C1001000000_b = I_ERI_S_Pz_Dxz_S_C1001000000_b+CDX*I_ERI_S_Pz_Pz_S_C1001000000_b;
  Double I_ERI_S_Px_Px_Py_C1001000000_b = I_ERI_S_Px_Dxy_S_C1001000000_b+CDY*I_ERI_S_Px_Px_S_C1001000000_b;
  Double I_ERI_S_Py_Px_Py_C1001000000_b = I_ERI_S_Py_Dxy_S_C1001000000_b+CDY*I_ERI_S_Py_Px_S_C1001000000_b;
  Double I_ERI_S_Pz_Px_Py_C1001000000_b = I_ERI_S_Pz_Dxy_S_C1001000000_b+CDY*I_ERI_S_Pz_Px_S_C1001000000_b;
  Double I_ERI_S_Px_Py_Py_C1001000000_b = I_ERI_S_Px_D2y_S_C1001000000_b+CDY*I_ERI_S_Px_Py_S_C1001000000_b;
  Double I_ERI_S_Py_Py_Py_C1001000000_b = I_ERI_S_Py_D2y_S_C1001000000_b+CDY*I_ERI_S_Py_Py_S_C1001000000_b;
  Double I_ERI_S_Pz_Py_Py_C1001000000_b = I_ERI_S_Pz_D2y_S_C1001000000_b+CDY*I_ERI_S_Pz_Py_S_C1001000000_b;
  Double I_ERI_S_Px_Pz_Py_C1001000000_b = I_ERI_S_Px_Dyz_S_C1001000000_b+CDY*I_ERI_S_Px_Pz_S_C1001000000_b;
  Double I_ERI_S_Py_Pz_Py_C1001000000_b = I_ERI_S_Py_Dyz_S_C1001000000_b+CDY*I_ERI_S_Py_Pz_S_C1001000000_b;
  Double I_ERI_S_Pz_Pz_Py_C1001000000_b = I_ERI_S_Pz_Dyz_S_C1001000000_b+CDY*I_ERI_S_Pz_Pz_S_C1001000000_b;
  Double I_ERI_S_Px_Px_Pz_C1001000000_b = I_ERI_S_Px_Dxz_S_C1001000000_b+CDZ*I_ERI_S_Px_Px_S_C1001000000_b;
  Double I_ERI_S_Py_Px_Pz_C1001000000_b = I_ERI_S_Py_Dxz_S_C1001000000_b+CDZ*I_ERI_S_Py_Px_S_C1001000000_b;
  Double I_ERI_S_Pz_Px_Pz_C1001000000_b = I_ERI_S_Pz_Dxz_S_C1001000000_b+CDZ*I_ERI_S_Pz_Px_S_C1001000000_b;
  Double I_ERI_S_Px_Py_Pz_C1001000000_b = I_ERI_S_Px_Dyz_S_C1001000000_b+CDZ*I_ERI_S_Px_Py_S_C1001000000_b;
  Double I_ERI_S_Py_Py_Pz_C1001000000_b = I_ERI_S_Py_Dyz_S_C1001000000_b+CDZ*I_ERI_S_Py_Py_S_C1001000000_b;
  Double I_ERI_S_Pz_Py_Pz_C1001000000_b = I_ERI_S_Pz_Dyz_S_C1001000000_b+CDZ*I_ERI_S_Pz_Py_S_C1001000000_b;
  Double I_ERI_S_Px_Pz_Pz_C1001000000_b = I_ERI_S_Px_D2z_S_C1001000000_b+CDZ*I_ERI_S_Px_Pz_S_C1001000000_b;
  Double I_ERI_S_Py_Pz_Pz_C1001000000_b = I_ERI_S_Py_D2z_S_C1001000000_b+CDZ*I_ERI_S_Py_Pz_S_C1001000000_b;
  Double I_ERI_S_Pz_Pz_Pz_C1001000000_b = I_ERI_S_Pz_D2z_S_C1001000000_b+CDZ*I_ERI_S_Pz_Pz_S_C1001000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001000001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001000001_b
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001000001_b = I_ERI_Px_S_D2x_S_C1001000001_b+CDX*I_ERI_Px_S_Px_S_C1001000001_b;
  Double I_ERI_Py_S_Px_Px_C1001000001_b = I_ERI_Py_S_D2x_S_C1001000001_b+CDX*I_ERI_Py_S_Px_S_C1001000001_b;
  Double I_ERI_Pz_S_Px_Px_C1001000001_b = I_ERI_Pz_S_D2x_S_C1001000001_b+CDX*I_ERI_Pz_S_Px_S_C1001000001_b;
  Double I_ERI_Px_S_Py_Px_C1001000001_b = I_ERI_Px_S_Dxy_S_C1001000001_b+CDX*I_ERI_Px_S_Py_S_C1001000001_b;
  Double I_ERI_Py_S_Py_Px_C1001000001_b = I_ERI_Py_S_Dxy_S_C1001000001_b+CDX*I_ERI_Py_S_Py_S_C1001000001_b;
  Double I_ERI_Pz_S_Py_Px_C1001000001_b = I_ERI_Pz_S_Dxy_S_C1001000001_b+CDX*I_ERI_Pz_S_Py_S_C1001000001_b;
  Double I_ERI_Px_S_Pz_Px_C1001000001_b = I_ERI_Px_S_Dxz_S_C1001000001_b+CDX*I_ERI_Px_S_Pz_S_C1001000001_b;
  Double I_ERI_Py_S_Pz_Px_C1001000001_b = I_ERI_Py_S_Dxz_S_C1001000001_b+CDX*I_ERI_Py_S_Pz_S_C1001000001_b;
  Double I_ERI_Pz_S_Pz_Px_C1001000001_b = I_ERI_Pz_S_Dxz_S_C1001000001_b+CDX*I_ERI_Pz_S_Pz_S_C1001000001_b;
  Double I_ERI_Px_S_Px_Py_C1001000001_b = I_ERI_Px_S_Dxy_S_C1001000001_b+CDY*I_ERI_Px_S_Px_S_C1001000001_b;
  Double I_ERI_Py_S_Px_Py_C1001000001_b = I_ERI_Py_S_Dxy_S_C1001000001_b+CDY*I_ERI_Py_S_Px_S_C1001000001_b;
  Double I_ERI_Pz_S_Px_Py_C1001000001_b = I_ERI_Pz_S_Dxy_S_C1001000001_b+CDY*I_ERI_Pz_S_Px_S_C1001000001_b;
  Double I_ERI_Px_S_Py_Py_C1001000001_b = I_ERI_Px_S_D2y_S_C1001000001_b+CDY*I_ERI_Px_S_Py_S_C1001000001_b;
  Double I_ERI_Py_S_Py_Py_C1001000001_b = I_ERI_Py_S_D2y_S_C1001000001_b+CDY*I_ERI_Py_S_Py_S_C1001000001_b;
  Double I_ERI_Pz_S_Py_Py_C1001000001_b = I_ERI_Pz_S_D2y_S_C1001000001_b+CDY*I_ERI_Pz_S_Py_S_C1001000001_b;
  Double I_ERI_Px_S_Pz_Py_C1001000001_b = I_ERI_Px_S_Dyz_S_C1001000001_b+CDY*I_ERI_Px_S_Pz_S_C1001000001_b;
  Double I_ERI_Py_S_Pz_Py_C1001000001_b = I_ERI_Py_S_Dyz_S_C1001000001_b+CDY*I_ERI_Py_S_Pz_S_C1001000001_b;
  Double I_ERI_Pz_S_Pz_Py_C1001000001_b = I_ERI_Pz_S_Dyz_S_C1001000001_b+CDY*I_ERI_Pz_S_Pz_S_C1001000001_b;
  Double I_ERI_Px_S_Px_Pz_C1001000001_b = I_ERI_Px_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Px_S_Px_S_C1001000001_b;
  Double I_ERI_Py_S_Px_Pz_C1001000001_b = I_ERI_Py_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Py_S_Px_S_C1001000001_b;
  Double I_ERI_Pz_S_Px_Pz_C1001000001_b = I_ERI_Pz_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Pz_S_Px_S_C1001000001_b;
  Double I_ERI_Px_S_Py_Pz_C1001000001_b = I_ERI_Px_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Px_S_Py_S_C1001000001_b;
  Double I_ERI_Py_S_Py_Pz_C1001000001_b = I_ERI_Py_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Py_S_Py_S_C1001000001_b;
  Double I_ERI_Pz_S_Py_Pz_C1001000001_b = I_ERI_Pz_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Pz_S_Py_S_C1001000001_b;
  Double I_ERI_Px_S_Pz_Pz_C1001000001_b = I_ERI_Px_S_D2z_S_C1001000001_b+CDZ*I_ERI_Px_S_Pz_S_C1001000001_b;
  Double I_ERI_Py_S_Pz_Pz_C1001000001_b = I_ERI_Py_S_D2z_S_C1001000001_b+CDZ*I_ERI_Py_S_Pz_S_C1001000001_b;
  Double I_ERI_Pz_S_Pz_Pz_C1001000001_b = I_ERI_Pz_S_D2z_S_C1001000001_b+CDZ*I_ERI_Pz_S_Pz_S_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001000001_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001000001_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001000001_b
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001000001_b = I_ERI_D2x_S_D2x_S_C1001000001_b+CDX*I_ERI_D2x_S_Px_S_C1001000001_b;
  Double I_ERI_Dxy_S_Px_Px_C1001000001_b = I_ERI_Dxy_S_D2x_S_C1001000001_b+CDX*I_ERI_Dxy_S_Px_S_C1001000001_b;
  Double I_ERI_Dxz_S_Px_Px_C1001000001_b = I_ERI_Dxz_S_D2x_S_C1001000001_b+CDX*I_ERI_Dxz_S_Px_S_C1001000001_b;
  Double I_ERI_D2y_S_Px_Px_C1001000001_b = I_ERI_D2y_S_D2x_S_C1001000001_b+CDX*I_ERI_D2y_S_Px_S_C1001000001_b;
  Double I_ERI_Dyz_S_Px_Px_C1001000001_b = I_ERI_Dyz_S_D2x_S_C1001000001_b+CDX*I_ERI_Dyz_S_Px_S_C1001000001_b;
  Double I_ERI_D2z_S_Px_Px_C1001000001_b = I_ERI_D2z_S_D2x_S_C1001000001_b+CDX*I_ERI_D2z_S_Px_S_C1001000001_b;
  Double I_ERI_D2x_S_Py_Px_C1001000001_b = I_ERI_D2x_S_Dxy_S_C1001000001_b+CDX*I_ERI_D2x_S_Py_S_C1001000001_b;
  Double I_ERI_Dxy_S_Py_Px_C1001000001_b = I_ERI_Dxy_S_Dxy_S_C1001000001_b+CDX*I_ERI_Dxy_S_Py_S_C1001000001_b;
  Double I_ERI_Dxz_S_Py_Px_C1001000001_b = I_ERI_Dxz_S_Dxy_S_C1001000001_b+CDX*I_ERI_Dxz_S_Py_S_C1001000001_b;
  Double I_ERI_D2y_S_Py_Px_C1001000001_b = I_ERI_D2y_S_Dxy_S_C1001000001_b+CDX*I_ERI_D2y_S_Py_S_C1001000001_b;
  Double I_ERI_Dyz_S_Py_Px_C1001000001_b = I_ERI_Dyz_S_Dxy_S_C1001000001_b+CDX*I_ERI_Dyz_S_Py_S_C1001000001_b;
  Double I_ERI_D2z_S_Py_Px_C1001000001_b = I_ERI_D2z_S_Dxy_S_C1001000001_b+CDX*I_ERI_D2z_S_Py_S_C1001000001_b;
  Double I_ERI_D2x_S_Pz_Px_C1001000001_b = I_ERI_D2x_S_Dxz_S_C1001000001_b+CDX*I_ERI_D2x_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxy_S_Pz_Px_C1001000001_b = I_ERI_Dxy_S_Dxz_S_C1001000001_b+CDX*I_ERI_Dxy_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxz_S_Pz_Px_C1001000001_b = I_ERI_Dxz_S_Dxz_S_C1001000001_b+CDX*I_ERI_Dxz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2y_S_Pz_Px_C1001000001_b = I_ERI_D2y_S_Dxz_S_C1001000001_b+CDX*I_ERI_D2y_S_Pz_S_C1001000001_b;
  Double I_ERI_Dyz_S_Pz_Px_C1001000001_b = I_ERI_Dyz_S_Dxz_S_C1001000001_b+CDX*I_ERI_Dyz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2z_S_Pz_Px_C1001000001_b = I_ERI_D2z_S_Dxz_S_C1001000001_b+CDX*I_ERI_D2z_S_Pz_S_C1001000001_b;
  Double I_ERI_D2x_S_Px_Py_C1001000001_b = I_ERI_D2x_S_Dxy_S_C1001000001_b+CDY*I_ERI_D2x_S_Px_S_C1001000001_b;
  Double I_ERI_Dxy_S_Px_Py_C1001000001_b = I_ERI_Dxy_S_Dxy_S_C1001000001_b+CDY*I_ERI_Dxy_S_Px_S_C1001000001_b;
  Double I_ERI_Dxz_S_Px_Py_C1001000001_b = I_ERI_Dxz_S_Dxy_S_C1001000001_b+CDY*I_ERI_Dxz_S_Px_S_C1001000001_b;
  Double I_ERI_D2y_S_Px_Py_C1001000001_b = I_ERI_D2y_S_Dxy_S_C1001000001_b+CDY*I_ERI_D2y_S_Px_S_C1001000001_b;
  Double I_ERI_Dyz_S_Px_Py_C1001000001_b = I_ERI_Dyz_S_Dxy_S_C1001000001_b+CDY*I_ERI_Dyz_S_Px_S_C1001000001_b;
  Double I_ERI_D2z_S_Px_Py_C1001000001_b = I_ERI_D2z_S_Dxy_S_C1001000001_b+CDY*I_ERI_D2z_S_Px_S_C1001000001_b;
  Double I_ERI_D2x_S_Py_Py_C1001000001_b = I_ERI_D2x_S_D2y_S_C1001000001_b+CDY*I_ERI_D2x_S_Py_S_C1001000001_b;
  Double I_ERI_Dxy_S_Py_Py_C1001000001_b = I_ERI_Dxy_S_D2y_S_C1001000001_b+CDY*I_ERI_Dxy_S_Py_S_C1001000001_b;
  Double I_ERI_Dxz_S_Py_Py_C1001000001_b = I_ERI_Dxz_S_D2y_S_C1001000001_b+CDY*I_ERI_Dxz_S_Py_S_C1001000001_b;
  Double I_ERI_D2y_S_Py_Py_C1001000001_b = I_ERI_D2y_S_D2y_S_C1001000001_b+CDY*I_ERI_D2y_S_Py_S_C1001000001_b;
  Double I_ERI_Dyz_S_Py_Py_C1001000001_b = I_ERI_Dyz_S_D2y_S_C1001000001_b+CDY*I_ERI_Dyz_S_Py_S_C1001000001_b;
  Double I_ERI_D2z_S_Py_Py_C1001000001_b = I_ERI_D2z_S_D2y_S_C1001000001_b+CDY*I_ERI_D2z_S_Py_S_C1001000001_b;
  Double I_ERI_D2x_S_Pz_Py_C1001000001_b = I_ERI_D2x_S_Dyz_S_C1001000001_b+CDY*I_ERI_D2x_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxy_S_Pz_Py_C1001000001_b = I_ERI_Dxy_S_Dyz_S_C1001000001_b+CDY*I_ERI_Dxy_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxz_S_Pz_Py_C1001000001_b = I_ERI_Dxz_S_Dyz_S_C1001000001_b+CDY*I_ERI_Dxz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2y_S_Pz_Py_C1001000001_b = I_ERI_D2y_S_Dyz_S_C1001000001_b+CDY*I_ERI_D2y_S_Pz_S_C1001000001_b;
  Double I_ERI_Dyz_S_Pz_Py_C1001000001_b = I_ERI_Dyz_S_Dyz_S_C1001000001_b+CDY*I_ERI_Dyz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2z_S_Pz_Py_C1001000001_b = I_ERI_D2z_S_Dyz_S_C1001000001_b+CDY*I_ERI_D2z_S_Pz_S_C1001000001_b;
  Double I_ERI_D2x_S_Px_Pz_C1001000001_b = I_ERI_D2x_S_Dxz_S_C1001000001_b+CDZ*I_ERI_D2x_S_Px_S_C1001000001_b;
  Double I_ERI_Dxy_S_Px_Pz_C1001000001_b = I_ERI_Dxy_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Dxy_S_Px_S_C1001000001_b;
  Double I_ERI_Dxz_S_Px_Pz_C1001000001_b = I_ERI_Dxz_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Dxz_S_Px_S_C1001000001_b;
  Double I_ERI_D2y_S_Px_Pz_C1001000001_b = I_ERI_D2y_S_Dxz_S_C1001000001_b+CDZ*I_ERI_D2y_S_Px_S_C1001000001_b;
  Double I_ERI_Dyz_S_Px_Pz_C1001000001_b = I_ERI_Dyz_S_Dxz_S_C1001000001_b+CDZ*I_ERI_Dyz_S_Px_S_C1001000001_b;
  Double I_ERI_D2z_S_Px_Pz_C1001000001_b = I_ERI_D2z_S_Dxz_S_C1001000001_b+CDZ*I_ERI_D2z_S_Px_S_C1001000001_b;
  Double I_ERI_D2x_S_Py_Pz_C1001000001_b = I_ERI_D2x_S_Dyz_S_C1001000001_b+CDZ*I_ERI_D2x_S_Py_S_C1001000001_b;
  Double I_ERI_Dxy_S_Py_Pz_C1001000001_b = I_ERI_Dxy_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Dxy_S_Py_S_C1001000001_b;
  Double I_ERI_Dxz_S_Py_Pz_C1001000001_b = I_ERI_Dxz_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Dxz_S_Py_S_C1001000001_b;
  Double I_ERI_D2y_S_Py_Pz_C1001000001_b = I_ERI_D2y_S_Dyz_S_C1001000001_b+CDZ*I_ERI_D2y_S_Py_S_C1001000001_b;
  Double I_ERI_Dyz_S_Py_Pz_C1001000001_b = I_ERI_Dyz_S_Dyz_S_C1001000001_b+CDZ*I_ERI_Dyz_S_Py_S_C1001000001_b;
  Double I_ERI_D2z_S_Py_Pz_C1001000001_b = I_ERI_D2z_S_Dyz_S_C1001000001_b+CDZ*I_ERI_D2z_S_Py_S_C1001000001_b;
  Double I_ERI_D2x_S_Pz_Pz_C1001000001_b = I_ERI_D2x_S_D2z_S_C1001000001_b+CDZ*I_ERI_D2x_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxy_S_Pz_Pz_C1001000001_b = I_ERI_Dxy_S_D2z_S_C1001000001_b+CDZ*I_ERI_Dxy_S_Pz_S_C1001000001_b;
  Double I_ERI_Dxz_S_Pz_Pz_C1001000001_b = I_ERI_Dxz_S_D2z_S_C1001000001_b+CDZ*I_ERI_Dxz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2y_S_Pz_Pz_C1001000001_b = I_ERI_D2y_S_D2z_S_C1001000001_b+CDZ*I_ERI_D2y_S_Pz_S_C1001000001_b;
  Double I_ERI_Dyz_S_Pz_Pz_C1001000001_b = I_ERI_Dyz_S_D2z_S_C1001000001_b+CDZ*I_ERI_Dyz_S_Pz_S_C1001000001_b;
  Double I_ERI_D2z_S_Pz_Pz_C1001000001_b = I_ERI_D2z_S_D2z_S_C1001000001_b+CDZ*I_ERI_D2z_S_Pz_S_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_D_P_P_C1001001000_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_D_S_C1001001000_b
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1001001000_b
   ************************************************************/
  Double I_ERI_S_D2x_Px_Px_C1001001000_b = I_ERI_S_D2x_D2x_S_C1001001000_b+CDX*I_ERI_S_D2x_Px_S_C1001001000_b;
  Double I_ERI_S_Dxy_Px_Px_C1001001000_b = I_ERI_S_Dxy_D2x_S_C1001001000_b+CDX*I_ERI_S_Dxy_Px_S_C1001001000_b;
  Double I_ERI_S_Dxz_Px_Px_C1001001000_b = I_ERI_S_Dxz_D2x_S_C1001001000_b+CDX*I_ERI_S_Dxz_Px_S_C1001001000_b;
  Double I_ERI_S_D2y_Px_Px_C1001001000_b = I_ERI_S_D2y_D2x_S_C1001001000_b+CDX*I_ERI_S_D2y_Px_S_C1001001000_b;
  Double I_ERI_S_Dyz_Px_Px_C1001001000_b = I_ERI_S_Dyz_D2x_S_C1001001000_b+CDX*I_ERI_S_Dyz_Px_S_C1001001000_b;
  Double I_ERI_S_D2z_Px_Px_C1001001000_b = I_ERI_S_D2z_D2x_S_C1001001000_b+CDX*I_ERI_S_D2z_Px_S_C1001001000_b;
  Double I_ERI_S_D2x_Py_Px_C1001001000_b = I_ERI_S_D2x_Dxy_S_C1001001000_b+CDX*I_ERI_S_D2x_Py_S_C1001001000_b;
  Double I_ERI_S_Dxy_Py_Px_C1001001000_b = I_ERI_S_Dxy_Dxy_S_C1001001000_b+CDX*I_ERI_S_Dxy_Py_S_C1001001000_b;
  Double I_ERI_S_Dxz_Py_Px_C1001001000_b = I_ERI_S_Dxz_Dxy_S_C1001001000_b+CDX*I_ERI_S_Dxz_Py_S_C1001001000_b;
  Double I_ERI_S_D2y_Py_Px_C1001001000_b = I_ERI_S_D2y_Dxy_S_C1001001000_b+CDX*I_ERI_S_D2y_Py_S_C1001001000_b;
  Double I_ERI_S_Dyz_Py_Px_C1001001000_b = I_ERI_S_Dyz_Dxy_S_C1001001000_b+CDX*I_ERI_S_Dyz_Py_S_C1001001000_b;
  Double I_ERI_S_D2z_Py_Px_C1001001000_b = I_ERI_S_D2z_Dxy_S_C1001001000_b+CDX*I_ERI_S_D2z_Py_S_C1001001000_b;
  Double I_ERI_S_D2x_Pz_Px_C1001001000_b = I_ERI_S_D2x_Dxz_S_C1001001000_b+CDX*I_ERI_S_D2x_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxy_Pz_Px_C1001001000_b = I_ERI_S_Dxy_Dxz_S_C1001001000_b+CDX*I_ERI_S_Dxy_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxz_Pz_Px_C1001001000_b = I_ERI_S_Dxz_Dxz_S_C1001001000_b+CDX*I_ERI_S_Dxz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2y_Pz_Px_C1001001000_b = I_ERI_S_D2y_Dxz_S_C1001001000_b+CDX*I_ERI_S_D2y_Pz_S_C1001001000_b;
  Double I_ERI_S_Dyz_Pz_Px_C1001001000_b = I_ERI_S_Dyz_Dxz_S_C1001001000_b+CDX*I_ERI_S_Dyz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2z_Pz_Px_C1001001000_b = I_ERI_S_D2z_Dxz_S_C1001001000_b+CDX*I_ERI_S_D2z_Pz_S_C1001001000_b;
  Double I_ERI_S_D2x_Px_Py_C1001001000_b = I_ERI_S_D2x_Dxy_S_C1001001000_b+CDY*I_ERI_S_D2x_Px_S_C1001001000_b;
  Double I_ERI_S_Dxy_Px_Py_C1001001000_b = I_ERI_S_Dxy_Dxy_S_C1001001000_b+CDY*I_ERI_S_Dxy_Px_S_C1001001000_b;
  Double I_ERI_S_Dxz_Px_Py_C1001001000_b = I_ERI_S_Dxz_Dxy_S_C1001001000_b+CDY*I_ERI_S_Dxz_Px_S_C1001001000_b;
  Double I_ERI_S_D2y_Px_Py_C1001001000_b = I_ERI_S_D2y_Dxy_S_C1001001000_b+CDY*I_ERI_S_D2y_Px_S_C1001001000_b;
  Double I_ERI_S_Dyz_Px_Py_C1001001000_b = I_ERI_S_Dyz_Dxy_S_C1001001000_b+CDY*I_ERI_S_Dyz_Px_S_C1001001000_b;
  Double I_ERI_S_D2z_Px_Py_C1001001000_b = I_ERI_S_D2z_Dxy_S_C1001001000_b+CDY*I_ERI_S_D2z_Px_S_C1001001000_b;
  Double I_ERI_S_D2x_Py_Py_C1001001000_b = I_ERI_S_D2x_D2y_S_C1001001000_b+CDY*I_ERI_S_D2x_Py_S_C1001001000_b;
  Double I_ERI_S_Dxy_Py_Py_C1001001000_b = I_ERI_S_Dxy_D2y_S_C1001001000_b+CDY*I_ERI_S_Dxy_Py_S_C1001001000_b;
  Double I_ERI_S_Dxz_Py_Py_C1001001000_b = I_ERI_S_Dxz_D2y_S_C1001001000_b+CDY*I_ERI_S_Dxz_Py_S_C1001001000_b;
  Double I_ERI_S_D2y_Py_Py_C1001001000_b = I_ERI_S_D2y_D2y_S_C1001001000_b+CDY*I_ERI_S_D2y_Py_S_C1001001000_b;
  Double I_ERI_S_Dyz_Py_Py_C1001001000_b = I_ERI_S_Dyz_D2y_S_C1001001000_b+CDY*I_ERI_S_Dyz_Py_S_C1001001000_b;
  Double I_ERI_S_D2z_Py_Py_C1001001000_b = I_ERI_S_D2z_D2y_S_C1001001000_b+CDY*I_ERI_S_D2z_Py_S_C1001001000_b;
  Double I_ERI_S_D2x_Pz_Py_C1001001000_b = I_ERI_S_D2x_Dyz_S_C1001001000_b+CDY*I_ERI_S_D2x_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxy_Pz_Py_C1001001000_b = I_ERI_S_Dxy_Dyz_S_C1001001000_b+CDY*I_ERI_S_Dxy_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxz_Pz_Py_C1001001000_b = I_ERI_S_Dxz_Dyz_S_C1001001000_b+CDY*I_ERI_S_Dxz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2y_Pz_Py_C1001001000_b = I_ERI_S_D2y_Dyz_S_C1001001000_b+CDY*I_ERI_S_D2y_Pz_S_C1001001000_b;
  Double I_ERI_S_Dyz_Pz_Py_C1001001000_b = I_ERI_S_Dyz_Dyz_S_C1001001000_b+CDY*I_ERI_S_Dyz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2z_Pz_Py_C1001001000_b = I_ERI_S_D2z_Dyz_S_C1001001000_b+CDY*I_ERI_S_D2z_Pz_S_C1001001000_b;
  Double I_ERI_S_D2x_Px_Pz_C1001001000_b = I_ERI_S_D2x_Dxz_S_C1001001000_b+CDZ*I_ERI_S_D2x_Px_S_C1001001000_b;
  Double I_ERI_S_Dxy_Px_Pz_C1001001000_b = I_ERI_S_Dxy_Dxz_S_C1001001000_b+CDZ*I_ERI_S_Dxy_Px_S_C1001001000_b;
  Double I_ERI_S_Dxz_Px_Pz_C1001001000_b = I_ERI_S_Dxz_Dxz_S_C1001001000_b+CDZ*I_ERI_S_Dxz_Px_S_C1001001000_b;
  Double I_ERI_S_D2y_Px_Pz_C1001001000_b = I_ERI_S_D2y_Dxz_S_C1001001000_b+CDZ*I_ERI_S_D2y_Px_S_C1001001000_b;
  Double I_ERI_S_Dyz_Px_Pz_C1001001000_b = I_ERI_S_Dyz_Dxz_S_C1001001000_b+CDZ*I_ERI_S_Dyz_Px_S_C1001001000_b;
  Double I_ERI_S_D2z_Px_Pz_C1001001000_b = I_ERI_S_D2z_Dxz_S_C1001001000_b+CDZ*I_ERI_S_D2z_Px_S_C1001001000_b;
  Double I_ERI_S_D2x_Py_Pz_C1001001000_b = I_ERI_S_D2x_Dyz_S_C1001001000_b+CDZ*I_ERI_S_D2x_Py_S_C1001001000_b;
  Double I_ERI_S_Dxy_Py_Pz_C1001001000_b = I_ERI_S_Dxy_Dyz_S_C1001001000_b+CDZ*I_ERI_S_Dxy_Py_S_C1001001000_b;
  Double I_ERI_S_Dxz_Py_Pz_C1001001000_b = I_ERI_S_Dxz_Dyz_S_C1001001000_b+CDZ*I_ERI_S_Dxz_Py_S_C1001001000_b;
  Double I_ERI_S_D2y_Py_Pz_C1001001000_b = I_ERI_S_D2y_Dyz_S_C1001001000_b+CDZ*I_ERI_S_D2y_Py_S_C1001001000_b;
  Double I_ERI_S_Dyz_Py_Pz_C1001001000_b = I_ERI_S_Dyz_Dyz_S_C1001001000_b+CDZ*I_ERI_S_Dyz_Py_S_C1001001000_b;
  Double I_ERI_S_D2z_Py_Pz_C1001001000_b = I_ERI_S_D2z_Dyz_S_C1001001000_b+CDZ*I_ERI_S_D2z_Py_S_C1001001000_b;
  Double I_ERI_S_D2x_Pz_Pz_C1001001000_b = I_ERI_S_D2x_D2z_S_C1001001000_b+CDZ*I_ERI_S_D2x_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxy_Pz_Pz_C1001001000_b = I_ERI_S_Dxy_D2z_S_C1001001000_b+CDZ*I_ERI_S_Dxy_Pz_S_C1001001000_b;
  Double I_ERI_S_Dxz_Pz_Pz_C1001001000_b = I_ERI_S_Dxz_D2z_S_C1001001000_b+CDZ*I_ERI_S_Dxz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2y_Pz_Pz_C1001001000_b = I_ERI_S_D2y_D2z_S_C1001001000_b+CDZ*I_ERI_S_D2y_Pz_S_C1001001000_b;
  Double I_ERI_S_Dyz_Pz_Pz_C1001001000_b = I_ERI_S_Dyz_D2z_S_C1001001000_b+CDZ*I_ERI_S_Dyz_Pz_S_C1001001000_b;
  Double I_ERI_S_D2z_Pz_Pz_C1001001000_b = I_ERI_S_D2z_D2z_S_C1001001000_b+CDZ*I_ERI_S_D2z_Pz_S_C1001001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001001001_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001001_b
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001001001_b = I_ERI_Px_S_D2x_S_C1001001001_b+CDX*I_ERI_Px_S_Px_S_C1001001001_b;
  Double I_ERI_Py_S_Px_Px_C1001001001_b = I_ERI_Py_S_D2x_S_C1001001001_b+CDX*I_ERI_Py_S_Px_S_C1001001001_b;
  Double I_ERI_Pz_S_Px_Px_C1001001001_b = I_ERI_Pz_S_D2x_S_C1001001001_b+CDX*I_ERI_Pz_S_Px_S_C1001001001_b;
  Double I_ERI_Px_S_Py_Px_C1001001001_b = I_ERI_Px_S_Dxy_S_C1001001001_b+CDX*I_ERI_Px_S_Py_S_C1001001001_b;
  Double I_ERI_Py_S_Py_Px_C1001001001_b = I_ERI_Py_S_Dxy_S_C1001001001_b+CDX*I_ERI_Py_S_Py_S_C1001001001_b;
  Double I_ERI_Pz_S_Py_Px_C1001001001_b = I_ERI_Pz_S_Dxy_S_C1001001001_b+CDX*I_ERI_Pz_S_Py_S_C1001001001_b;
  Double I_ERI_Px_S_Pz_Px_C1001001001_b = I_ERI_Px_S_Dxz_S_C1001001001_b+CDX*I_ERI_Px_S_Pz_S_C1001001001_b;
  Double I_ERI_Py_S_Pz_Px_C1001001001_b = I_ERI_Py_S_Dxz_S_C1001001001_b+CDX*I_ERI_Py_S_Pz_S_C1001001001_b;
  Double I_ERI_Pz_S_Pz_Px_C1001001001_b = I_ERI_Pz_S_Dxz_S_C1001001001_b+CDX*I_ERI_Pz_S_Pz_S_C1001001001_b;
  Double I_ERI_Px_S_Px_Py_C1001001001_b = I_ERI_Px_S_Dxy_S_C1001001001_b+CDY*I_ERI_Px_S_Px_S_C1001001001_b;
  Double I_ERI_Py_S_Px_Py_C1001001001_b = I_ERI_Py_S_Dxy_S_C1001001001_b+CDY*I_ERI_Py_S_Px_S_C1001001001_b;
  Double I_ERI_Pz_S_Px_Py_C1001001001_b = I_ERI_Pz_S_Dxy_S_C1001001001_b+CDY*I_ERI_Pz_S_Px_S_C1001001001_b;
  Double I_ERI_Px_S_Py_Py_C1001001001_b = I_ERI_Px_S_D2y_S_C1001001001_b+CDY*I_ERI_Px_S_Py_S_C1001001001_b;
  Double I_ERI_Py_S_Py_Py_C1001001001_b = I_ERI_Py_S_D2y_S_C1001001001_b+CDY*I_ERI_Py_S_Py_S_C1001001001_b;
  Double I_ERI_Pz_S_Py_Py_C1001001001_b = I_ERI_Pz_S_D2y_S_C1001001001_b+CDY*I_ERI_Pz_S_Py_S_C1001001001_b;
  Double I_ERI_Px_S_Pz_Py_C1001001001_b = I_ERI_Px_S_Dyz_S_C1001001001_b+CDY*I_ERI_Px_S_Pz_S_C1001001001_b;
  Double I_ERI_Py_S_Pz_Py_C1001001001_b = I_ERI_Py_S_Dyz_S_C1001001001_b+CDY*I_ERI_Py_S_Pz_S_C1001001001_b;
  Double I_ERI_Pz_S_Pz_Py_C1001001001_b = I_ERI_Pz_S_Dyz_S_C1001001001_b+CDY*I_ERI_Pz_S_Pz_S_C1001001001_b;
  Double I_ERI_Px_S_Px_Pz_C1001001001_b = I_ERI_Px_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Px_S_Px_S_C1001001001_b;
  Double I_ERI_Py_S_Px_Pz_C1001001001_b = I_ERI_Py_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Py_S_Px_S_C1001001001_b;
  Double I_ERI_Pz_S_Px_Pz_C1001001001_b = I_ERI_Pz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Pz_S_Px_S_C1001001001_b;
  Double I_ERI_Px_S_Py_Pz_C1001001001_b = I_ERI_Px_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Px_S_Py_S_C1001001001_b;
  Double I_ERI_Py_S_Py_Pz_C1001001001_b = I_ERI_Py_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Py_S_Py_S_C1001001001_b;
  Double I_ERI_Pz_S_Py_Pz_C1001001001_b = I_ERI_Pz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Pz_S_Py_S_C1001001001_b;
  Double I_ERI_Px_S_Pz_Pz_C1001001001_b = I_ERI_Px_S_D2z_S_C1001001001_b+CDZ*I_ERI_Px_S_Pz_S_C1001001001_b;
  Double I_ERI_Py_S_Pz_Pz_C1001001001_b = I_ERI_Py_S_D2z_S_C1001001001_b+CDZ*I_ERI_Py_S_Pz_S_C1001001001_b;
  Double I_ERI_Pz_S_Pz_Pz_C1001001001_b = I_ERI_Pz_S_D2z_S_C1001001001_b+CDZ*I_ERI_Pz_S_Pz_S_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001001_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001001_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001001_b
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001001_b = I_ERI_D2x_S_D2x_S_C1001001001_b+CDX*I_ERI_D2x_S_Px_S_C1001001001_b;
  Double I_ERI_Dxy_S_Px_Px_C1001001001_b = I_ERI_Dxy_S_D2x_S_C1001001001_b+CDX*I_ERI_Dxy_S_Px_S_C1001001001_b;
  Double I_ERI_Dxz_S_Px_Px_C1001001001_b = I_ERI_Dxz_S_D2x_S_C1001001001_b+CDX*I_ERI_Dxz_S_Px_S_C1001001001_b;
  Double I_ERI_D2y_S_Px_Px_C1001001001_b = I_ERI_D2y_S_D2x_S_C1001001001_b+CDX*I_ERI_D2y_S_Px_S_C1001001001_b;
  Double I_ERI_Dyz_S_Px_Px_C1001001001_b = I_ERI_Dyz_S_D2x_S_C1001001001_b+CDX*I_ERI_Dyz_S_Px_S_C1001001001_b;
  Double I_ERI_D2z_S_Px_Px_C1001001001_b = I_ERI_D2z_S_D2x_S_C1001001001_b+CDX*I_ERI_D2z_S_Px_S_C1001001001_b;
  Double I_ERI_D2x_S_Py_Px_C1001001001_b = I_ERI_D2x_S_Dxy_S_C1001001001_b+CDX*I_ERI_D2x_S_Py_S_C1001001001_b;
  Double I_ERI_Dxy_S_Py_Px_C1001001001_b = I_ERI_Dxy_S_Dxy_S_C1001001001_b+CDX*I_ERI_Dxy_S_Py_S_C1001001001_b;
  Double I_ERI_Dxz_S_Py_Px_C1001001001_b = I_ERI_Dxz_S_Dxy_S_C1001001001_b+CDX*I_ERI_Dxz_S_Py_S_C1001001001_b;
  Double I_ERI_D2y_S_Py_Px_C1001001001_b = I_ERI_D2y_S_Dxy_S_C1001001001_b+CDX*I_ERI_D2y_S_Py_S_C1001001001_b;
  Double I_ERI_Dyz_S_Py_Px_C1001001001_b = I_ERI_Dyz_S_Dxy_S_C1001001001_b+CDX*I_ERI_Dyz_S_Py_S_C1001001001_b;
  Double I_ERI_D2z_S_Py_Px_C1001001001_b = I_ERI_D2z_S_Dxy_S_C1001001001_b+CDX*I_ERI_D2z_S_Py_S_C1001001001_b;
  Double I_ERI_D2x_S_Pz_Px_C1001001001_b = I_ERI_D2x_S_Dxz_S_C1001001001_b+CDX*I_ERI_D2x_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxy_S_Pz_Px_C1001001001_b = I_ERI_Dxy_S_Dxz_S_C1001001001_b+CDX*I_ERI_Dxy_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxz_S_Pz_Px_C1001001001_b = I_ERI_Dxz_S_Dxz_S_C1001001001_b+CDX*I_ERI_Dxz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2y_S_Pz_Px_C1001001001_b = I_ERI_D2y_S_Dxz_S_C1001001001_b+CDX*I_ERI_D2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Dyz_S_Pz_Px_C1001001001_b = I_ERI_Dyz_S_Dxz_S_C1001001001_b+CDX*I_ERI_Dyz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2z_S_Pz_Px_C1001001001_b = I_ERI_D2z_S_Dxz_S_C1001001001_b+CDX*I_ERI_D2z_S_Pz_S_C1001001001_b;
  Double I_ERI_D2x_S_Px_Py_C1001001001_b = I_ERI_D2x_S_Dxy_S_C1001001001_b+CDY*I_ERI_D2x_S_Px_S_C1001001001_b;
  Double I_ERI_Dxy_S_Px_Py_C1001001001_b = I_ERI_Dxy_S_Dxy_S_C1001001001_b+CDY*I_ERI_Dxy_S_Px_S_C1001001001_b;
  Double I_ERI_Dxz_S_Px_Py_C1001001001_b = I_ERI_Dxz_S_Dxy_S_C1001001001_b+CDY*I_ERI_Dxz_S_Px_S_C1001001001_b;
  Double I_ERI_D2y_S_Px_Py_C1001001001_b = I_ERI_D2y_S_Dxy_S_C1001001001_b+CDY*I_ERI_D2y_S_Px_S_C1001001001_b;
  Double I_ERI_Dyz_S_Px_Py_C1001001001_b = I_ERI_Dyz_S_Dxy_S_C1001001001_b+CDY*I_ERI_Dyz_S_Px_S_C1001001001_b;
  Double I_ERI_D2z_S_Px_Py_C1001001001_b = I_ERI_D2z_S_Dxy_S_C1001001001_b+CDY*I_ERI_D2z_S_Px_S_C1001001001_b;
  Double I_ERI_D2x_S_Py_Py_C1001001001_b = I_ERI_D2x_S_D2y_S_C1001001001_b+CDY*I_ERI_D2x_S_Py_S_C1001001001_b;
  Double I_ERI_Dxy_S_Py_Py_C1001001001_b = I_ERI_Dxy_S_D2y_S_C1001001001_b+CDY*I_ERI_Dxy_S_Py_S_C1001001001_b;
  Double I_ERI_Dxz_S_Py_Py_C1001001001_b = I_ERI_Dxz_S_D2y_S_C1001001001_b+CDY*I_ERI_Dxz_S_Py_S_C1001001001_b;
  Double I_ERI_D2y_S_Py_Py_C1001001001_b = I_ERI_D2y_S_D2y_S_C1001001001_b+CDY*I_ERI_D2y_S_Py_S_C1001001001_b;
  Double I_ERI_Dyz_S_Py_Py_C1001001001_b = I_ERI_Dyz_S_D2y_S_C1001001001_b+CDY*I_ERI_Dyz_S_Py_S_C1001001001_b;
  Double I_ERI_D2z_S_Py_Py_C1001001001_b = I_ERI_D2z_S_D2y_S_C1001001001_b+CDY*I_ERI_D2z_S_Py_S_C1001001001_b;
  Double I_ERI_D2x_S_Pz_Py_C1001001001_b = I_ERI_D2x_S_Dyz_S_C1001001001_b+CDY*I_ERI_D2x_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxy_S_Pz_Py_C1001001001_b = I_ERI_Dxy_S_Dyz_S_C1001001001_b+CDY*I_ERI_Dxy_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxz_S_Pz_Py_C1001001001_b = I_ERI_Dxz_S_Dyz_S_C1001001001_b+CDY*I_ERI_Dxz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2y_S_Pz_Py_C1001001001_b = I_ERI_D2y_S_Dyz_S_C1001001001_b+CDY*I_ERI_D2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Dyz_S_Pz_Py_C1001001001_b = I_ERI_Dyz_S_Dyz_S_C1001001001_b+CDY*I_ERI_Dyz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2z_S_Pz_Py_C1001001001_b = I_ERI_D2z_S_Dyz_S_C1001001001_b+CDY*I_ERI_D2z_S_Pz_S_C1001001001_b;
  Double I_ERI_D2x_S_Px_Pz_C1001001001_b = I_ERI_D2x_S_Dxz_S_C1001001001_b+CDZ*I_ERI_D2x_S_Px_S_C1001001001_b;
  Double I_ERI_Dxy_S_Px_Pz_C1001001001_b = I_ERI_Dxy_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Dxy_S_Px_S_C1001001001_b;
  Double I_ERI_Dxz_S_Px_Pz_C1001001001_b = I_ERI_Dxz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Dxz_S_Px_S_C1001001001_b;
  Double I_ERI_D2y_S_Px_Pz_C1001001001_b = I_ERI_D2y_S_Dxz_S_C1001001001_b+CDZ*I_ERI_D2y_S_Px_S_C1001001001_b;
  Double I_ERI_Dyz_S_Px_Pz_C1001001001_b = I_ERI_Dyz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Dyz_S_Px_S_C1001001001_b;
  Double I_ERI_D2z_S_Px_Pz_C1001001001_b = I_ERI_D2z_S_Dxz_S_C1001001001_b+CDZ*I_ERI_D2z_S_Px_S_C1001001001_b;
  Double I_ERI_D2x_S_Py_Pz_C1001001001_b = I_ERI_D2x_S_Dyz_S_C1001001001_b+CDZ*I_ERI_D2x_S_Py_S_C1001001001_b;
  Double I_ERI_Dxy_S_Py_Pz_C1001001001_b = I_ERI_Dxy_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Dxy_S_Py_S_C1001001001_b;
  Double I_ERI_Dxz_S_Py_Pz_C1001001001_b = I_ERI_Dxz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Dxz_S_Py_S_C1001001001_b;
  Double I_ERI_D2y_S_Py_Pz_C1001001001_b = I_ERI_D2y_S_Dyz_S_C1001001001_b+CDZ*I_ERI_D2y_S_Py_S_C1001001001_b;
  Double I_ERI_Dyz_S_Py_Pz_C1001001001_b = I_ERI_Dyz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Dyz_S_Py_S_C1001001001_b;
  Double I_ERI_D2z_S_Py_Pz_C1001001001_b = I_ERI_D2z_S_Dyz_S_C1001001001_b+CDZ*I_ERI_D2z_S_Py_S_C1001001001_b;
  Double I_ERI_D2x_S_Pz_Pz_C1001001001_b = I_ERI_D2x_S_D2z_S_C1001001001_b+CDZ*I_ERI_D2x_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001001_b = I_ERI_Dxy_S_D2z_S_C1001001001_b+CDZ*I_ERI_Dxy_S_Pz_S_C1001001001_b;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001001_b = I_ERI_Dxz_S_D2z_S_C1001001001_b+CDZ*I_ERI_Dxz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2y_S_Pz_Pz_C1001001001_b = I_ERI_D2y_S_D2z_S_C1001001001_b+CDZ*I_ERI_D2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001001_b = I_ERI_Dyz_S_D2z_S_C1001001001_b+CDZ*I_ERI_Dyz_S_Pz_S_C1001001001_b;
  Double I_ERI_D2z_S_Pz_Pz_C1001001001_b = I_ERI_D2z_S_D2z_S_C1001001001_b+CDZ*I_ERI_D2z_S_Pz_S_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001001001_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001001001_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001001001_b
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1001001001_b = I_ERI_F3x_S_D2x_S_C1001001001_b+CDX*I_ERI_F3x_S_Px_S_C1001001001_b;
  Double I_ERI_F2xy_S_Px_Px_C1001001001_b = I_ERI_F2xy_S_D2x_S_C1001001001_b+CDX*I_ERI_F2xy_S_Px_S_C1001001001_b;
  Double I_ERI_F2xz_S_Px_Px_C1001001001_b = I_ERI_F2xz_S_D2x_S_C1001001001_b+CDX*I_ERI_F2xz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Px_Px_C1001001001_b = I_ERI_Fx2y_S_D2x_S_C1001001001_b+CDX*I_ERI_Fx2y_S_Px_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Px_Px_C1001001001_b = I_ERI_Fxyz_S_D2x_S_C1001001001_b+CDX*I_ERI_Fxyz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Px_Px_C1001001001_b = I_ERI_Fx2z_S_D2x_S_C1001001001_b+CDX*I_ERI_Fx2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3y_S_Px_Px_C1001001001_b = I_ERI_F3y_S_D2x_S_C1001001001_b+CDX*I_ERI_F3y_S_Px_S_C1001001001_b;
  Double I_ERI_F2yz_S_Px_Px_C1001001001_b = I_ERI_F2yz_S_D2x_S_C1001001001_b+CDX*I_ERI_F2yz_S_Px_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Px_Px_C1001001001_b = I_ERI_Fy2z_S_D2x_S_C1001001001_b+CDX*I_ERI_Fy2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3z_S_Px_Px_C1001001001_b = I_ERI_F3z_S_D2x_S_C1001001001_b+CDX*I_ERI_F3z_S_Px_S_C1001001001_b;
  Double I_ERI_F3x_S_Py_Px_C1001001001_b = I_ERI_F3x_S_Dxy_S_C1001001001_b+CDX*I_ERI_F3x_S_Py_S_C1001001001_b;
  Double I_ERI_F2xy_S_Py_Px_C1001001001_b = I_ERI_F2xy_S_Dxy_S_C1001001001_b+CDX*I_ERI_F2xy_S_Py_S_C1001001001_b;
  Double I_ERI_F2xz_S_Py_Px_C1001001001_b = I_ERI_F2xz_S_Dxy_S_C1001001001_b+CDX*I_ERI_F2xz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Py_Px_C1001001001_b = I_ERI_Fx2y_S_Dxy_S_C1001001001_b+CDX*I_ERI_Fx2y_S_Py_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Py_Px_C1001001001_b = I_ERI_Fxyz_S_Dxy_S_C1001001001_b+CDX*I_ERI_Fxyz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Py_Px_C1001001001_b = I_ERI_Fx2z_S_Dxy_S_C1001001001_b+CDX*I_ERI_Fx2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3y_S_Py_Px_C1001001001_b = I_ERI_F3y_S_Dxy_S_C1001001001_b+CDX*I_ERI_F3y_S_Py_S_C1001001001_b;
  Double I_ERI_F2yz_S_Py_Px_C1001001001_b = I_ERI_F2yz_S_Dxy_S_C1001001001_b+CDX*I_ERI_F2yz_S_Py_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Py_Px_C1001001001_b = I_ERI_Fy2z_S_Dxy_S_C1001001001_b+CDX*I_ERI_Fy2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3z_S_Py_Px_C1001001001_b = I_ERI_F3z_S_Dxy_S_C1001001001_b+CDX*I_ERI_F3z_S_Py_S_C1001001001_b;
  Double I_ERI_F3x_S_Pz_Px_C1001001001_b = I_ERI_F3x_S_Dxz_S_C1001001001_b+CDX*I_ERI_F3x_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xy_S_Pz_Px_C1001001001_b = I_ERI_F2xy_S_Dxz_S_C1001001001_b+CDX*I_ERI_F2xy_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xz_S_Pz_Px_C1001001001_b = I_ERI_F2xz_S_Dxz_S_C1001001001_b+CDX*I_ERI_F2xz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Pz_Px_C1001001001_b = I_ERI_Fx2y_S_Dxz_S_C1001001001_b+CDX*I_ERI_Fx2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Pz_Px_C1001001001_b = I_ERI_Fxyz_S_Dxz_S_C1001001001_b+CDX*I_ERI_Fxyz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Pz_Px_C1001001001_b = I_ERI_Fx2z_S_Dxz_S_C1001001001_b+CDX*I_ERI_Fx2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3y_S_Pz_Px_C1001001001_b = I_ERI_F3y_S_Dxz_S_C1001001001_b+CDX*I_ERI_F3y_S_Pz_S_C1001001001_b;
  Double I_ERI_F2yz_S_Pz_Px_C1001001001_b = I_ERI_F2yz_S_Dxz_S_C1001001001_b+CDX*I_ERI_F2yz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Pz_Px_C1001001001_b = I_ERI_Fy2z_S_Dxz_S_C1001001001_b+CDX*I_ERI_Fy2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3z_S_Pz_Px_C1001001001_b = I_ERI_F3z_S_Dxz_S_C1001001001_b+CDX*I_ERI_F3z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3x_S_Px_Py_C1001001001_b = I_ERI_F3x_S_Dxy_S_C1001001001_b+CDY*I_ERI_F3x_S_Px_S_C1001001001_b;
  Double I_ERI_F2xy_S_Px_Py_C1001001001_b = I_ERI_F2xy_S_Dxy_S_C1001001001_b+CDY*I_ERI_F2xy_S_Px_S_C1001001001_b;
  Double I_ERI_F2xz_S_Px_Py_C1001001001_b = I_ERI_F2xz_S_Dxy_S_C1001001001_b+CDY*I_ERI_F2xz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Px_Py_C1001001001_b = I_ERI_Fx2y_S_Dxy_S_C1001001001_b+CDY*I_ERI_Fx2y_S_Px_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Px_Py_C1001001001_b = I_ERI_Fxyz_S_Dxy_S_C1001001001_b+CDY*I_ERI_Fxyz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Px_Py_C1001001001_b = I_ERI_Fx2z_S_Dxy_S_C1001001001_b+CDY*I_ERI_Fx2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3y_S_Px_Py_C1001001001_b = I_ERI_F3y_S_Dxy_S_C1001001001_b+CDY*I_ERI_F3y_S_Px_S_C1001001001_b;
  Double I_ERI_F2yz_S_Px_Py_C1001001001_b = I_ERI_F2yz_S_Dxy_S_C1001001001_b+CDY*I_ERI_F2yz_S_Px_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Px_Py_C1001001001_b = I_ERI_Fy2z_S_Dxy_S_C1001001001_b+CDY*I_ERI_Fy2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3z_S_Px_Py_C1001001001_b = I_ERI_F3z_S_Dxy_S_C1001001001_b+CDY*I_ERI_F3z_S_Px_S_C1001001001_b;
  Double I_ERI_F3x_S_Py_Py_C1001001001_b = I_ERI_F3x_S_D2y_S_C1001001001_b+CDY*I_ERI_F3x_S_Py_S_C1001001001_b;
  Double I_ERI_F2xy_S_Py_Py_C1001001001_b = I_ERI_F2xy_S_D2y_S_C1001001001_b+CDY*I_ERI_F2xy_S_Py_S_C1001001001_b;
  Double I_ERI_F2xz_S_Py_Py_C1001001001_b = I_ERI_F2xz_S_D2y_S_C1001001001_b+CDY*I_ERI_F2xz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Py_Py_C1001001001_b = I_ERI_Fx2y_S_D2y_S_C1001001001_b+CDY*I_ERI_Fx2y_S_Py_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Py_Py_C1001001001_b = I_ERI_Fxyz_S_D2y_S_C1001001001_b+CDY*I_ERI_Fxyz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Py_Py_C1001001001_b = I_ERI_Fx2z_S_D2y_S_C1001001001_b+CDY*I_ERI_Fx2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3y_S_Py_Py_C1001001001_b = I_ERI_F3y_S_D2y_S_C1001001001_b+CDY*I_ERI_F3y_S_Py_S_C1001001001_b;
  Double I_ERI_F2yz_S_Py_Py_C1001001001_b = I_ERI_F2yz_S_D2y_S_C1001001001_b+CDY*I_ERI_F2yz_S_Py_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Py_Py_C1001001001_b = I_ERI_Fy2z_S_D2y_S_C1001001001_b+CDY*I_ERI_Fy2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3z_S_Py_Py_C1001001001_b = I_ERI_F3z_S_D2y_S_C1001001001_b+CDY*I_ERI_F3z_S_Py_S_C1001001001_b;
  Double I_ERI_F3x_S_Pz_Py_C1001001001_b = I_ERI_F3x_S_Dyz_S_C1001001001_b+CDY*I_ERI_F3x_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xy_S_Pz_Py_C1001001001_b = I_ERI_F2xy_S_Dyz_S_C1001001001_b+CDY*I_ERI_F2xy_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xz_S_Pz_Py_C1001001001_b = I_ERI_F2xz_S_Dyz_S_C1001001001_b+CDY*I_ERI_F2xz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Pz_Py_C1001001001_b = I_ERI_Fx2y_S_Dyz_S_C1001001001_b+CDY*I_ERI_Fx2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Pz_Py_C1001001001_b = I_ERI_Fxyz_S_Dyz_S_C1001001001_b+CDY*I_ERI_Fxyz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Pz_Py_C1001001001_b = I_ERI_Fx2z_S_Dyz_S_C1001001001_b+CDY*I_ERI_Fx2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3y_S_Pz_Py_C1001001001_b = I_ERI_F3y_S_Dyz_S_C1001001001_b+CDY*I_ERI_F3y_S_Pz_S_C1001001001_b;
  Double I_ERI_F2yz_S_Pz_Py_C1001001001_b = I_ERI_F2yz_S_Dyz_S_C1001001001_b+CDY*I_ERI_F2yz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Pz_Py_C1001001001_b = I_ERI_Fy2z_S_Dyz_S_C1001001001_b+CDY*I_ERI_Fy2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3z_S_Pz_Py_C1001001001_b = I_ERI_F3z_S_Dyz_S_C1001001001_b+CDY*I_ERI_F3z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3x_S_Px_Pz_C1001001001_b = I_ERI_F3x_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F3x_S_Px_S_C1001001001_b;
  Double I_ERI_F2xy_S_Px_Pz_C1001001001_b = I_ERI_F2xy_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F2xy_S_Px_S_C1001001001_b;
  Double I_ERI_F2xz_S_Px_Pz_C1001001001_b = I_ERI_F2xz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F2xz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Px_Pz_C1001001001_b = I_ERI_Fx2y_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Fx2y_S_Px_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Px_Pz_C1001001001_b = I_ERI_Fxyz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Fxyz_S_Px_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Px_Pz_C1001001001_b = I_ERI_Fx2z_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Fx2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3y_S_Px_Pz_C1001001001_b = I_ERI_F3y_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F3y_S_Px_S_C1001001001_b;
  Double I_ERI_F2yz_S_Px_Pz_C1001001001_b = I_ERI_F2yz_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F2yz_S_Px_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Px_Pz_C1001001001_b = I_ERI_Fy2z_S_Dxz_S_C1001001001_b+CDZ*I_ERI_Fy2z_S_Px_S_C1001001001_b;
  Double I_ERI_F3z_S_Px_Pz_C1001001001_b = I_ERI_F3z_S_Dxz_S_C1001001001_b+CDZ*I_ERI_F3z_S_Px_S_C1001001001_b;
  Double I_ERI_F3x_S_Py_Pz_C1001001001_b = I_ERI_F3x_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F3x_S_Py_S_C1001001001_b;
  Double I_ERI_F2xy_S_Py_Pz_C1001001001_b = I_ERI_F2xy_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F2xy_S_Py_S_C1001001001_b;
  Double I_ERI_F2xz_S_Py_Pz_C1001001001_b = I_ERI_F2xz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F2xz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Py_Pz_C1001001001_b = I_ERI_Fx2y_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Fx2y_S_Py_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Py_Pz_C1001001001_b = I_ERI_Fxyz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Fxyz_S_Py_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Py_Pz_C1001001001_b = I_ERI_Fx2z_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Fx2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3y_S_Py_Pz_C1001001001_b = I_ERI_F3y_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F3y_S_Py_S_C1001001001_b;
  Double I_ERI_F2yz_S_Py_Pz_C1001001001_b = I_ERI_F2yz_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F2yz_S_Py_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Py_Pz_C1001001001_b = I_ERI_Fy2z_S_Dyz_S_C1001001001_b+CDZ*I_ERI_Fy2z_S_Py_S_C1001001001_b;
  Double I_ERI_F3z_S_Py_Pz_C1001001001_b = I_ERI_F3z_S_Dyz_S_C1001001001_b+CDZ*I_ERI_F3z_S_Py_S_C1001001001_b;
  Double I_ERI_F3x_S_Pz_Pz_C1001001001_b = I_ERI_F3x_S_D2z_S_C1001001001_b+CDZ*I_ERI_F3x_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xy_S_Pz_Pz_C1001001001_b = I_ERI_F2xy_S_D2z_S_C1001001001_b+CDZ*I_ERI_F2xy_S_Pz_S_C1001001001_b;
  Double I_ERI_F2xz_S_Pz_Pz_C1001001001_b = I_ERI_F2xz_S_D2z_S_C1001001001_b+CDZ*I_ERI_F2xz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2y_S_Pz_Pz_C1001001001_b = I_ERI_Fx2y_S_D2z_S_C1001001001_b+CDZ*I_ERI_Fx2y_S_Pz_S_C1001001001_b;
  Double I_ERI_Fxyz_S_Pz_Pz_C1001001001_b = I_ERI_Fxyz_S_D2z_S_C1001001001_b+CDZ*I_ERI_Fxyz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fx2z_S_Pz_Pz_C1001001001_b = I_ERI_Fx2z_S_D2z_S_C1001001001_b+CDZ*I_ERI_Fx2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3y_S_Pz_Pz_C1001001001_b = I_ERI_F3y_S_D2z_S_C1001001001_b+CDZ*I_ERI_F3y_S_Pz_S_C1001001001_b;
  Double I_ERI_F2yz_S_Pz_Pz_C1001001001_b = I_ERI_F2yz_S_D2z_S_C1001001001_b+CDZ*I_ERI_F2yz_S_Pz_S_C1001001001_b;
  Double I_ERI_Fy2z_S_Pz_Pz_C1001001001_b = I_ERI_Fy2z_S_D2z_S_C1001001001_b+CDZ*I_ERI_Fy2z_S_Pz_S_C1001001001_b;
  Double I_ERI_F3z_S_Pz_Pz_C1001001001_b = I_ERI_F3z_S_D2z_S_C1001001001_b+CDZ*I_ERI_F3z_S_Pz_S_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1000000000_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000000_c
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000000000_c
   ************************************************************/
  Double I_ERI_S_S_Px_Px_C1000000000_c = I_ERI_S_S_D2x_S_C1000000000_c+CDX*I_ERI_S_S_Px_S_C1000000000_c;
  Double I_ERI_S_S_Py_Px_C1000000000_c = I_ERI_S_S_Dxy_S_C1000000000_c+CDX*I_ERI_S_S_Py_S_C1000000000_c;
  Double I_ERI_S_S_Pz_Px_C1000000000_c = I_ERI_S_S_Dxz_S_C1000000000_c+CDX*I_ERI_S_S_Pz_S_C1000000000_c;
  Double I_ERI_S_S_Px_Py_C1000000000_c = I_ERI_S_S_Dxy_S_C1000000000_c+CDY*I_ERI_S_S_Px_S_C1000000000_c;
  Double I_ERI_S_S_Py_Py_C1000000000_c = I_ERI_S_S_D2y_S_C1000000000_c+CDY*I_ERI_S_S_Py_S_C1000000000_c;
  Double I_ERI_S_S_Pz_Py_C1000000000_c = I_ERI_S_S_Dyz_S_C1000000000_c+CDY*I_ERI_S_S_Pz_S_C1000000000_c;
  Double I_ERI_S_S_Px_Pz_C1000000000_c = I_ERI_S_S_Dxz_S_C1000000000_c+CDZ*I_ERI_S_S_Px_S_C1000000000_c;
  Double I_ERI_S_S_Py_Pz_C1000000000_c = I_ERI_S_S_Dyz_S_C1000000000_c+CDZ*I_ERI_S_S_Py_S_C1000000000_c;
  Double I_ERI_S_S_Pz_Pz_C1000000000_c = I_ERI_S_S_D2z_S_C1000000000_c+CDZ*I_ERI_S_S_Pz_S_C1000000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1000000001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000000001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000001_c
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1000000001_c = I_ERI_Px_S_D2x_S_C1000000001_c+CDX*I_ERI_Px_S_Px_S_C1000000001_c;
  Double I_ERI_Py_S_Px_Px_C1000000001_c = I_ERI_Py_S_D2x_S_C1000000001_c+CDX*I_ERI_Py_S_Px_S_C1000000001_c;
  Double I_ERI_Pz_S_Px_Px_C1000000001_c = I_ERI_Pz_S_D2x_S_C1000000001_c+CDX*I_ERI_Pz_S_Px_S_C1000000001_c;
  Double I_ERI_Px_S_Py_Px_C1000000001_c = I_ERI_Px_S_Dxy_S_C1000000001_c+CDX*I_ERI_Px_S_Py_S_C1000000001_c;
  Double I_ERI_Py_S_Py_Px_C1000000001_c = I_ERI_Py_S_Dxy_S_C1000000001_c+CDX*I_ERI_Py_S_Py_S_C1000000001_c;
  Double I_ERI_Pz_S_Py_Px_C1000000001_c = I_ERI_Pz_S_Dxy_S_C1000000001_c+CDX*I_ERI_Pz_S_Py_S_C1000000001_c;
  Double I_ERI_Px_S_Pz_Px_C1000000001_c = I_ERI_Px_S_Dxz_S_C1000000001_c+CDX*I_ERI_Px_S_Pz_S_C1000000001_c;
  Double I_ERI_Py_S_Pz_Px_C1000000001_c = I_ERI_Py_S_Dxz_S_C1000000001_c+CDX*I_ERI_Py_S_Pz_S_C1000000001_c;
  Double I_ERI_Pz_S_Pz_Px_C1000000001_c = I_ERI_Pz_S_Dxz_S_C1000000001_c+CDX*I_ERI_Pz_S_Pz_S_C1000000001_c;
  Double I_ERI_Px_S_Px_Py_C1000000001_c = I_ERI_Px_S_Dxy_S_C1000000001_c+CDY*I_ERI_Px_S_Px_S_C1000000001_c;
  Double I_ERI_Py_S_Px_Py_C1000000001_c = I_ERI_Py_S_Dxy_S_C1000000001_c+CDY*I_ERI_Py_S_Px_S_C1000000001_c;
  Double I_ERI_Pz_S_Px_Py_C1000000001_c = I_ERI_Pz_S_Dxy_S_C1000000001_c+CDY*I_ERI_Pz_S_Px_S_C1000000001_c;
  Double I_ERI_Px_S_Py_Py_C1000000001_c = I_ERI_Px_S_D2y_S_C1000000001_c+CDY*I_ERI_Px_S_Py_S_C1000000001_c;
  Double I_ERI_Py_S_Py_Py_C1000000001_c = I_ERI_Py_S_D2y_S_C1000000001_c+CDY*I_ERI_Py_S_Py_S_C1000000001_c;
  Double I_ERI_Pz_S_Py_Py_C1000000001_c = I_ERI_Pz_S_D2y_S_C1000000001_c+CDY*I_ERI_Pz_S_Py_S_C1000000001_c;
  Double I_ERI_Px_S_Pz_Py_C1000000001_c = I_ERI_Px_S_Dyz_S_C1000000001_c+CDY*I_ERI_Px_S_Pz_S_C1000000001_c;
  Double I_ERI_Py_S_Pz_Py_C1000000001_c = I_ERI_Py_S_Dyz_S_C1000000001_c+CDY*I_ERI_Py_S_Pz_S_C1000000001_c;
  Double I_ERI_Pz_S_Pz_Py_C1000000001_c = I_ERI_Pz_S_Dyz_S_C1000000001_c+CDY*I_ERI_Pz_S_Pz_S_C1000000001_c;
  Double I_ERI_Px_S_Px_Pz_C1000000001_c = I_ERI_Px_S_Dxz_S_C1000000001_c+CDZ*I_ERI_Px_S_Px_S_C1000000001_c;
  Double I_ERI_Py_S_Px_Pz_C1000000001_c = I_ERI_Py_S_Dxz_S_C1000000001_c+CDZ*I_ERI_Py_S_Px_S_C1000000001_c;
  Double I_ERI_Pz_S_Px_Pz_C1000000001_c = I_ERI_Pz_S_Dxz_S_C1000000001_c+CDZ*I_ERI_Pz_S_Px_S_C1000000001_c;
  Double I_ERI_Px_S_Py_Pz_C1000000001_c = I_ERI_Px_S_Dyz_S_C1000000001_c+CDZ*I_ERI_Px_S_Py_S_C1000000001_c;
  Double I_ERI_Py_S_Py_Pz_C1000000001_c = I_ERI_Py_S_Dyz_S_C1000000001_c+CDZ*I_ERI_Py_S_Py_S_C1000000001_c;
  Double I_ERI_Pz_S_Py_Pz_C1000000001_c = I_ERI_Pz_S_Dyz_S_C1000000001_c+CDZ*I_ERI_Pz_S_Py_S_C1000000001_c;
  Double I_ERI_Px_S_Pz_Pz_C1000000001_c = I_ERI_Px_S_D2z_S_C1000000001_c+CDZ*I_ERI_Px_S_Pz_S_C1000000001_c;
  Double I_ERI_Py_S_Pz_Pz_C1000000001_c = I_ERI_Py_S_D2z_S_C1000000001_c+CDZ*I_ERI_Py_S_Pz_S_C1000000001_c;
  Double I_ERI_Pz_S_Pz_Pz_C1000000001_c = I_ERI_Pz_S_D2z_S_C1000000001_c+CDZ*I_ERI_Pz_S_Pz_S_C1000000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1000001000_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1000001000_c
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000001000_c
   ************************************************************/
  Double I_ERI_S_Px_Px_Px_C1000001000_c = I_ERI_S_Px_D2x_S_C1000001000_c+CDX*I_ERI_S_Px_Px_S_C1000001000_c;
  Double I_ERI_S_Py_Px_Px_C1000001000_c = I_ERI_S_Py_D2x_S_C1000001000_c+CDX*I_ERI_S_Py_Px_S_C1000001000_c;
  Double I_ERI_S_Pz_Px_Px_C1000001000_c = I_ERI_S_Pz_D2x_S_C1000001000_c+CDX*I_ERI_S_Pz_Px_S_C1000001000_c;
  Double I_ERI_S_Px_Py_Px_C1000001000_c = I_ERI_S_Px_Dxy_S_C1000001000_c+CDX*I_ERI_S_Px_Py_S_C1000001000_c;
  Double I_ERI_S_Py_Py_Px_C1000001000_c = I_ERI_S_Py_Dxy_S_C1000001000_c+CDX*I_ERI_S_Py_Py_S_C1000001000_c;
  Double I_ERI_S_Pz_Py_Px_C1000001000_c = I_ERI_S_Pz_Dxy_S_C1000001000_c+CDX*I_ERI_S_Pz_Py_S_C1000001000_c;
  Double I_ERI_S_Px_Pz_Px_C1000001000_c = I_ERI_S_Px_Dxz_S_C1000001000_c+CDX*I_ERI_S_Px_Pz_S_C1000001000_c;
  Double I_ERI_S_Py_Pz_Px_C1000001000_c = I_ERI_S_Py_Dxz_S_C1000001000_c+CDX*I_ERI_S_Py_Pz_S_C1000001000_c;
  Double I_ERI_S_Pz_Pz_Px_C1000001000_c = I_ERI_S_Pz_Dxz_S_C1000001000_c+CDX*I_ERI_S_Pz_Pz_S_C1000001000_c;
  Double I_ERI_S_Px_Px_Py_C1000001000_c = I_ERI_S_Px_Dxy_S_C1000001000_c+CDY*I_ERI_S_Px_Px_S_C1000001000_c;
  Double I_ERI_S_Py_Px_Py_C1000001000_c = I_ERI_S_Py_Dxy_S_C1000001000_c+CDY*I_ERI_S_Py_Px_S_C1000001000_c;
  Double I_ERI_S_Pz_Px_Py_C1000001000_c = I_ERI_S_Pz_Dxy_S_C1000001000_c+CDY*I_ERI_S_Pz_Px_S_C1000001000_c;
  Double I_ERI_S_Px_Py_Py_C1000001000_c = I_ERI_S_Px_D2y_S_C1000001000_c+CDY*I_ERI_S_Px_Py_S_C1000001000_c;
  Double I_ERI_S_Py_Py_Py_C1000001000_c = I_ERI_S_Py_D2y_S_C1000001000_c+CDY*I_ERI_S_Py_Py_S_C1000001000_c;
  Double I_ERI_S_Pz_Py_Py_C1000001000_c = I_ERI_S_Pz_D2y_S_C1000001000_c+CDY*I_ERI_S_Pz_Py_S_C1000001000_c;
  Double I_ERI_S_Px_Pz_Py_C1000001000_c = I_ERI_S_Px_Dyz_S_C1000001000_c+CDY*I_ERI_S_Px_Pz_S_C1000001000_c;
  Double I_ERI_S_Py_Pz_Py_C1000001000_c = I_ERI_S_Py_Dyz_S_C1000001000_c+CDY*I_ERI_S_Py_Pz_S_C1000001000_c;
  Double I_ERI_S_Pz_Pz_Py_C1000001000_c = I_ERI_S_Pz_Dyz_S_C1000001000_c+CDY*I_ERI_S_Pz_Pz_S_C1000001000_c;
  Double I_ERI_S_Px_Px_Pz_C1000001000_c = I_ERI_S_Px_Dxz_S_C1000001000_c+CDZ*I_ERI_S_Px_Px_S_C1000001000_c;
  Double I_ERI_S_Py_Px_Pz_C1000001000_c = I_ERI_S_Py_Dxz_S_C1000001000_c+CDZ*I_ERI_S_Py_Px_S_C1000001000_c;
  Double I_ERI_S_Pz_Px_Pz_C1000001000_c = I_ERI_S_Pz_Dxz_S_C1000001000_c+CDZ*I_ERI_S_Pz_Px_S_C1000001000_c;
  Double I_ERI_S_Px_Py_Pz_C1000001000_c = I_ERI_S_Px_Dyz_S_C1000001000_c+CDZ*I_ERI_S_Px_Py_S_C1000001000_c;
  Double I_ERI_S_Py_Py_Pz_C1000001000_c = I_ERI_S_Py_Dyz_S_C1000001000_c+CDZ*I_ERI_S_Py_Py_S_C1000001000_c;
  Double I_ERI_S_Pz_Py_Pz_C1000001000_c = I_ERI_S_Pz_Dyz_S_C1000001000_c+CDZ*I_ERI_S_Pz_Py_S_C1000001000_c;
  Double I_ERI_S_Px_Pz_Pz_C1000001000_c = I_ERI_S_Px_D2z_S_C1000001000_c+CDZ*I_ERI_S_Px_Pz_S_C1000001000_c;
  Double I_ERI_S_Py_Pz_Pz_C1000001000_c = I_ERI_S_Py_D2z_S_C1000001000_c+CDZ*I_ERI_S_Py_Pz_S_C1000001000_c;
  Double I_ERI_S_Pz_Pz_Pz_C1000001000_c = I_ERI_S_Pz_D2z_S_C1000001000_c+CDZ*I_ERI_S_Pz_Pz_S_C1000001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1000001001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001001_c
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1000001001_c = I_ERI_Px_S_D2x_S_C1000001001_c+CDX*I_ERI_Px_S_Px_S_C1000001001_c;
  Double I_ERI_Py_S_Px_Px_C1000001001_c = I_ERI_Py_S_D2x_S_C1000001001_c+CDX*I_ERI_Py_S_Px_S_C1000001001_c;
  Double I_ERI_Pz_S_Px_Px_C1000001001_c = I_ERI_Pz_S_D2x_S_C1000001001_c+CDX*I_ERI_Pz_S_Px_S_C1000001001_c;
  Double I_ERI_Px_S_Py_Px_C1000001001_c = I_ERI_Px_S_Dxy_S_C1000001001_c+CDX*I_ERI_Px_S_Py_S_C1000001001_c;
  Double I_ERI_Py_S_Py_Px_C1000001001_c = I_ERI_Py_S_Dxy_S_C1000001001_c+CDX*I_ERI_Py_S_Py_S_C1000001001_c;
  Double I_ERI_Pz_S_Py_Px_C1000001001_c = I_ERI_Pz_S_Dxy_S_C1000001001_c+CDX*I_ERI_Pz_S_Py_S_C1000001001_c;
  Double I_ERI_Px_S_Pz_Px_C1000001001_c = I_ERI_Px_S_Dxz_S_C1000001001_c+CDX*I_ERI_Px_S_Pz_S_C1000001001_c;
  Double I_ERI_Py_S_Pz_Px_C1000001001_c = I_ERI_Py_S_Dxz_S_C1000001001_c+CDX*I_ERI_Py_S_Pz_S_C1000001001_c;
  Double I_ERI_Pz_S_Pz_Px_C1000001001_c = I_ERI_Pz_S_Dxz_S_C1000001001_c+CDX*I_ERI_Pz_S_Pz_S_C1000001001_c;
  Double I_ERI_Px_S_Px_Py_C1000001001_c = I_ERI_Px_S_Dxy_S_C1000001001_c+CDY*I_ERI_Px_S_Px_S_C1000001001_c;
  Double I_ERI_Py_S_Px_Py_C1000001001_c = I_ERI_Py_S_Dxy_S_C1000001001_c+CDY*I_ERI_Py_S_Px_S_C1000001001_c;
  Double I_ERI_Pz_S_Px_Py_C1000001001_c = I_ERI_Pz_S_Dxy_S_C1000001001_c+CDY*I_ERI_Pz_S_Px_S_C1000001001_c;
  Double I_ERI_Px_S_Py_Py_C1000001001_c = I_ERI_Px_S_D2y_S_C1000001001_c+CDY*I_ERI_Px_S_Py_S_C1000001001_c;
  Double I_ERI_Py_S_Py_Py_C1000001001_c = I_ERI_Py_S_D2y_S_C1000001001_c+CDY*I_ERI_Py_S_Py_S_C1000001001_c;
  Double I_ERI_Pz_S_Py_Py_C1000001001_c = I_ERI_Pz_S_D2y_S_C1000001001_c+CDY*I_ERI_Pz_S_Py_S_C1000001001_c;
  Double I_ERI_Px_S_Pz_Py_C1000001001_c = I_ERI_Px_S_Dyz_S_C1000001001_c+CDY*I_ERI_Px_S_Pz_S_C1000001001_c;
  Double I_ERI_Py_S_Pz_Py_C1000001001_c = I_ERI_Py_S_Dyz_S_C1000001001_c+CDY*I_ERI_Py_S_Pz_S_C1000001001_c;
  Double I_ERI_Pz_S_Pz_Py_C1000001001_c = I_ERI_Pz_S_Dyz_S_C1000001001_c+CDY*I_ERI_Pz_S_Pz_S_C1000001001_c;
  Double I_ERI_Px_S_Px_Pz_C1000001001_c = I_ERI_Px_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Px_S_Px_S_C1000001001_c;
  Double I_ERI_Py_S_Px_Pz_C1000001001_c = I_ERI_Py_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Py_S_Px_S_C1000001001_c;
  Double I_ERI_Pz_S_Px_Pz_C1000001001_c = I_ERI_Pz_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Pz_S_Px_S_C1000001001_c;
  Double I_ERI_Px_S_Py_Pz_C1000001001_c = I_ERI_Px_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Px_S_Py_S_C1000001001_c;
  Double I_ERI_Py_S_Py_Pz_C1000001001_c = I_ERI_Py_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Py_S_Py_S_C1000001001_c;
  Double I_ERI_Pz_S_Py_Pz_C1000001001_c = I_ERI_Pz_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Pz_S_Py_S_C1000001001_c;
  Double I_ERI_Px_S_Pz_Pz_C1000001001_c = I_ERI_Px_S_D2z_S_C1000001001_c+CDZ*I_ERI_Px_S_Pz_S_C1000001001_c;
  Double I_ERI_Py_S_Pz_Pz_C1000001001_c = I_ERI_Py_S_D2z_S_C1000001001_c+CDZ*I_ERI_Py_S_Pz_S_C1000001001_c;
  Double I_ERI_Pz_S_Pz_Pz_C1000001001_c = I_ERI_Pz_S_D2z_S_C1000001001_c+CDZ*I_ERI_Pz_S_Pz_S_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1000001001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000001001_c
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001001_c
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1000001001_c = I_ERI_D2x_S_D2x_S_C1000001001_c+CDX*I_ERI_D2x_S_Px_S_C1000001001_c;
  Double I_ERI_Dxy_S_Px_Px_C1000001001_c = I_ERI_Dxy_S_D2x_S_C1000001001_c+CDX*I_ERI_Dxy_S_Px_S_C1000001001_c;
  Double I_ERI_Dxz_S_Px_Px_C1000001001_c = I_ERI_Dxz_S_D2x_S_C1000001001_c+CDX*I_ERI_Dxz_S_Px_S_C1000001001_c;
  Double I_ERI_D2y_S_Px_Px_C1000001001_c = I_ERI_D2y_S_D2x_S_C1000001001_c+CDX*I_ERI_D2y_S_Px_S_C1000001001_c;
  Double I_ERI_Dyz_S_Px_Px_C1000001001_c = I_ERI_Dyz_S_D2x_S_C1000001001_c+CDX*I_ERI_Dyz_S_Px_S_C1000001001_c;
  Double I_ERI_D2z_S_Px_Px_C1000001001_c = I_ERI_D2z_S_D2x_S_C1000001001_c+CDX*I_ERI_D2z_S_Px_S_C1000001001_c;
  Double I_ERI_D2x_S_Py_Px_C1000001001_c = I_ERI_D2x_S_Dxy_S_C1000001001_c+CDX*I_ERI_D2x_S_Py_S_C1000001001_c;
  Double I_ERI_Dxy_S_Py_Px_C1000001001_c = I_ERI_Dxy_S_Dxy_S_C1000001001_c+CDX*I_ERI_Dxy_S_Py_S_C1000001001_c;
  Double I_ERI_Dxz_S_Py_Px_C1000001001_c = I_ERI_Dxz_S_Dxy_S_C1000001001_c+CDX*I_ERI_Dxz_S_Py_S_C1000001001_c;
  Double I_ERI_D2y_S_Py_Px_C1000001001_c = I_ERI_D2y_S_Dxy_S_C1000001001_c+CDX*I_ERI_D2y_S_Py_S_C1000001001_c;
  Double I_ERI_Dyz_S_Py_Px_C1000001001_c = I_ERI_Dyz_S_Dxy_S_C1000001001_c+CDX*I_ERI_Dyz_S_Py_S_C1000001001_c;
  Double I_ERI_D2z_S_Py_Px_C1000001001_c = I_ERI_D2z_S_Dxy_S_C1000001001_c+CDX*I_ERI_D2z_S_Py_S_C1000001001_c;
  Double I_ERI_D2x_S_Pz_Px_C1000001001_c = I_ERI_D2x_S_Dxz_S_C1000001001_c+CDX*I_ERI_D2x_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxy_S_Pz_Px_C1000001001_c = I_ERI_Dxy_S_Dxz_S_C1000001001_c+CDX*I_ERI_Dxy_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxz_S_Pz_Px_C1000001001_c = I_ERI_Dxz_S_Dxz_S_C1000001001_c+CDX*I_ERI_Dxz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2y_S_Pz_Px_C1000001001_c = I_ERI_D2y_S_Dxz_S_C1000001001_c+CDX*I_ERI_D2y_S_Pz_S_C1000001001_c;
  Double I_ERI_Dyz_S_Pz_Px_C1000001001_c = I_ERI_Dyz_S_Dxz_S_C1000001001_c+CDX*I_ERI_Dyz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2z_S_Pz_Px_C1000001001_c = I_ERI_D2z_S_Dxz_S_C1000001001_c+CDX*I_ERI_D2z_S_Pz_S_C1000001001_c;
  Double I_ERI_D2x_S_Px_Py_C1000001001_c = I_ERI_D2x_S_Dxy_S_C1000001001_c+CDY*I_ERI_D2x_S_Px_S_C1000001001_c;
  Double I_ERI_Dxy_S_Px_Py_C1000001001_c = I_ERI_Dxy_S_Dxy_S_C1000001001_c+CDY*I_ERI_Dxy_S_Px_S_C1000001001_c;
  Double I_ERI_Dxz_S_Px_Py_C1000001001_c = I_ERI_Dxz_S_Dxy_S_C1000001001_c+CDY*I_ERI_Dxz_S_Px_S_C1000001001_c;
  Double I_ERI_D2y_S_Px_Py_C1000001001_c = I_ERI_D2y_S_Dxy_S_C1000001001_c+CDY*I_ERI_D2y_S_Px_S_C1000001001_c;
  Double I_ERI_Dyz_S_Px_Py_C1000001001_c = I_ERI_Dyz_S_Dxy_S_C1000001001_c+CDY*I_ERI_Dyz_S_Px_S_C1000001001_c;
  Double I_ERI_D2z_S_Px_Py_C1000001001_c = I_ERI_D2z_S_Dxy_S_C1000001001_c+CDY*I_ERI_D2z_S_Px_S_C1000001001_c;
  Double I_ERI_D2x_S_Py_Py_C1000001001_c = I_ERI_D2x_S_D2y_S_C1000001001_c+CDY*I_ERI_D2x_S_Py_S_C1000001001_c;
  Double I_ERI_Dxy_S_Py_Py_C1000001001_c = I_ERI_Dxy_S_D2y_S_C1000001001_c+CDY*I_ERI_Dxy_S_Py_S_C1000001001_c;
  Double I_ERI_Dxz_S_Py_Py_C1000001001_c = I_ERI_Dxz_S_D2y_S_C1000001001_c+CDY*I_ERI_Dxz_S_Py_S_C1000001001_c;
  Double I_ERI_D2y_S_Py_Py_C1000001001_c = I_ERI_D2y_S_D2y_S_C1000001001_c+CDY*I_ERI_D2y_S_Py_S_C1000001001_c;
  Double I_ERI_Dyz_S_Py_Py_C1000001001_c = I_ERI_Dyz_S_D2y_S_C1000001001_c+CDY*I_ERI_Dyz_S_Py_S_C1000001001_c;
  Double I_ERI_D2z_S_Py_Py_C1000001001_c = I_ERI_D2z_S_D2y_S_C1000001001_c+CDY*I_ERI_D2z_S_Py_S_C1000001001_c;
  Double I_ERI_D2x_S_Pz_Py_C1000001001_c = I_ERI_D2x_S_Dyz_S_C1000001001_c+CDY*I_ERI_D2x_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxy_S_Pz_Py_C1000001001_c = I_ERI_Dxy_S_Dyz_S_C1000001001_c+CDY*I_ERI_Dxy_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxz_S_Pz_Py_C1000001001_c = I_ERI_Dxz_S_Dyz_S_C1000001001_c+CDY*I_ERI_Dxz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2y_S_Pz_Py_C1000001001_c = I_ERI_D2y_S_Dyz_S_C1000001001_c+CDY*I_ERI_D2y_S_Pz_S_C1000001001_c;
  Double I_ERI_Dyz_S_Pz_Py_C1000001001_c = I_ERI_Dyz_S_Dyz_S_C1000001001_c+CDY*I_ERI_Dyz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2z_S_Pz_Py_C1000001001_c = I_ERI_D2z_S_Dyz_S_C1000001001_c+CDY*I_ERI_D2z_S_Pz_S_C1000001001_c;
  Double I_ERI_D2x_S_Px_Pz_C1000001001_c = I_ERI_D2x_S_Dxz_S_C1000001001_c+CDZ*I_ERI_D2x_S_Px_S_C1000001001_c;
  Double I_ERI_Dxy_S_Px_Pz_C1000001001_c = I_ERI_Dxy_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Dxy_S_Px_S_C1000001001_c;
  Double I_ERI_Dxz_S_Px_Pz_C1000001001_c = I_ERI_Dxz_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Dxz_S_Px_S_C1000001001_c;
  Double I_ERI_D2y_S_Px_Pz_C1000001001_c = I_ERI_D2y_S_Dxz_S_C1000001001_c+CDZ*I_ERI_D2y_S_Px_S_C1000001001_c;
  Double I_ERI_Dyz_S_Px_Pz_C1000001001_c = I_ERI_Dyz_S_Dxz_S_C1000001001_c+CDZ*I_ERI_Dyz_S_Px_S_C1000001001_c;
  Double I_ERI_D2z_S_Px_Pz_C1000001001_c = I_ERI_D2z_S_Dxz_S_C1000001001_c+CDZ*I_ERI_D2z_S_Px_S_C1000001001_c;
  Double I_ERI_D2x_S_Py_Pz_C1000001001_c = I_ERI_D2x_S_Dyz_S_C1000001001_c+CDZ*I_ERI_D2x_S_Py_S_C1000001001_c;
  Double I_ERI_Dxy_S_Py_Pz_C1000001001_c = I_ERI_Dxy_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Dxy_S_Py_S_C1000001001_c;
  Double I_ERI_Dxz_S_Py_Pz_C1000001001_c = I_ERI_Dxz_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Dxz_S_Py_S_C1000001001_c;
  Double I_ERI_D2y_S_Py_Pz_C1000001001_c = I_ERI_D2y_S_Dyz_S_C1000001001_c+CDZ*I_ERI_D2y_S_Py_S_C1000001001_c;
  Double I_ERI_Dyz_S_Py_Pz_C1000001001_c = I_ERI_Dyz_S_Dyz_S_C1000001001_c+CDZ*I_ERI_Dyz_S_Py_S_C1000001001_c;
  Double I_ERI_D2z_S_Py_Pz_C1000001001_c = I_ERI_D2z_S_Dyz_S_C1000001001_c+CDZ*I_ERI_D2z_S_Py_S_C1000001001_c;
  Double I_ERI_D2x_S_Pz_Pz_C1000001001_c = I_ERI_D2x_S_D2z_S_C1000001001_c+CDZ*I_ERI_D2x_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxy_S_Pz_Pz_C1000001001_c = I_ERI_Dxy_S_D2z_S_C1000001001_c+CDZ*I_ERI_Dxy_S_Pz_S_C1000001001_c;
  Double I_ERI_Dxz_S_Pz_Pz_C1000001001_c = I_ERI_Dxz_S_D2z_S_C1000001001_c+CDZ*I_ERI_Dxz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2y_S_Pz_Pz_C1000001001_c = I_ERI_D2y_S_D2z_S_C1000001001_c+CDZ*I_ERI_D2y_S_Pz_S_C1000001001_c;
  Double I_ERI_Dyz_S_Pz_Pz_C1000001001_c = I_ERI_Dyz_S_D2z_S_C1000001001_c+CDZ*I_ERI_Dyz_S_Pz_S_C1000001001_c;
  Double I_ERI_D2z_S_Pz_Pz_C1000001001_c = I_ERI_D2z_S_D2z_S_C1000001001_c+CDZ*I_ERI_D2z_S_Pz_S_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_D_P_C1001000000_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_F_S_C1001000000_c
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1001000000_c
   ************************************************************/
  Double I_ERI_S_S_D2x_Px_C1001000000_c = I_ERI_S_S_F3x_S_C1001000000_c+CDX*I_ERI_S_S_D2x_S_C1001000000_c;
  Double I_ERI_S_S_Dxy_Px_C1001000000_c = I_ERI_S_S_F2xy_S_C1001000000_c+CDX*I_ERI_S_S_Dxy_S_C1001000000_c;
  Double I_ERI_S_S_Dxz_Px_C1001000000_c = I_ERI_S_S_F2xz_S_C1001000000_c+CDX*I_ERI_S_S_Dxz_S_C1001000000_c;
  Double I_ERI_S_S_D2y_Px_C1001000000_c = I_ERI_S_S_Fx2y_S_C1001000000_c+CDX*I_ERI_S_S_D2y_S_C1001000000_c;
  Double I_ERI_S_S_Dyz_Px_C1001000000_c = I_ERI_S_S_Fxyz_S_C1001000000_c+CDX*I_ERI_S_S_Dyz_S_C1001000000_c;
  Double I_ERI_S_S_D2z_Px_C1001000000_c = I_ERI_S_S_Fx2z_S_C1001000000_c+CDX*I_ERI_S_S_D2z_S_C1001000000_c;
  Double I_ERI_S_S_D2x_Py_C1001000000_c = I_ERI_S_S_F2xy_S_C1001000000_c+CDY*I_ERI_S_S_D2x_S_C1001000000_c;
  Double I_ERI_S_S_Dxy_Py_C1001000000_c = I_ERI_S_S_Fx2y_S_C1001000000_c+CDY*I_ERI_S_S_Dxy_S_C1001000000_c;
  Double I_ERI_S_S_Dxz_Py_C1001000000_c = I_ERI_S_S_Fxyz_S_C1001000000_c+CDY*I_ERI_S_S_Dxz_S_C1001000000_c;
  Double I_ERI_S_S_D2y_Py_C1001000000_c = I_ERI_S_S_F3y_S_C1001000000_c+CDY*I_ERI_S_S_D2y_S_C1001000000_c;
  Double I_ERI_S_S_Dyz_Py_C1001000000_c = I_ERI_S_S_F2yz_S_C1001000000_c+CDY*I_ERI_S_S_Dyz_S_C1001000000_c;
  Double I_ERI_S_S_D2z_Py_C1001000000_c = I_ERI_S_S_Fy2z_S_C1001000000_c+CDY*I_ERI_S_S_D2z_S_C1001000000_c;
  Double I_ERI_S_S_D2x_Pz_C1001000000_c = I_ERI_S_S_F2xz_S_C1001000000_c+CDZ*I_ERI_S_S_D2x_S_C1001000000_c;
  Double I_ERI_S_S_Dxy_Pz_C1001000000_c = I_ERI_S_S_Fxyz_S_C1001000000_c+CDZ*I_ERI_S_S_Dxy_S_C1001000000_c;
  Double I_ERI_S_S_Dxz_Pz_C1001000000_c = I_ERI_S_S_Fx2z_S_C1001000000_c+CDZ*I_ERI_S_S_Dxz_S_C1001000000_c;
  Double I_ERI_S_S_D2y_Pz_C1001000000_c = I_ERI_S_S_F2yz_S_C1001000000_c+CDZ*I_ERI_S_S_D2y_S_C1001000000_c;
  Double I_ERI_S_S_Dyz_Pz_C1001000000_c = I_ERI_S_S_Fy2z_S_C1001000000_c+CDZ*I_ERI_S_S_Dyz_S_C1001000000_c;
  Double I_ERI_S_S_D2z_Pz_C1001000000_c = I_ERI_S_S_F3z_S_C1001000000_c+CDZ*I_ERI_S_S_D2z_S_C1001000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1001000001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1001000001_c
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001000001_c
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1001000001_c = I_ERI_Px_S_F3x_S_C1001000001_c+CDX*I_ERI_Px_S_D2x_S_C1001000001_c;
  Double I_ERI_Py_S_D2x_Px_C1001000001_c = I_ERI_Py_S_F3x_S_C1001000001_c+CDX*I_ERI_Py_S_D2x_S_C1001000001_c;
  Double I_ERI_Pz_S_D2x_Px_C1001000001_c = I_ERI_Pz_S_F3x_S_C1001000001_c+CDX*I_ERI_Pz_S_D2x_S_C1001000001_c;
  Double I_ERI_Px_S_Dxy_Px_C1001000001_c = I_ERI_Px_S_F2xy_S_C1001000001_c+CDX*I_ERI_Px_S_Dxy_S_C1001000001_c;
  Double I_ERI_Py_S_Dxy_Px_C1001000001_c = I_ERI_Py_S_F2xy_S_C1001000001_c+CDX*I_ERI_Py_S_Dxy_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxy_Px_C1001000001_c = I_ERI_Pz_S_F2xy_S_C1001000001_c+CDX*I_ERI_Pz_S_Dxy_S_C1001000001_c;
  Double I_ERI_Px_S_Dxz_Px_C1001000001_c = I_ERI_Px_S_F2xz_S_C1001000001_c+CDX*I_ERI_Px_S_Dxz_S_C1001000001_c;
  Double I_ERI_Py_S_Dxz_Px_C1001000001_c = I_ERI_Py_S_F2xz_S_C1001000001_c+CDX*I_ERI_Py_S_Dxz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxz_Px_C1001000001_c = I_ERI_Pz_S_F2xz_S_C1001000001_c+CDX*I_ERI_Pz_S_Dxz_S_C1001000001_c;
  Double I_ERI_Px_S_D2y_Px_C1001000001_c = I_ERI_Px_S_Fx2y_S_C1001000001_c+CDX*I_ERI_Px_S_D2y_S_C1001000001_c;
  Double I_ERI_Py_S_D2y_Px_C1001000001_c = I_ERI_Py_S_Fx2y_S_C1001000001_c+CDX*I_ERI_Py_S_D2y_S_C1001000001_c;
  Double I_ERI_Pz_S_D2y_Px_C1001000001_c = I_ERI_Pz_S_Fx2y_S_C1001000001_c+CDX*I_ERI_Pz_S_D2y_S_C1001000001_c;
  Double I_ERI_Px_S_Dyz_Px_C1001000001_c = I_ERI_Px_S_Fxyz_S_C1001000001_c+CDX*I_ERI_Px_S_Dyz_S_C1001000001_c;
  Double I_ERI_Py_S_Dyz_Px_C1001000001_c = I_ERI_Py_S_Fxyz_S_C1001000001_c+CDX*I_ERI_Py_S_Dyz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dyz_Px_C1001000001_c = I_ERI_Pz_S_Fxyz_S_C1001000001_c+CDX*I_ERI_Pz_S_Dyz_S_C1001000001_c;
  Double I_ERI_Px_S_D2z_Px_C1001000001_c = I_ERI_Px_S_Fx2z_S_C1001000001_c+CDX*I_ERI_Px_S_D2z_S_C1001000001_c;
  Double I_ERI_Py_S_D2z_Px_C1001000001_c = I_ERI_Py_S_Fx2z_S_C1001000001_c+CDX*I_ERI_Py_S_D2z_S_C1001000001_c;
  Double I_ERI_Pz_S_D2z_Px_C1001000001_c = I_ERI_Pz_S_Fx2z_S_C1001000001_c+CDX*I_ERI_Pz_S_D2z_S_C1001000001_c;
  Double I_ERI_Px_S_D2x_Py_C1001000001_c = I_ERI_Px_S_F2xy_S_C1001000001_c+CDY*I_ERI_Px_S_D2x_S_C1001000001_c;
  Double I_ERI_Py_S_D2x_Py_C1001000001_c = I_ERI_Py_S_F2xy_S_C1001000001_c+CDY*I_ERI_Py_S_D2x_S_C1001000001_c;
  Double I_ERI_Pz_S_D2x_Py_C1001000001_c = I_ERI_Pz_S_F2xy_S_C1001000001_c+CDY*I_ERI_Pz_S_D2x_S_C1001000001_c;
  Double I_ERI_Px_S_Dxy_Py_C1001000001_c = I_ERI_Px_S_Fx2y_S_C1001000001_c+CDY*I_ERI_Px_S_Dxy_S_C1001000001_c;
  Double I_ERI_Py_S_Dxy_Py_C1001000001_c = I_ERI_Py_S_Fx2y_S_C1001000001_c+CDY*I_ERI_Py_S_Dxy_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxy_Py_C1001000001_c = I_ERI_Pz_S_Fx2y_S_C1001000001_c+CDY*I_ERI_Pz_S_Dxy_S_C1001000001_c;
  Double I_ERI_Px_S_Dxz_Py_C1001000001_c = I_ERI_Px_S_Fxyz_S_C1001000001_c+CDY*I_ERI_Px_S_Dxz_S_C1001000001_c;
  Double I_ERI_Py_S_Dxz_Py_C1001000001_c = I_ERI_Py_S_Fxyz_S_C1001000001_c+CDY*I_ERI_Py_S_Dxz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxz_Py_C1001000001_c = I_ERI_Pz_S_Fxyz_S_C1001000001_c+CDY*I_ERI_Pz_S_Dxz_S_C1001000001_c;
  Double I_ERI_Px_S_D2y_Py_C1001000001_c = I_ERI_Px_S_F3y_S_C1001000001_c+CDY*I_ERI_Px_S_D2y_S_C1001000001_c;
  Double I_ERI_Py_S_D2y_Py_C1001000001_c = I_ERI_Py_S_F3y_S_C1001000001_c+CDY*I_ERI_Py_S_D2y_S_C1001000001_c;
  Double I_ERI_Pz_S_D2y_Py_C1001000001_c = I_ERI_Pz_S_F3y_S_C1001000001_c+CDY*I_ERI_Pz_S_D2y_S_C1001000001_c;
  Double I_ERI_Px_S_Dyz_Py_C1001000001_c = I_ERI_Px_S_F2yz_S_C1001000001_c+CDY*I_ERI_Px_S_Dyz_S_C1001000001_c;
  Double I_ERI_Py_S_Dyz_Py_C1001000001_c = I_ERI_Py_S_F2yz_S_C1001000001_c+CDY*I_ERI_Py_S_Dyz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dyz_Py_C1001000001_c = I_ERI_Pz_S_F2yz_S_C1001000001_c+CDY*I_ERI_Pz_S_Dyz_S_C1001000001_c;
  Double I_ERI_Px_S_D2z_Py_C1001000001_c = I_ERI_Px_S_Fy2z_S_C1001000001_c+CDY*I_ERI_Px_S_D2z_S_C1001000001_c;
  Double I_ERI_Py_S_D2z_Py_C1001000001_c = I_ERI_Py_S_Fy2z_S_C1001000001_c+CDY*I_ERI_Py_S_D2z_S_C1001000001_c;
  Double I_ERI_Pz_S_D2z_Py_C1001000001_c = I_ERI_Pz_S_Fy2z_S_C1001000001_c+CDY*I_ERI_Pz_S_D2z_S_C1001000001_c;
  Double I_ERI_Px_S_D2x_Pz_C1001000001_c = I_ERI_Px_S_F2xz_S_C1001000001_c+CDZ*I_ERI_Px_S_D2x_S_C1001000001_c;
  Double I_ERI_Py_S_D2x_Pz_C1001000001_c = I_ERI_Py_S_F2xz_S_C1001000001_c+CDZ*I_ERI_Py_S_D2x_S_C1001000001_c;
  Double I_ERI_Pz_S_D2x_Pz_C1001000001_c = I_ERI_Pz_S_F2xz_S_C1001000001_c+CDZ*I_ERI_Pz_S_D2x_S_C1001000001_c;
  Double I_ERI_Px_S_Dxy_Pz_C1001000001_c = I_ERI_Px_S_Fxyz_S_C1001000001_c+CDZ*I_ERI_Px_S_Dxy_S_C1001000001_c;
  Double I_ERI_Py_S_Dxy_Pz_C1001000001_c = I_ERI_Py_S_Fxyz_S_C1001000001_c+CDZ*I_ERI_Py_S_Dxy_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxy_Pz_C1001000001_c = I_ERI_Pz_S_Fxyz_S_C1001000001_c+CDZ*I_ERI_Pz_S_Dxy_S_C1001000001_c;
  Double I_ERI_Px_S_Dxz_Pz_C1001000001_c = I_ERI_Px_S_Fx2z_S_C1001000001_c+CDZ*I_ERI_Px_S_Dxz_S_C1001000001_c;
  Double I_ERI_Py_S_Dxz_Pz_C1001000001_c = I_ERI_Py_S_Fx2z_S_C1001000001_c+CDZ*I_ERI_Py_S_Dxz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dxz_Pz_C1001000001_c = I_ERI_Pz_S_Fx2z_S_C1001000001_c+CDZ*I_ERI_Pz_S_Dxz_S_C1001000001_c;
  Double I_ERI_Px_S_D2y_Pz_C1001000001_c = I_ERI_Px_S_F2yz_S_C1001000001_c+CDZ*I_ERI_Px_S_D2y_S_C1001000001_c;
  Double I_ERI_Py_S_D2y_Pz_C1001000001_c = I_ERI_Py_S_F2yz_S_C1001000001_c+CDZ*I_ERI_Py_S_D2y_S_C1001000001_c;
  Double I_ERI_Pz_S_D2y_Pz_C1001000001_c = I_ERI_Pz_S_F2yz_S_C1001000001_c+CDZ*I_ERI_Pz_S_D2y_S_C1001000001_c;
  Double I_ERI_Px_S_Dyz_Pz_C1001000001_c = I_ERI_Px_S_Fy2z_S_C1001000001_c+CDZ*I_ERI_Px_S_Dyz_S_C1001000001_c;
  Double I_ERI_Py_S_Dyz_Pz_C1001000001_c = I_ERI_Py_S_Fy2z_S_C1001000001_c+CDZ*I_ERI_Py_S_Dyz_S_C1001000001_c;
  Double I_ERI_Pz_S_Dyz_Pz_C1001000001_c = I_ERI_Pz_S_Fy2z_S_C1001000001_c+CDZ*I_ERI_Pz_S_Dyz_S_C1001000001_c;
  Double I_ERI_Px_S_D2z_Pz_C1001000001_c = I_ERI_Px_S_F3z_S_C1001000001_c+CDZ*I_ERI_Px_S_D2z_S_C1001000001_c;
  Double I_ERI_Py_S_D2z_Pz_C1001000001_c = I_ERI_Py_S_F3z_S_C1001000001_c+CDZ*I_ERI_Py_S_D2z_S_C1001000001_c;
  Double I_ERI_Pz_S_D2z_Pz_C1001000001_c = I_ERI_Pz_S_F3z_S_C1001000001_c+CDZ*I_ERI_Pz_S_D2z_S_C1001000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_D_P_C1001001000_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_F_S_C1001001000_c
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001001000_c
   ************************************************************/
  Double I_ERI_S_Px_D2x_Px_C1001001000_c = I_ERI_S_Px_F3x_S_C1001001000_c+CDX*I_ERI_S_Px_D2x_S_C1001001000_c;
  Double I_ERI_S_Py_D2x_Px_C1001001000_c = I_ERI_S_Py_F3x_S_C1001001000_c+CDX*I_ERI_S_Py_D2x_S_C1001001000_c;
  Double I_ERI_S_Pz_D2x_Px_C1001001000_c = I_ERI_S_Pz_F3x_S_C1001001000_c+CDX*I_ERI_S_Pz_D2x_S_C1001001000_c;
  Double I_ERI_S_Px_Dxy_Px_C1001001000_c = I_ERI_S_Px_F2xy_S_C1001001000_c+CDX*I_ERI_S_Px_Dxy_S_C1001001000_c;
  Double I_ERI_S_Py_Dxy_Px_C1001001000_c = I_ERI_S_Py_F2xy_S_C1001001000_c+CDX*I_ERI_S_Py_Dxy_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxy_Px_C1001001000_c = I_ERI_S_Pz_F2xy_S_C1001001000_c+CDX*I_ERI_S_Pz_Dxy_S_C1001001000_c;
  Double I_ERI_S_Px_Dxz_Px_C1001001000_c = I_ERI_S_Px_F2xz_S_C1001001000_c+CDX*I_ERI_S_Px_Dxz_S_C1001001000_c;
  Double I_ERI_S_Py_Dxz_Px_C1001001000_c = I_ERI_S_Py_F2xz_S_C1001001000_c+CDX*I_ERI_S_Py_Dxz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxz_Px_C1001001000_c = I_ERI_S_Pz_F2xz_S_C1001001000_c+CDX*I_ERI_S_Pz_Dxz_S_C1001001000_c;
  Double I_ERI_S_Px_D2y_Px_C1001001000_c = I_ERI_S_Px_Fx2y_S_C1001001000_c+CDX*I_ERI_S_Px_D2y_S_C1001001000_c;
  Double I_ERI_S_Py_D2y_Px_C1001001000_c = I_ERI_S_Py_Fx2y_S_C1001001000_c+CDX*I_ERI_S_Py_D2y_S_C1001001000_c;
  Double I_ERI_S_Pz_D2y_Px_C1001001000_c = I_ERI_S_Pz_Fx2y_S_C1001001000_c+CDX*I_ERI_S_Pz_D2y_S_C1001001000_c;
  Double I_ERI_S_Px_Dyz_Px_C1001001000_c = I_ERI_S_Px_Fxyz_S_C1001001000_c+CDX*I_ERI_S_Px_Dyz_S_C1001001000_c;
  Double I_ERI_S_Py_Dyz_Px_C1001001000_c = I_ERI_S_Py_Fxyz_S_C1001001000_c+CDX*I_ERI_S_Py_Dyz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dyz_Px_C1001001000_c = I_ERI_S_Pz_Fxyz_S_C1001001000_c+CDX*I_ERI_S_Pz_Dyz_S_C1001001000_c;
  Double I_ERI_S_Px_D2z_Px_C1001001000_c = I_ERI_S_Px_Fx2z_S_C1001001000_c+CDX*I_ERI_S_Px_D2z_S_C1001001000_c;
  Double I_ERI_S_Py_D2z_Px_C1001001000_c = I_ERI_S_Py_Fx2z_S_C1001001000_c+CDX*I_ERI_S_Py_D2z_S_C1001001000_c;
  Double I_ERI_S_Pz_D2z_Px_C1001001000_c = I_ERI_S_Pz_Fx2z_S_C1001001000_c+CDX*I_ERI_S_Pz_D2z_S_C1001001000_c;
  Double I_ERI_S_Px_D2x_Py_C1001001000_c = I_ERI_S_Px_F2xy_S_C1001001000_c+CDY*I_ERI_S_Px_D2x_S_C1001001000_c;
  Double I_ERI_S_Py_D2x_Py_C1001001000_c = I_ERI_S_Py_F2xy_S_C1001001000_c+CDY*I_ERI_S_Py_D2x_S_C1001001000_c;
  Double I_ERI_S_Pz_D2x_Py_C1001001000_c = I_ERI_S_Pz_F2xy_S_C1001001000_c+CDY*I_ERI_S_Pz_D2x_S_C1001001000_c;
  Double I_ERI_S_Px_Dxy_Py_C1001001000_c = I_ERI_S_Px_Fx2y_S_C1001001000_c+CDY*I_ERI_S_Px_Dxy_S_C1001001000_c;
  Double I_ERI_S_Py_Dxy_Py_C1001001000_c = I_ERI_S_Py_Fx2y_S_C1001001000_c+CDY*I_ERI_S_Py_Dxy_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxy_Py_C1001001000_c = I_ERI_S_Pz_Fx2y_S_C1001001000_c+CDY*I_ERI_S_Pz_Dxy_S_C1001001000_c;
  Double I_ERI_S_Px_Dxz_Py_C1001001000_c = I_ERI_S_Px_Fxyz_S_C1001001000_c+CDY*I_ERI_S_Px_Dxz_S_C1001001000_c;
  Double I_ERI_S_Py_Dxz_Py_C1001001000_c = I_ERI_S_Py_Fxyz_S_C1001001000_c+CDY*I_ERI_S_Py_Dxz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxz_Py_C1001001000_c = I_ERI_S_Pz_Fxyz_S_C1001001000_c+CDY*I_ERI_S_Pz_Dxz_S_C1001001000_c;
  Double I_ERI_S_Px_D2y_Py_C1001001000_c = I_ERI_S_Px_F3y_S_C1001001000_c+CDY*I_ERI_S_Px_D2y_S_C1001001000_c;
  Double I_ERI_S_Py_D2y_Py_C1001001000_c = I_ERI_S_Py_F3y_S_C1001001000_c+CDY*I_ERI_S_Py_D2y_S_C1001001000_c;
  Double I_ERI_S_Pz_D2y_Py_C1001001000_c = I_ERI_S_Pz_F3y_S_C1001001000_c+CDY*I_ERI_S_Pz_D2y_S_C1001001000_c;
  Double I_ERI_S_Px_Dyz_Py_C1001001000_c = I_ERI_S_Px_F2yz_S_C1001001000_c+CDY*I_ERI_S_Px_Dyz_S_C1001001000_c;
  Double I_ERI_S_Py_Dyz_Py_C1001001000_c = I_ERI_S_Py_F2yz_S_C1001001000_c+CDY*I_ERI_S_Py_Dyz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dyz_Py_C1001001000_c = I_ERI_S_Pz_F2yz_S_C1001001000_c+CDY*I_ERI_S_Pz_Dyz_S_C1001001000_c;
  Double I_ERI_S_Px_D2z_Py_C1001001000_c = I_ERI_S_Px_Fy2z_S_C1001001000_c+CDY*I_ERI_S_Px_D2z_S_C1001001000_c;
  Double I_ERI_S_Py_D2z_Py_C1001001000_c = I_ERI_S_Py_Fy2z_S_C1001001000_c+CDY*I_ERI_S_Py_D2z_S_C1001001000_c;
  Double I_ERI_S_Pz_D2z_Py_C1001001000_c = I_ERI_S_Pz_Fy2z_S_C1001001000_c+CDY*I_ERI_S_Pz_D2z_S_C1001001000_c;
  Double I_ERI_S_Px_D2x_Pz_C1001001000_c = I_ERI_S_Px_F2xz_S_C1001001000_c+CDZ*I_ERI_S_Px_D2x_S_C1001001000_c;
  Double I_ERI_S_Py_D2x_Pz_C1001001000_c = I_ERI_S_Py_F2xz_S_C1001001000_c+CDZ*I_ERI_S_Py_D2x_S_C1001001000_c;
  Double I_ERI_S_Pz_D2x_Pz_C1001001000_c = I_ERI_S_Pz_F2xz_S_C1001001000_c+CDZ*I_ERI_S_Pz_D2x_S_C1001001000_c;
  Double I_ERI_S_Px_Dxy_Pz_C1001001000_c = I_ERI_S_Px_Fxyz_S_C1001001000_c+CDZ*I_ERI_S_Px_Dxy_S_C1001001000_c;
  Double I_ERI_S_Py_Dxy_Pz_C1001001000_c = I_ERI_S_Py_Fxyz_S_C1001001000_c+CDZ*I_ERI_S_Py_Dxy_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxy_Pz_C1001001000_c = I_ERI_S_Pz_Fxyz_S_C1001001000_c+CDZ*I_ERI_S_Pz_Dxy_S_C1001001000_c;
  Double I_ERI_S_Px_Dxz_Pz_C1001001000_c = I_ERI_S_Px_Fx2z_S_C1001001000_c+CDZ*I_ERI_S_Px_Dxz_S_C1001001000_c;
  Double I_ERI_S_Py_Dxz_Pz_C1001001000_c = I_ERI_S_Py_Fx2z_S_C1001001000_c+CDZ*I_ERI_S_Py_Dxz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dxz_Pz_C1001001000_c = I_ERI_S_Pz_Fx2z_S_C1001001000_c+CDZ*I_ERI_S_Pz_Dxz_S_C1001001000_c;
  Double I_ERI_S_Px_D2y_Pz_C1001001000_c = I_ERI_S_Px_F2yz_S_C1001001000_c+CDZ*I_ERI_S_Px_D2y_S_C1001001000_c;
  Double I_ERI_S_Py_D2y_Pz_C1001001000_c = I_ERI_S_Py_F2yz_S_C1001001000_c+CDZ*I_ERI_S_Py_D2y_S_C1001001000_c;
  Double I_ERI_S_Pz_D2y_Pz_C1001001000_c = I_ERI_S_Pz_F2yz_S_C1001001000_c+CDZ*I_ERI_S_Pz_D2y_S_C1001001000_c;
  Double I_ERI_S_Px_Dyz_Pz_C1001001000_c = I_ERI_S_Px_Fy2z_S_C1001001000_c+CDZ*I_ERI_S_Px_Dyz_S_C1001001000_c;
  Double I_ERI_S_Py_Dyz_Pz_C1001001000_c = I_ERI_S_Py_Fy2z_S_C1001001000_c+CDZ*I_ERI_S_Py_Dyz_S_C1001001000_c;
  Double I_ERI_S_Pz_Dyz_Pz_C1001001000_c = I_ERI_S_Pz_Fy2z_S_C1001001000_c+CDZ*I_ERI_S_Pz_Dyz_S_C1001001000_c;
  Double I_ERI_S_Px_D2z_Pz_C1001001000_c = I_ERI_S_Px_F3z_S_C1001001000_c+CDZ*I_ERI_S_Px_D2z_S_C1001001000_c;
  Double I_ERI_S_Py_D2z_Pz_C1001001000_c = I_ERI_S_Py_F3z_S_C1001001000_c+CDZ*I_ERI_S_Py_D2z_S_C1001001000_c;
  Double I_ERI_S_Pz_D2z_Pz_C1001001000_c = I_ERI_S_Pz_F3z_S_C1001001000_c+CDZ*I_ERI_S_Pz_D2z_S_C1001001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_C1001001001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_C1001001001_c
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001001_c
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_C1001001001_c = I_ERI_Px_S_F3x_S_C1001001001_c+CDX*I_ERI_Px_S_D2x_S_C1001001001_c;
  Double I_ERI_Py_S_D2x_Px_C1001001001_c = I_ERI_Py_S_F3x_S_C1001001001_c+CDX*I_ERI_Py_S_D2x_S_C1001001001_c;
  Double I_ERI_Pz_S_D2x_Px_C1001001001_c = I_ERI_Pz_S_F3x_S_C1001001001_c+CDX*I_ERI_Pz_S_D2x_S_C1001001001_c;
  Double I_ERI_Px_S_Dxy_Px_C1001001001_c = I_ERI_Px_S_F2xy_S_C1001001001_c+CDX*I_ERI_Px_S_Dxy_S_C1001001001_c;
  Double I_ERI_Py_S_Dxy_Px_C1001001001_c = I_ERI_Py_S_F2xy_S_C1001001001_c+CDX*I_ERI_Py_S_Dxy_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxy_Px_C1001001001_c = I_ERI_Pz_S_F2xy_S_C1001001001_c+CDX*I_ERI_Pz_S_Dxy_S_C1001001001_c;
  Double I_ERI_Px_S_Dxz_Px_C1001001001_c = I_ERI_Px_S_F2xz_S_C1001001001_c+CDX*I_ERI_Px_S_Dxz_S_C1001001001_c;
  Double I_ERI_Py_S_Dxz_Px_C1001001001_c = I_ERI_Py_S_F2xz_S_C1001001001_c+CDX*I_ERI_Py_S_Dxz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxz_Px_C1001001001_c = I_ERI_Pz_S_F2xz_S_C1001001001_c+CDX*I_ERI_Pz_S_Dxz_S_C1001001001_c;
  Double I_ERI_Px_S_D2y_Px_C1001001001_c = I_ERI_Px_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Px_S_D2y_S_C1001001001_c;
  Double I_ERI_Py_S_D2y_Px_C1001001001_c = I_ERI_Py_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Py_S_D2y_S_C1001001001_c;
  Double I_ERI_Pz_S_D2y_Px_C1001001001_c = I_ERI_Pz_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Pz_S_D2y_S_C1001001001_c;
  Double I_ERI_Px_S_Dyz_Px_C1001001001_c = I_ERI_Px_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Px_S_Dyz_S_C1001001001_c;
  Double I_ERI_Py_S_Dyz_Px_C1001001001_c = I_ERI_Py_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Py_S_Dyz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dyz_Px_C1001001001_c = I_ERI_Pz_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Pz_S_Dyz_S_C1001001001_c;
  Double I_ERI_Px_S_D2z_Px_C1001001001_c = I_ERI_Px_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Px_S_D2z_S_C1001001001_c;
  Double I_ERI_Py_S_D2z_Px_C1001001001_c = I_ERI_Py_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Py_S_D2z_S_C1001001001_c;
  Double I_ERI_Pz_S_D2z_Px_C1001001001_c = I_ERI_Pz_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Pz_S_D2z_S_C1001001001_c;
  Double I_ERI_Px_S_D2x_Py_C1001001001_c = I_ERI_Px_S_F2xy_S_C1001001001_c+CDY*I_ERI_Px_S_D2x_S_C1001001001_c;
  Double I_ERI_Py_S_D2x_Py_C1001001001_c = I_ERI_Py_S_F2xy_S_C1001001001_c+CDY*I_ERI_Py_S_D2x_S_C1001001001_c;
  Double I_ERI_Pz_S_D2x_Py_C1001001001_c = I_ERI_Pz_S_F2xy_S_C1001001001_c+CDY*I_ERI_Pz_S_D2x_S_C1001001001_c;
  Double I_ERI_Px_S_Dxy_Py_C1001001001_c = I_ERI_Px_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Px_S_Dxy_S_C1001001001_c;
  Double I_ERI_Py_S_Dxy_Py_C1001001001_c = I_ERI_Py_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Py_S_Dxy_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxy_Py_C1001001001_c = I_ERI_Pz_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Pz_S_Dxy_S_C1001001001_c;
  Double I_ERI_Px_S_Dxz_Py_C1001001001_c = I_ERI_Px_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Px_S_Dxz_S_C1001001001_c;
  Double I_ERI_Py_S_Dxz_Py_C1001001001_c = I_ERI_Py_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Py_S_Dxz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxz_Py_C1001001001_c = I_ERI_Pz_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Pz_S_Dxz_S_C1001001001_c;
  Double I_ERI_Px_S_D2y_Py_C1001001001_c = I_ERI_Px_S_F3y_S_C1001001001_c+CDY*I_ERI_Px_S_D2y_S_C1001001001_c;
  Double I_ERI_Py_S_D2y_Py_C1001001001_c = I_ERI_Py_S_F3y_S_C1001001001_c+CDY*I_ERI_Py_S_D2y_S_C1001001001_c;
  Double I_ERI_Pz_S_D2y_Py_C1001001001_c = I_ERI_Pz_S_F3y_S_C1001001001_c+CDY*I_ERI_Pz_S_D2y_S_C1001001001_c;
  Double I_ERI_Px_S_Dyz_Py_C1001001001_c = I_ERI_Px_S_F2yz_S_C1001001001_c+CDY*I_ERI_Px_S_Dyz_S_C1001001001_c;
  Double I_ERI_Py_S_Dyz_Py_C1001001001_c = I_ERI_Py_S_F2yz_S_C1001001001_c+CDY*I_ERI_Py_S_Dyz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dyz_Py_C1001001001_c = I_ERI_Pz_S_F2yz_S_C1001001001_c+CDY*I_ERI_Pz_S_Dyz_S_C1001001001_c;
  Double I_ERI_Px_S_D2z_Py_C1001001001_c = I_ERI_Px_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Px_S_D2z_S_C1001001001_c;
  Double I_ERI_Py_S_D2z_Py_C1001001001_c = I_ERI_Py_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Py_S_D2z_S_C1001001001_c;
  Double I_ERI_Pz_S_D2z_Py_C1001001001_c = I_ERI_Pz_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Pz_S_D2z_S_C1001001001_c;
  Double I_ERI_Px_S_D2x_Pz_C1001001001_c = I_ERI_Px_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Px_S_D2x_S_C1001001001_c;
  Double I_ERI_Py_S_D2x_Pz_C1001001001_c = I_ERI_Py_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Py_S_D2x_S_C1001001001_c;
  Double I_ERI_Pz_S_D2x_Pz_C1001001001_c = I_ERI_Pz_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Pz_S_D2x_S_C1001001001_c;
  Double I_ERI_Px_S_Dxy_Pz_C1001001001_c = I_ERI_Px_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Px_S_Dxy_S_C1001001001_c;
  Double I_ERI_Py_S_Dxy_Pz_C1001001001_c = I_ERI_Py_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Py_S_Dxy_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxy_Pz_C1001001001_c = I_ERI_Pz_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Pz_S_Dxy_S_C1001001001_c;
  Double I_ERI_Px_S_Dxz_Pz_C1001001001_c = I_ERI_Px_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Px_S_Dxz_S_C1001001001_c;
  Double I_ERI_Py_S_Dxz_Pz_C1001001001_c = I_ERI_Py_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Py_S_Dxz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dxz_Pz_C1001001001_c = I_ERI_Pz_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Pz_S_Dxz_S_C1001001001_c;
  Double I_ERI_Px_S_D2y_Pz_C1001001001_c = I_ERI_Px_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Px_S_D2y_S_C1001001001_c;
  Double I_ERI_Py_S_D2y_Pz_C1001001001_c = I_ERI_Py_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Py_S_D2y_S_C1001001001_c;
  Double I_ERI_Pz_S_D2y_Pz_C1001001001_c = I_ERI_Pz_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Pz_S_D2y_S_C1001001001_c;
  Double I_ERI_Px_S_Dyz_Pz_C1001001001_c = I_ERI_Px_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Px_S_Dyz_S_C1001001001_c;
  Double I_ERI_Py_S_Dyz_Pz_C1001001001_c = I_ERI_Py_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Py_S_Dyz_S_C1001001001_c;
  Double I_ERI_Pz_S_Dyz_Pz_C1001001001_c = I_ERI_Pz_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Pz_S_Dyz_S_C1001001001_c;
  Double I_ERI_Px_S_D2z_Pz_C1001001001_c = I_ERI_Px_S_F3z_S_C1001001001_c+CDZ*I_ERI_Px_S_D2z_S_C1001001001_c;
  Double I_ERI_Py_S_D2z_Pz_C1001001001_c = I_ERI_Py_S_F3z_S_C1001001001_c+CDZ*I_ERI_Py_S_D2z_S_C1001001001_c;
  Double I_ERI_Pz_S_D2z_Pz_C1001001001_c = I_ERI_Pz_S_F3z_S_C1001001001_c+CDZ*I_ERI_Pz_S_D2z_S_C1001001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_D_P_C1001001001_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_C1001001001_c
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001001_c
   ************************************************************/
  Double I_ERI_D2x_S_D2x_Px_C1001001001_c = I_ERI_D2x_S_F3x_S_C1001001001_c+CDX*I_ERI_D2x_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2x_Px_C1001001001_c = I_ERI_Dxy_S_F3x_S_C1001001001_c+CDX*I_ERI_Dxy_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2x_Px_C1001001001_c = I_ERI_Dxz_S_F3x_S_C1001001001_c+CDX*I_ERI_Dxz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2y_S_D2x_Px_C1001001001_c = I_ERI_D2y_S_F3x_S_C1001001001_c+CDX*I_ERI_D2y_S_D2x_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2x_Px_C1001001001_c = I_ERI_Dyz_S_F3x_S_C1001001001_c+CDX*I_ERI_Dyz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2z_S_D2x_Px_C1001001001_c = I_ERI_D2z_S_F3x_S_C1001001001_c+CDX*I_ERI_D2z_S_D2x_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxy_Px_C1001001001_c = I_ERI_D2x_S_F2xy_S_C1001001001_c+CDX*I_ERI_D2x_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxy_Px_C1001001001_c = I_ERI_Dxy_S_F2xy_S_C1001001001_c+CDX*I_ERI_Dxy_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxy_Px_C1001001001_c = I_ERI_Dxz_S_F2xy_S_C1001001001_c+CDX*I_ERI_Dxz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxy_Px_C1001001001_c = I_ERI_D2y_S_F2xy_S_C1001001001_c+CDX*I_ERI_D2y_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxy_Px_C1001001001_c = I_ERI_Dyz_S_F2xy_S_C1001001001_c+CDX*I_ERI_Dyz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxy_Px_C1001001001_c = I_ERI_D2z_S_F2xy_S_C1001001001_c+CDX*I_ERI_D2z_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxz_Px_C1001001001_c = I_ERI_D2x_S_F2xz_S_C1001001001_c+CDX*I_ERI_D2x_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxz_Px_C1001001001_c = I_ERI_Dxy_S_F2xz_S_C1001001001_c+CDX*I_ERI_Dxy_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxz_Px_C1001001001_c = I_ERI_Dxz_S_F2xz_S_C1001001001_c+CDX*I_ERI_Dxz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxz_Px_C1001001001_c = I_ERI_D2y_S_F2xz_S_C1001001001_c+CDX*I_ERI_D2y_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxz_Px_C1001001001_c = I_ERI_Dyz_S_F2xz_S_C1001001001_c+CDX*I_ERI_Dyz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxz_Px_C1001001001_c = I_ERI_D2z_S_F2xz_S_C1001001001_c+CDX*I_ERI_D2z_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2y_Px_C1001001001_c = I_ERI_D2x_S_Fx2y_S_C1001001001_c+CDX*I_ERI_D2x_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2y_Px_C1001001001_c = I_ERI_Dxy_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Dxy_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2y_Px_C1001001001_c = I_ERI_Dxz_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Dxz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2y_S_D2y_Px_C1001001001_c = I_ERI_D2y_S_Fx2y_S_C1001001001_c+CDX*I_ERI_D2y_S_D2y_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2y_Px_C1001001001_c = I_ERI_Dyz_S_Fx2y_S_C1001001001_c+CDX*I_ERI_Dyz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2z_S_D2y_Px_C1001001001_c = I_ERI_D2z_S_Fx2y_S_C1001001001_c+CDX*I_ERI_D2z_S_D2y_S_C1001001001_c;
  Double I_ERI_D2x_S_Dyz_Px_C1001001001_c = I_ERI_D2x_S_Fxyz_S_C1001001001_c+CDX*I_ERI_D2x_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dyz_Px_C1001001001_c = I_ERI_Dxy_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Dxy_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dyz_Px_C1001001001_c = I_ERI_Dxz_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Dxz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dyz_Px_C1001001001_c = I_ERI_D2y_S_Fxyz_S_C1001001001_c+CDX*I_ERI_D2y_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dyz_Px_C1001001001_c = I_ERI_Dyz_S_Fxyz_S_C1001001001_c+CDX*I_ERI_Dyz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dyz_Px_C1001001001_c = I_ERI_D2z_S_Fxyz_S_C1001001001_c+CDX*I_ERI_D2z_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2z_Px_C1001001001_c = I_ERI_D2x_S_Fx2z_S_C1001001001_c+CDX*I_ERI_D2x_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2z_Px_C1001001001_c = I_ERI_Dxy_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Dxy_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2z_Px_C1001001001_c = I_ERI_Dxz_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Dxz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2y_S_D2z_Px_C1001001001_c = I_ERI_D2y_S_Fx2z_S_C1001001001_c+CDX*I_ERI_D2y_S_D2z_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2z_Px_C1001001001_c = I_ERI_Dyz_S_Fx2z_S_C1001001001_c+CDX*I_ERI_Dyz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2z_S_D2z_Px_C1001001001_c = I_ERI_D2z_S_Fx2z_S_C1001001001_c+CDX*I_ERI_D2z_S_D2z_S_C1001001001_c;
  Double I_ERI_D2x_S_D2x_Py_C1001001001_c = I_ERI_D2x_S_F2xy_S_C1001001001_c+CDY*I_ERI_D2x_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2x_Py_C1001001001_c = I_ERI_Dxy_S_F2xy_S_C1001001001_c+CDY*I_ERI_Dxy_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2x_Py_C1001001001_c = I_ERI_Dxz_S_F2xy_S_C1001001001_c+CDY*I_ERI_Dxz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2y_S_D2x_Py_C1001001001_c = I_ERI_D2y_S_F2xy_S_C1001001001_c+CDY*I_ERI_D2y_S_D2x_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2x_Py_C1001001001_c = I_ERI_Dyz_S_F2xy_S_C1001001001_c+CDY*I_ERI_Dyz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2z_S_D2x_Py_C1001001001_c = I_ERI_D2z_S_F2xy_S_C1001001001_c+CDY*I_ERI_D2z_S_D2x_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxy_Py_C1001001001_c = I_ERI_D2x_S_Fx2y_S_C1001001001_c+CDY*I_ERI_D2x_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxy_Py_C1001001001_c = I_ERI_Dxy_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Dxy_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxy_Py_C1001001001_c = I_ERI_Dxz_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Dxz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxy_Py_C1001001001_c = I_ERI_D2y_S_Fx2y_S_C1001001001_c+CDY*I_ERI_D2y_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxy_Py_C1001001001_c = I_ERI_Dyz_S_Fx2y_S_C1001001001_c+CDY*I_ERI_Dyz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxy_Py_C1001001001_c = I_ERI_D2z_S_Fx2y_S_C1001001001_c+CDY*I_ERI_D2z_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxz_Py_C1001001001_c = I_ERI_D2x_S_Fxyz_S_C1001001001_c+CDY*I_ERI_D2x_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxz_Py_C1001001001_c = I_ERI_Dxy_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Dxy_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxz_Py_C1001001001_c = I_ERI_Dxz_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Dxz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxz_Py_C1001001001_c = I_ERI_D2y_S_Fxyz_S_C1001001001_c+CDY*I_ERI_D2y_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxz_Py_C1001001001_c = I_ERI_Dyz_S_Fxyz_S_C1001001001_c+CDY*I_ERI_Dyz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxz_Py_C1001001001_c = I_ERI_D2z_S_Fxyz_S_C1001001001_c+CDY*I_ERI_D2z_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2y_Py_C1001001001_c = I_ERI_D2x_S_F3y_S_C1001001001_c+CDY*I_ERI_D2x_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2y_Py_C1001001001_c = I_ERI_Dxy_S_F3y_S_C1001001001_c+CDY*I_ERI_Dxy_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2y_Py_C1001001001_c = I_ERI_Dxz_S_F3y_S_C1001001001_c+CDY*I_ERI_Dxz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2y_S_D2y_Py_C1001001001_c = I_ERI_D2y_S_F3y_S_C1001001001_c+CDY*I_ERI_D2y_S_D2y_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2y_Py_C1001001001_c = I_ERI_Dyz_S_F3y_S_C1001001001_c+CDY*I_ERI_Dyz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2z_S_D2y_Py_C1001001001_c = I_ERI_D2z_S_F3y_S_C1001001001_c+CDY*I_ERI_D2z_S_D2y_S_C1001001001_c;
  Double I_ERI_D2x_S_Dyz_Py_C1001001001_c = I_ERI_D2x_S_F2yz_S_C1001001001_c+CDY*I_ERI_D2x_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dyz_Py_C1001001001_c = I_ERI_Dxy_S_F2yz_S_C1001001001_c+CDY*I_ERI_Dxy_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dyz_Py_C1001001001_c = I_ERI_Dxz_S_F2yz_S_C1001001001_c+CDY*I_ERI_Dxz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dyz_Py_C1001001001_c = I_ERI_D2y_S_F2yz_S_C1001001001_c+CDY*I_ERI_D2y_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dyz_Py_C1001001001_c = I_ERI_Dyz_S_F2yz_S_C1001001001_c+CDY*I_ERI_Dyz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dyz_Py_C1001001001_c = I_ERI_D2z_S_F2yz_S_C1001001001_c+CDY*I_ERI_D2z_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2z_Py_C1001001001_c = I_ERI_D2x_S_Fy2z_S_C1001001001_c+CDY*I_ERI_D2x_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2z_Py_C1001001001_c = I_ERI_Dxy_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Dxy_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2z_Py_C1001001001_c = I_ERI_Dxz_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Dxz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2y_S_D2z_Py_C1001001001_c = I_ERI_D2y_S_Fy2z_S_C1001001001_c+CDY*I_ERI_D2y_S_D2z_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2z_Py_C1001001001_c = I_ERI_Dyz_S_Fy2z_S_C1001001001_c+CDY*I_ERI_Dyz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2z_S_D2z_Py_C1001001001_c = I_ERI_D2z_S_Fy2z_S_C1001001001_c+CDY*I_ERI_D2z_S_D2z_S_C1001001001_c;
  Double I_ERI_D2x_S_D2x_Pz_C1001001001_c = I_ERI_D2x_S_F2xz_S_C1001001001_c+CDZ*I_ERI_D2x_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2x_Pz_C1001001001_c = I_ERI_Dxy_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Dxy_S_D2x_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2x_Pz_C1001001001_c = I_ERI_Dxz_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Dxz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2y_S_D2x_Pz_C1001001001_c = I_ERI_D2y_S_F2xz_S_C1001001001_c+CDZ*I_ERI_D2y_S_D2x_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2x_Pz_C1001001001_c = I_ERI_Dyz_S_F2xz_S_C1001001001_c+CDZ*I_ERI_Dyz_S_D2x_S_C1001001001_c;
  Double I_ERI_D2z_S_D2x_Pz_C1001001001_c = I_ERI_D2z_S_F2xz_S_C1001001001_c+CDZ*I_ERI_D2z_S_D2x_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxy_Pz_C1001001001_c = I_ERI_D2x_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_D2x_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxy_Pz_C1001001001_c = I_ERI_Dxy_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Dxy_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxy_Pz_C1001001001_c = I_ERI_Dxz_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Dxz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxy_Pz_C1001001001_c = I_ERI_D2y_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_D2y_S_Dxy_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxy_Pz_C1001001001_c = I_ERI_Dyz_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_Dyz_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxy_Pz_C1001001001_c = I_ERI_D2z_S_Fxyz_S_C1001001001_c+CDZ*I_ERI_D2z_S_Dxy_S_C1001001001_c;
  Double I_ERI_D2x_S_Dxz_Pz_C1001001001_c = I_ERI_D2x_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_D2x_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dxz_Pz_C1001001001_c = I_ERI_Dxy_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Dxy_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dxz_Pz_C1001001001_c = I_ERI_Dxz_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Dxz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dxz_Pz_C1001001001_c = I_ERI_D2y_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_D2y_S_Dxz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dxz_Pz_C1001001001_c = I_ERI_Dyz_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_Dyz_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dxz_Pz_C1001001001_c = I_ERI_D2z_S_Fx2z_S_C1001001001_c+CDZ*I_ERI_D2z_S_Dxz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2y_Pz_C1001001001_c = I_ERI_D2x_S_F2yz_S_C1001001001_c+CDZ*I_ERI_D2x_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2y_Pz_C1001001001_c = I_ERI_Dxy_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Dxy_S_D2y_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2y_Pz_C1001001001_c = I_ERI_Dxz_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Dxz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2y_S_D2y_Pz_C1001001001_c = I_ERI_D2y_S_F2yz_S_C1001001001_c+CDZ*I_ERI_D2y_S_D2y_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2y_Pz_C1001001001_c = I_ERI_Dyz_S_F2yz_S_C1001001001_c+CDZ*I_ERI_Dyz_S_D2y_S_C1001001001_c;
  Double I_ERI_D2z_S_D2y_Pz_C1001001001_c = I_ERI_D2z_S_F2yz_S_C1001001001_c+CDZ*I_ERI_D2z_S_D2y_S_C1001001001_c;
  Double I_ERI_D2x_S_Dyz_Pz_C1001001001_c = I_ERI_D2x_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_D2x_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxy_S_Dyz_Pz_C1001001001_c = I_ERI_Dxy_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Dxy_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dxz_S_Dyz_Pz_C1001001001_c = I_ERI_Dxz_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Dxz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2y_S_Dyz_Pz_C1001001001_c = I_ERI_D2y_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_D2y_S_Dyz_S_C1001001001_c;
  Double I_ERI_Dyz_S_Dyz_Pz_C1001001001_c = I_ERI_Dyz_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_Dyz_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2z_S_Dyz_Pz_C1001001001_c = I_ERI_D2z_S_Fy2z_S_C1001001001_c+CDZ*I_ERI_D2z_S_Dyz_S_C1001001001_c;
  Double I_ERI_D2x_S_D2z_Pz_C1001001001_c = I_ERI_D2x_S_F3z_S_C1001001001_c+CDZ*I_ERI_D2x_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxy_S_D2z_Pz_C1001001001_c = I_ERI_Dxy_S_F3z_S_C1001001001_c+CDZ*I_ERI_Dxy_S_D2z_S_C1001001001_c;
  Double I_ERI_Dxz_S_D2z_Pz_C1001001001_c = I_ERI_Dxz_S_F3z_S_C1001001001_c+CDZ*I_ERI_Dxz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2y_S_D2z_Pz_C1001001001_c = I_ERI_D2y_S_F3z_S_C1001001001_c+CDZ*I_ERI_D2y_S_D2z_S_C1001001001_c;
  Double I_ERI_Dyz_S_D2z_Pz_C1001001001_c = I_ERI_Dyz_S_F3z_S_C1001001001_c+CDZ*I_ERI_Dyz_S_D2z_S_C1001001001_c;
  Double I_ERI_D2z_S_D2z_Pz_C1001001001_c = I_ERI_D2z_S_F3z_S_C1001001001_c+CDZ*I_ERI_D2z_S_D2z_S_C1001001001_c;

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
   * shell quartet name: SQ_ERI_P_P_S_P_C1001001001
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1001001001
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1001001001
   ************************************************************/
  Double I_ERI_Px_Px_S_Px_C1001001001 = I_ERI_D2x_S_S_Px_C1001001001+ABX*I_ERI_Px_S_S_Px_C1001001001;
  Double I_ERI_Py_Px_S_Px_C1001001001 = I_ERI_Dxy_S_S_Px_C1001001001+ABX*I_ERI_Py_S_S_Px_C1001001001;
  Double I_ERI_Pz_Px_S_Px_C1001001001 = I_ERI_Dxz_S_S_Px_C1001001001+ABX*I_ERI_Pz_S_S_Px_C1001001001;
  Double I_ERI_Px_Py_S_Px_C1001001001 = I_ERI_Dxy_S_S_Px_C1001001001+ABY*I_ERI_Px_S_S_Px_C1001001001;
  Double I_ERI_Py_Py_S_Px_C1001001001 = I_ERI_D2y_S_S_Px_C1001001001+ABY*I_ERI_Py_S_S_Px_C1001001001;
  Double I_ERI_Pz_Py_S_Px_C1001001001 = I_ERI_Dyz_S_S_Px_C1001001001+ABY*I_ERI_Pz_S_S_Px_C1001001001;
  Double I_ERI_Px_Pz_S_Px_C1001001001 = I_ERI_Dxz_S_S_Px_C1001001001+ABZ*I_ERI_Px_S_S_Px_C1001001001;
  Double I_ERI_Py_Pz_S_Px_C1001001001 = I_ERI_Dyz_S_S_Px_C1001001001+ABZ*I_ERI_Py_S_S_Px_C1001001001;
  Double I_ERI_Pz_Pz_S_Px_C1001001001 = I_ERI_D2z_S_S_Px_C1001001001+ABZ*I_ERI_Pz_S_S_Px_C1001001001;
  Double I_ERI_Px_Px_S_Py_C1001001001 = I_ERI_D2x_S_S_Py_C1001001001+ABX*I_ERI_Px_S_S_Py_C1001001001;
  Double I_ERI_Py_Px_S_Py_C1001001001 = I_ERI_Dxy_S_S_Py_C1001001001+ABX*I_ERI_Py_S_S_Py_C1001001001;
  Double I_ERI_Pz_Px_S_Py_C1001001001 = I_ERI_Dxz_S_S_Py_C1001001001+ABX*I_ERI_Pz_S_S_Py_C1001001001;
  Double I_ERI_Px_Py_S_Py_C1001001001 = I_ERI_Dxy_S_S_Py_C1001001001+ABY*I_ERI_Px_S_S_Py_C1001001001;
  Double I_ERI_Py_Py_S_Py_C1001001001 = I_ERI_D2y_S_S_Py_C1001001001+ABY*I_ERI_Py_S_S_Py_C1001001001;
  Double I_ERI_Pz_Py_S_Py_C1001001001 = I_ERI_Dyz_S_S_Py_C1001001001+ABY*I_ERI_Pz_S_S_Py_C1001001001;
  Double I_ERI_Px_Pz_S_Py_C1001001001 = I_ERI_Dxz_S_S_Py_C1001001001+ABZ*I_ERI_Px_S_S_Py_C1001001001;
  Double I_ERI_Py_Pz_S_Py_C1001001001 = I_ERI_Dyz_S_S_Py_C1001001001+ABZ*I_ERI_Py_S_S_Py_C1001001001;
  Double I_ERI_Pz_Pz_S_Py_C1001001001 = I_ERI_D2z_S_S_Py_C1001001001+ABZ*I_ERI_Pz_S_S_Py_C1001001001;
  Double I_ERI_Px_Px_S_Pz_C1001001001 = I_ERI_D2x_S_S_Pz_C1001001001+ABX*I_ERI_Px_S_S_Pz_C1001001001;
  Double I_ERI_Py_Px_S_Pz_C1001001001 = I_ERI_Dxy_S_S_Pz_C1001001001+ABX*I_ERI_Py_S_S_Pz_C1001001001;
  Double I_ERI_Pz_Px_S_Pz_C1001001001 = I_ERI_Dxz_S_S_Pz_C1001001001+ABX*I_ERI_Pz_S_S_Pz_C1001001001;
  Double I_ERI_Px_Py_S_Pz_C1001001001 = I_ERI_Dxy_S_S_Pz_C1001001001+ABY*I_ERI_Px_S_S_Pz_C1001001001;
  Double I_ERI_Py_Py_S_Pz_C1001001001 = I_ERI_D2y_S_S_Pz_C1001001001+ABY*I_ERI_Py_S_S_Pz_C1001001001;
  Double I_ERI_Pz_Py_S_Pz_C1001001001 = I_ERI_Dyz_S_S_Pz_C1001001001+ABY*I_ERI_Pz_S_S_Pz_C1001001001;
  Double I_ERI_Px_Pz_S_Pz_C1001001001 = I_ERI_Dxz_S_S_Pz_C1001001001+ABZ*I_ERI_Px_S_S_Pz_C1001001001;
  Double I_ERI_Py_Pz_S_Pz_C1001001001 = I_ERI_Dyz_S_S_Pz_C1001001001+ABZ*I_ERI_Py_S_S_Pz_C1001001001;
  Double I_ERI_Pz_Pz_S_Pz_C1001001001 = I_ERI_D2z_S_S_Pz_C1001001001+ABZ*I_ERI_Pz_S_S_Pz_C1001001001;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1000_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000_a
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1000_a = I_ERI_D2x_S_S_S_C1000_a+ABX*I_ERI_Px_S_S_S_C1000_a;
  Double I_ERI_Py_Px_S_S_C1000_a = I_ERI_Dxy_S_S_S_C1000_a+ABX*I_ERI_Py_S_S_S_C1000_a;
  Double I_ERI_Pz_Px_S_S_C1000_a = I_ERI_Dxz_S_S_S_C1000_a+ABX*I_ERI_Pz_S_S_S_C1000_a;
  Double I_ERI_Px_Py_S_S_C1000_a = I_ERI_Dxy_S_S_S_C1000_a+ABY*I_ERI_Px_S_S_S_C1000_a;
  Double I_ERI_Py_Py_S_S_C1000_a = I_ERI_D2y_S_S_S_C1000_a+ABY*I_ERI_Py_S_S_S_C1000_a;
  Double I_ERI_Pz_Py_S_S_C1000_a = I_ERI_Dyz_S_S_S_C1000_a+ABY*I_ERI_Pz_S_S_S_C1000_a;
  Double I_ERI_Px_Pz_S_S_C1000_a = I_ERI_Dxz_S_S_S_C1000_a+ABZ*I_ERI_Px_S_S_S_C1000_a;
  Double I_ERI_Py_Pz_S_S_C1000_a = I_ERI_Dyz_S_S_S_C1000_a+ABZ*I_ERI_Py_S_S_S_C1000_a;
  Double I_ERI_Pz_Pz_S_S_C1000_a = I_ERI_D2z_S_S_S_C1000_a+ABZ*I_ERI_Pz_S_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_a
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001_a = I_ERI_F3x_S_S_S_C1001_a+ABX*I_ERI_D2x_S_S_S_C1001_a;
  Double I_ERI_Dxy_Px_S_S_C1001_a = I_ERI_F2xy_S_S_S_C1001_a+ABX*I_ERI_Dxy_S_S_S_C1001_a;
  Double I_ERI_Dxz_Px_S_S_C1001_a = I_ERI_F2xz_S_S_S_C1001_a+ABX*I_ERI_Dxz_S_S_S_C1001_a;
  Double I_ERI_D2y_Px_S_S_C1001_a = I_ERI_Fx2y_S_S_S_C1001_a+ABX*I_ERI_D2y_S_S_S_C1001_a;
  Double I_ERI_Dyz_Px_S_S_C1001_a = I_ERI_Fxyz_S_S_S_C1001_a+ABX*I_ERI_Dyz_S_S_S_C1001_a;
  Double I_ERI_D2z_Px_S_S_C1001_a = I_ERI_Fx2z_S_S_S_C1001_a+ABX*I_ERI_D2z_S_S_S_C1001_a;
  Double I_ERI_D2x_Py_S_S_C1001_a = I_ERI_F2xy_S_S_S_C1001_a+ABY*I_ERI_D2x_S_S_S_C1001_a;
  Double I_ERI_Dxy_Py_S_S_C1001_a = I_ERI_Fx2y_S_S_S_C1001_a+ABY*I_ERI_Dxy_S_S_S_C1001_a;
  Double I_ERI_Dxz_Py_S_S_C1001_a = I_ERI_Fxyz_S_S_S_C1001_a+ABY*I_ERI_Dxz_S_S_S_C1001_a;
  Double I_ERI_D2y_Py_S_S_C1001_a = I_ERI_F3y_S_S_S_C1001_a+ABY*I_ERI_D2y_S_S_S_C1001_a;
  Double I_ERI_Dyz_Py_S_S_C1001_a = I_ERI_F2yz_S_S_S_C1001_a+ABY*I_ERI_Dyz_S_S_S_C1001_a;
  Double I_ERI_D2z_Py_S_S_C1001_a = I_ERI_Fy2z_S_S_S_C1001_a+ABY*I_ERI_D2z_S_S_S_C1001_a;
  Double I_ERI_D2x_Pz_S_S_C1001_a = I_ERI_F2xz_S_S_S_C1001_a+ABZ*I_ERI_D2x_S_S_S_C1001_a;
  Double I_ERI_Dxy_Pz_S_S_C1001_a = I_ERI_Fxyz_S_S_S_C1001_a+ABZ*I_ERI_Dxy_S_S_S_C1001_a;
  Double I_ERI_Dxz_Pz_S_S_C1001_a = I_ERI_Fx2z_S_S_S_C1001_a+ABZ*I_ERI_Dxz_S_S_S_C1001_a;
  Double I_ERI_D2y_Pz_S_S_C1001_a = I_ERI_F2yz_S_S_S_C1001_a+ABZ*I_ERI_D2y_S_S_S_C1001_a;
  Double I_ERI_Dyz_Pz_S_S_C1001_a = I_ERI_Fy2z_S_S_S_C1001_a+ABZ*I_ERI_Dyz_S_S_S_C1001_a;
  Double I_ERI_D2z_Pz_S_S_C1001_a = I_ERI_F3z_S_S_S_C1001_a+ABZ*I_ERI_D2z_S_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001000_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001000_a
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001000_a = I_ERI_D2x_S_Px_S_C1001000_a+ABX*I_ERI_Px_S_Px_S_C1001000_a;
  Double I_ERI_Py_Px_Px_S_C1001000_a = I_ERI_Dxy_S_Px_S_C1001000_a+ABX*I_ERI_Py_S_Px_S_C1001000_a;
  Double I_ERI_Pz_Px_Px_S_C1001000_a = I_ERI_Dxz_S_Px_S_C1001000_a+ABX*I_ERI_Pz_S_Px_S_C1001000_a;
  Double I_ERI_Px_Py_Px_S_C1001000_a = I_ERI_Dxy_S_Px_S_C1001000_a+ABY*I_ERI_Px_S_Px_S_C1001000_a;
  Double I_ERI_Py_Py_Px_S_C1001000_a = I_ERI_D2y_S_Px_S_C1001000_a+ABY*I_ERI_Py_S_Px_S_C1001000_a;
  Double I_ERI_Pz_Py_Px_S_C1001000_a = I_ERI_Dyz_S_Px_S_C1001000_a+ABY*I_ERI_Pz_S_Px_S_C1001000_a;
  Double I_ERI_Px_Pz_Px_S_C1001000_a = I_ERI_Dxz_S_Px_S_C1001000_a+ABZ*I_ERI_Px_S_Px_S_C1001000_a;
  Double I_ERI_Py_Pz_Px_S_C1001000_a = I_ERI_Dyz_S_Px_S_C1001000_a+ABZ*I_ERI_Py_S_Px_S_C1001000_a;
  Double I_ERI_Pz_Pz_Px_S_C1001000_a = I_ERI_D2z_S_Px_S_C1001000_a+ABZ*I_ERI_Pz_S_Px_S_C1001000_a;
  Double I_ERI_Px_Px_Py_S_C1001000_a = I_ERI_D2x_S_Py_S_C1001000_a+ABX*I_ERI_Px_S_Py_S_C1001000_a;
  Double I_ERI_Py_Px_Py_S_C1001000_a = I_ERI_Dxy_S_Py_S_C1001000_a+ABX*I_ERI_Py_S_Py_S_C1001000_a;
  Double I_ERI_Pz_Px_Py_S_C1001000_a = I_ERI_Dxz_S_Py_S_C1001000_a+ABX*I_ERI_Pz_S_Py_S_C1001000_a;
  Double I_ERI_Px_Py_Py_S_C1001000_a = I_ERI_Dxy_S_Py_S_C1001000_a+ABY*I_ERI_Px_S_Py_S_C1001000_a;
  Double I_ERI_Py_Py_Py_S_C1001000_a = I_ERI_D2y_S_Py_S_C1001000_a+ABY*I_ERI_Py_S_Py_S_C1001000_a;
  Double I_ERI_Pz_Py_Py_S_C1001000_a = I_ERI_Dyz_S_Py_S_C1001000_a+ABY*I_ERI_Pz_S_Py_S_C1001000_a;
  Double I_ERI_Px_Pz_Py_S_C1001000_a = I_ERI_Dxz_S_Py_S_C1001000_a+ABZ*I_ERI_Px_S_Py_S_C1001000_a;
  Double I_ERI_Py_Pz_Py_S_C1001000_a = I_ERI_Dyz_S_Py_S_C1001000_a+ABZ*I_ERI_Py_S_Py_S_C1001000_a;
  Double I_ERI_Pz_Pz_Py_S_C1001000_a = I_ERI_D2z_S_Py_S_C1001000_a+ABZ*I_ERI_Pz_S_Py_S_C1001000_a;
  Double I_ERI_Px_Px_Pz_S_C1001000_a = I_ERI_D2x_S_Pz_S_C1001000_a+ABX*I_ERI_Px_S_Pz_S_C1001000_a;
  Double I_ERI_Py_Px_Pz_S_C1001000_a = I_ERI_Dxy_S_Pz_S_C1001000_a+ABX*I_ERI_Py_S_Pz_S_C1001000_a;
  Double I_ERI_Pz_Px_Pz_S_C1001000_a = I_ERI_Dxz_S_Pz_S_C1001000_a+ABX*I_ERI_Pz_S_Pz_S_C1001000_a;
  Double I_ERI_Px_Py_Pz_S_C1001000_a = I_ERI_Dxy_S_Pz_S_C1001000_a+ABY*I_ERI_Px_S_Pz_S_C1001000_a;
  Double I_ERI_Py_Py_Pz_S_C1001000_a = I_ERI_D2y_S_Pz_S_C1001000_a+ABY*I_ERI_Py_S_Pz_S_C1001000_a;
  Double I_ERI_Pz_Py_Pz_S_C1001000_a = I_ERI_Dyz_S_Pz_S_C1001000_a+ABY*I_ERI_Pz_S_Pz_S_C1001000_a;
  Double I_ERI_Px_Pz_Pz_S_C1001000_a = I_ERI_Dxz_S_Pz_S_C1001000_a+ABZ*I_ERI_Px_S_Pz_S_C1001000_a;
  Double I_ERI_Py_Pz_Pz_S_C1001000_a = I_ERI_Dyz_S_Pz_S_C1001000_a+ABZ*I_ERI_Py_S_Pz_S_C1001000_a;
  Double I_ERI_Pz_Pz_Pz_S_C1001000_a = I_ERI_D2z_S_Pz_S_C1001000_a+ABZ*I_ERI_Pz_S_Pz_S_C1001000_a;

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
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001000_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000001000_a
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001000_a
   ************************************************************/
  Double I_ERI_Px_Px_S_Px_C1000001000_a = I_ERI_D2x_S_S_Px_C1000001000_a+ABX*I_ERI_Px_S_S_Px_C1000001000_a;
  Double I_ERI_Py_Px_S_Px_C1000001000_a = I_ERI_Dxy_S_S_Px_C1000001000_a+ABX*I_ERI_Py_S_S_Px_C1000001000_a;
  Double I_ERI_Pz_Px_S_Px_C1000001000_a = I_ERI_Dxz_S_S_Px_C1000001000_a+ABX*I_ERI_Pz_S_S_Px_C1000001000_a;
  Double I_ERI_Px_Py_S_Px_C1000001000_a = I_ERI_Dxy_S_S_Px_C1000001000_a+ABY*I_ERI_Px_S_S_Px_C1000001000_a;
  Double I_ERI_Py_Py_S_Px_C1000001000_a = I_ERI_D2y_S_S_Px_C1000001000_a+ABY*I_ERI_Py_S_S_Px_C1000001000_a;
  Double I_ERI_Pz_Py_S_Px_C1000001000_a = I_ERI_Dyz_S_S_Px_C1000001000_a+ABY*I_ERI_Pz_S_S_Px_C1000001000_a;
  Double I_ERI_Px_Pz_S_Px_C1000001000_a = I_ERI_Dxz_S_S_Px_C1000001000_a+ABZ*I_ERI_Px_S_S_Px_C1000001000_a;
  Double I_ERI_Py_Pz_S_Px_C1000001000_a = I_ERI_Dyz_S_S_Px_C1000001000_a+ABZ*I_ERI_Py_S_S_Px_C1000001000_a;
  Double I_ERI_Pz_Pz_S_Px_C1000001000_a = I_ERI_D2z_S_S_Px_C1000001000_a+ABZ*I_ERI_Pz_S_S_Px_C1000001000_a;
  Double I_ERI_Px_Px_S_Py_C1000001000_a = I_ERI_D2x_S_S_Py_C1000001000_a+ABX*I_ERI_Px_S_S_Py_C1000001000_a;
  Double I_ERI_Py_Px_S_Py_C1000001000_a = I_ERI_Dxy_S_S_Py_C1000001000_a+ABX*I_ERI_Py_S_S_Py_C1000001000_a;
  Double I_ERI_Pz_Px_S_Py_C1000001000_a = I_ERI_Dxz_S_S_Py_C1000001000_a+ABX*I_ERI_Pz_S_S_Py_C1000001000_a;
  Double I_ERI_Px_Py_S_Py_C1000001000_a = I_ERI_Dxy_S_S_Py_C1000001000_a+ABY*I_ERI_Px_S_S_Py_C1000001000_a;
  Double I_ERI_Py_Py_S_Py_C1000001000_a = I_ERI_D2y_S_S_Py_C1000001000_a+ABY*I_ERI_Py_S_S_Py_C1000001000_a;
  Double I_ERI_Pz_Py_S_Py_C1000001000_a = I_ERI_Dyz_S_S_Py_C1000001000_a+ABY*I_ERI_Pz_S_S_Py_C1000001000_a;
  Double I_ERI_Px_Pz_S_Py_C1000001000_a = I_ERI_Dxz_S_S_Py_C1000001000_a+ABZ*I_ERI_Px_S_S_Py_C1000001000_a;
  Double I_ERI_Py_Pz_S_Py_C1000001000_a = I_ERI_Dyz_S_S_Py_C1000001000_a+ABZ*I_ERI_Py_S_S_Py_C1000001000_a;
  Double I_ERI_Pz_Pz_S_Py_C1000001000_a = I_ERI_D2z_S_S_Py_C1000001000_a+ABZ*I_ERI_Pz_S_S_Py_C1000001000_a;
  Double I_ERI_Px_Px_S_Pz_C1000001000_a = I_ERI_D2x_S_S_Pz_C1000001000_a+ABX*I_ERI_Px_S_S_Pz_C1000001000_a;
  Double I_ERI_Py_Px_S_Pz_C1000001000_a = I_ERI_Dxy_S_S_Pz_C1000001000_a+ABX*I_ERI_Py_S_S_Pz_C1000001000_a;
  Double I_ERI_Pz_Px_S_Pz_C1000001000_a = I_ERI_Dxz_S_S_Pz_C1000001000_a+ABX*I_ERI_Pz_S_S_Pz_C1000001000_a;
  Double I_ERI_Px_Py_S_Pz_C1000001000_a = I_ERI_Dxy_S_S_Pz_C1000001000_a+ABY*I_ERI_Px_S_S_Pz_C1000001000_a;
  Double I_ERI_Py_Py_S_Pz_C1000001000_a = I_ERI_D2y_S_S_Pz_C1000001000_a+ABY*I_ERI_Py_S_S_Pz_C1000001000_a;
  Double I_ERI_Pz_Py_S_Pz_C1000001000_a = I_ERI_Dyz_S_S_Pz_C1000001000_a+ABY*I_ERI_Pz_S_S_Pz_C1000001000_a;
  Double I_ERI_Px_Pz_S_Pz_C1000001000_a = I_ERI_Dxz_S_S_Pz_C1000001000_a+ABZ*I_ERI_Px_S_S_Pz_C1000001000_a;
  Double I_ERI_Py_Pz_S_Pz_C1000001000_a = I_ERI_Dyz_S_S_Pz_C1000001000_a+ABZ*I_ERI_Py_S_S_Pz_C1000001000_a;
  Double I_ERI_Pz_Pz_S_Pz_C1000001000_a = I_ERI_D2z_S_S_Pz_C1000001000_a+ABZ*I_ERI_Pz_S_S_Pz_C1000001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_P_C1000001001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1000001001_a
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000001001_a
   ************************************************************/
  Double I_ERI_D2x_Px_S_Px_C1000001001_a = I_ERI_F3x_S_S_Px_C1000001001_a+ABX*I_ERI_D2x_S_S_Px_C1000001001_a;
  Double I_ERI_Dxy_Px_S_Px_C1000001001_a = I_ERI_F2xy_S_S_Px_C1000001001_a+ABX*I_ERI_Dxy_S_S_Px_C1000001001_a;
  Double I_ERI_Dxz_Px_S_Px_C1000001001_a = I_ERI_F2xz_S_S_Px_C1000001001_a+ABX*I_ERI_Dxz_S_S_Px_C1000001001_a;
  Double I_ERI_D2y_Px_S_Px_C1000001001_a = I_ERI_Fx2y_S_S_Px_C1000001001_a+ABX*I_ERI_D2y_S_S_Px_C1000001001_a;
  Double I_ERI_Dyz_Px_S_Px_C1000001001_a = I_ERI_Fxyz_S_S_Px_C1000001001_a+ABX*I_ERI_Dyz_S_S_Px_C1000001001_a;
  Double I_ERI_D2z_Px_S_Px_C1000001001_a = I_ERI_Fx2z_S_S_Px_C1000001001_a+ABX*I_ERI_D2z_S_S_Px_C1000001001_a;
  Double I_ERI_D2x_Py_S_Px_C1000001001_a = I_ERI_F2xy_S_S_Px_C1000001001_a+ABY*I_ERI_D2x_S_S_Px_C1000001001_a;
  Double I_ERI_Dxy_Py_S_Px_C1000001001_a = I_ERI_Fx2y_S_S_Px_C1000001001_a+ABY*I_ERI_Dxy_S_S_Px_C1000001001_a;
  Double I_ERI_Dxz_Py_S_Px_C1000001001_a = I_ERI_Fxyz_S_S_Px_C1000001001_a+ABY*I_ERI_Dxz_S_S_Px_C1000001001_a;
  Double I_ERI_D2y_Py_S_Px_C1000001001_a = I_ERI_F3y_S_S_Px_C1000001001_a+ABY*I_ERI_D2y_S_S_Px_C1000001001_a;
  Double I_ERI_Dyz_Py_S_Px_C1000001001_a = I_ERI_F2yz_S_S_Px_C1000001001_a+ABY*I_ERI_Dyz_S_S_Px_C1000001001_a;
  Double I_ERI_D2z_Py_S_Px_C1000001001_a = I_ERI_Fy2z_S_S_Px_C1000001001_a+ABY*I_ERI_D2z_S_S_Px_C1000001001_a;
  Double I_ERI_D2x_Pz_S_Px_C1000001001_a = I_ERI_F2xz_S_S_Px_C1000001001_a+ABZ*I_ERI_D2x_S_S_Px_C1000001001_a;
  Double I_ERI_Dxy_Pz_S_Px_C1000001001_a = I_ERI_Fxyz_S_S_Px_C1000001001_a+ABZ*I_ERI_Dxy_S_S_Px_C1000001001_a;
  Double I_ERI_Dxz_Pz_S_Px_C1000001001_a = I_ERI_Fx2z_S_S_Px_C1000001001_a+ABZ*I_ERI_Dxz_S_S_Px_C1000001001_a;
  Double I_ERI_D2y_Pz_S_Px_C1000001001_a = I_ERI_F2yz_S_S_Px_C1000001001_a+ABZ*I_ERI_D2y_S_S_Px_C1000001001_a;
  Double I_ERI_Dyz_Pz_S_Px_C1000001001_a = I_ERI_Fy2z_S_S_Px_C1000001001_a+ABZ*I_ERI_Dyz_S_S_Px_C1000001001_a;
  Double I_ERI_D2z_Pz_S_Px_C1000001001_a = I_ERI_F3z_S_S_Px_C1000001001_a+ABZ*I_ERI_D2z_S_S_Px_C1000001001_a;
  Double I_ERI_D2x_Px_S_Py_C1000001001_a = I_ERI_F3x_S_S_Py_C1000001001_a+ABX*I_ERI_D2x_S_S_Py_C1000001001_a;
  Double I_ERI_Dxy_Px_S_Py_C1000001001_a = I_ERI_F2xy_S_S_Py_C1000001001_a+ABX*I_ERI_Dxy_S_S_Py_C1000001001_a;
  Double I_ERI_Dxz_Px_S_Py_C1000001001_a = I_ERI_F2xz_S_S_Py_C1000001001_a+ABX*I_ERI_Dxz_S_S_Py_C1000001001_a;
  Double I_ERI_D2y_Px_S_Py_C1000001001_a = I_ERI_Fx2y_S_S_Py_C1000001001_a+ABX*I_ERI_D2y_S_S_Py_C1000001001_a;
  Double I_ERI_Dyz_Px_S_Py_C1000001001_a = I_ERI_Fxyz_S_S_Py_C1000001001_a+ABX*I_ERI_Dyz_S_S_Py_C1000001001_a;
  Double I_ERI_D2z_Px_S_Py_C1000001001_a = I_ERI_Fx2z_S_S_Py_C1000001001_a+ABX*I_ERI_D2z_S_S_Py_C1000001001_a;
  Double I_ERI_D2x_Py_S_Py_C1000001001_a = I_ERI_F2xy_S_S_Py_C1000001001_a+ABY*I_ERI_D2x_S_S_Py_C1000001001_a;
  Double I_ERI_Dxy_Py_S_Py_C1000001001_a = I_ERI_Fx2y_S_S_Py_C1000001001_a+ABY*I_ERI_Dxy_S_S_Py_C1000001001_a;
  Double I_ERI_Dxz_Py_S_Py_C1000001001_a = I_ERI_Fxyz_S_S_Py_C1000001001_a+ABY*I_ERI_Dxz_S_S_Py_C1000001001_a;
  Double I_ERI_D2y_Py_S_Py_C1000001001_a = I_ERI_F3y_S_S_Py_C1000001001_a+ABY*I_ERI_D2y_S_S_Py_C1000001001_a;
  Double I_ERI_Dyz_Py_S_Py_C1000001001_a = I_ERI_F2yz_S_S_Py_C1000001001_a+ABY*I_ERI_Dyz_S_S_Py_C1000001001_a;
  Double I_ERI_D2z_Py_S_Py_C1000001001_a = I_ERI_Fy2z_S_S_Py_C1000001001_a+ABY*I_ERI_D2z_S_S_Py_C1000001001_a;
  Double I_ERI_D2x_Pz_S_Py_C1000001001_a = I_ERI_F2xz_S_S_Py_C1000001001_a+ABZ*I_ERI_D2x_S_S_Py_C1000001001_a;
  Double I_ERI_Dxy_Pz_S_Py_C1000001001_a = I_ERI_Fxyz_S_S_Py_C1000001001_a+ABZ*I_ERI_Dxy_S_S_Py_C1000001001_a;
  Double I_ERI_Dxz_Pz_S_Py_C1000001001_a = I_ERI_Fx2z_S_S_Py_C1000001001_a+ABZ*I_ERI_Dxz_S_S_Py_C1000001001_a;
  Double I_ERI_D2y_Pz_S_Py_C1000001001_a = I_ERI_F2yz_S_S_Py_C1000001001_a+ABZ*I_ERI_D2y_S_S_Py_C1000001001_a;
  Double I_ERI_Dyz_Pz_S_Py_C1000001001_a = I_ERI_Fy2z_S_S_Py_C1000001001_a+ABZ*I_ERI_Dyz_S_S_Py_C1000001001_a;
  Double I_ERI_D2z_Pz_S_Py_C1000001001_a = I_ERI_F3z_S_S_Py_C1000001001_a+ABZ*I_ERI_D2z_S_S_Py_C1000001001_a;
  Double I_ERI_D2x_Px_S_Pz_C1000001001_a = I_ERI_F3x_S_S_Pz_C1000001001_a+ABX*I_ERI_D2x_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxy_Px_S_Pz_C1000001001_a = I_ERI_F2xy_S_S_Pz_C1000001001_a+ABX*I_ERI_Dxy_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxz_Px_S_Pz_C1000001001_a = I_ERI_F2xz_S_S_Pz_C1000001001_a+ABX*I_ERI_Dxz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2y_Px_S_Pz_C1000001001_a = I_ERI_Fx2y_S_S_Pz_C1000001001_a+ABX*I_ERI_D2y_S_S_Pz_C1000001001_a;
  Double I_ERI_Dyz_Px_S_Pz_C1000001001_a = I_ERI_Fxyz_S_S_Pz_C1000001001_a+ABX*I_ERI_Dyz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2z_Px_S_Pz_C1000001001_a = I_ERI_Fx2z_S_S_Pz_C1000001001_a+ABX*I_ERI_D2z_S_S_Pz_C1000001001_a;
  Double I_ERI_D2x_Py_S_Pz_C1000001001_a = I_ERI_F2xy_S_S_Pz_C1000001001_a+ABY*I_ERI_D2x_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxy_Py_S_Pz_C1000001001_a = I_ERI_Fx2y_S_S_Pz_C1000001001_a+ABY*I_ERI_Dxy_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxz_Py_S_Pz_C1000001001_a = I_ERI_Fxyz_S_S_Pz_C1000001001_a+ABY*I_ERI_Dxz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2y_Py_S_Pz_C1000001001_a = I_ERI_F3y_S_S_Pz_C1000001001_a+ABY*I_ERI_D2y_S_S_Pz_C1000001001_a;
  Double I_ERI_Dyz_Py_S_Pz_C1000001001_a = I_ERI_F2yz_S_S_Pz_C1000001001_a+ABY*I_ERI_Dyz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2z_Py_S_Pz_C1000001001_a = I_ERI_Fy2z_S_S_Pz_C1000001001_a+ABY*I_ERI_D2z_S_S_Pz_C1000001001_a;
  Double I_ERI_D2x_Pz_S_Pz_C1000001001_a = I_ERI_F2xz_S_S_Pz_C1000001001_a+ABZ*I_ERI_D2x_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxy_Pz_S_Pz_C1000001001_a = I_ERI_Fxyz_S_S_Pz_C1000001001_a+ABZ*I_ERI_Dxy_S_S_Pz_C1000001001_a;
  Double I_ERI_Dxz_Pz_S_Pz_C1000001001_a = I_ERI_Fx2z_S_S_Pz_C1000001001_a+ABZ*I_ERI_Dxz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2y_Pz_S_Pz_C1000001001_a = I_ERI_F2yz_S_S_Pz_C1000001001_a+ABZ*I_ERI_D2y_S_S_Pz_C1000001001_a;
  Double I_ERI_Dyz_Pz_S_Pz_C1000001001_a = I_ERI_Fy2z_S_S_Pz_C1000001001_a+ABZ*I_ERI_Dyz_S_S_Pz_C1000001001_a;
  Double I_ERI_D2z_Pz_S_Pz_C1000001001_a = I_ERI_F3z_S_S_Pz_C1000001001_a+ABZ*I_ERI_D2z_S_S_Pz_C1000001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001000_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001000_a
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001000_a
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_C1001001000_a = I_ERI_D2x_S_Px_Px_C1001001000_a+ABX*I_ERI_Px_S_Px_Px_C1001001000_a;
  Double I_ERI_Py_Px_Px_Px_C1001001000_a = I_ERI_Dxy_S_Px_Px_C1001001000_a+ABX*I_ERI_Py_S_Px_Px_C1001001000_a;
  Double I_ERI_Pz_Px_Px_Px_C1001001000_a = I_ERI_Dxz_S_Px_Px_C1001001000_a+ABX*I_ERI_Pz_S_Px_Px_C1001001000_a;
  Double I_ERI_Px_Py_Px_Px_C1001001000_a = I_ERI_Dxy_S_Px_Px_C1001001000_a+ABY*I_ERI_Px_S_Px_Px_C1001001000_a;
  Double I_ERI_Py_Py_Px_Px_C1001001000_a = I_ERI_D2y_S_Px_Px_C1001001000_a+ABY*I_ERI_Py_S_Px_Px_C1001001000_a;
  Double I_ERI_Pz_Py_Px_Px_C1001001000_a = I_ERI_Dyz_S_Px_Px_C1001001000_a+ABY*I_ERI_Pz_S_Px_Px_C1001001000_a;
  Double I_ERI_Px_Pz_Px_Px_C1001001000_a = I_ERI_Dxz_S_Px_Px_C1001001000_a+ABZ*I_ERI_Px_S_Px_Px_C1001001000_a;
  Double I_ERI_Py_Pz_Px_Px_C1001001000_a = I_ERI_Dyz_S_Px_Px_C1001001000_a+ABZ*I_ERI_Py_S_Px_Px_C1001001000_a;
  Double I_ERI_Pz_Pz_Px_Px_C1001001000_a = I_ERI_D2z_S_Px_Px_C1001001000_a+ABZ*I_ERI_Pz_S_Px_Px_C1001001000_a;
  Double I_ERI_Px_Px_Py_Px_C1001001000_a = I_ERI_D2x_S_Py_Px_C1001001000_a+ABX*I_ERI_Px_S_Py_Px_C1001001000_a;
  Double I_ERI_Py_Px_Py_Px_C1001001000_a = I_ERI_Dxy_S_Py_Px_C1001001000_a+ABX*I_ERI_Py_S_Py_Px_C1001001000_a;
  Double I_ERI_Pz_Px_Py_Px_C1001001000_a = I_ERI_Dxz_S_Py_Px_C1001001000_a+ABX*I_ERI_Pz_S_Py_Px_C1001001000_a;
  Double I_ERI_Px_Py_Py_Px_C1001001000_a = I_ERI_Dxy_S_Py_Px_C1001001000_a+ABY*I_ERI_Px_S_Py_Px_C1001001000_a;
  Double I_ERI_Py_Py_Py_Px_C1001001000_a = I_ERI_D2y_S_Py_Px_C1001001000_a+ABY*I_ERI_Py_S_Py_Px_C1001001000_a;
  Double I_ERI_Pz_Py_Py_Px_C1001001000_a = I_ERI_Dyz_S_Py_Px_C1001001000_a+ABY*I_ERI_Pz_S_Py_Px_C1001001000_a;
  Double I_ERI_Px_Pz_Py_Px_C1001001000_a = I_ERI_Dxz_S_Py_Px_C1001001000_a+ABZ*I_ERI_Px_S_Py_Px_C1001001000_a;
  Double I_ERI_Py_Pz_Py_Px_C1001001000_a = I_ERI_Dyz_S_Py_Px_C1001001000_a+ABZ*I_ERI_Py_S_Py_Px_C1001001000_a;
  Double I_ERI_Pz_Pz_Py_Px_C1001001000_a = I_ERI_D2z_S_Py_Px_C1001001000_a+ABZ*I_ERI_Pz_S_Py_Px_C1001001000_a;
  Double I_ERI_Px_Px_Pz_Px_C1001001000_a = I_ERI_D2x_S_Pz_Px_C1001001000_a+ABX*I_ERI_Px_S_Pz_Px_C1001001000_a;
  Double I_ERI_Py_Px_Pz_Px_C1001001000_a = I_ERI_Dxy_S_Pz_Px_C1001001000_a+ABX*I_ERI_Py_S_Pz_Px_C1001001000_a;
  Double I_ERI_Pz_Px_Pz_Px_C1001001000_a = I_ERI_Dxz_S_Pz_Px_C1001001000_a+ABX*I_ERI_Pz_S_Pz_Px_C1001001000_a;
  Double I_ERI_Px_Py_Pz_Px_C1001001000_a = I_ERI_Dxy_S_Pz_Px_C1001001000_a+ABY*I_ERI_Px_S_Pz_Px_C1001001000_a;
  Double I_ERI_Py_Py_Pz_Px_C1001001000_a = I_ERI_D2y_S_Pz_Px_C1001001000_a+ABY*I_ERI_Py_S_Pz_Px_C1001001000_a;
  Double I_ERI_Pz_Py_Pz_Px_C1001001000_a = I_ERI_Dyz_S_Pz_Px_C1001001000_a+ABY*I_ERI_Pz_S_Pz_Px_C1001001000_a;
  Double I_ERI_Px_Pz_Pz_Px_C1001001000_a = I_ERI_Dxz_S_Pz_Px_C1001001000_a+ABZ*I_ERI_Px_S_Pz_Px_C1001001000_a;
  Double I_ERI_Py_Pz_Pz_Px_C1001001000_a = I_ERI_Dyz_S_Pz_Px_C1001001000_a+ABZ*I_ERI_Py_S_Pz_Px_C1001001000_a;
  Double I_ERI_Pz_Pz_Pz_Px_C1001001000_a = I_ERI_D2z_S_Pz_Px_C1001001000_a+ABZ*I_ERI_Pz_S_Pz_Px_C1001001000_a;
  Double I_ERI_Px_Px_Px_Py_C1001001000_a = I_ERI_D2x_S_Px_Py_C1001001000_a+ABX*I_ERI_Px_S_Px_Py_C1001001000_a;
  Double I_ERI_Py_Px_Px_Py_C1001001000_a = I_ERI_Dxy_S_Px_Py_C1001001000_a+ABX*I_ERI_Py_S_Px_Py_C1001001000_a;
  Double I_ERI_Pz_Px_Px_Py_C1001001000_a = I_ERI_Dxz_S_Px_Py_C1001001000_a+ABX*I_ERI_Pz_S_Px_Py_C1001001000_a;
  Double I_ERI_Px_Py_Px_Py_C1001001000_a = I_ERI_Dxy_S_Px_Py_C1001001000_a+ABY*I_ERI_Px_S_Px_Py_C1001001000_a;
  Double I_ERI_Py_Py_Px_Py_C1001001000_a = I_ERI_D2y_S_Px_Py_C1001001000_a+ABY*I_ERI_Py_S_Px_Py_C1001001000_a;
  Double I_ERI_Pz_Py_Px_Py_C1001001000_a = I_ERI_Dyz_S_Px_Py_C1001001000_a+ABY*I_ERI_Pz_S_Px_Py_C1001001000_a;
  Double I_ERI_Px_Pz_Px_Py_C1001001000_a = I_ERI_Dxz_S_Px_Py_C1001001000_a+ABZ*I_ERI_Px_S_Px_Py_C1001001000_a;
  Double I_ERI_Py_Pz_Px_Py_C1001001000_a = I_ERI_Dyz_S_Px_Py_C1001001000_a+ABZ*I_ERI_Py_S_Px_Py_C1001001000_a;
  Double I_ERI_Pz_Pz_Px_Py_C1001001000_a = I_ERI_D2z_S_Px_Py_C1001001000_a+ABZ*I_ERI_Pz_S_Px_Py_C1001001000_a;
  Double I_ERI_Px_Px_Py_Py_C1001001000_a = I_ERI_D2x_S_Py_Py_C1001001000_a+ABX*I_ERI_Px_S_Py_Py_C1001001000_a;
  Double I_ERI_Py_Px_Py_Py_C1001001000_a = I_ERI_Dxy_S_Py_Py_C1001001000_a+ABX*I_ERI_Py_S_Py_Py_C1001001000_a;
  Double I_ERI_Pz_Px_Py_Py_C1001001000_a = I_ERI_Dxz_S_Py_Py_C1001001000_a+ABX*I_ERI_Pz_S_Py_Py_C1001001000_a;
  Double I_ERI_Px_Py_Py_Py_C1001001000_a = I_ERI_Dxy_S_Py_Py_C1001001000_a+ABY*I_ERI_Px_S_Py_Py_C1001001000_a;
  Double I_ERI_Py_Py_Py_Py_C1001001000_a = I_ERI_D2y_S_Py_Py_C1001001000_a+ABY*I_ERI_Py_S_Py_Py_C1001001000_a;
  Double I_ERI_Pz_Py_Py_Py_C1001001000_a = I_ERI_Dyz_S_Py_Py_C1001001000_a+ABY*I_ERI_Pz_S_Py_Py_C1001001000_a;
  Double I_ERI_Px_Pz_Py_Py_C1001001000_a = I_ERI_Dxz_S_Py_Py_C1001001000_a+ABZ*I_ERI_Px_S_Py_Py_C1001001000_a;
  Double I_ERI_Py_Pz_Py_Py_C1001001000_a = I_ERI_Dyz_S_Py_Py_C1001001000_a+ABZ*I_ERI_Py_S_Py_Py_C1001001000_a;
  Double I_ERI_Pz_Pz_Py_Py_C1001001000_a = I_ERI_D2z_S_Py_Py_C1001001000_a+ABZ*I_ERI_Pz_S_Py_Py_C1001001000_a;
  Double I_ERI_Px_Px_Pz_Py_C1001001000_a = I_ERI_D2x_S_Pz_Py_C1001001000_a+ABX*I_ERI_Px_S_Pz_Py_C1001001000_a;
  Double I_ERI_Py_Px_Pz_Py_C1001001000_a = I_ERI_Dxy_S_Pz_Py_C1001001000_a+ABX*I_ERI_Py_S_Pz_Py_C1001001000_a;
  Double I_ERI_Pz_Px_Pz_Py_C1001001000_a = I_ERI_Dxz_S_Pz_Py_C1001001000_a+ABX*I_ERI_Pz_S_Pz_Py_C1001001000_a;
  Double I_ERI_Px_Py_Pz_Py_C1001001000_a = I_ERI_Dxy_S_Pz_Py_C1001001000_a+ABY*I_ERI_Px_S_Pz_Py_C1001001000_a;
  Double I_ERI_Py_Py_Pz_Py_C1001001000_a = I_ERI_D2y_S_Pz_Py_C1001001000_a+ABY*I_ERI_Py_S_Pz_Py_C1001001000_a;
  Double I_ERI_Pz_Py_Pz_Py_C1001001000_a = I_ERI_Dyz_S_Pz_Py_C1001001000_a+ABY*I_ERI_Pz_S_Pz_Py_C1001001000_a;
  Double I_ERI_Px_Pz_Pz_Py_C1001001000_a = I_ERI_Dxz_S_Pz_Py_C1001001000_a+ABZ*I_ERI_Px_S_Pz_Py_C1001001000_a;
  Double I_ERI_Py_Pz_Pz_Py_C1001001000_a = I_ERI_Dyz_S_Pz_Py_C1001001000_a+ABZ*I_ERI_Py_S_Pz_Py_C1001001000_a;
  Double I_ERI_Pz_Pz_Pz_Py_C1001001000_a = I_ERI_D2z_S_Pz_Py_C1001001000_a+ABZ*I_ERI_Pz_S_Pz_Py_C1001001000_a;
  Double I_ERI_Px_Px_Px_Pz_C1001001000_a = I_ERI_D2x_S_Px_Pz_C1001001000_a+ABX*I_ERI_Px_S_Px_Pz_C1001001000_a;
  Double I_ERI_Py_Px_Px_Pz_C1001001000_a = I_ERI_Dxy_S_Px_Pz_C1001001000_a+ABX*I_ERI_Py_S_Px_Pz_C1001001000_a;
  Double I_ERI_Pz_Px_Px_Pz_C1001001000_a = I_ERI_Dxz_S_Px_Pz_C1001001000_a+ABX*I_ERI_Pz_S_Px_Pz_C1001001000_a;
  Double I_ERI_Px_Py_Px_Pz_C1001001000_a = I_ERI_Dxy_S_Px_Pz_C1001001000_a+ABY*I_ERI_Px_S_Px_Pz_C1001001000_a;
  Double I_ERI_Py_Py_Px_Pz_C1001001000_a = I_ERI_D2y_S_Px_Pz_C1001001000_a+ABY*I_ERI_Py_S_Px_Pz_C1001001000_a;
  Double I_ERI_Pz_Py_Px_Pz_C1001001000_a = I_ERI_Dyz_S_Px_Pz_C1001001000_a+ABY*I_ERI_Pz_S_Px_Pz_C1001001000_a;
  Double I_ERI_Px_Pz_Px_Pz_C1001001000_a = I_ERI_Dxz_S_Px_Pz_C1001001000_a+ABZ*I_ERI_Px_S_Px_Pz_C1001001000_a;
  Double I_ERI_Py_Pz_Px_Pz_C1001001000_a = I_ERI_Dyz_S_Px_Pz_C1001001000_a+ABZ*I_ERI_Py_S_Px_Pz_C1001001000_a;
  Double I_ERI_Pz_Pz_Px_Pz_C1001001000_a = I_ERI_D2z_S_Px_Pz_C1001001000_a+ABZ*I_ERI_Pz_S_Px_Pz_C1001001000_a;
  Double I_ERI_Px_Px_Py_Pz_C1001001000_a = I_ERI_D2x_S_Py_Pz_C1001001000_a+ABX*I_ERI_Px_S_Py_Pz_C1001001000_a;
  Double I_ERI_Py_Px_Py_Pz_C1001001000_a = I_ERI_Dxy_S_Py_Pz_C1001001000_a+ABX*I_ERI_Py_S_Py_Pz_C1001001000_a;
  Double I_ERI_Pz_Px_Py_Pz_C1001001000_a = I_ERI_Dxz_S_Py_Pz_C1001001000_a+ABX*I_ERI_Pz_S_Py_Pz_C1001001000_a;
  Double I_ERI_Px_Py_Py_Pz_C1001001000_a = I_ERI_Dxy_S_Py_Pz_C1001001000_a+ABY*I_ERI_Px_S_Py_Pz_C1001001000_a;
  Double I_ERI_Py_Py_Py_Pz_C1001001000_a = I_ERI_D2y_S_Py_Pz_C1001001000_a+ABY*I_ERI_Py_S_Py_Pz_C1001001000_a;
  Double I_ERI_Pz_Py_Py_Pz_C1001001000_a = I_ERI_Dyz_S_Py_Pz_C1001001000_a+ABY*I_ERI_Pz_S_Py_Pz_C1001001000_a;
  Double I_ERI_Px_Pz_Py_Pz_C1001001000_a = I_ERI_Dxz_S_Py_Pz_C1001001000_a+ABZ*I_ERI_Px_S_Py_Pz_C1001001000_a;
  Double I_ERI_Py_Pz_Py_Pz_C1001001000_a = I_ERI_Dyz_S_Py_Pz_C1001001000_a+ABZ*I_ERI_Py_S_Py_Pz_C1001001000_a;
  Double I_ERI_Pz_Pz_Py_Pz_C1001001000_a = I_ERI_D2z_S_Py_Pz_C1001001000_a+ABZ*I_ERI_Pz_S_Py_Pz_C1001001000_a;
  Double I_ERI_Px_Px_Pz_Pz_C1001001000_a = I_ERI_D2x_S_Pz_Pz_C1001001000_a+ABX*I_ERI_Px_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Py_Px_Pz_Pz_C1001001000_a = I_ERI_Dxy_S_Pz_Pz_C1001001000_a+ABX*I_ERI_Py_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Pz_Px_Pz_Pz_C1001001000_a = I_ERI_Dxz_S_Pz_Pz_C1001001000_a+ABX*I_ERI_Pz_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Px_Py_Pz_Pz_C1001001000_a = I_ERI_Dxy_S_Pz_Pz_C1001001000_a+ABY*I_ERI_Px_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Py_Py_Pz_Pz_C1001001000_a = I_ERI_D2y_S_Pz_Pz_C1001001000_a+ABY*I_ERI_Py_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Pz_Py_Pz_Pz_C1001001000_a = I_ERI_Dyz_S_Pz_Pz_C1001001000_a+ABY*I_ERI_Pz_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Px_Pz_Pz_Pz_C1001001000_a = I_ERI_Dxz_S_Pz_Pz_C1001001000_a+ABZ*I_ERI_Px_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Py_Pz_Pz_Pz_C1001001000_a = I_ERI_Dyz_S_Pz_Pz_C1001001000_a+ABZ*I_ERI_Py_S_Pz_Pz_C1001001000_a;
  Double I_ERI_Pz_Pz_Pz_Pz_C1001001000_a = I_ERI_D2z_S_Pz_Pz_C1001001000_a+ABZ*I_ERI_Pz_S_Pz_Pz_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_P_C1001001001_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1001001001_a
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001001_a
   ************************************************************/
  Double I_ERI_D2x_Px_Px_Px_C1001001001_a = I_ERI_F3x_S_Px_Px_C1001001001_a+ABX*I_ERI_D2x_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxy_Px_Px_Px_C1001001001_a = I_ERI_F2xy_S_Px_Px_C1001001001_a+ABX*I_ERI_Dxy_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxz_Px_Px_Px_C1001001001_a = I_ERI_F2xz_S_Px_Px_C1001001001_a+ABX*I_ERI_Dxz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2y_Px_Px_Px_C1001001001_a = I_ERI_Fx2y_S_Px_Px_C1001001001_a+ABX*I_ERI_D2y_S_Px_Px_C1001001001_a;
  Double I_ERI_Dyz_Px_Px_Px_C1001001001_a = I_ERI_Fxyz_S_Px_Px_C1001001001_a+ABX*I_ERI_Dyz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2z_Px_Px_Px_C1001001001_a = I_ERI_Fx2z_S_Px_Px_C1001001001_a+ABX*I_ERI_D2z_S_Px_Px_C1001001001_a;
  Double I_ERI_D2x_Py_Px_Px_C1001001001_a = I_ERI_F2xy_S_Px_Px_C1001001001_a+ABY*I_ERI_D2x_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxy_Py_Px_Px_C1001001001_a = I_ERI_Fx2y_S_Px_Px_C1001001001_a+ABY*I_ERI_Dxy_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxz_Py_Px_Px_C1001001001_a = I_ERI_Fxyz_S_Px_Px_C1001001001_a+ABY*I_ERI_Dxz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2y_Py_Px_Px_C1001001001_a = I_ERI_F3y_S_Px_Px_C1001001001_a+ABY*I_ERI_D2y_S_Px_Px_C1001001001_a;
  Double I_ERI_Dyz_Py_Px_Px_C1001001001_a = I_ERI_F2yz_S_Px_Px_C1001001001_a+ABY*I_ERI_Dyz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2z_Py_Px_Px_C1001001001_a = I_ERI_Fy2z_S_Px_Px_C1001001001_a+ABY*I_ERI_D2z_S_Px_Px_C1001001001_a;
  Double I_ERI_D2x_Pz_Px_Px_C1001001001_a = I_ERI_F2xz_S_Px_Px_C1001001001_a+ABZ*I_ERI_D2x_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxy_Pz_Px_Px_C1001001001_a = I_ERI_Fxyz_S_Px_Px_C1001001001_a+ABZ*I_ERI_Dxy_S_Px_Px_C1001001001_a;
  Double I_ERI_Dxz_Pz_Px_Px_C1001001001_a = I_ERI_Fx2z_S_Px_Px_C1001001001_a+ABZ*I_ERI_Dxz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2y_Pz_Px_Px_C1001001001_a = I_ERI_F2yz_S_Px_Px_C1001001001_a+ABZ*I_ERI_D2y_S_Px_Px_C1001001001_a;
  Double I_ERI_Dyz_Pz_Px_Px_C1001001001_a = I_ERI_Fy2z_S_Px_Px_C1001001001_a+ABZ*I_ERI_Dyz_S_Px_Px_C1001001001_a;
  Double I_ERI_D2z_Pz_Px_Px_C1001001001_a = I_ERI_F3z_S_Px_Px_C1001001001_a+ABZ*I_ERI_D2z_S_Px_Px_C1001001001_a;
  Double I_ERI_D2x_Px_Py_Px_C1001001001_a = I_ERI_F3x_S_Py_Px_C1001001001_a+ABX*I_ERI_D2x_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxy_Px_Py_Px_C1001001001_a = I_ERI_F2xy_S_Py_Px_C1001001001_a+ABX*I_ERI_Dxy_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxz_Px_Py_Px_C1001001001_a = I_ERI_F2xz_S_Py_Px_C1001001001_a+ABX*I_ERI_Dxz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2y_Px_Py_Px_C1001001001_a = I_ERI_Fx2y_S_Py_Px_C1001001001_a+ABX*I_ERI_D2y_S_Py_Px_C1001001001_a;
  Double I_ERI_Dyz_Px_Py_Px_C1001001001_a = I_ERI_Fxyz_S_Py_Px_C1001001001_a+ABX*I_ERI_Dyz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2z_Px_Py_Px_C1001001001_a = I_ERI_Fx2z_S_Py_Px_C1001001001_a+ABX*I_ERI_D2z_S_Py_Px_C1001001001_a;
  Double I_ERI_D2x_Py_Py_Px_C1001001001_a = I_ERI_F2xy_S_Py_Px_C1001001001_a+ABY*I_ERI_D2x_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxy_Py_Py_Px_C1001001001_a = I_ERI_Fx2y_S_Py_Px_C1001001001_a+ABY*I_ERI_Dxy_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxz_Py_Py_Px_C1001001001_a = I_ERI_Fxyz_S_Py_Px_C1001001001_a+ABY*I_ERI_Dxz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2y_Py_Py_Px_C1001001001_a = I_ERI_F3y_S_Py_Px_C1001001001_a+ABY*I_ERI_D2y_S_Py_Px_C1001001001_a;
  Double I_ERI_Dyz_Py_Py_Px_C1001001001_a = I_ERI_F2yz_S_Py_Px_C1001001001_a+ABY*I_ERI_Dyz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2z_Py_Py_Px_C1001001001_a = I_ERI_Fy2z_S_Py_Px_C1001001001_a+ABY*I_ERI_D2z_S_Py_Px_C1001001001_a;
  Double I_ERI_D2x_Pz_Py_Px_C1001001001_a = I_ERI_F2xz_S_Py_Px_C1001001001_a+ABZ*I_ERI_D2x_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxy_Pz_Py_Px_C1001001001_a = I_ERI_Fxyz_S_Py_Px_C1001001001_a+ABZ*I_ERI_Dxy_S_Py_Px_C1001001001_a;
  Double I_ERI_Dxz_Pz_Py_Px_C1001001001_a = I_ERI_Fx2z_S_Py_Px_C1001001001_a+ABZ*I_ERI_Dxz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2y_Pz_Py_Px_C1001001001_a = I_ERI_F2yz_S_Py_Px_C1001001001_a+ABZ*I_ERI_D2y_S_Py_Px_C1001001001_a;
  Double I_ERI_Dyz_Pz_Py_Px_C1001001001_a = I_ERI_Fy2z_S_Py_Px_C1001001001_a+ABZ*I_ERI_Dyz_S_Py_Px_C1001001001_a;
  Double I_ERI_D2z_Pz_Py_Px_C1001001001_a = I_ERI_F3z_S_Py_Px_C1001001001_a+ABZ*I_ERI_D2z_S_Py_Px_C1001001001_a;
  Double I_ERI_D2x_Px_Pz_Px_C1001001001_a = I_ERI_F3x_S_Pz_Px_C1001001001_a+ABX*I_ERI_D2x_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxy_Px_Pz_Px_C1001001001_a = I_ERI_F2xy_S_Pz_Px_C1001001001_a+ABX*I_ERI_Dxy_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxz_Px_Pz_Px_C1001001001_a = I_ERI_F2xz_S_Pz_Px_C1001001001_a+ABX*I_ERI_Dxz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2y_Px_Pz_Px_C1001001001_a = I_ERI_Fx2y_S_Pz_Px_C1001001001_a+ABX*I_ERI_D2y_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dyz_Px_Pz_Px_C1001001001_a = I_ERI_Fxyz_S_Pz_Px_C1001001001_a+ABX*I_ERI_Dyz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2z_Px_Pz_Px_C1001001001_a = I_ERI_Fx2z_S_Pz_Px_C1001001001_a+ABX*I_ERI_D2z_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2x_Py_Pz_Px_C1001001001_a = I_ERI_F2xy_S_Pz_Px_C1001001001_a+ABY*I_ERI_D2x_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxy_Py_Pz_Px_C1001001001_a = I_ERI_Fx2y_S_Pz_Px_C1001001001_a+ABY*I_ERI_Dxy_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxz_Py_Pz_Px_C1001001001_a = I_ERI_Fxyz_S_Pz_Px_C1001001001_a+ABY*I_ERI_Dxz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2y_Py_Pz_Px_C1001001001_a = I_ERI_F3y_S_Pz_Px_C1001001001_a+ABY*I_ERI_D2y_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dyz_Py_Pz_Px_C1001001001_a = I_ERI_F2yz_S_Pz_Px_C1001001001_a+ABY*I_ERI_Dyz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2z_Py_Pz_Px_C1001001001_a = I_ERI_Fy2z_S_Pz_Px_C1001001001_a+ABY*I_ERI_D2z_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2x_Pz_Pz_Px_C1001001001_a = I_ERI_F2xz_S_Pz_Px_C1001001001_a+ABZ*I_ERI_D2x_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxy_Pz_Pz_Px_C1001001001_a = I_ERI_Fxyz_S_Pz_Px_C1001001001_a+ABZ*I_ERI_Dxy_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dxz_Pz_Pz_Px_C1001001001_a = I_ERI_Fx2z_S_Pz_Px_C1001001001_a+ABZ*I_ERI_Dxz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2y_Pz_Pz_Px_C1001001001_a = I_ERI_F2yz_S_Pz_Px_C1001001001_a+ABZ*I_ERI_D2y_S_Pz_Px_C1001001001_a;
  Double I_ERI_Dyz_Pz_Pz_Px_C1001001001_a = I_ERI_Fy2z_S_Pz_Px_C1001001001_a+ABZ*I_ERI_Dyz_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2z_Pz_Pz_Px_C1001001001_a = I_ERI_F3z_S_Pz_Px_C1001001001_a+ABZ*I_ERI_D2z_S_Pz_Px_C1001001001_a;
  Double I_ERI_D2x_Px_Px_Py_C1001001001_a = I_ERI_F3x_S_Px_Py_C1001001001_a+ABX*I_ERI_D2x_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxy_Px_Px_Py_C1001001001_a = I_ERI_F2xy_S_Px_Py_C1001001001_a+ABX*I_ERI_Dxy_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxz_Px_Px_Py_C1001001001_a = I_ERI_F2xz_S_Px_Py_C1001001001_a+ABX*I_ERI_Dxz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2y_Px_Px_Py_C1001001001_a = I_ERI_Fx2y_S_Px_Py_C1001001001_a+ABX*I_ERI_D2y_S_Px_Py_C1001001001_a;
  Double I_ERI_Dyz_Px_Px_Py_C1001001001_a = I_ERI_Fxyz_S_Px_Py_C1001001001_a+ABX*I_ERI_Dyz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2z_Px_Px_Py_C1001001001_a = I_ERI_Fx2z_S_Px_Py_C1001001001_a+ABX*I_ERI_D2z_S_Px_Py_C1001001001_a;
  Double I_ERI_D2x_Py_Px_Py_C1001001001_a = I_ERI_F2xy_S_Px_Py_C1001001001_a+ABY*I_ERI_D2x_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxy_Py_Px_Py_C1001001001_a = I_ERI_Fx2y_S_Px_Py_C1001001001_a+ABY*I_ERI_Dxy_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxz_Py_Px_Py_C1001001001_a = I_ERI_Fxyz_S_Px_Py_C1001001001_a+ABY*I_ERI_Dxz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2y_Py_Px_Py_C1001001001_a = I_ERI_F3y_S_Px_Py_C1001001001_a+ABY*I_ERI_D2y_S_Px_Py_C1001001001_a;
  Double I_ERI_Dyz_Py_Px_Py_C1001001001_a = I_ERI_F2yz_S_Px_Py_C1001001001_a+ABY*I_ERI_Dyz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2z_Py_Px_Py_C1001001001_a = I_ERI_Fy2z_S_Px_Py_C1001001001_a+ABY*I_ERI_D2z_S_Px_Py_C1001001001_a;
  Double I_ERI_D2x_Pz_Px_Py_C1001001001_a = I_ERI_F2xz_S_Px_Py_C1001001001_a+ABZ*I_ERI_D2x_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxy_Pz_Px_Py_C1001001001_a = I_ERI_Fxyz_S_Px_Py_C1001001001_a+ABZ*I_ERI_Dxy_S_Px_Py_C1001001001_a;
  Double I_ERI_Dxz_Pz_Px_Py_C1001001001_a = I_ERI_Fx2z_S_Px_Py_C1001001001_a+ABZ*I_ERI_Dxz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2y_Pz_Px_Py_C1001001001_a = I_ERI_F2yz_S_Px_Py_C1001001001_a+ABZ*I_ERI_D2y_S_Px_Py_C1001001001_a;
  Double I_ERI_Dyz_Pz_Px_Py_C1001001001_a = I_ERI_Fy2z_S_Px_Py_C1001001001_a+ABZ*I_ERI_Dyz_S_Px_Py_C1001001001_a;
  Double I_ERI_D2z_Pz_Px_Py_C1001001001_a = I_ERI_F3z_S_Px_Py_C1001001001_a+ABZ*I_ERI_D2z_S_Px_Py_C1001001001_a;
  Double I_ERI_D2x_Px_Py_Py_C1001001001_a = I_ERI_F3x_S_Py_Py_C1001001001_a+ABX*I_ERI_D2x_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxy_Px_Py_Py_C1001001001_a = I_ERI_F2xy_S_Py_Py_C1001001001_a+ABX*I_ERI_Dxy_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxz_Px_Py_Py_C1001001001_a = I_ERI_F2xz_S_Py_Py_C1001001001_a+ABX*I_ERI_Dxz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2y_Px_Py_Py_C1001001001_a = I_ERI_Fx2y_S_Py_Py_C1001001001_a+ABX*I_ERI_D2y_S_Py_Py_C1001001001_a;
  Double I_ERI_Dyz_Px_Py_Py_C1001001001_a = I_ERI_Fxyz_S_Py_Py_C1001001001_a+ABX*I_ERI_Dyz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2z_Px_Py_Py_C1001001001_a = I_ERI_Fx2z_S_Py_Py_C1001001001_a+ABX*I_ERI_D2z_S_Py_Py_C1001001001_a;
  Double I_ERI_D2x_Py_Py_Py_C1001001001_a = I_ERI_F2xy_S_Py_Py_C1001001001_a+ABY*I_ERI_D2x_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxy_Py_Py_Py_C1001001001_a = I_ERI_Fx2y_S_Py_Py_C1001001001_a+ABY*I_ERI_Dxy_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxz_Py_Py_Py_C1001001001_a = I_ERI_Fxyz_S_Py_Py_C1001001001_a+ABY*I_ERI_Dxz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2y_Py_Py_Py_C1001001001_a = I_ERI_F3y_S_Py_Py_C1001001001_a+ABY*I_ERI_D2y_S_Py_Py_C1001001001_a;
  Double I_ERI_Dyz_Py_Py_Py_C1001001001_a = I_ERI_F2yz_S_Py_Py_C1001001001_a+ABY*I_ERI_Dyz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2z_Py_Py_Py_C1001001001_a = I_ERI_Fy2z_S_Py_Py_C1001001001_a+ABY*I_ERI_D2z_S_Py_Py_C1001001001_a;
  Double I_ERI_D2x_Pz_Py_Py_C1001001001_a = I_ERI_F2xz_S_Py_Py_C1001001001_a+ABZ*I_ERI_D2x_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxy_Pz_Py_Py_C1001001001_a = I_ERI_Fxyz_S_Py_Py_C1001001001_a+ABZ*I_ERI_Dxy_S_Py_Py_C1001001001_a;
  Double I_ERI_Dxz_Pz_Py_Py_C1001001001_a = I_ERI_Fx2z_S_Py_Py_C1001001001_a+ABZ*I_ERI_Dxz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2y_Pz_Py_Py_C1001001001_a = I_ERI_F2yz_S_Py_Py_C1001001001_a+ABZ*I_ERI_D2y_S_Py_Py_C1001001001_a;
  Double I_ERI_Dyz_Pz_Py_Py_C1001001001_a = I_ERI_Fy2z_S_Py_Py_C1001001001_a+ABZ*I_ERI_Dyz_S_Py_Py_C1001001001_a;
  Double I_ERI_D2z_Pz_Py_Py_C1001001001_a = I_ERI_F3z_S_Py_Py_C1001001001_a+ABZ*I_ERI_D2z_S_Py_Py_C1001001001_a;
  Double I_ERI_D2x_Px_Pz_Py_C1001001001_a = I_ERI_F3x_S_Pz_Py_C1001001001_a+ABX*I_ERI_D2x_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxy_Px_Pz_Py_C1001001001_a = I_ERI_F2xy_S_Pz_Py_C1001001001_a+ABX*I_ERI_Dxy_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxz_Px_Pz_Py_C1001001001_a = I_ERI_F2xz_S_Pz_Py_C1001001001_a+ABX*I_ERI_Dxz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2y_Px_Pz_Py_C1001001001_a = I_ERI_Fx2y_S_Pz_Py_C1001001001_a+ABX*I_ERI_D2y_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dyz_Px_Pz_Py_C1001001001_a = I_ERI_Fxyz_S_Pz_Py_C1001001001_a+ABX*I_ERI_Dyz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2z_Px_Pz_Py_C1001001001_a = I_ERI_Fx2z_S_Pz_Py_C1001001001_a+ABX*I_ERI_D2z_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2x_Py_Pz_Py_C1001001001_a = I_ERI_F2xy_S_Pz_Py_C1001001001_a+ABY*I_ERI_D2x_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxy_Py_Pz_Py_C1001001001_a = I_ERI_Fx2y_S_Pz_Py_C1001001001_a+ABY*I_ERI_Dxy_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxz_Py_Pz_Py_C1001001001_a = I_ERI_Fxyz_S_Pz_Py_C1001001001_a+ABY*I_ERI_Dxz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2y_Py_Pz_Py_C1001001001_a = I_ERI_F3y_S_Pz_Py_C1001001001_a+ABY*I_ERI_D2y_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dyz_Py_Pz_Py_C1001001001_a = I_ERI_F2yz_S_Pz_Py_C1001001001_a+ABY*I_ERI_Dyz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2z_Py_Pz_Py_C1001001001_a = I_ERI_Fy2z_S_Pz_Py_C1001001001_a+ABY*I_ERI_D2z_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2x_Pz_Pz_Py_C1001001001_a = I_ERI_F2xz_S_Pz_Py_C1001001001_a+ABZ*I_ERI_D2x_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxy_Pz_Pz_Py_C1001001001_a = I_ERI_Fxyz_S_Pz_Py_C1001001001_a+ABZ*I_ERI_Dxy_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dxz_Pz_Pz_Py_C1001001001_a = I_ERI_Fx2z_S_Pz_Py_C1001001001_a+ABZ*I_ERI_Dxz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2y_Pz_Pz_Py_C1001001001_a = I_ERI_F2yz_S_Pz_Py_C1001001001_a+ABZ*I_ERI_D2y_S_Pz_Py_C1001001001_a;
  Double I_ERI_Dyz_Pz_Pz_Py_C1001001001_a = I_ERI_Fy2z_S_Pz_Py_C1001001001_a+ABZ*I_ERI_Dyz_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2z_Pz_Pz_Py_C1001001001_a = I_ERI_F3z_S_Pz_Py_C1001001001_a+ABZ*I_ERI_D2z_S_Pz_Py_C1001001001_a;
  Double I_ERI_D2x_Px_Px_Pz_C1001001001_a = I_ERI_F3x_S_Px_Pz_C1001001001_a+ABX*I_ERI_D2x_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxy_Px_Px_Pz_C1001001001_a = I_ERI_F2xy_S_Px_Pz_C1001001001_a+ABX*I_ERI_Dxy_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxz_Px_Px_Pz_C1001001001_a = I_ERI_F2xz_S_Px_Pz_C1001001001_a+ABX*I_ERI_Dxz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2y_Px_Px_Pz_C1001001001_a = I_ERI_Fx2y_S_Px_Pz_C1001001001_a+ABX*I_ERI_D2y_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dyz_Px_Px_Pz_C1001001001_a = I_ERI_Fxyz_S_Px_Pz_C1001001001_a+ABX*I_ERI_Dyz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2z_Px_Px_Pz_C1001001001_a = I_ERI_Fx2z_S_Px_Pz_C1001001001_a+ABX*I_ERI_D2z_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2x_Py_Px_Pz_C1001001001_a = I_ERI_F2xy_S_Px_Pz_C1001001001_a+ABY*I_ERI_D2x_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxy_Py_Px_Pz_C1001001001_a = I_ERI_Fx2y_S_Px_Pz_C1001001001_a+ABY*I_ERI_Dxy_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxz_Py_Px_Pz_C1001001001_a = I_ERI_Fxyz_S_Px_Pz_C1001001001_a+ABY*I_ERI_Dxz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2y_Py_Px_Pz_C1001001001_a = I_ERI_F3y_S_Px_Pz_C1001001001_a+ABY*I_ERI_D2y_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dyz_Py_Px_Pz_C1001001001_a = I_ERI_F2yz_S_Px_Pz_C1001001001_a+ABY*I_ERI_Dyz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2z_Py_Px_Pz_C1001001001_a = I_ERI_Fy2z_S_Px_Pz_C1001001001_a+ABY*I_ERI_D2z_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2x_Pz_Px_Pz_C1001001001_a = I_ERI_F2xz_S_Px_Pz_C1001001001_a+ABZ*I_ERI_D2x_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxy_Pz_Px_Pz_C1001001001_a = I_ERI_Fxyz_S_Px_Pz_C1001001001_a+ABZ*I_ERI_Dxy_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dxz_Pz_Px_Pz_C1001001001_a = I_ERI_Fx2z_S_Px_Pz_C1001001001_a+ABZ*I_ERI_Dxz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2y_Pz_Px_Pz_C1001001001_a = I_ERI_F2yz_S_Px_Pz_C1001001001_a+ABZ*I_ERI_D2y_S_Px_Pz_C1001001001_a;
  Double I_ERI_Dyz_Pz_Px_Pz_C1001001001_a = I_ERI_Fy2z_S_Px_Pz_C1001001001_a+ABZ*I_ERI_Dyz_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2z_Pz_Px_Pz_C1001001001_a = I_ERI_F3z_S_Px_Pz_C1001001001_a+ABZ*I_ERI_D2z_S_Px_Pz_C1001001001_a;
  Double I_ERI_D2x_Px_Py_Pz_C1001001001_a = I_ERI_F3x_S_Py_Pz_C1001001001_a+ABX*I_ERI_D2x_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxy_Px_Py_Pz_C1001001001_a = I_ERI_F2xy_S_Py_Pz_C1001001001_a+ABX*I_ERI_Dxy_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxz_Px_Py_Pz_C1001001001_a = I_ERI_F2xz_S_Py_Pz_C1001001001_a+ABX*I_ERI_Dxz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2y_Px_Py_Pz_C1001001001_a = I_ERI_Fx2y_S_Py_Pz_C1001001001_a+ABX*I_ERI_D2y_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dyz_Px_Py_Pz_C1001001001_a = I_ERI_Fxyz_S_Py_Pz_C1001001001_a+ABX*I_ERI_Dyz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2z_Px_Py_Pz_C1001001001_a = I_ERI_Fx2z_S_Py_Pz_C1001001001_a+ABX*I_ERI_D2z_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2x_Py_Py_Pz_C1001001001_a = I_ERI_F2xy_S_Py_Pz_C1001001001_a+ABY*I_ERI_D2x_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxy_Py_Py_Pz_C1001001001_a = I_ERI_Fx2y_S_Py_Pz_C1001001001_a+ABY*I_ERI_Dxy_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxz_Py_Py_Pz_C1001001001_a = I_ERI_Fxyz_S_Py_Pz_C1001001001_a+ABY*I_ERI_Dxz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2y_Py_Py_Pz_C1001001001_a = I_ERI_F3y_S_Py_Pz_C1001001001_a+ABY*I_ERI_D2y_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dyz_Py_Py_Pz_C1001001001_a = I_ERI_F2yz_S_Py_Pz_C1001001001_a+ABY*I_ERI_Dyz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2z_Py_Py_Pz_C1001001001_a = I_ERI_Fy2z_S_Py_Pz_C1001001001_a+ABY*I_ERI_D2z_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2x_Pz_Py_Pz_C1001001001_a = I_ERI_F2xz_S_Py_Pz_C1001001001_a+ABZ*I_ERI_D2x_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxy_Pz_Py_Pz_C1001001001_a = I_ERI_Fxyz_S_Py_Pz_C1001001001_a+ABZ*I_ERI_Dxy_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dxz_Pz_Py_Pz_C1001001001_a = I_ERI_Fx2z_S_Py_Pz_C1001001001_a+ABZ*I_ERI_Dxz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2y_Pz_Py_Pz_C1001001001_a = I_ERI_F2yz_S_Py_Pz_C1001001001_a+ABZ*I_ERI_D2y_S_Py_Pz_C1001001001_a;
  Double I_ERI_Dyz_Pz_Py_Pz_C1001001001_a = I_ERI_Fy2z_S_Py_Pz_C1001001001_a+ABZ*I_ERI_Dyz_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2z_Pz_Py_Pz_C1001001001_a = I_ERI_F3z_S_Py_Pz_C1001001001_a+ABZ*I_ERI_D2z_S_Py_Pz_C1001001001_a;
  Double I_ERI_D2x_Px_Pz_Pz_C1001001001_a = I_ERI_F3x_S_Pz_Pz_C1001001001_a+ABX*I_ERI_D2x_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxy_Px_Pz_Pz_C1001001001_a = I_ERI_F2xy_S_Pz_Pz_C1001001001_a+ABX*I_ERI_Dxy_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxz_Px_Pz_Pz_C1001001001_a = I_ERI_F2xz_S_Pz_Pz_C1001001001_a+ABX*I_ERI_Dxz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2y_Px_Pz_Pz_C1001001001_a = I_ERI_Fx2y_S_Pz_Pz_C1001001001_a+ABX*I_ERI_D2y_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dyz_Px_Pz_Pz_C1001001001_a = I_ERI_Fxyz_S_Pz_Pz_C1001001001_a+ABX*I_ERI_Dyz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2z_Px_Pz_Pz_C1001001001_a = I_ERI_Fx2z_S_Pz_Pz_C1001001001_a+ABX*I_ERI_D2z_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2x_Py_Pz_Pz_C1001001001_a = I_ERI_F2xy_S_Pz_Pz_C1001001001_a+ABY*I_ERI_D2x_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxy_Py_Pz_Pz_C1001001001_a = I_ERI_Fx2y_S_Pz_Pz_C1001001001_a+ABY*I_ERI_Dxy_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxz_Py_Pz_Pz_C1001001001_a = I_ERI_Fxyz_S_Pz_Pz_C1001001001_a+ABY*I_ERI_Dxz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2y_Py_Pz_Pz_C1001001001_a = I_ERI_F3y_S_Pz_Pz_C1001001001_a+ABY*I_ERI_D2y_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dyz_Py_Pz_Pz_C1001001001_a = I_ERI_F2yz_S_Pz_Pz_C1001001001_a+ABY*I_ERI_Dyz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2z_Py_Pz_Pz_C1001001001_a = I_ERI_Fy2z_S_Pz_Pz_C1001001001_a+ABY*I_ERI_D2z_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2x_Pz_Pz_Pz_C1001001001_a = I_ERI_F2xz_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_D2x_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxy_Pz_Pz_Pz_C1001001001_a = I_ERI_Fxyz_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_Dxy_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dxz_Pz_Pz_Pz_C1001001001_a = I_ERI_Fx2z_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_Dxz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2y_Pz_Pz_Pz_C1001001001_a = I_ERI_F2yz_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_D2y_S_Pz_Pz_C1001001001_a;
  Double I_ERI_Dyz_Pz_Pz_Pz_C1001001001_a = I_ERI_Fy2z_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_Dyz_S_Pz_Pz_C1001001001_a;
  Double I_ERI_D2z_Pz_Pz_Pz_C1001001001_a = I_ERI_F3z_S_Pz_Pz_C1001001001_a+ABZ*I_ERI_D2z_S_Pz_Pz_C1001001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1_b
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1_b = I_ERI_D2x_S_S_S_C1_b+ABX*I_ERI_Px_S_S_S_C1_b;
  Double I_ERI_Py_Px_S_S_C1_b = I_ERI_Dxy_S_S_S_C1_b+ABX*I_ERI_Py_S_S_S_C1_b;
  Double I_ERI_Pz_Px_S_S_C1_b = I_ERI_Dxz_S_S_S_C1_b+ABX*I_ERI_Pz_S_S_S_C1_b;
  Double I_ERI_Px_Py_S_S_C1_b = I_ERI_Dxy_S_S_S_C1_b+ABY*I_ERI_Px_S_S_S_C1_b;
  Double I_ERI_Py_Py_S_S_C1_b = I_ERI_D2y_S_S_S_C1_b+ABY*I_ERI_Py_S_S_S_C1_b;
  Double I_ERI_Pz_Py_S_S_C1_b = I_ERI_Dyz_S_S_S_C1_b+ABY*I_ERI_Pz_S_S_S_C1_b;
  Double I_ERI_Px_Pz_S_S_C1_b = I_ERI_Dxz_S_S_S_C1_b+ABZ*I_ERI_Px_S_S_S_C1_b;
  Double I_ERI_Py_Pz_S_S_C1_b = I_ERI_Dyz_S_S_S_C1_b+ABZ*I_ERI_Py_S_S_S_C1_b;
  Double I_ERI_Pz_Pz_S_S_C1_b = I_ERI_D2z_S_S_S_C1_b+ABZ*I_ERI_Pz_S_S_S_C1_b;

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
   * shell quartet name: SQ_ERI_D_P_S_S_C1001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 4 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001_b
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001_b = I_ERI_F3x_S_S_S_C1001_b+ABX*I_ERI_D2x_S_S_S_C1001_b;
  Double I_ERI_Dxy_Px_S_S_C1001_b = I_ERI_F2xy_S_S_S_C1001_b+ABX*I_ERI_Dxy_S_S_S_C1001_b;
  Double I_ERI_Dxz_Px_S_S_C1001_b = I_ERI_F2xz_S_S_S_C1001_b+ABX*I_ERI_Dxz_S_S_S_C1001_b;
  Double I_ERI_D2y_Px_S_S_C1001_b = I_ERI_Fx2y_S_S_S_C1001_b+ABX*I_ERI_D2y_S_S_S_C1001_b;
  Double I_ERI_Dyz_Px_S_S_C1001_b = I_ERI_Fxyz_S_S_S_C1001_b+ABX*I_ERI_Dyz_S_S_S_C1001_b;
  Double I_ERI_D2z_Px_S_S_C1001_b = I_ERI_Fx2z_S_S_S_C1001_b+ABX*I_ERI_D2z_S_S_S_C1001_b;
  Double I_ERI_Dxy_Py_S_S_C1001_b = I_ERI_Fx2y_S_S_S_C1001_b+ABY*I_ERI_Dxy_S_S_S_C1001_b;
  Double I_ERI_Dxz_Py_S_S_C1001_b = I_ERI_Fxyz_S_S_S_C1001_b+ABY*I_ERI_Dxz_S_S_S_C1001_b;
  Double I_ERI_D2y_Py_S_S_C1001_b = I_ERI_F3y_S_S_S_C1001_b+ABY*I_ERI_D2y_S_S_S_C1001_b;
  Double I_ERI_Dyz_Py_S_S_C1001_b = I_ERI_F2yz_S_S_S_C1001_b+ABY*I_ERI_Dyz_S_S_S_C1001_b;
  Double I_ERI_D2z_Py_S_S_C1001_b = I_ERI_Fy2z_S_S_S_C1001_b+ABY*I_ERI_D2z_S_S_S_C1001_b;
  Double I_ERI_Dxz_Pz_S_S_C1001_b = I_ERI_Fx2z_S_S_S_C1001_b+ABZ*I_ERI_Dxz_S_S_S_C1001_b;
  Double I_ERI_Dyz_Pz_S_S_C1001_b = I_ERI_Fy2z_S_S_S_C1001_b+ABZ*I_ERI_Dyz_S_S_S_C1001_b;
  Double I_ERI_D2z_Pz_S_S_C1001_b = I_ERI_F3z_S_S_S_C1001_b+ABZ*I_ERI_D2z_S_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S_C1001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001_b
   ************************************************************/
  Double I_ERI_Px_D2x_S_S_C1001_b = I_ERI_D2x_Px_S_S_C1001_b+ABX*I_ERI_Px_Px_S_S_C1001_b;
  Double I_ERI_Py_D2x_S_S_C1001_b = I_ERI_Dxy_Px_S_S_C1001_b+ABX*I_ERI_Py_Px_S_S_C1001_b;
  Double I_ERI_Pz_D2x_S_S_C1001_b = I_ERI_Dxz_Px_S_S_C1001_b+ABX*I_ERI_Pz_Px_S_S_C1001_b;
  Double I_ERI_Px_Dxy_S_S_C1001_b = I_ERI_Dxy_Px_S_S_C1001_b+ABY*I_ERI_Px_Px_S_S_C1001_b;
  Double I_ERI_Py_Dxy_S_S_C1001_b = I_ERI_D2y_Px_S_S_C1001_b+ABY*I_ERI_Py_Px_S_S_C1001_b;
  Double I_ERI_Pz_Dxy_S_S_C1001_b = I_ERI_Dyz_Px_S_S_C1001_b+ABY*I_ERI_Pz_Px_S_S_C1001_b;
  Double I_ERI_Px_Dxz_S_S_C1001_b = I_ERI_Dxz_Px_S_S_C1001_b+ABZ*I_ERI_Px_Px_S_S_C1001_b;
  Double I_ERI_Py_Dxz_S_S_C1001_b = I_ERI_Dyz_Px_S_S_C1001_b+ABZ*I_ERI_Py_Px_S_S_C1001_b;
  Double I_ERI_Pz_Dxz_S_S_C1001_b = I_ERI_D2z_Px_S_S_C1001_b+ABZ*I_ERI_Pz_Px_S_S_C1001_b;
  Double I_ERI_Px_D2y_S_S_C1001_b = I_ERI_Dxy_Py_S_S_C1001_b+ABY*I_ERI_Px_Py_S_S_C1001_b;
  Double I_ERI_Py_D2y_S_S_C1001_b = I_ERI_D2y_Py_S_S_C1001_b+ABY*I_ERI_Py_Py_S_S_C1001_b;
  Double I_ERI_Pz_D2y_S_S_C1001_b = I_ERI_Dyz_Py_S_S_C1001_b+ABY*I_ERI_Pz_Py_S_S_C1001_b;
  Double I_ERI_Px_Dyz_S_S_C1001_b = I_ERI_Dxz_Py_S_S_C1001_b+ABZ*I_ERI_Px_Py_S_S_C1001_b;
  Double I_ERI_Py_Dyz_S_S_C1001_b = I_ERI_Dyz_Py_S_S_C1001_b+ABZ*I_ERI_Py_Py_S_S_C1001_b;
  Double I_ERI_Pz_Dyz_S_S_C1001_b = I_ERI_D2z_Py_S_S_C1001_b+ABZ*I_ERI_Pz_Py_S_S_C1001_b;
  Double I_ERI_Px_D2z_S_S_C1001_b = I_ERI_Dxz_Pz_S_S_C1001_b+ABZ*I_ERI_Px_Pz_S_S_C1001_b;
  Double I_ERI_Py_D2z_S_S_C1001_b = I_ERI_Dyz_Pz_S_S_C1001_b+ABZ*I_ERI_Py_Pz_S_S_C1001_b;
  Double I_ERI_Pz_D2z_S_S_C1001_b = I_ERI_D2z_Pz_S_S_C1001_b+ABZ*I_ERI_Pz_Pz_S_S_C1001_b;

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
   * shell quartet name: SQ_ERI_P_P_S_P_C1000000001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000000001_b
   ************************************************************/
  Double I_ERI_Px_Px_S_Px_C1000000001_b = I_ERI_D2x_S_S_Px_C1000000001_b+ABX*I_ERI_Px_S_S_Px_C1000000001_b;
  Double I_ERI_Py_Px_S_Px_C1000000001_b = I_ERI_Dxy_S_S_Px_C1000000001_b+ABX*I_ERI_Py_S_S_Px_C1000000001_b;
  Double I_ERI_Pz_Px_S_Px_C1000000001_b = I_ERI_Dxz_S_S_Px_C1000000001_b+ABX*I_ERI_Pz_S_S_Px_C1000000001_b;
  Double I_ERI_Px_Py_S_Px_C1000000001_b = I_ERI_Dxy_S_S_Px_C1000000001_b+ABY*I_ERI_Px_S_S_Px_C1000000001_b;
  Double I_ERI_Py_Py_S_Px_C1000000001_b = I_ERI_D2y_S_S_Px_C1000000001_b+ABY*I_ERI_Py_S_S_Px_C1000000001_b;
  Double I_ERI_Pz_Py_S_Px_C1000000001_b = I_ERI_Dyz_S_S_Px_C1000000001_b+ABY*I_ERI_Pz_S_S_Px_C1000000001_b;
  Double I_ERI_Px_Pz_S_Px_C1000000001_b = I_ERI_Dxz_S_S_Px_C1000000001_b+ABZ*I_ERI_Px_S_S_Px_C1000000001_b;
  Double I_ERI_Py_Pz_S_Px_C1000000001_b = I_ERI_Dyz_S_S_Px_C1000000001_b+ABZ*I_ERI_Py_S_S_Px_C1000000001_b;
  Double I_ERI_Pz_Pz_S_Px_C1000000001_b = I_ERI_D2z_S_S_Px_C1000000001_b+ABZ*I_ERI_Pz_S_S_Px_C1000000001_b;
  Double I_ERI_Px_Px_S_Py_C1000000001_b = I_ERI_D2x_S_S_Py_C1000000001_b+ABX*I_ERI_Px_S_S_Py_C1000000001_b;
  Double I_ERI_Py_Px_S_Py_C1000000001_b = I_ERI_Dxy_S_S_Py_C1000000001_b+ABX*I_ERI_Py_S_S_Py_C1000000001_b;
  Double I_ERI_Pz_Px_S_Py_C1000000001_b = I_ERI_Dxz_S_S_Py_C1000000001_b+ABX*I_ERI_Pz_S_S_Py_C1000000001_b;
  Double I_ERI_Px_Py_S_Py_C1000000001_b = I_ERI_Dxy_S_S_Py_C1000000001_b+ABY*I_ERI_Px_S_S_Py_C1000000001_b;
  Double I_ERI_Py_Py_S_Py_C1000000001_b = I_ERI_D2y_S_S_Py_C1000000001_b+ABY*I_ERI_Py_S_S_Py_C1000000001_b;
  Double I_ERI_Pz_Py_S_Py_C1000000001_b = I_ERI_Dyz_S_S_Py_C1000000001_b+ABY*I_ERI_Pz_S_S_Py_C1000000001_b;
  Double I_ERI_Px_Pz_S_Py_C1000000001_b = I_ERI_Dxz_S_S_Py_C1000000001_b+ABZ*I_ERI_Px_S_S_Py_C1000000001_b;
  Double I_ERI_Py_Pz_S_Py_C1000000001_b = I_ERI_Dyz_S_S_Py_C1000000001_b+ABZ*I_ERI_Py_S_S_Py_C1000000001_b;
  Double I_ERI_Pz_Pz_S_Py_C1000000001_b = I_ERI_D2z_S_S_Py_C1000000001_b+ABZ*I_ERI_Pz_S_S_Py_C1000000001_b;
  Double I_ERI_Px_Px_S_Pz_C1000000001_b = I_ERI_D2x_S_S_Pz_C1000000001_b+ABX*I_ERI_Px_S_S_Pz_C1000000001_b;
  Double I_ERI_Py_Px_S_Pz_C1000000001_b = I_ERI_Dxy_S_S_Pz_C1000000001_b+ABX*I_ERI_Py_S_S_Pz_C1000000001_b;
  Double I_ERI_Pz_Px_S_Pz_C1000000001_b = I_ERI_Dxz_S_S_Pz_C1000000001_b+ABX*I_ERI_Pz_S_S_Pz_C1000000001_b;
  Double I_ERI_Px_Py_S_Pz_C1000000001_b = I_ERI_Dxy_S_S_Pz_C1000000001_b+ABY*I_ERI_Px_S_S_Pz_C1000000001_b;
  Double I_ERI_Py_Py_S_Pz_C1000000001_b = I_ERI_D2y_S_S_Pz_C1000000001_b+ABY*I_ERI_Py_S_S_Pz_C1000000001_b;
  Double I_ERI_Pz_Py_S_Pz_C1000000001_b = I_ERI_Dyz_S_S_Pz_C1000000001_b+ABY*I_ERI_Pz_S_S_Pz_C1000000001_b;
  Double I_ERI_Px_Pz_S_Pz_C1000000001_b = I_ERI_Dxz_S_S_Pz_C1000000001_b+ABZ*I_ERI_Px_S_S_Pz_C1000000001_b;
  Double I_ERI_Py_Pz_S_Pz_C1000000001_b = I_ERI_Dyz_S_S_Pz_C1000000001_b+ABZ*I_ERI_Py_S_S_Pz_C1000000001_b;
  Double I_ERI_Pz_Pz_S_Pz_C1000000001_b = I_ERI_D2z_S_S_Pz_C1000000001_b+ABZ*I_ERI_Pz_S_S_Pz_C1000000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001001_b
   ************************************************************/
  Double I_ERI_Px_Px_S_Px_C1000001001_b = I_ERI_D2x_S_S_Px_C1000001001_b+ABX*I_ERI_Px_S_S_Px_C1000001001_b;
  Double I_ERI_Py_Px_S_Px_C1000001001_b = I_ERI_Dxy_S_S_Px_C1000001001_b+ABX*I_ERI_Py_S_S_Px_C1000001001_b;
  Double I_ERI_Pz_Px_S_Px_C1000001001_b = I_ERI_Dxz_S_S_Px_C1000001001_b+ABX*I_ERI_Pz_S_S_Px_C1000001001_b;
  Double I_ERI_Px_Py_S_Px_C1000001001_b = I_ERI_Dxy_S_S_Px_C1000001001_b+ABY*I_ERI_Px_S_S_Px_C1000001001_b;
  Double I_ERI_Py_Py_S_Px_C1000001001_b = I_ERI_D2y_S_S_Px_C1000001001_b+ABY*I_ERI_Py_S_S_Px_C1000001001_b;
  Double I_ERI_Pz_Py_S_Px_C1000001001_b = I_ERI_Dyz_S_S_Px_C1000001001_b+ABY*I_ERI_Pz_S_S_Px_C1000001001_b;
  Double I_ERI_Px_Pz_S_Px_C1000001001_b = I_ERI_Dxz_S_S_Px_C1000001001_b+ABZ*I_ERI_Px_S_S_Px_C1000001001_b;
  Double I_ERI_Py_Pz_S_Px_C1000001001_b = I_ERI_Dyz_S_S_Px_C1000001001_b+ABZ*I_ERI_Py_S_S_Px_C1000001001_b;
  Double I_ERI_Pz_Pz_S_Px_C1000001001_b = I_ERI_D2z_S_S_Px_C1000001001_b+ABZ*I_ERI_Pz_S_S_Px_C1000001001_b;
  Double I_ERI_Px_Px_S_Py_C1000001001_b = I_ERI_D2x_S_S_Py_C1000001001_b+ABX*I_ERI_Px_S_S_Py_C1000001001_b;
  Double I_ERI_Py_Px_S_Py_C1000001001_b = I_ERI_Dxy_S_S_Py_C1000001001_b+ABX*I_ERI_Py_S_S_Py_C1000001001_b;
  Double I_ERI_Pz_Px_S_Py_C1000001001_b = I_ERI_Dxz_S_S_Py_C1000001001_b+ABX*I_ERI_Pz_S_S_Py_C1000001001_b;
  Double I_ERI_Px_Py_S_Py_C1000001001_b = I_ERI_Dxy_S_S_Py_C1000001001_b+ABY*I_ERI_Px_S_S_Py_C1000001001_b;
  Double I_ERI_Py_Py_S_Py_C1000001001_b = I_ERI_D2y_S_S_Py_C1000001001_b+ABY*I_ERI_Py_S_S_Py_C1000001001_b;
  Double I_ERI_Pz_Py_S_Py_C1000001001_b = I_ERI_Dyz_S_S_Py_C1000001001_b+ABY*I_ERI_Pz_S_S_Py_C1000001001_b;
  Double I_ERI_Px_Pz_S_Py_C1000001001_b = I_ERI_Dxz_S_S_Py_C1000001001_b+ABZ*I_ERI_Px_S_S_Py_C1000001001_b;
  Double I_ERI_Py_Pz_S_Py_C1000001001_b = I_ERI_Dyz_S_S_Py_C1000001001_b+ABZ*I_ERI_Py_S_S_Py_C1000001001_b;
  Double I_ERI_Pz_Pz_S_Py_C1000001001_b = I_ERI_D2z_S_S_Py_C1000001001_b+ABZ*I_ERI_Pz_S_S_Py_C1000001001_b;
  Double I_ERI_Px_Px_S_Pz_C1000001001_b = I_ERI_D2x_S_S_Pz_C1000001001_b+ABX*I_ERI_Px_S_S_Pz_C1000001001_b;
  Double I_ERI_Py_Px_S_Pz_C1000001001_b = I_ERI_Dxy_S_S_Pz_C1000001001_b+ABX*I_ERI_Py_S_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Px_S_Pz_C1000001001_b = I_ERI_Dxz_S_S_Pz_C1000001001_b+ABX*I_ERI_Pz_S_S_Pz_C1000001001_b;
  Double I_ERI_Px_Py_S_Pz_C1000001001_b = I_ERI_Dxy_S_S_Pz_C1000001001_b+ABY*I_ERI_Px_S_S_Pz_C1000001001_b;
  Double I_ERI_Py_Py_S_Pz_C1000001001_b = I_ERI_D2y_S_S_Pz_C1000001001_b+ABY*I_ERI_Py_S_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Py_S_Pz_C1000001001_b = I_ERI_Dyz_S_S_Pz_C1000001001_b+ABY*I_ERI_Pz_S_S_Pz_C1000001001_b;
  Double I_ERI_Px_Pz_S_Pz_C1000001001_b = I_ERI_Dxz_S_S_Pz_C1000001001_b+ABZ*I_ERI_Px_S_S_Pz_C1000001001_b;
  Double I_ERI_Py_Pz_S_Pz_C1000001001_b = I_ERI_Dyz_S_S_Pz_C1000001001_b+ABZ*I_ERI_Py_S_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Pz_S_Pz_C1000001001_b = I_ERI_D2z_S_S_Pz_C1000001001_b+ABZ*I_ERI_Pz_S_S_Pz_C1000001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_P_C1000001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000001001_b
   ************************************************************/
  Double I_ERI_D2x_Px_S_Px_C1000001001_b = I_ERI_F3x_S_S_Px_C1000001001_b+ABX*I_ERI_D2x_S_S_Px_C1000001001_b;
  Double I_ERI_Dxy_Px_S_Px_C1000001001_b = I_ERI_F2xy_S_S_Px_C1000001001_b+ABX*I_ERI_Dxy_S_S_Px_C1000001001_b;
  Double I_ERI_Dxz_Px_S_Px_C1000001001_b = I_ERI_F2xz_S_S_Px_C1000001001_b+ABX*I_ERI_Dxz_S_S_Px_C1000001001_b;
  Double I_ERI_D2y_Px_S_Px_C1000001001_b = I_ERI_Fx2y_S_S_Px_C1000001001_b+ABX*I_ERI_D2y_S_S_Px_C1000001001_b;
  Double I_ERI_Dyz_Px_S_Px_C1000001001_b = I_ERI_Fxyz_S_S_Px_C1000001001_b+ABX*I_ERI_Dyz_S_S_Px_C1000001001_b;
  Double I_ERI_D2z_Px_S_Px_C1000001001_b = I_ERI_Fx2z_S_S_Px_C1000001001_b+ABX*I_ERI_D2z_S_S_Px_C1000001001_b;
  Double I_ERI_Dxy_Py_S_Px_C1000001001_b = I_ERI_Fx2y_S_S_Px_C1000001001_b+ABY*I_ERI_Dxy_S_S_Px_C1000001001_b;
  Double I_ERI_Dxz_Py_S_Px_C1000001001_b = I_ERI_Fxyz_S_S_Px_C1000001001_b+ABY*I_ERI_Dxz_S_S_Px_C1000001001_b;
  Double I_ERI_D2y_Py_S_Px_C1000001001_b = I_ERI_F3y_S_S_Px_C1000001001_b+ABY*I_ERI_D2y_S_S_Px_C1000001001_b;
  Double I_ERI_Dyz_Py_S_Px_C1000001001_b = I_ERI_F2yz_S_S_Px_C1000001001_b+ABY*I_ERI_Dyz_S_S_Px_C1000001001_b;
  Double I_ERI_D2z_Py_S_Px_C1000001001_b = I_ERI_Fy2z_S_S_Px_C1000001001_b+ABY*I_ERI_D2z_S_S_Px_C1000001001_b;
  Double I_ERI_Dxz_Pz_S_Px_C1000001001_b = I_ERI_Fx2z_S_S_Px_C1000001001_b+ABZ*I_ERI_Dxz_S_S_Px_C1000001001_b;
  Double I_ERI_Dyz_Pz_S_Px_C1000001001_b = I_ERI_Fy2z_S_S_Px_C1000001001_b+ABZ*I_ERI_Dyz_S_S_Px_C1000001001_b;
  Double I_ERI_D2z_Pz_S_Px_C1000001001_b = I_ERI_F3z_S_S_Px_C1000001001_b+ABZ*I_ERI_D2z_S_S_Px_C1000001001_b;
  Double I_ERI_D2x_Px_S_Py_C1000001001_b = I_ERI_F3x_S_S_Py_C1000001001_b+ABX*I_ERI_D2x_S_S_Py_C1000001001_b;
  Double I_ERI_Dxy_Px_S_Py_C1000001001_b = I_ERI_F2xy_S_S_Py_C1000001001_b+ABX*I_ERI_Dxy_S_S_Py_C1000001001_b;
  Double I_ERI_Dxz_Px_S_Py_C1000001001_b = I_ERI_F2xz_S_S_Py_C1000001001_b+ABX*I_ERI_Dxz_S_S_Py_C1000001001_b;
  Double I_ERI_D2y_Px_S_Py_C1000001001_b = I_ERI_Fx2y_S_S_Py_C1000001001_b+ABX*I_ERI_D2y_S_S_Py_C1000001001_b;
  Double I_ERI_Dyz_Px_S_Py_C1000001001_b = I_ERI_Fxyz_S_S_Py_C1000001001_b+ABX*I_ERI_Dyz_S_S_Py_C1000001001_b;
  Double I_ERI_D2z_Px_S_Py_C1000001001_b = I_ERI_Fx2z_S_S_Py_C1000001001_b+ABX*I_ERI_D2z_S_S_Py_C1000001001_b;
  Double I_ERI_Dxy_Py_S_Py_C1000001001_b = I_ERI_Fx2y_S_S_Py_C1000001001_b+ABY*I_ERI_Dxy_S_S_Py_C1000001001_b;
  Double I_ERI_Dxz_Py_S_Py_C1000001001_b = I_ERI_Fxyz_S_S_Py_C1000001001_b+ABY*I_ERI_Dxz_S_S_Py_C1000001001_b;
  Double I_ERI_D2y_Py_S_Py_C1000001001_b = I_ERI_F3y_S_S_Py_C1000001001_b+ABY*I_ERI_D2y_S_S_Py_C1000001001_b;
  Double I_ERI_Dyz_Py_S_Py_C1000001001_b = I_ERI_F2yz_S_S_Py_C1000001001_b+ABY*I_ERI_Dyz_S_S_Py_C1000001001_b;
  Double I_ERI_D2z_Py_S_Py_C1000001001_b = I_ERI_Fy2z_S_S_Py_C1000001001_b+ABY*I_ERI_D2z_S_S_Py_C1000001001_b;
  Double I_ERI_Dxz_Pz_S_Py_C1000001001_b = I_ERI_Fx2z_S_S_Py_C1000001001_b+ABZ*I_ERI_Dxz_S_S_Py_C1000001001_b;
  Double I_ERI_Dyz_Pz_S_Py_C1000001001_b = I_ERI_Fy2z_S_S_Py_C1000001001_b+ABZ*I_ERI_Dyz_S_S_Py_C1000001001_b;
  Double I_ERI_D2z_Pz_S_Py_C1000001001_b = I_ERI_F3z_S_S_Py_C1000001001_b+ABZ*I_ERI_D2z_S_S_Py_C1000001001_b;
  Double I_ERI_D2x_Px_S_Pz_C1000001001_b = I_ERI_F3x_S_S_Pz_C1000001001_b+ABX*I_ERI_D2x_S_S_Pz_C1000001001_b;
  Double I_ERI_Dxy_Px_S_Pz_C1000001001_b = I_ERI_F2xy_S_S_Pz_C1000001001_b+ABX*I_ERI_Dxy_S_S_Pz_C1000001001_b;
  Double I_ERI_Dxz_Px_S_Pz_C1000001001_b = I_ERI_F2xz_S_S_Pz_C1000001001_b+ABX*I_ERI_Dxz_S_S_Pz_C1000001001_b;
  Double I_ERI_D2y_Px_S_Pz_C1000001001_b = I_ERI_Fx2y_S_S_Pz_C1000001001_b+ABX*I_ERI_D2y_S_S_Pz_C1000001001_b;
  Double I_ERI_Dyz_Px_S_Pz_C1000001001_b = I_ERI_Fxyz_S_S_Pz_C1000001001_b+ABX*I_ERI_Dyz_S_S_Pz_C1000001001_b;
  Double I_ERI_D2z_Px_S_Pz_C1000001001_b = I_ERI_Fx2z_S_S_Pz_C1000001001_b+ABX*I_ERI_D2z_S_S_Pz_C1000001001_b;
  Double I_ERI_Dxy_Py_S_Pz_C1000001001_b = I_ERI_Fx2y_S_S_Pz_C1000001001_b+ABY*I_ERI_Dxy_S_S_Pz_C1000001001_b;
  Double I_ERI_Dxz_Py_S_Pz_C1000001001_b = I_ERI_Fxyz_S_S_Pz_C1000001001_b+ABY*I_ERI_Dxz_S_S_Pz_C1000001001_b;
  Double I_ERI_D2y_Py_S_Pz_C1000001001_b = I_ERI_F3y_S_S_Pz_C1000001001_b+ABY*I_ERI_D2y_S_S_Pz_C1000001001_b;
  Double I_ERI_Dyz_Py_S_Pz_C1000001001_b = I_ERI_F2yz_S_S_Pz_C1000001001_b+ABY*I_ERI_Dyz_S_S_Pz_C1000001001_b;
  Double I_ERI_D2z_Py_S_Pz_C1000001001_b = I_ERI_Fy2z_S_S_Pz_C1000001001_b+ABY*I_ERI_D2z_S_S_Pz_C1000001001_b;
  Double I_ERI_Dxz_Pz_S_Pz_C1000001001_b = I_ERI_Fx2z_S_S_Pz_C1000001001_b+ABZ*I_ERI_Dxz_S_S_Pz_C1000001001_b;
  Double I_ERI_Dyz_Pz_S_Pz_C1000001001_b = I_ERI_Fy2z_S_S_Pz_C1000001001_b+ABZ*I_ERI_Dyz_S_S_Pz_C1000001001_b;
  Double I_ERI_D2z_Pz_S_Pz_C1000001001_b = I_ERI_F3z_S_S_Pz_C1000001001_b+ABZ*I_ERI_D2z_S_S_Pz_C1000001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_P_C1000001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000001001_b
   ************************************************************/
  Double I_ERI_Px_D2x_S_Px_C1000001001_b = I_ERI_D2x_Px_S_Px_C1000001001_b+ABX*I_ERI_Px_Px_S_Px_C1000001001_b;
  Double I_ERI_Py_D2x_S_Px_C1000001001_b = I_ERI_Dxy_Px_S_Px_C1000001001_b+ABX*I_ERI_Py_Px_S_Px_C1000001001_b;
  Double I_ERI_Pz_D2x_S_Px_C1000001001_b = I_ERI_Dxz_Px_S_Px_C1000001001_b+ABX*I_ERI_Pz_Px_S_Px_C1000001001_b;
  Double I_ERI_Px_Dxy_S_Px_C1000001001_b = I_ERI_Dxy_Px_S_Px_C1000001001_b+ABY*I_ERI_Px_Px_S_Px_C1000001001_b;
  Double I_ERI_Py_Dxy_S_Px_C1000001001_b = I_ERI_D2y_Px_S_Px_C1000001001_b+ABY*I_ERI_Py_Px_S_Px_C1000001001_b;
  Double I_ERI_Pz_Dxy_S_Px_C1000001001_b = I_ERI_Dyz_Px_S_Px_C1000001001_b+ABY*I_ERI_Pz_Px_S_Px_C1000001001_b;
  Double I_ERI_Px_Dxz_S_Px_C1000001001_b = I_ERI_Dxz_Px_S_Px_C1000001001_b+ABZ*I_ERI_Px_Px_S_Px_C1000001001_b;
  Double I_ERI_Py_Dxz_S_Px_C1000001001_b = I_ERI_Dyz_Px_S_Px_C1000001001_b+ABZ*I_ERI_Py_Px_S_Px_C1000001001_b;
  Double I_ERI_Pz_Dxz_S_Px_C1000001001_b = I_ERI_D2z_Px_S_Px_C1000001001_b+ABZ*I_ERI_Pz_Px_S_Px_C1000001001_b;
  Double I_ERI_Px_D2y_S_Px_C1000001001_b = I_ERI_Dxy_Py_S_Px_C1000001001_b+ABY*I_ERI_Px_Py_S_Px_C1000001001_b;
  Double I_ERI_Py_D2y_S_Px_C1000001001_b = I_ERI_D2y_Py_S_Px_C1000001001_b+ABY*I_ERI_Py_Py_S_Px_C1000001001_b;
  Double I_ERI_Pz_D2y_S_Px_C1000001001_b = I_ERI_Dyz_Py_S_Px_C1000001001_b+ABY*I_ERI_Pz_Py_S_Px_C1000001001_b;
  Double I_ERI_Px_Dyz_S_Px_C1000001001_b = I_ERI_Dxz_Py_S_Px_C1000001001_b+ABZ*I_ERI_Px_Py_S_Px_C1000001001_b;
  Double I_ERI_Py_Dyz_S_Px_C1000001001_b = I_ERI_Dyz_Py_S_Px_C1000001001_b+ABZ*I_ERI_Py_Py_S_Px_C1000001001_b;
  Double I_ERI_Pz_Dyz_S_Px_C1000001001_b = I_ERI_D2z_Py_S_Px_C1000001001_b+ABZ*I_ERI_Pz_Py_S_Px_C1000001001_b;
  Double I_ERI_Px_D2z_S_Px_C1000001001_b = I_ERI_Dxz_Pz_S_Px_C1000001001_b+ABZ*I_ERI_Px_Pz_S_Px_C1000001001_b;
  Double I_ERI_Py_D2z_S_Px_C1000001001_b = I_ERI_Dyz_Pz_S_Px_C1000001001_b+ABZ*I_ERI_Py_Pz_S_Px_C1000001001_b;
  Double I_ERI_Pz_D2z_S_Px_C1000001001_b = I_ERI_D2z_Pz_S_Px_C1000001001_b+ABZ*I_ERI_Pz_Pz_S_Px_C1000001001_b;
  Double I_ERI_Px_D2x_S_Py_C1000001001_b = I_ERI_D2x_Px_S_Py_C1000001001_b+ABX*I_ERI_Px_Px_S_Py_C1000001001_b;
  Double I_ERI_Py_D2x_S_Py_C1000001001_b = I_ERI_Dxy_Px_S_Py_C1000001001_b+ABX*I_ERI_Py_Px_S_Py_C1000001001_b;
  Double I_ERI_Pz_D2x_S_Py_C1000001001_b = I_ERI_Dxz_Px_S_Py_C1000001001_b+ABX*I_ERI_Pz_Px_S_Py_C1000001001_b;
  Double I_ERI_Px_Dxy_S_Py_C1000001001_b = I_ERI_Dxy_Px_S_Py_C1000001001_b+ABY*I_ERI_Px_Px_S_Py_C1000001001_b;
  Double I_ERI_Py_Dxy_S_Py_C1000001001_b = I_ERI_D2y_Px_S_Py_C1000001001_b+ABY*I_ERI_Py_Px_S_Py_C1000001001_b;
  Double I_ERI_Pz_Dxy_S_Py_C1000001001_b = I_ERI_Dyz_Px_S_Py_C1000001001_b+ABY*I_ERI_Pz_Px_S_Py_C1000001001_b;
  Double I_ERI_Px_Dxz_S_Py_C1000001001_b = I_ERI_Dxz_Px_S_Py_C1000001001_b+ABZ*I_ERI_Px_Px_S_Py_C1000001001_b;
  Double I_ERI_Py_Dxz_S_Py_C1000001001_b = I_ERI_Dyz_Px_S_Py_C1000001001_b+ABZ*I_ERI_Py_Px_S_Py_C1000001001_b;
  Double I_ERI_Pz_Dxz_S_Py_C1000001001_b = I_ERI_D2z_Px_S_Py_C1000001001_b+ABZ*I_ERI_Pz_Px_S_Py_C1000001001_b;
  Double I_ERI_Px_D2y_S_Py_C1000001001_b = I_ERI_Dxy_Py_S_Py_C1000001001_b+ABY*I_ERI_Px_Py_S_Py_C1000001001_b;
  Double I_ERI_Py_D2y_S_Py_C1000001001_b = I_ERI_D2y_Py_S_Py_C1000001001_b+ABY*I_ERI_Py_Py_S_Py_C1000001001_b;
  Double I_ERI_Pz_D2y_S_Py_C1000001001_b = I_ERI_Dyz_Py_S_Py_C1000001001_b+ABY*I_ERI_Pz_Py_S_Py_C1000001001_b;
  Double I_ERI_Px_Dyz_S_Py_C1000001001_b = I_ERI_Dxz_Py_S_Py_C1000001001_b+ABZ*I_ERI_Px_Py_S_Py_C1000001001_b;
  Double I_ERI_Py_Dyz_S_Py_C1000001001_b = I_ERI_Dyz_Py_S_Py_C1000001001_b+ABZ*I_ERI_Py_Py_S_Py_C1000001001_b;
  Double I_ERI_Pz_Dyz_S_Py_C1000001001_b = I_ERI_D2z_Py_S_Py_C1000001001_b+ABZ*I_ERI_Pz_Py_S_Py_C1000001001_b;
  Double I_ERI_Px_D2z_S_Py_C1000001001_b = I_ERI_Dxz_Pz_S_Py_C1000001001_b+ABZ*I_ERI_Px_Pz_S_Py_C1000001001_b;
  Double I_ERI_Py_D2z_S_Py_C1000001001_b = I_ERI_Dyz_Pz_S_Py_C1000001001_b+ABZ*I_ERI_Py_Pz_S_Py_C1000001001_b;
  Double I_ERI_Pz_D2z_S_Py_C1000001001_b = I_ERI_D2z_Pz_S_Py_C1000001001_b+ABZ*I_ERI_Pz_Pz_S_Py_C1000001001_b;
  Double I_ERI_Px_D2x_S_Pz_C1000001001_b = I_ERI_D2x_Px_S_Pz_C1000001001_b+ABX*I_ERI_Px_Px_S_Pz_C1000001001_b;
  Double I_ERI_Py_D2x_S_Pz_C1000001001_b = I_ERI_Dxy_Px_S_Pz_C1000001001_b+ABX*I_ERI_Py_Px_S_Pz_C1000001001_b;
  Double I_ERI_Pz_D2x_S_Pz_C1000001001_b = I_ERI_Dxz_Px_S_Pz_C1000001001_b+ABX*I_ERI_Pz_Px_S_Pz_C1000001001_b;
  Double I_ERI_Px_Dxy_S_Pz_C1000001001_b = I_ERI_Dxy_Px_S_Pz_C1000001001_b+ABY*I_ERI_Px_Px_S_Pz_C1000001001_b;
  Double I_ERI_Py_Dxy_S_Pz_C1000001001_b = I_ERI_D2y_Px_S_Pz_C1000001001_b+ABY*I_ERI_Py_Px_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Dxy_S_Pz_C1000001001_b = I_ERI_Dyz_Px_S_Pz_C1000001001_b+ABY*I_ERI_Pz_Px_S_Pz_C1000001001_b;
  Double I_ERI_Px_Dxz_S_Pz_C1000001001_b = I_ERI_Dxz_Px_S_Pz_C1000001001_b+ABZ*I_ERI_Px_Px_S_Pz_C1000001001_b;
  Double I_ERI_Py_Dxz_S_Pz_C1000001001_b = I_ERI_Dyz_Px_S_Pz_C1000001001_b+ABZ*I_ERI_Py_Px_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Dxz_S_Pz_C1000001001_b = I_ERI_D2z_Px_S_Pz_C1000001001_b+ABZ*I_ERI_Pz_Px_S_Pz_C1000001001_b;
  Double I_ERI_Px_D2y_S_Pz_C1000001001_b = I_ERI_Dxy_Py_S_Pz_C1000001001_b+ABY*I_ERI_Px_Py_S_Pz_C1000001001_b;
  Double I_ERI_Py_D2y_S_Pz_C1000001001_b = I_ERI_D2y_Py_S_Pz_C1000001001_b+ABY*I_ERI_Py_Py_S_Pz_C1000001001_b;
  Double I_ERI_Pz_D2y_S_Pz_C1000001001_b = I_ERI_Dyz_Py_S_Pz_C1000001001_b+ABY*I_ERI_Pz_Py_S_Pz_C1000001001_b;
  Double I_ERI_Px_Dyz_S_Pz_C1000001001_b = I_ERI_Dxz_Py_S_Pz_C1000001001_b+ABZ*I_ERI_Px_Py_S_Pz_C1000001001_b;
  Double I_ERI_Py_Dyz_S_Pz_C1000001001_b = I_ERI_Dyz_Py_S_Pz_C1000001001_b+ABZ*I_ERI_Py_Py_S_Pz_C1000001001_b;
  Double I_ERI_Pz_Dyz_S_Pz_C1000001001_b = I_ERI_D2z_Py_S_Pz_C1000001001_b+ABZ*I_ERI_Pz_Py_S_Pz_C1000001001_b;
  Double I_ERI_Px_D2z_S_Pz_C1000001001_b = I_ERI_Dxz_Pz_S_Pz_C1000001001_b+ABZ*I_ERI_Px_Pz_S_Pz_C1000001001_b;
  Double I_ERI_Py_D2z_S_Pz_C1000001001_b = I_ERI_Dyz_Pz_S_Pz_C1000001001_b+ABZ*I_ERI_Py_Pz_S_Pz_C1000001001_b;
  Double I_ERI_Pz_D2z_S_Pz_C1000001001_b = I_ERI_D2z_Pz_S_Pz_C1000001001_b+ABZ*I_ERI_Pz_Pz_S_Pz_C1000001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001000001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001000001_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_C1001000001_b = I_ERI_D2x_S_Px_Px_C1001000001_b+ABX*I_ERI_Px_S_Px_Px_C1001000001_b;
  Double I_ERI_Py_Px_Px_Px_C1001000001_b = I_ERI_Dxy_S_Px_Px_C1001000001_b+ABX*I_ERI_Py_S_Px_Px_C1001000001_b;
  Double I_ERI_Pz_Px_Px_Px_C1001000001_b = I_ERI_Dxz_S_Px_Px_C1001000001_b+ABX*I_ERI_Pz_S_Px_Px_C1001000001_b;
  Double I_ERI_Px_Py_Px_Px_C1001000001_b = I_ERI_Dxy_S_Px_Px_C1001000001_b+ABY*I_ERI_Px_S_Px_Px_C1001000001_b;
  Double I_ERI_Py_Py_Px_Px_C1001000001_b = I_ERI_D2y_S_Px_Px_C1001000001_b+ABY*I_ERI_Py_S_Px_Px_C1001000001_b;
  Double I_ERI_Pz_Py_Px_Px_C1001000001_b = I_ERI_Dyz_S_Px_Px_C1001000001_b+ABY*I_ERI_Pz_S_Px_Px_C1001000001_b;
  Double I_ERI_Px_Pz_Px_Px_C1001000001_b = I_ERI_Dxz_S_Px_Px_C1001000001_b+ABZ*I_ERI_Px_S_Px_Px_C1001000001_b;
  Double I_ERI_Py_Pz_Px_Px_C1001000001_b = I_ERI_Dyz_S_Px_Px_C1001000001_b+ABZ*I_ERI_Py_S_Px_Px_C1001000001_b;
  Double I_ERI_Pz_Pz_Px_Px_C1001000001_b = I_ERI_D2z_S_Px_Px_C1001000001_b+ABZ*I_ERI_Pz_S_Px_Px_C1001000001_b;
  Double I_ERI_Px_Px_Py_Px_C1001000001_b = I_ERI_D2x_S_Py_Px_C1001000001_b+ABX*I_ERI_Px_S_Py_Px_C1001000001_b;
  Double I_ERI_Py_Px_Py_Px_C1001000001_b = I_ERI_Dxy_S_Py_Px_C1001000001_b+ABX*I_ERI_Py_S_Py_Px_C1001000001_b;
  Double I_ERI_Pz_Px_Py_Px_C1001000001_b = I_ERI_Dxz_S_Py_Px_C1001000001_b+ABX*I_ERI_Pz_S_Py_Px_C1001000001_b;
  Double I_ERI_Px_Py_Py_Px_C1001000001_b = I_ERI_Dxy_S_Py_Px_C1001000001_b+ABY*I_ERI_Px_S_Py_Px_C1001000001_b;
  Double I_ERI_Py_Py_Py_Px_C1001000001_b = I_ERI_D2y_S_Py_Px_C1001000001_b+ABY*I_ERI_Py_S_Py_Px_C1001000001_b;
  Double I_ERI_Pz_Py_Py_Px_C1001000001_b = I_ERI_Dyz_S_Py_Px_C1001000001_b+ABY*I_ERI_Pz_S_Py_Px_C1001000001_b;
  Double I_ERI_Px_Pz_Py_Px_C1001000001_b = I_ERI_Dxz_S_Py_Px_C1001000001_b+ABZ*I_ERI_Px_S_Py_Px_C1001000001_b;
  Double I_ERI_Py_Pz_Py_Px_C1001000001_b = I_ERI_Dyz_S_Py_Px_C1001000001_b+ABZ*I_ERI_Py_S_Py_Px_C1001000001_b;
  Double I_ERI_Pz_Pz_Py_Px_C1001000001_b = I_ERI_D2z_S_Py_Px_C1001000001_b+ABZ*I_ERI_Pz_S_Py_Px_C1001000001_b;
  Double I_ERI_Px_Px_Pz_Px_C1001000001_b = I_ERI_D2x_S_Pz_Px_C1001000001_b+ABX*I_ERI_Px_S_Pz_Px_C1001000001_b;
  Double I_ERI_Py_Px_Pz_Px_C1001000001_b = I_ERI_Dxy_S_Pz_Px_C1001000001_b+ABX*I_ERI_Py_S_Pz_Px_C1001000001_b;
  Double I_ERI_Pz_Px_Pz_Px_C1001000001_b = I_ERI_Dxz_S_Pz_Px_C1001000001_b+ABX*I_ERI_Pz_S_Pz_Px_C1001000001_b;
  Double I_ERI_Px_Py_Pz_Px_C1001000001_b = I_ERI_Dxy_S_Pz_Px_C1001000001_b+ABY*I_ERI_Px_S_Pz_Px_C1001000001_b;
  Double I_ERI_Py_Py_Pz_Px_C1001000001_b = I_ERI_D2y_S_Pz_Px_C1001000001_b+ABY*I_ERI_Py_S_Pz_Px_C1001000001_b;
  Double I_ERI_Pz_Py_Pz_Px_C1001000001_b = I_ERI_Dyz_S_Pz_Px_C1001000001_b+ABY*I_ERI_Pz_S_Pz_Px_C1001000001_b;
  Double I_ERI_Px_Pz_Pz_Px_C1001000001_b = I_ERI_Dxz_S_Pz_Px_C1001000001_b+ABZ*I_ERI_Px_S_Pz_Px_C1001000001_b;
  Double I_ERI_Py_Pz_Pz_Px_C1001000001_b = I_ERI_Dyz_S_Pz_Px_C1001000001_b+ABZ*I_ERI_Py_S_Pz_Px_C1001000001_b;
  Double I_ERI_Pz_Pz_Pz_Px_C1001000001_b = I_ERI_D2z_S_Pz_Px_C1001000001_b+ABZ*I_ERI_Pz_S_Pz_Px_C1001000001_b;
  Double I_ERI_Px_Px_Px_Py_C1001000001_b = I_ERI_D2x_S_Px_Py_C1001000001_b+ABX*I_ERI_Px_S_Px_Py_C1001000001_b;
  Double I_ERI_Py_Px_Px_Py_C1001000001_b = I_ERI_Dxy_S_Px_Py_C1001000001_b+ABX*I_ERI_Py_S_Px_Py_C1001000001_b;
  Double I_ERI_Pz_Px_Px_Py_C1001000001_b = I_ERI_Dxz_S_Px_Py_C1001000001_b+ABX*I_ERI_Pz_S_Px_Py_C1001000001_b;
  Double I_ERI_Px_Py_Px_Py_C1001000001_b = I_ERI_Dxy_S_Px_Py_C1001000001_b+ABY*I_ERI_Px_S_Px_Py_C1001000001_b;
  Double I_ERI_Py_Py_Px_Py_C1001000001_b = I_ERI_D2y_S_Px_Py_C1001000001_b+ABY*I_ERI_Py_S_Px_Py_C1001000001_b;
  Double I_ERI_Pz_Py_Px_Py_C1001000001_b = I_ERI_Dyz_S_Px_Py_C1001000001_b+ABY*I_ERI_Pz_S_Px_Py_C1001000001_b;
  Double I_ERI_Px_Pz_Px_Py_C1001000001_b = I_ERI_Dxz_S_Px_Py_C1001000001_b+ABZ*I_ERI_Px_S_Px_Py_C1001000001_b;
  Double I_ERI_Py_Pz_Px_Py_C1001000001_b = I_ERI_Dyz_S_Px_Py_C1001000001_b+ABZ*I_ERI_Py_S_Px_Py_C1001000001_b;
  Double I_ERI_Pz_Pz_Px_Py_C1001000001_b = I_ERI_D2z_S_Px_Py_C1001000001_b+ABZ*I_ERI_Pz_S_Px_Py_C1001000001_b;
  Double I_ERI_Px_Px_Py_Py_C1001000001_b = I_ERI_D2x_S_Py_Py_C1001000001_b+ABX*I_ERI_Px_S_Py_Py_C1001000001_b;
  Double I_ERI_Py_Px_Py_Py_C1001000001_b = I_ERI_Dxy_S_Py_Py_C1001000001_b+ABX*I_ERI_Py_S_Py_Py_C1001000001_b;
  Double I_ERI_Pz_Px_Py_Py_C1001000001_b = I_ERI_Dxz_S_Py_Py_C1001000001_b+ABX*I_ERI_Pz_S_Py_Py_C1001000001_b;
  Double I_ERI_Px_Py_Py_Py_C1001000001_b = I_ERI_Dxy_S_Py_Py_C1001000001_b+ABY*I_ERI_Px_S_Py_Py_C1001000001_b;
  Double I_ERI_Py_Py_Py_Py_C1001000001_b = I_ERI_D2y_S_Py_Py_C1001000001_b+ABY*I_ERI_Py_S_Py_Py_C1001000001_b;
  Double I_ERI_Pz_Py_Py_Py_C1001000001_b = I_ERI_Dyz_S_Py_Py_C1001000001_b+ABY*I_ERI_Pz_S_Py_Py_C1001000001_b;
  Double I_ERI_Px_Pz_Py_Py_C1001000001_b = I_ERI_Dxz_S_Py_Py_C1001000001_b+ABZ*I_ERI_Px_S_Py_Py_C1001000001_b;
  Double I_ERI_Py_Pz_Py_Py_C1001000001_b = I_ERI_Dyz_S_Py_Py_C1001000001_b+ABZ*I_ERI_Py_S_Py_Py_C1001000001_b;
  Double I_ERI_Pz_Pz_Py_Py_C1001000001_b = I_ERI_D2z_S_Py_Py_C1001000001_b+ABZ*I_ERI_Pz_S_Py_Py_C1001000001_b;
  Double I_ERI_Px_Px_Pz_Py_C1001000001_b = I_ERI_D2x_S_Pz_Py_C1001000001_b+ABX*I_ERI_Px_S_Pz_Py_C1001000001_b;
  Double I_ERI_Py_Px_Pz_Py_C1001000001_b = I_ERI_Dxy_S_Pz_Py_C1001000001_b+ABX*I_ERI_Py_S_Pz_Py_C1001000001_b;
  Double I_ERI_Pz_Px_Pz_Py_C1001000001_b = I_ERI_Dxz_S_Pz_Py_C1001000001_b+ABX*I_ERI_Pz_S_Pz_Py_C1001000001_b;
  Double I_ERI_Px_Py_Pz_Py_C1001000001_b = I_ERI_Dxy_S_Pz_Py_C1001000001_b+ABY*I_ERI_Px_S_Pz_Py_C1001000001_b;
  Double I_ERI_Py_Py_Pz_Py_C1001000001_b = I_ERI_D2y_S_Pz_Py_C1001000001_b+ABY*I_ERI_Py_S_Pz_Py_C1001000001_b;
  Double I_ERI_Pz_Py_Pz_Py_C1001000001_b = I_ERI_Dyz_S_Pz_Py_C1001000001_b+ABY*I_ERI_Pz_S_Pz_Py_C1001000001_b;
  Double I_ERI_Px_Pz_Pz_Py_C1001000001_b = I_ERI_Dxz_S_Pz_Py_C1001000001_b+ABZ*I_ERI_Px_S_Pz_Py_C1001000001_b;
  Double I_ERI_Py_Pz_Pz_Py_C1001000001_b = I_ERI_Dyz_S_Pz_Py_C1001000001_b+ABZ*I_ERI_Py_S_Pz_Py_C1001000001_b;
  Double I_ERI_Pz_Pz_Pz_Py_C1001000001_b = I_ERI_D2z_S_Pz_Py_C1001000001_b+ABZ*I_ERI_Pz_S_Pz_Py_C1001000001_b;
  Double I_ERI_Px_Px_Px_Pz_C1001000001_b = I_ERI_D2x_S_Px_Pz_C1001000001_b+ABX*I_ERI_Px_S_Px_Pz_C1001000001_b;
  Double I_ERI_Py_Px_Px_Pz_C1001000001_b = I_ERI_Dxy_S_Px_Pz_C1001000001_b+ABX*I_ERI_Py_S_Px_Pz_C1001000001_b;
  Double I_ERI_Pz_Px_Px_Pz_C1001000001_b = I_ERI_Dxz_S_Px_Pz_C1001000001_b+ABX*I_ERI_Pz_S_Px_Pz_C1001000001_b;
  Double I_ERI_Px_Py_Px_Pz_C1001000001_b = I_ERI_Dxy_S_Px_Pz_C1001000001_b+ABY*I_ERI_Px_S_Px_Pz_C1001000001_b;
  Double I_ERI_Py_Py_Px_Pz_C1001000001_b = I_ERI_D2y_S_Px_Pz_C1001000001_b+ABY*I_ERI_Py_S_Px_Pz_C1001000001_b;
  Double I_ERI_Pz_Py_Px_Pz_C1001000001_b = I_ERI_Dyz_S_Px_Pz_C1001000001_b+ABY*I_ERI_Pz_S_Px_Pz_C1001000001_b;
  Double I_ERI_Px_Pz_Px_Pz_C1001000001_b = I_ERI_Dxz_S_Px_Pz_C1001000001_b+ABZ*I_ERI_Px_S_Px_Pz_C1001000001_b;
  Double I_ERI_Py_Pz_Px_Pz_C1001000001_b = I_ERI_Dyz_S_Px_Pz_C1001000001_b+ABZ*I_ERI_Py_S_Px_Pz_C1001000001_b;
  Double I_ERI_Pz_Pz_Px_Pz_C1001000001_b = I_ERI_D2z_S_Px_Pz_C1001000001_b+ABZ*I_ERI_Pz_S_Px_Pz_C1001000001_b;
  Double I_ERI_Px_Px_Py_Pz_C1001000001_b = I_ERI_D2x_S_Py_Pz_C1001000001_b+ABX*I_ERI_Px_S_Py_Pz_C1001000001_b;
  Double I_ERI_Py_Px_Py_Pz_C1001000001_b = I_ERI_Dxy_S_Py_Pz_C1001000001_b+ABX*I_ERI_Py_S_Py_Pz_C1001000001_b;
  Double I_ERI_Pz_Px_Py_Pz_C1001000001_b = I_ERI_Dxz_S_Py_Pz_C1001000001_b+ABX*I_ERI_Pz_S_Py_Pz_C1001000001_b;
  Double I_ERI_Px_Py_Py_Pz_C1001000001_b = I_ERI_Dxy_S_Py_Pz_C1001000001_b+ABY*I_ERI_Px_S_Py_Pz_C1001000001_b;
  Double I_ERI_Py_Py_Py_Pz_C1001000001_b = I_ERI_D2y_S_Py_Pz_C1001000001_b+ABY*I_ERI_Py_S_Py_Pz_C1001000001_b;
  Double I_ERI_Pz_Py_Py_Pz_C1001000001_b = I_ERI_Dyz_S_Py_Pz_C1001000001_b+ABY*I_ERI_Pz_S_Py_Pz_C1001000001_b;
  Double I_ERI_Px_Pz_Py_Pz_C1001000001_b = I_ERI_Dxz_S_Py_Pz_C1001000001_b+ABZ*I_ERI_Px_S_Py_Pz_C1001000001_b;
  Double I_ERI_Py_Pz_Py_Pz_C1001000001_b = I_ERI_Dyz_S_Py_Pz_C1001000001_b+ABZ*I_ERI_Py_S_Py_Pz_C1001000001_b;
  Double I_ERI_Pz_Pz_Py_Pz_C1001000001_b = I_ERI_D2z_S_Py_Pz_C1001000001_b+ABZ*I_ERI_Pz_S_Py_Pz_C1001000001_b;
  Double I_ERI_Px_Px_Pz_Pz_C1001000001_b = I_ERI_D2x_S_Pz_Pz_C1001000001_b+ABX*I_ERI_Px_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Py_Px_Pz_Pz_C1001000001_b = I_ERI_Dxy_S_Pz_Pz_C1001000001_b+ABX*I_ERI_Py_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Pz_Px_Pz_Pz_C1001000001_b = I_ERI_Dxz_S_Pz_Pz_C1001000001_b+ABX*I_ERI_Pz_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Px_Py_Pz_Pz_C1001000001_b = I_ERI_Dxy_S_Pz_Pz_C1001000001_b+ABY*I_ERI_Px_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Py_Py_Pz_Pz_C1001000001_b = I_ERI_D2y_S_Pz_Pz_C1001000001_b+ABY*I_ERI_Py_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Pz_Py_Pz_Pz_C1001000001_b = I_ERI_Dyz_S_Pz_Pz_C1001000001_b+ABY*I_ERI_Pz_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Px_Pz_Pz_Pz_C1001000001_b = I_ERI_Dxz_S_Pz_Pz_C1001000001_b+ABZ*I_ERI_Px_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Py_Pz_Pz_Pz_C1001000001_b = I_ERI_Dyz_S_Pz_Pz_C1001000001_b+ABZ*I_ERI_Py_S_Pz_Pz_C1001000001_b;
  Double I_ERI_Pz_Pz_Pz_Pz_C1001000001_b = I_ERI_D2z_S_Pz_Pz_C1001000001_b+ABZ*I_ERI_Pz_S_Pz_Pz_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001001_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_C1001001001_b = I_ERI_D2x_S_Px_Px_C1001001001_b+ABX*I_ERI_Px_S_Px_Px_C1001001001_b;
  Double I_ERI_Py_Px_Px_Px_C1001001001_b = I_ERI_Dxy_S_Px_Px_C1001001001_b+ABX*I_ERI_Py_S_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Px_Px_Px_C1001001001_b = I_ERI_Dxz_S_Px_Px_C1001001001_b+ABX*I_ERI_Pz_S_Px_Px_C1001001001_b;
  Double I_ERI_Px_Py_Px_Px_C1001001001_b = I_ERI_Dxy_S_Px_Px_C1001001001_b+ABY*I_ERI_Px_S_Px_Px_C1001001001_b;
  Double I_ERI_Py_Py_Px_Px_C1001001001_b = I_ERI_D2y_S_Px_Px_C1001001001_b+ABY*I_ERI_Py_S_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Py_Px_Px_C1001001001_b = I_ERI_Dyz_S_Px_Px_C1001001001_b+ABY*I_ERI_Pz_S_Px_Px_C1001001001_b;
  Double I_ERI_Px_Pz_Px_Px_C1001001001_b = I_ERI_Dxz_S_Px_Px_C1001001001_b+ABZ*I_ERI_Px_S_Px_Px_C1001001001_b;
  Double I_ERI_Py_Pz_Px_Px_C1001001001_b = I_ERI_Dyz_S_Px_Px_C1001001001_b+ABZ*I_ERI_Py_S_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Pz_Px_Px_C1001001001_b = I_ERI_D2z_S_Px_Px_C1001001001_b+ABZ*I_ERI_Pz_S_Px_Px_C1001001001_b;
  Double I_ERI_Px_Px_Py_Px_C1001001001_b = I_ERI_D2x_S_Py_Px_C1001001001_b+ABX*I_ERI_Px_S_Py_Px_C1001001001_b;
  Double I_ERI_Py_Px_Py_Px_C1001001001_b = I_ERI_Dxy_S_Py_Px_C1001001001_b+ABX*I_ERI_Py_S_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Px_Py_Px_C1001001001_b = I_ERI_Dxz_S_Py_Px_C1001001001_b+ABX*I_ERI_Pz_S_Py_Px_C1001001001_b;
  Double I_ERI_Px_Py_Py_Px_C1001001001_b = I_ERI_Dxy_S_Py_Px_C1001001001_b+ABY*I_ERI_Px_S_Py_Px_C1001001001_b;
  Double I_ERI_Py_Py_Py_Px_C1001001001_b = I_ERI_D2y_S_Py_Px_C1001001001_b+ABY*I_ERI_Py_S_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Py_Py_Px_C1001001001_b = I_ERI_Dyz_S_Py_Px_C1001001001_b+ABY*I_ERI_Pz_S_Py_Px_C1001001001_b;
  Double I_ERI_Px_Pz_Py_Px_C1001001001_b = I_ERI_Dxz_S_Py_Px_C1001001001_b+ABZ*I_ERI_Px_S_Py_Px_C1001001001_b;
  Double I_ERI_Py_Pz_Py_Px_C1001001001_b = I_ERI_Dyz_S_Py_Px_C1001001001_b+ABZ*I_ERI_Py_S_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Pz_Py_Px_C1001001001_b = I_ERI_D2z_S_Py_Px_C1001001001_b+ABZ*I_ERI_Pz_S_Py_Px_C1001001001_b;
  Double I_ERI_Px_Px_Pz_Px_C1001001001_b = I_ERI_D2x_S_Pz_Px_C1001001001_b+ABX*I_ERI_Px_S_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Px_Pz_Px_C1001001001_b = I_ERI_Dxy_S_Pz_Px_C1001001001_b+ABX*I_ERI_Py_S_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Px_Pz_Px_C1001001001_b = I_ERI_Dxz_S_Pz_Px_C1001001001_b+ABX*I_ERI_Pz_S_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Py_Pz_Px_C1001001001_b = I_ERI_Dxy_S_Pz_Px_C1001001001_b+ABY*I_ERI_Px_S_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Py_Pz_Px_C1001001001_b = I_ERI_D2y_S_Pz_Px_C1001001001_b+ABY*I_ERI_Py_S_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Py_Pz_Px_C1001001001_b = I_ERI_Dyz_S_Pz_Px_C1001001001_b+ABY*I_ERI_Pz_S_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Pz_Pz_Px_C1001001001_b = I_ERI_Dxz_S_Pz_Px_C1001001001_b+ABZ*I_ERI_Px_S_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Pz_Pz_Px_C1001001001_b = I_ERI_Dyz_S_Pz_Px_C1001001001_b+ABZ*I_ERI_Py_S_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Pz_Pz_Px_C1001001001_b = I_ERI_D2z_S_Pz_Px_C1001001001_b+ABZ*I_ERI_Pz_S_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Px_Px_Py_C1001001001_b = I_ERI_D2x_S_Px_Py_C1001001001_b+ABX*I_ERI_Px_S_Px_Py_C1001001001_b;
  Double I_ERI_Py_Px_Px_Py_C1001001001_b = I_ERI_Dxy_S_Px_Py_C1001001001_b+ABX*I_ERI_Py_S_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Px_Px_Py_C1001001001_b = I_ERI_Dxz_S_Px_Py_C1001001001_b+ABX*I_ERI_Pz_S_Px_Py_C1001001001_b;
  Double I_ERI_Px_Py_Px_Py_C1001001001_b = I_ERI_Dxy_S_Px_Py_C1001001001_b+ABY*I_ERI_Px_S_Px_Py_C1001001001_b;
  Double I_ERI_Py_Py_Px_Py_C1001001001_b = I_ERI_D2y_S_Px_Py_C1001001001_b+ABY*I_ERI_Py_S_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Py_Px_Py_C1001001001_b = I_ERI_Dyz_S_Px_Py_C1001001001_b+ABY*I_ERI_Pz_S_Px_Py_C1001001001_b;
  Double I_ERI_Px_Pz_Px_Py_C1001001001_b = I_ERI_Dxz_S_Px_Py_C1001001001_b+ABZ*I_ERI_Px_S_Px_Py_C1001001001_b;
  Double I_ERI_Py_Pz_Px_Py_C1001001001_b = I_ERI_Dyz_S_Px_Py_C1001001001_b+ABZ*I_ERI_Py_S_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Pz_Px_Py_C1001001001_b = I_ERI_D2z_S_Px_Py_C1001001001_b+ABZ*I_ERI_Pz_S_Px_Py_C1001001001_b;
  Double I_ERI_Px_Px_Py_Py_C1001001001_b = I_ERI_D2x_S_Py_Py_C1001001001_b+ABX*I_ERI_Px_S_Py_Py_C1001001001_b;
  Double I_ERI_Py_Px_Py_Py_C1001001001_b = I_ERI_Dxy_S_Py_Py_C1001001001_b+ABX*I_ERI_Py_S_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Px_Py_Py_C1001001001_b = I_ERI_Dxz_S_Py_Py_C1001001001_b+ABX*I_ERI_Pz_S_Py_Py_C1001001001_b;
  Double I_ERI_Px_Py_Py_Py_C1001001001_b = I_ERI_Dxy_S_Py_Py_C1001001001_b+ABY*I_ERI_Px_S_Py_Py_C1001001001_b;
  Double I_ERI_Py_Py_Py_Py_C1001001001_b = I_ERI_D2y_S_Py_Py_C1001001001_b+ABY*I_ERI_Py_S_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Py_Py_Py_C1001001001_b = I_ERI_Dyz_S_Py_Py_C1001001001_b+ABY*I_ERI_Pz_S_Py_Py_C1001001001_b;
  Double I_ERI_Px_Pz_Py_Py_C1001001001_b = I_ERI_Dxz_S_Py_Py_C1001001001_b+ABZ*I_ERI_Px_S_Py_Py_C1001001001_b;
  Double I_ERI_Py_Pz_Py_Py_C1001001001_b = I_ERI_Dyz_S_Py_Py_C1001001001_b+ABZ*I_ERI_Py_S_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Pz_Py_Py_C1001001001_b = I_ERI_D2z_S_Py_Py_C1001001001_b+ABZ*I_ERI_Pz_S_Py_Py_C1001001001_b;
  Double I_ERI_Px_Px_Pz_Py_C1001001001_b = I_ERI_D2x_S_Pz_Py_C1001001001_b+ABX*I_ERI_Px_S_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Px_Pz_Py_C1001001001_b = I_ERI_Dxy_S_Pz_Py_C1001001001_b+ABX*I_ERI_Py_S_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Px_Pz_Py_C1001001001_b = I_ERI_Dxz_S_Pz_Py_C1001001001_b+ABX*I_ERI_Pz_S_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Py_Pz_Py_C1001001001_b = I_ERI_Dxy_S_Pz_Py_C1001001001_b+ABY*I_ERI_Px_S_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Py_Pz_Py_C1001001001_b = I_ERI_D2y_S_Pz_Py_C1001001001_b+ABY*I_ERI_Py_S_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Py_Pz_Py_C1001001001_b = I_ERI_Dyz_S_Pz_Py_C1001001001_b+ABY*I_ERI_Pz_S_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Pz_Pz_Py_C1001001001_b = I_ERI_Dxz_S_Pz_Py_C1001001001_b+ABZ*I_ERI_Px_S_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Pz_Pz_Py_C1001001001_b = I_ERI_Dyz_S_Pz_Py_C1001001001_b+ABZ*I_ERI_Py_S_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Pz_Pz_Py_C1001001001_b = I_ERI_D2z_S_Pz_Py_C1001001001_b+ABZ*I_ERI_Pz_S_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Px_Px_Pz_C1001001001_b = I_ERI_D2x_S_Px_Pz_C1001001001_b+ABX*I_ERI_Px_S_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Px_Px_Pz_C1001001001_b = I_ERI_Dxy_S_Px_Pz_C1001001001_b+ABX*I_ERI_Py_S_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Px_Px_Pz_C1001001001_b = I_ERI_Dxz_S_Px_Pz_C1001001001_b+ABX*I_ERI_Pz_S_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Py_Px_Pz_C1001001001_b = I_ERI_Dxy_S_Px_Pz_C1001001001_b+ABY*I_ERI_Px_S_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Py_Px_Pz_C1001001001_b = I_ERI_D2y_S_Px_Pz_C1001001001_b+ABY*I_ERI_Py_S_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Py_Px_Pz_C1001001001_b = I_ERI_Dyz_S_Px_Pz_C1001001001_b+ABY*I_ERI_Pz_S_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Pz_Px_Pz_C1001001001_b = I_ERI_Dxz_S_Px_Pz_C1001001001_b+ABZ*I_ERI_Px_S_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Pz_Px_Pz_C1001001001_b = I_ERI_Dyz_S_Px_Pz_C1001001001_b+ABZ*I_ERI_Py_S_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Pz_Px_Pz_C1001001001_b = I_ERI_D2z_S_Px_Pz_C1001001001_b+ABZ*I_ERI_Pz_S_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Px_Py_Pz_C1001001001_b = I_ERI_D2x_S_Py_Pz_C1001001001_b+ABX*I_ERI_Px_S_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Px_Py_Pz_C1001001001_b = I_ERI_Dxy_S_Py_Pz_C1001001001_b+ABX*I_ERI_Py_S_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Px_Py_Pz_C1001001001_b = I_ERI_Dxz_S_Py_Pz_C1001001001_b+ABX*I_ERI_Pz_S_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Py_Py_Pz_C1001001001_b = I_ERI_Dxy_S_Py_Pz_C1001001001_b+ABY*I_ERI_Px_S_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Py_Py_Pz_C1001001001_b = I_ERI_D2y_S_Py_Pz_C1001001001_b+ABY*I_ERI_Py_S_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Py_Py_Pz_C1001001001_b = I_ERI_Dyz_S_Py_Pz_C1001001001_b+ABY*I_ERI_Pz_S_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Pz_Py_Pz_C1001001001_b = I_ERI_Dxz_S_Py_Pz_C1001001001_b+ABZ*I_ERI_Px_S_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Pz_Py_Pz_C1001001001_b = I_ERI_Dyz_S_Py_Pz_C1001001001_b+ABZ*I_ERI_Py_S_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Pz_Py_Pz_C1001001001_b = I_ERI_D2z_S_Py_Pz_C1001001001_b+ABZ*I_ERI_Pz_S_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Px_Pz_Pz_C1001001001_b = I_ERI_D2x_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Px_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Px_Pz_Pz_C1001001001_b = I_ERI_Dxy_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Py_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Px_Pz_Pz_C1001001001_b = I_ERI_Dxz_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Pz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_Py_Pz_Pz_C1001001001_b = I_ERI_Dxy_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Px_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Py_Pz_Pz_C1001001001_b = I_ERI_D2y_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Py_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Py_Pz_Pz_C1001001001_b = I_ERI_Dyz_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Pz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_Pz_Pz_Pz_C1001001001_b = I_ERI_Dxz_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_Px_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Pz_Pz_Pz_C1001001001_b = I_ERI_Dyz_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_Py_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Pz_Pz_Pz_C1001001001_b = I_ERI_D2z_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_Pz_S_Pz_Pz_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_P_C1001001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001001_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_Px_C1001001001_b = I_ERI_F3x_S_Px_Px_C1001001001_b+ABX*I_ERI_D2x_S_Px_Px_C1001001001_b;
  Double I_ERI_Dxy_Px_Px_Px_C1001001001_b = I_ERI_F2xy_S_Px_Px_C1001001001_b+ABX*I_ERI_Dxy_S_Px_Px_C1001001001_b;
  Double I_ERI_Dxz_Px_Px_Px_C1001001001_b = I_ERI_F2xz_S_Px_Px_C1001001001_b+ABX*I_ERI_Dxz_S_Px_Px_C1001001001_b;
  Double I_ERI_D2y_Px_Px_Px_C1001001001_b = I_ERI_Fx2y_S_Px_Px_C1001001001_b+ABX*I_ERI_D2y_S_Px_Px_C1001001001_b;
  Double I_ERI_Dyz_Px_Px_Px_C1001001001_b = I_ERI_Fxyz_S_Px_Px_C1001001001_b+ABX*I_ERI_Dyz_S_Px_Px_C1001001001_b;
  Double I_ERI_D2z_Px_Px_Px_C1001001001_b = I_ERI_Fx2z_S_Px_Px_C1001001001_b+ABX*I_ERI_D2z_S_Px_Px_C1001001001_b;
  Double I_ERI_Dxy_Py_Px_Px_C1001001001_b = I_ERI_Fx2y_S_Px_Px_C1001001001_b+ABY*I_ERI_Dxy_S_Px_Px_C1001001001_b;
  Double I_ERI_Dxz_Py_Px_Px_C1001001001_b = I_ERI_Fxyz_S_Px_Px_C1001001001_b+ABY*I_ERI_Dxz_S_Px_Px_C1001001001_b;
  Double I_ERI_D2y_Py_Px_Px_C1001001001_b = I_ERI_F3y_S_Px_Px_C1001001001_b+ABY*I_ERI_D2y_S_Px_Px_C1001001001_b;
  Double I_ERI_Dyz_Py_Px_Px_C1001001001_b = I_ERI_F2yz_S_Px_Px_C1001001001_b+ABY*I_ERI_Dyz_S_Px_Px_C1001001001_b;
  Double I_ERI_D2z_Py_Px_Px_C1001001001_b = I_ERI_Fy2z_S_Px_Px_C1001001001_b+ABY*I_ERI_D2z_S_Px_Px_C1001001001_b;
  Double I_ERI_Dxz_Pz_Px_Px_C1001001001_b = I_ERI_Fx2z_S_Px_Px_C1001001001_b+ABZ*I_ERI_Dxz_S_Px_Px_C1001001001_b;
  Double I_ERI_Dyz_Pz_Px_Px_C1001001001_b = I_ERI_Fy2z_S_Px_Px_C1001001001_b+ABZ*I_ERI_Dyz_S_Px_Px_C1001001001_b;
  Double I_ERI_D2z_Pz_Px_Px_C1001001001_b = I_ERI_F3z_S_Px_Px_C1001001001_b+ABZ*I_ERI_D2z_S_Px_Px_C1001001001_b;
  Double I_ERI_D2x_Px_Py_Px_C1001001001_b = I_ERI_F3x_S_Py_Px_C1001001001_b+ABX*I_ERI_D2x_S_Py_Px_C1001001001_b;
  Double I_ERI_Dxy_Px_Py_Px_C1001001001_b = I_ERI_F2xy_S_Py_Px_C1001001001_b+ABX*I_ERI_Dxy_S_Py_Px_C1001001001_b;
  Double I_ERI_Dxz_Px_Py_Px_C1001001001_b = I_ERI_F2xz_S_Py_Px_C1001001001_b+ABX*I_ERI_Dxz_S_Py_Px_C1001001001_b;
  Double I_ERI_D2y_Px_Py_Px_C1001001001_b = I_ERI_Fx2y_S_Py_Px_C1001001001_b+ABX*I_ERI_D2y_S_Py_Px_C1001001001_b;
  Double I_ERI_Dyz_Px_Py_Px_C1001001001_b = I_ERI_Fxyz_S_Py_Px_C1001001001_b+ABX*I_ERI_Dyz_S_Py_Px_C1001001001_b;
  Double I_ERI_D2z_Px_Py_Px_C1001001001_b = I_ERI_Fx2z_S_Py_Px_C1001001001_b+ABX*I_ERI_D2z_S_Py_Px_C1001001001_b;
  Double I_ERI_Dxy_Py_Py_Px_C1001001001_b = I_ERI_Fx2y_S_Py_Px_C1001001001_b+ABY*I_ERI_Dxy_S_Py_Px_C1001001001_b;
  Double I_ERI_Dxz_Py_Py_Px_C1001001001_b = I_ERI_Fxyz_S_Py_Px_C1001001001_b+ABY*I_ERI_Dxz_S_Py_Px_C1001001001_b;
  Double I_ERI_D2y_Py_Py_Px_C1001001001_b = I_ERI_F3y_S_Py_Px_C1001001001_b+ABY*I_ERI_D2y_S_Py_Px_C1001001001_b;
  Double I_ERI_Dyz_Py_Py_Px_C1001001001_b = I_ERI_F2yz_S_Py_Px_C1001001001_b+ABY*I_ERI_Dyz_S_Py_Px_C1001001001_b;
  Double I_ERI_D2z_Py_Py_Px_C1001001001_b = I_ERI_Fy2z_S_Py_Px_C1001001001_b+ABY*I_ERI_D2z_S_Py_Px_C1001001001_b;
  Double I_ERI_Dxz_Pz_Py_Px_C1001001001_b = I_ERI_Fx2z_S_Py_Px_C1001001001_b+ABZ*I_ERI_Dxz_S_Py_Px_C1001001001_b;
  Double I_ERI_Dyz_Pz_Py_Px_C1001001001_b = I_ERI_Fy2z_S_Py_Px_C1001001001_b+ABZ*I_ERI_Dyz_S_Py_Px_C1001001001_b;
  Double I_ERI_D2z_Pz_Py_Px_C1001001001_b = I_ERI_F3z_S_Py_Px_C1001001001_b+ABZ*I_ERI_D2z_S_Py_Px_C1001001001_b;
  Double I_ERI_D2x_Px_Pz_Px_C1001001001_b = I_ERI_F3x_S_Pz_Px_C1001001001_b+ABX*I_ERI_D2x_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dxy_Px_Pz_Px_C1001001001_b = I_ERI_F2xy_S_Pz_Px_C1001001001_b+ABX*I_ERI_Dxy_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dxz_Px_Pz_Px_C1001001001_b = I_ERI_F2xz_S_Pz_Px_C1001001001_b+ABX*I_ERI_Dxz_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2y_Px_Pz_Px_C1001001001_b = I_ERI_Fx2y_S_Pz_Px_C1001001001_b+ABX*I_ERI_D2y_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dyz_Px_Pz_Px_C1001001001_b = I_ERI_Fxyz_S_Pz_Px_C1001001001_b+ABX*I_ERI_Dyz_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2z_Px_Pz_Px_C1001001001_b = I_ERI_Fx2z_S_Pz_Px_C1001001001_b+ABX*I_ERI_D2z_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dxy_Py_Pz_Px_C1001001001_b = I_ERI_Fx2y_S_Pz_Px_C1001001001_b+ABY*I_ERI_Dxy_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dxz_Py_Pz_Px_C1001001001_b = I_ERI_Fxyz_S_Pz_Px_C1001001001_b+ABY*I_ERI_Dxz_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2y_Py_Pz_Px_C1001001001_b = I_ERI_F3y_S_Pz_Px_C1001001001_b+ABY*I_ERI_D2y_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dyz_Py_Pz_Px_C1001001001_b = I_ERI_F2yz_S_Pz_Px_C1001001001_b+ABY*I_ERI_Dyz_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2z_Py_Pz_Px_C1001001001_b = I_ERI_Fy2z_S_Pz_Px_C1001001001_b+ABY*I_ERI_D2z_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dxz_Pz_Pz_Px_C1001001001_b = I_ERI_Fx2z_S_Pz_Px_C1001001001_b+ABZ*I_ERI_Dxz_S_Pz_Px_C1001001001_b;
  Double I_ERI_Dyz_Pz_Pz_Px_C1001001001_b = I_ERI_Fy2z_S_Pz_Px_C1001001001_b+ABZ*I_ERI_Dyz_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2z_Pz_Pz_Px_C1001001001_b = I_ERI_F3z_S_Pz_Px_C1001001001_b+ABZ*I_ERI_D2z_S_Pz_Px_C1001001001_b;
  Double I_ERI_D2x_Px_Px_Py_C1001001001_b = I_ERI_F3x_S_Px_Py_C1001001001_b+ABX*I_ERI_D2x_S_Px_Py_C1001001001_b;
  Double I_ERI_Dxy_Px_Px_Py_C1001001001_b = I_ERI_F2xy_S_Px_Py_C1001001001_b+ABX*I_ERI_Dxy_S_Px_Py_C1001001001_b;
  Double I_ERI_Dxz_Px_Px_Py_C1001001001_b = I_ERI_F2xz_S_Px_Py_C1001001001_b+ABX*I_ERI_Dxz_S_Px_Py_C1001001001_b;
  Double I_ERI_D2y_Px_Px_Py_C1001001001_b = I_ERI_Fx2y_S_Px_Py_C1001001001_b+ABX*I_ERI_D2y_S_Px_Py_C1001001001_b;
  Double I_ERI_Dyz_Px_Px_Py_C1001001001_b = I_ERI_Fxyz_S_Px_Py_C1001001001_b+ABX*I_ERI_Dyz_S_Px_Py_C1001001001_b;
  Double I_ERI_D2z_Px_Px_Py_C1001001001_b = I_ERI_Fx2z_S_Px_Py_C1001001001_b+ABX*I_ERI_D2z_S_Px_Py_C1001001001_b;
  Double I_ERI_Dxy_Py_Px_Py_C1001001001_b = I_ERI_Fx2y_S_Px_Py_C1001001001_b+ABY*I_ERI_Dxy_S_Px_Py_C1001001001_b;
  Double I_ERI_Dxz_Py_Px_Py_C1001001001_b = I_ERI_Fxyz_S_Px_Py_C1001001001_b+ABY*I_ERI_Dxz_S_Px_Py_C1001001001_b;
  Double I_ERI_D2y_Py_Px_Py_C1001001001_b = I_ERI_F3y_S_Px_Py_C1001001001_b+ABY*I_ERI_D2y_S_Px_Py_C1001001001_b;
  Double I_ERI_Dyz_Py_Px_Py_C1001001001_b = I_ERI_F2yz_S_Px_Py_C1001001001_b+ABY*I_ERI_Dyz_S_Px_Py_C1001001001_b;
  Double I_ERI_D2z_Py_Px_Py_C1001001001_b = I_ERI_Fy2z_S_Px_Py_C1001001001_b+ABY*I_ERI_D2z_S_Px_Py_C1001001001_b;
  Double I_ERI_Dxz_Pz_Px_Py_C1001001001_b = I_ERI_Fx2z_S_Px_Py_C1001001001_b+ABZ*I_ERI_Dxz_S_Px_Py_C1001001001_b;
  Double I_ERI_Dyz_Pz_Px_Py_C1001001001_b = I_ERI_Fy2z_S_Px_Py_C1001001001_b+ABZ*I_ERI_Dyz_S_Px_Py_C1001001001_b;
  Double I_ERI_D2z_Pz_Px_Py_C1001001001_b = I_ERI_F3z_S_Px_Py_C1001001001_b+ABZ*I_ERI_D2z_S_Px_Py_C1001001001_b;
  Double I_ERI_D2x_Px_Py_Py_C1001001001_b = I_ERI_F3x_S_Py_Py_C1001001001_b+ABX*I_ERI_D2x_S_Py_Py_C1001001001_b;
  Double I_ERI_Dxy_Px_Py_Py_C1001001001_b = I_ERI_F2xy_S_Py_Py_C1001001001_b+ABX*I_ERI_Dxy_S_Py_Py_C1001001001_b;
  Double I_ERI_Dxz_Px_Py_Py_C1001001001_b = I_ERI_F2xz_S_Py_Py_C1001001001_b+ABX*I_ERI_Dxz_S_Py_Py_C1001001001_b;
  Double I_ERI_D2y_Px_Py_Py_C1001001001_b = I_ERI_Fx2y_S_Py_Py_C1001001001_b+ABX*I_ERI_D2y_S_Py_Py_C1001001001_b;
  Double I_ERI_Dyz_Px_Py_Py_C1001001001_b = I_ERI_Fxyz_S_Py_Py_C1001001001_b+ABX*I_ERI_Dyz_S_Py_Py_C1001001001_b;
  Double I_ERI_D2z_Px_Py_Py_C1001001001_b = I_ERI_Fx2z_S_Py_Py_C1001001001_b+ABX*I_ERI_D2z_S_Py_Py_C1001001001_b;
  Double I_ERI_Dxy_Py_Py_Py_C1001001001_b = I_ERI_Fx2y_S_Py_Py_C1001001001_b+ABY*I_ERI_Dxy_S_Py_Py_C1001001001_b;
  Double I_ERI_Dxz_Py_Py_Py_C1001001001_b = I_ERI_Fxyz_S_Py_Py_C1001001001_b+ABY*I_ERI_Dxz_S_Py_Py_C1001001001_b;
  Double I_ERI_D2y_Py_Py_Py_C1001001001_b = I_ERI_F3y_S_Py_Py_C1001001001_b+ABY*I_ERI_D2y_S_Py_Py_C1001001001_b;
  Double I_ERI_Dyz_Py_Py_Py_C1001001001_b = I_ERI_F2yz_S_Py_Py_C1001001001_b+ABY*I_ERI_Dyz_S_Py_Py_C1001001001_b;
  Double I_ERI_D2z_Py_Py_Py_C1001001001_b = I_ERI_Fy2z_S_Py_Py_C1001001001_b+ABY*I_ERI_D2z_S_Py_Py_C1001001001_b;
  Double I_ERI_Dxz_Pz_Py_Py_C1001001001_b = I_ERI_Fx2z_S_Py_Py_C1001001001_b+ABZ*I_ERI_Dxz_S_Py_Py_C1001001001_b;
  Double I_ERI_Dyz_Pz_Py_Py_C1001001001_b = I_ERI_Fy2z_S_Py_Py_C1001001001_b+ABZ*I_ERI_Dyz_S_Py_Py_C1001001001_b;
  Double I_ERI_D2z_Pz_Py_Py_C1001001001_b = I_ERI_F3z_S_Py_Py_C1001001001_b+ABZ*I_ERI_D2z_S_Py_Py_C1001001001_b;
  Double I_ERI_D2x_Px_Pz_Py_C1001001001_b = I_ERI_F3x_S_Pz_Py_C1001001001_b+ABX*I_ERI_D2x_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dxy_Px_Pz_Py_C1001001001_b = I_ERI_F2xy_S_Pz_Py_C1001001001_b+ABX*I_ERI_Dxy_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dxz_Px_Pz_Py_C1001001001_b = I_ERI_F2xz_S_Pz_Py_C1001001001_b+ABX*I_ERI_Dxz_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2y_Px_Pz_Py_C1001001001_b = I_ERI_Fx2y_S_Pz_Py_C1001001001_b+ABX*I_ERI_D2y_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dyz_Px_Pz_Py_C1001001001_b = I_ERI_Fxyz_S_Pz_Py_C1001001001_b+ABX*I_ERI_Dyz_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2z_Px_Pz_Py_C1001001001_b = I_ERI_Fx2z_S_Pz_Py_C1001001001_b+ABX*I_ERI_D2z_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dxy_Py_Pz_Py_C1001001001_b = I_ERI_Fx2y_S_Pz_Py_C1001001001_b+ABY*I_ERI_Dxy_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dxz_Py_Pz_Py_C1001001001_b = I_ERI_Fxyz_S_Pz_Py_C1001001001_b+ABY*I_ERI_Dxz_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2y_Py_Pz_Py_C1001001001_b = I_ERI_F3y_S_Pz_Py_C1001001001_b+ABY*I_ERI_D2y_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dyz_Py_Pz_Py_C1001001001_b = I_ERI_F2yz_S_Pz_Py_C1001001001_b+ABY*I_ERI_Dyz_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2z_Py_Pz_Py_C1001001001_b = I_ERI_Fy2z_S_Pz_Py_C1001001001_b+ABY*I_ERI_D2z_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dxz_Pz_Pz_Py_C1001001001_b = I_ERI_Fx2z_S_Pz_Py_C1001001001_b+ABZ*I_ERI_Dxz_S_Pz_Py_C1001001001_b;
  Double I_ERI_Dyz_Pz_Pz_Py_C1001001001_b = I_ERI_Fy2z_S_Pz_Py_C1001001001_b+ABZ*I_ERI_Dyz_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2z_Pz_Pz_Py_C1001001001_b = I_ERI_F3z_S_Pz_Py_C1001001001_b+ABZ*I_ERI_D2z_S_Pz_Py_C1001001001_b;
  Double I_ERI_D2x_Px_Px_Pz_C1001001001_b = I_ERI_F3x_S_Px_Pz_C1001001001_b+ABX*I_ERI_D2x_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dxy_Px_Px_Pz_C1001001001_b = I_ERI_F2xy_S_Px_Pz_C1001001001_b+ABX*I_ERI_Dxy_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dxz_Px_Px_Pz_C1001001001_b = I_ERI_F2xz_S_Px_Pz_C1001001001_b+ABX*I_ERI_Dxz_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2y_Px_Px_Pz_C1001001001_b = I_ERI_Fx2y_S_Px_Pz_C1001001001_b+ABX*I_ERI_D2y_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dyz_Px_Px_Pz_C1001001001_b = I_ERI_Fxyz_S_Px_Pz_C1001001001_b+ABX*I_ERI_Dyz_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2z_Px_Px_Pz_C1001001001_b = I_ERI_Fx2z_S_Px_Pz_C1001001001_b+ABX*I_ERI_D2z_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dxy_Py_Px_Pz_C1001001001_b = I_ERI_Fx2y_S_Px_Pz_C1001001001_b+ABY*I_ERI_Dxy_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dxz_Py_Px_Pz_C1001001001_b = I_ERI_Fxyz_S_Px_Pz_C1001001001_b+ABY*I_ERI_Dxz_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2y_Py_Px_Pz_C1001001001_b = I_ERI_F3y_S_Px_Pz_C1001001001_b+ABY*I_ERI_D2y_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dyz_Py_Px_Pz_C1001001001_b = I_ERI_F2yz_S_Px_Pz_C1001001001_b+ABY*I_ERI_Dyz_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2z_Py_Px_Pz_C1001001001_b = I_ERI_Fy2z_S_Px_Pz_C1001001001_b+ABY*I_ERI_D2z_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dxz_Pz_Px_Pz_C1001001001_b = I_ERI_Fx2z_S_Px_Pz_C1001001001_b+ABZ*I_ERI_Dxz_S_Px_Pz_C1001001001_b;
  Double I_ERI_Dyz_Pz_Px_Pz_C1001001001_b = I_ERI_Fy2z_S_Px_Pz_C1001001001_b+ABZ*I_ERI_Dyz_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2z_Pz_Px_Pz_C1001001001_b = I_ERI_F3z_S_Px_Pz_C1001001001_b+ABZ*I_ERI_D2z_S_Px_Pz_C1001001001_b;
  Double I_ERI_D2x_Px_Py_Pz_C1001001001_b = I_ERI_F3x_S_Py_Pz_C1001001001_b+ABX*I_ERI_D2x_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dxy_Px_Py_Pz_C1001001001_b = I_ERI_F2xy_S_Py_Pz_C1001001001_b+ABX*I_ERI_Dxy_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dxz_Px_Py_Pz_C1001001001_b = I_ERI_F2xz_S_Py_Pz_C1001001001_b+ABX*I_ERI_Dxz_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2y_Px_Py_Pz_C1001001001_b = I_ERI_Fx2y_S_Py_Pz_C1001001001_b+ABX*I_ERI_D2y_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dyz_Px_Py_Pz_C1001001001_b = I_ERI_Fxyz_S_Py_Pz_C1001001001_b+ABX*I_ERI_Dyz_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2z_Px_Py_Pz_C1001001001_b = I_ERI_Fx2z_S_Py_Pz_C1001001001_b+ABX*I_ERI_D2z_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dxy_Py_Py_Pz_C1001001001_b = I_ERI_Fx2y_S_Py_Pz_C1001001001_b+ABY*I_ERI_Dxy_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dxz_Py_Py_Pz_C1001001001_b = I_ERI_Fxyz_S_Py_Pz_C1001001001_b+ABY*I_ERI_Dxz_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2y_Py_Py_Pz_C1001001001_b = I_ERI_F3y_S_Py_Pz_C1001001001_b+ABY*I_ERI_D2y_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dyz_Py_Py_Pz_C1001001001_b = I_ERI_F2yz_S_Py_Pz_C1001001001_b+ABY*I_ERI_Dyz_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2z_Py_Py_Pz_C1001001001_b = I_ERI_Fy2z_S_Py_Pz_C1001001001_b+ABY*I_ERI_D2z_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dxz_Pz_Py_Pz_C1001001001_b = I_ERI_Fx2z_S_Py_Pz_C1001001001_b+ABZ*I_ERI_Dxz_S_Py_Pz_C1001001001_b;
  Double I_ERI_Dyz_Pz_Py_Pz_C1001001001_b = I_ERI_Fy2z_S_Py_Pz_C1001001001_b+ABZ*I_ERI_Dyz_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2z_Pz_Py_Pz_C1001001001_b = I_ERI_F3z_S_Py_Pz_C1001001001_b+ABZ*I_ERI_D2z_S_Py_Pz_C1001001001_b;
  Double I_ERI_D2x_Px_Pz_Pz_C1001001001_b = I_ERI_F3x_S_Pz_Pz_C1001001001_b+ABX*I_ERI_D2x_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dxy_Px_Pz_Pz_C1001001001_b = I_ERI_F2xy_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Dxy_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dxz_Px_Pz_Pz_C1001001001_b = I_ERI_F2xz_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Dxz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_D2y_Px_Pz_Pz_C1001001001_b = I_ERI_Fx2y_S_Pz_Pz_C1001001001_b+ABX*I_ERI_D2y_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dyz_Px_Pz_Pz_C1001001001_b = I_ERI_Fxyz_S_Pz_Pz_C1001001001_b+ABX*I_ERI_Dyz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_D2z_Px_Pz_Pz_C1001001001_b = I_ERI_Fx2z_S_Pz_Pz_C1001001001_b+ABX*I_ERI_D2z_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dxy_Py_Pz_Pz_C1001001001_b = I_ERI_Fx2y_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Dxy_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dxz_Py_Pz_Pz_C1001001001_b = I_ERI_Fxyz_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Dxz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_D2y_Py_Pz_Pz_C1001001001_b = I_ERI_F3y_S_Pz_Pz_C1001001001_b+ABY*I_ERI_D2y_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dyz_Py_Pz_Pz_C1001001001_b = I_ERI_F2yz_S_Pz_Pz_C1001001001_b+ABY*I_ERI_Dyz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_D2z_Py_Pz_Pz_C1001001001_b = I_ERI_Fy2z_S_Pz_Pz_C1001001001_b+ABY*I_ERI_D2z_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dxz_Pz_Pz_Pz_C1001001001_b = I_ERI_Fx2z_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_Dxz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_Dyz_Pz_Pz_Pz_C1001001001_b = I_ERI_Fy2z_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_Dyz_S_Pz_Pz_C1001001001_b;
  Double I_ERI_D2z_Pz_Pz_Pz_C1001001001_b = I_ERI_F3z_S_Pz_Pz_C1001001001_b+ABZ*I_ERI_D2z_S_Pz_Pz_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_P_P_C1001001001_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001001_b
   ************************************************************/
  Double I_ERI_Px_D2x_Px_Px_C1001001001_b = I_ERI_D2x_Px_Px_Px_C1001001001_b+ABX*I_ERI_Px_Px_Px_Px_C1001001001_b;
  Double I_ERI_Py_D2x_Px_Px_C1001001001_b = I_ERI_Dxy_Px_Px_Px_C1001001001_b+ABX*I_ERI_Py_Px_Px_Px_C1001001001_b;
  Double I_ERI_Pz_D2x_Px_Px_C1001001001_b = I_ERI_Dxz_Px_Px_Px_C1001001001_b+ABX*I_ERI_Pz_Px_Px_Px_C1001001001_b;
  Double I_ERI_Px_Dxy_Px_Px_C1001001001_b = I_ERI_Dxy_Px_Px_Px_C1001001001_b+ABY*I_ERI_Px_Px_Px_Px_C1001001001_b;
  Double I_ERI_Py_Dxy_Px_Px_C1001001001_b = I_ERI_D2y_Px_Px_Px_C1001001001_b+ABY*I_ERI_Py_Px_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Dxy_Px_Px_C1001001001_b = I_ERI_Dyz_Px_Px_Px_C1001001001_b+ABY*I_ERI_Pz_Px_Px_Px_C1001001001_b;
  Double I_ERI_Px_Dxz_Px_Px_C1001001001_b = I_ERI_Dxz_Px_Px_Px_C1001001001_b+ABZ*I_ERI_Px_Px_Px_Px_C1001001001_b;
  Double I_ERI_Py_Dxz_Px_Px_C1001001001_b = I_ERI_Dyz_Px_Px_Px_C1001001001_b+ABZ*I_ERI_Py_Px_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Dxz_Px_Px_C1001001001_b = I_ERI_D2z_Px_Px_Px_C1001001001_b+ABZ*I_ERI_Pz_Px_Px_Px_C1001001001_b;
  Double I_ERI_Px_D2y_Px_Px_C1001001001_b = I_ERI_Dxy_Py_Px_Px_C1001001001_b+ABY*I_ERI_Px_Py_Px_Px_C1001001001_b;
  Double I_ERI_Py_D2y_Px_Px_C1001001001_b = I_ERI_D2y_Py_Px_Px_C1001001001_b+ABY*I_ERI_Py_Py_Px_Px_C1001001001_b;
  Double I_ERI_Pz_D2y_Px_Px_C1001001001_b = I_ERI_Dyz_Py_Px_Px_C1001001001_b+ABY*I_ERI_Pz_Py_Px_Px_C1001001001_b;
  Double I_ERI_Px_Dyz_Px_Px_C1001001001_b = I_ERI_Dxz_Py_Px_Px_C1001001001_b+ABZ*I_ERI_Px_Py_Px_Px_C1001001001_b;
  Double I_ERI_Py_Dyz_Px_Px_C1001001001_b = I_ERI_Dyz_Py_Px_Px_C1001001001_b+ABZ*I_ERI_Py_Py_Px_Px_C1001001001_b;
  Double I_ERI_Pz_Dyz_Px_Px_C1001001001_b = I_ERI_D2z_Py_Px_Px_C1001001001_b+ABZ*I_ERI_Pz_Py_Px_Px_C1001001001_b;
  Double I_ERI_Px_D2z_Px_Px_C1001001001_b = I_ERI_Dxz_Pz_Px_Px_C1001001001_b+ABZ*I_ERI_Px_Pz_Px_Px_C1001001001_b;
  Double I_ERI_Py_D2z_Px_Px_C1001001001_b = I_ERI_Dyz_Pz_Px_Px_C1001001001_b+ABZ*I_ERI_Py_Pz_Px_Px_C1001001001_b;
  Double I_ERI_Pz_D2z_Px_Px_C1001001001_b = I_ERI_D2z_Pz_Px_Px_C1001001001_b+ABZ*I_ERI_Pz_Pz_Px_Px_C1001001001_b;
  Double I_ERI_Px_D2x_Py_Px_C1001001001_b = I_ERI_D2x_Px_Py_Px_C1001001001_b+ABX*I_ERI_Px_Px_Py_Px_C1001001001_b;
  Double I_ERI_Py_D2x_Py_Px_C1001001001_b = I_ERI_Dxy_Px_Py_Px_C1001001001_b+ABX*I_ERI_Py_Px_Py_Px_C1001001001_b;
  Double I_ERI_Pz_D2x_Py_Px_C1001001001_b = I_ERI_Dxz_Px_Py_Px_C1001001001_b+ABX*I_ERI_Pz_Px_Py_Px_C1001001001_b;
  Double I_ERI_Px_Dxy_Py_Px_C1001001001_b = I_ERI_Dxy_Px_Py_Px_C1001001001_b+ABY*I_ERI_Px_Px_Py_Px_C1001001001_b;
  Double I_ERI_Py_Dxy_Py_Px_C1001001001_b = I_ERI_D2y_Px_Py_Px_C1001001001_b+ABY*I_ERI_Py_Px_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Dxy_Py_Px_C1001001001_b = I_ERI_Dyz_Px_Py_Px_C1001001001_b+ABY*I_ERI_Pz_Px_Py_Px_C1001001001_b;
  Double I_ERI_Px_Dxz_Py_Px_C1001001001_b = I_ERI_Dxz_Px_Py_Px_C1001001001_b+ABZ*I_ERI_Px_Px_Py_Px_C1001001001_b;
  Double I_ERI_Py_Dxz_Py_Px_C1001001001_b = I_ERI_Dyz_Px_Py_Px_C1001001001_b+ABZ*I_ERI_Py_Px_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Dxz_Py_Px_C1001001001_b = I_ERI_D2z_Px_Py_Px_C1001001001_b+ABZ*I_ERI_Pz_Px_Py_Px_C1001001001_b;
  Double I_ERI_Px_D2y_Py_Px_C1001001001_b = I_ERI_Dxy_Py_Py_Px_C1001001001_b+ABY*I_ERI_Px_Py_Py_Px_C1001001001_b;
  Double I_ERI_Py_D2y_Py_Px_C1001001001_b = I_ERI_D2y_Py_Py_Px_C1001001001_b+ABY*I_ERI_Py_Py_Py_Px_C1001001001_b;
  Double I_ERI_Pz_D2y_Py_Px_C1001001001_b = I_ERI_Dyz_Py_Py_Px_C1001001001_b+ABY*I_ERI_Pz_Py_Py_Px_C1001001001_b;
  Double I_ERI_Px_Dyz_Py_Px_C1001001001_b = I_ERI_Dxz_Py_Py_Px_C1001001001_b+ABZ*I_ERI_Px_Py_Py_Px_C1001001001_b;
  Double I_ERI_Py_Dyz_Py_Px_C1001001001_b = I_ERI_Dyz_Py_Py_Px_C1001001001_b+ABZ*I_ERI_Py_Py_Py_Px_C1001001001_b;
  Double I_ERI_Pz_Dyz_Py_Px_C1001001001_b = I_ERI_D2z_Py_Py_Px_C1001001001_b+ABZ*I_ERI_Pz_Py_Py_Px_C1001001001_b;
  Double I_ERI_Px_D2z_Py_Px_C1001001001_b = I_ERI_Dxz_Pz_Py_Px_C1001001001_b+ABZ*I_ERI_Px_Pz_Py_Px_C1001001001_b;
  Double I_ERI_Py_D2z_Py_Px_C1001001001_b = I_ERI_Dyz_Pz_Py_Px_C1001001001_b+ABZ*I_ERI_Py_Pz_Py_Px_C1001001001_b;
  Double I_ERI_Pz_D2z_Py_Px_C1001001001_b = I_ERI_D2z_Pz_Py_Px_C1001001001_b+ABZ*I_ERI_Pz_Pz_Py_Px_C1001001001_b;
  Double I_ERI_Px_D2x_Pz_Px_C1001001001_b = I_ERI_D2x_Px_Pz_Px_C1001001001_b+ABX*I_ERI_Px_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Py_D2x_Pz_Px_C1001001001_b = I_ERI_Dxy_Px_Pz_Px_C1001001001_b+ABX*I_ERI_Py_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_D2x_Pz_Px_C1001001001_b = I_ERI_Dxz_Px_Pz_Px_C1001001001_b+ABX*I_ERI_Pz_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Dxy_Pz_Px_C1001001001_b = I_ERI_Dxy_Px_Pz_Px_C1001001001_b+ABY*I_ERI_Px_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Dxy_Pz_Px_C1001001001_b = I_ERI_D2y_Px_Pz_Px_C1001001001_b+ABY*I_ERI_Py_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Dxy_Pz_Px_C1001001001_b = I_ERI_Dyz_Px_Pz_Px_C1001001001_b+ABY*I_ERI_Pz_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Dxz_Pz_Px_C1001001001_b = I_ERI_Dxz_Px_Pz_Px_C1001001001_b+ABZ*I_ERI_Px_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Dxz_Pz_Px_C1001001001_b = I_ERI_Dyz_Px_Pz_Px_C1001001001_b+ABZ*I_ERI_Py_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Dxz_Pz_Px_C1001001001_b = I_ERI_D2z_Px_Pz_Px_C1001001001_b+ABZ*I_ERI_Pz_Px_Pz_Px_C1001001001_b;
  Double I_ERI_Px_D2y_Pz_Px_C1001001001_b = I_ERI_Dxy_Py_Pz_Px_C1001001001_b+ABY*I_ERI_Px_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Py_D2y_Pz_Px_C1001001001_b = I_ERI_D2y_Py_Pz_Px_C1001001001_b+ABY*I_ERI_Py_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_D2y_Pz_Px_C1001001001_b = I_ERI_Dyz_Py_Pz_Px_C1001001001_b+ABY*I_ERI_Pz_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Px_Dyz_Pz_Px_C1001001001_b = I_ERI_Dxz_Py_Pz_Px_C1001001001_b+ABZ*I_ERI_Px_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Py_Dyz_Pz_Px_C1001001001_b = I_ERI_Dyz_Py_Pz_Px_C1001001001_b+ABZ*I_ERI_Py_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_Dyz_Pz_Px_C1001001001_b = I_ERI_D2z_Py_Pz_Px_C1001001001_b+ABZ*I_ERI_Pz_Py_Pz_Px_C1001001001_b;
  Double I_ERI_Px_D2z_Pz_Px_C1001001001_b = I_ERI_Dxz_Pz_Pz_Px_C1001001001_b+ABZ*I_ERI_Px_Pz_Pz_Px_C1001001001_b;
  Double I_ERI_Py_D2z_Pz_Px_C1001001001_b = I_ERI_Dyz_Pz_Pz_Px_C1001001001_b+ABZ*I_ERI_Py_Pz_Pz_Px_C1001001001_b;
  Double I_ERI_Pz_D2z_Pz_Px_C1001001001_b = I_ERI_D2z_Pz_Pz_Px_C1001001001_b+ABZ*I_ERI_Pz_Pz_Pz_Px_C1001001001_b;
  Double I_ERI_Px_D2x_Px_Py_C1001001001_b = I_ERI_D2x_Px_Px_Py_C1001001001_b+ABX*I_ERI_Px_Px_Px_Py_C1001001001_b;
  Double I_ERI_Py_D2x_Px_Py_C1001001001_b = I_ERI_Dxy_Px_Px_Py_C1001001001_b+ABX*I_ERI_Py_Px_Px_Py_C1001001001_b;
  Double I_ERI_Pz_D2x_Px_Py_C1001001001_b = I_ERI_Dxz_Px_Px_Py_C1001001001_b+ABX*I_ERI_Pz_Px_Px_Py_C1001001001_b;
  Double I_ERI_Px_Dxy_Px_Py_C1001001001_b = I_ERI_Dxy_Px_Px_Py_C1001001001_b+ABY*I_ERI_Px_Px_Px_Py_C1001001001_b;
  Double I_ERI_Py_Dxy_Px_Py_C1001001001_b = I_ERI_D2y_Px_Px_Py_C1001001001_b+ABY*I_ERI_Py_Px_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Dxy_Px_Py_C1001001001_b = I_ERI_Dyz_Px_Px_Py_C1001001001_b+ABY*I_ERI_Pz_Px_Px_Py_C1001001001_b;
  Double I_ERI_Px_Dxz_Px_Py_C1001001001_b = I_ERI_Dxz_Px_Px_Py_C1001001001_b+ABZ*I_ERI_Px_Px_Px_Py_C1001001001_b;
  Double I_ERI_Py_Dxz_Px_Py_C1001001001_b = I_ERI_Dyz_Px_Px_Py_C1001001001_b+ABZ*I_ERI_Py_Px_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Dxz_Px_Py_C1001001001_b = I_ERI_D2z_Px_Px_Py_C1001001001_b+ABZ*I_ERI_Pz_Px_Px_Py_C1001001001_b;
  Double I_ERI_Px_D2y_Px_Py_C1001001001_b = I_ERI_Dxy_Py_Px_Py_C1001001001_b+ABY*I_ERI_Px_Py_Px_Py_C1001001001_b;
  Double I_ERI_Py_D2y_Px_Py_C1001001001_b = I_ERI_D2y_Py_Px_Py_C1001001001_b+ABY*I_ERI_Py_Py_Px_Py_C1001001001_b;
  Double I_ERI_Pz_D2y_Px_Py_C1001001001_b = I_ERI_Dyz_Py_Px_Py_C1001001001_b+ABY*I_ERI_Pz_Py_Px_Py_C1001001001_b;
  Double I_ERI_Px_Dyz_Px_Py_C1001001001_b = I_ERI_Dxz_Py_Px_Py_C1001001001_b+ABZ*I_ERI_Px_Py_Px_Py_C1001001001_b;
  Double I_ERI_Py_Dyz_Px_Py_C1001001001_b = I_ERI_Dyz_Py_Px_Py_C1001001001_b+ABZ*I_ERI_Py_Py_Px_Py_C1001001001_b;
  Double I_ERI_Pz_Dyz_Px_Py_C1001001001_b = I_ERI_D2z_Py_Px_Py_C1001001001_b+ABZ*I_ERI_Pz_Py_Px_Py_C1001001001_b;
  Double I_ERI_Px_D2z_Px_Py_C1001001001_b = I_ERI_Dxz_Pz_Px_Py_C1001001001_b+ABZ*I_ERI_Px_Pz_Px_Py_C1001001001_b;
  Double I_ERI_Py_D2z_Px_Py_C1001001001_b = I_ERI_Dyz_Pz_Px_Py_C1001001001_b+ABZ*I_ERI_Py_Pz_Px_Py_C1001001001_b;
  Double I_ERI_Pz_D2z_Px_Py_C1001001001_b = I_ERI_D2z_Pz_Px_Py_C1001001001_b+ABZ*I_ERI_Pz_Pz_Px_Py_C1001001001_b;
  Double I_ERI_Px_D2x_Py_Py_C1001001001_b = I_ERI_D2x_Px_Py_Py_C1001001001_b+ABX*I_ERI_Px_Px_Py_Py_C1001001001_b;
  Double I_ERI_Py_D2x_Py_Py_C1001001001_b = I_ERI_Dxy_Px_Py_Py_C1001001001_b+ABX*I_ERI_Py_Px_Py_Py_C1001001001_b;
  Double I_ERI_Pz_D2x_Py_Py_C1001001001_b = I_ERI_Dxz_Px_Py_Py_C1001001001_b+ABX*I_ERI_Pz_Px_Py_Py_C1001001001_b;
  Double I_ERI_Px_Dxy_Py_Py_C1001001001_b = I_ERI_Dxy_Px_Py_Py_C1001001001_b+ABY*I_ERI_Px_Px_Py_Py_C1001001001_b;
  Double I_ERI_Py_Dxy_Py_Py_C1001001001_b = I_ERI_D2y_Px_Py_Py_C1001001001_b+ABY*I_ERI_Py_Px_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Dxy_Py_Py_C1001001001_b = I_ERI_Dyz_Px_Py_Py_C1001001001_b+ABY*I_ERI_Pz_Px_Py_Py_C1001001001_b;
  Double I_ERI_Px_Dxz_Py_Py_C1001001001_b = I_ERI_Dxz_Px_Py_Py_C1001001001_b+ABZ*I_ERI_Px_Px_Py_Py_C1001001001_b;
  Double I_ERI_Py_Dxz_Py_Py_C1001001001_b = I_ERI_Dyz_Px_Py_Py_C1001001001_b+ABZ*I_ERI_Py_Px_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Dxz_Py_Py_C1001001001_b = I_ERI_D2z_Px_Py_Py_C1001001001_b+ABZ*I_ERI_Pz_Px_Py_Py_C1001001001_b;
  Double I_ERI_Px_D2y_Py_Py_C1001001001_b = I_ERI_Dxy_Py_Py_Py_C1001001001_b+ABY*I_ERI_Px_Py_Py_Py_C1001001001_b;
  Double I_ERI_Py_D2y_Py_Py_C1001001001_b = I_ERI_D2y_Py_Py_Py_C1001001001_b+ABY*I_ERI_Py_Py_Py_Py_C1001001001_b;
  Double I_ERI_Pz_D2y_Py_Py_C1001001001_b = I_ERI_Dyz_Py_Py_Py_C1001001001_b+ABY*I_ERI_Pz_Py_Py_Py_C1001001001_b;
  Double I_ERI_Px_Dyz_Py_Py_C1001001001_b = I_ERI_Dxz_Py_Py_Py_C1001001001_b+ABZ*I_ERI_Px_Py_Py_Py_C1001001001_b;
  Double I_ERI_Py_Dyz_Py_Py_C1001001001_b = I_ERI_Dyz_Py_Py_Py_C1001001001_b+ABZ*I_ERI_Py_Py_Py_Py_C1001001001_b;
  Double I_ERI_Pz_Dyz_Py_Py_C1001001001_b = I_ERI_D2z_Py_Py_Py_C1001001001_b+ABZ*I_ERI_Pz_Py_Py_Py_C1001001001_b;
  Double I_ERI_Px_D2z_Py_Py_C1001001001_b = I_ERI_Dxz_Pz_Py_Py_C1001001001_b+ABZ*I_ERI_Px_Pz_Py_Py_C1001001001_b;
  Double I_ERI_Py_D2z_Py_Py_C1001001001_b = I_ERI_Dyz_Pz_Py_Py_C1001001001_b+ABZ*I_ERI_Py_Pz_Py_Py_C1001001001_b;
  Double I_ERI_Pz_D2z_Py_Py_C1001001001_b = I_ERI_D2z_Pz_Py_Py_C1001001001_b+ABZ*I_ERI_Pz_Pz_Py_Py_C1001001001_b;
  Double I_ERI_Px_D2x_Pz_Py_C1001001001_b = I_ERI_D2x_Px_Pz_Py_C1001001001_b+ABX*I_ERI_Px_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Py_D2x_Pz_Py_C1001001001_b = I_ERI_Dxy_Px_Pz_Py_C1001001001_b+ABX*I_ERI_Py_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_D2x_Pz_Py_C1001001001_b = I_ERI_Dxz_Px_Pz_Py_C1001001001_b+ABX*I_ERI_Pz_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Dxy_Pz_Py_C1001001001_b = I_ERI_Dxy_Px_Pz_Py_C1001001001_b+ABY*I_ERI_Px_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Dxy_Pz_Py_C1001001001_b = I_ERI_D2y_Px_Pz_Py_C1001001001_b+ABY*I_ERI_Py_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Dxy_Pz_Py_C1001001001_b = I_ERI_Dyz_Px_Pz_Py_C1001001001_b+ABY*I_ERI_Pz_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Dxz_Pz_Py_C1001001001_b = I_ERI_Dxz_Px_Pz_Py_C1001001001_b+ABZ*I_ERI_Px_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Dxz_Pz_Py_C1001001001_b = I_ERI_Dyz_Px_Pz_Py_C1001001001_b+ABZ*I_ERI_Py_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Dxz_Pz_Py_C1001001001_b = I_ERI_D2z_Px_Pz_Py_C1001001001_b+ABZ*I_ERI_Pz_Px_Pz_Py_C1001001001_b;
  Double I_ERI_Px_D2y_Pz_Py_C1001001001_b = I_ERI_Dxy_Py_Pz_Py_C1001001001_b+ABY*I_ERI_Px_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Py_D2y_Pz_Py_C1001001001_b = I_ERI_D2y_Py_Pz_Py_C1001001001_b+ABY*I_ERI_Py_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_D2y_Pz_Py_C1001001001_b = I_ERI_Dyz_Py_Pz_Py_C1001001001_b+ABY*I_ERI_Pz_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Px_Dyz_Pz_Py_C1001001001_b = I_ERI_Dxz_Py_Pz_Py_C1001001001_b+ABZ*I_ERI_Px_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Py_Dyz_Pz_Py_C1001001001_b = I_ERI_Dyz_Py_Pz_Py_C1001001001_b+ABZ*I_ERI_Py_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_Dyz_Pz_Py_C1001001001_b = I_ERI_D2z_Py_Pz_Py_C1001001001_b+ABZ*I_ERI_Pz_Py_Pz_Py_C1001001001_b;
  Double I_ERI_Px_D2z_Pz_Py_C1001001001_b = I_ERI_Dxz_Pz_Pz_Py_C1001001001_b+ABZ*I_ERI_Px_Pz_Pz_Py_C1001001001_b;
  Double I_ERI_Py_D2z_Pz_Py_C1001001001_b = I_ERI_Dyz_Pz_Pz_Py_C1001001001_b+ABZ*I_ERI_Py_Pz_Pz_Py_C1001001001_b;
  Double I_ERI_Pz_D2z_Pz_Py_C1001001001_b = I_ERI_D2z_Pz_Pz_Py_C1001001001_b+ABZ*I_ERI_Pz_Pz_Pz_Py_C1001001001_b;
  Double I_ERI_Px_D2x_Px_Pz_C1001001001_b = I_ERI_D2x_Px_Px_Pz_C1001001001_b+ABX*I_ERI_Px_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Py_D2x_Px_Pz_C1001001001_b = I_ERI_Dxy_Px_Px_Pz_C1001001001_b+ABX*I_ERI_Py_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_D2x_Px_Pz_C1001001001_b = I_ERI_Dxz_Px_Px_Pz_C1001001001_b+ABX*I_ERI_Pz_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Dxy_Px_Pz_C1001001001_b = I_ERI_Dxy_Px_Px_Pz_C1001001001_b+ABY*I_ERI_Px_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Dxy_Px_Pz_C1001001001_b = I_ERI_D2y_Px_Px_Pz_C1001001001_b+ABY*I_ERI_Py_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxy_Px_Pz_C1001001001_b = I_ERI_Dyz_Px_Px_Pz_C1001001001_b+ABY*I_ERI_Pz_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Dxz_Px_Pz_C1001001001_b = I_ERI_Dxz_Px_Px_Pz_C1001001001_b+ABZ*I_ERI_Px_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Dxz_Px_Pz_C1001001001_b = I_ERI_Dyz_Px_Px_Pz_C1001001001_b+ABZ*I_ERI_Py_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxz_Px_Pz_C1001001001_b = I_ERI_D2z_Px_Px_Pz_C1001001001_b+ABZ*I_ERI_Pz_Px_Px_Pz_C1001001001_b;
  Double I_ERI_Px_D2y_Px_Pz_C1001001001_b = I_ERI_Dxy_Py_Px_Pz_C1001001001_b+ABY*I_ERI_Px_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Py_D2y_Px_Pz_C1001001001_b = I_ERI_D2y_Py_Px_Pz_C1001001001_b+ABY*I_ERI_Py_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_D2y_Px_Pz_C1001001001_b = I_ERI_Dyz_Py_Px_Pz_C1001001001_b+ABY*I_ERI_Pz_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Px_Dyz_Px_Pz_C1001001001_b = I_ERI_Dxz_Py_Px_Pz_C1001001001_b+ABZ*I_ERI_Px_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Py_Dyz_Px_Pz_C1001001001_b = I_ERI_Dyz_Py_Px_Pz_C1001001001_b+ABZ*I_ERI_Py_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_Dyz_Px_Pz_C1001001001_b = I_ERI_D2z_Py_Px_Pz_C1001001001_b+ABZ*I_ERI_Pz_Py_Px_Pz_C1001001001_b;
  Double I_ERI_Px_D2z_Px_Pz_C1001001001_b = I_ERI_Dxz_Pz_Px_Pz_C1001001001_b+ABZ*I_ERI_Px_Pz_Px_Pz_C1001001001_b;
  Double I_ERI_Py_D2z_Px_Pz_C1001001001_b = I_ERI_Dyz_Pz_Px_Pz_C1001001001_b+ABZ*I_ERI_Py_Pz_Px_Pz_C1001001001_b;
  Double I_ERI_Pz_D2z_Px_Pz_C1001001001_b = I_ERI_D2z_Pz_Px_Pz_C1001001001_b+ABZ*I_ERI_Pz_Pz_Px_Pz_C1001001001_b;
  Double I_ERI_Px_D2x_Py_Pz_C1001001001_b = I_ERI_D2x_Px_Py_Pz_C1001001001_b+ABX*I_ERI_Px_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Py_D2x_Py_Pz_C1001001001_b = I_ERI_Dxy_Px_Py_Pz_C1001001001_b+ABX*I_ERI_Py_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_D2x_Py_Pz_C1001001001_b = I_ERI_Dxz_Px_Py_Pz_C1001001001_b+ABX*I_ERI_Pz_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Dxy_Py_Pz_C1001001001_b = I_ERI_Dxy_Px_Py_Pz_C1001001001_b+ABY*I_ERI_Px_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Dxy_Py_Pz_C1001001001_b = I_ERI_D2y_Px_Py_Pz_C1001001001_b+ABY*I_ERI_Py_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxy_Py_Pz_C1001001001_b = I_ERI_Dyz_Px_Py_Pz_C1001001001_b+ABY*I_ERI_Pz_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Dxz_Py_Pz_C1001001001_b = I_ERI_Dxz_Px_Py_Pz_C1001001001_b+ABZ*I_ERI_Px_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Dxz_Py_Pz_C1001001001_b = I_ERI_Dyz_Px_Py_Pz_C1001001001_b+ABZ*I_ERI_Py_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxz_Py_Pz_C1001001001_b = I_ERI_D2z_Px_Py_Pz_C1001001001_b+ABZ*I_ERI_Pz_Px_Py_Pz_C1001001001_b;
  Double I_ERI_Px_D2y_Py_Pz_C1001001001_b = I_ERI_Dxy_Py_Py_Pz_C1001001001_b+ABY*I_ERI_Px_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Py_D2y_Py_Pz_C1001001001_b = I_ERI_D2y_Py_Py_Pz_C1001001001_b+ABY*I_ERI_Py_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_D2y_Py_Pz_C1001001001_b = I_ERI_Dyz_Py_Py_Pz_C1001001001_b+ABY*I_ERI_Pz_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Px_Dyz_Py_Pz_C1001001001_b = I_ERI_Dxz_Py_Py_Pz_C1001001001_b+ABZ*I_ERI_Px_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Py_Dyz_Py_Pz_C1001001001_b = I_ERI_Dyz_Py_Py_Pz_C1001001001_b+ABZ*I_ERI_Py_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_Dyz_Py_Pz_C1001001001_b = I_ERI_D2z_Py_Py_Pz_C1001001001_b+ABZ*I_ERI_Pz_Py_Py_Pz_C1001001001_b;
  Double I_ERI_Px_D2z_Py_Pz_C1001001001_b = I_ERI_Dxz_Pz_Py_Pz_C1001001001_b+ABZ*I_ERI_Px_Pz_Py_Pz_C1001001001_b;
  Double I_ERI_Py_D2z_Py_Pz_C1001001001_b = I_ERI_Dyz_Pz_Py_Pz_C1001001001_b+ABZ*I_ERI_Py_Pz_Py_Pz_C1001001001_b;
  Double I_ERI_Pz_D2z_Py_Pz_C1001001001_b = I_ERI_D2z_Pz_Py_Pz_C1001001001_b+ABZ*I_ERI_Pz_Pz_Py_Pz_C1001001001_b;
  Double I_ERI_Px_D2x_Pz_Pz_C1001001001_b = I_ERI_D2x_Px_Pz_Pz_C1001001001_b+ABX*I_ERI_Px_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_D2x_Pz_Pz_C1001001001_b = I_ERI_Dxy_Px_Pz_Pz_C1001001001_b+ABX*I_ERI_Py_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_D2x_Pz_Pz_C1001001001_b = I_ERI_Dxz_Px_Pz_Pz_C1001001001_b+ABX*I_ERI_Pz_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_Dxy_Pz_Pz_C1001001001_b = I_ERI_Dxy_Px_Pz_Pz_C1001001001_b+ABY*I_ERI_Px_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Dxy_Pz_Pz_C1001001001_b = I_ERI_D2y_Px_Pz_Pz_C1001001001_b+ABY*I_ERI_Py_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxy_Pz_Pz_C1001001001_b = I_ERI_Dyz_Px_Pz_Pz_C1001001001_b+ABY*I_ERI_Pz_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_Dxz_Pz_Pz_C1001001001_b = I_ERI_Dxz_Px_Pz_Pz_C1001001001_b+ABZ*I_ERI_Px_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Dxz_Pz_Pz_C1001001001_b = I_ERI_Dyz_Px_Pz_Pz_C1001001001_b+ABZ*I_ERI_Py_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Dxz_Pz_Pz_C1001001001_b = I_ERI_D2z_Px_Pz_Pz_C1001001001_b+ABZ*I_ERI_Pz_Px_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_D2y_Pz_Pz_C1001001001_b = I_ERI_Dxy_Py_Pz_Pz_C1001001001_b+ABY*I_ERI_Px_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_D2y_Pz_Pz_C1001001001_b = I_ERI_D2y_Py_Pz_Pz_C1001001001_b+ABY*I_ERI_Py_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_D2y_Pz_Pz_C1001001001_b = I_ERI_Dyz_Py_Pz_Pz_C1001001001_b+ABY*I_ERI_Pz_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_Dyz_Pz_Pz_C1001001001_b = I_ERI_Dxz_Py_Pz_Pz_C1001001001_b+ABZ*I_ERI_Px_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_Dyz_Pz_Pz_C1001001001_b = I_ERI_Dyz_Py_Pz_Pz_C1001001001_b+ABZ*I_ERI_Py_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_Dyz_Pz_Pz_C1001001001_b = I_ERI_D2z_Py_Pz_Pz_C1001001001_b+ABZ*I_ERI_Pz_Py_Pz_Pz_C1001001001_b;
  Double I_ERI_Px_D2z_Pz_Pz_C1001001001_b = I_ERI_Dxz_Pz_Pz_Pz_C1001001001_b+ABZ*I_ERI_Px_Pz_Pz_Pz_C1001001001_b;
  Double I_ERI_Py_D2z_Pz_Pz_C1001001001_b = I_ERI_Dyz_Pz_Pz_Pz_C1001001001_b+ABZ*I_ERI_Py_Pz_Pz_Pz_C1001001001_b;
  Double I_ERI_Pz_D2z_Pz_Pz_C1001001001_b = I_ERI_D2z_Pz_Pz_Pz_C1001001001_b+ABZ*I_ERI_Pz_Pz_Pz_Pz_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001_c
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001_c = I_ERI_D2x_S_Px_S_C1001_c+ABX*I_ERI_Px_S_Px_S_C1001_c;
  Double I_ERI_Py_Px_Px_S_C1001_c = I_ERI_Dxy_S_Px_S_C1001_c+ABX*I_ERI_Py_S_Px_S_C1001_c;
  Double I_ERI_Pz_Px_Px_S_C1001_c = I_ERI_Dxz_S_Px_S_C1001_c+ABX*I_ERI_Pz_S_Px_S_C1001_c;
  Double I_ERI_Px_Py_Px_S_C1001_c = I_ERI_Dxy_S_Px_S_C1001_c+ABY*I_ERI_Px_S_Px_S_C1001_c;
  Double I_ERI_Py_Py_Px_S_C1001_c = I_ERI_D2y_S_Px_S_C1001_c+ABY*I_ERI_Py_S_Px_S_C1001_c;
  Double I_ERI_Pz_Py_Px_S_C1001_c = I_ERI_Dyz_S_Px_S_C1001_c+ABY*I_ERI_Pz_S_Px_S_C1001_c;
  Double I_ERI_Px_Pz_Px_S_C1001_c = I_ERI_Dxz_S_Px_S_C1001_c+ABZ*I_ERI_Px_S_Px_S_C1001_c;
  Double I_ERI_Py_Pz_Px_S_C1001_c = I_ERI_Dyz_S_Px_S_C1001_c+ABZ*I_ERI_Py_S_Px_S_C1001_c;
  Double I_ERI_Pz_Pz_Px_S_C1001_c = I_ERI_D2z_S_Px_S_C1001_c+ABZ*I_ERI_Pz_S_Px_S_C1001_c;
  Double I_ERI_Px_Px_Py_S_C1001_c = I_ERI_D2x_S_Py_S_C1001_c+ABX*I_ERI_Px_S_Py_S_C1001_c;
  Double I_ERI_Py_Px_Py_S_C1001_c = I_ERI_Dxy_S_Py_S_C1001_c+ABX*I_ERI_Py_S_Py_S_C1001_c;
  Double I_ERI_Pz_Px_Py_S_C1001_c = I_ERI_Dxz_S_Py_S_C1001_c+ABX*I_ERI_Pz_S_Py_S_C1001_c;
  Double I_ERI_Px_Py_Py_S_C1001_c = I_ERI_Dxy_S_Py_S_C1001_c+ABY*I_ERI_Px_S_Py_S_C1001_c;
  Double I_ERI_Py_Py_Py_S_C1001_c = I_ERI_D2y_S_Py_S_C1001_c+ABY*I_ERI_Py_S_Py_S_C1001_c;
  Double I_ERI_Pz_Py_Py_S_C1001_c = I_ERI_Dyz_S_Py_S_C1001_c+ABY*I_ERI_Pz_S_Py_S_C1001_c;
  Double I_ERI_Px_Pz_Py_S_C1001_c = I_ERI_Dxz_S_Py_S_C1001_c+ABZ*I_ERI_Px_S_Py_S_C1001_c;
  Double I_ERI_Py_Pz_Py_S_C1001_c = I_ERI_Dyz_S_Py_S_C1001_c+ABZ*I_ERI_Py_S_Py_S_C1001_c;
  Double I_ERI_Pz_Pz_Py_S_C1001_c = I_ERI_D2z_S_Py_S_C1001_c+ABZ*I_ERI_Pz_S_Py_S_C1001_c;
  Double I_ERI_Px_Px_Pz_S_C1001_c = I_ERI_D2x_S_Pz_S_C1001_c+ABX*I_ERI_Px_S_Pz_S_C1001_c;
  Double I_ERI_Py_Px_Pz_S_C1001_c = I_ERI_Dxy_S_Pz_S_C1001_c+ABX*I_ERI_Py_S_Pz_S_C1001_c;
  Double I_ERI_Pz_Px_Pz_S_C1001_c = I_ERI_Dxz_S_Pz_S_C1001_c+ABX*I_ERI_Pz_S_Pz_S_C1001_c;
  Double I_ERI_Px_Py_Pz_S_C1001_c = I_ERI_Dxy_S_Pz_S_C1001_c+ABY*I_ERI_Px_S_Pz_S_C1001_c;
  Double I_ERI_Py_Py_Pz_S_C1001_c = I_ERI_D2y_S_Pz_S_C1001_c+ABY*I_ERI_Py_S_Pz_S_C1001_c;
  Double I_ERI_Pz_Py_Pz_S_C1001_c = I_ERI_Dyz_S_Pz_S_C1001_c+ABY*I_ERI_Pz_S_Pz_S_C1001_c;
  Double I_ERI_Px_Pz_Pz_S_C1001_c = I_ERI_Dxz_S_Pz_S_C1001_c+ABZ*I_ERI_Px_S_Pz_S_C1001_c;
  Double I_ERI_Py_Pz_Pz_S_C1001_c = I_ERI_Dyz_S_Pz_S_C1001_c+ABZ*I_ERI_Py_S_Pz_S_C1001_c;
  Double I_ERI_Pz_Pz_Pz_S_C1001_c = I_ERI_D2z_S_Pz_S_C1001_c+ABZ*I_ERI_Pz_S_Pz_S_C1001_c;

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
   * shell quartet name: SQ_ERI_P_P_P_P_C1000001001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1000001001_c
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000001001_c
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_C1000001001_c = I_ERI_D2x_S_Px_Px_C1000001001_c+ABX*I_ERI_Px_S_Px_Px_C1000001001_c;
  Double I_ERI_Py_Px_Px_Px_C1000001001_c = I_ERI_Dxy_S_Px_Px_C1000001001_c+ABX*I_ERI_Py_S_Px_Px_C1000001001_c;
  Double I_ERI_Pz_Px_Px_Px_C1000001001_c = I_ERI_Dxz_S_Px_Px_C1000001001_c+ABX*I_ERI_Pz_S_Px_Px_C1000001001_c;
  Double I_ERI_Px_Py_Px_Px_C1000001001_c = I_ERI_Dxy_S_Px_Px_C1000001001_c+ABY*I_ERI_Px_S_Px_Px_C1000001001_c;
  Double I_ERI_Py_Py_Px_Px_C1000001001_c = I_ERI_D2y_S_Px_Px_C1000001001_c+ABY*I_ERI_Py_S_Px_Px_C1000001001_c;
  Double I_ERI_Pz_Py_Px_Px_C1000001001_c = I_ERI_Dyz_S_Px_Px_C1000001001_c+ABY*I_ERI_Pz_S_Px_Px_C1000001001_c;
  Double I_ERI_Px_Pz_Px_Px_C1000001001_c = I_ERI_Dxz_S_Px_Px_C1000001001_c+ABZ*I_ERI_Px_S_Px_Px_C1000001001_c;
  Double I_ERI_Py_Pz_Px_Px_C1000001001_c = I_ERI_Dyz_S_Px_Px_C1000001001_c+ABZ*I_ERI_Py_S_Px_Px_C1000001001_c;
  Double I_ERI_Pz_Pz_Px_Px_C1000001001_c = I_ERI_D2z_S_Px_Px_C1000001001_c+ABZ*I_ERI_Pz_S_Px_Px_C1000001001_c;
  Double I_ERI_Px_Px_Py_Px_C1000001001_c = I_ERI_D2x_S_Py_Px_C1000001001_c+ABX*I_ERI_Px_S_Py_Px_C1000001001_c;
  Double I_ERI_Py_Px_Py_Px_C1000001001_c = I_ERI_Dxy_S_Py_Px_C1000001001_c+ABX*I_ERI_Py_S_Py_Px_C1000001001_c;
  Double I_ERI_Pz_Px_Py_Px_C1000001001_c = I_ERI_Dxz_S_Py_Px_C1000001001_c+ABX*I_ERI_Pz_S_Py_Px_C1000001001_c;
  Double I_ERI_Px_Py_Py_Px_C1000001001_c = I_ERI_Dxy_S_Py_Px_C1000001001_c+ABY*I_ERI_Px_S_Py_Px_C1000001001_c;
  Double I_ERI_Py_Py_Py_Px_C1000001001_c = I_ERI_D2y_S_Py_Px_C1000001001_c+ABY*I_ERI_Py_S_Py_Px_C1000001001_c;
  Double I_ERI_Pz_Py_Py_Px_C1000001001_c = I_ERI_Dyz_S_Py_Px_C1000001001_c+ABY*I_ERI_Pz_S_Py_Px_C1000001001_c;
  Double I_ERI_Px_Pz_Py_Px_C1000001001_c = I_ERI_Dxz_S_Py_Px_C1000001001_c+ABZ*I_ERI_Px_S_Py_Px_C1000001001_c;
  Double I_ERI_Py_Pz_Py_Px_C1000001001_c = I_ERI_Dyz_S_Py_Px_C1000001001_c+ABZ*I_ERI_Py_S_Py_Px_C1000001001_c;
  Double I_ERI_Pz_Pz_Py_Px_C1000001001_c = I_ERI_D2z_S_Py_Px_C1000001001_c+ABZ*I_ERI_Pz_S_Py_Px_C1000001001_c;
  Double I_ERI_Px_Px_Pz_Px_C1000001001_c = I_ERI_D2x_S_Pz_Px_C1000001001_c+ABX*I_ERI_Px_S_Pz_Px_C1000001001_c;
  Double I_ERI_Py_Px_Pz_Px_C1000001001_c = I_ERI_Dxy_S_Pz_Px_C1000001001_c+ABX*I_ERI_Py_S_Pz_Px_C1000001001_c;
  Double I_ERI_Pz_Px_Pz_Px_C1000001001_c = I_ERI_Dxz_S_Pz_Px_C1000001001_c+ABX*I_ERI_Pz_S_Pz_Px_C1000001001_c;
  Double I_ERI_Px_Py_Pz_Px_C1000001001_c = I_ERI_Dxy_S_Pz_Px_C1000001001_c+ABY*I_ERI_Px_S_Pz_Px_C1000001001_c;
  Double I_ERI_Py_Py_Pz_Px_C1000001001_c = I_ERI_D2y_S_Pz_Px_C1000001001_c+ABY*I_ERI_Py_S_Pz_Px_C1000001001_c;
  Double I_ERI_Pz_Py_Pz_Px_C1000001001_c = I_ERI_Dyz_S_Pz_Px_C1000001001_c+ABY*I_ERI_Pz_S_Pz_Px_C1000001001_c;
  Double I_ERI_Px_Pz_Pz_Px_C1000001001_c = I_ERI_Dxz_S_Pz_Px_C1000001001_c+ABZ*I_ERI_Px_S_Pz_Px_C1000001001_c;
  Double I_ERI_Py_Pz_Pz_Px_C1000001001_c = I_ERI_Dyz_S_Pz_Px_C1000001001_c+ABZ*I_ERI_Py_S_Pz_Px_C1000001001_c;
  Double I_ERI_Pz_Pz_Pz_Px_C1000001001_c = I_ERI_D2z_S_Pz_Px_C1000001001_c+ABZ*I_ERI_Pz_S_Pz_Px_C1000001001_c;
  Double I_ERI_Px_Px_Px_Py_C1000001001_c = I_ERI_D2x_S_Px_Py_C1000001001_c+ABX*I_ERI_Px_S_Px_Py_C1000001001_c;
  Double I_ERI_Py_Px_Px_Py_C1000001001_c = I_ERI_Dxy_S_Px_Py_C1000001001_c+ABX*I_ERI_Py_S_Px_Py_C1000001001_c;
  Double I_ERI_Pz_Px_Px_Py_C1000001001_c = I_ERI_Dxz_S_Px_Py_C1000001001_c+ABX*I_ERI_Pz_S_Px_Py_C1000001001_c;
  Double I_ERI_Px_Py_Px_Py_C1000001001_c = I_ERI_Dxy_S_Px_Py_C1000001001_c+ABY*I_ERI_Px_S_Px_Py_C1000001001_c;
  Double I_ERI_Py_Py_Px_Py_C1000001001_c = I_ERI_D2y_S_Px_Py_C1000001001_c+ABY*I_ERI_Py_S_Px_Py_C1000001001_c;
  Double I_ERI_Pz_Py_Px_Py_C1000001001_c = I_ERI_Dyz_S_Px_Py_C1000001001_c+ABY*I_ERI_Pz_S_Px_Py_C1000001001_c;
  Double I_ERI_Px_Pz_Px_Py_C1000001001_c = I_ERI_Dxz_S_Px_Py_C1000001001_c+ABZ*I_ERI_Px_S_Px_Py_C1000001001_c;
  Double I_ERI_Py_Pz_Px_Py_C1000001001_c = I_ERI_Dyz_S_Px_Py_C1000001001_c+ABZ*I_ERI_Py_S_Px_Py_C1000001001_c;
  Double I_ERI_Pz_Pz_Px_Py_C1000001001_c = I_ERI_D2z_S_Px_Py_C1000001001_c+ABZ*I_ERI_Pz_S_Px_Py_C1000001001_c;
  Double I_ERI_Px_Px_Py_Py_C1000001001_c = I_ERI_D2x_S_Py_Py_C1000001001_c+ABX*I_ERI_Px_S_Py_Py_C1000001001_c;
  Double I_ERI_Py_Px_Py_Py_C1000001001_c = I_ERI_Dxy_S_Py_Py_C1000001001_c+ABX*I_ERI_Py_S_Py_Py_C1000001001_c;
  Double I_ERI_Pz_Px_Py_Py_C1000001001_c = I_ERI_Dxz_S_Py_Py_C1000001001_c+ABX*I_ERI_Pz_S_Py_Py_C1000001001_c;
  Double I_ERI_Px_Py_Py_Py_C1000001001_c = I_ERI_Dxy_S_Py_Py_C1000001001_c+ABY*I_ERI_Px_S_Py_Py_C1000001001_c;
  Double I_ERI_Py_Py_Py_Py_C1000001001_c = I_ERI_D2y_S_Py_Py_C1000001001_c+ABY*I_ERI_Py_S_Py_Py_C1000001001_c;
  Double I_ERI_Pz_Py_Py_Py_C1000001001_c = I_ERI_Dyz_S_Py_Py_C1000001001_c+ABY*I_ERI_Pz_S_Py_Py_C1000001001_c;
  Double I_ERI_Px_Pz_Py_Py_C1000001001_c = I_ERI_Dxz_S_Py_Py_C1000001001_c+ABZ*I_ERI_Px_S_Py_Py_C1000001001_c;
  Double I_ERI_Py_Pz_Py_Py_C1000001001_c = I_ERI_Dyz_S_Py_Py_C1000001001_c+ABZ*I_ERI_Py_S_Py_Py_C1000001001_c;
  Double I_ERI_Pz_Pz_Py_Py_C1000001001_c = I_ERI_D2z_S_Py_Py_C1000001001_c+ABZ*I_ERI_Pz_S_Py_Py_C1000001001_c;
  Double I_ERI_Px_Px_Pz_Py_C1000001001_c = I_ERI_D2x_S_Pz_Py_C1000001001_c+ABX*I_ERI_Px_S_Pz_Py_C1000001001_c;
  Double I_ERI_Py_Px_Pz_Py_C1000001001_c = I_ERI_Dxy_S_Pz_Py_C1000001001_c+ABX*I_ERI_Py_S_Pz_Py_C1000001001_c;
  Double I_ERI_Pz_Px_Pz_Py_C1000001001_c = I_ERI_Dxz_S_Pz_Py_C1000001001_c+ABX*I_ERI_Pz_S_Pz_Py_C1000001001_c;
  Double I_ERI_Px_Py_Pz_Py_C1000001001_c = I_ERI_Dxy_S_Pz_Py_C1000001001_c+ABY*I_ERI_Px_S_Pz_Py_C1000001001_c;
  Double I_ERI_Py_Py_Pz_Py_C1000001001_c = I_ERI_D2y_S_Pz_Py_C1000001001_c+ABY*I_ERI_Py_S_Pz_Py_C1000001001_c;
  Double I_ERI_Pz_Py_Pz_Py_C1000001001_c = I_ERI_Dyz_S_Pz_Py_C1000001001_c+ABY*I_ERI_Pz_S_Pz_Py_C1000001001_c;
  Double I_ERI_Px_Pz_Pz_Py_C1000001001_c = I_ERI_Dxz_S_Pz_Py_C1000001001_c+ABZ*I_ERI_Px_S_Pz_Py_C1000001001_c;
  Double I_ERI_Py_Pz_Pz_Py_C1000001001_c = I_ERI_Dyz_S_Pz_Py_C1000001001_c+ABZ*I_ERI_Py_S_Pz_Py_C1000001001_c;
  Double I_ERI_Pz_Pz_Pz_Py_C1000001001_c = I_ERI_D2z_S_Pz_Py_C1000001001_c+ABZ*I_ERI_Pz_S_Pz_Py_C1000001001_c;
  Double I_ERI_Px_Px_Px_Pz_C1000001001_c = I_ERI_D2x_S_Px_Pz_C1000001001_c+ABX*I_ERI_Px_S_Px_Pz_C1000001001_c;
  Double I_ERI_Py_Px_Px_Pz_C1000001001_c = I_ERI_Dxy_S_Px_Pz_C1000001001_c+ABX*I_ERI_Py_S_Px_Pz_C1000001001_c;
  Double I_ERI_Pz_Px_Px_Pz_C1000001001_c = I_ERI_Dxz_S_Px_Pz_C1000001001_c+ABX*I_ERI_Pz_S_Px_Pz_C1000001001_c;
  Double I_ERI_Px_Py_Px_Pz_C1000001001_c = I_ERI_Dxy_S_Px_Pz_C1000001001_c+ABY*I_ERI_Px_S_Px_Pz_C1000001001_c;
  Double I_ERI_Py_Py_Px_Pz_C1000001001_c = I_ERI_D2y_S_Px_Pz_C1000001001_c+ABY*I_ERI_Py_S_Px_Pz_C1000001001_c;
  Double I_ERI_Pz_Py_Px_Pz_C1000001001_c = I_ERI_Dyz_S_Px_Pz_C1000001001_c+ABY*I_ERI_Pz_S_Px_Pz_C1000001001_c;
  Double I_ERI_Px_Pz_Px_Pz_C1000001001_c = I_ERI_Dxz_S_Px_Pz_C1000001001_c+ABZ*I_ERI_Px_S_Px_Pz_C1000001001_c;
  Double I_ERI_Py_Pz_Px_Pz_C1000001001_c = I_ERI_Dyz_S_Px_Pz_C1000001001_c+ABZ*I_ERI_Py_S_Px_Pz_C1000001001_c;
  Double I_ERI_Pz_Pz_Px_Pz_C1000001001_c = I_ERI_D2z_S_Px_Pz_C1000001001_c+ABZ*I_ERI_Pz_S_Px_Pz_C1000001001_c;
  Double I_ERI_Px_Px_Py_Pz_C1000001001_c = I_ERI_D2x_S_Py_Pz_C1000001001_c+ABX*I_ERI_Px_S_Py_Pz_C1000001001_c;
  Double I_ERI_Py_Px_Py_Pz_C1000001001_c = I_ERI_Dxy_S_Py_Pz_C1000001001_c+ABX*I_ERI_Py_S_Py_Pz_C1000001001_c;
  Double I_ERI_Pz_Px_Py_Pz_C1000001001_c = I_ERI_Dxz_S_Py_Pz_C1000001001_c+ABX*I_ERI_Pz_S_Py_Pz_C1000001001_c;
  Double I_ERI_Px_Py_Py_Pz_C1000001001_c = I_ERI_Dxy_S_Py_Pz_C1000001001_c+ABY*I_ERI_Px_S_Py_Pz_C1000001001_c;
  Double I_ERI_Py_Py_Py_Pz_C1000001001_c = I_ERI_D2y_S_Py_Pz_C1000001001_c+ABY*I_ERI_Py_S_Py_Pz_C1000001001_c;
  Double I_ERI_Pz_Py_Py_Pz_C1000001001_c = I_ERI_Dyz_S_Py_Pz_C1000001001_c+ABY*I_ERI_Pz_S_Py_Pz_C1000001001_c;
  Double I_ERI_Px_Pz_Py_Pz_C1000001001_c = I_ERI_Dxz_S_Py_Pz_C1000001001_c+ABZ*I_ERI_Px_S_Py_Pz_C1000001001_c;
  Double I_ERI_Py_Pz_Py_Pz_C1000001001_c = I_ERI_Dyz_S_Py_Pz_C1000001001_c+ABZ*I_ERI_Py_S_Py_Pz_C1000001001_c;
  Double I_ERI_Pz_Pz_Py_Pz_C1000001001_c = I_ERI_D2z_S_Py_Pz_C1000001001_c+ABZ*I_ERI_Pz_S_Py_Pz_C1000001001_c;
  Double I_ERI_Px_Px_Pz_Pz_C1000001001_c = I_ERI_D2x_S_Pz_Pz_C1000001001_c+ABX*I_ERI_Px_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Py_Px_Pz_Pz_C1000001001_c = I_ERI_Dxy_S_Pz_Pz_C1000001001_c+ABX*I_ERI_Py_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Pz_Px_Pz_Pz_C1000001001_c = I_ERI_Dxz_S_Pz_Pz_C1000001001_c+ABX*I_ERI_Pz_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Px_Py_Pz_Pz_C1000001001_c = I_ERI_Dxy_S_Pz_Pz_C1000001001_c+ABY*I_ERI_Px_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Py_Py_Pz_Pz_C1000001001_c = I_ERI_D2y_S_Pz_Pz_C1000001001_c+ABY*I_ERI_Py_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Pz_Py_Pz_Pz_C1000001001_c = I_ERI_Dyz_S_Pz_Pz_C1000001001_c+ABY*I_ERI_Pz_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Px_Pz_Pz_Pz_C1000001001_c = I_ERI_Dxz_S_Pz_Pz_C1000001001_c+ABZ*I_ERI_Px_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Py_Pz_Pz_Pz_C1000001001_c = I_ERI_Dyz_S_Pz_Pz_C1000001001_c+ABZ*I_ERI_Py_S_Pz_Pz_C1000001001_c;
  Double I_ERI_Pz_Pz_Pz_Pz_C1000001001_c = I_ERI_D2z_S_Pz_Pz_C1000001001_c+ABZ*I_ERI_Pz_S_Pz_Pz_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_P_C1001001001_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_P_C1001001001_c
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001001001_c
   ************************************************************/
  Double I_ERI_Px_Px_D2x_Px_C1001001001_c = I_ERI_D2x_S_D2x_Px_C1001001001_c+ABX*I_ERI_Px_S_D2x_Px_C1001001001_c;
  Double I_ERI_Py_Px_D2x_Px_C1001001001_c = I_ERI_Dxy_S_D2x_Px_C1001001001_c+ABX*I_ERI_Py_S_D2x_Px_C1001001001_c;
  Double I_ERI_Pz_Px_D2x_Px_C1001001001_c = I_ERI_Dxz_S_D2x_Px_C1001001001_c+ABX*I_ERI_Pz_S_D2x_Px_C1001001001_c;
  Double I_ERI_Px_Py_D2x_Px_C1001001001_c = I_ERI_Dxy_S_D2x_Px_C1001001001_c+ABY*I_ERI_Px_S_D2x_Px_C1001001001_c;
  Double I_ERI_Py_Py_D2x_Px_C1001001001_c = I_ERI_D2y_S_D2x_Px_C1001001001_c+ABY*I_ERI_Py_S_D2x_Px_C1001001001_c;
  Double I_ERI_Pz_Py_D2x_Px_C1001001001_c = I_ERI_Dyz_S_D2x_Px_C1001001001_c+ABY*I_ERI_Pz_S_D2x_Px_C1001001001_c;
  Double I_ERI_Px_Pz_D2x_Px_C1001001001_c = I_ERI_Dxz_S_D2x_Px_C1001001001_c+ABZ*I_ERI_Px_S_D2x_Px_C1001001001_c;
  Double I_ERI_Py_Pz_D2x_Px_C1001001001_c = I_ERI_Dyz_S_D2x_Px_C1001001001_c+ABZ*I_ERI_Py_S_D2x_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_D2x_Px_C1001001001_c = I_ERI_D2z_S_D2x_Px_C1001001001_c+ABZ*I_ERI_Pz_S_D2x_Px_C1001001001_c;
  Double I_ERI_Px_Px_Dxy_Px_C1001001001_c = I_ERI_D2x_S_Dxy_Px_C1001001001_c+ABX*I_ERI_Px_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Py_Px_Dxy_Px_C1001001001_c = I_ERI_Dxy_S_Dxy_Px_C1001001001_c+ABX*I_ERI_Py_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Pz_Px_Dxy_Px_C1001001001_c = I_ERI_Dxz_S_Dxy_Px_C1001001001_c+ABX*I_ERI_Pz_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Px_Py_Dxy_Px_C1001001001_c = I_ERI_Dxy_S_Dxy_Px_C1001001001_c+ABY*I_ERI_Px_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Py_Py_Dxy_Px_C1001001001_c = I_ERI_D2y_S_Dxy_Px_C1001001001_c+ABY*I_ERI_Py_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Pz_Py_Dxy_Px_C1001001001_c = I_ERI_Dyz_S_Dxy_Px_C1001001001_c+ABY*I_ERI_Pz_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Px_Pz_Dxy_Px_C1001001001_c = I_ERI_Dxz_S_Dxy_Px_C1001001001_c+ABZ*I_ERI_Px_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Py_Pz_Dxy_Px_C1001001001_c = I_ERI_Dyz_S_Dxy_Px_C1001001001_c+ABZ*I_ERI_Py_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxy_Px_C1001001001_c = I_ERI_D2z_S_Dxy_Px_C1001001001_c+ABZ*I_ERI_Pz_S_Dxy_Px_C1001001001_c;
  Double I_ERI_Px_Px_Dxz_Px_C1001001001_c = I_ERI_D2x_S_Dxz_Px_C1001001001_c+ABX*I_ERI_Px_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Py_Px_Dxz_Px_C1001001001_c = I_ERI_Dxy_S_Dxz_Px_C1001001001_c+ABX*I_ERI_Py_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Pz_Px_Dxz_Px_C1001001001_c = I_ERI_Dxz_S_Dxz_Px_C1001001001_c+ABX*I_ERI_Pz_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Px_Py_Dxz_Px_C1001001001_c = I_ERI_Dxy_S_Dxz_Px_C1001001001_c+ABY*I_ERI_Px_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Py_Py_Dxz_Px_C1001001001_c = I_ERI_D2y_S_Dxz_Px_C1001001001_c+ABY*I_ERI_Py_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Pz_Py_Dxz_Px_C1001001001_c = I_ERI_Dyz_S_Dxz_Px_C1001001001_c+ABY*I_ERI_Pz_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Px_Pz_Dxz_Px_C1001001001_c = I_ERI_Dxz_S_Dxz_Px_C1001001001_c+ABZ*I_ERI_Px_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Py_Pz_Dxz_Px_C1001001001_c = I_ERI_Dyz_S_Dxz_Px_C1001001001_c+ABZ*I_ERI_Py_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxz_Px_C1001001001_c = I_ERI_D2z_S_Dxz_Px_C1001001001_c+ABZ*I_ERI_Pz_S_Dxz_Px_C1001001001_c;
  Double I_ERI_Px_Px_D2y_Px_C1001001001_c = I_ERI_D2x_S_D2y_Px_C1001001001_c+ABX*I_ERI_Px_S_D2y_Px_C1001001001_c;
  Double I_ERI_Py_Px_D2y_Px_C1001001001_c = I_ERI_Dxy_S_D2y_Px_C1001001001_c+ABX*I_ERI_Py_S_D2y_Px_C1001001001_c;
  Double I_ERI_Pz_Px_D2y_Px_C1001001001_c = I_ERI_Dxz_S_D2y_Px_C1001001001_c+ABX*I_ERI_Pz_S_D2y_Px_C1001001001_c;
  Double I_ERI_Px_Py_D2y_Px_C1001001001_c = I_ERI_Dxy_S_D2y_Px_C1001001001_c+ABY*I_ERI_Px_S_D2y_Px_C1001001001_c;
  Double I_ERI_Py_Py_D2y_Px_C1001001001_c = I_ERI_D2y_S_D2y_Px_C1001001001_c+ABY*I_ERI_Py_S_D2y_Px_C1001001001_c;
  Double I_ERI_Pz_Py_D2y_Px_C1001001001_c = I_ERI_Dyz_S_D2y_Px_C1001001001_c+ABY*I_ERI_Pz_S_D2y_Px_C1001001001_c;
  Double I_ERI_Px_Pz_D2y_Px_C1001001001_c = I_ERI_Dxz_S_D2y_Px_C1001001001_c+ABZ*I_ERI_Px_S_D2y_Px_C1001001001_c;
  Double I_ERI_Py_Pz_D2y_Px_C1001001001_c = I_ERI_Dyz_S_D2y_Px_C1001001001_c+ABZ*I_ERI_Py_S_D2y_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_D2y_Px_C1001001001_c = I_ERI_D2z_S_D2y_Px_C1001001001_c+ABZ*I_ERI_Pz_S_D2y_Px_C1001001001_c;
  Double I_ERI_Px_Px_Dyz_Px_C1001001001_c = I_ERI_D2x_S_Dyz_Px_C1001001001_c+ABX*I_ERI_Px_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Py_Px_Dyz_Px_C1001001001_c = I_ERI_Dxy_S_Dyz_Px_C1001001001_c+ABX*I_ERI_Py_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Pz_Px_Dyz_Px_C1001001001_c = I_ERI_Dxz_S_Dyz_Px_C1001001001_c+ABX*I_ERI_Pz_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Px_Py_Dyz_Px_C1001001001_c = I_ERI_Dxy_S_Dyz_Px_C1001001001_c+ABY*I_ERI_Px_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Py_Py_Dyz_Px_C1001001001_c = I_ERI_D2y_S_Dyz_Px_C1001001001_c+ABY*I_ERI_Py_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Pz_Py_Dyz_Px_C1001001001_c = I_ERI_Dyz_S_Dyz_Px_C1001001001_c+ABY*I_ERI_Pz_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Px_Pz_Dyz_Px_C1001001001_c = I_ERI_Dxz_S_Dyz_Px_C1001001001_c+ABZ*I_ERI_Px_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Py_Pz_Dyz_Px_C1001001001_c = I_ERI_Dyz_S_Dyz_Px_C1001001001_c+ABZ*I_ERI_Py_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_Dyz_Px_C1001001001_c = I_ERI_D2z_S_Dyz_Px_C1001001001_c+ABZ*I_ERI_Pz_S_Dyz_Px_C1001001001_c;
  Double I_ERI_Px_Px_D2z_Px_C1001001001_c = I_ERI_D2x_S_D2z_Px_C1001001001_c+ABX*I_ERI_Px_S_D2z_Px_C1001001001_c;
  Double I_ERI_Py_Px_D2z_Px_C1001001001_c = I_ERI_Dxy_S_D2z_Px_C1001001001_c+ABX*I_ERI_Py_S_D2z_Px_C1001001001_c;
  Double I_ERI_Pz_Px_D2z_Px_C1001001001_c = I_ERI_Dxz_S_D2z_Px_C1001001001_c+ABX*I_ERI_Pz_S_D2z_Px_C1001001001_c;
  Double I_ERI_Px_Py_D2z_Px_C1001001001_c = I_ERI_Dxy_S_D2z_Px_C1001001001_c+ABY*I_ERI_Px_S_D2z_Px_C1001001001_c;
  Double I_ERI_Py_Py_D2z_Px_C1001001001_c = I_ERI_D2y_S_D2z_Px_C1001001001_c+ABY*I_ERI_Py_S_D2z_Px_C1001001001_c;
  Double I_ERI_Pz_Py_D2z_Px_C1001001001_c = I_ERI_Dyz_S_D2z_Px_C1001001001_c+ABY*I_ERI_Pz_S_D2z_Px_C1001001001_c;
  Double I_ERI_Px_Pz_D2z_Px_C1001001001_c = I_ERI_Dxz_S_D2z_Px_C1001001001_c+ABZ*I_ERI_Px_S_D2z_Px_C1001001001_c;
  Double I_ERI_Py_Pz_D2z_Px_C1001001001_c = I_ERI_Dyz_S_D2z_Px_C1001001001_c+ABZ*I_ERI_Py_S_D2z_Px_C1001001001_c;
  Double I_ERI_Pz_Pz_D2z_Px_C1001001001_c = I_ERI_D2z_S_D2z_Px_C1001001001_c+ABZ*I_ERI_Pz_S_D2z_Px_C1001001001_c;
  Double I_ERI_Px_Px_D2x_Py_C1001001001_c = I_ERI_D2x_S_D2x_Py_C1001001001_c+ABX*I_ERI_Px_S_D2x_Py_C1001001001_c;
  Double I_ERI_Py_Px_D2x_Py_C1001001001_c = I_ERI_Dxy_S_D2x_Py_C1001001001_c+ABX*I_ERI_Py_S_D2x_Py_C1001001001_c;
  Double I_ERI_Pz_Px_D2x_Py_C1001001001_c = I_ERI_Dxz_S_D2x_Py_C1001001001_c+ABX*I_ERI_Pz_S_D2x_Py_C1001001001_c;
  Double I_ERI_Px_Py_D2x_Py_C1001001001_c = I_ERI_Dxy_S_D2x_Py_C1001001001_c+ABY*I_ERI_Px_S_D2x_Py_C1001001001_c;
  Double I_ERI_Py_Py_D2x_Py_C1001001001_c = I_ERI_D2y_S_D2x_Py_C1001001001_c+ABY*I_ERI_Py_S_D2x_Py_C1001001001_c;
  Double I_ERI_Pz_Py_D2x_Py_C1001001001_c = I_ERI_Dyz_S_D2x_Py_C1001001001_c+ABY*I_ERI_Pz_S_D2x_Py_C1001001001_c;
  Double I_ERI_Px_Pz_D2x_Py_C1001001001_c = I_ERI_Dxz_S_D2x_Py_C1001001001_c+ABZ*I_ERI_Px_S_D2x_Py_C1001001001_c;
  Double I_ERI_Py_Pz_D2x_Py_C1001001001_c = I_ERI_Dyz_S_D2x_Py_C1001001001_c+ABZ*I_ERI_Py_S_D2x_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_D2x_Py_C1001001001_c = I_ERI_D2z_S_D2x_Py_C1001001001_c+ABZ*I_ERI_Pz_S_D2x_Py_C1001001001_c;
  Double I_ERI_Px_Px_Dxy_Py_C1001001001_c = I_ERI_D2x_S_Dxy_Py_C1001001001_c+ABX*I_ERI_Px_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Py_Px_Dxy_Py_C1001001001_c = I_ERI_Dxy_S_Dxy_Py_C1001001001_c+ABX*I_ERI_Py_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Pz_Px_Dxy_Py_C1001001001_c = I_ERI_Dxz_S_Dxy_Py_C1001001001_c+ABX*I_ERI_Pz_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Px_Py_Dxy_Py_C1001001001_c = I_ERI_Dxy_S_Dxy_Py_C1001001001_c+ABY*I_ERI_Px_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Py_Py_Dxy_Py_C1001001001_c = I_ERI_D2y_S_Dxy_Py_C1001001001_c+ABY*I_ERI_Py_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Pz_Py_Dxy_Py_C1001001001_c = I_ERI_Dyz_S_Dxy_Py_C1001001001_c+ABY*I_ERI_Pz_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Px_Pz_Dxy_Py_C1001001001_c = I_ERI_Dxz_S_Dxy_Py_C1001001001_c+ABZ*I_ERI_Px_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Py_Pz_Dxy_Py_C1001001001_c = I_ERI_Dyz_S_Dxy_Py_C1001001001_c+ABZ*I_ERI_Py_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxy_Py_C1001001001_c = I_ERI_D2z_S_Dxy_Py_C1001001001_c+ABZ*I_ERI_Pz_S_Dxy_Py_C1001001001_c;
  Double I_ERI_Px_Px_Dxz_Py_C1001001001_c = I_ERI_D2x_S_Dxz_Py_C1001001001_c+ABX*I_ERI_Px_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Py_Px_Dxz_Py_C1001001001_c = I_ERI_Dxy_S_Dxz_Py_C1001001001_c+ABX*I_ERI_Py_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Pz_Px_Dxz_Py_C1001001001_c = I_ERI_Dxz_S_Dxz_Py_C1001001001_c+ABX*I_ERI_Pz_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Px_Py_Dxz_Py_C1001001001_c = I_ERI_Dxy_S_Dxz_Py_C1001001001_c+ABY*I_ERI_Px_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Py_Py_Dxz_Py_C1001001001_c = I_ERI_D2y_S_Dxz_Py_C1001001001_c+ABY*I_ERI_Py_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Pz_Py_Dxz_Py_C1001001001_c = I_ERI_Dyz_S_Dxz_Py_C1001001001_c+ABY*I_ERI_Pz_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Px_Pz_Dxz_Py_C1001001001_c = I_ERI_Dxz_S_Dxz_Py_C1001001001_c+ABZ*I_ERI_Px_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Py_Pz_Dxz_Py_C1001001001_c = I_ERI_Dyz_S_Dxz_Py_C1001001001_c+ABZ*I_ERI_Py_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxz_Py_C1001001001_c = I_ERI_D2z_S_Dxz_Py_C1001001001_c+ABZ*I_ERI_Pz_S_Dxz_Py_C1001001001_c;
  Double I_ERI_Px_Px_D2y_Py_C1001001001_c = I_ERI_D2x_S_D2y_Py_C1001001001_c+ABX*I_ERI_Px_S_D2y_Py_C1001001001_c;
  Double I_ERI_Py_Px_D2y_Py_C1001001001_c = I_ERI_Dxy_S_D2y_Py_C1001001001_c+ABX*I_ERI_Py_S_D2y_Py_C1001001001_c;
  Double I_ERI_Pz_Px_D2y_Py_C1001001001_c = I_ERI_Dxz_S_D2y_Py_C1001001001_c+ABX*I_ERI_Pz_S_D2y_Py_C1001001001_c;
  Double I_ERI_Px_Py_D2y_Py_C1001001001_c = I_ERI_Dxy_S_D2y_Py_C1001001001_c+ABY*I_ERI_Px_S_D2y_Py_C1001001001_c;
  Double I_ERI_Py_Py_D2y_Py_C1001001001_c = I_ERI_D2y_S_D2y_Py_C1001001001_c+ABY*I_ERI_Py_S_D2y_Py_C1001001001_c;
  Double I_ERI_Pz_Py_D2y_Py_C1001001001_c = I_ERI_Dyz_S_D2y_Py_C1001001001_c+ABY*I_ERI_Pz_S_D2y_Py_C1001001001_c;
  Double I_ERI_Px_Pz_D2y_Py_C1001001001_c = I_ERI_Dxz_S_D2y_Py_C1001001001_c+ABZ*I_ERI_Px_S_D2y_Py_C1001001001_c;
  Double I_ERI_Py_Pz_D2y_Py_C1001001001_c = I_ERI_Dyz_S_D2y_Py_C1001001001_c+ABZ*I_ERI_Py_S_D2y_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_D2y_Py_C1001001001_c = I_ERI_D2z_S_D2y_Py_C1001001001_c+ABZ*I_ERI_Pz_S_D2y_Py_C1001001001_c;
  Double I_ERI_Px_Px_Dyz_Py_C1001001001_c = I_ERI_D2x_S_Dyz_Py_C1001001001_c+ABX*I_ERI_Px_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Py_Px_Dyz_Py_C1001001001_c = I_ERI_Dxy_S_Dyz_Py_C1001001001_c+ABX*I_ERI_Py_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Pz_Px_Dyz_Py_C1001001001_c = I_ERI_Dxz_S_Dyz_Py_C1001001001_c+ABX*I_ERI_Pz_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Px_Py_Dyz_Py_C1001001001_c = I_ERI_Dxy_S_Dyz_Py_C1001001001_c+ABY*I_ERI_Px_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Py_Py_Dyz_Py_C1001001001_c = I_ERI_D2y_S_Dyz_Py_C1001001001_c+ABY*I_ERI_Py_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Pz_Py_Dyz_Py_C1001001001_c = I_ERI_Dyz_S_Dyz_Py_C1001001001_c+ABY*I_ERI_Pz_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Px_Pz_Dyz_Py_C1001001001_c = I_ERI_Dxz_S_Dyz_Py_C1001001001_c+ABZ*I_ERI_Px_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Py_Pz_Dyz_Py_C1001001001_c = I_ERI_Dyz_S_Dyz_Py_C1001001001_c+ABZ*I_ERI_Py_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_Dyz_Py_C1001001001_c = I_ERI_D2z_S_Dyz_Py_C1001001001_c+ABZ*I_ERI_Pz_S_Dyz_Py_C1001001001_c;
  Double I_ERI_Px_Px_D2z_Py_C1001001001_c = I_ERI_D2x_S_D2z_Py_C1001001001_c+ABX*I_ERI_Px_S_D2z_Py_C1001001001_c;
  Double I_ERI_Py_Px_D2z_Py_C1001001001_c = I_ERI_Dxy_S_D2z_Py_C1001001001_c+ABX*I_ERI_Py_S_D2z_Py_C1001001001_c;
  Double I_ERI_Pz_Px_D2z_Py_C1001001001_c = I_ERI_Dxz_S_D2z_Py_C1001001001_c+ABX*I_ERI_Pz_S_D2z_Py_C1001001001_c;
  Double I_ERI_Px_Py_D2z_Py_C1001001001_c = I_ERI_Dxy_S_D2z_Py_C1001001001_c+ABY*I_ERI_Px_S_D2z_Py_C1001001001_c;
  Double I_ERI_Py_Py_D2z_Py_C1001001001_c = I_ERI_D2y_S_D2z_Py_C1001001001_c+ABY*I_ERI_Py_S_D2z_Py_C1001001001_c;
  Double I_ERI_Pz_Py_D2z_Py_C1001001001_c = I_ERI_Dyz_S_D2z_Py_C1001001001_c+ABY*I_ERI_Pz_S_D2z_Py_C1001001001_c;
  Double I_ERI_Px_Pz_D2z_Py_C1001001001_c = I_ERI_Dxz_S_D2z_Py_C1001001001_c+ABZ*I_ERI_Px_S_D2z_Py_C1001001001_c;
  Double I_ERI_Py_Pz_D2z_Py_C1001001001_c = I_ERI_Dyz_S_D2z_Py_C1001001001_c+ABZ*I_ERI_Py_S_D2z_Py_C1001001001_c;
  Double I_ERI_Pz_Pz_D2z_Py_C1001001001_c = I_ERI_D2z_S_D2z_Py_C1001001001_c+ABZ*I_ERI_Pz_S_D2z_Py_C1001001001_c;
  Double I_ERI_Px_Px_D2x_Pz_C1001001001_c = I_ERI_D2x_S_D2x_Pz_C1001001001_c+ABX*I_ERI_Px_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Py_Px_D2x_Pz_C1001001001_c = I_ERI_Dxy_S_D2x_Pz_C1001001001_c+ABX*I_ERI_Py_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_D2x_Pz_C1001001001_c = I_ERI_Dxz_S_D2x_Pz_C1001001001_c+ABX*I_ERI_Pz_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Px_Py_D2x_Pz_C1001001001_c = I_ERI_Dxy_S_D2x_Pz_C1001001001_c+ABY*I_ERI_Px_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Py_Py_D2x_Pz_C1001001001_c = I_ERI_D2y_S_D2x_Pz_C1001001001_c+ABY*I_ERI_Py_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_D2x_Pz_C1001001001_c = I_ERI_Dyz_S_D2x_Pz_C1001001001_c+ABY*I_ERI_Pz_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_D2x_Pz_C1001001001_c = I_ERI_Dxz_S_D2x_Pz_C1001001001_c+ABZ*I_ERI_Px_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_D2x_Pz_C1001001001_c = I_ERI_Dyz_S_D2x_Pz_C1001001001_c+ABZ*I_ERI_Py_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_D2x_Pz_C1001001001_c = I_ERI_D2z_S_D2x_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_D2x_Pz_C1001001001_c;
  Double I_ERI_Px_Px_Dxy_Pz_C1001001001_c = I_ERI_D2x_S_Dxy_Pz_C1001001001_c+ABX*I_ERI_Px_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Py_Px_Dxy_Pz_C1001001001_c = I_ERI_Dxy_S_Dxy_Pz_C1001001001_c+ABX*I_ERI_Py_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_Dxy_Pz_C1001001001_c = I_ERI_Dxz_S_Dxy_Pz_C1001001001_c+ABX*I_ERI_Pz_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Px_Py_Dxy_Pz_C1001001001_c = I_ERI_Dxy_S_Dxy_Pz_C1001001001_c+ABY*I_ERI_Px_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Py_Py_Dxy_Pz_C1001001001_c = I_ERI_D2y_S_Dxy_Pz_C1001001001_c+ABY*I_ERI_Py_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_Dxy_Pz_C1001001001_c = I_ERI_Dyz_S_Dxy_Pz_C1001001001_c+ABY*I_ERI_Pz_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_Dxy_Pz_C1001001001_c = I_ERI_Dxz_S_Dxy_Pz_C1001001001_c+ABZ*I_ERI_Px_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_Dxy_Pz_C1001001001_c = I_ERI_Dyz_S_Dxy_Pz_C1001001001_c+ABZ*I_ERI_Py_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxy_Pz_C1001001001_c = I_ERI_D2z_S_Dxy_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_Dxy_Pz_C1001001001_c;
  Double I_ERI_Px_Px_Dxz_Pz_C1001001001_c = I_ERI_D2x_S_Dxz_Pz_C1001001001_c+ABX*I_ERI_Px_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Py_Px_Dxz_Pz_C1001001001_c = I_ERI_Dxy_S_Dxz_Pz_C1001001001_c+ABX*I_ERI_Py_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_Dxz_Pz_C1001001001_c = I_ERI_Dxz_S_Dxz_Pz_C1001001001_c+ABX*I_ERI_Pz_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Px_Py_Dxz_Pz_C1001001001_c = I_ERI_Dxy_S_Dxz_Pz_C1001001001_c+ABY*I_ERI_Px_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Py_Py_Dxz_Pz_C1001001001_c = I_ERI_D2y_S_Dxz_Pz_C1001001001_c+ABY*I_ERI_Py_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_Dxz_Pz_C1001001001_c = I_ERI_Dyz_S_Dxz_Pz_C1001001001_c+ABY*I_ERI_Pz_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_Dxz_Pz_C1001001001_c = I_ERI_Dxz_S_Dxz_Pz_C1001001001_c+ABZ*I_ERI_Px_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_Dxz_Pz_C1001001001_c = I_ERI_Dyz_S_Dxz_Pz_C1001001001_c+ABZ*I_ERI_Py_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_Dxz_Pz_C1001001001_c = I_ERI_D2z_S_Dxz_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_Dxz_Pz_C1001001001_c;
  Double I_ERI_Px_Px_D2y_Pz_C1001001001_c = I_ERI_D2x_S_D2y_Pz_C1001001001_c+ABX*I_ERI_Px_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Py_Px_D2y_Pz_C1001001001_c = I_ERI_Dxy_S_D2y_Pz_C1001001001_c+ABX*I_ERI_Py_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_D2y_Pz_C1001001001_c = I_ERI_Dxz_S_D2y_Pz_C1001001001_c+ABX*I_ERI_Pz_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Px_Py_D2y_Pz_C1001001001_c = I_ERI_Dxy_S_D2y_Pz_C1001001001_c+ABY*I_ERI_Px_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Py_Py_D2y_Pz_C1001001001_c = I_ERI_D2y_S_D2y_Pz_C1001001001_c+ABY*I_ERI_Py_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_D2y_Pz_C1001001001_c = I_ERI_Dyz_S_D2y_Pz_C1001001001_c+ABY*I_ERI_Pz_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_D2y_Pz_C1001001001_c = I_ERI_Dxz_S_D2y_Pz_C1001001001_c+ABZ*I_ERI_Px_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_D2y_Pz_C1001001001_c = I_ERI_Dyz_S_D2y_Pz_C1001001001_c+ABZ*I_ERI_Py_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_D2y_Pz_C1001001001_c = I_ERI_D2z_S_D2y_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_D2y_Pz_C1001001001_c;
  Double I_ERI_Px_Px_Dyz_Pz_C1001001001_c = I_ERI_D2x_S_Dyz_Pz_C1001001001_c+ABX*I_ERI_Px_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Py_Px_Dyz_Pz_C1001001001_c = I_ERI_Dxy_S_Dyz_Pz_C1001001001_c+ABX*I_ERI_Py_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_Dyz_Pz_C1001001001_c = I_ERI_Dxz_S_Dyz_Pz_C1001001001_c+ABX*I_ERI_Pz_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Px_Py_Dyz_Pz_C1001001001_c = I_ERI_Dxy_S_Dyz_Pz_C1001001001_c+ABY*I_ERI_Px_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Py_Py_Dyz_Pz_C1001001001_c = I_ERI_D2y_S_Dyz_Pz_C1001001001_c+ABY*I_ERI_Py_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_Dyz_Pz_C1001001001_c = I_ERI_Dyz_S_Dyz_Pz_C1001001001_c+ABY*I_ERI_Pz_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_Dyz_Pz_C1001001001_c = I_ERI_Dxz_S_Dyz_Pz_C1001001001_c+ABZ*I_ERI_Px_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_Dyz_Pz_C1001001001_c = I_ERI_Dyz_S_Dyz_Pz_C1001001001_c+ABZ*I_ERI_Py_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_Dyz_Pz_C1001001001_c = I_ERI_D2z_S_Dyz_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_Dyz_Pz_C1001001001_c;
  Double I_ERI_Px_Px_D2z_Pz_C1001001001_c = I_ERI_D2x_S_D2z_Pz_C1001001001_c+ABX*I_ERI_Px_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Py_Px_D2z_Pz_C1001001001_c = I_ERI_Dxy_S_D2z_Pz_C1001001001_c+ABX*I_ERI_Py_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Pz_Px_D2z_Pz_C1001001001_c = I_ERI_Dxz_S_D2z_Pz_C1001001001_c+ABX*I_ERI_Pz_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Px_Py_D2z_Pz_C1001001001_c = I_ERI_Dxy_S_D2z_Pz_C1001001001_c+ABY*I_ERI_Px_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Py_Py_D2z_Pz_C1001001001_c = I_ERI_D2y_S_D2z_Pz_C1001001001_c+ABY*I_ERI_Py_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Pz_Py_D2z_Pz_C1001001001_c = I_ERI_Dyz_S_D2z_Pz_C1001001001_c+ABY*I_ERI_Pz_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Px_Pz_D2z_Pz_C1001001001_c = I_ERI_Dxz_S_D2z_Pz_C1001001001_c+ABZ*I_ERI_Px_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Py_Pz_D2z_Pz_C1001001001_c = I_ERI_Dyz_S_D2z_Pz_C1001001001_c+ABZ*I_ERI_Py_S_D2z_Pz_C1001001001_c;
  Double I_ERI_Pz_Pz_D2z_Pz_C1001001001_c = I_ERI_D2z_S_D2z_Pz_C1001001001_c+ABZ*I_ERI_Pz_S_D2z_Pz_C1001001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C0_a
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_Px_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1
   ************************************************************/
  abcd[1] = 2.0E0*I_ERI_D2x_S_S_S_C1_a-1*I_ERI_S_S_S_S_C1;
  abcd[2] = 2.0E0*I_ERI_Dxy_S_S_S_C1_a;
  abcd[3] = 2.0E0*I_ERI_Dxz_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1000_a
   ************************************************************/
  abcd[4] = 2.0E0*I_ERI_Px_Px_S_S_C1000_a;
  abcd[8] = 2.0E0*I_ERI_Px_Py_S_S_C1000_a;
  abcd[12] = 2.0E0*I_ERI_Px_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001
   ************************************************************/
  abcd[5] = 2.0E0*I_ERI_D2x_Px_S_S_C1001_a-1*I_ERI_S_Px_S_S_C1001;
  abcd[6] = 2.0E0*I_ERI_Dxy_Px_S_S_C1001_a;
  abcd[7] = 2.0E0*I_ERI_Dxz_Px_S_S_C1001_a;
  abcd[9] = 2.0E0*I_ERI_D2x_Py_S_S_C1001_a-1*I_ERI_S_Py_S_S_C1001;
  abcd[10] = 2.0E0*I_ERI_Dxy_Py_S_S_C1001_a;
  abcd[11] = 2.0E0*I_ERI_Dxz_Py_S_S_C1001_a;
  abcd[13] = 2.0E0*I_ERI_D2x_Pz_S_S_C1001_a-1*I_ERI_S_Pz_S_S_C1001;
  abcd[14] = 2.0E0*I_ERI_Dxy_Pz_S_S_C1001_a;
  abcd[15] = 2.0E0*I_ERI_Dxz_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[16] = 2.0E0*I_ERI_Px_S_Px_S_C1000000_a;
  abcd[32] = 2.0E0*I_ERI_Px_S_Py_S_C1000000_a;
  abcd[48] = 2.0E0*I_ERI_Px_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[17] = 2.0E0*I_ERI_D2x_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[18] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[19] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[33] = 2.0E0*I_ERI_D2x_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[34] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[35] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[49] = 2.0E0*I_ERI_D2x_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[50] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[51] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[20] = 2.0E0*I_ERI_Px_Px_Px_S_C1001000_a;
  abcd[24] = 2.0E0*I_ERI_Px_Py_Px_S_C1001000_a;
  abcd[28] = 2.0E0*I_ERI_Px_Pz_Px_S_C1001000_a;
  abcd[36] = 2.0E0*I_ERI_Px_Px_Py_S_C1001000_a;
  abcd[40] = 2.0E0*I_ERI_Px_Py_Py_S_C1001000_a;
  abcd[44] = 2.0E0*I_ERI_Px_Pz_Py_S_C1001000_a;
  abcd[52] = 2.0E0*I_ERI_Px_Px_Pz_S_C1001000_a;
  abcd[56] = 2.0E0*I_ERI_Px_Py_Pz_S_C1001000_a;
  abcd[60] = 2.0E0*I_ERI_Px_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[21] = 2.0E0*I_ERI_D2x_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[22] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[23] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[25] = 2.0E0*I_ERI_D2x_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[26] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[27] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[29] = 2.0E0*I_ERI_D2x_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[30] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[31] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[37] = 2.0E0*I_ERI_D2x_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[38] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[39] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[41] = 2.0E0*I_ERI_D2x_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[42] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[43] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[45] = 2.0E0*I_ERI_D2x_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[46] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[47] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[53] = 2.0E0*I_ERI_D2x_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[54] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[55] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[57] = 2.0E0*I_ERI_D2x_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[58] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[59] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[61] = 2.0E0*I_ERI_D2x_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[62] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[63] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000000000_a
   ************************************************************/
  abcd[64] = 2.0E0*I_ERI_Px_S_S_Px_C1000000000_a;
  abcd[128] = 2.0E0*I_ERI_Px_S_S_Py_C1000000000_a;
  abcd[192] = 2.0E0*I_ERI_Px_S_S_Pz_C1000000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000001_a
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000000001
   ************************************************************/
  abcd[65] = 2.0E0*I_ERI_D2x_S_S_Px_C1000000001_a-1*I_ERI_S_S_S_Px_C1000000001;
  abcd[66] = 2.0E0*I_ERI_Dxy_S_S_Px_C1000000001_a;
  abcd[67] = 2.0E0*I_ERI_Dxz_S_S_Px_C1000000001_a;
  abcd[129] = 2.0E0*I_ERI_D2x_S_S_Py_C1000000001_a-1*I_ERI_S_S_S_Py_C1000000001;
  abcd[130] = 2.0E0*I_ERI_Dxy_S_S_Py_C1000000001_a;
  abcd[131] = 2.0E0*I_ERI_Dxz_S_S_Py_C1000000001_a;
  abcd[193] = 2.0E0*I_ERI_D2x_S_S_Pz_C1000000001_a-1*I_ERI_S_S_S_Pz_C1000000001;
  abcd[194] = 2.0E0*I_ERI_Dxy_S_S_Pz_C1000000001_a;
  abcd[195] = 2.0E0*I_ERI_Dxz_S_S_Pz_C1000000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000001000_a
   ************************************************************/
  abcd[68] = 2.0E0*I_ERI_Px_Px_S_Px_C1000001000_a;
  abcd[72] = 2.0E0*I_ERI_Px_Py_S_Px_C1000001000_a;
  abcd[76] = 2.0E0*I_ERI_Px_Pz_S_Px_C1000001000_a;
  abcd[132] = 2.0E0*I_ERI_Px_Px_S_Py_C1000001000_a;
  abcd[136] = 2.0E0*I_ERI_Px_Py_S_Py_C1000001000_a;
  abcd[140] = 2.0E0*I_ERI_Px_Pz_S_Py_C1000001000_a;
  abcd[196] = 2.0E0*I_ERI_Px_Px_S_Pz_C1000001000_a;
  abcd[200] = 2.0E0*I_ERI_Px_Py_S_Pz_C1000001000_a;
  abcd[204] = 2.0E0*I_ERI_Px_Pz_S_Pz_C1000001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_P_C1000001001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000001001
   ************************************************************/
  abcd[69] = 2.0E0*I_ERI_D2x_Px_S_Px_C1000001001_a-1*I_ERI_S_Px_S_Px_C1000001001;
  abcd[70] = 2.0E0*I_ERI_Dxy_Px_S_Px_C1000001001_a;
  abcd[71] = 2.0E0*I_ERI_Dxz_Px_S_Px_C1000001001_a;
  abcd[73] = 2.0E0*I_ERI_D2x_Py_S_Px_C1000001001_a-1*I_ERI_S_Py_S_Px_C1000001001;
  abcd[74] = 2.0E0*I_ERI_Dxy_Py_S_Px_C1000001001_a;
  abcd[75] = 2.0E0*I_ERI_Dxz_Py_S_Px_C1000001001_a;
  abcd[77] = 2.0E0*I_ERI_D2x_Pz_S_Px_C1000001001_a-1*I_ERI_S_Pz_S_Px_C1000001001;
  abcd[78] = 2.0E0*I_ERI_Dxy_Pz_S_Px_C1000001001_a;
  abcd[79] = 2.0E0*I_ERI_Dxz_Pz_S_Px_C1000001001_a;
  abcd[133] = 2.0E0*I_ERI_D2x_Px_S_Py_C1000001001_a-1*I_ERI_S_Px_S_Py_C1000001001;
  abcd[134] = 2.0E0*I_ERI_Dxy_Px_S_Py_C1000001001_a;
  abcd[135] = 2.0E0*I_ERI_Dxz_Px_S_Py_C1000001001_a;
  abcd[137] = 2.0E0*I_ERI_D2x_Py_S_Py_C1000001001_a-1*I_ERI_S_Py_S_Py_C1000001001;
  abcd[138] = 2.0E0*I_ERI_Dxy_Py_S_Py_C1000001001_a;
  abcd[139] = 2.0E0*I_ERI_Dxz_Py_S_Py_C1000001001_a;
  abcd[141] = 2.0E0*I_ERI_D2x_Pz_S_Py_C1000001001_a-1*I_ERI_S_Pz_S_Py_C1000001001;
  abcd[142] = 2.0E0*I_ERI_Dxy_Pz_S_Py_C1000001001_a;
  abcd[143] = 2.0E0*I_ERI_Dxz_Pz_S_Py_C1000001001_a;
  abcd[197] = 2.0E0*I_ERI_D2x_Px_S_Pz_C1000001001_a-1*I_ERI_S_Px_S_Pz_C1000001001;
  abcd[198] = 2.0E0*I_ERI_Dxy_Px_S_Pz_C1000001001_a;
  abcd[199] = 2.0E0*I_ERI_Dxz_Px_S_Pz_C1000001001_a;
  abcd[201] = 2.0E0*I_ERI_D2x_Py_S_Pz_C1000001001_a-1*I_ERI_S_Py_S_Pz_C1000001001;
  abcd[202] = 2.0E0*I_ERI_Dxy_Py_S_Pz_C1000001001_a;
  abcd[203] = 2.0E0*I_ERI_Dxz_Py_S_Pz_C1000001001_a;
  abcd[205] = 2.0E0*I_ERI_D2x_Pz_S_Pz_C1000001001_a-1*I_ERI_S_Pz_S_Pz_C1000001001;
  abcd[206] = 2.0E0*I_ERI_Dxy_Pz_S_Pz_C1000001001_a;
  abcd[207] = 2.0E0*I_ERI_Dxz_Pz_S_Pz_C1000001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001000000_a
   ************************************************************/
  abcd[80] = 2.0E0*I_ERI_Px_S_Px_Px_C1001000000_a;
  abcd[96] = 2.0E0*I_ERI_Px_S_Py_Px_C1001000000_a;
  abcd[112] = 2.0E0*I_ERI_Px_S_Pz_Px_C1001000000_a;
  abcd[144] = 2.0E0*I_ERI_Px_S_Px_Py_C1001000000_a;
  abcd[160] = 2.0E0*I_ERI_Px_S_Py_Py_C1001000000_a;
  abcd[176] = 2.0E0*I_ERI_Px_S_Pz_Py_C1001000000_a;
  abcd[208] = 2.0E0*I_ERI_Px_S_Px_Pz_C1001000000_a;
  abcd[224] = 2.0E0*I_ERI_Px_S_Py_Pz_C1001000000_a;
  abcd[240] = 2.0E0*I_ERI_Px_S_Pz_Pz_C1001000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001000001
   ************************************************************/
  abcd[81] = 2.0E0*I_ERI_D2x_S_Px_Px_C1001000001_a-1*I_ERI_S_S_Px_Px_C1001000001;
  abcd[82] = 2.0E0*I_ERI_Dxy_S_Px_Px_C1001000001_a;
  abcd[83] = 2.0E0*I_ERI_Dxz_S_Px_Px_C1001000001_a;
  abcd[97] = 2.0E0*I_ERI_D2x_S_Py_Px_C1001000001_a-1*I_ERI_S_S_Py_Px_C1001000001;
  abcd[98] = 2.0E0*I_ERI_Dxy_S_Py_Px_C1001000001_a;
  abcd[99] = 2.0E0*I_ERI_Dxz_S_Py_Px_C1001000001_a;
  abcd[113] = 2.0E0*I_ERI_D2x_S_Pz_Px_C1001000001_a-1*I_ERI_S_S_Pz_Px_C1001000001;
  abcd[114] = 2.0E0*I_ERI_Dxy_S_Pz_Px_C1001000001_a;
  abcd[115] = 2.0E0*I_ERI_Dxz_S_Pz_Px_C1001000001_a;
  abcd[145] = 2.0E0*I_ERI_D2x_S_Px_Py_C1001000001_a-1*I_ERI_S_S_Px_Py_C1001000001;
  abcd[146] = 2.0E0*I_ERI_Dxy_S_Px_Py_C1001000001_a;
  abcd[147] = 2.0E0*I_ERI_Dxz_S_Px_Py_C1001000001_a;
  abcd[161] = 2.0E0*I_ERI_D2x_S_Py_Py_C1001000001_a-1*I_ERI_S_S_Py_Py_C1001000001;
  abcd[162] = 2.0E0*I_ERI_Dxy_S_Py_Py_C1001000001_a;
  abcd[163] = 2.0E0*I_ERI_Dxz_S_Py_Py_C1001000001_a;
  abcd[177] = 2.0E0*I_ERI_D2x_S_Pz_Py_C1001000001_a-1*I_ERI_S_S_Pz_Py_C1001000001;
  abcd[178] = 2.0E0*I_ERI_Dxy_S_Pz_Py_C1001000001_a;
  abcd[179] = 2.0E0*I_ERI_Dxz_S_Pz_Py_C1001000001_a;
  abcd[209] = 2.0E0*I_ERI_D2x_S_Px_Pz_C1001000001_a-1*I_ERI_S_S_Px_Pz_C1001000001;
  abcd[210] = 2.0E0*I_ERI_Dxy_S_Px_Pz_C1001000001_a;
  abcd[211] = 2.0E0*I_ERI_Dxz_S_Px_Pz_C1001000001_a;
  abcd[225] = 2.0E0*I_ERI_D2x_S_Py_Pz_C1001000001_a-1*I_ERI_S_S_Py_Pz_C1001000001;
  abcd[226] = 2.0E0*I_ERI_Dxy_S_Py_Pz_C1001000001_a;
  abcd[227] = 2.0E0*I_ERI_Dxz_S_Py_Pz_C1001000001_a;
  abcd[241] = 2.0E0*I_ERI_D2x_S_Pz_Pz_C1001000001_a-1*I_ERI_S_S_Pz_Pz_C1001000001;
  abcd[242] = 2.0E0*I_ERI_Dxy_S_Pz_Pz_C1001000001_a;
  abcd[243] = 2.0E0*I_ERI_Dxz_S_Pz_Pz_C1001000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001000_a
   ************************************************************/
  abcd[84] = 2.0E0*I_ERI_Px_Px_Px_Px_C1001001000_a;
  abcd[88] = 2.0E0*I_ERI_Px_Py_Px_Px_C1001001000_a;
  abcd[92] = 2.0E0*I_ERI_Px_Pz_Px_Px_C1001001000_a;
  abcd[100] = 2.0E0*I_ERI_Px_Px_Py_Px_C1001001000_a;
  abcd[104] = 2.0E0*I_ERI_Px_Py_Py_Px_C1001001000_a;
  abcd[108] = 2.0E0*I_ERI_Px_Pz_Py_Px_C1001001000_a;
  abcd[116] = 2.0E0*I_ERI_Px_Px_Pz_Px_C1001001000_a;
  abcd[120] = 2.0E0*I_ERI_Px_Py_Pz_Px_C1001001000_a;
  abcd[124] = 2.0E0*I_ERI_Px_Pz_Pz_Px_C1001001000_a;
  abcd[148] = 2.0E0*I_ERI_Px_Px_Px_Py_C1001001000_a;
  abcd[152] = 2.0E0*I_ERI_Px_Py_Px_Py_C1001001000_a;
  abcd[156] = 2.0E0*I_ERI_Px_Pz_Px_Py_C1001001000_a;
  abcd[164] = 2.0E0*I_ERI_Px_Px_Py_Py_C1001001000_a;
  abcd[168] = 2.0E0*I_ERI_Px_Py_Py_Py_C1001001000_a;
  abcd[172] = 2.0E0*I_ERI_Px_Pz_Py_Py_C1001001000_a;
  abcd[180] = 2.0E0*I_ERI_Px_Px_Pz_Py_C1001001000_a;
  abcd[184] = 2.0E0*I_ERI_Px_Py_Pz_Py_C1001001000_a;
  abcd[188] = 2.0E0*I_ERI_Px_Pz_Pz_Py_C1001001000_a;
  abcd[212] = 2.0E0*I_ERI_Px_Px_Px_Pz_C1001001000_a;
  abcd[216] = 2.0E0*I_ERI_Px_Py_Px_Pz_C1001001000_a;
  abcd[220] = 2.0E0*I_ERI_Px_Pz_Px_Pz_C1001001000_a;
  abcd[228] = 2.0E0*I_ERI_Px_Px_Py_Pz_C1001001000_a;
  abcd[232] = 2.0E0*I_ERI_Px_Py_Py_Pz_C1001001000_a;
  abcd[236] = 2.0E0*I_ERI_Px_Pz_Py_Pz_C1001001000_a;
  abcd[244] = 2.0E0*I_ERI_Px_Px_Pz_Pz_C1001001000_a;
  abcd[248] = 2.0E0*I_ERI_Px_Py_Pz_Pz_C1001001000_a;
  abcd[252] = 2.0E0*I_ERI_Px_Pz_Pz_Pz_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001001
   ************************************************************/
  abcd[85] = 2.0E0*I_ERI_D2x_Px_Px_Px_C1001001001_a-1*I_ERI_S_Px_Px_Px_C1001001001;
  abcd[86] = 2.0E0*I_ERI_Dxy_Px_Px_Px_C1001001001_a;
  abcd[87] = 2.0E0*I_ERI_Dxz_Px_Px_Px_C1001001001_a;
  abcd[89] = 2.0E0*I_ERI_D2x_Py_Px_Px_C1001001001_a-1*I_ERI_S_Py_Px_Px_C1001001001;
  abcd[90] = 2.0E0*I_ERI_Dxy_Py_Px_Px_C1001001001_a;
  abcd[91] = 2.0E0*I_ERI_Dxz_Py_Px_Px_C1001001001_a;
  abcd[93] = 2.0E0*I_ERI_D2x_Pz_Px_Px_C1001001001_a-1*I_ERI_S_Pz_Px_Px_C1001001001;
  abcd[94] = 2.0E0*I_ERI_Dxy_Pz_Px_Px_C1001001001_a;
  abcd[95] = 2.0E0*I_ERI_Dxz_Pz_Px_Px_C1001001001_a;
  abcd[101] = 2.0E0*I_ERI_D2x_Px_Py_Px_C1001001001_a-1*I_ERI_S_Px_Py_Px_C1001001001;
  abcd[102] = 2.0E0*I_ERI_Dxy_Px_Py_Px_C1001001001_a;
  abcd[103] = 2.0E0*I_ERI_Dxz_Px_Py_Px_C1001001001_a;
  abcd[105] = 2.0E0*I_ERI_D2x_Py_Py_Px_C1001001001_a-1*I_ERI_S_Py_Py_Px_C1001001001;
  abcd[106] = 2.0E0*I_ERI_Dxy_Py_Py_Px_C1001001001_a;
  abcd[107] = 2.0E0*I_ERI_Dxz_Py_Py_Px_C1001001001_a;
  abcd[109] = 2.0E0*I_ERI_D2x_Pz_Py_Px_C1001001001_a-1*I_ERI_S_Pz_Py_Px_C1001001001;
  abcd[110] = 2.0E0*I_ERI_Dxy_Pz_Py_Px_C1001001001_a;
  abcd[111] = 2.0E0*I_ERI_Dxz_Pz_Py_Px_C1001001001_a;
  abcd[117] = 2.0E0*I_ERI_D2x_Px_Pz_Px_C1001001001_a-1*I_ERI_S_Px_Pz_Px_C1001001001;
  abcd[118] = 2.0E0*I_ERI_Dxy_Px_Pz_Px_C1001001001_a;
  abcd[119] = 2.0E0*I_ERI_Dxz_Px_Pz_Px_C1001001001_a;
  abcd[121] = 2.0E0*I_ERI_D2x_Py_Pz_Px_C1001001001_a-1*I_ERI_S_Py_Pz_Px_C1001001001;
  abcd[122] = 2.0E0*I_ERI_Dxy_Py_Pz_Px_C1001001001_a;
  abcd[123] = 2.0E0*I_ERI_Dxz_Py_Pz_Px_C1001001001_a;
  abcd[125] = 2.0E0*I_ERI_D2x_Pz_Pz_Px_C1001001001_a-1*I_ERI_S_Pz_Pz_Px_C1001001001;
  abcd[126] = 2.0E0*I_ERI_Dxy_Pz_Pz_Px_C1001001001_a;
  abcd[127] = 2.0E0*I_ERI_Dxz_Pz_Pz_Px_C1001001001_a;
  abcd[149] = 2.0E0*I_ERI_D2x_Px_Px_Py_C1001001001_a-1*I_ERI_S_Px_Px_Py_C1001001001;
  abcd[150] = 2.0E0*I_ERI_Dxy_Px_Px_Py_C1001001001_a;
  abcd[151] = 2.0E0*I_ERI_Dxz_Px_Px_Py_C1001001001_a;
  abcd[153] = 2.0E0*I_ERI_D2x_Py_Px_Py_C1001001001_a-1*I_ERI_S_Py_Px_Py_C1001001001;
  abcd[154] = 2.0E0*I_ERI_Dxy_Py_Px_Py_C1001001001_a;
  abcd[155] = 2.0E0*I_ERI_Dxz_Py_Px_Py_C1001001001_a;
  abcd[157] = 2.0E0*I_ERI_D2x_Pz_Px_Py_C1001001001_a-1*I_ERI_S_Pz_Px_Py_C1001001001;
  abcd[158] = 2.0E0*I_ERI_Dxy_Pz_Px_Py_C1001001001_a;
  abcd[159] = 2.0E0*I_ERI_Dxz_Pz_Px_Py_C1001001001_a;
  abcd[165] = 2.0E0*I_ERI_D2x_Px_Py_Py_C1001001001_a-1*I_ERI_S_Px_Py_Py_C1001001001;
  abcd[166] = 2.0E0*I_ERI_Dxy_Px_Py_Py_C1001001001_a;
  abcd[167] = 2.0E0*I_ERI_Dxz_Px_Py_Py_C1001001001_a;
  abcd[169] = 2.0E0*I_ERI_D2x_Py_Py_Py_C1001001001_a-1*I_ERI_S_Py_Py_Py_C1001001001;
  abcd[170] = 2.0E0*I_ERI_Dxy_Py_Py_Py_C1001001001_a;
  abcd[171] = 2.0E0*I_ERI_Dxz_Py_Py_Py_C1001001001_a;
  abcd[173] = 2.0E0*I_ERI_D2x_Pz_Py_Py_C1001001001_a-1*I_ERI_S_Pz_Py_Py_C1001001001;
  abcd[174] = 2.0E0*I_ERI_Dxy_Pz_Py_Py_C1001001001_a;
  abcd[175] = 2.0E0*I_ERI_Dxz_Pz_Py_Py_C1001001001_a;
  abcd[181] = 2.0E0*I_ERI_D2x_Px_Pz_Py_C1001001001_a-1*I_ERI_S_Px_Pz_Py_C1001001001;
  abcd[182] = 2.0E0*I_ERI_Dxy_Px_Pz_Py_C1001001001_a;
  abcd[183] = 2.0E0*I_ERI_Dxz_Px_Pz_Py_C1001001001_a;
  abcd[185] = 2.0E0*I_ERI_D2x_Py_Pz_Py_C1001001001_a-1*I_ERI_S_Py_Pz_Py_C1001001001;
  abcd[186] = 2.0E0*I_ERI_Dxy_Py_Pz_Py_C1001001001_a;
  abcd[187] = 2.0E0*I_ERI_Dxz_Py_Pz_Py_C1001001001_a;
  abcd[189] = 2.0E0*I_ERI_D2x_Pz_Pz_Py_C1001001001_a-1*I_ERI_S_Pz_Pz_Py_C1001001001;
  abcd[190] = 2.0E0*I_ERI_Dxy_Pz_Pz_Py_C1001001001_a;
  abcd[191] = 2.0E0*I_ERI_Dxz_Pz_Pz_Py_C1001001001_a;
  abcd[213] = 2.0E0*I_ERI_D2x_Px_Px_Pz_C1001001001_a-1*I_ERI_S_Px_Px_Pz_C1001001001;
  abcd[214] = 2.0E0*I_ERI_Dxy_Px_Px_Pz_C1001001001_a;
  abcd[215] = 2.0E0*I_ERI_Dxz_Px_Px_Pz_C1001001001_a;
  abcd[217] = 2.0E0*I_ERI_D2x_Py_Px_Pz_C1001001001_a-1*I_ERI_S_Py_Px_Pz_C1001001001;
  abcd[218] = 2.0E0*I_ERI_Dxy_Py_Px_Pz_C1001001001_a;
  abcd[219] = 2.0E0*I_ERI_Dxz_Py_Px_Pz_C1001001001_a;
  abcd[221] = 2.0E0*I_ERI_D2x_Pz_Px_Pz_C1001001001_a-1*I_ERI_S_Pz_Px_Pz_C1001001001;
  abcd[222] = 2.0E0*I_ERI_Dxy_Pz_Px_Pz_C1001001001_a;
  abcd[223] = 2.0E0*I_ERI_Dxz_Pz_Px_Pz_C1001001001_a;
  abcd[229] = 2.0E0*I_ERI_D2x_Px_Py_Pz_C1001001001_a-1*I_ERI_S_Px_Py_Pz_C1001001001;
  abcd[230] = 2.0E0*I_ERI_Dxy_Px_Py_Pz_C1001001001_a;
  abcd[231] = 2.0E0*I_ERI_Dxz_Px_Py_Pz_C1001001001_a;
  abcd[233] = 2.0E0*I_ERI_D2x_Py_Py_Pz_C1001001001_a-1*I_ERI_S_Py_Py_Pz_C1001001001;
  abcd[234] = 2.0E0*I_ERI_Dxy_Py_Py_Pz_C1001001001_a;
  abcd[235] = 2.0E0*I_ERI_Dxz_Py_Py_Pz_C1001001001_a;
  abcd[237] = 2.0E0*I_ERI_D2x_Pz_Py_Pz_C1001001001_a-1*I_ERI_S_Pz_Py_Pz_C1001001001;
  abcd[238] = 2.0E0*I_ERI_Dxy_Pz_Py_Pz_C1001001001_a;
  abcd[239] = 2.0E0*I_ERI_Dxz_Pz_Py_Pz_C1001001001_a;
  abcd[245] = 2.0E0*I_ERI_D2x_Px_Pz_Pz_C1001001001_a-1*I_ERI_S_Px_Pz_Pz_C1001001001;
  abcd[246] = 2.0E0*I_ERI_Dxy_Px_Pz_Pz_C1001001001_a;
  abcd[247] = 2.0E0*I_ERI_Dxz_Px_Pz_Pz_C1001001001_a;
  abcd[249] = 2.0E0*I_ERI_D2x_Py_Pz_Pz_C1001001001_a-1*I_ERI_S_Py_Pz_Pz_C1001001001;
  abcd[250] = 2.0E0*I_ERI_Dxy_Py_Pz_Pz_C1001001001_a;
  abcd[251] = 2.0E0*I_ERI_Dxz_Py_Pz_Pz_C1001001001_a;
  abcd[253] = 2.0E0*I_ERI_D2x_Pz_Pz_Pz_C1001001001_a-1*I_ERI_S_Pz_Pz_Pz_C1001001001;
  abcd[254] = 2.0E0*I_ERI_Dxy_Pz_Pz_Pz_C1001001001_a;
  abcd[255] = 2.0E0*I_ERI_Dxz_Pz_Pz_Pz_C1001001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C0_a
   ************************************************************/
  abcd[256] = 2.0E0*I_ERI_Py_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1
   ************************************************************/
  abcd[257] = 2.0E0*I_ERI_Dxy_S_S_S_C1_a;
  abcd[258] = 2.0E0*I_ERI_D2y_S_S_S_C1_a-1*I_ERI_S_S_S_S_C1;
  abcd[259] = 2.0E0*I_ERI_Dyz_S_S_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1000_a
   ************************************************************/
  abcd[260] = 2.0E0*I_ERI_Py_Px_S_S_C1000_a;
  abcd[264] = 2.0E0*I_ERI_Py_Py_S_S_C1000_a;
  abcd[268] = 2.0E0*I_ERI_Py_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001
   ************************************************************/
  abcd[261] = 2.0E0*I_ERI_Dxy_Px_S_S_C1001_a;
  abcd[262] = 2.0E0*I_ERI_D2y_Px_S_S_C1001_a-1*I_ERI_S_Px_S_S_C1001;
  abcd[263] = 2.0E0*I_ERI_Dyz_Px_S_S_C1001_a;
  abcd[265] = 2.0E0*I_ERI_Dxy_Py_S_S_C1001_a;
  abcd[266] = 2.0E0*I_ERI_D2y_Py_S_S_C1001_a-1*I_ERI_S_Py_S_S_C1001;
  abcd[267] = 2.0E0*I_ERI_Dyz_Py_S_S_C1001_a;
  abcd[269] = 2.0E0*I_ERI_Dxy_Pz_S_S_C1001_a;
  abcd[270] = 2.0E0*I_ERI_D2y_Pz_S_S_C1001_a-1*I_ERI_S_Pz_S_S_C1001;
  abcd[271] = 2.0E0*I_ERI_Dyz_Pz_S_S_C1001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[272] = 2.0E0*I_ERI_Py_S_Px_S_C1000000_a;
  abcd[288] = 2.0E0*I_ERI_Py_S_Py_S_C1000000_a;
  abcd[304] = 2.0E0*I_ERI_Py_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[273] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[274] = 2.0E0*I_ERI_D2y_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[275] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[289] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[290] = 2.0E0*I_ERI_D2y_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[291] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[305] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[306] = 2.0E0*I_ERI_D2y_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[307] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[276] = 2.0E0*I_ERI_Py_Px_Px_S_C1001000_a;
  abcd[280] = 2.0E0*I_ERI_Py_Py_Px_S_C1001000_a;
  abcd[284] = 2.0E0*I_ERI_Py_Pz_Px_S_C1001000_a;
  abcd[292] = 2.0E0*I_ERI_Py_Px_Py_S_C1001000_a;
  abcd[296] = 2.0E0*I_ERI_Py_Py_Py_S_C1001000_a;
  abcd[300] = 2.0E0*I_ERI_Py_Pz_Py_S_C1001000_a;
  abcd[308] = 2.0E0*I_ERI_Py_Px_Pz_S_C1001000_a;
  abcd[312] = 2.0E0*I_ERI_Py_Py_Pz_S_C1001000_a;
  abcd[316] = 2.0E0*I_ERI_Py_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[277] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[278] = 2.0E0*I_ERI_D2y_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[279] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[281] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[282] = 2.0E0*I_ERI_D2y_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[283] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[285] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[286] = 2.0E0*I_ERI_D2y_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[287] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[293] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[294] = 2.0E0*I_ERI_D2y_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[295] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[297] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[298] = 2.0E0*I_ERI_D2y_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[299] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[301] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[302] = 2.0E0*I_ERI_D2y_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[303] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[309] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[310] = 2.0E0*I_ERI_D2y_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[311] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[313] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[314] = 2.0E0*I_ERI_D2y_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[315] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[317] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[318] = 2.0E0*I_ERI_D2y_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[319] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000000000_a
   ************************************************************/
  abcd[320] = 2.0E0*I_ERI_Py_S_S_Px_C1000000000_a;
  abcd[384] = 2.0E0*I_ERI_Py_S_S_Py_C1000000000_a;
  abcd[448] = 2.0E0*I_ERI_Py_S_S_Pz_C1000000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000001_a
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000000001
   ************************************************************/
  abcd[321] = 2.0E0*I_ERI_Dxy_S_S_Px_C1000000001_a;
  abcd[322] = 2.0E0*I_ERI_D2y_S_S_Px_C1000000001_a-1*I_ERI_S_S_S_Px_C1000000001;
  abcd[323] = 2.0E0*I_ERI_Dyz_S_S_Px_C1000000001_a;
  abcd[385] = 2.0E0*I_ERI_Dxy_S_S_Py_C1000000001_a;
  abcd[386] = 2.0E0*I_ERI_D2y_S_S_Py_C1000000001_a-1*I_ERI_S_S_S_Py_C1000000001;
  abcd[387] = 2.0E0*I_ERI_Dyz_S_S_Py_C1000000001_a;
  abcd[449] = 2.0E0*I_ERI_Dxy_S_S_Pz_C1000000001_a;
  abcd[450] = 2.0E0*I_ERI_D2y_S_S_Pz_C1000000001_a-1*I_ERI_S_S_S_Pz_C1000000001;
  abcd[451] = 2.0E0*I_ERI_Dyz_S_S_Pz_C1000000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000001000_a
   ************************************************************/
  abcd[324] = 2.0E0*I_ERI_Py_Px_S_Px_C1000001000_a;
  abcd[328] = 2.0E0*I_ERI_Py_Py_S_Px_C1000001000_a;
  abcd[332] = 2.0E0*I_ERI_Py_Pz_S_Px_C1000001000_a;
  abcd[388] = 2.0E0*I_ERI_Py_Px_S_Py_C1000001000_a;
  abcd[392] = 2.0E0*I_ERI_Py_Py_S_Py_C1000001000_a;
  abcd[396] = 2.0E0*I_ERI_Py_Pz_S_Py_C1000001000_a;
  abcd[452] = 2.0E0*I_ERI_Py_Px_S_Pz_C1000001000_a;
  abcd[456] = 2.0E0*I_ERI_Py_Py_S_Pz_C1000001000_a;
  abcd[460] = 2.0E0*I_ERI_Py_Pz_S_Pz_C1000001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_P_C1000001001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000001001
   ************************************************************/
  abcd[325] = 2.0E0*I_ERI_Dxy_Px_S_Px_C1000001001_a;
  abcd[326] = 2.0E0*I_ERI_D2y_Px_S_Px_C1000001001_a-1*I_ERI_S_Px_S_Px_C1000001001;
  abcd[327] = 2.0E0*I_ERI_Dyz_Px_S_Px_C1000001001_a;
  abcd[329] = 2.0E0*I_ERI_Dxy_Py_S_Px_C1000001001_a;
  abcd[330] = 2.0E0*I_ERI_D2y_Py_S_Px_C1000001001_a-1*I_ERI_S_Py_S_Px_C1000001001;
  abcd[331] = 2.0E0*I_ERI_Dyz_Py_S_Px_C1000001001_a;
  abcd[333] = 2.0E0*I_ERI_Dxy_Pz_S_Px_C1000001001_a;
  abcd[334] = 2.0E0*I_ERI_D2y_Pz_S_Px_C1000001001_a-1*I_ERI_S_Pz_S_Px_C1000001001;
  abcd[335] = 2.0E0*I_ERI_Dyz_Pz_S_Px_C1000001001_a;
  abcd[389] = 2.0E0*I_ERI_Dxy_Px_S_Py_C1000001001_a;
  abcd[390] = 2.0E0*I_ERI_D2y_Px_S_Py_C1000001001_a-1*I_ERI_S_Px_S_Py_C1000001001;
  abcd[391] = 2.0E0*I_ERI_Dyz_Px_S_Py_C1000001001_a;
  abcd[393] = 2.0E0*I_ERI_Dxy_Py_S_Py_C1000001001_a;
  abcd[394] = 2.0E0*I_ERI_D2y_Py_S_Py_C1000001001_a-1*I_ERI_S_Py_S_Py_C1000001001;
  abcd[395] = 2.0E0*I_ERI_Dyz_Py_S_Py_C1000001001_a;
  abcd[397] = 2.0E0*I_ERI_Dxy_Pz_S_Py_C1000001001_a;
  abcd[398] = 2.0E0*I_ERI_D2y_Pz_S_Py_C1000001001_a-1*I_ERI_S_Pz_S_Py_C1000001001;
  abcd[399] = 2.0E0*I_ERI_Dyz_Pz_S_Py_C1000001001_a;
  abcd[453] = 2.0E0*I_ERI_Dxy_Px_S_Pz_C1000001001_a;
  abcd[454] = 2.0E0*I_ERI_D2y_Px_S_Pz_C1000001001_a-1*I_ERI_S_Px_S_Pz_C1000001001;
  abcd[455] = 2.0E0*I_ERI_Dyz_Px_S_Pz_C1000001001_a;
  abcd[457] = 2.0E0*I_ERI_Dxy_Py_S_Pz_C1000001001_a;
  abcd[458] = 2.0E0*I_ERI_D2y_Py_S_Pz_C1000001001_a-1*I_ERI_S_Py_S_Pz_C1000001001;
  abcd[459] = 2.0E0*I_ERI_Dyz_Py_S_Pz_C1000001001_a;
  abcd[461] = 2.0E0*I_ERI_Dxy_Pz_S_Pz_C1000001001_a;
  abcd[462] = 2.0E0*I_ERI_D2y_Pz_S_Pz_C1000001001_a-1*I_ERI_S_Pz_S_Pz_C1000001001;
  abcd[463] = 2.0E0*I_ERI_Dyz_Pz_S_Pz_C1000001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001000000_a
   ************************************************************/
  abcd[336] = 2.0E0*I_ERI_Py_S_Px_Px_C1001000000_a;
  abcd[352] = 2.0E0*I_ERI_Py_S_Py_Px_C1001000000_a;
  abcd[368] = 2.0E0*I_ERI_Py_S_Pz_Px_C1001000000_a;
  abcd[400] = 2.0E0*I_ERI_Py_S_Px_Py_C1001000000_a;
  abcd[416] = 2.0E0*I_ERI_Py_S_Py_Py_C1001000000_a;
  abcd[432] = 2.0E0*I_ERI_Py_S_Pz_Py_C1001000000_a;
  abcd[464] = 2.0E0*I_ERI_Py_S_Px_Pz_C1001000000_a;
  abcd[480] = 2.0E0*I_ERI_Py_S_Py_Pz_C1001000000_a;
  abcd[496] = 2.0E0*I_ERI_Py_S_Pz_Pz_C1001000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001000001
   ************************************************************/
  abcd[337] = 2.0E0*I_ERI_Dxy_S_Px_Px_C1001000001_a;
  abcd[338] = 2.0E0*I_ERI_D2y_S_Px_Px_C1001000001_a-1*I_ERI_S_S_Px_Px_C1001000001;
  abcd[339] = 2.0E0*I_ERI_Dyz_S_Px_Px_C1001000001_a;
  abcd[353] = 2.0E0*I_ERI_Dxy_S_Py_Px_C1001000001_a;
  abcd[354] = 2.0E0*I_ERI_D2y_S_Py_Px_C1001000001_a-1*I_ERI_S_S_Py_Px_C1001000001;
  abcd[355] = 2.0E0*I_ERI_Dyz_S_Py_Px_C1001000001_a;
  abcd[369] = 2.0E0*I_ERI_Dxy_S_Pz_Px_C1001000001_a;
  abcd[370] = 2.0E0*I_ERI_D2y_S_Pz_Px_C1001000001_a-1*I_ERI_S_S_Pz_Px_C1001000001;
  abcd[371] = 2.0E0*I_ERI_Dyz_S_Pz_Px_C1001000001_a;
  abcd[401] = 2.0E0*I_ERI_Dxy_S_Px_Py_C1001000001_a;
  abcd[402] = 2.0E0*I_ERI_D2y_S_Px_Py_C1001000001_a-1*I_ERI_S_S_Px_Py_C1001000001;
  abcd[403] = 2.0E0*I_ERI_Dyz_S_Px_Py_C1001000001_a;
  abcd[417] = 2.0E0*I_ERI_Dxy_S_Py_Py_C1001000001_a;
  abcd[418] = 2.0E0*I_ERI_D2y_S_Py_Py_C1001000001_a-1*I_ERI_S_S_Py_Py_C1001000001;
  abcd[419] = 2.0E0*I_ERI_Dyz_S_Py_Py_C1001000001_a;
  abcd[433] = 2.0E0*I_ERI_Dxy_S_Pz_Py_C1001000001_a;
  abcd[434] = 2.0E0*I_ERI_D2y_S_Pz_Py_C1001000001_a-1*I_ERI_S_S_Pz_Py_C1001000001;
  abcd[435] = 2.0E0*I_ERI_Dyz_S_Pz_Py_C1001000001_a;
  abcd[465] = 2.0E0*I_ERI_Dxy_S_Px_Pz_C1001000001_a;
  abcd[466] = 2.0E0*I_ERI_D2y_S_Px_Pz_C1001000001_a-1*I_ERI_S_S_Px_Pz_C1001000001;
  abcd[467] = 2.0E0*I_ERI_Dyz_S_Px_Pz_C1001000001_a;
  abcd[481] = 2.0E0*I_ERI_Dxy_S_Py_Pz_C1001000001_a;
  abcd[482] = 2.0E0*I_ERI_D2y_S_Py_Pz_C1001000001_a-1*I_ERI_S_S_Py_Pz_C1001000001;
  abcd[483] = 2.0E0*I_ERI_Dyz_S_Py_Pz_C1001000001_a;
  abcd[497] = 2.0E0*I_ERI_Dxy_S_Pz_Pz_C1001000001_a;
  abcd[498] = 2.0E0*I_ERI_D2y_S_Pz_Pz_C1001000001_a-1*I_ERI_S_S_Pz_Pz_C1001000001;
  abcd[499] = 2.0E0*I_ERI_Dyz_S_Pz_Pz_C1001000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001000_a
   ************************************************************/
  abcd[340] = 2.0E0*I_ERI_Py_Px_Px_Px_C1001001000_a;
  abcd[344] = 2.0E0*I_ERI_Py_Py_Px_Px_C1001001000_a;
  abcd[348] = 2.0E0*I_ERI_Py_Pz_Px_Px_C1001001000_a;
  abcd[356] = 2.0E0*I_ERI_Py_Px_Py_Px_C1001001000_a;
  abcd[360] = 2.0E0*I_ERI_Py_Py_Py_Px_C1001001000_a;
  abcd[364] = 2.0E0*I_ERI_Py_Pz_Py_Px_C1001001000_a;
  abcd[372] = 2.0E0*I_ERI_Py_Px_Pz_Px_C1001001000_a;
  abcd[376] = 2.0E0*I_ERI_Py_Py_Pz_Px_C1001001000_a;
  abcd[380] = 2.0E0*I_ERI_Py_Pz_Pz_Px_C1001001000_a;
  abcd[404] = 2.0E0*I_ERI_Py_Px_Px_Py_C1001001000_a;
  abcd[408] = 2.0E0*I_ERI_Py_Py_Px_Py_C1001001000_a;
  abcd[412] = 2.0E0*I_ERI_Py_Pz_Px_Py_C1001001000_a;
  abcd[420] = 2.0E0*I_ERI_Py_Px_Py_Py_C1001001000_a;
  abcd[424] = 2.0E0*I_ERI_Py_Py_Py_Py_C1001001000_a;
  abcd[428] = 2.0E0*I_ERI_Py_Pz_Py_Py_C1001001000_a;
  abcd[436] = 2.0E0*I_ERI_Py_Px_Pz_Py_C1001001000_a;
  abcd[440] = 2.0E0*I_ERI_Py_Py_Pz_Py_C1001001000_a;
  abcd[444] = 2.0E0*I_ERI_Py_Pz_Pz_Py_C1001001000_a;
  abcd[468] = 2.0E0*I_ERI_Py_Px_Px_Pz_C1001001000_a;
  abcd[472] = 2.0E0*I_ERI_Py_Py_Px_Pz_C1001001000_a;
  abcd[476] = 2.0E0*I_ERI_Py_Pz_Px_Pz_C1001001000_a;
  abcd[484] = 2.0E0*I_ERI_Py_Px_Py_Pz_C1001001000_a;
  abcd[488] = 2.0E0*I_ERI_Py_Py_Py_Pz_C1001001000_a;
  abcd[492] = 2.0E0*I_ERI_Py_Pz_Py_Pz_C1001001000_a;
  abcd[500] = 2.0E0*I_ERI_Py_Px_Pz_Pz_C1001001000_a;
  abcd[504] = 2.0E0*I_ERI_Py_Py_Pz_Pz_C1001001000_a;
  abcd[508] = 2.0E0*I_ERI_Py_Pz_Pz_Pz_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001001
   ************************************************************/
  abcd[341] = 2.0E0*I_ERI_Dxy_Px_Px_Px_C1001001001_a;
  abcd[342] = 2.0E0*I_ERI_D2y_Px_Px_Px_C1001001001_a-1*I_ERI_S_Px_Px_Px_C1001001001;
  abcd[343] = 2.0E0*I_ERI_Dyz_Px_Px_Px_C1001001001_a;
  abcd[345] = 2.0E0*I_ERI_Dxy_Py_Px_Px_C1001001001_a;
  abcd[346] = 2.0E0*I_ERI_D2y_Py_Px_Px_C1001001001_a-1*I_ERI_S_Py_Px_Px_C1001001001;
  abcd[347] = 2.0E0*I_ERI_Dyz_Py_Px_Px_C1001001001_a;
  abcd[349] = 2.0E0*I_ERI_Dxy_Pz_Px_Px_C1001001001_a;
  abcd[350] = 2.0E0*I_ERI_D2y_Pz_Px_Px_C1001001001_a-1*I_ERI_S_Pz_Px_Px_C1001001001;
  abcd[351] = 2.0E0*I_ERI_Dyz_Pz_Px_Px_C1001001001_a;
  abcd[357] = 2.0E0*I_ERI_Dxy_Px_Py_Px_C1001001001_a;
  abcd[358] = 2.0E0*I_ERI_D2y_Px_Py_Px_C1001001001_a-1*I_ERI_S_Px_Py_Px_C1001001001;
  abcd[359] = 2.0E0*I_ERI_Dyz_Px_Py_Px_C1001001001_a;
  abcd[361] = 2.0E0*I_ERI_Dxy_Py_Py_Px_C1001001001_a;
  abcd[362] = 2.0E0*I_ERI_D2y_Py_Py_Px_C1001001001_a-1*I_ERI_S_Py_Py_Px_C1001001001;
  abcd[363] = 2.0E0*I_ERI_Dyz_Py_Py_Px_C1001001001_a;
  abcd[365] = 2.0E0*I_ERI_Dxy_Pz_Py_Px_C1001001001_a;
  abcd[366] = 2.0E0*I_ERI_D2y_Pz_Py_Px_C1001001001_a-1*I_ERI_S_Pz_Py_Px_C1001001001;
  abcd[367] = 2.0E0*I_ERI_Dyz_Pz_Py_Px_C1001001001_a;
  abcd[373] = 2.0E0*I_ERI_Dxy_Px_Pz_Px_C1001001001_a;
  abcd[374] = 2.0E0*I_ERI_D2y_Px_Pz_Px_C1001001001_a-1*I_ERI_S_Px_Pz_Px_C1001001001;
  abcd[375] = 2.0E0*I_ERI_Dyz_Px_Pz_Px_C1001001001_a;
  abcd[377] = 2.0E0*I_ERI_Dxy_Py_Pz_Px_C1001001001_a;
  abcd[378] = 2.0E0*I_ERI_D2y_Py_Pz_Px_C1001001001_a-1*I_ERI_S_Py_Pz_Px_C1001001001;
  abcd[379] = 2.0E0*I_ERI_Dyz_Py_Pz_Px_C1001001001_a;
  abcd[381] = 2.0E0*I_ERI_Dxy_Pz_Pz_Px_C1001001001_a;
  abcd[382] = 2.0E0*I_ERI_D2y_Pz_Pz_Px_C1001001001_a-1*I_ERI_S_Pz_Pz_Px_C1001001001;
  abcd[383] = 2.0E0*I_ERI_Dyz_Pz_Pz_Px_C1001001001_a;
  abcd[405] = 2.0E0*I_ERI_Dxy_Px_Px_Py_C1001001001_a;
  abcd[406] = 2.0E0*I_ERI_D2y_Px_Px_Py_C1001001001_a-1*I_ERI_S_Px_Px_Py_C1001001001;
  abcd[407] = 2.0E0*I_ERI_Dyz_Px_Px_Py_C1001001001_a;
  abcd[409] = 2.0E0*I_ERI_Dxy_Py_Px_Py_C1001001001_a;
  abcd[410] = 2.0E0*I_ERI_D2y_Py_Px_Py_C1001001001_a-1*I_ERI_S_Py_Px_Py_C1001001001;
  abcd[411] = 2.0E0*I_ERI_Dyz_Py_Px_Py_C1001001001_a;
  abcd[413] = 2.0E0*I_ERI_Dxy_Pz_Px_Py_C1001001001_a;
  abcd[414] = 2.0E0*I_ERI_D2y_Pz_Px_Py_C1001001001_a-1*I_ERI_S_Pz_Px_Py_C1001001001;
  abcd[415] = 2.0E0*I_ERI_Dyz_Pz_Px_Py_C1001001001_a;
  abcd[421] = 2.0E0*I_ERI_Dxy_Px_Py_Py_C1001001001_a;
  abcd[422] = 2.0E0*I_ERI_D2y_Px_Py_Py_C1001001001_a-1*I_ERI_S_Px_Py_Py_C1001001001;
  abcd[423] = 2.0E0*I_ERI_Dyz_Px_Py_Py_C1001001001_a;
  abcd[425] = 2.0E0*I_ERI_Dxy_Py_Py_Py_C1001001001_a;
  abcd[426] = 2.0E0*I_ERI_D2y_Py_Py_Py_C1001001001_a-1*I_ERI_S_Py_Py_Py_C1001001001;
  abcd[427] = 2.0E0*I_ERI_Dyz_Py_Py_Py_C1001001001_a;
  abcd[429] = 2.0E0*I_ERI_Dxy_Pz_Py_Py_C1001001001_a;
  abcd[430] = 2.0E0*I_ERI_D2y_Pz_Py_Py_C1001001001_a-1*I_ERI_S_Pz_Py_Py_C1001001001;
  abcd[431] = 2.0E0*I_ERI_Dyz_Pz_Py_Py_C1001001001_a;
  abcd[437] = 2.0E0*I_ERI_Dxy_Px_Pz_Py_C1001001001_a;
  abcd[438] = 2.0E0*I_ERI_D2y_Px_Pz_Py_C1001001001_a-1*I_ERI_S_Px_Pz_Py_C1001001001;
  abcd[439] = 2.0E0*I_ERI_Dyz_Px_Pz_Py_C1001001001_a;
  abcd[441] = 2.0E0*I_ERI_Dxy_Py_Pz_Py_C1001001001_a;
  abcd[442] = 2.0E0*I_ERI_D2y_Py_Pz_Py_C1001001001_a-1*I_ERI_S_Py_Pz_Py_C1001001001;
  abcd[443] = 2.0E0*I_ERI_Dyz_Py_Pz_Py_C1001001001_a;
  abcd[445] = 2.0E0*I_ERI_Dxy_Pz_Pz_Py_C1001001001_a;
  abcd[446] = 2.0E0*I_ERI_D2y_Pz_Pz_Py_C1001001001_a-1*I_ERI_S_Pz_Pz_Py_C1001001001;
  abcd[447] = 2.0E0*I_ERI_Dyz_Pz_Pz_Py_C1001001001_a;
  abcd[469] = 2.0E0*I_ERI_Dxy_Px_Px_Pz_C1001001001_a;
  abcd[470] = 2.0E0*I_ERI_D2y_Px_Px_Pz_C1001001001_a-1*I_ERI_S_Px_Px_Pz_C1001001001;
  abcd[471] = 2.0E0*I_ERI_Dyz_Px_Px_Pz_C1001001001_a;
  abcd[473] = 2.0E0*I_ERI_Dxy_Py_Px_Pz_C1001001001_a;
  abcd[474] = 2.0E0*I_ERI_D2y_Py_Px_Pz_C1001001001_a-1*I_ERI_S_Py_Px_Pz_C1001001001;
  abcd[475] = 2.0E0*I_ERI_Dyz_Py_Px_Pz_C1001001001_a;
  abcd[477] = 2.0E0*I_ERI_Dxy_Pz_Px_Pz_C1001001001_a;
  abcd[478] = 2.0E0*I_ERI_D2y_Pz_Px_Pz_C1001001001_a-1*I_ERI_S_Pz_Px_Pz_C1001001001;
  abcd[479] = 2.0E0*I_ERI_Dyz_Pz_Px_Pz_C1001001001_a;
  abcd[485] = 2.0E0*I_ERI_Dxy_Px_Py_Pz_C1001001001_a;
  abcd[486] = 2.0E0*I_ERI_D2y_Px_Py_Pz_C1001001001_a-1*I_ERI_S_Px_Py_Pz_C1001001001;
  abcd[487] = 2.0E0*I_ERI_Dyz_Px_Py_Pz_C1001001001_a;
  abcd[489] = 2.0E0*I_ERI_Dxy_Py_Py_Pz_C1001001001_a;
  abcd[490] = 2.0E0*I_ERI_D2y_Py_Py_Pz_C1001001001_a-1*I_ERI_S_Py_Py_Pz_C1001001001;
  abcd[491] = 2.0E0*I_ERI_Dyz_Py_Py_Pz_C1001001001_a;
  abcd[493] = 2.0E0*I_ERI_Dxy_Pz_Py_Pz_C1001001001_a;
  abcd[494] = 2.0E0*I_ERI_D2y_Pz_Py_Pz_C1001001001_a-1*I_ERI_S_Pz_Py_Pz_C1001001001;
  abcd[495] = 2.0E0*I_ERI_Dyz_Pz_Py_Pz_C1001001001_a;
  abcd[501] = 2.0E0*I_ERI_Dxy_Px_Pz_Pz_C1001001001_a;
  abcd[502] = 2.0E0*I_ERI_D2y_Px_Pz_Pz_C1001001001_a-1*I_ERI_S_Px_Pz_Pz_C1001001001;
  abcd[503] = 2.0E0*I_ERI_Dyz_Px_Pz_Pz_C1001001001_a;
  abcd[505] = 2.0E0*I_ERI_Dxy_Py_Pz_Pz_C1001001001_a;
  abcd[506] = 2.0E0*I_ERI_D2y_Py_Pz_Pz_C1001001001_a-1*I_ERI_S_Py_Pz_Pz_C1001001001;
  abcd[507] = 2.0E0*I_ERI_Dyz_Py_Pz_Pz_C1001001001_a;
  abcd[509] = 2.0E0*I_ERI_Dxy_Pz_Pz_Pz_C1001001001_a;
  abcd[510] = 2.0E0*I_ERI_D2y_Pz_Pz_Pz_C1001001001_a-1*I_ERI_S_Pz_Pz_Pz_C1001001001;
  abcd[511] = 2.0E0*I_ERI_Dyz_Pz_Pz_Pz_C1001001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C0_a
   ************************************************************/
  abcd[512] = 2.0E0*I_ERI_Pz_S_S_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1
   ************************************************************/
  abcd[513] = 2.0E0*I_ERI_Dxz_S_S_S_C1_a;
  abcd[514] = 2.0E0*I_ERI_Dyz_S_S_S_C1_a;
  abcd[515] = 2.0E0*I_ERI_D2z_S_S_S_C1_a-1*I_ERI_S_S_S_S_C1;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1000_a
   ************************************************************/
  abcd[516] = 2.0E0*I_ERI_Pz_Px_S_S_C1000_a;
  abcd[520] = 2.0E0*I_ERI_Pz_Py_S_S_C1000_a;
  abcd[524] = 2.0E0*I_ERI_Pz_Pz_S_S_C1000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001
   ************************************************************/
  abcd[517] = 2.0E0*I_ERI_Dxz_Px_S_S_C1001_a;
  abcd[518] = 2.0E0*I_ERI_Dyz_Px_S_S_C1001_a;
  abcd[519] = 2.0E0*I_ERI_D2z_Px_S_S_C1001_a-1*I_ERI_S_Px_S_S_C1001;
  abcd[521] = 2.0E0*I_ERI_Dxz_Py_S_S_C1001_a;
  abcd[522] = 2.0E0*I_ERI_Dyz_Py_S_S_C1001_a;
  abcd[523] = 2.0E0*I_ERI_D2z_Py_S_S_C1001_a-1*I_ERI_S_Py_S_S_C1001;
  abcd[525] = 2.0E0*I_ERI_Dxz_Pz_S_S_C1001_a;
  abcd[526] = 2.0E0*I_ERI_Dyz_Pz_S_S_C1001_a;
  abcd[527] = 2.0E0*I_ERI_D2z_Pz_S_S_C1001_a-1*I_ERI_S_Pz_S_S_C1001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[528] = 2.0E0*I_ERI_Pz_S_Px_S_C1000000_a;
  abcd[544] = 2.0E0*I_ERI_Pz_S_Py_S_C1000000_a;
  abcd[560] = 2.0E0*I_ERI_Pz_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[529] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[530] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[531] = 2.0E0*I_ERI_D2z_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[545] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[546] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[547] = 2.0E0*I_ERI_D2z_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[561] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;
  abcd[562] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;
  abcd[563] = 2.0E0*I_ERI_D2z_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[532] = 2.0E0*I_ERI_Pz_Px_Px_S_C1001000_a;
  abcd[536] = 2.0E0*I_ERI_Pz_Py_Px_S_C1001000_a;
  abcd[540] = 2.0E0*I_ERI_Pz_Pz_Px_S_C1001000_a;
  abcd[548] = 2.0E0*I_ERI_Pz_Px_Py_S_C1001000_a;
  abcd[552] = 2.0E0*I_ERI_Pz_Py_Py_S_C1001000_a;
  abcd[556] = 2.0E0*I_ERI_Pz_Pz_Py_S_C1001000_a;
  abcd[564] = 2.0E0*I_ERI_Pz_Px_Pz_S_C1001000_a;
  abcd[568] = 2.0E0*I_ERI_Pz_Py_Pz_S_C1001000_a;
  abcd[572] = 2.0E0*I_ERI_Pz_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[533] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[534] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[535] = 2.0E0*I_ERI_D2z_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[537] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[538] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[539] = 2.0E0*I_ERI_D2z_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[541] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[542] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[543] = 2.0E0*I_ERI_D2z_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[549] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[550] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[551] = 2.0E0*I_ERI_D2z_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[553] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[554] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[555] = 2.0E0*I_ERI_D2z_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[557] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[558] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[559] = 2.0E0*I_ERI_D2z_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[565] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[566] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[567] = 2.0E0*I_ERI_D2z_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[569] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[570] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[571] = 2.0E0*I_ERI_D2z_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[573] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;
  abcd[574] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;
  abcd[575] = 2.0E0*I_ERI_D2z_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000000000_a
   ************************************************************/
  abcd[576] = 2.0E0*I_ERI_Pz_S_S_Px_C1000000000_a;
  abcd[640] = 2.0E0*I_ERI_Pz_S_S_Py_C1000000000_a;
  abcd[704] = 2.0E0*I_ERI_Pz_S_S_Pz_C1000000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000001_a
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000000001
   ************************************************************/
  abcd[577] = 2.0E0*I_ERI_Dxz_S_S_Px_C1000000001_a;
  abcd[578] = 2.0E0*I_ERI_Dyz_S_S_Px_C1000000001_a;
  abcd[579] = 2.0E0*I_ERI_D2z_S_S_Px_C1000000001_a-1*I_ERI_S_S_S_Px_C1000000001;
  abcd[641] = 2.0E0*I_ERI_Dxz_S_S_Py_C1000000001_a;
  abcd[642] = 2.0E0*I_ERI_Dyz_S_S_Py_C1000000001_a;
  abcd[643] = 2.0E0*I_ERI_D2z_S_S_Py_C1000000001_a-1*I_ERI_S_S_S_Py_C1000000001;
  abcd[705] = 2.0E0*I_ERI_Dxz_S_S_Pz_C1000000001_a;
  abcd[706] = 2.0E0*I_ERI_Dyz_S_S_Pz_C1000000001_a;
  abcd[707] = 2.0E0*I_ERI_D2z_S_S_Pz_C1000000001_a-1*I_ERI_S_S_S_Pz_C1000000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000001000_a
   ************************************************************/
  abcd[580] = 2.0E0*I_ERI_Pz_Px_S_Px_C1000001000_a;
  abcd[584] = 2.0E0*I_ERI_Pz_Py_S_Px_C1000001000_a;
  abcd[588] = 2.0E0*I_ERI_Pz_Pz_S_Px_C1000001000_a;
  abcd[644] = 2.0E0*I_ERI_Pz_Px_S_Py_C1000001000_a;
  abcd[648] = 2.0E0*I_ERI_Pz_Py_S_Py_C1000001000_a;
  abcd[652] = 2.0E0*I_ERI_Pz_Pz_S_Py_C1000001000_a;
  abcd[708] = 2.0E0*I_ERI_Pz_Px_S_Pz_C1000001000_a;
  abcd[712] = 2.0E0*I_ERI_Pz_Py_S_Pz_C1000001000_a;
  abcd[716] = 2.0E0*I_ERI_Pz_Pz_S_Pz_C1000001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_P_C1000001001_a
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000001001
   ************************************************************/
  abcd[581] = 2.0E0*I_ERI_Dxz_Px_S_Px_C1000001001_a;
  abcd[582] = 2.0E0*I_ERI_Dyz_Px_S_Px_C1000001001_a;
  abcd[583] = 2.0E0*I_ERI_D2z_Px_S_Px_C1000001001_a-1*I_ERI_S_Px_S_Px_C1000001001;
  abcd[585] = 2.0E0*I_ERI_Dxz_Py_S_Px_C1000001001_a;
  abcd[586] = 2.0E0*I_ERI_Dyz_Py_S_Px_C1000001001_a;
  abcd[587] = 2.0E0*I_ERI_D2z_Py_S_Px_C1000001001_a-1*I_ERI_S_Py_S_Px_C1000001001;
  abcd[589] = 2.0E0*I_ERI_Dxz_Pz_S_Px_C1000001001_a;
  abcd[590] = 2.0E0*I_ERI_Dyz_Pz_S_Px_C1000001001_a;
  abcd[591] = 2.0E0*I_ERI_D2z_Pz_S_Px_C1000001001_a-1*I_ERI_S_Pz_S_Px_C1000001001;
  abcd[645] = 2.0E0*I_ERI_Dxz_Px_S_Py_C1000001001_a;
  abcd[646] = 2.0E0*I_ERI_Dyz_Px_S_Py_C1000001001_a;
  abcd[647] = 2.0E0*I_ERI_D2z_Px_S_Py_C1000001001_a-1*I_ERI_S_Px_S_Py_C1000001001;
  abcd[649] = 2.0E0*I_ERI_Dxz_Py_S_Py_C1000001001_a;
  abcd[650] = 2.0E0*I_ERI_Dyz_Py_S_Py_C1000001001_a;
  abcd[651] = 2.0E0*I_ERI_D2z_Py_S_Py_C1000001001_a-1*I_ERI_S_Py_S_Py_C1000001001;
  abcd[653] = 2.0E0*I_ERI_Dxz_Pz_S_Py_C1000001001_a;
  abcd[654] = 2.0E0*I_ERI_Dyz_Pz_S_Py_C1000001001_a;
  abcd[655] = 2.0E0*I_ERI_D2z_Pz_S_Py_C1000001001_a-1*I_ERI_S_Pz_S_Py_C1000001001;
  abcd[709] = 2.0E0*I_ERI_Dxz_Px_S_Pz_C1000001001_a;
  abcd[710] = 2.0E0*I_ERI_Dyz_Px_S_Pz_C1000001001_a;
  abcd[711] = 2.0E0*I_ERI_D2z_Px_S_Pz_C1000001001_a-1*I_ERI_S_Px_S_Pz_C1000001001;
  abcd[713] = 2.0E0*I_ERI_Dxz_Py_S_Pz_C1000001001_a;
  abcd[714] = 2.0E0*I_ERI_Dyz_Py_S_Pz_C1000001001_a;
  abcd[715] = 2.0E0*I_ERI_D2z_Py_S_Pz_C1000001001_a-1*I_ERI_S_Py_S_Pz_C1000001001;
  abcd[717] = 2.0E0*I_ERI_Dxz_Pz_S_Pz_C1000001001_a;
  abcd[718] = 2.0E0*I_ERI_Dyz_Pz_S_Pz_C1000001001_a;
  abcd[719] = 2.0E0*I_ERI_D2z_Pz_S_Pz_C1000001001_a-1*I_ERI_S_Pz_S_Pz_C1000001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001000000_a
   ************************************************************/
  abcd[592] = 2.0E0*I_ERI_Pz_S_Px_Px_C1001000000_a;
  abcd[608] = 2.0E0*I_ERI_Pz_S_Py_Px_C1001000000_a;
  abcd[624] = 2.0E0*I_ERI_Pz_S_Pz_Px_C1001000000_a;
  abcd[656] = 2.0E0*I_ERI_Pz_S_Px_Py_C1001000000_a;
  abcd[672] = 2.0E0*I_ERI_Pz_S_Py_Py_C1001000000_a;
  abcd[688] = 2.0E0*I_ERI_Pz_S_Pz_Py_C1001000000_a;
  abcd[720] = 2.0E0*I_ERI_Pz_S_Px_Pz_C1001000000_a;
  abcd[736] = 2.0E0*I_ERI_Pz_S_Py_Pz_C1001000000_a;
  abcd[752] = 2.0E0*I_ERI_Pz_S_Pz_Pz_C1001000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001000001
   ************************************************************/
  abcd[593] = 2.0E0*I_ERI_Dxz_S_Px_Px_C1001000001_a;
  abcd[594] = 2.0E0*I_ERI_Dyz_S_Px_Px_C1001000001_a;
  abcd[595] = 2.0E0*I_ERI_D2z_S_Px_Px_C1001000001_a-1*I_ERI_S_S_Px_Px_C1001000001;
  abcd[609] = 2.0E0*I_ERI_Dxz_S_Py_Px_C1001000001_a;
  abcd[610] = 2.0E0*I_ERI_Dyz_S_Py_Px_C1001000001_a;
  abcd[611] = 2.0E0*I_ERI_D2z_S_Py_Px_C1001000001_a-1*I_ERI_S_S_Py_Px_C1001000001;
  abcd[625] = 2.0E0*I_ERI_Dxz_S_Pz_Px_C1001000001_a;
  abcd[626] = 2.0E0*I_ERI_Dyz_S_Pz_Px_C1001000001_a;
  abcd[627] = 2.0E0*I_ERI_D2z_S_Pz_Px_C1001000001_a-1*I_ERI_S_S_Pz_Px_C1001000001;
  abcd[657] = 2.0E0*I_ERI_Dxz_S_Px_Py_C1001000001_a;
  abcd[658] = 2.0E0*I_ERI_Dyz_S_Px_Py_C1001000001_a;
  abcd[659] = 2.0E0*I_ERI_D2z_S_Px_Py_C1001000001_a-1*I_ERI_S_S_Px_Py_C1001000001;
  abcd[673] = 2.0E0*I_ERI_Dxz_S_Py_Py_C1001000001_a;
  abcd[674] = 2.0E0*I_ERI_Dyz_S_Py_Py_C1001000001_a;
  abcd[675] = 2.0E0*I_ERI_D2z_S_Py_Py_C1001000001_a-1*I_ERI_S_S_Py_Py_C1001000001;
  abcd[689] = 2.0E0*I_ERI_Dxz_S_Pz_Py_C1001000001_a;
  abcd[690] = 2.0E0*I_ERI_Dyz_S_Pz_Py_C1001000001_a;
  abcd[691] = 2.0E0*I_ERI_D2z_S_Pz_Py_C1001000001_a-1*I_ERI_S_S_Pz_Py_C1001000001;
  abcd[721] = 2.0E0*I_ERI_Dxz_S_Px_Pz_C1001000001_a;
  abcd[722] = 2.0E0*I_ERI_Dyz_S_Px_Pz_C1001000001_a;
  abcd[723] = 2.0E0*I_ERI_D2z_S_Px_Pz_C1001000001_a-1*I_ERI_S_S_Px_Pz_C1001000001;
  abcd[737] = 2.0E0*I_ERI_Dxz_S_Py_Pz_C1001000001_a;
  abcd[738] = 2.0E0*I_ERI_Dyz_S_Py_Pz_C1001000001_a;
  abcd[739] = 2.0E0*I_ERI_D2z_S_Py_Pz_C1001000001_a-1*I_ERI_S_S_Py_Pz_C1001000001;
  abcd[753] = 2.0E0*I_ERI_Dxz_S_Pz_Pz_C1001000001_a;
  abcd[754] = 2.0E0*I_ERI_Dyz_S_Pz_Pz_C1001000001_a;
  abcd[755] = 2.0E0*I_ERI_D2z_S_Pz_Pz_C1001000001_a-1*I_ERI_S_S_Pz_Pz_C1001000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001000_a
   ************************************************************/
  abcd[596] = 2.0E0*I_ERI_Pz_Px_Px_Px_C1001001000_a;
  abcd[600] = 2.0E0*I_ERI_Pz_Py_Px_Px_C1001001000_a;
  abcd[604] = 2.0E0*I_ERI_Pz_Pz_Px_Px_C1001001000_a;
  abcd[612] = 2.0E0*I_ERI_Pz_Px_Py_Px_C1001001000_a;
  abcd[616] = 2.0E0*I_ERI_Pz_Py_Py_Px_C1001001000_a;
  abcd[620] = 2.0E0*I_ERI_Pz_Pz_Py_Px_C1001001000_a;
  abcd[628] = 2.0E0*I_ERI_Pz_Px_Pz_Px_C1001001000_a;
  abcd[632] = 2.0E0*I_ERI_Pz_Py_Pz_Px_C1001001000_a;
  abcd[636] = 2.0E0*I_ERI_Pz_Pz_Pz_Px_C1001001000_a;
  abcd[660] = 2.0E0*I_ERI_Pz_Px_Px_Py_C1001001000_a;
  abcd[664] = 2.0E0*I_ERI_Pz_Py_Px_Py_C1001001000_a;
  abcd[668] = 2.0E0*I_ERI_Pz_Pz_Px_Py_C1001001000_a;
  abcd[676] = 2.0E0*I_ERI_Pz_Px_Py_Py_C1001001000_a;
  abcd[680] = 2.0E0*I_ERI_Pz_Py_Py_Py_C1001001000_a;
  abcd[684] = 2.0E0*I_ERI_Pz_Pz_Py_Py_C1001001000_a;
  abcd[692] = 2.0E0*I_ERI_Pz_Px_Pz_Py_C1001001000_a;
  abcd[696] = 2.0E0*I_ERI_Pz_Py_Pz_Py_C1001001000_a;
  abcd[700] = 2.0E0*I_ERI_Pz_Pz_Pz_Py_C1001001000_a;
  abcd[724] = 2.0E0*I_ERI_Pz_Px_Px_Pz_C1001001000_a;
  abcd[728] = 2.0E0*I_ERI_Pz_Py_Px_Pz_C1001001000_a;
  abcd[732] = 2.0E0*I_ERI_Pz_Pz_Px_Pz_C1001001000_a;
  abcd[740] = 2.0E0*I_ERI_Pz_Px_Py_Pz_C1001001000_a;
  abcd[744] = 2.0E0*I_ERI_Pz_Py_Py_Pz_C1001001000_a;
  abcd[748] = 2.0E0*I_ERI_Pz_Pz_Py_Pz_C1001001000_a;
  abcd[756] = 2.0E0*I_ERI_Pz_Px_Pz_Pz_C1001001000_a;
  abcd[760] = 2.0E0*I_ERI_Pz_Py_Pz_Pz_C1001001000_a;
  abcd[764] = 2.0E0*I_ERI_Pz_Pz_Pz_Pz_C1001001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_C1001001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001001001
   ************************************************************/
  abcd[597] = 2.0E0*I_ERI_Dxz_Px_Px_Px_C1001001001_a;
  abcd[598] = 2.0E0*I_ERI_Dyz_Px_Px_Px_C1001001001_a;
  abcd[599] = 2.0E0*I_ERI_D2z_Px_Px_Px_C1001001001_a-1*I_ERI_S_Px_Px_Px_C1001001001;
  abcd[601] = 2.0E0*I_ERI_Dxz_Py_Px_Px_C1001001001_a;
  abcd[602] = 2.0E0*I_ERI_Dyz_Py_Px_Px_C1001001001_a;
  abcd[603] = 2.0E0*I_ERI_D2z_Py_Px_Px_C1001001001_a-1*I_ERI_S_Py_Px_Px_C1001001001;
  abcd[605] = 2.0E0*I_ERI_Dxz_Pz_Px_Px_C1001001001_a;
  abcd[606] = 2.0E0*I_ERI_Dyz_Pz_Px_Px_C1001001001_a;
  abcd[607] = 2.0E0*I_ERI_D2z_Pz_Px_Px_C1001001001_a-1*I_ERI_S_Pz_Px_Px_C1001001001;
  abcd[613] = 2.0E0*I_ERI_Dxz_Px_Py_Px_C1001001001_a;
  abcd[614] = 2.0E0*I_ERI_Dyz_Px_Py_Px_C1001001001_a;
  abcd[615] = 2.0E0*I_ERI_D2z_Px_Py_Px_C1001001001_a-1*I_ERI_S_Px_Py_Px_C1001001001;
  abcd[617] = 2.0E0*I_ERI_Dxz_Py_Py_Px_C1001001001_a;
  abcd[618] = 2.0E0*I_ERI_Dyz_Py_Py_Px_C1001001001_a;
  abcd[619] = 2.0E0*I_ERI_D2z_Py_Py_Px_C1001001001_a-1*I_ERI_S_Py_Py_Px_C1001001001;
  abcd[621] = 2.0E0*I_ERI_Dxz_Pz_Py_Px_C1001001001_a;
  abcd[622] = 2.0E0*I_ERI_Dyz_Pz_Py_Px_C1001001001_a;
  abcd[623] = 2.0E0*I_ERI_D2z_Pz_Py_Px_C1001001001_a-1*I_ERI_S_Pz_Py_Px_C1001001001;
  abcd[629] = 2.0E0*I_ERI_Dxz_Px_Pz_Px_C1001001001_a;
  abcd[630] = 2.0E0*I_ERI_Dyz_Px_Pz_Px_C1001001001_a;
  abcd[631] = 2.0E0*I_ERI_D2z_Px_Pz_Px_C1001001001_a-1*I_ERI_S_Px_Pz_Px_C1001001001;
  abcd[633] = 2.0E0*I_ERI_Dxz_Py_Pz_Px_C1001001001_a;
  abcd[634] = 2.0E0*I_ERI_Dyz_Py_Pz_Px_C1001001001_a;
  abcd[635] = 2.0E0*I_ERI_D2z_Py_Pz_Px_C1001001001_a-1*I_ERI_S_Py_Pz_Px_C1001001001;
  abcd[637] = 2.0E0*I_ERI_Dxz_Pz_Pz_Px_C1001001001_a;
  abcd[638] = 2.0E0*I_ERI_Dyz_Pz_Pz_Px_C1001001001_a;
  abcd[639] = 2.0E0*I_ERI_D2z_Pz_Pz_Px_C1001001001_a-1*I_ERI_S_Pz_Pz_Px_C1001001001;
  abcd[661] = 2.0E0*I_ERI_Dxz_Px_Px_Py_C1001001001_a;
  abcd[662] = 2.0E0*I_ERI_Dyz_Px_Px_Py_C1001001001_a;
  abcd[663] = 2.0E0*I_ERI_D2z_Px_Px_Py_C1001001001_a-1*I_ERI_S_Px_Px_Py_C1001001001;
  abcd[665] = 2.0E0*I_ERI_Dxz_Py_Px_Py_C1001001001_a;
  abcd[666] = 2.0E0*I_ERI_Dyz_Py_Px_Py_C1001001001_a;
  abcd[667] = 2.0E0*I_ERI_D2z_Py_Px_Py_C1001001001_a-1*I_ERI_S_Py_Px_Py_C1001001001;
  abcd[669] = 2.0E0*I_ERI_Dxz_Pz_Px_Py_C1001001001_a;
  abcd[670] = 2.0E0*I_ERI_Dyz_Pz_Px_Py_C1001001001_a;
  abcd[671] = 2.0E0*I_ERI_D2z_Pz_Px_Py_C1001001001_a-1*I_ERI_S_Pz_Px_Py_C1001001001;
  abcd[677] = 2.0E0*I_ERI_Dxz_Px_Py_Py_C1001001001_a;
  abcd[678] = 2.0E0*I_ERI_Dyz_Px_Py_Py_C1001001001_a;
  abcd[679] = 2.0E0*I_ERI_D2z_Px_Py_Py_C1001001001_a-1*I_ERI_S_Px_Py_Py_C1001001001;
  abcd[681] = 2.0E0*I_ERI_Dxz_Py_Py_Py_C1001001001_a;
  abcd[682] = 2.0E0*I_ERI_Dyz_Py_Py_Py_C1001001001_a;
  abcd[683] = 2.0E0*I_ERI_D2z_Py_Py_Py_C1001001001_a-1*I_ERI_S_Py_Py_Py_C1001001001;
  abcd[685] = 2.0E0*I_ERI_Dxz_Pz_Py_Py_C1001001001_a;
  abcd[686] = 2.0E0*I_ERI_Dyz_Pz_Py_Py_C1001001001_a;
  abcd[687] = 2.0E0*I_ERI_D2z_Pz_Py_Py_C1001001001_a-1*I_ERI_S_Pz_Py_Py_C1001001001;
  abcd[693] = 2.0E0*I_ERI_Dxz_Px_Pz_Py_C1001001001_a;
  abcd[694] = 2.0E0*I_ERI_Dyz_Px_Pz_Py_C1001001001_a;
  abcd[695] = 2.0E0*I_ERI_D2z_Px_Pz_Py_C1001001001_a-1*I_ERI_S_Px_Pz_Py_C1001001001;
  abcd[697] = 2.0E0*I_ERI_Dxz_Py_Pz_Py_C1001001001_a;
  abcd[698] = 2.0E0*I_ERI_Dyz_Py_Pz_Py_C1001001001_a;
  abcd[699] = 2.0E0*I_ERI_D2z_Py_Pz_Py_C1001001001_a-1*I_ERI_S_Py_Pz_Py_C1001001001;
  abcd[701] = 2.0E0*I_ERI_Dxz_Pz_Pz_Py_C1001001001_a;
  abcd[702] = 2.0E0*I_ERI_Dyz_Pz_Pz_Py_C1001001001_a;
  abcd[703] = 2.0E0*I_ERI_D2z_Pz_Pz_Py_C1001001001_a-1*I_ERI_S_Pz_Pz_Py_C1001001001;
  abcd[725] = 2.0E0*I_ERI_Dxz_Px_Px_Pz_C1001001001_a;
  abcd[726] = 2.0E0*I_ERI_Dyz_Px_Px_Pz_C1001001001_a;
  abcd[727] = 2.0E0*I_ERI_D2z_Px_Px_Pz_C1001001001_a-1*I_ERI_S_Px_Px_Pz_C1001001001;
  abcd[729] = 2.0E0*I_ERI_Dxz_Py_Px_Pz_C1001001001_a;
  abcd[730] = 2.0E0*I_ERI_Dyz_Py_Px_Pz_C1001001001_a;
  abcd[731] = 2.0E0*I_ERI_D2z_Py_Px_Pz_C1001001001_a-1*I_ERI_S_Py_Px_Pz_C1001001001;
  abcd[733] = 2.0E0*I_ERI_Dxz_Pz_Px_Pz_C1001001001_a;
  abcd[734] = 2.0E0*I_ERI_Dyz_Pz_Px_Pz_C1001001001_a;
  abcd[735] = 2.0E0*I_ERI_D2z_Pz_Px_Pz_C1001001001_a-1*I_ERI_S_Pz_Px_Pz_C1001001001;
  abcd[741] = 2.0E0*I_ERI_Dxz_Px_Py_Pz_C1001001001_a;
  abcd[742] = 2.0E0*I_ERI_Dyz_Px_Py_Pz_C1001001001_a;
  abcd[743] = 2.0E0*I_ERI_D2z_Px_Py_Pz_C1001001001_a-1*I_ERI_S_Px_Py_Pz_C1001001001;
  abcd[745] = 2.0E0*I_ERI_Dxz_Py_Py_Pz_C1001001001_a;
  abcd[746] = 2.0E0*I_ERI_Dyz_Py_Py_Pz_C1001001001_a;
  abcd[747] = 2.0E0*I_ERI_D2z_Py_Py_Pz_C1001001001_a-1*I_ERI_S_Py_Py_Pz_C1001001001;
  abcd[749] = 2.0E0*I_ERI_Dxz_Pz_Py_Pz_C1001001001_a;
  abcd[750] = 2.0E0*I_ERI_Dyz_Pz_Py_Pz_C1001001001_a;
  abcd[751] = 2.0E0*I_ERI_D2z_Pz_Py_Pz_C1001001001_a-1*I_ERI_S_Pz_Py_Pz_C1001001001;
  abcd[757] = 2.0E0*I_ERI_Dxz_Px_Pz_Pz_C1001001001_a;
  abcd[758] = 2.0E0*I_ERI_Dyz_Px_Pz_Pz_C1001001001_a;
  abcd[759] = 2.0E0*I_ERI_D2z_Px_Pz_Pz_C1001001001_a-1*I_ERI_S_Px_Pz_Pz_C1001001001;
  abcd[761] = 2.0E0*I_ERI_Dxz_Py_Pz_Pz_C1001001001_a;
  abcd[762] = 2.0E0*I_ERI_Dyz_Py_Pz_Pz_C1001001001_a;
  abcd[763] = 2.0E0*I_ERI_D2z_Py_Pz_Pz_C1001001001_a-1*I_ERI_S_Py_Pz_Pz_C1001001001;
  abcd[765] = 2.0E0*I_ERI_Dxz_Pz_Pz_Pz_C1001001001_a;
  abcd[766] = 2.0E0*I_ERI_Dyz_Pz_Pz_Pz_C1001001001_a;
  abcd[767] = 2.0E0*I_ERI_D2z_Pz_Pz_Pz_C1001001001_a-1*I_ERI_S_Pz_Pz_Pz_C1001001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C0_b
   ************************************************************/
  abcd[768] = 2.0E0*I_ERI_S_Px_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1_b
   ************************************************************/
  abcd[769] = 2.0E0*I_ERI_Px_Px_S_S_C1_b;
  abcd[770] = 2.0E0*I_ERI_Py_Px_S_S_C1_b;
  abcd[771] = 2.0E0*I_ERI_Pz_Px_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000
   ************************************************************/
  abcd[772] = 2.0E0*I_ERI_S_D2x_S_S_C1000_b-1*I_ERI_S_S_S_S_C1000;
  abcd[776] = 2.0E0*I_ERI_S_Dxy_S_S_C1000_b;
  abcd[780] = 2.0E0*I_ERI_S_Dxz_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001
   ************************************************************/
  abcd[773] = 2.0E0*I_ERI_Px_D2x_S_S_C1001_b-1*I_ERI_Px_S_S_S_C1001;
  abcd[774] = 2.0E0*I_ERI_Py_D2x_S_S_C1001_b-1*I_ERI_Py_S_S_S_C1001;
  abcd[775] = 2.0E0*I_ERI_Pz_D2x_S_S_C1001_b-1*I_ERI_Pz_S_S_S_C1001;
  abcd[777] = 2.0E0*I_ERI_Px_Dxy_S_S_C1001_b;
  abcd[778] = 2.0E0*I_ERI_Py_Dxy_S_S_C1001_b;
  abcd[779] = 2.0E0*I_ERI_Pz_Dxy_S_S_C1001_b;
  abcd[781] = 2.0E0*I_ERI_Px_Dxz_S_S_C1001_b;
  abcd[782] = 2.0E0*I_ERI_Py_Dxz_S_S_C1001_b;
  abcd[783] = 2.0E0*I_ERI_Pz_Dxz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000000_b
   ************************************************************/
  abcd[784] = 2.0E0*I_ERI_S_Px_Px_S_C1000000_b;
  abcd[800] = 2.0E0*I_ERI_S_Px_Py_S_C1000000_b;
  abcd[816] = 2.0E0*I_ERI_S_Px_Pz_S_C1000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[785] = 2.0E0*I_ERI_Px_Px_Px_S_C1000001_b;
  abcd[786] = 2.0E0*I_ERI_Py_Px_Px_S_C1000001_b;
  abcd[787] = 2.0E0*I_ERI_Pz_Px_Px_S_C1000001_b;
  abcd[801] = 2.0E0*I_ERI_Px_Px_Py_S_C1000001_b;
  abcd[802] = 2.0E0*I_ERI_Py_Px_Py_S_C1000001_b;
  abcd[803] = 2.0E0*I_ERI_Pz_Px_Py_S_C1000001_b;
  abcd[817] = 2.0E0*I_ERI_Px_Px_Pz_S_C1000001_b;
  abcd[818] = 2.0E0*I_ERI_Py_Px_Pz_S_C1000001_b;
  abcd[819] = 2.0E0*I_ERI_Pz_Px_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1001000
   ************************************************************/
  abcd[788] = 2.0E0*I_ERI_S_D2x_Px_S_C1001000_b-1*I_ERI_S_S_Px_S_C1001000;
  abcd[792] = 2.0E0*I_ERI_S_Dxy_Px_S_C1001000_b;
  abcd[796] = 2.0E0*I_ERI_S_Dxz_Px_S_C1001000_b;
  abcd[804] = 2.0E0*I_ERI_S_D2x_Py_S_C1001000_b-1*I_ERI_S_S_Py_S_C1001000;
  abcd[808] = 2.0E0*I_ERI_S_Dxy_Py_S_C1001000_b;
  abcd[812] = 2.0E0*I_ERI_S_Dxz_Py_S_C1001000_b;
  abcd[820] = 2.0E0*I_ERI_S_D2x_Pz_S_C1001000_b-1*I_ERI_S_S_Pz_S_C1001000;
  abcd[824] = 2.0E0*I_ERI_S_Dxy_Pz_S_C1001000_b;
  abcd[828] = 2.0E0*I_ERI_S_Dxz_Pz_S_C1001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[789] = 2.0E0*I_ERI_Px_D2x_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[790] = 2.0E0*I_ERI_Py_D2x_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[791] = 2.0E0*I_ERI_Pz_D2x_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[793] = 2.0E0*I_ERI_Px_Dxy_Px_S_C1001001_b;
  abcd[794] = 2.0E0*I_ERI_Py_Dxy_Px_S_C1001001_b;
  abcd[795] = 2.0E0*I_ERI_Pz_Dxy_Px_S_C1001001_b;
  abcd[797] = 2.0E0*I_ERI_Px_Dxz_Px_S_C1001001_b;
  abcd[798] = 2.0E0*I_ERI_Py_Dxz_Px_S_C1001001_b;
  abcd[799] = 2.0E0*I_ERI_Pz_Dxz_Px_S_C1001001_b;
  abcd[805] = 2.0E0*I_ERI_Px_D2x_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[806] = 2.0E0*I_ERI_Py_D2x_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[807] = 2.0E0*I_ERI_Pz_D2x_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[809] = 2.0E0*I_ERI_Px_Dxy_Py_S_C1001001_b;
  abcd[810] = 2.0E0*I_ERI_Py_Dxy_Py_S_C1001001_b;
  abcd[811] = 2.0E0*I_ERI_Pz_Dxy_Py_S_C1001001_b;
  abcd[813] = 2.0E0*I_ERI_Px_Dxz_Py_S_C1001001_b;
  abcd[814] = 2.0E0*I_ERI_Py_Dxz_Py_S_C1001001_b;
  abcd[815] = 2.0E0*I_ERI_Pz_Dxz_Py_S_C1001001_b;
  abcd[821] = 2.0E0*I_ERI_Px_D2x_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[822] = 2.0E0*I_ERI_Py_D2x_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[823] = 2.0E0*I_ERI_Pz_D2x_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;
  abcd[825] = 2.0E0*I_ERI_Px_Dxy_Pz_S_C1001001_b;
  abcd[826] = 2.0E0*I_ERI_Py_Dxy_Pz_S_C1001001_b;
  abcd[827] = 2.0E0*I_ERI_Pz_Dxy_Pz_S_C1001001_b;
  abcd[829] = 2.0E0*I_ERI_Px_Dxz_Pz_S_C1001001_b;
  abcd[830] = 2.0E0*I_ERI_Py_Dxz_Pz_S_C1001001_b;
  abcd[831] = 2.0E0*I_ERI_Pz_Dxz_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000000000_b
   ************************************************************/
  abcd[832] = 2.0E0*I_ERI_S_Px_S_Px_C1000000000_b;
  abcd[896] = 2.0E0*I_ERI_S_Px_S_Py_C1000000000_b;
  abcd[960] = 2.0E0*I_ERI_S_Px_S_Pz_C1000000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000000001_b
   ************************************************************/
  abcd[833] = 2.0E0*I_ERI_Px_Px_S_Px_C1000000001_b;
  abcd[834] = 2.0E0*I_ERI_Py_Px_S_Px_C1000000001_b;
  abcd[835] = 2.0E0*I_ERI_Pz_Px_S_Px_C1000000001_b;
  abcd[897] = 2.0E0*I_ERI_Px_Px_S_Py_C1000000001_b;
  abcd[898] = 2.0E0*I_ERI_Py_Px_S_Py_C1000000001_b;
  abcd[899] = 2.0E0*I_ERI_Pz_Px_S_Py_C1000000001_b;
  abcd[961] = 2.0E0*I_ERI_Px_Px_S_Pz_C1000000001_b;
  abcd[962] = 2.0E0*I_ERI_Py_Px_S_Pz_C1000000001_b;
  abcd[963] = 2.0E0*I_ERI_Pz_Px_S_Pz_C1000000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_P_C1000001000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000001000
   ************************************************************/
  abcd[836] = 2.0E0*I_ERI_S_D2x_S_Px_C1000001000_b-1*I_ERI_S_S_S_Px_C1000001000;
  abcd[840] = 2.0E0*I_ERI_S_Dxy_S_Px_C1000001000_b;
  abcd[844] = 2.0E0*I_ERI_S_Dxz_S_Px_C1000001000_b;
  abcd[900] = 2.0E0*I_ERI_S_D2x_S_Py_C1000001000_b-1*I_ERI_S_S_S_Py_C1000001000;
  abcd[904] = 2.0E0*I_ERI_S_Dxy_S_Py_C1000001000_b;
  abcd[908] = 2.0E0*I_ERI_S_Dxz_S_Py_C1000001000_b;
  abcd[964] = 2.0E0*I_ERI_S_D2x_S_Pz_C1000001000_b-1*I_ERI_S_S_S_Pz_C1000001000;
  abcd[968] = 2.0E0*I_ERI_S_Dxy_S_Pz_C1000001000_b;
  abcd[972] = 2.0E0*I_ERI_S_Dxz_S_Pz_C1000001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001001
   ************************************************************/
  abcd[837] = 2.0E0*I_ERI_Px_D2x_S_Px_C1000001001_b-1*I_ERI_Px_S_S_Px_C1000001001;
  abcd[838] = 2.0E0*I_ERI_Py_D2x_S_Px_C1000001001_b-1*I_ERI_Py_S_S_Px_C1000001001;
  abcd[839] = 2.0E0*I_ERI_Pz_D2x_S_Px_C1000001001_b-1*I_ERI_Pz_S_S_Px_C1000001001;
  abcd[841] = 2.0E0*I_ERI_Px_Dxy_S_Px_C1000001001_b;
  abcd[842] = 2.0E0*I_ERI_Py_Dxy_S_Px_C1000001001_b;
  abcd[843] = 2.0E0*I_ERI_Pz_Dxy_S_Px_C1000001001_b;
  abcd[845] = 2.0E0*I_ERI_Px_Dxz_S_Px_C1000001001_b;
  abcd[846] = 2.0E0*I_ERI_Py_Dxz_S_Px_C1000001001_b;
  abcd[847] = 2.0E0*I_ERI_Pz_Dxz_S_Px_C1000001001_b;
  abcd[901] = 2.0E0*I_ERI_Px_D2x_S_Py_C1000001001_b-1*I_ERI_Px_S_S_Py_C1000001001;
  abcd[902] = 2.0E0*I_ERI_Py_D2x_S_Py_C1000001001_b-1*I_ERI_Py_S_S_Py_C1000001001;
  abcd[903] = 2.0E0*I_ERI_Pz_D2x_S_Py_C1000001001_b-1*I_ERI_Pz_S_S_Py_C1000001001;
  abcd[905] = 2.0E0*I_ERI_Px_Dxy_S_Py_C1000001001_b;
  abcd[906] = 2.0E0*I_ERI_Py_Dxy_S_Py_C1000001001_b;
  abcd[907] = 2.0E0*I_ERI_Pz_Dxy_S_Py_C1000001001_b;
  abcd[909] = 2.0E0*I_ERI_Px_Dxz_S_Py_C1000001001_b;
  abcd[910] = 2.0E0*I_ERI_Py_Dxz_S_Py_C1000001001_b;
  abcd[911] = 2.0E0*I_ERI_Pz_Dxz_S_Py_C1000001001_b;
  abcd[965] = 2.0E0*I_ERI_Px_D2x_S_Pz_C1000001001_b-1*I_ERI_Px_S_S_Pz_C1000001001;
  abcd[966] = 2.0E0*I_ERI_Py_D2x_S_Pz_C1000001001_b-1*I_ERI_Py_S_S_Pz_C1000001001;
  abcd[967] = 2.0E0*I_ERI_Pz_D2x_S_Pz_C1000001001_b-1*I_ERI_Pz_S_S_Pz_C1000001001;
  abcd[969] = 2.0E0*I_ERI_Px_Dxy_S_Pz_C1000001001_b;
  abcd[970] = 2.0E0*I_ERI_Py_Dxy_S_Pz_C1000001001_b;
  abcd[971] = 2.0E0*I_ERI_Pz_Dxy_S_Pz_C1000001001_b;
  abcd[973] = 2.0E0*I_ERI_Px_Dxz_S_Pz_C1000001001_b;
  abcd[974] = 2.0E0*I_ERI_Py_Dxz_S_Pz_C1000001001_b;
  abcd[975] = 2.0E0*I_ERI_Pz_Dxz_S_Pz_C1000001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000000_b
   ************************************************************/
  abcd[848] = 2.0E0*I_ERI_S_Px_Px_Px_C1001000000_b;
  abcd[864] = 2.0E0*I_ERI_S_Px_Py_Px_C1001000000_b;
  abcd[880] = 2.0E0*I_ERI_S_Px_Pz_Px_C1001000000_b;
  abcd[912] = 2.0E0*I_ERI_S_Px_Px_Py_C1001000000_b;
  abcd[928] = 2.0E0*I_ERI_S_Px_Py_Py_C1001000000_b;
  abcd[944] = 2.0E0*I_ERI_S_Px_Pz_Py_C1001000000_b;
  abcd[976] = 2.0E0*I_ERI_S_Px_Px_Pz_C1001000000_b;
  abcd[992] = 2.0E0*I_ERI_S_Px_Py_Pz_C1001000000_b;
  abcd[1008] = 2.0E0*I_ERI_S_Px_Pz_Pz_C1001000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001000001_b
   ************************************************************/
  abcd[849] = 2.0E0*I_ERI_Px_Px_Px_Px_C1001000001_b;
  abcd[850] = 2.0E0*I_ERI_Py_Px_Px_Px_C1001000001_b;
  abcd[851] = 2.0E0*I_ERI_Pz_Px_Px_Px_C1001000001_b;
  abcd[865] = 2.0E0*I_ERI_Px_Px_Py_Px_C1001000001_b;
  abcd[866] = 2.0E0*I_ERI_Py_Px_Py_Px_C1001000001_b;
  abcd[867] = 2.0E0*I_ERI_Pz_Px_Py_Px_C1001000001_b;
  abcd[881] = 2.0E0*I_ERI_Px_Px_Pz_Px_C1001000001_b;
  abcd[882] = 2.0E0*I_ERI_Py_Px_Pz_Px_C1001000001_b;
  abcd[883] = 2.0E0*I_ERI_Pz_Px_Pz_Px_C1001000001_b;
  abcd[913] = 2.0E0*I_ERI_Px_Px_Px_Py_C1001000001_b;
  abcd[914] = 2.0E0*I_ERI_Py_Px_Px_Py_C1001000001_b;
  abcd[915] = 2.0E0*I_ERI_Pz_Px_Px_Py_C1001000001_b;
  abcd[929] = 2.0E0*I_ERI_Px_Px_Py_Py_C1001000001_b;
  abcd[930] = 2.0E0*I_ERI_Py_Px_Py_Py_C1001000001_b;
  abcd[931] = 2.0E0*I_ERI_Pz_Px_Py_Py_C1001000001_b;
  abcd[945] = 2.0E0*I_ERI_Px_Px_Pz_Py_C1001000001_b;
  abcd[946] = 2.0E0*I_ERI_Py_Px_Pz_Py_C1001000001_b;
  abcd[947] = 2.0E0*I_ERI_Pz_Px_Pz_Py_C1001000001_b;
  abcd[977] = 2.0E0*I_ERI_Px_Px_Px_Pz_C1001000001_b;
  abcd[978] = 2.0E0*I_ERI_Py_Px_Px_Pz_C1001000001_b;
  abcd[979] = 2.0E0*I_ERI_Pz_Px_Px_Pz_C1001000001_b;
  abcd[993] = 2.0E0*I_ERI_Px_Px_Py_Pz_C1001000001_b;
  abcd[994] = 2.0E0*I_ERI_Py_Px_Py_Pz_C1001000001_b;
  abcd[995] = 2.0E0*I_ERI_Pz_Px_Py_Pz_C1001000001_b;
  abcd[1009] = 2.0E0*I_ERI_Px_Px_Pz_Pz_C1001000001_b;
  abcd[1010] = 2.0E0*I_ERI_Py_Px_Pz_Pz_C1001000001_b;
  abcd[1011] = 2.0E0*I_ERI_Pz_Px_Pz_Pz_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_P_C1001001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001001000
   ************************************************************/
  abcd[852] = 2.0E0*I_ERI_S_D2x_Px_Px_C1001001000_b-1*I_ERI_S_S_Px_Px_C1001001000;
  abcd[856] = 2.0E0*I_ERI_S_Dxy_Px_Px_C1001001000_b;
  abcd[860] = 2.0E0*I_ERI_S_Dxz_Px_Px_C1001001000_b;
  abcd[868] = 2.0E0*I_ERI_S_D2x_Py_Px_C1001001000_b-1*I_ERI_S_S_Py_Px_C1001001000;
  abcd[872] = 2.0E0*I_ERI_S_Dxy_Py_Px_C1001001000_b;
  abcd[876] = 2.0E0*I_ERI_S_Dxz_Py_Px_C1001001000_b;
  abcd[884] = 2.0E0*I_ERI_S_D2x_Pz_Px_C1001001000_b-1*I_ERI_S_S_Pz_Px_C1001001000;
  abcd[888] = 2.0E0*I_ERI_S_Dxy_Pz_Px_C1001001000_b;
  abcd[892] = 2.0E0*I_ERI_S_Dxz_Pz_Px_C1001001000_b;
  abcd[916] = 2.0E0*I_ERI_S_D2x_Px_Py_C1001001000_b-1*I_ERI_S_S_Px_Py_C1001001000;
  abcd[920] = 2.0E0*I_ERI_S_Dxy_Px_Py_C1001001000_b;
  abcd[924] = 2.0E0*I_ERI_S_Dxz_Px_Py_C1001001000_b;
  abcd[932] = 2.0E0*I_ERI_S_D2x_Py_Py_C1001001000_b-1*I_ERI_S_S_Py_Py_C1001001000;
  abcd[936] = 2.0E0*I_ERI_S_Dxy_Py_Py_C1001001000_b;
  abcd[940] = 2.0E0*I_ERI_S_Dxz_Py_Py_C1001001000_b;
  abcd[948] = 2.0E0*I_ERI_S_D2x_Pz_Py_C1001001000_b-1*I_ERI_S_S_Pz_Py_C1001001000;
  abcd[952] = 2.0E0*I_ERI_S_Dxy_Pz_Py_C1001001000_b;
  abcd[956] = 2.0E0*I_ERI_S_Dxz_Pz_Py_C1001001000_b;
  abcd[980] = 2.0E0*I_ERI_S_D2x_Px_Pz_C1001001000_b-1*I_ERI_S_S_Px_Pz_C1001001000;
  abcd[984] = 2.0E0*I_ERI_S_Dxy_Px_Pz_C1001001000_b;
  abcd[988] = 2.0E0*I_ERI_S_Dxz_Px_Pz_C1001001000_b;
  abcd[996] = 2.0E0*I_ERI_S_D2x_Py_Pz_C1001001000_b-1*I_ERI_S_S_Py_Pz_C1001001000;
  abcd[1000] = 2.0E0*I_ERI_S_Dxy_Py_Pz_C1001001000_b;
  abcd[1004] = 2.0E0*I_ERI_S_Dxz_Py_Pz_C1001001000_b;
  abcd[1012] = 2.0E0*I_ERI_S_D2x_Pz_Pz_C1001001000_b-1*I_ERI_S_S_Pz_Pz_C1001001000;
  abcd[1016] = 2.0E0*I_ERI_S_Dxy_Pz_Pz_C1001001000_b;
  abcd[1020] = 2.0E0*I_ERI_S_Dxz_Pz_Pz_C1001001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001001
   ************************************************************/
  abcd[853] = 2.0E0*I_ERI_Px_D2x_Px_Px_C1001001001_b-1*I_ERI_Px_S_Px_Px_C1001001001;
  abcd[854] = 2.0E0*I_ERI_Py_D2x_Px_Px_C1001001001_b-1*I_ERI_Py_S_Px_Px_C1001001001;
  abcd[855] = 2.0E0*I_ERI_Pz_D2x_Px_Px_C1001001001_b-1*I_ERI_Pz_S_Px_Px_C1001001001;
  abcd[857] = 2.0E0*I_ERI_Px_Dxy_Px_Px_C1001001001_b;
  abcd[858] = 2.0E0*I_ERI_Py_Dxy_Px_Px_C1001001001_b;
  abcd[859] = 2.0E0*I_ERI_Pz_Dxy_Px_Px_C1001001001_b;
  abcd[861] = 2.0E0*I_ERI_Px_Dxz_Px_Px_C1001001001_b;
  abcd[862] = 2.0E0*I_ERI_Py_Dxz_Px_Px_C1001001001_b;
  abcd[863] = 2.0E0*I_ERI_Pz_Dxz_Px_Px_C1001001001_b;
  abcd[869] = 2.0E0*I_ERI_Px_D2x_Py_Px_C1001001001_b-1*I_ERI_Px_S_Py_Px_C1001001001;
  abcd[870] = 2.0E0*I_ERI_Py_D2x_Py_Px_C1001001001_b-1*I_ERI_Py_S_Py_Px_C1001001001;
  abcd[871] = 2.0E0*I_ERI_Pz_D2x_Py_Px_C1001001001_b-1*I_ERI_Pz_S_Py_Px_C1001001001;
  abcd[873] = 2.0E0*I_ERI_Px_Dxy_Py_Px_C1001001001_b;
  abcd[874] = 2.0E0*I_ERI_Py_Dxy_Py_Px_C1001001001_b;
  abcd[875] = 2.0E0*I_ERI_Pz_Dxy_Py_Px_C1001001001_b;
  abcd[877] = 2.0E0*I_ERI_Px_Dxz_Py_Px_C1001001001_b;
  abcd[878] = 2.0E0*I_ERI_Py_Dxz_Py_Px_C1001001001_b;
  abcd[879] = 2.0E0*I_ERI_Pz_Dxz_Py_Px_C1001001001_b;
  abcd[885] = 2.0E0*I_ERI_Px_D2x_Pz_Px_C1001001001_b-1*I_ERI_Px_S_Pz_Px_C1001001001;
  abcd[886] = 2.0E0*I_ERI_Py_D2x_Pz_Px_C1001001001_b-1*I_ERI_Py_S_Pz_Px_C1001001001;
  abcd[887] = 2.0E0*I_ERI_Pz_D2x_Pz_Px_C1001001001_b-1*I_ERI_Pz_S_Pz_Px_C1001001001;
  abcd[889] = 2.0E0*I_ERI_Px_Dxy_Pz_Px_C1001001001_b;
  abcd[890] = 2.0E0*I_ERI_Py_Dxy_Pz_Px_C1001001001_b;
  abcd[891] = 2.0E0*I_ERI_Pz_Dxy_Pz_Px_C1001001001_b;
  abcd[893] = 2.0E0*I_ERI_Px_Dxz_Pz_Px_C1001001001_b;
  abcd[894] = 2.0E0*I_ERI_Py_Dxz_Pz_Px_C1001001001_b;
  abcd[895] = 2.0E0*I_ERI_Pz_Dxz_Pz_Px_C1001001001_b;
  abcd[917] = 2.0E0*I_ERI_Px_D2x_Px_Py_C1001001001_b-1*I_ERI_Px_S_Px_Py_C1001001001;
  abcd[918] = 2.0E0*I_ERI_Py_D2x_Px_Py_C1001001001_b-1*I_ERI_Py_S_Px_Py_C1001001001;
  abcd[919] = 2.0E0*I_ERI_Pz_D2x_Px_Py_C1001001001_b-1*I_ERI_Pz_S_Px_Py_C1001001001;
  abcd[921] = 2.0E0*I_ERI_Px_Dxy_Px_Py_C1001001001_b;
  abcd[922] = 2.0E0*I_ERI_Py_Dxy_Px_Py_C1001001001_b;
  abcd[923] = 2.0E0*I_ERI_Pz_Dxy_Px_Py_C1001001001_b;
  abcd[925] = 2.0E0*I_ERI_Px_Dxz_Px_Py_C1001001001_b;
  abcd[926] = 2.0E0*I_ERI_Py_Dxz_Px_Py_C1001001001_b;
  abcd[927] = 2.0E0*I_ERI_Pz_Dxz_Px_Py_C1001001001_b;
  abcd[933] = 2.0E0*I_ERI_Px_D2x_Py_Py_C1001001001_b-1*I_ERI_Px_S_Py_Py_C1001001001;
  abcd[934] = 2.0E0*I_ERI_Py_D2x_Py_Py_C1001001001_b-1*I_ERI_Py_S_Py_Py_C1001001001;
  abcd[935] = 2.0E0*I_ERI_Pz_D2x_Py_Py_C1001001001_b-1*I_ERI_Pz_S_Py_Py_C1001001001;
  abcd[937] = 2.0E0*I_ERI_Px_Dxy_Py_Py_C1001001001_b;
  abcd[938] = 2.0E0*I_ERI_Py_Dxy_Py_Py_C1001001001_b;
  abcd[939] = 2.0E0*I_ERI_Pz_Dxy_Py_Py_C1001001001_b;
  abcd[941] = 2.0E0*I_ERI_Px_Dxz_Py_Py_C1001001001_b;
  abcd[942] = 2.0E0*I_ERI_Py_Dxz_Py_Py_C1001001001_b;
  abcd[943] = 2.0E0*I_ERI_Pz_Dxz_Py_Py_C1001001001_b;
  abcd[949] = 2.0E0*I_ERI_Px_D2x_Pz_Py_C1001001001_b-1*I_ERI_Px_S_Pz_Py_C1001001001;
  abcd[950] = 2.0E0*I_ERI_Py_D2x_Pz_Py_C1001001001_b-1*I_ERI_Py_S_Pz_Py_C1001001001;
  abcd[951] = 2.0E0*I_ERI_Pz_D2x_Pz_Py_C1001001001_b-1*I_ERI_Pz_S_Pz_Py_C1001001001;
  abcd[953] = 2.0E0*I_ERI_Px_Dxy_Pz_Py_C1001001001_b;
  abcd[954] = 2.0E0*I_ERI_Py_Dxy_Pz_Py_C1001001001_b;
  abcd[955] = 2.0E0*I_ERI_Pz_Dxy_Pz_Py_C1001001001_b;
  abcd[957] = 2.0E0*I_ERI_Px_Dxz_Pz_Py_C1001001001_b;
  abcd[958] = 2.0E0*I_ERI_Py_Dxz_Pz_Py_C1001001001_b;
  abcd[959] = 2.0E0*I_ERI_Pz_Dxz_Pz_Py_C1001001001_b;
  abcd[981] = 2.0E0*I_ERI_Px_D2x_Px_Pz_C1001001001_b-1*I_ERI_Px_S_Px_Pz_C1001001001;
  abcd[982] = 2.0E0*I_ERI_Py_D2x_Px_Pz_C1001001001_b-1*I_ERI_Py_S_Px_Pz_C1001001001;
  abcd[983] = 2.0E0*I_ERI_Pz_D2x_Px_Pz_C1001001001_b-1*I_ERI_Pz_S_Px_Pz_C1001001001;
  abcd[985] = 2.0E0*I_ERI_Px_Dxy_Px_Pz_C1001001001_b;
  abcd[986] = 2.0E0*I_ERI_Py_Dxy_Px_Pz_C1001001001_b;
  abcd[987] = 2.0E0*I_ERI_Pz_Dxy_Px_Pz_C1001001001_b;
  abcd[989] = 2.0E0*I_ERI_Px_Dxz_Px_Pz_C1001001001_b;
  abcd[990] = 2.0E0*I_ERI_Py_Dxz_Px_Pz_C1001001001_b;
  abcd[991] = 2.0E0*I_ERI_Pz_Dxz_Px_Pz_C1001001001_b;
  abcd[997] = 2.0E0*I_ERI_Px_D2x_Py_Pz_C1001001001_b-1*I_ERI_Px_S_Py_Pz_C1001001001;
  abcd[998] = 2.0E0*I_ERI_Py_D2x_Py_Pz_C1001001001_b-1*I_ERI_Py_S_Py_Pz_C1001001001;
  abcd[999] = 2.0E0*I_ERI_Pz_D2x_Py_Pz_C1001001001_b-1*I_ERI_Pz_S_Py_Pz_C1001001001;
  abcd[1001] = 2.0E0*I_ERI_Px_Dxy_Py_Pz_C1001001001_b;
  abcd[1002] = 2.0E0*I_ERI_Py_Dxy_Py_Pz_C1001001001_b;
  abcd[1003] = 2.0E0*I_ERI_Pz_Dxy_Py_Pz_C1001001001_b;
  abcd[1005] = 2.0E0*I_ERI_Px_Dxz_Py_Pz_C1001001001_b;
  abcd[1006] = 2.0E0*I_ERI_Py_Dxz_Py_Pz_C1001001001_b;
  abcd[1007] = 2.0E0*I_ERI_Pz_Dxz_Py_Pz_C1001001001_b;
  abcd[1013] = 2.0E0*I_ERI_Px_D2x_Pz_Pz_C1001001001_b-1*I_ERI_Px_S_Pz_Pz_C1001001001;
  abcd[1014] = 2.0E0*I_ERI_Py_D2x_Pz_Pz_C1001001001_b-1*I_ERI_Py_S_Pz_Pz_C1001001001;
  abcd[1015] = 2.0E0*I_ERI_Pz_D2x_Pz_Pz_C1001001001_b-1*I_ERI_Pz_S_Pz_Pz_C1001001001;
  abcd[1017] = 2.0E0*I_ERI_Px_Dxy_Pz_Pz_C1001001001_b;
  abcd[1018] = 2.0E0*I_ERI_Py_Dxy_Pz_Pz_C1001001001_b;
  abcd[1019] = 2.0E0*I_ERI_Pz_Dxy_Pz_Pz_C1001001001_b;
  abcd[1021] = 2.0E0*I_ERI_Px_Dxz_Pz_Pz_C1001001001_b;
  abcd[1022] = 2.0E0*I_ERI_Py_Dxz_Pz_Pz_C1001001001_b;
  abcd[1023] = 2.0E0*I_ERI_Pz_Dxz_Pz_Pz_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C0_b
   ************************************************************/
  abcd[1024] = 2.0E0*I_ERI_S_Py_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1_b
   ************************************************************/
  abcd[1025] = 2.0E0*I_ERI_Px_Py_S_S_C1_b;
  abcd[1026] = 2.0E0*I_ERI_Py_Py_S_S_C1_b;
  abcd[1027] = 2.0E0*I_ERI_Pz_Py_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000
   ************************************************************/
  abcd[1028] = 2.0E0*I_ERI_S_Dxy_S_S_C1000_b;
  abcd[1032] = 2.0E0*I_ERI_S_D2y_S_S_C1000_b-1*I_ERI_S_S_S_S_C1000;
  abcd[1036] = 2.0E0*I_ERI_S_Dyz_S_S_C1000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001
   ************************************************************/
  abcd[1029] = 2.0E0*I_ERI_Px_Dxy_S_S_C1001_b;
  abcd[1030] = 2.0E0*I_ERI_Py_Dxy_S_S_C1001_b;
  abcd[1031] = 2.0E0*I_ERI_Pz_Dxy_S_S_C1001_b;
  abcd[1033] = 2.0E0*I_ERI_Px_D2y_S_S_C1001_b-1*I_ERI_Px_S_S_S_C1001;
  abcd[1034] = 2.0E0*I_ERI_Py_D2y_S_S_C1001_b-1*I_ERI_Py_S_S_S_C1001;
  abcd[1035] = 2.0E0*I_ERI_Pz_D2y_S_S_C1001_b-1*I_ERI_Pz_S_S_S_C1001;
  abcd[1037] = 2.0E0*I_ERI_Px_Dyz_S_S_C1001_b;
  abcd[1038] = 2.0E0*I_ERI_Py_Dyz_S_S_C1001_b;
  abcd[1039] = 2.0E0*I_ERI_Pz_Dyz_S_S_C1001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000000_b
   ************************************************************/
  abcd[1040] = 2.0E0*I_ERI_S_Py_Px_S_C1000000_b;
  abcd[1056] = 2.0E0*I_ERI_S_Py_Py_S_C1000000_b;
  abcd[1072] = 2.0E0*I_ERI_S_Py_Pz_S_C1000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[1041] = 2.0E0*I_ERI_Px_Py_Px_S_C1000001_b;
  abcd[1042] = 2.0E0*I_ERI_Py_Py_Px_S_C1000001_b;
  abcd[1043] = 2.0E0*I_ERI_Pz_Py_Px_S_C1000001_b;
  abcd[1057] = 2.0E0*I_ERI_Px_Py_Py_S_C1000001_b;
  abcd[1058] = 2.0E0*I_ERI_Py_Py_Py_S_C1000001_b;
  abcd[1059] = 2.0E0*I_ERI_Pz_Py_Py_S_C1000001_b;
  abcd[1073] = 2.0E0*I_ERI_Px_Py_Pz_S_C1000001_b;
  abcd[1074] = 2.0E0*I_ERI_Py_Py_Pz_S_C1000001_b;
  abcd[1075] = 2.0E0*I_ERI_Pz_Py_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1001000
   ************************************************************/
  abcd[1044] = 2.0E0*I_ERI_S_Dxy_Px_S_C1001000_b;
  abcd[1048] = 2.0E0*I_ERI_S_D2y_Px_S_C1001000_b-1*I_ERI_S_S_Px_S_C1001000;
  abcd[1052] = 2.0E0*I_ERI_S_Dyz_Px_S_C1001000_b;
  abcd[1060] = 2.0E0*I_ERI_S_Dxy_Py_S_C1001000_b;
  abcd[1064] = 2.0E0*I_ERI_S_D2y_Py_S_C1001000_b-1*I_ERI_S_S_Py_S_C1001000;
  abcd[1068] = 2.0E0*I_ERI_S_Dyz_Py_S_C1001000_b;
  abcd[1076] = 2.0E0*I_ERI_S_Dxy_Pz_S_C1001000_b;
  abcd[1080] = 2.0E0*I_ERI_S_D2y_Pz_S_C1001000_b-1*I_ERI_S_S_Pz_S_C1001000;
  abcd[1084] = 2.0E0*I_ERI_S_Dyz_Pz_S_C1001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[1045] = 2.0E0*I_ERI_Px_Dxy_Px_S_C1001001_b;
  abcd[1046] = 2.0E0*I_ERI_Py_Dxy_Px_S_C1001001_b;
  abcd[1047] = 2.0E0*I_ERI_Pz_Dxy_Px_S_C1001001_b;
  abcd[1049] = 2.0E0*I_ERI_Px_D2y_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[1050] = 2.0E0*I_ERI_Py_D2y_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[1051] = 2.0E0*I_ERI_Pz_D2y_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[1053] = 2.0E0*I_ERI_Px_Dyz_Px_S_C1001001_b;
  abcd[1054] = 2.0E0*I_ERI_Py_Dyz_Px_S_C1001001_b;
  abcd[1055] = 2.0E0*I_ERI_Pz_Dyz_Px_S_C1001001_b;
  abcd[1061] = 2.0E0*I_ERI_Px_Dxy_Py_S_C1001001_b;
  abcd[1062] = 2.0E0*I_ERI_Py_Dxy_Py_S_C1001001_b;
  abcd[1063] = 2.0E0*I_ERI_Pz_Dxy_Py_S_C1001001_b;
  abcd[1065] = 2.0E0*I_ERI_Px_D2y_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[1066] = 2.0E0*I_ERI_Py_D2y_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[1067] = 2.0E0*I_ERI_Pz_D2y_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[1069] = 2.0E0*I_ERI_Px_Dyz_Py_S_C1001001_b;
  abcd[1070] = 2.0E0*I_ERI_Py_Dyz_Py_S_C1001001_b;
  abcd[1071] = 2.0E0*I_ERI_Pz_Dyz_Py_S_C1001001_b;
  abcd[1077] = 2.0E0*I_ERI_Px_Dxy_Pz_S_C1001001_b;
  abcd[1078] = 2.0E0*I_ERI_Py_Dxy_Pz_S_C1001001_b;
  abcd[1079] = 2.0E0*I_ERI_Pz_Dxy_Pz_S_C1001001_b;
  abcd[1081] = 2.0E0*I_ERI_Px_D2y_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[1082] = 2.0E0*I_ERI_Py_D2y_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[1083] = 2.0E0*I_ERI_Pz_D2y_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;
  abcd[1085] = 2.0E0*I_ERI_Px_Dyz_Pz_S_C1001001_b;
  abcd[1086] = 2.0E0*I_ERI_Py_Dyz_Pz_S_C1001001_b;
  abcd[1087] = 2.0E0*I_ERI_Pz_Dyz_Pz_S_C1001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000000000_b
   ************************************************************/
  abcd[1088] = 2.0E0*I_ERI_S_Py_S_Px_C1000000000_b;
  abcd[1152] = 2.0E0*I_ERI_S_Py_S_Py_C1000000000_b;
  abcd[1216] = 2.0E0*I_ERI_S_Py_S_Pz_C1000000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000000001_b
   ************************************************************/
  abcd[1089] = 2.0E0*I_ERI_Px_Py_S_Px_C1000000001_b;
  abcd[1090] = 2.0E0*I_ERI_Py_Py_S_Px_C1000000001_b;
  abcd[1091] = 2.0E0*I_ERI_Pz_Py_S_Px_C1000000001_b;
  abcd[1153] = 2.0E0*I_ERI_Px_Py_S_Py_C1000000001_b;
  abcd[1154] = 2.0E0*I_ERI_Py_Py_S_Py_C1000000001_b;
  abcd[1155] = 2.0E0*I_ERI_Pz_Py_S_Py_C1000000001_b;
  abcd[1217] = 2.0E0*I_ERI_Px_Py_S_Pz_C1000000001_b;
  abcd[1218] = 2.0E0*I_ERI_Py_Py_S_Pz_C1000000001_b;
  abcd[1219] = 2.0E0*I_ERI_Pz_Py_S_Pz_C1000000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_P_C1000001000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000001000
   ************************************************************/
  abcd[1092] = 2.0E0*I_ERI_S_Dxy_S_Px_C1000001000_b;
  abcd[1096] = 2.0E0*I_ERI_S_D2y_S_Px_C1000001000_b-1*I_ERI_S_S_S_Px_C1000001000;
  abcd[1100] = 2.0E0*I_ERI_S_Dyz_S_Px_C1000001000_b;
  abcd[1156] = 2.0E0*I_ERI_S_Dxy_S_Py_C1000001000_b;
  abcd[1160] = 2.0E0*I_ERI_S_D2y_S_Py_C1000001000_b-1*I_ERI_S_S_S_Py_C1000001000;
  abcd[1164] = 2.0E0*I_ERI_S_Dyz_S_Py_C1000001000_b;
  abcd[1220] = 2.0E0*I_ERI_S_Dxy_S_Pz_C1000001000_b;
  abcd[1224] = 2.0E0*I_ERI_S_D2y_S_Pz_C1000001000_b-1*I_ERI_S_S_S_Pz_C1000001000;
  abcd[1228] = 2.0E0*I_ERI_S_Dyz_S_Pz_C1000001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001001
   ************************************************************/
  abcd[1093] = 2.0E0*I_ERI_Px_Dxy_S_Px_C1000001001_b;
  abcd[1094] = 2.0E0*I_ERI_Py_Dxy_S_Px_C1000001001_b;
  abcd[1095] = 2.0E0*I_ERI_Pz_Dxy_S_Px_C1000001001_b;
  abcd[1097] = 2.0E0*I_ERI_Px_D2y_S_Px_C1000001001_b-1*I_ERI_Px_S_S_Px_C1000001001;
  abcd[1098] = 2.0E0*I_ERI_Py_D2y_S_Px_C1000001001_b-1*I_ERI_Py_S_S_Px_C1000001001;
  abcd[1099] = 2.0E0*I_ERI_Pz_D2y_S_Px_C1000001001_b-1*I_ERI_Pz_S_S_Px_C1000001001;
  abcd[1101] = 2.0E0*I_ERI_Px_Dyz_S_Px_C1000001001_b;
  abcd[1102] = 2.0E0*I_ERI_Py_Dyz_S_Px_C1000001001_b;
  abcd[1103] = 2.0E0*I_ERI_Pz_Dyz_S_Px_C1000001001_b;
  abcd[1157] = 2.0E0*I_ERI_Px_Dxy_S_Py_C1000001001_b;
  abcd[1158] = 2.0E0*I_ERI_Py_Dxy_S_Py_C1000001001_b;
  abcd[1159] = 2.0E0*I_ERI_Pz_Dxy_S_Py_C1000001001_b;
  abcd[1161] = 2.0E0*I_ERI_Px_D2y_S_Py_C1000001001_b-1*I_ERI_Px_S_S_Py_C1000001001;
  abcd[1162] = 2.0E0*I_ERI_Py_D2y_S_Py_C1000001001_b-1*I_ERI_Py_S_S_Py_C1000001001;
  abcd[1163] = 2.0E0*I_ERI_Pz_D2y_S_Py_C1000001001_b-1*I_ERI_Pz_S_S_Py_C1000001001;
  abcd[1165] = 2.0E0*I_ERI_Px_Dyz_S_Py_C1000001001_b;
  abcd[1166] = 2.0E0*I_ERI_Py_Dyz_S_Py_C1000001001_b;
  abcd[1167] = 2.0E0*I_ERI_Pz_Dyz_S_Py_C1000001001_b;
  abcd[1221] = 2.0E0*I_ERI_Px_Dxy_S_Pz_C1000001001_b;
  abcd[1222] = 2.0E0*I_ERI_Py_Dxy_S_Pz_C1000001001_b;
  abcd[1223] = 2.0E0*I_ERI_Pz_Dxy_S_Pz_C1000001001_b;
  abcd[1225] = 2.0E0*I_ERI_Px_D2y_S_Pz_C1000001001_b-1*I_ERI_Px_S_S_Pz_C1000001001;
  abcd[1226] = 2.0E0*I_ERI_Py_D2y_S_Pz_C1000001001_b-1*I_ERI_Py_S_S_Pz_C1000001001;
  abcd[1227] = 2.0E0*I_ERI_Pz_D2y_S_Pz_C1000001001_b-1*I_ERI_Pz_S_S_Pz_C1000001001;
  abcd[1229] = 2.0E0*I_ERI_Px_Dyz_S_Pz_C1000001001_b;
  abcd[1230] = 2.0E0*I_ERI_Py_Dyz_S_Pz_C1000001001_b;
  abcd[1231] = 2.0E0*I_ERI_Pz_Dyz_S_Pz_C1000001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000000_b
   ************************************************************/
  abcd[1104] = 2.0E0*I_ERI_S_Py_Px_Px_C1001000000_b;
  abcd[1120] = 2.0E0*I_ERI_S_Py_Py_Px_C1001000000_b;
  abcd[1136] = 2.0E0*I_ERI_S_Py_Pz_Px_C1001000000_b;
  abcd[1168] = 2.0E0*I_ERI_S_Py_Px_Py_C1001000000_b;
  abcd[1184] = 2.0E0*I_ERI_S_Py_Py_Py_C1001000000_b;
  abcd[1200] = 2.0E0*I_ERI_S_Py_Pz_Py_C1001000000_b;
  abcd[1232] = 2.0E0*I_ERI_S_Py_Px_Pz_C1001000000_b;
  abcd[1248] = 2.0E0*I_ERI_S_Py_Py_Pz_C1001000000_b;
  abcd[1264] = 2.0E0*I_ERI_S_Py_Pz_Pz_C1001000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001000001_b
   ************************************************************/
  abcd[1105] = 2.0E0*I_ERI_Px_Py_Px_Px_C1001000001_b;
  abcd[1106] = 2.0E0*I_ERI_Py_Py_Px_Px_C1001000001_b;
  abcd[1107] = 2.0E0*I_ERI_Pz_Py_Px_Px_C1001000001_b;
  abcd[1121] = 2.0E0*I_ERI_Px_Py_Py_Px_C1001000001_b;
  abcd[1122] = 2.0E0*I_ERI_Py_Py_Py_Px_C1001000001_b;
  abcd[1123] = 2.0E0*I_ERI_Pz_Py_Py_Px_C1001000001_b;
  abcd[1137] = 2.0E0*I_ERI_Px_Py_Pz_Px_C1001000001_b;
  abcd[1138] = 2.0E0*I_ERI_Py_Py_Pz_Px_C1001000001_b;
  abcd[1139] = 2.0E0*I_ERI_Pz_Py_Pz_Px_C1001000001_b;
  abcd[1169] = 2.0E0*I_ERI_Px_Py_Px_Py_C1001000001_b;
  abcd[1170] = 2.0E0*I_ERI_Py_Py_Px_Py_C1001000001_b;
  abcd[1171] = 2.0E0*I_ERI_Pz_Py_Px_Py_C1001000001_b;
  abcd[1185] = 2.0E0*I_ERI_Px_Py_Py_Py_C1001000001_b;
  abcd[1186] = 2.0E0*I_ERI_Py_Py_Py_Py_C1001000001_b;
  abcd[1187] = 2.0E0*I_ERI_Pz_Py_Py_Py_C1001000001_b;
  abcd[1201] = 2.0E0*I_ERI_Px_Py_Pz_Py_C1001000001_b;
  abcd[1202] = 2.0E0*I_ERI_Py_Py_Pz_Py_C1001000001_b;
  abcd[1203] = 2.0E0*I_ERI_Pz_Py_Pz_Py_C1001000001_b;
  abcd[1233] = 2.0E0*I_ERI_Px_Py_Px_Pz_C1001000001_b;
  abcd[1234] = 2.0E0*I_ERI_Py_Py_Px_Pz_C1001000001_b;
  abcd[1235] = 2.0E0*I_ERI_Pz_Py_Px_Pz_C1001000001_b;
  abcd[1249] = 2.0E0*I_ERI_Px_Py_Py_Pz_C1001000001_b;
  abcd[1250] = 2.0E0*I_ERI_Py_Py_Py_Pz_C1001000001_b;
  abcd[1251] = 2.0E0*I_ERI_Pz_Py_Py_Pz_C1001000001_b;
  abcd[1265] = 2.0E0*I_ERI_Px_Py_Pz_Pz_C1001000001_b;
  abcd[1266] = 2.0E0*I_ERI_Py_Py_Pz_Pz_C1001000001_b;
  abcd[1267] = 2.0E0*I_ERI_Pz_Py_Pz_Pz_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_P_C1001001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001001000
   ************************************************************/
  abcd[1108] = 2.0E0*I_ERI_S_Dxy_Px_Px_C1001001000_b;
  abcd[1112] = 2.0E0*I_ERI_S_D2y_Px_Px_C1001001000_b-1*I_ERI_S_S_Px_Px_C1001001000;
  abcd[1116] = 2.0E0*I_ERI_S_Dyz_Px_Px_C1001001000_b;
  abcd[1124] = 2.0E0*I_ERI_S_Dxy_Py_Px_C1001001000_b;
  abcd[1128] = 2.0E0*I_ERI_S_D2y_Py_Px_C1001001000_b-1*I_ERI_S_S_Py_Px_C1001001000;
  abcd[1132] = 2.0E0*I_ERI_S_Dyz_Py_Px_C1001001000_b;
  abcd[1140] = 2.0E0*I_ERI_S_Dxy_Pz_Px_C1001001000_b;
  abcd[1144] = 2.0E0*I_ERI_S_D2y_Pz_Px_C1001001000_b-1*I_ERI_S_S_Pz_Px_C1001001000;
  abcd[1148] = 2.0E0*I_ERI_S_Dyz_Pz_Px_C1001001000_b;
  abcd[1172] = 2.0E0*I_ERI_S_Dxy_Px_Py_C1001001000_b;
  abcd[1176] = 2.0E0*I_ERI_S_D2y_Px_Py_C1001001000_b-1*I_ERI_S_S_Px_Py_C1001001000;
  abcd[1180] = 2.0E0*I_ERI_S_Dyz_Px_Py_C1001001000_b;
  abcd[1188] = 2.0E0*I_ERI_S_Dxy_Py_Py_C1001001000_b;
  abcd[1192] = 2.0E0*I_ERI_S_D2y_Py_Py_C1001001000_b-1*I_ERI_S_S_Py_Py_C1001001000;
  abcd[1196] = 2.0E0*I_ERI_S_Dyz_Py_Py_C1001001000_b;
  abcd[1204] = 2.0E0*I_ERI_S_Dxy_Pz_Py_C1001001000_b;
  abcd[1208] = 2.0E0*I_ERI_S_D2y_Pz_Py_C1001001000_b-1*I_ERI_S_S_Pz_Py_C1001001000;
  abcd[1212] = 2.0E0*I_ERI_S_Dyz_Pz_Py_C1001001000_b;
  abcd[1236] = 2.0E0*I_ERI_S_Dxy_Px_Pz_C1001001000_b;
  abcd[1240] = 2.0E0*I_ERI_S_D2y_Px_Pz_C1001001000_b-1*I_ERI_S_S_Px_Pz_C1001001000;
  abcd[1244] = 2.0E0*I_ERI_S_Dyz_Px_Pz_C1001001000_b;
  abcd[1252] = 2.0E0*I_ERI_S_Dxy_Py_Pz_C1001001000_b;
  abcd[1256] = 2.0E0*I_ERI_S_D2y_Py_Pz_C1001001000_b-1*I_ERI_S_S_Py_Pz_C1001001000;
  abcd[1260] = 2.0E0*I_ERI_S_Dyz_Py_Pz_C1001001000_b;
  abcd[1268] = 2.0E0*I_ERI_S_Dxy_Pz_Pz_C1001001000_b;
  abcd[1272] = 2.0E0*I_ERI_S_D2y_Pz_Pz_C1001001000_b-1*I_ERI_S_S_Pz_Pz_C1001001000;
  abcd[1276] = 2.0E0*I_ERI_S_Dyz_Pz_Pz_C1001001000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001001
   ************************************************************/
  abcd[1109] = 2.0E0*I_ERI_Px_Dxy_Px_Px_C1001001001_b;
  abcd[1110] = 2.0E0*I_ERI_Py_Dxy_Px_Px_C1001001001_b;
  abcd[1111] = 2.0E0*I_ERI_Pz_Dxy_Px_Px_C1001001001_b;
  abcd[1113] = 2.0E0*I_ERI_Px_D2y_Px_Px_C1001001001_b-1*I_ERI_Px_S_Px_Px_C1001001001;
  abcd[1114] = 2.0E0*I_ERI_Py_D2y_Px_Px_C1001001001_b-1*I_ERI_Py_S_Px_Px_C1001001001;
  abcd[1115] = 2.0E0*I_ERI_Pz_D2y_Px_Px_C1001001001_b-1*I_ERI_Pz_S_Px_Px_C1001001001;
  abcd[1117] = 2.0E0*I_ERI_Px_Dyz_Px_Px_C1001001001_b;
  abcd[1118] = 2.0E0*I_ERI_Py_Dyz_Px_Px_C1001001001_b;
  abcd[1119] = 2.0E0*I_ERI_Pz_Dyz_Px_Px_C1001001001_b;
  abcd[1125] = 2.0E0*I_ERI_Px_Dxy_Py_Px_C1001001001_b;
  abcd[1126] = 2.0E0*I_ERI_Py_Dxy_Py_Px_C1001001001_b;
  abcd[1127] = 2.0E0*I_ERI_Pz_Dxy_Py_Px_C1001001001_b;
  abcd[1129] = 2.0E0*I_ERI_Px_D2y_Py_Px_C1001001001_b-1*I_ERI_Px_S_Py_Px_C1001001001;
  abcd[1130] = 2.0E0*I_ERI_Py_D2y_Py_Px_C1001001001_b-1*I_ERI_Py_S_Py_Px_C1001001001;
  abcd[1131] = 2.0E0*I_ERI_Pz_D2y_Py_Px_C1001001001_b-1*I_ERI_Pz_S_Py_Px_C1001001001;
  abcd[1133] = 2.0E0*I_ERI_Px_Dyz_Py_Px_C1001001001_b;
  abcd[1134] = 2.0E0*I_ERI_Py_Dyz_Py_Px_C1001001001_b;
  abcd[1135] = 2.0E0*I_ERI_Pz_Dyz_Py_Px_C1001001001_b;
  abcd[1141] = 2.0E0*I_ERI_Px_Dxy_Pz_Px_C1001001001_b;
  abcd[1142] = 2.0E0*I_ERI_Py_Dxy_Pz_Px_C1001001001_b;
  abcd[1143] = 2.0E0*I_ERI_Pz_Dxy_Pz_Px_C1001001001_b;
  abcd[1145] = 2.0E0*I_ERI_Px_D2y_Pz_Px_C1001001001_b-1*I_ERI_Px_S_Pz_Px_C1001001001;
  abcd[1146] = 2.0E0*I_ERI_Py_D2y_Pz_Px_C1001001001_b-1*I_ERI_Py_S_Pz_Px_C1001001001;
  abcd[1147] = 2.0E0*I_ERI_Pz_D2y_Pz_Px_C1001001001_b-1*I_ERI_Pz_S_Pz_Px_C1001001001;
  abcd[1149] = 2.0E0*I_ERI_Px_Dyz_Pz_Px_C1001001001_b;
  abcd[1150] = 2.0E0*I_ERI_Py_Dyz_Pz_Px_C1001001001_b;
  abcd[1151] = 2.0E0*I_ERI_Pz_Dyz_Pz_Px_C1001001001_b;
  abcd[1173] = 2.0E0*I_ERI_Px_Dxy_Px_Py_C1001001001_b;
  abcd[1174] = 2.0E0*I_ERI_Py_Dxy_Px_Py_C1001001001_b;
  abcd[1175] = 2.0E0*I_ERI_Pz_Dxy_Px_Py_C1001001001_b;
  abcd[1177] = 2.0E0*I_ERI_Px_D2y_Px_Py_C1001001001_b-1*I_ERI_Px_S_Px_Py_C1001001001;
  abcd[1178] = 2.0E0*I_ERI_Py_D2y_Px_Py_C1001001001_b-1*I_ERI_Py_S_Px_Py_C1001001001;
  abcd[1179] = 2.0E0*I_ERI_Pz_D2y_Px_Py_C1001001001_b-1*I_ERI_Pz_S_Px_Py_C1001001001;
  abcd[1181] = 2.0E0*I_ERI_Px_Dyz_Px_Py_C1001001001_b;
  abcd[1182] = 2.0E0*I_ERI_Py_Dyz_Px_Py_C1001001001_b;
  abcd[1183] = 2.0E0*I_ERI_Pz_Dyz_Px_Py_C1001001001_b;
  abcd[1189] = 2.0E0*I_ERI_Px_Dxy_Py_Py_C1001001001_b;
  abcd[1190] = 2.0E0*I_ERI_Py_Dxy_Py_Py_C1001001001_b;
  abcd[1191] = 2.0E0*I_ERI_Pz_Dxy_Py_Py_C1001001001_b;
  abcd[1193] = 2.0E0*I_ERI_Px_D2y_Py_Py_C1001001001_b-1*I_ERI_Px_S_Py_Py_C1001001001;
  abcd[1194] = 2.0E0*I_ERI_Py_D2y_Py_Py_C1001001001_b-1*I_ERI_Py_S_Py_Py_C1001001001;
  abcd[1195] = 2.0E0*I_ERI_Pz_D2y_Py_Py_C1001001001_b-1*I_ERI_Pz_S_Py_Py_C1001001001;
  abcd[1197] = 2.0E0*I_ERI_Px_Dyz_Py_Py_C1001001001_b;
  abcd[1198] = 2.0E0*I_ERI_Py_Dyz_Py_Py_C1001001001_b;
  abcd[1199] = 2.0E0*I_ERI_Pz_Dyz_Py_Py_C1001001001_b;
  abcd[1205] = 2.0E0*I_ERI_Px_Dxy_Pz_Py_C1001001001_b;
  abcd[1206] = 2.0E0*I_ERI_Py_Dxy_Pz_Py_C1001001001_b;
  abcd[1207] = 2.0E0*I_ERI_Pz_Dxy_Pz_Py_C1001001001_b;
  abcd[1209] = 2.0E0*I_ERI_Px_D2y_Pz_Py_C1001001001_b-1*I_ERI_Px_S_Pz_Py_C1001001001;
  abcd[1210] = 2.0E0*I_ERI_Py_D2y_Pz_Py_C1001001001_b-1*I_ERI_Py_S_Pz_Py_C1001001001;
  abcd[1211] = 2.0E0*I_ERI_Pz_D2y_Pz_Py_C1001001001_b-1*I_ERI_Pz_S_Pz_Py_C1001001001;
  abcd[1213] = 2.0E0*I_ERI_Px_Dyz_Pz_Py_C1001001001_b;
  abcd[1214] = 2.0E0*I_ERI_Py_Dyz_Pz_Py_C1001001001_b;
  abcd[1215] = 2.0E0*I_ERI_Pz_Dyz_Pz_Py_C1001001001_b;
  abcd[1237] = 2.0E0*I_ERI_Px_Dxy_Px_Pz_C1001001001_b;
  abcd[1238] = 2.0E0*I_ERI_Py_Dxy_Px_Pz_C1001001001_b;
  abcd[1239] = 2.0E0*I_ERI_Pz_Dxy_Px_Pz_C1001001001_b;
  abcd[1241] = 2.0E0*I_ERI_Px_D2y_Px_Pz_C1001001001_b-1*I_ERI_Px_S_Px_Pz_C1001001001;
  abcd[1242] = 2.0E0*I_ERI_Py_D2y_Px_Pz_C1001001001_b-1*I_ERI_Py_S_Px_Pz_C1001001001;
  abcd[1243] = 2.0E0*I_ERI_Pz_D2y_Px_Pz_C1001001001_b-1*I_ERI_Pz_S_Px_Pz_C1001001001;
  abcd[1245] = 2.0E0*I_ERI_Px_Dyz_Px_Pz_C1001001001_b;
  abcd[1246] = 2.0E0*I_ERI_Py_Dyz_Px_Pz_C1001001001_b;
  abcd[1247] = 2.0E0*I_ERI_Pz_Dyz_Px_Pz_C1001001001_b;
  abcd[1253] = 2.0E0*I_ERI_Px_Dxy_Py_Pz_C1001001001_b;
  abcd[1254] = 2.0E0*I_ERI_Py_Dxy_Py_Pz_C1001001001_b;
  abcd[1255] = 2.0E0*I_ERI_Pz_Dxy_Py_Pz_C1001001001_b;
  abcd[1257] = 2.0E0*I_ERI_Px_D2y_Py_Pz_C1001001001_b-1*I_ERI_Px_S_Py_Pz_C1001001001;
  abcd[1258] = 2.0E0*I_ERI_Py_D2y_Py_Pz_C1001001001_b-1*I_ERI_Py_S_Py_Pz_C1001001001;
  abcd[1259] = 2.0E0*I_ERI_Pz_D2y_Py_Pz_C1001001001_b-1*I_ERI_Pz_S_Py_Pz_C1001001001;
  abcd[1261] = 2.0E0*I_ERI_Px_Dyz_Py_Pz_C1001001001_b;
  abcd[1262] = 2.0E0*I_ERI_Py_Dyz_Py_Pz_C1001001001_b;
  abcd[1263] = 2.0E0*I_ERI_Pz_Dyz_Py_Pz_C1001001001_b;
  abcd[1269] = 2.0E0*I_ERI_Px_Dxy_Pz_Pz_C1001001001_b;
  abcd[1270] = 2.0E0*I_ERI_Py_Dxy_Pz_Pz_C1001001001_b;
  abcd[1271] = 2.0E0*I_ERI_Pz_Dxy_Pz_Pz_C1001001001_b;
  abcd[1273] = 2.0E0*I_ERI_Px_D2y_Pz_Pz_C1001001001_b-1*I_ERI_Px_S_Pz_Pz_C1001001001;
  abcd[1274] = 2.0E0*I_ERI_Py_D2y_Pz_Pz_C1001001001_b-1*I_ERI_Py_S_Pz_Pz_C1001001001;
  abcd[1275] = 2.0E0*I_ERI_Pz_D2y_Pz_Pz_C1001001001_b-1*I_ERI_Pz_S_Pz_Pz_C1001001001;
  abcd[1277] = 2.0E0*I_ERI_Px_Dyz_Pz_Pz_C1001001001_b;
  abcd[1278] = 2.0E0*I_ERI_Py_Dyz_Pz_Pz_C1001001001_b;
  abcd[1279] = 2.0E0*I_ERI_Pz_Dyz_Pz_Pz_C1001001001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C0_b
   ************************************************************/
  abcd[1280] = 2.0E0*I_ERI_S_Pz_S_S_C0_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1_b
   ************************************************************/
  abcd[1281] = 2.0E0*I_ERI_Px_Pz_S_S_C1_b;
  abcd[1282] = 2.0E0*I_ERI_Py_Pz_S_S_C1_b;
  abcd[1283] = 2.0E0*I_ERI_Pz_Pz_S_S_C1_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_S_C1000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000
   ************************************************************/
  abcd[1284] = 2.0E0*I_ERI_S_Dxz_S_S_C1000_b;
  abcd[1288] = 2.0E0*I_ERI_S_Dyz_S_S_C1000_b;
  abcd[1292] = 2.0E0*I_ERI_S_D2z_S_S_C1000_b-1*I_ERI_S_S_S_S_C1000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_S_C1001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1001
   ************************************************************/
  abcd[1285] = 2.0E0*I_ERI_Px_Dxz_S_S_C1001_b;
  abcd[1286] = 2.0E0*I_ERI_Py_Dxz_S_S_C1001_b;
  abcd[1287] = 2.0E0*I_ERI_Pz_Dxz_S_S_C1001_b;
  abcd[1289] = 2.0E0*I_ERI_Px_Dyz_S_S_C1001_b;
  abcd[1290] = 2.0E0*I_ERI_Py_Dyz_S_S_C1001_b;
  abcd[1291] = 2.0E0*I_ERI_Pz_Dyz_S_S_C1001_b;
  abcd[1293] = 2.0E0*I_ERI_Px_D2z_S_S_C1001_b-1*I_ERI_Px_S_S_S_C1001;
  abcd[1294] = 2.0E0*I_ERI_Py_D2z_S_S_C1001_b-1*I_ERI_Py_S_S_S_C1001;
  abcd[1295] = 2.0E0*I_ERI_Pz_D2z_S_S_C1001_b-1*I_ERI_Pz_S_S_S_C1001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000000_b
   ************************************************************/
  abcd[1296] = 2.0E0*I_ERI_S_Pz_Px_S_C1000000_b;
  abcd[1312] = 2.0E0*I_ERI_S_Pz_Py_S_C1000000_b;
  abcd[1328] = 2.0E0*I_ERI_S_Pz_Pz_S_C1000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1000001_b
   ************************************************************/
  abcd[1297] = 2.0E0*I_ERI_Px_Pz_Px_S_C1000001_b;
  abcd[1298] = 2.0E0*I_ERI_Py_Pz_Px_S_C1000001_b;
  abcd[1299] = 2.0E0*I_ERI_Pz_Pz_Px_S_C1000001_b;
  abcd[1313] = 2.0E0*I_ERI_Px_Pz_Py_S_C1000001_b;
  abcd[1314] = 2.0E0*I_ERI_Py_Pz_Py_S_C1000001_b;
  abcd[1315] = 2.0E0*I_ERI_Pz_Pz_Py_S_C1000001_b;
  abcd[1329] = 2.0E0*I_ERI_Px_Pz_Pz_S_C1000001_b;
  abcd[1330] = 2.0E0*I_ERI_Py_Pz_Pz_S_C1000001_b;
  abcd[1331] = 2.0E0*I_ERI_Pz_Pz_Pz_S_C1000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_S_C1001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1001000
   ************************************************************/
  abcd[1300] = 2.0E0*I_ERI_S_Dxz_Px_S_C1001000_b;
  abcd[1304] = 2.0E0*I_ERI_S_Dyz_Px_S_C1001000_b;
  abcd[1308] = 2.0E0*I_ERI_S_D2z_Px_S_C1001000_b-1*I_ERI_S_S_Px_S_C1001000;
  abcd[1316] = 2.0E0*I_ERI_S_Dxz_Py_S_C1001000_b;
  abcd[1320] = 2.0E0*I_ERI_S_Dyz_Py_S_C1001000_b;
  abcd[1324] = 2.0E0*I_ERI_S_D2z_Py_S_C1001000_b-1*I_ERI_S_S_Py_S_C1001000;
  abcd[1332] = 2.0E0*I_ERI_S_Dxz_Pz_S_C1001000_b;
  abcd[1336] = 2.0E0*I_ERI_S_Dyz_Pz_S_C1001000_b;
  abcd[1340] = 2.0E0*I_ERI_S_D2z_Pz_S_C1001000_b-1*I_ERI_S_S_Pz_S_C1001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_S_C1001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001
   ************************************************************/
  abcd[1301] = 2.0E0*I_ERI_Px_Dxz_Px_S_C1001001_b;
  abcd[1302] = 2.0E0*I_ERI_Py_Dxz_Px_S_C1001001_b;
  abcd[1303] = 2.0E0*I_ERI_Pz_Dxz_Px_S_C1001001_b;
  abcd[1305] = 2.0E0*I_ERI_Px_Dyz_Px_S_C1001001_b;
  abcd[1306] = 2.0E0*I_ERI_Py_Dyz_Px_S_C1001001_b;
  abcd[1307] = 2.0E0*I_ERI_Pz_Dyz_Px_S_C1001001_b;
  abcd[1309] = 2.0E0*I_ERI_Px_D2z_Px_S_C1001001_b-1*I_ERI_Px_S_Px_S_C1001001;
  abcd[1310] = 2.0E0*I_ERI_Py_D2z_Px_S_C1001001_b-1*I_ERI_Py_S_Px_S_C1001001;
  abcd[1311] = 2.0E0*I_ERI_Pz_D2z_Px_S_C1001001_b-1*I_ERI_Pz_S_Px_S_C1001001;
  abcd[1317] = 2.0E0*I_ERI_Px_Dxz_Py_S_C1001001_b;
  abcd[1318] = 2.0E0*I_ERI_Py_Dxz_Py_S_C1001001_b;
  abcd[1319] = 2.0E0*I_ERI_Pz_Dxz_Py_S_C1001001_b;
  abcd[1321] = 2.0E0*I_ERI_Px_Dyz_Py_S_C1001001_b;
  abcd[1322] = 2.0E0*I_ERI_Py_Dyz_Py_S_C1001001_b;
  abcd[1323] = 2.0E0*I_ERI_Pz_Dyz_Py_S_C1001001_b;
  abcd[1325] = 2.0E0*I_ERI_Px_D2z_Py_S_C1001001_b-1*I_ERI_Px_S_Py_S_C1001001;
  abcd[1326] = 2.0E0*I_ERI_Py_D2z_Py_S_C1001001_b-1*I_ERI_Py_S_Py_S_C1001001;
  abcd[1327] = 2.0E0*I_ERI_Pz_D2z_Py_S_C1001001_b-1*I_ERI_Pz_S_Py_S_C1001001;
  abcd[1333] = 2.0E0*I_ERI_Px_Dxz_Pz_S_C1001001_b;
  abcd[1334] = 2.0E0*I_ERI_Py_Dxz_Pz_S_C1001001_b;
  abcd[1335] = 2.0E0*I_ERI_Pz_Dxz_Pz_S_C1001001_b;
  abcd[1337] = 2.0E0*I_ERI_Px_Dyz_Pz_S_C1001001_b;
  abcd[1338] = 2.0E0*I_ERI_Py_Dyz_Pz_S_C1001001_b;
  abcd[1339] = 2.0E0*I_ERI_Pz_Dyz_Pz_S_C1001001_b;
  abcd[1341] = 2.0E0*I_ERI_Px_D2z_Pz_S_C1001001_b-1*I_ERI_Px_S_Pz_S_C1001001;
  abcd[1342] = 2.0E0*I_ERI_Py_D2z_Pz_S_C1001001_b-1*I_ERI_Py_S_Pz_S_C1001001;
  abcd[1343] = 2.0E0*I_ERI_Pz_D2z_Pz_S_C1001001_b-1*I_ERI_Pz_S_Pz_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1000000000_b
   ************************************************************/
  abcd[1344] = 2.0E0*I_ERI_S_Pz_S_Px_C1000000000_b;
  abcd[1408] = 2.0E0*I_ERI_S_Pz_S_Py_C1000000000_b;
  abcd[1472] = 2.0E0*I_ERI_S_Pz_S_Pz_C1000000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1000000001_b
   ************************************************************/
  abcd[1345] = 2.0E0*I_ERI_Px_Pz_S_Px_C1000000001_b;
  abcd[1346] = 2.0E0*I_ERI_Py_Pz_S_Px_C1000000001_b;
  abcd[1347] = 2.0E0*I_ERI_Pz_Pz_S_Px_C1000000001_b;
  abcd[1409] = 2.0E0*I_ERI_Px_Pz_S_Py_C1000000001_b;
  abcd[1410] = 2.0E0*I_ERI_Py_Pz_S_Py_C1000000001_b;
  abcd[1411] = 2.0E0*I_ERI_Pz_Pz_S_Py_C1000000001_b;
  abcd[1473] = 2.0E0*I_ERI_Px_Pz_S_Pz_C1000000001_b;
  abcd[1474] = 2.0E0*I_ERI_Py_Pz_S_Pz_C1000000001_b;
  abcd[1475] = 2.0E0*I_ERI_Pz_Pz_S_Pz_C1000000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_S_P_C1000001000_b
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1000001000
   ************************************************************/
  abcd[1348] = 2.0E0*I_ERI_S_Dxz_S_Px_C1000001000_b;
  abcd[1352] = 2.0E0*I_ERI_S_Dyz_S_Px_C1000001000_b;
  abcd[1356] = 2.0E0*I_ERI_S_D2z_S_Px_C1000001000_b-1*I_ERI_S_S_S_Px_C1000001000;
  abcd[1412] = 2.0E0*I_ERI_S_Dxz_S_Py_C1000001000_b;
  abcd[1416] = 2.0E0*I_ERI_S_Dyz_S_Py_C1000001000_b;
  abcd[1420] = 2.0E0*I_ERI_S_D2z_S_Py_C1000001000_b-1*I_ERI_S_S_S_Py_C1000001000;
  abcd[1476] = 2.0E0*I_ERI_S_Dxz_S_Pz_C1000001000_b;
  abcd[1480] = 2.0E0*I_ERI_S_Dyz_S_Pz_C1000001000_b;
  abcd[1484] = 2.0E0*I_ERI_S_D2z_S_Pz_C1000001000_b-1*I_ERI_S_S_S_Pz_C1000001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_S_P_C1000001001_b
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1000001001
   ************************************************************/
  abcd[1349] = 2.0E0*I_ERI_Px_Dxz_S_Px_C1000001001_b;
  abcd[1350] = 2.0E0*I_ERI_Py_Dxz_S_Px_C1000001001_b;
  abcd[1351] = 2.0E0*I_ERI_Pz_Dxz_S_Px_C1000001001_b;
  abcd[1353] = 2.0E0*I_ERI_Px_Dyz_S_Px_C1000001001_b;
  abcd[1354] = 2.0E0*I_ERI_Py_Dyz_S_Px_C1000001001_b;
  abcd[1355] = 2.0E0*I_ERI_Pz_Dyz_S_Px_C1000001001_b;
  abcd[1357] = 2.0E0*I_ERI_Px_D2z_S_Px_C1000001001_b-1*I_ERI_Px_S_S_Px_C1000001001;
  abcd[1358] = 2.0E0*I_ERI_Py_D2z_S_Px_C1000001001_b-1*I_ERI_Py_S_S_Px_C1000001001;
  abcd[1359] = 2.0E0*I_ERI_Pz_D2z_S_Px_C1000001001_b-1*I_ERI_Pz_S_S_Px_C1000001001;
  abcd[1413] = 2.0E0*I_ERI_Px_Dxz_S_Py_C1000001001_b;
  abcd[1414] = 2.0E0*I_ERI_Py_Dxz_S_Py_C1000001001_b;
  abcd[1415] = 2.0E0*I_ERI_Pz_Dxz_S_Py_C1000001001_b;
  abcd[1417] = 2.0E0*I_ERI_Px_Dyz_S_Py_C1000001001_b;
  abcd[1418] = 2.0E0*I_ERI_Py_Dyz_S_Py_C1000001001_b;
  abcd[1419] = 2.0E0*I_ERI_Pz_Dyz_S_Py_C1000001001_b;
  abcd[1421] = 2.0E0*I_ERI_Px_D2z_S_Py_C1000001001_b-1*I_ERI_Px_S_S_Py_C1000001001;
  abcd[1422] = 2.0E0*I_ERI_Py_D2z_S_Py_C1000001001_b-1*I_ERI_Py_S_S_Py_C1000001001;
  abcd[1423] = 2.0E0*I_ERI_Pz_D2z_S_Py_C1000001001_b-1*I_ERI_Pz_S_S_Py_C1000001001;
  abcd[1477] = 2.0E0*I_ERI_Px_Dxz_S_Pz_C1000001001_b;
  abcd[1478] = 2.0E0*I_ERI_Py_Dxz_S_Pz_C1000001001_b;
  abcd[1479] = 2.0E0*I_ERI_Pz_Dxz_S_Pz_C1000001001_b;
  abcd[1481] = 2.0E0*I_ERI_Px_Dyz_S_Pz_C1000001001_b;
  abcd[1482] = 2.0E0*I_ERI_Py_Dyz_S_Pz_C1000001001_b;
  abcd[1483] = 2.0E0*I_ERI_Pz_Dyz_S_Pz_C1000001001_b;
  abcd[1485] = 2.0E0*I_ERI_Px_D2z_S_Pz_C1000001001_b-1*I_ERI_Px_S_S_Pz_C1000001001;
  abcd[1486] = 2.0E0*I_ERI_Py_D2z_S_Pz_C1000001001_b-1*I_ERI_Py_S_S_Pz_C1000001001;
  abcd[1487] = 2.0E0*I_ERI_Pz_D2z_S_Pz_C1000001001_b-1*I_ERI_Pz_S_S_Pz_C1000001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000000_b
   ************************************************************/
  abcd[1360] = 2.0E0*I_ERI_S_Pz_Px_Px_C1001000000_b;
  abcd[1376] = 2.0E0*I_ERI_S_Pz_Py_Px_C1001000000_b;
  abcd[1392] = 2.0E0*I_ERI_S_Pz_Pz_Px_C1001000000_b;
  abcd[1424] = 2.0E0*I_ERI_S_Pz_Px_Py_C1001000000_b;
  abcd[1440] = 2.0E0*I_ERI_S_Pz_Py_Py_C1001000000_b;
  abcd[1456] = 2.0E0*I_ERI_S_Pz_Pz_Py_C1001000000_b;
  abcd[1488] = 2.0E0*I_ERI_S_Pz_Px_Pz_C1001000000_b;
  abcd[1504] = 2.0E0*I_ERI_S_Pz_Py_Pz_C1001000000_b;
  abcd[1520] = 2.0E0*I_ERI_S_Pz_Pz_Pz_C1001000000_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001000001_b
   ************************************************************/
  abcd[1361] = 2.0E0*I_ERI_Px_Pz_Px_Px_C1001000001_b;
  abcd[1362] = 2.0E0*I_ERI_Py_Pz_Px_Px_C1001000001_b;
  abcd[1363] = 2.0E0*I_ERI_Pz_Pz_Px_Px_C1001000001_b;
  abcd[1377] = 2.0E0*I_ERI_Px_Pz_Py_Px_C1001000001_b;
  abcd[1378] = 2.0E0*I_ERI_Py_Pz_Py_Px_C1001000001_b;
  abcd[1379] = 2.0E0*I_ERI_Pz_Pz_Py_Px_C1001000001_b;
  abcd[1393] = 2.0E0*I_ERI_Px_Pz_Pz_Px_C1001000001_b;
  abcd[1394] = 2.0E0*I_ERI_Py_Pz_Pz_Px_C1001000001_b;
  abcd[1395] = 2.0E0*I_ERI_Pz_Pz_Pz_Px_C1001000001_b;
  abcd[1425] = 2.0E0*I_ERI_Px_Pz_Px_Py_C1001000001_b;
  abcd[1426] = 2.0E0*I_ERI_Py_Pz_Px_Py_C1001000001_b;
  abcd[1427] = 2.0E0*I_ERI_Pz_Pz_Px_Py_C1001000001_b;
  abcd[1441] = 2.0E0*I_ERI_Px_Pz_Py_Py_C1001000001_b;
  abcd[1442] = 2.0E0*I_ERI_Py_Pz_Py_Py_C1001000001_b;
  abcd[1443] = 2.0E0*I_ERI_Pz_Pz_Py_Py_C1001000001_b;
  abcd[1457] = 2.0E0*I_ERI_Px_Pz_Pz_Py_C1001000001_b;
  abcd[1458] = 2.0E0*I_ERI_Py_Pz_Pz_Py_C1001000001_b;
  abcd[1459] = 2.0E0*I_ERI_Pz_Pz_Pz_Py_C1001000001_b;
  abcd[1489] = 2.0E0*I_ERI_Px_Pz_Px_Pz_C1001000001_b;
  abcd[1490] = 2.0E0*I_ERI_Py_Pz_Px_Pz_C1001000001_b;
  abcd[1491] = 2.0E0*I_ERI_Pz_Pz_Px_Pz_C1001000001_b;
  abcd[1505] = 2.0E0*I_ERI_Px_Pz_Py_Pz_C1001000001_b;
  abcd[1506] = 2.0E0*I_ERI_Py_Pz_Py_Pz_C1001000001_b;
  abcd[1507] = 2.0E0*I_ERI_Pz_Pz_Py_Pz_C1001000001_b;
  abcd[1521] = 2.0E0*I_ERI_Px_Pz_Pz_Pz_C1001000001_b;
  abcd[1522] = 2.0E0*I_ERI_Py_Pz_Pz_Pz_C1001000001_b;
  abcd[1523] = 2.0E0*I_ERI_Pz_Pz_Pz_Pz_C1001000001_b;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_D_P_P_C1001001000_b
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1001001000
   ************************************************************/
  abcd[1364] = 2.0E0*I_ERI_S_Dxz_Px_Px_C1001001000_b;
  abcd[1368] = 2.0E0*I_ERI_S_Dyz_Px_Px_C1001001000_b;
  abcd[1372] = 2.0E0*I_ERI_S_D2z_Px_Px_C1001001000_b-1*I_ERI_S_S_Px_Px_C1001001000;
  abcd[1380] = 2.0E0*I_ERI_S_Dxz_Py_Px_C1001001000_b;
  abcd[1384] = 2.0E0*I_ERI_S_Dyz_Py_Px_C1001001000_b;
  abcd[1388] = 2.0E0*I_ERI_S_D2z_Py_Px_C1001001000_b-1*I_ERI_S_S_Py_Px_C1001001000;
  abcd[1396] = 2.0E0*I_ERI_S_Dxz_Pz_Px_C1001001000_b;
  abcd[1400] = 2.0E0*I_ERI_S_Dyz_Pz_Px_C1001001000_b;
  abcd[1404] = 2.0E0*I_ERI_S_D2z_Pz_Px_C1001001000_b-1*I_ERI_S_S_Pz_Px_C1001001000;
  abcd[1428] = 2.0E0*I_ERI_S_Dxz_Px_Py_C1001001000_b;
  abcd[1432] = 2.0E0*I_ERI_S_Dyz_Px_Py_C1001001000_b;
  abcd[1436] = 2.0E0*I_ERI_S_D2z_Px_Py_C1001001000_b-1*I_ERI_S_S_Px_Py_C1001001000;
  abcd[1444] = 2.0E0*I_ERI_S_Dxz_Py_Py_C1001001000_b;
  abcd[1448] = 2.0E0*I_ERI_S_Dyz_Py_Py_C1001001000_b;
  abcd[1452] = 2.0E0*I_ERI_S_D2z_Py_Py_C1001001000_b-1*I_ERI_S_S_Py_Py_C1001001000;
  abcd[1460] = 2.0E0*I_ERI_S_Dxz_Pz_Py_C1001001000_b;
  abcd[1464] = 2.0E0*I_ERI_S_Dyz_Pz_Py_C1001001000_b;
  abcd[1468] = 2.0E0*I_ERI_S_D2z_Pz_Py_C1001001000_b-1*I_ERI_S_S_Pz_Py_C1001001000;
  abcd[1492] = 2.0E0*I_ERI_S_Dxz_Px_Pz_C1001001000_b;
  abcd[1496] = 2.0E0*I_ERI_S_Dyz_Px_Pz_C1001001000_b;
  abcd[1500] = 2.0E0*I_ERI_S_D2z_Px_Pz_C1001001000_b-1*I_ERI_S_S_Px_Pz_C1001001000;
  abcd[1508] = 2.0E0*I_ERI_S_Dxz_Py_Pz_C1001001000_b;
  abcd[1512] = 2.0E0*I_ERI_S_Dyz_Py_Pz_C1001001000_b;
  abcd[1516] = 2.0E0*I_ERI_S_D2z_Py_Pz_C1001001000_b-1*I_ERI_S_S_Py_Pz_C1001001000;
  abcd[1524] = 2.0E0*I_ERI_S_Dxz_Pz_Pz_C1001001000_b;
  abcd[1528] = 2.0E0*I_ERI_S_Dyz_Pz_Pz_C1001001000_b;
  abcd[1532] = 2.0E0*I_ERI_S_D2z_Pz_Pz_C1001001000_b-1*I_ERI_S_S_Pz_Pz_C1001001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_C1001001001_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001001
   ************************************************************/
  abcd[1365] = 2.0E0*I_ERI_Px_Dxz_Px_Px_C1001001001_b;
  abcd[1366] = 2.0E0*I_ERI_Py_Dxz_Px_Px_C1001001001_b;
  abcd[1367] = 2.0E0*I_ERI_Pz_Dxz_Px_Px_C1001001001_b;
  abcd[1369] = 2.0E0*I_ERI_Px_Dyz_Px_Px_C1001001001_b;
  abcd[1370] = 2.0E0*I_ERI_Py_Dyz_Px_Px_C1001001001_b;
  abcd[1371] = 2.0E0*I_ERI_Pz_Dyz_Px_Px_C1001001001_b;
  abcd[1373] = 2.0E0*I_ERI_Px_D2z_Px_Px_C1001001001_b-1*I_ERI_Px_S_Px_Px_C1001001001;
  abcd[1374] = 2.0E0*I_ERI_Py_D2z_Px_Px_C1001001001_b-1*I_ERI_Py_S_Px_Px_C1001001001;
  abcd[1375] = 2.0E0*I_ERI_Pz_D2z_Px_Px_C1001001001_b-1*I_ERI_Pz_S_Px_Px_C1001001001;
  abcd[1381] = 2.0E0*I_ERI_Px_Dxz_Py_Px_C1001001001_b;
  abcd[1382] = 2.0E0*I_ERI_Py_Dxz_Py_Px_C1001001001_b;
  abcd[1383] = 2.0E0*I_ERI_Pz_Dxz_Py_Px_C1001001001_b;
  abcd[1385] = 2.0E0*I_ERI_Px_Dyz_Py_Px_C1001001001_b;
  abcd[1386] = 2.0E0*I_ERI_Py_Dyz_Py_Px_C1001001001_b;
  abcd[1387] = 2.0E0*I_ERI_Pz_Dyz_Py_Px_C1001001001_b;
  abcd[1389] = 2.0E0*I_ERI_Px_D2z_Py_Px_C1001001001_b-1*I_ERI_Px_S_Py_Px_C1001001001;
  abcd[1390] = 2.0E0*I_ERI_Py_D2z_Py_Px_C1001001001_b-1*I_ERI_Py_S_Py_Px_C1001001001;
  abcd[1391] = 2.0E0*I_ERI_Pz_D2z_Py_Px_C1001001001_b-1*I_ERI_Pz_S_Py_Px_C1001001001;
  abcd[1397] = 2.0E0*I_ERI_Px_Dxz_Pz_Px_C1001001001_b;
  abcd[1398] = 2.0E0*I_ERI_Py_Dxz_Pz_Px_C1001001001_b;
  abcd[1399] = 2.0E0*I_ERI_Pz_Dxz_Pz_Px_C1001001001_b;
  abcd[1401] = 2.0E0*I_ERI_Px_Dyz_Pz_Px_C1001001001_b;
  abcd[1402] = 2.0E0*I_ERI_Py_Dyz_Pz_Px_C1001001001_b;
  abcd[1403] = 2.0E0*I_ERI_Pz_Dyz_Pz_Px_C1001001001_b;
  abcd[1405] = 2.0E0*I_ERI_Px_D2z_Pz_Px_C1001001001_b-1*I_ERI_Px_S_Pz_Px_C1001001001;
  abcd[1406] = 2.0E0*I_ERI_Py_D2z_Pz_Px_C1001001001_b-1*I_ERI_Py_S_Pz_Px_C1001001001;
  abcd[1407] = 2.0E0*I_ERI_Pz_D2z_Pz_Px_C1001001001_b-1*I_ERI_Pz_S_Pz_Px_C1001001001;
  abcd[1429] = 2.0E0*I_ERI_Px_Dxz_Px_Py_C1001001001_b;
  abcd[1430] = 2.0E0*I_ERI_Py_Dxz_Px_Py_C1001001001_b;
  abcd[1431] = 2.0E0*I_ERI_Pz_Dxz_Px_Py_C1001001001_b;
  abcd[1433] = 2.0E0*I_ERI_Px_Dyz_Px_Py_C1001001001_b;
  abcd[1434] = 2.0E0*I_ERI_Py_Dyz_Px_Py_C1001001001_b;
  abcd[1435] = 2.0E0*I_ERI_Pz_Dyz_Px_Py_C1001001001_b;
  abcd[1437] = 2.0E0*I_ERI_Px_D2z_Px_Py_C1001001001_b-1*I_ERI_Px_S_Px_Py_C1001001001;
  abcd[1438] = 2.0E0*I_ERI_Py_D2z_Px_Py_C1001001001_b-1*I_ERI_Py_S_Px_Py_C1001001001;
  abcd[1439] = 2.0E0*I_ERI_Pz_D2z_Px_Py_C1001001001_b-1*I_ERI_Pz_S_Px_Py_C1001001001;
  abcd[1445] = 2.0E0*I_ERI_Px_Dxz_Py_Py_C1001001001_b;
  abcd[1446] = 2.0E0*I_ERI_Py_Dxz_Py_Py_C1001001001_b;
  abcd[1447] = 2.0E0*I_ERI_Pz_Dxz_Py_Py_C1001001001_b;
  abcd[1449] = 2.0E0*I_ERI_Px_Dyz_Py_Py_C1001001001_b;
  abcd[1450] = 2.0E0*I_ERI_Py_Dyz_Py_Py_C1001001001_b;
  abcd[1451] = 2.0E0*I_ERI_Pz_Dyz_Py_Py_C1001001001_b;
  abcd[1453] = 2.0E0*I_ERI_Px_D2z_Py_Py_C1001001001_b-1*I_ERI_Px_S_Py_Py_C1001001001;
  abcd[1454] = 2.0E0*I_ERI_Py_D2z_Py_Py_C1001001001_b-1*I_ERI_Py_S_Py_Py_C1001001001;
  abcd[1455] = 2.0E0*I_ERI_Pz_D2z_Py_Py_C1001001001_b-1*I_ERI_Pz_S_Py_Py_C1001001001;
  abcd[1461] = 2.0E0*I_ERI_Px_Dxz_Pz_Py_C1001001001_b;
  abcd[1462] = 2.0E0*I_ERI_Py_Dxz_Pz_Py_C1001001001_b;
  abcd[1463] = 2.0E0*I_ERI_Pz_Dxz_Pz_Py_C1001001001_b;
  abcd[1465] = 2.0E0*I_ERI_Px_Dyz_Pz_Py_C1001001001_b;
  abcd[1466] = 2.0E0*I_ERI_Py_Dyz_Pz_Py_C1001001001_b;
  abcd[1467] = 2.0E0*I_ERI_Pz_Dyz_Pz_Py_C1001001001_b;
  abcd[1469] = 2.0E0*I_ERI_Px_D2z_Pz_Py_C1001001001_b-1*I_ERI_Px_S_Pz_Py_C1001001001;
  abcd[1470] = 2.0E0*I_ERI_Py_D2z_Pz_Py_C1001001001_b-1*I_ERI_Py_S_Pz_Py_C1001001001;
  abcd[1471] = 2.0E0*I_ERI_Pz_D2z_Pz_Py_C1001001001_b-1*I_ERI_Pz_S_Pz_Py_C1001001001;
  abcd[1493] = 2.0E0*I_ERI_Px_Dxz_Px_Pz_C1001001001_b;
  abcd[1494] = 2.0E0*I_ERI_Py_Dxz_Px_Pz_C1001001001_b;
  abcd[1495] = 2.0E0*I_ERI_Pz_Dxz_Px_Pz_C1001001001_b;
  abcd[1497] = 2.0E0*I_ERI_Px_Dyz_Px_Pz_C1001001001_b;
  abcd[1498] = 2.0E0*I_ERI_Py_Dyz_Px_Pz_C1001001001_b;
  abcd[1499] = 2.0E0*I_ERI_Pz_Dyz_Px_Pz_C1001001001_b;
  abcd[1501] = 2.0E0*I_ERI_Px_D2z_Px_Pz_C1001001001_b-1*I_ERI_Px_S_Px_Pz_C1001001001;
  abcd[1502] = 2.0E0*I_ERI_Py_D2z_Px_Pz_C1001001001_b-1*I_ERI_Py_S_Px_Pz_C1001001001;
  abcd[1503] = 2.0E0*I_ERI_Pz_D2z_Px_Pz_C1001001001_b-1*I_ERI_Pz_S_Px_Pz_C1001001001;
  abcd[1509] = 2.0E0*I_ERI_Px_Dxz_Py_Pz_C1001001001_b;
  abcd[1510] = 2.0E0*I_ERI_Py_Dxz_Py_Pz_C1001001001_b;
  abcd[1511] = 2.0E0*I_ERI_Pz_Dxz_Py_Pz_C1001001001_b;
  abcd[1513] = 2.0E0*I_ERI_Px_Dyz_Py_Pz_C1001001001_b;
  abcd[1514] = 2.0E0*I_ERI_Py_Dyz_Py_Pz_C1001001001_b;
  abcd[1515] = 2.0E0*I_ERI_Pz_Dyz_Py_Pz_C1001001001_b;
  abcd[1517] = 2.0E0*I_ERI_Px_D2z_Py_Pz_C1001001001_b-1*I_ERI_Px_S_Py_Pz_C1001001001;
  abcd[1518] = 2.0E0*I_ERI_Py_D2z_Py_Pz_C1001001001_b-1*I_ERI_Py_S_Py_Pz_C1001001001;
  abcd[1519] = 2.0E0*I_ERI_Pz_D2z_Py_Pz_C1001001001_b-1*I_ERI_Pz_S_Py_Pz_C1001001001;
  abcd[1525] = 2.0E0*I_ERI_Px_Dxz_Pz_Pz_C1001001001_b;
  abcd[1526] = 2.0E0*I_ERI_Py_Dxz_Pz_Pz_C1001001001_b;
  abcd[1527] = 2.0E0*I_ERI_Pz_Dxz_Pz_Pz_C1001001001_b;
  abcd[1529] = 2.0E0*I_ERI_Px_Dyz_Pz_Pz_C1001001001_b;
  abcd[1530] = 2.0E0*I_ERI_Py_Dyz_Pz_Pz_C1001001001_b;
  abcd[1531] = 2.0E0*I_ERI_Pz_Dyz_Pz_Pz_C1001001001_b;
  abcd[1533] = 2.0E0*I_ERI_Px_D2z_Pz_Pz_C1001001001_b-1*I_ERI_Px_S_Pz_Pz_C1001001001;
  abcd[1534] = 2.0E0*I_ERI_Py_D2z_Pz_Pz_C1001001001_b-1*I_ERI_Py_S_Pz_Pz_C1001001001;
  abcd[1535] = 2.0E0*I_ERI_Pz_D2z_Pz_Pz_C1001001001_b-1*I_ERI_Pz_S_Pz_Pz_C1001001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C0_c
   ************************************************************/
  abcd[1536] = 2.0E0*I_ERI_S_S_Px_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1_c
   ************************************************************/
  abcd[1537] = 2.0E0*I_ERI_Px_S_Px_S_C1_c;
  abcd[1538] = 2.0E0*I_ERI_Py_S_Px_S_C1_c;
  abcd[1539] = 2.0E0*I_ERI_Pz_S_Px_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000_c
   ************************************************************/
  abcd[1540] = 2.0E0*I_ERI_S_Px_Px_S_C1000_c;
  abcd[1544] = 2.0E0*I_ERI_S_Py_Px_S_C1000_c;
  abcd[1548] = 2.0E0*I_ERI_S_Pz_Px_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001_c
   ************************************************************/
  abcd[1541] = 2.0E0*I_ERI_Px_Px_Px_S_C1001_c;
  abcd[1542] = 2.0E0*I_ERI_Py_Px_Px_S_C1001_c;
  abcd[1543] = 2.0E0*I_ERI_Pz_Px_Px_S_C1001_c;
  abcd[1545] = 2.0E0*I_ERI_Px_Py_Px_S_C1001_c;
  abcd[1546] = 2.0E0*I_ERI_Py_Py_Px_S_C1001_c;
  abcd[1547] = 2.0E0*I_ERI_Pz_Py_Px_S_C1001_c;
  abcd[1549] = 2.0E0*I_ERI_Px_Pz_Px_S_C1001_c;
  abcd[1550] = 2.0E0*I_ERI_Py_Pz_Px_S_C1001_c;
  abcd[1551] = 2.0E0*I_ERI_Pz_Pz_Px_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[1552] = 2.0E0*I_ERI_S_S_D2x_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;
  abcd[1568] = 2.0E0*I_ERI_S_S_Dxy_S_C1000000_c;
  abcd[1584] = 2.0E0*I_ERI_S_S_Dxz_S_C1000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[1553] = 2.0E0*I_ERI_Px_S_D2x_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[1554] = 2.0E0*I_ERI_Py_S_D2x_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[1555] = 2.0E0*I_ERI_Pz_S_D2x_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[1569] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[1570] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[1571] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[1585] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[1586] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[1587] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[1556] = 2.0E0*I_ERI_S_Px_D2x_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[1560] = 2.0E0*I_ERI_S_Py_D2x_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[1564] = 2.0E0*I_ERI_S_Pz_D2x_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;
  abcd[1572] = 2.0E0*I_ERI_S_Px_Dxy_S_C1001000_c;
  abcd[1576] = 2.0E0*I_ERI_S_Py_Dxy_S_C1001000_c;
  abcd[1580] = 2.0E0*I_ERI_S_Pz_Dxy_S_C1001000_c;
  abcd[1588] = 2.0E0*I_ERI_S_Px_Dxz_S_C1001000_c;
  abcd[1592] = 2.0E0*I_ERI_S_Py_Dxz_S_C1001000_c;
  abcd[1596] = 2.0E0*I_ERI_S_Pz_Dxz_S_C1001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[1557] = 2.0E0*I_ERI_Px_Px_D2x_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[1558] = 2.0E0*I_ERI_Py_Px_D2x_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[1559] = 2.0E0*I_ERI_Pz_Px_D2x_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[1561] = 2.0E0*I_ERI_Px_Py_D2x_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[1562] = 2.0E0*I_ERI_Py_Py_D2x_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[1563] = 2.0E0*I_ERI_Pz_Py_D2x_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[1565] = 2.0E0*I_ERI_Px_Pz_D2x_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[1566] = 2.0E0*I_ERI_Py_Pz_D2x_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[1567] = 2.0E0*I_ERI_Pz_Pz_D2x_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[1573] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[1574] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[1575] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[1577] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[1578] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[1579] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[1581] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[1582] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[1583] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[1589] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[1590] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[1591] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[1593] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[1594] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[1595] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[1597] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[1598] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[1599] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000000_c
   ************************************************************/
  abcd[1600] = 2.0E0*I_ERI_S_S_Px_Px_C1000000000_c;
  abcd[1664] = 2.0E0*I_ERI_S_S_Px_Py_C1000000000_c;
  abcd[1728] = 2.0E0*I_ERI_S_S_Px_Pz_C1000000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000000001_c
   ************************************************************/
  abcd[1601] = 2.0E0*I_ERI_Px_S_Px_Px_C1000000001_c;
  abcd[1602] = 2.0E0*I_ERI_Py_S_Px_Px_C1000000001_c;
  abcd[1603] = 2.0E0*I_ERI_Pz_S_Px_Px_C1000000001_c;
  abcd[1665] = 2.0E0*I_ERI_Px_S_Px_Py_C1000000001_c;
  abcd[1666] = 2.0E0*I_ERI_Py_S_Px_Py_C1000000001_c;
  abcd[1667] = 2.0E0*I_ERI_Pz_S_Px_Py_C1000000001_c;
  abcd[1729] = 2.0E0*I_ERI_Px_S_Px_Pz_C1000000001_c;
  abcd[1730] = 2.0E0*I_ERI_Py_S_Px_Pz_C1000000001_c;
  abcd[1731] = 2.0E0*I_ERI_Pz_S_Px_Pz_C1000000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1000001000_c
   ************************************************************/
  abcd[1604] = 2.0E0*I_ERI_S_Px_Px_Px_C1000001000_c;
  abcd[1608] = 2.0E0*I_ERI_S_Py_Px_Px_C1000001000_c;
  abcd[1612] = 2.0E0*I_ERI_S_Pz_Px_Px_C1000001000_c;
  abcd[1668] = 2.0E0*I_ERI_S_Px_Px_Py_C1000001000_c;
  abcd[1672] = 2.0E0*I_ERI_S_Py_Px_Py_C1000001000_c;
  abcd[1676] = 2.0E0*I_ERI_S_Pz_Px_Py_C1000001000_c;
  abcd[1732] = 2.0E0*I_ERI_S_Px_Px_Pz_C1000001000_c;
  abcd[1736] = 2.0E0*I_ERI_S_Py_Px_Pz_C1000001000_c;
  abcd[1740] = 2.0E0*I_ERI_S_Pz_Px_Pz_C1000001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1000001001_c
   ************************************************************/
  abcd[1605] = 2.0E0*I_ERI_Px_Px_Px_Px_C1000001001_c;
  abcd[1606] = 2.0E0*I_ERI_Py_Px_Px_Px_C1000001001_c;
  abcd[1607] = 2.0E0*I_ERI_Pz_Px_Px_Px_C1000001001_c;
  abcd[1609] = 2.0E0*I_ERI_Px_Py_Px_Px_C1000001001_c;
  abcd[1610] = 2.0E0*I_ERI_Py_Py_Px_Px_C1000001001_c;
  abcd[1611] = 2.0E0*I_ERI_Pz_Py_Px_Px_C1000001001_c;
  abcd[1613] = 2.0E0*I_ERI_Px_Pz_Px_Px_C1000001001_c;
  abcd[1614] = 2.0E0*I_ERI_Py_Pz_Px_Px_C1000001001_c;
  abcd[1615] = 2.0E0*I_ERI_Pz_Pz_Px_Px_C1000001001_c;
  abcd[1669] = 2.0E0*I_ERI_Px_Px_Px_Py_C1000001001_c;
  abcd[1670] = 2.0E0*I_ERI_Py_Px_Px_Py_C1000001001_c;
  abcd[1671] = 2.0E0*I_ERI_Pz_Px_Px_Py_C1000001001_c;
  abcd[1673] = 2.0E0*I_ERI_Px_Py_Px_Py_C1000001001_c;
  abcd[1674] = 2.0E0*I_ERI_Py_Py_Px_Py_C1000001001_c;
  abcd[1675] = 2.0E0*I_ERI_Pz_Py_Px_Py_C1000001001_c;
  abcd[1677] = 2.0E0*I_ERI_Px_Pz_Px_Py_C1000001001_c;
  abcd[1678] = 2.0E0*I_ERI_Py_Pz_Px_Py_C1000001001_c;
  abcd[1679] = 2.0E0*I_ERI_Pz_Pz_Px_Py_C1000001001_c;
  abcd[1733] = 2.0E0*I_ERI_Px_Px_Px_Pz_C1000001001_c;
  abcd[1734] = 2.0E0*I_ERI_Py_Px_Px_Pz_C1000001001_c;
  abcd[1735] = 2.0E0*I_ERI_Pz_Px_Px_Pz_C1000001001_c;
  abcd[1737] = 2.0E0*I_ERI_Px_Py_Px_Pz_C1000001001_c;
  abcd[1738] = 2.0E0*I_ERI_Py_Py_Px_Pz_C1000001001_c;
  abcd[1739] = 2.0E0*I_ERI_Pz_Py_Px_Pz_C1000001001_c;
  abcd[1741] = 2.0E0*I_ERI_Px_Pz_Px_Pz_C1000001001_c;
  abcd[1742] = 2.0E0*I_ERI_Py_Pz_Px_Pz_C1000001001_c;
  abcd[1743] = 2.0E0*I_ERI_Pz_Pz_Px_Pz_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_P_C1001000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1001000000
   ************************************************************/
  abcd[1616] = 2.0E0*I_ERI_S_S_D2x_Px_C1001000000_c-1*I_ERI_S_S_S_Px_C1001000000;
  abcd[1632] = 2.0E0*I_ERI_S_S_Dxy_Px_C1001000000_c;
  abcd[1648] = 2.0E0*I_ERI_S_S_Dxz_Px_C1001000000_c;
  abcd[1680] = 2.0E0*I_ERI_S_S_D2x_Py_C1001000000_c-1*I_ERI_S_S_S_Py_C1001000000;
  abcd[1696] = 2.0E0*I_ERI_S_S_Dxy_Py_C1001000000_c;
  abcd[1712] = 2.0E0*I_ERI_S_S_Dxz_Py_C1001000000_c;
  abcd[1744] = 2.0E0*I_ERI_S_S_D2x_Pz_C1001000000_c-1*I_ERI_S_S_S_Pz_C1001000000;
  abcd[1760] = 2.0E0*I_ERI_S_S_Dxy_Pz_C1001000000_c;
  abcd[1776] = 2.0E0*I_ERI_S_S_Dxz_Pz_C1001000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1001000001
   ************************************************************/
  abcd[1617] = 2.0E0*I_ERI_Px_S_D2x_Px_C1001000001_c-1*I_ERI_Px_S_S_Px_C1001000001;
  abcd[1618] = 2.0E0*I_ERI_Py_S_D2x_Px_C1001000001_c-1*I_ERI_Py_S_S_Px_C1001000001;
  abcd[1619] = 2.0E0*I_ERI_Pz_S_D2x_Px_C1001000001_c-1*I_ERI_Pz_S_S_Px_C1001000001;
  abcd[1633] = 2.0E0*I_ERI_Px_S_Dxy_Px_C1001000001_c;
  abcd[1634] = 2.0E0*I_ERI_Py_S_Dxy_Px_C1001000001_c;
  abcd[1635] = 2.0E0*I_ERI_Pz_S_Dxy_Px_C1001000001_c;
  abcd[1649] = 2.0E0*I_ERI_Px_S_Dxz_Px_C1001000001_c;
  abcd[1650] = 2.0E0*I_ERI_Py_S_Dxz_Px_C1001000001_c;
  abcd[1651] = 2.0E0*I_ERI_Pz_S_Dxz_Px_C1001000001_c;
  abcd[1681] = 2.0E0*I_ERI_Px_S_D2x_Py_C1001000001_c-1*I_ERI_Px_S_S_Py_C1001000001;
  abcd[1682] = 2.0E0*I_ERI_Py_S_D2x_Py_C1001000001_c-1*I_ERI_Py_S_S_Py_C1001000001;
  abcd[1683] = 2.0E0*I_ERI_Pz_S_D2x_Py_C1001000001_c-1*I_ERI_Pz_S_S_Py_C1001000001;
  abcd[1697] = 2.0E0*I_ERI_Px_S_Dxy_Py_C1001000001_c;
  abcd[1698] = 2.0E0*I_ERI_Py_S_Dxy_Py_C1001000001_c;
  abcd[1699] = 2.0E0*I_ERI_Pz_S_Dxy_Py_C1001000001_c;
  abcd[1713] = 2.0E0*I_ERI_Px_S_Dxz_Py_C1001000001_c;
  abcd[1714] = 2.0E0*I_ERI_Py_S_Dxz_Py_C1001000001_c;
  abcd[1715] = 2.0E0*I_ERI_Pz_S_Dxz_Py_C1001000001_c;
  abcd[1745] = 2.0E0*I_ERI_Px_S_D2x_Pz_C1001000001_c-1*I_ERI_Px_S_S_Pz_C1001000001;
  abcd[1746] = 2.0E0*I_ERI_Py_S_D2x_Pz_C1001000001_c-1*I_ERI_Py_S_S_Pz_C1001000001;
  abcd[1747] = 2.0E0*I_ERI_Pz_S_D2x_Pz_C1001000001_c-1*I_ERI_Pz_S_S_Pz_C1001000001;
  abcd[1761] = 2.0E0*I_ERI_Px_S_Dxy_Pz_C1001000001_c;
  abcd[1762] = 2.0E0*I_ERI_Py_S_Dxy_Pz_C1001000001_c;
  abcd[1763] = 2.0E0*I_ERI_Pz_S_Dxy_Pz_C1001000001_c;
  abcd[1777] = 2.0E0*I_ERI_Px_S_Dxz_Pz_C1001000001_c;
  abcd[1778] = 2.0E0*I_ERI_Py_S_Dxz_Pz_C1001000001_c;
  abcd[1779] = 2.0E0*I_ERI_Pz_S_Dxz_Pz_C1001000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_P_C1001001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1001001000
   ************************************************************/
  abcd[1620] = 2.0E0*I_ERI_S_Px_D2x_Px_C1001001000_c-1*I_ERI_S_Px_S_Px_C1001001000;
  abcd[1624] = 2.0E0*I_ERI_S_Py_D2x_Px_C1001001000_c-1*I_ERI_S_Py_S_Px_C1001001000;
  abcd[1628] = 2.0E0*I_ERI_S_Pz_D2x_Px_C1001001000_c-1*I_ERI_S_Pz_S_Px_C1001001000;
  abcd[1636] = 2.0E0*I_ERI_S_Px_Dxy_Px_C1001001000_c;
  abcd[1640] = 2.0E0*I_ERI_S_Py_Dxy_Px_C1001001000_c;
  abcd[1644] = 2.0E0*I_ERI_S_Pz_Dxy_Px_C1001001000_c;
  abcd[1652] = 2.0E0*I_ERI_S_Px_Dxz_Px_C1001001000_c;
  abcd[1656] = 2.0E0*I_ERI_S_Py_Dxz_Px_C1001001000_c;
  abcd[1660] = 2.0E0*I_ERI_S_Pz_Dxz_Px_C1001001000_c;
  abcd[1684] = 2.0E0*I_ERI_S_Px_D2x_Py_C1001001000_c-1*I_ERI_S_Px_S_Py_C1001001000;
  abcd[1688] = 2.0E0*I_ERI_S_Py_D2x_Py_C1001001000_c-1*I_ERI_S_Py_S_Py_C1001001000;
  abcd[1692] = 2.0E0*I_ERI_S_Pz_D2x_Py_C1001001000_c-1*I_ERI_S_Pz_S_Py_C1001001000;
  abcd[1700] = 2.0E0*I_ERI_S_Px_Dxy_Py_C1001001000_c;
  abcd[1704] = 2.0E0*I_ERI_S_Py_Dxy_Py_C1001001000_c;
  abcd[1708] = 2.0E0*I_ERI_S_Pz_Dxy_Py_C1001001000_c;
  abcd[1716] = 2.0E0*I_ERI_S_Px_Dxz_Py_C1001001000_c;
  abcd[1720] = 2.0E0*I_ERI_S_Py_Dxz_Py_C1001001000_c;
  abcd[1724] = 2.0E0*I_ERI_S_Pz_Dxz_Py_C1001001000_c;
  abcd[1748] = 2.0E0*I_ERI_S_Px_D2x_Pz_C1001001000_c-1*I_ERI_S_Px_S_Pz_C1001001000;
  abcd[1752] = 2.0E0*I_ERI_S_Py_D2x_Pz_C1001001000_c-1*I_ERI_S_Py_S_Pz_C1001001000;
  abcd[1756] = 2.0E0*I_ERI_S_Pz_D2x_Pz_C1001001000_c-1*I_ERI_S_Pz_S_Pz_C1001001000;
  abcd[1764] = 2.0E0*I_ERI_S_Px_Dxy_Pz_C1001001000_c;
  abcd[1768] = 2.0E0*I_ERI_S_Py_Dxy_Pz_C1001001000_c;
  abcd[1772] = 2.0E0*I_ERI_S_Pz_Dxy_Pz_C1001001000_c;
  abcd[1780] = 2.0E0*I_ERI_S_Px_Dxz_Pz_C1001001000_c;
  abcd[1784] = 2.0E0*I_ERI_S_Py_Dxz_Pz_C1001001000_c;
  abcd[1788] = 2.0E0*I_ERI_S_Pz_Dxz_Pz_C1001001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001001
   ************************************************************/
  abcd[1621] = 2.0E0*I_ERI_Px_Px_D2x_Px_C1001001001_c-1*I_ERI_Px_Px_S_Px_C1001001001;
  abcd[1622] = 2.0E0*I_ERI_Py_Px_D2x_Px_C1001001001_c-1*I_ERI_Py_Px_S_Px_C1001001001;
  abcd[1623] = 2.0E0*I_ERI_Pz_Px_D2x_Px_C1001001001_c-1*I_ERI_Pz_Px_S_Px_C1001001001;
  abcd[1625] = 2.0E0*I_ERI_Px_Py_D2x_Px_C1001001001_c-1*I_ERI_Px_Py_S_Px_C1001001001;
  abcd[1626] = 2.0E0*I_ERI_Py_Py_D2x_Px_C1001001001_c-1*I_ERI_Py_Py_S_Px_C1001001001;
  abcd[1627] = 2.0E0*I_ERI_Pz_Py_D2x_Px_C1001001001_c-1*I_ERI_Pz_Py_S_Px_C1001001001;
  abcd[1629] = 2.0E0*I_ERI_Px_Pz_D2x_Px_C1001001001_c-1*I_ERI_Px_Pz_S_Px_C1001001001;
  abcd[1630] = 2.0E0*I_ERI_Py_Pz_D2x_Px_C1001001001_c-1*I_ERI_Py_Pz_S_Px_C1001001001;
  abcd[1631] = 2.0E0*I_ERI_Pz_Pz_D2x_Px_C1001001001_c-1*I_ERI_Pz_Pz_S_Px_C1001001001;
  abcd[1637] = 2.0E0*I_ERI_Px_Px_Dxy_Px_C1001001001_c;
  abcd[1638] = 2.0E0*I_ERI_Py_Px_Dxy_Px_C1001001001_c;
  abcd[1639] = 2.0E0*I_ERI_Pz_Px_Dxy_Px_C1001001001_c;
  abcd[1641] = 2.0E0*I_ERI_Px_Py_Dxy_Px_C1001001001_c;
  abcd[1642] = 2.0E0*I_ERI_Py_Py_Dxy_Px_C1001001001_c;
  abcd[1643] = 2.0E0*I_ERI_Pz_Py_Dxy_Px_C1001001001_c;
  abcd[1645] = 2.0E0*I_ERI_Px_Pz_Dxy_Px_C1001001001_c;
  abcd[1646] = 2.0E0*I_ERI_Py_Pz_Dxy_Px_C1001001001_c;
  abcd[1647] = 2.0E0*I_ERI_Pz_Pz_Dxy_Px_C1001001001_c;
  abcd[1653] = 2.0E0*I_ERI_Px_Px_Dxz_Px_C1001001001_c;
  abcd[1654] = 2.0E0*I_ERI_Py_Px_Dxz_Px_C1001001001_c;
  abcd[1655] = 2.0E0*I_ERI_Pz_Px_Dxz_Px_C1001001001_c;
  abcd[1657] = 2.0E0*I_ERI_Px_Py_Dxz_Px_C1001001001_c;
  abcd[1658] = 2.0E0*I_ERI_Py_Py_Dxz_Px_C1001001001_c;
  abcd[1659] = 2.0E0*I_ERI_Pz_Py_Dxz_Px_C1001001001_c;
  abcd[1661] = 2.0E0*I_ERI_Px_Pz_Dxz_Px_C1001001001_c;
  abcd[1662] = 2.0E0*I_ERI_Py_Pz_Dxz_Px_C1001001001_c;
  abcd[1663] = 2.0E0*I_ERI_Pz_Pz_Dxz_Px_C1001001001_c;
  abcd[1685] = 2.0E0*I_ERI_Px_Px_D2x_Py_C1001001001_c-1*I_ERI_Px_Px_S_Py_C1001001001;
  abcd[1686] = 2.0E0*I_ERI_Py_Px_D2x_Py_C1001001001_c-1*I_ERI_Py_Px_S_Py_C1001001001;
  abcd[1687] = 2.0E0*I_ERI_Pz_Px_D2x_Py_C1001001001_c-1*I_ERI_Pz_Px_S_Py_C1001001001;
  abcd[1689] = 2.0E0*I_ERI_Px_Py_D2x_Py_C1001001001_c-1*I_ERI_Px_Py_S_Py_C1001001001;
  abcd[1690] = 2.0E0*I_ERI_Py_Py_D2x_Py_C1001001001_c-1*I_ERI_Py_Py_S_Py_C1001001001;
  abcd[1691] = 2.0E0*I_ERI_Pz_Py_D2x_Py_C1001001001_c-1*I_ERI_Pz_Py_S_Py_C1001001001;
  abcd[1693] = 2.0E0*I_ERI_Px_Pz_D2x_Py_C1001001001_c-1*I_ERI_Px_Pz_S_Py_C1001001001;
  abcd[1694] = 2.0E0*I_ERI_Py_Pz_D2x_Py_C1001001001_c-1*I_ERI_Py_Pz_S_Py_C1001001001;
  abcd[1695] = 2.0E0*I_ERI_Pz_Pz_D2x_Py_C1001001001_c-1*I_ERI_Pz_Pz_S_Py_C1001001001;
  abcd[1701] = 2.0E0*I_ERI_Px_Px_Dxy_Py_C1001001001_c;
  abcd[1702] = 2.0E0*I_ERI_Py_Px_Dxy_Py_C1001001001_c;
  abcd[1703] = 2.0E0*I_ERI_Pz_Px_Dxy_Py_C1001001001_c;
  abcd[1705] = 2.0E0*I_ERI_Px_Py_Dxy_Py_C1001001001_c;
  abcd[1706] = 2.0E0*I_ERI_Py_Py_Dxy_Py_C1001001001_c;
  abcd[1707] = 2.0E0*I_ERI_Pz_Py_Dxy_Py_C1001001001_c;
  abcd[1709] = 2.0E0*I_ERI_Px_Pz_Dxy_Py_C1001001001_c;
  abcd[1710] = 2.0E0*I_ERI_Py_Pz_Dxy_Py_C1001001001_c;
  abcd[1711] = 2.0E0*I_ERI_Pz_Pz_Dxy_Py_C1001001001_c;
  abcd[1717] = 2.0E0*I_ERI_Px_Px_Dxz_Py_C1001001001_c;
  abcd[1718] = 2.0E0*I_ERI_Py_Px_Dxz_Py_C1001001001_c;
  abcd[1719] = 2.0E0*I_ERI_Pz_Px_Dxz_Py_C1001001001_c;
  abcd[1721] = 2.0E0*I_ERI_Px_Py_Dxz_Py_C1001001001_c;
  abcd[1722] = 2.0E0*I_ERI_Py_Py_Dxz_Py_C1001001001_c;
  abcd[1723] = 2.0E0*I_ERI_Pz_Py_Dxz_Py_C1001001001_c;
  abcd[1725] = 2.0E0*I_ERI_Px_Pz_Dxz_Py_C1001001001_c;
  abcd[1726] = 2.0E0*I_ERI_Py_Pz_Dxz_Py_C1001001001_c;
  abcd[1727] = 2.0E0*I_ERI_Pz_Pz_Dxz_Py_C1001001001_c;
  abcd[1749] = 2.0E0*I_ERI_Px_Px_D2x_Pz_C1001001001_c-1*I_ERI_Px_Px_S_Pz_C1001001001;
  abcd[1750] = 2.0E0*I_ERI_Py_Px_D2x_Pz_C1001001001_c-1*I_ERI_Py_Px_S_Pz_C1001001001;
  abcd[1751] = 2.0E0*I_ERI_Pz_Px_D2x_Pz_C1001001001_c-1*I_ERI_Pz_Px_S_Pz_C1001001001;
  abcd[1753] = 2.0E0*I_ERI_Px_Py_D2x_Pz_C1001001001_c-1*I_ERI_Px_Py_S_Pz_C1001001001;
  abcd[1754] = 2.0E0*I_ERI_Py_Py_D2x_Pz_C1001001001_c-1*I_ERI_Py_Py_S_Pz_C1001001001;
  abcd[1755] = 2.0E0*I_ERI_Pz_Py_D2x_Pz_C1001001001_c-1*I_ERI_Pz_Py_S_Pz_C1001001001;
  abcd[1757] = 2.0E0*I_ERI_Px_Pz_D2x_Pz_C1001001001_c-1*I_ERI_Px_Pz_S_Pz_C1001001001;
  abcd[1758] = 2.0E0*I_ERI_Py_Pz_D2x_Pz_C1001001001_c-1*I_ERI_Py_Pz_S_Pz_C1001001001;
  abcd[1759] = 2.0E0*I_ERI_Pz_Pz_D2x_Pz_C1001001001_c-1*I_ERI_Pz_Pz_S_Pz_C1001001001;
  abcd[1765] = 2.0E0*I_ERI_Px_Px_Dxy_Pz_C1001001001_c;
  abcd[1766] = 2.0E0*I_ERI_Py_Px_Dxy_Pz_C1001001001_c;
  abcd[1767] = 2.0E0*I_ERI_Pz_Px_Dxy_Pz_C1001001001_c;
  abcd[1769] = 2.0E0*I_ERI_Px_Py_Dxy_Pz_C1001001001_c;
  abcd[1770] = 2.0E0*I_ERI_Py_Py_Dxy_Pz_C1001001001_c;
  abcd[1771] = 2.0E0*I_ERI_Pz_Py_Dxy_Pz_C1001001001_c;
  abcd[1773] = 2.0E0*I_ERI_Px_Pz_Dxy_Pz_C1001001001_c;
  abcd[1774] = 2.0E0*I_ERI_Py_Pz_Dxy_Pz_C1001001001_c;
  abcd[1775] = 2.0E0*I_ERI_Pz_Pz_Dxy_Pz_C1001001001_c;
  abcd[1781] = 2.0E0*I_ERI_Px_Px_Dxz_Pz_C1001001001_c;
  abcd[1782] = 2.0E0*I_ERI_Py_Px_Dxz_Pz_C1001001001_c;
  abcd[1783] = 2.0E0*I_ERI_Pz_Px_Dxz_Pz_C1001001001_c;
  abcd[1785] = 2.0E0*I_ERI_Px_Py_Dxz_Pz_C1001001001_c;
  abcd[1786] = 2.0E0*I_ERI_Py_Py_Dxz_Pz_C1001001001_c;
  abcd[1787] = 2.0E0*I_ERI_Pz_Py_Dxz_Pz_C1001001001_c;
  abcd[1789] = 2.0E0*I_ERI_Px_Pz_Dxz_Pz_C1001001001_c;
  abcd[1790] = 2.0E0*I_ERI_Py_Pz_Dxz_Pz_C1001001001_c;
  abcd[1791] = 2.0E0*I_ERI_Pz_Pz_Dxz_Pz_C1001001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C0_c
   ************************************************************/
  abcd[1792] = 2.0E0*I_ERI_S_S_Py_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1_c
   ************************************************************/
  abcd[1793] = 2.0E0*I_ERI_Px_S_Py_S_C1_c;
  abcd[1794] = 2.0E0*I_ERI_Py_S_Py_S_C1_c;
  abcd[1795] = 2.0E0*I_ERI_Pz_S_Py_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000_c
   ************************************************************/
  abcd[1796] = 2.0E0*I_ERI_S_Px_Py_S_C1000_c;
  abcd[1800] = 2.0E0*I_ERI_S_Py_Py_S_C1000_c;
  abcd[1804] = 2.0E0*I_ERI_S_Pz_Py_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001_c
   ************************************************************/
  abcd[1797] = 2.0E0*I_ERI_Px_Px_Py_S_C1001_c;
  abcd[1798] = 2.0E0*I_ERI_Py_Px_Py_S_C1001_c;
  abcd[1799] = 2.0E0*I_ERI_Pz_Px_Py_S_C1001_c;
  abcd[1801] = 2.0E0*I_ERI_Px_Py_Py_S_C1001_c;
  abcd[1802] = 2.0E0*I_ERI_Py_Py_Py_S_C1001_c;
  abcd[1803] = 2.0E0*I_ERI_Pz_Py_Py_S_C1001_c;
  abcd[1805] = 2.0E0*I_ERI_Px_Pz_Py_S_C1001_c;
  abcd[1806] = 2.0E0*I_ERI_Py_Pz_Py_S_C1001_c;
  abcd[1807] = 2.0E0*I_ERI_Pz_Pz_Py_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[1808] = 2.0E0*I_ERI_S_S_Dxy_S_C1000000_c;
  abcd[1824] = 2.0E0*I_ERI_S_S_D2y_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;
  abcd[1840] = 2.0E0*I_ERI_S_S_Dyz_S_C1000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[1809] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[1810] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[1811] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[1825] = 2.0E0*I_ERI_Px_S_D2y_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[1826] = 2.0E0*I_ERI_Py_S_D2y_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[1827] = 2.0E0*I_ERI_Pz_S_D2y_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[1841] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[1842] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[1843] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[1812] = 2.0E0*I_ERI_S_Px_Dxy_S_C1001000_c;
  abcd[1816] = 2.0E0*I_ERI_S_Py_Dxy_S_C1001000_c;
  abcd[1820] = 2.0E0*I_ERI_S_Pz_Dxy_S_C1001000_c;
  abcd[1828] = 2.0E0*I_ERI_S_Px_D2y_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[1832] = 2.0E0*I_ERI_S_Py_D2y_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[1836] = 2.0E0*I_ERI_S_Pz_D2y_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;
  abcd[1844] = 2.0E0*I_ERI_S_Px_Dyz_S_C1001000_c;
  abcd[1848] = 2.0E0*I_ERI_S_Py_Dyz_S_C1001000_c;
  abcd[1852] = 2.0E0*I_ERI_S_Pz_Dyz_S_C1001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[1813] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[1814] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[1815] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[1817] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[1818] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[1819] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[1821] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[1822] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[1823] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[1829] = 2.0E0*I_ERI_Px_Px_D2y_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[1830] = 2.0E0*I_ERI_Py_Px_D2y_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[1831] = 2.0E0*I_ERI_Pz_Px_D2y_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[1833] = 2.0E0*I_ERI_Px_Py_D2y_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[1834] = 2.0E0*I_ERI_Py_Py_D2y_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[1835] = 2.0E0*I_ERI_Pz_Py_D2y_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[1837] = 2.0E0*I_ERI_Px_Pz_D2y_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[1838] = 2.0E0*I_ERI_Py_Pz_D2y_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[1839] = 2.0E0*I_ERI_Pz_Pz_D2y_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[1845] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[1846] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[1847] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[1849] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[1850] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[1851] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[1853] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[1854] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[1855] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000000_c
   ************************************************************/
  abcd[1856] = 2.0E0*I_ERI_S_S_Py_Px_C1000000000_c;
  abcd[1920] = 2.0E0*I_ERI_S_S_Py_Py_C1000000000_c;
  abcd[1984] = 2.0E0*I_ERI_S_S_Py_Pz_C1000000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000000001_c
   ************************************************************/
  abcd[1857] = 2.0E0*I_ERI_Px_S_Py_Px_C1000000001_c;
  abcd[1858] = 2.0E0*I_ERI_Py_S_Py_Px_C1000000001_c;
  abcd[1859] = 2.0E0*I_ERI_Pz_S_Py_Px_C1000000001_c;
  abcd[1921] = 2.0E0*I_ERI_Px_S_Py_Py_C1000000001_c;
  abcd[1922] = 2.0E0*I_ERI_Py_S_Py_Py_C1000000001_c;
  abcd[1923] = 2.0E0*I_ERI_Pz_S_Py_Py_C1000000001_c;
  abcd[1985] = 2.0E0*I_ERI_Px_S_Py_Pz_C1000000001_c;
  abcd[1986] = 2.0E0*I_ERI_Py_S_Py_Pz_C1000000001_c;
  abcd[1987] = 2.0E0*I_ERI_Pz_S_Py_Pz_C1000000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1000001000_c
   ************************************************************/
  abcd[1860] = 2.0E0*I_ERI_S_Px_Py_Px_C1000001000_c;
  abcd[1864] = 2.0E0*I_ERI_S_Py_Py_Px_C1000001000_c;
  abcd[1868] = 2.0E0*I_ERI_S_Pz_Py_Px_C1000001000_c;
  abcd[1924] = 2.0E0*I_ERI_S_Px_Py_Py_C1000001000_c;
  abcd[1928] = 2.0E0*I_ERI_S_Py_Py_Py_C1000001000_c;
  abcd[1932] = 2.0E0*I_ERI_S_Pz_Py_Py_C1000001000_c;
  abcd[1988] = 2.0E0*I_ERI_S_Px_Py_Pz_C1000001000_c;
  abcd[1992] = 2.0E0*I_ERI_S_Py_Py_Pz_C1000001000_c;
  abcd[1996] = 2.0E0*I_ERI_S_Pz_Py_Pz_C1000001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1000001001_c
   ************************************************************/
  abcd[1861] = 2.0E0*I_ERI_Px_Px_Py_Px_C1000001001_c;
  abcd[1862] = 2.0E0*I_ERI_Py_Px_Py_Px_C1000001001_c;
  abcd[1863] = 2.0E0*I_ERI_Pz_Px_Py_Px_C1000001001_c;
  abcd[1865] = 2.0E0*I_ERI_Px_Py_Py_Px_C1000001001_c;
  abcd[1866] = 2.0E0*I_ERI_Py_Py_Py_Px_C1000001001_c;
  abcd[1867] = 2.0E0*I_ERI_Pz_Py_Py_Px_C1000001001_c;
  abcd[1869] = 2.0E0*I_ERI_Px_Pz_Py_Px_C1000001001_c;
  abcd[1870] = 2.0E0*I_ERI_Py_Pz_Py_Px_C1000001001_c;
  abcd[1871] = 2.0E0*I_ERI_Pz_Pz_Py_Px_C1000001001_c;
  abcd[1925] = 2.0E0*I_ERI_Px_Px_Py_Py_C1000001001_c;
  abcd[1926] = 2.0E0*I_ERI_Py_Px_Py_Py_C1000001001_c;
  abcd[1927] = 2.0E0*I_ERI_Pz_Px_Py_Py_C1000001001_c;
  abcd[1929] = 2.0E0*I_ERI_Px_Py_Py_Py_C1000001001_c;
  abcd[1930] = 2.0E0*I_ERI_Py_Py_Py_Py_C1000001001_c;
  abcd[1931] = 2.0E0*I_ERI_Pz_Py_Py_Py_C1000001001_c;
  abcd[1933] = 2.0E0*I_ERI_Px_Pz_Py_Py_C1000001001_c;
  abcd[1934] = 2.0E0*I_ERI_Py_Pz_Py_Py_C1000001001_c;
  abcd[1935] = 2.0E0*I_ERI_Pz_Pz_Py_Py_C1000001001_c;
  abcd[1989] = 2.0E0*I_ERI_Px_Px_Py_Pz_C1000001001_c;
  abcd[1990] = 2.0E0*I_ERI_Py_Px_Py_Pz_C1000001001_c;
  abcd[1991] = 2.0E0*I_ERI_Pz_Px_Py_Pz_C1000001001_c;
  abcd[1993] = 2.0E0*I_ERI_Px_Py_Py_Pz_C1000001001_c;
  abcd[1994] = 2.0E0*I_ERI_Py_Py_Py_Pz_C1000001001_c;
  abcd[1995] = 2.0E0*I_ERI_Pz_Py_Py_Pz_C1000001001_c;
  abcd[1997] = 2.0E0*I_ERI_Px_Pz_Py_Pz_C1000001001_c;
  abcd[1998] = 2.0E0*I_ERI_Py_Pz_Py_Pz_C1000001001_c;
  abcd[1999] = 2.0E0*I_ERI_Pz_Pz_Py_Pz_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_P_C1001000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1001000000
   ************************************************************/
  abcd[1872] = 2.0E0*I_ERI_S_S_Dxy_Px_C1001000000_c;
  abcd[1888] = 2.0E0*I_ERI_S_S_D2y_Px_C1001000000_c-1*I_ERI_S_S_S_Px_C1001000000;
  abcd[1904] = 2.0E0*I_ERI_S_S_Dyz_Px_C1001000000_c;
  abcd[1936] = 2.0E0*I_ERI_S_S_Dxy_Py_C1001000000_c;
  abcd[1952] = 2.0E0*I_ERI_S_S_D2y_Py_C1001000000_c-1*I_ERI_S_S_S_Py_C1001000000;
  abcd[1968] = 2.0E0*I_ERI_S_S_Dyz_Py_C1001000000_c;
  abcd[2000] = 2.0E0*I_ERI_S_S_Dxy_Pz_C1001000000_c;
  abcd[2016] = 2.0E0*I_ERI_S_S_D2y_Pz_C1001000000_c-1*I_ERI_S_S_S_Pz_C1001000000;
  abcd[2032] = 2.0E0*I_ERI_S_S_Dyz_Pz_C1001000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1001000001
   ************************************************************/
  abcd[1873] = 2.0E0*I_ERI_Px_S_Dxy_Px_C1001000001_c;
  abcd[1874] = 2.0E0*I_ERI_Py_S_Dxy_Px_C1001000001_c;
  abcd[1875] = 2.0E0*I_ERI_Pz_S_Dxy_Px_C1001000001_c;
  abcd[1889] = 2.0E0*I_ERI_Px_S_D2y_Px_C1001000001_c-1*I_ERI_Px_S_S_Px_C1001000001;
  abcd[1890] = 2.0E0*I_ERI_Py_S_D2y_Px_C1001000001_c-1*I_ERI_Py_S_S_Px_C1001000001;
  abcd[1891] = 2.0E0*I_ERI_Pz_S_D2y_Px_C1001000001_c-1*I_ERI_Pz_S_S_Px_C1001000001;
  abcd[1905] = 2.0E0*I_ERI_Px_S_Dyz_Px_C1001000001_c;
  abcd[1906] = 2.0E0*I_ERI_Py_S_Dyz_Px_C1001000001_c;
  abcd[1907] = 2.0E0*I_ERI_Pz_S_Dyz_Px_C1001000001_c;
  abcd[1937] = 2.0E0*I_ERI_Px_S_Dxy_Py_C1001000001_c;
  abcd[1938] = 2.0E0*I_ERI_Py_S_Dxy_Py_C1001000001_c;
  abcd[1939] = 2.0E0*I_ERI_Pz_S_Dxy_Py_C1001000001_c;
  abcd[1953] = 2.0E0*I_ERI_Px_S_D2y_Py_C1001000001_c-1*I_ERI_Px_S_S_Py_C1001000001;
  abcd[1954] = 2.0E0*I_ERI_Py_S_D2y_Py_C1001000001_c-1*I_ERI_Py_S_S_Py_C1001000001;
  abcd[1955] = 2.0E0*I_ERI_Pz_S_D2y_Py_C1001000001_c-1*I_ERI_Pz_S_S_Py_C1001000001;
  abcd[1969] = 2.0E0*I_ERI_Px_S_Dyz_Py_C1001000001_c;
  abcd[1970] = 2.0E0*I_ERI_Py_S_Dyz_Py_C1001000001_c;
  abcd[1971] = 2.0E0*I_ERI_Pz_S_Dyz_Py_C1001000001_c;
  abcd[2001] = 2.0E0*I_ERI_Px_S_Dxy_Pz_C1001000001_c;
  abcd[2002] = 2.0E0*I_ERI_Py_S_Dxy_Pz_C1001000001_c;
  abcd[2003] = 2.0E0*I_ERI_Pz_S_Dxy_Pz_C1001000001_c;
  abcd[2017] = 2.0E0*I_ERI_Px_S_D2y_Pz_C1001000001_c-1*I_ERI_Px_S_S_Pz_C1001000001;
  abcd[2018] = 2.0E0*I_ERI_Py_S_D2y_Pz_C1001000001_c-1*I_ERI_Py_S_S_Pz_C1001000001;
  abcd[2019] = 2.0E0*I_ERI_Pz_S_D2y_Pz_C1001000001_c-1*I_ERI_Pz_S_S_Pz_C1001000001;
  abcd[2033] = 2.0E0*I_ERI_Px_S_Dyz_Pz_C1001000001_c;
  abcd[2034] = 2.0E0*I_ERI_Py_S_Dyz_Pz_C1001000001_c;
  abcd[2035] = 2.0E0*I_ERI_Pz_S_Dyz_Pz_C1001000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_P_C1001001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1001001000
   ************************************************************/
  abcd[1876] = 2.0E0*I_ERI_S_Px_Dxy_Px_C1001001000_c;
  abcd[1880] = 2.0E0*I_ERI_S_Py_Dxy_Px_C1001001000_c;
  abcd[1884] = 2.0E0*I_ERI_S_Pz_Dxy_Px_C1001001000_c;
  abcd[1892] = 2.0E0*I_ERI_S_Px_D2y_Px_C1001001000_c-1*I_ERI_S_Px_S_Px_C1001001000;
  abcd[1896] = 2.0E0*I_ERI_S_Py_D2y_Px_C1001001000_c-1*I_ERI_S_Py_S_Px_C1001001000;
  abcd[1900] = 2.0E0*I_ERI_S_Pz_D2y_Px_C1001001000_c-1*I_ERI_S_Pz_S_Px_C1001001000;
  abcd[1908] = 2.0E0*I_ERI_S_Px_Dyz_Px_C1001001000_c;
  abcd[1912] = 2.0E0*I_ERI_S_Py_Dyz_Px_C1001001000_c;
  abcd[1916] = 2.0E0*I_ERI_S_Pz_Dyz_Px_C1001001000_c;
  abcd[1940] = 2.0E0*I_ERI_S_Px_Dxy_Py_C1001001000_c;
  abcd[1944] = 2.0E0*I_ERI_S_Py_Dxy_Py_C1001001000_c;
  abcd[1948] = 2.0E0*I_ERI_S_Pz_Dxy_Py_C1001001000_c;
  abcd[1956] = 2.0E0*I_ERI_S_Px_D2y_Py_C1001001000_c-1*I_ERI_S_Px_S_Py_C1001001000;
  abcd[1960] = 2.0E0*I_ERI_S_Py_D2y_Py_C1001001000_c-1*I_ERI_S_Py_S_Py_C1001001000;
  abcd[1964] = 2.0E0*I_ERI_S_Pz_D2y_Py_C1001001000_c-1*I_ERI_S_Pz_S_Py_C1001001000;
  abcd[1972] = 2.0E0*I_ERI_S_Px_Dyz_Py_C1001001000_c;
  abcd[1976] = 2.0E0*I_ERI_S_Py_Dyz_Py_C1001001000_c;
  abcd[1980] = 2.0E0*I_ERI_S_Pz_Dyz_Py_C1001001000_c;
  abcd[2004] = 2.0E0*I_ERI_S_Px_Dxy_Pz_C1001001000_c;
  abcd[2008] = 2.0E0*I_ERI_S_Py_Dxy_Pz_C1001001000_c;
  abcd[2012] = 2.0E0*I_ERI_S_Pz_Dxy_Pz_C1001001000_c;
  abcd[2020] = 2.0E0*I_ERI_S_Px_D2y_Pz_C1001001000_c-1*I_ERI_S_Px_S_Pz_C1001001000;
  abcd[2024] = 2.0E0*I_ERI_S_Py_D2y_Pz_C1001001000_c-1*I_ERI_S_Py_S_Pz_C1001001000;
  abcd[2028] = 2.0E0*I_ERI_S_Pz_D2y_Pz_C1001001000_c-1*I_ERI_S_Pz_S_Pz_C1001001000;
  abcd[2036] = 2.0E0*I_ERI_S_Px_Dyz_Pz_C1001001000_c;
  abcd[2040] = 2.0E0*I_ERI_S_Py_Dyz_Pz_C1001001000_c;
  abcd[2044] = 2.0E0*I_ERI_S_Pz_Dyz_Pz_C1001001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001001
   ************************************************************/
  abcd[1877] = 2.0E0*I_ERI_Px_Px_Dxy_Px_C1001001001_c;
  abcd[1878] = 2.0E0*I_ERI_Py_Px_Dxy_Px_C1001001001_c;
  abcd[1879] = 2.0E0*I_ERI_Pz_Px_Dxy_Px_C1001001001_c;
  abcd[1881] = 2.0E0*I_ERI_Px_Py_Dxy_Px_C1001001001_c;
  abcd[1882] = 2.0E0*I_ERI_Py_Py_Dxy_Px_C1001001001_c;
  abcd[1883] = 2.0E0*I_ERI_Pz_Py_Dxy_Px_C1001001001_c;
  abcd[1885] = 2.0E0*I_ERI_Px_Pz_Dxy_Px_C1001001001_c;
  abcd[1886] = 2.0E0*I_ERI_Py_Pz_Dxy_Px_C1001001001_c;
  abcd[1887] = 2.0E0*I_ERI_Pz_Pz_Dxy_Px_C1001001001_c;
  abcd[1893] = 2.0E0*I_ERI_Px_Px_D2y_Px_C1001001001_c-1*I_ERI_Px_Px_S_Px_C1001001001;
  abcd[1894] = 2.0E0*I_ERI_Py_Px_D2y_Px_C1001001001_c-1*I_ERI_Py_Px_S_Px_C1001001001;
  abcd[1895] = 2.0E0*I_ERI_Pz_Px_D2y_Px_C1001001001_c-1*I_ERI_Pz_Px_S_Px_C1001001001;
  abcd[1897] = 2.0E0*I_ERI_Px_Py_D2y_Px_C1001001001_c-1*I_ERI_Px_Py_S_Px_C1001001001;
  abcd[1898] = 2.0E0*I_ERI_Py_Py_D2y_Px_C1001001001_c-1*I_ERI_Py_Py_S_Px_C1001001001;
  abcd[1899] = 2.0E0*I_ERI_Pz_Py_D2y_Px_C1001001001_c-1*I_ERI_Pz_Py_S_Px_C1001001001;
  abcd[1901] = 2.0E0*I_ERI_Px_Pz_D2y_Px_C1001001001_c-1*I_ERI_Px_Pz_S_Px_C1001001001;
  abcd[1902] = 2.0E0*I_ERI_Py_Pz_D2y_Px_C1001001001_c-1*I_ERI_Py_Pz_S_Px_C1001001001;
  abcd[1903] = 2.0E0*I_ERI_Pz_Pz_D2y_Px_C1001001001_c-1*I_ERI_Pz_Pz_S_Px_C1001001001;
  abcd[1909] = 2.0E0*I_ERI_Px_Px_Dyz_Px_C1001001001_c;
  abcd[1910] = 2.0E0*I_ERI_Py_Px_Dyz_Px_C1001001001_c;
  abcd[1911] = 2.0E0*I_ERI_Pz_Px_Dyz_Px_C1001001001_c;
  abcd[1913] = 2.0E0*I_ERI_Px_Py_Dyz_Px_C1001001001_c;
  abcd[1914] = 2.0E0*I_ERI_Py_Py_Dyz_Px_C1001001001_c;
  abcd[1915] = 2.0E0*I_ERI_Pz_Py_Dyz_Px_C1001001001_c;
  abcd[1917] = 2.0E0*I_ERI_Px_Pz_Dyz_Px_C1001001001_c;
  abcd[1918] = 2.0E0*I_ERI_Py_Pz_Dyz_Px_C1001001001_c;
  abcd[1919] = 2.0E0*I_ERI_Pz_Pz_Dyz_Px_C1001001001_c;
  abcd[1941] = 2.0E0*I_ERI_Px_Px_Dxy_Py_C1001001001_c;
  abcd[1942] = 2.0E0*I_ERI_Py_Px_Dxy_Py_C1001001001_c;
  abcd[1943] = 2.0E0*I_ERI_Pz_Px_Dxy_Py_C1001001001_c;
  abcd[1945] = 2.0E0*I_ERI_Px_Py_Dxy_Py_C1001001001_c;
  abcd[1946] = 2.0E0*I_ERI_Py_Py_Dxy_Py_C1001001001_c;
  abcd[1947] = 2.0E0*I_ERI_Pz_Py_Dxy_Py_C1001001001_c;
  abcd[1949] = 2.0E0*I_ERI_Px_Pz_Dxy_Py_C1001001001_c;
  abcd[1950] = 2.0E0*I_ERI_Py_Pz_Dxy_Py_C1001001001_c;
  abcd[1951] = 2.0E0*I_ERI_Pz_Pz_Dxy_Py_C1001001001_c;
  abcd[1957] = 2.0E0*I_ERI_Px_Px_D2y_Py_C1001001001_c-1*I_ERI_Px_Px_S_Py_C1001001001;
  abcd[1958] = 2.0E0*I_ERI_Py_Px_D2y_Py_C1001001001_c-1*I_ERI_Py_Px_S_Py_C1001001001;
  abcd[1959] = 2.0E0*I_ERI_Pz_Px_D2y_Py_C1001001001_c-1*I_ERI_Pz_Px_S_Py_C1001001001;
  abcd[1961] = 2.0E0*I_ERI_Px_Py_D2y_Py_C1001001001_c-1*I_ERI_Px_Py_S_Py_C1001001001;
  abcd[1962] = 2.0E0*I_ERI_Py_Py_D2y_Py_C1001001001_c-1*I_ERI_Py_Py_S_Py_C1001001001;
  abcd[1963] = 2.0E0*I_ERI_Pz_Py_D2y_Py_C1001001001_c-1*I_ERI_Pz_Py_S_Py_C1001001001;
  abcd[1965] = 2.0E0*I_ERI_Px_Pz_D2y_Py_C1001001001_c-1*I_ERI_Px_Pz_S_Py_C1001001001;
  abcd[1966] = 2.0E0*I_ERI_Py_Pz_D2y_Py_C1001001001_c-1*I_ERI_Py_Pz_S_Py_C1001001001;
  abcd[1967] = 2.0E0*I_ERI_Pz_Pz_D2y_Py_C1001001001_c-1*I_ERI_Pz_Pz_S_Py_C1001001001;
  abcd[1973] = 2.0E0*I_ERI_Px_Px_Dyz_Py_C1001001001_c;
  abcd[1974] = 2.0E0*I_ERI_Py_Px_Dyz_Py_C1001001001_c;
  abcd[1975] = 2.0E0*I_ERI_Pz_Px_Dyz_Py_C1001001001_c;
  abcd[1977] = 2.0E0*I_ERI_Px_Py_Dyz_Py_C1001001001_c;
  abcd[1978] = 2.0E0*I_ERI_Py_Py_Dyz_Py_C1001001001_c;
  abcd[1979] = 2.0E0*I_ERI_Pz_Py_Dyz_Py_C1001001001_c;
  abcd[1981] = 2.0E0*I_ERI_Px_Pz_Dyz_Py_C1001001001_c;
  abcd[1982] = 2.0E0*I_ERI_Py_Pz_Dyz_Py_C1001001001_c;
  abcd[1983] = 2.0E0*I_ERI_Pz_Pz_Dyz_Py_C1001001001_c;
  abcd[2005] = 2.0E0*I_ERI_Px_Px_Dxy_Pz_C1001001001_c;
  abcd[2006] = 2.0E0*I_ERI_Py_Px_Dxy_Pz_C1001001001_c;
  abcd[2007] = 2.0E0*I_ERI_Pz_Px_Dxy_Pz_C1001001001_c;
  abcd[2009] = 2.0E0*I_ERI_Px_Py_Dxy_Pz_C1001001001_c;
  abcd[2010] = 2.0E0*I_ERI_Py_Py_Dxy_Pz_C1001001001_c;
  abcd[2011] = 2.0E0*I_ERI_Pz_Py_Dxy_Pz_C1001001001_c;
  abcd[2013] = 2.0E0*I_ERI_Px_Pz_Dxy_Pz_C1001001001_c;
  abcd[2014] = 2.0E0*I_ERI_Py_Pz_Dxy_Pz_C1001001001_c;
  abcd[2015] = 2.0E0*I_ERI_Pz_Pz_Dxy_Pz_C1001001001_c;
  abcd[2021] = 2.0E0*I_ERI_Px_Px_D2y_Pz_C1001001001_c-1*I_ERI_Px_Px_S_Pz_C1001001001;
  abcd[2022] = 2.0E0*I_ERI_Py_Px_D2y_Pz_C1001001001_c-1*I_ERI_Py_Px_S_Pz_C1001001001;
  abcd[2023] = 2.0E0*I_ERI_Pz_Px_D2y_Pz_C1001001001_c-1*I_ERI_Pz_Px_S_Pz_C1001001001;
  abcd[2025] = 2.0E0*I_ERI_Px_Py_D2y_Pz_C1001001001_c-1*I_ERI_Px_Py_S_Pz_C1001001001;
  abcd[2026] = 2.0E0*I_ERI_Py_Py_D2y_Pz_C1001001001_c-1*I_ERI_Py_Py_S_Pz_C1001001001;
  abcd[2027] = 2.0E0*I_ERI_Pz_Py_D2y_Pz_C1001001001_c-1*I_ERI_Pz_Py_S_Pz_C1001001001;
  abcd[2029] = 2.0E0*I_ERI_Px_Pz_D2y_Pz_C1001001001_c-1*I_ERI_Px_Pz_S_Pz_C1001001001;
  abcd[2030] = 2.0E0*I_ERI_Py_Pz_D2y_Pz_C1001001001_c-1*I_ERI_Py_Pz_S_Pz_C1001001001;
  abcd[2031] = 2.0E0*I_ERI_Pz_Pz_D2y_Pz_C1001001001_c-1*I_ERI_Pz_Pz_S_Pz_C1001001001;
  abcd[2037] = 2.0E0*I_ERI_Px_Px_Dyz_Pz_C1001001001_c;
  abcd[2038] = 2.0E0*I_ERI_Py_Px_Dyz_Pz_C1001001001_c;
  abcd[2039] = 2.0E0*I_ERI_Pz_Px_Dyz_Pz_C1001001001_c;
  abcd[2041] = 2.0E0*I_ERI_Px_Py_Dyz_Pz_C1001001001_c;
  abcd[2042] = 2.0E0*I_ERI_Py_Py_Dyz_Pz_C1001001001_c;
  abcd[2043] = 2.0E0*I_ERI_Pz_Py_Dyz_Pz_C1001001001_c;
  abcd[2045] = 2.0E0*I_ERI_Px_Pz_Dyz_Pz_C1001001001_c;
  abcd[2046] = 2.0E0*I_ERI_Py_Pz_Dyz_Pz_C1001001001_c;
  abcd[2047] = 2.0E0*I_ERI_Pz_Pz_Dyz_Pz_C1001001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_S_C0_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C0_c
   ************************************************************/
  abcd[2048] = 2.0E0*I_ERI_S_S_Pz_S_C0_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_S_C1_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1_c
   ************************************************************/
  abcd[2049] = 2.0E0*I_ERI_Px_S_Pz_S_C1_c;
  abcd[2050] = 2.0E0*I_ERI_Py_S_Pz_S_C1_c;
  abcd[2051] = 2.0E0*I_ERI_Pz_S_Pz_S_C1_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_S_C1000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1000_c
   ************************************************************/
  abcd[2052] = 2.0E0*I_ERI_S_Px_Pz_S_C1000_c;
  abcd[2056] = 2.0E0*I_ERI_S_Py_Pz_S_C1000_c;
  abcd[2060] = 2.0E0*I_ERI_S_Pz_Pz_S_C1000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_S_C1001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001_c
   ************************************************************/
  abcd[2053] = 2.0E0*I_ERI_Px_Px_Pz_S_C1001_c;
  abcd[2054] = 2.0E0*I_ERI_Py_Px_Pz_S_C1001_c;
  abcd[2055] = 2.0E0*I_ERI_Pz_Px_Pz_S_C1001_c;
  abcd[2057] = 2.0E0*I_ERI_Px_Py_Pz_S_C1001_c;
  abcd[2058] = 2.0E0*I_ERI_Py_Py_Pz_S_C1001_c;
  abcd[2059] = 2.0E0*I_ERI_Pz_Py_Pz_S_C1001_c;
  abcd[2061] = 2.0E0*I_ERI_Px_Pz_Pz_S_C1001_c;
  abcd[2062] = 2.0E0*I_ERI_Py_Pz_Pz_S_C1001_c;
  abcd[2063] = 2.0E0*I_ERI_Pz_Pz_Pz_S_C1001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[2064] = 2.0E0*I_ERI_S_S_Dxz_S_C1000000_c;
  abcd[2080] = 2.0E0*I_ERI_S_S_Dyz_S_C1000000_c;
  abcd[2096] = 2.0E0*I_ERI_S_S_D2z_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[2065] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[2066] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[2067] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;
  abcd[2081] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[2082] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[2083] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;
  abcd[2097] = 2.0E0*I_ERI_Px_S_D2z_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[2098] = 2.0E0*I_ERI_Py_S_D2z_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[2099] = 2.0E0*I_ERI_Pz_S_D2z_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[2068] = 2.0E0*I_ERI_S_Px_Dxz_S_C1001000_c;
  abcd[2072] = 2.0E0*I_ERI_S_Py_Dxz_S_C1001000_c;
  abcd[2076] = 2.0E0*I_ERI_S_Pz_Dxz_S_C1001000_c;
  abcd[2084] = 2.0E0*I_ERI_S_Px_Dyz_S_C1001000_c;
  abcd[2088] = 2.0E0*I_ERI_S_Py_Dyz_S_C1001000_c;
  abcd[2092] = 2.0E0*I_ERI_S_Pz_Dyz_S_C1001000_c;
  abcd[2100] = 2.0E0*I_ERI_S_Px_D2z_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[2104] = 2.0E0*I_ERI_S_Py_D2z_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[2108] = 2.0E0*I_ERI_S_Pz_D2z_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[2069] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[2070] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[2071] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[2073] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[2074] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[2075] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[2077] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[2078] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[2079] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;
  abcd[2085] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[2086] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[2087] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[2089] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[2090] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[2091] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[2093] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[2094] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[2095] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;
  abcd[2101] = 2.0E0*I_ERI_Px_Px_D2z_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[2102] = 2.0E0*I_ERI_Py_Px_D2z_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[2103] = 2.0E0*I_ERI_Pz_Px_D2z_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[2105] = 2.0E0*I_ERI_Px_Py_D2z_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[2106] = 2.0E0*I_ERI_Py_Py_D2z_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[2107] = 2.0E0*I_ERI_Pz_Py_D2z_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[2109] = 2.0E0*I_ERI_Px_Pz_D2z_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[2110] = 2.0E0*I_ERI_Py_Pz_D2z_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[2111] = 2.0E0*I_ERI_Pz_Pz_D2z_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_S_P_C1000000000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000000_c
   ************************************************************/
  abcd[2112] = 2.0E0*I_ERI_S_S_Pz_Px_C1000000000_c;
  abcd[2176] = 2.0E0*I_ERI_S_S_Pz_Py_C1000000000_c;
  abcd[2240] = 2.0E0*I_ERI_S_S_Pz_Pz_C1000000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_S_P_C1000000001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000000001_c
   ************************************************************/
  abcd[2113] = 2.0E0*I_ERI_Px_S_Pz_Px_C1000000001_c;
  abcd[2114] = 2.0E0*I_ERI_Py_S_Pz_Px_C1000000001_c;
  abcd[2115] = 2.0E0*I_ERI_Pz_S_Pz_Px_C1000000001_c;
  abcd[2177] = 2.0E0*I_ERI_Px_S_Pz_Py_C1000000001_c;
  abcd[2178] = 2.0E0*I_ERI_Py_S_Pz_Py_C1000000001_c;
  abcd[2179] = 2.0E0*I_ERI_Pz_S_Pz_Py_C1000000001_c;
  abcd[2241] = 2.0E0*I_ERI_Px_S_Pz_Pz_C1000000001_c;
  abcd[2242] = 2.0E0*I_ERI_Py_S_Pz_Pz_C1000000001_c;
  abcd[2243] = 2.0E0*I_ERI_Pz_S_Pz_Pz_C1000000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_S_P_C1000001000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1000001000_c
   ************************************************************/
  abcd[2116] = 2.0E0*I_ERI_S_Px_Pz_Px_C1000001000_c;
  abcd[2120] = 2.0E0*I_ERI_S_Py_Pz_Px_C1000001000_c;
  abcd[2124] = 2.0E0*I_ERI_S_Pz_Pz_Px_C1000001000_c;
  abcd[2180] = 2.0E0*I_ERI_S_Px_Pz_Py_C1000001000_c;
  abcd[2184] = 2.0E0*I_ERI_S_Py_Pz_Py_C1000001000_c;
  abcd[2188] = 2.0E0*I_ERI_S_Pz_Pz_Py_C1000001000_c;
  abcd[2244] = 2.0E0*I_ERI_S_Px_Pz_Pz_C1000001000_c;
  abcd[2248] = 2.0E0*I_ERI_S_Py_Pz_Pz_C1000001000_c;
  abcd[2252] = 2.0E0*I_ERI_S_Pz_Pz_Pz_C1000001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_S_P_C1000001001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1000001001_c
   ************************************************************/
  abcd[2117] = 2.0E0*I_ERI_Px_Px_Pz_Px_C1000001001_c;
  abcd[2118] = 2.0E0*I_ERI_Py_Px_Pz_Px_C1000001001_c;
  abcd[2119] = 2.0E0*I_ERI_Pz_Px_Pz_Px_C1000001001_c;
  abcd[2121] = 2.0E0*I_ERI_Px_Py_Pz_Px_C1000001001_c;
  abcd[2122] = 2.0E0*I_ERI_Py_Py_Pz_Px_C1000001001_c;
  abcd[2123] = 2.0E0*I_ERI_Pz_Py_Pz_Px_C1000001001_c;
  abcd[2125] = 2.0E0*I_ERI_Px_Pz_Pz_Px_C1000001001_c;
  abcd[2126] = 2.0E0*I_ERI_Py_Pz_Pz_Px_C1000001001_c;
  abcd[2127] = 2.0E0*I_ERI_Pz_Pz_Pz_Px_C1000001001_c;
  abcd[2181] = 2.0E0*I_ERI_Px_Px_Pz_Py_C1000001001_c;
  abcd[2182] = 2.0E0*I_ERI_Py_Px_Pz_Py_C1000001001_c;
  abcd[2183] = 2.0E0*I_ERI_Pz_Px_Pz_Py_C1000001001_c;
  abcd[2185] = 2.0E0*I_ERI_Px_Py_Pz_Py_C1000001001_c;
  abcd[2186] = 2.0E0*I_ERI_Py_Py_Pz_Py_C1000001001_c;
  abcd[2187] = 2.0E0*I_ERI_Pz_Py_Pz_Py_C1000001001_c;
  abcd[2189] = 2.0E0*I_ERI_Px_Pz_Pz_Py_C1000001001_c;
  abcd[2190] = 2.0E0*I_ERI_Py_Pz_Pz_Py_C1000001001_c;
  abcd[2191] = 2.0E0*I_ERI_Pz_Pz_Pz_Py_C1000001001_c;
  abcd[2245] = 2.0E0*I_ERI_Px_Px_Pz_Pz_C1000001001_c;
  abcd[2246] = 2.0E0*I_ERI_Py_Px_Pz_Pz_C1000001001_c;
  abcd[2247] = 2.0E0*I_ERI_Pz_Px_Pz_Pz_C1000001001_c;
  abcd[2249] = 2.0E0*I_ERI_Px_Py_Pz_Pz_C1000001001_c;
  abcd[2250] = 2.0E0*I_ERI_Py_Py_Pz_Pz_C1000001001_c;
  abcd[2251] = 2.0E0*I_ERI_Pz_Py_Pz_Pz_C1000001001_c;
  abcd[2253] = 2.0E0*I_ERI_Px_Pz_Pz_Pz_C1000001001_c;
  abcd[2254] = 2.0E0*I_ERI_Py_Pz_Pz_Pz_C1000001001_c;
  abcd[2255] = 2.0E0*I_ERI_Pz_Pz_Pz_Pz_C1000001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_P_C1001000000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_P_C1001000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_P_C1001000000
   ************************************************************/
  abcd[2128] = 2.0E0*I_ERI_S_S_Dxz_Px_C1001000000_c;
  abcd[2144] = 2.0E0*I_ERI_S_S_Dyz_Px_C1001000000_c;
  abcd[2160] = 2.0E0*I_ERI_S_S_D2z_Px_C1001000000_c-1*I_ERI_S_S_S_Px_C1001000000;
  abcd[2192] = 2.0E0*I_ERI_S_S_Dxz_Py_C1001000000_c;
  abcd[2208] = 2.0E0*I_ERI_S_S_Dyz_Py_C1001000000_c;
  abcd[2224] = 2.0E0*I_ERI_S_S_D2z_Py_C1001000000_c-1*I_ERI_S_S_S_Py_C1001000000;
  abcd[2256] = 2.0E0*I_ERI_S_S_Dxz_Pz_C1001000000_c;
  abcd[2272] = 2.0E0*I_ERI_S_S_Dyz_Pz_C1001000000_c;
  abcd[2288] = 2.0E0*I_ERI_S_S_D2z_Pz_C1001000000_c-1*I_ERI_S_S_S_Pz_C1001000000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001000001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_P_C1001000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_P_C1001000001
   ************************************************************/
  abcd[2129] = 2.0E0*I_ERI_Px_S_Dxz_Px_C1001000001_c;
  abcd[2130] = 2.0E0*I_ERI_Py_S_Dxz_Px_C1001000001_c;
  abcd[2131] = 2.0E0*I_ERI_Pz_S_Dxz_Px_C1001000001_c;
  abcd[2145] = 2.0E0*I_ERI_Px_S_Dyz_Px_C1001000001_c;
  abcd[2146] = 2.0E0*I_ERI_Py_S_Dyz_Px_C1001000001_c;
  abcd[2147] = 2.0E0*I_ERI_Pz_S_Dyz_Px_C1001000001_c;
  abcd[2161] = 2.0E0*I_ERI_Px_S_D2z_Px_C1001000001_c-1*I_ERI_Px_S_S_Px_C1001000001;
  abcd[2162] = 2.0E0*I_ERI_Py_S_D2z_Px_C1001000001_c-1*I_ERI_Py_S_S_Px_C1001000001;
  abcd[2163] = 2.0E0*I_ERI_Pz_S_D2z_Px_C1001000001_c-1*I_ERI_Pz_S_S_Px_C1001000001;
  abcd[2193] = 2.0E0*I_ERI_Px_S_Dxz_Py_C1001000001_c;
  abcd[2194] = 2.0E0*I_ERI_Py_S_Dxz_Py_C1001000001_c;
  abcd[2195] = 2.0E0*I_ERI_Pz_S_Dxz_Py_C1001000001_c;
  abcd[2209] = 2.0E0*I_ERI_Px_S_Dyz_Py_C1001000001_c;
  abcd[2210] = 2.0E0*I_ERI_Py_S_Dyz_Py_C1001000001_c;
  abcd[2211] = 2.0E0*I_ERI_Pz_S_Dyz_Py_C1001000001_c;
  abcd[2225] = 2.0E0*I_ERI_Px_S_D2z_Py_C1001000001_c-1*I_ERI_Px_S_S_Py_C1001000001;
  abcd[2226] = 2.0E0*I_ERI_Py_S_D2z_Py_C1001000001_c-1*I_ERI_Py_S_S_Py_C1001000001;
  abcd[2227] = 2.0E0*I_ERI_Pz_S_D2z_Py_C1001000001_c-1*I_ERI_Pz_S_S_Py_C1001000001;
  abcd[2257] = 2.0E0*I_ERI_Px_S_Dxz_Pz_C1001000001_c;
  abcd[2258] = 2.0E0*I_ERI_Py_S_Dxz_Pz_C1001000001_c;
  abcd[2259] = 2.0E0*I_ERI_Pz_S_Dxz_Pz_C1001000001_c;
  abcd[2273] = 2.0E0*I_ERI_Px_S_Dyz_Pz_C1001000001_c;
  abcd[2274] = 2.0E0*I_ERI_Py_S_Dyz_Pz_C1001000001_c;
  abcd[2275] = 2.0E0*I_ERI_Pz_S_Dyz_Pz_C1001000001_c;
  abcd[2289] = 2.0E0*I_ERI_Px_S_D2z_Pz_C1001000001_c-1*I_ERI_Px_S_S_Pz_C1001000001;
  abcd[2290] = 2.0E0*I_ERI_Py_S_D2z_Pz_C1001000001_c-1*I_ERI_Py_S_S_Pz_C1001000001;
  abcd[2291] = 2.0E0*I_ERI_Pz_S_D2z_Pz_C1001000001_c-1*I_ERI_Pz_S_S_Pz_C1001000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001001000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_P_C1001001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_P_C1001001000
   ************************************************************/
  abcd[2132] = 2.0E0*I_ERI_S_Px_Dxz_Px_C1001001000_c;
  abcd[2136] = 2.0E0*I_ERI_S_Py_Dxz_Px_C1001001000_c;
  abcd[2140] = 2.0E0*I_ERI_S_Pz_Dxz_Px_C1001001000_c;
  abcd[2148] = 2.0E0*I_ERI_S_Px_Dyz_Px_C1001001000_c;
  abcd[2152] = 2.0E0*I_ERI_S_Py_Dyz_Px_C1001001000_c;
  abcd[2156] = 2.0E0*I_ERI_S_Pz_Dyz_Px_C1001001000_c;
  abcd[2164] = 2.0E0*I_ERI_S_Px_D2z_Px_C1001001000_c-1*I_ERI_S_Px_S_Px_C1001001000;
  abcd[2168] = 2.0E0*I_ERI_S_Py_D2z_Px_C1001001000_c-1*I_ERI_S_Py_S_Px_C1001001000;
  abcd[2172] = 2.0E0*I_ERI_S_Pz_D2z_Px_C1001001000_c-1*I_ERI_S_Pz_S_Px_C1001001000;
  abcd[2196] = 2.0E0*I_ERI_S_Px_Dxz_Py_C1001001000_c;
  abcd[2200] = 2.0E0*I_ERI_S_Py_Dxz_Py_C1001001000_c;
  abcd[2204] = 2.0E0*I_ERI_S_Pz_Dxz_Py_C1001001000_c;
  abcd[2212] = 2.0E0*I_ERI_S_Px_Dyz_Py_C1001001000_c;
  abcd[2216] = 2.0E0*I_ERI_S_Py_Dyz_Py_C1001001000_c;
  abcd[2220] = 2.0E0*I_ERI_S_Pz_Dyz_Py_C1001001000_c;
  abcd[2228] = 2.0E0*I_ERI_S_Px_D2z_Py_C1001001000_c-1*I_ERI_S_Px_S_Py_C1001001000;
  abcd[2232] = 2.0E0*I_ERI_S_Py_D2z_Py_C1001001000_c-1*I_ERI_S_Py_S_Py_C1001001000;
  abcd[2236] = 2.0E0*I_ERI_S_Pz_D2z_Py_C1001001000_c-1*I_ERI_S_Pz_S_Py_C1001001000;
  abcd[2260] = 2.0E0*I_ERI_S_Px_Dxz_Pz_C1001001000_c;
  abcd[2264] = 2.0E0*I_ERI_S_Py_Dxz_Pz_C1001001000_c;
  abcd[2268] = 2.0E0*I_ERI_S_Pz_Dxz_Pz_C1001001000_c;
  abcd[2276] = 2.0E0*I_ERI_S_Px_Dyz_Pz_C1001001000_c;
  abcd[2280] = 2.0E0*I_ERI_S_Py_Dyz_Pz_C1001001000_c;
  abcd[2284] = 2.0E0*I_ERI_S_Pz_Dyz_Pz_C1001001000_c;
  abcd[2292] = 2.0E0*I_ERI_S_Px_D2z_Pz_C1001001000_c-1*I_ERI_S_Px_S_Pz_C1001001000;
  abcd[2296] = 2.0E0*I_ERI_S_Py_D2z_Pz_C1001001000_c-1*I_ERI_S_Py_S_Pz_C1001001000;
  abcd[2300] = 2.0E0*I_ERI_S_Pz_D2z_Pz_C1001001000_c-1*I_ERI_S_Pz_S_Pz_C1001001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_C1001001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P_C1001001001
   ************************************************************/
  abcd[2133] = 2.0E0*I_ERI_Px_Px_Dxz_Px_C1001001001_c;
  abcd[2134] = 2.0E0*I_ERI_Py_Px_Dxz_Px_C1001001001_c;
  abcd[2135] = 2.0E0*I_ERI_Pz_Px_Dxz_Px_C1001001001_c;
  abcd[2137] = 2.0E0*I_ERI_Px_Py_Dxz_Px_C1001001001_c;
  abcd[2138] = 2.0E0*I_ERI_Py_Py_Dxz_Px_C1001001001_c;
  abcd[2139] = 2.0E0*I_ERI_Pz_Py_Dxz_Px_C1001001001_c;
  abcd[2141] = 2.0E0*I_ERI_Px_Pz_Dxz_Px_C1001001001_c;
  abcd[2142] = 2.0E0*I_ERI_Py_Pz_Dxz_Px_C1001001001_c;
  abcd[2143] = 2.0E0*I_ERI_Pz_Pz_Dxz_Px_C1001001001_c;
  abcd[2149] = 2.0E0*I_ERI_Px_Px_Dyz_Px_C1001001001_c;
  abcd[2150] = 2.0E0*I_ERI_Py_Px_Dyz_Px_C1001001001_c;
  abcd[2151] = 2.0E0*I_ERI_Pz_Px_Dyz_Px_C1001001001_c;
  abcd[2153] = 2.0E0*I_ERI_Px_Py_Dyz_Px_C1001001001_c;
  abcd[2154] = 2.0E0*I_ERI_Py_Py_Dyz_Px_C1001001001_c;
  abcd[2155] = 2.0E0*I_ERI_Pz_Py_Dyz_Px_C1001001001_c;
  abcd[2157] = 2.0E0*I_ERI_Px_Pz_Dyz_Px_C1001001001_c;
  abcd[2158] = 2.0E0*I_ERI_Py_Pz_Dyz_Px_C1001001001_c;
  abcd[2159] = 2.0E0*I_ERI_Pz_Pz_Dyz_Px_C1001001001_c;
  abcd[2165] = 2.0E0*I_ERI_Px_Px_D2z_Px_C1001001001_c-1*I_ERI_Px_Px_S_Px_C1001001001;
  abcd[2166] = 2.0E0*I_ERI_Py_Px_D2z_Px_C1001001001_c-1*I_ERI_Py_Px_S_Px_C1001001001;
  abcd[2167] = 2.0E0*I_ERI_Pz_Px_D2z_Px_C1001001001_c-1*I_ERI_Pz_Px_S_Px_C1001001001;
  abcd[2169] = 2.0E0*I_ERI_Px_Py_D2z_Px_C1001001001_c-1*I_ERI_Px_Py_S_Px_C1001001001;
  abcd[2170] = 2.0E0*I_ERI_Py_Py_D2z_Px_C1001001001_c-1*I_ERI_Py_Py_S_Px_C1001001001;
  abcd[2171] = 2.0E0*I_ERI_Pz_Py_D2z_Px_C1001001001_c-1*I_ERI_Pz_Py_S_Px_C1001001001;
  abcd[2173] = 2.0E0*I_ERI_Px_Pz_D2z_Px_C1001001001_c-1*I_ERI_Px_Pz_S_Px_C1001001001;
  abcd[2174] = 2.0E0*I_ERI_Py_Pz_D2z_Px_C1001001001_c-1*I_ERI_Py_Pz_S_Px_C1001001001;
  abcd[2175] = 2.0E0*I_ERI_Pz_Pz_D2z_Px_C1001001001_c-1*I_ERI_Pz_Pz_S_Px_C1001001001;
  abcd[2197] = 2.0E0*I_ERI_Px_Px_Dxz_Py_C1001001001_c;
  abcd[2198] = 2.0E0*I_ERI_Py_Px_Dxz_Py_C1001001001_c;
  abcd[2199] = 2.0E0*I_ERI_Pz_Px_Dxz_Py_C1001001001_c;
  abcd[2201] = 2.0E0*I_ERI_Px_Py_Dxz_Py_C1001001001_c;
  abcd[2202] = 2.0E0*I_ERI_Py_Py_Dxz_Py_C1001001001_c;
  abcd[2203] = 2.0E0*I_ERI_Pz_Py_Dxz_Py_C1001001001_c;
  abcd[2205] = 2.0E0*I_ERI_Px_Pz_Dxz_Py_C1001001001_c;
  abcd[2206] = 2.0E0*I_ERI_Py_Pz_Dxz_Py_C1001001001_c;
  abcd[2207] = 2.0E0*I_ERI_Pz_Pz_Dxz_Py_C1001001001_c;
  abcd[2213] = 2.0E0*I_ERI_Px_Px_Dyz_Py_C1001001001_c;
  abcd[2214] = 2.0E0*I_ERI_Py_Px_Dyz_Py_C1001001001_c;
  abcd[2215] = 2.0E0*I_ERI_Pz_Px_Dyz_Py_C1001001001_c;
  abcd[2217] = 2.0E0*I_ERI_Px_Py_Dyz_Py_C1001001001_c;
  abcd[2218] = 2.0E0*I_ERI_Py_Py_Dyz_Py_C1001001001_c;
  abcd[2219] = 2.0E0*I_ERI_Pz_Py_Dyz_Py_C1001001001_c;
  abcd[2221] = 2.0E0*I_ERI_Px_Pz_Dyz_Py_C1001001001_c;
  abcd[2222] = 2.0E0*I_ERI_Py_Pz_Dyz_Py_C1001001001_c;
  abcd[2223] = 2.0E0*I_ERI_Pz_Pz_Dyz_Py_C1001001001_c;
  abcd[2229] = 2.0E0*I_ERI_Px_Px_D2z_Py_C1001001001_c-1*I_ERI_Px_Px_S_Py_C1001001001;
  abcd[2230] = 2.0E0*I_ERI_Py_Px_D2z_Py_C1001001001_c-1*I_ERI_Py_Px_S_Py_C1001001001;
  abcd[2231] = 2.0E0*I_ERI_Pz_Px_D2z_Py_C1001001001_c-1*I_ERI_Pz_Px_S_Py_C1001001001;
  abcd[2233] = 2.0E0*I_ERI_Px_Py_D2z_Py_C1001001001_c-1*I_ERI_Px_Py_S_Py_C1001001001;
  abcd[2234] = 2.0E0*I_ERI_Py_Py_D2z_Py_C1001001001_c-1*I_ERI_Py_Py_S_Py_C1001001001;
  abcd[2235] = 2.0E0*I_ERI_Pz_Py_D2z_Py_C1001001001_c-1*I_ERI_Pz_Py_S_Py_C1001001001;
  abcd[2237] = 2.0E0*I_ERI_Px_Pz_D2z_Py_C1001001001_c-1*I_ERI_Px_Pz_S_Py_C1001001001;
  abcd[2238] = 2.0E0*I_ERI_Py_Pz_D2z_Py_C1001001001_c-1*I_ERI_Py_Pz_S_Py_C1001001001;
  abcd[2239] = 2.0E0*I_ERI_Pz_Pz_D2z_Py_C1001001001_c-1*I_ERI_Pz_Pz_S_Py_C1001001001;
  abcd[2261] = 2.0E0*I_ERI_Px_Px_Dxz_Pz_C1001001001_c;
  abcd[2262] = 2.0E0*I_ERI_Py_Px_Dxz_Pz_C1001001001_c;
  abcd[2263] = 2.0E0*I_ERI_Pz_Px_Dxz_Pz_C1001001001_c;
  abcd[2265] = 2.0E0*I_ERI_Px_Py_Dxz_Pz_C1001001001_c;
  abcd[2266] = 2.0E0*I_ERI_Py_Py_Dxz_Pz_C1001001001_c;
  abcd[2267] = 2.0E0*I_ERI_Pz_Py_Dxz_Pz_C1001001001_c;
  abcd[2269] = 2.0E0*I_ERI_Px_Pz_Dxz_Pz_C1001001001_c;
  abcd[2270] = 2.0E0*I_ERI_Py_Pz_Dxz_Pz_C1001001001_c;
  abcd[2271] = 2.0E0*I_ERI_Pz_Pz_Dxz_Pz_C1001001001_c;
  abcd[2277] = 2.0E0*I_ERI_Px_Px_Dyz_Pz_C1001001001_c;
  abcd[2278] = 2.0E0*I_ERI_Py_Px_Dyz_Pz_C1001001001_c;
  abcd[2279] = 2.0E0*I_ERI_Pz_Px_Dyz_Pz_C1001001001_c;
  abcd[2281] = 2.0E0*I_ERI_Px_Py_Dyz_Pz_C1001001001_c;
  abcd[2282] = 2.0E0*I_ERI_Py_Py_Dyz_Pz_C1001001001_c;
  abcd[2283] = 2.0E0*I_ERI_Pz_Py_Dyz_Pz_C1001001001_c;
  abcd[2285] = 2.0E0*I_ERI_Px_Pz_Dyz_Pz_C1001001001_c;
  abcd[2286] = 2.0E0*I_ERI_Py_Pz_Dyz_Pz_C1001001001_c;
  abcd[2287] = 2.0E0*I_ERI_Pz_Pz_Dyz_Pz_C1001001001_c;
  abcd[2293] = 2.0E0*I_ERI_Px_Px_D2z_Pz_C1001001001_c-1*I_ERI_Px_Px_S_Pz_C1001001001;
  abcd[2294] = 2.0E0*I_ERI_Py_Px_D2z_Pz_C1001001001_c-1*I_ERI_Py_Px_S_Pz_C1001001001;
  abcd[2295] = 2.0E0*I_ERI_Pz_Px_D2z_Pz_C1001001001_c-1*I_ERI_Pz_Px_S_Pz_C1001001001;
  abcd[2297] = 2.0E0*I_ERI_Px_Py_D2z_Pz_C1001001001_c-1*I_ERI_Px_Py_S_Pz_C1001001001;
  abcd[2298] = 2.0E0*I_ERI_Py_Py_D2z_Pz_C1001001001_c-1*I_ERI_Py_Py_S_Pz_C1001001001;
  abcd[2299] = 2.0E0*I_ERI_Pz_Py_D2z_Pz_C1001001001_c-1*I_ERI_Pz_Py_S_Pz_C1001001001;
  abcd[2301] = 2.0E0*I_ERI_Px_Pz_D2z_Pz_C1001001001_c-1*I_ERI_Px_Pz_S_Pz_C1001001001;
  abcd[2302] = 2.0E0*I_ERI_Py_Pz_D2z_Pz_C1001001001_c-1*I_ERI_Py_Pz_S_Pz_C1001001001;
  abcd[2303] = 2.0E0*I_ERI_Pz_Pz_D2z_Pz_C1001001001_c-1*I_ERI_Pz_Pz_S_Pz_C1001001001;
}
