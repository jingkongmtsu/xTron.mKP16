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
// BRA1 as redundant position, total RHS integrals evaluated as: 50664
// BRA2 as redundant position, total RHS integrals evaluated as: 50664
// KET1 as redundant position, total RHS integrals evaluated as: 57387
// KET2 as redundant position, total RHS integrals evaluated as: 58569
// the redundant position is: BRA1
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####
// KET1
// X
// Y
// Z
// ####
// KET2
// X
// Y
// Z
// ####

void hgp_os_eri_sp_sp_p_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
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
  Double I_ERI_S_S_D2x_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Dxy_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Dxz_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_D2y_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Dyz_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_D2z_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Px_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Py_S_C1000000_d = 0.0E0;
  Double I_ERI_S_S_Pz_S_C1000000_d = 0.0E0;
  Double I_ERI_Px_S_D2x_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000001_d = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000001_d = 0.0E0;
  Double I_ERI_S_Px_D2x_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_D2x_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_D2x_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Dxy_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Dxy_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Dxz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Dxz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_D2y_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_D2y_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_D2y_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Dyz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Dyz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_D2z_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_D2z_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_D2z_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Px_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Px_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Px_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Py_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Py_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Py_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Px_Pz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Py_Pz_S_C1001000_d = 0.0E0;
  Double I_ERI_S_Pz_Pz_S_C1001000_d = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001001_d = 0.0E0;
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
  Double I_ERI_Px_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001001_d = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001001_d = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001001_d = 0.0E0;

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
       * shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000000_a_coefs = ic2*jc2*alpha;
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
      Double SQ_ERI_D_S_P_S_C1000001_a_coefs = ic2_1*jc2*alpha;
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
      Double SQ_ERI_S_S_P_S_C1000001_coefs = ic2_1*jc2;
      I_ERI_S_S_Px_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000001 += SQ_ERI_S_S_P_S_C1000001_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001001
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001001_coefs = ic2_3*jc2;
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
       * shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000000_c_coefs = ic2*jc2*gamma;
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
      Double SQ_ERI_S_S_S_S_C1000000_coefs = ic2*jc2;
      I_ERI_S_S_S_S_C1000000 += SQ_ERI_S_S_S_S_C1000000_coefs*I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_c_coefs = ic2_1*jc2*gamma;
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
      Double SQ_ERI_P_S_S_S_C1000001_coefs = ic2_1*jc2;
      I_ERI_Px_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1000001 += SQ_ERI_P_S_S_S_C1000001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001000_c_coefs = ic2_2*jc2*gamma;
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
      Double SQ_ERI_S_P_S_S_C1001000_coefs = ic2_2*jc2;
      I_ERI_S_Px_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Px_S_S_vrr;
      I_ERI_S_Py_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Py_S_S_vrr;
      I_ERI_S_Pz_S_S_C1001000 += SQ_ERI_S_P_S_S_C1001000_coefs*I_ERI_S_Pz_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001000_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001000_a_coefs = ic2_2*jc2*alpha;
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
      Double SQ_ERI_P_S_P_S_C1001000_a_coefs = ic2_2*jc2*alpha;
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
      Double SQ_ERI_F_S_P_S_C1001001_a_coefs = ic2_3*jc2*alpha;
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
      Double SQ_ERI_D_S_P_S_C1001001_a_coefs = ic2_3*jc2*alpha;
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
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_c_coefs = ic2_3*jc2*gamma;
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
      Double SQ_ERI_P_S_D_S_C1001001_c_coefs = ic2_3*jc2*gamma;
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
      Double SQ_ERI_D_S_S_S_C1001001_coefs = ic2_3*jc2;
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
      Double SQ_ERI_P_S_S_S_C1001001_coefs = ic2_3*jc2;
      I_ERI_Px_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C1001001 += SQ_ERI_P_S_S_S_C1001001_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_C1000000_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_D_S_C1000000_d_coefs = ic2*jc2*delta;
      I_ERI_S_S_D2x_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_D2x_S_vrr;
      I_ERI_S_S_Dxy_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_Dxy_S_vrr;
      I_ERI_S_S_Dxz_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_Dxz_S_vrr;
      I_ERI_S_S_D2y_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_D2y_S_vrr;
      I_ERI_S_S_Dyz_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_Dyz_S_vrr;
      I_ERI_S_S_D2z_S_C1000000_d += SQ_ERI_S_S_D_S_C1000000_d_coefs*I_ERI_S_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_C1000000_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_S_P_S_C1000000_d_coefs = ic2*jc2*delta;
      I_ERI_S_S_Px_S_C1000000_d += SQ_ERI_S_S_P_S_C1000000_d_coefs*I_ERI_S_S_Px_S_vrr;
      I_ERI_S_S_Py_S_C1000000_d += SQ_ERI_S_S_P_S_C1000000_d_coefs*I_ERI_S_S_Py_S_vrr;
      I_ERI_S_S_Pz_S_C1000000_d += SQ_ERI_S_S_P_S_C1000000_d_coefs*I_ERI_S_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1000001_d_coefs = ic2_1*jc2*delta;
      I_ERI_Px_S_D2x_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1000001_d += SQ_ERI_P_S_D_S_C1000001_d_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000001_d_coefs = ic2_1*jc2*delta;
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
       * shell quartet name: SQ_ERI_S_P_D_S_C1001000_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_D_S_C1001000_d_coefs = ic2_2*jc2*delta;
      I_ERI_S_Px_D2x_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S_C1001000_d += SQ_ERI_S_P_D_S_C1001000_d_coefs*I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S_C1001000_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_S_P_P_S_C1001000_d_coefs = ic2_2*jc2*delta;
      I_ERI_S_Px_Px_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S_C1001000_d += SQ_ERI_S_P_P_S_C1001000_d_coefs*I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001001_d_coefs = ic2_3*jc2*delta;
      I_ERI_D2x_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001001_d += SQ_ERI_D_S_D_S_C1001001_d_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001001_d_coefs = ic2_3*jc2*delta;
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
       * shell quartet name: SQ_ERI_P_S_D_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_C1001001_d_coefs = ic2_3*jc2*delta;
      I_ERI_Px_S_D2x_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_C1001001_d += SQ_ERI_P_S_D_S_C1001001_d_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001001_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001001_d_coefs = ic2_3*jc2*delta;
      I_ERI_Px_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001001_d += SQ_ERI_P_S_P_S_C1001001_d_coefs*I_ERI_Pz_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_S_S_P_P_C1000000_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_d
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000000_d
   ************************************************************/
  Double I_ERI_S_S_Px_Px_C1000000_d = I_ERI_S_S_D2x_S_C1000000_d+CDX*I_ERI_S_S_Px_S_C1000000_d;
  Double I_ERI_S_S_Py_Px_C1000000_d = I_ERI_S_S_Dxy_S_C1000000_d+CDX*I_ERI_S_S_Py_S_C1000000_d;
  Double I_ERI_S_S_Pz_Px_C1000000_d = I_ERI_S_S_Dxz_S_C1000000_d+CDX*I_ERI_S_S_Pz_S_C1000000_d;
  Double I_ERI_S_S_Px_Py_C1000000_d = I_ERI_S_S_Dxy_S_C1000000_d+CDY*I_ERI_S_S_Px_S_C1000000_d;
  Double I_ERI_S_S_Py_Py_C1000000_d = I_ERI_S_S_D2y_S_C1000000_d+CDY*I_ERI_S_S_Py_S_C1000000_d;
  Double I_ERI_S_S_Pz_Py_C1000000_d = I_ERI_S_S_Dyz_S_C1000000_d+CDY*I_ERI_S_S_Pz_S_C1000000_d;
  Double I_ERI_S_S_Px_Pz_C1000000_d = I_ERI_S_S_Dxz_S_C1000000_d+CDZ*I_ERI_S_S_Px_S_C1000000_d;
  Double I_ERI_S_S_Py_Pz_C1000000_d = I_ERI_S_S_Dyz_S_C1000000_d+CDZ*I_ERI_S_S_Py_S_C1000000_d;
  Double I_ERI_S_S_Pz_Pz_C1000000_d = I_ERI_S_S_D2z_S_C1000000_d+CDZ*I_ERI_S_S_Pz_S_C1000000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1000001_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_d
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000001_d
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1000001_d = I_ERI_Px_S_D2x_S_C1000001_d+CDX*I_ERI_Px_S_Px_S_C1000001_d;
  Double I_ERI_Py_S_Px_Px_C1000001_d = I_ERI_Py_S_D2x_S_C1000001_d+CDX*I_ERI_Py_S_Px_S_C1000001_d;
  Double I_ERI_Pz_S_Px_Px_C1000001_d = I_ERI_Pz_S_D2x_S_C1000001_d+CDX*I_ERI_Pz_S_Px_S_C1000001_d;
  Double I_ERI_Px_S_Py_Px_C1000001_d = I_ERI_Px_S_Dxy_S_C1000001_d+CDX*I_ERI_Px_S_Py_S_C1000001_d;
  Double I_ERI_Py_S_Py_Px_C1000001_d = I_ERI_Py_S_Dxy_S_C1000001_d+CDX*I_ERI_Py_S_Py_S_C1000001_d;
  Double I_ERI_Pz_S_Py_Px_C1000001_d = I_ERI_Pz_S_Dxy_S_C1000001_d+CDX*I_ERI_Pz_S_Py_S_C1000001_d;
  Double I_ERI_Px_S_Pz_Px_C1000001_d = I_ERI_Px_S_Dxz_S_C1000001_d+CDX*I_ERI_Px_S_Pz_S_C1000001_d;
  Double I_ERI_Py_S_Pz_Px_C1000001_d = I_ERI_Py_S_Dxz_S_C1000001_d+CDX*I_ERI_Py_S_Pz_S_C1000001_d;
  Double I_ERI_Pz_S_Pz_Px_C1000001_d = I_ERI_Pz_S_Dxz_S_C1000001_d+CDX*I_ERI_Pz_S_Pz_S_C1000001_d;
  Double I_ERI_Px_S_Px_Py_C1000001_d = I_ERI_Px_S_Dxy_S_C1000001_d+CDY*I_ERI_Px_S_Px_S_C1000001_d;
  Double I_ERI_Py_S_Px_Py_C1000001_d = I_ERI_Py_S_Dxy_S_C1000001_d+CDY*I_ERI_Py_S_Px_S_C1000001_d;
  Double I_ERI_Pz_S_Px_Py_C1000001_d = I_ERI_Pz_S_Dxy_S_C1000001_d+CDY*I_ERI_Pz_S_Px_S_C1000001_d;
  Double I_ERI_Px_S_Py_Py_C1000001_d = I_ERI_Px_S_D2y_S_C1000001_d+CDY*I_ERI_Px_S_Py_S_C1000001_d;
  Double I_ERI_Py_S_Py_Py_C1000001_d = I_ERI_Py_S_D2y_S_C1000001_d+CDY*I_ERI_Py_S_Py_S_C1000001_d;
  Double I_ERI_Pz_S_Py_Py_C1000001_d = I_ERI_Pz_S_D2y_S_C1000001_d+CDY*I_ERI_Pz_S_Py_S_C1000001_d;
  Double I_ERI_Px_S_Pz_Py_C1000001_d = I_ERI_Px_S_Dyz_S_C1000001_d+CDY*I_ERI_Px_S_Pz_S_C1000001_d;
  Double I_ERI_Py_S_Pz_Py_C1000001_d = I_ERI_Py_S_Dyz_S_C1000001_d+CDY*I_ERI_Py_S_Pz_S_C1000001_d;
  Double I_ERI_Pz_S_Pz_Py_C1000001_d = I_ERI_Pz_S_Dyz_S_C1000001_d+CDY*I_ERI_Pz_S_Pz_S_C1000001_d;
  Double I_ERI_Px_S_Px_Pz_C1000001_d = I_ERI_Px_S_Dxz_S_C1000001_d+CDZ*I_ERI_Px_S_Px_S_C1000001_d;
  Double I_ERI_Py_S_Px_Pz_C1000001_d = I_ERI_Py_S_Dxz_S_C1000001_d+CDZ*I_ERI_Py_S_Px_S_C1000001_d;
  Double I_ERI_Pz_S_Px_Pz_C1000001_d = I_ERI_Pz_S_Dxz_S_C1000001_d+CDZ*I_ERI_Pz_S_Px_S_C1000001_d;
  Double I_ERI_Px_S_Py_Pz_C1000001_d = I_ERI_Px_S_Dyz_S_C1000001_d+CDZ*I_ERI_Px_S_Py_S_C1000001_d;
  Double I_ERI_Py_S_Py_Pz_C1000001_d = I_ERI_Py_S_Dyz_S_C1000001_d+CDZ*I_ERI_Py_S_Py_S_C1000001_d;
  Double I_ERI_Pz_S_Py_Pz_C1000001_d = I_ERI_Pz_S_Dyz_S_C1000001_d+CDZ*I_ERI_Pz_S_Py_S_C1000001_d;
  Double I_ERI_Px_S_Pz_Pz_C1000001_d = I_ERI_Px_S_D2z_S_C1000001_d+CDZ*I_ERI_Px_S_Pz_S_C1000001_d;
  Double I_ERI_Py_S_Pz_Pz_C1000001_d = I_ERI_Py_S_D2z_S_C1000001_d+CDZ*I_ERI_Py_S_Pz_S_C1000001_d;
  Double I_ERI_Pz_S_Pz_Pz_C1000001_d = I_ERI_Pz_S_D2z_S_C1000001_d+CDZ*I_ERI_Pz_S_Pz_S_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P_C1001000_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_d
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001000_d
   ************************************************************/
  Double I_ERI_S_Px_Px_Px_C1001000_d = I_ERI_S_Px_D2x_S_C1001000_d+CDX*I_ERI_S_Px_Px_S_C1001000_d;
  Double I_ERI_S_Py_Px_Px_C1001000_d = I_ERI_S_Py_D2x_S_C1001000_d+CDX*I_ERI_S_Py_Px_S_C1001000_d;
  Double I_ERI_S_Pz_Px_Px_C1001000_d = I_ERI_S_Pz_D2x_S_C1001000_d+CDX*I_ERI_S_Pz_Px_S_C1001000_d;
  Double I_ERI_S_Px_Py_Px_C1001000_d = I_ERI_S_Px_Dxy_S_C1001000_d+CDX*I_ERI_S_Px_Py_S_C1001000_d;
  Double I_ERI_S_Py_Py_Px_C1001000_d = I_ERI_S_Py_Dxy_S_C1001000_d+CDX*I_ERI_S_Py_Py_S_C1001000_d;
  Double I_ERI_S_Pz_Py_Px_C1001000_d = I_ERI_S_Pz_Dxy_S_C1001000_d+CDX*I_ERI_S_Pz_Py_S_C1001000_d;
  Double I_ERI_S_Px_Pz_Px_C1001000_d = I_ERI_S_Px_Dxz_S_C1001000_d+CDX*I_ERI_S_Px_Pz_S_C1001000_d;
  Double I_ERI_S_Py_Pz_Px_C1001000_d = I_ERI_S_Py_Dxz_S_C1001000_d+CDX*I_ERI_S_Py_Pz_S_C1001000_d;
  Double I_ERI_S_Pz_Pz_Px_C1001000_d = I_ERI_S_Pz_Dxz_S_C1001000_d+CDX*I_ERI_S_Pz_Pz_S_C1001000_d;
  Double I_ERI_S_Px_Px_Py_C1001000_d = I_ERI_S_Px_Dxy_S_C1001000_d+CDY*I_ERI_S_Px_Px_S_C1001000_d;
  Double I_ERI_S_Py_Px_Py_C1001000_d = I_ERI_S_Py_Dxy_S_C1001000_d+CDY*I_ERI_S_Py_Px_S_C1001000_d;
  Double I_ERI_S_Pz_Px_Py_C1001000_d = I_ERI_S_Pz_Dxy_S_C1001000_d+CDY*I_ERI_S_Pz_Px_S_C1001000_d;
  Double I_ERI_S_Px_Py_Py_C1001000_d = I_ERI_S_Px_D2y_S_C1001000_d+CDY*I_ERI_S_Px_Py_S_C1001000_d;
  Double I_ERI_S_Py_Py_Py_C1001000_d = I_ERI_S_Py_D2y_S_C1001000_d+CDY*I_ERI_S_Py_Py_S_C1001000_d;
  Double I_ERI_S_Pz_Py_Py_C1001000_d = I_ERI_S_Pz_D2y_S_C1001000_d+CDY*I_ERI_S_Pz_Py_S_C1001000_d;
  Double I_ERI_S_Px_Pz_Py_C1001000_d = I_ERI_S_Px_Dyz_S_C1001000_d+CDY*I_ERI_S_Px_Pz_S_C1001000_d;
  Double I_ERI_S_Py_Pz_Py_C1001000_d = I_ERI_S_Py_Dyz_S_C1001000_d+CDY*I_ERI_S_Py_Pz_S_C1001000_d;
  Double I_ERI_S_Pz_Pz_Py_C1001000_d = I_ERI_S_Pz_Dyz_S_C1001000_d+CDY*I_ERI_S_Pz_Pz_S_C1001000_d;
  Double I_ERI_S_Px_Px_Pz_C1001000_d = I_ERI_S_Px_Dxz_S_C1001000_d+CDZ*I_ERI_S_Px_Px_S_C1001000_d;
  Double I_ERI_S_Py_Px_Pz_C1001000_d = I_ERI_S_Py_Dxz_S_C1001000_d+CDZ*I_ERI_S_Py_Px_S_C1001000_d;
  Double I_ERI_S_Pz_Px_Pz_C1001000_d = I_ERI_S_Pz_Dxz_S_C1001000_d+CDZ*I_ERI_S_Pz_Px_S_C1001000_d;
  Double I_ERI_S_Px_Py_Pz_C1001000_d = I_ERI_S_Px_Dyz_S_C1001000_d+CDZ*I_ERI_S_Px_Py_S_C1001000_d;
  Double I_ERI_S_Py_Py_Pz_C1001000_d = I_ERI_S_Py_Dyz_S_C1001000_d+CDZ*I_ERI_S_Py_Py_S_C1001000_d;
  Double I_ERI_S_Pz_Py_Pz_C1001000_d = I_ERI_S_Pz_Dyz_S_C1001000_d+CDZ*I_ERI_S_Pz_Py_S_C1001000_d;
  Double I_ERI_S_Px_Pz_Pz_C1001000_d = I_ERI_S_Px_D2z_S_C1001000_d+CDZ*I_ERI_S_Px_Pz_S_C1001000_d;
  Double I_ERI_S_Py_Pz_Pz_C1001000_d = I_ERI_S_Py_D2z_S_C1001000_d+CDZ*I_ERI_S_Py_Pz_S_C1001000_d;
  Double I_ERI_S_Pz_Pz_Pz_C1001000_d = I_ERI_S_Pz_D2z_S_C1001000_d+CDZ*I_ERI_S_Pz_Pz_S_C1001000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_C1001001_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1001001_d
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001001_d
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_C1001001_d = I_ERI_Px_S_D2x_S_C1001001_d+CDX*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_S_Px_Px_C1001001_d = I_ERI_Py_S_D2x_S_C1001001_d+CDX*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_S_Px_Px_C1001001_d = I_ERI_Pz_S_D2x_S_C1001001_d+CDX*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_S_Py_Px_C1001001_d = I_ERI_Px_S_Dxy_S_C1001001_d+CDX*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_S_Py_Px_C1001001_d = I_ERI_Py_S_Dxy_S_C1001001_d+CDX*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_S_Py_Px_C1001001_d = I_ERI_Pz_S_Dxy_S_C1001001_d+CDX*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_S_Pz_Px_C1001001_d = I_ERI_Px_S_Dxz_S_C1001001_d+CDX*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_S_Pz_Px_C1001001_d = I_ERI_Py_S_Dxz_S_C1001001_d+CDX*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_S_Pz_Px_C1001001_d = I_ERI_Pz_S_Dxz_S_C1001001_d+CDX*I_ERI_Pz_S_Pz_S_C1001001_d;
  Double I_ERI_Px_S_Px_Py_C1001001_d = I_ERI_Px_S_Dxy_S_C1001001_d+CDY*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_S_Px_Py_C1001001_d = I_ERI_Py_S_Dxy_S_C1001001_d+CDY*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_S_Px_Py_C1001001_d = I_ERI_Pz_S_Dxy_S_C1001001_d+CDY*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_S_Py_Py_C1001001_d = I_ERI_Px_S_D2y_S_C1001001_d+CDY*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_S_Py_Py_C1001001_d = I_ERI_Py_S_D2y_S_C1001001_d+CDY*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_S_Py_Py_C1001001_d = I_ERI_Pz_S_D2y_S_C1001001_d+CDY*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_S_Pz_Py_C1001001_d = I_ERI_Px_S_Dyz_S_C1001001_d+CDY*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_S_Pz_Py_C1001001_d = I_ERI_Py_S_Dyz_S_C1001001_d+CDY*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_S_Pz_Py_C1001001_d = I_ERI_Pz_S_Dyz_S_C1001001_d+CDY*I_ERI_Pz_S_Pz_S_C1001001_d;
  Double I_ERI_Px_S_Px_Pz_C1001001_d = I_ERI_Px_S_Dxz_S_C1001001_d+CDZ*I_ERI_Px_S_Px_S_C1001001_d;
  Double I_ERI_Py_S_Px_Pz_C1001001_d = I_ERI_Py_S_Dxz_S_C1001001_d+CDZ*I_ERI_Py_S_Px_S_C1001001_d;
  Double I_ERI_Pz_S_Px_Pz_C1001001_d = I_ERI_Pz_S_Dxz_S_C1001001_d+CDZ*I_ERI_Pz_S_Px_S_C1001001_d;
  Double I_ERI_Px_S_Py_Pz_C1001001_d = I_ERI_Px_S_Dyz_S_C1001001_d+CDZ*I_ERI_Px_S_Py_S_C1001001_d;
  Double I_ERI_Py_S_Py_Pz_C1001001_d = I_ERI_Py_S_Dyz_S_C1001001_d+CDZ*I_ERI_Py_S_Py_S_C1001001_d;
  Double I_ERI_Pz_S_Py_Pz_C1001001_d = I_ERI_Pz_S_Dyz_S_C1001001_d+CDZ*I_ERI_Pz_S_Py_S_C1001001_d;
  Double I_ERI_Px_S_Pz_Pz_C1001001_d = I_ERI_Px_S_D2z_S_C1001001_d+CDZ*I_ERI_Px_S_Pz_S_C1001001_d;
  Double I_ERI_Py_S_Pz_Pz_C1001001_d = I_ERI_Py_S_D2z_S_C1001001_d+CDZ*I_ERI_Py_S_Pz_S_C1001001_d;
  Double I_ERI_Pz_S_Pz_Pz_C1001001_d = I_ERI_Pz_S_D2z_S_C1001001_d+CDZ*I_ERI_Pz_S_Pz_S_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_C1001001_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001001_d
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001001_d
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001001_d = I_ERI_D2x_S_D2x_S_C1001001_d+CDX*I_ERI_D2x_S_Px_S_C1001001_d;
  Double I_ERI_Dxy_S_Px_Px_C1001001_d = I_ERI_Dxy_S_D2x_S_C1001001_d+CDX*I_ERI_Dxy_S_Px_S_C1001001_d;
  Double I_ERI_Dxz_S_Px_Px_C1001001_d = I_ERI_Dxz_S_D2x_S_C1001001_d+CDX*I_ERI_Dxz_S_Px_S_C1001001_d;
  Double I_ERI_D2y_S_Px_Px_C1001001_d = I_ERI_D2y_S_D2x_S_C1001001_d+CDX*I_ERI_D2y_S_Px_S_C1001001_d;
  Double I_ERI_Dyz_S_Px_Px_C1001001_d = I_ERI_Dyz_S_D2x_S_C1001001_d+CDX*I_ERI_Dyz_S_Px_S_C1001001_d;
  Double I_ERI_D2z_S_Px_Px_C1001001_d = I_ERI_D2z_S_D2x_S_C1001001_d+CDX*I_ERI_D2z_S_Px_S_C1001001_d;
  Double I_ERI_D2x_S_Py_Px_C1001001_d = I_ERI_D2x_S_Dxy_S_C1001001_d+CDX*I_ERI_D2x_S_Py_S_C1001001_d;
  Double I_ERI_Dxy_S_Py_Px_C1001001_d = I_ERI_Dxy_S_Dxy_S_C1001001_d+CDX*I_ERI_Dxy_S_Py_S_C1001001_d;
  Double I_ERI_Dxz_S_Py_Px_C1001001_d = I_ERI_Dxz_S_Dxy_S_C1001001_d+CDX*I_ERI_Dxz_S_Py_S_C1001001_d;
  Double I_ERI_D2y_S_Py_Px_C1001001_d = I_ERI_D2y_S_Dxy_S_C1001001_d+CDX*I_ERI_D2y_S_Py_S_C1001001_d;
  Double I_ERI_Dyz_S_Py_Px_C1001001_d = I_ERI_Dyz_S_Dxy_S_C1001001_d+CDX*I_ERI_Dyz_S_Py_S_C1001001_d;
  Double I_ERI_D2z_S_Py_Px_C1001001_d = I_ERI_D2z_S_Dxy_S_C1001001_d+CDX*I_ERI_D2z_S_Py_S_C1001001_d;
  Double I_ERI_D2x_S_Pz_Px_C1001001_d = I_ERI_D2x_S_Dxz_S_C1001001_d+CDX*I_ERI_D2x_S_Pz_S_C1001001_d;
  Double I_ERI_Dxy_S_Pz_Px_C1001001_d = I_ERI_Dxy_S_Dxz_S_C1001001_d+CDX*I_ERI_Dxy_S_Pz_S_C1001001_d;
  Double I_ERI_Dxz_S_Pz_Px_C1001001_d = I_ERI_Dxz_S_Dxz_S_C1001001_d+CDX*I_ERI_Dxz_S_Pz_S_C1001001_d;
  Double I_ERI_D2y_S_Pz_Px_C1001001_d = I_ERI_D2y_S_Dxz_S_C1001001_d+CDX*I_ERI_D2y_S_Pz_S_C1001001_d;
  Double I_ERI_Dyz_S_Pz_Px_C1001001_d = I_ERI_Dyz_S_Dxz_S_C1001001_d+CDX*I_ERI_Dyz_S_Pz_S_C1001001_d;
  Double I_ERI_D2z_S_Pz_Px_C1001001_d = I_ERI_D2z_S_Dxz_S_C1001001_d+CDX*I_ERI_D2z_S_Pz_S_C1001001_d;
  Double I_ERI_D2x_S_Px_Py_C1001001_d = I_ERI_D2x_S_Dxy_S_C1001001_d+CDY*I_ERI_D2x_S_Px_S_C1001001_d;
  Double I_ERI_Dxy_S_Px_Py_C1001001_d = I_ERI_Dxy_S_Dxy_S_C1001001_d+CDY*I_ERI_Dxy_S_Px_S_C1001001_d;
  Double I_ERI_Dxz_S_Px_Py_C1001001_d = I_ERI_Dxz_S_Dxy_S_C1001001_d+CDY*I_ERI_Dxz_S_Px_S_C1001001_d;
  Double I_ERI_D2y_S_Px_Py_C1001001_d = I_ERI_D2y_S_Dxy_S_C1001001_d+CDY*I_ERI_D2y_S_Px_S_C1001001_d;
  Double I_ERI_Dyz_S_Px_Py_C1001001_d = I_ERI_Dyz_S_Dxy_S_C1001001_d+CDY*I_ERI_Dyz_S_Px_S_C1001001_d;
  Double I_ERI_D2z_S_Px_Py_C1001001_d = I_ERI_D2z_S_Dxy_S_C1001001_d+CDY*I_ERI_D2z_S_Px_S_C1001001_d;
  Double I_ERI_D2x_S_Py_Py_C1001001_d = I_ERI_D2x_S_D2y_S_C1001001_d+CDY*I_ERI_D2x_S_Py_S_C1001001_d;
  Double I_ERI_Dxy_S_Py_Py_C1001001_d = I_ERI_Dxy_S_D2y_S_C1001001_d+CDY*I_ERI_Dxy_S_Py_S_C1001001_d;
  Double I_ERI_Dxz_S_Py_Py_C1001001_d = I_ERI_Dxz_S_D2y_S_C1001001_d+CDY*I_ERI_Dxz_S_Py_S_C1001001_d;
  Double I_ERI_D2y_S_Py_Py_C1001001_d = I_ERI_D2y_S_D2y_S_C1001001_d+CDY*I_ERI_D2y_S_Py_S_C1001001_d;
  Double I_ERI_Dyz_S_Py_Py_C1001001_d = I_ERI_Dyz_S_D2y_S_C1001001_d+CDY*I_ERI_Dyz_S_Py_S_C1001001_d;
  Double I_ERI_D2z_S_Py_Py_C1001001_d = I_ERI_D2z_S_D2y_S_C1001001_d+CDY*I_ERI_D2z_S_Py_S_C1001001_d;
  Double I_ERI_D2x_S_Pz_Py_C1001001_d = I_ERI_D2x_S_Dyz_S_C1001001_d+CDY*I_ERI_D2x_S_Pz_S_C1001001_d;
  Double I_ERI_Dxy_S_Pz_Py_C1001001_d = I_ERI_Dxy_S_Dyz_S_C1001001_d+CDY*I_ERI_Dxy_S_Pz_S_C1001001_d;
  Double I_ERI_Dxz_S_Pz_Py_C1001001_d = I_ERI_Dxz_S_Dyz_S_C1001001_d+CDY*I_ERI_Dxz_S_Pz_S_C1001001_d;
  Double I_ERI_D2y_S_Pz_Py_C1001001_d = I_ERI_D2y_S_Dyz_S_C1001001_d+CDY*I_ERI_D2y_S_Pz_S_C1001001_d;
  Double I_ERI_Dyz_S_Pz_Py_C1001001_d = I_ERI_Dyz_S_Dyz_S_C1001001_d+CDY*I_ERI_Dyz_S_Pz_S_C1001001_d;
  Double I_ERI_D2z_S_Pz_Py_C1001001_d = I_ERI_D2z_S_Dyz_S_C1001001_d+CDY*I_ERI_D2z_S_Pz_S_C1001001_d;
  Double I_ERI_D2x_S_Px_Pz_C1001001_d = I_ERI_D2x_S_Dxz_S_C1001001_d+CDZ*I_ERI_D2x_S_Px_S_C1001001_d;
  Double I_ERI_Dxy_S_Px_Pz_C1001001_d = I_ERI_Dxy_S_Dxz_S_C1001001_d+CDZ*I_ERI_Dxy_S_Px_S_C1001001_d;
  Double I_ERI_Dxz_S_Px_Pz_C1001001_d = I_ERI_Dxz_S_Dxz_S_C1001001_d+CDZ*I_ERI_Dxz_S_Px_S_C1001001_d;
  Double I_ERI_D2y_S_Px_Pz_C1001001_d = I_ERI_D2y_S_Dxz_S_C1001001_d+CDZ*I_ERI_D2y_S_Px_S_C1001001_d;
  Double I_ERI_Dyz_S_Px_Pz_C1001001_d = I_ERI_Dyz_S_Dxz_S_C1001001_d+CDZ*I_ERI_Dyz_S_Px_S_C1001001_d;
  Double I_ERI_D2z_S_Px_Pz_C1001001_d = I_ERI_D2z_S_Dxz_S_C1001001_d+CDZ*I_ERI_D2z_S_Px_S_C1001001_d;
  Double I_ERI_D2x_S_Py_Pz_C1001001_d = I_ERI_D2x_S_Dyz_S_C1001001_d+CDZ*I_ERI_D2x_S_Py_S_C1001001_d;
  Double I_ERI_Dxy_S_Py_Pz_C1001001_d = I_ERI_Dxy_S_Dyz_S_C1001001_d+CDZ*I_ERI_Dxy_S_Py_S_C1001001_d;
  Double I_ERI_Dxz_S_Py_Pz_C1001001_d = I_ERI_Dxz_S_Dyz_S_C1001001_d+CDZ*I_ERI_Dxz_S_Py_S_C1001001_d;
  Double I_ERI_D2y_S_Py_Pz_C1001001_d = I_ERI_D2y_S_Dyz_S_C1001001_d+CDZ*I_ERI_D2y_S_Py_S_C1001001_d;
  Double I_ERI_Dyz_S_Py_Pz_C1001001_d = I_ERI_Dyz_S_Dyz_S_C1001001_d+CDZ*I_ERI_Dyz_S_Py_S_C1001001_d;
  Double I_ERI_D2z_S_Py_Pz_C1001001_d = I_ERI_D2z_S_Dyz_S_C1001001_d+CDZ*I_ERI_D2z_S_Py_S_C1001001_d;
  Double I_ERI_D2x_S_Pz_Pz_C1001001_d = I_ERI_D2x_S_D2z_S_C1001001_d+CDZ*I_ERI_D2x_S_Pz_S_C1001001_d;
  Double I_ERI_Dxy_S_Pz_Pz_C1001001_d = I_ERI_Dxy_S_D2z_S_C1001001_d+CDZ*I_ERI_Dxy_S_Pz_S_C1001001_d;
  Double I_ERI_Dxz_S_Pz_Pz_C1001001_d = I_ERI_Dxz_S_D2z_S_C1001001_d+CDZ*I_ERI_Dxz_S_Pz_S_C1001001_d;
  Double I_ERI_D2y_S_Pz_Pz_C1001001_d = I_ERI_D2y_S_D2z_S_C1001001_d+CDZ*I_ERI_D2y_S_Pz_S_C1001001_d;
  Double I_ERI_Dyz_S_Pz_Pz_C1001001_d = I_ERI_Dyz_S_D2z_S_C1001001_d+CDZ*I_ERI_Dyz_S_Pz_S_C1001001_d;
  Double I_ERI_D2z_S_Pz_Pz_C1001001_d = I_ERI_D2z_S_D2z_S_C1001001_d+CDZ*I_ERI_D2z_S_Pz_S_C1001001_d;

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
   * shell quartet name: SQ_ERI_P_P_P_P_C1001001_d
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001001_d
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1001001_d
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_C1001001_d = I_ERI_D2x_S_Px_Px_C1001001_d+ABX*I_ERI_Px_S_Px_Px_C1001001_d;
  Double I_ERI_Py_Px_Px_Px_C1001001_d = I_ERI_Dxy_S_Px_Px_C1001001_d+ABX*I_ERI_Py_S_Px_Px_C1001001_d;
  Double I_ERI_Pz_Px_Px_Px_C1001001_d = I_ERI_Dxz_S_Px_Px_C1001001_d+ABX*I_ERI_Pz_S_Px_Px_C1001001_d;
  Double I_ERI_Px_Py_Px_Px_C1001001_d = I_ERI_Dxy_S_Px_Px_C1001001_d+ABY*I_ERI_Px_S_Px_Px_C1001001_d;
  Double I_ERI_Py_Py_Px_Px_C1001001_d = I_ERI_D2y_S_Px_Px_C1001001_d+ABY*I_ERI_Py_S_Px_Px_C1001001_d;
  Double I_ERI_Pz_Py_Px_Px_C1001001_d = I_ERI_Dyz_S_Px_Px_C1001001_d+ABY*I_ERI_Pz_S_Px_Px_C1001001_d;
  Double I_ERI_Px_Pz_Px_Px_C1001001_d = I_ERI_Dxz_S_Px_Px_C1001001_d+ABZ*I_ERI_Px_S_Px_Px_C1001001_d;
  Double I_ERI_Py_Pz_Px_Px_C1001001_d = I_ERI_Dyz_S_Px_Px_C1001001_d+ABZ*I_ERI_Py_S_Px_Px_C1001001_d;
  Double I_ERI_Pz_Pz_Px_Px_C1001001_d = I_ERI_D2z_S_Px_Px_C1001001_d+ABZ*I_ERI_Pz_S_Px_Px_C1001001_d;
  Double I_ERI_Px_Px_Py_Px_C1001001_d = I_ERI_D2x_S_Py_Px_C1001001_d+ABX*I_ERI_Px_S_Py_Px_C1001001_d;
  Double I_ERI_Py_Px_Py_Px_C1001001_d = I_ERI_Dxy_S_Py_Px_C1001001_d+ABX*I_ERI_Py_S_Py_Px_C1001001_d;
  Double I_ERI_Pz_Px_Py_Px_C1001001_d = I_ERI_Dxz_S_Py_Px_C1001001_d+ABX*I_ERI_Pz_S_Py_Px_C1001001_d;
  Double I_ERI_Px_Py_Py_Px_C1001001_d = I_ERI_Dxy_S_Py_Px_C1001001_d+ABY*I_ERI_Px_S_Py_Px_C1001001_d;
  Double I_ERI_Py_Py_Py_Px_C1001001_d = I_ERI_D2y_S_Py_Px_C1001001_d+ABY*I_ERI_Py_S_Py_Px_C1001001_d;
  Double I_ERI_Pz_Py_Py_Px_C1001001_d = I_ERI_Dyz_S_Py_Px_C1001001_d+ABY*I_ERI_Pz_S_Py_Px_C1001001_d;
  Double I_ERI_Px_Pz_Py_Px_C1001001_d = I_ERI_Dxz_S_Py_Px_C1001001_d+ABZ*I_ERI_Px_S_Py_Px_C1001001_d;
  Double I_ERI_Py_Pz_Py_Px_C1001001_d = I_ERI_Dyz_S_Py_Px_C1001001_d+ABZ*I_ERI_Py_S_Py_Px_C1001001_d;
  Double I_ERI_Pz_Pz_Py_Px_C1001001_d = I_ERI_D2z_S_Py_Px_C1001001_d+ABZ*I_ERI_Pz_S_Py_Px_C1001001_d;
  Double I_ERI_Px_Px_Pz_Px_C1001001_d = I_ERI_D2x_S_Pz_Px_C1001001_d+ABX*I_ERI_Px_S_Pz_Px_C1001001_d;
  Double I_ERI_Py_Px_Pz_Px_C1001001_d = I_ERI_Dxy_S_Pz_Px_C1001001_d+ABX*I_ERI_Py_S_Pz_Px_C1001001_d;
  Double I_ERI_Pz_Px_Pz_Px_C1001001_d = I_ERI_Dxz_S_Pz_Px_C1001001_d+ABX*I_ERI_Pz_S_Pz_Px_C1001001_d;
  Double I_ERI_Px_Py_Pz_Px_C1001001_d = I_ERI_Dxy_S_Pz_Px_C1001001_d+ABY*I_ERI_Px_S_Pz_Px_C1001001_d;
  Double I_ERI_Py_Py_Pz_Px_C1001001_d = I_ERI_D2y_S_Pz_Px_C1001001_d+ABY*I_ERI_Py_S_Pz_Px_C1001001_d;
  Double I_ERI_Pz_Py_Pz_Px_C1001001_d = I_ERI_Dyz_S_Pz_Px_C1001001_d+ABY*I_ERI_Pz_S_Pz_Px_C1001001_d;
  Double I_ERI_Px_Pz_Pz_Px_C1001001_d = I_ERI_Dxz_S_Pz_Px_C1001001_d+ABZ*I_ERI_Px_S_Pz_Px_C1001001_d;
  Double I_ERI_Py_Pz_Pz_Px_C1001001_d = I_ERI_Dyz_S_Pz_Px_C1001001_d+ABZ*I_ERI_Py_S_Pz_Px_C1001001_d;
  Double I_ERI_Pz_Pz_Pz_Px_C1001001_d = I_ERI_D2z_S_Pz_Px_C1001001_d+ABZ*I_ERI_Pz_S_Pz_Px_C1001001_d;
  Double I_ERI_Px_Px_Px_Py_C1001001_d = I_ERI_D2x_S_Px_Py_C1001001_d+ABX*I_ERI_Px_S_Px_Py_C1001001_d;
  Double I_ERI_Py_Px_Px_Py_C1001001_d = I_ERI_Dxy_S_Px_Py_C1001001_d+ABX*I_ERI_Py_S_Px_Py_C1001001_d;
  Double I_ERI_Pz_Px_Px_Py_C1001001_d = I_ERI_Dxz_S_Px_Py_C1001001_d+ABX*I_ERI_Pz_S_Px_Py_C1001001_d;
  Double I_ERI_Px_Py_Px_Py_C1001001_d = I_ERI_Dxy_S_Px_Py_C1001001_d+ABY*I_ERI_Px_S_Px_Py_C1001001_d;
  Double I_ERI_Py_Py_Px_Py_C1001001_d = I_ERI_D2y_S_Px_Py_C1001001_d+ABY*I_ERI_Py_S_Px_Py_C1001001_d;
  Double I_ERI_Pz_Py_Px_Py_C1001001_d = I_ERI_Dyz_S_Px_Py_C1001001_d+ABY*I_ERI_Pz_S_Px_Py_C1001001_d;
  Double I_ERI_Px_Pz_Px_Py_C1001001_d = I_ERI_Dxz_S_Px_Py_C1001001_d+ABZ*I_ERI_Px_S_Px_Py_C1001001_d;
  Double I_ERI_Py_Pz_Px_Py_C1001001_d = I_ERI_Dyz_S_Px_Py_C1001001_d+ABZ*I_ERI_Py_S_Px_Py_C1001001_d;
  Double I_ERI_Pz_Pz_Px_Py_C1001001_d = I_ERI_D2z_S_Px_Py_C1001001_d+ABZ*I_ERI_Pz_S_Px_Py_C1001001_d;
  Double I_ERI_Px_Px_Py_Py_C1001001_d = I_ERI_D2x_S_Py_Py_C1001001_d+ABX*I_ERI_Px_S_Py_Py_C1001001_d;
  Double I_ERI_Py_Px_Py_Py_C1001001_d = I_ERI_Dxy_S_Py_Py_C1001001_d+ABX*I_ERI_Py_S_Py_Py_C1001001_d;
  Double I_ERI_Pz_Px_Py_Py_C1001001_d = I_ERI_Dxz_S_Py_Py_C1001001_d+ABX*I_ERI_Pz_S_Py_Py_C1001001_d;
  Double I_ERI_Px_Py_Py_Py_C1001001_d = I_ERI_Dxy_S_Py_Py_C1001001_d+ABY*I_ERI_Px_S_Py_Py_C1001001_d;
  Double I_ERI_Py_Py_Py_Py_C1001001_d = I_ERI_D2y_S_Py_Py_C1001001_d+ABY*I_ERI_Py_S_Py_Py_C1001001_d;
  Double I_ERI_Pz_Py_Py_Py_C1001001_d = I_ERI_Dyz_S_Py_Py_C1001001_d+ABY*I_ERI_Pz_S_Py_Py_C1001001_d;
  Double I_ERI_Px_Pz_Py_Py_C1001001_d = I_ERI_Dxz_S_Py_Py_C1001001_d+ABZ*I_ERI_Px_S_Py_Py_C1001001_d;
  Double I_ERI_Py_Pz_Py_Py_C1001001_d = I_ERI_Dyz_S_Py_Py_C1001001_d+ABZ*I_ERI_Py_S_Py_Py_C1001001_d;
  Double I_ERI_Pz_Pz_Py_Py_C1001001_d = I_ERI_D2z_S_Py_Py_C1001001_d+ABZ*I_ERI_Pz_S_Py_Py_C1001001_d;
  Double I_ERI_Px_Px_Pz_Py_C1001001_d = I_ERI_D2x_S_Pz_Py_C1001001_d+ABX*I_ERI_Px_S_Pz_Py_C1001001_d;
  Double I_ERI_Py_Px_Pz_Py_C1001001_d = I_ERI_Dxy_S_Pz_Py_C1001001_d+ABX*I_ERI_Py_S_Pz_Py_C1001001_d;
  Double I_ERI_Pz_Px_Pz_Py_C1001001_d = I_ERI_Dxz_S_Pz_Py_C1001001_d+ABX*I_ERI_Pz_S_Pz_Py_C1001001_d;
  Double I_ERI_Px_Py_Pz_Py_C1001001_d = I_ERI_Dxy_S_Pz_Py_C1001001_d+ABY*I_ERI_Px_S_Pz_Py_C1001001_d;
  Double I_ERI_Py_Py_Pz_Py_C1001001_d = I_ERI_D2y_S_Pz_Py_C1001001_d+ABY*I_ERI_Py_S_Pz_Py_C1001001_d;
  Double I_ERI_Pz_Py_Pz_Py_C1001001_d = I_ERI_Dyz_S_Pz_Py_C1001001_d+ABY*I_ERI_Pz_S_Pz_Py_C1001001_d;
  Double I_ERI_Px_Pz_Pz_Py_C1001001_d = I_ERI_Dxz_S_Pz_Py_C1001001_d+ABZ*I_ERI_Px_S_Pz_Py_C1001001_d;
  Double I_ERI_Py_Pz_Pz_Py_C1001001_d = I_ERI_Dyz_S_Pz_Py_C1001001_d+ABZ*I_ERI_Py_S_Pz_Py_C1001001_d;
  Double I_ERI_Pz_Pz_Pz_Py_C1001001_d = I_ERI_D2z_S_Pz_Py_C1001001_d+ABZ*I_ERI_Pz_S_Pz_Py_C1001001_d;
  Double I_ERI_Px_Px_Px_Pz_C1001001_d = I_ERI_D2x_S_Px_Pz_C1001001_d+ABX*I_ERI_Px_S_Px_Pz_C1001001_d;
  Double I_ERI_Py_Px_Px_Pz_C1001001_d = I_ERI_Dxy_S_Px_Pz_C1001001_d+ABX*I_ERI_Py_S_Px_Pz_C1001001_d;
  Double I_ERI_Pz_Px_Px_Pz_C1001001_d = I_ERI_Dxz_S_Px_Pz_C1001001_d+ABX*I_ERI_Pz_S_Px_Pz_C1001001_d;
  Double I_ERI_Px_Py_Px_Pz_C1001001_d = I_ERI_Dxy_S_Px_Pz_C1001001_d+ABY*I_ERI_Px_S_Px_Pz_C1001001_d;
  Double I_ERI_Py_Py_Px_Pz_C1001001_d = I_ERI_D2y_S_Px_Pz_C1001001_d+ABY*I_ERI_Py_S_Px_Pz_C1001001_d;
  Double I_ERI_Pz_Py_Px_Pz_C1001001_d = I_ERI_Dyz_S_Px_Pz_C1001001_d+ABY*I_ERI_Pz_S_Px_Pz_C1001001_d;
  Double I_ERI_Px_Pz_Px_Pz_C1001001_d = I_ERI_Dxz_S_Px_Pz_C1001001_d+ABZ*I_ERI_Px_S_Px_Pz_C1001001_d;
  Double I_ERI_Py_Pz_Px_Pz_C1001001_d = I_ERI_Dyz_S_Px_Pz_C1001001_d+ABZ*I_ERI_Py_S_Px_Pz_C1001001_d;
  Double I_ERI_Pz_Pz_Px_Pz_C1001001_d = I_ERI_D2z_S_Px_Pz_C1001001_d+ABZ*I_ERI_Pz_S_Px_Pz_C1001001_d;
  Double I_ERI_Px_Px_Py_Pz_C1001001_d = I_ERI_D2x_S_Py_Pz_C1001001_d+ABX*I_ERI_Px_S_Py_Pz_C1001001_d;
  Double I_ERI_Py_Px_Py_Pz_C1001001_d = I_ERI_Dxy_S_Py_Pz_C1001001_d+ABX*I_ERI_Py_S_Py_Pz_C1001001_d;
  Double I_ERI_Pz_Px_Py_Pz_C1001001_d = I_ERI_Dxz_S_Py_Pz_C1001001_d+ABX*I_ERI_Pz_S_Py_Pz_C1001001_d;
  Double I_ERI_Px_Py_Py_Pz_C1001001_d = I_ERI_Dxy_S_Py_Pz_C1001001_d+ABY*I_ERI_Px_S_Py_Pz_C1001001_d;
  Double I_ERI_Py_Py_Py_Pz_C1001001_d = I_ERI_D2y_S_Py_Pz_C1001001_d+ABY*I_ERI_Py_S_Py_Pz_C1001001_d;
  Double I_ERI_Pz_Py_Py_Pz_C1001001_d = I_ERI_Dyz_S_Py_Pz_C1001001_d+ABY*I_ERI_Pz_S_Py_Pz_C1001001_d;
  Double I_ERI_Px_Pz_Py_Pz_C1001001_d = I_ERI_Dxz_S_Py_Pz_C1001001_d+ABZ*I_ERI_Px_S_Py_Pz_C1001001_d;
  Double I_ERI_Py_Pz_Py_Pz_C1001001_d = I_ERI_Dyz_S_Py_Pz_C1001001_d+ABZ*I_ERI_Py_S_Py_Pz_C1001001_d;
  Double I_ERI_Pz_Pz_Py_Pz_C1001001_d = I_ERI_D2z_S_Py_Pz_C1001001_d+ABZ*I_ERI_Pz_S_Py_Pz_C1001001_d;
  Double I_ERI_Px_Px_Pz_Pz_C1001001_d = I_ERI_D2x_S_Pz_Pz_C1001001_d+ABX*I_ERI_Px_S_Pz_Pz_C1001001_d;
  Double I_ERI_Py_Px_Pz_Pz_C1001001_d = I_ERI_Dxy_S_Pz_Pz_C1001001_d+ABX*I_ERI_Py_S_Pz_Pz_C1001001_d;
  Double I_ERI_Pz_Px_Pz_Pz_C1001001_d = I_ERI_Dxz_S_Pz_Pz_C1001001_d+ABX*I_ERI_Pz_S_Pz_Pz_C1001001_d;
  Double I_ERI_Px_Py_Pz_Pz_C1001001_d = I_ERI_Dxy_S_Pz_Pz_C1001001_d+ABY*I_ERI_Px_S_Pz_Pz_C1001001_d;
  Double I_ERI_Py_Py_Pz_Pz_C1001001_d = I_ERI_D2y_S_Pz_Pz_C1001001_d+ABY*I_ERI_Py_S_Pz_Pz_C1001001_d;
  Double I_ERI_Pz_Py_Pz_Pz_C1001001_d = I_ERI_Dyz_S_Pz_Pz_C1001001_d+ABY*I_ERI_Pz_S_Pz_Pz_C1001001_d;
  Double I_ERI_Px_Pz_Pz_Pz_C1001001_d = I_ERI_Dxz_S_Pz_Pz_C1001001_d+ABZ*I_ERI_Px_S_Pz_Pz_C1001001_d;
  Double I_ERI_Py_Pz_Pz_Pz_C1001001_d = I_ERI_Dyz_S_Pz_Pz_C1001001_d+ABZ*I_ERI_Py_S_Pz_Pz_C1001001_d;
  Double I_ERI_Pz_Pz_Pz_Pz_C1001001_d = I_ERI_D2z_S_Pz_Pz_C1001001_d+ABZ*I_ERI_Pz_S_Pz_Pz_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_Px_S_Px_S_C1000000_a;
  abcd[16] = 2.0E0*I_ERI_Px_S_Py_S_C1000000_a;
  abcd[32] = 2.0E0*I_ERI_Px_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[1] = 2.0E0*I_ERI_D2x_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[2] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[3] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[17] = 2.0E0*I_ERI_D2x_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[18] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[19] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[33] = 2.0E0*I_ERI_D2x_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[34] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[35] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[4] = 2.0E0*I_ERI_Px_Px_Px_S_C1001000_a;
  abcd[8] = 2.0E0*I_ERI_Px_Py_Px_S_C1001000_a;
  abcd[12] = 2.0E0*I_ERI_Px_Pz_Px_S_C1001000_a;
  abcd[20] = 2.0E0*I_ERI_Px_Px_Py_S_C1001000_a;
  abcd[24] = 2.0E0*I_ERI_Px_Py_Py_S_C1001000_a;
  abcd[28] = 2.0E0*I_ERI_Px_Pz_Py_S_C1001000_a;
  abcd[36] = 2.0E0*I_ERI_Px_Px_Pz_S_C1001000_a;
  abcd[40] = 2.0E0*I_ERI_Px_Py_Pz_S_C1001000_a;
  abcd[44] = 2.0E0*I_ERI_Px_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[5] = 2.0E0*I_ERI_D2x_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[6] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[7] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[9] = 2.0E0*I_ERI_D2x_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[10] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[11] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[13] = 2.0E0*I_ERI_D2x_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[14] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[15] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[21] = 2.0E0*I_ERI_D2x_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[22] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[23] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[25] = 2.0E0*I_ERI_D2x_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[26] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[27] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[29] = 2.0E0*I_ERI_D2x_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[30] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[31] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[37] = 2.0E0*I_ERI_D2x_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[38] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[39] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[41] = 2.0E0*I_ERI_D2x_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[42] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[43] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[45] = 2.0E0*I_ERI_D2x_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[46] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[47] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[48] = 2.0E0*I_ERI_Py_S_Px_S_C1000000_a;
  abcd[64] = 2.0E0*I_ERI_Py_S_Py_S_C1000000_a;
  abcd[80] = 2.0E0*I_ERI_Py_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[49] = 2.0E0*I_ERI_Dxy_S_Px_S_C1000001_a;
  abcd[50] = 2.0E0*I_ERI_D2y_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[51] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[65] = 2.0E0*I_ERI_Dxy_S_Py_S_C1000001_a;
  abcd[66] = 2.0E0*I_ERI_D2y_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[67] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[81] = 2.0E0*I_ERI_Dxy_S_Pz_S_C1000001_a;
  abcd[82] = 2.0E0*I_ERI_D2y_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;
  abcd[83] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[52] = 2.0E0*I_ERI_Py_Px_Px_S_C1001000_a;
  abcd[56] = 2.0E0*I_ERI_Py_Py_Px_S_C1001000_a;
  abcd[60] = 2.0E0*I_ERI_Py_Pz_Px_S_C1001000_a;
  abcd[68] = 2.0E0*I_ERI_Py_Px_Py_S_C1001000_a;
  abcd[72] = 2.0E0*I_ERI_Py_Py_Py_S_C1001000_a;
  abcd[76] = 2.0E0*I_ERI_Py_Pz_Py_S_C1001000_a;
  abcd[84] = 2.0E0*I_ERI_Py_Px_Pz_S_C1001000_a;
  abcd[88] = 2.0E0*I_ERI_Py_Py_Pz_S_C1001000_a;
  abcd[92] = 2.0E0*I_ERI_Py_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[53] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1001001_a;
  abcd[54] = 2.0E0*I_ERI_D2y_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[55] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[57] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1001001_a;
  abcd[58] = 2.0E0*I_ERI_D2y_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[59] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[61] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1001001_a;
  abcd[62] = 2.0E0*I_ERI_D2y_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[63] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[69] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1001001_a;
  abcd[70] = 2.0E0*I_ERI_D2y_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[71] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[73] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1001001_a;
  abcd[74] = 2.0E0*I_ERI_D2y_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[75] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[77] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1001001_a;
  abcd[78] = 2.0E0*I_ERI_D2y_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[79] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[85] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1001001_a;
  abcd[86] = 2.0E0*I_ERI_D2y_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[87] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[89] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1001001_a;
  abcd[90] = 2.0E0*I_ERI_D2y_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[91] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[93] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1001001_a;
  abcd[94] = 2.0E0*I_ERI_D2y_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;
  abcd[95] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000000_a
   ************************************************************/
  abcd[96] = 2.0E0*I_ERI_Pz_S_Px_S_C1000000_a;
  abcd[112] = 2.0E0*I_ERI_Pz_S_Py_S_C1000000_a;
  abcd[128] = 2.0E0*I_ERI_Pz_S_Pz_S_C1000000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000001_a
   * RHS shell quartet name: SQ_ERI_S_S_P_S_C1000001
   ************************************************************/
  abcd[97] = 2.0E0*I_ERI_Dxz_S_Px_S_C1000001_a;
  abcd[98] = 2.0E0*I_ERI_Dyz_S_Px_S_C1000001_a;
  abcd[99] = 2.0E0*I_ERI_D2z_S_Px_S_C1000001_a-1*I_ERI_S_S_Px_S_C1000001;
  abcd[113] = 2.0E0*I_ERI_Dxz_S_Py_S_C1000001_a;
  abcd[114] = 2.0E0*I_ERI_Dyz_S_Py_S_C1000001_a;
  abcd[115] = 2.0E0*I_ERI_D2z_S_Py_S_C1000001_a-1*I_ERI_S_S_Py_S_C1000001;
  abcd[129] = 2.0E0*I_ERI_Dxz_S_Pz_S_C1000001_a;
  abcd[130] = 2.0E0*I_ERI_Dyz_S_Pz_S_C1000001_a;
  abcd[131] = 2.0E0*I_ERI_D2z_S_Pz_S_C1000001_a-1*I_ERI_S_S_Pz_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001000_a
   ************************************************************/
  abcd[100] = 2.0E0*I_ERI_Pz_Px_Px_S_C1001000_a;
  abcd[104] = 2.0E0*I_ERI_Pz_Py_Px_S_C1001000_a;
  abcd[108] = 2.0E0*I_ERI_Pz_Pz_Px_S_C1001000_a;
  abcd[116] = 2.0E0*I_ERI_Pz_Px_Py_S_C1001000_a;
  abcd[120] = 2.0E0*I_ERI_Pz_Py_Py_S_C1001000_a;
  abcd[124] = 2.0E0*I_ERI_Pz_Pz_Py_S_C1001000_a;
  abcd[132] = 2.0E0*I_ERI_Pz_Px_Pz_S_C1001000_a;
  abcd[136] = 2.0E0*I_ERI_Pz_Py_Pz_S_C1001000_a;
  abcd[140] = 2.0E0*I_ERI_Pz_Pz_Pz_S_C1001000_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001001_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S_C1001001
   ************************************************************/
  abcd[101] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1001001_a;
  abcd[102] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1001001_a;
  abcd[103] = 2.0E0*I_ERI_D2z_Px_Px_S_C1001001_a-1*I_ERI_S_Px_Px_S_C1001001;
  abcd[105] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1001001_a;
  abcd[106] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1001001_a;
  abcd[107] = 2.0E0*I_ERI_D2z_Py_Px_S_C1001001_a-1*I_ERI_S_Py_Px_S_C1001001;
  abcd[109] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1001001_a;
  abcd[110] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1001001_a;
  abcd[111] = 2.0E0*I_ERI_D2z_Pz_Px_S_C1001001_a-1*I_ERI_S_Pz_Px_S_C1001001;
  abcd[117] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1001001_a;
  abcd[118] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1001001_a;
  abcd[119] = 2.0E0*I_ERI_D2z_Px_Py_S_C1001001_a-1*I_ERI_S_Px_Py_S_C1001001;
  abcd[121] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1001001_a;
  abcd[122] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1001001_a;
  abcd[123] = 2.0E0*I_ERI_D2z_Py_Py_S_C1001001_a-1*I_ERI_S_Py_Py_S_C1001001;
  abcd[125] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1001001_a;
  abcd[126] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1001001_a;
  abcd[127] = 2.0E0*I_ERI_D2z_Pz_Py_S_C1001001_a-1*I_ERI_S_Pz_Py_S_C1001001;
  abcd[133] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1001001_a;
  abcd[134] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1001001_a;
  abcd[135] = 2.0E0*I_ERI_D2z_Px_Pz_S_C1001001_a-1*I_ERI_S_Px_Pz_S_C1001001;
  abcd[137] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1001001_a;
  abcd[138] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1001001_a;
  abcd[139] = 2.0E0*I_ERI_D2z_Py_Pz_S_C1001001_a-1*I_ERI_S_Py_Pz_S_C1001001;
  abcd[141] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1001001_a;
  abcd[142] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1001001_a;
  abcd[143] = 2.0E0*I_ERI_D2z_Pz_Pz_S_C1001001_a-1*I_ERI_S_Pz_Pz_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[144] = 2.0E0*I_ERI_S_S_D2x_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;
  abcd[160] = 2.0E0*I_ERI_S_S_Dxy_S_C1000000_c;
  abcd[176] = 2.0E0*I_ERI_S_S_Dxz_S_C1000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[145] = 2.0E0*I_ERI_Px_S_D2x_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[146] = 2.0E0*I_ERI_Py_S_D2x_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[147] = 2.0E0*I_ERI_Pz_S_D2x_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[161] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[162] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[163] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[177] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[178] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[179] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[148] = 2.0E0*I_ERI_S_Px_D2x_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[152] = 2.0E0*I_ERI_S_Py_D2x_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[156] = 2.0E0*I_ERI_S_Pz_D2x_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;
  abcd[164] = 2.0E0*I_ERI_S_Px_Dxy_S_C1001000_c;
  abcd[168] = 2.0E0*I_ERI_S_Py_Dxy_S_C1001000_c;
  abcd[172] = 2.0E0*I_ERI_S_Pz_Dxy_S_C1001000_c;
  abcd[180] = 2.0E0*I_ERI_S_Px_Dxz_S_C1001000_c;
  abcd[184] = 2.0E0*I_ERI_S_Py_Dxz_S_C1001000_c;
  abcd[188] = 2.0E0*I_ERI_S_Pz_Dxz_S_C1001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[149] = 2.0E0*I_ERI_Px_Px_D2x_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[150] = 2.0E0*I_ERI_Py_Px_D2x_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[151] = 2.0E0*I_ERI_Pz_Px_D2x_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[153] = 2.0E0*I_ERI_Px_Py_D2x_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[154] = 2.0E0*I_ERI_Py_Py_D2x_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[155] = 2.0E0*I_ERI_Pz_Py_D2x_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[157] = 2.0E0*I_ERI_Px_Pz_D2x_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[158] = 2.0E0*I_ERI_Py_Pz_D2x_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[159] = 2.0E0*I_ERI_Pz_Pz_D2x_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[165] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[166] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[167] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[169] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[170] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[171] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[173] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[174] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[175] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[181] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[182] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[183] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[185] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[186] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[187] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[189] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[190] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[191] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[192] = 2.0E0*I_ERI_S_S_Dxy_S_C1000000_c;
  abcd[208] = 2.0E0*I_ERI_S_S_D2y_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;
  abcd[224] = 2.0E0*I_ERI_S_S_Dyz_S_C1000000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[193] = 2.0E0*I_ERI_Px_S_Dxy_S_C1000001_c;
  abcd[194] = 2.0E0*I_ERI_Py_S_Dxy_S_C1000001_c;
  abcd[195] = 2.0E0*I_ERI_Pz_S_Dxy_S_C1000001_c;
  abcd[209] = 2.0E0*I_ERI_Px_S_D2y_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[210] = 2.0E0*I_ERI_Py_S_D2y_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[211] = 2.0E0*I_ERI_Pz_S_D2y_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;
  abcd[225] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[226] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[227] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[196] = 2.0E0*I_ERI_S_Px_Dxy_S_C1001000_c;
  abcd[200] = 2.0E0*I_ERI_S_Py_Dxy_S_C1001000_c;
  abcd[204] = 2.0E0*I_ERI_S_Pz_Dxy_S_C1001000_c;
  abcd[212] = 2.0E0*I_ERI_S_Px_D2y_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[216] = 2.0E0*I_ERI_S_Py_D2y_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[220] = 2.0E0*I_ERI_S_Pz_D2y_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;
  abcd[228] = 2.0E0*I_ERI_S_Px_Dyz_S_C1001000_c;
  abcd[232] = 2.0E0*I_ERI_S_Py_Dyz_S_C1001000_c;
  abcd[236] = 2.0E0*I_ERI_S_Pz_Dyz_S_C1001000_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[197] = 2.0E0*I_ERI_Px_Px_Dxy_S_C1001001_c;
  abcd[198] = 2.0E0*I_ERI_Py_Px_Dxy_S_C1001001_c;
  abcd[199] = 2.0E0*I_ERI_Pz_Px_Dxy_S_C1001001_c;
  abcd[201] = 2.0E0*I_ERI_Px_Py_Dxy_S_C1001001_c;
  abcd[202] = 2.0E0*I_ERI_Py_Py_Dxy_S_C1001001_c;
  abcd[203] = 2.0E0*I_ERI_Pz_Py_Dxy_S_C1001001_c;
  abcd[205] = 2.0E0*I_ERI_Px_Pz_Dxy_S_C1001001_c;
  abcd[206] = 2.0E0*I_ERI_Py_Pz_Dxy_S_C1001001_c;
  abcd[207] = 2.0E0*I_ERI_Pz_Pz_Dxy_S_C1001001_c;
  abcd[213] = 2.0E0*I_ERI_Px_Px_D2y_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[214] = 2.0E0*I_ERI_Py_Px_D2y_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[215] = 2.0E0*I_ERI_Pz_Px_D2y_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[217] = 2.0E0*I_ERI_Px_Py_D2y_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[218] = 2.0E0*I_ERI_Py_Py_D2y_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[219] = 2.0E0*I_ERI_Pz_Py_D2y_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[221] = 2.0E0*I_ERI_Px_Pz_D2y_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[222] = 2.0E0*I_ERI_Py_Pz_D2y_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[223] = 2.0E0*I_ERI_Pz_Pz_D2y_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;
  abcd[229] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[230] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[231] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[233] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[234] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[235] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[237] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[238] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[239] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_D_S_C1000000_c
   * RHS shell quartet name: SQ_ERI_S_S_S_S_C1000000
   ************************************************************/
  abcd[240] = 2.0E0*I_ERI_S_S_Dxz_S_C1000000_c;
  abcd[256] = 2.0E0*I_ERI_S_S_Dyz_S_C1000000_c;
  abcd[272] = 2.0E0*I_ERI_S_S_D2z_S_C1000000_c-1*I_ERI_S_S_S_S_C1000000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_C1000001_c
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1000001
   ************************************************************/
  abcd[241] = 2.0E0*I_ERI_Px_S_Dxz_S_C1000001_c;
  abcd[242] = 2.0E0*I_ERI_Py_S_Dxz_S_C1000001_c;
  abcd[243] = 2.0E0*I_ERI_Pz_S_Dxz_S_C1000001_c;
  abcd[257] = 2.0E0*I_ERI_Px_S_Dyz_S_C1000001_c;
  abcd[258] = 2.0E0*I_ERI_Py_S_Dyz_S_C1000001_c;
  abcd[259] = 2.0E0*I_ERI_Pz_S_Dyz_S_C1000001_c;
  abcd[273] = 2.0E0*I_ERI_Px_S_D2z_S_C1000001_c-1*I_ERI_Px_S_S_S_C1000001;
  abcd[274] = 2.0E0*I_ERI_Py_S_D2z_S_C1000001_c-1*I_ERI_Py_S_S_S_C1000001;
  abcd[275] = 2.0E0*I_ERI_Pz_S_D2z_S_C1000001_c-1*I_ERI_Pz_S_S_S_C1000001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S_C1001000_c
   * RHS shell quartet name: SQ_ERI_S_P_S_S_C1001000
   ************************************************************/
  abcd[244] = 2.0E0*I_ERI_S_Px_Dxz_S_C1001000_c;
  abcd[248] = 2.0E0*I_ERI_S_Py_Dxz_S_C1001000_c;
  abcd[252] = 2.0E0*I_ERI_S_Pz_Dxz_S_C1001000_c;
  abcd[260] = 2.0E0*I_ERI_S_Px_Dyz_S_C1001000_c;
  abcd[264] = 2.0E0*I_ERI_S_Py_Dyz_S_C1001000_c;
  abcd[268] = 2.0E0*I_ERI_S_Pz_Dyz_S_C1001000_c;
  abcd[276] = 2.0E0*I_ERI_S_Px_D2z_S_C1001000_c-1*I_ERI_S_Px_S_S_C1001000;
  abcd[280] = 2.0E0*I_ERI_S_Py_D2z_S_C1001000_c-1*I_ERI_S_Py_S_S_C1001000;
  abcd[284] = 2.0E0*I_ERI_S_Pz_D2z_S_C1001000_c-1*I_ERI_S_Pz_S_S_C1001000;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_S_C1001001_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1001001
   ************************************************************/
  abcd[245] = 2.0E0*I_ERI_Px_Px_Dxz_S_C1001001_c;
  abcd[246] = 2.0E0*I_ERI_Py_Px_Dxz_S_C1001001_c;
  abcd[247] = 2.0E0*I_ERI_Pz_Px_Dxz_S_C1001001_c;
  abcd[249] = 2.0E0*I_ERI_Px_Py_Dxz_S_C1001001_c;
  abcd[250] = 2.0E0*I_ERI_Py_Py_Dxz_S_C1001001_c;
  abcd[251] = 2.0E0*I_ERI_Pz_Py_Dxz_S_C1001001_c;
  abcd[253] = 2.0E0*I_ERI_Px_Pz_Dxz_S_C1001001_c;
  abcd[254] = 2.0E0*I_ERI_Py_Pz_Dxz_S_C1001001_c;
  abcd[255] = 2.0E0*I_ERI_Pz_Pz_Dxz_S_C1001001_c;
  abcd[261] = 2.0E0*I_ERI_Px_Px_Dyz_S_C1001001_c;
  abcd[262] = 2.0E0*I_ERI_Py_Px_Dyz_S_C1001001_c;
  abcd[263] = 2.0E0*I_ERI_Pz_Px_Dyz_S_C1001001_c;
  abcd[265] = 2.0E0*I_ERI_Px_Py_Dyz_S_C1001001_c;
  abcd[266] = 2.0E0*I_ERI_Py_Py_Dyz_S_C1001001_c;
  abcd[267] = 2.0E0*I_ERI_Pz_Py_Dyz_S_C1001001_c;
  abcd[269] = 2.0E0*I_ERI_Px_Pz_Dyz_S_C1001001_c;
  abcd[270] = 2.0E0*I_ERI_Py_Pz_Dyz_S_C1001001_c;
  abcd[271] = 2.0E0*I_ERI_Pz_Pz_Dyz_S_C1001001_c;
  abcd[277] = 2.0E0*I_ERI_Px_Px_D2z_S_C1001001_c-1*I_ERI_Px_Px_S_S_C1001001;
  abcd[278] = 2.0E0*I_ERI_Py_Px_D2z_S_C1001001_c-1*I_ERI_Py_Px_S_S_C1001001;
  abcd[279] = 2.0E0*I_ERI_Pz_Px_D2z_S_C1001001_c-1*I_ERI_Pz_Px_S_S_C1001001;
  abcd[281] = 2.0E0*I_ERI_Px_Py_D2z_S_C1001001_c-1*I_ERI_Px_Py_S_S_C1001001;
  abcd[282] = 2.0E0*I_ERI_Py_Py_D2z_S_C1001001_c-1*I_ERI_Py_Py_S_S_C1001001;
  abcd[283] = 2.0E0*I_ERI_Pz_Py_D2z_S_C1001001_c-1*I_ERI_Pz_Py_S_S_C1001001;
  abcd[285] = 2.0E0*I_ERI_Px_Pz_D2z_S_C1001001_c-1*I_ERI_Px_Pz_S_S_C1001001;
  abcd[286] = 2.0E0*I_ERI_Py_Pz_D2z_S_C1001001_c-1*I_ERI_Py_Pz_S_S_C1001001;
  abcd[287] = 2.0E0*I_ERI_Pz_Pz_D2z_S_C1001001_c-1*I_ERI_Pz_Pz_S_S_C1001001;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000_d
   ************************************************************/
  abcd[288] = 2.0E0*I_ERI_S_S_Px_Px_C1000000_d;
  abcd[304] = 2.0E0*I_ERI_S_S_Py_Px_C1000000_d;
  abcd[320] = 2.0E0*I_ERI_S_S_Pz_Px_C1000000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000001_d
   ************************************************************/
  abcd[289] = 2.0E0*I_ERI_Px_S_Px_Px_C1000001_d;
  abcd[290] = 2.0E0*I_ERI_Py_S_Px_Px_C1000001_d;
  abcd[291] = 2.0E0*I_ERI_Pz_S_Px_Px_C1000001_d;
  abcd[305] = 2.0E0*I_ERI_Px_S_Py_Px_C1000001_d;
  abcd[306] = 2.0E0*I_ERI_Py_S_Py_Px_C1000001_d;
  abcd[307] = 2.0E0*I_ERI_Pz_S_Py_Px_C1000001_d;
  abcd[321] = 2.0E0*I_ERI_Px_S_Pz_Px_C1000001_d;
  abcd[322] = 2.0E0*I_ERI_Py_S_Pz_Px_C1000001_d;
  abcd[323] = 2.0E0*I_ERI_Pz_S_Pz_Px_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000_d
   ************************************************************/
  abcd[292] = 2.0E0*I_ERI_S_Px_Px_Px_C1001000_d;
  abcd[296] = 2.0E0*I_ERI_S_Py_Px_Px_C1001000_d;
  abcd[300] = 2.0E0*I_ERI_S_Pz_Px_Px_C1001000_d;
  abcd[308] = 2.0E0*I_ERI_S_Px_Py_Px_C1001000_d;
  abcd[312] = 2.0E0*I_ERI_S_Py_Py_Px_C1001000_d;
  abcd[316] = 2.0E0*I_ERI_S_Pz_Py_Px_C1001000_d;
  abcd[324] = 2.0E0*I_ERI_S_Px_Pz_Px_C1001000_d;
  abcd[328] = 2.0E0*I_ERI_S_Py_Pz_Px_C1001000_d;
  abcd[332] = 2.0E0*I_ERI_S_Pz_Pz_Px_C1001000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001_d
   ************************************************************/
  abcd[293] = 2.0E0*I_ERI_Px_Px_Px_Px_C1001001_d;
  abcd[294] = 2.0E0*I_ERI_Py_Px_Px_Px_C1001001_d;
  abcd[295] = 2.0E0*I_ERI_Pz_Px_Px_Px_C1001001_d;
  abcd[297] = 2.0E0*I_ERI_Px_Py_Px_Px_C1001001_d;
  abcd[298] = 2.0E0*I_ERI_Py_Py_Px_Px_C1001001_d;
  abcd[299] = 2.0E0*I_ERI_Pz_Py_Px_Px_C1001001_d;
  abcd[301] = 2.0E0*I_ERI_Px_Pz_Px_Px_C1001001_d;
  abcd[302] = 2.0E0*I_ERI_Py_Pz_Px_Px_C1001001_d;
  abcd[303] = 2.0E0*I_ERI_Pz_Pz_Px_Px_C1001001_d;
  abcd[309] = 2.0E0*I_ERI_Px_Px_Py_Px_C1001001_d;
  abcd[310] = 2.0E0*I_ERI_Py_Px_Py_Px_C1001001_d;
  abcd[311] = 2.0E0*I_ERI_Pz_Px_Py_Px_C1001001_d;
  abcd[313] = 2.0E0*I_ERI_Px_Py_Py_Px_C1001001_d;
  abcd[314] = 2.0E0*I_ERI_Py_Py_Py_Px_C1001001_d;
  abcd[315] = 2.0E0*I_ERI_Pz_Py_Py_Px_C1001001_d;
  abcd[317] = 2.0E0*I_ERI_Px_Pz_Py_Px_C1001001_d;
  abcd[318] = 2.0E0*I_ERI_Py_Pz_Py_Px_C1001001_d;
  abcd[319] = 2.0E0*I_ERI_Pz_Pz_Py_Px_C1001001_d;
  abcd[325] = 2.0E0*I_ERI_Px_Px_Pz_Px_C1001001_d;
  abcd[326] = 2.0E0*I_ERI_Py_Px_Pz_Px_C1001001_d;
  abcd[327] = 2.0E0*I_ERI_Pz_Px_Pz_Px_C1001001_d;
  abcd[329] = 2.0E0*I_ERI_Px_Py_Pz_Px_C1001001_d;
  abcd[330] = 2.0E0*I_ERI_Py_Py_Pz_Px_C1001001_d;
  abcd[331] = 2.0E0*I_ERI_Pz_Py_Pz_Px_C1001001_d;
  abcd[333] = 2.0E0*I_ERI_Px_Pz_Pz_Px_C1001001_d;
  abcd[334] = 2.0E0*I_ERI_Py_Pz_Pz_Px_C1001001_d;
  abcd[335] = 2.0E0*I_ERI_Pz_Pz_Pz_Px_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000_d
   ************************************************************/
  abcd[336] = 2.0E0*I_ERI_S_S_Px_Py_C1000000_d;
  abcd[352] = 2.0E0*I_ERI_S_S_Py_Py_C1000000_d;
  abcd[368] = 2.0E0*I_ERI_S_S_Pz_Py_C1000000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000001_d
   ************************************************************/
  abcd[337] = 2.0E0*I_ERI_Px_S_Px_Py_C1000001_d;
  abcd[338] = 2.0E0*I_ERI_Py_S_Px_Py_C1000001_d;
  abcd[339] = 2.0E0*I_ERI_Pz_S_Px_Py_C1000001_d;
  abcd[353] = 2.0E0*I_ERI_Px_S_Py_Py_C1000001_d;
  abcd[354] = 2.0E0*I_ERI_Py_S_Py_Py_C1000001_d;
  abcd[355] = 2.0E0*I_ERI_Pz_S_Py_Py_C1000001_d;
  abcd[369] = 2.0E0*I_ERI_Px_S_Pz_Py_C1000001_d;
  abcd[370] = 2.0E0*I_ERI_Py_S_Pz_Py_C1000001_d;
  abcd[371] = 2.0E0*I_ERI_Pz_S_Pz_Py_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000_d
   ************************************************************/
  abcd[340] = 2.0E0*I_ERI_S_Px_Px_Py_C1001000_d;
  abcd[344] = 2.0E0*I_ERI_S_Py_Px_Py_C1001000_d;
  abcd[348] = 2.0E0*I_ERI_S_Pz_Px_Py_C1001000_d;
  abcd[356] = 2.0E0*I_ERI_S_Px_Py_Py_C1001000_d;
  abcd[360] = 2.0E0*I_ERI_S_Py_Py_Py_C1001000_d;
  abcd[364] = 2.0E0*I_ERI_S_Pz_Py_Py_C1001000_d;
  abcd[372] = 2.0E0*I_ERI_S_Px_Pz_Py_C1001000_d;
  abcd[376] = 2.0E0*I_ERI_S_Py_Pz_Py_C1001000_d;
  abcd[380] = 2.0E0*I_ERI_S_Pz_Pz_Py_C1001000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001_d
   ************************************************************/
  abcd[341] = 2.0E0*I_ERI_Px_Px_Px_Py_C1001001_d;
  abcd[342] = 2.0E0*I_ERI_Py_Px_Px_Py_C1001001_d;
  abcd[343] = 2.0E0*I_ERI_Pz_Px_Px_Py_C1001001_d;
  abcd[345] = 2.0E0*I_ERI_Px_Py_Px_Py_C1001001_d;
  abcd[346] = 2.0E0*I_ERI_Py_Py_Px_Py_C1001001_d;
  abcd[347] = 2.0E0*I_ERI_Pz_Py_Px_Py_C1001001_d;
  abcd[349] = 2.0E0*I_ERI_Px_Pz_Px_Py_C1001001_d;
  abcd[350] = 2.0E0*I_ERI_Py_Pz_Px_Py_C1001001_d;
  abcd[351] = 2.0E0*I_ERI_Pz_Pz_Px_Py_C1001001_d;
  abcd[357] = 2.0E0*I_ERI_Px_Px_Py_Py_C1001001_d;
  abcd[358] = 2.0E0*I_ERI_Py_Px_Py_Py_C1001001_d;
  abcd[359] = 2.0E0*I_ERI_Pz_Px_Py_Py_C1001001_d;
  abcd[361] = 2.0E0*I_ERI_Px_Py_Py_Py_C1001001_d;
  abcd[362] = 2.0E0*I_ERI_Py_Py_Py_Py_C1001001_d;
  abcd[363] = 2.0E0*I_ERI_Pz_Py_Py_Py_C1001001_d;
  abcd[365] = 2.0E0*I_ERI_Px_Pz_Py_Py_C1001001_d;
  abcd[366] = 2.0E0*I_ERI_Py_Pz_Py_Py_C1001001_d;
  abcd[367] = 2.0E0*I_ERI_Pz_Pz_Py_Py_C1001001_d;
  abcd[373] = 2.0E0*I_ERI_Px_Px_Pz_Py_C1001001_d;
  abcd[374] = 2.0E0*I_ERI_Py_Px_Pz_Py_C1001001_d;
  abcd[375] = 2.0E0*I_ERI_Pz_Px_Pz_Py_C1001001_d;
  abcd[377] = 2.0E0*I_ERI_Px_Py_Pz_Py_C1001001_d;
  abcd[378] = 2.0E0*I_ERI_Py_Py_Pz_Py_C1001001_d;
  abcd[379] = 2.0E0*I_ERI_Pz_Py_Pz_Py_C1001001_d;
  abcd[381] = 2.0E0*I_ERI_Px_Pz_Pz_Py_C1001001_d;
  abcd[382] = 2.0E0*I_ERI_Py_Pz_Pz_Py_C1001001_d;
  abcd[383] = 2.0E0*I_ERI_Pz_Pz_Pz_Py_C1001001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_S_P_S_C1000000_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_S_P_P_C1000000_d
   ************************************************************/
  abcd[384] = 2.0E0*I_ERI_S_S_Px_Pz_C1000000_d;
  abcd[400] = 2.0E0*I_ERI_S_S_Py_Pz_C1000000_d;
  abcd[416] = 2.0E0*I_ERI_S_S_Pz_Pz_C1000000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_S_C1000001_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_P_P_C1000001_d
   ************************************************************/
  abcd[385] = 2.0E0*I_ERI_Px_S_Px_Pz_C1000001_d;
  abcd[386] = 2.0E0*I_ERI_Py_S_Px_Pz_C1000001_d;
  abcd[387] = 2.0E0*I_ERI_Pz_S_Px_Pz_C1000001_d;
  abcd[401] = 2.0E0*I_ERI_Px_S_Py_Pz_C1000001_d;
  abcd[402] = 2.0E0*I_ERI_Py_S_Py_Pz_C1000001_d;
  abcd[403] = 2.0E0*I_ERI_Pz_S_Py_Pz_C1000001_d;
  abcd[417] = 2.0E0*I_ERI_Px_S_Pz_Pz_C1000001_d;
  abcd[418] = 2.0E0*I_ERI_Py_S_Pz_Pz_C1000001_d;
  abcd[419] = 2.0E0*I_ERI_Pz_S_Pz_Pz_C1000001_d;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_S_C1001000_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_P_P_C1001000_d
   ************************************************************/
  abcd[388] = 2.0E0*I_ERI_S_Px_Px_Pz_C1001000_d;
  abcd[392] = 2.0E0*I_ERI_S_Py_Px_Pz_C1001000_d;
  abcd[396] = 2.0E0*I_ERI_S_Pz_Px_Pz_C1001000_d;
  abcd[404] = 2.0E0*I_ERI_S_Px_Py_Pz_C1001000_d;
  abcd[408] = 2.0E0*I_ERI_S_Py_Py_Pz_C1001000_d;
  abcd[412] = 2.0E0*I_ERI_S_Pz_Py_Pz_C1001000_d;
  abcd[420] = 2.0E0*I_ERI_S_Px_Pz_Pz_C1001000_d;
  abcd[424] = 2.0E0*I_ERI_S_Py_Pz_Pz_C1001000_d;
  abcd[428] = 2.0E0*I_ERI_S_Pz_Pz_Pz_C1001000_d;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001001_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_P_P_C1001001_d
   ************************************************************/
  abcd[389] = 2.0E0*I_ERI_Px_Px_Px_Pz_C1001001_d;
  abcd[390] = 2.0E0*I_ERI_Py_Px_Px_Pz_C1001001_d;
  abcd[391] = 2.0E0*I_ERI_Pz_Px_Px_Pz_C1001001_d;
  abcd[393] = 2.0E0*I_ERI_Px_Py_Px_Pz_C1001001_d;
  abcd[394] = 2.0E0*I_ERI_Py_Py_Px_Pz_C1001001_d;
  abcd[395] = 2.0E0*I_ERI_Pz_Py_Px_Pz_C1001001_d;
  abcd[397] = 2.0E0*I_ERI_Px_Pz_Px_Pz_C1001001_d;
  abcd[398] = 2.0E0*I_ERI_Py_Pz_Px_Pz_C1001001_d;
  abcd[399] = 2.0E0*I_ERI_Pz_Pz_Px_Pz_C1001001_d;
  abcd[405] = 2.0E0*I_ERI_Px_Px_Py_Pz_C1001001_d;
  abcd[406] = 2.0E0*I_ERI_Py_Px_Py_Pz_C1001001_d;
  abcd[407] = 2.0E0*I_ERI_Pz_Px_Py_Pz_C1001001_d;
  abcd[409] = 2.0E0*I_ERI_Px_Py_Py_Pz_C1001001_d;
  abcd[410] = 2.0E0*I_ERI_Py_Py_Py_Pz_C1001001_d;
  abcd[411] = 2.0E0*I_ERI_Pz_Py_Py_Pz_C1001001_d;
  abcd[413] = 2.0E0*I_ERI_Px_Pz_Py_Pz_C1001001_d;
  abcd[414] = 2.0E0*I_ERI_Py_Pz_Py_Pz_C1001001_d;
  abcd[415] = 2.0E0*I_ERI_Pz_Pz_Py_Pz_C1001001_d;
  abcd[421] = 2.0E0*I_ERI_Px_Px_Pz_Pz_C1001001_d;
  abcd[422] = 2.0E0*I_ERI_Py_Px_Pz_Pz_C1001001_d;
  abcd[423] = 2.0E0*I_ERI_Pz_Px_Pz_Pz_C1001001_d;
  abcd[425] = 2.0E0*I_ERI_Px_Py_Pz_Pz_C1001001_d;
  abcd[426] = 2.0E0*I_ERI_Py_Py_Pz_Pz_C1001001_d;
  abcd[427] = 2.0E0*I_ERI_Pz_Py_Pz_Pz_C1001001_d;
  abcd[429] = 2.0E0*I_ERI_Px_Pz_Pz_Pz_C1001001_d;
  abcd[430] = 2.0E0*I_ERI_Py_Pz_Pz_Pz_C1001001_d;
  abcd[431] = 2.0E0*I_ERI_Pz_Pz_Pz_Pz_C1001001_d;
}
