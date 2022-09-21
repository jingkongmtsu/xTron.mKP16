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

void hgp_os_eri_p_p_p_p_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_D2x_S_S_Px = 0.0E0;
  Double I_ERI_Dxy_S_S_Px = 0.0E0;
  Double I_ERI_Dxz_S_S_Px = 0.0E0;
  Double I_ERI_D2y_S_S_Px = 0.0E0;
  Double I_ERI_Dyz_S_S_Px = 0.0E0;
  Double I_ERI_D2z_S_S_Px = 0.0E0;
  Double I_ERI_D2x_S_S_Py = 0.0E0;
  Double I_ERI_Dxy_S_S_Py = 0.0E0;
  Double I_ERI_Dxz_S_S_Py = 0.0E0;
  Double I_ERI_D2y_S_S_Py = 0.0E0;
  Double I_ERI_Dyz_S_S_Py = 0.0E0;
  Double I_ERI_D2z_S_S_Py = 0.0E0;
  Double I_ERI_D2x_S_S_Pz = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz = 0.0E0;
  Double I_ERI_D2y_S_S_Pz = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz = 0.0E0;
  Double I_ERI_D2z_S_S_Pz = 0.0E0;
  Double I_ERI_Px_S_S_Px = 0.0E0;
  Double I_ERI_Py_S_S_Px = 0.0E0;
  Double I_ERI_Pz_S_S_Px = 0.0E0;
  Double I_ERI_Px_S_S_Py = 0.0E0;
  Double I_ERI_Py_S_S_Py = 0.0E0;
  Double I_ERI_Pz_S_S_Py = 0.0E0;
  Double I_ERI_Px_S_S_Pz = 0.0E0;
  Double I_ERI_Py_S_S_Pz = 0.0E0;
  Double I_ERI_Pz_S_S_Pz = 0.0E0;
  Double I_ERI_S_Px_D2x_S = 0.0E0;
  Double I_ERI_S_Py_D2x_S = 0.0E0;
  Double I_ERI_S_Pz_D2x_S = 0.0E0;
  Double I_ERI_S_Px_Dxy_S = 0.0E0;
  Double I_ERI_S_Py_Dxy_S = 0.0E0;
  Double I_ERI_S_Pz_Dxy_S = 0.0E0;
  Double I_ERI_S_Px_Dxz_S = 0.0E0;
  Double I_ERI_S_Py_Dxz_S = 0.0E0;
  Double I_ERI_S_Pz_Dxz_S = 0.0E0;
  Double I_ERI_S_Px_D2y_S = 0.0E0;
  Double I_ERI_S_Py_D2y_S = 0.0E0;
  Double I_ERI_S_Pz_D2y_S = 0.0E0;
  Double I_ERI_S_Px_Dyz_S = 0.0E0;
  Double I_ERI_S_Py_Dyz_S = 0.0E0;
  Double I_ERI_S_Pz_Dyz_S = 0.0E0;
  Double I_ERI_S_Px_D2z_S = 0.0E0;
  Double I_ERI_S_Py_D2z_S = 0.0E0;
  Double I_ERI_S_Pz_D2z_S = 0.0E0;
  Double I_ERI_S_Px_Px_S = 0.0E0;
  Double I_ERI_S_Py_Px_S = 0.0E0;
  Double I_ERI_S_Pz_Px_S = 0.0E0;
  Double I_ERI_S_Px_Py_S = 0.0E0;
  Double I_ERI_S_Py_Py_S = 0.0E0;
  Double I_ERI_S_Pz_Py_S = 0.0E0;
  Double I_ERI_S_Px_Pz_S = 0.0E0;
  Double I_ERI_S_Py_Pz_S = 0.0E0;
  Double I_ERI_S_Pz_Pz_S = 0.0E0;
  Double I_ERI_Px_S_D2x_S = 0.0E0;
  Double I_ERI_Py_S_D2x_S = 0.0E0;
  Double I_ERI_Pz_S_D2x_S = 0.0E0;
  Double I_ERI_Px_S_Dxy_S = 0.0E0;
  Double I_ERI_Py_S_Dxy_S = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S = 0.0E0;
  Double I_ERI_Px_S_Dxz_S = 0.0E0;
  Double I_ERI_Py_S_Dxz_S = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S = 0.0E0;
  Double I_ERI_Px_S_D2y_S = 0.0E0;
  Double I_ERI_Py_S_D2y_S = 0.0E0;
  Double I_ERI_Pz_S_D2y_S = 0.0E0;
  Double I_ERI_Px_S_Dyz_S = 0.0E0;
  Double I_ERI_Py_S_Dyz_S = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S = 0.0E0;
  Double I_ERI_Px_S_D2z_S = 0.0E0;
  Double I_ERI_Py_S_D2z_S = 0.0E0;
  Double I_ERI_Pz_S_D2z_S = 0.0E0;
  Double I_ERI_Px_S_Px_S = 0.0E0;
  Double I_ERI_Py_S_Px_S = 0.0E0;
  Double I_ERI_Pz_S_Px_S = 0.0E0;
  Double I_ERI_Px_S_Py_S = 0.0E0;
  Double I_ERI_Py_S_Py_S = 0.0E0;
  Double I_ERI_Pz_S_Py_S = 0.0E0;
  Double I_ERI_Px_S_Pz_S = 0.0E0;
  Double I_ERI_Py_S_Pz_S = 0.0E0;
  Double I_ERI_Pz_S_Pz_S = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_a = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_a = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_a = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_a = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_a = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_a = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_a = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_a = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_a = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_F3x_S_Px_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_a = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_a = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_a = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_a = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_a = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_a = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_a = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_a = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_a = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_a = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_a = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_a = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_a = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_a = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_a = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_a = 0.0E0;
  Double I_ERI_D2y_S_Px_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_a = 0.0E0;
  Double I_ERI_D2z_S_Px_S_a = 0.0E0;
  Double I_ERI_D2x_S_Py_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_a = 0.0E0;
  Double I_ERI_D2y_S_Py_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_a = 0.0E0;
  Double I_ERI_D2z_S_Py_S_a = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_a = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_a = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_a = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_a = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_a = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_a = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_c = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_c = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_c = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_c = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_c = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_c = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_c = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_c = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_c = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_c = 0.0E0;
  Double I_ERI_Px_S_F3x_S_c = 0.0E0;
  Double I_ERI_Py_S_F3x_S_c = 0.0E0;
  Double I_ERI_Pz_S_F3x_S_c = 0.0E0;
  Double I_ERI_Px_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Py_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Px_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Py_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Px_S_F3y_S_c = 0.0E0;
  Double I_ERI_Py_S_F3y_S_c = 0.0E0;
  Double I_ERI_Pz_S_F3y_S_c = 0.0E0;
  Double I_ERI_Px_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Py_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Px_S_F3z_S_c = 0.0E0;
  Double I_ERI_Py_S_F3z_S_c = 0.0E0;
  Double I_ERI_Pz_S_F3z_S_c = 0.0E0;
  Double I_ERI_Px_S_D2x_S_c = 0.0E0;
  Double I_ERI_Py_S_D2x_S_c = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_c = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Px_S_D2y_S_c = 0.0E0;
  Double I_ERI_Py_S_D2y_S_c = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_c = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Px_S_D2z_S_c = 0.0E0;
  Double I_ERI_Py_S_D2z_S_c = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_c = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_b = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_b = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_b = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_b = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_b = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_b = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_b = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_b = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_b = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_b = 0.0E0;
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
  Double I_ERI_D2x_S_D2x_S_b = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_b = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_b = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_b = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_b = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_b = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_b = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_b = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_b = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_b = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_b = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_b = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_b = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_b = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_b = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_Px_S_D2x_S_b = 0.0E0;
  Double I_ERI_Py_S_D2x_S_b = 0.0E0;
  Double I_ERI_Pz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Px_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Py_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Pz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Px_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Py_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Pz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Px_S_D2y_S_b = 0.0E0;
  Double I_ERI_Py_S_D2y_S_b = 0.0E0;
  Double I_ERI_Pz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Px_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Py_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Pz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Px_S_D2z_S_b = 0.0E0;
  Double I_ERI_Py_S_D2z_S_b = 0.0E0;
  Double I_ERI_Pz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Px_S_Px_S_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_b = 0.0E0;

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
       * shell quartet name: SQ_ERI_D_S_S_P
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_S_Px += I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px += I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px += I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px += I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px += I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px += I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py += I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py += I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py += I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py += I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py += I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py += I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz += I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz += I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz += I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz += I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz += I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz += I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_P
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_S_Px += I_ERI_Px_S_S_Px_vrr;
      I_ERI_Py_S_S_Px += I_ERI_Py_S_S_Px_vrr;
      I_ERI_Pz_S_S_Px += I_ERI_Pz_S_S_Px_vrr;
      I_ERI_Px_S_S_Py += I_ERI_Px_S_S_Py_vrr;
      I_ERI_Py_S_S_Py += I_ERI_Py_S_S_Py_vrr;
      I_ERI_Pz_S_S_Py += I_ERI_Pz_S_S_Py_vrr;
      I_ERI_Px_S_S_Pz += I_ERI_Px_S_S_Pz_vrr;
      I_ERI_Py_S_S_Pz += I_ERI_Py_S_S_Pz_vrr;
      I_ERI_Pz_S_S_Pz += I_ERI_Pz_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_S_Px_D2x_S += I_ERI_S_Px_D2x_S_vrr;
      I_ERI_S_Py_D2x_S += I_ERI_S_Py_D2x_S_vrr;
      I_ERI_S_Pz_D2x_S += I_ERI_S_Pz_D2x_S_vrr;
      I_ERI_S_Px_Dxy_S += I_ERI_S_Px_Dxy_S_vrr;
      I_ERI_S_Py_Dxy_S += I_ERI_S_Py_Dxy_S_vrr;
      I_ERI_S_Pz_Dxy_S += I_ERI_S_Pz_Dxy_S_vrr;
      I_ERI_S_Px_Dxz_S += I_ERI_S_Px_Dxz_S_vrr;
      I_ERI_S_Py_Dxz_S += I_ERI_S_Py_Dxz_S_vrr;
      I_ERI_S_Pz_Dxz_S += I_ERI_S_Pz_Dxz_S_vrr;
      I_ERI_S_Px_D2y_S += I_ERI_S_Px_D2y_S_vrr;
      I_ERI_S_Py_D2y_S += I_ERI_S_Py_D2y_S_vrr;
      I_ERI_S_Pz_D2y_S += I_ERI_S_Pz_D2y_S_vrr;
      I_ERI_S_Px_Dyz_S += I_ERI_S_Px_Dyz_S_vrr;
      I_ERI_S_Py_Dyz_S += I_ERI_S_Py_Dyz_S_vrr;
      I_ERI_S_Pz_Dyz_S += I_ERI_S_Pz_Dyz_S_vrr;
      I_ERI_S_Px_D2z_S += I_ERI_S_Px_D2z_S_vrr;
      I_ERI_S_Py_D2z_S += I_ERI_S_Py_D2z_S_vrr;
      I_ERI_S_Pz_D2z_S += I_ERI_S_Pz_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_P_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_S_Px_Px_S += I_ERI_S_Px_Px_S_vrr;
      I_ERI_S_Py_Px_S += I_ERI_S_Py_Px_S_vrr;
      I_ERI_S_Pz_Px_S += I_ERI_S_Pz_Px_S_vrr;
      I_ERI_S_Px_Py_S += I_ERI_S_Px_Py_S_vrr;
      I_ERI_S_Py_Py_S += I_ERI_S_Py_Py_S_vrr;
      I_ERI_S_Pz_Py_S += I_ERI_S_Pz_Py_S_vrr;
      I_ERI_S_Px_Pz_S += I_ERI_S_Px_Pz_S_vrr;
      I_ERI_S_Py_Pz_S += I_ERI_S_Py_Pz_S_vrr;
      I_ERI_S_Pz_Pz_S += I_ERI_S_Pz_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_D2x_S += I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S += I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S += I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S += I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S += I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S += I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S += I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S += I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S += I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S += I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S += I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S += I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S += I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S += I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S += I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S += I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S += I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S += I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_Px_S += I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S += I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S += I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S += I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S += I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S += I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S += I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S += I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S += I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_a_coefs = alpha;
      I_ERI_F3x_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_a += SQ_ERI_F_S_D_S_a_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_a_coefs = alpha;
      I_ERI_F3x_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_a += SQ_ERI_F_S_P_S_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_a_coefs = alpha;
      I_ERI_D2x_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_a += SQ_ERI_D_S_D_S_a_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_a_coefs = alpha;
      I_ERI_D2x_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_a += SQ_ERI_D_S_P_S_a_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_c_coefs = gamma;
      I_ERI_D2x_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_c += SQ_ERI_D_S_F_S_c_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_c_coefs = gamma;
      I_ERI_D2x_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_c += SQ_ERI_D_S_D_S_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_F_S_c_coefs = gamma;
      I_ERI_Px_S_F3x_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S_c += SQ_ERI_P_S_F_S_c_coefs*I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_c_coefs = gamma;
      I_ERI_Px_S_D2x_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_c += SQ_ERI_P_S_D_S_c_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_b_coefs = beta;
      I_ERI_F3x_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_b += SQ_ERI_F_S_D_S_b_coefs*I_ERI_F3z_S_D2z_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_b_coefs = beta;
      I_ERI_D2x_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_b += SQ_ERI_D_S_D_S_b_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_b_coefs = beta;
      I_ERI_D2x_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_b += SQ_ERI_D_S_P_S_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_D_S_b_coefs = beta;
      I_ERI_Px_S_D2x_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_D2x_S_vrr;
      I_ERI_Py_S_D2x_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_D2x_S_vrr;
      I_ERI_Pz_S_D2x_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_D2x_S_vrr;
      I_ERI_Px_S_Dxy_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_Dxy_S_vrr;
      I_ERI_Py_S_Dxy_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_Dxy_S_vrr;
      I_ERI_Pz_S_Dxy_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_Dxy_S_vrr;
      I_ERI_Px_S_Dxz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_Dxz_S_vrr;
      I_ERI_Py_S_Dxz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_Dxz_S_vrr;
      I_ERI_Pz_S_Dxz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_Dxz_S_vrr;
      I_ERI_Px_S_D2y_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_D2y_S_vrr;
      I_ERI_Py_S_D2y_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_D2y_S_vrr;
      I_ERI_Pz_S_D2y_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_D2y_S_vrr;
      I_ERI_Px_S_Dyz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_Dyz_S_vrr;
      I_ERI_Py_S_Dyz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_Dyz_S_vrr;
      I_ERI_Pz_S_Dyz_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_Dyz_S_vrr;
      I_ERI_Px_S_D2z_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Px_S_D2z_S_vrr;
      I_ERI_Py_S_D2z_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Py_S_D2z_S_vrr;
      I_ERI_Pz_S_D2z_S_b += SQ_ERI_P_S_D_S_b_coefs*I_ERI_Pz_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_b_coefs = beta;
      I_ERI_Px_S_Px_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_b += SQ_ERI_P_S_P_S_b_coefs*I_ERI_Pz_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_P_S_P_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  Double I_ERI_Px_S_Px_Px = I_ERI_Px_S_D2x_S+CDX*I_ERI_Px_S_Px_S;
  Double I_ERI_Py_S_Px_Px = I_ERI_Py_S_D2x_S+CDX*I_ERI_Py_S_Px_S;
  Double I_ERI_Pz_S_Px_Px = I_ERI_Pz_S_D2x_S+CDX*I_ERI_Pz_S_Px_S;
  Double I_ERI_Px_S_Py_Px = I_ERI_Px_S_Dxy_S+CDX*I_ERI_Px_S_Py_S;
  Double I_ERI_Py_S_Py_Px = I_ERI_Py_S_Dxy_S+CDX*I_ERI_Py_S_Py_S;
  Double I_ERI_Pz_S_Py_Px = I_ERI_Pz_S_Dxy_S+CDX*I_ERI_Pz_S_Py_S;
  Double I_ERI_Px_S_Pz_Px = I_ERI_Px_S_Dxz_S+CDX*I_ERI_Px_S_Pz_S;
  Double I_ERI_Py_S_Pz_Px = I_ERI_Py_S_Dxz_S+CDX*I_ERI_Py_S_Pz_S;
  Double I_ERI_Pz_S_Pz_Px = I_ERI_Pz_S_Dxz_S+CDX*I_ERI_Pz_S_Pz_S;
  Double I_ERI_Px_S_Px_Py = I_ERI_Px_S_Dxy_S+CDY*I_ERI_Px_S_Px_S;
  Double I_ERI_Py_S_Px_Py = I_ERI_Py_S_Dxy_S+CDY*I_ERI_Py_S_Px_S;
  Double I_ERI_Pz_S_Px_Py = I_ERI_Pz_S_Dxy_S+CDY*I_ERI_Pz_S_Px_S;
  Double I_ERI_Px_S_Py_Py = I_ERI_Px_S_D2y_S+CDY*I_ERI_Px_S_Py_S;
  Double I_ERI_Py_S_Py_Py = I_ERI_Py_S_D2y_S+CDY*I_ERI_Py_S_Py_S;
  Double I_ERI_Pz_S_Py_Py = I_ERI_Pz_S_D2y_S+CDY*I_ERI_Pz_S_Py_S;
  Double I_ERI_Px_S_Pz_Py = I_ERI_Px_S_Dyz_S+CDY*I_ERI_Px_S_Pz_S;
  Double I_ERI_Py_S_Pz_Py = I_ERI_Py_S_Dyz_S+CDY*I_ERI_Py_S_Pz_S;
  Double I_ERI_Pz_S_Pz_Py = I_ERI_Pz_S_Dyz_S+CDY*I_ERI_Pz_S_Pz_S;
  Double I_ERI_Px_S_Px_Pz = I_ERI_Px_S_Dxz_S+CDZ*I_ERI_Px_S_Px_S;
  Double I_ERI_Py_S_Px_Pz = I_ERI_Py_S_Dxz_S+CDZ*I_ERI_Py_S_Px_S;
  Double I_ERI_Pz_S_Px_Pz = I_ERI_Pz_S_Dxz_S+CDZ*I_ERI_Pz_S_Px_S;
  Double I_ERI_Px_S_Py_Pz = I_ERI_Px_S_Dyz_S+CDZ*I_ERI_Px_S_Py_S;
  Double I_ERI_Py_S_Py_Pz = I_ERI_Py_S_Dyz_S+CDZ*I_ERI_Py_S_Py_S;
  Double I_ERI_Pz_S_Py_Pz = I_ERI_Pz_S_Dyz_S+CDZ*I_ERI_Pz_S_Py_S;
  Double I_ERI_Px_S_Pz_Pz = I_ERI_Px_S_D2z_S+CDZ*I_ERI_Px_S_Pz_S;
  Double I_ERI_Py_S_Pz_Pz = I_ERI_Py_S_D2z_S+CDZ*I_ERI_Py_S_Pz_S;
  Double I_ERI_Pz_S_Pz_Pz = I_ERI_Pz_S_D2z_S+CDZ*I_ERI_Pz_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_S_P_P_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_S_P_D_S
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  Double I_ERI_S_Px_Px_Px = I_ERI_S_Px_D2x_S+CDX*I_ERI_S_Px_Px_S;
  Double I_ERI_S_Py_Px_Px = I_ERI_S_Py_D2x_S+CDX*I_ERI_S_Py_Px_S;
  Double I_ERI_S_Pz_Px_Px = I_ERI_S_Pz_D2x_S+CDX*I_ERI_S_Pz_Px_S;
  Double I_ERI_S_Px_Py_Px = I_ERI_S_Px_Dxy_S+CDX*I_ERI_S_Px_Py_S;
  Double I_ERI_S_Py_Py_Px = I_ERI_S_Py_Dxy_S+CDX*I_ERI_S_Py_Py_S;
  Double I_ERI_S_Pz_Py_Px = I_ERI_S_Pz_Dxy_S+CDX*I_ERI_S_Pz_Py_S;
  Double I_ERI_S_Px_Pz_Px = I_ERI_S_Px_Dxz_S+CDX*I_ERI_S_Px_Pz_S;
  Double I_ERI_S_Py_Pz_Px = I_ERI_S_Py_Dxz_S+CDX*I_ERI_S_Py_Pz_S;
  Double I_ERI_S_Pz_Pz_Px = I_ERI_S_Pz_Dxz_S+CDX*I_ERI_S_Pz_Pz_S;
  Double I_ERI_S_Px_Px_Py = I_ERI_S_Px_Dxy_S+CDY*I_ERI_S_Px_Px_S;
  Double I_ERI_S_Py_Px_Py = I_ERI_S_Py_Dxy_S+CDY*I_ERI_S_Py_Px_S;
  Double I_ERI_S_Pz_Px_Py = I_ERI_S_Pz_Dxy_S+CDY*I_ERI_S_Pz_Px_S;
  Double I_ERI_S_Px_Py_Py = I_ERI_S_Px_D2y_S+CDY*I_ERI_S_Px_Py_S;
  Double I_ERI_S_Py_Py_Py = I_ERI_S_Py_D2y_S+CDY*I_ERI_S_Py_Py_S;
  Double I_ERI_S_Pz_Py_Py = I_ERI_S_Pz_D2y_S+CDY*I_ERI_S_Pz_Py_S;
  Double I_ERI_S_Px_Pz_Py = I_ERI_S_Px_Dyz_S+CDY*I_ERI_S_Px_Pz_S;
  Double I_ERI_S_Py_Pz_Py = I_ERI_S_Py_Dyz_S+CDY*I_ERI_S_Py_Pz_S;
  Double I_ERI_S_Pz_Pz_Py = I_ERI_S_Pz_Dyz_S+CDY*I_ERI_S_Pz_Pz_S;
  Double I_ERI_S_Px_Px_Pz = I_ERI_S_Px_Dxz_S+CDZ*I_ERI_S_Px_Px_S;
  Double I_ERI_S_Py_Px_Pz = I_ERI_S_Py_Dxz_S+CDZ*I_ERI_S_Py_Px_S;
  Double I_ERI_S_Pz_Px_Pz = I_ERI_S_Pz_Dxz_S+CDZ*I_ERI_S_Pz_Px_S;
  Double I_ERI_S_Px_Py_Pz = I_ERI_S_Px_Dyz_S+CDZ*I_ERI_S_Px_Py_S;
  Double I_ERI_S_Py_Py_Pz = I_ERI_S_Py_Dyz_S+CDZ*I_ERI_S_Py_Py_S;
  Double I_ERI_S_Pz_Py_Pz = I_ERI_S_Pz_Dyz_S+CDZ*I_ERI_S_Pz_Py_S;
  Double I_ERI_S_Px_Pz_Pz = I_ERI_S_Px_D2z_S+CDZ*I_ERI_S_Px_Pz_S;
  Double I_ERI_S_Py_Pz_Pz = I_ERI_S_Py_D2z_S+CDZ*I_ERI_S_Py_Pz_S;
  Double I_ERI_S_Pz_Pz_Pz = I_ERI_S_Pz_D2z_S+CDZ*I_ERI_S_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_a
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_a = I_ERI_D2x_S_D2x_S_a+CDX*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_S_Px_Px_a = I_ERI_Dxy_S_D2x_S_a+CDX*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_S_Px_Px_a = I_ERI_Dxz_S_D2x_S_a+CDX*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_S_Px_Px_a = I_ERI_D2y_S_D2x_S_a+CDX*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_S_Px_Px_a = I_ERI_Dyz_S_D2x_S_a+CDX*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_S_Px_Px_a = I_ERI_D2z_S_D2x_S_a+CDX*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_S_Py_Px_a = I_ERI_D2x_S_Dxy_S_a+CDX*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_S_Py_Px_a = I_ERI_Dxy_S_Dxy_S_a+CDX*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_S_Py_Px_a = I_ERI_Dxz_S_Dxy_S_a+CDX*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_S_Py_Px_a = I_ERI_D2y_S_Dxy_S_a+CDX*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_S_Py_Px_a = I_ERI_Dyz_S_Dxy_S_a+CDX*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_S_Py_Px_a = I_ERI_D2z_S_Dxy_S_a+CDX*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_S_Pz_Px_a = I_ERI_D2x_S_Dxz_S_a+CDX*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_S_Pz_Px_a = I_ERI_Dxy_S_Dxz_S_a+CDX*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_S_Pz_Px_a = I_ERI_Dxz_S_Dxz_S_a+CDX*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_S_Pz_Px_a = I_ERI_D2y_S_Dxz_S_a+CDX*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_S_Pz_Px_a = I_ERI_Dyz_S_Dxz_S_a+CDX*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_S_Pz_Px_a = I_ERI_D2z_S_Dxz_S_a+CDX*I_ERI_D2z_S_Pz_S_a;
  Double I_ERI_D2x_S_Px_Py_a = I_ERI_D2x_S_Dxy_S_a+CDY*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_S_Px_Py_a = I_ERI_Dxy_S_Dxy_S_a+CDY*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_S_Px_Py_a = I_ERI_Dxz_S_Dxy_S_a+CDY*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_S_Px_Py_a = I_ERI_D2y_S_Dxy_S_a+CDY*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_S_Px_Py_a = I_ERI_Dyz_S_Dxy_S_a+CDY*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_S_Px_Py_a = I_ERI_D2z_S_Dxy_S_a+CDY*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_S_Py_Py_a = I_ERI_D2x_S_D2y_S_a+CDY*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_S_Py_Py_a = I_ERI_Dxy_S_D2y_S_a+CDY*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_S_Py_Py_a = I_ERI_Dxz_S_D2y_S_a+CDY*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_S_Py_Py_a = I_ERI_D2y_S_D2y_S_a+CDY*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_S_Py_Py_a = I_ERI_Dyz_S_D2y_S_a+CDY*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_S_Py_Py_a = I_ERI_D2z_S_D2y_S_a+CDY*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_S_Pz_Py_a = I_ERI_D2x_S_Dyz_S_a+CDY*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_S_Pz_Py_a = I_ERI_Dxy_S_Dyz_S_a+CDY*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_S_Pz_Py_a = I_ERI_Dxz_S_Dyz_S_a+CDY*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_S_Pz_Py_a = I_ERI_D2y_S_Dyz_S_a+CDY*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_S_Pz_Py_a = I_ERI_Dyz_S_Dyz_S_a+CDY*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_S_Pz_Py_a = I_ERI_D2z_S_Dyz_S_a+CDY*I_ERI_D2z_S_Pz_S_a;
  Double I_ERI_D2x_S_Px_Pz_a = I_ERI_D2x_S_Dxz_S_a+CDZ*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_S_Px_Pz_a = I_ERI_Dxy_S_Dxz_S_a+CDZ*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_S_Px_Pz_a = I_ERI_Dxz_S_Dxz_S_a+CDZ*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_S_Px_Pz_a = I_ERI_D2y_S_Dxz_S_a+CDZ*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_S_Px_Pz_a = I_ERI_Dyz_S_Dxz_S_a+CDZ*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_S_Px_Pz_a = I_ERI_D2z_S_Dxz_S_a+CDZ*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_S_Py_Pz_a = I_ERI_D2x_S_Dyz_S_a+CDZ*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_S_Py_Pz_a = I_ERI_Dxy_S_Dyz_S_a+CDZ*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_S_Py_Pz_a = I_ERI_Dxz_S_Dyz_S_a+CDZ*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_S_Py_Pz_a = I_ERI_D2y_S_Dyz_S_a+CDZ*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_S_Py_Pz_a = I_ERI_Dyz_S_Dyz_S_a+CDZ*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_S_Py_Pz_a = I_ERI_D2z_S_Dyz_S_a+CDZ*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_S_Pz_Pz_a = I_ERI_D2x_S_D2z_S_a+CDZ*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_S_Pz_Pz_a = I_ERI_Dxy_S_D2z_S_a+CDZ*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_S_Pz_Pz_a = I_ERI_Dxz_S_D2z_S_a+CDZ*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_S_Pz_Pz_a = I_ERI_D2y_S_D2z_S_a+CDZ*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_S_Pz_Pz_a = I_ERI_Dyz_S_D2z_S_a+CDZ*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_S_Pz_Pz_a = I_ERI_D2z_S_D2z_S_a+CDZ*I_ERI_D2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_a
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_a = I_ERI_F3x_S_D2x_S_a+CDX*I_ERI_F3x_S_Px_S_a;
  Double I_ERI_F2xy_S_Px_Px_a = I_ERI_F2xy_S_D2x_S_a+CDX*I_ERI_F2xy_S_Px_S_a;
  Double I_ERI_F2xz_S_Px_Px_a = I_ERI_F2xz_S_D2x_S_a+CDX*I_ERI_F2xz_S_Px_S_a;
  Double I_ERI_Fx2y_S_Px_Px_a = I_ERI_Fx2y_S_D2x_S_a+CDX*I_ERI_Fx2y_S_Px_S_a;
  Double I_ERI_Fxyz_S_Px_Px_a = I_ERI_Fxyz_S_D2x_S_a+CDX*I_ERI_Fxyz_S_Px_S_a;
  Double I_ERI_Fx2z_S_Px_Px_a = I_ERI_Fx2z_S_D2x_S_a+CDX*I_ERI_Fx2z_S_Px_S_a;
  Double I_ERI_F3y_S_Px_Px_a = I_ERI_F3y_S_D2x_S_a+CDX*I_ERI_F3y_S_Px_S_a;
  Double I_ERI_F2yz_S_Px_Px_a = I_ERI_F2yz_S_D2x_S_a+CDX*I_ERI_F2yz_S_Px_S_a;
  Double I_ERI_Fy2z_S_Px_Px_a = I_ERI_Fy2z_S_D2x_S_a+CDX*I_ERI_Fy2z_S_Px_S_a;
  Double I_ERI_F3z_S_Px_Px_a = I_ERI_F3z_S_D2x_S_a+CDX*I_ERI_F3z_S_Px_S_a;
  Double I_ERI_F3x_S_Py_Px_a = I_ERI_F3x_S_Dxy_S_a+CDX*I_ERI_F3x_S_Py_S_a;
  Double I_ERI_F2xy_S_Py_Px_a = I_ERI_F2xy_S_Dxy_S_a+CDX*I_ERI_F2xy_S_Py_S_a;
  Double I_ERI_F2xz_S_Py_Px_a = I_ERI_F2xz_S_Dxy_S_a+CDX*I_ERI_F2xz_S_Py_S_a;
  Double I_ERI_Fx2y_S_Py_Px_a = I_ERI_Fx2y_S_Dxy_S_a+CDX*I_ERI_Fx2y_S_Py_S_a;
  Double I_ERI_Fxyz_S_Py_Px_a = I_ERI_Fxyz_S_Dxy_S_a+CDX*I_ERI_Fxyz_S_Py_S_a;
  Double I_ERI_Fx2z_S_Py_Px_a = I_ERI_Fx2z_S_Dxy_S_a+CDX*I_ERI_Fx2z_S_Py_S_a;
  Double I_ERI_F3y_S_Py_Px_a = I_ERI_F3y_S_Dxy_S_a+CDX*I_ERI_F3y_S_Py_S_a;
  Double I_ERI_F2yz_S_Py_Px_a = I_ERI_F2yz_S_Dxy_S_a+CDX*I_ERI_F2yz_S_Py_S_a;
  Double I_ERI_Fy2z_S_Py_Px_a = I_ERI_Fy2z_S_Dxy_S_a+CDX*I_ERI_Fy2z_S_Py_S_a;
  Double I_ERI_F3z_S_Py_Px_a = I_ERI_F3z_S_Dxy_S_a+CDX*I_ERI_F3z_S_Py_S_a;
  Double I_ERI_F3x_S_Pz_Px_a = I_ERI_F3x_S_Dxz_S_a+CDX*I_ERI_F3x_S_Pz_S_a;
  Double I_ERI_F2xy_S_Pz_Px_a = I_ERI_F2xy_S_Dxz_S_a+CDX*I_ERI_F2xy_S_Pz_S_a;
  Double I_ERI_F2xz_S_Pz_Px_a = I_ERI_F2xz_S_Dxz_S_a+CDX*I_ERI_F2xz_S_Pz_S_a;
  Double I_ERI_Fx2y_S_Pz_Px_a = I_ERI_Fx2y_S_Dxz_S_a+CDX*I_ERI_Fx2y_S_Pz_S_a;
  Double I_ERI_Fxyz_S_Pz_Px_a = I_ERI_Fxyz_S_Dxz_S_a+CDX*I_ERI_Fxyz_S_Pz_S_a;
  Double I_ERI_Fx2z_S_Pz_Px_a = I_ERI_Fx2z_S_Dxz_S_a+CDX*I_ERI_Fx2z_S_Pz_S_a;
  Double I_ERI_F3y_S_Pz_Px_a = I_ERI_F3y_S_Dxz_S_a+CDX*I_ERI_F3y_S_Pz_S_a;
  Double I_ERI_F2yz_S_Pz_Px_a = I_ERI_F2yz_S_Dxz_S_a+CDX*I_ERI_F2yz_S_Pz_S_a;
  Double I_ERI_Fy2z_S_Pz_Px_a = I_ERI_Fy2z_S_Dxz_S_a+CDX*I_ERI_Fy2z_S_Pz_S_a;
  Double I_ERI_F3z_S_Pz_Px_a = I_ERI_F3z_S_Dxz_S_a+CDX*I_ERI_F3z_S_Pz_S_a;
  Double I_ERI_F3x_S_Px_Py_a = I_ERI_F3x_S_Dxy_S_a+CDY*I_ERI_F3x_S_Px_S_a;
  Double I_ERI_F2xy_S_Px_Py_a = I_ERI_F2xy_S_Dxy_S_a+CDY*I_ERI_F2xy_S_Px_S_a;
  Double I_ERI_F2xz_S_Px_Py_a = I_ERI_F2xz_S_Dxy_S_a+CDY*I_ERI_F2xz_S_Px_S_a;
  Double I_ERI_Fx2y_S_Px_Py_a = I_ERI_Fx2y_S_Dxy_S_a+CDY*I_ERI_Fx2y_S_Px_S_a;
  Double I_ERI_Fxyz_S_Px_Py_a = I_ERI_Fxyz_S_Dxy_S_a+CDY*I_ERI_Fxyz_S_Px_S_a;
  Double I_ERI_Fx2z_S_Px_Py_a = I_ERI_Fx2z_S_Dxy_S_a+CDY*I_ERI_Fx2z_S_Px_S_a;
  Double I_ERI_F3y_S_Px_Py_a = I_ERI_F3y_S_Dxy_S_a+CDY*I_ERI_F3y_S_Px_S_a;
  Double I_ERI_F2yz_S_Px_Py_a = I_ERI_F2yz_S_Dxy_S_a+CDY*I_ERI_F2yz_S_Px_S_a;
  Double I_ERI_Fy2z_S_Px_Py_a = I_ERI_Fy2z_S_Dxy_S_a+CDY*I_ERI_Fy2z_S_Px_S_a;
  Double I_ERI_F3z_S_Px_Py_a = I_ERI_F3z_S_Dxy_S_a+CDY*I_ERI_F3z_S_Px_S_a;
  Double I_ERI_F3x_S_Py_Py_a = I_ERI_F3x_S_D2y_S_a+CDY*I_ERI_F3x_S_Py_S_a;
  Double I_ERI_F2xy_S_Py_Py_a = I_ERI_F2xy_S_D2y_S_a+CDY*I_ERI_F2xy_S_Py_S_a;
  Double I_ERI_F2xz_S_Py_Py_a = I_ERI_F2xz_S_D2y_S_a+CDY*I_ERI_F2xz_S_Py_S_a;
  Double I_ERI_Fx2y_S_Py_Py_a = I_ERI_Fx2y_S_D2y_S_a+CDY*I_ERI_Fx2y_S_Py_S_a;
  Double I_ERI_Fxyz_S_Py_Py_a = I_ERI_Fxyz_S_D2y_S_a+CDY*I_ERI_Fxyz_S_Py_S_a;
  Double I_ERI_Fx2z_S_Py_Py_a = I_ERI_Fx2z_S_D2y_S_a+CDY*I_ERI_Fx2z_S_Py_S_a;
  Double I_ERI_F3y_S_Py_Py_a = I_ERI_F3y_S_D2y_S_a+CDY*I_ERI_F3y_S_Py_S_a;
  Double I_ERI_F2yz_S_Py_Py_a = I_ERI_F2yz_S_D2y_S_a+CDY*I_ERI_F2yz_S_Py_S_a;
  Double I_ERI_Fy2z_S_Py_Py_a = I_ERI_Fy2z_S_D2y_S_a+CDY*I_ERI_Fy2z_S_Py_S_a;
  Double I_ERI_F3z_S_Py_Py_a = I_ERI_F3z_S_D2y_S_a+CDY*I_ERI_F3z_S_Py_S_a;
  Double I_ERI_F3x_S_Pz_Py_a = I_ERI_F3x_S_Dyz_S_a+CDY*I_ERI_F3x_S_Pz_S_a;
  Double I_ERI_F2xy_S_Pz_Py_a = I_ERI_F2xy_S_Dyz_S_a+CDY*I_ERI_F2xy_S_Pz_S_a;
  Double I_ERI_F2xz_S_Pz_Py_a = I_ERI_F2xz_S_Dyz_S_a+CDY*I_ERI_F2xz_S_Pz_S_a;
  Double I_ERI_Fx2y_S_Pz_Py_a = I_ERI_Fx2y_S_Dyz_S_a+CDY*I_ERI_Fx2y_S_Pz_S_a;
  Double I_ERI_Fxyz_S_Pz_Py_a = I_ERI_Fxyz_S_Dyz_S_a+CDY*I_ERI_Fxyz_S_Pz_S_a;
  Double I_ERI_Fx2z_S_Pz_Py_a = I_ERI_Fx2z_S_Dyz_S_a+CDY*I_ERI_Fx2z_S_Pz_S_a;
  Double I_ERI_F3y_S_Pz_Py_a = I_ERI_F3y_S_Dyz_S_a+CDY*I_ERI_F3y_S_Pz_S_a;
  Double I_ERI_F2yz_S_Pz_Py_a = I_ERI_F2yz_S_Dyz_S_a+CDY*I_ERI_F2yz_S_Pz_S_a;
  Double I_ERI_Fy2z_S_Pz_Py_a = I_ERI_Fy2z_S_Dyz_S_a+CDY*I_ERI_Fy2z_S_Pz_S_a;
  Double I_ERI_F3z_S_Pz_Py_a = I_ERI_F3z_S_Dyz_S_a+CDY*I_ERI_F3z_S_Pz_S_a;
  Double I_ERI_F3x_S_Px_Pz_a = I_ERI_F3x_S_Dxz_S_a+CDZ*I_ERI_F3x_S_Px_S_a;
  Double I_ERI_F2xy_S_Px_Pz_a = I_ERI_F2xy_S_Dxz_S_a+CDZ*I_ERI_F2xy_S_Px_S_a;
  Double I_ERI_F2xz_S_Px_Pz_a = I_ERI_F2xz_S_Dxz_S_a+CDZ*I_ERI_F2xz_S_Px_S_a;
  Double I_ERI_Fx2y_S_Px_Pz_a = I_ERI_Fx2y_S_Dxz_S_a+CDZ*I_ERI_Fx2y_S_Px_S_a;
  Double I_ERI_Fxyz_S_Px_Pz_a = I_ERI_Fxyz_S_Dxz_S_a+CDZ*I_ERI_Fxyz_S_Px_S_a;
  Double I_ERI_Fx2z_S_Px_Pz_a = I_ERI_Fx2z_S_Dxz_S_a+CDZ*I_ERI_Fx2z_S_Px_S_a;
  Double I_ERI_F3y_S_Px_Pz_a = I_ERI_F3y_S_Dxz_S_a+CDZ*I_ERI_F3y_S_Px_S_a;
  Double I_ERI_F2yz_S_Px_Pz_a = I_ERI_F2yz_S_Dxz_S_a+CDZ*I_ERI_F2yz_S_Px_S_a;
  Double I_ERI_Fy2z_S_Px_Pz_a = I_ERI_Fy2z_S_Dxz_S_a+CDZ*I_ERI_Fy2z_S_Px_S_a;
  Double I_ERI_F3z_S_Px_Pz_a = I_ERI_F3z_S_Dxz_S_a+CDZ*I_ERI_F3z_S_Px_S_a;
  Double I_ERI_F3x_S_Py_Pz_a = I_ERI_F3x_S_Dyz_S_a+CDZ*I_ERI_F3x_S_Py_S_a;
  Double I_ERI_F2xy_S_Py_Pz_a = I_ERI_F2xy_S_Dyz_S_a+CDZ*I_ERI_F2xy_S_Py_S_a;
  Double I_ERI_F2xz_S_Py_Pz_a = I_ERI_F2xz_S_Dyz_S_a+CDZ*I_ERI_F2xz_S_Py_S_a;
  Double I_ERI_Fx2y_S_Py_Pz_a = I_ERI_Fx2y_S_Dyz_S_a+CDZ*I_ERI_Fx2y_S_Py_S_a;
  Double I_ERI_Fxyz_S_Py_Pz_a = I_ERI_Fxyz_S_Dyz_S_a+CDZ*I_ERI_Fxyz_S_Py_S_a;
  Double I_ERI_Fx2z_S_Py_Pz_a = I_ERI_Fx2z_S_Dyz_S_a+CDZ*I_ERI_Fx2z_S_Py_S_a;
  Double I_ERI_F3y_S_Py_Pz_a = I_ERI_F3y_S_Dyz_S_a+CDZ*I_ERI_F3y_S_Py_S_a;
  Double I_ERI_F2yz_S_Py_Pz_a = I_ERI_F2yz_S_Dyz_S_a+CDZ*I_ERI_F2yz_S_Py_S_a;
  Double I_ERI_Fy2z_S_Py_Pz_a = I_ERI_Fy2z_S_Dyz_S_a+CDZ*I_ERI_Fy2z_S_Py_S_a;
  Double I_ERI_F3z_S_Py_Pz_a = I_ERI_F3z_S_Dyz_S_a+CDZ*I_ERI_F3z_S_Py_S_a;
  Double I_ERI_F3x_S_Pz_Pz_a = I_ERI_F3x_S_D2z_S_a+CDZ*I_ERI_F3x_S_Pz_S_a;
  Double I_ERI_F2xy_S_Pz_Pz_a = I_ERI_F2xy_S_D2z_S_a+CDZ*I_ERI_F2xy_S_Pz_S_a;
  Double I_ERI_F2xz_S_Pz_Pz_a = I_ERI_F2xz_S_D2z_S_a+CDZ*I_ERI_F2xz_S_Pz_S_a;
  Double I_ERI_Fx2y_S_Pz_Pz_a = I_ERI_Fx2y_S_D2z_S_a+CDZ*I_ERI_Fx2y_S_Pz_S_a;
  Double I_ERI_Fxyz_S_Pz_Pz_a = I_ERI_Fxyz_S_D2z_S_a+CDZ*I_ERI_Fxyz_S_Pz_S_a;
  Double I_ERI_Fx2z_S_Pz_Pz_a = I_ERI_Fx2z_S_D2z_S_a+CDZ*I_ERI_Fx2z_S_Pz_S_a;
  Double I_ERI_F3y_S_Pz_Pz_a = I_ERI_F3y_S_D2z_S_a+CDZ*I_ERI_F3y_S_Pz_S_a;
  Double I_ERI_F2yz_S_Pz_Pz_a = I_ERI_F2yz_S_D2z_S_a+CDZ*I_ERI_F2yz_S_Pz_S_a;
  Double I_ERI_Fy2z_S_Pz_Pz_a = I_ERI_Fy2z_S_D2z_S_a+CDZ*I_ERI_Fy2z_S_Pz_S_a;
  Double I_ERI_F3z_S_Pz_Pz_a = I_ERI_F3z_S_D2z_S_a+CDZ*I_ERI_F3z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_P_P_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_D_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_b
   ************************************************************/
  Double I_ERI_Px_S_Px_Px_b = I_ERI_Px_S_D2x_S_b+CDX*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_S_Px_Px_b = I_ERI_Py_S_D2x_S_b+CDX*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_S_Px_Px_b = I_ERI_Pz_S_D2x_S_b+CDX*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_S_Py_Px_b = I_ERI_Px_S_Dxy_S_b+CDX*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_S_Py_Px_b = I_ERI_Py_S_Dxy_S_b+CDX*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_S_Py_Px_b = I_ERI_Pz_S_Dxy_S_b+CDX*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_S_Pz_Px_b = I_ERI_Px_S_Dxz_S_b+CDX*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_S_Pz_Px_b = I_ERI_Py_S_Dxz_S_b+CDX*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_S_Pz_Px_b = I_ERI_Pz_S_Dxz_S_b+CDX*I_ERI_Pz_S_Pz_S_b;
  Double I_ERI_Px_S_Px_Py_b = I_ERI_Px_S_Dxy_S_b+CDY*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_S_Px_Py_b = I_ERI_Py_S_Dxy_S_b+CDY*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_S_Px_Py_b = I_ERI_Pz_S_Dxy_S_b+CDY*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_S_Py_Py_b = I_ERI_Px_S_D2y_S_b+CDY*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_S_Py_Py_b = I_ERI_Py_S_D2y_S_b+CDY*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_S_Py_Py_b = I_ERI_Pz_S_D2y_S_b+CDY*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_S_Pz_Py_b = I_ERI_Px_S_Dyz_S_b+CDY*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_S_Pz_Py_b = I_ERI_Py_S_Dyz_S_b+CDY*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_S_Pz_Py_b = I_ERI_Pz_S_Dyz_S_b+CDY*I_ERI_Pz_S_Pz_S_b;
  Double I_ERI_Px_S_Px_Pz_b = I_ERI_Px_S_Dxz_S_b+CDZ*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_S_Px_Pz_b = I_ERI_Py_S_Dxz_S_b+CDZ*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_S_Px_Pz_b = I_ERI_Pz_S_Dxz_S_b+CDZ*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_S_Py_Pz_b = I_ERI_Px_S_Dyz_S_b+CDZ*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_S_Py_Pz_b = I_ERI_Py_S_Dyz_S_b+CDZ*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_S_Py_Pz_b = I_ERI_Pz_S_Dyz_S_b+CDZ*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_S_Pz_Pz_b = I_ERI_Px_S_D2z_S_b+CDZ*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_S_Pz_Pz_b = I_ERI_Py_S_D2z_S_b+CDZ*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_S_Pz_Pz_b = I_ERI_Pz_S_D2z_S_b+CDZ*I_ERI_Pz_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_P_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_b
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_b = I_ERI_D2x_S_D2x_S_b+CDX*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_S_Px_Px_b = I_ERI_Dxy_S_D2x_S_b+CDX*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_S_Px_Px_b = I_ERI_Dxz_S_D2x_S_b+CDX*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_S_Px_Px_b = I_ERI_D2y_S_D2x_S_b+CDX*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_S_Px_Px_b = I_ERI_Dyz_S_D2x_S_b+CDX*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_S_Px_Px_b = I_ERI_D2z_S_D2x_S_b+CDX*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_S_Py_Px_b = I_ERI_D2x_S_Dxy_S_b+CDX*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_S_Py_Px_b = I_ERI_Dxy_S_Dxy_S_b+CDX*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_S_Py_Px_b = I_ERI_Dxz_S_Dxy_S_b+CDX*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_S_Py_Px_b = I_ERI_D2y_S_Dxy_S_b+CDX*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_S_Py_Px_b = I_ERI_Dyz_S_Dxy_S_b+CDX*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_S_Py_Px_b = I_ERI_D2z_S_Dxy_S_b+CDX*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_S_Pz_Px_b = I_ERI_D2x_S_Dxz_S_b+CDX*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_S_Pz_Px_b = I_ERI_Dxy_S_Dxz_S_b+CDX*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_S_Pz_Px_b = I_ERI_Dxz_S_Dxz_S_b+CDX*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_S_Pz_Px_b = I_ERI_D2y_S_Dxz_S_b+CDX*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_S_Pz_Px_b = I_ERI_Dyz_S_Dxz_S_b+CDX*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_S_Pz_Px_b = I_ERI_D2z_S_Dxz_S_b+CDX*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_S_Px_Py_b = I_ERI_D2x_S_Dxy_S_b+CDY*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_S_Px_Py_b = I_ERI_Dxy_S_Dxy_S_b+CDY*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_S_Px_Py_b = I_ERI_Dxz_S_Dxy_S_b+CDY*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_S_Px_Py_b = I_ERI_D2y_S_Dxy_S_b+CDY*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_S_Px_Py_b = I_ERI_Dyz_S_Dxy_S_b+CDY*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_S_Px_Py_b = I_ERI_D2z_S_Dxy_S_b+CDY*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_S_Py_Py_b = I_ERI_D2x_S_D2y_S_b+CDY*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_S_Py_Py_b = I_ERI_Dxy_S_D2y_S_b+CDY*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_S_Py_Py_b = I_ERI_Dxz_S_D2y_S_b+CDY*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_S_Py_Py_b = I_ERI_D2y_S_D2y_S_b+CDY*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_S_Py_Py_b = I_ERI_Dyz_S_D2y_S_b+CDY*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_S_Py_Py_b = I_ERI_D2z_S_D2y_S_b+CDY*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_S_Pz_Py_b = I_ERI_D2x_S_Dyz_S_b+CDY*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_S_Pz_Py_b = I_ERI_Dxy_S_Dyz_S_b+CDY*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_S_Pz_Py_b = I_ERI_Dxz_S_Dyz_S_b+CDY*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_S_Pz_Py_b = I_ERI_D2y_S_Dyz_S_b+CDY*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_S_Pz_Py_b = I_ERI_Dyz_S_Dyz_S_b+CDY*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_S_Pz_Py_b = I_ERI_D2z_S_Dyz_S_b+CDY*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_S_Px_Pz_b = I_ERI_D2x_S_Dxz_S_b+CDZ*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_S_Px_Pz_b = I_ERI_Dxy_S_Dxz_S_b+CDZ*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_S_Px_Pz_b = I_ERI_Dxz_S_Dxz_S_b+CDZ*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_S_Px_Pz_b = I_ERI_D2y_S_Dxz_S_b+CDZ*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_S_Px_Pz_b = I_ERI_Dyz_S_Dxz_S_b+CDZ*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_S_Px_Pz_b = I_ERI_D2z_S_Dxz_S_b+CDZ*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_S_Py_Pz_b = I_ERI_D2x_S_Dyz_S_b+CDZ*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_S_Py_Pz_b = I_ERI_Dxy_S_Dyz_S_b+CDZ*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_S_Py_Pz_b = I_ERI_Dxz_S_Dyz_S_b+CDZ*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_S_Py_Pz_b = I_ERI_D2y_S_Dyz_S_b+CDZ*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_S_Py_Pz_b = I_ERI_Dyz_S_Dyz_S_b+CDZ*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_S_Py_Pz_b = I_ERI_D2z_S_Dyz_S_b+CDZ*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_S_Pz_Pz_b = I_ERI_D2x_S_D2z_S_b+CDZ*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_S_Pz_Pz_b = I_ERI_Dxy_S_D2z_S_b+CDZ*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_S_Pz_Pz_b = I_ERI_Dxz_S_D2z_S_b+CDZ*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_S_Pz_Pz_b = I_ERI_D2y_S_D2z_S_b+CDZ*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_S_Pz_Pz_b = I_ERI_Dyz_S_D2z_S_b+CDZ*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_S_Pz_Pz_b = I_ERI_D2z_S_D2z_S_b+CDZ*I_ERI_D2z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_b = I_ERI_F3x_S_D2x_S_b+CDX*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_S_Px_Px_b = I_ERI_F2xy_S_D2x_S_b+CDX*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_S_Px_Px_b = I_ERI_F2xz_S_D2x_S_b+CDX*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_S_Px_Px_b = I_ERI_Fx2y_S_D2x_S_b+CDX*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_S_Px_Px_b = I_ERI_Fxyz_S_D2x_S_b+CDX*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_S_Px_Px_b = I_ERI_Fx2z_S_D2x_S_b+CDX*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_S_Px_Px_b = I_ERI_F3y_S_D2x_S_b+CDX*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_S_Px_Px_b = I_ERI_F2yz_S_D2x_S_b+CDX*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_S_Px_Px_b = I_ERI_Fy2z_S_D2x_S_b+CDX*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_S_Px_Px_b = I_ERI_F3z_S_D2x_S_b+CDX*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_S_Py_Px_b = I_ERI_F3x_S_Dxy_S_b+CDX*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_S_Py_Px_b = I_ERI_F2xy_S_Dxy_S_b+CDX*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_S_Py_Px_b = I_ERI_F2xz_S_Dxy_S_b+CDX*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_S_Py_Px_b = I_ERI_Fx2y_S_Dxy_S_b+CDX*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_S_Py_Px_b = I_ERI_Fxyz_S_Dxy_S_b+CDX*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_S_Py_Px_b = I_ERI_Fx2z_S_Dxy_S_b+CDX*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_S_Py_Px_b = I_ERI_F3y_S_Dxy_S_b+CDX*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_S_Py_Px_b = I_ERI_F2yz_S_Dxy_S_b+CDX*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_S_Py_Px_b = I_ERI_Fy2z_S_Dxy_S_b+CDX*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_S_Py_Px_b = I_ERI_F3z_S_Dxy_S_b+CDX*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_S_Pz_Px_b = I_ERI_F3x_S_Dxz_S_b+CDX*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_S_Pz_Px_b = I_ERI_F2xy_S_Dxz_S_b+CDX*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_S_Pz_Px_b = I_ERI_F2xz_S_Dxz_S_b+CDX*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_S_Pz_Px_b = I_ERI_Fx2y_S_Dxz_S_b+CDX*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_S_Pz_Px_b = I_ERI_Fxyz_S_Dxz_S_b+CDX*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_S_Pz_Px_b = I_ERI_Fx2z_S_Dxz_S_b+CDX*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_S_Pz_Px_b = I_ERI_F3y_S_Dxz_S_b+CDX*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_S_Pz_Px_b = I_ERI_F2yz_S_Dxz_S_b+CDX*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_S_Pz_Px_b = I_ERI_Fy2z_S_Dxz_S_b+CDX*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_S_Pz_Px_b = I_ERI_F3z_S_Dxz_S_b+CDX*I_ERI_F3z_S_Pz_S_b;
  Double I_ERI_F3x_S_Px_Py_b = I_ERI_F3x_S_Dxy_S_b+CDY*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_S_Px_Py_b = I_ERI_F2xy_S_Dxy_S_b+CDY*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_S_Px_Py_b = I_ERI_F2xz_S_Dxy_S_b+CDY*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_S_Px_Py_b = I_ERI_Fx2y_S_Dxy_S_b+CDY*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_S_Px_Py_b = I_ERI_Fxyz_S_Dxy_S_b+CDY*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_S_Px_Py_b = I_ERI_Fx2z_S_Dxy_S_b+CDY*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_S_Px_Py_b = I_ERI_F3y_S_Dxy_S_b+CDY*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_S_Px_Py_b = I_ERI_F2yz_S_Dxy_S_b+CDY*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_S_Px_Py_b = I_ERI_Fy2z_S_Dxy_S_b+CDY*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_S_Px_Py_b = I_ERI_F3z_S_Dxy_S_b+CDY*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_S_Py_Py_b = I_ERI_F3x_S_D2y_S_b+CDY*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_S_Py_Py_b = I_ERI_F2xy_S_D2y_S_b+CDY*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_S_Py_Py_b = I_ERI_F2xz_S_D2y_S_b+CDY*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_S_Py_Py_b = I_ERI_Fx2y_S_D2y_S_b+CDY*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_S_Py_Py_b = I_ERI_Fxyz_S_D2y_S_b+CDY*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_S_Py_Py_b = I_ERI_Fx2z_S_D2y_S_b+CDY*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_S_Py_Py_b = I_ERI_F3y_S_D2y_S_b+CDY*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_S_Py_Py_b = I_ERI_F2yz_S_D2y_S_b+CDY*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_S_Py_Py_b = I_ERI_Fy2z_S_D2y_S_b+CDY*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_S_Py_Py_b = I_ERI_F3z_S_D2y_S_b+CDY*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_S_Pz_Py_b = I_ERI_F3x_S_Dyz_S_b+CDY*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_S_Pz_Py_b = I_ERI_F2xy_S_Dyz_S_b+CDY*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_S_Pz_Py_b = I_ERI_F2xz_S_Dyz_S_b+CDY*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_S_Pz_Py_b = I_ERI_Fx2y_S_Dyz_S_b+CDY*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_S_Pz_Py_b = I_ERI_Fxyz_S_Dyz_S_b+CDY*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_S_Pz_Py_b = I_ERI_Fx2z_S_Dyz_S_b+CDY*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_S_Pz_Py_b = I_ERI_F3y_S_Dyz_S_b+CDY*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_S_Pz_Py_b = I_ERI_F2yz_S_Dyz_S_b+CDY*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_S_Pz_Py_b = I_ERI_Fy2z_S_Dyz_S_b+CDY*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_S_Pz_Py_b = I_ERI_F3z_S_Dyz_S_b+CDY*I_ERI_F3z_S_Pz_S_b;
  Double I_ERI_F3x_S_Px_Pz_b = I_ERI_F3x_S_Dxz_S_b+CDZ*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_S_Px_Pz_b = I_ERI_F2xy_S_Dxz_S_b+CDZ*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_S_Px_Pz_b = I_ERI_F2xz_S_Dxz_S_b+CDZ*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_S_Px_Pz_b = I_ERI_Fx2y_S_Dxz_S_b+CDZ*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_S_Px_Pz_b = I_ERI_Fxyz_S_Dxz_S_b+CDZ*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_S_Px_Pz_b = I_ERI_Fx2z_S_Dxz_S_b+CDZ*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_S_Px_Pz_b = I_ERI_F3y_S_Dxz_S_b+CDZ*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_S_Px_Pz_b = I_ERI_F2yz_S_Dxz_S_b+CDZ*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_S_Px_Pz_b = I_ERI_Fy2z_S_Dxz_S_b+CDZ*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_S_Px_Pz_b = I_ERI_F3z_S_Dxz_S_b+CDZ*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_S_Py_Pz_b = I_ERI_F3x_S_Dyz_S_b+CDZ*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_S_Py_Pz_b = I_ERI_F2xy_S_Dyz_S_b+CDZ*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_S_Py_Pz_b = I_ERI_F2xz_S_Dyz_S_b+CDZ*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_S_Py_Pz_b = I_ERI_Fx2y_S_Dyz_S_b+CDZ*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_S_Py_Pz_b = I_ERI_Fxyz_S_Dyz_S_b+CDZ*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_S_Py_Pz_b = I_ERI_Fx2z_S_Dyz_S_b+CDZ*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_S_Py_Pz_b = I_ERI_F3y_S_Dyz_S_b+CDZ*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_S_Py_Pz_b = I_ERI_F2yz_S_Dyz_S_b+CDZ*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_S_Py_Pz_b = I_ERI_Fy2z_S_Dyz_S_b+CDZ*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_S_Py_Pz_b = I_ERI_F3z_S_Dyz_S_b+CDZ*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_S_Pz_Pz_b = I_ERI_F3x_S_D2z_S_b+CDZ*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_S_Pz_Pz_b = I_ERI_F2xy_S_D2z_S_b+CDZ*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_S_Pz_Pz_b = I_ERI_F2xz_S_D2z_S_b+CDZ*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_S_Pz_Pz_b = I_ERI_Fx2y_S_D2z_S_b+CDZ*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_S_Pz_Pz_b = I_ERI_Fxyz_S_D2z_S_b+CDZ*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_S_Pz_Pz_b = I_ERI_Fx2z_S_D2z_S_b+CDZ*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_S_Pz_Pz_b = I_ERI_F3y_S_D2z_S_b+CDZ*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_S_Pz_Pz_b = I_ERI_F2yz_S_D2z_S_b+CDZ*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_S_Pz_Pz_b = I_ERI_Fy2z_S_D2z_S_b+CDZ*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_S_Pz_Pz_b = I_ERI_F3z_S_D2z_S_b+CDZ*I_ERI_F3z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_S_D_P_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_S_F_S_c
   * RHS shell quartet name: SQ_ERI_P_S_D_S_c
   ************************************************************/
  Double I_ERI_Px_S_D2x_Px_c = I_ERI_Px_S_F3x_S_c+CDX*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_S_D2x_Px_c = I_ERI_Py_S_F3x_S_c+CDX*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_S_D2x_Px_c = I_ERI_Pz_S_F3x_S_c+CDX*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_S_Dxy_Px_c = I_ERI_Px_S_F2xy_S_c+CDX*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_S_Dxy_Px_c = I_ERI_Py_S_F2xy_S_c+CDX*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_S_Dxy_Px_c = I_ERI_Pz_S_F2xy_S_c+CDX*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_S_Dxz_Px_c = I_ERI_Px_S_F2xz_S_c+CDX*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_S_Dxz_Px_c = I_ERI_Py_S_F2xz_S_c+CDX*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_S_Dxz_Px_c = I_ERI_Pz_S_F2xz_S_c+CDX*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_S_D2y_Px_c = I_ERI_Px_S_Fx2y_S_c+CDX*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_S_D2y_Px_c = I_ERI_Py_S_Fx2y_S_c+CDX*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_S_D2y_Px_c = I_ERI_Pz_S_Fx2y_S_c+CDX*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_S_Dyz_Px_c = I_ERI_Px_S_Fxyz_S_c+CDX*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_S_Dyz_Px_c = I_ERI_Py_S_Fxyz_S_c+CDX*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_S_Dyz_Px_c = I_ERI_Pz_S_Fxyz_S_c+CDX*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_S_D2z_Px_c = I_ERI_Px_S_Fx2z_S_c+CDX*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_S_D2z_Px_c = I_ERI_Py_S_Fx2z_S_c+CDX*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_S_D2z_Px_c = I_ERI_Pz_S_Fx2z_S_c+CDX*I_ERI_Pz_S_D2z_S_c;
  Double I_ERI_Px_S_D2x_Py_c = I_ERI_Px_S_F2xy_S_c+CDY*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_S_D2x_Py_c = I_ERI_Py_S_F2xy_S_c+CDY*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_S_D2x_Py_c = I_ERI_Pz_S_F2xy_S_c+CDY*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_S_Dxy_Py_c = I_ERI_Px_S_Fx2y_S_c+CDY*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_S_Dxy_Py_c = I_ERI_Py_S_Fx2y_S_c+CDY*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_S_Dxy_Py_c = I_ERI_Pz_S_Fx2y_S_c+CDY*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_S_Dxz_Py_c = I_ERI_Px_S_Fxyz_S_c+CDY*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_S_Dxz_Py_c = I_ERI_Py_S_Fxyz_S_c+CDY*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_S_Dxz_Py_c = I_ERI_Pz_S_Fxyz_S_c+CDY*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_S_D2y_Py_c = I_ERI_Px_S_F3y_S_c+CDY*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_S_D2y_Py_c = I_ERI_Py_S_F3y_S_c+CDY*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_S_D2y_Py_c = I_ERI_Pz_S_F3y_S_c+CDY*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_S_Dyz_Py_c = I_ERI_Px_S_F2yz_S_c+CDY*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_S_Dyz_Py_c = I_ERI_Py_S_F2yz_S_c+CDY*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_S_Dyz_Py_c = I_ERI_Pz_S_F2yz_S_c+CDY*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_S_D2z_Py_c = I_ERI_Px_S_Fy2z_S_c+CDY*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_S_D2z_Py_c = I_ERI_Py_S_Fy2z_S_c+CDY*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_S_D2z_Py_c = I_ERI_Pz_S_Fy2z_S_c+CDY*I_ERI_Pz_S_D2z_S_c;
  Double I_ERI_Px_S_D2x_Pz_c = I_ERI_Px_S_F2xz_S_c+CDZ*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_S_D2x_Pz_c = I_ERI_Py_S_F2xz_S_c+CDZ*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_S_D2x_Pz_c = I_ERI_Pz_S_F2xz_S_c+CDZ*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_S_Dxy_Pz_c = I_ERI_Px_S_Fxyz_S_c+CDZ*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_S_Dxy_Pz_c = I_ERI_Py_S_Fxyz_S_c+CDZ*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_S_Dxy_Pz_c = I_ERI_Pz_S_Fxyz_S_c+CDZ*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_S_Dxz_Pz_c = I_ERI_Px_S_Fx2z_S_c+CDZ*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_S_Dxz_Pz_c = I_ERI_Py_S_Fx2z_S_c+CDZ*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_S_Dxz_Pz_c = I_ERI_Pz_S_Fx2z_S_c+CDZ*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_S_D2y_Pz_c = I_ERI_Px_S_F2yz_S_c+CDZ*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_S_D2y_Pz_c = I_ERI_Py_S_F2yz_S_c+CDZ*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_S_D2y_Pz_c = I_ERI_Pz_S_F2yz_S_c+CDZ*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_S_Dyz_Pz_c = I_ERI_Px_S_Fy2z_S_c+CDZ*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_S_Dyz_Pz_c = I_ERI_Py_S_Fy2z_S_c+CDZ*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_S_Dyz_Pz_c = I_ERI_Pz_S_Fy2z_S_c+CDZ*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_S_D2z_Pz_c = I_ERI_Px_S_F3z_S_c+CDZ*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_S_D2z_Pz_c = I_ERI_Py_S_F3z_S_c+CDZ*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_S_D2z_Pz_c = I_ERI_Pz_S_F3z_S_c+CDZ*I_ERI_Pz_S_D2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_D_P_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_c
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   ************************************************************/
  Double I_ERI_D2x_S_D2x_Px_c = I_ERI_D2x_S_F3x_S_c+CDX*I_ERI_D2x_S_D2x_S_c;
  Double I_ERI_Dxy_S_D2x_Px_c = I_ERI_Dxy_S_F3x_S_c+CDX*I_ERI_Dxy_S_D2x_S_c;
  Double I_ERI_Dxz_S_D2x_Px_c = I_ERI_Dxz_S_F3x_S_c+CDX*I_ERI_Dxz_S_D2x_S_c;
  Double I_ERI_D2y_S_D2x_Px_c = I_ERI_D2y_S_F3x_S_c+CDX*I_ERI_D2y_S_D2x_S_c;
  Double I_ERI_Dyz_S_D2x_Px_c = I_ERI_Dyz_S_F3x_S_c+CDX*I_ERI_Dyz_S_D2x_S_c;
  Double I_ERI_D2z_S_D2x_Px_c = I_ERI_D2z_S_F3x_S_c+CDX*I_ERI_D2z_S_D2x_S_c;
  Double I_ERI_D2x_S_Dxy_Px_c = I_ERI_D2x_S_F2xy_S_c+CDX*I_ERI_D2x_S_Dxy_S_c;
  Double I_ERI_Dxy_S_Dxy_Px_c = I_ERI_Dxy_S_F2xy_S_c+CDX*I_ERI_Dxy_S_Dxy_S_c;
  Double I_ERI_Dxz_S_Dxy_Px_c = I_ERI_Dxz_S_F2xy_S_c+CDX*I_ERI_Dxz_S_Dxy_S_c;
  Double I_ERI_D2y_S_Dxy_Px_c = I_ERI_D2y_S_F2xy_S_c+CDX*I_ERI_D2y_S_Dxy_S_c;
  Double I_ERI_Dyz_S_Dxy_Px_c = I_ERI_Dyz_S_F2xy_S_c+CDX*I_ERI_Dyz_S_Dxy_S_c;
  Double I_ERI_D2z_S_Dxy_Px_c = I_ERI_D2z_S_F2xy_S_c+CDX*I_ERI_D2z_S_Dxy_S_c;
  Double I_ERI_D2x_S_Dxz_Px_c = I_ERI_D2x_S_F2xz_S_c+CDX*I_ERI_D2x_S_Dxz_S_c;
  Double I_ERI_Dxy_S_Dxz_Px_c = I_ERI_Dxy_S_F2xz_S_c+CDX*I_ERI_Dxy_S_Dxz_S_c;
  Double I_ERI_Dxz_S_Dxz_Px_c = I_ERI_Dxz_S_F2xz_S_c+CDX*I_ERI_Dxz_S_Dxz_S_c;
  Double I_ERI_D2y_S_Dxz_Px_c = I_ERI_D2y_S_F2xz_S_c+CDX*I_ERI_D2y_S_Dxz_S_c;
  Double I_ERI_Dyz_S_Dxz_Px_c = I_ERI_Dyz_S_F2xz_S_c+CDX*I_ERI_Dyz_S_Dxz_S_c;
  Double I_ERI_D2z_S_Dxz_Px_c = I_ERI_D2z_S_F2xz_S_c+CDX*I_ERI_D2z_S_Dxz_S_c;
  Double I_ERI_D2x_S_D2y_Px_c = I_ERI_D2x_S_Fx2y_S_c+CDX*I_ERI_D2x_S_D2y_S_c;
  Double I_ERI_Dxy_S_D2y_Px_c = I_ERI_Dxy_S_Fx2y_S_c+CDX*I_ERI_Dxy_S_D2y_S_c;
  Double I_ERI_Dxz_S_D2y_Px_c = I_ERI_Dxz_S_Fx2y_S_c+CDX*I_ERI_Dxz_S_D2y_S_c;
  Double I_ERI_D2y_S_D2y_Px_c = I_ERI_D2y_S_Fx2y_S_c+CDX*I_ERI_D2y_S_D2y_S_c;
  Double I_ERI_Dyz_S_D2y_Px_c = I_ERI_Dyz_S_Fx2y_S_c+CDX*I_ERI_Dyz_S_D2y_S_c;
  Double I_ERI_D2z_S_D2y_Px_c = I_ERI_D2z_S_Fx2y_S_c+CDX*I_ERI_D2z_S_D2y_S_c;
  Double I_ERI_D2x_S_Dyz_Px_c = I_ERI_D2x_S_Fxyz_S_c+CDX*I_ERI_D2x_S_Dyz_S_c;
  Double I_ERI_Dxy_S_Dyz_Px_c = I_ERI_Dxy_S_Fxyz_S_c+CDX*I_ERI_Dxy_S_Dyz_S_c;
  Double I_ERI_Dxz_S_Dyz_Px_c = I_ERI_Dxz_S_Fxyz_S_c+CDX*I_ERI_Dxz_S_Dyz_S_c;
  Double I_ERI_D2y_S_Dyz_Px_c = I_ERI_D2y_S_Fxyz_S_c+CDX*I_ERI_D2y_S_Dyz_S_c;
  Double I_ERI_Dyz_S_Dyz_Px_c = I_ERI_Dyz_S_Fxyz_S_c+CDX*I_ERI_Dyz_S_Dyz_S_c;
  Double I_ERI_D2z_S_Dyz_Px_c = I_ERI_D2z_S_Fxyz_S_c+CDX*I_ERI_D2z_S_Dyz_S_c;
  Double I_ERI_D2x_S_D2z_Px_c = I_ERI_D2x_S_Fx2z_S_c+CDX*I_ERI_D2x_S_D2z_S_c;
  Double I_ERI_Dxy_S_D2z_Px_c = I_ERI_Dxy_S_Fx2z_S_c+CDX*I_ERI_Dxy_S_D2z_S_c;
  Double I_ERI_Dxz_S_D2z_Px_c = I_ERI_Dxz_S_Fx2z_S_c+CDX*I_ERI_Dxz_S_D2z_S_c;
  Double I_ERI_D2y_S_D2z_Px_c = I_ERI_D2y_S_Fx2z_S_c+CDX*I_ERI_D2y_S_D2z_S_c;
  Double I_ERI_Dyz_S_D2z_Px_c = I_ERI_Dyz_S_Fx2z_S_c+CDX*I_ERI_Dyz_S_D2z_S_c;
  Double I_ERI_D2z_S_D2z_Px_c = I_ERI_D2z_S_Fx2z_S_c+CDX*I_ERI_D2z_S_D2z_S_c;
  Double I_ERI_D2x_S_D2x_Py_c = I_ERI_D2x_S_F2xy_S_c+CDY*I_ERI_D2x_S_D2x_S_c;
  Double I_ERI_Dxy_S_D2x_Py_c = I_ERI_Dxy_S_F2xy_S_c+CDY*I_ERI_Dxy_S_D2x_S_c;
  Double I_ERI_Dxz_S_D2x_Py_c = I_ERI_Dxz_S_F2xy_S_c+CDY*I_ERI_Dxz_S_D2x_S_c;
  Double I_ERI_D2y_S_D2x_Py_c = I_ERI_D2y_S_F2xy_S_c+CDY*I_ERI_D2y_S_D2x_S_c;
  Double I_ERI_Dyz_S_D2x_Py_c = I_ERI_Dyz_S_F2xy_S_c+CDY*I_ERI_Dyz_S_D2x_S_c;
  Double I_ERI_D2z_S_D2x_Py_c = I_ERI_D2z_S_F2xy_S_c+CDY*I_ERI_D2z_S_D2x_S_c;
  Double I_ERI_D2x_S_Dxy_Py_c = I_ERI_D2x_S_Fx2y_S_c+CDY*I_ERI_D2x_S_Dxy_S_c;
  Double I_ERI_Dxy_S_Dxy_Py_c = I_ERI_Dxy_S_Fx2y_S_c+CDY*I_ERI_Dxy_S_Dxy_S_c;
  Double I_ERI_Dxz_S_Dxy_Py_c = I_ERI_Dxz_S_Fx2y_S_c+CDY*I_ERI_Dxz_S_Dxy_S_c;
  Double I_ERI_D2y_S_Dxy_Py_c = I_ERI_D2y_S_Fx2y_S_c+CDY*I_ERI_D2y_S_Dxy_S_c;
  Double I_ERI_Dyz_S_Dxy_Py_c = I_ERI_Dyz_S_Fx2y_S_c+CDY*I_ERI_Dyz_S_Dxy_S_c;
  Double I_ERI_D2z_S_Dxy_Py_c = I_ERI_D2z_S_Fx2y_S_c+CDY*I_ERI_D2z_S_Dxy_S_c;
  Double I_ERI_D2x_S_Dxz_Py_c = I_ERI_D2x_S_Fxyz_S_c+CDY*I_ERI_D2x_S_Dxz_S_c;
  Double I_ERI_Dxy_S_Dxz_Py_c = I_ERI_Dxy_S_Fxyz_S_c+CDY*I_ERI_Dxy_S_Dxz_S_c;
  Double I_ERI_Dxz_S_Dxz_Py_c = I_ERI_Dxz_S_Fxyz_S_c+CDY*I_ERI_Dxz_S_Dxz_S_c;
  Double I_ERI_D2y_S_Dxz_Py_c = I_ERI_D2y_S_Fxyz_S_c+CDY*I_ERI_D2y_S_Dxz_S_c;
  Double I_ERI_Dyz_S_Dxz_Py_c = I_ERI_Dyz_S_Fxyz_S_c+CDY*I_ERI_Dyz_S_Dxz_S_c;
  Double I_ERI_D2z_S_Dxz_Py_c = I_ERI_D2z_S_Fxyz_S_c+CDY*I_ERI_D2z_S_Dxz_S_c;
  Double I_ERI_D2x_S_D2y_Py_c = I_ERI_D2x_S_F3y_S_c+CDY*I_ERI_D2x_S_D2y_S_c;
  Double I_ERI_Dxy_S_D2y_Py_c = I_ERI_Dxy_S_F3y_S_c+CDY*I_ERI_Dxy_S_D2y_S_c;
  Double I_ERI_Dxz_S_D2y_Py_c = I_ERI_Dxz_S_F3y_S_c+CDY*I_ERI_Dxz_S_D2y_S_c;
  Double I_ERI_D2y_S_D2y_Py_c = I_ERI_D2y_S_F3y_S_c+CDY*I_ERI_D2y_S_D2y_S_c;
  Double I_ERI_Dyz_S_D2y_Py_c = I_ERI_Dyz_S_F3y_S_c+CDY*I_ERI_Dyz_S_D2y_S_c;
  Double I_ERI_D2z_S_D2y_Py_c = I_ERI_D2z_S_F3y_S_c+CDY*I_ERI_D2z_S_D2y_S_c;
  Double I_ERI_D2x_S_Dyz_Py_c = I_ERI_D2x_S_F2yz_S_c+CDY*I_ERI_D2x_S_Dyz_S_c;
  Double I_ERI_Dxy_S_Dyz_Py_c = I_ERI_Dxy_S_F2yz_S_c+CDY*I_ERI_Dxy_S_Dyz_S_c;
  Double I_ERI_Dxz_S_Dyz_Py_c = I_ERI_Dxz_S_F2yz_S_c+CDY*I_ERI_Dxz_S_Dyz_S_c;
  Double I_ERI_D2y_S_Dyz_Py_c = I_ERI_D2y_S_F2yz_S_c+CDY*I_ERI_D2y_S_Dyz_S_c;
  Double I_ERI_Dyz_S_Dyz_Py_c = I_ERI_Dyz_S_F2yz_S_c+CDY*I_ERI_Dyz_S_Dyz_S_c;
  Double I_ERI_D2z_S_Dyz_Py_c = I_ERI_D2z_S_F2yz_S_c+CDY*I_ERI_D2z_S_Dyz_S_c;
  Double I_ERI_D2x_S_D2z_Py_c = I_ERI_D2x_S_Fy2z_S_c+CDY*I_ERI_D2x_S_D2z_S_c;
  Double I_ERI_Dxy_S_D2z_Py_c = I_ERI_Dxy_S_Fy2z_S_c+CDY*I_ERI_Dxy_S_D2z_S_c;
  Double I_ERI_Dxz_S_D2z_Py_c = I_ERI_Dxz_S_Fy2z_S_c+CDY*I_ERI_Dxz_S_D2z_S_c;
  Double I_ERI_D2y_S_D2z_Py_c = I_ERI_D2y_S_Fy2z_S_c+CDY*I_ERI_D2y_S_D2z_S_c;
  Double I_ERI_Dyz_S_D2z_Py_c = I_ERI_Dyz_S_Fy2z_S_c+CDY*I_ERI_Dyz_S_D2z_S_c;
  Double I_ERI_D2z_S_D2z_Py_c = I_ERI_D2z_S_Fy2z_S_c+CDY*I_ERI_D2z_S_D2z_S_c;
  Double I_ERI_D2x_S_D2x_Pz_c = I_ERI_D2x_S_F2xz_S_c+CDZ*I_ERI_D2x_S_D2x_S_c;
  Double I_ERI_Dxy_S_D2x_Pz_c = I_ERI_Dxy_S_F2xz_S_c+CDZ*I_ERI_Dxy_S_D2x_S_c;
  Double I_ERI_Dxz_S_D2x_Pz_c = I_ERI_Dxz_S_F2xz_S_c+CDZ*I_ERI_Dxz_S_D2x_S_c;
  Double I_ERI_D2y_S_D2x_Pz_c = I_ERI_D2y_S_F2xz_S_c+CDZ*I_ERI_D2y_S_D2x_S_c;
  Double I_ERI_Dyz_S_D2x_Pz_c = I_ERI_Dyz_S_F2xz_S_c+CDZ*I_ERI_Dyz_S_D2x_S_c;
  Double I_ERI_D2z_S_D2x_Pz_c = I_ERI_D2z_S_F2xz_S_c+CDZ*I_ERI_D2z_S_D2x_S_c;
  Double I_ERI_D2x_S_Dxy_Pz_c = I_ERI_D2x_S_Fxyz_S_c+CDZ*I_ERI_D2x_S_Dxy_S_c;
  Double I_ERI_Dxy_S_Dxy_Pz_c = I_ERI_Dxy_S_Fxyz_S_c+CDZ*I_ERI_Dxy_S_Dxy_S_c;
  Double I_ERI_Dxz_S_Dxy_Pz_c = I_ERI_Dxz_S_Fxyz_S_c+CDZ*I_ERI_Dxz_S_Dxy_S_c;
  Double I_ERI_D2y_S_Dxy_Pz_c = I_ERI_D2y_S_Fxyz_S_c+CDZ*I_ERI_D2y_S_Dxy_S_c;
  Double I_ERI_Dyz_S_Dxy_Pz_c = I_ERI_Dyz_S_Fxyz_S_c+CDZ*I_ERI_Dyz_S_Dxy_S_c;
  Double I_ERI_D2z_S_Dxy_Pz_c = I_ERI_D2z_S_Fxyz_S_c+CDZ*I_ERI_D2z_S_Dxy_S_c;
  Double I_ERI_D2x_S_Dxz_Pz_c = I_ERI_D2x_S_Fx2z_S_c+CDZ*I_ERI_D2x_S_Dxz_S_c;
  Double I_ERI_Dxy_S_Dxz_Pz_c = I_ERI_Dxy_S_Fx2z_S_c+CDZ*I_ERI_Dxy_S_Dxz_S_c;
  Double I_ERI_Dxz_S_Dxz_Pz_c = I_ERI_Dxz_S_Fx2z_S_c+CDZ*I_ERI_Dxz_S_Dxz_S_c;
  Double I_ERI_D2y_S_Dxz_Pz_c = I_ERI_D2y_S_Fx2z_S_c+CDZ*I_ERI_D2y_S_Dxz_S_c;
  Double I_ERI_Dyz_S_Dxz_Pz_c = I_ERI_Dyz_S_Fx2z_S_c+CDZ*I_ERI_Dyz_S_Dxz_S_c;
  Double I_ERI_D2z_S_Dxz_Pz_c = I_ERI_D2z_S_Fx2z_S_c+CDZ*I_ERI_D2z_S_Dxz_S_c;
  Double I_ERI_D2x_S_D2y_Pz_c = I_ERI_D2x_S_F2yz_S_c+CDZ*I_ERI_D2x_S_D2y_S_c;
  Double I_ERI_Dxy_S_D2y_Pz_c = I_ERI_Dxy_S_F2yz_S_c+CDZ*I_ERI_Dxy_S_D2y_S_c;
  Double I_ERI_Dxz_S_D2y_Pz_c = I_ERI_Dxz_S_F2yz_S_c+CDZ*I_ERI_Dxz_S_D2y_S_c;
  Double I_ERI_D2y_S_D2y_Pz_c = I_ERI_D2y_S_F2yz_S_c+CDZ*I_ERI_D2y_S_D2y_S_c;
  Double I_ERI_Dyz_S_D2y_Pz_c = I_ERI_Dyz_S_F2yz_S_c+CDZ*I_ERI_Dyz_S_D2y_S_c;
  Double I_ERI_D2z_S_D2y_Pz_c = I_ERI_D2z_S_F2yz_S_c+CDZ*I_ERI_D2z_S_D2y_S_c;
  Double I_ERI_D2x_S_Dyz_Pz_c = I_ERI_D2x_S_Fy2z_S_c+CDZ*I_ERI_D2x_S_Dyz_S_c;
  Double I_ERI_Dxy_S_Dyz_Pz_c = I_ERI_Dxy_S_Fy2z_S_c+CDZ*I_ERI_Dxy_S_Dyz_S_c;
  Double I_ERI_Dxz_S_Dyz_Pz_c = I_ERI_Dxz_S_Fy2z_S_c+CDZ*I_ERI_Dxz_S_Dyz_S_c;
  Double I_ERI_D2y_S_Dyz_Pz_c = I_ERI_D2y_S_Fy2z_S_c+CDZ*I_ERI_D2y_S_Dyz_S_c;
  Double I_ERI_Dyz_S_Dyz_Pz_c = I_ERI_Dyz_S_Fy2z_S_c+CDZ*I_ERI_Dyz_S_Dyz_S_c;
  Double I_ERI_D2z_S_Dyz_Pz_c = I_ERI_D2z_S_Fy2z_S_c+CDZ*I_ERI_D2z_S_Dyz_S_c;
  Double I_ERI_D2x_S_D2z_Pz_c = I_ERI_D2x_S_F3z_S_c+CDZ*I_ERI_D2x_S_D2z_S_c;
  Double I_ERI_Dxy_S_D2z_Pz_c = I_ERI_Dxy_S_F3z_S_c+CDZ*I_ERI_Dxy_S_D2z_S_c;
  Double I_ERI_Dxz_S_D2z_Pz_c = I_ERI_Dxz_S_F3z_S_c+CDZ*I_ERI_Dxz_S_D2z_S_c;
  Double I_ERI_D2y_S_D2z_Pz_c = I_ERI_D2y_S_F3z_S_c+CDZ*I_ERI_D2y_S_D2z_S_c;
  Double I_ERI_Dyz_S_D2z_Pz_c = I_ERI_Dyz_S_F3z_S_c+CDZ*I_ERI_Dyz_S_D2z_S_c;
  Double I_ERI_D2z_S_D2z_Pz_c = I_ERI_D2z_S_F3z_S_c+CDZ*I_ERI_D2z_S_D2z_S_c;

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
   * shell quartet name: SQ_ERI_P_P_S_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_P
   * RHS shell quartet name: SQ_ERI_P_S_S_P
   ************************************************************/
  Double I_ERI_Px_Px_S_Px = I_ERI_D2x_S_S_Px+ABX*I_ERI_Px_S_S_Px;
  Double I_ERI_Py_Px_S_Px = I_ERI_Dxy_S_S_Px+ABX*I_ERI_Py_S_S_Px;
  Double I_ERI_Pz_Px_S_Px = I_ERI_Dxz_S_S_Px+ABX*I_ERI_Pz_S_S_Px;
  Double I_ERI_Px_Py_S_Px = I_ERI_Dxy_S_S_Px+ABY*I_ERI_Px_S_S_Px;
  Double I_ERI_Py_Py_S_Px = I_ERI_D2y_S_S_Px+ABY*I_ERI_Py_S_S_Px;
  Double I_ERI_Pz_Py_S_Px = I_ERI_Dyz_S_S_Px+ABY*I_ERI_Pz_S_S_Px;
  Double I_ERI_Px_Pz_S_Px = I_ERI_Dxz_S_S_Px+ABZ*I_ERI_Px_S_S_Px;
  Double I_ERI_Py_Pz_S_Px = I_ERI_Dyz_S_S_Px+ABZ*I_ERI_Py_S_S_Px;
  Double I_ERI_Pz_Pz_S_Px = I_ERI_D2z_S_S_Px+ABZ*I_ERI_Pz_S_S_Px;
  Double I_ERI_Px_Px_S_Py = I_ERI_D2x_S_S_Py+ABX*I_ERI_Px_S_S_Py;
  Double I_ERI_Py_Px_S_Py = I_ERI_Dxy_S_S_Py+ABX*I_ERI_Py_S_S_Py;
  Double I_ERI_Pz_Px_S_Py = I_ERI_Dxz_S_S_Py+ABX*I_ERI_Pz_S_S_Py;
  Double I_ERI_Px_Py_S_Py = I_ERI_Dxy_S_S_Py+ABY*I_ERI_Px_S_S_Py;
  Double I_ERI_Py_Py_S_Py = I_ERI_D2y_S_S_Py+ABY*I_ERI_Py_S_S_Py;
  Double I_ERI_Pz_Py_S_Py = I_ERI_Dyz_S_S_Py+ABY*I_ERI_Pz_S_S_Py;
  Double I_ERI_Px_Pz_S_Py = I_ERI_Dxz_S_S_Py+ABZ*I_ERI_Px_S_S_Py;
  Double I_ERI_Py_Pz_S_Py = I_ERI_Dyz_S_S_Py+ABZ*I_ERI_Py_S_S_Py;
  Double I_ERI_Pz_Pz_S_Py = I_ERI_D2z_S_S_Py+ABZ*I_ERI_Pz_S_S_Py;
  Double I_ERI_Px_Px_S_Pz = I_ERI_D2x_S_S_Pz+ABX*I_ERI_Px_S_S_Pz;
  Double I_ERI_Py_Px_S_Pz = I_ERI_Dxy_S_S_Pz+ABX*I_ERI_Py_S_S_Pz;
  Double I_ERI_Pz_Px_S_Pz = I_ERI_Dxz_S_S_Pz+ABX*I_ERI_Pz_S_S_Pz;
  Double I_ERI_Px_Py_S_Pz = I_ERI_Dxy_S_S_Pz+ABY*I_ERI_Px_S_S_Pz;
  Double I_ERI_Py_Py_S_Pz = I_ERI_D2y_S_S_Pz+ABY*I_ERI_Py_S_S_Pz;
  Double I_ERI_Pz_Py_S_Pz = I_ERI_Dyz_S_S_Pz+ABY*I_ERI_Pz_S_S_Pz;
  Double I_ERI_Px_Pz_S_Pz = I_ERI_Dxz_S_S_Pz+ABZ*I_ERI_Px_S_S_Pz;
  Double I_ERI_Py_Pz_S_Pz = I_ERI_Dyz_S_S_Pz+ABZ*I_ERI_Py_S_S_Pz;
  Double I_ERI_Pz_Pz_S_Pz = I_ERI_D2z_S_S_Pz+ABZ*I_ERI_Pz_S_S_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_a
   * RHS shell quartet name: SQ_ERI_D_S_P_P_a
   ************************************************************/
  Double I_ERI_D2x_Px_Px_Px_a = I_ERI_F3x_S_Px_Px_a+ABX*I_ERI_D2x_S_Px_Px_a;
  Double I_ERI_Dxy_Px_Px_Px_a = I_ERI_F2xy_S_Px_Px_a+ABX*I_ERI_Dxy_S_Px_Px_a;
  Double I_ERI_Dxz_Px_Px_Px_a = I_ERI_F2xz_S_Px_Px_a+ABX*I_ERI_Dxz_S_Px_Px_a;
  Double I_ERI_D2y_Px_Px_Px_a = I_ERI_Fx2y_S_Px_Px_a+ABX*I_ERI_D2y_S_Px_Px_a;
  Double I_ERI_Dyz_Px_Px_Px_a = I_ERI_Fxyz_S_Px_Px_a+ABX*I_ERI_Dyz_S_Px_Px_a;
  Double I_ERI_D2z_Px_Px_Px_a = I_ERI_Fx2z_S_Px_Px_a+ABX*I_ERI_D2z_S_Px_Px_a;
  Double I_ERI_D2x_Py_Px_Px_a = I_ERI_F2xy_S_Px_Px_a+ABY*I_ERI_D2x_S_Px_Px_a;
  Double I_ERI_Dxy_Py_Px_Px_a = I_ERI_Fx2y_S_Px_Px_a+ABY*I_ERI_Dxy_S_Px_Px_a;
  Double I_ERI_Dxz_Py_Px_Px_a = I_ERI_Fxyz_S_Px_Px_a+ABY*I_ERI_Dxz_S_Px_Px_a;
  Double I_ERI_D2y_Py_Px_Px_a = I_ERI_F3y_S_Px_Px_a+ABY*I_ERI_D2y_S_Px_Px_a;
  Double I_ERI_Dyz_Py_Px_Px_a = I_ERI_F2yz_S_Px_Px_a+ABY*I_ERI_Dyz_S_Px_Px_a;
  Double I_ERI_D2z_Py_Px_Px_a = I_ERI_Fy2z_S_Px_Px_a+ABY*I_ERI_D2z_S_Px_Px_a;
  Double I_ERI_D2x_Pz_Px_Px_a = I_ERI_F2xz_S_Px_Px_a+ABZ*I_ERI_D2x_S_Px_Px_a;
  Double I_ERI_Dxy_Pz_Px_Px_a = I_ERI_Fxyz_S_Px_Px_a+ABZ*I_ERI_Dxy_S_Px_Px_a;
  Double I_ERI_Dxz_Pz_Px_Px_a = I_ERI_Fx2z_S_Px_Px_a+ABZ*I_ERI_Dxz_S_Px_Px_a;
  Double I_ERI_D2y_Pz_Px_Px_a = I_ERI_F2yz_S_Px_Px_a+ABZ*I_ERI_D2y_S_Px_Px_a;
  Double I_ERI_Dyz_Pz_Px_Px_a = I_ERI_Fy2z_S_Px_Px_a+ABZ*I_ERI_Dyz_S_Px_Px_a;
  Double I_ERI_D2z_Pz_Px_Px_a = I_ERI_F3z_S_Px_Px_a+ABZ*I_ERI_D2z_S_Px_Px_a;
  Double I_ERI_D2x_Px_Py_Px_a = I_ERI_F3x_S_Py_Px_a+ABX*I_ERI_D2x_S_Py_Px_a;
  Double I_ERI_Dxy_Px_Py_Px_a = I_ERI_F2xy_S_Py_Px_a+ABX*I_ERI_Dxy_S_Py_Px_a;
  Double I_ERI_Dxz_Px_Py_Px_a = I_ERI_F2xz_S_Py_Px_a+ABX*I_ERI_Dxz_S_Py_Px_a;
  Double I_ERI_D2y_Px_Py_Px_a = I_ERI_Fx2y_S_Py_Px_a+ABX*I_ERI_D2y_S_Py_Px_a;
  Double I_ERI_Dyz_Px_Py_Px_a = I_ERI_Fxyz_S_Py_Px_a+ABX*I_ERI_Dyz_S_Py_Px_a;
  Double I_ERI_D2z_Px_Py_Px_a = I_ERI_Fx2z_S_Py_Px_a+ABX*I_ERI_D2z_S_Py_Px_a;
  Double I_ERI_D2x_Py_Py_Px_a = I_ERI_F2xy_S_Py_Px_a+ABY*I_ERI_D2x_S_Py_Px_a;
  Double I_ERI_Dxy_Py_Py_Px_a = I_ERI_Fx2y_S_Py_Px_a+ABY*I_ERI_Dxy_S_Py_Px_a;
  Double I_ERI_Dxz_Py_Py_Px_a = I_ERI_Fxyz_S_Py_Px_a+ABY*I_ERI_Dxz_S_Py_Px_a;
  Double I_ERI_D2y_Py_Py_Px_a = I_ERI_F3y_S_Py_Px_a+ABY*I_ERI_D2y_S_Py_Px_a;
  Double I_ERI_Dyz_Py_Py_Px_a = I_ERI_F2yz_S_Py_Px_a+ABY*I_ERI_Dyz_S_Py_Px_a;
  Double I_ERI_D2z_Py_Py_Px_a = I_ERI_Fy2z_S_Py_Px_a+ABY*I_ERI_D2z_S_Py_Px_a;
  Double I_ERI_D2x_Pz_Py_Px_a = I_ERI_F2xz_S_Py_Px_a+ABZ*I_ERI_D2x_S_Py_Px_a;
  Double I_ERI_Dxy_Pz_Py_Px_a = I_ERI_Fxyz_S_Py_Px_a+ABZ*I_ERI_Dxy_S_Py_Px_a;
  Double I_ERI_Dxz_Pz_Py_Px_a = I_ERI_Fx2z_S_Py_Px_a+ABZ*I_ERI_Dxz_S_Py_Px_a;
  Double I_ERI_D2y_Pz_Py_Px_a = I_ERI_F2yz_S_Py_Px_a+ABZ*I_ERI_D2y_S_Py_Px_a;
  Double I_ERI_Dyz_Pz_Py_Px_a = I_ERI_Fy2z_S_Py_Px_a+ABZ*I_ERI_Dyz_S_Py_Px_a;
  Double I_ERI_D2z_Pz_Py_Px_a = I_ERI_F3z_S_Py_Px_a+ABZ*I_ERI_D2z_S_Py_Px_a;
  Double I_ERI_D2x_Px_Pz_Px_a = I_ERI_F3x_S_Pz_Px_a+ABX*I_ERI_D2x_S_Pz_Px_a;
  Double I_ERI_Dxy_Px_Pz_Px_a = I_ERI_F2xy_S_Pz_Px_a+ABX*I_ERI_Dxy_S_Pz_Px_a;
  Double I_ERI_Dxz_Px_Pz_Px_a = I_ERI_F2xz_S_Pz_Px_a+ABX*I_ERI_Dxz_S_Pz_Px_a;
  Double I_ERI_D2y_Px_Pz_Px_a = I_ERI_Fx2y_S_Pz_Px_a+ABX*I_ERI_D2y_S_Pz_Px_a;
  Double I_ERI_Dyz_Px_Pz_Px_a = I_ERI_Fxyz_S_Pz_Px_a+ABX*I_ERI_Dyz_S_Pz_Px_a;
  Double I_ERI_D2z_Px_Pz_Px_a = I_ERI_Fx2z_S_Pz_Px_a+ABX*I_ERI_D2z_S_Pz_Px_a;
  Double I_ERI_D2x_Py_Pz_Px_a = I_ERI_F2xy_S_Pz_Px_a+ABY*I_ERI_D2x_S_Pz_Px_a;
  Double I_ERI_Dxy_Py_Pz_Px_a = I_ERI_Fx2y_S_Pz_Px_a+ABY*I_ERI_Dxy_S_Pz_Px_a;
  Double I_ERI_Dxz_Py_Pz_Px_a = I_ERI_Fxyz_S_Pz_Px_a+ABY*I_ERI_Dxz_S_Pz_Px_a;
  Double I_ERI_D2y_Py_Pz_Px_a = I_ERI_F3y_S_Pz_Px_a+ABY*I_ERI_D2y_S_Pz_Px_a;
  Double I_ERI_Dyz_Py_Pz_Px_a = I_ERI_F2yz_S_Pz_Px_a+ABY*I_ERI_Dyz_S_Pz_Px_a;
  Double I_ERI_D2z_Py_Pz_Px_a = I_ERI_Fy2z_S_Pz_Px_a+ABY*I_ERI_D2z_S_Pz_Px_a;
  Double I_ERI_D2x_Pz_Pz_Px_a = I_ERI_F2xz_S_Pz_Px_a+ABZ*I_ERI_D2x_S_Pz_Px_a;
  Double I_ERI_Dxy_Pz_Pz_Px_a = I_ERI_Fxyz_S_Pz_Px_a+ABZ*I_ERI_Dxy_S_Pz_Px_a;
  Double I_ERI_Dxz_Pz_Pz_Px_a = I_ERI_Fx2z_S_Pz_Px_a+ABZ*I_ERI_Dxz_S_Pz_Px_a;
  Double I_ERI_D2y_Pz_Pz_Px_a = I_ERI_F2yz_S_Pz_Px_a+ABZ*I_ERI_D2y_S_Pz_Px_a;
  Double I_ERI_Dyz_Pz_Pz_Px_a = I_ERI_Fy2z_S_Pz_Px_a+ABZ*I_ERI_Dyz_S_Pz_Px_a;
  Double I_ERI_D2z_Pz_Pz_Px_a = I_ERI_F3z_S_Pz_Px_a+ABZ*I_ERI_D2z_S_Pz_Px_a;
  Double I_ERI_D2x_Px_Px_Py_a = I_ERI_F3x_S_Px_Py_a+ABX*I_ERI_D2x_S_Px_Py_a;
  Double I_ERI_Dxy_Px_Px_Py_a = I_ERI_F2xy_S_Px_Py_a+ABX*I_ERI_Dxy_S_Px_Py_a;
  Double I_ERI_Dxz_Px_Px_Py_a = I_ERI_F2xz_S_Px_Py_a+ABX*I_ERI_Dxz_S_Px_Py_a;
  Double I_ERI_D2y_Px_Px_Py_a = I_ERI_Fx2y_S_Px_Py_a+ABX*I_ERI_D2y_S_Px_Py_a;
  Double I_ERI_Dyz_Px_Px_Py_a = I_ERI_Fxyz_S_Px_Py_a+ABX*I_ERI_Dyz_S_Px_Py_a;
  Double I_ERI_D2z_Px_Px_Py_a = I_ERI_Fx2z_S_Px_Py_a+ABX*I_ERI_D2z_S_Px_Py_a;
  Double I_ERI_D2x_Py_Px_Py_a = I_ERI_F2xy_S_Px_Py_a+ABY*I_ERI_D2x_S_Px_Py_a;
  Double I_ERI_Dxy_Py_Px_Py_a = I_ERI_Fx2y_S_Px_Py_a+ABY*I_ERI_Dxy_S_Px_Py_a;
  Double I_ERI_Dxz_Py_Px_Py_a = I_ERI_Fxyz_S_Px_Py_a+ABY*I_ERI_Dxz_S_Px_Py_a;
  Double I_ERI_D2y_Py_Px_Py_a = I_ERI_F3y_S_Px_Py_a+ABY*I_ERI_D2y_S_Px_Py_a;
  Double I_ERI_Dyz_Py_Px_Py_a = I_ERI_F2yz_S_Px_Py_a+ABY*I_ERI_Dyz_S_Px_Py_a;
  Double I_ERI_D2z_Py_Px_Py_a = I_ERI_Fy2z_S_Px_Py_a+ABY*I_ERI_D2z_S_Px_Py_a;
  Double I_ERI_D2x_Pz_Px_Py_a = I_ERI_F2xz_S_Px_Py_a+ABZ*I_ERI_D2x_S_Px_Py_a;
  Double I_ERI_Dxy_Pz_Px_Py_a = I_ERI_Fxyz_S_Px_Py_a+ABZ*I_ERI_Dxy_S_Px_Py_a;
  Double I_ERI_Dxz_Pz_Px_Py_a = I_ERI_Fx2z_S_Px_Py_a+ABZ*I_ERI_Dxz_S_Px_Py_a;
  Double I_ERI_D2y_Pz_Px_Py_a = I_ERI_F2yz_S_Px_Py_a+ABZ*I_ERI_D2y_S_Px_Py_a;
  Double I_ERI_Dyz_Pz_Px_Py_a = I_ERI_Fy2z_S_Px_Py_a+ABZ*I_ERI_Dyz_S_Px_Py_a;
  Double I_ERI_D2z_Pz_Px_Py_a = I_ERI_F3z_S_Px_Py_a+ABZ*I_ERI_D2z_S_Px_Py_a;
  Double I_ERI_D2x_Px_Py_Py_a = I_ERI_F3x_S_Py_Py_a+ABX*I_ERI_D2x_S_Py_Py_a;
  Double I_ERI_Dxy_Px_Py_Py_a = I_ERI_F2xy_S_Py_Py_a+ABX*I_ERI_Dxy_S_Py_Py_a;
  Double I_ERI_Dxz_Px_Py_Py_a = I_ERI_F2xz_S_Py_Py_a+ABX*I_ERI_Dxz_S_Py_Py_a;
  Double I_ERI_D2y_Px_Py_Py_a = I_ERI_Fx2y_S_Py_Py_a+ABX*I_ERI_D2y_S_Py_Py_a;
  Double I_ERI_Dyz_Px_Py_Py_a = I_ERI_Fxyz_S_Py_Py_a+ABX*I_ERI_Dyz_S_Py_Py_a;
  Double I_ERI_D2z_Px_Py_Py_a = I_ERI_Fx2z_S_Py_Py_a+ABX*I_ERI_D2z_S_Py_Py_a;
  Double I_ERI_D2x_Py_Py_Py_a = I_ERI_F2xy_S_Py_Py_a+ABY*I_ERI_D2x_S_Py_Py_a;
  Double I_ERI_Dxy_Py_Py_Py_a = I_ERI_Fx2y_S_Py_Py_a+ABY*I_ERI_Dxy_S_Py_Py_a;
  Double I_ERI_Dxz_Py_Py_Py_a = I_ERI_Fxyz_S_Py_Py_a+ABY*I_ERI_Dxz_S_Py_Py_a;
  Double I_ERI_D2y_Py_Py_Py_a = I_ERI_F3y_S_Py_Py_a+ABY*I_ERI_D2y_S_Py_Py_a;
  Double I_ERI_Dyz_Py_Py_Py_a = I_ERI_F2yz_S_Py_Py_a+ABY*I_ERI_Dyz_S_Py_Py_a;
  Double I_ERI_D2z_Py_Py_Py_a = I_ERI_Fy2z_S_Py_Py_a+ABY*I_ERI_D2z_S_Py_Py_a;
  Double I_ERI_D2x_Pz_Py_Py_a = I_ERI_F2xz_S_Py_Py_a+ABZ*I_ERI_D2x_S_Py_Py_a;
  Double I_ERI_Dxy_Pz_Py_Py_a = I_ERI_Fxyz_S_Py_Py_a+ABZ*I_ERI_Dxy_S_Py_Py_a;
  Double I_ERI_Dxz_Pz_Py_Py_a = I_ERI_Fx2z_S_Py_Py_a+ABZ*I_ERI_Dxz_S_Py_Py_a;
  Double I_ERI_D2y_Pz_Py_Py_a = I_ERI_F2yz_S_Py_Py_a+ABZ*I_ERI_D2y_S_Py_Py_a;
  Double I_ERI_Dyz_Pz_Py_Py_a = I_ERI_Fy2z_S_Py_Py_a+ABZ*I_ERI_Dyz_S_Py_Py_a;
  Double I_ERI_D2z_Pz_Py_Py_a = I_ERI_F3z_S_Py_Py_a+ABZ*I_ERI_D2z_S_Py_Py_a;
  Double I_ERI_D2x_Px_Pz_Py_a = I_ERI_F3x_S_Pz_Py_a+ABX*I_ERI_D2x_S_Pz_Py_a;
  Double I_ERI_Dxy_Px_Pz_Py_a = I_ERI_F2xy_S_Pz_Py_a+ABX*I_ERI_Dxy_S_Pz_Py_a;
  Double I_ERI_Dxz_Px_Pz_Py_a = I_ERI_F2xz_S_Pz_Py_a+ABX*I_ERI_Dxz_S_Pz_Py_a;
  Double I_ERI_D2y_Px_Pz_Py_a = I_ERI_Fx2y_S_Pz_Py_a+ABX*I_ERI_D2y_S_Pz_Py_a;
  Double I_ERI_Dyz_Px_Pz_Py_a = I_ERI_Fxyz_S_Pz_Py_a+ABX*I_ERI_Dyz_S_Pz_Py_a;
  Double I_ERI_D2z_Px_Pz_Py_a = I_ERI_Fx2z_S_Pz_Py_a+ABX*I_ERI_D2z_S_Pz_Py_a;
  Double I_ERI_D2x_Py_Pz_Py_a = I_ERI_F2xy_S_Pz_Py_a+ABY*I_ERI_D2x_S_Pz_Py_a;
  Double I_ERI_Dxy_Py_Pz_Py_a = I_ERI_Fx2y_S_Pz_Py_a+ABY*I_ERI_Dxy_S_Pz_Py_a;
  Double I_ERI_Dxz_Py_Pz_Py_a = I_ERI_Fxyz_S_Pz_Py_a+ABY*I_ERI_Dxz_S_Pz_Py_a;
  Double I_ERI_D2y_Py_Pz_Py_a = I_ERI_F3y_S_Pz_Py_a+ABY*I_ERI_D2y_S_Pz_Py_a;
  Double I_ERI_Dyz_Py_Pz_Py_a = I_ERI_F2yz_S_Pz_Py_a+ABY*I_ERI_Dyz_S_Pz_Py_a;
  Double I_ERI_D2z_Py_Pz_Py_a = I_ERI_Fy2z_S_Pz_Py_a+ABY*I_ERI_D2z_S_Pz_Py_a;
  Double I_ERI_D2x_Pz_Pz_Py_a = I_ERI_F2xz_S_Pz_Py_a+ABZ*I_ERI_D2x_S_Pz_Py_a;
  Double I_ERI_Dxy_Pz_Pz_Py_a = I_ERI_Fxyz_S_Pz_Py_a+ABZ*I_ERI_Dxy_S_Pz_Py_a;
  Double I_ERI_Dxz_Pz_Pz_Py_a = I_ERI_Fx2z_S_Pz_Py_a+ABZ*I_ERI_Dxz_S_Pz_Py_a;
  Double I_ERI_D2y_Pz_Pz_Py_a = I_ERI_F2yz_S_Pz_Py_a+ABZ*I_ERI_D2y_S_Pz_Py_a;
  Double I_ERI_Dyz_Pz_Pz_Py_a = I_ERI_Fy2z_S_Pz_Py_a+ABZ*I_ERI_Dyz_S_Pz_Py_a;
  Double I_ERI_D2z_Pz_Pz_Py_a = I_ERI_F3z_S_Pz_Py_a+ABZ*I_ERI_D2z_S_Pz_Py_a;
  Double I_ERI_D2x_Px_Px_Pz_a = I_ERI_F3x_S_Px_Pz_a+ABX*I_ERI_D2x_S_Px_Pz_a;
  Double I_ERI_Dxy_Px_Px_Pz_a = I_ERI_F2xy_S_Px_Pz_a+ABX*I_ERI_Dxy_S_Px_Pz_a;
  Double I_ERI_Dxz_Px_Px_Pz_a = I_ERI_F2xz_S_Px_Pz_a+ABX*I_ERI_Dxz_S_Px_Pz_a;
  Double I_ERI_D2y_Px_Px_Pz_a = I_ERI_Fx2y_S_Px_Pz_a+ABX*I_ERI_D2y_S_Px_Pz_a;
  Double I_ERI_Dyz_Px_Px_Pz_a = I_ERI_Fxyz_S_Px_Pz_a+ABX*I_ERI_Dyz_S_Px_Pz_a;
  Double I_ERI_D2z_Px_Px_Pz_a = I_ERI_Fx2z_S_Px_Pz_a+ABX*I_ERI_D2z_S_Px_Pz_a;
  Double I_ERI_D2x_Py_Px_Pz_a = I_ERI_F2xy_S_Px_Pz_a+ABY*I_ERI_D2x_S_Px_Pz_a;
  Double I_ERI_Dxy_Py_Px_Pz_a = I_ERI_Fx2y_S_Px_Pz_a+ABY*I_ERI_Dxy_S_Px_Pz_a;
  Double I_ERI_Dxz_Py_Px_Pz_a = I_ERI_Fxyz_S_Px_Pz_a+ABY*I_ERI_Dxz_S_Px_Pz_a;
  Double I_ERI_D2y_Py_Px_Pz_a = I_ERI_F3y_S_Px_Pz_a+ABY*I_ERI_D2y_S_Px_Pz_a;
  Double I_ERI_Dyz_Py_Px_Pz_a = I_ERI_F2yz_S_Px_Pz_a+ABY*I_ERI_Dyz_S_Px_Pz_a;
  Double I_ERI_D2z_Py_Px_Pz_a = I_ERI_Fy2z_S_Px_Pz_a+ABY*I_ERI_D2z_S_Px_Pz_a;
  Double I_ERI_D2x_Pz_Px_Pz_a = I_ERI_F2xz_S_Px_Pz_a+ABZ*I_ERI_D2x_S_Px_Pz_a;
  Double I_ERI_Dxy_Pz_Px_Pz_a = I_ERI_Fxyz_S_Px_Pz_a+ABZ*I_ERI_Dxy_S_Px_Pz_a;
  Double I_ERI_Dxz_Pz_Px_Pz_a = I_ERI_Fx2z_S_Px_Pz_a+ABZ*I_ERI_Dxz_S_Px_Pz_a;
  Double I_ERI_D2y_Pz_Px_Pz_a = I_ERI_F2yz_S_Px_Pz_a+ABZ*I_ERI_D2y_S_Px_Pz_a;
  Double I_ERI_Dyz_Pz_Px_Pz_a = I_ERI_Fy2z_S_Px_Pz_a+ABZ*I_ERI_Dyz_S_Px_Pz_a;
  Double I_ERI_D2z_Pz_Px_Pz_a = I_ERI_F3z_S_Px_Pz_a+ABZ*I_ERI_D2z_S_Px_Pz_a;
  Double I_ERI_D2x_Px_Py_Pz_a = I_ERI_F3x_S_Py_Pz_a+ABX*I_ERI_D2x_S_Py_Pz_a;
  Double I_ERI_Dxy_Px_Py_Pz_a = I_ERI_F2xy_S_Py_Pz_a+ABX*I_ERI_Dxy_S_Py_Pz_a;
  Double I_ERI_Dxz_Px_Py_Pz_a = I_ERI_F2xz_S_Py_Pz_a+ABX*I_ERI_Dxz_S_Py_Pz_a;
  Double I_ERI_D2y_Px_Py_Pz_a = I_ERI_Fx2y_S_Py_Pz_a+ABX*I_ERI_D2y_S_Py_Pz_a;
  Double I_ERI_Dyz_Px_Py_Pz_a = I_ERI_Fxyz_S_Py_Pz_a+ABX*I_ERI_Dyz_S_Py_Pz_a;
  Double I_ERI_D2z_Px_Py_Pz_a = I_ERI_Fx2z_S_Py_Pz_a+ABX*I_ERI_D2z_S_Py_Pz_a;
  Double I_ERI_D2x_Py_Py_Pz_a = I_ERI_F2xy_S_Py_Pz_a+ABY*I_ERI_D2x_S_Py_Pz_a;
  Double I_ERI_Dxy_Py_Py_Pz_a = I_ERI_Fx2y_S_Py_Pz_a+ABY*I_ERI_Dxy_S_Py_Pz_a;
  Double I_ERI_Dxz_Py_Py_Pz_a = I_ERI_Fxyz_S_Py_Pz_a+ABY*I_ERI_Dxz_S_Py_Pz_a;
  Double I_ERI_D2y_Py_Py_Pz_a = I_ERI_F3y_S_Py_Pz_a+ABY*I_ERI_D2y_S_Py_Pz_a;
  Double I_ERI_Dyz_Py_Py_Pz_a = I_ERI_F2yz_S_Py_Pz_a+ABY*I_ERI_Dyz_S_Py_Pz_a;
  Double I_ERI_D2z_Py_Py_Pz_a = I_ERI_Fy2z_S_Py_Pz_a+ABY*I_ERI_D2z_S_Py_Pz_a;
  Double I_ERI_D2x_Pz_Py_Pz_a = I_ERI_F2xz_S_Py_Pz_a+ABZ*I_ERI_D2x_S_Py_Pz_a;
  Double I_ERI_Dxy_Pz_Py_Pz_a = I_ERI_Fxyz_S_Py_Pz_a+ABZ*I_ERI_Dxy_S_Py_Pz_a;
  Double I_ERI_Dxz_Pz_Py_Pz_a = I_ERI_Fx2z_S_Py_Pz_a+ABZ*I_ERI_Dxz_S_Py_Pz_a;
  Double I_ERI_D2y_Pz_Py_Pz_a = I_ERI_F2yz_S_Py_Pz_a+ABZ*I_ERI_D2y_S_Py_Pz_a;
  Double I_ERI_Dyz_Pz_Py_Pz_a = I_ERI_Fy2z_S_Py_Pz_a+ABZ*I_ERI_Dyz_S_Py_Pz_a;
  Double I_ERI_D2z_Pz_Py_Pz_a = I_ERI_F3z_S_Py_Pz_a+ABZ*I_ERI_D2z_S_Py_Pz_a;
  Double I_ERI_D2x_Px_Pz_Pz_a = I_ERI_F3x_S_Pz_Pz_a+ABX*I_ERI_D2x_S_Pz_Pz_a;
  Double I_ERI_Dxy_Px_Pz_Pz_a = I_ERI_F2xy_S_Pz_Pz_a+ABX*I_ERI_Dxy_S_Pz_Pz_a;
  Double I_ERI_Dxz_Px_Pz_Pz_a = I_ERI_F2xz_S_Pz_Pz_a+ABX*I_ERI_Dxz_S_Pz_Pz_a;
  Double I_ERI_D2y_Px_Pz_Pz_a = I_ERI_Fx2y_S_Pz_Pz_a+ABX*I_ERI_D2y_S_Pz_Pz_a;
  Double I_ERI_Dyz_Px_Pz_Pz_a = I_ERI_Fxyz_S_Pz_Pz_a+ABX*I_ERI_Dyz_S_Pz_Pz_a;
  Double I_ERI_D2z_Px_Pz_Pz_a = I_ERI_Fx2z_S_Pz_Pz_a+ABX*I_ERI_D2z_S_Pz_Pz_a;
  Double I_ERI_D2x_Py_Pz_Pz_a = I_ERI_F2xy_S_Pz_Pz_a+ABY*I_ERI_D2x_S_Pz_Pz_a;
  Double I_ERI_Dxy_Py_Pz_Pz_a = I_ERI_Fx2y_S_Pz_Pz_a+ABY*I_ERI_Dxy_S_Pz_Pz_a;
  Double I_ERI_Dxz_Py_Pz_Pz_a = I_ERI_Fxyz_S_Pz_Pz_a+ABY*I_ERI_Dxz_S_Pz_Pz_a;
  Double I_ERI_D2y_Py_Pz_Pz_a = I_ERI_F3y_S_Pz_Pz_a+ABY*I_ERI_D2y_S_Pz_Pz_a;
  Double I_ERI_Dyz_Py_Pz_Pz_a = I_ERI_F2yz_S_Pz_Pz_a+ABY*I_ERI_Dyz_S_Pz_Pz_a;
  Double I_ERI_D2z_Py_Pz_Pz_a = I_ERI_Fy2z_S_Pz_Pz_a+ABY*I_ERI_D2z_S_Pz_Pz_a;
  Double I_ERI_D2x_Pz_Pz_Pz_a = I_ERI_F2xz_S_Pz_Pz_a+ABZ*I_ERI_D2x_S_Pz_Pz_a;
  Double I_ERI_Dxy_Pz_Pz_Pz_a = I_ERI_Fxyz_S_Pz_Pz_a+ABZ*I_ERI_Dxy_S_Pz_Pz_a;
  Double I_ERI_Dxz_Pz_Pz_Pz_a = I_ERI_Fx2z_S_Pz_Pz_a+ABZ*I_ERI_Dxz_S_Pz_Pz_a;
  Double I_ERI_D2y_Pz_Pz_Pz_a = I_ERI_F2yz_S_Pz_Pz_a+ABZ*I_ERI_D2y_S_Pz_Pz_a;
  Double I_ERI_Dyz_Pz_Pz_Pz_a = I_ERI_Fy2z_S_Pz_Pz_a+ABZ*I_ERI_Dyz_S_Pz_Pz_a;
  Double I_ERI_D2z_Pz_Pz_Pz_a = I_ERI_F3z_S_Pz_Pz_a+ABZ*I_ERI_D2z_S_Pz_Pz_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_P_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_Px_b = I_ERI_D2x_S_Px_Px_b+ABX*I_ERI_Px_S_Px_Px_b;
  Double I_ERI_Py_Px_Px_Px_b = I_ERI_Dxy_S_Px_Px_b+ABX*I_ERI_Py_S_Px_Px_b;
  Double I_ERI_Pz_Px_Px_Px_b = I_ERI_Dxz_S_Px_Px_b+ABX*I_ERI_Pz_S_Px_Px_b;
  Double I_ERI_Px_Py_Px_Px_b = I_ERI_Dxy_S_Px_Px_b+ABY*I_ERI_Px_S_Px_Px_b;
  Double I_ERI_Py_Py_Px_Px_b = I_ERI_D2y_S_Px_Px_b+ABY*I_ERI_Py_S_Px_Px_b;
  Double I_ERI_Pz_Py_Px_Px_b = I_ERI_Dyz_S_Px_Px_b+ABY*I_ERI_Pz_S_Px_Px_b;
  Double I_ERI_Px_Pz_Px_Px_b = I_ERI_Dxz_S_Px_Px_b+ABZ*I_ERI_Px_S_Px_Px_b;
  Double I_ERI_Py_Pz_Px_Px_b = I_ERI_Dyz_S_Px_Px_b+ABZ*I_ERI_Py_S_Px_Px_b;
  Double I_ERI_Pz_Pz_Px_Px_b = I_ERI_D2z_S_Px_Px_b+ABZ*I_ERI_Pz_S_Px_Px_b;
  Double I_ERI_Px_Px_Py_Px_b = I_ERI_D2x_S_Py_Px_b+ABX*I_ERI_Px_S_Py_Px_b;
  Double I_ERI_Py_Px_Py_Px_b = I_ERI_Dxy_S_Py_Px_b+ABX*I_ERI_Py_S_Py_Px_b;
  Double I_ERI_Pz_Px_Py_Px_b = I_ERI_Dxz_S_Py_Px_b+ABX*I_ERI_Pz_S_Py_Px_b;
  Double I_ERI_Px_Py_Py_Px_b = I_ERI_Dxy_S_Py_Px_b+ABY*I_ERI_Px_S_Py_Px_b;
  Double I_ERI_Py_Py_Py_Px_b = I_ERI_D2y_S_Py_Px_b+ABY*I_ERI_Py_S_Py_Px_b;
  Double I_ERI_Pz_Py_Py_Px_b = I_ERI_Dyz_S_Py_Px_b+ABY*I_ERI_Pz_S_Py_Px_b;
  Double I_ERI_Px_Pz_Py_Px_b = I_ERI_Dxz_S_Py_Px_b+ABZ*I_ERI_Px_S_Py_Px_b;
  Double I_ERI_Py_Pz_Py_Px_b = I_ERI_Dyz_S_Py_Px_b+ABZ*I_ERI_Py_S_Py_Px_b;
  Double I_ERI_Pz_Pz_Py_Px_b = I_ERI_D2z_S_Py_Px_b+ABZ*I_ERI_Pz_S_Py_Px_b;
  Double I_ERI_Px_Px_Pz_Px_b = I_ERI_D2x_S_Pz_Px_b+ABX*I_ERI_Px_S_Pz_Px_b;
  Double I_ERI_Py_Px_Pz_Px_b = I_ERI_Dxy_S_Pz_Px_b+ABX*I_ERI_Py_S_Pz_Px_b;
  Double I_ERI_Pz_Px_Pz_Px_b = I_ERI_Dxz_S_Pz_Px_b+ABX*I_ERI_Pz_S_Pz_Px_b;
  Double I_ERI_Px_Py_Pz_Px_b = I_ERI_Dxy_S_Pz_Px_b+ABY*I_ERI_Px_S_Pz_Px_b;
  Double I_ERI_Py_Py_Pz_Px_b = I_ERI_D2y_S_Pz_Px_b+ABY*I_ERI_Py_S_Pz_Px_b;
  Double I_ERI_Pz_Py_Pz_Px_b = I_ERI_Dyz_S_Pz_Px_b+ABY*I_ERI_Pz_S_Pz_Px_b;
  Double I_ERI_Px_Pz_Pz_Px_b = I_ERI_Dxz_S_Pz_Px_b+ABZ*I_ERI_Px_S_Pz_Px_b;
  Double I_ERI_Py_Pz_Pz_Px_b = I_ERI_Dyz_S_Pz_Px_b+ABZ*I_ERI_Py_S_Pz_Px_b;
  Double I_ERI_Pz_Pz_Pz_Px_b = I_ERI_D2z_S_Pz_Px_b+ABZ*I_ERI_Pz_S_Pz_Px_b;
  Double I_ERI_Px_Px_Px_Py_b = I_ERI_D2x_S_Px_Py_b+ABX*I_ERI_Px_S_Px_Py_b;
  Double I_ERI_Py_Px_Px_Py_b = I_ERI_Dxy_S_Px_Py_b+ABX*I_ERI_Py_S_Px_Py_b;
  Double I_ERI_Pz_Px_Px_Py_b = I_ERI_Dxz_S_Px_Py_b+ABX*I_ERI_Pz_S_Px_Py_b;
  Double I_ERI_Px_Py_Px_Py_b = I_ERI_Dxy_S_Px_Py_b+ABY*I_ERI_Px_S_Px_Py_b;
  Double I_ERI_Py_Py_Px_Py_b = I_ERI_D2y_S_Px_Py_b+ABY*I_ERI_Py_S_Px_Py_b;
  Double I_ERI_Pz_Py_Px_Py_b = I_ERI_Dyz_S_Px_Py_b+ABY*I_ERI_Pz_S_Px_Py_b;
  Double I_ERI_Px_Pz_Px_Py_b = I_ERI_Dxz_S_Px_Py_b+ABZ*I_ERI_Px_S_Px_Py_b;
  Double I_ERI_Py_Pz_Px_Py_b = I_ERI_Dyz_S_Px_Py_b+ABZ*I_ERI_Py_S_Px_Py_b;
  Double I_ERI_Pz_Pz_Px_Py_b = I_ERI_D2z_S_Px_Py_b+ABZ*I_ERI_Pz_S_Px_Py_b;
  Double I_ERI_Px_Px_Py_Py_b = I_ERI_D2x_S_Py_Py_b+ABX*I_ERI_Px_S_Py_Py_b;
  Double I_ERI_Py_Px_Py_Py_b = I_ERI_Dxy_S_Py_Py_b+ABX*I_ERI_Py_S_Py_Py_b;
  Double I_ERI_Pz_Px_Py_Py_b = I_ERI_Dxz_S_Py_Py_b+ABX*I_ERI_Pz_S_Py_Py_b;
  Double I_ERI_Px_Py_Py_Py_b = I_ERI_Dxy_S_Py_Py_b+ABY*I_ERI_Px_S_Py_Py_b;
  Double I_ERI_Py_Py_Py_Py_b = I_ERI_D2y_S_Py_Py_b+ABY*I_ERI_Py_S_Py_Py_b;
  Double I_ERI_Pz_Py_Py_Py_b = I_ERI_Dyz_S_Py_Py_b+ABY*I_ERI_Pz_S_Py_Py_b;
  Double I_ERI_Px_Pz_Py_Py_b = I_ERI_Dxz_S_Py_Py_b+ABZ*I_ERI_Px_S_Py_Py_b;
  Double I_ERI_Py_Pz_Py_Py_b = I_ERI_Dyz_S_Py_Py_b+ABZ*I_ERI_Py_S_Py_Py_b;
  Double I_ERI_Pz_Pz_Py_Py_b = I_ERI_D2z_S_Py_Py_b+ABZ*I_ERI_Pz_S_Py_Py_b;
  Double I_ERI_Px_Px_Pz_Py_b = I_ERI_D2x_S_Pz_Py_b+ABX*I_ERI_Px_S_Pz_Py_b;
  Double I_ERI_Py_Px_Pz_Py_b = I_ERI_Dxy_S_Pz_Py_b+ABX*I_ERI_Py_S_Pz_Py_b;
  Double I_ERI_Pz_Px_Pz_Py_b = I_ERI_Dxz_S_Pz_Py_b+ABX*I_ERI_Pz_S_Pz_Py_b;
  Double I_ERI_Px_Py_Pz_Py_b = I_ERI_Dxy_S_Pz_Py_b+ABY*I_ERI_Px_S_Pz_Py_b;
  Double I_ERI_Py_Py_Pz_Py_b = I_ERI_D2y_S_Pz_Py_b+ABY*I_ERI_Py_S_Pz_Py_b;
  Double I_ERI_Pz_Py_Pz_Py_b = I_ERI_Dyz_S_Pz_Py_b+ABY*I_ERI_Pz_S_Pz_Py_b;
  Double I_ERI_Px_Pz_Pz_Py_b = I_ERI_Dxz_S_Pz_Py_b+ABZ*I_ERI_Px_S_Pz_Py_b;
  Double I_ERI_Py_Pz_Pz_Py_b = I_ERI_Dyz_S_Pz_Py_b+ABZ*I_ERI_Py_S_Pz_Py_b;
  Double I_ERI_Pz_Pz_Pz_Py_b = I_ERI_D2z_S_Pz_Py_b+ABZ*I_ERI_Pz_S_Pz_Py_b;
  Double I_ERI_Px_Px_Px_Pz_b = I_ERI_D2x_S_Px_Pz_b+ABX*I_ERI_Px_S_Px_Pz_b;
  Double I_ERI_Py_Px_Px_Pz_b = I_ERI_Dxy_S_Px_Pz_b+ABX*I_ERI_Py_S_Px_Pz_b;
  Double I_ERI_Pz_Px_Px_Pz_b = I_ERI_Dxz_S_Px_Pz_b+ABX*I_ERI_Pz_S_Px_Pz_b;
  Double I_ERI_Px_Py_Px_Pz_b = I_ERI_Dxy_S_Px_Pz_b+ABY*I_ERI_Px_S_Px_Pz_b;
  Double I_ERI_Py_Py_Px_Pz_b = I_ERI_D2y_S_Px_Pz_b+ABY*I_ERI_Py_S_Px_Pz_b;
  Double I_ERI_Pz_Py_Px_Pz_b = I_ERI_Dyz_S_Px_Pz_b+ABY*I_ERI_Pz_S_Px_Pz_b;
  Double I_ERI_Px_Pz_Px_Pz_b = I_ERI_Dxz_S_Px_Pz_b+ABZ*I_ERI_Px_S_Px_Pz_b;
  Double I_ERI_Py_Pz_Px_Pz_b = I_ERI_Dyz_S_Px_Pz_b+ABZ*I_ERI_Py_S_Px_Pz_b;
  Double I_ERI_Pz_Pz_Px_Pz_b = I_ERI_D2z_S_Px_Pz_b+ABZ*I_ERI_Pz_S_Px_Pz_b;
  Double I_ERI_Px_Px_Py_Pz_b = I_ERI_D2x_S_Py_Pz_b+ABX*I_ERI_Px_S_Py_Pz_b;
  Double I_ERI_Py_Px_Py_Pz_b = I_ERI_Dxy_S_Py_Pz_b+ABX*I_ERI_Py_S_Py_Pz_b;
  Double I_ERI_Pz_Px_Py_Pz_b = I_ERI_Dxz_S_Py_Pz_b+ABX*I_ERI_Pz_S_Py_Pz_b;
  Double I_ERI_Px_Py_Py_Pz_b = I_ERI_Dxy_S_Py_Pz_b+ABY*I_ERI_Px_S_Py_Pz_b;
  Double I_ERI_Py_Py_Py_Pz_b = I_ERI_D2y_S_Py_Pz_b+ABY*I_ERI_Py_S_Py_Pz_b;
  Double I_ERI_Pz_Py_Py_Pz_b = I_ERI_Dyz_S_Py_Pz_b+ABY*I_ERI_Pz_S_Py_Pz_b;
  Double I_ERI_Px_Pz_Py_Pz_b = I_ERI_Dxz_S_Py_Pz_b+ABZ*I_ERI_Px_S_Py_Pz_b;
  Double I_ERI_Py_Pz_Py_Pz_b = I_ERI_Dyz_S_Py_Pz_b+ABZ*I_ERI_Py_S_Py_Pz_b;
  Double I_ERI_Pz_Pz_Py_Pz_b = I_ERI_D2z_S_Py_Pz_b+ABZ*I_ERI_Pz_S_Py_Pz_b;
  Double I_ERI_Px_Px_Pz_Pz_b = I_ERI_D2x_S_Pz_Pz_b+ABX*I_ERI_Px_S_Pz_Pz_b;
  Double I_ERI_Py_Px_Pz_Pz_b = I_ERI_Dxy_S_Pz_Pz_b+ABX*I_ERI_Py_S_Pz_Pz_b;
  Double I_ERI_Pz_Px_Pz_Pz_b = I_ERI_Dxz_S_Pz_Pz_b+ABX*I_ERI_Pz_S_Pz_Pz_b;
  Double I_ERI_Px_Py_Pz_Pz_b = I_ERI_Dxy_S_Pz_Pz_b+ABY*I_ERI_Px_S_Pz_Pz_b;
  Double I_ERI_Py_Py_Pz_Pz_b = I_ERI_D2y_S_Pz_Pz_b+ABY*I_ERI_Py_S_Pz_Pz_b;
  Double I_ERI_Pz_Py_Pz_Pz_b = I_ERI_Dyz_S_Pz_Pz_b+ABY*I_ERI_Pz_S_Pz_Pz_b;
  Double I_ERI_Px_Pz_Pz_Pz_b = I_ERI_Dxz_S_Pz_Pz_b+ABZ*I_ERI_Px_S_Pz_Pz_b;
  Double I_ERI_Py_Pz_Pz_Pz_b = I_ERI_Dyz_S_Pz_Pz_b+ABZ*I_ERI_Py_S_Pz_Pz_b;
  Double I_ERI_Pz_Pz_Pz_Pz_b = I_ERI_D2z_S_Pz_Pz_b+ABZ*I_ERI_Pz_S_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_b
   * RHS shell quartet name: SQ_ERI_D_S_P_P_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_Px_b = I_ERI_F3x_S_Px_Px_b+ABX*I_ERI_D2x_S_Px_Px_b;
  Double I_ERI_Dxy_Px_Px_Px_b = I_ERI_F2xy_S_Px_Px_b+ABX*I_ERI_Dxy_S_Px_Px_b;
  Double I_ERI_Dxz_Px_Px_Px_b = I_ERI_F2xz_S_Px_Px_b+ABX*I_ERI_Dxz_S_Px_Px_b;
  Double I_ERI_D2y_Px_Px_Px_b = I_ERI_Fx2y_S_Px_Px_b+ABX*I_ERI_D2y_S_Px_Px_b;
  Double I_ERI_Dyz_Px_Px_Px_b = I_ERI_Fxyz_S_Px_Px_b+ABX*I_ERI_Dyz_S_Px_Px_b;
  Double I_ERI_D2z_Px_Px_Px_b = I_ERI_Fx2z_S_Px_Px_b+ABX*I_ERI_D2z_S_Px_Px_b;
  Double I_ERI_Dxy_Py_Px_Px_b = I_ERI_Fx2y_S_Px_Px_b+ABY*I_ERI_Dxy_S_Px_Px_b;
  Double I_ERI_Dxz_Py_Px_Px_b = I_ERI_Fxyz_S_Px_Px_b+ABY*I_ERI_Dxz_S_Px_Px_b;
  Double I_ERI_D2y_Py_Px_Px_b = I_ERI_F3y_S_Px_Px_b+ABY*I_ERI_D2y_S_Px_Px_b;
  Double I_ERI_Dyz_Py_Px_Px_b = I_ERI_F2yz_S_Px_Px_b+ABY*I_ERI_Dyz_S_Px_Px_b;
  Double I_ERI_D2z_Py_Px_Px_b = I_ERI_Fy2z_S_Px_Px_b+ABY*I_ERI_D2z_S_Px_Px_b;
  Double I_ERI_Dxz_Pz_Px_Px_b = I_ERI_Fx2z_S_Px_Px_b+ABZ*I_ERI_Dxz_S_Px_Px_b;
  Double I_ERI_Dyz_Pz_Px_Px_b = I_ERI_Fy2z_S_Px_Px_b+ABZ*I_ERI_Dyz_S_Px_Px_b;
  Double I_ERI_D2z_Pz_Px_Px_b = I_ERI_F3z_S_Px_Px_b+ABZ*I_ERI_D2z_S_Px_Px_b;
  Double I_ERI_D2x_Px_Py_Px_b = I_ERI_F3x_S_Py_Px_b+ABX*I_ERI_D2x_S_Py_Px_b;
  Double I_ERI_Dxy_Px_Py_Px_b = I_ERI_F2xy_S_Py_Px_b+ABX*I_ERI_Dxy_S_Py_Px_b;
  Double I_ERI_Dxz_Px_Py_Px_b = I_ERI_F2xz_S_Py_Px_b+ABX*I_ERI_Dxz_S_Py_Px_b;
  Double I_ERI_D2y_Px_Py_Px_b = I_ERI_Fx2y_S_Py_Px_b+ABX*I_ERI_D2y_S_Py_Px_b;
  Double I_ERI_Dyz_Px_Py_Px_b = I_ERI_Fxyz_S_Py_Px_b+ABX*I_ERI_Dyz_S_Py_Px_b;
  Double I_ERI_D2z_Px_Py_Px_b = I_ERI_Fx2z_S_Py_Px_b+ABX*I_ERI_D2z_S_Py_Px_b;
  Double I_ERI_Dxy_Py_Py_Px_b = I_ERI_Fx2y_S_Py_Px_b+ABY*I_ERI_Dxy_S_Py_Px_b;
  Double I_ERI_Dxz_Py_Py_Px_b = I_ERI_Fxyz_S_Py_Px_b+ABY*I_ERI_Dxz_S_Py_Px_b;
  Double I_ERI_D2y_Py_Py_Px_b = I_ERI_F3y_S_Py_Px_b+ABY*I_ERI_D2y_S_Py_Px_b;
  Double I_ERI_Dyz_Py_Py_Px_b = I_ERI_F2yz_S_Py_Px_b+ABY*I_ERI_Dyz_S_Py_Px_b;
  Double I_ERI_D2z_Py_Py_Px_b = I_ERI_Fy2z_S_Py_Px_b+ABY*I_ERI_D2z_S_Py_Px_b;
  Double I_ERI_Dxz_Pz_Py_Px_b = I_ERI_Fx2z_S_Py_Px_b+ABZ*I_ERI_Dxz_S_Py_Px_b;
  Double I_ERI_Dyz_Pz_Py_Px_b = I_ERI_Fy2z_S_Py_Px_b+ABZ*I_ERI_Dyz_S_Py_Px_b;
  Double I_ERI_D2z_Pz_Py_Px_b = I_ERI_F3z_S_Py_Px_b+ABZ*I_ERI_D2z_S_Py_Px_b;
  Double I_ERI_D2x_Px_Pz_Px_b = I_ERI_F3x_S_Pz_Px_b+ABX*I_ERI_D2x_S_Pz_Px_b;
  Double I_ERI_Dxy_Px_Pz_Px_b = I_ERI_F2xy_S_Pz_Px_b+ABX*I_ERI_Dxy_S_Pz_Px_b;
  Double I_ERI_Dxz_Px_Pz_Px_b = I_ERI_F2xz_S_Pz_Px_b+ABX*I_ERI_Dxz_S_Pz_Px_b;
  Double I_ERI_D2y_Px_Pz_Px_b = I_ERI_Fx2y_S_Pz_Px_b+ABX*I_ERI_D2y_S_Pz_Px_b;
  Double I_ERI_Dyz_Px_Pz_Px_b = I_ERI_Fxyz_S_Pz_Px_b+ABX*I_ERI_Dyz_S_Pz_Px_b;
  Double I_ERI_D2z_Px_Pz_Px_b = I_ERI_Fx2z_S_Pz_Px_b+ABX*I_ERI_D2z_S_Pz_Px_b;
  Double I_ERI_Dxy_Py_Pz_Px_b = I_ERI_Fx2y_S_Pz_Px_b+ABY*I_ERI_Dxy_S_Pz_Px_b;
  Double I_ERI_Dxz_Py_Pz_Px_b = I_ERI_Fxyz_S_Pz_Px_b+ABY*I_ERI_Dxz_S_Pz_Px_b;
  Double I_ERI_D2y_Py_Pz_Px_b = I_ERI_F3y_S_Pz_Px_b+ABY*I_ERI_D2y_S_Pz_Px_b;
  Double I_ERI_Dyz_Py_Pz_Px_b = I_ERI_F2yz_S_Pz_Px_b+ABY*I_ERI_Dyz_S_Pz_Px_b;
  Double I_ERI_D2z_Py_Pz_Px_b = I_ERI_Fy2z_S_Pz_Px_b+ABY*I_ERI_D2z_S_Pz_Px_b;
  Double I_ERI_Dxz_Pz_Pz_Px_b = I_ERI_Fx2z_S_Pz_Px_b+ABZ*I_ERI_Dxz_S_Pz_Px_b;
  Double I_ERI_Dyz_Pz_Pz_Px_b = I_ERI_Fy2z_S_Pz_Px_b+ABZ*I_ERI_Dyz_S_Pz_Px_b;
  Double I_ERI_D2z_Pz_Pz_Px_b = I_ERI_F3z_S_Pz_Px_b+ABZ*I_ERI_D2z_S_Pz_Px_b;
  Double I_ERI_D2x_Px_Px_Py_b = I_ERI_F3x_S_Px_Py_b+ABX*I_ERI_D2x_S_Px_Py_b;
  Double I_ERI_Dxy_Px_Px_Py_b = I_ERI_F2xy_S_Px_Py_b+ABX*I_ERI_Dxy_S_Px_Py_b;
  Double I_ERI_Dxz_Px_Px_Py_b = I_ERI_F2xz_S_Px_Py_b+ABX*I_ERI_Dxz_S_Px_Py_b;
  Double I_ERI_D2y_Px_Px_Py_b = I_ERI_Fx2y_S_Px_Py_b+ABX*I_ERI_D2y_S_Px_Py_b;
  Double I_ERI_Dyz_Px_Px_Py_b = I_ERI_Fxyz_S_Px_Py_b+ABX*I_ERI_Dyz_S_Px_Py_b;
  Double I_ERI_D2z_Px_Px_Py_b = I_ERI_Fx2z_S_Px_Py_b+ABX*I_ERI_D2z_S_Px_Py_b;
  Double I_ERI_Dxy_Py_Px_Py_b = I_ERI_Fx2y_S_Px_Py_b+ABY*I_ERI_Dxy_S_Px_Py_b;
  Double I_ERI_Dxz_Py_Px_Py_b = I_ERI_Fxyz_S_Px_Py_b+ABY*I_ERI_Dxz_S_Px_Py_b;
  Double I_ERI_D2y_Py_Px_Py_b = I_ERI_F3y_S_Px_Py_b+ABY*I_ERI_D2y_S_Px_Py_b;
  Double I_ERI_Dyz_Py_Px_Py_b = I_ERI_F2yz_S_Px_Py_b+ABY*I_ERI_Dyz_S_Px_Py_b;
  Double I_ERI_D2z_Py_Px_Py_b = I_ERI_Fy2z_S_Px_Py_b+ABY*I_ERI_D2z_S_Px_Py_b;
  Double I_ERI_Dxz_Pz_Px_Py_b = I_ERI_Fx2z_S_Px_Py_b+ABZ*I_ERI_Dxz_S_Px_Py_b;
  Double I_ERI_Dyz_Pz_Px_Py_b = I_ERI_Fy2z_S_Px_Py_b+ABZ*I_ERI_Dyz_S_Px_Py_b;
  Double I_ERI_D2z_Pz_Px_Py_b = I_ERI_F3z_S_Px_Py_b+ABZ*I_ERI_D2z_S_Px_Py_b;
  Double I_ERI_D2x_Px_Py_Py_b = I_ERI_F3x_S_Py_Py_b+ABX*I_ERI_D2x_S_Py_Py_b;
  Double I_ERI_Dxy_Px_Py_Py_b = I_ERI_F2xy_S_Py_Py_b+ABX*I_ERI_Dxy_S_Py_Py_b;
  Double I_ERI_Dxz_Px_Py_Py_b = I_ERI_F2xz_S_Py_Py_b+ABX*I_ERI_Dxz_S_Py_Py_b;
  Double I_ERI_D2y_Px_Py_Py_b = I_ERI_Fx2y_S_Py_Py_b+ABX*I_ERI_D2y_S_Py_Py_b;
  Double I_ERI_Dyz_Px_Py_Py_b = I_ERI_Fxyz_S_Py_Py_b+ABX*I_ERI_Dyz_S_Py_Py_b;
  Double I_ERI_D2z_Px_Py_Py_b = I_ERI_Fx2z_S_Py_Py_b+ABX*I_ERI_D2z_S_Py_Py_b;
  Double I_ERI_Dxy_Py_Py_Py_b = I_ERI_Fx2y_S_Py_Py_b+ABY*I_ERI_Dxy_S_Py_Py_b;
  Double I_ERI_Dxz_Py_Py_Py_b = I_ERI_Fxyz_S_Py_Py_b+ABY*I_ERI_Dxz_S_Py_Py_b;
  Double I_ERI_D2y_Py_Py_Py_b = I_ERI_F3y_S_Py_Py_b+ABY*I_ERI_D2y_S_Py_Py_b;
  Double I_ERI_Dyz_Py_Py_Py_b = I_ERI_F2yz_S_Py_Py_b+ABY*I_ERI_Dyz_S_Py_Py_b;
  Double I_ERI_D2z_Py_Py_Py_b = I_ERI_Fy2z_S_Py_Py_b+ABY*I_ERI_D2z_S_Py_Py_b;
  Double I_ERI_Dxz_Pz_Py_Py_b = I_ERI_Fx2z_S_Py_Py_b+ABZ*I_ERI_Dxz_S_Py_Py_b;
  Double I_ERI_Dyz_Pz_Py_Py_b = I_ERI_Fy2z_S_Py_Py_b+ABZ*I_ERI_Dyz_S_Py_Py_b;
  Double I_ERI_D2z_Pz_Py_Py_b = I_ERI_F3z_S_Py_Py_b+ABZ*I_ERI_D2z_S_Py_Py_b;
  Double I_ERI_D2x_Px_Pz_Py_b = I_ERI_F3x_S_Pz_Py_b+ABX*I_ERI_D2x_S_Pz_Py_b;
  Double I_ERI_Dxy_Px_Pz_Py_b = I_ERI_F2xy_S_Pz_Py_b+ABX*I_ERI_Dxy_S_Pz_Py_b;
  Double I_ERI_Dxz_Px_Pz_Py_b = I_ERI_F2xz_S_Pz_Py_b+ABX*I_ERI_Dxz_S_Pz_Py_b;
  Double I_ERI_D2y_Px_Pz_Py_b = I_ERI_Fx2y_S_Pz_Py_b+ABX*I_ERI_D2y_S_Pz_Py_b;
  Double I_ERI_Dyz_Px_Pz_Py_b = I_ERI_Fxyz_S_Pz_Py_b+ABX*I_ERI_Dyz_S_Pz_Py_b;
  Double I_ERI_D2z_Px_Pz_Py_b = I_ERI_Fx2z_S_Pz_Py_b+ABX*I_ERI_D2z_S_Pz_Py_b;
  Double I_ERI_Dxy_Py_Pz_Py_b = I_ERI_Fx2y_S_Pz_Py_b+ABY*I_ERI_Dxy_S_Pz_Py_b;
  Double I_ERI_Dxz_Py_Pz_Py_b = I_ERI_Fxyz_S_Pz_Py_b+ABY*I_ERI_Dxz_S_Pz_Py_b;
  Double I_ERI_D2y_Py_Pz_Py_b = I_ERI_F3y_S_Pz_Py_b+ABY*I_ERI_D2y_S_Pz_Py_b;
  Double I_ERI_Dyz_Py_Pz_Py_b = I_ERI_F2yz_S_Pz_Py_b+ABY*I_ERI_Dyz_S_Pz_Py_b;
  Double I_ERI_D2z_Py_Pz_Py_b = I_ERI_Fy2z_S_Pz_Py_b+ABY*I_ERI_D2z_S_Pz_Py_b;
  Double I_ERI_Dxz_Pz_Pz_Py_b = I_ERI_Fx2z_S_Pz_Py_b+ABZ*I_ERI_Dxz_S_Pz_Py_b;
  Double I_ERI_Dyz_Pz_Pz_Py_b = I_ERI_Fy2z_S_Pz_Py_b+ABZ*I_ERI_Dyz_S_Pz_Py_b;
  Double I_ERI_D2z_Pz_Pz_Py_b = I_ERI_F3z_S_Pz_Py_b+ABZ*I_ERI_D2z_S_Pz_Py_b;
  Double I_ERI_D2x_Px_Px_Pz_b = I_ERI_F3x_S_Px_Pz_b+ABX*I_ERI_D2x_S_Px_Pz_b;
  Double I_ERI_Dxy_Px_Px_Pz_b = I_ERI_F2xy_S_Px_Pz_b+ABX*I_ERI_Dxy_S_Px_Pz_b;
  Double I_ERI_Dxz_Px_Px_Pz_b = I_ERI_F2xz_S_Px_Pz_b+ABX*I_ERI_Dxz_S_Px_Pz_b;
  Double I_ERI_D2y_Px_Px_Pz_b = I_ERI_Fx2y_S_Px_Pz_b+ABX*I_ERI_D2y_S_Px_Pz_b;
  Double I_ERI_Dyz_Px_Px_Pz_b = I_ERI_Fxyz_S_Px_Pz_b+ABX*I_ERI_Dyz_S_Px_Pz_b;
  Double I_ERI_D2z_Px_Px_Pz_b = I_ERI_Fx2z_S_Px_Pz_b+ABX*I_ERI_D2z_S_Px_Pz_b;
  Double I_ERI_Dxy_Py_Px_Pz_b = I_ERI_Fx2y_S_Px_Pz_b+ABY*I_ERI_Dxy_S_Px_Pz_b;
  Double I_ERI_Dxz_Py_Px_Pz_b = I_ERI_Fxyz_S_Px_Pz_b+ABY*I_ERI_Dxz_S_Px_Pz_b;
  Double I_ERI_D2y_Py_Px_Pz_b = I_ERI_F3y_S_Px_Pz_b+ABY*I_ERI_D2y_S_Px_Pz_b;
  Double I_ERI_Dyz_Py_Px_Pz_b = I_ERI_F2yz_S_Px_Pz_b+ABY*I_ERI_Dyz_S_Px_Pz_b;
  Double I_ERI_D2z_Py_Px_Pz_b = I_ERI_Fy2z_S_Px_Pz_b+ABY*I_ERI_D2z_S_Px_Pz_b;
  Double I_ERI_Dxz_Pz_Px_Pz_b = I_ERI_Fx2z_S_Px_Pz_b+ABZ*I_ERI_Dxz_S_Px_Pz_b;
  Double I_ERI_Dyz_Pz_Px_Pz_b = I_ERI_Fy2z_S_Px_Pz_b+ABZ*I_ERI_Dyz_S_Px_Pz_b;
  Double I_ERI_D2z_Pz_Px_Pz_b = I_ERI_F3z_S_Px_Pz_b+ABZ*I_ERI_D2z_S_Px_Pz_b;
  Double I_ERI_D2x_Px_Py_Pz_b = I_ERI_F3x_S_Py_Pz_b+ABX*I_ERI_D2x_S_Py_Pz_b;
  Double I_ERI_Dxy_Px_Py_Pz_b = I_ERI_F2xy_S_Py_Pz_b+ABX*I_ERI_Dxy_S_Py_Pz_b;
  Double I_ERI_Dxz_Px_Py_Pz_b = I_ERI_F2xz_S_Py_Pz_b+ABX*I_ERI_Dxz_S_Py_Pz_b;
  Double I_ERI_D2y_Px_Py_Pz_b = I_ERI_Fx2y_S_Py_Pz_b+ABX*I_ERI_D2y_S_Py_Pz_b;
  Double I_ERI_Dyz_Px_Py_Pz_b = I_ERI_Fxyz_S_Py_Pz_b+ABX*I_ERI_Dyz_S_Py_Pz_b;
  Double I_ERI_D2z_Px_Py_Pz_b = I_ERI_Fx2z_S_Py_Pz_b+ABX*I_ERI_D2z_S_Py_Pz_b;
  Double I_ERI_Dxy_Py_Py_Pz_b = I_ERI_Fx2y_S_Py_Pz_b+ABY*I_ERI_Dxy_S_Py_Pz_b;
  Double I_ERI_Dxz_Py_Py_Pz_b = I_ERI_Fxyz_S_Py_Pz_b+ABY*I_ERI_Dxz_S_Py_Pz_b;
  Double I_ERI_D2y_Py_Py_Pz_b = I_ERI_F3y_S_Py_Pz_b+ABY*I_ERI_D2y_S_Py_Pz_b;
  Double I_ERI_Dyz_Py_Py_Pz_b = I_ERI_F2yz_S_Py_Pz_b+ABY*I_ERI_Dyz_S_Py_Pz_b;
  Double I_ERI_D2z_Py_Py_Pz_b = I_ERI_Fy2z_S_Py_Pz_b+ABY*I_ERI_D2z_S_Py_Pz_b;
  Double I_ERI_Dxz_Pz_Py_Pz_b = I_ERI_Fx2z_S_Py_Pz_b+ABZ*I_ERI_Dxz_S_Py_Pz_b;
  Double I_ERI_Dyz_Pz_Py_Pz_b = I_ERI_Fy2z_S_Py_Pz_b+ABZ*I_ERI_Dyz_S_Py_Pz_b;
  Double I_ERI_D2z_Pz_Py_Pz_b = I_ERI_F3z_S_Py_Pz_b+ABZ*I_ERI_D2z_S_Py_Pz_b;
  Double I_ERI_D2x_Px_Pz_Pz_b = I_ERI_F3x_S_Pz_Pz_b+ABX*I_ERI_D2x_S_Pz_Pz_b;
  Double I_ERI_Dxy_Px_Pz_Pz_b = I_ERI_F2xy_S_Pz_Pz_b+ABX*I_ERI_Dxy_S_Pz_Pz_b;
  Double I_ERI_Dxz_Px_Pz_Pz_b = I_ERI_F2xz_S_Pz_Pz_b+ABX*I_ERI_Dxz_S_Pz_Pz_b;
  Double I_ERI_D2y_Px_Pz_Pz_b = I_ERI_Fx2y_S_Pz_Pz_b+ABX*I_ERI_D2y_S_Pz_Pz_b;
  Double I_ERI_Dyz_Px_Pz_Pz_b = I_ERI_Fxyz_S_Pz_Pz_b+ABX*I_ERI_Dyz_S_Pz_Pz_b;
  Double I_ERI_D2z_Px_Pz_Pz_b = I_ERI_Fx2z_S_Pz_Pz_b+ABX*I_ERI_D2z_S_Pz_Pz_b;
  Double I_ERI_Dxy_Py_Pz_Pz_b = I_ERI_Fx2y_S_Pz_Pz_b+ABY*I_ERI_Dxy_S_Pz_Pz_b;
  Double I_ERI_Dxz_Py_Pz_Pz_b = I_ERI_Fxyz_S_Pz_Pz_b+ABY*I_ERI_Dxz_S_Pz_Pz_b;
  Double I_ERI_D2y_Py_Pz_Pz_b = I_ERI_F3y_S_Pz_Pz_b+ABY*I_ERI_D2y_S_Pz_Pz_b;
  Double I_ERI_Dyz_Py_Pz_Pz_b = I_ERI_F2yz_S_Pz_Pz_b+ABY*I_ERI_Dyz_S_Pz_Pz_b;
  Double I_ERI_D2z_Py_Pz_Pz_b = I_ERI_Fy2z_S_Pz_Pz_b+ABY*I_ERI_D2z_S_Pz_Pz_b;
  Double I_ERI_Dxz_Pz_Pz_Pz_b = I_ERI_Fx2z_S_Pz_Pz_b+ABZ*I_ERI_Dxz_S_Pz_Pz_b;
  Double I_ERI_Dyz_Pz_Pz_Pz_b = I_ERI_Fy2z_S_Pz_Pz_b+ABZ*I_ERI_Dyz_S_Pz_Pz_b;
  Double I_ERI_D2z_Pz_Pz_Pz_b = I_ERI_F3z_S_Pz_Pz_b+ABZ*I_ERI_D2z_S_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_P_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_b
   * RHS shell quartet name: SQ_ERI_P_P_P_P_b
   ************************************************************/
  Double I_ERI_Px_D2x_Px_Px_b = I_ERI_D2x_Px_Px_Px_b+ABX*I_ERI_Px_Px_Px_Px_b;
  Double I_ERI_Py_D2x_Px_Px_b = I_ERI_Dxy_Px_Px_Px_b+ABX*I_ERI_Py_Px_Px_Px_b;
  Double I_ERI_Pz_D2x_Px_Px_b = I_ERI_Dxz_Px_Px_Px_b+ABX*I_ERI_Pz_Px_Px_Px_b;
  Double I_ERI_Px_Dxy_Px_Px_b = I_ERI_Dxy_Px_Px_Px_b+ABY*I_ERI_Px_Px_Px_Px_b;
  Double I_ERI_Py_Dxy_Px_Px_b = I_ERI_D2y_Px_Px_Px_b+ABY*I_ERI_Py_Px_Px_Px_b;
  Double I_ERI_Pz_Dxy_Px_Px_b = I_ERI_Dyz_Px_Px_Px_b+ABY*I_ERI_Pz_Px_Px_Px_b;
  Double I_ERI_Px_Dxz_Px_Px_b = I_ERI_Dxz_Px_Px_Px_b+ABZ*I_ERI_Px_Px_Px_Px_b;
  Double I_ERI_Py_Dxz_Px_Px_b = I_ERI_Dyz_Px_Px_Px_b+ABZ*I_ERI_Py_Px_Px_Px_b;
  Double I_ERI_Pz_Dxz_Px_Px_b = I_ERI_D2z_Px_Px_Px_b+ABZ*I_ERI_Pz_Px_Px_Px_b;
  Double I_ERI_Px_D2y_Px_Px_b = I_ERI_Dxy_Py_Px_Px_b+ABY*I_ERI_Px_Py_Px_Px_b;
  Double I_ERI_Py_D2y_Px_Px_b = I_ERI_D2y_Py_Px_Px_b+ABY*I_ERI_Py_Py_Px_Px_b;
  Double I_ERI_Pz_D2y_Px_Px_b = I_ERI_Dyz_Py_Px_Px_b+ABY*I_ERI_Pz_Py_Px_Px_b;
  Double I_ERI_Px_Dyz_Px_Px_b = I_ERI_Dxz_Py_Px_Px_b+ABZ*I_ERI_Px_Py_Px_Px_b;
  Double I_ERI_Py_Dyz_Px_Px_b = I_ERI_Dyz_Py_Px_Px_b+ABZ*I_ERI_Py_Py_Px_Px_b;
  Double I_ERI_Pz_Dyz_Px_Px_b = I_ERI_D2z_Py_Px_Px_b+ABZ*I_ERI_Pz_Py_Px_Px_b;
  Double I_ERI_Px_D2z_Px_Px_b = I_ERI_Dxz_Pz_Px_Px_b+ABZ*I_ERI_Px_Pz_Px_Px_b;
  Double I_ERI_Py_D2z_Px_Px_b = I_ERI_Dyz_Pz_Px_Px_b+ABZ*I_ERI_Py_Pz_Px_Px_b;
  Double I_ERI_Pz_D2z_Px_Px_b = I_ERI_D2z_Pz_Px_Px_b+ABZ*I_ERI_Pz_Pz_Px_Px_b;
  Double I_ERI_Px_D2x_Py_Px_b = I_ERI_D2x_Px_Py_Px_b+ABX*I_ERI_Px_Px_Py_Px_b;
  Double I_ERI_Py_D2x_Py_Px_b = I_ERI_Dxy_Px_Py_Px_b+ABX*I_ERI_Py_Px_Py_Px_b;
  Double I_ERI_Pz_D2x_Py_Px_b = I_ERI_Dxz_Px_Py_Px_b+ABX*I_ERI_Pz_Px_Py_Px_b;
  Double I_ERI_Px_Dxy_Py_Px_b = I_ERI_Dxy_Px_Py_Px_b+ABY*I_ERI_Px_Px_Py_Px_b;
  Double I_ERI_Py_Dxy_Py_Px_b = I_ERI_D2y_Px_Py_Px_b+ABY*I_ERI_Py_Px_Py_Px_b;
  Double I_ERI_Pz_Dxy_Py_Px_b = I_ERI_Dyz_Px_Py_Px_b+ABY*I_ERI_Pz_Px_Py_Px_b;
  Double I_ERI_Px_Dxz_Py_Px_b = I_ERI_Dxz_Px_Py_Px_b+ABZ*I_ERI_Px_Px_Py_Px_b;
  Double I_ERI_Py_Dxz_Py_Px_b = I_ERI_Dyz_Px_Py_Px_b+ABZ*I_ERI_Py_Px_Py_Px_b;
  Double I_ERI_Pz_Dxz_Py_Px_b = I_ERI_D2z_Px_Py_Px_b+ABZ*I_ERI_Pz_Px_Py_Px_b;
  Double I_ERI_Px_D2y_Py_Px_b = I_ERI_Dxy_Py_Py_Px_b+ABY*I_ERI_Px_Py_Py_Px_b;
  Double I_ERI_Py_D2y_Py_Px_b = I_ERI_D2y_Py_Py_Px_b+ABY*I_ERI_Py_Py_Py_Px_b;
  Double I_ERI_Pz_D2y_Py_Px_b = I_ERI_Dyz_Py_Py_Px_b+ABY*I_ERI_Pz_Py_Py_Px_b;
  Double I_ERI_Px_Dyz_Py_Px_b = I_ERI_Dxz_Py_Py_Px_b+ABZ*I_ERI_Px_Py_Py_Px_b;
  Double I_ERI_Py_Dyz_Py_Px_b = I_ERI_Dyz_Py_Py_Px_b+ABZ*I_ERI_Py_Py_Py_Px_b;
  Double I_ERI_Pz_Dyz_Py_Px_b = I_ERI_D2z_Py_Py_Px_b+ABZ*I_ERI_Pz_Py_Py_Px_b;
  Double I_ERI_Px_D2z_Py_Px_b = I_ERI_Dxz_Pz_Py_Px_b+ABZ*I_ERI_Px_Pz_Py_Px_b;
  Double I_ERI_Py_D2z_Py_Px_b = I_ERI_Dyz_Pz_Py_Px_b+ABZ*I_ERI_Py_Pz_Py_Px_b;
  Double I_ERI_Pz_D2z_Py_Px_b = I_ERI_D2z_Pz_Py_Px_b+ABZ*I_ERI_Pz_Pz_Py_Px_b;
  Double I_ERI_Px_D2x_Pz_Px_b = I_ERI_D2x_Px_Pz_Px_b+ABX*I_ERI_Px_Px_Pz_Px_b;
  Double I_ERI_Py_D2x_Pz_Px_b = I_ERI_Dxy_Px_Pz_Px_b+ABX*I_ERI_Py_Px_Pz_Px_b;
  Double I_ERI_Pz_D2x_Pz_Px_b = I_ERI_Dxz_Px_Pz_Px_b+ABX*I_ERI_Pz_Px_Pz_Px_b;
  Double I_ERI_Px_Dxy_Pz_Px_b = I_ERI_Dxy_Px_Pz_Px_b+ABY*I_ERI_Px_Px_Pz_Px_b;
  Double I_ERI_Py_Dxy_Pz_Px_b = I_ERI_D2y_Px_Pz_Px_b+ABY*I_ERI_Py_Px_Pz_Px_b;
  Double I_ERI_Pz_Dxy_Pz_Px_b = I_ERI_Dyz_Px_Pz_Px_b+ABY*I_ERI_Pz_Px_Pz_Px_b;
  Double I_ERI_Px_Dxz_Pz_Px_b = I_ERI_Dxz_Px_Pz_Px_b+ABZ*I_ERI_Px_Px_Pz_Px_b;
  Double I_ERI_Py_Dxz_Pz_Px_b = I_ERI_Dyz_Px_Pz_Px_b+ABZ*I_ERI_Py_Px_Pz_Px_b;
  Double I_ERI_Pz_Dxz_Pz_Px_b = I_ERI_D2z_Px_Pz_Px_b+ABZ*I_ERI_Pz_Px_Pz_Px_b;
  Double I_ERI_Px_D2y_Pz_Px_b = I_ERI_Dxy_Py_Pz_Px_b+ABY*I_ERI_Px_Py_Pz_Px_b;
  Double I_ERI_Py_D2y_Pz_Px_b = I_ERI_D2y_Py_Pz_Px_b+ABY*I_ERI_Py_Py_Pz_Px_b;
  Double I_ERI_Pz_D2y_Pz_Px_b = I_ERI_Dyz_Py_Pz_Px_b+ABY*I_ERI_Pz_Py_Pz_Px_b;
  Double I_ERI_Px_Dyz_Pz_Px_b = I_ERI_Dxz_Py_Pz_Px_b+ABZ*I_ERI_Px_Py_Pz_Px_b;
  Double I_ERI_Py_Dyz_Pz_Px_b = I_ERI_Dyz_Py_Pz_Px_b+ABZ*I_ERI_Py_Py_Pz_Px_b;
  Double I_ERI_Pz_Dyz_Pz_Px_b = I_ERI_D2z_Py_Pz_Px_b+ABZ*I_ERI_Pz_Py_Pz_Px_b;
  Double I_ERI_Px_D2z_Pz_Px_b = I_ERI_Dxz_Pz_Pz_Px_b+ABZ*I_ERI_Px_Pz_Pz_Px_b;
  Double I_ERI_Py_D2z_Pz_Px_b = I_ERI_Dyz_Pz_Pz_Px_b+ABZ*I_ERI_Py_Pz_Pz_Px_b;
  Double I_ERI_Pz_D2z_Pz_Px_b = I_ERI_D2z_Pz_Pz_Px_b+ABZ*I_ERI_Pz_Pz_Pz_Px_b;
  Double I_ERI_Px_D2x_Px_Py_b = I_ERI_D2x_Px_Px_Py_b+ABX*I_ERI_Px_Px_Px_Py_b;
  Double I_ERI_Py_D2x_Px_Py_b = I_ERI_Dxy_Px_Px_Py_b+ABX*I_ERI_Py_Px_Px_Py_b;
  Double I_ERI_Pz_D2x_Px_Py_b = I_ERI_Dxz_Px_Px_Py_b+ABX*I_ERI_Pz_Px_Px_Py_b;
  Double I_ERI_Px_Dxy_Px_Py_b = I_ERI_Dxy_Px_Px_Py_b+ABY*I_ERI_Px_Px_Px_Py_b;
  Double I_ERI_Py_Dxy_Px_Py_b = I_ERI_D2y_Px_Px_Py_b+ABY*I_ERI_Py_Px_Px_Py_b;
  Double I_ERI_Pz_Dxy_Px_Py_b = I_ERI_Dyz_Px_Px_Py_b+ABY*I_ERI_Pz_Px_Px_Py_b;
  Double I_ERI_Px_Dxz_Px_Py_b = I_ERI_Dxz_Px_Px_Py_b+ABZ*I_ERI_Px_Px_Px_Py_b;
  Double I_ERI_Py_Dxz_Px_Py_b = I_ERI_Dyz_Px_Px_Py_b+ABZ*I_ERI_Py_Px_Px_Py_b;
  Double I_ERI_Pz_Dxz_Px_Py_b = I_ERI_D2z_Px_Px_Py_b+ABZ*I_ERI_Pz_Px_Px_Py_b;
  Double I_ERI_Px_D2y_Px_Py_b = I_ERI_Dxy_Py_Px_Py_b+ABY*I_ERI_Px_Py_Px_Py_b;
  Double I_ERI_Py_D2y_Px_Py_b = I_ERI_D2y_Py_Px_Py_b+ABY*I_ERI_Py_Py_Px_Py_b;
  Double I_ERI_Pz_D2y_Px_Py_b = I_ERI_Dyz_Py_Px_Py_b+ABY*I_ERI_Pz_Py_Px_Py_b;
  Double I_ERI_Px_Dyz_Px_Py_b = I_ERI_Dxz_Py_Px_Py_b+ABZ*I_ERI_Px_Py_Px_Py_b;
  Double I_ERI_Py_Dyz_Px_Py_b = I_ERI_Dyz_Py_Px_Py_b+ABZ*I_ERI_Py_Py_Px_Py_b;
  Double I_ERI_Pz_Dyz_Px_Py_b = I_ERI_D2z_Py_Px_Py_b+ABZ*I_ERI_Pz_Py_Px_Py_b;
  Double I_ERI_Px_D2z_Px_Py_b = I_ERI_Dxz_Pz_Px_Py_b+ABZ*I_ERI_Px_Pz_Px_Py_b;
  Double I_ERI_Py_D2z_Px_Py_b = I_ERI_Dyz_Pz_Px_Py_b+ABZ*I_ERI_Py_Pz_Px_Py_b;
  Double I_ERI_Pz_D2z_Px_Py_b = I_ERI_D2z_Pz_Px_Py_b+ABZ*I_ERI_Pz_Pz_Px_Py_b;
  Double I_ERI_Px_D2x_Py_Py_b = I_ERI_D2x_Px_Py_Py_b+ABX*I_ERI_Px_Px_Py_Py_b;
  Double I_ERI_Py_D2x_Py_Py_b = I_ERI_Dxy_Px_Py_Py_b+ABX*I_ERI_Py_Px_Py_Py_b;
  Double I_ERI_Pz_D2x_Py_Py_b = I_ERI_Dxz_Px_Py_Py_b+ABX*I_ERI_Pz_Px_Py_Py_b;
  Double I_ERI_Px_Dxy_Py_Py_b = I_ERI_Dxy_Px_Py_Py_b+ABY*I_ERI_Px_Px_Py_Py_b;
  Double I_ERI_Py_Dxy_Py_Py_b = I_ERI_D2y_Px_Py_Py_b+ABY*I_ERI_Py_Px_Py_Py_b;
  Double I_ERI_Pz_Dxy_Py_Py_b = I_ERI_Dyz_Px_Py_Py_b+ABY*I_ERI_Pz_Px_Py_Py_b;
  Double I_ERI_Px_Dxz_Py_Py_b = I_ERI_Dxz_Px_Py_Py_b+ABZ*I_ERI_Px_Px_Py_Py_b;
  Double I_ERI_Py_Dxz_Py_Py_b = I_ERI_Dyz_Px_Py_Py_b+ABZ*I_ERI_Py_Px_Py_Py_b;
  Double I_ERI_Pz_Dxz_Py_Py_b = I_ERI_D2z_Px_Py_Py_b+ABZ*I_ERI_Pz_Px_Py_Py_b;
  Double I_ERI_Px_D2y_Py_Py_b = I_ERI_Dxy_Py_Py_Py_b+ABY*I_ERI_Px_Py_Py_Py_b;
  Double I_ERI_Py_D2y_Py_Py_b = I_ERI_D2y_Py_Py_Py_b+ABY*I_ERI_Py_Py_Py_Py_b;
  Double I_ERI_Pz_D2y_Py_Py_b = I_ERI_Dyz_Py_Py_Py_b+ABY*I_ERI_Pz_Py_Py_Py_b;
  Double I_ERI_Px_Dyz_Py_Py_b = I_ERI_Dxz_Py_Py_Py_b+ABZ*I_ERI_Px_Py_Py_Py_b;
  Double I_ERI_Py_Dyz_Py_Py_b = I_ERI_Dyz_Py_Py_Py_b+ABZ*I_ERI_Py_Py_Py_Py_b;
  Double I_ERI_Pz_Dyz_Py_Py_b = I_ERI_D2z_Py_Py_Py_b+ABZ*I_ERI_Pz_Py_Py_Py_b;
  Double I_ERI_Px_D2z_Py_Py_b = I_ERI_Dxz_Pz_Py_Py_b+ABZ*I_ERI_Px_Pz_Py_Py_b;
  Double I_ERI_Py_D2z_Py_Py_b = I_ERI_Dyz_Pz_Py_Py_b+ABZ*I_ERI_Py_Pz_Py_Py_b;
  Double I_ERI_Pz_D2z_Py_Py_b = I_ERI_D2z_Pz_Py_Py_b+ABZ*I_ERI_Pz_Pz_Py_Py_b;
  Double I_ERI_Px_D2x_Pz_Py_b = I_ERI_D2x_Px_Pz_Py_b+ABX*I_ERI_Px_Px_Pz_Py_b;
  Double I_ERI_Py_D2x_Pz_Py_b = I_ERI_Dxy_Px_Pz_Py_b+ABX*I_ERI_Py_Px_Pz_Py_b;
  Double I_ERI_Pz_D2x_Pz_Py_b = I_ERI_Dxz_Px_Pz_Py_b+ABX*I_ERI_Pz_Px_Pz_Py_b;
  Double I_ERI_Px_Dxy_Pz_Py_b = I_ERI_Dxy_Px_Pz_Py_b+ABY*I_ERI_Px_Px_Pz_Py_b;
  Double I_ERI_Py_Dxy_Pz_Py_b = I_ERI_D2y_Px_Pz_Py_b+ABY*I_ERI_Py_Px_Pz_Py_b;
  Double I_ERI_Pz_Dxy_Pz_Py_b = I_ERI_Dyz_Px_Pz_Py_b+ABY*I_ERI_Pz_Px_Pz_Py_b;
  Double I_ERI_Px_Dxz_Pz_Py_b = I_ERI_Dxz_Px_Pz_Py_b+ABZ*I_ERI_Px_Px_Pz_Py_b;
  Double I_ERI_Py_Dxz_Pz_Py_b = I_ERI_Dyz_Px_Pz_Py_b+ABZ*I_ERI_Py_Px_Pz_Py_b;
  Double I_ERI_Pz_Dxz_Pz_Py_b = I_ERI_D2z_Px_Pz_Py_b+ABZ*I_ERI_Pz_Px_Pz_Py_b;
  Double I_ERI_Px_D2y_Pz_Py_b = I_ERI_Dxy_Py_Pz_Py_b+ABY*I_ERI_Px_Py_Pz_Py_b;
  Double I_ERI_Py_D2y_Pz_Py_b = I_ERI_D2y_Py_Pz_Py_b+ABY*I_ERI_Py_Py_Pz_Py_b;
  Double I_ERI_Pz_D2y_Pz_Py_b = I_ERI_Dyz_Py_Pz_Py_b+ABY*I_ERI_Pz_Py_Pz_Py_b;
  Double I_ERI_Px_Dyz_Pz_Py_b = I_ERI_Dxz_Py_Pz_Py_b+ABZ*I_ERI_Px_Py_Pz_Py_b;
  Double I_ERI_Py_Dyz_Pz_Py_b = I_ERI_Dyz_Py_Pz_Py_b+ABZ*I_ERI_Py_Py_Pz_Py_b;
  Double I_ERI_Pz_Dyz_Pz_Py_b = I_ERI_D2z_Py_Pz_Py_b+ABZ*I_ERI_Pz_Py_Pz_Py_b;
  Double I_ERI_Px_D2z_Pz_Py_b = I_ERI_Dxz_Pz_Pz_Py_b+ABZ*I_ERI_Px_Pz_Pz_Py_b;
  Double I_ERI_Py_D2z_Pz_Py_b = I_ERI_Dyz_Pz_Pz_Py_b+ABZ*I_ERI_Py_Pz_Pz_Py_b;
  Double I_ERI_Pz_D2z_Pz_Py_b = I_ERI_D2z_Pz_Pz_Py_b+ABZ*I_ERI_Pz_Pz_Pz_Py_b;
  Double I_ERI_Px_D2x_Px_Pz_b = I_ERI_D2x_Px_Px_Pz_b+ABX*I_ERI_Px_Px_Px_Pz_b;
  Double I_ERI_Py_D2x_Px_Pz_b = I_ERI_Dxy_Px_Px_Pz_b+ABX*I_ERI_Py_Px_Px_Pz_b;
  Double I_ERI_Pz_D2x_Px_Pz_b = I_ERI_Dxz_Px_Px_Pz_b+ABX*I_ERI_Pz_Px_Px_Pz_b;
  Double I_ERI_Px_Dxy_Px_Pz_b = I_ERI_Dxy_Px_Px_Pz_b+ABY*I_ERI_Px_Px_Px_Pz_b;
  Double I_ERI_Py_Dxy_Px_Pz_b = I_ERI_D2y_Px_Px_Pz_b+ABY*I_ERI_Py_Px_Px_Pz_b;
  Double I_ERI_Pz_Dxy_Px_Pz_b = I_ERI_Dyz_Px_Px_Pz_b+ABY*I_ERI_Pz_Px_Px_Pz_b;
  Double I_ERI_Px_Dxz_Px_Pz_b = I_ERI_Dxz_Px_Px_Pz_b+ABZ*I_ERI_Px_Px_Px_Pz_b;
  Double I_ERI_Py_Dxz_Px_Pz_b = I_ERI_Dyz_Px_Px_Pz_b+ABZ*I_ERI_Py_Px_Px_Pz_b;
  Double I_ERI_Pz_Dxz_Px_Pz_b = I_ERI_D2z_Px_Px_Pz_b+ABZ*I_ERI_Pz_Px_Px_Pz_b;
  Double I_ERI_Px_D2y_Px_Pz_b = I_ERI_Dxy_Py_Px_Pz_b+ABY*I_ERI_Px_Py_Px_Pz_b;
  Double I_ERI_Py_D2y_Px_Pz_b = I_ERI_D2y_Py_Px_Pz_b+ABY*I_ERI_Py_Py_Px_Pz_b;
  Double I_ERI_Pz_D2y_Px_Pz_b = I_ERI_Dyz_Py_Px_Pz_b+ABY*I_ERI_Pz_Py_Px_Pz_b;
  Double I_ERI_Px_Dyz_Px_Pz_b = I_ERI_Dxz_Py_Px_Pz_b+ABZ*I_ERI_Px_Py_Px_Pz_b;
  Double I_ERI_Py_Dyz_Px_Pz_b = I_ERI_Dyz_Py_Px_Pz_b+ABZ*I_ERI_Py_Py_Px_Pz_b;
  Double I_ERI_Pz_Dyz_Px_Pz_b = I_ERI_D2z_Py_Px_Pz_b+ABZ*I_ERI_Pz_Py_Px_Pz_b;
  Double I_ERI_Px_D2z_Px_Pz_b = I_ERI_Dxz_Pz_Px_Pz_b+ABZ*I_ERI_Px_Pz_Px_Pz_b;
  Double I_ERI_Py_D2z_Px_Pz_b = I_ERI_Dyz_Pz_Px_Pz_b+ABZ*I_ERI_Py_Pz_Px_Pz_b;
  Double I_ERI_Pz_D2z_Px_Pz_b = I_ERI_D2z_Pz_Px_Pz_b+ABZ*I_ERI_Pz_Pz_Px_Pz_b;
  Double I_ERI_Px_D2x_Py_Pz_b = I_ERI_D2x_Px_Py_Pz_b+ABX*I_ERI_Px_Px_Py_Pz_b;
  Double I_ERI_Py_D2x_Py_Pz_b = I_ERI_Dxy_Px_Py_Pz_b+ABX*I_ERI_Py_Px_Py_Pz_b;
  Double I_ERI_Pz_D2x_Py_Pz_b = I_ERI_Dxz_Px_Py_Pz_b+ABX*I_ERI_Pz_Px_Py_Pz_b;
  Double I_ERI_Px_Dxy_Py_Pz_b = I_ERI_Dxy_Px_Py_Pz_b+ABY*I_ERI_Px_Px_Py_Pz_b;
  Double I_ERI_Py_Dxy_Py_Pz_b = I_ERI_D2y_Px_Py_Pz_b+ABY*I_ERI_Py_Px_Py_Pz_b;
  Double I_ERI_Pz_Dxy_Py_Pz_b = I_ERI_Dyz_Px_Py_Pz_b+ABY*I_ERI_Pz_Px_Py_Pz_b;
  Double I_ERI_Px_Dxz_Py_Pz_b = I_ERI_Dxz_Px_Py_Pz_b+ABZ*I_ERI_Px_Px_Py_Pz_b;
  Double I_ERI_Py_Dxz_Py_Pz_b = I_ERI_Dyz_Px_Py_Pz_b+ABZ*I_ERI_Py_Px_Py_Pz_b;
  Double I_ERI_Pz_Dxz_Py_Pz_b = I_ERI_D2z_Px_Py_Pz_b+ABZ*I_ERI_Pz_Px_Py_Pz_b;
  Double I_ERI_Px_D2y_Py_Pz_b = I_ERI_Dxy_Py_Py_Pz_b+ABY*I_ERI_Px_Py_Py_Pz_b;
  Double I_ERI_Py_D2y_Py_Pz_b = I_ERI_D2y_Py_Py_Pz_b+ABY*I_ERI_Py_Py_Py_Pz_b;
  Double I_ERI_Pz_D2y_Py_Pz_b = I_ERI_Dyz_Py_Py_Pz_b+ABY*I_ERI_Pz_Py_Py_Pz_b;
  Double I_ERI_Px_Dyz_Py_Pz_b = I_ERI_Dxz_Py_Py_Pz_b+ABZ*I_ERI_Px_Py_Py_Pz_b;
  Double I_ERI_Py_Dyz_Py_Pz_b = I_ERI_Dyz_Py_Py_Pz_b+ABZ*I_ERI_Py_Py_Py_Pz_b;
  Double I_ERI_Pz_Dyz_Py_Pz_b = I_ERI_D2z_Py_Py_Pz_b+ABZ*I_ERI_Pz_Py_Py_Pz_b;
  Double I_ERI_Px_D2z_Py_Pz_b = I_ERI_Dxz_Pz_Py_Pz_b+ABZ*I_ERI_Px_Pz_Py_Pz_b;
  Double I_ERI_Py_D2z_Py_Pz_b = I_ERI_Dyz_Pz_Py_Pz_b+ABZ*I_ERI_Py_Pz_Py_Pz_b;
  Double I_ERI_Pz_D2z_Py_Pz_b = I_ERI_D2z_Pz_Py_Pz_b+ABZ*I_ERI_Pz_Pz_Py_Pz_b;
  Double I_ERI_Px_D2x_Pz_Pz_b = I_ERI_D2x_Px_Pz_Pz_b+ABX*I_ERI_Px_Px_Pz_Pz_b;
  Double I_ERI_Py_D2x_Pz_Pz_b = I_ERI_Dxy_Px_Pz_Pz_b+ABX*I_ERI_Py_Px_Pz_Pz_b;
  Double I_ERI_Pz_D2x_Pz_Pz_b = I_ERI_Dxz_Px_Pz_Pz_b+ABX*I_ERI_Pz_Px_Pz_Pz_b;
  Double I_ERI_Px_Dxy_Pz_Pz_b = I_ERI_Dxy_Px_Pz_Pz_b+ABY*I_ERI_Px_Px_Pz_Pz_b;
  Double I_ERI_Py_Dxy_Pz_Pz_b = I_ERI_D2y_Px_Pz_Pz_b+ABY*I_ERI_Py_Px_Pz_Pz_b;
  Double I_ERI_Pz_Dxy_Pz_Pz_b = I_ERI_Dyz_Px_Pz_Pz_b+ABY*I_ERI_Pz_Px_Pz_Pz_b;
  Double I_ERI_Px_Dxz_Pz_Pz_b = I_ERI_Dxz_Px_Pz_Pz_b+ABZ*I_ERI_Px_Px_Pz_Pz_b;
  Double I_ERI_Py_Dxz_Pz_Pz_b = I_ERI_Dyz_Px_Pz_Pz_b+ABZ*I_ERI_Py_Px_Pz_Pz_b;
  Double I_ERI_Pz_Dxz_Pz_Pz_b = I_ERI_D2z_Px_Pz_Pz_b+ABZ*I_ERI_Pz_Px_Pz_Pz_b;
  Double I_ERI_Px_D2y_Pz_Pz_b = I_ERI_Dxy_Py_Pz_Pz_b+ABY*I_ERI_Px_Py_Pz_Pz_b;
  Double I_ERI_Py_D2y_Pz_Pz_b = I_ERI_D2y_Py_Pz_Pz_b+ABY*I_ERI_Py_Py_Pz_Pz_b;
  Double I_ERI_Pz_D2y_Pz_Pz_b = I_ERI_Dyz_Py_Pz_Pz_b+ABY*I_ERI_Pz_Py_Pz_Pz_b;
  Double I_ERI_Px_Dyz_Pz_Pz_b = I_ERI_Dxz_Py_Pz_Pz_b+ABZ*I_ERI_Px_Py_Pz_Pz_b;
  Double I_ERI_Py_Dyz_Pz_Pz_b = I_ERI_Dyz_Py_Pz_Pz_b+ABZ*I_ERI_Py_Py_Pz_Pz_b;
  Double I_ERI_Pz_Dyz_Pz_Pz_b = I_ERI_D2z_Py_Pz_Pz_b+ABZ*I_ERI_Pz_Py_Pz_Pz_b;
  Double I_ERI_Px_D2z_Pz_Pz_b = I_ERI_Dxz_Pz_Pz_Pz_b+ABZ*I_ERI_Px_Pz_Pz_Pz_b;
  Double I_ERI_Py_D2z_Pz_Pz_b = I_ERI_Dyz_Pz_Pz_Pz_b+ABZ*I_ERI_Py_Pz_Pz_Pz_b;
  Double I_ERI_Pz_D2z_Pz_Pz_b = I_ERI_D2z_Pz_Pz_Pz_b+ABZ*I_ERI_Pz_Pz_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_P_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_P_c
   * RHS shell quartet name: SQ_ERI_P_S_D_P_c
   ************************************************************/
  Double I_ERI_Px_Px_D2x_Px_c = I_ERI_D2x_S_D2x_Px_c+ABX*I_ERI_Px_S_D2x_Px_c;
  Double I_ERI_Py_Px_D2x_Px_c = I_ERI_Dxy_S_D2x_Px_c+ABX*I_ERI_Py_S_D2x_Px_c;
  Double I_ERI_Pz_Px_D2x_Px_c = I_ERI_Dxz_S_D2x_Px_c+ABX*I_ERI_Pz_S_D2x_Px_c;
  Double I_ERI_Px_Py_D2x_Px_c = I_ERI_Dxy_S_D2x_Px_c+ABY*I_ERI_Px_S_D2x_Px_c;
  Double I_ERI_Py_Py_D2x_Px_c = I_ERI_D2y_S_D2x_Px_c+ABY*I_ERI_Py_S_D2x_Px_c;
  Double I_ERI_Pz_Py_D2x_Px_c = I_ERI_Dyz_S_D2x_Px_c+ABY*I_ERI_Pz_S_D2x_Px_c;
  Double I_ERI_Px_Pz_D2x_Px_c = I_ERI_Dxz_S_D2x_Px_c+ABZ*I_ERI_Px_S_D2x_Px_c;
  Double I_ERI_Py_Pz_D2x_Px_c = I_ERI_Dyz_S_D2x_Px_c+ABZ*I_ERI_Py_S_D2x_Px_c;
  Double I_ERI_Pz_Pz_D2x_Px_c = I_ERI_D2z_S_D2x_Px_c+ABZ*I_ERI_Pz_S_D2x_Px_c;
  Double I_ERI_Px_Px_Dxy_Px_c = I_ERI_D2x_S_Dxy_Px_c+ABX*I_ERI_Px_S_Dxy_Px_c;
  Double I_ERI_Py_Px_Dxy_Px_c = I_ERI_Dxy_S_Dxy_Px_c+ABX*I_ERI_Py_S_Dxy_Px_c;
  Double I_ERI_Pz_Px_Dxy_Px_c = I_ERI_Dxz_S_Dxy_Px_c+ABX*I_ERI_Pz_S_Dxy_Px_c;
  Double I_ERI_Px_Py_Dxy_Px_c = I_ERI_Dxy_S_Dxy_Px_c+ABY*I_ERI_Px_S_Dxy_Px_c;
  Double I_ERI_Py_Py_Dxy_Px_c = I_ERI_D2y_S_Dxy_Px_c+ABY*I_ERI_Py_S_Dxy_Px_c;
  Double I_ERI_Pz_Py_Dxy_Px_c = I_ERI_Dyz_S_Dxy_Px_c+ABY*I_ERI_Pz_S_Dxy_Px_c;
  Double I_ERI_Px_Pz_Dxy_Px_c = I_ERI_Dxz_S_Dxy_Px_c+ABZ*I_ERI_Px_S_Dxy_Px_c;
  Double I_ERI_Py_Pz_Dxy_Px_c = I_ERI_Dyz_S_Dxy_Px_c+ABZ*I_ERI_Py_S_Dxy_Px_c;
  Double I_ERI_Pz_Pz_Dxy_Px_c = I_ERI_D2z_S_Dxy_Px_c+ABZ*I_ERI_Pz_S_Dxy_Px_c;
  Double I_ERI_Px_Px_Dxz_Px_c = I_ERI_D2x_S_Dxz_Px_c+ABX*I_ERI_Px_S_Dxz_Px_c;
  Double I_ERI_Py_Px_Dxz_Px_c = I_ERI_Dxy_S_Dxz_Px_c+ABX*I_ERI_Py_S_Dxz_Px_c;
  Double I_ERI_Pz_Px_Dxz_Px_c = I_ERI_Dxz_S_Dxz_Px_c+ABX*I_ERI_Pz_S_Dxz_Px_c;
  Double I_ERI_Px_Py_Dxz_Px_c = I_ERI_Dxy_S_Dxz_Px_c+ABY*I_ERI_Px_S_Dxz_Px_c;
  Double I_ERI_Py_Py_Dxz_Px_c = I_ERI_D2y_S_Dxz_Px_c+ABY*I_ERI_Py_S_Dxz_Px_c;
  Double I_ERI_Pz_Py_Dxz_Px_c = I_ERI_Dyz_S_Dxz_Px_c+ABY*I_ERI_Pz_S_Dxz_Px_c;
  Double I_ERI_Px_Pz_Dxz_Px_c = I_ERI_Dxz_S_Dxz_Px_c+ABZ*I_ERI_Px_S_Dxz_Px_c;
  Double I_ERI_Py_Pz_Dxz_Px_c = I_ERI_Dyz_S_Dxz_Px_c+ABZ*I_ERI_Py_S_Dxz_Px_c;
  Double I_ERI_Pz_Pz_Dxz_Px_c = I_ERI_D2z_S_Dxz_Px_c+ABZ*I_ERI_Pz_S_Dxz_Px_c;
  Double I_ERI_Px_Px_D2y_Px_c = I_ERI_D2x_S_D2y_Px_c+ABX*I_ERI_Px_S_D2y_Px_c;
  Double I_ERI_Py_Px_D2y_Px_c = I_ERI_Dxy_S_D2y_Px_c+ABX*I_ERI_Py_S_D2y_Px_c;
  Double I_ERI_Pz_Px_D2y_Px_c = I_ERI_Dxz_S_D2y_Px_c+ABX*I_ERI_Pz_S_D2y_Px_c;
  Double I_ERI_Px_Py_D2y_Px_c = I_ERI_Dxy_S_D2y_Px_c+ABY*I_ERI_Px_S_D2y_Px_c;
  Double I_ERI_Py_Py_D2y_Px_c = I_ERI_D2y_S_D2y_Px_c+ABY*I_ERI_Py_S_D2y_Px_c;
  Double I_ERI_Pz_Py_D2y_Px_c = I_ERI_Dyz_S_D2y_Px_c+ABY*I_ERI_Pz_S_D2y_Px_c;
  Double I_ERI_Px_Pz_D2y_Px_c = I_ERI_Dxz_S_D2y_Px_c+ABZ*I_ERI_Px_S_D2y_Px_c;
  Double I_ERI_Py_Pz_D2y_Px_c = I_ERI_Dyz_S_D2y_Px_c+ABZ*I_ERI_Py_S_D2y_Px_c;
  Double I_ERI_Pz_Pz_D2y_Px_c = I_ERI_D2z_S_D2y_Px_c+ABZ*I_ERI_Pz_S_D2y_Px_c;
  Double I_ERI_Px_Px_Dyz_Px_c = I_ERI_D2x_S_Dyz_Px_c+ABX*I_ERI_Px_S_Dyz_Px_c;
  Double I_ERI_Py_Px_Dyz_Px_c = I_ERI_Dxy_S_Dyz_Px_c+ABX*I_ERI_Py_S_Dyz_Px_c;
  Double I_ERI_Pz_Px_Dyz_Px_c = I_ERI_Dxz_S_Dyz_Px_c+ABX*I_ERI_Pz_S_Dyz_Px_c;
  Double I_ERI_Px_Py_Dyz_Px_c = I_ERI_Dxy_S_Dyz_Px_c+ABY*I_ERI_Px_S_Dyz_Px_c;
  Double I_ERI_Py_Py_Dyz_Px_c = I_ERI_D2y_S_Dyz_Px_c+ABY*I_ERI_Py_S_Dyz_Px_c;
  Double I_ERI_Pz_Py_Dyz_Px_c = I_ERI_Dyz_S_Dyz_Px_c+ABY*I_ERI_Pz_S_Dyz_Px_c;
  Double I_ERI_Px_Pz_Dyz_Px_c = I_ERI_Dxz_S_Dyz_Px_c+ABZ*I_ERI_Px_S_Dyz_Px_c;
  Double I_ERI_Py_Pz_Dyz_Px_c = I_ERI_Dyz_S_Dyz_Px_c+ABZ*I_ERI_Py_S_Dyz_Px_c;
  Double I_ERI_Pz_Pz_Dyz_Px_c = I_ERI_D2z_S_Dyz_Px_c+ABZ*I_ERI_Pz_S_Dyz_Px_c;
  Double I_ERI_Px_Px_D2z_Px_c = I_ERI_D2x_S_D2z_Px_c+ABX*I_ERI_Px_S_D2z_Px_c;
  Double I_ERI_Py_Px_D2z_Px_c = I_ERI_Dxy_S_D2z_Px_c+ABX*I_ERI_Py_S_D2z_Px_c;
  Double I_ERI_Pz_Px_D2z_Px_c = I_ERI_Dxz_S_D2z_Px_c+ABX*I_ERI_Pz_S_D2z_Px_c;
  Double I_ERI_Px_Py_D2z_Px_c = I_ERI_Dxy_S_D2z_Px_c+ABY*I_ERI_Px_S_D2z_Px_c;
  Double I_ERI_Py_Py_D2z_Px_c = I_ERI_D2y_S_D2z_Px_c+ABY*I_ERI_Py_S_D2z_Px_c;
  Double I_ERI_Pz_Py_D2z_Px_c = I_ERI_Dyz_S_D2z_Px_c+ABY*I_ERI_Pz_S_D2z_Px_c;
  Double I_ERI_Px_Pz_D2z_Px_c = I_ERI_Dxz_S_D2z_Px_c+ABZ*I_ERI_Px_S_D2z_Px_c;
  Double I_ERI_Py_Pz_D2z_Px_c = I_ERI_Dyz_S_D2z_Px_c+ABZ*I_ERI_Py_S_D2z_Px_c;
  Double I_ERI_Pz_Pz_D2z_Px_c = I_ERI_D2z_S_D2z_Px_c+ABZ*I_ERI_Pz_S_D2z_Px_c;
  Double I_ERI_Px_Px_D2x_Py_c = I_ERI_D2x_S_D2x_Py_c+ABX*I_ERI_Px_S_D2x_Py_c;
  Double I_ERI_Py_Px_D2x_Py_c = I_ERI_Dxy_S_D2x_Py_c+ABX*I_ERI_Py_S_D2x_Py_c;
  Double I_ERI_Pz_Px_D2x_Py_c = I_ERI_Dxz_S_D2x_Py_c+ABX*I_ERI_Pz_S_D2x_Py_c;
  Double I_ERI_Px_Py_D2x_Py_c = I_ERI_Dxy_S_D2x_Py_c+ABY*I_ERI_Px_S_D2x_Py_c;
  Double I_ERI_Py_Py_D2x_Py_c = I_ERI_D2y_S_D2x_Py_c+ABY*I_ERI_Py_S_D2x_Py_c;
  Double I_ERI_Pz_Py_D2x_Py_c = I_ERI_Dyz_S_D2x_Py_c+ABY*I_ERI_Pz_S_D2x_Py_c;
  Double I_ERI_Px_Pz_D2x_Py_c = I_ERI_Dxz_S_D2x_Py_c+ABZ*I_ERI_Px_S_D2x_Py_c;
  Double I_ERI_Py_Pz_D2x_Py_c = I_ERI_Dyz_S_D2x_Py_c+ABZ*I_ERI_Py_S_D2x_Py_c;
  Double I_ERI_Pz_Pz_D2x_Py_c = I_ERI_D2z_S_D2x_Py_c+ABZ*I_ERI_Pz_S_D2x_Py_c;
  Double I_ERI_Px_Px_Dxy_Py_c = I_ERI_D2x_S_Dxy_Py_c+ABX*I_ERI_Px_S_Dxy_Py_c;
  Double I_ERI_Py_Px_Dxy_Py_c = I_ERI_Dxy_S_Dxy_Py_c+ABX*I_ERI_Py_S_Dxy_Py_c;
  Double I_ERI_Pz_Px_Dxy_Py_c = I_ERI_Dxz_S_Dxy_Py_c+ABX*I_ERI_Pz_S_Dxy_Py_c;
  Double I_ERI_Px_Py_Dxy_Py_c = I_ERI_Dxy_S_Dxy_Py_c+ABY*I_ERI_Px_S_Dxy_Py_c;
  Double I_ERI_Py_Py_Dxy_Py_c = I_ERI_D2y_S_Dxy_Py_c+ABY*I_ERI_Py_S_Dxy_Py_c;
  Double I_ERI_Pz_Py_Dxy_Py_c = I_ERI_Dyz_S_Dxy_Py_c+ABY*I_ERI_Pz_S_Dxy_Py_c;
  Double I_ERI_Px_Pz_Dxy_Py_c = I_ERI_Dxz_S_Dxy_Py_c+ABZ*I_ERI_Px_S_Dxy_Py_c;
  Double I_ERI_Py_Pz_Dxy_Py_c = I_ERI_Dyz_S_Dxy_Py_c+ABZ*I_ERI_Py_S_Dxy_Py_c;
  Double I_ERI_Pz_Pz_Dxy_Py_c = I_ERI_D2z_S_Dxy_Py_c+ABZ*I_ERI_Pz_S_Dxy_Py_c;
  Double I_ERI_Px_Px_Dxz_Py_c = I_ERI_D2x_S_Dxz_Py_c+ABX*I_ERI_Px_S_Dxz_Py_c;
  Double I_ERI_Py_Px_Dxz_Py_c = I_ERI_Dxy_S_Dxz_Py_c+ABX*I_ERI_Py_S_Dxz_Py_c;
  Double I_ERI_Pz_Px_Dxz_Py_c = I_ERI_Dxz_S_Dxz_Py_c+ABX*I_ERI_Pz_S_Dxz_Py_c;
  Double I_ERI_Px_Py_Dxz_Py_c = I_ERI_Dxy_S_Dxz_Py_c+ABY*I_ERI_Px_S_Dxz_Py_c;
  Double I_ERI_Py_Py_Dxz_Py_c = I_ERI_D2y_S_Dxz_Py_c+ABY*I_ERI_Py_S_Dxz_Py_c;
  Double I_ERI_Pz_Py_Dxz_Py_c = I_ERI_Dyz_S_Dxz_Py_c+ABY*I_ERI_Pz_S_Dxz_Py_c;
  Double I_ERI_Px_Pz_Dxz_Py_c = I_ERI_Dxz_S_Dxz_Py_c+ABZ*I_ERI_Px_S_Dxz_Py_c;
  Double I_ERI_Py_Pz_Dxz_Py_c = I_ERI_Dyz_S_Dxz_Py_c+ABZ*I_ERI_Py_S_Dxz_Py_c;
  Double I_ERI_Pz_Pz_Dxz_Py_c = I_ERI_D2z_S_Dxz_Py_c+ABZ*I_ERI_Pz_S_Dxz_Py_c;
  Double I_ERI_Px_Px_D2y_Py_c = I_ERI_D2x_S_D2y_Py_c+ABX*I_ERI_Px_S_D2y_Py_c;
  Double I_ERI_Py_Px_D2y_Py_c = I_ERI_Dxy_S_D2y_Py_c+ABX*I_ERI_Py_S_D2y_Py_c;
  Double I_ERI_Pz_Px_D2y_Py_c = I_ERI_Dxz_S_D2y_Py_c+ABX*I_ERI_Pz_S_D2y_Py_c;
  Double I_ERI_Px_Py_D2y_Py_c = I_ERI_Dxy_S_D2y_Py_c+ABY*I_ERI_Px_S_D2y_Py_c;
  Double I_ERI_Py_Py_D2y_Py_c = I_ERI_D2y_S_D2y_Py_c+ABY*I_ERI_Py_S_D2y_Py_c;
  Double I_ERI_Pz_Py_D2y_Py_c = I_ERI_Dyz_S_D2y_Py_c+ABY*I_ERI_Pz_S_D2y_Py_c;
  Double I_ERI_Px_Pz_D2y_Py_c = I_ERI_Dxz_S_D2y_Py_c+ABZ*I_ERI_Px_S_D2y_Py_c;
  Double I_ERI_Py_Pz_D2y_Py_c = I_ERI_Dyz_S_D2y_Py_c+ABZ*I_ERI_Py_S_D2y_Py_c;
  Double I_ERI_Pz_Pz_D2y_Py_c = I_ERI_D2z_S_D2y_Py_c+ABZ*I_ERI_Pz_S_D2y_Py_c;
  Double I_ERI_Px_Px_Dyz_Py_c = I_ERI_D2x_S_Dyz_Py_c+ABX*I_ERI_Px_S_Dyz_Py_c;
  Double I_ERI_Py_Px_Dyz_Py_c = I_ERI_Dxy_S_Dyz_Py_c+ABX*I_ERI_Py_S_Dyz_Py_c;
  Double I_ERI_Pz_Px_Dyz_Py_c = I_ERI_Dxz_S_Dyz_Py_c+ABX*I_ERI_Pz_S_Dyz_Py_c;
  Double I_ERI_Px_Py_Dyz_Py_c = I_ERI_Dxy_S_Dyz_Py_c+ABY*I_ERI_Px_S_Dyz_Py_c;
  Double I_ERI_Py_Py_Dyz_Py_c = I_ERI_D2y_S_Dyz_Py_c+ABY*I_ERI_Py_S_Dyz_Py_c;
  Double I_ERI_Pz_Py_Dyz_Py_c = I_ERI_Dyz_S_Dyz_Py_c+ABY*I_ERI_Pz_S_Dyz_Py_c;
  Double I_ERI_Px_Pz_Dyz_Py_c = I_ERI_Dxz_S_Dyz_Py_c+ABZ*I_ERI_Px_S_Dyz_Py_c;
  Double I_ERI_Py_Pz_Dyz_Py_c = I_ERI_Dyz_S_Dyz_Py_c+ABZ*I_ERI_Py_S_Dyz_Py_c;
  Double I_ERI_Pz_Pz_Dyz_Py_c = I_ERI_D2z_S_Dyz_Py_c+ABZ*I_ERI_Pz_S_Dyz_Py_c;
  Double I_ERI_Px_Px_D2z_Py_c = I_ERI_D2x_S_D2z_Py_c+ABX*I_ERI_Px_S_D2z_Py_c;
  Double I_ERI_Py_Px_D2z_Py_c = I_ERI_Dxy_S_D2z_Py_c+ABX*I_ERI_Py_S_D2z_Py_c;
  Double I_ERI_Pz_Px_D2z_Py_c = I_ERI_Dxz_S_D2z_Py_c+ABX*I_ERI_Pz_S_D2z_Py_c;
  Double I_ERI_Px_Py_D2z_Py_c = I_ERI_Dxy_S_D2z_Py_c+ABY*I_ERI_Px_S_D2z_Py_c;
  Double I_ERI_Py_Py_D2z_Py_c = I_ERI_D2y_S_D2z_Py_c+ABY*I_ERI_Py_S_D2z_Py_c;
  Double I_ERI_Pz_Py_D2z_Py_c = I_ERI_Dyz_S_D2z_Py_c+ABY*I_ERI_Pz_S_D2z_Py_c;
  Double I_ERI_Px_Pz_D2z_Py_c = I_ERI_Dxz_S_D2z_Py_c+ABZ*I_ERI_Px_S_D2z_Py_c;
  Double I_ERI_Py_Pz_D2z_Py_c = I_ERI_Dyz_S_D2z_Py_c+ABZ*I_ERI_Py_S_D2z_Py_c;
  Double I_ERI_Pz_Pz_D2z_Py_c = I_ERI_D2z_S_D2z_Py_c+ABZ*I_ERI_Pz_S_D2z_Py_c;
  Double I_ERI_Px_Px_D2x_Pz_c = I_ERI_D2x_S_D2x_Pz_c+ABX*I_ERI_Px_S_D2x_Pz_c;
  Double I_ERI_Py_Px_D2x_Pz_c = I_ERI_Dxy_S_D2x_Pz_c+ABX*I_ERI_Py_S_D2x_Pz_c;
  Double I_ERI_Pz_Px_D2x_Pz_c = I_ERI_Dxz_S_D2x_Pz_c+ABX*I_ERI_Pz_S_D2x_Pz_c;
  Double I_ERI_Px_Py_D2x_Pz_c = I_ERI_Dxy_S_D2x_Pz_c+ABY*I_ERI_Px_S_D2x_Pz_c;
  Double I_ERI_Py_Py_D2x_Pz_c = I_ERI_D2y_S_D2x_Pz_c+ABY*I_ERI_Py_S_D2x_Pz_c;
  Double I_ERI_Pz_Py_D2x_Pz_c = I_ERI_Dyz_S_D2x_Pz_c+ABY*I_ERI_Pz_S_D2x_Pz_c;
  Double I_ERI_Px_Pz_D2x_Pz_c = I_ERI_Dxz_S_D2x_Pz_c+ABZ*I_ERI_Px_S_D2x_Pz_c;
  Double I_ERI_Py_Pz_D2x_Pz_c = I_ERI_Dyz_S_D2x_Pz_c+ABZ*I_ERI_Py_S_D2x_Pz_c;
  Double I_ERI_Pz_Pz_D2x_Pz_c = I_ERI_D2z_S_D2x_Pz_c+ABZ*I_ERI_Pz_S_D2x_Pz_c;
  Double I_ERI_Px_Px_Dxy_Pz_c = I_ERI_D2x_S_Dxy_Pz_c+ABX*I_ERI_Px_S_Dxy_Pz_c;
  Double I_ERI_Py_Px_Dxy_Pz_c = I_ERI_Dxy_S_Dxy_Pz_c+ABX*I_ERI_Py_S_Dxy_Pz_c;
  Double I_ERI_Pz_Px_Dxy_Pz_c = I_ERI_Dxz_S_Dxy_Pz_c+ABX*I_ERI_Pz_S_Dxy_Pz_c;
  Double I_ERI_Px_Py_Dxy_Pz_c = I_ERI_Dxy_S_Dxy_Pz_c+ABY*I_ERI_Px_S_Dxy_Pz_c;
  Double I_ERI_Py_Py_Dxy_Pz_c = I_ERI_D2y_S_Dxy_Pz_c+ABY*I_ERI_Py_S_Dxy_Pz_c;
  Double I_ERI_Pz_Py_Dxy_Pz_c = I_ERI_Dyz_S_Dxy_Pz_c+ABY*I_ERI_Pz_S_Dxy_Pz_c;
  Double I_ERI_Px_Pz_Dxy_Pz_c = I_ERI_Dxz_S_Dxy_Pz_c+ABZ*I_ERI_Px_S_Dxy_Pz_c;
  Double I_ERI_Py_Pz_Dxy_Pz_c = I_ERI_Dyz_S_Dxy_Pz_c+ABZ*I_ERI_Py_S_Dxy_Pz_c;
  Double I_ERI_Pz_Pz_Dxy_Pz_c = I_ERI_D2z_S_Dxy_Pz_c+ABZ*I_ERI_Pz_S_Dxy_Pz_c;
  Double I_ERI_Px_Px_Dxz_Pz_c = I_ERI_D2x_S_Dxz_Pz_c+ABX*I_ERI_Px_S_Dxz_Pz_c;
  Double I_ERI_Py_Px_Dxz_Pz_c = I_ERI_Dxy_S_Dxz_Pz_c+ABX*I_ERI_Py_S_Dxz_Pz_c;
  Double I_ERI_Pz_Px_Dxz_Pz_c = I_ERI_Dxz_S_Dxz_Pz_c+ABX*I_ERI_Pz_S_Dxz_Pz_c;
  Double I_ERI_Px_Py_Dxz_Pz_c = I_ERI_Dxy_S_Dxz_Pz_c+ABY*I_ERI_Px_S_Dxz_Pz_c;
  Double I_ERI_Py_Py_Dxz_Pz_c = I_ERI_D2y_S_Dxz_Pz_c+ABY*I_ERI_Py_S_Dxz_Pz_c;
  Double I_ERI_Pz_Py_Dxz_Pz_c = I_ERI_Dyz_S_Dxz_Pz_c+ABY*I_ERI_Pz_S_Dxz_Pz_c;
  Double I_ERI_Px_Pz_Dxz_Pz_c = I_ERI_Dxz_S_Dxz_Pz_c+ABZ*I_ERI_Px_S_Dxz_Pz_c;
  Double I_ERI_Py_Pz_Dxz_Pz_c = I_ERI_Dyz_S_Dxz_Pz_c+ABZ*I_ERI_Py_S_Dxz_Pz_c;
  Double I_ERI_Pz_Pz_Dxz_Pz_c = I_ERI_D2z_S_Dxz_Pz_c+ABZ*I_ERI_Pz_S_Dxz_Pz_c;
  Double I_ERI_Px_Px_D2y_Pz_c = I_ERI_D2x_S_D2y_Pz_c+ABX*I_ERI_Px_S_D2y_Pz_c;
  Double I_ERI_Py_Px_D2y_Pz_c = I_ERI_Dxy_S_D2y_Pz_c+ABX*I_ERI_Py_S_D2y_Pz_c;
  Double I_ERI_Pz_Px_D2y_Pz_c = I_ERI_Dxz_S_D2y_Pz_c+ABX*I_ERI_Pz_S_D2y_Pz_c;
  Double I_ERI_Px_Py_D2y_Pz_c = I_ERI_Dxy_S_D2y_Pz_c+ABY*I_ERI_Px_S_D2y_Pz_c;
  Double I_ERI_Py_Py_D2y_Pz_c = I_ERI_D2y_S_D2y_Pz_c+ABY*I_ERI_Py_S_D2y_Pz_c;
  Double I_ERI_Pz_Py_D2y_Pz_c = I_ERI_Dyz_S_D2y_Pz_c+ABY*I_ERI_Pz_S_D2y_Pz_c;
  Double I_ERI_Px_Pz_D2y_Pz_c = I_ERI_Dxz_S_D2y_Pz_c+ABZ*I_ERI_Px_S_D2y_Pz_c;
  Double I_ERI_Py_Pz_D2y_Pz_c = I_ERI_Dyz_S_D2y_Pz_c+ABZ*I_ERI_Py_S_D2y_Pz_c;
  Double I_ERI_Pz_Pz_D2y_Pz_c = I_ERI_D2z_S_D2y_Pz_c+ABZ*I_ERI_Pz_S_D2y_Pz_c;
  Double I_ERI_Px_Px_Dyz_Pz_c = I_ERI_D2x_S_Dyz_Pz_c+ABX*I_ERI_Px_S_Dyz_Pz_c;
  Double I_ERI_Py_Px_Dyz_Pz_c = I_ERI_Dxy_S_Dyz_Pz_c+ABX*I_ERI_Py_S_Dyz_Pz_c;
  Double I_ERI_Pz_Px_Dyz_Pz_c = I_ERI_Dxz_S_Dyz_Pz_c+ABX*I_ERI_Pz_S_Dyz_Pz_c;
  Double I_ERI_Px_Py_Dyz_Pz_c = I_ERI_Dxy_S_Dyz_Pz_c+ABY*I_ERI_Px_S_Dyz_Pz_c;
  Double I_ERI_Py_Py_Dyz_Pz_c = I_ERI_D2y_S_Dyz_Pz_c+ABY*I_ERI_Py_S_Dyz_Pz_c;
  Double I_ERI_Pz_Py_Dyz_Pz_c = I_ERI_Dyz_S_Dyz_Pz_c+ABY*I_ERI_Pz_S_Dyz_Pz_c;
  Double I_ERI_Px_Pz_Dyz_Pz_c = I_ERI_Dxz_S_Dyz_Pz_c+ABZ*I_ERI_Px_S_Dyz_Pz_c;
  Double I_ERI_Py_Pz_Dyz_Pz_c = I_ERI_Dyz_S_Dyz_Pz_c+ABZ*I_ERI_Py_S_Dyz_Pz_c;
  Double I_ERI_Pz_Pz_Dyz_Pz_c = I_ERI_D2z_S_Dyz_Pz_c+ABZ*I_ERI_Pz_S_Dyz_Pz_c;
  Double I_ERI_Px_Px_D2z_Pz_c = I_ERI_D2x_S_D2z_Pz_c+ABX*I_ERI_Px_S_D2z_Pz_c;
  Double I_ERI_Py_Px_D2z_Pz_c = I_ERI_Dxy_S_D2z_Pz_c+ABX*I_ERI_Py_S_D2z_Pz_c;
  Double I_ERI_Pz_Px_D2z_Pz_c = I_ERI_Dxz_S_D2z_Pz_c+ABX*I_ERI_Pz_S_D2z_Pz_c;
  Double I_ERI_Px_Py_D2z_Pz_c = I_ERI_Dxy_S_D2z_Pz_c+ABY*I_ERI_Px_S_D2z_Pz_c;
  Double I_ERI_Py_Py_D2z_Pz_c = I_ERI_D2y_S_D2z_Pz_c+ABY*I_ERI_Py_S_D2z_Pz_c;
  Double I_ERI_Pz_Py_D2z_Pz_c = I_ERI_Dyz_S_D2z_Pz_c+ABY*I_ERI_Pz_S_D2z_Pz_c;
  Double I_ERI_Px_Pz_D2z_Pz_c = I_ERI_Dxz_S_D2z_Pz_c+ABZ*I_ERI_Px_S_D2z_Pz_c;
  Double I_ERI_Py_Pz_D2z_Pz_c = I_ERI_Dyz_S_D2z_Pz_c+ABZ*I_ERI_Py_S_D2z_Pz_c;
  Double I_ERI_Pz_Pz_D2z_Pz_c = I_ERI_D2z_S_D2z_Pz_c+ABZ*I_ERI_Pz_S_D2z_Pz_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_D2x_Px_Px_Px_a-1*I_ERI_S_Px_Px_Px;
  abcd[1] = 2.0E0*I_ERI_Dxy_Px_Px_Px_a;
  abcd[2] = 2.0E0*I_ERI_Dxz_Px_Px_Px_a;
  abcd[3] = 2.0E0*I_ERI_D2x_Py_Px_Px_a-1*I_ERI_S_Py_Px_Px;
  abcd[4] = 2.0E0*I_ERI_Dxy_Py_Px_Px_a;
  abcd[5] = 2.0E0*I_ERI_Dxz_Py_Px_Px_a;
  abcd[6] = 2.0E0*I_ERI_D2x_Pz_Px_Px_a-1*I_ERI_S_Pz_Px_Px;
  abcd[7] = 2.0E0*I_ERI_Dxy_Pz_Px_Px_a;
  abcd[8] = 2.0E0*I_ERI_Dxz_Pz_Px_Px_a;
  abcd[9] = 2.0E0*I_ERI_D2x_Px_Py_Px_a-1*I_ERI_S_Px_Py_Px;
  abcd[10] = 2.0E0*I_ERI_Dxy_Px_Py_Px_a;
  abcd[11] = 2.0E0*I_ERI_Dxz_Px_Py_Px_a;
  abcd[12] = 2.0E0*I_ERI_D2x_Py_Py_Px_a-1*I_ERI_S_Py_Py_Px;
  abcd[13] = 2.0E0*I_ERI_Dxy_Py_Py_Px_a;
  abcd[14] = 2.0E0*I_ERI_Dxz_Py_Py_Px_a;
  abcd[15] = 2.0E0*I_ERI_D2x_Pz_Py_Px_a-1*I_ERI_S_Pz_Py_Px;
  abcd[16] = 2.0E0*I_ERI_Dxy_Pz_Py_Px_a;
  abcd[17] = 2.0E0*I_ERI_Dxz_Pz_Py_Px_a;
  abcd[18] = 2.0E0*I_ERI_D2x_Px_Pz_Px_a-1*I_ERI_S_Px_Pz_Px;
  abcd[19] = 2.0E0*I_ERI_Dxy_Px_Pz_Px_a;
  abcd[20] = 2.0E0*I_ERI_Dxz_Px_Pz_Px_a;
  abcd[21] = 2.0E0*I_ERI_D2x_Py_Pz_Px_a-1*I_ERI_S_Py_Pz_Px;
  abcd[22] = 2.0E0*I_ERI_Dxy_Py_Pz_Px_a;
  abcd[23] = 2.0E0*I_ERI_Dxz_Py_Pz_Px_a;
  abcd[24] = 2.0E0*I_ERI_D2x_Pz_Pz_Px_a-1*I_ERI_S_Pz_Pz_Px;
  abcd[25] = 2.0E0*I_ERI_Dxy_Pz_Pz_Px_a;
  abcd[26] = 2.0E0*I_ERI_Dxz_Pz_Pz_Px_a;
  abcd[27] = 2.0E0*I_ERI_D2x_Px_Px_Py_a-1*I_ERI_S_Px_Px_Py;
  abcd[28] = 2.0E0*I_ERI_Dxy_Px_Px_Py_a;
  abcd[29] = 2.0E0*I_ERI_Dxz_Px_Px_Py_a;
  abcd[30] = 2.0E0*I_ERI_D2x_Py_Px_Py_a-1*I_ERI_S_Py_Px_Py;
  abcd[31] = 2.0E0*I_ERI_Dxy_Py_Px_Py_a;
  abcd[32] = 2.0E0*I_ERI_Dxz_Py_Px_Py_a;
  abcd[33] = 2.0E0*I_ERI_D2x_Pz_Px_Py_a-1*I_ERI_S_Pz_Px_Py;
  abcd[34] = 2.0E0*I_ERI_Dxy_Pz_Px_Py_a;
  abcd[35] = 2.0E0*I_ERI_Dxz_Pz_Px_Py_a;
  abcd[36] = 2.0E0*I_ERI_D2x_Px_Py_Py_a-1*I_ERI_S_Px_Py_Py;
  abcd[37] = 2.0E0*I_ERI_Dxy_Px_Py_Py_a;
  abcd[38] = 2.0E0*I_ERI_Dxz_Px_Py_Py_a;
  abcd[39] = 2.0E0*I_ERI_D2x_Py_Py_Py_a-1*I_ERI_S_Py_Py_Py;
  abcd[40] = 2.0E0*I_ERI_Dxy_Py_Py_Py_a;
  abcd[41] = 2.0E0*I_ERI_Dxz_Py_Py_Py_a;
  abcd[42] = 2.0E0*I_ERI_D2x_Pz_Py_Py_a-1*I_ERI_S_Pz_Py_Py;
  abcd[43] = 2.0E0*I_ERI_Dxy_Pz_Py_Py_a;
  abcd[44] = 2.0E0*I_ERI_Dxz_Pz_Py_Py_a;
  abcd[45] = 2.0E0*I_ERI_D2x_Px_Pz_Py_a-1*I_ERI_S_Px_Pz_Py;
  abcd[46] = 2.0E0*I_ERI_Dxy_Px_Pz_Py_a;
  abcd[47] = 2.0E0*I_ERI_Dxz_Px_Pz_Py_a;
  abcd[48] = 2.0E0*I_ERI_D2x_Py_Pz_Py_a-1*I_ERI_S_Py_Pz_Py;
  abcd[49] = 2.0E0*I_ERI_Dxy_Py_Pz_Py_a;
  abcd[50] = 2.0E0*I_ERI_Dxz_Py_Pz_Py_a;
  abcd[51] = 2.0E0*I_ERI_D2x_Pz_Pz_Py_a-1*I_ERI_S_Pz_Pz_Py;
  abcd[52] = 2.0E0*I_ERI_Dxy_Pz_Pz_Py_a;
  abcd[53] = 2.0E0*I_ERI_Dxz_Pz_Pz_Py_a;
  abcd[54] = 2.0E0*I_ERI_D2x_Px_Px_Pz_a-1*I_ERI_S_Px_Px_Pz;
  abcd[55] = 2.0E0*I_ERI_Dxy_Px_Px_Pz_a;
  abcd[56] = 2.0E0*I_ERI_Dxz_Px_Px_Pz_a;
  abcd[57] = 2.0E0*I_ERI_D2x_Py_Px_Pz_a-1*I_ERI_S_Py_Px_Pz;
  abcd[58] = 2.0E0*I_ERI_Dxy_Py_Px_Pz_a;
  abcd[59] = 2.0E0*I_ERI_Dxz_Py_Px_Pz_a;
  abcd[60] = 2.0E0*I_ERI_D2x_Pz_Px_Pz_a-1*I_ERI_S_Pz_Px_Pz;
  abcd[61] = 2.0E0*I_ERI_Dxy_Pz_Px_Pz_a;
  abcd[62] = 2.0E0*I_ERI_Dxz_Pz_Px_Pz_a;
  abcd[63] = 2.0E0*I_ERI_D2x_Px_Py_Pz_a-1*I_ERI_S_Px_Py_Pz;
  abcd[64] = 2.0E0*I_ERI_Dxy_Px_Py_Pz_a;
  abcd[65] = 2.0E0*I_ERI_Dxz_Px_Py_Pz_a;
  abcd[66] = 2.0E0*I_ERI_D2x_Py_Py_Pz_a-1*I_ERI_S_Py_Py_Pz;
  abcd[67] = 2.0E0*I_ERI_Dxy_Py_Py_Pz_a;
  abcd[68] = 2.0E0*I_ERI_Dxz_Py_Py_Pz_a;
  abcd[69] = 2.0E0*I_ERI_D2x_Pz_Py_Pz_a-1*I_ERI_S_Pz_Py_Pz;
  abcd[70] = 2.0E0*I_ERI_Dxy_Pz_Py_Pz_a;
  abcd[71] = 2.0E0*I_ERI_Dxz_Pz_Py_Pz_a;
  abcd[72] = 2.0E0*I_ERI_D2x_Px_Pz_Pz_a-1*I_ERI_S_Px_Pz_Pz;
  abcd[73] = 2.0E0*I_ERI_Dxy_Px_Pz_Pz_a;
  abcd[74] = 2.0E0*I_ERI_Dxz_Px_Pz_Pz_a;
  abcd[75] = 2.0E0*I_ERI_D2x_Py_Pz_Pz_a-1*I_ERI_S_Py_Pz_Pz;
  abcd[76] = 2.0E0*I_ERI_Dxy_Py_Pz_Pz_a;
  abcd[77] = 2.0E0*I_ERI_Dxz_Py_Pz_Pz_a;
  abcd[78] = 2.0E0*I_ERI_D2x_Pz_Pz_Pz_a-1*I_ERI_S_Pz_Pz_Pz;
  abcd[79] = 2.0E0*I_ERI_Dxy_Pz_Pz_Pz_a;
  abcd[80] = 2.0E0*I_ERI_Dxz_Pz_Pz_Pz_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P
   ************************************************************/
  abcd[81] = 2.0E0*I_ERI_Dxy_Px_Px_Px_a;
  abcd[82] = 2.0E0*I_ERI_D2y_Px_Px_Px_a-1*I_ERI_S_Px_Px_Px;
  abcd[83] = 2.0E0*I_ERI_Dyz_Px_Px_Px_a;
  abcd[84] = 2.0E0*I_ERI_Dxy_Py_Px_Px_a;
  abcd[85] = 2.0E0*I_ERI_D2y_Py_Px_Px_a-1*I_ERI_S_Py_Px_Px;
  abcd[86] = 2.0E0*I_ERI_Dyz_Py_Px_Px_a;
  abcd[87] = 2.0E0*I_ERI_Dxy_Pz_Px_Px_a;
  abcd[88] = 2.0E0*I_ERI_D2y_Pz_Px_Px_a-1*I_ERI_S_Pz_Px_Px;
  abcd[89] = 2.0E0*I_ERI_Dyz_Pz_Px_Px_a;
  abcd[90] = 2.0E0*I_ERI_Dxy_Px_Py_Px_a;
  abcd[91] = 2.0E0*I_ERI_D2y_Px_Py_Px_a-1*I_ERI_S_Px_Py_Px;
  abcd[92] = 2.0E0*I_ERI_Dyz_Px_Py_Px_a;
  abcd[93] = 2.0E0*I_ERI_Dxy_Py_Py_Px_a;
  abcd[94] = 2.0E0*I_ERI_D2y_Py_Py_Px_a-1*I_ERI_S_Py_Py_Px;
  abcd[95] = 2.0E0*I_ERI_Dyz_Py_Py_Px_a;
  abcd[96] = 2.0E0*I_ERI_Dxy_Pz_Py_Px_a;
  abcd[97] = 2.0E0*I_ERI_D2y_Pz_Py_Px_a-1*I_ERI_S_Pz_Py_Px;
  abcd[98] = 2.0E0*I_ERI_Dyz_Pz_Py_Px_a;
  abcd[99] = 2.0E0*I_ERI_Dxy_Px_Pz_Px_a;
  abcd[100] = 2.0E0*I_ERI_D2y_Px_Pz_Px_a-1*I_ERI_S_Px_Pz_Px;
  abcd[101] = 2.0E0*I_ERI_Dyz_Px_Pz_Px_a;
  abcd[102] = 2.0E0*I_ERI_Dxy_Py_Pz_Px_a;
  abcd[103] = 2.0E0*I_ERI_D2y_Py_Pz_Px_a-1*I_ERI_S_Py_Pz_Px;
  abcd[104] = 2.0E0*I_ERI_Dyz_Py_Pz_Px_a;
  abcd[105] = 2.0E0*I_ERI_Dxy_Pz_Pz_Px_a;
  abcd[106] = 2.0E0*I_ERI_D2y_Pz_Pz_Px_a-1*I_ERI_S_Pz_Pz_Px;
  abcd[107] = 2.0E0*I_ERI_Dyz_Pz_Pz_Px_a;
  abcd[108] = 2.0E0*I_ERI_Dxy_Px_Px_Py_a;
  abcd[109] = 2.0E0*I_ERI_D2y_Px_Px_Py_a-1*I_ERI_S_Px_Px_Py;
  abcd[110] = 2.0E0*I_ERI_Dyz_Px_Px_Py_a;
  abcd[111] = 2.0E0*I_ERI_Dxy_Py_Px_Py_a;
  abcd[112] = 2.0E0*I_ERI_D2y_Py_Px_Py_a-1*I_ERI_S_Py_Px_Py;
  abcd[113] = 2.0E0*I_ERI_Dyz_Py_Px_Py_a;
  abcd[114] = 2.0E0*I_ERI_Dxy_Pz_Px_Py_a;
  abcd[115] = 2.0E0*I_ERI_D2y_Pz_Px_Py_a-1*I_ERI_S_Pz_Px_Py;
  abcd[116] = 2.0E0*I_ERI_Dyz_Pz_Px_Py_a;
  abcd[117] = 2.0E0*I_ERI_Dxy_Px_Py_Py_a;
  abcd[118] = 2.0E0*I_ERI_D2y_Px_Py_Py_a-1*I_ERI_S_Px_Py_Py;
  abcd[119] = 2.0E0*I_ERI_Dyz_Px_Py_Py_a;
  abcd[120] = 2.0E0*I_ERI_Dxy_Py_Py_Py_a;
  abcd[121] = 2.0E0*I_ERI_D2y_Py_Py_Py_a-1*I_ERI_S_Py_Py_Py;
  abcd[122] = 2.0E0*I_ERI_Dyz_Py_Py_Py_a;
  abcd[123] = 2.0E0*I_ERI_Dxy_Pz_Py_Py_a;
  abcd[124] = 2.0E0*I_ERI_D2y_Pz_Py_Py_a-1*I_ERI_S_Pz_Py_Py;
  abcd[125] = 2.0E0*I_ERI_Dyz_Pz_Py_Py_a;
  abcd[126] = 2.0E0*I_ERI_Dxy_Px_Pz_Py_a;
  abcd[127] = 2.0E0*I_ERI_D2y_Px_Pz_Py_a-1*I_ERI_S_Px_Pz_Py;
  abcd[128] = 2.0E0*I_ERI_Dyz_Px_Pz_Py_a;
  abcd[129] = 2.0E0*I_ERI_Dxy_Py_Pz_Py_a;
  abcd[130] = 2.0E0*I_ERI_D2y_Py_Pz_Py_a-1*I_ERI_S_Py_Pz_Py;
  abcd[131] = 2.0E0*I_ERI_Dyz_Py_Pz_Py_a;
  abcd[132] = 2.0E0*I_ERI_Dxy_Pz_Pz_Py_a;
  abcd[133] = 2.0E0*I_ERI_D2y_Pz_Pz_Py_a-1*I_ERI_S_Pz_Pz_Py;
  abcd[134] = 2.0E0*I_ERI_Dyz_Pz_Pz_Py_a;
  abcd[135] = 2.0E0*I_ERI_Dxy_Px_Px_Pz_a;
  abcd[136] = 2.0E0*I_ERI_D2y_Px_Px_Pz_a-1*I_ERI_S_Px_Px_Pz;
  abcd[137] = 2.0E0*I_ERI_Dyz_Px_Px_Pz_a;
  abcd[138] = 2.0E0*I_ERI_Dxy_Py_Px_Pz_a;
  abcd[139] = 2.0E0*I_ERI_D2y_Py_Px_Pz_a-1*I_ERI_S_Py_Px_Pz;
  abcd[140] = 2.0E0*I_ERI_Dyz_Py_Px_Pz_a;
  abcd[141] = 2.0E0*I_ERI_Dxy_Pz_Px_Pz_a;
  abcd[142] = 2.0E0*I_ERI_D2y_Pz_Px_Pz_a-1*I_ERI_S_Pz_Px_Pz;
  abcd[143] = 2.0E0*I_ERI_Dyz_Pz_Px_Pz_a;
  abcd[144] = 2.0E0*I_ERI_Dxy_Px_Py_Pz_a;
  abcd[145] = 2.0E0*I_ERI_D2y_Px_Py_Pz_a-1*I_ERI_S_Px_Py_Pz;
  abcd[146] = 2.0E0*I_ERI_Dyz_Px_Py_Pz_a;
  abcd[147] = 2.0E0*I_ERI_Dxy_Py_Py_Pz_a;
  abcd[148] = 2.0E0*I_ERI_D2y_Py_Py_Pz_a-1*I_ERI_S_Py_Py_Pz;
  abcd[149] = 2.0E0*I_ERI_Dyz_Py_Py_Pz_a;
  abcd[150] = 2.0E0*I_ERI_Dxy_Pz_Py_Pz_a;
  abcd[151] = 2.0E0*I_ERI_D2y_Pz_Py_Pz_a-1*I_ERI_S_Pz_Py_Pz;
  abcd[152] = 2.0E0*I_ERI_Dyz_Pz_Py_Pz_a;
  abcd[153] = 2.0E0*I_ERI_Dxy_Px_Pz_Pz_a;
  abcd[154] = 2.0E0*I_ERI_D2y_Px_Pz_Pz_a-1*I_ERI_S_Px_Pz_Pz;
  abcd[155] = 2.0E0*I_ERI_Dyz_Px_Pz_Pz_a;
  abcd[156] = 2.0E0*I_ERI_Dxy_Py_Pz_Pz_a;
  abcd[157] = 2.0E0*I_ERI_D2y_Py_Pz_Pz_a-1*I_ERI_S_Py_Pz_Pz;
  abcd[158] = 2.0E0*I_ERI_Dyz_Py_Pz_Pz_a;
  abcd[159] = 2.0E0*I_ERI_Dxy_Pz_Pz_Pz_a;
  abcd[160] = 2.0E0*I_ERI_D2y_Pz_Pz_Pz_a-1*I_ERI_S_Pz_Pz_Pz;
  abcd[161] = 2.0E0*I_ERI_Dyz_Pz_Pz_Pz_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_P_a
   * RHS shell quartet name: SQ_ERI_S_P_P_P
   ************************************************************/
  abcd[162] = 2.0E0*I_ERI_Dxz_Px_Px_Px_a;
  abcd[163] = 2.0E0*I_ERI_Dyz_Px_Px_Px_a;
  abcd[164] = 2.0E0*I_ERI_D2z_Px_Px_Px_a-1*I_ERI_S_Px_Px_Px;
  abcd[165] = 2.0E0*I_ERI_Dxz_Py_Px_Px_a;
  abcd[166] = 2.0E0*I_ERI_Dyz_Py_Px_Px_a;
  abcd[167] = 2.0E0*I_ERI_D2z_Py_Px_Px_a-1*I_ERI_S_Py_Px_Px;
  abcd[168] = 2.0E0*I_ERI_Dxz_Pz_Px_Px_a;
  abcd[169] = 2.0E0*I_ERI_Dyz_Pz_Px_Px_a;
  abcd[170] = 2.0E0*I_ERI_D2z_Pz_Px_Px_a-1*I_ERI_S_Pz_Px_Px;
  abcd[171] = 2.0E0*I_ERI_Dxz_Px_Py_Px_a;
  abcd[172] = 2.0E0*I_ERI_Dyz_Px_Py_Px_a;
  abcd[173] = 2.0E0*I_ERI_D2z_Px_Py_Px_a-1*I_ERI_S_Px_Py_Px;
  abcd[174] = 2.0E0*I_ERI_Dxz_Py_Py_Px_a;
  abcd[175] = 2.0E0*I_ERI_Dyz_Py_Py_Px_a;
  abcd[176] = 2.0E0*I_ERI_D2z_Py_Py_Px_a-1*I_ERI_S_Py_Py_Px;
  abcd[177] = 2.0E0*I_ERI_Dxz_Pz_Py_Px_a;
  abcd[178] = 2.0E0*I_ERI_Dyz_Pz_Py_Px_a;
  abcd[179] = 2.0E0*I_ERI_D2z_Pz_Py_Px_a-1*I_ERI_S_Pz_Py_Px;
  abcd[180] = 2.0E0*I_ERI_Dxz_Px_Pz_Px_a;
  abcd[181] = 2.0E0*I_ERI_Dyz_Px_Pz_Px_a;
  abcd[182] = 2.0E0*I_ERI_D2z_Px_Pz_Px_a-1*I_ERI_S_Px_Pz_Px;
  abcd[183] = 2.0E0*I_ERI_Dxz_Py_Pz_Px_a;
  abcd[184] = 2.0E0*I_ERI_Dyz_Py_Pz_Px_a;
  abcd[185] = 2.0E0*I_ERI_D2z_Py_Pz_Px_a-1*I_ERI_S_Py_Pz_Px;
  abcd[186] = 2.0E0*I_ERI_Dxz_Pz_Pz_Px_a;
  abcd[187] = 2.0E0*I_ERI_Dyz_Pz_Pz_Px_a;
  abcd[188] = 2.0E0*I_ERI_D2z_Pz_Pz_Px_a-1*I_ERI_S_Pz_Pz_Px;
  abcd[189] = 2.0E0*I_ERI_Dxz_Px_Px_Py_a;
  abcd[190] = 2.0E0*I_ERI_Dyz_Px_Px_Py_a;
  abcd[191] = 2.0E0*I_ERI_D2z_Px_Px_Py_a-1*I_ERI_S_Px_Px_Py;
  abcd[192] = 2.0E0*I_ERI_Dxz_Py_Px_Py_a;
  abcd[193] = 2.0E0*I_ERI_Dyz_Py_Px_Py_a;
  abcd[194] = 2.0E0*I_ERI_D2z_Py_Px_Py_a-1*I_ERI_S_Py_Px_Py;
  abcd[195] = 2.0E0*I_ERI_Dxz_Pz_Px_Py_a;
  abcd[196] = 2.0E0*I_ERI_Dyz_Pz_Px_Py_a;
  abcd[197] = 2.0E0*I_ERI_D2z_Pz_Px_Py_a-1*I_ERI_S_Pz_Px_Py;
  abcd[198] = 2.0E0*I_ERI_Dxz_Px_Py_Py_a;
  abcd[199] = 2.0E0*I_ERI_Dyz_Px_Py_Py_a;
  abcd[200] = 2.0E0*I_ERI_D2z_Px_Py_Py_a-1*I_ERI_S_Px_Py_Py;
  abcd[201] = 2.0E0*I_ERI_Dxz_Py_Py_Py_a;
  abcd[202] = 2.0E0*I_ERI_Dyz_Py_Py_Py_a;
  abcd[203] = 2.0E0*I_ERI_D2z_Py_Py_Py_a-1*I_ERI_S_Py_Py_Py;
  abcd[204] = 2.0E0*I_ERI_Dxz_Pz_Py_Py_a;
  abcd[205] = 2.0E0*I_ERI_Dyz_Pz_Py_Py_a;
  abcd[206] = 2.0E0*I_ERI_D2z_Pz_Py_Py_a-1*I_ERI_S_Pz_Py_Py;
  abcd[207] = 2.0E0*I_ERI_Dxz_Px_Pz_Py_a;
  abcd[208] = 2.0E0*I_ERI_Dyz_Px_Pz_Py_a;
  abcd[209] = 2.0E0*I_ERI_D2z_Px_Pz_Py_a-1*I_ERI_S_Px_Pz_Py;
  abcd[210] = 2.0E0*I_ERI_Dxz_Py_Pz_Py_a;
  abcd[211] = 2.0E0*I_ERI_Dyz_Py_Pz_Py_a;
  abcd[212] = 2.0E0*I_ERI_D2z_Py_Pz_Py_a-1*I_ERI_S_Py_Pz_Py;
  abcd[213] = 2.0E0*I_ERI_Dxz_Pz_Pz_Py_a;
  abcd[214] = 2.0E0*I_ERI_Dyz_Pz_Pz_Py_a;
  abcd[215] = 2.0E0*I_ERI_D2z_Pz_Pz_Py_a-1*I_ERI_S_Pz_Pz_Py;
  abcd[216] = 2.0E0*I_ERI_Dxz_Px_Px_Pz_a;
  abcd[217] = 2.0E0*I_ERI_Dyz_Px_Px_Pz_a;
  abcd[218] = 2.0E0*I_ERI_D2z_Px_Px_Pz_a-1*I_ERI_S_Px_Px_Pz;
  abcd[219] = 2.0E0*I_ERI_Dxz_Py_Px_Pz_a;
  abcd[220] = 2.0E0*I_ERI_Dyz_Py_Px_Pz_a;
  abcd[221] = 2.0E0*I_ERI_D2z_Py_Px_Pz_a-1*I_ERI_S_Py_Px_Pz;
  abcd[222] = 2.0E0*I_ERI_Dxz_Pz_Px_Pz_a;
  abcd[223] = 2.0E0*I_ERI_Dyz_Pz_Px_Pz_a;
  abcd[224] = 2.0E0*I_ERI_D2z_Pz_Px_Pz_a-1*I_ERI_S_Pz_Px_Pz;
  abcd[225] = 2.0E0*I_ERI_Dxz_Px_Py_Pz_a;
  abcd[226] = 2.0E0*I_ERI_Dyz_Px_Py_Pz_a;
  abcd[227] = 2.0E0*I_ERI_D2z_Px_Py_Pz_a-1*I_ERI_S_Px_Py_Pz;
  abcd[228] = 2.0E0*I_ERI_Dxz_Py_Py_Pz_a;
  abcd[229] = 2.0E0*I_ERI_Dyz_Py_Py_Pz_a;
  abcd[230] = 2.0E0*I_ERI_D2z_Py_Py_Pz_a-1*I_ERI_S_Py_Py_Pz;
  abcd[231] = 2.0E0*I_ERI_Dxz_Pz_Py_Pz_a;
  abcd[232] = 2.0E0*I_ERI_Dyz_Pz_Py_Pz_a;
  abcd[233] = 2.0E0*I_ERI_D2z_Pz_Py_Pz_a-1*I_ERI_S_Pz_Py_Pz;
  abcd[234] = 2.0E0*I_ERI_Dxz_Px_Pz_Pz_a;
  abcd[235] = 2.0E0*I_ERI_Dyz_Px_Pz_Pz_a;
  abcd[236] = 2.0E0*I_ERI_D2z_Px_Pz_Pz_a-1*I_ERI_S_Px_Pz_Pz;
  abcd[237] = 2.0E0*I_ERI_Dxz_Py_Pz_Pz_a;
  abcd[238] = 2.0E0*I_ERI_Dyz_Py_Pz_Pz_a;
  abcd[239] = 2.0E0*I_ERI_D2z_Py_Pz_Pz_a-1*I_ERI_S_Py_Pz_Pz;
  abcd[240] = 2.0E0*I_ERI_Dxz_Pz_Pz_Pz_a;
  abcd[241] = 2.0E0*I_ERI_Dyz_Pz_Pz_Pz_a;
  abcd[242] = 2.0E0*I_ERI_D2z_Pz_Pz_Pz_a-1*I_ERI_S_Pz_Pz_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P
   ************************************************************/
  abcd[243] = 2.0E0*I_ERI_Px_D2x_Px_Px_b-1*I_ERI_Px_S_Px_Px;
  abcd[244] = 2.0E0*I_ERI_Py_D2x_Px_Px_b-1*I_ERI_Py_S_Px_Px;
  abcd[245] = 2.0E0*I_ERI_Pz_D2x_Px_Px_b-1*I_ERI_Pz_S_Px_Px;
  abcd[246] = 2.0E0*I_ERI_Px_Dxy_Px_Px_b;
  abcd[247] = 2.0E0*I_ERI_Py_Dxy_Px_Px_b;
  abcd[248] = 2.0E0*I_ERI_Pz_Dxy_Px_Px_b;
  abcd[249] = 2.0E0*I_ERI_Px_Dxz_Px_Px_b;
  abcd[250] = 2.0E0*I_ERI_Py_Dxz_Px_Px_b;
  abcd[251] = 2.0E0*I_ERI_Pz_Dxz_Px_Px_b;
  abcd[252] = 2.0E0*I_ERI_Px_D2x_Py_Px_b-1*I_ERI_Px_S_Py_Px;
  abcd[253] = 2.0E0*I_ERI_Py_D2x_Py_Px_b-1*I_ERI_Py_S_Py_Px;
  abcd[254] = 2.0E0*I_ERI_Pz_D2x_Py_Px_b-1*I_ERI_Pz_S_Py_Px;
  abcd[255] = 2.0E0*I_ERI_Px_Dxy_Py_Px_b;
  abcd[256] = 2.0E0*I_ERI_Py_Dxy_Py_Px_b;
  abcd[257] = 2.0E0*I_ERI_Pz_Dxy_Py_Px_b;
  abcd[258] = 2.0E0*I_ERI_Px_Dxz_Py_Px_b;
  abcd[259] = 2.0E0*I_ERI_Py_Dxz_Py_Px_b;
  abcd[260] = 2.0E0*I_ERI_Pz_Dxz_Py_Px_b;
  abcd[261] = 2.0E0*I_ERI_Px_D2x_Pz_Px_b-1*I_ERI_Px_S_Pz_Px;
  abcd[262] = 2.0E0*I_ERI_Py_D2x_Pz_Px_b-1*I_ERI_Py_S_Pz_Px;
  abcd[263] = 2.0E0*I_ERI_Pz_D2x_Pz_Px_b-1*I_ERI_Pz_S_Pz_Px;
  abcd[264] = 2.0E0*I_ERI_Px_Dxy_Pz_Px_b;
  abcd[265] = 2.0E0*I_ERI_Py_Dxy_Pz_Px_b;
  abcd[266] = 2.0E0*I_ERI_Pz_Dxy_Pz_Px_b;
  abcd[267] = 2.0E0*I_ERI_Px_Dxz_Pz_Px_b;
  abcd[268] = 2.0E0*I_ERI_Py_Dxz_Pz_Px_b;
  abcd[269] = 2.0E0*I_ERI_Pz_Dxz_Pz_Px_b;
  abcd[270] = 2.0E0*I_ERI_Px_D2x_Px_Py_b-1*I_ERI_Px_S_Px_Py;
  abcd[271] = 2.0E0*I_ERI_Py_D2x_Px_Py_b-1*I_ERI_Py_S_Px_Py;
  abcd[272] = 2.0E0*I_ERI_Pz_D2x_Px_Py_b-1*I_ERI_Pz_S_Px_Py;
  abcd[273] = 2.0E0*I_ERI_Px_Dxy_Px_Py_b;
  abcd[274] = 2.0E0*I_ERI_Py_Dxy_Px_Py_b;
  abcd[275] = 2.0E0*I_ERI_Pz_Dxy_Px_Py_b;
  abcd[276] = 2.0E0*I_ERI_Px_Dxz_Px_Py_b;
  abcd[277] = 2.0E0*I_ERI_Py_Dxz_Px_Py_b;
  abcd[278] = 2.0E0*I_ERI_Pz_Dxz_Px_Py_b;
  abcd[279] = 2.0E0*I_ERI_Px_D2x_Py_Py_b-1*I_ERI_Px_S_Py_Py;
  abcd[280] = 2.0E0*I_ERI_Py_D2x_Py_Py_b-1*I_ERI_Py_S_Py_Py;
  abcd[281] = 2.0E0*I_ERI_Pz_D2x_Py_Py_b-1*I_ERI_Pz_S_Py_Py;
  abcd[282] = 2.0E0*I_ERI_Px_Dxy_Py_Py_b;
  abcd[283] = 2.0E0*I_ERI_Py_Dxy_Py_Py_b;
  abcd[284] = 2.0E0*I_ERI_Pz_Dxy_Py_Py_b;
  abcd[285] = 2.0E0*I_ERI_Px_Dxz_Py_Py_b;
  abcd[286] = 2.0E0*I_ERI_Py_Dxz_Py_Py_b;
  abcd[287] = 2.0E0*I_ERI_Pz_Dxz_Py_Py_b;
  abcd[288] = 2.0E0*I_ERI_Px_D2x_Pz_Py_b-1*I_ERI_Px_S_Pz_Py;
  abcd[289] = 2.0E0*I_ERI_Py_D2x_Pz_Py_b-1*I_ERI_Py_S_Pz_Py;
  abcd[290] = 2.0E0*I_ERI_Pz_D2x_Pz_Py_b-1*I_ERI_Pz_S_Pz_Py;
  abcd[291] = 2.0E0*I_ERI_Px_Dxy_Pz_Py_b;
  abcd[292] = 2.0E0*I_ERI_Py_Dxy_Pz_Py_b;
  abcd[293] = 2.0E0*I_ERI_Pz_Dxy_Pz_Py_b;
  abcd[294] = 2.0E0*I_ERI_Px_Dxz_Pz_Py_b;
  abcd[295] = 2.0E0*I_ERI_Py_Dxz_Pz_Py_b;
  abcd[296] = 2.0E0*I_ERI_Pz_Dxz_Pz_Py_b;
  abcd[297] = 2.0E0*I_ERI_Px_D2x_Px_Pz_b-1*I_ERI_Px_S_Px_Pz;
  abcd[298] = 2.0E0*I_ERI_Py_D2x_Px_Pz_b-1*I_ERI_Py_S_Px_Pz;
  abcd[299] = 2.0E0*I_ERI_Pz_D2x_Px_Pz_b-1*I_ERI_Pz_S_Px_Pz;
  abcd[300] = 2.0E0*I_ERI_Px_Dxy_Px_Pz_b;
  abcd[301] = 2.0E0*I_ERI_Py_Dxy_Px_Pz_b;
  abcd[302] = 2.0E0*I_ERI_Pz_Dxy_Px_Pz_b;
  abcd[303] = 2.0E0*I_ERI_Px_Dxz_Px_Pz_b;
  abcd[304] = 2.0E0*I_ERI_Py_Dxz_Px_Pz_b;
  abcd[305] = 2.0E0*I_ERI_Pz_Dxz_Px_Pz_b;
  abcd[306] = 2.0E0*I_ERI_Px_D2x_Py_Pz_b-1*I_ERI_Px_S_Py_Pz;
  abcd[307] = 2.0E0*I_ERI_Py_D2x_Py_Pz_b-1*I_ERI_Py_S_Py_Pz;
  abcd[308] = 2.0E0*I_ERI_Pz_D2x_Py_Pz_b-1*I_ERI_Pz_S_Py_Pz;
  abcd[309] = 2.0E0*I_ERI_Px_Dxy_Py_Pz_b;
  abcd[310] = 2.0E0*I_ERI_Py_Dxy_Py_Pz_b;
  abcd[311] = 2.0E0*I_ERI_Pz_Dxy_Py_Pz_b;
  abcd[312] = 2.0E0*I_ERI_Px_Dxz_Py_Pz_b;
  abcd[313] = 2.0E0*I_ERI_Py_Dxz_Py_Pz_b;
  abcd[314] = 2.0E0*I_ERI_Pz_Dxz_Py_Pz_b;
  abcd[315] = 2.0E0*I_ERI_Px_D2x_Pz_Pz_b-1*I_ERI_Px_S_Pz_Pz;
  abcd[316] = 2.0E0*I_ERI_Py_D2x_Pz_Pz_b-1*I_ERI_Py_S_Pz_Pz;
  abcd[317] = 2.0E0*I_ERI_Pz_D2x_Pz_Pz_b-1*I_ERI_Pz_S_Pz_Pz;
  abcd[318] = 2.0E0*I_ERI_Px_Dxy_Pz_Pz_b;
  abcd[319] = 2.0E0*I_ERI_Py_Dxy_Pz_Pz_b;
  abcd[320] = 2.0E0*I_ERI_Pz_Dxy_Pz_Pz_b;
  abcd[321] = 2.0E0*I_ERI_Px_Dxz_Pz_Pz_b;
  abcd[322] = 2.0E0*I_ERI_Py_Dxz_Pz_Pz_b;
  abcd[323] = 2.0E0*I_ERI_Pz_Dxz_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P
   ************************************************************/
  abcd[324] = 2.0E0*I_ERI_Px_Dxy_Px_Px_b;
  abcd[325] = 2.0E0*I_ERI_Py_Dxy_Px_Px_b;
  abcd[326] = 2.0E0*I_ERI_Pz_Dxy_Px_Px_b;
  abcd[327] = 2.0E0*I_ERI_Px_D2y_Px_Px_b-1*I_ERI_Px_S_Px_Px;
  abcd[328] = 2.0E0*I_ERI_Py_D2y_Px_Px_b-1*I_ERI_Py_S_Px_Px;
  abcd[329] = 2.0E0*I_ERI_Pz_D2y_Px_Px_b-1*I_ERI_Pz_S_Px_Px;
  abcd[330] = 2.0E0*I_ERI_Px_Dyz_Px_Px_b;
  abcd[331] = 2.0E0*I_ERI_Py_Dyz_Px_Px_b;
  abcd[332] = 2.0E0*I_ERI_Pz_Dyz_Px_Px_b;
  abcd[333] = 2.0E0*I_ERI_Px_Dxy_Py_Px_b;
  abcd[334] = 2.0E0*I_ERI_Py_Dxy_Py_Px_b;
  abcd[335] = 2.0E0*I_ERI_Pz_Dxy_Py_Px_b;
  abcd[336] = 2.0E0*I_ERI_Px_D2y_Py_Px_b-1*I_ERI_Px_S_Py_Px;
  abcd[337] = 2.0E0*I_ERI_Py_D2y_Py_Px_b-1*I_ERI_Py_S_Py_Px;
  abcd[338] = 2.0E0*I_ERI_Pz_D2y_Py_Px_b-1*I_ERI_Pz_S_Py_Px;
  abcd[339] = 2.0E0*I_ERI_Px_Dyz_Py_Px_b;
  abcd[340] = 2.0E0*I_ERI_Py_Dyz_Py_Px_b;
  abcd[341] = 2.0E0*I_ERI_Pz_Dyz_Py_Px_b;
  abcd[342] = 2.0E0*I_ERI_Px_Dxy_Pz_Px_b;
  abcd[343] = 2.0E0*I_ERI_Py_Dxy_Pz_Px_b;
  abcd[344] = 2.0E0*I_ERI_Pz_Dxy_Pz_Px_b;
  abcd[345] = 2.0E0*I_ERI_Px_D2y_Pz_Px_b-1*I_ERI_Px_S_Pz_Px;
  abcd[346] = 2.0E0*I_ERI_Py_D2y_Pz_Px_b-1*I_ERI_Py_S_Pz_Px;
  abcd[347] = 2.0E0*I_ERI_Pz_D2y_Pz_Px_b-1*I_ERI_Pz_S_Pz_Px;
  abcd[348] = 2.0E0*I_ERI_Px_Dyz_Pz_Px_b;
  abcd[349] = 2.0E0*I_ERI_Py_Dyz_Pz_Px_b;
  abcd[350] = 2.0E0*I_ERI_Pz_Dyz_Pz_Px_b;
  abcd[351] = 2.0E0*I_ERI_Px_Dxy_Px_Py_b;
  abcd[352] = 2.0E0*I_ERI_Py_Dxy_Px_Py_b;
  abcd[353] = 2.0E0*I_ERI_Pz_Dxy_Px_Py_b;
  abcd[354] = 2.0E0*I_ERI_Px_D2y_Px_Py_b-1*I_ERI_Px_S_Px_Py;
  abcd[355] = 2.0E0*I_ERI_Py_D2y_Px_Py_b-1*I_ERI_Py_S_Px_Py;
  abcd[356] = 2.0E0*I_ERI_Pz_D2y_Px_Py_b-1*I_ERI_Pz_S_Px_Py;
  abcd[357] = 2.0E0*I_ERI_Px_Dyz_Px_Py_b;
  abcd[358] = 2.0E0*I_ERI_Py_Dyz_Px_Py_b;
  abcd[359] = 2.0E0*I_ERI_Pz_Dyz_Px_Py_b;
  abcd[360] = 2.0E0*I_ERI_Px_Dxy_Py_Py_b;
  abcd[361] = 2.0E0*I_ERI_Py_Dxy_Py_Py_b;
  abcd[362] = 2.0E0*I_ERI_Pz_Dxy_Py_Py_b;
  abcd[363] = 2.0E0*I_ERI_Px_D2y_Py_Py_b-1*I_ERI_Px_S_Py_Py;
  abcd[364] = 2.0E0*I_ERI_Py_D2y_Py_Py_b-1*I_ERI_Py_S_Py_Py;
  abcd[365] = 2.0E0*I_ERI_Pz_D2y_Py_Py_b-1*I_ERI_Pz_S_Py_Py;
  abcd[366] = 2.0E0*I_ERI_Px_Dyz_Py_Py_b;
  abcd[367] = 2.0E0*I_ERI_Py_Dyz_Py_Py_b;
  abcd[368] = 2.0E0*I_ERI_Pz_Dyz_Py_Py_b;
  abcd[369] = 2.0E0*I_ERI_Px_Dxy_Pz_Py_b;
  abcd[370] = 2.0E0*I_ERI_Py_Dxy_Pz_Py_b;
  abcd[371] = 2.0E0*I_ERI_Pz_Dxy_Pz_Py_b;
  abcd[372] = 2.0E0*I_ERI_Px_D2y_Pz_Py_b-1*I_ERI_Px_S_Pz_Py;
  abcd[373] = 2.0E0*I_ERI_Py_D2y_Pz_Py_b-1*I_ERI_Py_S_Pz_Py;
  abcd[374] = 2.0E0*I_ERI_Pz_D2y_Pz_Py_b-1*I_ERI_Pz_S_Pz_Py;
  abcd[375] = 2.0E0*I_ERI_Px_Dyz_Pz_Py_b;
  abcd[376] = 2.0E0*I_ERI_Py_Dyz_Pz_Py_b;
  abcd[377] = 2.0E0*I_ERI_Pz_Dyz_Pz_Py_b;
  abcd[378] = 2.0E0*I_ERI_Px_Dxy_Px_Pz_b;
  abcd[379] = 2.0E0*I_ERI_Py_Dxy_Px_Pz_b;
  abcd[380] = 2.0E0*I_ERI_Pz_Dxy_Px_Pz_b;
  abcd[381] = 2.0E0*I_ERI_Px_D2y_Px_Pz_b-1*I_ERI_Px_S_Px_Pz;
  abcd[382] = 2.0E0*I_ERI_Py_D2y_Px_Pz_b-1*I_ERI_Py_S_Px_Pz;
  abcd[383] = 2.0E0*I_ERI_Pz_D2y_Px_Pz_b-1*I_ERI_Pz_S_Px_Pz;
  abcd[384] = 2.0E0*I_ERI_Px_Dyz_Px_Pz_b;
  abcd[385] = 2.0E0*I_ERI_Py_Dyz_Px_Pz_b;
  abcd[386] = 2.0E0*I_ERI_Pz_Dyz_Px_Pz_b;
  abcd[387] = 2.0E0*I_ERI_Px_Dxy_Py_Pz_b;
  abcd[388] = 2.0E0*I_ERI_Py_Dxy_Py_Pz_b;
  abcd[389] = 2.0E0*I_ERI_Pz_Dxy_Py_Pz_b;
  abcd[390] = 2.0E0*I_ERI_Px_D2y_Py_Pz_b-1*I_ERI_Px_S_Py_Pz;
  abcd[391] = 2.0E0*I_ERI_Py_D2y_Py_Pz_b-1*I_ERI_Py_S_Py_Pz;
  abcd[392] = 2.0E0*I_ERI_Pz_D2y_Py_Pz_b-1*I_ERI_Pz_S_Py_Pz;
  abcd[393] = 2.0E0*I_ERI_Px_Dyz_Py_Pz_b;
  abcd[394] = 2.0E0*I_ERI_Py_Dyz_Py_Pz_b;
  abcd[395] = 2.0E0*I_ERI_Pz_Dyz_Py_Pz_b;
  abcd[396] = 2.0E0*I_ERI_Px_Dxy_Pz_Pz_b;
  abcd[397] = 2.0E0*I_ERI_Py_Dxy_Pz_Pz_b;
  abcd[398] = 2.0E0*I_ERI_Pz_Dxy_Pz_Pz_b;
  abcd[399] = 2.0E0*I_ERI_Px_D2y_Pz_Pz_b-1*I_ERI_Px_S_Pz_Pz;
  abcd[400] = 2.0E0*I_ERI_Py_D2y_Pz_Pz_b-1*I_ERI_Py_S_Pz_Pz;
  abcd[401] = 2.0E0*I_ERI_Pz_D2y_Pz_Pz_b-1*I_ERI_Pz_S_Pz_Pz;
  abcd[402] = 2.0E0*I_ERI_Px_Dyz_Pz_Pz_b;
  abcd[403] = 2.0E0*I_ERI_Py_Dyz_Pz_Pz_b;
  abcd[404] = 2.0E0*I_ERI_Pz_Dyz_Pz_Pz_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_D_P_P_b
   * RHS shell quartet name: SQ_ERI_P_S_P_P
   ************************************************************/
  abcd[405] = 2.0E0*I_ERI_Px_Dxz_Px_Px_b;
  abcd[406] = 2.0E0*I_ERI_Py_Dxz_Px_Px_b;
  abcd[407] = 2.0E0*I_ERI_Pz_Dxz_Px_Px_b;
  abcd[408] = 2.0E0*I_ERI_Px_Dyz_Px_Px_b;
  abcd[409] = 2.0E0*I_ERI_Py_Dyz_Px_Px_b;
  abcd[410] = 2.0E0*I_ERI_Pz_Dyz_Px_Px_b;
  abcd[411] = 2.0E0*I_ERI_Px_D2z_Px_Px_b-1*I_ERI_Px_S_Px_Px;
  abcd[412] = 2.0E0*I_ERI_Py_D2z_Px_Px_b-1*I_ERI_Py_S_Px_Px;
  abcd[413] = 2.0E0*I_ERI_Pz_D2z_Px_Px_b-1*I_ERI_Pz_S_Px_Px;
  abcd[414] = 2.0E0*I_ERI_Px_Dxz_Py_Px_b;
  abcd[415] = 2.0E0*I_ERI_Py_Dxz_Py_Px_b;
  abcd[416] = 2.0E0*I_ERI_Pz_Dxz_Py_Px_b;
  abcd[417] = 2.0E0*I_ERI_Px_Dyz_Py_Px_b;
  abcd[418] = 2.0E0*I_ERI_Py_Dyz_Py_Px_b;
  abcd[419] = 2.0E0*I_ERI_Pz_Dyz_Py_Px_b;
  abcd[420] = 2.0E0*I_ERI_Px_D2z_Py_Px_b-1*I_ERI_Px_S_Py_Px;
  abcd[421] = 2.0E0*I_ERI_Py_D2z_Py_Px_b-1*I_ERI_Py_S_Py_Px;
  abcd[422] = 2.0E0*I_ERI_Pz_D2z_Py_Px_b-1*I_ERI_Pz_S_Py_Px;
  abcd[423] = 2.0E0*I_ERI_Px_Dxz_Pz_Px_b;
  abcd[424] = 2.0E0*I_ERI_Py_Dxz_Pz_Px_b;
  abcd[425] = 2.0E0*I_ERI_Pz_Dxz_Pz_Px_b;
  abcd[426] = 2.0E0*I_ERI_Px_Dyz_Pz_Px_b;
  abcd[427] = 2.0E0*I_ERI_Py_Dyz_Pz_Px_b;
  abcd[428] = 2.0E0*I_ERI_Pz_Dyz_Pz_Px_b;
  abcd[429] = 2.0E0*I_ERI_Px_D2z_Pz_Px_b-1*I_ERI_Px_S_Pz_Px;
  abcd[430] = 2.0E0*I_ERI_Py_D2z_Pz_Px_b-1*I_ERI_Py_S_Pz_Px;
  abcd[431] = 2.0E0*I_ERI_Pz_D2z_Pz_Px_b-1*I_ERI_Pz_S_Pz_Px;
  abcd[432] = 2.0E0*I_ERI_Px_Dxz_Px_Py_b;
  abcd[433] = 2.0E0*I_ERI_Py_Dxz_Px_Py_b;
  abcd[434] = 2.0E0*I_ERI_Pz_Dxz_Px_Py_b;
  abcd[435] = 2.0E0*I_ERI_Px_Dyz_Px_Py_b;
  abcd[436] = 2.0E0*I_ERI_Py_Dyz_Px_Py_b;
  abcd[437] = 2.0E0*I_ERI_Pz_Dyz_Px_Py_b;
  abcd[438] = 2.0E0*I_ERI_Px_D2z_Px_Py_b-1*I_ERI_Px_S_Px_Py;
  abcd[439] = 2.0E0*I_ERI_Py_D2z_Px_Py_b-1*I_ERI_Py_S_Px_Py;
  abcd[440] = 2.0E0*I_ERI_Pz_D2z_Px_Py_b-1*I_ERI_Pz_S_Px_Py;
  abcd[441] = 2.0E0*I_ERI_Px_Dxz_Py_Py_b;
  abcd[442] = 2.0E0*I_ERI_Py_Dxz_Py_Py_b;
  abcd[443] = 2.0E0*I_ERI_Pz_Dxz_Py_Py_b;
  abcd[444] = 2.0E0*I_ERI_Px_Dyz_Py_Py_b;
  abcd[445] = 2.0E0*I_ERI_Py_Dyz_Py_Py_b;
  abcd[446] = 2.0E0*I_ERI_Pz_Dyz_Py_Py_b;
  abcd[447] = 2.0E0*I_ERI_Px_D2z_Py_Py_b-1*I_ERI_Px_S_Py_Py;
  abcd[448] = 2.0E0*I_ERI_Py_D2z_Py_Py_b-1*I_ERI_Py_S_Py_Py;
  abcd[449] = 2.0E0*I_ERI_Pz_D2z_Py_Py_b-1*I_ERI_Pz_S_Py_Py;
  abcd[450] = 2.0E0*I_ERI_Px_Dxz_Pz_Py_b;
  abcd[451] = 2.0E0*I_ERI_Py_Dxz_Pz_Py_b;
  abcd[452] = 2.0E0*I_ERI_Pz_Dxz_Pz_Py_b;
  abcd[453] = 2.0E0*I_ERI_Px_Dyz_Pz_Py_b;
  abcd[454] = 2.0E0*I_ERI_Py_Dyz_Pz_Py_b;
  abcd[455] = 2.0E0*I_ERI_Pz_Dyz_Pz_Py_b;
  abcd[456] = 2.0E0*I_ERI_Px_D2z_Pz_Py_b-1*I_ERI_Px_S_Pz_Py;
  abcd[457] = 2.0E0*I_ERI_Py_D2z_Pz_Py_b-1*I_ERI_Py_S_Pz_Py;
  abcd[458] = 2.0E0*I_ERI_Pz_D2z_Pz_Py_b-1*I_ERI_Pz_S_Pz_Py;
  abcd[459] = 2.0E0*I_ERI_Px_Dxz_Px_Pz_b;
  abcd[460] = 2.0E0*I_ERI_Py_Dxz_Px_Pz_b;
  abcd[461] = 2.0E0*I_ERI_Pz_Dxz_Px_Pz_b;
  abcd[462] = 2.0E0*I_ERI_Px_Dyz_Px_Pz_b;
  abcd[463] = 2.0E0*I_ERI_Py_Dyz_Px_Pz_b;
  abcd[464] = 2.0E0*I_ERI_Pz_Dyz_Px_Pz_b;
  abcd[465] = 2.0E0*I_ERI_Px_D2z_Px_Pz_b-1*I_ERI_Px_S_Px_Pz;
  abcd[466] = 2.0E0*I_ERI_Py_D2z_Px_Pz_b-1*I_ERI_Py_S_Px_Pz;
  abcd[467] = 2.0E0*I_ERI_Pz_D2z_Px_Pz_b-1*I_ERI_Pz_S_Px_Pz;
  abcd[468] = 2.0E0*I_ERI_Px_Dxz_Py_Pz_b;
  abcd[469] = 2.0E0*I_ERI_Py_Dxz_Py_Pz_b;
  abcd[470] = 2.0E0*I_ERI_Pz_Dxz_Py_Pz_b;
  abcd[471] = 2.0E0*I_ERI_Px_Dyz_Py_Pz_b;
  abcd[472] = 2.0E0*I_ERI_Py_Dyz_Py_Pz_b;
  abcd[473] = 2.0E0*I_ERI_Pz_Dyz_Py_Pz_b;
  abcd[474] = 2.0E0*I_ERI_Px_D2z_Py_Pz_b-1*I_ERI_Px_S_Py_Pz;
  abcd[475] = 2.0E0*I_ERI_Py_D2z_Py_Pz_b-1*I_ERI_Py_S_Py_Pz;
  abcd[476] = 2.0E0*I_ERI_Pz_D2z_Py_Pz_b-1*I_ERI_Pz_S_Py_Pz;
  abcd[477] = 2.0E0*I_ERI_Px_Dxz_Pz_Pz_b;
  abcd[478] = 2.0E0*I_ERI_Py_Dxz_Pz_Pz_b;
  abcd[479] = 2.0E0*I_ERI_Pz_Dxz_Pz_Pz_b;
  abcd[480] = 2.0E0*I_ERI_Px_Dyz_Pz_Pz_b;
  abcd[481] = 2.0E0*I_ERI_Py_Dyz_Pz_Pz_b;
  abcd[482] = 2.0E0*I_ERI_Pz_Dyz_Pz_Pz_b;
  abcd[483] = 2.0E0*I_ERI_Px_D2z_Pz_Pz_b-1*I_ERI_Px_S_Pz_Pz;
  abcd[484] = 2.0E0*I_ERI_Py_D2z_Pz_Pz_b-1*I_ERI_Py_S_Pz_Pz;
  abcd[485] = 2.0E0*I_ERI_Pz_D2z_Pz_Pz_b-1*I_ERI_Pz_S_Pz_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P
   ************************************************************/
  abcd[486] = 2.0E0*I_ERI_Px_Px_D2x_Px_c-1*I_ERI_Px_Px_S_Px;
  abcd[487] = 2.0E0*I_ERI_Py_Px_D2x_Px_c-1*I_ERI_Py_Px_S_Px;
  abcd[488] = 2.0E0*I_ERI_Pz_Px_D2x_Px_c-1*I_ERI_Pz_Px_S_Px;
  abcd[489] = 2.0E0*I_ERI_Px_Py_D2x_Px_c-1*I_ERI_Px_Py_S_Px;
  abcd[490] = 2.0E0*I_ERI_Py_Py_D2x_Px_c-1*I_ERI_Py_Py_S_Px;
  abcd[491] = 2.0E0*I_ERI_Pz_Py_D2x_Px_c-1*I_ERI_Pz_Py_S_Px;
  abcd[492] = 2.0E0*I_ERI_Px_Pz_D2x_Px_c-1*I_ERI_Px_Pz_S_Px;
  abcd[493] = 2.0E0*I_ERI_Py_Pz_D2x_Px_c-1*I_ERI_Py_Pz_S_Px;
  abcd[494] = 2.0E0*I_ERI_Pz_Pz_D2x_Px_c-1*I_ERI_Pz_Pz_S_Px;
  abcd[495] = 2.0E0*I_ERI_Px_Px_Dxy_Px_c;
  abcd[496] = 2.0E0*I_ERI_Py_Px_Dxy_Px_c;
  abcd[497] = 2.0E0*I_ERI_Pz_Px_Dxy_Px_c;
  abcd[498] = 2.0E0*I_ERI_Px_Py_Dxy_Px_c;
  abcd[499] = 2.0E0*I_ERI_Py_Py_Dxy_Px_c;
  abcd[500] = 2.0E0*I_ERI_Pz_Py_Dxy_Px_c;
  abcd[501] = 2.0E0*I_ERI_Px_Pz_Dxy_Px_c;
  abcd[502] = 2.0E0*I_ERI_Py_Pz_Dxy_Px_c;
  abcd[503] = 2.0E0*I_ERI_Pz_Pz_Dxy_Px_c;
  abcd[504] = 2.0E0*I_ERI_Px_Px_Dxz_Px_c;
  abcd[505] = 2.0E0*I_ERI_Py_Px_Dxz_Px_c;
  abcd[506] = 2.0E0*I_ERI_Pz_Px_Dxz_Px_c;
  abcd[507] = 2.0E0*I_ERI_Px_Py_Dxz_Px_c;
  abcd[508] = 2.0E0*I_ERI_Py_Py_Dxz_Px_c;
  abcd[509] = 2.0E0*I_ERI_Pz_Py_Dxz_Px_c;
  abcd[510] = 2.0E0*I_ERI_Px_Pz_Dxz_Px_c;
  abcd[511] = 2.0E0*I_ERI_Py_Pz_Dxz_Px_c;
  abcd[512] = 2.0E0*I_ERI_Pz_Pz_Dxz_Px_c;
  abcd[513] = 2.0E0*I_ERI_Px_Px_D2x_Py_c-1*I_ERI_Px_Px_S_Py;
  abcd[514] = 2.0E0*I_ERI_Py_Px_D2x_Py_c-1*I_ERI_Py_Px_S_Py;
  abcd[515] = 2.0E0*I_ERI_Pz_Px_D2x_Py_c-1*I_ERI_Pz_Px_S_Py;
  abcd[516] = 2.0E0*I_ERI_Px_Py_D2x_Py_c-1*I_ERI_Px_Py_S_Py;
  abcd[517] = 2.0E0*I_ERI_Py_Py_D2x_Py_c-1*I_ERI_Py_Py_S_Py;
  abcd[518] = 2.0E0*I_ERI_Pz_Py_D2x_Py_c-1*I_ERI_Pz_Py_S_Py;
  abcd[519] = 2.0E0*I_ERI_Px_Pz_D2x_Py_c-1*I_ERI_Px_Pz_S_Py;
  abcd[520] = 2.0E0*I_ERI_Py_Pz_D2x_Py_c-1*I_ERI_Py_Pz_S_Py;
  abcd[521] = 2.0E0*I_ERI_Pz_Pz_D2x_Py_c-1*I_ERI_Pz_Pz_S_Py;
  abcd[522] = 2.0E0*I_ERI_Px_Px_Dxy_Py_c;
  abcd[523] = 2.0E0*I_ERI_Py_Px_Dxy_Py_c;
  abcd[524] = 2.0E0*I_ERI_Pz_Px_Dxy_Py_c;
  abcd[525] = 2.0E0*I_ERI_Px_Py_Dxy_Py_c;
  abcd[526] = 2.0E0*I_ERI_Py_Py_Dxy_Py_c;
  abcd[527] = 2.0E0*I_ERI_Pz_Py_Dxy_Py_c;
  abcd[528] = 2.0E0*I_ERI_Px_Pz_Dxy_Py_c;
  abcd[529] = 2.0E0*I_ERI_Py_Pz_Dxy_Py_c;
  abcd[530] = 2.0E0*I_ERI_Pz_Pz_Dxy_Py_c;
  abcd[531] = 2.0E0*I_ERI_Px_Px_Dxz_Py_c;
  abcd[532] = 2.0E0*I_ERI_Py_Px_Dxz_Py_c;
  abcd[533] = 2.0E0*I_ERI_Pz_Px_Dxz_Py_c;
  abcd[534] = 2.0E0*I_ERI_Px_Py_Dxz_Py_c;
  abcd[535] = 2.0E0*I_ERI_Py_Py_Dxz_Py_c;
  abcd[536] = 2.0E0*I_ERI_Pz_Py_Dxz_Py_c;
  abcd[537] = 2.0E0*I_ERI_Px_Pz_Dxz_Py_c;
  abcd[538] = 2.0E0*I_ERI_Py_Pz_Dxz_Py_c;
  abcd[539] = 2.0E0*I_ERI_Pz_Pz_Dxz_Py_c;
  abcd[540] = 2.0E0*I_ERI_Px_Px_D2x_Pz_c-1*I_ERI_Px_Px_S_Pz;
  abcd[541] = 2.0E0*I_ERI_Py_Px_D2x_Pz_c-1*I_ERI_Py_Px_S_Pz;
  abcd[542] = 2.0E0*I_ERI_Pz_Px_D2x_Pz_c-1*I_ERI_Pz_Px_S_Pz;
  abcd[543] = 2.0E0*I_ERI_Px_Py_D2x_Pz_c-1*I_ERI_Px_Py_S_Pz;
  abcd[544] = 2.0E0*I_ERI_Py_Py_D2x_Pz_c-1*I_ERI_Py_Py_S_Pz;
  abcd[545] = 2.0E0*I_ERI_Pz_Py_D2x_Pz_c-1*I_ERI_Pz_Py_S_Pz;
  abcd[546] = 2.0E0*I_ERI_Px_Pz_D2x_Pz_c-1*I_ERI_Px_Pz_S_Pz;
  abcd[547] = 2.0E0*I_ERI_Py_Pz_D2x_Pz_c-1*I_ERI_Py_Pz_S_Pz;
  abcd[548] = 2.0E0*I_ERI_Pz_Pz_D2x_Pz_c-1*I_ERI_Pz_Pz_S_Pz;
  abcd[549] = 2.0E0*I_ERI_Px_Px_Dxy_Pz_c;
  abcd[550] = 2.0E0*I_ERI_Py_Px_Dxy_Pz_c;
  abcd[551] = 2.0E0*I_ERI_Pz_Px_Dxy_Pz_c;
  abcd[552] = 2.0E0*I_ERI_Px_Py_Dxy_Pz_c;
  abcd[553] = 2.0E0*I_ERI_Py_Py_Dxy_Pz_c;
  abcd[554] = 2.0E0*I_ERI_Pz_Py_Dxy_Pz_c;
  abcd[555] = 2.0E0*I_ERI_Px_Pz_Dxy_Pz_c;
  abcd[556] = 2.0E0*I_ERI_Py_Pz_Dxy_Pz_c;
  abcd[557] = 2.0E0*I_ERI_Pz_Pz_Dxy_Pz_c;
  abcd[558] = 2.0E0*I_ERI_Px_Px_Dxz_Pz_c;
  abcd[559] = 2.0E0*I_ERI_Py_Px_Dxz_Pz_c;
  abcd[560] = 2.0E0*I_ERI_Pz_Px_Dxz_Pz_c;
  abcd[561] = 2.0E0*I_ERI_Px_Py_Dxz_Pz_c;
  abcd[562] = 2.0E0*I_ERI_Py_Py_Dxz_Pz_c;
  abcd[563] = 2.0E0*I_ERI_Pz_Py_Dxz_Pz_c;
  abcd[564] = 2.0E0*I_ERI_Px_Pz_Dxz_Pz_c;
  abcd[565] = 2.0E0*I_ERI_Py_Pz_Dxz_Pz_c;
  abcd[566] = 2.0E0*I_ERI_Pz_Pz_Dxz_Pz_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P
   ************************************************************/
  abcd[567] = 2.0E0*I_ERI_Px_Px_Dxy_Px_c;
  abcd[568] = 2.0E0*I_ERI_Py_Px_Dxy_Px_c;
  abcd[569] = 2.0E0*I_ERI_Pz_Px_Dxy_Px_c;
  abcd[570] = 2.0E0*I_ERI_Px_Py_Dxy_Px_c;
  abcd[571] = 2.0E0*I_ERI_Py_Py_Dxy_Px_c;
  abcd[572] = 2.0E0*I_ERI_Pz_Py_Dxy_Px_c;
  abcd[573] = 2.0E0*I_ERI_Px_Pz_Dxy_Px_c;
  abcd[574] = 2.0E0*I_ERI_Py_Pz_Dxy_Px_c;
  abcd[575] = 2.0E0*I_ERI_Pz_Pz_Dxy_Px_c;
  abcd[576] = 2.0E0*I_ERI_Px_Px_D2y_Px_c-1*I_ERI_Px_Px_S_Px;
  abcd[577] = 2.0E0*I_ERI_Py_Px_D2y_Px_c-1*I_ERI_Py_Px_S_Px;
  abcd[578] = 2.0E0*I_ERI_Pz_Px_D2y_Px_c-1*I_ERI_Pz_Px_S_Px;
  abcd[579] = 2.0E0*I_ERI_Px_Py_D2y_Px_c-1*I_ERI_Px_Py_S_Px;
  abcd[580] = 2.0E0*I_ERI_Py_Py_D2y_Px_c-1*I_ERI_Py_Py_S_Px;
  abcd[581] = 2.0E0*I_ERI_Pz_Py_D2y_Px_c-1*I_ERI_Pz_Py_S_Px;
  abcd[582] = 2.0E0*I_ERI_Px_Pz_D2y_Px_c-1*I_ERI_Px_Pz_S_Px;
  abcd[583] = 2.0E0*I_ERI_Py_Pz_D2y_Px_c-1*I_ERI_Py_Pz_S_Px;
  abcd[584] = 2.0E0*I_ERI_Pz_Pz_D2y_Px_c-1*I_ERI_Pz_Pz_S_Px;
  abcd[585] = 2.0E0*I_ERI_Px_Px_Dyz_Px_c;
  abcd[586] = 2.0E0*I_ERI_Py_Px_Dyz_Px_c;
  abcd[587] = 2.0E0*I_ERI_Pz_Px_Dyz_Px_c;
  abcd[588] = 2.0E0*I_ERI_Px_Py_Dyz_Px_c;
  abcd[589] = 2.0E0*I_ERI_Py_Py_Dyz_Px_c;
  abcd[590] = 2.0E0*I_ERI_Pz_Py_Dyz_Px_c;
  abcd[591] = 2.0E0*I_ERI_Px_Pz_Dyz_Px_c;
  abcd[592] = 2.0E0*I_ERI_Py_Pz_Dyz_Px_c;
  abcd[593] = 2.0E0*I_ERI_Pz_Pz_Dyz_Px_c;
  abcd[594] = 2.0E0*I_ERI_Px_Px_Dxy_Py_c;
  abcd[595] = 2.0E0*I_ERI_Py_Px_Dxy_Py_c;
  abcd[596] = 2.0E0*I_ERI_Pz_Px_Dxy_Py_c;
  abcd[597] = 2.0E0*I_ERI_Px_Py_Dxy_Py_c;
  abcd[598] = 2.0E0*I_ERI_Py_Py_Dxy_Py_c;
  abcd[599] = 2.0E0*I_ERI_Pz_Py_Dxy_Py_c;
  abcd[600] = 2.0E0*I_ERI_Px_Pz_Dxy_Py_c;
  abcd[601] = 2.0E0*I_ERI_Py_Pz_Dxy_Py_c;
  abcd[602] = 2.0E0*I_ERI_Pz_Pz_Dxy_Py_c;
  abcd[603] = 2.0E0*I_ERI_Px_Px_D2y_Py_c-1*I_ERI_Px_Px_S_Py;
  abcd[604] = 2.0E0*I_ERI_Py_Px_D2y_Py_c-1*I_ERI_Py_Px_S_Py;
  abcd[605] = 2.0E0*I_ERI_Pz_Px_D2y_Py_c-1*I_ERI_Pz_Px_S_Py;
  abcd[606] = 2.0E0*I_ERI_Px_Py_D2y_Py_c-1*I_ERI_Px_Py_S_Py;
  abcd[607] = 2.0E0*I_ERI_Py_Py_D2y_Py_c-1*I_ERI_Py_Py_S_Py;
  abcd[608] = 2.0E0*I_ERI_Pz_Py_D2y_Py_c-1*I_ERI_Pz_Py_S_Py;
  abcd[609] = 2.0E0*I_ERI_Px_Pz_D2y_Py_c-1*I_ERI_Px_Pz_S_Py;
  abcd[610] = 2.0E0*I_ERI_Py_Pz_D2y_Py_c-1*I_ERI_Py_Pz_S_Py;
  abcd[611] = 2.0E0*I_ERI_Pz_Pz_D2y_Py_c-1*I_ERI_Pz_Pz_S_Py;
  abcd[612] = 2.0E0*I_ERI_Px_Px_Dyz_Py_c;
  abcd[613] = 2.0E0*I_ERI_Py_Px_Dyz_Py_c;
  abcd[614] = 2.0E0*I_ERI_Pz_Px_Dyz_Py_c;
  abcd[615] = 2.0E0*I_ERI_Px_Py_Dyz_Py_c;
  abcd[616] = 2.0E0*I_ERI_Py_Py_Dyz_Py_c;
  abcd[617] = 2.0E0*I_ERI_Pz_Py_Dyz_Py_c;
  abcd[618] = 2.0E0*I_ERI_Px_Pz_Dyz_Py_c;
  abcd[619] = 2.0E0*I_ERI_Py_Pz_Dyz_Py_c;
  abcd[620] = 2.0E0*I_ERI_Pz_Pz_Dyz_Py_c;
  abcd[621] = 2.0E0*I_ERI_Px_Px_Dxy_Pz_c;
  abcd[622] = 2.0E0*I_ERI_Py_Px_Dxy_Pz_c;
  abcd[623] = 2.0E0*I_ERI_Pz_Px_Dxy_Pz_c;
  abcd[624] = 2.0E0*I_ERI_Px_Py_Dxy_Pz_c;
  abcd[625] = 2.0E0*I_ERI_Py_Py_Dxy_Pz_c;
  abcd[626] = 2.0E0*I_ERI_Pz_Py_Dxy_Pz_c;
  abcd[627] = 2.0E0*I_ERI_Px_Pz_Dxy_Pz_c;
  abcd[628] = 2.0E0*I_ERI_Py_Pz_Dxy_Pz_c;
  abcd[629] = 2.0E0*I_ERI_Pz_Pz_Dxy_Pz_c;
  abcd[630] = 2.0E0*I_ERI_Px_Px_D2y_Pz_c-1*I_ERI_Px_Px_S_Pz;
  abcd[631] = 2.0E0*I_ERI_Py_Px_D2y_Pz_c-1*I_ERI_Py_Px_S_Pz;
  abcd[632] = 2.0E0*I_ERI_Pz_Px_D2y_Pz_c-1*I_ERI_Pz_Px_S_Pz;
  abcd[633] = 2.0E0*I_ERI_Px_Py_D2y_Pz_c-1*I_ERI_Px_Py_S_Pz;
  abcd[634] = 2.0E0*I_ERI_Py_Py_D2y_Pz_c-1*I_ERI_Py_Py_S_Pz;
  abcd[635] = 2.0E0*I_ERI_Pz_Py_D2y_Pz_c-1*I_ERI_Pz_Py_S_Pz;
  abcd[636] = 2.0E0*I_ERI_Px_Pz_D2y_Pz_c-1*I_ERI_Px_Pz_S_Pz;
  abcd[637] = 2.0E0*I_ERI_Py_Pz_D2y_Pz_c-1*I_ERI_Py_Pz_S_Pz;
  abcd[638] = 2.0E0*I_ERI_Pz_Pz_D2y_Pz_c-1*I_ERI_Pz_Pz_S_Pz;
  abcd[639] = 2.0E0*I_ERI_Px_Px_Dyz_Pz_c;
  abcd[640] = 2.0E0*I_ERI_Py_Px_Dyz_Pz_c;
  abcd[641] = 2.0E0*I_ERI_Pz_Px_Dyz_Pz_c;
  abcd[642] = 2.0E0*I_ERI_Px_Py_Dyz_Pz_c;
  abcd[643] = 2.0E0*I_ERI_Py_Py_Dyz_Pz_c;
  abcd[644] = 2.0E0*I_ERI_Pz_Py_Dyz_Pz_c;
  abcd[645] = 2.0E0*I_ERI_Px_Pz_Dyz_Pz_c;
  abcd[646] = 2.0E0*I_ERI_Py_Pz_Dyz_Pz_c;
  abcd[647] = 2.0E0*I_ERI_Pz_Pz_Dyz_Pz_c;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_P_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_P_P_D_P_c
   * RHS shell quartet name: SQ_ERI_P_P_S_P
   ************************************************************/
  abcd[648] = 2.0E0*I_ERI_Px_Px_Dxz_Px_c;
  abcd[649] = 2.0E0*I_ERI_Py_Px_Dxz_Px_c;
  abcd[650] = 2.0E0*I_ERI_Pz_Px_Dxz_Px_c;
  abcd[651] = 2.0E0*I_ERI_Px_Py_Dxz_Px_c;
  abcd[652] = 2.0E0*I_ERI_Py_Py_Dxz_Px_c;
  abcd[653] = 2.0E0*I_ERI_Pz_Py_Dxz_Px_c;
  abcd[654] = 2.0E0*I_ERI_Px_Pz_Dxz_Px_c;
  abcd[655] = 2.0E0*I_ERI_Py_Pz_Dxz_Px_c;
  abcd[656] = 2.0E0*I_ERI_Pz_Pz_Dxz_Px_c;
  abcd[657] = 2.0E0*I_ERI_Px_Px_Dyz_Px_c;
  abcd[658] = 2.0E0*I_ERI_Py_Px_Dyz_Px_c;
  abcd[659] = 2.0E0*I_ERI_Pz_Px_Dyz_Px_c;
  abcd[660] = 2.0E0*I_ERI_Px_Py_Dyz_Px_c;
  abcd[661] = 2.0E0*I_ERI_Py_Py_Dyz_Px_c;
  abcd[662] = 2.0E0*I_ERI_Pz_Py_Dyz_Px_c;
  abcd[663] = 2.0E0*I_ERI_Px_Pz_Dyz_Px_c;
  abcd[664] = 2.0E0*I_ERI_Py_Pz_Dyz_Px_c;
  abcd[665] = 2.0E0*I_ERI_Pz_Pz_Dyz_Px_c;
  abcd[666] = 2.0E0*I_ERI_Px_Px_D2z_Px_c-1*I_ERI_Px_Px_S_Px;
  abcd[667] = 2.0E0*I_ERI_Py_Px_D2z_Px_c-1*I_ERI_Py_Px_S_Px;
  abcd[668] = 2.0E0*I_ERI_Pz_Px_D2z_Px_c-1*I_ERI_Pz_Px_S_Px;
  abcd[669] = 2.0E0*I_ERI_Px_Py_D2z_Px_c-1*I_ERI_Px_Py_S_Px;
  abcd[670] = 2.0E0*I_ERI_Py_Py_D2z_Px_c-1*I_ERI_Py_Py_S_Px;
  abcd[671] = 2.0E0*I_ERI_Pz_Py_D2z_Px_c-1*I_ERI_Pz_Py_S_Px;
  abcd[672] = 2.0E0*I_ERI_Px_Pz_D2z_Px_c-1*I_ERI_Px_Pz_S_Px;
  abcd[673] = 2.0E0*I_ERI_Py_Pz_D2z_Px_c-1*I_ERI_Py_Pz_S_Px;
  abcd[674] = 2.0E0*I_ERI_Pz_Pz_D2z_Px_c-1*I_ERI_Pz_Pz_S_Px;
  abcd[675] = 2.0E0*I_ERI_Px_Px_Dxz_Py_c;
  abcd[676] = 2.0E0*I_ERI_Py_Px_Dxz_Py_c;
  abcd[677] = 2.0E0*I_ERI_Pz_Px_Dxz_Py_c;
  abcd[678] = 2.0E0*I_ERI_Px_Py_Dxz_Py_c;
  abcd[679] = 2.0E0*I_ERI_Py_Py_Dxz_Py_c;
  abcd[680] = 2.0E0*I_ERI_Pz_Py_Dxz_Py_c;
  abcd[681] = 2.0E0*I_ERI_Px_Pz_Dxz_Py_c;
  abcd[682] = 2.0E0*I_ERI_Py_Pz_Dxz_Py_c;
  abcd[683] = 2.0E0*I_ERI_Pz_Pz_Dxz_Py_c;
  abcd[684] = 2.0E0*I_ERI_Px_Px_Dyz_Py_c;
  abcd[685] = 2.0E0*I_ERI_Py_Px_Dyz_Py_c;
  abcd[686] = 2.0E0*I_ERI_Pz_Px_Dyz_Py_c;
  abcd[687] = 2.0E0*I_ERI_Px_Py_Dyz_Py_c;
  abcd[688] = 2.0E0*I_ERI_Py_Py_Dyz_Py_c;
  abcd[689] = 2.0E0*I_ERI_Pz_Py_Dyz_Py_c;
  abcd[690] = 2.0E0*I_ERI_Px_Pz_Dyz_Py_c;
  abcd[691] = 2.0E0*I_ERI_Py_Pz_Dyz_Py_c;
  abcd[692] = 2.0E0*I_ERI_Pz_Pz_Dyz_Py_c;
  abcd[693] = 2.0E0*I_ERI_Px_Px_D2z_Py_c-1*I_ERI_Px_Px_S_Py;
  abcd[694] = 2.0E0*I_ERI_Py_Px_D2z_Py_c-1*I_ERI_Py_Px_S_Py;
  abcd[695] = 2.0E0*I_ERI_Pz_Px_D2z_Py_c-1*I_ERI_Pz_Px_S_Py;
  abcd[696] = 2.0E0*I_ERI_Px_Py_D2z_Py_c-1*I_ERI_Px_Py_S_Py;
  abcd[697] = 2.0E0*I_ERI_Py_Py_D2z_Py_c-1*I_ERI_Py_Py_S_Py;
  abcd[698] = 2.0E0*I_ERI_Pz_Py_D2z_Py_c-1*I_ERI_Pz_Py_S_Py;
  abcd[699] = 2.0E0*I_ERI_Px_Pz_D2z_Py_c-1*I_ERI_Px_Pz_S_Py;
  abcd[700] = 2.0E0*I_ERI_Py_Pz_D2z_Py_c-1*I_ERI_Py_Pz_S_Py;
  abcd[701] = 2.0E0*I_ERI_Pz_Pz_D2z_Py_c-1*I_ERI_Pz_Pz_S_Py;
  abcd[702] = 2.0E0*I_ERI_Px_Px_Dxz_Pz_c;
  abcd[703] = 2.0E0*I_ERI_Py_Px_Dxz_Pz_c;
  abcd[704] = 2.0E0*I_ERI_Pz_Px_Dxz_Pz_c;
  abcd[705] = 2.0E0*I_ERI_Px_Py_Dxz_Pz_c;
  abcd[706] = 2.0E0*I_ERI_Py_Py_Dxz_Pz_c;
  abcd[707] = 2.0E0*I_ERI_Pz_Py_Dxz_Pz_c;
  abcd[708] = 2.0E0*I_ERI_Px_Pz_Dxz_Pz_c;
  abcd[709] = 2.0E0*I_ERI_Py_Pz_Dxz_Pz_c;
  abcd[710] = 2.0E0*I_ERI_Pz_Pz_Dxz_Pz_c;
  abcd[711] = 2.0E0*I_ERI_Px_Px_Dyz_Pz_c;
  abcd[712] = 2.0E0*I_ERI_Py_Px_Dyz_Pz_c;
  abcd[713] = 2.0E0*I_ERI_Pz_Px_Dyz_Pz_c;
  abcd[714] = 2.0E0*I_ERI_Px_Py_Dyz_Pz_c;
  abcd[715] = 2.0E0*I_ERI_Py_Py_Dyz_Pz_c;
  abcd[716] = 2.0E0*I_ERI_Pz_Py_Dyz_Pz_c;
  abcd[717] = 2.0E0*I_ERI_Px_Pz_Dyz_Pz_c;
  abcd[718] = 2.0E0*I_ERI_Py_Pz_Dyz_Pz_c;
  abcd[719] = 2.0E0*I_ERI_Pz_Pz_Dyz_Pz_c;
  abcd[720] = 2.0E0*I_ERI_Px_Px_D2z_Pz_c-1*I_ERI_Px_Px_S_Pz;
  abcd[721] = 2.0E0*I_ERI_Py_Px_D2z_Pz_c-1*I_ERI_Py_Px_S_Pz;
  abcd[722] = 2.0E0*I_ERI_Pz_Px_D2z_Pz_c-1*I_ERI_Pz_Px_S_Pz;
  abcd[723] = 2.0E0*I_ERI_Px_Py_D2z_Pz_c-1*I_ERI_Px_Py_S_Pz;
  abcd[724] = 2.0E0*I_ERI_Py_Py_D2z_Pz_c-1*I_ERI_Py_Py_S_Pz;
  abcd[725] = 2.0E0*I_ERI_Pz_Py_D2z_Pz_c-1*I_ERI_Pz_Py_S_Pz;
  abcd[726] = 2.0E0*I_ERI_Px_Pz_D2z_Pz_c-1*I_ERI_Px_Pz_S_Pz;
  abcd[727] = 2.0E0*I_ERI_Py_Pz_D2z_Pz_c-1*I_ERI_Py_Pz_S_Pz;
  abcd[728] = 2.0E0*I_ERI_Pz_Pz_D2z_Pz_c-1*I_ERI_Pz_Pz_S_Pz;
}
