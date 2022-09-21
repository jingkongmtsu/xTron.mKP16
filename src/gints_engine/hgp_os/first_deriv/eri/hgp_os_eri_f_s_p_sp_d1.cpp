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
// BRA1 as redundant position, total RHS integrals evaluated as: 59521
// BRA2 as redundant position, total RHS integrals evaluated as: 61483
// KET1 as redundant position, total RHS integrals evaluated as: 68436
// KET2 as redundant position, total RHS integrals evaluated as: 61563
// the redundant position is: BRA1
//

//
// @@@@ derivative position-direction information
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
// KET2
// X
// Y
// Z
// ####

void hgp_os_eri_f_s_p_sp_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1000003_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1000003 = 0.0E0;
  Double I_ERI_F3x_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_S_Px_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_S_Py_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_S_Pz_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000003_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000003_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000003_b = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000003_d = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000003_d = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001000003_b = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001000003_d = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001000003_d = 0.0E0;

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
      Double jc2_1 = jcoe[jp2+1*jnp2];
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
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 6 integrals are omitted 
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
      Double I_ERI_D2y_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2y_S_Px_S_M1_vrr+WQZ*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2z_S_Px_S_M1_vrr+WQZ*I_ERI_D2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_D2y_S_M1_vrr = QCY*I_ERI_D2x_S_Py_S_M1_vrr+WQY*I_ERI_D2x_S_Py_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = QCY*I_ERI_Dxy_S_Py_S_M1_vrr+WQY*I_ERI_Dxy_S_Py_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = QCY*I_ERI_Dxz_S_Py_S_M1_vrr+WQY*I_ERI_Dxz_S_Py_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = QCY*I_ERI_D2y_S_Py_S_M1_vrr+WQY*I_ERI_D2y_S_Py_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = QCY*I_ERI_Dyz_S_Py_S_M1_vrr+WQY*I_ERI_Dyz_S_Py_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = QCY*I_ERI_D2z_S_Py_S_M1_vrr+WQY*I_ERI_D2z_S_Py_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Dyz_S_M1_vrr = QCZ*I_ERI_D2x_S_Py_S_M1_vrr+WQZ*I_ERI_D2x_S_Py_S_M2_vrr;
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
       * totally 0 integrals are omitted 
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
      Double I_ERI_Fxyz_S_Px_S_M1_vrr = PAZ*I_ERI_Dxy_S_Px_S_M1_vrr+WPZ*I_ERI_Dxy_S_Px_S_M2_vrr;
      Double I_ERI_Fx2z_S_Px_S_M1_vrr = PAX*I_ERI_D2z_S_Px_S_M1_vrr+WPX*I_ERI_D2z_S_Px_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_Px_S_M1_vrr = PAY*I_ERI_D2y_S_Px_S_M1_vrr+WPY*I_ERI_D2y_S_Px_S_M2_vrr+2*oned2z*I_ERI_Py_S_Px_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Px_S_M1_vrr = PAZ*I_ERI_D2y_S_Px_S_M1_vrr+WPZ*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Fy2z_S_Px_S_M1_vrr = PAY*I_ERI_D2z_S_Px_S_M1_vrr+WPY*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Px_S_M1_vrr = PAZ*I_ERI_D2z_S_Px_S_M1_vrr+WPZ*I_ERI_D2z_S_Px_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Px_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_Py_S_M1_vrr = PAX*I_ERI_D2x_S_Py_S_M1_vrr+WPX*I_ERI_D2x_S_Py_S_M2_vrr+2*oned2z*I_ERI_Px_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_Py_S_M1_vrr = PAY*I_ERI_D2x_S_Py_S_M1_vrr+WPY*I_ERI_D2x_S_Py_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_Py_S_M1_vrr = PAZ*I_ERI_D2x_S_Py_S_M1_vrr+WPZ*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_Py_S_M1_vrr = PAX*I_ERI_D2y_S_Py_S_M1_vrr+WPX*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Fxyz_S_Py_S_M1_vrr = PAZ*I_ERI_Dxy_S_Py_S_M1_vrr+WPZ*I_ERI_Dxy_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_Py_S_M1_vrr = PAX*I_ERI_D2z_S_Py_S_M1_vrr+WPX*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_Py_S_M1_vrr = PAY*I_ERI_D2y_S_Py_S_M1_vrr+WPY*I_ERI_D2y_S_Py_S_M2_vrr+2*oned2z*I_ERI_Py_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_Py_S_M1_vrr = PAZ*I_ERI_D2y_S_Py_S_M1_vrr+WPZ*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Fy2z_S_Py_S_M1_vrr = PAY*I_ERI_D2z_S_Py_S_M1_vrr+WPY*I_ERI_D2z_S_Py_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_Py_S_M1_vrr = PAZ*I_ERI_D2z_S_Py_S_M1_vrr+WPZ*I_ERI_D2z_S_Py_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Py_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_F3x_S_Pz_S_M1_vrr = PAX*I_ERI_D2x_S_Pz_S_M1_vrr+WPX*I_ERI_D2x_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_F2xy_S_Pz_S_M1_vrr = PAY*I_ERI_D2x_S_Pz_S_M1_vrr+WPY*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xz_S_Pz_S_M1_vrr = PAZ*I_ERI_D2x_S_Pz_S_M1_vrr+WPZ*I_ERI_D2x_S_Pz_S_M2_vrr+oned2k*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M1_vrr = PAX*I_ERI_D2y_S_Pz_S_M1_vrr+WPX*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fxyz_S_Pz_S_M1_vrr = PAZ*I_ERI_Dxy_S_Pz_S_M1_vrr+WPZ*I_ERI_Dxy_S_Pz_S_M2_vrr+oned2k*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M1_vrr = PAX*I_ERI_D2z_S_Pz_S_M1_vrr+WPX*I_ERI_D2z_S_Pz_S_M2_vrr;
      Double I_ERI_F3y_S_Pz_S_M1_vrr = PAY*I_ERI_D2y_S_Pz_S_M1_vrr+WPY*I_ERI_D2y_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_F2yz_S_Pz_S_M1_vrr = PAZ*I_ERI_D2y_S_Pz_S_M1_vrr+WPZ*I_ERI_D2y_S_Pz_S_M2_vrr+oned2k*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fy2z_S_Pz_S_M1_vrr = PAY*I_ERI_D2z_S_Pz_S_M1_vrr+WPY*I_ERI_D2z_S_Pz_S_M2_vrr;
      Double I_ERI_F3z_S_Pz_S_M1_vrr = PAZ*I_ERI_D2z_S_Pz_S_M1_vrr+WPZ*I_ERI_D2z_S_Pz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M2_vrr+oned2k*I_ERI_D2z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 4 integrals are omitted 
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
      Double I_ERI_F3x_S_Dxz_S_M1_vrr = QCZ*I_ERI_F3x_S_Px_S_M1_vrr+WQZ*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M1_vrr = QCZ*I_ERI_F2xy_S_Px_S_M1_vrr+WQZ*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxz_S_M1_vrr = QCZ*I_ERI_F2xz_S_Px_S_M1_vrr+WQZ*I_ERI_F2xz_S_Px_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_M1_vrr = QCZ*I_ERI_Fx2y_S_Px_S_M1_vrr+WQZ*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_Fx2z_S_Px_S_M1_vrr+WQZ*I_ERI_Fx2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_Px_S_M2_vrr;
      Double I_ERI_F3y_S_Dxz_S_M1_vrr = QCZ*I_ERI_F3y_S_Px_S_M1_vrr+WQZ*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxz_S_M1_vrr = QCZ*I_ERI_F2yz_S_Px_S_M1_vrr+WQZ*I_ERI_F2yz_S_Px_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Dxz_S_M1_vrr = QCZ*I_ERI_F3z_S_Px_S_M1_vrr+WQZ*I_ERI_F3z_S_Px_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
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
      Double I_ERI_F3x_S_Dyz_S_M1_vrr = QCZ*I_ERI_F3x_S_Py_S_M1_vrr+WQZ*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M1_vrr = QCZ*I_ERI_F2xy_S_Py_S_M1_vrr+WQZ*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_F2xz_S_Dyz_S_M1_vrr = QCZ*I_ERI_F2xz_S_Py_S_M1_vrr+WQZ*I_ERI_F2xz_S_Py_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_M1_vrr = QCZ*I_ERI_Fx2y_S_Py_S_M1_vrr+WQZ*I_ERI_Fx2y_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_M1_vrr = QCZ*I_ERI_Fx2z_S_Py_S_M1_vrr+WQZ*I_ERI_Fx2z_S_Py_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_Dyz_S_M1_vrr = QCZ*I_ERI_F3y_S_Py_S_M1_vrr+WQZ*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_F2yz_S_Dyz_S_M1_vrr = QCZ*I_ERI_F2yz_S_Py_S_M1_vrr+WQZ*I_ERI_F2yz_S_Py_S_M2_vrr+oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_Dyz_S_M1_vrr = QCZ*I_ERI_F3z_S_Py_S_M1_vrr+WQZ*I_ERI_F3z_S_Py_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
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
       * shell quartet name: SQ_ERI_D_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_vrr = QCX*I_ERI_D2x_S_S_S_vrr+WQX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Px_S_vrr = QCX*I_ERI_D2y_S_S_S_vrr+WQX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Px_S_vrr = QCX*I_ERI_D2z_S_S_S_vrr+WQX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Py_S_vrr = QCY*I_ERI_D2x_S_S_S_vrr+WQY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Py_S_vrr = QCY*I_ERI_D2y_S_S_S_vrr+WQY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_Py_S_vrr = QCY*I_ERI_D2z_S_S_S_vrr+WQY*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Pz_S_vrr = QCZ*I_ERI_D2x_S_S_S_vrr+WQZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_Pz_S_vrr = QCZ*I_ERI_D2y_S_S_S_vrr+WQZ*I_ERI_D2y_S_S_S_M1_vrr;
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
       * totally 18 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_vrr = QCX*I_ERI_D2x_S_Px_S_vrr+WQX*I_ERI_D2x_S_Px_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr+2*oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_D2x_S_vrr = QCX*I_ERI_D2y_S_Px_S_vrr+WQX*I_ERI_D2y_S_Px_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_D2x_S_vrr = QCX*I_ERI_D2z_S_Px_S_vrr+WQX*I_ERI_D2z_S_Px_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dxy_S_vrr = QCY*I_ERI_D2x_S_Px_S_vrr+WQY*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxy_S_vrr = QCY*I_ERI_D2y_S_Px_S_vrr+WQY*I_ERI_D2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxy_S_vrr = QCY*I_ERI_D2z_S_Px_S_vrr+WQY*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_Dxz_S_vrr = QCZ*I_ERI_D2x_S_Px_S_vrr+WQZ*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxz_S_vrr = QCZ*I_ERI_D2y_S_Px_S_vrr+WQZ*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxz_S_vrr = QCZ*I_ERI_D2z_S_Px_S_vrr+WQZ*I_ERI_D2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_D2y_S_vrr = QCY*I_ERI_D2x_S_Py_S_vrr+WQY*I_ERI_D2x_S_Py_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = QCY*I_ERI_D2y_S_Py_S_vrr+WQY*I_ERI_D2y_S_Py_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = QCY*I_ERI_D2z_S_Py_S_vrr+WQY*I_ERI_D2z_S_Py_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_D2x_S_Dyz_S_vrr = QCZ*I_ERI_D2x_S_Py_S_vrr+WQZ*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Dyz_S_vrr = QCZ*I_ERI_D2y_S_Py_S_vrr+WQZ*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Dyz_S_vrr = QCZ*I_ERI_D2z_S_Py_S_vrr+WQZ*I_ERI_D2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = QCZ*I_ERI_D2x_S_Pz_S_vrr+WQZ*I_ERI_D2x_S_Pz_S_M1_vrr+oned2e*I_ERI_D2x_S_S_S_vrr-rhod2esq*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = QCZ*I_ERI_D2y_S_Pz_S_vrr+WQZ*I_ERI_D2y_S_Pz_S_M1_vrr+oned2e*I_ERI_D2y_S_S_S_vrr-rhod2esq*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = QCZ*I_ERI_D2z_S_Pz_S_vrr+WQZ*I_ERI_D2z_S_Pz_S_M1_vrr+oned2e*I_ERI_D2z_S_S_S_vrr-rhod2esq*I_ERI_D2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

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
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_D_S
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_vrr = PAX*I_ERI_F3x_S_D2x_S_vrr+WPX*I_ERI_F3x_S_D2x_S_M1_vrr+3*oned2z*I_ERI_D2x_S_D2x_S_vrr-3*rhod2zsq*I_ERI_D2x_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G3xy_S_D2x_S_vrr = PAY*I_ERI_F3x_S_D2x_S_vrr+WPY*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_D2x_S_vrr = PAZ*I_ERI_F3x_S_D2x_S_vrr+WPZ*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2x_S_vrr = PAY*I_ERI_F2xy_S_D2x_S_vrr+WPY*I_ERI_F2xy_S_D2x_S_M1_vrr+oned2z*I_ERI_D2x_S_D2x_S_vrr-rhod2zsq*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2x_S_vrr = PAZ*I_ERI_F2xy_S_D2x_S_vrr+WPZ*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2x_S_vrr = PAZ*I_ERI_F2xz_S_D2x_S_vrr+WPZ*I_ERI_F2xz_S_D2x_S_M1_vrr+oned2z*I_ERI_D2x_S_D2x_S_vrr-rhod2zsq*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2x_S_vrr = PAX*I_ERI_F3y_S_D2x_S_vrr+WPX*I_ERI_F3y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2x_S_vrr = PAZ*I_ERI_Fx2y_S_D2x_S_vrr+WPZ*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2x_S_vrr = PAY*I_ERI_Fx2z_S_D2x_S_vrr+WPY*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2x_S_vrr = PAX*I_ERI_F3z_S_D2x_S_vrr+WPX*I_ERI_F3z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4y_S_D2x_S_vrr = PAY*I_ERI_F3y_S_D2x_S_vrr+WPY*I_ERI_F3y_S_D2x_S_M1_vrr+3*oned2z*I_ERI_D2y_S_D2x_S_vrr-3*rhod2zsq*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_D2x_S_vrr = PAZ*I_ERI_F3y_S_D2x_S_vrr+WPZ*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2x_S_vrr = PAZ*I_ERI_F2yz_S_D2x_S_vrr+WPZ*I_ERI_F2yz_S_D2x_S_M1_vrr+oned2z*I_ERI_D2y_S_D2x_S_vrr-rhod2zsq*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2x_S_vrr = PAY*I_ERI_F3z_S_D2x_S_vrr+WPY*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_D2x_S_vrr = PAZ*I_ERI_F3z_S_D2x_S_vrr+WPZ*I_ERI_F3z_S_D2x_S_M1_vrr+3*oned2z*I_ERI_D2z_S_D2x_S_vrr-3*rhod2zsq*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_Dxy_S_vrr = PAX*I_ERI_F3x_S_Dxy_S_vrr+WPX*I_ERI_F3x_S_Dxy_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Dxy_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G3xy_S_Dxy_S_vrr = PAY*I_ERI_F3x_S_Dxy_S_vrr+WPY*I_ERI_F3x_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G3xz_S_Dxy_S_vrr = PAZ*I_ERI_F3x_S_Dxy_S_vrr+WPZ*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_vrr = PAY*I_ERI_F2xy_S_Dxy_S_vrr+WPY*I_ERI_F2xy_S_Dxy_S_M1_vrr+oned2z*I_ERI_D2x_S_Dxy_S_vrr-rhod2zsq*I_ERI_D2x_S_Dxy_S_M1_vrr+oned2k*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dxy_S_vrr = PAZ*I_ERI_F2xy_S_Dxy_S_vrr+WPZ*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dxy_S_vrr = PAZ*I_ERI_F2xz_S_Dxy_S_vrr+WPZ*I_ERI_F2xz_S_Dxy_S_M1_vrr+oned2z*I_ERI_D2x_S_Dxy_S_vrr-rhod2zsq*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_vrr = PAX*I_ERI_F3y_S_Dxy_S_vrr+WPX*I_ERI_F3y_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dxy_S_vrr = PAZ*I_ERI_Fx2y_S_Dxy_S_vrr+WPZ*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dxy_S_vrr = PAY*I_ERI_Fx2z_S_Dxy_S_vrr+WPY*I_ERI_Fx2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Px_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_vrr = PAX*I_ERI_F3z_S_Dxy_S_vrr+WPX*I_ERI_F3z_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4y_S_Dxy_S_vrr = PAY*I_ERI_F3y_S_Dxy_S_vrr+WPY*I_ERI_F3y_S_Dxy_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Dxy_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_Dxy_S_vrr = PAZ*I_ERI_F3y_S_Dxy_S_vrr+WPZ*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dxy_S_vrr = PAZ*I_ERI_F2yz_S_Dxy_S_vrr+WPZ*I_ERI_F2yz_S_Dxy_S_M1_vrr+oned2z*I_ERI_D2y_S_Dxy_S_vrr-rhod2zsq*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_vrr = PAY*I_ERI_F3z_S_Dxy_S_vrr+WPY*I_ERI_F3z_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_Dxy_S_vrr = PAZ*I_ERI_F3z_S_Dxy_S_vrr+WPZ*I_ERI_F3z_S_Dxy_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Dxy_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4x_S_Dxz_S_vrr = PAX*I_ERI_F3x_S_Dxz_S_vrr+WPX*I_ERI_F3x_S_Dxz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Dxz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xy_S_Dxz_S_vrr = PAY*I_ERI_F3x_S_Dxz_S_vrr+WPY*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G3xz_S_Dxz_S_vrr = PAZ*I_ERI_F3x_S_Dxz_S_vrr+WPZ*I_ERI_F3x_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dxz_S_vrr = PAY*I_ERI_F2xy_S_Dxz_S_vrr+WPY*I_ERI_F2xy_S_Dxz_S_M1_vrr+oned2z*I_ERI_D2x_S_Dxz_S_vrr-rhod2zsq*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dxz_S_vrr = PAZ*I_ERI_F2xy_S_Dxz_S_vrr+WPZ*I_ERI_F2xy_S_Dxz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Px_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dxz_S_vrr = PAZ*I_ERI_F2xz_S_Dxz_S_vrr+WPZ*I_ERI_F2xz_S_Dxz_S_M1_vrr+oned2z*I_ERI_D2x_S_Dxz_S_vrr-rhod2zsq*I_ERI_D2x_S_Dxz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Px_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dxz_S_vrr = PAX*I_ERI_F3y_S_Dxz_S_vrr+WPX*I_ERI_F3y_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dxz_S_vrr = PAZ*I_ERI_Fx2y_S_Dxz_S_vrr+WPZ*I_ERI_Fx2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Px_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dxz_S_vrr = PAY*I_ERI_Fx2z_S_Dxz_S_vrr+WPY*I_ERI_Fx2z_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dxz_S_vrr = PAX*I_ERI_F3z_S_Dxz_S_vrr+WPX*I_ERI_F3z_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_G4y_S_Dxz_S_vrr = PAY*I_ERI_F3y_S_Dxz_S_vrr+WPY*I_ERI_F3y_S_Dxz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Dxz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_G3yz_S_Dxz_S_vrr = PAZ*I_ERI_F3y_S_Dxz_S_vrr+WPZ*I_ERI_F3y_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dxz_S_vrr = PAZ*I_ERI_F2yz_S_Dxz_S_vrr+WPZ*I_ERI_F2yz_S_Dxz_S_M1_vrr+oned2z*I_ERI_D2y_S_Dxz_S_vrr-rhod2zsq*I_ERI_D2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dxz_S_vrr = PAY*I_ERI_F3z_S_Dxz_S_vrr+WPY*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4z_S_Dxz_S_vrr = PAZ*I_ERI_F3z_S_Dxz_S_vrr+WPZ*I_ERI_F3z_S_Dxz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Dxz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxz_S_M1_vrr+oned2k*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_D2y_S_vrr = PAX*I_ERI_F3x_S_D2y_S_vrr+WPX*I_ERI_F3x_S_D2y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_D2y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_D2y_S_vrr = PAY*I_ERI_F3x_S_D2y_S_vrr+WPY*I_ERI_F3x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G3xz_S_D2y_S_vrr = PAZ*I_ERI_F3x_S_D2y_S_vrr+WPZ*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2y_S_vrr = PAY*I_ERI_F2xy_S_D2y_S_vrr+WPY*I_ERI_F2xy_S_D2y_S_M1_vrr+oned2z*I_ERI_D2x_S_D2y_S_vrr-rhod2zsq*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2y_S_vrr = PAZ*I_ERI_F2xy_S_D2y_S_vrr+WPZ*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2y_S_vrr = PAZ*I_ERI_F2xz_S_D2y_S_vrr+WPZ*I_ERI_F2xz_S_D2y_S_M1_vrr+oned2z*I_ERI_D2x_S_D2y_S_vrr-rhod2zsq*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2y_S_vrr = PAX*I_ERI_F3y_S_D2y_S_vrr+WPX*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2y_S_vrr = PAZ*I_ERI_Fx2y_S_D2y_S_vrr+WPZ*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2y_S_vrr = PAY*I_ERI_Fx2z_S_D2y_S_vrr+WPY*I_ERI_Fx2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Py_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2y_S_vrr = PAX*I_ERI_F3z_S_D2y_S_vrr+WPX*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_D2y_S_vrr = PAY*I_ERI_F3y_S_D2y_S_vrr+WPY*I_ERI_F3y_S_D2y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_D2y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_G3yz_S_D2y_S_vrr = PAZ*I_ERI_F3y_S_D2y_S_vrr+WPZ*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2y_S_vrr = PAZ*I_ERI_F2yz_S_D2y_S_vrr+WPZ*I_ERI_F2yz_S_D2y_S_M1_vrr+oned2z*I_ERI_D2y_S_D2y_S_vrr-rhod2zsq*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2y_S_vrr = PAY*I_ERI_F3z_S_D2y_S_vrr+WPY*I_ERI_F3z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4z_S_D2y_S_vrr = PAZ*I_ERI_F3z_S_D2y_S_vrr+WPZ*I_ERI_F3z_S_D2y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_D2y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Dyz_S_vrr = PAX*I_ERI_F3x_S_Dyz_S_vrr+WPX*I_ERI_F3x_S_Dyz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Dyz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_G3xy_S_Dyz_S_vrr = PAY*I_ERI_F3x_S_Dyz_S_vrr+WPY*I_ERI_F3x_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xz_S_Dyz_S_vrr = PAZ*I_ERI_F3x_S_Dyz_S_vrr+WPZ*I_ERI_F3x_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_G2x2y_S_Dyz_S_vrr = PAY*I_ERI_F2xy_S_Dyz_S_vrr+WPY*I_ERI_F2xy_S_Dyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Dyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Dyz_S_vrr = PAZ*I_ERI_F2xy_S_Dyz_S_vrr+WPZ*I_ERI_F2xy_S_Dyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Py_S_M1_vrr;
      Double I_ERI_G2x2z_S_Dyz_S_vrr = PAZ*I_ERI_F2xz_S_Dyz_S_vrr+WPZ*I_ERI_F2xz_S_Dyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Dyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Py_S_M1_vrr;
      Double I_ERI_Gx3y_S_Dyz_S_vrr = PAX*I_ERI_F3y_S_Dyz_S_vrr+WPX*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Dyz_S_vrr = PAZ*I_ERI_Fx2y_S_Dyz_S_vrr+WPZ*I_ERI_Fx2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Py_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Dyz_S_vrr = PAY*I_ERI_Fx2z_S_Dyz_S_vrr+WPY*I_ERI_Fx2z_S_Dyz_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Pz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Dyz_S_vrr = PAX*I_ERI_F3z_S_Dyz_S_vrr+WPX*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4y_S_Dyz_S_vrr = PAY*I_ERI_F3y_S_Dyz_S_vrr+WPY*I_ERI_F3y_S_Dyz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Dyz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_G3yz_S_Dyz_S_vrr = PAZ*I_ERI_F3y_S_Dyz_S_vrr+WPZ*I_ERI_F3y_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_G2y2z_S_Dyz_S_vrr = PAZ*I_ERI_F2yz_S_Dyz_S_vrr+WPZ*I_ERI_F2yz_S_Dyz_S_M1_vrr+oned2z*I_ERI_D2y_S_Dyz_S_vrr-rhod2zsq*I_ERI_D2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Py_S_M1_vrr;
      Double I_ERI_Gy3z_S_Dyz_S_vrr = PAY*I_ERI_F3z_S_Dyz_S_vrr+WPY*I_ERI_F3z_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_G4z_S_Dyz_S_vrr = PAZ*I_ERI_F3z_S_Dyz_S_vrr+WPZ*I_ERI_F3z_S_Dyz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Dyz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Dyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_G4x_S_D2z_S_vrr = PAX*I_ERI_F3x_S_D2z_S_vrr+WPX*I_ERI_F3x_S_D2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_D2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_D2z_S_vrr = PAY*I_ERI_F3x_S_D2z_S_vrr+WPY*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_D2z_S_vrr = PAZ*I_ERI_F3x_S_D2z_S_vrr+WPZ*I_ERI_F3x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_G2x2y_S_D2z_S_vrr = PAY*I_ERI_F2xy_S_D2z_S_vrr+WPY*I_ERI_F2xy_S_D2z_S_M1_vrr+oned2z*I_ERI_D2x_S_D2z_S_vrr-rhod2zsq*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_D2z_S_vrr = PAZ*I_ERI_F2xy_S_D2z_S_vrr+WPZ*I_ERI_F2xy_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Pz_S_M1_vrr;
      Double I_ERI_G2x2z_S_D2z_S_vrr = PAZ*I_ERI_F2xz_S_D2z_S_vrr+WPZ*I_ERI_F2xz_S_D2z_S_M1_vrr+oned2z*I_ERI_D2x_S_D2z_S_vrr-rhod2zsq*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Pz_S_M1_vrr;
      Double I_ERI_Gx3y_S_D2z_S_vrr = PAX*I_ERI_F3y_S_D2z_S_vrr+WPX*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_D2z_S_vrr = PAZ*I_ERI_Fx2y_S_D2z_S_vrr+WPZ*I_ERI_Fx2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Pz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_D2z_S_vrr = PAY*I_ERI_Fx2z_S_D2z_S_vrr+WPY*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_D2z_S_vrr = PAX*I_ERI_F3z_S_D2z_S_vrr+WPX*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_D2z_S_vrr = PAY*I_ERI_F3y_S_D2z_S_vrr+WPY*I_ERI_F3y_S_D2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_D2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_D2z_S_vrr = PAZ*I_ERI_F3y_S_D2z_S_vrr+WPZ*I_ERI_F3y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_G2y2z_S_D2z_S_vrr = PAZ*I_ERI_F2yz_S_D2z_S_vrr+WPZ*I_ERI_F2yz_S_D2z_S_M1_vrr+oned2z*I_ERI_D2y_S_D2z_S_vrr-rhod2zsq*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Pz_S_M1_vrr;
      Double I_ERI_Gy3z_S_D2z_S_vrr = PAY*I_ERI_F3z_S_D2z_S_vrr+WPY*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_D2z_S_vrr = PAZ*I_ERI_F3z_S_D2z_S_vrr+WPZ*I_ERI_F3z_S_D2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_D2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1000003_c_coefs = ic2*jc2*gamma;
      I_ERI_F3x_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1000003_c += SQ_ERI_F_S_D_S_C1000003_c_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1000003_coefs = ic2*jc2;
      I_ERI_F3x_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1000003 += SQ_ERI_F_S_S_S_C1000003_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_P_C1001000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_P_C1001000003_coefs = ic2*jc2_1;
      I_ERI_F3x_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3x_S_S_Px_vrr;
      I_ERI_F2xy_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xy_S_S_Px_vrr;
      I_ERI_F2xz_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xz_S_S_Px_vrr;
      I_ERI_Fx2y_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2y_S_S_Px_vrr;
      I_ERI_Fxyz_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fxyz_S_S_Px_vrr;
      I_ERI_Fx2z_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2z_S_S_Px_vrr;
      I_ERI_F3y_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3y_S_S_Px_vrr;
      I_ERI_F2yz_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2yz_S_S_Px_vrr;
      I_ERI_Fy2z_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fy2z_S_S_Px_vrr;
      I_ERI_F3z_S_S_Px_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3z_S_S_Px_vrr;
      I_ERI_F3x_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3x_S_S_Py_vrr;
      I_ERI_F2xy_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xy_S_S_Py_vrr;
      I_ERI_F2xz_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xz_S_S_Py_vrr;
      I_ERI_Fx2y_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2y_S_S_Py_vrr;
      I_ERI_Fxyz_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fxyz_S_S_Py_vrr;
      I_ERI_Fx2z_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2z_S_S_Py_vrr;
      I_ERI_F3y_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3y_S_S_Py_vrr;
      I_ERI_F2yz_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2yz_S_S_Py_vrr;
      I_ERI_Fy2z_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fy2z_S_S_Py_vrr;
      I_ERI_F3z_S_S_Py_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3z_S_S_Py_vrr;
      I_ERI_F3x_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3x_S_S_Pz_vrr;
      I_ERI_F2xy_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xy_S_S_Pz_vrr;
      I_ERI_F2xz_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2xz_S_S_Pz_vrr;
      I_ERI_Fx2y_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2y_S_S_Pz_vrr;
      I_ERI_Fxyz_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fxyz_S_S_Pz_vrr;
      I_ERI_Fx2z_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fx2z_S_S_Pz_vrr;
      I_ERI_F3y_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3y_S_S_Pz_vrr;
      I_ERI_F2yz_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F2yz_S_S_Pz_vrr;
      I_ERI_Fy2z_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_Fy2z_S_S_Pz_vrr;
      I_ERI_F3z_S_S_Pz_C1001000003 += SQ_ERI_F_S_S_P_C1001000003_coefs*I_ERI_F3z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001000003_coefs = ic2*jc2_1;
      I_ERI_F3x_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001000003 += SQ_ERI_F_S_P_S_C1001000003_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1000003_b_coefs = ic2*jc2*beta;
      I_ERI_G4x_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1000003_b += SQ_ERI_G_S_P_S_C1000003_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000003_b_coefs = ic2*jc2*beta;
      I_ERI_F3x_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000003_b += SQ_ERI_F_S_P_S_C1000003_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1001000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1001000003_c_coefs = ic2*jc2_1*gamma;
      I_ERI_F3x_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_C1001000003_c += SQ_ERI_F_S_F_S_C1001000003_c_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001000003_c_coefs = ic2*jc2_1*gamma;
      I_ERI_F3x_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001000003_c += SQ_ERI_F_S_D_S_C1001000003_c_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1000003_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1000003_d_coefs = ic2*jc2*delta;
      I_ERI_F3x_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1000003_d += SQ_ERI_F_S_D_S_C1000003_d_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1000003_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000003_d_coefs = ic2*jc2*delta;
      I_ERI_F3x_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000003_d += SQ_ERI_F_S_P_S_C1000003_d_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_C1001000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_C1001000003_b_coefs = ic2*jc2_1*beta;
      I_ERI_G4x_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_C1001000003_b += SQ_ERI_G_S_D_S_C1001000003_b_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1001000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1001000003_b_coefs = ic2*jc2_1*beta;
      I_ERI_G4x_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1001000003_b += SQ_ERI_G_S_P_S_C1001000003_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001000003_b_coefs = ic2*jc2_1*beta;
      I_ERI_F3x_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001000003_b += SQ_ERI_F_S_D_S_C1001000003_b_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001000003_b_coefs = ic2*jc2_1*beta;
      I_ERI_F3x_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001000003_b += SQ_ERI_F_S_P_S_C1001000003_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1001000003_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1001000003_d_coefs = ic2*jc2_1*delta;
      I_ERI_F3x_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_C1001000003_d += SQ_ERI_F_S_F_S_C1001000003_d_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1001000003_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001000003_d_coefs = ic2*jc2_1*delta;
      I_ERI_F3x_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001000003_d += SQ_ERI_F_S_D_S_C1001000003_d_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001000003_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001000003_d_coefs = ic2*jc2_1*delta;
      I_ERI_F3x_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001000003_d += SQ_ERI_F_S_P_S_C1001000003_d_coefs*I_ERI_F3z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001000003_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001000003_b
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1001000003_b = I_ERI_F3x_S_D2x_S_C1001000003_b+CDX*I_ERI_F3x_S_Px_S_C1001000003_b;
  Double I_ERI_F2xy_S_Px_Px_C1001000003_b = I_ERI_F2xy_S_D2x_S_C1001000003_b+CDX*I_ERI_F2xy_S_Px_S_C1001000003_b;
  Double I_ERI_F2xz_S_Px_Px_C1001000003_b = I_ERI_F2xz_S_D2x_S_C1001000003_b+CDX*I_ERI_F2xz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Px_Px_C1001000003_b = I_ERI_Fx2y_S_D2x_S_C1001000003_b+CDX*I_ERI_Fx2y_S_Px_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Px_Px_C1001000003_b = I_ERI_Fxyz_S_D2x_S_C1001000003_b+CDX*I_ERI_Fxyz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Px_Px_C1001000003_b = I_ERI_Fx2z_S_D2x_S_C1001000003_b+CDX*I_ERI_Fx2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3y_S_Px_Px_C1001000003_b = I_ERI_F3y_S_D2x_S_C1001000003_b+CDX*I_ERI_F3y_S_Px_S_C1001000003_b;
  Double I_ERI_F2yz_S_Px_Px_C1001000003_b = I_ERI_F2yz_S_D2x_S_C1001000003_b+CDX*I_ERI_F2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Px_Px_C1001000003_b = I_ERI_Fy2z_S_D2x_S_C1001000003_b+CDX*I_ERI_Fy2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3z_S_Px_Px_C1001000003_b = I_ERI_F3z_S_D2x_S_C1001000003_b+CDX*I_ERI_F3z_S_Px_S_C1001000003_b;
  Double I_ERI_F3x_S_Py_Px_C1001000003_b = I_ERI_F3x_S_Dxy_S_C1001000003_b+CDX*I_ERI_F3x_S_Py_S_C1001000003_b;
  Double I_ERI_F2xy_S_Py_Px_C1001000003_b = I_ERI_F2xy_S_Dxy_S_C1001000003_b+CDX*I_ERI_F2xy_S_Py_S_C1001000003_b;
  Double I_ERI_F2xz_S_Py_Px_C1001000003_b = I_ERI_F2xz_S_Dxy_S_C1001000003_b+CDX*I_ERI_F2xz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Py_Px_C1001000003_b = I_ERI_Fx2y_S_Dxy_S_C1001000003_b+CDX*I_ERI_Fx2y_S_Py_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Py_Px_C1001000003_b = I_ERI_Fxyz_S_Dxy_S_C1001000003_b+CDX*I_ERI_Fxyz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Py_Px_C1001000003_b = I_ERI_Fx2z_S_Dxy_S_C1001000003_b+CDX*I_ERI_Fx2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3y_S_Py_Px_C1001000003_b = I_ERI_F3y_S_Dxy_S_C1001000003_b+CDX*I_ERI_F3y_S_Py_S_C1001000003_b;
  Double I_ERI_F2yz_S_Py_Px_C1001000003_b = I_ERI_F2yz_S_Dxy_S_C1001000003_b+CDX*I_ERI_F2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Py_Px_C1001000003_b = I_ERI_Fy2z_S_Dxy_S_C1001000003_b+CDX*I_ERI_Fy2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3z_S_Py_Px_C1001000003_b = I_ERI_F3z_S_Dxy_S_C1001000003_b+CDX*I_ERI_F3z_S_Py_S_C1001000003_b;
  Double I_ERI_F3x_S_Pz_Px_C1001000003_b = I_ERI_F3x_S_Dxz_S_C1001000003_b+CDX*I_ERI_F3x_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xy_S_Pz_Px_C1001000003_b = I_ERI_F2xy_S_Dxz_S_C1001000003_b+CDX*I_ERI_F2xy_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xz_S_Pz_Px_C1001000003_b = I_ERI_F2xz_S_Dxz_S_C1001000003_b+CDX*I_ERI_F2xz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Pz_Px_C1001000003_b = I_ERI_Fx2y_S_Dxz_S_C1001000003_b+CDX*I_ERI_Fx2y_S_Pz_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Pz_Px_C1001000003_b = I_ERI_Fxyz_S_Dxz_S_C1001000003_b+CDX*I_ERI_Fxyz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Pz_Px_C1001000003_b = I_ERI_Fx2z_S_Dxz_S_C1001000003_b+CDX*I_ERI_Fx2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3y_S_Pz_Px_C1001000003_b = I_ERI_F3y_S_Dxz_S_C1001000003_b+CDX*I_ERI_F3y_S_Pz_S_C1001000003_b;
  Double I_ERI_F2yz_S_Pz_Px_C1001000003_b = I_ERI_F2yz_S_Dxz_S_C1001000003_b+CDX*I_ERI_F2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Pz_Px_C1001000003_b = I_ERI_Fy2z_S_Dxz_S_C1001000003_b+CDX*I_ERI_Fy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3z_S_Pz_Px_C1001000003_b = I_ERI_F3z_S_Dxz_S_C1001000003_b+CDX*I_ERI_F3z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3x_S_Px_Py_C1001000003_b = I_ERI_F3x_S_Dxy_S_C1001000003_b+CDY*I_ERI_F3x_S_Px_S_C1001000003_b;
  Double I_ERI_F2xy_S_Px_Py_C1001000003_b = I_ERI_F2xy_S_Dxy_S_C1001000003_b+CDY*I_ERI_F2xy_S_Px_S_C1001000003_b;
  Double I_ERI_F2xz_S_Px_Py_C1001000003_b = I_ERI_F2xz_S_Dxy_S_C1001000003_b+CDY*I_ERI_F2xz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Px_Py_C1001000003_b = I_ERI_Fx2y_S_Dxy_S_C1001000003_b+CDY*I_ERI_Fx2y_S_Px_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Px_Py_C1001000003_b = I_ERI_Fxyz_S_Dxy_S_C1001000003_b+CDY*I_ERI_Fxyz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Px_Py_C1001000003_b = I_ERI_Fx2z_S_Dxy_S_C1001000003_b+CDY*I_ERI_Fx2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3y_S_Px_Py_C1001000003_b = I_ERI_F3y_S_Dxy_S_C1001000003_b+CDY*I_ERI_F3y_S_Px_S_C1001000003_b;
  Double I_ERI_F2yz_S_Px_Py_C1001000003_b = I_ERI_F2yz_S_Dxy_S_C1001000003_b+CDY*I_ERI_F2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Px_Py_C1001000003_b = I_ERI_Fy2z_S_Dxy_S_C1001000003_b+CDY*I_ERI_Fy2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3z_S_Px_Py_C1001000003_b = I_ERI_F3z_S_Dxy_S_C1001000003_b+CDY*I_ERI_F3z_S_Px_S_C1001000003_b;
  Double I_ERI_F3x_S_Py_Py_C1001000003_b = I_ERI_F3x_S_D2y_S_C1001000003_b+CDY*I_ERI_F3x_S_Py_S_C1001000003_b;
  Double I_ERI_F2xy_S_Py_Py_C1001000003_b = I_ERI_F2xy_S_D2y_S_C1001000003_b+CDY*I_ERI_F2xy_S_Py_S_C1001000003_b;
  Double I_ERI_F2xz_S_Py_Py_C1001000003_b = I_ERI_F2xz_S_D2y_S_C1001000003_b+CDY*I_ERI_F2xz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Py_Py_C1001000003_b = I_ERI_Fx2y_S_D2y_S_C1001000003_b+CDY*I_ERI_Fx2y_S_Py_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Py_Py_C1001000003_b = I_ERI_Fxyz_S_D2y_S_C1001000003_b+CDY*I_ERI_Fxyz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Py_Py_C1001000003_b = I_ERI_Fx2z_S_D2y_S_C1001000003_b+CDY*I_ERI_Fx2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3y_S_Py_Py_C1001000003_b = I_ERI_F3y_S_D2y_S_C1001000003_b+CDY*I_ERI_F3y_S_Py_S_C1001000003_b;
  Double I_ERI_F2yz_S_Py_Py_C1001000003_b = I_ERI_F2yz_S_D2y_S_C1001000003_b+CDY*I_ERI_F2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Py_Py_C1001000003_b = I_ERI_Fy2z_S_D2y_S_C1001000003_b+CDY*I_ERI_Fy2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3z_S_Py_Py_C1001000003_b = I_ERI_F3z_S_D2y_S_C1001000003_b+CDY*I_ERI_F3z_S_Py_S_C1001000003_b;
  Double I_ERI_F3x_S_Pz_Py_C1001000003_b = I_ERI_F3x_S_Dyz_S_C1001000003_b+CDY*I_ERI_F3x_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xy_S_Pz_Py_C1001000003_b = I_ERI_F2xy_S_Dyz_S_C1001000003_b+CDY*I_ERI_F2xy_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xz_S_Pz_Py_C1001000003_b = I_ERI_F2xz_S_Dyz_S_C1001000003_b+CDY*I_ERI_F2xz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Pz_Py_C1001000003_b = I_ERI_Fx2y_S_Dyz_S_C1001000003_b+CDY*I_ERI_Fx2y_S_Pz_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Pz_Py_C1001000003_b = I_ERI_Fxyz_S_Dyz_S_C1001000003_b+CDY*I_ERI_Fxyz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Pz_Py_C1001000003_b = I_ERI_Fx2z_S_Dyz_S_C1001000003_b+CDY*I_ERI_Fx2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3y_S_Pz_Py_C1001000003_b = I_ERI_F3y_S_Dyz_S_C1001000003_b+CDY*I_ERI_F3y_S_Pz_S_C1001000003_b;
  Double I_ERI_F2yz_S_Pz_Py_C1001000003_b = I_ERI_F2yz_S_Dyz_S_C1001000003_b+CDY*I_ERI_F2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Pz_Py_C1001000003_b = I_ERI_Fy2z_S_Dyz_S_C1001000003_b+CDY*I_ERI_Fy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3z_S_Pz_Py_C1001000003_b = I_ERI_F3z_S_Dyz_S_C1001000003_b+CDY*I_ERI_F3z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3x_S_Px_Pz_C1001000003_b = I_ERI_F3x_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F3x_S_Px_S_C1001000003_b;
  Double I_ERI_F2xy_S_Px_Pz_C1001000003_b = I_ERI_F2xy_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F2xy_S_Px_S_C1001000003_b;
  Double I_ERI_F2xz_S_Px_Pz_C1001000003_b = I_ERI_F2xz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F2xz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Px_Pz_C1001000003_b = I_ERI_Fx2y_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Fx2y_S_Px_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Px_Pz_C1001000003_b = I_ERI_Fxyz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Fxyz_S_Px_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Px_Pz_C1001000003_b = I_ERI_Fx2z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Fx2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3y_S_Px_Pz_C1001000003_b = I_ERI_F3y_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F3y_S_Px_S_C1001000003_b;
  Double I_ERI_F2yz_S_Px_Pz_C1001000003_b = I_ERI_F2yz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Px_Pz_C1001000003_b = I_ERI_Fy2z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Fy2z_S_Px_S_C1001000003_b;
  Double I_ERI_F3z_S_Px_Pz_C1001000003_b = I_ERI_F3z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_F3z_S_Px_S_C1001000003_b;
  Double I_ERI_F3x_S_Py_Pz_C1001000003_b = I_ERI_F3x_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F3x_S_Py_S_C1001000003_b;
  Double I_ERI_F2xy_S_Py_Pz_C1001000003_b = I_ERI_F2xy_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F2xy_S_Py_S_C1001000003_b;
  Double I_ERI_F2xz_S_Py_Pz_C1001000003_b = I_ERI_F2xz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F2xz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Py_Pz_C1001000003_b = I_ERI_Fx2y_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Fx2y_S_Py_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Py_Pz_C1001000003_b = I_ERI_Fxyz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Fxyz_S_Py_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Py_Pz_C1001000003_b = I_ERI_Fx2z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Fx2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3y_S_Py_Pz_C1001000003_b = I_ERI_F3y_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F3y_S_Py_S_C1001000003_b;
  Double I_ERI_F2yz_S_Py_Pz_C1001000003_b = I_ERI_F2yz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Py_Pz_C1001000003_b = I_ERI_Fy2z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Fy2z_S_Py_S_C1001000003_b;
  Double I_ERI_F3z_S_Py_Pz_C1001000003_b = I_ERI_F3z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_F3z_S_Py_S_C1001000003_b;
  Double I_ERI_F3x_S_Pz_Pz_C1001000003_b = I_ERI_F3x_S_D2z_S_C1001000003_b+CDZ*I_ERI_F3x_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xy_S_Pz_Pz_C1001000003_b = I_ERI_F2xy_S_D2z_S_C1001000003_b+CDZ*I_ERI_F2xy_S_Pz_S_C1001000003_b;
  Double I_ERI_F2xz_S_Pz_Pz_C1001000003_b = I_ERI_F2xz_S_D2z_S_C1001000003_b+CDZ*I_ERI_F2xz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2y_S_Pz_Pz_C1001000003_b = I_ERI_Fx2y_S_D2z_S_C1001000003_b+CDZ*I_ERI_Fx2y_S_Pz_S_C1001000003_b;
  Double I_ERI_Fxyz_S_Pz_Pz_C1001000003_b = I_ERI_Fxyz_S_D2z_S_C1001000003_b+CDZ*I_ERI_Fxyz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fx2z_S_Pz_Pz_C1001000003_b = I_ERI_Fx2z_S_D2z_S_C1001000003_b+CDZ*I_ERI_Fx2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3y_S_Pz_Pz_C1001000003_b = I_ERI_F3y_S_D2z_S_C1001000003_b+CDZ*I_ERI_F3y_S_Pz_S_C1001000003_b;
  Double I_ERI_F2yz_S_Pz_Pz_C1001000003_b = I_ERI_F2yz_S_D2z_S_C1001000003_b+CDZ*I_ERI_F2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Fy2z_S_Pz_Pz_C1001000003_b = I_ERI_Fy2z_S_D2z_S_C1001000003_b+CDZ*I_ERI_Fy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_F3z_S_Pz_Pz_C1001000003_b = I_ERI_F3z_S_D2z_S_C1001000003_b+CDZ*I_ERI_F3z_S_Pz_S_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_P_C1001000003_b
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_C1001000003_b
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1001000003_b
   ************************************************************/
  Double I_ERI_G4x_S_Px_Px_C1001000003_b = I_ERI_G4x_S_D2x_S_C1001000003_b+CDX*I_ERI_G4x_S_Px_S_C1001000003_b;
  Double I_ERI_G3xy_S_Px_Px_C1001000003_b = I_ERI_G3xy_S_D2x_S_C1001000003_b+CDX*I_ERI_G3xy_S_Px_S_C1001000003_b;
  Double I_ERI_G3xz_S_Px_Px_C1001000003_b = I_ERI_G3xz_S_D2x_S_C1001000003_b+CDX*I_ERI_G3xz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Px_Px_C1001000003_b = I_ERI_G2x2y_S_D2x_S_C1001000003_b+CDX*I_ERI_G2x2y_S_Px_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Px_Px_C1001000003_b = I_ERI_G2xyz_S_D2x_S_C1001000003_b+CDX*I_ERI_G2xyz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Px_Px_C1001000003_b = I_ERI_G2x2z_S_D2x_S_C1001000003_b+CDX*I_ERI_G2x2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Px_Px_C1001000003_b = I_ERI_Gx3y_S_D2x_S_C1001000003_b+CDX*I_ERI_Gx3y_S_Px_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Px_Px_C1001000003_b = I_ERI_Gx2yz_S_D2x_S_C1001000003_b+CDX*I_ERI_Gx2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Px_Px_C1001000003_b = I_ERI_Gxy2z_S_D2x_S_C1001000003_b+CDX*I_ERI_Gxy2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Px_Px_C1001000003_b = I_ERI_Gx3z_S_D2x_S_C1001000003_b+CDX*I_ERI_Gx3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4y_S_Px_Px_C1001000003_b = I_ERI_G4y_S_D2x_S_C1001000003_b+CDX*I_ERI_G4y_S_Px_S_C1001000003_b;
  Double I_ERI_G3yz_S_Px_Px_C1001000003_b = I_ERI_G3yz_S_D2x_S_C1001000003_b+CDX*I_ERI_G3yz_S_Px_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Px_Px_C1001000003_b = I_ERI_G2y2z_S_D2x_S_C1001000003_b+CDX*I_ERI_G2y2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Px_Px_C1001000003_b = I_ERI_Gy3z_S_D2x_S_C1001000003_b+CDX*I_ERI_Gy3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4z_S_Px_Px_C1001000003_b = I_ERI_G4z_S_D2x_S_C1001000003_b+CDX*I_ERI_G4z_S_Px_S_C1001000003_b;
  Double I_ERI_G4x_S_Py_Px_C1001000003_b = I_ERI_G4x_S_Dxy_S_C1001000003_b+CDX*I_ERI_G4x_S_Py_S_C1001000003_b;
  Double I_ERI_G3xy_S_Py_Px_C1001000003_b = I_ERI_G3xy_S_Dxy_S_C1001000003_b+CDX*I_ERI_G3xy_S_Py_S_C1001000003_b;
  Double I_ERI_G3xz_S_Py_Px_C1001000003_b = I_ERI_G3xz_S_Dxy_S_C1001000003_b+CDX*I_ERI_G3xz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Py_Px_C1001000003_b = I_ERI_G2x2y_S_Dxy_S_C1001000003_b+CDX*I_ERI_G2x2y_S_Py_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Py_Px_C1001000003_b = I_ERI_G2xyz_S_Dxy_S_C1001000003_b+CDX*I_ERI_G2xyz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Py_Px_C1001000003_b = I_ERI_G2x2z_S_Dxy_S_C1001000003_b+CDX*I_ERI_G2x2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Py_Px_C1001000003_b = I_ERI_Gx3y_S_Dxy_S_C1001000003_b+CDX*I_ERI_Gx3y_S_Py_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Py_Px_C1001000003_b = I_ERI_Gx2yz_S_Dxy_S_C1001000003_b+CDX*I_ERI_Gx2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Py_Px_C1001000003_b = I_ERI_Gxy2z_S_Dxy_S_C1001000003_b+CDX*I_ERI_Gxy2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Py_Px_C1001000003_b = I_ERI_Gx3z_S_Dxy_S_C1001000003_b+CDX*I_ERI_Gx3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4y_S_Py_Px_C1001000003_b = I_ERI_G4y_S_Dxy_S_C1001000003_b+CDX*I_ERI_G4y_S_Py_S_C1001000003_b;
  Double I_ERI_G3yz_S_Py_Px_C1001000003_b = I_ERI_G3yz_S_Dxy_S_C1001000003_b+CDX*I_ERI_G3yz_S_Py_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Py_Px_C1001000003_b = I_ERI_G2y2z_S_Dxy_S_C1001000003_b+CDX*I_ERI_G2y2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Py_Px_C1001000003_b = I_ERI_Gy3z_S_Dxy_S_C1001000003_b+CDX*I_ERI_Gy3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4z_S_Py_Px_C1001000003_b = I_ERI_G4z_S_Dxy_S_C1001000003_b+CDX*I_ERI_G4z_S_Py_S_C1001000003_b;
  Double I_ERI_G4x_S_Pz_Px_C1001000003_b = I_ERI_G4x_S_Dxz_S_C1001000003_b+CDX*I_ERI_G4x_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xy_S_Pz_Px_C1001000003_b = I_ERI_G3xy_S_Dxz_S_C1001000003_b+CDX*I_ERI_G3xy_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xz_S_Pz_Px_C1001000003_b = I_ERI_G3xz_S_Dxz_S_C1001000003_b+CDX*I_ERI_G3xz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Pz_Px_C1001000003_b = I_ERI_G2x2y_S_Dxz_S_C1001000003_b+CDX*I_ERI_G2x2y_S_Pz_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Pz_Px_C1001000003_b = I_ERI_G2xyz_S_Dxz_S_C1001000003_b+CDX*I_ERI_G2xyz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Pz_Px_C1001000003_b = I_ERI_G2x2z_S_Dxz_S_C1001000003_b+CDX*I_ERI_G2x2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Pz_Px_C1001000003_b = I_ERI_Gx3y_S_Dxz_S_C1001000003_b+CDX*I_ERI_Gx3y_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Pz_Px_C1001000003_b = I_ERI_Gx2yz_S_Dxz_S_C1001000003_b+CDX*I_ERI_Gx2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Pz_Px_C1001000003_b = I_ERI_Gxy2z_S_Dxz_S_C1001000003_b+CDX*I_ERI_Gxy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Pz_Px_C1001000003_b = I_ERI_Gx3z_S_Dxz_S_C1001000003_b+CDX*I_ERI_Gx3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4y_S_Pz_Px_C1001000003_b = I_ERI_G4y_S_Dxz_S_C1001000003_b+CDX*I_ERI_G4y_S_Pz_S_C1001000003_b;
  Double I_ERI_G3yz_S_Pz_Px_C1001000003_b = I_ERI_G3yz_S_Dxz_S_C1001000003_b+CDX*I_ERI_G3yz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Pz_Px_C1001000003_b = I_ERI_G2y2z_S_Dxz_S_C1001000003_b+CDX*I_ERI_G2y2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Pz_Px_C1001000003_b = I_ERI_Gy3z_S_Dxz_S_C1001000003_b+CDX*I_ERI_Gy3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4z_S_Pz_Px_C1001000003_b = I_ERI_G4z_S_Dxz_S_C1001000003_b+CDX*I_ERI_G4z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4x_S_Px_Py_C1001000003_b = I_ERI_G4x_S_Dxy_S_C1001000003_b+CDY*I_ERI_G4x_S_Px_S_C1001000003_b;
  Double I_ERI_G3xy_S_Px_Py_C1001000003_b = I_ERI_G3xy_S_Dxy_S_C1001000003_b+CDY*I_ERI_G3xy_S_Px_S_C1001000003_b;
  Double I_ERI_G3xz_S_Px_Py_C1001000003_b = I_ERI_G3xz_S_Dxy_S_C1001000003_b+CDY*I_ERI_G3xz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Px_Py_C1001000003_b = I_ERI_G2x2y_S_Dxy_S_C1001000003_b+CDY*I_ERI_G2x2y_S_Px_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Px_Py_C1001000003_b = I_ERI_G2xyz_S_Dxy_S_C1001000003_b+CDY*I_ERI_G2xyz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Px_Py_C1001000003_b = I_ERI_G2x2z_S_Dxy_S_C1001000003_b+CDY*I_ERI_G2x2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Px_Py_C1001000003_b = I_ERI_Gx3y_S_Dxy_S_C1001000003_b+CDY*I_ERI_Gx3y_S_Px_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Px_Py_C1001000003_b = I_ERI_Gx2yz_S_Dxy_S_C1001000003_b+CDY*I_ERI_Gx2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Px_Py_C1001000003_b = I_ERI_Gxy2z_S_Dxy_S_C1001000003_b+CDY*I_ERI_Gxy2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Px_Py_C1001000003_b = I_ERI_Gx3z_S_Dxy_S_C1001000003_b+CDY*I_ERI_Gx3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4y_S_Px_Py_C1001000003_b = I_ERI_G4y_S_Dxy_S_C1001000003_b+CDY*I_ERI_G4y_S_Px_S_C1001000003_b;
  Double I_ERI_G3yz_S_Px_Py_C1001000003_b = I_ERI_G3yz_S_Dxy_S_C1001000003_b+CDY*I_ERI_G3yz_S_Px_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Px_Py_C1001000003_b = I_ERI_G2y2z_S_Dxy_S_C1001000003_b+CDY*I_ERI_G2y2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Px_Py_C1001000003_b = I_ERI_Gy3z_S_Dxy_S_C1001000003_b+CDY*I_ERI_Gy3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4z_S_Px_Py_C1001000003_b = I_ERI_G4z_S_Dxy_S_C1001000003_b+CDY*I_ERI_G4z_S_Px_S_C1001000003_b;
  Double I_ERI_G4x_S_Py_Py_C1001000003_b = I_ERI_G4x_S_D2y_S_C1001000003_b+CDY*I_ERI_G4x_S_Py_S_C1001000003_b;
  Double I_ERI_G3xy_S_Py_Py_C1001000003_b = I_ERI_G3xy_S_D2y_S_C1001000003_b+CDY*I_ERI_G3xy_S_Py_S_C1001000003_b;
  Double I_ERI_G3xz_S_Py_Py_C1001000003_b = I_ERI_G3xz_S_D2y_S_C1001000003_b+CDY*I_ERI_G3xz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Py_Py_C1001000003_b = I_ERI_G2x2y_S_D2y_S_C1001000003_b+CDY*I_ERI_G2x2y_S_Py_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Py_Py_C1001000003_b = I_ERI_G2xyz_S_D2y_S_C1001000003_b+CDY*I_ERI_G2xyz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Py_Py_C1001000003_b = I_ERI_G2x2z_S_D2y_S_C1001000003_b+CDY*I_ERI_G2x2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Py_Py_C1001000003_b = I_ERI_Gx3y_S_D2y_S_C1001000003_b+CDY*I_ERI_Gx3y_S_Py_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Py_Py_C1001000003_b = I_ERI_Gx2yz_S_D2y_S_C1001000003_b+CDY*I_ERI_Gx2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Py_Py_C1001000003_b = I_ERI_Gxy2z_S_D2y_S_C1001000003_b+CDY*I_ERI_Gxy2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Py_Py_C1001000003_b = I_ERI_Gx3z_S_D2y_S_C1001000003_b+CDY*I_ERI_Gx3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4y_S_Py_Py_C1001000003_b = I_ERI_G4y_S_D2y_S_C1001000003_b+CDY*I_ERI_G4y_S_Py_S_C1001000003_b;
  Double I_ERI_G3yz_S_Py_Py_C1001000003_b = I_ERI_G3yz_S_D2y_S_C1001000003_b+CDY*I_ERI_G3yz_S_Py_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Py_Py_C1001000003_b = I_ERI_G2y2z_S_D2y_S_C1001000003_b+CDY*I_ERI_G2y2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Py_Py_C1001000003_b = I_ERI_Gy3z_S_D2y_S_C1001000003_b+CDY*I_ERI_Gy3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4z_S_Py_Py_C1001000003_b = I_ERI_G4z_S_D2y_S_C1001000003_b+CDY*I_ERI_G4z_S_Py_S_C1001000003_b;
  Double I_ERI_G4x_S_Pz_Py_C1001000003_b = I_ERI_G4x_S_Dyz_S_C1001000003_b+CDY*I_ERI_G4x_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xy_S_Pz_Py_C1001000003_b = I_ERI_G3xy_S_Dyz_S_C1001000003_b+CDY*I_ERI_G3xy_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xz_S_Pz_Py_C1001000003_b = I_ERI_G3xz_S_Dyz_S_C1001000003_b+CDY*I_ERI_G3xz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Pz_Py_C1001000003_b = I_ERI_G2x2y_S_Dyz_S_C1001000003_b+CDY*I_ERI_G2x2y_S_Pz_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Pz_Py_C1001000003_b = I_ERI_G2xyz_S_Dyz_S_C1001000003_b+CDY*I_ERI_G2xyz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Pz_Py_C1001000003_b = I_ERI_G2x2z_S_Dyz_S_C1001000003_b+CDY*I_ERI_G2x2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Pz_Py_C1001000003_b = I_ERI_Gx3y_S_Dyz_S_C1001000003_b+CDY*I_ERI_Gx3y_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Pz_Py_C1001000003_b = I_ERI_Gx2yz_S_Dyz_S_C1001000003_b+CDY*I_ERI_Gx2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Pz_Py_C1001000003_b = I_ERI_Gxy2z_S_Dyz_S_C1001000003_b+CDY*I_ERI_Gxy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Pz_Py_C1001000003_b = I_ERI_Gx3z_S_Dyz_S_C1001000003_b+CDY*I_ERI_Gx3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4y_S_Pz_Py_C1001000003_b = I_ERI_G4y_S_Dyz_S_C1001000003_b+CDY*I_ERI_G4y_S_Pz_S_C1001000003_b;
  Double I_ERI_G3yz_S_Pz_Py_C1001000003_b = I_ERI_G3yz_S_Dyz_S_C1001000003_b+CDY*I_ERI_G3yz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Pz_Py_C1001000003_b = I_ERI_G2y2z_S_Dyz_S_C1001000003_b+CDY*I_ERI_G2y2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Pz_Py_C1001000003_b = I_ERI_Gy3z_S_Dyz_S_C1001000003_b+CDY*I_ERI_Gy3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4z_S_Pz_Py_C1001000003_b = I_ERI_G4z_S_Dyz_S_C1001000003_b+CDY*I_ERI_G4z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4x_S_Px_Pz_C1001000003_b = I_ERI_G4x_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G4x_S_Px_S_C1001000003_b;
  Double I_ERI_G3xy_S_Px_Pz_C1001000003_b = I_ERI_G3xy_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G3xy_S_Px_S_C1001000003_b;
  Double I_ERI_G3xz_S_Px_Pz_C1001000003_b = I_ERI_G3xz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G3xz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Px_Pz_C1001000003_b = I_ERI_G2x2y_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G2x2y_S_Px_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Px_Pz_C1001000003_b = I_ERI_G2xyz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G2xyz_S_Px_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Px_Pz_C1001000003_b = I_ERI_G2x2z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G2x2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Px_Pz_C1001000003_b = I_ERI_Gx3y_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Gx3y_S_Px_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Px_Pz_C1001000003_b = I_ERI_Gx2yz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Gx2yz_S_Px_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Px_Pz_C1001000003_b = I_ERI_Gxy2z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Gxy2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Px_Pz_C1001000003_b = I_ERI_Gx3z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Gx3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4y_S_Px_Pz_C1001000003_b = I_ERI_G4y_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G4y_S_Px_S_C1001000003_b;
  Double I_ERI_G3yz_S_Px_Pz_C1001000003_b = I_ERI_G3yz_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G3yz_S_Px_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Px_Pz_C1001000003_b = I_ERI_G2y2z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G2y2z_S_Px_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Px_Pz_C1001000003_b = I_ERI_Gy3z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_Gy3z_S_Px_S_C1001000003_b;
  Double I_ERI_G4z_S_Px_Pz_C1001000003_b = I_ERI_G4z_S_Dxz_S_C1001000003_b+CDZ*I_ERI_G4z_S_Px_S_C1001000003_b;
  Double I_ERI_G4x_S_Py_Pz_C1001000003_b = I_ERI_G4x_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G4x_S_Py_S_C1001000003_b;
  Double I_ERI_G3xy_S_Py_Pz_C1001000003_b = I_ERI_G3xy_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G3xy_S_Py_S_C1001000003_b;
  Double I_ERI_G3xz_S_Py_Pz_C1001000003_b = I_ERI_G3xz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G3xz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Py_Pz_C1001000003_b = I_ERI_G2x2y_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G2x2y_S_Py_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Py_Pz_C1001000003_b = I_ERI_G2xyz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G2xyz_S_Py_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Py_Pz_C1001000003_b = I_ERI_G2x2z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G2x2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Py_Pz_C1001000003_b = I_ERI_Gx3y_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Gx3y_S_Py_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Py_Pz_C1001000003_b = I_ERI_Gx2yz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Gx2yz_S_Py_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Py_Pz_C1001000003_b = I_ERI_Gxy2z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Gxy2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Py_Pz_C1001000003_b = I_ERI_Gx3z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Gx3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4y_S_Py_Pz_C1001000003_b = I_ERI_G4y_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G4y_S_Py_S_C1001000003_b;
  Double I_ERI_G3yz_S_Py_Pz_C1001000003_b = I_ERI_G3yz_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G3yz_S_Py_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Py_Pz_C1001000003_b = I_ERI_G2y2z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G2y2z_S_Py_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Py_Pz_C1001000003_b = I_ERI_Gy3z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_Gy3z_S_Py_S_C1001000003_b;
  Double I_ERI_G4z_S_Py_Pz_C1001000003_b = I_ERI_G4z_S_Dyz_S_C1001000003_b+CDZ*I_ERI_G4z_S_Py_S_C1001000003_b;
  Double I_ERI_G4x_S_Pz_Pz_C1001000003_b = I_ERI_G4x_S_D2z_S_C1001000003_b+CDZ*I_ERI_G4x_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xy_S_Pz_Pz_C1001000003_b = I_ERI_G3xy_S_D2z_S_C1001000003_b+CDZ*I_ERI_G3xy_S_Pz_S_C1001000003_b;
  Double I_ERI_G3xz_S_Pz_Pz_C1001000003_b = I_ERI_G3xz_S_D2z_S_C1001000003_b+CDZ*I_ERI_G3xz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2y_S_Pz_Pz_C1001000003_b = I_ERI_G2x2y_S_D2z_S_C1001000003_b+CDZ*I_ERI_G2x2y_S_Pz_S_C1001000003_b;
  Double I_ERI_G2xyz_S_Pz_Pz_C1001000003_b = I_ERI_G2xyz_S_D2z_S_C1001000003_b+CDZ*I_ERI_G2xyz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2x2z_S_Pz_Pz_C1001000003_b = I_ERI_G2x2z_S_D2z_S_C1001000003_b+CDZ*I_ERI_G2x2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3y_S_Pz_Pz_C1001000003_b = I_ERI_Gx3y_S_D2z_S_C1001000003_b+CDZ*I_ERI_Gx3y_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx2yz_S_Pz_Pz_C1001000003_b = I_ERI_Gx2yz_S_D2z_S_C1001000003_b+CDZ*I_ERI_Gx2yz_S_Pz_S_C1001000003_b;
  Double I_ERI_Gxy2z_S_Pz_Pz_C1001000003_b = I_ERI_Gxy2z_S_D2z_S_C1001000003_b+CDZ*I_ERI_Gxy2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gx3z_S_Pz_Pz_C1001000003_b = I_ERI_Gx3z_S_D2z_S_C1001000003_b+CDZ*I_ERI_Gx3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4y_S_Pz_Pz_C1001000003_b = I_ERI_G4y_S_D2z_S_C1001000003_b+CDZ*I_ERI_G4y_S_Pz_S_C1001000003_b;
  Double I_ERI_G3yz_S_Pz_Pz_C1001000003_b = I_ERI_G3yz_S_D2z_S_C1001000003_b+CDZ*I_ERI_G3yz_S_Pz_S_C1001000003_b;
  Double I_ERI_G2y2z_S_Pz_Pz_C1001000003_b = I_ERI_G2y2z_S_D2z_S_C1001000003_b+CDZ*I_ERI_G2y2z_S_Pz_S_C1001000003_b;
  Double I_ERI_Gy3z_S_Pz_Pz_C1001000003_b = I_ERI_Gy3z_S_D2z_S_C1001000003_b+CDZ*I_ERI_Gy3z_S_Pz_S_C1001000003_b;
  Double I_ERI_G4z_S_Pz_Pz_C1001000003_b = I_ERI_G4z_S_D2z_S_C1001000003_b+CDZ*I_ERI_G4z_S_Pz_S_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001000003_c
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_C1001000003_c = I_ERI_F3x_S_F3x_S_C1001000003_c+CDX*I_ERI_F3x_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2x_Px_C1001000003_c = I_ERI_F2xy_S_F3x_S_C1001000003_c+CDX*I_ERI_F2xy_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2x_Px_C1001000003_c = I_ERI_F2xz_S_F3x_S_C1001000003_c+CDX*I_ERI_F2xz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2x_Px_C1001000003_c = I_ERI_Fx2y_S_F3x_S_C1001000003_c+CDX*I_ERI_Fx2y_S_D2x_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2x_Px_C1001000003_c = I_ERI_Fxyz_S_F3x_S_C1001000003_c+CDX*I_ERI_Fxyz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2x_Px_C1001000003_c = I_ERI_Fx2z_S_F3x_S_C1001000003_c+CDX*I_ERI_Fx2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3y_S_D2x_Px_C1001000003_c = I_ERI_F3y_S_F3x_S_C1001000003_c+CDX*I_ERI_F3y_S_D2x_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2x_Px_C1001000003_c = I_ERI_F2yz_S_F3x_S_C1001000003_c+CDX*I_ERI_F2yz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2x_Px_C1001000003_c = I_ERI_Fy2z_S_F3x_S_C1001000003_c+CDX*I_ERI_Fy2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3z_S_D2x_Px_C1001000003_c = I_ERI_F3z_S_F3x_S_C1001000003_c+CDX*I_ERI_F3z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxy_Px_C1001000003_c = I_ERI_F3x_S_F2xy_S_C1001000003_c+CDX*I_ERI_F3x_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxy_Px_C1001000003_c = I_ERI_F2xy_S_F2xy_S_C1001000003_c+CDX*I_ERI_F2xy_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxy_Px_C1001000003_c = I_ERI_F2xz_S_F2xy_S_C1001000003_c+CDX*I_ERI_F2xz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxy_Px_C1001000003_c = I_ERI_Fx2y_S_F2xy_S_C1001000003_c+CDX*I_ERI_Fx2y_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxy_Px_C1001000003_c = I_ERI_Fxyz_S_F2xy_S_C1001000003_c+CDX*I_ERI_Fxyz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxy_Px_C1001000003_c = I_ERI_Fx2z_S_F2xy_S_C1001000003_c+CDX*I_ERI_Fx2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxy_Px_C1001000003_c = I_ERI_F3y_S_F2xy_S_C1001000003_c+CDX*I_ERI_F3y_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxy_Px_C1001000003_c = I_ERI_F2yz_S_F2xy_S_C1001000003_c+CDX*I_ERI_F2yz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxy_Px_C1001000003_c = I_ERI_Fy2z_S_F2xy_S_C1001000003_c+CDX*I_ERI_Fy2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxy_Px_C1001000003_c = I_ERI_F3z_S_F2xy_S_C1001000003_c+CDX*I_ERI_F3z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxz_Px_C1001000003_c = I_ERI_F3x_S_F2xz_S_C1001000003_c+CDX*I_ERI_F3x_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxz_Px_C1001000003_c = I_ERI_F2xy_S_F2xz_S_C1001000003_c+CDX*I_ERI_F2xy_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxz_Px_C1001000003_c = I_ERI_F2xz_S_F2xz_S_C1001000003_c+CDX*I_ERI_F2xz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxz_Px_C1001000003_c = I_ERI_Fx2y_S_F2xz_S_C1001000003_c+CDX*I_ERI_Fx2y_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxz_Px_C1001000003_c = I_ERI_Fxyz_S_F2xz_S_C1001000003_c+CDX*I_ERI_Fxyz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxz_Px_C1001000003_c = I_ERI_Fx2z_S_F2xz_S_C1001000003_c+CDX*I_ERI_Fx2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxz_Px_C1001000003_c = I_ERI_F3y_S_F2xz_S_C1001000003_c+CDX*I_ERI_F3y_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxz_Px_C1001000003_c = I_ERI_F2yz_S_F2xz_S_C1001000003_c+CDX*I_ERI_F2yz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxz_Px_C1001000003_c = I_ERI_Fy2z_S_F2xz_S_C1001000003_c+CDX*I_ERI_Fy2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxz_Px_C1001000003_c = I_ERI_F3z_S_F2xz_S_C1001000003_c+CDX*I_ERI_F3z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2y_Px_C1001000003_c = I_ERI_F3x_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F3x_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2y_Px_C1001000003_c = I_ERI_F2xy_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F2xy_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2y_Px_C1001000003_c = I_ERI_F2xz_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F2xz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2y_Px_C1001000003_c = I_ERI_Fx2y_S_Fx2y_S_C1001000003_c+CDX*I_ERI_Fx2y_S_D2y_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2y_Px_C1001000003_c = I_ERI_Fxyz_S_Fx2y_S_C1001000003_c+CDX*I_ERI_Fxyz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2y_Px_C1001000003_c = I_ERI_Fx2z_S_Fx2y_S_C1001000003_c+CDX*I_ERI_Fx2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3y_S_D2y_Px_C1001000003_c = I_ERI_F3y_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F3y_S_D2y_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2y_Px_C1001000003_c = I_ERI_F2yz_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F2yz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2y_Px_C1001000003_c = I_ERI_Fy2z_S_Fx2y_S_C1001000003_c+CDX*I_ERI_Fy2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3z_S_D2y_Px_C1001000003_c = I_ERI_F3z_S_Fx2y_S_C1001000003_c+CDX*I_ERI_F3z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3x_S_Dyz_Px_C1001000003_c = I_ERI_F3x_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F3x_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dyz_Px_C1001000003_c = I_ERI_F2xy_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F2xy_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dyz_Px_C1001000003_c = I_ERI_F2xz_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F2xz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dyz_Px_C1001000003_c = I_ERI_Fx2y_S_Fxyz_S_C1001000003_c+CDX*I_ERI_Fx2y_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dyz_Px_C1001000003_c = I_ERI_Fxyz_S_Fxyz_S_C1001000003_c+CDX*I_ERI_Fxyz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dyz_Px_C1001000003_c = I_ERI_Fx2z_S_Fxyz_S_C1001000003_c+CDX*I_ERI_Fx2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dyz_Px_C1001000003_c = I_ERI_F3y_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F3y_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dyz_Px_C1001000003_c = I_ERI_F2yz_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F2yz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dyz_Px_C1001000003_c = I_ERI_Fy2z_S_Fxyz_S_C1001000003_c+CDX*I_ERI_Fy2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dyz_Px_C1001000003_c = I_ERI_F3z_S_Fxyz_S_C1001000003_c+CDX*I_ERI_F3z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2z_Px_C1001000003_c = I_ERI_F3x_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F3x_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2z_Px_C1001000003_c = I_ERI_F2xy_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F2xy_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2z_Px_C1001000003_c = I_ERI_F2xz_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F2xz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2z_Px_C1001000003_c = I_ERI_Fx2y_S_Fx2z_S_C1001000003_c+CDX*I_ERI_Fx2y_S_D2z_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2z_Px_C1001000003_c = I_ERI_Fxyz_S_Fx2z_S_C1001000003_c+CDX*I_ERI_Fxyz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2z_Px_C1001000003_c = I_ERI_Fx2z_S_Fx2z_S_C1001000003_c+CDX*I_ERI_Fx2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3y_S_D2z_Px_C1001000003_c = I_ERI_F3y_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F3y_S_D2z_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2z_Px_C1001000003_c = I_ERI_F2yz_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F2yz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2z_Px_C1001000003_c = I_ERI_Fy2z_S_Fx2z_S_C1001000003_c+CDX*I_ERI_Fy2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3z_S_D2z_Px_C1001000003_c = I_ERI_F3z_S_Fx2z_S_C1001000003_c+CDX*I_ERI_F3z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3x_S_D2x_Py_C1001000003_c = I_ERI_F3x_S_F2xy_S_C1001000003_c+CDY*I_ERI_F3x_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2x_Py_C1001000003_c = I_ERI_F2xy_S_F2xy_S_C1001000003_c+CDY*I_ERI_F2xy_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2x_Py_C1001000003_c = I_ERI_F2xz_S_F2xy_S_C1001000003_c+CDY*I_ERI_F2xz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2x_Py_C1001000003_c = I_ERI_Fx2y_S_F2xy_S_C1001000003_c+CDY*I_ERI_Fx2y_S_D2x_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2x_Py_C1001000003_c = I_ERI_Fxyz_S_F2xy_S_C1001000003_c+CDY*I_ERI_Fxyz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2x_Py_C1001000003_c = I_ERI_Fx2z_S_F2xy_S_C1001000003_c+CDY*I_ERI_Fx2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3y_S_D2x_Py_C1001000003_c = I_ERI_F3y_S_F2xy_S_C1001000003_c+CDY*I_ERI_F3y_S_D2x_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2x_Py_C1001000003_c = I_ERI_F2yz_S_F2xy_S_C1001000003_c+CDY*I_ERI_F2yz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2x_Py_C1001000003_c = I_ERI_Fy2z_S_F2xy_S_C1001000003_c+CDY*I_ERI_Fy2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3z_S_D2x_Py_C1001000003_c = I_ERI_F3z_S_F2xy_S_C1001000003_c+CDY*I_ERI_F3z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxy_Py_C1001000003_c = I_ERI_F3x_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F3x_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxy_Py_C1001000003_c = I_ERI_F2xy_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F2xy_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxy_Py_C1001000003_c = I_ERI_F2xz_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F2xz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxy_Py_C1001000003_c = I_ERI_Fx2y_S_Fx2y_S_C1001000003_c+CDY*I_ERI_Fx2y_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxy_Py_C1001000003_c = I_ERI_Fxyz_S_Fx2y_S_C1001000003_c+CDY*I_ERI_Fxyz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxy_Py_C1001000003_c = I_ERI_Fx2z_S_Fx2y_S_C1001000003_c+CDY*I_ERI_Fx2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxy_Py_C1001000003_c = I_ERI_F3y_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F3y_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxy_Py_C1001000003_c = I_ERI_F2yz_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F2yz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxy_Py_C1001000003_c = I_ERI_Fy2z_S_Fx2y_S_C1001000003_c+CDY*I_ERI_Fy2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxy_Py_C1001000003_c = I_ERI_F3z_S_Fx2y_S_C1001000003_c+CDY*I_ERI_F3z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxz_Py_C1001000003_c = I_ERI_F3x_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F3x_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxz_Py_C1001000003_c = I_ERI_F2xy_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F2xy_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxz_Py_C1001000003_c = I_ERI_F2xz_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F2xz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxz_Py_C1001000003_c = I_ERI_Fx2y_S_Fxyz_S_C1001000003_c+CDY*I_ERI_Fx2y_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxz_Py_C1001000003_c = I_ERI_Fxyz_S_Fxyz_S_C1001000003_c+CDY*I_ERI_Fxyz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxz_Py_C1001000003_c = I_ERI_Fx2z_S_Fxyz_S_C1001000003_c+CDY*I_ERI_Fx2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxz_Py_C1001000003_c = I_ERI_F3y_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F3y_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxz_Py_C1001000003_c = I_ERI_F2yz_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F2yz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxz_Py_C1001000003_c = I_ERI_Fy2z_S_Fxyz_S_C1001000003_c+CDY*I_ERI_Fy2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxz_Py_C1001000003_c = I_ERI_F3z_S_Fxyz_S_C1001000003_c+CDY*I_ERI_F3z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2y_Py_C1001000003_c = I_ERI_F3x_S_F3y_S_C1001000003_c+CDY*I_ERI_F3x_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2y_Py_C1001000003_c = I_ERI_F2xy_S_F3y_S_C1001000003_c+CDY*I_ERI_F2xy_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2y_Py_C1001000003_c = I_ERI_F2xz_S_F3y_S_C1001000003_c+CDY*I_ERI_F2xz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2y_Py_C1001000003_c = I_ERI_Fx2y_S_F3y_S_C1001000003_c+CDY*I_ERI_Fx2y_S_D2y_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2y_Py_C1001000003_c = I_ERI_Fxyz_S_F3y_S_C1001000003_c+CDY*I_ERI_Fxyz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2y_Py_C1001000003_c = I_ERI_Fx2z_S_F3y_S_C1001000003_c+CDY*I_ERI_Fx2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3y_S_D2y_Py_C1001000003_c = I_ERI_F3y_S_F3y_S_C1001000003_c+CDY*I_ERI_F3y_S_D2y_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2y_Py_C1001000003_c = I_ERI_F2yz_S_F3y_S_C1001000003_c+CDY*I_ERI_F2yz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2y_Py_C1001000003_c = I_ERI_Fy2z_S_F3y_S_C1001000003_c+CDY*I_ERI_Fy2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3z_S_D2y_Py_C1001000003_c = I_ERI_F3z_S_F3y_S_C1001000003_c+CDY*I_ERI_F3z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3x_S_Dyz_Py_C1001000003_c = I_ERI_F3x_S_F2yz_S_C1001000003_c+CDY*I_ERI_F3x_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dyz_Py_C1001000003_c = I_ERI_F2xy_S_F2yz_S_C1001000003_c+CDY*I_ERI_F2xy_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dyz_Py_C1001000003_c = I_ERI_F2xz_S_F2yz_S_C1001000003_c+CDY*I_ERI_F2xz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dyz_Py_C1001000003_c = I_ERI_Fx2y_S_F2yz_S_C1001000003_c+CDY*I_ERI_Fx2y_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dyz_Py_C1001000003_c = I_ERI_Fxyz_S_F2yz_S_C1001000003_c+CDY*I_ERI_Fxyz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dyz_Py_C1001000003_c = I_ERI_Fx2z_S_F2yz_S_C1001000003_c+CDY*I_ERI_Fx2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dyz_Py_C1001000003_c = I_ERI_F3y_S_F2yz_S_C1001000003_c+CDY*I_ERI_F3y_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dyz_Py_C1001000003_c = I_ERI_F2yz_S_F2yz_S_C1001000003_c+CDY*I_ERI_F2yz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dyz_Py_C1001000003_c = I_ERI_Fy2z_S_F2yz_S_C1001000003_c+CDY*I_ERI_Fy2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dyz_Py_C1001000003_c = I_ERI_F3z_S_F2yz_S_C1001000003_c+CDY*I_ERI_F3z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2z_Py_C1001000003_c = I_ERI_F3x_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F3x_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2z_Py_C1001000003_c = I_ERI_F2xy_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F2xy_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2z_Py_C1001000003_c = I_ERI_F2xz_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F2xz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2z_Py_C1001000003_c = I_ERI_Fx2y_S_Fy2z_S_C1001000003_c+CDY*I_ERI_Fx2y_S_D2z_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2z_Py_C1001000003_c = I_ERI_Fxyz_S_Fy2z_S_C1001000003_c+CDY*I_ERI_Fxyz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2z_Py_C1001000003_c = I_ERI_Fx2z_S_Fy2z_S_C1001000003_c+CDY*I_ERI_Fx2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3y_S_D2z_Py_C1001000003_c = I_ERI_F3y_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F3y_S_D2z_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2z_Py_C1001000003_c = I_ERI_F2yz_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F2yz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2z_Py_C1001000003_c = I_ERI_Fy2z_S_Fy2z_S_C1001000003_c+CDY*I_ERI_Fy2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3z_S_D2z_Py_C1001000003_c = I_ERI_F3z_S_Fy2z_S_C1001000003_c+CDY*I_ERI_F3z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3x_S_D2x_Pz_C1001000003_c = I_ERI_F3x_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F3x_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2x_Pz_C1001000003_c = I_ERI_F2xy_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F2xy_S_D2x_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2x_Pz_C1001000003_c = I_ERI_F2xz_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F2xz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2x_Pz_C1001000003_c = I_ERI_Fx2y_S_F2xz_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_D2x_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2x_Pz_C1001000003_c = I_ERI_Fxyz_S_F2xz_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2x_Pz_C1001000003_c = I_ERI_Fx2z_S_F2xz_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3y_S_D2x_Pz_C1001000003_c = I_ERI_F3y_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F3y_S_D2x_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2x_Pz_C1001000003_c = I_ERI_F2yz_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F2yz_S_D2x_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2x_Pz_C1001000003_c = I_ERI_Fy2z_S_F2xz_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3z_S_D2x_Pz_C1001000003_c = I_ERI_F3z_S_F2xz_S_C1001000003_c+CDZ*I_ERI_F3z_S_D2x_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxy_Pz_C1001000003_c = I_ERI_F3x_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F3x_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxy_Pz_C1001000003_c = I_ERI_F2xy_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F2xy_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxy_Pz_C1001000003_c = I_ERI_F2xz_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F2xz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxy_Pz_C1001000003_c = I_ERI_Fx2y_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxy_Pz_C1001000003_c = I_ERI_Fxyz_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxy_Pz_C1001000003_c = I_ERI_Fx2z_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxy_Pz_C1001000003_c = I_ERI_F3y_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F3y_S_Dxy_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxy_Pz_C1001000003_c = I_ERI_F2yz_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F2yz_S_Dxy_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxy_Pz_C1001000003_c = I_ERI_Fy2z_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxy_Pz_C1001000003_c = I_ERI_F3z_S_Fxyz_S_C1001000003_c+CDZ*I_ERI_F3z_S_Dxy_S_C1001000003_c;
  Double I_ERI_F3x_S_Dxz_Pz_C1001000003_c = I_ERI_F3x_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F3x_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dxz_Pz_C1001000003_c = I_ERI_F2xy_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F2xy_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dxz_Pz_C1001000003_c = I_ERI_F2xz_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F2xz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dxz_Pz_C1001000003_c = I_ERI_Fx2y_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dxz_Pz_C1001000003_c = I_ERI_Fxyz_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dxz_Pz_C1001000003_c = I_ERI_Fx2z_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dxz_Pz_C1001000003_c = I_ERI_F3y_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F3y_S_Dxz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dxz_Pz_C1001000003_c = I_ERI_F2yz_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F2yz_S_Dxz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dxz_Pz_C1001000003_c = I_ERI_Fy2z_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dxz_Pz_C1001000003_c = I_ERI_F3z_S_Fx2z_S_C1001000003_c+CDZ*I_ERI_F3z_S_Dxz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2y_Pz_C1001000003_c = I_ERI_F3x_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F3x_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2y_Pz_C1001000003_c = I_ERI_F2xy_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F2xy_S_D2y_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2y_Pz_C1001000003_c = I_ERI_F2xz_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F2xz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2y_Pz_C1001000003_c = I_ERI_Fx2y_S_F2yz_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_D2y_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2y_Pz_C1001000003_c = I_ERI_Fxyz_S_F2yz_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2y_Pz_C1001000003_c = I_ERI_Fx2z_S_F2yz_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3y_S_D2y_Pz_C1001000003_c = I_ERI_F3y_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F3y_S_D2y_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2y_Pz_C1001000003_c = I_ERI_F2yz_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F2yz_S_D2y_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2y_Pz_C1001000003_c = I_ERI_Fy2z_S_F2yz_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3z_S_D2y_Pz_C1001000003_c = I_ERI_F3z_S_F2yz_S_C1001000003_c+CDZ*I_ERI_F3z_S_D2y_S_C1001000003_c;
  Double I_ERI_F3x_S_Dyz_Pz_C1001000003_c = I_ERI_F3x_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F3x_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xy_S_Dyz_Pz_C1001000003_c = I_ERI_F2xy_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F2xy_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2xz_S_Dyz_Pz_C1001000003_c = I_ERI_F2xz_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F2xz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2y_S_Dyz_Pz_C1001000003_c = I_ERI_Fx2y_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fxyz_S_Dyz_Pz_C1001000003_c = I_ERI_Fxyz_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fx2z_S_Dyz_Pz_C1001000003_c = I_ERI_Fx2z_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3y_S_Dyz_Pz_C1001000003_c = I_ERI_F3y_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F3y_S_Dyz_S_C1001000003_c;
  Double I_ERI_F2yz_S_Dyz_Pz_C1001000003_c = I_ERI_F2yz_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F2yz_S_Dyz_S_C1001000003_c;
  Double I_ERI_Fy2z_S_Dyz_Pz_C1001000003_c = I_ERI_Fy2z_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3z_S_Dyz_Pz_C1001000003_c = I_ERI_F3z_S_Fy2z_S_C1001000003_c+CDZ*I_ERI_F3z_S_Dyz_S_C1001000003_c;
  Double I_ERI_F3x_S_D2z_Pz_C1001000003_c = I_ERI_F3x_S_F3z_S_C1001000003_c+CDZ*I_ERI_F3x_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xy_S_D2z_Pz_C1001000003_c = I_ERI_F2xy_S_F3z_S_C1001000003_c+CDZ*I_ERI_F2xy_S_D2z_S_C1001000003_c;
  Double I_ERI_F2xz_S_D2z_Pz_C1001000003_c = I_ERI_F2xz_S_F3z_S_C1001000003_c+CDZ*I_ERI_F2xz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2y_S_D2z_Pz_C1001000003_c = I_ERI_Fx2y_S_F3z_S_C1001000003_c+CDZ*I_ERI_Fx2y_S_D2z_S_C1001000003_c;
  Double I_ERI_Fxyz_S_D2z_Pz_C1001000003_c = I_ERI_Fxyz_S_F3z_S_C1001000003_c+CDZ*I_ERI_Fxyz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fx2z_S_D2z_Pz_C1001000003_c = I_ERI_Fx2z_S_F3z_S_C1001000003_c+CDZ*I_ERI_Fx2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3y_S_D2z_Pz_C1001000003_c = I_ERI_F3y_S_F3z_S_C1001000003_c+CDZ*I_ERI_F3y_S_D2z_S_C1001000003_c;
  Double I_ERI_F2yz_S_D2z_Pz_C1001000003_c = I_ERI_F2yz_S_F3z_S_C1001000003_c+CDZ*I_ERI_F2yz_S_D2z_S_C1001000003_c;
  Double I_ERI_Fy2z_S_D2z_Pz_C1001000003_c = I_ERI_Fy2z_S_F3z_S_C1001000003_c+CDZ*I_ERI_Fy2z_S_D2z_S_C1001000003_c;
  Double I_ERI_F3z_S_D2z_Pz_C1001000003_c = I_ERI_F3z_S_F3z_S_C1001000003_c+CDZ*I_ERI_F3z_S_D2z_S_C1001000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1000003_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000003_d
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1000003_d = I_ERI_F3x_S_D2x_S_C1000003_d+CDX*I_ERI_F3x_S_Px_S_C1000003_d;
  Double I_ERI_F2xy_S_Px_Px_C1000003_d = I_ERI_F2xy_S_D2x_S_C1000003_d+CDX*I_ERI_F2xy_S_Px_S_C1000003_d;
  Double I_ERI_F2xz_S_Px_Px_C1000003_d = I_ERI_F2xz_S_D2x_S_C1000003_d+CDX*I_ERI_F2xz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2y_S_Px_Px_C1000003_d = I_ERI_Fx2y_S_D2x_S_C1000003_d+CDX*I_ERI_Fx2y_S_Px_S_C1000003_d;
  Double I_ERI_Fxyz_S_Px_Px_C1000003_d = I_ERI_Fxyz_S_D2x_S_C1000003_d+CDX*I_ERI_Fxyz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2z_S_Px_Px_C1000003_d = I_ERI_Fx2z_S_D2x_S_C1000003_d+CDX*I_ERI_Fx2z_S_Px_S_C1000003_d;
  Double I_ERI_F3y_S_Px_Px_C1000003_d = I_ERI_F3y_S_D2x_S_C1000003_d+CDX*I_ERI_F3y_S_Px_S_C1000003_d;
  Double I_ERI_F2yz_S_Px_Px_C1000003_d = I_ERI_F2yz_S_D2x_S_C1000003_d+CDX*I_ERI_F2yz_S_Px_S_C1000003_d;
  Double I_ERI_Fy2z_S_Px_Px_C1000003_d = I_ERI_Fy2z_S_D2x_S_C1000003_d+CDX*I_ERI_Fy2z_S_Px_S_C1000003_d;
  Double I_ERI_F3z_S_Px_Px_C1000003_d = I_ERI_F3z_S_D2x_S_C1000003_d+CDX*I_ERI_F3z_S_Px_S_C1000003_d;
  Double I_ERI_F3x_S_Py_Px_C1000003_d = I_ERI_F3x_S_Dxy_S_C1000003_d+CDX*I_ERI_F3x_S_Py_S_C1000003_d;
  Double I_ERI_F2xy_S_Py_Px_C1000003_d = I_ERI_F2xy_S_Dxy_S_C1000003_d+CDX*I_ERI_F2xy_S_Py_S_C1000003_d;
  Double I_ERI_F2xz_S_Py_Px_C1000003_d = I_ERI_F2xz_S_Dxy_S_C1000003_d+CDX*I_ERI_F2xz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2y_S_Py_Px_C1000003_d = I_ERI_Fx2y_S_Dxy_S_C1000003_d+CDX*I_ERI_Fx2y_S_Py_S_C1000003_d;
  Double I_ERI_Fxyz_S_Py_Px_C1000003_d = I_ERI_Fxyz_S_Dxy_S_C1000003_d+CDX*I_ERI_Fxyz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2z_S_Py_Px_C1000003_d = I_ERI_Fx2z_S_Dxy_S_C1000003_d+CDX*I_ERI_Fx2z_S_Py_S_C1000003_d;
  Double I_ERI_F3y_S_Py_Px_C1000003_d = I_ERI_F3y_S_Dxy_S_C1000003_d+CDX*I_ERI_F3y_S_Py_S_C1000003_d;
  Double I_ERI_F2yz_S_Py_Px_C1000003_d = I_ERI_F2yz_S_Dxy_S_C1000003_d+CDX*I_ERI_F2yz_S_Py_S_C1000003_d;
  Double I_ERI_Fy2z_S_Py_Px_C1000003_d = I_ERI_Fy2z_S_Dxy_S_C1000003_d+CDX*I_ERI_Fy2z_S_Py_S_C1000003_d;
  Double I_ERI_F3z_S_Py_Px_C1000003_d = I_ERI_F3z_S_Dxy_S_C1000003_d+CDX*I_ERI_F3z_S_Py_S_C1000003_d;
  Double I_ERI_F3x_S_Pz_Px_C1000003_d = I_ERI_F3x_S_Dxz_S_C1000003_d+CDX*I_ERI_F3x_S_Pz_S_C1000003_d;
  Double I_ERI_F2xy_S_Pz_Px_C1000003_d = I_ERI_F2xy_S_Dxz_S_C1000003_d+CDX*I_ERI_F2xy_S_Pz_S_C1000003_d;
  Double I_ERI_F2xz_S_Pz_Px_C1000003_d = I_ERI_F2xz_S_Dxz_S_C1000003_d+CDX*I_ERI_F2xz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2y_S_Pz_Px_C1000003_d = I_ERI_Fx2y_S_Dxz_S_C1000003_d+CDX*I_ERI_Fx2y_S_Pz_S_C1000003_d;
  Double I_ERI_Fxyz_S_Pz_Px_C1000003_d = I_ERI_Fxyz_S_Dxz_S_C1000003_d+CDX*I_ERI_Fxyz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2z_S_Pz_Px_C1000003_d = I_ERI_Fx2z_S_Dxz_S_C1000003_d+CDX*I_ERI_Fx2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3y_S_Pz_Px_C1000003_d = I_ERI_F3y_S_Dxz_S_C1000003_d+CDX*I_ERI_F3y_S_Pz_S_C1000003_d;
  Double I_ERI_F2yz_S_Pz_Px_C1000003_d = I_ERI_F2yz_S_Dxz_S_C1000003_d+CDX*I_ERI_F2yz_S_Pz_S_C1000003_d;
  Double I_ERI_Fy2z_S_Pz_Px_C1000003_d = I_ERI_Fy2z_S_Dxz_S_C1000003_d+CDX*I_ERI_Fy2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3z_S_Pz_Px_C1000003_d = I_ERI_F3z_S_Dxz_S_C1000003_d+CDX*I_ERI_F3z_S_Pz_S_C1000003_d;
  Double I_ERI_F3x_S_Px_Py_C1000003_d = I_ERI_F3x_S_Dxy_S_C1000003_d+CDY*I_ERI_F3x_S_Px_S_C1000003_d;
  Double I_ERI_F2xy_S_Px_Py_C1000003_d = I_ERI_F2xy_S_Dxy_S_C1000003_d+CDY*I_ERI_F2xy_S_Px_S_C1000003_d;
  Double I_ERI_F2xz_S_Px_Py_C1000003_d = I_ERI_F2xz_S_Dxy_S_C1000003_d+CDY*I_ERI_F2xz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2y_S_Px_Py_C1000003_d = I_ERI_Fx2y_S_Dxy_S_C1000003_d+CDY*I_ERI_Fx2y_S_Px_S_C1000003_d;
  Double I_ERI_Fxyz_S_Px_Py_C1000003_d = I_ERI_Fxyz_S_Dxy_S_C1000003_d+CDY*I_ERI_Fxyz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2z_S_Px_Py_C1000003_d = I_ERI_Fx2z_S_Dxy_S_C1000003_d+CDY*I_ERI_Fx2z_S_Px_S_C1000003_d;
  Double I_ERI_F3y_S_Px_Py_C1000003_d = I_ERI_F3y_S_Dxy_S_C1000003_d+CDY*I_ERI_F3y_S_Px_S_C1000003_d;
  Double I_ERI_F2yz_S_Px_Py_C1000003_d = I_ERI_F2yz_S_Dxy_S_C1000003_d+CDY*I_ERI_F2yz_S_Px_S_C1000003_d;
  Double I_ERI_Fy2z_S_Px_Py_C1000003_d = I_ERI_Fy2z_S_Dxy_S_C1000003_d+CDY*I_ERI_Fy2z_S_Px_S_C1000003_d;
  Double I_ERI_F3z_S_Px_Py_C1000003_d = I_ERI_F3z_S_Dxy_S_C1000003_d+CDY*I_ERI_F3z_S_Px_S_C1000003_d;
  Double I_ERI_F3x_S_Py_Py_C1000003_d = I_ERI_F3x_S_D2y_S_C1000003_d+CDY*I_ERI_F3x_S_Py_S_C1000003_d;
  Double I_ERI_F2xy_S_Py_Py_C1000003_d = I_ERI_F2xy_S_D2y_S_C1000003_d+CDY*I_ERI_F2xy_S_Py_S_C1000003_d;
  Double I_ERI_F2xz_S_Py_Py_C1000003_d = I_ERI_F2xz_S_D2y_S_C1000003_d+CDY*I_ERI_F2xz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2y_S_Py_Py_C1000003_d = I_ERI_Fx2y_S_D2y_S_C1000003_d+CDY*I_ERI_Fx2y_S_Py_S_C1000003_d;
  Double I_ERI_Fxyz_S_Py_Py_C1000003_d = I_ERI_Fxyz_S_D2y_S_C1000003_d+CDY*I_ERI_Fxyz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2z_S_Py_Py_C1000003_d = I_ERI_Fx2z_S_D2y_S_C1000003_d+CDY*I_ERI_Fx2z_S_Py_S_C1000003_d;
  Double I_ERI_F3y_S_Py_Py_C1000003_d = I_ERI_F3y_S_D2y_S_C1000003_d+CDY*I_ERI_F3y_S_Py_S_C1000003_d;
  Double I_ERI_F2yz_S_Py_Py_C1000003_d = I_ERI_F2yz_S_D2y_S_C1000003_d+CDY*I_ERI_F2yz_S_Py_S_C1000003_d;
  Double I_ERI_Fy2z_S_Py_Py_C1000003_d = I_ERI_Fy2z_S_D2y_S_C1000003_d+CDY*I_ERI_Fy2z_S_Py_S_C1000003_d;
  Double I_ERI_F3z_S_Py_Py_C1000003_d = I_ERI_F3z_S_D2y_S_C1000003_d+CDY*I_ERI_F3z_S_Py_S_C1000003_d;
  Double I_ERI_F3x_S_Pz_Py_C1000003_d = I_ERI_F3x_S_Dyz_S_C1000003_d+CDY*I_ERI_F3x_S_Pz_S_C1000003_d;
  Double I_ERI_F2xy_S_Pz_Py_C1000003_d = I_ERI_F2xy_S_Dyz_S_C1000003_d+CDY*I_ERI_F2xy_S_Pz_S_C1000003_d;
  Double I_ERI_F2xz_S_Pz_Py_C1000003_d = I_ERI_F2xz_S_Dyz_S_C1000003_d+CDY*I_ERI_F2xz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2y_S_Pz_Py_C1000003_d = I_ERI_Fx2y_S_Dyz_S_C1000003_d+CDY*I_ERI_Fx2y_S_Pz_S_C1000003_d;
  Double I_ERI_Fxyz_S_Pz_Py_C1000003_d = I_ERI_Fxyz_S_Dyz_S_C1000003_d+CDY*I_ERI_Fxyz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2z_S_Pz_Py_C1000003_d = I_ERI_Fx2z_S_Dyz_S_C1000003_d+CDY*I_ERI_Fx2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3y_S_Pz_Py_C1000003_d = I_ERI_F3y_S_Dyz_S_C1000003_d+CDY*I_ERI_F3y_S_Pz_S_C1000003_d;
  Double I_ERI_F2yz_S_Pz_Py_C1000003_d = I_ERI_F2yz_S_Dyz_S_C1000003_d+CDY*I_ERI_F2yz_S_Pz_S_C1000003_d;
  Double I_ERI_Fy2z_S_Pz_Py_C1000003_d = I_ERI_Fy2z_S_Dyz_S_C1000003_d+CDY*I_ERI_Fy2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3z_S_Pz_Py_C1000003_d = I_ERI_F3z_S_Dyz_S_C1000003_d+CDY*I_ERI_F3z_S_Pz_S_C1000003_d;
  Double I_ERI_F3x_S_Px_Pz_C1000003_d = I_ERI_F3x_S_Dxz_S_C1000003_d+CDZ*I_ERI_F3x_S_Px_S_C1000003_d;
  Double I_ERI_F2xy_S_Px_Pz_C1000003_d = I_ERI_F2xy_S_Dxz_S_C1000003_d+CDZ*I_ERI_F2xy_S_Px_S_C1000003_d;
  Double I_ERI_F2xz_S_Px_Pz_C1000003_d = I_ERI_F2xz_S_Dxz_S_C1000003_d+CDZ*I_ERI_F2xz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2y_S_Px_Pz_C1000003_d = I_ERI_Fx2y_S_Dxz_S_C1000003_d+CDZ*I_ERI_Fx2y_S_Px_S_C1000003_d;
  Double I_ERI_Fxyz_S_Px_Pz_C1000003_d = I_ERI_Fxyz_S_Dxz_S_C1000003_d+CDZ*I_ERI_Fxyz_S_Px_S_C1000003_d;
  Double I_ERI_Fx2z_S_Px_Pz_C1000003_d = I_ERI_Fx2z_S_Dxz_S_C1000003_d+CDZ*I_ERI_Fx2z_S_Px_S_C1000003_d;
  Double I_ERI_F3y_S_Px_Pz_C1000003_d = I_ERI_F3y_S_Dxz_S_C1000003_d+CDZ*I_ERI_F3y_S_Px_S_C1000003_d;
  Double I_ERI_F2yz_S_Px_Pz_C1000003_d = I_ERI_F2yz_S_Dxz_S_C1000003_d+CDZ*I_ERI_F2yz_S_Px_S_C1000003_d;
  Double I_ERI_Fy2z_S_Px_Pz_C1000003_d = I_ERI_Fy2z_S_Dxz_S_C1000003_d+CDZ*I_ERI_Fy2z_S_Px_S_C1000003_d;
  Double I_ERI_F3z_S_Px_Pz_C1000003_d = I_ERI_F3z_S_Dxz_S_C1000003_d+CDZ*I_ERI_F3z_S_Px_S_C1000003_d;
  Double I_ERI_F3x_S_Py_Pz_C1000003_d = I_ERI_F3x_S_Dyz_S_C1000003_d+CDZ*I_ERI_F3x_S_Py_S_C1000003_d;
  Double I_ERI_F2xy_S_Py_Pz_C1000003_d = I_ERI_F2xy_S_Dyz_S_C1000003_d+CDZ*I_ERI_F2xy_S_Py_S_C1000003_d;
  Double I_ERI_F2xz_S_Py_Pz_C1000003_d = I_ERI_F2xz_S_Dyz_S_C1000003_d+CDZ*I_ERI_F2xz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2y_S_Py_Pz_C1000003_d = I_ERI_Fx2y_S_Dyz_S_C1000003_d+CDZ*I_ERI_Fx2y_S_Py_S_C1000003_d;
  Double I_ERI_Fxyz_S_Py_Pz_C1000003_d = I_ERI_Fxyz_S_Dyz_S_C1000003_d+CDZ*I_ERI_Fxyz_S_Py_S_C1000003_d;
  Double I_ERI_Fx2z_S_Py_Pz_C1000003_d = I_ERI_Fx2z_S_Dyz_S_C1000003_d+CDZ*I_ERI_Fx2z_S_Py_S_C1000003_d;
  Double I_ERI_F3y_S_Py_Pz_C1000003_d = I_ERI_F3y_S_Dyz_S_C1000003_d+CDZ*I_ERI_F3y_S_Py_S_C1000003_d;
  Double I_ERI_F2yz_S_Py_Pz_C1000003_d = I_ERI_F2yz_S_Dyz_S_C1000003_d+CDZ*I_ERI_F2yz_S_Py_S_C1000003_d;
  Double I_ERI_Fy2z_S_Py_Pz_C1000003_d = I_ERI_Fy2z_S_Dyz_S_C1000003_d+CDZ*I_ERI_Fy2z_S_Py_S_C1000003_d;
  Double I_ERI_F3z_S_Py_Pz_C1000003_d = I_ERI_F3z_S_Dyz_S_C1000003_d+CDZ*I_ERI_F3z_S_Py_S_C1000003_d;
  Double I_ERI_F3x_S_Pz_Pz_C1000003_d = I_ERI_F3x_S_D2z_S_C1000003_d+CDZ*I_ERI_F3x_S_Pz_S_C1000003_d;
  Double I_ERI_F2xy_S_Pz_Pz_C1000003_d = I_ERI_F2xy_S_D2z_S_C1000003_d+CDZ*I_ERI_F2xy_S_Pz_S_C1000003_d;
  Double I_ERI_F2xz_S_Pz_Pz_C1000003_d = I_ERI_F2xz_S_D2z_S_C1000003_d+CDZ*I_ERI_F2xz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2y_S_Pz_Pz_C1000003_d = I_ERI_Fx2y_S_D2z_S_C1000003_d+CDZ*I_ERI_Fx2y_S_Pz_S_C1000003_d;
  Double I_ERI_Fxyz_S_Pz_Pz_C1000003_d = I_ERI_Fxyz_S_D2z_S_C1000003_d+CDZ*I_ERI_Fxyz_S_Pz_S_C1000003_d;
  Double I_ERI_Fx2z_S_Pz_Pz_C1000003_d = I_ERI_Fx2z_S_D2z_S_C1000003_d+CDZ*I_ERI_Fx2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3y_S_Pz_Pz_C1000003_d = I_ERI_F3y_S_D2z_S_C1000003_d+CDZ*I_ERI_F3y_S_Pz_S_C1000003_d;
  Double I_ERI_F2yz_S_Pz_Pz_C1000003_d = I_ERI_F2yz_S_D2z_S_C1000003_d+CDZ*I_ERI_F2yz_S_Pz_S_C1000003_d;
  Double I_ERI_Fy2z_S_Pz_Pz_C1000003_d = I_ERI_Fy2z_S_D2z_S_C1000003_d+CDZ*I_ERI_Fy2z_S_Pz_S_C1000003_d;
  Double I_ERI_F3z_S_Pz_Pz_C1000003_d = I_ERI_F3z_S_D2z_S_C1000003_d+CDZ*I_ERI_F3z_S_Pz_S_C1000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001000003_d
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1001000003_d = I_ERI_F3x_S_D2x_S_C1001000003_d+CDX*I_ERI_F3x_S_Px_S_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Px_C1001000003_d = I_ERI_F2xy_S_D2x_S_C1001000003_d+CDX*I_ERI_F2xy_S_Px_S_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Px_C1001000003_d = I_ERI_F2xz_S_D2x_S_C1001000003_d+CDX*I_ERI_F2xz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Px_C1001000003_d = I_ERI_Fx2y_S_D2x_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Px_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Px_C1001000003_d = I_ERI_Fxyz_S_D2x_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Px_C1001000003_d = I_ERI_Fx2z_S_D2x_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3y_S_Px_Px_C1001000003_d = I_ERI_F3y_S_D2x_S_C1001000003_d+CDX*I_ERI_F3y_S_Px_S_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Px_C1001000003_d = I_ERI_F2yz_S_D2x_S_C1001000003_d+CDX*I_ERI_F2yz_S_Px_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Px_C1001000003_d = I_ERI_Fy2z_S_D2x_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3z_S_Px_Px_C1001000003_d = I_ERI_F3z_S_D2x_S_C1001000003_d+CDX*I_ERI_F3z_S_Px_S_C1001000003_d;
  Double I_ERI_F3x_S_Py_Px_C1001000003_d = I_ERI_F3x_S_Dxy_S_C1001000003_d+CDX*I_ERI_F3x_S_Py_S_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Px_C1001000003_d = I_ERI_F2xy_S_Dxy_S_C1001000003_d+CDX*I_ERI_F2xy_S_Py_S_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Px_C1001000003_d = I_ERI_F2xz_S_Dxy_S_C1001000003_d+CDX*I_ERI_F2xz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Px_C1001000003_d = I_ERI_Fx2y_S_Dxy_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Py_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Px_C1001000003_d = I_ERI_Fxyz_S_Dxy_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Px_C1001000003_d = I_ERI_Fx2z_S_Dxy_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3y_S_Py_Px_C1001000003_d = I_ERI_F3y_S_Dxy_S_C1001000003_d+CDX*I_ERI_F3y_S_Py_S_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Px_C1001000003_d = I_ERI_F2yz_S_Dxy_S_C1001000003_d+CDX*I_ERI_F2yz_S_Py_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Px_C1001000003_d = I_ERI_Fy2z_S_Dxy_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3z_S_Py_Px_C1001000003_d = I_ERI_F3z_S_Dxy_S_C1001000003_d+CDX*I_ERI_F3z_S_Py_S_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Px_C1001000003_d = I_ERI_F3x_S_Dxz_S_C1001000003_d+CDX*I_ERI_F3x_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Px_C1001000003_d = I_ERI_F2xy_S_Dxz_S_C1001000003_d+CDX*I_ERI_F2xy_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Px_C1001000003_d = I_ERI_F2xz_S_Dxz_S_C1001000003_d+CDX*I_ERI_F2xz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Px_C1001000003_d = I_ERI_Fx2y_S_Dxz_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Pz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Px_C1001000003_d = I_ERI_Fxyz_S_Dxz_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Px_C1001000003_d = I_ERI_Fx2z_S_Dxz_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Px_C1001000003_d = I_ERI_F3y_S_Dxz_S_C1001000003_d+CDX*I_ERI_F3y_S_Pz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Px_C1001000003_d = I_ERI_F2yz_S_Dxz_S_C1001000003_d+CDX*I_ERI_F2yz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Px_C1001000003_d = I_ERI_Fy2z_S_Dxz_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Px_C1001000003_d = I_ERI_F3z_S_Dxz_S_C1001000003_d+CDX*I_ERI_F3z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3x_S_Px_Py_C1001000003_d = I_ERI_F3x_S_Dxy_S_C1001000003_d+CDY*I_ERI_F3x_S_Px_S_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Py_C1001000003_d = I_ERI_F2xy_S_Dxy_S_C1001000003_d+CDY*I_ERI_F2xy_S_Px_S_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Py_C1001000003_d = I_ERI_F2xz_S_Dxy_S_C1001000003_d+CDY*I_ERI_F2xz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Py_C1001000003_d = I_ERI_Fx2y_S_Dxy_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Px_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Py_C1001000003_d = I_ERI_Fxyz_S_Dxy_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Py_C1001000003_d = I_ERI_Fx2z_S_Dxy_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3y_S_Px_Py_C1001000003_d = I_ERI_F3y_S_Dxy_S_C1001000003_d+CDY*I_ERI_F3y_S_Px_S_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Py_C1001000003_d = I_ERI_F2yz_S_Dxy_S_C1001000003_d+CDY*I_ERI_F2yz_S_Px_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Py_C1001000003_d = I_ERI_Fy2z_S_Dxy_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3z_S_Px_Py_C1001000003_d = I_ERI_F3z_S_Dxy_S_C1001000003_d+CDY*I_ERI_F3z_S_Px_S_C1001000003_d;
  Double I_ERI_F3x_S_Py_Py_C1001000003_d = I_ERI_F3x_S_D2y_S_C1001000003_d+CDY*I_ERI_F3x_S_Py_S_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Py_C1001000003_d = I_ERI_F2xy_S_D2y_S_C1001000003_d+CDY*I_ERI_F2xy_S_Py_S_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Py_C1001000003_d = I_ERI_F2xz_S_D2y_S_C1001000003_d+CDY*I_ERI_F2xz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Py_C1001000003_d = I_ERI_Fx2y_S_D2y_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Py_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Py_C1001000003_d = I_ERI_Fxyz_S_D2y_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Py_C1001000003_d = I_ERI_Fx2z_S_D2y_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3y_S_Py_Py_C1001000003_d = I_ERI_F3y_S_D2y_S_C1001000003_d+CDY*I_ERI_F3y_S_Py_S_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Py_C1001000003_d = I_ERI_F2yz_S_D2y_S_C1001000003_d+CDY*I_ERI_F2yz_S_Py_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Py_C1001000003_d = I_ERI_Fy2z_S_D2y_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3z_S_Py_Py_C1001000003_d = I_ERI_F3z_S_D2y_S_C1001000003_d+CDY*I_ERI_F3z_S_Py_S_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Py_C1001000003_d = I_ERI_F3x_S_Dyz_S_C1001000003_d+CDY*I_ERI_F3x_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Py_C1001000003_d = I_ERI_F2xy_S_Dyz_S_C1001000003_d+CDY*I_ERI_F2xy_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Py_C1001000003_d = I_ERI_F2xz_S_Dyz_S_C1001000003_d+CDY*I_ERI_F2xz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Py_C1001000003_d = I_ERI_Fx2y_S_Dyz_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Pz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Py_C1001000003_d = I_ERI_Fxyz_S_Dyz_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Py_C1001000003_d = I_ERI_Fx2z_S_Dyz_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Py_C1001000003_d = I_ERI_F3y_S_Dyz_S_C1001000003_d+CDY*I_ERI_F3y_S_Pz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Py_C1001000003_d = I_ERI_F2yz_S_Dyz_S_C1001000003_d+CDY*I_ERI_F2yz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Py_C1001000003_d = I_ERI_Fy2z_S_Dyz_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Py_C1001000003_d = I_ERI_F3z_S_Dyz_S_C1001000003_d+CDY*I_ERI_F3z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3x_S_Px_Pz_C1001000003_d = I_ERI_F3x_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F3x_S_Px_S_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Pz_C1001000003_d = I_ERI_F2xy_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F2xy_S_Px_S_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Pz_C1001000003_d = I_ERI_F2xz_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F2xz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Pz_C1001000003_d = I_ERI_Fx2y_S_Dxz_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_Px_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Pz_C1001000003_d = I_ERI_Fxyz_S_Dxz_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_Px_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Pz_C1001000003_d = I_ERI_Fx2z_S_Dxz_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3y_S_Px_Pz_C1001000003_d = I_ERI_F3y_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F3y_S_Px_S_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Pz_C1001000003_d = I_ERI_F2yz_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F2yz_S_Px_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Pz_C1001000003_d = I_ERI_Fy2z_S_Dxz_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_Px_S_C1001000003_d;
  Double I_ERI_F3z_S_Px_Pz_C1001000003_d = I_ERI_F3z_S_Dxz_S_C1001000003_d+CDZ*I_ERI_F3z_S_Px_S_C1001000003_d;
  Double I_ERI_F3x_S_Py_Pz_C1001000003_d = I_ERI_F3x_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F3x_S_Py_S_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Pz_C1001000003_d = I_ERI_F2xy_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F2xy_S_Py_S_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Pz_C1001000003_d = I_ERI_F2xz_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F2xz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Pz_C1001000003_d = I_ERI_Fx2y_S_Dyz_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_Py_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Pz_C1001000003_d = I_ERI_Fxyz_S_Dyz_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_Py_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Pz_C1001000003_d = I_ERI_Fx2z_S_Dyz_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3y_S_Py_Pz_C1001000003_d = I_ERI_F3y_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F3y_S_Py_S_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Pz_C1001000003_d = I_ERI_F2yz_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F2yz_S_Py_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Pz_C1001000003_d = I_ERI_Fy2z_S_Dyz_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_Py_S_C1001000003_d;
  Double I_ERI_F3z_S_Py_Pz_C1001000003_d = I_ERI_F3z_S_Dyz_S_C1001000003_d+CDZ*I_ERI_F3z_S_Py_S_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Pz_C1001000003_d = I_ERI_F3x_S_D2z_S_C1001000003_d+CDZ*I_ERI_F3x_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Pz_C1001000003_d = I_ERI_F2xy_S_D2z_S_C1001000003_d+CDZ*I_ERI_F2xy_S_Pz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Pz_C1001000003_d = I_ERI_F2xz_S_D2z_S_C1001000003_d+CDZ*I_ERI_F2xz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Pz_C1001000003_d = I_ERI_Fx2y_S_D2z_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_Pz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Pz_C1001000003_d = I_ERI_Fxyz_S_D2z_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Pz_C1001000003_d = I_ERI_Fx2z_S_D2z_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Pz_C1001000003_d = I_ERI_F3y_S_D2z_S_C1001000003_d+CDZ*I_ERI_F3y_S_Pz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Pz_C1001000003_d = I_ERI_F2yz_S_D2z_S_C1001000003_d+CDZ*I_ERI_F2yz_S_Pz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Pz_C1001000003_d = I_ERI_Fy2z_S_D2z_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_Pz_S_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Pz_C1001000003_d = I_ERI_F3z_S_D2z_S_C1001000003_d+CDZ*I_ERI_F3z_S_Pz_S_C1001000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_C1001000003_d
   * expanding position: KET2
   * code section is: HRR
   * totally 40 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001000003_d
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_C1001000003_d = I_ERI_F3x_S_F3x_S_C1001000003_d+CDX*I_ERI_F3x_S_D2x_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2x_Px_C1001000003_d = I_ERI_F2xy_S_F3x_S_C1001000003_d+CDX*I_ERI_F2xy_S_D2x_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2x_Px_C1001000003_d = I_ERI_F2xz_S_F3x_S_C1001000003_d+CDX*I_ERI_F2xz_S_D2x_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2x_Px_C1001000003_d = I_ERI_Fx2y_S_F3x_S_C1001000003_d+CDX*I_ERI_Fx2y_S_D2x_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2x_Px_C1001000003_d = I_ERI_Fxyz_S_F3x_S_C1001000003_d+CDX*I_ERI_Fxyz_S_D2x_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2x_Px_C1001000003_d = I_ERI_Fx2z_S_F3x_S_C1001000003_d+CDX*I_ERI_Fx2z_S_D2x_S_C1001000003_d;
  Double I_ERI_F3y_S_D2x_Px_C1001000003_d = I_ERI_F3y_S_F3x_S_C1001000003_d+CDX*I_ERI_F3y_S_D2x_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2x_Px_C1001000003_d = I_ERI_F2yz_S_F3x_S_C1001000003_d+CDX*I_ERI_F2yz_S_D2x_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2x_Px_C1001000003_d = I_ERI_Fy2z_S_F3x_S_C1001000003_d+CDX*I_ERI_Fy2z_S_D2x_S_C1001000003_d;
  Double I_ERI_F3z_S_D2x_Px_C1001000003_d = I_ERI_F3z_S_F3x_S_C1001000003_d+CDX*I_ERI_F3z_S_D2x_S_C1001000003_d;
  Double I_ERI_F3x_S_Dxy_Px_C1001000003_d = I_ERI_F3x_S_F2xy_S_C1001000003_d+CDX*I_ERI_F3x_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dxy_Px_C1001000003_d = I_ERI_F2xy_S_F2xy_S_C1001000003_d+CDX*I_ERI_F2xy_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dxy_Px_C1001000003_d = I_ERI_F2xz_S_F2xy_S_C1001000003_d+CDX*I_ERI_F2xz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dxy_Px_C1001000003_d = I_ERI_Fx2y_S_F2xy_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dxy_Px_C1001000003_d = I_ERI_Fxyz_S_F2xy_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dxy_Px_C1001000003_d = I_ERI_Fx2z_S_F2xy_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3y_S_Dxy_Px_C1001000003_d = I_ERI_F3y_S_F2xy_S_C1001000003_d+CDX*I_ERI_F3y_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dxy_Px_C1001000003_d = I_ERI_F2yz_S_F2xy_S_C1001000003_d+CDX*I_ERI_F2yz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dxy_Px_C1001000003_d = I_ERI_Fy2z_S_F2xy_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3z_S_Dxy_Px_C1001000003_d = I_ERI_F3z_S_F2xy_S_C1001000003_d+CDX*I_ERI_F3z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3x_S_Dxz_Px_C1001000003_d = I_ERI_F3x_S_F2xz_S_C1001000003_d+CDX*I_ERI_F3x_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dxz_Px_C1001000003_d = I_ERI_F2xy_S_F2xz_S_C1001000003_d+CDX*I_ERI_F2xy_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dxz_Px_C1001000003_d = I_ERI_F2xz_S_F2xz_S_C1001000003_d+CDX*I_ERI_F2xz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dxz_Px_C1001000003_d = I_ERI_Fx2y_S_F2xz_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dxz_Px_C1001000003_d = I_ERI_Fxyz_S_F2xz_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dxz_Px_C1001000003_d = I_ERI_Fx2z_S_F2xz_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dxz_Px_C1001000003_d = I_ERI_F3y_S_F2xz_S_C1001000003_d+CDX*I_ERI_F3y_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dxz_Px_C1001000003_d = I_ERI_F2yz_S_F2xz_S_C1001000003_d+CDX*I_ERI_F2yz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dxz_Px_C1001000003_d = I_ERI_Fy2z_S_F2xz_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dxz_Px_C1001000003_d = I_ERI_F3z_S_F2xz_S_C1001000003_d+CDX*I_ERI_F3z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3x_S_D2y_Px_C1001000003_d = I_ERI_F3x_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F3x_S_D2y_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2y_Px_C1001000003_d = I_ERI_F2xy_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F2xy_S_D2y_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2y_Px_C1001000003_d = I_ERI_F2xz_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F2xz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2y_Px_C1001000003_d = I_ERI_Fx2y_S_Fx2y_S_C1001000003_d+CDX*I_ERI_Fx2y_S_D2y_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2y_Px_C1001000003_d = I_ERI_Fxyz_S_Fx2y_S_C1001000003_d+CDX*I_ERI_Fxyz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2y_Px_C1001000003_d = I_ERI_Fx2z_S_Fx2y_S_C1001000003_d+CDX*I_ERI_Fx2z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3y_S_D2y_Px_C1001000003_d = I_ERI_F3y_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F3y_S_D2y_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2y_Px_C1001000003_d = I_ERI_F2yz_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F2yz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2y_Px_C1001000003_d = I_ERI_Fy2z_S_Fx2y_S_C1001000003_d+CDX*I_ERI_Fy2z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3z_S_D2y_Px_C1001000003_d = I_ERI_F3z_S_Fx2y_S_C1001000003_d+CDX*I_ERI_F3z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3x_S_Dyz_Px_C1001000003_d = I_ERI_F3x_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F3x_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dyz_Px_C1001000003_d = I_ERI_F2xy_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F2xy_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dyz_Px_C1001000003_d = I_ERI_F2xz_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F2xz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dyz_Px_C1001000003_d = I_ERI_Fx2y_S_Fxyz_S_C1001000003_d+CDX*I_ERI_Fx2y_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dyz_Px_C1001000003_d = I_ERI_Fxyz_S_Fxyz_S_C1001000003_d+CDX*I_ERI_Fxyz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dyz_Px_C1001000003_d = I_ERI_Fx2z_S_Fxyz_S_C1001000003_d+CDX*I_ERI_Fx2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dyz_Px_C1001000003_d = I_ERI_F3y_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F3y_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dyz_Px_C1001000003_d = I_ERI_F2yz_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F2yz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dyz_Px_C1001000003_d = I_ERI_Fy2z_S_Fxyz_S_C1001000003_d+CDX*I_ERI_Fy2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dyz_Px_C1001000003_d = I_ERI_F3z_S_Fxyz_S_C1001000003_d+CDX*I_ERI_F3z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3x_S_D2z_Px_C1001000003_d = I_ERI_F3x_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F3x_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2z_Px_C1001000003_d = I_ERI_F2xy_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F2xy_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2z_Px_C1001000003_d = I_ERI_F2xz_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F2xz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2z_Px_C1001000003_d = I_ERI_Fx2y_S_Fx2z_S_C1001000003_d+CDX*I_ERI_Fx2y_S_D2z_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2z_Px_C1001000003_d = I_ERI_Fxyz_S_Fx2z_S_C1001000003_d+CDX*I_ERI_Fxyz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2z_Px_C1001000003_d = I_ERI_Fx2z_S_Fx2z_S_C1001000003_d+CDX*I_ERI_Fx2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3y_S_D2z_Px_C1001000003_d = I_ERI_F3y_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F3y_S_D2z_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2z_Px_C1001000003_d = I_ERI_F2yz_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F2yz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2z_Px_C1001000003_d = I_ERI_Fy2z_S_Fx2z_S_C1001000003_d+CDX*I_ERI_Fy2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3z_S_D2z_Px_C1001000003_d = I_ERI_F3z_S_Fx2z_S_C1001000003_d+CDX*I_ERI_F3z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3x_S_Dxy_Py_C1001000003_d = I_ERI_F3x_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F3x_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dxy_Py_C1001000003_d = I_ERI_F2xy_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F2xy_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dxy_Py_C1001000003_d = I_ERI_F2xz_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F2xz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dxy_Py_C1001000003_d = I_ERI_Fx2y_S_Fx2y_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dxy_Py_C1001000003_d = I_ERI_Fxyz_S_Fx2y_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dxy_Py_C1001000003_d = I_ERI_Fx2z_S_Fx2y_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3y_S_Dxy_Py_C1001000003_d = I_ERI_F3y_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F3y_S_Dxy_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dxy_Py_C1001000003_d = I_ERI_F2yz_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F2yz_S_Dxy_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dxy_Py_C1001000003_d = I_ERI_Fy2z_S_Fx2y_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3z_S_Dxy_Py_C1001000003_d = I_ERI_F3z_S_Fx2y_S_C1001000003_d+CDY*I_ERI_F3z_S_Dxy_S_C1001000003_d;
  Double I_ERI_F3x_S_Dxz_Py_C1001000003_d = I_ERI_F3x_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F3x_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dxz_Py_C1001000003_d = I_ERI_F2xy_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F2xy_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dxz_Py_C1001000003_d = I_ERI_F2xz_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F2xz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dxz_Py_C1001000003_d = I_ERI_Fx2y_S_Fxyz_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dxz_Py_C1001000003_d = I_ERI_Fxyz_S_Fxyz_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dxz_Py_C1001000003_d = I_ERI_Fx2z_S_Fxyz_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dxz_Py_C1001000003_d = I_ERI_F3y_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F3y_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dxz_Py_C1001000003_d = I_ERI_F2yz_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F2yz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dxz_Py_C1001000003_d = I_ERI_Fy2z_S_Fxyz_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dxz_Py_C1001000003_d = I_ERI_F3z_S_Fxyz_S_C1001000003_d+CDY*I_ERI_F3z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3x_S_D2y_Py_C1001000003_d = I_ERI_F3x_S_F3y_S_C1001000003_d+CDY*I_ERI_F3x_S_D2y_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2y_Py_C1001000003_d = I_ERI_F2xy_S_F3y_S_C1001000003_d+CDY*I_ERI_F2xy_S_D2y_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2y_Py_C1001000003_d = I_ERI_F2xz_S_F3y_S_C1001000003_d+CDY*I_ERI_F2xz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2y_Py_C1001000003_d = I_ERI_Fx2y_S_F3y_S_C1001000003_d+CDY*I_ERI_Fx2y_S_D2y_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2y_Py_C1001000003_d = I_ERI_Fxyz_S_F3y_S_C1001000003_d+CDY*I_ERI_Fxyz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2y_Py_C1001000003_d = I_ERI_Fx2z_S_F3y_S_C1001000003_d+CDY*I_ERI_Fx2z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3y_S_D2y_Py_C1001000003_d = I_ERI_F3y_S_F3y_S_C1001000003_d+CDY*I_ERI_F3y_S_D2y_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2y_Py_C1001000003_d = I_ERI_F2yz_S_F3y_S_C1001000003_d+CDY*I_ERI_F2yz_S_D2y_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2y_Py_C1001000003_d = I_ERI_Fy2z_S_F3y_S_C1001000003_d+CDY*I_ERI_Fy2z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3z_S_D2y_Py_C1001000003_d = I_ERI_F3z_S_F3y_S_C1001000003_d+CDY*I_ERI_F3z_S_D2y_S_C1001000003_d;
  Double I_ERI_F3x_S_Dyz_Py_C1001000003_d = I_ERI_F3x_S_F2yz_S_C1001000003_d+CDY*I_ERI_F3x_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dyz_Py_C1001000003_d = I_ERI_F2xy_S_F2yz_S_C1001000003_d+CDY*I_ERI_F2xy_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dyz_Py_C1001000003_d = I_ERI_F2xz_S_F2yz_S_C1001000003_d+CDY*I_ERI_F2xz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dyz_Py_C1001000003_d = I_ERI_Fx2y_S_F2yz_S_C1001000003_d+CDY*I_ERI_Fx2y_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dyz_Py_C1001000003_d = I_ERI_Fxyz_S_F2yz_S_C1001000003_d+CDY*I_ERI_Fxyz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dyz_Py_C1001000003_d = I_ERI_Fx2z_S_F2yz_S_C1001000003_d+CDY*I_ERI_Fx2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dyz_Py_C1001000003_d = I_ERI_F3y_S_F2yz_S_C1001000003_d+CDY*I_ERI_F3y_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dyz_Py_C1001000003_d = I_ERI_F2yz_S_F2yz_S_C1001000003_d+CDY*I_ERI_F2yz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dyz_Py_C1001000003_d = I_ERI_Fy2z_S_F2yz_S_C1001000003_d+CDY*I_ERI_Fy2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dyz_Py_C1001000003_d = I_ERI_F3z_S_F2yz_S_C1001000003_d+CDY*I_ERI_F3z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3x_S_D2z_Py_C1001000003_d = I_ERI_F3x_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F3x_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2z_Py_C1001000003_d = I_ERI_F2xy_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F2xy_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2z_Py_C1001000003_d = I_ERI_F2xz_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F2xz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2z_Py_C1001000003_d = I_ERI_Fx2y_S_Fy2z_S_C1001000003_d+CDY*I_ERI_Fx2y_S_D2z_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2z_Py_C1001000003_d = I_ERI_Fxyz_S_Fy2z_S_C1001000003_d+CDY*I_ERI_Fxyz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2z_Py_C1001000003_d = I_ERI_Fx2z_S_Fy2z_S_C1001000003_d+CDY*I_ERI_Fx2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3y_S_D2z_Py_C1001000003_d = I_ERI_F3y_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F3y_S_D2z_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2z_Py_C1001000003_d = I_ERI_F2yz_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F2yz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2z_Py_C1001000003_d = I_ERI_Fy2z_S_Fy2z_S_C1001000003_d+CDY*I_ERI_Fy2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3z_S_D2z_Py_C1001000003_d = I_ERI_F3z_S_Fy2z_S_C1001000003_d+CDY*I_ERI_F3z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3x_S_Dxz_Pz_C1001000003_d = I_ERI_F3x_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F3x_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dxz_Pz_C1001000003_d = I_ERI_F2xy_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F2xy_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dxz_Pz_C1001000003_d = I_ERI_F2xz_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F2xz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dxz_Pz_C1001000003_d = I_ERI_Fx2y_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dxz_Pz_C1001000003_d = I_ERI_Fxyz_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dxz_Pz_C1001000003_d = I_ERI_Fx2z_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dxz_Pz_C1001000003_d = I_ERI_F3y_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F3y_S_Dxz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dxz_Pz_C1001000003_d = I_ERI_F2yz_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F2yz_S_Dxz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dxz_Pz_C1001000003_d = I_ERI_Fy2z_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dxz_Pz_C1001000003_d = I_ERI_F3z_S_Fx2z_S_C1001000003_d+CDZ*I_ERI_F3z_S_Dxz_S_C1001000003_d;
  Double I_ERI_F3x_S_Dyz_Pz_C1001000003_d = I_ERI_F3x_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F3x_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xy_S_Dyz_Pz_C1001000003_d = I_ERI_F2xy_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F2xy_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2xz_S_Dyz_Pz_C1001000003_d = I_ERI_F2xz_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F2xz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2y_S_Dyz_Pz_C1001000003_d = I_ERI_Fx2y_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fxyz_S_Dyz_Pz_C1001000003_d = I_ERI_Fxyz_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fx2z_S_Dyz_Pz_C1001000003_d = I_ERI_Fx2z_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3y_S_Dyz_Pz_C1001000003_d = I_ERI_F3y_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F3y_S_Dyz_S_C1001000003_d;
  Double I_ERI_F2yz_S_Dyz_Pz_C1001000003_d = I_ERI_F2yz_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F2yz_S_Dyz_S_C1001000003_d;
  Double I_ERI_Fy2z_S_Dyz_Pz_C1001000003_d = I_ERI_Fy2z_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3z_S_Dyz_Pz_C1001000003_d = I_ERI_F3z_S_Fy2z_S_C1001000003_d+CDZ*I_ERI_F3z_S_Dyz_S_C1001000003_d;
  Double I_ERI_F3x_S_D2z_Pz_C1001000003_d = I_ERI_F3x_S_F3z_S_C1001000003_d+CDZ*I_ERI_F3x_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xy_S_D2z_Pz_C1001000003_d = I_ERI_F2xy_S_F3z_S_C1001000003_d+CDZ*I_ERI_F2xy_S_D2z_S_C1001000003_d;
  Double I_ERI_F2xz_S_D2z_Pz_C1001000003_d = I_ERI_F2xz_S_F3z_S_C1001000003_d+CDZ*I_ERI_F2xz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2y_S_D2z_Pz_C1001000003_d = I_ERI_Fx2y_S_F3z_S_C1001000003_d+CDZ*I_ERI_Fx2y_S_D2z_S_C1001000003_d;
  Double I_ERI_Fxyz_S_D2z_Pz_C1001000003_d = I_ERI_Fxyz_S_F3z_S_C1001000003_d+CDZ*I_ERI_Fxyz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fx2z_S_D2z_Pz_C1001000003_d = I_ERI_Fx2z_S_F3z_S_C1001000003_d+CDZ*I_ERI_Fx2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3y_S_D2z_Pz_C1001000003_d = I_ERI_F3y_S_F3z_S_C1001000003_d+CDZ*I_ERI_F3y_S_D2z_S_C1001000003_d;
  Double I_ERI_F2yz_S_D2z_Pz_C1001000003_d = I_ERI_F2yz_S_F3z_S_C1001000003_d+CDZ*I_ERI_F2yz_S_D2z_S_C1001000003_d;
  Double I_ERI_Fy2z_S_D2z_Pz_C1001000003_d = I_ERI_Fy2z_S_F3z_S_C1001000003_d+CDZ*I_ERI_Fy2z_S_D2z_S_C1001000003_d;
  Double I_ERI_F3z_S_D2z_Pz_C1001000003_d = I_ERI_F3z_S_F3z_S_C1001000003_d+CDZ*I_ERI_F3z_S_D2z_S_C1001000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_D_C1001000003_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1001000003_d
   ************************************************************/
  Double I_ERI_F3x_S_Px_D2x_C1001000003_d = I_ERI_F3x_S_D2x_Px_C1001000003_d+CDX*I_ERI_F3x_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Px_D2x_C1001000003_d = I_ERI_F2xy_S_D2x_Px_C1001000003_d+CDX*I_ERI_F2xy_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Px_D2x_C1001000003_d = I_ERI_F2xz_S_D2x_Px_C1001000003_d+CDX*I_ERI_F2xz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_D2x_C1001000003_d = I_ERI_Fx2y_S_D2x_Px_C1001000003_d+CDX*I_ERI_Fx2y_S_Px_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_D2x_C1001000003_d = I_ERI_Fxyz_S_D2x_Px_C1001000003_d+CDX*I_ERI_Fxyz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_D2x_C1001000003_d = I_ERI_Fx2z_S_D2x_Px_C1001000003_d+CDX*I_ERI_Fx2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3y_S_Px_D2x_C1001000003_d = I_ERI_F3y_S_D2x_Px_C1001000003_d+CDX*I_ERI_F3y_S_Px_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Px_D2x_C1001000003_d = I_ERI_F2yz_S_D2x_Px_C1001000003_d+CDX*I_ERI_F2yz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_D2x_C1001000003_d = I_ERI_Fy2z_S_D2x_Px_C1001000003_d+CDX*I_ERI_Fy2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3z_S_Px_D2x_C1001000003_d = I_ERI_F3z_S_D2x_Px_C1001000003_d+CDX*I_ERI_F3z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3x_S_Py_D2x_C1001000003_d = I_ERI_F3x_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F3x_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Py_D2x_C1001000003_d = I_ERI_F2xy_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F2xy_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Py_D2x_C1001000003_d = I_ERI_F2xz_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F2xz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_D2x_C1001000003_d = I_ERI_Fx2y_S_Dxy_Px_C1001000003_d+CDX*I_ERI_Fx2y_S_Py_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_D2x_C1001000003_d = I_ERI_Fxyz_S_Dxy_Px_C1001000003_d+CDX*I_ERI_Fxyz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_D2x_C1001000003_d = I_ERI_Fx2z_S_Dxy_Px_C1001000003_d+CDX*I_ERI_Fx2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3y_S_Py_D2x_C1001000003_d = I_ERI_F3y_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F3y_S_Py_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Py_D2x_C1001000003_d = I_ERI_F2yz_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F2yz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_D2x_C1001000003_d = I_ERI_Fy2z_S_Dxy_Px_C1001000003_d+CDX*I_ERI_Fy2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3z_S_Py_D2x_C1001000003_d = I_ERI_F3z_S_Dxy_Px_C1001000003_d+CDX*I_ERI_F3z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3x_S_Pz_D2x_C1001000003_d = I_ERI_F3x_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F3x_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_D2x_C1001000003_d = I_ERI_F2xy_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F2xy_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_D2x_C1001000003_d = I_ERI_F2xz_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F2xz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_D2x_C1001000003_d = I_ERI_Fx2y_S_Dxz_Px_C1001000003_d+CDX*I_ERI_Fx2y_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_D2x_C1001000003_d = I_ERI_Fxyz_S_Dxz_Px_C1001000003_d+CDX*I_ERI_Fxyz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_D2x_C1001000003_d = I_ERI_Fx2z_S_Dxz_Px_C1001000003_d+CDX*I_ERI_Fx2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3y_S_Pz_D2x_C1001000003_d = I_ERI_F3y_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F3y_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_D2x_C1001000003_d = I_ERI_F2yz_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F2yz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_D2x_C1001000003_d = I_ERI_Fy2z_S_Dxz_Px_C1001000003_d+CDX*I_ERI_Fy2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3z_S_Pz_D2x_C1001000003_d = I_ERI_F3z_S_Dxz_Px_C1001000003_d+CDX*I_ERI_F3z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3x_S_Px_Dxy_C1001000003_d = I_ERI_F3x_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F3x_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Dxy_C1001000003_d = I_ERI_F2xy_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F2xy_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Dxy_C1001000003_d = I_ERI_F2xz_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F2xz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Dxy_C1001000003_d = I_ERI_Fx2y_S_Dxy_Px_C1001000003_d+CDY*I_ERI_Fx2y_S_Px_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Dxy_C1001000003_d = I_ERI_Fxyz_S_Dxy_Px_C1001000003_d+CDY*I_ERI_Fxyz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Dxy_C1001000003_d = I_ERI_Fx2z_S_Dxy_Px_C1001000003_d+CDY*I_ERI_Fx2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3y_S_Px_Dxy_C1001000003_d = I_ERI_F3y_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F3y_S_Px_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Dxy_C1001000003_d = I_ERI_F2yz_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F2yz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Dxy_C1001000003_d = I_ERI_Fy2z_S_Dxy_Px_C1001000003_d+CDY*I_ERI_Fy2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3z_S_Px_Dxy_C1001000003_d = I_ERI_F3z_S_Dxy_Px_C1001000003_d+CDY*I_ERI_F3z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3x_S_Py_Dxy_C1001000003_d = I_ERI_F3x_S_D2y_Px_C1001000003_d+CDY*I_ERI_F3x_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Dxy_C1001000003_d = I_ERI_F2xy_S_D2y_Px_C1001000003_d+CDY*I_ERI_F2xy_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Dxy_C1001000003_d = I_ERI_F2xz_S_D2y_Px_C1001000003_d+CDY*I_ERI_F2xz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Dxy_C1001000003_d = I_ERI_Fx2y_S_D2y_Px_C1001000003_d+CDY*I_ERI_Fx2y_S_Py_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Dxy_C1001000003_d = I_ERI_Fxyz_S_D2y_Px_C1001000003_d+CDY*I_ERI_Fxyz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Dxy_C1001000003_d = I_ERI_Fx2z_S_D2y_Px_C1001000003_d+CDY*I_ERI_Fx2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3y_S_Py_Dxy_C1001000003_d = I_ERI_F3y_S_D2y_Px_C1001000003_d+CDY*I_ERI_F3y_S_Py_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Dxy_C1001000003_d = I_ERI_F2yz_S_D2y_Px_C1001000003_d+CDY*I_ERI_F2yz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Dxy_C1001000003_d = I_ERI_Fy2z_S_D2y_Px_C1001000003_d+CDY*I_ERI_Fy2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3z_S_Py_Dxy_C1001000003_d = I_ERI_F3z_S_D2y_Px_C1001000003_d+CDY*I_ERI_F3z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Dxy_C1001000003_d = I_ERI_F3x_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F3x_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Dxy_C1001000003_d = I_ERI_F2xy_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F2xy_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Dxy_C1001000003_d = I_ERI_F2xz_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F2xz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Dxy_C1001000003_d = I_ERI_Fx2y_S_Dyz_Px_C1001000003_d+CDY*I_ERI_Fx2y_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Dxy_C1001000003_d = I_ERI_Fxyz_S_Dyz_Px_C1001000003_d+CDY*I_ERI_Fxyz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Dxy_C1001000003_d = I_ERI_Fx2z_S_Dyz_Px_C1001000003_d+CDY*I_ERI_Fx2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Dxy_C1001000003_d = I_ERI_F3y_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F3y_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Dxy_C1001000003_d = I_ERI_F2yz_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F2yz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Dxy_C1001000003_d = I_ERI_Fy2z_S_Dyz_Px_C1001000003_d+CDY*I_ERI_Fy2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Dxy_C1001000003_d = I_ERI_F3z_S_Dyz_Px_C1001000003_d+CDY*I_ERI_F3z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3x_S_Px_Dxz_C1001000003_d = I_ERI_F3x_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F3x_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Dxz_C1001000003_d = I_ERI_F2xy_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F2xy_S_Px_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Dxz_C1001000003_d = I_ERI_F2xz_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F2xz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Dxz_C1001000003_d = I_ERI_Fx2y_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_Fx2y_S_Px_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Dxz_C1001000003_d = I_ERI_Fxyz_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_Fxyz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Dxz_C1001000003_d = I_ERI_Fx2z_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_Fx2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3y_S_Px_Dxz_C1001000003_d = I_ERI_F3y_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F3y_S_Px_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Dxz_C1001000003_d = I_ERI_F2yz_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F2yz_S_Px_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Dxz_C1001000003_d = I_ERI_Fy2z_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_Fy2z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3z_S_Px_Dxz_C1001000003_d = I_ERI_F3z_S_Dxz_Px_C1001000003_d+CDZ*I_ERI_F3z_S_Px_Px_C1001000003_d;
  Double I_ERI_F3x_S_Py_Dxz_C1001000003_d = I_ERI_F3x_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F3x_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Dxz_C1001000003_d = I_ERI_F2xy_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F2xy_S_Py_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Dxz_C1001000003_d = I_ERI_F2xz_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F2xz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Dxz_C1001000003_d = I_ERI_Fx2y_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_Fx2y_S_Py_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Dxz_C1001000003_d = I_ERI_Fxyz_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_Fxyz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Dxz_C1001000003_d = I_ERI_Fx2z_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_Fx2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3y_S_Py_Dxz_C1001000003_d = I_ERI_F3y_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F3y_S_Py_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Dxz_C1001000003_d = I_ERI_F2yz_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F2yz_S_Py_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Dxz_C1001000003_d = I_ERI_Fy2z_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_Fy2z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3z_S_Py_Dxz_C1001000003_d = I_ERI_F3z_S_Dyz_Px_C1001000003_d+CDZ*I_ERI_F3z_S_Py_Px_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Dxz_C1001000003_d = I_ERI_F3x_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F3x_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Dxz_C1001000003_d = I_ERI_F2xy_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F2xy_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Dxz_C1001000003_d = I_ERI_F2xz_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F2xz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Dxz_C1001000003_d = I_ERI_Fx2y_S_D2z_Px_C1001000003_d+CDZ*I_ERI_Fx2y_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Dxz_C1001000003_d = I_ERI_Fxyz_S_D2z_Px_C1001000003_d+CDZ*I_ERI_Fxyz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Dxz_C1001000003_d = I_ERI_Fx2z_S_D2z_Px_C1001000003_d+CDZ*I_ERI_Fx2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Dxz_C1001000003_d = I_ERI_F3y_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F3y_S_Pz_Px_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Dxz_C1001000003_d = I_ERI_F2yz_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F2yz_S_Pz_Px_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Dxz_C1001000003_d = I_ERI_Fy2z_S_D2z_Px_C1001000003_d+CDZ*I_ERI_Fy2z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Dxz_C1001000003_d = I_ERI_F3z_S_D2z_Px_C1001000003_d+CDZ*I_ERI_F3z_S_Pz_Px_C1001000003_d;
  Double I_ERI_F3x_S_Px_D2y_C1001000003_d = I_ERI_F3x_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F3x_S_Px_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Px_D2y_C1001000003_d = I_ERI_F2xy_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F2xy_S_Px_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Px_D2y_C1001000003_d = I_ERI_F2xz_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F2xz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_D2y_C1001000003_d = I_ERI_Fx2y_S_Dxy_Py_C1001000003_d+CDY*I_ERI_Fx2y_S_Px_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_D2y_C1001000003_d = I_ERI_Fxyz_S_Dxy_Py_C1001000003_d+CDY*I_ERI_Fxyz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_D2y_C1001000003_d = I_ERI_Fx2z_S_Dxy_Py_C1001000003_d+CDY*I_ERI_Fx2z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3y_S_Px_D2y_C1001000003_d = I_ERI_F3y_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F3y_S_Px_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Px_D2y_C1001000003_d = I_ERI_F2yz_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F2yz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_D2y_C1001000003_d = I_ERI_Fy2z_S_Dxy_Py_C1001000003_d+CDY*I_ERI_Fy2z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3z_S_Px_D2y_C1001000003_d = I_ERI_F3z_S_Dxy_Py_C1001000003_d+CDY*I_ERI_F3z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3x_S_Py_D2y_C1001000003_d = I_ERI_F3x_S_D2y_Py_C1001000003_d+CDY*I_ERI_F3x_S_Py_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Py_D2y_C1001000003_d = I_ERI_F2xy_S_D2y_Py_C1001000003_d+CDY*I_ERI_F2xy_S_Py_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Py_D2y_C1001000003_d = I_ERI_F2xz_S_D2y_Py_C1001000003_d+CDY*I_ERI_F2xz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_D2y_C1001000003_d = I_ERI_Fx2y_S_D2y_Py_C1001000003_d+CDY*I_ERI_Fx2y_S_Py_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_D2y_C1001000003_d = I_ERI_Fxyz_S_D2y_Py_C1001000003_d+CDY*I_ERI_Fxyz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_D2y_C1001000003_d = I_ERI_Fx2z_S_D2y_Py_C1001000003_d+CDY*I_ERI_Fx2z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3y_S_Py_D2y_C1001000003_d = I_ERI_F3y_S_D2y_Py_C1001000003_d+CDY*I_ERI_F3y_S_Py_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Py_D2y_C1001000003_d = I_ERI_F2yz_S_D2y_Py_C1001000003_d+CDY*I_ERI_F2yz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_D2y_C1001000003_d = I_ERI_Fy2z_S_D2y_Py_C1001000003_d+CDY*I_ERI_Fy2z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3z_S_Py_D2y_C1001000003_d = I_ERI_F3z_S_D2y_Py_C1001000003_d+CDY*I_ERI_F3z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3x_S_Pz_D2y_C1001000003_d = I_ERI_F3x_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F3x_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_D2y_C1001000003_d = I_ERI_F2xy_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F2xy_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_D2y_C1001000003_d = I_ERI_F2xz_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F2xz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_D2y_C1001000003_d = I_ERI_Fx2y_S_Dyz_Py_C1001000003_d+CDY*I_ERI_Fx2y_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_D2y_C1001000003_d = I_ERI_Fxyz_S_Dyz_Py_C1001000003_d+CDY*I_ERI_Fxyz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_D2y_C1001000003_d = I_ERI_Fx2z_S_Dyz_Py_C1001000003_d+CDY*I_ERI_Fx2z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3y_S_Pz_D2y_C1001000003_d = I_ERI_F3y_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F3y_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_D2y_C1001000003_d = I_ERI_F2yz_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F2yz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_D2y_C1001000003_d = I_ERI_Fy2z_S_Dyz_Py_C1001000003_d+CDY*I_ERI_Fy2z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3z_S_Pz_D2y_C1001000003_d = I_ERI_F3z_S_Dyz_Py_C1001000003_d+CDY*I_ERI_F3z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3x_S_Px_Dyz_C1001000003_d = I_ERI_F3x_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F3x_S_Px_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Px_Dyz_C1001000003_d = I_ERI_F2xy_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F2xy_S_Px_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Px_Dyz_C1001000003_d = I_ERI_F2xz_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F2xz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_Dyz_C1001000003_d = I_ERI_Fx2y_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_Fx2y_S_Px_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_Dyz_C1001000003_d = I_ERI_Fxyz_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_Fxyz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_Dyz_C1001000003_d = I_ERI_Fx2z_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_Fx2z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3y_S_Px_Dyz_C1001000003_d = I_ERI_F3y_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F3y_S_Px_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Px_Dyz_C1001000003_d = I_ERI_F2yz_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F2yz_S_Px_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_Dyz_C1001000003_d = I_ERI_Fy2z_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_Fy2z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3z_S_Px_Dyz_C1001000003_d = I_ERI_F3z_S_Dxz_Py_C1001000003_d+CDZ*I_ERI_F3z_S_Px_Py_C1001000003_d;
  Double I_ERI_F3x_S_Py_Dyz_C1001000003_d = I_ERI_F3x_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F3x_S_Py_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Py_Dyz_C1001000003_d = I_ERI_F2xy_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F2xy_S_Py_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Py_Dyz_C1001000003_d = I_ERI_F2xz_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F2xz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_Dyz_C1001000003_d = I_ERI_Fx2y_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_Fx2y_S_Py_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_Dyz_C1001000003_d = I_ERI_Fxyz_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_Fxyz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_Dyz_C1001000003_d = I_ERI_Fx2z_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_Fx2z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3y_S_Py_Dyz_C1001000003_d = I_ERI_F3y_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F3y_S_Py_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Py_Dyz_C1001000003_d = I_ERI_F2yz_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F2yz_S_Py_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_Dyz_C1001000003_d = I_ERI_Fy2z_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_Fy2z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3z_S_Py_Dyz_C1001000003_d = I_ERI_F3z_S_Dyz_Py_C1001000003_d+CDZ*I_ERI_F3z_S_Py_Py_C1001000003_d;
  Double I_ERI_F3x_S_Pz_Dyz_C1001000003_d = I_ERI_F3x_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F3x_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_Dyz_C1001000003_d = I_ERI_F2xy_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F2xy_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_Dyz_C1001000003_d = I_ERI_F2xz_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F2xz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_Dyz_C1001000003_d = I_ERI_Fx2y_S_D2z_Py_C1001000003_d+CDZ*I_ERI_Fx2y_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_Dyz_C1001000003_d = I_ERI_Fxyz_S_D2z_Py_C1001000003_d+CDZ*I_ERI_Fxyz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_Dyz_C1001000003_d = I_ERI_Fx2z_S_D2z_Py_C1001000003_d+CDZ*I_ERI_Fx2z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3y_S_Pz_Dyz_C1001000003_d = I_ERI_F3y_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F3y_S_Pz_Py_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_Dyz_C1001000003_d = I_ERI_F2yz_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F2yz_S_Pz_Py_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_Dyz_C1001000003_d = I_ERI_Fy2z_S_D2z_Py_C1001000003_d+CDZ*I_ERI_Fy2z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3z_S_Pz_Dyz_C1001000003_d = I_ERI_F3z_S_D2z_Py_C1001000003_d+CDZ*I_ERI_F3z_S_Pz_Py_C1001000003_d;
  Double I_ERI_F3x_S_Px_D2z_C1001000003_d = I_ERI_F3x_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F3x_S_Px_Pz_C1001000003_d;
  Double I_ERI_F2xy_S_Px_D2z_C1001000003_d = I_ERI_F2xy_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F2xy_S_Px_Pz_C1001000003_d;
  Double I_ERI_F2xz_S_Px_D2z_C1001000003_d = I_ERI_F2xz_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F2xz_S_Px_Pz_C1001000003_d;
  Double I_ERI_Fx2y_S_Px_D2z_C1001000003_d = I_ERI_Fx2y_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_Fx2y_S_Px_Pz_C1001000003_d;
  Double I_ERI_Fxyz_S_Px_D2z_C1001000003_d = I_ERI_Fxyz_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_Fxyz_S_Px_Pz_C1001000003_d;
  Double I_ERI_Fx2z_S_Px_D2z_C1001000003_d = I_ERI_Fx2z_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_Fx2z_S_Px_Pz_C1001000003_d;
  Double I_ERI_F3y_S_Px_D2z_C1001000003_d = I_ERI_F3y_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F3y_S_Px_Pz_C1001000003_d;
  Double I_ERI_F2yz_S_Px_D2z_C1001000003_d = I_ERI_F2yz_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F2yz_S_Px_Pz_C1001000003_d;
  Double I_ERI_Fy2z_S_Px_D2z_C1001000003_d = I_ERI_Fy2z_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_Fy2z_S_Px_Pz_C1001000003_d;
  Double I_ERI_F3z_S_Px_D2z_C1001000003_d = I_ERI_F3z_S_Dxz_Pz_C1001000003_d+CDZ*I_ERI_F3z_S_Px_Pz_C1001000003_d;
  Double I_ERI_F3x_S_Py_D2z_C1001000003_d = I_ERI_F3x_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F3x_S_Py_Pz_C1001000003_d;
  Double I_ERI_F2xy_S_Py_D2z_C1001000003_d = I_ERI_F2xy_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F2xy_S_Py_Pz_C1001000003_d;
  Double I_ERI_F2xz_S_Py_D2z_C1001000003_d = I_ERI_F2xz_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F2xz_S_Py_Pz_C1001000003_d;
  Double I_ERI_Fx2y_S_Py_D2z_C1001000003_d = I_ERI_Fx2y_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_Fx2y_S_Py_Pz_C1001000003_d;
  Double I_ERI_Fxyz_S_Py_D2z_C1001000003_d = I_ERI_Fxyz_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_Fxyz_S_Py_Pz_C1001000003_d;
  Double I_ERI_Fx2z_S_Py_D2z_C1001000003_d = I_ERI_Fx2z_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_Fx2z_S_Py_Pz_C1001000003_d;
  Double I_ERI_F3y_S_Py_D2z_C1001000003_d = I_ERI_F3y_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F3y_S_Py_Pz_C1001000003_d;
  Double I_ERI_F2yz_S_Py_D2z_C1001000003_d = I_ERI_F2yz_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F2yz_S_Py_Pz_C1001000003_d;
  Double I_ERI_Fy2z_S_Py_D2z_C1001000003_d = I_ERI_Fy2z_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_Fy2z_S_Py_Pz_C1001000003_d;
  Double I_ERI_F3z_S_Py_D2z_C1001000003_d = I_ERI_F3z_S_Dyz_Pz_C1001000003_d+CDZ*I_ERI_F3z_S_Py_Pz_C1001000003_d;
  Double I_ERI_F3x_S_Pz_D2z_C1001000003_d = I_ERI_F3x_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F3x_S_Pz_Pz_C1001000003_d;
  Double I_ERI_F2xy_S_Pz_D2z_C1001000003_d = I_ERI_F2xy_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F2xy_S_Pz_Pz_C1001000003_d;
  Double I_ERI_F2xz_S_Pz_D2z_C1001000003_d = I_ERI_F2xz_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F2xz_S_Pz_Pz_C1001000003_d;
  Double I_ERI_Fx2y_S_Pz_D2z_C1001000003_d = I_ERI_Fx2y_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_Fx2y_S_Pz_Pz_C1001000003_d;
  Double I_ERI_Fxyz_S_Pz_D2z_C1001000003_d = I_ERI_Fxyz_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_Fxyz_S_Pz_Pz_C1001000003_d;
  Double I_ERI_Fx2z_S_Pz_D2z_C1001000003_d = I_ERI_Fx2z_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_Fx2z_S_Pz_Pz_C1001000003_d;
  Double I_ERI_F3y_S_Pz_D2z_C1001000003_d = I_ERI_F3y_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F3y_S_Pz_Pz_C1001000003_d;
  Double I_ERI_F2yz_S_Pz_D2z_C1001000003_d = I_ERI_F2yz_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F2yz_S_Pz_Pz_C1001000003_d;
  Double I_ERI_Fy2z_S_Pz_D2z_C1001000003_d = I_ERI_Fy2z_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_Fy2z_S_Pz_Pz_C1001000003_d;
  Double I_ERI_F3z_S_Pz_D2z_C1001000003_d = I_ERI_F3z_S_D2z_Pz_C1001000003_d+CDZ*I_ERI_F3z_S_Pz_Pz_C1001000003_d;

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
   * shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1000003_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000003_b
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1000003_b = I_ERI_G4x_S_Px_S_C1000003_b+ABX*I_ERI_F3x_S_Px_S_C1000003_b;
  Double I_ERI_F2xy_Px_Px_S_C1000003_b = I_ERI_G3xy_S_Px_S_C1000003_b+ABX*I_ERI_F2xy_S_Px_S_C1000003_b;
  Double I_ERI_F2xz_Px_Px_S_C1000003_b = I_ERI_G3xz_S_Px_S_C1000003_b+ABX*I_ERI_F2xz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2y_Px_Px_S_C1000003_b = I_ERI_G2x2y_S_Px_S_C1000003_b+ABX*I_ERI_Fx2y_S_Px_S_C1000003_b;
  Double I_ERI_Fxyz_Px_Px_S_C1000003_b = I_ERI_G2xyz_S_Px_S_C1000003_b+ABX*I_ERI_Fxyz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2z_Px_Px_S_C1000003_b = I_ERI_G2x2z_S_Px_S_C1000003_b+ABX*I_ERI_Fx2z_S_Px_S_C1000003_b;
  Double I_ERI_F3y_Px_Px_S_C1000003_b = I_ERI_Gx3y_S_Px_S_C1000003_b+ABX*I_ERI_F3y_S_Px_S_C1000003_b;
  Double I_ERI_F2yz_Px_Px_S_C1000003_b = I_ERI_Gx2yz_S_Px_S_C1000003_b+ABX*I_ERI_F2yz_S_Px_S_C1000003_b;
  Double I_ERI_Fy2z_Px_Px_S_C1000003_b = I_ERI_Gxy2z_S_Px_S_C1000003_b+ABX*I_ERI_Fy2z_S_Px_S_C1000003_b;
  Double I_ERI_F3z_Px_Px_S_C1000003_b = I_ERI_Gx3z_S_Px_S_C1000003_b+ABX*I_ERI_F3z_S_Px_S_C1000003_b;
  Double I_ERI_F3x_Py_Px_S_C1000003_b = I_ERI_G3xy_S_Px_S_C1000003_b+ABY*I_ERI_F3x_S_Px_S_C1000003_b;
  Double I_ERI_F2xy_Py_Px_S_C1000003_b = I_ERI_G2x2y_S_Px_S_C1000003_b+ABY*I_ERI_F2xy_S_Px_S_C1000003_b;
  Double I_ERI_F2xz_Py_Px_S_C1000003_b = I_ERI_G2xyz_S_Px_S_C1000003_b+ABY*I_ERI_F2xz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2y_Py_Px_S_C1000003_b = I_ERI_Gx3y_S_Px_S_C1000003_b+ABY*I_ERI_Fx2y_S_Px_S_C1000003_b;
  Double I_ERI_Fxyz_Py_Px_S_C1000003_b = I_ERI_Gx2yz_S_Px_S_C1000003_b+ABY*I_ERI_Fxyz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2z_Py_Px_S_C1000003_b = I_ERI_Gxy2z_S_Px_S_C1000003_b+ABY*I_ERI_Fx2z_S_Px_S_C1000003_b;
  Double I_ERI_F3y_Py_Px_S_C1000003_b = I_ERI_G4y_S_Px_S_C1000003_b+ABY*I_ERI_F3y_S_Px_S_C1000003_b;
  Double I_ERI_F2yz_Py_Px_S_C1000003_b = I_ERI_G3yz_S_Px_S_C1000003_b+ABY*I_ERI_F2yz_S_Px_S_C1000003_b;
  Double I_ERI_Fy2z_Py_Px_S_C1000003_b = I_ERI_G2y2z_S_Px_S_C1000003_b+ABY*I_ERI_Fy2z_S_Px_S_C1000003_b;
  Double I_ERI_F3z_Py_Px_S_C1000003_b = I_ERI_Gy3z_S_Px_S_C1000003_b+ABY*I_ERI_F3z_S_Px_S_C1000003_b;
  Double I_ERI_F3x_Pz_Px_S_C1000003_b = I_ERI_G3xz_S_Px_S_C1000003_b+ABZ*I_ERI_F3x_S_Px_S_C1000003_b;
  Double I_ERI_F2xy_Pz_Px_S_C1000003_b = I_ERI_G2xyz_S_Px_S_C1000003_b+ABZ*I_ERI_F2xy_S_Px_S_C1000003_b;
  Double I_ERI_F2xz_Pz_Px_S_C1000003_b = I_ERI_G2x2z_S_Px_S_C1000003_b+ABZ*I_ERI_F2xz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2y_Pz_Px_S_C1000003_b = I_ERI_Gx2yz_S_Px_S_C1000003_b+ABZ*I_ERI_Fx2y_S_Px_S_C1000003_b;
  Double I_ERI_Fxyz_Pz_Px_S_C1000003_b = I_ERI_Gxy2z_S_Px_S_C1000003_b+ABZ*I_ERI_Fxyz_S_Px_S_C1000003_b;
  Double I_ERI_Fx2z_Pz_Px_S_C1000003_b = I_ERI_Gx3z_S_Px_S_C1000003_b+ABZ*I_ERI_Fx2z_S_Px_S_C1000003_b;
  Double I_ERI_F3y_Pz_Px_S_C1000003_b = I_ERI_G3yz_S_Px_S_C1000003_b+ABZ*I_ERI_F3y_S_Px_S_C1000003_b;
  Double I_ERI_F2yz_Pz_Px_S_C1000003_b = I_ERI_G2y2z_S_Px_S_C1000003_b+ABZ*I_ERI_F2yz_S_Px_S_C1000003_b;
  Double I_ERI_Fy2z_Pz_Px_S_C1000003_b = I_ERI_Gy3z_S_Px_S_C1000003_b+ABZ*I_ERI_Fy2z_S_Px_S_C1000003_b;
  Double I_ERI_F3z_Pz_Px_S_C1000003_b = I_ERI_G4z_S_Px_S_C1000003_b+ABZ*I_ERI_F3z_S_Px_S_C1000003_b;
  Double I_ERI_F3x_Px_Py_S_C1000003_b = I_ERI_G4x_S_Py_S_C1000003_b+ABX*I_ERI_F3x_S_Py_S_C1000003_b;
  Double I_ERI_F2xy_Px_Py_S_C1000003_b = I_ERI_G3xy_S_Py_S_C1000003_b+ABX*I_ERI_F2xy_S_Py_S_C1000003_b;
  Double I_ERI_F2xz_Px_Py_S_C1000003_b = I_ERI_G3xz_S_Py_S_C1000003_b+ABX*I_ERI_F2xz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2y_Px_Py_S_C1000003_b = I_ERI_G2x2y_S_Py_S_C1000003_b+ABX*I_ERI_Fx2y_S_Py_S_C1000003_b;
  Double I_ERI_Fxyz_Px_Py_S_C1000003_b = I_ERI_G2xyz_S_Py_S_C1000003_b+ABX*I_ERI_Fxyz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2z_Px_Py_S_C1000003_b = I_ERI_G2x2z_S_Py_S_C1000003_b+ABX*I_ERI_Fx2z_S_Py_S_C1000003_b;
  Double I_ERI_F3y_Px_Py_S_C1000003_b = I_ERI_Gx3y_S_Py_S_C1000003_b+ABX*I_ERI_F3y_S_Py_S_C1000003_b;
  Double I_ERI_F2yz_Px_Py_S_C1000003_b = I_ERI_Gx2yz_S_Py_S_C1000003_b+ABX*I_ERI_F2yz_S_Py_S_C1000003_b;
  Double I_ERI_Fy2z_Px_Py_S_C1000003_b = I_ERI_Gxy2z_S_Py_S_C1000003_b+ABX*I_ERI_Fy2z_S_Py_S_C1000003_b;
  Double I_ERI_F3z_Px_Py_S_C1000003_b = I_ERI_Gx3z_S_Py_S_C1000003_b+ABX*I_ERI_F3z_S_Py_S_C1000003_b;
  Double I_ERI_F3x_Py_Py_S_C1000003_b = I_ERI_G3xy_S_Py_S_C1000003_b+ABY*I_ERI_F3x_S_Py_S_C1000003_b;
  Double I_ERI_F2xy_Py_Py_S_C1000003_b = I_ERI_G2x2y_S_Py_S_C1000003_b+ABY*I_ERI_F2xy_S_Py_S_C1000003_b;
  Double I_ERI_F2xz_Py_Py_S_C1000003_b = I_ERI_G2xyz_S_Py_S_C1000003_b+ABY*I_ERI_F2xz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2y_Py_Py_S_C1000003_b = I_ERI_Gx3y_S_Py_S_C1000003_b+ABY*I_ERI_Fx2y_S_Py_S_C1000003_b;
  Double I_ERI_Fxyz_Py_Py_S_C1000003_b = I_ERI_Gx2yz_S_Py_S_C1000003_b+ABY*I_ERI_Fxyz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2z_Py_Py_S_C1000003_b = I_ERI_Gxy2z_S_Py_S_C1000003_b+ABY*I_ERI_Fx2z_S_Py_S_C1000003_b;
  Double I_ERI_F3y_Py_Py_S_C1000003_b = I_ERI_G4y_S_Py_S_C1000003_b+ABY*I_ERI_F3y_S_Py_S_C1000003_b;
  Double I_ERI_F2yz_Py_Py_S_C1000003_b = I_ERI_G3yz_S_Py_S_C1000003_b+ABY*I_ERI_F2yz_S_Py_S_C1000003_b;
  Double I_ERI_Fy2z_Py_Py_S_C1000003_b = I_ERI_G2y2z_S_Py_S_C1000003_b+ABY*I_ERI_Fy2z_S_Py_S_C1000003_b;
  Double I_ERI_F3z_Py_Py_S_C1000003_b = I_ERI_Gy3z_S_Py_S_C1000003_b+ABY*I_ERI_F3z_S_Py_S_C1000003_b;
  Double I_ERI_F3x_Pz_Py_S_C1000003_b = I_ERI_G3xz_S_Py_S_C1000003_b+ABZ*I_ERI_F3x_S_Py_S_C1000003_b;
  Double I_ERI_F2xy_Pz_Py_S_C1000003_b = I_ERI_G2xyz_S_Py_S_C1000003_b+ABZ*I_ERI_F2xy_S_Py_S_C1000003_b;
  Double I_ERI_F2xz_Pz_Py_S_C1000003_b = I_ERI_G2x2z_S_Py_S_C1000003_b+ABZ*I_ERI_F2xz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2y_Pz_Py_S_C1000003_b = I_ERI_Gx2yz_S_Py_S_C1000003_b+ABZ*I_ERI_Fx2y_S_Py_S_C1000003_b;
  Double I_ERI_Fxyz_Pz_Py_S_C1000003_b = I_ERI_Gxy2z_S_Py_S_C1000003_b+ABZ*I_ERI_Fxyz_S_Py_S_C1000003_b;
  Double I_ERI_Fx2z_Pz_Py_S_C1000003_b = I_ERI_Gx3z_S_Py_S_C1000003_b+ABZ*I_ERI_Fx2z_S_Py_S_C1000003_b;
  Double I_ERI_F3y_Pz_Py_S_C1000003_b = I_ERI_G3yz_S_Py_S_C1000003_b+ABZ*I_ERI_F3y_S_Py_S_C1000003_b;
  Double I_ERI_F2yz_Pz_Py_S_C1000003_b = I_ERI_G2y2z_S_Py_S_C1000003_b+ABZ*I_ERI_F2yz_S_Py_S_C1000003_b;
  Double I_ERI_Fy2z_Pz_Py_S_C1000003_b = I_ERI_Gy3z_S_Py_S_C1000003_b+ABZ*I_ERI_Fy2z_S_Py_S_C1000003_b;
  Double I_ERI_F3z_Pz_Py_S_C1000003_b = I_ERI_G4z_S_Py_S_C1000003_b+ABZ*I_ERI_F3z_S_Py_S_C1000003_b;
  Double I_ERI_F3x_Px_Pz_S_C1000003_b = I_ERI_G4x_S_Pz_S_C1000003_b+ABX*I_ERI_F3x_S_Pz_S_C1000003_b;
  Double I_ERI_F2xy_Px_Pz_S_C1000003_b = I_ERI_G3xy_S_Pz_S_C1000003_b+ABX*I_ERI_F2xy_S_Pz_S_C1000003_b;
  Double I_ERI_F2xz_Px_Pz_S_C1000003_b = I_ERI_G3xz_S_Pz_S_C1000003_b+ABX*I_ERI_F2xz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2y_Px_Pz_S_C1000003_b = I_ERI_G2x2y_S_Pz_S_C1000003_b+ABX*I_ERI_Fx2y_S_Pz_S_C1000003_b;
  Double I_ERI_Fxyz_Px_Pz_S_C1000003_b = I_ERI_G2xyz_S_Pz_S_C1000003_b+ABX*I_ERI_Fxyz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2z_Px_Pz_S_C1000003_b = I_ERI_G2x2z_S_Pz_S_C1000003_b+ABX*I_ERI_Fx2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3y_Px_Pz_S_C1000003_b = I_ERI_Gx3y_S_Pz_S_C1000003_b+ABX*I_ERI_F3y_S_Pz_S_C1000003_b;
  Double I_ERI_F2yz_Px_Pz_S_C1000003_b = I_ERI_Gx2yz_S_Pz_S_C1000003_b+ABX*I_ERI_F2yz_S_Pz_S_C1000003_b;
  Double I_ERI_Fy2z_Px_Pz_S_C1000003_b = I_ERI_Gxy2z_S_Pz_S_C1000003_b+ABX*I_ERI_Fy2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3z_Px_Pz_S_C1000003_b = I_ERI_Gx3z_S_Pz_S_C1000003_b+ABX*I_ERI_F3z_S_Pz_S_C1000003_b;
  Double I_ERI_F3x_Py_Pz_S_C1000003_b = I_ERI_G3xy_S_Pz_S_C1000003_b+ABY*I_ERI_F3x_S_Pz_S_C1000003_b;
  Double I_ERI_F2xy_Py_Pz_S_C1000003_b = I_ERI_G2x2y_S_Pz_S_C1000003_b+ABY*I_ERI_F2xy_S_Pz_S_C1000003_b;
  Double I_ERI_F2xz_Py_Pz_S_C1000003_b = I_ERI_G2xyz_S_Pz_S_C1000003_b+ABY*I_ERI_F2xz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2y_Py_Pz_S_C1000003_b = I_ERI_Gx3y_S_Pz_S_C1000003_b+ABY*I_ERI_Fx2y_S_Pz_S_C1000003_b;
  Double I_ERI_Fxyz_Py_Pz_S_C1000003_b = I_ERI_Gx2yz_S_Pz_S_C1000003_b+ABY*I_ERI_Fxyz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2z_Py_Pz_S_C1000003_b = I_ERI_Gxy2z_S_Pz_S_C1000003_b+ABY*I_ERI_Fx2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3y_Py_Pz_S_C1000003_b = I_ERI_G4y_S_Pz_S_C1000003_b+ABY*I_ERI_F3y_S_Pz_S_C1000003_b;
  Double I_ERI_F2yz_Py_Pz_S_C1000003_b = I_ERI_G3yz_S_Pz_S_C1000003_b+ABY*I_ERI_F2yz_S_Pz_S_C1000003_b;
  Double I_ERI_Fy2z_Py_Pz_S_C1000003_b = I_ERI_G2y2z_S_Pz_S_C1000003_b+ABY*I_ERI_Fy2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3z_Py_Pz_S_C1000003_b = I_ERI_Gy3z_S_Pz_S_C1000003_b+ABY*I_ERI_F3z_S_Pz_S_C1000003_b;
  Double I_ERI_F3x_Pz_Pz_S_C1000003_b = I_ERI_G3xz_S_Pz_S_C1000003_b+ABZ*I_ERI_F3x_S_Pz_S_C1000003_b;
  Double I_ERI_F2xy_Pz_Pz_S_C1000003_b = I_ERI_G2xyz_S_Pz_S_C1000003_b+ABZ*I_ERI_F2xy_S_Pz_S_C1000003_b;
  Double I_ERI_F2xz_Pz_Pz_S_C1000003_b = I_ERI_G2x2z_S_Pz_S_C1000003_b+ABZ*I_ERI_F2xz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2y_Pz_Pz_S_C1000003_b = I_ERI_Gx2yz_S_Pz_S_C1000003_b+ABZ*I_ERI_Fx2y_S_Pz_S_C1000003_b;
  Double I_ERI_Fxyz_Pz_Pz_S_C1000003_b = I_ERI_Gxy2z_S_Pz_S_C1000003_b+ABZ*I_ERI_Fxyz_S_Pz_S_C1000003_b;
  Double I_ERI_Fx2z_Pz_Pz_S_C1000003_b = I_ERI_Gx3z_S_Pz_S_C1000003_b+ABZ*I_ERI_Fx2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3y_Pz_Pz_S_C1000003_b = I_ERI_G3yz_S_Pz_S_C1000003_b+ABZ*I_ERI_F3y_S_Pz_S_C1000003_b;
  Double I_ERI_F2yz_Pz_Pz_S_C1000003_b = I_ERI_G2y2z_S_Pz_S_C1000003_b+ABZ*I_ERI_F2yz_S_Pz_S_C1000003_b;
  Double I_ERI_Fy2z_Pz_Pz_S_C1000003_b = I_ERI_Gy3z_S_Pz_S_C1000003_b+ABZ*I_ERI_Fy2z_S_Pz_S_C1000003_b;
  Double I_ERI_F3z_Pz_Pz_S_C1000003_b = I_ERI_G4z_S_Pz_S_C1000003_b+ABZ*I_ERI_F3z_S_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P_C1001000003_b
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1001000003_b
   ************************************************************/
  Double I_ERI_F3x_Px_Px_Px_C1001000003_b = I_ERI_G4x_S_Px_Px_C1001000003_b+ABX*I_ERI_F3x_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xy_Px_Px_Px_C1001000003_b = I_ERI_G3xy_S_Px_Px_C1001000003_b+ABX*I_ERI_F2xy_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xz_Px_Px_Px_C1001000003_b = I_ERI_G3xz_S_Px_Px_C1001000003_b+ABX*I_ERI_F2xz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2y_Px_Px_Px_C1001000003_b = I_ERI_G2x2y_S_Px_Px_C1001000003_b+ABX*I_ERI_Fx2y_S_Px_Px_C1001000003_b;
  Double I_ERI_Fxyz_Px_Px_Px_C1001000003_b = I_ERI_G2xyz_S_Px_Px_C1001000003_b+ABX*I_ERI_Fxyz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2z_Px_Px_Px_C1001000003_b = I_ERI_G2x2z_S_Px_Px_C1001000003_b+ABX*I_ERI_Fx2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3y_Px_Px_Px_C1001000003_b = I_ERI_Gx3y_S_Px_Px_C1001000003_b+ABX*I_ERI_F3y_S_Px_Px_C1001000003_b;
  Double I_ERI_F2yz_Px_Px_Px_C1001000003_b = I_ERI_Gx2yz_S_Px_Px_C1001000003_b+ABX*I_ERI_F2yz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fy2z_Px_Px_Px_C1001000003_b = I_ERI_Gxy2z_S_Px_Px_C1001000003_b+ABX*I_ERI_Fy2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3z_Px_Px_Px_C1001000003_b = I_ERI_Gx3z_S_Px_Px_C1001000003_b+ABX*I_ERI_F3z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3x_Py_Px_Px_C1001000003_b = I_ERI_G3xy_S_Px_Px_C1001000003_b+ABY*I_ERI_F3x_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xy_Py_Px_Px_C1001000003_b = I_ERI_G2x2y_S_Px_Px_C1001000003_b+ABY*I_ERI_F2xy_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xz_Py_Px_Px_C1001000003_b = I_ERI_G2xyz_S_Px_Px_C1001000003_b+ABY*I_ERI_F2xz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2y_Py_Px_Px_C1001000003_b = I_ERI_Gx3y_S_Px_Px_C1001000003_b+ABY*I_ERI_Fx2y_S_Px_Px_C1001000003_b;
  Double I_ERI_Fxyz_Py_Px_Px_C1001000003_b = I_ERI_Gx2yz_S_Px_Px_C1001000003_b+ABY*I_ERI_Fxyz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2z_Py_Px_Px_C1001000003_b = I_ERI_Gxy2z_S_Px_Px_C1001000003_b+ABY*I_ERI_Fx2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3y_Py_Px_Px_C1001000003_b = I_ERI_G4y_S_Px_Px_C1001000003_b+ABY*I_ERI_F3y_S_Px_Px_C1001000003_b;
  Double I_ERI_F2yz_Py_Px_Px_C1001000003_b = I_ERI_G3yz_S_Px_Px_C1001000003_b+ABY*I_ERI_F2yz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fy2z_Py_Px_Px_C1001000003_b = I_ERI_G2y2z_S_Px_Px_C1001000003_b+ABY*I_ERI_Fy2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3z_Py_Px_Px_C1001000003_b = I_ERI_Gy3z_S_Px_Px_C1001000003_b+ABY*I_ERI_F3z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3x_Pz_Px_Px_C1001000003_b = I_ERI_G3xz_S_Px_Px_C1001000003_b+ABZ*I_ERI_F3x_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xy_Pz_Px_Px_C1001000003_b = I_ERI_G2xyz_S_Px_Px_C1001000003_b+ABZ*I_ERI_F2xy_S_Px_Px_C1001000003_b;
  Double I_ERI_F2xz_Pz_Px_Px_C1001000003_b = I_ERI_G2x2z_S_Px_Px_C1001000003_b+ABZ*I_ERI_F2xz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Px_Px_C1001000003_b = I_ERI_Gx2yz_S_Px_Px_C1001000003_b+ABZ*I_ERI_Fx2y_S_Px_Px_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Px_Px_C1001000003_b = I_ERI_Gxy2z_S_Px_Px_C1001000003_b+ABZ*I_ERI_Fxyz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Px_Px_C1001000003_b = I_ERI_Gx3z_S_Px_Px_C1001000003_b+ABZ*I_ERI_Fx2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3y_Pz_Px_Px_C1001000003_b = I_ERI_G3yz_S_Px_Px_C1001000003_b+ABZ*I_ERI_F3y_S_Px_Px_C1001000003_b;
  Double I_ERI_F2yz_Pz_Px_Px_C1001000003_b = I_ERI_G2y2z_S_Px_Px_C1001000003_b+ABZ*I_ERI_F2yz_S_Px_Px_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Px_Px_C1001000003_b = I_ERI_Gy3z_S_Px_Px_C1001000003_b+ABZ*I_ERI_Fy2z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3z_Pz_Px_Px_C1001000003_b = I_ERI_G4z_S_Px_Px_C1001000003_b+ABZ*I_ERI_F3z_S_Px_Px_C1001000003_b;
  Double I_ERI_F3x_Px_Py_Px_C1001000003_b = I_ERI_G4x_S_Py_Px_C1001000003_b+ABX*I_ERI_F3x_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xy_Px_Py_Px_C1001000003_b = I_ERI_G3xy_S_Py_Px_C1001000003_b+ABX*I_ERI_F2xy_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xz_Px_Py_Px_C1001000003_b = I_ERI_G3xz_S_Py_Px_C1001000003_b+ABX*I_ERI_F2xz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2y_Px_Py_Px_C1001000003_b = I_ERI_G2x2y_S_Py_Px_C1001000003_b+ABX*I_ERI_Fx2y_S_Py_Px_C1001000003_b;
  Double I_ERI_Fxyz_Px_Py_Px_C1001000003_b = I_ERI_G2xyz_S_Py_Px_C1001000003_b+ABX*I_ERI_Fxyz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2z_Px_Py_Px_C1001000003_b = I_ERI_G2x2z_S_Py_Px_C1001000003_b+ABX*I_ERI_Fx2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3y_Px_Py_Px_C1001000003_b = I_ERI_Gx3y_S_Py_Px_C1001000003_b+ABX*I_ERI_F3y_S_Py_Px_C1001000003_b;
  Double I_ERI_F2yz_Px_Py_Px_C1001000003_b = I_ERI_Gx2yz_S_Py_Px_C1001000003_b+ABX*I_ERI_F2yz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fy2z_Px_Py_Px_C1001000003_b = I_ERI_Gxy2z_S_Py_Px_C1001000003_b+ABX*I_ERI_Fy2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3z_Px_Py_Px_C1001000003_b = I_ERI_Gx3z_S_Py_Px_C1001000003_b+ABX*I_ERI_F3z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3x_Py_Py_Px_C1001000003_b = I_ERI_G3xy_S_Py_Px_C1001000003_b+ABY*I_ERI_F3x_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xy_Py_Py_Px_C1001000003_b = I_ERI_G2x2y_S_Py_Px_C1001000003_b+ABY*I_ERI_F2xy_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xz_Py_Py_Px_C1001000003_b = I_ERI_G2xyz_S_Py_Px_C1001000003_b+ABY*I_ERI_F2xz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2y_Py_Py_Px_C1001000003_b = I_ERI_Gx3y_S_Py_Px_C1001000003_b+ABY*I_ERI_Fx2y_S_Py_Px_C1001000003_b;
  Double I_ERI_Fxyz_Py_Py_Px_C1001000003_b = I_ERI_Gx2yz_S_Py_Px_C1001000003_b+ABY*I_ERI_Fxyz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2z_Py_Py_Px_C1001000003_b = I_ERI_Gxy2z_S_Py_Px_C1001000003_b+ABY*I_ERI_Fx2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3y_Py_Py_Px_C1001000003_b = I_ERI_G4y_S_Py_Px_C1001000003_b+ABY*I_ERI_F3y_S_Py_Px_C1001000003_b;
  Double I_ERI_F2yz_Py_Py_Px_C1001000003_b = I_ERI_G3yz_S_Py_Px_C1001000003_b+ABY*I_ERI_F2yz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fy2z_Py_Py_Px_C1001000003_b = I_ERI_G2y2z_S_Py_Px_C1001000003_b+ABY*I_ERI_Fy2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3z_Py_Py_Px_C1001000003_b = I_ERI_Gy3z_S_Py_Px_C1001000003_b+ABY*I_ERI_F3z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3x_Pz_Py_Px_C1001000003_b = I_ERI_G3xz_S_Py_Px_C1001000003_b+ABZ*I_ERI_F3x_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xy_Pz_Py_Px_C1001000003_b = I_ERI_G2xyz_S_Py_Px_C1001000003_b+ABZ*I_ERI_F2xy_S_Py_Px_C1001000003_b;
  Double I_ERI_F2xz_Pz_Py_Px_C1001000003_b = I_ERI_G2x2z_S_Py_Px_C1001000003_b+ABZ*I_ERI_F2xz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Py_Px_C1001000003_b = I_ERI_Gx2yz_S_Py_Px_C1001000003_b+ABZ*I_ERI_Fx2y_S_Py_Px_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Py_Px_C1001000003_b = I_ERI_Gxy2z_S_Py_Px_C1001000003_b+ABZ*I_ERI_Fxyz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Py_Px_C1001000003_b = I_ERI_Gx3z_S_Py_Px_C1001000003_b+ABZ*I_ERI_Fx2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3y_Pz_Py_Px_C1001000003_b = I_ERI_G3yz_S_Py_Px_C1001000003_b+ABZ*I_ERI_F3y_S_Py_Px_C1001000003_b;
  Double I_ERI_F2yz_Pz_Py_Px_C1001000003_b = I_ERI_G2y2z_S_Py_Px_C1001000003_b+ABZ*I_ERI_F2yz_S_Py_Px_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Py_Px_C1001000003_b = I_ERI_Gy3z_S_Py_Px_C1001000003_b+ABZ*I_ERI_Fy2z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3z_Pz_Py_Px_C1001000003_b = I_ERI_G4z_S_Py_Px_C1001000003_b+ABZ*I_ERI_F3z_S_Py_Px_C1001000003_b;
  Double I_ERI_F3x_Px_Pz_Px_C1001000003_b = I_ERI_G4x_S_Pz_Px_C1001000003_b+ABX*I_ERI_F3x_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xy_Px_Pz_Px_C1001000003_b = I_ERI_G3xy_S_Pz_Px_C1001000003_b+ABX*I_ERI_F2xy_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xz_Px_Pz_Px_C1001000003_b = I_ERI_G3xz_S_Pz_Px_C1001000003_b+ABX*I_ERI_F2xz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2y_Px_Pz_Px_C1001000003_b = I_ERI_G2x2y_S_Pz_Px_C1001000003_b+ABX*I_ERI_Fx2y_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fxyz_Px_Pz_Px_C1001000003_b = I_ERI_G2xyz_S_Pz_Px_C1001000003_b+ABX*I_ERI_Fxyz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2z_Px_Pz_Px_C1001000003_b = I_ERI_G2x2z_S_Pz_Px_C1001000003_b+ABX*I_ERI_Fx2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3y_Px_Pz_Px_C1001000003_b = I_ERI_Gx3y_S_Pz_Px_C1001000003_b+ABX*I_ERI_F3y_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2yz_Px_Pz_Px_C1001000003_b = I_ERI_Gx2yz_S_Pz_Px_C1001000003_b+ABX*I_ERI_F2yz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fy2z_Px_Pz_Px_C1001000003_b = I_ERI_Gxy2z_S_Pz_Px_C1001000003_b+ABX*I_ERI_Fy2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3z_Px_Pz_Px_C1001000003_b = I_ERI_Gx3z_S_Pz_Px_C1001000003_b+ABX*I_ERI_F3z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3x_Py_Pz_Px_C1001000003_b = I_ERI_G3xy_S_Pz_Px_C1001000003_b+ABY*I_ERI_F3x_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xy_Py_Pz_Px_C1001000003_b = I_ERI_G2x2y_S_Pz_Px_C1001000003_b+ABY*I_ERI_F2xy_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xz_Py_Pz_Px_C1001000003_b = I_ERI_G2xyz_S_Pz_Px_C1001000003_b+ABY*I_ERI_F2xz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2y_Py_Pz_Px_C1001000003_b = I_ERI_Gx3y_S_Pz_Px_C1001000003_b+ABY*I_ERI_Fx2y_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fxyz_Py_Pz_Px_C1001000003_b = I_ERI_Gx2yz_S_Pz_Px_C1001000003_b+ABY*I_ERI_Fxyz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2z_Py_Pz_Px_C1001000003_b = I_ERI_Gxy2z_S_Pz_Px_C1001000003_b+ABY*I_ERI_Fx2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3y_Py_Pz_Px_C1001000003_b = I_ERI_G4y_S_Pz_Px_C1001000003_b+ABY*I_ERI_F3y_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2yz_Py_Pz_Px_C1001000003_b = I_ERI_G3yz_S_Pz_Px_C1001000003_b+ABY*I_ERI_F2yz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fy2z_Py_Pz_Px_C1001000003_b = I_ERI_G2y2z_S_Pz_Px_C1001000003_b+ABY*I_ERI_Fy2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3z_Py_Pz_Px_C1001000003_b = I_ERI_Gy3z_S_Pz_Px_C1001000003_b+ABY*I_ERI_F3z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3x_Pz_Pz_Px_C1001000003_b = I_ERI_G3xz_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F3x_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xy_Pz_Pz_Px_C1001000003_b = I_ERI_G2xyz_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F2xy_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2xz_Pz_Pz_Px_C1001000003_b = I_ERI_G2x2z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F2xz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Pz_Px_C1001000003_b = I_ERI_Gx2yz_S_Pz_Px_C1001000003_b+ABZ*I_ERI_Fx2y_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Pz_Px_C1001000003_b = I_ERI_Gxy2z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_Fxyz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Pz_Px_C1001000003_b = I_ERI_Gx3z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_Fx2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3y_Pz_Pz_Px_C1001000003_b = I_ERI_G3yz_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F3y_S_Pz_Px_C1001000003_b;
  Double I_ERI_F2yz_Pz_Pz_Px_C1001000003_b = I_ERI_G2y2z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F2yz_S_Pz_Px_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Pz_Px_C1001000003_b = I_ERI_Gy3z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_Fy2z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3z_Pz_Pz_Px_C1001000003_b = I_ERI_G4z_S_Pz_Px_C1001000003_b+ABZ*I_ERI_F3z_S_Pz_Px_C1001000003_b;
  Double I_ERI_F3x_Px_Px_Py_C1001000003_b = I_ERI_G4x_S_Px_Py_C1001000003_b+ABX*I_ERI_F3x_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xy_Px_Px_Py_C1001000003_b = I_ERI_G3xy_S_Px_Py_C1001000003_b+ABX*I_ERI_F2xy_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xz_Px_Px_Py_C1001000003_b = I_ERI_G3xz_S_Px_Py_C1001000003_b+ABX*I_ERI_F2xz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2y_Px_Px_Py_C1001000003_b = I_ERI_G2x2y_S_Px_Py_C1001000003_b+ABX*I_ERI_Fx2y_S_Px_Py_C1001000003_b;
  Double I_ERI_Fxyz_Px_Px_Py_C1001000003_b = I_ERI_G2xyz_S_Px_Py_C1001000003_b+ABX*I_ERI_Fxyz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2z_Px_Px_Py_C1001000003_b = I_ERI_G2x2z_S_Px_Py_C1001000003_b+ABX*I_ERI_Fx2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3y_Px_Px_Py_C1001000003_b = I_ERI_Gx3y_S_Px_Py_C1001000003_b+ABX*I_ERI_F3y_S_Px_Py_C1001000003_b;
  Double I_ERI_F2yz_Px_Px_Py_C1001000003_b = I_ERI_Gx2yz_S_Px_Py_C1001000003_b+ABX*I_ERI_F2yz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fy2z_Px_Px_Py_C1001000003_b = I_ERI_Gxy2z_S_Px_Py_C1001000003_b+ABX*I_ERI_Fy2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3z_Px_Px_Py_C1001000003_b = I_ERI_Gx3z_S_Px_Py_C1001000003_b+ABX*I_ERI_F3z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3x_Py_Px_Py_C1001000003_b = I_ERI_G3xy_S_Px_Py_C1001000003_b+ABY*I_ERI_F3x_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xy_Py_Px_Py_C1001000003_b = I_ERI_G2x2y_S_Px_Py_C1001000003_b+ABY*I_ERI_F2xy_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xz_Py_Px_Py_C1001000003_b = I_ERI_G2xyz_S_Px_Py_C1001000003_b+ABY*I_ERI_F2xz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2y_Py_Px_Py_C1001000003_b = I_ERI_Gx3y_S_Px_Py_C1001000003_b+ABY*I_ERI_Fx2y_S_Px_Py_C1001000003_b;
  Double I_ERI_Fxyz_Py_Px_Py_C1001000003_b = I_ERI_Gx2yz_S_Px_Py_C1001000003_b+ABY*I_ERI_Fxyz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2z_Py_Px_Py_C1001000003_b = I_ERI_Gxy2z_S_Px_Py_C1001000003_b+ABY*I_ERI_Fx2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3y_Py_Px_Py_C1001000003_b = I_ERI_G4y_S_Px_Py_C1001000003_b+ABY*I_ERI_F3y_S_Px_Py_C1001000003_b;
  Double I_ERI_F2yz_Py_Px_Py_C1001000003_b = I_ERI_G3yz_S_Px_Py_C1001000003_b+ABY*I_ERI_F2yz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fy2z_Py_Px_Py_C1001000003_b = I_ERI_G2y2z_S_Px_Py_C1001000003_b+ABY*I_ERI_Fy2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3z_Py_Px_Py_C1001000003_b = I_ERI_Gy3z_S_Px_Py_C1001000003_b+ABY*I_ERI_F3z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3x_Pz_Px_Py_C1001000003_b = I_ERI_G3xz_S_Px_Py_C1001000003_b+ABZ*I_ERI_F3x_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xy_Pz_Px_Py_C1001000003_b = I_ERI_G2xyz_S_Px_Py_C1001000003_b+ABZ*I_ERI_F2xy_S_Px_Py_C1001000003_b;
  Double I_ERI_F2xz_Pz_Px_Py_C1001000003_b = I_ERI_G2x2z_S_Px_Py_C1001000003_b+ABZ*I_ERI_F2xz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Px_Py_C1001000003_b = I_ERI_Gx2yz_S_Px_Py_C1001000003_b+ABZ*I_ERI_Fx2y_S_Px_Py_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Px_Py_C1001000003_b = I_ERI_Gxy2z_S_Px_Py_C1001000003_b+ABZ*I_ERI_Fxyz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Px_Py_C1001000003_b = I_ERI_Gx3z_S_Px_Py_C1001000003_b+ABZ*I_ERI_Fx2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3y_Pz_Px_Py_C1001000003_b = I_ERI_G3yz_S_Px_Py_C1001000003_b+ABZ*I_ERI_F3y_S_Px_Py_C1001000003_b;
  Double I_ERI_F2yz_Pz_Px_Py_C1001000003_b = I_ERI_G2y2z_S_Px_Py_C1001000003_b+ABZ*I_ERI_F2yz_S_Px_Py_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Px_Py_C1001000003_b = I_ERI_Gy3z_S_Px_Py_C1001000003_b+ABZ*I_ERI_Fy2z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3z_Pz_Px_Py_C1001000003_b = I_ERI_G4z_S_Px_Py_C1001000003_b+ABZ*I_ERI_F3z_S_Px_Py_C1001000003_b;
  Double I_ERI_F3x_Px_Py_Py_C1001000003_b = I_ERI_G4x_S_Py_Py_C1001000003_b+ABX*I_ERI_F3x_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xy_Px_Py_Py_C1001000003_b = I_ERI_G3xy_S_Py_Py_C1001000003_b+ABX*I_ERI_F2xy_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xz_Px_Py_Py_C1001000003_b = I_ERI_G3xz_S_Py_Py_C1001000003_b+ABX*I_ERI_F2xz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2y_Px_Py_Py_C1001000003_b = I_ERI_G2x2y_S_Py_Py_C1001000003_b+ABX*I_ERI_Fx2y_S_Py_Py_C1001000003_b;
  Double I_ERI_Fxyz_Px_Py_Py_C1001000003_b = I_ERI_G2xyz_S_Py_Py_C1001000003_b+ABX*I_ERI_Fxyz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2z_Px_Py_Py_C1001000003_b = I_ERI_G2x2z_S_Py_Py_C1001000003_b+ABX*I_ERI_Fx2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3y_Px_Py_Py_C1001000003_b = I_ERI_Gx3y_S_Py_Py_C1001000003_b+ABX*I_ERI_F3y_S_Py_Py_C1001000003_b;
  Double I_ERI_F2yz_Px_Py_Py_C1001000003_b = I_ERI_Gx2yz_S_Py_Py_C1001000003_b+ABX*I_ERI_F2yz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fy2z_Px_Py_Py_C1001000003_b = I_ERI_Gxy2z_S_Py_Py_C1001000003_b+ABX*I_ERI_Fy2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3z_Px_Py_Py_C1001000003_b = I_ERI_Gx3z_S_Py_Py_C1001000003_b+ABX*I_ERI_F3z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3x_Py_Py_Py_C1001000003_b = I_ERI_G3xy_S_Py_Py_C1001000003_b+ABY*I_ERI_F3x_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xy_Py_Py_Py_C1001000003_b = I_ERI_G2x2y_S_Py_Py_C1001000003_b+ABY*I_ERI_F2xy_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xz_Py_Py_Py_C1001000003_b = I_ERI_G2xyz_S_Py_Py_C1001000003_b+ABY*I_ERI_F2xz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2y_Py_Py_Py_C1001000003_b = I_ERI_Gx3y_S_Py_Py_C1001000003_b+ABY*I_ERI_Fx2y_S_Py_Py_C1001000003_b;
  Double I_ERI_Fxyz_Py_Py_Py_C1001000003_b = I_ERI_Gx2yz_S_Py_Py_C1001000003_b+ABY*I_ERI_Fxyz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2z_Py_Py_Py_C1001000003_b = I_ERI_Gxy2z_S_Py_Py_C1001000003_b+ABY*I_ERI_Fx2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3y_Py_Py_Py_C1001000003_b = I_ERI_G4y_S_Py_Py_C1001000003_b+ABY*I_ERI_F3y_S_Py_Py_C1001000003_b;
  Double I_ERI_F2yz_Py_Py_Py_C1001000003_b = I_ERI_G3yz_S_Py_Py_C1001000003_b+ABY*I_ERI_F2yz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fy2z_Py_Py_Py_C1001000003_b = I_ERI_G2y2z_S_Py_Py_C1001000003_b+ABY*I_ERI_Fy2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3z_Py_Py_Py_C1001000003_b = I_ERI_Gy3z_S_Py_Py_C1001000003_b+ABY*I_ERI_F3z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3x_Pz_Py_Py_C1001000003_b = I_ERI_G3xz_S_Py_Py_C1001000003_b+ABZ*I_ERI_F3x_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xy_Pz_Py_Py_C1001000003_b = I_ERI_G2xyz_S_Py_Py_C1001000003_b+ABZ*I_ERI_F2xy_S_Py_Py_C1001000003_b;
  Double I_ERI_F2xz_Pz_Py_Py_C1001000003_b = I_ERI_G2x2z_S_Py_Py_C1001000003_b+ABZ*I_ERI_F2xz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Py_Py_C1001000003_b = I_ERI_Gx2yz_S_Py_Py_C1001000003_b+ABZ*I_ERI_Fx2y_S_Py_Py_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Py_Py_C1001000003_b = I_ERI_Gxy2z_S_Py_Py_C1001000003_b+ABZ*I_ERI_Fxyz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Py_Py_C1001000003_b = I_ERI_Gx3z_S_Py_Py_C1001000003_b+ABZ*I_ERI_Fx2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3y_Pz_Py_Py_C1001000003_b = I_ERI_G3yz_S_Py_Py_C1001000003_b+ABZ*I_ERI_F3y_S_Py_Py_C1001000003_b;
  Double I_ERI_F2yz_Pz_Py_Py_C1001000003_b = I_ERI_G2y2z_S_Py_Py_C1001000003_b+ABZ*I_ERI_F2yz_S_Py_Py_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Py_Py_C1001000003_b = I_ERI_Gy3z_S_Py_Py_C1001000003_b+ABZ*I_ERI_Fy2z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3z_Pz_Py_Py_C1001000003_b = I_ERI_G4z_S_Py_Py_C1001000003_b+ABZ*I_ERI_F3z_S_Py_Py_C1001000003_b;
  Double I_ERI_F3x_Px_Pz_Py_C1001000003_b = I_ERI_G4x_S_Pz_Py_C1001000003_b+ABX*I_ERI_F3x_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xy_Px_Pz_Py_C1001000003_b = I_ERI_G3xy_S_Pz_Py_C1001000003_b+ABX*I_ERI_F2xy_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xz_Px_Pz_Py_C1001000003_b = I_ERI_G3xz_S_Pz_Py_C1001000003_b+ABX*I_ERI_F2xz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2y_Px_Pz_Py_C1001000003_b = I_ERI_G2x2y_S_Pz_Py_C1001000003_b+ABX*I_ERI_Fx2y_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fxyz_Px_Pz_Py_C1001000003_b = I_ERI_G2xyz_S_Pz_Py_C1001000003_b+ABX*I_ERI_Fxyz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2z_Px_Pz_Py_C1001000003_b = I_ERI_G2x2z_S_Pz_Py_C1001000003_b+ABX*I_ERI_Fx2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3y_Px_Pz_Py_C1001000003_b = I_ERI_Gx3y_S_Pz_Py_C1001000003_b+ABX*I_ERI_F3y_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2yz_Px_Pz_Py_C1001000003_b = I_ERI_Gx2yz_S_Pz_Py_C1001000003_b+ABX*I_ERI_F2yz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fy2z_Px_Pz_Py_C1001000003_b = I_ERI_Gxy2z_S_Pz_Py_C1001000003_b+ABX*I_ERI_Fy2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3z_Px_Pz_Py_C1001000003_b = I_ERI_Gx3z_S_Pz_Py_C1001000003_b+ABX*I_ERI_F3z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3x_Py_Pz_Py_C1001000003_b = I_ERI_G3xy_S_Pz_Py_C1001000003_b+ABY*I_ERI_F3x_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xy_Py_Pz_Py_C1001000003_b = I_ERI_G2x2y_S_Pz_Py_C1001000003_b+ABY*I_ERI_F2xy_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xz_Py_Pz_Py_C1001000003_b = I_ERI_G2xyz_S_Pz_Py_C1001000003_b+ABY*I_ERI_F2xz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2y_Py_Pz_Py_C1001000003_b = I_ERI_Gx3y_S_Pz_Py_C1001000003_b+ABY*I_ERI_Fx2y_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fxyz_Py_Pz_Py_C1001000003_b = I_ERI_Gx2yz_S_Pz_Py_C1001000003_b+ABY*I_ERI_Fxyz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2z_Py_Pz_Py_C1001000003_b = I_ERI_Gxy2z_S_Pz_Py_C1001000003_b+ABY*I_ERI_Fx2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3y_Py_Pz_Py_C1001000003_b = I_ERI_G4y_S_Pz_Py_C1001000003_b+ABY*I_ERI_F3y_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2yz_Py_Pz_Py_C1001000003_b = I_ERI_G3yz_S_Pz_Py_C1001000003_b+ABY*I_ERI_F2yz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fy2z_Py_Pz_Py_C1001000003_b = I_ERI_G2y2z_S_Pz_Py_C1001000003_b+ABY*I_ERI_Fy2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3z_Py_Pz_Py_C1001000003_b = I_ERI_Gy3z_S_Pz_Py_C1001000003_b+ABY*I_ERI_F3z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3x_Pz_Pz_Py_C1001000003_b = I_ERI_G3xz_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F3x_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xy_Pz_Pz_Py_C1001000003_b = I_ERI_G2xyz_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F2xy_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2xz_Pz_Pz_Py_C1001000003_b = I_ERI_G2x2z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F2xz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Pz_Py_C1001000003_b = I_ERI_Gx2yz_S_Pz_Py_C1001000003_b+ABZ*I_ERI_Fx2y_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Pz_Py_C1001000003_b = I_ERI_Gxy2z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_Fxyz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Pz_Py_C1001000003_b = I_ERI_Gx3z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_Fx2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3y_Pz_Pz_Py_C1001000003_b = I_ERI_G3yz_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F3y_S_Pz_Py_C1001000003_b;
  Double I_ERI_F2yz_Pz_Pz_Py_C1001000003_b = I_ERI_G2y2z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F2yz_S_Pz_Py_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Pz_Py_C1001000003_b = I_ERI_Gy3z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_Fy2z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3z_Pz_Pz_Py_C1001000003_b = I_ERI_G4z_S_Pz_Py_C1001000003_b+ABZ*I_ERI_F3z_S_Pz_Py_C1001000003_b;
  Double I_ERI_F3x_Px_Px_Pz_C1001000003_b = I_ERI_G4x_S_Px_Pz_C1001000003_b+ABX*I_ERI_F3x_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xy_Px_Px_Pz_C1001000003_b = I_ERI_G3xy_S_Px_Pz_C1001000003_b+ABX*I_ERI_F2xy_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xz_Px_Px_Pz_C1001000003_b = I_ERI_G3xz_S_Px_Pz_C1001000003_b+ABX*I_ERI_F2xz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Px_Px_Pz_C1001000003_b = I_ERI_G2x2y_S_Px_Pz_C1001000003_b+ABX*I_ERI_Fx2y_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Px_Px_Pz_C1001000003_b = I_ERI_G2xyz_S_Px_Pz_C1001000003_b+ABX*I_ERI_Fxyz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Px_Px_Pz_C1001000003_b = I_ERI_G2x2z_S_Px_Pz_C1001000003_b+ABX*I_ERI_Fx2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3y_Px_Px_Pz_C1001000003_b = I_ERI_Gx3y_S_Px_Pz_C1001000003_b+ABX*I_ERI_F3y_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2yz_Px_Px_Pz_C1001000003_b = I_ERI_Gx2yz_S_Px_Pz_C1001000003_b+ABX*I_ERI_F2yz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Px_Px_Pz_C1001000003_b = I_ERI_Gxy2z_S_Px_Pz_C1001000003_b+ABX*I_ERI_Fy2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3z_Px_Px_Pz_C1001000003_b = I_ERI_Gx3z_S_Px_Pz_C1001000003_b+ABX*I_ERI_F3z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3x_Py_Px_Pz_C1001000003_b = I_ERI_G3xy_S_Px_Pz_C1001000003_b+ABY*I_ERI_F3x_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xy_Py_Px_Pz_C1001000003_b = I_ERI_G2x2y_S_Px_Pz_C1001000003_b+ABY*I_ERI_F2xy_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xz_Py_Px_Pz_C1001000003_b = I_ERI_G2xyz_S_Px_Pz_C1001000003_b+ABY*I_ERI_F2xz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Py_Px_Pz_C1001000003_b = I_ERI_Gx3y_S_Px_Pz_C1001000003_b+ABY*I_ERI_Fx2y_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Py_Px_Pz_C1001000003_b = I_ERI_Gx2yz_S_Px_Pz_C1001000003_b+ABY*I_ERI_Fxyz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Py_Px_Pz_C1001000003_b = I_ERI_Gxy2z_S_Px_Pz_C1001000003_b+ABY*I_ERI_Fx2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3y_Py_Px_Pz_C1001000003_b = I_ERI_G4y_S_Px_Pz_C1001000003_b+ABY*I_ERI_F3y_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2yz_Py_Px_Pz_C1001000003_b = I_ERI_G3yz_S_Px_Pz_C1001000003_b+ABY*I_ERI_F2yz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Py_Px_Pz_C1001000003_b = I_ERI_G2y2z_S_Px_Pz_C1001000003_b+ABY*I_ERI_Fy2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3z_Py_Px_Pz_C1001000003_b = I_ERI_Gy3z_S_Px_Pz_C1001000003_b+ABY*I_ERI_F3z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3x_Pz_Px_Pz_C1001000003_b = I_ERI_G3xz_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F3x_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xy_Pz_Px_Pz_C1001000003_b = I_ERI_G2xyz_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F2xy_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2xz_Pz_Px_Pz_C1001000003_b = I_ERI_G2x2z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F2xz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Px_Pz_C1001000003_b = I_ERI_Gx2yz_S_Px_Pz_C1001000003_b+ABZ*I_ERI_Fx2y_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Px_Pz_C1001000003_b = I_ERI_Gxy2z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_Fxyz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Px_Pz_C1001000003_b = I_ERI_Gx3z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_Fx2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3y_Pz_Px_Pz_C1001000003_b = I_ERI_G3yz_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F3y_S_Px_Pz_C1001000003_b;
  Double I_ERI_F2yz_Pz_Px_Pz_C1001000003_b = I_ERI_G2y2z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F2yz_S_Px_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Px_Pz_C1001000003_b = I_ERI_Gy3z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_Fy2z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3z_Pz_Px_Pz_C1001000003_b = I_ERI_G4z_S_Px_Pz_C1001000003_b+ABZ*I_ERI_F3z_S_Px_Pz_C1001000003_b;
  Double I_ERI_F3x_Px_Py_Pz_C1001000003_b = I_ERI_G4x_S_Py_Pz_C1001000003_b+ABX*I_ERI_F3x_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xy_Px_Py_Pz_C1001000003_b = I_ERI_G3xy_S_Py_Pz_C1001000003_b+ABX*I_ERI_F2xy_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xz_Px_Py_Pz_C1001000003_b = I_ERI_G3xz_S_Py_Pz_C1001000003_b+ABX*I_ERI_F2xz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Px_Py_Pz_C1001000003_b = I_ERI_G2x2y_S_Py_Pz_C1001000003_b+ABX*I_ERI_Fx2y_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Px_Py_Pz_C1001000003_b = I_ERI_G2xyz_S_Py_Pz_C1001000003_b+ABX*I_ERI_Fxyz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Px_Py_Pz_C1001000003_b = I_ERI_G2x2z_S_Py_Pz_C1001000003_b+ABX*I_ERI_Fx2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3y_Px_Py_Pz_C1001000003_b = I_ERI_Gx3y_S_Py_Pz_C1001000003_b+ABX*I_ERI_F3y_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2yz_Px_Py_Pz_C1001000003_b = I_ERI_Gx2yz_S_Py_Pz_C1001000003_b+ABX*I_ERI_F2yz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Px_Py_Pz_C1001000003_b = I_ERI_Gxy2z_S_Py_Pz_C1001000003_b+ABX*I_ERI_Fy2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3z_Px_Py_Pz_C1001000003_b = I_ERI_Gx3z_S_Py_Pz_C1001000003_b+ABX*I_ERI_F3z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3x_Py_Py_Pz_C1001000003_b = I_ERI_G3xy_S_Py_Pz_C1001000003_b+ABY*I_ERI_F3x_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xy_Py_Py_Pz_C1001000003_b = I_ERI_G2x2y_S_Py_Pz_C1001000003_b+ABY*I_ERI_F2xy_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xz_Py_Py_Pz_C1001000003_b = I_ERI_G2xyz_S_Py_Pz_C1001000003_b+ABY*I_ERI_F2xz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Py_Py_Pz_C1001000003_b = I_ERI_Gx3y_S_Py_Pz_C1001000003_b+ABY*I_ERI_Fx2y_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Py_Py_Pz_C1001000003_b = I_ERI_Gx2yz_S_Py_Pz_C1001000003_b+ABY*I_ERI_Fxyz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Py_Py_Pz_C1001000003_b = I_ERI_Gxy2z_S_Py_Pz_C1001000003_b+ABY*I_ERI_Fx2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3y_Py_Py_Pz_C1001000003_b = I_ERI_G4y_S_Py_Pz_C1001000003_b+ABY*I_ERI_F3y_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2yz_Py_Py_Pz_C1001000003_b = I_ERI_G3yz_S_Py_Pz_C1001000003_b+ABY*I_ERI_F2yz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Py_Py_Pz_C1001000003_b = I_ERI_G2y2z_S_Py_Pz_C1001000003_b+ABY*I_ERI_Fy2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3z_Py_Py_Pz_C1001000003_b = I_ERI_Gy3z_S_Py_Pz_C1001000003_b+ABY*I_ERI_F3z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3x_Pz_Py_Pz_C1001000003_b = I_ERI_G3xz_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F3x_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xy_Pz_Py_Pz_C1001000003_b = I_ERI_G2xyz_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F2xy_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2xz_Pz_Py_Pz_C1001000003_b = I_ERI_G2x2z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F2xz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Py_Pz_C1001000003_b = I_ERI_Gx2yz_S_Py_Pz_C1001000003_b+ABZ*I_ERI_Fx2y_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Py_Pz_C1001000003_b = I_ERI_Gxy2z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_Fxyz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Py_Pz_C1001000003_b = I_ERI_Gx3z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_Fx2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3y_Pz_Py_Pz_C1001000003_b = I_ERI_G3yz_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F3y_S_Py_Pz_C1001000003_b;
  Double I_ERI_F2yz_Pz_Py_Pz_C1001000003_b = I_ERI_G2y2z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F2yz_S_Py_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Py_Pz_C1001000003_b = I_ERI_Gy3z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_Fy2z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3z_Pz_Py_Pz_C1001000003_b = I_ERI_G4z_S_Py_Pz_C1001000003_b+ABZ*I_ERI_F3z_S_Py_Pz_C1001000003_b;
  Double I_ERI_F3x_Px_Pz_Pz_C1001000003_b = I_ERI_G4x_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F3x_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xy_Px_Pz_Pz_C1001000003_b = I_ERI_G3xy_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F2xy_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xz_Px_Pz_Pz_C1001000003_b = I_ERI_G3xz_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F2xz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Px_Pz_Pz_C1001000003_b = I_ERI_G2x2y_S_Pz_Pz_C1001000003_b+ABX*I_ERI_Fx2y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Px_Pz_Pz_C1001000003_b = I_ERI_G2xyz_S_Pz_Pz_C1001000003_b+ABX*I_ERI_Fxyz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Px_Pz_Pz_C1001000003_b = I_ERI_G2x2z_S_Pz_Pz_C1001000003_b+ABX*I_ERI_Fx2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3y_Px_Pz_Pz_C1001000003_b = I_ERI_Gx3y_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F3y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2yz_Px_Pz_Pz_C1001000003_b = I_ERI_Gx2yz_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F2yz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Px_Pz_Pz_C1001000003_b = I_ERI_Gxy2z_S_Pz_Pz_C1001000003_b+ABX*I_ERI_Fy2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3z_Px_Pz_Pz_C1001000003_b = I_ERI_Gx3z_S_Pz_Pz_C1001000003_b+ABX*I_ERI_F3z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3x_Py_Pz_Pz_C1001000003_b = I_ERI_G3xy_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F3x_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xy_Py_Pz_Pz_C1001000003_b = I_ERI_G2x2y_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F2xy_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xz_Py_Pz_Pz_C1001000003_b = I_ERI_G2xyz_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F2xz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Py_Pz_Pz_C1001000003_b = I_ERI_Gx3y_S_Pz_Pz_C1001000003_b+ABY*I_ERI_Fx2y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Py_Pz_Pz_C1001000003_b = I_ERI_Gx2yz_S_Pz_Pz_C1001000003_b+ABY*I_ERI_Fxyz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Py_Pz_Pz_C1001000003_b = I_ERI_Gxy2z_S_Pz_Pz_C1001000003_b+ABY*I_ERI_Fx2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3y_Py_Pz_Pz_C1001000003_b = I_ERI_G4y_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F3y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2yz_Py_Pz_Pz_C1001000003_b = I_ERI_G3yz_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F2yz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Py_Pz_Pz_C1001000003_b = I_ERI_G2y2z_S_Pz_Pz_C1001000003_b+ABY*I_ERI_Fy2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3z_Py_Pz_Pz_C1001000003_b = I_ERI_Gy3z_S_Pz_Pz_C1001000003_b+ABY*I_ERI_F3z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3x_Pz_Pz_Pz_C1001000003_b = I_ERI_G3xz_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F3x_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xy_Pz_Pz_Pz_C1001000003_b = I_ERI_G2xyz_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F2xy_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2xz_Pz_Pz_Pz_C1001000003_b = I_ERI_G2x2z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F2xz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2y_Pz_Pz_Pz_C1001000003_b = I_ERI_Gx2yz_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_Fx2y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fxyz_Pz_Pz_Pz_C1001000003_b = I_ERI_Gxy2z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_Fxyz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fx2z_Pz_Pz_Pz_C1001000003_b = I_ERI_Gx3z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_Fx2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3y_Pz_Pz_Pz_C1001000003_b = I_ERI_G3yz_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F3y_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F2yz_Pz_Pz_Pz_C1001000003_b = I_ERI_G2y2z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F2yz_S_Pz_Pz_C1001000003_b;
  Double I_ERI_Fy2z_Pz_Pz_Pz_C1001000003_b = I_ERI_Gy3z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_Fy2z_S_Pz_Pz_C1001000003_b;
  Double I_ERI_F3z_Pz_Pz_Pz_C1001000003_b = I_ERI_G4z_S_Pz_Pz_C1001000003_b+ABZ*I_ERI_F3z_S_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_F3x_Px_Px_S_C1000003_b;
  abcd[1] = 2.0E0*I_ERI_F2xy_Px_Px_S_C1000003_b;
  abcd[2] = 2.0E0*I_ERI_F2xz_Px_Px_S_C1000003_b;
  abcd[3] = 2.0E0*I_ERI_Fx2y_Px_Px_S_C1000003_b;
  abcd[4] = 2.0E0*I_ERI_Fxyz_Px_Px_S_C1000003_b;
  abcd[5] = 2.0E0*I_ERI_Fx2z_Px_Px_S_C1000003_b;
  abcd[6] = 2.0E0*I_ERI_F3y_Px_Px_S_C1000003_b;
  abcd[7] = 2.0E0*I_ERI_F2yz_Px_Px_S_C1000003_b;
  abcd[8] = 2.0E0*I_ERI_Fy2z_Px_Px_S_C1000003_b;
  abcd[9] = 2.0E0*I_ERI_F3z_Px_Px_S_C1000003_b;
  abcd[10] = 2.0E0*I_ERI_F3x_Px_Py_S_C1000003_b;
  abcd[11] = 2.0E0*I_ERI_F2xy_Px_Py_S_C1000003_b;
  abcd[12] = 2.0E0*I_ERI_F2xz_Px_Py_S_C1000003_b;
  abcd[13] = 2.0E0*I_ERI_Fx2y_Px_Py_S_C1000003_b;
  abcd[14] = 2.0E0*I_ERI_Fxyz_Px_Py_S_C1000003_b;
  abcd[15] = 2.0E0*I_ERI_Fx2z_Px_Py_S_C1000003_b;
  abcd[16] = 2.0E0*I_ERI_F3y_Px_Py_S_C1000003_b;
  abcd[17] = 2.0E0*I_ERI_F2yz_Px_Py_S_C1000003_b;
  abcd[18] = 2.0E0*I_ERI_Fy2z_Px_Py_S_C1000003_b;
  abcd[19] = 2.0E0*I_ERI_F3z_Px_Py_S_C1000003_b;
  abcd[20] = 2.0E0*I_ERI_F3x_Px_Pz_S_C1000003_b;
  abcd[21] = 2.0E0*I_ERI_F2xy_Px_Pz_S_C1000003_b;
  abcd[22] = 2.0E0*I_ERI_F2xz_Px_Pz_S_C1000003_b;
  abcd[23] = 2.0E0*I_ERI_Fx2y_Px_Pz_S_C1000003_b;
  abcd[24] = 2.0E0*I_ERI_Fxyz_Px_Pz_S_C1000003_b;
  abcd[25] = 2.0E0*I_ERI_Fx2z_Px_Pz_S_C1000003_b;
  abcd[26] = 2.0E0*I_ERI_F3y_Px_Pz_S_C1000003_b;
  abcd[27] = 2.0E0*I_ERI_F2yz_Px_Pz_S_C1000003_b;
  abcd[28] = 2.0E0*I_ERI_Fy2z_Px_Pz_S_C1000003_b;
  abcd[29] = 2.0E0*I_ERI_F3z_Px_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[30] = 2.0E0*I_ERI_F3x_Px_Px_Px_C1001000003_b;
  abcd[31] = 2.0E0*I_ERI_F2xy_Px_Px_Px_C1001000003_b;
  abcd[32] = 2.0E0*I_ERI_F2xz_Px_Px_Px_C1001000003_b;
  abcd[33] = 2.0E0*I_ERI_Fx2y_Px_Px_Px_C1001000003_b;
  abcd[34] = 2.0E0*I_ERI_Fxyz_Px_Px_Px_C1001000003_b;
  abcd[35] = 2.0E0*I_ERI_Fx2z_Px_Px_Px_C1001000003_b;
  abcd[36] = 2.0E0*I_ERI_F3y_Px_Px_Px_C1001000003_b;
  abcd[37] = 2.0E0*I_ERI_F2yz_Px_Px_Px_C1001000003_b;
  abcd[38] = 2.0E0*I_ERI_Fy2z_Px_Px_Px_C1001000003_b;
  abcd[39] = 2.0E0*I_ERI_F3z_Px_Px_Px_C1001000003_b;
  abcd[40] = 2.0E0*I_ERI_F3x_Px_Py_Px_C1001000003_b;
  abcd[41] = 2.0E0*I_ERI_F2xy_Px_Py_Px_C1001000003_b;
  abcd[42] = 2.0E0*I_ERI_F2xz_Px_Py_Px_C1001000003_b;
  abcd[43] = 2.0E0*I_ERI_Fx2y_Px_Py_Px_C1001000003_b;
  abcd[44] = 2.0E0*I_ERI_Fxyz_Px_Py_Px_C1001000003_b;
  abcd[45] = 2.0E0*I_ERI_Fx2z_Px_Py_Px_C1001000003_b;
  abcd[46] = 2.0E0*I_ERI_F3y_Px_Py_Px_C1001000003_b;
  abcd[47] = 2.0E0*I_ERI_F2yz_Px_Py_Px_C1001000003_b;
  abcd[48] = 2.0E0*I_ERI_Fy2z_Px_Py_Px_C1001000003_b;
  abcd[49] = 2.0E0*I_ERI_F3z_Px_Py_Px_C1001000003_b;
  abcd[50] = 2.0E0*I_ERI_F3x_Px_Pz_Px_C1001000003_b;
  abcd[51] = 2.0E0*I_ERI_F2xy_Px_Pz_Px_C1001000003_b;
  abcd[52] = 2.0E0*I_ERI_F2xz_Px_Pz_Px_C1001000003_b;
  abcd[53] = 2.0E0*I_ERI_Fx2y_Px_Pz_Px_C1001000003_b;
  abcd[54] = 2.0E0*I_ERI_Fxyz_Px_Pz_Px_C1001000003_b;
  abcd[55] = 2.0E0*I_ERI_Fx2z_Px_Pz_Px_C1001000003_b;
  abcd[56] = 2.0E0*I_ERI_F3y_Px_Pz_Px_C1001000003_b;
  abcd[57] = 2.0E0*I_ERI_F2yz_Px_Pz_Px_C1001000003_b;
  abcd[58] = 2.0E0*I_ERI_Fy2z_Px_Pz_Px_C1001000003_b;
  abcd[59] = 2.0E0*I_ERI_F3z_Px_Pz_Px_C1001000003_b;
  abcd[60] = 2.0E0*I_ERI_F3x_Px_Px_Py_C1001000003_b;
  abcd[61] = 2.0E0*I_ERI_F2xy_Px_Px_Py_C1001000003_b;
  abcd[62] = 2.0E0*I_ERI_F2xz_Px_Px_Py_C1001000003_b;
  abcd[63] = 2.0E0*I_ERI_Fx2y_Px_Px_Py_C1001000003_b;
  abcd[64] = 2.0E0*I_ERI_Fxyz_Px_Px_Py_C1001000003_b;
  abcd[65] = 2.0E0*I_ERI_Fx2z_Px_Px_Py_C1001000003_b;
  abcd[66] = 2.0E0*I_ERI_F3y_Px_Px_Py_C1001000003_b;
  abcd[67] = 2.0E0*I_ERI_F2yz_Px_Px_Py_C1001000003_b;
  abcd[68] = 2.0E0*I_ERI_Fy2z_Px_Px_Py_C1001000003_b;
  abcd[69] = 2.0E0*I_ERI_F3z_Px_Px_Py_C1001000003_b;
  abcd[70] = 2.0E0*I_ERI_F3x_Px_Py_Py_C1001000003_b;
  abcd[71] = 2.0E0*I_ERI_F2xy_Px_Py_Py_C1001000003_b;
  abcd[72] = 2.0E0*I_ERI_F2xz_Px_Py_Py_C1001000003_b;
  abcd[73] = 2.0E0*I_ERI_Fx2y_Px_Py_Py_C1001000003_b;
  abcd[74] = 2.0E0*I_ERI_Fxyz_Px_Py_Py_C1001000003_b;
  abcd[75] = 2.0E0*I_ERI_Fx2z_Px_Py_Py_C1001000003_b;
  abcd[76] = 2.0E0*I_ERI_F3y_Px_Py_Py_C1001000003_b;
  abcd[77] = 2.0E0*I_ERI_F2yz_Px_Py_Py_C1001000003_b;
  abcd[78] = 2.0E0*I_ERI_Fy2z_Px_Py_Py_C1001000003_b;
  abcd[79] = 2.0E0*I_ERI_F3z_Px_Py_Py_C1001000003_b;
  abcd[80] = 2.0E0*I_ERI_F3x_Px_Pz_Py_C1001000003_b;
  abcd[81] = 2.0E0*I_ERI_F2xy_Px_Pz_Py_C1001000003_b;
  abcd[82] = 2.0E0*I_ERI_F2xz_Px_Pz_Py_C1001000003_b;
  abcd[83] = 2.0E0*I_ERI_Fx2y_Px_Pz_Py_C1001000003_b;
  abcd[84] = 2.0E0*I_ERI_Fxyz_Px_Pz_Py_C1001000003_b;
  abcd[85] = 2.0E0*I_ERI_Fx2z_Px_Pz_Py_C1001000003_b;
  abcd[86] = 2.0E0*I_ERI_F3y_Px_Pz_Py_C1001000003_b;
  abcd[87] = 2.0E0*I_ERI_F2yz_Px_Pz_Py_C1001000003_b;
  abcd[88] = 2.0E0*I_ERI_Fy2z_Px_Pz_Py_C1001000003_b;
  abcd[89] = 2.0E0*I_ERI_F3z_Px_Pz_Py_C1001000003_b;
  abcd[90] = 2.0E0*I_ERI_F3x_Px_Px_Pz_C1001000003_b;
  abcd[91] = 2.0E0*I_ERI_F2xy_Px_Px_Pz_C1001000003_b;
  abcd[92] = 2.0E0*I_ERI_F2xz_Px_Px_Pz_C1001000003_b;
  abcd[93] = 2.0E0*I_ERI_Fx2y_Px_Px_Pz_C1001000003_b;
  abcd[94] = 2.0E0*I_ERI_Fxyz_Px_Px_Pz_C1001000003_b;
  abcd[95] = 2.0E0*I_ERI_Fx2z_Px_Px_Pz_C1001000003_b;
  abcd[96] = 2.0E0*I_ERI_F3y_Px_Px_Pz_C1001000003_b;
  abcd[97] = 2.0E0*I_ERI_F2yz_Px_Px_Pz_C1001000003_b;
  abcd[98] = 2.0E0*I_ERI_Fy2z_Px_Px_Pz_C1001000003_b;
  abcd[99] = 2.0E0*I_ERI_F3z_Px_Px_Pz_C1001000003_b;
  abcd[100] = 2.0E0*I_ERI_F3x_Px_Py_Pz_C1001000003_b;
  abcd[101] = 2.0E0*I_ERI_F2xy_Px_Py_Pz_C1001000003_b;
  abcd[102] = 2.0E0*I_ERI_F2xz_Px_Py_Pz_C1001000003_b;
  abcd[103] = 2.0E0*I_ERI_Fx2y_Px_Py_Pz_C1001000003_b;
  abcd[104] = 2.0E0*I_ERI_Fxyz_Px_Py_Pz_C1001000003_b;
  abcd[105] = 2.0E0*I_ERI_Fx2z_Px_Py_Pz_C1001000003_b;
  abcd[106] = 2.0E0*I_ERI_F3y_Px_Py_Pz_C1001000003_b;
  abcd[107] = 2.0E0*I_ERI_F2yz_Px_Py_Pz_C1001000003_b;
  abcd[108] = 2.0E0*I_ERI_Fy2z_Px_Py_Pz_C1001000003_b;
  abcd[109] = 2.0E0*I_ERI_F3z_Px_Py_Pz_C1001000003_b;
  abcd[110] = 2.0E0*I_ERI_F3x_Px_Pz_Pz_C1001000003_b;
  abcd[111] = 2.0E0*I_ERI_F2xy_Px_Pz_Pz_C1001000003_b;
  abcd[112] = 2.0E0*I_ERI_F2xz_Px_Pz_Pz_C1001000003_b;
  abcd[113] = 2.0E0*I_ERI_Fx2y_Px_Pz_Pz_C1001000003_b;
  abcd[114] = 2.0E0*I_ERI_Fxyz_Px_Pz_Pz_C1001000003_b;
  abcd[115] = 2.0E0*I_ERI_Fx2z_Px_Pz_Pz_C1001000003_b;
  abcd[116] = 2.0E0*I_ERI_F3y_Px_Pz_Pz_C1001000003_b;
  abcd[117] = 2.0E0*I_ERI_F2yz_Px_Pz_Pz_C1001000003_b;
  abcd[118] = 2.0E0*I_ERI_Fy2z_Px_Pz_Pz_C1001000003_b;
  abcd[119] = 2.0E0*I_ERI_F3z_Px_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[120] = 2.0E0*I_ERI_F3x_Py_Px_S_C1000003_b;
  abcd[121] = 2.0E0*I_ERI_F2xy_Py_Px_S_C1000003_b;
  abcd[122] = 2.0E0*I_ERI_F2xz_Py_Px_S_C1000003_b;
  abcd[123] = 2.0E0*I_ERI_Fx2y_Py_Px_S_C1000003_b;
  abcd[124] = 2.0E0*I_ERI_Fxyz_Py_Px_S_C1000003_b;
  abcd[125] = 2.0E0*I_ERI_Fx2z_Py_Px_S_C1000003_b;
  abcd[126] = 2.0E0*I_ERI_F3y_Py_Px_S_C1000003_b;
  abcd[127] = 2.0E0*I_ERI_F2yz_Py_Px_S_C1000003_b;
  abcd[128] = 2.0E0*I_ERI_Fy2z_Py_Px_S_C1000003_b;
  abcd[129] = 2.0E0*I_ERI_F3z_Py_Px_S_C1000003_b;
  abcd[130] = 2.0E0*I_ERI_F3x_Py_Py_S_C1000003_b;
  abcd[131] = 2.0E0*I_ERI_F2xy_Py_Py_S_C1000003_b;
  abcd[132] = 2.0E0*I_ERI_F2xz_Py_Py_S_C1000003_b;
  abcd[133] = 2.0E0*I_ERI_Fx2y_Py_Py_S_C1000003_b;
  abcd[134] = 2.0E0*I_ERI_Fxyz_Py_Py_S_C1000003_b;
  abcd[135] = 2.0E0*I_ERI_Fx2z_Py_Py_S_C1000003_b;
  abcd[136] = 2.0E0*I_ERI_F3y_Py_Py_S_C1000003_b;
  abcd[137] = 2.0E0*I_ERI_F2yz_Py_Py_S_C1000003_b;
  abcd[138] = 2.0E0*I_ERI_Fy2z_Py_Py_S_C1000003_b;
  abcd[139] = 2.0E0*I_ERI_F3z_Py_Py_S_C1000003_b;
  abcd[140] = 2.0E0*I_ERI_F3x_Py_Pz_S_C1000003_b;
  abcd[141] = 2.0E0*I_ERI_F2xy_Py_Pz_S_C1000003_b;
  abcd[142] = 2.0E0*I_ERI_F2xz_Py_Pz_S_C1000003_b;
  abcd[143] = 2.0E0*I_ERI_Fx2y_Py_Pz_S_C1000003_b;
  abcd[144] = 2.0E0*I_ERI_Fxyz_Py_Pz_S_C1000003_b;
  abcd[145] = 2.0E0*I_ERI_Fx2z_Py_Pz_S_C1000003_b;
  abcd[146] = 2.0E0*I_ERI_F3y_Py_Pz_S_C1000003_b;
  abcd[147] = 2.0E0*I_ERI_F2yz_Py_Pz_S_C1000003_b;
  abcd[148] = 2.0E0*I_ERI_Fy2z_Py_Pz_S_C1000003_b;
  abcd[149] = 2.0E0*I_ERI_F3z_Py_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[150] = 2.0E0*I_ERI_F3x_Py_Px_Px_C1001000003_b;
  abcd[151] = 2.0E0*I_ERI_F2xy_Py_Px_Px_C1001000003_b;
  abcd[152] = 2.0E0*I_ERI_F2xz_Py_Px_Px_C1001000003_b;
  abcd[153] = 2.0E0*I_ERI_Fx2y_Py_Px_Px_C1001000003_b;
  abcd[154] = 2.0E0*I_ERI_Fxyz_Py_Px_Px_C1001000003_b;
  abcd[155] = 2.0E0*I_ERI_Fx2z_Py_Px_Px_C1001000003_b;
  abcd[156] = 2.0E0*I_ERI_F3y_Py_Px_Px_C1001000003_b;
  abcd[157] = 2.0E0*I_ERI_F2yz_Py_Px_Px_C1001000003_b;
  abcd[158] = 2.0E0*I_ERI_Fy2z_Py_Px_Px_C1001000003_b;
  abcd[159] = 2.0E0*I_ERI_F3z_Py_Px_Px_C1001000003_b;
  abcd[160] = 2.0E0*I_ERI_F3x_Py_Py_Px_C1001000003_b;
  abcd[161] = 2.0E0*I_ERI_F2xy_Py_Py_Px_C1001000003_b;
  abcd[162] = 2.0E0*I_ERI_F2xz_Py_Py_Px_C1001000003_b;
  abcd[163] = 2.0E0*I_ERI_Fx2y_Py_Py_Px_C1001000003_b;
  abcd[164] = 2.0E0*I_ERI_Fxyz_Py_Py_Px_C1001000003_b;
  abcd[165] = 2.0E0*I_ERI_Fx2z_Py_Py_Px_C1001000003_b;
  abcd[166] = 2.0E0*I_ERI_F3y_Py_Py_Px_C1001000003_b;
  abcd[167] = 2.0E0*I_ERI_F2yz_Py_Py_Px_C1001000003_b;
  abcd[168] = 2.0E0*I_ERI_Fy2z_Py_Py_Px_C1001000003_b;
  abcd[169] = 2.0E0*I_ERI_F3z_Py_Py_Px_C1001000003_b;
  abcd[170] = 2.0E0*I_ERI_F3x_Py_Pz_Px_C1001000003_b;
  abcd[171] = 2.0E0*I_ERI_F2xy_Py_Pz_Px_C1001000003_b;
  abcd[172] = 2.0E0*I_ERI_F2xz_Py_Pz_Px_C1001000003_b;
  abcd[173] = 2.0E0*I_ERI_Fx2y_Py_Pz_Px_C1001000003_b;
  abcd[174] = 2.0E0*I_ERI_Fxyz_Py_Pz_Px_C1001000003_b;
  abcd[175] = 2.0E0*I_ERI_Fx2z_Py_Pz_Px_C1001000003_b;
  abcd[176] = 2.0E0*I_ERI_F3y_Py_Pz_Px_C1001000003_b;
  abcd[177] = 2.0E0*I_ERI_F2yz_Py_Pz_Px_C1001000003_b;
  abcd[178] = 2.0E0*I_ERI_Fy2z_Py_Pz_Px_C1001000003_b;
  abcd[179] = 2.0E0*I_ERI_F3z_Py_Pz_Px_C1001000003_b;
  abcd[180] = 2.0E0*I_ERI_F3x_Py_Px_Py_C1001000003_b;
  abcd[181] = 2.0E0*I_ERI_F2xy_Py_Px_Py_C1001000003_b;
  abcd[182] = 2.0E0*I_ERI_F2xz_Py_Px_Py_C1001000003_b;
  abcd[183] = 2.0E0*I_ERI_Fx2y_Py_Px_Py_C1001000003_b;
  abcd[184] = 2.0E0*I_ERI_Fxyz_Py_Px_Py_C1001000003_b;
  abcd[185] = 2.0E0*I_ERI_Fx2z_Py_Px_Py_C1001000003_b;
  abcd[186] = 2.0E0*I_ERI_F3y_Py_Px_Py_C1001000003_b;
  abcd[187] = 2.0E0*I_ERI_F2yz_Py_Px_Py_C1001000003_b;
  abcd[188] = 2.0E0*I_ERI_Fy2z_Py_Px_Py_C1001000003_b;
  abcd[189] = 2.0E0*I_ERI_F3z_Py_Px_Py_C1001000003_b;
  abcd[190] = 2.0E0*I_ERI_F3x_Py_Py_Py_C1001000003_b;
  abcd[191] = 2.0E0*I_ERI_F2xy_Py_Py_Py_C1001000003_b;
  abcd[192] = 2.0E0*I_ERI_F2xz_Py_Py_Py_C1001000003_b;
  abcd[193] = 2.0E0*I_ERI_Fx2y_Py_Py_Py_C1001000003_b;
  abcd[194] = 2.0E0*I_ERI_Fxyz_Py_Py_Py_C1001000003_b;
  abcd[195] = 2.0E0*I_ERI_Fx2z_Py_Py_Py_C1001000003_b;
  abcd[196] = 2.0E0*I_ERI_F3y_Py_Py_Py_C1001000003_b;
  abcd[197] = 2.0E0*I_ERI_F2yz_Py_Py_Py_C1001000003_b;
  abcd[198] = 2.0E0*I_ERI_Fy2z_Py_Py_Py_C1001000003_b;
  abcd[199] = 2.0E0*I_ERI_F3z_Py_Py_Py_C1001000003_b;
  abcd[200] = 2.0E0*I_ERI_F3x_Py_Pz_Py_C1001000003_b;
  abcd[201] = 2.0E0*I_ERI_F2xy_Py_Pz_Py_C1001000003_b;
  abcd[202] = 2.0E0*I_ERI_F2xz_Py_Pz_Py_C1001000003_b;
  abcd[203] = 2.0E0*I_ERI_Fx2y_Py_Pz_Py_C1001000003_b;
  abcd[204] = 2.0E0*I_ERI_Fxyz_Py_Pz_Py_C1001000003_b;
  abcd[205] = 2.0E0*I_ERI_Fx2z_Py_Pz_Py_C1001000003_b;
  abcd[206] = 2.0E0*I_ERI_F3y_Py_Pz_Py_C1001000003_b;
  abcd[207] = 2.0E0*I_ERI_F2yz_Py_Pz_Py_C1001000003_b;
  abcd[208] = 2.0E0*I_ERI_Fy2z_Py_Pz_Py_C1001000003_b;
  abcd[209] = 2.0E0*I_ERI_F3z_Py_Pz_Py_C1001000003_b;
  abcd[210] = 2.0E0*I_ERI_F3x_Py_Px_Pz_C1001000003_b;
  abcd[211] = 2.0E0*I_ERI_F2xy_Py_Px_Pz_C1001000003_b;
  abcd[212] = 2.0E0*I_ERI_F2xz_Py_Px_Pz_C1001000003_b;
  abcd[213] = 2.0E0*I_ERI_Fx2y_Py_Px_Pz_C1001000003_b;
  abcd[214] = 2.0E0*I_ERI_Fxyz_Py_Px_Pz_C1001000003_b;
  abcd[215] = 2.0E0*I_ERI_Fx2z_Py_Px_Pz_C1001000003_b;
  abcd[216] = 2.0E0*I_ERI_F3y_Py_Px_Pz_C1001000003_b;
  abcd[217] = 2.0E0*I_ERI_F2yz_Py_Px_Pz_C1001000003_b;
  abcd[218] = 2.0E0*I_ERI_Fy2z_Py_Px_Pz_C1001000003_b;
  abcd[219] = 2.0E0*I_ERI_F3z_Py_Px_Pz_C1001000003_b;
  abcd[220] = 2.0E0*I_ERI_F3x_Py_Py_Pz_C1001000003_b;
  abcd[221] = 2.0E0*I_ERI_F2xy_Py_Py_Pz_C1001000003_b;
  abcd[222] = 2.0E0*I_ERI_F2xz_Py_Py_Pz_C1001000003_b;
  abcd[223] = 2.0E0*I_ERI_Fx2y_Py_Py_Pz_C1001000003_b;
  abcd[224] = 2.0E0*I_ERI_Fxyz_Py_Py_Pz_C1001000003_b;
  abcd[225] = 2.0E0*I_ERI_Fx2z_Py_Py_Pz_C1001000003_b;
  abcd[226] = 2.0E0*I_ERI_F3y_Py_Py_Pz_C1001000003_b;
  abcd[227] = 2.0E0*I_ERI_F2yz_Py_Py_Pz_C1001000003_b;
  abcd[228] = 2.0E0*I_ERI_Fy2z_Py_Py_Pz_C1001000003_b;
  abcd[229] = 2.0E0*I_ERI_F3z_Py_Py_Pz_C1001000003_b;
  abcd[230] = 2.0E0*I_ERI_F3x_Py_Pz_Pz_C1001000003_b;
  abcd[231] = 2.0E0*I_ERI_F2xy_Py_Pz_Pz_C1001000003_b;
  abcd[232] = 2.0E0*I_ERI_F2xz_Py_Pz_Pz_C1001000003_b;
  abcd[233] = 2.0E0*I_ERI_Fx2y_Py_Pz_Pz_C1001000003_b;
  abcd[234] = 2.0E0*I_ERI_Fxyz_Py_Pz_Pz_C1001000003_b;
  abcd[235] = 2.0E0*I_ERI_Fx2z_Py_Pz_Pz_C1001000003_b;
  abcd[236] = 2.0E0*I_ERI_F3y_Py_Pz_Pz_C1001000003_b;
  abcd[237] = 2.0E0*I_ERI_F2yz_Py_Pz_Pz_C1001000003_b;
  abcd[238] = 2.0E0*I_ERI_Fy2z_Py_Pz_Pz_C1001000003_b;
  abcd[239] = 2.0E0*I_ERI_F3z_Py_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[240] = 2.0E0*I_ERI_F3x_Pz_Px_S_C1000003_b;
  abcd[241] = 2.0E0*I_ERI_F2xy_Pz_Px_S_C1000003_b;
  abcd[242] = 2.0E0*I_ERI_F2xz_Pz_Px_S_C1000003_b;
  abcd[243] = 2.0E0*I_ERI_Fx2y_Pz_Px_S_C1000003_b;
  abcd[244] = 2.0E0*I_ERI_Fxyz_Pz_Px_S_C1000003_b;
  abcd[245] = 2.0E0*I_ERI_Fx2z_Pz_Px_S_C1000003_b;
  abcd[246] = 2.0E0*I_ERI_F3y_Pz_Px_S_C1000003_b;
  abcd[247] = 2.0E0*I_ERI_F2yz_Pz_Px_S_C1000003_b;
  abcd[248] = 2.0E0*I_ERI_Fy2z_Pz_Px_S_C1000003_b;
  abcd[249] = 2.0E0*I_ERI_F3z_Pz_Px_S_C1000003_b;
  abcd[250] = 2.0E0*I_ERI_F3x_Pz_Py_S_C1000003_b;
  abcd[251] = 2.0E0*I_ERI_F2xy_Pz_Py_S_C1000003_b;
  abcd[252] = 2.0E0*I_ERI_F2xz_Pz_Py_S_C1000003_b;
  abcd[253] = 2.0E0*I_ERI_Fx2y_Pz_Py_S_C1000003_b;
  abcd[254] = 2.0E0*I_ERI_Fxyz_Pz_Py_S_C1000003_b;
  abcd[255] = 2.0E0*I_ERI_Fx2z_Pz_Py_S_C1000003_b;
  abcd[256] = 2.0E0*I_ERI_F3y_Pz_Py_S_C1000003_b;
  abcd[257] = 2.0E0*I_ERI_F2yz_Pz_Py_S_C1000003_b;
  abcd[258] = 2.0E0*I_ERI_Fy2z_Pz_Py_S_C1000003_b;
  abcd[259] = 2.0E0*I_ERI_F3z_Pz_Py_S_C1000003_b;
  abcd[260] = 2.0E0*I_ERI_F3x_Pz_Pz_S_C1000003_b;
  abcd[261] = 2.0E0*I_ERI_F2xy_Pz_Pz_S_C1000003_b;
  abcd[262] = 2.0E0*I_ERI_F2xz_Pz_Pz_S_C1000003_b;
  abcd[263] = 2.0E0*I_ERI_Fx2y_Pz_Pz_S_C1000003_b;
  abcd[264] = 2.0E0*I_ERI_Fxyz_Pz_Pz_S_C1000003_b;
  abcd[265] = 2.0E0*I_ERI_Fx2z_Pz_Pz_S_C1000003_b;
  abcd[266] = 2.0E0*I_ERI_F3y_Pz_Pz_S_C1000003_b;
  abcd[267] = 2.0E0*I_ERI_F2yz_Pz_Pz_S_C1000003_b;
  abcd[268] = 2.0E0*I_ERI_Fy2z_Pz_Pz_S_C1000003_b;
  abcd[269] = 2.0E0*I_ERI_F3z_Pz_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[270] = 2.0E0*I_ERI_F3x_Pz_Px_Px_C1001000003_b;
  abcd[271] = 2.0E0*I_ERI_F2xy_Pz_Px_Px_C1001000003_b;
  abcd[272] = 2.0E0*I_ERI_F2xz_Pz_Px_Px_C1001000003_b;
  abcd[273] = 2.0E0*I_ERI_Fx2y_Pz_Px_Px_C1001000003_b;
  abcd[274] = 2.0E0*I_ERI_Fxyz_Pz_Px_Px_C1001000003_b;
  abcd[275] = 2.0E0*I_ERI_Fx2z_Pz_Px_Px_C1001000003_b;
  abcd[276] = 2.0E0*I_ERI_F3y_Pz_Px_Px_C1001000003_b;
  abcd[277] = 2.0E0*I_ERI_F2yz_Pz_Px_Px_C1001000003_b;
  abcd[278] = 2.0E0*I_ERI_Fy2z_Pz_Px_Px_C1001000003_b;
  abcd[279] = 2.0E0*I_ERI_F3z_Pz_Px_Px_C1001000003_b;
  abcd[280] = 2.0E0*I_ERI_F3x_Pz_Py_Px_C1001000003_b;
  abcd[281] = 2.0E0*I_ERI_F2xy_Pz_Py_Px_C1001000003_b;
  abcd[282] = 2.0E0*I_ERI_F2xz_Pz_Py_Px_C1001000003_b;
  abcd[283] = 2.0E0*I_ERI_Fx2y_Pz_Py_Px_C1001000003_b;
  abcd[284] = 2.0E0*I_ERI_Fxyz_Pz_Py_Px_C1001000003_b;
  abcd[285] = 2.0E0*I_ERI_Fx2z_Pz_Py_Px_C1001000003_b;
  abcd[286] = 2.0E0*I_ERI_F3y_Pz_Py_Px_C1001000003_b;
  abcd[287] = 2.0E0*I_ERI_F2yz_Pz_Py_Px_C1001000003_b;
  abcd[288] = 2.0E0*I_ERI_Fy2z_Pz_Py_Px_C1001000003_b;
  abcd[289] = 2.0E0*I_ERI_F3z_Pz_Py_Px_C1001000003_b;
  abcd[290] = 2.0E0*I_ERI_F3x_Pz_Pz_Px_C1001000003_b;
  abcd[291] = 2.0E0*I_ERI_F2xy_Pz_Pz_Px_C1001000003_b;
  abcd[292] = 2.0E0*I_ERI_F2xz_Pz_Pz_Px_C1001000003_b;
  abcd[293] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Px_C1001000003_b;
  abcd[294] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Px_C1001000003_b;
  abcd[295] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Px_C1001000003_b;
  abcd[296] = 2.0E0*I_ERI_F3y_Pz_Pz_Px_C1001000003_b;
  abcd[297] = 2.0E0*I_ERI_F2yz_Pz_Pz_Px_C1001000003_b;
  abcd[298] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Px_C1001000003_b;
  abcd[299] = 2.0E0*I_ERI_F3z_Pz_Pz_Px_C1001000003_b;
  abcd[300] = 2.0E0*I_ERI_F3x_Pz_Px_Py_C1001000003_b;
  abcd[301] = 2.0E0*I_ERI_F2xy_Pz_Px_Py_C1001000003_b;
  abcd[302] = 2.0E0*I_ERI_F2xz_Pz_Px_Py_C1001000003_b;
  abcd[303] = 2.0E0*I_ERI_Fx2y_Pz_Px_Py_C1001000003_b;
  abcd[304] = 2.0E0*I_ERI_Fxyz_Pz_Px_Py_C1001000003_b;
  abcd[305] = 2.0E0*I_ERI_Fx2z_Pz_Px_Py_C1001000003_b;
  abcd[306] = 2.0E0*I_ERI_F3y_Pz_Px_Py_C1001000003_b;
  abcd[307] = 2.0E0*I_ERI_F2yz_Pz_Px_Py_C1001000003_b;
  abcd[308] = 2.0E0*I_ERI_Fy2z_Pz_Px_Py_C1001000003_b;
  abcd[309] = 2.0E0*I_ERI_F3z_Pz_Px_Py_C1001000003_b;
  abcd[310] = 2.0E0*I_ERI_F3x_Pz_Py_Py_C1001000003_b;
  abcd[311] = 2.0E0*I_ERI_F2xy_Pz_Py_Py_C1001000003_b;
  abcd[312] = 2.0E0*I_ERI_F2xz_Pz_Py_Py_C1001000003_b;
  abcd[313] = 2.0E0*I_ERI_Fx2y_Pz_Py_Py_C1001000003_b;
  abcd[314] = 2.0E0*I_ERI_Fxyz_Pz_Py_Py_C1001000003_b;
  abcd[315] = 2.0E0*I_ERI_Fx2z_Pz_Py_Py_C1001000003_b;
  abcd[316] = 2.0E0*I_ERI_F3y_Pz_Py_Py_C1001000003_b;
  abcd[317] = 2.0E0*I_ERI_F2yz_Pz_Py_Py_C1001000003_b;
  abcd[318] = 2.0E0*I_ERI_Fy2z_Pz_Py_Py_C1001000003_b;
  abcd[319] = 2.0E0*I_ERI_F3z_Pz_Py_Py_C1001000003_b;
  abcd[320] = 2.0E0*I_ERI_F3x_Pz_Pz_Py_C1001000003_b;
  abcd[321] = 2.0E0*I_ERI_F2xy_Pz_Pz_Py_C1001000003_b;
  abcd[322] = 2.0E0*I_ERI_F2xz_Pz_Pz_Py_C1001000003_b;
  abcd[323] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Py_C1001000003_b;
  abcd[324] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Py_C1001000003_b;
  abcd[325] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Py_C1001000003_b;
  abcd[326] = 2.0E0*I_ERI_F3y_Pz_Pz_Py_C1001000003_b;
  abcd[327] = 2.0E0*I_ERI_F2yz_Pz_Pz_Py_C1001000003_b;
  abcd[328] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Py_C1001000003_b;
  abcd[329] = 2.0E0*I_ERI_F3z_Pz_Pz_Py_C1001000003_b;
  abcd[330] = 2.0E0*I_ERI_F3x_Pz_Px_Pz_C1001000003_b;
  abcd[331] = 2.0E0*I_ERI_F2xy_Pz_Px_Pz_C1001000003_b;
  abcd[332] = 2.0E0*I_ERI_F2xz_Pz_Px_Pz_C1001000003_b;
  abcd[333] = 2.0E0*I_ERI_Fx2y_Pz_Px_Pz_C1001000003_b;
  abcd[334] = 2.0E0*I_ERI_Fxyz_Pz_Px_Pz_C1001000003_b;
  abcd[335] = 2.0E0*I_ERI_Fx2z_Pz_Px_Pz_C1001000003_b;
  abcd[336] = 2.0E0*I_ERI_F3y_Pz_Px_Pz_C1001000003_b;
  abcd[337] = 2.0E0*I_ERI_F2yz_Pz_Px_Pz_C1001000003_b;
  abcd[338] = 2.0E0*I_ERI_Fy2z_Pz_Px_Pz_C1001000003_b;
  abcd[339] = 2.0E0*I_ERI_F3z_Pz_Px_Pz_C1001000003_b;
  abcd[340] = 2.0E0*I_ERI_F3x_Pz_Py_Pz_C1001000003_b;
  abcd[341] = 2.0E0*I_ERI_F2xy_Pz_Py_Pz_C1001000003_b;
  abcd[342] = 2.0E0*I_ERI_F2xz_Pz_Py_Pz_C1001000003_b;
  abcd[343] = 2.0E0*I_ERI_Fx2y_Pz_Py_Pz_C1001000003_b;
  abcd[344] = 2.0E0*I_ERI_Fxyz_Pz_Py_Pz_C1001000003_b;
  abcd[345] = 2.0E0*I_ERI_Fx2z_Pz_Py_Pz_C1001000003_b;
  abcd[346] = 2.0E0*I_ERI_F3y_Pz_Py_Pz_C1001000003_b;
  abcd[347] = 2.0E0*I_ERI_F2yz_Pz_Py_Pz_C1001000003_b;
  abcd[348] = 2.0E0*I_ERI_Fy2z_Pz_Py_Pz_C1001000003_b;
  abcd[349] = 2.0E0*I_ERI_F3z_Pz_Py_Pz_C1001000003_b;
  abcd[350] = 2.0E0*I_ERI_F3x_Pz_Pz_Pz_C1001000003_b;
  abcd[351] = 2.0E0*I_ERI_F2xy_Pz_Pz_Pz_C1001000003_b;
  abcd[352] = 2.0E0*I_ERI_F2xz_Pz_Pz_Pz_C1001000003_b;
  abcd[353] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Pz_C1001000003_b;
  abcd[354] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Pz_C1001000003_b;
  abcd[355] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Pz_C1001000003_b;
  abcd[356] = 2.0E0*I_ERI_F3y_Pz_Pz_Pz_C1001000003_b;
  abcd[357] = 2.0E0*I_ERI_F2yz_Pz_Pz_Pz_C1001000003_b;
  abcd[358] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Pz_C1001000003_b;
  abcd[359] = 2.0E0*I_ERI_F3z_Pz_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[360] = 2.0E0*I_ERI_F3x_S_D2x_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[361] = 2.0E0*I_ERI_F2xy_S_D2x_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[362] = 2.0E0*I_ERI_F2xz_S_D2x_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[363] = 2.0E0*I_ERI_Fx2y_S_D2x_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[364] = 2.0E0*I_ERI_Fxyz_S_D2x_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[365] = 2.0E0*I_ERI_Fx2z_S_D2x_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[366] = 2.0E0*I_ERI_F3y_S_D2x_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[367] = 2.0E0*I_ERI_F2yz_S_D2x_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[368] = 2.0E0*I_ERI_Fy2z_S_D2x_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[369] = 2.0E0*I_ERI_F3z_S_D2x_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;
  abcd[370] = 2.0E0*I_ERI_F3x_S_Dxy_S_C1000003_c;
  abcd[371] = 2.0E0*I_ERI_F2xy_S_Dxy_S_C1000003_c;
  abcd[372] = 2.0E0*I_ERI_F2xz_S_Dxy_S_C1000003_c;
  abcd[373] = 2.0E0*I_ERI_Fx2y_S_Dxy_S_C1000003_c;
  abcd[374] = 2.0E0*I_ERI_Fxyz_S_Dxy_S_C1000003_c;
  abcd[375] = 2.0E0*I_ERI_Fx2z_S_Dxy_S_C1000003_c;
  abcd[376] = 2.0E0*I_ERI_F3y_S_Dxy_S_C1000003_c;
  abcd[377] = 2.0E0*I_ERI_F2yz_S_Dxy_S_C1000003_c;
  abcd[378] = 2.0E0*I_ERI_Fy2z_S_Dxy_S_C1000003_c;
  abcd[379] = 2.0E0*I_ERI_F3z_S_Dxy_S_C1000003_c;
  abcd[380] = 2.0E0*I_ERI_F3x_S_Dxz_S_C1000003_c;
  abcd[381] = 2.0E0*I_ERI_F2xy_S_Dxz_S_C1000003_c;
  abcd[382] = 2.0E0*I_ERI_F2xz_S_Dxz_S_C1000003_c;
  abcd[383] = 2.0E0*I_ERI_Fx2y_S_Dxz_S_C1000003_c;
  abcd[384] = 2.0E0*I_ERI_Fxyz_S_Dxz_S_C1000003_c;
  abcd[385] = 2.0E0*I_ERI_Fx2z_S_Dxz_S_C1000003_c;
  abcd[386] = 2.0E0*I_ERI_F3y_S_Dxz_S_C1000003_c;
  abcd[387] = 2.0E0*I_ERI_F2yz_S_Dxz_S_C1000003_c;
  abcd[388] = 2.0E0*I_ERI_Fy2z_S_Dxz_S_C1000003_c;
  abcd[389] = 2.0E0*I_ERI_F3z_S_Dxz_S_C1000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[390] = 2.0E0*I_ERI_F3x_S_D2x_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[391] = 2.0E0*I_ERI_F2xy_S_D2x_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[392] = 2.0E0*I_ERI_F2xz_S_D2x_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[393] = 2.0E0*I_ERI_Fx2y_S_D2x_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[394] = 2.0E0*I_ERI_Fxyz_S_D2x_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[395] = 2.0E0*I_ERI_Fx2z_S_D2x_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[396] = 2.0E0*I_ERI_F3y_S_D2x_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[397] = 2.0E0*I_ERI_F2yz_S_D2x_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[398] = 2.0E0*I_ERI_Fy2z_S_D2x_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[399] = 2.0E0*I_ERI_F3z_S_D2x_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[400] = 2.0E0*I_ERI_F3x_S_Dxy_Px_C1001000003_c;
  abcd[401] = 2.0E0*I_ERI_F2xy_S_Dxy_Px_C1001000003_c;
  abcd[402] = 2.0E0*I_ERI_F2xz_S_Dxy_Px_C1001000003_c;
  abcd[403] = 2.0E0*I_ERI_Fx2y_S_Dxy_Px_C1001000003_c;
  abcd[404] = 2.0E0*I_ERI_Fxyz_S_Dxy_Px_C1001000003_c;
  abcd[405] = 2.0E0*I_ERI_Fx2z_S_Dxy_Px_C1001000003_c;
  abcd[406] = 2.0E0*I_ERI_F3y_S_Dxy_Px_C1001000003_c;
  abcd[407] = 2.0E0*I_ERI_F2yz_S_Dxy_Px_C1001000003_c;
  abcd[408] = 2.0E0*I_ERI_Fy2z_S_Dxy_Px_C1001000003_c;
  abcd[409] = 2.0E0*I_ERI_F3z_S_Dxy_Px_C1001000003_c;
  abcd[410] = 2.0E0*I_ERI_F3x_S_Dxz_Px_C1001000003_c;
  abcd[411] = 2.0E0*I_ERI_F2xy_S_Dxz_Px_C1001000003_c;
  abcd[412] = 2.0E0*I_ERI_F2xz_S_Dxz_Px_C1001000003_c;
  abcd[413] = 2.0E0*I_ERI_Fx2y_S_Dxz_Px_C1001000003_c;
  abcd[414] = 2.0E0*I_ERI_Fxyz_S_Dxz_Px_C1001000003_c;
  abcd[415] = 2.0E0*I_ERI_Fx2z_S_Dxz_Px_C1001000003_c;
  abcd[416] = 2.0E0*I_ERI_F3y_S_Dxz_Px_C1001000003_c;
  abcd[417] = 2.0E0*I_ERI_F2yz_S_Dxz_Px_C1001000003_c;
  abcd[418] = 2.0E0*I_ERI_Fy2z_S_Dxz_Px_C1001000003_c;
  abcd[419] = 2.0E0*I_ERI_F3z_S_Dxz_Px_C1001000003_c;
  abcd[420] = 2.0E0*I_ERI_F3x_S_D2x_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[421] = 2.0E0*I_ERI_F2xy_S_D2x_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[422] = 2.0E0*I_ERI_F2xz_S_D2x_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[423] = 2.0E0*I_ERI_Fx2y_S_D2x_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[424] = 2.0E0*I_ERI_Fxyz_S_D2x_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[425] = 2.0E0*I_ERI_Fx2z_S_D2x_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[426] = 2.0E0*I_ERI_F3y_S_D2x_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[427] = 2.0E0*I_ERI_F2yz_S_D2x_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[428] = 2.0E0*I_ERI_Fy2z_S_D2x_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[429] = 2.0E0*I_ERI_F3z_S_D2x_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[430] = 2.0E0*I_ERI_F3x_S_Dxy_Py_C1001000003_c;
  abcd[431] = 2.0E0*I_ERI_F2xy_S_Dxy_Py_C1001000003_c;
  abcd[432] = 2.0E0*I_ERI_F2xz_S_Dxy_Py_C1001000003_c;
  abcd[433] = 2.0E0*I_ERI_Fx2y_S_Dxy_Py_C1001000003_c;
  abcd[434] = 2.0E0*I_ERI_Fxyz_S_Dxy_Py_C1001000003_c;
  abcd[435] = 2.0E0*I_ERI_Fx2z_S_Dxy_Py_C1001000003_c;
  abcd[436] = 2.0E0*I_ERI_F3y_S_Dxy_Py_C1001000003_c;
  abcd[437] = 2.0E0*I_ERI_F2yz_S_Dxy_Py_C1001000003_c;
  abcd[438] = 2.0E0*I_ERI_Fy2z_S_Dxy_Py_C1001000003_c;
  abcd[439] = 2.0E0*I_ERI_F3z_S_Dxy_Py_C1001000003_c;
  abcd[440] = 2.0E0*I_ERI_F3x_S_Dxz_Py_C1001000003_c;
  abcd[441] = 2.0E0*I_ERI_F2xy_S_Dxz_Py_C1001000003_c;
  abcd[442] = 2.0E0*I_ERI_F2xz_S_Dxz_Py_C1001000003_c;
  abcd[443] = 2.0E0*I_ERI_Fx2y_S_Dxz_Py_C1001000003_c;
  abcd[444] = 2.0E0*I_ERI_Fxyz_S_Dxz_Py_C1001000003_c;
  abcd[445] = 2.0E0*I_ERI_Fx2z_S_Dxz_Py_C1001000003_c;
  abcd[446] = 2.0E0*I_ERI_F3y_S_Dxz_Py_C1001000003_c;
  abcd[447] = 2.0E0*I_ERI_F2yz_S_Dxz_Py_C1001000003_c;
  abcd[448] = 2.0E0*I_ERI_Fy2z_S_Dxz_Py_C1001000003_c;
  abcd[449] = 2.0E0*I_ERI_F3z_S_Dxz_Py_C1001000003_c;
  abcd[450] = 2.0E0*I_ERI_F3x_S_D2x_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[451] = 2.0E0*I_ERI_F2xy_S_D2x_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[452] = 2.0E0*I_ERI_F2xz_S_D2x_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[453] = 2.0E0*I_ERI_Fx2y_S_D2x_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[454] = 2.0E0*I_ERI_Fxyz_S_D2x_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[455] = 2.0E0*I_ERI_Fx2z_S_D2x_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[456] = 2.0E0*I_ERI_F3y_S_D2x_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[457] = 2.0E0*I_ERI_F2yz_S_D2x_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[458] = 2.0E0*I_ERI_Fy2z_S_D2x_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[459] = 2.0E0*I_ERI_F3z_S_D2x_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;
  abcd[460] = 2.0E0*I_ERI_F3x_S_Dxy_Pz_C1001000003_c;
  abcd[461] = 2.0E0*I_ERI_F2xy_S_Dxy_Pz_C1001000003_c;
  abcd[462] = 2.0E0*I_ERI_F2xz_S_Dxy_Pz_C1001000003_c;
  abcd[463] = 2.0E0*I_ERI_Fx2y_S_Dxy_Pz_C1001000003_c;
  abcd[464] = 2.0E0*I_ERI_Fxyz_S_Dxy_Pz_C1001000003_c;
  abcd[465] = 2.0E0*I_ERI_Fx2z_S_Dxy_Pz_C1001000003_c;
  abcd[466] = 2.0E0*I_ERI_F3y_S_Dxy_Pz_C1001000003_c;
  abcd[467] = 2.0E0*I_ERI_F2yz_S_Dxy_Pz_C1001000003_c;
  abcd[468] = 2.0E0*I_ERI_Fy2z_S_Dxy_Pz_C1001000003_c;
  abcd[469] = 2.0E0*I_ERI_F3z_S_Dxy_Pz_C1001000003_c;
  abcd[470] = 2.0E0*I_ERI_F3x_S_Dxz_Pz_C1001000003_c;
  abcd[471] = 2.0E0*I_ERI_F2xy_S_Dxz_Pz_C1001000003_c;
  abcd[472] = 2.0E0*I_ERI_F2xz_S_Dxz_Pz_C1001000003_c;
  abcd[473] = 2.0E0*I_ERI_Fx2y_S_Dxz_Pz_C1001000003_c;
  abcd[474] = 2.0E0*I_ERI_Fxyz_S_Dxz_Pz_C1001000003_c;
  abcd[475] = 2.0E0*I_ERI_Fx2z_S_Dxz_Pz_C1001000003_c;
  abcd[476] = 2.0E0*I_ERI_F3y_S_Dxz_Pz_C1001000003_c;
  abcd[477] = 2.0E0*I_ERI_F2yz_S_Dxz_Pz_C1001000003_c;
  abcd[478] = 2.0E0*I_ERI_Fy2z_S_Dxz_Pz_C1001000003_c;
  abcd[479] = 2.0E0*I_ERI_F3z_S_Dxz_Pz_C1001000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[480] = 2.0E0*I_ERI_F3x_S_Dxy_S_C1000003_c;
  abcd[481] = 2.0E0*I_ERI_F2xy_S_Dxy_S_C1000003_c;
  abcd[482] = 2.0E0*I_ERI_F2xz_S_Dxy_S_C1000003_c;
  abcd[483] = 2.0E0*I_ERI_Fx2y_S_Dxy_S_C1000003_c;
  abcd[484] = 2.0E0*I_ERI_Fxyz_S_Dxy_S_C1000003_c;
  abcd[485] = 2.0E0*I_ERI_Fx2z_S_Dxy_S_C1000003_c;
  abcd[486] = 2.0E0*I_ERI_F3y_S_Dxy_S_C1000003_c;
  abcd[487] = 2.0E0*I_ERI_F2yz_S_Dxy_S_C1000003_c;
  abcd[488] = 2.0E0*I_ERI_Fy2z_S_Dxy_S_C1000003_c;
  abcd[489] = 2.0E0*I_ERI_F3z_S_Dxy_S_C1000003_c;
  abcd[490] = 2.0E0*I_ERI_F3x_S_D2y_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[491] = 2.0E0*I_ERI_F2xy_S_D2y_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[492] = 2.0E0*I_ERI_F2xz_S_D2y_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[493] = 2.0E0*I_ERI_Fx2y_S_D2y_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[494] = 2.0E0*I_ERI_Fxyz_S_D2y_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[495] = 2.0E0*I_ERI_Fx2z_S_D2y_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[496] = 2.0E0*I_ERI_F3y_S_D2y_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[497] = 2.0E0*I_ERI_F2yz_S_D2y_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[498] = 2.0E0*I_ERI_Fy2z_S_D2y_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[499] = 2.0E0*I_ERI_F3z_S_D2y_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;
  abcd[500] = 2.0E0*I_ERI_F3x_S_Dyz_S_C1000003_c;
  abcd[501] = 2.0E0*I_ERI_F2xy_S_Dyz_S_C1000003_c;
  abcd[502] = 2.0E0*I_ERI_F2xz_S_Dyz_S_C1000003_c;
  abcd[503] = 2.0E0*I_ERI_Fx2y_S_Dyz_S_C1000003_c;
  abcd[504] = 2.0E0*I_ERI_Fxyz_S_Dyz_S_C1000003_c;
  abcd[505] = 2.0E0*I_ERI_Fx2z_S_Dyz_S_C1000003_c;
  abcd[506] = 2.0E0*I_ERI_F3y_S_Dyz_S_C1000003_c;
  abcd[507] = 2.0E0*I_ERI_F2yz_S_Dyz_S_C1000003_c;
  abcd[508] = 2.0E0*I_ERI_Fy2z_S_Dyz_S_C1000003_c;
  abcd[509] = 2.0E0*I_ERI_F3z_S_Dyz_S_C1000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[510] = 2.0E0*I_ERI_F3x_S_Dxy_Px_C1001000003_c;
  abcd[511] = 2.0E0*I_ERI_F2xy_S_Dxy_Px_C1001000003_c;
  abcd[512] = 2.0E0*I_ERI_F2xz_S_Dxy_Px_C1001000003_c;
  abcd[513] = 2.0E0*I_ERI_Fx2y_S_Dxy_Px_C1001000003_c;
  abcd[514] = 2.0E0*I_ERI_Fxyz_S_Dxy_Px_C1001000003_c;
  abcd[515] = 2.0E0*I_ERI_Fx2z_S_Dxy_Px_C1001000003_c;
  abcd[516] = 2.0E0*I_ERI_F3y_S_Dxy_Px_C1001000003_c;
  abcd[517] = 2.0E0*I_ERI_F2yz_S_Dxy_Px_C1001000003_c;
  abcd[518] = 2.0E0*I_ERI_Fy2z_S_Dxy_Px_C1001000003_c;
  abcd[519] = 2.0E0*I_ERI_F3z_S_Dxy_Px_C1001000003_c;
  abcd[520] = 2.0E0*I_ERI_F3x_S_D2y_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[521] = 2.0E0*I_ERI_F2xy_S_D2y_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[522] = 2.0E0*I_ERI_F2xz_S_D2y_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[523] = 2.0E0*I_ERI_Fx2y_S_D2y_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[524] = 2.0E0*I_ERI_Fxyz_S_D2y_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[525] = 2.0E0*I_ERI_Fx2z_S_D2y_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[526] = 2.0E0*I_ERI_F3y_S_D2y_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[527] = 2.0E0*I_ERI_F2yz_S_D2y_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[528] = 2.0E0*I_ERI_Fy2z_S_D2y_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[529] = 2.0E0*I_ERI_F3z_S_D2y_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[530] = 2.0E0*I_ERI_F3x_S_Dyz_Px_C1001000003_c;
  abcd[531] = 2.0E0*I_ERI_F2xy_S_Dyz_Px_C1001000003_c;
  abcd[532] = 2.0E0*I_ERI_F2xz_S_Dyz_Px_C1001000003_c;
  abcd[533] = 2.0E0*I_ERI_Fx2y_S_Dyz_Px_C1001000003_c;
  abcd[534] = 2.0E0*I_ERI_Fxyz_S_Dyz_Px_C1001000003_c;
  abcd[535] = 2.0E0*I_ERI_Fx2z_S_Dyz_Px_C1001000003_c;
  abcd[536] = 2.0E0*I_ERI_F3y_S_Dyz_Px_C1001000003_c;
  abcd[537] = 2.0E0*I_ERI_F2yz_S_Dyz_Px_C1001000003_c;
  abcd[538] = 2.0E0*I_ERI_Fy2z_S_Dyz_Px_C1001000003_c;
  abcd[539] = 2.0E0*I_ERI_F3z_S_Dyz_Px_C1001000003_c;
  abcd[540] = 2.0E0*I_ERI_F3x_S_Dxy_Py_C1001000003_c;
  abcd[541] = 2.0E0*I_ERI_F2xy_S_Dxy_Py_C1001000003_c;
  abcd[542] = 2.0E0*I_ERI_F2xz_S_Dxy_Py_C1001000003_c;
  abcd[543] = 2.0E0*I_ERI_Fx2y_S_Dxy_Py_C1001000003_c;
  abcd[544] = 2.0E0*I_ERI_Fxyz_S_Dxy_Py_C1001000003_c;
  abcd[545] = 2.0E0*I_ERI_Fx2z_S_Dxy_Py_C1001000003_c;
  abcd[546] = 2.0E0*I_ERI_F3y_S_Dxy_Py_C1001000003_c;
  abcd[547] = 2.0E0*I_ERI_F2yz_S_Dxy_Py_C1001000003_c;
  abcd[548] = 2.0E0*I_ERI_Fy2z_S_Dxy_Py_C1001000003_c;
  abcd[549] = 2.0E0*I_ERI_F3z_S_Dxy_Py_C1001000003_c;
  abcd[550] = 2.0E0*I_ERI_F3x_S_D2y_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[551] = 2.0E0*I_ERI_F2xy_S_D2y_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[552] = 2.0E0*I_ERI_F2xz_S_D2y_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[553] = 2.0E0*I_ERI_Fx2y_S_D2y_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[554] = 2.0E0*I_ERI_Fxyz_S_D2y_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[555] = 2.0E0*I_ERI_Fx2z_S_D2y_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[556] = 2.0E0*I_ERI_F3y_S_D2y_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[557] = 2.0E0*I_ERI_F2yz_S_D2y_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[558] = 2.0E0*I_ERI_Fy2z_S_D2y_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[559] = 2.0E0*I_ERI_F3z_S_D2y_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[560] = 2.0E0*I_ERI_F3x_S_Dyz_Py_C1001000003_c;
  abcd[561] = 2.0E0*I_ERI_F2xy_S_Dyz_Py_C1001000003_c;
  abcd[562] = 2.0E0*I_ERI_F2xz_S_Dyz_Py_C1001000003_c;
  abcd[563] = 2.0E0*I_ERI_Fx2y_S_Dyz_Py_C1001000003_c;
  abcd[564] = 2.0E0*I_ERI_Fxyz_S_Dyz_Py_C1001000003_c;
  abcd[565] = 2.0E0*I_ERI_Fx2z_S_Dyz_Py_C1001000003_c;
  abcd[566] = 2.0E0*I_ERI_F3y_S_Dyz_Py_C1001000003_c;
  abcd[567] = 2.0E0*I_ERI_F2yz_S_Dyz_Py_C1001000003_c;
  abcd[568] = 2.0E0*I_ERI_Fy2z_S_Dyz_Py_C1001000003_c;
  abcd[569] = 2.0E0*I_ERI_F3z_S_Dyz_Py_C1001000003_c;
  abcd[570] = 2.0E0*I_ERI_F3x_S_Dxy_Pz_C1001000003_c;
  abcd[571] = 2.0E0*I_ERI_F2xy_S_Dxy_Pz_C1001000003_c;
  abcd[572] = 2.0E0*I_ERI_F2xz_S_Dxy_Pz_C1001000003_c;
  abcd[573] = 2.0E0*I_ERI_Fx2y_S_Dxy_Pz_C1001000003_c;
  abcd[574] = 2.0E0*I_ERI_Fxyz_S_Dxy_Pz_C1001000003_c;
  abcd[575] = 2.0E0*I_ERI_Fx2z_S_Dxy_Pz_C1001000003_c;
  abcd[576] = 2.0E0*I_ERI_F3y_S_Dxy_Pz_C1001000003_c;
  abcd[577] = 2.0E0*I_ERI_F2yz_S_Dxy_Pz_C1001000003_c;
  abcd[578] = 2.0E0*I_ERI_Fy2z_S_Dxy_Pz_C1001000003_c;
  abcd[579] = 2.0E0*I_ERI_F3z_S_Dxy_Pz_C1001000003_c;
  abcd[580] = 2.0E0*I_ERI_F3x_S_D2y_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[581] = 2.0E0*I_ERI_F2xy_S_D2y_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[582] = 2.0E0*I_ERI_F2xz_S_D2y_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[583] = 2.0E0*I_ERI_Fx2y_S_D2y_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[584] = 2.0E0*I_ERI_Fxyz_S_D2y_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[585] = 2.0E0*I_ERI_Fx2z_S_D2y_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[586] = 2.0E0*I_ERI_F3y_S_D2y_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[587] = 2.0E0*I_ERI_F2yz_S_D2y_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[588] = 2.0E0*I_ERI_Fy2z_S_D2y_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[589] = 2.0E0*I_ERI_F3z_S_D2y_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;
  abcd[590] = 2.0E0*I_ERI_F3x_S_Dyz_Pz_C1001000003_c;
  abcd[591] = 2.0E0*I_ERI_F2xy_S_Dyz_Pz_C1001000003_c;
  abcd[592] = 2.0E0*I_ERI_F2xz_S_Dyz_Pz_C1001000003_c;
  abcd[593] = 2.0E0*I_ERI_Fx2y_S_Dyz_Pz_C1001000003_c;
  abcd[594] = 2.0E0*I_ERI_Fxyz_S_Dyz_Pz_C1001000003_c;
  abcd[595] = 2.0E0*I_ERI_Fx2z_S_Dyz_Pz_C1001000003_c;
  abcd[596] = 2.0E0*I_ERI_F3y_S_Dyz_Pz_C1001000003_c;
  abcd[597] = 2.0E0*I_ERI_F2yz_S_Dyz_Pz_C1001000003_c;
  abcd[598] = 2.0E0*I_ERI_Fy2z_S_Dyz_Pz_C1001000003_c;
  abcd[599] = 2.0E0*I_ERI_F3z_S_Dyz_Pz_C1001000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[600] = 2.0E0*I_ERI_F3x_S_Dxz_S_C1000003_c;
  abcd[601] = 2.0E0*I_ERI_F2xy_S_Dxz_S_C1000003_c;
  abcd[602] = 2.0E0*I_ERI_F2xz_S_Dxz_S_C1000003_c;
  abcd[603] = 2.0E0*I_ERI_Fx2y_S_Dxz_S_C1000003_c;
  abcd[604] = 2.0E0*I_ERI_Fxyz_S_Dxz_S_C1000003_c;
  abcd[605] = 2.0E0*I_ERI_Fx2z_S_Dxz_S_C1000003_c;
  abcd[606] = 2.0E0*I_ERI_F3y_S_Dxz_S_C1000003_c;
  abcd[607] = 2.0E0*I_ERI_F2yz_S_Dxz_S_C1000003_c;
  abcd[608] = 2.0E0*I_ERI_Fy2z_S_Dxz_S_C1000003_c;
  abcd[609] = 2.0E0*I_ERI_F3z_S_Dxz_S_C1000003_c;
  abcd[610] = 2.0E0*I_ERI_F3x_S_Dyz_S_C1000003_c;
  abcd[611] = 2.0E0*I_ERI_F2xy_S_Dyz_S_C1000003_c;
  abcd[612] = 2.0E0*I_ERI_F2xz_S_Dyz_S_C1000003_c;
  abcd[613] = 2.0E0*I_ERI_Fx2y_S_Dyz_S_C1000003_c;
  abcd[614] = 2.0E0*I_ERI_Fxyz_S_Dyz_S_C1000003_c;
  abcd[615] = 2.0E0*I_ERI_Fx2z_S_Dyz_S_C1000003_c;
  abcd[616] = 2.0E0*I_ERI_F3y_S_Dyz_S_C1000003_c;
  abcd[617] = 2.0E0*I_ERI_F2yz_S_Dyz_S_C1000003_c;
  abcd[618] = 2.0E0*I_ERI_Fy2z_S_Dyz_S_C1000003_c;
  abcd[619] = 2.0E0*I_ERI_F3z_S_Dyz_S_C1000003_c;
  abcd[620] = 2.0E0*I_ERI_F3x_S_D2z_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[621] = 2.0E0*I_ERI_F2xy_S_D2z_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[622] = 2.0E0*I_ERI_F2xz_S_D2z_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[623] = 2.0E0*I_ERI_Fx2y_S_D2z_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[624] = 2.0E0*I_ERI_Fxyz_S_D2z_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[625] = 2.0E0*I_ERI_Fx2z_S_D2z_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[626] = 2.0E0*I_ERI_F3y_S_D2z_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[627] = 2.0E0*I_ERI_F2yz_S_D2z_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[628] = 2.0E0*I_ERI_Fy2z_S_D2z_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[629] = 2.0E0*I_ERI_F3z_S_D2z_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[630] = 2.0E0*I_ERI_F3x_S_Dxz_Px_C1001000003_c;
  abcd[631] = 2.0E0*I_ERI_F2xy_S_Dxz_Px_C1001000003_c;
  abcd[632] = 2.0E0*I_ERI_F2xz_S_Dxz_Px_C1001000003_c;
  abcd[633] = 2.0E0*I_ERI_Fx2y_S_Dxz_Px_C1001000003_c;
  abcd[634] = 2.0E0*I_ERI_Fxyz_S_Dxz_Px_C1001000003_c;
  abcd[635] = 2.0E0*I_ERI_Fx2z_S_Dxz_Px_C1001000003_c;
  abcd[636] = 2.0E0*I_ERI_F3y_S_Dxz_Px_C1001000003_c;
  abcd[637] = 2.0E0*I_ERI_F2yz_S_Dxz_Px_C1001000003_c;
  abcd[638] = 2.0E0*I_ERI_Fy2z_S_Dxz_Px_C1001000003_c;
  abcd[639] = 2.0E0*I_ERI_F3z_S_Dxz_Px_C1001000003_c;
  abcd[640] = 2.0E0*I_ERI_F3x_S_Dyz_Px_C1001000003_c;
  abcd[641] = 2.0E0*I_ERI_F2xy_S_Dyz_Px_C1001000003_c;
  abcd[642] = 2.0E0*I_ERI_F2xz_S_Dyz_Px_C1001000003_c;
  abcd[643] = 2.0E0*I_ERI_Fx2y_S_Dyz_Px_C1001000003_c;
  abcd[644] = 2.0E0*I_ERI_Fxyz_S_Dyz_Px_C1001000003_c;
  abcd[645] = 2.0E0*I_ERI_Fx2z_S_Dyz_Px_C1001000003_c;
  abcd[646] = 2.0E0*I_ERI_F3y_S_Dyz_Px_C1001000003_c;
  abcd[647] = 2.0E0*I_ERI_F2yz_S_Dyz_Px_C1001000003_c;
  abcd[648] = 2.0E0*I_ERI_Fy2z_S_Dyz_Px_C1001000003_c;
  abcd[649] = 2.0E0*I_ERI_F3z_S_Dyz_Px_C1001000003_c;
  abcd[650] = 2.0E0*I_ERI_F3x_S_D2z_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[651] = 2.0E0*I_ERI_F2xy_S_D2z_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[652] = 2.0E0*I_ERI_F2xz_S_D2z_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[653] = 2.0E0*I_ERI_Fx2y_S_D2z_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[654] = 2.0E0*I_ERI_Fxyz_S_D2z_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[655] = 2.0E0*I_ERI_Fx2z_S_D2z_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[656] = 2.0E0*I_ERI_F3y_S_D2z_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[657] = 2.0E0*I_ERI_F2yz_S_D2z_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[658] = 2.0E0*I_ERI_Fy2z_S_D2z_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[659] = 2.0E0*I_ERI_F3z_S_D2z_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[660] = 2.0E0*I_ERI_F3x_S_Dxz_Py_C1001000003_c;
  abcd[661] = 2.0E0*I_ERI_F2xy_S_Dxz_Py_C1001000003_c;
  abcd[662] = 2.0E0*I_ERI_F2xz_S_Dxz_Py_C1001000003_c;
  abcd[663] = 2.0E0*I_ERI_Fx2y_S_Dxz_Py_C1001000003_c;
  abcd[664] = 2.0E0*I_ERI_Fxyz_S_Dxz_Py_C1001000003_c;
  abcd[665] = 2.0E0*I_ERI_Fx2z_S_Dxz_Py_C1001000003_c;
  abcd[666] = 2.0E0*I_ERI_F3y_S_Dxz_Py_C1001000003_c;
  abcd[667] = 2.0E0*I_ERI_F2yz_S_Dxz_Py_C1001000003_c;
  abcd[668] = 2.0E0*I_ERI_Fy2z_S_Dxz_Py_C1001000003_c;
  abcd[669] = 2.0E0*I_ERI_F3z_S_Dxz_Py_C1001000003_c;
  abcd[670] = 2.0E0*I_ERI_F3x_S_Dyz_Py_C1001000003_c;
  abcd[671] = 2.0E0*I_ERI_F2xy_S_Dyz_Py_C1001000003_c;
  abcd[672] = 2.0E0*I_ERI_F2xz_S_Dyz_Py_C1001000003_c;
  abcd[673] = 2.0E0*I_ERI_Fx2y_S_Dyz_Py_C1001000003_c;
  abcd[674] = 2.0E0*I_ERI_Fxyz_S_Dyz_Py_C1001000003_c;
  abcd[675] = 2.0E0*I_ERI_Fx2z_S_Dyz_Py_C1001000003_c;
  abcd[676] = 2.0E0*I_ERI_F3y_S_Dyz_Py_C1001000003_c;
  abcd[677] = 2.0E0*I_ERI_F2yz_S_Dyz_Py_C1001000003_c;
  abcd[678] = 2.0E0*I_ERI_Fy2z_S_Dyz_Py_C1001000003_c;
  abcd[679] = 2.0E0*I_ERI_F3z_S_Dyz_Py_C1001000003_c;
  abcd[680] = 2.0E0*I_ERI_F3x_S_D2z_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[681] = 2.0E0*I_ERI_F2xy_S_D2z_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[682] = 2.0E0*I_ERI_F2xz_S_D2z_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[683] = 2.0E0*I_ERI_Fx2y_S_D2z_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[684] = 2.0E0*I_ERI_Fxyz_S_D2z_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[685] = 2.0E0*I_ERI_Fx2z_S_D2z_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[686] = 2.0E0*I_ERI_F3y_S_D2z_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[687] = 2.0E0*I_ERI_F2yz_S_D2z_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[688] = 2.0E0*I_ERI_Fy2z_S_D2z_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[689] = 2.0E0*I_ERI_F3z_S_D2z_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[690] = 2.0E0*I_ERI_F3x_S_Dxz_Pz_C1001000003_c;
  abcd[691] = 2.0E0*I_ERI_F2xy_S_Dxz_Pz_C1001000003_c;
  abcd[692] = 2.0E0*I_ERI_F2xz_S_Dxz_Pz_C1001000003_c;
  abcd[693] = 2.0E0*I_ERI_Fx2y_S_Dxz_Pz_C1001000003_c;
  abcd[694] = 2.0E0*I_ERI_Fxyz_S_Dxz_Pz_C1001000003_c;
  abcd[695] = 2.0E0*I_ERI_Fx2z_S_Dxz_Pz_C1001000003_c;
  abcd[696] = 2.0E0*I_ERI_F3y_S_Dxz_Pz_C1001000003_c;
  abcd[697] = 2.0E0*I_ERI_F2yz_S_Dxz_Pz_C1001000003_c;
  abcd[698] = 2.0E0*I_ERI_Fy2z_S_Dxz_Pz_C1001000003_c;
  abcd[699] = 2.0E0*I_ERI_F3z_S_Dxz_Pz_C1001000003_c;
  abcd[700] = 2.0E0*I_ERI_F3x_S_Dyz_Pz_C1001000003_c;
  abcd[701] = 2.0E0*I_ERI_F2xy_S_Dyz_Pz_C1001000003_c;
  abcd[702] = 2.0E0*I_ERI_F2xz_S_Dyz_Pz_C1001000003_c;
  abcd[703] = 2.0E0*I_ERI_Fx2y_S_Dyz_Pz_C1001000003_c;
  abcd[704] = 2.0E0*I_ERI_Fxyz_S_Dyz_Pz_C1001000003_c;
  abcd[705] = 2.0E0*I_ERI_Fx2z_S_Dyz_Pz_C1001000003_c;
  abcd[706] = 2.0E0*I_ERI_F3y_S_Dyz_Pz_C1001000003_c;
  abcd[707] = 2.0E0*I_ERI_F2yz_S_Dyz_Pz_C1001000003_c;
  abcd[708] = 2.0E0*I_ERI_Fy2z_S_Dyz_Pz_C1001000003_c;
  abcd[709] = 2.0E0*I_ERI_F3z_S_Dyz_Pz_C1001000003_c;
  abcd[710] = 2.0E0*I_ERI_F3x_S_D2z_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[711] = 2.0E0*I_ERI_F2xy_S_D2z_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[712] = 2.0E0*I_ERI_F2xz_S_D2z_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[713] = 2.0E0*I_ERI_Fx2y_S_D2z_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[714] = 2.0E0*I_ERI_Fxyz_S_D2z_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[715] = 2.0E0*I_ERI_Fx2z_S_D2z_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[716] = 2.0E0*I_ERI_F3y_S_D2z_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[717] = 2.0E0*I_ERI_F2yz_S_D2z_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[718] = 2.0E0*I_ERI_Fy2z_S_D2z_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[719] = 2.0E0*I_ERI_F3z_S_D2z_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000003_d
   ************************************************************/
  abcd[720] = 2.0E0*I_ERI_F3x_S_Px_Px_C1000003_d;
  abcd[721] = 2.0E0*I_ERI_F2xy_S_Px_Px_C1000003_d;
  abcd[722] = 2.0E0*I_ERI_F2xz_S_Px_Px_C1000003_d;
  abcd[723] = 2.0E0*I_ERI_Fx2y_S_Px_Px_C1000003_d;
  abcd[724] = 2.0E0*I_ERI_Fxyz_S_Px_Px_C1000003_d;
  abcd[725] = 2.0E0*I_ERI_Fx2z_S_Px_Px_C1000003_d;
  abcd[726] = 2.0E0*I_ERI_F3y_S_Px_Px_C1000003_d;
  abcd[727] = 2.0E0*I_ERI_F2yz_S_Px_Px_C1000003_d;
  abcd[728] = 2.0E0*I_ERI_Fy2z_S_Px_Px_C1000003_d;
  abcd[729] = 2.0E0*I_ERI_F3z_S_Px_Px_C1000003_d;
  abcd[730] = 2.0E0*I_ERI_F3x_S_Py_Px_C1000003_d;
  abcd[731] = 2.0E0*I_ERI_F2xy_S_Py_Px_C1000003_d;
  abcd[732] = 2.0E0*I_ERI_F2xz_S_Py_Px_C1000003_d;
  abcd[733] = 2.0E0*I_ERI_Fx2y_S_Py_Px_C1000003_d;
  abcd[734] = 2.0E0*I_ERI_Fxyz_S_Py_Px_C1000003_d;
  abcd[735] = 2.0E0*I_ERI_Fx2z_S_Py_Px_C1000003_d;
  abcd[736] = 2.0E0*I_ERI_F3y_S_Py_Px_C1000003_d;
  abcd[737] = 2.0E0*I_ERI_F2yz_S_Py_Px_C1000003_d;
  abcd[738] = 2.0E0*I_ERI_Fy2z_S_Py_Px_C1000003_d;
  abcd[739] = 2.0E0*I_ERI_F3z_S_Py_Px_C1000003_d;
  abcd[740] = 2.0E0*I_ERI_F3x_S_Pz_Px_C1000003_d;
  abcd[741] = 2.0E0*I_ERI_F2xy_S_Pz_Px_C1000003_d;
  abcd[742] = 2.0E0*I_ERI_F2xz_S_Pz_Px_C1000003_d;
  abcd[743] = 2.0E0*I_ERI_Fx2y_S_Pz_Px_C1000003_d;
  abcd[744] = 2.0E0*I_ERI_Fxyz_S_Pz_Px_C1000003_d;
  abcd[745] = 2.0E0*I_ERI_Fx2z_S_Pz_Px_C1000003_d;
  abcd[746] = 2.0E0*I_ERI_F3y_S_Pz_Px_C1000003_d;
  abcd[747] = 2.0E0*I_ERI_F2yz_S_Pz_Px_C1000003_d;
  abcd[748] = 2.0E0*I_ERI_Fy2z_S_Pz_Px_C1000003_d;
  abcd[749] = 2.0E0*I_ERI_F3z_S_Pz_Px_C1000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001000003
   ************************************************************/
  abcd[750] = 2.0E0*I_ERI_F3x_S_Px_D2x_C1001000003_d-1*I_ERI_F3x_S_Px_S_C1001000003;
  abcd[751] = 2.0E0*I_ERI_F2xy_S_Px_D2x_C1001000003_d-1*I_ERI_F2xy_S_Px_S_C1001000003;
  abcd[752] = 2.0E0*I_ERI_F2xz_S_Px_D2x_C1001000003_d-1*I_ERI_F2xz_S_Px_S_C1001000003;
  abcd[753] = 2.0E0*I_ERI_Fx2y_S_Px_D2x_C1001000003_d-1*I_ERI_Fx2y_S_Px_S_C1001000003;
  abcd[754] = 2.0E0*I_ERI_Fxyz_S_Px_D2x_C1001000003_d-1*I_ERI_Fxyz_S_Px_S_C1001000003;
  abcd[755] = 2.0E0*I_ERI_Fx2z_S_Px_D2x_C1001000003_d-1*I_ERI_Fx2z_S_Px_S_C1001000003;
  abcd[756] = 2.0E0*I_ERI_F3y_S_Px_D2x_C1001000003_d-1*I_ERI_F3y_S_Px_S_C1001000003;
  abcd[757] = 2.0E0*I_ERI_F2yz_S_Px_D2x_C1001000003_d-1*I_ERI_F2yz_S_Px_S_C1001000003;
  abcd[758] = 2.0E0*I_ERI_Fy2z_S_Px_D2x_C1001000003_d-1*I_ERI_Fy2z_S_Px_S_C1001000003;
  abcd[759] = 2.0E0*I_ERI_F3z_S_Px_D2x_C1001000003_d-1*I_ERI_F3z_S_Px_S_C1001000003;
  abcd[760] = 2.0E0*I_ERI_F3x_S_Py_D2x_C1001000003_d-1*I_ERI_F3x_S_Py_S_C1001000003;
  abcd[761] = 2.0E0*I_ERI_F2xy_S_Py_D2x_C1001000003_d-1*I_ERI_F2xy_S_Py_S_C1001000003;
  abcd[762] = 2.0E0*I_ERI_F2xz_S_Py_D2x_C1001000003_d-1*I_ERI_F2xz_S_Py_S_C1001000003;
  abcd[763] = 2.0E0*I_ERI_Fx2y_S_Py_D2x_C1001000003_d-1*I_ERI_Fx2y_S_Py_S_C1001000003;
  abcd[764] = 2.0E0*I_ERI_Fxyz_S_Py_D2x_C1001000003_d-1*I_ERI_Fxyz_S_Py_S_C1001000003;
  abcd[765] = 2.0E0*I_ERI_Fx2z_S_Py_D2x_C1001000003_d-1*I_ERI_Fx2z_S_Py_S_C1001000003;
  abcd[766] = 2.0E0*I_ERI_F3y_S_Py_D2x_C1001000003_d-1*I_ERI_F3y_S_Py_S_C1001000003;
  abcd[767] = 2.0E0*I_ERI_F2yz_S_Py_D2x_C1001000003_d-1*I_ERI_F2yz_S_Py_S_C1001000003;
  abcd[768] = 2.0E0*I_ERI_Fy2z_S_Py_D2x_C1001000003_d-1*I_ERI_Fy2z_S_Py_S_C1001000003;
  abcd[769] = 2.0E0*I_ERI_F3z_S_Py_D2x_C1001000003_d-1*I_ERI_F3z_S_Py_S_C1001000003;
  abcd[770] = 2.0E0*I_ERI_F3x_S_Pz_D2x_C1001000003_d-1*I_ERI_F3x_S_Pz_S_C1001000003;
  abcd[771] = 2.0E0*I_ERI_F2xy_S_Pz_D2x_C1001000003_d-1*I_ERI_F2xy_S_Pz_S_C1001000003;
  abcd[772] = 2.0E0*I_ERI_F2xz_S_Pz_D2x_C1001000003_d-1*I_ERI_F2xz_S_Pz_S_C1001000003;
  abcd[773] = 2.0E0*I_ERI_Fx2y_S_Pz_D2x_C1001000003_d-1*I_ERI_Fx2y_S_Pz_S_C1001000003;
  abcd[774] = 2.0E0*I_ERI_Fxyz_S_Pz_D2x_C1001000003_d-1*I_ERI_Fxyz_S_Pz_S_C1001000003;
  abcd[775] = 2.0E0*I_ERI_Fx2z_S_Pz_D2x_C1001000003_d-1*I_ERI_Fx2z_S_Pz_S_C1001000003;
  abcd[776] = 2.0E0*I_ERI_F3y_S_Pz_D2x_C1001000003_d-1*I_ERI_F3y_S_Pz_S_C1001000003;
  abcd[777] = 2.0E0*I_ERI_F2yz_S_Pz_D2x_C1001000003_d-1*I_ERI_F2yz_S_Pz_S_C1001000003;
  abcd[778] = 2.0E0*I_ERI_Fy2z_S_Pz_D2x_C1001000003_d-1*I_ERI_Fy2z_S_Pz_S_C1001000003;
  abcd[779] = 2.0E0*I_ERI_F3z_S_Pz_D2x_C1001000003_d-1*I_ERI_F3z_S_Pz_S_C1001000003;
  abcd[780] = 2.0E0*I_ERI_F3x_S_Px_Dxy_C1001000003_d;
  abcd[781] = 2.0E0*I_ERI_F2xy_S_Px_Dxy_C1001000003_d;
  abcd[782] = 2.0E0*I_ERI_F2xz_S_Px_Dxy_C1001000003_d;
  abcd[783] = 2.0E0*I_ERI_Fx2y_S_Px_Dxy_C1001000003_d;
  abcd[784] = 2.0E0*I_ERI_Fxyz_S_Px_Dxy_C1001000003_d;
  abcd[785] = 2.0E0*I_ERI_Fx2z_S_Px_Dxy_C1001000003_d;
  abcd[786] = 2.0E0*I_ERI_F3y_S_Px_Dxy_C1001000003_d;
  abcd[787] = 2.0E0*I_ERI_F2yz_S_Px_Dxy_C1001000003_d;
  abcd[788] = 2.0E0*I_ERI_Fy2z_S_Px_Dxy_C1001000003_d;
  abcd[789] = 2.0E0*I_ERI_F3z_S_Px_Dxy_C1001000003_d;
  abcd[790] = 2.0E0*I_ERI_F3x_S_Py_Dxy_C1001000003_d;
  abcd[791] = 2.0E0*I_ERI_F2xy_S_Py_Dxy_C1001000003_d;
  abcd[792] = 2.0E0*I_ERI_F2xz_S_Py_Dxy_C1001000003_d;
  abcd[793] = 2.0E0*I_ERI_Fx2y_S_Py_Dxy_C1001000003_d;
  abcd[794] = 2.0E0*I_ERI_Fxyz_S_Py_Dxy_C1001000003_d;
  abcd[795] = 2.0E0*I_ERI_Fx2z_S_Py_Dxy_C1001000003_d;
  abcd[796] = 2.0E0*I_ERI_F3y_S_Py_Dxy_C1001000003_d;
  abcd[797] = 2.0E0*I_ERI_F2yz_S_Py_Dxy_C1001000003_d;
  abcd[798] = 2.0E0*I_ERI_Fy2z_S_Py_Dxy_C1001000003_d;
  abcd[799] = 2.0E0*I_ERI_F3z_S_Py_Dxy_C1001000003_d;
  abcd[800] = 2.0E0*I_ERI_F3x_S_Pz_Dxy_C1001000003_d;
  abcd[801] = 2.0E0*I_ERI_F2xy_S_Pz_Dxy_C1001000003_d;
  abcd[802] = 2.0E0*I_ERI_F2xz_S_Pz_Dxy_C1001000003_d;
  abcd[803] = 2.0E0*I_ERI_Fx2y_S_Pz_Dxy_C1001000003_d;
  abcd[804] = 2.0E0*I_ERI_Fxyz_S_Pz_Dxy_C1001000003_d;
  abcd[805] = 2.0E0*I_ERI_Fx2z_S_Pz_Dxy_C1001000003_d;
  abcd[806] = 2.0E0*I_ERI_F3y_S_Pz_Dxy_C1001000003_d;
  abcd[807] = 2.0E0*I_ERI_F2yz_S_Pz_Dxy_C1001000003_d;
  abcd[808] = 2.0E0*I_ERI_Fy2z_S_Pz_Dxy_C1001000003_d;
  abcd[809] = 2.0E0*I_ERI_F3z_S_Pz_Dxy_C1001000003_d;
  abcd[810] = 2.0E0*I_ERI_F3x_S_Px_Dxz_C1001000003_d;
  abcd[811] = 2.0E0*I_ERI_F2xy_S_Px_Dxz_C1001000003_d;
  abcd[812] = 2.0E0*I_ERI_F2xz_S_Px_Dxz_C1001000003_d;
  abcd[813] = 2.0E0*I_ERI_Fx2y_S_Px_Dxz_C1001000003_d;
  abcd[814] = 2.0E0*I_ERI_Fxyz_S_Px_Dxz_C1001000003_d;
  abcd[815] = 2.0E0*I_ERI_Fx2z_S_Px_Dxz_C1001000003_d;
  abcd[816] = 2.0E0*I_ERI_F3y_S_Px_Dxz_C1001000003_d;
  abcd[817] = 2.0E0*I_ERI_F2yz_S_Px_Dxz_C1001000003_d;
  abcd[818] = 2.0E0*I_ERI_Fy2z_S_Px_Dxz_C1001000003_d;
  abcd[819] = 2.0E0*I_ERI_F3z_S_Px_Dxz_C1001000003_d;
  abcd[820] = 2.0E0*I_ERI_F3x_S_Py_Dxz_C1001000003_d;
  abcd[821] = 2.0E0*I_ERI_F2xy_S_Py_Dxz_C1001000003_d;
  abcd[822] = 2.0E0*I_ERI_F2xz_S_Py_Dxz_C1001000003_d;
  abcd[823] = 2.0E0*I_ERI_Fx2y_S_Py_Dxz_C1001000003_d;
  abcd[824] = 2.0E0*I_ERI_Fxyz_S_Py_Dxz_C1001000003_d;
  abcd[825] = 2.0E0*I_ERI_Fx2z_S_Py_Dxz_C1001000003_d;
  abcd[826] = 2.0E0*I_ERI_F3y_S_Py_Dxz_C1001000003_d;
  abcd[827] = 2.0E0*I_ERI_F2yz_S_Py_Dxz_C1001000003_d;
  abcd[828] = 2.0E0*I_ERI_Fy2z_S_Py_Dxz_C1001000003_d;
  abcd[829] = 2.0E0*I_ERI_F3z_S_Py_Dxz_C1001000003_d;
  abcd[830] = 2.0E0*I_ERI_F3x_S_Pz_Dxz_C1001000003_d;
  abcd[831] = 2.0E0*I_ERI_F2xy_S_Pz_Dxz_C1001000003_d;
  abcd[832] = 2.0E0*I_ERI_F2xz_S_Pz_Dxz_C1001000003_d;
  abcd[833] = 2.0E0*I_ERI_Fx2y_S_Pz_Dxz_C1001000003_d;
  abcd[834] = 2.0E0*I_ERI_Fxyz_S_Pz_Dxz_C1001000003_d;
  abcd[835] = 2.0E0*I_ERI_Fx2z_S_Pz_Dxz_C1001000003_d;
  abcd[836] = 2.0E0*I_ERI_F3y_S_Pz_Dxz_C1001000003_d;
  abcd[837] = 2.0E0*I_ERI_F2yz_S_Pz_Dxz_C1001000003_d;
  abcd[838] = 2.0E0*I_ERI_Fy2z_S_Pz_Dxz_C1001000003_d;
  abcd[839] = 2.0E0*I_ERI_F3z_S_Pz_Dxz_C1001000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000003_d
   ************************************************************/
  abcd[840] = 2.0E0*I_ERI_F3x_S_Px_Py_C1000003_d;
  abcd[841] = 2.0E0*I_ERI_F2xy_S_Px_Py_C1000003_d;
  abcd[842] = 2.0E0*I_ERI_F2xz_S_Px_Py_C1000003_d;
  abcd[843] = 2.0E0*I_ERI_Fx2y_S_Px_Py_C1000003_d;
  abcd[844] = 2.0E0*I_ERI_Fxyz_S_Px_Py_C1000003_d;
  abcd[845] = 2.0E0*I_ERI_Fx2z_S_Px_Py_C1000003_d;
  abcd[846] = 2.0E0*I_ERI_F3y_S_Px_Py_C1000003_d;
  abcd[847] = 2.0E0*I_ERI_F2yz_S_Px_Py_C1000003_d;
  abcd[848] = 2.0E0*I_ERI_Fy2z_S_Px_Py_C1000003_d;
  abcd[849] = 2.0E0*I_ERI_F3z_S_Px_Py_C1000003_d;
  abcd[850] = 2.0E0*I_ERI_F3x_S_Py_Py_C1000003_d;
  abcd[851] = 2.0E0*I_ERI_F2xy_S_Py_Py_C1000003_d;
  abcd[852] = 2.0E0*I_ERI_F2xz_S_Py_Py_C1000003_d;
  abcd[853] = 2.0E0*I_ERI_Fx2y_S_Py_Py_C1000003_d;
  abcd[854] = 2.0E0*I_ERI_Fxyz_S_Py_Py_C1000003_d;
  abcd[855] = 2.0E0*I_ERI_Fx2z_S_Py_Py_C1000003_d;
  abcd[856] = 2.0E0*I_ERI_F3y_S_Py_Py_C1000003_d;
  abcd[857] = 2.0E0*I_ERI_F2yz_S_Py_Py_C1000003_d;
  abcd[858] = 2.0E0*I_ERI_Fy2z_S_Py_Py_C1000003_d;
  abcd[859] = 2.0E0*I_ERI_F3z_S_Py_Py_C1000003_d;
  abcd[860] = 2.0E0*I_ERI_F3x_S_Pz_Py_C1000003_d;
  abcd[861] = 2.0E0*I_ERI_F2xy_S_Pz_Py_C1000003_d;
  abcd[862] = 2.0E0*I_ERI_F2xz_S_Pz_Py_C1000003_d;
  abcd[863] = 2.0E0*I_ERI_Fx2y_S_Pz_Py_C1000003_d;
  abcd[864] = 2.0E0*I_ERI_Fxyz_S_Pz_Py_C1000003_d;
  abcd[865] = 2.0E0*I_ERI_Fx2z_S_Pz_Py_C1000003_d;
  abcd[866] = 2.0E0*I_ERI_F3y_S_Pz_Py_C1000003_d;
  abcd[867] = 2.0E0*I_ERI_F2yz_S_Pz_Py_C1000003_d;
  abcd[868] = 2.0E0*I_ERI_Fy2z_S_Pz_Py_C1000003_d;
  abcd[869] = 2.0E0*I_ERI_F3z_S_Pz_Py_C1000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001000003
   ************************************************************/
  abcd[870] = 2.0E0*I_ERI_F3x_S_Px_Dxy_C1001000003_d;
  abcd[871] = 2.0E0*I_ERI_F2xy_S_Px_Dxy_C1001000003_d;
  abcd[872] = 2.0E0*I_ERI_F2xz_S_Px_Dxy_C1001000003_d;
  abcd[873] = 2.0E0*I_ERI_Fx2y_S_Px_Dxy_C1001000003_d;
  abcd[874] = 2.0E0*I_ERI_Fxyz_S_Px_Dxy_C1001000003_d;
  abcd[875] = 2.0E0*I_ERI_Fx2z_S_Px_Dxy_C1001000003_d;
  abcd[876] = 2.0E0*I_ERI_F3y_S_Px_Dxy_C1001000003_d;
  abcd[877] = 2.0E0*I_ERI_F2yz_S_Px_Dxy_C1001000003_d;
  abcd[878] = 2.0E0*I_ERI_Fy2z_S_Px_Dxy_C1001000003_d;
  abcd[879] = 2.0E0*I_ERI_F3z_S_Px_Dxy_C1001000003_d;
  abcd[880] = 2.0E0*I_ERI_F3x_S_Py_Dxy_C1001000003_d;
  abcd[881] = 2.0E0*I_ERI_F2xy_S_Py_Dxy_C1001000003_d;
  abcd[882] = 2.0E0*I_ERI_F2xz_S_Py_Dxy_C1001000003_d;
  abcd[883] = 2.0E0*I_ERI_Fx2y_S_Py_Dxy_C1001000003_d;
  abcd[884] = 2.0E0*I_ERI_Fxyz_S_Py_Dxy_C1001000003_d;
  abcd[885] = 2.0E0*I_ERI_Fx2z_S_Py_Dxy_C1001000003_d;
  abcd[886] = 2.0E0*I_ERI_F3y_S_Py_Dxy_C1001000003_d;
  abcd[887] = 2.0E0*I_ERI_F2yz_S_Py_Dxy_C1001000003_d;
  abcd[888] = 2.0E0*I_ERI_Fy2z_S_Py_Dxy_C1001000003_d;
  abcd[889] = 2.0E0*I_ERI_F3z_S_Py_Dxy_C1001000003_d;
  abcd[890] = 2.0E0*I_ERI_F3x_S_Pz_Dxy_C1001000003_d;
  abcd[891] = 2.0E0*I_ERI_F2xy_S_Pz_Dxy_C1001000003_d;
  abcd[892] = 2.0E0*I_ERI_F2xz_S_Pz_Dxy_C1001000003_d;
  abcd[893] = 2.0E0*I_ERI_Fx2y_S_Pz_Dxy_C1001000003_d;
  abcd[894] = 2.0E0*I_ERI_Fxyz_S_Pz_Dxy_C1001000003_d;
  abcd[895] = 2.0E0*I_ERI_Fx2z_S_Pz_Dxy_C1001000003_d;
  abcd[896] = 2.0E0*I_ERI_F3y_S_Pz_Dxy_C1001000003_d;
  abcd[897] = 2.0E0*I_ERI_F2yz_S_Pz_Dxy_C1001000003_d;
  abcd[898] = 2.0E0*I_ERI_Fy2z_S_Pz_Dxy_C1001000003_d;
  abcd[899] = 2.0E0*I_ERI_F3z_S_Pz_Dxy_C1001000003_d;
  abcd[900] = 2.0E0*I_ERI_F3x_S_Px_D2y_C1001000003_d-1*I_ERI_F3x_S_Px_S_C1001000003;
  abcd[901] = 2.0E0*I_ERI_F2xy_S_Px_D2y_C1001000003_d-1*I_ERI_F2xy_S_Px_S_C1001000003;
  abcd[902] = 2.0E0*I_ERI_F2xz_S_Px_D2y_C1001000003_d-1*I_ERI_F2xz_S_Px_S_C1001000003;
  abcd[903] = 2.0E0*I_ERI_Fx2y_S_Px_D2y_C1001000003_d-1*I_ERI_Fx2y_S_Px_S_C1001000003;
  abcd[904] = 2.0E0*I_ERI_Fxyz_S_Px_D2y_C1001000003_d-1*I_ERI_Fxyz_S_Px_S_C1001000003;
  abcd[905] = 2.0E0*I_ERI_Fx2z_S_Px_D2y_C1001000003_d-1*I_ERI_Fx2z_S_Px_S_C1001000003;
  abcd[906] = 2.0E0*I_ERI_F3y_S_Px_D2y_C1001000003_d-1*I_ERI_F3y_S_Px_S_C1001000003;
  abcd[907] = 2.0E0*I_ERI_F2yz_S_Px_D2y_C1001000003_d-1*I_ERI_F2yz_S_Px_S_C1001000003;
  abcd[908] = 2.0E0*I_ERI_Fy2z_S_Px_D2y_C1001000003_d-1*I_ERI_Fy2z_S_Px_S_C1001000003;
  abcd[909] = 2.0E0*I_ERI_F3z_S_Px_D2y_C1001000003_d-1*I_ERI_F3z_S_Px_S_C1001000003;
  abcd[910] = 2.0E0*I_ERI_F3x_S_Py_D2y_C1001000003_d-1*I_ERI_F3x_S_Py_S_C1001000003;
  abcd[911] = 2.0E0*I_ERI_F2xy_S_Py_D2y_C1001000003_d-1*I_ERI_F2xy_S_Py_S_C1001000003;
  abcd[912] = 2.0E0*I_ERI_F2xz_S_Py_D2y_C1001000003_d-1*I_ERI_F2xz_S_Py_S_C1001000003;
  abcd[913] = 2.0E0*I_ERI_Fx2y_S_Py_D2y_C1001000003_d-1*I_ERI_Fx2y_S_Py_S_C1001000003;
  abcd[914] = 2.0E0*I_ERI_Fxyz_S_Py_D2y_C1001000003_d-1*I_ERI_Fxyz_S_Py_S_C1001000003;
  abcd[915] = 2.0E0*I_ERI_Fx2z_S_Py_D2y_C1001000003_d-1*I_ERI_Fx2z_S_Py_S_C1001000003;
  abcd[916] = 2.0E0*I_ERI_F3y_S_Py_D2y_C1001000003_d-1*I_ERI_F3y_S_Py_S_C1001000003;
  abcd[917] = 2.0E0*I_ERI_F2yz_S_Py_D2y_C1001000003_d-1*I_ERI_F2yz_S_Py_S_C1001000003;
  abcd[918] = 2.0E0*I_ERI_Fy2z_S_Py_D2y_C1001000003_d-1*I_ERI_Fy2z_S_Py_S_C1001000003;
  abcd[919] = 2.0E0*I_ERI_F3z_S_Py_D2y_C1001000003_d-1*I_ERI_F3z_S_Py_S_C1001000003;
  abcd[920] = 2.0E0*I_ERI_F3x_S_Pz_D2y_C1001000003_d-1*I_ERI_F3x_S_Pz_S_C1001000003;
  abcd[921] = 2.0E0*I_ERI_F2xy_S_Pz_D2y_C1001000003_d-1*I_ERI_F2xy_S_Pz_S_C1001000003;
  abcd[922] = 2.0E0*I_ERI_F2xz_S_Pz_D2y_C1001000003_d-1*I_ERI_F2xz_S_Pz_S_C1001000003;
  abcd[923] = 2.0E0*I_ERI_Fx2y_S_Pz_D2y_C1001000003_d-1*I_ERI_Fx2y_S_Pz_S_C1001000003;
  abcd[924] = 2.0E0*I_ERI_Fxyz_S_Pz_D2y_C1001000003_d-1*I_ERI_Fxyz_S_Pz_S_C1001000003;
  abcd[925] = 2.0E0*I_ERI_Fx2z_S_Pz_D2y_C1001000003_d-1*I_ERI_Fx2z_S_Pz_S_C1001000003;
  abcd[926] = 2.0E0*I_ERI_F3y_S_Pz_D2y_C1001000003_d-1*I_ERI_F3y_S_Pz_S_C1001000003;
  abcd[927] = 2.0E0*I_ERI_F2yz_S_Pz_D2y_C1001000003_d-1*I_ERI_F2yz_S_Pz_S_C1001000003;
  abcd[928] = 2.0E0*I_ERI_Fy2z_S_Pz_D2y_C1001000003_d-1*I_ERI_Fy2z_S_Pz_S_C1001000003;
  abcd[929] = 2.0E0*I_ERI_F3z_S_Pz_D2y_C1001000003_d-1*I_ERI_F3z_S_Pz_S_C1001000003;
  abcd[930] = 2.0E0*I_ERI_F3x_S_Px_Dyz_C1001000003_d;
  abcd[931] = 2.0E0*I_ERI_F2xy_S_Px_Dyz_C1001000003_d;
  abcd[932] = 2.0E0*I_ERI_F2xz_S_Px_Dyz_C1001000003_d;
  abcd[933] = 2.0E0*I_ERI_Fx2y_S_Px_Dyz_C1001000003_d;
  abcd[934] = 2.0E0*I_ERI_Fxyz_S_Px_Dyz_C1001000003_d;
  abcd[935] = 2.0E0*I_ERI_Fx2z_S_Px_Dyz_C1001000003_d;
  abcd[936] = 2.0E0*I_ERI_F3y_S_Px_Dyz_C1001000003_d;
  abcd[937] = 2.0E0*I_ERI_F2yz_S_Px_Dyz_C1001000003_d;
  abcd[938] = 2.0E0*I_ERI_Fy2z_S_Px_Dyz_C1001000003_d;
  abcd[939] = 2.0E0*I_ERI_F3z_S_Px_Dyz_C1001000003_d;
  abcd[940] = 2.0E0*I_ERI_F3x_S_Py_Dyz_C1001000003_d;
  abcd[941] = 2.0E0*I_ERI_F2xy_S_Py_Dyz_C1001000003_d;
  abcd[942] = 2.0E0*I_ERI_F2xz_S_Py_Dyz_C1001000003_d;
  abcd[943] = 2.0E0*I_ERI_Fx2y_S_Py_Dyz_C1001000003_d;
  abcd[944] = 2.0E0*I_ERI_Fxyz_S_Py_Dyz_C1001000003_d;
  abcd[945] = 2.0E0*I_ERI_Fx2z_S_Py_Dyz_C1001000003_d;
  abcd[946] = 2.0E0*I_ERI_F3y_S_Py_Dyz_C1001000003_d;
  abcd[947] = 2.0E0*I_ERI_F2yz_S_Py_Dyz_C1001000003_d;
  abcd[948] = 2.0E0*I_ERI_Fy2z_S_Py_Dyz_C1001000003_d;
  abcd[949] = 2.0E0*I_ERI_F3z_S_Py_Dyz_C1001000003_d;
  abcd[950] = 2.0E0*I_ERI_F3x_S_Pz_Dyz_C1001000003_d;
  abcd[951] = 2.0E0*I_ERI_F2xy_S_Pz_Dyz_C1001000003_d;
  abcd[952] = 2.0E0*I_ERI_F2xz_S_Pz_Dyz_C1001000003_d;
  abcd[953] = 2.0E0*I_ERI_Fx2y_S_Pz_Dyz_C1001000003_d;
  abcd[954] = 2.0E0*I_ERI_Fxyz_S_Pz_Dyz_C1001000003_d;
  abcd[955] = 2.0E0*I_ERI_Fx2z_S_Pz_Dyz_C1001000003_d;
  abcd[956] = 2.0E0*I_ERI_F3y_S_Pz_Dyz_C1001000003_d;
  abcd[957] = 2.0E0*I_ERI_F2yz_S_Pz_Dyz_C1001000003_d;
  abcd[958] = 2.0E0*I_ERI_Fy2z_S_Pz_Dyz_C1001000003_d;
  abcd[959] = 2.0E0*I_ERI_F3z_S_Pz_Dyz_C1001000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000003_d
   ************************************************************/
  abcd[960] = 2.0E0*I_ERI_F3x_S_Px_Pz_C1000003_d;
  abcd[961] = 2.0E0*I_ERI_F2xy_S_Px_Pz_C1000003_d;
  abcd[962] = 2.0E0*I_ERI_F2xz_S_Px_Pz_C1000003_d;
  abcd[963] = 2.0E0*I_ERI_Fx2y_S_Px_Pz_C1000003_d;
  abcd[964] = 2.0E0*I_ERI_Fxyz_S_Px_Pz_C1000003_d;
  abcd[965] = 2.0E0*I_ERI_Fx2z_S_Px_Pz_C1000003_d;
  abcd[966] = 2.0E0*I_ERI_F3y_S_Px_Pz_C1000003_d;
  abcd[967] = 2.0E0*I_ERI_F2yz_S_Px_Pz_C1000003_d;
  abcd[968] = 2.0E0*I_ERI_Fy2z_S_Px_Pz_C1000003_d;
  abcd[969] = 2.0E0*I_ERI_F3z_S_Px_Pz_C1000003_d;
  abcd[970] = 2.0E0*I_ERI_F3x_S_Py_Pz_C1000003_d;
  abcd[971] = 2.0E0*I_ERI_F2xy_S_Py_Pz_C1000003_d;
  abcd[972] = 2.0E0*I_ERI_F2xz_S_Py_Pz_C1000003_d;
  abcd[973] = 2.0E0*I_ERI_Fx2y_S_Py_Pz_C1000003_d;
  abcd[974] = 2.0E0*I_ERI_Fxyz_S_Py_Pz_C1000003_d;
  abcd[975] = 2.0E0*I_ERI_Fx2z_S_Py_Pz_C1000003_d;
  abcd[976] = 2.0E0*I_ERI_F3y_S_Py_Pz_C1000003_d;
  abcd[977] = 2.0E0*I_ERI_F2yz_S_Py_Pz_C1000003_d;
  abcd[978] = 2.0E0*I_ERI_Fy2z_S_Py_Pz_C1000003_d;
  abcd[979] = 2.0E0*I_ERI_F3z_S_Py_Pz_C1000003_d;
  abcd[980] = 2.0E0*I_ERI_F3x_S_Pz_Pz_C1000003_d;
  abcd[981] = 2.0E0*I_ERI_F2xy_S_Pz_Pz_C1000003_d;
  abcd[982] = 2.0E0*I_ERI_F2xz_S_Pz_Pz_C1000003_d;
  abcd[983] = 2.0E0*I_ERI_Fx2y_S_Pz_Pz_C1000003_d;
  abcd[984] = 2.0E0*I_ERI_Fxyz_S_Pz_Pz_C1000003_d;
  abcd[985] = 2.0E0*I_ERI_Fx2z_S_Pz_Pz_C1000003_d;
  abcd[986] = 2.0E0*I_ERI_F3y_S_Pz_Pz_C1000003_d;
  abcd[987] = 2.0E0*I_ERI_F2yz_S_Pz_Pz_C1000003_d;
  abcd[988] = 2.0E0*I_ERI_Fy2z_S_Pz_Pz_C1000003_d;
  abcd[989] = 2.0E0*I_ERI_F3z_S_Pz_Pz_C1000003_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_D_C1001000003_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001000003
   ************************************************************/
  abcd[990] = 2.0E0*I_ERI_F3x_S_Px_Dxz_C1001000003_d;
  abcd[991] = 2.0E0*I_ERI_F2xy_S_Px_Dxz_C1001000003_d;
  abcd[992] = 2.0E0*I_ERI_F2xz_S_Px_Dxz_C1001000003_d;
  abcd[993] = 2.0E0*I_ERI_Fx2y_S_Px_Dxz_C1001000003_d;
  abcd[994] = 2.0E0*I_ERI_Fxyz_S_Px_Dxz_C1001000003_d;
  abcd[995] = 2.0E0*I_ERI_Fx2z_S_Px_Dxz_C1001000003_d;
  abcd[996] = 2.0E0*I_ERI_F3y_S_Px_Dxz_C1001000003_d;
  abcd[997] = 2.0E0*I_ERI_F2yz_S_Px_Dxz_C1001000003_d;
  abcd[998] = 2.0E0*I_ERI_Fy2z_S_Px_Dxz_C1001000003_d;
  abcd[999] = 2.0E0*I_ERI_F3z_S_Px_Dxz_C1001000003_d;
  abcd[1000] = 2.0E0*I_ERI_F3x_S_Py_Dxz_C1001000003_d;
  abcd[1001] = 2.0E0*I_ERI_F2xy_S_Py_Dxz_C1001000003_d;
  abcd[1002] = 2.0E0*I_ERI_F2xz_S_Py_Dxz_C1001000003_d;
  abcd[1003] = 2.0E0*I_ERI_Fx2y_S_Py_Dxz_C1001000003_d;
  abcd[1004] = 2.0E0*I_ERI_Fxyz_S_Py_Dxz_C1001000003_d;
  abcd[1005] = 2.0E0*I_ERI_Fx2z_S_Py_Dxz_C1001000003_d;
  abcd[1006] = 2.0E0*I_ERI_F3y_S_Py_Dxz_C1001000003_d;
  abcd[1007] = 2.0E0*I_ERI_F2yz_S_Py_Dxz_C1001000003_d;
  abcd[1008] = 2.0E0*I_ERI_Fy2z_S_Py_Dxz_C1001000003_d;
  abcd[1009] = 2.0E0*I_ERI_F3z_S_Py_Dxz_C1001000003_d;
  abcd[1010] = 2.0E0*I_ERI_F3x_S_Pz_Dxz_C1001000003_d;
  abcd[1011] = 2.0E0*I_ERI_F2xy_S_Pz_Dxz_C1001000003_d;
  abcd[1012] = 2.0E0*I_ERI_F2xz_S_Pz_Dxz_C1001000003_d;
  abcd[1013] = 2.0E0*I_ERI_Fx2y_S_Pz_Dxz_C1001000003_d;
  abcd[1014] = 2.0E0*I_ERI_Fxyz_S_Pz_Dxz_C1001000003_d;
  abcd[1015] = 2.0E0*I_ERI_Fx2z_S_Pz_Dxz_C1001000003_d;
  abcd[1016] = 2.0E0*I_ERI_F3y_S_Pz_Dxz_C1001000003_d;
  abcd[1017] = 2.0E0*I_ERI_F2yz_S_Pz_Dxz_C1001000003_d;
  abcd[1018] = 2.0E0*I_ERI_Fy2z_S_Pz_Dxz_C1001000003_d;
  abcd[1019] = 2.0E0*I_ERI_F3z_S_Pz_Dxz_C1001000003_d;
  abcd[1020] = 2.0E0*I_ERI_F3x_S_Px_Dyz_C1001000003_d;
  abcd[1021] = 2.0E0*I_ERI_F2xy_S_Px_Dyz_C1001000003_d;
  abcd[1022] = 2.0E0*I_ERI_F2xz_S_Px_Dyz_C1001000003_d;
  abcd[1023] = 2.0E0*I_ERI_Fx2y_S_Px_Dyz_C1001000003_d;
  abcd[1024] = 2.0E0*I_ERI_Fxyz_S_Px_Dyz_C1001000003_d;
  abcd[1025] = 2.0E0*I_ERI_Fx2z_S_Px_Dyz_C1001000003_d;
  abcd[1026] = 2.0E0*I_ERI_F3y_S_Px_Dyz_C1001000003_d;
  abcd[1027] = 2.0E0*I_ERI_F2yz_S_Px_Dyz_C1001000003_d;
  abcd[1028] = 2.0E0*I_ERI_Fy2z_S_Px_Dyz_C1001000003_d;
  abcd[1029] = 2.0E0*I_ERI_F3z_S_Px_Dyz_C1001000003_d;
  abcd[1030] = 2.0E0*I_ERI_F3x_S_Py_Dyz_C1001000003_d;
  abcd[1031] = 2.0E0*I_ERI_F2xy_S_Py_Dyz_C1001000003_d;
  abcd[1032] = 2.0E0*I_ERI_F2xz_S_Py_Dyz_C1001000003_d;
  abcd[1033] = 2.0E0*I_ERI_Fx2y_S_Py_Dyz_C1001000003_d;
  abcd[1034] = 2.0E0*I_ERI_Fxyz_S_Py_Dyz_C1001000003_d;
  abcd[1035] = 2.0E0*I_ERI_Fx2z_S_Py_Dyz_C1001000003_d;
  abcd[1036] = 2.0E0*I_ERI_F3y_S_Py_Dyz_C1001000003_d;
  abcd[1037] = 2.0E0*I_ERI_F2yz_S_Py_Dyz_C1001000003_d;
  abcd[1038] = 2.0E0*I_ERI_Fy2z_S_Py_Dyz_C1001000003_d;
  abcd[1039] = 2.0E0*I_ERI_F3z_S_Py_Dyz_C1001000003_d;
  abcd[1040] = 2.0E0*I_ERI_F3x_S_Pz_Dyz_C1001000003_d;
  abcd[1041] = 2.0E0*I_ERI_F2xy_S_Pz_Dyz_C1001000003_d;
  abcd[1042] = 2.0E0*I_ERI_F2xz_S_Pz_Dyz_C1001000003_d;
  abcd[1043] = 2.0E0*I_ERI_Fx2y_S_Pz_Dyz_C1001000003_d;
  abcd[1044] = 2.0E0*I_ERI_Fxyz_S_Pz_Dyz_C1001000003_d;
  abcd[1045] = 2.0E0*I_ERI_Fx2z_S_Pz_Dyz_C1001000003_d;
  abcd[1046] = 2.0E0*I_ERI_F3y_S_Pz_Dyz_C1001000003_d;
  abcd[1047] = 2.0E0*I_ERI_F2yz_S_Pz_Dyz_C1001000003_d;
  abcd[1048] = 2.0E0*I_ERI_Fy2z_S_Pz_Dyz_C1001000003_d;
  abcd[1049] = 2.0E0*I_ERI_F3z_S_Pz_Dyz_C1001000003_d;
  abcd[1050] = 2.0E0*I_ERI_F3x_S_Px_D2z_C1001000003_d-1*I_ERI_F3x_S_Px_S_C1001000003;
  abcd[1051] = 2.0E0*I_ERI_F2xy_S_Px_D2z_C1001000003_d-1*I_ERI_F2xy_S_Px_S_C1001000003;
  abcd[1052] = 2.0E0*I_ERI_F2xz_S_Px_D2z_C1001000003_d-1*I_ERI_F2xz_S_Px_S_C1001000003;
  abcd[1053] = 2.0E0*I_ERI_Fx2y_S_Px_D2z_C1001000003_d-1*I_ERI_Fx2y_S_Px_S_C1001000003;
  abcd[1054] = 2.0E0*I_ERI_Fxyz_S_Px_D2z_C1001000003_d-1*I_ERI_Fxyz_S_Px_S_C1001000003;
  abcd[1055] = 2.0E0*I_ERI_Fx2z_S_Px_D2z_C1001000003_d-1*I_ERI_Fx2z_S_Px_S_C1001000003;
  abcd[1056] = 2.0E0*I_ERI_F3y_S_Px_D2z_C1001000003_d-1*I_ERI_F3y_S_Px_S_C1001000003;
  abcd[1057] = 2.0E0*I_ERI_F2yz_S_Px_D2z_C1001000003_d-1*I_ERI_F2yz_S_Px_S_C1001000003;
  abcd[1058] = 2.0E0*I_ERI_Fy2z_S_Px_D2z_C1001000003_d-1*I_ERI_Fy2z_S_Px_S_C1001000003;
  abcd[1059] = 2.0E0*I_ERI_F3z_S_Px_D2z_C1001000003_d-1*I_ERI_F3z_S_Px_S_C1001000003;
  abcd[1060] = 2.0E0*I_ERI_F3x_S_Py_D2z_C1001000003_d-1*I_ERI_F3x_S_Py_S_C1001000003;
  abcd[1061] = 2.0E0*I_ERI_F2xy_S_Py_D2z_C1001000003_d-1*I_ERI_F2xy_S_Py_S_C1001000003;
  abcd[1062] = 2.0E0*I_ERI_F2xz_S_Py_D2z_C1001000003_d-1*I_ERI_F2xz_S_Py_S_C1001000003;
  abcd[1063] = 2.0E0*I_ERI_Fx2y_S_Py_D2z_C1001000003_d-1*I_ERI_Fx2y_S_Py_S_C1001000003;
  abcd[1064] = 2.0E0*I_ERI_Fxyz_S_Py_D2z_C1001000003_d-1*I_ERI_Fxyz_S_Py_S_C1001000003;
  abcd[1065] = 2.0E0*I_ERI_Fx2z_S_Py_D2z_C1001000003_d-1*I_ERI_Fx2z_S_Py_S_C1001000003;
  abcd[1066] = 2.0E0*I_ERI_F3y_S_Py_D2z_C1001000003_d-1*I_ERI_F3y_S_Py_S_C1001000003;
  abcd[1067] = 2.0E0*I_ERI_F2yz_S_Py_D2z_C1001000003_d-1*I_ERI_F2yz_S_Py_S_C1001000003;
  abcd[1068] = 2.0E0*I_ERI_Fy2z_S_Py_D2z_C1001000003_d-1*I_ERI_Fy2z_S_Py_S_C1001000003;
  abcd[1069] = 2.0E0*I_ERI_F3z_S_Py_D2z_C1001000003_d-1*I_ERI_F3z_S_Py_S_C1001000003;
  abcd[1070] = 2.0E0*I_ERI_F3x_S_Pz_D2z_C1001000003_d-1*I_ERI_F3x_S_Pz_S_C1001000003;
  abcd[1071] = 2.0E0*I_ERI_F2xy_S_Pz_D2z_C1001000003_d-1*I_ERI_F2xy_S_Pz_S_C1001000003;
  abcd[1072] = 2.0E0*I_ERI_F2xz_S_Pz_D2z_C1001000003_d-1*I_ERI_F2xz_S_Pz_S_C1001000003;
  abcd[1073] = 2.0E0*I_ERI_Fx2y_S_Pz_D2z_C1001000003_d-1*I_ERI_Fx2y_S_Pz_S_C1001000003;
  abcd[1074] = 2.0E0*I_ERI_Fxyz_S_Pz_D2z_C1001000003_d-1*I_ERI_Fxyz_S_Pz_S_C1001000003;
  abcd[1075] = 2.0E0*I_ERI_Fx2z_S_Pz_D2z_C1001000003_d-1*I_ERI_Fx2z_S_Pz_S_C1001000003;
  abcd[1076] = 2.0E0*I_ERI_F3y_S_Pz_D2z_C1001000003_d-1*I_ERI_F3y_S_Pz_S_C1001000003;
  abcd[1077] = 2.0E0*I_ERI_F2yz_S_Pz_D2z_C1001000003_d-1*I_ERI_F2yz_S_Pz_S_C1001000003;
  abcd[1078] = 2.0E0*I_ERI_Fy2z_S_Pz_D2z_C1001000003_d-1*I_ERI_Fy2z_S_Pz_S_C1001000003;
  abcd[1079] = 2.0E0*I_ERI_F3z_S_Pz_D2z_C1001000003_d-1*I_ERI_F3z_S_Pz_S_C1001000003;
}
