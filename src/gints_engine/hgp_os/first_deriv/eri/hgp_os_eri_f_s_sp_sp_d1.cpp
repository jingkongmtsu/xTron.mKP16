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
// BRA1 as redundant position, total RHS integrals evaluated as: 74992
// BRA2 as redundant position, total RHS integrals evaluated as: 77434
// KET1 as redundant position, total RHS integrals evaluated as: 69138
// KET2 as redundant position, total RHS integrals evaluated as: 69138
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

void hgp_os_eri_f_s_sp_sp_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_G4x_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C3_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C3_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G4y_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_G4z_S_S_S_C3_a = 0.0E0;
  Double I_ERI_D2x_S_S_S_C3 = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C3 = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C3 = 0.0E0;
  Double I_ERI_D2y_S_S_S_C3 = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C3 = 0.0E0;
  Double I_ERI_D2z_S_S_S_C3 = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1000003_a = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1000003_a = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1000003_a = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000003 = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000003 = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000003 = 0.0E0;
  Double I_ERI_G4x_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G3xy_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G3xz_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G4y_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G3yz_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G4z_S_S_Px_C1000000003_a = 0.0E0;
  Double I_ERI_G4x_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G3xy_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G3xz_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G4y_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G3yz_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G4z_S_S_Py_C1000000003_a = 0.0E0;
  Double I_ERI_G4x_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G3xy_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G3xz_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G4y_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G3yz_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_G4z_S_S_Pz_C1000000003_a = 0.0E0;
  Double I_ERI_D2x_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_Dxy_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_Dxz_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_D2y_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_Dyz_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_D2z_S_S_Px_C1000000003 = 0.0E0;
  Double I_ERI_D2x_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_Dxy_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_Dxz_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_D2y_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_Dyz_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_D2z_S_S_Py_C1000000003 = 0.0E0;
  Double I_ERI_D2x_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_Dxy_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_Dxz_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_D2y_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_Dyz_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_D2z_S_S_Pz_C1000000003 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C3_c = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C3_c = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C3_c = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C3_c = 0.0E0;
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
  Double I_ERI_G4x_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G4y_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_G4z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F3x_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C3_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C3_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_C3_b = 0.0E0;
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
  Double I_ERI_G4x_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G3xy_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G3xz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G4y_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G3yz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G4z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_G4x_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G3xy_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G3xz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G4y_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G3yz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G4z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_G4x_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G3xy_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G3xz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G4y_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G3yz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G4z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F3x_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F3y_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F3z_S_S_Px_C1000000003_b = 0.0E0;
  Double I_ERI_F3x_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F3y_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F3z_S_S_Py_C1000000003_b = 0.0E0;
  Double I_ERI_F3x_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F2xy_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F2xz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F3y_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F2yz_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_F3z_S_S_Pz_C1000000003_b = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1001000003_a = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001000003 = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001000003 = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000000003_c = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000000003_c = 0.0E0;
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
       * shell quartet name: SQ_ERI_G_S_S_P
       * expanding position: KET2
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_S_Px_vrr = QDX*I_ERI_G4x_S_S_S_vrr+WQX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_S_Px_vrr = QDX*I_ERI_G3xy_S_S_S_vrr+WQX*I_ERI_G3xy_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_S_Px_vrr = QDX*I_ERI_G3xz_S_S_S_vrr+WQX*I_ERI_G3xz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_S_Px_vrr = QDX*I_ERI_G2x2y_S_S_S_vrr+WQX*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_S_Px_vrr = QDX*I_ERI_G2xyz_S_S_S_vrr+WQX*I_ERI_G2xyz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_S_Px_vrr = QDX*I_ERI_G2x2z_S_S_S_vrr+WQX*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_S_Px_vrr = QDX*I_ERI_Gx3y_S_S_S_vrr+WQX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_S_Px_vrr = QDX*I_ERI_Gx2yz_S_S_S_vrr+WQX*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_S_Px_vrr = QDX*I_ERI_Gxy2z_S_S_S_vrr+WQX*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_S_Px_vrr = QDX*I_ERI_Gx3z_S_S_S_vrr+WQX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_S_Px_vrr = QDX*I_ERI_G4y_S_S_S_vrr+WQX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_S_Px_vrr = QDX*I_ERI_G3yz_S_S_S_vrr+WQX*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_S_Px_vrr = QDX*I_ERI_G2y2z_S_S_S_vrr+WQX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_S_Px_vrr = QDX*I_ERI_Gy3z_S_S_S_vrr+WQX*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_S_Px_vrr = QDX*I_ERI_G4z_S_S_S_vrr+WQX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_S_Py_vrr = QDY*I_ERI_G4x_S_S_S_vrr+WQY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_S_Py_vrr = QDY*I_ERI_G3xy_S_S_S_vrr+WQY*I_ERI_G3xy_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_S_Py_vrr = QDY*I_ERI_G3xz_S_S_S_vrr+WQY*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_S_Py_vrr = QDY*I_ERI_G2x2y_S_S_S_vrr+WQY*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_S_Py_vrr = QDY*I_ERI_G2xyz_S_S_S_vrr+WQY*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_S_Py_vrr = QDY*I_ERI_G2x2z_S_S_S_vrr+WQY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_S_Py_vrr = QDY*I_ERI_Gx3y_S_S_S_vrr+WQY*I_ERI_Gx3y_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_S_Py_vrr = QDY*I_ERI_Gx2yz_S_S_S_vrr+WQY*I_ERI_Gx2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_S_Py_vrr = QDY*I_ERI_Gxy2z_S_S_S_vrr+WQY*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_S_Py_vrr = QDY*I_ERI_Gx3z_S_S_S_vrr+WQY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_S_Py_vrr = QDY*I_ERI_G4y_S_S_S_vrr+WQY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_S_Py_vrr = QDY*I_ERI_G3yz_S_S_S_vrr+WQY*I_ERI_G3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_S_Py_vrr = QDY*I_ERI_G2y2z_S_S_S_vrr+WQY*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_S_Py_vrr = QDY*I_ERI_Gy3z_S_S_S_vrr+WQY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_S_Py_vrr = QDY*I_ERI_G4z_S_S_S_vrr+WQY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_S_Pz_vrr = QDZ*I_ERI_G4x_S_S_S_vrr+WQZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_S_Pz_vrr = QDZ*I_ERI_G3xy_S_S_S_vrr+WQZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_S_Pz_vrr = QDZ*I_ERI_G3xz_S_S_S_vrr+WQZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_S_Pz_vrr = QDZ*I_ERI_G2x2y_S_S_S_vrr+WQZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_S_Pz_vrr = QDZ*I_ERI_G2xyz_S_S_S_vrr+WQZ*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_S_Pz_vrr = QDZ*I_ERI_G2x2z_S_S_S_vrr+WQZ*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_S_Pz_vrr = QDZ*I_ERI_Gx3y_S_S_S_vrr+WQZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_S_Pz_vrr = QDZ*I_ERI_Gx2yz_S_S_S_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_S_Pz_vrr = QDZ*I_ERI_Gxy2z_S_S_S_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_S_Pz_vrr = QDZ*I_ERI_Gx3z_S_S_S_vrr+WQZ*I_ERI_Gx3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_S_Pz_vrr = QDZ*I_ERI_G4y_S_S_S_vrr+WQZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_S_Pz_vrr = QDZ*I_ERI_G3yz_S_S_S_vrr+WQZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_S_Pz_vrr = QDZ*I_ERI_G2y2z_S_S_S_vrr+WQZ*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_S_Pz_vrr = QDZ*I_ERI_Gy3z_S_S_S_vrr+WQZ*I_ERI_Gy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_S_Pz_vrr = QDZ*I_ERI_G4z_S_S_S_vrr+WQZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2k*I_ERI_F3z_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S_C3_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C3_a_coefs = ic2*jc2*alpha;
      I_ERI_G4x_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C3_a += SQ_ERI_G_S_S_S_C3_a_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C3
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C3_coefs = ic2*jc2;
      I_ERI_D2x_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C3 += SQ_ERI_D_S_S_S_C3_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1000003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1000003_a_coefs = ic2*jc2_1*alpha;
      I_ERI_G4x_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1000003_a += SQ_ERI_G_S_P_S_C1000003_a_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000003_coefs = ic2*jc2_1;
      I_ERI_D2x_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000003 += SQ_ERI_D_S_P_S_C1000003_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_P_C1000000003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_P_C1000000003_a_coefs = ic2*jc2_2*alpha;
      I_ERI_G4x_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4x_S_S_Px_vrr;
      I_ERI_G3xy_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xy_S_S_Px_vrr;
      I_ERI_G3xz_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xz_S_S_Px_vrr;
      I_ERI_G2x2y_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2y_S_S_Px_vrr;
      I_ERI_G2xyz_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2xyz_S_S_Px_vrr;
      I_ERI_G2x2z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2z_S_S_Px_vrr;
      I_ERI_Gx3y_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3y_S_S_Px_vrr;
      I_ERI_Gx2yz_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
      I_ERI_Gxy2z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
      I_ERI_Gx3z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3z_S_S_Px_vrr;
      I_ERI_G4y_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4y_S_S_Px_vrr;
      I_ERI_G3yz_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3yz_S_S_Px_vrr;
      I_ERI_G2y2z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2y2z_S_S_Px_vrr;
      I_ERI_Gy3z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gy3z_S_S_Px_vrr;
      I_ERI_G4z_S_S_Px_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4z_S_S_Px_vrr;
      I_ERI_G4x_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4x_S_S_Py_vrr;
      I_ERI_G3xy_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xy_S_S_Py_vrr;
      I_ERI_G3xz_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xz_S_S_Py_vrr;
      I_ERI_G2x2y_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2y_S_S_Py_vrr;
      I_ERI_G2xyz_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2xyz_S_S_Py_vrr;
      I_ERI_G2x2z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2z_S_S_Py_vrr;
      I_ERI_Gx3y_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3y_S_S_Py_vrr;
      I_ERI_Gx2yz_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
      I_ERI_Gxy2z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
      I_ERI_Gx3z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3z_S_S_Py_vrr;
      I_ERI_G4y_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4y_S_S_Py_vrr;
      I_ERI_G3yz_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3yz_S_S_Py_vrr;
      I_ERI_G2y2z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2y2z_S_S_Py_vrr;
      I_ERI_Gy3z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gy3z_S_S_Py_vrr;
      I_ERI_G4z_S_S_Py_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4z_S_S_Py_vrr;
      I_ERI_G4x_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4x_S_S_Pz_vrr;
      I_ERI_G3xy_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xy_S_S_Pz_vrr;
      I_ERI_G3xz_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3xz_S_S_Pz_vrr;
      I_ERI_G2x2y_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
      I_ERI_G2xyz_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
      I_ERI_G2x2z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
      I_ERI_Gx3y_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
      I_ERI_Gx2yz_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
      I_ERI_Gxy2z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
      I_ERI_Gx3z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
      I_ERI_G4y_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4y_S_S_Pz_vrr;
      I_ERI_G3yz_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G3yz_S_S_Pz_vrr;
      I_ERI_G2y2z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
      I_ERI_Gy3z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
      I_ERI_G4z_S_S_Pz_C1000000003_a += SQ_ERI_G_S_S_P_C1000000003_a_coefs*I_ERI_G4z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_P_C1000000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_P_C1000000003_coefs = ic2*jc2_2;
      I_ERI_D2x_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2x_S_S_Px_vrr;
      I_ERI_Dxy_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxy_S_S_Px_vrr;
      I_ERI_Dxz_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxz_S_S_Px_vrr;
      I_ERI_D2y_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2y_S_S_Px_vrr;
      I_ERI_Dyz_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dyz_S_S_Px_vrr;
      I_ERI_D2z_S_S_Px_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2z_S_S_Px_vrr;
      I_ERI_D2x_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2x_S_S_Py_vrr;
      I_ERI_Dxy_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxy_S_S_Py_vrr;
      I_ERI_Dxz_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxz_S_S_Py_vrr;
      I_ERI_D2y_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2y_S_S_Py_vrr;
      I_ERI_Dyz_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dyz_S_S_Py_vrr;
      I_ERI_D2z_S_S_Py_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2z_S_S_Py_vrr;
      I_ERI_D2x_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2x_S_S_Pz_vrr;
      I_ERI_Dxy_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxy_S_S_Pz_vrr;
      I_ERI_Dxz_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dxz_S_S_Pz_vrr;
      I_ERI_D2y_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2y_S_S_Pz_vrr;
      I_ERI_Dyz_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_Dyz_S_S_Pz_vrr;
      I_ERI_D2z_S_S_Pz_C1000000003 += SQ_ERI_D_S_S_P_C1000000003_coefs*I_ERI_D2z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C3_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C3_c_coefs = ic2*jc2*gamma;
      I_ERI_F3x_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C3_c += SQ_ERI_F_S_P_S_C3_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1000003_c_coefs = ic2*jc2_1*gamma;
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
      Double SQ_ERI_F_S_S_S_C1000003_coefs = ic2*jc2_1;
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
      Double SQ_ERI_F_S_S_P_C1001000003_coefs = ic2*jc2_3;
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
       * shell quartet name: SQ_ERI_G_S_S_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C3_b_coefs = ic2*jc2*beta;
      I_ERI_G4x_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C3_b += SQ_ERI_G_S_S_S_C3_b_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C3_b_coefs = ic2*jc2*beta;
      I_ERI_F3x_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C3_b += SQ_ERI_F_S_S_S_C3_b_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1000003_b_coefs = ic2*jc2_1*beta;
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
      Double SQ_ERI_F_S_P_S_C1000003_b_coefs = ic2*jc2_1*beta;
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
       * shell quartet name: SQ_ERI_G_S_S_P_C1000000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_P_C1000000003_b_coefs = ic2*jc2_2*beta;
      I_ERI_G4x_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4x_S_S_Px_vrr;
      I_ERI_G3xy_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xy_S_S_Px_vrr;
      I_ERI_G3xz_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xz_S_S_Px_vrr;
      I_ERI_G2x2y_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2y_S_S_Px_vrr;
      I_ERI_G2xyz_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2xyz_S_S_Px_vrr;
      I_ERI_G2x2z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2z_S_S_Px_vrr;
      I_ERI_Gx3y_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3y_S_S_Px_vrr;
      I_ERI_Gx2yz_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
      I_ERI_Gxy2z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
      I_ERI_Gx3z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3z_S_S_Px_vrr;
      I_ERI_G4y_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4y_S_S_Px_vrr;
      I_ERI_G3yz_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3yz_S_S_Px_vrr;
      I_ERI_G2y2z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2y2z_S_S_Px_vrr;
      I_ERI_Gy3z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gy3z_S_S_Px_vrr;
      I_ERI_G4z_S_S_Px_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4z_S_S_Px_vrr;
      I_ERI_G4x_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4x_S_S_Py_vrr;
      I_ERI_G3xy_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xy_S_S_Py_vrr;
      I_ERI_G3xz_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xz_S_S_Py_vrr;
      I_ERI_G2x2y_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2y_S_S_Py_vrr;
      I_ERI_G2xyz_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2xyz_S_S_Py_vrr;
      I_ERI_G2x2z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2z_S_S_Py_vrr;
      I_ERI_Gx3y_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3y_S_S_Py_vrr;
      I_ERI_Gx2yz_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
      I_ERI_Gxy2z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
      I_ERI_Gx3z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3z_S_S_Py_vrr;
      I_ERI_G4y_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4y_S_S_Py_vrr;
      I_ERI_G3yz_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3yz_S_S_Py_vrr;
      I_ERI_G2y2z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2y2z_S_S_Py_vrr;
      I_ERI_Gy3z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gy3z_S_S_Py_vrr;
      I_ERI_G4z_S_S_Py_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4z_S_S_Py_vrr;
      I_ERI_G4x_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4x_S_S_Pz_vrr;
      I_ERI_G3xy_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xy_S_S_Pz_vrr;
      I_ERI_G3xz_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3xz_S_S_Pz_vrr;
      I_ERI_G2x2y_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
      I_ERI_G2xyz_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
      I_ERI_G2x2z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
      I_ERI_Gx3y_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
      I_ERI_Gx2yz_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
      I_ERI_Gxy2z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
      I_ERI_Gx3z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
      I_ERI_G4y_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4y_S_S_Pz_vrr;
      I_ERI_G3yz_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G3yz_S_S_Pz_vrr;
      I_ERI_G2y2z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
      I_ERI_Gy3z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
      I_ERI_G4z_S_S_Pz_C1000000003_b += SQ_ERI_G_S_S_P_C1000000003_b_coefs*I_ERI_G4z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_P_C1000000003_b_coefs = ic2*jc2_2*beta;
      I_ERI_F3x_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3x_S_S_Px_vrr;
      I_ERI_F2xy_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xy_S_S_Px_vrr;
      I_ERI_F2xz_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xz_S_S_Px_vrr;
      I_ERI_Fx2y_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2y_S_S_Px_vrr;
      I_ERI_Fxyz_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fxyz_S_S_Px_vrr;
      I_ERI_Fx2z_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2z_S_S_Px_vrr;
      I_ERI_F3y_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3y_S_S_Px_vrr;
      I_ERI_F2yz_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2yz_S_S_Px_vrr;
      I_ERI_Fy2z_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fy2z_S_S_Px_vrr;
      I_ERI_F3z_S_S_Px_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3z_S_S_Px_vrr;
      I_ERI_F3x_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3x_S_S_Py_vrr;
      I_ERI_F2xy_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xy_S_S_Py_vrr;
      I_ERI_F2xz_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xz_S_S_Py_vrr;
      I_ERI_Fx2y_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2y_S_S_Py_vrr;
      I_ERI_Fxyz_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fxyz_S_S_Py_vrr;
      I_ERI_Fx2z_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2z_S_S_Py_vrr;
      I_ERI_F3y_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3y_S_S_Py_vrr;
      I_ERI_F2yz_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2yz_S_S_Py_vrr;
      I_ERI_Fy2z_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fy2z_S_S_Py_vrr;
      I_ERI_F3z_S_S_Py_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3z_S_S_Py_vrr;
      I_ERI_F3x_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3x_S_S_Pz_vrr;
      I_ERI_F2xy_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xy_S_S_Pz_vrr;
      I_ERI_F2xz_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2xz_S_S_Pz_vrr;
      I_ERI_Fx2y_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2y_S_S_Pz_vrr;
      I_ERI_Fxyz_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fxyz_S_S_Pz_vrr;
      I_ERI_Fx2z_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fx2z_S_S_Pz_vrr;
      I_ERI_F3y_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3y_S_S_Pz_vrr;
      I_ERI_F2yz_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F2yz_S_S_Pz_vrr;
      I_ERI_Fy2z_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_Fy2z_S_S_Pz_vrr;
      I_ERI_F3z_S_S_Pz_C1000000003_b += SQ_ERI_F_S_S_P_C1000000003_b_coefs*I_ERI_F3z_S_S_Pz_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_C1001000003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_C1001000003_a_coefs = ic2*jc2_3*alpha;
      I_ERI_G4x_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_C1001000003_a += SQ_ERI_G_S_D_S_C1001000003_a_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_C1001000003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1001000003_a_coefs = ic2*jc2_3*alpha;
      I_ERI_G4x_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1001000003_a += SQ_ERI_G_S_P_S_C1001000003_a_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001000003_coefs = ic2*jc2_3;
      I_ERI_D2x_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001000003 += SQ_ERI_D_S_D_S_C1001000003_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001000003_coefs = ic2*jc2_3;
      I_ERI_D2x_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001000003 += SQ_ERI_D_S_P_S_C1001000003_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1000000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1000000003_c_coefs = ic2*jc2_2*gamma;
      I_ERI_F3x_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1000000003_c += SQ_ERI_F_S_D_S_C1000000003_c_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1000000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000000003_c_coefs = ic2*jc2_2*gamma;
      I_ERI_F3x_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000000003_c += SQ_ERI_F_S_P_S_C1000000003_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1001000003_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1001000003_c_coefs = ic2*jc2_3*gamma;
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
      Double SQ_ERI_F_S_D_S_C1001000003_c_coefs = ic2*jc2_3*gamma;
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
       * shell quartet name: SQ_ERI_G_S_D_S_C1001000003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_C1001000003_b_coefs = ic2*jc2_3*beta;
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
      Double SQ_ERI_G_S_P_S_C1001000003_b_coefs = ic2*jc2_3*beta;
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
      Double SQ_ERI_F_S_D_S_C1001000003_b_coefs = ic2*jc2_3*beta;
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
      Double SQ_ERI_F_S_P_S_C1001000003_b_coefs = ic2*jc2_3*beta;
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
   * shell quartet name: SQ_ERI_D_S_P_P_C1001000003
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001000003
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001000003
   ************************************************************/
  Double I_ERI_D2x_S_Px_Px_C1001000003 = I_ERI_D2x_S_D2x_S_C1001000003+CDX*I_ERI_D2x_S_Px_S_C1001000003;
  Double I_ERI_Dxy_S_Px_Px_C1001000003 = I_ERI_Dxy_S_D2x_S_C1001000003+CDX*I_ERI_Dxy_S_Px_S_C1001000003;
  Double I_ERI_Dxz_S_Px_Px_C1001000003 = I_ERI_Dxz_S_D2x_S_C1001000003+CDX*I_ERI_Dxz_S_Px_S_C1001000003;
  Double I_ERI_D2y_S_Px_Px_C1001000003 = I_ERI_D2y_S_D2x_S_C1001000003+CDX*I_ERI_D2y_S_Px_S_C1001000003;
  Double I_ERI_Dyz_S_Px_Px_C1001000003 = I_ERI_Dyz_S_D2x_S_C1001000003+CDX*I_ERI_Dyz_S_Px_S_C1001000003;
  Double I_ERI_D2z_S_Px_Px_C1001000003 = I_ERI_D2z_S_D2x_S_C1001000003+CDX*I_ERI_D2z_S_Px_S_C1001000003;
  Double I_ERI_D2x_S_Py_Px_C1001000003 = I_ERI_D2x_S_Dxy_S_C1001000003+CDX*I_ERI_D2x_S_Py_S_C1001000003;
  Double I_ERI_Dxy_S_Py_Px_C1001000003 = I_ERI_Dxy_S_Dxy_S_C1001000003+CDX*I_ERI_Dxy_S_Py_S_C1001000003;
  Double I_ERI_Dxz_S_Py_Px_C1001000003 = I_ERI_Dxz_S_Dxy_S_C1001000003+CDX*I_ERI_Dxz_S_Py_S_C1001000003;
  Double I_ERI_D2y_S_Py_Px_C1001000003 = I_ERI_D2y_S_Dxy_S_C1001000003+CDX*I_ERI_D2y_S_Py_S_C1001000003;
  Double I_ERI_Dyz_S_Py_Px_C1001000003 = I_ERI_Dyz_S_Dxy_S_C1001000003+CDX*I_ERI_Dyz_S_Py_S_C1001000003;
  Double I_ERI_D2z_S_Py_Px_C1001000003 = I_ERI_D2z_S_Dxy_S_C1001000003+CDX*I_ERI_D2z_S_Py_S_C1001000003;
  Double I_ERI_D2x_S_Pz_Px_C1001000003 = I_ERI_D2x_S_Dxz_S_C1001000003+CDX*I_ERI_D2x_S_Pz_S_C1001000003;
  Double I_ERI_Dxy_S_Pz_Px_C1001000003 = I_ERI_Dxy_S_Dxz_S_C1001000003+CDX*I_ERI_Dxy_S_Pz_S_C1001000003;
  Double I_ERI_Dxz_S_Pz_Px_C1001000003 = I_ERI_Dxz_S_Dxz_S_C1001000003+CDX*I_ERI_Dxz_S_Pz_S_C1001000003;
  Double I_ERI_D2y_S_Pz_Px_C1001000003 = I_ERI_D2y_S_Dxz_S_C1001000003+CDX*I_ERI_D2y_S_Pz_S_C1001000003;
  Double I_ERI_Dyz_S_Pz_Px_C1001000003 = I_ERI_Dyz_S_Dxz_S_C1001000003+CDX*I_ERI_Dyz_S_Pz_S_C1001000003;
  Double I_ERI_D2z_S_Pz_Px_C1001000003 = I_ERI_D2z_S_Dxz_S_C1001000003+CDX*I_ERI_D2z_S_Pz_S_C1001000003;
  Double I_ERI_D2x_S_Px_Py_C1001000003 = I_ERI_D2x_S_Dxy_S_C1001000003+CDY*I_ERI_D2x_S_Px_S_C1001000003;
  Double I_ERI_Dxy_S_Px_Py_C1001000003 = I_ERI_Dxy_S_Dxy_S_C1001000003+CDY*I_ERI_Dxy_S_Px_S_C1001000003;
  Double I_ERI_Dxz_S_Px_Py_C1001000003 = I_ERI_Dxz_S_Dxy_S_C1001000003+CDY*I_ERI_Dxz_S_Px_S_C1001000003;
  Double I_ERI_D2y_S_Px_Py_C1001000003 = I_ERI_D2y_S_Dxy_S_C1001000003+CDY*I_ERI_D2y_S_Px_S_C1001000003;
  Double I_ERI_Dyz_S_Px_Py_C1001000003 = I_ERI_Dyz_S_Dxy_S_C1001000003+CDY*I_ERI_Dyz_S_Px_S_C1001000003;
  Double I_ERI_D2z_S_Px_Py_C1001000003 = I_ERI_D2z_S_Dxy_S_C1001000003+CDY*I_ERI_D2z_S_Px_S_C1001000003;
  Double I_ERI_D2x_S_Py_Py_C1001000003 = I_ERI_D2x_S_D2y_S_C1001000003+CDY*I_ERI_D2x_S_Py_S_C1001000003;
  Double I_ERI_Dxy_S_Py_Py_C1001000003 = I_ERI_Dxy_S_D2y_S_C1001000003+CDY*I_ERI_Dxy_S_Py_S_C1001000003;
  Double I_ERI_Dxz_S_Py_Py_C1001000003 = I_ERI_Dxz_S_D2y_S_C1001000003+CDY*I_ERI_Dxz_S_Py_S_C1001000003;
  Double I_ERI_D2y_S_Py_Py_C1001000003 = I_ERI_D2y_S_D2y_S_C1001000003+CDY*I_ERI_D2y_S_Py_S_C1001000003;
  Double I_ERI_Dyz_S_Py_Py_C1001000003 = I_ERI_Dyz_S_D2y_S_C1001000003+CDY*I_ERI_Dyz_S_Py_S_C1001000003;
  Double I_ERI_D2z_S_Py_Py_C1001000003 = I_ERI_D2z_S_D2y_S_C1001000003+CDY*I_ERI_D2z_S_Py_S_C1001000003;
  Double I_ERI_D2x_S_Pz_Py_C1001000003 = I_ERI_D2x_S_Dyz_S_C1001000003+CDY*I_ERI_D2x_S_Pz_S_C1001000003;
  Double I_ERI_Dxy_S_Pz_Py_C1001000003 = I_ERI_Dxy_S_Dyz_S_C1001000003+CDY*I_ERI_Dxy_S_Pz_S_C1001000003;
  Double I_ERI_Dxz_S_Pz_Py_C1001000003 = I_ERI_Dxz_S_Dyz_S_C1001000003+CDY*I_ERI_Dxz_S_Pz_S_C1001000003;
  Double I_ERI_D2y_S_Pz_Py_C1001000003 = I_ERI_D2y_S_Dyz_S_C1001000003+CDY*I_ERI_D2y_S_Pz_S_C1001000003;
  Double I_ERI_Dyz_S_Pz_Py_C1001000003 = I_ERI_Dyz_S_Dyz_S_C1001000003+CDY*I_ERI_Dyz_S_Pz_S_C1001000003;
  Double I_ERI_D2z_S_Pz_Py_C1001000003 = I_ERI_D2z_S_Dyz_S_C1001000003+CDY*I_ERI_D2z_S_Pz_S_C1001000003;
  Double I_ERI_D2x_S_Px_Pz_C1001000003 = I_ERI_D2x_S_Dxz_S_C1001000003+CDZ*I_ERI_D2x_S_Px_S_C1001000003;
  Double I_ERI_Dxy_S_Px_Pz_C1001000003 = I_ERI_Dxy_S_Dxz_S_C1001000003+CDZ*I_ERI_Dxy_S_Px_S_C1001000003;
  Double I_ERI_Dxz_S_Px_Pz_C1001000003 = I_ERI_Dxz_S_Dxz_S_C1001000003+CDZ*I_ERI_Dxz_S_Px_S_C1001000003;
  Double I_ERI_D2y_S_Px_Pz_C1001000003 = I_ERI_D2y_S_Dxz_S_C1001000003+CDZ*I_ERI_D2y_S_Px_S_C1001000003;
  Double I_ERI_Dyz_S_Px_Pz_C1001000003 = I_ERI_Dyz_S_Dxz_S_C1001000003+CDZ*I_ERI_Dyz_S_Px_S_C1001000003;
  Double I_ERI_D2z_S_Px_Pz_C1001000003 = I_ERI_D2z_S_Dxz_S_C1001000003+CDZ*I_ERI_D2z_S_Px_S_C1001000003;
  Double I_ERI_D2x_S_Py_Pz_C1001000003 = I_ERI_D2x_S_Dyz_S_C1001000003+CDZ*I_ERI_D2x_S_Py_S_C1001000003;
  Double I_ERI_Dxy_S_Py_Pz_C1001000003 = I_ERI_Dxy_S_Dyz_S_C1001000003+CDZ*I_ERI_Dxy_S_Py_S_C1001000003;
  Double I_ERI_Dxz_S_Py_Pz_C1001000003 = I_ERI_Dxz_S_Dyz_S_C1001000003+CDZ*I_ERI_Dxz_S_Py_S_C1001000003;
  Double I_ERI_D2y_S_Py_Pz_C1001000003 = I_ERI_D2y_S_Dyz_S_C1001000003+CDZ*I_ERI_D2y_S_Py_S_C1001000003;
  Double I_ERI_Dyz_S_Py_Pz_C1001000003 = I_ERI_Dyz_S_Dyz_S_C1001000003+CDZ*I_ERI_Dyz_S_Py_S_C1001000003;
  Double I_ERI_D2z_S_Py_Pz_C1001000003 = I_ERI_D2z_S_Dyz_S_C1001000003+CDZ*I_ERI_D2z_S_Py_S_C1001000003;
  Double I_ERI_D2x_S_Pz_Pz_C1001000003 = I_ERI_D2x_S_D2z_S_C1001000003+CDZ*I_ERI_D2x_S_Pz_S_C1001000003;
  Double I_ERI_Dxy_S_Pz_Pz_C1001000003 = I_ERI_Dxy_S_D2z_S_C1001000003+CDZ*I_ERI_Dxy_S_Pz_S_C1001000003;
  Double I_ERI_Dxz_S_Pz_Pz_C1001000003 = I_ERI_Dxz_S_D2z_S_C1001000003+CDZ*I_ERI_Dxz_S_Pz_S_C1001000003;
  Double I_ERI_D2y_S_Pz_Pz_C1001000003 = I_ERI_D2y_S_D2z_S_C1001000003+CDZ*I_ERI_D2y_S_Pz_S_C1001000003;
  Double I_ERI_Dyz_S_Pz_Pz_C1001000003 = I_ERI_Dyz_S_D2z_S_C1001000003+CDZ*I_ERI_Dyz_S_Pz_S_C1001000003;
  Double I_ERI_D2z_S_Pz_Pz_C1001000003 = I_ERI_D2z_S_D2z_S_C1001000003+CDZ*I_ERI_D2z_S_Pz_S_C1001000003;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_P_C1001000003_a
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_C1001000003_a
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1001000003_a
   ************************************************************/
  Double I_ERI_G4x_S_Px_Px_C1001000003_a = I_ERI_G4x_S_D2x_S_C1001000003_a+CDX*I_ERI_G4x_S_Px_S_C1001000003_a;
  Double I_ERI_G3xy_S_Px_Px_C1001000003_a = I_ERI_G3xy_S_D2x_S_C1001000003_a+CDX*I_ERI_G3xy_S_Px_S_C1001000003_a;
  Double I_ERI_G3xz_S_Px_Px_C1001000003_a = I_ERI_G3xz_S_D2x_S_C1001000003_a+CDX*I_ERI_G3xz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Px_Px_C1001000003_a = I_ERI_G2x2y_S_D2x_S_C1001000003_a+CDX*I_ERI_G2x2y_S_Px_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Px_Px_C1001000003_a = I_ERI_G2xyz_S_D2x_S_C1001000003_a+CDX*I_ERI_G2xyz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Px_Px_C1001000003_a = I_ERI_G2x2z_S_D2x_S_C1001000003_a+CDX*I_ERI_G2x2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Px_Px_C1001000003_a = I_ERI_Gx3y_S_D2x_S_C1001000003_a+CDX*I_ERI_Gx3y_S_Px_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Px_Px_C1001000003_a = I_ERI_Gx2yz_S_D2x_S_C1001000003_a+CDX*I_ERI_Gx2yz_S_Px_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Px_Px_C1001000003_a = I_ERI_Gxy2z_S_D2x_S_C1001000003_a+CDX*I_ERI_Gxy2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Px_Px_C1001000003_a = I_ERI_Gx3z_S_D2x_S_C1001000003_a+CDX*I_ERI_Gx3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4y_S_Px_Px_C1001000003_a = I_ERI_G4y_S_D2x_S_C1001000003_a+CDX*I_ERI_G4y_S_Px_S_C1001000003_a;
  Double I_ERI_G3yz_S_Px_Px_C1001000003_a = I_ERI_G3yz_S_D2x_S_C1001000003_a+CDX*I_ERI_G3yz_S_Px_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Px_Px_C1001000003_a = I_ERI_G2y2z_S_D2x_S_C1001000003_a+CDX*I_ERI_G2y2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Px_Px_C1001000003_a = I_ERI_Gy3z_S_D2x_S_C1001000003_a+CDX*I_ERI_Gy3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4z_S_Px_Px_C1001000003_a = I_ERI_G4z_S_D2x_S_C1001000003_a+CDX*I_ERI_G4z_S_Px_S_C1001000003_a;
  Double I_ERI_G4x_S_Py_Px_C1001000003_a = I_ERI_G4x_S_Dxy_S_C1001000003_a+CDX*I_ERI_G4x_S_Py_S_C1001000003_a;
  Double I_ERI_G3xy_S_Py_Px_C1001000003_a = I_ERI_G3xy_S_Dxy_S_C1001000003_a+CDX*I_ERI_G3xy_S_Py_S_C1001000003_a;
  Double I_ERI_G3xz_S_Py_Px_C1001000003_a = I_ERI_G3xz_S_Dxy_S_C1001000003_a+CDX*I_ERI_G3xz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Py_Px_C1001000003_a = I_ERI_G2x2y_S_Dxy_S_C1001000003_a+CDX*I_ERI_G2x2y_S_Py_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Py_Px_C1001000003_a = I_ERI_G2xyz_S_Dxy_S_C1001000003_a+CDX*I_ERI_G2xyz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Py_Px_C1001000003_a = I_ERI_G2x2z_S_Dxy_S_C1001000003_a+CDX*I_ERI_G2x2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Py_Px_C1001000003_a = I_ERI_Gx3y_S_Dxy_S_C1001000003_a+CDX*I_ERI_Gx3y_S_Py_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Py_Px_C1001000003_a = I_ERI_Gx2yz_S_Dxy_S_C1001000003_a+CDX*I_ERI_Gx2yz_S_Py_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Py_Px_C1001000003_a = I_ERI_Gxy2z_S_Dxy_S_C1001000003_a+CDX*I_ERI_Gxy2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Py_Px_C1001000003_a = I_ERI_Gx3z_S_Dxy_S_C1001000003_a+CDX*I_ERI_Gx3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4y_S_Py_Px_C1001000003_a = I_ERI_G4y_S_Dxy_S_C1001000003_a+CDX*I_ERI_G4y_S_Py_S_C1001000003_a;
  Double I_ERI_G3yz_S_Py_Px_C1001000003_a = I_ERI_G3yz_S_Dxy_S_C1001000003_a+CDX*I_ERI_G3yz_S_Py_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Py_Px_C1001000003_a = I_ERI_G2y2z_S_Dxy_S_C1001000003_a+CDX*I_ERI_G2y2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Py_Px_C1001000003_a = I_ERI_Gy3z_S_Dxy_S_C1001000003_a+CDX*I_ERI_Gy3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4z_S_Py_Px_C1001000003_a = I_ERI_G4z_S_Dxy_S_C1001000003_a+CDX*I_ERI_G4z_S_Py_S_C1001000003_a;
  Double I_ERI_G4x_S_Pz_Px_C1001000003_a = I_ERI_G4x_S_Dxz_S_C1001000003_a+CDX*I_ERI_G4x_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xy_S_Pz_Px_C1001000003_a = I_ERI_G3xy_S_Dxz_S_C1001000003_a+CDX*I_ERI_G3xy_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xz_S_Pz_Px_C1001000003_a = I_ERI_G3xz_S_Dxz_S_C1001000003_a+CDX*I_ERI_G3xz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Pz_Px_C1001000003_a = I_ERI_G2x2y_S_Dxz_S_C1001000003_a+CDX*I_ERI_G2x2y_S_Pz_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Pz_Px_C1001000003_a = I_ERI_G2xyz_S_Dxz_S_C1001000003_a+CDX*I_ERI_G2xyz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Pz_Px_C1001000003_a = I_ERI_G2x2z_S_Dxz_S_C1001000003_a+CDX*I_ERI_G2x2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Pz_Px_C1001000003_a = I_ERI_Gx3y_S_Dxz_S_C1001000003_a+CDX*I_ERI_Gx3y_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Pz_Px_C1001000003_a = I_ERI_Gx2yz_S_Dxz_S_C1001000003_a+CDX*I_ERI_Gx2yz_S_Pz_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Pz_Px_C1001000003_a = I_ERI_Gxy2z_S_Dxz_S_C1001000003_a+CDX*I_ERI_Gxy2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Pz_Px_C1001000003_a = I_ERI_Gx3z_S_Dxz_S_C1001000003_a+CDX*I_ERI_Gx3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4y_S_Pz_Px_C1001000003_a = I_ERI_G4y_S_Dxz_S_C1001000003_a+CDX*I_ERI_G4y_S_Pz_S_C1001000003_a;
  Double I_ERI_G3yz_S_Pz_Px_C1001000003_a = I_ERI_G3yz_S_Dxz_S_C1001000003_a+CDX*I_ERI_G3yz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Pz_Px_C1001000003_a = I_ERI_G2y2z_S_Dxz_S_C1001000003_a+CDX*I_ERI_G2y2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Pz_Px_C1001000003_a = I_ERI_Gy3z_S_Dxz_S_C1001000003_a+CDX*I_ERI_Gy3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4z_S_Pz_Px_C1001000003_a = I_ERI_G4z_S_Dxz_S_C1001000003_a+CDX*I_ERI_G4z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4x_S_Px_Py_C1001000003_a = I_ERI_G4x_S_Dxy_S_C1001000003_a+CDY*I_ERI_G4x_S_Px_S_C1001000003_a;
  Double I_ERI_G3xy_S_Px_Py_C1001000003_a = I_ERI_G3xy_S_Dxy_S_C1001000003_a+CDY*I_ERI_G3xy_S_Px_S_C1001000003_a;
  Double I_ERI_G3xz_S_Px_Py_C1001000003_a = I_ERI_G3xz_S_Dxy_S_C1001000003_a+CDY*I_ERI_G3xz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Px_Py_C1001000003_a = I_ERI_G2x2y_S_Dxy_S_C1001000003_a+CDY*I_ERI_G2x2y_S_Px_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Px_Py_C1001000003_a = I_ERI_G2xyz_S_Dxy_S_C1001000003_a+CDY*I_ERI_G2xyz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Px_Py_C1001000003_a = I_ERI_G2x2z_S_Dxy_S_C1001000003_a+CDY*I_ERI_G2x2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Px_Py_C1001000003_a = I_ERI_Gx3y_S_Dxy_S_C1001000003_a+CDY*I_ERI_Gx3y_S_Px_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Px_Py_C1001000003_a = I_ERI_Gx2yz_S_Dxy_S_C1001000003_a+CDY*I_ERI_Gx2yz_S_Px_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Px_Py_C1001000003_a = I_ERI_Gxy2z_S_Dxy_S_C1001000003_a+CDY*I_ERI_Gxy2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Px_Py_C1001000003_a = I_ERI_Gx3z_S_Dxy_S_C1001000003_a+CDY*I_ERI_Gx3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4y_S_Px_Py_C1001000003_a = I_ERI_G4y_S_Dxy_S_C1001000003_a+CDY*I_ERI_G4y_S_Px_S_C1001000003_a;
  Double I_ERI_G3yz_S_Px_Py_C1001000003_a = I_ERI_G3yz_S_Dxy_S_C1001000003_a+CDY*I_ERI_G3yz_S_Px_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Px_Py_C1001000003_a = I_ERI_G2y2z_S_Dxy_S_C1001000003_a+CDY*I_ERI_G2y2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Px_Py_C1001000003_a = I_ERI_Gy3z_S_Dxy_S_C1001000003_a+CDY*I_ERI_Gy3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4z_S_Px_Py_C1001000003_a = I_ERI_G4z_S_Dxy_S_C1001000003_a+CDY*I_ERI_G4z_S_Px_S_C1001000003_a;
  Double I_ERI_G4x_S_Py_Py_C1001000003_a = I_ERI_G4x_S_D2y_S_C1001000003_a+CDY*I_ERI_G4x_S_Py_S_C1001000003_a;
  Double I_ERI_G3xy_S_Py_Py_C1001000003_a = I_ERI_G3xy_S_D2y_S_C1001000003_a+CDY*I_ERI_G3xy_S_Py_S_C1001000003_a;
  Double I_ERI_G3xz_S_Py_Py_C1001000003_a = I_ERI_G3xz_S_D2y_S_C1001000003_a+CDY*I_ERI_G3xz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Py_Py_C1001000003_a = I_ERI_G2x2y_S_D2y_S_C1001000003_a+CDY*I_ERI_G2x2y_S_Py_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Py_Py_C1001000003_a = I_ERI_G2xyz_S_D2y_S_C1001000003_a+CDY*I_ERI_G2xyz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Py_Py_C1001000003_a = I_ERI_G2x2z_S_D2y_S_C1001000003_a+CDY*I_ERI_G2x2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Py_Py_C1001000003_a = I_ERI_Gx3y_S_D2y_S_C1001000003_a+CDY*I_ERI_Gx3y_S_Py_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Py_Py_C1001000003_a = I_ERI_Gx2yz_S_D2y_S_C1001000003_a+CDY*I_ERI_Gx2yz_S_Py_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Py_Py_C1001000003_a = I_ERI_Gxy2z_S_D2y_S_C1001000003_a+CDY*I_ERI_Gxy2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Py_Py_C1001000003_a = I_ERI_Gx3z_S_D2y_S_C1001000003_a+CDY*I_ERI_Gx3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4y_S_Py_Py_C1001000003_a = I_ERI_G4y_S_D2y_S_C1001000003_a+CDY*I_ERI_G4y_S_Py_S_C1001000003_a;
  Double I_ERI_G3yz_S_Py_Py_C1001000003_a = I_ERI_G3yz_S_D2y_S_C1001000003_a+CDY*I_ERI_G3yz_S_Py_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Py_Py_C1001000003_a = I_ERI_G2y2z_S_D2y_S_C1001000003_a+CDY*I_ERI_G2y2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Py_Py_C1001000003_a = I_ERI_Gy3z_S_D2y_S_C1001000003_a+CDY*I_ERI_Gy3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4z_S_Py_Py_C1001000003_a = I_ERI_G4z_S_D2y_S_C1001000003_a+CDY*I_ERI_G4z_S_Py_S_C1001000003_a;
  Double I_ERI_G4x_S_Pz_Py_C1001000003_a = I_ERI_G4x_S_Dyz_S_C1001000003_a+CDY*I_ERI_G4x_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xy_S_Pz_Py_C1001000003_a = I_ERI_G3xy_S_Dyz_S_C1001000003_a+CDY*I_ERI_G3xy_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xz_S_Pz_Py_C1001000003_a = I_ERI_G3xz_S_Dyz_S_C1001000003_a+CDY*I_ERI_G3xz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Pz_Py_C1001000003_a = I_ERI_G2x2y_S_Dyz_S_C1001000003_a+CDY*I_ERI_G2x2y_S_Pz_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Pz_Py_C1001000003_a = I_ERI_G2xyz_S_Dyz_S_C1001000003_a+CDY*I_ERI_G2xyz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Pz_Py_C1001000003_a = I_ERI_G2x2z_S_Dyz_S_C1001000003_a+CDY*I_ERI_G2x2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Pz_Py_C1001000003_a = I_ERI_Gx3y_S_Dyz_S_C1001000003_a+CDY*I_ERI_Gx3y_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Pz_Py_C1001000003_a = I_ERI_Gx2yz_S_Dyz_S_C1001000003_a+CDY*I_ERI_Gx2yz_S_Pz_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Pz_Py_C1001000003_a = I_ERI_Gxy2z_S_Dyz_S_C1001000003_a+CDY*I_ERI_Gxy2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Pz_Py_C1001000003_a = I_ERI_Gx3z_S_Dyz_S_C1001000003_a+CDY*I_ERI_Gx3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4y_S_Pz_Py_C1001000003_a = I_ERI_G4y_S_Dyz_S_C1001000003_a+CDY*I_ERI_G4y_S_Pz_S_C1001000003_a;
  Double I_ERI_G3yz_S_Pz_Py_C1001000003_a = I_ERI_G3yz_S_Dyz_S_C1001000003_a+CDY*I_ERI_G3yz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Pz_Py_C1001000003_a = I_ERI_G2y2z_S_Dyz_S_C1001000003_a+CDY*I_ERI_G2y2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Pz_Py_C1001000003_a = I_ERI_Gy3z_S_Dyz_S_C1001000003_a+CDY*I_ERI_Gy3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4z_S_Pz_Py_C1001000003_a = I_ERI_G4z_S_Dyz_S_C1001000003_a+CDY*I_ERI_G4z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4x_S_Px_Pz_C1001000003_a = I_ERI_G4x_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G4x_S_Px_S_C1001000003_a;
  Double I_ERI_G3xy_S_Px_Pz_C1001000003_a = I_ERI_G3xy_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G3xy_S_Px_S_C1001000003_a;
  Double I_ERI_G3xz_S_Px_Pz_C1001000003_a = I_ERI_G3xz_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G3xz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Px_Pz_C1001000003_a = I_ERI_G2x2y_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G2x2y_S_Px_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Px_Pz_C1001000003_a = I_ERI_G2xyz_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G2xyz_S_Px_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Px_Pz_C1001000003_a = I_ERI_G2x2z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G2x2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Px_Pz_C1001000003_a = I_ERI_Gx3y_S_Dxz_S_C1001000003_a+CDZ*I_ERI_Gx3y_S_Px_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Px_Pz_C1001000003_a = I_ERI_Gx2yz_S_Dxz_S_C1001000003_a+CDZ*I_ERI_Gx2yz_S_Px_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Px_Pz_C1001000003_a = I_ERI_Gxy2z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_Gxy2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Px_Pz_C1001000003_a = I_ERI_Gx3z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_Gx3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4y_S_Px_Pz_C1001000003_a = I_ERI_G4y_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G4y_S_Px_S_C1001000003_a;
  Double I_ERI_G3yz_S_Px_Pz_C1001000003_a = I_ERI_G3yz_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G3yz_S_Px_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Px_Pz_C1001000003_a = I_ERI_G2y2z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G2y2z_S_Px_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Px_Pz_C1001000003_a = I_ERI_Gy3z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_Gy3z_S_Px_S_C1001000003_a;
  Double I_ERI_G4z_S_Px_Pz_C1001000003_a = I_ERI_G4z_S_Dxz_S_C1001000003_a+CDZ*I_ERI_G4z_S_Px_S_C1001000003_a;
  Double I_ERI_G4x_S_Py_Pz_C1001000003_a = I_ERI_G4x_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G4x_S_Py_S_C1001000003_a;
  Double I_ERI_G3xy_S_Py_Pz_C1001000003_a = I_ERI_G3xy_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G3xy_S_Py_S_C1001000003_a;
  Double I_ERI_G3xz_S_Py_Pz_C1001000003_a = I_ERI_G3xz_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G3xz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Py_Pz_C1001000003_a = I_ERI_G2x2y_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G2x2y_S_Py_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Py_Pz_C1001000003_a = I_ERI_G2xyz_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G2xyz_S_Py_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Py_Pz_C1001000003_a = I_ERI_G2x2z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G2x2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Py_Pz_C1001000003_a = I_ERI_Gx3y_S_Dyz_S_C1001000003_a+CDZ*I_ERI_Gx3y_S_Py_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Py_Pz_C1001000003_a = I_ERI_Gx2yz_S_Dyz_S_C1001000003_a+CDZ*I_ERI_Gx2yz_S_Py_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Py_Pz_C1001000003_a = I_ERI_Gxy2z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_Gxy2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Py_Pz_C1001000003_a = I_ERI_Gx3z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_Gx3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4y_S_Py_Pz_C1001000003_a = I_ERI_G4y_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G4y_S_Py_S_C1001000003_a;
  Double I_ERI_G3yz_S_Py_Pz_C1001000003_a = I_ERI_G3yz_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G3yz_S_Py_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Py_Pz_C1001000003_a = I_ERI_G2y2z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G2y2z_S_Py_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Py_Pz_C1001000003_a = I_ERI_Gy3z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_Gy3z_S_Py_S_C1001000003_a;
  Double I_ERI_G4z_S_Py_Pz_C1001000003_a = I_ERI_G4z_S_Dyz_S_C1001000003_a+CDZ*I_ERI_G4z_S_Py_S_C1001000003_a;
  Double I_ERI_G4x_S_Pz_Pz_C1001000003_a = I_ERI_G4x_S_D2z_S_C1001000003_a+CDZ*I_ERI_G4x_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xy_S_Pz_Pz_C1001000003_a = I_ERI_G3xy_S_D2z_S_C1001000003_a+CDZ*I_ERI_G3xy_S_Pz_S_C1001000003_a;
  Double I_ERI_G3xz_S_Pz_Pz_C1001000003_a = I_ERI_G3xz_S_D2z_S_C1001000003_a+CDZ*I_ERI_G3xz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2y_S_Pz_Pz_C1001000003_a = I_ERI_G2x2y_S_D2z_S_C1001000003_a+CDZ*I_ERI_G2x2y_S_Pz_S_C1001000003_a;
  Double I_ERI_G2xyz_S_Pz_Pz_C1001000003_a = I_ERI_G2xyz_S_D2z_S_C1001000003_a+CDZ*I_ERI_G2xyz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2x2z_S_Pz_Pz_C1001000003_a = I_ERI_G2x2z_S_D2z_S_C1001000003_a+CDZ*I_ERI_G2x2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3y_S_Pz_Pz_C1001000003_a = I_ERI_Gx3y_S_D2z_S_C1001000003_a+CDZ*I_ERI_Gx3y_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx2yz_S_Pz_Pz_C1001000003_a = I_ERI_Gx2yz_S_D2z_S_C1001000003_a+CDZ*I_ERI_Gx2yz_S_Pz_S_C1001000003_a;
  Double I_ERI_Gxy2z_S_Pz_Pz_C1001000003_a = I_ERI_Gxy2z_S_D2z_S_C1001000003_a+CDZ*I_ERI_Gxy2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gx3z_S_Pz_Pz_C1001000003_a = I_ERI_Gx3z_S_D2z_S_C1001000003_a+CDZ*I_ERI_Gx3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4y_S_Pz_Pz_C1001000003_a = I_ERI_G4y_S_D2z_S_C1001000003_a+CDZ*I_ERI_G4y_S_Pz_S_C1001000003_a;
  Double I_ERI_G3yz_S_Pz_Pz_C1001000003_a = I_ERI_G3yz_S_D2z_S_C1001000003_a+CDZ*I_ERI_G3yz_S_Pz_S_C1001000003_a;
  Double I_ERI_G2y2z_S_Pz_Pz_C1001000003_a = I_ERI_G2y2z_S_D2z_S_C1001000003_a+CDZ*I_ERI_G2y2z_S_Pz_S_C1001000003_a;
  Double I_ERI_Gy3z_S_Pz_Pz_C1001000003_a = I_ERI_Gy3z_S_D2z_S_C1001000003_a+CDZ*I_ERI_Gy3z_S_Pz_S_C1001000003_a;
  Double I_ERI_G4z_S_Pz_Pz_C1001000003_a = I_ERI_G4z_S_D2z_S_C1001000003_a+CDZ*I_ERI_G4z_S_Pz_S_C1001000003_a;

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
   * shell quartet name: SQ_ERI_F_S_P_P_C1000000003_c
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000000003_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000000003_c
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_C1000000003_c = I_ERI_F3x_S_D2x_S_C1000000003_c+CDX*I_ERI_F3x_S_Px_S_C1000000003_c;
  Double I_ERI_F2xy_S_Px_Px_C1000000003_c = I_ERI_F2xy_S_D2x_S_C1000000003_c+CDX*I_ERI_F2xy_S_Px_S_C1000000003_c;
  Double I_ERI_F2xz_S_Px_Px_C1000000003_c = I_ERI_F2xz_S_D2x_S_C1000000003_c+CDX*I_ERI_F2xz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Px_Px_C1000000003_c = I_ERI_Fx2y_S_D2x_S_C1000000003_c+CDX*I_ERI_Fx2y_S_Px_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Px_Px_C1000000003_c = I_ERI_Fxyz_S_D2x_S_C1000000003_c+CDX*I_ERI_Fxyz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Px_Px_C1000000003_c = I_ERI_Fx2z_S_D2x_S_C1000000003_c+CDX*I_ERI_Fx2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3y_S_Px_Px_C1000000003_c = I_ERI_F3y_S_D2x_S_C1000000003_c+CDX*I_ERI_F3y_S_Px_S_C1000000003_c;
  Double I_ERI_F2yz_S_Px_Px_C1000000003_c = I_ERI_F2yz_S_D2x_S_C1000000003_c+CDX*I_ERI_F2yz_S_Px_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Px_Px_C1000000003_c = I_ERI_Fy2z_S_D2x_S_C1000000003_c+CDX*I_ERI_Fy2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3z_S_Px_Px_C1000000003_c = I_ERI_F3z_S_D2x_S_C1000000003_c+CDX*I_ERI_F3z_S_Px_S_C1000000003_c;
  Double I_ERI_F3x_S_Py_Px_C1000000003_c = I_ERI_F3x_S_Dxy_S_C1000000003_c+CDX*I_ERI_F3x_S_Py_S_C1000000003_c;
  Double I_ERI_F2xy_S_Py_Px_C1000000003_c = I_ERI_F2xy_S_Dxy_S_C1000000003_c+CDX*I_ERI_F2xy_S_Py_S_C1000000003_c;
  Double I_ERI_F2xz_S_Py_Px_C1000000003_c = I_ERI_F2xz_S_Dxy_S_C1000000003_c+CDX*I_ERI_F2xz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Py_Px_C1000000003_c = I_ERI_Fx2y_S_Dxy_S_C1000000003_c+CDX*I_ERI_Fx2y_S_Py_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Py_Px_C1000000003_c = I_ERI_Fxyz_S_Dxy_S_C1000000003_c+CDX*I_ERI_Fxyz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Py_Px_C1000000003_c = I_ERI_Fx2z_S_Dxy_S_C1000000003_c+CDX*I_ERI_Fx2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3y_S_Py_Px_C1000000003_c = I_ERI_F3y_S_Dxy_S_C1000000003_c+CDX*I_ERI_F3y_S_Py_S_C1000000003_c;
  Double I_ERI_F2yz_S_Py_Px_C1000000003_c = I_ERI_F2yz_S_Dxy_S_C1000000003_c+CDX*I_ERI_F2yz_S_Py_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Py_Px_C1000000003_c = I_ERI_Fy2z_S_Dxy_S_C1000000003_c+CDX*I_ERI_Fy2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3z_S_Py_Px_C1000000003_c = I_ERI_F3z_S_Dxy_S_C1000000003_c+CDX*I_ERI_F3z_S_Py_S_C1000000003_c;
  Double I_ERI_F3x_S_Pz_Px_C1000000003_c = I_ERI_F3x_S_Dxz_S_C1000000003_c+CDX*I_ERI_F3x_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xy_S_Pz_Px_C1000000003_c = I_ERI_F2xy_S_Dxz_S_C1000000003_c+CDX*I_ERI_F2xy_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xz_S_Pz_Px_C1000000003_c = I_ERI_F2xz_S_Dxz_S_C1000000003_c+CDX*I_ERI_F2xz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Pz_Px_C1000000003_c = I_ERI_Fx2y_S_Dxz_S_C1000000003_c+CDX*I_ERI_Fx2y_S_Pz_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Pz_Px_C1000000003_c = I_ERI_Fxyz_S_Dxz_S_C1000000003_c+CDX*I_ERI_Fxyz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Pz_Px_C1000000003_c = I_ERI_Fx2z_S_Dxz_S_C1000000003_c+CDX*I_ERI_Fx2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3y_S_Pz_Px_C1000000003_c = I_ERI_F3y_S_Dxz_S_C1000000003_c+CDX*I_ERI_F3y_S_Pz_S_C1000000003_c;
  Double I_ERI_F2yz_S_Pz_Px_C1000000003_c = I_ERI_F2yz_S_Dxz_S_C1000000003_c+CDX*I_ERI_F2yz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Pz_Px_C1000000003_c = I_ERI_Fy2z_S_Dxz_S_C1000000003_c+CDX*I_ERI_Fy2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3z_S_Pz_Px_C1000000003_c = I_ERI_F3z_S_Dxz_S_C1000000003_c+CDX*I_ERI_F3z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3x_S_Px_Py_C1000000003_c = I_ERI_F3x_S_Dxy_S_C1000000003_c+CDY*I_ERI_F3x_S_Px_S_C1000000003_c;
  Double I_ERI_F2xy_S_Px_Py_C1000000003_c = I_ERI_F2xy_S_Dxy_S_C1000000003_c+CDY*I_ERI_F2xy_S_Px_S_C1000000003_c;
  Double I_ERI_F2xz_S_Px_Py_C1000000003_c = I_ERI_F2xz_S_Dxy_S_C1000000003_c+CDY*I_ERI_F2xz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Px_Py_C1000000003_c = I_ERI_Fx2y_S_Dxy_S_C1000000003_c+CDY*I_ERI_Fx2y_S_Px_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Px_Py_C1000000003_c = I_ERI_Fxyz_S_Dxy_S_C1000000003_c+CDY*I_ERI_Fxyz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Px_Py_C1000000003_c = I_ERI_Fx2z_S_Dxy_S_C1000000003_c+CDY*I_ERI_Fx2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3y_S_Px_Py_C1000000003_c = I_ERI_F3y_S_Dxy_S_C1000000003_c+CDY*I_ERI_F3y_S_Px_S_C1000000003_c;
  Double I_ERI_F2yz_S_Px_Py_C1000000003_c = I_ERI_F2yz_S_Dxy_S_C1000000003_c+CDY*I_ERI_F2yz_S_Px_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Px_Py_C1000000003_c = I_ERI_Fy2z_S_Dxy_S_C1000000003_c+CDY*I_ERI_Fy2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3z_S_Px_Py_C1000000003_c = I_ERI_F3z_S_Dxy_S_C1000000003_c+CDY*I_ERI_F3z_S_Px_S_C1000000003_c;
  Double I_ERI_F3x_S_Py_Py_C1000000003_c = I_ERI_F3x_S_D2y_S_C1000000003_c+CDY*I_ERI_F3x_S_Py_S_C1000000003_c;
  Double I_ERI_F2xy_S_Py_Py_C1000000003_c = I_ERI_F2xy_S_D2y_S_C1000000003_c+CDY*I_ERI_F2xy_S_Py_S_C1000000003_c;
  Double I_ERI_F2xz_S_Py_Py_C1000000003_c = I_ERI_F2xz_S_D2y_S_C1000000003_c+CDY*I_ERI_F2xz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Py_Py_C1000000003_c = I_ERI_Fx2y_S_D2y_S_C1000000003_c+CDY*I_ERI_Fx2y_S_Py_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Py_Py_C1000000003_c = I_ERI_Fxyz_S_D2y_S_C1000000003_c+CDY*I_ERI_Fxyz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Py_Py_C1000000003_c = I_ERI_Fx2z_S_D2y_S_C1000000003_c+CDY*I_ERI_Fx2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3y_S_Py_Py_C1000000003_c = I_ERI_F3y_S_D2y_S_C1000000003_c+CDY*I_ERI_F3y_S_Py_S_C1000000003_c;
  Double I_ERI_F2yz_S_Py_Py_C1000000003_c = I_ERI_F2yz_S_D2y_S_C1000000003_c+CDY*I_ERI_F2yz_S_Py_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Py_Py_C1000000003_c = I_ERI_Fy2z_S_D2y_S_C1000000003_c+CDY*I_ERI_Fy2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3z_S_Py_Py_C1000000003_c = I_ERI_F3z_S_D2y_S_C1000000003_c+CDY*I_ERI_F3z_S_Py_S_C1000000003_c;
  Double I_ERI_F3x_S_Pz_Py_C1000000003_c = I_ERI_F3x_S_Dyz_S_C1000000003_c+CDY*I_ERI_F3x_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xy_S_Pz_Py_C1000000003_c = I_ERI_F2xy_S_Dyz_S_C1000000003_c+CDY*I_ERI_F2xy_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xz_S_Pz_Py_C1000000003_c = I_ERI_F2xz_S_Dyz_S_C1000000003_c+CDY*I_ERI_F2xz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Pz_Py_C1000000003_c = I_ERI_Fx2y_S_Dyz_S_C1000000003_c+CDY*I_ERI_Fx2y_S_Pz_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Pz_Py_C1000000003_c = I_ERI_Fxyz_S_Dyz_S_C1000000003_c+CDY*I_ERI_Fxyz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Pz_Py_C1000000003_c = I_ERI_Fx2z_S_Dyz_S_C1000000003_c+CDY*I_ERI_Fx2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3y_S_Pz_Py_C1000000003_c = I_ERI_F3y_S_Dyz_S_C1000000003_c+CDY*I_ERI_F3y_S_Pz_S_C1000000003_c;
  Double I_ERI_F2yz_S_Pz_Py_C1000000003_c = I_ERI_F2yz_S_Dyz_S_C1000000003_c+CDY*I_ERI_F2yz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Pz_Py_C1000000003_c = I_ERI_Fy2z_S_Dyz_S_C1000000003_c+CDY*I_ERI_Fy2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3z_S_Pz_Py_C1000000003_c = I_ERI_F3z_S_Dyz_S_C1000000003_c+CDY*I_ERI_F3z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3x_S_Px_Pz_C1000000003_c = I_ERI_F3x_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F3x_S_Px_S_C1000000003_c;
  Double I_ERI_F2xy_S_Px_Pz_C1000000003_c = I_ERI_F2xy_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F2xy_S_Px_S_C1000000003_c;
  Double I_ERI_F2xz_S_Px_Pz_C1000000003_c = I_ERI_F2xz_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F2xz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Px_Pz_C1000000003_c = I_ERI_Fx2y_S_Dxz_S_C1000000003_c+CDZ*I_ERI_Fx2y_S_Px_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Px_Pz_C1000000003_c = I_ERI_Fxyz_S_Dxz_S_C1000000003_c+CDZ*I_ERI_Fxyz_S_Px_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Px_Pz_C1000000003_c = I_ERI_Fx2z_S_Dxz_S_C1000000003_c+CDZ*I_ERI_Fx2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3y_S_Px_Pz_C1000000003_c = I_ERI_F3y_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F3y_S_Px_S_C1000000003_c;
  Double I_ERI_F2yz_S_Px_Pz_C1000000003_c = I_ERI_F2yz_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F2yz_S_Px_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Px_Pz_C1000000003_c = I_ERI_Fy2z_S_Dxz_S_C1000000003_c+CDZ*I_ERI_Fy2z_S_Px_S_C1000000003_c;
  Double I_ERI_F3z_S_Px_Pz_C1000000003_c = I_ERI_F3z_S_Dxz_S_C1000000003_c+CDZ*I_ERI_F3z_S_Px_S_C1000000003_c;
  Double I_ERI_F3x_S_Py_Pz_C1000000003_c = I_ERI_F3x_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F3x_S_Py_S_C1000000003_c;
  Double I_ERI_F2xy_S_Py_Pz_C1000000003_c = I_ERI_F2xy_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F2xy_S_Py_S_C1000000003_c;
  Double I_ERI_F2xz_S_Py_Pz_C1000000003_c = I_ERI_F2xz_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F2xz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Py_Pz_C1000000003_c = I_ERI_Fx2y_S_Dyz_S_C1000000003_c+CDZ*I_ERI_Fx2y_S_Py_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Py_Pz_C1000000003_c = I_ERI_Fxyz_S_Dyz_S_C1000000003_c+CDZ*I_ERI_Fxyz_S_Py_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Py_Pz_C1000000003_c = I_ERI_Fx2z_S_Dyz_S_C1000000003_c+CDZ*I_ERI_Fx2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3y_S_Py_Pz_C1000000003_c = I_ERI_F3y_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F3y_S_Py_S_C1000000003_c;
  Double I_ERI_F2yz_S_Py_Pz_C1000000003_c = I_ERI_F2yz_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F2yz_S_Py_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Py_Pz_C1000000003_c = I_ERI_Fy2z_S_Dyz_S_C1000000003_c+CDZ*I_ERI_Fy2z_S_Py_S_C1000000003_c;
  Double I_ERI_F3z_S_Py_Pz_C1000000003_c = I_ERI_F3z_S_Dyz_S_C1000000003_c+CDZ*I_ERI_F3z_S_Py_S_C1000000003_c;
  Double I_ERI_F3x_S_Pz_Pz_C1000000003_c = I_ERI_F3x_S_D2z_S_C1000000003_c+CDZ*I_ERI_F3x_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xy_S_Pz_Pz_C1000000003_c = I_ERI_F2xy_S_D2z_S_C1000000003_c+CDZ*I_ERI_F2xy_S_Pz_S_C1000000003_c;
  Double I_ERI_F2xz_S_Pz_Pz_C1000000003_c = I_ERI_F2xz_S_D2z_S_C1000000003_c+CDZ*I_ERI_F2xz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2y_S_Pz_Pz_C1000000003_c = I_ERI_Fx2y_S_D2z_S_C1000000003_c+CDZ*I_ERI_Fx2y_S_Pz_S_C1000000003_c;
  Double I_ERI_Fxyz_S_Pz_Pz_C1000000003_c = I_ERI_Fxyz_S_D2z_S_C1000000003_c+CDZ*I_ERI_Fxyz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fx2z_S_Pz_Pz_C1000000003_c = I_ERI_Fx2z_S_D2z_S_C1000000003_c+CDZ*I_ERI_Fx2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3y_S_Pz_Pz_C1000000003_c = I_ERI_F3y_S_D2z_S_C1000000003_c+CDZ*I_ERI_F3y_S_Pz_S_C1000000003_c;
  Double I_ERI_F2yz_S_Pz_Pz_C1000000003_c = I_ERI_F2yz_S_D2z_S_C1000000003_c+CDZ*I_ERI_F2yz_S_Pz_S_C1000000003_c;
  Double I_ERI_Fy2z_S_Pz_Pz_C1000000003_c = I_ERI_Fy2z_S_D2z_S_C1000000003_c+CDZ*I_ERI_Fy2z_S_Pz_S_C1000000003_c;
  Double I_ERI_F3z_S_Pz_Pz_C1000000003_c = I_ERI_F3z_S_D2z_S_C1000000003_c+CDZ*I_ERI_F3z_S_Pz_S_C1000000003_c;

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
   * declare the HRR2 result shell quartets in array form
   ************************************************************/

  /************************************************************
   * initilize the HRR steps : build the AB/CD variables
   ************************************************************/
  Double ABX = A[0] - B[0];
  Double ABY = A[1] - B[1];
  Double ABZ = A[2] - B[2];

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C3_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C3_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C3_b
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C3_b = I_ERI_G4x_S_S_S_C3_b+ABX*I_ERI_F3x_S_S_S_C3_b;
  Double I_ERI_F2xy_Px_S_S_C3_b = I_ERI_G3xy_S_S_S_C3_b+ABX*I_ERI_F2xy_S_S_S_C3_b;
  Double I_ERI_F2xz_Px_S_S_C3_b = I_ERI_G3xz_S_S_S_C3_b+ABX*I_ERI_F2xz_S_S_S_C3_b;
  Double I_ERI_Fx2y_Px_S_S_C3_b = I_ERI_G2x2y_S_S_S_C3_b+ABX*I_ERI_Fx2y_S_S_S_C3_b;
  Double I_ERI_Fxyz_Px_S_S_C3_b = I_ERI_G2xyz_S_S_S_C3_b+ABX*I_ERI_Fxyz_S_S_S_C3_b;
  Double I_ERI_Fx2z_Px_S_S_C3_b = I_ERI_G2x2z_S_S_S_C3_b+ABX*I_ERI_Fx2z_S_S_S_C3_b;
  Double I_ERI_F3y_Px_S_S_C3_b = I_ERI_Gx3y_S_S_S_C3_b+ABX*I_ERI_F3y_S_S_S_C3_b;
  Double I_ERI_F2yz_Px_S_S_C3_b = I_ERI_Gx2yz_S_S_S_C3_b+ABX*I_ERI_F2yz_S_S_S_C3_b;
  Double I_ERI_Fy2z_Px_S_S_C3_b = I_ERI_Gxy2z_S_S_S_C3_b+ABX*I_ERI_Fy2z_S_S_S_C3_b;
  Double I_ERI_F3z_Px_S_S_C3_b = I_ERI_Gx3z_S_S_S_C3_b+ABX*I_ERI_F3z_S_S_S_C3_b;
  Double I_ERI_F3x_Py_S_S_C3_b = I_ERI_G3xy_S_S_S_C3_b+ABY*I_ERI_F3x_S_S_S_C3_b;
  Double I_ERI_F2xy_Py_S_S_C3_b = I_ERI_G2x2y_S_S_S_C3_b+ABY*I_ERI_F2xy_S_S_S_C3_b;
  Double I_ERI_F2xz_Py_S_S_C3_b = I_ERI_G2xyz_S_S_S_C3_b+ABY*I_ERI_F2xz_S_S_S_C3_b;
  Double I_ERI_Fx2y_Py_S_S_C3_b = I_ERI_Gx3y_S_S_S_C3_b+ABY*I_ERI_Fx2y_S_S_S_C3_b;
  Double I_ERI_Fxyz_Py_S_S_C3_b = I_ERI_Gx2yz_S_S_S_C3_b+ABY*I_ERI_Fxyz_S_S_S_C3_b;
  Double I_ERI_Fx2z_Py_S_S_C3_b = I_ERI_Gxy2z_S_S_S_C3_b+ABY*I_ERI_Fx2z_S_S_S_C3_b;
  Double I_ERI_F3y_Py_S_S_C3_b = I_ERI_G4y_S_S_S_C3_b+ABY*I_ERI_F3y_S_S_S_C3_b;
  Double I_ERI_F2yz_Py_S_S_C3_b = I_ERI_G3yz_S_S_S_C3_b+ABY*I_ERI_F2yz_S_S_S_C3_b;
  Double I_ERI_Fy2z_Py_S_S_C3_b = I_ERI_G2y2z_S_S_S_C3_b+ABY*I_ERI_Fy2z_S_S_S_C3_b;
  Double I_ERI_F3z_Py_S_S_C3_b = I_ERI_Gy3z_S_S_S_C3_b+ABY*I_ERI_F3z_S_S_S_C3_b;
  Double I_ERI_F3x_Pz_S_S_C3_b = I_ERI_G3xz_S_S_S_C3_b+ABZ*I_ERI_F3x_S_S_S_C3_b;
  Double I_ERI_F2xy_Pz_S_S_C3_b = I_ERI_G2xyz_S_S_S_C3_b+ABZ*I_ERI_F2xy_S_S_S_C3_b;
  Double I_ERI_F2xz_Pz_S_S_C3_b = I_ERI_G2x2z_S_S_S_C3_b+ABZ*I_ERI_F2xz_S_S_S_C3_b;
  Double I_ERI_Fx2y_Pz_S_S_C3_b = I_ERI_Gx2yz_S_S_S_C3_b+ABZ*I_ERI_Fx2y_S_S_S_C3_b;
  Double I_ERI_Fxyz_Pz_S_S_C3_b = I_ERI_Gxy2z_S_S_S_C3_b+ABZ*I_ERI_Fxyz_S_S_S_C3_b;
  Double I_ERI_Fx2z_Pz_S_S_C3_b = I_ERI_Gx3z_S_S_S_C3_b+ABZ*I_ERI_Fx2z_S_S_S_C3_b;
  Double I_ERI_F3y_Pz_S_S_C3_b = I_ERI_G3yz_S_S_S_C3_b+ABZ*I_ERI_F3y_S_S_S_C3_b;
  Double I_ERI_F2yz_Pz_S_S_C3_b = I_ERI_G2y2z_S_S_S_C3_b+ABZ*I_ERI_F2yz_S_S_S_C3_b;
  Double I_ERI_Fy2z_Pz_S_S_C3_b = I_ERI_Gy3z_S_S_S_C3_b+ABZ*I_ERI_Fy2z_S_S_S_C3_b;
  Double I_ERI_F3z_Pz_S_S_C3_b = I_ERI_G4z_S_S_S_C3_b+ABZ*I_ERI_F3z_S_S_S_C3_b;

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
   * shell quartet name: SQ_ERI_F_P_S_P_C1000000003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_P_C1000000003_b
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1000000003_b
   ************************************************************/
  Double I_ERI_F3x_Px_S_Px_C1000000003_b = I_ERI_G4x_S_S_Px_C1000000003_b+ABX*I_ERI_F3x_S_S_Px_C1000000003_b;
  Double I_ERI_F2xy_Px_S_Px_C1000000003_b = I_ERI_G3xy_S_S_Px_C1000000003_b+ABX*I_ERI_F2xy_S_S_Px_C1000000003_b;
  Double I_ERI_F2xz_Px_S_Px_C1000000003_b = I_ERI_G3xz_S_S_Px_C1000000003_b+ABX*I_ERI_F2xz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2y_Px_S_Px_C1000000003_b = I_ERI_G2x2y_S_S_Px_C1000000003_b+ABX*I_ERI_Fx2y_S_S_Px_C1000000003_b;
  Double I_ERI_Fxyz_Px_S_Px_C1000000003_b = I_ERI_G2xyz_S_S_Px_C1000000003_b+ABX*I_ERI_Fxyz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2z_Px_S_Px_C1000000003_b = I_ERI_G2x2z_S_S_Px_C1000000003_b+ABX*I_ERI_Fx2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3y_Px_S_Px_C1000000003_b = I_ERI_Gx3y_S_S_Px_C1000000003_b+ABX*I_ERI_F3y_S_S_Px_C1000000003_b;
  Double I_ERI_F2yz_Px_S_Px_C1000000003_b = I_ERI_Gx2yz_S_S_Px_C1000000003_b+ABX*I_ERI_F2yz_S_S_Px_C1000000003_b;
  Double I_ERI_Fy2z_Px_S_Px_C1000000003_b = I_ERI_Gxy2z_S_S_Px_C1000000003_b+ABX*I_ERI_Fy2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3z_Px_S_Px_C1000000003_b = I_ERI_Gx3z_S_S_Px_C1000000003_b+ABX*I_ERI_F3z_S_S_Px_C1000000003_b;
  Double I_ERI_F3x_Py_S_Px_C1000000003_b = I_ERI_G3xy_S_S_Px_C1000000003_b+ABY*I_ERI_F3x_S_S_Px_C1000000003_b;
  Double I_ERI_F2xy_Py_S_Px_C1000000003_b = I_ERI_G2x2y_S_S_Px_C1000000003_b+ABY*I_ERI_F2xy_S_S_Px_C1000000003_b;
  Double I_ERI_F2xz_Py_S_Px_C1000000003_b = I_ERI_G2xyz_S_S_Px_C1000000003_b+ABY*I_ERI_F2xz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2y_Py_S_Px_C1000000003_b = I_ERI_Gx3y_S_S_Px_C1000000003_b+ABY*I_ERI_Fx2y_S_S_Px_C1000000003_b;
  Double I_ERI_Fxyz_Py_S_Px_C1000000003_b = I_ERI_Gx2yz_S_S_Px_C1000000003_b+ABY*I_ERI_Fxyz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2z_Py_S_Px_C1000000003_b = I_ERI_Gxy2z_S_S_Px_C1000000003_b+ABY*I_ERI_Fx2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3y_Py_S_Px_C1000000003_b = I_ERI_G4y_S_S_Px_C1000000003_b+ABY*I_ERI_F3y_S_S_Px_C1000000003_b;
  Double I_ERI_F2yz_Py_S_Px_C1000000003_b = I_ERI_G3yz_S_S_Px_C1000000003_b+ABY*I_ERI_F2yz_S_S_Px_C1000000003_b;
  Double I_ERI_Fy2z_Py_S_Px_C1000000003_b = I_ERI_G2y2z_S_S_Px_C1000000003_b+ABY*I_ERI_Fy2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3z_Py_S_Px_C1000000003_b = I_ERI_Gy3z_S_S_Px_C1000000003_b+ABY*I_ERI_F3z_S_S_Px_C1000000003_b;
  Double I_ERI_F3x_Pz_S_Px_C1000000003_b = I_ERI_G3xz_S_S_Px_C1000000003_b+ABZ*I_ERI_F3x_S_S_Px_C1000000003_b;
  Double I_ERI_F2xy_Pz_S_Px_C1000000003_b = I_ERI_G2xyz_S_S_Px_C1000000003_b+ABZ*I_ERI_F2xy_S_S_Px_C1000000003_b;
  Double I_ERI_F2xz_Pz_S_Px_C1000000003_b = I_ERI_G2x2z_S_S_Px_C1000000003_b+ABZ*I_ERI_F2xz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2y_Pz_S_Px_C1000000003_b = I_ERI_Gx2yz_S_S_Px_C1000000003_b+ABZ*I_ERI_Fx2y_S_S_Px_C1000000003_b;
  Double I_ERI_Fxyz_Pz_S_Px_C1000000003_b = I_ERI_Gxy2z_S_S_Px_C1000000003_b+ABZ*I_ERI_Fxyz_S_S_Px_C1000000003_b;
  Double I_ERI_Fx2z_Pz_S_Px_C1000000003_b = I_ERI_Gx3z_S_S_Px_C1000000003_b+ABZ*I_ERI_Fx2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3y_Pz_S_Px_C1000000003_b = I_ERI_G3yz_S_S_Px_C1000000003_b+ABZ*I_ERI_F3y_S_S_Px_C1000000003_b;
  Double I_ERI_F2yz_Pz_S_Px_C1000000003_b = I_ERI_G2y2z_S_S_Px_C1000000003_b+ABZ*I_ERI_F2yz_S_S_Px_C1000000003_b;
  Double I_ERI_Fy2z_Pz_S_Px_C1000000003_b = I_ERI_Gy3z_S_S_Px_C1000000003_b+ABZ*I_ERI_Fy2z_S_S_Px_C1000000003_b;
  Double I_ERI_F3z_Pz_S_Px_C1000000003_b = I_ERI_G4z_S_S_Px_C1000000003_b+ABZ*I_ERI_F3z_S_S_Px_C1000000003_b;
  Double I_ERI_F3x_Px_S_Py_C1000000003_b = I_ERI_G4x_S_S_Py_C1000000003_b+ABX*I_ERI_F3x_S_S_Py_C1000000003_b;
  Double I_ERI_F2xy_Px_S_Py_C1000000003_b = I_ERI_G3xy_S_S_Py_C1000000003_b+ABX*I_ERI_F2xy_S_S_Py_C1000000003_b;
  Double I_ERI_F2xz_Px_S_Py_C1000000003_b = I_ERI_G3xz_S_S_Py_C1000000003_b+ABX*I_ERI_F2xz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2y_Px_S_Py_C1000000003_b = I_ERI_G2x2y_S_S_Py_C1000000003_b+ABX*I_ERI_Fx2y_S_S_Py_C1000000003_b;
  Double I_ERI_Fxyz_Px_S_Py_C1000000003_b = I_ERI_G2xyz_S_S_Py_C1000000003_b+ABX*I_ERI_Fxyz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2z_Px_S_Py_C1000000003_b = I_ERI_G2x2z_S_S_Py_C1000000003_b+ABX*I_ERI_Fx2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3y_Px_S_Py_C1000000003_b = I_ERI_Gx3y_S_S_Py_C1000000003_b+ABX*I_ERI_F3y_S_S_Py_C1000000003_b;
  Double I_ERI_F2yz_Px_S_Py_C1000000003_b = I_ERI_Gx2yz_S_S_Py_C1000000003_b+ABX*I_ERI_F2yz_S_S_Py_C1000000003_b;
  Double I_ERI_Fy2z_Px_S_Py_C1000000003_b = I_ERI_Gxy2z_S_S_Py_C1000000003_b+ABX*I_ERI_Fy2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3z_Px_S_Py_C1000000003_b = I_ERI_Gx3z_S_S_Py_C1000000003_b+ABX*I_ERI_F3z_S_S_Py_C1000000003_b;
  Double I_ERI_F3x_Py_S_Py_C1000000003_b = I_ERI_G3xy_S_S_Py_C1000000003_b+ABY*I_ERI_F3x_S_S_Py_C1000000003_b;
  Double I_ERI_F2xy_Py_S_Py_C1000000003_b = I_ERI_G2x2y_S_S_Py_C1000000003_b+ABY*I_ERI_F2xy_S_S_Py_C1000000003_b;
  Double I_ERI_F2xz_Py_S_Py_C1000000003_b = I_ERI_G2xyz_S_S_Py_C1000000003_b+ABY*I_ERI_F2xz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2y_Py_S_Py_C1000000003_b = I_ERI_Gx3y_S_S_Py_C1000000003_b+ABY*I_ERI_Fx2y_S_S_Py_C1000000003_b;
  Double I_ERI_Fxyz_Py_S_Py_C1000000003_b = I_ERI_Gx2yz_S_S_Py_C1000000003_b+ABY*I_ERI_Fxyz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2z_Py_S_Py_C1000000003_b = I_ERI_Gxy2z_S_S_Py_C1000000003_b+ABY*I_ERI_Fx2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3y_Py_S_Py_C1000000003_b = I_ERI_G4y_S_S_Py_C1000000003_b+ABY*I_ERI_F3y_S_S_Py_C1000000003_b;
  Double I_ERI_F2yz_Py_S_Py_C1000000003_b = I_ERI_G3yz_S_S_Py_C1000000003_b+ABY*I_ERI_F2yz_S_S_Py_C1000000003_b;
  Double I_ERI_Fy2z_Py_S_Py_C1000000003_b = I_ERI_G2y2z_S_S_Py_C1000000003_b+ABY*I_ERI_Fy2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3z_Py_S_Py_C1000000003_b = I_ERI_Gy3z_S_S_Py_C1000000003_b+ABY*I_ERI_F3z_S_S_Py_C1000000003_b;
  Double I_ERI_F3x_Pz_S_Py_C1000000003_b = I_ERI_G3xz_S_S_Py_C1000000003_b+ABZ*I_ERI_F3x_S_S_Py_C1000000003_b;
  Double I_ERI_F2xy_Pz_S_Py_C1000000003_b = I_ERI_G2xyz_S_S_Py_C1000000003_b+ABZ*I_ERI_F2xy_S_S_Py_C1000000003_b;
  Double I_ERI_F2xz_Pz_S_Py_C1000000003_b = I_ERI_G2x2z_S_S_Py_C1000000003_b+ABZ*I_ERI_F2xz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2y_Pz_S_Py_C1000000003_b = I_ERI_Gx2yz_S_S_Py_C1000000003_b+ABZ*I_ERI_Fx2y_S_S_Py_C1000000003_b;
  Double I_ERI_Fxyz_Pz_S_Py_C1000000003_b = I_ERI_Gxy2z_S_S_Py_C1000000003_b+ABZ*I_ERI_Fxyz_S_S_Py_C1000000003_b;
  Double I_ERI_Fx2z_Pz_S_Py_C1000000003_b = I_ERI_Gx3z_S_S_Py_C1000000003_b+ABZ*I_ERI_Fx2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3y_Pz_S_Py_C1000000003_b = I_ERI_G3yz_S_S_Py_C1000000003_b+ABZ*I_ERI_F3y_S_S_Py_C1000000003_b;
  Double I_ERI_F2yz_Pz_S_Py_C1000000003_b = I_ERI_G2y2z_S_S_Py_C1000000003_b+ABZ*I_ERI_F2yz_S_S_Py_C1000000003_b;
  Double I_ERI_Fy2z_Pz_S_Py_C1000000003_b = I_ERI_Gy3z_S_S_Py_C1000000003_b+ABZ*I_ERI_Fy2z_S_S_Py_C1000000003_b;
  Double I_ERI_F3z_Pz_S_Py_C1000000003_b = I_ERI_G4z_S_S_Py_C1000000003_b+ABZ*I_ERI_F3z_S_S_Py_C1000000003_b;
  Double I_ERI_F3x_Px_S_Pz_C1000000003_b = I_ERI_G4x_S_S_Pz_C1000000003_b+ABX*I_ERI_F3x_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xy_Px_S_Pz_C1000000003_b = I_ERI_G3xy_S_S_Pz_C1000000003_b+ABX*I_ERI_F2xy_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xz_Px_S_Pz_C1000000003_b = I_ERI_G3xz_S_S_Pz_C1000000003_b+ABX*I_ERI_F2xz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2y_Px_S_Pz_C1000000003_b = I_ERI_G2x2y_S_S_Pz_C1000000003_b+ABX*I_ERI_Fx2y_S_S_Pz_C1000000003_b;
  Double I_ERI_Fxyz_Px_S_Pz_C1000000003_b = I_ERI_G2xyz_S_S_Pz_C1000000003_b+ABX*I_ERI_Fxyz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2z_Px_S_Pz_C1000000003_b = I_ERI_G2x2z_S_S_Pz_C1000000003_b+ABX*I_ERI_Fx2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3y_Px_S_Pz_C1000000003_b = I_ERI_Gx3y_S_S_Pz_C1000000003_b+ABX*I_ERI_F3y_S_S_Pz_C1000000003_b;
  Double I_ERI_F2yz_Px_S_Pz_C1000000003_b = I_ERI_Gx2yz_S_S_Pz_C1000000003_b+ABX*I_ERI_F2yz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fy2z_Px_S_Pz_C1000000003_b = I_ERI_Gxy2z_S_S_Pz_C1000000003_b+ABX*I_ERI_Fy2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3z_Px_S_Pz_C1000000003_b = I_ERI_Gx3z_S_S_Pz_C1000000003_b+ABX*I_ERI_F3z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3x_Py_S_Pz_C1000000003_b = I_ERI_G3xy_S_S_Pz_C1000000003_b+ABY*I_ERI_F3x_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xy_Py_S_Pz_C1000000003_b = I_ERI_G2x2y_S_S_Pz_C1000000003_b+ABY*I_ERI_F2xy_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xz_Py_S_Pz_C1000000003_b = I_ERI_G2xyz_S_S_Pz_C1000000003_b+ABY*I_ERI_F2xz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2y_Py_S_Pz_C1000000003_b = I_ERI_Gx3y_S_S_Pz_C1000000003_b+ABY*I_ERI_Fx2y_S_S_Pz_C1000000003_b;
  Double I_ERI_Fxyz_Py_S_Pz_C1000000003_b = I_ERI_Gx2yz_S_S_Pz_C1000000003_b+ABY*I_ERI_Fxyz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2z_Py_S_Pz_C1000000003_b = I_ERI_Gxy2z_S_S_Pz_C1000000003_b+ABY*I_ERI_Fx2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3y_Py_S_Pz_C1000000003_b = I_ERI_G4y_S_S_Pz_C1000000003_b+ABY*I_ERI_F3y_S_S_Pz_C1000000003_b;
  Double I_ERI_F2yz_Py_S_Pz_C1000000003_b = I_ERI_G3yz_S_S_Pz_C1000000003_b+ABY*I_ERI_F2yz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fy2z_Py_S_Pz_C1000000003_b = I_ERI_G2y2z_S_S_Pz_C1000000003_b+ABY*I_ERI_Fy2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3z_Py_S_Pz_C1000000003_b = I_ERI_Gy3z_S_S_Pz_C1000000003_b+ABY*I_ERI_F3z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3x_Pz_S_Pz_C1000000003_b = I_ERI_G3xz_S_S_Pz_C1000000003_b+ABZ*I_ERI_F3x_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xy_Pz_S_Pz_C1000000003_b = I_ERI_G2xyz_S_S_Pz_C1000000003_b+ABZ*I_ERI_F2xy_S_S_Pz_C1000000003_b;
  Double I_ERI_F2xz_Pz_S_Pz_C1000000003_b = I_ERI_G2x2z_S_S_Pz_C1000000003_b+ABZ*I_ERI_F2xz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2y_Pz_S_Pz_C1000000003_b = I_ERI_Gx2yz_S_S_Pz_C1000000003_b+ABZ*I_ERI_Fx2y_S_S_Pz_C1000000003_b;
  Double I_ERI_Fxyz_Pz_S_Pz_C1000000003_b = I_ERI_Gxy2z_S_S_Pz_C1000000003_b+ABZ*I_ERI_Fxyz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fx2z_Pz_S_Pz_C1000000003_b = I_ERI_Gx3z_S_S_Pz_C1000000003_b+ABZ*I_ERI_Fx2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3y_Pz_S_Pz_C1000000003_b = I_ERI_G3yz_S_S_Pz_C1000000003_b+ABZ*I_ERI_F3y_S_S_Pz_C1000000003_b;
  Double I_ERI_F2yz_Pz_S_Pz_C1000000003_b = I_ERI_G2y2z_S_S_Pz_C1000000003_b+ABZ*I_ERI_F2yz_S_S_Pz_C1000000003_b;
  Double I_ERI_Fy2z_Pz_S_Pz_C1000000003_b = I_ERI_Gy3z_S_S_Pz_C1000000003_b+ABZ*I_ERI_Fy2z_S_S_Pz_C1000000003_b;
  Double I_ERI_F3z_Pz_S_Pz_C1000000003_b = I_ERI_G4z_S_S_Pz_C1000000003_b+ABZ*I_ERI_F3z_S_S_Pz_C1000000003_b;

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
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C3_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C3
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_G4x_S_S_S_C3_a-3*I_ERI_D2x_S_S_S_C3;
  abcd[1] = 2.0E0*I_ERI_G3xy_S_S_S_C3_a-2*I_ERI_Dxy_S_S_S_C3;
  abcd[2] = 2.0E0*I_ERI_G3xz_S_S_S_C3_a-2*I_ERI_Dxz_S_S_S_C3;
  abcd[3] = 2.0E0*I_ERI_G2x2y_S_S_S_C3_a-1*I_ERI_D2y_S_S_S_C3;
  abcd[4] = 2.0E0*I_ERI_G2xyz_S_S_S_C3_a-1*I_ERI_Dyz_S_S_S_C3;
  abcd[5] = 2.0E0*I_ERI_G2x2z_S_S_S_C3_a-1*I_ERI_D2z_S_S_S_C3;
  abcd[6] = 2.0E0*I_ERI_Gx3y_S_S_S_C3_a;
  abcd[7] = 2.0E0*I_ERI_Gx2yz_S_S_S_C3_a;
  abcd[8] = 2.0E0*I_ERI_Gxy2z_S_S_S_C3_a;
  abcd[9] = 2.0E0*I_ERI_Gx3z_S_S_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000003
   ************************************************************/
  abcd[10] = 2.0E0*I_ERI_G4x_S_Px_S_C1000003_a-3*I_ERI_D2x_S_Px_S_C1000003;
  abcd[11] = 2.0E0*I_ERI_G3xy_S_Px_S_C1000003_a-2*I_ERI_Dxy_S_Px_S_C1000003;
  abcd[12] = 2.0E0*I_ERI_G3xz_S_Px_S_C1000003_a-2*I_ERI_Dxz_S_Px_S_C1000003;
  abcd[13] = 2.0E0*I_ERI_G2x2y_S_Px_S_C1000003_a-1*I_ERI_D2y_S_Px_S_C1000003;
  abcd[14] = 2.0E0*I_ERI_G2xyz_S_Px_S_C1000003_a-1*I_ERI_Dyz_S_Px_S_C1000003;
  abcd[15] = 2.0E0*I_ERI_G2x2z_S_Px_S_C1000003_a-1*I_ERI_D2z_S_Px_S_C1000003;
  abcd[16] = 2.0E0*I_ERI_Gx3y_S_Px_S_C1000003_a;
  abcd[17] = 2.0E0*I_ERI_Gx2yz_S_Px_S_C1000003_a;
  abcd[18] = 2.0E0*I_ERI_Gxy2z_S_Px_S_C1000003_a;
  abcd[19] = 2.0E0*I_ERI_Gx3z_S_Px_S_C1000003_a;
  abcd[20] = 2.0E0*I_ERI_G4x_S_Py_S_C1000003_a-3*I_ERI_D2x_S_Py_S_C1000003;
  abcd[21] = 2.0E0*I_ERI_G3xy_S_Py_S_C1000003_a-2*I_ERI_Dxy_S_Py_S_C1000003;
  abcd[22] = 2.0E0*I_ERI_G3xz_S_Py_S_C1000003_a-2*I_ERI_Dxz_S_Py_S_C1000003;
  abcd[23] = 2.0E0*I_ERI_G2x2y_S_Py_S_C1000003_a-1*I_ERI_D2y_S_Py_S_C1000003;
  abcd[24] = 2.0E0*I_ERI_G2xyz_S_Py_S_C1000003_a-1*I_ERI_Dyz_S_Py_S_C1000003;
  abcd[25] = 2.0E0*I_ERI_G2x2z_S_Py_S_C1000003_a-1*I_ERI_D2z_S_Py_S_C1000003;
  abcd[26] = 2.0E0*I_ERI_Gx3y_S_Py_S_C1000003_a;
  abcd[27] = 2.0E0*I_ERI_Gx2yz_S_Py_S_C1000003_a;
  abcd[28] = 2.0E0*I_ERI_Gxy2z_S_Py_S_C1000003_a;
  abcd[29] = 2.0E0*I_ERI_Gx3z_S_Py_S_C1000003_a;
  abcd[30] = 2.0E0*I_ERI_G4x_S_Pz_S_C1000003_a-3*I_ERI_D2x_S_Pz_S_C1000003;
  abcd[31] = 2.0E0*I_ERI_G3xy_S_Pz_S_C1000003_a-2*I_ERI_Dxy_S_Pz_S_C1000003;
  abcd[32] = 2.0E0*I_ERI_G3xz_S_Pz_S_C1000003_a-2*I_ERI_Dxz_S_Pz_S_C1000003;
  abcd[33] = 2.0E0*I_ERI_G2x2y_S_Pz_S_C1000003_a-1*I_ERI_D2y_S_Pz_S_C1000003;
  abcd[34] = 2.0E0*I_ERI_G2xyz_S_Pz_S_C1000003_a-1*I_ERI_Dyz_S_Pz_S_C1000003;
  abcd[35] = 2.0E0*I_ERI_G2x2z_S_Pz_S_C1000003_a-1*I_ERI_D2z_S_Pz_S_C1000003;
  abcd[36] = 2.0E0*I_ERI_Gx3y_S_Pz_S_C1000003_a;
  abcd[37] = 2.0E0*I_ERI_Gx2yz_S_Pz_S_C1000003_a;
  abcd[38] = 2.0E0*I_ERI_Gxy2z_S_Pz_S_C1000003_a;
  abcd[39] = 2.0E0*I_ERI_Gx3z_S_Pz_S_C1000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_P_C1000000003_a
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000003
   ************************************************************/
  abcd[40] = 2.0E0*I_ERI_G4x_S_S_Px_C1000000003_a-3*I_ERI_D2x_S_S_Px_C1000000003;
  abcd[41] = 2.0E0*I_ERI_G3xy_S_S_Px_C1000000003_a-2*I_ERI_Dxy_S_S_Px_C1000000003;
  abcd[42] = 2.0E0*I_ERI_G3xz_S_S_Px_C1000000003_a-2*I_ERI_Dxz_S_S_Px_C1000000003;
  abcd[43] = 2.0E0*I_ERI_G2x2y_S_S_Px_C1000000003_a-1*I_ERI_D2y_S_S_Px_C1000000003;
  abcd[44] = 2.0E0*I_ERI_G2xyz_S_S_Px_C1000000003_a-1*I_ERI_Dyz_S_S_Px_C1000000003;
  abcd[45] = 2.0E0*I_ERI_G2x2z_S_S_Px_C1000000003_a-1*I_ERI_D2z_S_S_Px_C1000000003;
  abcd[46] = 2.0E0*I_ERI_Gx3y_S_S_Px_C1000000003_a;
  abcd[47] = 2.0E0*I_ERI_Gx2yz_S_S_Px_C1000000003_a;
  abcd[48] = 2.0E0*I_ERI_Gxy2z_S_S_Px_C1000000003_a;
  abcd[49] = 2.0E0*I_ERI_Gx3z_S_S_Px_C1000000003_a;
  abcd[80] = 2.0E0*I_ERI_G4x_S_S_Py_C1000000003_a-3*I_ERI_D2x_S_S_Py_C1000000003;
  abcd[81] = 2.0E0*I_ERI_G3xy_S_S_Py_C1000000003_a-2*I_ERI_Dxy_S_S_Py_C1000000003;
  abcd[82] = 2.0E0*I_ERI_G3xz_S_S_Py_C1000000003_a-2*I_ERI_Dxz_S_S_Py_C1000000003;
  abcd[83] = 2.0E0*I_ERI_G2x2y_S_S_Py_C1000000003_a-1*I_ERI_D2y_S_S_Py_C1000000003;
  abcd[84] = 2.0E0*I_ERI_G2xyz_S_S_Py_C1000000003_a-1*I_ERI_Dyz_S_S_Py_C1000000003;
  abcd[85] = 2.0E0*I_ERI_G2x2z_S_S_Py_C1000000003_a-1*I_ERI_D2z_S_S_Py_C1000000003;
  abcd[86] = 2.0E0*I_ERI_Gx3y_S_S_Py_C1000000003_a;
  abcd[87] = 2.0E0*I_ERI_Gx2yz_S_S_Py_C1000000003_a;
  abcd[88] = 2.0E0*I_ERI_Gxy2z_S_S_Py_C1000000003_a;
  abcd[89] = 2.0E0*I_ERI_Gx3z_S_S_Py_C1000000003_a;
  abcd[120] = 2.0E0*I_ERI_G4x_S_S_Pz_C1000000003_a-3*I_ERI_D2x_S_S_Pz_C1000000003;
  abcd[121] = 2.0E0*I_ERI_G3xy_S_S_Pz_C1000000003_a-2*I_ERI_Dxy_S_S_Pz_C1000000003;
  abcd[122] = 2.0E0*I_ERI_G3xz_S_S_Pz_C1000000003_a-2*I_ERI_Dxz_S_S_Pz_C1000000003;
  abcd[123] = 2.0E0*I_ERI_G2x2y_S_S_Pz_C1000000003_a-1*I_ERI_D2y_S_S_Pz_C1000000003;
  abcd[124] = 2.0E0*I_ERI_G2xyz_S_S_Pz_C1000000003_a-1*I_ERI_Dyz_S_S_Pz_C1000000003;
  abcd[125] = 2.0E0*I_ERI_G2x2z_S_S_Pz_C1000000003_a-1*I_ERI_D2z_S_S_Pz_C1000000003;
  abcd[126] = 2.0E0*I_ERI_Gx3y_S_S_Pz_C1000000003_a;
  abcd[127] = 2.0E0*I_ERI_Gx2yz_S_S_Pz_C1000000003_a;
  abcd[128] = 2.0E0*I_ERI_Gxy2z_S_S_Pz_C1000000003_a;
  abcd[129] = 2.0E0*I_ERI_Gx3z_S_S_Pz_C1000000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P_C1001000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000003
   ************************************************************/
  abcd[50] = 2.0E0*I_ERI_G4x_S_Px_Px_C1001000003_a-3*I_ERI_D2x_S_Px_Px_C1001000003;
  abcd[51] = 2.0E0*I_ERI_G3xy_S_Px_Px_C1001000003_a-2*I_ERI_Dxy_S_Px_Px_C1001000003;
  abcd[52] = 2.0E0*I_ERI_G3xz_S_Px_Px_C1001000003_a-2*I_ERI_Dxz_S_Px_Px_C1001000003;
  abcd[53] = 2.0E0*I_ERI_G2x2y_S_Px_Px_C1001000003_a-1*I_ERI_D2y_S_Px_Px_C1001000003;
  abcd[54] = 2.0E0*I_ERI_G2xyz_S_Px_Px_C1001000003_a-1*I_ERI_Dyz_S_Px_Px_C1001000003;
  abcd[55] = 2.0E0*I_ERI_G2x2z_S_Px_Px_C1001000003_a-1*I_ERI_D2z_S_Px_Px_C1001000003;
  abcd[56] = 2.0E0*I_ERI_Gx3y_S_Px_Px_C1001000003_a;
  abcd[57] = 2.0E0*I_ERI_Gx2yz_S_Px_Px_C1001000003_a;
  abcd[58] = 2.0E0*I_ERI_Gxy2z_S_Px_Px_C1001000003_a;
  abcd[59] = 2.0E0*I_ERI_Gx3z_S_Px_Px_C1001000003_a;
  abcd[60] = 2.0E0*I_ERI_G4x_S_Py_Px_C1001000003_a-3*I_ERI_D2x_S_Py_Px_C1001000003;
  abcd[61] = 2.0E0*I_ERI_G3xy_S_Py_Px_C1001000003_a-2*I_ERI_Dxy_S_Py_Px_C1001000003;
  abcd[62] = 2.0E0*I_ERI_G3xz_S_Py_Px_C1001000003_a-2*I_ERI_Dxz_S_Py_Px_C1001000003;
  abcd[63] = 2.0E0*I_ERI_G2x2y_S_Py_Px_C1001000003_a-1*I_ERI_D2y_S_Py_Px_C1001000003;
  abcd[64] = 2.0E0*I_ERI_G2xyz_S_Py_Px_C1001000003_a-1*I_ERI_Dyz_S_Py_Px_C1001000003;
  abcd[65] = 2.0E0*I_ERI_G2x2z_S_Py_Px_C1001000003_a-1*I_ERI_D2z_S_Py_Px_C1001000003;
  abcd[66] = 2.0E0*I_ERI_Gx3y_S_Py_Px_C1001000003_a;
  abcd[67] = 2.0E0*I_ERI_Gx2yz_S_Py_Px_C1001000003_a;
  abcd[68] = 2.0E0*I_ERI_Gxy2z_S_Py_Px_C1001000003_a;
  abcd[69] = 2.0E0*I_ERI_Gx3z_S_Py_Px_C1001000003_a;
  abcd[70] = 2.0E0*I_ERI_G4x_S_Pz_Px_C1001000003_a-3*I_ERI_D2x_S_Pz_Px_C1001000003;
  abcd[71] = 2.0E0*I_ERI_G3xy_S_Pz_Px_C1001000003_a-2*I_ERI_Dxy_S_Pz_Px_C1001000003;
  abcd[72] = 2.0E0*I_ERI_G3xz_S_Pz_Px_C1001000003_a-2*I_ERI_Dxz_S_Pz_Px_C1001000003;
  abcd[73] = 2.0E0*I_ERI_G2x2y_S_Pz_Px_C1001000003_a-1*I_ERI_D2y_S_Pz_Px_C1001000003;
  abcd[74] = 2.0E0*I_ERI_G2xyz_S_Pz_Px_C1001000003_a-1*I_ERI_Dyz_S_Pz_Px_C1001000003;
  abcd[75] = 2.0E0*I_ERI_G2x2z_S_Pz_Px_C1001000003_a-1*I_ERI_D2z_S_Pz_Px_C1001000003;
  abcd[76] = 2.0E0*I_ERI_Gx3y_S_Pz_Px_C1001000003_a;
  abcd[77] = 2.0E0*I_ERI_Gx2yz_S_Pz_Px_C1001000003_a;
  abcd[78] = 2.0E0*I_ERI_Gxy2z_S_Pz_Px_C1001000003_a;
  abcd[79] = 2.0E0*I_ERI_Gx3z_S_Pz_Px_C1001000003_a;
  abcd[90] = 2.0E0*I_ERI_G4x_S_Px_Py_C1001000003_a-3*I_ERI_D2x_S_Px_Py_C1001000003;
  abcd[91] = 2.0E0*I_ERI_G3xy_S_Px_Py_C1001000003_a-2*I_ERI_Dxy_S_Px_Py_C1001000003;
  abcd[92] = 2.0E0*I_ERI_G3xz_S_Px_Py_C1001000003_a-2*I_ERI_Dxz_S_Px_Py_C1001000003;
  abcd[93] = 2.0E0*I_ERI_G2x2y_S_Px_Py_C1001000003_a-1*I_ERI_D2y_S_Px_Py_C1001000003;
  abcd[94] = 2.0E0*I_ERI_G2xyz_S_Px_Py_C1001000003_a-1*I_ERI_Dyz_S_Px_Py_C1001000003;
  abcd[95] = 2.0E0*I_ERI_G2x2z_S_Px_Py_C1001000003_a-1*I_ERI_D2z_S_Px_Py_C1001000003;
  abcd[96] = 2.0E0*I_ERI_Gx3y_S_Px_Py_C1001000003_a;
  abcd[97] = 2.0E0*I_ERI_Gx2yz_S_Px_Py_C1001000003_a;
  abcd[98] = 2.0E0*I_ERI_Gxy2z_S_Px_Py_C1001000003_a;
  abcd[99] = 2.0E0*I_ERI_Gx3z_S_Px_Py_C1001000003_a;
  abcd[100] = 2.0E0*I_ERI_G4x_S_Py_Py_C1001000003_a-3*I_ERI_D2x_S_Py_Py_C1001000003;
  abcd[101] = 2.0E0*I_ERI_G3xy_S_Py_Py_C1001000003_a-2*I_ERI_Dxy_S_Py_Py_C1001000003;
  abcd[102] = 2.0E0*I_ERI_G3xz_S_Py_Py_C1001000003_a-2*I_ERI_Dxz_S_Py_Py_C1001000003;
  abcd[103] = 2.0E0*I_ERI_G2x2y_S_Py_Py_C1001000003_a-1*I_ERI_D2y_S_Py_Py_C1001000003;
  abcd[104] = 2.0E0*I_ERI_G2xyz_S_Py_Py_C1001000003_a-1*I_ERI_Dyz_S_Py_Py_C1001000003;
  abcd[105] = 2.0E0*I_ERI_G2x2z_S_Py_Py_C1001000003_a-1*I_ERI_D2z_S_Py_Py_C1001000003;
  abcd[106] = 2.0E0*I_ERI_Gx3y_S_Py_Py_C1001000003_a;
  abcd[107] = 2.0E0*I_ERI_Gx2yz_S_Py_Py_C1001000003_a;
  abcd[108] = 2.0E0*I_ERI_Gxy2z_S_Py_Py_C1001000003_a;
  abcd[109] = 2.0E0*I_ERI_Gx3z_S_Py_Py_C1001000003_a;
  abcd[110] = 2.0E0*I_ERI_G4x_S_Pz_Py_C1001000003_a-3*I_ERI_D2x_S_Pz_Py_C1001000003;
  abcd[111] = 2.0E0*I_ERI_G3xy_S_Pz_Py_C1001000003_a-2*I_ERI_Dxy_S_Pz_Py_C1001000003;
  abcd[112] = 2.0E0*I_ERI_G3xz_S_Pz_Py_C1001000003_a-2*I_ERI_Dxz_S_Pz_Py_C1001000003;
  abcd[113] = 2.0E0*I_ERI_G2x2y_S_Pz_Py_C1001000003_a-1*I_ERI_D2y_S_Pz_Py_C1001000003;
  abcd[114] = 2.0E0*I_ERI_G2xyz_S_Pz_Py_C1001000003_a-1*I_ERI_Dyz_S_Pz_Py_C1001000003;
  abcd[115] = 2.0E0*I_ERI_G2x2z_S_Pz_Py_C1001000003_a-1*I_ERI_D2z_S_Pz_Py_C1001000003;
  abcd[116] = 2.0E0*I_ERI_Gx3y_S_Pz_Py_C1001000003_a;
  abcd[117] = 2.0E0*I_ERI_Gx2yz_S_Pz_Py_C1001000003_a;
  abcd[118] = 2.0E0*I_ERI_Gxy2z_S_Pz_Py_C1001000003_a;
  abcd[119] = 2.0E0*I_ERI_Gx3z_S_Pz_Py_C1001000003_a;
  abcd[130] = 2.0E0*I_ERI_G4x_S_Px_Pz_C1001000003_a-3*I_ERI_D2x_S_Px_Pz_C1001000003;
  abcd[131] = 2.0E0*I_ERI_G3xy_S_Px_Pz_C1001000003_a-2*I_ERI_Dxy_S_Px_Pz_C1001000003;
  abcd[132] = 2.0E0*I_ERI_G3xz_S_Px_Pz_C1001000003_a-2*I_ERI_Dxz_S_Px_Pz_C1001000003;
  abcd[133] = 2.0E0*I_ERI_G2x2y_S_Px_Pz_C1001000003_a-1*I_ERI_D2y_S_Px_Pz_C1001000003;
  abcd[134] = 2.0E0*I_ERI_G2xyz_S_Px_Pz_C1001000003_a-1*I_ERI_Dyz_S_Px_Pz_C1001000003;
  abcd[135] = 2.0E0*I_ERI_G2x2z_S_Px_Pz_C1001000003_a-1*I_ERI_D2z_S_Px_Pz_C1001000003;
  abcd[136] = 2.0E0*I_ERI_Gx3y_S_Px_Pz_C1001000003_a;
  abcd[137] = 2.0E0*I_ERI_Gx2yz_S_Px_Pz_C1001000003_a;
  abcd[138] = 2.0E0*I_ERI_Gxy2z_S_Px_Pz_C1001000003_a;
  abcd[139] = 2.0E0*I_ERI_Gx3z_S_Px_Pz_C1001000003_a;
  abcd[140] = 2.0E0*I_ERI_G4x_S_Py_Pz_C1001000003_a-3*I_ERI_D2x_S_Py_Pz_C1001000003;
  abcd[141] = 2.0E0*I_ERI_G3xy_S_Py_Pz_C1001000003_a-2*I_ERI_Dxy_S_Py_Pz_C1001000003;
  abcd[142] = 2.0E0*I_ERI_G3xz_S_Py_Pz_C1001000003_a-2*I_ERI_Dxz_S_Py_Pz_C1001000003;
  abcd[143] = 2.0E0*I_ERI_G2x2y_S_Py_Pz_C1001000003_a-1*I_ERI_D2y_S_Py_Pz_C1001000003;
  abcd[144] = 2.0E0*I_ERI_G2xyz_S_Py_Pz_C1001000003_a-1*I_ERI_Dyz_S_Py_Pz_C1001000003;
  abcd[145] = 2.0E0*I_ERI_G2x2z_S_Py_Pz_C1001000003_a-1*I_ERI_D2z_S_Py_Pz_C1001000003;
  abcd[146] = 2.0E0*I_ERI_Gx3y_S_Py_Pz_C1001000003_a;
  abcd[147] = 2.0E0*I_ERI_Gx2yz_S_Py_Pz_C1001000003_a;
  abcd[148] = 2.0E0*I_ERI_Gxy2z_S_Py_Pz_C1001000003_a;
  abcd[149] = 2.0E0*I_ERI_Gx3z_S_Py_Pz_C1001000003_a;
  abcd[150] = 2.0E0*I_ERI_G4x_S_Pz_Pz_C1001000003_a-3*I_ERI_D2x_S_Pz_Pz_C1001000003;
  abcd[151] = 2.0E0*I_ERI_G3xy_S_Pz_Pz_C1001000003_a-2*I_ERI_Dxy_S_Pz_Pz_C1001000003;
  abcd[152] = 2.0E0*I_ERI_G3xz_S_Pz_Pz_C1001000003_a-2*I_ERI_Dxz_S_Pz_Pz_C1001000003;
  abcd[153] = 2.0E0*I_ERI_G2x2y_S_Pz_Pz_C1001000003_a-1*I_ERI_D2y_S_Pz_Pz_C1001000003;
  abcd[154] = 2.0E0*I_ERI_G2xyz_S_Pz_Pz_C1001000003_a-1*I_ERI_Dyz_S_Pz_Pz_C1001000003;
  abcd[155] = 2.0E0*I_ERI_G2x2z_S_Pz_Pz_C1001000003_a-1*I_ERI_D2z_S_Pz_Pz_C1001000003;
  abcd[156] = 2.0E0*I_ERI_Gx3y_S_Pz_Pz_C1001000003_a;
  abcd[157] = 2.0E0*I_ERI_Gx2yz_S_Pz_Pz_C1001000003_a;
  abcd[158] = 2.0E0*I_ERI_Gxy2z_S_Pz_Pz_C1001000003_a;
  abcd[159] = 2.0E0*I_ERI_Gx3z_S_Pz_Pz_C1001000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C3_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C3
   ************************************************************/
  abcd[160] = 2.0E0*I_ERI_G3xy_S_S_S_C3_a;
  abcd[161] = 2.0E0*I_ERI_G2x2y_S_S_S_C3_a-1*I_ERI_D2x_S_S_S_C3;
  abcd[162] = 2.0E0*I_ERI_G2xyz_S_S_S_C3_a;
  abcd[163] = 2.0E0*I_ERI_Gx3y_S_S_S_C3_a-2*I_ERI_Dxy_S_S_S_C3;
  abcd[164] = 2.0E0*I_ERI_Gx2yz_S_S_S_C3_a-1*I_ERI_Dxz_S_S_S_C3;
  abcd[165] = 2.0E0*I_ERI_Gxy2z_S_S_S_C3_a;
  abcd[166] = 2.0E0*I_ERI_G4y_S_S_S_C3_a-3*I_ERI_D2y_S_S_S_C3;
  abcd[167] = 2.0E0*I_ERI_G3yz_S_S_S_C3_a-2*I_ERI_Dyz_S_S_S_C3;
  abcd[168] = 2.0E0*I_ERI_G2y2z_S_S_S_C3_a-1*I_ERI_D2z_S_S_S_C3;
  abcd[169] = 2.0E0*I_ERI_Gy3z_S_S_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000003
   ************************************************************/
  abcd[170] = 2.0E0*I_ERI_G3xy_S_Px_S_C1000003_a;
  abcd[171] = 2.0E0*I_ERI_G2x2y_S_Px_S_C1000003_a-1*I_ERI_D2x_S_Px_S_C1000003;
  abcd[172] = 2.0E0*I_ERI_G2xyz_S_Px_S_C1000003_a;
  abcd[173] = 2.0E0*I_ERI_Gx3y_S_Px_S_C1000003_a-2*I_ERI_Dxy_S_Px_S_C1000003;
  abcd[174] = 2.0E0*I_ERI_Gx2yz_S_Px_S_C1000003_a-1*I_ERI_Dxz_S_Px_S_C1000003;
  abcd[175] = 2.0E0*I_ERI_Gxy2z_S_Px_S_C1000003_a;
  abcd[176] = 2.0E0*I_ERI_G4y_S_Px_S_C1000003_a-3*I_ERI_D2y_S_Px_S_C1000003;
  abcd[177] = 2.0E0*I_ERI_G3yz_S_Px_S_C1000003_a-2*I_ERI_Dyz_S_Px_S_C1000003;
  abcd[178] = 2.0E0*I_ERI_G2y2z_S_Px_S_C1000003_a-1*I_ERI_D2z_S_Px_S_C1000003;
  abcd[179] = 2.0E0*I_ERI_Gy3z_S_Px_S_C1000003_a;
  abcd[180] = 2.0E0*I_ERI_G3xy_S_Py_S_C1000003_a;
  abcd[181] = 2.0E0*I_ERI_G2x2y_S_Py_S_C1000003_a-1*I_ERI_D2x_S_Py_S_C1000003;
  abcd[182] = 2.0E0*I_ERI_G2xyz_S_Py_S_C1000003_a;
  abcd[183] = 2.0E0*I_ERI_Gx3y_S_Py_S_C1000003_a-2*I_ERI_Dxy_S_Py_S_C1000003;
  abcd[184] = 2.0E0*I_ERI_Gx2yz_S_Py_S_C1000003_a-1*I_ERI_Dxz_S_Py_S_C1000003;
  abcd[185] = 2.0E0*I_ERI_Gxy2z_S_Py_S_C1000003_a;
  abcd[186] = 2.0E0*I_ERI_G4y_S_Py_S_C1000003_a-3*I_ERI_D2y_S_Py_S_C1000003;
  abcd[187] = 2.0E0*I_ERI_G3yz_S_Py_S_C1000003_a-2*I_ERI_Dyz_S_Py_S_C1000003;
  abcd[188] = 2.0E0*I_ERI_G2y2z_S_Py_S_C1000003_a-1*I_ERI_D2z_S_Py_S_C1000003;
  abcd[189] = 2.0E0*I_ERI_Gy3z_S_Py_S_C1000003_a;
  abcd[190] = 2.0E0*I_ERI_G3xy_S_Pz_S_C1000003_a;
  abcd[191] = 2.0E0*I_ERI_G2x2y_S_Pz_S_C1000003_a-1*I_ERI_D2x_S_Pz_S_C1000003;
  abcd[192] = 2.0E0*I_ERI_G2xyz_S_Pz_S_C1000003_a;
  abcd[193] = 2.0E0*I_ERI_Gx3y_S_Pz_S_C1000003_a-2*I_ERI_Dxy_S_Pz_S_C1000003;
  abcd[194] = 2.0E0*I_ERI_Gx2yz_S_Pz_S_C1000003_a-1*I_ERI_Dxz_S_Pz_S_C1000003;
  abcd[195] = 2.0E0*I_ERI_Gxy2z_S_Pz_S_C1000003_a;
  abcd[196] = 2.0E0*I_ERI_G4y_S_Pz_S_C1000003_a-3*I_ERI_D2y_S_Pz_S_C1000003;
  abcd[197] = 2.0E0*I_ERI_G3yz_S_Pz_S_C1000003_a-2*I_ERI_Dyz_S_Pz_S_C1000003;
  abcd[198] = 2.0E0*I_ERI_G2y2z_S_Pz_S_C1000003_a-1*I_ERI_D2z_S_Pz_S_C1000003;
  abcd[199] = 2.0E0*I_ERI_Gy3z_S_Pz_S_C1000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_P_C1000000003_a
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000003
   ************************************************************/
  abcd[200] = 2.0E0*I_ERI_G3xy_S_S_Px_C1000000003_a;
  abcd[201] = 2.0E0*I_ERI_G2x2y_S_S_Px_C1000000003_a-1*I_ERI_D2x_S_S_Px_C1000000003;
  abcd[202] = 2.0E0*I_ERI_G2xyz_S_S_Px_C1000000003_a;
  abcd[203] = 2.0E0*I_ERI_Gx3y_S_S_Px_C1000000003_a-2*I_ERI_Dxy_S_S_Px_C1000000003;
  abcd[204] = 2.0E0*I_ERI_Gx2yz_S_S_Px_C1000000003_a-1*I_ERI_Dxz_S_S_Px_C1000000003;
  abcd[205] = 2.0E0*I_ERI_Gxy2z_S_S_Px_C1000000003_a;
  abcd[206] = 2.0E0*I_ERI_G4y_S_S_Px_C1000000003_a-3*I_ERI_D2y_S_S_Px_C1000000003;
  abcd[207] = 2.0E0*I_ERI_G3yz_S_S_Px_C1000000003_a-2*I_ERI_Dyz_S_S_Px_C1000000003;
  abcd[208] = 2.0E0*I_ERI_G2y2z_S_S_Px_C1000000003_a-1*I_ERI_D2z_S_S_Px_C1000000003;
  abcd[209] = 2.0E0*I_ERI_Gy3z_S_S_Px_C1000000003_a;
  abcd[240] = 2.0E0*I_ERI_G3xy_S_S_Py_C1000000003_a;
  abcd[241] = 2.0E0*I_ERI_G2x2y_S_S_Py_C1000000003_a-1*I_ERI_D2x_S_S_Py_C1000000003;
  abcd[242] = 2.0E0*I_ERI_G2xyz_S_S_Py_C1000000003_a;
  abcd[243] = 2.0E0*I_ERI_Gx3y_S_S_Py_C1000000003_a-2*I_ERI_Dxy_S_S_Py_C1000000003;
  abcd[244] = 2.0E0*I_ERI_Gx2yz_S_S_Py_C1000000003_a-1*I_ERI_Dxz_S_S_Py_C1000000003;
  abcd[245] = 2.0E0*I_ERI_Gxy2z_S_S_Py_C1000000003_a;
  abcd[246] = 2.0E0*I_ERI_G4y_S_S_Py_C1000000003_a-3*I_ERI_D2y_S_S_Py_C1000000003;
  abcd[247] = 2.0E0*I_ERI_G3yz_S_S_Py_C1000000003_a-2*I_ERI_Dyz_S_S_Py_C1000000003;
  abcd[248] = 2.0E0*I_ERI_G2y2z_S_S_Py_C1000000003_a-1*I_ERI_D2z_S_S_Py_C1000000003;
  abcd[249] = 2.0E0*I_ERI_Gy3z_S_S_Py_C1000000003_a;
  abcd[280] = 2.0E0*I_ERI_G3xy_S_S_Pz_C1000000003_a;
  abcd[281] = 2.0E0*I_ERI_G2x2y_S_S_Pz_C1000000003_a-1*I_ERI_D2x_S_S_Pz_C1000000003;
  abcd[282] = 2.0E0*I_ERI_G2xyz_S_S_Pz_C1000000003_a;
  abcd[283] = 2.0E0*I_ERI_Gx3y_S_S_Pz_C1000000003_a-2*I_ERI_Dxy_S_S_Pz_C1000000003;
  abcd[284] = 2.0E0*I_ERI_Gx2yz_S_S_Pz_C1000000003_a-1*I_ERI_Dxz_S_S_Pz_C1000000003;
  abcd[285] = 2.0E0*I_ERI_Gxy2z_S_S_Pz_C1000000003_a;
  abcd[286] = 2.0E0*I_ERI_G4y_S_S_Pz_C1000000003_a-3*I_ERI_D2y_S_S_Pz_C1000000003;
  abcd[287] = 2.0E0*I_ERI_G3yz_S_S_Pz_C1000000003_a-2*I_ERI_Dyz_S_S_Pz_C1000000003;
  abcd[288] = 2.0E0*I_ERI_G2y2z_S_S_Pz_C1000000003_a-1*I_ERI_D2z_S_S_Pz_C1000000003;
  abcd[289] = 2.0E0*I_ERI_Gy3z_S_S_Pz_C1000000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P_C1001000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000003
   ************************************************************/
  abcd[210] = 2.0E0*I_ERI_G3xy_S_Px_Px_C1001000003_a;
  abcd[211] = 2.0E0*I_ERI_G2x2y_S_Px_Px_C1001000003_a-1*I_ERI_D2x_S_Px_Px_C1001000003;
  abcd[212] = 2.0E0*I_ERI_G2xyz_S_Px_Px_C1001000003_a;
  abcd[213] = 2.0E0*I_ERI_Gx3y_S_Px_Px_C1001000003_a-2*I_ERI_Dxy_S_Px_Px_C1001000003;
  abcd[214] = 2.0E0*I_ERI_Gx2yz_S_Px_Px_C1001000003_a-1*I_ERI_Dxz_S_Px_Px_C1001000003;
  abcd[215] = 2.0E0*I_ERI_Gxy2z_S_Px_Px_C1001000003_a;
  abcd[216] = 2.0E0*I_ERI_G4y_S_Px_Px_C1001000003_a-3*I_ERI_D2y_S_Px_Px_C1001000003;
  abcd[217] = 2.0E0*I_ERI_G3yz_S_Px_Px_C1001000003_a-2*I_ERI_Dyz_S_Px_Px_C1001000003;
  abcd[218] = 2.0E0*I_ERI_G2y2z_S_Px_Px_C1001000003_a-1*I_ERI_D2z_S_Px_Px_C1001000003;
  abcd[219] = 2.0E0*I_ERI_Gy3z_S_Px_Px_C1001000003_a;
  abcd[220] = 2.0E0*I_ERI_G3xy_S_Py_Px_C1001000003_a;
  abcd[221] = 2.0E0*I_ERI_G2x2y_S_Py_Px_C1001000003_a-1*I_ERI_D2x_S_Py_Px_C1001000003;
  abcd[222] = 2.0E0*I_ERI_G2xyz_S_Py_Px_C1001000003_a;
  abcd[223] = 2.0E0*I_ERI_Gx3y_S_Py_Px_C1001000003_a-2*I_ERI_Dxy_S_Py_Px_C1001000003;
  abcd[224] = 2.0E0*I_ERI_Gx2yz_S_Py_Px_C1001000003_a-1*I_ERI_Dxz_S_Py_Px_C1001000003;
  abcd[225] = 2.0E0*I_ERI_Gxy2z_S_Py_Px_C1001000003_a;
  abcd[226] = 2.0E0*I_ERI_G4y_S_Py_Px_C1001000003_a-3*I_ERI_D2y_S_Py_Px_C1001000003;
  abcd[227] = 2.0E0*I_ERI_G3yz_S_Py_Px_C1001000003_a-2*I_ERI_Dyz_S_Py_Px_C1001000003;
  abcd[228] = 2.0E0*I_ERI_G2y2z_S_Py_Px_C1001000003_a-1*I_ERI_D2z_S_Py_Px_C1001000003;
  abcd[229] = 2.0E0*I_ERI_Gy3z_S_Py_Px_C1001000003_a;
  abcd[230] = 2.0E0*I_ERI_G3xy_S_Pz_Px_C1001000003_a;
  abcd[231] = 2.0E0*I_ERI_G2x2y_S_Pz_Px_C1001000003_a-1*I_ERI_D2x_S_Pz_Px_C1001000003;
  abcd[232] = 2.0E0*I_ERI_G2xyz_S_Pz_Px_C1001000003_a;
  abcd[233] = 2.0E0*I_ERI_Gx3y_S_Pz_Px_C1001000003_a-2*I_ERI_Dxy_S_Pz_Px_C1001000003;
  abcd[234] = 2.0E0*I_ERI_Gx2yz_S_Pz_Px_C1001000003_a-1*I_ERI_Dxz_S_Pz_Px_C1001000003;
  abcd[235] = 2.0E0*I_ERI_Gxy2z_S_Pz_Px_C1001000003_a;
  abcd[236] = 2.0E0*I_ERI_G4y_S_Pz_Px_C1001000003_a-3*I_ERI_D2y_S_Pz_Px_C1001000003;
  abcd[237] = 2.0E0*I_ERI_G3yz_S_Pz_Px_C1001000003_a-2*I_ERI_Dyz_S_Pz_Px_C1001000003;
  abcd[238] = 2.0E0*I_ERI_G2y2z_S_Pz_Px_C1001000003_a-1*I_ERI_D2z_S_Pz_Px_C1001000003;
  abcd[239] = 2.0E0*I_ERI_Gy3z_S_Pz_Px_C1001000003_a;
  abcd[250] = 2.0E0*I_ERI_G3xy_S_Px_Py_C1001000003_a;
  abcd[251] = 2.0E0*I_ERI_G2x2y_S_Px_Py_C1001000003_a-1*I_ERI_D2x_S_Px_Py_C1001000003;
  abcd[252] = 2.0E0*I_ERI_G2xyz_S_Px_Py_C1001000003_a;
  abcd[253] = 2.0E0*I_ERI_Gx3y_S_Px_Py_C1001000003_a-2*I_ERI_Dxy_S_Px_Py_C1001000003;
  abcd[254] = 2.0E0*I_ERI_Gx2yz_S_Px_Py_C1001000003_a-1*I_ERI_Dxz_S_Px_Py_C1001000003;
  abcd[255] = 2.0E0*I_ERI_Gxy2z_S_Px_Py_C1001000003_a;
  abcd[256] = 2.0E0*I_ERI_G4y_S_Px_Py_C1001000003_a-3*I_ERI_D2y_S_Px_Py_C1001000003;
  abcd[257] = 2.0E0*I_ERI_G3yz_S_Px_Py_C1001000003_a-2*I_ERI_Dyz_S_Px_Py_C1001000003;
  abcd[258] = 2.0E0*I_ERI_G2y2z_S_Px_Py_C1001000003_a-1*I_ERI_D2z_S_Px_Py_C1001000003;
  abcd[259] = 2.0E0*I_ERI_Gy3z_S_Px_Py_C1001000003_a;
  abcd[260] = 2.0E0*I_ERI_G3xy_S_Py_Py_C1001000003_a;
  abcd[261] = 2.0E0*I_ERI_G2x2y_S_Py_Py_C1001000003_a-1*I_ERI_D2x_S_Py_Py_C1001000003;
  abcd[262] = 2.0E0*I_ERI_G2xyz_S_Py_Py_C1001000003_a;
  abcd[263] = 2.0E0*I_ERI_Gx3y_S_Py_Py_C1001000003_a-2*I_ERI_Dxy_S_Py_Py_C1001000003;
  abcd[264] = 2.0E0*I_ERI_Gx2yz_S_Py_Py_C1001000003_a-1*I_ERI_Dxz_S_Py_Py_C1001000003;
  abcd[265] = 2.0E0*I_ERI_Gxy2z_S_Py_Py_C1001000003_a;
  abcd[266] = 2.0E0*I_ERI_G4y_S_Py_Py_C1001000003_a-3*I_ERI_D2y_S_Py_Py_C1001000003;
  abcd[267] = 2.0E0*I_ERI_G3yz_S_Py_Py_C1001000003_a-2*I_ERI_Dyz_S_Py_Py_C1001000003;
  abcd[268] = 2.0E0*I_ERI_G2y2z_S_Py_Py_C1001000003_a-1*I_ERI_D2z_S_Py_Py_C1001000003;
  abcd[269] = 2.0E0*I_ERI_Gy3z_S_Py_Py_C1001000003_a;
  abcd[270] = 2.0E0*I_ERI_G3xy_S_Pz_Py_C1001000003_a;
  abcd[271] = 2.0E0*I_ERI_G2x2y_S_Pz_Py_C1001000003_a-1*I_ERI_D2x_S_Pz_Py_C1001000003;
  abcd[272] = 2.0E0*I_ERI_G2xyz_S_Pz_Py_C1001000003_a;
  abcd[273] = 2.0E0*I_ERI_Gx3y_S_Pz_Py_C1001000003_a-2*I_ERI_Dxy_S_Pz_Py_C1001000003;
  abcd[274] = 2.0E0*I_ERI_Gx2yz_S_Pz_Py_C1001000003_a-1*I_ERI_Dxz_S_Pz_Py_C1001000003;
  abcd[275] = 2.0E0*I_ERI_Gxy2z_S_Pz_Py_C1001000003_a;
  abcd[276] = 2.0E0*I_ERI_G4y_S_Pz_Py_C1001000003_a-3*I_ERI_D2y_S_Pz_Py_C1001000003;
  abcd[277] = 2.0E0*I_ERI_G3yz_S_Pz_Py_C1001000003_a-2*I_ERI_Dyz_S_Pz_Py_C1001000003;
  abcd[278] = 2.0E0*I_ERI_G2y2z_S_Pz_Py_C1001000003_a-1*I_ERI_D2z_S_Pz_Py_C1001000003;
  abcd[279] = 2.0E0*I_ERI_Gy3z_S_Pz_Py_C1001000003_a;
  abcd[290] = 2.0E0*I_ERI_G3xy_S_Px_Pz_C1001000003_a;
  abcd[291] = 2.0E0*I_ERI_G2x2y_S_Px_Pz_C1001000003_a-1*I_ERI_D2x_S_Px_Pz_C1001000003;
  abcd[292] = 2.0E0*I_ERI_G2xyz_S_Px_Pz_C1001000003_a;
  abcd[293] = 2.0E0*I_ERI_Gx3y_S_Px_Pz_C1001000003_a-2*I_ERI_Dxy_S_Px_Pz_C1001000003;
  abcd[294] = 2.0E0*I_ERI_Gx2yz_S_Px_Pz_C1001000003_a-1*I_ERI_Dxz_S_Px_Pz_C1001000003;
  abcd[295] = 2.0E0*I_ERI_Gxy2z_S_Px_Pz_C1001000003_a;
  abcd[296] = 2.0E0*I_ERI_G4y_S_Px_Pz_C1001000003_a-3*I_ERI_D2y_S_Px_Pz_C1001000003;
  abcd[297] = 2.0E0*I_ERI_G3yz_S_Px_Pz_C1001000003_a-2*I_ERI_Dyz_S_Px_Pz_C1001000003;
  abcd[298] = 2.0E0*I_ERI_G2y2z_S_Px_Pz_C1001000003_a-1*I_ERI_D2z_S_Px_Pz_C1001000003;
  abcd[299] = 2.0E0*I_ERI_Gy3z_S_Px_Pz_C1001000003_a;
  abcd[300] = 2.0E0*I_ERI_G3xy_S_Py_Pz_C1001000003_a;
  abcd[301] = 2.0E0*I_ERI_G2x2y_S_Py_Pz_C1001000003_a-1*I_ERI_D2x_S_Py_Pz_C1001000003;
  abcd[302] = 2.0E0*I_ERI_G2xyz_S_Py_Pz_C1001000003_a;
  abcd[303] = 2.0E0*I_ERI_Gx3y_S_Py_Pz_C1001000003_a-2*I_ERI_Dxy_S_Py_Pz_C1001000003;
  abcd[304] = 2.0E0*I_ERI_Gx2yz_S_Py_Pz_C1001000003_a-1*I_ERI_Dxz_S_Py_Pz_C1001000003;
  abcd[305] = 2.0E0*I_ERI_Gxy2z_S_Py_Pz_C1001000003_a;
  abcd[306] = 2.0E0*I_ERI_G4y_S_Py_Pz_C1001000003_a-3*I_ERI_D2y_S_Py_Pz_C1001000003;
  abcd[307] = 2.0E0*I_ERI_G3yz_S_Py_Pz_C1001000003_a-2*I_ERI_Dyz_S_Py_Pz_C1001000003;
  abcd[308] = 2.0E0*I_ERI_G2y2z_S_Py_Pz_C1001000003_a-1*I_ERI_D2z_S_Py_Pz_C1001000003;
  abcd[309] = 2.0E0*I_ERI_Gy3z_S_Py_Pz_C1001000003_a;
  abcd[310] = 2.0E0*I_ERI_G3xy_S_Pz_Pz_C1001000003_a;
  abcd[311] = 2.0E0*I_ERI_G2x2y_S_Pz_Pz_C1001000003_a-1*I_ERI_D2x_S_Pz_Pz_C1001000003;
  abcd[312] = 2.0E0*I_ERI_G2xyz_S_Pz_Pz_C1001000003_a;
  abcd[313] = 2.0E0*I_ERI_Gx3y_S_Pz_Pz_C1001000003_a-2*I_ERI_Dxy_S_Pz_Pz_C1001000003;
  abcd[314] = 2.0E0*I_ERI_Gx2yz_S_Pz_Pz_C1001000003_a-1*I_ERI_Dxz_S_Pz_Pz_C1001000003;
  abcd[315] = 2.0E0*I_ERI_Gxy2z_S_Pz_Pz_C1001000003_a;
  abcd[316] = 2.0E0*I_ERI_G4y_S_Pz_Pz_C1001000003_a-3*I_ERI_D2y_S_Pz_Pz_C1001000003;
  abcd[317] = 2.0E0*I_ERI_G3yz_S_Pz_Pz_C1001000003_a-2*I_ERI_Dyz_S_Pz_Pz_C1001000003;
  abcd[318] = 2.0E0*I_ERI_G2y2z_S_Pz_Pz_C1001000003_a-1*I_ERI_D2z_S_Pz_Pz_C1001000003;
  abcd[319] = 2.0E0*I_ERI_Gy3z_S_Pz_Pz_C1001000003_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C3_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C3
   ************************************************************/
  abcd[320] = 2.0E0*I_ERI_G3xz_S_S_S_C3_a;
  abcd[321] = 2.0E0*I_ERI_G2xyz_S_S_S_C3_a;
  abcd[322] = 2.0E0*I_ERI_G2x2z_S_S_S_C3_a-1*I_ERI_D2x_S_S_S_C3;
  abcd[323] = 2.0E0*I_ERI_Gx2yz_S_S_S_C3_a;
  abcd[324] = 2.0E0*I_ERI_Gxy2z_S_S_S_C3_a-1*I_ERI_Dxy_S_S_S_C3;
  abcd[325] = 2.0E0*I_ERI_Gx3z_S_S_S_C3_a-2*I_ERI_Dxz_S_S_S_C3;
  abcd[326] = 2.0E0*I_ERI_G3yz_S_S_S_C3_a;
  abcd[327] = 2.0E0*I_ERI_G2y2z_S_S_S_C3_a-1*I_ERI_D2y_S_S_S_C3;
  abcd[328] = 2.0E0*I_ERI_Gy3z_S_S_S_C3_a-2*I_ERI_Dyz_S_S_S_C3;
  abcd[329] = 2.0E0*I_ERI_G4z_S_S_S_C3_a-3*I_ERI_D2z_S_S_S_C3;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000003
   ************************************************************/
  abcd[330] = 2.0E0*I_ERI_G3xz_S_Px_S_C1000003_a;
  abcd[331] = 2.0E0*I_ERI_G2xyz_S_Px_S_C1000003_a;
  abcd[332] = 2.0E0*I_ERI_G2x2z_S_Px_S_C1000003_a-1*I_ERI_D2x_S_Px_S_C1000003;
  abcd[333] = 2.0E0*I_ERI_Gx2yz_S_Px_S_C1000003_a;
  abcd[334] = 2.0E0*I_ERI_Gxy2z_S_Px_S_C1000003_a-1*I_ERI_Dxy_S_Px_S_C1000003;
  abcd[335] = 2.0E0*I_ERI_Gx3z_S_Px_S_C1000003_a-2*I_ERI_Dxz_S_Px_S_C1000003;
  abcd[336] = 2.0E0*I_ERI_G3yz_S_Px_S_C1000003_a;
  abcd[337] = 2.0E0*I_ERI_G2y2z_S_Px_S_C1000003_a-1*I_ERI_D2y_S_Px_S_C1000003;
  abcd[338] = 2.0E0*I_ERI_Gy3z_S_Px_S_C1000003_a-2*I_ERI_Dyz_S_Px_S_C1000003;
  abcd[339] = 2.0E0*I_ERI_G4z_S_Px_S_C1000003_a-3*I_ERI_D2z_S_Px_S_C1000003;
  abcd[340] = 2.0E0*I_ERI_G3xz_S_Py_S_C1000003_a;
  abcd[341] = 2.0E0*I_ERI_G2xyz_S_Py_S_C1000003_a;
  abcd[342] = 2.0E0*I_ERI_G2x2z_S_Py_S_C1000003_a-1*I_ERI_D2x_S_Py_S_C1000003;
  abcd[343] = 2.0E0*I_ERI_Gx2yz_S_Py_S_C1000003_a;
  abcd[344] = 2.0E0*I_ERI_Gxy2z_S_Py_S_C1000003_a-1*I_ERI_Dxy_S_Py_S_C1000003;
  abcd[345] = 2.0E0*I_ERI_Gx3z_S_Py_S_C1000003_a-2*I_ERI_Dxz_S_Py_S_C1000003;
  abcd[346] = 2.0E0*I_ERI_G3yz_S_Py_S_C1000003_a;
  abcd[347] = 2.0E0*I_ERI_G2y2z_S_Py_S_C1000003_a-1*I_ERI_D2y_S_Py_S_C1000003;
  abcd[348] = 2.0E0*I_ERI_Gy3z_S_Py_S_C1000003_a-2*I_ERI_Dyz_S_Py_S_C1000003;
  abcd[349] = 2.0E0*I_ERI_G4z_S_Py_S_C1000003_a-3*I_ERI_D2z_S_Py_S_C1000003;
  abcd[350] = 2.0E0*I_ERI_G3xz_S_Pz_S_C1000003_a;
  abcd[351] = 2.0E0*I_ERI_G2xyz_S_Pz_S_C1000003_a;
  abcd[352] = 2.0E0*I_ERI_G2x2z_S_Pz_S_C1000003_a-1*I_ERI_D2x_S_Pz_S_C1000003;
  abcd[353] = 2.0E0*I_ERI_Gx2yz_S_Pz_S_C1000003_a;
  abcd[354] = 2.0E0*I_ERI_Gxy2z_S_Pz_S_C1000003_a-1*I_ERI_Dxy_S_Pz_S_C1000003;
  abcd[355] = 2.0E0*I_ERI_Gx3z_S_Pz_S_C1000003_a-2*I_ERI_Dxz_S_Pz_S_C1000003;
  abcd[356] = 2.0E0*I_ERI_G3yz_S_Pz_S_C1000003_a;
  abcd[357] = 2.0E0*I_ERI_G2y2z_S_Pz_S_C1000003_a-1*I_ERI_D2y_S_Pz_S_C1000003;
  abcd[358] = 2.0E0*I_ERI_Gy3z_S_Pz_S_C1000003_a-2*I_ERI_Dyz_S_Pz_S_C1000003;
  abcd[359] = 2.0E0*I_ERI_G4z_S_Pz_S_C1000003_a-3*I_ERI_D2z_S_Pz_S_C1000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_P_C1000000003_a
   * RHS shell quartet name: SQ_ERI_D_S_S_P_C1000000003
   ************************************************************/
  abcd[360] = 2.0E0*I_ERI_G3xz_S_S_Px_C1000000003_a;
  abcd[361] = 2.0E0*I_ERI_G2xyz_S_S_Px_C1000000003_a;
  abcd[362] = 2.0E0*I_ERI_G2x2z_S_S_Px_C1000000003_a-1*I_ERI_D2x_S_S_Px_C1000000003;
  abcd[363] = 2.0E0*I_ERI_Gx2yz_S_S_Px_C1000000003_a;
  abcd[364] = 2.0E0*I_ERI_Gxy2z_S_S_Px_C1000000003_a-1*I_ERI_Dxy_S_S_Px_C1000000003;
  abcd[365] = 2.0E0*I_ERI_Gx3z_S_S_Px_C1000000003_a-2*I_ERI_Dxz_S_S_Px_C1000000003;
  abcd[366] = 2.0E0*I_ERI_G3yz_S_S_Px_C1000000003_a;
  abcd[367] = 2.0E0*I_ERI_G2y2z_S_S_Px_C1000000003_a-1*I_ERI_D2y_S_S_Px_C1000000003;
  abcd[368] = 2.0E0*I_ERI_Gy3z_S_S_Px_C1000000003_a-2*I_ERI_Dyz_S_S_Px_C1000000003;
  abcd[369] = 2.0E0*I_ERI_G4z_S_S_Px_C1000000003_a-3*I_ERI_D2z_S_S_Px_C1000000003;
  abcd[400] = 2.0E0*I_ERI_G3xz_S_S_Py_C1000000003_a;
  abcd[401] = 2.0E0*I_ERI_G2xyz_S_S_Py_C1000000003_a;
  abcd[402] = 2.0E0*I_ERI_G2x2z_S_S_Py_C1000000003_a-1*I_ERI_D2x_S_S_Py_C1000000003;
  abcd[403] = 2.0E0*I_ERI_Gx2yz_S_S_Py_C1000000003_a;
  abcd[404] = 2.0E0*I_ERI_Gxy2z_S_S_Py_C1000000003_a-1*I_ERI_Dxy_S_S_Py_C1000000003;
  abcd[405] = 2.0E0*I_ERI_Gx3z_S_S_Py_C1000000003_a-2*I_ERI_Dxz_S_S_Py_C1000000003;
  abcd[406] = 2.0E0*I_ERI_G3yz_S_S_Py_C1000000003_a;
  abcd[407] = 2.0E0*I_ERI_G2y2z_S_S_Py_C1000000003_a-1*I_ERI_D2y_S_S_Py_C1000000003;
  abcd[408] = 2.0E0*I_ERI_Gy3z_S_S_Py_C1000000003_a-2*I_ERI_Dyz_S_S_Py_C1000000003;
  abcd[409] = 2.0E0*I_ERI_G4z_S_S_Py_C1000000003_a-3*I_ERI_D2z_S_S_Py_C1000000003;
  abcd[440] = 2.0E0*I_ERI_G3xz_S_S_Pz_C1000000003_a;
  abcd[441] = 2.0E0*I_ERI_G2xyz_S_S_Pz_C1000000003_a;
  abcd[442] = 2.0E0*I_ERI_G2x2z_S_S_Pz_C1000000003_a-1*I_ERI_D2x_S_S_Pz_C1000000003;
  abcd[443] = 2.0E0*I_ERI_Gx2yz_S_S_Pz_C1000000003_a;
  abcd[444] = 2.0E0*I_ERI_Gxy2z_S_S_Pz_C1000000003_a-1*I_ERI_Dxy_S_S_Pz_C1000000003;
  abcd[445] = 2.0E0*I_ERI_Gx3z_S_S_Pz_C1000000003_a-2*I_ERI_Dxz_S_S_Pz_C1000000003;
  abcd[446] = 2.0E0*I_ERI_G3yz_S_S_Pz_C1000000003_a;
  abcd[447] = 2.0E0*I_ERI_G2y2z_S_S_Pz_C1000000003_a-1*I_ERI_D2y_S_S_Pz_C1000000003;
  abcd[448] = 2.0E0*I_ERI_Gy3z_S_S_Pz_C1000000003_a-2*I_ERI_Dyz_S_S_Pz_C1000000003;
  abcd[449] = 2.0E0*I_ERI_G4z_S_S_Pz_C1000000003_a-3*I_ERI_D2z_S_S_Pz_C1000000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P_C1001000003_a
   * RHS shell quartet name: SQ_ERI_D_S_P_P_C1001000003
   ************************************************************/
  abcd[370] = 2.0E0*I_ERI_G3xz_S_Px_Px_C1001000003_a;
  abcd[371] = 2.0E0*I_ERI_G2xyz_S_Px_Px_C1001000003_a;
  abcd[372] = 2.0E0*I_ERI_G2x2z_S_Px_Px_C1001000003_a-1*I_ERI_D2x_S_Px_Px_C1001000003;
  abcd[373] = 2.0E0*I_ERI_Gx2yz_S_Px_Px_C1001000003_a;
  abcd[374] = 2.0E0*I_ERI_Gxy2z_S_Px_Px_C1001000003_a-1*I_ERI_Dxy_S_Px_Px_C1001000003;
  abcd[375] = 2.0E0*I_ERI_Gx3z_S_Px_Px_C1001000003_a-2*I_ERI_Dxz_S_Px_Px_C1001000003;
  abcd[376] = 2.0E0*I_ERI_G3yz_S_Px_Px_C1001000003_a;
  abcd[377] = 2.0E0*I_ERI_G2y2z_S_Px_Px_C1001000003_a-1*I_ERI_D2y_S_Px_Px_C1001000003;
  abcd[378] = 2.0E0*I_ERI_Gy3z_S_Px_Px_C1001000003_a-2*I_ERI_Dyz_S_Px_Px_C1001000003;
  abcd[379] = 2.0E0*I_ERI_G4z_S_Px_Px_C1001000003_a-3*I_ERI_D2z_S_Px_Px_C1001000003;
  abcd[380] = 2.0E0*I_ERI_G3xz_S_Py_Px_C1001000003_a;
  abcd[381] = 2.0E0*I_ERI_G2xyz_S_Py_Px_C1001000003_a;
  abcd[382] = 2.0E0*I_ERI_G2x2z_S_Py_Px_C1001000003_a-1*I_ERI_D2x_S_Py_Px_C1001000003;
  abcd[383] = 2.0E0*I_ERI_Gx2yz_S_Py_Px_C1001000003_a;
  abcd[384] = 2.0E0*I_ERI_Gxy2z_S_Py_Px_C1001000003_a-1*I_ERI_Dxy_S_Py_Px_C1001000003;
  abcd[385] = 2.0E0*I_ERI_Gx3z_S_Py_Px_C1001000003_a-2*I_ERI_Dxz_S_Py_Px_C1001000003;
  abcd[386] = 2.0E0*I_ERI_G3yz_S_Py_Px_C1001000003_a;
  abcd[387] = 2.0E0*I_ERI_G2y2z_S_Py_Px_C1001000003_a-1*I_ERI_D2y_S_Py_Px_C1001000003;
  abcd[388] = 2.0E0*I_ERI_Gy3z_S_Py_Px_C1001000003_a-2*I_ERI_Dyz_S_Py_Px_C1001000003;
  abcd[389] = 2.0E0*I_ERI_G4z_S_Py_Px_C1001000003_a-3*I_ERI_D2z_S_Py_Px_C1001000003;
  abcd[390] = 2.0E0*I_ERI_G3xz_S_Pz_Px_C1001000003_a;
  abcd[391] = 2.0E0*I_ERI_G2xyz_S_Pz_Px_C1001000003_a;
  abcd[392] = 2.0E0*I_ERI_G2x2z_S_Pz_Px_C1001000003_a-1*I_ERI_D2x_S_Pz_Px_C1001000003;
  abcd[393] = 2.0E0*I_ERI_Gx2yz_S_Pz_Px_C1001000003_a;
  abcd[394] = 2.0E0*I_ERI_Gxy2z_S_Pz_Px_C1001000003_a-1*I_ERI_Dxy_S_Pz_Px_C1001000003;
  abcd[395] = 2.0E0*I_ERI_Gx3z_S_Pz_Px_C1001000003_a-2*I_ERI_Dxz_S_Pz_Px_C1001000003;
  abcd[396] = 2.0E0*I_ERI_G3yz_S_Pz_Px_C1001000003_a;
  abcd[397] = 2.0E0*I_ERI_G2y2z_S_Pz_Px_C1001000003_a-1*I_ERI_D2y_S_Pz_Px_C1001000003;
  abcd[398] = 2.0E0*I_ERI_Gy3z_S_Pz_Px_C1001000003_a-2*I_ERI_Dyz_S_Pz_Px_C1001000003;
  abcd[399] = 2.0E0*I_ERI_G4z_S_Pz_Px_C1001000003_a-3*I_ERI_D2z_S_Pz_Px_C1001000003;
  abcd[410] = 2.0E0*I_ERI_G3xz_S_Px_Py_C1001000003_a;
  abcd[411] = 2.0E0*I_ERI_G2xyz_S_Px_Py_C1001000003_a;
  abcd[412] = 2.0E0*I_ERI_G2x2z_S_Px_Py_C1001000003_a-1*I_ERI_D2x_S_Px_Py_C1001000003;
  abcd[413] = 2.0E0*I_ERI_Gx2yz_S_Px_Py_C1001000003_a;
  abcd[414] = 2.0E0*I_ERI_Gxy2z_S_Px_Py_C1001000003_a-1*I_ERI_Dxy_S_Px_Py_C1001000003;
  abcd[415] = 2.0E0*I_ERI_Gx3z_S_Px_Py_C1001000003_a-2*I_ERI_Dxz_S_Px_Py_C1001000003;
  abcd[416] = 2.0E0*I_ERI_G3yz_S_Px_Py_C1001000003_a;
  abcd[417] = 2.0E0*I_ERI_G2y2z_S_Px_Py_C1001000003_a-1*I_ERI_D2y_S_Px_Py_C1001000003;
  abcd[418] = 2.0E0*I_ERI_Gy3z_S_Px_Py_C1001000003_a-2*I_ERI_Dyz_S_Px_Py_C1001000003;
  abcd[419] = 2.0E0*I_ERI_G4z_S_Px_Py_C1001000003_a-3*I_ERI_D2z_S_Px_Py_C1001000003;
  abcd[420] = 2.0E0*I_ERI_G3xz_S_Py_Py_C1001000003_a;
  abcd[421] = 2.0E0*I_ERI_G2xyz_S_Py_Py_C1001000003_a;
  abcd[422] = 2.0E0*I_ERI_G2x2z_S_Py_Py_C1001000003_a-1*I_ERI_D2x_S_Py_Py_C1001000003;
  abcd[423] = 2.0E0*I_ERI_Gx2yz_S_Py_Py_C1001000003_a;
  abcd[424] = 2.0E0*I_ERI_Gxy2z_S_Py_Py_C1001000003_a-1*I_ERI_Dxy_S_Py_Py_C1001000003;
  abcd[425] = 2.0E0*I_ERI_Gx3z_S_Py_Py_C1001000003_a-2*I_ERI_Dxz_S_Py_Py_C1001000003;
  abcd[426] = 2.0E0*I_ERI_G3yz_S_Py_Py_C1001000003_a;
  abcd[427] = 2.0E0*I_ERI_G2y2z_S_Py_Py_C1001000003_a-1*I_ERI_D2y_S_Py_Py_C1001000003;
  abcd[428] = 2.0E0*I_ERI_Gy3z_S_Py_Py_C1001000003_a-2*I_ERI_Dyz_S_Py_Py_C1001000003;
  abcd[429] = 2.0E0*I_ERI_G4z_S_Py_Py_C1001000003_a-3*I_ERI_D2z_S_Py_Py_C1001000003;
  abcd[430] = 2.0E0*I_ERI_G3xz_S_Pz_Py_C1001000003_a;
  abcd[431] = 2.0E0*I_ERI_G2xyz_S_Pz_Py_C1001000003_a;
  abcd[432] = 2.0E0*I_ERI_G2x2z_S_Pz_Py_C1001000003_a-1*I_ERI_D2x_S_Pz_Py_C1001000003;
  abcd[433] = 2.0E0*I_ERI_Gx2yz_S_Pz_Py_C1001000003_a;
  abcd[434] = 2.0E0*I_ERI_Gxy2z_S_Pz_Py_C1001000003_a-1*I_ERI_Dxy_S_Pz_Py_C1001000003;
  abcd[435] = 2.0E0*I_ERI_Gx3z_S_Pz_Py_C1001000003_a-2*I_ERI_Dxz_S_Pz_Py_C1001000003;
  abcd[436] = 2.0E0*I_ERI_G3yz_S_Pz_Py_C1001000003_a;
  abcd[437] = 2.0E0*I_ERI_G2y2z_S_Pz_Py_C1001000003_a-1*I_ERI_D2y_S_Pz_Py_C1001000003;
  abcd[438] = 2.0E0*I_ERI_Gy3z_S_Pz_Py_C1001000003_a-2*I_ERI_Dyz_S_Pz_Py_C1001000003;
  abcd[439] = 2.0E0*I_ERI_G4z_S_Pz_Py_C1001000003_a-3*I_ERI_D2z_S_Pz_Py_C1001000003;
  abcd[450] = 2.0E0*I_ERI_G3xz_S_Px_Pz_C1001000003_a;
  abcd[451] = 2.0E0*I_ERI_G2xyz_S_Px_Pz_C1001000003_a;
  abcd[452] = 2.0E0*I_ERI_G2x2z_S_Px_Pz_C1001000003_a-1*I_ERI_D2x_S_Px_Pz_C1001000003;
  abcd[453] = 2.0E0*I_ERI_Gx2yz_S_Px_Pz_C1001000003_a;
  abcd[454] = 2.0E0*I_ERI_Gxy2z_S_Px_Pz_C1001000003_a-1*I_ERI_Dxy_S_Px_Pz_C1001000003;
  abcd[455] = 2.0E0*I_ERI_Gx3z_S_Px_Pz_C1001000003_a-2*I_ERI_Dxz_S_Px_Pz_C1001000003;
  abcd[456] = 2.0E0*I_ERI_G3yz_S_Px_Pz_C1001000003_a;
  abcd[457] = 2.0E0*I_ERI_G2y2z_S_Px_Pz_C1001000003_a-1*I_ERI_D2y_S_Px_Pz_C1001000003;
  abcd[458] = 2.0E0*I_ERI_Gy3z_S_Px_Pz_C1001000003_a-2*I_ERI_Dyz_S_Px_Pz_C1001000003;
  abcd[459] = 2.0E0*I_ERI_G4z_S_Px_Pz_C1001000003_a-3*I_ERI_D2z_S_Px_Pz_C1001000003;
  abcd[460] = 2.0E0*I_ERI_G3xz_S_Py_Pz_C1001000003_a;
  abcd[461] = 2.0E0*I_ERI_G2xyz_S_Py_Pz_C1001000003_a;
  abcd[462] = 2.0E0*I_ERI_G2x2z_S_Py_Pz_C1001000003_a-1*I_ERI_D2x_S_Py_Pz_C1001000003;
  abcd[463] = 2.0E0*I_ERI_Gx2yz_S_Py_Pz_C1001000003_a;
  abcd[464] = 2.0E0*I_ERI_Gxy2z_S_Py_Pz_C1001000003_a-1*I_ERI_Dxy_S_Py_Pz_C1001000003;
  abcd[465] = 2.0E0*I_ERI_Gx3z_S_Py_Pz_C1001000003_a-2*I_ERI_Dxz_S_Py_Pz_C1001000003;
  abcd[466] = 2.0E0*I_ERI_G3yz_S_Py_Pz_C1001000003_a;
  abcd[467] = 2.0E0*I_ERI_G2y2z_S_Py_Pz_C1001000003_a-1*I_ERI_D2y_S_Py_Pz_C1001000003;
  abcd[468] = 2.0E0*I_ERI_Gy3z_S_Py_Pz_C1001000003_a-2*I_ERI_Dyz_S_Py_Pz_C1001000003;
  abcd[469] = 2.0E0*I_ERI_G4z_S_Py_Pz_C1001000003_a-3*I_ERI_D2z_S_Py_Pz_C1001000003;
  abcd[470] = 2.0E0*I_ERI_G3xz_S_Pz_Pz_C1001000003_a;
  abcd[471] = 2.0E0*I_ERI_G2xyz_S_Pz_Pz_C1001000003_a;
  abcd[472] = 2.0E0*I_ERI_G2x2z_S_Pz_Pz_C1001000003_a-1*I_ERI_D2x_S_Pz_Pz_C1001000003;
  abcd[473] = 2.0E0*I_ERI_Gx2yz_S_Pz_Pz_C1001000003_a;
  abcd[474] = 2.0E0*I_ERI_Gxy2z_S_Pz_Pz_C1001000003_a-1*I_ERI_Dxy_S_Pz_Pz_C1001000003;
  abcd[475] = 2.0E0*I_ERI_Gx3z_S_Pz_Pz_C1001000003_a-2*I_ERI_Dxz_S_Pz_Pz_C1001000003;
  abcd[476] = 2.0E0*I_ERI_G3yz_S_Pz_Pz_C1001000003_a;
  abcd[477] = 2.0E0*I_ERI_G2y2z_S_Pz_Pz_C1001000003_a-1*I_ERI_D2y_S_Pz_Pz_C1001000003;
  abcd[478] = 2.0E0*I_ERI_Gy3z_S_Pz_Pz_C1001000003_a-2*I_ERI_Dyz_S_Pz_Pz_C1001000003;
  abcd[479] = 2.0E0*I_ERI_G4z_S_Pz_Pz_C1001000003_a-3*I_ERI_D2z_S_Pz_Pz_C1001000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C3_b
   ************************************************************/
  abcd[480] = 2.0E0*I_ERI_F3x_Px_S_S_C3_b;
  abcd[481] = 2.0E0*I_ERI_F2xy_Px_S_S_C3_b;
  abcd[482] = 2.0E0*I_ERI_F2xz_Px_S_S_C3_b;
  abcd[483] = 2.0E0*I_ERI_Fx2y_Px_S_S_C3_b;
  abcd[484] = 2.0E0*I_ERI_Fxyz_Px_S_S_C3_b;
  abcd[485] = 2.0E0*I_ERI_Fx2z_Px_S_S_C3_b;
  abcd[486] = 2.0E0*I_ERI_F3y_Px_S_S_C3_b;
  abcd[487] = 2.0E0*I_ERI_F2yz_Px_S_S_C3_b;
  abcd[488] = 2.0E0*I_ERI_Fy2z_Px_S_S_C3_b;
  abcd[489] = 2.0E0*I_ERI_F3z_Px_S_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[490] = 2.0E0*I_ERI_F3x_Px_Px_S_C1000003_b;
  abcd[491] = 2.0E0*I_ERI_F2xy_Px_Px_S_C1000003_b;
  abcd[492] = 2.0E0*I_ERI_F2xz_Px_Px_S_C1000003_b;
  abcd[493] = 2.0E0*I_ERI_Fx2y_Px_Px_S_C1000003_b;
  abcd[494] = 2.0E0*I_ERI_Fxyz_Px_Px_S_C1000003_b;
  abcd[495] = 2.0E0*I_ERI_Fx2z_Px_Px_S_C1000003_b;
  abcd[496] = 2.0E0*I_ERI_F3y_Px_Px_S_C1000003_b;
  abcd[497] = 2.0E0*I_ERI_F2yz_Px_Px_S_C1000003_b;
  abcd[498] = 2.0E0*I_ERI_Fy2z_Px_Px_S_C1000003_b;
  abcd[499] = 2.0E0*I_ERI_F3z_Px_Px_S_C1000003_b;
  abcd[500] = 2.0E0*I_ERI_F3x_Px_Py_S_C1000003_b;
  abcd[501] = 2.0E0*I_ERI_F2xy_Px_Py_S_C1000003_b;
  abcd[502] = 2.0E0*I_ERI_F2xz_Px_Py_S_C1000003_b;
  abcd[503] = 2.0E0*I_ERI_Fx2y_Px_Py_S_C1000003_b;
  abcd[504] = 2.0E0*I_ERI_Fxyz_Px_Py_S_C1000003_b;
  abcd[505] = 2.0E0*I_ERI_Fx2z_Px_Py_S_C1000003_b;
  abcd[506] = 2.0E0*I_ERI_F3y_Px_Py_S_C1000003_b;
  abcd[507] = 2.0E0*I_ERI_F2yz_Px_Py_S_C1000003_b;
  abcd[508] = 2.0E0*I_ERI_Fy2z_Px_Py_S_C1000003_b;
  abcd[509] = 2.0E0*I_ERI_F3z_Px_Py_S_C1000003_b;
  abcd[510] = 2.0E0*I_ERI_F3x_Px_Pz_S_C1000003_b;
  abcd[511] = 2.0E0*I_ERI_F2xy_Px_Pz_S_C1000003_b;
  abcd[512] = 2.0E0*I_ERI_F2xz_Px_Pz_S_C1000003_b;
  abcd[513] = 2.0E0*I_ERI_Fx2y_Px_Pz_S_C1000003_b;
  abcd[514] = 2.0E0*I_ERI_Fxyz_Px_Pz_S_C1000003_b;
  abcd[515] = 2.0E0*I_ERI_Fx2z_Px_Pz_S_C1000003_b;
  abcd[516] = 2.0E0*I_ERI_F3y_Px_Pz_S_C1000003_b;
  abcd[517] = 2.0E0*I_ERI_F2yz_Px_Pz_S_C1000003_b;
  abcd[518] = 2.0E0*I_ERI_Fy2z_Px_Pz_S_C1000003_b;
  abcd[519] = 2.0E0*I_ERI_F3z_Px_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_P_C1000000003_b
   ************************************************************/
  abcd[520] = 2.0E0*I_ERI_F3x_Px_S_Px_C1000000003_b;
  abcd[521] = 2.0E0*I_ERI_F2xy_Px_S_Px_C1000000003_b;
  abcd[522] = 2.0E0*I_ERI_F2xz_Px_S_Px_C1000000003_b;
  abcd[523] = 2.0E0*I_ERI_Fx2y_Px_S_Px_C1000000003_b;
  abcd[524] = 2.0E0*I_ERI_Fxyz_Px_S_Px_C1000000003_b;
  abcd[525] = 2.0E0*I_ERI_Fx2z_Px_S_Px_C1000000003_b;
  abcd[526] = 2.0E0*I_ERI_F3y_Px_S_Px_C1000000003_b;
  abcd[527] = 2.0E0*I_ERI_F2yz_Px_S_Px_C1000000003_b;
  abcd[528] = 2.0E0*I_ERI_Fy2z_Px_S_Px_C1000000003_b;
  abcd[529] = 2.0E0*I_ERI_F3z_Px_S_Px_C1000000003_b;
  abcd[560] = 2.0E0*I_ERI_F3x_Px_S_Py_C1000000003_b;
  abcd[561] = 2.0E0*I_ERI_F2xy_Px_S_Py_C1000000003_b;
  abcd[562] = 2.0E0*I_ERI_F2xz_Px_S_Py_C1000000003_b;
  abcd[563] = 2.0E0*I_ERI_Fx2y_Px_S_Py_C1000000003_b;
  abcd[564] = 2.0E0*I_ERI_Fxyz_Px_S_Py_C1000000003_b;
  abcd[565] = 2.0E0*I_ERI_Fx2z_Px_S_Py_C1000000003_b;
  abcd[566] = 2.0E0*I_ERI_F3y_Px_S_Py_C1000000003_b;
  abcd[567] = 2.0E0*I_ERI_F2yz_Px_S_Py_C1000000003_b;
  abcd[568] = 2.0E0*I_ERI_Fy2z_Px_S_Py_C1000000003_b;
  abcd[569] = 2.0E0*I_ERI_F3z_Px_S_Py_C1000000003_b;
  abcd[600] = 2.0E0*I_ERI_F3x_Px_S_Pz_C1000000003_b;
  abcd[601] = 2.0E0*I_ERI_F2xy_Px_S_Pz_C1000000003_b;
  abcd[602] = 2.0E0*I_ERI_F2xz_Px_S_Pz_C1000000003_b;
  abcd[603] = 2.0E0*I_ERI_Fx2y_Px_S_Pz_C1000000003_b;
  abcd[604] = 2.0E0*I_ERI_Fxyz_Px_S_Pz_C1000000003_b;
  abcd[605] = 2.0E0*I_ERI_Fx2z_Px_S_Pz_C1000000003_b;
  abcd[606] = 2.0E0*I_ERI_F3y_Px_S_Pz_C1000000003_b;
  abcd[607] = 2.0E0*I_ERI_F2yz_Px_S_Pz_C1000000003_b;
  abcd[608] = 2.0E0*I_ERI_Fy2z_Px_S_Pz_C1000000003_b;
  abcd[609] = 2.0E0*I_ERI_F3z_Px_S_Pz_C1000000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[530] = 2.0E0*I_ERI_F3x_Px_Px_Px_C1001000003_b;
  abcd[531] = 2.0E0*I_ERI_F2xy_Px_Px_Px_C1001000003_b;
  abcd[532] = 2.0E0*I_ERI_F2xz_Px_Px_Px_C1001000003_b;
  abcd[533] = 2.0E0*I_ERI_Fx2y_Px_Px_Px_C1001000003_b;
  abcd[534] = 2.0E0*I_ERI_Fxyz_Px_Px_Px_C1001000003_b;
  abcd[535] = 2.0E0*I_ERI_Fx2z_Px_Px_Px_C1001000003_b;
  abcd[536] = 2.0E0*I_ERI_F3y_Px_Px_Px_C1001000003_b;
  abcd[537] = 2.0E0*I_ERI_F2yz_Px_Px_Px_C1001000003_b;
  abcd[538] = 2.0E0*I_ERI_Fy2z_Px_Px_Px_C1001000003_b;
  abcd[539] = 2.0E0*I_ERI_F3z_Px_Px_Px_C1001000003_b;
  abcd[540] = 2.0E0*I_ERI_F3x_Px_Py_Px_C1001000003_b;
  abcd[541] = 2.0E0*I_ERI_F2xy_Px_Py_Px_C1001000003_b;
  abcd[542] = 2.0E0*I_ERI_F2xz_Px_Py_Px_C1001000003_b;
  abcd[543] = 2.0E0*I_ERI_Fx2y_Px_Py_Px_C1001000003_b;
  abcd[544] = 2.0E0*I_ERI_Fxyz_Px_Py_Px_C1001000003_b;
  abcd[545] = 2.0E0*I_ERI_Fx2z_Px_Py_Px_C1001000003_b;
  abcd[546] = 2.0E0*I_ERI_F3y_Px_Py_Px_C1001000003_b;
  abcd[547] = 2.0E0*I_ERI_F2yz_Px_Py_Px_C1001000003_b;
  abcd[548] = 2.0E0*I_ERI_Fy2z_Px_Py_Px_C1001000003_b;
  abcd[549] = 2.0E0*I_ERI_F3z_Px_Py_Px_C1001000003_b;
  abcd[550] = 2.0E0*I_ERI_F3x_Px_Pz_Px_C1001000003_b;
  abcd[551] = 2.0E0*I_ERI_F2xy_Px_Pz_Px_C1001000003_b;
  abcd[552] = 2.0E0*I_ERI_F2xz_Px_Pz_Px_C1001000003_b;
  abcd[553] = 2.0E0*I_ERI_Fx2y_Px_Pz_Px_C1001000003_b;
  abcd[554] = 2.0E0*I_ERI_Fxyz_Px_Pz_Px_C1001000003_b;
  abcd[555] = 2.0E0*I_ERI_Fx2z_Px_Pz_Px_C1001000003_b;
  abcd[556] = 2.0E0*I_ERI_F3y_Px_Pz_Px_C1001000003_b;
  abcd[557] = 2.0E0*I_ERI_F2yz_Px_Pz_Px_C1001000003_b;
  abcd[558] = 2.0E0*I_ERI_Fy2z_Px_Pz_Px_C1001000003_b;
  abcd[559] = 2.0E0*I_ERI_F3z_Px_Pz_Px_C1001000003_b;
  abcd[570] = 2.0E0*I_ERI_F3x_Px_Px_Py_C1001000003_b;
  abcd[571] = 2.0E0*I_ERI_F2xy_Px_Px_Py_C1001000003_b;
  abcd[572] = 2.0E0*I_ERI_F2xz_Px_Px_Py_C1001000003_b;
  abcd[573] = 2.0E0*I_ERI_Fx2y_Px_Px_Py_C1001000003_b;
  abcd[574] = 2.0E0*I_ERI_Fxyz_Px_Px_Py_C1001000003_b;
  abcd[575] = 2.0E0*I_ERI_Fx2z_Px_Px_Py_C1001000003_b;
  abcd[576] = 2.0E0*I_ERI_F3y_Px_Px_Py_C1001000003_b;
  abcd[577] = 2.0E0*I_ERI_F2yz_Px_Px_Py_C1001000003_b;
  abcd[578] = 2.0E0*I_ERI_Fy2z_Px_Px_Py_C1001000003_b;
  abcd[579] = 2.0E0*I_ERI_F3z_Px_Px_Py_C1001000003_b;
  abcd[580] = 2.0E0*I_ERI_F3x_Px_Py_Py_C1001000003_b;
  abcd[581] = 2.0E0*I_ERI_F2xy_Px_Py_Py_C1001000003_b;
  abcd[582] = 2.0E0*I_ERI_F2xz_Px_Py_Py_C1001000003_b;
  abcd[583] = 2.0E0*I_ERI_Fx2y_Px_Py_Py_C1001000003_b;
  abcd[584] = 2.0E0*I_ERI_Fxyz_Px_Py_Py_C1001000003_b;
  abcd[585] = 2.0E0*I_ERI_Fx2z_Px_Py_Py_C1001000003_b;
  abcd[586] = 2.0E0*I_ERI_F3y_Px_Py_Py_C1001000003_b;
  abcd[587] = 2.0E0*I_ERI_F2yz_Px_Py_Py_C1001000003_b;
  abcd[588] = 2.0E0*I_ERI_Fy2z_Px_Py_Py_C1001000003_b;
  abcd[589] = 2.0E0*I_ERI_F3z_Px_Py_Py_C1001000003_b;
  abcd[590] = 2.0E0*I_ERI_F3x_Px_Pz_Py_C1001000003_b;
  abcd[591] = 2.0E0*I_ERI_F2xy_Px_Pz_Py_C1001000003_b;
  abcd[592] = 2.0E0*I_ERI_F2xz_Px_Pz_Py_C1001000003_b;
  abcd[593] = 2.0E0*I_ERI_Fx2y_Px_Pz_Py_C1001000003_b;
  abcd[594] = 2.0E0*I_ERI_Fxyz_Px_Pz_Py_C1001000003_b;
  abcd[595] = 2.0E0*I_ERI_Fx2z_Px_Pz_Py_C1001000003_b;
  abcd[596] = 2.0E0*I_ERI_F3y_Px_Pz_Py_C1001000003_b;
  abcd[597] = 2.0E0*I_ERI_F2yz_Px_Pz_Py_C1001000003_b;
  abcd[598] = 2.0E0*I_ERI_Fy2z_Px_Pz_Py_C1001000003_b;
  abcd[599] = 2.0E0*I_ERI_F3z_Px_Pz_Py_C1001000003_b;
  abcd[610] = 2.0E0*I_ERI_F3x_Px_Px_Pz_C1001000003_b;
  abcd[611] = 2.0E0*I_ERI_F2xy_Px_Px_Pz_C1001000003_b;
  abcd[612] = 2.0E0*I_ERI_F2xz_Px_Px_Pz_C1001000003_b;
  abcd[613] = 2.0E0*I_ERI_Fx2y_Px_Px_Pz_C1001000003_b;
  abcd[614] = 2.0E0*I_ERI_Fxyz_Px_Px_Pz_C1001000003_b;
  abcd[615] = 2.0E0*I_ERI_Fx2z_Px_Px_Pz_C1001000003_b;
  abcd[616] = 2.0E0*I_ERI_F3y_Px_Px_Pz_C1001000003_b;
  abcd[617] = 2.0E0*I_ERI_F2yz_Px_Px_Pz_C1001000003_b;
  abcd[618] = 2.0E0*I_ERI_Fy2z_Px_Px_Pz_C1001000003_b;
  abcd[619] = 2.0E0*I_ERI_F3z_Px_Px_Pz_C1001000003_b;
  abcd[620] = 2.0E0*I_ERI_F3x_Px_Py_Pz_C1001000003_b;
  abcd[621] = 2.0E0*I_ERI_F2xy_Px_Py_Pz_C1001000003_b;
  abcd[622] = 2.0E0*I_ERI_F2xz_Px_Py_Pz_C1001000003_b;
  abcd[623] = 2.0E0*I_ERI_Fx2y_Px_Py_Pz_C1001000003_b;
  abcd[624] = 2.0E0*I_ERI_Fxyz_Px_Py_Pz_C1001000003_b;
  abcd[625] = 2.0E0*I_ERI_Fx2z_Px_Py_Pz_C1001000003_b;
  abcd[626] = 2.0E0*I_ERI_F3y_Px_Py_Pz_C1001000003_b;
  abcd[627] = 2.0E0*I_ERI_F2yz_Px_Py_Pz_C1001000003_b;
  abcd[628] = 2.0E0*I_ERI_Fy2z_Px_Py_Pz_C1001000003_b;
  abcd[629] = 2.0E0*I_ERI_F3z_Px_Py_Pz_C1001000003_b;
  abcd[630] = 2.0E0*I_ERI_F3x_Px_Pz_Pz_C1001000003_b;
  abcd[631] = 2.0E0*I_ERI_F2xy_Px_Pz_Pz_C1001000003_b;
  abcd[632] = 2.0E0*I_ERI_F2xz_Px_Pz_Pz_C1001000003_b;
  abcd[633] = 2.0E0*I_ERI_Fx2y_Px_Pz_Pz_C1001000003_b;
  abcd[634] = 2.0E0*I_ERI_Fxyz_Px_Pz_Pz_C1001000003_b;
  abcd[635] = 2.0E0*I_ERI_Fx2z_Px_Pz_Pz_C1001000003_b;
  abcd[636] = 2.0E0*I_ERI_F3y_Px_Pz_Pz_C1001000003_b;
  abcd[637] = 2.0E0*I_ERI_F2yz_Px_Pz_Pz_C1001000003_b;
  abcd[638] = 2.0E0*I_ERI_Fy2z_Px_Pz_Pz_C1001000003_b;
  abcd[639] = 2.0E0*I_ERI_F3z_Px_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C3_b
   ************************************************************/
  abcd[640] = 2.0E0*I_ERI_F3x_Py_S_S_C3_b;
  abcd[641] = 2.0E0*I_ERI_F2xy_Py_S_S_C3_b;
  abcd[642] = 2.0E0*I_ERI_F2xz_Py_S_S_C3_b;
  abcd[643] = 2.0E0*I_ERI_Fx2y_Py_S_S_C3_b;
  abcd[644] = 2.0E0*I_ERI_Fxyz_Py_S_S_C3_b;
  abcd[645] = 2.0E0*I_ERI_Fx2z_Py_S_S_C3_b;
  abcd[646] = 2.0E0*I_ERI_F3y_Py_S_S_C3_b;
  abcd[647] = 2.0E0*I_ERI_F2yz_Py_S_S_C3_b;
  abcd[648] = 2.0E0*I_ERI_Fy2z_Py_S_S_C3_b;
  abcd[649] = 2.0E0*I_ERI_F3z_Py_S_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[650] = 2.0E0*I_ERI_F3x_Py_Px_S_C1000003_b;
  abcd[651] = 2.0E0*I_ERI_F2xy_Py_Px_S_C1000003_b;
  abcd[652] = 2.0E0*I_ERI_F2xz_Py_Px_S_C1000003_b;
  abcd[653] = 2.0E0*I_ERI_Fx2y_Py_Px_S_C1000003_b;
  abcd[654] = 2.0E0*I_ERI_Fxyz_Py_Px_S_C1000003_b;
  abcd[655] = 2.0E0*I_ERI_Fx2z_Py_Px_S_C1000003_b;
  abcd[656] = 2.0E0*I_ERI_F3y_Py_Px_S_C1000003_b;
  abcd[657] = 2.0E0*I_ERI_F2yz_Py_Px_S_C1000003_b;
  abcd[658] = 2.0E0*I_ERI_Fy2z_Py_Px_S_C1000003_b;
  abcd[659] = 2.0E0*I_ERI_F3z_Py_Px_S_C1000003_b;
  abcd[660] = 2.0E0*I_ERI_F3x_Py_Py_S_C1000003_b;
  abcd[661] = 2.0E0*I_ERI_F2xy_Py_Py_S_C1000003_b;
  abcd[662] = 2.0E0*I_ERI_F2xz_Py_Py_S_C1000003_b;
  abcd[663] = 2.0E0*I_ERI_Fx2y_Py_Py_S_C1000003_b;
  abcd[664] = 2.0E0*I_ERI_Fxyz_Py_Py_S_C1000003_b;
  abcd[665] = 2.0E0*I_ERI_Fx2z_Py_Py_S_C1000003_b;
  abcd[666] = 2.0E0*I_ERI_F3y_Py_Py_S_C1000003_b;
  abcd[667] = 2.0E0*I_ERI_F2yz_Py_Py_S_C1000003_b;
  abcd[668] = 2.0E0*I_ERI_Fy2z_Py_Py_S_C1000003_b;
  abcd[669] = 2.0E0*I_ERI_F3z_Py_Py_S_C1000003_b;
  abcd[670] = 2.0E0*I_ERI_F3x_Py_Pz_S_C1000003_b;
  abcd[671] = 2.0E0*I_ERI_F2xy_Py_Pz_S_C1000003_b;
  abcd[672] = 2.0E0*I_ERI_F2xz_Py_Pz_S_C1000003_b;
  abcd[673] = 2.0E0*I_ERI_Fx2y_Py_Pz_S_C1000003_b;
  abcd[674] = 2.0E0*I_ERI_Fxyz_Py_Pz_S_C1000003_b;
  abcd[675] = 2.0E0*I_ERI_Fx2z_Py_Pz_S_C1000003_b;
  abcd[676] = 2.0E0*I_ERI_F3y_Py_Pz_S_C1000003_b;
  abcd[677] = 2.0E0*I_ERI_F2yz_Py_Pz_S_C1000003_b;
  abcd[678] = 2.0E0*I_ERI_Fy2z_Py_Pz_S_C1000003_b;
  abcd[679] = 2.0E0*I_ERI_F3z_Py_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_P_C1000000003_b
   ************************************************************/
  abcd[680] = 2.0E0*I_ERI_F3x_Py_S_Px_C1000000003_b;
  abcd[681] = 2.0E0*I_ERI_F2xy_Py_S_Px_C1000000003_b;
  abcd[682] = 2.0E0*I_ERI_F2xz_Py_S_Px_C1000000003_b;
  abcd[683] = 2.0E0*I_ERI_Fx2y_Py_S_Px_C1000000003_b;
  abcd[684] = 2.0E0*I_ERI_Fxyz_Py_S_Px_C1000000003_b;
  abcd[685] = 2.0E0*I_ERI_Fx2z_Py_S_Px_C1000000003_b;
  abcd[686] = 2.0E0*I_ERI_F3y_Py_S_Px_C1000000003_b;
  abcd[687] = 2.0E0*I_ERI_F2yz_Py_S_Px_C1000000003_b;
  abcd[688] = 2.0E0*I_ERI_Fy2z_Py_S_Px_C1000000003_b;
  abcd[689] = 2.0E0*I_ERI_F3z_Py_S_Px_C1000000003_b;
  abcd[720] = 2.0E0*I_ERI_F3x_Py_S_Py_C1000000003_b;
  abcd[721] = 2.0E0*I_ERI_F2xy_Py_S_Py_C1000000003_b;
  abcd[722] = 2.0E0*I_ERI_F2xz_Py_S_Py_C1000000003_b;
  abcd[723] = 2.0E0*I_ERI_Fx2y_Py_S_Py_C1000000003_b;
  abcd[724] = 2.0E0*I_ERI_Fxyz_Py_S_Py_C1000000003_b;
  abcd[725] = 2.0E0*I_ERI_Fx2z_Py_S_Py_C1000000003_b;
  abcd[726] = 2.0E0*I_ERI_F3y_Py_S_Py_C1000000003_b;
  abcd[727] = 2.0E0*I_ERI_F2yz_Py_S_Py_C1000000003_b;
  abcd[728] = 2.0E0*I_ERI_Fy2z_Py_S_Py_C1000000003_b;
  abcd[729] = 2.0E0*I_ERI_F3z_Py_S_Py_C1000000003_b;
  abcd[760] = 2.0E0*I_ERI_F3x_Py_S_Pz_C1000000003_b;
  abcd[761] = 2.0E0*I_ERI_F2xy_Py_S_Pz_C1000000003_b;
  abcd[762] = 2.0E0*I_ERI_F2xz_Py_S_Pz_C1000000003_b;
  abcd[763] = 2.0E0*I_ERI_Fx2y_Py_S_Pz_C1000000003_b;
  abcd[764] = 2.0E0*I_ERI_Fxyz_Py_S_Pz_C1000000003_b;
  abcd[765] = 2.0E0*I_ERI_Fx2z_Py_S_Pz_C1000000003_b;
  abcd[766] = 2.0E0*I_ERI_F3y_Py_S_Pz_C1000000003_b;
  abcd[767] = 2.0E0*I_ERI_F2yz_Py_S_Pz_C1000000003_b;
  abcd[768] = 2.0E0*I_ERI_Fy2z_Py_S_Pz_C1000000003_b;
  abcd[769] = 2.0E0*I_ERI_F3z_Py_S_Pz_C1000000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[690] = 2.0E0*I_ERI_F3x_Py_Px_Px_C1001000003_b;
  abcd[691] = 2.0E0*I_ERI_F2xy_Py_Px_Px_C1001000003_b;
  abcd[692] = 2.0E0*I_ERI_F2xz_Py_Px_Px_C1001000003_b;
  abcd[693] = 2.0E0*I_ERI_Fx2y_Py_Px_Px_C1001000003_b;
  abcd[694] = 2.0E0*I_ERI_Fxyz_Py_Px_Px_C1001000003_b;
  abcd[695] = 2.0E0*I_ERI_Fx2z_Py_Px_Px_C1001000003_b;
  abcd[696] = 2.0E0*I_ERI_F3y_Py_Px_Px_C1001000003_b;
  abcd[697] = 2.0E0*I_ERI_F2yz_Py_Px_Px_C1001000003_b;
  abcd[698] = 2.0E0*I_ERI_Fy2z_Py_Px_Px_C1001000003_b;
  abcd[699] = 2.0E0*I_ERI_F3z_Py_Px_Px_C1001000003_b;
  abcd[700] = 2.0E0*I_ERI_F3x_Py_Py_Px_C1001000003_b;
  abcd[701] = 2.0E0*I_ERI_F2xy_Py_Py_Px_C1001000003_b;
  abcd[702] = 2.0E0*I_ERI_F2xz_Py_Py_Px_C1001000003_b;
  abcd[703] = 2.0E0*I_ERI_Fx2y_Py_Py_Px_C1001000003_b;
  abcd[704] = 2.0E0*I_ERI_Fxyz_Py_Py_Px_C1001000003_b;
  abcd[705] = 2.0E0*I_ERI_Fx2z_Py_Py_Px_C1001000003_b;
  abcd[706] = 2.0E0*I_ERI_F3y_Py_Py_Px_C1001000003_b;
  abcd[707] = 2.0E0*I_ERI_F2yz_Py_Py_Px_C1001000003_b;
  abcd[708] = 2.0E0*I_ERI_Fy2z_Py_Py_Px_C1001000003_b;
  abcd[709] = 2.0E0*I_ERI_F3z_Py_Py_Px_C1001000003_b;
  abcd[710] = 2.0E0*I_ERI_F3x_Py_Pz_Px_C1001000003_b;
  abcd[711] = 2.0E0*I_ERI_F2xy_Py_Pz_Px_C1001000003_b;
  abcd[712] = 2.0E0*I_ERI_F2xz_Py_Pz_Px_C1001000003_b;
  abcd[713] = 2.0E0*I_ERI_Fx2y_Py_Pz_Px_C1001000003_b;
  abcd[714] = 2.0E0*I_ERI_Fxyz_Py_Pz_Px_C1001000003_b;
  abcd[715] = 2.0E0*I_ERI_Fx2z_Py_Pz_Px_C1001000003_b;
  abcd[716] = 2.0E0*I_ERI_F3y_Py_Pz_Px_C1001000003_b;
  abcd[717] = 2.0E0*I_ERI_F2yz_Py_Pz_Px_C1001000003_b;
  abcd[718] = 2.0E0*I_ERI_Fy2z_Py_Pz_Px_C1001000003_b;
  abcd[719] = 2.0E0*I_ERI_F3z_Py_Pz_Px_C1001000003_b;
  abcd[730] = 2.0E0*I_ERI_F3x_Py_Px_Py_C1001000003_b;
  abcd[731] = 2.0E0*I_ERI_F2xy_Py_Px_Py_C1001000003_b;
  abcd[732] = 2.0E0*I_ERI_F2xz_Py_Px_Py_C1001000003_b;
  abcd[733] = 2.0E0*I_ERI_Fx2y_Py_Px_Py_C1001000003_b;
  abcd[734] = 2.0E0*I_ERI_Fxyz_Py_Px_Py_C1001000003_b;
  abcd[735] = 2.0E0*I_ERI_Fx2z_Py_Px_Py_C1001000003_b;
  abcd[736] = 2.0E0*I_ERI_F3y_Py_Px_Py_C1001000003_b;
  abcd[737] = 2.0E0*I_ERI_F2yz_Py_Px_Py_C1001000003_b;
  abcd[738] = 2.0E0*I_ERI_Fy2z_Py_Px_Py_C1001000003_b;
  abcd[739] = 2.0E0*I_ERI_F3z_Py_Px_Py_C1001000003_b;
  abcd[740] = 2.0E0*I_ERI_F3x_Py_Py_Py_C1001000003_b;
  abcd[741] = 2.0E0*I_ERI_F2xy_Py_Py_Py_C1001000003_b;
  abcd[742] = 2.0E0*I_ERI_F2xz_Py_Py_Py_C1001000003_b;
  abcd[743] = 2.0E0*I_ERI_Fx2y_Py_Py_Py_C1001000003_b;
  abcd[744] = 2.0E0*I_ERI_Fxyz_Py_Py_Py_C1001000003_b;
  abcd[745] = 2.0E0*I_ERI_Fx2z_Py_Py_Py_C1001000003_b;
  abcd[746] = 2.0E0*I_ERI_F3y_Py_Py_Py_C1001000003_b;
  abcd[747] = 2.0E0*I_ERI_F2yz_Py_Py_Py_C1001000003_b;
  abcd[748] = 2.0E0*I_ERI_Fy2z_Py_Py_Py_C1001000003_b;
  abcd[749] = 2.0E0*I_ERI_F3z_Py_Py_Py_C1001000003_b;
  abcd[750] = 2.0E0*I_ERI_F3x_Py_Pz_Py_C1001000003_b;
  abcd[751] = 2.0E0*I_ERI_F2xy_Py_Pz_Py_C1001000003_b;
  abcd[752] = 2.0E0*I_ERI_F2xz_Py_Pz_Py_C1001000003_b;
  abcd[753] = 2.0E0*I_ERI_Fx2y_Py_Pz_Py_C1001000003_b;
  abcd[754] = 2.0E0*I_ERI_Fxyz_Py_Pz_Py_C1001000003_b;
  abcd[755] = 2.0E0*I_ERI_Fx2z_Py_Pz_Py_C1001000003_b;
  abcd[756] = 2.0E0*I_ERI_F3y_Py_Pz_Py_C1001000003_b;
  abcd[757] = 2.0E0*I_ERI_F2yz_Py_Pz_Py_C1001000003_b;
  abcd[758] = 2.0E0*I_ERI_Fy2z_Py_Pz_Py_C1001000003_b;
  abcd[759] = 2.0E0*I_ERI_F3z_Py_Pz_Py_C1001000003_b;
  abcd[770] = 2.0E0*I_ERI_F3x_Py_Px_Pz_C1001000003_b;
  abcd[771] = 2.0E0*I_ERI_F2xy_Py_Px_Pz_C1001000003_b;
  abcd[772] = 2.0E0*I_ERI_F2xz_Py_Px_Pz_C1001000003_b;
  abcd[773] = 2.0E0*I_ERI_Fx2y_Py_Px_Pz_C1001000003_b;
  abcd[774] = 2.0E0*I_ERI_Fxyz_Py_Px_Pz_C1001000003_b;
  abcd[775] = 2.0E0*I_ERI_Fx2z_Py_Px_Pz_C1001000003_b;
  abcd[776] = 2.0E0*I_ERI_F3y_Py_Px_Pz_C1001000003_b;
  abcd[777] = 2.0E0*I_ERI_F2yz_Py_Px_Pz_C1001000003_b;
  abcd[778] = 2.0E0*I_ERI_Fy2z_Py_Px_Pz_C1001000003_b;
  abcd[779] = 2.0E0*I_ERI_F3z_Py_Px_Pz_C1001000003_b;
  abcd[780] = 2.0E0*I_ERI_F3x_Py_Py_Pz_C1001000003_b;
  abcd[781] = 2.0E0*I_ERI_F2xy_Py_Py_Pz_C1001000003_b;
  abcd[782] = 2.0E0*I_ERI_F2xz_Py_Py_Pz_C1001000003_b;
  abcd[783] = 2.0E0*I_ERI_Fx2y_Py_Py_Pz_C1001000003_b;
  abcd[784] = 2.0E0*I_ERI_Fxyz_Py_Py_Pz_C1001000003_b;
  abcd[785] = 2.0E0*I_ERI_Fx2z_Py_Py_Pz_C1001000003_b;
  abcd[786] = 2.0E0*I_ERI_F3y_Py_Py_Pz_C1001000003_b;
  abcd[787] = 2.0E0*I_ERI_F2yz_Py_Py_Pz_C1001000003_b;
  abcd[788] = 2.0E0*I_ERI_Fy2z_Py_Py_Pz_C1001000003_b;
  abcd[789] = 2.0E0*I_ERI_F3z_Py_Py_Pz_C1001000003_b;
  abcd[790] = 2.0E0*I_ERI_F3x_Py_Pz_Pz_C1001000003_b;
  abcd[791] = 2.0E0*I_ERI_F2xy_Py_Pz_Pz_C1001000003_b;
  abcd[792] = 2.0E0*I_ERI_F2xz_Py_Pz_Pz_C1001000003_b;
  abcd[793] = 2.0E0*I_ERI_Fx2y_Py_Pz_Pz_C1001000003_b;
  abcd[794] = 2.0E0*I_ERI_Fxyz_Py_Pz_Pz_C1001000003_b;
  abcd[795] = 2.0E0*I_ERI_Fx2z_Py_Pz_Pz_C1001000003_b;
  abcd[796] = 2.0E0*I_ERI_F3y_Py_Pz_Pz_C1001000003_b;
  abcd[797] = 2.0E0*I_ERI_F2yz_Py_Pz_Pz_C1001000003_b;
  abcd[798] = 2.0E0*I_ERI_Fy2z_Py_Pz_Pz_C1001000003_b;
  abcd[799] = 2.0E0*I_ERI_F3z_Py_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C3_b
   ************************************************************/
  abcd[800] = 2.0E0*I_ERI_F3x_Pz_S_S_C3_b;
  abcd[801] = 2.0E0*I_ERI_F2xy_Pz_S_S_C3_b;
  abcd[802] = 2.0E0*I_ERI_F2xz_Pz_S_S_C3_b;
  abcd[803] = 2.0E0*I_ERI_Fx2y_Pz_S_S_C3_b;
  abcd[804] = 2.0E0*I_ERI_Fxyz_Pz_S_S_C3_b;
  abcd[805] = 2.0E0*I_ERI_Fx2z_Pz_S_S_C3_b;
  abcd[806] = 2.0E0*I_ERI_F3y_Pz_S_S_C3_b;
  abcd[807] = 2.0E0*I_ERI_F2yz_Pz_S_S_C3_b;
  abcd[808] = 2.0E0*I_ERI_Fy2z_Pz_S_S_C3_b;
  abcd[809] = 2.0E0*I_ERI_F3z_Pz_S_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1000003_b
   ************************************************************/
  abcd[810] = 2.0E0*I_ERI_F3x_Pz_Px_S_C1000003_b;
  abcd[811] = 2.0E0*I_ERI_F2xy_Pz_Px_S_C1000003_b;
  abcd[812] = 2.0E0*I_ERI_F2xz_Pz_Px_S_C1000003_b;
  abcd[813] = 2.0E0*I_ERI_Fx2y_Pz_Px_S_C1000003_b;
  abcd[814] = 2.0E0*I_ERI_Fxyz_Pz_Px_S_C1000003_b;
  abcd[815] = 2.0E0*I_ERI_Fx2z_Pz_Px_S_C1000003_b;
  abcd[816] = 2.0E0*I_ERI_F3y_Pz_Px_S_C1000003_b;
  abcd[817] = 2.0E0*I_ERI_F2yz_Pz_Px_S_C1000003_b;
  abcd[818] = 2.0E0*I_ERI_Fy2z_Pz_Px_S_C1000003_b;
  abcd[819] = 2.0E0*I_ERI_F3z_Pz_Px_S_C1000003_b;
  abcd[820] = 2.0E0*I_ERI_F3x_Pz_Py_S_C1000003_b;
  abcd[821] = 2.0E0*I_ERI_F2xy_Pz_Py_S_C1000003_b;
  abcd[822] = 2.0E0*I_ERI_F2xz_Pz_Py_S_C1000003_b;
  abcd[823] = 2.0E0*I_ERI_Fx2y_Pz_Py_S_C1000003_b;
  abcd[824] = 2.0E0*I_ERI_Fxyz_Pz_Py_S_C1000003_b;
  abcd[825] = 2.0E0*I_ERI_Fx2z_Pz_Py_S_C1000003_b;
  abcd[826] = 2.0E0*I_ERI_F3y_Pz_Py_S_C1000003_b;
  abcd[827] = 2.0E0*I_ERI_F2yz_Pz_Py_S_C1000003_b;
  abcd[828] = 2.0E0*I_ERI_Fy2z_Pz_Py_S_C1000003_b;
  abcd[829] = 2.0E0*I_ERI_F3z_Pz_Py_S_C1000003_b;
  abcd[830] = 2.0E0*I_ERI_F3x_Pz_Pz_S_C1000003_b;
  abcd[831] = 2.0E0*I_ERI_F2xy_Pz_Pz_S_C1000003_b;
  abcd[832] = 2.0E0*I_ERI_F2xz_Pz_Pz_S_C1000003_b;
  abcd[833] = 2.0E0*I_ERI_Fx2y_Pz_Pz_S_C1000003_b;
  abcd[834] = 2.0E0*I_ERI_Fxyz_Pz_Pz_S_C1000003_b;
  abcd[835] = 2.0E0*I_ERI_Fx2z_Pz_Pz_S_C1000003_b;
  abcd[836] = 2.0E0*I_ERI_F3y_Pz_Pz_S_C1000003_b;
  abcd[837] = 2.0E0*I_ERI_F2yz_Pz_Pz_S_C1000003_b;
  abcd[838] = 2.0E0*I_ERI_Fy2z_Pz_Pz_S_C1000003_b;
  abcd[839] = 2.0E0*I_ERI_F3z_Pz_Pz_S_C1000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_P_C1000000003_b
   ************************************************************/
  abcd[840] = 2.0E0*I_ERI_F3x_Pz_S_Px_C1000000003_b;
  abcd[841] = 2.0E0*I_ERI_F2xy_Pz_S_Px_C1000000003_b;
  abcd[842] = 2.0E0*I_ERI_F2xz_Pz_S_Px_C1000000003_b;
  abcd[843] = 2.0E0*I_ERI_Fx2y_Pz_S_Px_C1000000003_b;
  abcd[844] = 2.0E0*I_ERI_Fxyz_Pz_S_Px_C1000000003_b;
  abcd[845] = 2.0E0*I_ERI_Fx2z_Pz_S_Px_C1000000003_b;
  abcd[846] = 2.0E0*I_ERI_F3y_Pz_S_Px_C1000000003_b;
  abcd[847] = 2.0E0*I_ERI_F2yz_Pz_S_Px_C1000000003_b;
  abcd[848] = 2.0E0*I_ERI_Fy2z_Pz_S_Px_C1000000003_b;
  abcd[849] = 2.0E0*I_ERI_F3z_Pz_S_Px_C1000000003_b;
  abcd[880] = 2.0E0*I_ERI_F3x_Pz_S_Py_C1000000003_b;
  abcd[881] = 2.0E0*I_ERI_F2xy_Pz_S_Py_C1000000003_b;
  abcd[882] = 2.0E0*I_ERI_F2xz_Pz_S_Py_C1000000003_b;
  abcd[883] = 2.0E0*I_ERI_Fx2y_Pz_S_Py_C1000000003_b;
  abcd[884] = 2.0E0*I_ERI_Fxyz_Pz_S_Py_C1000000003_b;
  abcd[885] = 2.0E0*I_ERI_Fx2z_Pz_S_Py_C1000000003_b;
  abcd[886] = 2.0E0*I_ERI_F3y_Pz_S_Py_C1000000003_b;
  abcd[887] = 2.0E0*I_ERI_F2yz_Pz_S_Py_C1000000003_b;
  abcd[888] = 2.0E0*I_ERI_Fy2z_Pz_S_Py_C1000000003_b;
  abcd[889] = 2.0E0*I_ERI_F3z_Pz_S_Py_C1000000003_b;
  abcd[920] = 2.0E0*I_ERI_F3x_Pz_S_Pz_C1000000003_b;
  abcd[921] = 2.0E0*I_ERI_F2xy_Pz_S_Pz_C1000000003_b;
  abcd[922] = 2.0E0*I_ERI_F2xz_Pz_S_Pz_C1000000003_b;
  abcd[923] = 2.0E0*I_ERI_Fx2y_Pz_S_Pz_C1000000003_b;
  abcd[924] = 2.0E0*I_ERI_Fxyz_Pz_S_Pz_C1000000003_b;
  abcd[925] = 2.0E0*I_ERI_Fx2z_Pz_S_Pz_C1000000003_b;
  abcd[926] = 2.0E0*I_ERI_F3y_Pz_S_Pz_C1000000003_b;
  abcd[927] = 2.0E0*I_ERI_F2yz_Pz_S_Pz_C1000000003_b;
  abcd[928] = 2.0E0*I_ERI_Fy2z_Pz_S_Pz_C1000000003_b;
  abcd[929] = 2.0E0*I_ERI_F3z_Pz_S_Pz_C1000000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_P_C1001000003_b
   ************************************************************/
  abcd[850] = 2.0E0*I_ERI_F3x_Pz_Px_Px_C1001000003_b;
  abcd[851] = 2.0E0*I_ERI_F2xy_Pz_Px_Px_C1001000003_b;
  abcd[852] = 2.0E0*I_ERI_F2xz_Pz_Px_Px_C1001000003_b;
  abcd[853] = 2.0E0*I_ERI_Fx2y_Pz_Px_Px_C1001000003_b;
  abcd[854] = 2.0E0*I_ERI_Fxyz_Pz_Px_Px_C1001000003_b;
  abcd[855] = 2.0E0*I_ERI_Fx2z_Pz_Px_Px_C1001000003_b;
  abcd[856] = 2.0E0*I_ERI_F3y_Pz_Px_Px_C1001000003_b;
  abcd[857] = 2.0E0*I_ERI_F2yz_Pz_Px_Px_C1001000003_b;
  abcd[858] = 2.0E0*I_ERI_Fy2z_Pz_Px_Px_C1001000003_b;
  abcd[859] = 2.0E0*I_ERI_F3z_Pz_Px_Px_C1001000003_b;
  abcd[860] = 2.0E0*I_ERI_F3x_Pz_Py_Px_C1001000003_b;
  abcd[861] = 2.0E0*I_ERI_F2xy_Pz_Py_Px_C1001000003_b;
  abcd[862] = 2.0E0*I_ERI_F2xz_Pz_Py_Px_C1001000003_b;
  abcd[863] = 2.0E0*I_ERI_Fx2y_Pz_Py_Px_C1001000003_b;
  abcd[864] = 2.0E0*I_ERI_Fxyz_Pz_Py_Px_C1001000003_b;
  abcd[865] = 2.0E0*I_ERI_Fx2z_Pz_Py_Px_C1001000003_b;
  abcd[866] = 2.0E0*I_ERI_F3y_Pz_Py_Px_C1001000003_b;
  abcd[867] = 2.0E0*I_ERI_F2yz_Pz_Py_Px_C1001000003_b;
  abcd[868] = 2.0E0*I_ERI_Fy2z_Pz_Py_Px_C1001000003_b;
  abcd[869] = 2.0E0*I_ERI_F3z_Pz_Py_Px_C1001000003_b;
  abcd[870] = 2.0E0*I_ERI_F3x_Pz_Pz_Px_C1001000003_b;
  abcd[871] = 2.0E0*I_ERI_F2xy_Pz_Pz_Px_C1001000003_b;
  abcd[872] = 2.0E0*I_ERI_F2xz_Pz_Pz_Px_C1001000003_b;
  abcd[873] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Px_C1001000003_b;
  abcd[874] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Px_C1001000003_b;
  abcd[875] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Px_C1001000003_b;
  abcd[876] = 2.0E0*I_ERI_F3y_Pz_Pz_Px_C1001000003_b;
  abcd[877] = 2.0E0*I_ERI_F2yz_Pz_Pz_Px_C1001000003_b;
  abcd[878] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Px_C1001000003_b;
  abcd[879] = 2.0E0*I_ERI_F3z_Pz_Pz_Px_C1001000003_b;
  abcd[890] = 2.0E0*I_ERI_F3x_Pz_Px_Py_C1001000003_b;
  abcd[891] = 2.0E0*I_ERI_F2xy_Pz_Px_Py_C1001000003_b;
  abcd[892] = 2.0E0*I_ERI_F2xz_Pz_Px_Py_C1001000003_b;
  abcd[893] = 2.0E0*I_ERI_Fx2y_Pz_Px_Py_C1001000003_b;
  abcd[894] = 2.0E0*I_ERI_Fxyz_Pz_Px_Py_C1001000003_b;
  abcd[895] = 2.0E0*I_ERI_Fx2z_Pz_Px_Py_C1001000003_b;
  abcd[896] = 2.0E0*I_ERI_F3y_Pz_Px_Py_C1001000003_b;
  abcd[897] = 2.0E0*I_ERI_F2yz_Pz_Px_Py_C1001000003_b;
  abcd[898] = 2.0E0*I_ERI_Fy2z_Pz_Px_Py_C1001000003_b;
  abcd[899] = 2.0E0*I_ERI_F3z_Pz_Px_Py_C1001000003_b;
  abcd[900] = 2.0E0*I_ERI_F3x_Pz_Py_Py_C1001000003_b;
  abcd[901] = 2.0E0*I_ERI_F2xy_Pz_Py_Py_C1001000003_b;
  abcd[902] = 2.0E0*I_ERI_F2xz_Pz_Py_Py_C1001000003_b;
  abcd[903] = 2.0E0*I_ERI_Fx2y_Pz_Py_Py_C1001000003_b;
  abcd[904] = 2.0E0*I_ERI_Fxyz_Pz_Py_Py_C1001000003_b;
  abcd[905] = 2.0E0*I_ERI_Fx2z_Pz_Py_Py_C1001000003_b;
  abcd[906] = 2.0E0*I_ERI_F3y_Pz_Py_Py_C1001000003_b;
  abcd[907] = 2.0E0*I_ERI_F2yz_Pz_Py_Py_C1001000003_b;
  abcd[908] = 2.0E0*I_ERI_Fy2z_Pz_Py_Py_C1001000003_b;
  abcd[909] = 2.0E0*I_ERI_F3z_Pz_Py_Py_C1001000003_b;
  abcd[910] = 2.0E0*I_ERI_F3x_Pz_Pz_Py_C1001000003_b;
  abcd[911] = 2.0E0*I_ERI_F2xy_Pz_Pz_Py_C1001000003_b;
  abcd[912] = 2.0E0*I_ERI_F2xz_Pz_Pz_Py_C1001000003_b;
  abcd[913] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Py_C1001000003_b;
  abcd[914] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Py_C1001000003_b;
  abcd[915] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Py_C1001000003_b;
  abcd[916] = 2.0E0*I_ERI_F3y_Pz_Pz_Py_C1001000003_b;
  abcd[917] = 2.0E0*I_ERI_F2yz_Pz_Pz_Py_C1001000003_b;
  abcd[918] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Py_C1001000003_b;
  abcd[919] = 2.0E0*I_ERI_F3z_Pz_Pz_Py_C1001000003_b;
  abcd[930] = 2.0E0*I_ERI_F3x_Pz_Px_Pz_C1001000003_b;
  abcd[931] = 2.0E0*I_ERI_F2xy_Pz_Px_Pz_C1001000003_b;
  abcd[932] = 2.0E0*I_ERI_F2xz_Pz_Px_Pz_C1001000003_b;
  abcd[933] = 2.0E0*I_ERI_Fx2y_Pz_Px_Pz_C1001000003_b;
  abcd[934] = 2.0E0*I_ERI_Fxyz_Pz_Px_Pz_C1001000003_b;
  abcd[935] = 2.0E0*I_ERI_Fx2z_Pz_Px_Pz_C1001000003_b;
  abcd[936] = 2.0E0*I_ERI_F3y_Pz_Px_Pz_C1001000003_b;
  abcd[937] = 2.0E0*I_ERI_F2yz_Pz_Px_Pz_C1001000003_b;
  abcd[938] = 2.0E0*I_ERI_Fy2z_Pz_Px_Pz_C1001000003_b;
  abcd[939] = 2.0E0*I_ERI_F3z_Pz_Px_Pz_C1001000003_b;
  abcd[940] = 2.0E0*I_ERI_F3x_Pz_Py_Pz_C1001000003_b;
  abcd[941] = 2.0E0*I_ERI_F2xy_Pz_Py_Pz_C1001000003_b;
  abcd[942] = 2.0E0*I_ERI_F2xz_Pz_Py_Pz_C1001000003_b;
  abcd[943] = 2.0E0*I_ERI_Fx2y_Pz_Py_Pz_C1001000003_b;
  abcd[944] = 2.0E0*I_ERI_Fxyz_Pz_Py_Pz_C1001000003_b;
  abcd[945] = 2.0E0*I_ERI_Fx2z_Pz_Py_Pz_C1001000003_b;
  abcd[946] = 2.0E0*I_ERI_F3y_Pz_Py_Pz_C1001000003_b;
  abcd[947] = 2.0E0*I_ERI_F2yz_Pz_Py_Pz_C1001000003_b;
  abcd[948] = 2.0E0*I_ERI_Fy2z_Pz_Py_Pz_C1001000003_b;
  abcd[949] = 2.0E0*I_ERI_F3z_Pz_Py_Pz_C1001000003_b;
  abcd[950] = 2.0E0*I_ERI_F3x_Pz_Pz_Pz_C1001000003_b;
  abcd[951] = 2.0E0*I_ERI_F2xy_Pz_Pz_Pz_C1001000003_b;
  abcd[952] = 2.0E0*I_ERI_F2xz_Pz_Pz_Pz_C1001000003_b;
  abcd[953] = 2.0E0*I_ERI_Fx2y_Pz_Pz_Pz_C1001000003_b;
  abcd[954] = 2.0E0*I_ERI_Fxyz_Pz_Pz_Pz_C1001000003_b;
  abcd[955] = 2.0E0*I_ERI_Fx2z_Pz_Pz_Pz_C1001000003_b;
  abcd[956] = 2.0E0*I_ERI_F3y_Pz_Pz_Pz_C1001000003_b;
  abcd[957] = 2.0E0*I_ERI_F2yz_Pz_Pz_Pz_C1001000003_b;
  abcd[958] = 2.0E0*I_ERI_Fy2z_Pz_Pz_Pz_C1001000003_b;
  abcd[959] = 2.0E0*I_ERI_F3z_Pz_Pz_Pz_C1001000003_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C3_c
   ************************************************************/
  abcd[960] = 2.0E0*I_ERI_F3x_S_Px_S_C3_c;
  abcd[961] = 2.0E0*I_ERI_F2xy_S_Px_S_C3_c;
  abcd[962] = 2.0E0*I_ERI_F2xz_S_Px_S_C3_c;
  abcd[963] = 2.0E0*I_ERI_Fx2y_S_Px_S_C3_c;
  abcd[964] = 2.0E0*I_ERI_Fxyz_S_Px_S_C3_c;
  abcd[965] = 2.0E0*I_ERI_Fx2z_S_Px_S_C3_c;
  abcd[966] = 2.0E0*I_ERI_F3y_S_Px_S_C3_c;
  abcd[967] = 2.0E0*I_ERI_F2yz_S_Px_S_C3_c;
  abcd[968] = 2.0E0*I_ERI_Fy2z_S_Px_S_C3_c;
  abcd[969] = 2.0E0*I_ERI_F3z_S_Px_S_C3_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[970] = 2.0E0*I_ERI_F3x_S_D2x_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[971] = 2.0E0*I_ERI_F2xy_S_D2x_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[972] = 2.0E0*I_ERI_F2xz_S_D2x_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[973] = 2.0E0*I_ERI_Fx2y_S_D2x_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[974] = 2.0E0*I_ERI_Fxyz_S_D2x_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[975] = 2.0E0*I_ERI_Fx2z_S_D2x_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[976] = 2.0E0*I_ERI_F3y_S_D2x_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[977] = 2.0E0*I_ERI_F2yz_S_D2x_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[978] = 2.0E0*I_ERI_Fy2z_S_D2x_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[979] = 2.0E0*I_ERI_F3z_S_D2x_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;
  abcd[980] = 2.0E0*I_ERI_F3x_S_Dxy_S_C1000003_c;
  abcd[981] = 2.0E0*I_ERI_F2xy_S_Dxy_S_C1000003_c;
  abcd[982] = 2.0E0*I_ERI_F2xz_S_Dxy_S_C1000003_c;
  abcd[983] = 2.0E0*I_ERI_Fx2y_S_Dxy_S_C1000003_c;
  abcd[984] = 2.0E0*I_ERI_Fxyz_S_Dxy_S_C1000003_c;
  abcd[985] = 2.0E0*I_ERI_Fx2z_S_Dxy_S_C1000003_c;
  abcd[986] = 2.0E0*I_ERI_F3y_S_Dxy_S_C1000003_c;
  abcd[987] = 2.0E0*I_ERI_F2yz_S_Dxy_S_C1000003_c;
  abcd[988] = 2.0E0*I_ERI_Fy2z_S_Dxy_S_C1000003_c;
  abcd[989] = 2.0E0*I_ERI_F3z_S_Dxy_S_C1000003_c;
  abcd[990] = 2.0E0*I_ERI_F3x_S_Dxz_S_C1000003_c;
  abcd[991] = 2.0E0*I_ERI_F2xy_S_Dxz_S_C1000003_c;
  abcd[992] = 2.0E0*I_ERI_F2xz_S_Dxz_S_C1000003_c;
  abcd[993] = 2.0E0*I_ERI_Fx2y_S_Dxz_S_C1000003_c;
  abcd[994] = 2.0E0*I_ERI_Fxyz_S_Dxz_S_C1000003_c;
  abcd[995] = 2.0E0*I_ERI_Fx2z_S_Dxz_S_C1000003_c;
  abcd[996] = 2.0E0*I_ERI_F3y_S_Dxz_S_C1000003_c;
  abcd[997] = 2.0E0*I_ERI_F2yz_S_Dxz_S_C1000003_c;
  abcd[998] = 2.0E0*I_ERI_Fy2z_S_Dxz_S_C1000003_c;
  abcd[999] = 2.0E0*I_ERI_F3z_S_Dxz_S_C1000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000000003_c
   ************************************************************/
  abcd[1000] = 2.0E0*I_ERI_F3x_S_Px_Px_C1000000003_c;
  abcd[1001] = 2.0E0*I_ERI_F2xy_S_Px_Px_C1000000003_c;
  abcd[1002] = 2.0E0*I_ERI_F2xz_S_Px_Px_C1000000003_c;
  abcd[1003] = 2.0E0*I_ERI_Fx2y_S_Px_Px_C1000000003_c;
  abcd[1004] = 2.0E0*I_ERI_Fxyz_S_Px_Px_C1000000003_c;
  abcd[1005] = 2.0E0*I_ERI_Fx2z_S_Px_Px_C1000000003_c;
  abcd[1006] = 2.0E0*I_ERI_F3y_S_Px_Px_C1000000003_c;
  abcd[1007] = 2.0E0*I_ERI_F2yz_S_Px_Px_C1000000003_c;
  abcd[1008] = 2.0E0*I_ERI_Fy2z_S_Px_Px_C1000000003_c;
  abcd[1009] = 2.0E0*I_ERI_F3z_S_Px_Px_C1000000003_c;
  abcd[1040] = 2.0E0*I_ERI_F3x_S_Px_Py_C1000000003_c;
  abcd[1041] = 2.0E0*I_ERI_F2xy_S_Px_Py_C1000000003_c;
  abcd[1042] = 2.0E0*I_ERI_F2xz_S_Px_Py_C1000000003_c;
  abcd[1043] = 2.0E0*I_ERI_Fx2y_S_Px_Py_C1000000003_c;
  abcd[1044] = 2.0E0*I_ERI_Fxyz_S_Px_Py_C1000000003_c;
  abcd[1045] = 2.0E0*I_ERI_Fx2z_S_Px_Py_C1000000003_c;
  abcd[1046] = 2.0E0*I_ERI_F3y_S_Px_Py_C1000000003_c;
  abcd[1047] = 2.0E0*I_ERI_F2yz_S_Px_Py_C1000000003_c;
  abcd[1048] = 2.0E0*I_ERI_Fy2z_S_Px_Py_C1000000003_c;
  abcd[1049] = 2.0E0*I_ERI_F3z_S_Px_Py_C1000000003_c;
  abcd[1080] = 2.0E0*I_ERI_F3x_S_Px_Pz_C1000000003_c;
  abcd[1081] = 2.0E0*I_ERI_F2xy_S_Px_Pz_C1000000003_c;
  abcd[1082] = 2.0E0*I_ERI_F2xz_S_Px_Pz_C1000000003_c;
  abcd[1083] = 2.0E0*I_ERI_Fx2y_S_Px_Pz_C1000000003_c;
  abcd[1084] = 2.0E0*I_ERI_Fxyz_S_Px_Pz_C1000000003_c;
  abcd[1085] = 2.0E0*I_ERI_Fx2z_S_Px_Pz_C1000000003_c;
  abcd[1086] = 2.0E0*I_ERI_F3y_S_Px_Pz_C1000000003_c;
  abcd[1087] = 2.0E0*I_ERI_F2yz_S_Px_Pz_C1000000003_c;
  abcd[1088] = 2.0E0*I_ERI_Fy2z_S_Px_Pz_C1000000003_c;
  abcd[1089] = 2.0E0*I_ERI_F3z_S_Px_Pz_C1000000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[1010] = 2.0E0*I_ERI_F3x_S_D2x_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[1011] = 2.0E0*I_ERI_F2xy_S_D2x_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[1012] = 2.0E0*I_ERI_F2xz_S_D2x_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[1013] = 2.0E0*I_ERI_Fx2y_S_D2x_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[1014] = 2.0E0*I_ERI_Fxyz_S_D2x_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[1015] = 2.0E0*I_ERI_Fx2z_S_D2x_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[1016] = 2.0E0*I_ERI_F3y_S_D2x_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[1017] = 2.0E0*I_ERI_F2yz_S_D2x_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[1018] = 2.0E0*I_ERI_Fy2z_S_D2x_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[1019] = 2.0E0*I_ERI_F3z_S_D2x_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[1020] = 2.0E0*I_ERI_F3x_S_Dxy_Px_C1001000003_c;
  abcd[1021] = 2.0E0*I_ERI_F2xy_S_Dxy_Px_C1001000003_c;
  abcd[1022] = 2.0E0*I_ERI_F2xz_S_Dxy_Px_C1001000003_c;
  abcd[1023] = 2.0E0*I_ERI_Fx2y_S_Dxy_Px_C1001000003_c;
  abcd[1024] = 2.0E0*I_ERI_Fxyz_S_Dxy_Px_C1001000003_c;
  abcd[1025] = 2.0E0*I_ERI_Fx2z_S_Dxy_Px_C1001000003_c;
  abcd[1026] = 2.0E0*I_ERI_F3y_S_Dxy_Px_C1001000003_c;
  abcd[1027] = 2.0E0*I_ERI_F2yz_S_Dxy_Px_C1001000003_c;
  abcd[1028] = 2.0E0*I_ERI_Fy2z_S_Dxy_Px_C1001000003_c;
  abcd[1029] = 2.0E0*I_ERI_F3z_S_Dxy_Px_C1001000003_c;
  abcd[1030] = 2.0E0*I_ERI_F3x_S_Dxz_Px_C1001000003_c;
  abcd[1031] = 2.0E0*I_ERI_F2xy_S_Dxz_Px_C1001000003_c;
  abcd[1032] = 2.0E0*I_ERI_F2xz_S_Dxz_Px_C1001000003_c;
  abcd[1033] = 2.0E0*I_ERI_Fx2y_S_Dxz_Px_C1001000003_c;
  abcd[1034] = 2.0E0*I_ERI_Fxyz_S_Dxz_Px_C1001000003_c;
  abcd[1035] = 2.0E0*I_ERI_Fx2z_S_Dxz_Px_C1001000003_c;
  abcd[1036] = 2.0E0*I_ERI_F3y_S_Dxz_Px_C1001000003_c;
  abcd[1037] = 2.0E0*I_ERI_F2yz_S_Dxz_Px_C1001000003_c;
  abcd[1038] = 2.0E0*I_ERI_Fy2z_S_Dxz_Px_C1001000003_c;
  abcd[1039] = 2.0E0*I_ERI_F3z_S_Dxz_Px_C1001000003_c;
  abcd[1050] = 2.0E0*I_ERI_F3x_S_D2x_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[1051] = 2.0E0*I_ERI_F2xy_S_D2x_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[1052] = 2.0E0*I_ERI_F2xz_S_D2x_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[1053] = 2.0E0*I_ERI_Fx2y_S_D2x_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[1054] = 2.0E0*I_ERI_Fxyz_S_D2x_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[1055] = 2.0E0*I_ERI_Fx2z_S_D2x_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[1056] = 2.0E0*I_ERI_F3y_S_D2x_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[1057] = 2.0E0*I_ERI_F2yz_S_D2x_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[1058] = 2.0E0*I_ERI_Fy2z_S_D2x_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[1059] = 2.0E0*I_ERI_F3z_S_D2x_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[1060] = 2.0E0*I_ERI_F3x_S_Dxy_Py_C1001000003_c;
  abcd[1061] = 2.0E0*I_ERI_F2xy_S_Dxy_Py_C1001000003_c;
  abcd[1062] = 2.0E0*I_ERI_F2xz_S_Dxy_Py_C1001000003_c;
  abcd[1063] = 2.0E0*I_ERI_Fx2y_S_Dxy_Py_C1001000003_c;
  abcd[1064] = 2.0E0*I_ERI_Fxyz_S_Dxy_Py_C1001000003_c;
  abcd[1065] = 2.0E0*I_ERI_Fx2z_S_Dxy_Py_C1001000003_c;
  abcd[1066] = 2.0E0*I_ERI_F3y_S_Dxy_Py_C1001000003_c;
  abcd[1067] = 2.0E0*I_ERI_F2yz_S_Dxy_Py_C1001000003_c;
  abcd[1068] = 2.0E0*I_ERI_Fy2z_S_Dxy_Py_C1001000003_c;
  abcd[1069] = 2.0E0*I_ERI_F3z_S_Dxy_Py_C1001000003_c;
  abcd[1070] = 2.0E0*I_ERI_F3x_S_Dxz_Py_C1001000003_c;
  abcd[1071] = 2.0E0*I_ERI_F2xy_S_Dxz_Py_C1001000003_c;
  abcd[1072] = 2.0E0*I_ERI_F2xz_S_Dxz_Py_C1001000003_c;
  abcd[1073] = 2.0E0*I_ERI_Fx2y_S_Dxz_Py_C1001000003_c;
  abcd[1074] = 2.0E0*I_ERI_Fxyz_S_Dxz_Py_C1001000003_c;
  abcd[1075] = 2.0E0*I_ERI_Fx2z_S_Dxz_Py_C1001000003_c;
  abcd[1076] = 2.0E0*I_ERI_F3y_S_Dxz_Py_C1001000003_c;
  abcd[1077] = 2.0E0*I_ERI_F2yz_S_Dxz_Py_C1001000003_c;
  abcd[1078] = 2.0E0*I_ERI_Fy2z_S_Dxz_Py_C1001000003_c;
  abcd[1079] = 2.0E0*I_ERI_F3z_S_Dxz_Py_C1001000003_c;
  abcd[1090] = 2.0E0*I_ERI_F3x_S_D2x_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[1091] = 2.0E0*I_ERI_F2xy_S_D2x_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[1092] = 2.0E0*I_ERI_F2xz_S_D2x_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[1093] = 2.0E0*I_ERI_Fx2y_S_D2x_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[1094] = 2.0E0*I_ERI_Fxyz_S_D2x_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[1095] = 2.0E0*I_ERI_Fx2z_S_D2x_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[1096] = 2.0E0*I_ERI_F3y_S_D2x_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[1097] = 2.0E0*I_ERI_F2yz_S_D2x_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[1098] = 2.0E0*I_ERI_Fy2z_S_D2x_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[1099] = 2.0E0*I_ERI_F3z_S_D2x_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;
  abcd[1100] = 2.0E0*I_ERI_F3x_S_Dxy_Pz_C1001000003_c;
  abcd[1101] = 2.0E0*I_ERI_F2xy_S_Dxy_Pz_C1001000003_c;
  abcd[1102] = 2.0E0*I_ERI_F2xz_S_Dxy_Pz_C1001000003_c;
  abcd[1103] = 2.0E0*I_ERI_Fx2y_S_Dxy_Pz_C1001000003_c;
  abcd[1104] = 2.0E0*I_ERI_Fxyz_S_Dxy_Pz_C1001000003_c;
  abcd[1105] = 2.0E0*I_ERI_Fx2z_S_Dxy_Pz_C1001000003_c;
  abcd[1106] = 2.0E0*I_ERI_F3y_S_Dxy_Pz_C1001000003_c;
  abcd[1107] = 2.0E0*I_ERI_F2yz_S_Dxy_Pz_C1001000003_c;
  abcd[1108] = 2.0E0*I_ERI_Fy2z_S_Dxy_Pz_C1001000003_c;
  abcd[1109] = 2.0E0*I_ERI_F3z_S_Dxy_Pz_C1001000003_c;
  abcd[1110] = 2.0E0*I_ERI_F3x_S_Dxz_Pz_C1001000003_c;
  abcd[1111] = 2.0E0*I_ERI_F2xy_S_Dxz_Pz_C1001000003_c;
  abcd[1112] = 2.0E0*I_ERI_F2xz_S_Dxz_Pz_C1001000003_c;
  abcd[1113] = 2.0E0*I_ERI_Fx2y_S_Dxz_Pz_C1001000003_c;
  abcd[1114] = 2.0E0*I_ERI_Fxyz_S_Dxz_Pz_C1001000003_c;
  abcd[1115] = 2.0E0*I_ERI_Fx2z_S_Dxz_Pz_C1001000003_c;
  abcd[1116] = 2.0E0*I_ERI_F3y_S_Dxz_Pz_C1001000003_c;
  abcd[1117] = 2.0E0*I_ERI_F2yz_S_Dxz_Pz_C1001000003_c;
  abcd[1118] = 2.0E0*I_ERI_Fy2z_S_Dxz_Pz_C1001000003_c;
  abcd[1119] = 2.0E0*I_ERI_F3z_S_Dxz_Pz_C1001000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C3_c
   ************************************************************/
  abcd[1120] = 2.0E0*I_ERI_F3x_S_Py_S_C3_c;
  abcd[1121] = 2.0E0*I_ERI_F2xy_S_Py_S_C3_c;
  abcd[1122] = 2.0E0*I_ERI_F2xz_S_Py_S_C3_c;
  abcd[1123] = 2.0E0*I_ERI_Fx2y_S_Py_S_C3_c;
  abcd[1124] = 2.0E0*I_ERI_Fxyz_S_Py_S_C3_c;
  abcd[1125] = 2.0E0*I_ERI_Fx2z_S_Py_S_C3_c;
  abcd[1126] = 2.0E0*I_ERI_F3y_S_Py_S_C3_c;
  abcd[1127] = 2.0E0*I_ERI_F2yz_S_Py_S_C3_c;
  abcd[1128] = 2.0E0*I_ERI_Fy2z_S_Py_S_C3_c;
  abcd[1129] = 2.0E0*I_ERI_F3z_S_Py_S_C3_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[1130] = 2.0E0*I_ERI_F3x_S_Dxy_S_C1000003_c;
  abcd[1131] = 2.0E0*I_ERI_F2xy_S_Dxy_S_C1000003_c;
  abcd[1132] = 2.0E0*I_ERI_F2xz_S_Dxy_S_C1000003_c;
  abcd[1133] = 2.0E0*I_ERI_Fx2y_S_Dxy_S_C1000003_c;
  abcd[1134] = 2.0E0*I_ERI_Fxyz_S_Dxy_S_C1000003_c;
  abcd[1135] = 2.0E0*I_ERI_Fx2z_S_Dxy_S_C1000003_c;
  abcd[1136] = 2.0E0*I_ERI_F3y_S_Dxy_S_C1000003_c;
  abcd[1137] = 2.0E0*I_ERI_F2yz_S_Dxy_S_C1000003_c;
  abcd[1138] = 2.0E0*I_ERI_Fy2z_S_Dxy_S_C1000003_c;
  abcd[1139] = 2.0E0*I_ERI_F3z_S_Dxy_S_C1000003_c;
  abcd[1140] = 2.0E0*I_ERI_F3x_S_D2y_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[1141] = 2.0E0*I_ERI_F2xy_S_D2y_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[1142] = 2.0E0*I_ERI_F2xz_S_D2y_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[1143] = 2.0E0*I_ERI_Fx2y_S_D2y_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[1144] = 2.0E0*I_ERI_Fxyz_S_D2y_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[1145] = 2.0E0*I_ERI_Fx2z_S_D2y_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[1146] = 2.0E0*I_ERI_F3y_S_D2y_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[1147] = 2.0E0*I_ERI_F2yz_S_D2y_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[1148] = 2.0E0*I_ERI_Fy2z_S_D2y_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[1149] = 2.0E0*I_ERI_F3z_S_D2y_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;
  abcd[1150] = 2.0E0*I_ERI_F3x_S_Dyz_S_C1000003_c;
  abcd[1151] = 2.0E0*I_ERI_F2xy_S_Dyz_S_C1000003_c;
  abcd[1152] = 2.0E0*I_ERI_F2xz_S_Dyz_S_C1000003_c;
  abcd[1153] = 2.0E0*I_ERI_Fx2y_S_Dyz_S_C1000003_c;
  abcd[1154] = 2.0E0*I_ERI_Fxyz_S_Dyz_S_C1000003_c;
  abcd[1155] = 2.0E0*I_ERI_Fx2z_S_Dyz_S_C1000003_c;
  abcd[1156] = 2.0E0*I_ERI_F3y_S_Dyz_S_C1000003_c;
  abcd[1157] = 2.0E0*I_ERI_F2yz_S_Dyz_S_C1000003_c;
  abcd[1158] = 2.0E0*I_ERI_Fy2z_S_Dyz_S_C1000003_c;
  abcd[1159] = 2.0E0*I_ERI_F3z_S_Dyz_S_C1000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000000003_c
   ************************************************************/
  abcd[1160] = 2.0E0*I_ERI_F3x_S_Py_Px_C1000000003_c;
  abcd[1161] = 2.0E0*I_ERI_F2xy_S_Py_Px_C1000000003_c;
  abcd[1162] = 2.0E0*I_ERI_F2xz_S_Py_Px_C1000000003_c;
  abcd[1163] = 2.0E0*I_ERI_Fx2y_S_Py_Px_C1000000003_c;
  abcd[1164] = 2.0E0*I_ERI_Fxyz_S_Py_Px_C1000000003_c;
  abcd[1165] = 2.0E0*I_ERI_Fx2z_S_Py_Px_C1000000003_c;
  abcd[1166] = 2.0E0*I_ERI_F3y_S_Py_Px_C1000000003_c;
  abcd[1167] = 2.0E0*I_ERI_F2yz_S_Py_Px_C1000000003_c;
  abcd[1168] = 2.0E0*I_ERI_Fy2z_S_Py_Px_C1000000003_c;
  abcd[1169] = 2.0E0*I_ERI_F3z_S_Py_Px_C1000000003_c;
  abcd[1200] = 2.0E0*I_ERI_F3x_S_Py_Py_C1000000003_c;
  abcd[1201] = 2.0E0*I_ERI_F2xy_S_Py_Py_C1000000003_c;
  abcd[1202] = 2.0E0*I_ERI_F2xz_S_Py_Py_C1000000003_c;
  abcd[1203] = 2.0E0*I_ERI_Fx2y_S_Py_Py_C1000000003_c;
  abcd[1204] = 2.0E0*I_ERI_Fxyz_S_Py_Py_C1000000003_c;
  abcd[1205] = 2.0E0*I_ERI_Fx2z_S_Py_Py_C1000000003_c;
  abcd[1206] = 2.0E0*I_ERI_F3y_S_Py_Py_C1000000003_c;
  abcd[1207] = 2.0E0*I_ERI_F2yz_S_Py_Py_C1000000003_c;
  abcd[1208] = 2.0E0*I_ERI_Fy2z_S_Py_Py_C1000000003_c;
  abcd[1209] = 2.0E0*I_ERI_F3z_S_Py_Py_C1000000003_c;
  abcd[1240] = 2.0E0*I_ERI_F3x_S_Py_Pz_C1000000003_c;
  abcd[1241] = 2.0E0*I_ERI_F2xy_S_Py_Pz_C1000000003_c;
  abcd[1242] = 2.0E0*I_ERI_F2xz_S_Py_Pz_C1000000003_c;
  abcd[1243] = 2.0E0*I_ERI_Fx2y_S_Py_Pz_C1000000003_c;
  abcd[1244] = 2.0E0*I_ERI_Fxyz_S_Py_Pz_C1000000003_c;
  abcd[1245] = 2.0E0*I_ERI_Fx2z_S_Py_Pz_C1000000003_c;
  abcd[1246] = 2.0E0*I_ERI_F3y_S_Py_Pz_C1000000003_c;
  abcd[1247] = 2.0E0*I_ERI_F2yz_S_Py_Pz_C1000000003_c;
  abcd[1248] = 2.0E0*I_ERI_Fy2z_S_Py_Pz_C1000000003_c;
  abcd[1249] = 2.0E0*I_ERI_F3z_S_Py_Pz_C1000000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[1170] = 2.0E0*I_ERI_F3x_S_Dxy_Px_C1001000003_c;
  abcd[1171] = 2.0E0*I_ERI_F2xy_S_Dxy_Px_C1001000003_c;
  abcd[1172] = 2.0E0*I_ERI_F2xz_S_Dxy_Px_C1001000003_c;
  abcd[1173] = 2.0E0*I_ERI_Fx2y_S_Dxy_Px_C1001000003_c;
  abcd[1174] = 2.0E0*I_ERI_Fxyz_S_Dxy_Px_C1001000003_c;
  abcd[1175] = 2.0E0*I_ERI_Fx2z_S_Dxy_Px_C1001000003_c;
  abcd[1176] = 2.0E0*I_ERI_F3y_S_Dxy_Px_C1001000003_c;
  abcd[1177] = 2.0E0*I_ERI_F2yz_S_Dxy_Px_C1001000003_c;
  abcd[1178] = 2.0E0*I_ERI_Fy2z_S_Dxy_Px_C1001000003_c;
  abcd[1179] = 2.0E0*I_ERI_F3z_S_Dxy_Px_C1001000003_c;
  abcd[1180] = 2.0E0*I_ERI_F3x_S_D2y_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[1181] = 2.0E0*I_ERI_F2xy_S_D2y_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[1182] = 2.0E0*I_ERI_F2xz_S_D2y_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[1183] = 2.0E0*I_ERI_Fx2y_S_D2y_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[1184] = 2.0E0*I_ERI_Fxyz_S_D2y_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[1185] = 2.0E0*I_ERI_Fx2z_S_D2y_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[1186] = 2.0E0*I_ERI_F3y_S_D2y_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[1187] = 2.0E0*I_ERI_F2yz_S_D2y_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[1188] = 2.0E0*I_ERI_Fy2z_S_D2y_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[1189] = 2.0E0*I_ERI_F3z_S_D2y_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[1190] = 2.0E0*I_ERI_F3x_S_Dyz_Px_C1001000003_c;
  abcd[1191] = 2.0E0*I_ERI_F2xy_S_Dyz_Px_C1001000003_c;
  abcd[1192] = 2.0E0*I_ERI_F2xz_S_Dyz_Px_C1001000003_c;
  abcd[1193] = 2.0E0*I_ERI_Fx2y_S_Dyz_Px_C1001000003_c;
  abcd[1194] = 2.0E0*I_ERI_Fxyz_S_Dyz_Px_C1001000003_c;
  abcd[1195] = 2.0E0*I_ERI_Fx2z_S_Dyz_Px_C1001000003_c;
  abcd[1196] = 2.0E0*I_ERI_F3y_S_Dyz_Px_C1001000003_c;
  abcd[1197] = 2.0E0*I_ERI_F2yz_S_Dyz_Px_C1001000003_c;
  abcd[1198] = 2.0E0*I_ERI_Fy2z_S_Dyz_Px_C1001000003_c;
  abcd[1199] = 2.0E0*I_ERI_F3z_S_Dyz_Px_C1001000003_c;
  abcd[1210] = 2.0E0*I_ERI_F3x_S_Dxy_Py_C1001000003_c;
  abcd[1211] = 2.0E0*I_ERI_F2xy_S_Dxy_Py_C1001000003_c;
  abcd[1212] = 2.0E0*I_ERI_F2xz_S_Dxy_Py_C1001000003_c;
  abcd[1213] = 2.0E0*I_ERI_Fx2y_S_Dxy_Py_C1001000003_c;
  abcd[1214] = 2.0E0*I_ERI_Fxyz_S_Dxy_Py_C1001000003_c;
  abcd[1215] = 2.0E0*I_ERI_Fx2z_S_Dxy_Py_C1001000003_c;
  abcd[1216] = 2.0E0*I_ERI_F3y_S_Dxy_Py_C1001000003_c;
  abcd[1217] = 2.0E0*I_ERI_F2yz_S_Dxy_Py_C1001000003_c;
  abcd[1218] = 2.0E0*I_ERI_Fy2z_S_Dxy_Py_C1001000003_c;
  abcd[1219] = 2.0E0*I_ERI_F3z_S_Dxy_Py_C1001000003_c;
  abcd[1220] = 2.0E0*I_ERI_F3x_S_D2y_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[1221] = 2.0E0*I_ERI_F2xy_S_D2y_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[1222] = 2.0E0*I_ERI_F2xz_S_D2y_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[1223] = 2.0E0*I_ERI_Fx2y_S_D2y_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[1224] = 2.0E0*I_ERI_Fxyz_S_D2y_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[1225] = 2.0E0*I_ERI_Fx2z_S_D2y_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[1226] = 2.0E0*I_ERI_F3y_S_D2y_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[1227] = 2.0E0*I_ERI_F2yz_S_D2y_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[1228] = 2.0E0*I_ERI_Fy2z_S_D2y_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[1229] = 2.0E0*I_ERI_F3z_S_D2y_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[1230] = 2.0E0*I_ERI_F3x_S_Dyz_Py_C1001000003_c;
  abcd[1231] = 2.0E0*I_ERI_F2xy_S_Dyz_Py_C1001000003_c;
  abcd[1232] = 2.0E0*I_ERI_F2xz_S_Dyz_Py_C1001000003_c;
  abcd[1233] = 2.0E0*I_ERI_Fx2y_S_Dyz_Py_C1001000003_c;
  abcd[1234] = 2.0E0*I_ERI_Fxyz_S_Dyz_Py_C1001000003_c;
  abcd[1235] = 2.0E0*I_ERI_Fx2z_S_Dyz_Py_C1001000003_c;
  abcd[1236] = 2.0E0*I_ERI_F3y_S_Dyz_Py_C1001000003_c;
  abcd[1237] = 2.0E0*I_ERI_F2yz_S_Dyz_Py_C1001000003_c;
  abcd[1238] = 2.0E0*I_ERI_Fy2z_S_Dyz_Py_C1001000003_c;
  abcd[1239] = 2.0E0*I_ERI_F3z_S_Dyz_Py_C1001000003_c;
  abcd[1250] = 2.0E0*I_ERI_F3x_S_Dxy_Pz_C1001000003_c;
  abcd[1251] = 2.0E0*I_ERI_F2xy_S_Dxy_Pz_C1001000003_c;
  abcd[1252] = 2.0E0*I_ERI_F2xz_S_Dxy_Pz_C1001000003_c;
  abcd[1253] = 2.0E0*I_ERI_Fx2y_S_Dxy_Pz_C1001000003_c;
  abcd[1254] = 2.0E0*I_ERI_Fxyz_S_Dxy_Pz_C1001000003_c;
  abcd[1255] = 2.0E0*I_ERI_Fx2z_S_Dxy_Pz_C1001000003_c;
  abcd[1256] = 2.0E0*I_ERI_F3y_S_Dxy_Pz_C1001000003_c;
  abcd[1257] = 2.0E0*I_ERI_F2yz_S_Dxy_Pz_C1001000003_c;
  abcd[1258] = 2.0E0*I_ERI_Fy2z_S_Dxy_Pz_C1001000003_c;
  abcd[1259] = 2.0E0*I_ERI_F3z_S_Dxy_Pz_C1001000003_c;
  abcd[1260] = 2.0E0*I_ERI_F3x_S_D2y_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[1261] = 2.0E0*I_ERI_F2xy_S_D2y_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[1262] = 2.0E0*I_ERI_F2xz_S_D2y_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[1263] = 2.0E0*I_ERI_Fx2y_S_D2y_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[1264] = 2.0E0*I_ERI_Fxyz_S_D2y_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[1265] = 2.0E0*I_ERI_Fx2z_S_D2y_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[1266] = 2.0E0*I_ERI_F3y_S_D2y_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[1267] = 2.0E0*I_ERI_F2yz_S_D2y_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[1268] = 2.0E0*I_ERI_Fy2z_S_D2y_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[1269] = 2.0E0*I_ERI_F3z_S_D2y_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;
  abcd[1270] = 2.0E0*I_ERI_F3x_S_Dyz_Pz_C1001000003_c;
  abcd[1271] = 2.0E0*I_ERI_F2xy_S_Dyz_Pz_C1001000003_c;
  abcd[1272] = 2.0E0*I_ERI_F2xz_S_Dyz_Pz_C1001000003_c;
  abcd[1273] = 2.0E0*I_ERI_Fx2y_S_Dyz_Pz_C1001000003_c;
  abcd[1274] = 2.0E0*I_ERI_Fxyz_S_Dyz_Pz_C1001000003_c;
  abcd[1275] = 2.0E0*I_ERI_Fx2z_S_Dyz_Pz_C1001000003_c;
  abcd[1276] = 2.0E0*I_ERI_F3y_S_Dyz_Pz_C1001000003_c;
  abcd[1277] = 2.0E0*I_ERI_F2yz_S_Dyz_Pz_C1001000003_c;
  abcd[1278] = 2.0E0*I_ERI_Fy2z_S_Dyz_Pz_C1001000003_c;
  abcd[1279] = 2.0E0*I_ERI_F3z_S_Dyz_Pz_C1001000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C3_c
   ************************************************************/
  abcd[1280] = 2.0E0*I_ERI_F3x_S_Pz_S_C3_c;
  abcd[1281] = 2.0E0*I_ERI_F2xy_S_Pz_S_C3_c;
  abcd[1282] = 2.0E0*I_ERI_F2xz_S_Pz_S_C3_c;
  abcd[1283] = 2.0E0*I_ERI_Fx2y_S_Pz_S_C3_c;
  abcd[1284] = 2.0E0*I_ERI_Fxyz_S_Pz_S_C3_c;
  abcd[1285] = 2.0E0*I_ERI_Fx2z_S_Pz_S_C3_c;
  abcd[1286] = 2.0E0*I_ERI_F3y_S_Pz_S_C3_c;
  abcd[1287] = 2.0E0*I_ERI_F2yz_S_Pz_S_C3_c;
  abcd[1288] = 2.0E0*I_ERI_Fy2z_S_Pz_S_C3_c;
  abcd[1289] = 2.0E0*I_ERI_F3z_S_Pz_S_C3_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1000003
   ************************************************************/
  abcd[1290] = 2.0E0*I_ERI_F3x_S_Dxz_S_C1000003_c;
  abcd[1291] = 2.0E0*I_ERI_F2xy_S_Dxz_S_C1000003_c;
  abcd[1292] = 2.0E0*I_ERI_F2xz_S_Dxz_S_C1000003_c;
  abcd[1293] = 2.0E0*I_ERI_Fx2y_S_Dxz_S_C1000003_c;
  abcd[1294] = 2.0E0*I_ERI_Fxyz_S_Dxz_S_C1000003_c;
  abcd[1295] = 2.0E0*I_ERI_Fx2z_S_Dxz_S_C1000003_c;
  abcd[1296] = 2.0E0*I_ERI_F3y_S_Dxz_S_C1000003_c;
  abcd[1297] = 2.0E0*I_ERI_F2yz_S_Dxz_S_C1000003_c;
  abcd[1298] = 2.0E0*I_ERI_Fy2z_S_Dxz_S_C1000003_c;
  abcd[1299] = 2.0E0*I_ERI_F3z_S_Dxz_S_C1000003_c;
  abcd[1300] = 2.0E0*I_ERI_F3x_S_Dyz_S_C1000003_c;
  abcd[1301] = 2.0E0*I_ERI_F2xy_S_Dyz_S_C1000003_c;
  abcd[1302] = 2.0E0*I_ERI_F2xz_S_Dyz_S_C1000003_c;
  abcd[1303] = 2.0E0*I_ERI_Fx2y_S_Dyz_S_C1000003_c;
  abcd[1304] = 2.0E0*I_ERI_Fxyz_S_Dyz_S_C1000003_c;
  abcd[1305] = 2.0E0*I_ERI_Fx2z_S_Dyz_S_C1000003_c;
  abcd[1306] = 2.0E0*I_ERI_F3y_S_Dyz_S_C1000003_c;
  abcd[1307] = 2.0E0*I_ERI_F2yz_S_Dyz_S_C1000003_c;
  abcd[1308] = 2.0E0*I_ERI_Fy2z_S_Dyz_S_C1000003_c;
  abcd[1309] = 2.0E0*I_ERI_F3z_S_Dyz_S_C1000003_c;
  abcd[1310] = 2.0E0*I_ERI_F3x_S_D2z_S_C1000003_c-1*I_ERI_F3x_S_S_S_C1000003;
  abcd[1311] = 2.0E0*I_ERI_F2xy_S_D2z_S_C1000003_c-1*I_ERI_F2xy_S_S_S_C1000003;
  abcd[1312] = 2.0E0*I_ERI_F2xz_S_D2z_S_C1000003_c-1*I_ERI_F2xz_S_S_S_C1000003;
  abcd[1313] = 2.0E0*I_ERI_Fx2y_S_D2z_S_C1000003_c-1*I_ERI_Fx2y_S_S_S_C1000003;
  abcd[1314] = 2.0E0*I_ERI_Fxyz_S_D2z_S_C1000003_c-1*I_ERI_Fxyz_S_S_S_C1000003;
  abcd[1315] = 2.0E0*I_ERI_Fx2z_S_D2z_S_C1000003_c-1*I_ERI_Fx2z_S_S_S_C1000003;
  abcd[1316] = 2.0E0*I_ERI_F3y_S_D2z_S_C1000003_c-1*I_ERI_F3y_S_S_S_C1000003;
  abcd[1317] = 2.0E0*I_ERI_F2yz_S_D2z_S_C1000003_c-1*I_ERI_F2yz_S_S_S_C1000003;
  abcd[1318] = 2.0E0*I_ERI_Fy2z_S_D2z_S_C1000003_c-1*I_ERI_Fy2z_S_S_S_C1000003;
  abcd[1319] = 2.0E0*I_ERI_F3z_S_D2z_S_C1000003_c-1*I_ERI_F3z_S_S_S_C1000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_P_C1000000003_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_P_C1000000003_c
   ************************************************************/
  abcd[1320] = 2.0E0*I_ERI_F3x_S_Pz_Px_C1000000003_c;
  abcd[1321] = 2.0E0*I_ERI_F2xy_S_Pz_Px_C1000000003_c;
  abcd[1322] = 2.0E0*I_ERI_F2xz_S_Pz_Px_C1000000003_c;
  abcd[1323] = 2.0E0*I_ERI_Fx2y_S_Pz_Px_C1000000003_c;
  abcd[1324] = 2.0E0*I_ERI_Fxyz_S_Pz_Px_C1000000003_c;
  abcd[1325] = 2.0E0*I_ERI_Fx2z_S_Pz_Px_C1000000003_c;
  abcd[1326] = 2.0E0*I_ERI_F3y_S_Pz_Px_C1000000003_c;
  abcd[1327] = 2.0E0*I_ERI_F2yz_S_Pz_Px_C1000000003_c;
  abcd[1328] = 2.0E0*I_ERI_Fy2z_S_Pz_Px_C1000000003_c;
  abcd[1329] = 2.0E0*I_ERI_F3z_S_Pz_Px_C1000000003_c;
  abcd[1360] = 2.0E0*I_ERI_F3x_S_Pz_Py_C1000000003_c;
  abcd[1361] = 2.0E0*I_ERI_F2xy_S_Pz_Py_C1000000003_c;
  abcd[1362] = 2.0E0*I_ERI_F2xz_S_Pz_Py_C1000000003_c;
  abcd[1363] = 2.0E0*I_ERI_Fx2y_S_Pz_Py_C1000000003_c;
  abcd[1364] = 2.0E0*I_ERI_Fxyz_S_Pz_Py_C1000000003_c;
  abcd[1365] = 2.0E0*I_ERI_Fx2z_S_Pz_Py_C1000000003_c;
  abcd[1366] = 2.0E0*I_ERI_F3y_S_Pz_Py_C1000000003_c;
  abcd[1367] = 2.0E0*I_ERI_F2yz_S_Pz_Py_C1000000003_c;
  abcd[1368] = 2.0E0*I_ERI_Fy2z_S_Pz_Py_C1000000003_c;
  abcd[1369] = 2.0E0*I_ERI_F3z_S_Pz_Py_C1000000003_c;
  abcd[1400] = 2.0E0*I_ERI_F3x_S_Pz_Pz_C1000000003_c;
  abcd[1401] = 2.0E0*I_ERI_F2xy_S_Pz_Pz_C1000000003_c;
  abcd[1402] = 2.0E0*I_ERI_F2xz_S_Pz_Pz_C1000000003_c;
  abcd[1403] = 2.0E0*I_ERI_Fx2y_S_Pz_Pz_C1000000003_c;
  abcd[1404] = 2.0E0*I_ERI_Fxyz_S_Pz_Pz_C1000000003_c;
  abcd[1405] = 2.0E0*I_ERI_Fx2z_S_Pz_Pz_C1000000003_c;
  abcd[1406] = 2.0E0*I_ERI_F3y_S_Pz_Pz_C1000000003_c;
  abcd[1407] = 2.0E0*I_ERI_F2yz_S_Pz_Pz_C1000000003_c;
  abcd[1408] = 2.0E0*I_ERI_Fy2z_S_Pz_Pz_C1000000003_c;
  abcd[1409] = 2.0E0*I_ERI_F3z_S_Pz_Pz_C1000000003_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_P_C1001000003_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1001000003_c
   * RHS shell quartet name: SQ_ERI_F_S_S_P_C1001000003
   ************************************************************/
  abcd[1330] = 2.0E0*I_ERI_F3x_S_Dxz_Px_C1001000003_c;
  abcd[1331] = 2.0E0*I_ERI_F2xy_S_Dxz_Px_C1001000003_c;
  abcd[1332] = 2.0E0*I_ERI_F2xz_S_Dxz_Px_C1001000003_c;
  abcd[1333] = 2.0E0*I_ERI_Fx2y_S_Dxz_Px_C1001000003_c;
  abcd[1334] = 2.0E0*I_ERI_Fxyz_S_Dxz_Px_C1001000003_c;
  abcd[1335] = 2.0E0*I_ERI_Fx2z_S_Dxz_Px_C1001000003_c;
  abcd[1336] = 2.0E0*I_ERI_F3y_S_Dxz_Px_C1001000003_c;
  abcd[1337] = 2.0E0*I_ERI_F2yz_S_Dxz_Px_C1001000003_c;
  abcd[1338] = 2.0E0*I_ERI_Fy2z_S_Dxz_Px_C1001000003_c;
  abcd[1339] = 2.0E0*I_ERI_F3z_S_Dxz_Px_C1001000003_c;
  abcd[1340] = 2.0E0*I_ERI_F3x_S_Dyz_Px_C1001000003_c;
  abcd[1341] = 2.0E0*I_ERI_F2xy_S_Dyz_Px_C1001000003_c;
  abcd[1342] = 2.0E0*I_ERI_F2xz_S_Dyz_Px_C1001000003_c;
  abcd[1343] = 2.0E0*I_ERI_Fx2y_S_Dyz_Px_C1001000003_c;
  abcd[1344] = 2.0E0*I_ERI_Fxyz_S_Dyz_Px_C1001000003_c;
  abcd[1345] = 2.0E0*I_ERI_Fx2z_S_Dyz_Px_C1001000003_c;
  abcd[1346] = 2.0E0*I_ERI_F3y_S_Dyz_Px_C1001000003_c;
  abcd[1347] = 2.0E0*I_ERI_F2yz_S_Dyz_Px_C1001000003_c;
  abcd[1348] = 2.0E0*I_ERI_Fy2z_S_Dyz_Px_C1001000003_c;
  abcd[1349] = 2.0E0*I_ERI_F3z_S_Dyz_Px_C1001000003_c;
  abcd[1350] = 2.0E0*I_ERI_F3x_S_D2z_Px_C1001000003_c-1*I_ERI_F3x_S_S_Px_C1001000003;
  abcd[1351] = 2.0E0*I_ERI_F2xy_S_D2z_Px_C1001000003_c-1*I_ERI_F2xy_S_S_Px_C1001000003;
  abcd[1352] = 2.0E0*I_ERI_F2xz_S_D2z_Px_C1001000003_c-1*I_ERI_F2xz_S_S_Px_C1001000003;
  abcd[1353] = 2.0E0*I_ERI_Fx2y_S_D2z_Px_C1001000003_c-1*I_ERI_Fx2y_S_S_Px_C1001000003;
  abcd[1354] = 2.0E0*I_ERI_Fxyz_S_D2z_Px_C1001000003_c-1*I_ERI_Fxyz_S_S_Px_C1001000003;
  abcd[1355] = 2.0E0*I_ERI_Fx2z_S_D2z_Px_C1001000003_c-1*I_ERI_Fx2z_S_S_Px_C1001000003;
  abcd[1356] = 2.0E0*I_ERI_F3y_S_D2z_Px_C1001000003_c-1*I_ERI_F3y_S_S_Px_C1001000003;
  abcd[1357] = 2.0E0*I_ERI_F2yz_S_D2z_Px_C1001000003_c-1*I_ERI_F2yz_S_S_Px_C1001000003;
  abcd[1358] = 2.0E0*I_ERI_Fy2z_S_D2z_Px_C1001000003_c-1*I_ERI_Fy2z_S_S_Px_C1001000003;
  abcd[1359] = 2.0E0*I_ERI_F3z_S_D2z_Px_C1001000003_c-1*I_ERI_F3z_S_S_Px_C1001000003;
  abcd[1370] = 2.0E0*I_ERI_F3x_S_Dxz_Py_C1001000003_c;
  abcd[1371] = 2.0E0*I_ERI_F2xy_S_Dxz_Py_C1001000003_c;
  abcd[1372] = 2.0E0*I_ERI_F2xz_S_Dxz_Py_C1001000003_c;
  abcd[1373] = 2.0E0*I_ERI_Fx2y_S_Dxz_Py_C1001000003_c;
  abcd[1374] = 2.0E0*I_ERI_Fxyz_S_Dxz_Py_C1001000003_c;
  abcd[1375] = 2.0E0*I_ERI_Fx2z_S_Dxz_Py_C1001000003_c;
  abcd[1376] = 2.0E0*I_ERI_F3y_S_Dxz_Py_C1001000003_c;
  abcd[1377] = 2.0E0*I_ERI_F2yz_S_Dxz_Py_C1001000003_c;
  abcd[1378] = 2.0E0*I_ERI_Fy2z_S_Dxz_Py_C1001000003_c;
  abcd[1379] = 2.0E0*I_ERI_F3z_S_Dxz_Py_C1001000003_c;
  abcd[1380] = 2.0E0*I_ERI_F3x_S_Dyz_Py_C1001000003_c;
  abcd[1381] = 2.0E0*I_ERI_F2xy_S_Dyz_Py_C1001000003_c;
  abcd[1382] = 2.0E0*I_ERI_F2xz_S_Dyz_Py_C1001000003_c;
  abcd[1383] = 2.0E0*I_ERI_Fx2y_S_Dyz_Py_C1001000003_c;
  abcd[1384] = 2.0E0*I_ERI_Fxyz_S_Dyz_Py_C1001000003_c;
  abcd[1385] = 2.0E0*I_ERI_Fx2z_S_Dyz_Py_C1001000003_c;
  abcd[1386] = 2.0E0*I_ERI_F3y_S_Dyz_Py_C1001000003_c;
  abcd[1387] = 2.0E0*I_ERI_F2yz_S_Dyz_Py_C1001000003_c;
  abcd[1388] = 2.0E0*I_ERI_Fy2z_S_Dyz_Py_C1001000003_c;
  abcd[1389] = 2.0E0*I_ERI_F3z_S_Dyz_Py_C1001000003_c;
  abcd[1390] = 2.0E0*I_ERI_F3x_S_D2z_Py_C1001000003_c-1*I_ERI_F3x_S_S_Py_C1001000003;
  abcd[1391] = 2.0E0*I_ERI_F2xy_S_D2z_Py_C1001000003_c-1*I_ERI_F2xy_S_S_Py_C1001000003;
  abcd[1392] = 2.0E0*I_ERI_F2xz_S_D2z_Py_C1001000003_c-1*I_ERI_F2xz_S_S_Py_C1001000003;
  abcd[1393] = 2.0E0*I_ERI_Fx2y_S_D2z_Py_C1001000003_c-1*I_ERI_Fx2y_S_S_Py_C1001000003;
  abcd[1394] = 2.0E0*I_ERI_Fxyz_S_D2z_Py_C1001000003_c-1*I_ERI_Fxyz_S_S_Py_C1001000003;
  abcd[1395] = 2.0E0*I_ERI_Fx2z_S_D2z_Py_C1001000003_c-1*I_ERI_Fx2z_S_S_Py_C1001000003;
  abcd[1396] = 2.0E0*I_ERI_F3y_S_D2z_Py_C1001000003_c-1*I_ERI_F3y_S_S_Py_C1001000003;
  abcd[1397] = 2.0E0*I_ERI_F2yz_S_D2z_Py_C1001000003_c-1*I_ERI_F2yz_S_S_Py_C1001000003;
  abcd[1398] = 2.0E0*I_ERI_Fy2z_S_D2z_Py_C1001000003_c-1*I_ERI_Fy2z_S_S_Py_C1001000003;
  abcd[1399] = 2.0E0*I_ERI_F3z_S_D2z_Py_C1001000003_c-1*I_ERI_F3z_S_S_Py_C1001000003;
  abcd[1410] = 2.0E0*I_ERI_F3x_S_Dxz_Pz_C1001000003_c;
  abcd[1411] = 2.0E0*I_ERI_F2xy_S_Dxz_Pz_C1001000003_c;
  abcd[1412] = 2.0E0*I_ERI_F2xz_S_Dxz_Pz_C1001000003_c;
  abcd[1413] = 2.0E0*I_ERI_Fx2y_S_Dxz_Pz_C1001000003_c;
  abcd[1414] = 2.0E0*I_ERI_Fxyz_S_Dxz_Pz_C1001000003_c;
  abcd[1415] = 2.0E0*I_ERI_Fx2z_S_Dxz_Pz_C1001000003_c;
  abcd[1416] = 2.0E0*I_ERI_F3y_S_Dxz_Pz_C1001000003_c;
  abcd[1417] = 2.0E0*I_ERI_F2yz_S_Dxz_Pz_C1001000003_c;
  abcd[1418] = 2.0E0*I_ERI_Fy2z_S_Dxz_Pz_C1001000003_c;
  abcd[1419] = 2.0E0*I_ERI_F3z_S_Dxz_Pz_C1001000003_c;
  abcd[1420] = 2.0E0*I_ERI_F3x_S_Dyz_Pz_C1001000003_c;
  abcd[1421] = 2.0E0*I_ERI_F2xy_S_Dyz_Pz_C1001000003_c;
  abcd[1422] = 2.0E0*I_ERI_F2xz_S_Dyz_Pz_C1001000003_c;
  abcd[1423] = 2.0E0*I_ERI_Fx2y_S_Dyz_Pz_C1001000003_c;
  abcd[1424] = 2.0E0*I_ERI_Fxyz_S_Dyz_Pz_C1001000003_c;
  abcd[1425] = 2.0E0*I_ERI_Fx2z_S_Dyz_Pz_C1001000003_c;
  abcd[1426] = 2.0E0*I_ERI_F3y_S_Dyz_Pz_C1001000003_c;
  abcd[1427] = 2.0E0*I_ERI_F2yz_S_Dyz_Pz_C1001000003_c;
  abcd[1428] = 2.0E0*I_ERI_Fy2z_S_Dyz_Pz_C1001000003_c;
  abcd[1429] = 2.0E0*I_ERI_F3z_S_Dyz_Pz_C1001000003_c;
  abcd[1430] = 2.0E0*I_ERI_F3x_S_D2z_Pz_C1001000003_c-1*I_ERI_F3x_S_S_Pz_C1001000003;
  abcd[1431] = 2.0E0*I_ERI_F2xy_S_D2z_Pz_C1001000003_c-1*I_ERI_F2xy_S_S_Pz_C1001000003;
  abcd[1432] = 2.0E0*I_ERI_F2xz_S_D2z_Pz_C1001000003_c-1*I_ERI_F2xz_S_S_Pz_C1001000003;
  abcd[1433] = 2.0E0*I_ERI_Fx2y_S_D2z_Pz_C1001000003_c-1*I_ERI_Fx2y_S_S_Pz_C1001000003;
  abcd[1434] = 2.0E0*I_ERI_Fxyz_S_D2z_Pz_C1001000003_c-1*I_ERI_Fxyz_S_S_Pz_C1001000003;
  abcd[1435] = 2.0E0*I_ERI_Fx2z_S_D2z_Pz_C1001000003_c-1*I_ERI_Fx2z_S_S_Pz_C1001000003;
  abcd[1436] = 2.0E0*I_ERI_F3y_S_D2z_Pz_C1001000003_c-1*I_ERI_F3y_S_S_Pz_C1001000003;
  abcd[1437] = 2.0E0*I_ERI_F2yz_S_D2z_Pz_C1001000003_c-1*I_ERI_F2yz_S_S_Pz_C1001000003;
  abcd[1438] = 2.0E0*I_ERI_Fy2z_S_D2z_Pz_C1001000003_c-1*I_ERI_Fy2z_S_S_Pz_C1001000003;
  abcd[1439] = 2.0E0*I_ERI_F3z_S_D2z_Pz_C1001000003_c-1*I_ERI_F3z_S_S_Pz_C1001000003;
}
