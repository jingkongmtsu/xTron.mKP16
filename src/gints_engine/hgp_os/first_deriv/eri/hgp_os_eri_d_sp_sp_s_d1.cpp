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
// BRA1 as redundant position, total RHS integrals evaluated as: 32753
// BRA2 as redundant position, total RHS integrals evaluated as: 33360
// KET1 as redundant position, total RHS integrals evaluated as: 33776
// KET2 as redundant position, total RHS integrals evaluated as: 30059
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

void hgp_os_eri_d_sp_sp_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_S_S_C2_a = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C2_a = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C2_a = 0.0E0;
  Double I_ERI_F3y_S_S_S_C2_a = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C2_a = 0.0E0;
  Double I_ERI_F3z_S_S_S_C2_a = 0.0E0;
  Double I_ERI_Px_S_S_S_C2 = 0.0E0;
  Double I_ERI_Py_S_S_S_C2 = 0.0E0;
  Double I_ERI_Pz_S_S_S_C2 = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000002_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000002_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000002_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1000002 = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1000002 = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1000002 = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1000002 = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1000002 = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1000002 = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1000002 = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1000002 = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1000002 = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1002 = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1002 = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1002 = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1002 = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1002 = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1002 = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C2_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C2_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C2_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1000002_c = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1000002 = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1002_a = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1002_a = 0.0E0;
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
  Double I_ERI_G4x_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001002_a = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001002_a = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001002_a = 0.0E0;
  Double I_ERI_Px_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Py_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Pz_S_Px_S_C1001002 = 0.0E0;
  Double I_ERI_Px_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Py_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Pz_S_Py_S_C1001002 = 0.0E0;
  Double I_ERI_Px_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_Py_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_C1001002 = 0.0E0;
  Double I_ERI_F3x_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2x_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_C2_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C2_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_C2_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1000002_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1000002_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1000002_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1002_c = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1002_c = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1002_c = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1002_c = 0.0E0;
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
  Double I_ERI_F3x_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1001002_c = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1001002_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_F2xy_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_F2xz_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_F3y_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_F2yz_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_F3z_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_D2x_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Dxy_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Dxz_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_D2y_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_Dyz_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_D2z_S_S_S_C1001002 = 0.0E0;
  Double I_ERI_G4x_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G3xy_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G3xz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G4y_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G3yz_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_C1002_b = 0.0E0;
  Double I_ERI_G4z_S_S_S_C1002_b = 0.0E0;
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
  Double I_ERI_G4x_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G4y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G4z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_G4x_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G4y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G4z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F3x_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F3y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F3z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_F3x_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F3y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F3z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_D2x_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_D2y_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_D2z_S_Px_S_C1001002_b = 0.0E0;
  Double I_ERI_D2x_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_D2y_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_D2z_S_Py_S_C1001002_b = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_C1001002_b = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_C1001002_b = 0.0E0;

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
       * shell quartet name: SQ_ERI_F_S_S_S_C2_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C2_a_coefs = ic2*jc2*alpha;
      I_ERI_F3x_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C2_a += SQ_ERI_F_S_S_S_C2_a_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_C2
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_C2_coefs = ic2*jc2;
      I_ERI_Px_S_S_S_C2 += SQ_ERI_P_S_S_S_C2_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_C2 += SQ_ERI_P_S_S_S_C2_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_C2 += SQ_ERI_P_S_S_S_C2_coefs*I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1000002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000002_a_coefs = ic2*jc2_1*alpha;
      I_ERI_F3x_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000002_a += SQ_ERI_F_S_P_S_C1000002_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1000002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1000002_coefs = ic2*jc2_1;
      I_ERI_Px_S_Px_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1000002 += SQ_ERI_P_S_P_S_C1000002_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1002_coefs = ic2_1*jc2;
      I_ERI_D2x_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1002 += SQ_ERI_D_S_S_S_C1002_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001002_coefs = ic2_1*jc2_1;
      I_ERI_D2x_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001002 += SQ_ERI_D_S_P_S_C1001002_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C2_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C2_c_coefs = ic2*jc2*gamma;
      I_ERI_D2x_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C2_c += SQ_ERI_D_S_P_S_C2_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1000002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1000002_c_coefs = ic2*jc2_1*gamma;
      I_ERI_D2x_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1000002_c += SQ_ERI_D_S_D_S_C1000002_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1000002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1000002_coefs = ic2*jc2_1;
      I_ERI_D2x_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1000002 += SQ_ERI_D_S_S_S_C1000002_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1002_a_coefs = ic2_1*jc2*alpha;
      I_ERI_G4x_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1002_a += SQ_ERI_G_S_S_S_C1002_a_coefs*I_ERI_G4z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_C1001002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1001002_a_coefs = ic2_1*jc2_1*alpha;
      I_ERI_G4x_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1001002_a += SQ_ERI_G_S_P_S_C1001002_a_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001002_a_coefs = ic2_1*jc2_1*alpha;
      I_ERI_F3x_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001002_a += SQ_ERI_F_S_P_S_C1001002_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_C1001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_C1001002_coefs = ic2_1*jc2_1;
      I_ERI_Px_S_Px_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_C1001002 += SQ_ERI_P_S_P_S_C1001002_coefs*I_ERI_Pz_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C2_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C2_b_coefs = ic2*jc2*beta;
      I_ERI_F3x_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C2_b += SQ_ERI_F_S_S_S_C2_b_coefs*I_ERI_F3z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_P_S_C1000002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1000002_b_coefs = ic2*jc2_1*beta;
      I_ERI_F3x_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1000002_b += SQ_ERI_F_S_P_S_C1000002_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1000002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1000002_b_coefs = ic2*jc2_1*beta;
      I_ERI_D2x_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1000002_b += SQ_ERI_D_S_P_S_C1000002_b_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1002_c_coefs = ic2_1*jc2*gamma;
      I_ERI_F3x_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1002_c += SQ_ERI_F_S_P_S_C1002_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_D_S_C1001002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1001002_c_coefs = ic2_1*jc2_1*gamma;
      I_ERI_F3x_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1001002_c += SQ_ERI_F_S_D_S_C1001002_c_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1001002_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1001002_c_coefs = ic2_1*jc2_1*gamma;
      I_ERI_D2x_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1001002_c += SQ_ERI_D_S_D_S_C1001002_c_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_C1001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_C1001002_coefs = ic2_1*jc2_1;
      I_ERI_F3x_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_C1001002 += SQ_ERI_F_S_S_S_C1001002_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_C1001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_C1001002_coefs = ic2_1*jc2_1;
      I_ERI_D2x_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_C1001002 += SQ_ERI_D_S_S_S_C1001002_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_C1002_b_coefs = ic2_1*jc2*beta;
      I_ERI_G4x_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_C1002_b += SQ_ERI_G_S_S_S_C1002_b_coefs*I_ERI_G4z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_C1001002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_C1001002_b_coefs = ic2_1*jc2_1*beta;
      I_ERI_G4x_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_C1001002_b += SQ_ERI_G_S_P_S_C1001002_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_C1001002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_C1001002_b_coefs = ic2_1*jc2_1*beta;
      I_ERI_F3x_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_C1001002_b += SQ_ERI_F_S_P_S_C1001002_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_C1001002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_C1001002_b_coefs = ic2_1*jc2_1*beta;
      I_ERI_D2x_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_C1001002_b += SQ_ERI_D_S_P_S_C1001002_b_coefs*I_ERI_D2z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_P_P_S_S_C1002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C1002
   ************************************************************/
  Double I_ERI_Px_Px_S_S_C1002 = I_ERI_D2x_S_S_S_C1002+ABX*I_ERI_Px_S_S_S_C1002;
  Double I_ERI_Py_Px_S_S_C1002 = I_ERI_Dxy_S_S_S_C1002+ABX*I_ERI_Py_S_S_S_C1002;
  Double I_ERI_Pz_Px_S_S_C1002 = I_ERI_Dxz_S_S_S_C1002+ABX*I_ERI_Pz_S_S_S_C1002;
  Double I_ERI_Px_Py_S_S_C1002 = I_ERI_Dxy_S_S_S_C1002+ABY*I_ERI_Px_S_S_S_C1002;
  Double I_ERI_Py_Py_S_S_C1002 = I_ERI_D2y_S_S_S_C1002+ABY*I_ERI_Py_S_S_S_C1002;
  Double I_ERI_Pz_Py_S_S_C1002 = I_ERI_Dyz_S_S_S_C1002+ABY*I_ERI_Pz_S_S_S_C1002;
  Double I_ERI_Px_Pz_S_S_C1002 = I_ERI_Dxz_S_S_S_C1002+ABZ*I_ERI_Px_S_S_S_C1002;
  Double I_ERI_Py_Pz_S_S_C1002 = I_ERI_Dyz_S_S_S_C1002+ABZ*I_ERI_Py_S_S_S_C1002;
  Double I_ERI_Pz_Pz_S_S_C1002 = I_ERI_D2z_S_S_S_C1002+ABZ*I_ERI_Pz_S_S_S_C1002;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_C1001002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001002
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1001002
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_C1001002 = I_ERI_D2x_S_Px_S_C1001002+ABX*I_ERI_Px_S_Px_S_C1001002;
  Double I_ERI_Py_Px_Px_S_C1001002 = I_ERI_Dxy_S_Px_S_C1001002+ABX*I_ERI_Py_S_Px_S_C1001002;
  Double I_ERI_Pz_Px_Px_S_C1001002 = I_ERI_Dxz_S_Px_S_C1001002+ABX*I_ERI_Pz_S_Px_S_C1001002;
  Double I_ERI_Px_Py_Px_S_C1001002 = I_ERI_Dxy_S_Px_S_C1001002+ABY*I_ERI_Px_S_Px_S_C1001002;
  Double I_ERI_Py_Py_Px_S_C1001002 = I_ERI_D2y_S_Px_S_C1001002+ABY*I_ERI_Py_S_Px_S_C1001002;
  Double I_ERI_Pz_Py_Px_S_C1001002 = I_ERI_Dyz_S_Px_S_C1001002+ABY*I_ERI_Pz_S_Px_S_C1001002;
  Double I_ERI_Px_Pz_Px_S_C1001002 = I_ERI_Dxz_S_Px_S_C1001002+ABZ*I_ERI_Px_S_Px_S_C1001002;
  Double I_ERI_Py_Pz_Px_S_C1001002 = I_ERI_Dyz_S_Px_S_C1001002+ABZ*I_ERI_Py_S_Px_S_C1001002;
  Double I_ERI_Pz_Pz_Px_S_C1001002 = I_ERI_D2z_S_Px_S_C1001002+ABZ*I_ERI_Pz_S_Px_S_C1001002;
  Double I_ERI_Px_Px_Py_S_C1001002 = I_ERI_D2x_S_Py_S_C1001002+ABX*I_ERI_Px_S_Py_S_C1001002;
  Double I_ERI_Py_Px_Py_S_C1001002 = I_ERI_Dxy_S_Py_S_C1001002+ABX*I_ERI_Py_S_Py_S_C1001002;
  Double I_ERI_Pz_Px_Py_S_C1001002 = I_ERI_Dxz_S_Py_S_C1001002+ABX*I_ERI_Pz_S_Py_S_C1001002;
  Double I_ERI_Px_Py_Py_S_C1001002 = I_ERI_Dxy_S_Py_S_C1001002+ABY*I_ERI_Px_S_Py_S_C1001002;
  Double I_ERI_Py_Py_Py_S_C1001002 = I_ERI_D2y_S_Py_S_C1001002+ABY*I_ERI_Py_S_Py_S_C1001002;
  Double I_ERI_Pz_Py_Py_S_C1001002 = I_ERI_Dyz_S_Py_S_C1001002+ABY*I_ERI_Pz_S_Py_S_C1001002;
  Double I_ERI_Px_Pz_Py_S_C1001002 = I_ERI_Dxz_S_Py_S_C1001002+ABZ*I_ERI_Px_S_Py_S_C1001002;
  Double I_ERI_Py_Pz_Py_S_C1001002 = I_ERI_Dyz_S_Py_S_C1001002+ABZ*I_ERI_Py_S_Py_S_C1001002;
  Double I_ERI_Pz_Pz_Py_S_C1001002 = I_ERI_D2z_S_Py_S_C1001002+ABZ*I_ERI_Pz_S_Py_S_C1001002;
  Double I_ERI_Px_Px_Pz_S_C1001002 = I_ERI_D2x_S_Pz_S_C1001002+ABX*I_ERI_Px_S_Pz_S_C1001002;
  Double I_ERI_Py_Px_Pz_S_C1001002 = I_ERI_Dxy_S_Pz_S_C1001002+ABX*I_ERI_Py_S_Pz_S_C1001002;
  Double I_ERI_Pz_Px_Pz_S_C1001002 = I_ERI_Dxz_S_Pz_S_C1001002+ABX*I_ERI_Pz_S_Pz_S_C1001002;
  Double I_ERI_Px_Py_Pz_S_C1001002 = I_ERI_Dxy_S_Pz_S_C1001002+ABY*I_ERI_Px_S_Pz_S_C1001002;
  Double I_ERI_Py_Py_Pz_S_C1001002 = I_ERI_D2y_S_Pz_S_C1001002+ABY*I_ERI_Py_S_Pz_S_C1001002;
  Double I_ERI_Pz_Py_Pz_S_C1001002 = I_ERI_Dyz_S_Pz_S_C1001002+ABY*I_ERI_Pz_S_Pz_S_C1001002;
  Double I_ERI_Px_Pz_Pz_S_C1001002 = I_ERI_Dxz_S_Pz_S_C1001002+ABZ*I_ERI_Px_S_Pz_S_C1001002;
  Double I_ERI_Py_Pz_Pz_S_C1001002 = I_ERI_Dyz_S_Pz_S_C1001002+ABZ*I_ERI_Py_S_Pz_S_C1001002;
  Double I_ERI_Pz_Pz_Pz_S_C1001002 = I_ERI_D2z_S_Pz_S_C1001002+ABZ*I_ERI_Pz_S_Pz_S_C1001002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1001002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1001002
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1001002
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C1001002 = I_ERI_F3x_S_S_S_C1001002+ABX*I_ERI_D2x_S_S_S_C1001002;
  Double I_ERI_Dxy_Px_S_S_C1001002 = I_ERI_F2xy_S_S_S_C1001002+ABX*I_ERI_Dxy_S_S_S_C1001002;
  Double I_ERI_Dxz_Px_S_S_C1001002 = I_ERI_F2xz_S_S_S_C1001002+ABX*I_ERI_Dxz_S_S_S_C1001002;
  Double I_ERI_D2y_Px_S_S_C1001002 = I_ERI_Fx2y_S_S_S_C1001002+ABX*I_ERI_D2y_S_S_S_C1001002;
  Double I_ERI_Dyz_Px_S_S_C1001002 = I_ERI_Fxyz_S_S_S_C1001002+ABX*I_ERI_Dyz_S_S_S_C1001002;
  Double I_ERI_D2z_Px_S_S_C1001002 = I_ERI_Fx2z_S_S_S_C1001002+ABX*I_ERI_D2z_S_S_S_C1001002;
  Double I_ERI_D2x_Py_S_S_C1001002 = I_ERI_F2xy_S_S_S_C1001002+ABY*I_ERI_D2x_S_S_S_C1001002;
  Double I_ERI_Dxy_Py_S_S_C1001002 = I_ERI_Fx2y_S_S_S_C1001002+ABY*I_ERI_Dxy_S_S_S_C1001002;
  Double I_ERI_Dxz_Py_S_S_C1001002 = I_ERI_Fxyz_S_S_S_C1001002+ABY*I_ERI_Dxz_S_S_S_C1001002;
  Double I_ERI_D2y_Py_S_S_C1001002 = I_ERI_F3y_S_S_S_C1001002+ABY*I_ERI_D2y_S_S_S_C1001002;
  Double I_ERI_Dyz_Py_S_S_C1001002 = I_ERI_F2yz_S_S_S_C1001002+ABY*I_ERI_Dyz_S_S_S_C1001002;
  Double I_ERI_D2z_Py_S_S_C1001002 = I_ERI_Fy2z_S_S_S_C1001002+ABY*I_ERI_D2z_S_S_S_C1001002;
  Double I_ERI_D2x_Pz_S_S_C1001002 = I_ERI_F2xz_S_S_S_C1001002+ABZ*I_ERI_D2x_S_S_S_C1001002;
  Double I_ERI_Dxy_Pz_S_S_C1001002 = I_ERI_Fxyz_S_S_S_C1001002+ABZ*I_ERI_Dxy_S_S_S_C1001002;
  Double I_ERI_Dxz_Pz_S_S_C1001002 = I_ERI_Fx2z_S_S_S_C1001002+ABZ*I_ERI_Dxz_S_S_S_C1001002;
  Double I_ERI_D2y_Pz_S_S_C1001002 = I_ERI_F2yz_S_S_S_C1001002+ABZ*I_ERI_D2y_S_S_S_C1001002;
  Double I_ERI_Dyz_Pz_S_S_C1001002 = I_ERI_Fy2z_S_S_S_C1001002+ABZ*I_ERI_Dyz_S_S_S_C1001002;
  Double I_ERI_D2z_Pz_S_S_C1001002 = I_ERI_F3z_S_S_S_C1001002+ABZ*I_ERI_D2z_S_S_S_C1001002;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_C1002_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_a
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1002_a = I_ERI_G4x_S_S_S_C1002_a+ABX*I_ERI_F3x_S_S_S_C1002_a;
  Double I_ERI_F2xy_Px_S_S_C1002_a = I_ERI_G3xy_S_S_S_C1002_a+ABX*I_ERI_F2xy_S_S_S_C1002_a;
  Double I_ERI_F2xz_Px_S_S_C1002_a = I_ERI_G3xz_S_S_S_C1002_a+ABX*I_ERI_F2xz_S_S_S_C1002_a;
  Double I_ERI_Fx2y_Px_S_S_C1002_a = I_ERI_G2x2y_S_S_S_C1002_a+ABX*I_ERI_Fx2y_S_S_S_C1002_a;
  Double I_ERI_Fxyz_Px_S_S_C1002_a = I_ERI_G2xyz_S_S_S_C1002_a+ABX*I_ERI_Fxyz_S_S_S_C1002_a;
  Double I_ERI_Fx2z_Px_S_S_C1002_a = I_ERI_G2x2z_S_S_S_C1002_a+ABX*I_ERI_Fx2z_S_S_S_C1002_a;
  Double I_ERI_F3y_Px_S_S_C1002_a = I_ERI_Gx3y_S_S_S_C1002_a+ABX*I_ERI_F3y_S_S_S_C1002_a;
  Double I_ERI_F2yz_Px_S_S_C1002_a = I_ERI_Gx2yz_S_S_S_C1002_a+ABX*I_ERI_F2yz_S_S_S_C1002_a;
  Double I_ERI_Fy2z_Px_S_S_C1002_a = I_ERI_Gxy2z_S_S_S_C1002_a+ABX*I_ERI_Fy2z_S_S_S_C1002_a;
  Double I_ERI_F3z_Px_S_S_C1002_a = I_ERI_Gx3z_S_S_S_C1002_a+ABX*I_ERI_F3z_S_S_S_C1002_a;
  Double I_ERI_F3x_Py_S_S_C1002_a = I_ERI_G3xy_S_S_S_C1002_a+ABY*I_ERI_F3x_S_S_S_C1002_a;
  Double I_ERI_F2xy_Py_S_S_C1002_a = I_ERI_G2x2y_S_S_S_C1002_a+ABY*I_ERI_F2xy_S_S_S_C1002_a;
  Double I_ERI_F2xz_Py_S_S_C1002_a = I_ERI_G2xyz_S_S_S_C1002_a+ABY*I_ERI_F2xz_S_S_S_C1002_a;
  Double I_ERI_Fx2y_Py_S_S_C1002_a = I_ERI_Gx3y_S_S_S_C1002_a+ABY*I_ERI_Fx2y_S_S_S_C1002_a;
  Double I_ERI_Fxyz_Py_S_S_C1002_a = I_ERI_Gx2yz_S_S_S_C1002_a+ABY*I_ERI_Fxyz_S_S_S_C1002_a;
  Double I_ERI_Fx2z_Py_S_S_C1002_a = I_ERI_Gxy2z_S_S_S_C1002_a+ABY*I_ERI_Fx2z_S_S_S_C1002_a;
  Double I_ERI_F3y_Py_S_S_C1002_a = I_ERI_G4y_S_S_S_C1002_a+ABY*I_ERI_F3y_S_S_S_C1002_a;
  Double I_ERI_F2yz_Py_S_S_C1002_a = I_ERI_G3yz_S_S_S_C1002_a+ABY*I_ERI_F2yz_S_S_S_C1002_a;
  Double I_ERI_Fy2z_Py_S_S_C1002_a = I_ERI_G2y2z_S_S_S_C1002_a+ABY*I_ERI_Fy2z_S_S_S_C1002_a;
  Double I_ERI_F3z_Py_S_S_C1002_a = I_ERI_Gy3z_S_S_S_C1002_a+ABY*I_ERI_F3z_S_S_S_C1002_a;
  Double I_ERI_F3x_Pz_S_S_C1002_a = I_ERI_G3xz_S_S_S_C1002_a+ABZ*I_ERI_F3x_S_S_S_C1002_a;
  Double I_ERI_F2xy_Pz_S_S_C1002_a = I_ERI_G2xyz_S_S_S_C1002_a+ABZ*I_ERI_F2xy_S_S_S_C1002_a;
  Double I_ERI_F2xz_Pz_S_S_C1002_a = I_ERI_G2x2z_S_S_S_C1002_a+ABZ*I_ERI_F2xz_S_S_S_C1002_a;
  Double I_ERI_Fx2y_Pz_S_S_C1002_a = I_ERI_Gx2yz_S_S_S_C1002_a+ABZ*I_ERI_Fx2y_S_S_S_C1002_a;
  Double I_ERI_Fxyz_Pz_S_S_C1002_a = I_ERI_Gxy2z_S_S_S_C1002_a+ABZ*I_ERI_Fxyz_S_S_S_C1002_a;
  Double I_ERI_Fx2z_Pz_S_S_C1002_a = I_ERI_Gx3z_S_S_S_C1002_a+ABZ*I_ERI_Fx2z_S_S_S_C1002_a;
  Double I_ERI_F3y_Pz_S_S_C1002_a = I_ERI_G3yz_S_S_S_C1002_a+ABZ*I_ERI_F3y_S_S_S_C1002_a;
  Double I_ERI_F2yz_Pz_S_S_C1002_a = I_ERI_G2y2z_S_S_S_C1002_a+ABZ*I_ERI_F2yz_S_S_S_C1002_a;
  Double I_ERI_Fy2z_Pz_S_S_C1002_a = I_ERI_Gy3z_S_S_S_C1002_a+ABZ*I_ERI_Fy2z_S_S_S_C1002_a;
  Double I_ERI_F3z_Pz_S_S_C1002_a = I_ERI_G4z_S_S_S_C1002_a+ABZ*I_ERI_F3z_S_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_C1001002_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1001002_a
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001002_a
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1001002_a = I_ERI_G4x_S_Px_S_C1001002_a+ABX*I_ERI_F3x_S_Px_S_C1001002_a;
  Double I_ERI_F2xy_Px_Px_S_C1001002_a = I_ERI_G3xy_S_Px_S_C1001002_a+ABX*I_ERI_F2xy_S_Px_S_C1001002_a;
  Double I_ERI_F2xz_Px_Px_S_C1001002_a = I_ERI_G3xz_S_Px_S_C1001002_a+ABX*I_ERI_F2xz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2y_Px_Px_S_C1001002_a = I_ERI_G2x2y_S_Px_S_C1001002_a+ABX*I_ERI_Fx2y_S_Px_S_C1001002_a;
  Double I_ERI_Fxyz_Px_Px_S_C1001002_a = I_ERI_G2xyz_S_Px_S_C1001002_a+ABX*I_ERI_Fxyz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2z_Px_Px_S_C1001002_a = I_ERI_G2x2z_S_Px_S_C1001002_a+ABX*I_ERI_Fx2z_S_Px_S_C1001002_a;
  Double I_ERI_F3y_Px_Px_S_C1001002_a = I_ERI_Gx3y_S_Px_S_C1001002_a+ABX*I_ERI_F3y_S_Px_S_C1001002_a;
  Double I_ERI_F2yz_Px_Px_S_C1001002_a = I_ERI_Gx2yz_S_Px_S_C1001002_a+ABX*I_ERI_F2yz_S_Px_S_C1001002_a;
  Double I_ERI_Fy2z_Px_Px_S_C1001002_a = I_ERI_Gxy2z_S_Px_S_C1001002_a+ABX*I_ERI_Fy2z_S_Px_S_C1001002_a;
  Double I_ERI_F3z_Px_Px_S_C1001002_a = I_ERI_Gx3z_S_Px_S_C1001002_a+ABX*I_ERI_F3z_S_Px_S_C1001002_a;
  Double I_ERI_F3x_Py_Px_S_C1001002_a = I_ERI_G3xy_S_Px_S_C1001002_a+ABY*I_ERI_F3x_S_Px_S_C1001002_a;
  Double I_ERI_F2xy_Py_Px_S_C1001002_a = I_ERI_G2x2y_S_Px_S_C1001002_a+ABY*I_ERI_F2xy_S_Px_S_C1001002_a;
  Double I_ERI_F2xz_Py_Px_S_C1001002_a = I_ERI_G2xyz_S_Px_S_C1001002_a+ABY*I_ERI_F2xz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2y_Py_Px_S_C1001002_a = I_ERI_Gx3y_S_Px_S_C1001002_a+ABY*I_ERI_Fx2y_S_Px_S_C1001002_a;
  Double I_ERI_Fxyz_Py_Px_S_C1001002_a = I_ERI_Gx2yz_S_Px_S_C1001002_a+ABY*I_ERI_Fxyz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2z_Py_Px_S_C1001002_a = I_ERI_Gxy2z_S_Px_S_C1001002_a+ABY*I_ERI_Fx2z_S_Px_S_C1001002_a;
  Double I_ERI_F3y_Py_Px_S_C1001002_a = I_ERI_G4y_S_Px_S_C1001002_a+ABY*I_ERI_F3y_S_Px_S_C1001002_a;
  Double I_ERI_F2yz_Py_Px_S_C1001002_a = I_ERI_G3yz_S_Px_S_C1001002_a+ABY*I_ERI_F2yz_S_Px_S_C1001002_a;
  Double I_ERI_Fy2z_Py_Px_S_C1001002_a = I_ERI_G2y2z_S_Px_S_C1001002_a+ABY*I_ERI_Fy2z_S_Px_S_C1001002_a;
  Double I_ERI_F3z_Py_Px_S_C1001002_a = I_ERI_Gy3z_S_Px_S_C1001002_a+ABY*I_ERI_F3z_S_Px_S_C1001002_a;
  Double I_ERI_F3x_Pz_Px_S_C1001002_a = I_ERI_G3xz_S_Px_S_C1001002_a+ABZ*I_ERI_F3x_S_Px_S_C1001002_a;
  Double I_ERI_F2xy_Pz_Px_S_C1001002_a = I_ERI_G2xyz_S_Px_S_C1001002_a+ABZ*I_ERI_F2xy_S_Px_S_C1001002_a;
  Double I_ERI_F2xz_Pz_Px_S_C1001002_a = I_ERI_G2x2z_S_Px_S_C1001002_a+ABZ*I_ERI_F2xz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2y_Pz_Px_S_C1001002_a = I_ERI_Gx2yz_S_Px_S_C1001002_a+ABZ*I_ERI_Fx2y_S_Px_S_C1001002_a;
  Double I_ERI_Fxyz_Pz_Px_S_C1001002_a = I_ERI_Gxy2z_S_Px_S_C1001002_a+ABZ*I_ERI_Fxyz_S_Px_S_C1001002_a;
  Double I_ERI_Fx2z_Pz_Px_S_C1001002_a = I_ERI_Gx3z_S_Px_S_C1001002_a+ABZ*I_ERI_Fx2z_S_Px_S_C1001002_a;
  Double I_ERI_F3y_Pz_Px_S_C1001002_a = I_ERI_G3yz_S_Px_S_C1001002_a+ABZ*I_ERI_F3y_S_Px_S_C1001002_a;
  Double I_ERI_F2yz_Pz_Px_S_C1001002_a = I_ERI_G2y2z_S_Px_S_C1001002_a+ABZ*I_ERI_F2yz_S_Px_S_C1001002_a;
  Double I_ERI_Fy2z_Pz_Px_S_C1001002_a = I_ERI_Gy3z_S_Px_S_C1001002_a+ABZ*I_ERI_Fy2z_S_Px_S_C1001002_a;
  Double I_ERI_F3z_Pz_Px_S_C1001002_a = I_ERI_G4z_S_Px_S_C1001002_a+ABZ*I_ERI_F3z_S_Px_S_C1001002_a;
  Double I_ERI_F3x_Px_Py_S_C1001002_a = I_ERI_G4x_S_Py_S_C1001002_a+ABX*I_ERI_F3x_S_Py_S_C1001002_a;
  Double I_ERI_F2xy_Px_Py_S_C1001002_a = I_ERI_G3xy_S_Py_S_C1001002_a+ABX*I_ERI_F2xy_S_Py_S_C1001002_a;
  Double I_ERI_F2xz_Px_Py_S_C1001002_a = I_ERI_G3xz_S_Py_S_C1001002_a+ABX*I_ERI_F2xz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2y_Px_Py_S_C1001002_a = I_ERI_G2x2y_S_Py_S_C1001002_a+ABX*I_ERI_Fx2y_S_Py_S_C1001002_a;
  Double I_ERI_Fxyz_Px_Py_S_C1001002_a = I_ERI_G2xyz_S_Py_S_C1001002_a+ABX*I_ERI_Fxyz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2z_Px_Py_S_C1001002_a = I_ERI_G2x2z_S_Py_S_C1001002_a+ABX*I_ERI_Fx2z_S_Py_S_C1001002_a;
  Double I_ERI_F3y_Px_Py_S_C1001002_a = I_ERI_Gx3y_S_Py_S_C1001002_a+ABX*I_ERI_F3y_S_Py_S_C1001002_a;
  Double I_ERI_F2yz_Px_Py_S_C1001002_a = I_ERI_Gx2yz_S_Py_S_C1001002_a+ABX*I_ERI_F2yz_S_Py_S_C1001002_a;
  Double I_ERI_Fy2z_Px_Py_S_C1001002_a = I_ERI_Gxy2z_S_Py_S_C1001002_a+ABX*I_ERI_Fy2z_S_Py_S_C1001002_a;
  Double I_ERI_F3z_Px_Py_S_C1001002_a = I_ERI_Gx3z_S_Py_S_C1001002_a+ABX*I_ERI_F3z_S_Py_S_C1001002_a;
  Double I_ERI_F3x_Py_Py_S_C1001002_a = I_ERI_G3xy_S_Py_S_C1001002_a+ABY*I_ERI_F3x_S_Py_S_C1001002_a;
  Double I_ERI_F2xy_Py_Py_S_C1001002_a = I_ERI_G2x2y_S_Py_S_C1001002_a+ABY*I_ERI_F2xy_S_Py_S_C1001002_a;
  Double I_ERI_F2xz_Py_Py_S_C1001002_a = I_ERI_G2xyz_S_Py_S_C1001002_a+ABY*I_ERI_F2xz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2y_Py_Py_S_C1001002_a = I_ERI_Gx3y_S_Py_S_C1001002_a+ABY*I_ERI_Fx2y_S_Py_S_C1001002_a;
  Double I_ERI_Fxyz_Py_Py_S_C1001002_a = I_ERI_Gx2yz_S_Py_S_C1001002_a+ABY*I_ERI_Fxyz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2z_Py_Py_S_C1001002_a = I_ERI_Gxy2z_S_Py_S_C1001002_a+ABY*I_ERI_Fx2z_S_Py_S_C1001002_a;
  Double I_ERI_F3y_Py_Py_S_C1001002_a = I_ERI_G4y_S_Py_S_C1001002_a+ABY*I_ERI_F3y_S_Py_S_C1001002_a;
  Double I_ERI_F2yz_Py_Py_S_C1001002_a = I_ERI_G3yz_S_Py_S_C1001002_a+ABY*I_ERI_F2yz_S_Py_S_C1001002_a;
  Double I_ERI_Fy2z_Py_Py_S_C1001002_a = I_ERI_G2y2z_S_Py_S_C1001002_a+ABY*I_ERI_Fy2z_S_Py_S_C1001002_a;
  Double I_ERI_F3z_Py_Py_S_C1001002_a = I_ERI_Gy3z_S_Py_S_C1001002_a+ABY*I_ERI_F3z_S_Py_S_C1001002_a;
  Double I_ERI_F3x_Pz_Py_S_C1001002_a = I_ERI_G3xz_S_Py_S_C1001002_a+ABZ*I_ERI_F3x_S_Py_S_C1001002_a;
  Double I_ERI_F2xy_Pz_Py_S_C1001002_a = I_ERI_G2xyz_S_Py_S_C1001002_a+ABZ*I_ERI_F2xy_S_Py_S_C1001002_a;
  Double I_ERI_F2xz_Pz_Py_S_C1001002_a = I_ERI_G2x2z_S_Py_S_C1001002_a+ABZ*I_ERI_F2xz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2y_Pz_Py_S_C1001002_a = I_ERI_Gx2yz_S_Py_S_C1001002_a+ABZ*I_ERI_Fx2y_S_Py_S_C1001002_a;
  Double I_ERI_Fxyz_Pz_Py_S_C1001002_a = I_ERI_Gxy2z_S_Py_S_C1001002_a+ABZ*I_ERI_Fxyz_S_Py_S_C1001002_a;
  Double I_ERI_Fx2z_Pz_Py_S_C1001002_a = I_ERI_Gx3z_S_Py_S_C1001002_a+ABZ*I_ERI_Fx2z_S_Py_S_C1001002_a;
  Double I_ERI_F3y_Pz_Py_S_C1001002_a = I_ERI_G3yz_S_Py_S_C1001002_a+ABZ*I_ERI_F3y_S_Py_S_C1001002_a;
  Double I_ERI_F2yz_Pz_Py_S_C1001002_a = I_ERI_G2y2z_S_Py_S_C1001002_a+ABZ*I_ERI_F2yz_S_Py_S_C1001002_a;
  Double I_ERI_Fy2z_Pz_Py_S_C1001002_a = I_ERI_Gy3z_S_Py_S_C1001002_a+ABZ*I_ERI_Fy2z_S_Py_S_C1001002_a;
  Double I_ERI_F3z_Pz_Py_S_C1001002_a = I_ERI_G4z_S_Py_S_C1001002_a+ABZ*I_ERI_F3z_S_Py_S_C1001002_a;
  Double I_ERI_F3x_Px_Pz_S_C1001002_a = I_ERI_G4x_S_Pz_S_C1001002_a+ABX*I_ERI_F3x_S_Pz_S_C1001002_a;
  Double I_ERI_F2xy_Px_Pz_S_C1001002_a = I_ERI_G3xy_S_Pz_S_C1001002_a+ABX*I_ERI_F2xy_S_Pz_S_C1001002_a;
  Double I_ERI_F2xz_Px_Pz_S_C1001002_a = I_ERI_G3xz_S_Pz_S_C1001002_a+ABX*I_ERI_F2xz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2y_Px_Pz_S_C1001002_a = I_ERI_G2x2y_S_Pz_S_C1001002_a+ABX*I_ERI_Fx2y_S_Pz_S_C1001002_a;
  Double I_ERI_Fxyz_Px_Pz_S_C1001002_a = I_ERI_G2xyz_S_Pz_S_C1001002_a+ABX*I_ERI_Fxyz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2z_Px_Pz_S_C1001002_a = I_ERI_G2x2z_S_Pz_S_C1001002_a+ABX*I_ERI_Fx2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3y_Px_Pz_S_C1001002_a = I_ERI_Gx3y_S_Pz_S_C1001002_a+ABX*I_ERI_F3y_S_Pz_S_C1001002_a;
  Double I_ERI_F2yz_Px_Pz_S_C1001002_a = I_ERI_Gx2yz_S_Pz_S_C1001002_a+ABX*I_ERI_F2yz_S_Pz_S_C1001002_a;
  Double I_ERI_Fy2z_Px_Pz_S_C1001002_a = I_ERI_Gxy2z_S_Pz_S_C1001002_a+ABX*I_ERI_Fy2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3z_Px_Pz_S_C1001002_a = I_ERI_Gx3z_S_Pz_S_C1001002_a+ABX*I_ERI_F3z_S_Pz_S_C1001002_a;
  Double I_ERI_F3x_Py_Pz_S_C1001002_a = I_ERI_G3xy_S_Pz_S_C1001002_a+ABY*I_ERI_F3x_S_Pz_S_C1001002_a;
  Double I_ERI_F2xy_Py_Pz_S_C1001002_a = I_ERI_G2x2y_S_Pz_S_C1001002_a+ABY*I_ERI_F2xy_S_Pz_S_C1001002_a;
  Double I_ERI_F2xz_Py_Pz_S_C1001002_a = I_ERI_G2xyz_S_Pz_S_C1001002_a+ABY*I_ERI_F2xz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2y_Py_Pz_S_C1001002_a = I_ERI_Gx3y_S_Pz_S_C1001002_a+ABY*I_ERI_Fx2y_S_Pz_S_C1001002_a;
  Double I_ERI_Fxyz_Py_Pz_S_C1001002_a = I_ERI_Gx2yz_S_Pz_S_C1001002_a+ABY*I_ERI_Fxyz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2z_Py_Pz_S_C1001002_a = I_ERI_Gxy2z_S_Pz_S_C1001002_a+ABY*I_ERI_Fx2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3y_Py_Pz_S_C1001002_a = I_ERI_G4y_S_Pz_S_C1001002_a+ABY*I_ERI_F3y_S_Pz_S_C1001002_a;
  Double I_ERI_F2yz_Py_Pz_S_C1001002_a = I_ERI_G3yz_S_Pz_S_C1001002_a+ABY*I_ERI_F2yz_S_Pz_S_C1001002_a;
  Double I_ERI_Fy2z_Py_Pz_S_C1001002_a = I_ERI_G2y2z_S_Pz_S_C1001002_a+ABY*I_ERI_Fy2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3z_Py_Pz_S_C1001002_a = I_ERI_Gy3z_S_Pz_S_C1001002_a+ABY*I_ERI_F3z_S_Pz_S_C1001002_a;
  Double I_ERI_F3x_Pz_Pz_S_C1001002_a = I_ERI_G3xz_S_Pz_S_C1001002_a+ABZ*I_ERI_F3x_S_Pz_S_C1001002_a;
  Double I_ERI_F2xy_Pz_Pz_S_C1001002_a = I_ERI_G2xyz_S_Pz_S_C1001002_a+ABZ*I_ERI_F2xy_S_Pz_S_C1001002_a;
  Double I_ERI_F2xz_Pz_Pz_S_C1001002_a = I_ERI_G2x2z_S_Pz_S_C1001002_a+ABZ*I_ERI_F2xz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2y_Pz_Pz_S_C1001002_a = I_ERI_Gx2yz_S_Pz_S_C1001002_a+ABZ*I_ERI_Fx2y_S_Pz_S_C1001002_a;
  Double I_ERI_Fxyz_Pz_Pz_S_C1001002_a = I_ERI_Gxy2z_S_Pz_S_C1001002_a+ABZ*I_ERI_Fxyz_S_Pz_S_C1001002_a;
  Double I_ERI_Fx2z_Pz_Pz_S_C1001002_a = I_ERI_Gx3z_S_Pz_S_C1001002_a+ABZ*I_ERI_Fx2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3y_Pz_Pz_S_C1001002_a = I_ERI_G3yz_S_Pz_S_C1001002_a+ABZ*I_ERI_F3y_S_Pz_S_C1001002_a;
  Double I_ERI_F2yz_Pz_Pz_S_C1001002_a = I_ERI_G2y2z_S_Pz_S_C1001002_a+ABZ*I_ERI_F2yz_S_Pz_S_C1001002_a;
  Double I_ERI_Fy2z_Pz_Pz_S_C1001002_a = I_ERI_Gy3z_S_Pz_S_C1001002_a+ABZ*I_ERI_Fy2z_S_Pz_S_C1001002_a;
  Double I_ERI_F3z_Pz_Pz_S_C1001002_a = I_ERI_G4z_S_Pz_S_C1001002_a+ABZ*I_ERI_F3z_S_Pz_S_C1001002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C2_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C2_b
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_C2_b = I_ERI_F3x_S_S_S_C2_b+ABX*I_ERI_D2x_S_S_S_C2_b;
  Double I_ERI_Dxy_Px_S_S_C2_b = I_ERI_F2xy_S_S_S_C2_b+ABX*I_ERI_Dxy_S_S_S_C2_b;
  Double I_ERI_Dxz_Px_S_S_C2_b = I_ERI_F2xz_S_S_S_C2_b+ABX*I_ERI_Dxz_S_S_S_C2_b;
  Double I_ERI_D2y_Px_S_S_C2_b = I_ERI_Fx2y_S_S_S_C2_b+ABX*I_ERI_D2y_S_S_S_C2_b;
  Double I_ERI_Dyz_Px_S_S_C2_b = I_ERI_Fxyz_S_S_S_C2_b+ABX*I_ERI_Dyz_S_S_S_C2_b;
  Double I_ERI_D2z_Px_S_S_C2_b = I_ERI_Fx2z_S_S_S_C2_b+ABX*I_ERI_D2z_S_S_S_C2_b;
  Double I_ERI_D2x_Py_S_S_C2_b = I_ERI_F2xy_S_S_S_C2_b+ABY*I_ERI_D2x_S_S_S_C2_b;
  Double I_ERI_Dxy_Py_S_S_C2_b = I_ERI_Fx2y_S_S_S_C2_b+ABY*I_ERI_Dxy_S_S_S_C2_b;
  Double I_ERI_Dxz_Py_S_S_C2_b = I_ERI_Fxyz_S_S_S_C2_b+ABY*I_ERI_Dxz_S_S_S_C2_b;
  Double I_ERI_D2y_Py_S_S_C2_b = I_ERI_F3y_S_S_S_C2_b+ABY*I_ERI_D2y_S_S_S_C2_b;
  Double I_ERI_Dyz_Py_S_S_C2_b = I_ERI_F2yz_S_S_S_C2_b+ABY*I_ERI_Dyz_S_S_S_C2_b;
  Double I_ERI_D2z_Py_S_S_C2_b = I_ERI_Fy2z_S_S_S_C2_b+ABY*I_ERI_D2z_S_S_S_C2_b;
  Double I_ERI_D2x_Pz_S_S_C2_b = I_ERI_F2xz_S_S_S_C2_b+ABZ*I_ERI_D2x_S_S_S_C2_b;
  Double I_ERI_Dxy_Pz_S_S_C2_b = I_ERI_Fxyz_S_S_S_C2_b+ABZ*I_ERI_Dxy_S_S_S_C2_b;
  Double I_ERI_Dxz_Pz_S_S_C2_b = I_ERI_Fx2z_S_S_S_C2_b+ABZ*I_ERI_Dxz_S_S_S_C2_b;
  Double I_ERI_D2y_Pz_S_S_C2_b = I_ERI_F2yz_S_S_S_C2_b+ABZ*I_ERI_D2y_S_S_S_C2_b;
  Double I_ERI_Dyz_Pz_S_S_C2_b = I_ERI_Fy2z_S_S_S_C2_b+ABZ*I_ERI_Dyz_S_S_S_C2_b;
  Double I_ERI_D2z_Pz_S_S_C2_b = I_ERI_F3z_S_S_S_C2_b+ABZ*I_ERI_D2z_S_S_S_C2_b;

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
   * shell quartet name: SQ_ERI_F_P_S_S_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C1002_b
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_C1002_b = I_ERI_G4x_S_S_S_C1002_b+ABX*I_ERI_F3x_S_S_S_C1002_b;
  Double I_ERI_F2xy_Px_S_S_C1002_b = I_ERI_G3xy_S_S_S_C1002_b+ABX*I_ERI_F2xy_S_S_S_C1002_b;
  Double I_ERI_F2xz_Px_S_S_C1002_b = I_ERI_G3xz_S_S_S_C1002_b+ABX*I_ERI_F2xz_S_S_S_C1002_b;
  Double I_ERI_Fx2y_Px_S_S_C1002_b = I_ERI_G2x2y_S_S_S_C1002_b+ABX*I_ERI_Fx2y_S_S_S_C1002_b;
  Double I_ERI_Fxyz_Px_S_S_C1002_b = I_ERI_G2xyz_S_S_S_C1002_b+ABX*I_ERI_Fxyz_S_S_S_C1002_b;
  Double I_ERI_Fx2z_Px_S_S_C1002_b = I_ERI_G2x2z_S_S_S_C1002_b+ABX*I_ERI_Fx2z_S_S_S_C1002_b;
  Double I_ERI_F3y_Px_S_S_C1002_b = I_ERI_Gx3y_S_S_S_C1002_b+ABX*I_ERI_F3y_S_S_S_C1002_b;
  Double I_ERI_F2yz_Px_S_S_C1002_b = I_ERI_Gx2yz_S_S_S_C1002_b+ABX*I_ERI_F2yz_S_S_S_C1002_b;
  Double I_ERI_Fy2z_Px_S_S_C1002_b = I_ERI_Gxy2z_S_S_S_C1002_b+ABX*I_ERI_Fy2z_S_S_S_C1002_b;
  Double I_ERI_F3z_Px_S_S_C1002_b = I_ERI_Gx3z_S_S_S_C1002_b+ABX*I_ERI_F3z_S_S_S_C1002_b;
  Double I_ERI_F2xy_Py_S_S_C1002_b = I_ERI_G2x2y_S_S_S_C1002_b+ABY*I_ERI_F2xy_S_S_S_C1002_b;
  Double I_ERI_F2xz_Py_S_S_C1002_b = I_ERI_G2xyz_S_S_S_C1002_b+ABY*I_ERI_F2xz_S_S_S_C1002_b;
  Double I_ERI_Fx2y_Py_S_S_C1002_b = I_ERI_Gx3y_S_S_S_C1002_b+ABY*I_ERI_Fx2y_S_S_S_C1002_b;
  Double I_ERI_Fxyz_Py_S_S_C1002_b = I_ERI_Gx2yz_S_S_S_C1002_b+ABY*I_ERI_Fxyz_S_S_S_C1002_b;
  Double I_ERI_Fx2z_Py_S_S_C1002_b = I_ERI_Gxy2z_S_S_S_C1002_b+ABY*I_ERI_Fx2z_S_S_S_C1002_b;
  Double I_ERI_F3y_Py_S_S_C1002_b = I_ERI_G4y_S_S_S_C1002_b+ABY*I_ERI_F3y_S_S_S_C1002_b;
  Double I_ERI_F2yz_Py_S_S_C1002_b = I_ERI_G3yz_S_S_S_C1002_b+ABY*I_ERI_F2yz_S_S_S_C1002_b;
  Double I_ERI_Fy2z_Py_S_S_C1002_b = I_ERI_G2y2z_S_S_S_C1002_b+ABY*I_ERI_Fy2z_S_S_S_C1002_b;
  Double I_ERI_F3z_Py_S_S_C1002_b = I_ERI_Gy3z_S_S_S_C1002_b+ABY*I_ERI_F3z_S_S_S_C1002_b;
  Double I_ERI_F2xz_Pz_S_S_C1002_b = I_ERI_G2x2z_S_S_S_C1002_b+ABZ*I_ERI_F2xz_S_S_S_C1002_b;
  Double I_ERI_Fxyz_Pz_S_S_C1002_b = I_ERI_Gxy2z_S_S_S_C1002_b+ABZ*I_ERI_Fxyz_S_S_S_C1002_b;
  Double I_ERI_Fx2z_Pz_S_S_C1002_b = I_ERI_Gx3z_S_S_S_C1002_b+ABZ*I_ERI_Fx2z_S_S_S_C1002_b;
  Double I_ERI_F2yz_Pz_S_S_C1002_b = I_ERI_G2y2z_S_S_S_C1002_b+ABZ*I_ERI_F2yz_S_S_S_C1002_b;
  Double I_ERI_Fy2z_Pz_S_S_C1002_b = I_ERI_Gy3z_S_S_S_C1002_b+ABZ*I_ERI_Fy2z_S_S_S_C1002_b;
  Double I_ERI_F3z_Pz_S_S_C1002_b = I_ERI_G4z_S_S_S_C1002_b+ABZ*I_ERI_F3z_S_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1002_b
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_C1002_b = I_ERI_F3x_Px_S_S_C1002_b+ABX*I_ERI_D2x_Px_S_S_C1002_b;
  Double I_ERI_Dxy_D2x_S_S_C1002_b = I_ERI_F2xy_Px_S_S_C1002_b+ABX*I_ERI_Dxy_Px_S_S_C1002_b;
  Double I_ERI_Dxz_D2x_S_S_C1002_b = I_ERI_F2xz_Px_S_S_C1002_b+ABX*I_ERI_Dxz_Px_S_S_C1002_b;
  Double I_ERI_D2y_D2x_S_S_C1002_b = I_ERI_Fx2y_Px_S_S_C1002_b+ABX*I_ERI_D2y_Px_S_S_C1002_b;
  Double I_ERI_Dyz_D2x_S_S_C1002_b = I_ERI_Fxyz_Px_S_S_C1002_b+ABX*I_ERI_Dyz_Px_S_S_C1002_b;
  Double I_ERI_D2z_D2x_S_S_C1002_b = I_ERI_Fx2z_Px_S_S_C1002_b+ABX*I_ERI_D2z_Px_S_S_C1002_b;
  Double I_ERI_D2x_Dxy_S_S_C1002_b = I_ERI_F2xy_Px_S_S_C1002_b+ABY*I_ERI_D2x_Px_S_S_C1002_b;
  Double I_ERI_Dxy_Dxy_S_S_C1002_b = I_ERI_Fx2y_Px_S_S_C1002_b+ABY*I_ERI_Dxy_Px_S_S_C1002_b;
  Double I_ERI_Dxz_Dxy_S_S_C1002_b = I_ERI_Fxyz_Px_S_S_C1002_b+ABY*I_ERI_Dxz_Px_S_S_C1002_b;
  Double I_ERI_D2y_Dxy_S_S_C1002_b = I_ERI_F3y_Px_S_S_C1002_b+ABY*I_ERI_D2y_Px_S_S_C1002_b;
  Double I_ERI_Dyz_Dxy_S_S_C1002_b = I_ERI_F2yz_Px_S_S_C1002_b+ABY*I_ERI_Dyz_Px_S_S_C1002_b;
  Double I_ERI_D2z_Dxy_S_S_C1002_b = I_ERI_Fy2z_Px_S_S_C1002_b+ABY*I_ERI_D2z_Px_S_S_C1002_b;
  Double I_ERI_D2x_Dxz_S_S_C1002_b = I_ERI_F2xz_Px_S_S_C1002_b+ABZ*I_ERI_D2x_Px_S_S_C1002_b;
  Double I_ERI_Dxy_Dxz_S_S_C1002_b = I_ERI_Fxyz_Px_S_S_C1002_b+ABZ*I_ERI_Dxy_Px_S_S_C1002_b;
  Double I_ERI_Dxz_Dxz_S_S_C1002_b = I_ERI_Fx2z_Px_S_S_C1002_b+ABZ*I_ERI_Dxz_Px_S_S_C1002_b;
  Double I_ERI_D2y_Dxz_S_S_C1002_b = I_ERI_F2yz_Px_S_S_C1002_b+ABZ*I_ERI_D2y_Px_S_S_C1002_b;
  Double I_ERI_Dyz_Dxz_S_S_C1002_b = I_ERI_Fy2z_Px_S_S_C1002_b+ABZ*I_ERI_Dyz_Px_S_S_C1002_b;
  Double I_ERI_D2z_Dxz_S_S_C1002_b = I_ERI_F3z_Px_S_S_C1002_b+ABZ*I_ERI_D2z_Px_S_S_C1002_b;
  Double I_ERI_D2x_D2y_S_S_C1002_b = I_ERI_F2xy_Py_S_S_C1002_b+ABY*I_ERI_D2x_Py_S_S_C1002_b;
  Double I_ERI_Dxy_D2y_S_S_C1002_b = I_ERI_Fx2y_Py_S_S_C1002_b+ABY*I_ERI_Dxy_Py_S_S_C1002_b;
  Double I_ERI_Dxz_D2y_S_S_C1002_b = I_ERI_Fxyz_Py_S_S_C1002_b+ABY*I_ERI_Dxz_Py_S_S_C1002_b;
  Double I_ERI_D2y_D2y_S_S_C1002_b = I_ERI_F3y_Py_S_S_C1002_b+ABY*I_ERI_D2y_Py_S_S_C1002_b;
  Double I_ERI_Dyz_D2y_S_S_C1002_b = I_ERI_F2yz_Py_S_S_C1002_b+ABY*I_ERI_Dyz_Py_S_S_C1002_b;
  Double I_ERI_D2z_D2y_S_S_C1002_b = I_ERI_Fy2z_Py_S_S_C1002_b+ABY*I_ERI_D2z_Py_S_S_C1002_b;
  Double I_ERI_D2x_Dyz_S_S_C1002_b = I_ERI_F2xz_Py_S_S_C1002_b+ABZ*I_ERI_D2x_Py_S_S_C1002_b;
  Double I_ERI_Dxy_Dyz_S_S_C1002_b = I_ERI_Fxyz_Py_S_S_C1002_b+ABZ*I_ERI_Dxy_Py_S_S_C1002_b;
  Double I_ERI_Dxz_Dyz_S_S_C1002_b = I_ERI_Fx2z_Py_S_S_C1002_b+ABZ*I_ERI_Dxz_Py_S_S_C1002_b;
  Double I_ERI_D2y_Dyz_S_S_C1002_b = I_ERI_F2yz_Py_S_S_C1002_b+ABZ*I_ERI_D2y_Py_S_S_C1002_b;
  Double I_ERI_Dyz_Dyz_S_S_C1002_b = I_ERI_Fy2z_Py_S_S_C1002_b+ABZ*I_ERI_Dyz_Py_S_S_C1002_b;
  Double I_ERI_D2z_Dyz_S_S_C1002_b = I_ERI_F3z_Py_S_S_C1002_b+ABZ*I_ERI_D2z_Py_S_S_C1002_b;
  Double I_ERI_D2x_D2z_S_S_C1002_b = I_ERI_F2xz_Pz_S_S_C1002_b+ABZ*I_ERI_D2x_Pz_S_S_C1002_b;
  Double I_ERI_Dxy_D2z_S_S_C1002_b = I_ERI_Fxyz_Pz_S_S_C1002_b+ABZ*I_ERI_Dxy_Pz_S_S_C1002_b;
  Double I_ERI_Dxz_D2z_S_S_C1002_b = I_ERI_Fx2z_Pz_S_S_C1002_b+ABZ*I_ERI_Dxz_Pz_S_S_C1002_b;
  Double I_ERI_D2y_D2z_S_S_C1002_b = I_ERI_F2yz_Pz_S_S_C1002_b+ABZ*I_ERI_D2y_Pz_S_S_C1002_b;
  Double I_ERI_Dyz_D2z_S_S_C1002_b = I_ERI_Fy2z_Pz_S_S_C1002_b+ABZ*I_ERI_Dyz_Pz_S_S_C1002_b;
  Double I_ERI_D2z_D2z_S_S_C1002_b = I_ERI_F3z_Pz_S_S_C1002_b+ABZ*I_ERI_D2z_Pz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1000002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000002_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1000002_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1000002_b = I_ERI_F3x_S_Px_S_C1000002_b+ABX*I_ERI_D2x_S_Px_S_C1000002_b;
  Double I_ERI_Dxy_Px_Px_S_C1000002_b = I_ERI_F2xy_S_Px_S_C1000002_b+ABX*I_ERI_Dxy_S_Px_S_C1000002_b;
  Double I_ERI_Dxz_Px_Px_S_C1000002_b = I_ERI_F2xz_S_Px_S_C1000002_b+ABX*I_ERI_Dxz_S_Px_S_C1000002_b;
  Double I_ERI_D2y_Px_Px_S_C1000002_b = I_ERI_Fx2y_S_Px_S_C1000002_b+ABX*I_ERI_D2y_S_Px_S_C1000002_b;
  Double I_ERI_Dyz_Px_Px_S_C1000002_b = I_ERI_Fxyz_S_Px_S_C1000002_b+ABX*I_ERI_Dyz_S_Px_S_C1000002_b;
  Double I_ERI_D2z_Px_Px_S_C1000002_b = I_ERI_Fx2z_S_Px_S_C1000002_b+ABX*I_ERI_D2z_S_Px_S_C1000002_b;
  Double I_ERI_D2x_Py_Px_S_C1000002_b = I_ERI_F2xy_S_Px_S_C1000002_b+ABY*I_ERI_D2x_S_Px_S_C1000002_b;
  Double I_ERI_Dxy_Py_Px_S_C1000002_b = I_ERI_Fx2y_S_Px_S_C1000002_b+ABY*I_ERI_Dxy_S_Px_S_C1000002_b;
  Double I_ERI_Dxz_Py_Px_S_C1000002_b = I_ERI_Fxyz_S_Px_S_C1000002_b+ABY*I_ERI_Dxz_S_Px_S_C1000002_b;
  Double I_ERI_D2y_Py_Px_S_C1000002_b = I_ERI_F3y_S_Px_S_C1000002_b+ABY*I_ERI_D2y_S_Px_S_C1000002_b;
  Double I_ERI_Dyz_Py_Px_S_C1000002_b = I_ERI_F2yz_S_Px_S_C1000002_b+ABY*I_ERI_Dyz_S_Px_S_C1000002_b;
  Double I_ERI_D2z_Py_Px_S_C1000002_b = I_ERI_Fy2z_S_Px_S_C1000002_b+ABY*I_ERI_D2z_S_Px_S_C1000002_b;
  Double I_ERI_D2x_Pz_Px_S_C1000002_b = I_ERI_F2xz_S_Px_S_C1000002_b+ABZ*I_ERI_D2x_S_Px_S_C1000002_b;
  Double I_ERI_Dxy_Pz_Px_S_C1000002_b = I_ERI_Fxyz_S_Px_S_C1000002_b+ABZ*I_ERI_Dxy_S_Px_S_C1000002_b;
  Double I_ERI_Dxz_Pz_Px_S_C1000002_b = I_ERI_Fx2z_S_Px_S_C1000002_b+ABZ*I_ERI_Dxz_S_Px_S_C1000002_b;
  Double I_ERI_D2y_Pz_Px_S_C1000002_b = I_ERI_F2yz_S_Px_S_C1000002_b+ABZ*I_ERI_D2y_S_Px_S_C1000002_b;
  Double I_ERI_Dyz_Pz_Px_S_C1000002_b = I_ERI_Fy2z_S_Px_S_C1000002_b+ABZ*I_ERI_Dyz_S_Px_S_C1000002_b;
  Double I_ERI_D2z_Pz_Px_S_C1000002_b = I_ERI_F3z_S_Px_S_C1000002_b+ABZ*I_ERI_D2z_S_Px_S_C1000002_b;
  Double I_ERI_D2x_Px_Py_S_C1000002_b = I_ERI_F3x_S_Py_S_C1000002_b+ABX*I_ERI_D2x_S_Py_S_C1000002_b;
  Double I_ERI_Dxy_Px_Py_S_C1000002_b = I_ERI_F2xy_S_Py_S_C1000002_b+ABX*I_ERI_Dxy_S_Py_S_C1000002_b;
  Double I_ERI_Dxz_Px_Py_S_C1000002_b = I_ERI_F2xz_S_Py_S_C1000002_b+ABX*I_ERI_Dxz_S_Py_S_C1000002_b;
  Double I_ERI_D2y_Px_Py_S_C1000002_b = I_ERI_Fx2y_S_Py_S_C1000002_b+ABX*I_ERI_D2y_S_Py_S_C1000002_b;
  Double I_ERI_Dyz_Px_Py_S_C1000002_b = I_ERI_Fxyz_S_Py_S_C1000002_b+ABX*I_ERI_Dyz_S_Py_S_C1000002_b;
  Double I_ERI_D2z_Px_Py_S_C1000002_b = I_ERI_Fx2z_S_Py_S_C1000002_b+ABX*I_ERI_D2z_S_Py_S_C1000002_b;
  Double I_ERI_D2x_Py_Py_S_C1000002_b = I_ERI_F2xy_S_Py_S_C1000002_b+ABY*I_ERI_D2x_S_Py_S_C1000002_b;
  Double I_ERI_Dxy_Py_Py_S_C1000002_b = I_ERI_Fx2y_S_Py_S_C1000002_b+ABY*I_ERI_Dxy_S_Py_S_C1000002_b;
  Double I_ERI_Dxz_Py_Py_S_C1000002_b = I_ERI_Fxyz_S_Py_S_C1000002_b+ABY*I_ERI_Dxz_S_Py_S_C1000002_b;
  Double I_ERI_D2y_Py_Py_S_C1000002_b = I_ERI_F3y_S_Py_S_C1000002_b+ABY*I_ERI_D2y_S_Py_S_C1000002_b;
  Double I_ERI_Dyz_Py_Py_S_C1000002_b = I_ERI_F2yz_S_Py_S_C1000002_b+ABY*I_ERI_Dyz_S_Py_S_C1000002_b;
  Double I_ERI_D2z_Py_Py_S_C1000002_b = I_ERI_Fy2z_S_Py_S_C1000002_b+ABY*I_ERI_D2z_S_Py_S_C1000002_b;
  Double I_ERI_D2x_Pz_Py_S_C1000002_b = I_ERI_F2xz_S_Py_S_C1000002_b+ABZ*I_ERI_D2x_S_Py_S_C1000002_b;
  Double I_ERI_Dxy_Pz_Py_S_C1000002_b = I_ERI_Fxyz_S_Py_S_C1000002_b+ABZ*I_ERI_Dxy_S_Py_S_C1000002_b;
  Double I_ERI_Dxz_Pz_Py_S_C1000002_b = I_ERI_Fx2z_S_Py_S_C1000002_b+ABZ*I_ERI_Dxz_S_Py_S_C1000002_b;
  Double I_ERI_D2y_Pz_Py_S_C1000002_b = I_ERI_F2yz_S_Py_S_C1000002_b+ABZ*I_ERI_D2y_S_Py_S_C1000002_b;
  Double I_ERI_Dyz_Pz_Py_S_C1000002_b = I_ERI_Fy2z_S_Py_S_C1000002_b+ABZ*I_ERI_Dyz_S_Py_S_C1000002_b;
  Double I_ERI_D2z_Pz_Py_S_C1000002_b = I_ERI_F3z_S_Py_S_C1000002_b+ABZ*I_ERI_D2z_S_Py_S_C1000002_b;
  Double I_ERI_D2x_Px_Pz_S_C1000002_b = I_ERI_F3x_S_Pz_S_C1000002_b+ABX*I_ERI_D2x_S_Pz_S_C1000002_b;
  Double I_ERI_Dxy_Px_Pz_S_C1000002_b = I_ERI_F2xy_S_Pz_S_C1000002_b+ABX*I_ERI_Dxy_S_Pz_S_C1000002_b;
  Double I_ERI_Dxz_Px_Pz_S_C1000002_b = I_ERI_F2xz_S_Pz_S_C1000002_b+ABX*I_ERI_Dxz_S_Pz_S_C1000002_b;
  Double I_ERI_D2y_Px_Pz_S_C1000002_b = I_ERI_Fx2y_S_Pz_S_C1000002_b+ABX*I_ERI_D2y_S_Pz_S_C1000002_b;
  Double I_ERI_Dyz_Px_Pz_S_C1000002_b = I_ERI_Fxyz_S_Pz_S_C1000002_b+ABX*I_ERI_Dyz_S_Pz_S_C1000002_b;
  Double I_ERI_D2z_Px_Pz_S_C1000002_b = I_ERI_Fx2z_S_Pz_S_C1000002_b+ABX*I_ERI_D2z_S_Pz_S_C1000002_b;
  Double I_ERI_D2x_Py_Pz_S_C1000002_b = I_ERI_F2xy_S_Pz_S_C1000002_b+ABY*I_ERI_D2x_S_Pz_S_C1000002_b;
  Double I_ERI_Dxy_Py_Pz_S_C1000002_b = I_ERI_Fx2y_S_Pz_S_C1000002_b+ABY*I_ERI_Dxy_S_Pz_S_C1000002_b;
  Double I_ERI_Dxz_Py_Pz_S_C1000002_b = I_ERI_Fxyz_S_Pz_S_C1000002_b+ABY*I_ERI_Dxz_S_Pz_S_C1000002_b;
  Double I_ERI_D2y_Py_Pz_S_C1000002_b = I_ERI_F3y_S_Pz_S_C1000002_b+ABY*I_ERI_D2y_S_Pz_S_C1000002_b;
  Double I_ERI_Dyz_Py_Pz_S_C1000002_b = I_ERI_F2yz_S_Pz_S_C1000002_b+ABY*I_ERI_Dyz_S_Pz_S_C1000002_b;
  Double I_ERI_D2z_Py_Pz_S_C1000002_b = I_ERI_Fy2z_S_Pz_S_C1000002_b+ABY*I_ERI_D2z_S_Pz_S_C1000002_b;
  Double I_ERI_D2x_Pz_Pz_S_C1000002_b = I_ERI_F2xz_S_Pz_S_C1000002_b+ABZ*I_ERI_D2x_S_Pz_S_C1000002_b;
  Double I_ERI_Dxy_Pz_Pz_S_C1000002_b = I_ERI_Fxyz_S_Pz_S_C1000002_b+ABZ*I_ERI_Dxy_S_Pz_S_C1000002_b;
  Double I_ERI_Dxz_Pz_Pz_S_C1000002_b = I_ERI_Fx2z_S_Pz_S_C1000002_b+ABZ*I_ERI_Dxz_S_Pz_S_C1000002_b;
  Double I_ERI_D2y_Pz_Pz_S_C1000002_b = I_ERI_F2yz_S_Pz_S_C1000002_b+ABZ*I_ERI_D2y_S_Pz_S_C1000002_b;
  Double I_ERI_Dyz_Pz_Pz_S_C1000002_b = I_ERI_Fy2z_S_Pz_S_C1000002_b+ABZ*I_ERI_Dyz_S_Pz_S_C1000002_b;
  Double I_ERI_D2z_Pz_Pz_S_C1000002_b = I_ERI_F3z_S_Pz_S_C1000002_b+ABZ*I_ERI_D2z_S_Pz_S_C1000002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001002_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1001002_b = I_ERI_F3x_S_Px_S_C1001002_b+ABX*I_ERI_D2x_S_Px_S_C1001002_b;
  Double I_ERI_Dxy_Px_Px_S_C1001002_b = I_ERI_F2xy_S_Px_S_C1001002_b+ABX*I_ERI_Dxy_S_Px_S_C1001002_b;
  Double I_ERI_Dxz_Px_Px_S_C1001002_b = I_ERI_F2xz_S_Px_S_C1001002_b+ABX*I_ERI_Dxz_S_Px_S_C1001002_b;
  Double I_ERI_D2y_Px_Px_S_C1001002_b = I_ERI_Fx2y_S_Px_S_C1001002_b+ABX*I_ERI_D2y_S_Px_S_C1001002_b;
  Double I_ERI_Dyz_Px_Px_S_C1001002_b = I_ERI_Fxyz_S_Px_S_C1001002_b+ABX*I_ERI_Dyz_S_Px_S_C1001002_b;
  Double I_ERI_D2z_Px_Px_S_C1001002_b = I_ERI_Fx2z_S_Px_S_C1001002_b+ABX*I_ERI_D2z_S_Px_S_C1001002_b;
  Double I_ERI_D2x_Py_Px_S_C1001002_b = I_ERI_F2xy_S_Px_S_C1001002_b+ABY*I_ERI_D2x_S_Px_S_C1001002_b;
  Double I_ERI_Dxy_Py_Px_S_C1001002_b = I_ERI_Fx2y_S_Px_S_C1001002_b+ABY*I_ERI_Dxy_S_Px_S_C1001002_b;
  Double I_ERI_Dxz_Py_Px_S_C1001002_b = I_ERI_Fxyz_S_Px_S_C1001002_b+ABY*I_ERI_Dxz_S_Px_S_C1001002_b;
  Double I_ERI_D2y_Py_Px_S_C1001002_b = I_ERI_F3y_S_Px_S_C1001002_b+ABY*I_ERI_D2y_S_Px_S_C1001002_b;
  Double I_ERI_Dyz_Py_Px_S_C1001002_b = I_ERI_F2yz_S_Px_S_C1001002_b+ABY*I_ERI_Dyz_S_Px_S_C1001002_b;
  Double I_ERI_D2z_Py_Px_S_C1001002_b = I_ERI_Fy2z_S_Px_S_C1001002_b+ABY*I_ERI_D2z_S_Px_S_C1001002_b;
  Double I_ERI_D2x_Pz_Px_S_C1001002_b = I_ERI_F2xz_S_Px_S_C1001002_b+ABZ*I_ERI_D2x_S_Px_S_C1001002_b;
  Double I_ERI_Dxy_Pz_Px_S_C1001002_b = I_ERI_Fxyz_S_Px_S_C1001002_b+ABZ*I_ERI_Dxy_S_Px_S_C1001002_b;
  Double I_ERI_Dxz_Pz_Px_S_C1001002_b = I_ERI_Fx2z_S_Px_S_C1001002_b+ABZ*I_ERI_Dxz_S_Px_S_C1001002_b;
  Double I_ERI_D2y_Pz_Px_S_C1001002_b = I_ERI_F2yz_S_Px_S_C1001002_b+ABZ*I_ERI_D2y_S_Px_S_C1001002_b;
  Double I_ERI_Dyz_Pz_Px_S_C1001002_b = I_ERI_Fy2z_S_Px_S_C1001002_b+ABZ*I_ERI_Dyz_S_Px_S_C1001002_b;
  Double I_ERI_D2z_Pz_Px_S_C1001002_b = I_ERI_F3z_S_Px_S_C1001002_b+ABZ*I_ERI_D2z_S_Px_S_C1001002_b;
  Double I_ERI_D2x_Px_Py_S_C1001002_b = I_ERI_F3x_S_Py_S_C1001002_b+ABX*I_ERI_D2x_S_Py_S_C1001002_b;
  Double I_ERI_Dxy_Px_Py_S_C1001002_b = I_ERI_F2xy_S_Py_S_C1001002_b+ABX*I_ERI_Dxy_S_Py_S_C1001002_b;
  Double I_ERI_Dxz_Px_Py_S_C1001002_b = I_ERI_F2xz_S_Py_S_C1001002_b+ABX*I_ERI_Dxz_S_Py_S_C1001002_b;
  Double I_ERI_D2y_Px_Py_S_C1001002_b = I_ERI_Fx2y_S_Py_S_C1001002_b+ABX*I_ERI_D2y_S_Py_S_C1001002_b;
  Double I_ERI_Dyz_Px_Py_S_C1001002_b = I_ERI_Fxyz_S_Py_S_C1001002_b+ABX*I_ERI_Dyz_S_Py_S_C1001002_b;
  Double I_ERI_D2z_Px_Py_S_C1001002_b = I_ERI_Fx2z_S_Py_S_C1001002_b+ABX*I_ERI_D2z_S_Py_S_C1001002_b;
  Double I_ERI_D2x_Py_Py_S_C1001002_b = I_ERI_F2xy_S_Py_S_C1001002_b+ABY*I_ERI_D2x_S_Py_S_C1001002_b;
  Double I_ERI_Dxy_Py_Py_S_C1001002_b = I_ERI_Fx2y_S_Py_S_C1001002_b+ABY*I_ERI_Dxy_S_Py_S_C1001002_b;
  Double I_ERI_Dxz_Py_Py_S_C1001002_b = I_ERI_Fxyz_S_Py_S_C1001002_b+ABY*I_ERI_Dxz_S_Py_S_C1001002_b;
  Double I_ERI_D2y_Py_Py_S_C1001002_b = I_ERI_F3y_S_Py_S_C1001002_b+ABY*I_ERI_D2y_S_Py_S_C1001002_b;
  Double I_ERI_Dyz_Py_Py_S_C1001002_b = I_ERI_F2yz_S_Py_S_C1001002_b+ABY*I_ERI_Dyz_S_Py_S_C1001002_b;
  Double I_ERI_D2z_Py_Py_S_C1001002_b = I_ERI_Fy2z_S_Py_S_C1001002_b+ABY*I_ERI_D2z_S_Py_S_C1001002_b;
  Double I_ERI_D2x_Pz_Py_S_C1001002_b = I_ERI_F2xz_S_Py_S_C1001002_b+ABZ*I_ERI_D2x_S_Py_S_C1001002_b;
  Double I_ERI_Dxy_Pz_Py_S_C1001002_b = I_ERI_Fxyz_S_Py_S_C1001002_b+ABZ*I_ERI_Dxy_S_Py_S_C1001002_b;
  Double I_ERI_Dxz_Pz_Py_S_C1001002_b = I_ERI_Fx2z_S_Py_S_C1001002_b+ABZ*I_ERI_Dxz_S_Py_S_C1001002_b;
  Double I_ERI_D2y_Pz_Py_S_C1001002_b = I_ERI_F2yz_S_Py_S_C1001002_b+ABZ*I_ERI_D2y_S_Py_S_C1001002_b;
  Double I_ERI_Dyz_Pz_Py_S_C1001002_b = I_ERI_Fy2z_S_Py_S_C1001002_b+ABZ*I_ERI_Dyz_S_Py_S_C1001002_b;
  Double I_ERI_D2z_Pz_Py_S_C1001002_b = I_ERI_F3z_S_Py_S_C1001002_b+ABZ*I_ERI_D2z_S_Py_S_C1001002_b;
  Double I_ERI_D2x_Px_Pz_S_C1001002_b = I_ERI_F3x_S_Pz_S_C1001002_b+ABX*I_ERI_D2x_S_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Px_Pz_S_C1001002_b = I_ERI_F2xy_S_Pz_S_C1001002_b+ABX*I_ERI_Dxy_S_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Px_Pz_S_C1001002_b = I_ERI_F2xz_S_Pz_S_C1001002_b+ABX*I_ERI_Dxz_S_Pz_S_C1001002_b;
  Double I_ERI_D2y_Px_Pz_S_C1001002_b = I_ERI_Fx2y_S_Pz_S_C1001002_b+ABX*I_ERI_D2y_S_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Px_Pz_S_C1001002_b = I_ERI_Fxyz_S_Pz_S_C1001002_b+ABX*I_ERI_Dyz_S_Pz_S_C1001002_b;
  Double I_ERI_D2z_Px_Pz_S_C1001002_b = I_ERI_Fx2z_S_Pz_S_C1001002_b+ABX*I_ERI_D2z_S_Pz_S_C1001002_b;
  Double I_ERI_D2x_Py_Pz_S_C1001002_b = I_ERI_F2xy_S_Pz_S_C1001002_b+ABY*I_ERI_D2x_S_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Py_Pz_S_C1001002_b = I_ERI_Fx2y_S_Pz_S_C1001002_b+ABY*I_ERI_Dxy_S_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Py_Pz_S_C1001002_b = I_ERI_Fxyz_S_Pz_S_C1001002_b+ABY*I_ERI_Dxz_S_Pz_S_C1001002_b;
  Double I_ERI_D2y_Py_Pz_S_C1001002_b = I_ERI_F3y_S_Pz_S_C1001002_b+ABY*I_ERI_D2y_S_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Py_Pz_S_C1001002_b = I_ERI_F2yz_S_Pz_S_C1001002_b+ABY*I_ERI_Dyz_S_Pz_S_C1001002_b;
  Double I_ERI_D2z_Py_Pz_S_C1001002_b = I_ERI_Fy2z_S_Pz_S_C1001002_b+ABY*I_ERI_D2z_S_Pz_S_C1001002_b;
  Double I_ERI_D2x_Pz_Pz_S_C1001002_b = I_ERI_F2xz_S_Pz_S_C1001002_b+ABZ*I_ERI_D2x_S_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Pz_Pz_S_C1001002_b = I_ERI_Fxyz_S_Pz_S_C1001002_b+ABZ*I_ERI_Dxy_S_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Pz_Pz_S_C1001002_b = I_ERI_Fx2z_S_Pz_S_C1001002_b+ABZ*I_ERI_Dxz_S_Pz_S_C1001002_b;
  Double I_ERI_D2y_Pz_Pz_S_C1001002_b = I_ERI_F2yz_S_Pz_S_C1001002_b+ABZ*I_ERI_D2y_S_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Pz_Pz_S_C1001002_b = I_ERI_Fy2z_S_Pz_S_C1001002_b+ABZ*I_ERI_Dyz_S_Pz_S_C1001002_b;
  Double I_ERI_D2z_Pz_Pz_S_C1001002_b = I_ERI_F3z_S_Pz_S_C1001002_b+ABZ*I_ERI_D2z_S_Pz_S_C1001002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_C1001002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1001002_b
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_C1001002_b = I_ERI_G4x_S_Px_S_C1001002_b+ABX*I_ERI_F3x_S_Px_S_C1001002_b;
  Double I_ERI_F2xy_Px_Px_S_C1001002_b = I_ERI_G3xy_S_Px_S_C1001002_b+ABX*I_ERI_F2xy_S_Px_S_C1001002_b;
  Double I_ERI_F2xz_Px_Px_S_C1001002_b = I_ERI_G3xz_S_Px_S_C1001002_b+ABX*I_ERI_F2xz_S_Px_S_C1001002_b;
  Double I_ERI_Fx2y_Px_Px_S_C1001002_b = I_ERI_G2x2y_S_Px_S_C1001002_b+ABX*I_ERI_Fx2y_S_Px_S_C1001002_b;
  Double I_ERI_Fxyz_Px_Px_S_C1001002_b = I_ERI_G2xyz_S_Px_S_C1001002_b+ABX*I_ERI_Fxyz_S_Px_S_C1001002_b;
  Double I_ERI_Fx2z_Px_Px_S_C1001002_b = I_ERI_G2x2z_S_Px_S_C1001002_b+ABX*I_ERI_Fx2z_S_Px_S_C1001002_b;
  Double I_ERI_F3y_Px_Px_S_C1001002_b = I_ERI_Gx3y_S_Px_S_C1001002_b+ABX*I_ERI_F3y_S_Px_S_C1001002_b;
  Double I_ERI_F2yz_Px_Px_S_C1001002_b = I_ERI_Gx2yz_S_Px_S_C1001002_b+ABX*I_ERI_F2yz_S_Px_S_C1001002_b;
  Double I_ERI_Fy2z_Px_Px_S_C1001002_b = I_ERI_Gxy2z_S_Px_S_C1001002_b+ABX*I_ERI_Fy2z_S_Px_S_C1001002_b;
  Double I_ERI_F3z_Px_Px_S_C1001002_b = I_ERI_Gx3z_S_Px_S_C1001002_b+ABX*I_ERI_F3z_S_Px_S_C1001002_b;
  Double I_ERI_F2xy_Py_Px_S_C1001002_b = I_ERI_G2x2y_S_Px_S_C1001002_b+ABY*I_ERI_F2xy_S_Px_S_C1001002_b;
  Double I_ERI_F2xz_Py_Px_S_C1001002_b = I_ERI_G2xyz_S_Px_S_C1001002_b+ABY*I_ERI_F2xz_S_Px_S_C1001002_b;
  Double I_ERI_Fx2y_Py_Px_S_C1001002_b = I_ERI_Gx3y_S_Px_S_C1001002_b+ABY*I_ERI_Fx2y_S_Px_S_C1001002_b;
  Double I_ERI_Fxyz_Py_Px_S_C1001002_b = I_ERI_Gx2yz_S_Px_S_C1001002_b+ABY*I_ERI_Fxyz_S_Px_S_C1001002_b;
  Double I_ERI_Fx2z_Py_Px_S_C1001002_b = I_ERI_Gxy2z_S_Px_S_C1001002_b+ABY*I_ERI_Fx2z_S_Px_S_C1001002_b;
  Double I_ERI_F3y_Py_Px_S_C1001002_b = I_ERI_G4y_S_Px_S_C1001002_b+ABY*I_ERI_F3y_S_Px_S_C1001002_b;
  Double I_ERI_F2yz_Py_Px_S_C1001002_b = I_ERI_G3yz_S_Px_S_C1001002_b+ABY*I_ERI_F2yz_S_Px_S_C1001002_b;
  Double I_ERI_Fy2z_Py_Px_S_C1001002_b = I_ERI_G2y2z_S_Px_S_C1001002_b+ABY*I_ERI_Fy2z_S_Px_S_C1001002_b;
  Double I_ERI_F3z_Py_Px_S_C1001002_b = I_ERI_Gy3z_S_Px_S_C1001002_b+ABY*I_ERI_F3z_S_Px_S_C1001002_b;
  Double I_ERI_F2xz_Pz_Px_S_C1001002_b = I_ERI_G2x2z_S_Px_S_C1001002_b+ABZ*I_ERI_F2xz_S_Px_S_C1001002_b;
  Double I_ERI_Fxyz_Pz_Px_S_C1001002_b = I_ERI_Gxy2z_S_Px_S_C1001002_b+ABZ*I_ERI_Fxyz_S_Px_S_C1001002_b;
  Double I_ERI_Fx2z_Pz_Px_S_C1001002_b = I_ERI_Gx3z_S_Px_S_C1001002_b+ABZ*I_ERI_Fx2z_S_Px_S_C1001002_b;
  Double I_ERI_F2yz_Pz_Px_S_C1001002_b = I_ERI_G2y2z_S_Px_S_C1001002_b+ABZ*I_ERI_F2yz_S_Px_S_C1001002_b;
  Double I_ERI_Fy2z_Pz_Px_S_C1001002_b = I_ERI_Gy3z_S_Px_S_C1001002_b+ABZ*I_ERI_Fy2z_S_Px_S_C1001002_b;
  Double I_ERI_F3z_Pz_Px_S_C1001002_b = I_ERI_G4z_S_Px_S_C1001002_b+ABZ*I_ERI_F3z_S_Px_S_C1001002_b;
  Double I_ERI_F3x_Px_Py_S_C1001002_b = I_ERI_G4x_S_Py_S_C1001002_b+ABX*I_ERI_F3x_S_Py_S_C1001002_b;
  Double I_ERI_F2xy_Px_Py_S_C1001002_b = I_ERI_G3xy_S_Py_S_C1001002_b+ABX*I_ERI_F2xy_S_Py_S_C1001002_b;
  Double I_ERI_F2xz_Px_Py_S_C1001002_b = I_ERI_G3xz_S_Py_S_C1001002_b+ABX*I_ERI_F2xz_S_Py_S_C1001002_b;
  Double I_ERI_Fx2y_Px_Py_S_C1001002_b = I_ERI_G2x2y_S_Py_S_C1001002_b+ABX*I_ERI_Fx2y_S_Py_S_C1001002_b;
  Double I_ERI_Fxyz_Px_Py_S_C1001002_b = I_ERI_G2xyz_S_Py_S_C1001002_b+ABX*I_ERI_Fxyz_S_Py_S_C1001002_b;
  Double I_ERI_Fx2z_Px_Py_S_C1001002_b = I_ERI_G2x2z_S_Py_S_C1001002_b+ABX*I_ERI_Fx2z_S_Py_S_C1001002_b;
  Double I_ERI_F3y_Px_Py_S_C1001002_b = I_ERI_Gx3y_S_Py_S_C1001002_b+ABX*I_ERI_F3y_S_Py_S_C1001002_b;
  Double I_ERI_F2yz_Px_Py_S_C1001002_b = I_ERI_Gx2yz_S_Py_S_C1001002_b+ABX*I_ERI_F2yz_S_Py_S_C1001002_b;
  Double I_ERI_Fy2z_Px_Py_S_C1001002_b = I_ERI_Gxy2z_S_Py_S_C1001002_b+ABX*I_ERI_Fy2z_S_Py_S_C1001002_b;
  Double I_ERI_F3z_Px_Py_S_C1001002_b = I_ERI_Gx3z_S_Py_S_C1001002_b+ABX*I_ERI_F3z_S_Py_S_C1001002_b;
  Double I_ERI_F2xy_Py_Py_S_C1001002_b = I_ERI_G2x2y_S_Py_S_C1001002_b+ABY*I_ERI_F2xy_S_Py_S_C1001002_b;
  Double I_ERI_F2xz_Py_Py_S_C1001002_b = I_ERI_G2xyz_S_Py_S_C1001002_b+ABY*I_ERI_F2xz_S_Py_S_C1001002_b;
  Double I_ERI_Fx2y_Py_Py_S_C1001002_b = I_ERI_Gx3y_S_Py_S_C1001002_b+ABY*I_ERI_Fx2y_S_Py_S_C1001002_b;
  Double I_ERI_Fxyz_Py_Py_S_C1001002_b = I_ERI_Gx2yz_S_Py_S_C1001002_b+ABY*I_ERI_Fxyz_S_Py_S_C1001002_b;
  Double I_ERI_Fx2z_Py_Py_S_C1001002_b = I_ERI_Gxy2z_S_Py_S_C1001002_b+ABY*I_ERI_Fx2z_S_Py_S_C1001002_b;
  Double I_ERI_F3y_Py_Py_S_C1001002_b = I_ERI_G4y_S_Py_S_C1001002_b+ABY*I_ERI_F3y_S_Py_S_C1001002_b;
  Double I_ERI_F2yz_Py_Py_S_C1001002_b = I_ERI_G3yz_S_Py_S_C1001002_b+ABY*I_ERI_F2yz_S_Py_S_C1001002_b;
  Double I_ERI_Fy2z_Py_Py_S_C1001002_b = I_ERI_G2y2z_S_Py_S_C1001002_b+ABY*I_ERI_Fy2z_S_Py_S_C1001002_b;
  Double I_ERI_F3z_Py_Py_S_C1001002_b = I_ERI_Gy3z_S_Py_S_C1001002_b+ABY*I_ERI_F3z_S_Py_S_C1001002_b;
  Double I_ERI_F2xz_Pz_Py_S_C1001002_b = I_ERI_G2x2z_S_Py_S_C1001002_b+ABZ*I_ERI_F2xz_S_Py_S_C1001002_b;
  Double I_ERI_Fxyz_Pz_Py_S_C1001002_b = I_ERI_Gxy2z_S_Py_S_C1001002_b+ABZ*I_ERI_Fxyz_S_Py_S_C1001002_b;
  Double I_ERI_Fx2z_Pz_Py_S_C1001002_b = I_ERI_Gx3z_S_Py_S_C1001002_b+ABZ*I_ERI_Fx2z_S_Py_S_C1001002_b;
  Double I_ERI_F2yz_Pz_Py_S_C1001002_b = I_ERI_G2y2z_S_Py_S_C1001002_b+ABZ*I_ERI_F2yz_S_Py_S_C1001002_b;
  Double I_ERI_Fy2z_Pz_Py_S_C1001002_b = I_ERI_Gy3z_S_Py_S_C1001002_b+ABZ*I_ERI_Fy2z_S_Py_S_C1001002_b;
  Double I_ERI_F3z_Pz_Py_S_C1001002_b = I_ERI_G4z_S_Py_S_C1001002_b+ABZ*I_ERI_F3z_S_Py_S_C1001002_b;
  Double I_ERI_F3x_Px_Pz_S_C1001002_b = I_ERI_G4x_S_Pz_S_C1001002_b+ABX*I_ERI_F3x_S_Pz_S_C1001002_b;
  Double I_ERI_F2xy_Px_Pz_S_C1001002_b = I_ERI_G3xy_S_Pz_S_C1001002_b+ABX*I_ERI_F2xy_S_Pz_S_C1001002_b;
  Double I_ERI_F2xz_Px_Pz_S_C1001002_b = I_ERI_G3xz_S_Pz_S_C1001002_b+ABX*I_ERI_F2xz_S_Pz_S_C1001002_b;
  Double I_ERI_Fx2y_Px_Pz_S_C1001002_b = I_ERI_G2x2y_S_Pz_S_C1001002_b+ABX*I_ERI_Fx2y_S_Pz_S_C1001002_b;
  Double I_ERI_Fxyz_Px_Pz_S_C1001002_b = I_ERI_G2xyz_S_Pz_S_C1001002_b+ABX*I_ERI_Fxyz_S_Pz_S_C1001002_b;
  Double I_ERI_Fx2z_Px_Pz_S_C1001002_b = I_ERI_G2x2z_S_Pz_S_C1001002_b+ABX*I_ERI_Fx2z_S_Pz_S_C1001002_b;
  Double I_ERI_F3y_Px_Pz_S_C1001002_b = I_ERI_Gx3y_S_Pz_S_C1001002_b+ABX*I_ERI_F3y_S_Pz_S_C1001002_b;
  Double I_ERI_F2yz_Px_Pz_S_C1001002_b = I_ERI_Gx2yz_S_Pz_S_C1001002_b+ABX*I_ERI_F2yz_S_Pz_S_C1001002_b;
  Double I_ERI_Fy2z_Px_Pz_S_C1001002_b = I_ERI_Gxy2z_S_Pz_S_C1001002_b+ABX*I_ERI_Fy2z_S_Pz_S_C1001002_b;
  Double I_ERI_F3z_Px_Pz_S_C1001002_b = I_ERI_Gx3z_S_Pz_S_C1001002_b+ABX*I_ERI_F3z_S_Pz_S_C1001002_b;
  Double I_ERI_F2xy_Py_Pz_S_C1001002_b = I_ERI_G2x2y_S_Pz_S_C1001002_b+ABY*I_ERI_F2xy_S_Pz_S_C1001002_b;
  Double I_ERI_F2xz_Py_Pz_S_C1001002_b = I_ERI_G2xyz_S_Pz_S_C1001002_b+ABY*I_ERI_F2xz_S_Pz_S_C1001002_b;
  Double I_ERI_Fx2y_Py_Pz_S_C1001002_b = I_ERI_Gx3y_S_Pz_S_C1001002_b+ABY*I_ERI_Fx2y_S_Pz_S_C1001002_b;
  Double I_ERI_Fxyz_Py_Pz_S_C1001002_b = I_ERI_Gx2yz_S_Pz_S_C1001002_b+ABY*I_ERI_Fxyz_S_Pz_S_C1001002_b;
  Double I_ERI_Fx2z_Py_Pz_S_C1001002_b = I_ERI_Gxy2z_S_Pz_S_C1001002_b+ABY*I_ERI_Fx2z_S_Pz_S_C1001002_b;
  Double I_ERI_F3y_Py_Pz_S_C1001002_b = I_ERI_G4y_S_Pz_S_C1001002_b+ABY*I_ERI_F3y_S_Pz_S_C1001002_b;
  Double I_ERI_F2yz_Py_Pz_S_C1001002_b = I_ERI_G3yz_S_Pz_S_C1001002_b+ABY*I_ERI_F2yz_S_Pz_S_C1001002_b;
  Double I_ERI_Fy2z_Py_Pz_S_C1001002_b = I_ERI_G2y2z_S_Pz_S_C1001002_b+ABY*I_ERI_Fy2z_S_Pz_S_C1001002_b;
  Double I_ERI_F3z_Py_Pz_S_C1001002_b = I_ERI_Gy3z_S_Pz_S_C1001002_b+ABY*I_ERI_F3z_S_Pz_S_C1001002_b;
  Double I_ERI_F2xz_Pz_Pz_S_C1001002_b = I_ERI_G2x2z_S_Pz_S_C1001002_b+ABZ*I_ERI_F2xz_S_Pz_S_C1001002_b;
  Double I_ERI_Fxyz_Pz_Pz_S_C1001002_b = I_ERI_Gxy2z_S_Pz_S_C1001002_b+ABZ*I_ERI_Fxyz_S_Pz_S_C1001002_b;
  Double I_ERI_Fx2z_Pz_Pz_S_C1001002_b = I_ERI_Gx3z_S_Pz_S_C1001002_b+ABZ*I_ERI_Fx2z_S_Pz_S_C1001002_b;
  Double I_ERI_F2yz_Pz_Pz_S_C1001002_b = I_ERI_G2y2z_S_Pz_S_C1001002_b+ABZ*I_ERI_F2yz_S_Pz_S_C1001002_b;
  Double I_ERI_Fy2z_Pz_Pz_S_C1001002_b = I_ERI_Gy3z_S_Pz_S_C1001002_b+ABZ*I_ERI_Fy2z_S_Pz_S_C1001002_b;
  Double I_ERI_F3z_Pz_Pz_S_C1001002_b = I_ERI_G4z_S_Pz_S_C1001002_b+ABZ*I_ERI_F3z_S_Pz_S_C1001002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_P_S_C1001002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1001002_b
   ************************************************************/
  Double I_ERI_D2x_D2x_Px_S_C1001002_b = I_ERI_F3x_Px_Px_S_C1001002_b+ABX*I_ERI_D2x_Px_Px_S_C1001002_b;
  Double I_ERI_Dxy_D2x_Px_S_C1001002_b = I_ERI_F2xy_Px_Px_S_C1001002_b+ABX*I_ERI_Dxy_Px_Px_S_C1001002_b;
  Double I_ERI_Dxz_D2x_Px_S_C1001002_b = I_ERI_F2xz_Px_Px_S_C1001002_b+ABX*I_ERI_Dxz_Px_Px_S_C1001002_b;
  Double I_ERI_D2y_D2x_Px_S_C1001002_b = I_ERI_Fx2y_Px_Px_S_C1001002_b+ABX*I_ERI_D2y_Px_Px_S_C1001002_b;
  Double I_ERI_Dyz_D2x_Px_S_C1001002_b = I_ERI_Fxyz_Px_Px_S_C1001002_b+ABX*I_ERI_Dyz_Px_Px_S_C1001002_b;
  Double I_ERI_D2z_D2x_Px_S_C1001002_b = I_ERI_Fx2z_Px_Px_S_C1001002_b+ABX*I_ERI_D2z_Px_Px_S_C1001002_b;
  Double I_ERI_D2x_Dxy_Px_S_C1001002_b = I_ERI_F2xy_Px_Px_S_C1001002_b+ABY*I_ERI_D2x_Px_Px_S_C1001002_b;
  Double I_ERI_Dxy_Dxy_Px_S_C1001002_b = I_ERI_Fx2y_Px_Px_S_C1001002_b+ABY*I_ERI_Dxy_Px_Px_S_C1001002_b;
  Double I_ERI_Dxz_Dxy_Px_S_C1001002_b = I_ERI_Fxyz_Px_Px_S_C1001002_b+ABY*I_ERI_Dxz_Px_Px_S_C1001002_b;
  Double I_ERI_D2y_Dxy_Px_S_C1001002_b = I_ERI_F3y_Px_Px_S_C1001002_b+ABY*I_ERI_D2y_Px_Px_S_C1001002_b;
  Double I_ERI_Dyz_Dxy_Px_S_C1001002_b = I_ERI_F2yz_Px_Px_S_C1001002_b+ABY*I_ERI_Dyz_Px_Px_S_C1001002_b;
  Double I_ERI_D2z_Dxy_Px_S_C1001002_b = I_ERI_Fy2z_Px_Px_S_C1001002_b+ABY*I_ERI_D2z_Px_Px_S_C1001002_b;
  Double I_ERI_D2x_Dxz_Px_S_C1001002_b = I_ERI_F2xz_Px_Px_S_C1001002_b+ABZ*I_ERI_D2x_Px_Px_S_C1001002_b;
  Double I_ERI_Dxy_Dxz_Px_S_C1001002_b = I_ERI_Fxyz_Px_Px_S_C1001002_b+ABZ*I_ERI_Dxy_Px_Px_S_C1001002_b;
  Double I_ERI_Dxz_Dxz_Px_S_C1001002_b = I_ERI_Fx2z_Px_Px_S_C1001002_b+ABZ*I_ERI_Dxz_Px_Px_S_C1001002_b;
  Double I_ERI_D2y_Dxz_Px_S_C1001002_b = I_ERI_F2yz_Px_Px_S_C1001002_b+ABZ*I_ERI_D2y_Px_Px_S_C1001002_b;
  Double I_ERI_Dyz_Dxz_Px_S_C1001002_b = I_ERI_Fy2z_Px_Px_S_C1001002_b+ABZ*I_ERI_Dyz_Px_Px_S_C1001002_b;
  Double I_ERI_D2z_Dxz_Px_S_C1001002_b = I_ERI_F3z_Px_Px_S_C1001002_b+ABZ*I_ERI_D2z_Px_Px_S_C1001002_b;
  Double I_ERI_D2x_D2y_Px_S_C1001002_b = I_ERI_F2xy_Py_Px_S_C1001002_b+ABY*I_ERI_D2x_Py_Px_S_C1001002_b;
  Double I_ERI_Dxy_D2y_Px_S_C1001002_b = I_ERI_Fx2y_Py_Px_S_C1001002_b+ABY*I_ERI_Dxy_Py_Px_S_C1001002_b;
  Double I_ERI_Dxz_D2y_Px_S_C1001002_b = I_ERI_Fxyz_Py_Px_S_C1001002_b+ABY*I_ERI_Dxz_Py_Px_S_C1001002_b;
  Double I_ERI_D2y_D2y_Px_S_C1001002_b = I_ERI_F3y_Py_Px_S_C1001002_b+ABY*I_ERI_D2y_Py_Px_S_C1001002_b;
  Double I_ERI_Dyz_D2y_Px_S_C1001002_b = I_ERI_F2yz_Py_Px_S_C1001002_b+ABY*I_ERI_Dyz_Py_Px_S_C1001002_b;
  Double I_ERI_D2z_D2y_Px_S_C1001002_b = I_ERI_Fy2z_Py_Px_S_C1001002_b+ABY*I_ERI_D2z_Py_Px_S_C1001002_b;
  Double I_ERI_D2x_Dyz_Px_S_C1001002_b = I_ERI_F2xz_Py_Px_S_C1001002_b+ABZ*I_ERI_D2x_Py_Px_S_C1001002_b;
  Double I_ERI_Dxy_Dyz_Px_S_C1001002_b = I_ERI_Fxyz_Py_Px_S_C1001002_b+ABZ*I_ERI_Dxy_Py_Px_S_C1001002_b;
  Double I_ERI_Dxz_Dyz_Px_S_C1001002_b = I_ERI_Fx2z_Py_Px_S_C1001002_b+ABZ*I_ERI_Dxz_Py_Px_S_C1001002_b;
  Double I_ERI_D2y_Dyz_Px_S_C1001002_b = I_ERI_F2yz_Py_Px_S_C1001002_b+ABZ*I_ERI_D2y_Py_Px_S_C1001002_b;
  Double I_ERI_Dyz_Dyz_Px_S_C1001002_b = I_ERI_Fy2z_Py_Px_S_C1001002_b+ABZ*I_ERI_Dyz_Py_Px_S_C1001002_b;
  Double I_ERI_D2z_Dyz_Px_S_C1001002_b = I_ERI_F3z_Py_Px_S_C1001002_b+ABZ*I_ERI_D2z_Py_Px_S_C1001002_b;
  Double I_ERI_D2x_D2z_Px_S_C1001002_b = I_ERI_F2xz_Pz_Px_S_C1001002_b+ABZ*I_ERI_D2x_Pz_Px_S_C1001002_b;
  Double I_ERI_Dxy_D2z_Px_S_C1001002_b = I_ERI_Fxyz_Pz_Px_S_C1001002_b+ABZ*I_ERI_Dxy_Pz_Px_S_C1001002_b;
  Double I_ERI_Dxz_D2z_Px_S_C1001002_b = I_ERI_Fx2z_Pz_Px_S_C1001002_b+ABZ*I_ERI_Dxz_Pz_Px_S_C1001002_b;
  Double I_ERI_D2y_D2z_Px_S_C1001002_b = I_ERI_F2yz_Pz_Px_S_C1001002_b+ABZ*I_ERI_D2y_Pz_Px_S_C1001002_b;
  Double I_ERI_Dyz_D2z_Px_S_C1001002_b = I_ERI_Fy2z_Pz_Px_S_C1001002_b+ABZ*I_ERI_Dyz_Pz_Px_S_C1001002_b;
  Double I_ERI_D2z_D2z_Px_S_C1001002_b = I_ERI_F3z_Pz_Px_S_C1001002_b+ABZ*I_ERI_D2z_Pz_Px_S_C1001002_b;
  Double I_ERI_D2x_D2x_Py_S_C1001002_b = I_ERI_F3x_Px_Py_S_C1001002_b+ABX*I_ERI_D2x_Px_Py_S_C1001002_b;
  Double I_ERI_Dxy_D2x_Py_S_C1001002_b = I_ERI_F2xy_Px_Py_S_C1001002_b+ABX*I_ERI_Dxy_Px_Py_S_C1001002_b;
  Double I_ERI_Dxz_D2x_Py_S_C1001002_b = I_ERI_F2xz_Px_Py_S_C1001002_b+ABX*I_ERI_Dxz_Px_Py_S_C1001002_b;
  Double I_ERI_D2y_D2x_Py_S_C1001002_b = I_ERI_Fx2y_Px_Py_S_C1001002_b+ABX*I_ERI_D2y_Px_Py_S_C1001002_b;
  Double I_ERI_Dyz_D2x_Py_S_C1001002_b = I_ERI_Fxyz_Px_Py_S_C1001002_b+ABX*I_ERI_Dyz_Px_Py_S_C1001002_b;
  Double I_ERI_D2z_D2x_Py_S_C1001002_b = I_ERI_Fx2z_Px_Py_S_C1001002_b+ABX*I_ERI_D2z_Px_Py_S_C1001002_b;
  Double I_ERI_D2x_Dxy_Py_S_C1001002_b = I_ERI_F2xy_Px_Py_S_C1001002_b+ABY*I_ERI_D2x_Px_Py_S_C1001002_b;
  Double I_ERI_Dxy_Dxy_Py_S_C1001002_b = I_ERI_Fx2y_Px_Py_S_C1001002_b+ABY*I_ERI_Dxy_Px_Py_S_C1001002_b;
  Double I_ERI_Dxz_Dxy_Py_S_C1001002_b = I_ERI_Fxyz_Px_Py_S_C1001002_b+ABY*I_ERI_Dxz_Px_Py_S_C1001002_b;
  Double I_ERI_D2y_Dxy_Py_S_C1001002_b = I_ERI_F3y_Px_Py_S_C1001002_b+ABY*I_ERI_D2y_Px_Py_S_C1001002_b;
  Double I_ERI_Dyz_Dxy_Py_S_C1001002_b = I_ERI_F2yz_Px_Py_S_C1001002_b+ABY*I_ERI_Dyz_Px_Py_S_C1001002_b;
  Double I_ERI_D2z_Dxy_Py_S_C1001002_b = I_ERI_Fy2z_Px_Py_S_C1001002_b+ABY*I_ERI_D2z_Px_Py_S_C1001002_b;
  Double I_ERI_D2x_Dxz_Py_S_C1001002_b = I_ERI_F2xz_Px_Py_S_C1001002_b+ABZ*I_ERI_D2x_Px_Py_S_C1001002_b;
  Double I_ERI_Dxy_Dxz_Py_S_C1001002_b = I_ERI_Fxyz_Px_Py_S_C1001002_b+ABZ*I_ERI_Dxy_Px_Py_S_C1001002_b;
  Double I_ERI_Dxz_Dxz_Py_S_C1001002_b = I_ERI_Fx2z_Px_Py_S_C1001002_b+ABZ*I_ERI_Dxz_Px_Py_S_C1001002_b;
  Double I_ERI_D2y_Dxz_Py_S_C1001002_b = I_ERI_F2yz_Px_Py_S_C1001002_b+ABZ*I_ERI_D2y_Px_Py_S_C1001002_b;
  Double I_ERI_Dyz_Dxz_Py_S_C1001002_b = I_ERI_Fy2z_Px_Py_S_C1001002_b+ABZ*I_ERI_Dyz_Px_Py_S_C1001002_b;
  Double I_ERI_D2z_Dxz_Py_S_C1001002_b = I_ERI_F3z_Px_Py_S_C1001002_b+ABZ*I_ERI_D2z_Px_Py_S_C1001002_b;
  Double I_ERI_D2x_D2y_Py_S_C1001002_b = I_ERI_F2xy_Py_Py_S_C1001002_b+ABY*I_ERI_D2x_Py_Py_S_C1001002_b;
  Double I_ERI_Dxy_D2y_Py_S_C1001002_b = I_ERI_Fx2y_Py_Py_S_C1001002_b+ABY*I_ERI_Dxy_Py_Py_S_C1001002_b;
  Double I_ERI_Dxz_D2y_Py_S_C1001002_b = I_ERI_Fxyz_Py_Py_S_C1001002_b+ABY*I_ERI_Dxz_Py_Py_S_C1001002_b;
  Double I_ERI_D2y_D2y_Py_S_C1001002_b = I_ERI_F3y_Py_Py_S_C1001002_b+ABY*I_ERI_D2y_Py_Py_S_C1001002_b;
  Double I_ERI_Dyz_D2y_Py_S_C1001002_b = I_ERI_F2yz_Py_Py_S_C1001002_b+ABY*I_ERI_Dyz_Py_Py_S_C1001002_b;
  Double I_ERI_D2z_D2y_Py_S_C1001002_b = I_ERI_Fy2z_Py_Py_S_C1001002_b+ABY*I_ERI_D2z_Py_Py_S_C1001002_b;
  Double I_ERI_D2x_Dyz_Py_S_C1001002_b = I_ERI_F2xz_Py_Py_S_C1001002_b+ABZ*I_ERI_D2x_Py_Py_S_C1001002_b;
  Double I_ERI_Dxy_Dyz_Py_S_C1001002_b = I_ERI_Fxyz_Py_Py_S_C1001002_b+ABZ*I_ERI_Dxy_Py_Py_S_C1001002_b;
  Double I_ERI_Dxz_Dyz_Py_S_C1001002_b = I_ERI_Fx2z_Py_Py_S_C1001002_b+ABZ*I_ERI_Dxz_Py_Py_S_C1001002_b;
  Double I_ERI_D2y_Dyz_Py_S_C1001002_b = I_ERI_F2yz_Py_Py_S_C1001002_b+ABZ*I_ERI_D2y_Py_Py_S_C1001002_b;
  Double I_ERI_Dyz_Dyz_Py_S_C1001002_b = I_ERI_Fy2z_Py_Py_S_C1001002_b+ABZ*I_ERI_Dyz_Py_Py_S_C1001002_b;
  Double I_ERI_D2z_Dyz_Py_S_C1001002_b = I_ERI_F3z_Py_Py_S_C1001002_b+ABZ*I_ERI_D2z_Py_Py_S_C1001002_b;
  Double I_ERI_D2x_D2z_Py_S_C1001002_b = I_ERI_F2xz_Pz_Py_S_C1001002_b+ABZ*I_ERI_D2x_Pz_Py_S_C1001002_b;
  Double I_ERI_Dxy_D2z_Py_S_C1001002_b = I_ERI_Fxyz_Pz_Py_S_C1001002_b+ABZ*I_ERI_Dxy_Pz_Py_S_C1001002_b;
  Double I_ERI_Dxz_D2z_Py_S_C1001002_b = I_ERI_Fx2z_Pz_Py_S_C1001002_b+ABZ*I_ERI_Dxz_Pz_Py_S_C1001002_b;
  Double I_ERI_D2y_D2z_Py_S_C1001002_b = I_ERI_F2yz_Pz_Py_S_C1001002_b+ABZ*I_ERI_D2y_Pz_Py_S_C1001002_b;
  Double I_ERI_Dyz_D2z_Py_S_C1001002_b = I_ERI_Fy2z_Pz_Py_S_C1001002_b+ABZ*I_ERI_Dyz_Pz_Py_S_C1001002_b;
  Double I_ERI_D2z_D2z_Py_S_C1001002_b = I_ERI_F3z_Pz_Py_S_C1001002_b+ABZ*I_ERI_D2z_Pz_Py_S_C1001002_b;
  Double I_ERI_D2x_D2x_Pz_S_C1001002_b = I_ERI_F3x_Px_Pz_S_C1001002_b+ABX*I_ERI_D2x_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxy_D2x_Pz_S_C1001002_b = I_ERI_F2xy_Px_Pz_S_C1001002_b+ABX*I_ERI_Dxy_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxz_D2x_Pz_S_C1001002_b = I_ERI_F2xz_Px_Pz_S_C1001002_b+ABX*I_ERI_Dxz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2y_D2x_Pz_S_C1001002_b = I_ERI_Fx2y_Px_Pz_S_C1001002_b+ABX*I_ERI_D2y_Px_Pz_S_C1001002_b;
  Double I_ERI_Dyz_D2x_Pz_S_C1001002_b = I_ERI_Fxyz_Px_Pz_S_C1001002_b+ABX*I_ERI_Dyz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2z_D2x_Pz_S_C1001002_b = I_ERI_Fx2z_Px_Pz_S_C1001002_b+ABX*I_ERI_D2z_Px_Pz_S_C1001002_b;
  Double I_ERI_D2x_Dxy_Pz_S_C1001002_b = I_ERI_F2xy_Px_Pz_S_C1001002_b+ABY*I_ERI_D2x_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Dxy_Pz_S_C1001002_b = I_ERI_Fx2y_Px_Pz_S_C1001002_b+ABY*I_ERI_Dxy_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Dxy_Pz_S_C1001002_b = I_ERI_Fxyz_Px_Pz_S_C1001002_b+ABY*I_ERI_Dxz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2y_Dxy_Pz_S_C1001002_b = I_ERI_F3y_Px_Pz_S_C1001002_b+ABY*I_ERI_D2y_Px_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Dxy_Pz_S_C1001002_b = I_ERI_F2yz_Px_Pz_S_C1001002_b+ABY*I_ERI_Dyz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2z_Dxy_Pz_S_C1001002_b = I_ERI_Fy2z_Px_Pz_S_C1001002_b+ABY*I_ERI_D2z_Px_Pz_S_C1001002_b;
  Double I_ERI_D2x_Dxz_Pz_S_C1001002_b = I_ERI_F2xz_Px_Pz_S_C1001002_b+ABZ*I_ERI_D2x_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Dxz_Pz_S_C1001002_b = I_ERI_Fxyz_Px_Pz_S_C1001002_b+ABZ*I_ERI_Dxy_Px_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Dxz_Pz_S_C1001002_b = I_ERI_Fx2z_Px_Pz_S_C1001002_b+ABZ*I_ERI_Dxz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2y_Dxz_Pz_S_C1001002_b = I_ERI_F2yz_Px_Pz_S_C1001002_b+ABZ*I_ERI_D2y_Px_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Dxz_Pz_S_C1001002_b = I_ERI_Fy2z_Px_Pz_S_C1001002_b+ABZ*I_ERI_Dyz_Px_Pz_S_C1001002_b;
  Double I_ERI_D2z_Dxz_Pz_S_C1001002_b = I_ERI_F3z_Px_Pz_S_C1001002_b+ABZ*I_ERI_D2z_Px_Pz_S_C1001002_b;
  Double I_ERI_D2x_D2y_Pz_S_C1001002_b = I_ERI_F2xy_Py_Pz_S_C1001002_b+ABY*I_ERI_D2x_Py_Pz_S_C1001002_b;
  Double I_ERI_Dxy_D2y_Pz_S_C1001002_b = I_ERI_Fx2y_Py_Pz_S_C1001002_b+ABY*I_ERI_Dxy_Py_Pz_S_C1001002_b;
  Double I_ERI_Dxz_D2y_Pz_S_C1001002_b = I_ERI_Fxyz_Py_Pz_S_C1001002_b+ABY*I_ERI_Dxz_Py_Pz_S_C1001002_b;
  Double I_ERI_D2y_D2y_Pz_S_C1001002_b = I_ERI_F3y_Py_Pz_S_C1001002_b+ABY*I_ERI_D2y_Py_Pz_S_C1001002_b;
  Double I_ERI_Dyz_D2y_Pz_S_C1001002_b = I_ERI_F2yz_Py_Pz_S_C1001002_b+ABY*I_ERI_Dyz_Py_Pz_S_C1001002_b;
  Double I_ERI_D2z_D2y_Pz_S_C1001002_b = I_ERI_Fy2z_Py_Pz_S_C1001002_b+ABY*I_ERI_D2z_Py_Pz_S_C1001002_b;
  Double I_ERI_D2x_Dyz_Pz_S_C1001002_b = I_ERI_F2xz_Py_Pz_S_C1001002_b+ABZ*I_ERI_D2x_Py_Pz_S_C1001002_b;
  Double I_ERI_Dxy_Dyz_Pz_S_C1001002_b = I_ERI_Fxyz_Py_Pz_S_C1001002_b+ABZ*I_ERI_Dxy_Py_Pz_S_C1001002_b;
  Double I_ERI_Dxz_Dyz_Pz_S_C1001002_b = I_ERI_Fx2z_Py_Pz_S_C1001002_b+ABZ*I_ERI_Dxz_Py_Pz_S_C1001002_b;
  Double I_ERI_D2y_Dyz_Pz_S_C1001002_b = I_ERI_F2yz_Py_Pz_S_C1001002_b+ABZ*I_ERI_D2y_Py_Pz_S_C1001002_b;
  Double I_ERI_Dyz_Dyz_Pz_S_C1001002_b = I_ERI_Fy2z_Py_Pz_S_C1001002_b+ABZ*I_ERI_Dyz_Py_Pz_S_C1001002_b;
  Double I_ERI_D2z_Dyz_Pz_S_C1001002_b = I_ERI_F3z_Py_Pz_S_C1001002_b+ABZ*I_ERI_D2z_Py_Pz_S_C1001002_b;
  Double I_ERI_D2x_D2z_Pz_S_C1001002_b = I_ERI_F2xz_Pz_Pz_S_C1001002_b+ABZ*I_ERI_D2x_Pz_Pz_S_C1001002_b;
  Double I_ERI_Dxy_D2z_Pz_S_C1001002_b = I_ERI_Fxyz_Pz_Pz_S_C1001002_b+ABZ*I_ERI_Dxy_Pz_Pz_S_C1001002_b;
  Double I_ERI_Dxz_D2z_Pz_S_C1001002_b = I_ERI_Fx2z_Pz_Pz_S_C1001002_b+ABZ*I_ERI_Dxz_Pz_Pz_S_C1001002_b;
  Double I_ERI_D2y_D2z_Pz_S_C1001002_b = I_ERI_F2yz_Pz_Pz_S_C1001002_b+ABZ*I_ERI_D2y_Pz_Pz_S_C1001002_b;
  Double I_ERI_Dyz_D2z_Pz_S_C1001002_b = I_ERI_Fy2z_Pz_Pz_S_C1001002_b+ABZ*I_ERI_Dyz_Pz_Pz_S_C1001002_b;
  Double I_ERI_D2z_D2z_Pz_S_C1001002_b = I_ERI_F3z_Pz_Pz_S_C1001002_b+ABZ*I_ERI_D2z_Pz_Pz_S_C1001002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1002_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1002_c
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1002_c
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_C1002_c = I_ERI_F3x_S_Px_S_C1002_c+ABX*I_ERI_D2x_S_Px_S_C1002_c;
  Double I_ERI_Dxy_Px_Px_S_C1002_c = I_ERI_F2xy_S_Px_S_C1002_c+ABX*I_ERI_Dxy_S_Px_S_C1002_c;
  Double I_ERI_Dxz_Px_Px_S_C1002_c = I_ERI_F2xz_S_Px_S_C1002_c+ABX*I_ERI_Dxz_S_Px_S_C1002_c;
  Double I_ERI_D2y_Px_Px_S_C1002_c = I_ERI_Fx2y_S_Px_S_C1002_c+ABX*I_ERI_D2y_S_Px_S_C1002_c;
  Double I_ERI_Dyz_Px_Px_S_C1002_c = I_ERI_Fxyz_S_Px_S_C1002_c+ABX*I_ERI_Dyz_S_Px_S_C1002_c;
  Double I_ERI_D2z_Px_Px_S_C1002_c = I_ERI_Fx2z_S_Px_S_C1002_c+ABX*I_ERI_D2z_S_Px_S_C1002_c;
  Double I_ERI_D2x_Py_Px_S_C1002_c = I_ERI_F2xy_S_Px_S_C1002_c+ABY*I_ERI_D2x_S_Px_S_C1002_c;
  Double I_ERI_Dxy_Py_Px_S_C1002_c = I_ERI_Fx2y_S_Px_S_C1002_c+ABY*I_ERI_Dxy_S_Px_S_C1002_c;
  Double I_ERI_Dxz_Py_Px_S_C1002_c = I_ERI_Fxyz_S_Px_S_C1002_c+ABY*I_ERI_Dxz_S_Px_S_C1002_c;
  Double I_ERI_D2y_Py_Px_S_C1002_c = I_ERI_F3y_S_Px_S_C1002_c+ABY*I_ERI_D2y_S_Px_S_C1002_c;
  Double I_ERI_Dyz_Py_Px_S_C1002_c = I_ERI_F2yz_S_Px_S_C1002_c+ABY*I_ERI_Dyz_S_Px_S_C1002_c;
  Double I_ERI_D2z_Py_Px_S_C1002_c = I_ERI_Fy2z_S_Px_S_C1002_c+ABY*I_ERI_D2z_S_Px_S_C1002_c;
  Double I_ERI_D2x_Pz_Px_S_C1002_c = I_ERI_F2xz_S_Px_S_C1002_c+ABZ*I_ERI_D2x_S_Px_S_C1002_c;
  Double I_ERI_Dxy_Pz_Px_S_C1002_c = I_ERI_Fxyz_S_Px_S_C1002_c+ABZ*I_ERI_Dxy_S_Px_S_C1002_c;
  Double I_ERI_Dxz_Pz_Px_S_C1002_c = I_ERI_Fx2z_S_Px_S_C1002_c+ABZ*I_ERI_Dxz_S_Px_S_C1002_c;
  Double I_ERI_D2y_Pz_Px_S_C1002_c = I_ERI_F2yz_S_Px_S_C1002_c+ABZ*I_ERI_D2y_S_Px_S_C1002_c;
  Double I_ERI_Dyz_Pz_Px_S_C1002_c = I_ERI_Fy2z_S_Px_S_C1002_c+ABZ*I_ERI_Dyz_S_Px_S_C1002_c;
  Double I_ERI_D2z_Pz_Px_S_C1002_c = I_ERI_F3z_S_Px_S_C1002_c+ABZ*I_ERI_D2z_S_Px_S_C1002_c;
  Double I_ERI_D2x_Px_Py_S_C1002_c = I_ERI_F3x_S_Py_S_C1002_c+ABX*I_ERI_D2x_S_Py_S_C1002_c;
  Double I_ERI_Dxy_Px_Py_S_C1002_c = I_ERI_F2xy_S_Py_S_C1002_c+ABX*I_ERI_Dxy_S_Py_S_C1002_c;
  Double I_ERI_Dxz_Px_Py_S_C1002_c = I_ERI_F2xz_S_Py_S_C1002_c+ABX*I_ERI_Dxz_S_Py_S_C1002_c;
  Double I_ERI_D2y_Px_Py_S_C1002_c = I_ERI_Fx2y_S_Py_S_C1002_c+ABX*I_ERI_D2y_S_Py_S_C1002_c;
  Double I_ERI_Dyz_Px_Py_S_C1002_c = I_ERI_Fxyz_S_Py_S_C1002_c+ABX*I_ERI_Dyz_S_Py_S_C1002_c;
  Double I_ERI_D2z_Px_Py_S_C1002_c = I_ERI_Fx2z_S_Py_S_C1002_c+ABX*I_ERI_D2z_S_Py_S_C1002_c;
  Double I_ERI_D2x_Py_Py_S_C1002_c = I_ERI_F2xy_S_Py_S_C1002_c+ABY*I_ERI_D2x_S_Py_S_C1002_c;
  Double I_ERI_Dxy_Py_Py_S_C1002_c = I_ERI_Fx2y_S_Py_S_C1002_c+ABY*I_ERI_Dxy_S_Py_S_C1002_c;
  Double I_ERI_Dxz_Py_Py_S_C1002_c = I_ERI_Fxyz_S_Py_S_C1002_c+ABY*I_ERI_Dxz_S_Py_S_C1002_c;
  Double I_ERI_D2y_Py_Py_S_C1002_c = I_ERI_F3y_S_Py_S_C1002_c+ABY*I_ERI_D2y_S_Py_S_C1002_c;
  Double I_ERI_Dyz_Py_Py_S_C1002_c = I_ERI_F2yz_S_Py_S_C1002_c+ABY*I_ERI_Dyz_S_Py_S_C1002_c;
  Double I_ERI_D2z_Py_Py_S_C1002_c = I_ERI_Fy2z_S_Py_S_C1002_c+ABY*I_ERI_D2z_S_Py_S_C1002_c;
  Double I_ERI_D2x_Pz_Py_S_C1002_c = I_ERI_F2xz_S_Py_S_C1002_c+ABZ*I_ERI_D2x_S_Py_S_C1002_c;
  Double I_ERI_Dxy_Pz_Py_S_C1002_c = I_ERI_Fxyz_S_Py_S_C1002_c+ABZ*I_ERI_Dxy_S_Py_S_C1002_c;
  Double I_ERI_Dxz_Pz_Py_S_C1002_c = I_ERI_Fx2z_S_Py_S_C1002_c+ABZ*I_ERI_Dxz_S_Py_S_C1002_c;
  Double I_ERI_D2y_Pz_Py_S_C1002_c = I_ERI_F2yz_S_Py_S_C1002_c+ABZ*I_ERI_D2y_S_Py_S_C1002_c;
  Double I_ERI_Dyz_Pz_Py_S_C1002_c = I_ERI_Fy2z_S_Py_S_C1002_c+ABZ*I_ERI_Dyz_S_Py_S_C1002_c;
  Double I_ERI_D2z_Pz_Py_S_C1002_c = I_ERI_F3z_S_Py_S_C1002_c+ABZ*I_ERI_D2z_S_Py_S_C1002_c;
  Double I_ERI_D2x_Px_Pz_S_C1002_c = I_ERI_F3x_S_Pz_S_C1002_c+ABX*I_ERI_D2x_S_Pz_S_C1002_c;
  Double I_ERI_Dxy_Px_Pz_S_C1002_c = I_ERI_F2xy_S_Pz_S_C1002_c+ABX*I_ERI_Dxy_S_Pz_S_C1002_c;
  Double I_ERI_Dxz_Px_Pz_S_C1002_c = I_ERI_F2xz_S_Pz_S_C1002_c+ABX*I_ERI_Dxz_S_Pz_S_C1002_c;
  Double I_ERI_D2y_Px_Pz_S_C1002_c = I_ERI_Fx2y_S_Pz_S_C1002_c+ABX*I_ERI_D2y_S_Pz_S_C1002_c;
  Double I_ERI_Dyz_Px_Pz_S_C1002_c = I_ERI_Fxyz_S_Pz_S_C1002_c+ABX*I_ERI_Dyz_S_Pz_S_C1002_c;
  Double I_ERI_D2z_Px_Pz_S_C1002_c = I_ERI_Fx2z_S_Pz_S_C1002_c+ABX*I_ERI_D2z_S_Pz_S_C1002_c;
  Double I_ERI_D2x_Py_Pz_S_C1002_c = I_ERI_F2xy_S_Pz_S_C1002_c+ABY*I_ERI_D2x_S_Pz_S_C1002_c;
  Double I_ERI_Dxy_Py_Pz_S_C1002_c = I_ERI_Fx2y_S_Pz_S_C1002_c+ABY*I_ERI_Dxy_S_Pz_S_C1002_c;
  Double I_ERI_Dxz_Py_Pz_S_C1002_c = I_ERI_Fxyz_S_Pz_S_C1002_c+ABY*I_ERI_Dxz_S_Pz_S_C1002_c;
  Double I_ERI_D2y_Py_Pz_S_C1002_c = I_ERI_F3y_S_Pz_S_C1002_c+ABY*I_ERI_D2y_S_Pz_S_C1002_c;
  Double I_ERI_Dyz_Py_Pz_S_C1002_c = I_ERI_F2yz_S_Pz_S_C1002_c+ABY*I_ERI_Dyz_S_Pz_S_C1002_c;
  Double I_ERI_D2z_Py_Pz_S_C1002_c = I_ERI_Fy2z_S_Pz_S_C1002_c+ABY*I_ERI_D2z_S_Pz_S_C1002_c;
  Double I_ERI_D2x_Pz_Pz_S_C1002_c = I_ERI_F2xz_S_Pz_S_C1002_c+ABZ*I_ERI_D2x_S_Pz_S_C1002_c;
  Double I_ERI_Dxy_Pz_Pz_S_C1002_c = I_ERI_Fxyz_S_Pz_S_C1002_c+ABZ*I_ERI_Dxy_S_Pz_S_C1002_c;
  Double I_ERI_Dxz_Pz_Pz_S_C1002_c = I_ERI_Fx2z_S_Pz_S_C1002_c+ABZ*I_ERI_Dxz_S_Pz_S_C1002_c;
  Double I_ERI_D2y_Pz_Pz_S_C1002_c = I_ERI_F2yz_S_Pz_S_C1002_c+ABZ*I_ERI_D2y_S_Pz_S_C1002_c;
  Double I_ERI_Dyz_Pz_Pz_S_C1002_c = I_ERI_Fy2z_S_Pz_S_C1002_c+ABZ*I_ERI_Dyz_S_Pz_S_C1002_c;
  Double I_ERI_D2z_Pz_Pz_S_C1002_c = I_ERI_F3z_S_Pz_S_C1002_c+ABZ*I_ERI_D2z_S_Pz_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S_C1001002_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1001002_c
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1001002_c
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S_C1001002_c = I_ERI_F3x_S_D2x_S_C1001002_c+ABX*I_ERI_D2x_S_D2x_S_C1001002_c;
  Double I_ERI_Dxy_Px_D2x_S_C1001002_c = I_ERI_F2xy_S_D2x_S_C1001002_c+ABX*I_ERI_Dxy_S_D2x_S_C1001002_c;
  Double I_ERI_Dxz_Px_D2x_S_C1001002_c = I_ERI_F2xz_S_D2x_S_C1001002_c+ABX*I_ERI_Dxz_S_D2x_S_C1001002_c;
  Double I_ERI_D2y_Px_D2x_S_C1001002_c = I_ERI_Fx2y_S_D2x_S_C1001002_c+ABX*I_ERI_D2y_S_D2x_S_C1001002_c;
  Double I_ERI_Dyz_Px_D2x_S_C1001002_c = I_ERI_Fxyz_S_D2x_S_C1001002_c+ABX*I_ERI_Dyz_S_D2x_S_C1001002_c;
  Double I_ERI_D2z_Px_D2x_S_C1001002_c = I_ERI_Fx2z_S_D2x_S_C1001002_c+ABX*I_ERI_D2z_S_D2x_S_C1001002_c;
  Double I_ERI_D2x_Py_D2x_S_C1001002_c = I_ERI_F2xy_S_D2x_S_C1001002_c+ABY*I_ERI_D2x_S_D2x_S_C1001002_c;
  Double I_ERI_Dxy_Py_D2x_S_C1001002_c = I_ERI_Fx2y_S_D2x_S_C1001002_c+ABY*I_ERI_Dxy_S_D2x_S_C1001002_c;
  Double I_ERI_Dxz_Py_D2x_S_C1001002_c = I_ERI_Fxyz_S_D2x_S_C1001002_c+ABY*I_ERI_Dxz_S_D2x_S_C1001002_c;
  Double I_ERI_D2y_Py_D2x_S_C1001002_c = I_ERI_F3y_S_D2x_S_C1001002_c+ABY*I_ERI_D2y_S_D2x_S_C1001002_c;
  Double I_ERI_Dyz_Py_D2x_S_C1001002_c = I_ERI_F2yz_S_D2x_S_C1001002_c+ABY*I_ERI_Dyz_S_D2x_S_C1001002_c;
  Double I_ERI_D2z_Py_D2x_S_C1001002_c = I_ERI_Fy2z_S_D2x_S_C1001002_c+ABY*I_ERI_D2z_S_D2x_S_C1001002_c;
  Double I_ERI_D2x_Pz_D2x_S_C1001002_c = I_ERI_F2xz_S_D2x_S_C1001002_c+ABZ*I_ERI_D2x_S_D2x_S_C1001002_c;
  Double I_ERI_Dxy_Pz_D2x_S_C1001002_c = I_ERI_Fxyz_S_D2x_S_C1001002_c+ABZ*I_ERI_Dxy_S_D2x_S_C1001002_c;
  Double I_ERI_Dxz_Pz_D2x_S_C1001002_c = I_ERI_Fx2z_S_D2x_S_C1001002_c+ABZ*I_ERI_Dxz_S_D2x_S_C1001002_c;
  Double I_ERI_D2y_Pz_D2x_S_C1001002_c = I_ERI_F2yz_S_D2x_S_C1001002_c+ABZ*I_ERI_D2y_S_D2x_S_C1001002_c;
  Double I_ERI_Dyz_Pz_D2x_S_C1001002_c = I_ERI_Fy2z_S_D2x_S_C1001002_c+ABZ*I_ERI_Dyz_S_D2x_S_C1001002_c;
  Double I_ERI_D2z_Pz_D2x_S_C1001002_c = I_ERI_F3z_S_D2x_S_C1001002_c+ABZ*I_ERI_D2z_S_D2x_S_C1001002_c;
  Double I_ERI_D2x_Px_Dxy_S_C1001002_c = I_ERI_F3x_S_Dxy_S_C1001002_c+ABX*I_ERI_D2x_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxy_Px_Dxy_S_C1001002_c = I_ERI_F2xy_S_Dxy_S_C1001002_c+ABX*I_ERI_Dxy_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxz_Px_Dxy_S_C1001002_c = I_ERI_F2xz_S_Dxy_S_C1001002_c+ABX*I_ERI_Dxz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2y_Px_Dxy_S_C1001002_c = I_ERI_Fx2y_S_Dxy_S_C1001002_c+ABX*I_ERI_D2y_S_Dxy_S_C1001002_c;
  Double I_ERI_Dyz_Px_Dxy_S_C1001002_c = I_ERI_Fxyz_S_Dxy_S_C1001002_c+ABX*I_ERI_Dyz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2z_Px_Dxy_S_C1001002_c = I_ERI_Fx2z_S_Dxy_S_C1001002_c+ABX*I_ERI_D2z_S_Dxy_S_C1001002_c;
  Double I_ERI_D2x_Py_Dxy_S_C1001002_c = I_ERI_F2xy_S_Dxy_S_C1001002_c+ABY*I_ERI_D2x_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxy_Py_Dxy_S_C1001002_c = I_ERI_Fx2y_S_Dxy_S_C1001002_c+ABY*I_ERI_Dxy_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxz_Py_Dxy_S_C1001002_c = I_ERI_Fxyz_S_Dxy_S_C1001002_c+ABY*I_ERI_Dxz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2y_Py_Dxy_S_C1001002_c = I_ERI_F3y_S_Dxy_S_C1001002_c+ABY*I_ERI_D2y_S_Dxy_S_C1001002_c;
  Double I_ERI_Dyz_Py_Dxy_S_C1001002_c = I_ERI_F2yz_S_Dxy_S_C1001002_c+ABY*I_ERI_Dyz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2z_Py_Dxy_S_C1001002_c = I_ERI_Fy2z_S_Dxy_S_C1001002_c+ABY*I_ERI_D2z_S_Dxy_S_C1001002_c;
  Double I_ERI_D2x_Pz_Dxy_S_C1001002_c = I_ERI_F2xz_S_Dxy_S_C1001002_c+ABZ*I_ERI_D2x_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxy_Pz_Dxy_S_C1001002_c = I_ERI_Fxyz_S_Dxy_S_C1001002_c+ABZ*I_ERI_Dxy_S_Dxy_S_C1001002_c;
  Double I_ERI_Dxz_Pz_Dxy_S_C1001002_c = I_ERI_Fx2z_S_Dxy_S_C1001002_c+ABZ*I_ERI_Dxz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2y_Pz_Dxy_S_C1001002_c = I_ERI_F2yz_S_Dxy_S_C1001002_c+ABZ*I_ERI_D2y_S_Dxy_S_C1001002_c;
  Double I_ERI_Dyz_Pz_Dxy_S_C1001002_c = I_ERI_Fy2z_S_Dxy_S_C1001002_c+ABZ*I_ERI_Dyz_S_Dxy_S_C1001002_c;
  Double I_ERI_D2z_Pz_Dxy_S_C1001002_c = I_ERI_F3z_S_Dxy_S_C1001002_c+ABZ*I_ERI_D2z_S_Dxy_S_C1001002_c;
  Double I_ERI_D2x_Px_Dxz_S_C1001002_c = I_ERI_F3x_S_Dxz_S_C1001002_c+ABX*I_ERI_D2x_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxy_Px_Dxz_S_C1001002_c = I_ERI_F2xy_S_Dxz_S_C1001002_c+ABX*I_ERI_Dxy_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxz_Px_Dxz_S_C1001002_c = I_ERI_F2xz_S_Dxz_S_C1001002_c+ABX*I_ERI_Dxz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2y_Px_Dxz_S_C1001002_c = I_ERI_Fx2y_S_Dxz_S_C1001002_c+ABX*I_ERI_D2y_S_Dxz_S_C1001002_c;
  Double I_ERI_Dyz_Px_Dxz_S_C1001002_c = I_ERI_Fxyz_S_Dxz_S_C1001002_c+ABX*I_ERI_Dyz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2z_Px_Dxz_S_C1001002_c = I_ERI_Fx2z_S_Dxz_S_C1001002_c+ABX*I_ERI_D2z_S_Dxz_S_C1001002_c;
  Double I_ERI_D2x_Py_Dxz_S_C1001002_c = I_ERI_F2xy_S_Dxz_S_C1001002_c+ABY*I_ERI_D2x_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxy_Py_Dxz_S_C1001002_c = I_ERI_Fx2y_S_Dxz_S_C1001002_c+ABY*I_ERI_Dxy_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxz_Py_Dxz_S_C1001002_c = I_ERI_Fxyz_S_Dxz_S_C1001002_c+ABY*I_ERI_Dxz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2y_Py_Dxz_S_C1001002_c = I_ERI_F3y_S_Dxz_S_C1001002_c+ABY*I_ERI_D2y_S_Dxz_S_C1001002_c;
  Double I_ERI_Dyz_Py_Dxz_S_C1001002_c = I_ERI_F2yz_S_Dxz_S_C1001002_c+ABY*I_ERI_Dyz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2z_Py_Dxz_S_C1001002_c = I_ERI_Fy2z_S_Dxz_S_C1001002_c+ABY*I_ERI_D2z_S_Dxz_S_C1001002_c;
  Double I_ERI_D2x_Pz_Dxz_S_C1001002_c = I_ERI_F2xz_S_Dxz_S_C1001002_c+ABZ*I_ERI_D2x_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxy_Pz_Dxz_S_C1001002_c = I_ERI_Fxyz_S_Dxz_S_C1001002_c+ABZ*I_ERI_Dxy_S_Dxz_S_C1001002_c;
  Double I_ERI_Dxz_Pz_Dxz_S_C1001002_c = I_ERI_Fx2z_S_Dxz_S_C1001002_c+ABZ*I_ERI_Dxz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2y_Pz_Dxz_S_C1001002_c = I_ERI_F2yz_S_Dxz_S_C1001002_c+ABZ*I_ERI_D2y_S_Dxz_S_C1001002_c;
  Double I_ERI_Dyz_Pz_Dxz_S_C1001002_c = I_ERI_Fy2z_S_Dxz_S_C1001002_c+ABZ*I_ERI_Dyz_S_Dxz_S_C1001002_c;
  Double I_ERI_D2z_Pz_Dxz_S_C1001002_c = I_ERI_F3z_S_Dxz_S_C1001002_c+ABZ*I_ERI_D2z_S_Dxz_S_C1001002_c;
  Double I_ERI_D2x_Px_D2y_S_C1001002_c = I_ERI_F3x_S_D2y_S_C1001002_c+ABX*I_ERI_D2x_S_D2y_S_C1001002_c;
  Double I_ERI_Dxy_Px_D2y_S_C1001002_c = I_ERI_F2xy_S_D2y_S_C1001002_c+ABX*I_ERI_Dxy_S_D2y_S_C1001002_c;
  Double I_ERI_Dxz_Px_D2y_S_C1001002_c = I_ERI_F2xz_S_D2y_S_C1001002_c+ABX*I_ERI_Dxz_S_D2y_S_C1001002_c;
  Double I_ERI_D2y_Px_D2y_S_C1001002_c = I_ERI_Fx2y_S_D2y_S_C1001002_c+ABX*I_ERI_D2y_S_D2y_S_C1001002_c;
  Double I_ERI_Dyz_Px_D2y_S_C1001002_c = I_ERI_Fxyz_S_D2y_S_C1001002_c+ABX*I_ERI_Dyz_S_D2y_S_C1001002_c;
  Double I_ERI_D2z_Px_D2y_S_C1001002_c = I_ERI_Fx2z_S_D2y_S_C1001002_c+ABX*I_ERI_D2z_S_D2y_S_C1001002_c;
  Double I_ERI_D2x_Py_D2y_S_C1001002_c = I_ERI_F2xy_S_D2y_S_C1001002_c+ABY*I_ERI_D2x_S_D2y_S_C1001002_c;
  Double I_ERI_Dxy_Py_D2y_S_C1001002_c = I_ERI_Fx2y_S_D2y_S_C1001002_c+ABY*I_ERI_Dxy_S_D2y_S_C1001002_c;
  Double I_ERI_Dxz_Py_D2y_S_C1001002_c = I_ERI_Fxyz_S_D2y_S_C1001002_c+ABY*I_ERI_Dxz_S_D2y_S_C1001002_c;
  Double I_ERI_D2y_Py_D2y_S_C1001002_c = I_ERI_F3y_S_D2y_S_C1001002_c+ABY*I_ERI_D2y_S_D2y_S_C1001002_c;
  Double I_ERI_Dyz_Py_D2y_S_C1001002_c = I_ERI_F2yz_S_D2y_S_C1001002_c+ABY*I_ERI_Dyz_S_D2y_S_C1001002_c;
  Double I_ERI_D2z_Py_D2y_S_C1001002_c = I_ERI_Fy2z_S_D2y_S_C1001002_c+ABY*I_ERI_D2z_S_D2y_S_C1001002_c;
  Double I_ERI_D2x_Pz_D2y_S_C1001002_c = I_ERI_F2xz_S_D2y_S_C1001002_c+ABZ*I_ERI_D2x_S_D2y_S_C1001002_c;
  Double I_ERI_Dxy_Pz_D2y_S_C1001002_c = I_ERI_Fxyz_S_D2y_S_C1001002_c+ABZ*I_ERI_Dxy_S_D2y_S_C1001002_c;
  Double I_ERI_Dxz_Pz_D2y_S_C1001002_c = I_ERI_Fx2z_S_D2y_S_C1001002_c+ABZ*I_ERI_Dxz_S_D2y_S_C1001002_c;
  Double I_ERI_D2y_Pz_D2y_S_C1001002_c = I_ERI_F2yz_S_D2y_S_C1001002_c+ABZ*I_ERI_D2y_S_D2y_S_C1001002_c;
  Double I_ERI_Dyz_Pz_D2y_S_C1001002_c = I_ERI_Fy2z_S_D2y_S_C1001002_c+ABZ*I_ERI_Dyz_S_D2y_S_C1001002_c;
  Double I_ERI_D2z_Pz_D2y_S_C1001002_c = I_ERI_F3z_S_D2y_S_C1001002_c+ABZ*I_ERI_D2z_S_D2y_S_C1001002_c;
  Double I_ERI_D2x_Px_Dyz_S_C1001002_c = I_ERI_F3x_S_Dyz_S_C1001002_c+ABX*I_ERI_D2x_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxy_Px_Dyz_S_C1001002_c = I_ERI_F2xy_S_Dyz_S_C1001002_c+ABX*I_ERI_Dxy_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxz_Px_Dyz_S_C1001002_c = I_ERI_F2xz_S_Dyz_S_C1001002_c+ABX*I_ERI_Dxz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2y_Px_Dyz_S_C1001002_c = I_ERI_Fx2y_S_Dyz_S_C1001002_c+ABX*I_ERI_D2y_S_Dyz_S_C1001002_c;
  Double I_ERI_Dyz_Px_Dyz_S_C1001002_c = I_ERI_Fxyz_S_Dyz_S_C1001002_c+ABX*I_ERI_Dyz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2z_Px_Dyz_S_C1001002_c = I_ERI_Fx2z_S_Dyz_S_C1001002_c+ABX*I_ERI_D2z_S_Dyz_S_C1001002_c;
  Double I_ERI_D2x_Py_Dyz_S_C1001002_c = I_ERI_F2xy_S_Dyz_S_C1001002_c+ABY*I_ERI_D2x_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxy_Py_Dyz_S_C1001002_c = I_ERI_Fx2y_S_Dyz_S_C1001002_c+ABY*I_ERI_Dxy_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxz_Py_Dyz_S_C1001002_c = I_ERI_Fxyz_S_Dyz_S_C1001002_c+ABY*I_ERI_Dxz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2y_Py_Dyz_S_C1001002_c = I_ERI_F3y_S_Dyz_S_C1001002_c+ABY*I_ERI_D2y_S_Dyz_S_C1001002_c;
  Double I_ERI_Dyz_Py_Dyz_S_C1001002_c = I_ERI_F2yz_S_Dyz_S_C1001002_c+ABY*I_ERI_Dyz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2z_Py_Dyz_S_C1001002_c = I_ERI_Fy2z_S_Dyz_S_C1001002_c+ABY*I_ERI_D2z_S_Dyz_S_C1001002_c;
  Double I_ERI_D2x_Pz_Dyz_S_C1001002_c = I_ERI_F2xz_S_Dyz_S_C1001002_c+ABZ*I_ERI_D2x_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxy_Pz_Dyz_S_C1001002_c = I_ERI_Fxyz_S_Dyz_S_C1001002_c+ABZ*I_ERI_Dxy_S_Dyz_S_C1001002_c;
  Double I_ERI_Dxz_Pz_Dyz_S_C1001002_c = I_ERI_Fx2z_S_Dyz_S_C1001002_c+ABZ*I_ERI_Dxz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2y_Pz_Dyz_S_C1001002_c = I_ERI_F2yz_S_Dyz_S_C1001002_c+ABZ*I_ERI_D2y_S_Dyz_S_C1001002_c;
  Double I_ERI_Dyz_Pz_Dyz_S_C1001002_c = I_ERI_Fy2z_S_Dyz_S_C1001002_c+ABZ*I_ERI_Dyz_S_Dyz_S_C1001002_c;
  Double I_ERI_D2z_Pz_Dyz_S_C1001002_c = I_ERI_F3z_S_Dyz_S_C1001002_c+ABZ*I_ERI_D2z_S_Dyz_S_C1001002_c;
  Double I_ERI_D2x_Px_D2z_S_C1001002_c = I_ERI_F3x_S_D2z_S_C1001002_c+ABX*I_ERI_D2x_S_D2z_S_C1001002_c;
  Double I_ERI_Dxy_Px_D2z_S_C1001002_c = I_ERI_F2xy_S_D2z_S_C1001002_c+ABX*I_ERI_Dxy_S_D2z_S_C1001002_c;
  Double I_ERI_Dxz_Px_D2z_S_C1001002_c = I_ERI_F2xz_S_D2z_S_C1001002_c+ABX*I_ERI_Dxz_S_D2z_S_C1001002_c;
  Double I_ERI_D2y_Px_D2z_S_C1001002_c = I_ERI_Fx2y_S_D2z_S_C1001002_c+ABX*I_ERI_D2y_S_D2z_S_C1001002_c;
  Double I_ERI_Dyz_Px_D2z_S_C1001002_c = I_ERI_Fxyz_S_D2z_S_C1001002_c+ABX*I_ERI_Dyz_S_D2z_S_C1001002_c;
  Double I_ERI_D2z_Px_D2z_S_C1001002_c = I_ERI_Fx2z_S_D2z_S_C1001002_c+ABX*I_ERI_D2z_S_D2z_S_C1001002_c;
  Double I_ERI_D2x_Py_D2z_S_C1001002_c = I_ERI_F2xy_S_D2z_S_C1001002_c+ABY*I_ERI_D2x_S_D2z_S_C1001002_c;
  Double I_ERI_Dxy_Py_D2z_S_C1001002_c = I_ERI_Fx2y_S_D2z_S_C1001002_c+ABY*I_ERI_Dxy_S_D2z_S_C1001002_c;
  Double I_ERI_Dxz_Py_D2z_S_C1001002_c = I_ERI_Fxyz_S_D2z_S_C1001002_c+ABY*I_ERI_Dxz_S_D2z_S_C1001002_c;
  Double I_ERI_D2y_Py_D2z_S_C1001002_c = I_ERI_F3y_S_D2z_S_C1001002_c+ABY*I_ERI_D2y_S_D2z_S_C1001002_c;
  Double I_ERI_Dyz_Py_D2z_S_C1001002_c = I_ERI_F2yz_S_D2z_S_C1001002_c+ABY*I_ERI_Dyz_S_D2z_S_C1001002_c;
  Double I_ERI_D2z_Py_D2z_S_C1001002_c = I_ERI_Fy2z_S_D2z_S_C1001002_c+ABY*I_ERI_D2z_S_D2z_S_C1001002_c;
  Double I_ERI_D2x_Pz_D2z_S_C1001002_c = I_ERI_F2xz_S_D2z_S_C1001002_c+ABZ*I_ERI_D2x_S_D2z_S_C1001002_c;
  Double I_ERI_Dxy_Pz_D2z_S_C1001002_c = I_ERI_Fxyz_S_D2z_S_C1001002_c+ABZ*I_ERI_Dxy_S_D2z_S_C1001002_c;
  Double I_ERI_Dxz_Pz_D2z_S_C1001002_c = I_ERI_Fx2z_S_D2z_S_C1001002_c+ABZ*I_ERI_Dxz_S_D2z_S_C1001002_c;
  Double I_ERI_D2y_Pz_D2z_S_C1001002_c = I_ERI_F2yz_S_D2z_S_C1001002_c+ABZ*I_ERI_D2y_S_D2z_S_C1001002_c;
  Double I_ERI_Dyz_Pz_D2z_S_C1001002_c = I_ERI_Fy2z_S_D2z_S_C1001002_c+ABZ*I_ERI_Dyz_S_D2z_S_C1001002_c;
  Double I_ERI_D2z_Pz_D2z_S_C1001002_c = I_ERI_F3z_S_D2z_S_C1001002_c+ABZ*I_ERI_D2z_S_D2z_S_C1001002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C2
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_F3x_S_S_S_C2_a-2*I_ERI_Px_S_S_S_C2;
  abcd[1] = 2.0E0*I_ERI_F2xy_S_S_S_C2_a-1*I_ERI_Py_S_S_S_C2;
  abcd[2] = 2.0E0*I_ERI_F2xz_S_S_S_C2_a-1*I_ERI_Pz_S_S_S_C2;
  abcd[3] = 2.0E0*I_ERI_Fx2y_S_S_S_C2_a;
  abcd[4] = 2.0E0*I_ERI_Fxyz_S_S_S_C2_a;
  abcd[5] = 2.0E0*I_ERI_Fx2z_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1002
   ************************************************************/
  abcd[6] = 2.0E0*I_ERI_F3x_Px_S_S_C1002_a-2*I_ERI_Px_Px_S_S_C1002;
  abcd[7] = 2.0E0*I_ERI_F2xy_Px_S_S_C1002_a-1*I_ERI_Py_Px_S_S_C1002;
  abcd[8] = 2.0E0*I_ERI_F2xz_Px_S_S_C1002_a-1*I_ERI_Pz_Px_S_S_C1002;
  abcd[9] = 2.0E0*I_ERI_Fx2y_Px_S_S_C1002_a;
  abcd[10] = 2.0E0*I_ERI_Fxyz_Px_S_S_C1002_a;
  abcd[11] = 2.0E0*I_ERI_Fx2z_Px_S_S_C1002_a;
  abcd[12] = 2.0E0*I_ERI_F3x_Py_S_S_C1002_a-2*I_ERI_Px_Py_S_S_C1002;
  abcd[13] = 2.0E0*I_ERI_F2xy_Py_S_S_C1002_a-1*I_ERI_Py_Py_S_S_C1002;
  abcd[14] = 2.0E0*I_ERI_F2xz_Py_S_S_C1002_a-1*I_ERI_Pz_Py_S_S_C1002;
  abcd[15] = 2.0E0*I_ERI_Fx2y_Py_S_S_C1002_a;
  abcd[16] = 2.0E0*I_ERI_Fxyz_Py_S_S_C1002_a;
  abcd[17] = 2.0E0*I_ERI_Fx2z_Py_S_S_C1002_a;
  abcd[18] = 2.0E0*I_ERI_F3x_Pz_S_S_C1002_a-2*I_ERI_Px_Pz_S_S_C1002;
  abcd[19] = 2.0E0*I_ERI_F2xy_Pz_S_S_C1002_a-1*I_ERI_Py_Pz_S_S_C1002;
  abcd[20] = 2.0E0*I_ERI_F2xz_Pz_S_S_C1002_a-1*I_ERI_Pz_Pz_S_S_C1002;
  abcd[21] = 2.0E0*I_ERI_Fx2y_Pz_S_S_C1002_a;
  abcd[22] = 2.0E0*I_ERI_Fxyz_Pz_S_S_C1002_a;
  abcd[23] = 2.0E0*I_ERI_Fx2z_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000002_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000002
   ************************************************************/
  abcd[24] = 2.0E0*I_ERI_F3x_S_Px_S_C1000002_a-2*I_ERI_Px_S_Px_S_C1000002;
  abcd[25] = 2.0E0*I_ERI_F2xy_S_Px_S_C1000002_a-1*I_ERI_Py_S_Px_S_C1000002;
  abcd[26] = 2.0E0*I_ERI_F2xz_S_Px_S_C1000002_a-1*I_ERI_Pz_S_Px_S_C1000002;
  abcd[27] = 2.0E0*I_ERI_Fx2y_S_Px_S_C1000002_a;
  abcd[28] = 2.0E0*I_ERI_Fxyz_S_Px_S_C1000002_a;
  abcd[29] = 2.0E0*I_ERI_Fx2z_S_Px_S_C1000002_a;
  abcd[48] = 2.0E0*I_ERI_F3x_S_Py_S_C1000002_a-2*I_ERI_Px_S_Py_S_C1000002;
  abcd[49] = 2.0E0*I_ERI_F2xy_S_Py_S_C1000002_a-1*I_ERI_Py_S_Py_S_C1000002;
  abcd[50] = 2.0E0*I_ERI_F2xz_S_Py_S_C1000002_a-1*I_ERI_Pz_S_Py_S_C1000002;
  abcd[51] = 2.0E0*I_ERI_Fx2y_S_Py_S_C1000002_a;
  abcd[52] = 2.0E0*I_ERI_Fxyz_S_Py_S_C1000002_a;
  abcd[53] = 2.0E0*I_ERI_Fx2z_S_Py_S_C1000002_a;
  abcd[72] = 2.0E0*I_ERI_F3x_S_Pz_S_C1000002_a-2*I_ERI_Px_S_Pz_S_C1000002;
  abcd[73] = 2.0E0*I_ERI_F2xy_S_Pz_S_C1000002_a-1*I_ERI_Py_S_Pz_S_C1000002;
  abcd[74] = 2.0E0*I_ERI_F2xz_S_Pz_S_C1000002_a-1*I_ERI_Pz_S_Pz_S_C1000002;
  abcd[75] = 2.0E0*I_ERI_Fx2y_S_Pz_S_C1000002_a;
  abcd[76] = 2.0E0*I_ERI_Fxyz_S_Pz_S_C1000002_a;
  abcd[77] = 2.0E0*I_ERI_Fx2z_S_Pz_S_C1000002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001002_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001002
   ************************************************************/
  abcd[30] = 2.0E0*I_ERI_F3x_Px_Px_S_C1001002_a-2*I_ERI_Px_Px_Px_S_C1001002;
  abcd[31] = 2.0E0*I_ERI_F2xy_Px_Px_S_C1001002_a-1*I_ERI_Py_Px_Px_S_C1001002;
  abcd[32] = 2.0E0*I_ERI_F2xz_Px_Px_S_C1001002_a-1*I_ERI_Pz_Px_Px_S_C1001002;
  abcd[33] = 2.0E0*I_ERI_Fx2y_Px_Px_S_C1001002_a;
  abcd[34] = 2.0E0*I_ERI_Fxyz_Px_Px_S_C1001002_a;
  abcd[35] = 2.0E0*I_ERI_Fx2z_Px_Px_S_C1001002_a;
  abcd[36] = 2.0E0*I_ERI_F3x_Py_Px_S_C1001002_a-2*I_ERI_Px_Py_Px_S_C1001002;
  abcd[37] = 2.0E0*I_ERI_F2xy_Py_Px_S_C1001002_a-1*I_ERI_Py_Py_Px_S_C1001002;
  abcd[38] = 2.0E0*I_ERI_F2xz_Py_Px_S_C1001002_a-1*I_ERI_Pz_Py_Px_S_C1001002;
  abcd[39] = 2.0E0*I_ERI_Fx2y_Py_Px_S_C1001002_a;
  abcd[40] = 2.0E0*I_ERI_Fxyz_Py_Px_S_C1001002_a;
  abcd[41] = 2.0E0*I_ERI_Fx2z_Py_Px_S_C1001002_a;
  abcd[42] = 2.0E0*I_ERI_F3x_Pz_Px_S_C1001002_a-2*I_ERI_Px_Pz_Px_S_C1001002;
  abcd[43] = 2.0E0*I_ERI_F2xy_Pz_Px_S_C1001002_a-1*I_ERI_Py_Pz_Px_S_C1001002;
  abcd[44] = 2.0E0*I_ERI_F2xz_Pz_Px_S_C1001002_a-1*I_ERI_Pz_Pz_Px_S_C1001002;
  abcd[45] = 2.0E0*I_ERI_Fx2y_Pz_Px_S_C1001002_a;
  abcd[46] = 2.0E0*I_ERI_Fxyz_Pz_Px_S_C1001002_a;
  abcd[47] = 2.0E0*I_ERI_Fx2z_Pz_Px_S_C1001002_a;
  abcd[54] = 2.0E0*I_ERI_F3x_Px_Py_S_C1001002_a-2*I_ERI_Px_Px_Py_S_C1001002;
  abcd[55] = 2.0E0*I_ERI_F2xy_Px_Py_S_C1001002_a-1*I_ERI_Py_Px_Py_S_C1001002;
  abcd[56] = 2.0E0*I_ERI_F2xz_Px_Py_S_C1001002_a-1*I_ERI_Pz_Px_Py_S_C1001002;
  abcd[57] = 2.0E0*I_ERI_Fx2y_Px_Py_S_C1001002_a;
  abcd[58] = 2.0E0*I_ERI_Fxyz_Px_Py_S_C1001002_a;
  abcd[59] = 2.0E0*I_ERI_Fx2z_Px_Py_S_C1001002_a;
  abcd[60] = 2.0E0*I_ERI_F3x_Py_Py_S_C1001002_a-2*I_ERI_Px_Py_Py_S_C1001002;
  abcd[61] = 2.0E0*I_ERI_F2xy_Py_Py_S_C1001002_a-1*I_ERI_Py_Py_Py_S_C1001002;
  abcd[62] = 2.0E0*I_ERI_F2xz_Py_Py_S_C1001002_a-1*I_ERI_Pz_Py_Py_S_C1001002;
  abcd[63] = 2.0E0*I_ERI_Fx2y_Py_Py_S_C1001002_a;
  abcd[64] = 2.0E0*I_ERI_Fxyz_Py_Py_S_C1001002_a;
  abcd[65] = 2.0E0*I_ERI_Fx2z_Py_Py_S_C1001002_a;
  abcd[66] = 2.0E0*I_ERI_F3x_Pz_Py_S_C1001002_a-2*I_ERI_Px_Pz_Py_S_C1001002;
  abcd[67] = 2.0E0*I_ERI_F2xy_Pz_Py_S_C1001002_a-1*I_ERI_Py_Pz_Py_S_C1001002;
  abcd[68] = 2.0E0*I_ERI_F2xz_Pz_Py_S_C1001002_a-1*I_ERI_Pz_Pz_Py_S_C1001002;
  abcd[69] = 2.0E0*I_ERI_Fx2y_Pz_Py_S_C1001002_a;
  abcd[70] = 2.0E0*I_ERI_Fxyz_Pz_Py_S_C1001002_a;
  abcd[71] = 2.0E0*I_ERI_Fx2z_Pz_Py_S_C1001002_a;
  abcd[78] = 2.0E0*I_ERI_F3x_Px_Pz_S_C1001002_a-2*I_ERI_Px_Px_Pz_S_C1001002;
  abcd[79] = 2.0E0*I_ERI_F2xy_Px_Pz_S_C1001002_a-1*I_ERI_Py_Px_Pz_S_C1001002;
  abcd[80] = 2.0E0*I_ERI_F2xz_Px_Pz_S_C1001002_a-1*I_ERI_Pz_Px_Pz_S_C1001002;
  abcd[81] = 2.0E0*I_ERI_Fx2y_Px_Pz_S_C1001002_a;
  abcd[82] = 2.0E0*I_ERI_Fxyz_Px_Pz_S_C1001002_a;
  abcd[83] = 2.0E0*I_ERI_Fx2z_Px_Pz_S_C1001002_a;
  abcd[84] = 2.0E0*I_ERI_F3x_Py_Pz_S_C1001002_a-2*I_ERI_Px_Py_Pz_S_C1001002;
  abcd[85] = 2.0E0*I_ERI_F2xy_Py_Pz_S_C1001002_a-1*I_ERI_Py_Py_Pz_S_C1001002;
  abcd[86] = 2.0E0*I_ERI_F2xz_Py_Pz_S_C1001002_a-1*I_ERI_Pz_Py_Pz_S_C1001002;
  abcd[87] = 2.0E0*I_ERI_Fx2y_Py_Pz_S_C1001002_a;
  abcd[88] = 2.0E0*I_ERI_Fxyz_Py_Pz_S_C1001002_a;
  abcd[89] = 2.0E0*I_ERI_Fx2z_Py_Pz_S_C1001002_a;
  abcd[90] = 2.0E0*I_ERI_F3x_Pz_Pz_S_C1001002_a-2*I_ERI_Px_Pz_Pz_S_C1001002;
  abcd[91] = 2.0E0*I_ERI_F2xy_Pz_Pz_S_C1001002_a-1*I_ERI_Py_Pz_Pz_S_C1001002;
  abcd[92] = 2.0E0*I_ERI_F2xz_Pz_Pz_S_C1001002_a-1*I_ERI_Pz_Pz_Pz_S_C1001002;
  abcd[93] = 2.0E0*I_ERI_Fx2y_Pz_Pz_S_C1001002_a;
  abcd[94] = 2.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001002_a;
  abcd[95] = 2.0E0*I_ERI_Fx2z_Pz_Pz_S_C1001002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C2
   ************************************************************/
  abcd[96] = 2.0E0*I_ERI_F2xy_S_S_S_C2_a;
  abcd[97] = 2.0E0*I_ERI_Fx2y_S_S_S_C2_a-1*I_ERI_Px_S_S_S_C2;
  abcd[98] = 2.0E0*I_ERI_Fxyz_S_S_S_C2_a;
  abcd[99] = 2.0E0*I_ERI_F3y_S_S_S_C2_a-2*I_ERI_Py_S_S_S_C2;
  abcd[100] = 2.0E0*I_ERI_F2yz_S_S_S_C2_a-1*I_ERI_Pz_S_S_S_C2;
  abcd[101] = 2.0E0*I_ERI_Fy2z_S_S_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1002
   ************************************************************/
  abcd[102] = 2.0E0*I_ERI_F2xy_Px_S_S_C1002_a;
  abcd[103] = 2.0E0*I_ERI_Fx2y_Px_S_S_C1002_a-1*I_ERI_Px_Px_S_S_C1002;
  abcd[104] = 2.0E0*I_ERI_Fxyz_Px_S_S_C1002_a;
  abcd[105] = 2.0E0*I_ERI_F3y_Px_S_S_C1002_a-2*I_ERI_Py_Px_S_S_C1002;
  abcd[106] = 2.0E0*I_ERI_F2yz_Px_S_S_C1002_a-1*I_ERI_Pz_Px_S_S_C1002;
  abcd[107] = 2.0E0*I_ERI_Fy2z_Px_S_S_C1002_a;
  abcd[108] = 2.0E0*I_ERI_F2xy_Py_S_S_C1002_a;
  abcd[109] = 2.0E0*I_ERI_Fx2y_Py_S_S_C1002_a-1*I_ERI_Px_Py_S_S_C1002;
  abcd[110] = 2.0E0*I_ERI_Fxyz_Py_S_S_C1002_a;
  abcd[111] = 2.0E0*I_ERI_F3y_Py_S_S_C1002_a-2*I_ERI_Py_Py_S_S_C1002;
  abcd[112] = 2.0E0*I_ERI_F2yz_Py_S_S_C1002_a-1*I_ERI_Pz_Py_S_S_C1002;
  abcd[113] = 2.0E0*I_ERI_Fy2z_Py_S_S_C1002_a;
  abcd[114] = 2.0E0*I_ERI_F2xy_Pz_S_S_C1002_a;
  abcd[115] = 2.0E0*I_ERI_Fx2y_Pz_S_S_C1002_a-1*I_ERI_Px_Pz_S_S_C1002;
  abcd[116] = 2.0E0*I_ERI_Fxyz_Pz_S_S_C1002_a;
  abcd[117] = 2.0E0*I_ERI_F3y_Pz_S_S_C1002_a-2*I_ERI_Py_Pz_S_S_C1002;
  abcd[118] = 2.0E0*I_ERI_F2yz_Pz_S_S_C1002_a-1*I_ERI_Pz_Pz_S_S_C1002;
  abcd[119] = 2.0E0*I_ERI_Fy2z_Pz_S_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000002_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000002
   ************************************************************/
  abcd[120] = 2.0E0*I_ERI_F2xy_S_Px_S_C1000002_a;
  abcd[121] = 2.0E0*I_ERI_Fx2y_S_Px_S_C1000002_a-1*I_ERI_Px_S_Px_S_C1000002;
  abcd[122] = 2.0E0*I_ERI_Fxyz_S_Px_S_C1000002_a;
  abcd[123] = 2.0E0*I_ERI_F3y_S_Px_S_C1000002_a-2*I_ERI_Py_S_Px_S_C1000002;
  abcd[124] = 2.0E0*I_ERI_F2yz_S_Px_S_C1000002_a-1*I_ERI_Pz_S_Px_S_C1000002;
  abcd[125] = 2.0E0*I_ERI_Fy2z_S_Px_S_C1000002_a;
  abcd[144] = 2.0E0*I_ERI_F2xy_S_Py_S_C1000002_a;
  abcd[145] = 2.0E0*I_ERI_Fx2y_S_Py_S_C1000002_a-1*I_ERI_Px_S_Py_S_C1000002;
  abcd[146] = 2.0E0*I_ERI_Fxyz_S_Py_S_C1000002_a;
  abcd[147] = 2.0E0*I_ERI_F3y_S_Py_S_C1000002_a-2*I_ERI_Py_S_Py_S_C1000002;
  abcd[148] = 2.0E0*I_ERI_F2yz_S_Py_S_C1000002_a-1*I_ERI_Pz_S_Py_S_C1000002;
  abcd[149] = 2.0E0*I_ERI_Fy2z_S_Py_S_C1000002_a;
  abcd[168] = 2.0E0*I_ERI_F2xy_S_Pz_S_C1000002_a;
  abcd[169] = 2.0E0*I_ERI_Fx2y_S_Pz_S_C1000002_a-1*I_ERI_Px_S_Pz_S_C1000002;
  abcd[170] = 2.0E0*I_ERI_Fxyz_S_Pz_S_C1000002_a;
  abcd[171] = 2.0E0*I_ERI_F3y_S_Pz_S_C1000002_a-2*I_ERI_Py_S_Pz_S_C1000002;
  abcd[172] = 2.0E0*I_ERI_F2yz_S_Pz_S_C1000002_a-1*I_ERI_Pz_S_Pz_S_C1000002;
  abcd[173] = 2.0E0*I_ERI_Fy2z_S_Pz_S_C1000002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001002_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001002
   ************************************************************/
  abcd[126] = 2.0E0*I_ERI_F2xy_Px_Px_S_C1001002_a;
  abcd[127] = 2.0E0*I_ERI_Fx2y_Px_Px_S_C1001002_a-1*I_ERI_Px_Px_Px_S_C1001002;
  abcd[128] = 2.0E0*I_ERI_Fxyz_Px_Px_S_C1001002_a;
  abcd[129] = 2.0E0*I_ERI_F3y_Px_Px_S_C1001002_a-2*I_ERI_Py_Px_Px_S_C1001002;
  abcd[130] = 2.0E0*I_ERI_F2yz_Px_Px_S_C1001002_a-1*I_ERI_Pz_Px_Px_S_C1001002;
  abcd[131] = 2.0E0*I_ERI_Fy2z_Px_Px_S_C1001002_a;
  abcd[132] = 2.0E0*I_ERI_F2xy_Py_Px_S_C1001002_a;
  abcd[133] = 2.0E0*I_ERI_Fx2y_Py_Px_S_C1001002_a-1*I_ERI_Px_Py_Px_S_C1001002;
  abcd[134] = 2.0E0*I_ERI_Fxyz_Py_Px_S_C1001002_a;
  abcd[135] = 2.0E0*I_ERI_F3y_Py_Px_S_C1001002_a-2*I_ERI_Py_Py_Px_S_C1001002;
  abcd[136] = 2.0E0*I_ERI_F2yz_Py_Px_S_C1001002_a-1*I_ERI_Pz_Py_Px_S_C1001002;
  abcd[137] = 2.0E0*I_ERI_Fy2z_Py_Px_S_C1001002_a;
  abcd[138] = 2.0E0*I_ERI_F2xy_Pz_Px_S_C1001002_a;
  abcd[139] = 2.0E0*I_ERI_Fx2y_Pz_Px_S_C1001002_a-1*I_ERI_Px_Pz_Px_S_C1001002;
  abcd[140] = 2.0E0*I_ERI_Fxyz_Pz_Px_S_C1001002_a;
  abcd[141] = 2.0E0*I_ERI_F3y_Pz_Px_S_C1001002_a-2*I_ERI_Py_Pz_Px_S_C1001002;
  abcd[142] = 2.0E0*I_ERI_F2yz_Pz_Px_S_C1001002_a-1*I_ERI_Pz_Pz_Px_S_C1001002;
  abcd[143] = 2.0E0*I_ERI_Fy2z_Pz_Px_S_C1001002_a;
  abcd[150] = 2.0E0*I_ERI_F2xy_Px_Py_S_C1001002_a;
  abcd[151] = 2.0E0*I_ERI_Fx2y_Px_Py_S_C1001002_a-1*I_ERI_Px_Px_Py_S_C1001002;
  abcd[152] = 2.0E0*I_ERI_Fxyz_Px_Py_S_C1001002_a;
  abcd[153] = 2.0E0*I_ERI_F3y_Px_Py_S_C1001002_a-2*I_ERI_Py_Px_Py_S_C1001002;
  abcd[154] = 2.0E0*I_ERI_F2yz_Px_Py_S_C1001002_a-1*I_ERI_Pz_Px_Py_S_C1001002;
  abcd[155] = 2.0E0*I_ERI_Fy2z_Px_Py_S_C1001002_a;
  abcd[156] = 2.0E0*I_ERI_F2xy_Py_Py_S_C1001002_a;
  abcd[157] = 2.0E0*I_ERI_Fx2y_Py_Py_S_C1001002_a-1*I_ERI_Px_Py_Py_S_C1001002;
  abcd[158] = 2.0E0*I_ERI_Fxyz_Py_Py_S_C1001002_a;
  abcd[159] = 2.0E0*I_ERI_F3y_Py_Py_S_C1001002_a-2*I_ERI_Py_Py_Py_S_C1001002;
  abcd[160] = 2.0E0*I_ERI_F2yz_Py_Py_S_C1001002_a-1*I_ERI_Pz_Py_Py_S_C1001002;
  abcd[161] = 2.0E0*I_ERI_Fy2z_Py_Py_S_C1001002_a;
  abcd[162] = 2.0E0*I_ERI_F2xy_Pz_Py_S_C1001002_a;
  abcd[163] = 2.0E0*I_ERI_Fx2y_Pz_Py_S_C1001002_a-1*I_ERI_Px_Pz_Py_S_C1001002;
  abcd[164] = 2.0E0*I_ERI_Fxyz_Pz_Py_S_C1001002_a;
  abcd[165] = 2.0E0*I_ERI_F3y_Pz_Py_S_C1001002_a-2*I_ERI_Py_Pz_Py_S_C1001002;
  abcd[166] = 2.0E0*I_ERI_F2yz_Pz_Py_S_C1001002_a-1*I_ERI_Pz_Pz_Py_S_C1001002;
  abcd[167] = 2.0E0*I_ERI_Fy2z_Pz_Py_S_C1001002_a;
  abcd[174] = 2.0E0*I_ERI_F2xy_Px_Pz_S_C1001002_a;
  abcd[175] = 2.0E0*I_ERI_Fx2y_Px_Pz_S_C1001002_a-1*I_ERI_Px_Px_Pz_S_C1001002;
  abcd[176] = 2.0E0*I_ERI_Fxyz_Px_Pz_S_C1001002_a;
  abcd[177] = 2.0E0*I_ERI_F3y_Px_Pz_S_C1001002_a-2*I_ERI_Py_Px_Pz_S_C1001002;
  abcd[178] = 2.0E0*I_ERI_F2yz_Px_Pz_S_C1001002_a-1*I_ERI_Pz_Px_Pz_S_C1001002;
  abcd[179] = 2.0E0*I_ERI_Fy2z_Px_Pz_S_C1001002_a;
  abcd[180] = 2.0E0*I_ERI_F2xy_Py_Pz_S_C1001002_a;
  abcd[181] = 2.0E0*I_ERI_Fx2y_Py_Pz_S_C1001002_a-1*I_ERI_Px_Py_Pz_S_C1001002;
  abcd[182] = 2.0E0*I_ERI_Fxyz_Py_Pz_S_C1001002_a;
  abcd[183] = 2.0E0*I_ERI_F3y_Py_Pz_S_C1001002_a-2*I_ERI_Py_Py_Pz_S_C1001002;
  abcd[184] = 2.0E0*I_ERI_F2yz_Py_Pz_S_C1001002_a-1*I_ERI_Pz_Py_Pz_S_C1001002;
  abcd[185] = 2.0E0*I_ERI_Fy2z_Py_Pz_S_C1001002_a;
  abcd[186] = 2.0E0*I_ERI_F2xy_Pz_Pz_S_C1001002_a;
  abcd[187] = 2.0E0*I_ERI_Fx2y_Pz_Pz_S_C1001002_a-1*I_ERI_Px_Pz_Pz_S_C1001002;
  abcd[188] = 2.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001002_a;
  abcd[189] = 2.0E0*I_ERI_F3y_Pz_Pz_S_C1001002_a-2*I_ERI_Py_Pz_Pz_S_C1001002;
  abcd[190] = 2.0E0*I_ERI_F2yz_Pz_Pz_S_C1001002_a-1*I_ERI_Pz_Pz_Pz_S_C1001002;
  abcd[191] = 2.0E0*I_ERI_Fy2z_Pz_Pz_S_C1001002_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_C2_a
   * RHS shell quartet name: SQ_ERI_P_S_S_S_C2
   ************************************************************/
  abcd[192] = 2.0E0*I_ERI_F2xz_S_S_S_C2_a;
  abcd[193] = 2.0E0*I_ERI_Fxyz_S_S_S_C2_a;
  abcd[194] = 2.0E0*I_ERI_Fx2z_S_S_S_C2_a-1*I_ERI_Px_S_S_S_C2;
  abcd[195] = 2.0E0*I_ERI_F2yz_S_S_S_C2_a;
  abcd[196] = 2.0E0*I_ERI_Fy2z_S_S_S_C2_a-1*I_ERI_Py_S_S_S_C2;
  abcd[197] = 2.0E0*I_ERI_F3z_S_S_S_C2_a-2*I_ERI_Pz_S_S_S_C2;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_C1002_a
   * RHS shell quartet name: SQ_ERI_P_P_S_S_C1002
   ************************************************************/
  abcd[198] = 2.0E0*I_ERI_F2xz_Px_S_S_C1002_a;
  abcd[199] = 2.0E0*I_ERI_Fxyz_Px_S_S_C1002_a;
  abcd[200] = 2.0E0*I_ERI_Fx2z_Px_S_S_C1002_a-1*I_ERI_Px_Px_S_S_C1002;
  abcd[201] = 2.0E0*I_ERI_F2yz_Px_S_S_C1002_a;
  abcd[202] = 2.0E0*I_ERI_Fy2z_Px_S_S_C1002_a-1*I_ERI_Py_Px_S_S_C1002;
  abcd[203] = 2.0E0*I_ERI_F3z_Px_S_S_C1002_a-2*I_ERI_Pz_Px_S_S_C1002;
  abcd[204] = 2.0E0*I_ERI_F2xz_Py_S_S_C1002_a;
  abcd[205] = 2.0E0*I_ERI_Fxyz_Py_S_S_C1002_a;
  abcd[206] = 2.0E0*I_ERI_Fx2z_Py_S_S_C1002_a-1*I_ERI_Px_Py_S_S_C1002;
  abcd[207] = 2.0E0*I_ERI_F2yz_Py_S_S_C1002_a;
  abcd[208] = 2.0E0*I_ERI_Fy2z_Py_S_S_C1002_a-1*I_ERI_Py_Py_S_S_C1002;
  abcd[209] = 2.0E0*I_ERI_F3z_Py_S_S_C1002_a-2*I_ERI_Pz_Py_S_S_C1002;
  abcd[210] = 2.0E0*I_ERI_F2xz_Pz_S_S_C1002_a;
  abcd[211] = 2.0E0*I_ERI_Fxyz_Pz_S_S_C1002_a;
  abcd[212] = 2.0E0*I_ERI_Fx2z_Pz_S_S_C1002_a-1*I_ERI_Px_Pz_S_S_C1002;
  abcd[213] = 2.0E0*I_ERI_F2yz_Pz_S_S_C1002_a;
  abcd[214] = 2.0E0*I_ERI_Fy2z_Pz_S_S_C1002_a-1*I_ERI_Py_Pz_S_S_C1002;
  abcd[215] = 2.0E0*I_ERI_F3z_Pz_S_S_C1002_a-2*I_ERI_Pz_Pz_S_S_C1002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_C1000002_a
   * RHS shell quartet name: SQ_ERI_P_S_P_S_C1000002
   ************************************************************/
  abcd[216] = 2.0E0*I_ERI_F2xz_S_Px_S_C1000002_a;
  abcd[217] = 2.0E0*I_ERI_Fxyz_S_Px_S_C1000002_a;
  abcd[218] = 2.0E0*I_ERI_Fx2z_S_Px_S_C1000002_a-1*I_ERI_Px_S_Px_S_C1000002;
  abcd[219] = 2.0E0*I_ERI_F2yz_S_Px_S_C1000002_a;
  abcd[220] = 2.0E0*I_ERI_Fy2z_S_Px_S_C1000002_a-1*I_ERI_Py_S_Px_S_C1000002;
  abcd[221] = 2.0E0*I_ERI_F3z_S_Px_S_C1000002_a-2*I_ERI_Pz_S_Px_S_C1000002;
  abcd[240] = 2.0E0*I_ERI_F2xz_S_Py_S_C1000002_a;
  abcd[241] = 2.0E0*I_ERI_Fxyz_S_Py_S_C1000002_a;
  abcd[242] = 2.0E0*I_ERI_Fx2z_S_Py_S_C1000002_a-1*I_ERI_Px_S_Py_S_C1000002;
  abcd[243] = 2.0E0*I_ERI_F2yz_S_Py_S_C1000002_a;
  abcd[244] = 2.0E0*I_ERI_Fy2z_S_Py_S_C1000002_a-1*I_ERI_Py_S_Py_S_C1000002;
  abcd[245] = 2.0E0*I_ERI_F3z_S_Py_S_C1000002_a-2*I_ERI_Pz_S_Py_S_C1000002;
  abcd[264] = 2.0E0*I_ERI_F2xz_S_Pz_S_C1000002_a;
  abcd[265] = 2.0E0*I_ERI_Fxyz_S_Pz_S_C1000002_a;
  abcd[266] = 2.0E0*I_ERI_Fx2z_S_Pz_S_C1000002_a-1*I_ERI_Px_S_Pz_S_C1000002;
  abcd[267] = 2.0E0*I_ERI_F2yz_S_Pz_S_C1000002_a;
  abcd[268] = 2.0E0*I_ERI_Fy2z_S_Pz_S_C1000002_a-1*I_ERI_Py_S_Pz_S_C1000002;
  abcd[269] = 2.0E0*I_ERI_F3z_S_Pz_S_C1000002_a-2*I_ERI_Pz_S_Pz_S_C1000002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_C1001002_a
   * RHS shell quartet name: SQ_ERI_P_P_P_S_C1001002
   ************************************************************/
  abcd[222] = 2.0E0*I_ERI_F2xz_Px_Px_S_C1001002_a;
  abcd[223] = 2.0E0*I_ERI_Fxyz_Px_Px_S_C1001002_a;
  abcd[224] = 2.0E0*I_ERI_Fx2z_Px_Px_S_C1001002_a-1*I_ERI_Px_Px_Px_S_C1001002;
  abcd[225] = 2.0E0*I_ERI_F2yz_Px_Px_S_C1001002_a;
  abcd[226] = 2.0E0*I_ERI_Fy2z_Px_Px_S_C1001002_a-1*I_ERI_Py_Px_Px_S_C1001002;
  abcd[227] = 2.0E0*I_ERI_F3z_Px_Px_S_C1001002_a-2*I_ERI_Pz_Px_Px_S_C1001002;
  abcd[228] = 2.0E0*I_ERI_F2xz_Py_Px_S_C1001002_a;
  abcd[229] = 2.0E0*I_ERI_Fxyz_Py_Px_S_C1001002_a;
  abcd[230] = 2.0E0*I_ERI_Fx2z_Py_Px_S_C1001002_a-1*I_ERI_Px_Py_Px_S_C1001002;
  abcd[231] = 2.0E0*I_ERI_F2yz_Py_Px_S_C1001002_a;
  abcd[232] = 2.0E0*I_ERI_Fy2z_Py_Px_S_C1001002_a-1*I_ERI_Py_Py_Px_S_C1001002;
  abcd[233] = 2.0E0*I_ERI_F3z_Py_Px_S_C1001002_a-2*I_ERI_Pz_Py_Px_S_C1001002;
  abcd[234] = 2.0E0*I_ERI_F2xz_Pz_Px_S_C1001002_a;
  abcd[235] = 2.0E0*I_ERI_Fxyz_Pz_Px_S_C1001002_a;
  abcd[236] = 2.0E0*I_ERI_Fx2z_Pz_Px_S_C1001002_a-1*I_ERI_Px_Pz_Px_S_C1001002;
  abcd[237] = 2.0E0*I_ERI_F2yz_Pz_Px_S_C1001002_a;
  abcd[238] = 2.0E0*I_ERI_Fy2z_Pz_Px_S_C1001002_a-1*I_ERI_Py_Pz_Px_S_C1001002;
  abcd[239] = 2.0E0*I_ERI_F3z_Pz_Px_S_C1001002_a-2*I_ERI_Pz_Pz_Px_S_C1001002;
  abcd[246] = 2.0E0*I_ERI_F2xz_Px_Py_S_C1001002_a;
  abcd[247] = 2.0E0*I_ERI_Fxyz_Px_Py_S_C1001002_a;
  abcd[248] = 2.0E0*I_ERI_Fx2z_Px_Py_S_C1001002_a-1*I_ERI_Px_Px_Py_S_C1001002;
  abcd[249] = 2.0E0*I_ERI_F2yz_Px_Py_S_C1001002_a;
  abcd[250] = 2.0E0*I_ERI_Fy2z_Px_Py_S_C1001002_a-1*I_ERI_Py_Px_Py_S_C1001002;
  abcd[251] = 2.0E0*I_ERI_F3z_Px_Py_S_C1001002_a-2*I_ERI_Pz_Px_Py_S_C1001002;
  abcd[252] = 2.0E0*I_ERI_F2xz_Py_Py_S_C1001002_a;
  abcd[253] = 2.0E0*I_ERI_Fxyz_Py_Py_S_C1001002_a;
  abcd[254] = 2.0E0*I_ERI_Fx2z_Py_Py_S_C1001002_a-1*I_ERI_Px_Py_Py_S_C1001002;
  abcd[255] = 2.0E0*I_ERI_F2yz_Py_Py_S_C1001002_a;
  abcd[256] = 2.0E0*I_ERI_Fy2z_Py_Py_S_C1001002_a-1*I_ERI_Py_Py_Py_S_C1001002;
  abcd[257] = 2.0E0*I_ERI_F3z_Py_Py_S_C1001002_a-2*I_ERI_Pz_Py_Py_S_C1001002;
  abcd[258] = 2.0E0*I_ERI_F2xz_Pz_Py_S_C1001002_a;
  abcd[259] = 2.0E0*I_ERI_Fxyz_Pz_Py_S_C1001002_a;
  abcd[260] = 2.0E0*I_ERI_Fx2z_Pz_Py_S_C1001002_a-1*I_ERI_Px_Pz_Py_S_C1001002;
  abcd[261] = 2.0E0*I_ERI_F2yz_Pz_Py_S_C1001002_a;
  abcd[262] = 2.0E0*I_ERI_Fy2z_Pz_Py_S_C1001002_a-1*I_ERI_Py_Pz_Py_S_C1001002;
  abcd[263] = 2.0E0*I_ERI_F3z_Pz_Py_S_C1001002_a-2*I_ERI_Pz_Pz_Py_S_C1001002;
  abcd[270] = 2.0E0*I_ERI_F2xz_Px_Pz_S_C1001002_a;
  abcd[271] = 2.0E0*I_ERI_Fxyz_Px_Pz_S_C1001002_a;
  abcd[272] = 2.0E0*I_ERI_Fx2z_Px_Pz_S_C1001002_a-1*I_ERI_Px_Px_Pz_S_C1001002;
  abcd[273] = 2.0E0*I_ERI_F2yz_Px_Pz_S_C1001002_a;
  abcd[274] = 2.0E0*I_ERI_Fy2z_Px_Pz_S_C1001002_a-1*I_ERI_Py_Px_Pz_S_C1001002;
  abcd[275] = 2.0E0*I_ERI_F3z_Px_Pz_S_C1001002_a-2*I_ERI_Pz_Px_Pz_S_C1001002;
  abcd[276] = 2.0E0*I_ERI_F2xz_Py_Pz_S_C1001002_a;
  abcd[277] = 2.0E0*I_ERI_Fxyz_Py_Pz_S_C1001002_a;
  abcd[278] = 2.0E0*I_ERI_Fx2z_Py_Pz_S_C1001002_a-1*I_ERI_Px_Py_Pz_S_C1001002;
  abcd[279] = 2.0E0*I_ERI_F2yz_Py_Pz_S_C1001002_a;
  abcd[280] = 2.0E0*I_ERI_Fy2z_Py_Pz_S_C1001002_a-1*I_ERI_Py_Py_Pz_S_C1001002;
  abcd[281] = 2.0E0*I_ERI_F3z_Py_Pz_S_C1001002_a-2*I_ERI_Pz_Py_Pz_S_C1001002;
  abcd[282] = 2.0E0*I_ERI_F2xz_Pz_Pz_S_C1001002_a;
  abcd[283] = 2.0E0*I_ERI_Fxyz_Pz_Pz_S_C1001002_a;
  abcd[284] = 2.0E0*I_ERI_Fx2z_Pz_Pz_S_C1001002_a-1*I_ERI_Px_Pz_Pz_S_C1001002;
  abcd[285] = 2.0E0*I_ERI_F2yz_Pz_Pz_S_C1001002_a;
  abcd[286] = 2.0E0*I_ERI_Fy2z_Pz_Pz_S_C1001002_a-1*I_ERI_Py_Pz_Pz_S_C1001002;
  abcd[287] = 2.0E0*I_ERI_F3z_Pz_Pz_S_C1001002_a-2*I_ERI_Pz_Pz_Pz_S_C1001002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C2_b
   ************************************************************/
  abcd[288] = 2.0E0*I_ERI_D2x_Px_S_S_C2_b;
  abcd[289] = 2.0E0*I_ERI_Dxy_Px_S_S_C2_b;
  abcd[290] = 2.0E0*I_ERI_Dxz_Px_S_S_C2_b;
  abcd[291] = 2.0E0*I_ERI_D2y_Px_S_S_C2_b;
  abcd[292] = 2.0E0*I_ERI_Dyz_Px_S_S_C2_b;
  abcd[293] = 2.0E0*I_ERI_D2z_Px_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002
   ************************************************************/
  abcd[294] = 2.0E0*I_ERI_D2x_D2x_S_S_C1002_b-1*I_ERI_D2x_S_S_S_C1002;
  abcd[295] = 2.0E0*I_ERI_Dxy_D2x_S_S_C1002_b-1*I_ERI_Dxy_S_S_S_C1002;
  abcd[296] = 2.0E0*I_ERI_Dxz_D2x_S_S_C1002_b-1*I_ERI_Dxz_S_S_S_C1002;
  abcd[297] = 2.0E0*I_ERI_D2y_D2x_S_S_C1002_b-1*I_ERI_D2y_S_S_S_C1002;
  abcd[298] = 2.0E0*I_ERI_Dyz_D2x_S_S_C1002_b-1*I_ERI_Dyz_S_S_S_C1002;
  abcd[299] = 2.0E0*I_ERI_D2z_D2x_S_S_C1002_b-1*I_ERI_D2z_S_S_S_C1002;
  abcd[300] = 2.0E0*I_ERI_D2x_Dxy_S_S_C1002_b;
  abcd[301] = 2.0E0*I_ERI_Dxy_Dxy_S_S_C1002_b;
  abcd[302] = 2.0E0*I_ERI_Dxz_Dxy_S_S_C1002_b;
  abcd[303] = 2.0E0*I_ERI_D2y_Dxy_S_S_C1002_b;
  abcd[304] = 2.0E0*I_ERI_Dyz_Dxy_S_S_C1002_b;
  abcd[305] = 2.0E0*I_ERI_D2z_Dxy_S_S_C1002_b;
  abcd[306] = 2.0E0*I_ERI_D2x_Dxz_S_S_C1002_b;
  abcd[307] = 2.0E0*I_ERI_Dxy_Dxz_S_S_C1002_b;
  abcd[308] = 2.0E0*I_ERI_Dxz_Dxz_S_S_C1002_b;
  abcd[309] = 2.0E0*I_ERI_D2y_Dxz_S_S_C1002_b;
  abcd[310] = 2.0E0*I_ERI_Dyz_Dxz_S_S_C1002_b;
  abcd[311] = 2.0E0*I_ERI_D2z_Dxz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1000002_b
   ************************************************************/
  abcd[312] = 2.0E0*I_ERI_D2x_Px_Px_S_C1000002_b;
  abcd[313] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1000002_b;
  abcd[314] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1000002_b;
  abcd[315] = 2.0E0*I_ERI_D2y_Px_Px_S_C1000002_b;
  abcd[316] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1000002_b;
  abcd[317] = 2.0E0*I_ERI_D2z_Px_Px_S_C1000002_b;
  abcd[336] = 2.0E0*I_ERI_D2x_Px_Py_S_C1000002_b;
  abcd[337] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1000002_b;
  abcd[338] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1000002_b;
  abcd[339] = 2.0E0*I_ERI_D2y_Px_Py_S_C1000002_b;
  abcd[340] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1000002_b;
  abcd[341] = 2.0E0*I_ERI_D2z_Px_Py_S_C1000002_b;
  abcd[360] = 2.0E0*I_ERI_D2x_Px_Pz_S_C1000002_b;
  abcd[361] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1000002_b;
  abcd[362] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1000002_b;
  abcd[363] = 2.0E0*I_ERI_D2y_Px_Pz_S_C1000002_b;
  abcd[364] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1000002_b;
  abcd[365] = 2.0E0*I_ERI_D2z_Px_Pz_S_C1000002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001002
   ************************************************************/
  abcd[318] = 2.0E0*I_ERI_D2x_D2x_Px_S_C1001002_b-1*I_ERI_D2x_S_Px_S_C1001002;
  abcd[319] = 2.0E0*I_ERI_Dxy_D2x_Px_S_C1001002_b-1*I_ERI_Dxy_S_Px_S_C1001002;
  abcd[320] = 2.0E0*I_ERI_Dxz_D2x_Px_S_C1001002_b-1*I_ERI_Dxz_S_Px_S_C1001002;
  abcd[321] = 2.0E0*I_ERI_D2y_D2x_Px_S_C1001002_b-1*I_ERI_D2y_S_Px_S_C1001002;
  abcd[322] = 2.0E0*I_ERI_Dyz_D2x_Px_S_C1001002_b-1*I_ERI_Dyz_S_Px_S_C1001002;
  abcd[323] = 2.0E0*I_ERI_D2z_D2x_Px_S_C1001002_b-1*I_ERI_D2z_S_Px_S_C1001002;
  abcd[324] = 2.0E0*I_ERI_D2x_Dxy_Px_S_C1001002_b;
  abcd[325] = 2.0E0*I_ERI_Dxy_Dxy_Px_S_C1001002_b;
  abcd[326] = 2.0E0*I_ERI_Dxz_Dxy_Px_S_C1001002_b;
  abcd[327] = 2.0E0*I_ERI_D2y_Dxy_Px_S_C1001002_b;
  abcd[328] = 2.0E0*I_ERI_Dyz_Dxy_Px_S_C1001002_b;
  abcd[329] = 2.0E0*I_ERI_D2z_Dxy_Px_S_C1001002_b;
  abcd[330] = 2.0E0*I_ERI_D2x_Dxz_Px_S_C1001002_b;
  abcd[331] = 2.0E0*I_ERI_Dxy_Dxz_Px_S_C1001002_b;
  abcd[332] = 2.0E0*I_ERI_Dxz_Dxz_Px_S_C1001002_b;
  abcd[333] = 2.0E0*I_ERI_D2y_Dxz_Px_S_C1001002_b;
  abcd[334] = 2.0E0*I_ERI_Dyz_Dxz_Px_S_C1001002_b;
  abcd[335] = 2.0E0*I_ERI_D2z_Dxz_Px_S_C1001002_b;
  abcd[342] = 2.0E0*I_ERI_D2x_D2x_Py_S_C1001002_b-1*I_ERI_D2x_S_Py_S_C1001002;
  abcd[343] = 2.0E0*I_ERI_Dxy_D2x_Py_S_C1001002_b-1*I_ERI_Dxy_S_Py_S_C1001002;
  abcd[344] = 2.0E0*I_ERI_Dxz_D2x_Py_S_C1001002_b-1*I_ERI_Dxz_S_Py_S_C1001002;
  abcd[345] = 2.0E0*I_ERI_D2y_D2x_Py_S_C1001002_b-1*I_ERI_D2y_S_Py_S_C1001002;
  abcd[346] = 2.0E0*I_ERI_Dyz_D2x_Py_S_C1001002_b-1*I_ERI_Dyz_S_Py_S_C1001002;
  abcd[347] = 2.0E0*I_ERI_D2z_D2x_Py_S_C1001002_b-1*I_ERI_D2z_S_Py_S_C1001002;
  abcd[348] = 2.0E0*I_ERI_D2x_Dxy_Py_S_C1001002_b;
  abcd[349] = 2.0E0*I_ERI_Dxy_Dxy_Py_S_C1001002_b;
  abcd[350] = 2.0E0*I_ERI_Dxz_Dxy_Py_S_C1001002_b;
  abcd[351] = 2.0E0*I_ERI_D2y_Dxy_Py_S_C1001002_b;
  abcd[352] = 2.0E0*I_ERI_Dyz_Dxy_Py_S_C1001002_b;
  abcd[353] = 2.0E0*I_ERI_D2z_Dxy_Py_S_C1001002_b;
  abcd[354] = 2.0E0*I_ERI_D2x_Dxz_Py_S_C1001002_b;
  abcd[355] = 2.0E0*I_ERI_Dxy_Dxz_Py_S_C1001002_b;
  abcd[356] = 2.0E0*I_ERI_Dxz_Dxz_Py_S_C1001002_b;
  abcd[357] = 2.0E0*I_ERI_D2y_Dxz_Py_S_C1001002_b;
  abcd[358] = 2.0E0*I_ERI_Dyz_Dxz_Py_S_C1001002_b;
  abcd[359] = 2.0E0*I_ERI_D2z_Dxz_Py_S_C1001002_b;
  abcd[366] = 2.0E0*I_ERI_D2x_D2x_Pz_S_C1001002_b-1*I_ERI_D2x_S_Pz_S_C1001002;
  abcd[367] = 2.0E0*I_ERI_Dxy_D2x_Pz_S_C1001002_b-1*I_ERI_Dxy_S_Pz_S_C1001002;
  abcd[368] = 2.0E0*I_ERI_Dxz_D2x_Pz_S_C1001002_b-1*I_ERI_Dxz_S_Pz_S_C1001002;
  abcd[369] = 2.0E0*I_ERI_D2y_D2x_Pz_S_C1001002_b-1*I_ERI_D2y_S_Pz_S_C1001002;
  abcd[370] = 2.0E0*I_ERI_Dyz_D2x_Pz_S_C1001002_b-1*I_ERI_Dyz_S_Pz_S_C1001002;
  abcd[371] = 2.0E0*I_ERI_D2z_D2x_Pz_S_C1001002_b-1*I_ERI_D2z_S_Pz_S_C1001002;
  abcd[372] = 2.0E0*I_ERI_D2x_Dxy_Pz_S_C1001002_b;
  abcd[373] = 2.0E0*I_ERI_Dxy_Dxy_Pz_S_C1001002_b;
  abcd[374] = 2.0E0*I_ERI_Dxz_Dxy_Pz_S_C1001002_b;
  abcd[375] = 2.0E0*I_ERI_D2y_Dxy_Pz_S_C1001002_b;
  abcd[376] = 2.0E0*I_ERI_Dyz_Dxy_Pz_S_C1001002_b;
  abcd[377] = 2.0E0*I_ERI_D2z_Dxy_Pz_S_C1001002_b;
  abcd[378] = 2.0E0*I_ERI_D2x_Dxz_Pz_S_C1001002_b;
  abcd[379] = 2.0E0*I_ERI_Dxy_Dxz_Pz_S_C1001002_b;
  abcd[380] = 2.0E0*I_ERI_Dxz_Dxz_Pz_S_C1001002_b;
  abcd[381] = 2.0E0*I_ERI_D2y_Dxz_Pz_S_C1001002_b;
  abcd[382] = 2.0E0*I_ERI_Dyz_Dxz_Pz_S_C1001002_b;
  abcd[383] = 2.0E0*I_ERI_D2z_Dxz_Pz_S_C1001002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C2_b
   ************************************************************/
  abcd[384] = 2.0E0*I_ERI_D2x_Py_S_S_C2_b;
  abcd[385] = 2.0E0*I_ERI_Dxy_Py_S_S_C2_b;
  abcd[386] = 2.0E0*I_ERI_Dxz_Py_S_S_C2_b;
  abcd[387] = 2.0E0*I_ERI_D2y_Py_S_S_C2_b;
  abcd[388] = 2.0E0*I_ERI_Dyz_Py_S_S_C2_b;
  abcd[389] = 2.0E0*I_ERI_D2z_Py_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002
   ************************************************************/
  abcd[390] = 2.0E0*I_ERI_D2x_Dxy_S_S_C1002_b;
  abcd[391] = 2.0E0*I_ERI_Dxy_Dxy_S_S_C1002_b;
  abcd[392] = 2.0E0*I_ERI_Dxz_Dxy_S_S_C1002_b;
  abcd[393] = 2.0E0*I_ERI_D2y_Dxy_S_S_C1002_b;
  abcd[394] = 2.0E0*I_ERI_Dyz_Dxy_S_S_C1002_b;
  abcd[395] = 2.0E0*I_ERI_D2z_Dxy_S_S_C1002_b;
  abcd[396] = 2.0E0*I_ERI_D2x_D2y_S_S_C1002_b-1*I_ERI_D2x_S_S_S_C1002;
  abcd[397] = 2.0E0*I_ERI_Dxy_D2y_S_S_C1002_b-1*I_ERI_Dxy_S_S_S_C1002;
  abcd[398] = 2.0E0*I_ERI_Dxz_D2y_S_S_C1002_b-1*I_ERI_Dxz_S_S_S_C1002;
  abcd[399] = 2.0E0*I_ERI_D2y_D2y_S_S_C1002_b-1*I_ERI_D2y_S_S_S_C1002;
  abcd[400] = 2.0E0*I_ERI_Dyz_D2y_S_S_C1002_b-1*I_ERI_Dyz_S_S_S_C1002;
  abcd[401] = 2.0E0*I_ERI_D2z_D2y_S_S_C1002_b-1*I_ERI_D2z_S_S_S_C1002;
  abcd[402] = 2.0E0*I_ERI_D2x_Dyz_S_S_C1002_b;
  abcd[403] = 2.0E0*I_ERI_Dxy_Dyz_S_S_C1002_b;
  abcd[404] = 2.0E0*I_ERI_Dxz_Dyz_S_S_C1002_b;
  abcd[405] = 2.0E0*I_ERI_D2y_Dyz_S_S_C1002_b;
  abcd[406] = 2.0E0*I_ERI_Dyz_Dyz_S_S_C1002_b;
  abcd[407] = 2.0E0*I_ERI_D2z_Dyz_S_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1000002_b
   ************************************************************/
  abcd[408] = 2.0E0*I_ERI_D2x_Py_Px_S_C1000002_b;
  abcd[409] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1000002_b;
  abcd[410] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1000002_b;
  abcd[411] = 2.0E0*I_ERI_D2y_Py_Px_S_C1000002_b;
  abcd[412] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1000002_b;
  abcd[413] = 2.0E0*I_ERI_D2z_Py_Px_S_C1000002_b;
  abcd[432] = 2.0E0*I_ERI_D2x_Py_Py_S_C1000002_b;
  abcd[433] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1000002_b;
  abcd[434] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1000002_b;
  abcd[435] = 2.0E0*I_ERI_D2y_Py_Py_S_C1000002_b;
  abcd[436] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1000002_b;
  abcd[437] = 2.0E0*I_ERI_D2z_Py_Py_S_C1000002_b;
  abcd[456] = 2.0E0*I_ERI_D2x_Py_Pz_S_C1000002_b;
  abcd[457] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1000002_b;
  abcd[458] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1000002_b;
  abcd[459] = 2.0E0*I_ERI_D2y_Py_Pz_S_C1000002_b;
  abcd[460] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1000002_b;
  abcd[461] = 2.0E0*I_ERI_D2z_Py_Pz_S_C1000002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001002
   ************************************************************/
  abcd[414] = 2.0E0*I_ERI_D2x_Dxy_Px_S_C1001002_b;
  abcd[415] = 2.0E0*I_ERI_Dxy_Dxy_Px_S_C1001002_b;
  abcd[416] = 2.0E0*I_ERI_Dxz_Dxy_Px_S_C1001002_b;
  abcd[417] = 2.0E0*I_ERI_D2y_Dxy_Px_S_C1001002_b;
  abcd[418] = 2.0E0*I_ERI_Dyz_Dxy_Px_S_C1001002_b;
  abcd[419] = 2.0E0*I_ERI_D2z_Dxy_Px_S_C1001002_b;
  abcd[420] = 2.0E0*I_ERI_D2x_D2y_Px_S_C1001002_b-1*I_ERI_D2x_S_Px_S_C1001002;
  abcd[421] = 2.0E0*I_ERI_Dxy_D2y_Px_S_C1001002_b-1*I_ERI_Dxy_S_Px_S_C1001002;
  abcd[422] = 2.0E0*I_ERI_Dxz_D2y_Px_S_C1001002_b-1*I_ERI_Dxz_S_Px_S_C1001002;
  abcd[423] = 2.0E0*I_ERI_D2y_D2y_Px_S_C1001002_b-1*I_ERI_D2y_S_Px_S_C1001002;
  abcd[424] = 2.0E0*I_ERI_Dyz_D2y_Px_S_C1001002_b-1*I_ERI_Dyz_S_Px_S_C1001002;
  abcd[425] = 2.0E0*I_ERI_D2z_D2y_Px_S_C1001002_b-1*I_ERI_D2z_S_Px_S_C1001002;
  abcd[426] = 2.0E0*I_ERI_D2x_Dyz_Px_S_C1001002_b;
  abcd[427] = 2.0E0*I_ERI_Dxy_Dyz_Px_S_C1001002_b;
  abcd[428] = 2.0E0*I_ERI_Dxz_Dyz_Px_S_C1001002_b;
  abcd[429] = 2.0E0*I_ERI_D2y_Dyz_Px_S_C1001002_b;
  abcd[430] = 2.0E0*I_ERI_Dyz_Dyz_Px_S_C1001002_b;
  abcd[431] = 2.0E0*I_ERI_D2z_Dyz_Px_S_C1001002_b;
  abcd[438] = 2.0E0*I_ERI_D2x_Dxy_Py_S_C1001002_b;
  abcd[439] = 2.0E0*I_ERI_Dxy_Dxy_Py_S_C1001002_b;
  abcd[440] = 2.0E0*I_ERI_Dxz_Dxy_Py_S_C1001002_b;
  abcd[441] = 2.0E0*I_ERI_D2y_Dxy_Py_S_C1001002_b;
  abcd[442] = 2.0E0*I_ERI_Dyz_Dxy_Py_S_C1001002_b;
  abcd[443] = 2.0E0*I_ERI_D2z_Dxy_Py_S_C1001002_b;
  abcd[444] = 2.0E0*I_ERI_D2x_D2y_Py_S_C1001002_b-1*I_ERI_D2x_S_Py_S_C1001002;
  abcd[445] = 2.0E0*I_ERI_Dxy_D2y_Py_S_C1001002_b-1*I_ERI_Dxy_S_Py_S_C1001002;
  abcd[446] = 2.0E0*I_ERI_Dxz_D2y_Py_S_C1001002_b-1*I_ERI_Dxz_S_Py_S_C1001002;
  abcd[447] = 2.0E0*I_ERI_D2y_D2y_Py_S_C1001002_b-1*I_ERI_D2y_S_Py_S_C1001002;
  abcd[448] = 2.0E0*I_ERI_Dyz_D2y_Py_S_C1001002_b-1*I_ERI_Dyz_S_Py_S_C1001002;
  abcd[449] = 2.0E0*I_ERI_D2z_D2y_Py_S_C1001002_b-1*I_ERI_D2z_S_Py_S_C1001002;
  abcd[450] = 2.0E0*I_ERI_D2x_Dyz_Py_S_C1001002_b;
  abcd[451] = 2.0E0*I_ERI_Dxy_Dyz_Py_S_C1001002_b;
  abcd[452] = 2.0E0*I_ERI_Dxz_Dyz_Py_S_C1001002_b;
  abcd[453] = 2.0E0*I_ERI_D2y_Dyz_Py_S_C1001002_b;
  abcd[454] = 2.0E0*I_ERI_Dyz_Dyz_Py_S_C1001002_b;
  abcd[455] = 2.0E0*I_ERI_D2z_Dyz_Py_S_C1001002_b;
  abcd[462] = 2.0E0*I_ERI_D2x_Dxy_Pz_S_C1001002_b;
  abcd[463] = 2.0E0*I_ERI_Dxy_Dxy_Pz_S_C1001002_b;
  abcd[464] = 2.0E0*I_ERI_Dxz_Dxy_Pz_S_C1001002_b;
  abcd[465] = 2.0E0*I_ERI_D2y_Dxy_Pz_S_C1001002_b;
  abcd[466] = 2.0E0*I_ERI_Dyz_Dxy_Pz_S_C1001002_b;
  abcd[467] = 2.0E0*I_ERI_D2z_Dxy_Pz_S_C1001002_b;
  abcd[468] = 2.0E0*I_ERI_D2x_D2y_Pz_S_C1001002_b-1*I_ERI_D2x_S_Pz_S_C1001002;
  abcd[469] = 2.0E0*I_ERI_Dxy_D2y_Pz_S_C1001002_b-1*I_ERI_Dxy_S_Pz_S_C1001002;
  abcd[470] = 2.0E0*I_ERI_Dxz_D2y_Pz_S_C1001002_b-1*I_ERI_Dxz_S_Pz_S_C1001002;
  abcd[471] = 2.0E0*I_ERI_D2y_D2y_Pz_S_C1001002_b-1*I_ERI_D2y_S_Pz_S_C1001002;
  abcd[472] = 2.0E0*I_ERI_Dyz_D2y_Pz_S_C1001002_b-1*I_ERI_Dyz_S_Pz_S_C1001002;
  abcd[473] = 2.0E0*I_ERI_D2z_D2y_Pz_S_C1001002_b-1*I_ERI_D2z_S_Pz_S_C1001002;
  abcd[474] = 2.0E0*I_ERI_D2x_Dyz_Pz_S_C1001002_b;
  abcd[475] = 2.0E0*I_ERI_Dxy_Dyz_Pz_S_C1001002_b;
  abcd[476] = 2.0E0*I_ERI_Dxz_Dyz_Pz_S_C1001002_b;
  abcd[477] = 2.0E0*I_ERI_D2y_Dyz_Pz_S_C1001002_b;
  abcd[478] = 2.0E0*I_ERI_Dyz_Dyz_Pz_S_C1001002_b;
  abcd[479] = 2.0E0*I_ERI_D2z_Dyz_Pz_S_C1001002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C2_b
   ************************************************************/
  abcd[480] = 2.0E0*I_ERI_D2x_Pz_S_S_C2_b;
  abcd[481] = 2.0E0*I_ERI_Dxy_Pz_S_S_C2_b;
  abcd[482] = 2.0E0*I_ERI_Dxz_Pz_S_S_C2_b;
  abcd[483] = 2.0E0*I_ERI_D2y_Pz_S_S_C2_b;
  abcd[484] = 2.0E0*I_ERI_Dyz_Pz_S_S_C2_b;
  abcd[485] = 2.0E0*I_ERI_D2z_Pz_S_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_C1002_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1002
   ************************************************************/
  abcd[486] = 2.0E0*I_ERI_D2x_Dxz_S_S_C1002_b;
  abcd[487] = 2.0E0*I_ERI_Dxy_Dxz_S_S_C1002_b;
  abcd[488] = 2.0E0*I_ERI_Dxz_Dxz_S_S_C1002_b;
  abcd[489] = 2.0E0*I_ERI_D2y_Dxz_S_S_C1002_b;
  abcd[490] = 2.0E0*I_ERI_Dyz_Dxz_S_S_C1002_b;
  abcd[491] = 2.0E0*I_ERI_D2z_Dxz_S_S_C1002_b;
  abcd[492] = 2.0E0*I_ERI_D2x_Dyz_S_S_C1002_b;
  abcd[493] = 2.0E0*I_ERI_Dxy_Dyz_S_S_C1002_b;
  abcd[494] = 2.0E0*I_ERI_Dxz_Dyz_S_S_C1002_b;
  abcd[495] = 2.0E0*I_ERI_D2y_Dyz_S_S_C1002_b;
  abcd[496] = 2.0E0*I_ERI_Dyz_Dyz_S_S_C1002_b;
  abcd[497] = 2.0E0*I_ERI_D2z_Dyz_S_S_C1002_b;
  abcd[498] = 2.0E0*I_ERI_D2x_D2z_S_S_C1002_b-1*I_ERI_D2x_S_S_S_C1002;
  abcd[499] = 2.0E0*I_ERI_Dxy_D2z_S_S_C1002_b-1*I_ERI_Dxy_S_S_S_C1002;
  abcd[500] = 2.0E0*I_ERI_Dxz_D2z_S_S_C1002_b-1*I_ERI_Dxz_S_S_S_C1002;
  abcd[501] = 2.0E0*I_ERI_D2y_D2z_S_S_C1002_b-1*I_ERI_D2y_S_S_S_C1002;
  abcd[502] = 2.0E0*I_ERI_Dyz_D2z_S_S_C1002_b-1*I_ERI_Dyz_S_S_S_C1002;
  abcd[503] = 2.0E0*I_ERI_D2z_D2z_S_S_C1002_b-1*I_ERI_D2z_S_S_S_C1002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1000002_b
   ************************************************************/
  abcd[504] = 2.0E0*I_ERI_D2x_Pz_Px_S_C1000002_b;
  abcd[505] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1000002_b;
  abcd[506] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1000002_b;
  abcd[507] = 2.0E0*I_ERI_D2y_Pz_Px_S_C1000002_b;
  abcd[508] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1000002_b;
  abcd[509] = 2.0E0*I_ERI_D2z_Pz_Px_S_C1000002_b;
  abcd[528] = 2.0E0*I_ERI_D2x_Pz_Py_S_C1000002_b;
  abcd[529] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1000002_b;
  abcd[530] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1000002_b;
  abcd[531] = 2.0E0*I_ERI_D2y_Pz_Py_S_C1000002_b;
  abcd[532] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1000002_b;
  abcd[533] = 2.0E0*I_ERI_D2z_Pz_Py_S_C1000002_b;
  abcd[552] = 2.0E0*I_ERI_D2x_Pz_Pz_S_C1000002_b;
  abcd[553] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1000002_b;
  abcd[554] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1000002_b;
  abcd[555] = 2.0E0*I_ERI_D2y_Pz_Pz_S_C1000002_b;
  abcd[556] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1000002_b;
  abcd[557] = 2.0E0*I_ERI_D2z_Pz_Pz_S_C1000002_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_P_S_C1001002_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C1001002
   ************************************************************/
  abcd[510] = 2.0E0*I_ERI_D2x_Dxz_Px_S_C1001002_b;
  abcd[511] = 2.0E0*I_ERI_Dxy_Dxz_Px_S_C1001002_b;
  abcd[512] = 2.0E0*I_ERI_Dxz_Dxz_Px_S_C1001002_b;
  abcd[513] = 2.0E0*I_ERI_D2y_Dxz_Px_S_C1001002_b;
  abcd[514] = 2.0E0*I_ERI_Dyz_Dxz_Px_S_C1001002_b;
  abcd[515] = 2.0E0*I_ERI_D2z_Dxz_Px_S_C1001002_b;
  abcd[516] = 2.0E0*I_ERI_D2x_Dyz_Px_S_C1001002_b;
  abcd[517] = 2.0E0*I_ERI_Dxy_Dyz_Px_S_C1001002_b;
  abcd[518] = 2.0E0*I_ERI_Dxz_Dyz_Px_S_C1001002_b;
  abcd[519] = 2.0E0*I_ERI_D2y_Dyz_Px_S_C1001002_b;
  abcd[520] = 2.0E0*I_ERI_Dyz_Dyz_Px_S_C1001002_b;
  abcd[521] = 2.0E0*I_ERI_D2z_Dyz_Px_S_C1001002_b;
  abcd[522] = 2.0E0*I_ERI_D2x_D2z_Px_S_C1001002_b-1*I_ERI_D2x_S_Px_S_C1001002;
  abcd[523] = 2.0E0*I_ERI_Dxy_D2z_Px_S_C1001002_b-1*I_ERI_Dxy_S_Px_S_C1001002;
  abcd[524] = 2.0E0*I_ERI_Dxz_D2z_Px_S_C1001002_b-1*I_ERI_Dxz_S_Px_S_C1001002;
  abcd[525] = 2.0E0*I_ERI_D2y_D2z_Px_S_C1001002_b-1*I_ERI_D2y_S_Px_S_C1001002;
  abcd[526] = 2.0E0*I_ERI_Dyz_D2z_Px_S_C1001002_b-1*I_ERI_Dyz_S_Px_S_C1001002;
  abcd[527] = 2.0E0*I_ERI_D2z_D2z_Px_S_C1001002_b-1*I_ERI_D2z_S_Px_S_C1001002;
  abcd[534] = 2.0E0*I_ERI_D2x_Dxz_Py_S_C1001002_b;
  abcd[535] = 2.0E0*I_ERI_Dxy_Dxz_Py_S_C1001002_b;
  abcd[536] = 2.0E0*I_ERI_Dxz_Dxz_Py_S_C1001002_b;
  abcd[537] = 2.0E0*I_ERI_D2y_Dxz_Py_S_C1001002_b;
  abcd[538] = 2.0E0*I_ERI_Dyz_Dxz_Py_S_C1001002_b;
  abcd[539] = 2.0E0*I_ERI_D2z_Dxz_Py_S_C1001002_b;
  abcd[540] = 2.0E0*I_ERI_D2x_Dyz_Py_S_C1001002_b;
  abcd[541] = 2.0E0*I_ERI_Dxy_Dyz_Py_S_C1001002_b;
  abcd[542] = 2.0E0*I_ERI_Dxz_Dyz_Py_S_C1001002_b;
  abcd[543] = 2.0E0*I_ERI_D2y_Dyz_Py_S_C1001002_b;
  abcd[544] = 2.0E0*I_ERI_Dyz_Dyz_Py_S_C1001002_b;
  abcd[545] = 2.0E0*I_ERI_D2z_Dyz_Py_S_C1001002_b;
  abcd[546] = 2.0E0*I_ERI_D2x_D2z_Py_S_C1001002_b-1*I_ERI_D2x_S_Py_S_C1001002;
  abcd[547] = 2.0E0*I_ERI_Dxy_D2z_Py_S_C1001002_b-1*I_ERI_Dxy_S_Py_S_C1001002;
  abcd[548] = 2.0E0*I_ERI_Dxz_D2z_Py_S_C1001002_b-1*I_ERI_Dxz_S_Py_S_C1001002;
  abcd[549] = 2.0E0*I_ERI_D2y_D2z_Py_S_C1001002_b-1*I_ERI_D2y_S_Py_S_C1001002;
  abcd[550] = 2.0E0*I_ERI_Dyz_D2z_Py_S_C1001002_b-1*I_ERI_Dyz_S_Py_S_C1001002;
  abcd[551] = 2.0E0*I_ERI_D2z_D2z_Py_S_C1001002_b-1*I_ERI_D2z_S_Py_S_C1001002;
  abcd[558] = 2.0E0*I_ERI_D2x_Dxz_Pz_S_C1001002_b;
  abcd[559] = 2.0E0*I_ERI_Dxy_Dxz_Pz_S_C1001002_b;
  abcd[560] = 2.0E0*I_ERI_Dxz_Dxz_Pz_S_C1001002_b;
  abcd[561] = 2.0E0*I_ERI_D2y_Dxz_Pz_S_C1001002_b;
  abcd[562] = 2.0E0*I_ERI_Dyz_Dxz_Pz_S_C1001002_b;
  abcd[563] = 2.0E0*I_ERI_D2z_Dxz_Pz_S_C1001002_b;
  abcd[564] = 2.0E0*I_ERI_D2x_Dyz_Pz_S_C1001002_b;
  abcd[565] = 2.0E0*I_ERI_Dxy_Dyz_Pz_S_C1001002_b;
  abcd[566] = 2.0E0*I_ERI_Dxz_Dyz_Pz_S_C1001002_b;
  abcd[567] = 2.0E0*I_ERI_D2y_Dyz_Pz_S_C1001002_b;
  abcd[568] = 2.0E0*I_ERI_Dyz_Dyz_Pz_S_C1001002_b;
  abcd[569] = 2.0E0*I_ERI_D2z_Dyz_Pz_S_C1001002_b;
  abcd[570] = 2.0E0*I_ERI_D2x_D2z_Pz_S_C1001002_b-1*I_ERI_D2x_S_Pz_S_C1001002;
  abcd[571] = 2.0E0*I_ERI_Dxy_D2z_Pz_S_C1001002_b-1*I_ERI_Dxy_S_Pz_S_C1001002;
  abcd[572] = 2.0E0*I_ERI_Dxz_D2z_Pz_S_C1001002_b-1*I_ERI_Dxz_S_Pz_S_C1001002;
  abcd[573] = 2.0E0*I_ERI_D2y_D2z_Pz_S_C1001002_b-1*I_ERI_D2y_S_Pz_S_C1001002;
  abcd[574] = 2.0E0*I_ERI_Dyz_D2z_Pz_S_C1001002_b-1*I_ERI_Dyz_S_Pz_S_C1001002;
  abcd[575] = 2.0E0*I_ERI_D2z_D2z_Pz_S_C1001002_b-1*I_ERI_D2z_S_Pz_S_C1001002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C2_c
   ************************************************************/
  abcd[576] = 2.0E0*I_ERI_D2x_S_Px_S_C2_c;
  abcd[577] = 2.0E0*I_ERI_Dxy_S_Px_S_C2_c;
  abcd[578] = 2.0E0*I_ERI_Dxz_S_Px_S_C2_c;
  abcd[579] = 2.0E0*I_ERI_D2y_S_Px_S_C2_c;
  abcd[580] = 2.0E0*I_ERI_Dyz_S_Px_S_C2_c;
  abcd[581] = 2.0E0*I_ERI_D2z_S_Px_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1002_c
   ************************************************************/
  abcd[582] = 2.0E0*I_ERI_D2x_Px_Px_S_C1002_c;
  abcd[583] = 2.0E0*I_ERI_Dxy_Px_Px_S_C1002_c;
  abcd[584] = 2.0E0*I_ERI_Dxz_Px_Px_S_C1002_c;
  abcd[585] = 2.0E0*I_ERI_D2y_Px_Px_S_C1002_c;
  abcd[586] = 2.0E0*I_ERI_Dyz_Px_Px_S_C1002_c;
  abcd[587] = 2.0E0*I_ERI_D2z_Px_Px_S_C1002_c;
  abcd[588] = 2.0E0*I_ERI_D2x_Py_Px_S_C1002_c;
  abcd[589] = 2.0E0*I_ERI_Dxy_Py_Px_S_C1002_c;
  abcd[590] = 2.0E0*I_ERI_Dxz_Py_Px_S_C1002_c;
  abcd[591] = 2.0E0*I_ERI_D2y_Py_Px_S_C1002_c;
  abcd[592] = 2.0E0*I_ERI_Dyz_Py_Px_S_C1002_c;
  abcd[593] = 2.0E0*I_ERI_D2z_Py_Px_S_C1002_c;
  abcd[594] = 2.0E0*I_ERI_D2x_Pz_Px_S_C1002_c;
  abcd[595] = 2.0E0*I_ERI_Dxy_Pz_Px_S_C1002_c;
  abcd[596] = 2.0E0*I_ERI_Dxz_Pz_Px_S_C1002_c;
  abcd[597] = 2.0E0*I_ERI_D2y_Pz_Px_S_C1002_c;
  abcd[598] = 2.0E0*I_ERI_Dyz_Pz_Px_S_C1002_c;
  abcd[599] = 2.0E0*I_ERI_D2z_Pz_Px_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000002_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000002
   ************************************************************/
  abcd[600] = 2.0E0*I_ERI_D2x_S_D2x_S_C1000002_c-1*I_ERI_D2x_S_S_S_C1000002;
  abcd[601] = 2.0E0*I_ERI_Dxy_S_D2x_S_C1000002_c-1*I_ERI_Dxy_S_S_S_C1000002;
  abcd[602] = 2.0E0*I_ERI_Dxz_S_D2x_S_C1000002_c-1*I_ERI_Dxz_S_S_S_C1000002;
  abcd[603] = 2.0E0*I_ERI_D2y_S_D2x_S_C1000002_c-1*I_ERI_D2y_S_S_S_C1000002;
  abcd[604] = 2.0E0*I_ERI_Dyz_S_D2x_S_C1000002_c-1*I_ERI_Dyz_S_S_S_C1000002;
  abcd[605] = 2.0E0*I_ERI_D2z_S_D2x_S_C1000002_c-1*I_ERI_D2z_S_S_S_C1000002;
  abcd[624] = 2.0E0*I_ERI_D2x_S_Dxy_S_C1000002_c;
  abcd[625] = 2.0E0*I_ERI_Dxy_S_Dxy_S_C1000002_c;
  abcd[626] = 2.0E0*I_ERI_Dxz_S_Dxy_S_C1000002_c;
  abcd[627] = 2.0E0*I_ERI_D2y_S_Dxy_S_C1000002_c;
  abcd[628] = 2.0E0*I_ERI_Dyz_S_Dxy_S_C1000002_c;
  abcd[629] = 2.0E0*I_ERI_D2z_S_Dxy_S_C1000002_c;
  abcd[648] = 2.0E0*I_ERI_D2x_S_Dxz_S_C1000002_c;
  abcd[649] = 2.0E0*I_ERI_Dxy_S_Dxz_S_C1000002_c;
  abcd[650] = 2.0E0*I_ERI_Dxz_S_Dxz_S_C1000002_c;
  abcd[651] = 2.0E0*I_ERI_D2y_S_Dxz_S_C1000002_c;
  abcd[652] = 2.0E0*I_ERI_Dyz_S_Dxz_S_C1000002_c;
  abcd[653] = 2.0E0*I_ERI_D2z_S_Dxz_S_C1000002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001002_c
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001002
   ************************************************************/
  abcd[606] = 2.0E0*I_ERI_D2x_Px_D2x_S_C1001002_c-1*I_ERI_D2x_Px_S_S_C1001002;
  abcd[607] = 2.0E0*I_ERI_Dxy_Px_D2x_S_C1001002_c-1*I_ERI_Dxy_Px_S_S_C1001002;
  abcd[608] = 2.0E0*I_ERI_Dxz_Px_D2x_S_C1001002_c-1*I_ERI_Dxz_Px_S_S_C1001002;
  abcd[609] = 2.0E0*I_ERI_D2y_Px_D2x_S_C1001002_c-1*I_ERI_D2y_Px_S_S_C1001002;
  abcd[610] = 2.0E0*I_ERI_Dyz_Px_D2x_S_C1001002_c-1*I_ERI_Dyz_Px_S_S_C1001002;
  abcd[611] = 2.0E0*I_ERI_D2z_Px_D2x_S_C1001002_c-1*I_ERI_D2z_Px_S_S_C1001002;
  abcd[612] = 2.0E0*I_ERI_D2x_Py_D2x_S_C1001002_c-1*I_ERI_D2x_Py_S_S_C1001002;
  abcd[613] = 2.0E0*I_ERI_Dxy_Py_D2x_S_C1001002_c-1*I_ERI_Dxy_Py_S_S_C1001002;
  abcd[614] = 2.0E0*I_ERI_Dxz_Py_D2x_S_C1001002_c-1*I_ERI_Dxz_Py_S_S_C1001002;
  abcd[615] = 2.0E0*I_ERI_D2y_Py_D2x_S_C1001002_c-1*I_ERI_D2y_Py_S_S_C1001002;
  abcd[616] = 2.0E0*I_ERI_Dyz_Py_D2x_S_C1001002_c-1*I_ERI_Dyz_Py_S_S_C1001002;
  abcd[617] = 2.0E0*I_ERI_D2z_Py_D2x_S_C1001002_c-1*I_ERI_D2z_Py_S_S_C1001002;
  abcd[618] = 2.0E0*I_ERI_D2x_Pz_D2x_S_C1001002_c-1*I_ERI_D2x_Pz_S_S_C1001002;
  abcd[619] = 2.0E0*I_ERI_Dxy_Pz_D2x_S_C1001002_c-1*I_ERI_Dxy_Pz_S_S_C1001002;
  abcd[620] = 2.0E0*I_ERI_Dxz_Pz_D2x_S_C1001002_c-1*I_ERI_Dxz_Pz_S_S_C1001002;
  abcd[621] = 2.0E0*I_ERI_D2y_Pz_D2x_S_C1001002_c-1*I_ERI_D2y_Pz_S_S_C1001002;
  abcd[622] = 2.0E0*I_ERI_Dyz_Pz_D2x_S_C1001002_c-1*I_ERI_Dyz_Pz_S_S_C1001002;
  abcd[623] = 2.0E0*I_ERI_D2z_Pz_D2x_S_C1001002_c-1*I_ERI_D2z_Pz_S_S_C1001002;
  abcd[630] = 2.0E0*I_ERI_D2x_Px_Dxy_S_C1001002_c;
  abcd[631] = 2.0E0*I_ERI_Dxy_Px_Dxy_S_C1001002_c;
  abcd[632] = 2.0E0*I_ERI_Dxz_Px_Dxy_S_C1001002_c;
  abcd[633] = 2.0E0*I_ERI_D2y_Px_Dxy_S_C1001002_c;
  abcd[634] = 2.0E0*I_ERI_Dyz_Px_Dxy_S_C1001002_c;
  abcd[635] = 2.0E0*I_ERI_D2z_Px_Dxy_S_C1001002_c;
  abcd[636] = 2.0E0*I_ERI_D2x_Py_Dxy_S_C1001002_c;
  abcd[637] = 2.0E0*I_ERI_Dxy_Py_Dxy_S_C1001002_c;
  abcd[638] = 2.0E0*I_ERI_Dxz_Py_Dxy_S_C1001002_c;
  abcd[639] = 2.0E0*I_ERI_D2y_Py_Dxy_S_C1001002_c;
  abcd[640] = 2.0E0*I_ERI_Dyz_Py_Dxy_S_C1001002_c;
  abcd[641] = 2.0E0*I_ERI_D2z_Py_Dxy_S_C1001002_c;
  abcd[642] = 2.0E0*I_ERI_D2x_Pz_Dxy_S_C1001002_c;
  abcd[643] = 2.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001002_c;
  abcd[644] = 2.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001002_c;
  abcd[645] = 2.0E0*I_ERI_D2y_Pz_Dxy_S_C1001002_c;
  abcd[646] = 2.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001002_c;
  abcd[647] = 2.0E0*I_ERI_D2z_Pz_Dxy_S_C1001002_c;
  abcd[654] = 2.0E0*I_ERI_D2x_Px_Dxz_S_C1001002_c;
  abcd[655] = 2.0E0*I_ERI_Dxy_Px_Dxz_S_C1001002_c;
  abcd[656] = 2.0E0*I_ERI_Dxz_Px_Dxz_S_C1001002_c;
  abcd[657] = 2.0E0*I_ERI_D2y_Px_Dxz_S_C1001002_c;
  abcd[658] = 2.0E0*I_ERI_Dyz_Px_Dxz_S_C1001002_c;
  abcd[659] = 2.0E0*I_ERI_D2z_Px_Dxz_S_C1001002_c;
  abcd[660] = 2.0E0*I_ERI_D2x_Py_Dxz_S_C1001002_c;
  abcd[661] = 2.0E0*I_ERI_Dxy_Py_Dxz_S_C1001002_c;
  abcd[662] = 2.0E0*I_ERI_Dxz_Py_Dxz_S_C1001002_c;
  abcd[663] = 2.0E0*I_ERI_D2y_Py_Dxz_S_C1001002_c;
  abcd[664] = 2.0E0*I_ERI_Dyz_Py_Dxz_S_C1001002_c;
  abcd[665] = 2.0E0*I_ERI_D2z_Py_Dxz_S_C1001002_c;
  abcd[666] = 2.0E0*I_ERI_D2x_Pz_Dxz_S_C1001002_c;
  abcd[667] = 2.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001002_c;
  abcd[668] = 2.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001002_c;
  abcd[669] = 2.0E0*I_ERI_D2y_Pz_Dxz_S_C1001002_c;
  abcd[670] = 2.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001002_c;
  abcd[671] = 2.0E0*I_ERI_D2z_Pz_Dxz_S_C1001002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C2_c
   ************************************************************/
  abcd[672] = 2.0E0*I_ERI_D2x_S_Py_S_C2_c;
  abcd[673] = 2.0E0*I_ERI_Dxy_S_Py_S_C2_c;
  abcd[674] = 2.0E0*I_ERI_Dxz_S_Py_S_C2_c;
  abcd[675] = 2.0E0*I_ERI_D2y_S_Py_S_C2_c;
  abcd[676] = 2.0E0*I_ERI_Dyz_S_Py_S_C2_c;
  abcd[677] = 2.0E0*I_ERI_D2z_S_Py_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1002_c
   ************************************************************/
  abcd[678] = 2.0E0*I_ERI_D2x_Px_Py_S_C1002_c;
  abcd[679] = 2.0E0*I_ERI_Dxy_Px_Py_S_C1002_c;
  abcd[680] = 2.0E0*I_ERI_Dxz_Px_Py_S_C1002_c;
  abcd[681] = 2.0E0*I_ERI_D2y_Px_Py_S_C1002_c;
  abcd[682] = 2.0E0*I_ERI_Dyz_Px_Py_S_C1002_c;
  abcd[683] = 2.0E0*I_ERI_D2z_Px_Py_S_C1002_c;
  abcd[684] = 2.0E0*I_ERI_D2x_Py_Py_S_C1002_c;
  abcd[685] = 2.0E0*I_ERI_Dxy_Py_Py_S_C1002_c;
  abcd[686] = 2.0E0*I_ERI_Dxz_Py_Py_S_C1002_c;
  abcd[687] = 2.0E0*I_ERI_D2y_Py_Py_S_C1002_c;
  abcd[688] = 2.0E0*I_ERI_Dyz_Py_Py_S_C1002_c;
  abcd[689] = 2.0E0*I_ERI_D2z_Py_Py_S_C1002_c;
  abcd[690] = 2.0E0*I_ERI_D2x_Pz_Py_S_C1002_c;
  abcd[691] = 2.0E0*I_ERI_Dxy_Pz_Py_S_C1002_c;
  abcd[692] = 2.0E0*I_ERI_Dxz_Pz_Py_S_C1002_c;
  abcd[693] = 2.0E0*I_ERI_D2y_Pz_Py_S_C1002_c;
  abcd[694] = 2.0E0*I_ERI_Dyz_Pz_Py_S_C1002_c;
  abcd[695] = 2.0E0*I_ERI_D2z_Pz_Py_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000002_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000002
   ************************************************************/
  abcd[696] = 2.0E0*I_ERI_D2x_S_Dxy_S_C1000002_c;
  abcd[697] = 2.0E0*I_ERI_Dxy_S_Dxy_S_C1000002_c;
  abcd[698] = 2.0E0*I_ERI_Dxz_S_Dxy_S_C1000002_c;
  abcd[699] = 2.0E0*I_ERI_D2y_S_Dxy_S_C1000002_c;
  abcd[700] = 2.0E0*I_ERI_Dyz_S_Dxy_S_C1000002_c;
  abcd[701] = 2.0E0*I_ERI_D2z_S_Dxy_S_C1000002_c;
  abcd[720] = 2.0E0*I_ERI_D2x_S_D2y_S_C1000002_c-1*I_ERI_D2x_S_S_S_C1000002;
  abcd[721] = 2.0E0*I_ERI_Dxy_S_D2y_S_C1000002_c-1*I_ERI_Dxy_S_S_S_C1000002;
  abcd[722] = 2.0E0*I_ERI_Dxz_S_D2y_S_C1000002_c-1*I_ERI_Dxz_S_S_S_C1000002;
  abcd[723] = 2.0E0*I_ERI_D2y_S_D2y_S_C1000002_c-1*I_ERI_D2y_S_S_S_C1000002;
  abcd[724] = 2.0E0*I_ERI_Dyz_S_D2y_S_C1000002_c-1*I_ERI_Dyz_S_S_S_C1000002;
  abcd[725] = 2.0E0*I_ERI_D2z_S_D2y_S_C1000002_c-1*I_ERI_D2z_S_S_S_C1000002;
  abcd[744] = 2.0E0*I_ERI_D2x_S_Dyz_S_C1000002_c;
  abcd[745] = 2.0E0*I_ERI_Dxy_S_Dyz_S_C1000002_c;
  abcd[746] = 2.0E0*I_ERI_Dxz_S_Dyz_S_C1000002_c;
  abcd[747] = 2.0E0*I_ERI_D2y_S_Dyz_S_C1000002_c;
  abcd[748] = 2.0E0*I_ERI_Dyz_S_Dyz_S_C1000002_c;
  abcd[749] = 2.0E0*I_ERI_D2z_S_Dyz_S_C1000002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001002_c
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001002
   ************************************************************/
  abcd[702] = 2.0E0*I_ERI_D2x_Px_Dxy_S_C1001002_c;
  abcd[703] = 2.0E0*I_ERI_Dxy_Px_Dxy_S_C1001002_c;
  abcd[704] = 2.0E0*I_ERI_Dxz_Px_Dxy_S_C1001002_c;
  abcd[705] = 2.0E0*I_ERI_D2y_Px_Dxy_S_C1001002_c;
  abcd[706] = 2.0E0*I_ERI_Dyz_Px_Dxy_S_C1001002_c;
  abcd[707] = 2.0E0*I_ERI_D2z_Px_Dxy_S_C1001002_c;
  abcd[708] = 2.0E0*I_ERI_D2x_Py_Dxy_S_C1001002_c;
  abcd[709] = 2.0E0*I_ERI_Dxy_Py_Dxy_S_C1001002_c;
  abcd[710] = 2.0E0*I_ERI_Dxz_Py_Dxy_S_C1001002_c;
  abcd[711] = 2.0E0*I_ERI_D2y_Py_Dxy_S_C1001002_c;
  abcd[712] = 2.0E0*I_ERI_Dyz_Py_Dxy_S_C1001002_c;
  abcd[713] = 2.0E0*I_ERI_D2z_Py_Dxy_S_C1001002_c;
  abcd[714] = 2.0E0*I_ERI_D2x_Pz_Dxy_S_C1001002_c;
  abcd[715] = 2.0E0*I_ERI_Dxy_Pz_Dxy_S_C1001002_c;
  abcd[716] = 2.0E0*I_ERI_Dxz_Pz_Dxy_S_C1001002_c;
  abcd[717] = 2.0E0*I_ERI_D2y_Pz_Dxy_S_C1001002_c;
  abcd[718] = 2.0E0*I_ERI_Dyz_Pz_Dxy_S_C1001002_c;
  abcd[719] = 2.0E0*I_ERI_D2z_Pz_Dxy_S_C1001002_c;
  abcd[726] = 2.0E0*I_ERI_D2x_Px_D2y_S_C1001002_c-1*I_ERI_D2x_Px_S_S_C1001002;
  abcd[727] = 2.0E0*I_ERI_Dxy_Px_D2y_S_C1001002_c-1*I_ERI_Dxy_Px_S_S_C1001002;
  abcd[728] = 2.0E0*I_ERI_Dxz_Px_D2y_S_C1001002_c-1*I_ERI_Dxz_Px_S_S_C1001002;
  abcd[729] = 2.0E0*I_ERI_D2y_Px_D2y_S_C1001002_c-1*I_ERI_D2y_Px_S_S_C1001002;
  abcd[730] = 2.0E0*I_ERI_Dyz_Px_D2y_S_C1001002_c-1*I_ERI_Dyz_Px_S_S_C1001002;
  abcd[731] = 2.0E0*I_ERI_D2z_Px_D2y_S_C1001002_c-1*I_ERI_D2z_Px_S_S_C1001002;
  abcd[732] = 2.0E0*I_ERI_D2x_Py_D2y_S_C1001002_c-1*I_ERI_D2x_Py_S_S_C1001002;
  abcd[733] = 2.0E0*I_ERI_Dxy_Py_D2y_S_C1001002_c-1*I_ERI_Dxy_Py_S_S_C1001002;
  abcd[734] = 2.0E0*I_ERI_Dxz_Py_D2y_S_C1001002_c-1*I_ERI_Dxz_Py_S_S_C1001002;
  abcd[735] = 2.0E0*I_ERI_D2y_Py_D2y_S_C1001002_c-1*I_ERI_D2y_Py_S_S_C1001002;
  abcd[736] = 2.0E0*I_ERI_Dyz_Py_D2y_S_C1001002_c-1*I_ERI_Dyz_Py_S_S_C1001002;
  abcd[737] = 2.0E0*I_ERI_D2z_Py_D2y_S_C1001002_c-1*I_ERI_D2z_Py_S_S_C1001002;
  abcd[738] = 2.0E0*I_ERI_D2x_Pz_D2y_S_C1001002_c-1*I_ERI_D2x_Pz_S_S_C1001002;
  abcd[739] = 2.0E0*I_ERI_Dxy_Pz_D2y_S_C1001002_c-1*I_ERI_Dxy_Pz_S_S_C1001002;
  abcd[740] = 2.0E0*I_ERI_Dxz_Pz_D2y_S_C1001002_c-1*I_ERI_Dxz_Pz_S_S_C1001002;
  abcd[741] = 2.0E0*I_ERI_D2y_Pz_D2y_S_C1001002_c-1*I_ERI_D2y_Pz_S_S_C1001002;
  abcd[742] = 2.0E0*I_ERI_Dyz_Pz_D2y_S_C1001002_c-1*I_ERI_Dyz_Pz_S_S_C1001002;
  abcd[743] = 2.0E0*I_ERI_D2z_Pz_D2y_S_C1001002_c-1*I_ERI_D2z_Pz_S_S_C1001002;
  abcd[750] = 2.0E0*I_ERI_D2x_Px_Dyz_S_C1001002_c;
  abcd[751] = 2.0E0*I_ERI_Dxy_Px_Dyz_S_C1001002_c;
  abcd[752] = 2.0E0*I_ERI_Dxz_Px_Dyz_S_C1001002_c;
  abcd[753] = 2.0E0*I_ERI_D2y_Px_Dyz_S_C1001002_c;
  abcd[754] = 2.0E0*I_ERI_Dyz_Px_Dyz_S_C1001002_c;
  abcd[755] = 2.0E0*I_ERI_D2z_Px_Dyz_S_C1001002_c;
  abcd[756] = 2.0E0*I_ERI_D2x_Py_Dyz_S_C1001002_c;
  abcd[757] = 2.0E0*I_ERI_Dxy_Py_Dyz_S_C1001002_c;
  abcd[758] = 2.0E0*I_ERI_Dxz_Py_Dyz_S_C1001002_c;
  abcd[759] = 2.0E0*I_ERI_D2y_Py_Dyz_S_C1001002_c;
  abcd[760] = 2.0E0*I_ERI_Dyz_Py_Dyz_S_C1001002_c;
  abcd[761] = 2.0E0*I_ERI_D2z_Py_Dyz_S_C1001002_c;
  abcd[762] = 2.0E0*I_ERI_D2x_Pz_Dyz_S_C1001002_c;
  abcd[763] = 2.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001002_c;
  abcd[764] = 2.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001002_c;
  abcd[765] = 2.0E0*I_ERI_D2y_Pz_Dyz_S_C1001002_c;
  abcd[766] = 2.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001002_c;
  abcd[767] = 2.0E0*I_ERI_D2z_Pz_Dyz_S_C1001002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_C2_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_C2_c
   ************************************************************/
  abcd[768] = 2.0E0*I_ERI_D2x_S_Pz_S_C2_c;
  abcd[769] = 2.0E0*I_ERI_Dxy_S_Pz_S_C2_c;
  abcd[770] = 2.0E0*I_ERI_Dxz_S_Pz_S_C2_c;
  abcd[771] = 2.0E0*I_ERI_D2y_S_Pz_S_C2_c;
  abcd[772] = 2.0E0*I_ERI_Dyz_S_Pz_S_C2_c;
  abcd[773] = 2.0E0*I_ERI_D2z_S_Pz_S_C2_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_C1002_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_C1002_c
   ************************************************************/
  abcd[774] = 2.0E0*I_ERI_D2x_Px_Pz_S_C1002_c;
  abcd[775] = 2.0E0*I_ERI_Dxy_Px_Pz_S_C1002_c;
  abcd[776] = 2.0E0*I_ERI_Dxz_Px_Pz_S_C1002_c;
  abcd[777] = 2.0E0*I_ERI_D2y_Px_Pz_S_C1002_c;
  abcd[778] = 2.0E0*I_ERI_Dyz_Px_Pz_S_C1002_c;
  abcd[779] = 2.0E0*I_ERI_D2z_Px_Pz_S_C1002_c;
  abcd[780] = 2.0E0*I_ERI_D2x_Py_Pz_S_C1002_c;
  abcd[781] = 2.0E0*I_ERI_Dxy_Py_Pz_S_C1002_c;
  abcd[782] = 2.0E0*I_ERI_Dxz_Py_Pz_S_C1002_c;
  abcd[783] = 2.0E0*I_ERI_D2y_Py_Pz_S_C1002_c;
  abcd[784] = 2.0E0*I_ERI_Dyz_Py_Pz_S_C1002_c;
  abcd[785] = 2.0E0*I_ERI_D2z_Py_Pz_S_C1002_c;
  abcd[786] = 2.0E0*I_ERI_D2x_Pz_Pz_S_C1002_c;
  abcd[787] = 2.0E0*I_ERI_Dxy_Pz_Pz_S_C1002_c;
  abcd[788] = 2.0E0*I_ERI_Dxz_Pz_Pz_S_C1002_c;
  abcd[789] = 2.0E0*I_ERI_D2y_Pz_Pz_S_C1002_c;
  abcd[790] = 2.0E0*I_ERI_Dyz_Pz_Pz_S_C1002_c;
  abcd[791] = 2.0E0*I_ERI_D2z_Pz_Pz_S_C1002_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_P_S_C1000002_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1000002_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S_C1000002
   ************************************************************/
  abcd[792] = 2.0E0*I_ERI_D2x_S_Dxz_S_C1000002_c;
  abcd[793] = 2.0E0*I_ERI_Dxy_S_Dxz_S_C1000002_c;
  abcd[794] = 2.0E0*I_ERI_Dxz_S_Dxz_S_C1000002_c;
  abcd[795] = 2.0E0*I_ERI_D2y_S_Dxz_S_C1000002_c;
  abcd[796] = 2.0E0*I_ERI_Dyz_S_Dxz_S_C1000002_c;
  abcd[797] = 2.0E0*I_ERI_D2z_S_Dxz_S_C1000002_c;
  abcd[816] = 2.0E0*I_ERI_D2x_S_Dyz_S_C1000002_c;
  abcd[817] = 2.0E0*I_ERI_Dxy_S_Dyz_S_C1000002_c;
  abcd[818] = 2.0E0*I_ERI_Dxz_S_Dyz_S_C1000002_c;
  abcd[819] = 2.0E0*I_ERI_D2y_S_Dyz_S_C1000002_c;
  abcd[820] = 2.0E0*I_ERI_Dyz_S_Dyz_S_C1000002_c;
  abcd[821] = 2.0E0*I_ERI_D2z_S_Dyz_S_C1000002_c;
  abcd[840] = 2.0E0*I_ERI_D2x_S_D2z_S_C1000002_c-1*I_ERI_D2x_S_S_S_C1000002;
  abcd[841] = 2.0E0*I_ERI_Dxy_S_D2z_S_C1000002_c-1*I_ERI_Dxy_S_S_S_C1000002;
  abcd[842] = 2.0E0*I_ERI_Dxz_S_D2z_S_C1000002_c-1*I_ERI_Dxz_S_S_S_C1000002;
  abcd[843] = 2.0E0*I_ERI_D2y_S_D2z_S_C1000002_c-1*I_ERI_D2y_S_S_S_C1000002;
  abcd[844] = 2.0E0*I_ERI_Dyz_S_D2z_S_C1000002_c-1*I_ERI_Dyz_S_S_S_C1000002;
  abcd[845] = 2.0E0*I_ERI_D2z_S_D2z_S_C1000002_c-1*I_ERI_D2z_S_S_S_C1000002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_C1001002_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S_C1001002_c
   * RHS shell quartet name: SQ_ERI_D_P_S_S_C1001002
   ************************************************************/
  abcd[798] = 2.0E0*I_ERI_D2x_Px_Dxz_S_C1001002_c;
  abcd[799] = 2.0E0*I_ERI_Dxy_Px_Dxz_S_C1001002_c;
  abcd[800] = 2.0E0*I_ERI_Dxz_Px_Dxz_S_C1001002_c;
  abcd[801] = 2.0E0*I_ERI_D2y_Px_Dxz_S_C1001002_c;
  abcd[802] = 2.0E0*I_ERI_Dyz_Px_Dxz_S_C1001002_c;
  abcd[803] = 2.0E0*I_ERI_D2z_Px_Dxz_S_C1001002_c;
  abcd[804] = 2.0E0*I_ERI_D2x_Py_Dxz_S_C1001002_c;
  abcd[805] = 2.0E0*I_ERI_Dxy_Py_Dxz_S_C1001002_c;
  abcd[806] = 2.0E0*I_ERI_Dxz_Py_Dxz_S_C1001002_c;
  abcd[807] = 2.0E0*I_ERI_D2y_Py_Dxz_S_C1001002_c;
  abcd[808] = 2.0E0*I_ERI_Dyz_Py_Dxz_S_C1001002_c;
  abcd[809] = 2.0E0*I_ERI_D2z_Py_Dxz_S_C1001002_c;
  abcd[810] = 2.0E0*I_ERI_D2x_Pz_Dxz_S_C1001002_c;
  abcd[811] = 2.0E0*I_ERI_Dxy_Pz_Dxz_S_C1001002_c;
  abcd[812] = 2.0E0*I_ERI_Dxz_Pz_Dxz_S_C1001002_c;
  abcd[813] = 2.0E0*I_ERI_D2y_Pz_Dxz_S_C1001002_c;
  abcd[814] = 2.0E0*I_ERI_Dyz_Pz_Dxz_S_C1001002_c;
  abcd[815] = 2.0E0*I_ERI_D2z_Pz_Dxz_S_C1001002_c;
  abcd[822] = 2.0E0*I_ERI_D2x_Px_Dyz_S_C1001002_c;
  abcd[823] = 2.0E0*I_ERI_Dxy_Px_Dyz_S_C1001002_c;
  abcd[824] = 2.0E0*I_ERI_Dxz_Px_Dyz_S_C1001002_c;
  abcd[825] = 2.0E0*I_ERI_D2y_Px_Dyz_S_C1001002_c;
  abcd[826] = 2.0E0*I_ERI_Dyz_Px_Dyz_S_C1001002_c;
  abcd[827] = 2.0E0*I_ERI_D2z_Px_Dyz_S_C1001002_c;
  abcd[828] = 2.0E0*I_ERI_D2x_Py_Dyz_S_C1001002_c;
  abcd[829] = 2.0E0*I_ERI_Dxy_Py_Dyz_S_C1001002_c;
  abcd[830] = 2.0E0*I_ERI_Dxz_Py_Dyz_S_C1001002_c;
  abcd[831] = 2.0E0*I_ERI_D2y_Py_Dyz_S_C1001002_c;
  abcd[832] = 2.0E0*I_ERI_Dyz_Py_Dyz_S_C1001002_c;
  abcd[833] = 2.0E0*I_ERI_D2z_Py_Dyz_S_C1001002_c;
  abcd[834] = 2.0E0*I_ERI_D2x_Pz_Dyz_S_C1001002_c;
  abcd[835] = 2.0E0*I_ERI_Dxy_Pz_Dyz_S_C1001002_c;
  abcd[836] = 2.0E0*I_ERI_Dxz_Pz_Dyz_S_C1001002_c;
  abcd[837] = 2.0E0*I_ERI_D2y_Pz_Dyz_S_C1001002_c;
  abcd[838] = 2.0E0*I_ERI_Dyz_Pz_Dyz_S_C1001002_c;
  abcd[839] = 2.0E0*I_ERI_D2z_Pz_Dyz_S_C1001002_c;
  abcd[846] = 2.0E0*I_ERI_D2x_Px_D2z_S_C1001002_c-1*I_ERI_D2x_Px_S_S_C1001002;
  abcd[847] = 2.0E0*I_ERI_Dxy_Px_D2z_S_C1001002_c-1*I_ERI_Dxy_Px_S_S_C1001002;
  abcd[848] = 2.0E0*I_ERI_Dxz_Px_D2z_S_C1001002_c-1*I_ERI_Dxz_Px_S_S_C1001002;
  abcd[849] = 2.0E0*I_ERI_D2y_Px_D2z_S_C1001002_c-1*I_ERI_D2y_Px_S_S_C1001002;
  abcd[850] = 2.0E0*I_ERI_Dyz_Px_D2z_S_C1001002_c-1*I_ERI_Dyz_Px_S_S_C1001002;
  abcd[851] = 2.0E0*I_ERI_D2z_Px_D2z_S_C1001002_c-1*I_ERI_D2z_Px_S_S_C1001002;
  abcd[852] = 2.0E0*I_ERI_D2x_Py_D2z_S_C1001002_c-1*I_ERI_D2x_Py_S_S_C1001002;
  abcd[853] = 2.0E0*I_ERI_Dxy_Py_D2z_S_C1001002_c-1*I_ERI_Dxy_Py_S_S_C1001002;
  abcd[854] = 2.0E0*I_ERI_Dxz_Py_D2z_S_C1001002_c-1*I_ERI_Dxz_Py_S_S_C1001002;
  abcd[855] = 2.0E0*I_ERI_D2y_Py_D2z_S_C1001002_c-1*I_ERI_D2y_Py_S_S_C1001002;
  abcd[856] = 2.0E0*I_ERI_Dyz_Py_D2z_S_C1001002_c-1*I_ERI_Dyz_Py_S_S_C1001002;
  abcd[857] = 2.0E0*I_ERI_D2z_Py_D2z_S_C1001002_c-1*I_ERI_D2z_Py_S_S_C1001002;
  abcd[858] = 2.0E0*I_ERI_D2x_Pz_D2z_S_C1001002_c-1*I_ERI_D2x_Pz_S_S_C1001002;
  abcd[859] = 2.0E0*I_ERI_Dxy_Pz_D2z_S_C1001002_c-1*I_ERI_Dxy_Pz_S_S_C1001002;
  abcd[860] = 2.0E0*I_ERI_Dxz_Pz_D2z_S_C1001002_c-1*I_ERI_Dxz_Pz_S_S_C1001002;
  abcd[861] = 2.0E0*I_ERI_D2y_Pz_D2z_S_C1001002_c-1*I_ERI_D2y_Pz_S_S_C1001002;
  abcd[862] = 2.0E0*I_ERI_Dyz_Pz_D2z_S_C1001002_c-1*I_ERI_Dyz_Pz_S_S_C1001002;
  abcd[863] = 2.0E0*I_ERI_D2z_Pz_D2z_S_C1001002_c-1*I_ERI_D2z_Pz_S_S_C1001002;
}
