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
// BRA1 as redundant position, total RHS integrals evaluated as: 53618
// BRA2 as redundant position, total RHS integrals evaluated as: 55673
// KET1 as redundant position, total RHS integrals evaluated as: 56498
// KET2 as redundant position, total RHS integrals evaluated as: 55641
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

void hgp_os_eri_f_s_d_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
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
  Double I_ERI_F3x_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G4x_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G4x_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G4x_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G4x_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G4x_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G3xy_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G3xz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G2x2y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G2xyz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G2x2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Gx3y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Gx2yz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Gxy2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Gx3z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G4y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G4y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G4y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G4y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G4y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G3yz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G2y2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Gy3z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_G4z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_G4z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_G4z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_G4z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_G4z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_c = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_c = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_c = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_c = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_c = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_c = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_c = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_c = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_c = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_c = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_c = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_c = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_c = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_c = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_c = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_c = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_c = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_c = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_c = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_c = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_c = 0.0E0;
  Double I_ERI_F3x_S_S_S = 0.0E0;
  Double I_ERI_F2xy_S_S_S = 0.0E0;
  Double I_ERI_F2xz_S_S_S = 0.0E0;
  Double I_ERI_Fx2y_S_S_S = 0.0E0;
  Double I_ERI_Fxyz_S_S_S = 0.0E0;
  Double I_ERI_Fx2z_S_S_S = 0.0E0;
  Double I_ERI_F3y_S_S_S = 0.0E0;
  Double I_ERI_F2yz_S_S_S = 0.0E0;
  Double I_ERI_Fy2z_S_S_S = 0.0E0;
  Double I_ERI_F3z_S_S_S = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_d = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_d = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_d = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_d = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_d = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_d = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_d = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_d = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_d = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_d = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_d = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_d = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_d = 0.0E0;
  Double I_ERI_G4x_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_bc = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Px_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_b = 0.0E0;
  Double I_ERI_G4y_S_Px_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_b = 0.0E0;
  Double I_ERI_G4z_S_Px_S_b = 0.0E0;
  Double I_ERI_G4x_S_Py_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_b = 0.0E0;
  Double I_ERI_G4y_S_Py_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_b = 0.0E0;
  Double I_ERI_G4z_S_Py_S_b = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_b = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_b = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_b = 0.0E0;
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
  Double I_ERI_H5x_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H5y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H5z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_H5x_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H5y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H5z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_H5x_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H5y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_H5z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_bb = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_bb = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_bb = 0.0E0;
  Double I_ERI_F3x_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G4x_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G4x_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G4x_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G4x_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G4x_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G3xy_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G3xz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G2x2y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G2xyz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G2x2z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Gx3y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Gx2yz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Gxy2z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Gx3z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G4y_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G4y_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G4y_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G4y_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G4y_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G3yz_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G2y2z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_Gy3z_S_cd = 0.0E0;
  Double I_ERI_F3x_S_G4z_S_cd = 0.0E0;
  Double I_ERI_F2xy_S_G4z_S_cd = 0.0E0;
  Double I_ERI_F2xz_S_G4z_S_cd = 0.0E0;
  Double I_ERI_Fx2y_S_G4z_S_cd = 0.0E0;
  Double I_ERI_Fxyz_S_G4z_S_cd = 0.0E0;
  Double I_ERI_Fx2z_S_G4z_S_cd = 0.0E0;
  Double I_ERI_F3y_S_G4z_S_cd = 0.0E0;
  Double I_ERI_F2yz_S_G4z_S_cd = 0.0E0;
  Double I_ERI_Fy2z_S_G4z_S_cd = 0.0E0;
  Double I_ERI_F3z_S_G4z_S_cd = 0.0E0;
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
  Double I_ERI_G4x_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_bd = 0.0E0;
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
  Double I_ERI_F3x_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_bd = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_bd = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_bd = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_bd = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_bd = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_bd = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_bd = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_bd = 0.0E0;
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
  Double I_ERI_F3x_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G4x_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G4x_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G4x_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G4x_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G4x_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G3xy_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G3xz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G2x2y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G2xyz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G2x2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Gx3y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Gx2yz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Gxy2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Gx3z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G4y_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G4y_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G4y_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G4y_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G4y_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G3yz_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G2y2z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_Gy3z_S_dd = 0.0E0;
  Double I_ERI_F3x_S_G4z_S_dd = 0.0E0;
  Double I_ERI_F2xy_S_G4z_S_dd = 0.0E0;
  Double I_ERI_F2xz_S_G4z_S_dd = 0.0E0;
  Double I_ERI_Fx2y_S_G4z_S_dd = 0.0E0;
  Double I_ERI_Fxyz_S_G4z_S_dd = 0.0E0;
  Double I_ERI_Fx2z_S_G4z_S_dd = 0.0E0;
  Double I_ERI_F3y_S_G4z_S_dd = 0.0E0;
  Double I_ERI_F2yz_S_G4z_S_dd = 0.0E0;
  Double I_ERI_Fy2z_S_G4z_S_dd = 0.0E0;
  Double I_ERI_F3z_S_G4z_S_dd = 0.0E0;
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
      Double I_ERI_S_S_S_S_M7_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ERI_S_S_S_S_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M1_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M2_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M3_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M4_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M5_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M6_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M7_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER49;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = 1.0E0+u2*ONEOVER17*I_ERI_S_S_S_S_M7_vrr;
        I_ERI_S_S_S_S_M7_vrr = ONEOVER15*I_ERI_S_S_S_S_M7_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M7_vrr  = f*I_ERI_S_S_S_S_M7_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M6_vrr  = ONEOVER13*(u2*I_ERI_S_S_S_S_M7_vrr+f);
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
        I_ERI_S_S_S_S_M7_vrr_d = oneO2u*(13.0E0*I_ERI_S_S_S_S_M6_vrr_d-f);

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);
        I_ERI_S_S_S_S_M5_vrr = static_cast<Double>(I_ERI_S_S_S_S_M5_vrr_d);
        I_ERI_S_S_S_S_M6_vrr = static_cast<Double>(I_ERI_S_S_S_S_M6_vrr_d);
        I_ERI_S_S_S_S_M7_vrr = static_cast<Double>(I_ERI_S_S_S_S_M7_vrr_d);

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
        I_ERI_S_S_S_S_M7_vrr = oneO2u*(13.0E0*I_ERI_S_S_S_S_M6_vrr-f);

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
        Double erfPref_15 = erfPref_13*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_Px_S_S_S_M6_vrr = PAX*I_ERI_S_S_S_S_M6_vrr+WPX*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Py_S_S_S_M6_vrr = PAY*I_ERI_S_S_S_S_M6_vrr+WPY*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Pz_S_S_S_M6_vrr = PAZ*I_ERI_S_S_S_S_M6_vrr+WPZ*I_ERI_S_S_S_S_M7_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M5_vrr = PAX*I_ERI_Px_S_S_S_M5_vrr+WPX*I_ERI_Px_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Dxy_S_S_S_M5_vrr = PAY*I_ERI_Px_S_S_S_M5_vrr+WPY*I_ERI_Px_S_S_S_M6_vrr;
      Double I_ERI_D2y_S_S_S_M5_vrr = PAY*I_ERI_Py_S_S_S_M5_vrr+WPY*I_ERI_Py_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_D2z_S_S_S_M5_vrr = PAZ*I_ERI_Pz_S_S_S_M5_vrr+WPZ*I_ERI_Pz_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M4_vrr = PAX*I_ERI_Px_S_S_S_M4_vrr+WPX*I_ERI_Px_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Dxy_S_S_S_M4_vrr = PAY*I_ERI_Px_S_S_S_M4_vrr+WPY*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_Dxz_S_S_S_M4_vrr = PAZ*I_ERI_Px_S_S_S_M4_vrr+WPZ*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_D2y_S_S_S_M4_vrr = PAY*I_ERI_Py_S_S_S_M4_vrr+WPY*I_ERI_Py_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Dyz_S_S_S_M4_vrr = PAZ*I_ERI_Py_S_S_S_M4_vrr+WPZ*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_D2z_S_S_S_M4_vrr = PAZ*I_ERI_Pz_S_S_S_M4_vrr+WPZ*I_ERI_Pz_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F2xz_S_S_S_M4_vrr = PAZ*I_ERI_D2x_S_S_S_M4_vrr+WPZ*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_Fx2y_S_S_S_M4_vrr = PAX*I_ERI_D2y_S_S_S_M4_vrr+WPX*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_Fxyz_S_S_S_M4_vrr = PAZ*I_ERI_Dxy_S_S_S_M4_vrr+WPZ*I_ERI_Dxy_S_S_S_M5_vrr;
      Double I_ERI_Fx2z_S_S_S_M4_vrr = PAX*I_ERI_D2z_S_S_S_M4_vrr+WPX*I_ERI_D2z_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F2yz_S_S_S_M4_vrr = PAZ*I_ERI_D2y_S_S_S_M4_vrr+WPZ*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_Fy2z_S_S_S_M4_vrr = PAY*I_ERI_D2z_S_S_S_M4_vrr+WPY*I_ERI_D2z_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M3_vrr = PAX*I_ERI_S_S_Px_S_M3_vrr+WPX*I_ERI_S_S_Px_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Py_S_Px_S_M3_vrr = PAY*I_ERI_S_S_Px_S_M3_vrr+WPY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Px_S_M3_vrr = PAZ*I_ERI_S_S_Px_S_M3_vrr+WPZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_Py_S_M3_vrr = PAX*I_ERI_S_S_Py_S_M3_vrr+WPX*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Py_S_M3_vrr = PAY*I_ERI_S_S_Py_S_M3_vrr+WPY*I_ERI_S_S_Py_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Pz_S_Py_S_M3_vrr = PAZ*I_ERI_S_S_Py_S_M3_vrr+WPZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_Pz_S_M3_vrr = PAX*I_ERI_S_S_Pz_S_M3_vrr+WPX*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Pz_S_M3_vrr = PAY*I_ERI_S_S_Pz_S_M3_vrr+WPY*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Pz_S_M3_vrr = PAZ*I_ERI_S_S_Pz_S_M3_vrr+WPZ*I_ERI_S_S_Pz_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M3_vrr = QCX*I_ERI_D2x_S_S_S_M3_vrr+WQX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_Px_S_M3_vrr = QCX*I_ERI_Dxy_S_S_S_M3_vrr+WQX*I_ERI_Dxy_S_S_S_M4_vrr+oned2k*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_Dxz_S_Px_S_M3_vrr = QCX*I_ERI_Dxz_S_S_S_M3_vrr+WQX*I_ERI_Dxz_S_S_S_M4_vrr+oned2k*I_ERI_Pz_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_Px_S_M3_vrr = QCX*I_ERI_D2y_S_S_S_M3_vrr+WQX*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Dyz_S_Px_S_M3_vrr = QCX*I_ERI_Dyz_S_S_S_M3_vrr+WQX*I_ERI_Dyz_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_Px_S_M3_vrr = QCX*I_ERI_D2z_S_S_S_M3_vrr+WQX*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_D2x_S_Py_S_M3_vrr = QCY*I_ERI_D2x_S_S_S_M3_vrr+WQY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_Py_S_M3_vrr = QCY*I_ERI_Dxy_S_S_S_M3_vrr+WQY*I_ERI_Dxy_S_S_S_M4_vrr+oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_Dxz_S_Py_S_M3_vrr = QCY*I_ERI_Dxz_S_S_S_M3_vrr+WQY*I_ERI_Dxz_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_Py_S_M3_vrr = QCY*I_ERI_D2y_S_S_S_M3_vrr+WQY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2k*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_Dyz_S_Py_S_M3_vrr = QCY*I_ERI_Dyz_S_S_S_M3_vrr+WQY*I_ERI_Dyz_S_S_S_M4_vrr+oned2k*I_ERI_Pz_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_Py_S_M3_vrr = QCY*I_ERI_D2z_S_S_S_M3_vrr+WQY*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_D2x_S_Pz_S_M3_vrr = QCZ*I_ERI_D2x_S_S_S_M3_vrr+WQZ*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_Pz_S_M3_vrr = QCZ*I_ERI_Dxy_S_S_S_M3_vrr+WQZ*I_ERI_Dxy_S_S_S_M4_vrr;
      Double I_ERI_Dxz_S_Pz_S_M3_vrr = QCZ*I_ERI_Dxz_S_S_S_M3_vrr+WQZ*I_ERI_Dxz_S_S_S_M4_vrr+oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_Pz_S_M3_vrr = QCZ*I_ERI_D2y_S_S_S_M3_vrr+WQZ*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Dyz_S_Pz_S_M3_vrr = QCZ*I_ERI_Dyz_S_S_S_M3_vrr+WQZ*I_ERI_Dyz_S_S_S_M4_vrr+oned2k*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_Pz_S_M3_vrr = QCZ*I_ERI_D2z_S_S_S_M3_vrr+WQZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2k*I_ERI_Pz_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_F_S_P_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M3_vrr = QCX*I_ERI_F3x_S_S_S_M3_vrr+WQX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2k*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_Px_S_M3_vrr = QCX*I_ERI_F2xy_S_S_S_M3_vrr+WQX*I_ERI_F2xy_S_S_S_M4_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_Px_S_M3_vrr = QCX*I_ERI_F2xz_S_S_S_M3_vrr+WQX*I_ERI_F2xz_S_S_S_M4_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_Px_S_M3_vrr = QCX*I_ERI_Fx2y_S_S_S_M3_vrr+WQX*I_ERI_Fx2y_S_S_S_M4_vrr+oned2k*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fxyz_S_Px_S_M3_vrr = QCX*I_ERI_Fxyz_S_S_S_M3_vrr+WQX*I_ERI_Fxyz_S_S_S_M4_vrr+oned2k*I_ERI_Dyz_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_Px_S_M3_vrr = QCX*I_ERI_Fx2z_S_S_S_M3_vrr+WQX*I_ERI_Fx2z_S_S_S_M4_vrr+oned2k*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_Px_S_M3_vrr = QCX*I_ERI_F3y_S_S_S_M3_vrr+WQX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_Px_S_M3_vrr = QCX*I_ERI_F2yz_S_S_S_M3_vrr+WQX*I_ERI_F2yz_S_S_S_M4_vrr;
      Double I_ERI_Fy2z_S_Px_S_M3_vrr = QCX*I_ERI_Fy2z_S_S_S_M3_vrr+WQX*I_ERI_Fy2z_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_Px_S_M3_vrr = QCX*I_ERI_F3z_S_S_S_M3_vrr+WQX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_F3x_S_Py_S_M3_vrr = QCY*I_ERI_F3x_S_S_S_M3_vrr+WQY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_Py_S_M3_vrr = QCY*I_ERI_F2xy_S_S_S_M3_vrr+WQY*I_ERI_F2xy_S_S_S_M4_vrr+oned2k*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_Py_S_M3_vrr = QCY*I_ERI_F2xz_S_S_S_M3_vrr+WQY*I_ERI_F2xz_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_Py_S_M3_vrr = QCY*I_ERI_Fx2y_S_S_S_M3_vrr+WQY*I_ERI_Fx2y_S_S_S_M4_vrr+2*oned2k*I_ERI_Dxy_S_S_S_M4_vrr;
      Double I_ERI_Fxyz_S_Py_S_M3_vrr = QCY*I_ERI_Fxyz_S_S_S_M3_vrr+WQY*I_ERI_Fxyz_S_S_S_M4_vrr+oned2k*I_ERI_Dxz_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_Py_S_M3_vrr = QCY*I_ERI_Fx2z_S_S_S_M3_vrr+WQY*I_ERI_Fx2z_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_Py_S_M3_vrr = QCY*I_ERI_F3y_S_S_S_M3_vrr+WQY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2k*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_Py_S_M3_vrr = QCY*I_ERI_F2yz_S_S_S_M3_vrr+WQY*I_ERI_F2yz_S_S_S_M4_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M4_vrr;
      Double I_ERI_Fy2z_S_Py_S_M3_vrr = QCY*I_ERI_Fy2z_S_S_S_M3_vrr+WQY*I_ERI_Fy2z_S_S_S_M4_vrr+oned2k*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_Py_S_M3_vrr = QCY*I_ERI_F3z_S_S_S_M3_vrr+WQY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_F3x_S_Pz_S_M3_vrr = QCZ*I_ERI_F3x_S_S_S_M3_vrr+WQZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_Pz_S_M3_vrr = QCZ*I_ERI_F2xy_S_S_S_M3_vrr+WQZ*I_ERI_F2xy_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_Pz_S_M3_vrr = QCZ*I_ERI_F2xz_S_S_S_M3_vrr+WQZ*I_ERI_F2xz_S_S_S_M4_vrr+oned2k*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M3_vrr = QCZ*I_ERI_Fx2y_S_S_S_M3_vrr+WQZ*I_ERI_Fx2y_S_S_S_M4_vrr;
      Double I_ERI_Fxyz_S_Pz_S_M3_vrr = QCZ*I_ERI_Fxyz_S_S_S_M3_vrr+WQZ*I_ERI_Fxyz_S_S_S_M4_vrr+oned2k*I_ERI_Dxy_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M3_vrr = QCZ*I_ERI_Fx2z_S_S_S_M3_vrr+WQZ*I_ERI_Fx2z_S_S_S_M4_vrr+2*oned2k*I_ERI_Dxz_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_Pz_S_M3_vrr = QCZ*I_ERI_F3y_S_S_S_M3_vrr+WQZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_Pz_S_M3_vrr = QCZ*I_ERI_F2yz_S_S_S_M3_vrr+WQZ*I_ERI_F2yz_S_S_S_M4_vrr+oned2k*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fy2z_S_Pz_S_M3_vrr = QCZ*I_ERI_Fy2z_S_S_S_M3_vrr+WQZ*I_ERI_Fy2z_S_S_S_M4_vrr+2*oned2k*I_ERI_Dyz_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_Pz_S_M3_vrr = QCZ*I_ERI_F3z_S_S_S_M3_vrr+WQZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2k*I_ERI_D2z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_P_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M2_vrr = QCX*I_ERI_Px_S_Px_S_M2_vrr+WQX*I_ERI_Px_S_Px_S_M3_vrr+oned2e*I_ERI_Px_S_S_S_M2_vrr-rhod2esq*I_ERI_Px_S_S_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Py_S_D2x_S_M2_vrr = QCX*I_ERI_Py_S_Px_S_M2_vrr+WQX*I_ERI_Py_S_Px_S_M3_vrr+oned2e*I_ERI_Py_S_S_S_M2_vrr-rhod2esq*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_Pz_S_D2x_S_M2_vrr = QCX*I_ERI_Pz_S_Px_S_M2_vrr+WQX*I_ERI_Pz_S_Px_S_M3_vrr+oned2e*I_ERI_Pz_S_S_S_M2_vrr-rhod2esq*I_ERI_Pz_S_S_S_M3_vrr;
      Double I_ERI_Px_S_D2y_S_M2_vrr = QCY*I_ERI_Px_S_Py_S_M2_vrr+WQY*I_ERI_Px_S_Py_S_M3_vrr+oned2e*I_ERI_Px_S_S_S_M2_vrr-rhod2esq*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Py_S_D2y_S_M2_vrr = QCY*I_ERI_Py_S_Py_S_M2_vrr+WQY*I_ERI_Py_S_Py_S_M3_vrr+oned2e*I_ERI_Py_S_S_S_M2_vrr-rhod2esq*I_ERI_Py_S_S_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Pz_S_D2y_S_M2_vrr = QCY*I_ERI_Pz_S_Py_S_M2_vrr+WQY*I_ERI_Pz_S_Py_S_M3_vrr+oned2e*I_ERI_Pz_S_S_S_M2_vrr-rhod2esq*I_ERI_Pz_S_S_S_M3_vrr;
      Double I_ERI_Px_S_D2z_S_M2_vrr = QCZ*I_ERI_Px_S_Pz_S_M2_vrr+WQZ*I_ERI_Px_S_Pz_S_M3_vrr+oned2e*I_ERI_Px_S_S_S_M2_vrr-rhod2esq*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_Py_S_D2z_S_M2_vrr = QCZ*I_ERI_Py_S_Pz_S_M2_vrr+WQZ*I_ERI_Py_S_Pz_S_M3_vrr+oned2e*I_ERI_Py_S_S_S_M2_vrr-rhod2esq*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_Pz_S_D2z_S_M2_vrr = QCZ*I_ERI_Pz_S_Pz_S_M2_vrr+WQZ*I_ERI_Pz_S_Pz_S_M3_vrr+oned2e*I_ERI_Pz_S_S_S_M2_vrr-rhod2esq*I_ERI_Pz_S_S_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M2_vrr = QCX*I_ERI_D2x_S_Px_S_M2_vrr+WQX*I_ERI_D2x_S_Px_S_M3_vrr+oned2e*I_ERI_D2x_S_S_S_M2_vrr-rhod2esq*I_ERI_D2x_S_S_S_M3_vrr+2*oned2k*I_ERI_Px_S_Px_S_M3_vrr;
      Double I_ERI_Dxy_S_D2x_S_M2_vrr = QCX*I_ERI_Dxy_S_Px_S_M2_vrr+WQX*I_ERI_Dxy_S_Px_S_M3_vrr+oned2e*I_ERI_Dxy_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M3_vrr+oned2k*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_Dxz_S_D2x_S_M2_vrr = QCX*I_ERI_Dxz_S_Px_S_M2_vrr+WQX*I_ERI_Dxz_S_Px_S_M3_vrr+oned2e*I_ERI_Dxz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M3_vrr+oned2k*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_D2y_S_D2x_S_M2_vrr = QCX*I_ERI_D2y_S_Px_S_M2_vrr+WQX*I_ERI_D2y_S_Px_S_M3_vrr+oned2e*I_ERI_D2y_S_S_S_M2_vrr-rhod2esq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_D2x_S_M2_vrr = QCX*I_ERI_Dyz_S_Px_S_M2_vrr+WQX*I_ERI_Dyz_S_Px_S_M3_vrr+oned2e*I_ERI_Dyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_D2x_S_M2_vrr = QCX*I_ERI_D2z_S_Px_S_M2_vrr+WQX*I_ERI_D2z_S_Px_S_M3_vrr+oned2e*I_ERI_D2z_S_S_S_M2_vrr-rhod2esq*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_D2x_S_Dxy_S_M2_vrr = QCY*I_ERI_D2x_S_Px_S_M2_vrr+WQY*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_D2y_S_Dxy_S_M2_vrr = QCY*I_ERI_D2y_S_Px_S_M2_vrr+WQY*I_ERI_D2y_S_Px_S_M3_vrr+2*oned2k*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Dxy_S_M2_vrr = QCY*I_ERI_D2z_S_Px_S_M2_vrr+WQY*I_ERI_D2z_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_Dxz_S_M2_vrr = QCZ*I_ERI_D2x_S_Px_S_M2_vrr+WQZ*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_D2y_S_Dxz_S_M2_vrr = QCZ*I_ERI_D2y_S_Px_S_M2_vrr+WQZ*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Dxz_S_M2_vrr = QCZ*I_ERI_D2z_S_Px_S_M2_vrr+WQZ*I_ERI_D2z_S_Px_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_D2y_S_M2_vrr = QCY*I_ERI_D2x_S_Py_S_M2_vrr+WQY*I_ERI_D2x_S_Py_S_M3_vrr+oned2e*I_ERI_D2x_S_S_S_M2_vrr-rhod2esq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_D2y_S_M2_vrr = QCY*I_ERI_Dxy_S_Py_S_M2_vrr+WQY*I_ERI_Dxy_S_Py_S_M3_vrr+oned2e*I_ERI_Dxy_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M3_vrr+oned2k*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_Dxz_S_D2y_S_M2_vrr = QCY*I_ERI_Dxz_S_Py_S_M2_vrr+WQY*I_ERI_Dxz_S_Py_S_M3_vrr+oned2e*I_ERI_Dxz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_D2y_S_M2_vrr = QCY*I_ERI_D2y_S_Py_S_M2_vrr+WQY*I_ERI_D2y_S_Py_S_M3_vrr+oned2e*I_ERI_D2y_S_S_S_M2_vrr-rhod2esq*I_ERI_D2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Py_S_Py_S_M3_vrr;
      Double I_ERI_Dyz_S_D2y_S_M2_vrr = QCY*I_ERI_Dyz_S_Py_S_M2_vrr+WQY*I_ERI_Dyz_S_Py_S_M3_vrr+oned2e*I_ERI_Dyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M3_vrr+oned2k*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_D2z_S_D2y_S_M2_vrr = QCY*I_ERI_D2z_S_Py_S_M2_vrr+WQY*I_ERI_D2z_S_Py_S_M3_vrr+oned2e*I_ERI_D2z_S_S_S_M2_vrr-rhod2esq*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_D2x_S_Dyz_S_M2_vrr = QCZ*I_ERI_D2x_S_Py_S_M2_vrr+WQZ*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Dyz_S_M2_vrr = QCZ*I_ERI_D2y_S_Py_S_M2_vrr+WQZ*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_D2z_S_Dyz_S_M2_vrr = QCZ*I_ERI_D2z_S_Py_S_M2_vrr+WQZ*I_ERI_D2z_S_Py_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_D2z_S_M2_vrr = QCZ*I_ERI_D2x_S_Pz_S_M2_vrr+WQZ*I_ERI_D2x_S_Pz_S_M3_vrr+oned2e*I_ERI_D2x_S_S_S_M2_vrr-rhod2esq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_D2z_S_M2_vrr = QCZ*I_ERI_Dxy_S_Pz_S_M2_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M3_vrr+oned2e*I_ERI_Dxy_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Dxz_S_D2z_S_M2_vrr = QCZ*I_ERI_Dxz_S_Pz_S_M2_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M3_vrr+oned2e*I_ERI_Dxz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M3_vrr+oned2k*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_D2y_S_D2z_S_M2_vrr = QCZ*I_ERI_D2y_S_Pz_S_M2_vrr+WQZ*I_ERI_D2y_S_Pz_S_M3_vrr+oned2e*I_ERI_D2y_S_S_S_M2_vrr-rhod2esq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Dyz_S_D2z_S_M2_vrr = QCZ*I_ERI_Dyz_S_Pz_S_M2_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M3_vrr+oned2e*I_ERI_Dyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M3_vrr+oned2k*I_ERI_Py_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_D2z_S_M2_vrr = QCZ*I_ERI_D2z_S_Pz_S_M2_vrr+WQZ*I_ERI_D2z_S_Pz_S_M3_vrr+oned2e*I_ERI_D2z_S_S_S_M2_vrr-rhod2esq*I_ERI_D2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M2_vrr = PAX*I_ERI_D2x_S_Px_S_M2_vrr+WPX*I_ERI_D2x_S_Px_S_M3_vrr+2*oned2z*I_ERI_Px_S_Px_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Px_S_M3_vrr+oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_Px_S_M2_vrr = PAY*I_ERI_D2x_S_Px_S_M2_vrr+WPY*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F2xz_S_Px_S_M2_vrr = PAZ*I_ERI_D2x_S_Px_S_M2_vrr+WPZ*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_Fx2y_S_Px_S_M2_vrr = PAX*I_ERI_D2y_S_Px_S_M2_vrr+WPX*I_ERI_D2y_S_Px_S_M3_vrr+oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_Px_S_M2_vrr = PAZ*I_ERI_Dxy_S_Px_S_M2_vrr+WPZ*I_ERI_Dxy_S_Px_S_M3_vrr;
      Double I_ERI_Fx2z_S_Px_S_M2_vrr = PAX*I_ERI_D2z_S_Px_S_M2_vrr+WPX*I_ERI_D2z_S_Px_S_M3_vrr+oned2k*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_Px_S_M2_vrr = PAY*I_ERI_D2y_S_Px_S_M2_vrr+WPY*I_ERI_D2y_S_Px_S_M3_vrr+2*oned2z*I_ERI_Py_S_Px_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_F2yz_S_Px_S_M2_vrr = PAZ*I_ERI_D2y_S_Px_S_M2_vrr+WPZ*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_Fy2z_S_Px_S_M2_vrr = PAY*I_ERI_D2z_S_Px_S_M2_vrr+WPY*I_ERI_D2z_S_Px_S_M3_vrr;
      Double I_ERI_F3z_S_Px_S_M2_vrr = PAZ*I_ERI_D2z_S_Px_S_M2_vrr+WPZ*I_ERI_D2z_S_Px_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Px_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_F3x_S_Py_S_M2_vrr = PAX*I_ERI_D2x_S_Py_S_M2_vrr+WPX*I_ERI_D2x_S_Py_S_M3_vrr+2*oned2z*I_ERI_Px_S_Py_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_F2xy_S_Py_S_M2_vrr = PAY*I_ERI_D2x_S_Py_S_M2_vrr+WPY*I_ERI_D2x_S_Py_S_M3_vrr+oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_Py_S_M2_vrr = PAZ*I_ERI_D2x_S_Py_S_M2_vrr+WPZ*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_Fx2y_S_Py_S_M2_vrr = PAX*I_ERI_D2y_S_Py_S_M2_vrr+WPX*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_Fxyz_S_Py_S_M2_vrr = PAZ*I_ERI_Dxy_S_Py_S_M2_vrr+WPZ*I_ERI_Dxy_S_Py_S_M3_vrr;
      Double I_ERI_Fx2z_S_Py_S_M2_vrr = PAX*I_ERI_D2z_S_Py_S_M2_vrr+WPX*I_ERI_D2z_S_Py_S_M3_vrr;
      Double I_ERI_F3y_S_Py_S_M2_vrr = PAY*I_ERI_D2y_S_Py_S_M2_vrr+WPY*I_ERI_D2y_S_Py_S_M3_vrr+2*oned2z*I_ERI_Py_S_Py_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M3_vrr+oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_Py_S_M2_vrr = PAZ*I_ERI_D2y_S_Py_S_M2_vrr+WPZ*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_Fy2z_S_Py_S_M2_vrr = PAY*I_ERI_D2z_S_Py_S_M2_vrr+WPY*I_ERI_D2z_S_Py_S_M3_vrr+oned2k*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_Py_S_M2_vrr = PAZ*I_ERI_D2z_S_Py_S_M2_vrr+WPZ*I_ERI_D2z_S_Py_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Py_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_F3x_S_Pz_S_M2_vrr = PAX*I_ERI_D2x_S_Pz_S_M2_vrr+WPX*I_ERI_D2x_S_Pz_S_M3_vrr+2*oned2z*I_ERI_Px_S_Pz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_F2xy_S_Pz_S_M2_vrr = PAY*I_ERI_D2x_S_Pz_S_M2_vrr+WPY*I_ERI_D2x_S_Pz_S_M3_vrr;
      Double I_ERI_F2xz_S_Pz_S_M2_vrr = PAZ*I_ERI_D2x_S_Pz_S_M2_vrr+WPZ*I_ERI_D2x_S_Pz_S_M3_vrr+oned2k*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_Pz_S_M2_vrr = PAX*I_ERI_D2y_S_Pz_S_M2_vrr+WPX*I_ERI_D2y_S_Pz_S_M3_vrr;
      Double I_ERI_Fxyz_S_Pz_S_M2_vrr = PAZ*I_ERI_Dxy_S_Pz_S_M2_vrr+WPZ*I_ERI_Dxy_S_Pz_S_M3_vrr+oned2k*I_ERI_Dxy_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_Pz_S_M2_vrr = PAX*I_ERI_D2z_S_Pz_S_M2_vrr+WPX*I_ERI_D2z_S_Pz_S_M3_vrr;
      Double I_ERI_F3y_S_Pz_S_M2_vrr = PAY*I_ERI_D2y_S_Pz_S_M2_vrr+WPY*I_ERI_D2y_S_Pz_S_M3_vrr+2*oned2z*I_ERI_Py_S_Pz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M3_vrr;
      Double I_ERI_F2yz_S_Pz_S_M2_vrr = PAZ*I_ERI_D2y_S_Pz_S_M2_vrr+WPZ*I_ERI_D2y_S_Pz_S_M3_vrr+oned2k*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_Pz_S_M2_vrr = PAY*I_ERI_D2z_S_Pz_S_M2_vrr+WPY*I_ERI_D2z_S_Pz_S_M3_vrr;
      Double I_ERI_F3z_S_Pz_S_M2_vrr = PAZ*I_ERI_D2z_S_Pz_S_M2_vrr+WPZ*I_ERI_D2z_S_Pz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M3_vrr+oned2k*I_ERI_D2z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 10 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M2_vrr = QCX*I_ERI_F3x_S_Px_S_M2_vrr+WQX*I_ERI_F3x_S_Px_S_M3_vrr+oned2e*I_ERI_F3x_S_S_S_M2_vrr-rhod2esq*I_ERI_F3x_S_S_S_M3_vrr+3*oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F2xy_S_D2x_S_M2_vrr = QCX*I_ERI_F2xy_S_Px_S_M2_vrr+WQX*I_ERI_F2xy_S_Px_S_M3_vrr+oned2e*I_ERI_F2xy_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M3_vrr;
      Double I_ERI_F2xz_S_D2x_S_M2_vrr = QCX*I_ERI_F2xz_S_Px_S_M2_vrr+WQX*I_ERI_F2xz_S_Px_S_M3_vrr+oned2e*I_ERI_F2xz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_Px_S_M3_vrr;
      Double I_ERI_Fx2y_S_D2x_S_M2_vrr = QCX*I_ERI_Fx2y_S_Px_S_M2_vrr+WQX*I_ERI_Fx2y_S_Px_S_M3_vrr+oned2e*I_ERI_Fx2y_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M3_vrr+oned2k*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_Fxyz_S_D2x_S_M2_vrr = QCX*I_ERI_Fxyz_S_Px_S_M2_vrr+WQX*I_ERI_Fxyz_S_Px_S_M3_vrr+oned2e*I_ERI_Fxyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dyz_S_Px_S_M3_vrr;
      Double I_ERI_Fx2z_S_D2x_S_M2_vrr = QCX*I_ERI_Fx2z_S_Px_S_M2_vrr+WQX*I_ERI_Fx2z_S_Px_S_M3_vrr+oned2e*I_ERI_Fx2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M3_vrr+oned2k*I_ERI_D2z_S_Px_S_M3_vrr;
      Double I_ERI_F3y_S_D2x_S_M2_vrr = QCX*I_ERI_F3y_S_Px_S_M2_vrr+WQX*I_ERI_F3y_S_Px_S_M3_vrr+oned2e*I_ERI_F3y_S_S_S_M2_vrr-rhod2esq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_D2x_S_M2_vrr = QCX*I_ERI_F2yz_S_Px_S_M2_vrr+WQX*I_ERI_F2yz_S_Px_S_M3_vrr+oned2e*I_ERI_F2yz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Fy2z_S_D2x_S_M2_vrr = QCX*I_ERI_Fy2z_S_Px_S_M2_vrr+WQX*I_ERI_Fy2z_S_Px_S_M3_vrr+oned2e*I_ERI_Fy2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_D2x_S_M2_vrr = QCX*I_ERI_F3z_S_Px_S_M2_vrr+WQX*I_ERI_F3z_S_Px_S_M3_vrr+oned2e*I_ERI_F3z_S_S_S_M2_vrr-rhod2esq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_F3x_S_Dxy_S_M2_vrr = QCY*I_ERI_F3x_S_Px_S_M2_vrr+WQY*I_ERI_F3x_S_Px_S_M3_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M2_vrr = QCY*I_ERI_F2xy_S_Px_S_M2_vrr+WQY*I_ERI_F2xy_S_Px_S_M3_vrr+oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F2xz_S_Dxy_S_M2_vrr = QCY*I_ERI_F2xz_S_Px_S_M2_vrr+WQY*I_ERI_F2xz_S_Px_S_M3_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_M2_vrr = QCY*I_ERI_Fx2y_S_Px_S_M2_vrr+WQY*I_ERI_Fx2y_S_Px_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_Px_S_M3_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_M2_vrr = QCY*I_ERI_Fx2z_S_Px_S_M2_vrr+WQY*I_ERI_Fx2z_S_Px_S_M3_vrr;
      Double I_ERI_F3y_S_Dxy_S_M2_vrr = QCY*I_ERI_F3y_S_Px_S_M2_vrr+WQY*I_ERI_F3y_S_Px_S_M3_vrr+3*oned2k*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_F2yz_S_Dxy_S_M2_vrr = QCY*I_ERI_F2yz_S_Px_S_M2_vrr+WQY*I_ERI_F2yz_S_Px_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_Px_S_M3_vrr;
      Double I_ERI_F3z_S_Dxy_S_M2_vrr = QCY*I_ERI_F3z_S_Px_S_M2_vrr+WQY*I_ERI_F3z_S_Px_S_M3_vrr;
      Double I_ERI_F3x_S_Dxz_S_M2_vrr = QCZ*I_ERI_F3x_S_Px_S_M2_vrr+WQZ*I_ERI_F3x_S_Px_S_M3_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M2_vrr = QCZ*I_ERI_F2xy_S_Px_S_M2_vrr+WQZ*I_ERI_F2xy_S_Px_S_M3_vrr;
      Double I_ERI_F2xz_S_Dxz_S_M2_vrr = QCZ*I_ERI_F2xz_S_Px_S_M2_vrr+WQZ*I_ERI_F2xz_S_Px_S_M3_vrr+oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F3y_S_Dxz_S_M2_vrr = QCZ*I_ERI_F3y_S_Px_S_M2_vrr+WQZ*I_ERI_F3y_S_Px_S_M3_vrr;
      Double I_ERI_F2yz_S_Dxz_S_M2_vrr = QCZ*I_ERI_F2yz_S_Px_S_M2_vrr+WQZ*I_ERI_F2yz_S_Px_S_M3_vrr+oned2k*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_F3z_S_Dxz_S_M2_vrr = QCZ*I_ERI_F3z_S_Px_S_M2_vrr+WQZ*I_ERI_F3z_S_Px_S_M3_vrr+3*oned2k*I_ERI_D2z_S_Px_S_M3_vrr;
      Double I_ERI_F3x_S_D2y_S_M2_vrr = QCY*I_ERI_F3x_S_Py_S_M2_vrr+WQY*I_ERI_F3x_S_Py_S_M3_vrr+oned2e*I_ERI_F3x_S_S_S_M2_vrr-rhod2esq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_D2y_S_M2_vrr = QCY*I_ERI_F2xy_S_Py_S_M2_vrr+WQY*I_ERI_F2xy_S_Py_S_M3_vrr+oned2e*I_ERI_F2xy_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M3_vrr+oned2k*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_F2xz_S_D2y_S_M2_vrr = QCY*I_ERI_F2xz_S_Py_S_M2_vrr+WQY*I_ERI_F2xz_S_Py_S_M3_vrr+oned2e*I_ERI_F2xz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_D2y_S_M2_vrr = QCY*I_ERI_Fx2y_S_Py_S_M2_vrr+WQY*I_ERI_Fx2y_S_Py_S_M3_vrr+oned2e*I_ERI_Fx2y_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_Py_S_M3_vrr;
      Double I_ERI_Fxyz_S_D2y_S_M2_vrr = QCY*I_ERI_Fxyz_S_Py_S_M2_vrr+WQY*I_ERI_Fxyz_S_Py_S_M3_vrr+oned2e*I_ERI_Fxyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dxz_S_Py_S_M3_vrr;
      Double I_ERI_Fx2z_S_D2y_S_M2_vrr = QCY*I_ERI_Fx2z_S_Py_S_M2_vrr+WQY*I_ERI_Fx2z_S_Py_S_M3_vrr+oned2e*I_ERI_Fx2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_D2y_S_M2_vrr = QCY*I_ERI_F3y_S_Py_S_M2_vrr+WQY*I_ERI_F3y_S_Py_S_M3_vrr+oned2e*I_ERI_F3y_S_S_S_M2_vrr-rhod2esq*I_ERI_F3y_S_S_S_M3_vrr+3*oned2k*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_F2yz_S_D2y_S_M2_vrr = QCY*I_ERI_F2yz_S_Py_S_M2_vrr+WQY*I_ERI_F2yz_S_Py_S_M3_vrr+oned2e*I_ERI_F2yz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_Py_S_M3_vrr;
      Double I_ERI_Fy2z_S_D2y_S_M2_vrr = QCY*I_ERI_Fy2z_S_Py_S_M2_vrr+WQY*I_ERI_Fy2z_S_Py_S_M3_vrr+oned2e*I_ERI_Fy2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M3_vrr+oned2k*I_ERI_D2z_S_Py_S_M3_vrr;
      Double I_ERI_F3z_S_D2y_S_M2_vrr = QCY*I_ERI_F3z_S_Py_S_M2_vrr+WQY*I_ERI_F3z_S_Py_S_M3_vrr+oned2e*I_ERI_F3z_S_S_S_M2_vrr-rhod2esq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_F3x_S_Dyz_S_M2_vrr = QCZ*I_ERI_F3x_S_Py_S_M2_vrr+WQZ*I_ERI_F3x_S_Py_S_M3_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M2_vrr = QCZ*I_ERI_F2xy_S_Py_S_M2_vrr+WQZ*I_ERI_F2xy_S_Py_S_M3_vrr;
      Double I_ERI_F2xz_S_Dyz_S_M2_vrr = QCZ*I_ERI_F2xz_S_Py_S_M2_vrr+WQZ*I_ERI_F2xz_S_Py_S_M3_vrr+oned2k*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_F3y_S_Dyz_S_M2_vrr = QCZ*I_ERI_F3y_S_Py_S_M2_vrr+WQZ*I_ERI_F3y_S_Py_S_M3_vrr;
      Double I_ERI_F2yz_S_Dyz_S_M2_vrr = QCZ*I_ERI_F2yz_S_Py_S_M2_vrr+WQZ*I_ERI_F2yz_S_Py_S_M3_vrr+oned2k*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_F3z_S_Dyz_S_M2_vrr = QCZ*I_ERI_F3z_S_Py_S_M2_vrr+WQZ*I_ERI_F3z_S_Py_S_M3_vrr+3*oned2k*I_ERI_D2z_S_Py_S_M3_vrr;
      Double I_ERI_F3x_S_D2z_S_M2_vrr = QCZ*I_ERI_F3x_S_Pz_S_M2_vrr+WQZ*I_ERI_F3x_S_Pz_S_M3_vrr+oned2e*I_ERI_F3x_S_S_S_M2_vrr-rhod2esq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_D2z_S_M2_vrr = QCZ*I_ERI_F2xy_S_Pz_S_M2_vrr+WQZ*I_ERI_F2xy_S_Pz_S_M3_vrr+oned2e*I_ERI_F2xy_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_D2z_S_M2_vrr = QCZ*I_ERI_F2xz_S_Pz_S_M2_vrr+WQZ*I_ERI_F2xz_S_Pz_S_M3_vrr+oned2e*I_ERI_F2xz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_S_S_M3_vrr+oned2k*I_ERI_D2x_S_Pz_S_M3_vrr;
      Double I_ERI_Fx2y_S_D2z_S_M2_vrr = QCZ*I_ERI_Fx2y_S_Pz_S_M2_vrr+WQZ*I_ERI_Fx2y_S_Pz_S_M3_vrr+oned2e*I_ERI_Fx2y_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Fxyz_S_D2z_S_M2_vrr = QCZ*I_ERI_Fxyz_S_Pz_S_M2_vrr+WQZ*I_ERI_Fxyz_S_Pz_S_M3_vrr+oned2e*I_ERI_Fxyz_S_S_S_M2_vrr-rhod2esq*I_ERI_Fxyz_S_S_S_M3_vrr+oned2k*I_ERI_Dxy_S_Pz_S_M3_vrr;
      Double I_ERI_Fx2z_S_D2z_S_M2_vrr = QCZ*I_ERI_Fx2z_S_Pz_S_M2_vrr+WQZ*I_ERI_Fx2z_S_Pz_S_M3_vrr+oned2e*I_ERI_Fx2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fx2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_Pz_S_M3_vrr;
      Double I_ERI_F3y_S_D2z_S_M2_vrr = QCZ*I_ERI_F3y_S_Pz_S_M2_vrr+WQZ*I_ERI_F3y_S_Pz_S_M3_vrr+oned2e*I_ERI_F3y_S_S_S_M2_vrr-rhod2esq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_D2z_S_M2_vrr = QCZ*I_ERI_F2yz_S_Pz_S_M2_vrr+WQZ*I_ERI_F2yz_S_Pz_S_M3_vrr+oned2e*I_ERI_F2yz_S_S_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_S_S_M3_vrr+oned2k*I_ERI_D2y_S_Pz_S_M3_vrr;
      Double I_ERI_Fy2z_S_D2z_S_M2_vrr = QCZ*I_ERI_Fy2z_S_Pz_S_M2_vrr+WQZ*I_ERI_Fy2z_S_Pz_S_M3_vrr+oned2e*I_ERI_Fy2z_S_S_S_M2_vrr-rhod2esq*I_ERI_Fy2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_Pz_S_M3_vrr;
      Double I_ERI_F3z_S_D2z_S_M2_vrr = QCZ*I_ERI_F3z_S_Pz_S_M2_vrr+WQZ*I_ERI_F3z_S_Pz_S_M3_vrr+oned2e*I_ERI_F3z_S_S_S_M2_vrr-rhod2esq*I_ERI_F3z_S_S_S_M3_vrr+3*oned2k*I_ERI_D2z_S_Pz_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_F_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 12 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M1_vrr = QCX*I_ERI_D2x_S_D2x_S_M1_vrr+WQX*I_ERI_D2x_S_D2x_S_M2_vrr+2*oned2e*I_ERI_D2x_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_D2x_S_Px_S_M2_vrr+2*oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxy_S_F3x_S_M1_vrr = QCX*I_ERI_Dxy_S_D2x_S_M1_vrr+WQX*I_ERI_Dxy_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Dxy_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Dxy_S_Px_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_Dxz_S_F3x_S_M1_vrr = QCX*I_ERI_Dxz_S_D2x_S_M1_vrr+WQX*I_ERI_Dxz_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Dxz_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Dxz_S_Px_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F3x_S_M1_vrr = QCX*I_ERI_D2y_S_D2x_S_M1_vrr+WQX*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2e*I_ERI_D2y_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Dyz_S_F3x_S_M1_vrr = QCX*I_ERI_Dyz_S_D2x_S_M1_vrr+WQX*I_ERI_Dyz_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Dyz_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Dyz_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_F3x_S_M1_vrr = QCX*I_ERI_D2z_S_D2x_S_M1_vrr+WQX*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2e*I_ERI_D2z_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_F2xy_S_M1_vrr = QCY*I_ERI_D2x_S_D2x_S_M1_vrr+WQY*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xy_S_M1_vrr = QCY*I_ERI_Dxy_S_D2x_S_M1_vrr+WQY*I_ERI_Dxy_S_D2x_S_M2_vrr+oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxz_S_F2xy_S_M1_vrr = QCY*I_ERI_Dxz_S_D2x_S_M1_vrr+WQY*I_ERI_Dxz_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F2xy_S_M1_vrr = QCY*I_ERI_D2y_S_D2x_S_M1_vrr+WQY*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_Dyz_S_F2xy_S_M1_vrr = QCY*I_ERI_Dyz_S_D2x_S_M1_vrr+WQY*I_ERI_Dyz_S_D2x_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xy_S_M1_vrr = QCY*I_ERI_D2z_S_D2x_S_M1_vrr+WQY*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_F2xz_S_M1_vrr = QCZ*I_ERI_D2x_S_D2x_S_M1_vrr+WQZ*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xz_S_M1_vrr = QCZ*I_ERI_Dxy_S_D2x_S_M1_vrr+WQZ*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_Dxz_S_F2xz_S_M1_vrr = QCZ*I_ERI_Dxz_S_D2x_S_M1_vrr+WQZ*I_ERI_Dxz_S_D2x_S_M2_vrr+oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F2xz_S_M1_vrr = QCZ*I_ERI_D2y_S_D2x_S_M1_vrr+WQZ*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Dyz_S_F2xz_S_M1_vrr = QCZ*I_ERI_Dyz_S_D2x_S_M1_vrr+WQZ*I_ERI_Dyz_S_D2x_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xz_S_M1_vrr = QCZ*I_ERI_D2z_S_D2x_S_M1_vrr+WQZ*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M1_vrr = QCX*I_ERI_D2x_S_D2y_S_M1_vrr+WQX*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_M1_vrr = QCX*I_ERI_Dxy_S_D2y_S_M1_vrr+WQX*I_ERI_Dxy_S_D2y_S_M2_vrr+oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_Dxz_S_Fx2y_S_M1_vrr = QCX*I_ERI_Dxz_S_D2y_S_M1_vrr+WQX*I_ERI_Dxz_S_D2y_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M1_vrr = QCX*I_ERI_D2y_S_D2y_S_M1_vrr+WQX*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Dyz_S_Fx2y_S_M1_vrr = QCX*I_ERI_Dyz_S_D2y_S_M1_vrr+WQX*I_ERI_Dyz_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M1_vrr = QCX*I_ERI_D2z_S_D2y_S_M1_vrr+WQX*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M1_vrr = QCX*I_ERI_D2x_S_D2z_S_M1_vrr+WQX*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_M1_vrr = QCX*I_ERI_Dxy_S_D2z_S_M1_vrr+WQX*I_ERI_Dxy_S_D2z_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_Dxz_S_Fx2z_S_M1_vrr = QCX*I_ERI_Dxz_S_D2z_S_M1_vrr+WQX*I_ERI_Dxz_S_D2z_S_M2_vrr+oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M1_vrr = QCX*I_ERI_D2y_S_D2z_S_M1_vrr+WQX*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_M1_vrr = QCX*I_ERI_Dyz_S_D2z_S_M1_vrr+WQX*I_ERI_Dyz_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M1_vrr = QCX*I_ERI_D2z_S_D2z_S_M1_vrr+WQX*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_D2x_S_F3y_S_M1_vrr = QCY*I_ERI_D2x_S_D2y_S_M1_vrr+WQY*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2e*I_ERI_D2x_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_F3y_S_M1_vrr = QCY*I_ERI_Dxy_S_D2y_S_M1_vrr+WQY*I_ERI_Dxy_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Dxy_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Dxy_S_Py_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_Dxz_S_F3y_S_M1_vrr = QCY*I_ERI_Dxz_S_D2y_S_M1_vrr+WQY*I_ERI_Dxz_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Dxz_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Dxz_S_Py_S_M2_vrr;
      Double I_ERI_D2y_S_F3y_S_M1_vrr = QCY*I_ERI_D2y_S_D2y_S_M1_vrr+WQY*I_ERI_D2y_S_D2y_S_M2_vrr+2*oned2e*I_ERI_D2y_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_D2y_S_Py_S_M2_vrr+2*oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_Dyz_S_F3y_S_M1_vrr = QCY*I_ERI_Dyz_S_D2y_S_M1_vrr+WQY*I_ERI_Dyz_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Dyz_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Dyz_S_Py_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F3y_S_M1_vrr = QCY*I_ERI_D2z_S_D2y_S_M1_vrr+WQY*I_ERI_D2z_S_D2y_S_M2_vrr+2*oned2e*I_ERI_D2z_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_F2yz_S_M1_vrr = QCZ*I_ERI_D2x_S_D2y_S_M1_vrr+WQZ*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_F2yz_S_M1_vrr = QCZ*I_ERI_Dxy_S_D2y_S_M1_vrr+WQZ*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_Dxz_S_F2yz_S_M1_vrr = QCZ*I_ERI_Dxz_S_D2y_S_M1_vrr+WQZ*I_ERI_Dxz_S_D2y_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_F2yz_S_M1_vrr = QCZ*I_ERI_D2y_S_D2y_S_M1_vrr+WQZ*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Dyz_S_F2yz_S_M1_vrr = QCZ*I_ERI_Dyz_S_D2y_S_M1_vrr+WQZ*I_ERI_Dyz_S_D2y_S_M2_vrr+oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F2yz_S_M1_vrr = QCZ*I_ERI_D2z_S_D2y_S_M1_vrr+WQZ*I_ERI_D2z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_F3z_S_M1_vrr = QCZ*I_ERI_D2x_S_D2z_S_M1_vrr+WQZ*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2e*I_ERI_D2x_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_F3z_S_M1_vrr = QCZ*I_ERI_Dxy_S_D2z_S_M1_vrr+WQZ*I_ERI_Dxy_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Dxy_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Dxy_S_Pz_S_M2_vrr;
      Double I_ERI_Dxz_S_F3z_S_M1_vrr = QCZ*I_ERI_Dxz_S_D2z_S_M1_vrr+WQZ*I_ERI_Dxz_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Dxz_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Dxz_S_Pz_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_F3z_S_M1_vrr = QCZ*I_ERI_D2y_S_D2z_S_M1_vrr+WQZ*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2e*I_ERI_D2y_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Dyz_S_F3z_S_M1_vrr = QCZ*I_ERI_Dyz_S_D2z_S_M1_vrr+WQZ*I_ERI_Dyz_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Dyz_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Dyz_S_Pz_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_F3z_S_M1_vrr = QCZ*I_ERI_D2z_S_D2z_S_M1_vrr+WQZ*I_ERI_D2z_S_D2z_S_M2_vrr+2*oned2e*I_ERI_D2z_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_D2z_S_Pz_S_M2_vrr+2*oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 8 integrals are omitted 
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
       * shell quartet name: SQ_ERI_F_S_F_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 20 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M1_vrr = QCX*I_ERI_F3x_S_D2x_S_M1_vrr+WQX*I_ERI_F3x_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F3x_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F3x_S_Px_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F3x_S_M1_vrr = QCX*I_ERI_F2xy_S_D2x_S_M1_vrr+WQX*I_ERI_F2xy_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F2xy_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F2xy_S_Px_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F3x_S_M1_vrr = QCX*I_ERI_F2xz_S_D2x_S_M1_vrr+WQX*I_ERI_F2xz_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F2xz_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F2xz_S_Px_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3x_S_M1_vrr = QCX*I_ERI_Fx2y_S_D2x_S_M1_vrr+WQX*I_ERI_Fx2y_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Fx2y_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Fx2y_S_Px_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3x_S_M1_vrr = QCX*I_ERI_Fxyz_S_D2x_S_M1_vrr+WQX*I_ERI_Fxyz_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Fxyz_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Fxyz_S_Px_S_M2_vrr+oned2k*I_ERI_Dyz_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3x_S_M1_vrr = QCX*I_ERI_Fx2z_S_D2x_S_M1_vrr+WQX*I_ERI_Fx2z_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Fx2z_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Fx2z_S_Px_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F3x_S_M1_vrr = QCX*I_ERI_F3y_S_D2x_S_M1_vrr+WQX*I_ERI_F3y_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F3y_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_F3x_S_M1_vrr = QCX*I_ERI_F2yz_S_D2x_S_M1_vrr+WQX*I_ERI_F2yz_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F2yz_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3x_S_M1_vrr = QCX*I_ERI_Fy2z_S_D2x_S_M1_vrr+WQX*I_ERI_Fy2z_S_D2x_S_M2_vrr+2*oned2e*I_ERI_Fy2z_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_F3x_S_M1_vrr = QCX*I_ERI_F3z_S_D2x_S_M1_vrr+WQX*I_ERI_F3z_S_D2x_S_M2_vrr+2*oned2e*I_ERI_F3z_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_F2xy_S_M1_vrr = QCY*I_ERI_F3x_S_D2x_S_M1_vrr+WQY*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M1_vrr = QCY*I_ERI_F2xy_S_D2x_S_M1_vrr+WQY*I_ERI_F2xy_S_D2x_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M1_vrr = QCY*I_ERI_F2xz_S_D2x_S_M1_vrr+WQY*I_ERI_F2xz_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_M1_vrr = QCY*I_ERI_Fx2y_S_D2x_S_M1_vrr+WQY*I_ERI_Fx2y_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2xy_S_M1_vrr = QCY*I_ERI_Fxyz_S_D2x_S_M1_vrr+WQY*I_ERI_Fxyz_S_D2x_S_M2_vrr+oned2k*I_ERI_Dxz_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_M1_vrr = QCY*I_ERI_Fx2z_S_D2x_S_M1_vrr+WQY*I_ERI_Fx2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F2xy_S_M1_vrr = QCY*I_ERI_F3y_S_D2x_S_M1_vrr+WQY*I_ERI_F3y_S_D2x_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M1_vrr = QCY*I_ERI_F2yz_S_D2x_S_M1_vrr+WQY*I_ERI_F2yz_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_D2x_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2xy_S_M1_vrr = QCY*I_ERI_Fy2z_S_D2x_S_M1_vrr+WQY*I_ERI_Fy2z_S_D2x_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_F2xy_S_M1_vrr = QCY*I_ERI_F3z_S_D2x_S_M1_vrr+WQY*I_ERI_F3z_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_F2xz_S_M1_vrr = QCZ*I_ERI_F3x_S_D2x_S_M1_vrr+WQZ*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M1_vrr = QCZ*I_ERI_F2xy_S_D2x_S_M1_vrr+WQZ*I_ERI_F2xy_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M1_vrr = QCZ*I_ERI_F2xz_S_D2x_S_M1_vrr+WQZ*I_ERI_F2xz_S_D2x_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_M1_vrr = QCZ*I_ERI_Fx2y_S_D2x_S_M1_vrr+WQZ*I_ERI_Fx2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2xz_S_M1_vrr = QCZ*I_ERI_Fxyz_S_D2x_S_M1_vrr+WQZ*I_ERI_Fxyz_S_D2x_S_M2_vrr+oned2k*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_M1_vrr = QCZ*I_ERI_Fx2z_S_D2x_S_M1_vrr+WQZ*I_ERI_Fx2z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F2xz_S_M1_vrr = QCZ*I_ERI_F3y_S_D2x_S_M1_vrr+WQZ*I_ERI_F3y_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M1_vrr = QCZ*I_ERI_F2yz_S_D2x_S_M1_vrr+WQZ*I_ERI_F2yz_S_D2x_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2xz_S_M1_vrr = QCZ*I_ERI_Fy2z_S_D2x_S_M1_vrr+WQZ*I_ERI_Fy2z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_F2xz_S_M1_vrr = QCZ*I_ERI_F3z_S_D2x_S_M1_vrr+WQZ*I_ERI_F3z_S_D2x_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M1_vrr = QCX*I_ERI_F3x_S_D2y_S_M1_vrr+WQX*I_ERI_F3x_S_D2y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M1_vrr = QCX*I_ERI_F2xy_S_D2y_S_M1_vrr+WQX*I_ERI_F2xy_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M1_vrr = QCX*I_ERI_F2xz_S_D2y_S_M1_vrr+WQX*I_ERI_F2xz_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_M1_vrr = QCX*I_ERI_Fx2y_S_D2y_S_M1_vrr+WQX*I_ERI_Fx2y_S_D2y_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fxyz_S_Fx2y_S_M1_vrr = QCX*I_ERI_Fxyz_S_D2y_S_M1_vrr+WQX*I_ERI_Fxyz_S_D2y_S_M2_vrr+oned2k*I_ERI_Dyz_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_M1_vrr = QCX*I_ERI_Fx2z_S_D2y_S_M1_vrr+WQX*I_ERI_Fx2z_S_D2y_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M1_vrr = QCX*I_ERI_F3y_S_D2y_S_M1_vrr+WQX*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M1_vrr = QCX*I_ERI_F2yz_S_D2y_S_M1_vrr+WQX*I_ERI_F2yz_S_D2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_Fx2y_S_M1_vrr = QCX*I_ERI_Fy2z_S_D2y_S_M1_vrr+WQX*I_ERI_Fy2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M1_vrr = QCX*I_ERI_F3z_S_D2y_S_M1_vrr+WQX*I_ERI_F3z_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M1_vrr = QCX*I_ERI_F3x_S_D2z_S_M1_vrr+WQX*I_ERI_F3x_S_D2z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M1_vrr = QCX*I_ERI_F2xy_S_D2z_S_M1_vrr+WQX*I_ERI_F2xy_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M1_vrr = QCX*I_ERI_F2xz_S_D2z_S_M1_vrr+WQX*I_ERI_F2xz_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_M1_vrr = QCX*I_ERI_Fx2y_S_D2z_S_M1_vrr+WQX*I_ERI_Fx2y_S_D2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fxyz_S_Fx2z_S_M1_vrr = QCX*I_ERI_Fxyz_S_D2z_S_M1_vrr+WQX*I_ERI_Fxyz_S_D2z_S_M2_vrr+oned2k*I_ERI_Dyz_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_M1_vrr = QCX*I_ERI_Fx2z_S_D2z_S_M1_vrr+WQX*I_ERI_Fx2z_S_D2z_S_M2_vrr+oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M1_vrr = QCX*I_ERI_F3y_S_D2z_S_M1_vrr+WQX*I_ERI_F3y_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M1_vrr = QCX*I_ERI_F2yz_S_D2z_S_M1_vrr+WQX*I_ERI_F2yz_S_D2z_S_M2_vrr;
      Double I_ERI_Fy2z_S_Fx2z_S_M1_vrr = QCX*I_ERI_Fy2z_S_D2z_S_M1_vrr+WQX*I_ERI_Fy2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M1_vrr = QCX*I_ERI_F3z_S_D2z_S_M1_vrr+WQX*I_ERI_F3z_S_D2z_S_M2_vrr;
      Double I_ERI_F3x_S_F3y_S_M1_vrr = QCY*I_ERI_F3x_S_D2y_S_M1_vrr+WQY*I_ERI_F3x_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F3x_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_F3y_S_M1_vrr = QCY*I_ERI_F2xy_S_D2y_S_M1_vrr+WQY*I_ERI_F2xy_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F2xy_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F2xy_S_Py_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_F3y_S_M1_vrr = QCY*I_ERI_F2xz_S_D2y_S_M1_vrr+WQY*I_ERI_F2xz_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F2xz_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F2xz_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3y_S_M1_vrr = QCY*I_ERI_Fx2y_S_D2y_S_M1_vrr+WQY*I_ERI_Fx2y_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Fx2y_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Fx2y_S_Py_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3y_S_M1_vrr = QCY*I_ERI_Fxyz_S_D2y_S_M1_vrr+WQY*I_ERI_Fxyz_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Fxyz_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Fxyz_S_Py_S_M2_vrr+oned2k*I_ERI_Dxz_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3y_S_M1_vrr = QCY*I_ERI_Fx2z_S_D2y_S_M1_vrr+WQY*I_ERI_Fx2z_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Fx2z_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Fx2z_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_F3y_S_M1_vrr = QCY*I_ERI_F3y_S_D2y_S_M1_vrr+WQY*I_ERI_F3y_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F3y_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F3y_S_Py_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_F3y_S_M1_vrr = QCY*I_ERI_F2yz_S_D2y_S_M1_vrr+WQY*I_ERI_F2yz_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F2yz_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F2yz_S_Py_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_D2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3y_S_M1_vrr = QCY*I_ERI_Fy2z_S_D2y_S_M1_vrr+WQY*I_ERI_Fy2z_S_D2y_S_M2_vrr+2*oned2e*I_ERI_Fy2z_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_Fy2z_S_Py_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_F3y_S_M1_vrr = QCY*I_ERI_F3z_S_D2y_S_M1_vrr+WQY*I_ERI_F3z_S_D2y_S_M2_vrr+2*oned2e*I_ERI_F3z_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_F3x_S_F2yz_S_M1_vrr = QCZ*I_ERI_F3x_S_D2y_S_M1_vrr+WQZ*I_ERI_F3x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M1_vrr = QCZ*I_ERI_F2xy_S_D2y_S_M1_vrr+WQZ*I_ERI_F2xy_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M1_vrr = QCZ*I_ERI_F2xz_S_D2y_S_M1_vrr+WQZ*I_ERI_F2xz_S_D2y_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_M1_vrr = QCZ*I_ERI_Fx2y_S_D2y_S_M1_vrr+WQZ*I_ERI_Fx2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2yz_S_M1_vrr = QCZ*I_ERI_Fxyz_S_D2y_S_M1_vrr+WQZ*I_ERI_Fxyz_S_D2y_S_M2_vrr+oned2k*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_M1_vrr = QCZ*I_ERI_Fx2z_S_D2y_S_M1_vrr+WQZ*I_ERI_Fx2z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_F2yz_S_M1_vrr = QCZ*I_ERI_F3y_S_D2y_S_M1_vrr+WQZ*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M1_vrr = QCZ*I_ERI_F2yz_S_D2y_S_M1_vrr+WQZ*I_ERI_F2yz_S_D2y_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2yz_S_M1_vrr = QCZ*I_ERI_Fy2z_S_D2y_S_M1_vrr+WQZ*I_ERI_Fy2z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_F2yz_S_M1_vrr = QCZ*I_ERI_F3z_S_D2y_S_M1_vrr+WQZ*I_ERI_F3z_S_D2y_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_F3z_S_M1_vrr = QCZ*I_ERI_F3x_S_D2z_S_M1_vrr+WQZ*I_ERI_F3x_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F3x_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xy_S_F3z_S_M1_vrr = QCZ*I_ERI_F2xy_S_D2z_S_M1_vrr+WQZ*I_ERI_F2xy_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F2xy_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F2xy_S_Pz_S_M2_vrr;
      Double I_ERI_F2xz_S_F3z_S_M1_vrr = QCZ*I_ERI_F2xz_S_D2z_S_M1_vrr+WQZ*I_ERI_F2xz_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F2xz_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F2xz_S_Pz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3z_S_M1_vrr = QCZ*I_ERI_Fx2y_S_D2z_S_M1_vrr+WQZ*I_ERI_Fx2y_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Fx2y_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Fx2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3z_S_M1_vrr = QCZ*I_ERI_Fxyz_S_D2z_S_M1_vrr+WQZ*I_ERI_Fxyz_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Fxyz_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Fxyz_S_Pz_S_M2_vrr+oned2k*I_ERI_Dxy_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3z_S_M1_vrr = QCZ*I_ERI_Fx2z_S_D2z_S_M1_vrr+WQZ*I_ERI_Fx2z_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Fx2z_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Fx2z_S_Pz_S_M2_vrr+2*oned2k*I_ERI_Dxz_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_F3z_S_M1_vrr = QCZ*I_ERI_F3y_S_D2z_S_M1_vrr+WQZ*I_ERI_F3y_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F3y_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_F2yz_S_F3z_S_M1_vrr = QCZ*I_ERI_F2yz_S_D2z_S_M1_vrr+WQZ*I_ERI_F2yz_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F2yz_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F2yz_S_Pz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3z_S_M1_vrr = QCZ*I_ERI_Fy2z_S_D2z_S_M1_vrr+WQZ*I_ERI_Fy2z_S_D2z_S_M2_vrr+2*oned2e*I_ERI_Fy2z_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_Fy2z_S_Pz_S_M2_vrr+2*oned2k*I_ERI_Dyz_S_D2z_S_M2_vrr;
      Double I_ERI_F3z_S_F3z_S_M1_vrr = QCZ*I_ERI_F3z_S_D2z_S_M1_vrr+WQZ*I_ERI_F3z_S_D2z_S_M2_vrr+2*oned2e*I_ERI_F3z_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_F3z_S_Pz_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_M1_vrr = PAX*I_ERI_F3x_S_D2x_S_M1_vrr+WPX*I_ERI_F3x_S_D2x_S_M2_vrr+3*oned2z*I_ERI_D2x_S_D2x_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_D2x_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_D2x_S_M1_vrr = PAY*I_ERI_F3x_S_D2x_S_M1_vrr+WPY*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_G3xz_S_D2x_S_M1_vrr = PAZ*I_ERI_F3x_S_D2x_S_M1_vrr+WPZ*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2x_S_M1_vrr = PAY*I_ERI_F2xy_S_D2x_S_M1_vrr+WPY*I_ERI_F2xy_S_D2x_S_M2_vrr+oned2z*I_ERI_D2x_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2x_S_M1_vrr = PAZ*I_ERI_F2xy_S_D2x_S_M1_vrr+WPZ*I_ERI_F2xy_S_D2x_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2x_S_M1_vrr = PAZ*I_ERI_F2xz_S_D2x_S_M1_vrr+WPZ*I_ERI_F2xz_S_D2x_S_M2_vrr+oned2z*I_ERI_D2x_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2x_S_M1_vrr = PAX*I_ERI_F3y_S_D2x_S_M1_vrr+WPX*I_ERI_F3y_S_D2x_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2x_S_M1_vrr = PAZ*I_ERI_Fx2y_S_D2x_S_M1_vrr+WPZ*I_ERI_Fx2y_S_D2x_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2x_S_M1_vrr = PAY*I_ERI_Fx2z_S_D2x_S_M1_vrr+WPY*I_ERI_Fx2z_S_D2x_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2x_S_M1_vrr = PAX*I_ERI_F3z_S_D2x_S_M1_vrr+WPX*I_ERI_F3z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_D2x_S_M1_vrr = PAY*I_ERI_F3y_S_D2x_S_M1_vrr+WPY*I_ERI_F3y_S_D2x_S_M2_vrr+3*oned2z*I_ERI_D2y_S_D2x_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_G3yz_S_D2x_S_M1_vrr = PAZ*I_ERI_F3y_S_D2x_S_M1_vrr+WPZ*I_ERI_F3y_S_D2x_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2x_S_M1_vrr = PAZ*I_ERI_F2yz_S_D2x_S_M1_vrr+WPZ*I_ERI_F2yz_S_D2x_S_M2_vrr+oned2z*I_ERI_D2y_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2x_S_M1_vrr = PAY*I_ERI_F3z_S_D2x_S_M1_vrr+WPY*I_ERI_F3z_S_D2x_S_M2_vrr;
      Double I_ERI_G4z_S_D2x_S_M1_vrr = PAZ*I_ERI_F3z_S_D2x_S_M1_vrr+WPZ*I_ERI_F3z_S_D2x_S_M2_vrr+3*oned2z*I_ERI_D2z_S_D2x_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_G4x_S_Dxy_S_M1_vrr = PAX*I_ERI_F3x_S_Dxy_S_M1_vrr+WPX*I_ERI_F3x_S_Dxy_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Dxy_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G3xy_S_Dxy_S_M1_vrr = PAY*I_ERI_F3x_S_Dxy_S_M1_vrr+WPY*I_ERI_F3x_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Dxy_S_M1_vrr = PAZ*I_ERI_F3x_S_Dxy_S_M1_vrr+WPZ*I_ERI_F3x_S_Dxy_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_M1_vrr = PAY*I_ERI_F2xy_S_Dxy_S_M1_vrr+WPY*I_ERI_F2xy_S_Dxy_S_M2_vrr+oned2z*I_ERI_D2x_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dxy_S_M2_vrr+oned2k*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_Dxy_S_M1_vrr = PAZ*I_ERI_F2xy_S_Dxy_S_M1_vrr+WPZ*I_ERI_F2xy_S_Dxy_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dxy_S_M1_vrr = PAZ*I_ERI_F2xz_S_Dxy_S_M1_vrr+WPZ*I_ERI_F2xz_S_Dxy_S_M2_vrr+oned2z*I_ERI_D2x_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_M1_vrr = PAX*I_ERI_F3y_S_Dxy_S_M1_vrr+WPX*I_ERI_F3y_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Dxy_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Dxy_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Dxy_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Dxy_S_M1_vrr = PAY*I_ERI_Fx2z_S_Dxy_S_M1_vrr+WPY*I_ERI_Fx2z_S_Dxy_S_M2_vrr+oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_M1_vrr = PAX*I_ERI_F3z_S_Dxy_S_M1_vrr+WPX*I_ERI_F3z_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4y_S_Dxy_S_M1_vrr = PAY*I_ERI_F3y_S_Dxy_S_M1_vrr+WPY*I_ERI_F3y_S_Dxy_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Dxy_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Dxy_S_M1_vrr = PAZ*I_ERI_F3y_S_Dxy_S_M1_vrr+WPZ*I_ERI_F3y_S_Dxy_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dxy_S_M1_vrr = PAZ*I_ERI_F2yz_S_Dxy_S_M1_vrr+WPZ*I_ERI_F2yz_S_Dxy_S_M2_vrr+oned2z*I_ERI_D2y_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_M1_vrr = PAY*I_ERI_F3z_S_Dxy_S_M1_vrr+WPY*I_ERI_F3z_S_Dxy_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Dxy_S_M1_vrr = PAZ*I_ERI_F3z_S_Dxy_S_M1_vrr+WPZ*I_ERI_F3z_S_Dxy_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Dxy_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_G4x_S_Dxz_S_M1_vrr = PAX*I_ERI_F3x_S_Dxz_S_M1_vrr+WPX*I_ERI_F3x_S_Dxz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Dxz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xy_S_Dxz_S_M1_vrr = PAY*I_ERI_F3x_S_Dxz_S_M1_vrr+WPY*I_ERI_F3x_S_Dxz_S_M2_vrr;
      Double I_ERI_G3xz_S_Dxz_S_M1_vrr = PAZ*I_ERI_F3x_S_Dxz_S_M1_vrr+WPZ*I_ERI_F3x_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dxz_S_M1_vrr = PAY*I_ERI_F2xy_S_Dxz_S_M1_vrr+WPY*I_ERI_F2xy_S_Dxz_S_M2_vrr+oned2z*I_ERI_D2x_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dxz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Dxz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Dxz_S_M2_vrr+oned2z*I_ERI_D2x_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dxz_S_M2_vrr+oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dxz_S_M1_vrr = PAX*I_ERI_F3y_S_Dxz_S_M1_vrr+WPX*I_ERI_F3y_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dxz_S_M1_vrr = PAX*I_ERI_F3z_S_Dxz_S_M1_vrr+WPX*I_ERI_F3z_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4y_S_Dxz_S_M1_vrr = PAY*I_ERI_F3y_S_Dxz_S_M1_vrr+WPY*I_ERI_F3y_S_Dxz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Dxz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_G3yz_S_Dxz_S_M1_vrr = PAZ*I_ERI_F3y_S_Dxz_S_M1_vrr+WPZ*I_ERI_F3y_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dxz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Dxz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Dxz_S_M2_vrr+oned2z*I_ERI_D2y_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Dxz_S_M2_vrr+oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dxz_S_M1_vrr = PAY*I_ERI_F3z_S_Dxz_S_M1_vrr+WPY*I_ERI_F3z_S_Dxz_S_M2_vrr;
      Double I_ERI_G4z_S_Dxz_S_M1_vrr = PAZ*I_ERI_F3z_S_Dxz_S_M1_vrr+WPZ*I_ERI_F3z_S_Dxz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Dxz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxz_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_D2y_S_M1_vrr = PAX*I_ERI_F3x_S_D2y_S_M1_vrr+WPX*I_ERI_F3x_S_D2y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_D2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_G3xy_S_D2y_S_M1_vrr = PAY*I_ERI_F3x_S_D2y_S_M1_vrr+WPY*I_ERI_F3x_S_D2y_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G3xz_S_D2y_S_M1_vrr = PAZ*I_ERI_F3x_S_D2y_S_M1_vrr+WPZ*I_ERI_F3x_S_D2y_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2y_S_M1_vrr = PAY*I_ERI_F2xy_S_D2y_S_M1_vrr+WPY*I_ERI_F2xy_S_D2y_S_M2_vrr+oned2z*I_ERI_D2x_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2y_S_M1_vrr = PAZ*I_ERI_F2xy_S_D2y_S_M1_vrr+WPZ*I_ERI_F2xy_S_D2y_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2y_S_M1_vrr = PAZ*I_ERI_F2xz_S_D2y_S_M1_vrr+WPZ*I_ERI_F2xz_S_D2y_S_M2_vrr+oned2z*I_ERI_D2x_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2y_S_M1_vrr = PAX*I_ERI_F3y_S_D2y_S_M1_vrr+WPX*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2y_S_M1_vrr = PAZ*I_ERI_Fx2y_S_D2y_S_M1_vrr+WPZ*I_ERI_Fx2y_S_D2y_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2y_S_M1_vrr = PAY*I_ERI_Fx2z_S_D2y_S_M1_vrr+WPY*I_ERI_Fx2z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Fx2z_S_Py_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2y_S_M1_vrr = PAX*I_ERI_F3z_S_D2y_S_M1_vrr+WPX*I_ERI_F3z_S_D2y_S_M2_vrr;
      Double I_ERI_G4y_S_D2y_S_M1_vrr = PAY*I_ERI_F3y_S_D2y_S_M1_vrr+WPY*I_ERI_F3y_S_D2y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_D2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_D2y_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G3yz_S_D2y_S_M1_vrr = PAZ*I_ERI_F3y_S_D2y_S_M1_vrr+WPZ*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2y_S_M1_vrr = PAZ*I_ERI_F2yz_S_D2y_S_M1_vrr+WPZ*I_ERI_F2yz_S_D2y_S_M2_vrr+oned2z*I_ERI_D2y_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2y_S_M1_vrr = PAY*I_ERI_F3z_S_D2y_S_M1_vrr+WPY*I_ERI_F3z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4z_S_D2y_S_M1_vrr = PAZ*I_ERI_F3z_S_D2y_S_M1_vrr+WPZ*I_ERI_F3z_S_D2y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_D2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_G4x_S_Dyz_S_M1_vrr = PAX*I_ERI_F3x_S_Dyz_S_M1_vrr+WPX*I_ERI_F3x_S_Dyz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Dyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_G3xy_S_Dyz_S_M1_vrr = PAY*I_ERI_F3x_S_Dyz_S_M1_vrr+WPY*I_ERI_F3x_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xz_S_Dyz_S_M1_vrr = PAZ*I_ERI_F3x_S_Dyz_S_M1_vrr+WPZ*I_ERI_F3x_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dyz_S_M1_vrr = PAY*I_ERI_F2xy_S_Dyz_S_M1_vrr+WPY*I_ERI_F2xy_S_Dyz_S_M2_vrr+oned2z*I_ERI_D2x_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dyz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Dyz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Dyz_S_M2_vrr+oned2z*I_ERI_D2x_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_F2xz_S_Py_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dyz_S_M1_vrr = PAX*I_ERI_F3y_S_Dyz_S_M1_vrr+WPX*I_ERI_F3y_S_Dyz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dyz_S_M1_vrr = PAX*I_ERI_F3z_S_Dyz_S_M1_vrr+WPX*I_ERI_F3z_S_Dyz_S_M2_vrr;
      Double I_ERI_G4y_S_Dyz_S_M1_vrr = PAY*I_ERI_F3y_S_Dyz_S_M1_vrr+WPY*I_ERI_F3y_S_Dyz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Dyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_G3yz_S_Dyz_S_M1_vrr = PAZ*I_ERI_F3y_S_Dyz_S_M1_vrr+WPZ*I_ERI_F3y_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dyz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Dyz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Dyz_S_M2_vrr+oned2z*I_ERI_D2y_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Dyz_S_M2_vrr+oned2k*I_ERI_F2yz_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dyz_S_M1_vrr = PAY*I_ERI_F3z_S_Dyz_S_M1_vrr+WPY*I_ERI_F3z_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4z_S_Dyz_S_M1_vrr = PAZ*I_ERI_F3z_S_Dyz_S_M1_vrr+WPZ*I_ERI_F3z_S_Dyz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Dyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Dyz_S_M2_vrr+oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4x_S_D2z_S_M1_vrr = PAX*I_ERI_F3x_S_D2z_S_M1_vrr+WPX*I_ERI_F3x_S_D2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_D2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_G3xy_S_D2z_S_M1_vrr = PAY*I_ERI_F3x_S_D2z_S_M1_vrr+WPY*I_ERI_F3x_S_D2z_S_M2_vrr;
      Double I_ERI_G3xz_S_D2z_S_M1_vrr = PAZ*I_ERI_F3x_S_D2z_S_M1_vrr+WPZ*I_ERI_F3x_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2z_S_M1_vrr = PAY*I_ERI_F2xy_S_D2z_S_M1_vrr+WPY*I_ERI_F2xy_S_D2z_S_M2_vrr+oned2z*I_ERI_D2x_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2z_S_M1_vrr = PAZ*I_ERI_F2xy_S_D2z_S_M1_vrr+WPZ*I_ERI_F2xy_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_D2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_D2z_S_M2_vrr+oned2z*I_ERI_D2x_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2z_S_M1_vrr = PAX*I_ERI_F3y_S_D2z_S_M1_vrr+WPX*I_ERI_F3y_S_D2z_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2z_S_M1_vrr = PAZ*I_ERI_Fx2y_S_D2z_S_M1_vrr+WPZ*I_ERI_Fx2y_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Fx2y_S_Pz_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2z_S_M1_vrr = PAY*I_ERI_Fx2z_S_D2z_S_M1_vrr+WPY*I_ERI_Fx2z_S_D2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2z_S_M1_vrr = PAX*I_ERI_F3z_S_D2z_S_M1_vrr+WPX*I_ERI_F3z_S_D2z_S_M2_vrr;
      Double I_ERI_G4y_S_D2z_S_M1_vrr = PAY*I_ERI_F3y_S_D2z_S_M1_vrr+WPY*I_ERI_F3y_S_D2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_D2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_G3yz_S_D2z_S_M1_vrr = PAZ*I_ERI_F3y_S_D2z_S_M1_vrr+WPZ*I_ERI_F3y_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_D2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_D2z_S_M2_vrr+oned2z*I_ERI_D2y_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Pz_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2z_S_M1_vrr = PAY*I_ERI_F3z_S_D2z_S_M1_vrr+WPY*I_ERI_F3z_S_D2z_S_M2_vrr;
      Double I_ERI_G4z_S_D2z_S_M1_vrr = PAZ*I_ERI_F3z_S_D2z_S_M1_vrr+WPZ*I_ERI_F3z_S_D2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_D2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_D2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Pz_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_F_S_G_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_vrr = QCX*I_ERI_F3x_S_F3x_S_vrr+WQX*I_ERI_F3x_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F3x_S_D2x_S_vrr-3*rhod2esq*I_ERI_F3x_S_D2x_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xy_S_G4x_S_vrr = QCX*I_ERI_F2xy_S_F3x_S_vrr+WQX*I_ERI_F2xy_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F2xy_S_D2x_S_vrr-3*rhod2esq*I_ERI_F2xy_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G4x_S_vrr = QCX*I_ERI_F2xz_S_F3x_S_vrr+WQX*I_ERI_F2xz_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F2xz_S_D2x_S_vrr-3*rhod2esq*I_ERI_F2xz_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4x_S_vrr = QCX*I_ERI_Fx2y_S_F3x_S_vrr+WQX*I_ERI_Fx2y_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Fx2y_S_D2x_S_vrr-3*rhod2esq*I_ERI_Fx2y_S_D2x_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fxyz_S_G4x_S_vrr = QCX*I_ERI_Fxyz_S_F3x_S_vrr+WQX*I_ERI_Fxyz_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Fxyz_S_D2x_S_vrr-3*rhod2esq*I_ERI_Fxyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Dyz_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4x_S_vrr = QCX*I_ERI_Fx2z_S_F3x_S_vrr+WQX*I_ERI_Fx2z_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Fx2z_S_D2x_S_vrr-3*rhod2esq*I_ERI_Fx2z_S_D2x_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G4x_S_vrr = QCX*I_ERI_F3y_S_F3x_S_vrr+WQX*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F3y_S_D2x_S_vrr-3*rhod2esq*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_G4x_S_vrr = QCX*I_ERI_F2yz_S_F3x_S_vrr+WQX*I_ERI_F2yz_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F2yz_S_D2x_S_vrr-3*rhod2esq*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4x_S_vrr = QCX*I_ERI_Fy2z_S_F3x_S_vrr+WQX*I_ERI_Fy2z_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Fy2z_S_D2x_S_vrr-3*rhod2esq*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_G4x_S_vrr = QCX*I_ERI_F3z_S_F3x_S_vrr+WQX*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2e*I_ERI_F3z_S_D2x_S_vrr-3*rhod2esq*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_G3xy_S_vrr = QCY*I_ERI_F3x_S_F3x_S_vrr+WQY*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xy_S_vrr = QCY*I_ERI_F2xy_S_F3x_S_vrr+WQY*I_ERI_F2xy_S_F3x_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xy_S_vrr = QCY*I_ERI_F2xz_S_F3x_S_vrr+WQY*I_ERI_F2xz_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_vrr = QCY*I_ERI_Fx2y_S_F3x_S_vrr+WQY*I_ERI_Fx2y_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3xy_S_vrr = QCY*I_ERI_Fxyz_S_F3x_S_vrr+WQY*I_ERI_Fxyz_S_F3x_S_M1_vrr+oned2k*I_ERI_Dxz_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_vrr = QCY*I_ERI_Fx2z_S_F3x_S_vrr+WQY*I_ERI_Fx2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G3xy_S_vrr = QCY*I_ERI_F3y_S_F3x_S_vrr+WQY*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xy_S_vrr = QCY*I_ERI_F2yz_S_F3x_S_vrr+WQY*I_ERI_F2yz_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3x_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3xy_S_vrr = QCY*I_ERI_Fy2z_S_F3x_S_vrr+WQY*I_ERI_Fy2z_S_F3x_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_G3xy_S_vrr = QCY*I_ERI_F3z_S_F3x_S_vrr+WQY*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_G3xz_S_vrr = QCZ*I_ERI_F3x_S_F3x_S_vrr+WQZ*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xz_S_vrr = QCZ*I_ERI_F2xy_S_F3x_S_vrr+WQZ*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xz_S_vrr = QCZ*I_ERI_F2xz_S_F3x_S_vrr+WQZ*I_ERI_F2xz_S_F3x_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_vrr = QCZ*I_ERI_Fx2y_S_F3x_S_vrr+WQZ*I_ERI_Fx2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3xz_S_vrr = QCZ*I_ERI_Fxyz_S_F3x_S_vrr+WQZ*I_ERI_Fxyz_S_F3x_S_M1_vrr+oned2k*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_vrr = QCZ*I_ERI_Fx2z_S_F3x_S_vrr+WQZ*I_ERI_Fx2z_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G3xz_S_vrr = QCZ*I_ERI_F3y_S_F3x_S_vrr+WQZ*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xz_S_vrr = QCZ*I_ERI_F2yz_S_F3x_S_vrr+WQZ*I_ERI_F2yz_S_F3x_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3xz_S_vrr = QCZ*I_ERI_Fy2z_S_F3x_S_vrr+WQZ*I_ERI_Fy2z_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_G3xz_S_vrr = QCZ*I_ERI_F3z_S_F3x_S_vrr+WQZ*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2y_S_vrr = QCY*I_ERI_F3x_S_F2xy_S_vrr+WQY*I_ERI_F3x_S_F2xy_S_M1_vrr+oned2e*I_ERI_F3x_S_D2x_S_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_vrr = QCY*I_ERI_F2xy_S_F2xy_S_vrr+WQY*I_ERI_F2xy_S_F2xy_S_M1_vrr+oned2e*I_ERI_F2xy_S_D2x_S_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_vrr = QCY*I_ERI_F2xz_S_F2xy_S_vrr+WQY*I_ERI_F2xz_S_F2xy_S_M1_vrr+oned2e*I_ERI_F2xz_S_D2x_S_vrr-rhod2esq*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_vrr = QCY*I_ERI_Fx2y_S_F2xy_S_vrr+WQY*I_ERI_Fx2y_S_F2xy_S_M1_vrr+oned2e*I_ERI_Fx2y_S_D2x_S_vrr-rhod2esq*I_ERI_Fx2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2x2y_S_vrr = QCY*I_ERI_Fxyz_S_F2xy_S_vrr+WQY*I_ERI_Fxyz_S_F2xy_S_M1_vrr+oned2e*I_ERI_Fxyz_S_D2x_S_vrr-rhod2esq*I_ERI_Fxyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Dxz_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_vrr = QCY*I_ERI_Fx2z_S_F2xy_S_vrr+WQY*I_ERI_Fx2z_S_F2xy_S_M1_vrr+oned2e*I_ERI_Fx2z_S_D2x_S_vrr-rhod2esq*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2y_S_vrr = QCY*I_ERI_F3y_S_F2xy_S_vrr+WQY*I_ERI_F3y_S_F2xy_S_M1_vrr+oned2e*I_ERI_F3y_S_D2x_S_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_vrr = QCY*I_ERI_F2yz_S_F2xy_S_vrr+WQY*I_ERI_F2yz_S_F2xy_S_M1_vrr+oned2e*I_ERI_F2yz_S_D2x_S_vrr-rhod2esq*I_ERI_F2yz_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F2xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2x2y_S_vrr = QCY*I_ERI_Fy2z_S_F2xy_S_vrr+WQY*I_ERI_Fy2z_S_F2xy_S_M1_vrr+oned2e*I_ERI_Fy2z_S_D2x_S_vrr-rhod2esq*I_ERI_Fy2z_S_D2x_S_M1_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2y_S_vrr = QCY*I_ERI_F3z_S_F2xy_S_vrr+WQY*I_ERI_F3z_S_F2xy_S_M1_vrr+oned2e*I_ERI_F3z_S_D2x_S_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_G2xyz_S_vrr = QCZ*I_ERI_F3x_S_F2xy_S_vrr+WQZ*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_vrr = QCZ*I_ERI_F2xy_S_F2xy_S_vrr+WQZ*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_vrr = QCZ*I_ERI_F2xz_S_F2xy_S_vrr+WQZ*I_ERI_F2xz_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_vrr = QCZ*I_ERI_Fx2y_S_F2xy_S_vrr+WQZ*I_ERI_Fx2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2xyz_S_vrr = QCZ*I_ERI_Fxyz_S_F2xy_S_vrr+WQZ*I_ERI_Fxyz_S_F2xy_S_M1_vrr+oned2k*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_vrr = QCZ*I_ERI_Fx2z_S_F2xy_S_vrr+WQZ*I_ERI_Fx2z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F2xy_S_M1_vrr;
      Double I_ERI_F3y_S_G2xyz_S_vrr = QCZ*I_ERI_F3y_S_F2xy_S_vrr+WQZ*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_vrr = QCZ*I_ERI_F2yz_S_F2xy_S_vrr+WQZ*I_ERI_F2yz_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2xyz_S_vrr = QCZ*I_ERI_Fy2z_S_F2xy_S_vrr+WQZ*I_ERI_Fy2z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F2xy_S_M1_vrr;
      Double I_ERI_F3z_S_G2xyz_S_vrr = QCZ*I_ERI_F3z_S_F2xy_S_vrr+WQZ*I_ERI_F3z_S_F2xy_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2z_S_vrr = QCZ*I_ERI_F3x_S_F2xz_S_vrr+WQZ*I_ERI_F3x_S_F2xz_S_M1_vrr+oned2e*I_ERI_F3x_S_D2x_S_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_vrr = QCZ*I_ERI_F2xy_S_F2xz_S_vrr+WQZ*I_ERI_F2xy_S_F2xz_S_M1_vrr+oned2e*I_ERI_F2xy_S_D2x_S_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_vrr = QCZ*I_ERI_F2xz_S_F2xz_S_vrr+WQZ*I_ERI_F2xz_S_F2xz_S_M1_vrr+oned2e*I_ERI_F2xz_S_D2x_S_vrr-rhod2esq*I_ERI_F2xz_S_D2x_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_vrr = QCZ*I_ERI_Fx2y_S_F2xz_S_vrr+WQZ*I_ERI_Fx2y_S_F2xz_S_M1_vrr+oned2e*I_ERI_Fx2y_S_D2x_S_vrr-rhod2esq*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2x2z_S_vrr = QCZ*I_ERI_Fxyz_S_F2xz_S_vrr+WQZ*I_ERI_Fxyz_S_F2xz_S_M1_vrr+oned2e*I_ERI_Fxyz_S_D2x_S_vrr-rhod2esq*I_ERI_Fxyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Dxy_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_vrr = QCZ*I_ERI_Fx2z_S_F2xz_S_vrr+WQZ*I_ERI_Fx2z_S_F2xz_S_M1_vrr+oned2e*I_ERI_Fx2z_S_D2x_S_vrr-rhod2esq*I_ERI_Fx2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F2xz_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2z_S_vrr = QCZ*I_ERI_F3y_S_F2xz_S_vrr+WQZ*I_ERI_F3y_S_F2xz_S_M1_vrr+oned2e*I_ERI_F3y_S_D2x_S_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_vrr = QCZ*I_ERI_F2yz_S_F2xz_S_vrr+WQZ*I_ERI_F2yz_S_F2xz_S_M1_vrr+oned2e*I_ERI_F2yz_S_D2x_S_vrr-rhod2esq*I_ERI_F2yz_S_D2x_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2x2z_S_vrr = QCZ*I_ERI_Fy2z_S_F2xz_S_vrr+WQZ*I_ERI_Fy2z_S_F2xz_S_M1_vrr+oned2e*I_ERI_Fy2z_S_D2x_S_vrr-rhod2esq*I_ERI_Fy2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2z_S_vrr = QCZ*I_ERI_F3z_S_F2xz_S_vrr+WQZ*I_ERI_F3z_S_F2xz_S_M1_vrr+oned2e*I_ERI_F3z_S_D2x_S_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3y_S_vrr = QCX*I_ERI_F3x_S_F3y_S_vrr+WQX*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_vrr = QCX*I_ERI_F2xy_S_F3y_S_vrr+WQX*I_ERI_F2xy_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_vrr = QCX*I_ERI_F2xz_S_F3y_S_vrr+WQX*I_ERI_F2xz_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_vrr = QCX*I_ERI_Fx2y_S_F3y_S_vrr+WQX*I_ERI_Fx2y_S_F3y_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx3y_S_vrr = QCX*I_ERI_Fxyz_S_F3y_S_vrr+WQX*I_ERI_Fxyz_S_F3y_S_M1_vrr+oned2k*I_ERI_Dyz_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_vrr = QCX*I_ERI_Fx2z_S_F3y_S_vrr+WQX*I_ERI_Fx2z_S_F3y_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3y_S_vrr = QCX*I_ERI_F3y_S_F3y_S_vrr+WQX*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_vrr = QCX*I_ERI_F2yz_S_F3y_S_vrr+WQX*I_ERI_F2yz_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx3y_S_vrr = QCX*I_ERI_Fy2z_S_F3y_S_vrr+WQX*I_ERI_Fy2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3y_S_vrr = QCX*I_ERI_F3z_S_F3y_S_vrr+WQX*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_vrr = QCZ*I_ERI_F3x_S_Fx2y_S_vrr+WQZ*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_vrr = QCZ*I_ERI_F2xy_S_Fx2y_S_vrr+WQZ*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_vrr = QCZ*I_ERI_F2xz_S_Fx2y_S_vrr+WQZ*I_ERI_F2xz_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_vrr = QCZ*I_ERI_Fx2y_S_Fx2y_S_vrr+WQZ*I_ERI_Fx2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx2yz_S_vrr = QCZ*I_ERI_Fxyz_S_Fx2y_S_vrr+WQZ*I_ERI_Fxyz_S_Fx2y_S_M1_vrr+oned2k*I_ERI_Dxy_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_vrr = QCZ*I_ERI_Fx2z_S_Fx2y_S_vrr+WQZ*I_ERI_Fx2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_vrr = QCZ*I_ERI_F3y_S_Fx2y_S_vrr+WQZ*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_vrr = QCZ*I_ERI_F2yz_S_Fx2y_S_vrr+WQZ*I_ERI_F2yz_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx2yz_S_vrr = QCZ*I_ERI_Fy2z_S_Fx2y_S_vrr+WQZ*I_ERI_Fy2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_vrr = QCZ*I_ERI_F3z_S_Fx2y_S_vrr+WQZ*I_ERI_F3z_S_Fx2y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_vrr = QCY*I_ERI_F3x_S_Fx2z_S_vrr+WQY*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_vrr = QCY*I_ERI_F2xy_S_Fx2z_S_vrr+WQY*I_ERI_F2xy_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_vrr = QCY*I_ERI_F2xz_S_Fx2z_S_vrr+WQY*I_ERI_F2xz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_vrr = QCY*I_ERI_Fx2y_S_Fx2z_S_vrr+WQY*I_ERI_Fx2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gxy2z_S_vrr = QCY*I_ERI_Fxyz_S_Fx2z_S_vrr+WQY*I_ERI_Fxyz_S_Fx2z_S_M1_vrr+oned2k*I_ERI_Dxz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_vrr = QCY*I_ERI_Fx2z_S_Fx2z_S_vrr+WQY*I_ERI_Fx2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_vrr = QCY*I_ERI_F3y_S_Fx2z_S_vrr+WQY*I_ERI_F3y_S_Fx2z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_vrr = QCY*I_ERI_F2yz_S_Fx2z_S_vrr+WQY*I_ERI_F2yz_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gxy2z_S_vrr = QCY*I_ERI_Fy2z_S_Fx2z_S_vrr+WQY*I_ERI_Fy2z_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_vrr = QCY*I_ERI_F3z_S_Fx2z_S_vrr+WQY*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3z_S_vrr = QCX*I_ERI_F3x_S_F3z_S_vrr+WQX*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_vrr = QCX*I_ERI_F2xy_S_F3z_S_vrr+WQX*I_ERI_F2xy_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_vrr = QCX*I_ERI_F2xz_S_F3z_S_vrr+WQX*I_ERI_F2xz_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_vrr = QCX*I_ERI_Fx2y_S_F3z_S_vrr+WQX*I_ERI_Fx2y_S_F3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx3z_S_vrr = QCX*I_ERI_Fxyz_S_F3z_S_vrr+WQX*I_ERI_Fxyz_S_F3z_S_M1_vrr+oned2k*I_ERI_Dyz_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_vrr = QCX*I_ERI_Fx2z_S_F3z_S_vrr+WQX*I_ERI_Fx2z_S_F3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3z_S_vrr = QCX*I_ERI_F3y_S_F3z_S_vrr+WQX*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_vrr = QCX*I_ERI_F2yz_S_F3z_S_vrr+WQX*I_ERI_F2yz_S_F3z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx3z_S_vrr = QCX*I_ERI_Fy2z_S_F3z_S_vrr+WQX*I_ERI_Fy2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3z_S_vrr = QCX*I_ERI_F3z_S_F3z_S_vrr+WQX*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_F3x_S_G4y_S_vrr = QCY*I_ERI_F3x_S_F3y_S_vrr+WQY*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F3x_S_D2y_S_vrr-3*rhod2esq*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_G4y_S_vrr = QCY*I_ERI_F2xy_S_F3y_S_vrr+WQY*I_ERI_F2xy_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F2xy_S_D2y_S_vrr-3*rhod2esq*I_ERI_F2xy_S_D2y_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_G4y_S_vrr = QCY*I_ERI_F2xz_S_F3y_S_vrr+WQY*I_ERI_F2xz_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F2xz_S_D2y_S_vrr-3*rhod2esq*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4y_S_vrr = QCY*I_ERI_Fx2y_S_F3y_S_vrr+WQY*I_ERI_Fx2y_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Fx2y_S_D2y_S_vrr-3*rhod2esq*I_ERI_Fx2y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_G4y_S_vrr = QCY*I_ERI_Fxyz_S_F3y_S_vrr+WQY*I_ERI_Fxyz_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Fxyz_S_D2y_S_vrr-3*rhod2esq*I_ERI_Fxyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Dxz_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4y_S_vrr = QCY*I_ERI_Fx2z_S_F3y_S_vrr+WQY*I_ERI_Fx2z_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Fx2z_S_D2y_S_vrr-3*rhod2esq*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_G4y_S_vrr = QCY*I_ERI_F3y_S_F3y_S_vrr+WQY*I_ERI_F3y_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F3y_S_D2y_S_vrr-3*rhod2esq*I_ERI_F3y_S_D2y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_G4y_S_vrr = QCY*I_ERI_F2yz_S_F3y_S_vrr+WQY*I_ERI_F2yz_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F2yz_S_D2y_S_vrr-3*rhod2esq*I_ERI_F2yz_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4y_S_vrr = QCY*I_ERI_Fy2z_S_F3y_S_vrr+WQY*I_ERI_Fy2z_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Fy2z_S_D2y_S_vrr-3*rhod2esq*I_ERI_Fy2z_S_D2y_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_G4y_S_vrr = QCY*I_ERI_F3z_S_F3y_S_vrr+WQY*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2e*I_ERI_F3z_S_D2y_S_vrr-3*rhod2esq*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_G3yz_S_vrr = QCZ*I_ERI_F3x_S_F3y_S_vrr+WQZ*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_G3yz_S_vrr = QCZ*I_ERI_F2xy_S_F3y_S_vrr+WQZ*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_G3yz_S_vrr = QCZ*I_ERI_F2xz_S_F3y_S_vrr+WQZ*I_ERI_F2xz_S_F3y_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_vrr = QCZ*I_ERI_Fx2y_S_F3y_S_vrr+WQZ*I_ERI_Fx2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3yz_S_vrr = QCZ*I_ERI_Fxyz_S_F3y_S_vrr+WQZ*I_ERI_Fxyz_S_F3y_S_M1_vrr+oned2k*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_vrr = QCZ*I_ERI_Fx2z_S_F3y_S_vrr+WQZ*I_ERI_Fx2z_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_G3yz_S_vrr = QCZ*I_ERI_F3y_S_F3y_S_vrr+WQZ*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_G3yz_S_vrr = QCZ*I_ERI_F2yz_S_F3y_S_vrr+WQZ*I_ERI_F2yz_S_F3y_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3yz_S_vrr = QCZ*I_ERI_Fy2z_S_F3y_S_vrr+WQZ*I_ERI_Fy2z_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_G3yz_S_vrr = QCZ*I_ERI_F3z_S_F3y_S_vrr+WQZ*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_G2y2z_S_vrr = QCZ*I_ERI_F3x_S_F2yz_S_vrr+WQZ*I_ERI_F3x_S_F2yz_S_M1_vrr+oned2e*I_ERI_F3x_S_D2y_S_vrr-rhod2esq*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_vrr = QCZ*I_ERI_F2xy_S_F2yz_S_vrr+WQZ*I_ERI_F2xy_S_F2yz_S_M1_vrr+oned2e*I_ERI_F2xy_S_D2y_S_vrr-rhod2esq*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_vrr = QCZ*I_ERI_F2xz_S_F2yz_S_vrr+WQZ*I_ERI_F2xz_S_F2yz_S_M1_vrr+oned2e*I_ERI_F2xz_S_D2y_S_vrr-rhod2esq*I_ERI_F2xz_S_D2y_S_M1_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_vrr = QCZ*I_ERI_Fx2y_S_F2yz_S_vrr+WQZ*I_ERI_Fx2y_S_F2yz_S_M1_vrr+oned2e*I_ERI_Fx2y_S_D2y_S_vrr-rhod2esq*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2y2z_S_vrr = QCZ*I_ERI_Fxyz_S_F2yz_S_vrr+WQZ*I_ERI_Fxyz_S_F2yz_S_M1_vrr+oned2e*I_ERI_Fxyz_S_D2y_S_vrr-rhod2esq*I_ERI_Fxyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Dxy_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_vrr = QCZ*I_ERI_Fx2z_S_F2yz_S_vrr+WQZ*I_ERI_Fx2z_S_F2yz_S_M1_vrr+oned2e*I_ERI_Fx2z_S_D2y_S_vrr-rhod2esq*I_ERI_Fx2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_G2y2z_S_vrr = QCZ*I_ERI_F3y_S_F2yz_S_vrr+WQZ*I_ERI_F3y_S_F2yz_S_M1_vrr+oned2e*I_ERI_F3y_S_D2y_S_vrr-rhod2esq*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_vrr = QCZ*I_ERI_F2yz_S_F2yz_S_vrr+WQZ*I_ERI_F2yz_S_F2yz_S_M1_vrr+oned2e*I_ERI_F2yz_S_D2y_S_vrr-rhod2esq*I_ERI_F2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2y2z_S_vrr = QCZ*I_ERI_Fy2z_S_F2yz_S_vrr+WQZ*I_ERI_Fy2z_S_F2yz_S_M1_vrr+oned2e*I_ERI_Fy2z_S_D2y_S_vrr-rhod2esq*I_ERI_Fy2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F2yz_S_M1_vrr;
      Double I_ERI_F3z_S_G2y2z_S_vrr = QCZ*I_ERI_F3z_S_F2yz_S_vrr+WQZ*I_ERI_F3z_S_F2yz_S_M1_vrr+oned2e*I_ERI_F3z_S_D2y_S_vrr-rhod2esq*I_ERI_F3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3x_S_Gy3z_S_vrr = QCY*I_ERI_F3x_S_F3z_S_vrr+WQY*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_vrr = QCY*I_ERI_F2xy_S_F3z_S_vrr+WQY*I_ERI_F2xy_S_F3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_vrr = QCY*I_ERI_F2xz_S_F3z_S_vrr+WQY*I_ERI_F2xz_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_vrr = QCY*I_ERI_Fx2y_S_F3z_S_vrr+WQY*I_ERI_Fx2y_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gy3z_S_vrr = QCY*I_ERI_Fxyz_S_F3z_S_vrr+WQY*I_ERI_Fxyz_S_F3z_S_M1_vrr+oned2k*I_ERI_Dxz_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_vrr = QCY*I_ERI_Fx2z_S_F3z_S_vrr+WQY*I_ERI_Fx2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gy3z_S_vrr = QCY*I_ERI_F3y_S_F3z_S_vrr+WQY*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_vrr = QCY*I_ERI_F2yz_S_F3z_S_vrr+WQY*I_ERI_F2yz_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gy3z_S_vrr = QCY*I_ERI_Fy2z_S_F3z_S_vrr+WQY*I_ERI_Fy2z_S_F3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_Gy3z_S_vrr = QCY*I_ERI_F3z_S_F3z_S_vrr+WQY*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_F3x_S_G4z_S_vrr = QCZ*I_ERI_F3x_S_F3z_S_vrr+WQZ*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F3x_S_D2z_S_vrr-3*rhod2esq*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G4z_S_vrr = QCZ*I_ERI_F2xy_S_F3z_S_vrr+WQZ*I_ERI_F2xy_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F2xy_S_D2z_S_vrr-3*rhod2esq*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G4z_S_vrr = QCZ*I_ERI_F2xz_S_F3z_S_vrr+WQZ*I_ERI_F2xz_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F2xz_S_D2z_S_vrr-3*rhod2esq*I_ERI_F2xz_S_D2z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4z_S_vrr = QCZ*I_ERI_Fx2y_S_F3z_S_vrr+WQZ*I_ERI_Fx2y_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Fx2y_S_D2z_S_vrr-3*rhod2esq*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_G4z_S_vrr = QCZ*I_ERI_Fxyz_S_F3z_S_vrr+WQZ*I_ERI_Fxyz_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Fxyz_S_D2z_S_vrr-3*rhod2esq*I_ERI_Fxyz_S_D2z_S_M1_vrr+oned2k*I_ERI_Dxy_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4z_S_vrr = QCZ*I_ERI_Fx2z_S_F3z_S_vrr+WQZ*I_ERI_Fx2z_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Fx2z_S_D2z_S_vrr-3*rhod2esq*I_ERI_Fx2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxz_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_G4z_S_vrr = QCZ*I_ERI_F3y_S_F3z_S_vrr+WQZ*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F3y_S_D2z_S_vrr-3*rhod2esq*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G4z_S_vrr = QCZ*I_ERI_F2yz_S_F3z_S_vrr+WQZ*I_ERI_F2yz_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F2yz_S_D2z_S_vrr-3*rhod2esq*I_ERI_F2yz_S_D2z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4z_S_vrr = QCZ*I_ERI_Fy2z_S_F3z_S_vrr+WQZ*I_ERI_Fy2z_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Fy2z_S_D2z_S_vrr-3*rhod2esq*I_ERI_Fy2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dyz_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_G4z_S_vrr = QCZ*I_ERI_F3z_S_F3z_S_vrr+WQZ*I_ERI_F3z_S_F3z_S_M1_vrr+3*oned2e*I_ERI_F3z_S_D2z_S_vrr-3*rhod2esq*I_ERI_F3z_S_D2z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_D_S
       * RHS shell quartet name: SQ_ERI_G_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_vrr = QCX*I_ERI_G4x_S_D2x_S_vrr+WQX*I_ERI_G4x_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Px_S_vrr-2*rhod2esq*I_ERI_G4x_S_Px_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F3x_S_vrr = QCX*I_ERI_G3xy_S_D2x_S_vrr+WQX*I_ERI_G3xy_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Px_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Px_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F3x_S_vrr = QCX*I_ERI_G3xz_S_D2x_S_vrr+WQX*I_ERI_G3xz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Px_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Px_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3x_S_vrr = QCX*I_ERI_G2x2y_S_D2x_S_vrr+WQX*I_ERI_G2x2y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Px_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3x_S_vrr = QCX*I_ERI_G2xyz_S_D2x_S_vrr+WQX*I_ERI_G2xyz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Px_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3x_S_vrr = QCX*I_ERI_G2x2z_S_D2x_S_vrr+WQX*I_ERI_G2x2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Px_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3x_S_vrr = QCX*I_ERI_Gx3y_S_D2x_S_vrr+WQX*I_ERI_Gx3y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3x_S_vrr = QCX*I_ERI_Gx2yz_S_D2x_S_vrr+WQX*I_ERI_Gx2yz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Px_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3x_S_vrr = QCX*I_ERI_Gxy2z_S_D2x_S_vrr+WQX*I_ERI_Gxy2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Px_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3x_S_vrr = QCX*I_ERI_Gx3z_S_D2x_S_vrr+WQX*I_ERI_Gx3z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F3x_S_vrr = QCX*I_ERI_G4y_S_D2x_S_vrr+WQX*I_ERI_G4y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Px_S_vrr-2*rhod2esq*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_F3x_S_vrr = QCX*I_ERI_G3yz_S_D2x_S_vrr+WQX*I_ERI_G3yz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Px_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3x_S_vrr = QCX*I_ERI_G2y2z_S_D2x_S_vrr+WQX*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Px_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3x_S_vrr = QCX*I_ERI_Gy3z_S_D2x_S_vrr+WQX*I_ERI_Gy3z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_F3x_S_vrr = QCX*I_ERI_G4z_S_D2x_S_vrr+WQX*I_ERI_G4z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Px_S_vrr-2*rhod2esq*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_F2xy_S_vrr = QCY*I_ERI_G4x_S_D2x_S_vrr+WQY*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xy_S_vrr = QCY*I_ERI_G3xy_S_D2x_S_vrr+WQY*I_ERI_G3xy_S_D2x_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xy_S_vrr = QCY*I_ERI_G3xz_S_D2x_S_vrr+WQY*I_ERI_G3xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_vrr = QCY*I_ERI_G2x2y_S_D2x_S_vrr+WQY*I_ERI_G2x2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xy_S_vrr = QCY*I_ERI_G2xyz_S_D2x_S_vrr+WQY*I_ERI_G2xyz_S_D2x_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xy_S_vrr = QCY*I_ERI_G2x2z_S_D2x_S_vrr+WQY*I_ERI_G2x2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_vrr = QCY*I_ERI_Gx3y_S_D2x_S_vrr+WQY*I_ERI_Gx3y_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xy_S_vrr = QCY*I_ERI_Gx2yz_S_D2x_S_vrr+WQY*I_ERI_Gx2yz_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xy_S_vrr = QCY*I_ERI_Gxy2z_S_D2x_S_vrr+WQY*I_ERI_Gxy2z_S_D2x_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_vrr = QCY*I_ERI_Gx3z_S_D2x_S_vrr+WQY*I_ERI_Gx3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F2xy_S_vrr = QCY*I_ERI_G4y_S_D2x_S_vrr+WQY*I_ERI_G4y_S_D2x_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xy_S_vrr = QCY*I_ERI_G3yz_S_D2x_S_vrr+WQY*I_ERI_G3yz_S_D2x_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xy_S_vrr = QCY*I_ERI_G2y2z_S_D2x_S_vrr+WQY*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_vrr = QCY*I_ERI_Gy3z_S_D2x_S_vrr+WQY*I_ERI_Gy3z_S_D2x_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xy_S_vrr = QCY*I_ERI_G4z_S_D2x_S_vrr+WQY*I_ERI_G4z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_F2xz_S_vrr = QCZ*I_ERI_G4x_S_D2x_S_vrr+WQZ*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xz_S_vrr = QCZ*I_ERI_G3xy_S_D2x_S_vrr+WQZ*I_ERI_G3xy_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xz_S_vrr = QCZ*I_ERI_G3xz_S_D2x_S_vrr+WQZ*I_ERI_G3xz_S_D2x_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_vrr = QCZ*I_ERI_G2x2y_S_D2x_S_vrr+WQZ*I_ERI_G2x2y_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xz_S_vrr = QCZ*I_ERI_G2xyz_S_D2x_S_vrr+WQZ*I_ERI_G2xyz_S_D2x_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xz_S_vrr = QCZ*I_ERI_G2x2z_S_D2x_S_vrr+WQZ*I_ERI_G2x2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_vrr = QCZ*I_ERI_Gx3y_S_D2x_S_vrr+WQZ*I_ERI_Gx3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xz_S_vrr = QCZ*I_ERI_Gx2yz_S_D2x_S_vrr+WQZ*I_ERI_Gx2yz_S_D2x_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xz_S_vrr = QCZ*I_ERI_Gxy2z_S_D2x_S_vrr+WQZ*I_ERI_Gxy2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_vrr = QCZ*I_ERI_Gx3z_S_D2x_S_vrr+WQZ*I_ERI_Gx3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F2xz_S_vrr = QCZ*I_ERI_G4y_S_D2x_S_vrr+WQZ*I_ERI_G4y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xz_S_vrr = QCZ*I_ERI_G3yz_S_D2x_S_vrr+WQZ*I_ERI_G3yz_S_D2x_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xz_S_vrr = QCZ*I_ERI_G2y2z_S_D2x_S_vrr+WQZ*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_vrr = QCZ*I_ERI_Gy3z_S_D2x_S_vrr+WQZ*I_ERI_Gy3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xz_S_vrr = QCZ*I_ERI_G4z_S_D2x_S_vrr+WQZ*I_ERI_G4z_S_D2x_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2y_S_vrr = QCX*I_ERI_G4x_S_D2y_S_vrr+WQX*I_ERI_G4x_S_D2y_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_vrr = QCX*I_ERI_G3xy_S_D2y_S_vrr+WQX*I_ERI_G3xy_S_D2y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_vrr = QCX*I_ERI_G3xz_S_D2y_S_vrr+WQX*I_ERI_G3xz_S_D2y_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_vrr = QCX*I_ERI_G2x2y_S_D2y_S_vrr+WQX*I_ERI_G2x2y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2y_S_vrr = QCX*I_ERI_G2xyz_S_D2y_S_vrr+WQX*I_ERI_G2xyz_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2y_S_vrr = QCX*I_ERI_G2x2z_S_D2y_S_vrr+WQX*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_vrr = QCX*I_ERI_Gx3y_S_D2y_S_vrr+WQX*I_ERI_Gx3y_S_D2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2y_S_vrr = QCX*I_ERI_Gx2yz_S_D2y_S_vrr+WQX*I_ERI_Gx2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2y_S_vrr = QCX*I_ERI_Gxy2z_S_D2y_S_vrr+WQX*I_ERI_Gxy2z_S_D2y_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_vrr = QCX*I_ERI_Gx3z_S_D2y_S_vrr+WQX*I_ERI_Gx3z_S_D2y_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2y_S_vrr = QCX*I_ERI_G4y_S_D2y_S_vrr+WQX*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_vrr = QCX*I_ERI_G3yz_S_D2y_S_vrr+WQX*I_ERI_G3yz_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2y_S_vrr = QCX*I_ERI_G2y2z_S_D2y_S_vrr+WQX*I_ERI_G2y2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_vrr = QCX*I_ERI_Gy3z_S_D2y_S_vrr+WQX*I_ERI_Gy3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2y_S_vrr = QCX*I_ERI_G4z_S_D2y_S_vrr+WQX*I_ERI_G4z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fxyz_S_vrr = QCZ*I_ERI_G4x_S_Dxy_S_vrr+WQZ*I_ERI_G4x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_vrr = QCZ*I_ERI_G3xy_S_Dxy_S_vrr+WQZ*I_ERI_G3xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_vrr = QCZ*I_ERI_G3xz_S_Dxy_S_vrr+WQZ*I_ERI_G3xz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_vrr = QCZ*I_ERI_G2x2y_S_Dxy_S_vrr+WQZ*I_ERI_G2x2y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fxyz_S_vrr = QCZ*I_ERI_G2xyz_S_Dxy_S_vrr+WQZ*I_ERI_G2xyz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fxyz_S_vrr = QCZ*I_ERI_G2x2z_S_Dxy_S_vrr+WQZ*I_ERI_G2x2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_vrr = QCZ*I_ERI_Gx3y_S_Dxy_S_vrr+WQZ*I_ERI_Gx3y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fxyz_S_vrr = QCZ*I_ERI_Gx2yz_S_Dxy_S_vrr+WQZ*I_ERI_Gx2yz_S_Dxy_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fxyz_S_vrr = QCZ*I_ERI_Gxy2z_S_Dxy_S_vrr+WQZ*I_ERI_Gxy2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_vrr = QCZ*I_ERI_Gx3z_S_Dxy_S_vrr+WQZ*I_ERI_Gx3z_S_Dxy_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4y_S_Fxyz_S_vrr = QCZ*I_ERI_G4y_S_Dxy_S_vrr+WQZ*I_ERI_G4y_S_Dxy_S_M1_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_vrr = QCZ*I_ERI_G3yz_S_Dxy_S_vrr+WQZ*I_ERI_G3yz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fxyz_S_vrr = QCZ*I_ERI_G2y2z_S_Dxy_S_vrr+WQZ*I_ERI_G2y2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_vrr = QCZ*I_ERI_Gy3z_S_Dxy_S_vrr+WQZ*I_ERI_Gy3z_S_Dxy_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4z_S_Fxyz_S_vrr = QCZ*I_ERI_G4z_S_Dxy_S_vrr+WQZ*I_ERI_G4z_S_Dxy_S_M1_vrr+4*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2z_S_vrr = QCX*I_ERI_G4x_S_D2z_S_vrr+WQX*I_ERI_G4x_S_D2z_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_vrr = QCX*I_ERI_G3xy_S_D2z_S_vrr+WQX*I_ERI_G3xy_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_vrr = QCX*I_ERI_G3xz_S_D2z_S_vrr+WQX*I_ERI_G3xz_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_vrr = QCX*I_ERI_G2x2y_S_D2z_S_vrr+WQX*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2z_S_vrr = QCX*I_ERI_G2xyz_S_D2z_S_vrr+WQX*I_ERI_G2xyz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2z_S_vrr = QCX*I_ERI_G2x2z_S_D2z_S_vrr+WQX*I_ERI_G2x2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_vrr = QCX*I_ERI_Gx3y_S_D2z_S_vrr+WQX*I_ERI_Gx3y_S_D2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2z_S_vrr = QCX*I_ERI_Gx2yz_S_D2z_S_vrr+WQX*I_ERI_Gx2yz_S_D2z_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2z_S_vrr = QCX*I_ERI_Gxy2z_S_D2z_S_vrr+WQX*I_ERI_Gxy2z_S_D2z_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_vrr = QCX*I_ERI_Gx3z_S_D2z_S_vrr+WQX*I_ERI_Gx3z_S_D2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2z_S_vrr = QCX*I_ERI_G4y_S_D2z_S_vrr+WQX*I_ERI_G4y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_vrr = QCX*I_ERI_G3yz_S_D2z_S_vrr+WQX*I_ERI_G3yz_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2z_S_vrr = QCX*I_ERI_G2y2z_S_D2z_S_vrr+WQX*I_ERI_G2y2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_vrr = QCX*I_ERI_Gy3z_S_D2z_S_vrr+WQX*I_ERI_Gy3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2z_S_vrr = QCX*I_ERI_G4z_S_D2z_S_vrr+WQX*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_G4x_S_F3y_S_vrr = QCY*I_ERI_G4x_S_D2y_S_vrr+WQY*I_ERI_G4x_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Py_S_vrr-2*rhod2esq*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_G3xy_S_F3y_S_vrr = QCY*I_ERI_G3xy_S_D2y_S_vrr+WQY*I_ERI_G3xy_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Py_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Py_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F3y_S_vrr = QCY*I_ERI_G3xz_S_D2y_S_vrr+WQY*I_ERI_G3xz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Py_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3y_S_vrr = QCY*I_ERI_G2x2y_S_D2y_S_vrr+WQY*I_ERI_G2x2y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Py_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Py_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3y_S_vrr = QCY*I_ERI_G2xyz_S_D2y_S_vrr+WQY*I_ERI_G2xyz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Py_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Py_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3y_S_vrr = QCY*I_ERI_G2x2z_S_D2y_S_vrr+WQY*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Py_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3y_S_vrr = QCY*I_ERI_Gx3y_S_D2y_S_vrr+WQY*I_ERI_Gx3y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Py_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3y_S_vrr = QCY*I_ERI_Gx2yz_S_D2y_S_vrr+WQY*I_ERI_Gx2yz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Py_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3y_S_vrr = QCY*I_ERI_Gxy2z_S_D2y_S_vrr+WQY*I_ERI_Gxy2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Py_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3y_S_vrr = QCY*I_ERI_Gx3z_S_D2y_S_vrr+WQY*I_ERI_Gx3z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_G4y_S_F3y_S_vrr = QCY*I_ERI_G4y_S_D2y_S_vrr+WQY*I_ERI_G4y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Py_S_vrr-2*rhod2esq*I_ERI_G4y_S_Py_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F3y_S_vrr = QCY*I_ERI_G3yz_S_D2y_S_vrr+WQY*I_ERI_G3yz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Py_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Py_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3y_S_vrr = QCY*I_ERI_G2y2z_S_D2y_S_vrr+WQY*I_ERI_G2y2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Py_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3y_S_vrr = QCY*I_ERI_Gy3z_S_D2y_S_vrr+WQY*I_ERI_Gy3z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F3y_S_vrr = QCY*I_ERI_G4z_S_D2y_S_vrr+WQY*I_ERI_G4z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Py_S_vrr-2*rhod2esq*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_G4x_S_F2yz_S_vrr = QCZ*I_ERI_G4x_S_D2y_S_vrr+WQZ*I_ERI_G4x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_F2yz_S_vrr = QCZ*I_ERI_G3xy_S_D2y_S_vrr+WQZ*I_ERI_G3xy_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F2yz_S_vrr = QCZ*I_ERI_G3xz_S_D2y_S_vrr+WQZ*I_ERI_G3xz_S_D2y_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_vrr = QCZ*I_ERI_G2x2y_S_D2y_S_vrr+WQZ*I_ERI_G2x2y_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2yz_S_vrr = QCZ*I_ERI_G2xyz_S_D2y_S_vrr+WQZ*I_ERI_G2xyz_S_D2y_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2yz_S_vrr = QCZ*I_ERI_G2x2z_S_D2y_S_vrr+WQZ*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_vrr = QCZ*I_ERI_Gx3y_S_D2y_S_vrr+WQZ*I_ERI_Gx3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2yz_S_vrr = QCZ*I_ERI_Gx2yz_S_D2y_S_vrr+WQZ*I_ERI_Gx2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2yz_S_vrr = QCZ*I_ERI_Gxy2z_S_D2y_S_vrr+WQZ*I_ERI_Gxy2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_vrr = QCZ*I_ERI_Gx3z_S_D2y_S_vrr+WQZ*I_ERI_Gx3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_F2yz_S_vrr = QCZ*I_ERI_G4y_S_D2y_S_vrr+WQZ*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F2yz_S_vrr = QCZ*I_ERI_G3yz_S_D2y_S_vrr+WQZ*I_ERI_G3yz_S_D2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2yz_S_vrr = QCZ*I_ERI_G2y2z_S_D2y_S_vrr+WQZ*I_ERI_G2y2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_vrr = QCZ*I_ERI_Gy3z_S_D2y_S_vrr+WQZ*I_ERI_Gy3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F2yz_S_vrr = QCZ*I_ERI_G4z_S_D2y_S_vrr+WQZ*I_ERI_G4z_S_D2y_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fy2z_S_vrr = QCY*I_ERI_G4x_S_D2z_S_vrr+WQY*I_ERI_G4x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_vrr = QCY*I_ERI_G3xy_S_D2z_S_vrr+WQY*I_ERI_G3xy_S_D2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_vrr = QCY*I_ERI_G3xz_S_D2z_S_vrr+WQY*I_ERI_G3xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_vrr = QCY*I_ERI_G2x2y_S_D2z_S_vrr+WQY*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fy2z_S_vrr = QCY*I_ERI_G2xyz_S_D2z_S_vrr+WQY*I_ERI_G2xyz_S_D2z_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fy2z_S_vrr = QCY*I_ERI_G2x2z_S_D2z_S_vrr+WQY*I_ERI_G2x2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_vrr = QCY*I_ERI_Gx3y_S_D2z_S_vrr+WQY*I_ERI_Gx3y_S_D2z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fy2z_S_vrr = QCY*I_ERI_Gx2yz_S_D2z_S_vrr+WQY*I_ERI_Gx2yz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fy2z_S_vrr = QCY*I_ERI_Gxy2z_S_D2z_S_vrr+WQY*I_ERI_Gxy2z_S_D2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_vrr = QCY*I_ERI_Gx3z_S_D2z_S_vrr+WQY*I_ERI_Gx3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fy2z_S_vrr = QCY*I_ERI_G4y_S_D2z_S_vrr+WQY*I_ERI_G4y_S_D2z_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_vrr = QCY*I_ERI_G3yz_S_D2z_S_vrr+WQY*I_ERI_G3yz_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fy2z_S_vrr = QCY*I_ERI_G2y2z_S_D2z_S_vrr+WQY*I_ERI_G2y2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_vrr = QCY*I_ERI_Gy3z_S_D2z_S_vrr+WQY*I_ERI_Gy3z_S_D2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fy2z_S_vrr = QCY*I_ERI_G4z_S_D2z_S_vrr+WQY*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_G4x_S_F3z_S_vrr = QCZ*I_ERI_G4x_S_D2z_S_vrr+WQZ*I_ERI_G4x_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xy_S_F3z_S_vrr = QCZ*I_ERI_G3xy_S_D2z_S_vrr+WQZ*I_ERI_G3xy_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_G3xz_S_F3z_S_vrr = QCZ*I_ERI_G3xz_S_D2z_S_vrr+WQZ*I_ERI_G3xz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3z_S_vrr = QCZ*I_ERI_G2x2y_S_D2z_S_vrr+WQZ*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Pz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3z_S_vrr = QCZ*I_ERI_G2xyz_S_D2z_S_vrr+WQZ*I_ERI_G2xyz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Pz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3z_S_vrr = QCZ*I_ERI_G2x2z_S_D2z_S_vrr+WQZ*I_ERI_G2x2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3z_S_vrr = QCZ*I_ERI_Gx3y_S_D2z_S_vrr+WQZ*I_ERI_Gx3y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3z_S_vrr = QCZ*I_ERI_Gx2yz_S_D2z_S_vrr+WQZ*I_ERI_Gx2yz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Pz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3z_S_vrr = QCZ*I_ERI_Gxy2z_S_D2z_S_vrr+WQZ*I_ERI_Gxy2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3z_S_vrr = QCZ*I_ERI_Gx3z_S_D2z_S_vrr+WQZ*I_ERI_Gx3z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Pz_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_F3z_S_vrr = QCZ*I_ERI_G4y_S_D2z_S_vrr+WQZ*I_ERI_G4y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_G3yz_S_F3z_S_vrr = QCZ*I_ERI_G3yz_S_D2z_S_vrr+WQZ*I_ERI_G3yz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3z_S_vrr = QCZ*I_ERI_G2y2z_S_D2z_S_vrr+WQZ*I_ERI_G2y2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3z_S_vrr = QCZ*I_ERI_Gy3z_S_D2z_S_vrr+WQZ*I_ERI_Gy3z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Pz_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_F3z_S_vrr = QCZ*I_ERI_G4z_S_D2z_S_vrr+WQZ*I_ERI_G4z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4z_S_Pz_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_D_S
       * RHS shell quartet name: SQ_ERI_G_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_D2x_S_vrr = PAX*I_ERI_G4x_S_D2x_S_vrr+WPX*I_ERI_G4x_S_D2x_S_M1_vrr+4*oned2z*I_ERI_F3x_S_D2x_S_vrr-4*rhod2zsq*I_ERI_F3x_S_D2x_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_D2x_S_vrr = PAY*I_ERI_G4x_S_D2x_S_vrr+WPY*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_H4xz_S_D2x_S_vrr = PAZ*I_ERI_G4x_S_D2x_S_vrr+WPZ*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2x_S_vrr = PAY*I_ERI_G3xy_S_D2x_S_vrr+WPY*I_ERI_G3xy_S_D2x_S_M1_vrr+oned2z*I_ERI_F3x_S_D2x_S_vrr-rhod2zsq*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2x_S_vrr = PAZ*I_ERI_G3xy_S_D2x_S_vrr+WPZ*I_ERI_G3xy_S_D2x_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2x_S_vrr = PAZ*I_ERI_G3xz_S_D2x_S_vrr+WPZ*I_ERI_G3xz_S_D2x_S_M1_vrr+oned2z*I_ERI_F3x_S_D2x_S_vrr-rhod2zsq*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2x_S_vrr = PAX*I_ERI_Gx3y_S_D2x_S_vrr+WPX*I_ERI_Gx3y_S_D2x_S_M1_vrr+oned2z*I_ERI_F3y_S_D2x_S_vrr-rhod2zsq*I_ERI_F3y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2x_S_vrr = PAZ*I_ERI_G2x2y_S_D2x_S_vrr+WPZ*I_ERI_G2x2y_S_D2x_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2x_S_vrr = PAY*I_ERI_G2x2z_S_D2x_S_vrr+WPY*I_ERI_G2x2z_S_D2x_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2x_S_vrr = PAX*I_ERI_Gx3z_S_D2x_S_vrr+WPX*I_ERI_Gx3z_S_D2x_S_M1_vrr+oned2z*I_ERI_F3z_S_D2x_S_vrr-rhod2zsq*I_ERI_F3z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2x_S_vrr = PAX*I_ERI_G4y_S_D2x_S_vrr+WPX*I_ERI_G4y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2x_S_vrr = PAZ*I_ERI_Gx3y_S_D2x_S_vrr+WPZ*I_ERI_Gx3y_S_D2x_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2x_S_vrr = PAX*I_ERI_G2y2z_S_D2x_S_vrr+WPX*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2x_S_vrr = PAY*I_ERI_Gx3z_S_D2x_S_vrr+WPY*I_ERI_Gx3z_S_D2x_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2x_S_vrr = PAX*I_ERI_G4z_S_D2x_S_vrr+WPX*I_ERI_G4z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_D2x_S_vrr = PAY*I_ERI_G4y_S_D2x_S_vrr+WPY*I_ERI_G4y_S_D2x_S_M1_vrr+4*oned2z*I_ERI_F3y_S_D2x_S_vrr-4*rhod2zsq*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_H4yz_S_D2x_S_vrr = PAZ*I_ERI_G4y_S_D2x_S_vrr+WPZ*I_ERI_G4y_S_D2x_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2x_S_vrr = PAZ*I_ERI_G3yz_S_D2x_S_vrr+WPZ*I_ERI_G3yz_S_D2x_S_M1_vrr+oned2z*I_ERI_F3y_S_D2x_S_vrr-rhod2zsq*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2x_S_vrr = PAY*I_ERI_Gy3z_S_D2x_S_vrr+WPY*I_ERI_Gy3z_S_D2x_S_M1_vrr+oned2z*I_ERI_F3z_S_D2x_S_vrr-rhod2zsq*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2x_S_vrr = PAY*I_ERI_G4z_S_D2x_S_vrr+WPY*I_ERI_G4z_S_D2x_S_M1_vrr;
      Double I_ERI_H5z_S_D2x_S_vrr = PAZ*I_ERI_G4z_S_D2x_S_vrr+WPZ*I_ERI_G4z_S_D2x_S_M1_vrr+4*oned2z*I_ERI_F3z_S_D2x_S_vrr-4*rhod2zsq*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_H5x_S_Dxy_S_vrr = PAX*I_ERI_G4x_S_Dxy_S_vrr+WPX*I_ERI_G4x_S_Dxy_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Dxy_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxy_S_vrr = PAY*I_ERI_G4x_S_Dxy_S_vrr+WPY*I_ERI_G4x_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxy_S_vrr = PAZ*I_ERI_G4x_S_Dxy_S_vrr+WPZ*I_ERI_G4x_S_Dxy_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxy_S_vrr = PAY*I_ERI_G3xy_S_Dxy_S_vrr+WPY*I_ERI_G3xy_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3x_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3x_S_Dxy_S_M1_vrr+oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxy_S_vrr = PAZ*I_ERI_G3xy_S_Dxy_S_vrr+WPZ*I_ERI_G3xy_S_Dxy_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxy_S_vrr = PAZ*I_ERI_G3xz_S_Dxy_S_vrr+WPZ*I_ERI_G3xz_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3x_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxy_S_vrr = PAX*I_ERI_Gx3y_S_Dxy_S_vrr+WPX*I_ERI_Gx3y_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3y_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3y_S_Dxy_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxy_S_vrr = PAZ*I_ERI_G2x2y_S_Dxy_S_vrr+WPZ*I_ERI_G2x2y_S_Dxy_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxy_S_vrr = PAY*I_ERI_G2x2z_S_Dxy_S_vrr+WPY*I_ERI_G2x2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxy_S_vrr = PAX*I_ERI_Gx3z_S_Dxy_S_vrr+WPX*I_ERI_Gx3z_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3z_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3z_S_Dxy_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxy_S_vrr = PAX*I_ERI_G4y_S_Dxy_S_vrr+WPX*I_ERI_G4y_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxy_S_vrr = PAZ*I_ERI_Gx3y_S_Dxy_S_vrr+WPZ*I_ERI_Gx3y_S_Dxy_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxy_S_vrr = PAX*I_ERI_G2y2z_S_Dxy_S_vrr+WPX*I_ERI_G2y2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxy_S_vrr = PAY*I_ERI_Gx3z_S_Dxy_S_vrr+WPY*I_ERI_Gx3z_S_Dxy_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxy_S_vrr = PAX*I_ERI_G4z_S_Dxy_S_vrr+WPX*I_ERI_G4z_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Dxy_S_vrr = PAY*I_ERI_G4y_S_Dxy_S_vrr+WPY*I_ERI_G4y_S_Dxy_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Dxy_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxy_S_vrr = PAZ*I_ERI_G4y_S_Dxy_S_vrr+WPZ*I_ERI_G4y_S_Dxy_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxy_S_vrr = PAZ*I_ERI_G3yz_S_Dxy_S_vrr+WPZ*I_ERI_G3yz_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3y_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxy_S_vrr = PAY*I_ERI_Gy3z_S_Dxy_S_vrr+WPY*I_ERI_Gy3z_S_Dxy_S_M1_vrr+oned2z*I_ERI_F3z_S_Dxy_S_vrr-rhod2zsq*I_ERI_F3z_S_Dxy_S_M1_vrr+oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxy_S_vrr = PAY*I_ERI_G4z_S_Dxy_S_vrr+WPY*I_ERI_G4z_S_Dxy_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Dxy_S_vrr = PAZ*I_ERI_G4z_S_Dxy_S_vrr+WPZ*I_ERI_G4z_S_Dxy_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Dxy_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_H5x_S_Dxz_S_vrr = PAX*I_ERI_G4x_S_Dxz_S_vrr+WPX*I_ERI_G4x_S_Dxz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Dxz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxz_S_vrr = PAY*I_ERI_G4x_S_Dxz_S_vrr+WPY*I_ERI_G4x_S_Dxz_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxz_S_vrr = PAZ*I_ERI_G4x_S_Dxz_S_vrr+WPZ*I_ERI_G4x_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxz_S_vrr = PAY*I_ERI_G3xy_S_Dxz_S_vrr+WPY*I_ERI_G3xy_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3x_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxz_S_vrr = PAZ*I_ERI_G3xy_S_Dxz_S_vrr+WPZ*I_ERI_G3xy_S_Dxz_S_M1_vrr+oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxz_S_vrr = PAZ*I_ERI_G3xz_S_Dxz_S_vrr+WPZ*I_ERI_G3xz_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3x_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3x_S_Dxz_S_M1_vrr+oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxz_S_vrr = PAX*I_ERI_Gx3y_S_Dxz_S_vrr+WPX*I_ERI_Gx3y_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3y_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3y_S_Dxz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxz_S_vrr = PAZ*I_ERI_G2x2y_S_Dxz_S_vrr+WPZ*I_ERI_G2x2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxz_S_vrr = PAY*I_ERI_G2x2z_S_Dxz_S_vrr+WPY*I_ERI_G2x2z_S_Dxz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxz_S_vrr = PAX*I_ERI_Gx3z_S_Dxz_S_vrr+WPX*I_ERI_Gx3z_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3z_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3z_S_Dxz_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxz_S_vrr = PAX*I_ERI_G4y_S_Dxz_S_vrr+WPX*I_ERI_G4y_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxz_S_vrr = PAZ*I_ERI_Gx3y_S_Dxz_S_vrr+WPZ*I_ERI_Gx3y_S_Dxz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxz_S_vrr = PAX*I_ERI_G2y2z_S_Dxz_S_vrr+WPX*I_ERI_G2y2z_S_Dxz_S_M1_vrr+oned2k*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxz_S_vrr = PAY*I_ERI_Gx3z_S_Dxz_S_vrr+WPY*I_ERI_Gx3z_S_Dxz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxz_S_vrr = PAX*I_ERI_G4z_S_Dxz_S_vrr+WPX*I_ERI_G4z_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_Dxz_S_vrr = PAY*I_ERI_G4y_S_Dxz_S_vrr+WPY*I_ERI_G4y_S_Dxz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Dxz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxz_S_vrr = PAZ*I_ERI_G4y_S_Dxz_S_vrr+WPZ*I_ERI_G4y_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxz_S_vrr = PAZ*I_ERI_G3yz_S_Dxz_S_vrr+WPZ*I_ERI_G3yz_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3y_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3y_S_Dxz_S_M1_vrr+oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxz_S_vrr = PAY*I_ERI_Gy3z_S_Dxz_S_vrr+WPY*I_ERI_Gy3z_S_Dxz_S_M1_vrr+oned2z*I_ERI_F3z_S_Dxz_S_vrr-rhod2zsq*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxz_S_vrr = PAY*I_ERI_G4z_S_Dxz_S_vrr+WPY*I_ERI_G4z_S_Dxz_S_M1_vrr;
      Double I_ERI_H5z_S_Dxz_S_vrr = PAZ*I_ERI_G4z_S_Dxz_S_vrr+WPZ*I_ERI_G4z_S_Dxz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Dxz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Dxz_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_D2y_S_vrr = PAX*I_ERI_G4x_S_D2y_S_vrr+WPX*I_ERI_G4x_S_D2y_S_M1_vrr+4*oned2z*I_ERI_F3x_S_D2y_S_vrr-4*rhod2zsq*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_H4xy_S_D2y_S_vrr = PAY*I_ERI_G4x_S_D2y_S_vrr+WPY*I_ERI_G4x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H4xz_S_D2y_S_vrr = PAZ*I_ERI_G4x_S_D2y_S_vrr+WPZ*I_ERI_G4x_S_D2y_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2y_S_vrr = PAY*I_ERI_G3xy_S_D2y_S_vrr+WPY*I_ERI_G3xy_S_D2y_S_M1_vrr+oned2z*I_ERI_F3x_S_D2y_S_vrr-rhod2zsq*I_ERI_F3x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2y_S_vrr = PAZ*I_ERI_G3xy_S_D2y_S_vrr+WPZ*I_ERI_G3xy_S_D2y_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2y_S_vrr = PAZ*I_ERI_G3xz_S_D2y_S_vrr+WPZ*I_ERI_G3xz_S_D2y_S_M1_vrr+oned2z*I_ERI_F3x_S_D2y_S_vrr-rhod2zsq*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2y_S_vrr = PAX*I_ERI_Gx3y_S_D2y_S_vrr+WPX*I_ERI_Gx3y_S_D2y_S_M1_vrr+oned2z*I_ERI_F3y_S_D2y_S_vrr-rhod2zsq*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2y_S_vrr = PAZ*I_ERI_G2x2y_S_D2y_S_vrr+WPZ*I_ERI_G2x2y_S_D2y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2y_S_vrr = PAY*I_ERI_G2x2z_S_D2y_S_vrr+WPY*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2y_S_vrr = PAX*I_ERI_Gx3z_S_D2y_S_vrr+WPX*I_ERI_Gx3z_S_D2y_S_M1_vrr+oned2z*I_ERI_F3z_S_D2y_S_vrr-rhod2zsq*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2y_S_vrr = PAX*I_ERI_G4y_S_D2y_S_vrr+WPX*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2y_S_vrr = PAZ*I_ERI_Gx3y_S_D2y_S_vrr+WPZ*I_ERI_Gx3y_S_D2y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2y_S_vrr = PAX*I_ERI_G2y2z_S_D2y_S_vrr+WPX*I_ERI_G2y2z_S_D2y_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2y_S_vrr = PAY*I_ERI_Gx3z_S_D2y_S_vrr+WPY*I_ERI_Gx3z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2y_S_vrr = PAX*I_ERI_G4z_S_D2y_S_vrr+WPX*I_ERI_G4z_S_D2y_S_M1_vrr;
      Double I_ERI_H5y_S_D2y_S_vrr = PAY*I_ERI_G4y_S_D2y_S_vrr+WPY*I_ERI_G4y_S_D2y_S_M1_vrr+4*oned2z*I_ERI_F3y_S_D2y_S_vrr-4*rhod2zsq*I_ERI_F3y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H4yz_S_D2y_S_vrr = PAZ*I_ERI_G4y_S_D2y_S_vrr+WPZ*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2y_S_vrr = PAZ*I_ERI_G3yz_S_D2y_S_vrr+WPZ*I_ERI_G3yz_S_D2y_S_M1_vrr+oned2z*I_ERI_F3y_S_D2y_S_vrr-rhod2zsq*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2y_S_vrr = PAY*I_ERI_Gy3z_S_D2y_S_vrr+WPY*I_ERI_Gy3z_S_D2y_S_M1_vrr+oned2z*I_ERI_F3z_S_D2y_S_vrr-rhod2zsq*I_ERI_F3z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Py_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2y_S_vrr = PAY*I_ERI_G4z_S_D2y_S_vrr+WPY*I_ERI_G4z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5z_S_D2y_S_vrr = PAZ*I_ERI_G4z_S_D2y_S_vrr+WPZ*I_ERI_G4z_S_D2y_S_M1_vrr+4*oned2z*I_ERI_F3z_S_D2y_S_vrr-4*rhod2zsq*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_H5x_S_Dyz_S_vrr = PAX*I_ERI_G4x_S_Dyz_S_vrr+WPX*I_ERI_G4x_S_Dyz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Dyz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_H4xy_S_Dyz_S_vrr = PAY*I_ERI_G4x_S_Dyz_S_vrr+WPY*I_ERI_G4x_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xz_S_Dyz_S_vrr = PAZ*I_ERI_G4x_S_Dyz_S_vrr+WPZ*I_ERI_G4x_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dyz_S_vrr = PAY*I_ERI_G3xy_S_Dyz_S_vrr+WPY*I_ERI_G3xy_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3x_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3x_S_Dyz_S_M1_vrr+oned2k*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dyz_S_vrr = PAZ*I_ERI_G3xy_S_Dyz_S_vrr+WPZ*I_ERI_G3xy_S_Dyz_S_M1_vrr+oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dyz_S_vrr = PAZ*I_ERI_G3xz_S_Dyz_S_vrr+WPZ*I_ERI_G3xz_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3x_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3x_S_Dyz_S_M1_vrr+oned2k*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dyz_S_vrr = PAX*I_ERI_Gx3y_S_Dyz_S_vrr+WPX*I_ERI_Gx3y_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3y_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dyz_S_vrr = PAZ*I_ERI_G2x2y_S_Dyz_S_vrr+WPZ*I_ERI_G2x2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dyz_S_vrr = PAY*I_ERI_G2x2z_S_Dyz_S_vrr+WPY*I_ERI_G2x2z_S_Dyz_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dyz_S_vrr = PAX*I_ERI_Gx3z_S_Dyz_S_vrr+WPX*I_ERI_Gx3z_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3z_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dyz_S_vrr = PAX*I_ERI_G4y_S_Dyz_S_vrr+WPX*I_ERI_G4y_S_Dyz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dyz_S_vrr = PAZ*I_ERI_Gx3y_S_Dyz_S_vrr+WPZ*I_ERI_Gx3y_S_Dyz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dyz_S_vrr = PAX*I_ERI_G2y2z_S_Dyz_S_vrr+WPX*I_ERI_G2y2z_S_Dyz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dyz_S_vrr = PAY*I_ERI_Gx3z_S_Dyz_S_vrr+WPY*I_ERI_Gx3z_S_Dyz_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dyz_S_vrr = PAX*I_ERI_G4z_S_Dyz_S_vrr+WPX*I_ERI_G4z_S_Dyz_S_M1_vrr;
      Double I_ERI_H5y_S_Dyz_S_vrr = PAY*I_ERI_G4y_S_Dyz_S_vrr+WPY*I_ERI_G4y_S_Dyz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Dyz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_H4yz_S_Dyz_S_vrr = PAZ*I_ERI_G4y_S_Dyz_S_vrr+WPZ*I_ERI_G4y_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dyz_S_vrr = PAZ*I_ERI_G3yz_S_Dyz_S_vrr+WPZ*I_ERI_G3yz_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3y_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3y_S_Dyz_S_M1_vrr+oned2k*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dyz_S_vrr = PAY*I_ERI_Gy3z_S_Dyz_S_vrr+WPY*I_ERI_Gy3z_S_Dyz_S_M1_vrr+oned2z*I_ERI_F3z_S_Dyz_S_vrr-rhod2zsq*I_ERI_F3z_S_Dyz_S_M1_vrr+oned2k*I_ERI_Gy3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dyz_S_vrr = PAY*I_ERI_G4z_S_Dyz_S_vrr+WPY*I_ERI_G4z_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_Dyz_S_vrr = PAZ*I_ERI_G4z_S_Dyz_S_vrr+WPZ*I_ERI_G4z_S_Dyz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Dyz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Dyz_S_M1_vrr+oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_D2z_S_vrr = PAX*I_ERI_G4x_S_D2z_S_vrr+WPX*I_ERI_G4x_S_D2z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_D2z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_H4xy_S_D2z_S_vrr = PAY*I_ERI_G4x_S_D2z_S_vrr+WPY*I_ERI_G4x_S_D2z_S_M1_vrr;
      Double I_ERI_H4xz_S_D2z_S_vrr = PAZ*I_ERI_G4x_S_D2z_S_vrr+WPZ*I_ERI_G4x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2z_S_vrr = PAY*I_ERI_G3xy_S_D2z_S_vrr+WPY*I_ERI_G3xy_S_D2z_S_M1_vrr+oned2z*I_ERI_F3x_S_D2z_S_vrr-rhod2zsq*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2z_S_vrr = PAZ*I_ERI_G3xy_S_D2z_S_vrr+WPZ*I_ERI_G3xy_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2z_S_vrr = PAZ*I_ERI_G3xz_S_D2z_S_vrr+WPZ*I_ERI_G3xz_S_D2z_S_M1_vrr+oned2z*I_ERI_F3x_S_D2z_S_vrr-rhod2zsq*I_ERI_F3x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2z_S_vrr = PAX*I_ERI_Gx3y_S_D2z_S_vrr+WPX*I_ERI_Gx3y_S_D2z_S_M1_vrr+oned2z*I_ERI_F3y_S_D2z_S_vrr-rhod2zsq*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2z_S_vrr = PAZ*I_ERI_G2x2y_S_D2z_S_vrr+WPZ*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G2x2y_S_Pz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2z_S_vrr = PAY*I_ERI_G2x2z_S_D2z_S_vrr+WPY*I_ERI_G2x2z_S_D2z_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2z_S_vrr = PAX*I_ERI_Gx3z_S_D2z_S_vrr+WPX*I_ERI_Gx3z_S_D2z_S_M1_vrr+oned2z*I_ERI_F3z_S_D2z_S_vrr-rhod2zsq*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2z_S_vrr = PAX*I_ERI_G4y_S_D2z_S_vrr+WPX*I_ERI_G4y_S_D2z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2z_S_vrr = PAZ*I_ERI_Gx3y_S_D2z_S_vrr+WPZ*I_ERI_Gx3y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2z_S_vrr = PAX*I_ERI_G2y2z_S_D2z_S_vrr+WPX*I_ERI_G2y2z_S_D2z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2z_S_vrr = PAY*I_ERI_Gx3z_S_D2z_S_vrr+WPY*I_ERI_Gx3z_S_D2z_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2z_S_vrr = PAX*I_ERI_G4z_S_D2z_S_vrr+WPX*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_H5y_S_D2z_S_vrr = PAY*I_ERI_G4y_S_D2z_S_vrr+WPY*I_ERI_G4y_S_D2z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_D2z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_H4yz_S_D2z_S_vrr = PAZ*I_ERI_G4y_S_D2z_S_vrr+WPZ*I_ERI_G4y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2z_S_vrr = PAZ*I_ERI_G3yz_S_D2z_S_vrr+WPZ*I_ERI_G3yz_S_D2z_S_M1_vrr+oned2z*I_ERI_F3y_S_D2z_S_vrr-rhod2zsq*I_ERI_F3y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Pz_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2z_S_vrr = PAY*I_ERI_Gy3z_S_D2z_S_vrr+WPY*I_ERI_Gy3z_S_D2z_S_M1_vrr+oned2z*I_ERI_F3z_S_D2z_S_vrr-rhod2zsq*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2z_S_vrr = PAY*I_ERI_G4z_S_D2z_S_vrr+WPY*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_H5z_S_D2z_S_vrr = PAZ*I_ERI_G4z_S_D2z_S_vrr+WPZ*I_ERI_G4z_S_D2z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_D2z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_G_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_G_S_cc_coefs = gamma*gamma;
      I_ERI_F3x_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G4x_S_vrr;
      I_ERI_F2xy_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G4x_S_vrr;
      I_ERI_F2xz_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G4x_S_vrr;
      I_ERI_Fx2y_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G4x_S_vrr;
      I_ERI_Fxyz_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G4x_S_vrr;
      I_ERI_Fx2z_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G4x_S_vrr;
      I_ERI_F3y_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G4x_S_vrr;
      I_ERI_F2yz_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G4x_S_vrr;
      I_ERI_Fy2z_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G4x_S_vrr;
      I_ERI_F3z_S_G4x_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G4x_S_vrr;
      I_ERI_F3x_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G3xy_S_vrr;
      I_ERI_F2xy_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G3xy_S_vrr;
      I_ERI_F2xz_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G3xy_S_vrr;
      I_ERI_Fx2y_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G3xy_S_vrr;
      I_ERI_Fxyz_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G3xy_S_vrr;
      I_ERI_Fx2z_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G3xy_S_vrr;
      I_ERI_F3y_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G3xy_S_vrr;
      I_ERI_F2yz_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G3xy_S_vrr;
      I_ERI_Fy2z_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G3xy_S_vrr;
      I_ERI_F3z_S_G3xy_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G3xy_S_vrr;
      I_ERI_F3x_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G3xz_S_vrr;
      I_ERI_F2xy_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G3xz_S_vrr;
      I_ERI_F2xz_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G3xz_S_vrr;
      I_ERI_Fx2y_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G3xz_S_vrr;
      I_ERI_Fxyz_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G3xz_S_vrr;
      I_ERI_Fx2z_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G3xz_S_vrr;
      I_ERI_F3y_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G3xz_S_vrr;
      I_ERI_F2yz_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G3xz_S_vrr;
      I_ERI_Fy2z_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G3xz_S_vrr;
      I_ERI_F3z_S_G3xz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G3xz_S_vrr;
      I_ERI_F3x_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G2x2y_S_vrr;
      I_ERI_F2xy_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G2x2y_S_vrr;
      I_ERI_F2xz_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G2x2y_S_vrr;
      I_ERI_Fx2y_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G2x2y_S_vrr;
      I_ERI_Fxyz_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G2x2y_S_vrr;
      I_ERI_Fx2z_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G2x2y_S_vrr;
      I_ERI_F3y_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G2x2y_S_vrr;
      I_ERI_F2yz_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G2x2y_S_vrr;
      I_ERI_Fy2z_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G2x2y_S_vrr;
      I_ERI_F3z_S_G2x2y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G2x2y_S_vrr;
      I_ERI_F3x_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G2xyz_S_vrr;
      I_ERI_F2xy_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G2xyz_S_vrr;
      I_ERI_F2xz_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G2xyz_S_vrr;
      I_ERI_Fx2y_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G2xyz_S_vrr;
      I_ERI_Fxyz_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G2xyz_S_vrr;
      I_ERI_Fx2z_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G2xyz_S_vrr;
      I_ERI_F3y_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G2xyz_S_vrr;
      I_ERI_F2yz_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G2xyz_S_vrr;
      I_ERI_Fy2z_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G2xyz_S_vrr;
      I_ERI_F3z_S_G2xyz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G2xyz_S_vrr;
      I_ERI_F3x_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G2x2z_S_vrr;
      I_ERI_F2xy_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G2x2z_S_vrr;
      I_ERI_F2xz_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G2x2z_S_vrr;
      I_ERI_Fx2y_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G2x2z_S_vrr;
      I_ERI_Fxyz_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G2x2z_S_vrr;
      I_ERI_Fx2z_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G2x2z_S_vrr;
      I_ERI_F3y_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G2x2z_S_vrr;
      I_ERI_F2yz_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G2x2z_S_vrr;
      I_ERI_Fy2z_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G2x2z_S_vrr;
      I_ERI_F3z_S_G2x2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G2x2z_S_vrr;
      I_ERI_F3x_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_Gx3y_S_vrr;
      I_ERI_F2xy_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_Gx3y_S_vrr;
      I_ERI_F2xz_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_Gx3y_S_vrr;
      I_ERI_Fx2y_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_Gx3y_S_vrr;
      I_ERI_Fxyz_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_Gx3y_S_vrr;
      I_ERI_Fx2z_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_Gx3y_S_vrr;
      I_ERI_F3y_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_Gx3y_S_vrr;
      I_ERI_F2yz_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_Gx3y_S_vrr;
      I_ERI_Fy2z_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_Gx3y_S_vrr;
      I_ERI_F3z_S_Gx3y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_Gx3y_S_vrr;
      I_ERI_F3x_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_Gx2yz_S_vrr;
      I_ERI_F2xy_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_Gx2yz_S_vrr;
      I_ERI_F2xz_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_Gx2yz_S_vrr;
      I_ERI_Fx2y_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_Gx2yz_S_vrr;
      I_ERI_Fxyz_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_Gx2yz_S_vrr;
      I_ERI_Fx2z_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_Gx2yz_S_vrr;
      I_ERI_F3y_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_Gx2yz_S_vrr;
      I_ERI_F2yz_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_Gx2yz_S_vrr;
      I_ERI_Fy2z_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_Gx2yz_S_vrr;
      I_ERI_F3z_S_Gx2yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_Gx2yz_S_vrr;
      I_ERI_F3x_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_Gxy2z_S_vrr;
      I_ERI_F2xy_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_Gxy2z_S_vrr;
      I_ERI_F2xz_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_Gxy2z_S_vrr;
      I_ERI_Fx2y_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_Gxy2z_S_vrr;
      I_ERI_Fxyz_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_Gxy2z_S_vrr;
      I_ERI_Fx2z_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_Gxy2z_S_vrr;
      I_ERI_F3y_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_Gxy2z_S_vrr;
      I_ERI_F2yz_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_Gxy2z_S_vrr;
      I_ERI_Fy2z_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_Gxy2z_S_vrr;
      I_ERI_F3z_S_Gxy2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_Gxy2z_S_vrr;
      I_ERI_F3x_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_Gx3z_S_vrr;
      I_ERI_F2xy_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_Gx3z_S_vrr;
      I_ERI_F2xz_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_Gx3z_S_vrr;
      I_ERI_Fx2y_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_Gx3z_S_vrr;
      I_ERI_Fxyz_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_Gx3z_S_vrr;
      I_ERI_Fx2z_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_Gx3z_S_vrr;
      I_ERI_F3y_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_Gx3z_S_vrr;
      I_ERI_F2yz_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_Gx3z_S_vrr;
      I_ERI_Fy2z_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_Gx3z_S_vrr;
      I_ERI_F3z_S_Gx3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_Gx3z_S_vrr;
      I_ERI_F3x_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G4y_S_vrr;
      I_ERI_F2xy_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G4y_S_vrr;
      I_ERI_F2xz_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G4y_S_vrr;
      I_ERI_Fx2y_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G4y_S_vrr;
      I_ERI_Fxyz_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G4y_S_vrr;
      I_ERI_Fx2z_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G4y_S_vrr;
      I_ERI_F3y_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G4y_S_vrr;
      I_ERI_F2yz_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G4y_S_vrr;
      I_ERI_Fy2z_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G4y_S_vrr;
      I_ERI_F3z_S_G4y_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G4y_S_vrr;
      I_ERI_F3x_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G3yz_S_vrr;
      I_ERI_F2xy_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G3yz_S_vrr;
      I_ERI_F2xz_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G3yz_S_vrr;
      I_ERI_Fx2y_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G3yz_S_vrr;
      I_ERI_Fxyz_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G3yz_S_vrr;
      I_ERI_Fx2z_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G3yz_S_vrr;
      I_ERI_F3y_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G3yz_S_vrr;
      I_ERI_F2yz_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G3yz_S_vrr;
      I_ERI_Fy2z_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G3yz_S_vrr;
      I_ERI_F3z_S_G3yz_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G3yz_S_vrr;
      I_ERI_F3x_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G2y2z_S_vrr;
      I_ERI_F2xy_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G2y2z_S_vrr;
      I_ERI_F2xz_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G2y2z_S_vrr;
      I_ERI_Fx2y_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G2y2z_S_vrr;
      I_ERI_Fxyz_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G2y2z_S_vrr;
      I_ERI_Fx2z_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G2y2z_S_vrr;
      I_ERI_F3y_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G2y2z_S_vrr;
      I_ERI_F2yz_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G2y2z_S_vrr;
      I_ERI_Fy2z_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G2y2z_S_vrr;
      I_ERI_F3z_S_G2y2z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G2y2z_S_vrr;
      I_ERI_F3x_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_Gy3z_S_vrr;
      I_ERI_F2xy_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_Gy3z_S_vrr;
      I_ERI_F2xz_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_Gy3z_S_vrr;
      I_ERI_Fx2y_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_Gy3z_S_vrr;
      I_ERI_Fxyz_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_Gy3z_S_vrr;
      I_ERI_Fx2z_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_Gy3z_S_vrr;
      I_ERI_F3y_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_Gy3z_S_vrr;
      I_ERI_F2yz_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_Gy3z_S_vrr;
      I_ERI_Fy2z_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_Gy3z_S_vrr;
      I_ERI_F3z_S_Gy3z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_Gy3z_S_vrr;
      I_ERI_F3x_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3x_S_G4z_S_vrr;
      I_ERI_F2xy_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xy_S_G4z_S_vrr;
      I_ERI_F2xz_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2xz_S_G4z_S_vrr;
      I_ERI_Fx2y_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2y_S_G4z_S_vrr;
      I_ERI_Fxyz_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fxyz_S_G4z_S_vrr;
      I_ERI_Fx2z_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fx2z_S_G4z_S_vrr;
      I_ERI_F3y_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3y_S_G4z_S_vrr;
      I_ERI_F2yz_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F2yz_S_G4z_S_vrr;
      I_ERI_Fy2z_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_Fy2z_S_G4z_S_vrr;
      I_ERI_F3z_S_G4z_S_cc += SQ_ERI_F_S_G_S_cc_coefs*I_ERI_F3z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_c_coefs = gamma;
      I_ERI_F3x_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_c += SQ_ERI_F_S_D_S_c_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_S_S += I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S += I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S += I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S += I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S += I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S += I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S += I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S += I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S += I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S += I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_d_coefs = delta;
      I_ERI_F3x_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_d += SQ_ERI_F_S_D_S_d_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_bc_coefs = beta*gamma;
      I_ERI_G4x_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_bc += SQ_ERI_G_S_F_S_bc_coefs*I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_bc_coefs = beta*gamma;
      I_ERI_F3x_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_bc += SQ_ERI_F_S_F_S_bc_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_b_coefs = beta;
      I_ERI_G4x_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_b += SQ_ERI_G_S_P_S_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_H_S_D_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_D_S_bb_coefs = beta*beta;
      I_ERI_H5x_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_D2x_S_vrr;
      I_ERI_H4xy_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_D2x_S_vrr;
      I_ERI_H4xz_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_D2x_S_vrr;
      I_ERI_H3x2y_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_D2x_S_vrr;
      I_ERI_H3xyz_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_D2x_S_vrr;
      I_ERI_H3x2z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_D2x_S_vrr;
      I_ERI_H2x3y_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_D2x_S_vrr;
      I_ERI_H2x2yz_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_D2x_S_vrr;
      I_ERI_H2xy2z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_D2x_S_vrr;
      I_ERI_H2x3z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_D2x_S_vrr;
      I_ERI_Hx4y_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_D2x_S_vrr;
      I_ERI_Hx3yz_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_D2x_S_vrr;
      I_ERI_Hx2y2z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_D2x_S_vrr;
      I_ERI_Hxy3z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_D2x_S_vrr;
      I_ERI_Hx4z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_D2x_S_vrr;
      I_ERI_H5y_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_D2x_S_vrr;
      I_ERI_H4yz_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_D2x_S_vrr;
      I_ERI_H3y2z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_D2x_S_vrr;
      I_ERI_H2y3z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_D2x_S_vrr;
      I_ERI_Hy4z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_D2x_S_vrr;
      I_ERI_H5z_S_D2x_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_D2x_S_vrr;
      I_ERI_H5x_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_Dxy_S_vrr;
      I_ERI_H4xy_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_Dxy_S_vrr;
      I_ERI_H4xz_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_Dxy_S_vrr;
      I_ERI_H3x2y_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_Dxy_S_vrr;
      I_ERI_H3xyz_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_Dxy_S_vrr;
      I_ERI_H3x2z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_Dxy_S_vrr;
      I_ERI_H2x3y_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_Dxy_S_vrr;
      I_ERI_H2x2yz_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_Dxy_S_vrr;
      I_ERI_H2xy2z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_Dxy_S_vrr;
      I_ERI_H2x3z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_Dxy_S_vrr;
      I_ERI_Hx4y_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_Dxy_S_vrr;
      I_ERI_Hx3yz_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_Dxy_S_vrr;
      I_ERI_Hx2y2z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_Dxy_S_vrr;
      I_ERI_Hxy3z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_Dxy_S_vrr;
      I_ERI_Hx4z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_Dxy_S_vrr;
      I_ERI_H5y_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_Dxy_S_vrr;
      I_ERI_H4yz_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_Dxy_S_vrr;
      I_ERI_H3y2z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_Dxy_S_vrr;
      I_ERI_H2y3z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_Dxy_S_vrr;
      I_ERI_Hy4z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_Dxy_S_vrr;
      I_ERI_H5z_S_Dxy_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_Dxy_S_vrr;
      I_ERI_H5x_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_Dxz_S_vrr;
      I_ERI_H4xy_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_Dxz_S_vrr;
      I_ERI_H4xz_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_Dxz_S_vrr;
      I_ERI_H3x2y_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_Dxz_S_vrr;
      I_ERI_H3xyz_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_Dxz_S_vrr;
      I_ERI_H3x2z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_Dxz_S_vrr;
      I_ERI_H2x3y_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_Dxz_S_vrr;
      I_ERI_H2x2yz_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_Dxz_S_vrr;
      I_ERI_H2xy2z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_Dxz_S_vrr;
      I_ERI_H2x3z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_Dxz_S_vrr;
      I_ERI_Hx4y_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_Dxz_S_vrr;
      I_ERI_Hx3yz_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_Dxz_S_vrr;
      I_ERI_Hx2y2z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_Dxz_S_vrr;
      I_ERI_Hxy3z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_Dxz_S_vrr;
      I_ERI_Hx4z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_Dxz_S_vrr;
      I_ERI_H5y_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_Dxz_S_vrr;
      I_ERI_H4yz_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_Dxz_S_vrr;
      I_ERI_H3y2z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_Dxz_S_vrr;
      I_ERI_H2y3z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_Dxz_S_vrr;
      I_ERI_Hy4z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_Dxz_S_vrr;
      I_ERI_H5z_S_Dxz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_Dxz_S_vrr;
      I_ERI_H5x_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_D2y_S_vrr;
      I_ERI_H4xy_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_D2y_S_vrr;
      I_ERI_H4xz_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_D2y_S_vrr;
      I_ERI_H3x2y_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_D2y_S_vrr;
      I_ERI_H3xyz_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_D2y_S_vrr;
      I_ERI_H3x2z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_D2y_S_vrr;
      I_ERI_H2x3y_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_D2y_S_vrr;
      I_ERI_H2x2yz_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_D2y_S_vrr;
      I_ERI_H2xy2z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_D2y_S_vrr;
      I_ERI_H2x3z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_D2y_S_vrr;
      I_ERI_Hx4y_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_D2y_S_vrr;
      I_ERI_Hx3yz_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_D2y_S_vrr;
      I_ERI_Hx2y2z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_D2y_S_vrr;
      I_ERI_Hxy3z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_D2y_S_vrr;
      I_ERI_Hx4z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_D2y_S_vrr;
      I_ERI_H5y_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_D2y_S_vrr;
      I_ERI_H4yz_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_D2y_S_vrr;
      I_ERI_H3y2z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_D2y_S_vrr;
      I_ERI_H2y3z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_D2y_S_vrr;
      I_ERI_Hy4z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_D2y_S_vrr;
      I_ERI_H5z_S_D2y_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_D2y_S_vrr;
      I_ERI_H5x_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_Dyz_S_vrr;
      I_ERI_H4xy_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_Dyz_S_vrr;
      I_ERI_H4xz_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_Dyz_S_vrr;
      I_ERI_H3x2y_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_Dyz_S_vrr;
      I_ERI_H3xyz_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_Dyz_S_vrr;
      I_ERI_H3x2z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_Dyz_S_vrr;
      I_ERI_H2x3y_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_Dyz_S_vrr;
      I_ERI_H2x2yz_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_Dyz_S_vrr;
      I_ERI_H2xy2z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_Dyz_S_vrr;
      I_ERI_H2x3z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_Dyz_S_vrr;
      I_ERI_Hx4y_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_Dyz_S_vrr;
      I_ERI_Hx3yz_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_Dyz_S_vrr;
      I_ERI_Hx2y2z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_Dyz_S_vrr;
      I_ERI_Hxy3z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_Dyz_S_vrr;
      I_ERI_Hx4z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_Dyz_S_vrr;
      I_ERI_H5y_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_Dyz_S_vrr;
      I_ERI_H4yz_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_Dyz_S_vrr;
      I_ERI_H3y2z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_Dyz_S_vrr;
      I_ERI_H2y3z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_Dyz_S_vrr;
      I_ERI_Hy4z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_Dyz_S_vrr;
      I_ERI_H5z_S_Dyz_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_Dyz_S_vrr;
      I_ERI_H5x_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5x_S_D2z_S_vrr;
      I_ERI_H4xy_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xy_S_D2z_S_vrr;
      I_ERI_H4xz_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4xz_S_D2z_S_vrr;
      I_ERI_H3x2y_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2y_S_D2z_S_vrr;
      I_ERI_H3xyz_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3xyz_S_D2z_S_vrr;
      I_ERI_H3x2z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3x2z_S_D2z_S_vrr;
      I_ERI_H2x3y_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3y_S_D2z_S_vrr;
      I_ERI_H2x2yz_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x2yz_S_D2z_S_vrr;
      I_ERI_H2xy2z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2xy2z_S_D2z_S_vrr;
      I_ERI_H2x3z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2x3z_S_D2z_S_vrr;
      I_ERI_Hx4y_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4y_S_D2z_S_vrr;
      I_ERI_Hx3yz_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx3yz_S_D2z_S_vrr;
      I_ERI_Hx2y2z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx2y2z_S_D2z_S_vrr;
      I_ERI_Hxy3z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hxy3z_S_D2z_S_vrr;
      I_ERI_Hx4z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hx4z_S_D2z_S_vrr;
      I_ERI_H5y_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5y_S_D2z_S_vrr;
      I_ERI_H4yz_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H4yz_S_D2z_S_vrr;
      I_ERI_H3y2z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H3y2z_S_D2z_S_vrr;
      I_ERI_H2y3z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H2y3z_S_D2z_S_vrr;
      I_ERI_Hy4z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_Hy4z_S_D2z_S_vrr;
      I_ERI_H5z_S_D2z_S_bb += SQ_ERI_H_S_D_S_bb_coefs*I_ERI_H5z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_bb_coefs = beta*beta;
      I_ERI_G4x_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_bb += SQ_ERI_G_S_D_S_bb_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_bb_coefs = beta*beta;
      I_ERI_F3x_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_bb += SQ_ERI_F_S_D_S_bb_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_cd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_G_S_cd_coefs = gamma*delta;
      I_ERI_F3x_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G4x_S_vrr;
      I_ERI_F2xy_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G4x_S_vrr;
      I_ERI_F2xz_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G4x_S_vrr;
      I_ERI_Fx2y_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G4x_S_vrr;
      I_ERI_Fxyz_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G4x_S_vrr;
      I_ERI_Fx2z_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G4x_S_vrr;
      I_ERI_F3y_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G4x_S_vrr;
      I_ERI_F2yz_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G4x_S_vrr;
      I_ERI_Fy2z_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G4x_S_vrr;
      I_ERI_F3z_S_G4x_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G4x_S_vrr;
      I_ERI_F3x_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G3xy_S_vrr;
      I_ERI_F2xy_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G3xy_S_vrr;
      I_ERI_F2xz_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G3xy_S_vrr;
      I_ERI_Fx2y_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G3xy_S_vrr;
      I_ERI_Fxyz_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G3xy_S_vrr;
      I_ERI_Fx2z_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G3xy_S_vrr;
      I_ERI_F3y_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G3xy_S_vrr;
      I_ERI_F2yz_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G3xy_S_vrr;
      I_ERI_Fy2z_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G3xy_S_vrr;
      I_ERI_F3z_S_G3xy_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G3xy_S_vrr;
      I_ERI_F3x_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G3xz_S_vrr;
      I_ERI_F2xy_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G3xz_S_vrr;
      I_ERI_F2xz_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G3xz_S_vrr;
      I_ERI_Fx2y_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G3xz_S_vrr;
      I_ERI_Fxyz_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G3xz_S_vrr;
      I_ERI_Fx2z_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G3xz_S_vrr;
      I_ERI_F3y_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G3xz_S_vrr;
      I_ERI_F2yz_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G3xz_S_vrr;
      I_ERI_Fy2z_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G3xz_S_vrr;
      I_ERI_F3z_S_G3xz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G3xz_S_vrr;
      I_ERI_F3x_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G2x2y_S_vrr;
      I_ERI_F2xy_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G2x2y_S_vrr;
      I_ERI_F2xz_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G2x2y_S_vrr;
      I_ERI_Fx2y_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G2x2y_S_vrr;
      I_ERI_Fxyz_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G2x2y_S_vrr;
      I_ERI_Fx2z_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G2x2y_S_vrr;
      I_ERI_F3y_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G2x2y_S_vrr;
      I_ERI_F2yz_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G2x2y_S_vrr;
      I_ERI_Fy2z_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G2x2y_S_vrr;
      I_ERI_F3z_S_G2x2y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G2x2y_S_vrr;
      I_ERI_F3x_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G2xyz_S_vrr;
      I_ERI_F2xy_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G2xyz_S_vrr;
      I_ERI_F2xz_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G2xyz_S_vrr;
      I_ERI_Fx2y_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G2xyz_S_vrr;
      I_ERI_Fxyz_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G2xyz_S_vrr;
      I_ERI_Fx2z_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G2xyz_S_vrr;
      I_ERI_F3y_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G2xyz_S_vrr;
      I_ERI_F2yz_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G2xyz_S_vrr;
      I_ERI_Fy2z_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G2xyz_S_vrr;
      I_ERI_F3z_S_G2xyz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G2xyz_S_vrr;
      I_ERI_F3x_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G2x2z_S_vrr;
      I_ERI_F2xy_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G2x2z_S_vrr;
      I_ERI_F2xz_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G2x2z_S_vrr;
      I_ERI_Fx2y_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G2x2z_S_vrr;
      I_ERI_Fxyz_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G2x2z_S_vrr;
      I_ERI_Fx2z_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G2x2z_S_vrr;
      I_ERI_F3y_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G2x2z_S_vrr;
      I_ERI_F2yz_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G2x2z_S_vrr;
      I_ERI_Fy2z_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G2x2z_S_vrr;
      I_ERI_F3z_S_G2x2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G2x2z_S_vrr;
      I_ERI_F3x_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_Gx3y_S_vrr;
      I_ERI_F2xy_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_Gx3y_S_vrr;
      I_ERI_F2xz_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_Gx3y_S_vrr;
      I_ERI_Fx2y_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_Gx3y_S_vrr;
      I_ERI_Fxyz_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_Gx3y_S_vrr;
      I_ERI_Fx2z_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_Gx3y_S_vrr;
      I_ERI_F3y_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_Gx3y_S_vrr;
      I_ERI_F2yz_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_Gx3y_S_vrr;
      I_ERI_Fy2z_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_Gx3y_S_vrr;
      I_ERI_F3z_S_Gx3y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_Gx3y_S_vrr;
      I_ERI_F3x_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_Gx2yz_S_vrr;
      I_ERI_F2xy_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_Gx2yz_S_vrr;
      I_ERI_F2xz_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_Gx2yz_S_vrr;
      I_ERI_Fx2y_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_Gx2yz_S_vrr;
      I_ERI_Fxyz_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_Gx2yz_S_vrr;
      I_ERI_Fx2z_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_Gx2yz_S_vrr;
      I_ERI_F3y_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_Gx2yz_S_vrr;
      I_ERI_F2yz_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_Gx2yz_S_vrr;
      I_ERI_Fy2z_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_Gx2yz_S_vrr;
      I_ERI_F3z_S_Gx2yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_Gx2yz_S_vrr;
      I_ERI_F3x_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_Gxy2z_S_vrr;
      I_ERI_F2xy_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_Gxy2z_S_vrr;
      I_ERI_F2xz_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_Gxy2z_S_vrr;
      I_ERI_Fx2y_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_Gxy2z_S_vrr;
      I_ERI_Fxyz_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_Gxy2z_S_vrr;
      I_ERI_Fx2z_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_Gxy2z_S_vrr;
      I_ERI_F3y_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_Gxy2z_S_vrr;
      I_ERI_F2yz_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_Gxy2z_S_vrr;
      I_ERI_Fy2z_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_Gxy2z_S_vrr;
      I_ERI_F3z_S_Gxy2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_Gxy2z_S_vrr;
      I_ERI_F3x_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_Gx3z_S_vrr;
      I_ERI_F2xy_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_Gx3z_S_vrr;
      I_ERI_F2xz_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_Gx3z_S_vrr;
      I_ERI_Fx2y_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_Gx3z_S_vrr;
      I_ERI_Fxyz_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_Gx3z_S_vrr;
      I_ERI_Fx2z_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_Gx3z_S_vrr;
      I_ERI_F3y_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_Gx3z_S_vrr;
      I_ERI_F2yz_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_Gx3z_S_vrr;
      I_ERI_Fy2z_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_Gx3z_S_vrr;
      I_ERI_F3z_S_Gx3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_Gx3z_S_vrr;
      I_ERI_F3x_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G4y_S_vrr;
      I_ERI_F2xy_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G4y_S_vrr;
      I_ERI_F2xz_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G4y_S_vrr;
      I_ERI_Fx2y_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G4y_S_vrr;
      I_ERI_Fxyz_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G4y_S_vrr;
      I_ERI_Fx2z_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G4y_S_vrr;
      I_ERI_F3y_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G4y_S_vrr;
      I_ERI_F2yz_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G4y_S_vrr;
      I_ERI_Fy2z_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G4y_S_vrr;
      I_ERI_F3z_S_G4y_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G4y_S_vrr;
      I_ERI_F3x_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G3yz_S_vrr;
      I_ERI_F2xy_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G3yz_S_vrr;
      I_ERI_F2xz_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G3yz_S_vrr;
      I_ERI_Fx2y_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G3yz_S_vrr;
      I_ERI_Fxyz_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G3yz_S_vrr;
      I_ERI_Fx2z_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G3yz_S_vrr;
      I_ERI_F3y_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G3yz_S_vrr;
      I_ERI_F2yz_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G3yz_S_vrr;
      I_ERI_Fy2z_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G3yz_S_vrr;
      I_ERI_F3z_S_G3yz_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G3yz_S_vrr;
      I_ERI_F3x_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G2y2z_S_vrr;
      I_ERI_F2xy_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G2y2z_S_vrr;
      I_ERI_F2xz_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G2y2z_S_vrr;
      I_ERI_Fx2y_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G2y2z_S_vrr;
      I_ERI_Fxyz_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G2y2z_S_vrr;
      I_ERI_Fx2z_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G2y2z_S_vrr;
      I_ERI_F3y_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G2y2z_S_vrr;
      I_ERI_F2yz_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G2y2z_S_vrr;
      I_ERI_Fy2z_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G2y2z_S_vrr;
      I_ERI_F3z_S_G2y2z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G2y2z_S_vrr;
      I_ERI_F3x_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_Gy3z_S_vrr;
      I_ERI_F2xy_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_Gy3z_S_vrr;
      I_ERI_F2xz_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_Gy3z_S_vrr;
      I_ERI_Fx2y_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_Gy3z_S_vrr;
      I_ERI_Fxyz_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_Gy3z_S_vrr;
      I_ERI_Fx2z_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_Gy3z_S_vrr;
      I_ERI_F3y_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_Gy3z_S_vrr;
      I_ERI_F2yz_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_Gy3z_S_vrr;
      I_ERI_Fy2z_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_Gy3z_S_vrr;
      I_ERI_F3z_S_Gy3z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_Gy3z_S_vrr;
      I_ERI_F3x_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3x_S_G4z_S_vrr;
      I_ERI_F2xy_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xy_S_G4z_S_vrr;
      I_ERI_F2xz_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2xz_S_G4z_S_vrr;
      I_ERI_Fx2y_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2y_S_G4z_S_vrr;
      I_ERI_Fxyz_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fxyz_S_G4z_S_vrr;
      I_ERI_Fx2z_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fx2z_S_G4z_S_vrr;
      I_ERI_F3y_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3y_S_G4z_S_vrr;
      I_ERI_F2yz_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F2yz_S_G4z_S_vrr;
      I_ERI_Fy2z_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_Fy2z_S_G4z_S_vrr;
      I_ERI_F3z_S_G4z_S_cd += SQ_ERI_F_S_G_S_cd_coefs*I_ERI_F3z_S_G4z_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_F_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_bd_coefs = beta*delta;
      I_ERI_G4x_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_bd += SQ_ERI_G_S_F_S_bd_coefs*I_ERI_G4z_S_F3z_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_F_S_bd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_bd_coefs = beta*delta;
      I_ERI_F3x_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_bd += SQ_ERI_F_S_F_S_bd_coefs*I_ERI_F3z_S_F3z_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_G_S_dd
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_G_S_dd_coefs = delta*delta;
      I_ERI_F3x_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G4x_S_vrr;
      I_ERI_F2xy_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G4x_S_vrr;
      I_ERI_F2xz_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G4x_S_vrr;
      I_ERI_Fx2y_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G4x_S_vrr;
      I_ERI_Fxyz_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G4x_S_vrr;
      I_ERI_Fx2z_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G4x_S_vrr;
      I_ERI_F3y_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G4x_S_vrr;
      I_ERI_F2yz_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G4x_S_vrr;
      I_ERI_Fy2z_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G4x_S_vrr;
      I_ERI_F3z_S_G4x_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G4x_S_vrr;
      I_ERI_F3x_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G3xy_S_vrr;
      I_ERI_F2xy_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G3xy_S_vrr;
      I_ERI_F2xz_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G3xy_S_vrr;
      I_ERI_Fx2y_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G3xy_S_vrr;
      I_ERI_Fxyz_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G3xy_S_vrr;
      I_ERI_Fx2z_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G3xy_S_vrr;
      I_ERI_F3y_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G3xy_S_vrr;
      I_ERI_F2yz_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G3xy_S_vrr;
      I_ERI_Fy2z_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G3xy_S_vrr;
      I_ERI_F3z_S_G3xy_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G3xy_S_vrr;
      I_ERI_F3x_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G3xz_S_vrr;
      I_ERI_F2xy_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G3xz_S_vrr;
      I_ERI_F2xz_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G3xz_S_vrr;
      I_ERI_Fx2y_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G3xz_S_vrr;
      I_ERI_Fxyz_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G3xz_S_vrr;
      I_ERI_Fx2z_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G3xz_S_vrr;
      I_ERI_F3y_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G3xz_S_vrr;
      I_ERI_F2yz_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G3xz_S_vrr;
      I_ERI_Fy2z_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G3xz_S_vrr;
      I_ERI_F3z_S_G3xz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G3xz_S_vrr;
      I_ERI_F3x_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G2x2y_S_vrr;
      I_ERI_F2xy_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G2x2y_S_vrr;
      I_ERI_F2xz_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G2x2y_S_vrr;
      I_ERI_Fx2y_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G2x2y_S_vrr;
      I_ERI_Fxyz_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G2x2y_S_vrr;
      I_ERI_Fx2z_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G2x2y_S_vrr;
      I_ERI_F3y_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G2x2y_S_vrr;
      I_ERI_F2yz_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G2x2y_S_vrr;
      I_ERI_Fy2z_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G2x2y_S_vrr;
      I_ERI_F3z_S_G2x2y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G2x2y_S_vrr;
      I_ERI_F3x_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G2xyz_S_vrr;
      I_ERI_F2xy_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G2xyz_S_vrr;
      I_ERI_F2xz_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G2xyz_S_vrr;
      I_ERI_Fx2y_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G2xyz_S_vrr;
      I_ERI_Fxyz_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G2xyz_S_vrr;
      I_ERI_Fx2z_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G2xyz_S_vrr;
      I_ERI_F3y_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G2xyz_S_vrr;
      I_ERI_F2yz_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G2xyz_S_vrr;
      I_ERI_Fy2z_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G2xyz_S_vrr;
      I_ERI_F3z_S_G2xyz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G2xyz_S_vrr;
      I_ERI_F3x_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G2x2z_S_vrr;
      I_ERI_F2xy_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G2x2z_S_vrr;
      I_ERI_F2xz_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G2x2z_S_vrr;
      I_ERI_Fx2y_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G2x2z_S_vrr;
      I_ERI_Fxyz_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G2x2z_S_vrr;
      I_ERI_Fx2z_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G2x2z_S_vrr;
      I_ERI_F3y_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G2x2z_S_vrr;
      I_ERI_F2yz_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G2x2z_S_vrr;
      I_ERI_Fy2z_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G2x2z_S_vrr;
      I_ERI_F3z_S_G2x2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G2x2z_S_vrr;
      I_ERI_F3x_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_Gx3y_S_vrr;
      I_ERI_F2xy_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_Gx3y_S_vrr;
      I_ERI_F2xz_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_Gx3y_S_vrr;
      I_ERI_Fx2y_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_Gx3y_S_vrr;
      I_ERI_Fxyz_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_Gx3y_S_vrr;
      I_ERI_Fx2z_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_Gx3y_S_vrr;
      I_ERI_F3y_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_Gx3y_S_vrr;
      I_ERI_F2yz_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_Gx3y_S_vrr;
      I_ERI_Fy2z_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_Gx3y_S_vrr;
      I_ERI_F3z_S_Gx3y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_Gx3y_S_vrr;
      I_ERI_F3x_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_Gx2yz_S_vrr;
      I_ERI_F2xy_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_Gx2yz_S_vrr;
      I_ERI_F2xz_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_Gx2yz_S_vrr;
      I_ERI_Fx2y_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_Gx2yz_S_vrr;
      I_ERI_Fxyz_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_Gx2yz_S_vrr;
      I_ERI_Fx2z_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_Gx2yz_S_vrr;
      I_ERI_F3y_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_Gx2yz_S_vrr;
      I_ERI_F2yz_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_Gx2yz_S_vrr;
      I_ERI_Fy2z_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_Gx2yz_S_vrr;
      I_ERI_F3z_S_Gx2yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_Gx2yz_S_vrr;
      I_ERI_F3x_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_Gxy2z_S_vrr;
      I_ERI_F2xy_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_Gxy2z_S_vrr;
      I_ERI_F2xz_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_Gxy2z_S_vrr;
      I_ERI_Fx2y_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_Gxy2z_S_vrr;
      I_ERI_Fxyz_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_Gxy2z_S_vrr;
      I_ERI_Fx2z_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_Gxy2z_S_vrr;
      I_ERI_F3y_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_Gxy2z_S_vrr;
      I_ERI_F2yz_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_Gxy2z_S_vrr;
      I_ERI_Fy2z_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_Gxy2z_S_vrr;
      I_ERI_F3z_S_Gxy2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_Gxy2z_S_vrr;
      I_ERI_F3x_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_Gx3z_S_vrr;
      I_ERI_F2xy_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_Gx3z_S_vrr;
      I_ERI_F2xz_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_Gx3z_S_vrr;
      I_ERI_Fx2y_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_Gx3z_S_vrr;
      I_ERI_Fxyz_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_Gx3z_S_vrr;
      I_ERI_Fx2z_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_Gx3z_S_vrr;
      I_ERI_F3y_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_Gx3z_S_vrr;
      I_ERI_F2yz_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_Gx3z_S_vrr;
      I_ERI_Fy2z_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_Gx3z_S_vrr;
      I_ERI_F3z_S_Gx3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_Gx3z_S_vrr;
      I_ERI_F3x_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G4y_S_vrr;
      I_ERI_F2xy_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G4y_S_vrr;
      I_ERI_F2xz_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G4y_S_vrr;
      I_ERI_Fx2y_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G4y_S_vrr;
      I_ERI_Fxyz_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G4y_S_vrr;
      I_ERI_Fx2z_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G4y_S_vrr;
      I_ERI_F3y_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G4y_S_vrr;
      I_ERI_F2yz_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G4y_S_vrr;
      I_ERI_Fy2z_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G4y_S_vrr;
      I_ERI_F3z_S_G4y_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G4y_S_vrr;
      I_ERI_F3x_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G3yz_S_vrr;
      I_ERI_F2xy_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G3yz_S_vrr;
      I_ERI_F2xz_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G3yz_S_vrr;
      I_ERI_Fx2y_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G3yz_S_vrr;
      I_ERI_Fxyz_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G3yz_S_vrr;
      I_ERI_Fx2z_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G3yz_S_vrr;
      I_ERI_F3y_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G3yz_S_vrr;
      I_ERI_F2yz_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G3yz_S_vrr;
      I_ERI_Fy2z_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G3yz_S_vrr;
      I_ERI_F3z_S_G3yz_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G3yz_S_vrr;
      I_ERI_F3x_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G2y2z_S_vrr;
      I_ERI_F2xy_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G2y2z_S_vrr;
      I_ERI_F2xz_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G2y2z_S_vrr;
      I_ERI_Fx2y_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G2y2z_S_vrr;
      I_ERI_Fxyz_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G2y2z_S_vrr;
      I_ERI_Fx2z_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G2y2z_S_vrr;
      I_ERI_F3y_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G2y2z_S_vrr;
      I_ERI_F2yz_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G2y2z_S_vrr;
      I_ERI_Fy2z_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G2y2z_S_vrr;
      I_ERI_F3z_S_G2y2z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G2y2z_S_vrr;
      I_ERI_F3x_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_Gy3z_S_vrr;
      I_ERI_F2xy_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_Gy3z_S_vrr;
      I_ERI_F2xz_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_Gy3z_S_vrr;
      I_ERI_Fx2y_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_Gy3z_S_vrr;
      I_ERI_Fxyz_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_Gy3z_S_vrr;
      I_ERI_Fx2z_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_Gy3z_S_vrr;
      I_ERI_F3y_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_Gy3z_S_vrr;
      I_ERI_F2yz_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_Gy3z_S_vrr;
      I_ERI_Fy2z_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_Gy3z_S_vrr;
      I_ERI_F3z_S_Gy3z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_Gy3z_S_vrr;
      I_ERI_F3x_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3x_S_G4z_S_vrr;
      I_ERI_F2xy_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xy_S_G4z_S_vrr;
      I_ERI_F2xz_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2xz_S_G4z_S_vrr;
      I_ERI_Fx2y_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2y_S_G4z_S_vrr;
      I_ERI_Fxyz_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fxyz_S_G4z_S_vrr;
      I_ERI_Fx2z_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fx2z_S_G4z_S_vrr;
      I_ERI_F3y_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3y_S_G4z_S_vrr;
      I_ERI_F2yz_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F2yz_S_G4z_S_vrr;
      I_ERI_Fy2z_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_Fy2z_S_G4z_S_vrr;
      I_ERI_F3z_S_G4z_S_dd += SQ_ERI_F_S_G_S_dd_coefs*I_ERI_F3z_S_G4z_S_vrr;

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
   * shell quartet name: SQ_ERI_F_S_P_P_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   * RHS shell quartet name: SQ_ERI_F_S_P_S_d
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px_d = I_ERI_F3x_S_D2x_S_d+CDX*I_ERI_F3x_S_Px_S_d;
  Double I_ERI_F2xy_S_Px_Px_d = I_ERI_F2xy_S_D2x_S_d+CDX*I_ERI_F2xy_S_Px_S_d;
  Double I_ERI_F2xz_S_Px_Px_d = I_ERI_F2xz_S_D2x_S_d+CDX*I_ERI_F2xz_S_Px_S_d;
  Double I_ERI_Fx2y_S_Px_Px_d = I_ERI_Fx2y_S_D2x_S_d+CDX*I_ERI_Fx2y_S_Px_S_d;
  Double I_ERI_Fxyz_S_Px_Px_d = I_ERI_Fxyz_S_D2x_S_d+CDX*I_ERI_Fxyz_S_Px_S_d;
  Double I_ERI_Fx2z_S_Px_Px_d = I_ERI_Fx2z_S_D2x_S_d+CDX*I_ERI_Fx2z_S_Px_S_d;
  Double I_ERI_F3y_S_Px_Px_d = I_ERI_F3y_S_D2x_S_d+CDX*I_ERI_F3y_S_Px_S_d;
  Double I_ERI_F2yz_S_Px_Px_d = I_ERI_F2yz_S_D2x_S_d+CDX*I_ERI_F2yz_S_Px_S_d;
  Double I_ERI_Fy2z_S_Px_Px_d = I_ERI_Fy2z_S_D2x_S_d+CDX*I_ERI_Fy2z_S_Px_S_d;
  Double I_ERI_F3z_S_Px_Px_d = I_ERI_F3z_S_D2x_S_d+CDX*I_ERI_F3z_S_Px_S_d;
  Double I_ERI_F3x_S_Py_Px_d = I_ERI_F3x_S_Dxy_S_d+CDX*I_ERI_F3x_S_Py_S_d;
  Double I_ERI_F2xy_S_Py_Px_d = I_ERI_F2xy_S_Dxy_S_d+CDX*I_ERI_F2xy_S_Py_S_d;
  Double I_ERI_F2xz_S_Py_Px_d = I_ERI_F2xz_S_Dxy_S_d+CDX*I_ERI_F2xz_S_Py_S_d;
  Double I_ERI_Fx2y_S_Py_Px_d = I_ERI_Fx2y_S_Dxy_S_d+CDX*I_ERI_Fx2y_S_Py_S_d;
  Double I_ERI_Fxyz_S_Py_Px_d = I_ERI_Fxyz_S_Dxy_S_d+CDX*I_ERI_Fxyz_S_Py_S_d;
  Double I_ERI_Fx2z_S_Py_Px_d = I_ERI_Fx2z_S_Dxy_S_d+CDX*I_ERI_Fx2z_S_Py_S_d;
  Double I_ERI_F3y_S_Py_Px_d = I_ERI_F3y_S_Dxy_S_d+CDX*I_ERI_F3y_S_Py_S_d;
  Double I_ERI_F2yz_S_Py_Px_d = I_ERI_F2yz_S_Dxy_S_d+CDX*I_ERI_F2yz_S_Py_S_d;
  Double I_ERI_Fy2z_S_Py_Px_d = I_ERI_Fy2z_S_Dxy_S_d+CDX*I_ERI_Fy2z_S_Py_S_d;
  Double I_ERI_F3z_S_Py_Px_d = I_ERI_F3z_S_Dxy_S_d+CDX*I_ERI_F3z_S_Py_S_d;
  Double I_ERI_F3x_S_Pz_Px_d = I_ERI_F3x_S_Dxz_S_d+CDX*I_ERI_F3x_S_Pz_S_d;
  Double I_ERI_F2xy_S_Pz_Px_d = I_ERI_F2xy_S_Dxz_S_d+CDX*I_ERI_F2xy_S_Pz_S_d;
  Double I_ERI_F2xz_S_Pz_Px_d = I_ERI_F2xz_S_Dxz_S_d+CDX*I_ERI_F2xz_S_Pz_S_d;
  Double I_ERI_Fx2y_S_Pz_Px_d = I_ERI_Fx2y_S_Dxz_S_d+CDX*I_ERI_Fx2y_S_Pz_S_d;
  Double I_ERI_Fxyz_S_Pz_Px_d = I_ERI_Fxyz_S_Dxz_S_d+CDX*I_ERI_Fxyz_S_Pz_S_d;
  Double I_ERI_Fx2z_S_Pz_Px_d = I_ERI_Fx2z_S_Dxz_S_d+CDX*I_ERI_Fx2z_S_Pz_S_d;
  Double I_ERI_F3y_S_Pz_Px_d = I_ERI_F3y_S_Dxz_S_d+CDX*I_ERI_F3y_S_Pz_S_d;
  Double I_ERI_F2yz_S_Pz_Px_d = I_ERI_F2yz_S_Dxz_S_d+CDX*I_ERI_F2yz_S_Pz_S_d;
  Double I_ERI_Fy2z_S_Pz_Px_d = I_ERI_Fy2z_S_Dxz_S_d+CDX*I_ERI_Fy2z_S_Pz_S_d;
  Double I_ERI_F3z_S_Pz_Px_d = I_ERI_F3z_S_Dxz_S_d+CDX*I_ERI_F3z_S_Pz_S_d;
  Double I_ERI_F3x_S_Px_Py_d = I_ERI_F3x_S_Dxy_S_d+CDY*I_ERI_F3x_S_Px_S_d;
  Double I_ERI_F2xy_S_Px_Py_d = I_ERI_F2xy_S_Dxy_S_d+CDY*I_ERI_F2xy_S_Px_S_d;
  Double I_ERI_F2xz_S_Px_Py_d = I_ERI_F2xz_S_Dxy_S_d+CDY*I_ERI_F2xz_S_Px_S_d;
  Double I_ERI_Fx2y_S_Px_Py_d = I_ERI_Fx2y_S_Dxy_S_d+CDY*I_ERI_Fx2y_S_Px_S_d;
  Double I_ERI_Fxyz_S_Px_Py_d = I_ERI_Fxyz_S_Dxy_S_d+CDY*I_ERI_Fxyz_S_Px_S_d;
  Double I_ERI_Fx2z_S_Px_Py_d = I_ERI_Fx2z_S_Dxy_S_d+CDY*I_ERI_Fx2z_S_Px_S_d;
  Double I_ERI_F3y_S_Px_Py_d = I_ERI_F3y_S_Dxy_S_d+CDY*I_ERI_F3y_S_Px_S_d;
  Double I_ERI_F2yz_S_Px_Py_d = I_ERI_F2yz_S_Dxy_S_d+CDY*I_ERI_F2yz_S_Px_S_d;
  Double I_ERI_Fy2z_S_Px_Py_d = I_ERI_Fy2z_S_Dxy_S_d+CDY*I_ERI_Fy2z_S_Px_S_d;
  Double I_ERI_F3z_S_Px_Py_d = I_ERI_F3z_S_Dxy_S_d+CDY*I_ERI_F3z_S_Px_S_d;
  Double I_ERI_F3x_S_Py_Py_d = I_ERI_F3x_S_D2y_S_d+CDY*I_ERI_F3x_S_Py_S_d;
  Double I_ERI_F2xy_S_Py_Py_d = I_ERI_F2xy_S_D2y_S_d+CDY*I_ERI_F2xy_S_Py_S_d;
  Double I_ERI_F2xz_S_Py_Py_d = I_ERI_F2xz_S_D2y_S_d+CDY*I_ERI_F2xz_S_Py_S_d;
  Double I_ERI_Fx2y_S_Py_Py_d = I_ERI_Fx2y_S_D2y_S_d+CDY*I_ERI_Fx2y_S_Py_S_d;
  Double I_ERI_Fxyz_S_Py_Py_d = I_ERI_Fxyz_S_D2y_S_d+CDY*I_ERI_Fxyz_S_Py_S_d;
  Double I_ERI_Fx2z_S_Py_Py_d = I_ERI_Fx2z_S_D2y_S_d+CDY*I_ERI_Fx2z_S_Py_S_d;
  Double I_ERI_F3y_S_Py_Py_d = I_ERI_F3y_S_D2y_S_d+CDY*I_ERI_F3y_S_Py_S_d;
  Double I_ERI_F2yz_S_Py_Py_d = I_ERI_F2yz_S_D2y_S_d+CDY*I_ERI_F2yz_S_Py_S_d;
  Double I_ERI_Fy2z_S_Py_Py_d = I_ERI_Fy2z_S_D2y_S_d+CDY*I_ERI_Fy2z_S_Py_S_d;
  Double I_ERI_F3z_S_Py_Py_d = I_ERI_F3z_S_D2y_S_d+CDY*I_ERI_F3z_S_Py_S_d;
  Double I_ERI_F3x_S_Pz_Py_d = I_ERI_F3x_S_Dyz_S_d+CDY*I_ERI_F3x_S_Pz_S_d;
  Double I_ERI_F2xy_S_Pz_Py_d = I_ERI_F2xy_S_Dyz_S_d+CDY*I_ERI_F2xy_S_Pz_S_d;
  Double I_ERI_F2xz_S_Pz_Py_d = I_ERI_F2xz_S_Dyz_S_d+CDY*I_ERI_F2xz_S_Pz_S_d;
  Double I_ERI_Fx2y_S_Pz_Py_d = I_ERI_Fx2y_S_Dyz_S_d+CDY*I_ERI_Fx2y_S_Pz_S_d;
  Double I_ERI_Fxyz_S_Pz_Py_d = I_ERI_Fxyz_S_Dyz_S_d+CDY*I_ERI_Fxyz_S_Pz_S_d;
  Double I_ERI_Fx2z_S_Pz_Py_d = I_ERI_Fx2z_S_Dyz_S_d+CDY*I_ERI_Fx2z_S_Pz_S_d;
  Double I_ERI_F3y_S_Pz_Py_d = I_ERI_F3y_S_Dyz_S_d+CDY*I_ERI_F3y_S_Pz_S_d;
  Double I_ERI_F2yz_S_Pz_Py_d = I_ERI_F2yz_S_Dyz_S_d+CDY*I_ERI_F2yz_S_Pz_S_d;
  Double I_ERI_Fy2z_S_Pz_Py_d = I_ERI_Fy2z_S_Dyz_S_d+CDY*I_ERI_Fy2z_S_Pz_S_d;
  Double I_ERI_F3z_S_Pz_Py_d = I_ERI_F3z_S_Dyz_S_d+CDY*I_ERI_F3z_S_Pz_S_d;
  Double I_ERI_F3x_S_Px_Pz_d = I_ERI_F3x_S_Dxz_S_d+CDZ*I_ERI_F3x_S_Px_S_d;
  Double I_ERI_F2xy_S_Px_Pz_d = I_ERI_F2xy_S_Dxz_S_d+CDZ*I_ERI_F2xy_S_Px_S_d;
  Double I_ERI_F2xz_S_Px_Pz_d = I_ERI_F2xz_S_Dxz_S_d+CDZ*I_ERI_F2xz_S_Px_S_d;
  Double I_ERI_Fx2y_S_Px_Pz_d = I_ERI_Fx2y_S_Dxz_S_d+CDZ*I_ERI_Fx2y_S_Px_S_d;
  Double I_ERI_Fxyz_S_Px_Pz_d = I_ERI_Fxyz_S_Dxz_S_d+CDZ*I_ERI_Fxyz_S_Px_S_d;
  Double I_ERI_Fx2z_S_Px_Pz_d = I_ERI_Fx2z_S_Dxz_S_d+CDZ*I_ERI_Fx2z_S_Px_S_d;
  Double I_ERI_F3y_S_Px_Pz_d = I_ERI_F3y_S_Dxz_S_d+CDZ*I_ERI_F3y_S_Px_S_d;
  Double I_ERI_F2yz_S_Px_Pz_d = I_ERI_F2yz_S_Dxz_S_d+CDZ*I_ERI_F2yz_S_Px_S_d;
  Double I_ERI_Fy2z_S_Px_Pz_d = I_ERI_Fy2z_S_Dxz_S_d+CDZ*I_ERI_Fy2z_S_Px_S_d;
  Double I_ERI_F3z_S_Px_Pz_d = I_ERI_F3z_S_Dxz_S_d+CDZ*I_ERI_F3z_S_Px_S_d;
  Double I_ERI_F3x_S_Py_Pz_d = I_ERI_F3x_S_Dyz_S_d+CDZ*I_ERI_F3x_S_Py_S_d;
  Double I_ERI_F2xy_S_Py_Pz_d = I_ERI_F2xy_S_Dyz_S_d+CDZ*I_ERI_F2xy_S_Py_S_d;
  Double I_ERI_F2xz_S_Py_Pz_d = I_ERI_F2xz_S_Dyz_S_d+CDZ*I_ERI_F2xz_S_Py_S_d;
  Double I_ERI_Fx2y_S_Py_Pz_d = I_ERI_Fx2y_S_Dyz_S_d+CDZ*I_ERI_Fx2y_S_Py_S_d;
  Double I_ERI_Fxyz_S_Py_Pz_d = I_ERI_Fxyz_S_Dyz_S_d+CDZ*I_ERI_Fxyz_S_Py_S_d;
  Double I_ERI_Fx2z_S_Py_Pz_d = I_ERI_Fx2z_S_Dyz_S_d+CDZ*I_ERI_Fx2z_S_Py_S_d;
  Double I_ERI_F3y_S_Py_Pz_d = I_ERI_F3y_S_Dyz_S_d+CDZ*I_ERI_F3y_S_Py_S_d;
  Double I_ERI_F2yz_S_Py_Pz_d = I_ERI_F2yz_S_Dyz_S_d+CDZ*I_ERI_F2yz_S_Py_S_d;
  Double I_ERI_Fy2z_S_Py_Pz_d = I_ERI_Fy2z_S_Dyz_S_d+CDZ*I_ERI_Fy2z_S_Py_S_d;
  Double I_ERI_F3z_S_Py_Pz_d = I_ERI_F3z_S_Dyz_S_d+CDZ*I_ERI_F3z_S_Py_S_d;
  Double I_ERI_F3x_S_Pz_Pz_d = I_ERI_F3x_S_D2z_S_d+CDZ*I_ERI_F3x_S_Pz_S_d;
  Double I_ERI_F2xy_S_Pz_Pz_d = I_ERI_F2xy_S_D2z_S_d+CDZ*I_ERI_F2xy_S_Pz_S_d;
  Double I_ERI_F2xz_S_Pz_Pz_d = I_ERI_F2xz_S_D2z_S_d+CDZ*I_ERI_F2xz_S_Pz_S_d;
  Double I_ERI_Fx2y_S_Pz_Pz_d = I_ERI_Fx2y_S_D2z_S_d+CDZ*I_ERI_Fx2y_S_Pz_S_d;
  Double I_ERI_Fxyz_S_Pz_Pz_d = I_ERI_Fxyz_S_D2z_S_d+CDZ*I_ERI_Fxyz_S_Pz_S_d;
  Double I_ERI_Fx2z_S_Pz_Pz_d = I_ERI_Fx2z_S_D2z_S_d+CDZ*I_ERI_Fx2z_S_Pz_S_d;
  Double I_ERI_F3y_S_Pz_Pz_d = I_ERI_F3y_S_D2z_S_d+CDZ*I_ERI_F3y_S_Pz_S_d;
  Double I_ERI_F2yz_S_Pz_Pz_d = I_ERI_F2yz_S_D2z_S_d+CDZ*I_ERI_F2yz_S_Pz_S_d;
  Double I_ERI_Fy2z_S_Pz_Pz_d = I_ERI_Fy2z_S_D2z_S_d+CDZ*I_ERI_Fy2z_S_Pz_S_d;
  Double I_ERI_F3z_S_Pz_Pz_d = I_ERI_F3z_S_D2z_S_d+CDZ*I_ERI_F3z_S_Pz_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_bd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_bd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bd
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_bd = I_ERI_F3x_S_F3x_S_bd+CDX*I_ERI_F3x_S_D2x_S_bd;
  Double I_ERI_F2xy_S_D2x_Px_bd = I_ERI_F2xy_S_F3x_S_bd+CDX*I_ERI_F2xy_S_D2x_S_bd;
  Double I_ERI_F2xz_S_D2x_Px_bd = I_ERI_F2xz_S_F3x_S_bd+CDX*I_ERI_F2xz_S_D2x_S_bd;
  Double I_ERI_Fx2y_S_D2x_Px_bd = I_ERI_Fx2y_S_F3x_S_bd+CDX*I_ERI_Fx2y_S_D2x_S_bd;
  Double I_ERI_Fxyz_S_D2x_Px_bd = I_ERI_Fxyz_S_F3x_S_bd+CDX*I_ERI_Fxyz_S_D2x_S_bd;
  Double I_ERI_Fx2z_S_D2x_Px_bd = I_ERI_Fx2z_S_F3x_S_bd+CDX*I_ERI_Fx2z_S_D2x_S_bd;
  Double I_ERI_F3y_S_D2x_Px_bd = I_ERI_F3y_S_F3x_S_bd+CDX*I_ERI_F3y_S_D2x_S_bd;
  Double I_ERI_F2yz_S_D2x_Px_bd = I_ERI_F2yz_S_F3x_S_bd+CDX*I_ERI_F2yz_S_D2x_S_bd;
  Double I_ERI_Fy2z_S_D2x_Px_bd = I_ERI_Fy2z_S_F3x_S_bd+CDX*I_ERI_Fy2z_S_D2x_S_bd;
  Double I_ERI_F3z_S_D2x_Px_bd = I_ERI_F3z_S_F3x_S_bd+CDX*I_ERI_F3z_S_D2x_S_bd;
  Double I_ERI_F3x_S_Dxy_Px_bd = I_ERI_F3x_S_F2xy_S_bd+CDX*I_ERI_F3x_S_Dxy_S_bd;
  Double I_ERI_F2xy_S_Dxy_Px_bd = I_ERI_F2xy_S_F2xy_S_bd+CDX*I_ERI_F2xy_S_Dxy_S_bd;
  Double I_ERI_F2xz_S_Dxy_Px_bd = I_ERI_F2xz_S_F2xy_S_bd+CDX*I_ERI_F2xz_S_Dxy_S_bd;
  Double I_ERI_Fx2y_S_Dxy_Px_bd = I_ERI_Fx2y_S_F2xy_S_bd+CDX*I_ERI_Fx2y_S_Dxy_S_bd;
  Double I_ERI_Fxyz_S_Dxy_Px_bd = I_ERI_Fxyz_S_F2xy_S_bd+CDX*I_ERI_Fxyz_S_Dxy_S_bd;
  Double I_ERI_Fx2z_S_Dxy_Px_bd = I_ERI_Fx2z_S_F2xy_S_bd+CDX*I_ERI_Fx2z_S_Dxy_S_bd;
  Double I_ERI_F3y_S_Dxy_Px_bd = I_ERI_F3y_S_F2xy_S_bd+CDX*I_ERI_F3y_S_Dxy_S_bd;
  Double I_ERI_F2yz_S_Dxy_Px_bd = I_ERI_F2yz_S_F2xy_S_bd+CDX*I_ERI_F2yz_S_Dxy_S_bd;
  Double I_ERI_Fy2z_S_Dxy_Px_bd = I_ERI_Fy2z_S_F2xy_S_bd+CDX*I_ERI_Fy2z_S_Dxy_S_bd;
  Double I_ERI_F3z_S_Dxy_Px_bd = I_ERI_F3z_S_F2xy_S_bd+CDX*I_ERI_F3z_S_Dxy_S_bd;
  Double I_ERI_F3x_S_Dxz_Px_bd = I_ERI_F3x_S_F2xz_S_bd+CDX*I_ERI_F3x_S_Dxz_S_bd;
  Double I_ERI_F2xy_S_Dxz_Px_bd = I_ERI_F2xy_S_F2xz_S_bd+CDX*I_ERI_F2xy_S_Dxz_S_bd;
  Double I_ERI_F2xz_S_Dxz_Px_bd = I_ERI_F2xz_S_F2xz_S_bd+CDX*I_ERI_F2xz_S_Dxz_S_bd;
  Double I_ERI_Fx2y_S_Dxz_Px_bd = I_ERI_Fx2y_S_F2xz_S_bd+CDX*I_ERI_Fx2y_S_Dxz_S_bd;
  Double I_ERI_Fxyz_S_Dxz_Px_bd = I_ERI_Fxyz_S_F2xz_S_bd+CDX*I_ERI_Fxyz_S_Dxz_S_bd;
  Double I_ERI_Fx2z_S_Dxz_Px_bd = I_ERI_Fx2z_S_F2xz_S_bd+CDX*I_ERI_Fx2z_S_Dxz_S_bd;
  Double I_ERI_F3y_S_Dxz_Px_bd = I_ERI_F3y_S_F2xz_S_bd+CDX*I_ERI_F3y_S_Dxz_S_bd;
  Double I_ERI_F2yz_S_Dxz_Px_bd = I_ERI_F2yz_S_F2xz_S_bd+CDX*I_ERI_F2yz_S_Dxz_S_bd;
  Double I_ERI_Fy2z_S_Dxz_Px_bd = I_ERI_Fy2z_S_F2xz_S_bd+CDX*I_ERI_Fy2z_S_Dxz_S_bd;
  Double I_ERI_F3z_S_Dxz_Px_bd = I_ERI_F3z_S_F2xz_S_bd+CDX*I_ERI_F3z_S_Dxz_S_bd;
  Double I_ERI_F3x_S_D2y_Px_bd = I_ERI_F3x_S_Fx2y_S_bd+CDX*I_ERI_F3x_S_D2y_S_bd;
  Double I_ERI_F2xy_S_D2y_Px_bd = I_ERI_F2xy_S_Fx2y_S_bd+CDX*I_ERI_F2xy_S_D2y_S_bd;
  Double I_ERI_F2xz_S_D2y_Px_bd = I_ERI_F2xz_S_Fx2y_S_bd+CDX*I_ERI_F2xz_S_D2y_S_bd;
  Double I_ERI_Fx2y_S_D2y_Px_bd = I_ERI_Fx2y_S_Fx2y_S_bd+CDX*I_ERI_Fx2y_S_D2y_S_bd;
  Double I_ERI_Fxyz_S_D2y_Px_bd = I_ERI_Fxyz_S_Fx2y_S_bd+CDX*I_ERI_Fxyz_S_D2y_S_bd;
  Double I_ERI_Fx2z_S_D2y_Px_bd = I_ERI_Fx2z_S_Fx2y_S_bd+CDX*I_ERI_Fx2z_S_D2y_S_bd;
  Double I_ERI_F3y_S_D2y_Px_bd = I_ERI_F3y_S_Fx2y_S_bd+CDX*I_ERI_F3y_S_D2y_S_bd;
  Double I_ERI_F2yz_S_D2y_Px_bd = I_ERI_F2yz_S_Fx2y_S_bd+CDX*I_ERI_F2yz_S_D2y_S_bd;
  Double I_ERI_Fy2z_S_D2y_Px_bd = I_ERI_Fy2z_S_Fx2y_S_bd+CDX*I_ERI_Fy2z_S_D2y_S_bd;
  Double I_ERI_F3z_S_D2y_Px_bd = I_ERI_F3z_S_Fx2y_S_bd+CDX*I_ERI_F3z_S_D2y_S_bd;
  Double I_ERI_F3x_S_Dyz_Px_bd = I_ERI_F3x_S_Fxyz_S_bd+CDX*I_ERI_F3x_S_Dyz_S_bd;
  Double I_ERI_F2xy_S_Dyz_Px_bd = I_ERI_F2xy_S_Fxyz_S_bd+CDX*I_ERI_F2xy_S_Dyz_S_bd;
  Double I_ERI_F2xz_S_Dyz_Px_bd = I_ERI_F2xz_S_Fxyz_S_bd+CDX*I_ERI_F2xz_S_Dyz_S_bd;
  Double I_ERI_Fx2y_S_Dyz_Px_bd = I_ERI_Fx2y_S_Fxyz_S_bd+CDX*I_ERI_Fx2y_S_Dyz_S_bd;
  Double I_ERI_Fxyz_S_Dyz_Px_bd = I_ERI_Fxyz_S_Fxyz_S_bd+CDX*I_ERI_Fxyz_S_Dyz_S_bd;
  Double I_ERI_Fx2z_S_Dyz_Px_bd = I_ERI_Fx2z_S_Fxyz_S_bd+CDX*I_ERI_Fx2z_S_Dyz_S_bd;
  Double I_ERI_F3y_S_Dyz_Px_bd = I_ERI_F3y_S_Fxyz_S_bd+CDX*I_ERI_F3y_S_Dyz_S_bd;
  Double I_ERI_F2yz_S_Dyz_Px_bd = I_ERI_F2yz_S_Fxyz_S_bd+CDX*I_ERI_F2yz_S_Dyz_S_bd;
  Double I_ERI_Fy2z_S_Dyz_Px_bd = I_ERI_Fy2z_S_Fxyz_S_bd+CDX*I_ERI_Fy2z_S_Dyz_S_bd;
  Double I_ERI_F3z_S_Dyz_Px_bd = I_ERI_F3z_S_Fxyz_S_bd+CDX*I_ERI_F3z_S_Dyz_S_bd;
  Double I_ERI_F3x_S_D2z_Px_bd = I_ERI_F3x_S_Fx2z_S_bd+CDX*I_ERI_F3x_S_D2z_S_bd;
  Double I_ERI_F2xy_S_D2z_Px_bd = I_ERI_F2xy_S_Fx2z_S_bd+CDX*I_ERI_F2xy_S_D2z_S_bd;
  Double I_ERI_F2xz_S_D2z_Px_bd = I_ERI_F2xz_S_Fx2z_S_bd+CDX*I_ERI_F2xz_S_D2z_S_bd;
  Double I_ERI_Fx2y_S_D2z_Px_bd = I_ERI_Fx2y_S_Fx2z_S_bd+CDX*I_ERI_Fx2y_S_D2z_S_bd;
  Double I_ERI_Fxyz_S_D2z_Px_bd = I_ERI_Fxyz_S_Fx2z_S_bd+CDX*I_ERI_Fxyz_S_D2z_S_bd;
  Double I_ERI_Fx2z_S_D2z_Px_bd = I_ERI_Fx2z_S_Fx2z_S_bd+CDX*I_ERI_Fx2z_S_D2z_S_bd;
  Double I_ERI_F3y_S_D2z_Px_bd = I_ERI_F3y_S_Fx2z_S_bd+CDX*I_ERI_F3y_S_D2z_S_bd;
  Double I_ERI_F2yz_S_D2z_Px_bd = I_ERI_F2yz_S_Fx2z_S_bd+CDX*I_ERI_F2yz_S_D2z_S_bd;
  Double I_ERI_Fy2z_S_D2z_Px_bd = I_ERI_Fy2z_S_Fx2z_S_bd+CDX*I_ERI_Fy2z_S_D2z_S_bd;
  Double I_ERI_F3z_S_D2z_Px_bd = I_ERI_F3z_S_Fx2z_S_bd+CDX*I_ERI_F3z_S_D2z_S_bd;
  Double I_ERI_F3x_S_D2x_Py_bd = I_ERI_F3x_S_F2xy_S_bd+CDY*I_ERI_F3x_S_D2x_S_bd;
  Double I_ERI_F2xy_S_D2x_Py_bd = I_ERI_F2xy_S_F2xy_S_bd+CDY*I_ERI_F2xy_S_D2x_S_bd;
  Double I_ERI_F2xz_S_D2x_Py_bd = I_ERI_F2xz_S_F2xy_S_bd+CDY*I_ERI_F2xz_S_D2x_S_bd;
  Double I_ERI_Fx2y_S_D2x_Py_bd = I_ERI_Fx2y_S_F2xy_S_bd+CDY*I_ERI_Fx2y_S_D2x_S_bd;
  Double I_ERI_Fxyz_S_D2x_Py_bd = I_ERI_Fxyz_S_F2xy_S_bd+CDY*I_ERI_Fxyz_S_D2x_S_bd;
  Double I_ERI_Fx2z_S_D2x_Py_bd = I_ERI_Fx2z_S_F2xy_S_bd+CDY*I_ERI_Fx2z_S_D2x_S_bd;
  Double I_ERI_F3y_S_D2x_Py_bd = I_ERI_F3y_S_F2xy_S_bd+CDY*I_ERI_F3y_S_D2x_S_bd;
  Double I_ERI_F2yz_S_D2x_Py_bd = I_ERI_F2yz_S_F2xy_S_bd+CDY*I_ERI_F2yz_S_D2x_S_bd;
  Double I_ERI_Fy2z_S_D2x_Py_bd = I_ERI_Fy2z_S_F2xy_S_bd+CDY*I_ERI_Fy2z_S_D2x_S_bd;
  Double I_ERI_F3z_S_D2x_Py_bd = I_ERI_F3z_S_F2xy_S_bd+CDY*I_ERI_F3z_S_D2x_S_bd;
  Double I_ERI_F3x_S_Dxy_Py_bd = I_ERI_F3x_S_Fx2y_S_bd+CDY*I_ERI_F3x_S_Dxy_S_bd;
  Double I_ERI_F2xy_S_Dxy_Py_bd = I_ERI_F2xy_S_Fx2y_S_bd+CDY*I_ERI_F2xy_S_Dxy_S_bd;
  Double I_ERI_F2xz_S_Dxy_Py_bd = I_ERI_F2xz_S_Fx2y_S_bd+CDY*I_ERI_F2xz_S_Dxy_S_bd;
  Double I_ERI_Fx2y_S_Dxy_Py_bd = I_ERI_Fx2y_S_Fx2y_S_bd+CDY*I_ERI_Fx2y_S_Dxy_S_bd;
  Double I_ERI_Fxyz_S_Dxy_Py_bd = I_ERI_Fxyz_S_Fx2y_S_bd+CDY*I_ERI_Fxyz_S_Dxy_S_bd;
  Double I_ERI_Fx2z_S_Dxy_Py_bd = I_ERI_Fx2z_S_Fx2y_S_bd+CDY*I_ERI_Fx2z_S_Dxy_S_bd;
  Double I_ERI_F3y_S_Dxy_Py_bd = I_ERI_F3y_S_Fx2y_S_bd+CDY*I_ERI_F3y_S_Dxy_S_bd;
  Double I_ERI_F2yz_S_Dxy_Py_bd = I_ERI_F2yz_S_Fx2y_S_bd+CDY*I_ERI_F2yz_S_Dxy_S_bd;
  Double I_ERI_Fy2z_S_Dxy_Py_bd = I_ERI_Fy2z_S_Fx2y_S_bd+CDY*I_ERI_Fy2z_S_Dxy_S_bd;
  Double I_ERI_F3z_S_Dxy_Py_bd = I_ERI_F3z_S_Fx2y_S_bd+CDY*I_ERI_F3z_S_Dxy_S_bd;
  Double I_ERI_F3x_S_Dxz_Py_bd = I_ERI_F3x_S_Fxyz_S_bd+CDY*I_ERI_F3x_S_Dxz_S_bd;
  Double I_ERI_F2xy_S_Dxz_Py_bd = I_ERI_F2xy_S_Fxyz_S_bd+CDY*I_ERI_F2xy_S_Dxz_S_bd;
  Double I_ERI_F2xz_S_Dxz_Py_bd = I_ERI_F2xz_S_Fxyz_S_bd+CDY*I_ERI_F2xz_S_Dxz_S_bd;
  Double I_ERI_Fx2y_S_Dxz_Py_bd = I_ERI_Fx2y_S_Fxyz_S_bd+CDY*I_ERI_Fx2y_S_Dxz_S_bd;
  Double I_ERI_Fxyz_S_Dxz_Py_bd = I_ERI_Fxyz_S_Fxyz_S_bd+CDY*I_ERI_Fxyz_S_Dxz_S_bd;
  Double I_ERI_Fx2z_S_Dxz_Py_bd = I_ERI_Fx2z_S_Fxyz_S_bd+CDY*I_ERI_Fx2z_S_Dxz_S_bd;
  Double I_ERI_F3y_S_Dxz_Py_bd = I_ERI_F3y_S_Fxyz_S_bd+CDY*I_ERI_F3y_S_Dxz_S_bd;
  Double I_ERI_F2yz_S_Dxz_Py_bd = I_ERI_F2yz_S_Fxyz_S_bd+CDY*I_ERI_F2yz_S_Dxz_S_bd;
  Double I_ERI_Fy2z_S_Dxz_Py_bd = I_ERI_Fy2z_S_Fxyz_S_bd+CDY*I_ERI_Fy2z_S_Dxz_S_bd;
  Double I_ERI_F3z_S_Dxz_Py_bd = I_ERI_F3z_S_Fxyz_S_bd+CDY*I_ERI_F3z_S_Dxz_S_bd;
  Double I_ERI_F3x_S_D2y_Py_bd = I_ERI_F3x_S_F3y_S_bd+CDY*I_ERI_F3x_S_D2y_S_bd;
  Double I_ERI_F2xy_S_D2y_Py_bd = I_ERI_F2xy_S_F3y_S_bd+CDY*I_ERI_F2xy_S_D2y_S_bd;
  Double I_ERI_F2xz_S_D2y_Py_bd = I_ERI_F2xz_S_F3y_S_bd+CDY*I_ERI_F2xz_S_D2y_S_bd;
  Double I_ERI_Fx2y_S_D2y_Py_bd = I_ERI_Fx2y_S_F3y_S_bd+CDY*I_ERI_Fx2y_S_D2y_S_bd;
  Double I_ERI_Fxyz_S_D2y_Py_bd = I_ERI_Fxyz_S_F3y_S_bd+CDY*I_ERI_Fxyz_S_D2y_S_bd;
  Double I_ERI_Fx2z_S_D2y_Py_bd = I_ERI_Fx2z_S_F3y_S_bd+CDY*I_ERI_Fx2z_S_D2y_S_bd;
  Double I_ERI_F3y_S_D2y_Py_bd = I_ERI_F3y_S_F3y_S_bd+CDY*I_ERI_F3y_S_D2y_S_bd;
  Double I_ERI_F2yz_S_D2y_Py_bd = I_ERI_F2yz_S_F3y_S_bd+CDY*I_ERI_F2yz_S_D2y_S_bd;
  Double I_ERI_Fy2z_S_D2y_Py_bd = I_ERI_Fy2z_S_F3y_S_bd+CDY*I_ERI_Fy2z_S_D2y_S_bd;
  Double I_ERI_F3z_S_D2y_Py_bd = I_ERI_F3z_S_F3y_S_bd+CDY*I_ERI_F3z_S_D2y_S_bd;
  Double I_ERI_F3x_S_Dyz_Py_bd = I_ERI_F3x_S_F2yz_S_bd+CDY*I_ERI_F3x_S_Dyz_S_bd;
  Double I_ERI_F2xy_S_Dyz_Py_bd = I_ERI_F2xy_S_F2yz_S_bd+CDY*I_ERI_F2xy_S_Dyz_S_bd;
  Double I_ERI_F2xz_S_Dyz_Py_bd = I_ERI_F2xz_S_F2yz_S_bd+CDY*I_ERI_F2xz_S_Dyz_S_bd;
  Double I_ERI_Fx2y_S_Dyz_Py_bd = I_ERI_Fx2y_S_F2yz_S_bd+CDY*I_ERI_Fx2y_S_Dyz_S_bd;
  Double I_ERI_Fxyz_S_Dyz_Py_bd = I_ERI_Fxyz_S_F2yz_S_bd+CDY*I_ERI_Fxyz_S_Dyz_S_bd;
  Double I_ERI_Fx2z_S_Dyz_Py_bd = I_ERI_Fx2z_S_F2yz_S_bd+CDY*I_ERI_Fx2z_S_Dyz_S_bd;
  Double I_ERI_F3y_S_Dyz_Py_bd = I_ERI_F3y_S_F2yz_S_bd+CDY*I_ERI_F3y_S_Dyz_S_bd;
  Double I_ERI_F2yz_S_Dyz_Py_bd = I_ERI_F2yz_S_F2yz_S_bd+CDY*I_ERI_F2yz_S_Dyz_S_bd;
  Double I_ERI_Fy2z_S_Dyz_Py_bd = I_ERI_Fy2z_S_F2yz_S_bd+CDY*I_ERI_Fy2z_S_Dyz_S_bd;
  Double I_ERI_F3z_S_Dyz_Py_bd = I_ERI_F3z_S_F2yz_S_bd+CDY*I_ERI_F3z_S_Dyz_S_bd;
  Double I_ERI_F3x_S_D2z_Py_bd = I_ERI_F3x_S_Fy2z_S_bd+CDY*I_ERI_F3x_S_D2z_S_bd;
  Double I_ERI_F2xy_S_D2z_Py_bd = I_ERI_F2xy_S_Fy2z_S_bd+CDY*I_ERI_F2xy_S_D2z_S_bd;
  Double I_ERI_F2xz_S_D2z_Py_bd = I_ERI_F2xz_S_Fy2z_S_bd+CDY*I_ERI_F2xz_S_D2z_S_bd;
  Double I_ERI_Fx2y_S_D2z_Py_bd = I_ERI_Fx2y_S_Fy2z_S_bd+CDY*I_ERI_Fx2y_S_D2z_S_bd;
  Double I_ERI_Fxyz_S_D2z_Py_bd = I_ERI_Fxyz_S_Fy2z_S_bd+CDY*I_ERI_Fxyz_S_D2z_S_bd;
  Double I_ERI_Fx2z_S_D2z_Py_bd = I_ERI_Fx2z_S_Fy2z_S_bd+CDY*I_ERI_Fx2z_S_D2z_S_bd;
  Double I_ERI_F3y_S_D2z_Py_bd = I_ERI_F3y_S_Fy2z_S_bd+CDY*I_ERI_F3y_S_D2z_S_bd;
  Double I_ERI_F2yz_S_D2z_Py_bd = I_ERI_F2yz_S_Fy2z_S_bd+CDY*I_ERI_F2yz_S_D2z_S_bd;
  Double I_ERI_Fy2z_S_D2z_Py_bd = I_ERI_Fy2z_S_Fy2z_S_bd+CDY*I_ERI_Fy2z_S_D2z_S_bd;
  Double I_ERI_F3z_S_D2z_Py_bd = I_ERI_F3z_S_Fy2z_S_bd+CDY*I_ERI_F3z_S_D2z_S_bd;
  Double I_ERI_F3x_S_D2x_Pz_bd = I_ERI_F3x_S_F2xz_S_bd+CDZ*I_ERI_F3x_S_D2x_S_bd;
  Double I_ERI_F2xy_S_D2x_Pz_bd = I_ERI_F2xy_S_F2xz_S_bd+CDZ*I_ERI_F2xy_S_D2x_S_bd;
  Double I_ERI_F2xz_S_D2x_Pz_bd = I_ERI_F2xz_S_F2xz_S_bd+CDZ*I_ERI_F2xz_S_D2x_S_bd;
  Double I_ERI_Fx2y_S_D2x_Pz_bd = I_ERI_Fx2y_S_F2xz_S_bd+CDZ*I_ERI_Fx2y_S_D2x_S_bd;
  Double I_ERI_Fxyz_S_D2x_Pz_bd = I_ERI_Fxyz_S_F2xz_S_bd+CDZ*I_ERI_Fxyz_S_D2x_S_bd;
  Double I_ERI_Fx2z_S_D2x_Pz_bd = I_ERI_Fx2z_S_F2xz_S_bd+CDZ*I_ERI_Fx2z_S_D2x_S_bd;
  Double I_ERI_F3y_S_D2x_Pz_bd = I_ERI_F3y_S_F2xz_S_bd+CDZ*I_ERI_F3y_S_D2x_S_bd;
  Double I_ERI_F2yz_S_D2x_Pz_bd = I_ERI_F2yz_S_F2xz_S_bd+CDZ*I_ERI_F2yz_S_D2x_S_bd;
  Double I_ERI_Fy2z_S_D2x_Pz_bd = I_ERI_Fy2z_S_F2xz_S_bd+CDZ*I_ERI_Fy2z_S_D2x_S_bd;
  Double I_ERI_F3z_S_D2x_Pz_bd = I_ERI_F3z_S_F2xz_S_bd+CDZ*I_ERI_F3z_S_D2x_S_bd;
  Double I_ERI_F3x_S_Dxy_Pz_bd = I_ERI_F3x_S_Fxyz_S_bd+CDZ*I_ERI_F3x_S_Dxy_S_bd;
  Double I_ERI_F2xy_S_Dxy_Pz_bd = I_ERI_F2xy_S_Fxyz_S_bd+CDZ*I_ERI_F2xy_S_Dxy_S_bd;
  Double I_ERI_F2xz_S_Dxy_Pz_bd = I_ERI_F2xz_S_Fxyz_S_bd+CDZ*I_ERI_F2xz_S_Dxy_S_bd;
  Double I_ERI_Fx2y_S_Dxy_Pz_bd = I_ERI_Fx2y_S_Fxyz_S_bd+CDZ*I_ERI_Fx2y_S_Dxy_S_bd;
  Double I_ERI_Fxyz_S_Dxy_Pz_bd = I_ERI_Fxyz_S_Fxyz_S_bd+CDZ*I_ERI_Fxyz_S_Dxy_S_bd;
  Double I_ERI_Fx2z_S_Dxy_Pz_bd = I_ERI_Fx2z_S_Fxyz_S_bd+CDZ*I_ERI_Fx2z_S_Dxy_S_bd;
  Double I_ERI_F3y_S_Dxy_Pz_bd = I_ERI_F3y_S_Fxyz_S_bd+CDZ*I_ERI_F3y_S_Dxy_S_bd;
  Double I_ERI_F2yz_S_Dxy_Pz_bd = I_ERI_F2yz_S_Fxyz_S_bd+CDZ*I_ERI_F2yz_S_Dxy_S_bd;
  Double I_ERI_Fy2z_S_Dxy_Pz_bd = I_ERI_Fy2z_S_Fxyz_S_bd+CDZ*I_ERI_Fy2z_S_Dxy_S_bd;
  Double I_ERI_F3z_S_Dxy_Pz_bd = I_ERI_F3z_S_Fxyz_S_bd+CDZ*I_ERI_F3z_S_Dxy_S_bd;
  Double I_ERI_F3x_S_Dxz_Pz_bd = I_ERI_F3x_S_Fx2z_S_bd+CDZ*I_ERI_F3x_S_Dxz_S_bd;
  Double I_ERI_F2xy_S_Dxz_Pz_bd = I_ERI_F2xy_S_Fx2z_S_bd+CDZ*I_ERI_F2xy_S_Dxz_S_bd;
  Double I_ERI_F2xz_S_Dxz_Pz_bd = I_ERI_F2xz_S_Fx2z_S_bd+CDZ*I_ERI_F2xz_S_Dxz_S_bd;
  Double I_ERI_Fx2y_S_Dxz_Pz_bd = I_ERI_Fx2y_S_Fx2z_S_bd+CDZ*I_ERI_Fx2y_S_Dxz_S_bd;
  Double I_ERI_Fxyz_S_Dxz_Pz_bd = I_ERI_Fxyz_S_Fx2z_S_bd+CDZ*I_ERI_Fxyz_S_Dxz_S_bd;
  Double I_ERI_Fx2z_S_Dxz_Pz_bd = I_ERI_Fx2z_S_Fx2z_S_bd+CDZ*I_ERI_Fx2z_S_Dxz_S_bd;
  Double I_ERI_F3y_S_Dxz_Pz_bd = I_ERI_F3y_S_Fx2z_S_bd+CDZ*I_ERI_F3y_S_Dxz_S_bd;
  Double I_ERI_F2yz_S_Dxz_Pz_bd = I_ERI_F2yz_S_Fx2z_S_bd+CDZ*I_ERI_F2yz_S_Dxz_S_bd;
  Double I_ERI_Fy2z_S_Dxz_Pz_bd = I_ERI_Fy2z_S_Fx2z_S_bd+CDZ*I_ERI_Fy2z_S_Dxz_S_bd;
  Double I_ERI_F3z_S_Dxz_Pz_bd = I_ERI_F3z_S_Fx2z_S_bd+CDZ*I_ERI_F3z_S_Dxz_S_bd;
  Double I_ERI_F3x_S_D2y_Pz_bd = I_ERI_F3x_S_F2yz_S_bd+CDZ*I_ERI_F3x_S_D2y_S_bd;
  Double I_ERI_F2xy_S_D2y_Pz_bd = I_ERI_F2xy_S_F2yz_S_bd+CDZ*I_ERI_F2xy_S_D2y_S_bd;
  Double I_ERI_F2xz_S_D2y_Pz_bd = I_ERI_F2xz_S_F2yz_S_bd+CDZ*I_ERI_F2xz_S_D2y_S_bd;
  Double I_ERI_Fx2y_S_D2y_Pz_bd = I_ERI_Fx2y_S_F2yz_S_bd+CDZ*I_ERI_Fx2y_S_D2y_S_bd;
  Double I_ERI_Fxyz_S_D2y_Pz_bd = I_ERI_Fxyz_S_F2yz_S_bd+CDZ*I_ERI_Fxyz_S_D2y_S_bd;
  Double I_ERI_Fx2z_S_D2y_Pz_bd = I_ERI_Fx2z_S_F2yz_S_bd+CDZ*I_ERI_Fx2z_S_D2y_S_bd;
  Double I_ERI_F3y_S_D2y_Pz_bd = I_ERI_F3y_S_F2yz_S_bd+CDZ*I_ERI_F3y_S_D2y_S_bd;
  Double I_ERI_F2yz_S_D2y_Pz_bd = I_ERI_F2yz_S_F2yz_S_bd+CDZ*I_ERI_F2yz_S_D2y_S_bd;
  Double I_ERI_Fy2z_S_D2y_Pz_bd = I_ERI_Fy2z_S_F2yz_S_bd+CDZ*I_ERI_Fy2z_S_D2y_S_bd;
  Double I_ERI_F3z_S_D2y_Pz_bd = I_ERI_F3z_S_F2yz_S_bd+CDZ*I_ERI_F3z_S_D2y_S_bd;
  Double I_ERI_F3x_S_Dyz_Pz_bd = I_ERI_F3x_S_Fy2z_S_bd+CDZ*I_ERI_F3x_S_Dyz_S_bd;
  Double I_ERI_F2xy_S_Dyz_Pz_bd = I_ERI_F2xy_S_Fy2z_S_bd+CDZ*I_ERI_F2xy_S_Dyz_S_bd;
  Double I_ERI_F2xz_S_Dyz_Pz_bd = I_ERI_F2xz_S_Fy2z_S_bd+CDZ*I_ERI_F2xz_S_Dyz_S_bd;
  Double I_ERI_Fx2y_S_Dyz_Pz_bd = I_ERI_Fx2y_S_Fy2z_S_bd+CDZ*I_ERI_Fx2y_S_Dyz_S_bd;
  Double I_ERI_Fxyz_S_Dyz_Pz_bd = I_ERI_Fxyz_S_Fy2z_S_bd+CDZ*I_ERI_Fxyz_S_Dyz_S_bd;
  Double I_ERI_Fx2z_S_Dyz_Pz_bd = I_ERI_Fx2z_S_Fy2z_S_bd+CDZ*I_ERI_Fx2z_S_Dyz_S_bd;
  Double I_ERI_F3y_S_Dyz_Pz_bd = I_ERI_F3y_S_Fy2z_S_bd+CDZ*I_ERI_F3y_S_Dyz_S_bd;
  Double I_ERI_F2yz_S_Dyz_Pz_bd = I_ERI_F2yz_S_Fy2z_S_bd+CDZ*I_ERI_F2yz_S_Dyz_S_bd;
  Double I_ERI_Fy2z_S_Dyz_Pz_bd = I_ERI_Fy2z_S_Fy2z_S_bd+CDZ*I_ERI_Fy2z_S_Dyz_S_bd;
  Double I_ERI_F3z_S_Dyz_Pz_bd = I_ERI_F3z_S_Fy2z_S_bd+CDZ*I_ERI_F3z_S_Dyz_S_bd;
  Double I_ERI_F3x_S_D2z_Pz_bd = I_ERI_F3x_S_F3z_S_bd+CDZ*I_ERI_F3x_S_D2z_S_bd;
  Double I_ERI_F2xy_S_D2z_Pz_bd = I_ERI_F2xy_S_F3z_S_bd+CDZ*I_ERI_F2xy_S_D2z_S_bd;
  Double I_ERI_F2xz_S_D2z_Pz_bd = I_ERI_F2xz_S_F3z_S_bd+CDZ*I_ERI_F2xz_S_D2z_S_bd;
  Double I_ERI_Fx2y_S_D2z_Pz_bd = I_ERI_Fx2y_S_F3z_S_bd+CDZ*I_ERI_Fx2y_S_D2z_S_bd;
  Double I_ERI_Fxyz_S_D2z_Pz_bd = I_ERI_Fxyz_S_F3z_S_bd+CDZ*I_ERI_Fxyz_S_D2z_S_bd;
  Double I_ERI_Fx2z_S_D2z_Pz_bd = I_ERI_Fx2z_S_F3z_S_bd+CDZ*I_ERI_Fx2z_S_D2z_S_bd;
  Double I_ERI_F3y_S_D2z_Pz_bd = I_ERI_F3y_S_F3z_S_bd+CDZ*I_ERI_F3y_S_D2z_S_bd;
  Double I_ERI_F2yz_S_D2z_Pz_bd = I_ERI_F2yz_S_F3z_S_bd+CDZ*I_ERI_F2yz_S_D2z_S_bd;
  Double I_ERI_Fy2z_S_D2z_Pz_bd = I_ERI_Fy2z_S_F3z_S_bd+CDZ*I_ERI_Fy2z_S_D2z_S_bd;
  Double I_ERI_F3z_S_D2z_Pz_bd = I_ERI_F3z_S_F3z_S_bd+CDZ*I_ERI_F3z_S_D2z_S_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_D_P_bd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_bd
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bd
   ************************************************************/
  Double I_ERI_G4x_S_D2x_Px_bd = I_ERI_G4x_S_F3x_S_bd+CDX*I_ERI_G4x_S_D2x_S_bd;
  Double I_ERI_G3xy_S_D2x_Px_bd = I_ERI_G3xy_S_F3x_S_bd+CDX*I_ERI_G3xy_S_D2x_S_bd;
  Double I_ERI_G3xz_S_D2x_Px_bd = I_ERI_G3xz_S_F3x_S_bd+CDX*I_ERI_G3xz_S_D2x_S_bd;
  Double I_ERI_G2x2y_S_D2x_Px_bd = I_ERI_G2x2y_S_F3x_S_bd+CDX*I_ERI_G2x2y_S_D2x_S_bd;
  Double I_ERI_G2xyz_S_D2x_Px_bd = I_ERI_G2xyz_S_F3x_S_bd+CDX*I_ERI_G2xyz_S_D2x_S_bd;
  Double I_ERI_G2x2z_S_D2x_Px_bd = I_ERI_G2x2z_S_F3x_S_bd+CDX*I_ERI_G2x2z_S_D2x_S_bd;
  Double I_ERI_Gx3y_S_D2x_Px_bd = I_ERI_Gx3y_S_F3x_S_bd+CDX*I_ERI_Gx3y_S_D2x_S_bd;
  Double I_ERI_Gx2yz_S_D2x_Px_bd = I_ERI_Gx2yz_S_F3x_S_bd+CDX*I_ERI_Gx2yz_S_D2x_S_bd;
  Double I_ERI_Gxy2z_S_D2x_Px_bd = I_ERI_Gxy2z_S_F3x_S_bd+CDX*I_ERI_Gxy2z_S_D2x_S_bd;
  Double I_ERI_Gx3z_S_D2x_Px_bd = I_ERI_Gx3z_S_F3x_S_bd+CDX*I_ERI_Gx3z_S_D2x_S_bd;
  Double I_ERI_G4y_S_D2x_Px_bd = I_ERI_G4y_S_F3x_S_bd+CDX*I_ERI_G4y_S_D2x_S_bd;
  Double I_ERI_G3yz_S_D2x_Px_bd = I_ERI_G3yz_S_F3x_S_bd+CDX*I_ERI_G3yz_S_D2x_S_bd;
  Double I_ERI_G2y2z_S_D2x_Px_bd = I_ERI_G2y2z_S_F3x_S_bd+CDX*I_ERI_G2y2z_S_D2x_S_bd;
  Double I_ERI_Gy3z_S_D2x_Px_bd = I_ERI_Gy3z_S_F3x_S_bd+CDX*I_ERI_Gy3z_S_D2x_S_bd;
  Double I_ERI_G4z_S_D2x_Px_bd = I_ERI_G4z_S_F3x_S_bd+CDX*I_ERI_G4z_S_D2x_S_bd;
  Double I_ERI_G4x_S_Dxy_Px_bd = I_ERI_G4x_S_F2xy_S_bd+CDX*I_ERI_G4x_S_Dxy_S_bd;
  Double I_ERI_G3xy_S_Dxy_Px_bd = I_ERI_G3xy_S_F2xy_S_bd+CDX*I_ERI_G3xy_S_Dxy_S_bd;
  Double I_ERI_G3xz_S_Dxy_Px_bd = I_ERI_G3xz_S_F2xy_S_bd+CDX*I_ERI_G3xz_S_Dxy_S_bd;
  Double I_ERI_G2x2y_S_Dxy_Px_bd = I_ERI_G2x2y_S_F2xy_S_bd+CDX*I_ERI_G2x2y_S_Dxy_S_bd;
  Double I_ERI_G2xyz_S_Dxy_Px_bd = I_ERI_G2xyz_S_F2xy_S_bd+CDX*I_ERI_G2xyz_S_Dxy_S_bd;
  Double I_ERI_G2x2z_S_Dxy_Px_bd = I_ERI_G2x2z_S_F2xy_S_bd+CDX*I_ERI_G2x2z_S_Dxy_S_bd;
  Double I_ERI_Gx3y_S_Dxy_Px_bd = I_ERI_Gx3y_S_F2xy_S_bd+CDX*I_ERI_Gx3y_S_Dxy_S_bd;
  Double I_ERI_Gx2yz_S_Dxy_Px_bd = I_ERI_Gx2yz_S_F2xy_S_bd+CDX*I_ERI_Gx2yz_S_Dxy_S_bd;
  Double I_ERI_Gxy2z_S_Dxy_Px_bd = I_ERI_Gxy2z_S_F2xy_S_bd+CDX*I_ERI_Gxy2z_S_Dxy_S_bd;
  Double I_ERI_Gx3z_S_Dxy_Px_bd = I_ERI_Gx3z_S_F2xy_S_bd+CDX*I_ERI_Gx3z_S_Dxy_S_bd;
  Double I_ERI_G4y_S_Dxy_Px_bd = I_ERI_G4y_S_F2xy_S_bd+CDX*I_ERI_G4y_S_Dxy_S_bd;
  Double I_ERI_G3yz_S_Dxy_Px_bd = I_ERI_G3yz_S_F2xy_S_bd+CDX*I_ERI_G3yz_S_Dxy_S_bd;
  Double I_ERI_G2y2z_S_Dxy_Px_bd = I_ERI_G2y2z_S_F2xy_S_bd+CDX*I_ERI_G2y2z_S_Dxy_S_bd;
  Double I_ERI_Gy3z_S_Dxy_Px_bd = I_ERI_Gy3z_S_F2xy_S_bd+CDX*I_ERI_Gy3z_S_Dxy_S_bd;
  Double I_ERI_G4z_S_Dxy_Px_bd = I_ERI_G4z_S_F2xy_S_bd+CDX*I_ERI_G4z_S_Dxy_S_bd;
  Double I_ERI_G4x_S_Dxz_Px_bd = I_ERI_G4x_S_F2xz_S_bd+CDX*I_ERI_G4x_S_Dxz_S_bd;
  Double I_ERI_G3xy_S_Dxz_Px_bd = I_ERI_G3xy_S_F2xz_S_bd+CDX*I_ERI_G3xy_S_Dxz_S_bd;
  Double I_ERI_G3xz_S_Dxz_Px_bd = I_ERI_G3xz_S_F2xz_S_bd+CDX*I_ERI_G3xz_S_Dxz_S_bd;
  Double I_ERI_G2x2y_S_Dxz_Px_bd = I_ERI_G2x2y_S_F2xz_S_bd+CDX*I_ERI_G2x2y_S_Dxz_S_bd;
  Double I_ERI_G2xyz_S_Dxz_Px_bd = I_ERI_G2xyz_S_F2xz_S_bd+CDX*I_ERI_G2xyz_S_Dxz_S_bd;
  Double I_ERI_G2x2z_S_Dxz_Px_bd = I_ERI_G2x2z_S_F2xz_S_bd+CDX*I_ERI_G2x2z_S_Dxz_S_bd;
  Double I_ERI_Gx3y_S_Dxz_Px_bd = I_ERI_Gx3y_S_F2xz_S_bd+CDX*I_ERI_Gx3y_S_Dxz_S_bd;
  Double I_ERI_Gx2yz_S_Dxz_Px_bd = I_ERI_Gx2yz_S_F2xz_S_bd+CDX*I_ERI_Gx2yz_S_Dxz_S_bd;
  Double I_ERI_Gxy2z_S_Dxz_Px_bd = I_ERI_Gxy2z_S_F2xz_S_bd+CDX*I_ERI_Gxy2z_S_Dxz_S_bd;
  Double I_ERI_Gx3z_S_Dxz_Px_bd = I_ERI_Gx3z_S_F2xz_S_bd+CDX*I_ERI_Gx3z_S_Dxz_S_bd;
  Double I_ERI_G4y_S_Dxz_Px_bd = I_ERI_G4y_S_F2xz_S_bd+CDX*I_ERI_G4y_S_Dxz_S_bd;
  Double I_ERI_G3yz_S_Dxz_Px_bd = I_ERI_G3yz_S_F2xz_S_bd+CDX*I_ERI_G3yz_S_Dxz_S_bd;
  Double I_ERI_G2y2z_S_Dxz_Px_bd = I_ERI_G2y2z_S_F2xz_S_bd+CDX*I_ERI_G2y2z_S_Dxz_S_bd;
  Double I_ERI_Gy3z_S_Dxz_Px_bd = I_ERI_Gy3z_S_F2xz_S_bd+CDX*I_ERI_Gy3z_S_Dxz_S_bd;
  Double I_ERI_G4z_S_Dxz_Px_bd = I_ERI_G4z_S_F2xz_S_bd+CDX*I_ERI_G4z_S_Dxz_S_bd;
  Double I_ERI_G4x_S_D2y_Px_bd = I_ERI_G4x_S_Fx2y_S_bd+CDX*I_ERI_G4x_S_D2y_S_bd;
  Double I_ERI_G3xy_S_D2y_Px_bd = I_ERI_G3xy_S_Fx2y_S_bd+CDX*I_ERI_G3xy_S_D2y_S_bd;
  Double I_ERI_G3xz_S_D2y_Px_bd = I_ERI_G3xz_S_Fx2y_S_bd+CDX*I_ERI_G3xz_S_D2y_S_bd;
  Double I_ERI_G2x2y_S_D2y_Px_bd = I_ERI_G2x2y_S_Fx2y_S_bd+CDX*I_ERI_G2x2y_S_D2y_S_bd;
  Double I_ERI_G2xyz_S_D2y_Px_bd = I_ERI_G2xyz_S_Fx2y_S_bd+CDX*I_ERI_G2xyz_S_D2y_S_bd;
  Double I_ERI_G2x2z_S_D2y_Px_bd = I_ERI_G2x2z_S_Fx2y_S_bd+CDX*I_ERI_G2x2z_S_D2y_S_bd;
  Double I_ERI_Gx3y_S_D2y_Px_bd = I_ERI_Gx3y_S_Fx2y_S_bd+CDX*I_ERI_Gx3y_S_D2y_S_bd;
  Double I_ERI_Gx2yz_S_D2y_Px_bd = I_ERI_Gx2yz_S_Fx2y_S_bd+CDX*I_ERI_Gx2yz_S_D2y_S_bd;
  Double I_ERI_Gxy2z_S_D2y_Px_bd = I_ERI_Gxy2z_S_Fx2y_S_bd+CDX*I_ERI_Gxy2z_S_D2y_S_bd;
  Double I_ERI_Gx3z_S_D2y_Px_bd = I_ERI_Gx3z_S_Fx2y_S_bd+CDX*I_ERI_Gx3z_S_D2y_S_bd;
  Double I_ERI_G4y_S_D2y_Px_bd = I_ERI_G4y_S_Fx2y_S_bd+CDX*I_ERI_G4y_S_D2y_S_bd;
  Double I_ERI_G3yz_S_D2y_Px_bd = I_ERI_G3yz_S_Fx2y_S_bd+CDX*I_ERI_G3yz_S_D2y_S_bd;
  Double I_ERI_G2y2z_S_D2y_Px_bd = I_ERI_G2y2z_S_Fx2y_S_bd+CDX*I_ERI_G2y2z_S_D2y_S_bd;
  Double I_ERI_Gy3z_S_D2y_Px_bd = I_ERI_Gy3z_S_Fx2y_S_bd+CDX*I_ERI_Gy3z_S_D2y_S_bd;
  Double I_ERI_G4z_S_D2y_Px_bd = I_ERI_G4z_S_Fx2y_S_bd+CDX*I_ERI_G4z_S_D2y_S_bd;
  Double I_ERI_G4x_S_Dyz_Px_bd = I_ERI_G4x_S_Fxyz_S_bd+CDX*I_ERI_G4x_S_Dyz_S_bd;
  Double I_ERI_G3xy_S_Dyz_Px_bd = I_ERI_G3xy_S_Fxyz_S_bd+CDX*I_ERI_G3xy_S_Dyz_S_bd;
  Double I_ERI_G3xz_S_Dyz_Px_bd = I_ERI_G3xz_S_Fxyz_S_bd+CDX*I_ERI_G3xz_S_Dyz_S_bd;
  Double I_ERI_G2x2y_S_Dyz_Px_bd = I_ERI_G2x2y_S_Fxyz_S_bd+CDX*I_ERI_G2x2y_S_Dyz_S_bd;
  Double I_ERI_G2xyz_S_Dyz_Px_bd = I_ERI_G2xyz_S_Fxyz_S_bd+CDX*I_ERI_G2xyz_S_Dyz_S_bd;
  Double I_ERI_G2x2z_S_Dyz_Px_bd = I_ERI_G2x2z_S_Fxyz_S_bd+CDX*I_ERI_G2x2z_S_Dyz_S_bd;
  Double I_ERI_Gx3y_S_Dyz_Px_bd = I_ERI_Gx3y_S_Fxyz_S_bd+CDX*I_ERI_Gx3y_S_Dyz_S_bd;
  Double I_ERI_Gx2yz_S_Dyz_Px_bd = I_ERI_Gx2yz_S_Fxyz_S_bd+CDX*I_ERI_Gx2yz_S_Dyz_S_bd;
  Double I_ERI_Gxy2z_S_Dyz_Px_bd = I_ERI_Gxy2z_S_Fxyz_S_bd+CDX*I_ERI_Gxy2z_S_Dyz_S_bd;
  Double I_ERI_Gx3z_S_Dyz_Px_bd = I_ERI_Gx3z_S_Fxyz_S_bd+CDX*I_ERI_Gx3z_S_Dyz_S_bd;
  Double I_ERI_G4y_S_Dyz_Px_bd = I_ERI_G4y_S_Fxyz_S_bd+CDX*I_ERI_G4y_S_Dyz_S_bd;
  Double I_ERI_G3yz_S_Dyz_Px_bd = I_ERI_G3yz_S_Fxyz_S_bd+CDX*I_ERI_G3yz_S_Dyz_S_bd;
  Double I_ERI_G2y2z_S_Dyz_Px_bd = I_ERI_G2y2z_S_Fxyz_S_bd+CDX*I_ERI_G2y2z_S_Dyz_S_bd;
  Double I_ERI_Gy3z_S_Dyz_Px_bd = I_ERI_Gy3z_S_Fxyz_S_bd+CDX*I_ERI_Gy3z_S_Dyz_S_bd;
  Double I_ERI_G4z_S_Dyz_Px_bd = I_ERI_G4z_S_Fxyz_S_bd+CDX*I_ERI_G4z_S_Dyz_S_bd;
  Double I_ERI_G4x_S_D2z_Px_bd = I_ERI_G4x_S_Fx2z_S_bd+CDX*I_ERI_G4x_S_D2z_S_bd;
  Double I_ERI_G3xy_S_D2z_Px_bd = I_ERI_G3xy_S_Fx2z_S_bd+CDX*I_ERI_G3xy_S_D2z_S_bd;
  Double I_ERI_G3xz_S_D2z_Px_bd = I_ERI_G3xz_S_Fx2z_S_bd+CDX*I_ERI_G3xz_S_D2z_S_bd;
  Double I_ERI_G2x2y_S_D2z_Px_bd = I_ERI_G2x2y_S_Fx2z_S_bd+CDX*I_ERI_G2x2y_S_D2z_S_bd;
  Double I_ERI_G2xyz_S_D2z_Px_bd = I_ERI_G2xyz_S_Fx2z_S_bd+CDX*I_ERI_G2xyz_S_D2z_S_bd;
  Double I_ERI_G2x2z_S_D2z_Px_bd = I_ERI_G2x2z_S_Fx2z_S_bd+CDX*I_ERI_G2x2z_S_D2z_S_bd;
  Double I_ERI_Gx3y_S_D2z_Px_bd = I_ERI_Gx3y_S_Fx2z_S_bd+CDX*I_ERI_Gx3y_S_D2z_S_bd;
  Double I_ERI_Gx2yz_S_D2z_Px_bd = I_ERI_Gx2yz_S_Fx2z_S_bd+CDX*I_ERI_Gx2yz_S_D2z_S_bd;
  Double I_ERI_Gxy2z_S_D2z_Px_bd = I_ERI_Gxy2z_S_Fx2z_S_bd+CDX*I_ERI_Gxy2z_S_D2z_S_bd;
  Double I_ERI_Gx3z_S_D2z_Px_bd = I_ERI_Gx3z_S_Fx2z_S_bd+CDX*I_ERI_Gx3z_S_D2z_S_bd;
  Double I_ERI_G4y_S_D2z_Px_bd = I_ERI_G4y_S_Fx2z_S_bd+CDX*I_ERI_G4y_S_D2z_S_bd;
  Double I_ERI_G3yz_S_D2z_Px_bd = I_ERI_G3yz_S_Fx2z_S_bd+CDX*I_ERI_G3yz_S_D2z_S_bd;
  Double I_ERI_G2y2z_S_D2z_Px_bd = I_ERI_G2y2z_S_Fx2z_S_bd+CDX*I_ERI_G2y2z_S_D2z_S_bd;
  Double I_ERI_Gy3z_S_D2z_Px_bd = I_ERI_Gy3z_S_Fx2z_S_bd+CDX*I_ERI_Gy3z_S_D2z_S_bd;
  Double I_ERI_G4z_S_D2z_Px_bd = I_ERI_G4z_S_Fx2z_S_bd+CDX*I_ERI_G4z_S_D2z_S_bd;
  Double I_ERI_G4x_S_D2x_Py_bd = I_ERI_G4x_S_F2xy_S_bd+CDY*I_ERI_G4x_S_D2x_S_bd;
  Double I_ERI_G3xy_S_D2x_Py_bd = I_ERI_G3xy_S_F2xy_S_bd+CDY*I_ERI_G3xy_S_D2x_S_bd;
  Double I_ERI_G3xz_S_D2x_Py_bd = I_ERI_G3xz_S_F2xy_S_bd+CDY*I_ERI_G3xz_S_D2x_S_bd;
  Double I_ERI_G2x2y_S_D2x_Py_bd = I_ERI_G2x2y_S_F2xy_S_bd+CDY*I_ERI_G2x2y_S_D2x_S_bd;
  Double I_ERI_G2xyz_S_D2x_Py_bd = I_ERI_G2xyz_S_F2xy_S_bd+CDY*I_ERI_G2xyz_S_D2x_S_bd;
  Double I_ERI_G2x2z_S_D2x_Py_bd = I_ERI_G2x2z_S_F2xy_S_bd+CDY*I_ERI_G2x2z_S_D2x_S_bd;
  Double I_ERI_Gx3y_S_D2x_Py_bd = I_ERI_Gx3y_S_F2xy_S_bd+CDY*I_ERI_Gx3y_S_D2x_S_bd;
  Double I_ERI_Gx2yz_S_D2x_Py_bd = I_ERI_Gx2yz_S_F2xy_S_bd+CDY*I_ERI_Gx2yz_S_D2x_S_bd;
  Double I_ERI_Gxy2z_S_D2x_Py_bd = I_ERI_Gxy2z_S_F2xy_S_bd+CDY*I_ERI_Gxy2z_S_D2x_S_bd;
  Double I_ERI_Gx3z_S_D2x_Py_bd = I_ERI_Gx3z_S_F2xy_S_bd+CDY*I_ERI_Gx3z_S_D2x_S_bd;
  Double I_ERI_G4y_S_D2x_Py_bd = I_ERI_G4y_S_F2xy_S_bd+CDY*I_ERI_G4y_S_D2x_S_bd;
  Double I_ERI_G3yz_S_D2x_Py_bd = I_ERI_G3yz_S_F2xy_S_bd+CDY*I_ERI_G3yz_S_D2x_S_bd;
  Double I_ERI_G2y2z_S_D2x_Py_bd = I_ERI_G2y2z_S_F2xy_S_bd+CDY*I_ERI_G2y2z_S_D2x_S_bd;
  Double I_ERI_Gy3z_S_D2x_Py_bd = I_ERI_Gy3z_S_F2xy_S_bd+CDY*I_ERI_Gy3z_S_D2x_S_bd;
  Double I_ERI_G4z_S_D2x_Py_bd = I_ERI_G4z_S_F2xy_S_bd+CDY*I_ERI_G4z_S_D2x_S_bd;
  Double I_ERI_G4x_S_Dxy_Py_bd = I_ERI_G4x_S_Fx2y_S_bd+CDY*I_ERI_G4x_S_Dxy_S_bd;
  Double I_ERI_G3xy_S_Dxy_Py_bd = I_ERI_G3xy_S_Fx2y_S_bd+CDY*I_ERI_G3xy_S_Dxy_S_bd;
  Double I_ERI_G3xz_S_Dxy_Py_bd = I_ERI_G3xz_S_Fx2y_S_bd+CDY*I_ERI_G3xz_S_Dxy_S_bd;
  Double I_ERI_G2x2y_S_Dxy_Py_bd = I_ERI_G2x2y_S_Fx2y_S_bd+CDY*I_ERI_G2x2y_S_Dxy_S_bd;
  Double I_ERI_G2xyz_S_Dxy_Py_bd = I_ERI_G2xyz_S_Fx2y_S_bd+CDY*I_ERI_G2xyz_S_Dxy_S_bd;
  Double I_ERI_G2x2z_S_Dxy_Py_bd = I_ERI_G2x2z_S_Fx2y_S_bd+CDY*I_ERI_G2x2z_S_Dxy_S_bd;
  Double I_ERI_Gx3y_S_Dxy_Py_bd = I_ERI_Gx3y_S_Fx2y_S_bd+CDY*I_ERI_Gx3y_S_Dxy_S_bd;
  Double I_ERI_Gx2yz_S_Dxy_Py_bd = I_ERI_Gx2yz_S_Fx2y_S_bd+CDY*I_ERI_Gx2yz_S_Dxy_S_bd;
  Double I_ERI_Gxy2z_S_Dxy_Py_bd = I_ERI_Gxy2z_S_Fx2y_S_bd+CDY*I_ERI_Gxy2z_S_Dxy_S_bd;
  Double I_ERI_Gx3z_S_Dxy_Py_bd = I_ERI_Gx3z_S_Fx2y_S_bd+CDY*I_ERI_Gx3z_S_Dxy_S_bd;
  Double I_ERI_G4y_S_Dxy_Py_bd = I_ERI_G4y_S_Fx2y_S_bd+CDY*I_ERI_G4y_S_Dxy_S_bd;
  Double I_ERI_G3yz_S_Dxy_Py_bd = I_ERI_G3yz_S_Fx2y_S_bd+CDY*I_ERI_G3yz_S_Dxy_S_bd;
  Double I_ERI_G2y2z_S_Dxy_Py_bd = I_ERI_G2y2z_S_Fx2y_S_bd+CDY*I_ERI_G2y2z_S_Dxy_S_bd;
  Double I_ERI_Gy3z_S_Dxy_Py_bd = I_ERI_Gy3z_S_Fx2y_S_bd+CDY*I_ERI_Gy3z_S_Dxy_S_bd;
  Double I_ERI_G4z_S_Dxy_Py_bd = I_ERI_G4z_S_Fx2y_S_bd+CDY*I_ERI_G4z_S_Dxy_S_bd;
  Double I_ERI_G4x_S_Dxz_Py_bd = I_ERI_G4x_S_Fxyz_S_bd+CDY*I_ERI_G4x_S_Dxz_S_bd;
  Double I_ERI_G3xy_S_Dxz_Py_bd = I_ERI_G3xy_S_Fxyz_S_bd+CDY*I_ERI_G3xy_S_Dxz_S_bd;
  Double I_ERI_G3xz_S_Dxz_Py_bd = I_ERI_G3xz_S_Fxyz_S_bd+CDY*I_ERI_G3xz_S_Dxz_S_bd;
  Double I_ERI_G2x2y_S_Dxz_Py_bd = I_ERI_G2x2y_S_Fxyz_S_bd+CDY*I_ERI_G2x2y_S_Dxz_S_bd;
  Double I_ERI_G2xyz_S_Dxz_Py_bd = I_ERI_G2xyz_S_Fxyz_S_bd+CDY*I_ERI_G2xyz_S_Dxz_S_bd;
  Double I_ERI_G2x2z_S_Dxz_Py_bd = I_ERI_G2x2z_S_Fxyz_S_bd+CDY*I_ERI_G2x2z_S_Dxz_S_bd;
  Double I_ERI_Gx3y_S_Dxz_Py_bd = I_ERI_Gx3y_S_Fxyz_S_bd+CDY*I_ERI_Gx3y_S_Dxz_S_bd;
  Double I_ERI_Gx2yz_S_Dxz_Py_bd = I_ERI_Gx2yz_S_Fxyz_S_bd+CDY*I_ERI_Gx2yz_S_Dxz_S_bd;
  Double I_ERI_Gxy2z_S_Dxz_Py_bd = I_ERI_Gxy2z_S_Fxyz_S_bd+CDY*I_ERI_Gxy2z_S_Dxz_S_bd;
  Double I_ERI_Gx3z_S_Dxz_Py_bd = I_ERI_Gx3z_S_Fxyz_S_bd+CDY*I_ERI_Gx3z_S_Dxz_S_bd;
  Double I_ERI_G4y_S_Dxz_Py_bd = I_ERI_G4y_S_Fxyz_S_bd+CDY*I_ERI_G4y_S_Dxz_S_bd;
  Double I_ERI_G3yz_S_Dxz_Py_bd = I_ERI_G3yz_S_Fxyz_S_bd+CDY*I_ERI_G3yz_S_Dxz_S_bd;
  Double I_ERI_G2y2z_S_Dxz_Py_bd = I_ERI_G2y2z_S_Fxyz_S_bd+CDY*I_ERI_G2y2z_S_Dxz_S_bd;
  Double I_ERI_Gy3z_S_Dxz_Py_bd = I_ERI_Gy3z_S_Fxyz_S_bd+CDY*I_ERI_Gy3z_S_Dxz_S_bd;
  Double I_ERI_G4z_S_Dxz_Py_bd = I_ERI_G4z_S_Fxyz_S_bd+CDY*I_ERI_G4z_S_Dxz_S_bd;
  Double I_ERI_G4x_S_D2y_Py_bd = I_ERI_G4x_S_F3y_S_bd+CDY*I_ERI_G4x_S_D2y_S_bd;
  Double I_ERI_G3xy_S_D2y_Py_bd = I_ERI_G3xy_S_F3y_S_bd+CDY*I_ERI_G3xy_S_D2y_S_bd;
  Double I_ERI_G3xz_S_D2y_Py_bd = I_ERI_G3xz_S_F3y_S_bd+CDY*I_ERI_G3xz_S_D2y_S_bd;
  Double I_ERI_G2x2y_S_D2y_Py_bd = I_ERI_G2x2y_S_F3y_S_bd+CDY*I_ERI_G2x2y_S_D2y_S_bd;
  Double I_ERI_G2xyz_S_D2y_Py_bd = I_ERI_G2xyz_S_F3y_S_bd+CDY*I_ERI_G2xyz_S_D2y_S_bd;
  Double I_ERI_G2x2z_S_D2y_Py_bd = I_ERI_G2x2z_S_F3y_S_bd+CDY*I_ERI_G2x2z_S_D2y_S_bd;
  Double I_ERI_Gx3y_S_D2y_Py_bd = I_ERI_Gx3y_S_F3y_S_bd+CDY*I_ERI_Gx3y_S_D2y_S_bd;
  Double I_ERI_Gx2yz_S_D2y_Py_bd = I_ERI_Gx2yz_S_F3y_S_bd+CDY*I_ERI_Gx2yz_S_D2y_S_bd;
  Double I_ERI_Gxy2z_S_D2y_Py_bd = I_ERI_Gxy2z_S_F3y_S_bd+CDY*I_ERI_Gxy2z_S_D2y_S_bd;
  Double I_ERI_Gx3z_S_D2y_Py_bd = I_ERI_Gx3z_S_F3y_S_bd+CDY*I_ERI_Gx3z_S_D2y_S_bd;
  Double I_ERI_G4y_S_D2y_Py_bd = I_ERI_G4y_S_F3y_S_bd+CDY*I_ERI_G4y_S_D2y_S_bd;
  Double I_ERI_G3yz_S_D2y_Py_bd = I_ERI_G3yz_S_F3y_S_bd+CDY*I_ERI_G3yz_S_D2y_S_bd;
  Double I_ERI_G2y2z_S_D2y_Py_bd = I_ERI_G2y2z_S_F3y_S_bd+CDY*I_ERI_G2y2z_S_D2y_S_bd;
  Double I_ERI_Gy3z_S_D2y_Py_bd = I_ERI_Gy3z_S_F3y_S_bd+CDY*I_ERI_Gy3z_S_D2y_S_bd;
  Double I_ERI_G4z_S_D2y_Py_bd = I_ERI_G4z_S_F3y_S_bd+CDY*I_ERI_G4z_S_D2y_S_bd;
  Double I_ERI_G4x_S_Dyz_Py_bd = I_ERI_G4x_S_F2yz_S_bd+CDY*I_ERI_G4x_S_Dyz_S_bd;
  Double I_ERI_G3xy_S_Dyz_Py_bd = I_ERI_G3xy_S_F2yz_S_bd+CDY*I_ERI_G3xy_S_Dyz_S_bd;
  Double I_ERI_G3xz_S_Dyz_Py_bd = I_ERI_G3xz_S_F2yz_S_bd+CDY*I_ERI_G3xz_S_Dyz_S_bd;
  Double I_ERI_G2x2y_S_Dyz_Py_bd = I_ERI_G2x2y_S_F2yz_S_bd+CDY*I_ERI_G2x2y_S_Dyz_S_bd;
  Double I_ERI_G2xyz_S_Dyz_Py_bd = I_ERI_G2xyz_S_F2yz_S_bd+CDY*I_ERI_G2xyz_S_Dyz_S_bd;
  Double I_ERI_G2x2z_S_Dyz_Py_bd = I_ERI_G2x2z_S_F2yz_S_bd+CDY*I_ERI_G2x2z_S_Dyz_S_bd;
  Double I_ERI_Gx3y_S_Dyz_Py_bd = I_ERI_Gx3y_S_F2yz_S_bd+CDY*I_ERI_Gx3y_S_Dyz_S_bd;
  Double I_ERI_Gx2yz_S_Dyz_Py_bd = I_ERI_Gx2yz_S_F2yz_S_bd+CDY*I_ERI_Gx2yz_S_Dyz_S_bd;
  Double I_ERI_Gxy2z_S_Dyz_Py_bd = I_ERI_Gxy2z_S_F2yz_S_bd+CDY*I_ERI_Gxy2z_S_Dyz_S_bd;
  Double I_ERI_Gx3z_S_Dyz_Py_bd = I_ERI_Gx3z_S_F2yz_S_bd+CDY*I_ERI_Gx3z_S_Dyz_S_bd;
  Double I_ERI_G4y_S_Dyz_Py_bd = I_ERI_G4y_S_F2yz_S_bd+CDY*I_ERI_G4y_S_Dyz_S_bd;
  Double I_ERI_G3yz_S_Dyz_Py_bd = I_ERI_G3yz_S_F2yz_S_bd+CDY*I_ERI_G3yz_S_Dyz_S_bd;
  Double I_ERI_G2y2z_S_Dyz_Py_bd = I_ERI_G2y2z_S_F2yz_S_bd+CDY*I_ERI_G2y2z_S_Dyz_S_bd;
  Double I_ERI_Gy3z_S_Dyz_Py_bd = I_ERI_Gy3z_S_F2yz_S_bd+CDY*I_ERI_Gy3z_S_Dyz_S_bd;
  Double I_ERI_G4z_S_Dyz_Py_bd = I_ERI_G4z_S_F2yz_S_bd+CDY*I_ERI_G4z_S_Dyz_S_bd;
  Double I_ERI_G4x_S_D2z_Py_bd = I_ERI_G4x_S_Fy2z_S_bd+CDY*I_ERI_G4x_S_D2z_S_bd;
  Double I_ERI_G3xy_S_D2z_Py_bd = I_ERI_G3xy_S_Fy2z_S_bd+CDY*I_ERI_G3xy_S_D2z_S_bd;
  Double I_ERI_G3xz_S_D2z_Py_bd = I_ERI_G3xz_S_Fy2z_S_bd+CDY*I_ERI_G3xz_S_D2z_S_bd;
  Double I_ERI_G2x2y_S_D2z_Py_bd = I_ERI_G2x2y_S_Fy2z_S_bd+CDY*I_ERI_G2x2y_S_D2z_S_bd;
  Double I_ERI_G2xyz_S_D2z_Py_bd = I_ERI_G2xyz_S_Fy2z_S_bd+CDY*I_ERI_G2xyz_S_D2z_S_bd;
  Double I_ERI_G2x2z_S_D2z_Py_bd = I_ERI_G2x2z_S_Fy2z_S_bd+CDY*I_ERI_G2x2z_S_D2z_S_bd;
  Double I_ERI_Gx3y_S_D2z_Py_bd = I_ERI_Gx3y_S_Fy2z_S_bd+CDY*I_ERI_Gx3y_S_D2z_S_bd;
  Double I_ERI_Gx2yz_S_D2z_Py_bd = I_ERI_Gx2yz_S_Fy2z_S_bd+CDY*I_ERI_Gx2yz_S_D2z_S_bd;
  Double I_ERI_Gxy2z_S_D2z_Py_bd = I_ERI_Gxy2z_S_Fy2z_S_bd+CDY*I_ERI_Gxy2z_S_D2z_S_bd;
  Double I_ERI_Gx3z_S_D2z_Py_bd = I_ERI_Gx3z_S_Fy2z_S_bd+CDY*I_ERI_Gx3z_S_D2z_S_bd;
  Double I_ERI_G4y_S_D2z_Py_bd = I_ERI_G4y_S_Fy2z_S_bd+CDY*I_ERI_G4y_S_D2z_S_bd;
  Double I_ERI_G3yz_S_D2z_Py_bd = I_ERI_G3yz_S_Fy2z_S_bd+CDY*I_ERI_G3yz_S_D2z_S_bd;
  Double I_ERI_G2y2z_S_D2z_Py_bd = I_ERI_G2y2z_S_Fy2z_S_bd+CDY*I_ERI_G2y2z_S_D2z_S_bd;
  Double I_ERI_Gy3z_S_D2z_Py_bd = I_ERI_Gy3z_S_Fy2z_S_bd+CDY*I_ERI_Gy3z_S_D2z_S_bd;
  Double I_ERI_G4z_S_D2z_Py_bd = I_ERI_G4z_S_Fy2z_S_bd+CDY*I_ERI_G4z_S_D2z_S_bd;
  Double I_ERI_G4x_S_D2x_Pz_bd = I_ERI_G4x_S_F2xz_S_bd+CDZ*I_ERI_G4x_S_D2x_S_bd;
  Double I_ERI_G3xy_S_D2x_Pz_bd = I_ERI_G3xy_S_F2xz_S_bd+CDZ*I_ERI_G3xy_S_D2x_S_bd;
  Double I_ERI_G3xz_S_D2x_Pz_bd = I_ERI_G3xz_S_F2xz_S_bd+CDZ*I_ERI_G3xz_S_D2x_S_bd;
  Double I_ERI_G2x2y_S_D2x_Pz_bd = I_ERI_G2x2y_S_F2xz_S_bd+CDZ*I_ERI_G2x2y_S_D2x_S_bd;
  Double I_ERI_G2xyz_S_D2x_Pz_bd = I_ERI_G2xyz_S_F2xz_S_bd+CDZ*I_ERI_G2xyz_S_D2x_S_bd;
  Double I_ERI_G2x2z_S_D2x_Pz_bd = I_ERI_G2x2z_S_F2xz_S_bd+CDZ*I_ERI_G2x2z_S_D2x_S_bd;
  Double I_ERI_Gx3y_S_D2x_Pz_bd = I_ERI_Gx3y_S_F2xz_S_bd+CDZ*I_ERI_Gx3y_S_D2x_S_bd;
  Double I_ERI_Gx2yz_S_D2x_Pz_bd = I_ERI_Gx2yz_S_F2xz_S_bd+CDZ*I_ERI_Gx2yz_S_D2x_S_bd;
  Double I_ERI_Gxy2z_S_D2x_Pz_bd = I_ERI_Gxy2z_S_F2xz_S_bd+CDZ*I_ERI_Gxy2z_S_D2x_S_bd;
  Double I_ERI_Gx3z_S_D2x_Pz_bd = I_ERI_Gx3z_S_F2xz_S_bd+CDZ*I_ERI_Gx3z_S_D2x_S_bd;
  Double I_ERI_G4y_S_D2x_Pz_bd = I_ERI_G4y_S_F2xz_S_bd+CDZ*I_ERI_G4y_S_D2x_S_bd;
  Double I_ERI_G3yz_S_D2x_Pz_bd = I_ERI_G3yz_S_F2xz_S_bd+CDZ*I_ERI_G3yz_S_D2x_S_bd;
  Double I_ERI_G2y2z_S_D2x_Pz_bd = I_ERI_G2y2z_S_F2xz_S_bd+CDZ*I_ERI_G2y2z_S_D2x_S_bd;
  Double I_ERI_Gy3z_S_D2x_Pz_bd = I_ERI_Gy3z_S_F2xz_S_bd+CDZ*I_ERI_Gy3z_S_D2x_S_bd;
  Double I_ERI_G4z_S_D2x_Pz_bd = I_ERI_G4z_S_F2xz_S_bd+CDZ*I_ERI_G4z_S_D2x_S_bd;
  Double I_ERI_G4x_S_Dxy_Pz_bd = I_ERI_G4x_S_Fxyz_S_bd+CDZ*I_ERI_G4x_S_Dxy_S_bd;
  Double I_ERI_G3xy_S_Dxy_Pz_bd = I_ERI_G3xy_S_Fxyz_S_bd+CDZ*I_ERI_G3xy_S_Dxy_S_bd;
  Double I_ERI_G3xz_S_Dxy_Pz_bd = I_ERI_G3xz_S_Fxyz_S_bd+CDZ*I_ERI_G3xz_S_Dxy_S_bd;
  Double I_ERI_G2x2y_S_Dxy_Pz_bd = I_ERI_G2x2y_S_Fxyz_S_bd+CDZ*I_ERI_G2x2y_S_Dxy_S_bd;
  Double I_ERI_G2xyz_S_Dxy_Pz_bd = I_ERI_G2xyz_S_Fxyz_S_bd+CDZ*I_ERI_G2xyz_S_Dxy_S_bd;
  Double I_ERI_G2x2z_S_Dxy_Pz_bd = I_ERI_G2x2z_S_Fxyz_S_bd+CDZ*I_ERI_G2x2z_S_Dxy_S_bd;
  Double I_ERI_Gx3y_S_Dxy_Pz_bd = I_ERI_Gx3y_S_Fxyz_S_bd+CDZ*I_ERI_Gx3y_S_Dxy_S_bd;
  Double I_ERI_Gx2yz_S_Dxy_Pz_bd = I_ERI_Gx2yz_S_Fxyz_S_bd+CDZ*I_ERI_Gx2yz_S_Dxy_S_bd;
  Double I_ERI_Gxy2z_S_Dxy_Pz_bd = I_ERI_Gxy2z_S_Fxyz_S_bd+CDZ*I_ERI_Gxy2z_S_Dxy_S_bd;
  Double I_ERI_Gx3z_S_Dxy_Pz_bd = I_ERI_Gx3z_S_Fxyz_S_bd+CDZ*I_ERI_Gx3z_S_Dxy_S_bd;
  Double I_ERI_G4y_S_Dxy_Pz_bd = I_ERI_G4y_S_Fxyz_S_bd+CDZ*I_ERI_G4y_S_Dxy_S_bd;
  Double I_ERI_G3yz_S_Dxy_Pz_bd = I_ERI_G3yz_S_Fxyz_S_bd+CDZ*I_ERI_G3yz_S_Dxy_S_bd;
  Double I_ERI_G2y2z_S_Dxy_Pz_bd = I_ERI_G2y2z_S_Fxyz_S_bd+CDZ*I_ERI_G2y2z_S_Dxy_S_bd;
  Double I_ERI_Gy3z_S_Dxy_Pz_bd = I_ERI_Gy3z_S_Fxyz_S_bd+CDZ*I_ERI_Gy3z_S_Dxy_S_bd;
  Double I_ERI_G4z_S_Dxy_Pz_bd = I_ERI_G4z_S_Fxyz_S_bd+CDZ*I_ERI_G4z_S_Dxy_S_bd;
  Double I_ERI_G4x_S_Dxz_Pz_bd = I_ERI_G4x_S_Fx2z_S_bd+CDZ*I_ERI_G4x_S_Dxz_S_bd;
  Double I_ERI_G3xy_S_Dxz_Pz_bd = I_ERI_G3xy_S_Fx2z_S_bd+CDZ*I_ERI_G3xy_S_Dxz_S_bd;
  Double I_ERI_G3xz_S_Dxz_Pz_bd = I_ERI_G3xz_S_Fx2z_S_bd+CDZ*I_ERI_G3xz_S_Dxz_S_bd;
  Double I_ERI_G2x2y_S_Dxz_Pz_bd = I_ERI_G2x2y_S_Fx2z_S_bd+CDZ*I_ERI_G2x2y_S_Dxz_S_bd;
  Double I_ERI_G2xyz_S_Dxz_Pz_bd = I_ERI_G2xyz_S_Fx2z_S_bd+CDZ*I_ERI_G2xyz_S_Dxz_S_bd;
  Double I_ERI_G2x2z_S_Dxz_Pz_bd = I_ERI_G2x2z_S_Fx2z_S_bd+CDZ*I_ERI_G2x2z_S_Dxz_S_bd;
  Double I_ERI_Gx3y_S_Dxz_Pz_bd = I_ERI_Gx3y_S_Fx2z_S_bd+CDZ*I_ERI_Gx3y_S_Dxz_S_bd;
  Double I_ERI_Gx2yz_S_Dxz_Pz_bd = I_ERI_Gx2yz_S_Fx2z_S_bd+CDZ*I_ERI_Gx2yz_S_Dxz_S_bd;
  Double I_ERI_Gxy2z_S_Dxz_Pz_bd = I_ERI_Gxy2z_S_Fx2z_S_bd+CDZ*I_ERI_Gxy2z_S_Dxz_S_bd;
  Double I_ERI_Gx3z_S_Dxz_Pz_bd = I_ERI_Gx3z_S_Fx2z_S_bd+CDZ*I_ERI_Gx3z_S_Dxz_S_bd;
  Double I_ERI_G4y_S_Dxz_Pz_bd = I_ERI_G4y_S_Fx2z_S_bd+CDZ*I_ERI_G4y_S_Dxz_S_bd;
  Double I_ERI_G3yz_S_Dxz_Pz_bd = I_ERI_G3yz_S_Fx2z_S_bd+CDZ*I_ERI_G3yz_S_Dxz_S_bd;
  Double I_ERI_G2y2z_S_Dxz_Pz_bd = I_ERI_G2y2z_S_Fx2z_S_bd+CDZ*I_ERI_G2y2z_S_Dxz_S_bd;
  Double I_ERI_Gy3z_S_Dxz_Pz_bd = I_ERI_Gy3z_S_Fx2z_S_bd+CDZ*I_ERI_Gy3z_S_Dxz_S_bd;
  Double I_ERI_G4z_S_Dxz_Pz_bd = I_ERI_G4z_S_Fx2z_S_bd+CDZ*I_ERI_G4z_S_Dxz_S_bd;
  Double I_ERI_G4x_S_D2y_Pz_bd = I_ERI_G4x_S_F2yz_S_bd+CDZ*I_ERI_G4x_S_D2y_S_bd;
  Double I_ERI_G3xy_S_D2y_Pz_bd = I_ERI_G3xy_S_F2yz_S_bd+CDZ*I_ERI_G3xy_S_D2y_S_bd;
  Double I_ERI_G3xz_S_D2y_Pz_bd = I_ERI_G3xz_S_F2yz_S_bd+CDZ*I_ERI_G3xz_S_D2y_S_bd;
  Double I_ERI_G2x2y_S_D2y_Pz_bd = I_ERI_G2x2y_S_F2yz_S_bd+CDZ*I_ERI_G2x2y_S_D2y_S_bd;
  Double I_ERI_G2xyz_S_D2y_Pz_bd = I_ERI_G2xyz_S_F2yz_S_bd+CDZ*I_ERI_G2xyz_S_D2y_S_bd;
  Double I_ERI_G2x2z_S_D2y_Pz_bd = I_ERI_G2x2z_S_F2yz_S_bd+CDZ*I_ERI_G2x2z_S_D2y_S_bd;
  Double I_ERI_Gx3y_S_D2y_Pz_bd = I_ERI_Gx3y_S_F2yz_S_bd+CDZ*I_ERI_Gx3y_S_D2y_S_bd;
  Double I_ERI_Gx2yz_S_D2y_Pz_bd = I_ERI_Gx2yz_S_F2yz_S_bd+CDZ*I_ERI_Gx2yz_S_D2y_S_bd;
  Double I_ERI_Gxy2z_S_D2y_Pz_bd = I_ERI_Gxy2z_S_F2yz_S_bd+CDZ*I_ERI_Gxy2z_S_D2y_S_bd;
  Double I_ERI_Gx3z_S_D2y_Pz_bd = I_ERI_Gx3z_S_F2yz_S_bd+CDZ*I_ERI_Gx3z_S_D2y_S_bd;
  Double I_ERI_G4y_S_D2y_Pz_bd = I_ERI_G4y_S_F2yz_S_bd+CDZ*I_ERI_G4y_S_D2y_S_bd;
  Double I_ERI_G3yz_S_D2y_Pz_bd = I_ERI_G3yz_S_F2yz_S_bd+CDZ*I_ERI_G3yz_S_D2y_S_bd;
  Double I_ERI_G2y2z_S_D2y_Pz_bd = I_ERI_G2y2z_S_F2yz_S_bd+CDZ*I_ERI_G2y2z_S_D2y_S_bd;
  Double I_ERI_Gy3z_S_D2y_Pz_bd = I_ERI_Gy3z_S_F2yz_S_bd+CDZ*I_ERI_Gy3z_S_D2y_S_bd;
  Double I_ERI_G4z_S_D2y_Pz_bd = I_ERI_G4z_S_F2yz_S_bd+CDZ*I_ERI_G4z_S_D2y_S_bd;
  Double I_ERI_G4x_S_Dyz_Pz_bd = I_ERI_G4x_S_Fy2z_S_bd+CDZ*I_ERI_G4x_S_Dyz_S_bd;
  Double I_ERI_G3xy_S_Dyz_Pz_bd = I_ERI_G3xy_S_Fy2z_S_bd+CDZ*I_ERI_G3xy_S_Dyz_S_bd;
  Double I_ERI_G3xz_S_Dyz_Pz_bd = I_ERI_G3xz_S_Fy2z_S_bd+CDZ*I_ERI_G3xz_S_Dyz_S_bd;
  Double I_ERI_G2x2y_S_Dyz_Pz_bd = I_ERI_G2x2y_S_Fy2z_S_bd+CDZ*I_ERI_G2x2y_S_Dyz_S_bd;
  Double I_ERI_G2xyz_S_Dyz_Pz_bd = I_ERI_G2xyz_S_Fy2z_S_bd+CDZ*I_ERI_G2xyz_S_Dyz_S_bd;
  Double I_ERI_G2x2z_S_Dyz_Pz_bd = I_ERI_G2x2z_S_Fy2z_S_bd+CDZ*I_ERI_G2x2z_S_Dyz_S_bd;
  Double I_ERI_Gx3y_S_Dyz_Pz_bd = I_ERI_Gx3y_S_Fy2z_S_bd+CDZ*I_ERI_Gx3y_S_Dyz_S_bd;
  Double I_ERI_Gx2yz_S_Dyz_Pz_bd = I_ERI_Gx2yz_S_Fy2z_S_bd+CDZ*I_ERI_Gx2yz_S_Dyz_S_bd;
  Double I_ERI_Gxy2z_S_Dyz_Pz_bd = I_ERI_Gxy2z_S_Fy2z_S_bd+CDZ*I_ERI_Gxy2z_S_Dyz_S_bd;
  Double I_ERI_Gx3z_S_Dyz_Pz_bd = I_ERI_Gx3z_S_Fy2z_S_bd+CDZ*I_ERI_Gx3z_S_Dyz_S_bd;
  Double I_ERI_G4y_S_Dyz_Pz_bd = I_ERI_G4y_S_Fy2z_S_bd+CDZ*I_ERI_G4y_S_Dyz_S_bd;
  Double I_ERI_G3yz_S_Dyz_Pz_bd = I_ERI_G3yz_S_Fy2z_S_bd+CDZ*I_ERI_G3yz_S_Dyz_S_bd;
  Double I_ERI_G2y2z_S_Dyz_Pz_bd = I_ERI_G2y2z_S_Fy2z_S_bd+CDZ*I_ERI_G2y2z_S_Dyz_S_bd;
  Double I_ERI_Gy3z_S_Dyz_Pz_bd = I_ERI_Gy3z_S_Fy2z_S_bd+CDZ*I_ERI_Gy3z_S_Dyz_S_bd;
  Double I_ERI_G4z_S_Dyz_Pz_bd = I_ERI_G4z_S_Fy2z_S_bd+CDZ*I_ERI_G4z_S_Dyz_S_bd;
  Double I_ERI_G4x_S_D2z_Pz_bd = I_ERI_G4x_S_F3z_S_bd+CDZ*I_ERI_G4x_S_D2z_S_bd;
  Double I_ERI_G3xy_S_D2z_Pz_bd = I_ERI_G3xy_S_F3z_S_bd+CDZ*I_ERI_G3xy_S_D2z_S_bd;
  Double I_ERI_G3xz_S_D2z_Pz_bd = I_ERI_G3xz_S_F3z_S_bd+CDZ*I_ERI_G3xz_S_D2z_S_bd;
  Double I_ERI_G2x2y_S_D2z_Pz_bd = I_ERI_G2x2y_S_F3z_S_bd+CDZ*I_ERI_G2x2y_S_D2z_S_bd;
  Double I_ERI_G2xyz_S_D2z_Pz_bd = I_ERI_G2xyz_S_F3z_S_bd+CDZ*I_ERI_G2xyz_S_D2z_S_bd;
  Double I_ERI_G2x2z_S_D2z_Pz_bd = I_ERI_G2x2z_S_F3z_S_bd+CDZ*I_ERI_G2x2z_S_D2z_S_bd;
  Double I_ERI_Gx3y_S_D2z_Pz_bd = I_ERI_Gx3y_S_F3z_S_bd+CDZ*I_ERI_Gx3y_S_D2z_S_bd;
  Double I_ERI_Gx2yz_S_D2z_Pz_bd = I_ERI_Gx2yz_S_F3z_S_bd+CDZ*I_ERI_Gx2yz_S_D2z_S_bd;
  Double I_ERI_Gxy2z_S_D2z_Pz_bd = I_ERI_Gxy2z_S_F3z_S_bd+CDZ*I_ERI_Gxy2z_S_D2z_S_bd;
  Double I_ERI_Gx3z_S_D2z_Pz_bd = I_ERI_Gx3z_S_F3z_S_bd+CDZ*I_ERI_Gx3z_S_D2z_S_bd;
  Double I_ERI_G4y_S_D2z_Pz_bd = I_ERI_G4y_S_F3z_S_bd+CDZ*I_ERI_G4y_S_D2z_S_bd;
  Double I_ERI_G3yz_S_D2z_Pz_bd = I_ERI_G3yz_S_F3z_S_bd+CDZ*I_ERI_G3yz_S_D2z_S_bd;
  Double I_ERI_G2y2z_S_D2z_Pz_bd = I_ERI_G2y2z_S_F3z_S_bd+CDZ*I_ERI_G2y2z_S_D2z_S_bd;
  Double I_ERI_Gy3z_S_D2z_Pz_bd = I_ERI_Gy3z_S_F3z_S_bd+CDZ*I_ERI_Gy3z_S_D2z_S_bd;
  Double I_ERI_G4z_S_D2z_Pz_bd = I_ERI_G4z_S_F3z_S_bd+CDZ*I_ERI_G4z_S_D2z_S_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_F_P_cd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cd
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cd
   ************************************************************/
  Double I_ERI_F3x_S_F3x_Px_cd = I_ERI_F3x_S_G4x_S_cd+CDX*I_ERI_F3x_S_F3x_S_cd;
  Double I_ERI_F2xy_S_F3x_Px_cd = I_ERI_F2xy_S_G4x_S_cd+CDX*I_ERI_F2xy_S_F3x_S_cd;
  Double I_ERI_F2xz_S_F3x_Px_cd = I_ERI_F2xz_S_G4x_S_cd+CDX*I_ERI_F2xz_S_F3x_S_cd;
  Double I_ERI_Fx2y_S_F3x_Px_cd = I_ERI_Fx2y_S_G4x_S_cd+CDX*I_ERI_Fx2y_S_F3x_S_cd;
  Double I_ERI_Fxyz_S_F3x_Px_cd = I_ERI_Fxyz_S_G4x_S_cd+CDX*I_ERI_Fxyz_S_F3x_S_cd;
  Double I_ERI_Fx2z_S_F3x_Px_cd = I_ERI_Fx2z_S_G4x_S_cd+CDX*I_ERI_Fx2z_S_F3x_S_cd;
  Double I_ERI_F3y_S_F3x_Px_cd = I_ERI_F3y_S_G4x_S_cd+CDX*I_ERI_F3y_S_F3x_S_cd;
  Double I_ERI_F2yz_S_F3x_Px_cd = I_ERI_F2yz_S_G4x_S_cd+CDX*I_ERI_F2yz_S_F3x_S_cd;
  Double I_ERI_Fy2z_S_F3x_Px_cd = I_ERI_Fy2z_S_G4x_S_cd+CDX*I_ERI_Fy2z_S_F3x_S_cd;
  Double I_ERI_F3z_S_F3x_Px_cd = I_ERI_F3z_S_G4x_S_cd+CDX*I_ERI_F3z_S_F3x_S_cd;
  Double I_ERI_F3x_S_F2xy_Px_cd = I_ERI_F3x_S_G3xy_S_cd+CDX*I_ERI_F3x_S_F2xy_S_cd;
  Double I_ERI_F2xy_S_F2xy_Px_cd = I_ERI_F2xy_S_G3xy_S_cd+CDX*I_ERI_F2xy_S_F2xy_S_cd;
  Double I_ERI_F2xz_S_F2xy_Px_cd = I_ERI_F2xz_S_G3xy_S_cd+CDX*I_ERI_F2xz_S_F2xy_S_cd;
  Double I_ERI_Fx2y_S_F2xy_Px_cd = I_ERI_Fx2y_S_G3xy_S_cd+CDX*I_ERI_Fx2y_S_F2xy_S_cd;
  Double I_ERI_Fxyz_S_F2xy_Px_cd = I_ERI_Fxyz_S_G3xy_S_cd+CDX*I_ERI_Fxyz_S_F2xy_S_cd;
  Double I_ERI_Fx2z_S_F2xy_Px_cd = I_ERI_Fx2z_S_G3xy_S_cd+CDX*I_ERI_Fx2z_S_F2xy_S_cd;
  Double I_ERI_F3y_S_F2xy_Px_cd = I_ERI_F3y_S_G3xy_S_cd+CDX*I_ERI_F3y_S_F2xy_S_cd;
  Double I_ERI_F2yz_S_F2xy_Px_cd = I_ERI_F2yz_S_G3xy_S_cd+CDX*I_ERI_F2yz_S_F2xy_S_cd;
  Double I_ERI_Fy2z_S_F2xy_Px_cd = I_ERI_Fy2z_S_G3xy_S_cd+CDX*I_ERI_Fy2z_S_F2xy_S_cd;
  Double I_ERI_F3z_S_F2xy_Px_cd = I_ERI_F3z_S_G3xy_S_cd+CDX*I_ERI_F3z_S_F2xy_S_cd;
  Double I_ERI_F3x_S_F2xz_Px_cd = I_ERI_F3x_S_G3xz_S_cd+CDX*I_ERI_F3x_S_F2xz_S_cd;
  Double I_ERI_F2xy_S_F2xz_Px_cd = I_ERI_F2xy_S_G3xz_S_cd+CDX*I_ERI_F2xy_S_F2xz_S_cd;
  Double I_ERI_F2xz_S_F2xz_Px_cd = I_ERI_F2xz_S_G3xz_S_cd+CDX*I_ERI_F2xz_S_F2xz_S_cd;
  Double I_ERI_Fx2y_S_F2xz_Px_cd = I_ERI_Fx2y_S_G3xz_S_cd+CDX*I_ERI_Fx2y_S_F2xz_S_cd;
  Double I_ERI_Fxyz_S_F2xz_Px_cd = I_ERI_Fxyz_S_G3xz_S_cd+CDX*I_ERI_Fxyz_S_F2xz_S_cd;
  Double I_ERI_Fx2z_S_F2xz_Px_cd = I_ERI_Fx2z_S_G3xz_S_cd+CDX*I_ERI_Fx2z_S_F2xz_S_cd;
  Double I_ERI_F3y_S_F2xz_Px_cd = I_ERI_F3y_S_G3xz_S_cd+CDX*I_ERI_F3y_S_F2xz_S_cd;
  Double I_ERI_F2yz_S_F2xz_Px_cd = I_ERI_F2yz_S_G3xz_S_cd+CDX*I_ERI_F2yz_S_F2xz_S_cd;
  Double I_ERI_Fy2z_S_F2xz_Px_cd = I_ERI_Fy2z_S_G3xz_S_cd+CDX*I_ERI_Fy2z_S_F2xz_S_cd;
  Double I_ERI_F3z_S_F2xz_Px_cd = I_ERI_F3z_S_G3xz_S_cd+CDX*I_ERI_F3z_S_F2xz_S_cd;
  Double I_ERI_F3x_S_Fx2y_Px_cd = I_ERI_F3x_S_G2x2y_S_cd+CDX*I_ERI_F3x_S_Fx2y_S_cd;
  Double I_ERI_F2xy_S_Fx2y_Px_cd = I_ERI_F2xy_S_G2x2y_S_cd+CDX*I_ERI_F2xy_S_Fx2y_S_cd;
  Double I_ERI_F2xz_S_Fx2y_Px_cd = I_ERI_F2xz_S_G2x2y_S_cd+CDX*I_ERI_F2xz_S_Fx2y_S_cd;
  Double I_ERI_Fx2y_S_Fx2y_Px_cd = I_ERI_Fx2y_S_G2x2y_S_cd+CDX*I_ERI_Fx2y_S_Fx2y_S_cd;
  Double I_ERI_Fxyz_S_Fx2y_Px_cd = I_ERI_Fxyz_S_G2x2y_S_cd+CDX*I_ERI_Fxyz_S_Fx2y_S_cd;
  Double I_ERI_Fx2z_S_Fx2y_Px_cd = I_ERI_Fx2z_S_G2x2y_S_cd+CDX*I_ERI_Fx2z_S_Fx2y_S_cd;
  Double I_ERI_F3y_S_Fx2y_Px_cd = I_ERI_F3y_S_G2x2y_S_cd+CDX*I_ERI_F3y_S_Fx2y_S_cd;
  Double I_ERI_F2yz_S_Fx2y_Px_cd = I_ERI_F2yz_S_G2x2y_S_cd+CDX*I_ERI_F2yz_S_Fx2y_S_cd;
  Double I_ERI_Fy2z_S_Fx2y_Px_cd = I_ERI_Fy2z_S_G2x2y_S_cd+CDX*I_ERI_Fy2z_S_Fx2y_S_cd;
  Double I_ERI_F3z_S_Fx2y_Px_cd = I_ERI_F3z_S_G2x2y_S_cd+CDX*I_ERI_F3z_S_Fx2y_S_cd;
  Double I_ERI_F3x_S_Fxyz_Px_cd = I_ERI_F3x_S_G2xyz_S_cd+CDX*I_ERI_F3x_S_Fxyz_S_cd;
  Double I_ERI_F2xy_S_Fxyz_Px_cd = I_ERI_F2xy_S_G2xyz_S_cd+CDX*I_ERI_F2xy_S_Fxyz_S_cd;
  Double I_ERI_F2xz_S_Fxyz_Px_cd = I_ERI_F2xz_S_G2xyz_S_cd+CDX*I_ERI_F2xz_S_Fxyz_S_cd;
  Double I_ERI_Fx2y_S_Fxyz_Px_cd = I_ERI_Fx2y_S_G2xyz_S_cd+CDX*I_ERI_Fx2y_S_Fxyz_S_cd;
  Double I_ERI_Fxyz_S_Fxyz_Px_cd = I_ERI_Fxyz_S_G2xyz_S_cd+CDX*I_ERI_Fxyz_S_Fxyz_S_cd;
  Double I_ERI_Fx2z_S_Fxyz_Px_cd = I_ERI_Fx2z_S_G2xyz_S_cd+CDX*I_ERI_Fx2z_S_Fxyz_S_cd;
  Double I_ERI_F3y_S_Fxyz_Px_cd = I_ERI_F3y_S_G2xyz_S_cd+CDX*I_ERI_F3y_S_Fxyz_S_cd;
  Double I_ERI_F2yz_S_Fxyz_Px_cd = I_ERI_F2yz_S_G2xyz_S_cd+CDX*I_ERI_F2yz_S_Fxyz_S_cd;
  Double I_ERI_Fy2z_S_Fxyz_Px_cd = I_ERI_Fy2z_S_G2xyz_S_cd+CDX*I_ERI_Fy2z_S_Fxyz_S_cd;
  Double I_ERI_F3z_S_Fxyz_Px_cd = I_ERI_F3z_S_G2xyz_S_cd+CDX*I_ERI_F3z_S_Fxyz_S_cd;
  Double I_ERI_F3x_S_Fx2z_Px_cd = I_ERI_F3x_S_G2x2z_S_cd+CDX*I_ERI_F3x_S_Fx2z_S_cd;
  Double I_ERI_F2xy_S_Fx2z_Px_cd = I_ERI_F2xy_S_G2x2z_S_cd+CDX*I_ERI_F2xy_S_Fx2z_S_cd;
  Double I_ERI_F2xz_S_Fx2z_Px_cd = I_ERI_F2xz_S_G2x2z_S_cd+CDX*I_ERI_F2xz_S_Fx2z_S_cd;
  Double I_ERI_Fx2y_S_Fx2z_Px_cd = I_ERI_Fx2y_S_G2x2z_S_cd+CDX*I_ERI_Fx2y_S_Fx2z_S_cd;
  Double I_ERI_Fxyz_S_Fx2z_Px_cd = I_ERI_Fxyz_S_G2x2z_S_cd+CDX*I_ERI_Fxyz_S_Fx2z_S_cd;
  Double I_ERI_Fx2z_S_Fx2z_Px_cd = I_ERI_Fx2z_S_G2x2z_S_cd+CDX*I_ERI_Fx2z_S_Fx2z_S_cd;
  Double I_ERI_F3y_S_Fx2z_Px_cd = I_ERI_F3y_S_G2x2z_S_cd+CDX*I_ERI_F3y_S_Fx2z_S_cd;
  Double I_ERI_F2yz_S_Fx2z_Px_cd = I_ERI_F2yz_S_G2x2z_S_cd+CDX*I_ERI_F2yz_S_Fx2z_S_cd;
  Double I_ERI_Fy2z_S_Fx2z_Px_cd = I_ERI_Fy2z_S_G2x2z_S_cd+CDX*I_ERI_Fy2z_S_Fx2z_S_cd;
  Double I_ERI_F3z_S_Fx2z_Px_cd = I_ERI_F3z_S_G2x2z_S_cd+CDX*I_ERI_F3z_S_Fx2z_S_cd;
  Double I_ERI_F3x_S_F3y_Px_cd = I_ERI_F3x_S_Gx3y_S_cd+CDX*I_ERI_F3x_S_F3y_S_cd;
  Double I_ERI_F2xy_S_F3y_Px_cd = I_ERI_F2xy_S_Gx3y_S_cd+CDX*I_ERI_F2xy_S_F3y_S_cd;
  Double I_ERI_F2xz_S_F3y_Px_cd = I_ERI_F2xz_S_Gx3y_S_cd+CDX*I_ERI_F2xz_S_F3y_S_cd;
  Double I_ERI_Fx2y_S_F3y_Px_cd = I_ERI_Fx2y_S_Gx3y_S_cd+CDX*I_ERI_Fx2y_S_F3y_S_cd;
  Double I_ERI_Fxyz_S_F3y_Px_cd = I_ERI_Fxyz_S_Gx3y_S_cd+CDX*I_ERI_Fxyz_S_F3y_S_cd;
  Double I_ERI_Fx2z_S_F3y_Px_cd = I_ERI_Fx2z_S_Gx3y_S_cd+CDX*I_ERI_Fx2z_S_F3y_S_cd;
  Double I_ERI_F3y_S_F3y_Px_cd = I_ERI_F3y_S_Gx3y_S_cd+CDX*I_ERI_F3y_S_F3y_S_cd;
  Double I_ERI_F2yz_S_F3y_Px_cd = I_ERI_F2yz_S_Gx3y_S_cd+CDX*I_ERI_F2yz_S_F3y_S_cd;
  Double I_ERI_Fy2z_S_F3y_Px_cd = I_ERI_Fy2z_S_Gx3y_S_cd+CDX*I_ERI_Fy2z_S_F3y_S_cd;
  Double I_ERI_F3z_S_F3y_Px_cd = I_ERI_F3z_S_Gx3y_S_cd+CDX*I_ERI_F3z_S_F3y_S_cd;
  Double I_ERI_F3x_S_F2yz_Px_cd = I_ERI_F3x_S_Gx2yz_S_cd+CDX*I_ERI_F3x_S_F2yz_S_cd;
  Double I_ERI_F2xy_S_F2yz_Px_cd = I_ERI_F2xy_S_Gx2yz_S_cd+CDX*I_ERI_F2xy_S_F2yz_S_cd;
  Double I_ERI_F2xz_S_F2yz_Px_cd = I_ERI_F2xz_S_Gx2yz_S_cd+CDX*I_ERI_F2xz_S_F2yz_S_cd;
  Double I_ERI_Fx2y_S_F2yz_Px_cd = I_ERI_Fx2y_S_Gx2yz_S_cd+CDX*I_ERI_Fx2y_S_F2yz_S_cd;
  Double I_ERI_Fxyz_S_F2yz_Px_cd = I_ERI_Fxyz_S_Gx2yz_S_cd+CDX*I_ERI_Fxyz_S_F2yz_S_cd;
  Double I_ERI_Fx2z_S_F2yz_Px_cd = I_ERI_Fx2z_S_Gx2yz_S_cd+CDX*I_ERI_Fx2z_S_F2yz_S_cd;
  Double I_ERI_F3y_S_F2yz_Px_cd = I_ERI_F3y_S_Gx2yz_S_cd+CDX*I_ERI_F3y_S_F2yz_S_cd;
  Double I_ERI_F2yz_S_F2yz_Px_cd = I_ERI_F2yz_S_Gx2yz_S_cd+CDX*I_ERI_F2yz_S_F2yz_S_cd;
  Double I_ERI_Fy2z_S_F2yz_Px_cd = I_ERI_Fy2z_S_Gx2yz_S_cd+CDX*I_ERI_Fy2z_S_F2yz_S_cd;
  Double I_ERI_F3z_S_F2yz_Px_cd = I_ERI_F3z_S_Gx2yz_S_cd+CDX*I_ERI_F3z_S_F2yz_S_cd;
  Double I_ERI_F3x_S_Fy2z_Px_cd = I_ERI_F3x_S_Gxy2z_S_cd+CDX*I_ERI_F3x_S_Fy2z_S_cd;
  Double I_ERI_F2xy_S_Fy2z_Px_cd = I_ERI_F2xy_S_Gxy2z_S_cd+CDX*I_ERI_F2xy_S_Fy2z_S_cd;
  Double I_ERI_F2xz_S_Fy2z_Px_cd = I_ERI_F2xz_S_Gxy2z_S_cd+CDX*I_ERI_F2xz_S_Fy2z_S_cd;
  Double I_ERI_Fx2y_S_Fy2z_Px_cd = I_ERI_Fx2y_S_Gxy2z_S_cd+CDX*I_ERI_Fx2y_S_Fy2z_S_cd;
  Double I_ERI_Fxyz_S_Fy2z_Px_cd = I_ERI_Fxyz_S_Gxy2z_S_cd+CDX*I_ERI_Fxyz_S_Fy2z_S_cd;
  Double I_ERI_Fx2z_S_Fy2z_Px_cd = I_ERI_Fx2z_S_Gxy2z_S_cd+CDX*I_ERI_Fx2z_S_Fy2z_S_cd;
  Double I_ERI_F3y_S_Fy2z_Px_cd = I_ERI_F3y_S_Gxy2z_S_cd+CDX*I_ERI_F3y_S_Fy2z_S_cd;
  Double I_ERI_F2yz_S_Fy2z_Px_cd = I_ERI_F2yz_S_Gxy2z_S_cd+CDX*I_ERI_F2yz_S_Fy2z_S_cd;
  Double I_ERI_Fy2z_S_Fy2z_Px_cd = I_ERI_Fy2z_S_Gxy2z_S_cd+CDX*I_ERI_Fy2z_S_Fy2z_S_cd;
  Double I_ERI_F3z_S_Fy2z_Px_cd = I_ERI_F3z_S_Gxy2z_S_cd+CDX*I_ERI_F3z_S_Fy2z_S_cd;
  Double I_ERI_F3x_S_F3z_Px_cd = I_ERI_F3x_S_Gx3z_S_cd+CDX*I_ERI_F3x_S_F3z_S_cd;
  Double I_ERI_F2xy_S_F3z_Px_cd = I_ERI_F2xy_S_Gx3z_S_cd+CDX*I_ERI_F2xy_S_F3z_S_cd;
  Double I_ERI_F2xz_S_F3z_Px_cd = I_ERI_F2xz_S_Gx3z_S_cd+CDX*I_ERI_F2xz_S_F3z_S_cd;
  Double I_ERI_Fx2y_S_F3z_Px_cd = I_ERI_Fx2y_S_Gx3z_S_cd+CDX*I_ERI_Fx2y_S_F3z_S_cd;
  Double I_ERI_Fxyz_S_F3z_Px_cd = I_ERI_Fxyz_S_Gx3z_S_cd+CDX*I_ERI_Fxyz_S_F3z_S_cd;
  Double I_ERI_Fx2z_S_F3z_Px_cd = I_ERI_Fx2z_S_Gx3z_S_cd+CDX*I_ERI_Fx2z_S_F3z_S_cd;
  Double I_ERI_F3y_S_F3z_Px_cd = I_ERI_F3y_S_Gx3z_S_cd+CDX*I_ERI_F3y_S_F3z_S_cd;
  Double I_ERI_F2yz_S_F3z_Px_cd = I_ERI_F2yz_S_Gx3z_S_cd+CDX*I_ERI_F2yz_S_F3z_S_cd;
  Double I_ERI_Fy2z_S_F3z_Px_cd = I_ERI_Fy2z_S_Gx3z_S_cd+CDX*I_ERI_Fy2z_S_F3z_S_cd;
  Double I_ERI_F3z_S_F3z_Px_cd = I_ERI_F3z_S_Gx3z_S_cd+CDX*I_ERI_F3z_S_F3z_S_cd;
  Double I_ERI_F3x_S_F3x_Py_cd = I_ERI_F3x_S_G3xy_S_cd+CDY*I_ERI_F3x_S_F3x_S_cd;
  Double I_ERI_F2xy_S_F3x_Py_cd = I_ERI_F2xy_S_G3xy_S_cd+CDY*I_ERI_F2xy_S_F3x_S_cd;
  Double I_ERI_F2xz_S_F3x_Py_cd = I_ERI_F2xz_S_G3xy_S_cd+CDY*I_ERI_F2xz_S_F3x_S_cd;
  Double I_ERI_Fx2y_S_F3x_Py_cd = I_ERI_Fx2y_S_G3xy_S_cd+CDY*I_ERI_Fx2y_S_F3x_S_cd;
  Double I_ERI_Fxyz_S_F3x_Py_cd = I_ERI_Fxyz_S_G3xy_S_cd+CDY*I_ERI_Fxyz_S_F3x_S_cd;
  Double I_ERI_Fx2z_S_F3x_Py_cd = I_ERI_Fx2z_S_G3xy_S_cd+CDY*I_ERI_Fx2z_S_F3x_S_cd;
  Double I_ERI_F3y_S_F3x_Py_cd = I_ERI_F3y_S_G3xy_S_cd+CDY*I_ERI_F3y_S_F3x_S_cd;
  Double I_ERI_F2yz_S_F3x_Py_cd = I_ERI_F2yz_S_G3xy_S_cd+CDY*I_ERI_F2yz_S_F3x_S_cd;
  Double I_ERI_Fy2z_S_F3x_Py_cd = I_ERI_Fy2z_S_G3xy_S_cd+CDY*I_ERI_Fy2z_S_F3x_S_cd;
  Double I_ERI_F3z_S_F3x_Py_cd = I_ERI_F3z_S_G3xy_S_cd+CDY*I_ERI_F3z_S_F3x_S_cd;
  Double I_ERI_F3x_S_F2xy_Py_cd = I_ERI_F3x_S_G2x2y_S_cd+CDY*I_ERI_F3x_S_F2xy_S_cd;
  Double I_ERI_F2xy_S_F2xy_Py_cd = I_ERI_F2xy_S_G2x2y_S_cd+CDY*I_ERI_F2xy_S_F2xy_S_cd;
  Double I_ERI_F2xz_S_F2xy_Py_cd = I_ERI_F2xz_S_G2x2y_S_cd+CDY*I_ERI_F2xz_S_F2xy_S_cd;
  Double I_ERI_Fx2y_S_F2xy_Py_cd = I_ERI_Fx2y_S_G2x2y_S_cd+CDY*I_ERI_Fx2y_S_F2xy_S_cd;
  Double I_ERI_Fxyz_S_F2xy_Py_cd = I_ERI_Fxyz_S_G2x2y_S_cd+CDY*I_ERI_Fxyz_S_F2xy_S_cd;
  Double I_ERI_Fx2z_S_F2xy_Py_cd = I_ERI_Fx2z_S_G2x2y_S_cd+CDY*I_ERI_Fx2z_S_F2xy_S_cd;
  Double I_ERI_F3y_S_F2xy_Py_cd = I_ERI_F3y_S_G2x2y_S_cd+CDY*I_ERI_F3y_S_F2xy_S_cd;
  Double I_ERI_F2yz_S_F2xy_Py_cd = I_ERI_F2yz_S_G2x2y_S_cd+CDY*I_ERI_F2yz_S_F2xy_S_cd;
  Double I_ERI_Fy2z_S_F2xy_Py_cd = I_ERI_Fy2z_S_G2x2y_S_cd+CDY*I_ERI_Fy2z_S_F2xy_S_cd;
  Double I_ERI_F3z_S_F2xy_Py_cd = I_ERI_F3z_S_G2x2y_S_cd+CDY*I_ERI_F3z_S_F2xy_S_cd;
  Double I_ERI_F3x_S_F2xz_Py_cd = I_ERI_F3x_S_G2xyz_S_cd+CDY*I_ERI_F3x_S_F2xz_S_cd;
  Double I_ERI_F2xy_S_F2xz_Py_cd = I_ERI_F2xy_S_G2xyz_S_cd+CDY*I_ERI_F2xy_S_F2xz_S_cd;
  Double I_ERI_F2xz_S_F2xz_Py_cd = I_ERI_F2xz_S_G2xyz_S_cd+CDY*I_ERI_F2xz_S_F2xz_S_cd;
  Double I_ERI_Fx2y_S_F2xz_Py_cd = I_ERI_Fx2y_S_G2xyz_S_cd+CDY*I_ERI_Fx2y_S_F2xz_S_cd;
  Double I_ERI_Fxyz_S_F2xz_Py_cd = I_ERI_Fxyz_S_G2xyz_S_cd+CDY*I_ERI_Fxyz_S_F2xz_S_cd;
  Double I_ERI_Fx2z_S_F2xz_Py_cd = I_ERI_Fx2z_S_G2xyz_S_cd+CDY*I_ERI_Fx2z_S_F2xz_S_cd;
  Double I_ERI_F3y_S_F2xz_Py_cd = I_ERI_F3y_S_G2xyz_S_cd+CDY*I_ERI_F3y_S_F2xz_S_cd;
  Double I_ERI_F2yz_S_F2xz_Py_cd = I_ERI_F2yz_S_G2xyz_S_cd+CDY*I_ERI_F2yz_S_F2xz_S_cd;
  Double I_ERI_Fy2z_S_F2xz_Py_cd = I_ERI_Fy2z_S_G2xyz_S_cd+CDY*I_ERI_Fy2z_S_F2xz_S_cd;
  Double I_ERI_F3z_S_F2xz_Py_cd = I_ERI_F3z_S_G2xyz_S_cd+CDY*I_ERI_F3z_S_F2xz_S_cd;
  Double I_ERI_F3x_S_Fx2y_Py_cd = I_ERI_F3x_S_Gx3y_S_cd+CDY*I_ERI_F3x_S_Fx2y_S_cd;
  Double I_ERI_F2xy_S_Fx2y_Py_cd = I_ERI_F2xy_S_Gx3y_S_cd+CDY*I_ERI_F2xy_S_Fx2y_S_cd;
  Double I_ERI_F2xz_S_Fx2y_Py_cd = I_ERI_F2xz_S_Gx3y_S_cd+CDY*I_ERI_F2xz_S_Fx2y_S_cd;
  Double I_ERI_Fx2y_S_Fx2y_Py_cd = I_ERI_Fx2y_S_Gx3y_S_cd+CDY*I_ERI_Fx2y_S_Fx2y_S_cd;
  Double I_ERI_Fxyz_S_Fx2y_Py_cd = I_ERI_Fxyz_S_Gx3y_S_cd+CDY*I_ERI_Fxyz_S_Fx2y_S_cd;
  Double I_ERI_Fx2z_S_Fx2y_Py_cd = I_ERI_Fx2z_S_Gx3y_S_cd+CDY*I_ERI_Fx2z_S_Fx2y_S_cd;
  Double I_ERI_F3y_S_Fx2y_Py_cd = I_ERI_F3y_S_Gx3y_S_cd+CDY*I_ERI_F3y_S_Fx2y_S_cd;
  Double I_ERI_F2yz_S_Fx2y_Py_cd = I_ERI_F2yz_S_Gx3y_S_cd+CDY*I_ERI_F2yz_S_Fx2y_S_cd;
  Double I_ERI_Fy2z_S_Fx2y_Py_cd = I_ERI_Fy2z_S_Gx3y_S_cd+CDY*I_ERI_Fy2z_S_Fx2y_S_cd;
  Double I_ERI_F3z_S_Fx2y_Py_cd = I_ERI_F3z_S_Gx3y_S_cd+CDY*I_ERI_F3z_S_Fx2y_S_cd;
  Double I_ERI_F3x_S_Fxyz_Py_cd = I_ERI_F3x_S_Gx2yz_S_cd+CDY*I_ERI_F3x_S_Fxyz_S_cd;
  Double I_ERI_F2xy_S_Fxyz_Py_cd = I_ERI_F2xy_S_Gx2yz_S_cd+CDY*I_ERI_F2xy_S_Fxyz_S_cd;
  Double I_ERI_F2xz_S_Fxyz_Py_cd = I_ERI_F2xz_S_Gx2yz_S_cd+CDY*I_ERI_F2xz_S_Fxyz_S_cd;
  Double I_ERI_Fx2y_S_Fxyz_Py_cd = I_ERI_Fx2y_S_Gx2yz_S_cd+CDY*I_ERI_Fx2y_S_Fxyz_S_cd;
  Double I_ERI_Fxyz_S_Fxyz_Py_cd = I_ERI_Fxyz_S_Gx2yz_S_cd+CDY*I_ERI_Fxyz_S_Fxyz_S_cd;
  Double I_ERI_Fx2z_S_Fxyz_Py_cd = I_ERI_Fx2z_S_Gx2yz_S_cd+CDY*I_ERI_Fx2z_S_Fxyz_S_cd;
  Double I_ERI_F3y_S_Fxyz_Py_cd = I_ERI_F3y_S_Gx2yz_S_cd+CDY*I_ERI_F3y_S_Fxyz_S_cd;
  Double I_ERI_F2yz_S_Fxyz_Py_cd = I_ERI_F2yz_S_Gx2yz_S_cd+CDY*I_ERI_F2yz_S_Fxyz_S_cd;
  Double I_ERI_Fy2z_S_Fxyz_Py_cd = I_ERI_Fy2z_S_Gx2yz_S_cd+CDY*I_ERI_Fy2z_S_Fxyz_S_cd;
  Double I_ERI_F3z_S_Fxyz_Py_cd = I_ERI_F3z_S_Gx2yz_S_cd+CDY*I_ERI_F3z_S_Fxyz_S_cd;
  Double I_ERI_F3x_S_Fx2z_Py_cd = I_ERI_F3x_S_Gxy2z_S_cd+CDY*I_ERI_F3x_S_Fx2z_S_cd;
  Double I_ERI_F2xy_S_Fx2z_Py_cd = I_ERI_F2xy_S_Gxy2z_S_cd+CDY*I_ERI_F2xy_S_Fx2z_S_cd;
  Double I_ERI_F2xz_S_Fx2z_Py_cd = I_ERI_F2xz_S_Gxy2z_S_cd+CDY*I_ERI_F2xz_S_Fx2z_S_cd;
  Double I_ERI_Fx2y_S_Fx2z_Py_cd = I_ERI_Fx2y_S_Gxy2z_S_cd+CDY*I_ERI_Fx2y_S_Fx2z_S_cd;
  Double I_ERI_Fxyz_S_Fx2z_Py_cd = I_ERI_Fxyz_S_Gxy2z_S_cd+CDY*I_ERI_Fxyz_S_Fx2z_S_cd;
  Double I_ERI_Fx2z_S_Fx2z_Py_cd = I_ERI_Fx2z_S_Gxy2z_S_cd+CDY*I_ERI_Fx2z_S_Fx2z_S_cd;
  Double I_ERI_F3y_S_Fx2z_Py_cd = I_ERI_F3y_S_Gxy2z_S_cd+CDY*I_ERI_F3y_S_Fx2z_S_cd;
  Double I_ERI_F2yz_S_Fx2z_Py_cd = I_ERI_F2yz_S_Gxy2z_S_cd+CDY*I_ERI_F2yz_S_Fx2z_S_cd;
  Double I_ERI_Fy2z_S_Fx2z_Py_cd = I_ERI_Fy2z_S_Gxy2z_S_cd+CDY*I_ERI_Fy2z_S_Fx2z_S_cd;
  Double I_ERI_F3z_S_Fx2z_Py_cd = I_ERI_F3z_S_Gxy2z_S_cd+CDY*I_ERI_F3z_S_Fx2z_S_cd;
  Double I_ERI_F3x_S_F3y_Py_cd = I_ERI_F3x_S_G4y_S_cd+CDY*I_ERI_F3x_S_F3y_S_cd;
  Double I_ERI_F2xy_S_F3y_Py_cd = I_ERI_F2xy_S_G4y_S_cd+CDY*I_ERI_F2xy_S_F3y_S_cd;
  Double I_ERI_F2xz_S_F3y_Py_cd = I_ERI_F2xz_S_G4y_S_cd+CDY*I_ERI_F2xz_S_F3y_S_cd;
  Double I_ERI_Fx2y_S_F3y_Py_cd = I_ERI_Fx2y_S_G4y_S_cd+CDY*I_ERI_Fx2y_S_F3y_S_cd;
  Double I_ERI_Fxyz_S_F3y_Py_cd = I_ERI_Fxyz_S_G4y_S_cd+CDY*I_ERI_Fxyz_S_F3y_S_cd;
  Double I_ERI_Fx2z_S_F3y_Py_cd = I_ERI_Fx2z_S_G4y_S_cd+CDY*I_ERI_Fx2z_S_F3y_S_cd;
  Double I_ERI_F3y_S_F3y_Py_cd = I_ERI_F3y_S_G4y_S_cd+CDY*I_ERI_F3y_S_F3y_S_cd;
  Double I_ERI_F2yz_S_F3y_Py_cd = I_ERI_F2yz_S_G4y_S_cd+CDY*I_ERI_F2yz_S_F3y_S_cd;
  Double I_ERI_Fy2z_S_F3y_Py_cd = I_ERI_Fy2z_S_G4y_S_cd+CDY*I_ERI_Fy2z_S_F3y_S_cd;
  Double I_ERI_F3z_S_F3y_Py_cd = I_ERI_F3z_S_G4y_S_cd+CDY*I_ERI_F3z_S_F3y_S_cd;
  Double I_ERI_F3x_S_F2yz_Py_cd = I_ERI_F3x_S_G3yz_S_cd+CDY*I_ERI_F3x_S_F2yz_S_cd;
  Double I_ERI_F2xy_S_F2yz_Py_cd = I_ERI_F2xy_S_G3yz_S_cd+CDY*I_ERI_F2xy_S_F2yz_S_cd;
  Double I_ERI_F2xz_S_F2yz_Py_cd = I_ERI_F2xz_S_G3yz_S_cd+CDY*I_ERI_F2xz_S_F2yz_S_cd;
  Double I_ERI_Fx2y_S_F2yz_Py_cd = I_ERI_Fx2y_S_G3yz_S_cd+CDY*I_ERI_Fx2y_S_F2yz_S_cd;
  Double I_ERI_Fxyz_S_F2yz_Py_cd = I_ERI_Fxyz_S_G3yz_S_cd+CDY*I_ERI_Fxyz_S_F2yz_S_cd;
  Double I_ERI_Fx2z_S_F2yz_Py_cd = I_ERI_Fx2z_S_G3yz_S_cd+CDY*I_ERI_Fx2z_S_F2yz_S_cd;
  Double I_ERI_F3y_S_F2yz_Py_cd = I_ERI_F3y_S_G3yz_S_cd+CDY*I_ERI_F3y_S_F2yz_S_cd;
  Double I_ERI_F2yz_S_F2yz_Py_cd = I_ERI_F2yz_S_G3yz_S_cd+CDY*I_ERI_F2yz_S_F2yz_S_cd;
  Double I_ERI_Fy2z_S_F2yz_Py_cd = I_ERI_Fy2z_S_G3yz_S_cd+CDY*I_ERI_Fy2z_S_F2yz_S_cd;
  Double I_ERI_F3z_S_F2yz_Py_cd = I_ERI_F3z_S_G3yz_S_cd+CDY*I_ERI_F3z_S_F2yz_S_cd;
  Double I_ERI_F3x_S_Fy2z_Py_cd = I_ERI_F3x_S_G2y2z_S_cd+CDY*I_ERI_F3x_S_Fy2z_S_cd;
  Double I_ERI_F2xy_S_Fy2z_Py_cd = I_ERI_F2xy_S_G2y2z_S_cd+CDY*I_ERI_F2xy_S_Fy2z_S_cd;
  Double I_ERI_F2xz_S_Fy2z_Py_cd = I_ERI_F2xz_S_G2y2z_S_cd+CDY*I_ERI_F2xz_S_Fy2z_S_cd;
  Double I_ERI_Fx2y_S_Fy2z_Py_cd = I_ERI_Fx2y_S_G2y2z_S_cd+CDY*I_ERI_Fx2y_S_Fy2z_S_cd;
  Double I_ERI_Fxyz_S_Fy2z_Py_cd = I_ERI_Fxyz_S_G2y2z_S_cd+CDY*I_ERI_Fxyz_S_Fy2z_S_cd;
  Double I_ERI_Fx2z_S_Fy2z_Py_cd = I_ERI_Fx2z_S_G2y2z_S_cd+CDY*I_ERI_Fx2z_S_Fy2z_S_cd;
  Double I_ERI_F3y_S_Fy2z_Py_cd = I_ERI_F3y_S_G2y2z_S_cd+CDY*I_ERI_F3y_S_Fy2z_S_cd;
  Double I_ERI_F2yz_S_Fy2z_Py_cd = I_ERI_F2yz_S_G2y2z_S_cd+CDY*I_ERI_F2yz_S_Fy2z_S_cd;
  Double I_ERI_Fy2z_S_Fy2z_Py_cd = I_ERI_Fy2z_S_G2y2z_S_cd+CDY*I_ERI_Fy2z_S_Fy2z_S_cd;
  Double I_ERI_F3z_S_Fy2z_Py_cd = I_ERI_F3z_S_G2y2z_S_cd+CDY*I_ERI_F3z_S_Fy2z_S_cd;
  Double I_ERI_F3x_S_F3z_Py_cd = I_ERI_F3x_S_Gy3z_S_cd+CDY*I_ERI_F3x_S_F3z_S_cd;
  Double I_ERI_F2xy_S_F3z_Py_cd = I_ERI_F2xy_S_Gy3z_S_cd+CDY*I_ERI_F2xy_S_F3z_S_cd;
  Double I_ERI_F2xz_S_F3z_Py_cd = I_ERI_F2xz_S_Gy3z_S_cd+CDY*I_ERI_F2xz_S_F3z_S_cd;
  Double I_ERI_Fx2y_S_F3z_Py_cd = I_ERI_Fx2y_S_Gy3z_S_cd+CDY*I_ERI_Fx2y_S_F3z_S_cd;
  Double I_ERI_Fxyz_S_F3z_Py_cd = I_ERI_Fxyz_S_Gy3z_S_cd+CDY*I_ERI_Fxyz_S_F3z_S_cd;
  Double I_ERI_Fx2z_S_F3z_Py_cd = I_ERI_Fx2z_S_Gy3z_S_cd+CDY*I_ERI_Fx2z_S_F3z_S_cd;
  Double I_ERI_F3y_S_F3z_Py_cd = I_ERI_F3y_S_Gy3z_S_cd+CDY*I_ERI_F3y_S_F3z_S_cd;
  Double I_ERI_F2yz_S_F3z_Py_cd = I_ERI_F2yz_S_Gy3z_S_cd+CDY*I_ERI_F2yz_S_F3z_S_cd;
  Double I_ERI_Fy2z_S_F3z_Py_cd = I_ERI_Fy2z_S_Gy3z_S_cd+CDY*I_ERI_Fy2z_S_F3z_S_cd;
  Double I_ERI_F3z_S_F3z_Py_cd = I_ERI_F3z_S_Gy3z_S_cd+CDY*I_ERI_F3z_S_F3z_S_cd;
  Double I_ERI_F3x_S_F3x_Pz_cd = I_ERI_F3x_S_G3xz_S_cd+CDZ*I_ERI_F3x_S_F3x_S_cd;
  Double I_ERI_F2xy_S_F3x_Pz_cd = I_ERI_F2xy_S_G3xz_S_cd+CDZ*I_ERI_F2xy_S_F3x_S_cd;
  Double I_ERI_F2xz_S_F3x_Pz_cd = I_ERI_F2xz_S_G3xz_S_cd+CDZ*I_ERI_F2xz_S_F3x_S_cd;
  Double I_ERI_Fx2y_S_F3x_Pz_cd = I_ERI_Fx2y_S_G3xz_S_cd+CDZ*I_ERI_Fx2y_S_F3x_S_cd;
  Double I_ERI_Fxyz_S_F3x_Pz_cd = I_ERI_Fxyz_S_G3xz_S_cd+CDZ*I_ERI_Fxyz_S_F3x_S_cd;
  Double I_ERI_Fx2z_S_F3x_Pz_cd = I_ERI_Fx2z_S_G3xz_S_cd+CDZ*I_ERI_Fx2z_S_F3x_S_cd;
  Double I_ERI_F3y_S_F3x_Pz_cd = I_ERI_F3y_S_G3xz_S_cd+CDZ*I_ERI_F3y_S_F3x_S_cd;
  Double I_ERI_F2yz_S_F3x_Pz_cd = I_ERI_F2yz_S_G3xz_S_cd+CDZ*I_ERI_F2yz_S_F3x_S_cd;
  Double I_ERI_Fy2z_S_F3x_Pz_cd = I_ERI_Fy2z_S_G3xz_S_cd+CDZ*I_ERI_Fy2z_S_F3x_S_cd;
  Double I_ERI_F3z_S_F3x_Pz_cd = I_ERI_F3z_S_G3xz_S_cd+CDZ*I_ERI_F3z_S_F3x_S_cd;
  Double I_ERI_F3x_S_F2xy_Pz_cd = I_ERI_F3x_S_G2xyz_S_cd+CDZ*I_ERI_F3x_S_F2xy_S_cd;
  Double I_ERI_F2xy_S_F2xy_Pz_cd = I_ERI_F2xy_S_G2xyz_S_cd+CDZ*I_ERI_F2xy_S_F2xy_S_cd;
  Double I_ERI_F2xz_S_F2xy_Pz_cd = I_ERI_F2xz_S_G2xyz_S_cd+CDZ*I_ERI_F2xz_S_F2xy_S_cd;
  Double I_ERI_Fx2y_S_F2xy_Pz_cd = I_ERI_Fx2y_S_G2xyz_S_cd+CDZ*I_ERI_Fx2y_S_F2xy_S_cd;
  Double I_ERI_Fxyz_S_F2xy_Pz_cd = I_ERI_Fxyz_S_G2xyz_S_cd+CDZ*I_ERI_Fxyz_S_F2xy_S_cd;
  Double I_ERI_Fx2z_S_F2xy_Pz_cd = I_ERI_Fx2z_S_G2xyz_S_cd+CDZ*I_ERI_Fx2z_S_F2xy_S_cd;
  Double I_ERI_F3y_S_F2xy_Pz_cd = I_ERI_F3y_S_G2xyz_S_cd+CDZ*I_ERI_F3y_S_F2xy_S_cd;
  Double I_ERI_F2yz_S_F2xy_Pz_cd = I_ERI_F2yz_S_G2xyz_S_cd+CDZ*I_ERI_F2yz_S_F2xy_S_cd;
  Double I_ERI_Fy2z_S_F2xy_Pz_cd = I_ERI_Fy2z_S_G2xyz_S_cd+CDZ*I_ERI_Fy2z_S_F2xy_S_cd;
  Double I_ERI_F3z_S_F2xy_Pz_cd = I_ERI_F3z_S_G2xyz_S_cd+CDZ*I_ERI_F3z_S_F2xy_S_cd;
  Double I_ERI_F3x_S_F2xz_Pz_cd = I_ERI_F3x_S_G2x2z_S_cd+CDZ*I_ERI_F3x_S_F2xz_S_cd;
  Double I_ERI_F2xy_S_F2xz_Pz_cd = I_ERI_F2xy_S_G2x2z_S_cd+CDZ*I_ERI_F2xy_S_F2xz_S_cd;
  Double I_ERI_F2xz_S_F2xz_Pz_cd = I_ERI_F2xz_S_G2x2z_S_cd+CDZ*I_ERI_F2xz_S_F2xz_S_cd;
  Double I_ERI_Fx2y_S_F2xz_Pz_cd = I_ERI_Fx2y_S_G2x2z_S_cd+CDZ*I_ERI_Fx2y_S_F2xz_S_cd;
  Double I_ERI_Fxyz_S_F2xz_Pz_cd = I_ERI_Fxyz_S_G2x2z_S_cd+CDZ*I_ERI_Fxyz_S_F2xz_S_cd;
  Double I_ERI_Fx2z_S_F2xz_Pz_cd = I_ERI_Fx2z_S_G2x2z_S_cd+CDZ*I_ERI_Fx2z_S_F2xz_S_cd;
  Double I_ERI_F3y_S_F2xz_Pz_cd = I_ERI_F3y_S_G2x2z_S_cd+CDZ*I_ERI_F3y_S_F2xz_S_cd;
  Double I_ERI_F2yz_S_F2xz_Pz_cd = I_ERI_F2yz_S_G2x2z_S_cd+CDZ*I_ERI_F2yz_S_F2xz_S_cd;
  Double I_ERI_Fy2z_S_F2xz_Pz_cd = I_ERI_Fy2z_S_G2x2z_S_cd+CDZ*I_ERI_Fy2z_S_F2xz_S_cd;
  Double I_ERI_F3z_S_F2xz_Pz_cd = I_ERI_F3z_S_G2x2z_S_cd+CDZ*I_ERI_F3z_S_F2xz_S_cd;
  Double I_ERI_F3x_S_Fx2y_Pz_cd = I_ERI_F3x_S_Gx2yz_S_cd+CDZ*I_ERI_F3x_S_Fx2y_S_cd;
  Double I_ERI_F2xy_S_Fx2y_Pz_cd = I_ERI_F2xy_S_Gx2yz_S_cd+CDZ*I_ERI_F2xy_S_Fx2y_S_cd;
  Double I_ERI_F2xz_S_Fx2y_Pz_cd = I_ERI_F2xz_S_Gx2yz_S_cd+CDZ*I_ERI_F2xz_S_Fx2y_S_cd;
  Double I_ERI_Fx2y_S_Fx2y_Pz_cd = I_ERI_Fx2y_S_Gx2yz_S_cd+CDZ*I_ERI_Fx2y_S_Fx2y_S_cd;
  Double I_ERI_Fxyz_S_Fx2y_Pz_cd = I_ERI_Fxyz_S_Gx2yz_S_cd+CDZ*I_ERI_Fxyz_S_Fx2y_S_cd;
  Double I_ERI_Fx2z_S_Fx2y_Pz_cd = I_ERI_Fx2z_S_Gx2yz_S_cd+CDZ*I_ERI_Fx2z_S_Fx2y_S_cd;
  Double I_ERI_F3y_S_Fx2y_Pz_cd = I_ERI_F3y_S_Gx2yz_S_cd+CDZ*I_ERI_F3y_S_Fx2y_S_cd;
  Double I_ERI_F2yz_S_Fx2y_Pz_cd = I_ERI_F2yz_S_Gx2yz_S_cd+CDZ*I_ERI_F2yz_S_Fx2y_S_cd;
  Double I_ERI_Fy2z_S_Fx2y_Pz_cd = I_ERI_Fy2z_S_Gx2yz_S_cd+CDZ*I_ERI_Fy2z_S_Fx2y_S_cd;
  Double I_ERI_F3z_S_Fx2y_Pz_cd = I_ERI_F3z_S_Gx2yz_S_cd+CDZ*I_ERI_F3z_S_Fx2y_S_cd;
  Double I_ERI_F3x_S_Fxyz_Pz_cd = I_ERI_F3x_S_Gxy2z_S_cd+CDZ*I_ERI_F3x_S_Fxyz_S_cd;
  Double I_ERI_F2xy_S_Fxyz_Pz_cd = I_ERI_F2xy_S_Gxy2z_S_cd+CDZ*I_ERI_F2xy_S_Fxyz_S_cd;
  Double I_ERI_F2xz_S_Fxyz_Pz_cd = I_ERI_F2xz_S_Gxy2z_S_cd+CDZ*I_ERI_F2xz_S_Fxyz_S_cd;
  Double I_ERI_Fx2y_S_Fxyz_Pz_cd = I_ERI_Fx2y_S_Gxy2z_S_cd+CDZ*I_ERI_Fx2y_S_Fxyz_S_cd;
  Double I_ERI_Fxyz_S_Fxyz_Pz_cd = I_ERI_Fxyz_S_Gxy2z_S_cd+CDZ*I_ERI_Fxyz_S_Fxyz_S_cd;
  Double I_ERI_Fx2z_S_Fxyz_Pz_cd = I_ERI_Fx2z_S_Gxy2z_S_cd+CDZ*I_ERI_Fx2z_S_Fxyz_S_cd;
  Double I_ERI_F3y_S_Fxyz_Pz_cd = I_ERI_F3y_S_Gxy2z_S_cd+CDZ*I_ERI_F3y_S_Fxyz_S_cd;
  Double I_ERI_F2yz_S_Fxyz_Pz_cd = I_ERI_F2yz_S_Gxy2z_S_cd+CDZ*I_ERI_F2yz_S_Fxyz_S_cd;
  Double I_ERI_Fy2z_S_Fxyz_Pz_cd = I_ERI_Fy2z_S_Gxy2z_S_cd+CDZ*I_ERI_Fy2z_S_Fxyz_S_cd;
  Double I_ERI_F3z_S_Fxyz_Pz_cd = I_ERI_F3z_S_Gxy2z_S_cd+CDZ*I_ERI_F3z_S_Fxyz_S_cd;
  Double I_ERI_F3x_S_Fx2z_Pz_cd = I_ERI_F3x_S_Gx3z_S_cd+CDZ*I_ERI_F3x_S_Fx2z_S_cd;
  Double I_ERI_F2xy_S_Fx2z_Pz_cd = I_ERI_F2xy_S_Gx3z_S_cd+CDZ*I_ERI_F2xy_S_Fx2z_S_cd;
  Double I_ERI_F2xz_S_Fx2z_Pz_cd = I_ERI_F2xz_S_Gx3z_S_cd+CDZ*I_ERI_F2xz_S_Fx2z_S_cd;
  Double I_ERI_Fx2y_S_Fx2z_Pz_cd = I_ERI_Fx2y_S_Gx3z_S_cd+CDZ*I_ERI_Fx2y_S_Fx2z_S_cd;
  Double I_ERI_Fxyz_S_Fx2z_Pz_cd = I_ERI_Fxyz_S_Gx3z_S_cd+CDZ*I_ERI_Fxyz_S_Fx2z_S_cd;
  Double I_ERI_Fx2z_S_Fx2z_Pz_cd = I_ERI_Fx2z_S_Gx3z_S_cd+CDZ*I_ERI_Fx2z_S_Fx2z_S_cd;
  Double I_ERI_F3y_S_Fx2z_Pz_cd = I_ERI_F3y_S_Gx3z_S_cd+CDZ*I_ERI_F3y_S_Fx2z_S_cd;
  Double I_ERI_F2yz_S_Fx2z_Pz_cd = I_ERI_F2yz_S_Gx3z_S_cd+CDZ*I_ERI_F2yz_S_Fx2z_S_cd;
  Double I_ERI_Fy2z_S_Fx2z_Pz_cd = I_ERI_Fy2z_S_Gx3z_S_cd+CDZ*I_ERI_Fy2z_S_Fx2z_S_cd;
  Double I_ERI_F3z_S_Fx2z_Pz_cd = I_ERI_F3z_S_Gx3z_S_cd+CDZ*I_ERI_F3z_S_Fx2z_S_cd;
  Double I_ERI_F3x_S_F3y_Pz_cd = I_ERI_F3x_S_G3yz_S_cd+CDZ*I_ERI_F3x_S_F3y_S_cd;
  Double I_ERI_F2xy_S_F3y_Pz_cd = I_ERI_F2xy_S_G3yz_S_cd+CDZ*I_ERI_F2xy_S_F3y_S_cd;
  Double I_ERI_F2xz_S_F3y_Pz_cd = I_ERI_F2xz_S_G3yz_S_cd+CDZ*I_ERI_F2xz_S_F3y_S_cd;
  Double I_ERI_Fx2y_S_F3y_Pz_cd = I_ERI_Fx2y_S_G3yz_S_cd+CDZ*I_ERI_Fx2y_S_F3y_S_cd;
  Double I_ERI_Fxyz_S_F3y_Pz_cd = I_ERI_Fxyz_S_G3yz_S_cd+CDZ*I_ERI_Fxyz_S_F3y_S_cd;
  Double I_ERI_Fx2z_S_F3y_Pz_cd = I_ERI_Fx2z_S_G3yz_S_cd+CDZ*I_ERI_Fx2z_S_F3y_S_cd;
  Double I_ERI_F3y_S_F3y_Pz_cd = I_ERI_F3y_S_G3yz_S_cd+CDZ*I_ERI_F3y_S_F3y_S_cd;
  Double I_ERI_F2yz_S_F3y_Pz_cd = I_ERI_F2yz_S_G3yz_S_cd+CDZ*I_ERI_F2yz_S_F3y_S_cd;
  Double I_ERI_Fy2z_S_F3y_Pz_cd = I_ERI_Fy2z_S_G3yz_S_cd+CDZ*I_ERI_Fy2z_S_F3y_S_cd;
  Double I_ERI_F3z_S_F3y_Pz_cd = I_ERI_F3z_S_G3yz_S_cd+CDZ*I_ERI_F3z_S_F3y_S_cd;
  Double I_ERI_F3x_S_F2yz_Pz_cd = I_ERI_F3x_S_G2y2z_S_cd+CDZ*I_ERI_F3x_S_F2yz_S_cd;
  Double I_ERI_F2xy_S_F2yz_Pz_cd = I_ERI_F2xy_S_G2y2z_S_cd+CDZ*I_ERI_F2xy_S_F2yz_S_cd;
  Double I_ERI_F2xz_S_F2yz_Pz_cd = I_ERI_F2xz_S_G2y2z_S_cd+CDZ*I_ERI_F2xz_S_F2yz_S_cd;
  Double I_ERI_Fx2y_S_F2yz_Pz_cd = I_ERI_Fx2y_S_G2y2z_S_cd+CDZ*I_ERI_Fx2y_S_F2yz_S_cd;
  Double I_ERI_Fxyz_S_F2yz_Pz_cd = I_ERI_Fxyz_S_G2y2z_S_cd+CDZ*I_ERI_Fxyz_S_F2yz_S_cd;
  Double I_ERI_Fx2z_S_F2yz_Pz_cd = I_ERI_Fx2z_S_G2y2z_S_cd+CDZ*I_ERI_Fx2z_S_F2yz_S_cd;
  Double I_ERI_F3y_S_F2yz_Pz_cd = I_ERI_F3y_S_G2y2z_S_cd+CDZ*I_ERI_F3y_S_F2yz_S_cd;
  Double I_ERI_F2yz_S_F2yz_Pz_cd = I_ERI_F2yz_S_G2y2z_S_cd+CDZ*I_ERI_F2yz_S_F2yz_S_cd;
  Double I_ERI_Fy2z_S_F2yz_Pz_cd = I_ERI_Fy2z_S_G2y2z_S_cd+CDZ*I_ERI_Fy2z_S_F2yz_S_cd;
  Double I_ERI_F3z_S_F2yz_Pz_cd = I_ERI_F3z_S_G2y2z_S_cd+CDZ*I_ERI_F3z_S_F2yz_S_cd;
  Double I_ERI_F3x_S_Fy2z_Pz_cd = I_ERI_F3x_S_Gy3z_S_cd+CDZ*I_ERI_F3x_S_Fy2z_S_cd;
  Double I_ERI_F2xy_S_Fy2z_Pz_cd = I_ERI_F2xy_S_Gy3z_S_cd+CDZ*I_ERI_F2xy_S_Fy2z_S_cd;
  Double I_ERI_F2xz_S_Fy2z_Pz_cd = I_ERI_F2xz_S_Gy3z_S_cd+CDZ*I_ERI_F2xz_S_Fy2z_S_cd;
  Double I_ERI_Fx2y_S_Fy2z_Pz_cd = I_ERI_Fx2y_S_Gy3z_S_cd+CDZ*I_ERI_Fx2y_S_Fy2z_S_cd;
  Double I_ERI_Fxyz_S_Fy2z_Pz_cd = I_ERI_Fxyz_S_Gy3z_S_cd+CDZ*I_ERI_Fxyz_S_Fy2z_S_cd;
  Double I_ERI_Fx2z_S_Fy2z_Pz_cd = I_ERI_Fx2z_S_Gy3z_S_cd+CDZ*I_ERI_Fx2z_S_Fy2z_S_cd;
  Double I_ERI_F3y_S_Fy2z_Pz_cd = I_ERI_F3y_S_Gy3z_S_cd+CDZ*I_ERI_F3y_S_Fy2z_S_cd;
  Double I_ERI_F2yz_S_Fy2z_Pz_cd = I_ERI_F2yz_S_Gy3z_S_cd+CDZ*I_ERI_F2yz_S_Fy2z_S_cd;
  Double I_ERI_Fy2z_S_Fy2z_Pz_cd = I_ERI_Fy2z_S_Gy3z_S_cd+CDZ*I_ERI_Fy2z_S_Fy2z_S_cd;
  Double I_ERI_F3z_S_Fy2z_Pz_cd = I_ERI_F3z_S_Gy3z_S_cd+CDZ*I_ERI_F3z_S_Fy2z_S_cd;
  Double I_ERI_F3x_S_F3z_Pz_cd = I_ERI_F3x_S_G4z_S_cd+CDZ*I_ERI_F3x_S_F3z_S_cd;
  Double I_ERI_F2xy_S_F3z_Pz_cd = I_ERI_F2xy_S_G4z_S_cd+CDZ*I_ERI_F2xy_S_F3z_S_cd;
  Double I_ERI_F2xz_S_F3z_Pz_cd = I_ERI_F2xz_S_G4z_S_cd+CDZ*I_ERI_F2xz_S_F3z_S_cd;
  Double I_ERI_Fx2y_S_F3z_Pz_cd = I_ERI_Fx2y_S_G4z_S_cd+CDZ*I_ERI_Fx2y_S_F3z_S_cd;
  Double I_ERI_Fxyz_S_F3z_Pz_cd = I_ERI_Fxyz_S_G4z_S_cd+CDZ*I_ERI_Fxyz_S_F3z_S_cd;
  Double I_ERI_Fx2z_S_F3z_Pz_cd = I_ERI_Fx2z_S_G4z_S_cd+CDZ*I_ERI_Fx2z_S_F3z_S_cd;
  Double I_ERI_F3y_S_F3z_Pz_cd = I_ERI_F3y_S_G4z_S_cd+CDZ*I_ERI_F3y_S_F3z_S_cd;
  Double I_ERI_F2yz_S_F3z_Pz_cd = I_ERI_F2yz_S_G4z_S_cd+CDZ*I_ERI_F2yz_S_F3z_S_cd;
  Double I_ERI_Fy2z_S_F3z_Pz_cd = I_ERI_Fy2z_S_G4z_S_cd+CDZ*I_ERI_Fy2z_S_F3z_S_cd;
  Double I_ERI_F3z_S_F3z_Pz_cd = I_ERI_F3z_S_G4z_S_cd+CDZ*I_ERI_F3z_S_F3z_S_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_ERI_F3x_S_D2x_Py_dd = I_ERI_F3x_S_F2xy_S_dd+CDY*I_ERI_F3x_S_D2x_S_dd;
  Double I_ERI_F2xy_S_D2x_Py_dd = I_ERI_F2xy_S_F2xy_S_dd+CDY*I_ERI_F2xy_S_D2x_S_dd;
  Double I_ERI_F2xz_S_D2x_Py_dd = I_ERI_F2xz_S_F2xy_S_dd+CDY*I_ERI_F2xz_S_D2x_S_dd;
  Double I_ERI_Fx2y_S_D2x_Py_dd = I_ERI_Fx2y_S_F2xy_S_dd+CDY*I_ERI_Fx2y_S_D2x_S_dd;
  Double I_ERI_Fxyz_S_D2x_Py_dd = I_ERI_Fxyz_S_F2xy_S_dd+CDY*I_ERI_Fxyz_S_D2x_S_dd;
  Double I_ERI_Fx2z_S_D2x_Py_dd = I_ERI_Fx2z_S_F2xy_S_dd+CDY*I_ERI_Fx2z_S_D2x_S_dd;
  Double I_ERI_F3y_S_D2x_Py_dd = I_ERI_F3y_S_F2xy_S_dd+CDY*I_ERI_F3y_S_D2x_S_dd;
  Double I_ERI_F2yz_S_D2x_Py_dd = I_ERI_F2yz_S_F2xy_S_dd+CDY*I_ERI_F2yz_S_D2x_S_dd;
  Double I_ERI_Fy2z_S_D2x_Py_dd = I_ERI_Fy2z_S_F2xy_S_dd+CDY*I_ERI_Fy2z_S_D2x_S_dd;
  Double I_ERI_F3z_S_D2x_Py_dd = I_ERI_F3z_S_F2xy_S_dd+CDY*I_ERI_F3z_S_D2x_S_dd;
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
  Double I_ERI_F3x_S_D2x_Pz_dd = I_ERI_F3x_S_F2xz_S_dd+CDZ*I_ERI_F3x_S_D2x_S_dd;
  Double I_ERI_F2xy_S_D2x_Pz_dd = I_ERI_F2xy_S_F2xz_S_dd+CDZ*I_ERI_F2xy_S_D2x_S_dd;
  Double I_ERI_F2xz_S_D2x_Pz_dd = I_ERI_F2xz_S_F2xz_S_dd+CDZ*I_ERI_F2xz_S_D2x_S_dd;
  Double I_ERI_Fx2y_S_D2x_Pz_dd = I_ERI_Fx2y_S_F2xz_S_dd+CDZ*I_ERI_Fx2y_S_D2x_S_dd;
  Double I_ERI_Fxyz_S_D2x_Pz_dd = I_ERI_Fxyz_S_F2xz_S_dd+CDZ*I_ERI_Fxyz_S_D2x_S_dd;
  Double I_ERI_Fx2z_S_D2x_Pz_dd = I_ERI_Fx2z_S_F2xz_S_dd+CDZ*I_ERI_Fx2z_S_D2x_S_dd;
  Double I_ERI_F3y_S_D2x_Pz_dd = I_ERI_F3y_S_F2xz_S_dd+CDZ*I_ERI_F3y_S_D2x_S_dd;
  Double I_ERI_F2yz_S_D2x_Pz_dd = I_ERI_F2yz_S_F2xz_S_dd+CDZ*I_ERI_F2yz_S_D2x_S_dd;
  Double I_ERI_Fy2z_S_D2x_Pz_dd = I_ERI_Fy2z_S_F2xz_S_dd+CDZ*I_ERI_Fy2z_S_D2x_S_dd;
  Double I_ERI_F3z_S_D2x_Pz_dd = I_ERI_F3z_S_F2xz_S_dd+CDZ*I_ERI_F3z_S_D2x_S_dd;
  Double I_ERI_F3x_S_Dxy_Pz_dd = I_ERI_F3x_S_Fxyz_S_dd+CDZ*I_ERI_F3x_S_Dxy_S_dd;
  Double I_ERI_F2xy_S_Dxy_Pz_dd = I_ERI_F2xy_S_Fxyz_S_dd+CDZ*I_ERI_F2xy_S_Dxy_S_dd;
  Double I_ERI_F2xz_S_Dxy_Pz_dd = I_ERI_F2xz_S_Fxyz_S_dd+CDZ*I_ERI_F2xz_S_Dxy_S_dd;
  Double I_ERI_Fx2y_S_Dxy_Pz_dd = I_ERI_Fx2y_S_Fxyz_S_dd+CDZ*I_ERI_Fx2y_S_Dxy_S_dd;
  Double I_ERI_Fxyz_S_Dxy_Pz_dd = I_ERI_Fxyz_S_Fxyz_S_dd+CDZ*I_ERI_Fxyz_S_Dxy_S_dd;
  Double I_ERI_Fx2z_S_Dxy_Pz_dd = I_ERI_Fx2z_S_Fxyz_S_dd+CDZ*I_ERI_Fx2z_S_Dxy_S_dd;
  Double I_ERI_F3y_S_Dxy_Pz_dd = I_ERI_F3y_S_Fxyz_S_dd+CDZ*I_ERI_F3y_S_Dxy_S_dd;
  Double I_ERI_F2yz_S_Dxy_Pz_dd = I_ERI_F2yz_S_Fxyz_S_dd+CDZ*I_ERI_F2yz_S_Dxy_S_dd;
  Double I_ERI_Fy2z_S_Dxy_Pz_dd = I_ERI_Fy2z_S_Fxyz_S_dd+CDZ*I_ERI_Fy2z_S_Dxy_S_dd;
  Double I_ERI_F3z_S_Dxy_Pz_dd = I_ERI_F3z_S_Fxyz_S_dd+CDZ*I_ERI_F3z_S_Dxy_S_dd;
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
  Double I_ERI_F3x_S_D2y_Pz_dd = I_ERI_F3x_S_F2yz_S_dd+CDZ*I_ERI_F3x_S_D2y_S_dd;
  Double I_ERI_F2xy_S_D2y_Pz_dd = I_ERI_F2xy_S_F2yz_S_dd+CDZ*I_ERI_F2xy_S_D2y_S_dd;
  Double I_ERI_F2xz_S_D2y_Pz_dd = I_ERI_F2xz_S_F2yz_S_dd+CDZ*I_ERI_F2xz_S_D2y_S_dd;
  Double I_ERI_Fx2y_S_D2y_Pz_dd = I_ERI_Fx2y_S_F2yz_S_dd+CDZ*I_ERI_Fx2y_S_D2y_S_dd;
  Double I_ERI_Fxyz_S_D2y_Pz_dd = I_ERI_Fxyz_S_F2yz_S_dd+CDZ*I_ERI_Fxyz_S_D2y_S_dd;
  Double I_ERI_Fx2z_S_D2y_Pz_dd = I_ERI_Fx2z_S_F2yz_S_dd+CDZ*I_ERI_Fx2z_S_D2y_S_dd;
  Double I_ERI_F3y_S_D2y_Pz_dd = I_ERI_F3y_S_F2yz_S_dd+CDZ*I_ERI_F3y_S_D2y_S_dd;
  Double I_ERI_F2yz_S_D2y_Pz_dd = I_ERI_F2yz_S_F2yz_S_dd+CDZ*I_ERI_F2yz_S_D2y_S_dd;
  Double I_ERI_Fy2z_S_D2y_Pz_dd = I_ERI_Fy2z_S_F2yz_S_dd+CDZ*I_ERI_Fy2z_S_D2y_S_dd;
  Double I_ERI_F3z_S_D2y_Pz_dd = I_ERI_F3z_S_F2yz_S_dd+CDZ*I_ERI_F3z_S_D2y_S_dd;
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
   * shell quartet name: SQ_ERI_F_S_F_P_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 50 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_dd
   * RHS shell quartet name: SQ_ERI_F_S_F_S_dd
   ************************************************************/
  Double I_ERI_F3x_S_F3x_Px_dd = I_ERI_F3x_S_G4x_S_dd+CDX*I_ERI_F3x_S_F3x_S_dd;
  Double I_ERI_F2xy_S_F3x_Px_dd = I_ERI_F2xy_S_G4x_S_dd+CDX*I_ERI_F2xy_S_F3x_S_dd;
  Double I_ERI_F2xz_S_F3x_Px_dd = I_ERI_F2xz_S_G4x_S_dd+CDX*I_ERI_F2xz_S_F3x_S_dd;
  Double I_ERI_Fx2y_S_F3x_Px_dd = I_ERI_Fx2y_S_G4x_S_dd+CDX*I_ERI_Fx2y_S_F3x_S_dd;
  Double I_ERI_Fxyz_S_F3x_Px_dd = I_ERI_Fxyz_S_G4x_S_dd+CDX*I_ERI_Fxyz_S_F3x_S_dd;
  Double I_ERI_Fx2z_S_F3x_Px_dd = I_ERI_Fx2z_S_G4x_S_dd+CDX*I_ERI_Fx2z_S_F3x_S_dd;
  Double I_ERI_F3y_S_F3x_Px_dd = I_ERI_F3y_S_G4x_S_dd+CDX*I_ERI_F3y_S_F3x_S_dd;
  Double I_ERI_F2yz_S_F3x_Px_dd = I_ERI_F2yz_S_G4x_S_dd+CDX*I_ERI_F2yz_S_F3x_S_dd;
  Double I_ERI_Fy2z_S_F3x_Px_dd = I_ERI_Fy2z_S_G4x_S_dd+CDX*I_ERI_Fy2z_S_F3x_S_dd;
  Double I_ERI_F3z_S_F3x_Px_dd = I_ERI_F3z_S_G4x_S_dd+CDX*I_ERI_F3z_S_F3x_S_dd;
  Double I_ERI_F3x_S_F2xy_Px_dd = I_ERI_F3x_S_G3xy_S_dd+CDX*I_ERI_F3x_S_F2xy_S_dd;
  Double I_ERI_F2xy_S_F2xy_Px_dd = I_ERI_F2xy_S_G3xy_S_dd+CDX*I_ERI_F2xy_S_F2xy_S_dd;
  Double I_ERI_F2xz_S_F2xy_Px_dd = I_ERI_F2xz_S_G3xy_S_dd+CDX*I_ERI_F2xz_S_F2xy_S_dd;
  Double I_ERI_Fx2y_S_F2xy_Px_dd = I_ERI_Fx2y_S_G3xy_S_dd+CDX*I_ERI_Fx2y_S_F2xy_S_dd;
  Double I_ERI_Fxyz_S_F2xy_Px_dd = I_ERI_Fxyz_S_G3xy_S_dd+CDX*I_ERI_Fxyz_S_F2xy_S_dd;
  Double I_ERI_Fx2z_S_F2xy_Px_dd = I_ERI_Fx2z_S_G3xy_S_dd+CDX*I_ERI_Fx2z_S_F2xy_S_dd;
  Double I_ERI_F3y_S_F2xy_Px_dd = I_ERI_F3y_S_G3xy_S_dd+CDX*I_ERI_F3y_S_F2xy_S_dd;
  Double I_ERI_F2yz_S_F2xy_Px_dd = I_ERI_F2yz_S_G3xy_S_dd+CDX*I_ERI_F2yz_S_F2xy_S_dd;
  Double I_ERI_Fy2z_S_F2xy_Px_dd = I_ERI_Fy2z_S_G3xy_S_dd+CDX*I_ERI_Fy2z_S_F2xy_S_dd;
  Double I_ERI_F3z_S_F2xy_Px_dd = I_ERI_F3z_S_G3xy_S_dd+CDX*I_ERI_F3z_S_F2xy_S_dd;
  Double I_ERI_F3x_S_F2xz_Px_dd = I_ERI_F3x_S_G3xz_S_dd+CDX*I_ERI_F3x_S_F2xz_S_dd;
  Double I_ERI_F2xy_S_F2xz_Px_dd = I_ERI_F2xy_S_G3xz_S_dd+CDX*I_ERI_F2xy_S_F2xz_S_dd;
  Double I_ERI_F2xz_S_F2xz_Px_dd = I_ERI_F2xz_S_G3xz_S_dd+CDX*I_ERI_F2xz_S_F2xz_S_dd;
  Double I_ERI_Fx2y_S_F2xz_Px_dd = I_ERI_Fx2y_S_G3xz_S_dd+CDX*I_ERI_Fx2y_S_F2xz_S_dd;
  Double I_ERI_Fxyz_S_F2xz_Px_dd = I_ERI_Fxyz_S_G3xz_S_dd+CDX*I_ERI_Fxyz_S_F2xz_S_dd;
  Double I_ERI_Fx2z_S_F2xz_Px_dd = I_ERI_Fx2z_S_G3xz_S_dd+CDX*I_ERI_Fx2z_S_F2xz_S_dd;
  Double I_ERI_F3y_S_F2xz_Px_dd = I_ERI_F3y_S_G3xz_S_dd+CDX*I_ERI_F3y_S_F2xz_S_dd;
  Double I_ERI_F2yz_S_F2xz_Px_dd = I_ERI_F2yz_S_G3xz_S_dd+CDX*I_ERI_F2yz_S_F2xz_S_dd;
  Double I_ERI_Fy2z_S_F2xz_Px_dd = I_ERI_Fy2z_S_G3xz_S_dd+CDX*I_ERI_Fy2z_S_F2xz_S_dd;
  Double I_ERI_F3z_S_F2xz_Px_dd = I_ERI_F3z_S_G3xz_S_dd+CDX*I_ERI_F3z_S_F2xz_S_dd;
  Double I_ERI_F3x_S_Fx2y_Px_dd = I_ERI_F3x_S_G2x2y_S_dd+CDX*I_ERI_F3x_S_Fx2y_S_dd;
  Double I_ERI_F2xy_S_Fx2y_Px_dd = I_ERI_F2xy_S_G2x2y_S_dd+CDX*I_ERI_F2xy_S_Fx2y_S_dd;
  Double I_ERI_F2xz_S_Fx2y_Px_dd = I_ERI_F2xz_S_G2x2y_S_dd+CDX*I_ERI_F2xz_S_Fx2y_S_dd;
  Double I_ERI_Fx2y_S_Fx2y_Px_dd = I_ERI_Fx2y_S_G2x2y_S_dd+CDX*I_ERI_Fx2y_S_Fx2y_S_dd;
  Double I_ERI_Fxyz_S_Fx2y_Px_dd = I_ERI_Fxyz_S_G2x2y_S_dd+CDX*I_ERI_Fxyz_S_Fx2y_S_dd;
  Double I_ERI_Fx2z_S_Fx2y_Px_dd = I_ERI_Fx2z_S_G2x2y_S_dd+CDX*I_ERI_Fx2z_S_Fx2y_S_dd;
  Double I_ERI_F3y_S_Fx2y_Px_dd = I_ERI_F3y_S_G2x2y_S_dd+CDX*I_ERI_F3y_S_Fx2y_S_dd;
  Double I_ERI_F2yz_S_Fx2y_Px_dd = I_ERI_F2yz_S_G2x2y_S_dd+CDX*I_ERI_F2yz_S_Fx2y_S_dd;
  Double I_ERI_Fy2z_S_Fx2y_Px_dd = I_ERI_Fy2z_S_G2x2y_S_dd+CDX*I_ERI_Fy2z_S_Fx2y_S_dd;
  Double I_ERI_F3z_S_Fx2y_Px_dd = I_ERI_F3z_S_G2x2y_S_dd+CDX*I_ERI_F3z_S_Fx2y_S_dd;
  Double I_ERI_F3x_S_Fxyz_Px_dd = I_ERI_F3x_S_G2xyz_S_dd+CDX*I_ERI_F3x_S_Fxyz_S_dd;
  Double I_ERI_F2xy_S_Fxyz_Px_dd = I_ERI_F2xy_S_G2xyz_S_dd+CDX*I_ERI_F2xy_S_Fxyz_S_dd;
  Double I_ERI_F2xz_S_Fxyz_Px_dd = I_ERI_F2xz_S_G2xyz_S_dd+CDX*I_ERI_F2xz_S_Fxyz_S_dd;
  Double I_ERI_Fx2y_S_Fxyz_Px_dd = I_ERI_Fx2y_S_G2xyz_S_dd+CDX*I_ERI_Fx2y_S_Fxyz_S_dd;
  Double I_ERI_Fxyz_S_Fxyz_Px_dd = I_ERI_Fxyz_S_G2xyz_S_dd+CDX*I_ERI_Fxyz_S_Fxyz_S_dd;
  Double I_ERI_Fx2z_S_Fxyz_Px_dd = I_ERI_Fx2z_S_G2xyz_S_dd+CDX*I_ERI_Fx2z_S_Fxyz_S_dd;
  Double I_ERI_F3y_S_Fxyz_Px_dd = I_ERI_F3y_S_G2xyz_S_dd+CDX*I_ERI_F3y_S_Fxyz_S_dd;
  Double I_ERI_F2yz_S_Fxyz_Px_dd = I_ERI_F2yz_S_G2xyz_S_dd+CDX*I_ERI_F2yz_S_Fxyz_S_dd;
  Double I_ERI_Fy2z_S_Fxyz_Px_dd = I_ERI_Fy2z_S_G2xyz_S_dd+CDX*I_ERI_Fy2z_S_Fxyz_S_dd;
  Double I_ERI_F3z_S_Fxyz_Px_dd = I_ERI_F3z_S_G2xyz_S_dd+CDX*I_ERI_F3z_S_Fxyz_S_dd;
  Double I_ERI_F3x_S_Fx2z_Px_dd = I_ERI_F3x_S_G2x2z_S_dd+CDX*I_ERI_F3x_S_Fx2z_S_dd;
  Double I_ERI_F2xy_S_Fx2z_Px_dd = I_ERI_F2xy_S_G2x2z_S_dd+CDX*I_ERI_F2xy_S_Fx2z_S_dd;
  Double I_ERI_F2xz_S_Fx2z_Px_dd = I_ERI_F2xz_S_G2x2z_S_dd+CDX*I_ERI_F2xz_S_Fx2z_S_dd;
  Double I_ERI_Fx2y_S_Fx2z_Px_dd = I_ERI_Fx2y_S_G2x2z_S_dd+CDX*I_ERI_Fx2y_S_Fx2z_S_dd;
  Double I_ERI_Fxyz_S_Fx2z_Px_dd = I_ERI_Fxyz_S_G2x2z_S_dd+CDX*I_ERI_Fxyz_S_Fx2z_S_dd;
  Double I_ERI_Fx2z_S_Fx2z_Px_dd = I_ERI_Fx2z_S_G2x2z_S_dd+CDX*I_ERI_Fx2z_S_Fx2z_S_dd;
  Double I_ERI_F3y_S_Fx2z_Px_dd = I_ERI_F3y_S_G2x2z_S_dd+CDX*I_ERI_F3y_S_Fx2z_S_dd;
  Double I_ERI_F2yz_S_Fx2z_Px_dd = I_ERI_F2yz_S_G2x2z_S_dd+CDX*I_ERI_F2yz_S_Fx2z_S_dd;
  Double I_ERI_Fy2z_S_Fx2z_Px_dd = I_ERI_Fy2z_S_G2x2z_S_dd+CDX*I_ERI_Fy2z_S_Fx2z_S_dd;
  Double I_ERI_F3z_S_Fx2z_Px_dd = I_ERI_F3z_S_G2x2z_S_dd+CDX*I_ERI_F3z_S_Fx2z_S_dd;
  Double I_ERI_F3x_S_F3y_Px_dd = I_ERI_F3x_S_Gx3y_S_dd+CDX*I_ERI_F3x_S_F3y_S_dd;
  Double I_ERI_F2xy_S_F3y_Px_dd = I_ERI_F2xy_S_Gx3y_S_dd+CDX*I_ERI_F2xy_S_F3y_S_dd;
  Double I_ERI_F2xz_S_F3y_Px_dd = I_ERI_F2xz_S_Gx3y_S_dd+CDX*I_ERI_F2xz_S_F3y_S_dd;
  Double I_ERI_Fx2y_S_F3y_Px_dd = I_ERI_Fx2y_S_Gx3y_S_dd+CDX*I_ERI_Fx2y_S_F3y_S_dd;
  Double I_ERI_Fxyz_S_F3y_Px_dd = I_ERI_Fxyz_S_Gx3y_S_dd+CDX*I_ERI_Fxyz_S_F3y_S_dd;
  Double I_ERI_Fx2z_S_F3y_Px_dd = I_ERI_Fx2z_S_Gx3y_S_dd+CDX*I_ERI_Fx2z_S_F3y_S_dd;
  Double I_ERI_F3y_S_F3y_Px_dd = I_ERI_F3y_S_Gx3y_S_dd+CDX*I_ERI_F3y_S_F3y_S_dd;
  Double I_ERI_F2yz_S_F3y_Px_dd = I_ERI_F2yz_S_Gx3y_S_dd+CDX*I_ERI_F2yz_S_F3y_S_dd;
  Double I_ERI_Fy2z_S_F3y_Px_dd = I_ERI_Fy2z_S_Gx3y_S_dd+CDX*I_ERI_Fy2z_S_F3y_S_dd;
  Double I_ERI_F3z_S_F3y_Px_dd = I_ERI_F3z_S_Gx3y_S_dd+CDX*I_ERI_F3z_S_F3y_S_dd;
  Double I_ERI_F3x_S_F2yz_Px_dd = I_ERI_F3x_S_Gx2yz_S_dd+CDX*I_ERI_F3x_S_F2yz_S_dd;
  Double I_ERI_F2xy_S_F2yz_Px_dd = I_ERI_F2xy_S_Gx2yz_S_dd+CDX*I_ERI_F2xy_S_F2yz_S_dd;
  Double I_ERI_F2xz_S_F2yz_Px_dd = I_ERI_F2xz_S_Gx2yz_S_dd+CDX*I_ERI_F2xz_S_F2yz_S_dd;
  Double I_ERI_Fx2y_S_F2yz_Px_dd = I_ERI_Fx2y_S_Gx2yz_S_dd+CDX*I_ERI_Fx2y_S_F2yz_S_dd;
  Double I_ERI_Fxyz_S_F2yz_Px_dd = I_ERI_Fxyz_S_Gx2yz_S_dd+CDX*I_ERI_Fxyz_S_F2yz_S_dd;
  Double I_ERI_Fx2z_S_F2yz_Px_dd = I_ERI_Fx2z_S_Gx2yz_S_dd+CDX*I_ERI_Fx2z_S_F2yz_S_dd;
  Double I_ERI_F3y_S_F2yz_Px_dd = I_ERI_F3y_S_Gx2yz_S_dd+CDX*I_ERI_F3y_S_F2yz_S_dd;
  Double I_ERI_F2yz_S_F2yz_Px_dd = I_ERI_F2yz_S_Gx2yz_S_dd+CDX*I_ERI_F2yz_S_F2yz_S_dd;
  Double I_ERI_Fy2z_S_F2yz_Px_dd = I_ERI_Fy2z_S_Gx2yz_S_dd+CDX*I_ERI_Fy2z_S_F2yz_S_dd;
  Double I_ERI_F3z_S_F2yz_Px_dd = I_ERI_F3z_S_Gx2yz_S_dd+CDX*I_ERI_F3z_S_F2yz_S_dd;
  Double I_ERI_F3x_S_Fy2z_Px_dd = I_ERI_F3x_S_Gxy2z_S_dd+CDX*I_ERI_F3x_S_Fy2z_S_dd;
  Double I_ERI_F2xy_S_Fy2z_Px_dd = I_ERI_F2xy_S_Gxy2z_S_dd+CDX*I_ERI_F2xy_S_Fy2z_S_dd;
  Double I_ERI_F2xz_S_Fy2z_Px_dd = I_ERI_F2xz_S_Gxy2z_S_dd+CDX*I_ERI_F2xz_S_Fy2z_S_dd;
  Double I_ERI_Fx2y_S_Fy2z_Px_dd = I_ERI_Fx2y_S_Gxy2z_S_dd+CDX*I_ERI_Fx2y_S_Fy2z_S_dd;
  Double I_ERI_Fxyz_S_Fy2z_Px_dd = I_ERI_Fxyz_S_Gxy2z_S_dd+CDX*I_ERI_Fxyz_S_Fy2z_S_dd;
  Double I_ERI_Fx2z_S_Fy2z_Px_dd = I_ERI_Fx2z_S_Gxy2z_S_dd+CDX*I_ERI_Fx2z_S_Fy2z_S_dd;
  Double I_ERI_F3y_S_Fy2z_Px_dd = I_ERI_F3y_S_Gxy2z_S_dd+CDX*I_ERI_F3y_S_Fy2z_S_dd;
  Double I_ERI_F2yz_S_Fy2z_Px_dd = I_ERI_F2yz_S_Gxy2z_S_dd+CDX*I_ERI_F2yz_S_Fy2z_S_dd;
  Double I_ERI_Fy2z_S_Fy2z_Px_dd = I_ERI_Fy2z_S_Gxy2z_S_dd+CDX*I_ERI_Fy2z_S_Fy2z_S_dd;
  Double I_ERI_F3z_S_Fy2z_Px_dd = I_ERI_F3z_S_Gxy2z_S_dd+CDX*I_ERI_F3z_S_Fy2z_S_dd;
  Double I_ERI_F3x_S_F3z_Px_dd = I_ERI_F3x_S_Gx3z_S_dd+CDX*I_ERI_F3x_S_F3z_S_dd;
  Double I_ERI_F2xy_S_F3z_Px_dd = I_ERI_F2xy_S_Gx3z_S_dd+CDX*I_ERI_F2xy_S_F3z_S_dd;
  Double I_ERI_F2xz_S_F3z_Px_dd = I_ERI_F2xz_S_Gx3z_S_dd+CDX*I_ERI_F2xz_S_F3z_S_dd;
  Double I_ERI_Fx2y_S_F3z_Px_dd = I_ERI_Fx2y_S_Gx3z_S_dd+CDX*I_ERI_Fx2y_S_F3z_S_dd;
  Double I_ERI_Fxyz_S_F3z_Px_dd = I_ERI_Fxyz_S_Gx3z_S_dd+CDX*I_ERI_Fxyz_S_F3z_S_dd;
  Double I_ERI_Fx2z_S_F3z_Px_dd = I_ERI_Fx2z_S_Gx3z_S_dd+CDX*I_ERI_Fx2z_S_F3z_S_dd;
  Double I_ERI_F3y_S_F3z_Px_dd = I_ERI_F3y_S_Gx3z_S_dd+CDX*I_ERI_F3y_S_F3z_S_dd;
  Double I_ERI_F2yz_S_F3z_Px_dd = I_ERI_F2yz_S_Gx3z_S_dd+CDX*I_ERI_F2yz_S_F3z_S_dd;
  Double I_ERI_Fy2z_S_F3z_Px_dd = I_ERI_Fy2z_S_Gx3z_S_dd+CDX*I_ERI_Fy2z_S_F3z_S_dd;
  Double I_ERI_F3z_S_F3z_Px_dd = I_ERI_F3z_S_Gx3z_S_dd+CDX*I_ERI_F3z_S_F3z_S_dd;
  Double I_ERI_F3x_S_F2xy_Py_dd = I_ERI_F3x_S_G2x2y_S_dd+CDY*I_ERI_F3x_S_F2xy_S_dd;
  Double I_ERI_F2xy_S_F2xy_Py_dd = I_ERI_F2xy_S_G2x2y_S_dd+CDY*I_ERI_F2xy_S_F2xy_S_dd;
  Double I_ERI_F2xz_S_F2xy_Py_dd = I_ERI_F2xz_S_G2x2y_S_dd+CDY*I_ERI_F2xz_S_F2xy_S_dd;
  Double I_ERI_Fx2y_S_F2xy_Py_dd = I_ERI_Fx2y_S_G2x2y_S_dd+CDY*I_ERI_Fx2y_S_F2xy_S_dd;
  Double I_ERI_Fxyz_S_F2xy_Py_dd = I_ERI_Fxyz_S_G2x2y_S_dd+CDY*I_ERI_Fxyz_S_F2xy_S_dd;
  Double I_ERI_Fx2z_S_F2xy_Py_dd = I_ERI_Fx2z_S_G2x2y_S_dd+CDY*I_ERI_Fx2z_S_F2xy_S_dd;
  Double I_ERI_F3y_S_F2xy_Py_dd = I_ERI_F3y_S_G2x2y_S_dd+CDY*I_ERI_F3y_S_F2xy_S_dd;
  Double I_ERI_F2yz_S_F2xy_Py_dd = I_ERI_F2yz_S_G2x2y_S_dd+CDY*I_ERI_F2yz_S_F2xy_S_dd;
  Double I_ERI_Fy2z_S_F2xy_Py_dd = I_ERI_Fy2z_S_G2x2y_S_dd+CDY*I_ERI_Fy2z_S_F2xy_S_dd;
  Double I_ERI_F3z_S_F2xy_Py_dd = I_ERI_F3z_S_G2x2y_S_dd+CDY*I_ERI_F3z_S_F2xy_S_dd;
  Double I_ERI_F3x_S_F2xz_Py_dd = I_ERI_F3x_S_G2xyz_S_dd+CDY*I_ERI_F3x_S_F2xz_S_dd;
  Double I_ERI_F2xy_S_F2xz_Py_dd = I_ERI_F2xy_S_G2xyz_S_dd+CDY*I_ERI_F2xy_S_F2xz_S_dd;
  Double I_ERI_F2xz_S_F2xz_Py_dd = I_ERI_F2xz_S_G2xyz_S_dd+CDY*I_ERI_F2xz_S_F2xz_S_dd;
  Double I_ERI_Fx2y_S_F2xz_Py_dd = I_ERI_Fx2y_S_G2xyz_S_dd+CDY*I_ERI_Fx2y_S_F2xz_S_dd;
  Double I_ERI_Fxyz_S_F2xz_Py_dd = I_ERI_Fxyz_S_G2xyz_S_dd+CDY*I_ERI_Fxyz_S_F2xz_S_dd;
  Double I_ERI_Fx2z_S_F2xz_Py_dd = I_ERI_Fx2z_S_G2xyz_S_dd+CDY*I_ERI_Fx2z_S_F2xz_S_dd;
  Double I_ERI_F3y_S_F2xz_Py_dd = I_ERI_F3y_S_G2xyz_S_dd+CDY*I_ERI_F3y_S_F2xz_S_dd;
  Double I_ERI_F2yz_S_F2xz_Py_dd = I_ERI_F2yz_S_G2xyz_S_dd+CDY*I_ERI_F2yz_S_F2xz_S_dd;
  Double I_ERI_Fy2z_S_F2xz_Py_dd = I_ERI_Fy2z_S_G2xyz_S_dd+CDY*I_ERI_Fy2z_S_F2xz_S_dd;
  Double I_ERI_F3z_S_F2xz_Py_dd = I_ERI_F3z_S_G2xyz_S_dd+CDY*I_ERI_F3z_S_F2xz_S_dd;
  Double I_ERI_F3x_S_Fx2y_Py_dd = I_ERI_F3x_S_Gx3y_S_dd+CDY*I_ERI_F3x_S_Fx2y_S_dd;
  Double I_ERI_F2xy_S_Fx2y_Py_dd = I_ERI_F2xy_S_Gx3y_S_dd+CDY*I_ERI_F2xy_S_Fx2y_S_dd;
  Double I_ERI_F2xz_S_Fx2y_Py_dd = I_ERI_F2xz_S_Gx3y_S_dd+CDY*I_ERI_F2xz_S_Fx2y_S_dd;
  Double I_ERI_Fx2y_S_Fx2y_Py_dd = I_ERI_Fx2y_S_Gx3y_S_dd+CDY*I_ERI_Fx2y_S_Fx2y_S_dd;
  Double I_ERI_Fxyz_S_Fx2y_Py_dd = I_ERI_Fxyz_S_Gx3y_S_dd+CDY*I_ERI_Fxyz_S_Fx2y_S_dd;
  Double I_ERI_Fx2z_S_Fx2y_Py_dd = I_ERI_Fx2z_S_Gx3y_S_dd+CDY*I_ERI_Fx2z_S_Fx2y_S_dd;
  Double I_ERI_F3y_S_Fx2y_Py_dd = I_ERI_F3y_S_Gx3y_S_dd+CDY*I_ERI_F3y_S_Fx2y_S_dd;
  Double I_ERI_F2yz_S_Fx2y_Py_dd = I_ERI_F2yz_S_Gx3y_S_dd+CDY*I_ERI_F2yz_S_Fx2y_S_dd;
  Double I_ERI_Fy2z_S_Fx2y_Py_dd = I_ERI_Fy2z_S_Gx3y_S_dd+CDY*I_ERI_Fy2z_S_Fx2y_S_dd;
  Double I_ERI_F3z_S_Fx2y_Py_dd = I_ERI_F3z_S_Gx3y_S_dd+CDY*I_ERI_F3z_S_Fx2y_S_dd;
  Double I_ERI_F3x_S_Fxyz_Py_dd = I_ERI_F3x_S_Gx2yz_S_dd+CDY*I_ERI_F3x_S_Fxyz_S_dd;
  Double I_ERI_F2xy_S_Fxyz_Py_dd = I_ERI_F2xy_S_Gx2yz_S_dd+CDY*I_ERI_F2xy_S_Fxyz_S_dd;
  Double I_ERI_F2xz_S_Fxyz_Py_dd = I_ERI_F2xz_S_Gx2yz_S_dd+CDY*I_ERI_F2xz_S_Fxyz_S_dd;
  Double I_ERI_Fx2y_S_Fxyz_Py_dd = I_ERI_Fx2y_S_Gx2yz_S_dd+CDY*I_ERI_Fx2y_S_Fxyz_S_dd;
  Double I_ERI_Fxyz_S_Fxyz_Py_dd = I_ERI_Fxyz_S_Gx2yz_S_dd+CDY*I_ERI_Fxyz_S_Fxyz_S_dd;
  Double I_ERI_Fx2z_S_Fxyz_Py_dd = I_ERI_Fx2z_S_Gx2yz_S_dd+CDY*I_ERI_Fx2z_S_Fxyz_S_dd;
  Double I_ERI_F3y_S_Fxyz_Py_dd = I_ERI_F3y_S_Gx2yz_S_dd+CDY*I_ERI_F3y_S_Fxyz_S_dd;
  Double I_ERI_F2yz_S_Fxyz_Py_dd = I_ERI_F2yz_S_Gx2yz_S_dd+CDY*I_ERI_F2yz_S_Fxyz_S_dd;
  Double I_ERI_Fy2z_S_Fxyz_Py_dd = I_ERI_Fy2z_S_Gx2yz_S_dd+CDY*I_ERI_Fy2z_S_Fxyz_S_dd;
  Double I_ERI_F3z_S_Fxyz_Py_dd = I_ERI_F3z_S_Gx2yz_S_dd+CDY*I_ERI_F3z_S_Fxyz_S_dd;
  Double I_ERI_F3x_S_Fx2z_Py_dd = I_ERI_F3x_S_Gxy2z_S_dd+CDY*I_ERI_F3x_S_Fx2z_S_dd;
  Double I_ERI_F2xy_S_Fx2z_Py_dd = I_ERI_F2xy_S_Gxy2z_S_dd+CDY*I_ERI_F2xy_S_Fx2z_S_dd;
  Double I_ERI_F2xz_S_Fx2z_Py_dd = I_ERI_F2xz_S_Gxy2z_S_dd+CDY*I_ERI_F2xz_S_Fx2z_S_dd;
  Double I_ERI_Fx2y_S_Fx2z_Py_dd = I_ERI_Fx2y_S_Gxy2z_S_dd+CDY*I_ERI_Fx2y_S_Fx2z_S_dd;
  Double I_ERI_Fxyz_S_Fx2z_Py_dd = I_ERI_Fxyz_S_Gxy2z_S_dd+CDY*I_ERI_Fxyz_S_Fx2z_S_dd;
  Double I_ERI_Fx2z_S_Fx2z_Py_dd = I_ERI_Fx2z_S_Gxy2z_S_dd+CDY*I_ERI_Fx2z_S_Fx2z_S_dd;
  Double I_ERI_F3y_S_Fx2z_Py_dd = I_ERI_F3y_S_Gxy2z_S_dd+CDY*I_ERI_F3y_S_Fx2z_S_dd;
  Double I_ERI_F2yz_S_Fx2z_Py_dd = I_ERI_F2yz_S_Gxy2z_S_dd+CDY*I_ERI_F2yz_S_Fx2z_S_dd;
  Double I_ERI_Fy2z_S_Fx2z_Py_dd = I_ERI_Fy2z_S_Gxy2z_S_dd+CDY*I_ERI_Fy2z_S_Fx2z_S_dd;
  Double I_ERI_F3z_S_Fx2z_Py_dd = I_ERI_F3z_S_Gxy2z_S_dd+CDY*I_ERI_F3z_S_Fx2z_S_dd;
  Double I_ERI_F3x_S_F3y_Py_dd = I_ERI_F3x_S_G4y_S_dd+CDY*I_ERI_F3x_S_F3y_S_dd;
  Double I_ERI_F2xy_S_F3y_Py_dd = I_ERI_F2xy_S_G4y_S_dd+CDY*I_ERI_F2xy_S_F3y_S_dd;
  Double I_ERI_F2xz_S_F3y_Py_dd = I_ERI_F2xz_S_G4y_S_dd+CDY*I_ERI_F2xz_S_F3y_S_dd;
  Double I_ERI_Fx2y_S_F3y_Py_dd = I_ERI_Fx2y_S_G4y_S_dd+CDY*I_ERI_Fx2y_S_F3y_S_dd;
  Double I_ERI_Fxyz_S_F3y_Py_dd = I_ERI_Fxyz_S_G4y_S_dd+CDY*I_ERI_Fxyz_S_F3y_S_dd;
  Double I_ERI_Fx2z_S_F3y_Py_dd = I_ERI_Fx2z_S_G4y_S_dd+CDY*I_ERI_Fx2z_S_F3y_S_dd;
  Double I_ERI_F3y_S_F3y_Py_dd = I_ERI_F3y_S_G4y_S_dd+CDY*I_ERI_F3y_S_F3y_S_dd;
  Double I_ERI_F2yz_S_F3y_Py_dd = I_ERI_F2yz_S_G4y_S_dd+CDY*I_ERI_F2yz_S_F3y_S_dd;
  Double I_ERI_Fy2z_S_F3y_Py_dd = I_ERI_Fy2z_S_G4y_S_dd+CDY*I_ERI_Fy2z_S_F3y_S_dd;
  Double I_ERI_F3z_S_F3y_Py_dd = I_ERI_F3z_S_G4y_S_dd+CDY*I_ERI_F3z_S_F3y_S_dd;
  Double I_ERI_F3x_S_F2yz_Py_dd = I_ERI_F3x_S_G3yz_S_dd+CDY*I_ERI_F3x_S_F2yz_S_dd;
  Double I_ERI_F2xy_S_F2yz_Py_dd = I_ERI_F2xy_S_G3yz_S_dd+CDY*I_ERI_F2xy_S_F2yz_S_dd;
  Double I_ERI_F2xz_S_F2yz_Py_dd = I_ERI_F2xz_S_G3yz_S_dd+CDY*I_ERI_F2xz_S_F2yz_S_dd;
  Double I_ERI_Fx2y_S_F2yz_Py_dd = I_ERI_Fx2y_S_G3yz_S_dd+CDY*I_ERI_Fx2y_S_F2yz_S_dd;
  Double I_ERI_Fxyz_S_F2yz_Py_dd = I_ERI_Fxyz_S_G3yz_S_dd+CDY*I_ERI_Fxyz_S_F2yz_S_dd;
  Double I_ERI_Fx2z_S_F2yz_Py_dd = I_ERI_Fx2z_S_G3yz_S_dd+CDY*I_ERI_Fx2z_S_F2yz_S_dd;
  Double I_ERI_F3y_S_F2yz_Py_dd = I_ERI_F3y_S_G3yz_S_dd+CDY*I_ERI_F3y_S_F2yz_S_dd;
  Double I_ERI_F2yz_S_F2yz_Py_dd = I_ERI_F2yz_S_G3yz_S_dd+CDY*I_ERI_F2yz_S_F2yz_S_dd;
  Double I_ERI_Fy2z_S_F2yz_Py_dd = I_ERI_Fy2z_S_G3yz_S_dd+CDY*I_ERI_Fy2z_S_F2yz_S_dd;
  Double I_ERI_F3z_S_F2yz_Py_dd = I_ERI_F3z_S_G3yz_S_dd+CDY*I_ERI_F3z_S_F2yz_S_dd;
  Double I_ERI_F3x_S_Fy2z_Py_dd = I_ERI_F3x_S_G2y2z_S_dd+CDY*I_ERI_F3x_S_Fy2z_S_dd;
  Double I_ERI_F2xy_S_Fy2z_Py_dd = I_ERI_F2xy_S_G2y2z_S_dd+CDY*I_ERI_F2xy_S_Fy2z_S_dd;
  Double I_ERI_F2xz_S_Fy2z_Py_dd = I_ERI_F2xz_S_G2y2z_S_dd+CDY*I_ERI_F2xz_S_Fy2z_S_dd;
  Double I_ERI_Fx2y_S_Fy2z_Py_dd = I_ERI_Fx2y_S_G2y2z_S_dd+CDY*I_ERI_Fx2y_S_Fy2z_S_dd;
  Double I_ERI_Fxyz_S_Fy2z_Py_dd = I_ERI_Fxyz_S_G2y2z_S_dd+CDY*I_ERI_Fxyz_S_Fy2z_S_dd;
  Double I_ERI_Fx2z_S_Fy2z_Py_dd = I_ERI_Fx2z_S_G2y2z_S_dd+CDY*I_ERI_Fx2z_S_Fy2z_S_dd;
  Double I_ERI_F3y_S_Fy2z_Py_dd = I_ERI_F3y_S_G2y2z_S_dd+CDY*I_ERI_F3y_S_Fy2z_S_dd;
  Double I_ERI_F2yz_S_Fy2z_Py_dd = I_ERI_F2yz_S_G2y2z_S_dd+CDY*I_ERI_F2yz_S_Fy2z_S_dd;
  Double I_ERI_Fy2z_S_Fy2z_Py_dd = I_ERI_Fy2z_S_G2y2z_S_dd+CDY*I_ERI_Fy2z_S_Fy2z_S_dd;
  Double I_ERI_F3z_S_Fy2z_Py_dd = I_ERI_F3z_S_G2y2z_S_dd+CDY*I_ERI_F3z_S_Fy2z_S_dd;
  Double I_ERI_F3x_S_F3z_Py_dd = I_ERI_F3x_S_Gy3z_S_dd+CDY*I_ERI_F3x_S_F3z_S_dd;
  Double I_ERI_F2xy_S_F3z_Py_dd = I_ERI_F2xy_S_Gy3z_S_dd+CDY*I_ERI_F2xy_S_F3z_S_dd;
  Double I_ERI_F2xz_S_F3z_Py_dd = I_ERI_F2xz_S_Gy3z_S_dd+CDY*I_ERI_F2xz_S_F3z_S_dd;
  Double I_ERI_Fx2y_S_F3z_Py_dd = I_ERI_Fx2y_S_Gy3z_S_dd+CDY*I_ERI_Fx2y_S_F3z_S_dd;
  Double I_ERI_Fxyz_S_F3z_Py_dd = I_ERI_Fxyz_S_Gy3z_S_dd+CDY*I_ERI_Fxyz_S_F3z_S_dd;
  Double I_ERI_Fx2z_S_F3z_Py_dd = I_ERI_Fx2z_S_Gy3z_S_dd+CDY*I_ERI_Fx2z_S_F3z_S_dd;
  Double I_ERI_F3y_S_F3z_Py_dd = I_ERI_F3y_S_Gy3z_S_dd+CDY*I_ERI_F3y_S_F3z_S_dd;
  Double I_ERI_F2yz_S_F3z_Py_dd = I_ERI_F2yz_S_Gy3z_S_dd+CDY*I_ERI_F2yz_S_F3z_S_dd;
  Double I_ERI_Fy2z_S_F3z_Py_dd = I_ERI_Fy2z_S_Gy3z_S_dd+CDY*I_ERI_Fy2z_S_F3z_S_dd;
  Double I_ERI_F3z_S_F3z_Py_dd = I_ERI_F3z_S_Gy3z_S_dd+CDY*I_ERI_F3z_S_F3z_S_dd;
  Double I_ERI_F3x_S_F2xz_Pz_dd = I_ERI_F3x_S_G2x2z_S_dd+CDZ*I_ERI_F3x_S_F2xz_S_dd;
  Double I_ERI_F2xy_S_F2xz_Pz_dd = I_ERI_F2xy_S_G2x2z_S_dd+CDZ*I_ERI_F2xy_S_F2xz_S_dd;
  Double I_ERI_F2xz_S_F2xz_Pz_dd = I_ERI_F2xz_S_G2x2z_S_dd+CDZ*I_ERI_F2xz_S_F2xz_S_dd;
  Double I_ERI_Fx2y_S_F2xz_Pz_dd = I_ERI_Fx2y_S_G2x2z_S_dd+CDZ*I_ERI_Fx2y_S_F2xz_S_dd;
  Double I_ERI_Fxyz_S_F2xz_Pz_dd = I_ERI_Fxyz_S_G2x2z_S_dd+CDZ*I_ERI_Fxyz_S_F2xz_S_dd;
  Double I_ERI_Fx2z_S_F2xz_Pz_dd = I_ERI_Fx2z_S_G2x2z_S_dd+CDZ*I_ERI_Fx2z_S_F2xz_S_dd;
  Double I_ERI_F3y_S_F2xz_Pz_dd = I_ERI_F3y_S_G2x2z_S_dd+CDZ*I_ERI_F3y_S_F2xz_S_dd;
  Double I_ERI_F2yz_S_F2xz_Pz_dd = I_ERI_F2yz_S_G2x2z_S_dd+CDZ*I_ERI_F2yz_S_F2xz_S_dd;
  Double I_ERI_Fy2z_S_F2xz_Pz_dd = I_ERI_Fy2z_S_G2x2z_S_dd+CDZ*I_ERI_Fy2z_S_F2xz_S_dd;
  Double I_ERI_F3z_S_F2xz_Pz_dd = I_ERI_F3z_S_G2x2z_S_dd+CDZ*I_ERI_F3z_S_F2xz_S_dd;
  Double I_ERI_F3x_S_Fxyz_Pz_dd = I_ERI_F3x_S_Gxy2z_S_dd+CDZ*I_ERI_F3x_S_Fxyz_S_dd;
  Double I_ERI_F2xy_S_Fxyz_Pz_dd = I_ERI_F2xy_S_Gxy2z_S_dd+CDZ*I_ERI_F2xy_S_Fxyz_S_dd;
  Double I_ERI_F2xz_S_Fxyz_Pz_dd = I_ERI_F2xz_S_Gxy2z_S_dd+CDZ*I_ERI_F2xz_S_Fxyz_S_dd;
  Double I_ERI_Fx2y_S_Fxyz_Pz_dd = I_ERI_Fx2y_S_Gxy2z_S_dd+CDZ*I_ERI_Fx2y_S_Fxyz_S_dd;
  Double I_ERI_Fxyz_S_Fxyz_Pz_dd = I_ERI_Fxyz_S_Gxy2z_S_dd+CDZ*I_ERI_Fxyz_S_Fxyz_S_dd;
  Double I_ERI_Fx2z_S_Fxyz_Pz_dd = I_ERI_Fx2z_S_Gxy2z_S_dd+CDZ*I_ERI_Fx2z_S_Fxyz_S_dd;
  Double I_ERI_F3y_S_Fxyz_Pz_dd = I_ERI_F3y_S_Gxy2z_S_dd+CDZ*I_ERI_F3y_S_Fxyz_S_dd;
  Double I_ERI_F2yz_S_Fxyz_Pz_dd = I_ERI_F2yz_S_Gxy2z_S_dd+CDZ*I_ERI_F2yz_S_Fxyz_S_dd;
  Double I_ERI_Fy2z_S_Fxyz_Pz_dd = I_ERI_Fy2z_S_Gxy2z_S_dd+CDZ*I_ERI_Fy2z_S_Fxyz_S_dd;
  Double I_ERI_F3z_S_Fxyz_Pz_dd = I_ERI_F3z_S_Gxy2z_S_dd+CDZ*I_ERI_F3z_S_Fxyz_S_dd;
  Double I_ERI_F3x_S_Fx2z_Pz_dd = I_ERI_F3x_S_Gx3z_S_dd+CDZ*I_ERI_F3x_S_Fx2z_S_dd;
  Double I_ERI_F2xy_S_Fx2z_Pz_dd = I_ERI_F2xy_S_Gx3z_S_dd+CDZ*I_ERI_F2xy_S_Fx2z_S_dd;
  Double I_ERI_F2xz_S_Fx2z_Pz_dd = I_ERI_F2xz_S_Gx3z_S_dd+CDZ*I_ERI_F2xz_S_Fx2z_S_dd;
  Double I_ERI_Fx2y_S_Fx2z_Pz_dd = I_ERI_Fx2y_S_Gx3z_S_dd+CDZ*I_ERI_Fx2y_S_Fx2z_S_dd;
  Double I_ERI_Fxyz_S_Fx2z_Pz_dd = I_ERI_Fxyz_S_Gx3z_S_dd+CDZ*I_ERI_Fxyz_S_Fx2z_S_dd;
  Double I_ERI_Fx2z_S_Fx2z_Pz_dd = I_ERI_Fx2z_S_Gx3z_S_dd+CDZ*I_ERI_Fx2z_S_Fx2z_S_dd;
  Double I_ERI_F3y_S_Fx2z_Pz_dd = I_ERI_F3y_S_Gx3z_S_dd+CDZ*I_ERI_F3y_S_Fx2z_S_dd;
  Double I_ERI_F2yz_S_Fx2z_Pz_dd = I_ERI_F2yz_S_Gx3z_S_dd+CDZ*I_ERI_F2yz_S_Fx2z_S_dd;
  Double I_ERI_Fy2z_S_Fx2z_Pz_dd = I_ERI_Fy2z_S_Gx3z_S_dd+CDZ*I_ERI_Fy2z_S_Fx2z_S_dd;
  Double I_ERI_F3z_S_Fx2z_Pz_dd = I_ERI_F3z_S_Gx3z_S_dd+CDZ*I_ERI_F3z_S_Fx2z_S_dd;
  Double I_ERI_F3x_S_F2yz_Pz_dd = I_ERI_F3x_S_G2y2z_S_dd+CDZ*I_ERI_F3x_S_F2yz_S_dd;
  Double I_ERI_F2xy_S_F2yz_Pz_dd = I_ERI_F2xy_S_G2y2z_S_dd+CDZ*I_ERI_F2xy_S_F2yz_S_dd;
  Double I_ERI_F2xz_S_F2yz_Pz_dd = I_ERI_F2xz_S_G2y2z_S_dd+CDZ*I_ERI_F2xz_S_F2yz_S_dd;
  Double I_ERI_Fx2y_S_F2yz_Pz_dd = I_ERI_Fx2y_S_G2y2z_S_dd+CDZ*I_ERI_Fx2y_S_F2yz_S_dd;
  Double I_ERI_Fxyz_S_F2yz_Pz_dd = I_ERI_Fxyz_S_G2y2z_S_dd+CDZ*I_ERI_Fxyz_S_F2yz_S_dd;
  Double I_ERI_Fx2z_S_F2yz_Pz_dd = I_ERI_Fx2z_S_G2y2z_S_dd+CDZ*I_ERI_Fx2z_S_F2yz_S_dd;
  Double I_ERI_F3y_S_F2yz_Pz_dd = I_ERI_F3y_S_G2y2z_S_dd+CDZ*I_ERI_F3y_S_F2yz_S_dd;
  Double I_ERI_F2yz_S_F2yz_Pz_dd = I_ERI_F2yz_S_G2y2z_S_dd+CDZ*I_ERI_F2yz_S_F2yz_S_dd;
  Double I_ERI_Fy2z_S_F2yz_Pz_dd = I_ERI_Fy2z_S_G2y2z_S_dd+CDZ*I_ERI_Fy2z_S_F2yz_S_dd;
  Double I_ERI_F3z_S_F2yz_Pz_dd = I_ERI_F3z_S_G2y2z_S_dd+CDZ*I_ERI_F3z_S_F2yz_S_dd;
  Double I_ERI_F3x_S_Fy2z_Pz_dd = I_ERI_F3x_S_Gy3z_S_dd+CDZ*I_ERI_F3x_S_Fy2z_S_dd;
  Double I_ERI_F2xy_S_Fy2z_Pz_dd = I_ERI_F2xy_S_Gy3z_S_dd+CDZ*I_ERI_F2xy_S_Fy2z_S_dd;
  Double I_ERI_F2xz_S_Fy2z_Pz_dd = I_ERI_F2xz_S_Gy3z_S_dd+CDZ*I_ERI_F2xz_S_Fy2z_S_dd;
  Double I_ERI_Fx2y_S_Fy2z_Pz_dd = I_ERI_Fx2y_S_Gy3z_S_dd+CDZ*I_ERI_Fx2y_S_Fy2z_S_dd;
  Double I_ERI_Fxyz_S_Fy2z_Pz_dd = I_ERI_Fxyz_S_Gy3z_S_dd+CDZ*I_ERI_Fxyz_S_Fy2z_S_dd;
  Double I_ERI_Fx2z_S_Fy2z_Pz_dd = I_ERI_Fx2z_S_Gy3z_S_dd+CDZ*I_ERI_Fx2z_S_Fy2z_S_dd;
  Double I_ERI_F3y_S_Fy2z_Pz_dd = I_ERI_F3y_S_Gy3z_S_dd+CDZ*I_ERI_F3y_S_Fy2z_S_dd;
  Double I_ERI_F2yz_S_Fy2z_Pz_dd = I_ERI_F2yz_S_Gy3z_S_dd+CDZ*I_ERI_F2yz_S_Fy2z_S_dd;
  Double I_ERI_Fy2z_S_Fy2z_Pz_dd = I_ERI_Fy2z_S_Gy3z_S_dd+CDZ*I_ERI_Fy2z_S_Fy2z_S_dd;
  Double I_ERI_F3z_S_Fy2z_Pz_dd = I_ERI_F3z_S_Gy3z_S_dd+CDZ*I_ERI_F3z_S_Fy2z_S_dd;
  Double I_ERI_F3x_S_F3z_Pz_dd = I_ERI_F3x_S_G4z_S_dd+CDZ*I_ERI_F3x_S_F3z_S_dd;
  Double I_ERI_F2xy_S_F3z_Pz_dd = I_ERI_F2xy_S_G4z_S_dd+CDZ*I_ERI_F2xy_S_F3z_S_dd;
  Double I_ERI_F2xz_S_F3z_Pz_dd = I_ERI_F2xz_S_G4z_S_dd+CDZ*I_ERI_F2xz_S_F3z_S_dd;
  Double I_ERI_Fx2y_S_F3z_Pz_dd = I_ERI_Fx2y_S_G4z_S_dd+CDZ*I_ERI_Fx2y_S_F3z_S_dd;
  Double I_ERI_Fxyz_S_F3z_Pz_dd = I_ERI_Fxyz_S_G4z_S_dd+CDZ*I_ERI_Fxyz_S_F3z_S_dd;
  Double I_ERI_Fx2z_S_F3z_Pz_dd = I_ERI_Fx2z_S_G4z_S_dd+CDZ*I_ERI_Fx2z_S_F3z_S_dd;
  Double I_ERI_F3y_S_F3z_Pz_dd = I_ERI_F3y_S_G4z_S_dd+CDZ*I_ERI_F3y_S_F3z_S_dd;
  Double I_ERI_F2yz_S_F3z_Pz_dd = I_ERI_F2yz_S_G4z_S_dd+CDZ*I_ERI_F2yz_S_F3z_S_dd;
  Double I_ERI_Fy2z_S_F3z_Pz_dd = I_ERI_Fy2z_S_G4z_S_dd+CDZ*I_ERI_Fy2z_S_F3z_S_dd;
  Double I_ERI_F3z_S_F3z_Pz_dd = I_ERI_F3z_S_G4z_S_dd+CDZ*I_ERI_F3z_S_F3z_S_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_D_dd
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_P_dd
   ************************************************************/
  Double I_ERI_F3x_S_D2x_D2x_dd = I_ERI_F3x_S_F3x_Px_dd+CDX*I_ERI_F3x_S_D2x_Px_dd;
  Double I_ERI_F2xy_S_D2x_D2x_dd = I_ERI_F2xy_S_F3x_Px_dd+CDX*I_ERI_F2xy_S_D2x_Px_dd;
  Double I_ERI_F2xz_S_D2x_D2x_dd = I_ERI_F2xz_S_F3x_Px_dd+CDX*I_ERI_F2xz_S_D2x_Px_dd;
  Double I_ERI_Fx2y_S_D2x_D2x_dd = I_ERI_Fx2y_S_F3x_Px_dd+CDX*I_ERI_Fx2y_S_D2x_Px_dd;
  Double I_ERI_Fxyz_S_D2x_D2x_dd = I_ERI_Fxyz_S_F3x_Px_dd+CDX*I_ERI_Fxyz_S_D2x_Px_dd;
  Double I_ERI_Fx2z_S_D2x_D2x_dd = I_ERI_Fx2z_S_F3x_Px_dd+CDX*I_ERI_Fx2z_S_D2x_Px_dd;
  Double I_ERI_F3y_S_D2x_D2x_dd = I_ERI_F3y_S_F3x_Px_dd+CDX*I_ERI_F3y_S_D2x_Px_dd;
  Double I_ERI_F2yz_S_D2x_D2x_dd = I_ERI_F2yz_S_F3x_Px_dd+CDX*I_ERI_F2yz_S_D2x_Px_dd;
  Double I_ERI_Fy2z_S_D2x_D2x_dd = I_ERI_Fy2z_S_F3x_Px_dd+CDX*I_ERI_Fy2z_S_D2x_Px_dd;
  Double I_ERI_F3z_S_D2x_D2x_dd = I_ERI_F3z_S_F3x_Px_dd+CDX*I_ERI_F3z_S_D2x_Px_dd;
  Double I_ERI_F3x_S_Dxy_D2x_dd = I_ERI_F3x_S_F2xy_Px_dd+CDX*I_ERI_F3x_S_Dxy_Px_dd;
  Double I_ERI_F2xy_S_Dxy_D2x_dd = I_ERI_F2xy_S_F2xy_Px_dd+CDX*I_ERI_F2xy_S_Dxy_Px_dd;
  Double I_ERI_F2xz_S_Dxy_D2x_dd = I_ERI_F2xz_S_F2xy_Px_dd+CDX*I_ERI_F2xz_S_Dxy_Px_dd;
  Double I_ERI_Fx2y_S_Dxy_D2x_dd = I_ERI_Fx2y_S_F2xy_Px_dd+CDX*I_ERI_Fx2y_S_Dxy_Px_dd;
  Double I_ERI_Fxyz_S_Dxy_D2x_dd = I_ERI_Fxyz_S_F2xy_Px_dd+CDX*I_ERI_Fxyz_S_Dxy_Px_dd;
  Double I_ERI_Fx2z_S_Dxy_D2x_dd = I_ERI_Fx2z_S_F2xy_Px_dd+CDX*I_ERI_Fx2z_S_Dxy_Px_dd;
  Double I_ERI_F3y_S_Dxy_D2x_dd = I_ERI_F3y_S_F2xy_Px_dd+CDX*I_ERI_F3y_S_Dxy_Px_dd;
  Double I_ERI_F2yz_S_Dxy_D2x_dd = I_ERI_F2yz_S_F2xy_Px_dd+CDX*I_ERI_F2yz_S_Dxy_Px_dd;
  Double I_ERI_Fy2z_S_Dxy_D2x_dd = I_ERI_Fy2z_S_F2xy_Px_dd+CDX*I_ERI_Fy2z_S_Dxy_Px_dd;
  Double I_ERI_F3z_S_Dxy_D2x_dd = I_ERI_F3z_S_F2xy_Px_dd+CDX*I_ERI_F3z_S_Dxy_Px_dd;
  Double I_ERI_F3x_S_Dxz_D2x_dd = I_ERI_F3x_S_F2xz_Px_dd+CDX*I_ERI_F3x_S_Dxz_Px_dd;
  Double I_ERI_F2xy_S_Dxz_D2x_dd = I_ERI_F2xy_S_F2xz_Px_dd+CDX*I_ERI_F2xy_S_Dxz_Px_dd;
  Double I_ERI_F2xz_S_Dxz_D2x_dd = I_ERI_F2xz_S_F2xz_Px_dd+CDX*I_ERI_F2xz_S_Dxz_Px_dd;
  Double I_ERI_Fx2y_S_Dxz_D2x_dd = I_ERI_Fx2y_S_F2xz_Px_dd+CDX*I_ERI_Fx2y_S_Dxz_Px_dd;
  Double I_ERI_Fxyz_S_Dxz_D2x_dd = I_ERI_Fxyz_S_F2xz_Px_dd+CDX*I_ERI_Fxyz_S_Dxz_Px_dd;
  Double I_ERI_Fx2z_S_Dxz_D2x_dd = I_ERI_Fx2z_S_F2xz_Px_dd+CDX*I_ERI_Fx2z_S_Dxz_Px_dd;
  Double I_ERI_F3y_S_Dxz_D2x_dd = I_ERI_F3y_S_F2xz_Px_dd+CDX*I_ERI_F3y_S_Dxz_Px_dd;
  Double I_ERI_F2yz_S_Dxz_D2x_dd = I_ERI_F2yz_S_F2xz_Px_dd+CDX*I_ERI_F2yz_S_Dxz_Px_dd;
  Double I_ERI_Fy2z_S_Dxz_D2x_dd = I_ERI_Fy2z_S_F2xz_Px_dd+CDX*I_ERI_Fy2z_S_Dxz_Px_dd;
  Double I_ERI_F3z_S_Dxz_D2x_dd = I_ERI_F3z_S_F2xz_Px_dd+CDX*I_ERI_F3z_S_Dxz_Px_dd;
  Double I_ERI_F3x_S_D2y_D2x_dd = I_ERI_F3x_S_Fx2y_Px_dd+CDX*I_ERI_F3x_S_D2y_Px_dd;
  Double I_ERI_F2xy_S_D2y_D2x_dd = I_ERI_F2xy_S_Fx2y_Px_dd+CDX*I_ERI_F2xy_S_D2y_Px_dd;
  Double I_ERI_F2xz_S_D2y_D2x_dd = I_ERI_F2xz_S_Fx2y_Px_dd+CDX*I_ERI_F2xz_S_D2y_Px_dd;
  Double I_ERI_Fx2y_S_D2y_D2x_dd = I_ERI_Fx2y_S_Fx2y_Px_dd+CDX*I_ERI_Fx2y_S_D2y_Px_dd;
  Double I_ERI_Fxyz_S_D2y_D2x_dd = I_ERI_Fxyz_S_Fx2y_Px_dd+CDX*I_ERI_Fxyz_S_D2y_Px_dd;
  Double I_ERI_Fx2z_S_D2y_D2x_dd = I_ERI_Fx2z_S_Fx2y_Px_dd+CDX*I_ERI_Fx2z_S_D2y_Px_dd;
  Double I_ERI_F3y_S_D2y_D2x_dd = I_ERI_F3y_S_Fx2y_Px_dd+CDX*I_ERI_F3y_S_D2y_Px_dd;
  Double I_ERI_F2yz_S_D2y_D2x_dd = I_ERI_F2yz_S_Fx2y_Px_dd+CDX*I_ERI_F2yz_S_D2y_Px_dd;
  Double I_ERI_Fy2z_S_D2y_D2x_dd = I_ERI_Fy2z_S_Fx2y_Px_dd+CDX*I_ERI_Fy2z_S_D2y_Px_dd;
  Double I_ERI_F3z_S_D2y_D2x_dd = I_ERI_F3z_S_Fx2y_Px_dd+CDX*I_ERI_F3z_S_D2y_Px_dd;
  Double I_ERI_F3x_S_Dyz_D2x_dd = I_ERI_F3x_S_Fxyz_Px_dd+CDX*I_ERI_F3x_S_Dyz_Px_dd;
  Double I_ERI_F2xy_S_Dyz_D2x_dd = I_ERI_F2xy_S_Fxyz_Px_dd+CDX*I_ERI_F2xy_S_Dyz_Px_dd;
  Double I_ERI_F2xz_S_Dyz_D2x_dd = I_ERI_F2xz_S_Fxyz_Px_dd+CDX*I_ERI_F2xz_S_Dyz_Px_dd;
  Double I_ERI_Fx2y_S_Dyz_D2x_dd = I_ERI_Fx2y_S_Fxyz_Px_dd+CDX*I_ERI_Fx2y_S_Dyz_Px_dd;
  Double I_ERI_Fxyz_S_Dyz_D2x_dd = I_ERI_Fxyz_S_Fxyz_Px_dd+CDX*I_ERI_Fxyz_S_Dyz_Px_dd;
  Double I_ERI_Fx2z_S_Dyz_D2x_dd = I_ERI_Fx2z_S_Fxyz_Px_dd+CDX*I_ERI_Fx2z_S_Dyz_Px_dd;
  Double I_ERI_F3y_S_Dyz_D2x_dd = I_ERI_F3y_S_Fxyz_Px_dd+CDX*I_ERI_F3y_S_Dyz_Px_dd;
  Double I_ERI_F2yz_S_Dyz_D2x_dd = I_ERI_F2yz_S_Fxyz_Px_dd+CDX*I_ERI_F2yz_S_Dyz_Px_dd;
  Double I_ERI_Fy2z_S_Dyz_D2x_dd = I_ERI_Fy2z_S_Fxyz_Px_dd+CDX*I_ERI_Fy2z_S_Dyz_Px_dd;
  Double I_ERI_F3z_S_Dyz_D2x_dd = I_ERI_F3z_S_Fxyz_Px_dd+CDX*I_ERI_F3z_S_Dyz_Px_dd;
  Double I_ERI_F3x_S_D2z_D2x_dd = I_ERI_F3x_S_Fx2z_Px_dd+CDX*I_ERI_F3x_S_D2z_Px_dd;
  Double I_ERI_F2xy_S_D2z_D2x_dd = I_ERI_F2xy_S_Fx2z_Px_dd+CDX*I_ERI_F2xy_S_D2z_Px_dd;
  Double I_ERI_F2xz_S_D2z_D2x_dd = I_ERI_F2xz_S_Fx2z_Px_dd+CDX*I_ERI_F2xz_S_D2z_Px_dd;
  Double I_ERI_Fx2y_S_D2z_D2x_dd = I_ERI_Fx2y_S_Fx2z_Px_dd+CDX*I_ERI_Fx2y_S_D2z_Px_dd;
  Double I_ERI_Fxyz_S_D2z_D2x_dd = I_ERI_Fxyz_S_Fx2z_Px_dd+CDX*I_ERI_Fxyz_S_D2z_Px_dd;
  Double I_ERI_Fx2z_S_D2z_D2x_dd = I_ERI_Fx2z_S_Fx2z_Px_dd+CDX*I_ERI_Fx2z_S_D2z_Px_dd;
  Double I_ERI_F3y_S_D2z_D2x_dd = I_ERI_F3y_S_Fx2z_Px_dd+CDX*I_ERI_F3y_S_D2z_Px_dd;
  Double I_ERI_F2yz_S_D2z_D2x_dd = I_ERI_F2yz_S_Fx2z_Px_dd+CDX*I_ERI_F2yz_S_D2z_Px_dd;
  Double I_ERI_Fy2z_S_D2z_D2x_dd = I_ERI_Fy2z_S_Fx2z_Px_dd+CDX*I_ERI_Fy2z_S_D2z_Px_dd;
  Double I_ERI_F3z_S_D2z_D2x_dd = I_ERI_F3z_S_Fx2z_Px_dd+CDX*I_ERI_F3z_S_D2z_Px_dd;
  Double I_ERI_F3x_S_D2x_Dxy_dd = I_ERI_F3x_S_F2xy_Px_dd+CDY*I_ERI_F3x_S_D2x_Px_dd;
  Double I_ERI_F2xy_S_D2x_Dxy_dd = I_ERI_F2xy_S_F2xy_Px_dd+CDY*I_ERI_F2xy_S_D2x_Px_dd;
  Double I_ERI_F2xz_S_D2x_Dxy_dd = I_ERI_F2xz_S_F2xy_Px_dd+CDY*I_ERI_F2xz_S_D2x_Px_dd;
  Double I_ERI_Fx2y_S_D2x_Dxy_dd = I_ERI_Fx2y_S_F2xy_Px_dd+CDY*I_ERI_Fx2y_S_D2x_Px_dd;
  Double I_ERI_Fxyz_S_D2x_Dxy_dd = I_ERI_Fxyz_S_F2xy_Px_dd+CDY*I_ERI_Fxyz_S_D2x_Px_dd;
  Double I_ERI_Fx2z_S_D2x_Dxy_dd = I_ERI_Fx2z_S_F2xy_Px_dd+CDY*I_ERI_Fx2z_S_D2x_Px_dd;
  Double I_ERI_F3y_S_D2x_Dxy_dd = I_ERI_F3y_S_F2xy_Px_dd+CDY*I_ERI_F3y_S_D2x_Px_dd;
  Double I_ERI_F2yz_S_D2x_Dxy_dd = I_ERI_F2yz_S_F2xy_Px_dd+CDY*I_ERI_F2yz_S_D2x_Px_dd;
  Double I_ERI_Fy2z_S_D2x_Dxy_dd = I_ERI_Fy2z_S_F2xy_Px_dd+CDY*I_ERI_Fy2z_S_D2x_Px_dd;
  Double I_ERI_F3z_S_D2x_Dxy_dd = I_ERI_F3z_S_F2xy_Px_dd+CDY*I_ERI_F3z_S_D2x_Px_dd;
  Double I_ERI_F3x_S_Dxy_Dxy_dd = I_ERI_F3x_S_Fx2y_Px_dd+CDY*I_ERI_F3x_S_Dxy_Px_dd;
  Double I_ERI_F2xy_S_Dxy_Dxy_dd = I_ERI_F2xy_S_Fx2y_Px_dd+CDY*I_ERI_F2xy_S_Dxy_Px_dd;
  Double I_ERI_F2xz_S_Dxy_Dxy_dd = I_ERI_F2xz_S_Fx2y_Px_dd+CDY*I_ERI_F2xz_S_Dxy_Px_dd;
  Double I_ERI_Fx2y_S_Dxy_Dxy_dd = I_ERI_Fx2y_S_Fx2y_Px_dd+CDY*I_ERI_Fx2y_S_Dxy_Px_dd;
  Double I_ERI_Fxyz_S_Dxy_Dxy_dd = I_ERI_Fxyz_S_Fx2y_Px_dd+CDY*I_ERI_Fxyz_S_Dxy_Px_dd;
  Double I_ERI_Fx2z_S_Dxy_Dxy_dd = I_ERI_Fx2z_S_Fx2y_Px_dd+CDY*I_ERI_Fx2z_S_Dxy_Px_dd;
  Double I_ERI_F3y_S_Dxy_Dxy_dd = I_ERI_F3y_S_Fx2y_Px_dd+CDY*I_ERI_F3y_S_Dxy_Px_dd;
  Double I_ERI_F2yz_S_Dxy_Dxy_dd = I_ERI_F2yz_S_Fx2y_Px_dd+CDY*I_ERI_F2yz_S_Dxy_Px_dd;
  Double I_ERI_Fy2z_S_Dxy_Dxy_dd = I_ERI_Fy2z_S_Fx2y_Px_dd+CDY*I_ERI_Fy2z_S_Dxy_Px_dd;
  Double I_ERI_F3z_S_Dxy_Dxy_dd = I_ERI_F3z_S_Fx2y_Px_dd+CDY*I_ERI_F3z_S_Dxy_Px_dd;
  Double I_ERI_F3x_S_Dxz_Dxy_dd = I_ERI_F3x_S_Fxyz_Px_dd+CDY*I_ERI_F3x_S_Dxz_Px_dd;
  Double I_ERI_F2xy_S_Dxz_Dxy_dd = I_ERI_F2xy_S_Fxyz_Px_dd+CDY*I_ERI_F2xy_S_Dxz_Px_dd;
  Double I_ERI_F2xz_S_Dxz_Dxy_dd = I_ERI_F2xz_S_Fxyz_Px_dd+CDY*I_ERI_F2xz_S_Dxz_Px_dd;
  Double I_ERI_Fx2y_S_Dxz_Dxy_dd = I_ERI_Fx2y_S_Fxyz_Px_dd+CDY*I_ERI_Fx2y_S_Dxz_Px_dd;
  Double I_ERI_Fxyz_S_Dxz_Dxy_dd = I_ERI_Fxyz_S_Fxyz_Px_dd+CDY*I_ERI_Fxyz_S_Dxz_Px_dd;
  Double I_ERI_Fx2z_S_Dxz_Dxy_dd = I_ERI_Fx2z_S_Fxyz_Px_dd+CDY*I_ERI_Fx2z_S_Dxz_Px_dd;
  Double I_ERI_F3y_S_Dxz_Dxy_dd = I_ERI_F3y_S_Fxyz_Px_dd+CDY*I_ERI_F3y_S_Dxz_Px_dd;
  Double I_ERI_F2yz_S_Dxz_Dxy_dd = I_ERI_F2yz_S_Fxyz_Px_dd+CDY*I_ERI_F2yz_S_Dxz_Px_dd;
  Double I_ERI_Fy2z_S_Dxz_Dxy_dd = I_ERI_Fy2z_S_Fxyz_Px_dd+CDY*I_ERI_Fy2z_S_Dxz_Px_dd;
  Double I_ERI_F3z_S_Dxz_Dxy_dd = I_ERI_F3z_S_Fxyz_Px_dd+CDY*I_ERI_F3z_S_Dxz_Px_dd;
  Double I_ERI_F3x_S_D2y_Dxy_dd = I_ERI_F3x_S_F3y_Px_dd+CDY*I_ERI_F3x_S_D2y_Px_dd;
  Double I_ERI_F2xy_S_D2y_Dxy_dd = I_ERI_F2xy_S_F3y_Px_dd+CDY*I_ERI_F2xy_S_D2y_Px_dd;
  Double I_ERI_F2xz_S_D2y_Dxy_dd = I_ERI_F2xz_S_F3y_Px_dd+CDY*I_ERI_F2xz_S_D2y_Px_dd;
  Double I_ERI_Fx2y_S_D2y_Dxy_dd = I_ERI_Fx2y_S_F3y_Px_dd+CDY*I_ERI_Fx2y_S_D2y_Px_dd;
  Double I_ERI_Fxyz_S_D2y_Dxy_dd = I_ERI_Fxyz_S_F3y_Px_dd+CDY*I_ERI_Fxyz_S_D2y_Px_dd;
  Double I_ERI_Fx2z_S_D2y_Dxy_dd = I_ERI_Fx2z_S_F3y_Px_dd+CDY*I_ERI_Fx2z_S_D2y_Px_dd;
  Double I_ERI_F3y_S_D2y_Dxy_dd = I_ERI_F3y_S_F3y_Px_dd+CDY*I_ERI_F3y_S_D2y_Px_dd;
  Double I_ERI_F2yz_S_D2y_Dxy_dd = I_ERI_F2yz_S_F3y_Px_dd+CDY*I_ERI_F2yz_S_D2y_Px_dd;
  Double I_ERI_Fy2z_S_D2y_Dxy_dd = I_ERI_Fy2z_S_F3y_Px_dd+CDY*I_ERI_Fy2z_S_D2y_Px_dd;
  Double I_ERI_F3z_S_D2y_Dxy_dd = I_ERI_F3z_S_F3y_Px_dd+CDY*I_ERI_F3z_S_D2y_Px_dd;
  Double I_ERI_F3x_S_Dyz_Dxy_dd = I_ERI_F3x_S_F2yz_Px_dd+CDY*I_ERI_F3x_S_Dyz_Px_dd;
  Double I_ERI_F2xy_S_Dyz_Dxy_dd = I_ERI_F2xy_S_F2yz_Px_dd+CDY*I_ERI_F2xy_S_Dyz_Px_dd;
  Double I_ERI_F2xz_S_Dyz_Dxy_dd = I_ERI_F2xz_S_F2yz_Px_dd+CDY*I_ERI_F2xz_S_Dyz_Px_dd;
  Double I_ERI_Fx2y_S_Dyz_Dxy_dd = I_ERI_Fx2y_S_F2yz_Px_dd+CDY*I_ERI_Fx2y_S_Dyz_Px_dd;
  Double I_ERI_Fxyz_S_Dyz_Dxy_dd = I_ERI_Fxyz_S_F2yz_Px_dd+CDY*I_ERI_Fxyz_S_Dyz_Px_dd;
  Double I_ERI_Fx2z_S_Dyz_Dxy_dd = I_ERI_Fx2z_S_F2yz_Px_dd+CDY*I_ERI_Fx2z_S_Dyz_Px_dd;
  Double I_ERI_F3y_S_Dyz_Dxy_dd = I_ERI_F3y_S_F2yz_Px_dd+CDY*I_ERI_F3y_S_Dyz_Px_dd;
  Double I_ERI_F2yz_S_Dyz_Dxy_dd = I_ERI_F2yz_S_F2yz_Px_dd+CDY*I_ERI_F2yz_S_Dyz_Px_dd;
  Double I_ERI_Fy2z_S_Dyz_Dxy_dd = I_ERI_Fy2z_S_F2yz_Px_dd+CDY*I_ERI_Fy2z_S_Dyz_Px_dd;
  Double I_ERI_F3z_S_Dyz_Dxy_dd = I_ERI_F3z_S_F2yz_Px_dd+CDY*I_ERI_F3z_S_Dyz_Px_dd;
  Double I_ERI_F3x_S_D2z_Dxy_dd = I_ERI_F3x_S_Fy2z_Px_dd+CDY*I_ERI_F3x_S_D2z_Px_dd;
  Double I_ERI_F2xy_S_D2z_Dxy_dd = I_ERI_F2xy_S_Fy2z_Px_dd+CDY*I_ERI_F2xy_S_D2z_Px_dd;
  Double I_ERI_F2xz_S_D2z_Dxy_dd = I_ERI_F2xz_S_Fy2z_Px_dd+CDY*I_ERI_F2xz_S_D2z_Px_dd;
  Double I_ERI_Fx2y_S_D2z_Dxy_dd = I_ERI_Fx2y_S_Fy2z_Px_dd+CDY*I_ERI_Fx2y_S_D2z_Px_dd;
  Double I_ERI_Fxyz_S_D2z_Dxy_dd = I_ERI_Fxyz_S_Fy2z_Px_dd+CDY*I_ERI_Fxyz_S_D2z_Px_dd;
  Double I_ERI_Fx2z_S_D2z_Dxy_dd = I_ERI_Fx2z_S_Fy2z_Px_dd+CDY*I_ERI_Fx2z_S_D2z_Px_dd;
  Double I_ERI_F3y_S_D2z_Dxy_dd = I_ERI_F3y_S_Fy2z_Px_dd+CDY*I_ERI_F3y_S_D2z_Px_dd;
  Double I_ERI_F2yz_S_D2z_Dxy_dd = I_ERI_F2yz_S_Fy2z_Px_dd+CDY*I_ERI_F2yz_S_D2z_Px_dd;
  Double I_ERI_Fy2z_S_D2z_Dxy_dd = I_ERI_Fy2z_S_Fy2z_Px_dd+CDY*I_ERI_Fy2z_S_D2z_Px_dd;
  Double I_ERI_F3z_S_D2z_Dxy_dd = I_ERI_F3z_S_Fy2z_Px_dd+CDY*I_ERI_F3z_S_D2z_Px_dd;
  Double I_ERI_F3x_S_D2x_Dxz_dd = I_ERI_F3x_S_F2xz_Px_dd+CDZ*I_ERI_F3x_S_D2x_Px_dd;
  Double I_ERI_F2xy_S_D2x_Dxz_dd = I_ERI_F2xy_S_F2xz_Px_dd+CDZ*I_ERI_F2xy_S_D2x_Px_dd;
  Double I_ERI_F2xz_S_D2x_Dxz_dd = I_ERI_F2xz_S_F2xz_Px_dd+CDZ*I_ERI_F2xz_S_D2x_Px_dd;
  Double I_ERI_Fx2y_S_D2x_Dxz_dd = I_ERI_Fx2y_S_F2xz_Px_dd+CDZ*I_ERI_Fx2y_S_D2x_Px_dd;
  Double I_ERI_Fxyz_S_D2x_Dxz_dd = I_ERI_Fxyz_S_F2xz_Px_dd+CDZ*I_ERI_Fxyz_S_D2x_Px_dd;
  Double I_ERI_Fx2z_S_D2x_Dxz_dd = I_ERI_Fx2z_S_F2xz_Px_dd+CDZ*I_ERI_Fx2z_S_D2x_Px_dd;
  Double I_ERI_F3y_S_D2x_Dxz_dd = I_ERI_F3y_S_F2xz_Px_dd+CDZ*I_ERI_F3y_S_D2x_Px_dd;
  Double I_ERI_F2yz_S_D2x_Dxz_dd = I_ERI_F2yz_S_F2xz_Px_dd+CDZ*I_ERI_F2yz_S_D2x_Px_dd;
  Double I_ERI_Fy2z_S_D2x_Dxz_dd = I_ERI_Fy2z_S_F2xz_Px_dd+CDZ*I_ERI_Fy2z_S_D2x_Px_dd;
  Double I_ERI_F3z_S_D2x_Dxz_dd = I_ERI_F3z_S_F2xz_Px_dd+CDZ*I_ERI_F3z_S_D2x_Px_dd;
  Double I_ERI_F3x_S_Dxy_Dxz_dd = I_ERI_F3x_S_Fxyz_Px_dd+CDZ*I_ERI_F3x_S_Dxy_Px_dd;
  Double I_ERI_F2xy_S_Dxy_Dxz_dd = I_ERI_F2xy_S_Fxyz_Px_dd+CDZ*I_ERI_F2xy_S_Dxy_Px_dd;
  Double I_ERI_F2xz_S_Dxy_Dxz_dd = I_ERI_F2xz_S_Fxyz_Px_dd+CDZ*I_ERI_F2xz_S_Dxy_Px_dd;
  Double I_ERI_Fx2y_S_Dxy_Dxz_dd = I_ERI_Fx2y_S_Fxyz_Px_dd+CDZ*I_ERI_Fx2y_S_Dxy_Px_dd;
  Double I_ERI_Fxyz_S_Dxy_Dxz_dd = I_ERI_Fxyz_S_Fxyz_Px_dd+CDZ*I_ERI_Fxyz_S_Dxy_Px_dd;
  Double I_ERI_Fx2z_S_Dxy_Dxz_dd = I_ERI_Fx2z_S_Fxyz_Px_dd+CDZ*I_ERI_Fx2z_S_Dxy_Px_dd;
  Double I_ERI_F3y_S_Dxy_Dxz_dd = I_ERI_F3y_S_Fxyz_Px_dd+CDZ*I_ERI_F3y_S_Dxy_Px_dd;
  Double I_ERI_F2yz_S_Dxy_Dxz_dd = I_ERI_F2yz_S_Fxyz_Px_dd+CDZ*I_ERI_F2yz_S_Dxy_Px_dd;
  Double I_ERI_Fy2z_S_Dxy_Dxz_dd = I_ERI_Fy2z_S_Fxyz_Px_dd+CDZ*I_ERI_Fy2z_S_Dxy_Px_dd;
  Double I_ERI_F3z_S_Dxy_Dxz_dd = I_ERI_F3z_S_Fxyz_Px_dd+CDZ*I_ERI_F3z_S_Dxy_Px_dd;
  Double I_ERI_F3x_S_Dxz_Dxz_dd = I_ERI_F3x_S_Fx2z_Px_dd+CDZ*I_ERI_F3x_S_Dxz_Px_dd;
  Double I_ERI_F2xy_S_Dxz_Dxz_dd = I_ERI_F2xy_S_Fx2z_Px_dd+CDZ*I_ERI_F2xy_S_Dxz_Px_dd;
  Double I_ERI_F2xz_S_Dxz_Dxz_dd = I_ERI_F2xz_S_Fx2z_Px_dd+CDZ*I_ERI_F2xz_S_Dxz_Px_dd;
  Double I_ERI_Fx2y_S_Dxz_Dxz_dd = I_ERI_Fx2y_S_Fx2z_Px_dd+CDZ*I_ERI_Fx2y_S_Dxz_Px_dd;
  Double I_ERI_Fxyz_S_Dxz_Dxz_dd = I_ERI_Fxyz_S_Fx2z_Px_dd+CDZ*I_ERI_Fxyz_S_Dxz_Px_dd;
  Double I_ERI_Fx2z_S_Dxz_Dxz_dd = I_ERI_Fx2z_S_Fx2z_Px_dd+CDZ*I_ERI_Fx2z_S_Dxz_Px_dd;
  Double I_ERI_F3y_S_Dxz_Dxz_dd = I_ERI_F3y_S_Fx2z_Px_dd+CDZ*I_ERI_F3y_S_Dxz_Px_dd;
  Double I_ERI_F2yz_S_Dxz_Dxz_dd = I_ERI_F2yz_S_Fx2z_Px_dd+CDZ*I_ERI_F2yz_S_Dxz_Px_dd;
  Double I_ERI_Fy2z_S_Dxz_Dxz_dd = I_ERI_Fy2z_S_Fx2z_Px_dd+CDZ*I_ERI_Fy2z_S_Dxz_Px_dd;
  Double I_ERI_F3z_S_Dxz_Dxz_dd = I_ERI_F3z_S_Fx2z_Px_dd+CDZ*I_ERI_F3z_S_Dxz_Px_dd;
  Double I_ERI_F3x_S_D2y_Dxz_dd = I_ERI_F3x_S_F2yz_Px_dd+CDZ*I_ERI_F3x_S_D2y_Px_dd;
  Double I_ERI_F2xy_S_D2y_Dxz_dd = I_ERI_F2xy_S_F2yz_Px_dd+CDZ*I_ERI_F2xy_S_D2y_Px_dd;
  Double I_ERI_F2xz_S_D2y_Dxz_dd = I_ERI_F2xz_S_F2yz_Px_dd+CDZ*I_ERI_F2xz_S_D2y_Px_dd;
  Double I_ERI_Fx2y_S_D2y_Dxz_dd = I_ERI_Fx2y_S_F2yz_Px_dd+CDZ*I_ERI_Fx2y_S_D2y_Px_dd;
  Double I_ERI_Fxyz_S_D2y_Dxz_dd = I_ERI_Fxyz_S_F2yz_Px_dd+CDZ*I_ERI_Fxyz_S_D2y_Px_dd;
  Double I_ERI_Fx2z_S_D2y_Dxz_dd = I_ERI_Fx2z_S_F2yz_Px_dd+CDZ*I_ERI_Fx2z_S_D2y_Px_dd;
  Double I_ERI_F3y_S_D2y_Dxz_dd = I_ERI_F3y_S_F2yz_Px_dd+CDZ*I_ERI_F3y_S_D2y_Px_dd;
  Double I_ERI_F2yz_S_D2y_Dxz_dd = I_ERI_F2yz_S_F2yz_Px_dd+CDZ*I_ERI_F2yz_S_D2y_Px_dd;
  Double I_ERI_Fy2z_S_D2y_Dxz_dd = I_ERI_Fy2z_S_F2yz_Px_dd+CDZ*I_ERI_Fy2z_S_D2y_Px_dd;
  Double I_ERI_F3z_S_D2y_Dxz_dd = I_ERI_F3z_S_F2yz_Px_dd+CDZ*I_ERI_F3z_S_D2y_Px_dd;
  Double I_ERI_F3x_S_Dyz_Dxz_dd = I_ERI_F3x_S_Fy2z_Px_dd+CDZ*I_ERI_F3x_S_Dyz_Px_dd;
  Double I_ERI_F2xy_S_Dyz_Dxz_dd = I_ERI_F2xy_S_Fy2z_Px_dd+CDZ*I_ERI_F2xy_S_Dyz_Px_dd;
  Double I_ERI_F2xz_S_Dyz_Dxz_dd = I_ERI_F2xz_S_Fy2z_Px_dd+CDZ*I_ERI_F2xz_S_Dyz_Px_dd;
  Double I_ERI_Fx2y_S_Dyz_Dxz_dd = I_ERI_Fx2y_S_Fy2z_Px_dd+CDZ*I_ERI_Fx2y_S_Dyz_Px_dd;
  Double I_ERI_Fxyz_S_Dyz_Dxz_dd = I_ERI_Fxyz_S_Fy2z_Px_dd+CDZ*I_ERI_Fxyz_S_Dyz_Px_dd;
  Double I_ERI_Fx2z_S_Dyz_Dxz_dd = I_ERI_Fx2z_S_Fy2z_Px_dd+CDZ*I_ERI_Fx2z_S_Dyz_Px_dd;
  Double I_ERI_F3y_S_Dyz_Dxz_dd = I_ERI_F3y_S_Fy2z_Px_dd+CDZ*I_ERI_F3y_S_Dyz_Px_dd;
  Double I_ERI_F2yz_S_Dyz_Dxz_dd = I_ERI_F2yz_S_Fy2z_Px_dd+CDZ*I_ERI_F2yz_S_Dyz_Px_dd;
  Double I_ERI_Fy2z_S_Dyz_Dxz_dd = I_ERI_Fy2z_S_Fy2z_Px_dd+CDZ*I_ERI_Fy2z_S_Dyz_Px_dd;
  Double I_ERI_F3z_S_Dyz_Dxz_dd = I_ERI_F3z_S_Fy2z_Px_dd+CDZ*I_ERI_F3z_S_Dyz_Px_dd;
  Double I_ERI_F3x_S_D2z_Dxz_dd = I_ERI_F3x_S_F3z_Px_dd+CDZ*I_ERI_F3x_S_D2z_Px_dd;
  Double I_ERI_F2xy_S_D2z_Dxz_dd = I_ERI_F2xy_S_F3z_Px_dd+CDZ*I_ERI_F2xy_S_D2z_Px_dd;
  Double I_ERI_F2xz_S_D2z_Dxz_dd = I_ERI_F2xz_S_F3z_Px_dd+CDZ*I_ERI_F2xz_S_D2z_Px_dd;
  Double I_ERI_Fx2y_S_D2z_Dxz_dd = I_ERI_Fx2y_S_F3z_Px_dd+CDZ*I_ERI_Fx2y_S_D2z_Px_dd;
  Double I_ERI_Fxyz_S_D2z_Dxz_dd = I_ERI_Fxyz_S_F3z_Px_dd+CDZ*I_ERI_Fxyz_S_D2z_Px_dd;
  Double I_ERI_Fx2z_S_D2z_Dxz_dd = I_ERI_Fx2z_S_F3z_Px_dd+CDZ*I_ERI_Fx2z_S_D2z_Px_dd;
  Double I_ERI_F3y_S_D2z_Dxz_dd = I_ERI_F3y_S_F3z_Px_dd+CDZ*I_ERI_F3y_S_D2z_Px_dd;
  Double I_ERI_F2yz_S_D2z_Dxz_dd = I_ERI_F2yz_S_F3z_Px_dd+CDZ*I_ERI_F2yz_S_D2z_Px_dd;
  Double I_ERI_Fy2z_S_D2z_Dxz_dd = I_ERI_Fy2z_S_F3z_Px_dd+CDZ*I_ERI_Fy2z_S_D2z_Px_dd;
  Double I_ERI_F3z_S_D2z_Dxz_dd = I_ERI_F3z_S_F3z_Px_dd+CDZ*I_ERI_F3z_S_D2z_Px_dd;
  Double I_ERI_F3x_S_D2x_D2y_dd = I_ERI_F3x_S_F2xy_Py_dd+CDY*I_ERI_F3x_S_D2x_Py_dd;
  Double I_ERI_F2xy_S_D2x_D2y_dd = I_ERI_F2xy_S_F2xy_Py_dd+CDY*I_ERI_F2xy_S_D2x_Py_dd;
  Double I_ERI_F2xz_S_D2x_D2y_dd = I_ERI_F2xz_S_F2xy_Py_dd+CDY*I_ERI_F2xz_S_D2x_Py_dd;
  Double I_ERI_Fx2y_S_D2x_D2y_dd = I_ERI_Fx2y_S_F2xy_Py_dd+CDY*I_ERI_Fx2y_S_D2x_Py_dd;
  Double I_ERI_Fxyz_S_D2x_D2y_dd = I_ERI_Fxyz_S_F2xy_Py_dd+CDY*I_ERI_Fxyz_S_D2x_Py_dd;
  Double I_ERI_Fx2z_S_D2x_D2y_dd = I_ERI_Fx2z_S_F2xy_Py_dd+CDY*I_ERI_Fx2z_S_D2x_Py_dd;
  Double I_ERI_F3y_S_D2x_D2y_dd = I_ERI_F3y_S_F2xy_Py_dd+CDY*I_ERI_F3y_S_D2x_Py_dd;
  Double I_ERI_F2yz_S_D2x_D2y_dd = I_ERI_F2yz_S_F2xy_Py_dd+CDY*I_ERI_F2yz_S_D2x_Py_dd;
  Double I_ERI_Fy2z_S_D2x_D2y_dd = I_ERI_Fy2z_S_F2xy_Py_dd+CDY*I_ERI_Fy2z_S_D2x_Py_dd;
  Double I_ERI_F3z_S_D2x_D2y_dd = I_ERI_F3z_S_F2xy_Py_dd+CDY*I_ERI_F3z_S_D2x_Py_dd;
  Double I_ERI_F3x_S_Dxy_D2y_dd = I_ERI_F3x_S_Fx2y_Py_dd+CDY*I_ERI_F3x_S_Dxy_Py_dd;
  Double I_ERI_F2xy_S_Dxy_D2y_dd = I_ERI_F2xy_S_Fx2y_Py_dd+CDY*I_ERI_F2xy_S_Dxy_Py_dd;
  Double I_ERI_F2xz_S_Dxy_D2y_dd = I_ERI_F2xz_S_Fx2y_Py_dd+CDY*I_ERI_F2xz_S_Dxy_Py_dd;
  Double I_ERI_Fx2y_S_Dxy_D2y_dd = I_ERI_Fx2y_S_Fx2y_Py_dd+CDY*I_ERI_Fx2y_S_Dxy_Py_dd;
  Double I_ERI_Fxyz_S_Dxy_D2y_dd = I_ERI_Fxyz_S_Fx2y_Py_dd+CDY*I_ERI_Fxyz_S_Dxy_Py_dd;
  Double I_ERI_Fx2z_S_Dxy_D2y_dd = I_ERI_Fx2z_S_Fx2y_Py_dd+CDY*I_ERI_Fx2z_S_Dxy_Py_dd;
  Double I_ERI_F3y_S_Dxy_D2y_dd = I_ERI_F3y_S_Fx2y_Py_dd+CDY*I_ERI_F3y_S_Dxy_Py_dd;
  Double I_ERI_F2yz_S_Dxy_D2y_dd = I_ERI_F2yz_S_Fx2y_Py_dd+CDY*I_ERI_F2yz_S_Dxy_Py_dd;
  Double I_ERI_Fy2z_S_Dxy_D2y_dd = I_ERI_Fy2z_S_Fx2y_Py_dd+CDY*I_ERI_Fy2z_S_Dxy_Py_dd;
  Double I_ERI_F3z_S_Dxy_D2y_dd = I_ERI_F3z_S_Fx2y_Py_dd+CDY*I_ERI_F3z_S_Dxy_Py_dd;
  Double I_ERI_F3x_S_Dxz_D2y_dd = I_ERI_F3x_S_Fxyz_Py_dd+CDY*I_ERI_F3x_S_Dxz_Py_dd;
  Double I_ERI_F2xy_S_Dxz_D2y_dd = I_ERI_F2xy_S_Fxyz_Py_dd+CDY*I_ERI_F2xy_S_Dxz_Py_dd;
  Double I_ERI_F2xz_S_Dxz_D2y_dd = I_ERI_F2xz_S_Fxyz_Py_dd+CDY*I_ERI_F2xz_S_Dxz_Py_dd;
  Double I_ERI_Fx2y_S_Dxz_D2y_dd = I_ERI_Fx2y_S_Fxyz_Py_dd+CDY*I_ERI_Fx2y_S_Dxz_Py_dd;
  Double I_ERI_Fxyz_S_Dxz_D2y_dd = I_ERI_Fxyz_S_Fxyz_Py_dd+CDY*I_ERI_Fxyz_S_Dxz_Py_dd;
  Double I_ERI_Fx2z_S_Dxz_D2y_dd = I_ERI_Fx2z_S_Fxyz_Py_dd+CDY*I_ERI_Fx2z_S_Dxz_Py_dd;
  Double I_ERI_F3y_S_Dxz_D2y_dd = I_ERI_F3y_S_Fxyz_Py_dd+CDY*I_ERI_F3y_S_Dxz_Py_dd;
  Double I_ERI_F2yz_S_Dxz_D2y_dd = I_ERI_F2yz_S_Fxyz_Py_dd+CDY*I_ERI_F2yz_S_Dxz_Py_dd;
  Double I_ERI_Fy2z_S_Dxz_D2y_dd = I_ERI_Fy2z_S_Fxyz_Py_dd+CDY*I_ERI_Fy2z_S_Dxz_Py_dd;
  Double I_ERI_F3z_S_Dxz_D2y_dd = I_ERI_F3z_S_Fxyz_Py_dd+CDY*I_ERI_F3z_S_Dxz_Py_dd;
  Double I_ERI_F3x_S_D2y_D2y_dd = I_ERI_F3x_S_F3y_Py_dd+CDY*I_ERI_F3x_S_D2y_Py_dd;
  Double I_ERI_F2xy_S_D2y_D2y_dd = I_ERI_F2xy_S_F3y_Py_dd+CDY*I_ERI_F2xy_S_D2y_Py_dd;
  Double I_ERI_F2xz_S_D2y_D2y_dd = I_ERI_F2xz_S_F3y_Py_dd+CDY*I_ERI_F2xz_S_D2y_Py_dd;
  Double I_ERI_Fx2y_S_D2y_D2y_dd = I_ERI_Fx2y_S_F3y_Py_dd+CDY*I_ERI_Fx2y_S_D2y_Py_dd;
  Double I_ERI_Fxyz_S_D2y_D2y_dd = I_ERI_Fxyz_S_F3y_Py_dd+CDY*I_ERI_Fxyz_S_D2y_Py_dd;
  Double I_ERI_Fx2z_S_D2y_D2y_dd = I_ERI_Fx2z_S_F3y_Py_dd+CDY*I_ERI_Fx2z_S_D2y_Py_dd;
  Double I_ERI_F3y_S_D2y_D2y_dd = I_ERI_F3y_S_F3y_Py_dd+CDY*I_ERI_F3y_S_D2y_Py_dd;
  Double I_ERI_F2yz_S_D2y_D2y_dd = I_ERI_F2yz_S_F3y_Py_dd+CDY*I_ERI_F2yz_S_D2y_Py_dd;
  Double I_ERI_Fy2z_S_D2y_D2y_dd = I_ERI_Fy2z_S_F3y_Py_dd+CDY*I_ERI_Fy2z_S_D2y_Py_dd;
  Double I_ERI_F3z_S_D2y_D2y_dd = I_ERI_F3z_S_F3y_Py_dd+CDY*I_ERI_F3z_S_D2y_Py_dd;
  Double I_ERI_F3x_S_Dyz_D2y_dd = I_ERI_F3x_S_F2yz_Py_dd+CDY*I_ERI_F3x_S_Dyz_Py_dd;
  Double I_ERI_F2xy_S_Dyz_D2y_dd = I_ERI_F2xy_S_F2yz_Py_dd+CDY*I_ERI_F2xy_S_Dyz_Py_dd;
  Double I_ERI_F2xz_S_Dyz_D2y_dd = I_ERI_F2xz_S_F2yz_Py_dd+CDY*I_ERI_F2xz_S_Dyz_Py_dd;
  Double I_ERI_Fx2y_S_Dyz_D2y_dd = I_ERI_Fx2y_S_F2yz_Py_dd+CDY*I_ERI_Fx2y_S_Dyz_Py_dd;
  Double I_ERI_Fxyz_S_Dyz_D2y_dd = I_ERI_Fxyz_S_F2yz_Py_dd+CDY*I_ERI_Fxyz_S_Dyz_Py_dd;
  Double I_ERI_Fx2z_S_Dyz_D2y_dd = I_ERI_Fx2z_S_F2yz_Py_dd+CDY*I_ERI_Fx2z_S_Dyz_Py_dd;
  Double I_ERI_F3y_S_Dyz_D2y_dd = I_ERI_F3y_S_F2yz_Py_dd+CDY*I_ERI_F3y_S_Dyz_Py_dd;
  Double I_ERI_F2yz_S_Dyz_D2y_dd = I_ERI_F2yz_S_F2yz_Py_dd+CDY*I_ERI_F2yz_S_Dyz_Py_dd;
  Double I_ERI_Fy2z_S_Dyz_D2y_dd = I_ERI_Fy2z_S_F2yz_Py_dd+CDY*I_ERI_Fy2z_S_Dyz_Py_dd;
  Double I_ERI_F3z_S_Dyz_D2y_dd = I_ERI_F3z_S_F2yz_Py_dd+CDY*I_ERI_F3z_S_Dyz_Py_dd;
  Double I_ERI_F3x_S_D2z_D2y_dd = I_ERI_F3x_S_Fy2z_Py_dd+CDY*I_ERI_F3x_S_D2z_Py_dd;
  Double I_ERI_F2xy_S_D2z_D2y_dd = I_ERI_F2xy_S_Fy2z_Py_dd+CDY*I_ERI_F2xy_S_D2z_Py_dd;
  Double I_ERI_F2xz_S_D2z_D2y_dd = I_ERI_F2xz_S_Fy2z_Py_dd+CDY*I_ERI_F2xz_S_D2z_Py_dd;
  Double I_ERI_Fx2y_S_D2z_D2y_dd = I_ERI_Fx2y_S_Fy2z_Py_dd+CDY*I_ERI_Fx2y_S_D2z_Py_dd;
  Double I_ERI_Fxyz_S_D2z_D2y_dd = I_ERI_Fxyz_S_Fy2z_Py_dd+CDY*I_ERI_Fxyz_S_D2z_Py_dd;
  Double I_ERI_Fx2z_S_D2z_D2y_dd = I_ERI_Fx2z_S_Fy2z_Py_dd+CDY*I_ERI_Fx2z_S_D2z_Py_dd;
  Double I_ERI_F3y_S_D2z_D2y_dd = I_ERI_F3y_S_Fy2z_Py_dd+CDY*I_ERI_F3y_S_D2z_Py_dd;
  Double I_ERI_F2yz_S_D2z_D2y_dd = I_ERI_F2yz_S_Fy2z_Py_dd+CDY*I_ERI_F2yz_S_D2z_Py_dd;
  Double I_ERI_Fy2z_S_D2z_D2y_dd = I_ERI_Fy2z_S_Fy2z_Py_dd+CDY*I_ERI_Fy2z_S_D2z_Py_dd;
  Double I_ERI_F3z_S_D2z_D2y_dd = I_ERI_F3z_S_Fy2z_Py_dd+CDY*I_ERI_F3z_S_D2z_Py_dd;
  Double I_ERI_F3x_S_D2x_Dyz_dd = I_ERI_F3x_S_F2xz_Py_dd+CDZ*I_ERI_F3x_S_D2x_Py_dd;
  Double I_ERI_F2xy_S_D2x_Dyz_dd = I_ERI_F2xy_S_F2xz_Py_dd+CDZ*I_ERI_F2xy_S_D2x_Py_dd;
  Double I_ERI_F2xz_S_D2x_Dyz_dd = I_ERI_F2xz_S_F2xz_Py_dd+CDZ*I_ERI_F2xz_S_D2x_Py_dd;
  Double I_ERI_Fx2y_S_D2x_Dyz_dd = I_ERI_Fx2y_S_F2xz_Py_dd+CDZ*I_ERI_Fx2y_S_D2x_Py_dd;
  Double I_ERI_Fxyz_S_D2x_Dyz_dd = I_ERI_Fxyz_S_F2xz_Py_dd+CDZ*I_ERI_Fxyz_S_D2x_Py_dd;
  Double I_ERI_Fx2z_S_D2x_Dyz_dd = I_ERI_Fx2z_S_F2xz_Py_dd+CDZ*I_ERI_Fx2z_S_D2x_Py_dd;
  Double I_ERI_F3y_S_D2x_Dyz_dd = I_ERI_F3y_S_F2xz_Py_dd+CDZ*I_ERI_F3y_S_D2x_Py_dd;
  Double I_ERI_F2yz_S_D2x_Dyz_dd = I_ERI_F2yz_S_F2xz_Py_dd+CDZ*I_ERI_F2yz_S_D2x_Py_dd;
  Double I_ERI_Fy2z_S_D2x_Dyz_dd = I_ERI_Fy2z_S_F2xz_Py_dd+CDZ*I_ERI_Fy2z_S_D2x_Py_dd;
  Double I_ERI_F3z_S_D2x_Dyz_dd = I_ERI_F3z_S_F2xz_Py_dd+CDZ*I_ERI_F3z_S_D2x_Py_dd;
  Double I_ERI_F3x_S_Dxy_Dyz_dd = I_ERI_F3x_S_Fxyz_Py_dd+CDZ*I_ERI_F3x_S_Dxy_Py_dd;
  Double I_ERI_F2xy_S_Dxy_Dyz_dd = I_ERI_F2xy_S_Fxyz_Py_dd+CDZ*I_ERI_F2xy_S_Dxy_Py_dd;
  Double I_ERI_F2xz_S_Dxy_Dyz_dd = I_ERI_F2xz_S_Fxyz_Py_dd+CDZ*I_ERI_F2xz_S_Dxy_Py_dd;
  Double I_ERI_Fx2y_S_Dxy_Dyz_dd = I_ERI_Fx2y_S_Fxyz_Py_dd+CDZ*I_ERI_Fx2y_S_Dxy_Py_dd;
  Double I_ERI_Fxyz_S_Dxy_Dyz_dd = I_ERI_Fxyz_S_Fxyz_Py_dd+CDZ*I_ERI_Fxyz_S_Dxy_Py_dd;
  Double I_ERI_Fx2z_S_Dxy_Dyz_dd = I_ERI_Fx2z_S_Fxyz_Py_dd+CDZ*I_ERI_Fx2z_S_Dxy_Py_dd;
  Double I_ERI_F3y_S_Dxy_Dyz_dd = I_ERI_F3y_S_Fxyz_Py_dd+CDZ*I_ERI_F3y_S_Dxy_Py_dd;
  Double I_ERI_F2yz_S_Dxy_Dyz_dd = I_ERI_F2yz_S_Fxyz_Py_dd+CDZ*I_ERI_F2yz_S_Dxy_Py_dd;
  Double I_ERI_Fy2z_S_Dxy_Dyz_dd = I_ERI_Fy2z_S_Fxyz_Py_dd+CDZ*I_ERI_Fy2z_S_Dxy_Py_dd;
  Double I_ERI_F3z_S_Dxy_Dyz_dd = I_ERI_F3z_S_Fxyz_Py_dd+CDZ*I_ERI_F3z_S_Dxy_Py_dd;
  Double I_ERI_F3x_S_Dxz_Dyz_dd = I_ERI_F3x_S_Fx2z_Py_dd+CDZ*I_ERI_F3x_S_Dxz_Py_dd;
  Double I_ERI_F2xy_S_Dxz_Dyz_dd = I_ERI_F2xy_S_Fx2z_Py_dd+CDZ*I_ERI_F2xy_S_Dxz_Py_dd;
  Double I_ERI_F2xz_S_Dxz_Dyz_dd = I_ERI_F2xz_S_Fx2z_Py_dd+CDZ*I_ERI_F2xz_S_Dxz_Py_dd;
  Double I_ERI_Fx2y_S_Dxz_Dyz_dd = I_ERI_Fx2y_S_Fx2z_Py_dd+CDZ*I_ERI_Fx2y_S_Dxz_Py_dd;
  Double I_ERI_Fxyz_S_Dxz_Dyz_dd = I_ERI_Fxyz_S_Fx2z_Py_dd+CDZ*I_ERI_Fxyz_S_Dxz_Py_dd;
  Double I_ERI_Fx2z_S_Dxz_Dyz_dd = I_ERI_Fx2z_S_Fx2z_Py_dd+CDZ*I_ERI_Fx2z_S_Dxz_Py_dd;
  Double I_ERI_F3y_S_Dxz_Dyz_dd = I_ERI_F3y_S_Fx2z_Py_dd+CDZ*I_ERI_F3y_S_Dxz_Py_dd;
  Double I_ERI_F2yz_S_Dxz_Dyz_dd = I_ERI_F2yz_S_Fx2z_Py_dd+CDZ*I_ERI_F2yz_S_Dxz_Py_dd;
  Double I_ERI_Fy2z_S_Dxz_Dyz_dd = I_ERI_Fy2z_S_Fx2z_Py_dd+CDZ*I_ERI_Fy2z_S_Dxz_Py_dd;
  Double I_ERI_F3z_S_Dxz_Dyz_dd = I_ERI_F3z_S_Fx2z_Py_dd+CDZ*I_ERI_F3z_S_Dxz_Py_dd;
  Double I_ERI_F3x_S_D2y_Dyz_dd = I_ERI_F3x_S_F2yz_Py_dd+CDZ*I_ERI_F3x_S_D2y_Py_dd;
  Double I_ERI_F2xy_S_D2y_Dyz_dd = I_ERI_F2xy_S_F2yz_Py_dd+CDZ*I_ERI_F2xy_S_D2y_Py_dd;
  Double I_ERI_F2xz_S_D2y_Dyz_dd = I_ERI_F2xz_S_F2yz_Py_dd+CDZ*I_ERI_F2xz_S_D2y_Py_dd;
  Double I_ERI_Fx2y_S_D2y_Dyz_dd = I_ERI_Fx2y_S_F2yz_Py_dd+CDZ*I_ERI_Fx2y_S_D2y_Py_dd;
  Double I_ERI_Fxyz_S_D2y_Dyz_dd = I_ERI_Fxyz_S_F2yz_Py_dd+CDZ*I_ERI_Fxyz_S_D2y_Py_dd;
  Double I_ERI_Fx2z_S_D2y_Dyz_dd = I_ERI_Fx2z_S_F2yz_Py_dd+CDZ*I_ERI_Fx2z_S_D2y_Py_dd;
  Double I_ERI_F3y_S_D2y_Dyz_dd = I_ERI_F3y_S_F2yz_Py_dd+CDZ*I_ERI_F3y_S_D2y_Py_dd;
  Double I_ERI_F2yz_S_D2y_Dyz_dd = I_ERI_F2yz_S_F2yz_Py_dd+CDZ*I_ERI_F2yz_S_D2y_Py_dd;
  Double I_ERI_Fy2z_S_D2y_Dyz_dd = I_ERI_Fy2z_S_F2yz_Py_dd+CDZ*I_ERI_Fy2z_S_D2y_Py_dd;
  Double I_ERI_F3z_S_D2y_Dyz_dd = I_ERI_F3z_S_F2yz_Py_dd+CDZ*I_ERI_F3z_S_D2y_Py_dd;
  Double I_ERI_F3x_S_Dyz_Dyz_dd = I_ERI_F3x_S_Fy2z_Py_dd+CDZ*I_ERI_F3x_S_Dyz_Py_dd;
  Double I_ERI_F2xy_S_Dyz_Dyz_dd = I_ERI_F2xy_S_Fy2z_Py_dd+CDZ*I_ERI_F2xy_S_Dyz_Py_dd;
  Double I_ERI_F2xz_S_Dyz_Dyz_dd = I_ERI_F2xz_S_Fy2z_Py_dd+CDZ*I_ERI_F2xz_S_Dyz_Py_dd;
  Double I_ERI_Fx2y_S_Dyz_Dyz_dd = I_ERI_Fx2y_S_Fy2z_Py_dd+CDZ*I_ERI_Fx2y_S_Dyz_Py_dd;
  Double I_ERI_Fxyz_S_Dyz_Dyz_dd = I_ERI_Fxyz_S_Fy2z_Py_dd+CDZ*I_ERI_Fxyz_S_Dyz_Py_dd;
  Double I_ERI_Fx2z_S_Dyz_Dyz_dd = I_ERI_Fx2z_S_Fy2z_Py_dd+CDZ*I_ERI_Fx2z_S_Dyz_Py_dd;
  Double I_ERI_F3y_S_Dyz_Dyz_dd = I_ERI_F3y_S_Fy2z_Py_dd+CDZ*I_ERI_F3y_S_Dyz_Py_dd;
  Double I_ERI_F2yz_S_Dyz_Dyz_dd = I_ERI_F2yz_S_Fy2z_Py_dd+CDZ*I_ERI_F2yz_S_Dyz_Py_dd;
  Double I_ERI_Fy2z_S_Dyz_Dyz_dd = I_ERI_Fy2z_S_Fy2z_Py_dd+CDZ*I_ERI_Fy2z_S_Dyz_Py_dd;
  Double I_ERI_F3z_S_Dyz_Dyz_dd = I_ERI_F3z_S_Fy2z_Py_dd+CDZ*I_ERI_F3z_S_Dyz_Py_dd;
  Double I_ERI_F3x_S_D2z_Dyz_dd = I_ERI_F3x_S_F3z_Py_dd+CDZ*I_ERI_F3x_S_D2z_Py_dd;
  Double I_ERI_F2xy_S_D2z_Dyz_dd = I_ERI_F2xy_S_F3z_Py_dd+CDZ*I_ERI_F2xy_S_D2z_Py_dd;
  Double I_ERI_F2xz_S_D2z_Dyz_dd = I_ERI_F2xz_S_F3z_Py_dd+CDZ*I_ERI_F2xz_S_D2z_Py_dd;
  Double I_ERI_Fx2y_S_D2z_Dyz_dd = I_ERI_Fx2y_S_F3z_Py_dd+CDZ*I_ERI_Fx2y_S_D2z_Py_dd;
  Double I_ERI_Fxyz_S_D2z_Dyz_dd = I_ERI_Fxyz_S_F3z_Py_dd+CDZ*I_ERI_Fxyz_S_D2z_Py_dd;
  Double I_ERI_Fx2z_S_D2z_Dyz_dd = I_ERI_Fx2z_S_F3z_Py_dd+CDZ*I_ERI_Fx2z_S_D2z_Py_dd;
  Double I_ERI_F3y_S_D2z_Dyz_dd = I_ERI_F3y_S_F3z_Py_dd+CDZ*I_ERI_F3y_S_D2z_Py_dd;
  Double I_ERI_F2yz_S_D2z_Dyz_dd = I_ERI_F2yz_S_F3z_Py_dd+CDZ*I_ERI_F2yz_S_D2z_Py_dd;
  Double I_ERI_Fy2z_S_D2z_Dyz_dd = I_ERI_Fy2z_S_F3z_Py_dd+CDZ*I_ERI_Fy2z_S_D2z_Py_dd;
  Double I_ERI_F3z_S_D2z_Dyz_dd = I_ERI_F3z_S_F3z_Py_dd+CDZ*I_ERI_F3z_S_D2z_Py_dd;
  Double I_ERI_F3x_S_D2x_D2z_dd = I_ERI_F3x_S_F2xz_Pz_dd+CDZ*I_ERI_F3x_S_D2x_Pz_dd;
  Double I_ERI_F2xy_S_D2x_D2z_dd = I_ERI_F2xy_S_F2xz_Pz_dd+CDZ*I_ERI_F2xy_S_D2x_Pz_dd;
  Double I_ERI_F2xz_S_D2x_D2z_dd = I_ERI_F2xz_S_F2xz_Pz_dd+CDZ*I_ERI_F2xz_S_D2x_Pz_dd;
  Double I_ERI_Fx2y_S_D2x_D2z_dd = I_ERI_Fx2y_S_F2xz_Pz_dd+CDZ*I_ERI_Fx2y_S_D2x_Pz_dd;
  Double I_ERI_Fxyz_S_D2x_D2z_dd = I_ERI_Fxyz_S_F2xz_Pz_dd+CDZ*I_ERI_Fxyz_S_D2x_Pz_dd;
  Double I_ERI_Fx2z_S_D2x_D2z_dd = I_ERI_Fx2z_S_F2xz_Pz_dd+CDZ*I_ERI_Fx2z_S_D2x_Pz_dd;
  Double I_ERI_F3y_S_D2x_D2z_dd = I_ERI_F3y_S_F2xz_Pz_dd+CDZ*I_ERI_F3y_S_D2x_Pz_dd;
  Double I_ERI_F2yz_S_D2x_D2z_dd = I_ERI_F2yz_S_F2xz_Pz_dd+CDZ*I_ERI_F2yz_S_D2x_Pz_dd;
  Double I_ERI_Fy2z_S_D2x_D2z_dd = I_ERI_Fy2z_S_F2xz_Pz_dd+CDZ*I_ERI_Fy2z_S_D2x_Pz_dd;
  Double I_ERI_F3z_S_D2x_D2z_dd = I_ERI_F3z_S_F2xz_Pz_dd+CDZ*I_ERI_F3z_S_D2x_Pz_dd;
  Double I_ERI_F3x_S_Dxy_D2z_dd = I_ERI_F3x_S_Fxyz_Pz_dd+CDZ*I_ERI_F3x_S_Dxy_Pz_dd;
  Double I_ERI_F2xy_S_Dxy_D2z_dd = I_ERI_F2xy_S_Fxyz_Pz_dd+CDZ*I_ERI_F2xy_S_Dxy_Pz_dd;
  Double I_ERI_F2xz_S_Dxy_D2z_dd = I_ERI_F2xz_S_Fxyz_Pz_dd+CDZ*I_ERI_F2xz_S_Dxy_Pz_dd;
  Double I_ERI_Fx2y_S_Dxy_D2z_dd = I_ERI_Fx2y_S_Fxyz_Pz_dd+CDZ*I_ERI_Fx2y_S_Dxy_Pz_dd;
  Double I_ERI_Fxyz_S_Dxy_D2z_dd = I_ERI_Fxyz_S_Fxyz_Pz_dd+CDZ*I_ERI_Fxyz_S_Dxy_Pz_dd;
  Double I_ERI_Fx2z_S_Dxy_D2z_dd = I_ERI_Fx2z_S_Fxyz_Pz_dd+CDZ*I_ERI_Fx2z_S_Dxy_Pz_dd;
  Double I_ERI_F3y_S_Dxy_D2z_dd = I_ERI_F3y_S_Fxyz_Pz_dd+CDZ*I_ERI_F3y_S_Dxy_Pz_dd;
  Double I_ERI_F2yz_S_Dxy_D2z_dd = I_ERI_F2yz_S_Fxyz_Pz_dd+CDZ*I_ERI_F2yz_S_Dxy_Pz_dd;
  Double I_ERI_Fy2z_S_Dxy_D2z_dd = I_ERI_Fy2z_S_Fxyz_Pz_dd+CDZ*I_ERI_Fy2z_S_Dxy_Pz_dd;
  Double I_ERI_F3z_S_Dxy_D2z_dd = I_ERI_F3z_S_Fxyz_Pz_dd+CDZ*I_ERI_F3z_S_Dxy_Pz_dd;
  Double I_ERI_F3x_S_Dxz_D2z_dd = I_ERI_F3x_S_Fx2z_Pz_dd+CDZ*I_ERI_F3x_S_Dxz_Pz_dd;
  Double I_ERI_F2xy_S_Dxz_D2z_dd = I_ERI_F2xy_S_Fx2z_Pz_dd+CDZ*I_ERI_F2xy_S_Dxz_Pz_dd;
  Double I_ERI_F2xz_S_Dxz_D2z_dd = I_ERI_F2xz_S_Fx2z_Pz_dd+CDZ*I_ERI_F2xz_S_Dxz_Pz_dd;
  Double I_ERI_Fx2y_S_Dxz_D2z_dd = I_ERI_Fx2y_S_Fx2z_Pz_dd+CDZ*I_ERI_Fx2y_S_Dxz_Pz_dd;
  Double I_ERI_Fxyz_S_Dxz_D2z_dd = I_ERI_Fxyz_S_Fx2z_Pz_dd+CDZ*I_ERI_Fxyz_S_Dxz_Pz_dd;
  Double I_ERI_Fx2z_S_Dxz_D2z_dd = I_ERI_Fx2z_S_Fx2z_Pz_dd+CDZ*I_ERI_Fx2z_S_Dxz_Pz_dd;
  Double I_ERI_F3y_S_Dxz_D2z_dd = I_ERI_F3y_S_Fx2z_Pz_dd+CDZ*I_ERI_F3y_S_Dxz_Pz_dd;
  Double I_ERI_F2yz_S_Dxz_D2z_dd = I_ERI_F2yz_S_Fx2z_Pz_dd+CDZ*I_ERI_F2yz_S_Dxz_Pz_dd;
  Double I_ERI_Fy2z_S_Dxz_D2z_dd = I_ERI_Fy2z_S_Fx2z_Pz_dd+CDZ*I_ERI_Fy2z_S_Dxz_Pz_dd;
  Double I_ERI_F3z_S_Dxz_D2z_dd = I_ERI_F3z_S_Fx2z_Pz_dd+CDZ*I_ERI_F3z_S_Dxz_Pz_dd;
  Double I_ERI_F3x_S_D2y_D2z_dd = I_ERI_F3x_S_F2yz_Pz_dd+CDZ*I_ERI_F3x_S_D2y_Pz_dd;
  Double I_ERI_F2xy_S_D2y_D2z_dd = I_ERI_F2xy_S_F2yz_Pz_dd+CDZ*I_ERI_F2xy_S_D2y_Pz_dd;
  Double I_ERI_F2xz_S_D2y_D2z_dd = I_ERI_F2xz_S_F2yz_Pz_dd+CDZ*I_ERI_F2xz_S_D2y_Pz_dd;
  Double I_ERI_Fx2y_S_D2y_D2z_dd = I_ERI_Fx2y_S_F2yz_Pz_dd+CDZ*I_ERI_Fx2y_S_D2y_Pz_dd;
  Double I_ERI_Fxyz_S_D2y_D2z_dd = I_ERI_Fxyz_S_F2yz_Pz_dd+CDZ*I_ERI_Fxyz_S_D2y_Pz_dd;
  Double I_ERI_Fx2z_S_D2y_D2z_dd = I_ERI_Fx2z_S_F2yz_Pz_dd+CDZ*I_ERI_Fx2z_S_D2y_Pz_dd;
  Double I_ERI_F3y_S_D2y_D2z_dd = I_ERI_F3y_S_F2yz_Pz_dd+CDZ*I_ERI_F3y_S_D2y_Pz_dd;
  Double I_ERI_F2yz_S_D2y_D2z_dd = I_ERI_F2yz_S_F2yz_Pz_dd+CDZ*I_ERI_F2yz_S_D2y_Pz_dd;
  Double I_ERI_Fy2z_S_D2y_D2z_dd = I_ERI_Fy2z_S_F2yz_Pz_dd+CDZ*I_ERI_Fy2z_S_D2y_Pz_dd;
  Double I_ERI_F3z_S_D2y_D2z_dd = I_ERI_F3z_S_F2yz_Pz_dd+CDZ*I_ERI_F3z_S_D2y_Pz_dd;
  Double I_ERI_F3x_S_Dyz_D2z_dd = I_ERI_F3x_S_Fy2z_Pz_dd+CDZ*I_ERI_F3x_S_Dyz_Pz_dd;
  Double I_ERI_F2xy_S_Dyz_D2z_dd = I_ERI_F2xy_S_Fy2z_Pz_dd+CDZ*I_ERI_F2xy_S_Dyz_Pz_dd;
  Double I_ERI_F2xz_S_Dyz_D2z_dd = I_ERI_F2xz_S_Fy2z_Pz_dd+CDZ*I_ERI_F2xz_S_Dyz_Pz_dd;
  Double I_ERI_Fx2y_S_Dyz_D2z_dd = I_ERI_Fx2y_S_Fy2z_Pz_dd+CDZ*I_ERI_Fx2y_S_Dyz_Pz_dd;
  Double I_ERI_Fxyz_S_Dyz_D2z_dd = I_ERI_Fxyz_S_Fy2z_Pz_dd+CDZ*I_ERI_Fxyz_S_Dyz_Pz_dd;
  Double I_ERI_Fx2z_S_Dyz_D2z_dd = I_ERI_Fx2z_S_Fy2z_Pz_dd+CDZ*I_ERI_Fx2z_S_Dyz_Pz_dd;
  Double I_ERI_F3y_S_Dyz_D2z_dd = I_ERI_F3y_S_Fy2z_Pz_dd+CDZ*I_ERI_F3y_S_Dyz_Pz_dd;
  Double I_ERI_F2yz_S_Dyz_D2z_dd = I_ERI_F2yz_S_Fy2z_Pz_dd+CDZ*I_ERI_F2yz_S_Dyz_Pz_dd;
  Double I_ERI_Fy2z_S_Dyz_D2z_dd = I_ERI_Fy2z_S_Fy2z_Pz_dd+CDZ*I_ERI_Fy2z_S_Dyz_Pz_dd;
  Double I_ERI_F3z_S_Dyz_D2z_dd = I_ERI_F3z_S_Fy2z_Pz_dd+CDZ*I_ERI_F3z_S_Dyz_Pz_dd;
  Double I_ERI_F3x_S_D2z_D2z_dd = I_ERI_F3x_S_F3z_Pz_dd+CDZ*I_ERI_F3x_S_D2z_Pz_dd;
  Double I_ERI_F2xy_S_D2z_D2z_dd = I_ERI_F2xy_S_F3z_Pz_dd+CDZ*I_ERI_F2xy_S_D2z_Pz_dd;
  Double I_ERI_F2xz_S_D2z_D2z_dd = I_ERI_F2xz_S_F3z_Pz_dd+CDZ*I_ERI_F2xz_S_D2z_Pz_dd;
  Double I_ERI_Fx2y_S_D2z_D2z_dd = I_ERI_Fx2y_S_F3z_Pz_dd+CDZ*I_ERI_Fx2y_S_D2z_Pz_dd;
  Double I_ERI_Fxyz_S_D2z_D2z_dd = I_ERI_Fxyz_S_F3z_Pz_dd+CDZ*I_ERI_Fxyz_S_D2z_Pz_dd;
  Double I_ERI_Fx2z_S_D2z_D2z_dd = I_ERI_Fx2z_S_F3z_Pz_dd+CDZ*I_ERI_Fx2z_S_D2z_Pz_dd;
  Double I_ERI_F3y_S_D2z_D2z_dd = I_ERI_F3y_S_F3z_Pz_dd+CDZ*I_ERI_F3y_S_D2z_Pz_dd;
  Double I_ERI_F2yz_S_D2z_D2z_dd = I_ERI_F2yz_S_F3z_Pz_dd+CDZ*I_ERI_F2yz_S_D2z_Pz_dd;
  Double I_ERI_Fy2z_S_D2z_D2z_dd = I_ERI_Fy2z_S_F3z_Pz_dd+CDZ*I_ERI_Fy2z_S_D2z_Pz_dd;
  Double I_ERI_F3z_S_D2z_D2z_dd = I_ERI_F3z_S_F3z_Pz_dd+CDZ*I_ERI_F3z_S_D2z_Pz_dd;

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
   * shell quartet name: SQ_ERI_F_P_P_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_b
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_b = I_ERI_G4x_S_Px_S_b+ABX*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_Px_Px_S_b = I_ERI_G3xy_S_Px_S_b+ABX*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_Px_Px_S_b = I_ERI_G3xz_S_Px_S_b+ABX*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_Px_Px_S_b = I_ERI_G2x2y_S_Px_S_b+ABX*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_Px_Px_S_b = I_ERI_G2xyz_S_Px_S_b+ABX*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_Px_Px_S_b = I_ERI_G2x2z_S_Px_S_b+ABX*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_Px_Px_S_b = I_ERI_Gx3y_S_Px_S_b+ABX*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_Px_Px_S_b = I_ERI_Gx2yz_S_Px_S_b+ABX*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_Px_Px_S_b = I_ERI_Gxy2z_S_Px_S_b+ABX*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_Px_Px_S_b = I_ERI_Gx3z_S_Px_S_b+ABX*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_Py_Px_S_b = I_ERI_G3xy_S_Px_S_b+ABY*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_Py_Px_S_b = I_ERI_G2x2y_S_Px_S_b+ABY*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_Py_Px_S_b = I_ERI_G2xyz_S_Px_S_b+ABY*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_Py_Px_S_b = I_ERI_Gx3y_S_Px_S_b+ABY*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_Py_Px_S_b = I_ERI_Gx2yz_S_Px_S_b+ABY*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_Py_Px_S_b = I_ERI_Gxy2z_S_Px_S_b+ABY*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_Py_Px_S_b = I_ERI_G4y_S_Px_S_b+ABY*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_Py_Px_S_b = I_ERI_G3yz_S_Px_S_b+ABY*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_Py_Px_S_b = I_ERI_G2y2z_S_Px_S_b+ABY*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_Py_Px_S_b = I_ERI_Gy3z_S_Px_S_b+ABY*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_Pz_Px_S_b = I_ERI_G3xz_S_Px_S_b+ABZ*I_ERI_F3x_S_Px_S_b;
  Double I_ERI_F2xy_Pz_Px_S_b = I_ERI_G2xyz_S_Px_S_b+ABZ*I_ERI_F2xy_S_Px_S_b;
  Double I_ERI_F2xz_Pz_Px_S_b = I_ERI_G2x2z_S_Px_S_b+ABZ*I_ERI_F2xz_S_Px_S_b;
  Double I_ERI_Fx2y_Pz_Px_S_b = I_ERI_Gx2yz_S_Px_S_b+ABZ*I_ERI_Fx2y_S_Px_S_b;
  Double I_ERI_Fxyz_Pz_Px_S_b = I_ERI_Gxy2z_S_Px_S_b+ABZ*I_ERI_Fxyz_S_Px_S_b;
  Double I_ERI_Fx2z_Pz_Px_S_b = I_ERI_Gx3z_S_Px_S_b+ABZ*I_ERI_Fx2z_S_Px_S_b;
  Double I_ERI_F3y_Pz_Px_S_b = I_ERI_G3yz_S_Px_S_b+ABZ*I_ERI_F3y_S_Px_S_b;
  Double I_ERI_F2yz_Pz_Px_S_b = I_ERI_G2y2z_S_Px_S_b+ABZ*I_ERI_F2yz_S_Px_S_b;
  Double I_ERI_Fy2z_Pz_Px_S_b = I_ERI_Gy3z_S_Px_S_b+ABZ*I_ERI_Fy2z_S_Px_S_b;
  Double I_ERI_F3z_Pz_Px_S_b = I_ERI_G4z_S_Px_S_b+ABZ*I_ERI_F3z_S_Px_S_b;
  Double I_ERI_F3x_Px_Py_S_b = I_ERI_G4x_S_Py_S_b+ABX*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_Px_Py_S_b = I_ERI_G3xy_S_Py_S_b+ABX*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_Px_Py_S_b = I_ERI_G3xz_S_Py_S_b+ABX*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_Px_Py_S_b = I_ERI_G2x2y_S_Py_S_b+ABX*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_Px_Py_S_b = I_ERI_G2xyz_S_Py_S_b+ABX*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_Px_Py_S_b = I_ERI_G2x2z_S_Py_S_b+ABX*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_Px_Py_S_b = I_ERI_Gx3y_S_Py_S_b+ABX*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_Px_Py_S_b = I_ERI_Gx2yz_S_Py_S_b+ABX*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_Px_Py_S_b = I_ERI_Gxy2z_S_Py_S_b+ABX*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_Px_Py_S_b = I_ERI_Gx3z_S_Py_S_b+ABX*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_Py_Py_S_b = I_ERI_G3xy_S_Py_S_b+ABY*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_Py_Py_S_b = I_ERI_G2x2y_S_Py_S_b+ABY*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_Py_Py_S_b = I_ERI_G2xyz_S_Py_S_b+ABY*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_Py_Py_S_b = I_ERI_Gx3y_S_Py_S_b+ABY*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_Py_Py_S_b = I_ERI_Gx2yz_S_Py_S_b+ABY*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_Py_Py_S_b = I_ERI_Gxy2z_S_Py_S_b+ABY*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_Py_Py_S_b = I_ERI_G4y_S_Py_S_b+ABY*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_Py_Py_S_b = I_ERI_G3yz_S_Py_S_b+ABY*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_Py_Py_S_b = I_ERI_G2y2z_S_Py_S_b+ABY*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_Py_Py_S_b = I_ERI_Gy3z_S_Py_S_b+ABY*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_Pz_Py_S_b = I_ERI_G3xz_S_Py_S_b+ABZ*I_ERI_F3x_S_Py_S_b;
  Double I_ERI_F2xy_Pz_Py_S_b = I_ERI_G2xyz_S_Py_S_b+ABZ*I_ERI_F2xy_S_Py_S_b;
  Double I_ERI_F2xz_Pz_Py_S_b = I_ERI_G2x2z_S_Py_S_b+ABZ*I_ERI_F2xz_S_Py_S_b;
  Double I_ERI_Fx2y_Pz_Py_S_b = I_ERI_Gx2yz_S_Py_S_b+ABZ*I_ERI_Fx2y_S_Py_S_b;
  Double I_ERI_Fxyz_Pz_Py_S_b = I_ERI_Gxy2z_S_Py_S_b+ABZ*I_ERI_Fxyz_S_Py_S_b;
  Double I_ERI_Fx2z_Pz_Py_S_b = I_ERI_Gx3z_S_Py_S_b+ABZ*I_ERI_Fx2z_S_Py_S_b;
  Double I_ERI_F3y_Pz_Py_S_b = I_ERI_G3yz_S_Py_S_b+ABZ*I_ERI_F3y_S_Py_S_b;
  Double I_ERI_F2yz_Pz_Py_S_b = I_ERI_G2y2z_S_Py_S_b+ABZ*I_ERI_F2yz_S_Py_S_b;
  Double I_ERI_Fy2z_Pz_Py_S_b = I_ERI_Gy3z_S_Py_S_b+ABZ*I_ERI_Fy2z_S_Py_S_b;
  Double I_ERI_F3z_Pz_Py_S_b = I_ERI_G4z_S_Py_S_b+ABZ*I_ERI_F3z_S_Py_S_b;
  Double I_ERI_F3x_Px_Pz_S_b = I_ERI_G4x_S_Pz_S_b+ABX*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_Px_Pz_S_b = I_ERI_G3xy_S_Pz_S_b+ABX*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_Px_Pz_S_b = I_ERI_G3xz_S_Pz_S_b+ABX*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_Px_Pz_S_b = I_ERI_G2x2y_S_Pz_S_b+ABX*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_Px_Pz_S_b = I_ERI_G2xyz_S_Pz_S_b+ABX*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_Px_Pz_S_b = I_ERI_G2x2z_S_Pz_S_b+ABX*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_Px_Pz_S_b = I_ERI_Gx3y_S_Pz_S_b+ABX*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_Px_Pz_S_b = I_ERI_Gx2yz_S_Pz_S_b+ABX*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_Px_Pz_S_b = I_ERI_Gxy2z_S_Pz_S_b+ABX*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_Px_Pz_S_b = I_ERI_Gx3z_S_Pz_S_b+ABX*I_ERI_F3z_S_Pz_S_b;
  Double I_ERI_F3x_Py_Pz_S_b = I_ERI_G3xy_S_Pz_S_b+ABY*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_Py_Pz_S_b = I_ERI_G2x2y_S_Pz_S_b+ABY*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_Py_Pz_S_b = I_ERI_G2xyz_S_Pz_S_b+ABY*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_Py_Pz_S_b = I_ERI_Gx3y_S_Pz_S_b+ABY*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_Py_Pz_S_b = I_ERI_Gx2yz_S_Pz_S_b+ABY*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_Py_Pz_S_b = I_ERI_Gxy2z_S_Pz_S_b+ABY*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_Py_Pz_S_b = I_ERI_G4y_S_Pz_S_b+ABY*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_Py_Pz_S_b = I_ERI_G3yz_S_Pz_S_b+ABY*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_Py_Pz_S_b = I_ERI_G2y2z_S_Pz_S_b+ABY*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_Py_Pz_S_b = I_ERI_Gy3z_S_Pz_S_b+ABY*I_ERI_F3z_S_Pz_S_b;
  Double I_ERI_F3x_Pz_Pz_S_b = I_ERI_G3xz_S_Pz_S_b+ABZ*I_ERI_F3x_S_Pz_S_b;
  Double I_ERI_F2xy_Pz_Pz_S_b = I_ERI_G2xyz_S_Pz_S_b+ABZ*I_ERI_F2xy_S_Pz_S_b;
  Double I_ERI_F2xz_Pz_Pz_S_b = I_ERI_G2x2z_S_Pz_S_b+ABZ*I_ERI_F2xz_S_Pz_S_b;
  Double I_ERI_Fx2y_Pz_Pz_S_b = I_ERI_Gx2yz_S_Pz_S_b+ABZ*I_ERI_Fx2y_S_Pz_S_b;
  Double I_ERI_Fxyz_Pz_Pz_S_b = I_ERI_Gxy2z_S_Pz_S_b+ABZ*I_ERI_Fxyz_S_Pz_S_b;
  Double I_ERI_Fx2z_Pz_Pz_S_b = I_ERI_Gx3z_S_Pz_S_b+ABZ*I_ERI_Fx2z_S_Pz_S_b;
  Double I_ERI_F3y_Pz_Pz_S_b = I_ERI_G3yz_S_Pz_S_b+ABZ*I_ERI_F3y_S_Pz_S_b;
  Double I_ERI_F2yz_Pz_Pz_S_b = I_ERI_G2y2z_S_Pz_S_b+ABZ*I_ERI_F2yz_S_Pz_S_b;
  Double I_ERI_Fy2z_Pz_Pz_S_b = I_ERI_Gy3z_S_Pz_S_b+ABZ*I_ERI_Fy2z_S_Pz_S_b;
  Double I_ERI_F3z_Pz_Pz_S_b = I_ERI_G4z_S_Pz_S_b+ABZ*I_ERI_F3z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bb
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_bb = I_ERI_G4x_S_D2x_S_bb+ABX*I_ERI_F3x_S_D2x_S_bb;
  Double I_ERI_F2xy_Px_D2x_S_bb = I_ERI_G3xy_S_D2x_S_bb+ABX*I_ERI_F2xy_S_D2x_S_bb;
  Double I_ERI_F2xz_Px_D2x_S_bb = I_ERI_G3xz_S_D2x_S_bb+ABX*I_ERI_F2xz_S_D2x_S_bb;
  Double I_ERI_Fx2y_Px_D2x_S_bb = I_ERI_G2x2y_S_D2x_S_bb+ABX*I_ERI_Fx2y_S_D2x_S_bb;
  Double I_ERI_Fxyz_Px_D2x_S_bb = I_ERI_G2xyz_S_D2x_S_bb+ABX*I_ERI_Fxyz_S_D2x_S_bb;
  Double I_ERI_Fx2z_Px_D2x_S_bb = I_ERI_G2x2z_S_D2x_S_bb+ABX*I_ERI_Fx2z_S_D2x_S_bb;
  Double I_ERI_F3y_Px_D2x_S_bb = I_ERI_Gx3y_S_D2x_S_bb+ABX*I_ERI_F3y_S_D2x_S_bb;
  Double I_ERI_F2yz_Px_D2x_S_bb = I_ERI_Gx2yz_S_D2x_S_bb+ABX*I_ERI_F2yz_S_D2x_S_bb;
  Double I_ERI_Fy2z_Px_D2x_S_bb = I_ERI_Gxy2z_S_D2x_S_bb+ABX*I_ERI_Fy2z_S_D2x_S_bb;
  Double I_ERI_F3z_Px_D2x_S_bb = I_ERI_Gx3z_S_D2x_S_bb+ABX*I_ERI_F3z_S_D2x_S_bb;
  Double I_ERI_F3x_Py_D2x_S_bb = I_ERI_G3xy_S_D2x_S_bb+ABY*I_ERI_F3x_S_D2x_S_bb;
  Double I_ERI_F2xy_Py_D2x_S_bb = I_ERI_G2x2y_S_D2x_S_bb+ABY*I_ERI_F2xy_S_D2x_S_bb;
  Double I_ERI_F2xz_Py_D2x_S_bb = I_ERI_G2xyz_S_D2x_S_bb+ABY*I_ERI_F2xz_S_D2x_S_bb;
  Double I_ERI_Fx2y_Py_D2x_S_bb = I_ERI_Gx3y_S_D2x_S_bb+ABY*I_ERI_Fx2y_S_D2x_S_bb;
  Double I_ERI_Fxyz_Py_D2x_S_bb = I_ERI_Gx2yz_S_D2x_S_bb+ABY*I_ERI_Fxyz_S_D2x_S_bb;
  Double I_ERI_Fx2z_Py_D2x_S_bb = I_ERI_Gxy2z_S_D2x_S_bb+ABY*I_ERI_Fx2z_S_D2x_S_bb;
  Double I_ERI_F3y_Py_D2x_S_bb = I_ERI_G4y_S_D2x_S_bb+ABY*I_ERI_F3y_S_D2x_S_bb;
  Double I_ERI_F2yz_Py_D2x_S_bb = I_ERI_G3yz_S_D2x_S_bb+ABY*I_ERI_F2yz_S_D2x_S_bb;
  Double I_ERI_Fy2z_Py_D2x_S_bb = I_ERI_G2y2z_S_D2x_S_bb+ABY*I_ERI_Fy2z_S_D2x_S_bb;
  Double I_ERI_F3z_Py_D2x_S_bb = I_ERI_Gy3z_S_D2x_S_bb+ABY*I_ERI_F3z_S_D2x_S_bb;
  Double I_ERI_F3x_Pz_D2x_S_bb = I_ERI_G3xz_S_D2x_S_bb+ABZ*I_ERI_F3x_S_D2x_S_bb;
  Double I_ERI_F2xy_Pz_D2x_S_bb = I_ERI_G2xyz_S_D2x_S_bb+ABZ*I_ERI_F2xy_S_D2x_S_bb;
  Double I_ERI_F2xz_Pz_D2x_S_bb = I_ERI_G2x2z_S_D2x_S_bb+ABZ*I_ERI_F2xz_S_D2x_S_bb;
  Double I_ERI_Fx2y_Pz_D2x_S_bb = I_ERI_Gx2yz_S_D2x_S_bb+ABZ*I_ERI_Fx2y_S_D2x_S_bb;
  Double I_ERI_Fxyz_Pz_D2x_S_bb = I_ERI_Gxy2z_S_D2x_S_bb+ABZ*I_ERI_Fxyz_S_D2x_S_bb;
  Double I_ERI_Fx2z_Pz_D2x_S_bb = I_ERI_Gx3z_S_D2x_S_bb+ABZ*I_ERI_Fx2z_S_D2x_S_bb;
  Double I_ERI_F3y_Pz_D2x_S_bb = I_ERI_G3yz_S_D2x_S_bb+ABZ*I_ERI_F3y_S_D2x_S_bb;
  Double I_ERI_F2yz_Pz_D2x_S_bb = I_ERI_G2y2z_S_D2x_S_bb+ABZ*I_ERI_F2yz_S_D2x_S_bb;
  Double I_ERI_Fy2z_Pz_D2x_S_bb = I_ERI_Gy3z_S_D2x_S_bb+ABZ*I_ERI_Fy2z_S_D2x_S_bb;
  Double I_ERI_F3z_Pz_D2x_S_bb = I_ERI_G4z_S_D2x_S_bb+ABZ*I_ERI_F3z_S_D2x_S_bb;
  Double I_ERI_F3x_Px_Dxy_S_bb = I_ERI_G4x_S_Dxy_S_bb+ABX*I_ERI_F3x_S_Dxy_S_bb;
  Double I_ERI_F2xy_Px_Dxy_S_bb = I_ERI_G3xy_S_Dxy_S_bb+ABX*I_ERI_F2xy_S_Dxy_S_bb;
  Double I_ERI_F2xz_Px_Dxy_S_bb = I_ERI_G3xz_S_Dxy_S_bb+ABX*I_ERI_F2xz_S_Dxy_S_bb;
  Double I_ERI_Fx2y_Px_Dxy_S_bb = I_ERI_G2x2y_S_Dxy_S_bb+ABX*I_ERI_Fx2y_S_Dxy_S_bb;
  Double I_ERI_Fxyz_Px_Dxy_S_bb = I_ERI_G2xyz_S_Dxy_S_bb+ABX*I_ERI_Fxyz_S_Dxy_S_bb;
  Double I_ERI_Fx2z_Px_Dxy_S_bb = I_ERI_G2x2z_S_Dxy_S_bb+ABX*I_ERI_Fx2z_S_Dxy_S_bb;
  Double I_ERI_F3y_Px_Dxy_S_bb = I_ERI_Gx3y_S_Dxy_S_bb+ABX*I_ERI_F3y_S_Dxy_S_bb;
  Double I_ERI_F2yz_Px_Dxy_S_bb = I_ERI_Gx2yz_S_Dxy_S_bb+ABX*I_ERI_F2yz_S_Dxy_S_bb;
  Double I_ERI_Fy2z_Px_Dxy_S_bb = I_ERI_Gxy2z_S_Dxy_S_bb+ABX*I_ERI_Fy2z_S_Dxy_S_bb;
  Double I_ERI_F3z_Px_Dxy_S_bb = I_ERI_Gx3z_S_Dxy_S_bb+ABX*I_ERI_F3z_S_Dxy_S_bb;
  Double I_ERI_F3x_Py_Dxy_S_bb = I_ERI_G3xy_S_Dxy_S_bb+ABY*I_ERI_F3x_S_Dxy_S_bb;
  Double I_ERI_F2xy_Py_Dxy_S_bb = I_ERI_G2x2y_S_Dxy_S_bb+ABY*I_ERI_F2xy_S_Dxy_S_bb;
  Double I_ERI_F2xz_Py_Dxy_S_bb = I_ERI_G2xyz_S_Dxy_S_bb+ABY*I_ERI_F2xz_S_Dxy_S_bb;
  Double I_ERI_Fx2y_Py_Dxy_S_bb = I_ERI_Gx3y_S_Dxy_S_bb+ABY*I_ERI_Fx2y_S_Dxy_S_bb;
  Double I_ERI_Fxyz_Py_Dxy_S_bb = I_ERI_Gx2yz_S_Dxy_S_bb+ABY*I_ERI_Fxyz_S_Dxy_S_bb;
  Double I_ERI_Fx2z_Py_Dxy_S_bb = I_ERI_Gxy2z_S_Dxy_S_bb+ABY*I_ERI_Fx2z_S_Dxy_S_bb;
  Double I_ERI_F3y_Py_Dxy_S_bb = I_ERI_G4y_S_Dxy_S_bb+ABY*I_ERI_F3y_S_Dxy_S_bb;
  Double I_ERI_F2yz_Py_Dxy_S_bb = I_ERI_G3yz_S_Dxy_S_bb+ABY*I_ERI_F2yz_S_Dxy_S_bb;
  Double I_ERI_Fy2z_Py_Dxy_S_bb = I_ERI_G2y2z_S_Dxy_S_bb+ABY*I_ERI_Fy2z_S_Dxy_S_bb;
  Double I_ERI_F3z_Py_Dxy_S_bb = I_ERI_Gy3z_S_Dxy_S_bb+ABY*I_ERI_F3z_S_Dxy_S_bb;
  Double I_ERI_F3x_Pz_Dxy_S_bb = I_ERI_G3xz_S_Dxy_S_bb+ABZ*I_ERI_F3x_S_Dxy_S_bb;
  Double I_ERI_F2xy_Pz_Dxy_S_bb = I_ERI_G2xyz_S_Dxy_S_bb+ABZ*I_ERI_F2xy_S_Dxy_S_bb;
  Double I_ERI_F2xz_Pz_Dxy_S_bb = I_ERI_G2x2z_S_Dxy_S_bb+ABZ*I_ERI_F2xz_S_Dxy_S_bb;
  Double I_ERI_Fx2y_Pz_Dxy_S_bb = I_ERI_Gx2yz_S_Dxy_S_bb+ABZ*I_ERI_Fx2y_S_Dxy_S_bb;
  Double I_ERI_Fxyz_Pz_Dxy_S_bb = I_ERI_Gxy2z_S_Dxy_S_bb+ABZ*I_ERI_Fxyz_S_Dxy_S_bb;
  Double I_ERI_Fx2z_Pz_Dxy_S_bb = I_ERI_Gx3z_S_Dxy_S_bb+ABZ*I_ERI_Fx2z_S_Dxy_S_bb;
  Double I_ERI_F3y_Pz_Dxy_S_bb = I_ERI_G3yz_S_Dxy_S_bb+ABZ*I_ERI_F3y_S_Dxy_S_bb;
  Double I_ERI_F2yz_Pz_Dxy_S_bb = I_ERI_G2y2z_S_Dxy_S_bb+ABZ*I_ERI_F2yz_S_Dxy_S_bb;
  Double I_ERI_Fy2z_Pz_Dxy_S_bb = I_ERI_Gy3z_S_Dxy_S_bb+ABZ*I_ERI_Fy2z_S_Dxy_S_bb;
  Double I_ERI_F3z_Pz_Dxy_S_bb = I_ERI_G4z_S_Dxy_S_bb+ABZ*I_ERI_F3z_S_Dxy_S_bb;
  Double I_ERI_F3x_Px_Dxz_S_bb = I_ERI_G4x_S_Dxz_S_bb+ABX*I_ERI_F3x_S_Dxz_S_bb;
  Double I_ERI_F2xy_Px_Dxz_S_bb = I_ERI_G3xy_S_Dxz_S_bb+ABX*I_ERI_F2xy_S_Dxz_S_bb;
  Double I_ERI_F2xz_Px_Dxz_S_bb = I_ERI_G3xz_S_Dxz_S_bb+ABX*I_ERI_F2xz_S_Dxz_S_bb;
  Double I_ERI_Fx2y_Px_Dxz_S_bb = I_ERI_G2x2y_S_Dxz_S_bb+ABX*I_ERI_Fx2y_S_Dxz_S_bb;
  Double I_ERI_Fxyz_Px_Dxz_S_bb = I_ERI_G2xyz_S_Dxz_S_bb+ABX*I_ERI_Fxyz_S_Dxz_S_bb;
  Double I_ERI_Fx2z_Px_Dxz_S_bb = I_ERI_G2x2z_S_Dxz_S_bb+ABX*I_ERI_Fx2z_S_Dxz_S_bb;
  Double I_ERI_F3y_Px_Dxz_S_bb = I_ERI_Gx3y_S_Dxz_S_bb+ABX*I_ERI_F3y_S_Dxz_S_bb;
  Double I_ERI_F2yz_Px_Dxz_S_bb = I_ERI_Gx2yz_S_Dxz_S_bb+ABX*I_ERI_F2yz_S_Dxz_S_bb;
  Double I_ERI_Fy2z_Px_Dxz_S_bb = I_ERI_Gxy2z_S_Dxz_S_bb+ABX*I_ERI_Fy2z_S_Dxz_S_bb;
  Double I_ERI_F3z_Px_Dxz_S_bb = I_ERI_Gx3z_S_Dxz_S_bb+ABX*I_ERI_F3z_S_Dxz_S_bb;
  Double I_ERI_F3x_Py_Dxz_S_bb = I_ERI_G3xy_S_Dxz_S_bb+ABY*I_ERI_F3x_S_Dxz_S_bb;
  Double I_ERI_F2xy_Py_Dxz_S_bb = I_ERI_G2x2y_S_Dxz_S_bb+ABY*I_ERI_F2xy_S_Dxz_S_bb;
  Double I_ERI_F2xz_Py_Dxz_S_bb = I_ERI_G2xyz_S_Dxz_S_bb+ABY*I_ERI_F2xz_S_Dxz_S_bb;
  Double I_ERI_Fx2y_Py_Dxz_S_bb = I_ERI_Gx3y_S_Dxz_S_bb+ABY*I_ERI_Fx2y_S_Dxz_S_bb;
  Double I_ERI_Fxyz_Py_Dxz_S_bb = I_ERI_Gx2yz_S_Dxz_S_bb+ABY*I_ERI_Fxyz_S_Dxz_S_bb;
  Double I_ERI_Fx2z_Py_Dxz_S_bb = I_ERI_Gxy2z_S_Dxz_S_bb+ABY*I_ERI_Fx2z_S_Dxz_S_bb;
  Double I_ERI_F3y_Py_Dxz_S_bb = I_ERI_G4y_S_Dxz_S_bb+ABY*I_ERI_F3y_S_Dxz_S_bb;
  Double I_ERI_F2yz_Py_Dxz_S_bb = I_ERI_G3yz_S_Dxz_S_bb+ABY*I_ERI_F2yz_S_Dxz_S_bb;
  Double I_ERI_Fy2z_Py_Dxz_S_bb = I_ERI_G2y2z_S_Dxz_S_bb+ABY*I_ERI_Fy2z_S_Dxz_S_bb;
  Double I_ERI_F3z_Py_Dxz_S_bb = I_ERI_Gy3z_S_Dxz_S_bb+ABY*I_ERI_F3z_S_Dxz_S_bb;
  Double I_ERI_F3x_Pz_Dxz_S_bb = I_ERI_G3xz_S_Dxz_S_bb+ABZ*I_ERI_F3x_S_Dxz_S_bb;
  Double I_ERI_F2xy_Pz_Dxz_S_bb = I_ERI_G2xyz_S_Dxz_S_bb+ABZ*I_ERI_F2xy_S_Dxz_S_bb;
  Double I_ERI_F2xz_Pz_Dxz_S_bb = I_ERI_G2x2z_S_Dxz_S_bb+ABZ*I_ERI_F2xz_S_Dxz_S_bb;
  Double I_ERI_Fx2y_Pz_Dxz_S_bb = I_ERI_Gx2yz_S_Dxz_S_bb+ABZ*I_ERI_Fx2y_S_Dxz_S_bb;
  Double I_ERI_Fxyz_Pz_Dxz_S_bb = I_ERI_Gxy2z_S_Dxz_S_bb+ABZ*I_ERI_Fxyz_S_Dxz_S_bb;
  Double I_ERI_Fx2z_Pz_Dxz_S_bb = I_ERI_Gx3z_S_Dxz_S_bb+ABZ*I_ERI_Fx2z_S_Dxz_S_bb;
  Double I_ERI_F3y_Pz_Dxz_S_bb = I_ERI_G3yz_S_Dxz_S_bb+ABZ*I_ERI_F3y_S_Dxz_S_bb;
  Double I_ERI_F2yz_Pz_Dxz_S_bb = I_ERI_G2y2z_S_Dxz_S_bb+ABZ*I_ERI_F2yz_S_Dxz_S_bb;
  Double I_ERI_Fy2z_Pz_Dxz_S_bb = I_ERI_Gy3z_S_Dxz_S_bb+ABZ*I_ERI_Fy2z_S_Dxz_S_bb;
  Double I_ERI_F3z_Pz_Dxz_S_bb = I_ERI_G4z_S_Dxz_S_bb+ABZ*I_ERI_F3z_S_Dxz_S_bb;
  Double I_ERI_F3x_Px_D2y_S_bb = I_ERI_G4x_S_D2y_S_bb+ABX*I_ERI_F3x_S_D2y_S_bb;
  Double I_ERI_F2xy_Px_D2y_S_bb = I_ERI_G3xy_S_D2y_S_bb+ABX*I_ERI_F2xy_S_D2y_S_bb;
  Double I_ERI_F2xz_Px_D2y_S_bb = I_ERI_G3xz_S_D2y_S_bb+ABX*I_ERI_F2xz_S_D2y_S_bb;
  Double I_ERI_Fx2y_Px_D2y_S_bb = I_ERI_G2x2y_S_D2y_S_bb+ABX*I_ERI_Fx2y_S_D2y_S_bb;
  Double I_ERI_Fxyz_Px_D2y_S_bb = I_ERI_G2xyz_S_D2y_S_bb+ABX*I_ERI_Fxyz_S_D2y_S_bb;
  Double I_ERI_Fx2z_Px_D2y_S_bb = I_ERI_G2x2z_S_D2y_S_bb+ABX*I_ERI_Fx2z_S_D2y_S_bb;
  Double I_ERI_F3y_Px_D2y_S_bb = I_ERI_Gx3y_S_D2y_S_bb+ABX*I_ERI_F3y_S_D2y_S_bb;
  Double I_ERI_F2yz_Px_D2y_S_bb = I_ERI_Gx2yz_S_D2y_S_bb+ABX*I_ERI_F2yz_S_D2y_S_bb;
  Double I_ERI_Fy2z_Px_D2y_S_bb = I_ERI_Gxy2z_S_D2y_S_bb+ABX*I_ERI_Fy2z_S_D2y_S_bb;
  Double I_ERI_F3z_Px_D2y_S_bb = I_ERI_Gx3z_S_D2y_S_bb+ABX*I_ERI_F3z_S_D2y_S_bb;
  Double I_ERI_F3x_Py_D2y_S_bb = I_ERI_G3xy_S_D2y_S_bb+ABY*I_ERI_F3x_S_D2y_S_bb;
  Double I_ERI_F2xy_Py_D2y_S_bb = I_ERI_G2x2y_S_D2y_S_bb+ABY*I_ERI_F2xy_S_D2y_S_bb;
  Double I_ERI_F2xz_Py_D2y_S_bb = I_ERI_G2xyz_S_D2y_S_bb+ABY*I_ERI_F2xz_S_D2y_S_bb;
  Double I_ERI_Fx2y_Py_D2y_S_bb = I_ERI_Gx3y_S_D2y_S_bb+ABY*I_ERI_Fx2y_S_D2y_S_bb;
  Double I_ERI_Fxyz_Py_D2y_S_bb = I_ERI_Gx2yz_S_D2y_S_bb+ABY*I_ERI_Fxyz_S_D2y_S_bb;
  Double I_ERI_Fx2z_Py_D2y_S_bb = I_ERI_Gxy2z_S_D2y_S_bb+ABY*I_ERI_Fx2z_S_D2y_S_bb;
  Double I_ERI_F3y_Py_D2y_S_bb = I_ERI_G4y_S_D2y_S_bb+ABY*I_ERI_F3y_S_D2y_S_bb;
  Double I_ERI_F2yz_Py_D2y_S_bb = I_ERI_G3yz_S_D2y_S_bb+ABY*I_ERI_F2yz_S_D2y_S_bb;
  Double I_ERI_Fy2z_Py_D2y_S_bb = I_ERI_G2y2z_S_D2y_S_bb+ABY*I_ERI_Fy2z_S_D2y_S_bb;
  Double I_ERI_F3z_Py_D2y_S_bb = I_ERI_Gy3z_S_D2y_S_bb+ABY*I_ERI_F3z_S_D2y_S_bb;
  Double I_ERI_F3x_Pz_D2y_S_bb = I_ERI_G3xz_S_D2y_S_bb+ABZ*I_ERI_F3x_S_D2y_S_bb;
  Double I_ERI_F2xy_Pz_D2y_S_bb = I_ERI_G2xyz_S_D2y_S_bb+ABZ*I_ERI_F2xy_S_D2y_S_bb;
  Double I_ERI_F2xz_Pz_D2y_S_bb = I_ERI_G2x2z_S_D2y_S_bb+ABZ*I_ERI_F2xz_S_D2y_S_bb;
  Double I_ERI_Fx2y_Pz_D2y_S_bb = I_ERI_Gx2yz_S_D2y_S_bb+ABZ*I_ERI_Fx2y_S_D2y_S_bb;
  Double I_ERI_Fxyz_Pz_D2y_S_bb = I_ERI_Gxy2z_S_D2y_S_bb+ABZ*I_ERI_Fxyz_S_D2y_S_bb;
  Double I_ERI_Fx2z_Pz_D2y_S_bb = I_ERI_Gx3z_S_D2y_S_bb+ABZ*I_ERI_Fx2z_S_D2y_S_bb;
  Double I_ERI_F3y_Pz_D2y_S_bb = I_ERI_G3yz_S_D2y_S_bb+ABZ*I_ERI_F3y_S_D2y_S_bb;
  Double I_ERI_F2yz_Pz_D2y_S_bb = I_ERI_G2y2z_S_D2y_S_bb+ABZ*I_ERI_F2yz_S_D2y_S_bb;
  Double I_ERI_Fy2z_Pz_D2y_S_bb = I_ERI_Gy3z_S_D2y_S_bb+ABZ*I_ERI_Fy2z_S_D2y_S_bb;
  Double I_ERI_F3z_Pz_D2y_S_bb = I_ERI_G4z_S_D2y_S_bb+ABZ*I_ERI_F3z_S_D2y_S_bb;
  Double I_ERI_F3x_Px_Dyz_S_bb = I_ERI_G4x_S_Dyz_S_bb+ABX*I_ERI_F3x_S_Dyz_S_bb;
  Double I_ERI_F2xy_Px_Dyz_S_bb = I_ERI_G3xy_S_Dyz_S_bb+ABX*I_ERI_F2xy_S_Dyz_S_bb;
  Double I_ERI_F2xz_Px_Dyz_S_bb = I_ERI_G3xz_S_Dyz_S_bb+ABX*I_ERI_F2xz_S_Dyz_S_bb;
  Double I_ERI_Fx2y_Px_Dyz_S_bb = I_ERI_G2x2y_S_Dyz_S_bb+ABX*I_ERI_Fx2y_S_Dyz_S_bb;
  Double I_ERI_Fxyz_Px_Dyz_S_bb = I_ERI_G2xyz_S_Dyz_S_bb+ABX*I_ERI_Fxyz_S_Dyz_S_bb;
  Double I_ERI_Fx2z_Px_Dyz_S_bb = I_ERI_G2x2z_S_Dyz_S_bb+ABX*I_ERI_Fx2z_S_Dyz_S_bb;
  Double I_ERI_F3y_Px_Dyz_S_bb = I_ERI_Gx3y_S_Dyz_S_bb+ABX*I_ERI_F3y_S_Dyz_S_bb;
  Double I_ERI_F2yz_Px_Dyz_S_bb = I_ERI_Gx2yz_S_Dyz_S_bb+ABX*I_ERI_F2yz_S_Dyz_S_bb;
  Double I_ERI_Fy2z_Px_Dyz_S_bb = I_ERI_Gxy2z_S_Dyz_S_bb+ABX*I_ERI_Fy2z_S_Dyz_S_bb;
  Double I_ERI_F3z_Px_Dyz_S_bb = I_ERI_Gx3z_S_Dyz_S_bb+ABX*I_ERI_F3z_S_Dyz_S_bb;
  Double I_ERI_F3x_Py_Dyz_S_bb = I_ERI_G3xy_S_Dyz_S_bb+ABY*I_ERI_F3x_S_Dyz_S_bb;
  Double I_ERI_F2xy_Py_Dyz_S_bb = I_ERI_G2x2y_S_Dyz_S_bb+ABY*I_ERI_F2xy_S_Dyz_S_bb;
  Double I_ERI_F2xz_Py_Dyz_S_bb = I_ERI_G2xyz_S_Dyz_S_bb+ABY*I_ERI_F2xz_S_Dyz_S_bb;
  Double I_ERI_Fx2y_Py_Dyz_S_bb = I_ERI_Gx3y_S_Dyz_S_bb+ABY*I_ERI_Fx2y_S_Dyz_S_bb;
  Double I_ERI_Fxyz_Py_Dyz_S_bb = I_ERI_Gx2yz_S_Dyz_S_bb+ABY*I_ERI_Fxyz_S_Dyz_S_bb;
  Double I_ERI_Fx2z_Py_Dyz_S_bb = I_ERI_Gxy2z_S_Dyz_S_bb+ABY*I_ERI_Fx2z_S_Dyz_S_bb;
  Double I_ERI_F3y_Py_Dyz_S_bb = I_ERI_G4y_S_Dyz_S_bb+ABY*I_ERI_F3y_S_Dyz_S_bb;
  Double I_ERI_F2yz_Py_Dyz_S_bb = I_ERI_G3yz_S_Dyz_S_bb+ABY*I_ERI_F2yz_S_Dyz_S_bb;
  Double I_ERI_Fy2z_Py_Dyz_S_bb = I_ERI_G2y2z_S_Dyz_S_bb+ABY*I_ERI_Fy2z_S_Dyz_S_bb;
  Double I_ERI_F3z_Py_Dyz_S_bb = I_ERI_Gy3z_S_Dyz_S_bb+ABY*I_ERI_F3z_S_Dyz_S_bb;
  Double I_ERI_F3x_Pz_Dyz_S_bb = I_ERI_G3xz_S_Dyz_S_bb+ABZ*I_ERI_F3x_S_Dyz_S_bb;
  Double I_ERI_F2xy_Pz_Dyz_S_bb = I_ERI_G2xyz_S_Dyz_S_bb+ABZ*I_ERI_F2xy_S_Dyz_S_bb;
  Double I_ERI_F2xz_Pz_Dyz_S_bb = I_ERI_G2x2z_S_Dyz_S_bb+ABZ*I_ERI_F2xz_S_Dyz_S_bb;
  Double I_ERI_Fx2y_Pz_Dyz_S_bb = I_ERI_Gx2yz_S_Dyz_S_bb+ABZ*I_ERI_Fx2y_S_Dyz_S_bb;
  Double I_ERI_Fxyz_Pz_Dyz_S_bb = I_ERI_Gxy2z_S_Dyz_S_bb+ABZ*I_ERI_Fxyz_S_Dyz_S_bb;
  Double I_ERI_Fx2z_Pz_Dyz_S_bb = I_ERI_Gx3z_S_Dyz_S_bb+ABZ*I_ERI_Fx2z_S_Dyz_S_bb;
  Double I_ERI_F3y_Pz_Dyz_S_bb = I_ERI_G3yz_S_Dyz_S_bb+ABZ*I_ERI_F3y_S_Dyz_S_bb;
  Double I_ERI_F2yz_Pz_Dyz_S_bb = I_ERI_G2y2z_S_Dyz_S_bb+ABZ*I_ERI_F2yz_S_Dyz_S_bb;
  Double I_ERI_Fy2z_Pz_Dyz_S_bb = I_ERI_Gy3z_S_Dyz_S_bb+ABZ*I_ERI_Fy2z_S_Dyz_S_bb;
  Double I_ERI_F3z_Pz_Dyz_S_bb = I_ERI_G4z_S_Dyz_S_bb+ABZ*I_ERI_F3z_S_Dyz_S_bb;
  Double I_ERI_F3x_Px_D2z_S_bb = I_ERI_G4x_S_D2z_S_bb+ABX*I_ERI_F3x_S_D2z_S_bb;
  Double I_ERI_F2xy_Px_D2z_S_bb = I_ERI_G3xy_S_D2z_S_bb+ABX*I_ERI_F2xy_S_D2z_S_bb;
  Double I_ERI_F2xz_Px_D2z_S_bb = I_ERI_G3xz_S_D2z_S_bb+ABX*I_ERI_F2xz_S_D2z_S_bb;
  Double I_ERI_Fx2y_Px_D2z_S_bb = I_ERI_G2x2y_S_D2z_S_bb+ABX*I_ERI_Fx2y_S_D2z_S_bb;
  Double I_ERI_Fxyz_Px_D2z_S_bb = I_ERI_G2xyz_S_D2z_S_bb+ABX*I_ERI_Fxyz_S_D2z_S_bb;
  Double I_ERI_Fx2z_Px_D2z_S_bb = I_ERI_G2x2z_S_D2z_S_bb+ABX*I_ERI_Fx2z_S_D2z_S_bb;
  Double I_ERI_F3y_Px_D2z_S_bb = I_ERI_Gx3y_S_D2z_S_bb+ABX*I_ERI_F3y_S_D2z_S_bb;
  Double I_ERI_F2yz_Px_D2z_S_bb = I_ERI_Gx2yz_S_D2z_S_bb+ABX*I_ERI_F2yz_S_D2z_S_bb;
  Double I_ERI_Fy2z_Px_D2z_S_bb = I_ERI_Gxy2z_S_D2z_S_bb+ABX*I_ERI_Fy2z_S_D2z_S_bb;
  Double I_ERI_F3z_Px_D2z_S_bb = I_ERI_Gx3z_S_D2z_S_bb+ABX*I_ERI_F3z_S_D2z_S_bb;
  Double I_ERI_F3x_Py_D2z_S_bb = I_ERI_G3xy_S_D2z_S_bb+ABY*I_ERI_F3x_S_D2z_S_bb;
  Double I_ERI_F2xy_Py_D2z_S_bb = I_ERI_G2x2y_S_D2z_S_bb+ABY*I_ERI_F2xy_S_D2z_S_bb;
  Double I_ERI_F2xz_Py_D2z_S_bb = I_ERI_G2xyz_S_D2z_S_bb+ABY*I_ERI_F2xz_S_D2z_S_bb;
  Double I_ERI_Fx2y_Py_D2z_S_bb = I_ERI_Gx3y_S_D2z_S_bb+ABY*I_ERI_Fx2y_S_D2z_S_bb;
  Double I_ERI_Fxyz_Py_D2z_S_bb = I_ERI_Gx2yz_S_D2z_S_bb+ABY*I_ERI_Fxyz_S_D2z_S_bb;
  Double I_ERI_Fx2z_Py_D2z_S_bb = I_ERI_Gxy2z_S_D2z_S_bb+ABY*I_ERI_Fx2z_S_D2z_S_bb;
  Double I_ERI_F3y_Py_D2z_S_bb = I_ERI_G4y_S_D2z_S_bb+ABY*I_ERI_F3y_S_D2z_S_bb;
  Double I_ERI_F2yz_Py_D2z_S_bb = I_ERI_G3yz_S_D2z_S_bb+ABY*I_ERI_F2yz_S_D2z_S_bb;
  Double I_ERI_Fy2z_Py_D2z_S_bb = I_ERI_G2y2z_S_D2z_S_bb+ABY*I_ERI_Fy2z_S_D2z_S_bb;
  Double I_ERI_F3z_Py_D2z_S_bb = I_ERI_Gy3z_S_D2z_S_bb+ABY*I_ERI_F3z_S_D2z_S_bb;
  Double I_ERI_F3x_Pz_D2z_S_bb = I_ERI_G3xz_S_D2z_S_bb+ABZ*I_ERI_F3x_S_D2z_S_bb;
  Double I_ERI_F2xy_Pz_D2z_S_bb = I_ERI_G2xyz_S_D2z_S_bb+ABZ*I_ERI_F2xy_S_D2z_S_bb;
  Double I_ERI_F2xz_Pz_D2z_S_bb = I_ERI_G2x2z_S_D2z_S_bb+ABZ*I_ERI_F2xz_S_D2z_S_bb;
  Double I_ERI_Fx2y_Pz_D2z_S_bb = I_ERI_Gx2yz_S_D2z_S_bb+ABZ*I_ERI_Fx2y_S_D2z_S_bb;
  Double I_ERI_Fxyz_Pz_D2z_S_bb = I_ERI_Gxy2z_S_D2z_S_bb+ABZ*I_ERI_Fxyz_S_D2z_S_bb;
  Double I_ERI_Fx2z_Pz_D2z_S_bb = I_ERI_Gx3z_S_D2z_S_bb+ABZ*I_ERI_Fx2z_S_D2z_S_bb;
  Double I_ERI_F3y_Pz_D2z_S_bb = I_ERI_G3yz_S_D2z_S_bb+ABZ*I_ERI_F3y_S_D2z_S_bb;
  Double I_ERI_F2yz_Pz_D2z_S_bb = I_ERI_G2y2z_S_D2z_S_bb+ABZ*I_ERI_F2yz_S_D2z_S_bb;
  Double I_ERI_Fy2z_Pz_D2z_S_bb = I_ERI_Gy3z_S_D2z_S_bb+ABZ*I_ERI_Fy2z_S_D2z_S_bb;
  Double I_ERI_F3z_Pz_D2z_S_bb = I_ERI_G4z_S_D2z_S_bb+ABZ*I_ERI_F3z_S_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_D_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_D_S_bb
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bb
   ************************************************************/
  Double I_ERI_G4x_Px_D2x_S_bb = I_ERI_H5x_S_D2x_S_bb+ABX*I_ERI_G4x_S_D2x_S_bb;
  Double I_ERI_G3xy_Px_D2x_S_bb = I_ERI_H4xy_S_D2x_S_bb+ABX*I_ERI_G3xy_S_D2x_S_bb;
  Double I_ERI_G3xz_Px_D2x_S_bb = I_ERI_H4xz_S_D2x_S_bb+ABX*I_ERI_G3xz_S_D2x_S_bb;
  Double I_ERI_G2x2y_Px_D2x_S_bb = I_ERI_H3x2y_S_D2x_S_bb+ABX*I_ERI_G2x2y_S_D2x_S_bb;
  Double I_ERI_G2xyz_Px_D2x_S_bb = I_ERI_H3xyz_S_D2x_S_bb+ABX*I_ERI_G2xyz_S_D2x_S_bb;
  Double I_ERI_G2x2z_Px_D2x_S_bb = I_ERI_H3x2z_S_D2x_S_bb+ABX*I_ERI_G2x2z_S_D2x_S_bb;
  Double I_ERI_Gx3y_Px_D2x_S_bb = I_ERI_H2x3y_S_D2x_S_bb+ABX*I_ERI_Gx3y_S_D2x_S_bb;
  Double I_ERI_Gx2yz_Px_D2x_S_bb = I_ERI_H2x2yz_S_D2x_S_bb+ABX*I_ERI_Gx2yz_S_D2x_S_bb;
  Double I_ERI_Gxy2z_Px_D2x_S_bb = I_ERI_H2xy2z_S_D2x_S_bb+ABX*I_ERI_Gxy2z_S_D2x_S_bb;
  Double I_ERI_Gx3z_Px_D2x_S_bb = I_ERI_H2x3z_S_D2x_S_bb+ABX*I_ERI_Gx3z_S_D2x_S_bb;
  Double I_ERI_G4y_Px_D2x_S_bb = I_ERI_Hx4y_S_D2x_S_bb+ABX*I_ERI_G4y_S_D2x_S_bb;
  Double I_ERI_G3yz_Px_D2x_S_bb = I_ERI_Hx3yz_S_D2x_S_bb+ABX*I_ERI_G3yz_S_D2x_S_bb;
  Double I_ERI_G2y2z_Px_D2x_S_bb = I_ERI_Hx2y2z_S_D2x_S_bb+ABX*I_ERI_G2y2z_S_D2x_S_bb;
  Double I_ERI_Gy3z_Px_D2x_S_bb = I_ERI_Hxy3z_S_D2x_S_bb+ABX*I_ERI_Gy3z_S_D2x_S_bb;
  Double I_ERI_G4z_Px_D2x_S_bb = I_ERI_Hx4z_S_D2x_S_bb+ABX*I_ERI_G4z_S_D2x_S_bb;
  Double I_ERI_G3xy_Py_D2x_S_bb = I_ERI_H3x2y_S_D2x_S_bb+ABY*I_ERI_G3xy_S_D2x_S_bb;
  Double I_ERI_G3xz_Py_D2x_S_bb = I_ERI_H3xyz_S_D2x_S_bb+ABY*I_ERI_G3xz_S_D2x_S_bb;
  Double I_ERI_G2x2y_Py_D2x_S_bb = I_ERI_H2x3y_S_D2x_S_bb+ABY*I_ERI_G2x2y_S_D2x_S_bb;
  Double I_ERI_G2xyz_Py_D2x_S_bb = I_ERI_H2x2yz_S_D2x_S_bb+ABY*I_ERI_G2xyz_S_D2x_S_bb;
  Double I_ERI_G2x2z_Py_D2x_S_bb = I_ERI_H2xy2z_S_D2x_S_bb+ABY*I_ERI_G2x2z_S_D2x_S_bb;
  Double I_ERI_Gx3y_Py_D2x_S_bb = I_ERI_Hx4y_S_D2x_S_bb+ABY*I_ERI_Gx3y_S_D2x_S_bb;
  Double I_ERI_Gx2yz_Py_D2x_S_bb = I_ERI_Hx3yz_S_D2x_S_bb+ABY*I_ERI_Gx2yz_S_D2x_S_bb;
  Double I_ERI_Gxy2z_Py_D2x_S_bb = I_ERI_Hx2y2z_S_D2x_S_bb+ABY*I_ERI_Gxy2z_S_D2x_S_bb;
  Double I_ERI_Gx3z_Py_D2x_S_bb = I_ERI_Hxy3z_S_D2x_S_bb+ABY*I_ERI_Gx3z_S_D2x_S_bb;
  Double I_ERI_G4y_Py_D2x_S_bb = I_ERI_H5y_S_D2x_S_bb+ABY*I_ERI_G4y_S_D2x_S_bb;
  Double I_ERI_G3yz_Py_D2x_S_bb = I_ERI_H4yz_S_D2x_S_bb+ABY*I_ERI_G3yz_S_D2x_S_bb;
  Double I_ERI_G2y2z_Py_D2x_S_bb = I_ERI_H3y2z_S_D2x_S_bb+ABY*I_ERI_G2y2z_S_D2x_S_bb;
  Double I_ERI_Gy3z_Py_D2x_S_bb = I_ERI_H2y3z_S_D2x_S_bb+ABY*I_ERI_Gy3z_S_D2x_S_bb;
  Double I_ERI_G4z_Py_D2x_S_bb = I_ERI_Hy4z_S_D2x_S_bb+ABY*I_ERI_G4z_S_D2x_S_bb;
  Double I_ERI_G3xz_Pz_D2x_S_bb = I_ERI_H3x2z_S_D2x_S_bb+ABZ*I_ERI_G3xz_S_D2x_S_bb;
  Double I_ERI_G2xyz_Pz_D2x_S_bb = I_ERI_H2xy2z_S_D2x_S_bb+ABZ*I_ERI_G2xyz_S_D2x_S_bb;
  Double I_ERI_G2x2z_Pz_D2x_S_bb = I_ERI_H2x3z_S_D2x_S_bb+ABZ*I_ERI_G2x2z_S_D2x_S_bb;
  Double I_ERI_Gx2yz_Pz_D2x_S_bb = I_ERI_Hx2y2z_S_D2x_S_bb+ABZ*I_ERI_Gx2yz_S_D2x_S_bb;
  Double I_ERI_Gxy2z_Pz_D2x_S_bb = I_ERI_Hxy3z_S_D2x_S_bb+ABZ*I_ERI_Gxy2z_S_D2x_S_bb;
  Double I_ERI_Gx3z_Pz_D2x_S_bb = I_ERI_Hx4z_S_D2x_S_bb+ABZ*I_ERI_Gx3z_S_D2x_S_bb;
  Double I_ERI_G3yz_Pz_D2x_S_bb = I_ERI_H3y2z_S_D2x_S_bb+ABZ*I_ERI_G3yz_S_D2x_S_bb;
  Double I_ERI_G2y2z_Pz_D2x_S_bb = I_ERI_H2y3z_S_D2x_S_bb+ABZ*I_ERI_G2y2z_S_D2x_S_bb;
  Double I_ERI_Gy3z_Pz_D2x_S_bb = I_ERI_Hy4z_S_D2x_S_bb+ABZ*I_ERI_Gy3z_S_D2x_S_bb;
  Double I_ERI_G4z_Pz_D2x_S_bb = I_ERI_H5z_S_D2x_S_bb+ABZ*I_ERI_G4z_S_D2x_S_bb;
  Double I_ERI_G4x_Px_Dxy_S_bb = I_ERI_H5x_S_Dxy_S_bb+ABX*I_ERI_G4x_S_Dxy_S_bb;
  Double I_ERI_G3xy_Px_Dxy_S_bb = I_ERI_H4xy_S_Dxy_S_bb+ABX*I_ERI_G3xy_S_Dxy_S_bb;
  Double I_ERI_G3xz_Px_Dxy_S_bb = I_ERI_H4xz_S_Dxy_S_bb+ABX*I_ERI_G3xz_S_Dxy_S_bb;
  Double I_ERI_G2x2y_Px_Dxy_S_bb = I_ERI_H3x2y_S_Dxy_S_bb+ABX*I_ERI_G2x2y_S_Dxy_S_bb;
  Double I_ERI_G2xyz_Px_Dxy_S_bb = I_ERI_H3xyz_S_Dxy_S_bb+ABX*I_ERI_G2xyz_S_Dxy_S_bb;
  Double I_ERI_G2x2z_Px_Dxy_S_bb = I_ERI_H3x2z_S_Dxy_S_bb+ABX*I_ERI_G2x2z_S_Dxy_S_bb;
  Double I_ERI_Gx3y_Px_Dxy_S_bb = I_ERI_H2x3y_S_Dxy_S_bb+ABX*I_ERI_Gx3y_S_Dxy_S_bb;
  Double I_ERI_Gx2yz_Px_Dxy_S_bb = I_ERI_H2x2yz_S_Dxy_S_bb+ABX*I_ERI_Gx2yz_S_Dxy_S_bb;
  Double I_ERI_Gxy2z_Px_Dxy_S_bb = I_ERI_H2xy2z_S_Dxy_S_bb+ABX*I_ERI_Gxy2z_S_Dxy_S_bb;
  Double I_ERI_Gx3z_Px_Dxy_S_bb = I_ERI_H2x3z_S_Dxy_S_bb+ABX*I_ERI_Gx3z_S_Dxy_S_bb;
  Double I_ERI_G4y_Px_Dxy_S_bb = I_ERI_Hx4y_S_Dxy_S_bb+ABX*I_ERI_G4y_S_Dxy_S_bb;
  Double I_ERI_G3yz_Px_Dxy_S_bb = I_ERI_Hx3yz_S_Dxy_S_bb+ABX*I_ERI_G3yz_S_Dxy_S_bb;
  Double I_ERI_G2y2z_Px_Dxy_S_bb = I_ERI_Hx2y2z_S_Dxy_S_bb+ABX*I_ERI_G2y2z_S_Dxy_S_bb;
  Double I_ERI_Gy3z_Px_Dxy_S_bb = I_ERI_Hxy3z_S_Dxy_S_bb+ABX*I_ERI_Gy3z_S_Dxy_S_bb;
  Double I_ERI_G4z_Px_Dxy_S_bb = I_ERI_Hx4z_S_Dxy_S_bb+ABX*I_ERI_G4z_S_Dxy_S_bb;
  Double I_ERI_G3xy_Py_Dxy_S_bb = I_ERI_H3x2y_S_Dxy_S_bb+ABY*I_ERI_G3xy_S_Dxy_S_bb;
  Double I_ERI_G3xz_Py_Dxy_S_bb = I_ERI_H3xyz_S_Dxy_S_bb+ABY*I_ERI_G3xz_S_Dxy_S_bb;
  Double I_ERI_G2x2y_Py_Dxy_S_bb = I_ERI_H2x3y_S_Dxy_S_bb+ABY*I_ERI_G2x2y_S_Dxy_S_bb;
  Double I_ERI_G2xyz_Py_Dxy_S_bb = I_ERI_H2x2yz_S_Dxy_S_bb+ABY*I_ERI_G2xyz_S_Dxy_S_bb;
  Double I_ERI_G2x2z_Py_Dxy_S_bb = I_ERI_H2xy2z_S_Dxy_S_bb+ABY*I_ERI_G2x2z_S_Dxy_S_bb;
  Double I_ERI_Gx3y_Py_Dxy_S_bb = I_ERI_Hx4y_S_Dxy_S_bb+ABY*I_ERI_Gx3y_S_Dxy_S_bb;
  Double I_ERI_Gx2yz_Py_Dxy_S_bb = I_ERI_Hx3yz_S_Dxy_S_bb+ABY*I_ERI_Gx2yz_S_Dxy_S_bb;
  Double I_ERI_Gxy2z_Py_Dxy_S_bb = I_ERI_Hx2y2z_S_Dxy_S_bb+ABY*I_ERI_Gxy2z_S_Dxy_S_bb;
  Double I_ERI_Gx3z_Py_Dxy_S_bb = I_ERI_Hxy3z_S_Dxy_S_bb+ABY*I_ERI_Gx3z_S_Dxy_S_bb;
  Double I_ERI_G4y_Py_Dxy_S_bb = I_ERI_H5y_S_Dxy_S_bb+ABY*I_ERI_G4y_S_Dxy_S_bb;
  Double I_ERI_G3yz_Py_Dxy_S_bb = I_ERI_H4yz_S_Dxy_S_bb+ABY*I_ERI_G3yz_S_Dxy_S_bb;
  Double I_ERI_G2y2z_Py_Dxy_S_bb = I_ERI_H3y2z_S_Dxy_S_bb+ABY*I_ERI_G2y2z_S_Dxy_S_bb;
  Double I_ERI_Gy3z_Py_Dxy_S_bb = I_ERI_H2y3z_S_Dxy_S_bb+ABY*I_ERI_Gy3z_S_Dxy_S_bb;
  Double I_ERI_G4z_Py_Dxy_S_bb = I_ERI_Hy4z_S_Dxy_S_bb+ABY*I_ERI_G4z_S_Dxy_S_bb;
  Double I_ERI_G3xz_Pz_Dxy_S_bb = I_ERI_H3x2z_S_Dxy_S_bb+ABZ*I_ERI_G3xz_S_Dxy_S_bb;
  Double I_ERI_G2xyz_Pz_Dxy_S_bb = I_ERI_H2xy2z_S_Dxy_S_bb+ABZ*I_ERI_G2xyz_S_Dxy_S_bb;
  Double I_ERI_G2x2z_Pz_Dxy_S_bb = I_ERI_H2x3z_S_Dxy_S_bb+ABZ*I_ERI_G2x2z_S_Dxy_S_bb;
  Double I_ERI_Gx2yz_Pz_Dxy_S_bb = I_ERI_Hx2y2z_S_Dxy_S_bb+ABZ*I_ERI_Gx2yz_S_Dxy_S_bb;
  Double I_ERI_Gxy2z_Pz_Dxy_S_bb = I_ERI_Hxy3z_S_Dxy_S_bb+ABZ*I_ERI_Gxy2z_S_Dxy_S_bb;
  Double I_ERI_Gx3z_Pz_Dxy_S_bb = I_ERI_Hx4z_S_Dxy_S_bb+ABZ*I_ERI_Gx3z_S_Dxy_S_bb;
  Double I_ERI_G3yz_Pz_Dxy_S_bb = I_ERI_H3y2z_S_Dxy_S_bb+ABZ*I_ERI_G3yz_S_Dxy_S_bb;
  Double I_ERI_G2y2z_Pz_Dxy_S_bb = I_ERI_H2y3z_S_Dxy_S_bb+ABZ*I_ERI_G2y2z_S_Dxy_S_bb;
  Double I_ERI_Gy3z_Pz_Dxy_S_bb = I_ERI_Hy4z_S_Dxy_S_bb+ABZ*I_ERI_Gy3z_S_Dxy_S_bb;
  Double I_ERI_G4z_Pz_Dxy_S_bb = I_ERI_H5z_S_Dxy_S_bb+ABZ*I_ERI_G4z_S_Dxy_S_bb;
  Double I_ERI_G4x_Px_Dxz_S_bb = I_ERI_H5x_S_Dxz_S_bb+ABX*I_ERI_G4x_S_Dxz_S_bb;
  Double I_ERI_G3xy_Px_Dxz_S_bb = I_ERI_H4xy_S_Dxz_S_bb+ABX*I_ERI_G3xy_S_Dxz_S_bb;
  Double I_ERI_G3xz_Px_Dxz_S_bb = I_ERI_H4xz_S_Dxz_S_bb+ABX*I_ERI_G3xz_S_Dxz_S_bb;
  Double I_ERI_G2x2y_Px_Dxz_S_bb = I_ERI_H3x2y_S_Dxz_S_bb+ABX*I_ERI_G2x2y_S_Dxz_S_bb;
  Double I_ERI_G2xyz_Px_Dxz_S_bb = I_ERI_H3xyz_S_Dxz_S_bb+ABX*I_ERI_G2xyz_S_Dxz_S_bb;
  Double I_ERI_G2x2z_Px_Dxz_S_bb = I_ERI_H3x2z_S_Dxz_S_bb+ABX*I_ERI_G2x2z_S_Dxz_S_bb;
  Double I_ERI_Gx3y_Px_Dxz_S_bb = I_ERI_H2x3y_S_Dxz_S_bb+ABX*I_ERI_Gx3y_S_Dxz_S_bb;
  Double I_ERI_Gx2yz_Px_Dxz_S_bb = I_ERI_H2x2yz_S_Dxz_S_bb+ABX*I_ERI_Gx2yz_S_Dxz_S_bb;
  Double I_ERI_Gxy2z_Px_Dxz_S_bb = I_ERI_H2xy2z_S_Dxz_S_bb+ABX*I_ERI_Gxy2z_S_Dxz_S_bb;
  Double I_ERI_Gx3z_Px_Dxz_S_bb = I_ERI_H2x3z_S_Dxz_S_bb+ABX*I_ERI_Gx3z_S_Dxz_S_bb;
  Double I_ERI_G4y_Px_Dxz_S_bb = I_ERI_Hx4y_S_Dxz_S_bb+ABX*I_ERI_G4y_S_Dxz_S_bb;
  Double I_ERI_G3yz_Px_Dxz_S_bb = I_ERI_Hx3yz_S_Dxz_S_bb+ABX*I_ERI_G3yz_S_Dxz_S_bb;
  Double I_ERI_G2y2z_Px_Dxz_S_bb = I_ERI_Hx2y2z_S_Dxz_S_bb+ABX*I_ERI_G2y2z_S_Dxz_S_bb;
  Double I_ERI_Gy3z_Px_Dxz_S_bb = I_ERI_Hxy3z_S_Dxz_S_bb+ABX*I_ERI_Gy3z_S_Dxz_S_bb;
  Double I_ERI_G4z_Px_Dxz_S_bb = I_ERI_Hx4z_S_Dxz_S_bb+ABX*I_ERI_G4z_S_Dxz_S_bb;
  Double I_ERI_G3xy_Py_Dxz_S_bb = I_ERI_H3x2y_S_Dxz_S_bb+ABY*I_ERI_G3xy_S_Dxz_S_bb;
  Double I_ERI_G3xz_Py_Dxz_S_bb = I_ERI_H3xyz_S_Dxz_S_bb+ABY*I_ERI_G3xz_S_Dxz_S_bb;
  Double I_ERI_G2x2y_Py_Dxz_S_bb = I_ERI_H2x3y_S_Dxz_S_bb+ABY*I_ERI_G2x2y_S_Dxz_S_bb;
  Double I_ERI_G2xyz_Py_Dxz_S_bb = I_ERI_H2x2yz_S_Dxz_S_bb+ABY*I_ERI_G2xyz_S_Dxz_S_bb;
  Double I_ERI_G2x2z_Py_Dxz_S_bb = I_ERI_H2xy2z_S_Dxz_S_bb+ABY*I_ERI_G2x2z_S_Dxz_S_bb;
  Double I_ERI_Gx3y_Py_Dxz_S_bb = I_ERI_Hx4y_S_Dxz_S_bb+ABY*I_ERI_Gx3y_S_Dxz_S_bb;
  Double I_ERI_Gx2yz_Py_Dxz_S_bb = I_ERI_Hx3yz_S_Dxz_S_bb+ABY*I_ERI_Gx2yz_S_Dxz_S_bb;
  Double I_ERI_Gxy2z_Py_Dxz_S_bb = I_ERI_Hx2y2z_S_Dxz_S_bb+ABY*I_ERI_Gxy2z_S_Dxz_S_bb;
  Double I_ERI_Gx3z_Py_Dxz_S_bb = I_ERI_Hxy3z_S_Dxz_S_bb+ABY*I_ERI_Gx3z_S_Dxz_S_bb;
  Double I_ERI_G4y_Py_Dxz_S_bb = I_ERI_H5y_S_Dxz_S_bb+ABY*I_ERI_G4y_S_Dxz_S_bb;
  Double I_ERI_G3yz_Py_Dxz_S_bb = I_ERI_H4yz_S_Dxz_S_bb+ABY*I_ERI_G3yz_S_Dxz_S_bb;
  Double I_ERI_G2y2z_Py_Dxz_S_bb = I_ERI_H3y2z_S_Dxz_S_bb+ABY*I_ERI_G2y2z_S_Dxz_S_bb;
  Double I_ERI_Gy3z_Py_Dxz_S_bb = I_ERI_H2y3z_S_Dxz_S_bb+ABY*I_ERI_Gy3z_S_Dxz_S_bb;
  Double I_ERI_G4z_Py_Dxz_S_bb = I_ERI_Hy4z_S_Dxz_S_bb+ABY*I_ERI_G4z_S_Dxz_S_bb;
  Double I_ERI_G3xz_Pz_Dxz_S_bb = I_ERI_H3x2z_S_Dxz_S_bb+ABZ*I_ERI_G3xz_S_Dxz_S_bb;
  Double I_ERI_G2xyz_Pz_Dxz_S_bb = I_ERI_H2xy2z_S_Dxz_S_bb+ABZ*I_ERI_G2xyz_S_Dxz_S_bb;
  Double I_ERI_G2x2z_Pz_Dxz_S_bb = I_ERI_H2x3z_S_Dxz_S_bb+ABZ*I_ERI_G2x2z_S_Dxz_S_bb;
  Double I_ERI_Gx2yz_Pz_Dxz_S_bb = I_ERI_Hx2y2z_S_Dxz_S_bb+ABZ*I_ERI_Gx2yz_S_Dxz_S_bb;
  Double I_ERI_Gxy2z_Pz_Dxz_S_bb = I_ERI_Hxy3z_S_Dxz_S_bb+ABZ*I_ERI_Gxy2z_S_Dxz_S_bb;
  Double I_ERI_Gx3z_Pz_Dxz_S_bb = I_ERI_Hx4z_S_Dxz_S_bb+ABZ*I_ERI_Gx3z_S_Dxz_S_bb;
  Double I_ERI_G3yz_Pz_Dxz_S_bb = I_ERI_H3y2z_S_Dxz_S_bb+ABZ*I_ERI_G3yz_S_Dxz_S_bb;
  Double I_ERI_G2y2z_Pz_Dxz_S_bb = I_ERI_H2y3z_S_Dxz_S_bb+ABZ*I_ERI_G2y2z_S_Dxz_S_bb;
  Double I_ERI_Gy3z_Pz_Dxz_S_bb = I_ERI_Hy4z_S_Dxz_S_bb+ABZ*I_ERI_Gy3z_S_Dxz_S_bb;
  Double I_ERI_G4z_Pz_Dxz_S_bb = I_ERI_H5z_S_Dxz_S_bb+ABZ*I_ERI_G4z_S_Dxz_S_bb;
  Double I_ERI_G4x_Px_D2y_S_bb = I_ERI_H5x_S_D2y_S_bb+ABX*I_ERI_G4x_S_D2y_S_bb;
  Double I_ERI_G3xy_Px_D2y_S_bb = I_ERI_H4xy_S_D2y_S_bb+ABX*I_ERI_G3xy_S_D2y_S_bb;
  Double I_ERI_G3xz_Px_D2y_S_bb = I_ERI_H4xz_S_D2y_S_bb+ABX*I_ERI_G3xz_S_D2y_S_bb;
  Double I_ERI_G2x2y_Px_D2y_S_bb = I_ERI_H3x2y_S_D2y_S_bb+ABX*I_ERI_G2x2y_S_D2y_S_bb;
  Double I_ERI_G2xyz_Px_D2y_S_bb = I_ERI_H3xyz_S_D2y_S_bb+ABX*I_ERI_G2xyz_S_D2y_S_bb;
  Double I_ERI_G2x2z_Px_D2y_S_bb = I_ERI_H3x2z_S_D2y_S_bb+ABX*I_ERI_G2x2z_S_D2y_S_bb;
  Double I_ERI_Gx3y_Px_D2y_S_bb = I_ERI_H2x3y_S_D2y_S_bb+ABX*I_ERI_Gx3y_S_D2y_S_bb;
  Double I_ERI_Gx2yz_Px_D2y_S_bb = I_ERI_H2x2yz_S_D2y_S_bb+ABX*I_ERI_Gx2yz_S_D2y_S_bb;
  Double I_ERI_Gxy2z_Px_D2y_S_bb = I_ERI_H2xy2z_S_D2y_S_bb+ABX*I_ERI_Gxy2z_S_D2y_S_bb;
  Double I_ERI_Gx3z_Px_D2y_S_bb = I_ERI_H2x3z_S_D2y_S_bb+ABX*I_ERI_Gx3z_S_D2y_S_bb;
  Double I_ERI_G4y_Px_D2y_S_bb = I_ERI_Hx4y_S_D2y_S_bb+ABX*I_ERI_G4y_S_D2y_S_bb;
  Double I_ERI_G3yz_Px_D2y_S_bb = I_ERI_Hx3yz_S_D2y_S_bb+ABX*I_ERI_G3yz_S_D2y_S_bb;
  Double I_ERI_G2y2z_Px_D2y_S_bb = I_ERI_Hx2y2z_S_D2y_S_bb+ABX*I_ERI_G2y2z_S_D2y_S_bb;
  Double I_ERI_Gy3z_Px_D2y_S_bb = I_ERI_Hxy3z_S_D2y_S_bb+ABX*I_ERI_Gy3z_S_D2y_S_bb;
  Double I_ERI_G4z_Px_D2y_S_bb = I_ERI_Hx4z_S_D2y_S_bb+ABX*I_ERI_G4z_S_D2y_S_bb;
  Double I_ERI_G3xy_Py_D2y_S_bb = I_ERI_H3x2y_S_D2y_S_bb+ABY*I_ERI_G3xy_S_D2y_S_bb;
  Double I_ERI_G3xz_Py_D2y_S_bb = I_ERI_H3xyz_S_D2y_S_bb+ABY*I_ERI_G3xz_S_D2y_S_bb;
  Double I_ERI_G2x2y_Py_D2y_S_bb = I_ERI_H2x3y_S_D2y_S_bb+ABY*I_ERI_G2x2y_S_D2y_S_bb;
  Double I_ERI_G2xyz_Py_D2y_S_bb = I_ERI_H2x2yz_S_D2y_S_bb+ABY*I_ERI_G2xyz_S_D2y_S_bb;
  Double I_ERI_G2x2z_Py_D2y_S_bb = I_ERI_H2xy2z_S_D2y_S_bb+ABY*I_ERI_G2x2z_S_D2y_S_bb;
  Double I_ERI_Gx3y_Py_D2y_S_bb = I_ERI_Hx4y_S_D2y_S_bb+ABY*I_ERI_Gx3y_S_D2y_S_bb;
  Double I_ERI_Gx2yz_Py_D2y_S_bb = I_ERI_Hx3yz_S_D2y_S_bb+ABY*I_ERI_Gx2yz_S_D2y_S_bb;
  Double I_ERI_Gxy2z_Py_D2y_S_bb = I_ERI_Hx2y2z_S_D2y_S_bb+ABY*I_ERI_Gxy2z_S_D2y_S_bb;
  Double I_ERI_Gx3z_Py_D2y_S_bb = I_ERI_Hxy3z_S_D2y_S_bb+ABY*I_ERI_Gx3z_S_D2y_S_bb;
  Double I_ERI_G4y_Py_D2y_S_bb = I_ERI_H5y_S_D2y_S_bb+ABY*I_ERI_G4y_S_D2y_S_bb;
  Double I_ERI_G3yz_Py_D2y_S_bb = I_ERI_H4yz_S_D2y_S_bb+ABY*I_ERI_G3yz_S_D2y_S_bb;
  Double I_ERI_G2y2z_Py_D2y_S_bb = I_ERI_H3y2z_S_D2y_S_bb+ABY*I_ERI_G2y2z_S_D2y_S_bb;
  Double I_ERI_Gy3z_Py_D2y_S_bb = I_ERI_H2y3z_S_D2y_S_bb+ABY*I_ERI_Gy3z_S_D2y_S_bb;
  Double I_ERI_G4z_Py_D2y_S_bb = I_ERI_Hy4z_S_D2y_S_bb+ABY*I_ERI_G4z_S_D2y_S_bb;
  Double I_ERI_G3xz_Pz_D2y_S_bb = I_ERI_H3x2z_S_D2y_S_bb+ABZ*I_ERI_G3xz_S_D2y_S_bb;
  Double I_ERI_G2xyz_Pz_D2y_S_bb = I_ERI_H2xy2z_S_D2y_S_bb+ABZ*I_ERI_G2xyz_S_D2y_S_bb;
  Double I_ERI_G2x2z_Pz_D2y_S_bb = I_ERI_H2x3z_S_D2y_S_bb+ABZ*I_ERI_G2x2z_S_D2y_S_bb;
  Double I_ERI_Gx2yz_Pz_D2y_S_bb = I_ERI_Hx2y2z_S_D2y_S_bb+ABZ*I_ERI_Gx2yz_S_D2y_S_bb;
  Double I_ERI_Gxy2z_Pz_D2y_S_bb = I_ERI_Hxy3z_S_D2y_S_bb+ABZ*I_ERI_Gxy2z_S_D2y_S_bb;
  Double I_ERI_Gx3z_Pz_D2y_S_bb = I_ERI_Hx4z_S_D2y_S_bb+ABZ*I_ERI_Gx3z_S_D2y_S_bb;
  Double I_ERI_G3yz_Pz_D2y_S_bb = I_ERI_H3y2z_S_D2y_S_bb+ABZ*I_ERI_G3yz_S_D2y_S_bb;
  Double I_ERI_G2y2z_Pz_D2y_S_bb = I_ERI_H2y3z_S_D2y_S_bb+ABZ*I_ERI_G2y2z_S_D2y_S_bb;
  Double I_ERI_Gy3z_Pz_D2y_S_bb = I_ERI_Hy4z_S_D2y_S_bb+ABZ*I_ERI_Gy3z_S_D2y_S_bb;
  Double I_ERI_G4z_Pz_D2y_S_bb = I_ERI_H5z_S_D2y_S_bb+ABZ*I_ERI_G4z_S_D2y_S_bb;
  Double I_ERI_G4x_Px_Dyz_S_bb = I_ERI_H5x_S_Dyz_S_bb+ABX*I_ERI_G4x_S_Dyz_S_bb;
  Double I_ERI_G3xy_Px_Dyz_S_bb = I_ERI_H4xy_S_Dyz_S_bb+ABX*I_ERI_G3xy_S_Dyz_S_bb;
  Double I_ERI_G3xz_Px_Dyz_S_bb = I_ERI_H4xz_S_Dyz_S_bb+ABX*I_ERI_G3xz_S_Dyz_S_bb;
  Double I_ERI_G2x2y_Px_Dyz_S_bb = I_ERI_H3x2y_S_Dyz_S_bb+ABX*I_ERI_G2x2y_S_Dyz_S_bb;
  Double I_ERI_G2xyz_Px_Dyz_S_bb = I_ERI_H3xyz_S_Dyz_S_bb+ABX*I_ERI_G2xyz_S_Dyz_S_bb;
  Double I_ERI_G2x2z_Px_Dyz_S_bb = I_ERI_H3x2z_S_Dyz_S_bb+ABX*I_ERI_G2x2z_S_Dyz_S_bb;
  Double I_ERI_Gx3y_Px_Dyz_S_bb = I_ERI_H2x3y_S_Dyz_S_bb+ABX*I_ERI_Gx3y_S_Dyz_S_bb;
  Double I_ERI_Gx2yz_Px_Dyz_S_bb = I_ERI_H2x2yz_S_Dyz_S_bb+ABX*I_ERI_Gx2yz_S_Dyz_S_bb;
  Double I_ERI_Gxy2z_Px_Dyz_S_bb = I_ERI_H2xy2z_S_Dyz_S_bb+ABX*I_ERI_Gxy2z_S_Dyz_S_bb;
  Double I_ERI_Gx3z_Px_Dyz_S_bb = I_ERI_H2x3z_S_Dyz_S_bb+ABX*I_ERI_Gx3z_S_Dyz_S_bb;
  Double I_ERI_G4y_Px_Dyz_S_bb = I_ERI_Hx4y_S_Dyz_S_bb+ABX*I_ERI_G4y_S_Dyz_S_bb;
  Double I_ERI_G3yz_Px_Dyz_S_bb = I_ERI_Hx3yz_S_Dyz_S_bb+ABX*I_ERI_G3yz_S_Dyz_S_bb;
  Double I_ERI_G2y2z_Px_Dyz_S_bb = I_ERI_Hx2y2z_S_Dyz_S_bb+ABX*I_ERI_G2y2z_S_Dyz_S_bb;
  Double I_ERI_Gy3z_Px_Dyz_S_bb = I_ERI_Hxy3z_S_Dyz_S_bb+ABX*I_ERI_Gy3z_S_Dyz_S_bb;
  Double I_ERI_G4z_Px_Dyz_S_bb = I_ERI_Hx4z_S_Dyz_S_bb+ABX*I_ERI_G4z_S_Dyz_S_bb;
  Double I_ERI_G3xy_Py_Dyz_S_bb = I_ERI_H3x2y_S_Dyz_S_bb+ABY*I_ERI_G3xy_S_Dyz_S_bb;
  Double I_ERI_G3xz_Py_Dyz_S_bb = I_ERI_H3xyz_S_Dyz_S_bb+ABY*I_ERI_G3xz_S_Dyz_S_bb;
  Double I_ERI_G2x2y_Py_Dyz_S_bb = I_ERI_H2x3y_S_Dyz_S_bb+ABY*I_ERI_G2x2y_S_Dyz_S_bb;
  Double I_ERI_G2xyz_Py_Dyz_S_bb = I_ERI_H2x2yz_S_Dyz_S_bb+ABY*I_ERI_G2xyz_S_Dyz_S_bb;
  Double I_ERI_G2x2z_Py_Dyz_S_bb = I_ERI_H2xy2z_S_Dyz_S_bb+ABY*I_ERI_G2x2z_S_Dyz_S_bb;
  Double I_ERI_Gx3y_Py_Dyz_S_bb = I_ERI_Hx4y_S_Dyz_S_bb+ABY*I_ERI_Gx3y_S_Dyz_S_bb;
  Double I_ERI_Gx2yz_Py_Dyz_S_bb = I_ERI_Hx3yz_S_Dyz_S_bb+ABY*I_ERI_Gx2yz_S_Dyz_S_bb;
  Double I_ERI_Gxy2z_Py_Dyz_S_bb = I_ERI_Hx2y2z_S_Dyz_S_bb+ABY*I_ERI_Gxy2z_S_Dyz_S_bb;
  Double I_ERI_Gx3z_Py_Dyz_S_bb = I_ERI_Hxy3z_S_Dyz_S_bb+ABY*I_ERI_Gx3z_S_Dyz_S_bb;
  Double I_ERI_G4y_Py_Dyz_S_bb = I_ERI_H5y_S_Dyz_S_bb+ABY*I_ERI_G4y_S_Dyz_S_bb;
  Double I_ERI_G3yz_Py_Dyz_S_bb = I_ERI_H4yz_S_Dyz_S_bb+ABY*I_ERI_G3yz_S_Dyz_S_bb;
  Double I_ERI_G2y2z_Py_Dyz_S_bb = I_ERI_H3y2z_S_Dyz_S_bb+ABY*I_ERI_G2y2z_S_Dyz_S_bb;
  Double I_ERI_Gy3z_Py_Dyz_S_bb = I_ERI_H2y3z_S_Dyz_S_bb+ABY*I_ERI_Gy3z_S_Dyz_S_bb;
  Double I_ERI_G4z_Py_Dyz_S_bb = I_ERI_Hy4z_S_Dyz_S_bb+ABY*I_ERI_G4z_S_Dyz_S_bb;
  Double I_ERI_G3xz_Pz_Dyz_S_bb = I_ERI_H3x2z_S_Dyz_S_bb+ABZ*I_ERI_G3xz_S_Dyz_S_bb;
  Double I_ERI_G2xyz_Pz_Dyz_S_bb = I_ERI_H2xy2z_S_Dyz_S_bb+ABZ*I_ERI_G2xyz_S_Dyz_S_bb;
  Double I_ERI_G2x2z_Pz_Dyz_S_bb = I_ERI_H2x3z_S_Dyz_S_bb+ABZ*I_ERI_G2x2z_S_Dyz_S_bb;
  Double I_ERI_Gx2yz_Pz_Dyz_S_bb = I_ERI_Hx2y2z_S_Dyz_S_bb+ABZ*I_ERI_Gx2yz_S_Dyz_S_bb;
  Double I_ERI_Gxy2z_Pz_Dyz_S_bb = I_ERI_Hxy3z_S_Dyz_S_bb+ABZ*I_ERI_Gxy2z_S_Dyz_S_bb;
  Double I_ERI_Gx3z_Pz_Dyz_S_bb = I_ERI_Hx4z_S_Dyz_S_bb+ABZ*I_ERI_Gx3z_S_Dyz_S_bb;
  Double I_ERI_G3yz_Pz_Dyz_S_bb = I_ERI_H3y2z_S_Dyz_S_bb+ABZ*I_ERI_G3yz_S_Dyz_S_bb;
  Double I_ERI_G2y2z_Pz_Dyz_S_bb = I_ERI_H2y3z_S_Dyz_S_bb+ABZ*I_ERI_G2y2z_S_Dyz_S_bb;
  Double I_ERI_Gy3z_Pz_Dyz_S_bb = I_ERI_Hy4z_S_Dyz_S_bb+ABZ*I_ERI_Gy3z_S_Dyz_S_bb;
  Double I_ERI_G4z_Pz_Dyz_S_bb = I_ERI_H5z_S_Dyz_S_bb+ABZ*I_ERI_G4z_S_Dyz_S_bb;
  Double I_ERI_G4x_Px_D2z_S_bb = I_ERI_H5x_S_D2z_S_bb+ABX*I_ERI_G4x_S_D2z_S_bb;
  Double I_ERI_G3xy_Px_D2z_S_bb = I_ERI_H4xy_S_D2z_S_bb+ABX*I_ERI_G3xy_S_D2z_S_bb;
  Double I_ERI_G3xz_Px_D2z_S_bb = I_ERI_H4xz_S_D2z_S_bb+ABX*I_ERI_G3xz_S_D2z_S_bb;
  Double I_ERI_G2x2y_Px_D2z_S_bb = I_ERI_H3x2y_S_D2z_S_bb+ABX*I_ERI_G2x2y_S_D2z_S_bb;
  Double I_ERI_G2xyz_Px_D2z_S_bb = I_ERI_H3xyz_S_D2z_S_bb+ABX*I_ERI_G2xyz_S_D2z_S_bb;
  Double I_ERI_G2x2z_Px_D2z_S_bb = I_ERI_H3x2z_S_D2z_S_bb+ABX*I_ERI_G2x2z_S_D2z_S_bb;
  Double I_ERI_Gx3y_Px_D2z_S_bb = I_ERI_H2x3y_S_D2z_S_bb+ABX*I_ERI_Gx3y_S_D2z_S_bb;
  Double I_ERI_Gx2yz_Px_D2z_S_bb = I_ERI_H2x2yz_S_D2z_S_bb+ABX*I_ERI_Gx2yz_S_D2z_S_bb;
  Double I_ERI_Gxy2z_Px_D2z_S_bb = I_ERI_H2xy2z_S_D2z_S_bb+ABX*I_ERI_Gxy2z_S_D2z_S_bb;
  Double I_ERI_Gx3z_Px_D2z_S_bb = I_ERI_H2x3z_S_D2z_S_bb+ABX*I_ERI_Gx3z_S_D2z_S_bb;
  Double I_ERI_G4y_Px_D2z_S_bb = I_ERI_Hx4y_S_D2z_S_bb+ABX*I_ERI_G4y_S_D2z_S_bb;
  Double I_ERI_G3yz_Px_D2z_S_bb = I_ERI_Hx3yz_S_D2z_S_bb+ABX*I_ERI_G3yz_S_D2z_S_bb;
  Double I_ERI_G2y2z_Px_D2z_S_bb = I_ERI_Hx2y2z_S_D2z_S_bb+ABX*I_ERI_G2y2z_S_D2z_S_bb;
  Double I_ERI_Gy3z_Px_D2z_S_bb = I_ERI_Hxy3z_S_D2z_S_bb+ABX*I_ERI_Gy3z_S_D2z_S_bb;
  Double I_ERI_G4z_Px_D2z_S_bb = I_ERI_Hx4z_S_D2z_S_bb+ABX*I_ERI_G4z_S_D2z_S_bb;
  Double I_ERI_G3xy_Py_D2z_S_bb = I_ERI_H3x2y_S_D2z_S_bb+ABY*I_ERI_G3xy_S_D2z_S_bb;
  Double I_ERI_G3xz_Py_D2z_S_bb = I_ERI_H3xyz_S_D2z_S_bb+ABY*I_ERI_G3xz_S_D2z_S_bb;
  Double I_ERI_G2x2y_Py_D2z_S_bb = I_ERI_H2x3y_S_D2z_S_bb+ABY*I_ERI_G2x2y_S_D2z_S_bb;
  Double I_ERI_G2xyz_Py_D2z_S_bb = I_ERI_H2x2yz_S_D2z_S_bb+ABY*I_ERI_G2xyz_S_D2z_S_bb;
  Double I_ERI_G2x2z_Py_D2z_S_bb = I_ERI_H2xy2z_S_D2z_S_bb+ABY*I_ERI_G2x2z_S_D2z_S_bb;
  Double I_ERI_Gx3y_Py_D2z_S_bb = I_ERI_Hx4y_S_D2z_S_bb+ABY*I_ERI_Gx3y_S_D2z_S_bb;
  Double I_ERI_Gx2yz_Py_D2z_S_bb = I_ERI_Hx3yz_S_D2z_S_bb+ABY*I_ERI_Gx2yz_S_D2z_S_bb;
  Double I_ERI_Gxy2z_Py_D2z_S_bb = I_ERI_Hx2y2z_S_D2z_S_bb+ABY*I_ERI_Gxy2z_S_D2z_S_bb;
  Double I_ERI_Gx3z_Py_D2z_S_bb = I_ERI_Hxy3z_S_D2z_S_bb+ABY*I_ERI_Gx3z_S_D2z_S_bb;
  Double I_ERI_G4y_Py_D2z_S_bb = I_ERI_H5y_S_D2z_S_bb+ABY*I_ERI_G4y_S_D2z_S_bb;
  Double I_ERI_G3yz_Py_D2z_S_bb = I_ERI_H4yz_S_D2z_S_bb+ABY*I_ERI_G3yz_S_D2z_S_bb;
  Double I_ERI_G2y2z_Py_D2z_S_bb = I_ERI_H3y2z_S_D2z_S_bb+ABY*I_ERI_G2y2z_S_D2z_S_bb;
  Double I_ERI_Gy3z_Py_D2z_S_bb = I_ERI_H2y3z_S_D2z_S_bb+ABY*I_ERI_Gy3z_S_D2z_S_bb;
  Double I_ERI_G4z_Py_D2z_S_bb = I_ERI_Hy4z_S_D2z_S_bb+ABY*I_ERI_G4z_S_D2z_S_bb;
  Double I_ERI_G3xz_Pz_D2z_S_bb = I_ERI_H3x2z_S_D2z_S_bb+ABZ*I_ERI_G3xz_S_D2z_S_bb;
  Double I_ERI_G2xyz_Pz_D2z_S_bb = I_ERI_H2xy2z_S_D2z_S_bb+ABZ*I_ERI_G2xyz_S_D2z_S_bb;
  Double I_ERI_G2x2z_Pz_D2z_S_bb = I_ERI_H2x3z_S_D2z_S_bb+ABZ*I_ERI_G2x2z_S_D2z_S_bb;
  Double I_ERI_Gx2yz_Pz_D2z_S_bb = I_ERI_Hx2y2z_S_D2z_S_bb+ABZ*I_ERI_Gx2yz_S_D2z_S_bb;
  Double I_ERI_Gxy2z_Pz_D2z_S_bb = I_ERI_Hxy3z_S_D2z_S_bb+ABZ*I_ERI_Gxy2z_S_D2z_S_bb;
  Double I_ERI_Gx3z_Pz_D2z_S_bb = I_ERI_Hx4z_S_D2z_S_bb+ABZ*I_ERI_Gx3z_S_D2z_S_bb;
  Double I_ERI_G3yz_Pz_D2z_S_bb = I_ERI_H3y2z_S_D2z_S_bb+ABZ*I_ERI_G3yz_S_D2z_S_bb;
  Double I_ERI_G2y2z_Pz_D2z_S_bb = I_ERI_H2y3z_S_D2z_S_bb+ABZ*I_ERI_G2y2z_S_D2z_S_bb;
  Double I_ERI_Gy3z_Pz_D2z_S_bb = I_ERI_Hy4z_S_D2z_S_bb+ABZ*I_ERI_Gy3z_S_D2z_S_bb;
  Double I_ERI_G4z_Pz_D2z_S_bb = I_ERI_H5z_S_D2z_S_bb+ABZ*I_ERI_G4z_S_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_D_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bb
   ************************************************************/
  Double I_ERI_F3x_D2x_D2x_S_bb = I_ERI_G4x_Px_D2x_S_bb+ABX*I_ERI_F3x_Px_D2x_S_bb;
  Double I_ERI_F2xy_D2x_D2x_S_bb = I_ERI_G3xy_Px_D2x_S_bb+ABX*I_ERI_F2xy_Px_D2x_S_bb;
  Double I_ERI_F2xz_D2x_D2x_S_bb = I_ERI_G3xz_Px_D2x_S_bb+ABX*I_ERI_F2xz_Px_D2x_S_bb;
  Double I_ERI_Fx2y_D2x_D2x_S_bb = I_ERI_G2x2y_Px_D2x_S_bb+ABX*I_ERI_Fx2y_Px_D2x_S_bb;
  Double I_ERI_Fxyz_D2x_D2x_S_bb = I_ERI_G2xyz_Px_D2x_S_bb+ABX*I_ERI_Fxyz_Px_D2x_S_bb;
  Double I_ERI_Fx2z_D2x_D2x_S_bb = I_ERI_G2x2z_Px_D2x_S_bb+ABX*I_ERI_Fx2z_Px_D2x_S_bb;
  Double I_ERI_F3y_D2x_D2x_S_bb = I_ERI_Gx3y_Px_D2x_S_bb+ABX*I_ERI_F3y_Px_D2x_S_bb;
  Double I_ERI_F2yz_D2x_D2x_S_bb = I_ERI_Gx2yz_Px_D2x_S_bb+ABX*I_ERI_F2yz_Px_D2x_S_bb;
  Double I_ERI_Fy2z_D2x_D2x_S_bb = I_ERI_Gxy2z_Px_D2x_S_bb+ABX*I_ERI_Fy2z_Px_D2x_S_bb;
  Double I_ERI_F3z_D2x_D2x_S_bb = I_ERI_Gx3z_Px_D2x_S_bb+ABX*I_ERI_F3z_Px_D2x_S_bb;
  Double I_ERI_F3x_Dxy_D2x_S_bb = I_ERI_G3xy_Px_D2x_S_bb+ABY*I_ERI_F3x_Px_D2x_S_bb;
  Double I_ERI_F2xy_Dxy_D2x_S_bb = I_ERI_G2x2y_Px_D2x_S_bb+ABY*I_ERI_F2xy_Px_D2x_S_bb;
  Double I_ERI_F2xz_Dxy_D2x_S_bb = I_ERI_G2xyz_Px_D2x_S_bb+ABY*I_ERI_F2xz_Px_D2x_S_bb;
  Double I_ERI_Fx2y_Dxy_D2x_S_bb = I_ERI_Gx3y_Px_D2x_S_bb+ABY*I_ERI_Fx2y_Px_D2x_S_bb;
  Double I_ERI_Fxyz_Dxy_D2x_S_bb = I_ERI_Gx2yz_Px_D2x_S_bb+ABY*I_ERI_Fxyz_Px_D2x_S_bb;
  Double I_ERI_Fx2z_Dxy_D2x_S_bb = I_ERI_Gxy2z_Px_D2x_S_bb+ABY*I_ERI_Fx2z_Px_D2x_S_bb;
  Double I_ERI_F3y_Dxy_D2x_S_bb = I_ERI_G4y_Px_D2x_S_bb+ABY*I_ERI_F3y_Px_D2x_S_bb;
  Double I_ERI_F2yz_Dxy_D2x_S_bb = I_ERI_G3yz_Px_D2x_S_bb+ABY*I_ERI_F2yz_Px_D2x_S_bb;
  Double I_ERI_Fy2z_Dxy_D2x_S_bb = I_ERI_G2y2z_Px_D2x_S_bb+ABY*I_ERI_Fy2z_Px_D2x_S_bb;
  Double I_ERI_F3z_Dxy_D2x_S_bb = I_ERI_Gy3z_Px_D2x_S_bb+ABY*I_ERI_F3z_Px_D2x_S_bb;
  Double I_ERI_F3x_Dxz_D2x_S_bb = I_ERI_G3xz_Px_D2x_S_bb+ABZ*I_ERI_F3x_Px_D2x_S_bb;
  Double I_ERI_F2xy_Dxz_D2x_S_bb = I_ERI_G2xyz_Px_D2x_S_bb+ABZ*I_ERI_F2xy_Px_D2x_S_bb;
  Double I_ERI_F2xz_Dxz_D2x_S_bb = I_ERI_G2x2z_Px_D2x_S_bb+ABZ*I_ERI_F2xz_Px_D2x_S_bb;
  Double I_ERI_Fx2y_Dxz_D2x_S_bb = I_ERI_Gx2yz_Px_D2x_S_bb+ABZ*I_ERI_Fx2y_Px_D2x_S_bb;
  Double I_ERI_Fxyz_Dxz_D2x_S_bb = I_ERI_Gxy2z_Px_D2x_S_bb+ABZ*I_ERI_Fxyz_Px_D2x_S_bb;
  Double I_ERI_Fx2z_Dxz_D2x_S_bb = I_ERI_Gx3z_Px_D2x_S_bb+ABZ*I_ERI_Fx2z_Px_D2x_S_bb;
  Double I_ERI_F3y_Dxz_D2x_S_bb = I_ERI_G3yz_Px_D2x_S_bb+ABZ*I_ERI_F3y_Px_D2x_S_bb;
  Double I_ERI_F2yz_Dxz_D2x_S_bb = I_ERI_G2y2z_Px_D2x_S_bb+ABZ*I_ERI_F2yz_Px_D2x_S_bb;
  Double I_ERI_Fy2z_Dxz_D2x_S_bb = I_ERI_Gy3z_Px_D2x_S_bb+ABZ*I_ERI_Fy2z_Px_D2x_S_bb;
  Double I_ERI_F3z_Dxz_D2x_S_bb = I_ERI_G4z_Px_D2x_S_bb+ABZ*I_ERI_F3z_Px_D2x_S_bb;
  Double I_ERI_F3x_D2y_D2x_S_bb = I_ERI_G3xy_Py_D2x_S_bb+ABY*I_ERI_F3x_Py_D2x_S_bb;
  Double I_ERI_F2xy_D2y_D2x_S_bb = I_ERI_G2x2y_Py_D2x_S_bb+ABY*I_ERI_F2xy_Py_D2x_S_bb;
  Double I_ERI_F2xz_D2y_D2x_S_bb = I_ERI_G2xyz_Py_D2x_S_bb+ABY*I_ERI_F2xz_Py_D2x_S_bb;
  Double I_ERI_Fx2y_D2y_D2x_S_bb = I_ERI_Gx3y_Py_D2x_S_bb+ABY*I_ERI_Fx2y_Py_D2x_S_bb;
  Double I_ERI_Fxyz_D2y_D2x_S_bb = I_ERI_Gx2yz_Py_D2x_S_bb+ABY*I_ERI_Fxyz_Py_D2x_S_bb;
  Double I_ERI_Fx2z_D2y_D2x_S_bb = I_ERI_Gxy2z_Py_D2x_S_bb+ABY*I_ERI_Fx2z_Py_D2x_S_bb;
  Double I_ERI_F3y_D2y_D2x_S_bb = I_ERI_G4y_Py_D2x_S_bb+ABY*I_ERI_F3y_Py_D2x_S_bb;
  Double I_ERI_F2yz_D2y_D2x_S_bb = I_ERI_G3yz_Py_D2x_S_bb+ABY*I_ERI_F2yz_Py_D2x_S_bb;
  Double I_ERI_Fy2z_D2y_D2x_S_bb = I_ERI_G2y2z_Py_D2x_S_bb+ABY*I_ERI_Fy2z_Py_D2x_S_bb;
  Double I_ERI_F3z_D2y_D2x_S_bb = I_ERI_Gy3z_Py_D2x_S_bb+ABY*I_ERI_F3z_Py_D2x_S_bb;
  Double I_ERI_F3x_Dyz_D2x_S_bb = I_ERI_G3xz_Py_D2x_S_bb+ABZ*I_ERI_F3x_Py_D2x_S_bb;
  Double I_ERI_F2xy_Dyz_D2x_S_bb = I_ERI_G2xyz_Py_D2x_S_bb+ABZ*I_ERI_F2xy_Py_D2x_S_bb;
  Double I_ERI_F2xz_Dyz_D2x_S_bb = I_ERI_G2x2z_Py_D2x_S_bb+ABZ*I_ERI_F2xz_Py_D2x_S_bb;
  Double I_ERI_Fx2y_Dyz_D2x_S_bb = I_ERI_Gx2yz_Py_D2x_S_bb+ABZ*I_ERI_Fx2y_Py_D2x_S_bb;
  Double I_ERI_Fxyz_Dyz_D2x_S_bb = I_ERI_Gxy2z_Py_D2x_S_bb+ABZ*I_ERI_Fxyz_Py_D2x_S_bb;
  Double I_ERI_Fx2z_Dyz_D2x_S_bb = I_ERI_Gx3z_Py_D2x_S_bb+ABZ*I_ERI_Fx2z_Py_D2x_S_bb;
  Double I_ERI_F3y_Dyz_D2x_S_bb = I_ERI_G3yz_Py_D2x_S_bb+ABZ*I_ERI_F3y_Py_D2x_S_bb;
  Double I_ERI_F2yz_Dyz_D2x_S_bb = I_ERI_G2y2z_Py_D2x_S_bb+ABZ*I_ERI_F2yz_Py_D2x_S_bb;
  Double I_ERI_Fy2z_Dyz_D2x_S_bb = I_ERI_Gy3z_Py_D2x_S_bb+ABZ*I_ERI_Fy2z_Py_D2x_S_bb;
  Double I_ERI_F3z_Dyz_D2x_S_bb = I_ERI_G4z_Py_D2x_S_bb+ABZ*I_ERI_F3z_Py_D2x_S_bb;
  Double I_ERI_F3x_D2z_D2x_S_bb = I_ERI_G3xz_Pz_D2x_S_bb+ABZ*I_ERI_F3x_Pz_D2x_S_bb;
  Double I_ERI_F2xy_D2z_D2x_S_bb = I_ERI_G2xyz_Pz_D2x_S_bb+ABZ*I_ERI_F2xy_Pz_D2x_S_bb;
  Double I_ERI_F2xz_D2z_D2x_S_bb = I_ERI_G2x2z_Pz_D2x_S_bb+ABZ*I_ERI_F2xz_Pz_D2x_S_bb;
  Double I_ERI_Fx2y_D2z_D2x_S_bb = I_ERI_Gx2yz_Pz_D2x_S_bb+ABZ*I_ERI_Fx2y_Pz_D2x_S_bb;
  Double I_ERI_Fxyz_D2z_D2x_S_bb = I_ERI_Gxy2z_Pz_D2x_S_bb+ABZ*I_ERI_Fxyz_Pz_D2x_S_bb;
  Double I_ERI_Fx2z_D2z_D2x_S_bb = I_ERI_Gx3z_Pz_D2x_S_bb+ABZ*I_ERI_Fx2z_Pz_D2x_S_bb;
  Double I_ERI_F3y_D2z_D2x_S_bb = I_ERI_G3yz_Pz_D2x_S_bb+ABZ*I_ERI_F3y_Pz_D2x_S_bb;
  Double I_ERI_F2yz_D2z_D2x_S_bb = I_ERI_G2y2z_Pz_D2x_S_bb+ABZ*I_ERI_F2yz_Pz_D2x_S_bb;
  Double I_ERI_Fy2z_D2z_D2x_S_bb = I_ERI_Gy3z_Pz_D2x_S_bb+ABZ*I_ERI_Fy2z_Pz_D2x_S_bb;
  Double I_ERI_F3z_D2z_D2x_S_bb = I_ERI_G4z_Pz_D2x_S_bb+ABZ*I_ERI_F3z_Pz_D2x_S_bb;
  Double I_ERI_F3x_D2x_Dxy_S_bb = I_ERI_G4x_Px_Dxy_S_bb+ABX*I_ERI_F3x_Px_Dxy_S_bb;
  Double I_ERI_F2xy_D2x_Dxy_S_bb = I_ERI_G3xy_Px_Dxy_S_bb+ABX*I_ERI_F2xy_Px_Dxy_S_bb;
  Double I_ERI_F2xz_D2x_Dxy_S_bb = I_ERI_G3xz_Px_Dxy_S_bb+ABX*I_ERI_F2xz_Px_Dxy_S_bb;
  Double I_ERI_Fx2y_D2x_Dxy_S_bb = I_ERI_G2x2y_Px_Dxy_S_bb+ABX*I_ERI_Fx2y_Px_Dxy_S_bb;
  Double I_ERI_Fxyz_D2x_Dxy_S_bb = I_ERI_G2xyz_Px_Dxy_S_bb+ABX*I_ERI_Fxyz_Px_Dxy_S_bb;
  Double I_ERI_Fx2z_D2x_Dxy_S_bb = I_ERI_G2x2z_Px_Dxy_S_bb+ABX*I_ERI_Fx2z_Px_Dxy_S_bb;
  Double I_ERI_F3y_D2x_Dxy_S_bb = I_ERI_Gx3y_Px_Dxy_S_bb+ABX*I_ERI_F3y_Px_Dxy_S_bb;
  Double I_ERI_F2yz_D2x_Dxy_S_bb = I_ERI_Gx2yz_Px_Dxy_S_bb+ABX*I_ERI_F2yz_Px_Dxy_S_bb;
  Double I_ERI_Fy2z_D2x_Dxy_S_bb = I_ERI_Gxy2z_Px_Dxy_S_bb+ABX*I_ERI_Fy2z_Px_Dxy_S_bb;
  Double I_ERI_F3z_D2x_Dxy_S_bb = I_ERI_Gx3z_Px_Dxy_S_bb+ABX*I_ERI_F3z_Px_Dxy_S_bb;
  Double I_ERI_F3x_Dxy_Dxy_S_bb = I_ERI_G3xy_Px_Dxy_S_bb+ABY*I_ERI_F3x_Px_Dxy_S_bb;
  Double I_ERI_F2xy_Dxy_Dxy_S_bb = I_ERI_G2x2y_Px_Dxy_S_bb+ABY*I_ERI_F2xy_Px_Dxy_S_bb;
  Double I_ERI_F2xz_Dxy_Dxy_S_bb = I_ERI_G2xyz_Px_Dxy_S_bb+ABY*I_ERI_F2xz_Px_Dxy_S_bb;
  Double I_ERI_Fx2y_Dxy_Dxy_S_bb = I_ERI_Gx3y_Px_Dxy_S_bb+ABY*I_ERI_Fx2y_Px_Dxy_S_bb;
  Double I_ERI_Fxyz_Dxy_Dxy_S_bb = I_ERI_Gx2yz_Px_Dxy_S_bb+ABY*I_ERI_Fxyz_Px_Dxy_S_bb;
  Double I_ERI_Fx2z_Dxy_Dxy_S_bb = I_ERI_Gxy2z_Px_Dxy_S_bb+ABY*I_ERI_Fx2z_Px_Dxy_S_bb;
  Double I_ERI_F3y_Dxy_Dxy_S_bb = I_ERI_G4y_Px_Dxy_S_bb+ABY*I_ERI_F3y_Px_Dxy_S_bb;
  Double I_ERI_F2yz_Dxy_Dxy_S_bb = I_ERI_G3yz_Px_Dxy_S_bb+ABY*I_ERI_F2yz_Px_Dxy_S_bb;
  Double I_ERI_Fy2z_Dxy_Dxy_S_bb = I_ERI_G2y2z_Px_Dxy_S_bb+ABY*I_ERI_Fy2z_Px_Dxy_S_bb;
  Double I_ERI_F3z_Dxy_Dxy_S_bb = I_ERI_Gy3z_Px_Dxy_S_bb+ABY*I_ERI_F3z_Px_Dxy_S_bb;
  Double I_ERI_F3x_Dxz_Dxy_S_bb = I_ERI_G3xz_Px_Dxy_S_bb+ABZ*I_ERI_F3x_Px_Dxy_S_bb;
  Double I_ERI_F2xy_Dxz_Dxy_S_bb = I_ERI_G2xyz_Px_Dxy_S_bb+ABZ*I_ERI_F2xy_Px_Dxy_S_bb;
  Double I_ERI_F2xz_Dxz_Dxy_S_bb = I_ERI_G2x2z_Px_Dxy_S_bb+ABZ*I_ERI_F2xz_Px_Dxy_S_bb;
  Double I_ERI_Fx2y_Dxz_Dxy_S_bb = I_ERI_Gx2yz_Px_Dxy_S_bb+ABZ*I_ERI_Fx2y_Px_Dxy_S_bb;
  Double I_ERI_Fxyz_Dxz_Dxy_S_bb = I_ERI_Gxy2z_Px_Dxy_S_bb+ABZ*I_ERI_Fxyz_Px_Dxy_S_bb;
  Double I_ERI_Fx2z_Dxz_Dxy_S_bb = I_ERI_Gx3z_Px_Dxy_S_bb+ABZ*I_ERI_Fx2z_Px_Dxy_S_bb;
  Double I_ERI_F3y_Dxz_Dxy_S_bb = I_ERI_G3yz_Px_Dxy_S_bb+ABZ*I_ERI_F3y_Px_Dxy_S_bb;
  Double I_ERI_F2yz_Dxz_Dxy_S_bb = I_ERI_G2y2z_Px_Dxy_S_bb+ABZ*I_ERI_F2yz_Px_Dxy_S_bb;
  Double I_ERI_Fy2z_Dxz_Dxy_S_bb = I_ERI_Gy3z_Px_Dxy_S_bb+ABZ*I_ERI_Fy2z_Px_Dxy_S_bb;
  Double I_ERI_F3z_Dxz_Dxy_S_bb = I_ERI_G4z_Px_Dxy_S_bb+ABZ*I_ERI_F3z_Px_Dxy_S_bb;
  Double I_ERI_F3x_D2y_Dxy_S_bb = I_ERI_G3xy_Py_Dxy_S_bb+ABY*I_ERI_F3x_Py_Dxy_S_bb;
  Double I_ERI_F2xy_D2y_Dxy_S_bb = I_ERI_G2x2y_Py_Dxy_S_bb+ABY*I_ERI_F2xy_Py_Dxy_S_bb;
  Double I_ERI_F2xz_D2y_Dxy_S_bb = I_ERI_G2xyz_Py_Dxy_S_bb+ABY*I_ERI_F2xz_Py_Dxy_S_bb;
  Double I_ERI_Fx2y_D2y_Dxy_S_bb = I_ERI_Gx3y_Py_Dxy_S_bb+ABY*I_ERI_Fx2y_Py_Dxy_S_bb;
  Double I_ERI_Fxyz_D2y_Dxy_S_bb = I_ERI_Gx2yz_Py_Dxy_S_bb+ABY*I_ERI_Fxyz_Py_Dxy_S_bb;
  Double I_ERI_Fx2z_D2y_Dxy_S_bb = I_ERI_Gxy2z_Py_Dxy_S_bb+ABY*I_ERI_Fx2z_Py_Dxy_S_bb;
  Double I_ERI_F3y_D2y_Dxy_S_bb = I_ERI_G4y_Py_Dxy_S_bb+ABY*I_ERI_F3y_Py_Dxy_S_bb;
  Double I_ERI_F2yz_D2y_Dxy_S_bb = I_ERI_G3yz_Py_Dxy_S_bb+ABY*I_ERI_F2yz_Py_Dxy_S_bb;
  Double I_ERI_Fy2z_D2y_Dxy_S_bb = I_ERI_G2y2z_Py_Dxy_S_bb+ABY*I_ERI_Fy2z_Py_Dxy_S_bb;
  Double I_ERI_F3z_D2y_Dxy_S_bb = I_ERI_Gy3z_Py_Dxy_S_bb+ABY*I_ERI_F3z_Py_Dxy_S_bb;
  Double I_ERI_F3x_Dyz_Dxy_S_bb = I_ERI_G3xz_Py_Dxy_S_bb+ABZ*I_ERI_F3x_Py_Dxy_S_bb;
  Double I_ERI_F2xy_Dyz_Dxy_S_bb = I_ERI_G2xyz_Py_Dxy_S_bb+ABZ*I_ERI_F2xy_Py_Dxy_S_bb;
  Double I_ERI_F2xz_Dyz_Dxy_S_bb = I_ERI_G2x2z_Py_Dxy_S_bb+ABZ*I_ERI_F2xz_Py_Dxy_S_bb;
  Double I_ERI_Fx2y_Dyz_Dxy_S_bb = I_ERI_Gx2yz_Py_Dxy_S_bb+ABZ*I_ERI_Fx2y_Py_Dxy_S_bb;
  Double I_ERI_Fxyz_Dyz_Dxy_S_bb = I_ERI_Gxy2z_Py_Dxy_S_bb+ABZ*I_ERI_Fxyz_Py_Dxy_S_bb;
  Double I_ERI_Fx2z_Dyz_Dxy_S_bb = I_ERI_Gx3z_Py_Dxy_S_bb+ABZ*I_ERI_Fx2z_Py_Dxy_S_bb;
  Double I_ERI_F3y_Dyz_Dxy_S_bb = I_ERI_G3yz_Py_Dxy_S_bb+ABZ*I_ERI_F3y_Py_Dxy_S_bb;
  Double I_ERI_F2yz_Dyz_Dxy_S_bb = I_ERI_G2y2z_Py_Dxy_S_bb+ABZ*I_ERI_F2yz_Py_Dxy_S_bb;
  Double I_ERI_Fy2z_Dyz_Dxy_S_bb = I_ERI_Gy3z_Py_Dxy_S_bb+ABZ*I_ERI_Fy2z_Py_Dxy_S_bb;
  Double I_ERI_F3z_Dyz_Dxy_S_bb = I_ERI_G4z_Py_Dxy_S_bb+ABZ*I_ERI_F3z_Py_Dxy_S_bb;
  Double I_ERI_F3x_D2z_Dxy_S_bb = I_ERI_G3xz_Pz_Dxy_S_bb+ABZ*I_ERI_F3x_Pz_Dxy_S_bb;
  Double I_ERI_F2xy_D2z_Dxy_S_bb = I_ERI_G2xyz_Pz_Dxy_S_bb+ABZ*I_ERI_F2xy_Pz_Dxy_S_bb;
  Double I_ERI_F2xz_D2z_Dxy_S_bb = I_ERI_G2x2z_Pz_Dxy_S_bb+ABZ*I_ERI_F2xz_Pz_Dxy_S_bb;
  Double I_ERI_Fx2y_D2z_Dxy_S_bb = I_ERI_Gx2yz_Pz_Dxy_S_bb+ABZ*I_ERI_Fx2y_Pz_Dxy_S_bb;
  Double I_ERI_Fxyz_D2z_Dxy_S_bb = I_ERI_Gxy2z_Pz_Dxy_S_bb+ABZ*I_ERI_Fxyz_Pz_Dxy_S_bb;
  Double I_ERI_Fx2z_D2z_Dxy_S_bb = I_ERI_Gx3z_Pz_Dxy_S_bb+ABZ*I_ERI_Fx2z_Pz_Dxy_S_bb;
  Double I_ERI_F3y_D2z_Dxy_S_bb = I_ERI_G3yz_Pz_Dxy_S_bb+ABZ*I_ERI_F3y_Pz_Dxy_S_bb;
  Double I_ERI_F2yz_D2z_Dxy_S_bb = I_ERI_G2y2z_Pz_Dxy_S_bb+ABZ*I_ERI_F2yz_Pz_Dxy_S_bb;
  Double I_ERI_Fy2z_D2z_Dxy_S_bb = I_ERI_Gy3z_Pz_Dxy_S_bb+ABZ*I_ERI_Fy2z_Pz_Dxy_S_bb;
  Double I_ERI_F3z_D2z_Dxy_S_bb = I_ERI_G4z_Pz_Dxy_S_bb+ABZ*I_ERI_F3z_Pz_Dxy_S_bb;
  Double I_ERI_F3x_D2x_Dxz_S_bb = I_ERI_G4x_Px_Dxz_S_bb+ABX*I_ERI_F3x_Px_Dxz_S_bb;
  Double I_ERI_F2xy_D2x_Dxz_S_bb = I_ERI_G3xy_Px_Dxz_S_bb+ABX*I_ERI_F2xy_Px_Dxz_S_bb;
  Double I_ERI_F2xz_D2x_Dxz_S_bb = I_ERI_G3xz_Px_Dxz_S_bb+ABX*I_ERI_F2xz_Px_Dxz_S_bb;
  Double I_ERI_Fx2y_D2x_Dxz_S_bb = I_ERI_G2x2y_Px_Dxz_S_bb+ABX*I_ERI_Fx2y_Px_Dxz_S_bb;
  Double I_ERI_Fxyz_D2x_Dxz_S_bb = I_ERI_G2xyz_Px_Dxz_S_bb+ABX*I_ERI_Fxyz_Px_Dxz_S_bb;
  Double I_ERI_Fx2z_D2x_Dxz_S_bb = I_ERI_G2x2z_Px_Dxz_S_bb+ABX*I_ERI_Fx2z_Px_Dxz_S_bb;
  Double I_ERI_F3y_D2x_Dxz_S_bb = I_ERI_Gx3y_Px_Dxz_S_bb+ABX*I_ERI_F3y_Px_Dxz_S_bb;
  Double I_ERI_F2yz_D2x_Dxz_S_bb = I_ERI_Gx2yz_Px_Dxz_S_bb+ABX*I_ERI_F2yz_Px_Dxz_S_bb;
  Double I_ERI_Fy2z_D2x_Dxz_S_bb = I_ERI_Gxy2z_Px_Dxz_S_bb+ABX*I_ERI_Fy2z_Px_Dxz_S_bb;
  Double I_ERI_F3z_D2x_Dxz_S_bb = I_ERI_Gx3z_Px_Dxz_S_bb+ABX*I_ERI_F3z_Px_Dxz_S_bb;
  Double I_ERI_F3x_Dxy_Dxz_S_bb = I_ERI_G3xy_Px_Dxz_S_bb+ABY*I_ERI_F3x_Px_Dxz_S_bb;
  Double I_ERI_F2xy_Dxy_Dxz_S_bb = I_ERI_G2x2y_Px_Dxz_S_bb+ABY*I_ERI_F2xy_Px_Dxz_S_bb;
  Double I_ERI_F2xz_Dxy_Dxz_S_bb = I_ERI_G2xyz_Px_Dxz_S_bb+ABY*I_ERI_F2xz_Px_Dxz_S_bb;
  Double I_ERI_Fx2y_Dxy_Dxz_S_bb = I_ERI_Gx3y_Px_Dxz_S_bb+ABY*I_ERI_Fx2y_Px_Dxz_S_bb;
  Double I_ERI_Fxyz_Dxy_Dxz_S_bb = I_ERI_Gx2yz_Px_Dxz_S_bb+ABY*I_ERI_Fxyz_Px_Dxz_S_bb;
  Double I_ERI_Fx2z_Dxy_Dxz_S_bb = I_ERI_Gxy2z_Px_Dxz_S_bb+ABY*I_ERI_Fx2z_Px_Dxz_S_bb;
  Double I_ERI_F3y_Dxy_Dxz_S_bb = I_ERI_G4y_Px_Dxz_S_bb+ABY*I_ERI_F3y_Px_Dxz_S_bb;
  Double I_ERI_F2yz_Dxy_Dxz_S_bb = I_ERI_G3yz_Px_Dxz_S_bb+ABY*I_ERI_F2yz_Px_Dxz_S_bb;
  Double I_ERI_Fy2z_Dxy_Dxz_S_bb = I_ERI_G2y2z_Px_Dxz_S_bb+ABY*I_ERI_Fy2z_Px_Dxz_S_bb;
  Double I_ERI_F3z_Dxy_Dxz_S_bb = I_ERI_Gy3z_Px_Dxz_S_bb+ABY*I_ERI_F3z_Px_Dxz_S_bb;
  Double I_ERI_F3x_Dxz_Dxz_S_bb = I_ERI_G3xz_Px_Dxz_S_bb+ABZ*I_ERI_F3x_Px_Dxz_S_bb;
  Double I_ERI_F2xy_Dxz_Dxz_S_bb = I_ERI_G2xyz_Px_Dxz_S_bb+ABZ*I_ERI_F2xy_Px_Dxz_S_bb;
  Double I_ERI_F2xz_Dxz_Dxz_S_bb = I_ERI_G2x2z_Px_Dxz_S_bb+ABZ*I_ERI_F2xz_Px_Dxz_S_bb;
  Double I_ERI_Fx2y_Dxz_Dxz_S_bb = I_ERI_Gx2yz_Px_Dxz_S_bb+ABZ*I_ERI_Fx2y_Px_Dxz_S_bb;
  Double I_ERI_Fxyz_Dxz_Dxz_S_bb = I_ERI_Gxy2z_Px_Dxz_S_bb+ABZ*I_ERI_Fxyz_Px_Dxz_S_bb;
  Double I_ERI_Fx2z_Dxz_Dxz_S_bb = I_ERI_Gx3z_Px_Dxz_S_bb+ABZ*I_ERI_Fx2z_Px_Dxz_S_bb;
  Double I_ERI_F3y_Dxz_Dxz_S_bb = I_ERI_G3yz_Px_Dxz_S_bb+ABZ*I_ERI_F3y_Px_Dxz_S_bb;
  Double I_ERI_F2yz_Dxz_Dxz_S_bb = I_ERI_G2y2z_Px_Dxz_S_bb+ABZ*I_ERI_F2yz_Px_Dxz_S_bb;
  Double I_ERI_Fy2z_Dxz_Dxz_S_bb = I_ERI_Gy3z_Px_Dxz_S_bb+ABZ*I_ERI_Fy2z_Px_Dxz_S_bb;
  Double I_ERI_F3z_Dxz_Dxz_S_bb = I_ERI_G4z_Px_Dxz_S_bb+ABZ*I_ERI_F3z_Px_Dxz_S_bb;
  Double I_ERI_F3x_D2y_Dxz_S_bb = I_ERI_G3xy_Py_Dxz_S_bb+ABY*I_ERI_F3x_Py_Dxz_S_bb;
  Double I_ERI_F2xy_D2y_Dxz_S_bb = I_ERI_G2x2y_Py_Dxz_S_bb+ABY*I_ERI_F2xy_Py_Dxz_S_bb;
  Double I_ERI_F2xz_D2y_Dxz_S_bb = I_ERI_G2xyz_Py_Dxz_S_bb+ABY*I_ERI_F2xz_Py_Dxz_S_bb;
  Double I_ERI_Fx2y_D2y_Dxz_S_bb = I_ERI_Gx3y_Py_Dxz_S_bb+ABY*I_ERI_Fx2y_Py_Dxz_S_bb;
  Double I_ERI_Fxyz_D2y_Dxz_S_bb = I_ERI_Gx2yz_Py_Dxz_S_bb+ABY*I_ERI_Fxyz_Py_Dxz_S_bb;
  Double I_ERI_Fx2z_D2y_Dxz_S_bb = I_ERI_Gxy2z_Py_Dxz_S_bb+ABY*I_ERI_Fx2z_Py_Dxz_S_bb;
  Double I_ERI_F3y_D2y_Dxz_S_bb = I_ERI_G4y_Py_Dxz_S_bb+ABY*I_ERI_F3y_Py_Dxz_S_bb;
  Double I_ERI_F2yz_D2y_Dxz_S_bb = I_ERI_G3yz_Py_Dxz_S_bb+ABY*I_ERI_F2yz_Py_Dxz_S_bb;
  Double I_ERI_Fy2z_D2y_Dxz_S_bb = I_ERI_G2y2z_Py_Dxz_S_bb+ABY*I_ERI_Fy2z_Py_Dxz_S_bb;
  Double I_ERI_F3z_D2y_Dxz_S_bb = I_ERI_Gy3z_Py_Dxz_S_bb+ABY*I_ERI_F3z_Py_Dxz_S_bb;
  Double I_ERI_F3x_Dyz_Dxz_S_bb = I_ERI_G3xz_Py_Dxz_S_bb+ABZ*I_ERI_F3x_Py_Dxz_S_bb;
  Double I_ERI_F2xy_Dyz_Dxz_S_bb = I_ERI_G2xyz_Py_Dxz_S_bb+ABZ*I_ERI_F2xy_Py_Dxz_S_bb;
  Double I_ERI_F2xz_Dyz_Dxz_S_bb = I_ERI_G2x2z_Py_Dxz_S_bb+ABZ*I_ERI_F2xz_Py_Dxz_S_bb;
  Double I_ERI_Fx2y_Dyz_Dxz_S_bb = I_ERI_Gx2yz_Py_Dxz_S_bb+ABZ*I_ERI_Fx2y_Py_Dxz_S_bb;
  Double I_ERI_Fxyz_Dyz_Dxz_S_bb = I_ERI_Gxy2z_Py_Dxz_S_bb+ABZ*I_ERI_Fxyz_Py_Dxz_S_bb;
  Double I_ERI_Fx2z_Dyz_Dxz_S_bb = I_ERI_Gx3z_Py_Dxz_S_bb+ABZ*I_ERI_Fx2z_Py_Dxz_S_bb;
  Double I_ERI_F3y_Dyz_Dxz_S_bb = I_ERI_G3yz_Py_Dxz_S_bb+ABZ*I_ERI_F3y_Py_Dxz_S_bb;
  Double I_ERI_F2yz_Dyz_Dxz_S_bb = I_ERI_G2y2z_Py_Dxz_S_bb+ABZ*I_ERI_F2yz_Py_Dxz_S_bb;
  Double I_ERI_Fy2z_Dyz_Dxz_S_bb = I_ERI_Gy3z_Py_Dxz_S_bb+ABZ*I_ERI_Fy2z_Py_Dxz_S_bb;
  Double I_ERI_F3z_Dyz_Dxz_S_bb = I_ERI_G4z_Py_Dxz_S_bb+ABZ*I_ERI_F3z_Py_Dxz_S_bb;
  Double I_ERI_F3x_D2z_Dxz_S_bb = I_ERI_G3xz_Pz_Dxz_S_bb+ABZ*I_ERI_F3x_Pz_Dxz_S_bb;
  Double I_ERI_F2xy_D2z_Dxz_S_bb = I_ERI_G2xyz_Pz_Dxz_S_bb+ABZ*I_ERI_F2xy_Pz_Dxz_S_bb;
  Double I_ERI_F2xz_D2z_Dxz_S_bb = I_ERI_G2x2z_Pz_Dxz_S_bb+ABZ*I_ERI_F2xz_Pz_Dxz_S_bb;
  Double I_ERI_Fx2y_D2z_Dxz_S_bb = I_ERI_Gx2yz_Pz_Dxz_S_bb+ABZ*I_ERI_Fx2y_Pz_Dxz_S_bb;
  Double I_ERI_Fxyz_D2z_Dxz_S_bb = I_ERI_Gxy2z_Pz_Dxz_S_bb+ABZ*I_ERI_Fxyz_Pz_Dxz_S_bb;
  Double I_ERI_Fx2z_D2z_Dxz_S_bb = I_ERI_Gx3z_Pz_Dxz_S_bb+ABZ*I_ERI_Fx2z_Pz_Dxz_S_bb;
  Double I_ERI_F3y_D2z_Dxz_S_bb = I_ERI_G3yz_Pz_Dxz_S_bb+ABZ*I_ERI_F3y_Pz_Dxz_S_bb;
  Double I_ERI_F2yz_D2z_Dxz_S_bb = I_ERI_G2y2z_Pz_Dxz_S_bb+ABZ*I_ERI_F2yz_Pz_Dxz_S_bb;
  Double I_ERI_Fy2z_D2z_Dxz_S_bb = I_ERI_Gy3z_Pz_Dxz_S_bb+ABZ*I_ERI_Fy2z_Pz_Dxz_S_bb;
  Double I_ERI_F3z_D2z_Dxz_S_bb = I_ERI_G4z_Pz_Dxz_S_bb+ABZ*I_ERI_F3z_Pz_Dxz_S_bb;
  Double I_ERI_F3x_D2x_D2y_S_bb = I_ERI_G4x_Px_D2y_S_bb+ABX*I_ERI_F3x_Px_D2y_S_bb;
  Double I_ERI_F2xy_D2x_D2y_S_bb = I_ERI_G3xy_Px_D2y_S_bb+ABX*I_ERI_F2xy_Px_D2y_S_bb;
  Double I_ERI_F2xz_D2x_D2y_S_bb = I_ERI_G3xz_Px_D2y_S_bb+ABX*I_ERI_F2xz_Px_D2y_S_bb;
  Double I_ERI_Fx2y_D2x_D2y_S_bb = I_ERI_G2x2y_Px_D2y_S_bb+ABX*I_ERI_Fx2y_Px_D2y_S_bb;
  Double I_ERI_Fxyz_D2x_D2y_S_bb = I_ERI_G2xyz_Px_D2y_S_bb+ABX*I_ERI_Fxyz_Px_D2y_S_bb;
  Double I_ERI_Fx2z_D2x_D2y_S_bb = I_ERI_G2x2z_Px_D2y_S_bb+ABX*I_ERI_Fx2z_Px_D2y_S_bb;
  Double I_ERI_F3y_D2x_D2y_S_bb = I_ERI_Gx3y_Px_D2y_S_bb+ABX*I_ERI_F3y_Px_D2y_S_bb;
  Double I_ERI_F2yz_D2x_D2y_S_bb = I_ERI_Gx2yz_Px_D2y_S_bb+ABX*I_ERI_F2yz_Px_D2y_S_bb;
  Double I_ERI_Fy2z_D2x_D2y_S_bb = I_ERI_Gxy2z_Px_D2y_S_bb+ABX*I_ERI_Fy2z_Px_D2y_S_bb;
  Double I_ERI_F3z_D2x_D2y_S_bb = I_ERI_Gx3z_Px_D2y_S_bb+ABX*I_ERI_F3z_Px_D2y_S_bb;
  Double I_ERI_F3x_Dxy_D2y_S_bb = I_ERI_G3xy_Px_D2y_S_bb+ABY*I_ERI_F3x_Px_D2y_S_bb;
  Double I_ERI_F2xy_Dxy_D2y_S_bb = I_ERI_G2x2y_Px_D2y_S_bb+ABY*I_ERI_F2xy_Px_D2y_S_bb;
  Double I_ERI_F2xz_Dxy_D2y_S_bb = I_ERI_G2xyz_Px_D2y_S_bb+ABY*I_ERI_F2xz_Px_D2y_S_bb;
  Double I_ERI_Fx2y_Dxy_D2y_S_bb = I_ERI_Gx3y_Px_D2y_S_bb+ABY*I_ERI_Fx2y_Px_D2y_S_bb;
  Double I_ERI_Fxyz_Dxy_D2y_S_bb = I_ERI_Gx2yz_Px_D2y_S_bb+ABY*I_ERI_Fxyz_Px_D2y_S_bb;
  Double I_ERI_Fx2z_Dxy_D2y_S_bb = I_ERI_Gxy2z_Px_D2y_S_bb+ABY*I_ERI_Fx2z_Px_D2y_S_bb;
  Double I_ERI_F3y_Dxy_D2y_S_bb = I_ERI_G4y_Px_D2y_S_bb+ABY*I_ERI_F3y_Px_D2y_S_bb;
  Double I_ERI_F2yz_Dxy_D2y_S_bb = I_ERI_G3yz_Px_D2y_S_bb+ABY*I_ERI_F2yz_Px_D2y_S_bb;
  Double I_ERI_Fy2z_Dxy_D2y_S_bb = I_ERI_G2y2z_Px_D2y_S_bb+ABY*I_ERI_Fy2z_Px_D2y_S_bb;
  Double I_ERI_F3z_Dxy_D2y_S_bb = I_ERI_Gy3z_Px_D2y_S_bb+ABY*I_ERI_F3z_Px_D2y_S_bb;
  Double I_ERI_F3x_Dxz_D2y_S_bb = I_ERI_G3xz_Px_D2y_S_bb+ABZ*I_ERI_F3x_Px_D2y_S_bb;
  Double I_ERI_F2xy_Dxz_D2y_S_bb = I_ERI_G2xyz_Px_D2y_S_bb+ABZ*I_ERI_F2xy_Px_D2y_S_bb;
  Double I_ERI_F2xz_Dxz_D2y_S_bb = I_ERI_G2x2z_Px_D2y_S_bb+ABZ*I_ERI_F2xz_Px_D2y_S_bb;
  Double I_ERI_Fx2y_Dxz_D2y_S_bb = I_ERI_Gx2yz_Px_D2y_S_bb+ABZ*I_ERI_Fx2y_Px_D2y_S_bb;
  Double I_ERI_Fxyz_Dxz_D2y_S_bb = I_ERI_Gxy2z_Px_D2y_S_bb+ABZ*I_ERI_Fxyz_Px_D2y_S_bb;
  Double I_ERI_Fx2z_Dxz_D2y_S_bb = I_ERI_Gx3z_Px_D2y_S_bb+ABZ*I_ERI_Fx2z_Px_D2y_S_bb;
  Double I_ERI_F3y_Dxz_D2y_S_bb = I_ERI_G3yz_Px_D2y_S_bb+ABZ*I_ERI_F3y_Px_D2y_S_bb;
  Double I_ERI_F2yz_Dxz_D2y_S_bb = I_ERI_G2y2z_Px_D2y_S_bb+ABZ*I_ERI_F2yz_Px_D2y_S_bb;
  Double I_ERI_Fy2z_Dxz_D2y_S_bb = I_ERI_Gy3z_Px_D2y_S_bb+ABZ*I_ERI_Fy2z_Px_D2y_S_bb;
  Double I_ERI_F3z_Dxz_D2y_S_bb = I_ERI_G4z_Px_D2y_S_bb+ABZ*I_ERI_F3z_Px_D2y_S_bb;
  Double I_ERI_F3x_D2y_D2y_S_bb = I_ERI_G3xy_Py_D2y_S_bb+ABY*I_ERI_F3x_Py_D2y_S_bb;
  Double I_ERI_F2xy_D2y_D2y_S_bb = I_ERI_G2x2y_Py_D2y_S_bb+ABY*I_ERI_F2xy_Py_D2y_S_bb;
  Double I_ERI_F2xz_D2y_D2y_S_bb = I_ERI_G2xyz_Py_D2y_S_bb+ABY*I_ERI_F2xz_Py_D2y_S_bb;
  Double I_ERI_Fx2y_D2y_D2y_S_bb = I_ERI_Gx3y_Py_D2y_S_bb+ABY*I_ERI_Fx2y_Py_D2y_S_bb;
  Double I_ERI_Fxyz_D2y_D2y_S_bb = I_ERI_Gx2yz_Py_D2y_S_bb+ABY*I_ERI_Fxyz_Py_D2y_S_bb;
  Double I_ERI_Fx2z_D2y_D2y_S_bb = I_ERI_Gxy2z_Py_D2y_S_bb+ABY*I_ERI_Fx2z_Py_D2y_S_bb;
  Double I_ERI_F3y_D2y_D2y_S_bb = I_ERI_G4y_Py_D2y_S_bb+ABY*I_ERI_F3y_Py_D2y_S_bb;
  Double I_ERI_F2yz_D2y_D2y_S_bb = I_ERI_G3yz_Py_D2y_S_bb+ABY*I_ERI_F2yz_Py_D2y_S_bb;
  Double I_ERI_Fy2z_D2y_D2y_S_bb = I_ERI_G2y2z_Py_D2y_S_bb+ABY*I_ERI_Fy2z_Py_D2y_S_bb;
  Double I_ERI_F3z_D2y_D2y_S_bb = I_ERI_Gy3z_Py_D2y_S_bb+ABY*I_ERI_F3z_Py_D2y_S_bb;
  Double I_ERI_F3x_Dyz_D2y_S_bb = I_ERI_G3xz_Py_D2y_S_bb+ABZ*I_ERI_F3x_Py_D2y_S_bb;
  Double I_ERI_F2xy_Dyz_D2y_S_bb = I_ERI_G2xyz_Py_D2y_S_bb+ABZ*I_ERI_F2xy_Py_D2y_S_bb;
  Double I_ERI_F2xz_Dyz_D2y_S_bb = I_ERI_G2x2z_Py_D2y_S_bb+ABZ*I_ERI_F2xz_Py_D2y_S_bb;
  Double I_ERI_Fx2y_Dyz_D2y_S_bb = I_ERI_Gx2yz_Py_D2y_S_bb+ABZ*I_ERI_Fx2y_Py_D2y_S_bb;
  Double I_ERI_Fxyz_Dyz_D2y_S_bb = I_ERI_Gxy2z_Py_D2y_S_bb+ABZ*I_ERI_Fxyz_Py_D2y_S_bb;
  Double I_ERI_Fx2z_Dyz_D2y_S_bb = I_ERI_Gx3z_Py_D2y_S_bb+ABZ*I_ERI_Fx2z_Py_D2y_S_bb;
  Double I_ERI_F3y_Dyz_D2y_S_bb = I_ERI_G3yz_Py_D2y_S_bb+ABZ*I_ERI_F3y_Py_D2y_S_bb;
  Double I_ERI_F2yz_Dyz_D2y_S_bb = I_ERI_G2y2z_Py_D2y_S_bb+ABZ*I_ERI_F2yz_Py_D2y_S_bb;
  Double I_ERI_Fy2z_Dyz_D2y_S_bb = I_ERI_Gy3z_Py_D2y_S_bb+ABZ*I_ERI_Fy2z_Py_D2y_S_bb;
  Double I_ERI_F3z_Dyz_D2y_S_bb = I_ERI_G4z_Py_D2y_S_bb+ABZ*I_ERI_F3z_Py_D2y_S_bb;
  Double I_ERI_F3x_D2z_D2y_S_bb = I_ERI_G3xz_Pz_D2y_S_bb+ABZ*I_ERI_F3x_Pz_D2y_S_bb;
  Double I_ERI_F2xy_D2z_D2y_S_bb = I_ERI_G2xyz_Pz_D2y_S_bb+ABZ*I_ERI_F2xy_Pz_D2y_S_bb;
  Double I_ERI_F2xz_D2z_D2y_S_bb = I_ERI_G2x2z_Pz_D2y_S_bb+ABZ*I_ERI_F2xz_Pz_D2y_S_bb;
  Double I_ERI_Fx2y_D2z_D2y_S_bb = I_ERI_Gx2yz_Pz_D2y_S_bb+ABZ*I_ERI_Fx2y_Pz_D2y_S_bb;
  Double I_ERI_Fxyz_D2z_D2y_S_bb = I_ERI_Gxy2z_Pz_D2y_S_bb+ABZ*I_ERI_Fxyz_Pz_D2y_S_bb;
  Double I_ERI_Fx2z_D2z_D2y_S_bb = I_ERI_Gx3z_Pz_D2y_S_bb+ABZ*I_ERI_Fx2z_Pz_D2y_S_bb;
  Double I_ERI_F3y_D2z_D2y_S_bb = I_ERI_G3yz_Pz_D2y_S_bb+ABZ*I_ERI_F3y_Pz_D2y_S_bb;
  Double I_ERI_F2yz_D2z_D2y_S_bb = I_ERI_G2y2z_Pz_D2y_S_bb+ABZ*I_ERI_F2yz_Pz_D2y_S_bb;
  Double I_ERI_Fy2z_D2z_D2y_S_bb = I_ERI_Gy3z_Pz_D2y_S_bb+ABZ*I_ERI_Fy2z_Pz_D2y_S_bb;
  Double I_ERI_F3z_D2z_D2y_S_bb = I_ERI_G4z_Pz_D2y_S_bb+ABZ*I_ERI_F3z_Pz_D2y_S_bb;
  Double I_ERI_F3x_D2x_Dyz_S_bb = I_ERI_G4x_Px_Dyz_S_bb+ABX*I_ERI_F3x_Px_Dyz_S_bb;
  Double I_ERI_F2xy_D2x_Dyz_S_bb = I_ERI_G3xy_Px_Dyz_S_bb+ABX*I_ERI_F2xy_Px_Dyz_S_bb;
  Double I_ERI_F2xz_D2x_Dyz_S_bb = I_ERI_G3xz_Px_Dyz_S_bb+ABX*I_ERI_F2xz_Px_Dyz_S_bb;
  Double I_ERI_Fx2y_D2x_Dyz_S_bb = I_ERI_G2x2y_Px_Dyz_S_bb+ABX*I_ERI_Fx2y_Px_Dyz_S_bb;
  Double I_ERI_Fxyz_D2x_Dyz_S_bb = I_ERI_G2xyz_Px_Dyz_S_bb+ABX*I_ERI_Fxyz_Px_Dyz_S_bb;
  Double I_ERI_Fx2z_D2x_Dyz_S_bb = I_ERI_G2x2z_Px_Dyz_S_bb+ABX*I_ERI_Fx2z_Px_Dyz_S_bb;
  Double I_ERI_F3y_D2x_Dyz_S_bb = I_ERI_Gx3y_Px_Dyz_S_bb+ABX*I_ERI_F3y_Px_Dyz_S_bb;
  Double I_ERI_F2yz_D2x_Dyz_S_bb = I_ERI_Gx2yz_Px_Dyz_S_bb+ABX*I_ERI_F2yz_Px_Dyz_S_bb;
  Double I_ERI_Fy2z_D2x_Dyz_S_bb = I_ERI_Gxy2z_Px_Dyz_S_bb+ABX*I_ERI_Fy2z_Px_Dyz_S_bb;
  Double I_ERI_F3z_D2x_Dyz_S_bb = I_ERI_Gx3z_Px_Dyz_S_bb+ABX*I_ERI_F3z_Px_Dyz_S_bb;
  Double I_ERI_F3x_Dxy_Dyz_S_bb = I_ERI_G3xy_Px_Dyz_S_bb+ABY*I_ERI_F3x_Px_Dyz_S_bb;
  Double I_ERI_F2xy_Dxy_Dyz_S_bb = I_ERI_G2x2y_Px_Dyz_S_bb+ABY*I_ERI_F2xy_Px_Dyz_S_bb;
  Double I_ERI_F2xz_Dxy_Dyz_S_bb = I_ERI_G2xyz_Px_Dyz_S_bb+ABY*I_ERI_F2xz_Px_Dyz_S_bb;
  Double I_ERI_Fx2y_Dxy_Dyz_S_bb = I_ERI_Gx3y_Px_Dyz_S_bb+ABY*I_ERI_Fx2y_Px_Dyz_S_bb;
  Double I_ERI_Fxyz_Dxy_Dyz_S_bb = I_ERI_Gx2yz_Px_Dyz_S_bb+ABY*I_ERI_Fxyz_Px_Dyz_S_bb;
  Double I_ERI_Fx2z_Dxy_Dyz_S_bb = I_ERI_Gxy2z_Px_Dyz_S_bb+ABY*I_ERI_Fx2z_Px_Dyz_S_bb;
  Double I_ERI_F3y_Dxy_Dyz_S_bb = I_ERI_G4y_Px_Dyz_S_bb+ABY*I_ERI_F3y_Px_Dyz_S_bb;
  Double I_ERI_F2yz_Dxy_Dyz_S_bb = I_ERI_G3yz_Px_Dyz_S_bb+ABY*I_ERI_F2yz_Px_Dyz_S_bb;
  Double I_ERI_Fy2z_Dxy_Dyz_S_bb = I_ERI_G2y2z_Px_Dyz_S_bb+ABY*I_ERI_Fy2z_Px_Dyz_S_bb;
  Double I_ERI_F3z_Dxy_Dyz_S_bb = I_ERI_Gy3z_Px_Dyz_S_bb+ABY*I_ERI_F3z_Px_Dyz_S_bb;
  Double I_ERI_F3x_Dxz_Dyz_S_bb = I_ERI_G3xz_Px_Dyz_S_bb+ABZ*I_ERI_F3x_Px_Dyz_S_bb;
  Double I_ERI_F2xy_Dxz_Dyz_S_bb = I_ERI_G2xyz_Px_Dyz_S_bb+ABZ*I_ERI_F2xy_Px_Dyz_S_bb;
  Double I_ERI_F2xz_Dxz_Dyz_S_bb = I_ERI_G2x2z_Px_Dyz_S_bb+ABZ*I_ERI_F2xz_Px_Dyz_S_bb;
  Double I_ERI_Fx2y_Dxz_Dyz_S_bb = I_ERI_Gx2yz_Px_Dyz_S_bb+ABZ*I_ERI_Fx2y_Px_Dyz_S_bb;
  Double I_ERI_Fxyz_Dxz_Dyz_S_bb = I_ERI_Gxy2z_Px_Dyz_S_bb+ABZ*I_ERI_Fxyz_Px_Dyz_S_bb;
  Double I_ERI_Fx2z_Dxz_Dyz_S_bb = I_ERI_Gx3z_Px_Dyz_S_bb+ABZ*I_ERI_Fx2z_Px_Dyz_S_bb;
  Double I_ERI_F3y_Dxz_Dyz_S_bb = I_ERI_G3yz_Px_Dyz_S_bb+ABZ*I_ERI_F3y_Px_Dyz_S_bb;
  Double I_ERI_F2yz_Dxz_Dyz_S_bb = I_ERI_G2y2z_Px_Dyz_S_bb+ABZ*I_ERI_F2yz_Px_Dyz_S_bb;
  Double I_ERI_Fy2z_Dxz_Dyz_S_bb = I_ERI_Gy3z_Px_Dyz_S_bb+ABZ*I_ERI_Fy2z_Px_Dyz_S_bb;
  Double I_ERI_F3z_Dxz_Dyz_S_bb = I_ERI_G4z_Px_Dyz_S_bb+ABZ*I_ERI_F3z_Px_Dyz_S_bb;
  Double I_ERI_F3x_D2y_Dyz_S_bb = I_ERI_G3xy_Py_Dyz_S_bb+ABY*I_ERI_F3x_Py_Dyz_S_bb;
  Double I_ERI_F2xy_D2y_Dyz_S_bb = I_ERI_G2x2y_Py_Dyz_S_bb+ABY*I_ERI_F2xy_Py_Dyz_S_bb;
  Double I_ERI_F2xz_D2y_Dyz_S_bb = I_ERI_G2xyz_Py_Dyz_S_bb+ABY*I_ERI_F2xz_Py_Dyz_S_bb;
  Double I_ERI_Fx2y_D2y_Dyz_S_bb = I_ERI_Gx3y_Py_Dyz_S_bb+ABY*I_ERI_Fx2y_Py_Dyz_S_bb;
  Double I_ERI_Fxyz_D2y_Dyz_S_bb = I_ERI_Gx2yz_Py_Dyz_S_bb+ABY*I_ERI_Fxyz_Py_Dyz_S_bb;
  Double I_ERI_Fx2z_D2y_Dyz_S_bb = I_ERI_Gxy2z_Py_Dyz_S_bb+ABY*I_ERI_Fx2z_Py_Dyz_S_bb;
  Double I_ERI_F3y_D2y_Dyz_S_bb = I_ERI_G4y_Py_Dyz_S_bb+ABY*I_ERI_F3y_Py_Dyz_S_bb;
  Double I_ERI_F2yz_D2y_Dyz_S_bb = I_ERI_G3yz_Py_Dyz_S_bb+ABY*I_ERI_F2yz_Py_Dyz_S_bb;
  Double I_ERI_Fy2z_D2y_Dyz_S_bb = I_ERI_G2y2z_Py_Dyz_S_bb+ABY*I_ERI_Fy2z_Py_Dyz_S_bb;
  Double I_ERI_F3z_D2y_Dyz_S_bb = I_ERI_Gy3z_Py_Dyz_S_bb+ABY*I_ERI_F3z_Py_Dyz_S_bb;
  Double I_ERI_F3x_Dyz_Dyz_S_bb = I_ERI_G3xz_Py_Dyz_S_bb+ABZ*I_ERI_F3x_Py_Dyz_S_bb;
  Double I_ERI_F2xy_Dyz_Dyz_S_bb = I_ERI_G2xyz_Py_Dyz_S_bb+ABZ*I_ERI_F2xy_Py_Dyz_S_bb;
  Double I_ERI_F2xz_Dyz_Dyz_S_bb = I_ERI_G2x2z_Py_Dyz_S_bb+ABZ*I_ERI_F2xz_Py_Dyz_S_bb;
  Double I_ERI_Fx2y_Dyz_Dyz_S_bb = I_ERI_Gx2yz_Py_Dyz_S_bb+ABZ*I_ERI_Fx2y_Py_Dyz_S_bb;
  Double I_ERI_Fxyz_Dyz_Dyz_S_bb = I_ERI_Gxy2z_Py_Dyz_S_bb+ABZ*I_ERI_Fxyz_Py_Dyz_S_bb;
  Double I_ERI_Fx2z_Dyz_Dyz_S_bb = I_ERI_Gx3z_Py_Dyz_S_bb+ABZ*I_ERI_Fx2z_Py_Dyz_S_bb;
  Double I_ERI_F3y_Dyz_Dyz_S_bb = I_ERI_G3yz_Py_Dyz_S_bb+ABZ*I_ERI_F3y_Py_Dyz_S_bb;
  Double I_ERI_F2yz_Dyz_Dyz_S_bb = I_ERI_G2y2z_Py_Dyz_S_bb+ABZ*I_ERI_F2yz_Py_Dyz_S_bb;
  Double I_ERI_Fy2z_Dyz_Dyz_S_bb = I_ERI_Gy3z_Py_Dyz_S_bb+ABZ*I_ERI_Fy2z_Py_Dyz_S_bb;
  Double I_ERI_F3z_Dyz_Dyz_S_bb = I_ERI_G4z_Py_Dyz_S_bb+ABZ*I_ERI_F3z_Py_Dyz_S_bb;
  Double I_ERI_F3x_D2z_Dyz_S_bb = I_ERI_G3xz_Pz_Dyz_S_bb+ABZ*I_ERI_F3x_Pz_Dyz_S_bb;
  Double I_ERI_F2xy_D2z_Dyz_S_bb = I_ERI_G2xyz_Pz_Dyz_S_bb+ABZ*I_ERI_F2xy_Pz_Dyz_S_bb;
  Double I_ERI_F2xz_D2z_Dyz_S_bb = I_ERI_G2x2z_Pz_Dyz_S_bb+ABZ*I_ERI_F2xz_Pz_Dyz_S_bb;
  Double I_ERI_Fx2y_D2z_Dyz_S_bb = I_ERI_Gx2yz_Pz_Dyz_S_bb+ABZ*I_ERI_Fx2y_Pz_Dyz_S_bb;
  Double I_ERI_Fxyz_D2z_Dyz_S_bb = I_ERI_Gxy2z_Pz_Dyz_S_bb+ABZ*I_ERI_Fxyz_Pz_Dyz_S_bb;
  Double I_ERI_Fx2z_D2z_Dyz_S_bb = I_ERI_Gx3z_Pz_Dyz_S_bb+ABZ*I_ERI_Fx2z_Pz_Dyz_S_bb;
  Double I_ERI_F3y_D2z_Dyz_S_bb = I_ERI_G3yz_Pz_Dyz_S_bb+ABZ*I_ERI_F3y_Pz_Dyz_S_bb;
  Double I_ERI_F2yz_D2z_Dyz_S_bb = I_ERI_G2y2z_Pz_Dyz_S_bb+ABZ*I_ERI_F2yz_Pz_Dyz_S_bb;
  Double I_ERI_Fy2z_D2z_Dyz_S_bb = I_ERI_Gy3z_Pz_Dyz_S_bb+ABZ*I_ERI_Fy2z_Pz_Dyz_S_bb;
  Double I_ERI_F3z_D2z_Dyz_S_bb = I_ERI_G4z_Pz_Dyz_S_bb+ABZ*I_ERI_F3z_Pz_Dyz_S_bb;
  Double I_ERI_F3x_D2x_D2z_S_bb = I_ERI_G4x_Px_D2z_S_bb+ABX*I_ERI_F3x_Px_D2z_S_bb;
  Double I_ERI_F2xy_D2x_D2z_S_bb = I_ERI_G3xy_Px_D2z_S_bb+ABX*I_ERI_F2xy_Px_D2z_S_bb;
  Double I_ERI_F2xz_D2x_D2z_S_bb = I_ERI_G3xz_Px_D2z_S_bb+ABX*I_ERI_F2xz_Px_D2z_S_bb;
  Double I_ERI_Fx2y_D2x_D2z_S_bb = I_ERI_G2x2y_Px_D2z_S_bb+ABX*I_ERI_Fx2y_Px_D2z_S_bb;
  Double I_ERI_Fxyz_D2x_D2z_S_bb = I_ERI_G2xyz_Px_D2z_S_bb+ABX*I_ERI_Fxyz_Px_D2z_S_bb;
  Double I_ERI_Fx2z_D2x_D2z_S_bb = I_ERI_G2x2z_Px_D2z_S_bb+ABX*I_ERI_Fx2z_Px_D2z_S_bb;
  Double I_ERI_F3y_D2x_D2z_S_bb = I_ERI_Gx3y_Px_D2z_S_bb+ABX*I_ERI_F3y_Px_D2z_S_bb;
  Double I_ERI_F2yz_D2x_D2z_S_bb = I_ERI_Gx2yz_Px_D2z_S_bb+ABX*I_ERI_F2yz_Px_D2z_S_bb;
  Double I_ERI_Fy2z_D2x_D2z_S_bb = I_ERI_Gxy2z_Px_D2z_S_bb+ABX*I_ERI_Fy2z_Px_D2z_S_bb;
  Double I_ERI_F3z_D2x_D2z_S_bb = I_ERI_Gx3z_Px_D2z_S_bb+ABX*I_ERI_F3z_Px_D2z_S_bb;
  Double I_ERI_F3x_Dxy_D2z_S_bb = I_ERI_G3xy_Px_D2z_S_bb+ABY*I_ERI_F3x_Px_D2z_S_bb;
  Double I_ERI_F2xy_Dxy_D2z_S_bb = I_ERI_G2x2y_Px_D2z_S_bb+ABY*I_ERI_F2xy_Px_D2z_S_bb;
  Double I_ERI_F2xz_Dxy_D2z_S_bb = I_ERI_G2xyz_Px_D2z_S_bb+ABY*I_ERI_F2xz_Px_D2z_S_bb;
  Double I_ERI_Fx2y_Dxy_D2z_S_bb = I_ERI_Gx3y_Px_D2z_S_bb+ABY*I_ERI_Fx2y_Px_D2z_S_bb;
  Double I_ERI_Fxyz_Dxy_D2z_S_bb = I_ERI_Gx2yz_Px_D2z_S_bb+ABY*I_ERI_Fxyz_Px_D2z_S_bb;
  Double I_ERI_Fx2z_Dxy_D2z_S_bb = I_ERI_Gxy2z_Px_D2z_S_bb+ABY*I_ERI_Fx2z_Px_D2z_S_bb;
  Double I_ERI_F3y_Dxy_D2z_S_bb = I_ERI_G4y_Px_D2z_S_bb+ABY*I_ERI_F3y_Px_D2z_S_bb;
  Double I_ERI_F2yz_Dxy_D2z_S_bb = I_ERI_G3yz_Px_D2z_S_bb+ABY*I_ERI_F2yz_Px_D2z_S_bb;
  Double I_ERI_Fy2z_Dxy_D2z_S_bb = I_ERI_G2y2z_Px_D2z_S_bb+ABY*I_ERI_Fy2z_Px_D2z_S_bb;
  Double I_ERI_F3z_Dxy_D2z_S_bb = I_ERI_Gy3z_Px_D2z_S_bb+ABY*I_ERI_F3z_Px_D2z_S_bb;
  Double I_ERI_F3x_Dxz_D2z_S_bb = I_ERI_G3xz_Px_D2z_S_bb+ABZ*I_ERI_F3x_Px_D2z_S_bb;
  Double I_ERI_F2xy_Dxz_D2z_S_bb = I_ERI_G2xyz_Px_D2z_S_bb+ABZ*I_ERI_F2xy_Px_D2z_S_bb;
  Double I_ERI_F2xz_Dxz_D2z_S_bb = I_ERI_G2x2z_Px_D2z_S_bb+ABZ*I_ERI_F2xz_Px_D2z_S_bb;
  Double I_ERI_Fx2y_Dxz_D2z_S_bb = I_ERI_Gx2yz_Px_D2z_S_bb+ABZ*I_ERI_Fx2y_Px_D2z_S_bb;
  Double I_ERI_Fxyz_Dxz_D2z_S_bb = I_ERI_Gxy2z_Px_D2z_S_bb+ABZ*I_ERI_Fxyz_Px_D2z_S_bb;
  Double I_ERI_Fx2z_Dxz_D2z_S_bb = I_ERI_Gx3z_Px_D2z_S_bb+ABZ*I_ERI_Fx2z_Px_D2z_S_bb;
  Double I_ERI_F3y_Dxz_D2z_S_bb = I_ERI_G3yz_Px_D2z_S_bb+ABZ*I_ERI_F3y_Px_D2z_S_bb;
  Double I_ERI_F2yz_Dxz_D2z_S_bb = I_ERI_G2y2z_Px_D2z_S_bb+ABZ*I_ERI_F2yz_Px_D2z_S_bb;
  Double I_ERI_Fy2z_Dxz_D2z_S_bb = I_ERI_Gy3z_Px_D2z_S_bb+ABZ*I_ERI_Fy2z_Px_D2z_S_bb;
  Double I_ERI_F3z_Dxz_D2z_S_bb = I_ERI_G4z_Px_D2z_S_bb+ABZ*I_ERI_F3z_Px_D2z_S_bb;
  Double I_ERI_F3x_D2y_D2z_S_bb = I_ERI_G3xy_Py_D2z_S_bb+ABY*I_ERI_F3x_Py_D2z_S_bb;
  Double I_ERI_F2xy_D2y_D2z_S_bb = I_ERI_G2x2y_Py_D2z_S_bb+ABY*I_ERI_F2xy_Py_D2z_S_bb;
  Double I_ERI_F2xz_D2y_D2z_S_bb = I_ERI_G2xyz_Py_D2z_S_bb+ABY*I_ERI_F2xz_Py_D2z_S_bb;
  Double I_ERI_Fx2y_D2y_D2z_S_bb = I_ERI_Gx3y_Py_D2z_S_bb+ABY*I_ERI_Fx2y_Py_D2z_S_bb;
  Double I_ERI_Fxyz_D2y_D2z_S_bb = I_ERI_Gx2yz_Py_D2z_S_bb+ABY*I_ERI_Fxyz_Py_D2z_S_bb;
  Double I_ERI_Fx2z_D2y_D2z_S_bb = I_ERI_Gxy2z_Py_D2z_S_bb+ABY*I_ERI_Fx2z_Py_D2z_S_bb;
  Double I_ERI_F3y_D2y_D2z_S_bb = I_ERI_G4y_Py_D2z_S_bb+ABY*I_ERI_F3y_Py_D2z_S_bb;
  Double I_ERI_F2yz_D2y_D2z_S_bb = I_ERI_G3yz_Py_D2z_S_bb+ABY*I_ERI_F2yz_Py_D2z_S_bb;
  Double I_ERI_Fy2z_D2y_D2z_S_bb = I_ERI_G2y2z_Py_D2z_S_bb+ABY*I_ERI_Fy2z_Py_D2z_S_bb;
  Double I_ERI_F3z_D2y_D2z_S_bb = I_ERI_Gy3z_Py_D2z_S_bb+ABY*I_ERI_F3z_Py_D2z_S_bb;
  Double I_ERI_F3x_Dyz_D2z_S_bb = I_ERI_G3xz_Py_D2z_S_bb+ABZ*I_ERI_F3x_Py_D2z_S_bb;
  Double I_ERI_F2xy_Dyz_D2z_S_bb = I_ERI_G2xyz_Py_D2z_S_bb+ABZ*I_ERI_F2xy_Py_D2z_S_bb;
  Double I_ERI_F2xz_Dyz_D2z_S_bb = I_ERI_G2x2z_Py_D2z_S_bb+ABZ*I_ERI_F2xz_Py_D2z_S_bb;
  Double I_ERI_Fx2y_Dyz_D2z_S_bb = I_ERI_Gx2yz_Py_D2z_S_bb+ABZ*I_ERI_Fx2y_Py_D2z_S_bb;
  Double I_ERI_Fxyz_Dyz_D2z_S_bb = I_ERI_Gxy2z_Py_D2z_S_bb+ABZ*I_ERI_Fxyz_Py_D2z_S_bb;
  Double I_ERI_Fx2z_Dyz_D2z_S_bb = I_ERI_Gx3z_Py_D2z_S_bb+ABZ*I_ERI_Fx2z_Py_D2z_S_bb;
  Double I_ERI_F3y_Dyz_D2z_S_bb = I_ERI_G3yz_Py_D2z_S_bb+ABZ*I_ERI_F3y_Py_D2z_S_bb;
  Double I_ERI_F2yz_Dyz_D2z_S_bb = I_ERI_G2y2z_Py_D2z_S_bb+ABZ*I_ERI_F2yz_Py_D2z_S_bb;
  Double I_ERI_Fy2z_Dyz_D2z_S_bb = I_ERI_Gy3z_Py_D2z_S_bb+ABZ*I_ERI_Fy2z_Py_D2z_S_bb;
  Double I_ERI_F3z_Dyz_D2z_S_bb = I_ERI_G4z_Py_D2z_S_bb+ABZ*I_ERI_F3z_Py_D2z_S_bb;
  Double I_ERI_F3x_D2z_D2z_S_bb = I_ERI_G3xz_Pz_D2z_S_bb+ABZ*I_ERI_F3x_Pz_D2z_S_bb;
  Double I_ERI_F2xy_D2z_D2z_S_bb = I_ERI_G2xyz_Pz_D2z_S_bb+ABZ*I_ERI_F2xy_Pz_D2z_S_bb;
  Double I_ERI_F2xz_D2z_D2z_S_bb = I_ERI_G2x2z_Pz_D2z_S_bb+ABZ*I_ERI_F2xz_Pz_D2z_S_bb;
  Double I_ERI_Fx2y_D2z_D2z_S_bb = I_ERI_Gx2yz_Pz_D2z_S_bb+ABZ*I_ERI_Fx2y_Pz_D2z_S_bb;
  Double I_ERI_Fxyz_D2z_D2z_S_bb = I_ERI_Gxy2z_Pz_D2z_S_bb+ABZ*I_ERI_Fxyz_Pz_D2z_S_bb;
  Double I_ERI_Fx2z_D2z_D2z_S_bb = I_ERI_Gx3z_Pz_D2z_S_bb+ABZ*I_ERI_Fx2z_Pz_D2z_S_bb;
  Double I_ERI_F3y_D2z_D2z_S_bb = I_ERI_G3yz_Pz_D2z_S_bb+ABZ*I_ERI_F3y_Pz_D2z_S_bb;
  Double I_ERI_F2yz_D2z_D2z_S_bb = I_ERI_G2y2z_Pz_D2z_S_bb+ABZ*I_ERI_F2yz_Pz_D2z_S_bb;
  Double I_ERI_Fy2z_D2z_D2z_S_bb = I_ERI_Gy3z_Pz_D2z_S_bb+ABZ*I_ERI_Fy2z_Pz_D2z_S_bb;
  Double I_ERI_F3z_D2z_D2z_S_bb = I_ERI_G4z_Pz_D2z_S_bb+ABZ*I_ERI_F3z_Pz_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_F_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_S_F_S_bc
   ************************************************************/
  Double I_ERI_F3x_Px_F3x_S_bc = I_ERI_G4x_S_F3x_S_bc+ABX*I_ERI_F3x_S_F3x_S_bc;
  Double I_ERI_F2xy_Px_F3x_S_bc = I_ERI_G3xy_S_F3x_S_bc+ABX*I_ERI_F2xy_S_F3x_S_bc;
  Double I_ERI_F2xz_Px_F3x_S_bc = I_ERI_G3xz_S_F3x_S_bc+ABX*I_ERI_F2xz_S_F3x_S_bc;
  Double I_ERI_Fx2y_Px_F3x_S_bc = I_ERI_G2x2y_S_F3x_S_bc+ABX*I_ERI_Fx2y_S_F3x_S_bc;
  Double I_ERI_Fxyz_Px_F3x_S_bc = I_ERI_G2xyz_S_F3x_S_bc+ABX*I_ERI_Fxyz_S_F3x_S_bc;
  Double I_ERI_Fx2z_Px_F3x_S_bc = I_ERI_G2x2z_S_F3x_S_bc+ABX*I_ERI_Fx2z_S_F3x_S_bc;
  Double I_ERI_F3y_Px_F3x_S_bc = I_ERI_Gx3y_S_F3x_S_bc+ABX*I_ERI_F3y_S_F3x_S_bc;
  Double I_ERI_F2yz_Px_F3x_S_bc = I_ERI_Gx2yz_S_F3x_S_bc+ABX*I_ERI_F2yz_S_F3x_S_bc;
  Double I_ERI_Fy2z_Px_F3x_S_bc = I_ERI_Gxy2z_S_F3x_S_bc+ABX*I_ERI_Fy2z_S_F3x_S_bc;
  Double I_ERI_F3z_Px_F3x_S_bc = I_ERI_Gx3z_S_F3x_S_bc+ABX*I_ERI_F3z_S_F3x_S_bc;
  Double I_ERI_F3x_Py_F3x_S_bc = I_ERI_G3xy_S_F3x_S_bc+ABY*I_ERI_F3x_S_F3x_S_bc;
  Double I_ERI_F2xy_Py_F3x_S_bc = I_ERI_G2x2y_S_F3x_S_bc+ABY*I_ERI_F2xy_S_F3x_S_bc;
  Double I_ERI_F2xz_Py_F3x_S_bc = I_ERI_G2xyz_S_F3x_S_bc+ABY*I_ERI_F2xz_S_F3x_S_bc;
  Double I_ERI_Fx2y_Py_F3x_S_bc = I_ERI_Gx3y_S_F3x_S_bc+ABY*I_ERI_Fx2y_S_F3x_S_bc;
  Double I_ERI_Fxyz_Py_F3x_S_bc = I_ERI_Gx2yz_S_F3x_S_bc+ABY*I_ERI_Fxyz_S_F3x_S_bc;
  Double I_ERI_Fx2z_Py_F3x_S_bc = I_ERI_Gxy2z_S_F3x_S_bc+ABY*I_ERI_Fx2z_S_F3x_S_bc;
  Double I_ERI_F3y_Py_F3x_S_bc = I_ERI_G4y_S_F3x_S_bc+ABY*I_ERI_F3y_S_F3x_S_bc;
  Double I_ERI_F2yz_Py_F3x_S_bc = I_ERI_G3yz_S_F3x_S_bc+ABY*I_ERI_F2yz_S_F3x_S_bc;
  Double I_ERI_Fy2z_Py_F3x_S_bc = I_ERI_G2y2z_S_F3x_S_bc+ABY*I_ERI_Fy2z_S_F3x_S_bc;
  Double I_ERI_F3z_Py_F3x_S_bc = I_ERI_Gy3z_S_F3x_S_bc+ABY*I_ERI_F3z_S_F3x_S_bc;
  Double I_ERI_F3x_Pz_F3x_S_bc = I_ERI_G3xz_S_F3x_S_bc+ABZ*I_ERI_F3x_S_F3x_S_bc;
  Double I_ERI_F2xy_Pz_F3x_S_bc = I_ERI_G2xyz_S_F3x_S_bc+ABZ*I_ERI_F2xy_S_F3x_S_bc;
  Double I_ERI_F2xz_Pz_F3x_S_bc = I_ERI_G2x2z_S_F3x_S_bc+ABZ*I_ERI_F2xz_S_F3x_S_bc;
  Double I_ERI_Fx2y_Pz_F3x_S_bc = I_ERI_Gx2yz_S_F3x_S_bc+ABZ*I_ERI_Fx2y_S_F3x_S_bc;
  Double I_ERI_Fxyz_Pz_F3x_S_bc = I_ERI_Gxy2z_S_F3x_S_bc+ABZ*I_ERI_Fxyz_S_F3x_S_bc;
  Double I_ERI_Fx2z_Pz_F3x_S_bc = I_ERI_Gx3z_S_F3x_S_bc+ABZ*I_ERI_Fx2z_S_F3x_S_bc;
  Double I_ERI_F3y_Pz_F3x_S_bc = I_ERI_G3yz_S_F3x_S_bc+ABZ*I_ERI_F3y_S_F3x_S_bc;
  Double I_ERI_F2yz_Pz_F3x_S_bc = I_ERI_G2y2z_S_F3x_S_bc+ABZ*I_ERI_F2yz_S_F3x_S_bc;
  Double I_ERI_Fy2z_Pz_F3x_S_bc = I_ERI_Gy3z_S_F3x_S_bc+ABZ*I_ERI_Fy2z_S_F3x_S_bc;
  Double I_ERI_F3z_Pz_F3x_S_bc = I_ERI_G4z_S_F3x_S_bc+ABZ*I_ERI_F3z_S_F3x_S_bc;
  Double I_ERI_F3x_Px_F2xy_S_bc = I_ERI_G4x_S_F2xy_S_bc+ABX*I_ERI_F3x_S_F2xy_S_bc;
  Double I_ERI_F2xy_Px_F2xy_S_bc = I_ERI_G3xy_S_F2xy_S_bc+ABX*I_ERI_F2xy_S_F2xy_S_bc;
  Double I_ERI_F2xz_Px_F2xy_S_bc = I_ERI_G3xz_S_F2xy_S_bc+ABX*I_ERI_F2xz_S_F2xy_S_bc;
  Double I_ERI_Fx2y_Px_F2xy_S_bc = I_ERI_G2x2y_S_F2xy_S_bc+ABX*I_ERI_Fx2y_S_F2xy_S_bc;
  Double I_ERI_Fxyz_Px_F2xy_S_bc = I_ERI_G2xyz_S_F2xy_S_bc+ABX*I_ERI_Fxyz_S_F2xy_S_bc;
  Double I_ERI_Fx2z_Px_F2xy_S_bc = I_ERI_G2x2z_S_F2xy_S_bc+ABX*I_ERI_Fx2z_S_F2xy_S_bc;
  Double I_ERI_F3y_Px_F2xy_S_bc = I_ERI_Gx3y_S_F2xy_S_bc+ABX*I_ERI_F3y_S_F2xy_S_bc;
  Double I_ERI_F2yz_Px_F2xy_S_bc = I_ERI_Gx2yz_S_F2xy_S_bc+ABX*I_ERI_F2yz_S_F2xy_S_bc;
  Double I_ERI_Fy2z_Px_F2xy_S_bc = I_ERI_Gxy2z_S_F2xy_S_bc+ABX*I_ERI_Fy2z_S_F2xy_S_bc;
  Double I_ERI_F3z_Px_F2xy_S_bc = I_ERI_Gx3z_S_F2xy_S_bc+ABX*I_ERI_F3z_S_F2xy_S_bc;
  Double I_ERI_F3x_Py_F2xy_S_bc = I_ERI_G3xy_S_F2xy_S_bc+ABY*I_ERI_F3x_S_F2xy_S_bc;
  Double I_ERI_F2xy_Py_F2xy_S_bc = I_ERI_G2x2y_S_F2xy_S_bc+ABY*I_ERI_F2xy_S_F2xy_S_bc;
  Double I_ERI_F2xz_Py_F2xy_S_bc = I_ERI_G2xyz_S_F2xy_S_bc+ABY*I_ERI_F2xz_S_F2xy_S_bc;
  Double I_ERI_Fx2y_Py_F2xy_S_bc = I_ERI_Gx3y_S_F2xy_S_bc+ABY*I_ERI_Fx2y_S_F2xy_S_bc;
  Double I_ERI_Fxyz_Py_F2xy_S_bc = I_ERI_Gx2yz_S_F2xy_S_bc+ABY*I_ERI_Fxyz_S_F2xy_S_bc;
  Double I_ERI_Fx2z_Py_F2xy_S_bc = I_ERI_Gxy2z_S_F2xy_S_bc+ABY*I_ERI_Fx2z_S_F2xy_S_bc;
  Double I_ERI_F3y_Py_F2xy_S_bc = I_ERI_G4y_S_F2xy_S_bc+ABY*I_ERI_F3y_S_F2xy_S_bc;
  Double I_ERI_F2yz_Py_F2xy_S_bc = I_ERI_G3yz_S_F2xy_S_bc+ABY*I_ERI_F2yz_S_F2xy_S_bc;
  Double I_ERI_Fy2z_Py_F2xy_S_bc = I_ERI_G2y2z_S_F2xy_S_bc+ABY*I_ERI_Fy2z_S_F2xy_S_bc;
  Double I_ERI_F3z_Py_F2xy_S_bc = I_ERI_Gy3z_S_F2xy_S_bc+ABY*I_ERI_F3z_S_F2xy_S_bc;
  Double I_ERI_F3x_Pz_F2xy_S_bc = I_ERI_G3xz_S_F2xy_S_bc+ABZ*I_ERI_F3x_S_F2xy_S_bc;
  Double I_ERI_F2xy_Pz_F2xy_S_bc = I_ERI_G2xyz_S_F2xy_S_bc+ABZ*I_ERI_F2xy_S_F2xy_S_bc;
  Double I_ERI_F2xz_Pz_F2xy_S_bc = I_ERI_G2x2z_S_F2xy_S_bc+ABZ*I_ERI_F2xz_S_F2xy_S_bc;
  Double I_ERI_Fx2y_Pz_F2xy_S_bc = I_ERI_Gx2yz_S_F2xy_S_bc+ABZ*I_ERI_Fx2y_S_F2xy_S_bc;
  Double I_ERI_Fxyz_Pz_F2xy_S_bc = I_ERI_Gxy2z_S_F2xy_S_bc+ABZ*I_ERI_Fxyz_S_F2xy_S_bc;
  Double I_ERI_Fx2z_Pz_F2xy_S_bc = I_ERI_Gx3z_S_F2xy_S_bc+ABZ*I_ERI_Fx2z_S_F2xy_S_bc;
  Double I_ERI_F3y_Pz_F2xy_S_bc = I_ERI_G3yz_S_F2xy_S_bc+ABZ*I_ERI_F3y_S_F2xy_S_bc;
  Double I_ERI_F2yz_Pz_F2xy_S_bc = I_ERI_G2y2z_S_F2xy_S_bc+ABZ*I_ERI_F2yz_S_F2xy_S_bc;
  Double I_ERI_Fy2z_Pz_F2xy_S_bc = I_ERI_Gy3z_S_F2xy_S_bc+ABZ*I_ERI_Fy2z_S_F2xy_S_bc;
  Double I_ERI_F3z_Pz_F2xy_S_bc = I_ERI_G4z_S_F2xy_S_bc+ABZ*I_ERI_F3z_S_F2xy_S_bc;
  Double I_ERI_F3x_Px_F2xz_S_bc = I_ERI_G4x_S_F2xz_S_bc+ABX*I_ERI_F3x_S_F2xz_S_bc;
  Double I_ERI_F2xy_Px_F2xz_S_bc = I_ERI_G3xy_S_F2xz_S_bc+ABX*I_ERI_F2xy_S_F2xz_S_bc;
  Double I_ERI_F2xz_Px_F2xz_S_bc = I_ERI_G3xz_S_F2xz_S_bc+ABX*I_ERI_F2xz_S_F2xz_S_bc;
  Double I_ERI_Fx2y_Px_F2xz_S_bc = I_ERI_G2x2y_S_F2xz_S_bc+ABX*I_ERI_Fx2y_S_F2xz_S_bc;
  Double I_ERI_Fxyz_Px_F2xz_S_bc = I_ERI_G2xyz_S_F2xz_S_bc+ABX*I_ERI_Fxyz_S_F2xz_S_bc;
  Double I_ERI_Fx2z_Px_F2xz_S_bc = I_ERI_G2x2z_S_F2xz_S_bc+ABX*I_ERI_Fx2z_S_F2xz_S_bc;
  Double I_ERI_F3y_Px_F2xz_S_bc = I_ERI_Gx3y_S_F2xz_S_bc+ABX*I_ERI_F3y_S_F2xz_S_bc;
  Double I_ERI_F2yz_Px_F2xz_S_bc = I_ERI_Gx2yz_S_F2xz_S_bc+ABX*I_ERI_F2yz_S_F2xz_S_bc;
  Double I_ERI_Fy2z_Px_F2xz_S_bc = I_ERI_Gxy2z_S_F2xz_S_bc+ABX*I_ERI_Fy2z_S_F2xz_S_bc;
  Double I_ERI_F3z_Px_F2xz_S_bc = I_ERI_Gx3z_S_F2xz_S_bc+ABX*I_ERI_F3z_S_F2xz_S_bc;
  Double I_ERI_F3x_Py_F2xz_S_bc = I_ERI_G3xy_S_F2xz_S_bc+ABY*I_ERI_F3x_S_F2xz_S_bc;
  Double I_ERI_F2xy_Py_F2xz_S_bc = I_ERI_G2x2y_S_F2xz_S_bc+ABY*I_ERI_F2xy_S_F2xz_S_bc;
  Double I_ERI_F2xz_Py_F2xz_S_bc = I_ERI_G2xyz_S_F2xz_S_bc+ABY*I_ERI_F2xz_S_F2xz_S_bc;
  Double I_ERI_Fx2y_Py_F2xz_S_bc = I_ERI_Gx3y_S_F2xz_S_bc+ABY*I_ERI_Fx2y_S_F2xz_S_bc;
  Double I_ERI_Fxyz_Py_F2xz_S_bc = I_ERI_Gx2yz_S_F2xz_S_bc+ABY*I_ERI_Fxyz_S_F2xz_S_bc;
  Double I_ERI_Fx2z_Py_F2xz_S_bc = I_ERI_Gxy2z_S_F2xz_S_bc+ABY*I_ERI_Fx2z_S_F2xz_S_bc;
  Double I_ERI_F3y_Py_F2xz_S_bc = I_ERI_G4y_S_F2xz_S_bc+ABY*I_ERI_F3y_S_F2xz_S_bc;
  Double I_ERI_F2yz_Py_F2xz_S_bc = I_ERI_G3yz_S_F2xz_S_bc+ABY*I_ERI_F2yz_S_F2xz_S_bc;
  Double I_ERI_Fy2z_Py_F2xz_S_bc = I_ERI_G2y2z_S_F2xz_S_bc+ABY*I_ERI_Fy2z_S_F2xz_S_bc;
  Double I_ERI_F3z_Py_F2xz_S_bc = I_ERI_Gy3z_S_F2xz_S_bc+ABY*I_ERI_F3z_S_F2xz_S_bc;
  Double I_ERI_F3x_Pz_F2xz_S_bc = I_ERI_G3xz_S_F2xz_S_bc+ABZ*I_ERI_F3x_S_F2xz_S_bc;
  Double I_ERI_F2xy_Pz_F2xz_S_bc = I_ERI_G2xyz_S_F2xz_S_bc+ABZ*I_ERI_F2xy_S_F2xz_S_bc;
  Double I_ERI_F2xz_Pz_F2xz_S_bc = I_ERI_G2x2z_S_F2xz_S_bc+ABZ*I_ERI_F2xz_S_F2xz_S_bc;
  Double I_ERI_Fx2y_Pz_F2xz_S_bc = I_ERI_Gx2yz_S_F2xz_S_bc+ABZ*I_ERI_Fx2y_S_F2xz_S_bc;
  Double I_ERI_Fxyz_Pz_F2xz_S_bc = I_ERI_Gxy2z_S_F2xz_S_bc+ABZ*I_ERI_Fxyz_S_F2xz_S_bc;
  Double I_ERI_Fx2z_Pz_F2xz_S_bc = I_ERI_Gx3z_S_F2xz_S_bc+ABZ*I_ERI_Fx2z_S_F2xz_S_bc;
  Double I_ERI_F3y_Pz_F2xz_S_bc = I_ERI_G3yz_S_F2xz_S_bc+ABZ*I_ERI_F3y_S_F2xz_S_bc;
  Double I_ERI_F2yz_Pz_F2xz_S_bc = I_ERI_G2y2z_S_F2xz_S_bc+ABZ*I_ERI_F2yz_S_F2xz_S_bc;
  Double I_ERI_Fy2z_Pz_F2xz_S_bc = I_ERI_Gy3z_S_F2xz_S_bc+ABZ*I_ERI_Fy2z_S_F2xz_S_bc;
  Double I_ERI_F3z_Pz_F2xz_S_bc = I_ERI_G4z_S_F2xz_S_bc+ABZ*I_ERI_F3z_S_F2xz_S_bc;
  Double I_ERI_F3x_Px_Fx2y_S_bc = I_ERI_G4x_S_Fx2y_S_bc+ABX*I_ERI_F3x_S_Fx2y_S_bc;
  Double I_ERI_F2xy_Px_Fx2y_S_bc = I_ERI_G3xy_S_Fx2y_S_bc+ABX*I_ERI_F2xy_S_Fx2y_S_bc;
  Double I_ERI_F2xz_Px_Fx2y_S_bc = I_ERI_G3xz_S_Fx2y_S_bc+ABX*I_ERI_F2xz_S_Fx2y_S_bc;
  Double I_ERI_Fx2y_Px_Fx2y_S_bc = I_ERI_G2x2y_S_Fx2y_S_bc+ABX*I_ERI_Fx2y_S_Fx2y_S_bc;
  Double I_ERI_Fxyz_Px_Fx2y_S_bc = I_ERI_G2xyz_S_Fx2y_S_bc+ABX*I_ERI_Fxyz_S_Fx2y_S_bc;
  Double I_ERI_Fx2z_Px_Fx2y_S_bc = I_ERI_G2x2z_S_Fx2y_S_bc+ABX*I_ERI_Fx2z_S_Fx2y_S_bc;
  Double I_ERI_F3y_Px_Fx2y_S_bc = I_ERI_Gx3y_S_Fx2y_S_bc+ABX*I_ERI_F3y_S_Fx2y_S_bc;
  Double I_ERI_F2yz_Px_Fx2y_S_bc = I_ERI_Gx2yz_S_Fx2y_S_bc+ABX*I_ERI_F2yz_S_Fx2y_S_bc;
  Double I_ERI_Fy2z_Px_Fx2y_S_bc = I_ERI_Gxy2z_S_Fx2y_S_bc+ABX*I_ERI_Fy2z_S_Fx2y_S_bc;
  Double I_ERI_F3z_Px_Fx2y_S_bc = I_ERI_Gx3z_S_Fx2y_S_bc+ABX*I_ERI_F3z_S_Fx2y_S_bc;
  Double I_ERI_F3x_Py_Fx2y_S_bc = I_ERI_G3xy_S_Fx2y_S_bc+ABY*I_ERI_F3x_S_Fx2y_S_bc;
  Double I_ERI_F2xy_Py_Fx2y_S_bc = I_ERI_G2x2y_S_Fx2y_S_bc+ABY*I_ERI_F2xy_S_Fx2y_S_bc;
  Double I_ERI_F2xz_Py_Fx2y_S_bc = I_ERI_G2xyz_S_Fx2y_S_bc+ABY*I_ERI_F2xz_S_Fx2y_S_bc;
  Double I_ERI_Fx2y_Py_Fx2y_S_bc = I_ERI_Gx3y_S_Fx2y_S_bc+ABY*I_ERI_Fx2y_S_Fx2y_S_bc;
  Double I_ERI_Fxyz_Py_Fx2y_S_bc = I_ERI_Gx2yz_S_Fx2y_S_bc+ABY*I_ERI_Fxyz_S_Fx2y_S_bc;
  Double I_ERI_Fx2z_Py_Fx2y_S_bc = I_ERI_Gxy2z_S_Fx2y_S_bc+ABY*I_ERI_Fx2z_S_Fx2y_S_bc;
  Double I_ERI_F3y_Py_Fx2y_S_bc = I_ERI_G4y_S_Fx2y_S_bc+ABY*I_ERI_F3y_S_Fx2y_S_bc;
  Double I_ERI_F2yz_Py_Fx2y_S_bc = I_ERI_G3yz_S_Fx2y_S_bc+ABY*I_ERI_F2yz_S_Fx2y_S_bc;
  Double I_ERI_Fy2z_Py_Fx2y_S_bc = I_ERI_G2y2z_S_Fx2y_S_bc+ABY*I_ERI_Fy2z_S_Fx2y_S_bc;
  Double I_ERI_F3z_Py_Fx2y_S_bc = I_ERI_Gy3z_S_Fx2y_S_bc+ABY*I_ERI_F3z_S_Fx2y_S_bc;
  Double I_ERI_F3x_Pz_Fx2y_S_bc = I_ERI_G3xz_S_Fx2y_S_bc+ABZ*I_ERI_F3x_S_Fx2y_S_bc;
  Double I_ERI_F2xy_Pz_Fx2y_S_bc = I_ERI_G2xyz_S_Fx2y_S_bc+ABZ*I_ERI_F2xy_S_Fx2y_S_bc;
  Double I_ERI_F2xz_Pz_Fx2y_S_bc = I_ERI_G2x2z_S_Fx2y_S_bc+ABZ*I_ERI_F2xz_S_Fx2y_S_bc;
  Double I_ERI_Fx2y_Pz_Fx2y_S_bc = I_ERI_Gx2yz_S_Fx2y_S_bc+ABZ*I_ERI_Fx2y_S_Fx2y_S_bc;
  Double I_ERI_Fxyz_Pz_Fx2y_S_bc = I_ERI_Gxy2z_S_Fx2y_S_bc+ABZ*I_ERI_Fxyz_S_Fx2y_S_bc;
  Double I_ERI_Fx2z_Pz_Fx2y_S_bc = I_ERI_Gx3z_S_Fx2y_S_bc+ABZ*I_ERI_Fx2z_S_Fx2y_S_bc;
  Double I_ERI_F3y_Pz_Fx2y_S_bc = I_ERI_G3yz_S_Fx2y_S_bc+ABZ*I_ERI_F3y_S_Fx2y_S_bc;
  Double I_ERI_F2yz_Pz_Fx2y_S_bc = I_ERI_G2y2z_S_Fx2y_S_bc+ABZ*I_ERI_F2yz_S_Fx2y_S_bc;
  Double I_ERI_Fy2z_Pz_Fx2y_S_bc = I_ERI_Gy3z_S_Fx2y_S_bc+ABZ*I_ERI_Fy2z_S_Fx2y_S_bc;
  Double I_ERI_F3z_Pz_Fx2y_S_bc = I_ERI_G4z_S_Fx2y_S_bc+ABZ*I_ERI_F3z_S_Fx2y_S_bc;
  Double I_ERI_F3x_Px_Fxyz_S_bc = I_ERI_G4x_S_Fxyz_S_bc+ABX*I_ERI_F3x_S_Fxyz_S_bc;
  Double I_ERI_F2xy_Px_Fxyz_S_bc = I_ERI_G3xy_S_Fxyz_S_bc+ABX*I_ERI_F2xy_S_Fxyz_S_bc;
  Double I_ERI_F2xz_Px_Fxyz_S_bc = I_ERI_G3xz_S_Fxyz_S_bc+ABX*I_ERI_F2xz_S_Fxyz_S_bc;
  Double I_ERI_Fx2y_Px_Fxyz_S_bc = I_ERI_G2x2y_S_Fxyz_S_bc+ABX*I_ERI_Fx2y_S_Fxyz_S_bc;
  Double I_ERI_Fxyz_Px_Fxyz_S_bc = I_ERI_G2xyz_S_Fxyz_S_bc+ABX*I_ERI_Fxyz_S_Fxyz_S_bc;
  Double I_ERI_Fx2z_Px_Fxyz_S_bc = I_ERI_G2x2z_S_Fxyz_S_bc+ABX*I_ERI_Fx2z_S_Fxyz_S_bc;
  Double I_ERI_F3y_Px_Fxyz_S_bc = I_ERI_Gx3y_S_Fxyz_S_bc+ABX*I_ERI_F3y_S_Fxyz_S_bc;
  Double I_ERI_F2yz_Px_Fxyz_S_bc = I_ERI_Gx2yz_S_Fxyz_S_bc+ABX*I_ERI_F2yz_S_Fxyz_S_bc;
  Double I_ERI_Fy2z_Px_Fxyz_S_bc = I_ERI_Gxy2z_S_Fxyz_S_bc+ABX*I_ERI_Fy2z_S_Fxyz_S_bc;
  Double I_ERI_F3z_Px_Fxyz_S_bc = I_ERI_Gx3z_S_Fxyz_S_bc+ABX*I_ERI_F3z_S_Fxyz_S_bc;
  Double I_ERI_F3x_Py_Fxyz_S_bc = I_ERI_G3xy_S_Fxyz_S_bc+ABY*I_ERI_F3x_S_Fxyz_S_bc;
  Double I_ERI_F2xy_Py_Fxyz_S_bc = I_ERI_G2x2y_S_Fxyz_S_bc+ABY*I_ERI_F2xy_S_Fxyz_S_bc;
  Double I_ERI_F2xz_Py_Fxyz_S_bc = I_ERI_G2xyz_S_Fxyz_S_bc+ABY*I_ERI_F2xz_S_Fxyz_S_bc;
  Double I_ERI_Fx2y_Py_Fxyz_S_bc = I_ERI_Gx3y_S_Fxyz_S_bc+ABY*I_ERI_Fx2y_S_Fxyz_S_bc;
  Double I_ERI_Fxyz_Py_Fxyz_S_bc = I_ERI_Gx2yz_S_Fxyz_S_bc+ABY*I_ERI_Fxyz_S_Fxyz_S_bc;
  Double I_ERI_Fx2z_Py_Fxyz_S_bc = I_ERI_Gxy2z_S_Fxyz_S_bc+ABY*I_ERI_Fx2z_S_Fxyz_S_bc;
  Double I_ERI_F3y_Py_Fxyz_S_bc = I_ERI_G4y_S_Fxyz_S_bc+ABY*I_ERI_F3y_S_Fxyz_S_bc;
  Double I_ERI_F2yz_Py_Fxyz_S_bc = I_ERI_G3yz_S_Fxyz_S_bc+ABY*I_ERI_F2yz_S_Fxyz_S_bc;
  Double I_ERI_Fy2z_Py_Fxyz_S_bc = I_ERI_G2y2z_S_Fxyz_S_bc+ABY*I_ERI_Fy2z_S_Fxyz_S_bc;
  Double I_ERI_F3z_Py_Fxyz_S_bc = I_ERI_Gy3z_S_Fxyz_S_bc+ABY*I_ERI_F3z_S_Fxyz_S_bc;
  Double I_ERI_F3x_Pz_Fxyz_S_bc = I_ERI_G3xz_S_Fxyz_S_bc+ABZ*I_ERI_F3x_S_Fxyz_S_bc;
  Double I_ERI_F2xy_Pz_Fxyz_S_bc = I_ERI_G2xyz_S_Fxyz_S_bc+ABZ*I_ERI_F2xy_S_Fxyz_S_bc;
  Double I_ERI_F2xz_Pz_Fxyz_S_bc = I_ERI_G2x2z_S_Fxyz_S_bc+ABZ*I_ERI_F2xz_S_Fxyz_S_bc;
  Double I_ERI_Fx2y_Pz_Fxyz_S_bc = I_ERI_Gx2yz_S_Fxyz_S_bc+ABZ*I_ERI_Fx2y_S_Fxyz_S_bc;
  Double I_ERI_Fxyz_Pz_Fxyz_S_bc = I_ERI_Gxy2z_S_Fxyz_S_bc+ABZ*I_ERI_Fxyz_S_Fxyz_S_bc;
  Double I_ERI_Fx2z_Pz_Fxyz_S_bc = I_ERI_Gx3z_S_Fxyz_S_bc+ABZ*I_ERI_Fx2z_S_Fxyz_S_bc;
  Double I_ERI_F3y_Pz_Fxyz_S_bc = I_ERI_G3yz_S_Fxyz_S_bc+ABZ*I_ERI_F3y_S_Fxyz_S_bc;
  Double I_ERI_F2yz_Pz_Fxyz_S_bc = I_ERI_G2y2z_S_Fxyz_S_bc+ABZ*I_ERI_F2yz_S_Fxyz_S_bc;
  Double I_ERI_Fy2z_Pz_Fxyz_S_bc = I_ERI_Gy3z_S_Fxyz_S_bc+ABZ*I_ERI_Fy2z_S_Fxyz_S_bc;
  Double I_ERI_F3z_Pz_Fxyz_S_bc = I_ERI_G4z_S_Fxyz_S_bc+ABZ*I_ERI_F3z_S_Fxyz_S_bc;
  Double I_ERI_F3x_Px_Fx2z_S_bc = I_ERI_G4x_S_Fx2z_S_bc+ABX*I_ERI_F3x_S_Fx2z_S_bc;
  Double I_ERI_F2xy_Px_Fx2z_S_bc = I_ERI_G3xy_S_Fx2z_S_bc+ABX*I_ERI_F2xy_S_Fx2z_S_bc;
  Double I_ERI_F2xz_Px_Fx2z_S_bc = I_ERI_G3xz_S_Fx2z_S_bc+ABX*I_ERI_F2xz_S_Fx2z_S_bc;
  Double I_ERI_Fx2y_Px_Fx2z_S_bc = I_ERI_G2x2y_S_Fx2z_S_bc+ABX*I_ERI_Fx2y_S_Fx2z_S_bc;
  Double I_ERI_Fxyz_Px_Fx2z_S_bc = I_ERI_G2xyz_S_Fx2z_S_bc+ABX*I_ERI_Fxyz_S_Fx2z_S_bc;
  Double I_ERI_Fx2z_Px_Fx2z_S_bc = I_ERI_G2x2z_S_Fx2z_S_bc+ABX*I_ERI_Fx2z_S_Fx2z_S_bc;
  Double I_ERI_F3y_Px_Fx2z_S_bc = I_ERI_Gx3y_S_Fx2z_S_bc+ABX*I_ERI_F3y_S_Fx2z_S_bc;
  Double I_ERI_F2yz_Px_Fx2z_S_bc = I_ERI_Gx2yz_S_Fx2z_S_bc+ABX*I_ERI_F2yz_S_Fx2z_S_bc;
  Double I_ERI_Fy2z_Px_Fx2z_S_bc = I_ERI_Gxy2z_S_Fx2z_S_bc+ABX*I_ERI_Fy2z_S_Fx2z_S_bc;
  Double I_ERI_F3z_Px_Fx2z_S_bc = I_ERI_Gx3z_S_Fx2z_S_bc+ABX*I_ERI_F3z_S_Fx2z_S_bc;
  Double I_ERI_F3x_Py_Fx2z_S_bc = I_ERI_G3xy_S_Fx2z_S_bc+ABY*I_ERI_F3x_S_Fx2z_S_bc;
  Double I_ERI_F2xy_Py_Fx2z_S_bc = I_ERI_G2x2y_S_Fx2z_S_bc+ABY*I_ERI_F2xy_S_Fx2z_S_bc;
  Double I_ERI_F2xz_Py_Fx2z_S_bc = I_ERI_G2xyz_S_Fx2z_S_bc+ABY*I_ERI_F2xz_S_Fx2z_S_bc;
  Double I_ERI_Fx2y_Py_Fx2z_S_bc = I_ERI_Gx3y_S_Fx2z_S_bc+ABY*I_ERI_Fx2y_S_Fx2z_S_bc;
  Double I_ERI_Fxyz_Py_Fx2z_S_bc = I_ERI_Gx2yz_S_Fx2z_S_bc+ABY*I_ERI_Fxyz_S_Fx2z_S_bc;
  Double I_ERI_Fx2z_Py_Fx2z_S_bc = I_ERI_Gxy2z_S_Fx2z_S_bc+ABY*I_ERI_Fx2z_S_Fx2z_S_bc;
  Double I_ERI_F3y_Py_Fx2z_S_bc = I_ERI_G4y_S_Fx2z_S_bc+ABY*I_ERI_F3y_S_Fx2z_S_bc;
  Double I_ERI_F2yz_Py_Fx2z_S_bc = I_ERI_G3yz_S_Fx2z_S_bc+ABY*I_ERI_F2yz_S_Fx2z_S_bc;
  Double I_ERI_Fy2z_Py_Fx2z_S_bc = I_ERI_G2y2z_S_Fx2z_S_bc+ABY*I_ERI_Fy2z_S_Fx2z_S_bc;
  Double I_ERI_F3z_Py_Fx2z_S_bc = I_ERI_Gy3z_S_Fx2z_S_bc+ABY*I_ERI_F3z_S_Fx2z_S_bc;
  Double I_ERI_F3x_Pz_Fx2z_S_bc = I_ERI_G3xz_S_Fx2z_S_bc+ABZ*I_ERI_F3x_S_Fx2z_S_bc;
  Double I_ERI_F2xy_Pz_Fx2z_S_bc = I_ERI_G2xyz_S_Fx2z_S_bc+ABZ*I_ERI_F2xy_S_Fx2z_S_bc;
  Double I_ERI_F2xz_Pz_Fx2z_S_bc = I_ERI_G2x2z_S_Fx2z_S_bc+ABZ*I_ERI_F2xz_S_Fx2z_S_bc;
  Double I_ERI_Fx2y_Pz_Fx2z_S_bc = I_ERI_Gx2yz_S_Fx2z_S_bc+ABZ*I_ERI_Fx2y_S_Fx2z_S_bc;
  Double I_ERI_Fxyz_Pz_Fx2z_S_bc = I_ERI_Gxy2z_S_Fx2z_S_bc+ABZ*I_ERI_Fxyz_S_Fx2z_S_bc;
  Double I_ERI_Fx2z_Pz_Fx2z_S_bc = I_ERI_Gx3z_S_Fx2z_S_bc+ABZ*I_ERI_Fx2z_S_Fx2z_S_bc;
  Double I_ERI_F3y_Pz_Fx2z_S_bc = I_ERI_G3yz_S_Fx2z_S_bc+ABZ*I_ERI_F3y_S_Fx2z_S_bc;
  Double I_ERI_F2yz_Pz_Fx2z_S_bc = I_ERI_G2y2z_S_Fx2z_S_bc+ABZ*I_ERI_F2yz_S_Fx2z_S_bc;
  Double I_ERI_Fy2z_Pz_Fx2z_S_bc = I_ERI_Gy3z_S_Fx2z_S_bc+ABZ*I_ERI_Fy2z_S_Fx2z_S_bc;
  Double I_ERI_F3z_Pz_Fx2z_S_bc = I_ERI_G4z_S_Fx2z_S_bc+ABZ*I_ERI_F3z_S_Fx2z_S_bc;
  Double I_ERI_F3x_Px_F3y_S_bc = I_ERI_G4x_S_F3y_S_bc+ABX*I_ERI_F3x_S_F3y_S_bc;
  Double I_ERI_F2xy_Px_F3y_S_bc = I_ERI_G3xy_S_F3y_S_bc+ABX*I_ERI_F2xy_S_F3y_S_bc;
  Double I_ERI_F2xz_Px_F3y_S_bc = I_ERI_G3xz_S_F3y_S_bc+ABX*I_ERI_F2xz_S_F3y_S_bc;
  Double I_ERI_Fx2y_Px_F3y_S_bc = I_ERI_G2x2y_S_F3y_S_bc+ABX*I_ERI_Fx2y_S_F3y_S_bc;
  Double I_ERI_Fxyz_Px_F3y_S_bc = I_ERI_G2xyz_S_F3y_S_bc+ABX*I_ERI_Fxyz_S_F3y_S_bc;
  Double I_ERI_Fx2z_Px_F3y_S_bc = I_ERI_G2x2z_S_F3y_S_bc+ABX*I_ERI_Fx2z_S_F3y_S_bc;
  Double I_ERI_F3y_Px_F3y_S_bc = I_ERI_Gx3y_S_F3y_S_bc+ABX*I_ERI_F3y_S_F3y_S_bc;
  Double I_ERI_F2yz_Px_F3y_S_bc = I_ERI_Gx2yz_S_F3y_S_bc+ABX*I_ERI_F2yz_S_F3y_S_bc;
  Double I_ERI_Fy2z_Px_F3y_S_bc = I_ERI_Gxy2z_S_F3y_S_bc+ABX*I_ERI_Fy2z_S_F3y_S_bc;
  Double I_ERI_F3z_Px_F3y_S_bc = I_ERI_Gx3z_S_F3y_S_bc+ABX*I_ERI_F3z_S_F3y_S_bc;
  Double I_ERI_F3x_Py_F3y_S_bc = I_ERI_G3xy_S_F3y_S_bc+ABY*I_ERI_F3x_S_F3y_S_bc;
  Double I_ERI_F2xy_Py_F3y_S_bc = I_ERI_G2x2y_S_F3y_S_bc+ABY*I_ERI_F2xy_S_F3y_S_bc;
  Double I_ERI_F2xz_Py_F3y_S_bc = I_ERI_G2xyz_S_F3y_S_bc+ABY*I_ERI_F2xz_S_F3y_S_bc;
  Double I_ERI_Fx2y_Py_F3y_S_bc = I_ERI_Gx3y_S_F3y_S_bc+ABY*I_ERI_Fx2y_S_F3y_S_bc;
  Double I_ERI_Fxyz_Py_F3y_S_bc = I_ERI_Gx2yz_S_F3y_S_bc+ABY*I_ERI_Fxyz_S_F3y_S_bc;
  Double I_ERI_Fx2z_Py_F3y_S_bc = I_ERI_Gxy2z_S_F3y_S_bc+ABY*I_ERI_Fx2z_S_F3y_S_bc;
  Double I_ERI_F3y_Py_F3y_S_bc = I_ERI_G4y_S_F3y_S_bc+ABY*I_ERI_F3y_S_F3y_S_bc;
  Double I_ERI_F2yz_Py_F3y_S_bc = I_ERI_G3yz_S_F3y_S_bc+ABY*I_ERI_F2yz_S_F3y_S_bc;
  Double I_ERI_Fy2z_Py_F3y_S_bc = I_ERI_G2y2z_S_F3y_S_bc+ABY*I_ERI_Fy2z_S_F3y_S_bc;
  Double I_ERI_F3z_Py_F3y_S_bc = I_ERI_Gy3z_S_F3y_S_bc+ABY*I_ERI_F3z_S_F3y_S_bc;
  Double I_ERI_F3x_Pz_F3y_S_bc = I_ERI_G3xz_S_F3y_S_bc+ABZ*I_ERI_F3x_S_F3y_S_bc;
  Double I_ERI_F2xy_Pz_F3y_S_bc = I_ERI_G2xyz_S_F3y_S_bc+ABZ*I_ERI_F2xy_S_F3y_S_bc;
  Double I_ERI_F2xz_Pz_F3y_S_bc = I_ERI_G2x2z_S_F3y_S_bc+ABZ*I_ERI_F2xz_S_F3y_S_bc;
  Double I_ERI_Fx2y_Pz_F3y_S_bc = I_ERI_Gx2yz_S_F3y_S_bc+ABZ*I_ERI_Fx2y_S_F3y_S_bc;
  Double I_ERI_Fxyz_Pz_F3y_S_bc = I_ERI_Gxy2z_S_F3y_S_bc+ABZ*I_ERI_Fxyz_S_F3y_S_bc;
  Double I_ERI_Fx2z_Pz_F3y_S_bc = I_ERI_Gx3z_S_F3y_S_bc+ABZ*I_ERI_Fx2z_S_F3y_S_bc;
  Double I_ERI_F3y_Pz_F3y_S_bc = I_ERI_G3yz_S_F3y_S_bc+ABZ*I_ERI_F3y_S_F3y_S_bc;
  Double I_ERI_F2yz_Pz_F3y_S_bc = I_ERI_G2y2z_S_F3y_S_bc+ABZ*I_ERI_F2yz_S_F3y_S_bc;
  Double I_ERI_Fy2z_Pz_F3y_S_bc = I_ERI_Gy3z_S_F3y_S_bc+ABZ*I_ERI_Fy2z_S_F3y_S_bc;
  Double I_ERI_F3z_Pz_F3y_S_bc = I_ERI_G4z_S_F3y_S_bc+ABZ*I_ERI_F3z_S_F3y_S_bc;
  Double I_ERI_F3x_Px_F2yz_S_bc = I_ERI_G4x_S_F2yz_S_bc+ABX*I_ERI_F3x_S_F2yz_S_bc;
  Double I_ERI_F2xy_Px_F2yz_S_bc = I_ERI_G3xy_S_F2yz_S_bc+ABX*I_ERI_F2xy_S_F2yz_S_bc;
  Double I_ERI_F2xz_Px_F2yz_S_bc = I_ERI_G3xz_S_F2yz_S_bc+ABX*I_ERI_F2xz_S_F2yz_S_bc;
  Double I_ERI_Fx2y_Px_F2yz_S_bc = I_ERI_G2x2y_S_F2yz_S_bc+ABX*I_ERI_Fx2y_S_F2yz_S_bc;
  Double I_ERI_Fxyz_Px_F2yz_S_bc = I_ERI_G2xyz_S_F2yz_S_bc+ABX*I_ERI_Fxyz_S_F2yz_S_bc;
  Double I_ERI_Fx2z_Px_F2yz_S_bc = I_ERI_G2x2z_S_F2yz_S_bc+ABX*I_ERI_Fx2z_S_F2yz_S_bc;
  Double I_ERI_F3y_Px_F2yz_S_bc = I_ERI_Gx3y_S_F2yz_S_bc+ABX*I_ERI_F3y_S_F2yz_S_bc;
  Double I_ERI_F2yz_Px_F2yz_S_bc = I_ERI_Gx2yz_S_F2yz_S_bc+ABX*I_ERI_F2yz_S_F2yz_S_bc;
  Double I_ERI_Fy2z_Px_F2yz_S_bc = I_ERI_Gxy2z_S_F2yz_S_bc+ABX*I_ERI_Fy2z_S_F2yz_S_bc;
  Double I_ERI_F3z_Px_F2yz_S_bc = I_ERI_Gx3z_S_F2yz_S_bc+ABX*I_ERI_F3z_S_F2yz_S_bc;
  Double I_ERI_F3x_Py_F2yz_S_bc = I_ERI_G3xy_S_F2yz_S_bc+ABY*I_ERI_F3x_S_F2yz_S_bc;
  Double I_ERI_F2xy_Py_F2yz_S_bc = I_ERI_G2x2y_S_F2yz_S_bc+ABY*I_ERI_F2xy_S_F2yz_S_bc;
  Double I_ERI_F2xz_Py_F2yz_S_bc = I_ERI_G2xyz_S_F2yz_S_bc+ABY*I_ERI_F2xz_S_F2yz_S_bc;
  Double I_ERI_Fx2y_Py_F2yz_S_bc = I_ERI_Gx3y_S_F2yz_S_bc+ABY*I_ERI_Fx2y_S_F2yz_S_bc;
  Double I_ERI_Fxyz_Py_F2yz_S_bc = I_ERI_Gx2yz_S_F2yz_S_bc+ABY*I_ERI_Fxyz_S_F2yz_S_bc;
  Double I_ERI_Fx2z_Py_F2yz_S_bc = I_ERI_Gxy2z_S_F2yz_S_bc+ABY*I_ERI_Fx2z_S_F2yz_S_bc;
  Double I_ERI_F3y_Py_F2yz_S_bc = I_ERI_G4y_S_F2yz_S_bc+ABY*I_ERI_F3y_S_F2yz_S_bc;
  Double I_ERI_F2yz_Py_F2yz_S_bc = I_ERI_G3yz_S_F2yz_S_bc+ABY*I_ERI_F2yz_S_F2yz_S_bc;
  Double I_ERI_Fy2z_Py_F2yz_S_bc = I_ERI_G2y2z_S_F2yz_S_bc+ABY*I_ERI_Fy2z_S_F2yz_S_bc;
  Double I_ERI_F3z_Py_F2yz_S_bc = I_ERI_Gy3z_S_F2yz_S_bc+ABY*I_ERI_F3z_S_F2yz_S_bc;
  Double I_ERI_F3x_Pz_F2yz_S_bc = I_ERI_G3xz_S_F2yz_S_bc+ABZ*I_ERI_F3x_S_F2yz_S_bc;
  Double I_ERI_F2xy_Pz_F2yz_S_bc = I_ERI_G2xyz_S_F2yz_S_bc+ABZ*I_ERI_F2xy_S_F2yz_S_bc;
  Double I_ERI_F2xz_Pz_F2yz_S_bc = I_ERI_G2x2z_S_F2yz_S_bc+ABZ*I_ERI_F2xz_S_F2yz_S_bc;
  Double I_ERI_Fx2y_Pz_F2yz_S_bc = I_ERI_Gx2yz_S_F2yz_S_bc+ABZ*I_ERI_Fx2y_S_F2yz_S_bc;
  Double I_ERI_Fxyz_Pz_F2yz_S_bc = I_ERI_Gxy2z_S_F2yz_S_bc+ABZ*I_ERI_Fxyz_S_F2yz_S_bc;
  Double I_ERI_Fx2z_Pz_F2yz_S_bc = I_ERI_Gx3z_S_F2yz_S_bc+ABZ*I_ERI_Fx2z_S_F2yz_S_bc;
  Double I_ERI_F3y_Pz_F2yz_S_bc = I_ERI_G3yz_S_F2yz_S_bc+ABZ*I_ERI_F3y_S_F2yz_S_bc;
  Double I_ERI_F2yz_Pz_F2yz_S_bc = I_ERI_G2y2z_S_F2yz_S_bc+ABZ*I_ERI_F2yz_S_F2yz_S_bc;
  Double I_ERI_Fy2z_Pz_F2yz_S_bc = I_ERI_Gy3z_S_F2yz_S_bc+ABZ*I_ERI_Fy2z_S_F2yz_S_bc;
  Double I_ERI_F3z_Pz_F2yz_S_bc = I_ERI_G4z_S_F2yz_S_bc+ABZ*I_ERI_F3z_S_F2yz_S_bc;
  Double I_ERI_F3x_Px_Fy2z_S_bc = I_ERI_G4x_S_Fy2z_S_bc+ABX*I_ERI_F3x_S_Fy2z_S_bc;
  Double I_ERI_F2xy_Px_Fy2z_S_bc = I_ERI_G3xy_S_Fy2z_S_bc+ABX*I_ERI_F2xy_S_Fy2z_S_bc;
  Double I_ERI_F2xz_Px_Fy2z_S_bc = I_ERI_G3xz_S_Fy2z_S_bc+ABX*I_ERI_F2xz_S_Fy2z_S_bc;
  Double I_ERI_Fx2y_Px_Fy2z_S_bc = I_ERI_G2x2y_S_Fy2z_S_bc+ABX*I_ERI_Fx2y_S_Fy2z_S_bc;
  Double I_ERI_Fxyz_Px_Fy2z_S_bc = I_ERI_G2xyz_S_Fy2z_S_bc+ABX*I_ERI_Fxyz_S_Fy2z_S_bc;
  Double I_ERI_Fx2z_Px_Fy2z_S_bc = I_ERI_G2x2z_S_Fy2z_S_bc+ABX*I_ERI_Fx2z_S_Fy2z_S_bc;
  Double I_ERI_F3y_Px_Fy2z_S_bc = I_ERI_Gx3y_S_Fy2z_S_bc+ABX*I_ERI_F3y_S_Fy2z_S_bc;
  Double I_ERI_F2yz_Px_Fy2z_S_bc = I_ERI_Gx2yz_S_Fy2z_S_bc+ABX*I_ERI_F2yz_S_Fy2z_S_bc;
  Double I_ERI_Fy2z_Px_Fy2z_S_bc = I_ERI_Gxy2z_S_Fy2z_S_bc+ABX*I_ERI_Fy2z_S_Fy2z_S_bc;
  Double I_ERI_F3z_Px_Fy2z_S_bc = I_ERI_Gx3z_S_Fy2z_S_bc+ABX*I_ERI_F3z_S_Fy2z_S_bc;
  Double I_ERI_F3x_Py_Fy2z_S_bc = I_ERI_G3xy_S_Fy2z_S_bc+ABY*I_ERI_F3x_S_Fy2z_S_bc;
  Double I_ERI_F2xy_Py_Fy2z_S_bc = I_ERI_G2x2y_S_Fy2z_S_bc+ABY*I_ERI_F2xy_S_Fy2z_S_bc;
  Double I_ERI_F2xz_Py_Fy2z_S_bc = I_ERI_G2xyz_S_Fy2z_S_bc+ABY*I_ERI_F2xz_S_Fy2z_S_bc;
  Double I_ERI_Fx2y_Py_Fy2z_S_bc = I_ERI_Gx3y_S_Fy2z_S_bc+ABY*I_ERI_Fx2y_S_Fy2z_S_bc;
  Double I_ERI_Fxyz_Py_Fy2z_S_bc = I_ERI_Gx2yz_S_Fy2z_S_bc+ABY*I_ERI_Fxyz_S_Fy2z_S_bc;
  Double I_ERI_Fx2z_Py_Fy2z_S_bc = I_ERI_Gxy2z_S_Fy2z_S_bc+ABY*I_ERI_Fx2z_S_Fy2z_S_bc;
  Double I_ERI_F3y_Py_Fy2z_S_bc = I_ERI_G4y_S_Fy2z_S_bc+ABY*I_ERI_F3y_S_Fy2z_S_bc;
  Double I_ERI_F2yz_Py_Fy2z_S_bc = I_ERI_G3yz_S_Fy2z_S_bc+ABY*I_ERI_F2yz_S_Fy2z_S_bc;
  Double I_ERI_Fy2z_Py_Fy2z_S_bc = I_ERI_G2y2z_S_Fy2z_S_bc+ABY*I_ERI_Fy2z_S_Fy2z_S_bc;
  Double I_ERI_F3z_Py_Fy2z_S_bc = I_ERI_Gy3z_S_Fy2z_S_bc+ABY*I_ERI_F3z_S_Fy2z_S_bc;
  Double I_ERI_F3x_Pz_Fy2z_S_bc = I_ERI_G3xz_S_Fy2z_S_bc+ABZ*I_ERI_F3x_S_Fy2z_S_bc;
  Double I_ERI_F2xy_Pz_Fy2z_S_bc = I_ERI_G2xyz_S_Fy2z_S_bc+ABZ*I_ERI_F2xy_S_Fy2z_S_bc;
  Double I_ERI_F2xz_Pz_Fy2z_S_bc = I_ERI_G2x2z_S_Fy2z_S_bc+ABZ*I_ERI_F2xz_S_Fy2z_S_bc;
  Double I_ERI_Fx2y_Pz_Fy2z_S_bc = I_ERI_Gx2yz_S_Fy2z_S_bc+ABZ*I_ERI_Fx2y_S_Fy2z_S_bc;
  Double I_ERI_Fxyz_Pz_Fy2z_S_bc = I_ERI_Gxy2z_S_Fy2z_S_bc+ABZ*I_ERI_Fxyz_S_Fy2z_S_bc;
  Double I_ERI_Fx2z_Pz_Fy2z_S_bc = I_ERI_Gx3z_S_Fy2z_S_bc+ABZ*I_ERI_Fx2z_S_Fy2z_S_bc;
  Double I_ERI_F3y_Pz_Fy2z_S_bc = I_ERI_G3yz_S_Fy2z_S_bc+ABZ*I_ERI_F3y_S_Fy2z_S_bc;
  Double I_ERI_F2yz_Pz_Fy2z_S_bc = I_ERI_G2y2z_S_Fy2z_S_bc+ABZ*I_ERI_F2yz_S_Fy2z_S_bc;
  Double I_ERI_Fy2z_Pz_Fy2z_S_bc = I_ERI_Gy3z_S_Fy2z_S_bc+ABZ*I_ERI_Fy2z_S_Fy2z_S_bc;
  Double I_ERI_F3z_Pz_Fy2z_S_bc = I_ERI_G4z_S_Fy2z_S_bc+ABZ*I_ERI_F3z_S_Fy2z_S_bc;
  Double I_ERI_F3x_Px_F3z_S_bc = I_ERI_G4x_S_F3z_S_bc+ABX*I_ERI_F3x_S_F3z_S_bc;
  Double I_ERI_F2xy_Px_F3z_S_bc = I_ERI_G3xy_S_F3z_S_bc+ABX*I_ERI_F2xy_S_F3z_S_bc;
  Double I_ERI_F2xz_Px_F3z_S_bc = I_ERI_G3xz_S_F3z_S_bc+ABX*I_ERI_F2xz_S_F3z_S_bc;
  Double I_ERI_Fx2y_Px_F3z_S_bc = I_ERI_G2x2y_S_F3z_S_bc+ABX*I_ERI_Fx2y_S_F3z_S_bc;
  Double I_ERI_Fxyz_Px_F3z_S_bc = I_ERI_G2xyz_S_F3z_S_bc+ABX*I_ERI_Fxyz_S_F3z_S_bc;
  Double I_ERI_Fx2z_Px_F3z_S_bc = I_ERI_G2x2z_S_F3z_S_bc+ABX*I_ERI_Fx2z_S_F3z_S_bc;
  Double I_ERI_F3y_Px_F3z_S_bc = I_ERI_Gx3y_S_F3z_S_bc+ABX*I_ERI_F3y_S_F3z_S_bc;
  Double I_ERI_F2yz_Px_F3z_S_bc = I_ERI_Gx2yz_S_F3z_S_bc+ABX*I_ERI_F2yz_S_F3z_S_bc;
  Double I_ERI_Fy2z_Px_F3z_S_bc = I_ERI_Gxy2z_S_F3z_S_bc+ABX*I_ERI_Fy2z_S_F3z_S_bc;
  Double I_ERI_F3z_Px_F3z_S_bc = I_ERI_Gx3z_S_F3z_S_bc+ABX*I_ERI_F3z_S_F3z_S_bc;
  Double I_ERI_F3x_Py_F3z_S_bc = I_ERI_G3xy_S_F3z_S_bc+ABY*I_ERI_F3x_S_F3z_S_bc;
  Double I_ERI_F2xy_Py_F3z_S_bc = I_ERI_G2x2y_S_F3z_S_bc+ABY*I_ERI_F2xy_S_F3z_S_bc;
  Double I_ERI_F2xz_Py_F3z_S_bc = I_ERI_G2xyz_S_F3z_S_bc+ABY*I_ERI_F2xz_S_F3z_S_bc;
  Double I_ERI_Fx2y_Py_F3z_S_bc = I_ERI_Gx3y_S_F3z_S_bc+ABY*I_ERI_Fx2y_S_F3z_S_bc;
  Double I_ERI_Fxyz_Py_F3z_S_bc = I_ERI_Gx2yz_S_F3z_S_bc+ABY*I_ERI_Fxyz_S_F3z_S_bc;
  Double I_ERI_Fx2z_Py_F3z_S_bc = I_ERI_Gxy2z_S_F3z_S_bc+ABY*I_ERI_Fx2z_S_F3z_S_bc;
  Double I_ERI_F3y_Py_F3z_S_bc = I_ERI_G4y_S_F3z_S_bc+ABY*I_ERI_F3y_S_F3z_S_bc;
  Double I_ERI_F2yz_Py_F3z_S_bc = I_ERI_G3yz_S_F3z_S_bc+ABY*I_ERI_F2yz_S_F3z_S_bc;
  Double I_ERI_Fy2z_Py_F3z_S_bc = I_ERI_G2y2z_S_F3z_S_bc+ABY*I_ERI_Fy2z_S_F3z_S_bc;
  Double I_ERI_F3z_Py_F3z_S_bc = I_ERI_Gy3z_S_F3z_S_bc+ABY*I_ERI_F3z_S_F3z_S_bc;
  Double I_ERI_F3x_Pz_F3z_S_bc = I_ERI_G3xz_S_F3z_S_bc+ABZ*I_ERI_F3x_S_F3z_S_bc;
  Double I_ERI_F2xy_Pz_F3z_S_bc = I_ERI_G2xyz_S_F3z_S_bc+ABZ*I_ERI_F2xy_S_F3z_S_bc;
  Double I_ERI_F2xz_Pz_F3z_S_bc = I_ERI_G2x2z_S_F3z_S_bc+ABZ*I_ERI_F2xz_S_F3z_S_bc;
  Double I_ERI_Fx2y_Pz_F3z_S_bc = I_ERI_Gx2yz_S_F3z_S_bc+ABZ*I_ERI_Fx2y_S_F3z_S_bc;
  Double I_ERI_Fxyz_Pz_F3z_S_bc = I_ERI_Gxy2z_S_F3z_S_bc+ABZ*I_ERI_Fxyz_S_F3z_S_bc;
  Double I_ERI_Fx2z_Pz_F3z_S_bc = I_ERI_Gx3z_S_F3z_S_bc+ABZ*I_ERI_Fx2z_S_F3z_S_bc;
  Double I_ERI_F3y_Pz_F3z_S_bc = I_ERI_G3yz_S_F3z_S_bc+ABZ*I_ERI_F3y_S_F3z_S_bc;
  Double I_ERI_F2yz_Pz_F3z_S_bc = I_ERI_G2y2z_S_F3z_S_bc+ABZ*I_ERI_F2yz_S_F3z_S_bc;
  Double I_ERI_Fy2z_Pz_F3z_S_bc = I_ERI_Gy3z_S_F3z_S_bc+ABZ*I_ERI_Fy2z_S_F3z_S_bc;
  Double I_ERI_F3z_Pz_F3z_S_bc = I_ERI_G4z_S_F3z_S_bc+ABZ*I_ERI_F3z_S_F3z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_P_bd
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_P_bd
   * RHS shell quartet name: SQ_ERI_F_S_D_P_bd
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_Px_bd = I_ERI_G4x_S_D2x_Px_bd+ABX*I_ERI_F3x_S_D2x_Px_bd;
  Double I_ERI_F2xy_Px_D2x_Px_bd = I_ERI_G3xy_S_D2x_Px_bd+ABX*I_ERI_F2xy_S_D2x_Px_bd;
  Double I_ERI_F2xz_Px_D2x_Px_bd = I_ERI_G3xz_S_D2x_Px_bd+ABX*I_ERI_F2xz_S_D2x_Px_bd;
  Double I_ERI_Fx2y_Px_D2x_Px_bd = I_ERI_G2x2y_S_D2x_Px_bd+ABX*I_ERI_Fx2y_S_D2x_Px_bd;
  Double I_ERI_Fxyz_Px_D2x_Px_bd = I_ERI_G2xyz_S_D2x_Px_bd+ABX*I_ERI_Fxyz_S_D2x_Px_bd;
  Double I_ERI_Fx2z_Px_D2x_Px_bd = I_ERI_G2x2z_S_D2x_Px_bd+ABX*I_ERI_Fx2z_S_D2x_Px_bd;
  Double I_ERI_F3y_Px_D2x_Px_bd = I_ERI_Gx3y_S_D2x_Px_bd+ABX*I_ERI_F3y_S_D2x_Px_bd;
  Double I_ERI_F2yz_Px_D2x_Px_bd = I_ERI_Gx2yz_S_D2x_Px_bd+ABX*I_ERI_F2yz_S_D2x_Px_bd;
  Double I_ERI_Fy2z_Px_D2x_Px_bd = I_ERI_Gxy2z_S_D2x_Px_bd+ABX*I_ERI_Fy2z_S_D2x_Px_bd;
  Double I_ERI_F3z_Px_D2x_Px_bd = I_ERI_Gx3z_S_D2x_Px_bd+ABX*I_ERI_F3z_S_D2x_Px_bd;
  Double I_ERI_F3x_Py_D2x_Px_bd = I_ERI_G3xy_S_D2x_Px_bd+ABY*I_ERI_F3x_S_D2x_Px_bd;
  Double I_ERI_F2xy_Py_D2x_Px_bd = I_ERI_G2x2y_S_D2x_Px_bd+ABY*I_ERI_F2xy_S_D2x_Px_bd;
  Double I_ERI_F2xz_Py_D2x_Px_bd = I_ERI_G2xyz_S_D2x_Px_bd+ABY*I_ERI_F2xz_S_D2x_Px_bd;
  Double I_ERI_Fx2y_Py_D2x_Px_bd = I_ERI_Gx3y_S_D2x_Px_bd+ABY*I_ERI_Fx2y_S_D2x_Px_bd;
  Double I_ERI_Fxyz_Py_D2x_Px_bd = I_ERI_Gx2yz_S_D2x_Px_bd+ABY*I_ERI_Fxyz_S_D2x_Px_bd;
  Double I_ERI_Fx2z_Py_D2x_Px_bd = I_ERI_Gxy2z_S_D2x_Px_bd+ABY*I_ERI_Fx2z_S_D2x_Px_bd;
  Double I_ERI_F3y_Py_D2x_Px_bd = I_ERI_G4y_S_D2x_Px_bd+ABY*I_ERI_F3y_S_D2x_Px_bd;
  Double I_ERI_F2yz_Py_D2x_Px_bd = I_ERI_G3yz_S_D2x_Px_bd+ABY*I_ERI_F2yz_S_D2x_Px_bd;
  Double I_ERI_Fy2z_Py_D2x_Px_bd = I_ERI_G2y2z_S_D2x_Px_bd+ABY*I_ERI_Fy2z_S_D2x_Px_bd;
  Double I_ERI_F3z_Py_D2x_Px_bd = I_ERI_Gy3z_S_D2x_Px_bd+ABY*I_ERI_F3z_S_D2x_Px_bd;
  Double I_ERI_F3x_Pz_D2x_Px_bd = I_ERI_G3xz_S_D2x_Px_bd+ABZ*I_ERI_F3x_S_D2x_Px_bd;
  Double I_ERI_F2xy_Pz_D2x_Px_bd = I_ERI_G2xyz_S_D2x_Px_bd+ABZ*I_ERI_F2xy_S_D2x_Px_bd;
  Double I_ERI_F2xz_Pz_D2x_Px_bd = I_ERI_G2x2z_S_D2x_Px_bd+ABZ*I_ERI_F2xz_S_D2x_Px_bd;
  Double I_ERI_Fx2y_Pz_D2x_Px_bd = I_ERI_Gx2yz_S_D2x_Px_bd+ABZ*I_ERI_Fx2y_S_D2x_Px_bd;
  Double I_ERI_Fxyz_Pz_D2x_Px_bd = I_ERI_Gxy2z_S_D2x_Px_bd+ABZ*I_ERI_Fxyz_S_D2x_Px_bd;
  Double I_ERI_Fx2z_Pz_D2x_Px_bd = I_ERI_Gx3z_S_D2x_Px_bd+ABZ*I_ERI_Fx2z_S_D2x_Px_bd;
  Double I_ERI_F3y_Pz_D2x_Px_bd = I_ERI_G3yz_S_D2x_Px_bd+ABZ*I_ERI_F3y_S_D2x_Px_bd;
  Double I_ERI_F2yz_Pz_D2x_Px_bd = I_ERI_G2y2z_S_D2x_Px_bd+ABZ*I_ERI_F2yz_S_D2x_Px_bd;
  Double I_ERI_Fy2z_Pz_D2x_Px_bd = I_ERI_Gy3z_S_D2x_Px_bd+ABZ*I_ERI_Fy2z_S_D2x_Px_bd;
  Double I_ERI_F3z_Pz_D2x_Px_bd = I_ERI_G4z_S_D2x_Px_bd+ABZ*I_ERI_F3z_S_D2x_Px_bd;
  Double I_ERI_F3x_Px_Dxy_Px_bd = I_ERI_G4x_S_Dxy_Px_bd+ABX*I_ERI_F3x_S_Dxy_Px_bd;
  Double I_ERI_F2xy_Px_Dxy_Px_bd = I_ERI_G3xy_S_Dxy_Px_bd+ABX*I_ERI_F2xy_S_Dxy_Px_bd;
  Double I_ERI_F2xz_Px_Dxy_Px_bd = I_ERI_G3xz_S_Dxy_Px_bd+ABX*I_ERI_F2xz_S_Dxy_Px_bd;
  Double I_ERI_Fx2y_Px_Dxy_Px_bd = I_ERI_G2x2y_S_Dxy_Px_bd+ABX*I_ERI_Fx2y_S_Dxy_Px_bd;
  Double I_ERI_Fxyz_Px_Dxy_Px_bd = I_ERI_G2xyz_S_Dxy_Px_bd+ABX*I_ERI_Fxyz_S_Dxy_Px_bd;
  Double I_ERI_Fx2z_Px_Dxy_Px_bd = I_ERI_G2x2z_S_Dxy_Px_bd+ABX*I_ERI_Fx2z_S_Dxy_Px_bd;
  Double I_ERI_F3y_Px_Dxy_Px_bd = I_ERI_Gx3y_S_Dxy_Px_bd+ABX*I_ERI_F3y_S_Dxy_Px_bd;
  Double I_ERI_F2yz_Px_Dxy_Px_bd = I_ERI_Gx2yz_S_Dxy_Px_bd+ABX*I_ERI_F2yz_S_Dxy_Px_bd;
  Double I_ERI_Fy2z_Px_Dxy_Px_bd = I_ERI_Gxy2z_S_Dxy_Px_bd+ABX*I_ERI_Fy2z_S_Dxy_Px_bd;
  Double I_ERI_F3z_Px_Dxy_Px_bd = I_ERI_Gx3z_S_Dxy_Px_bd+ABX*I_ERI_F3z_S_Dxy_Px_bd;
  Double I_ERI_F3x_Py_Dxy_Px_bd = I_ERI_G3xy_S_Dxy_Px_bd+ABY*I_ERI_F3x_S_Dxy_Px_bd;
  Double I_ERI_F2xy_Py_Dxy_Px_bd = I_ERI_G2x2y_S_Dxy_Px_bd+ABY*I_ERI_F2xy_S_Dxy_Px_bd;
  Double I_ERI_F2xz_Py_Dxy_Px_bd = I_ERI_G2xyz_S_Dxy_Px_bd+ABY*I_ERI_F2xz_S_Dxy_Px_bd;
  Double I_ERI_Fx2y_Py_Dxy_Px_bd = I_ERI_Gx3y_S_Dxy_Px_bd+ABY*I_ERI_Fx2y_S_Dxy_Px_bd;
  Double I_ERI_Fxyz_Py_Dxy_Px_bd = I_ERI_Gx2yz_S_Dxy_Px_bd+ABY*I_ERI_Fxyz_S_Dxy_Px_bd;
  Double I_ERI_Fx2z_Py_Dxy_Px_bd = I_ERI_Gxy2z_S_Dxy_Px_bd+ABY*I_ERI_Fx2z_S_Dxy_Px_bd;
  Double I_ERI_F3y_Py_Dxy_Px_bd = I_ERI_G4y_S_Dxy_Px_bd+ABY*I_ERI_F3y_S_Dxy_Px_bd;
  Double I_ERI_F2yz_Py_Dxy_Px_bd = I_ERI_G3yz_S_Dxy_Px_bd+ABY*I_ERI_F2yz_S_Dxy_Px_bd;
  Double I_ERI_Fy2z_Py_Dxy_Px_bd = I_ERI_G2y2z_S_Dxy_Px_bd+ABY*I_ERI_Fy2z_S_Dxy_Px_bd;
  Double I_ERI_F3z_Py_Dxy_Px_bd = I_ERI_Gy3z_S_Dxy_Px_bd+ABY*I_ERI_F3z_S_Dxy_Px_bd;
  Double I_ERI_F3x_Pz_Dxy_Px_bd = I_ERI_G3xz_S_Dxy_Px_bd+ABZ*I_ERI_F3x_S_Dxy_Px_bd;
  Double I_ERI_F2xy_Pz_Dxy_Px_bd = I_ERI_G2xyz_S_Dxy_Px_bd+ABZ*I_ERI_F2xy_S_Dxy_Px_bd;
  Double I_ERI_F2xz_Pz_Dxy_Px_bd = I_ERI_G2x2z_S_Dxy_Px_bd+ABZ*I_ERI_F2xz_S_Dxy_Px_bd;
  Double I_ERI_Fx2y_Pz_Dxy_Px_bd = I_ERI_Gx2yz_S_Dxy_Px_bd+ABZ*I_ERI_Fx2y_S_Dxy_Px_bd;
  Double I_ERI_Fxyz_Pz_Dxy_Px_bd = I_ERI_Gxy2z_S_Dxy_Px_bd+ABZ*I_ERI_Fxyz_S_Dxy_Px_bd;
  Double I_ERI_Fx2z_Pz_Dxy_Px_bd = I_ERI_Gx3z_S_Dxy_Px_bd+ABZ*I_ERI_Fx2z_S_Dxy_Px_bd;
  Double I_ERI_F3y_Pz_Dxy_Px_bd = I_ERI_G3yz_S_Dxy_Px_bd+ABZ*I_ERI_F3y_S_Dxy_Px_bd;
  Double I_ERI_F2yz_Pz_Dxy_Px_bd = I_ERI_G2y2z_S_Dxy_Px_bd+ABZ*I_ERI_F2yz_S_Dxy_Px_bd;
  Double I_ERI_Fy2z_Pz_Dxy_Px_bd = I_ERI_Gy3z_S_Dxy_Px_bd+ABZ*I_ERI_Fy2z_S_Dxy_Px_bd;
  Double I_ERI_F3z_Pz_Dxy_Px_bd = I_ERI_G4z_S_Dxy_Px_bd+ABZ*I_ERI_F3z_S_Dxy_Px_bd;
  Double I_ERI_F3x_Px_Dxz_Px_bd = I_ERI_G4x_S_Dxz_Px_bd+ABX*I_ERI_F3x_S_Dxz_Px_bd;
  Double I_ERI_F2xy_Px_Dxz_Px_bd = I_ERI_G3xy_S_Dxz_Px_bd+ABX*I_ERI_F2xy_S_Dxz_Px_bd;
  Double I_ERI_F2xz_Px_Dxz_Px_bd = I_ERI_G3xz_S_Dxz_Px_bd+ABX*I_ERI_F2xz_S_Dxz_Px_bd;
  Double I_ERI_Fx2y_Px_Dxz_Px_bd = I_ERI_G2x2y_S_Dxz_Px_bd+ABX*I_ERI_Fx2y_S_Dxz_Px_bd;
  Double I_ERI_Fxyz_Px_Dxz_Px_bd = I_ERI_G2xyz_S_Dxz_Px_bd+ABX*I_ERI_Fxyz_S_Dxz_Px_bd;
  Double I_ERI_Fx2z_Px_Dxz_Px_bd = I_ERI_G2x2z_S_Dxz_Px_bd+ABX*I_ERI_Fx2z_S_Dxz_Px_bd;
  Double I_ERI_F3y_Px_Dxz_Px_bd = I_ERI_Gx3y_S_Dxz_Px_bd+ABX*I_ERI_F3y_S_Dxz_Px_bd;
  Double I_ERI_F2yz_Px_Dxz_Px_bd = I_ERI_Gx2yz_S_Dxz_Px_bd+ABX*I_ERI_F2yz_S_Dxz_Px_bd;
  Double I_ERI_Fy2z_Px_Dxz_Px_bd = I_ERI_Gxy2z_S_Dxz_Px_bd+ABX*I_ERI_Fy2z_S_Dxz_Px_bd;
  Double I_ERI_F3z_Px_Dxz_Px_bd = I_ERI_Gx3z_S_Dxz_Px_bd+ABX*I_ERI_F3z_S_Dxz_Px_bd;
  Double I_ERI_F3x_Py_Dxz_Px_bd = I_ERI_G3xy_S_Dxz_Px_bd+ABY*I_ERI_F3x_S_Dxz_Px_bd;
  Double I_ERI_F2xy_Py_Dxz_Px_bd = I_ERI_G2x2y_S_Dxz_Px_bd+ABY*I_ERI_F2xy_S_Dxz_Px_bd;
  Double I_ERI_F2xz_Py_Dxz_Px_bd = I_ERI_G2xyz_S_Dxz_Px_bd+ABY*I_ERI_F2xz_S_Dxz_Px_bd;
  Double I_ERI_Fx2y_Py_Dxz_Px_bd = I_ERI_Gx3y_S_Dxz_Px_bd+ABY*I_ERI_Fx2y_S_Dxz_Px_bd;
  Double I_ERI_Fxyz_Py_Dxz_Px_bd = I_ERI_Gx2yz_S_Dxz_Px_bd+ABY*I_ERI_Fxyz_S_Dxz_Px_bd;
  Double I_ERI_Fx2z_Py_Dxz_Px_bd = I_ERI_Gxy2z_S_Dxz_Px_bd+ABY*I_ERI_Fx2z_S_Dxz_Px_bd;
  Double I_ERI_F3y_Py_Dxz_Px_bd = I_ERI_G4y_S_Dxz_Px_bd+ABY*I_ERI_F3y_S_Dxz_Px_bd;
  Double I_ERI_F2yz_Py_Dxz_Px_bd = I_ERI_G3yz_S_Dxz_Px_bd+ABY*I_ERI_F2yz_S_Dxz_Px_bd;
  Double I_ERI_Fy2z_Py_Dxz_Px_bd = I_ERI_G2y2z_S_Dxz_Px_bd+ABY*I_ERI_Fy2z_S_Dxz_Px_bd;
  Double I_ERI_F3z_Py_Dxz_Px_bd = I_ERI_Gy3z_S_Dxz_Px_bd+ABY*I_ERI_F3z_S_Dxz_Px_bd;
  Double I_ERI_F3x_Pz_Dxz_Px_bd = I_ERI_G3xz_S_Dxz_Px_bd+ABZ*I_ERI_F3x_S_Dxz_Px_bd;
  Double I_ERI_F2xy_Pz_Dxz_Px_bd = I_ERI_G2xyz_S_Dxz_Px_bd+ABZ*I_ERI_F2xy_S_Dxz_Px_bd;
  Double I_ERI_F2xz_Pz_Dxz_Px_bd = I_ERI_G2x2z_S_Dxz_Px_bd+ABZ*I_ERI_F2xz_S_Dxz_Px_bd;
  Double I_ERI_Fx2y_Pz_Dxz_Px_bd = I_ERI_Gx2yz_S_Dxz_Px_bd+ABZ*I_ERI_Fx2y_S_Dxz_Px_bd;
  Double I_ERI_Fxyz_Pz_Dxz_Px_bd = I_ERI_Gxy2z_S_Dxz_Px_bd+ABZ*I_ERI_Fxyz_S_Dxz_Px_bd;
  Double I_ERI_Fx2z_Pz_Dxz_Px_bd = I_ERI_Gx3z_S_Dxz_Px_bd+ABZ*I_ERI_Fx2z_S_Dxz_Px_bd;
  Double I_ERI_F3y_Pz_Dxz_Px_bd = I_ERI_G3yz_S_Dxz_Px_bd+ABZ*I_ERI_F3y_S_Dxz_Px_bd;
  Double I_ERI_F2yz_Pz_Dxz_Px_bd = I_ERI_G2y2z_S_Dxz_Px_bd+ABZ*I_ERI_F2yz_S_Dxz_Px_bd;
  Double I_ERI_Fy2z_Pz_Dxz_Px_bd = I_ERI_Gy3z_S_Dxz_Px_bd+ABZ*I_ERI_Fy2z_S_Dxz_Px_bd;
  Double I_ERI_F3z_Pz_Dxz_Px_bd = I_ERI_G4z_S_Dxz_Px_bd+ABZ*I_ERI_F3z_S_Dxz_Px_bd;
  Double I_ERI_F3x_Px_D2y_Px_bd = I_ERI_G4x_S_D2y_Px_bd+ABX*I_ERI_F3x_S_D2y_Px_bd;
  Double I_ERI_F2xy_Px_D2y_Px_bd = I_ERI_G3xy_S_D2y_Px_bd+ABX*I_ERI_F2xy_S_D2y_Px_bd;
  Double I_ERI_F2xz_Px_D2y_Px_bd = I_ERI_G3xz_S_D2y_Px_bd+ABX*I_ERI_F2xz_S_D2y_Px_bd;
  Double I_ERI_Fx2y_Px_D2y_Px_bd = I_ERI_G2x2y_S_D2y_Px_bd+ABX*I_ERI_Fx2y_S_D2y_Px_bd;
  Double I_ERI_Fxyz_Px_D2y_Px_bd = I_ERI_G2xyz_S_D2y_Px_bd+ABX*I_ERI_Fxyz_S_D2y_Px_bd;
  Double I_ERI_Fx2z_Px_D2y_Px_bd = I_ERI_G2x2z_S_D2y_Px_bd+ABX*I_ERI_Fx2z_S_D2y_Px_bd;
  Double I_ERI_F3y_Px_D2y_Px_bd = I_ERI_Gx3y_S_D2y_Px_bd+ABX*I_ERI_F3y_S_D2y_Px_bd;
  Double I_ERI_F2yz_Px_D2y_Px_bd = I_ERI_Gx2yz_S_D2y_Px_bd+ABX*I_ERI_F2yz_S_D2y_Px_bd;
  Double I_ERI_Fy2z_Px_D2y_Px_bd = I_ERI_Gxy2z_S_D2y_Px_bd+ABX*I_ERI_Fy2z_S_D2y_Px_bd;
  Double I_ERI_F3z_Px_D2y_Px_bd = I_ERI_Gx3z_S_D2y_Px_bd+ABX*I_ERI_F3z_S_D2y_Px_bd;
  Double I_ERI_F3x_Py_D2y_Px_bd = I_ERI_G3xy_S_D2y_Px_bd+ABY*I_ERI_F3x_S_D2y_Px_bd;
  Double I_ERI_F2xy_Py_D2y_Px_bd = I_ERI_G2x2y_S_D2y_Px_bd+ABY*I_ERI_F2xy_S_D2y_Px_bd;
  Double I_ERI_F2xz_Py_D2y_Px_bd = I_ERI_G2xyz_S_D2y_Px_bd+ABY*I_ERI_F2xz_S_D2y_Px_bd;
  Double I_ERI_Fx2y_Py_D2y_Px_bd = I_ERI_Gx3y_S_D2y_Px_bd+ABY*I_ERI_Fx2y_S_D2y_Px_bd;
  Double I_ERI_Fxyz_Py_D2y_Px_bd = I_ERI_Gx2yz_S_D2y_Px_bd+ABY*I_ERI_Fxyz_S_D2y_Px_bd;
  Double I_ERI_Fx2z_Py_D2y_Px_bd = I_ERI_Gxy2z_S_D2y_Px_bd+ABY*I_ERI_Fx2z_S_D2y_Px_bd;
  Double I_ERI_F3y_Py_D2y_Px_bd = I_ERI_G4y_S_D2y_Px_bd+ABY*I_ERI_F3y_S_D2y_Px_bd;
  Double I_ERI_F2yz_Py_D2y_Px_bd = I_ERI_G3yz_S_D2y_Px_bd+ABY*I_ERI_F2yz_S_D2y_Px_bd;
  Double I_ERI_Fy2z_Py_D2y_Px_bd = I_ERI_G2y2z_S_D2y_Px_bd+ABY*I_ERI_Fy2z_S_D2y_Px_bd;
  Double I_ERI_F3z_Py_D2y_Px_bd = I_ERI_Gy3z_S_D2y_Px_bd+ABY*I_ERI_F3z_S_D2y_Px_bd;
  Double I_ERI_F3x_Pz_D2y_Px_bd = I_ERI_G3xz_S_D2y_Px_bd+ABZ*I_ERI_F3x_S_D2y_Px_bd;
  Double I_ERI_F2xy_Pz_D2y_Px_bd = I_ERI_G2xyz_S_D2y_Px_bd+ABZ*I_ERI_F2xy_S_D2y_Px_bd;
  Double I_ERI_F2xz_Pz_D2y_Px_bd = I_ERI_G2x2z_S_D2y_Px_bd+ABZ*I_ERI_F2xz_S_D2y_Px_bd;
  Double I_ERI_Fx2y_Pz_D2y_Px_bd = I_ERI_Gx2yz_S_D2y_Px_bd+ABZ*I_ERI_Fx2y_S_D2y_Px_bd;
  Double I_ERI_Fxyz_Pz_D2y_Px_bd = I_ERI_Gxy2z_S_D2y_Px_bd+ABZ*I_ERI_Fxyz_S_D2y_Px_bd;
  Double I_ERI_Fx2z_Pz_D2y_Px_bd = I_ERI_Gx3z_S_D2y_Px_bd+ABZ*I_ERI_Fx2z_S_D2y_Px_bd;
  Double I_ERI_F3y_Pz_D2y_Px_bd = I_ERI_G3yz_S_D2y_Px_bd+ABZ*I_ERI_F3y_S_D2y_Px_bd;
  Double I_ERI_F2yz_Pz_D2y_Px_bd = I_ERI_G2y2z_S_D2y_Px_bd+ABZ*I_ERI_F2yz_S_D2y_Px_bd;
  Double I_ERI_Fy2z_Pz_D2y_Px_bd = I_ERI_Gy3z_S_D2y_Px_bd+ABZ*I_ERI_Fy2z_S_D2y_Px_bd;
  Double I_ERI_F3z_Pz_D2y_Px_bd = I_ERI_G4z_S_D2y_Px_bd+ABZ*I_ERI_F3z_S_D2y_Px_bd;
  Double I_ERI_F3x_Px_Dyz_Px_bd = I_ERI_G4x_S_Dyz_Px_bd+ABX*I_ERI_F3x_S_Dyz_Px_bd;
  Double I_ERI_F2xy_Px_Dyz_Px_bd = I_ERI_G3xy_S_Dyz_Px_bd+ABX*I_ERI_F2xy_S_Dyz_Px_bd;
  Double I_ERI_F2xz_Px_Dyz_Px_bd = I_ERI_G3xz_S_Dyz_Px_bd+ABX*I_ERI_F2xz_S_Dyz_Px_bd;
  Double I_ERI_Fx2y_Px_Dyz_Px_bd = I_ERI_G2x2y_S_Dyz_Px_bd+ABX*I_ERI_Fx2y_S_Dyz_Px_bd;
  Double I_ERI_Fxyz_Px_Dyz_Px_bd = I_ERI_G2xyz_S_Dyz_Px_bd+ABX*I_ERI_Fxyz_S_Dyz_Px_bd;
  Double I_ERI_Fx2z_Px_Dyz_Px_bd = I_ERI_G2x2z_S_Dyz_Px_bd+ABX*I_ERI_Fx2z_S_Dyz_Px_bd;
  Double I_ERI_F3y_Px_Dyz_Px_bd = I_ERI_Gx3y_S_Dyz_Px_bd+ABX*I_ERI_F3y_S_Dyz_Px_bd;
  Double I_ERI_F2yz_Px_Dyz_Px_bd = I_ERI_Gx2yz_S_Dyz_Px_bd+ABX*I_ERI_F2yz_S_Dyz_Px_bd;
  Double I_ERI_Fy2z_Px_Dyz_Px_bd = I_ERI_Gxy2z_S_Dyz_Px_bd+ABX*I_ERI_Fy2z_S_Dyz_Px_bd;
  Double I_ERI_F3z_Px_Dyz_Px_bd = I_ERI_Gx3z_S_Dyz_Px_bd+ABX*I_ERI_F3z_S_Dyz_Px_bd;
  Double I_ERI_F3x_Py_Dyz_Px_bd = I_ERI_G3xy_S_Dyz_Px_bd+ABY*I_ERI_F3x_S_Dyz_Px_bd;
  Double I_ERI_F2xy_Py_Dyz_Px_bd = I_ERI_G2x2y_S_Dyz_Px_bd+ABY*I_ERI_F2xy_S_Dyz_Px_bd;
  Double I_ERI_F2xz_Py_Dyz_Px_bd = I_ERI_G2xyz_S_Dyz_Px_bd+ABY*I_ERI_F2xz_S_Dyz_Px_bd;
  Double I_ERI_Fx2y_Py_Dyz_Px_bd = I_ERI_Gx3y_S_Dyz_Px_bd+ABY*I_ERI_Fx2y_S_Dyz_Px_bd;
  Double I_ERI_Fxyz_Py_Dyz_Px_bd = I_ERI_Gx2yz_S_Dyz_Px_bd+ABY*I_ERI_Fxyz_S_Dyz_Px_bd;
  Double I_ERI_Fx2z_Py_Dyz_Px_bd = I_ERI_Gxy2z_S_Dyz_Px_bd+ABY*I_ERI_Fx2z_S_Dyz_Px_bd;
  Double I_ERI_F3y_Py_Dyz_Px_bd = I_ERI_G4y_S_Dyz_Px_bd+ABY*I_ERI_F3y_S_Dyz_Px_bd;
  Double I_ERI_F2yz_Py_Dyz_Px_bd = I_ERI_G3yz_S_Dyz_Px_bd+ABY*I_ERI_F2yz_S_Dyz_Px_bd;
  Double I_ERI_Fy2z_Py_Dyz_Px_bd = I_ERI_G2y2z_S_Dyz_Px_bd+ABY*I_ERI_Fy2z_S_Dyz_Px_bd;
  Double I_ERI_F3z_Py_Dyz_Px_bd = I_ERI_Gy3z_S_Dyz_Px_bd+ABY*I_ERI_F3z_S_Dyz_Px_bd;
  Double I_ERI_F3x_Pz_Dyz_Px_bd = I_ERI_G3xz_S_Dyz_Px_bd+ABZ*I_ERI_F3x_S_Dyz_Px_bd;
  Double I_ERI_F2xy_Pz_Dyz_Px_bd = I_ERI_G2xyz_S_Dyz_Px_bd+ABZ*I_ERI_F2xy_S_Dyz_Px_bd;
  Double I_ERI_F2xz_Pz_Dyz_Px_bd = I_ERI_G2x2z_S_Dyz_Px_bd+ABZ*I_ERI_F2xz_S_Dyz_Px_bd;
  Double I_ERI_Fx2y_Pz_Dyz_Px_bd = I_ERI_Gx2yz_S_Dyz_Px_bd+ABZ*I_ERI_Fx2y_S_Dyz_Px_bd;
  Double I_ERI_Fxyz_Pz_Dyz_Px_bd = I_ERI_Gxy2z_S_Dyz_Px_bd+ABZ*I_ERI_Fxyz_S_Dyz_Px_bd;
  Double I_ERI_Fx2z_Pz_Dyz_Px_bd = I_ERI_Gx3z_S_Dyz_Px_bd+ABZ*I_ERI_Fx2z_S_Dyz_Px_bd;
  Double I_ERI_F3y_Pz_Dyz_Px_bd = I_ERI_G3yz_S_Dyz_Px_bd+ABZ*I_ERI_F3y_S_Dyz_Px_bd;
  Double I_ERI_F2yz_Pz_Dyz_Px_bd = I_ERI_G2y2z_S_Dyz_Px_bd+ABZ*I_ERI_F2yz_S_Dyz_Px_bd;
  Double I_ERI_Fy2z_Pz_Dyz_Px_bd = I_ERI_Gy3z_S_Dyz_Px_bd+ABZ*I_ERI_Fy2z_S_Dyz_Px_bd;
  Double I_ERI_F3z_Pz_Dyz_Px_bd = I_ERI_G4z_S_Dyz_Px_bd+ABZ*I_ERI_F3z_S_Dyz_Px_bd;
  Double I_ERI_F3x_Px_D2z_Px_bd = I_ERI_G4x_S_D2z_Px_bd+ABX*I_ERI_F3x_S_D2z_Px_bd;
  Double I_ERI_F2xy_Px_D2z_Px_bd = I_ERI_G3xy_S_D2z_Px_bd+ABX*I_ERI_F2xy_S_D2z_Px_bd;
  Double I_ERI_F2xz_Px_D2z_Px_bd = I_ERI_G3xz_S_D2z_Px_bd+ABX*I_ERI_F2xz_S_D2z_Px_bd;
  Double I_ERI_Fx2y_Px_D2z_Px_bd = I_ERI_G2x2y_S_D2z_Px_bd+ABX*I_ERI_Fx2y_S_D2z_Px_bd;
  Double I_ERI_Fxyz_Px_D2z_Px_bd = I_ERI_G2xyz_S_D2z_Px_bd+ABX*I_ERI_Fxyz_S_D2z_Px_bd;
  Double I_ERI_Fx2z_Px_D2z_Px_bd = I_ERI_G2x2z_S_D2z_Px_bd+ABX*I_ERI_Fx2z_S_D2z_Px_bd;
  Double I_ERI_F3y_Px_D2z_Px_bd = I_ERI_Gx3y_S_D2z_Px_bd+ABX*I_ERI_F3y_S_D2z_Px_bd;
  Double I_ERI_F2yz_Px_D2z_Px_bd = I_ERI_Gx2yz_S_D2z_Px_bd+ABX*I_ERI_F2yz_S_D2z_Px_bd;
  Double I_ERI_Fy2z_Px_D2z_Px_bd = I_ERI_Gxy2z_S_D2z_Px_bd+ABX*I_ERI_Fy2z_S_D2z_Px_bd;
  Double I_ERI_F3z_Px_D2z_Px_bd = I_ERI_Gx3z_S_D2z_Px_bd+ABX*I_ERI_F3z_S_D2z_Px_bd;
  Double I_ERI_F3x_Py_D2z_Px_bd = I_ERI_G3xy_S_D2z_Px_bd+ABY*I_ERI_F3x_S_D2z_Px_bd;
  Double I_ERI_F2xy_Py_D2z_Px_bd = I_ERI_G2x2y_S_D2z_Px_bd+ABY*I_ERI_F2xy_S_D2z_Px_bd;
  Double I_ERI_F2xz_Py_D2z_Px_bd = I_ERI_G2xyz_S_D2z_Px_bd+ABY*I_ERI_F2xz_S_D2z_Px_bd;
  Double I_ERI_Fx2y_Py_D2z_Px_bd = I_ERI_Gx3y_S_D2z_Px_bd+ABY*I_ERI_Fx2y_S_D2z_Px_bd;
  Double I_ERI_Fxyz_Py_D2z_Px_bd = I_ERI_Gx2yz_S_D2z_Px_bd+ABY*I_ERI_Fxyz_S_D2z_Px_bd;
  Double I_ERI_Fx2z_Py_D2z_Px_bd = I_ERI_Gxy2z_S_D2z_Px_bd+ABY*I_ERI_Fx2z_S_D2z_Px_bd;
  Double I_ERI_F3y_Py_D2z_Px_bd = I_ERI_G4y_S_D2z_Px_bd+ABY*I_ERI_F3y_S_D2z_Px_bd;
  Double I_ERI_F2yz_Py_D2z_Px_bd = I_ERI_G3yz_S_D2z_Px_bd+ABY*I_ERI_F2yz_S_D2z_Px_bd;
  Double I_ERI_Fy2z_Py_D2z_Px_bd = I_ERI_G2y2z_S_D2z_Px_bd+ABY*I_ERI_Fy2z_S_D2z_Px_bd;
  Double I_ERI_F3z_Py_D2z_Px_bd = I_ERI_Gy3z_S_D2z_Px_bd+ABY*I_ERI_F3z_S_D2z_Px_bd;
  Double I_ERI_F3x_Pz_D2z_Px_bd = I_ERI_G3xz_S_D2z_Px_bd+ABZ*I_ERI_F3x_S_D2z_Px_bd;
  Double I_ERI_F2xy_Pz_D2z_Px_bd = I_ERI_G2xyz_S_D2z_Px_bd+ABZ*I_ERI_F2xy_S_D2z_Px_bd;
  Double I_ERI_F2xz_Pz_D2z_Px_bd = I_ERI_G2x2z_S_D2z_Px_bd+ABZ*I_ERI_F2xz_S_D2z_Px_bd;
  Double I_ERI_Fx2y_Pz_D2z_Px_bd = I_ERI_Gx2yz_S_D2z_Px_bd+ABZ*I_ERI_Fx2y_S_D2z_Px_bd;
  Double I_ERI_Fxyz_Pz_D2z_Px_bd = I_ERI_Gxy2z_S_D2z_Px_bd+ABZ*I_ERI_Fxyz_S_D2z_Px_bd;
  Double I_ERI_Fx2z_Pz_D2z_Px_bd = I_ERI_Gx3z_S_D2z_Px_bd+ABZ*I_ERI_Fx2z_S_D2z_Px_bd;
  Double I_ERI_F3y_Pz_D2z_Px_bd = I_ERI_G3yz_S_D2z_Px_bd+ABZ*I_ERI_F3y_S_D2z_Px_bd;
  Double I_ERI_F2yz_Pz_D2z_Px_bd = I_ERI_G2y2z_S_D2z_Px_bd+ABZ*I_ERI_F2yz_S_D2z_Px_bd;
  Double I_ERI_Fy2z_Pz_D2z_Px_bd = I_ERI_Gy3z_S_D2z_Px_bd+ABZ*I_ERI_Fy2z_S_D2z_Px_bd;
  Double I_ERI_F3z_Pz_D2z_Px_bd = I_ERI_G4z_S_D2z_Px_bd+ABZ*I_ERI_F3z_S_D2z_Px_bd;
  Double I_ERI_F3x_Px_D2x_Py_bd = I_ERI_G4x_S_D2x_Py_bd+ABX*I_ERI_F3x_S_D2x_Py_bd;
  Double I_ERI_F2xy_Px_D2x_Py_bd = I_ERI_G3xy_S_D2x_Py_bd+ABX*I_ERI_F2xy_S_D2x_Py_bd;
  Double I_ERI_F2xz_Px_D2x_Py_bd = I_ERI_G3xz_S_D2x_Py_bd+ABX*I_ERI_F2xz_S_D2x_Py_bd;
  Double I_ERI_Fx2y_Px_D2x_Py_bd = I_ERI_G2x2y_S_D2x_Py_bd+ABX*I_ERI_Fx2y_S_D2x_Py_bd;
  Double I_ERI_Fxyz_Px_D2x_Py_bd = I_ERI_G2xyz_S_D2x_Py_bd+ABX*I_ERI_Fxyz_S_D2x_Py_bd;
  Double I_ERI_Fx2z_Px_D2x_Py_bd = I_ERI_G2x2z_S_D2x_Py_bd+ABX*I_ERI_Fx2z_S_D2x_Py_bd;
  Double I_ERI_F3y_Px_D2x_Py_bd = I_ERI_Gx3y_S_D2x_Py_bd+ABX*I_ERI_F3y_S_D2x_Py_bd;
  Double I_ERI_F2yz_Px_D2x_Py_bd = I_ERI_Gx2yz_S_D2x_Py_bd+ABX*I_ERI_F2yz_S_D2x_Py_bd;
  Double I_ERI_Fy2z_Px_D2x_Py_bd = I_ERI_Gxy2z_S_D2x_Py_bd+ABX*I_ERI_Fy2z_S_D2x_Py_bd;
  Double I_ERI_F3z_Px_D2x_Py_bd = I_ERI_Gx3z_S_D2x_Py_bd+ABX*I_ERI_F3z_S_D2x_Py_bd;
  Double I_ERI_F3x_Py_D2x_Py_bd = I_ERI_G3xy_S_D2x_Py_bd+ABY*I_ERI_F3x_S_D2x_Py_bd;
  Double I_ERI_F2xy_Py_D2x_Py_bd = I_ERI_G2x2y_S_D2x_Py_bd+ABY*I_ERI_F2xy_S_D2x_Py_bd;
  Double I_ERI_F2xz_Py_D2x_Py_bd = I_ERI_G2xyz_S_D2x_Py_bd+ABY*I_ERI_F2xz_S_D2x_Py_bd;
  Double I_ERI_Fx2y_Py_D2x_Py_bd = I_ERI_Gx3y_S_D2x_Py_bd+ABY*I_ERI_Fx2y_S_D2x_Py_bd;
  Double I_ERI_Fxyz_Py_D2x_Py_bd = I_ERI_Gx2yz_S_D2x_Py_bd+ABY*I_ERI_Fxyz_S_D2x_Py_bd;
  Double I_ERI_Fx2z_Py_D2x_Py_bd = I_ERI_Gxy2z_S_D2x_Py_bd+ABY*I_ERI_Fx2z_S_D2x_Py_bd;
  Double I_ERI_F3y_Py_D2x_Py_bd = I_ERI_G4y_S_D2x_Py_bd+ABY*I_ERI_F3y_S_D2x_Py_bd;
  Double I_ERI_F2yz_Py_D2x_Py_bd = I_ERI_G3yz_S_D2x_Py_bd+ABY*I_ERI_F2yz_S_D2x_Py_bd;
  Double I_ERI_Fy2z_Py_D2x_Py_bd = I_ERI_G2y2z_S_D2x_Py_bd+ABY*I_ERI_Fy2z_S_D2x_Py_bd;
  Double I_ERI_F3z_Py_D2x_Py_bd = I_ERI_Gy3z_S_D2x_Py_bd+ABY*I_ERI_F3z_S_D2x_Py_bd;
  Double I_ERI_F3x_Pz_D2x_Py_bd = I_ERI_G3xz_S_D2x_Py_bd+ABZ*I_ERI_F3x_S_D2x_Py_bd;
  Double I_ERI_F2xy_Pz_D2x_Py_bd = I_ERI_G2xyz_S_D2x_Py_bd+ABZ*I_ERI_F2xy_S_D2x_Py_bd;
  Double I_ERI_F2xz_Pz_D2x_Py_bd = I_ERI_G2x2z_S_D2x_Py_bd+ABZ*I_ERI_F2xz_S_D2x_Py_bd;
  Double I_ERI_Fx2y_Pz_D2x_Py_bd = I_ERI_Gx2yz_S_D2x_Py_bd+ABZ*I_ERI_Fx2y_S_D2x_Py_bd;
  Double I_ERI_Fxyz_Pz_D2x_Py_bd = I_ERI_Gxy2z_S_D2x_Py_bd+ABZ*I_ERI_Fxyz_S_D2x_Py_bd;
  Double I_ERI_Fx2z_Pz_D2x_Py_bd = I_ERI_Gx3z_S_D2x_Py_bd+ABZ*I_ERI_Fx2z_S_D2x_Py_bd;
  Double I_ERI_F3y_Pz_D2x_Py_bd = I_ERI_G3yz_S_D2x_Py_bd+ABZ*I_ERI_F3y_S_D2x_Py_bd;
  Double I_ERI_F2yz_Pz_D2x_Py_bd = I_ERI_G2y2z_S_D2x_Py_bd+ABZ*I_ERI_F2yz_S_D2x_Py_bd;
  Double I_ERI_Fy2z_Pz_D2x_Py_bd = I_ERI_Gy3z_S_D2x_Py_bd+ABZ*I_ERI_Fy2z_S_D2x_Py_bd;
  Double I_ERI_F3z_Pz_D2x_Py_bd = I_ERI_G4z_S_D2x_Py_bd+ABZ*I_ERI_F3z_S_D2x_Py_bd;
  Double I_ERI_F3x_Px_Dxy_Py_bd = I_ERI_G4x_S_Dxy_Py_bd+ABX*I_ERI_F3x_S_Dxy_Py_bd;
  Double I_ERI_F2xy_Px_Dxy_Py_bd = I_ERI_G3xy_S_Dxy_Py_bd+ABX*I_ERI_F2xy_S_Dxy_Py_bd;
  Double I_ERI_F2xz_Px_Dxy_Py_bd = I_ERI_G3xz_S_Dxy_Py_bd+ABX*I_ERI_F2xz_S_Dxy_Py_bd;
  Double I_ERI_Fx2y_Px_Dxy_Py_bd = I_ERI_G2x2y_S_Dxy_Py_bd+ABX*I_ERI_Fx2y_S_Dxy_Py_bd;
  Double I_ERI_Fxyz_Px_Dxy_Py_bd = I_ERI_G2xyz_S_Dxy_Py_bd+ABX*I_ERI_Fxyz_S_Dxy_Py_bd;
  Double I_ERI_Fx2z_Px_Dxy_Py_bd = I_ERI_G2x2z_S_Dxy_Py_bd+ABX*I_ERI_Fx2z_S_Dxy_Py_bd;
  Double I_ERI_F3y_Px_Dxy_Py_bd = I_ERI_Gx3y_S_Dxy_Py_bd+ABX*I_ERI_F3y_S_Dxy_Py_bd;
  Double I_ERI_F2yz_Px_Dxy_Py_bd = I_ERI_Gx2yz_S_Dxy_Py_bd+ABX*I_ERI_F2yz_S_Dxy_Py_bd;
  Double I_ERI_Fy2z_Px_Dxy_Py_bd = I_ERI_Gxy2z_S_Dxy_Py_bd+ABX*I_ERI_Fy2z_S_Dxy_Py_bd;
  Double I_ERI_F3z_Px_Dxy_Py_bd = I_ERI_Gx3z_S_Dxy_Py_bd+ABX*I_ERI_F3z_S_Dxy_Py_bd;
  Double I_ERI_F3x_Py_Dxy_Py_bd = I_ERI_G3xy_S_Dxy_Py_bd+ABY*I_ERI_F3x_S_Dxy_Py_bd;
  Double I_ERI_F2xy_Py_Dxy_Py_bd = I_ERI_G2x2y_S_Dxy_Py_bd+ABY*I_ERI_F2xy_S_Dxy_Py_bd;
  Double I_ERI_F2xz_Py_Dxy_Py_bd = I_ERI_G2xyz_S_Dxy_Py_bd+ABY*I_ERI_F2xz_S_Dxy_Py_bd;
  Double I_ERI_Fx2y_Py_Dxy_Py_bd = I_ERI_Gx3y_S_Dxy_Py_bd+ABY*I_ERI_Fx2y_S_Dxy_Py_bd;
  Double I_ERI_Fxyz_Py_Dxy_Py_bd = I_ERI_Gx2yz_S_Dxy_Py_bd+ABY*I_ERI_Fxyz_S_Dxy_Py_bd;
  Double I_ERI_Fx2z_Py_Dxy_Py_bd = I_ERI_Gxy2z_S_Dxy_Py_bd+ABY*I_ERI_Fx2z_S_Dxy_Py_bd;
  Double I_ERI_F3y_Py_Dxy_Py_bd = I_ERI_G4y_S_Dxy_Py_bd+ABY*I_ERI_F3y_S_Dxy_Py_bd;
  Double I_ERI_F2yz_Py_Dxy_Py_bd = I_ERI_G3yz_S_Dxy_Py_bd+ABY*I_ERI_F2yz_S_Dxy_Py_bd;
  Double I_ERI_Fy2z_Py_Dxy_Py_bd = I_ERI_G2y2z_S_Dxy_Py_bd+ABY*I_ERI_Fy2z_S_Dxy_Py_bd;
  Double I_ERI_F3z_Py_Dxy_Py_bd = I_ERI_Gy3z_S_Dxy_Py_bd+ABY*I_ERI_F3z_S_Dxy_Py_bd;
  Double I_ERI_F3x_Pz_Dxy_Py_bd = I_ERI_G3xz_S_Dxy_Py_bd+ABZ*I_ERI_F3x_S_Dxy_Py_bd;
  Double I_ERI_F2xy_Pz_Dxy_Py_bd = I_ERI_G2xyz_S_Dxy_Py_bd+ABZ*I_ERI_F2xy_S_Dxy_Py_bd;
  Double I_ERI_F2xz_Pz_Dxy_Py_bd = I_ERI_G2x2z_S_Dxy_Py_bd+ABZ*I_ERI_F2xz_S_Dxy_Py_bd;
  Double I_ERI_Fx2y_Pz_Dxy_Py_bd = I_ERI_Gx2yz_S_Dxy_Py_bd+ABZ*I_ERI_Fx2y_S_Dxy_Py_bd;
  Double I_ERI_Fxyz_Pz_Dxy_Py_bd = I_ERI_Gxy2z_S_Dxy_Py_bd+ABZ*I_ERI_Fxyz_S_Dxy_Py_bd;
  Double I_ERI_Fx2z_Pz_Dxy_Py_bd = I_ERI_Gx3z_S_Dxy_Py_bd+ABZ*I_ERI_Fx2z_S_Dxy_Py_bd;
  Double I_ERI_F3y_Pz_Dxy_Py_bd = I_ERI_G3yz_S_Dxy_Py_bd+ABZ*I_ERI_F3y_S_Dxy_Py_bd;
  Double I_ERI_F2yz_Pz_Dxy_Py_bd = I_ERI_G2y2z_S_Dxy_Py_bd+ABZ*I_ERI_F2yz_S_Dxy_Py_bd;
  Double I_ERI_Fy2z_Pz_Dxy_Py_bd = I_ERI_Gy3z_S_Dxy_Py_bd+ABZ*I_ERI_Fy2z_S_Dxy_Py_bd;
  Double I_ERI_F3z_Pz_Dxy_Py_bd = I_ERI_G4z_S_Dxy_Py_bd+ABZ*I_ERI_F3z_S_Dxy_Py_bd;
  Double I_ERI_F3x_Px_Dxz_Py_bd = I_ERI_G4x_S_Dxz_Py_bd+ABX*I_ERI_F3x_S_Dxz_Py_bd;
  Double I_ERI_F2xy_Px_Dxz_Py_bd = I_ERI_G3xy_S_Dxz_Py_bd+ABX*I_ERI_F2xy_S_Dxz_Py_bd;
  Double I_ERI_F2xz_Px_Dxz_Py_bd = I_ERI_G3xz_S_Dxz_Py_bd+ABX*I_ERI_F2xz_S_Dxz_Py_bd;
  Double I_ERI_Fx2y_Px_Dxz_Py_bd = I_ERI_G2x2y_S_Dxz_Py_bd+ABX*I_ERI_Fx2y_S_Dxz_Py_bd;
  Double I_ERI_Fxyz_Px_Dxz_Py_bd = I_ERI_G2xyz_S_Dxz_Py_bd+ABX*I_ERI_Fxyz_S_Dxz_Py_bd;
  Double I_ERI_Fx2z_Px_Dxz_Py_bd = I_ERI_G2x2z_S_Dxz_Py_bd+ABX*I_ERI_Fx2z_S_Dxz_Py_bd;
  Double I_ERI_F3y_Px_Dxz_Py_bd = I_ERI_Gx3y_S_Dxz_Py_bd+ABX*I_ERI_F3y_S_Dxz_Py_bd;
  Double I_ERI_F2yz_Px_Dxz_Py_bd = I_ERI_Gx2yz_S_Dxz_Py_bd+ABX*I_ERI_F2yz_S_Dxz_Py_bd;
  Double I_ERI_Fy2z_Px_Dxz_Py_bd = I_ERI_Gxy2z_S_Dxz_Py_bd+ABX*I_ERI_Fy2z_S_Dxz_Py_bd;
  Double I_ERI_F3z_Px_Dxz_Py_bd = I_ERI_Gx3z_S_Dxz_Py_bd+ABX*I_ERI_F3z_S_Dxz_Py_bd;
  Double I_ERI_F3x_Py_Dxz_Py_bd = I_ERI_G3xy_S_Dxz_Py_bd+ABY*I_ERI_F3x_S_Dxz_Py_bd;
  Double I_ERI_F2xy_Py_Dxz_Py_bd = I_ERI_G2x2y_S_Dxz_Py_bd+ABY*I_ERI_F2xy_S_Dxz_Py_bd;
  Double I_ERI_F2xz_Py_Dxz_Py_bd = I_ERI_G2xyz_S_Dxz_Py_bd+ABY*I_ERI_F2xz_S_Dxz_Py_bd;
  Double I_ERI_Fx2y_Py_Dxz_Py_bd = I_ERI_Gx3y_S_Dxz_Py_bd+ABY*I_ERI_Fx2y_S_Dxz_Py_bd;
  Double I_ERI_Fxyz_Py_Dxz_Py_bd = I_ERI_Gx2yz_S_Dxz_Py_bd+ABY*I_ERI_Fxyz_S_Dxz_Py_bd;
  Double I_ERI_Fx2z_Py_Dxz_Py_bd = I_ERI_Gxy2z_S_Dxz_Py_bd+ABY*I_ERI_Fx2z_S_Dxz_Py_bd;
  Double I_ERI_F3y_Py_Dxz_Py_bd = I_ERI_G4y_S_Dxz_Py_bd+ABY*I_ERI_F3y_S_Dxz_Py_bd;
  Double I_ERI_F2yz_Py_Dxz_Py_bd = I_ERI_G3yz_S_Dxz_Py_bd+ABY*I_ERI_F2yz_S_Dxz_Py_bd;
  Double I_ERI_Fy2z_Py_Dxz_Py_bd = I_ERI_G2y2z_S_Dxz_Py_bd+ABY*I_ERI_Fy2z_S_Dxz_Py_bd;
  Double I_ERI_F3z_Py_Dxz_Py_bd = I_ERI_Gy3z_S_Dxz_Py_bd+ABY*I_ERI_F3z_S_Dxz_Py_bd;
  Double I_ERI_F3x_Pz_Dxz_Py_bd = I_ERI_G3xz_S_Dxz_Py_bd+ABZ*I_ERI_F3x_S_Dxz_Py_bd;
  Double I_ERI_F2xy_Pz_Dxz_Py_bd = I_ERI_G2xyz_S_Dxz_Py_bd+ABZ*I_ERI_F2xy_S_Dxz_Py_bd;
  Double I_ERI_F2xz_Pz_Dxz_Py_bd = I_ERI_G2x2z_S_Dxz_Py_bd+ABZ*I_ERI_F2xz_S_Dxz_Py_bd;
  Double I_ERI_Fx2y_Pz_Dxz_Py_bd = I_ERI_Gx2yz_S_Dxz_Py_bd+ABZ*I_ERI_Fx2y_S_Dxz_Py_bd;
  Double I_ERI_Fxyz_Pz_Dxz_Py_bd = I_ERI_Gxy2z_S_Dxz_Py_bd+ABZ*I_ERI_Fxyz_S_Dxz_Py_bd;
  Double I_ERI_Fx2z_Pz_Dxz_Py_bd = I_ERI_Gx3z_S_Dxz_Py_bd+ABZ*I_ERI_Fx2z_S_Dxz_Py_bd;
  Double I_ERI_F3y_Pz_Dxz_Py_bd = I_ERI_G3yz_S_Dxz_Py_bd+ABZ*I_ERI_F3y_S_Dxz_Py_bd;
  Double I_ERI_F2yz_Pz_Dxz_Py_bd = I_ERI_G2y2z_S_Dxz_Py_bd+ABZ*I_ERI_F2yz_S_Dxz_Py_bd;
  Double I_ERI_Fy2z_Pz_Dxz_Py_bd = I_ERI_Gy3z_S_Dxz_Py_bd+ABZ*I_ERI_Fy2z_S_Dxz_Py_bd;
  Double I_ERI_F3z_Pz_Dxz_Py_bd = I_ERI_G4z_S_Dxz_Py_bd+ABZ*I_ERI_F3z_S_Dxz_Py_bd;
  Double I_ERI_F3x_Px_D2y_Py_bd = I_ERI_G4x_S_D2y_Py_bd+ABX*I_ERI_F3x_S_D2y_Py_bd;
  Double I_ERI_F2xy_Px_D2y_Py_bd = I_ERI_G3xy_S_D2y_Py_bd+ABX*I_ERI_F2xy_S_D2y_Py_bd;
  Double I_ERI_F2xz_Px_D2y_Py_bd = I_ERI_G3xz_S_D2y_Py_bd+ABX*I_ERI_F2xz_S_D2y_Py_bd;
  Double I_ERI_Fx2y_Px_D2y_Py_bd = I_ERI_G2x2y_S_D2y_Py_bd+ABX*I_ERI_Fx2y_S_D2y_Py_bd;
  Double I_ERI_Fxyz_Px_D2y_Py_bd = I_ERI_G2xyz_S_D2y_Py_bd+ABX*I_ERI_Fxyz_S_D2y_Py_bd;
  Double I_ERI_Fx2z_Px_D2y_Py_bd = I_ERI_G2x2z_S_D2y_Py_bd+ABX*I_ERI_Fx2z_S_D2y_Py_bd;
  Double I_ERI_F3y_Px_D2y_Py_bd = I_ERI_Gx3y_S_D2y_Py_bd+ABX*I_ERI_F3y_S_D2y_Py_bd;
  Double I_ERI_F2yz_Px_D2y_Py_bd = I_ERI_Gx2yz_S_D2y_Py_bd+ABX*I_ERI_F2yz_S_D2y_Py_bd;
  Double I_ERI_Fy2z_Px_D2y_Py_bd = I_ERI_Gxy2z_S_D2y_Py_bd+ABX*I_ERI_Fy2z_S_D2y_Py_bd;
  Double I_ERI_F3z_Px_D2y_Py_bd = I_ERI_Gx3z_S_D2y_Py_bd+ABX*I_ERI_F3z_S_D2y_Py_bd;
  Double I_ERI_F3x_Py_D2y_Py_bd = I_ERI_G3xy_S_D2y_Py_bd+ABY*I_ERI_F3x_S_D2y_Py_bd;
  Double I_ERI_F2xy_Py_D2y_Py_bd = I_ERI_G2x2y_S_D2y_Py_bd+ABY*I_ERI_F2xy_S_D2y_Py_bd;
  Double I_ERI_F2xz_Py_D2y_Py_bd = I_ERI_G2xyz_S_D2y_Py_bd+ABY*I_ERI_F2xz_S_D2y_Py_bd;
  Double I_ERI_Fx2y_Py_D2y_Py_bd = I_ERI_Gx3y_S_D2y_Py_bd+ABY*I_ERI_Fx2y_S_D2y_Py_bd;
  Double I_ERI_Fxyz_Py_D2y_Py_bd = I_ERI_Gx2yz_S_D2y_Py_bd+ABY*I_ERI_Fxyz_S_D2y_Py_bd;
  Double I_ERI_Fx2z_Py_D2y_Py_bd = I_ERI_Gxy2z_S_D2y_Py_bd+ABY*I_ERI_Fx2z_S_D2y_Py_bd;
  Double I_ERI_F3y_Py_D2y_Py_bd = I_ERI_G4y_S_D2y_Py_bd+ABY*I_ERI_F3y_S_D2y_Py_bd;
  Double I_ERI_F2yz_Py_D2y_Py_bd = I_ERI_G3yz_S_D2y_Py_bd+ABY*I_ERI_F2yz_S_D2y_Py_bd;
  Double I_ERI_Fy2z_Py_D2y_Py_bd = I_ERI_G2y2z_S_D2y_Py_bd+ABY*I_ERI_Fy2z_S_D2y_Py_bd;
  Double I_ERI_F3z_Py_D2y_Py_bd = I_ERI_Gy3z_S_D2y_Py_bd+ABY*I_ERI_F3z_S_D2y_Py_bd;
  Double I_ERI_F3x_Pz_D2y_Py_bd = I_ERI_G3xz_S_D2y_Py_bd+ABZ*I_ERI_F3x_S_D2y_Py_bd;
  Double I_ERI_F2xy_Pz_D2y_Py_bd = I_ERI_G2xyz_S_D2y_Py_bd+ABZ*I_ERI_F2xy_S_D2y_Py_bd;
  Double I_ERI_F2xz_Pz_D2y_Py_bd = I_ERI_G2x2z_S_D2y_Py_bd+ABZ*I_ERI_F2xz_S_D2y_Py_bd;
  Double I_ERI_Fx2y_Pz_D2y_Py_bd = I_ERI_Gx2yz_S_D2y_Py_bd+ABZ*I_ERI_Fx2y_S_D2y_Py_bd;
  Double I_ERI_Fxyz_Pz_D2y_Py_bd = I_ERI_Gxy2z_S_D2y_Py_bd+ABZ*I_ERI_Fxyz_S_D2y_Py_bd;
  Double I_ERI_Fx2z_Pz_D2y_Py_bd = I_ERI_Gx3z_S_D2y_Py_bd+ABZ*I_ERI_Fx2z_S_D2y_Py_bd;
  Double I_ERI_F3y_Pz_D2y_Py_bd = I_ERI_G3yz_S_D2y_Py_bd+ABZ*I_ERI_F3y_S_D2y_Py_bd;
  Double I_ERI_F2yz_Pz_D2y_Py_bd = I_ERI_G2y2z_S_D2y_Py_bd+ABZ*I_ERI_F2yz_S_D2y_Py_bd;
  Double I_ERI_Fy2z_Pz_D2y_Py_bd = I_ERI_Gy3z_S_D2y_Py_bd+ABZ*I_ERI_Fy2z_S_D2y_Py_bd;
  Double I_ERI_F3z_Pz_D2y_Py_bd = I_ERI_G4z_S_D2y_Py_bd+ABZ*I_ERI_F3z_S_D2y_Py_bd;
  Double I_ERI_F3x_Px_Dyz_Py_bd = I_ERI_G4x_S_Dyz_Py_bd+ABX*I_ERI_F3x_S_Dyz_Py_bd;
  Double I_ERI_F2xy_Px_Dyz_Py_bd = I_ERI_G3xy_S_Dyz_Py_bd+ABX*I_ERI_F2xy_S_Dyz_Py_bd;
  Double I_ERI_F2xz_Px_Dyz_Py_bd = I_ERI_G3xz_S_Dyz_Py_bd+ABX*I_ERI_F2xz_S_Dyz_Py_bd;
  Double I_ERI_Fx2y_Px_Dyz_Py_bd = I_ERI_G2x2y_S_Dyz_Py_bd+ABX*I_ERI_Fx2y_S_Dyz_Py_bd;
  Double I_ERI_Fxyz_Px_Dyz_Py_bd = I_ERI_G2xyz_S_Dyz_Py_bd+ABX*I_ERI_Fxyz_S_Dyz_Py_bd;
  Double I_ERI_Fx2z_Px_Dyz_Py_bd = I_ERI_G2x2z_S_Dyz_Py_bd+ABX*I_ERI_Fx2z_S_Dyz_Py_bd;
  Double I_ERI_F3y_Px_Dyz_Py_bd = I_ERI_Gx3y_S_Dyz_Py_bd+ABX*I_ERI_F3y_S_Dyz_Py_bd;
  Double I_ERI_F2yz_Px_Dyz_Py_bd = I_ERI_Gx2yz_S_Dyz_Py_bd+ABX*I_ERI_F2yz_S_Dyz_Py_bd;
  Double I_ERI_Fy2z_Px_Dyz_Py_bd = I_ERI_Gxy2z_S_Dyz_Py_bd+ABX*I_ERI_Fy2z_S_Dyz_Py_bd;
  Double I_ERI_F3z_Px_Dyz_Py_bd = I_ERI_Gx3z_S_Dyz_Py_bd+ABX*I_ERI_F3z_S_Dyz_Py_bd;
  Double I_ERI_F3x_Py_Dyz_Py_bd = I_ERI_G3xy_S_Dyz_Py_bd+ABY*I_ERI_F3x_S_Dyz_Py_bd;
  Double I_ERI_F2xy_Py_Dyz_Py_bd = I_ERI_G2x2y_S_Dyz_Py_bd+ABY*I_ERI_F2xy_S_Dyz_Py_bd;
  Double I_ERI_F2xz_Py_Dyz_Py_bd = I_ERI_G2xyz_S_Dyz_Py_bd+ABY*I_ERI_F2xz_S_Dyz_Py_bd;
  Double I_ERI_Fx2y_Py_Dyz_Py_bd = I_ERI_Gx3y_S_Dyz_Py_bd+ABY*I_ERI_Fx2y_S_Dyz_Py_bd;
  Double I_ERI_Fxyz_Py_Dyz_Py_bd = I_ERI_Gx2yz_S_Dyz_Py_bd+ABY*I_ERI_Fxyz_S_Dyz_Py_bd;
  Double I_ERI_Fx2z_Py_Dyz_Py_bd = I_ERI_Gxy2z_S_Dyz_Py_bd+ABY*I_ERI_Fx2z_S_Dyz_Py_bd;
  Double I_ERI_F3y_Py_Dyz_Py_bd = I_ERI_G4y_S_Dyz_Py_bd+ABY*I_ERI_F3y_S_Dyz_Py_bd;
  Double I_ERI_F2yz_Py_Dyz_Py_bd = I_ERI_G3yz_S_Dyz_Py_bd+ABY*I_ERI_F2yz_S_Dyz_Py_bd;
  Double I_ERI_Fy2z_Py_Dyz_Py_bd = I_ERI_G2y2z_S_Dyz_Py_bd+ABY*I_ERI_Fy2z_S_Dyz_Py_bd;
  Double I_ERI_F3z_Py_Dyz_Py_bd = I_ERI_Gy3z_S_Dyz_Py_bd+ABY*I_ERI_F3z_S_Dyz_Py_bd;
  Double I_ERI_F3x_Pz_Dyz_Py_bd = I_ERI_G3xz_S_Dyz_Py_bd+ABZ*I_ERI_F3x_S_Dyz_Py_bd;
  Double I_ERI_F2xy_Pz_Dyz_Py_bd = I_ERI_G2xyz_S_Dyz_Py_bd+ABZ*I_ERI_F2xy_S_Dyz_Py_bd;
  Double I_ERI_F2xz_Pz_Dyz_Py_bd = I_ERI_G2x2z_S_Dyz_Py_bd+ABZ*I_ERI_F2xz_S_Dyz_Py_bd;
  Double I_ERI_Fx2y_Pz_Dyz_Py_bd = I_ERI_Gx2yz_S_Dyz_Py_bd+ABZ*I_ERI_Fx2y_S_Dyz_Py_bd;
  Double I_ERI_Fxyz_Pz_Dyz_Py_bd = I_ERI_Gxy2z_S_Dyz_Py_bd+ABZ*I_ERI_Fxyz_S_Dyz_Py_bd;
  Double I_ERI_Fx2z_Pz_Dyz_Py_bd = I_ERI_Gx3z_S_Dyz_Py_bd+ABZ*I_ERI_Fx2z_S_Dyz_Py_bd;
  Double I_ERI_F3y_Pz_Dyz_Py_bd = I_ERI_G3yz_S_Dyz_Py_bd+ABZ*I_ERI_F3y_S_Dyz_Py_bd;
  Double I_ERI_F2yz_Pz_Dyz_Py_bd = I_ERI_G2y2z_S_Dyz_Py_bd+ABZ*I_ERI_F2yz_S_Dyz_Py_bd;
  Double I_ERI_Fy2z_Pz_Dyz_Py_bd = I_ERI_Gy3z_S_Dyz_Py_bd+ABZ*I_ERI_Fy2z_S_Dyz_Py_bd;
  Double I_ERI_F3z_Pz_Dyz_Py_bd = I_ERI_G4z_S_Dyz_Py_bd+ABZ*I_ERI_F3z_S_Dyz_Py_bd;
  Double I_ERI_F3x_Px_D2z_Py_bd = I_ERI_G4x_S_D2z_Py_bd+ABX*I_ERI_F3x_S_D2z_Py_bd;
  Double I_ERI_F2xy_Px_D2z_Py_bd = I_ERI_G3xy_S_D2z_Py_bd+ABX*I_ERI_F2xy_S_D2z_Py_bd;
  Double I_ERI_F2xz_Px_D2z_Py_bd = I_ERI_G3xz_S_D2z_Py_bd+ABX*I_ERI_F2xz_S_D2z_Py_bd;
  Double I_ERI_Fx2y_Px_D2z_Py_bd = I_ERI_G2x2y_S_D2z_Py_bd+ABX*I_ERI_Fx2y_S_D2z_Py_bd;
  Double I_ERI_Fxyz_Px_D2z_Py_bd = I_ERI_G2xyz_S_D2z_Py_bd+ABX*I_ERI_Fxyz_S_D2z_Py_bd;
  Double I_ERI_Fx2z_Px_D2z_Py_bd = I_ERI_G2x2z_S_D2z_Py_bd+ABX*I_ERI_Fx2z_S_D2z_Py_bd;
  Double I_ERI_F3y_Px_D2z_Py_bd = I_ERI_Gx3y_S_D2z_Py_bd+ABX*I_ERI_F3y_S_D2z_Py_bd;
  Double I_ERI_F2yz_Px_D2z_Py_bd = I_ERI_Gx2yz_S_D2z_Py_bd+ABX*I_ERI_F2yz_S_D2z_Py_bd;
  Double I_ERI_Fy2z_Px_D2z_Py_bd = I_ERI_Gxy2z_S_D2z_Py_bd+ABX*I_ERI_Fy2z_S_D2z_Py_bd;
  Double I_ERI_F3z_Px_D2z_Py_bd = I_ERI_Gx3z_S_D2z_Py_bd+ABX*I_ERI_F3z_S_D2z_Py_bd;
  Double I_ERI_F3x_Py_D2z_Py_bd = I_ERI_G3xy_S_D2z_Py_bd+ABY*I_ERI_F3x_S_D2z_Py_bd;
  Double I_ERI_F2xy_Py_D2z_Py_bd = I_ERI_G2x2y_S_D2z_Py_bd+ABY*I_ERI_F2xy_S_D2z_Py_bd;
  Double I_ERI_F2xz_Py_D2z_Py_bd = I_ERI_G2xyz_S_D2z_Py_bd+ABY*I_ERI_F2xz_S_D2z_Py_bd;
  Double I_ERI_Fx2y_Py_D2z_Py_bd = I_ERI_Gx3y_S_D2z_Py_bd+ABY*I_ERI_Fx2y_S_D2z_Py_bd;
  Double I_ERI_Fxyz_Py_D2z_Py_bd = I_ERI_Gx2yz_S_D2z_Py_bd+ABY*I_ERI_Fxyz_S_D2z_Py_bd;
  Double I_ERI_Fx2z_Py_D2z_Py_bd = I_ERI_Gxy2z_S_D2z_Py_bd+ABY*I_ERI_Fx2z_S_D2z_Py_bd;
  Double I_ERI_F3y_Py_D2z_Py_bd = I_ERI_G4y_S_D2z_Py_bd+ABY*I_ERI_F3y_S_D2z_Py_bd;
  Double I_ERI_F2yz_Py_D2z_Py_bd = I_ERI_G3yz_S_D2z_Py_bd+ABY*I_ERI_F2yz_S_D2z_Py_bd;
  Double I_ERI_Fy2z_Py_D2z_Py_bd = I_ERI_G2y2z_S_D2z_Py_bd+ABY*I_ERI_Fy2z_S_D2z_Py_bd;
  Double I_ERI_F3z_Py_D2z_Py_bd = I_ERI_Gy3z_S_D2z_Py_bd+ABY*I_ERI_F3z_S_D2z_Py_bd;
  Double I_ERI_F3x_Pz_D2z_Py_bd = I_ERI_G3xz_S_D2z_Py_bd+ABZ*I_ERI_F3x_S_D2z_Py_bd;
  Double I_ERI_F2xy_Pz_D2z_Py_bd = I_ERI_G2xyz_S_D2z_Py_bd+ABZ*I_ERI_F2xy_S_D2z_Py_bd;
  Double I_ERI_F2xz_Pz_D2z_Py_bd = I_ERI_G2x2z_S_D2z_Py_bd+ABZ*I_ERI_F2xz_S_D2z_Py_bd;
  Double I_ERI_Fx2y_Pz_D2z_Py_bd = I_ERI_Gx2yz_S_D2z_Py_bd+ABZ*I_ERI_Fx2y_S_D2z_Py_bd;
  Double I_ERI_Fxyz_Pz_D2z_Py_bd = I_ERI_Gxy2z_S_D2z_Py_bd+ABZ*I_ERI_Fxyz_S_D2z_Py_bd;
  Double I_ERI_Fx2z_Pz_D2z_Py_bd = I_ERI_Gx3z_S_D2z_Py_bd+ABZ*I_ERI_Fx2z_S_D2z_Py_bd;
  Double I_ERI_F3y_Pz_D2z_Py_bd = I_ERI_G3yz_S_D2z_Py_bd+ABZ*I_ERI_F3y_S_D2z_Py_bd;
  Double I_ERI_F2yz_Pz_D2z_Py_bd = I_ERI_G2y2z_S_D2z_Py_bd+ABZ*I_ERI_F2yz_S_D2z_Py_bd;
  Double I_ERI_Fy2z_Pz_D2z_Py_bd = I_ERI_Gy3z_S_D2z_Py_bd+ABZ*I_ERI_Fy2z_S_D2z_Py_bd;
  Double I_ERI_F3z_Pz_D2z_Py_bd = I_ERI_G4z_S_D2z_Py_bd+ABZ*I_ERI_F3z_S_D2z_Py_bd;
  Double I_ERI_F3x_Px_D2x_Pz_bd = I_ERI_G4x_S_D2x_Pz_bd+ABX*I_ERI_F3x_S_D2x_Pz_bd;
  Double I_ERI_F2xy_Px_D2x_Pz_bd = I_ERI_G3xy_S_D2x_Pz_bd+ABX*I_ERI_F2xy_S_D2x_Pz_bd;
  Double I_ERI_F2xz_Px_D2x_Pz_bd = I_ERI_G3xz_S_D2x_Pz_bd+ABX*I_ERI_F2xz_S_D2x_Pz_bd;
  Double I_ERI_Fx2y_Px_D2x_Pz_bd = I_ERI_G2x2y_S_D2x_Pz_bd+ABX*I_ERI_Fx2y_S_D2x_Pz_bd;
  Double I_ERI_Fxyz_Px_D2x_Pz_bd = I_ERI_G2xyz_S_D2x_Pz_bd+ABX*I_ERI_Fxyz_S_D2x_Pz_bd;
  Double I_ERI_Fx2z_Px_D2x_Pz_bd = I_ERI_G2x2z_S_D2x_Pz_bd+ABX*I_ERI_Fx2z_S_D2x_Pz_bd;
  Double I_ERI_F3y_Px_D2x_Pz_bd = I_ERI_Gx3y_S_D2x_Pz_bd+ABX*I_ERI_F3y_S_D2x_Pz_bd;
  Double I_ERI_F2yz_Px_D2x_Pz_bd = I_ERI_Gx2yz_S_D2x_Pz_bd+ABX*I_ERI_F2yz_S_D2x_Pz_bd;
  Double I_ERI_Fy2z_Px_D2x_Pz_bd = I_ERI_Gxy2z_S_D2x_Pz_bd+ABX*I_ERI_Fy2z_S_D2x_Pz_bd;
  Double I_ERI_F3z_Px_D2x_Pz_bd = I_ERI_Gx3z_S_D2x_Pz_bd+ABX*I_ERI_F3z_S_D2x_Pz_bd;
  Double I_ERI_F3x_Py_D2x_Pz_bd = I_ERI_G3xy_S_D2x_Pz_bd+ABY*I_ERI_F3x_S_D2x_Pz_bd;
  Double I_ERI_F2xy_Py_D2x_Pz_bd = I_ERI_G2x2y_S_D2x_Pz_bd+ABY*I_ERI_F2xy_S_D2x_Pz_bd;
  Double I_ERI_F2xz_Py_D2x_Pz_bd = I_ERI_G2xyz_S_D2x_Pz_bd+ABY*I_ERI_F2xz_S_D2x_Pz_bd;
  Double I_ERI_Fx2y_Py_D2x_Pz_bd = I_ERI_Gx3y_S_D2x_Pz_bd+ABY*I_ERI_Fx2y_S_D2x_Pz_bd;
  Double I_ERI_Fxyz_Py_D2x_Pz_bd = I_ERI_Gx2yz_S_D2x_Pz_bd+ABY*I_ERI_Fxyz_S_D2x_Pz_bd;
  Double I_ERI_Fx2z_Py_D2x_Pz_bd = I_ERI_Gxy2z_S_D2x_Pz_bd+ABY*I_ERI_Fx2z_S_D2x_Pz_bd;
  Double I_ERI_F3y_Py_D2x_Pz_bd = I_ERI_G4y_S_D2x_Pz_bd+ABY*I_ERI_F3y_S_D2x_Pz_bd;
  Double I_ERI_F2yz_Py_D2x_Pz_bd = I_ERI_G3yz_S_D2x_Pz_bd+ABY*I_ERI_F2yz_S_D2x_Pz_bd;
  Double I_ERI_Fy2z_Py_D2x_Pz_bd = I_ERI_G2y2z_S_D2x_Pz_bd+ABY*I_ERI_Fy2z_S_D2x_Pz_bd;
  Double I_ERI_F3z_Py_D2x_Pz_bd = I_ERI_Gy3z_S_D2x_Pz_bd+ABY*I_ERI_F3z_S_D2x_Pz_bd;
  Double I_ERI_F3x_Pz_D2x_Pz_bd = I_ERI_G3xz_S_D2x_Pz_bd+ABZ*I_ERI_F3x_S_D2x_Pz_bd;
  Double I_ERI_F2xy_Pz_D2x_Pz_bd = I_ERI_G2xyz_S_D2x_Pz_bd+ABZ*I_ERI_F2xy_S_D2x_Pz_bd;
  Double I_ERI_F2xz_Pz_D2x_Pz_bd = I_ERI_G2x2z_S_D2x_Pz_bd+ABZ*I_ERI_F2xz_S_D2x_Pz_bd;
  Double I_ERI_Fx2y_Pz_D2x_Pz_bd = I_ERI_Gx2yz_S_D2x_Pz_bd+ABZ*I_ERI_Fx2y_S_D2x_Pz_bd;
  Double I_ERI_Fxyz_Pz_D2x_Pz_bd = I_ERI_Gxy2z_S_D2x_Pz_bd+ABZ*I_ERI_Fxyz_S_D2x_Pz_bd;
  Double I_ERI_Fx2z_Pz_D2x_Pz_bd = I_ERI_Gx3z_S_D2x_Pz_bd+ABZ*I_ERI_Fx2z_S_D2x_Pz_bd;
  Double I_ERI_F3y_Pz_D2x_Pz_bd = I_ERI_G3yz_S_D2x_Pz_bd+ABZ*I_ERI_F3y_S_D2x_Pz_bd;
  Double I_ERI_F2yz_Pz_D2x_Pz_bd = I_ERI_G2y2z_S_D2x_Pz_bd+ABZ*I_ERI_F2yz_S_D2x_Pz_bd;
  Double I_ERI_Fy2z_Pz_D2x_Pz_bd = I_ERI_Gy3z_S_D2x_Pz_bd+ABZ*I_ERI_Fy2z_S_D2x_Pz_bd;
  Double I_ERI_F3z_Pz_D2x_Pz_bd = I_ERI_G4z_S_D2x_Pz_bd+ABZ*I_ERI_F3z_S_D2x_Pz_bd;
  Double I_ERI_F3x_Px_Dxy_Pz_bd = I_ERI_G4x_S_Dxy_Pz_bd+ABX*I_ERI_F3x_S_Dxy_Pz_bd;
  Double I_ERI_F2xy_Px_Dxy_Pz_bd = I_ERI_G3xy_S_Dxy_Pz_bd+ABX*I_ERI_F2xy_S_Dxy_Pz_bd;
  Double I_ERI_F2xz_Px_Dxy_Pz_bd = I_ERI_G3xz_S_Dxy_Pz_bd+ABX*I_ERI_F2xz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2y_Px_Dxy_Pz_bd = I_ERI_G2x2y_S_Dxy_Pz_bd+ABX*I_ERI_Fx2y_S_Dxy_Pz_bd;
  Double I_ERI_Fxyz_Px_Dxy_Pz_bd = I_ERI_G2xyz_S_Dxy_Pz_bd+ABX*I_ERI_Fxyz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2z_Px_Dxy_Pz_bd = I_ERI_G2x2z_S_Dxy_Pz_bd+ABX*I_ERI_Fx2z_S_Dxy_Pz_bd;
  Double I_ERI_F3y_Px_Dxy_Pz_bd = I_ERI_Gx3y_S_Dxy_Pz_bd+ABX*I_ERI_F3y_S_Dxy_Pz_bd;
  Double I_ERI_F2yz_Px_Dxy_Pz_bd = I_ERI_Gx2yz_S_Dxy_Pz_bd+ABX*I_ERI_F2yz_S_Dxy_Pz_bd;
  Double I_ERI_Fy2z_Px_Dxy_Pz_bd = I_ERI_Gxy2z_S_Dxy_Pz_bd+ABX*I_ERI_Fy2z_S_Dxy_Pz_bd;
  Double I_ERI_F3z_Px_Dxy_Pz_bd = I_ERI_Gx3z_S_Dxy_Pz_bd+ABX*I_ERI_F3z_S_Dxy_Pz_bd;
  Double I_ERI_F3x_Py_Dxy_Pz_bd = I_ERI_G3xy_S_Dxy_Pz_bd+ABY*I_ERI_F3x_S_Dxy_Pz_bd;
  Double I_ERI_F2xy_Py_Dxy_Pz_bd = I_ERI_G2x2y_S_Dxy_Pz_bd+ABY*I_ERI_F2xy_S_Dxy_Pz_bd;
  Double I_ERI_F2xz_Py_Dxy_Pz_bd = I_ERI_G2xyz_S_Dxy_Pz_bd+ABY*I_ERI_F2xz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2y_Py_Dxy_Pz_bd = I_ERI_Gx3y_S_Dxy_Pz_bd+ABY*I_ERI_Fx2y_S_Dxy_Pz_bd;
  Double I_ERI_Fxyz_Py_Dxy_Pz_bd = I_ERI_Gx2yz_S_Dxy_Pz_bd+ABY*I_ERI_Fxyz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2z_Py_Dxy_Pz_bd = I_ERI_Gxy2z_S_Dxy_Pz_bd+ABY*I_ERI_Fx2z_S_Dxy_Pz_bd;
  Double I_ERI_F3y_Py_Dxy_Pz_bd = I_ERI_G4y_S_Dxy_Pz_bd+ABY*I_ERI_F3y_S_Dxy_Pz_bd;
  Double I_ERI_F2yz_Py_Dxy_Pz_bd = I_ERI_G3yz_S_Dxy_Pz_bd+ABY*I_ERI_F2yz_S_Dxy_Pz_bd;
  Double I_ERI_Fy2z_Py_Dxy_Pz_bd = I_ERI_G2y2z_S_Dxy_Pz_bd+ABY*I_ERI_Fy2z_S_Dxy_Pz_bd;
  Double I_ERI_F3z_Py_Dxy_Pz_bd = I_ERI_Gy3z_S_Dxy_Pz_bd+ABY*I_ERI_F3z_S_Dxy_Pz_bd;
  Double I_ERI_F3x_Pz_Dxy_Pz_bd = I_ERI_G3xz_S_Dxy_Pz_bd+ABZ*I_ERI_F3x_S_Dxy_Pz_bd;
  Double I_ERI_F2xy_Pz_Dxy_Pz_bd = I_ERI_G2xyz_S_Dxy_Pz_bd+ABZ*I_ERI_F2xy_S_Dxy_Pz_bd;
  Double I_ERI_F2xz_Pz_Dxy_Pz_bd = I_ERI_G2x2z_S_Dxy_Pz_bd+ABZ*I_ERI_F2xz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2y_Pz_Dxy_Pz_bd = I_ERI_Gx2yz_S_Dxy_Pz_bd+ABZ*I_ERI_Fx2y_S_Dxy_Pz_bd;
  Double I_ERI_Fxyz_Pz_Dxy_Pz_bd = I_ERI_Gxy2z_S_Dxy_Pz_bd+ABZ*I_ERI_Fxyz_S_Dxy_Pz_bd;
  Double I_ERI_Fx2z_Pz_Dxy_Pz_bd = I_ERI_Gx3z_S_Dxy_Pz_bd+ABZ*I_ERI_Fx2z_S_Dxy_Pz_bd;
  Double I_ERI_F3y_Pz_Dxy_Pz_bd = I_ERI_G3yz_S_Dxy_Pz_bd+ABZ*I_ERI_F3y_S_Dxy_Pz_bd;
  Double I_ERI_F2yz_Pz_Dxy_Pz_bd = I_ERI_G2y2z_S_Dxy_Pz_bd+ABZ*I_ERI_F2yz_S_Dxy_Pz_bd;
  Double I_ERI_Fy2z_Pz_Dxy_Pz_bd = I_ERI_Gy3z_S_Dxy_Pz_bd+ABZ*I_ERI_Fy2z_S_Dxy_Pz_bd;
  Double I_ERI_F3z_Pz_Dxy_Pz_bd = I_ERI_G4z_S_Dxy_Pz_bd+ABZ*I_ERI_F3z_S_Dxy_Pz_bd;
  Double I_ERI_F3x_Px_Dxz_Pz_bd = I_ERI_G4x_S_Dxz_Pz_bd+ABX*I_ERI_F3x_S_Dxz_Pz_bd;
  Double I_ERI_F2xy_Px_Dxz_Pz_bd = I_ERI_G3xy_S_Dxz_Pz_bd+ABX*I_ERI_F2xy_S_Dxz_Pz_bd;
  Double I_ERI_F2xz_Px_Dxz_Pz_bd = I_ERI_G3xz_S_Dxz_Pz_bd+ABX*I_ERI_F2xz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2y_Px_Dxz_Pz_bd = I_ERI_G2x2y_S_Dxz_Pz_bd+ABX*I_ERI_Fx2y_S_Dxz_Pz_bd;
  Double I_ERI_Fxyz_Px_Dxz_Pz_bd = I_ERI_G2xyz_S_Dxz_Pz_bd+ABX*I_ERI_Fxyz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2z_Px_Dxz_Pz_bd = I_ERI_G2x2z_S_Dxz_Pz_bd+ABX*I_ERI_Fx2z_S_Dxz_Pz_bd;
  Double I_ERI_F3y_Px_Dxz_Pz_bd = I_ERI_Gx3y_S_Dxz_Pz_bd+ABX*I_ERI_F3y_S_Dxz_Pz_bd;
  Double I_ERI_F2yz_Px_Dxz_Pz_bd = I_ERI_Gx2yz_S_Dxz_Pz_bd+ABX*I_ERI_F2yz_S_Dxz_Pz_bd;
  Double I_ERI_Fy2z_Px_Dxz_Pz_bd = I_ERI_Gxy2z_S_Dxz_Pz_bd+ABX*I_ERI_Fy2z_S_Dxz_Pz_bd;
  Double I_ERI_F3z_Px_Dxz_Pz_bd = I_ERI_Gx3z_S_Dxz_Pz_bd+ABX*I_ERI_F3z_S_Dxz_Pz_bd;
  Double I_ERI_F3x_Py_Dxz_Pz_bd = I_ERI_G3xy_S_Dxz_Pz_bd+ABY*I_ERI_F3x_S_Dxz_Pz_bd;
  Double I_ERI_F2xy_Py_Dxz_Pz_bd = I_ERI_G2x2y_S_Dxz_Pz_bd+ABY*I_ERI_F2xy_S_Dxz_Pz_bd;
  Double I_ERI_F2xz_Py_Dxz_Pz_bd = I_ERI_G2xyz_S_Dxz_Pz_bd+ABY*I_ERI_F2xz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2y_Py_Dxz_Pz_bd = I_ERI_Gx3y_S_Dxz_Pz_bd+ABY*I_ERI_Fx2y_S_Dxz_Pz_bd;
  Double I_ERI_Fxyz_Py_Dxz_Pz_bd = I_ERI_Gx2yz_S_Dxz_Pz_bd+ABY*I_ERI_Fxyz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2z_Py_Dxz_Pz_bd = I_ERI_Gxy2z_S_Dxz_Pz_bd+ABY*I_ERI_Fx2z_S_Dxz_Pz_bd;
  Double I_ERI_F3y_Py_Dxz_Pz_bd = I_ERI_G4y_S_Dxz_Pz_bd+ABY*I_ERI_F3y_S_Dxz_Pz_bd;
  Double I_ERI_F2yz_Py_Dxz_Pz_bd = I_ERI_G3yz_S_Dxz_Pz_bd+ABY*I_ERI_F2yz_S_Dxz_Pz_bd;
  Double I_ERI_Fy2z_Py_Dxz_Pz_bd = I_ERI_G2y2z_S_Dxz_Pz_bd+ABY*I_ERI_Fy2z_S_Dxz_Pz_bd;
  Double I_ERI_F3z_Py_Dxz_Pz_bd = I_ERI_Gy3z_S_Dxz_Pz_bd+ABY*I_ERI_F3z_S_Dxz_Pz_bd;
  Double I_ERI_F3x_Pz_Dxz_Pz_bd = I_ERI_G3xz_S_Dxz_Pz_bd+ABZ*I_ERI_F3x_S_Dxz_Pz_bd;
  Double I_ERI_F2xy_Pz_Dxz_Pz_bd = I_ERI_G2xyz_S_Dxz_Pz_bd+ABZ*I_ERI_F2xy_S_Dxz_Pz_bd;
  Double I_ERI_F2xz_Pz_Dxz_Pz_bd = I_ERI_G2x2z_S_Dxz_Pz_bd+ABZ*I_ERI_F2xz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2y_Pz_Dxz_Pz_bd = I_ERI_Gx2yz_S_Dxz_Pz_bd+ABZ*I_ERI_Fx2y_S_Dxz_Pz_bd;
  Double I_ERI_Fxyz_Pz_Dxz_Pz_bd = I_ERI_Gxy2z_S_Dxz_Pz_bd+ABZ*I_ERI_Fxyz_S_Dxz_Pz_bd;
  Double I_ERI_Fx2z_Pz_Dxz_Pz_bd = I_ERI_Gx3z_S_Dxz_Pz_bd+ABZ*I_ERI_Fx2z_S_Dxz_Pz_bd;
  Double I_ERI_F3y_Pz_Dxz_Pz_bd = I_ERI_G3yz_S_Dxz_Pz_bd+ABZ*I_ERI_F3y_S_Dxz_Pz_bd;
  Double I_ERI_F2yz_Pz_Dxz_Pz_bd = I_ERI_G2y2z_S_Dxz_Pz_bd+ABZ*I_ERI_F2yz_S_Dxz_Pz_bd;
  Double I_ERI_Fy2z_Pz_Dxz_Pz_bd = I_ERI_Gy3z_S_Dxz_Pz_bd+ABZ*I_ERI_Fy2z_S_Dxz_Pz_bd;
  Double I_ERI_F3z_Pz_Dxz_Pz_bd = I_ERI_G4z_S_Dxz_Pz_bd+ABZ*I_ERI_F3z_S_Dxz_Pz_bd;
  Double I_ERI_F3x_Px_D2y_Pz_bd = I_ERI_G4x_S_D2y_Pz_bd+ABX*I_ERI_F3x_S_D2y_Pz_bd;
  Double I_ERI_F2xy_Px_D2y_Pz_bd = I_ERI_G3xy_S_D2y_Pz_bd+ABX*I_ERI_F2xy_S_D2y_Pz_bd;
  Double I_ERI_F2xz_Px_D2y_Pz_bd = I_ERI_G3xz_S_D2y_Pz_bd+ABX*I_ERI_F2xz_S_D2y_Pz_bd;
  Double I_ERI_Fx2y_Px_D2y_Pz_bd = I_ERI_G2x2y_S_D2y_Pz_bd+ABX*I_ERI_Fx2y_S_D2y_Pz_bd;
  Double I_ERI_Fxyz_Px_D2y_Pz_bd = I_ERI_G2xyz_S_D2y_Pz_bd+ABX*I_ERI_Fxyz_S_D2y_Pz_bd;
  Double I_ERI_Fx2z_Px_D2y_Pz_bd = I_ERI_G2x2z_S_D2y_Pz_bd+ABX*I_ERI_Fx2z_S_D2y_Pz_bd;
  Double I_ERI_F3y_Px_D2y_Pz_bd = I_ERI_Gx3y_S_D2y_Pz_bd+ABX*I_ERI_F3y_S_D2y_Pz_bd;
  Double I_ERI_F2yz_Px_D2y_Pz_bd = I_ERI_Gx2yz_S_D2y_Pz_bd+ABX*I_ERI_F2yz_S_D2y_Pz_bd;
  Double I_ERI_Fy2z_Px_D2y_Pz_bd = I_ERI_Gxy2z_S_D2y_Pz_bd+ABX*I_ERI_Fy2z_S_D2y_Pz_bd;
  Double I_ERI_F3z_Px_D2y_Pz_bd = I_ERI_Gx3z_S_D2y_Pz_bd+ABX*I_ERI_F3z_S_D2y_Pz_bd;
  Double I_ERI_F3x_Py_D2y_Pz_bd = I_ERI_G3xy_S_D2y_Pz_bd+ABY*I_ERI_F3x_S_D2y_Pz_bd;
  Double I_ERI_F2xy_Py_D2y_Pz_bd = I_ERI_G2x2y_S_D2y_Pz_bd+ABY*I_ERI_F2xy_S_D2y_Pz_bd;
  Double I_ERI_F2xz_Py_D2y_Pz_bd = I_ERI_G2xyz_S_D2y_Pz_bd+ABY*I_ERI_F2xz_S_D2y_Pz_bd;
  Double I_ERI_Fx2y_Py_D2y_Pz_bd = I_ERI_Gx3y_S_D2y_Pz_bd+ABY*I_ERI_Fx2y_S_D2y_Pz_bd;
  Double I_ERI_Fxyz_Py_D2y_Pz_bd = I_ERI_Gx2yz_S_D2y_Pz_bd+ABY*I_ERI_Fxyz_S_D2y_Pz_bd;
  Double I_ERI_Fx2z_Py_D2y_Pz_bd = I_ERI_Gxy2z_S_D2y_Pz_bd+ABY*I_ERI_Fx2z_S_D2y_Pz_bd;
  Double I_ERI_F3y_Py_D2y_Pz_bd = I_ERI_G4y_S_D2y_Pz_bd+ABY*I_ERI_F3y_S_D2y_Pz_bd;
  Double I_ERI_F2yz_Py_D2y_Pz_bd = I_ERI_G3yz_S_D2y_Pz_bd+ABY*I_ERI_F2yz_S_D2y_Pz_bd;
  Double I_ERI_Fy2z_Py_D2y_Pz_bd = I_ERI_G2y2z_S_D2y_Pz_bd+ABY*I_ERI_Fy2z_S_D2y_Pz_bd;
  Double I_ERI_F3z_Py_D2y_Pz_bd = I_ERI_Gy3z_S_D2y_Pz_bd+ABY*I_ERI_F3z_S_D2y_Pz_bd;
  Double I_ERI_F3x_Pz_D2y_Pz_bd = I_ERI_G3xz_S_D2y_Pz_bd+ABZ*I_ERI_F3x_S_D2y_Pz_bd;
  Double I_ERI_F2xy_Pz_D2y_Pz_bd = I_ERI_G2xyz_S_D2y_Pz_bd+ABZ*I_ERI_F2xy_S_D2y_Pz_bd;
  Double I_ERI_F2xz_Pz_D2y_Pz_bd = I_ERI_G2x2z_S_D2y_Pz_bd+ABZ*I_ERI_F2xz_S_D2y_Pz_bd;
  Double I_ERI_Fx2y_Pz_D2y_Pz_bd = I_ERI_Gx2yz_S_D2y_Pz_bd+ABZ*I_ERI_Fx2y_S_D2y_Pz_bd;
  Double I_ERI_Fxyz_Pz_D2y_Pz_bd = I_ERI_Gxy2z_S_D2y_Pz_bd+ABZ*I_ERI_Fxyz_S_D2y_Pz_bd;
  Double I_ERI_Fx2z_Pz_D2y_Pz_bd = I_ERI_Gx3z_S_D2y_Pz_bd+ABZ*I_ERI_Fx2z_S_D2y_Pz_bd;
  Double I_ERI_F3y_Pz_D2y_Pz_bd = I_ERI_G3yz_S_D2y_Pz_bd+ABZ*I_ERI_F3y_S_D2y_Pz_bd;
  Double I_ERI_F2yz_Pz_D2y_Pz_bd = I_ERI_G2y2z_S_D2y_Pz_bd+ABZ*I_ERI_F2yz_S_D2y_Pz_bd;
  Double I_ERI_Fy2z_Pz_D2y_Pz_bd = I_ERI_Gy3z_S_D2y_Pz_bd+ABZ*I_ERI_Fy2z_S_D2y_Pz_bd;
  Double I_ERI_F3z_Pz_D2y_Pz_bd = I_ERI_G4z_S_D2y_Pz_bd+ABZ*I_ERI_F3z_S_D2y_Pz_bd;
  Double I_ERI_F3x_Px_Dyz_Pz_bd = I_ERI_G4x_S_Dyz_Pz_bd+ABX*I_ERI_F3x_S_Dyz_Pz_bd;
  Double I_ERI_F2xy_Px_Dyz_Pz_bd = I_ERI_G3xy_S_Dyz_Pz_bd+ABX*I_ERI_F2xy_S_Dyz_Pz_bd;
  Double I_ERI_F2xz_Px_Dyz_Pz_bd = I_ERI_G3xz_S_Dyz_Pz_bd+ABX*I_ERI_F2xz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2y_Px_Dyz_Pz_bd = I_ERI_G2x2y_S_Dyz_Pz_bd+ABX*I_ERI_Fx2y_S_Dyz_Pz_bd;
  Double I_ERI_Fxyz_Px_Dyz_Pz_bd = I_ERI_G2xyz_S_Dyz_Pz_bd+ABX*I_ERI_Fxyz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2z_Px_Dyz_Pz_bd = I_ERI_G2x2z_S_Dyz_Pz_bd+ABX*I_ERI_Fx2z_S_Dyz_Pz_bd;
  Double I_ERI_F3y_Px_Dyz_Pz_bd = I_ERI_Gx3y_S_Dyz_Pz_bd+ABX*I_ERI_F3y_S_Dyz_Pz_bd;
  Double I_ERI_F2yz_Px_Dyz_Pz_bd = I_ERI_Gx2yz_S_Dyz_Pz_bd+ABX*I_ERI_F2yz_S_Dyz_Pz_bd;
  Double I_ERI_Fy2z_Px_Dyz_Pz_bd = I_ERI_Gxy2z_S_Dyz_Pz_bd+ABX*I_ERI_Fy2z_S_Dyz_Pz_bd;
  Double I_ERI_F3z_Px_Dyz_Pz_bd = I_ERI_Gx3z_S_Dyz_Pz_bd+ABX*I_ERI_F3z_S_Dyz_Pz_bd;
  Double I_ERI_F3x_Py_Dyz_Pz_bd = I_ERI_G3xy_S_Dyz_Pz_bd+ABY*I_ERI_F3x_S_Dyz_Pz_bd;
  Double I_ERI_F2xy_Py_Dyz_Pz_bd = I_ERI_G2x2y_S_Dyz_Pz_bd+ABY*I_ERI_F2xy_S_Dyz_Pz_bd;
  Double I_ERI_F2xz_Py_Dyz_Pz_bd = I_ERI_G2xyz_S_Dyz_Pz_bd+ABY*I_ERI_F2xz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2y_Py_Dyz_Pz_bd = I_ERI_Gx3y_S_Dyz_Pz_bd+ABY*I_ERI_Fx2y_S_Dyz_Pz_bd;
  Double I_ERI_Fxyz_Py_Dyz_Pz_bd = I_ERI_Gx2yz_S_Dyz_Pz_bd+ABY*I_ERI_Fxyz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2z_Py_Dyz_Pz_bd = I_ERI_Gxy2z_S_Dyz_Pz_bd+ABY*I_ERI_Fx2z_S_Dyz_Pz_bd;
  Double I_ERI_F3y_Py_Dyz_Pz_bd = I_ERI_G4y_S_Dyz_Pz_bd+ABY*I_ERI_F3y_S_Dyz_Pz_bd;
  Double I_ERI_F2yz_Py_Dyz_Pz_bd = I_ERI_G3yz_S_Dyz_Pz_bd+ABY*I_ERI_F2yz_S_Dyz_Pz_bd;
  Double I_ERI_Fy2z_Py_Dyz_Pz_bd = I_ERI_G2y2z_S_Dyz_Pz_bd+ABY*I_ERI_Fy2z_S_Dyz_Pz_bd;
  Double I_ERI_F3z_Py_Dyz_Pz_bd = I_ERI_Gy3z_S_Dyz_Pz_bd+ABY*I_ERI_F3z_S_Dyz_Pz_bd;
  Double I_ERI_F3x_Pz_Dyz_Pz_bd = I_ERI_G3xz_S_Dyz_Pz_bd+ABZ*I_ERI_F3x_S_Dyz_Pz_bd;
  Double I_ERI_F2xy_Pz_Dyz_Pz_bd = I_ERI_G2xyz_S_Dyz_Pz_bd+ABZ*I_ERI_F2xy_S_Dyz_Pz_bd;
  Double I_ERI_F2xz_Pz_Dyz_Pz_bd = I_ERI_G2x2z_S_Dyz_Pz_bd+ABZ*I_ERI_F2xz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2y_Pz_Dyz_Pz_bd = I_ERI_Gx2yz_S_Dyz_Pz_bd+ABZ*I_ERI_Fx2y_S_Dyz_Pz_bd;
  Double I_ERI_Fxyz_Pz_Dyz_Pz_bd = I_ERI_Gxy2z_S_Dyz_Pz_bd+ABZ*I_ERI_Fxyz_S_Dyz_Pz_bd;
  Double I_ERI_Fx2z_Pz_Dyz_Pz_bd = I_ERI_Gx3z_S_Dyz_Pz_bd+ABZ*I_ERI_Fx2z_S_Dyz_Pz_bd;
  Double I_ERI_F3y_Pz_Dyz_Pz_bd = I_ERI_G3yz_S_Dyz_Pz_bd+ABZ*I_ERI_F3y_S_Dyz_Pz_bd;
  Double I_ERI_F2yz_Pz_Dyz_Pz_bd = I_ERI_G2y2z_S_Dyz_Pz_bd+ABZ*I_ERI_F2yz_S_Dyz_Pz_bd;
  Double I_ERI_Fy2z_Pz_Dyz_Pz_bd = I_ERI_Gy3z_S_Dyz_Pz_bd+ABZ*I_ERI_Fy2z_S_Dyz_Pz_bd;
  Double I_ERI_F3z_Pz_Dyz_Pz_bd = I_ERI_G4z_S_Dyz_Pz_bd+ABZ*I_ERI_F3z_S_Dyz_Pz_bd;
  Double I_ERI_F3x_Px_D2z_Pz_bd = I_ERI_G4x_S_D2z_Pz_bd+ABX*I_ERI_F3x_S_D2z_Pz_bd;
  Double I_ERI_F2xy_Px_D2z_Pz_bd = I_ERI_G3xy_S_D2z_Pz_bd+ABX*I_ERI_F2xy_S_D2z_Pz_bd;
  Double I_ERI_F2xz_Px_D2z_Pz_bd = I_ERI_G3xz_S_D2z_Pz_bd+ABX*I_ERI_F2xz_S_D2z_Pz_bd;
  Double I_ERI_Fx2y_Px_D2z_Pz_bd = I_ERI_G2x2y_S_D2z_Pz_bd+ABX*I_ERI_Fx2y_S_D2z_Pz_bd;
  Double I_ERI_Fxyz_Px_D2z_Pz_bd = I_ERI_G2xyz_S_D2z_Pz_bd+ABX*I_ERI_Fxyz_S_D2z_Pz_bd;
  Double I_ERI_Fx2z_Px_D2z_Pz_bd = I_ERI_G2x2z_S_D2z_Pz_bd+ABX*I_ERI_Fx2z_S_D2z_Pz_bd;
  Double I_ERI_F3y_Px_D2z_Pz_bd = I_ERI_Gx3y_S_D2z_Pz_bd+ABX*I_ERI_F3y_S_D2z_Pz_bd;
  Double I_ERI_F2yz_Px_D2z_Pz_bd = I_ERI_Gx2yz_S_D2z_Pz_bd+ABX*I_ERI_F2yz_S_D2z_Pz_bd;
  Double I_ERI_Fy2z_Px_D2z_Pz_bd = I_ERI_Gxy2z_S_D2z_Pz_bd+ABX*I_ERI_Fy2z_S_D2z_Pz_bd;
  Double I_ERI_F3z_Px_D2z_Pz_bd = I_ERI_Gx3z_S_D2z_Pz_bd+ABX*I_ERI_F3z_S_D2z_Pz_bd;
  Double I_ERI_F3x_Py_D2z_Pz_bd = I_ERI_G3xy_S_D2z_Pz_bd+ABY*I_ERI_F3x_S_D2z_Pz_bd;
  Double I_ERI_F2xy_Py_D2z_Pz_bd = I_ERI_G2x2y_S_D2z_Pz_bd+ABY*I_ERI_F2xy_S_D2z_Pz_bd;
  Double I_ERI_F2xz_Py_D2z_Pz_bd = I_ERI_G2xyz_S_D2z_Pz_bd+ABY*I_ERI_F2xz_S_D2z_Pz_bd;
  Double I_ERI_Fx2y_Py_D2z_Pz_bd = I_ERI_Gx3y_S_D2z_Pz_bd+ABY*I_ERI_Fx2y_S_D2z_Pz_bd;
  Double I_ERI_Fxyz_Py_D2z_Pz_bd = I_ERI_Gx2yz_S_D2z_Pz_bd+ABY*I_ERI_Fxyz_S_D2z_Pz_bd;
  Double I_ERI_Fx2z_Py_D2z_Pz_bd = I_ERI_Gxy2z_S_D2z_Pz_bd+ABY*I_ERI_Fx2z_S_D2z_Pz_bd;
  Double I_ERI_F3y_Py_D2z_Pz_bd = I_ERI_G4y_S_D2z_Pz_bd+ABY*I_ERI_F3y_S_D2z_Pz_bd;
  Double I_ERI_F2yz_Py_D2z_Pz_bd = I_ERI_G3yz_S_D2z_Pz_bd+ABY*I_ERI_F2yz_S_D2z_Pz_bd;
  Double I_ERI_Fy2z_Py_D2z_Pz_bd = I_ERI_G2y2z_S_D2z_Pz_bd+ABY*I_ERI_Fy2z_S_D2z_Pz_bd;
  Double I_ERI_F3z_Py_D2z_Pz_bd = I_ERI_Gy3z_S_D2z_Pz_bd+ABY*I_ERI_F3z_S_D2z_Pz_bd;
  Double I_ERI_F3x_Pz_D2z_Pz_bd = I_ERI_G3xz_S_D2z_Pz_bd+ABZ*I_ERI_F3x_S_D2z_Pz_bd;
  Double I_ERI_F2xy_Pz_D2z_Pz_bd = I_ERI_G2xyz_S_D2z_Pz_bd+ABZ*I_ERI_F2xy_S_D2z_Pz_bd;
  Double I_ERI_F2xz_Pz_D2z_Pz_bd = I_ERI_G2x2z_S_D2z_Pz_bd+ABZ*I_ERI_F2xz_S_D2z_Pz_bd;
  Double I_ERI_Fx2y_Pz_D2z_Pz_bd = I_ERI_Gx2yz_S_D2z_Pz_bd+ABZ*I_ERI_Fx2y_S_D2z_Pz_bd;
  Double I_ERI_Fxyz_Pz_D2z_Pz_bd = I_ERI_Gxy2z_S_D2z_Pz_bd+ABZ*I_ERI_Fxyz_S_D2z_Pz_bd;
  Double I_ERI_Fx2z_Pz_D2z_Pz_bd = I_ERI_Gx3z_S_D2z_Pz_bd+ABZ*I_ERI_Fx2z_S_D2z_Pz_bd;
  Double I_ERI_F3y_Pz_D2z_Pz_bd = I_ERI_G3yz_S_D2z_Pz_bd+ABZ*I_ERI_F3y_S_D2z_Pz_bd;
  Double I_ERI_F2yz_Pz_D2z_Pz_bd = I_ERI_G2y2z_S_D2z_Pz_bd+ABZ*I_ERI_F2yz_S_D2z_Pz_bd;
  Double I_ERI_Fy2z_Pz_D2z_Pz_bd = I_ERI_Gy3z_S_D2z_Pz_bd+ABZ*I_ERI_Fy2z_S_D2z_Pz_bd;
  Double I_ERI_F3z_Pz_D2z_Pz_bd = I_ERI_G4z_S_D2z_Pz_bd+ABZ*I_ERI_F3z_S_D2z_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_F3x_D2x_D2x_S_bb-2.0E0*1*I_ERI_F3x_S_D2x_S_b;
  abcd[1] = 4.0E0*I_ERI_F2xy_D2x_D2x_S_bb-2.0E0*1*I_ERI_F2xy_S_D2x_S_b;
  abcd[2] = 4.0E0*I_ERI_F2xz_D2x_D2x_S_bb-2.0E0*1*I_ERI_F2xz_S_D2x_S_b;
  abcd[3] = 4.0E0*I_ERI_Fx2y_D2x_D2x_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2x_S_b;
  abcd[4] = 4.0E0*I_ERI_Fxyz_D2x_D2x_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2x_S_b;
  abcd[5] = 4.0E0*I_ERI_Fx2z_D2x_D2x_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2x_S_b;
  abcd[6] = 4.0E0*I_ERI_F3y_D2x_D2x_S_bb-2.0E0*1*I_ERI_F3y_S_D2x_S_b;
  abcd[7] = 4.0E0*I_ERI_F2yz_D2x_D2x_S_bb-2.0E0*1*I_ERI_F2yz_S_D2x_S_b;
  abcd[8] = 4.0E0*I_ERI_Fy2z_D2x_D2x_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2x_S_b;
  abcd[9] = 4.0E0*I_ERI_F3z_D2x_D2x_S_bb-2.0E0*1*I_ERI_F3z_S_D2x_S_b;
  abcd[10] = 4.0E0*I_ERI_F3x_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F3x_S_Dxy_S_b;
  abcd[11] = 4.0E0*I_ERI_F2xy_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxy_S_b;
  abcd[12] = 4.0E0*I_ERI_F2xz_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxy_S_b;
  abcd[13] = 4.0E0*I_ERI_Fx2y_D2x_Dxy_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_b;
  abcd[14] = 4.0E0*I_ERI_Fxyz_D2x_Dxy_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_b;
  abcd[15] = 4.0E0*I_ERI_Fx2z_D2x_Dxy_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_b;
  abcd[16] = 4.0E0*I_ERI_F3y_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F3y_S_Dxy_S_b;
  abcd[17] = 4.0E0*I_ERI_F2yz_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxy_S_b;
  abcd[18] = 4.0E0*I_ERI_Fy2z_D2x_Dxy_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_b;
  abcd[19] = 4.0E0*I_ERI_F3z_D2x_Dxy_S_bb-2.0E0*1*I_ERI_F3z_S_Dxy_S_b;
  abcd[20] = 4.0E0*I_ERI_F3x_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F3x_S_Dxz_S_b;
  abcd[21] = 4.0E0*I_ERI_F2xy_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxz_S_b;
  abcd[22] = 4.0E0*I_ERI_F2xz_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxz_S_b;
  abcd[23] = 4.0E0*I_ERI_Fx2y_D2x_Dxz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_b;
  abcd[24] = 4.0E0*I_ERI_Fxyz_D2x_Dxz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_b;
  abcd[25] = 4.0E0*I_ERI_Fx2z_D2x_Dxz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_b;
  abcd[26] = 4.0E0*I_ERI_F3y_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F3y_S_Dxz_S_b;
  abcd[27] = 4.0E0*I_ERI_F2yz_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxz_S_b;
  abcd[28] = 4.0E0*I_ERI_Fy2z_D2x_Dxz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_b;
  abcd[29] = 4.0E0*I_ERI_F3z_D2x_Dxz_S_bb-2.0E0*1*I_ERI_F3z_S_Dxz_S_b;
  abcd[30] = 4.0E0*I_ERI_F3x_D2x_D2y_S_bb-2.0E0*1*I_ERI_F3x_S_D2y_S_b;
  abcd[31] = 4.0E0*I_ERI_F2xy_D2x_D2y_S_bb-2.0E0*1*I_ERI_F2xy_S_D2y_S_b;
  abcd[32] = 4.0E0*I_ERI_F2xz_D2x_D2y_S_bb-2.0E0*1*I_ERI_F2xz_S_D2y_S_b;
  abcd[33] = 4.0E0*I_ERI_Fx2y_D2x_D2y_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2y_S_b;
  abcd[34] = 4.0E0*I_ERI_Fxyz_D2x_D2y_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2y_S_b;
  abcd[35] = 4.0E0*I_ERI_Fx2z_D2x_D2y_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2y_S_b;
  abcd[36] = 4.0E0*I_ERI_F3y_D2x_D2y_S_bb-2.0E0*1*I_ERI_F3y_S_D2y_S_b;
  abcd[37] = 4.0E0*I_ERI_F2yz_D2x_D2y_S_bb-2.0E0*1*I_ERI_F2yz_S_D2y_S_b;
  abcd[38] = 4.0E0*I_ERI_Fy2z_D2x_D2y_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2y_S_b;
  abcd[39] = 4.0E0*I_ERI_F3z_D2x_D2y_S_bb-2.0E0*1*I_ERI_F3z_S_D2y_S_b;
  abcd[40] = 4.0E0*I_ERI_F3x_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F3x_S_Dyz_S_b;
  abcd[41] = 4.0E0*I_ERI_F2xy_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dyz_S_b;
  abcd[42] = 4.0E0*I_ERI_F2xz_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dyz_S_b;
  abcd[43] = 4.0E0*I_ERI_Fx2y_D2x_Dyz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_b;
  abcd[44] = 4.0E0*I_ERI_Fxyz_D2x_Dyz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_b;
  abcd[45] = 4.0E0*I_ERI_Fx2z_D2x_Dyz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_b;
  abcd[46] = 4.0E0*I_ERI_F3y_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F3y_S_Dyz_S_b;
  abcd[47] = 4.0E0*I_ERI_F2yz_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dyz_S_b;
  abcd[48] = 4.0E0*I_ERI_Fy2z_D2x_Dyz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_b;
  abcd[49] = 4.0E0*I_ERI_F3z_D2x_Dyz_S_bb-2.0E0*1*I_ERI_F3z_S_Dyz_S_b;
  abcd[50] = 4.0E0*I_ERI_F3x_D2x_D2z_S_bb-2.0E0*1*I_ERI_F3x_S_D2z_S_b;
  abcd[51] = 4.0E0*I_ERI_F2xy_D2x_D2z_S_bb-2.0E0*1*I_ERI_F2xy_S_D2z_S_b;
  abcd[52] = 4.0E0*I_ERI_F2xz_D2x_D2z_S_bb-2.0E0*1*I_ERI_F2xz_S_D2z_S_b;
  abcd[53] = 4.0E0*I_ERI_Fx2y_D2x_D2z_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2z_S_b;
  abcd[54] = 4.0E0*I_ERI_Fxyz_D2x_D2z_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2z_S_b;
  abcd[55] = 4.0E0*I_ERI_Fx2z_D2x_D2z_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2z_S_b;
  abcd[56] = 4.0E0*I_ERI_F3y_D2x_D2z_S_bb-2.0E0*1*I_ERI_F3y_S_D2z_S_b;
  abcd[57] = 4.0E0*I_ERI_F2yz_D2x_D2z_S_bb-2.0E0*1*I_ERI_F2yz_S_D2z_S_b;
  abcd[58] = 4.0E0*I_ERI_Fy2z_D2x_D2z_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2z_S_b;
  abcd[59] = 4.0E0*I_ERI_F3z_D2x_D2z_S_bb-2.0E0*1*I_ERI_F3z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[60] = 4.0E0*I_ERI_F3x_Dxy_D2x_S_bb;
  abcd[61] = 4.0E0*I_ERI_F2xy_Dxy_D2x_S_bb;
  abcd[62] = 4.0E0*I_ERI_F2xz_Dxy_D2x_S_bb;
  abcd[63] = 4.0E0*I_ERI_Fx2y_Dxy_D2x_S_bb;
  abcd[64] = 4.0E0*I_ERI_Fxyz_Dxy_D2x_S_bb;
  abcd[65] = 4.0E0*I_ERI_Fx2z_Dxy_D2x_S_bb;
  abcd[66] = 4.0E0*I_ERI_F3y_Dxy_D2x_S_bb;
  abcd[67] = 4.0E0*I_ERI_F2yz_Dxy_D2x_S_bb;
  abcd[68] = 4.0E0*I_ERI_Fy2z_Dxy_D2x_S_bb;
  abcd[69] = 4.0E0*I_ERI_F3z_Dxy_D2x_S_bb;
  abcd[70] = 4.0E0*I_ERI_F3x_Dxy_Dxy_S_bb;
  abcd[71] = 4.0E0*I_ERI_F2xy_Dxy_Dxy_S_bb;
  abcd[72] = 4.0E0*I_ERI_F2xz_Dxy_Dxy_S_bb;
  abcd[73] = 4.0E0*I_ERI_Fx2y_Dxy_Dxy_S_bb;
  abcd[74] = 4.0E0*I_ERI_Fxyz_Dxy_Dxy_S_bb;
  abcd[75] = 4.0E0*I_ERI_Fx2z_Dxy_Dxy_S_bb;
  abcd[76] = 4.0E0*I_ERI_F3y_Dxy_Dxy_S_bb;
  abcd[77] = 4.0E0*I_ERI_F2yz_Dxy_Dxy_S_bb;
  abcd[78] = 4.0E0*I_ERI_Fy2z_Dxy_Dxy_S_bb;
  abcd[79] = 4.0E0*I_ERI_F3z_Dxy_Dxy_S_bb;
  abcd[80] = 4.0E0*I_ERI_F3x_Dxy_Dxz_S_bb;
  abcd[81] = 4.0E0*I_ERI_F2xy_Dxy_Dxz_S_bb;
  abcd[82] = 4.0E0*I_ERI_F2xz_Dxy_Dxz_S_bb;
  abcd[83] = 4.0E0*I_ERI_Fx2y_Dxy_Dxz_S_bb;
  abcd[84] = 4.0E0*I_ERI_Fxyz_Dxy_Dxz_S_bb;
  abcd[85] = 4.0E0*I_ERI_Fx2z_Dxy_Dxz_S_bb;
  abcd[86] = 4.0E0*I_ERI_F3y_Dxy_Dxz_S_bb;
  abcd[87] = 4.0E0*I_ERI_F2yz_Dxy_Dxz_S_bb;
  abcd[88] = 4.0E0*I_ERI_Fy2z_Dxy_Dxz_S_bb;
  abcd[89] = 4.0E0*I_ERI_F3z_Dxy_Dxz_S_bb;
  abcd[90] = 4.0E0*I_ERI_F3x_Dxy_D2y_S_bb;
  abcd[91] = 4.0E0*I_ERI_F2xy_Dxy_D2y_S_bb;
  abcd[92] = 4.0E0*I_ERI_F2xz_Dxy_D2y_S_bb;
  abcd[93] = 4.0E0*I_ERI_Fx2y_Dxy_D2y_S_bb;
  abcd[94] = 4.0E0*I_ERI_Fxyz_Dxy_D2y_S_bb;
  abcd[95] = 4.0E0*I_ERI_Fx2z_Dxy_D2y_S_bb;
  abcd[96] = 4.0E0*I_ERI_F3y_Dxy_D2y_S_bb;
  abcd[97] = 4.0E0*I_ERI_F2yz_Dxy_D2y_S_bb;
  abcd[98] = 4.0E0*I_ERI_Fy2z_Dxy_D2y_S_bb;
  abcd[99] = 4.0E0*I_ERI_F3z_Dxy_D2y_S_bb;
  abcd[100] = 4.0E0*I_ERI_F3x_Dxy_Dyz_S_bb;
  abcd[101] = 4.0E0*I_ERI_F2xy_Dxy_Dyz_S_bb;
  abcd[102] = 4.0E0*I_ERI_F2xz_Dxy_Dyz_S_bb;
  abcd[103] = 4.0E0*I_ERI_Fx2y_Dxy_Dyz_S_bb;
  abcd[104] = 4.0E0*I_ERI_Fxyz_Dxy_Dyz_S_bb;
  abcd[105] = 4.0E0*I_ERI_Fx2z_Dxy_Dyz_S_bb;
  abcd[106] = 4.0E0*I_ERI_F3y_Dxy_Dyz_S_bb;
  abcd[107] = 4.0E0*I_ERI_F2yz_Dxy_Dyz_S_bb;
  abcd[108] = 4.0E0*I_ERI_Fy2z_Dxy_Dyz_S_bb;
  abcd[109] = 4.0E0*I_ERI_F3z_Dxy_Dyz_S_bb;
  abcd[110] = 4.0E0*I_ERI_F3x_Dxy_D2z_S_bb;
  abcd[111] = 4.0E0*I_ERI_F2xy_Dxy_D2z_S_bb;
  abcd[112] = 4.0E0*I_ERI_F2xz_Dxy_D2z_S_bb;
  abcd[113] = 4.0E0*I_ERI_Fx2y_Dxy_D2z_S_bb;
  abcd[114] = 4.0E0*I_ERI_Fxyz_Dxy_D2z_S_bb;
  abcd[115] = 4.0E0*I_ERI_Fx2z_Dxy_D2z_S_bb;
  abcd[116] = 4.0E0*I_ERI_F3y_Dxy_D2z_S_bb;
  abcd[117] = 4.0E0*I_ERI_F2yz_Dxy_D2z_S_bb;
  abcd[118] = 4.0E0*I_ERI_Fy2z_Dxy_D2z_S_bb;
  abcd[119] = 4.0E0*I_ERI_F3z_Dxy_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[120] = 4.0E0*I_ERI_F3x_Dxz_D2x_S_bb;
  abcd[121] = 4.0E0*I_ERI_F2xy_Dxz_D2x_S_bb;
  abcd[122] = 4.0E0*I_ERI_F2xz_Dxz_D2x_S_bb;
  abcd[123] = 4.0E0*I_ERI_Fx2y_Dxz_D2x_S_bb;
  abcd[124] = 4.0E0*I_ERI_Fxyz_Dxz_D2x_S_bb;
  abcd[125] = 4.0E0*I_ERI_Fx2z_Dxz_D2x_S_bb;
  abcd[126] = 4.0E0*I_ERI_F3y_Dxz_D2x_S_bb;
  abcd[127] = 4.0E0*I_ERI_F2yz_Dxz_D2x_S_bb;
  abcd[128] = 4.0E0*I_ERI_Fy2z_Dxz_D2x_S_bb;
  abcd[129] = 4.0E0*I_ERI_F3z_Dxz_D2x_S_bb;
  abcd[130] = 4.0E0*I_ERI_F3x_Dxz_Dxy_S_bb;
  abcd[131] = 4.0E0*I_ERI_F2xy_Dxz_Dxy_S_bb;
  abcd[132] = 4.0E0*I_ERI_F2xz_Dxz_Dxy_S_bb;
  abcd[133] = 4.0E0*I_ERI_Fx2y_Dxz_Dxy_S_bb;
  abcd[134] = 4.0E0*I_ERI_Fxyz_Dxz_Dxy_S_bb;
  abcd[135] = 4.0E0*I_ERI_Fx2z_Dxz_Dxy_S_bb;
  abcd[136] = 4.0E0*I_ERI_F3y_Dxz_Dxy_S_bb;
  abcd[137] = 4.0E0*I_ERI_F2yz_Dxz_Dxy_S_bb;
  abcd[138] = 4.0E0*I_ERI_Fy2z_Dxz_Dxy_S_bb;
  abcd[139] = 4.0E0*I_ERI_F3z_Dxz_Dxy_S_bb;
  abcd[140] = 4.0E0*I_ERI_F3x_Dxz_Dxz_S_bb;
  abcd[141] = 4.0E0*I_ERI_F2xy_Dxz_Dxz_S_bb;
  abcd[142] = 4.0E0*I_ERI_F2xz_Dxz_Dxz_S_bb;
  abcd[143] = 4.0E0*I_ERI_Fx2y_Dxz_Dxz_S_bb;
  abcd[144] = 4.0E0*I_ERI_Fxyz_Dxz_Dxz_S_bb;
  abcd[145] = 4.0E0*I_ERI_Fx2z_Dxz_Dxz_S_bb;
  abcd[146] = 4.0E0*I_ERI_F3y_Dxz_Dxz_S_bb;
  abcd[147] = 4.0E0*I_ERI_F2yz_Dxz_Dxz_S_bb;
  abcd[148] = 4.0E0*I_ERI_Fy2z_Dxz_Dxz_S_bb;
  abcd[149] = 4.0E0*I_ERI_F3z_Dxz_Dxz_S_bb;
  abcd[150] = 4.0E0*I_ERI_F3x_Dxz_D2y_S_bb;
  abcd[151] = 4.0E0*I_ERI_F2xy_Dxz_D2y_S_bb;
  abcd[152] = 4.0E0*I_ERI_F2xz_Dxz_D2y_S_bb;
  abcd[153] = 4.0E0*I_ERI_Fx2y_Dxz_D2y_S_bb;
  abcd[154] = 4.0E0*I_ERI_Fxyz_Dxz_D2y_S_bb;
  abcd[155] = 4.0E0*I_ERI_Fx2z_Dxz_D2y_S_bb;
  abcd[156] = 4.0E0*I_ERI_F3y_Dxz_D2y_S_bb;
  abcd[157] = 4.0E0*I_ERI_F2yz_Dxz_D2y_S_bb;
  abcd[158] = 4.0E0*I_ERI_Fy2z_Dxz_D2y_S_bb;
  abcd[159] = 4.0E0*I_ERI_F3z_Dxz_D2y_S_bb;
  abcd[160] = 4.0E0*I_ERI_F3x_Dxz_Dyz_S_bb;
  abcd[161] = 4.0E0*I_ERI_F2xy_Dxz_Dyz_S_bb;
  abcd[162] = 4.0E0*I_ERI_F2xz_Dxz_Dyz_S_bb;
  abcd[163] = 4.0E0*I_ERI_Fx2y_Dxz_Dyz_S_bb;
  abcd[164] = 4.0E0*I_ERI_Fxyz_Dxz_Dyz_S_bb;
  abcd[165] = 4.0E0*I_ERI_Fx2z_Dxz_Dyz_S_bb;
  abcd[166] = 4.0E0*I_ERI_F3y_Dxz_Dyz_S_bb;
  abcd[167] = 4.0E0*I_ERI_F2yz_Dxz_Dyz_S_bb;
  abcd[168] = 4.0E0*I_ERI_Fy2z_Dxz_Dyz_S_bb;
  abcd[169] = 4.0E0*I_ERI_F3z_Dxz_Dyz_S_bb;
  abcd[170] = 4.0E0*I_ERI_F3x_Dxz_D2z_S_bb;
  abcd[171] = 4.0E0*I_ERI_F2xy_Dxz_D2z_S_bb;
  abcd[172] = 4.0E0*I_ERI_F2xz_Dxz_D2z_S_bb;
  abcd[173] = 4.0E0*I_ERI_Fx2y_Dxz_D2z_S_bb;
  abcd[174] = 4.0E0*I_ERI_Fxyz_Dxz_D2z_S_bb;
  abcd[175] = 4.0E0*I_ERI_Fx2z_Dxz_D2z_S_bb;
  abcd[176] = 4.0E0*I_ERI_F3y_Dxz_D2z_S_bb;
  abcd[177] = 4.0E0*I_ERI_F2yz_Dxz_D2z_S_bb;
  abcd[178] = 4.0E0*I_ERI_Fy2z_Dxz_D2z_S_bb;
  abcd[179] = 4.0E0*I_ERI_F3z_Dxz_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_F3x_D2y_D2x_S_bb-2.0E0*1*I_ERI_F3x_S_D2x_S_b;
  abcd[181] = 4.0E0*I_ERI_F2xy_D2y_D2x_S_bb-2.0E0*1*I_ERI_F2xy_S_D2x_S_b;
  abcd[182] = 4.0E0*I_ERI_F2xz_D2y_D2x_S_bb-2.0E0*1*I_ERI_F2xz_S_D2x_S_b;
  abcd[183] = 4.0E0*I_ERI_Fx2y_D2y_D2x_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2x_S_b;
  abcd[184] = 4.0E0*I_ERI_Fxyz_D2y_D2x_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2x_S_b;
  abcd[185] = 4.0E0*I_ERI_Fx2z_D2y_D2x_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2x_S_b;
  abcd[186] = 4.0E0*I_ERI_F3y_D2y_D2x_S_bb-2.0E0*1*I_ERI_F3y_S_D2x_S_b;
  abcd[187] = 4.0E0*I_ERI_F2yz_D2y_D2x_S_bb-2.0E0*1*I_ERI_F2yz_S_D2x_S_b;
  abcd[188] = 4.0E0*I_ERI_Fy2z_D2y_D2x_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2x_S_b;
  abcd[189] = 4.0E0*I_ERI_F3z_D2y_D2x_S_bb-2.0E0*1*I_ERI_F3z_S_D2x_S_b;
  abcd[190] = 4.0E0*I_ERI_F3x_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F3x_S_Dxy_S_b;
  abcd[191] = 4.0E0*I_ERI_F2xy_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxy_S_b;
  abcd[192] = 4.0E0*I_ERI_F2xz_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxy_S_b;
  abcd[193] = 4.0E0*I_ERI_Fx2y_D2y_Dxy_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_b;
  abcd[194] = 4.0E0*I_ERI_Fxyz_D2y_Dxy_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_b;
  abcd[195] = 4.0E0*I_ERI_Fx2z_D2y_Dxy_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_b;
  abcd[196] = 4.0E0*I_ERI_F3y_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F3y_S_Dxy_S_b;
  abcd[197] = 4.0E0*I_ERI_F2yz_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxy_S_b;
  abcd[198] = 4.0E0*I_ERI_Fy2z_D2y_Dxy_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_b;
  abcd[199] = 4.0E0*I_ERI_F3z_D2y_Dxy_S_bb-2.0E0*1*I_ERI_F3z_S_Dxy_S_b;
  abcd[200] = 4.0E0*I_ERI_F3x_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F3x_S_Dxz_S_b;
  abcd[201] = 4.0E0*I_ERI_F2xy_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxz_S_b;
  abcd[202] = 4.0E0*I_ERI_F2xz_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxz_S_b;
  abcd[203] = 4.0E0*I_ERI_Fx2y_D2y_Dxz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_b;
  abcd[204] = 4.0E0*I_ERI_Fxyz_D2y_Dxz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_b;
  abcd[205] = 4.0E0*I_ERI_Fx2z_D2y_Dxz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_b;
  abcd[206] = 4.0E0*I_ERI_F3y_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F3y_S_Dxz_S_b;
  abcd[207] = 4.0E0*I_ERI_F2yz_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxz_S_b;
  abcd[208] = 4.0E0*I_ERI_Fy2z_D2y_Dxz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_b;
  abcd[209] = 4.0E0*I_ERI_F3z_D2y_Dxz_S_bb-2.0E0*1*I_ERI_F3z_S_Dxz_S_b;
  abcd[210] = 4.0E0*I_ERI_F3x_D2y_D2y_S_bb-2.0E0*1*I_ERI_F3x_S_D2y_S_b;
  abcd[211] = 4.0E0*I_ERI_F2xy_D2y_D2y_S_bb-2.0E0*1*I_ERI_F2xy_S_D2y_S_b;
  abcd[212] = 4.0E0*I_ERI_F2xz_D2y_D2y_S_bb-2.0E0*1*I_ERI_F2xz_S_D2y_S_b;
  abcd[213] = 4.0E0*I_ERI_Fx2y_D2y_D2y_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2y_S_b;
  abcd[214] = 4.0E0*I_ERI_Fxyz_D2y_D2y_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2y_S_b;
  abcd[215] = 4.0E0*I_ERI_Fx2z_D2y_D2y_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2y_S_b;
  abcd[216] = 4.0E0*I_ERI_F3y_D2y_D2y_S_bb-2.0E0*1*I_ERI_F3y_S_D2y_S_b;
  abcd[217] = 4.0E0*I_ERI_F2yz_D2y_D2y_S_bb-2.0E0*1*I_ERI_F2yz_S_D2y_S_b;
  abcd[218] = 4.0E0*I_ERI_Fy2z_D2y_D2y_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2y_S_b;
  abcd[219] = 4.0E0*I_ERI_F3z_D2y_D2y_S_bb-2.0E0*1*I_ERI_F3z_S_D2y_S_b;
  abcd[220] = 4.0E0*I_ERI_F3x_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F3x_S_Dyz_S_b;
  abcd[221] = 4.0E0*I_ERI_F2xy_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dyz_S_b;
  abcd[222] = 4.0E0*I_ERI_F2xz_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dyz_S_b;
  abcd[223] = 4.0E0*I_ERI_Fx2y_D2y_Dyz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_b;
  abcd[224] = 4.0E0*I_ERI_Fxyz_D2y_Dyz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_b;
  abcd[225] = 4.0E0*I_ERI_Fx2z_D2y_Dyz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_b;
  abcd[226] = 4.0E0*I_ERI_F3y_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F3y_S_Dyz_S_b;
  abcd[227] = 4.0E0*I_ERI_F2yz_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dyz_S_b;
  abcd[228] = 4.0E0*I_ERI_Fy2z_D2y_Dyz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_b;
  abcd[229] = 4.0E0*I_ERI_F3z_D2y_Dyz_S_bb-2.0E0*1*I_ERI_F3z_S_Dyz_S_b;
  abcd[230] = 4.0E0*I_ERI_F3x_D2y_D2z_S_bb-2.0E0*1*I_ERI_F3x_S_D2z_S_b;
  abcd[231] = 4.0E0*I_ERI_F2xy_D2y_D2z_S_bb-2.0E0*1*I_ERI_F2xy_S_D2z_S_b;
  abcd[232] = 4.0E0*I_ERI_F2xz_D2y_D2z_S_bb-2.0E0*1*I_ERI_F2xz_S_D2z_S_b;
  abcd[233] = 4.0E0*I_ERI_Fx2y_D2y_D2z_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2z_S_b;
  abcd[234] = 4.0E0*I_ERI_Fxyz_D2y_D2z_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2z_S_b;
  abcd[235] = 4.0E0*I_ERI_Fx2z_D2y_D2z_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2z_S_b;
  abcd[236] = 4.0E0*I_ERI_F3y_D2y_D2z_S_bb-2.0E0*1*I_ERI_F3y_S_D2z_S_b;
  abcd[237] = 4.0E0*I_ERI_F2yz_D2y_D2z_S_bb-2.0E0*1*I_ERI_F2yz_S_D2z_S_b;
  abcd[238] = 4.0E0*I_ERI_Fy2z_D2y_D2z_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2z_S_b;
  abcd[239] = 4.0E0*I_ERI_F3z_D2y_D2z_S_bb-2.0E0*1*I_ERI_F3z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_F3x_Dyz_D2x_S_bb;
  abcd[241] = 4.0E0*I_ERI_F2xy_Dyz_D2x_S_bb;
  abcd[242] = 4.0E0*I_ERI_F2xz_Dyz_D2x_S_bb;
  abcd[243] = 4.0E0*I_ERI_Fx2y_Dyz_D2x_S_bb;
  abcd[244] = 4.0E0*I_ERI_Fxyz_Dyz_D2x_S_bb;
  abcd[245] = 4.0E0*I_ERI_Fx2z_Dyz_D2x_S_bb;
  abcd[246] = 4.0E0*I_ERI_F3y_Dyz_D2x_S_bb;
  abcd[247] = 4.0E0*I_ERI_F2yz_Dyz_D2x_S_bb;
  abcd[248] = 4.0E0*I_ERI_Fy2z_Dyz_D2x_S_bb;
  abcd[249] = 4.0E0*I_ERI_F3z_Dyz_D2x_S_bb;
  abcd[250] = 4.0E0*I_ERI_F3x_Dyz_Dxy_S_bb;
  abcd[251] = 4.0E0*I_ERI_F2xy_Dyz_Dxy_S_bb;
  abcd[252] = 4.0E0*I_ERI_F2xz_Dyz_Dxy_S_bb;
  abcd[253] = 4.0E0*I_ERI_Fx2y_Dyz_Dxy_S_bb;
  abcd[254] = 4.0E0*I_ERI_Fxyz_Dyz_Dxy_S_bb;
  abcd[255] = 4.0E0*I_ERI_Fx2z_Dyz_Dxy_S_bb;
  abcd[256] = 4.0E0*I_ERI_F3y_Dyz_Dxy_S_bb;
  abcd[257] = 4.0E0*I_ERI_F2yz_Dyz_Dxy_S_bb;
  abcd[258] = 4.0E0*I_ERI_Fy2z_Dyz_Dxy_S_bb;
  abcd[259] = 4.0E0*I_ERI_F3z_Dyz_Dxy_S_bb;
  abcd[260] = 4.0E0*I_ERI_F3x_Dyz_Dxz_S_bb;
  abcd[261] = 4.0E0*I_ERI_F2xy_Dyz_Dxz_S_bb;
  abcd[262] = 4.0E0*I_ERI_F2xz_Dyz_Dxz_S_bb;
  abcd[263] = 4.0E0*I_ERI_Fx2y_Dyz_Dxz_S_bb;
  abcd[264] = 4.0E0*I_ERI_Fxyz_Dyz_Dxz_S_bb;
  abcd[265] = 4.0E0*I_ERI_Fx2z_Dyz_Dxz_S_bb;
  abcd[266] = 4.0E0*I_ERI_F3y_Dyz_Dxz_S_bb;
  abcd[267] = 4.0E0*I_ERI_F2yz_Dyz_Dxz_S_bb;
  abcd[268] = 4.0E0*I_ERI_Fy2z_Dyz_Dxz_S_bb;
  abcd[269] = 4.0E0*I_ERI_F3z_Dyz_Dxz_S_bb;
  abcd[270] = 4.0E0*I_ERI_F3x_Dyz_D2y_S_bb;
  abcd[271] = 4.0E0*I_ERI_F2xy_Dyz_D2y_S_bb;
  abcd[272] = 4.0E0*I_ERI_F2xz_Dyz_D2y_S_bb;
  abcd[273] = 4.0E0*I_ERI_Fx2y_Dyz_D2y_S_bb;
  abcd[274] = 4.0E0*I_ERI_Fxyz_Dyz_D2y_S_bb;
  abcd[275] = 4.0E0*I_ERI_Fx2z_Dyz_D2y_S_bb;
  abcd[276] = 4.0E0*I_ERI_F3y_Dyz_D2y_S_bb;
  abcd[277] = 4.0E0*I_ERI_F2yz_Dyz_D2y_S_bb;
  abcd[278] = 4.0E0*I_ERI_Fy2z_Dyz_D2y_S_bb;
  abcd[279] = 4.0E0*I_ERI_F3z_Dyz_D2y_S_bb;
  abcd[280] = 4.0E0*I_ERI_F3x_Dyz_Dyz_S_bb;
  abcd[281] = 4.0E0*I_ERI_F2xy_Dyz_Dyz_S_bb;
  abcd[282] = 4.0E0*I_ERI_F2xz_Dyz_Dyz_S_bb;
  abcd[283] = 4.0E0*I_ERI_Fx2y_Dyz_Dyz_S_bb;
  abcd[284] = 4.0E0*I_ERI_Fxyz_Dyz_Dyz_S_bb;
  abcd[285] = 4.0E0*I_ERI_Fx2z_Dyz_Dyz_S_bb;
  abcd[286] = 4.0E0*I_ERI_F3y_Dyz_Dyz_S_bb;
  abcd[287] = 4.0E0*I_ERI_F2yz_Dyz_Dyz_S_bb;
  abcd[288] = 4.0E0*I_ERI_Fy2z_Dyz_Dyz_S_bb;
  abcd[289] = 4.0E0*I_ERI_F3z_Dyz_Dyz_S_bb;
  abcd[290] = 4.0E0*I_ERI_F3x_Dyz_D2z_S_bb;
  abcd[291] = 4.0E0*I_ERI_F2xy_Dyz_D2z_S_bb;
  abcd[292] = 4.0E0*I_ERI_F2xz_Dyz_D2z_S_bb;
  abcd[293] = 4.0E0*I_ERI_Fx2y_Dyz_D2z_S_bb;
  abcd[294] = 4.0E0*I_ERI_Fxyz_Dyz_D2z_S_bb;
  abcd[295] = 4.0E0*I_ERI_Fx2z_Dyz_D2z_S_bb;
  abcd[296] = 4.0E0*I_ERI_F3y_Dyz_D2z_S_bb;
  abcd[297] = 4.0E0*I_ERI_F2yz_Dyz_D2z_S_bb;
  abcd[298] = 4.0E0*I_ERI_Fy2z_Dyz_D2z_S_bb;
  abcd[299] = 4.0E0*I_ERI_F3z_Dyz_D2z_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  abcd[300] = 4.0E0*I_ERI_F3x_D2z_D2x_S_bb-2.0E0*1*I_ERI_F3x_S_D2x_S_b;
  abcd[301] = 4.0E0*I_ERI_F2xy_D2z_D2x_S_bb-2.0E0*1*I_ERI_F2xy_S_D2x_S_b;
  abcd[302] = 4.0E0*I_ERI_F2xz_D2z_D2x_S_bb-2.0E0*1*I_ERI_F2xz_S_D2x_S_b;
  abcd[303] = 4.0E0*I_ERI_Fx2y_D2z_D2x_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2x_S_b;
  abcd[304] = 4.0E0*I_ERI_Fxyz_D2z_D2x_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2x_S_b;
  abcd[305] = 4.0E0*I_ERI_Fx2z_D2z_D2x_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2x_S_b;
  abcd[306] = 4.0E0*I_ERI_F3y_D2z_D2x_S_bb-2.0E0*1*I_ERI_F3y_S_D2x_S_b;
  abcd[307] = 4.0E0*I_ERI_F2yz_D2z_D2x_S_bb-2.0E0*1*I_ERI_F2yz_S_D2x_S_b;
  abcd[308] = 4.0E0*I_ERI_Fy2z_D2z_D2x_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2x_S_b;
  abcd[309] = 4.0E0*I_ERI_F3z_D2z_D2x_S_bb-2.0E0*1*I_ERI_F3z_S_D2x_S_b;
  abcd[310] = 4.0E0*I_ERI_F3x_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F3x_S_Dxy_S_b;
  abcd[311] = 4.0E0*I_ERI_F2xy_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxy_S_b;
  abcd[312] = 4.0E0*I_ERI_F2xz_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxy_S_b;
  abcd[313] = 4.0E0*I_ERI_Fx2y_D2z_Dxy_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_b;
  abcd[314] = 4.0E0*I_ERI_Fxyz_D2z_Dxy_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_b;
  abcd[315] = 4.0E0*I_ERI_Fx2z_D2z_Dxy_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_b;
  abcd[316] = 4.0E0*I_ERI_F3y_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F3y_S_Dxy_S_b;
  abcd[317] = 4.0E0*I_ERI_F2yz_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxy_S_b;
  abcd[318] = 4.0E0*I_ERI_Fy2z_D2z_Dxy_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_b;
  abcd[319] = 4.0E0*I_ERI_F3z_D2z_Dxy_S_bb-2.0E0*1*I_ERI_F3z_S_Dxy_S_b;
  abcd[320] = 4.0E0*I_ERI_F3x_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F3x_S_Dxz_S_b;
  abcd[321] = 4.0E0*I_ERI_F2xy_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dxz_S_b;
  abcd[322] = 4.0E0*I_ERI_F2xz_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dxz_S_b;
  abcd[323] = 4.0E0*I_ERI_Fx2y_D2z_Dxz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_b;
  abcd[324] = 4.0E0*I_ERI_Fxyz_D2z_Dxz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_b;
  abcd[325] = 4.0E0*I_ERI_Fx2z_D2z_Dxz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_b;
  abcd[326] = 4.0E0*I_ERI_F3y_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F3y_S_Dxz_S_b;
  abcd[327] = 4.0E0*I_ERI_F2yz_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dxz_S_b;
  abcd[328] = 4.0E0*I_ERI_Fy2z_D2z_Dxz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_b;
  abcd[329] = 4.0E0*I_ERI_F3z_D2z_Dxz_S_bb-2.0E0*1*I_ERI_F3z_S_Dxz_S_b;
  abcd[330] = 4.0E0*I_ERI_F3x_D2z_D2y_S_bb-2.0E0*1*I_ERI_F3x_S_D2y_S_b;
  abcd[331] = 4.0E0*I_ERI_F2xy_D2z_D2y_S_bb-2.0E0*1*I_ERI_F2xy_S_D2y_S_b;
  abcd[332] = 4.0E0*I_ERI_F2xz_D2z_D2y_S_bb-2.0E0*1*I_ERI_F2xz_S_D2y_S_b;
  abcd[333] = 4.0E0*I_ERI_Fx2y_D2z_D2y_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2y_S_b;
  abcd[334] = 4.0E0*I_ERI_Fxyz_D2z_D2y_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2y_S_b;
  abcd[335] = 4.0E0*I_ERI_Fx2z_D2z_D2y_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2y_S_b;
  abcd[336] = 4.0E0*I_ERI_F3y_D2z_D2y_S_bb-2.0E0*1*I_ERI_F3y_S_D2y_S_b;
  abcd[337] = 4.0E0*I_ERI_F2yz_D2z_D2y_S_bb-2.0E0*1*I_ERI_F2yz_S_D2y_S_b;
  abcd[338] = 4.0E0*I_ERI_Fy2z_D2z_D2y_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2y_S_b;
  abcd[339] = 4.0E0*I_ERI_F3z_D2z_D2y_S_bb-2.0E0*1*I_ERI_F3z_S_D2y_S_b;
  abcd[340] = 4.0E0*I_ERI_F3x_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F3x_S_Dyz_S_b;
  abcd[341] = 4.0E0*I_ERI_F2xy_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F2xy_S_Dyz_S_b;
  abcd[342] = 4.0E0*I_ERI_F2xz_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F2xz_S_Dyz_S_b;
  abcd[343] = 4.0E0*I_ERI_Fx2y_D2z_Dyz_S_bb-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_b;
  abcd[344] = 4.0E0*I_ERI_Fxyz_D2z_Dyz_S_bb-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_b;
  abcd[345] = 4.0E0*I_ERI_Fx2z_D2z_Dyz_S_bb-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_b;
  abcd[346] = 4.0E0*I_ERI_F3y_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F3y_S_Dyz_S_b;
  abcd[347] = 4.0E0*I_ERI_F2yz_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F2yz_S_Dyz_S_b;
  abcd[348] = 4.0E0*I_ERI_Fy2z_D2z_Dyz_S_bb-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_b;
  abcd[349] = 4.0E0*I_ERI_F3z_D2z_Dyz_S_bb-2.0E0*1*I_ERI_F3z_S_Dyz_S_b;
  abcd[350] = 4.0E0*I_ERI_F3x_D2z_D2z_S_bb-2.0E0*1*I_ERI_F3x_S_D2z_S_b;
  abcd[351] = 4.0E0*I_ERI_F2xy_D2z_D2z_S_bb-2.0E0*1*I_ERI_F2xy_S_D2z_S_b;
  abcd[352] = 4.0E0*I_ERI_F2xz_D2z_D2z_S_bb-2.0E0*1*I_ERI_F2xz_S_D2z_S_b;
  abcd[353] = 4.0E0*I_ERI_Fx2y_D2z_D2z_S_bb-2.0E0*1*I_ERI_Fx2y_S_D2z_S_b;
  abcd[354] = 4.0E0*I_ERI_Fxyz_D2z_D2z_S_bb-2.0E0*1*I_ERI_Fxyz_S_D2z_S_b;
  abcd[355] = 4.0E0*I_ERI_Fx2z_D2z_D2z_S_bb-2.0E0*1*I_ERI_Fx2z_S_D2z_S_b;
  abcd[356] = 4.0E0*I_ERI_F3y_D2z_D2z_S_bb-2.0E0*1*I_ERI_F3y_S_D2z_S_b;
  abcd[357] = 4.0E0*I_ERI_F2yz_D2z_D2z_S_bb-2.0E0*1*I_ERI_F2yz_S_D2z_S_b;
  abcd[358] = 4.0E0*I_ERI_Fy2z_D2z_D2z_S_bb-2.0E0*1*I_ERI_Fy2z_S_D2z_S_b;
  abcd[359] = 4.0E0*I_ERI_F3z_D2z_D2z_S_bb-2.0E0*1*I_ERI_F3z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[360] = 4.0E0*I_ERI_F3x_Px_F3x_S_bc-2.0E0*2*I_ERI_F3x_Px_Px_S_b;
  abcd[361] = 4.0E0*I_ERI_F2xy_Px_F3x_S_bc-2.0E0*2*I_ERI_F2xy_Px_Px_S_b;
  abcd[362] = 4.0E0*I_ERI_F2xz_Px_F3x_S_bc-2.0E0*2*I_ERI_F2xz_Px_Px_S_b;
  abcd[363] = 4.0E0*I_ERI_Fx2y_Px_F3x_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Px_S_b;
  abcd[364] = 4.0E0*I_ERI_Fxyz_Px_F3x_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Px_S_b;
  abcd[365] = 4.0E0*I_ERI_Fx2z_Px_F3x_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Px_S_b;
  abcd[366] = 4.0E0*I_ERI_F3y_Px_F3x_S_bc-2.0E0*2*I_ERI_F3y_Px_Px_S_b;
  abcd[367] = 4.0E0*I_ERI_F2yz_Px_F3x_S_bc-2.0E0*2*I_ERI_F2yz_Px_Px_S_b;
  abcd[368] = 4.0E0*I_ERI_Fy2z_Px_F3x_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Px_S_b;
  abcd[369] = 4.0E0*I_ERI_F3z_Px_F3x_S_bc-2.0E0*2*I_ERI_F3z_Px_Px_S_b;
  abcd[370] = 4.0E0*I_ERI_F3x_Px_F2xy_S_bc-2.0E0*1*I_ERI_F3x_Px_Py_S_b;
  abcd[371] = 4.0E0*I_ERI_F2xy_Px_F2xy_S_bc-2.0E0*1*I_ERI_F2xy_Px_Py_S_b;
  abcd[372] = 4.0E0*I_ERI_F2xz_Px_F2xy_S_bc-2.0E0*1*I_ERI_F2xz_Px_Py_S_b;
  abcd[373] = 4.0E0*I_ERI_Fx2y_Px_F2xy_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Py_S_b;
  abcd[374] = 4.0E0*I_ERI_Fxyz_Px_F2xy_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Py_S_b;
  abcd[375] = 4.0E0*I_ERI_Fx2z_Px_F2xy_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Py_S_b;
  abcd[376] = 4.0E0*I_ERI_F3y_Px_F2xy_S_bc-2.0E0*1*I_ERI_F3y_Px_Py_S_b;
  abcd[377] = 4.0E0*I_ERI_F2yz_Px_F2xy_S_bc-2.0E0*1*I_ERI_F2yz_Px_Py_S_b;
  abcd[378] = 4.0E0*I_ERI_Fy2z_Px_F2xy_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Py_S_b;
  abcd[379] = 4.0E0*I_ERI_F3z_Px_F2xy_S_bc-2.0E0*1*I_ERI_F3z_Px_Py_S_b;
  abcd[380] = 4.0E0*I_ERI_F3x_Px_F2xz_S_bc-2.0E0*1*I_ERI_F3x_Px_Pz_S_b;
  abcd[381] = 4.0E0*I_ERI_F2xy_Px_F2xz_S_bc-2.0E0*1*I_ERI_F2xy_Px_Pz_S_b;
  abcd[382] = 4.0E0*I_ERI_F2xz_Px_F2xz_S_bc-2.0E0*1*I_ERI_F2xz_Px_Pz_S_b;
  abcd[383] = 4.0E0*I_ERI_Fx2y_Px_F2xz_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Pz_S_b;
  abcd[384] = 4.0E0*I_ERI_Fxyz_Px_F2xz_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Pz_S_b;
  abcd[385] = 4.0E0*I_ERI_Fx2z_Px_F2xz_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Pz_S_b;
  abcd[386] = 4.0E0*I_ERI_F3y_Px_F2xz_S_bc-2.0E0*1*I_ERI_F3y_Px_Pz_S_b;
  abcd[387] = 4.0E0*I_ERI_F2yz_Px_F2xz_S_bc-2.0E0*1*I_ERI_F2yz_Px_Pz_S_b;
  abcd[388] = 4.0E0*I_ERI_Fy2z_Px_F2xz_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Pz_S_b;
  abcd[389] = 4.0E0*I_ERI_F3z_Px_F2xz_S_bc-2.0E0*1*I_ERI_F3z_Px_Pz_S_b;
  abcd[390] = 4.0E0*I_ERI_F3x_Px_Fx2y_S_bc;
  abcd[391] = 4.0E0*I_ERI_F2xy_Px_Fx2y_S_bc;
  abcd[392] = 4.0E0*I_ERI_F2xz_Px_Fx2y_S_bc;
  abcd[393] = 4.0E0*I_ERI_Fx2y_Px_Fx2y_S_bc;
  abcd[394] = 4.0E0*I_ERI_Fxyz_Px_Fx2y_S_bc;
  abcd[395] = 4.0E0*I_ERI_Fx2z_Px_Fx2y_S_bc;
  abcd[396] = 4.0E0*I_ERI_F3y_Px_Fx2y_S_bc;
  abcd[397] = 4.0E0*I_ERI_F2yz_Px_Fx2y_S_bc;
  abcd[398] = 4.0E0*I_ERI_Fy2z_Px_Fx2y_S_bc;
  abcd[399] = 4.0E0*I_ERI_F3z_Px_Fx2y_S_bc;
  abcd[400] = 4.0E0*I_ERI_F3x_Px_Fxyz_S_bc;
  abcd[401] = 4.0E0*I_ERI_F2xy_Px_Fxyz_S_bc;
  abcd[402] = 4.0E0*I_ERI_F2xz_Px_Fxyz_S_bc;
  abcd[403] = 4.0E0*I_ERI_Fx2y_Px_Fxyz_S_bc;
  abcd[404] = 4.0E0*I_ERI_Fxyz_Px_Fxyz_S_bc;
  abcd[405] = 4.0E0*I_ERI_Fx2z_Px_Fxyz_S_bc;
  abcd[406] = 4.0E0*I_ERI_F3y_Px_Fxyz_S_bc;
  abcd[407] = 4.0E0*I_ERI_F2yz_Px_Fxyz_S_bc;
  abcd[408] = 4.0E0*I_ERI_Fy2z_Px_Fxyz_S_bc;
  abcd[409] = 4.0E0*I_ERI_F3z_Px_Fxyz_S_bc;
  abcd[410] = 4.0E0*I_ERI_F3x_Px_Fx2z_S_bc;
  abcd[411] = 4.0E0*I_ERI_F2xy_Px_Fx2z_S_bc;
  abcd[412] = 4.0E0*I_ERI_F2xz_Px_Fx2z_S_bc;
  abcd[413] = 4.0E0*I_ERI_Fx2y_Px_Fx2z_S_bc;
  abcd[414] = 4.0E0*I_ERI_Fxyz_Px_Fx2z_S_bc;
  abcd[415] = 4.0E0*I_ERI_Fx2z_Px_Fx2z_S_bc;
  abcd[416] = 4.0E0*I_ERI_F3y_Px_Fx2z_S_bc;
  abcd[417] = 4.0E0*I_ERI_F2yz_Px_Fx2z_S_bc;
  abcd[418] = 4.0E0*I_ERI_Fy2z_Px_Fx2z_S_bc;
  abcd[419] = 4.0E0*I_ERI_F3z_Px_Fx2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[420] = 4.0E0*I_ERI_F3x_Px_F2xy_S_bc;
  abcd[421] = 4.0E0*I_ERI_F2xy_Px_F2xy_S_bc;
  abcd[422] = 4.0E0*I_ERI_F2xz_Px_F2xy_S_bc;
  abcd[423] = 4.0E0*I_ERI_Fx2y_Px_F2xy_S_bc;
  abcd[424] = 4.0E0*I_ERI_Fxyz_Px_F2xy_S_bc;
  abcd[425] = 4.0E0*I_ERI_Fx2z_Px_F2xy_S_bc;
  abcd[426] = 4.0E0*I_ERI_F3y_Px_F2xy_S_bc;
  abcd[427] = 4.0E0*I_ERI_F2yz_Px_F2xy_S_bc;
  abcd[428] = 4.0E0*I_ERI_Fy2z_Px_F2xy_S_bc;
  abcd[429] = 4.0E0*I_ERI_F3z_Px_F2xy_S_bc;
  abcd[430] = 4.0E0*I_ERI_F3x_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F3x_Px_Px_S_b;
  abcd[431] = 4.0E0*I_ERI_F2xy_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F2xy_Px_Px_S_b;
  abcd[432] = 4.0E0*I_ERI_F2xz_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F2xz_Px_Px_S_b;
  abcd[433] = 4.0E0*I_ERI_Fx2y_Px_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Px_S_b;
  abcd[434] = 4.0E0*I_ERI_Fxyz_Px_Fx2y_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Px_S_b;
  abcd[435] = 4.0E0*I_ERI_Fx2z_Px_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Px_S_b;
  abcd[436] = 4.0E0*I_ERI_F3y_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F3y_Px_Px_S_b;
  abcd[437] = 4.0E0*I_ERI_F2yz_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F2yz_Px_Px_S_b;
  abcd[438] = 4.0E0*I_ERI_Fy2z_Px_Fx2y_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Px_S_b;
  abcd[439] = 4.0E0*I_ERI_F3z_Px_Fx2y_S_bc-2.0E0*1*I_ERI_F3z_Px_Px_S_b;
  abcd[440] = 4.0E0*I_ERI_F3x_Px_Fxyz_S_bc;
  abcd[441] = 4.0E0*I_ERI_F2xy_Px_Fxyz_S_bc;
  abcd[442] = 4.0E0*I_ERI_F2xz_Px_Fxyz_S_bc;
  abcd[443] = 4.0E0*I_ERI_Fx2y_Px_Fxyz_S_bc;
  abcd[444] = 4.0E0*I_ERI_Fxyz_Px_Fxyz_S_bc;
  abcd[445] = 4.0E0*I_ERI_Fx2z_Px_Fxyz_S_bc;
  abcd[446] = 4.0E0*I_ERI_F3y_Px_Fxyz_S_bc;
  abcd[447] = 4.0E0*I_ERI_F2yz_Px_Fxyz_S_bc;
  abcd[448] = 4.0E0*I_ERI_Fy2z_Px_Fxyz_S_bc;
  abcd[449] = 4.0E0*I_ERI_F3z_Px_Fxyz_S_bc;
  abcd[450] = 4.0E0*I_ERI_F3x_Px_F3y_S_bc-2.0E0*2*I_ERI_F3x_Px_Py_S_b;
  abcd[451] = 4.0E0*I_ERI_F2xy_Px_F3y_S_bc-2.0E0*2*I_ERI_F2xy_Px_Py_S_b;
  abcd[452] = 4.0E0*I_ERI_F2xz_Px_F3y_S_bc-2.0E0*2*I_ERI_F2xz_Px_Py_S_b;
  abcd[453] = 4.0E0*I_ERI_Fx2y_Px_F3y_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Py_S_b;
  abcd[454] = 4.0E0*I_ERI_Fxyz_Px_F3y_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Py_S_b;
  abcd[455] = 4.0E0*I_ERI_Fx2z_Px_F3y_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Py_S_b;
  abcd[456] = 4.0E0*I_ERI_F3y_Px_F3y_S_bc-2.0E0*2*I_ERI_F3y_Px_Py_S_b;
  abcd[457] = 4.0E0*I_ERI_F2yz_Px_F3y_S_bc-2.0E0*2*I_ERI_F2yz_Px_Py_S_b;
  abcd[458] = 4.0E0*I_ERI_Fy2z_Px_F3y_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Py_S_b;
  abcd[459] = 4.0E0*I_ERI_F3z_Px_F3y_S_bc-2.0E0*2*I_ERI_F3z_Px_Py_S_b;
  abcd[460] = 4.0E0*I_ERI_F3x_Px_F2yz_S_bc-2.0E0*1*I_ERI_F3x_Px_Pz_S_b;
  abcd[461] = 4.0E0*I_ERI_F2xy_Px_F2yz_S_bc-2.0E0*1*I_ERI_F2xy_Px_Pz_S_b;
  abcd[462] = 4.0E0*I_ERI_F2xz_Px_F2yz_S_bc-2.0E0*1*I_ERI_F2xz_Px_Pz_S_b;
  abcd[463] = 4.0E0*I_ERI_Fx2y_Px_F2yz_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Pz_S_b;
  abcd[464] = 4.0E0*I_ERI_Fxyz_Px_F2yz_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Pz_S_b;
  abcd[465] = 4.0E0*I_ERI_Fx2z_Px_F2yz_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Pz_S_b;
  abcd[466] = 4.0E0*I_ERI_F3y_Px_F2yz_S_bc-2.0E0*1*I_ERI_F3y_Px_Pz_S_b;
  abcd[467] = 4.0E0*I_ERI_F2yz_Px_F2yz_S_bc-2.0E0*1*I_ERI_F2yz_Px_Pz_S_b;
  abcd[468] = 4.0E0*I_ERI_Fy2z_Px_F2yz_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Pz_S_b;
  abcd[469] = 4.0E0*I_ERI_F3z_Px_F2yz_S_bc-2.0E0*1*I_ERI_F3z_Px_Pz_S_b;
  abcd[470] = 4.0E0*I_ERI_F3x_Px_Fy2z_S_bc;
  abcd[471] = 4.0E0*I_ERI_F2xy_Px_Fy2z_S_bc;
  abcd[472] = 4.0E0*I_ERI_F2xz_Px_Fy2z_S_bc;
  abcd[473] = 4.0E0*I_ERI_Fx2y_Px_Fy2z_S_bc;
  abcd[474] = 4.0E0*I_ERI_Fxyz_Px_Fy2z_S_bc;
  abcd[475] = 4.0E0*I_ERI_Fx2z_Px_Fy2z_S_bc;
  abcd[476] = 4.0E0*I_ERI_F3y_Px_Fy2z_S_bc;
  abcd[477] = 4.0E0*I_ERI_F2yz_Px_Fy2z_S_bc;
  abcd[478] = 4.0E0*I_ERI_Fy2z_Px_Fy2z_S_bc;
  abcd[479] = 4.0E0*I_ERI_F3z_Px_Fy2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[480] = 4.0E0*I_ERI_F3x_Px_F2xz_S_bc;
  abcd[481] = 4.0E0*I_ERI_F2xy_Px_F2xz_S_bc;
  abcd[482] = 4.0E0*I_ERI_F2xz_Px_F2xz_S_bc;
  abcd[483] = 4.0E0*I_ERI_Fx2y_Px_F2xz_S_bc;
  abcd[484] = 4.0E0*I_ERI_Fxyz_Px_F2xz_S_bc;
  abcd[485] = 4.0E0*I_ERI_Fx2z_Px_F2xz_S_bc;
  abcd[486] = 4.0E0*I_ERI_F3y_Px_F2xz_S_bc;
  abcd[487] = 4.0E0*I_ERI_F2yz_Px_F2xz_S_bc;
  abcd[488] = 4.0E0*I_ERI_Fy2z_Px_F2xz_S_bc;
  abcd[489] = 4.0E0*I_ERI_F3z_Px_F2xz_S_bc;
  abcd[490] = 4.0E0*I_ERI_F3x_Px_Fxyz_S_bc;
  abcd[491] = 4.0E0*I_ERI_F2xy_Px_Fxyz_S_bc;
  abcd[492] = 4.0E0*I_ERI_F2xz_Px_Fxyz_S_bc;
  abcd[493] = 4.0E0*I_ERI_Fx2y_Px_Fxyz_S_bc;
  abcd[494] = 4.0E0*I_ERI_Fxyz_Px_Fxyz_S_bc;
  abcd[495] = 4.0E0*I_ERI_Fx2z_Px_Fxyz_S_bc;
  abcd[496] = 4.0E0*I_ERI_F3y_Px_Fxyz_S_bc;
  abcd[497] = 4.0E0*I_ERI_F2yz_Px_Fxyz_S_bc;
  abcd[498] = 4.0E0*I_ERI_Fy2z_Px_Fxyz_S_bc;
  abcd[499] = 4.0E0*I_ERI_F3z_Px_Fxyz_S_bc;
  abcd[500] = 4.0E0*I_ERI_F3x_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F3x_Px_Px_S_b;
  abcd[501] = 4.0E0*I_ERI_F2xy_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F2xy_Px_Px_S_b;
  abcd[502] = 4.0E0*I_ERI_F2xz_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F2xz_Px_Px_S_b;
  abcd[503] = 4.0E0*I_ERI_Fx2y_Px_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Px_S_b;
  abcd[504] = 4.0E0*I_ERI_Fxyz_Px_Fx2z_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Px_S_b;
  abcd[505] = 4.0E0*I_ERI_Fx2z_Px_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Px_S_b;
  abcd[506] = 4.0E0*I_ERI_F3y_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F3y_Px_Px_S_b;
  abcd[507] = 4.0E0*I_ERI_F2yz_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F2yz_Px_Px_S_b;
  abcd[508] = 4.0E0*I_ERI_Fy2z_Px_Fx2z_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Px_S_b;
  abcd[509] = 4.0E0*I_ERI_F3z_Px_Fx2z_S_bc-2.0E0*1*I_ERI_F3z_Px_Px_S_b;
  abcd[510] = 4.0E0*I_ERI_F3x_Px_F2yz_S_bc;
  abcd[511] = 4.0E0*I_ERI_F2xy_Px_F2yz_S_bc;
  abcd[512] = 4.0E0*I_ERI_F2xz_Px_F2yz_S_bc;
  abcd[513] = 4.0E0*I_ERI_Fx2y_Px_F2yz_S_bc;
  abcd[514] = 4.0E0*I_ERI_Fxyz_Px_F2yz_S_bc;
  abcd[515] = 4.0E0*I_ERI_Fx2z_Px_F2yz_S_bc;
  abcd[516] = 4.0E0*I_ERI_F3y_Px_F2yz_S_bc;
  abcd[517] = 4.0E0*I_ERI_F2yz_Px_F2yz_S_bc;
  abcd[518] = 4.0E0*I_ERI_Fy2z_Px_F2yz_S_bc;
  abcd[519] = 4.0E0*I_ERI_F3z_Px_F2yz_S_bc;
  abcd[520] = 4.0E0*I_ERI_F3x_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F3x_Px_Py_S_b;
  abcd[521] = 4.0E0*I_ERI_F2xy_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F2xy_Px_Py_S_b;
  abcd[522] = 4.0E0*I_ERI_F2xz_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F2xz_Px_Py_S_b;
  abcd[523] = 4.0E0*I_ERI_Fx2y_Px_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Py_S_b;
  abcd[524] = 4.0E0*I_ERI_Fxyz_Px_Fy2z_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Py_S_b;
  abcd[525] = 4.0E0*I_ERI_Fx2z_Px_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Py_S_b;
  abcd[526] = 4.0E0*I_ERI_F3y_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F3y_Px_Py_S_b;
  abcd[527] = 4.0E0*I_ERI_F2yz_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F2yz_Px_Py_S_b;
  abcd[528] = 4.0E0*I_ERI_Fy2z_Px_Fy2z_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Py_S_b;
  abcd[529] = 4.0E0*I_ERI_F3z_Px_Fy2z_S_bc-2.0E0*1*I_ERI_F3z_Px_Py_S_b;
  abcd[530] = 4.0E0*I_ERI_F3x_Px_F3z_S_bc-2.0E0*2*I_ERI_F3x_Px_Pz_S_b;
  abcd[531] = 4.0E0*I_ERI_F2xy_Px_F3z_S_bc-2.0E0*2*I_ERI_F2xy_Px_Pz_S_b;
  abcd[532] = 4.0E0*I_ERI_F2xz_Px_F3z_S_bc-2.0E0*2*I_ERI_F2xz_Px_Pz_S_b;
  abcd[533] = 4.0E0*I_ERI_Fx2y_Px_F3z_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Pz_S_b;
  abcd[534] = 4.0E0*I_ERI_Fxyz_Px_F3z_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Pz_S_b;
  abcd[535] = 4.0E0*I_ERI_Fx2z_Px_F3z_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Pz_S_b;
  abcd[536] = 4.0E0*I_ERI_F3y_Px_F3z_S_bc-2.0E0*2*I_ERI_F3y_Px_Pz_S_b;
  abcd[537] = 4.0E0*I_ERI_F2yz_Px_F3z_S_bc-2.0E0*2*I_ERI_F2yz_Px_Pz_S_b;
  abcd[538] = 4.0E0*I_ERI_Fy2z_Px_F3z_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Pz_S_b;
  abcd[539] = 4.0E0*I_ERI_F3z_Px_F3z_S_bc-2.0E0*2*I_ERI_F3z_Px_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[540] = 4.0E0*I_ERI_F3x_Py_F3x_S_bc-2.0E0*2*I_ERI_F3x_Py_Px_S_b;
  abcd[541] = 4.0E0*I_ERI_F2xy_Py_F3x_S_bc-2.0E0*2*I_ERI_F2xy_Py_Px_S_b;
  abcd[542] = 4.0E0*I_ERI_F2xz_Py_F3x_S_bc-2.0E0*2*I_ERI_F2xz_Py_Px_S_b;
  abcd[543] = 4.0E0*I_ERI_Fx2y_Py_F3x_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Px_S_b;
  abcd[544] = 4.0E0*I_ERI_Fxyz_Py_F3x_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Px_S_b;
  abcd[545] = 4.0E0*I_ERI_Fx2z_Py_F3x_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Px_S_b;
  abcd[546] = 4.0E0*I_ERI_F3y_Py_F3x_S_bc-2.0E0*2*I_ERI_F3y_Py_Px_S_b;
  abcd[547] = 4.0E0*I_ERI_F2yz_Py_F3x_S_bc-2.0E0*2*I_ERI_F2yz_Py_Px_S_b;
  abcd[548] = 4.0E0*I_ERI_Fy2z_Py_F3x_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Px_S_b;
  abcd[549] = 4.0E0*I_ERI_F3z_Py_F3x_S_bc-2.0E0*2*I_ERI_F3z_Py_Px_S_b;
  abcd[550] = 4.0E0*I_ERI_F3x_Py_F2xy_S_bc-2.0E0*1*I_ERI_F3x_Py_Py_S_b;
  abcd[551] = 4.0E0*I_ERI_F2xy_Py_F2xy_S_bc-2.0E0*1*I_ERI_F2xy_Py_Py_S_b;
  abcd[552] = 4.0E0*I_ERI_F2xz_Py_F2xy_S_bc-2.0E0*1*I_ERI_F2xz_Py_Py_S_b;
  abcd[553] = 4.0E0*I_ERI_Fx2y_Py_F2xy_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Py_S_b;
  abcd[554] = 4.0E0*I_ERI_Fxyz_Py_F2xy_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Py_S_b;
  abcd[555] = 4.0E0*I_ERI_Fx2z_Py_F2xy_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Py_S_b;
  abcd[556] = 4.0E0*I_ERI_F3y_Py_F2xy_S_bc-2.0E0*1*I_ERI_F3y_Py_Py_S_b;
  abcd[557] = 4.0E0*I_ERI_F2yz_Py_F2xy_S_bc-2.0E0*1*I_ERI_F2yz_Py_Py_S_b;
  abcd[558] = 4.0E0*I_ERI_Fy2z_Py_F2xy_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Py_S_b;
  abcd[559] = 4.0E0*I_ERI_F3z_Py_F2xy_S_bc-2.0E0*1*I_ERI_F3z_Py_Py_S_b;
  abcd[560] = 4.0E0*I_ERI_F3x_Py_F2xz_S_bc-2.0E0*1*I_ERI_F3x_Py_Pz_S_b;
  abcd[561] = 4.0E0*I_ERI_F2xy_Py_F2xz_S_bc-2.0E0*1*I_ERI_F2xy_Py_Pz_S_b;
  abcd[562] = 4.0E0*I_ERI_F2xz_Py_F2xz_S_bc-2.0E0*1*I_ERI_F2xz_Py_Pz_S_b;
  abcd[563] = 4.0E0*I_ERI_Fx2y_Py_F2xz_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Pz_S_b;
  abcd[564] = 4.0E0*I_ERI_Fxyz_Py_F2xz_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Pz_S_b;
  abcd[565] = 4.0E0*I_ERI_Fx2z_Py_F2xz_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Pz_S_b;
  abcd[566] = 4.0E0*I_ERI_F3y_Py_F2xz_S_bc-2.0E0*1*I_ERI_F3y_Py_Pz_S_b;
  abcd[567] = 4.0E0*I_ERI_F2yz_Py_F2xz_S_bc-2.0E0*1*I_ERI_F2yz_Py_Pz_S_b;
  abcd[568] = 4.0E0*I_ERI_Fy2z_Py_F2xz_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Pz_S_b;
  abcd[569] = 4.0E0*I_ERI_F3z_Py_F2xz_S_bc-2.0E0*1*I_ERI_F3z_Py_Pz_S_b;
  abcd[570] = 4.0E0*I_ERI_F3x_Py_Fx2y_S_bc;
  abcd[571] = 4.0E0*I_ERI_F2xy_Py_Fx2y_S_bc;
  abcd[572] = 4.0E0*I_ERI_F2xz_Py_Fx2y_S_bc;
  abcd[573] = 4.0E0*I_ERI_Fx2y_Py_Fx2y_S_bc;
  abcd[574] = 4.0E0*I_ERI_Fxyz_Py_Fx2y_S_bc;
  abcd[575] = 4.0E0*I_ERI_Fx2z_Py_Fx2y_S_bc;
  abcd[576] = 4.0E0*I_ERI_F3y_Py_Fx2y_S_bc;
  abcd[577] = 4.0E0*I_ERI_F2yz_Py_Fx2y_S_bc;
  abcd[578] = 4.0E0*I_ERI_Fy2z_Py_Fx2y_S_bc;
  abcd[579] = 4.0E0*I_ERI_F3z_Py_Fx2y_S_bc;
  abcd[580] = 4.0E0*I_ERI_F3x_Py_Fxyz_S_bc;
  abcd[581] = 4.0E0*I_ERI_F2xy_Py_Fxyz_S_bc;
  abcd[582] = 4.0E0*I_ERI_F2xz_Py_Fxyz_S_bc;
  abcd[583] = 4.0E0*I_ERI_Fx2y_Py_Fxyz_S_bc;
  abcd[584] = 4.0E0*I_ERI_Fxyz_Py_Fxyz_S_bc;
  abcd[585] = 4.0E0*I_ERI_Fx2z_Py_Fxyz_S_bc;
  abcd[586] = 4.0E0*I_ERI_F3y_Py_Fxyz_S_bc;
  abcd[587] = 4.0E0*I_ERI_F2yz_Py_Fxyz_S_bc;
  abcd[588] = 4.0E0*I_ERI_Fy2z_Py_Fxyz_S_bc;
  abcd[589] = 4.0E0*I_ERI_F3z_Py_Fxyz_S_bc;
  abcd[590] = 4.0E0*I_ERI_F3x_Py_Fx2z_S_bc;
  abcd[591] = 4.0E0*I_ERI_F2xy_Py_Fx2z_S_bc;
  abcd[592] = 4.0E0*I_ERI_F2xz_Py_Fx2z_S_bc;
  abcd[593] = 4.0E0*I_ERI_Fx2y_Py_Fx2z_S_bc;
  abcd[594] = 4.0E0*I_ERI_Fxyz_Py_Fx2z_S_bc;
  abcd[595] = 4.0E0*I_ERI_Fx2z_Py_Fx2z_S_bc;
  abcd[596] = 4.0E0*I_ERI_F3y_Py_Fx2z_S_bc;
  abcd[597] = 4.0E0*I_ERI_F2yz_Py_Fx2z_S_bc;
  abcd[598] = 4.0E0*I_ERI_Fy2z_Py_Fx2z_S_bc;
  abcd[599] = 4.0E0*I_ERI_F3z_Py_Fx2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[600] = 4.0E0*I_ERI_F3x_Py_F2xy_S_bc;
  abcd[601] = 4.0E0*I_ERI_F2xy_Py_F2xy_S_bc;
  abcd[602] = 4.0E0*I_ERI_F2xz_Py_F2xy_S_bc;
  abcd[603] = 4.0E0*I_ERI_Fx2y_Py_F2xy_S_bc;
  abcd[604] = 4.0E0*I_ERI_Fxyz_Py_F2xy_S_bc;
  abcd[605] = 4.0E0*I_ERI_Fx2z_Py_F2xy_S_bc;
  abcd[606] = 4.0E0*I_ERI_F3y_Py_F2xy_S_bc;
  abcd[607] = 4.0E0*I_ERI_F2yz_Py_F2xy_S_bc;
  abcd[608] = 4.0E0*I_ERI_Fy2z_Py_F2xy_S_bc;
  abcd[609] = 4.0E0*I_ERI_F3z_Py_F2xy_S_bc;
  abcd[610] = 4.0E0*I_ERI_F3x_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F3x_Py_Px_S_b;
  abcd[611] = 4.0E0*I_ERI_F2xy_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F2xy_Py_Px_S_b;
  abcd[612] = 4.0E0*I_ERI_F2xz_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F2xz_Py_Px_S_b;
  abcd[613] = 4.0E0*I_ERI_Fx2y_Py_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Px_S_b;
  abcd[614] = 4.0E0*I_ERI_Fxyz_Py_Fx2y_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Px_S_b;
  abcd[615] = 4.0E0*I_ERI_Fx2z_Py_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Px_S_b;
  abcd[616] = 4.0E0*I_ERI_F3y_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F3y_Py_Px_S_b;
  abcd[617] = 4.0E0*I_ERI_F2yz_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F2yz_Py_Px_S_b;
  abcd[618] = 4.0E0*I_ERI_Fy2z_Py_Fx2y_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Px_S_b;
  abcd[619] = 4.0E0*I_ERI_F3z_Py_Fx2y_S_bc-2.0E0*1*I_ERI_F3z_Py_Px_S_b;
  abcd[620] = 4.0E0*I_ERI_F3x_Py_Fxyz_S_bc;
  abcd[621] = 4.0E0*I_ERI_F2xy_Py_Fxyz_S_bc;
  abcd[622] = 4.0E0*I_ERI_F2xz_Py_Fxyz_S_bc;
  abcd[623] = 4.0E0*I_ERI_Fx2y_Py_Fxyz_S_bc;
  abcd[624] = 4.0E0*I_ERI_Fxyz_Py_Fxyz_S_bc;
  abcd[625] = 4.0E0*I_ERI_Fx2z_Py_Fxyz_S_bc;
  abcd[626] = 4.0E0*I_ERI_F3y_Py_Fxyz_S_bc;
  abcd[627] = 4.0E0*I_ERI_F2yz_Py_Fxyz_S_bc;
  abcd[628] = 4.0E0*I_ERI_Fy2z_Py_Fxyz_S_bc;
  abcd[629] = 4.0E0*I_ERI_F3z_Py_Fxyz_S_bc;
  abcd[630] = 4.0E0*I_ERI_F3x_Py_F3y_S_bc-2.0E0*2*I_ERI_F3x_Py_Py_S_b;
  abcd[631] = 4.0E0*I_ERI_F2xy_Py_F3y_S_bc-2.0E0*2*I_ERI_F2xy_Py_Py_S_b;
  abcd[632] = 4.0E0*I_ERI_F2xz_Py_F3y_S_bc-2.0E0*2*I_ERI_F2xz_Py_Py_S_b;
  abcd[633] = 4.0E0*I_ERI_Fx2y_Py_F3y_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Py_S_b;
  abcd[634] = 4.0E0*I_ERI_Fxyz_Py_F3y_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Py_S_b;
  abcd[635] = 4.0E0*I_ERI_Fx2z_Py_F3y_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Py_S_b;
  abcd[636] = 4.0E0*I_ERI_F3y_Py_F3y_S_bc-2.0E0*2*I_ERI_F3y_Py_Py_S_b;
  abcd[637] = 4.0E0*I_ERI_F2yz_Py_F3y_S_bc-2.0E0*2*I_ERI_F2yz_Py_Py_S_b;
  abcd[638] = 4.0E0*I_ERI_Fy2z_Py_F3y_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Py_S_b;
  abcd[639] = 4.0E0*I_ERI_F3z_Py_F3y_S_bc-2.0E0*2*I_ERI_F3z_Py_Py_S_b;
  abcd[640] = 4.0E0*I_ERI_F3x_Py_F2yz_S_bc-2.0E0*1*I_ERI_F3x_Py_Pz_S_b;
  abcd[641] = 4.0E0*I_ERI_F2xy_Py_F2yz_S_bc-2.0E0*1*I_ERI_F2xy_Py_Pz_S_b;
  abcd[642] = 4.0E0*I_ERI_F2xz_Py_F2yz_S_bc-2.0E0*1*I_ERI_F2xz_Py_Pz_S_b;
  abcd[643] = 4.0E0*I_ERI_Fx2y_Py_F2yz_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Pz_S_b;
  abcd[644] = 4.0E0*I_ERI_Fxyz_Py_F2yz_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Pz_S_b;
  abcd[645] = 4.0E0*I_ERI_Fx2z_Py_F2yz_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Pz_S_b;
  abcd[646] = 4.0E0*I_ERI_F3y_Py_F2yz_S_bc-2.0E0*1*I_ERI_F3y_Py_Pz_S_b;
  abcd[647] = 4.0E0*I_ERI_F2yz_Py_F2yz_S_bc-2.0E0*1*I_ERI_F2yz_Py_Pz_S_b;
  abcd[648] = 4.0E0*I_ERI_Fy2z_Py_F2yz_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Pz_S_b;
  abcd[649] = 4.0E0*I_ERI_F3z_Py_F2yz_S_bc-2.0E0*1*I_ERI_F3z_Py_Pz_S_b;
  abcd[650] = 4.0E0*I_ERI_F3x_Py_Fy2z_S_bc;
  abcd[651] = 4.0E0*I_ERI_F2xy_Py_Fy2z_S_bc;
  abcd[652] = 4.0E0*I_ERI_F2xz_Py_Fy2z_S_bc;
  abcd[653] = 4.0E0*I_ERI_Fx2y_Py_Fy2z_S_bc;
  abcd[654] = 4.0E0*I_ERI_Fxyz_Py_Fy2z_S_bc;
  abcd[655] = 4.0E0*I_ERI_Fx2z_Py_Fy2z_S_bc;
  abcd[656] = 4.0E0*I_ERI_F3y_Py_Fy2z_S_bc;
  abcd[657] = 4.0E0*I_ERI_F2yz_Py_Fy2z_S_bc;
  abcd[658] = 4.0E0*I_ERI_Fy2z_Py_Fy2z_S_bc;
  abcd[659] = 4.0E0*I_ERI_F3z_Py_Fy2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[660] = 4.0E0*I_ERI_F3x_Py_F2xz_S_bc;
  abcd[661] = 4.0E0*I_ERI_F2xy_Py_F2xz_S_bc;
  abcd[662] = 4.0E0*I_ERI_F2xz_Py_F2xz_S_bc;
  abcd[663] = 4.0E0*I_ERI_Fx2y_Py_F2xz_S_bc;
  abcd[664] = 4.0E0*I_ERI_Fxyz_Py_F2xz_S_bc;
  abcd[665] = 4.0E0*I_ERI_Fx2z_Py_F2xz_S_bc;
  abcd[666] = 4.0E0*I_ERI_F3y_Py_F2xz_S_bc;
  abcd[667] = 4.0E0*I_ERI_F2yz_Py_F2xz_S_bc;
  abcd[668] = 4.0E0*I_ERI_Fy2z_Py_F2xz_S_bc;
  abcd[669] = 4.0E0*I_ERI_F3z_Py_F2xz_S_bc;
  abcd[670] = 4.0E0*I_ERI_F3x_Py_Fxyz_S_bc;
  abcd[671] = 4.0E0*I_ERI_F2xy_Py_Fxyz_S_bc;
  abcd[672] = 4.0E0*I_ERI_F2xz_Py_Fxyz_S_bc;
  abcd[673] = 4.0E0*I_ERI_Fx2y_Py_Fxyz_S_bc;
  abcd[674] = 4.0E0*I_ERI_Fxyz_Py_Fxyz_S_bc;
  abcd[675] = 4.0E0*I_ERI_Fx2z_Py_Fxyz_S_bc;
  abcd[676] = 4.0E0*I_ERI_F3y_Py_Fxyz_S_bc;
  abcd[677] = 4.0E0*I_ERI_F2yz_Py_Fxyz_S_bc;
  abcd[678] = 4.0E0*I_ERI_Fy2z_Py_Fxyz_S_bc;
  abcd[679] = 4.0E0*I_ERI_F3z_Py_Fxyz_S_bc;
  abcd[680] = 4.0E0*I_ERI_F3x_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F3x_Py_Px_S_b;
  abcd[681] = 4.0E0*I_ERI_F2xy_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F2xy_Py_Px_S_b;
  abcd[682] = 4.0E0*I_ERI_F2xz_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F2xz_Py_Px_S_b;
  abcd[683] = 4.0E0*I_ERI_Fx2y_Py_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Px_S_b;
  abcd[684] = 4.0E0*I_ERI_Fxyz_Py_Fx2z_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Px_S_b;
  abcd[685] = 4.0E0*I_ERI_Fx2z_Py_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Px_S_b;
  abcd[686] = 4.0E0*I_ERI_F3y_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F3y_Py_Px_S_b;
  abcd[687] = 4.0E0*I_ERI_F2yz_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F2yz_Py_Px_S_b;
  abcd[688] = 4.0E0*I_ERI_Fy2z_Py_Fx2z_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Px_S_b;
  abcd[689] = 4.0E0*I_ERI_F3z_Py_Fx2z_S_bc-2.0E0*1*I_ERI_F3z_Py_Px_S_b;
  abcd[690] = 4.0E0*I_ERI_F3x_Py_F2yz_S_bc;
  abcd[691] = 4.0E0*I_ERI_F2xy_Py_F2yz_S_bc;
  abcd[692] = 4.0E0*I_ERI_F2xz_Py_F2yz_S_bc;
  abcd[693] = 4.0E0*I_ERI_Fx2y_Py_F2yz_S_bc;
  abcd[694] = 4.0E0*I_ERI_Fxyz_Py_F2yz_S_bc;
  abcd[695] = 4.0E0*I_ERI_Fx2z_Py_F2yz_S_bc;
  abcd[696] = 4.0E0*I_ERI_F3y_Py_F2yz_S_bc;
  abcd[697] = 4.0E0*I_ERI_F2yz_Py_F2yz_S_bc;
  abcd[698] = 4.0E0*I_ERI_Fy2z_Py_F2yz_S_bc;
  abcd[699] = 4.0E0*I_ERI_F3z_Py_F2yz_S_bc;
  abcd[700] = 4.0E0*I_ERI_F3x_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F3x_Py_Py_S_b;
  abcd[701] = 4.0E0*I_ERI_F2xy_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F2xy_Py_Py_S_b;
  abcd[702] = 4.0E0*I_ERI_F2xz_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F2xz_Py_Py_S_b;
  abcd[703] = 4.0E0*I_ERI_Fx2y_Py_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Py_S_b;
  abcd[704] = 4.0E0*I_ERI_Fxyz_Py_Fy2z_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Py_S_b;
  abcd[705] = 4.0E0*I_ERI_Fx2z_Py_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Py_S_b;
  abcd[706] = 4.0E0*I_ERI_F3y_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F3y_Py_Py_S_b;
  abcd[707] = 4.0E0*I_ERI_F2yz_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F2yz_Py_Py_S_b;
  abcd[708] = 4.0E0*I_ERI_Fy2z_Py_Fy2z_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Py_S_b;
  abcd[709] = 4.0E0*I_ERI_F3z_Py_Fy2z_S_bc-2.0E0*1*I_ERI_F3z_Py_Py_S_b;
  abcd[710] = 4.0E0*I_ERI_F3x_Py_F3z_S_bc-2.0E0*2*I_ERI_F3x_Py_Pz_S_b;
  abcd[711] = 4.0E0*I_ERI_F2xy_Py_F3z_S_bc-2.0E0*2*I_ERI_F2xy_Py_Pz_S_b;
  abcd[712] = 4.0E0*I_ERI_F2xz_Py_F3z_S_bc-2.0E0*2*I_ERI_F2xz_Py_Pz_S_b;
  abcd[713] = 4.0E0*I_ERI_Fx2y_Py_F3z_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Pz_S_b;
  abcd[714] = 4.0E0*I_ERI_Fxyz_Py_F3z_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Pz_S_b;
  abcd[715] = 4.0E0*I_ERI_Fx2z_Py_F3z_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Pz_S_b;
  abcd[716] = 4.0E0*I_ERI_F3y_Py_F3z_S_bc-2.0E0*2*I_ERI_F3y_Py_Pz_S_b;
  abcd[717] = 4.0E0*I_ERI_F2yz_Py_F3z_S_bc-2.0E0*2*I_ERI_F2yz_Py_Pz_S_b;
  abcd[718] = 4.0E0*I_ERI_Fy2z_Py_F3z_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Pz_S_b;
  abcd[719] = 4.0E0*I_ERI_F3z_Py_F3z_S_bc-2.0E0*2*I_ERI_F3z_Py_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[720] = 4.0E0*I_ERI_F3x_Pz_F3x_S_bc-2.0E0*2*I_ERI_F3x_Pz_Px_S_b;
  abcd[721] = 4.0E0*I_ERI_F2xy_Pz_F3x_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Px_S_b;
  abcd[722] = 4.0E0*I_ERI_F2xz_Pz_F3x_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Px_S_b;
  abcd[723] = 4.0E0*I_ERI_Fx2y_Pz_F3x_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Px_S_b;
  abcd[724] = 4.0E0*I_ERI_Fxyz_Pz_F3x_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Px_S_b;
  abcd[725] = 4.0E0*I_ERI_Fx2z_Pz_F3x_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Px_S_b;
  abcd[726] = 4.0E0*I_ERI_F3y_Pz_F3x_S_bc-2.0E0*2*I_ERI_F3y_Pz_Px_S_b;
  abcd[727] = 4.0E0*I_ERI_F2yz_Pz_F3x_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Px_S_b;
  abcd[728] = 4.0E0*I_ERI_Fy2z_Pz_F3x_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Px_S_b;
  abcd[729] = 4.0E0*I_ERI_F3z_Pz_F3x_S_bc-2.0E0*2*I_ERI_F3z_Pz_Px_S_b;
  abcd[730] = 4.0E0*I_ERI_F3x_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F3x_Pz_Py_S_b;
  abcd[731] = 4.0E0*I_ERI_F2xy_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Py_S_b;
  abcd[732] = 4.0E0*I_ERI_F2xz_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Py_S_b;
  abcd[733] = 4.0E0*I_ERI_Fx2y_Pz_F2xy_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Py_S_b;
  abcd[734] = 4.0E0*I_ERI_Fxyz_Pz_F2xy_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Py_S_b;
  abcd[735] = 4.0E0*I_ERI_Fx2z_Pz_F2xy_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Py_S_b;
  abcd[736] = 4.0E0*I_ERI_F3y_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F3y_Pz_Py_S_b;
  abcd[737] = 4.0E0*I_ERI_F2yz_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Py_S_b;
  abcd[738] = 4.0E0*I_ERI_Fy2z_Pz_F2xy_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Py_S_b;
  abcd[739] = 4.0E0*I_ERI_F3z_Pz_F2xy_S_bc-2.0E0*1*I_ERI_F3z_Pz_Py_S_b;
  abcd[740] = 4.0E0*I_ERI_F3x_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F3x_Pz_Pz_S_b;
  abcd[741] = 4.0E0*I_ERI_F2xy_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Pz_S_b;
  abcd[742] = 4.0E0*I_ERI_F2xz_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Pz_S_b;
  abcd[743] = 4.0E0*I_ERI_Fx2y_Pz_F2xz_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Pz_S_b;
  abcd[744] = 4.0E0*I_ERI_Fxyz_Pz_F2xz_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Pz_S_b;
  abcd[745] = 4.0E0*I_ERI_Fx2z_Pz_F2xz_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Pz_S_b;
  abcd[746] = 4.0E0*I_ERI_F3y_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F3y_Pz_Pz_S_b;
  abcd[747] = 4.0E0*I_ERI_F2yz_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Pz_S_b;
  abcd[748] = 4.0E0*I_ERI_Fy2z_Pz_F2xz_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Pz_S_b;
  abcd[749] = 4.0E0*I_ERI_F3z_Pz_F2xz_S_bc-2.0E0*1*I_ERI_F3z_Pz_Pz_S_b;
  abcd[750] = 4.0E0*I_ERI_F3x_Pz_Fx2y_S_bc;
  abcd[751] = 4.0E0*I_ERI_F2xy_Pz_Fx2y_S_bc;
  abcd[752] = 4.0E0*I_ERI_F2xz_Pz_Fx2y_S_bc;
  abcd[753] = 4.0E0*I_ERI_Fx2y_Pz_Fx2y_S_bc;
  abcd[754] = 4.0E0*I_ERI_Fxyz_Pz_Fx2y_S_bc;
  abcd[755] = 4.0E0*I_ERI_Fx2z_Pz_Fx2y_S_bc;
  abcd[756] = 4.0E0*I_ERI_F3y_Pz_Fx2y_S_bc;
  abcd[757] = 4.0E0*I_ERI_F2yz_Pz_Fx2y_S_bc;
  abcd[758] = 4.0E0*I_ERI_Fy2z_Pz_Fx2y_S_bc;
  abcd[759] = 4.0E0*I_ERI_F3z_Pz_Fx2y_S_bc;
  abcd[760] = 4.0E0*I_ERI_F3x_Pz_Fxyz_S_bc;
  abcd[761] = 4.0E0*I_ERI_F2xy_Pz_Fxyz_S_bc;
  abcd[762] = 4.0E0*I_ERI_F2xz_Pz_Fxyz_S_bc;
  abcd[763] = 4.0E0*I_ERI_Fx2y_Pz_Fxyz_S_bc;
  abcd[764] = 4.0E0*I_ERI_Fxyz_Pz_Fxyz_S_bc;
  abcd[765] = 4.0E0*I_ERI_Fx2z_Pz_Fxyz_S_bc;
  abcd[766] = 4.0E0*I_ERI_F3y_Pz_Fxyz_S_bc;
  abcd[767] = 4.0E0*I_ERI_F2yz_Pz_Fxyz_S_bc;
  abcd[768] = 4.0E0*I_ERI_Fy2z_Pz_Fxyz_S_bc;
  abcd[769] = 4.0E0*I_ERI_F3z_Pz_Fxyz_S_bc;
  abcd[770] = 4.0E0*I_ERI_F3x_Pz_Fx2z_S_bc;
  abcd[771] = 4.0E0*I_ERI_F2xy_Pz_Fx2z_S_bc;
  abcd[772] = 4.0E0*I_ERI_F2xz_Pz_Fx2z_S_bc;
  abcd[773] = 4.0E0*I_ERI_Fx2y_Pz_Fx2z_S_bc;
  abcd[774] = 4.0E0*I_ERI_Fxyz_Pz_Fx2z_S_bc;
  abcd[775] = 4.0E0*I_ERI_Fx2z_Pz_Fx2z_S_bc;
  abcd[776] = 4.0E0*I_ERI_F3y_Pz_Fx2z_S_bc;
  abcd[777] = 4.0E0*I_ERI_F2yz_Pz_Fx2z_S_bc;
  abcd[778] = 4.0E0*I_ERI_Fy2z_Pz_Fx2z_S_bc;
  abcd[779] = 4.0E0*I_ERI_F3z_Pz_Fx2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[780] = 4.0E0*I_ERI_F3x_Pz_F2xy_S_bc;
  abcd[781] = 4.0E0*I_ERI_F2xy_Pz_F2xy_S_bc;
  abcd[782] = 4.0E0*I_ERI_F2xz_Pz_F2xy_S_bc;
  abcd[783] = 4.0E0*I_ERI_Fx2y_Pz_F2xy_S_bc;
  abcd[784] = 4.0E0*I_ERI_Fxyz_Pz_F2xy_S_bc;
  abcd[785] = 4.0E0*I_ERI_Fx2z_Pz_F2xy_S_bc;
  abcd[786] = 4.0E0*I_ERI_F3y_Pz_F2xy_S_bc;
  abcd[787] = 4.0E0*I_ERI_F2yz_Pz_F2xy_S_bc;
  abcd[788] = 4.0E0*I_ERI_Fy2z_Pz_F2xy_S_bc;
  abcd[789] = 4.0E0*I_ERI_F3z_Pz_F2xy_S_bc;
  abcd[790] = 4.0E0*I_ERI_F3x_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F3x_Pz_Px_S_b;
  abcd[791] = 4.0E0*I_ERI_F2xy_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Px_S_b;
  abcd[792] = 4.0E0*I_ERI_F2xz_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Px_S_b;
  abcd[793] = 4.0E0*I_ERI_Fx2y_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Px_S_b;
  abcd[794] = 4.0E0*I_ERI_Fxyz_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Px_S_b;
  abcd[795] = 4.0E0*I_ERI_Fx2z_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Px_S_b;
  abcd[796] = 4.0E0*I_ERI_F3y_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F3y_Pz_Px_S_b;
  abcd[797] = 4.0E0*I_ERI_F2yz_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Px_S_b;
  abcd[798] = 4.0E0*I_ERI_Fy2z_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Px_S_b;
  abcd[799] = 4.0E0*I_ERI_F3z_Pz_Fx2y_S_bc-2.0E0*1*I_ERI_F3z_Pz_Px_S_b;
  abcd[800] = 4.0E0*I_ERI_F3x_Pz_Fxyz_S_bc;
  abcd[801] = 4.0E0*I_ERI_F2xy_Pz_Fxyz_S_bc;
  abcd[802] = 4.0E0*I_ERI_F2xz_Pz_Fxyz_S_bc;
  abcd[803] = 4.0E0*I_ERI_Fx2y_Pz_Fxyz_S_bc;
  abcd[804] = 4.0E0*I_ERI_Fxyz_Pz_Fxyz_S_bc;
  abcd[805] = 4.0E0*I_ERI_Fx2z_Pz_Fxyz_S_bc;
  abcd[806] = 4.0E0*I_ERI_F3y_Pz_Fxyz_S_bc;
  abcd[807] = 4.0E0*I_ERI_F2yz_Pz_Fxyz_S_bc;
  abcd[808] = 4.0E0*I_ERI_Fy2z_Pz_Fxyz_S_bc;
  abcd[809] = 4.0E0*I_ERI_F3z_Pz_Fxyz_S_bc;
  abcd[810] = 4.0E0*I_ERI_F3x_Pz_F3y_S_bc-2.0E0*2*I_ERI_F3x_Pz_Py_S_b;
  abcd[811] = 4.0E0*I_ERI_F2xy_Pz_F3y_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Py_S_b;
  abcd[812] = 4.0E0*I_ERI_F2xz_Pz_F3y_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Py_S_b;
  abcd[813] = 4.0E0*I_ERI_Fx2y_Pz_F3y_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Py_S_b;
  abcd[814] = 4.0E0*I_ERI_Fxyz_Pz_F3y_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Py_S_b;
  abcd[815] = 4.0E0*I_ERI_Fx2z_Pz_F3y_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Py_S_b;
  abcd[816] = 4.0E0*I_ERI_F3y_Pz_F3y_S_bc-2.0E0*2*I_ERI_F3y_Pz_Py_S_b;
  abcd[817] = 4.0E0*I_ERI_F2yz_Pz_F3y_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Py_S_b;
  abcd[818] = 4.0E0*I_ERI_Fy2z_Pz_F3y_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Py_S_b;
  abcd[819] = 4.0E0*I_ERI_F3z_Pz_F3y_S_bc-2.0E0*2*I_ERI_F3z_Pz_Py_S_b;
  abcd[820] = 4.0E0*I_ERI_F3x_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F3x_Pz_Pz_S_b;
  abcd[821] = 4.0E0*I_ERI_F2xy_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Pz_S_b;
  abcd[822] = 4.0E0*I_ERI_F2xz_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Pz_S_b;
  abcd[823] = 4.0E0*I_ERI_Fx2y_Pz_F2yz_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Pz_S_b;
  abcd[824] = 4.0E0*I_ERI_Fxyz_Pz_F2yz_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Pz_S_b;
  abcd[825] = 4.0E0*I_ERI_Fx2z_Pz_F2yz_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Pz_S_b;
  abcd[826] = 4.0E0*I_ERI_F3y_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F3y_Pz_Pz_S_b;
  abcd[827] = 4.0E0*I_ERI_F2yz_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Pz_S_b;
  abcd[828] = 4.0E0*I_ERI_Fy2z_Pz_F2yz_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Pz_S_b;
  abcd[829] = 4.0E0*I_ERI_F3z_Pz_F2yz_S_bc-2.0E0*1*I_ERI_F3z_Pz_Pz_S_b;
  abcd[830] = 4.0E0*I_ERI_F3x_Pz_Fy2z_S_bc;
  abcd[831] = 4.0E0*I_ERI_F2xy_Pz_Fy2z_S_bc;
  abcd[832] = 4.0E0*I_ERI_F2xz_Pz_Fy2z_S_bc;
  abcd[833] = 4.0E0*I_ERI_Fx2y_Pz_Fy2z_S_bc;
  abcd[834] = 4.0E0*I_ERI_Fxyz_Pz_Fy2z_S_bc;
  abcd[835] = 4.0E0*I_ERI_Fx2z_Pz_Fy2z_S_bc;
  abcd[836] = 4.0E0*I_ERI_F3y_Pz_Fy2z_S_bc;
  abcd[837] = 4.0E0*I_ERI_F2yz_Pz_Fy2z_S_bc;
  abcd[838] = 4.0E0*I_ERI_Fy2z_Pz_Fy2z_S_bc;
  abcd[839] = 4.0E0*I_ERI_F3z_Pz_Fy2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_b
   ************************************************************/
  abcd[840] = 4.0E0*I_ERI_F3x_Pz_F2xz_S_bc;
  abcd[841] = 4.0E0*I_ERI_F2xy_Pz_F2xz_S_bc;
  abcd[842] = 4.0E0*I_ERI_F2xz_Pz_F2xz_S_bc;
  abcd[843] = 4.0E0*I_ERI_Fx2y_Pz_F2xz_S_bc;
  abcd[844] = 4.0E0*I_ERI_Fxyz_Pz_F2xz_S_bc;
  abcd[845] = 4.0E0*I_ERI_Fx2z_Pz_F2xz_S_bc;
  abcd[846] = 4.0E0*I_ERI_F3y_Pz_F2xz_S_bc;
  abcd[847] = 4.0E0*I_ERI_F2yz_Pz_F2xz_S_bc;
  abcd[848] = 4.0E0*I_ERI_Fy2z_Pz_F2xz_S_bc;
  abcd[849] = 4.0E0*I_ERI_F3z_Pz_F2xz_S_bc;
  abcd[850] = 4.0E0*I_ERI_F3x_Pz_Fxyz_S_bc;
  abcd[851] = 4.0E0*I_ERI_F2xy_Pz_Fxyz_S_bc;
  abcd[852] = 4.0E0*I_ERI_F2xz_Pz_Fxyz_S_bc;
  abcd[853] = 4.0E0*I_ERI_Fx2y_Pz_Fxyz_S_bc;
  abcd[854] = 4.0E0*I_ERI_Fxyz_Pz_Fxyz_S_bc;
  abcd[855] = 4.0E0*I_ERI_Fx2z_Pz_Fxyz_S_bc;
  abcd[856] = 4.0E0*I_ERI_F3y_Pz_Fxyz_S_bc;
  abcd[857] = 4.0E0*I_ERI_F2yz_Pz_Fxyz_S_bc;
  abcd[858] = 4.0E0*I_ERI_Fy2z_Pz_Fxyz_S_bc;
  abcd[859] = 4.0E0*I_ERI_F3z_Pz_Fxyz_S_bc;
  abcd[860] = 4.0E0*I_ERI_F3x_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F3x_Pz_Px_S_b;
  abcd[861] = 4.0E0*I_ERI_F2xy_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Px_S_b;
  abcd[862] = 4.0E0*I_ERI_F2xz_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Px_S_b;
  abcd[863] = 4.0E0*I_ERI_Fx2y_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Px_S_b;
  abcd[864] = 4.0E0*I_ERI_Fxyz_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Px_S_b;
  abcd[865] = 4.0E0*I_ERI_Fx2z_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Px_S_b;
  abcd[866] = 4.0E0*I_ERI_F3y_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F3y_Pz_Px_S_b;
  abcd[867] = 4.0E0*I_ERI_F2yz_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Px_S_b;
  abcd[868] = 4.0E0*I_ERI_Fy2z_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Px_S_b;
  abcd[869] = 4.0E0*I_ERI_F3z_Pz_Fx2z_S_bc-2.0E0*1*I_ERI_F3z_Pz_Px_S_b;
  abcd[870] = 4.0E0*I_ERI_F3x_Pz_F2yz_S_bc;
  abcd[871] = 4.0E0*I_ERI_F2xy_Pz_F2yz_S_bc;
  abcd[872] = 4.0E0*I_ERI_F2xz_Pz_F2yz_S_bc;
  abcd[873] = 4.0E0*I_ERI_Fx2y_Pz_F2yz_S_bc;
  abcd[874] = 4.0E0*I_ERI_Fxyz_Pz_F2yz_S_bc;
  abcd[875] = 4.0E0*I_ERI_Fx2z_Pz_F2yz_S_bc;
  abcd[876] = 4.0E0*I_ERI_F3y_Pz_F2yz_S_bc;
  abcd[877] = 4.0E0*I_ERI_F2yz_Pz_F2yz_S_bc;
  abcd[878] = 4.0E0*I_ERI_Fy2z_Pz_F2yz_S_bc;
  abcd[879] = 4.0E0*I_ERI_F3z_Pz_F2yz_S_bc;
  abcd[880] = 4.0E0*I_ERI_F3x_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F3x_Pz_Py_S_b;
  abcd[881] = 4.0E0*I_ERI_F2xy_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Py_S_b;
  abcd[882] = 4.0E0*I_ERI_F2xz_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Py_S_b;
  abcd[883] = 4.0E0*I_ERI_Fx2y_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Py_S_b;
  abcd[884] = 4.0E0*I_ERI_Fxyz_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Py_S_b;
  abcd[885] = 4.0E0*I_ERI_Fx2z_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Py_S_b;
  abcd[886] = 4.0E0*I_ERI_F3y_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F3y_Pz_Py_S_b;
  abcd[887] = 4.0E0*I_ERI_F2yz_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Py_S_b;
  abcd[888] = 4.0E0*I_ERI_Fy2z_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Py_S_b;
  abcd[889] = 4.0E0*I_ERI_F3z_Pz_Fy2z_S_bc-2.0E0*1*I_ERI_F3z_Pz_Py_S_b;
  abcd[890] = 4.0E0*I_ERI_F3x_Pz_F3z_S_bc-2.0E0*2*I_ERI_F3x_Pz_Pz_S_b;
  abcd[891] = 4.0E0*I_ERI_F2xy_Pz_F3z_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Pz_S_b;
  abcd[892] = 4.0E0*I_ERI_F2xz_Pz_F3z_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Pz_S_b;
  abcd[893] = 4.0E0*I_ERI_Fx2y_Pz_F3z_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Pz_S_b;
  abcd[894] = 4.0E0*I_ERI_Fxyz_Pz_F3z_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Pz_S_b;
  abcd[895] = 4.0E0*I_ERI_Fx2z_Pz_F3z_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Pz_S_b;
  abcd[896] = 4.0E0*I_ERI_F3y_Pz_F3z_S_bc-2.0E0*2*I_ERI_F3y_Pz_Pz_S_b;
  abcd[897] = 4.0E0*I_ERI_F2yz_Pz_F3z_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Pz_S_b;
  abcd[898] = 4.0E0*I_ERI_Fy2z_Pz_F3z_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Pz_S_b;
  abcd[899] = 4.0E0*I_ERI_F3z_Pz_F3z_S_bc-2.0E0*2*I_ERI_F3z_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[900] = 4.0E0*I_ERI_F3x_Px_D2x_Px_bd;
  abcd[901] = 4.0E0*I_ERI_F2xy_Px_D2x_Px_bd;
  abcd[902] = 4.0E0*I_ERI_F2xz_Px_D2x_Px_bd;
  abcd[903] = 4.0E0*I_ERI_Fx2y_Px_D2x_Px_bd;
  abcd[904] = 4.0E0*I_ERI_Fxyz_Px_D2x_Px_bd;
  abcd[905] = 4.0E0*I_ERI_Fx2z_Px_D2x_Px_bd;
  abcd[906] = 4.0E0*I_ERI_F3y_Px_D2x_Px_bd;
  abcd[907] = 4.0E0*I_ERI_F2yz_Px_D2x_Px_bd;
  abcd[908] = 4.0E0*I_ERI_Fy2z_Px_D2x_Px_bd;
  abcd[909] = 4.0E0*I_ERI_F3z_Px_D2x_Px_bd;
  abcd[910] = 4.0E0*I_ERI_F3x_Px_Dxy_Px_bd;
  abcd[911] = 4.0E0*I_ERI_F2xy_Px_Dxy_Px_bd;
  abcd[912] = 4.0E0*I_ERI_F2xz_Px_Dxy_Px_bd;
  abcd[913] = 4.0E0*I_ERI_Fx2y_Px_Dxy_Px_bd;
  abcd[914] = 4.0E0*I_ERI_Fxyz_Px_Dxy_Px_bd;
  abcd[915] = 4.0E0*I_ERI_Fx2z_Px_Dxy_Px_bd;
  abcd[916] = 4.0E0*I_ERI_F3y_Px_Dxy_Px_bd;
  abcd[917] = 4.0E0*I_ERI_F2yz_Px_Dxy_Px_bd;
  abcd[918] = 4.0E0*I_ERI_Fy2z_Px_Dxy_Px_bd;
  abcd[919] = 4.0E0*I_ERI_F3z_Px_Dxy_Px_bd;
  abcd[920] = 4.0E0*I_ERI_F3x_Px_Dxz_Px_bd;
  abcd[921] = 4.0E0*I_ERI_F2xy_Px_Dxz_Px_bd;
  abcd[922] = 4.0E0*I_ERI_F2xz_Px_Dxz_Px_bd;
  abcd[923] = 4.0E0*I_ERI_Fx2y_Px_Dxz_Px_bd;
  abcd[924] = 4.0E0*I_ERI_Fxyz_Px_Dxz_Px_bd;
  abcd[925] = 4.0E0*I_ERI_Fx2z_Px_Dxz_Px_bd;
  abcd[926] = 4.0E0*I_ERI_F3y_Px_Dxz_Px_bd;
  abcd[927] = 4.0E0*I_ERI_F2yz_Px_Dxz_Px_bd;
  abcd[928] = 4.0E0*I_ERI_Fy2z_Px_Dxz_Px_bd;
  abcd[929] = 4.0E0*I_ERI_F3z_Px_Dxz_Px_bd;
  abcd[930] = 4.0E0*I_ERI_F3x_Px_D2y_Px_bd;
  abcd[931] = 4.0E0*I_ERI_F2xy_Px_D2y_Px_bd;
  abcd[932] = 4.0E0*I_ERI_F2xz_Px_D2y_Px_bd;
  abcd[933] = 4.0E0*I_ERI_Fx2y_Px_D2y_Px_bd;
  abcd[934] = 4.0E0*I_ERI_Fxyz_Px_D2y_Px_bd;
  abcd[935] = 4.0E0*I_ERI_Fx2z_Px_D2y_Px_bd;
  abcd[936] = 4.0E0*I_ERI_F3y_Px_D2y_Px_bd;
  abcd[937] = 4.0E0*I_ERI_F2yz_Px_D2y_Px_bd;
  abcd[938] = 4.0E0*I_ERI_Fy2z_Px_D2y_Px_bd;
  abcd[939] = 4.0E0*I_ERI_F3z_Px_D2y_Px_bd;
  abcd[940] = 4.0E0*I_ERI_F3x_Px_Dyz_Px_bd;
  abcd[941] = 4.0E0*I_ERI_F2xy_Px_Dyz_Px_bd;
  abcd[942] = 4.0E0*I_ERI_F2xz_Px_Dyz_Px_bd;
  abcd[943] = 4.0E0*I_ERI_Fx2y_Px_Dyz_Px_bd;
  abcd[944] = 4.0E0*I_ERI_Fxyz_Px_Dyz_Px_bd;
  abcd[945] = 4.0E0*I_ERI_Fx2z_Px_Dyz_Px_bd;
  abcd[946] = 4.0E0*I_ERI_F3y_Px_Dyz_Px_bd;
  abcd[947] = 4.0E0*I_ERI_F2yz_Px_Dyz_Px_bd;
  abcd[948] = 4.0E0*I_ERI_Fy2z_Px_Dyz_Px_bd;
  abcd[949] = 4.0E0*I_ERI_F3z_Px_Dyz_Px_bd;
  abcd[950] = 4.0E0*I_ERI_F3x_Px_D2z_Px_bd;
  abcd[951] = 4.0E0*I_ERI_F2xy_Px_D2z_Px_bd;
  abcd[952] = 4.0E0*I_ERI_F2xz_Px_D2z_Px_bd;
  abcd[953] = 4.0E0*I_ERI_Fx2y_Px_D2z_Px_bd;
  abcd[954] = 4.0E0*I_ERI_Fxyz_Px_D2z_Px_bd;
  abcd[955] = 4.0E0*I_ERI_Fx2z_Px_D2z_Px_bd;
  abcd[956] = 4.0E0*I_ERI_F3y_Px_D2z_Px_bd;
  abcd[957] = 4.0E0*I_ERI_F2yz_Px_D2z_Px_bd;
  abcd[958] = 4.0E0*I_ERI_Fy2z_Px_D2z_Px_bd;
  abcd[959] = 4.0E0*I_ERI_F3z_Px_D2z_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[960] = 4.0E0*I_ERI_F3x_Px_D2x_Py_bd;
  abcd[961] = 4.0E0*I_ERI_F2xy_Px_D2x_Py_bd;
  abcd[962] = 4.0E0*I_ERI_F2xz_Px_D2x_Py_bd;
  abcd[963] = 4.0E0*I_ERI_Fx2y_Px_D2x_Py_bd;
  abcd[964] = 4.0E0*I_ERI_Fxyz_Px_D2x_Py_bd;
  abcd[965] = 4.0E0*I_ERI_Fx2z_Px_D2x_Py_bd;
  abcd[966] = 4.0E0*I_ERI_F3y_Px_D2x_Py_bd;
  abcd[967] = 4.0E0*I_ERI_F2yz_Px_D2x_Py_bd;
  abcd[968] = 4.0E0*I_ERI_Fy2z_Px_D2x_Py_bd;
  abcd[969] = 4.0E0*I_ERI_F3z_Px_D2x_Py_bd;
  abcd[970] = 4.0E0*I_ERI_F3x_Px_Dxy_Py_bd;
  abcd[971] = 4.0E0*I_ERI_F2xy_Px_Dxy_Py_bd;
  abcd[972] = 4.0E0*I_ERI_F2xz_Px_Dxy_Py_bd;
  abcd[973] = 4.0E0*I_ERI_Fx2y_Px_Dxy_Py_bd;
  abcd[974] = 4.0E0*I_ERI_Fxyz_Px_Dxy_Py_bd;
  abcd[975] = 4.0E0*I_ERI_Fx2z_Px_Dxy_Py_bd;
  abcd[976] = 4.0E0*I_ERI_F3y_Px_Dxy_Py_bd;
  abcd[977] = 4.0E0*I_ERI_F2yz_Px_Dxy_Py_bd;
  abcd[978] = 4.0E0*I_ERI_Fy2z_Px_Dxy_Py_bd;
  abcd[979] = 4.0E0*I_ERI_F3z_Px_Dxy_Py_bd;
  abcd[980] = 4.0E0*I_ERI_F3x_Px_Dxz_Py_bd;
  abcd[981] = 4.0E0*I_ERI_F2xy_Px_Dxz_Py_bd;
  abcd[982] = 4.0E0*I_ERI_F2xz_Px_Dxz_Py_bd;
  abcd[983] = 4.0E0*I_ERI_Fx2y_Px_Dxz_Py_bd;
  abcd[984] = 4.0E0*I_ERI_Fxyz_Px_Dxz_Py_bd;
  abcd[985] = 4.0E0*I_ERI_Fx2z_Px_Dxz_Py_bd;
  abcd[986] = 4.0E0*I_ERI_F3y_Px_Dxz_Py_bd;
  abcd[987] = 4.0E0*I_ERI_F2yz_Px_Dxz_Py_bd;
  abcd[988] = 4.0E0*I_ERI_Fy2z_Px_Dxz_Py_bd;
  abcd[989] = 4.0E0*I_ERI_F3z_Px_Dxz_Py_bd;
  abcd[990] = 4.0E0*I_ERI_F3x_Px_D2y_Py_bd;
  abcd[991] = 4.0E0*I_ERI_F2xy_Px_D2y_Py_bd;
  abcd[992] = 4.0E0*I_ERI_F2xz_Px_D2y_Py_bd;
  abcd[993] = 4.0E0*I_ERI_Fx2y_Px_D2y_Py_bd;
  abcd[994] = 4.0E0*I_ERI_Fxyz_Px_D2y_Py_bd;
  abcd[995] = 4.0E0*I_ERI_Fx2z_Px_D2y_Py_bd;
  abcd[996] = 4.0E0*I_ERI_F3y_Px_D2y_Py_bd;
  abcd[997] = 4.0E0*I_ERI_F2yz_Px_D2y_Py_bd;
  abcd[998] = 4.0E0*I_ERI_Fy2z_Px_D2y_Py_bd;
  abcd[999] = 4.0E0*I_ERI_F3z_Px_D2y_Py_bd;
  abcd[1000] = 4.0E0*I_ERI_F3x_Px_Dyz_Py_bd;
  abcd[1001] = 4.0E0*I_ERI_F2xy_Px_Dyz_Py_bd;
  abcd[1002] = 4.0E0*I_ERI_F2xz_Px_Dyz_Py_bd;
  abcd[1003] = 4.0E0*I_ERI_Fx2y_Px_Dyz_Py_bd;
  abcd[1004] = 4.0E0*I_ERI_Fxyz_Px_Dyz_Py_bd;
  abcd[1005] = 4.0E0*I_ERI_Fx2z_Px_Dyz_Py_bd;
  abcd[1006] = 4.0E0*I_ERI_F3y_Px_Dyz_Py_bd;
  abcd[1007] = 4.0E0*I_ERI_F2yz_Px_Dyz_Py_bd;
  abcd[1008] = 4.0E0*I_ERI_Fy2z_Px_Dyz_Py_bd;
  abcd[1009] = 4.0E0*I_ERI_F3z_Px_Dyz_Py_bd;
  abcd[1010] = 4.0E0*I_ERI_F3x_Px_D2z_Py_bd;
  abcd[1011] = 4.0E0*I_ERI_F2xy_Px_D2z_Py_bd;
  abcd[1012] = 4.0E0*I_ERI_F2xz_Px_D2z_Py_bd;
  abcd[1013] = 4.0E0*I_ERI_Fx2y_Px_D2z_Py_bd;
  abcd[1014] = 4.0E0*I_ERI_Fxyz_Px_D2z_Py_bd;
  abcd[1015] = 4.0E0*I_ERI_Fx2z_Px_D2z_Py_bd;
  abcd[1016] = 4.0E0*I_ERI_F3y_Px_D2z_Py_bd;
  abcd[1017] = 4.0E0*I_ERI_F2yz_Px_D2z_Py_bd;
  abcd[1018] = 4.0E0*I_ERI_Fy2z_Px_D2z_Py_bd;
  abcd[1019] = 4.0E0*I_ERI_F3z_Px_D2z_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1020] = 4.0E0*I_ERI_F3x_Px_D2x_Pz_bd;
  abcd[1021] = 4.0E0*I_ERI_F2xy_Px_D2x_Pz_bd;
  abcd[1022] = 4.0E0*I_ERI_F2xz_Px_D2x_Pz_bd;
  abcd[1023] = 4.0E0*I_ERI_Fx2y_Px_D2x_Pz_bd;
  abcd[1024] = 4.0E0*I_ERI_Fxyz_Px_D2x_Pz_bd;
  abcd[1025] = 4.0E0*I_ERI_Fx2z_Px_D2x_Pz_bd;
  abcd[1026] = 4.0E0*I_ERI_F3y_Px_D2x_Pz_bd;
  abcd[1027] = 4.0E0*I_ERI_F2yz_Px_D2x_Pz_bd;
  abcd[1028] = 4.0E0*I_ERI_Fy2z_Px_D2x_Pz_bd;
  abcd[1029] = 4.0E0*I_ERI_F3z_Px_D2x_Pz_bd;
  abcd[1030] = 4.0E0*I_ERI_F3x_Px_Dxy_Pz_bd;
  abcd[1031] = 4.0E0*I_ERI_F2xy_Px_Dxy_Pz_bd;
  abcd[1032] = 4.0E0*I_ERI_F2xz_Px_Dxy_Pz_bd;
  abcd[1033] = 4.0E0*I_ERI_Fx2y_Px_Dxy_Pz_bd;
  abcd[1034] = 4.0E0*I_ERI_Fxyz_Px_Dxy_Pz_bd;
  abcd[1035] = 4.0E0*I_ERI_Fx2z_Px_Dxy_Pz_bd;
  abcd[1036] = 4.0E0*I_ERI_F3y_Px_Dxy_Pz_bd;
  abcd[1037] = 4.0E0*I_ERI_F2yz_Px_Dxy_Pz_bd;
  abcd[1038] = 4.0E0*I_ERI_Fy2z_Px_Dxy_Pz_bd;
  abcd[1039] = 4.0E0*I_ERI_F3z_Px_Dxy_Pz_bd;
  abcd[1040] = 4.0E0*I_ERI_F3x_Px_Dxz_Pz_bd;
  abcd[1041] = 4.0E0*I_ERI_F2xy_Px_Dxz_Pz_bd;
  abcd[1042] = 4.0E0*I_ERI_F2xz_Px_Dxz_Pz_bd;
  abcd[1043] = 4.0E0*I_ERI_Fx2y_Px_Dxz_Pz_bd;
  abcd[1044] = 4.0E0*I_ERI_Fxyz_Px_Dxz_Pz_bd;
  abcd[1045] = 4.0E0*I_ERI_Fx2z_Px_Dxz_Pz_bd;
  abcd[1046] = 4.0E0*I_ERI_F3y_Px_Dxz_Pz_bd;
  abcd[1047] = 4.0E0*I_ERI_F2yz_Px_Dxz_Pz_bd;
  abcd[1048] = 4.0E0*I_ERI_Fy2z_Px_Dxz_Pz_bd;
  abcd[1049] = 4.0E0*I_ERI_F3z_Px_Dxz_Pz_bd;
  abcd[1050] = 4.0E0*I_ERI_F3x_Px_D2y_Pz_bd;
  abcd[1051] = 4.0E0*I_ERI_F2xy_Px_D2y_Pz_bd;
  abcd[1052] = 4.0E0*I_ERI_F2xz_Px_D2y_Pz_bd;
  abcd[1053] = 4.0E0*I_ERI_Fx2y_Px_D2y_Pz_bd;
  abcd[1054] = 4.0E0*I_ERI_Fxyz_Px_D2y_Pz_bd;
  abcd[1055] = 4.0E0*I_ERI_Fx2z_Px_D2y_Pz_bd;
  abcd[1056] = 4.0E0*I_ERI_F3y_Px_D2y_Pz_bd;
  abcd[1057] = 4.0E0*I_ERI_F2yz_Px_D2y_Pz_bd;
  abcd[1058] = 4.0E0*I_ERI_Fy2z_Px_D2y_Pz_bd;
  abcd[1059] = 4.0E0*I_ERI_F3z_Px_D2y_Pz_bd;
  abcd[1060] = 4.0E0*I_ERI_F3x_Px_Dyz_Pz_bd;
  abcd[1061] = 4.0E0*I_ERI_F2xy_Px_Dyz_Pz_bd;
  abcd[1062] = 4.0E0*I_ERI_F2xz_Px_Dyz_Pz_bd;
  abcd[1063] = 4.0E0*I_ERI_Fx2y_Px_Dyz_Pz_bd;
  abcd[1064] = 4.0E0*I_ERI_Fxyz_Px_Dyz_Pz_bd;
  abcd[1065] = 4.0E0*I_ERI_Fx2z_Px_Dyz_Pz_bd;
  abcd[1066] = 4.0E0*I_ERI_F3y_Px_Dyz_Pz_bd;
  abcd[1067] = 4.0E0*I_ERI_F2yz_Px_Dyz_Pz_bd;
  abcd[1068] = 4.0E0*I_ERI_Fy2z_Px_Dyz_Pz_bd;
  abcd[1069] = 4.0E0*I_ERI_F3z_Px_Dyz_Pz_bd;
  abcd[1070] = 4.0E0*I_ERI_F3x_Px_D2z_Pz_bd;
  abcd[1071] = 4.0E0*I_ERI_F2xy_Px_D2z_Pz_bd;
  abcd[1072] = 4.0E0*I_ERI_F2xz_Px_D2z_Pz_bd;
  abcd[1073] = 4.0E0*I_ERI_Fx2y_Px_D2z_Pz_bd;
  abcd[1074] = 4.0E0*I_ERI_Fxyz_Px_D2z_Pz_bd;
  abcd[1075] = 4.0E0*I_ERI_Fx2z_Px_D2z_Pz_bd;
  abcd[1076] = 4.0E0*I_ERI_F3y_Px_D2z_Pz_bd;
  abcd[1077] = 4.0E0*I_ERI_F2yz_Px_D2z_Pz_bd;
  abcd[1078] = 4.0E0*I_ERI_Fy2z_Px_D2z_Pz_bd;
  abcd[1079] = 4.0E0*I_ERI_F3z_Px_D2z_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1080] = 4.0E0*I_ERI_F3x_Py_D2x_Px_bd;
  abcd[1081] = 4.0E0*I_ERI_F2xy_Py_D2x_Px_bd;
  abcd[1082] = 4.0E0*I_ERI_F2xz_Py_D2x_Px_bd;
  abcd[1083] = 4.0E0*I_ERI_Fx2y_Py_D2x_Px_bd;
  abcd[1084] = 4.0E0*I_ERI_Fxyz_Py_D2x_Px_bd;
  abcd[1085] = 4.0E0*I_ERI_Fx2z_Py_D2x_Px_bd;
  abcd[1086] = 4.0E0*I_ERI_F3y_Py_D2x_Px_bd;
  abcd[1087] = 4.0E0*I_ERI_F2yz_Py_D2x_Px_bd;
  abcd[1088] = 4.0E0*I_ERI_Fy2z_Py_D2x_Px_bd;
  abcd[1089] = 4.0E0*I_ERI_F3z_Py_D2x_Px_bd;
  abcd[1090] = 4.0E0*I_ERI_F3x_Py_Dxy_Px_bd;
  abcd[1091] = 4.0E0*I_ERI_F2xy_Py_Dxy_Px_bd;
  abcd[1092] = 4.0E0*I_ERI_F2xz_Py_Dxy_Px_bd;
  abcd[1093] = 4.0E0*I_ERI_Fx2y_Py_Dxy_Px_bd;
  abcd[1094] = 4.0E0*I_ERI_Fxyz_Py_Dxy_Px_bd;
  abcd[1095] = 4.0E0*I_ERI_Fx2z_Py_Dxy_Px_bd;
  abcd[1096] = 4.0E0*I_ERI_F3y_Py_Dxy_Px_bd;
  abcd[1097] = 4.0E0*I_ERI_F2yz_Py_Dxy_Px_bd;
  abcd[1098] = 4.0E0*I_ERI_Fy2z_Py_Dxy_Px_bd;
  abcd[1099] = 4.0E0*I_ERI_F3z_Py_Dxy_Px_bd;
  abcd[1100] = 4.0E0*I_ERI_F3x_Py_Dxz_Px_bd;
  abcd[1101] = 4.0E0*I_ERI_F2xy_Py_Dxz_Px_bd;
  abcd[1102] = 4.0E0*I_ERI_F2xz_Py_Dxz_Px_bd;
  abcd[1103] = 4.0E0*I_ERI_Fx2y_Py_Dxz_Px_bd;
  abcd[1104] = 4.0E0*I_ERI_Fxyz_Py_Dxz_Px_bd;
  abcd[1105] = 4.0E0*I_ERI_Fx2z_Py_Dxz_Px_bd;
  abcd[1106] = 4.0E0*I_ERI_F3y_Py_Dxz_Px_bd;
  abcd[1107] = 4.0E0*I_ERI_F2yz_Py_Dxz_Px_bd;
  abcd[1108] = 4.0E0*I_ERI_Fy2z_Py_Dxz_Px_bd;
  abcd[1109] = 4.0E0*I_ERI_F3z_Py_Dxz_Px_bd;
  abcd[1110] = 4.0E0*I_ERI_F3x_Py_D2y_Px_bd;
  abcd[1111] = 4.0E0*I_ERI_F2xy_Py_D2y_Px_bd;
  abcd[1112] = 4.0E0*I_ERI_F2xz_Py_D2y_Px_bd;
  abcd[1113] = 4.0E0*I_ERI_Fx2y_Py_D2y_Px_bd;
  abcd[1114] = 4.0E0*I_ERI_Fxyz_Py_D2y_Px_bd;
  abcd[1115] = 4.0E0*I_ERI_Fx2z_Py_D2y_Px_bd;
  abcd[1116] = 4.0E0*I_ERI_F3y_Py_D2y_Px_bd;
  abcd[1117] = 4.0E0*I_ERI_F2yz_Py_D2y_Px_bd;
  abcd[1118] = 4.0E0*I_ERI_Fy2z_Py_D2y_Px_bd;
  abcd[1119] = 4.0E0*I_ERI_F3z_Py_D2y_Px_bd;
  abcd[1120] = 4.0E0*I_ERI_F3x_Py_Dyz_Px_bd;
  abcd[1121] = 4.0E0*I_ERI_F2xy_Py_Dyz_Px_bd;
  abcd[1122] = 4.0E0*I_ERI_F2xz_Py_Dyz_Px_bd;
  abcd[1123] = 4.0E0*I_ERI_Fx2y_Py_Dyz_Px_bd;
  abcd[1124] = 4.0E0*I_ERI_Fxyz_Py_Dyz_Px_bd;
  abcd[1125] = 4.0E0*I_ERI_Fx2z_Py_Dyz_Px_bd;
  abcd[1126] = 4.0E0*I_ERI_F3y_Py_Dyz_Px_bd;
  abcd[1127] = 4.0E0*I_ERI_F2yz_Py_Dyz_Px_bd;
  abcd[1128] = 4.0E0*I_ERI_Fy2z_Py_Dyz_Px_bd;
  abcd[1129] = 4.0E0*I_ERI_F3z_Py_Dyz_Px_bd;
  abcd[1130] = 4.0E0*I_ERI_F3x_Py_D2z_Px_bd;
  abcd[1131] = 4.0E0*I_ERI_F2xy_Py_D2z_Px_bd;
  abcd[1132] = 4.0E0*I_ERI_F2xz_Py_D2z_Px_bd;
  abcd[1133] = 4.0E0*I_ERI_Fx2y_Py_D2z_Px_bd;
  abcd[1134] = 4.0E0*I_ERI_Fxyz_Py_D2z_Px_bd;
  abcd[1135] = 4.0E0*I_ERI_Fx2z_Py_D2z_Px_bd;
  abcd[1136] = 4.0E0*I_ERI_F3y_Py_D2z_Px_bd;
  abcd[1137] = 4.0E0*I_ERI_F2yz_Py_D2z_Px_bd;
  abcd[1138] = 4.0E0*I_ERI_Fy2z_Py_D2z_Px_bd;
  abcd[1139] = 4.0E0*I_ERI_F3z_Py_D2z_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1140] = 4.0E0*I_ERI_F3x_Py_D2x_Py_bd;
  abcd[1141] = 4.0E0*I_ERI_F2xy_Py_D2x_Py_bd;
  abcd[1142] = 4.0E0*I_ERI_F2xz_Py_D2x_Py_bd;
  abcd[1143] = 4.0E0*I_ERI_Fx2y_Py_D2x_Py_bd;
  abcd[1144] = 4.0E0*I_ERI_Fxyz_Py_D2x_Py_bd;
  abcd[1145] = 4.0E0*I_ERI_Fx2z_Py_D2x_Py_bd;
  abcd[1146] = 4.0E0*I_ERI_F3y_Py_D2x_Py_bd;
  abcd[1147] = 4.0E0*I_ERI_F2yz_Py_D2x_Py_bd;
  abcd[1148] = 4.0E0*I_ERI_Fy2z_Py_D2x_Py_bd;
  abcd[1149] = 4.0E0*I_ERI_F3z_Py_D2x_Py_bd;
  abcd[1150] = 4.0E0*I_ERI_F3x_Py_Dxy_Py_bd;
  abcd[1151] = 4.0E0*I_ERI_F2xy_Py_Dxy_Py_bd;
  abcd[1152] = 4.0E0*I_ERI_F2xz_Py_Dxy_Py_bd;
  abcd[1153] = 4.0E0*I_ERI_Fx2y_Py_Dxy_Py_bd;
  abcd[1154] = 4.0E0*I_ERI_Fxyz_Py_Dxy_Py_bd;
  abcd[1155] = 4.0E0*I_ERI_Fx2z_Py_Dxy_Py_bd;
  abcd[1156] = 4.0E0*I_ERI_F3y_Py_Dxy_Py_bd;
  abcd[1157] = 4.0E0*I_ERI_F2yz_Py_Dxy_Py_bd;
  abcd[1158] = 4.0E0*I_ERI_Fy2z_Py_Dxy_Py_bd;
  abcd[1159] = 4.0E0*I_ERI_F3z_Py_Dxy_Py_bd;
  abcd[1160] = 4.0E0*I_ERI_F3x_Py_Dxz_Py_bd;
  abcd[1161] = 4.0E0*I_ERI_F2xy_Py_Dxz_Py_bd;
  abcd[1162] = 4.0E0*I_ERI_F2xz_Py_Dxz_Py_bd;
  abcd[1163] = 4.0E0*I_ERI_Fx2y_Py_Dxz_Py_bd;
  abcd[1164] = 4.0E0*I_ERI_Fxyz_Py_Dxz_Py_bd;
  abcd[1165] = 4.0E0*I_ERI_Fx2z_Py_Dxz_Py_bd;
  abcd[1166] = 4.0E0*I_ERI_F3y_Py_Dxz_Py_bd;
  abcd[1167] = 4.0E0*I_ERI_F2yz_Py_Dxz_Py_bd;
  abcd[1168] = 4.0E0*I_ERI_Fy2z_Py_Dxz_Py_bd;
  abcd[1169] = 4.0E0*I_ERI_F3z_Py_Dxz_Py_bd;
  abcd[1170] = 4.0E0*I_ERI_F3x_Py_D2y_Py_bd;
  abcd[1171] = 4.0E0*I_ERI_F2xy_Py_D2y_Py_bd;
  abcd[1172] = 4.0E0*I_ERI_F2xz_Py_D2y_Py_bd;
  abcd[1173] = 4.0E0*I_ERI_Fx2y_Py_D2y_Py_bd;
  abcd[1174] = 4.0E0*I_ERI_Fxyz_Py_D2y_Py_bd;
  abcd[1175] = 4.0E0*I_ERI_Fx2z_Py_D2y_Py_bd;
  abcd[1176] = 4.0E0*I_ERI_F3y_Py_D2y_Py_bd;
  abcd[1177] = 4.0E0*I_ERI_F2yz_Py_D2y_Py_bd;
  abcd[1178] = 4.0E0*I_ERI_Fy2z_Py_D2y_Py_bd;
  abcd[1179] = 4.0E0*I_ERI_F3z_Py_D2y_Py_bd;
  abcd[1180] = 4.0E0*I_ERI_F3x_Py_Dyz_Py_bd;
  abcd[1181] = 4.0E0*I_ERI_F2xy_Py_Dyz_Py_bd;
  abcd[1182] = 4.0E0*I_ERI_F2xz_Py_Dyz_Py_bd;
  abcd[1183] = 4.0E0*I_ERI_Fx2y_Py_Dyz_Py_bd;
  abcd[1184] = 4.0E0*I_ERI_Fxyz_Py_Dyz_Py_bd;
  abcd[1185] = 4.0E0*I_ERI_Fx2z_Py_Dyz_Py_bd;
  abcd[1186] = 4.0E0*I_ERI_F3y_Py_Dyz_Py_bd;
  abcd[1187] = 4.0E0*I_ERI_F2yz_Py_Dyz_Py_bd;
  abcd[1188] = 4.0E0*I_ERI_Fy2z_Py_Dyz_Py_bd;
  abcd[1189] = 4.0E0*I_ERI_F3z_Py_Dyz_Py_bd;
  abcd[1190] = 4.0E0*I_ERI_F3x_Py_D2z_Py_bd;
  abcd[1191] = 4.0E0*I_ERI_F2xy_Py_D2z_Py_bd;
  abcd[1192] = 4.0E0*I_ERI_F2xz_Py_D2z_Py_bd;
  abcd[1193] = 4.0E0*I_ERI_Fx2y_Py_D2z_Py_bd;
  abcd[1194] = 4.0E0*I_ERI_Fxyz_Py_D2z_Py_bd;
  abcd[1195] = 4.0E0*I_ERI_Fx2z_Py_D2z_Py_bd;
  abcd[1196] = 4.0E0*I_ERI_F3y_Py_D2z_Py_bd;
  abcd[1197] = 4.0E0*I_ERI_F2yz_Py_D2z_Py_bd;
  abcd[1198] = 4.0E0*I_ERI_Fy2z_Py_D2z_Py_bd;
  abcd[1199] = 4.0E0*I_ERI_F3z_Py_D2z_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dby_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1200] = 4.0E0*I_ERI_F3x_Py_D2x_Pz_bd;
  abcd[1201] = 4.0E0*I_ERI_F2xy_Py_D2x_Pz_bd;
  abcd[1202] = 4.0E0*I_ERI_F2xz_Py_D2x_Pz_bd;
  abcd[1203] = 4.0E0*I_ERI_Fx2y_Py_D2x_Pz_bd;
  abcd[1204] = 4.0E0*I_ERI_Fxyz_Py_D2x_Pz_bd;
  abcd[1205] = 4.0E0*I_ERI_Fx2z_Py_D2x_Pz_bd;
  abcd[1206] = 4.0E0*I_ERI_F3y_Py_D2x_Pz_bd;
  abcd[1207] = 4.0E0*I_ERI_F2yz_Py_D2x_Pz_bd;
  abcd[1208] = 4.0E0*I_ERI_Fy2z_Py_D2x_Pz_bd;
  abcd[1209] = 4.0E0*I_ERI_F3z_Py_D2x_Pz_bd;
  abcd[1210] = 4.0E0*I_ERI_F3x_Py_Dxy_Pz_bd;
  abcd[1211] = 4.0E0*I_ERI_F2xy_Py_Dxy_Pz_bd;
  abcd[1212] = 4.0E0*I_ERI_F2xz_Py_Dxy_Pz_bd;
  abcd[1213] = 4.0E0*I_ERI_Fx2y_Py_Dxy_Pz_bd;
  abcd[1214] = 4.0E0*I_ERI_Fxyz_Py_Dxy_Pz_bd;
  abcd[1215] = 4.0E0*I_ERI_Fx2z_Py_Dxy_Pz_bd;
  abcd[1216] = 4.0E0*I_ERI_F3y_Py_Dxy_Pz_bd;
  abcd[1217] = 4.0E0*I_ERI_F2yz_Py_Dxy_Pz_bd;
  abcd[1218] = 4.0E0*I_ERI_Fy2z_Py_Dxy_Pz_bd;
  abcd[1219] = 4.0E0*I_ERI_F3z_Py_Dxy_Pz_bd;
  abcd[1220] = 4.0E0*I_ERI_F3x_Py_Dxz_Pz_bd;
  abcd[1221] = 4.0E0*I_ERI_F2xy_Py_Dxz_Pz_bd;
  abcd[1222] = 4.0E0*I_ERI_F2xz_Py_Dxz_Pz_bd;
  abcd[1223] = 4.0E0*I_ERI_Fx2y_Py_Dxz_Pz_bd;
  abcd[1224] = 4.0E0*I_ERI_Fxyz_Py_Dxz_Pz_bd;
  abcd[1225] = 4.0E0*I_ERI_Fx2z_Py_Dxz_Pz_bd;
  abcd[1226] = 4.0E0*I_ERI_F3y_Py_Dxz_Pz_bd;
  abcd[1227] = 4.0E0*I_ERI_F2yz_Py_Dxz_Pz_bd;
  abcd[1228] = 4.0E0*I_ERI_Fy2z_Py_Dxz_Pz_bd;
  abcd[1229] = 4.0E0*I_ERI_F3z_Py_Dxz_Pz_bd;
  abcd[1230] = 4.0E0*I_ERI_F3x_Py_D2y_Pz_bd;
  abcd[1231] = 4.0E0*I_ERI_F2xy_Py_D2y_Pz_bd;
  abcd[1232] = 4.0E0*I_ERI_F2xz_Py_D2y_Pz_bd;
  abcd[1233] = 4.0E0*I_ERI_Fx2y_Py_D2y_Pz_bd;
  abcd[1234] = 4.0E0*I_ERI_Fxyz_Py_D2y_Pz_bd;
  abcd[1235] = 4.0E0*I_ERI_Fx2z_Py_D2y_Pz_bd;
  abcd[1236] = 4.0E0*I_ERI_F3y_Py_D2y_Pz_bd;
  abcd[1237] = 4.0E0*I_ERI_F2yz_Py_D2y_Pz_bd;
  abcd[1238] = 4.0E0*I_ERI_Fy2z_Py_D2y_Pz_bd;
  abcd[1239] = 4.0E0*I_ERI_F3z_Py_D2y_Pz_bd;
  abcd[1240] = 4.0E0*I_ERI_F3x_Py_Dyz_Pz_bd;
  abcd[1241] = 4.0E0*I_ERI_F2xy_Py_Dyz_Pz_bd;
  abcd[1242] = 4.0E0*I_ERI_F2xz_Py_Dyz_Pz_bd;
  abcd[1243] = 4.0E0*I_ERI_Fx2y_Py_Dyz_Pz_bd;
  abcd[1244] = 4.0E0*I_ERI_Fxyz_Py_Dyz_Pz_bd;
  abcd[1245] = 4.0E0*I_ERI_Fx2z_Py_Dyz_Pz_bd;
  abcd[1246] = 4.0E0*I_ERI_F3y_Py_Dyz_Pz_bd;
  abcd[1247] = 4.0E0*I_ERI_F2yz_Py_Dyz_Pz_bd;
  abcd[1248] = 4.0E0*I_ERI_Fy2z_Py_Dyz_Pz_bd;
  abcd[1249] = 4.0E0*I_ERI_F3z_Py_Dyz_Pz_bd;
  abcd[1250] = 4.0E0*I_ERI_F3x_Py_D2z_Pz_bd;
  abcd[1251] = 4.0E0*I_ERI_F2xy_Py_D2z_Pz_bd;
  abcd[1252] = 4.0E0*I_ERI_F2xz_Py_D2z_Pz_bd;
  abcd[1253] = 4.0E0*I_ERI_Fx2y_Py_D2z_Pz_bd;
  abcd[1254] = 4.0E0*I_ERI_Fxyz_Py_D2z_Pz_bd;
  abcd[1255] = 4.0E0*I_ERI_Fx2z_Py_D2z_Pz_bd;
  abcd[1256] = 4.0E0*I_ERI_F3y_Py_D2z_Pz_bd;
  abcd[1257] = 4.0E0*I_ERI_F2yz_Py_D2z_Pz_bd;
  abcd[1258] = 4.0E0*I_ERI_Fy2z_Py_D2z_Pz_bd;
  abcd[1259] = 4.0E0*I_ERI_F3z_Py_D2z_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1260] = 4.0E0*I_ERI_F3x_Pz_D2x_Px_bd;
  abcd[1261] = 4.0E0*I_ERI_F2xy_Pz_D2x_Px_bd;
  abcd[1262] = 4.0E0*I_ERI_F2xz_Pz_D2x_Px_bd;
  abcd[1263] = 4.0E0*I_ERI_Fx2y_Pz_D2x_Px_bd;
  abcd[1264] = 4.0E0*I_ERI_Fxyz_Pz_D2x_Px_bd;
  abcd[1265] = 4.0E0*I_ERI_Fx2z_Pz_D2x_Px_bd;
  abcd[1266] = 4.0E0*I_ERI_F3y_Pz_D2x_Px_bd;
  abcd[1267] = 4.0E0*I_ERI_F2yz_Pz_D2x_Px_bd;
  abcd[1268] = 4.0E0*I_ERI_Fy2z_Pz_D2x_Px_bd;
  abcd[1269] = 4.0E0*I_ERI_F3z_Pz_D2x_Px_bd;
  abcd[1270] = 4.0E0*I_ERI_F3x_Pz_Dxy_Px_bd;
  abcd[1271] = 4.0E0*I_ERI_F2xy_Pz_Dxy_Px_bd;
  abcd[1272] = 4.0E0*I_ERI_F2xz_Pz_Dxy_Px_bd;
  abcd[1273] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_Px_bd;
  abcd[1274] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_Px_bd;
  abcd[1275] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_Px_bd;
  abcd[1276] = 4.0E0*I_ERI_F3y_Pz_Dxy_Px_bd;
  abcd[1277] = 4.0E0*I_ERI_F2yz_Pz_Dxy_Px_bd;
  abcd[1278] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_Px_bd;
  abcd[1279] = 4.0E0*I_ERI_F3z_Pz_Dxy_Px_bd;
  abcd[1280] = 4.0E0*I_ERI_F3x_Pz_Dxz_Px_bd;
  abcd[1281] = 4.0E0*I_ERI_F2xy_Pz_Dxz_Px_bd;
  abcd[1282] = 4.0E0*I_ERI_F2xz_Pz_Dxz_Px_bd;
  abcd[1283] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_Px_bd;
  abcd[1284] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_Px_bd;
  abcd[1285] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_Px_bd;
  abcd[1286] = 4.0E0*I_ERI_F3y_Pz_Dxz_Px_bd;
  abcd[1287] = 4.0E0*I_ERI_F2yz_Pz_Dxz_Px_bd;
  abcd[1288] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_Px_bd;
  abcd[1289] = 4.0E0*I_ERI_F3z_Pz_Dxz_Px_bd;
  abcd[1290] = 4.0E0*I_ERI_F3x_Pz_D2y_Px_bd;
  abcd[1291] = 4.0E0*I_ERI_F2xy_Pz_D2y_Px_bd;
  abcd[1292] = 4.0E0*I_ERI_F2xz_Pz_D2y_Px_bd;
  abcd[1293] = 4.0E0*I_ERI_Fx2y_Pz_D2y_Px_bd;
  abcd[1294] = 4.0E0*I_ERI_Fxyz_Pz_D2y_Px_bd;
  abcd[1295] = 4.0E0*I_ERI_Fx2z_Pz_D2y_Px_bd;
  abcd[1296] = 4.0E0*I_ERI_F3y_Pz_D2y_Px_bd;
  abcd[1297] = 4.0E0*I_ERI_F2yz_Pz_D2y_Px_bd;
  abcd[1298] = 4.0E0*I_ERI_Fy2z_Pz_D2y_Px_bd;
  abcd[1299] = 4.0E0*I_ERI_F3z_Pz_D2y_Px_bd;
  abcd[1300] = 4.0E0*I_ERI_F3x_Pz_Dyz_Px_bd;
  abcd[1301] = 4.0E0*I_ERI_F2xy_Pz_Dyz_Px_bd;
  abcd[1302] = 4.0E0*I_ERI_F2xz_Pz_Dyz_Px_bd;
  abcd[1303] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_Px_bd;
  abcd[1304] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_Px_bd;
  abcd[1305] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_Px_bd;
  abcd[1306] = 4.0E0*I_ERI_F3y_Pz_Dyz_Px_bd;
  abcd[1307] = 4.0E0*I_ERI_F2yz_Pz_Dyz_Px_bd;
  abcd[1308] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_Px_bd;
  abcd[1309] = 4.0E0*I_ERI_F3z_Pz_Dyz_Px_bd;
  abcd[1310] = 4.0E0*I_ERI_F3x_Pz_D2z_Px_bd;
  abcd[1311] = 4.0E0*I_ERI_F2xy_Pz_D2z_Px_bd;
  abcd[1312] = 4.0E0*I_ERI_F2xz_Pz_D2z_Px_bd;
  abcd[1313] = 4.0E0*I_ERI_Fx2y_Pz_D2z_Px_bd;
  abcd[1314] = 4.0E0*I_ERI_Fxyz_Pz_D2z_Px_bd;
  abcd[1315] = 4.0E0*I_ERI_Fx2z_Pz_D2z_Px_bd;
  abcd[1316] = 4.0E0*I_ERI_F3y_Pz_D2z_Px_bd;
  abcd[1317] = 4.0E0*I_ERI_F2yz_Pz_D2z_Px_bd;
  abcd[1318] = 4.0E0*I_ERI_Fy2z_Pz_D2z_Px_bd;
  abcd[1319] = 4.0E0*I_ERI_F3z_Pz_D2z_Px_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1320] = 4.0E0*I_ERI_F3x_Pz_D2x_Py_bd;
  abcd[1321] = 4.0E0*I_ERI_F2xy_Pz_D2x_Py_bd;
  abcd[1322] = 4.0E0*I_ERI_F2xz_Pz_D2x_Py_bd;
  abcd[1323] = 4.0E0*I_ERI_Fx2y_Pz_D2x_Py_bd;
  abcd[1324] = 4.0E0*I_ERI_Fxyz_Pz_D2x_Py_bd;
  abcd[1325] = 4.0E0*I_ERI_Fx2z_Pz_D2x_Py_bd;
  abcd[1326] = 4.0E0*I_ERI_F3y_Pz_D2x_Py_bd;
  abcd[1327] = 4.0E0*I_ERI_F2yz_Pz_D2x_Py_bd;
  abcd[1328] = 4.0E0*I_ERI_Fy2z_Pz_D2x_Py_bd;
  abcd[1329] = 4.0E0*I_ERI_F3z_Pz_D2x_Py_bd;
  abcd[1330] = 4.0E0*I_ERI_F3x_Pz_Dxy_Py_bd;
  abcd[1331] = 4.0E0*I_ERI_F2xy_Pz_Dxy_Py_bd;
  abcd[1332] = 4.0E0*I_ERI_F2xz_Pz_Dxy_Py_bd;
  abcd[1333] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_Py_bd;
  abcd[1334] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_Py_bd;
  abcd[1335] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_Py_bd;
  abcd[1336] = 4.0E0*I_ERI_F3y_Pz_Dxy_Py_bd;
  abcd[1337] = 4.0E0*I_ERI_F2yz_Pz_Dxy_Py_bd;
  abcd[1338] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_Py_bd;
  abcd[1339] = 4.0E0*I_ERI_F3z_Pz_Dxy_Py_bd;
  abcd[1340] = 4.0E0*I_ERI_F3x_Pz_Dxz_Py_bd;
  abcd[1341] = 4.0E0*I_ERI_F2xy_Pz_Dxz_Py_bd;
  abcd[1342] = 4.0E0*I_ERI_F2xz_Pz_Dxz_Py_bd;
  abcd[1343] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_Py_bd;
  abcd[1344] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_Py_bd;
  abcd[1345] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_Py_bd;
  abcd[1346] = 4.0E0*I_ERI_F3y_Pz_Dxz_Py_bd;
  abcd[1347] = 4.0E0*I_ERI_F2yz_Pz_Dxz_Py_bd;
  abcd[1348] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_Py_bd;
  abcd[1349] = 4.0E0*I_ERI_F3z_Pz_Dxz_Py_bd;
  abcd[1350] = 4.0E0*I_ERI_F3x_Pz_D2y_Py_bd;
  abcd[1351] = 4.0E0*I_ERI_F2xy_Pz_D2y_Py_bd;
  abcd[1352] = 4.0E0*I_ERI_F2xz_Pz_D2y_Py_bd;
  abcd[1353] = 4.0E0*I_ERI_Fx2y_Pz_D2y_Py_bd;
  abcd[1354] = 4.0E0*I_ERI_Fxyz_Pz_D2y_Py_bd;
  abcd[1355] = 4.0E0*I_ERI_Fx2z_Pz_D2y_Py_bd;
  abcd[1356] = 4.0E0*I_ERI_F3y_Pz_D2y_Py_bd;
  abcd[1357] = 4.0E0*I_ERI_F2yz_Pz_D2y_Py_bd;
  abcd[1358] = 4.0E0*I_ERI_Fy2z_Pz_D2y_Py_bd;
  abcd[1359] = 4.0E0*I_ERI_F3z_Pz_D2y_Py_bd;
  abcd[1360] = 4.0E0*I_ERI_F3x_Pz_Dyz_Py_bd;
  abcd[1361] = 4.0E0*I_ERI_F2xy_Pz_Dyz_Py_bd;
  abcd[1362] = 4.0E0*I_ERI_F2xz_Pz_Dyz_Py_bd;
  abcd[1363] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_Py_bd;
  abcd[1364] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_Py_bd;
  abcd[1365] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_Py_bd;
  abcd[1366] = 4.0E0*I_ERI_F3y_Pz_Dyz_Py_bd;
  abcd[1367] = 4.0E0*I_ERI_F2yz_Pz_Dyz_Py_bd;
  abcd[1368] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_Py_bd;
  abcd[1369] = 4.0E0*I_ERI_F3z_Pz_Dyz_Py_bd;
  abcd[1370] = 4.0E0*I_ERI_F3x_Pz_D2z_Py_bd;
  abcd[1371] = 4.0E0*I_ERI_F2xy_Pz_D2z_Py_bd;
  abcd[1372] = 4.0E0*I_ERI_F2xz_Pz_D2z_Py_bd;
  abcd[1373] = 4.0E0*I_ERI_Fx2y_Pz_D2z_Py_bd;
  abcd[1374] = 4.0E0*I_ERI_Fxyz_Pz_D2z_Py_bd;
  abcd[1375] = 4.0E0*I_ERI_Fx2z_Pz_D2z_Py_bd;
  abcd[1376] = 4.0E0*I_ERI_F3y_Pz_D2z_Py_bd;
  abcd[1377] = 4.0E0*I_ERI_F2yz_Pz_D2z_Py_bd;
  abcd[1378] = 4.0E0*I_ERI_Fy2z_Pz_D2z_Py_bd;
  abcd[1379] = 4.0E0*I_ERI_F3z_Pz_D2z_Py_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dbz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_P_bd
   ************************************************************/
  abcd[1380] = 4.0E0*I_ERI_F3x_Pz_D2x_Pz_bd;
  abcd[1381] = 4.0E0*I_ERI_F2xy_Pz_D2x_Pz_bd;
  abcd[1382] = 4.0E0*I_ERI_F2xz_Pz_D2x_Pz_bd;
  abcd[1383] = 4.0E0*I_ERI_Fx2y_Pz_D2x_Pz_bd;
  abcd[1384] = 4.0E0*I_ERI_Fxyz_Pz_D2x_Pz_bd;
  abcd[1385] = 4.0E0*I_ERI_Fx2z_Pz_D2x_Pz_bd;
  abcd[1386] = 4.0E0*I_ERI_F3y_Pz_D2x_Pz_bd;
  abcd[1387] = 4.0E0*I_ERI_F2yz_Pz_D2x_Pz_bd;
  abcd[1388] = 4.0E0*I_ERI_Fy2z_Pz_D2x_Pz_bd;
  abcd[1389] = 4.0E0*I_ERI_F3z_Pz_D2x_Pz_bd;
  abcd[1390] = 4.0E0*I_ERI_F3x_Pz_Dxy_Pz_bd;
  abcd[1391] = 4.0E0*I_ERI_F2xy_Pz_Dxy_Pz_bd;
  abcd[1392] = 4.0E0*I_ERI_F2xz_Pz_Dxy_Pz_bd;
  abcd[1393] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_Pz_bd;
  abcd[1394] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_Pz_bd;
  abcd[1395] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_Pz_bd;
  abcd[1396] = 4.0E0*I_ERI_F3y_Pz_Dxy_Pz_bd;
  abcd[1397] = 4.0E0*I_ERI_F2yz_Pz_Dxy_Pz_bd;
  abcd[1398] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_Pz_bd;
  abcd[1399] = 4.0E0*I_ERI_F3z_Pz_Dxy_Pz_bd;
  abcd[1400] = 4.0E0*I_ERI_F3x_Pz_Dxz_Pz_bd;
  abcd[1401] = 4.0E0*I_ERI_F2xy_Pz_Dxz_Pz_bd;
  abcd[1402] = 4.0E0*I_ERI_F2xz_Pz_Dxz_Pz_bd;
  abcd[1403] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_Pz_bd;
  abcd[1404] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_Pz_bd;
  abcd[1405] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_Pz_bd;
  abcd[1406] = 4.0E0*I_ERI_F3y_Pz_Dxz_Pz_bd;
  abcd[1407] = 4.0E0*I_ERI_F2yz_Pz_Dxz_Pz_bd;
  abcd[1408] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_Pz_bd;
  abcd[1409] = 4.0E0*I_ERI_F3z_Pz_Dxz_Pz_bd;
  abcd[1410] = 4.0E0*I_ERI_F3x_Pz_D2y_Pz_bd;
  abcd[1411] = 4.0E0*I_ERI_F2xy_Pz_D2y_Pz_bd;
  abcd[1412] = 4.0E0*I_ERI_F2xz_Pz_D2y_Pz_bd;
  abcd[1413] = 4.0E0*I_ERI_Fx2y_Pz_D2y_Pz_bd;
  abcd[1414] = 4.0E0*I_ERI_Fxyz_Pz_D2y_Pz_bd;
  abcd[1415] = 4.0E0*I_ERI_Fx2z_Pz_D2y_Pz_bd;
  abcd[1416] = 4.0E0*I_ERI_F3y_Pz_D2y_Pz_bd;
  abcd[1417] = 4.0E0*I_ERI_F2yz_Pz_D2y_Pz_bd;
  abcd[1418] = 4.0E0*I_ERI_Fy2z_Pz_D2y_Pz_bd;
  abcd[1419] = 4.0E0*I_ERI_F3z_Pz_D2y_Pz_bd;
  abcd[1420] = 4.0E0*I_ERI_F3x_Pz_Dyz_Pz_bd;
  abcd[1421] = 4.0E0*I_ERI_F2xy_Pz_Dyz_Pz_bd;
  abcd[1422] = 4.0E0*I_ERI_F2xz_Pz_Dyz_Pz_bd;
  abcd[1423] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_Pz_bd;
  abcd[1424] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_Pz_bd;
  abcd[1425] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_Pz_bd;
  abcd[1426] = 4.0E0*I_ERI_F3y_Pz_Dyz_Pz_bd;
  abcd[1427] = 4.0E0*I_ERI_F2yz_Pz_Dyz_Pz_bd;
  abcd[1428] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_Pz_bd;
  abcd[1429] = 4.0E0*I_ERI_F3z_Pz_Dyz_Pz_bd;
  abcd[1430] = 4.0E0*I_ERI_F3x_Pz_D2z_Pz_bd;
  abcd[1431] = 4.0E0*I_ERI_F2xy_Pz_D2z_Pz_bd;
  abcd[1432] = 4.0E0*I_ERI_F2xz_Pz_D2z_Pz_bd;
  abcd[1433] = 4.0E0*I_ERI_Fx2y_Pz_D2z_Pz_bd;
  abcd[1434] = 4.0E0*I_ERI_Fxyz_Pz_D2z_Pz_bd;
  abcd[1435] = 4.0E0*I_ERI_Fx2z_Pz_D2z_Pz_bd;
  abcd[1436] = 4.0E0*I_ERI_F3y_Pz_D2z_Pz_bd;
  abcd[1437] = 4.0E0*I_ERI_F2yz_Pz_D2z_Pz_bd;
  abcd[1438] = 4.0E0*I_ERI_Fy2z_Pz_D2z_Pz_bd;
  abcd[1439] = 4.0E0*I_ERI_F3z_Pz_D2z_Pz_bd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1440] = 4.0E0*I_ERI_F3x_S_G4x_S_cc-2.0E0*2*I_ERI_F3x_S_D2x_S_c-2.0E0*3*I_ERI_F3x_S_D2x_S_c+2*1*I_ERI_F3x_S_S_S;
  abcd[1441] = 4.0E0*I_ERI_F2xy_S_G4x_S_cc-2.0E0*2*I_ERI_F2xy_S_D2x_S_c-2.0E0*3*I_ERI_F2xy_S_D2x_S_c+2*1*I_ERI_F2xy_S_S_S;
  abcd[1442] = 4.0E0*I_ERI_F2xz_S_G4x_S_cc-2.0E0*2*I_ERI_F2xz_S_D2x_S_c-2.0E0*3*I_ERI_F2xz_S_D2x_S_c+2*1*I_ERI_F2xz_S_S_S;
  abcd[1443] = 4.0E0*I_ERI_Fx2y_S_G4x_S_cc-2.0E0*2*I_ERI_Fx2y_S_D2x_S_c-2.0E0*3*I_ERI_Fx2y_S_D2x_S_c+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1444] = 4.0E0*I_ERI_Fxyz_S_G4x_S_cc-2.0E0*2*I_ERI_Fxyz_S_D2x_S_c-2.0E0*3*I_ERI_Fxyz_S_D2x_S_c+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1445] = 4.0E0*I_ERI_Fx2z_S_G4x_S_cc-2.0E0*2*I_ERI_Fx2z_S_D2x_S_c-2.0E0*3*I_ERI_Fx2z_S_D2x_S_c+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1446] = 4.0E0*I_ERI_F3y_S_G4x_S_cc-2.0E0*2*I_ERI_F3y_S_D2x_S_c-2.0E0*3*I_ERI_F3y_S_D2x_S_c+2*1*I_ERI_F3y_S_S_S;
  abcd[1447] = 4.0E0*I_ERI_F2yz_S_G4x_S_cc-2.0E0*2*I_ERI_F2yz_S_D2x_S_c-2.0E0*3*I_ERI_F2yz_S_D2x_S_c+2*1*I_ERI_F2yz_S_S_S;
  abcd[1448] = 4.0E0*I_ERI_Fy2z_S_G4x_S_cc-2.0E0*2*I_ERI_Fy2z_S_D2x_S_c-2.0E0*3*I_ERI_Fy2z_S_D2x_S_c+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1449] = 4.0E0*I_ERI_F3z_S_G4x_S_cc-2.0E0*2*I_ERI_F3z_S_D2x_S_c-2.0E0*3*I_ERI_F3z_S_D2x_S_c+2*1*I_ERI_F3z_S_S_S;
  abcd[1450] = 4.0E0*I_ERI_F3x_S_G3xy_S_cc-2.0E0*1*I_ERI_F3x_S_Dxy_S_c-2.0E0*2*I_ERI_F3x_S_Dxy_S_c;
  abcd[1451] = 4.0E0*I_ERI_F2xy_S_G3xy_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxy_S_c-2.0E0*2*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1452] = 4.0E0*I_ERI_F2xz_S_G3xy_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxy_S_c-2.0E0*2*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1453] = 4.0E0*I_ERI_Fx2y_S_G3xy_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_c-2.0E0*2*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1454] = 4.0E0*I_ERI_Fxyz_S_G3xy_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_c-2.0E0*2*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1455] = 4.0E0*I_ERI_Fx2z_S_G3xy_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_c-2.0E0*2*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1456] = 4.0E0*I_ERI_F3y_S_G3xy_S_cc-2.0E0*1*I_ERI_F3y_S_Dxy_S_c-2.0E0*2*I_ERI_F3y_S_Dxy_S_c;
  abcd[1457] = 4.0E0*I_ERI_F2yz_S_G3xy_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxy_S_c-2.0E0*2*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1458] = 4.0E0*I_ERI_Fy2z_S_G3xy_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_c-2.0E0*2*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1459] = 4.0E0*I_ERI_F3z_S_G3xy_S_cc-2.0E0*1*I_ERI_F3z_S_Dxy_S_c-2.0E0*2*I_ERI_F3z_S_Dxy_S_c;
  abcd[1460] = 4.0E0*I_ERI_F3x_S_G3xz_S_cc-2.0E0*1*I_ERI_F3x_S_Dxz_S_c-2.0E0*2*I_ERI_F3x_S_Dxz_S_c;
  abcd[1461] = 4.0E0*I_ERI_F2xy_S_G3xz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxz_S_c-2.0E0*2*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1462] = 4.0E0*I_ERI_F2xz_S_G3xz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxz_S_c-2.0E0*2*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1463] = 4.0E0*I_ERI_Fx2y_S_G3xz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_c-2.0E0*2*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1464] = 4.0E0*I_ERI_Fxyz_S_G3xz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_c-2.0E0*2*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1465] = 4.0E0*I_ERI_Fx2z_S_G3xz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_c-2.0E0*2*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1466] = 4.0E0*I_ERI_F3y_S_G3xz_S_cc-2.0E0*1*I_ERI_F3y_S_Dxz_S_c-2.0E0*2*I_ERI_F3y_S_Dxz_S_c;
  abcd[1467] = 4.0E0*I_ERI_F2yz_S_G3xz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxz_S_c-2.0E0*2*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1468] = 4.0E0*I_ERI_Fy2z_S_G3xz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_c-2.0E0*2*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1469] = 4.0E0*I_ERI_F3z_S_G3xz_S_cc-2.0E0*1*I_ERI_F3z_S_Dxz_S_c-2.0E0*2*I_ERI_F3z_S_Dxz_S_c;
  abcd[1470] = 4.0E0*I_ERI_F3x_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3x_S_D2y_S_c;
  abcd[1471] = 4.0E0*I_ERI_F2xy_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xy_S_D2y_S_c;
  abcd[1472] = 4.0E0*I_ERI_F2xz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xz_S_D2y_S_c;
  abcd[1473] = 4.0E0*I_ERI_Fx2y_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2y_S_c;
  abcd[1474] = 4.0E0*I_ERI_Fxyz_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2y_S_c;
  abcd[1475] = 4.0E0*I_ERI_Fx2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2y_S_c;
  abcd[1476] = 4.0E0*I_ERI_F3y_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3y_S_D2y_S_c;
  abcd[1477] = 4.0E0*I_ERI_F2yz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2yz_S_D2y_S_c;
  abcd[1478] = 4.0E0*I_ERI_Fy2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2y_S_c;
  abcd[1479] = 4.0E0*I_ERI_F3z_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3z_S_D2y_S_c;
  abcd[1480] = 4.0E0*I_ERI_F3x_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3x_S_Dyz_S_c;
  abcd[1481] = 4.0E0*I_ERI_F2xy_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1482] = 4.0E0*I_ERI_F2xz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1483] = 4.0E0*I_ERI_Fx2y_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1484] = 4.0E0*I_ERI_Fxyz_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1485] = 4.0E0*I_ERI_Fx2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1486] = 4.0E0*I_ERI_F3y_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3y_S_Dyz_S_c;
  abcd[1487] = 4.0E0*I_ERI_F2yz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1488] = 4.0E0*I_ERI_Fy2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1489] = 4.0E0*I_ERI_F3z_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3z_S_Dyz_S_c;
  abcd[1490] = 4.0E0*I_ERI_F3x_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2z_S_c;
  abcd[1491] = 4.0E0*I_ERI_F2xy_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2z_S_c;
  abcd[1492] = 4.0E0*I_ERI_F2xz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2z_S_c;
  abcd[1493] = 4.0E0*I_ERI_Fx2y_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2z_S_c;
  abcd[1494] = 4.0E0*I_ERI_Fxyz_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2z_S_c;
  abcd[1495] = 4.0E0*I_ERI_Fx2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2z_S_c;
  abcd[1496] = 4.0E0*I_ERI_F3y_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2z_S_c;
  abcd[1497] = 4.0E0*I_ERI_F2yz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2z_S_c;
  abcd[1498] = 4.0E0*I_ERI_Fy2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2z_S_c;
  abcd[1499] = 4.0E0*I_ERI_F3z_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1500] = 4.0E0*I_ERI_F3x_S_G3xy_S_cc-2.0E0*2*I_ERI_F3x_S_Dxy_S_c;
  abcd[1501] = 4.0E0*I_ERI_F2xy_S_G3xy_S_cc-2.0E0*2*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1502] = 4.0E0*I_ERI_F2xz_S_G3xy_S_cc-2.0E0*2*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1503] = 4.0E0*I_ERI_Fx2y_S_G3xy_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1504] = 4.0E0*I_ERI_Fxyz_S_G3xy_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1505] = 4.0E0*I_ERI_Fx2z_S_G3xy_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1506] = 4.0E0*I_ERI_F3y_S_G3xy_S_cc-2.0E0*2*I_ERI_F3y_S_Dxy_S_c;
  abcd[1507] = 4.0E0*I_ERI_F2yz_S_G3xy_S_cc-2.0E0*2*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1508] = 4.0E0*I_ERI_Fy2z_S_G3xy_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1509] = 4.0E0*I_ERI_F3z_S_G3xy_S_cc-2.0E0*2*I_ERI_F3z_S_Dxy_S_c;
  abcd[1510] = 4.0E0*I_ERI_F3x_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3x_S_D2x_S_c-2.0E0*1*I_ERI_F3x_S_D2y_S_c+1*I_ERI_F3x_S_S_S;
  abcd[1511] = 4.0E0*I_ERI_F2xy_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xy_S_D2x_S_c-2.0E0*1*I_ERI_F2xy_S_D2y_S_c+1*I_ERI_F2xy_S_S_S;
  abcd[1512] = 4.0E0*I_ERI_F2xz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xz_S_D2x_S_c-2.0E0*1*I_ERI_F2xz_S_D2y_S_c+1*I_ERI_F2xz_S_S_S;
  abcd[1513] = 4.0E0*I_ERI_Fx2y_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2x_S_c-2.0E0*1*I_ERI_Fx2y_S_D2y_S_c+1*I_ERI_Fx2y_S_S_S;
  abcd[1514] = 4.0E0*I_ERI_Fxyz_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2x_S_c-2.0E0*1*I_ERI_Fxyz_S_D2y_S_c+1*I_ERI_Fxyz_S_S_S;
  abcd[1515] = 4.0E0*I_ERI_Fx2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2x_S_c-2.0E0*1*I_ERI_Fx2z_S_D2y_S_c+1*I_ERI_Fx2z_S_S_S;
  abcd[1516] = 4.0E0*I_ERI_F3y_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3y_S_D2x_S_c-2.0E0*1*I_ERI_F3y_S_D2y_S_c+1*I_ERI_F3y_S_S_S;
  abcd[1517] = 4.0E0*I_ERI_F2yz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2yz_S_D2x_S_c-2.0E0*1*I_ERI_F2yz_S_D2y_S_c+1*I_ERI_F2yz_S_S_S;
  abcd[1518] = 4.0E0*I_ERI_Fy2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2x_S_c-2.0E0*1*I_ERI_Fy2z_S_D2y_S_c+1*I_ERI_Fy2z_S_S_S;
  abcd[1519] = 4.0E0*I_ERI_F3z_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3z_S_D2x_S_c-2.0E0*1*I_ERI_F3z_S_D2y_S_c+1*I_ERI_F3z_S_S_S;
  abcd[1520] = 4.0E0*I_ERI_F3x_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3x_S_Dyz_S_c;
  abcd[1521] = 4.0E0*I_ERI_F2xy_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1522] = 4.0E0*I_ERI_F2xz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1523] = 4.0E0*I_ERI_Fx2y_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1524] = 4.0E0*I_ERI_Fxyz_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1525] = 4.0E0*I_ERI_Fx2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1526] = 4.0E0*I_ERI_F3y_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3y_S_Dyz_S_c;
  abcd[1527] = 4.0E0*I_ERI_F2yz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1528] = 4.0E0*I_ERI_Fy2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1529] = 4.0E0*I_ERI_F3z_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3z_S_Dyz_S_c;
  abcd[1530] = 4.0E0*I_ERI_F3x_S_Gx3y_S_cc-2.0E0*2*I_ERI_F3x_S_Dxy_S_c;
  abcd[1531] = 4.0E0*I_ERI_F2xy_S_Gx3y_S_cc-2.0E0*2*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1532] = 4.0E0*I_ERI_F2xz_S_Gx3y_S_cc-2.0E0*2*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1533] = 4.0E0*I_ERI_Fx2y_S_Gx3y_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1534] = 4.0E0*I_ERI_Fxyz_S_Gx3y_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1535] = 4.0E0*I_ERI_Fx2z_S_Gx3y_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1536] = 4.0E0*I_ERI_F3y_S_Gx3y_S_cc-2.0E0*2*I_ERI_F3y_S_Dxy_S_c;
  abcd[1537] = 4.0E0*I_ERI_F2yz_S_Gx3y_S_cc-2.0E0*2*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1538] = 4.0E0*I_ERI_Fy2z_S_Gx3y_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1539] = 4.0E0*I_ERI_F3z_S_Gx3y_S_cc-2.0E0*2*I_ERI_F3z_S_Dxy_S_c;
  abcd[1540] = 4.0E0*I_ERI_F3x_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3x_S_Dxz_S_c;
  abcd[1541] = 4.0E0*I_ERI_F2xy_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1542] = 4.0E0*I_ERI_F2xz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1543] = 4.0E0*I_ERI_Fx2y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1544] = 4.0E0*I_ERI_Fxyz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1545] = 4.0E0*I_ERI_Fx2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1546] = 4.0E0*I_ERI_F3y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3y_S_Dxz_S_c;
  abcd[1547] = 4.0E0*I_ERI_F2yz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1548] = 4.0E0*I_ERI_Fy2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1549] = 4.0E0*I_ERI_F3z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3z_S_Dxz_S_c;
  abcd[1550] = 4.0E0*I_ERI_F3x_S_Gxy2z_S_cc;
  abcd[1551] = 4.0E0*I_ERI_F2xy_S_Gxy2z_S_cc;
  abcd[1552] = 4.0E0*I_ERI_F2xz_S_Gxy2z_S_cc;
  abcd[1553] = 4.0E0*I_ERI_Fx2y_S_Gxy2z_S_cc;
  abcd[1554] = 4.0E0*I_ERI_Fxyz_S_Gxy2z_S_cc;
  abcd[1555] = 4.0E0*I_ERI_Fx2z_S_Gxy2z_S_cc;
  abcd[1556] = 4.0E0*I_ERI_F3y_S_Gxy2z_S_cc;
  abcd[1557] = 4.0E0*I_ERI_F2yz_S_Gxy2z_S_cc;
  abcd[1558] = 4.0E0*I_ERI_Fy2z_S_Gxy2z_S_cc;
  abcd[1559] = 4.0E0*I_ERI_F3z_S_Gxy2z_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1560] = 4.0E0*I_ERI_F3x_S_G3xz_S_cc-2.0E0*2*I_ERI_F3x_S_Dxz_S_c;
  abcd[1561] = 4.0E0*I_ERI_F2xy_S_G3xz_S_cc-2.0E0*2*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1562] = 4.0E0*I_ERI_F2xz_S_G3xz_S_cc-2.0E0*2*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1563] = 4.0E0*I_ERI_Fx2y_S_G3xz_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1564] = 4.0E0*I_ERI_Fxyz_S_G3xz_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1565] = 4.0E0*I_ERI_Fx2z_S_G3xz_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1566] = 4.0E0*I_ERI_F3y_S_G3xz_S_cc-2.0E0*2*I_ERI_F3y_S_Dxz_S_c;
  abcd[1567] = 4.0E0*I_ERI_F2yz_S_G3xz_S_cc-2.0E0*2*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1568] = 4.0E0*I_ERI_Fy2z_S_G3xz_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1569] = 4.0E0*I_ERI_F3z_S_G3xz_S_cc-2.0E0*2*I_ERI_F3z_S_Dxz_S_c;
  abcd[1570] = 4.0E0*I_ERI_F3x_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3x_S_Dyz_S_c;
  abcd[1571] = 4.0E0*I_ERI_F2xy_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1572] = 4.0E0*I_ERI_F2xz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1573] = 4.0E0*I_ERI_Fx2y_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1574] = 4.0E0*I_ERI_Fxyz_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1575] = 4.0E0*I_ERI_Fx2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1576] = 4.0E0*I_ERI_F3y_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3y_S_Dyz_S_c;
  abcd[1577] = 4.0E0*I_ERI_F2yz_S_G2xyz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1578] = 4.0E0*I_ERI_Fy2z_S_G2xyz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1579] = 4.0E0*I_ERI_F3z_S_G2xyz_S_cc-2.0E0*1*I_ERI_F3z_S_Dyz_S_c;
  abcd[1580] = 4.0E0*I_ERI_F3x_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2x_S_c-2.0E0*1*I_ERI_F3x_S_D2z_S_c+1*I_ERI_F3x_S_S_S;
  abcd[1581] = 4.0E0*I_ERI_F2xy_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2x_S_c-2.0E0*1*I_ERI_F2xy_S_D2z_S_c+1*I_ERI_F2xy_S_S_S;
  abcd[1582] = 4.0E0*I_ERI_F2xz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2x_S_c-2.0E0*1*I_ERI_F2xz_S_D2z_S_c+1*I_ERI_F2xz_S_S_S;
  abcd[1583] = 4.0E0*I_ERI_Fx2y_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2x_S_c-2.0E0*1*I_ERI_Fx2y_S_D2z_S_c+1*I_ERI_Fx2y_S_S_S;
  abcd[1584] = 4.0E0*I_ERI_Fxyz_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2x_S_c-2.0E0*1*I_ERI_Fxyz_S_D2z_S_c+1*I_ERI_Fxyz_S_S_S;
  abcd[1585] = 4.0E0*I_ERI_Fx2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2x_S_c-2.0E0*1*I_ERI_Fx2z_S_D2z_S_c+1*I_ERI_Fx2z_S_S_S;
  abcd[1586] = 4.0E0*I_ERI_F3y_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2x_S_c-2.0E0*1*I_ERI_F3y_S_D2z_S_c+1*I_ERI_F3y_S_S_S;
  abcd[1587] = 4.0E0*I_ERI_F2yz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2x_S_c-2.0E0*1*I_ERI_F2yz_S_D2z_S_c+1*I_ERI_F2yz_S_S_S;
  abcd[1588] = 4.0E0*I_ERI_Fy2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2x_S_c-2.0E0*1*I_ERI_Fy2z_S_D2z_S_c+1*I_ERI_Fy2z_S_S_S;
  abcd[1589] = 4.0E0*I_ERI_F3z_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2x_S_c-2.0E0*1*I_ERI_F3z_S_D2z_S_c+1*I_ERI_F3z_S_S_S;
  abcd[1590] = 4.0E0*I_ERI_F3x_S_Gx2yz_S_cc;
  abcd[1591] = 4.0E0*I_ERI_F2xy_S_Gx2yz_S_cc;
  abcd[1592] = 4.0E0*I_ERI_F2xz_S_Gx2yz_S_cc;
  abcd[1593] = 4.0E0*I_ERI_Fx2y_S_Gx2yz_S_cc;
  abcd[1594] = 4.0E0*I_ERI_Fxyz_S_Gx2yz_S_cc;
  abcd[1595] = 4.0E0*I_ERI_Fx2z_S_Gx2yz_S_cc;
  abcd[1596] = 4.0E0*I_ERI_F3y_S_Gx2yz_S_cc;
  abcd[1597] = 4.0E0*I_ERI_F2yz_S_Gx2yz_S_cc;
  abcd[1598] = 4.0E0*I_ERI_Fy2z_S_Gx2yz_S_cc;
  abcd[1599] = 4.0E0*I_ERI_F3z_S_Gx2yz_S_cc;
  abcd[1600] = 4.0E0*I_ERI_F3x_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3x_S_Dxy_S_c;
  abcd[1601] = 4.0E0*I_ERI_F2xy_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1602] = 4.0E0*I_ERI_F2xz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1603] = 4.0E0*I_ERI_Fx2y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1604] = 4.0E0*I_ERI_Fxyz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1605] = 4.0E0*I_ERI_Fx2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1606] = 4.0E0*I_ERI_F3y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3y_S_Dxy_S_c;
  abcd[1607] = 4.0E0*I_ERI_F2yz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1608] = 4.0E0*I_ERI_Fy2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1609] = 4.0E0*I_ERI_F3z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3z_S_Dxy_S_c;
  abcd[1610] = 4.0E0*I_ERI_F3x_S_Gx3z_S_cc-2.0E0*2*I_ERI_F3x_S_Dxz_S_c;
  abcd[1611] = 4.0E0*I_ERI_F2xy_S_Gx3z_S_cc-2.0E0*2*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1612] = 4.0E0*I_ERI_F2xz_S_Gx3z_S_cc-2.0E0*2*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1613] = 4.0E0*I_ERI_Fx2y_S_Gx3z_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1614] = 4.0E0*I_ERI_Fxyz_S_Gx3z_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1615] = 4.0E0*I_ERI_Fx2z_S_Gx3z_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1616] = 4.0E0*I_ERI_F3y_S_Gx3z_S_cc-2.0E0*2*I_ERI_F3y_S_Dxz_S_c;
  abcd[1617] = 4.0E0*I_ERI_F2yz_S_Gx3z_S_cc-2.0E0*2*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1618] = 4.0E0*I_ERI_Fy2z_S_Gx3z_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1619] = 4.0E0*I_ERI_F3z_S_Gx3z_S_cc-2.0E0*2*I_ERI_F3z_S_Dxz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1620] = 4.0E0*I_ERI_F3x_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3x_S_D2x_S_c;
  abcd[1621] = 4.0E0*I_ERI_F2xy_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xy_S_D2x_S_c;
  abcd[1622] = 4.0E0*I_ERI_F2xz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2xz_S_D2x_S_c;
  abcd[1623] = 4.0E0*I_ERI_Fx2y_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2x_S_c;
  abcd[1624] = 4.0E0*I_ERI_Fxyz_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2x_S_c;
  abcd[1625] = 4.0E0*I_ERI_Fx2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2x_S_c;
  abcd[1626] = 4.0E0*I_ERI_F3y_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3y_S_D2x_S_c;
  abcd[1627] = 4.0E0*I_ERI_F2yz_S_G2x2y_S_cc-2.0E0*1*I_ERI_F2yz_S_D2x_S_c;
  abcd[1628] = 4.0E0*I_ERI_Fy2z_S_G2x2y_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2x_S_c;
  abcd[1629] = 4.0E0*I_ERI_F3z_S_G2x2y_S_cc-2.0E0*1*I_ERI_F3z_S_D2x_S_c;
  abcd[1630] = 4.0E0*I_ERI_F3x_S_Gx3y_S_cc-2.0E0*1*I_ERI_F3x_S_Dxy_S_c-2.0E0*2*I_ERI_F3x_S_Dxy_S_c;
  abcd[1631] = 4.0E0*I_ERI_F2xy_S_Gx3y_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxy_S_c-2.0E0*2*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1632] = 4.0E0*I_ERI_F2xz_S_Gx3y_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxy_S_c-2.0E0*2*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1633] = 4.0E0*I_ERI_Fx2y_S_Gx3y_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_c-2.0E0*2*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1634] = 4.0E0*I_ERI_Fxyz_S_Gx3y_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_c-2.0E0*2*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1635] = 4.0E0*I_ERI_Fx2z_S_Gx3y_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_c-2.0E0*2*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1636] = 4.0E0*I_ERI_F3y_S_Gx3y_S_cc-2.0E0*1*I_ERI_F3y_S_Dxy_S_c-2.0E0*2*I_ERI_F3y_S_Dxy_S_c;
  abcd[1637] = 4.0E0*I_ERI_F2yz_S_Gx3y_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxy_S_c-2.0E0*2*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1638] = 4.0E0*I_ERI_Fy2z_S_Gx3y_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_c-2.0E0*2*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1639] = 4.0E0*I_ERI_F3z_S_Gx3y_S_cc-2.0E0*1*I_ERI_F3z_S_Dxy_S_c-2.0E0*2*I_ERI_F3z_S_Dxy_S_c;
  abcd[1640] = 4.0E0*I_ERI_F3x_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3x_S_Dxz_S_c;
  abcd[1641] = 4.0E0*I_ERI_F2xy_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1642] = 4.0E0*I_ERI_F2xz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1643] = 4.0E0*I_ERI_Fx2y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1644] = 4.0E0*I_ERI_Fxyz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1645] = 4.0E0*I_ERI_Fx2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1646] = 4.0E0*I_ERI_F3y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3y_S_Dxz_S_c;
  abcd[1647] = 4.0E0*I_ERI_F2yz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1648] = 4.0E0*I_ERI_Fy2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1649] = 4.0E0*I_ERI_F3z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3z_S_Dxz_S_c;
  abcd[1650] = 4.0E0*I_ERI_F3x_S_G4y_S_cc-2.0E0*2*I_ERI_F3x_S_D2y_S_c-2.0E0*3*I_ERI_F3x_S_D2y_S_c+2*1*I_ERI_F3x_S_S_S;
  abcd[1651] = 4.0E0*I_ERI_F2xy_S_G4y_S_cc-2.0E0*2*I_ERI_F2xy_S_D2y_S_c-2.0E0*3*I_ERI_F2xy_S_D2y_S_c+2*1*I_ERI_F2xy_S_S_S;
  abcd[1652] = 4.0E0*I_ERI_F2xz_S_G4y_S_cc-2.0E0*2*I_ERI_F2xz_S_D2y_S_c-2.0E0*3*I_ERI_F2xz_S_D2y_S_c+2*1*I_ERI_F2xz_S_S_S;
  abcd[1653] = 4.0E0*I_ERI_Fx2y_S_G4y_S_cc-2.0E0*2*I_ERI_Fx2y_S_D2y_S_c-2.0E0*3*I_ERI_Fx2y_S_D2y_S_c+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1654] = 4.0E0*I_ERI_Fxyz_S_G4y_S_cc-2.0E0*2*I_ERI_Fxyz_S_D2y_S_c-2.0E0*3*I_ERI_Fxyz_S_D2y_S_c+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1655] = 4.0E0*I_ERI_Fx2z_S_G4y_S_cc-2.0E0*2*I_ERI_Fx2z_S_D2y_S_c-2.0E0*3*I_ERI_Fx2z_S_D2y_S_c+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1656] = 4.0E0*I_ERI_F3y_S_G4y_S_cc-2.0E0*2*I_ERI_F3y_S_D2y_S_c-2.0E0*3*I_ERI_F3y_S_D2y_S_c+2*1*I_ERI_F3y_S_S_S;
  abcd[1657] = 4.0E0*I_ERI_F2yz_S_G4y_S_cc-2.0E0*2*I_ERI_F2yz_S_D2y_S_c-2.0E0*3*I_ERI_F2yz_S_D2y_S_c+2*1*I_ERI_F2yz_S_S_S;
  abcd[1658] = 4.0E0*I_ERI_Fy2z_S_G4y_S_cc-2.0E0*2*I_ERI_Fy2z_S_D2y_S_c-2.0E0*3*I_ERI_Fy2z_S_D2y_S_c+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1659] = 4.0E0*I_ERI_F3z_S_G4y_S_cc-2.0E0*2*I_ERI_F3z_S_D2y_S_c-2.0E0*3*I_ERI_F3z_S_D2y_S_c+2*1*I_ERI_F3z_S_S_S;
  abcd[1660] = 4.0E0*I_ERI_F3x_S_G3yz_S_cc-2.0E0*1*I_ERI_F3x_S_Dyz_S_c-2.0E0*2*I_ERI_F3x_S_Dyz_S_c;
  abcd[1661] = 4.0E0*I_ERI_F2xy_S_G3yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dyz_S_c-2.0E0*2*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1662] = 4.0E0*I_ERI_F2xz_S_G3yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dyz_S_c-2.0E0*2*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1663] = 4.0E0*I_ERI_Fx2y_S_G3yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_c-2.0E0*2*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1664] = 4.0E0*I_ERI_Fxyz_S_G3yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_c-2.0E0*2*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1665] = 4.0E0*I_ERI_Fx2z_S_G3yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_c-2.0E0*2*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1666] = 4.0E0*I_ERI_F3y_S_G3yz_S_cc-2.0E0*1*I_ERI_F3y_S_Dyz_S_c-2.0E0*2*I_ERI_F3y_S_Dyz_S_c;
  abcd[1667] = 4.0E0*I_ERI_F2yz_S_G3yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dyz_S_c-2.0E0*2*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1668] = 4.0E0*I_ERI_Fy2z_S_G3yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_c-2.0E0*2*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1669] = 4.0E0*I_ERI_F3z_S_G3yz_S_cc-2.0E0*1*I_ERI_F3z_S_Dyz_S_c-2.0E0*2*I_ERI_F3z_S_Dyz_S_c;
  abcd[1670] = 4.0E0*I_ERI_F3x_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2z_S_c;
  abcd[1671] = 4.0E0*I_ERI_F2xy_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2z_S_c;
  abcd[1672] = 4.0E0*I_ERI_F2xz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2z_S_c;
  abcd[1673] = 4.0E0*I_ERI_Fx2y_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2z_S_c;
  abcd[1674] = 4.0E0*I_ERI_Fxyz_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2z_S_c;
  abcd[1675] = 4.0E0*I_ERI_Fx2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2z_S_c;
  abcd[1676] = 4.0E0*I_ERI_F3y_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2z_S_c;
  abcd[1677] = 4.0E0*I_ERI_F2yz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2z_S_c;
  abcd[1678] = 4.0E0*I_ERI_Fy2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2z_S_c;
  abcd[1679] = 4.0E0*I_ERI_F3z_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1680] = 4.0E0*I_ERI_F3x_S_G2xyz_S_cc;
  abcd[1681] = 4.0E0*I_ERI_F2xy_S_G2xyz_S_cc;
  abcd[1682] = 4.0E0*I_ERI_F2xz_S_G2xyz_S_cc;
  abcd[1683] = 4.0E0*I_ERI_Fx2y_S_G2xyz_S_cc;
  abcd[1684] = 4.0E0*I_ERI_Fxyz_S_G2xyz_S_cc;
  abcd[1685] = 4.0E0*I_ERI_Fx2z_S_G2xyz_S_cc;
  abcd[1686] = 4.0E0*I_ERI_F3y_S_G2xyz_S_cc;
  abcd[1687] = 4.0E0*I_ERI_F2yz_S_G2xyz_S_cc;
  abcd[1688] = 4.0E0*I_ERI_Fy2z_S_G2xyz_S_cc;
  abcd[1689] = 4.0E0*I_ERI_F3z_S_G2xyz_S_cc;
  abcd[1690] = 4.0E0*I_ERI_F3x_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3x_S_Dxz_S_c;
  abcd[1691] = 4.0E0*I_ERI_F2xy_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1692] = 4.0E0*I_ERI_F2xz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1693] = 4.0E0*I_ERI_Fx2y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1694] = 4.0E0*I_ERI_Fxyz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1695] = 4.0E0*I_ERI_Fx2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1696] = 4.0E0*I_ERI_F3y_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3y_S_Dxz_S_c;
  abcd[1697] = 4.0E0*I_ERI_F2yz_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1698] = 4.0E0*I_ERI_Fy2z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1699] = 4.0E0*I_ERI_F3z_S_Gx2yz_S_cc-2.0E0*1*I_ERI_F3z_S_Dxz_S_c;
  abcd[1700] = 4.0E0*I_ERI_F3x_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3x_S_Dxy_S_c;
  abcd[1701] = 4.0E0*I_ERI_F2xy_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1702] = 4.0E0*I_ERI_F2xz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1703] = 4.0E0*I_ERI_Fx2y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1704] = 4.0E0*I_ERI_Fxyz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1705] = 4.0E0*I_ERI_Fx2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1706] = 4.0E0*I_ERI_F3y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3y_S_Dxy_S_c;
  abcd[1707] = 4.0E0*I_ERI_F2yz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1708] = 4.0E0*I_ERI_Fy2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1709] = 4.0E0*I_ERI_F3z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3z_S_Dxy_S_c;
  abcd[1710] = 4.0E0*I_ERI_F3x_S_G3yz_S_cc-2.0E0*2*I_ERI_F3x_S_Dyz_S_c;
  abcd[1711] = 4.0E0*I_ERI_F2xy_S_G3yz_S_cc-2.0E0*2*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1712] = 4.0E0*I_ERI_F2xz_S_G3yz_S_cc-2.0E0*2*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1713] = 4.0E0*I_ERI_Fx2y_S_G3yz_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1714] = 4.0E0*I_ERI_Fxyz_S_G3yz_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1715] = 4.0E0*I_ERI_Fx2z_S_G3yz_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1716] = 4.0E0*I_ERI_F3y_S_G3yz_S_cc-2.0E0*2*I_ERI_F3y_S_Dyz_S_c;
  abcd[1717] = 4.0E0*I_ERI_F2yz_S_G3yz_S_cc-2.0E0*2*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1718] = 4.0E0*I_ERI_Fy2z_S_G3yz_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1719] = 4.0E0*I_ERI_F3z_S_G3yz_S_cc-2.0E0*2*I_ERI_F3z_S_Dyz_S_c;
  abcd[1720] = 4.0E0*I_ERI_F3x_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2y_S_c-2.0E0*1*I_ERI_F3x_S_D2z_S_c+1*I_ERI_F3x_S_S_S;
  abcd[1721] = 4.0E0*I_ERI_F2xy_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2y_S_c-2.0E0*1*I_ERI_F2xy_S_D2z_S_c+1*I_ERI_F2xy_S_S_S;
  abcd[1722] = 4.0E0*I_ERI_F2xz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2y_S_c-2.0E0*1*I_ERI_F2xz_S_D2z_S_c+1*I_ERI_F2xz_S_S_S;
  abcd[1723] = 4.0E0*I_ERI_Fx2y_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2y_S_c-2.0E0*1*I_ERI_Fx2y_S_D2z_S_c+1*I_ERI_Fx2y_S_S_S;
  abcd[1724] = 4.0E0*I_ERI_Fxyz_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2y_S_c-2.0E0*1*I_ERI_Fxyz_S_D2z_S_c+1*I_ERI_Fxyz_S_S_S;
  abcd[1725] = 4.0E0*I_ERI_Fx2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2y_S_c-2.0E0*1*I_ERI_Fx2z_S_D2z_S_c+1*I_ERI_Fx2z_S_S_S;
  abcd[1726] = 4.0E0*I_ERI_F3y_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2y_S_c-2.0E0*1*I_ERI_F3y_S_D2z_S_c+1*I_ERI_F3y_S_S_S;
  abcd[1727] = 4.0E0*I_ERI_F2yz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2y_S_c-2.0E0*1*I_ERI_F2yz_S_D2z_S_c+1*I_ERI_F2yz_S_S_S;
  abcd[1728] = 4.0E0*I_ERI_Fy2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2y_S_c-2.0E0*1*I_ERI_Fy2z_S_D2z_S_c+1*I_ERI_Fy2z_S_S_S;
  abcd[1729] = 4.0E0*I_ERI_F3z_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2y_S_c-2.0E0*1*I_ERI_F3z_S_D2z_S_c+1*I_ERI_F3z_S_S_S;
  abcd[1730] = 4.0E0*I_ERI_F3x_S_Gy3z_S_cc-2.0E0*2*I_ERI_F3x_S_Dyz_S_c;
  abcd[1731] = 4.0E0*I_ERI_F2xy_S_Gy3z_S_cc-2.0E0*2*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1732] = 4.0E0*I_ERI_F2xz_S_Gy3z_S_cc-2.0E0*2*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1733] = 4.0E0*I_ERI_Fx2y_S_Gy3z_S_cc-2.0E0*2*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1734] = 4.0E0*I_ERI_Fxyz_S_Gy3z_S_cc-2.0E0*2*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1735] = 4.0E0*I_ERI_Fx2z_S_Gy3z_S_cc-2.0E0*2*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1736] = 4.0E0*I_ERI_F3y_S_Gy3z_S_cc-2.0E0*2*I_ERI_F3y_S_Dyz_S_c;
  abcd[1737] = 4.0E0*I_ERI_F2yz_S_Gy3z_S_cc-2.0E0*2*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1738] = 4.0E0*I_ERI_Fy2z_S_Gy3z_S_cc-2.0E0*2*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1739] = 4.0E0*I_ERI_F3z_S_Gy3z_S_cc-2.0E0*2*I_ERI_F3z_S_Dyz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_D_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1740] = 4.0E0*I_ERI_F3x_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2x_S_c;
  abcd[1741] = 4.0E0*I_ERI_F2xy_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2x_S_c;
  abcd[1742] = 4.0E0*I_ERI_F2xz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2x_S_c;
  abcd[1743] = 4.0E0*I_ERI_Fx2y_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2x_S_c;
  abcd[1744] = 4.0E0*I_ERI_Fxyz_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2x_S_c;
  abcd[1745] = 4.0E0*I_ERI_Fx2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2x_S_c;
  abcd[1746] = 4.0E0*I_ERI_F3y_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2x_S_c;
  abcd[1747] = 4.0E0*I_ERI_F2yz_S_G2x2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2x_S_c;
  abcd[1748] = 4.0E0*I_ERI_Fy2z_S_G2x2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2x_S_c;
  abcd[1749] = 4.0E0*I_ERI_F3z_S_G2x2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2x_S_c;
  abcd[1750] = 4.0E0*I_ERI_F3x_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3x_S_Dxy_S_c;
  abcd[1751] = 4.0E0*I_ERI_F2xy_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxy_S_c;
  abcd[1752] = 4.0E0*I_ERI_F2xz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxy_S_c;
  abcd[1753] = 4.0E0*I_ERI_Fx2y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_c;
  abcd[1754] = 4.0E0*I_ERI_Fxyz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_c;
  abcd[1755] = 4.0E0*I_ERI_Fx2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_c;
  abcd[1756] = 4.0E0*I_ERI_F3y_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3y_S_Dxy_S_c;
  abcd[1757] = 4.0E0*I_ERI_F2yz_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxy_S_c;
  abcd[1758] = 4.0E0*I_ERI_Fy2z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_c;
  abcd[1759] = 4.0E0*I_ERI_F3z_S_Gxy2z_S_cc-2.0E0*1*I_ERI_F3z_S_Dxy_S_c;
  abcd[1760] = 4.0E0*I_ERI_F3x_S_Gx3z_S_cc-2.0E0*1*I_ERI_F3x_S_Dxz_S_c-2.0E0*2*I_ERI_F3x_S_Dxz_S_c;
  abcd[1761] = 4.0E0*I_ERI_F2xy_S_Gx3z_S_cc-2.0E0*1*I_ERI_F2xy_S_Dxz_S_c-2.0E0*2*I_ERI_F2xy_S_Dxz_S_c;
  abcd[1762] = 4.0E0*I_ERI_F2xz_S_Gx3z_S_cc-2.0E0*1*I_ERI_F2xz_S_Dxz_S_c-2.0E0*2*I_ERI_F2xz_S_Dxz_S_c;
  abcd[1763] = 4.0E0*I_ERI_Fx2y_S_Gx3z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_c-2.0E0*2*I_ERI_Fx2y_S_Dxz_S_c;
  abcd[1764] = 4.0E0*I_ERI_Fxyz_S_Gx3z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_c-2.0E0*2*I_ERI_Fxyz_S_Dxz_S_c;
  abcd[1765] = 4.0E0*I_ERI_Fx2z_S_Gx3z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_c-2.0E0*2*I_ERI_Fx2z_S_Dxz_S_c;
  abcd[1766] = 4.0E0*I_ERI_F3y_S_Gx3z_S_cc-2.0E0*1*I_ERI_F3y_S_Dxz_S_c-2.0E0*2*I_ERI_F3y_S_Dxz_S_c;
  abcd[1767] = 4.0E0*I_ERI_F2yz_S_Gx3z_S_cc-2.0E0*1*I_ERI_F2yz_S_Dxz_S_c-2.0E0*2*I_ERI_F2yz_S_Dxz_S_c;
  abcd[1768] = 4.0E0*I_ERI_Fy2z_S_Gx3z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_c-2.0E0*2*I_ERI_Fy2z_S_Dxz_S_c;
  abcd[1769] = 4.0E0*I_ERI_F3z_S_Gx3z_S_cc-2.0E0*1*I_ERI_F3z_S_Dxz_S_c-2.0E0*2*I_ERI_F3z_S_Dxz_S_c;
  abcd[1770] = 4.0E0*I_ERI_F3x_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3x_S_D2y_S_c;
  abcd[1771] = 4.0E0*I_ERI_F2xy_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xy_S_D2y_S_c;
  abcd[1772] = 4.0E0*I_ERI_F2xz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2xz_S_D2y_S_c;
  abcd[1773] = 4.0E0*I_ERI_Fx2y_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2y_S_D2y_S_c;
  abcd[1774] = 4.0E0*I_ERI_Fxyz_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fxyz_S_D2y_S_c;
  abcd[1775] = 4.0E0*I_ERI_Fx2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fx2z_S_D2y_S_c;
  abcd[1776] = 4.0E0*I_ERI_F3y_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3y_S_D2y_S_c;
  abcd[1777] = 4.0E0*I_ERI_F2yz_S_G2y2z_S_cc-2.0E0*1*I_ERI_F2yz_S_D2y_S_c;
  abcd[1778] = 4.0E0*I_ERI_Fy2z_S_G2y2z_S_cc-2.0E0*1*I_ERI_Fy2z_S_D2y_S_c;
  abcd[1779] = 4.0E0*I_ERI_F3z_S_G2y2z_S_cc-2.0E0*1*I_ERI_F3z_S_D2y_S_c;
  abcd[1780] = 4.0E0*I_ERI_F3x_S_Gy3z_S_cc-2.0E0*1*I_ERI_F3x_S_Dyz_S_c-2.0E0*2*I_ERI_F3x_S_Dyz_S_c;
  abcd[1781] = 4.0E0*I_ERI_F2xy_S_Gy3z_S_cc-2.0E0*1*I_ERI_F2xy_S_Dyz_S_c-2.0E0*2*I_ERI_F2xy_S_Dyz_S_c;
  abcd[1782] = 4.0E0*I_ERI_F2xz_S_Gy3z_S_cc-2.0E0*1*I_ERI_F2xz_S_Dyz_S_c-2.0E0*2*I_ERI_F2xz_S_Dyz_S_c;
  abcd[1783] = 4.0E0*I_ERI_Fx2y_S_Gy3z_S_cc-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_c-2.0E0*2*I_ERI_Fx2y_S_Dyz_S_c;
  abcd[1784] = 4.0E0*I_ERI_Fxyz_S_Gy3z_S_cc-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_c-2.0E0*2*I_ERI_Fxyz_S_Dyz_S_c;
  abcd[1785] = 4.0E0*I_ERI_Fx2z_S_Gy3z_S_cc-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_c-2.0E0*2*I_ERI_Fx2z_S_Dyz_S_c;
  abcd[1786] = 4.0E0*I_ERI_F3y_S_Gy3z_S_cc-2.0E0*1*I_ERI_F3y_S_Dyz_S_c-2.0E0*2*I_ERI_F3y_S_Dyz_S_c;
  abcd[1787] = 4.0E0*I_ERI_F2yz_S_Gy3z_S_cc-2.0E0*1*I_ERI_F2yz_S_Dyz_S_c-2.0E0*2*I_ERI_F2yz_S_Dyz_S_c;
  abcd[1788] = 4.0E0*I_ERI_Fy2z_S_Gy3z_S_cc-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_c-2.0E0*2*I_ERI_Fy2z_S_Dyz_S_c;
  abcd[1789] = 4.0E0*I_ERI_F3z_S_Gy3z_S_cc-2.0E0*1*I_ERI_F3z_S_Dyz_S_c-2.0E0*2*I_ERI_F3z_S_Dyz_S_c;
  abcd[1790] = 4.0E0*I_ERI_F3x_S_G4z_S_cc-2.0E0*2*I_ERI_F3x_S_D2z_S_c-2.0E0*3*I_ERI_F3x_S_D2z_S_c+2*1*I_ERI_F3x_S_S_S;
  abcd[1791] = 4.0E0*I_ERI_F2xy_S_G4z_S_cc-2.0E0*2*I_ERI_F2xy_S_D2z_S_c-2.0E0*3*I_ERI_F2xy_S_D2z_S_c+2*1*I_ERI_F2xy_S_S_S;
  abcd[1792] = 4.0E0*I_ERI_F2xz_S_G4z_S_cc-2.0E0*2*I_ERI_F2xz_S_D2z_S_c-2.0E0*3*I_ERI_F2xz_S_D2z_S_c+2*1*I_ERI_F2xz_S_S_S;
  abcd[1793] = 4.0E0*I_ERI_Fx2y_S_G4z_S_cc-2.0E0*2*I_ERI_Fx2y_S_D2z_S_c-2.0E0*3*I_ERI_Fx2y_S_D2z_S_c+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1794] = 4.0E0*I_ERI_Fxyz_S_G4z_S_cc-2.0E0*2*I_ERI_Fxyz_S_D2z_S_c-2.0E0*3*I_ERI_Fxyz_S_D2z_S_c+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1795] = 4.0E0*I_ERI_Fx2z_S_G4z_S_cc-2.0E0*2*I_ERI_Fx2z_S_D2z_S_c-2.0E0*3*I_ERI_Fx2z_S_D2z_S_c+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1796] = 4.0E0*I_ERI_F3y_S_G4z_S_cc-2.0E0*2*I_ERI_F3y_S_D2z_S_c-2.0E0*3*I_ERI_F3y_S_D2z_S_c+2*1*I_ERI_F3y_S_S_S;
  abcd[1797] = 4.0E0*I_ERI_F2yz_S_G4z_S_cc-2.0E0*2*I_ERI_F2yz_S_D2z_S_c-2.0E0*3*I_ERI_F2yz_S_D2z_S_c+2*1*I_ERI_F2yz_S_S_S;
  abcd[1798] = 4.0E0*I_ERI_Fy2z_S_G4z_S_cc-2.0E0*2*I_ERI_Fy2z_S_D2z_S_c-2.0E0*3*I_ERI_Fy2z_S_D2z_S_c+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1799] = 4.0E0*I_ERI_F3z_S_G4z_S_cc-2.0E0*2*I_ERI_F3z_S_D2z_S_c-2.0E0*3*I_ERI_F3z_S_D2z_S_c+2*1*I_ERI_F3z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[1800] = 4.0E0*I_ERI_F3x_S_F3x_Px_cd-2.0E0*2*I_ERI_F3x_S_Px_Px_d;
  abcd[1801] = 4.0E0*I_ERI_F2xy_S_F3x_Px_cd-2.0E0*2*I_ERI_F2xy_S_Px_Px_d;
  abcd[1802] = 4.0E0*I_ERI_F2xz_S_F3x_Px_cd-2.0E0*2*I_ERI_F2xz_S_Px_Px_d;
  abcd[1803] = 4.0E0*I_ERI_Fx2y_S_F3x_Px_cd-2.0E0*2*I_ERI_Fx2y_S_Px_Px_d;
  abcd[1804] = 4.0E0*I_ERI_Fxyz_S_F3x_Px_cd-2.0E0*2*I_ERI_Fxyz_S_Px_Px_d;
  abcd[1805] = 4.0E0*I_ERI_Fx2z_S_F3x_Px_cd-2.0E0*2*I_ERI_Fx2z_S_Px_Px_d;
  abcd[1806] = 4.0E0*I_ERI_F3y_S_F3x_Px_cd-2.0E0*2*I_ERI_F3y_S_Px_Px_d;
  abcd[1807] = 4.0E0*I_ERI_F2yz_S_F3x_Px_cd-2.0E0*2*I_ERI_F2yz_S_Px_Px_d;
  abcd[1808] = 4.0E0*I_ERI_Fy2z_S_F3x_Px_cd-2.0E0*2*I_ERI_Fy2z_S_Px_Px_d;
  abcd[1809] = 4.0E0*I_ERI_F3z_S_F3x_Px_cd-2.0E0*2*I_ERI_F3z_S_Px_Px_d;
  abcd[1810] = 4.0E0*I_ERI_F3x_S_F2xy_Px_cd-2.0E0*1*I_ERI_F3x_S_Py_Px_d;
  abcd[1811] = 4.0E0*I_ERI_F2xy_S_F2xy_Px_cd-2.0E0*1*I_ERI_F2xy_S_Py_Px_d;
  abcd[1812] = 4.0E0*I_ERI_F2xz_S_F2xy_Px_cd-2.0E0*1*I_ERI_F2xz_S_Py_Px_d;
  abcd[1813] = 4.0E0*I_ERI_Fx2y_S_F2xy_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Px_d;
  abcd[1814] = 4.0E0*I_ERI_Fxyz_S_F2xy_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Px_d;
  abcd[1815] = 4.0E0*I_ERI_Fx2z_S_F2xy_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Px_d;
  abcd[1816] = 4.0E0*I_ERI_F3y_S_F2xy_Px_cd-2.0E0*1*I_ERI_F3y_S_Py_Px_d;
  abcd[1817] = 4.0E0*I_ERI_F2yz_S_F2xy_Px_cd-2.0E0*1*I_ERI_F2yz_S_Py_Px_d;
  abcd[1818] = 4.0E0*I_ERI_Fy2z_S_F2xy_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Px_d;
  abcd[1819] = 4.0E0*I_ERI_F3z_S_F2xy_Px_cd-2.0E0*1*I_ERI_F3z_S_Py_Px_d;
  abcd[1820] = 4.0E0*I_ERI_F3x_S_F2xz_Px_cd-2.0E0*1*I_ERI_F3x_S_Pz_Px_d;
  abcd[1821] = 4.0E0*I_ERI_F2xy_S_F2xz_Px_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Px_d;
  abcd[1822] = 4.0E0*I_ERI_F2xz_S_F2xz_Px_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Px_d;
  abcd[1823] = 4.0E0*I_ERI_Fx2y_S_F2xz_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Px_d;
  abcd[1824] = 4.0E0*I_ERI_Fxyz_S_F2xz_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Px_d;
  abcd[1825] = 4.0E0*I_ERI_Fx2z_S_F2xz_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Px_d;
  abcd[1826] = 4.0E0*I_ERI_F3y_S_F2xz_Px_cd-2.0E0*1*I_ERI_F3y_S_Pz_Px_d;
  abcd[1827] = 4.0E0*I_ERI_F2yz_S_F2xz_Px_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Px_d;
  abcd[1828] = 4.0E0*I_ERI_Fy2z_S_F2xz_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Px_d;
  abcd[1829] = 4.0E0*I_ERI_F3z_S_F2xz_Px_cd-2.0E0*1*I_ERI_F3z_S_Pz_Px_d;
  abcd[1830] = 4.0E0*I_ERI_F3x_S_Fx2y_Px_cd;
  abcd[1831] = 4.0E0*I_ERI_F2xy_S_Fx2y_Px_cd;
  abcd[1832] = 4.0E0*I_ERI_F2xz_S_Fx2y_Px_cd;
  abcd[1833] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Px_cd;
  abcd[1834] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Px_cd;
  abcd[1835] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Px_cd;
  abcd[1836] = 4.0E0*I_ERI_F3y_S_Fx2y_Px_cd;
  abcd[1837] = 4.0E0*I_ERI_F2yz_S_Fx2y_Px_cd;
  abcd[1838] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Px_cd;
  abcd[1839] = 4.0E0*I_ERI_F3z_S_Fx2y_Px_cd;
  abcd[1840] = 4.0E0*I_ERI_F3x_S_Fxyz_Px_cd;
  abcd[1841] = 4.0E0*I_ERI_F2xy_S_Fxyz_Px_cd;
  abcd[1842] = 4.0E0*I_ERI_F2xz_S_Fxyz_Px_cd;
  abcd[1843] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Px_cd;
  abcd[1844] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Px_cd;
  abcd[1845] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Px_cd;
  abcd[1846] = 4.0E0*I_ERI_F3y_S_Fxyz_Px_cd;
  abcd[1847] = 4.0E0*I_ERI_F2yz_S_Fxyz_Px_cd;
  abcd[1848] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Px_cd;
  abcd[1849] = 4.0E0*I_ERI_F3z_S_Fxyz_Px_cd;
  abcd[1850] = 4.0E0*I_ERI_F3x_S_Fx2z_Px_cd;
  abcd[1851] = 4.0E0*I_ERI_F2xy_S_Fx2z_Px_cd;
  abcd[1852] = 4.0E0*I_ERI_F2xz_S_Fx2z_Px_cd;
  abcd[1853] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Px_cd;
  abcd[1854] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Px_cd;
  abcd[1855] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Px_cd;
  abcd[1856] = 4.0E0*I_ERI_F3y_S_Fx2z_Px_cd;
  abcd[1857] = 4.0E0*I_ERI_F2yz_S_Fx2z_Px_cd;
  abcd[1858] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Px_cd;
  abcd[1859] = 4.0E0*I_ERI_F3z_S_Fx2z_Px_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[1860] = 4.0E0*I_ERI_F3x_S_F3x_Py_cd-2.0E0*2*I_ERI_F3x_S_Px_Py_d;
  abcd[1861] = 4.0E0*I_ERI_F2xy_S_F3x_Py_cd-2.0E0*2*I_ERI_F2xy_S_Px_Py_d;
  abcd[1862] = 4.0E0*I_ERI_F2xz_S_F3x_Py_cd-2.0E0*2*I_ERI_F2xz_S_Px_Py_d;
  abcd[1863] = 4.0E0*I_ERI_Fx2y_S_F3x_Py_cd-2.0E0*2*I_ERI_Fx2y_S_Px_Py_d;
  abcd[1864] = 4.0E0*I_ERI_Fxyz_S_F3x_Py_cd-2.0E0*2*I_ERI_Fxyz_S_Px_Py_d;
  abcd[1865] = 4.0E0*I_ERI_Fx2z_S_F3x_Py_cd-2.0E0*2*I_ERI_Fx2z_S_Px_Py_d;
  abcd[1866] = 4.0E0*I_ERI_F3y_S_F3x_Py_cd-2.0E0*2*I_ERI_F3y_S_Px_Py_d;
  abcd[1867] = 4.0E0*I_ERI_F2yz_S_F3x_Py_cd-2.0E0*2*I_ERI_F2yz_S_Px_Py_d;
  abcd[1868] = 4.0E0*I_ERI_Fy2z_S_F3x_Py_cd-2.0E0*2*I_ERI_Fy2z_S_Px_Py_d;
  abcd[1869] = 4.0E0*I_ERI_F3z_S_F3x_Py_cd-2.0E0*2*I_ERI_F3z_S_Px_Py_d;
  abcd[1870] = 4.0E0*I_ERI_F3x_S_F2xy_Py_cd-2.0E0*1*I_ERI_F3x_S_Py_Py_d;
  abcd[1871] = 4.0E0*I_ERI_F2xy_S_F2xy_Py_cd-2.0E0*1*I_ERI_F2xy_S_Py_Py_d;
  abcd[1872] = 4.0E0*I_ERI_F2xz_S_F2xy_Py_cd-2.0E0*1*I_ERI_F2xz_S_Py_Py_d;
  abcd[1873] = 4.0E0*I_ERI_Fx2y_S_F2xy_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Py_d;
  abcd[1874] = 4.0E0*I_ERI_Fxyz_S_F2xy_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Py_d;
  abcd[1875] = 4.0E0*I_ERI_Fx2z_S_F2xy_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Py_d;
  abcd[1876] = 4.0E0*I_ERI_F3y_S_F2xy_Py_cd-2.0E0*1*I_ERI_F3y_S_Py_Py_d;
  abcd[1877] = 4.0E0*I_ERI_F2yz_S_F2xy_Py_cd-2.0E0*1*I_ERI_F2yz_S_Py_Py_d;
  abcd[1878] = 4.0E0*I_ERI_Fy2z_S_F2xy_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Py_d;
  abcd[1879] = 4.0E0*I_ERI_F3z_S_F2xy_Py_cd-2.0E0*1*I_ERI_F3z_S_Py_Py_d;
  abcd[1880] = 4.0E0*I_ERI_F3x_S_F2xz_Py_cd-2.0E0*1*I_ERI_F3x_S_Pz_Py_d;
  abcd[1881] = 4.0E0*I_ERI_F2xy_S_F2xz_Py_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Py_d;
  abcd[1882] = 4.0E0*I_ERI_F2xz_S_F2xz_Py_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Py_d;
  abcd[1883] = 4.0E0*I_ERI_Fx2y_S_F2xz_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Py_d;
  abcd[1884] = 4.0E0*I_ERI_Fxyz_S_F2xz_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Py_d;
  abcd[1885] = 4.0E0*I_ERI_Fx2z_S_F2xz_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Py_d;
  abcd[1886] = 4.0E0*I_ERI_F3y_S_F2xz_Py_cd-2.0E0*1*I_ERI_F3y_S_Pz_Py_d;
  abcd[1887] = 4.0E0*I_ERI_F2yz_S_F2xz_Py_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Py_d;
  abcd[1888] = 4.0E0*I_ERI_Fy2z_S_F2xz_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Py_d;
  abcd[1889] = 4.0E0*I_ERI_F3z_S_F2xz_Py_cd-2.0E0*1*I_ERI_F3z_S_Pz_Py_d;
  abcd[1890] = 4.0E0*I_ERI_F3x_S_Fx2y_Py_cd;
  abcd[1891] = 4.0E0*I_ERI_F2xy_S_Fx2y_Py_cd;
  abcd[1892] = 4.0E0*I_ERI_F2xz_S_Fx2y_Py_cd;
  abcd[1893] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Py_cd;
  abcd[1894] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Py_cd;
  abcd[1895] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Py_cd;
  abcd[1896] = 4.0E0*I_ERI_F3y_S_Fx2y_Py_cd;
  abcd[1897] = 4.0E0*I_ERI_F2yz_S_Fx2y_Py_cd;
  abcd[1898] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Py_cd;
  abcd[1899] = 4.0E0*I_ERI_F3z_S_Fx2y_Py_cd;
  abcd[1900] = 4.0E0*I_ERI_F3x_S_Fxyz_Py_cd;
  abcd[1901] = 4.0E0*I_ERI_F2xy_S_Fxyz_Py_cd;
  abcd[1902] = 4.0E0*I_ERI_F2xz_S_Fxyz_Py_cd;
  abcd[1903] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Py_cd;
  abcd[1904] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Py_cd;
  abcd[1905] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Py_cd;
  abcd[1906] = 4.0E0*I_ERI_F3y_S_Fxyz_Py_cd;
  abcd[1907] = 4.0E0*I_ERI_F2yz_S_Fxyz_Py_cd;
  abcd[1908] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Py_cd;
  abcd[1909] = 4.0E0*I_ERI_F3z_S_Fxyz_Py_cd;
  abcd[1910] = 4.0E0*I_ERI_F3x_S_Fx2z_Py_cd;
  abcd[1911] = 4.0E0*I_ERI_F2xy_S_Fx2z_Py_cd;
  abcd[1912] = 4.0E0*I_ERI_F2xz_S_Fx2z_Py_cd;
  abcd[1913] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Py_cd;
  abcd[1914] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Py_cd;
  abcd[1915] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Py_cd;
  abcd[1916] = 4.0E0*I_ERI_F3y_S_Fx2z_Py_cd;
  abcd[1917] = 4.0E0*I_ERI_F2yz_S_Fx2z_Py_cd;
  abcd[1918] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Py_cd;
  abcd[1919] = 4.0E0*I_ERI_F3z_S_Fx2z_Py_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[1920] = 4.0E0*I_ERI_F3x_S_F3x_Pz_cd-2.0E0*2*I_ERI_F3x_S_Px_Pz_d;
  abcd[1921] = 4.0E0*I_ERI_F2xy_S_F3x_Pz_cd-2.0E0*2*I_ERI_F2xy_S_Px_Pz_d;
  abcd[1922] = 4.0E0*I_ERI_F2xz_S_F3x_Pz_cd-2.0E0*2*I_ERI_F2xz_S_Px_Pz_d;
  abcd[1923] = 4.0E0*I_ERI_Fx2y_S_F3x_Pz_cd-2.0E0*2*I_ERI_Fx2y_S_Px_Pz_d;
  abcd[1924] = 4.0E0*I_ERI_Fxyz_S_F3x_Pz_cd-2.0E0*2*I_ERI_Fxyz_S_Px_Pz_d;
  abcd[1925] = 4.0E0*I_ERI_Fx2z_S_F3x_Pz_cd-2.0E0*2*I_ERI_Fx2z_S_Px_Pz_d;
  abcd[1926] = 4.0E0*I_ERI_F3y_S_F3x_Pz_cd-2.0E0*2*I_ERI_F3y_S_Px_Pz_d;
  abcd[1927] = 4.0E0*I_ERI_F2yz_S_F3x_Pz_cd-2.0E0*2*I_ERI_F2yz_S_Px_Pz_d;
  abcd[1928] = 4.0E0*I_ERI_Fy2z_S_F3x_Pz_cd-2.0E0*2*I_ERI_Fy2z_S_Px_Pz_d;
  abcd[1929] = 4.0E0*I_ERI_F3z_S_F3x_Pz_cd-2.0E0*2*I_ERI_F3z_S_Px_Pz_d;
  abcd[1930] = 4.0E0*I_ERI_F3x_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F3x_S_Py_Pz_d;
  abcd[1931] = 4.0E0*I_ERI_F2xy_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Py_Pz_d;
  abcd[1932] = 4.0E0*I_ERI_F2xz_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Py_Pz_d;
  abcd[1933] = 4.0E0*I_ERI_Fx2y_S_F2xy_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Pz_d;
  abcd[1934] = 4.0E0*I_ERI_Fxyz_S_F2xy_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Pz_d;
  abcd[1935] = 4.0E0*I_ERI_Fx2z_S_F2xy_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Pz_d;
  abcd[1936] = 4.0E0*I_ERI_F3y_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F3y_S_Py_Pz_d;
  abcd[1937] = 4.0E0*I_ERI_F2yz_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Py_Pz_d;
  abcd[1938] = 4.0E0*I_ERI_Fy2z_S_F2xy_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Pz_d;
  abcd[1939] = 4.0E0*I_ERI_F3z_S_F2xy_Pz_cd-2.0E0*1*I_ERI_F3z_S_Py_Pz_d;
  abcd[1940] = 4.0E0*I_ERI_F3x_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F3x_S_Pz_Pz_d;
  abcd[1941] = 4.0E0*I_ERI_F2xy_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Pz_d;
  abcd[1942] = 4.0E0*I_ERI_F2xz_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Pz_d;
  abcd[1943] = 4.0E0*I_ERI_Fx2y_S_F2xz_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Pz_d;
  abcd[1944] = 4.0E0*I_ERI_Fxyz_S_F2xz_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Pz_d;
  abcd[1945] = 4.0E0*I_ERI_Fx2z_S_F2xz_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Pz_d;
  abcd[1946] = 4.0E0*I_ERI_F3y_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F3y_S_Pz_Pz_d;
  abcd[1947] = 4.0E0*I_ERI_F2yz_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Pz_d;
  abcd[1948] = 4.0E0*I_ERI_Fy2z_S_F2xz_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Pz_d;
  abcd[1949] = 4.0E0*I_ERI_F3z_S_F2xz_Pz_cd-2.0E0*1*I_ERI_F3z_S_Pz_Pz_d;
  abcd[1950] = 4.0E0*I_ERI_F3x_S_Fx2y_Pz_cd;
  abcd[1951] = 4.0E0*I_ERI_F2xy_S_Fx2y_Pz_cd;
  abcd[1952] = 4.0E0*I_ERI_F2xz_S_Fx2y_Pz_cd;
  abcd[1953] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Pz_cd;
  abcd[1954] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Pz_cd;
  abcd[1955] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Pz_cd;
  abcd[1956] = 4.0E0*I_ERI_F3y_S_Fx2y_Pz_cd;
  abcd[1957] = 4.0E0*I_ERI_F2yz_S_Fx2y_Pz_cd;
  abcd[1958] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Pz_cd;
  abcd[1959] = 4.0E0*I_ERI_F3z_S_Fx2y_Pz_cd;
  abcd[1960] = 4.0E0*I_ERI_F3x_S_Fxyz_Pz_cd;
  abcd[1961] = 4.0E0*I_ERI_F2xy_S_Fxyz_Pz_cd;
  abcd[1962] = 4.0E0*I_ERI_F2xz_S_Fxyz_Pz_cd;
  abcd[1963] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Pz_cd;
  abcd[1964] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Pz_cd;
  abcd[1965] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Pz_cd;
  abcd[1966] = 4.0E0*I_ERI_F3y_S_Fxyz_Pz_cd;
  abcd[1967] = 4.0E0*I_ERI_F2yz_S_Fxyz_Pz_cd;
  abcd[1968] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Pz_cd;
  abcd[1969] = 4.0E0*I_ERI_F3z_S_Fxyz_Pz_cd;
  abcd[1970] = 4.0E0*I_ERI_F3x_S_Fx2z_Pz_cd;
  abcd[1971] = 4.0E0*I_ERI_F2xy_S_Fx2z_Pz_cd;
  abcd[1972] = 4.0E0*I_ERI_F2xz_S_Fx2z_Pz_cd;
  abcd[1973] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Pz_cd;
  abcd[1974] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Pz_cd;
  abcd[1975] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Pz_cd;
  abcd[1976] = 4.0E0*I_ERI_F3y_S_Fx2z_Pz_cd;
  abcd[1977] = 4.0E0*I_ERI_F2yz_S_Fx2z_Pz_cd;
  abcd[1978] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Pz_cd;
  abcd[1979] = 4.0E0*I_ERI_F3z_S_Fx2z_Pz_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcy_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[1980] = 4.0E0*I_ERI_F3x_S_F2xy_Px_cd;
  abcd[1981] = 4.0E0*I_ERI_F2xy_S_F2xy_Px_cd;
  abcd[1982] = 4.0E0*I_ERI_F2xz_S_F2xy_Px_cd;
  abcd[1983] = 4.0E0*I_ERI_Fx2y_S_F2xy_Px_cd;
  abcd[1984] = 4.0E0*I_ERI_Fxyz_S_F2xy_Px_cd;
  abcd[1985] = 4.0E0*I_ERI_Fx2z_S_F2xy_Px_cd;
  abcd[1986] = 4.0E0*I_ERI_F3y_S_F2xy_Px_cd;
  abcd[1987] = 4.0E0*I_ERI_F2yz_S_F2xy_Px_cd;
  abcd[1988] = 4.0E0*I_ERI_Fy2z_S_F2xy_Px_cd;
  abcd[1989] = 4.0E0*I_ERI_F3z_S_F2xy_Px_cd;
  abcd[1990] = 4.0E0*I_ERI_F3x_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F3x_S_Px_Px_d;
  abcd[1991] = 4.0E0*I_ERI_F2xy_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F2xy_S_Px_Px_d;
  abcd[1992] = 4.0E0*I_ERI_F2xz_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F2xz_S_Px_Px_d;
  abcd[1993] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Px_d;
  abcd[1994] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Px_d;
  abcd[1995] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Px_d;
  abcd[1996] = 4.0E0*I_ERI_F3y_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F3y_S_Px_Px_d;
  abcd[1997] = 4.0E0*I_ERI_F2yz_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F2yz_S_Px_Px_d;
  abcd[1998] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Px_d;
  abcd[1999] = 4.0E0*I_ERI_F3z_S_Fx2y_Px_cd-2.0E0*1*I_ERI_F3z_S_Px_Px_d;
  abcd[2000] = 4.0E0*I_ERI_F3x_S_Fxyz_Px_cd;
  abcd[2001] = 4.0E0*I_ERI_F2xy_S_Fxyz_Px_cd;
  abcd[2002] = 4.0E0*I_ERI_F2xz_S_Fxyz_Px_cd;
  abcd[2003] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Px_cd;
  abcd[2004] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Px_cd;
  abcd[2005] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Px_cd;
  abcd[2006] = 4.0E0*I_ERI_F3y_S_Fxyz_Px_cd;
  abcd[2007] = 4.0E0*I_ERI_F2yz_S_Fxyz_Px_cd;
  abcd[2008] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Px_cd;
  abcd[2009] = 4.0E0*I_ERI_F3z_S_Fxyz_Px_cd;
  abcd[2010] = 4.0E0*I_ERI_F3x_S_F3y_Px_cd-2.0E0*2*I_ERI_F3x_S_Py_Px_d;
  abcd[2011] = 4.0E0*I_ERI_F2xy_S_F3y_Px_cd-2.0E0*2*I_ERI_F2xy_S_Py_Px_d;
  abcd[2012] = 4.0E0*I_ERI_F2xz_S_F3y_Px_cd-2.0E0*2*I_ERI_F2xz_S_Py_Px_d;
  abcd[2013] = 4.0E0*I_ERI_Fx2y_S_F3y_Px_cd-2.0E0*2*I_ERI_Fx2y_S_Py_Px_d;
  abcd[2014] = 4.0E0*I_ERI_Fxyz_S_F3y_Px_cd-2.0E0*2*I_ERI_Fxyz_S_Py_Px_d;
  abcd[2015] = 4.0E0*I_ERI_Fx2z_S_F3y_Px_cd-2.0E0*2*I_ERI_Fx2z_S_Py_Px_d;
  abcd[2016] = 4.0E0*I_ERI_F3y_S_F3y_Px_cd-2.0E0*2*I_ERI_F3y_S_Py_Px_d;
  abcd[2017] = 4.0E0*I_ERI_F2yz_S_F3y_Px_cd-2.0E0*2*I_ERI_F2yz_S_Py_Px_d;
  abcd[2018] = 4.0E0*I_ERI_Fy2z_S_F3y_Px_cd-2.0E0*2*I_ERI_Fy2z_S_Py_Px_d;
  abcd[2019] = 4.0E0*I_ERI_F3z_S_F3y_Px_cd-2.0E0*2*I_ERI_F3z_S_Py_Px_d;
  abcd[2020] = 4.0E0*I_ERI_F3x_S_F2yz_Px_cd-2.0E0*1*I_ERI_F3x_S_Pz_Px_d;
  abcd[2021] = 4.0E0*I_ERI_F2xy_S_F2yz_Px_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Px_d;
  abcd[2022] = 4.0E0*I_ERI_F2xz_S_F2yz_Px_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Px_d;
  abcd[2023] = 4.0E0*I_ERI_Fx2y_S_F2yz_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Px_d;
  abcd[2024] = 4.0E0*I_ERI_Fxyz_S_F2yz_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Px_d;
  abcd[2025] = 4.0E0*I_ERI_Fx2z_S_F2yz_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Px_d;
  abcd[2026] = 4.0E0*I_ERI_F3y_S_F2yz_Px_cd-2.0E0*1*I_ERI_F3y_S_Pz_Px_d;
  abcd[2027] = 4.0E0*I_ERI_F2yz_S_F2yz_Px_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Px_d;
  abcd[2028] = 4.0E0*I_ERI_Fy2z_S_F2yz_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Px_d;
  abcd[2029] = 4.0E0*I_ERI_F3z_S_F2yz_Px_cd-2.0E0*1*I_ERI_F3z_S_Pz_Px_d;
  abcd[2030] = 4.0E0*I_ERI_F3x_S_Fy2z_Px_cd;
  abcd[2031] = 4.0E0*I_ERI_F2xy_S_Fy2z_Px_cd;
  abcd[2032] = 4.0E0*I_ERI_F2xz_S_Fy2z_Px_cd;
  abcd[2033] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Px_cd;
  abcd[2034] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Px_cd;
  abcd[2035] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Px_cd;
  abcd[2036] = 4.0E0*I_ERI_F3y_S_Fy2z_Px_cd;
  abcd[2037] = 4.0E0*I_ERI_F2yz_S_Fy2z_Px_cd;
  abcd[2038] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Px_cd;
  abcd[2039] = 4.0E0*I_ERI_F3z_S_Fy2z_Px_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[2040] = 4.0E0*I_ERI_F3x_S_F2xy_Py_cd;
  abcd[2041] = 4.0E0*I_ERI_F2xy_S_F2xy_Py_cd;
  abcd[2042] = 4.0E0*I_ERI_F2xz_S_F2xy_Py_cd;
  abcd[2043] = 4.0E0*I_ERI_Fx2y_S_F2xy_Py_cd;
  abcd[2044] = 4.0E0*I_ERI_Fxyz_S_F2xy_Py_cd;
  abcd[2045] = 4.0E0*I_ERI_Fx2z_S_F2xy_Py_cd;
  abcd[2046] = 4.0E0*I_ERI_F3y_S_F2xy_Py_cd;
  abcd[2047] = 4.0E0*I_ERI_F2yz_S_F2xy_Py_cd;
  abcd[2048] = 4.0E0*I_ERI_Fy2z_S_F2xy_Py_cd;
  abcd[2049] = 4.0E0*I_ERI_F3z_S_F2xy_Py_cd;
  abcd[2050] = 4.0E0*I_ERI_F3x_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F3x_S_Px_Py_d;
  abcd[2051] = 4.0E0*I_ERI_F2xy_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F2xy_S_Px_Py_d;
  abcd[2052] = 4.0E0*I_ERI_F2xz_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F2xz_S_Px_Py_d;
  abcd[2053] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Py_d;
  abcd[2054] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Py_d;
  abcd[2055] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Py_d;
  abcd[2056] = 4.0E0*I_ERI_F3y_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F3y_S_Px_Py_d;
  abcd[2057] = 4.0E0*I_ERI_F2yz_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F2yz_S_Px_Py_d;
  abcd[2058] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Py_d;
  abcd[2059] = 4.0E0*I_ERI_F3z_S_Fx2y_Py_cd-2.0E0*1*I_ERI_F3z_S_Px_Py_d;
  abcd[2060] = 4.0E0*I_ERI_F3x_S_Fxyz_Py_cd;
  abcd[2061] = 4.0E0*I_ERI_F2xy_S_Fxyz_Py_cd;
  abcd[2062] = 4.0E0*I_ERI_F2xz_S_Fxyz_Py_cd;
  abcd[2063] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Py_cd;
  abcd[2064] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Py_cd;
  abcd[2065] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Py_cd;
  abcd[2066] = 4.0E0*I_ERI_F3y_S_Fxyz_Py_cd;
  abcd[2067] = 4.0E0*I_ERI_F2yz_S_Fxyz_Py_cd;
  abcd[2068] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Py_cd;
  abcd[2069] = 4.0E0*I_ERI_F3z_S_Fxyz_Py_cd;
  abcd[2070] = 4.0E0*I_ERI_F3x_S_F3y_Py_cd-2.0E0*2*I_ERI_F3x_S_Py_Py_d;
  abcd[2071] = 4.0E0*I_ERI_F2xy_S_F3y_Py_cd-2.0E0*2*I_ERI_F2xy_S_Py_Py_d;
  abcd[2072] = 4.0E0*I_ERI_F2xz_S_F3y_Py_cd-2.0E0*2*I_ERI_F2xz_S_Py_Py_d;
  abcd[2073] = 4.0E0*I_ERI_Fx2y_S_F3y_Py_cd-2.0E0*2*I_ERI_Fx2y_S_Py_Py_d;
  abcd[2074] = 4.0E0*I_ERI_Fxyz_S_F3y_Py_cd-2.0E0*2*I_ERI_Fxyz_S_Py_Py_d;
  abcd[2075] = 4.0E0*I_ERI_Fx2z_S_F3y_Py_cd-2.0E0*2*I_ERI_Fx2z_S_Py_Py_d;
  abcd[2076] = 4.0E0*I_ERI_F3y_S_F3y_Py_cd-2.0E0*2*I_ERI_F3y_S_Py_Py_d;
  abcd[2077] = 4.0E0*I_ERI_F2yz_S_F3y_Py_cd-2.0E0*2*I_ERI_F2yz_S_Py_Py_d;
  abcd[2078] = 4.0E0*I_ERI_Fy2z_S_F3y_Py_cd-2.0E0*2*I_ERI_Fy2z_S_Py_Py_d;
  abcd[2079] = 4.0E0*I_ERI_F3z_S_F3y_Py_cd-2.0E0*2*I_ERI_F3z_S_Py_Py_d;
  abcd[2080] = 4.0E0*I_ERI_F3x_S_F2yz_Py_cd-2.0E0*1*I_ERI_F3x_S_Pz_Py_d;
  abcd[2081] = 4.0E0*I_ERI_F2xy_S_F2yz_Py_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Py_d;
  abcd[2082] = 4.0E0*I_ERI_F2xz_S_F2yz_Py_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Py_d;
  abcd[2083] = 4.0E0*I_ERI_Fx2y_S_F2yz_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Py_d;
  abcd[2084] = 4.0E0*I_ERI_Fxyz_S_F2yz_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Py_d;
  abcd[2085] = 4.0E0*I_ERI_Fx2z_S_F2yz_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Py_d;
  abcd[2086] = 4.0E0*I_ERI_F3y_S_F2yz_Py_cd-2.0E0*1*I_ERI_F3y_S_Pz_Py_d;
  abcd[2087] = 4.0E0*I_ERI_F2yz_S_F2yz_Py_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Py_d;
  abcd[2088] = 4.0E0*I_ERI_Fy2z_S_F2yz_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Py_d;
  abcd[2089] = 4.0E0*I_ERI_F3z_S_F2yz_Py_cd-2.0E0*1*I_ERI_F3z_S_Pz_Py_d;
  abcd[2090] = 4.0E0*I_ERI_F3x_S_Fy2z_Py_cd;
  abcd[2091] = 4.0E0*I_ERI_F2xy_S_Fy2z_Py_cd;
  abcd[2092] = 4.0E0*I_ERI_F2xz_S_Fy2z_Py_cd;
  abcd[2093] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Py_cd;
  abcd[2094] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Py_cd;
  abcd[2095] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Py_cd;
  abcd[2096] = 4.0E0*I_ERI_F3y_S_Fy2z_Py_cd;
  abcd[2097] = 4.0E0*I_ERI_F2yz_S_Fy2z_Py_cd;
  abcd[2098] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Py_cd;
  abcd[2099] = 4.0E0*I_ERI_F3z_S_Fy2z_Py_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[2100] = 4.0E0*I_ERI_F3x_S_F2xy_Pz_cd;
  abcd[2101] = 4.0E0*I_ERI_F2xy_S_F2xy_Pz_cd;
  abcd[2102] = 4.0E0*I_ERI_F2xz_S_F2xy_Pz_cd;
  abcd[2103] = 4.0E0*I_ERI_Fx2y_S_F2xy_Pz_cd;
  abcd[2104] = 4.0E0*I_ERI_Fxyz_S_F2xy_Pz_cd;
  abcd[2105] = 4.0E0*I_ERI_Fx2z_S_F2xy_Pz_cd;
  abcd[2106] = 4.0E0*I_ERI_F3y_S_F2xy_Pz_cd;
  abcd[2107] = 4.0E0*I_ERI_F2yz_S_F2xy_Pz_cd;
  abcd[2108] = 4.0E0*I_ERI_Fy2z_S_F2xy_Pz_cd;
  abcd[2109] = 4.0E0*I_ERI_F3z_S_F2xy_Pz_cd;
  abcd[2110] = 4.0E0*I_ERI_F3x_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F3x_S_Px_Pz_d;
  abcd[2111] = 4.0E0*I_ERI_F2xy_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Px_Pz_d;
  abcd[2112] = 4.0E0*I_ERI_F2xz_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Px_Pz_d;
  abcd[2113] = 4.0E0*I_ERI_Fx2y_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Pz_d;
  abcd[2114] = 4.0E0*I_ERI_Fxyz_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Pz_d;
  abcd[2115] = 4.0E0*I_ERI_Fx2z_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Pz_d;
  abcd[2116] = 4.0E0*I_ERI_F3y_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F3y_S_Px_Pz_d;
  abcd[2117] = 4.0E0*I_ERI_F2yz_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Px_Pz_d;
  abcd[2118] = 4.0E0*I_ERI_Fy2z_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Pz_d;
  abcd[2119] = 4.0E0*I_ERI_F3z_S_Fx2y_Pz_cd-2.0E0*1*I_ERI_F3z_S_Px_Pz_d;
  abcd[2120] = 4.0E0*I_ERI_F3x_S_Fxyz_Pz_cd;
  abcd[2121] = 4.0E0*I_ERI_F2xy_S_Fxyz_Pz_cd;
  abcd[2122] = 4.0E0*I_ERI_F2xz_S_Fxyz_Pz_cd;
  abcd[2123] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Pz_cd;
  abcd[2124] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Pz_cd;
  abcd[2125] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Pz_cd;
  abcd[2126] = 4.0E0*I_ERI_F3y_S_Fxyz_Pz_cd;
  abcd[2127] = 4.0E0*I_ERI_F2yz_S_Fxyz_Pz_cd;
  abcd[2128] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Pz_cd;
  abcd[2129] = 4.0E0*I_ERI_F3z_S_Fxyz_Pz_cd;
  abcd[2130] = 4.0E0*I_ERI_F3x_S_F3y_Pz_cd-2.0E0*2*I_ERI_F3x_S_Py_Pz_d;
  abcd[2131] = 4.0E0*I_ERI_F2xy_S_F3y_Pz_cd-2.0E0*2*I_ERI_F2xy_S_Py_Pz_d;
  abcd[2132] = 4.0E0*I_ERI_F2xz_S_F3y_Pz_cd-2.0E0*2*I_ERI_F2xz_S_Py_Pz_d;
  abcd[2133] = 4.0E0*I_ERI_Fx2y_S_F3y_Pz_cd-2.0E0*2*I_ERI_Fx2y_S_Py_Pz_d;
  abcd[2134] = 4.0E0*I_ERI_Fxyz_S_F3y_Pz_cd-2.0E0*2*I_ERI_Fxyz_S_Py_Pz_d;
  abcd[2135] = 4.0E0*I_ERI_Fx2z_S_F3y_Pz_cd-2.0E0*2*I_ERI_Fx2z_S_Py_Pz_d;
  abcd[2136] = 4.0E0*I_ERI_F3y_S_F3y_Pz_cd-2.0E0*2*I_ERI_F3y_S_Py_Pz_d;
  abcd[2137] = 4.0E0*I_ERI_F2yz_S_F3y_Pz_cd-2.0E0*2*I_ERI_F2yz_S_Py_Pz_d;
  abcd[2138] = 4.0E0*I_ERI_Fy2z_S_F3y_Pz_cd-2.0E0*2*I_ERI_Fy2z_S_Py_Pz_d;
  abcd[2139] = 4.0E0*I_ERI_F3z_S_F3y_Pz_cd-2.0E0*2*I_ERI_F3z_S_Py_Pz_d;
  abcd[2140] = 4.0E0*I_ERI_F3x_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F3x_S_Pz_Pz_d;
  abcd[2141] = 4.0E0*I_ERI_F2xy_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Pz_Pz_d;
  abcd[2142] = 4.0E0*I_ERI_F2xz_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Pz_Pz_d;
  abcd[2143] = 4.0E0*I_ERI_Fx2y_S_F2yz_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Pz_Pz_d;
  abcd[2144] = 4.0E0*I_ERI_Fxyz_S_F2yz_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Pz_Pz_d;
  abcd[2145] = 4.0E0*I_ERI_Fx2z_S_F2yz_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Pz_Pz_d;
  abcd[2146] = 4.0E0*I_ERI_F3y_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F3y_S_Pz_Pz_d;
  abcd[2147] = 4.0E0*I_ERI_F2yz_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Pz_Pz_d;
  abcd[2148] = 4.0E0*I_ERI_Fy2z_S_F2yz_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Pz_Pz_d;
  abcd[2149] = 4.0E0*I_ERI_F3z_S_F2yz_Pz_cd-2.0E0*1*I_ERI_F3z_S_Pz_Pz_d;
  abcd[2150] = 4.0E0*I_ERI_F3x_S_Fy2z_Pz_cd;
  abcd[2151] = 4.0E0*I_ERI_F2xy_S_Fy2z_Pz_cd;
  abcd[2152] = 4.0E0*I_ERI_F2xz_S_Fy2z_Pz_cd;
  abcd[2153] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Pz_cd;
  abcd[2154] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Pz_cd;
  abcd[2155] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Pz_cd;
  abcd[2156] = 4.0E0*I_ERI_F3y_S_Fy2z_Pz_cd;
  abcd[2157] = 4.0E0*I_ERI_F2yz_S_Fy2z_Pz_cd;
  abcd[2158] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Pz_cd;
  abcd[2159] = 4.0E0*I_ERI_F3z_S_Fy2z_Pz_cd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcz_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[2160] = 4.0E0*I_ERI_F3x_S_F2xz_Px_cd;
  abcd[2161] = 4.0E0*I_ERI_F2xy_S_F2xz_Px_cd;
  abcd[2162] = 4.0E0*I_ERI_F2xz_S_F2xz_Px_cd;
  abcd[2163] = 4.0E0*I_ERI_Fx2y_S_F2xz_Px_cd;
  abcd[2164] = 4.0E0*I_ERI_Fxyz_S_F2xz_Px_cd;
  abcd[2165] = 4.0E0*I_ERI_Fx2z_S_F2xz_Px_cd;
  abcd[2166] = 4.0E0*I_ERI_F3y_S_F2xz_Px_cd;
  abcd[2167] = 4.0E0*I_ERI_F2yz_S_F2xz_Px_cd;
  abcd[2168] = 4.0E0*I_ERI_Fy2z_S_F2xz_Px_cd;
  abcd[2169] = 4.0E0*I_ERI_F3z_S_F2xz_Px_cd;
  abcd[2170] = 4.0E0*I_ERI_F3x_S_Fxyz_Px_cd;
  abcd[2171] = 4.0E0*I_ERI_F2xy_S_Fxyz_Px_cd;
  abcd[2172] = 4.0E0*I_ERI_F2xz_S_Fxyz_Px_cd;
  abcd[2173] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Px_cd;
  abcd[2174] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Px_cd;
  abcd[2175] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Px_cd;
  abcd[2176] = 4.0E0*I_ERI_F3y_S_Fxyz_Px_cd;
  abcd[2177] = 4.0E0*I_ERI_F2yz_S_Fxyz_Px_cd;
  abcd[2178] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Px_cd;
  abcd[2179] = 4.0E0*I_ERI_F3z_S_Fxyz_Px_cd;
  abcd[2180] = 4.0E0*I_ERI_F3x_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F3x_S_Px_Px_d;
  abcd[2181] = 4.0E0*I_ERI_F2xy_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F2xy_S_Px_Px_d;
  abcd[2182] = 4.0E0*I_ERI_F2xz_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F2xz_S_Px_Px_d;
  abcd[2183] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Px_d;
  abcd[2184] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Px_d;
  abcd[2185] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Px_d;
  abcd[2186] = 4.0E0*I_ERI_F3y_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F3y_S_Px_Px_d;
  abcd[2187] = 4.0E0*I_ERI_F2yz_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F2yz_S_Px_Px_d;
  abcd[2188] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Px_d;
  abcd[2189] = 4.0E0*I_ERI_F3z_S_Fx2z_Px_cd-2.0E0*1*I_ERI_F3z_S_Px_Px_d;
  abcd[2190] = 4.0E0*I_ERI_F3x_S_F2yz_Px_cd;
  abcd[2191] = 4.0E0*I_ERI_F2xy_S_F2yz_Px_cd;
  abcd[2192] = 4.0E0*I_ERI_F2xz_S_F2yz_Px_cd;
  abcd[2193] = 4.0E0*I_ERI_Fx2y_S_F2yz_Px_cd;
  abcd[2194] = 4.0E0*I_ERI_Fxyz_S_F2yz_Px_cd;
  abcd[2195] = 4.0E0*I_ERI_Fx2z_S_F2yz_Px_cd;
  abcd[2196] = 4.0E0*I_ERI_F3y_S_F2yz_Px_cd;
  abcd[2197] = 4.0E0*I_ERI_F2yz_S_F2yz_Px_cd;
  abcd[2198] = 4.0E0*I_ERI_Fy2z_S_F2yz_Px_cd;
  abcd[2199] = 4.0E0*I_ERI_F3z_S_F2yz_Px_cd;
  abcd[2200] = 4.0E0*I_ERI_F3x_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F3x_S_Py_Px_d;
  abcd[2201] = 4.0E0*I_ERI_F2xy_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F2xy_S_Py_Px_d;
  abcd[2202] = 4.0E0*I_ERI_F2xz_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F2xz_S_Py_Px_d;
  abcd[2203] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Px_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Px_d;
  abcd[2204] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Px_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Px_d;
  abcd[2205] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Px_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Px_d;
  abcd[2206] = 4.0E0*I_ERI_F3y_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F3y_S_Py_Px_d;
  abcd[2207] = 4.0E0*I_ERI_F2yz_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F2yz_S_Py_Px_d;
  abcd[2208] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Px_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Px_d;
  abcd[2209] = 4.0E0*I_ERI_F3z_S_Fy2z_Px_cd-2.0E0*1*I_ERI_F3z_S_Py_Px_d;
  abcd[2210] = 4.0E0*I_ERI_F3x_S_F3z_Px_cd-2.0E0*2*I_ERI_F3x_S_Pz_Px_d;
  abcd[2211] = 4.0E0*I_ERI_F2xy_S_F3z_Px_cd-2.0E0*2*I_ERI_F2xy_S_Pz_Px_d;
  abcd[2212] = 4.0E0*I_ERI_F2xz_S_F3z_Px_cd-2.0E0*2*I_ERI_F2xz_S_Pz_Px_d;
  abcd[2213] = 4.0E0*I_ERI_Fx2y_S_F3z_Px_cd-2.0E0*2*I_ERI_Fx2y_S_Pz_Px_d;
  abcd[2214] = 4.0E0*I_ERI_Fxyz_S_F3z_Px_cd-2.0E0*2*I_ERI_Fxyz_S_Pz_Px_d;
  abcd[2215] = 4.0E0*I_ERI_Fx2z_S_F3z_Px_cd-2.0E0*2*I_ERI_Fx2z_S_Pz_Px_d;
  abcd[2216] = 4.0E0*I_ERI_F3y_S_F3z_Px_cd-2.0E0*2*I_ERI_F3y_S_Pz_Px_d;
  abcd[2217] = 4.0E0*I_ERI_F2yz_S_F3z_Px_cd-2.0E0*2*I_ERI_F2yz_S_Pz_Px_d;
  abcd[2218] = 4.0E0*I_ERI_Fy2z_S_F3z_Px_cd-2.0E0*2*I_ERI_Fy2z_S_Pz_Px_d;
  abcd[2219] = 4.0E0*I_ERI_F3z_S_F3z_Px_cd-2.0E0*2*I_ERI_F3z_S_Pz_Px_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcz_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[2220] = 4.0E0*I_ERI_F3x_S_F2xz_Py_cd;
  abcd[2221] = 4.0E0*I_ERI_F2xy_S_F2xz_Py_cd;
  abcd[2222] = 4.0E0*I_ERI_F2xz_S_F2xz_Py_cd;
  abcd[2223] = 4.0E0*I_ERI_Fx2y_S_F2xz_Py_cd;
  abcd[2224] = 4.0E0*I_ERI_Fxyz_S_F2xz_Py_cd;
  abcd[2225] = 4.0E0*I_ERI_Fx2z_S_F2xz_Py_cd;
  abcd[2226] = 4.0E0*I_ERI_F3y_S_F2xz_Py_cd;
  abcd[2227] = 4.0E0*I_ERI_F2yz_S_F2xz_Py_cd;
  abcd[2228] = 4.0E0*I_ERI_Fy2z_S_F2xz_Py_cd;
  abcd[2229] = 4.0E0*I_ERI_F3z_S_F2xz_Py_cd;
  abcd[2230] = 4.0E0*I_ERI_F3x_S_Fxyz_Py_cd;
  abcd[2231] = 4.0E0*I_ERI_F2xy_S_Fxyz_Py_cd;
  abcd[2232] = 4.0E0*I_ERI_F2xz_S_Fxyz_Py_cd;
  abcd[2233] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Py_cd;
  abcd[2234] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Py_cd;
  abcd[2235] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Py_cd;
  abcd[2236] = 4.0E0*I_ERI_F3y_S_Fxyz_Py_cd;
  abcd[2237] = 4.0E0*I_ERI_F2yz_S_Fxyz_Py_cd;
  abcd[2238] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Py_cd;
  abcd[2239] = 4.0E0*I_ERI_F3z_S_Fxyz_Py_cd;
  abcd[2240] = 4.0E0*I_ERI_F3x_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F3x_S_Px_Py_d;
  abcd[2241] = 4.0E0*I_ERI_F2xy_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F2xy_S_Px_Py_d;
  abcd[2242] = 4.0E0*I_ERI_F2xz_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F2xz_S_Px_Py_d;
  abcd[2243] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Py_d;
  abcd[2244] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Py_d;
  abcd[2245] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Py_d;
  abcd[2246] = 4.0E0*I_ERI_F3y_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F3y_S_Px_Py_d;
  abcd[2247] = 4.0E0*I_ERI_F2yz_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F2yz_S_Px_Py_d;
  abcd[2248] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Py_d;
  abcd[2249] = 4.0E0*I_ERI_F3z_S_Fx2z_Py_cd-2.0E0*1*I_ERI_F3z_S_Px_Py_d;
  abcd[2250] = 4.0E0*I_ERI_F3x_S_F2yz_Py_cd;
  abcd[2251] = 4.0E0*I_ERI_F2xy_S_F2yz_Py_cd;
  abcd[2252] = 4.0E0*I_ERI_F2xz_S_F2yz_Py_cd;
  abcd[2253] = 4.0E0*I_ERI_Fx2y_S_F2yz_Py_cd;
  abcd[2254] = 4.0E0*I_ERI_Fxyz_S_F2yz_Py_cd;
  abcd[2255] = 4.0E0*I_ERI_Fx2z_S_F2yz_Py_cd;
  abcd[2256] = 4.0E0*I_ERI_F3y_S_F2yz_Py_cd;
  abcd[2257] = 4.0E0*I_ERI_F2yz_S_F2yz_Py_cd;
  abcd[2258] = 4.0E0*I_ERI_Fy2z_S_F2yz_Py_cd;
  abcd[2259] = 4.0E0*I_ERI_F3z_S_F2yz_Py_cd;
  abcd[2260] = 4.0E0*I_ERI_F3x_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F3x_S_Py_Py_d;
  abcd[2261] = 4.0E0*I_ERI_F2xy_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F2xy_S_Py_Py_d;
  abcd[2262] = 4.0E0*I_ERI_F2xz_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F2xz_S_Py_Py_d;
  abcd[2263] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Py_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Py_d;
  abcd[2264] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Py_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Py_d;
  abcd[2265] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Py_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Py_d;
  abcd[2266] = 4.0E0*I_ERI_F3y_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F3y_S_Py_Py_d;
  abcd[2267] = 4.0E0*I_ERI_F2yz_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F2yz_S_Py_Py_d;
  abcd[2268] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Py_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Py_d;
  abcd[2269] = 4.0E0*I_ERI_F3z_S_Fy2z_Py_cd-2.0E0*1*I_ERI_F3z_S_Py_Py_d;
  abcd[2270] = 4.0E0*I_ERI_F3x_S_F3z_Py_cd-2.0E0*2*I_ERI_F3x_S_Pz_Py_d;
  abcd[2271] = 4.0E0*I_ERI_F2xy_S_F3z_Py_cd-2.0E0*2*I_ERI_F2xy_S_Pz_Py_d;
  abcd[2272] = 4.0E0*I_ERI_F2xz_S_F3z_Py_cd-2.0E0*2*I_ERI_F2xz_S_Pz_Py_d;
  abcd[2273] = 4.0E0*I_ERI_Fx2y_S_F3z_Py_cd-2.0E0*2*I_ERI_Fx2y_S_Pz_Py_d;
  abcd[2274] = 4.0E0*I_ERI_Fxyz_S_F3z_Py_cd-2.0E0*2*I_ERI_Fxyz_S_Pz_Py_d;
  abcd[2275] = 4.0E0*I_ERI_Fx2z_S_F3z_Py_cd-2.0E0*2*I_ERI_Fx2z_S_Pz_Py_d;
  abcd[2276] = 4.0E0*I_ERI_F3y_S_F3z_Py_cd-2.0E0*2*I_ERI_F3y_S_Pz_Py_d;
  abcd[2277] = 4.0E0*I_ERI_F2yz_S_F3z_Py_cd-2.0E0*2*I_ERI_F2yz_S_Pz_Py_d;
  abcd[2278] = 4.0E0*I_ERI_Fy2z_S_F3z_Py_cd-2.0E0*2*I_ERI_Fy2z_S_Pz_Py_d;
  abcd[2279] = 4.0E0*I_ERI_F3z_S_F3z_Py_cd-2.0E0*2*I_ERI_F3z_S_Pz_Py_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_dcz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_cd
   * RHS shell quartet name: SQ_ERI_F_S_P_P_d
   ************************************************************/
  abcd[2280] = 4.0E0*I_ERI_F3x_S_F2xz_Pz_cd;
  abcd[2281] = 4.0E0*I_ERI_F2xy_S_F2xz_Pz_cd;
  abcd[2282] = 4.0E0*I_ERI_F2xz_S_F2xz_Pz_cd;
  abcd[2283] = 4.0E0*I_ERI_Fx2y_S_F2xz_Pz_cd;
  abcd[2284] = 4.0E0*I_ERI_Fxyz_S_F2xz_Pz_cd;
  abcd[2285] = 4.0E0*I_ERI_Fx2z_S_F2xz_Pz_cd;
  abcd[2286] = 4.0E0*I_ERI_F3y_S_F2xz_Pz_cd;
  abcd[2287] = 4.0E0*I_ERI_F2yz_S_F2xz_Pz_cd;
  abcd[2288] = 4.0E0*I_ERI_Fy2z_S_F2xz_Pz_cd;
  abcd[2289] = 4.0E0*I_ERI_F3z_S_F2xz_Pz_cd;
  abcd[2290] = 4.0E0*I_ERI_F3x_S_Fxyz_Pz_cd;
  abcd[2291] = 4.0E0*I_ERI_F2xy_S_Fxyz_Pz_cd;
  abcd[2292] = 4.0E0*I_ERI_F2xz_S_Fxyz_Pz_cd;
  abcd[2293] = 4.0E0*I_ERI_Fx2y_S_Fxyz_Pz_cd;
  abcd[2294] = 4.0E0*I_ERI_Fxyz_S_Fxyz_Pz_cd;
  abcd[2295] = 4.0E0*I_ERI_Fx2z_S_Fxyz_Pz_cd;
  abcd[2296] = 4.0E0*I_ERI_F3y_S_Fxyz_Pz_cd;
  abcd[2297] = 4.0E0*I_ERI_F2yz_S_Fxyz_Pz_cd;
  abcd[2298] = 4.0E0*I_ERI_Fy2z_S_Fxyz_Pz_cd;
  abcd[2299] = 4.0E0*I_ERI_F3z_S_Fxyz_Pz_cd;
  abcd[2300] = 4.0E0*I_ERI_F3x_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F3x_S_Px_Pz_d;
  abcd[2301] = 4.0E0*I_ERI_F2xy_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Px_Pz_d;
  abcd[2302] = 4.0E0*I_ERI_F2xz_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Px_Pz_d;
  abcd[2303] = 4.0E0*I_ERI_Fx2y_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Px_Pz_d;
  abcd[2304] = 4.0E0*I_ERI_Fxyz_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Px_Pz_d;
  abcd[2305] = 4.0E0*I_ERI_Fx2z_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Px_Pz_d;
  abcd[2306] = 4.0E0*I_ERI_F3y_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F3y_S_Px_Pz_d;
  abcd[2307] = 4.0E0*I_ERI_F2yz_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Px_Pz_d;
  abcd[2308] = 4.0E0*I_ERI_Fy2z_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Px_Pz_d;
  abcd[2309] = 4.0E0*I_ERI_F3z_S_Fx2z_Pz_cd-2.0E0*1*I_ERI_F3z_S_Px_Pz_d;
  abcd[2310] = 4.0E0*I_ERI_F3x_S_F2yz_Pz_cd;
  abcd[2311] = 4.0E0*I_ERI_F2xy_S_F2yz_Pz_cd;
  abcd[2312] = 4.0E0*I_ERI_F2xz_S_F2yz_Pz_cd;
  abcd[2313] = 4.0E0*I_ERI_Fx2y_S_F2yz_Pz_cd;
  abcd[2314] = 4.0E0*I_ERI_Fxyz_S_F2yz_Pz_cd;
  abcd[2315] = 4.0E0*I_ERI_Fx2z_S_F2yz_Pz_cd;
  abcd[2316] = 4.0E0*I_ERI_F3y_S_F2yz_Pz_cd;
  abcd[2317] = 4.0E0*I_ERI_F2yz_S_F2yz_Pz_cd;
  abcd[2318] = 4.0E0*I_ERI_Fy2z_S_F2yz_Pz_cd;
  abcd[2319] = 4.0E0*I_ERI_F3z_S_F2yz_Pz_cd;
  abcd[2320] = 4.0E0*I_ERI_F3x_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F3x_S_Py_Pz_d;
  abcd[2321] = 4.0E0*I_ERI_F2xy_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F2xy_S_Py_Pz_d;
  abcd[2322] = 4.0E0*I_ERI_F2xz_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F2xz_S_Py_Pz_d;
  abcd[2323] = 4.0E0*I_ERI_Fx2y_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_Fx2y_S_Py_Pz_d;
  abcd[2324] = 4.0E0*I_ERI_Fxyz_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_Fxyz_S_Py_Pz_d;
  abcd[2325] = 4.0E0*I_ERI_Fx2z_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_Fx2z_S_Py_Pz_d;
  abcd[2326] = 4.0E0*I_ERI_F3y_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F3y_S_Py_Pz_d;
  abcd[2327] = 4.0E0*I_ERI_F2yz_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F2yz_S_Py_Pz_d;
  abcd[2328] = 4.0E0*I_ERI_Fy2z_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_Fy2z_S_Py_Pz_d;
  abcd[2329] = 4.0E0*I_ERI_F3z_S_Fy2z_Pz_cd-2.0E0*1*I_ERI_F3z_S_Py_Pz_d;
  abcd[2330] = 4.0E0*I_ERI_F3x_S_F3z_Pz_cd-2.0E0*2*I_ERI_F3x_S_Pz_Pz_d;
  abcd[2331] = 4.0E0*I_ERI_F2xy_S_F3z_Pz_cd-2.0E0*2*I_ERI_F2xy_S_Pz_Pz_d;
  abcd[2332] = 4.0E0*I_ERI_F2xz_S_F3z_Pz_cd-2.0E0*2*I_ERI_F2xz_S_Pz_Pz_d;
  abcd[2333] = 4.0E0*I_ERI_Fx2y_S_F3z_Pz_cd-2.0E0*2*I_ERI_Fx2y_S_Pz_Pz_d;
  abcd[2334] = 4.0E0*I_ERI_Fxyz_S_F3z_Pz_cd-2.0E0*2*I_ERI_Fxyz_S_Pz_Pz_d;
  abcd[2335] = 4.0E0*I_ERI_Fx2z_S_F3z_Pz_cd-2.0E0*2*I_ERI_Fx2z_S_Pz_Pz_d;
  abcd[2336] = 4.0E0*I_ERI_F3y_S_F3z_Pz_cd-2.0E0*2*I_ERI_F3y_S_Pz_Pz_d;
  abcd[2337] = 4.0E0*I_ERI_F2yz_S_F3z_Pz_cd-2.0E0*2*I_ERI_F2yz_S_Pz_Pz_d;
  abcd[2338] = 4.0E0*I_ERI_Fy2z_S_F3z_Pz_cd-2.0E0*2*I_ERI_Fy2z_S_Pz_Pz_d;
  abcd[2339] = 4.0E0*I_ERI_F3z_S_F3z_Pz_cd-2.0E0*2*I_ERI_F3z_S_Pz_Pz_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddx_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2340] = 4.0E0*I_ERI_F3x_S_D2x_D2x_dd-2.0E0*1*I_ERI_F3x_S_D2x_S_d;
  abcd[2341] = 4.0E0*I_ERI_F2xy_S_D2x_D2x_dd-2.0E0*1*I_ERI_F2xy_S_D2x_S_d;
  abcd[2342] = 4.0E0*I_ERI_F2xz_S_D2x_D2x_dd-2.0E0*1*I_ERI_F2xz_S_D2x_S_d;
  abcd[2343] = 4.0E0*I_ERI_Fx2y_S_D2x_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_D2x_S_d;
  abcd[2344] = 4.0E0*I_ERI_Fxyz_S_D2x_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_D2x_S_d;
  abcd[2345] = 4.0E0*I_ERI_Fx2z_S_D2x_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_D2x_S_d;
  abcd[2346] = 4.0E0*I_ERI_F3y_S_D2x_D2x_dd-2.0E0*1*I_ERI_F3y_S_D2x_S_d;
  abcd[2347] = 4.0E0*I_ERI_F2yz_S_D2x_D2x_dd-2.0E0*1*I_ERI_F2yz_S_D2x_S_d;
  abcd[2348] = 4.0E0*I_ERI_Fy2z_S_D2x_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_D2x_S_d;
  abcd[2349] = 4.0E0*I_ERI_F3z_S_D2x_D2x_dd-2.0E0*1*I_ERI_F3z_S_D2x_S_d;
  abcd[2350] = 4.0E0*I_ERI_F3x_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F3x_S_Dxy_S_d;
  abcd[2351] = 4.0E0*I_ERI_F2xy_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Dxy_S_d;
  abcd[2352] = 4.0E0*I_ERI_F2xz_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Dxy_S_d;
  abcd[2353] = 4.0E0*I_ERI_Fx2y_S_Dxy_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_d;
  abcd[2354] = 4.0E0*I_ERI_Fxyz_S_Dxy_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_d;
  abcd[2355] = 4.0E0*I_ERI_Fx2z_S_Dxy_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_d;
  abcd[2356] = 4.0E0*I_ERI_F3y_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F3y_S_Dxy_S_d;
  abcd[2357] = 4.0E0*I_ERI_F2yz_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Dxy_S_d;
  abcd[2358] = 4.0E0*I_ERI_Fy2z_S_Dxy_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_d;
  abcd[2359] = 4.0E0*I_ERI_F3z_S_Dxy_D2x_dd-2.0E0*1*I_ERI_F3z_S_Dxy_S_d;
  abcd[2360] = 4.0E0*I_ERI_F3x_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F3x_S_Dxz_S_d;
  abcd[2361] = 4.0E0*I_ERI_F2xy_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Dxz_S_d;
  abcd[2362] = 4.0E0*I_ERI_F2xz_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Dxz_S_d;
  abcd[2363] = 4.0E0*I_ERI_Fx2y_S_Dxz_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_d;
  abcd[2364] = 4.0E0*I_ERI_Fxyz_S_Dxz_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_d;
  abcd[2365] = 4.0E0*I_ERI_Fx2z_S_Dxz_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_d;
  abcd[2366] = 4.0E0*I_ERI_F3y_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F3y_S_Dxz_S_d;
  abcd[2367] = 4.0E0*I_ERI_F2yz_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Dxz_S_d;
  abcd[2368] = 4.0E0*I_ERI_Fy2z_S_Dxz_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_d;
  abcd[2369] = 4.0E0*I_ERI_F3z_S_Dxz_D2x_dd-2.0E0*1*I_ERI_F3z_S_Dxz_S_d;
  abcd[2370] = 4.0E0*I_ERI_F3x_S_D2y_D2x_dd-2.0E0*1*I_ERI_F3x_S_D2y_S_d;
  abcd[2371] = 4.0E0*I_ERI_F2xy_S_D2y_D2x_dd-2.0E0*1*I_ERI_F2xy_S_D2y_S_d;
  abcd[2372] = 4.0E0*I_ERI_F2xz_S_D2y_D2x_dd-2.0E0*1*I_ERI_F2xz_S_D2y_S_d;
  abcd[2373] = 4.0E0*I_ERI_Fx2y_S_D2y_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_D2y_S_d;
  abcd[2374] = 4.0E0*I_ERI_Fxyz_S_D2y_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_D2y_S_d;
  abcd[2375] = 4.0E0*I_ERI_Fx2z_S_D2y_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_D2y_S_d;
  abcd[2376] = 4.0E0*I_ERI_F3y_S_D2y_D2x_dd-2.0E0*1*I_ERI_F3y_S_D2y_S_d;
  abcd[2377] = 4.0E0*I_ERI_F2yz_S_D2y_D2x_dd-2.0E0*1*I_ERI_F2yz_S_D2y_S_d;
  abcd[2378] = 4.0E0*I_ERI_Fy2z_S_D2y_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_D2y_S_d;
  abcd[2379] = 4.0E0*I_ERI_F3z_S_D2y_D2x_dd-2.0E0*1*I_ERI_F3z_S_D2y_S_d;
  abcd[2380] = 4.0E0*I_ERI_F3x_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F3x_S_Dyz_S_d;
  abcd[2381] = 4.0E0*I_ERI_F2xy_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F2xy_S_Dyz_S_d;
  abcd[2382] = 4.0E0*I_ERI_F2xz_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F2xz_S_Dyz_S_d;
  abcd[2383] = 4.0E0*I_ERI_Fx2y_S_Dyz_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_d;
  abcd[2384] = 4.0E0*I_ERI_Fxyz_S_Dyz_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_d;
  abcd[2385] = 4.0E0*I_ERI_Fx2z_S_Dyz_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_d;
  abcd[2386] = 4.0E0*I_ERI_F3y_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F3y_S_Dyz_S_d;
  abcd[2387] = 4.0E0*I_ERI_F2yz_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F2yz_S_Dyz_S_d;
  abcd[2388] = 4.0E0*I_ERI_Fy2z_S_Dyz_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_d;
  abcd[2389] = 4.0E0*I_ERI_F3z_S_Dyz_D2x_dd-2.0E0*1*I_ERI_F3z_S_Dyz_S_d;
  abcd[2390] = 4.0E0*I_ERI_F3x_S_D2z_D2x_dd-2.0E0*1*I_ERI_F3x_S_D2z_S_d;
  abcd[2391] = 4.0E0*I_ERI_F2xy_S_D2z_D2x_dd-2.0E0*1*I_ERI_F2xy_S_D2z_S_d;
  abcd[2392] = 4.0E0*I_ERI_F2xz_S_D2z_D2x_dd-2.0E0*1*I_ERI_F2xz_S_D2z_S_d;
  abcd[2393] = 4.0E0*I_ERI_Fx2y_S_D2z_D2x_dd-2.0E0*1*I_ERI_Fx2y_S_D2z_S_d;
  abcd[2394] = 4.0E0*I_ERI_Fxyz_S_D2z_D2x_dd-2.0E0*1*I_ERI_Fxyz_S_D2z_S_d;
  abcd[2395] = 4.0E0*I_ERI_Fx2z_S_D2z_D2x_dd-2.0E0*1*I_ERI_Fx2z_S_D2z_S_d;
  abcd[2396] = 4.0E0*I_ERI_F3y_S_D2z_D2x_dd-2.0E0*1*I_ERI_F3y_S_D2z_S_d;
  abcd[2397] = 4.0E0*I_ERI_F2yz_S_D2z_D2x_dd-2.0E0*1*I_ERI_F2yz_S_D2z_S_d;
  abcd[2398] = 4.0E0*I_ERI_Fy2z_S_D2z_D2x_dd-2.0E0*1*I_ERI_Fy2z_S_D2z_S_d;
  abcd[2399] = 4.0E0*I_ERI_F3z_S_D2z_D2x_dd-2.0E0*1*I_ERI_F3z_S_D2z_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddx_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2400] = 4.0E0*I_ERI_F3x_S_D2x_Dxy_dd;
  abcd[2401] = 4.0E0*I_ERI_F2xy_S_D2x_Dxy_dd;
  abcd[2402] = 4.0E0*I_ERI_F2xz_S_D2x_Dxy_dd;
  abcd[2403] = 4.0E0*I_ERI_Fx2y_S_D2x_Dxy_dd;
  abcd[2404] = 4.0E0*I_ERI_Fxyz_S_D2x_Dxy_dd;
  abcd[2405] = 4.0E0*I_ERI_Fx2z_S_D2x_Dxy_dd;
  abcd[2406] = 4.0E0*I_ERI_F3y_S_D2x_Dxy_dd;
  abcd[2407] = 4.0E0*I_ERI_F2yz_S_D2x_Dxy_dd;
  abcd[2408] = 4.0E0*I_ERI_Fy2z_S_D2x_Dxy_dd;
  abcd[2409] = 4.0E0*I_ERI_F3z_S_D2x_Dxy_dd;
  abcd[2410] = 4.0E0*I_ERI_F3x_S_Dxy_Dxy_dd;
  abcd[2411] = 4.0E0*I_ERI_F2xy_S_Dxy_Dxy_dd;
  abcd[2412] = 4.0E0*I_ERI_F2xz_S_Dxy_Dxy_dd;
  abcd[2413] = 4.0E0*I_ERI_Fx2y_S_Dxy_Dxy_dd;
  abcd[2414] = 4.0E0*I_ERI_Fxyz_S_Dxy_Dxy_dd;
  abcd[2415] = 4.0E0*I_ERI_Fx2z_S_Dxy_Dxy_dd;
  abcd[2416] = 4.0E0*I_ERI_F3y_S_Dxy_Dxy_dd;
  abcd[2417] = 4.0E0*I_ERI_F2yz_S_Dxy_Dxy_dd;
  abcd[2418] = 4.0E0*I_ERI_Fy2z_S_Dxy_Dxy_dd;
  abcd[2419] = 4.0E0*I_ERI_F3z_S_Dxy_Dxy_dd;
  abcd[2420] = 4.0E0*I_ERI_F3x_S_Dxz_Dxy_dd;
  abcd[2421] = 4.0E0*I_ERI_F2xy_S_Dxz_Dxy_dd;
  abcd[2422] = 4.0E0*I_ERI_F2xz_S_Dxz_Dxy_dd;
  abcd[2423] = 4.0E0*I_ERI_Fx2y_S_Dxz_Dxy_dd;
  abcd[2424] = 4.0E0*I_ERI_Fxyz_S_Dxz_Dxy_dd;
  abcd[2425] = 4.0E0*I_ERI_Fx2z_S_Dxz_Dxy_dd;
  abcd[2426] = 4.0E0*I_ERI_F3y_S_Dxz_Dxy_dd;
  abcd[2427] = 4.0E0*I_ERI_F2yz_S_Dxz_Dxy_dd;
  abcd[2428] = 4.0E0*I_ERI_Fy2z_S_Dxz_Dxy_dd;
  abcd[2429] = 4.0E0*I_ERI_F3z_S_Dxz_Dxy_dd;
  abcd[2430] = 4.0E0*I_ERI_F3x_S_D2y_Dxy_dd;
  abcd[2431] = 4.0E0*I_ERI_F2xy_S_D2y_Dxy_dd;
  abcd[2432] = 4.0E0*I_ERI_F2xz_S_D2y_Dxy_dd;
  abcd[2433] = 4.0E0*I_ERI_Fx2y_S_D2y_Dxy_dd;
  abcd[2434] = 4.0E0*I_ERI_Fxyz_S_D2y_Dxy_dd;
  abcd[2435] = 4.0E0*I_ERI_Fx2z_S_D2y_Dxy_dd;
  abcd[2436] = 4.0E0*I_ERI_F3y_S_D2y_Dxy_dd;
  abcd[2437] = 4.0E0*I_ERI_F2yz_S_D2y_Dxy_dd;
  abcd[2438] = 4.0E0*I_ERI_Fy2z_S_D2y_Dxy_dd;
  abcd[2439] = 4.0E0*I_ERI_F3z_S_D2y_Dxy_dd;
  abcd[2440] = 4.0E0*I_ERI_F3x_S_Dyz_Dxy_dd;
  abcd[2441] = 4.0E0*I_ERI_F2xy_S_Dyz_Dxy_dd;
  abcd[2442] = 4.0E0*I_ERI_F2xz_S_Dyz_Dxy_dd;
  abcd[2443] = 4.0E0*I_ERI_Fx2y_S_Dyz_Dxy_dd;
  abcd[2444] = 4.0E0*I_ERI_Fxyz_S_Dyz_Dxy_dd;
  abcd[2445] = 4.0E0*I_ERI_Fx2z_S_Dyz_Dxy_dd;
  abcd[2446] = 4.0E0*I_ERI_F3y_S_Dyz_Dxy_dd;
  abcd[2447] = 4.0E0*I_ERI_F2yz_S_Dyz_Dxy_dd;
  abcd[2448] = 4.0E0*I_ERI_Fy2z_S_Dyz_Dxy_dd;
  abcd[2449] = 4.0E0*I_ERI_F3z_S_Dyz_Dxy_dd;
  abcd[2450] = 4.0E0*I_ERI_F3x_S_D2z_Dxy_dd;
  abcd[2451] = 4.0E0*I_ERI_F2xy_S_D2z_Dxy_dd;
  abcd[2452] = 4.0E0*I_ERI_F2xz_S_D2z_Dxy_dd;
  abcd[2453] = 4.0E0*I_ERI_Fx2y_S_D2z_Dxy_dd;
  abcd[2454] = 4.0E0*I_ERI_Fxyz_S_D2z_Dxy_dd;
  abcd[2455] = 4.0E0*I_ERI_Fx2z_S_D2z_Dxy_dd;
  abcd[2456] = 4.0E0*I_ERI_F3y_S_D2z_Dxy_dd;
  abcd[2457] = 4.0E0*I_ERI_F2yz_S_D2z_Dxy_dd;
  abcd[2458] = 4.0E0*I_ERI_Fy2z_S_D2z_Dxy_dd;
  abcd[2459] = 4.0E0*I_ERI_F3z_S_D2z_Dxy_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddx_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2460] = 4.0E0*I_ERI_F3x_S_D2x_Dxz_dd;
  abcd[2461] = 4.0E0*I_ERI_F2xy_S_D2x_Dxz_dd;
  abcd[2462] = 4.0E0*I_ERI_F2xz_S_D2x_Dxz_dd;
  abcd[2463] = 4.0E0*I_ERI_Fx2y_S_D2x_Dxz_dd;
  abcd[2464] = 4.0E0*I_ERI_Fxyz_S_D2x_Dxz_dd;
  abcd[2465] = 4.0E0*I_ERI_Fx2z_S_D2x_Dxz_dd;
  abcd[2466] = 4.0E0*I_ERI_F3y_S_D2x_Dxz_dd;
  abcd[2467] = 4.0E0*I_ERI_F2yz_S_D2x_Dxz_dd;
  abcd[2468] = 4.0E0*I_ERI_Fy2z_S_D2x_Dxz_dd;
  abcd[2469] = 4.0E0*I_ERI_F3z_S_D2x_Dxz_dd;
  abcd[2470] = 4.0E0*I_ERI_F3x_S_Dxy_Dxz_dd;
  abcd[2471] = 4.0E0*I_ERI_F2xy_S_Dxy_Dxz_dd;
  abcd[2472] = 4.0E0*I_ERI_F2xz_S_Dxy_Dxz_dd;
  abcd[2473] = 4.0E0*I_ERI_Fx2y_S_Dxy_Dxz_dd;
  abcd[2474] = 4.0E0*I_ERI_Fxyz_S_Dxy_Dxz_dd;
  abcd[2475] = 4.0E0*I_ERI_Fx2z_S_Dxy_Dxz_dd;
  abcd[2476] = 4.0E0*I_ERI_F3y_S_Dxy_Dxz_dd;
  abcd[2477] = 4.0E0*I_ERI_F2yz_S_Dxy_Dxz_dd;
  abcd[2478] = 4.0E0*I_ERI_Fy2z_S_Dxy_Dxz_dd;
  abcd[2479] = 4.0E0*I_ERI_F3z_S_Dxy_Dxz_dd;
  abcd[2480] = 4.0E0*I_ERI_F3x_S_Dxz_Dxz_dd;
  abcd[2481] = 4.0E0*I_ERI_F2xy_S_Dxz_Dxz_dd;
  abcd[2482] = 4.0E0*I_ERI_F2xz_S_Dxz_Dxz_dd;
  abcd[2483] = 4.0E0*I_ERI_Fx2y_S_Dxz_Dxz_dd;
  abcd[2484] = 4.0E0*I_ERI_Fxyz_S_Dxz_Dxz_dd;
  abcd[2485] = 4.0E0*I_ERI_Fx2z_S_Dxz_Dxz_dd;
  abcd[2486] = 4.0E0*I_ERI_F3y_S_Dxz_Dxz_dd;
  abcd[2487] = 4.0E0*I_ERI_F2yz_S_Dxz_Dxz_dd;
  abcd[2488] = 4.0E0*I_ERI_Fy2z_S_Dxz_Dxz_dd;
  abcd[2489] = 4.0E0*I_ERI_F3z_S_Dxz_Dxz_dd;
  abcd[2490] = 4.0E0*I_ERI_F3x_S_D2y_Dxz_dd;
  abcd[2491] = 4.0E0*I_ERI_F2xy_S_D2y_Dxz_dd;
  abcd[2492] = 4.0E0*I_ERI_F2xz_S_D2y_Dxz_dd;
  abcd[2493] = 4.0E0*I_ERI_Fx2y_S_D2y_Dxz_dd;
  abcd[2494] = 4.0E0*I_ERI_Fxyz_S_D2y_Dxz_dd;
  abcd[2495] = 4.0E0*I_ERI_Fx2z_S_D2y_Dxz_dd;
  abcd[2496] = 4.0E0*I_ERI_F3y_S_D2y_Dxz_dd;
  abcd[2497] = 4.0E0*I_ERI_F2yz_S_D2y_Dxz_dd;
  abcd[2498] = 4.0E0*I_ERI_Fy2z_S_D2y_Dxz_dd;
  abcd[2499] = 4.0E0*I_ERI_F3z_S_D2y_Dxz_dd;
  abcd[2500] = 4.0E0*I_ERI_F3x_S_Dyz_Dxz_dd;
  abcd[2501] = 4.0E0*I_ERI_F2xy_S_Dyz_Dxz_dd;
  abcd[2502] = 4.0E0*I_ERI_F2xz_S_Dyz_Dxz_dd;
  abcd[2503] = 4.0E0*I_ERI_Fx2y_S_Dyz_Dxz_dd;
  abcd[2504] = 4.0E0*I_ERI_Fxyz_S_Dyz_Dxz_dd;
  abcd[2505] = 4.0E0*I_ERI_Fx2z_S_Dyz_Dxz_dd;
  abcd[2506] = 4.0E0*I_ERI_F3y_S_Dyz_Dxz_dd;
  abcd[2507] = 4.0E0*I_ERI_F2yz_S_Dyz_Dxz_dd;
  abcd[2508] = 4.0E0*I_ERI_Fy2z_S_Dyz_Dxz_dd;
  abcd[2509] = 4.0E0*I_ERI_F3z_S_Dyz_Dxz_dd;
  abcd[2510] = 4.0E0*I_ERI_F3x_S_D2z_Dxz_dd;
  abcd[2511] = 4.0E0*I_ERI_F2xy_S_D2z_Dxz_dd;
  abcd[2512] = 4.0E0*I_ERI_F2xz_S_D2z_Dxz_dd;
  abcd[2513] = 4.0E0*I_ERI_Fx2y_S_D2z_Dxz_dd;
  abcd[2514] = 4.0E0*I_ERI_Fxyz_S_D2z_Dxz_dd;
  abcd[2515] = 4.0E0*I_ERI_Fx2z_S_D2z_Dxz_dd;
  abcd[2516] = 4.0E0*I_ERI_F3y_S_D2z_Dxz_dd;
  abcd[2517] = 4.0E0*I_ERI_F2yz_S_D2z_Dxz_dd;
  abcd[2518] = 4.0E0*I_ERI_Fy2z_S_D2z_Dxz_dd;
  abcd[2519] = 4.0E0*I_ERI_F3z_S_D2z_Dxz_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddy_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2520] = 4.0E0*I_ERI_F3x_S_D2x_D2y_dd-2.0E0*1*I_ERI_F3x_S_D2x_S_d;
  abcd[2521] = 4.0E0*I_ERI_F2xy_S_D2x_D2y_dd-2.0E0*1*I_ERI_F2xy_S_D2x_S_d;
  abcd[2522] = 4.0E0*I_ERI_F2xz_S_D2x_D2y_dd-2.0E0*1*I_ERI_F2xz_S_D2x_S_d;
  abcd[2523] = 4.0E0*I_ERI_Fx2y_S_D2x_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_D2x_S_d;
  abcd[2524] = 4.0E0*I_ERI_Fxyz_S_D2x_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_D2x_S_d;
  abcd[2525] = 4.0E0*I_ERI_Fx2z_S_D2x_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_D2x_S_d;
  abcd[2526] = 4.0E0*I_ERI_F3y_S_D2x_D2y_dd-2.0E0*1*I_ERI_F3y_S_D2x_S_d;
  abcd[2527] = 4.0E0*I_ERI_F2yz_S_D2x_D2y_dd-2.0E0*1*I_ERI_F2yz_S_D2x_S_d;
  abcd[2528] = 4.0E0*I_ERI_Fy2z_S_D2x_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_D2x_S_d;
  abcd[2529] = 4.0E0*I_ERI_F3z_S_D2x_D2y_dd-2.0E0*1*I_ERI_F3z_S_D2x_S_d;
  abcd[2530] = 4.0E0*I_ERI_F3x_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F3x_S_Dxy_S_d;
  abcd[2531] = 4.0E0*I_ERI_F2xy_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Dxy_S_d;
  abcd[2532] = 4.0E0*I_ERI_F2xz_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Dxy_S_d;
  abcd[2533] = 4.0E0*I_ERI_Fx2y_S_Dxy_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_d;
  abcd[2534] = 4.0E0*I_ERI_Fxyz_S_Dxy_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_d;
  abcd[2535] = 4.0E0*I_ERI_Fx2z_S_Dxy_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_d;
  abcd[2536] = 4.0E0*I_ERI_F3y_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F3y_S_Dxy_S_d;
  abcd[2537] = 4.0E0*I_ERI_F2yz_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Dxy_S_d;
  abcd[2538] = 4.0E0*I_ERI_Fy2z_S_Dxy_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_d;
  abcd[2539] = 4.0E0*I_ERI_F3z_S_Dxy_D2y_dd-2.0E0*1*I_ERI_F3z_S_Dxy_S_d;
  abcd[2540] = 4.0E0*I_ERI_F3x_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F3x_S_Dxz_S_d;
  abcd[2541] = 4.0E0*I_ERI_F2xy_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Dxz_S_d;
  abcd[2542] = 4.0E0*I_ERI_F2xz_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Dxz_S_d;
  abcd[2543] = 4.0E0*I_ERI_Fx2y_S_Dxz_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_d;
  abcd[2544] = 4.0E0*I_ERI_Fxyz_S_Dxz_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_d;
  abcd[2545] = 4.0E0*I_ERI_Fx2z_S_Dxz_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_d;
  abcd[2546] = 4.0E0*I_ERI_F3y_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F3y_S_Dxz_S_d;
  abcd[2547] = 4.0E0*I_ERI_F2yz_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Dxz_S_d;
  abcd[2548] = 4.0E0*I_ERI_Fy2z_S_Dxz_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_d;
  abcd[2549] = 4.0E0*I_ERI_F3z_S_Dxz_D2y_dd-2.0E0*1*I_ERI_F3z_S_Dxz_S_d;
  abcd[2550] = 4.0E0*I_ERI_F3x_S_D2y_D2y_dd-2.0E0*1*I_ERI_F3x_S_D2y_S_d;
  abcd[2551] = 4.0E0*I_ERI_F2xy_S_D2y_D2y_dd-2.0E0*1*I_ERI_F2xy_S_D2y_S_d;
  abcd[2552] = 4.0E0*I_ERI_F2xz_S_D2y_D2y_dd-2.0E0*1*I_ERI_F2xz_S_D2y_S_d;
  abcd[2553] = 4.0E0*I_ERI_Fx2y_S_D2y_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_D2y_S_d;
  abcd[2554] = 4.0E0*I_ERI_Fxyz_S_D2y_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_D2y_S_d;
  abcd[2555] = 4.0E0*I_ERI_Fx2z_S_D2y_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_D2y_S_d;
  abcd[2556] = 4.0E0*I_ERI_F3y_S_D2y_D2y_dd-2.0E0*1*I_ERI_F3y_S_D2y_S_d;
  abcd[2557] = 4.0E0*I_ERI_F2yz_S_D2y_D2y_dd-2.0E0*1*I_ERI_F2yz_S_D2y_S_d;
  abcd[2558] = 4.0E0*I_ERI_Fy2z_S_D2y_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_D2y_S_d;
  abcd[2559] = 4.0E0*I_ERI_F3z_S_D2y_D2y_dd-2.0E0*1*I_ERI_F3z_S_D2y_S_d;
  abcd[2560] = 4.0E0*I_ERI_F3x_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F3x_S_Dyz_S_d;
  abcd[2561] = 4.0E0*I_ERI_F2xy_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F2xy_S_Dyz_S_d;
  abcd[2562] = 4.0E0*I_ERI_F2xz_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F2xz_S_Dyz_S_d;
  abcd[2563] = 4.0E0*I_ERI_Fx2y_S_Dyz_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_d;
  abcd[2564] = 4.0E0*I_ERI_Fxyz_S_Dyz_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_d;
  abcd[2565] = 4.0E0*I_ERI_Fx2z_S_Dyz_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_d;
  abcd[2566] = 4.0E0*I_ERI_F3y_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F3y_S_Dyz_S_d;
  abcd[2567] = 4.0E0*I_ERI_F2yz_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F2yz_S_Dyz_S_d;
  abcd[2568] = 4.0E0*I_ERI_Fy2z_S_Dyz_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_d;
  abcd[2569] = 4.0E0*I_ERI_F3z_S_Dyz_D2y_dd-2.0E0*1*I_ERI_F3z_S_Dyz_S_d;
  abcd[2570] = 4.0E0*I_ERI_F3x_S_D2z_D2y_dd-2.0E0*1*I_ERI_F3x_S_D2z_S_d;
  abcd[2571] = 4.0E0*I_ERI_F2xy_S_D2z_D2y_dd-2.0E0*1*I_ERI_F2xy_S_D2z_S_d;
  abcd[2572] = 4.0E0*I_ERI_F2xz_S_D2z_D2y_dd-2.0E0*1*I_ERI_F2xz_S_D2z_S_d;
  abcd[2573] = 4.0E0*I_ERI_Fx2y_S_D2z_D2y_dd-2.0E0*1*I_ERI_Fx2y_S_D2z_S_d;
  abcd[2574] = 4.0E0*I_ERI_Fxyz_S_D2z_D2y_dd-2.0E0*1*I_ERI_Fxyz_S_D2z_S_d;
  abcd[2575] = 4.0E0*I_ERI_Fx2z_S_D2z_D2y_dd-2.0E0*1*I_ERI_Fx2z_S_D2z_S_d;
  abcd[2576] = 4.0E0*I_ERI_F3y_S_D2z_D2y_dd-2.0E0*1*I_ERI_F3y_S_D2z_S_d;
  abcd[2577] = 4.0E0*I_ERI_F2yz_S_D2z_D2y_dd-2.0E0*1*I_ERI_F2yz_S_D2z_S_d;
  abcd[2578] = 4.0E0*I_ERI_Fy2z_S_D2z_D2y_dd-2.0E0*1*I_ERI_Fy2z_S_D2z_S_d;
  abcd[2579] = 4.0E0*I_ERI_F3z_S_D2z_D2y_dd-2.0E0*1*I_ERI_F3z_S_D2z_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddy_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2580] = 4.0E0*I_ERI_F3x_S_D2x_Dyz_dd;
  abcd[2581] = 4.0E0*I_ERI_F2xy_S_D2x_Dyz_dd;
  abcd[2582] = 4.0E0*I_ERI_F2xz_S_D2x_Dyz_dd;
  abcd[2583] = 4.0E0*I_ERI_Fx2y_S_D2x_Dyz_dd;
  abcd[2584] = 4.0E0*I_ERI_Fxyz_S_D2x_Dyz_dd;
  abcd[2585] = 4.0E0*I_ERI_Fx2z_S_D2x_Dyz_dd;
  abcd[2586] = 4.0E0*I_ERI_F3y_S_D2x_Dyz_dd;
  abcd[2587] = 4.0E0*I_ERI_F2yz_S_D2x_Dyz_dd;
  abcd[2588] = 4.0E0*I_ERI_Fy2z_S_D2x_Dyz_dd;
  abcd[2589] = 4.0E0*I_ERI_F3z_S_D2x_Dyz_dd;
  abcd[2590] = 4.0E0*I_ERI_F3x_S_Dxy_Dyz_dd;
  abcd[2591] = 4.0E0*I_ERI_F2xy_S_Dxy_Dyz_dd;
  abcd[2592] = 4.0E0*I_ERI_F2xz_S_Dxy_Dyz_dd;
  abcd[2593] = 4.0E0*I_ERI_Fx2y_S_Dxy_Dyz_dd;
  abcd[2594] = 4.0E0*I_ERI_Fxyz_S_Dxy_Dyz_dd;
  abcd[2595] = 4.0E0*I_ERI_Fx2z_S_Dxy_Dyz_dd;
  abcd[2596] = 4.0E0*I_ERI_F3y_S_Dxy_Dyz_dd;
  abcd[2597] = 4.0E0*I_ERI_F2yz_S_Dxy_Dyz_dd;
  abcd[2598] = 4.0E0*I_ERI_Fy2z_S_Dxy_Dyz_dd;
  abcd[2599] = 4.0E0*I_ERI_F3z_S_Dxy_Dyz_dd;
  abcd[2600] = 4.0E0*I_ERI_F3x_S_Dxz_Dyz_dd;
  abcd[2601] = 4.0E0*I_ERI_F2xy_S_Dxz_Dyz_dd;
  abcd[2602] = 4.0E0*I_ERI_F2xz_S_Dxz_Dyz_dd;
  abcd[2603] = 4.0E0*I_ERI_Fx2y_S_Dxz_Dyz_dd;
  abcd[2604] = 4.0E0*I_ERI_Fxyz_S_Dxz_Dyz_dd;
  abcd[2605] = 4.0E0*I_ERI_Fx2z_S_Dxz_Dyz_dd;
  abcd[2606] = 4.0E0*I_ERI_F3y_S_Dxz_Dyz_dd;
  abcd[2607] = 4.0E0*I_ERI_F2yz_S_Dxz_Dyz_dd;
  abcd[2608] = 4.0E0*I_ERI_Fy2z_S_Dxz_Dyz_dd;
  abcd[2609] = 4.0E0*I_ERI_F3z_S_Dxz_Dyz_dd;
  abcd[2610] = 4.0E0*I_ERI_F3x_S_D2y_Dyz_dd;
  abcd[2611] = 4.0E0*I_ERI_F2xy_S_D2y_Dyz_dd;
  abcd[2612] = 4.0E0*I_ERI_F2xz_S_D2y_Dyz_dd;
  abcd[2613] = 4.0E0*I_ERI_Fx2y_S_D2y_Dyz_dd;
  abcd[2614] = 4.0E0*I_ERI_Fxyz_S_D2y_Dyz_dd;
  abcd[2615] = 4.0E0*I_ERI_Fx2z_S_D2y_Dyz_dd;
  abcd[2616] = 4.0E0*I_ERI_F3y_S_D2y_Dyz_dd;
  abcd[2617] = 4.0E0*I_ERI_F2yz_S_D2y_Dyz_dd;
  abcd[2618] = 4.0E0*I_ERI_Fy2z_S_D2y_Dyz_dd;
  abcd[2619] = 4.0E0*I_ERI_F3z_S_D2y_Dyz_dd;
  abcd[2620] = 4.0E0*I_ERI_F3x_S_Dyz_Dyz_dd;
  abcd[2621] = 4.0E0*I_ERI_F2xy_S_Dyz_Dyz_dd;
  abcd[2622] = 4.0E0*I_ERI_F2xz_S_Dyz_Dyz_dd;
  abcd[2623] = 4.0E0*I_ERI_Fx2y_S_Dyz_Dyz_dd;
  abcd[2624] = 4.0E0*I_ERI_Fxyz_S_Dyz_Dyz_dd;
  abcd[2625] = 4.0E0*I_ERI_Fx2z_S_Dyz_Dyz_dd;
  abcd[2626] = 4.0E0*I_ERI_F3y_S_Dyz_Dyz_dd;
  abcd[2627] = 4.0E0*I_ERI_F2yz_S_Dyz_Dyz_dd;
  abcd[2628] = 4.0E0*I_ERI_Fy2z_S_Dyz_Dyz_dd;
  abcd[2629] = 4.0E0*I_ERI_F3z_S_Dyz_Dyz_dd;
  abcd[2630] = 4.0E0*I_ERI_F3x_S_D2z_Dyz_dd;
  abcd[2631] = 4.0E0*I_ERI_F2xy_S_D2z_Dyz_dd;
  abcd[2632] = 4.0E0*I_ERI_F2xz_S_D2z_Dyz_dd;
  abcd[2633] = 4.0E0*I_ERI_Fx2y_S_D2z_Dyz_dd;
  abcd[2634] = 4.0E0*I_ERI_Fxyz_S_D2z_Dyz_dd;
  abcd[2635] = 4.0E0*I_ERI_Fx2z_S_D2z_Dyz_dd;
  abcd[2636] = 4.0E0*I_ERI_F3y_S_D2z_Dyz_dd;
  abcd[2637] = 4.0E0*I_ERI_F2yz_S_D2z_Dyz_dd;
  abcd[2638] = 4.0E0*I_ERI_Fy2z_S_D2z_Dyz_dd;
  abcd[2639] = 4.0E0*I_ERI_F3z_S_D2z_Dyz_dd;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_S_ddz_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_D_dd
   * RHS shell quartet name: SQ_ERI_F_S_D_S_d
   ************************************************************/
  abcd[2640] = 4.0E0*I_ERI_F3x_S_D2x_D2z_dd-2.0E0*1*I_ERI_F3x_S_D2x_S_d;
  abcd[2641] = 4.0E0*I_ERI_F2xy_S_D2x_D2z_dd-2.0E0*1*I_ERI_F2xy_S_D2x_S_d;
  abcd[2642] = 4.0E0*I_ERI_F2xz_S_D2x_D2z_dd-2.0E0*1*I_ERI_F2xz_S_D2x_S_d;
  abcd[2643] = 4.0E0*I_ERI_Fx2y_S_D2x_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_D2x_S_d;
  abcd[2644] = 4.0E0*I_ERI_Fxyz_S_D2x_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_D2x_S_d;
  abcd[2645] = 4.0E0*I_ERI_Fx2z_S_D2x_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_D2x_S_d;
  abcd[2646] = 4.0E0*I_ERI_F3y_S_D2x_D2z_dd-2.0E0*1*I_ERI_F3y_S_D2x_S_d;
  abcd[2647] = 4.0E0*I_ERI_F2yz_S_D2x_D2z_dd-2.0E0*1*I_ERI_F2yz_S_D2x_S_d;
  abcd[2648] = 4.0E0*I_ERI_Fy2z_S_D2x_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_D2x_S_d;
  abcd[2649] = 4.0E0*I_ERI_F3z_S_D2x_D2z_dd-2.0E0*1*I_ERI_F3z_S_D2x_S_d;
  abcd[2650] = 4.0E0*I_ERI_F3x_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F3x_S_Dxy_S_d;
  abcd[2651] = 4.0E0*I_ERI_F2xy_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Dxy_S_d;
  abcd[2652] = 4.0E0*I_ERI_F2xz_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Dxy_S_d;
  abcd[2653] = 4.0E0*I_ERI_Fx2y_S_Dxy_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Dxy_S_d;
  abcd[2654] = 4.0E0*I_ERI_Fxyz_S_Dxy_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Dxy_S_d;
  abcd[2655] = 4.0E0*I_ERI_Fx2z_S_Dxy_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Dxy_S_d;
  abcd[2656] = 4.0E0*I_ERI_F3y_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F3y_S_Dxy_S_d;
  abcd[2657] = 4.0E0*I_ERI_F2yz_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Dxy_S_d;
  abcd[2658] = 4.0E0*I_ERI_Fy2z_S_Dxy_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Dxy_S_d;
  abcd[2659] = 4.0E0*I_ERI_F3z_S_Dxy_D2z_dd-2.0E0*1*I_ERI_F3z_S_Dxy_S_d;
  abcd[2660] = 4.0E0*I_ERI_F3x_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F3x_S_Dxz_S_d;
  abcd[2661] = 4.0E0*I_ERI_F2xy_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Dxz_S_d;
  abcd[2662] = 4.0E0*I_ERI_F2xz_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Dxz_S_d;
  abcd[2663] = 4.0E0*I_ERI_Fx2y_S_Dxz_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Dxz_S_d;
  abcd[2664] = 4.0E0*I_ERI_Fxyz_S_Dxz_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Dxz_S_d;
  abcd[2665] = 4.0E0*I_ERI_Fx2z_S_Dxz_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Dxz_S_d;
  abcd[2666] = 4.0E0*I_ERI_F3y_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F3y_S_Dxz_S_d;
  abcd[2667] = 4.0E0*I_ERI_F2yz_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Dxz_S_d;
  abcd[2668] = 4.0E0*I_ERI_Fy2z_S_Dxz_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Dxz_S_d;
  abcd[2669] = 4.0E0*I_ERI_F3z_S_Dxz_D2z_dd-2.0E0*1*I_ERI_F3z_S_Dxz_S_d;
  abcd[2670] = 4.0E0*I_ERI_F3x_S_D2y_D2z_dd-2.0E0*1*I_ERI_F3x_S_D2y_S_d;
  abcd[2671] = 4.0E0*I_ERI_F2xy_S_D2y_D2z_dd-2.0E0*1*I_ERI_F2xy_S_D2y_S_d;
  abcd[2672] = 4.0E0*I_ERI_F2xz_S_D2y_D2z_dd-2.0E0*1*I_ERI_F2xz_S_D2y_S_d;
  abcd[2673] = 4.0E0*I_ERI_Fx2y_S_D2y_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_D2y_S_d;
  abcd[2674] = 4.0E0*I_ERI_Fxyz_S_D2y_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_D2y_S_d;
  abcd[2675] = 4.0E0*I_ERI_Fx2z_S_D2y_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_D2y_S_d;
  abcd[2676] = 4.0E0*I_ERI_F3y_S_D2y_D2z_dd-2.0E0*1*I_ERI_F3y_S_D2y_S_d;
  abcd[2677] = 4.0E0*I_ERI_F2yz_S_D2y_D2z_dd-2.0E0*1*I_ERI_F2yz_S_D2y_S_d;
  abcd[2678] = 4.0E0*I_ERI_Fy2z_S_D2y_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_D2y_S_d;
  abcd[2679] = 4.0E0*I_ERI_F3z_S_D2y_D2z_dd-2.0E0*1*I_ERI_F3z_S_D2y_S_d;
  abcd[2680] = 4.0E0*I_ERI_F3x_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F3x_S_Dyz_S_d;
  abcd[2681] = 4.0E0*I_ERI_F2xy_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F2xy_S_Dyz_S_d;
  abcd[2682] = 4.0E0*I_ERI_F2xz_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F2xz_S_Dyz_S_d;
  abcd[2683] = 4.0E0*I_ERI_Fx2y_S_Dyz_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_Dyz_S_d;
  abcd[2684] = 4.0E0*I_ERI_Fxyz_S_Dyz_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_Dyz_S_d;
  abcd[2685] = 4.0E0*I_ERI_Fx2z_S_Dyz_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_Dyz_S_d;
  abcd[2686] = 4.0E0*I_ERI_F3y_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F3y_S_Dyz_S_d;
  abcd[2687] = 4.0E0*I_ERI_F2yz_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F2yz_S_Dyz_S_d;
  abcd[2688] = 4.0E0*I_ERI_Fy2z_S_Dyz_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_Dyz_S_d;
  abcd[2689] = 4.0E0*I_ERI_F3z_S_Dyz_D2z_dd-2.0E0*1*I_ERI_F3z_S_Dyz_S_d;
  abcd[2690] = 4.0E0*I_ERI_F3x_S_D2z_D2z_dd-2.0E0*1*I_ERI_F3x_S_D2z_S_d;
  abcd[2691] = 4.0E0*I_ERI_F2xy_S_D2z_D2z_dd-2.0E0*1*I_ERI_F2xy_S_D2z_S_d;
  abcd[2692] = 4.0E0*I_ERI_F2xz_S_D2z_D2z_dd-2.0E0*1*I_ERI_F2xz_S_D2z_S_d;
  abcd[2693] = 4.0E0*I_ERI_Fx2y_S_D2z_D2z_dd-2.0E0*1*I_ERI_Fx2y_S_D2z_S_d;
  abcd[2694] = 4.0E0*I_ERI_Fxyz_S_D2z_D2z_dd-2.0E0*1*I_ERI_Fxyz_S_D2z_S_d;
  abcd[2695] = 4.0E0*I_ERI_Fx2z_S_D2z_D2z_dd-2.0E0*1*I_ERI_Fx2z_S_D2z_S_d;
  abcd[2696] = 4.0E0*I_ERI_F3y_S_D2z_D2z_dd-2.0E0*1*I_ERI_F3y_S_D2z_S_d;
  abcd[2697] = 4.0E0*I_ERI_F2yz_S_D2z_D2z_dd-2.0E0*1*I_ERI_F2yz_S_D2z_S_d;
  abcd[2698] = 4.0E0*I_ERI_Fy2z_S_D2z_D2z_dd-2.0E0*1*I_ERI_Fy2z_S_D2z_S_d;
  abcd[2699] = 4.0E0*I_ERI_F3z_S_D2z_D2z_dd-2.0E0*1*I_ERI_F3z_S_D2z_S_d;
}
