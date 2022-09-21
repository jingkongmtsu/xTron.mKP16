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
// BRA1 as redundant position, total RHS integrals evaluated as: 43715
// BRA2 as redundant position, total RHS integrals evaluated as: 43005
// KET1 as redundant position, total RHS integrals evaluated as: 41912
// KET2 as redundant position, total RHS integrals evaluated as: 42725
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
// KET2
// X
// Y
// Z
// ####

void hgp_os_eri_d_p_f_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_D2x_S_F3x_S = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S = 0.0E0;
  Double I_ERI_D2y_S_F3x_S = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S = 0.0E0;
  Double I_ERI_D2z_S_F3x_S = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2x_S_F3y_S = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S = 0.0E0;
  Double I_ERI_D2y_S_F3y_S = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S = 0.0E0;
  Double I_ERI_D2z_S_F3y_S = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2x_S_F3z_S = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S = 0.0E0;
  Double I_ERI_D2y_S_F3z_S = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S = 0.0E0;
  Double I_ERI_D2z_S_F3z_S = 0.0E0;
  Double I_ERI_G4x_S_F3x_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_a = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_a = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_a = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_a = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_a = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_a = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_a = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_a = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_a = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_a = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_a = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_a = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_a = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_a = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_a = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_a = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_a = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_a = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_a = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_a = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_a = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_a = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_a = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_a = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_a = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_a = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_a = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_a = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_a = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_a = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_a = 0.0E0;
  Double I_ERI_Px_S_F3x_S = 0.0E0;
  Double I_ERI_Py_S_F3x_S = 0.0E0;
  Double I_ERI_Pz_S_F3x_S = 0.0E0;
  Double I_ERI_Px_S_F2xy_S = 0.0E0;
  Double I_ERI_Py_S_F2xy_S = 0.0E0;
  Double I_ERI_Pz_S_F2xy_S = 0.0E0;
  Double I_ERI_Px_S_F2xz_S = 0.0E0;
  Double I_ERI_Py_S_F2xz_S = 0.0E0;
  Double I_ERI_Pz_S_F2xz_S = 0.0E0;
  Double I_ERI_Px_S_Fx2y_S = 0.0E0;
  Double I_ERI_Py_S_Fx2y_S = 0.0E0;
  Double I_ERI_Pz_S_Fx2y_S = 0.0E0;
  Double I_ERI_Px_S_Fxyz_S = 0.0E0;
  Double I_ERI_Py_S_Fxyz_S = 0.0E0;
  Double I_ERI_Pz_S_Fxyz_S = 0.0E0;
  Double I_ERI_Px_S_Fx2z_S = 0.0E0;
  Double I_ERI_Py_S_Fx2z_S = 0.0E0;
  Double I_ERI_Pz_S_Fx2z_S = 0.0E0;
  Double I_ERI_Px_S_F3y_S = 0.0E0;
  Double I_ERI_Py_S_F3y_S = 0.0E0;
  Double I_ERI_Pz_S_F3y_S = 0.0E0;
  Double I_ERI_Px_S_F2yz_S = 0.0E0;
  Double I_ERI_Py_S_F2yz_S = 0.0E0;
  Double I_ERI_Pz_S_F2yz_S = 0.0E0;
  Double I_ERI_Px_S_Fy2z_S = 0.0E0;
  Double I_ERI_Py_S_Fy2z_S = 0.0E0;
  Double I_ERI_Pz_S_Fy2z_S = 0.0E0;
  Double I_ERI_Px_S_F3z_S = 0.0E0;
  Double I_ERI_Py_S_F3z_S = 0.0E0;
  Double I_ERI_Pz_S_F3z_S = 0.0E0;
  Double I_ERI_G4x_S_F3x_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_b = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_b = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_b = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_b = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_b = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_b = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_b = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_b = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_b = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_b = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_b = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_b = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_b = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_b = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_b = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_b = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_b = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_b = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_b = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_b = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_b = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_b = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_b = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_b = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_b = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_b = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_b = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_b = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_b = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_b = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_b = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_b = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_b = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_b = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_b = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_b = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_b = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_b = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_b = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_b = 0.0E0;
  Double I_ERI_F3x_S_G4x_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G4x_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G4x_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G4x_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G4x_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G4x_S_d = 0.0E0;
  Double I_ERI_F3y_S_G4x_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G4x_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G4x_S_d = 0.0E0;
  Double I_ERI_F3z_S_G4x_S_d = 0.0E0;
  Double I_ERI_F3x_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F3y_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F3z_S_G3xy_S_d = 0.0E0;
  Double I_ERI_F3x_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F3y_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F3z_S_G3xz_S_d = 0.0E0;
  Double I_ERI_F3x_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F3y_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F3z_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_F3x_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F3y_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F3z_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_F3x_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_F3x_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F3y_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F3z_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_F3x_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F3y_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F3z_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_F3x_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_F3x_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F3y_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F3z_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_F3x_S_G4y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G4y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G4y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G4y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G4y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G4y_S_d = 0.0E0;
  Double I_ERI_F3y_S_G4y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G4y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G4y_S_d = 0.0E0;
  Double I_ERI_F3z_S_G4y_S_d = 0.0E0;
  Double I_ERI_F3x_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F3y_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F3z_S_G3yz_S_d = 0.0E0;
  Double I_ERI_F3x_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_F3x_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F3y_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F3z_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_F3x_S_G4z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_G4z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_G4z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_G4z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_G4z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_G4z_S_d = 0.0E0;
  Double I_ERI_F3y_S_G4z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_G4z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_G4z_S_d = 0.0E0;
  Double I_ERI_F3z_S_G4z_S_d = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_d = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_d = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_d = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_d = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_d = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_d = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_d = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_d = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_d = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_d = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_d = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_d = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_d = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_d = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_d = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_d = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_d = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_d = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_d = 0.0E0;
  Double I_ERI_D2x_S_G4x_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G4x_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G4x_S_d = 0.0E0;
  Double I_ERI_D2y_S_G4x_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G4x_S_d = 0.0E0;
  Double I_ERI_D2z_S_G4x_S_d = 0.0E0;
  Double I_ERI_D2x_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G3xy_S_d = 0.0E0;
  Double I_ERI_D2y_S_G3xy_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G3xy_S_d = 0.0E0;
  Double I_ERI_D2z_S_G3xy_S_d = 0.0E0;
  Double I_ERI_D2x_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G3xz_S_d = 0.0E0;
  Double I_ERI_D2y_S_G3xz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G3xz_S_d = 0.0E0;
  Double I_ERI_D2z_S_G3xz_S_d = 0.0E0;
  Double I_ERI_D2x_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_D2y_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_D2z_S_G2x2y_S_d = 0.0E0;
  Double I_ERI_D2x_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_D2y_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_D2z_S_G2xyz_S_d = 0.0E0;
  Double I_ERI_D2x_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_D2y_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_D2z_S_G2x2z_S_d = 0.0E0;
  Double I_ERI_D2x_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_D2y_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_D2z_S_Gx3y_S_d = 0.0E0;
  Double I_ERI_D2x_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_D2y_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_D2z_S_Gx2yz_S_d = 0.0E0;
  Double I_ERI_D2x_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_D2y_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_D2z_S_Gxy2z_S_d = 0.0E0;
  Double I_ERI_D2x_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_D2y_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_D2z_S_Gx3z_S_d = 0.0E0;
  Double I_ERI_D2x_S_G4y_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G4y_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G4y_S_d = 0.0E0;
  Double I_ERI_D2y_S_G4y_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G4y_S_d = 0.0E0;
  Double I_ERI_D2z_S_G4y_S_d = 0.0E0;
  Double I_ERI_D2x_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G3yz_S_d = 0.0E0;
  Double I_ERI_D2y_S_G3yz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G3yz_S_d = 0.0E0;
  Double I_ERI_D2z_S_G3yz_S_d = 0.0E0;
  Double I_ERI_D2x_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_D2y_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_D2z_S_G2y2z_S_d = 0.0E0;
  Double I_ERI_D2x_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_D2y_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_D2z_S_Gy3z_S_d = 0.0E0;
  Double I_ERI_D2x_S_G4z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_G4z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_G4z_S_d = 0.0E0;
  Double I_ERI_D2y_S_G4z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_G4z_S_d = 0.0E0;
  Double I_ERI_D2z_S_G4z_S_d = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_d = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_d = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_d = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_d = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_d = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_d = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_d = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_d = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_d = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_d = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_d = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_d = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_d = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_d = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_d = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_d = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_d = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_d = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_d = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_d = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_d = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_d = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_d = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_d = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_d = 0.0E0;

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
       * shell quartet name: SQ_ERI_S_S_P_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_Px_S_M6_vrr = QCX*I_ERI_S_S_S_S_M6_vrr+WQX*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Py_S_M6_vrr = QCY*I_ERI_S_S_S_S_M6_vrr+WQY*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Pz_S_M6_vrr = QCZ*I_ERI_S_S_S_S_M6_vrr+WQZ*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_Px_S_M5_vrr = QCX*I_ERI_S_S_S_S_M5_vrr+WQX*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Py_S_M5_vrr = QCY*I_ERI_S_S_S_S_M5_vrr+WQY*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Pz_S_M5_vrr = QCZ*I_ERI_S_S_S_S_M5_vrr+WQZ*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M5_vrr = QCX*I_ERI_S_S_Px_S_M5_vrr+WQX*I_ERI_S_S_Px_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Dxy_S_M5_vrr = QCY*I_ERI_S_S_Px_S_M5_vrr+WQY*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_D2y_S_M5_vrr = QCY*I_ERI_S_S_Py_S_M5_vrr+WQY*I_ERI_S_S_Py_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_D2z_S_M5_vrr = QCZ*I_ERI_S_S_Pz_S_M5_vrr+WQZ*I_ERI_S_S_Pz_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M4_vrr = QCX*I_ERI_S_S_Px_S_M4_vrr+WQX*I_ERI_S_S_Px_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dxy_S_M4_vrr = QCY*I_ERI_S_S_Px_S_M4_vrr+WQY*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_Dxz_S_M4_vrr = QCZ*I_ERI_S_S_Px_S_M4_vrr+WQZ*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_D2y_S_M4_vrr = QCY*I_ERI_S_S_Py_S_M4_vrr+WQY*I_ERI_S_S_Py_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dyz_S_M4_vrr = QCZ*I_ERI_S_S_Py_S_M4_vrr+WQZ*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_D2z_S_M4_vrr = QCZ*I_ERI_S_S_Pz_S_M4_vrr+WQZ*I_ERI_S_S_Pz_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M4_vrr = QCX*I_ERI_S_S_D2x_S_M4_vrr+WQX*I_ERI_S_S_D2x_S_M5_vrr+2*oned2e*I_ERI_S_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_F2xy_S_M4_vrr = QCY*I_ERI_S_S_D2x_S_M4_vrr+WQY*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_F2xz_S_M4_vrr = QCZ*I_ERI_S_S_D2x_S_M4_vrr+WQZ*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_Fx2y_S_M4_vrr = QCX*I_ERI_S_S_D2y_S_M4_vrr+WQX*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fxyz_S_M4_vrr = QCZ*I_ERI_S_S_Dxy_S_M4_vrr+WQZ*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_S_S_Fx2z_S_M4_vrr = QCX*I_ERI_S_S_D2z_S_M4_vrr+WQX*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3y_S_M4_vrr = QCY*I_ERI_S_S_D2y_S_M4_vrr+WQY*I_ERI_S_S_D2y_S_M5_vrr+2*oned2e*I_ERI_S_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_F2yz_S_M4_vrr = QCZ*I_ERI_S_S_D2y_S_M4_vrr+WQZ*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fy2z_S_M4_vrr = QCY*I_ERI_S_S_D2z_S_M4_vrr+WQY*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3z_S_M4_vrr = QCZ*I_ERI_S_S_D2z_S_M4_vrr+WQZ*I_ERI_S_S_D2z_S_M5_vrr+2*oned2e*I_ERI_S_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M3_vrr = PAX*I_ERI_S_S_D2x_S_M3_vrr+WPX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Py_S_D2x_S_M3_vrr = PAY*I_ERI_S_S_D2x_S_M3_vrr+WPY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_D2x_S_M3_vrr = PAZ*I_ERI_S_S_D2x_S_M3_vrr+WPZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Dxy_S_M3_vrr = PAX*I_ERI_S_S_Dxy_S_M3_vrr+WPX*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Dxy_S_M3_vrr = PAY*I_ERI_S_S_Dxy_S_M3_vrr+WPY*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Dxy_S_M3_vrr = PAZ*I_ERI_S_S_Dxy_S_M3_vrr+WPZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Dxz_S_M3_vrr = PAX*I_ERI_S_S_Dxz_S_M3_vrr+WPX*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Dxz_S_M3_vrr = PAY*I_ERI_S_S_Dxz_S_M3_vrr+WPY*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Dxz_S_M3_vrr = PAZ*I_ERI_S_S_Dxz_S_M3_vrr+WPZ*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_D2y_S_M3_vrr = PAX*I_ERI_S_S_D2y_S_M3_vrr+WPX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_D2y_S_M3_vrr = PAY*I_ERI_S_S_D2y_S_M3_vrr+WPY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Pz_S_D2y_S_M3_vrr = PAZ*I_ERI_S_S_D2y_S_M3_vrr+WPZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Dyz_S_M3_vrr = PAX*I_ERI_S_S_Dyz_S_M3_vrr+WPX*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Dyz_S_M3_vrr = PAY*I_ERI_S_S_Dyz_S_M3_vrr+WPY*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Dyz_S_M3_vrr = PAZ*I_ERI_S_S_Dyz_S_M3_vrr+WPZ*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_D2z_S_M3_vrr = PAX*I_ERI_S_S_D2z_S_M3_vrr+WPX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_D2z_S_M3_vrr = PAY*I_ERI_S_S_D2z_S_M3_vrr+WPY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_D2z_S_M3_vrr = PAZ*I_ERI_S_S_D2z_S_M3_vrr+WPZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M3_vrr = QCX*I_ERI_S_S_D2x_S_M3_vrr+WQX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2e*I_ERI_S_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_F2xy_S_M3_vrr = QCY*I_ERI_S_S_D2x_S_M3_vrr+WQY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_F2xz_S_M3_vrr = QCZ*I_ERI_S_S_D2x_S_M3_vrr+WQZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_Fx2y_S_M3_vrr = QCX*I_ERI_S_S_D2y_S_M3_vrr+WQX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fxyz_S_M3_vrr = QCZ*I_ERI_S_S_Dxy_S_M3_vrr+WQZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_S_S_Fx2z_S_M3_vrr = QCX*I_ERI_S_S_D2z_S_M3_vrr+WQX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3y_S_M3_vrr = QCY*I_ERI_S_S_D2y_S_M3_vrr+WQY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2e*I_ERI_S_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_F2yz_S_M3_vrr = QCZ*I_ERI_S_S_D2y_S_M3_vrr+WQZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fy2z_S_M3_vrr = QCY*I_ERI_S_S_D2z_S_M3_vrr+WQY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3z_S_M3_vrr = QCZ*I_ERI_S_S_D2z_S_M3_vrr+WQZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2e*I_ERI_S_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M3_vrr = PAX*I_ERI_S_S_F3x_S_M3_vrr+WPX*I_ERI_S_S_F3x_S_M4_vrr+3*oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F3x_S_M3_vrr = PAY*I_ERI_S_S_F3x_S_M3_vrr+WPY*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Pz_S_F3x_S_M3_vrr = PAZ*I_ERI_S_S_F3x_S_M3_vrr+WPZ*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Px_S_F2xy_S_M3_vrr = PAX*I_ERI_S_S_F2xy_S_M3_vrr+WPX*I_ERI_S_S_F2xy_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Py_S_F2xy_S_M3_vrr = PAY*I_ERI_S_S_F2xy_S_M3_vrr+WPY*I_ERI_S_S_F2xy_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_F2xy_S_M3_vrr = PAZ*I_ERI_S_S_F2xy_S_M3_vrr+WPZ*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Px_S_F2xz_S_M3_vrr = PAX*I_ERI_S_S_F2xz_S_M3_vrr+WPX*I_ERI_S_S_F2xz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Py_S_F2xz_S_M3_vrr = PAY*I_ERI_S_S_F2xz_S_M3_vrr+WPY*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Pz_S_F2xz_S_M3_vrr = PAZ*I_ERI_S_S_F2xz_S_M3_vrr+WPZ*I_ERI_S_S_F2xz_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Fx2y_S_M3_vrr = PAX*I_ERI_S_S_Fx2y_S_M3_vrr+WPX*I_ERI_S_S_Fx2y_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_Fx2y_S_M3_vrr = PAY*I_ERI_S_S_Fx2y_S_M3_vrr+WPY*I_ERI_S_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M3_vrr = PAZ*I_ERI_S_S_Fx2y_S_M3_vrr+WPZ*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Px_S_Fxyz_S_M3_vrr = PAX*I_ERI_S_S_Fxyz_S_M3_vrr+WPX*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Fxyz_S_M3_vrr = PAY*I_ERI_S_S_Fxyz_S_M3_vrr+WPY*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M3_vrr = PAZ*I_ERI_S_S_Fxyz_S_M3_vrr+WPZ*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Fx2z_S_M3_vrr = PAX*I_ERI_S_S_Fx2z_S_M3_vrr+WPX*I_ERI_S_S_Fx2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_Fx2z_S_M3_vrr = PAY*I_ERI_S_S_Fx2z_S_M3_vrr+WPY*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M3_vrr = PAZ*I_ERI_S_S_Fx2z_S_M3_vrr+WPZ*I_ERI_S_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Px_S_F3y_S_M3_vrr = PAX*I_ERI_S_S_F3y_S_M3_vrr+WPX*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Py_S_F3y_S_M3_vrr = PAY*I_ERI_S_S_F3y_S_M3_vrr+WPY*I_ERI_S_S_F3y_S_M4_vrr+3*oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_F3y_S_M3_vrr = PAZ*I_ERI_S_S_F3y_S_M3_vrr+WPZ*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Px_S_F2yz_S_M3_vrr = PAX*I_ERI_S_S_F2yz_S_M3_vrr+WPX*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Py_S_F2yz_S_M3_vrr = PAY*I_ERI_S_S_F2yz_S_M3_vrr+WPY*I_ERI_S_S_F2yz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Pz_S_F2yz_S_M3_vrr = PAZ*I_ERI_S_S_F2yz_S_M3_vrr+WPZ*I_ERI_S_S_F2yz_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Fy2z_S_M3_vrr = PAX*I_ERI_S_S_Fy2z_S_M3_vrr+WPX*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Py_S_Fy2z_S_M3_vrr = PAY*I_ERI_S_S_Fy2z_S_M3_vrr+WPY*I_ERI_S_S_Fy2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M3_vrr = PAZ*I_ERI_S_S_Fy2z_S_M3_vrr+WPZ*I_ERI_S_S_Fy2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Px_S_F3z_S_M3_vrr = PAX*I_ERI_S_S_F3z_S_M3_vrr+WPX*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Py_S_F3z_S_M3_vrr = PAY*I_ERI_S_S_F3z_S_M3_vrr+WPY*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Pz_S_F3z_S_M3_vrr = PAZ*I_ERI_S_S_F3z_S_M3_vrr+WPZ*I_ERI_S_S_F3z_S_M4_vrr+3*oned2k*I_ERI_S_S_D2z_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 8 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M2_vrr = PAX*I_ERI_Px_S_Px_S_M2_vrr+WPX*I_ERI_Px_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr+oned2k*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_Px_S_M2_vrr = PAY*I_ERI_Py_S_Px_S_M2_vrr+WPY*I_ERI_Py_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Px_S_M2_vrr = PAZ*I_ERI_Pz_S_Px_S_M2_vrr+WPZ*I_ERI_Pz_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_Py_S_M2_vrr = PAX*I_ERI_Px_S_Py_S_M2_vrr+WPX*I_ERI_Px_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Py_S_M2_vrr = PAY*I_ERI_Py_S_Py_S_M2_vrr+WPY*I_ERI_Py_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr+oned2k*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_Py_S_M2_vrr = PAZ*I_ERI_Pz_S_Py_S_M2_vrr+WPZ*I_ERI_Pz_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_Pz_S_M2_vrr = PAX*I_ERI_Px_S_Pz_S_M2_vrr+WPX*I_ERI_Px_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Dxy_S_Pz_S_M2_vrr = PAY*I_ERI_Px_S_Pz_S_M2_vrr+WPY*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_D2y_S_Pz_S_M2_vrr = PAY*I_ERI_Py_S_Pz_S_M2_vrr+WPY*I_ERI_Py_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_Pz_S_M2_vrr = PAZ*I_ERI_Pz_S_Pz_S_M2_vrr+WPZ*I_ERI_Pz_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr+oned2k*I_ERI_Pz_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M2_vrr = QCX*I_ERI_S_S_D2x_S_M2_vrr+WQX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2e*I_ERI_S_S_Px_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_F2xy_S_M2_vrr = QCY*I_ERI_S_S_D2x_S_M2_vrr+WQY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_F2xz_S_M2_vrr = QCZ*I_ERI_S_S_D2x_S_M2_vrr+WQZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_Fx2y_S_M2_vrr = QCX*I_ERI_S_S_D2y_S_M2_vrr+WQX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fxyz_S_M2_vrr = QCZ*I_ERI_S_S_Dxy_S_M2_vrr+WQZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_S_S_Fx2z_S_M2_vrr = QCX*I_ERI_S_S_D2z_S_M2_vrr+WQX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3y_S_M2_vrr = QCY*I_ERI_S_S_D2y_S_M2_vrr+WQY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2e*I_ERI_S_S_Py_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_F2yz_S_M2_vrr = QCZ*I_ERI_S_S_D2y_S_M2_vrr+WQZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fy2z_S_M2_vrr = QCY*I_ERI_S_S_D2z_S_M2_vrr+WQY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3z_S_M2_vrr = QCZ*I_ERI_S_S_D2z_S_M2_vrr+WQZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2e*I_ERI_S_S_Pz_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 14 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M2_vrr = PAX*I_ERI_Px_S_D2x_S_M2_vrr+WPX*I_ERI_Px_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Px_S_Px_S_M3_vrr;
      Double I_ERI_Dxy_S_D2x_S_M2_vrr = PAY*I_ERI_Px_S_D2x_S_M2_vrr+WPY*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_D2x_S_M2_vrr = PAY*I_ERI_Py_S_D2x_S_M2_vrr+WPY*I_ERI_Py_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_D2x_S_M2_vrr = PAZ*I_ERI_Pz_S_D2x_S_M2_vrr+WPZ*I_ERI_Pz_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Dxy_S_M2_vrr = PAX*I_ERI_Px_S_Dxy_S_M2_vrr+WPX*I_ERI_Px_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Dxy_S_M2_vrr = PAY*I_ERI_Py_S_Dxy_S_M2_vrr+WPY*I_ERI_Py_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Dxy_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Dxz_S_M2_vrr = PAX*I_ERI_Px_S_Dxz_S_M2_vrr+WPX*I_ERI_Px_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_Dxy_S_Dxz_S_M2_vrr = PAY*I_ERI_Px_S_Dxz_S_M2_vrr+WPY*I_ERI_Px_S_Dxz_S_M3_vrr;
      Double I_ERI_D2y_S_Dxz_S_M2_vrr = PAY*I_ERI_Py_S_Dxz_S_M2_vrr+WPY*I_ERI_Py_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Dxz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_D2y_S_M2_vrr = PAX*I_ERI_Px_S_D2y_S_M2_vrr+WPX*I_ERI_Px_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Dxy_S_D2y_S_M2_vrr = PAY*I_ERI_Px_S_D2y_S_M2_vrr+WPY*I_ERI_Px_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_D2y_S_M2_vrr = PAY*I_ERI_Py_S_D2y_S_M2_vrr+WPY*I_ERI_Py_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Py_S_M3_vrr;
      Double I_ERI_D2z_S_D2y_S_M2_vrr = PAZ*I_ERI_Pz_S_D2y_S_M2_vrr+WPZ*I_ERI_Pz_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Dyz_S_M2_vrr = PAX*I_ERI_Px_S_Dyz_S_M2_vrr+WPX*I_ERI_Px_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Dyz_S_M2_vrr = PAY*I_ERI_Py_S_Dyz_S_M2_vrr+WPY*I_ERI_Py_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Py_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_Dyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_D2z_S_M2_vrr = PAX*I_ERI_Px_S_D2z_S_M2_vrr+WPX*I_ERI_Px_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Dxy_S_D2z_S_M2_vrr = PAY*I_ERI_Px_S_D2z_S_M2_vrr+WPY*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_D2z_S_M2_vrr = PAY*I_ERI_Py_S_D2z_S_M2_vrr+WPY*I_ERI_Py_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_D2z_S_M2_vrr = PAZ*I_ERI_Pz_S_D2z_S_M2_vrr+WPZ*I_ERI_Pz_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M2_vrr = PAX*I_ERI_S_S_F3x_S_M2_vrr+WPX*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Py_S_F3x_S_M2_vrr = PAY*I_ERI_S_S_F3x_S_M2_vrr+WPY*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_F3x_S_M2_vrr = PAZ*I_ERI_S_S_F3x_S_M2_vrr+WPZ*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_F2xy_S_M2_vrr = PAX*I_ERI_S_S_F2xy_S_M2_vrr+WPX*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Py_S_F2xy_S_M2_vrr = PAY*I_ERI_S_S_F2xy_S_M2_vrr+WPY*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_F2xy_S_M2_vrr = PAZ*I_ERI_S_S_F2xy_S_M2_vrr+WPZ*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_F2xz_S_M2_vrr = PAX*I_ERI_S_S_F2xz_S_M2_vrr+WPX*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Py_S_F2xz_S_M2_vrr = PAY*I_ERI_S_S_F2xz_S_M2_vrr+WPY*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Pz_S_F2xz_S_M2_vrr = PAZ*I_ERI_S_S_F2xz_S_M2_vrr+WPZ*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Fx2y_S_M2_vrr = PAX*I_ERI_S_S_Fx2y_S_M2_vrr+WPX*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_Fx2y_S_M2_vrr = PAY*I_ERI_S_S_Fx2y_S_M2_vrr+WPY*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_S_S_Fx2y_S_M2_vrr+WPZ*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Fxyz_S_M2_vrr = PAX*I_ERI_S_S_Fxyz_S_M2_vrr+WPX*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Fxyz_S_M2_vrr = PAY*I_ERI_S_S_Fxyz_S_M2_vrr+WPY*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_S_S_Fxyz_S_M2_vrr+WPZ*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Fx2z_S_M2_vrr = PAX*I_ERI_S_S_Fx2z_S_M2_vrr+WPX*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_Fx2z_S_M2_vrr = PAY*I_ERI_S_S_Fx2z_S_M2_vrr+WPY*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_S_S_Fx2z_S_M2_vrr+WPZ*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Px_S_F3y_S_M2_vrr = PAX*I_ERI_S_S_F3y_S_M2_vrr+WPX*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_F3y_S_M2_vrr = PAY*I_ERI_S_S_F3y_S_M2_vrr+WPY*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Pz_S_F3y_S_M2_vrr = PAZ*I_ERI_S_S_F3y_S_M2_vrr+WPZ*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_F2yz_S_M2_vrr = PAX*I_ERI_S_S_F2yz_S_M2_vrr+WPX*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Py_S_F2yz_S_M2_vrr = PAY*I_ERI_S_S_F2yz_S_M2_vrr+WPY*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Pz_S_F2yz_S_M2_vrr = PAZ*I_ERI_S_S_F2yz_S_M2_vrr+WPZ*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Fy2z_S_M2_vrr = PAX*I_ERI_S_S_Fy2z_S_M2_vrr+WPX*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Py_S_Fy2z_S_M2_vrr = PAY*I_ERI_S_S_Fy2z_S_M2_vrr+WPY*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_S_S_Fy2z_S_M2_vrr+WPZ*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Px_S_F3z_S_M2_vrr = PAX*I_ERI_S_S_F3z_S_M2_vrr+WPX*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_F3z_S_M2_vrr = PAY*I_ERI_S_S_F3z_S_M2_vrr+WPY*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_F3z_S_M2_vrr = PAZ*I_ERI_S_S_F3z_S_M2_vrr+WPZ*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_S_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 22 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M2_vrr = PAX*I_ERI_Px_S_F3x_S_M2_vrr+WPX*I_ERI_Px_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_Dxy_S_F3x_S_M2_vrr = PAY*I_ERI_Px_S_F3x_S_M2_vrr+WPY*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_F3x_S_M2_vrr = PAY*I_ERI_Py_S_F3x_S_M2_vrr+WPY*I_ERI_Py_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_F3x_S_M2_vrr = PAZ*I_ERI_Pz_S_F3x_S_M2_vrr+WPZ*I_ERI_Pz_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_F2xy_S_M2_vrr = PAX*I_ERI_Px_S_F2xy_S_M2_vrr+WPX*I_ERI_Px_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M3_vrr;
      Double I_ERI_Dxy_S_F2xy_S_M2_vrr = PAY*I_ERI_Px_S_F2xy_S_M2_vrr+WPY*I_ERI_Px_S_F2xy_S_M3_vrr+oned2k*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_F2xy_S_M2_vrr = PAY*I_ERI_Py_S_F2xy_S_M2_vrr+WPY*I_ERI_Py_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_F2xy_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_F2xz_S_M2_vrr = PAX*I_ERI_Px_S_F2xz_S_M2_vrr+WPX*I_ERI_Px_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M3_vrr;
      Double I_ERI_Dxy_S_F2xz_S_M2_vrr = PAY*I_ERI_Px_S_F2xz_S_M2_vrr+WPY*I_ERI_Px_S_F2xz_S_M3_vrr;
      Double I_ERI_D2y_S_F2xz_S_M2_vrr = PAY*I_ERI_Py_S_F2xz_S_M2_vrr+WPY*I_ERI_Py_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_D2z_S_F2xz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M2_vrr = PAX*I_ERI_Px_S_Fx2y_S_M2_vrr+WPX*I_ERI_Px_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_M2_vrr = PAY*I_ERI_Px_S_Fx2y_S_M2_vrr+WPY*I_ERI_Px_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M2_vrr = PAY*I_ERI_Py_S_Fx2y_S_M2_vrr+WPY*I_ERI_Py_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M2_vrr = PAX*I_ERI_Px_S_Fxyz_S_M2_vrr+WPX*I_ERI_Px_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Px_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M2_vrr = PAY*I_ERI_Py_S_Fxyz_S_M2_vrr+WPY*I_ERI_Py_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Py_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M2_vrr = PAX*I_ERI_Px_S_Fx2z_S_M2_vrr+WPX*I_ERI_Px_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_M2_vrr = PAY*I_ERI_Px_S_Fx2z_S_M2_vrr+WPY*I_ERI_Px_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M2_vrr = PAY*I_ERI_Py_S_Fx2z_S_M2_vrr+WPY*I_ERI_Py_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M3_vrr;
      Double I_ERI_D2x_S_F3y_S_M2_vrr = PAX*I_ERI_Px_S_F3y_S_M2_vrr+WPX*I_ERI_Px_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Dxy_S_F3y_S_M2_vrr = PAY*I_ERI_Px_S_F3y_S_M2_vrr+WPY*I_ERI_Px_S_F3y_S_M3_vrr+3*oned2k*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_F3y_S_M2_vrr = PAY*I_ERI_Py_S_F3y_S_M2_vrr+WPY*I_ERI_Py_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M3_vrr;
      Double I_ERI_D2z_S_F3y_S_M2_vrr = PAZ*I_ERI_Pz_S_F3y_S_M2_vrr+WPZ*I_ERI_Pz_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_F2yz_S_M2_vrr = PAX*I_ERI_Px_S_F2yz_S_M2_vrr+WPX*I_ERI_Px_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Dxy_S_F2yz_S_M2_vrr = PAY*I_ERI_Px_S_F2yz_S_M2_vrr+WPY*I_ERI_Px_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_F2yz_S_M2_vrr = PAY*I_ERI_Py_S_F2yz_S_M2_vrr+WPY*I_ERI_Py_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M3_vrr;
      Double I_ERI_D2z_S_F2yz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M2_vrr = PAX*I_ERI_Px_S_Fy2z_S_M2_vrr+WPX*I_ERI_Px_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M2_vrr = PAY*I_ERI_Py_S_Fy2z_S_M2_vrr+WPY*I_ERI_Py_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_Py_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M3_vrr;
      Double I_ERI_D2x_S_F3z_S_M2_vrr = PAX*I_ERI_Px_S_F3z_S_M2_vrr+WPX*I_ERI_Px_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Dxy_S_F3z_S_M2_vrr = PAY*I_ERI_Px_S_F3z_S_M2_vrr+WPY*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_F3z_S_M2_vrr = PAY*I_ERI_Py_S_F3z_S_M2_vrr+WPY*I_ERI_Py_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_F3z_S_M2_vrr = PAZ*I_ERI_Pz_S_F3z_S_M2_vrr+WPZ*I_ERI_Pz_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M3_vrr;

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
       * totally 8 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M1_vrr = PAX*I_ERI_Px_S_Px_S_M1_vrr+WPX*I_ERI_Px_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Px_S_M1_vrr = PAY*I_ERI_Py_S_Px_S_M1_vrr+WPY*I_ERI_Py_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Px_S_M1_vrr = PAZ*I_ERI_Pz_S_Px_S_M1_vrr+WPZ*I_ERI_Pz_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_Py_S_M1_vrr = PAX*I_ERI_Px_S_Py_S_M1_vrr+WPX*I_ERI_Px_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_D2y_S_Py_S_M1_vrr = PAY*I_ERI_Py_S_Py_S_M1_vrr+WPY*I_ERI_Py_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Py_S_M1_vrr = PAZ*I_ERI_Pz_S_Py_S_M1_vrr+WPZ*I_ERI_Pz_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_Pz_S_M1_vrr = PAX*I_ERI_Px_S_Pz_S_M1_vrr+WPX*I_ERI_Px_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_Pz_S_M1_vrr = PAY*I_ERI_Px_S_Pz_S_M1_vrr+WPY*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_Pz_S_M1_vrr = PAY*I_ERI_Py_S_Pz_S_M1_vrr+WPY*I_ERI_Py_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_Pz_S_M1_vrr = PAZ*I_ERI_Pz_S_Pz_S_M1_vrr+WPZ*I_ERI_Pz_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M1_vrr = QCX*I_ERI_S_S_D2x_S_M1_vrr+WQX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2e*I_ERI_S_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_F2xy_S_M1_vrr = QCY*I_ERI_S_S_D2x_S_M1_vrr+WQY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_F2xz_S_M1_vrr = QCZ*I_ERI_S_S_D2x_S_M1_vrr+WQZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_Fx2y_S_M1_vrr = QCX*I_ERI_S_S_D2y_S_M1_vrr+WQX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fxyz_S_M1_vrr = QCZ*I_ERI_S_S_Dxy_S_M1_vrr+WQZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_S_S_Fx2z_S_M1_vrr = QCX*I_ERI_S_S_D2z_S_M1_vrr+WQX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3y_S_M1_vrr = QCY*I_ERI_S_S_D2y_S_M1_vrr+WQY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2e*I_ERI_S_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_S_F2yz_S_M1_vrr = QCZ*I_ERI_S_S_D2y_S_M1_vrr+WQZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fy2z_S_M1_vrr = QCY*I_ERI_S_S_D2z_S_M1_vrr+WQY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3z_S_M1_vrr = QCZ*I_ERI_S_S_D2z_S_M1_vrr+WQZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2e*I_ERI_S_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
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
      Double I_ERI_D2y_S_Dxy_S_M1_vrr = PAY*I_ERI_Py_S_Dxy_S_M1_vrr+WPY*I_ERI_Py_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_Py_S_Px_S_M2_vrr;
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
       * shell quartet name: SQ_ERI_P_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M1_vrr = PAX*I_ERI_S_S_F3x_S_M1_vrr+WPX*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Py_S_F3x_S_M1_vrr = PAY*I_ERI_S_S_F3x_S_M1_vrr+WPY*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Pz_S_F3x_S_M1_vrr = PAZ*I_ERI_S_S_F3x_S_M1_vrr+WPZ*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Px_S_F2xy_S_M1_vrr = PAX*I_ERI_S_S_F2xy_S_M1_vrr+WPX*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Py_S_F2xy_S_M1_vrr = PAY*I_ERI_S_S_F2xy_S_M1_vrr+WPY*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Pz_S_F2xy_S_M1_vrr = PAZ*I_ERI_S_S_F2xy_S_M1_vrr+WPZ*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Px_S_F2xz_S_M1_vrr = PAX*I_ERI_S_S_F2xz_S_M1_vrr+WPX*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Py_S_F2xz_S_M1_vrr = PAY*I_ERI_S_S_F2xz_S_M1_vrr+WPY*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Pz_S_F2xz_S_M1_vrr = PAZ*I_ERI_S_S_F2xz_S_M1_vrr+WPZ*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Px_S_Fx2y_S_M1_vrr = PAX*I_ERI_S_S_Fx2y_S_M1_vrr+WPX*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Py_S_Fx2y_S_M1_vrr = PAY*I_ERI_S_S_Fx2y_S_M1_vrr+WPY*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_S_S_Fx2y_S_M1_vrr+WPZ*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Px_S_Fxyz_S_M1_vrr = PAX*I_ERI_S_S_Fxyz_S_M1_vrr+WPX*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Py_S_Fxyz_S_M1_vrr = PAY*I_ERI_S_S_Fxyz_S_M1_vrr+WPY*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_S_S_Fxyz_S_M1_vrr+WPZ*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Px_S_Fx2z_S_M1_vrr = PAX*I_ERI_S_S_Fx2z_S_M1_vrr+WPX*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Py_S_Fx2z_S_M1_vrr = PAY*I_ERI_S_S_Fx2z_S_M1_vrr+WPY*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_S_S_Fx2z_S_M1_vrr+WPZ*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Px_S_F3y_S_M1_vrr = PAX*I_ERI_S_S_F3y_S_M1_vrr+WPX*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Py_S_F3y_S_M1_vrr = PAY*I_ERI_S_S_F3y_S_M1_vrr+WPY*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Pz_S_F3y_S_M1_vrr = PAZ*I_ERI_S_S_F3y_S_M1_vrr+WPZ*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Px_S_F2yz_S_M1_vrr = PAX*I_ERI_S_S_F2yz_S_M1_vrr+WPX*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Py_S_F2yz_S_M1_vrr = PAY*I_ERI_S_S_F2yz_S_M1_vrr+WPY*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Pz_S_F2yz_S_M1_vrr = PAZ*I_ERI_S_S_F2yz_S_M1_vrr+WPZ*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Px_S_Fy2z_S_M1_vrr = PAX*I_ERI_S_S_Fy2z_S_M1_vrr+WPX*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Py_S_Fy2z_S_M1_vrr = PAY*I_ERI_S_S_Fy2z_S_M1_vrr+WPY*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_S_S_Fy2z_S_M1_vrr+WPZ*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Px_S_F3z_S_M1_vrr = PAX*I_ERI_S_S_F3z_S_M1_vrr+WPX*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Py_S_F3z_S_M1_vrr = PAY*I_ERI_S_S_F3z_S_M1_vrr+WPY*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Pz_S_F3z_S_M1_vrr = PAZ*I_ERI_S_S_F3z_S_M1_vrr+WPZ*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_S_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M1_vrr = PAX*I_ERI_Px_S_F3x_S_M1_vrr+WPX*I_ERI_Px_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxy_S_F3x_S_M1_vrr = PAY*I_ERI_Px_S_F3x_S_M1_vrr+WPY*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_Dxz_S_F3x_S_M1_vrr = PAZ*I_ERI_Px_S_F3x_S_M1_vrr+WPZ*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_D2y_S_F3x_S_M1_vrr = PAY*I_ERI_Py_S_F3x_S_M1_vrr+WPY*I_ERI_Py_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Dyz_S_F3x_S_M1_vrr = PAZ*I_ERI_Py_S_F3x_S_M1_vrr+WPZ*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_F3x_S_M1_vrr = PAZ*I_ERI_Pz_S_F3x_S_M1_vrr+WPZ*I_ERI_Pz_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_F2xy_S_M1_vrr = PAX*I_ERI_Px_S_F2xy_S_M1_vrr+WPX*I_ERI_Px_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xy_S_M1_vrr = PAY*I_ERI_Px_S_F2xy_S_M1_vrr+WPY*I_ERI_Px_S_F2xy_S_M2_vrr+oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxz_S_F2xy_S_M1_vrr = PAZ*I_ERI_Px_S_F2xy_S_M1_vrr+WPZ*I_ERI_Px_S_F2xy_S_M2_vrr;
      Double I_ERI_D2y_S_F2xy_S_M1_vrr = PAY*I_ERI_Py_S_F2xy_S_M1_vrr+WPY*I_ERI_Py_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_Dyz_S_F2xy_S_M1_vrr = PAZ*I_ERI_Py_S_F2xy_S_M1_vrr+WPZ*I_ERI_Py_S_F2xy_S_M2_vrr;
      Double I_ERI_D2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_F2xz_S_M1_vrr = PAX*I_ERI_Px_S_F2xz_S_M1_vrr+WPX*I_ERI_Px_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xz_S_M1_vrr = PAY*I_ERI_Px_S_F2xz_S_M1_vrr+WPY*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_Dxz_S_F2xz_S_M1_vrr = PAZ*I_ERI_Px_S_F2xz_S_M1_vrr+WPZ*I_ERI_Px_S_F2xz_S_M2_vrr+oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F2xz_S_M1_vrr = PAY*I_ERI_Py_S_F2xz_S_M1_vrr+WPY*I_ERI_Py_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Dyz_S_F2xz_S_M1_vrr = PAZ*I_ERI_Py_S_F2xz_S_M1_vrr+WPZ*I_ERI_Py_S_F2xz_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M1_vrr = PAX*I_ERI_Px_S_Fx2y_S_M1_vrr+WPX*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_M1_vrr = PAY*I_ERI_Px_S_Fx2y_S_M1_vrr+WPY*I_ERI_Px_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_Dxz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Px_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Px_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_Py_S_Fx2y_S_M1_vrr+WPY*I_ERI_Py_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M2_vrr;
      Double I_ERI_Dyz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Py_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Py_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M1_vrr = PAX*I_ERI_Px_S_Fxyz_S_M1_vrr+WPX*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_Dxy_S_Fxyz_S_M1_vrr = PAY*I_ERI_Px_S_Fxyz_S_M1_vrr+WPY*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_Py_S_Fxyz_S_M1_vrr+WPY*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M1_vrr = PAX*I_ERI_Px_S_Fx2z_S_M1_vrr+WPX*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_M1_vrr = PAY*I_ERI_Px_S_Fx2z_S_M1_vrr+WPY*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_Dxz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Px_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Px_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_Py_S_Fx2z_S_M1_vrr+WPY*I_ERI_Py_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Py_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Py_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M2_vrr;
      Double I_ERI_D2x_S_F3y_S_M1_vrr = PAX*I_ERI_Px_S_F3y_S_M1_vrr+WPX*I_ERI_Px_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Dxy_S_F3y_S_M1_vrr = PAY*I_ERI_Px_S_F3y_S_M1_vrr+WPY*I_ERI_Px_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_Dxz_S_F3y_S_M1_vrr = PAZ*I_ERI_Px_S_F3y_S_M1_vrr+WPZ*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_F3y_S_M1_vrr = PAY*I_ERI_Py_S_F3y_S_M1_vrr+WPY*I_ERI_Py_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_Dyz_S_F3y_S_M1_vrr = PAZ*I_ERI_Py_S_F3y_S_M1_vrr+WPZ*I_ERI_Py_S_F3y_S_M2_vrr;
      Double I_ERI_D2z_S_F3y_S_M1_vrr = PAZ*I_ERI_Pz_S_F3y_S_M1_vrr+WPZ*I_ERI_Pz_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_F2yz_S_M1_vrr = PAX*I_ERI_Px_S_F2yz_S_M1_vrr+WPX*I_ERI_Px_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Dxy_S_F2yz_S_M1_vrr = PAY*I_ERI_Px_S_F2yz_S_M1_vrr+WPY*I_ERI_Px_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_Dxz_S_F2yz_S_M1_vrr = PAZ*I_ERI_Px_S_F2yz_S_M1_vrr+WPZ*I_ERI_Px_S_F2yz_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_F2yz_S_M1_vrr = PAY*I_ERI_Py_S_F2yz_S_M1_vrr+WPY*I_ERI_Py_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M2_vrr;
      Double I_ERI_Dyz_S_F2yz_S_M1_vrr = PAZ*I_ERI_Py_S_F2yz_S_M1_vrr+WPZ*I_ERI_Py_S_F2yz_S_M2_vrr+oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M1_vrr = PAX*I_ERI_Px_S_Fy2z_S_M1_vrr+WPX*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Fy2z_S_M1_vrr = PAY*I_ERI_Px_S_Fy2z_S_M1_vrr+WPY*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_Py_S_Fy2z_S_M1_vrr+WPY*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M2_vrr;
      Double I_ERI_D2x_S_F3z_S_M1_vrr = PAX*I_ERI_Px_S_F3z_S_M1_vrr+WPX*I_ERI_Px_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Dxy_S_F3z_S_M1_vrr = PAY*I_ERI_Px_S_F3z_S_M1_vrr+WPY*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_Dxz_S_F3z_S_M1_vrr = PAZ*I_ERI_Px_S_F3z_S_M1_vrr+WPZ*I_ERI_Px_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_F3z_S_M1_vrr = PAY*I_ERI_Py_S_F3z_S_M1_vrr+WPY*I_ERI_Py_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Dyz_S_F3z_S_M1_vrr = PAZ*I_ERI_Py_S_F3z_S_M1_vrr+WPZ*I_ERI_Py_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_F3z_S_M1_vrr = PAZ*I_ERI_Pz_S_F3z_S_M1_vrr+WPZ*I_ERI_Pz_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M1_vrr = PAX*I_ERI_D2x_S_D2x_S_M1_vrr+WPX*I_ERI_D2x_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xy_S_D2x_S_M1_vrr = PAY*I_ERI_D2x_S_D2x_S_M1_vrr+WPY*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_D2x_S_M1_vrr = PAZ*I_ERI_D2x_S_D2x_S_M1_vrr+WPZ*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2x_S_M1_vrr = PAX*I_ERI_D2y_S_D2x_S_M1_vrr+WPX*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2x_S_M1_vrr = PAZ*I_ERI_Dxy_S_D2x_S_M1_vrr+WPZ*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2x_S_M1_vrr = PAX*I_ERI_D2z_S_D2x_S_M1_vrr+WPX*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3y_S_D2x_S_M1_vrr = PAY*I_ERI_D2y_S_D2x_S_M1_vrr+WPY*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_D2x_S_M1_vrr = PAZ*I_ERI_D2y_S_D2x_S_M1_vrr+WPZ*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2x_S_M1_vrr = PAY*I_ERI_D2z_S_D2x_S_M1_vrr+WPY*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_D2x_S_M1_vrr = PAZ*I_ERI_D2z_S_D2x_S_M1_vrr+WPZ*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Dxy_S_M1_vrr = PAX*I_ERI_D2x_S_Dxy_S_M1_vrr+WPX*I_ERI_D2x_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M1_vrr = PAY*I_ERI_D2x_S_Dxy_S_M1_vrr+WPY*I_ERI_D2x_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2x_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_M1_vrr = PAX*I_ERI_D2y_S_Dxy_S_M1_vrr+WPX*I_ERI_D2y_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_M1_vrr = PAX*I_ERI_D2z_S_Dxy_S_M1_vrr+WPX*I_ERI_D2z_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_Dxy_S_M1_vrr = PAY*I_ERI_D2y_S_Dxy_S_M1_vrr+WPY*I_ERI_D2y_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2y_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2z_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_F3x_S_Dxz_S_M1_vrr = PAX*I_ERI_D2x_S_Dxz_S_M1_vrr+WPX*I_ERI_D2x_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M1_vrr = PAY*I_ERI_D2x_S_Dxz_S_M1_vrr+WPY*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2x_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2x_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_M1_vrr = PAX*I_ERI_D2y_S_Dxz_S_M1_vrr+WPX*I_ERI_D2y_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_M1_vrr = PAX*I_ERI_D2z_S_Dxz_S_M1_vrr+WPX*I_ERI_D2z_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2z_S_Pz_S_M2_vrr;
      Double I_ERI_F3y_S_Dxz_S_M1_vrr = PAY*I_ERI_D2y_S_Dxz_S_M1_vrr+WPY*I_ERI_D2y_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2y_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2y_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2z_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_D2y_S_M1_vrr = PAX*I_ERI_D2x_S_D2y_S_M1_vrr+WPX*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_D2y_S_M1_vrr = PAY*I_ERI_D2x_S_D2y_S_M1_vrr+WPY*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_F2xz_S_D2y_S_M1_vrr = PAZ*I_ERI_D2x_S_D2y_S_M1_vrr+WPZ*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2y_S_M1_vrr = PAX*I_ERI_D2y_S_D2y_S_M1_vrr+WPX*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2y_S_M1_vrr = PAZ*I_ERI_Dxy_S_D2y_S_M1_vrr+WPZ*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2y_S_M1_vrr = PAX*I_ERI_D2z_S_D2y_S_M1_vrr+WPX*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_D2y_S_M1_vrr = PAY*I_ERI_D2y_S_D2y_S_M1_vrr+WPY*I_ERI_D2y_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F2yz_S_D2y_S_M1_vrr = PAZ*I_ERI_D2y_S_D2y_S_M1_vrr+WPZ*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2y_S_M1_vrr = PAY*I_ERI_D2z_S_D2y_S_M1_vrr+WPY*I_ERI_D2z_S_D2y_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_D2y_S_M1_vrr = PAZ*I_ERI_D2z_S_D2y_S_M1_vrr+WPZ*I_ERI_D2z_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Dyz_S_M1_vrr = PAX*I_ERI_D2x_S_Dyz_S_M1_vrr+WPX*I_ERI_D2x_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M1_vrr = PAY*I_ERI_D2x_S_Dyz_S_M1_vrr+WPY*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xz_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2x_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_M1_vrr = PAX*I_ERI_D2y_S_Dyz_S_M1_vrr+WPX*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_M1_vrr = PAX*I_ERI_D2z_S_Dyz_S_M1_vrr+WPX*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3y_S_Dyz_S_M1_vrr = PAY*I_ERI_D2y_S_Dyz_S_M1_vrr+WPY*I_ERI_D2y_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_F2yz_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2y_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2y_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2z_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3x_S_D2z_S_M1_vrr = PAX*I_ERI_D2x_S_D2z_S_M1_vrr+WPX*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_D2z_S_M1_vrr = PAY*I_ERI_D2x_S_D2z_S_M1_vrr+WPY*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_D2z_S_M1_vrr = PAZ*I_ERI_D2x_S_D2z_S_M1_vrr+WPZ*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2z_S_M1_vrr = PAX*I_ERI_D2y_S_D2z_S_M1_vrr+WPX*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fxyz_S_D2z_S_M1_vrr = PAZ*I_ERI_Dxy_S_D2z_S_M1_vrr+WPZ*I_ERI_Dxy_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2z_S_M1_vrr = PAX*I_ERI_D2z_S_D2z_S_M1_vrr+WPX*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_D2z_S_M1_vrr = PAY*I_ERI_D2y_S_D2z_S_M1_vrr+WPY*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_D2z_S_M1_vrr = PAZ*I_ERI_D2y_S_D2z_S_M1_vrr+WPZ*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fy2z_S_D2z_S_M1_vrr = PAY*I_ERI_D2z_S_D2z_S_M1_vrr+WPY*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3z_S_D2z_S_M1_vrr = PAZ*I_ERI_D2z_S_D2z_S_M1_vrr+WPZ*I_ERI_D2z_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M1_vrr = PAX*I_ERI_D2x_S_F3x_S_M1_vrr+WPX*I_ERI_D2x_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F3x_S_M1_vrr = PAY*I_ERI_D2x_S_F3x_S_M1_vrr+WPY*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2x_S_F3x_S_M1_vrr+WPZ*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3x_S_M1_vrr = PAX*I_ERI_D2y_S_F3x_S_M1_vrr+WPX*I_ERI_D2y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3x_S_M1_vrr = PAZ*I_ERI_Dxy_S_F3x_S_M1_vrr+WPZ*I_ERI_Dxy_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3x_S_M1_vrr = PAX*I_ERI_D2z_S_F3x_S_M1_vrr+WPX*I_ERI_D2z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F3x_S_M1_vrr = PAY*I_ERI_D2y_S_F3x_S_M1_vrr+WPY*I_ERI_D2y_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_F2yz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2y_S_F3x_S_M1_vrr+WPZ*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3x_S_M1_vrr = PAY*I_ERI_D2z_S_F3x_S_M1_vrr+WPY*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_F3z_S_F3x_S_M1_vrr = PAZ*I_ERI_D2z_S_F3x_S_M1_vrr+WPZ*I_ERI_D2z_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_F3x_S_F2xy_S_M1_vrr = PAX*I_ERI_D2x_S_F2xy_S_M1_vrr+WPX*I_ERI_D2x_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M1_vrr = PAY*I_ERI_D2x_S_F2xy_S_M1_vrr+WPY*I_ERI_D2x_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_M1_vrr = PAX*I_ERI_D2y_S_F2xy_S_M1_vrr+WPX*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2xy_S_M1_vrr = PAZ*I_ERI_Dxy_S_F2xy_S_M1_vrr+WPZ*I_ERI_Dxy_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_M1_vrr = PAX*I_ERI_D2z_S_F2xy_S_M1_vrr+WPX*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3y_S_F2xy_S_M1_vrr = PAY*I_ERI_D2y_S_F2xy_S_M1_vrr+WPY*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2xy_S_M1_vrr = PAY*I_ERI_D2z_S_F2xy_S_M1_vrr+WPY*I_ERI_D2z_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_F3x_S_F2xz_S_M1_vrr = PAX*I_ERI_D2x_S_F2xz_S_M1_vrr+WPX*I_ERI_D2x_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M1_vrr = PAY*I_ERI_D2x_S_F2xz_S_M1_vrr+WPY*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_M1_vrr = PAX*I_ERI_D2y_S_F2xz_S_M1_vrr+WPX*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2xz_S_M1_vrr = PAZ*I_ERI_Dxy_S_F2xz_S_M1_vrr+WPZ*I_ERI_Dxy_S_F2xz_S_M2_vrr+oned2k*I_ERI_Dxy_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_M1_vrr = PAX*I_ERI_D2z_S_F2xz_S_M1_vrr+WPX*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3y_S_F2xz_S_M1_vrr = PAY*I_ERI_D2y_S_F2xz_S_M1_vrr+WPY*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2xz_S_M1_vrr = PAY*I_ERI_D2z_S_F2xz_S_M1_vrr+WPY*I_ERI_D2z_S_F2xz_S_M2_vrr;
      Double I_ERI_F3z_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fxyz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Dxy_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Dxy_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2z_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3y_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fxyz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Dxy_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Dxy_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Dxy_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_Fy2z_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2z_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3x_S_F3y_S_M1_vrr = PAX*I_ERI_D2x_S_F3y_S_M1_vrr+WPX*I_ERI_D2x_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_F2xy_S_F3y_S_M1_vrr = PAY*I_ERI_D2x_S_F3y_S_M1_vrr+WPY*I_ERI_D2x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2x_S_F3y_S_M1_vrr+WPZ*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3y_S_M1_vrr = PAX*I_ERI_D2y_S_F3y_S_M1_vrr+WPX*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3y_S_M1_vrr = PAZ*I_ERI_Dxy_S_F3y_S_M1_vrr+WPZ*I_ERI_Dxy_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3y_S_M1_vrr = PAX*I_ERI_D2z_S_F3y_S_M1_vrr+WPX*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3y_S_F3y_S_M1_vrr = PAY*I_ERI_D2y_S_F3y_S_M1_vrr+WPY*I_ERI_D2y_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2y_S_F3y_S_M1_vrr+WPZ*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3y_S_M1_vrr = PAY*I_ERI_D2z_S_F3y_S_M1_vrr+WPY*I_ERI_D2z_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_F3y_S_M1_vrr = PAZ*I_ERI_D2z_S_F3y_S_M1_vrr+WPZ*I_ERI_D2z_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_F3x_S_F2yz_S_M1_vrr = PAX*I_ERI_D2x_S_F2yz_S_M1_vrr+WPX*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M1_vrr = PAY*I_ERI_D2x_S_F2yz_S_M1_vrr+WPY*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_M1_vrr = PAX*I_ERI_D2y_S_F2yz_S_M1_vrr+WPX*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_Fxyz_S_F2yz_S_M1_vrr = PAZ*I_ERI_Dxy_S_F2yz_S_M1_vrr+WPZ*I_ERI_Dxy_S_F2yz_S_M2_vrr+oned2k*I_ERI_Dxy_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_M1_vrr = PAX*I_ERI_D2z_S_F2yz_S_M1_vrr+WPX*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3y_S_F2yz_S_M1_vrr = PAY*I_ERI_D2y_S_F2yz_S_M1_vrr+WPY*I_ERI_D2y_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fy2z_S_F2yz_S_M1_vrr = PAY*I_ERI_D2z_S_F2yz_S_M1_vrr+WPY*I_ERI_D2z_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3z_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F3z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3x_S_F3z_S_M1_vrr = PAX*I_ERI_D2x_S_F3z_S_M1_vrr+WPX*I_ERI_D2x_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_F2xy_S_F3z_S_M1_vrr = PAY*I_ERI_D2x_S_F3z_S_M1_vrr+WPY*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2x_S_F3z_S_M1_vrr+WPZ*I_ERI_D2x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3z_S_M1_vrr = PAX*I_ERI_D2y_S_F3z_S_M1_vrr+WPX*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_Fxyz_S_F3z_S_M1_vrr = PAZ*I_ERI_Dxy_S_F3z_S_M1_vrr+WPZ*I_ERI_Dxy_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Dxy_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3z_S_M1_vrr = PAX*I_ERI_D2z_S_F3z_S_M1_vrr+WPX*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3y_S_F3z_S_M1_vrr = PAY*I_ERI_D2y_S_F3z_S_M1_vrr+WPY*I_ERI_D2y_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_F2yz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2y_S_F3z_S_M1_vrr+WPZ*I_ERI_D2y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fy2z_S_F3z_S_M1_vrr = PAY*I_ERI_D2z_S_F3z_S_M1_vrr+WPY*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3z_S_F3z_S_M1_vrr = PAZ*I_ERI_D2z_S_F3z_S_M1_vrr+WPZ*I_ERI_D2z_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dxy_S_vrr = QCY*I_ERI_S_S_Px_S_vrr+WQY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_vrr = PAX*I_ERI_S_S_D2x_S_vrr+WPX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_D2x_S_vrr = PAY*I_ERI_S_S_D2x_S_vrr+WPY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_D2x_S_vrr = PAZ*I_ERI_S_S_D2x_S_vrr+WPZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_D2y_S_vrr = PAX*I_ERI_S_S_D2y_S_vrr+WPX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_D2y_S_vrr = PAY*I_ERI_S_S_D2y_S_vrr+WPY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Pz_S_D2y_S_vrr = PAZ*I_ERI_S_S_D2y_S_vrr+WPZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_D2z_S_vrr = PAX*I_ERI_S_S_D2z_S_vrr+WPX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_D2z_S_vrr = PAY*I_ERI_S_S_D2z_S_vrr+WPY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_D2z_S_vrr = PAZ*I_ERI_S_S_D2z_S_vrr+WPZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 18 integrals are omitted 
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
      Double I_ERI_D2x_S_D2y_S_vrr = PAX*I_ERI_Px_S_D2y_S_vrr+WPX*I_ERI_Px_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_D2y_S_vrr = PAY*I_ERI_Px_S_D2y_S_vrr+WPY*I_ERI_Px_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_D2y_S_vrr = PAZ*I_ERI_Px_S_D2y_S_vrr+WPZ*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = PAY*I_ERI_Py_S_D2y_S_vrr+WPY*I_ERI_Py_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_D2y_S_vrr = PAZ*I_ERI_Py_S_D2y_S_vrr+WPZ*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = PAZ*I_ERI_Pz_S_D2y_S_vrr+WPZ*I_ERI_Pz_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = PAX*I_ERI_Px_S_D2z_S_vrr+WPX*I_ERI_Px_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_D2z_S_vrr = PAY*I_ERI_Px_S_D2z_S_vrr+WPY*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_D2z_S_vrr = PAZ*I_ERI_Px_S_D2z_S_vrr+WPZ*I_ERI_Px_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = PAY*I_ERI_Py_S_D2z_S_vrr+WPY*I_ERI_Py_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_D2z_S_vrr = PAZ*I_ERI_Py_S_D2z_S_vrr+WPZ*I_ERI_Py_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = PAZ*I_ERI_Pz_S_D2z_S_vrr+WPZ*I_ERI_Pz_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_vrr = PAX*I_ERI_S_S_F3x_S_vrr+WPX*I_ERI_S_S_F3x_S_M1_vrr+3*oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Py_S_F3x_S_vrr = PAY*I_ERI_S_S_F3x_S_vrr+WPY*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Pz_S_F3x_S_vrr = PAZ*I_ERI_S_S_F3x_S_vrr+WPZ*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Px_S_F2xy_S_vrr = PAX*I_ERI_S_S_F2xy_S_vrr+WPX*I_ERI_S_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Py_S_F2xy_S_vrr = PAY*I_ERI_S_S_F2xy_S_vrr+WPY*I_ERI_S_S_F2xy_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_F2xy_S_vrr = PAZ*I_ERI_S_S_F2xy_S_vrr+WPZ*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Px_S_F2xz_S_vrr = PAX*I_ERI_S_S_F2xz_S_vrr+WPX*I_ERI_S_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Py_S_F2xz_S_vrr = PAY*I_ERI_S_S_F2xz_S_vrr+WPY*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Pz_S_F2xz_S_vrr = PAZ*I_ERI_S_S_F2xz_S_vrr+WPZ*I_ERI_S_S_F2xz_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_Fx2y_S_vrr = PAX*I_ERI_S_S_Fx2y_S_vrr+WPX*I_ERI_S_S_Fx2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_Fx2y_S_vrr = PAY*I_ERI_S_S_Fx2y_S_vrr+WPY*I_ERI_S_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2y_S_vrr = PAZ*I_ERI_S_S_Fx2y_S_vrr+WPZ*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Px_S_Fxyz_S_vrr = PAX*I_ERI_S_S_Fxyz_S_vrr+WPX*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Py_S_Fxyz_S_vrr = PAY*I_ERI_S_S_Fxyz_S_vrr+WPY*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Pz_S_Fxyz_S_vrr = PAZ*I_ERI_S_S_Fxyz_S_vrr+WPZ*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Px_S_Fx2z_S_vrr = PAX*I_ERI_S_S_Fx2z_S_vrr+WPX*I_ERI_S_S_Fx2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_Fx2z_S_vrr = PAY*I_ERI_S_S_Fx2z_S_vrr+WPY*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2z_S_vrr = PAZ*I_ERI_S_S_Fx2z_S_vrr+WPZ*I_ERI_S_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Px_S_F3y_S_vrr = PAX*I_ERI_S_S_F3y_S_vrr+WPX*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Py_S_F3y_S_vrr = PAY*I_ERI_S_S_F3y_S_vrr+WPY*I_ERI_S_S_F3y_S_M1_vrr+3*oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Pz_S_F3y_S_vrr = PAZ*I_ERI_S_S_F3y_S_vrr+WPZ*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Px_S_F2yz_S_vrr = PAX*I_ERI_S_S_F2yz_S_vrr+WPX*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Py_S_F2yz_S_vrr = PAY*I_ERI_S_S_F2yz_S_vrr+WPY*I_ERI_S_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Pz_S_F2yz_S_vrr = PAZ*I_ERI_S_S_F2yz_S_vrr+WPZ*I_ERI_S_S_F2yz_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Fy2z_S_vrr = PAX*I_ERI_S_S_Fy2z_S_vrr+WPX*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Py_S_Fy2z_S_vrr = PAY*I_ERI_S_S_Fy2z_S_vrr+WPY*I_ERI_S_S_Fy2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fy2z_S_vrr = PAZ*I_ERI_S_S_Fy2z_S_vrr+WPZ*I_ERI_S_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Px_S_F3z_S_vrr = PAX*I_ERI_S_S_F3z_S_vrr+WPX*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Py_S_F3z_S_vrr = PAY*I_ERI_S_S_F3z_S_vrr+WPY*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Pz_S_F3z_S_vrr = PAZ*I_ERI_S_S_F3z_S_vrr+WPZ*I_ERI_S_S_F3z_S_M1_vrr+3*oned2k*I_ERI_S_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_vrr = PAX*I_ERI_Px_S_F3x_S_vrr+WPX*I_ERI_Px_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_F3x_S_vrr = PAY*I_ERI_Px_S_F3x_S_vrr+WPY*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_F3x_S_vrr = PAZ*I_ERI_Px_S_F3x_S_vrr+WPZ*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_F3x_S_vrr = PAY*I_ERI_Py_S_F3x_S_vrr+WPY*I_ERI_Py_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Dyz_S_F3x_S_vrr = PAZ*I_ERI_Py_S_F3x_S_vrr+WPZ*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_F3x_S_vrr = PAZ*I_ERI_Pz_S_F3x_S_vrr+WPZ*I_ERI_Pz_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_F2xy_S_vrr = PAX*I_ERI_Px_S_F2xy_S_vrr+WPX*I_ERI_Px_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xy_S_vrr = PAY*I_ERI_Px_S_F2xy_S_vrr+WPY*I_ERI_Px_S_F2xy_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xy_S_vrr = PAZ*I_ERI_Px_S_F2xy_S_vrr+WPZ*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_F2xy_S_vrr = PAY*I_ERI_Py_S_F2xy_S_vrr+WPY*I_ERI_Py_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xy_S_vrr = PAZ*I_ERI_Py_S_F2xy_S_vrr+WPZ*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_F2xy_S_vrr = PAZ*I_ERI_Pz_S_F2xy_S_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_F2xz_S_vrr = PAX*I_ERI_Px_S_F2xz_S_vrr+WPX*I_ERI_Px_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xz_S_vrr = PAY*I_ERI_Px_S_F2xz_S_vrr+WPY*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xz_S_vrr = PAZ*I_ERI_Px_S_F2xz_S_vrr+WPZ*I_ERI_Px_S_F2xz_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_F2xz_S_vrr = PAY*I_ERI_Py_S_F2xz_S_vrr+WPY*I_ERI_Py_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xz_S_vrr = PAZ*I_ERI_Py_S_F2xz_S_vrr+WPZ*I_ERI_Py_S_F2xz_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_F2xz_S_vrr = PAZ*I_ERI_Pz_S_F2xz_S_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr+oned2k*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2y_S_vrr = PAX*I_ERI_Px_S_Fx2y_S_vrr+WPX*I_ERI_Px_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_vrr = PAY*I_ERI_Px_S_Fx2y_S_vrr+WPY*I_ERI_Px_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2y_S_vrr = PAZ*I_ERI_Px_S_Fx2y_S_vrr+WPZ*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2y_S_vrr = PAY*I_ERI_Py_S_Fx2y_S_vrr+WPY*I_ERI_Py_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2y_S_vrr = PAZ*I_ERI_Py_S_Fx2y_S_vrr+WPZ*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2y_S_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fxyz_S_vrr = PAX*I_ERI_Px_S_Fxyz_S_vrr+WPX*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxy_S_Fxyz_S_vrr = PAY*I_ERI_Px_S_Fxyz_S_vrr+WPY*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxz_S_Fxyz_S_vrr = PAZ*I_ERI_Px_S_Fxyz_S_vrr+WPZ*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_D2y_S_Fxyz_S_vrr = PAY*I_ERI_Py_S_Fxyz_S_vrr+WPY*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_Dyz_S_Fxyz_S_vrr = PAZ*I_ERI_Py_S_Fxyz_S_vrr+WPZ*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_D2z_S_Fxyz_S_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2z_S_vrr = PAX*I_ERI_Px_S_Fx2z_S_vrr+WPX*I_ERI_Px_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_vrr = PAY*I_ERI_Px_S_Fx2z_S_vrr+WPY*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2z_S_vrr = PAZ*I_ERI_Px_S_Fx2z_S_vrr+WPZ*I_ERI_Px_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2z_S_vrr = PAY*I_ERI_Py_S_Fx2z_S_vrr+WPY*I_ERI_Py_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_vrr = PAZ*I_ERI_Py_S_Fx2z_S_vrr+WPZ*I_ERI_Py_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2z_S_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M1_vrr;
      Double I_ERI_D2x_S_F3y_S_vrr = PAX*I_ERI_Px_S_F3y_S_vrr+WPX*I_ERI_Px_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Dxy_S_F3y_S_vrr = PAY*I_ERI_Px_S_F3y_S_vrr+WPY*I_ERI_Px_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_F3y_S_vrr = PAZ*I_ERI_Px_S_F3y_S_vrr+WPZ*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_F3y_S_vrr = PAY*I_ERI_Py_S_F3y_S_vrr+WPY*I_ERI_Py_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_F3y_S_vrr = PAZ*I_ERI_Py_S_F3y_S_vrr+WPZ*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_F3y_S_vrr = PAZ*I_ERI_Pz_S_F3y_S_vrr+WPZ*I_ERI_Pz_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_F2yz_S_vrr = PAX*I_ERI_Px_S_F2yz_S_vrr+WPX*I_ERI_Px_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Dxy_S_F2yz_S_vrr = PAY*I_ERI_Px_S_F2yz_S_vrr+WPY*I_ERI_Px_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxz_S_F2yz_S_vrr = PAZ*I_ERI_Px_S_F2yz_S_vrr+WPZ*I_ERI_Px_S_F2yz_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_F2yz_S_vrr = PAY*I_ERI_Py_S_F2yz_S_vrr+WPY*I_ERI_Py_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M1_vrr;
      Double I_ERI_Dyz_S_F2yz_S_vrr = PAZ*I_ERI_Py_S_F2yz_S_vrr+WPZ*I_ERI_Py_S_F2yz_S_M1_vrr+oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_F2yz_S_vrr = PAZ*I_ERI_Pz_S_F2yz_S_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr+oned2k*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fy2z_S_vrr = PAX*I_ERI_Px_S_Fy2z_S_vrr+WPX*I_ERI_Px_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fy2z_S_vrr = PAY*I_ERI_Px_S_Fy2z_S_vrr+WPY*I_ERI_Px_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fy2z_S_vrr = PAZ*I_ERI_Px_S_Fy2z_S_vrr+WPZ*I_ERI_Px_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_D2y_S_Fy2z_S_vrr = PAY*I_ERI_Py_S_Fy2z_S_vrr+WPY*I_ERI_Py_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fy2z_S_vrr = PAZ*I_ERI_Py_S_Fy2z_S_vrr+WPZ*I_ERI_Py_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M1_vrr;
      Double I_ERI_D2z_S_Fy2z_S_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M1_vrr;
      Double I_ERI_D2x_S_F3z_S_vrr = PAX*I_ERI_Px_S_F3z_S_vrr+WPX*I_ERI_Px_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Dxy_S_F3z_S_vrr = PAY*I_ERI_Px_S_F3z_S_vrr+WPY*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxz_S_F3z_S_vrr = PAZ*I_ERI_Px_S_F3z_S_vrr+WPZ*I_ERI_Px_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_D2y_S_F3z_S_vrr = PAY*I_ERI_Py_S_F3z_S_vrr+WPY*I_ERI_Py_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Dyz_S_F3z_S_vrr = PAZ*I_ERI_Py_S_F3z_S_vrr+WPZ*I_ERI_Py_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_D2z_S_F3z_S_vrr = PAZ*I_ERI_Pz_S_F3z_S_vrr+WPZ*I_ERI_Pz_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
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
       * shell quartet name: SQ_ERI_F_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_vrr = PAX*I_ERI_D2x_S_F3x_S_vrr+WPX*I_ERI_D2x_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_F3x_S_vrr = PAY*I_ERI_D2x_S_F3x_S_vrr+WPY*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_F3x_S_vrr = PAZ*I_ERI_D2x_S_F3x_S_vrr+WPZ*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3x_S_vrr = PAX*I_ERI_D2y_S_F3x_S_vrr+WPX*I_ERI_D2y_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3x_S_vrr = PAZ*I_ERI_Dxy_S_F3x_S_vrr+WPZ*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3x_S_vrr = PAX*I_ERI_D2z_S_F3x_S_vrr+WPX*I_ERI_D2z_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_F3x_S_vrr = PAY*I_ERI_D2y_S_F3x_S_vrr+WPY*I_ERI_D2y_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_F3x_S_vrr = PAZ*I_ERI_D2y_S_F3x_S_vrr+WPZ*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3x_S_vrr = PAY*I_ERI_D2z_S_F3x_S_vrr+WPY*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_F3x_S_vrr = PAZ*I_ERI_D2z_S_F3x_S_vrr+WPZ*I_ERI_D2z_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_F2xy_S_vrr = PAX*I_ERI_D2x_S_F2xy_S_vrr+WPX*I_ERI_D2x_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xy_S_vrr = PAY*I_ERI_D2x_S_F2xy_S_vrr+WPY*I_ERI_D2x_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xy_S_vrr = PAZ*I_ERI_D2x_S_F2xy_S_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_vrr = PAX*I_ERI_D2y_S_F2xy_S_vrr+WPX*I_ERI_D2y_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xy_S_vrr = PAZ*I_ERI_Dxy_S_F2xy_S_vrr+WPZ*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_vrr = PAX*I_ERI_D2z_S_F2xy_S_vrr+WPX*I_ERI_D2z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3y_S_F2xy_S_vrr = PAY*I_ERI_D2y_S_F2xy_S_vrr+WPY*I_ERI_D2y_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xy_S_vrr = PAZ*I_ERI_D2y_S_F2xy_S_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xy_S_vrr = PAY*I_ERI_D2z_S_F2xy_S_vrr+WPY*I_ERI_D2z_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_F2xy_S_vrr = PAZ*I_ERI_D2z_S_F2xy_S_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_F2xz_S_vrr = PAX*I_ERI_D2x_S_F2xz_S_vrr+WPX*I_ERI_D2x_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xz_S_vrr = PAY*I_ERI_D2x_S_F2xz_S_vrr+WPY*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xz_S_vrr = PAZ*I_ERI_D2x_S_F2xz_S_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_vrr = PAX*I_ERI_D2y_S_F2xz_S_vrr+WPX*I_ERI_D2y_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xz_S_vrr = PAZ*I_ERI_Dxy_S_F2xz_S_vrr+WPZ*I_ERI_Dxy_S_F2xz_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_vrr = PAX*I_ERI_D2z_S_F2xz_S_vrr+WPX*I_ERI_D2z_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3y_S_F2xz_S_vrr = PAY*I_ERI_D2y_S_F2xz_S_vrr+WPY*I_ERI_D2y_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xz_S_vrr = PAZ*I_ERI_D2y_S_F2xz_S_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xz_S_vrr = PAY*I_ERI_D2z_S_F2xz_S_vrr+WPY*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_F2xz_S_vrr = PAZ*I_ERI_D2z_S_F2xz_S_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2y_S_vrr = PAX*I_ERI_D2x_S_Fx2y_S_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_vrr = PAY*I_ERI_D2x_S_Fx2y_S_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_vrr = PAX*I_ERI_D2y_S_Fx2y_S_vrr+WPX*I_ERI_D2y_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2y_S_vrr = PAZ*I_ERI_Dxy_S_Fx2y_S_vrr+WPZ*I_ERI_Dxy_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_vrr = PAX*I_ERI_D2z_S_Fx2y_S_vrr+WPX*I_ERI_D2z_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2y_S_vrr = PAY*I_ERI_D2y_S_Fx2y_S_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2y_S_vrr = PAY*I_ERI_D2z_S_Fx2y_S_vrr+WPY*I_ERI_D2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2y_S_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fxyz_S_vrr = PAX*I_ERI_D2x_S_Fxyz_S_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_vrr = PAY*I_ERI_D2x_S_Fxyz_S_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_vrr = PAX*I_ERI_D2y_S_Fxyz_S_vrr+WPX*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fxyz_S_vrr = PAZ*I_ERI_Dxy_S_Fxyz_S_vrr+WPZ*I_ERI_Dxy_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_vrr = PAX*I_ERI_D2z_S_Fxyz_S_vrr+WPX*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3y_S_Fxyz_S_vrr = PAY*I_ERI_D2y_S_Fxyz_S_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fxyz_S_vrr = PAY*I_ERI_D2z_S_Fxyz_S_vrr+WPY*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3z_S_Fxyz_S_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2z_S_vrr = PAX*I_ERI_D2x_S_Fx2z_S_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_vrr = PAY*I_ERI_D2x_S_Fx2z_S_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_vrr = PAX*I_ERI_D2y_S_Fx2z_S_vrr+WPX*I_ERI_D2y_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2z_S_vrr = PAZ*I_ERI_Dxy_S_Fx2z_S_vrr+WPZ*I_ERI_Dxy_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Dxz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_vrr = PAX*I_ERI_D2z_S_Fx2z_S_vrr+WPX*I_ERI_D2z_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2z_S_vrr = PAY*I_ERI_D2y_S_Fx2z_S_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2z_S_vrr = PAY*I_ERI_D2z_S_Fx2z_S_vrr+WPY*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2z_S_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3x_S_F3y_S_vrr = PAX*I_ERI_D2x_S_F3y_S_vrr+WPX*I_ERI_D2x_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_F3y_S_vrr = PAY*I_ERI_D2x_S_F3y_S_vrr+WPY*I_ERI_D2x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_F3y_S_vrr = PAZ*I_ERI_D2x_S_F3y_S_vrr+WPZ*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3y_S_vrr = PAX*I_ERI_D2y_S_F3y_S_vrr+WPX*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3y_S_vrr = PAZ*I_ERI_Dxy_S_F3y_S_vrr+WPZ*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3y_S_vrr = PAX*I_ERI_D2z_S_F3y_S_vrr+WPX*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_F3y_S_vrr = PAY*I_ERI_D2y_S_F3y_S_vrr+WPY*I_ERI_D2y_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_F3y_S_vrr = PAZ*I_ERI_D2y_S_F3y_S_vrr+WPZ*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3y_S_vrr = PAY*I_ERI_D2z_S_F3y_S_vrr+WPY*I_ERI_D2z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3z_S_F3y_S_vrr = PAZ*I_ERI_D2z_S_F3y_S_vrr+WPZ*I_ERI_D2z_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_F2yz_S_vrr = PAX*I_ERI_D2x_S_F2yz_S_vrr+WPX*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xy_S_F2yz_S_vrr = PAY*I_ERI_D2x_S_F2yz_S_vrr+WPY*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xz_S_F2yz_S_vrr = PAZ*I_ERI_D2x_S_F2yz_S_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_vrr = PAX*I_ERI_D2y_S_F2yz_S_vrr+WPX*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2yz_S_vrr = PAZ*I_ERI_Dxy_S_F2yz_S_vrr+WPZ*I_ERI_Dxy_S_F2yz_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_vrr = PAX*I_ERI_D2z_S_F2yz_S_vrr+WPX*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_F2yz_S_vrr = PAY*I_ERI_D2y_S_F2yz_S_vrr+WPY*I_ERI_D2y_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_F2yz_S_F2yz_S_vrr = PAZ*I_ERI_D2y_S_F2yz_S_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2yz_S_vrr = PAY*I_ERI_D2z_S_F2yz_S_vrr+WPY*I_ERI_D2z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3z_S_F2yz_S_vrr = PAZ*I_ERI_D2z_S_F2yz_S_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fy2z_S_vrr = PAX*I_ERI_D2x_S_Fy2z_S_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_vrr = PAY*I_ERI_D2x_S_Fy2z_S_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_vrr = PAX*I_ERI_D2y_S_Fy2z_S_vrr+WPX*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fy2z_S_vrr = PAZ*I_ERI_Dxy_S_Fy2z_S_vrr+WPZ*I_ERI_Dxy_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Dyz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_vrr = PAX*I_ERI_D2z_S_Fy2z_S_vrr+WPX*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fy2z_S_vrr = PAY*I_ERI_D2y_S_Fy2z_S_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fy2z_S_vrr = PAY*I_ERI_D2z_S_Fy2z_S_vrr+WPY*I_ERI_D2z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fy2z_S_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3x_S_F3z_S_vrr = PAX*I_ERI_D2x_S_F3z_S_vrr+WPX*I_ERI_D2x_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_F3z_S_vrr = PAY*I_ERI_D2x_S_F3z_S_vrr+WPY*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_F3z_S_vrr = PAZ*I_ERI_D2x_S_F3z_S_vrr+WPZ*I_ERI_D2x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3z_S_vrr = PAX*I_ERI_D2y_S_F3z_S_vrr+WPX*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3z_S_vrr = PAZ*I_ERI_Dxy_S_F3z_S_vrr+WPZ*I_ERI_Dxy_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3z_S_vrr = PAX*I_ERI_D2z_S_F3z_S_vrr+WPX*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_F3z_S_vrr = PAY*I_ERI_D2y_S_F3z_S_vrr+WPY*I_ERI_D2y_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_F3z_S_vrr = PAZ*I_ERI_D2y_S_F3z_S_vrr+WPZ*I_ERI_D2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3z_S_vrr = PAY*I_ERI_D2z_S_F3z_S_vrr+WPY*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_F3z_S_vrr = PAZ*I_ERI_D2z_S_F3z_S_vrr+WPZ*I_ERI_D2z_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_vrr = QCX*I_ERI_D2x_S_F3x_S_vrr+WQX*I_ERI_D2x_S_F3x_S_M1_vrr+3*oned2e*I_ERI_D2x_S_D2x_S_vrr-3*rhod2esq*I_ERI_D2x_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxy_S_G4x_S_vrr = QCX*I_ERI_Dxy_S_F3x_S_vrr+WQX*I_ERI_Dxy_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Dxy_S_D2x_S_vrr-3*rhod2esq*I_ERI_Dxy_S_D2x_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_G4x_S_vrr = QCX*I_ERI_Dxz_S_F3x_S_vrr+WQX*I_ERI_Dxz_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Dxz_S_D2x_S_vrr-3*rhod2esq*I_ERI_Dxz_S_D2x_S_M1_vrr+oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G4x_S_vrr = QCX*I_ERI_D2y_S_F3x_S_vrr+WQX*I_ERI_D2y_S_F3x_S_M1_vrr+3*oned2e*I_ERI_D2y_S_D2x_S_vrr-3*rhod2esq*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_G4x_S_vrr = QCX*I_ERI_Dyz_S_F3x_S_vrr+WQX*I_ERI_Dyz_S_F3x_S_M1_vrr+3*oned2e*I_ERI_Dyz_S_D2x_S_vrr-3*rhod2esq*I_ERI_Dyz_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_G4x_S_vrr = QCX*I_ERI_D2z_S_F3x_S_vrr+WQX*I_ERI_D2z_S_F3x_S_M1_vrr+3*oned2e*I_ERI_D2z_S_D2x_S_vrr-3*rhod2esq*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_G3xy_S_vrr = QCY*I_ERI_D2x_S_F3x_S_vrr+WQY*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Dxy_S_G3xy_S_vrr = QCY*I_ERI_Dxy_S_F3x_S_vrr+WQY*I_ERI_Dxy_S_F3x_S_M1_vrr+oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_G3xy_S_vrr = QCY*I_ERI_Dxz_S_F3x_S_vrr+WQY*I_ERI_Dxz_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G3xy_S_vrr = QCY*I_ERI_D2y_S_F3x_S_vrr+WQY*I_ERI_D2y_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_Dyz_S_G3xy_S_vrr = QCY*I_ERI_Dyz_S_F3x_S_vrr+WQY*I_ERI_Dyz_S_F3x_S_M1_vrr+oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_G3xy_S_vrr = QCY*I_ERI_D2z_S_F3x_S_vrr+WQY*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_G3xz_S_vrr = QCZ*I_ERI_D2x_S_F3x_S_vrr+WQZ*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Dxy_S_G3xz_S_vrr = QCZ*I_ERI_Dxy_S_F3x_S_vrr+WQZ*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_G3xz_S_vrr = QCZ*I_ERI_Dxz_S_F3x_S_vrr+WQZ*I_ERI_Dxz_S_F3x_S_M1_vrr+oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G3xz_S_vrr = QCZ*I_ERI_D2y_S_F3x_S_vrr+WQZ*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Dyz_S_G3xz_S_vrr = QCZ*I_ERI_Dyz_S_F3x_S_vrr+WQZ*I_ERI_Dyz_S_F3x_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_G3xz_S_vrr = QCZ*I_ERI_D2z_S_F3x_S_vrr+WQZ*I_ERI_D2z_S_F3x_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2y_S_vrr = QCY*I_ERI_D2x_S_F2xy_S_vrr+WQY*I_ERI_D2x_S_F2xy_S_M1_vrr+oned2e*I_ERI_D2x_S_D2x_S_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_G2x2y_S_vrr = QCY*I_ERI_Dxy_S_F2xy_S_vrr+WQY*I_ERI_Dxy_S_F2xy_S_M1_vrr+oned2e*I_ERI_Dxy_S_D2x_S_vrr-rhod2esq*I_ERI_Dxy_S_D2x_S_M1_vrr+oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_Dxz_S_G2x2y_S_vrr = QCY*I_ERI_Dxz_S_F2xy_S_vrr+WQY*I_ERI_Dxz_S_F2xy_S_M1_vrr+oned2e*I_ERI_Dxz_S_D2x_S_vrr-rhod2esq*I_ERI_Dxz_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2y_S_vrr = QCY*I_ERI_D2y_S_F2xy_S_vrr+WQY*I_ERI_D2y_S_F2xy_S_M1_vrr+oned2e*I_ERI_D2y_S_D2x_S_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_Dyz_S_G2x2y_S_vrr = QCY*I_ERI_Dyz_S_F2xy_S_vrr+WQY*I_ERI_Dyz_S_F2xy_S_M1_vrr+oned2e*I_ERI_Dyz_S_D2x_S_vrr-rhod2esq*I_ERI_Dyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2y_S_vrr = QCY*I_ERI_D2z_S_F2xy_S_vrr+WQY*I_ERI_D2z_S_F2xy_S_M1_vrr+oned2e*I_ERI_D2z_S_D2x_S_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_G2xyz_S_vrr = QCZ*I_ERI_D2x_S_F2xy_S_vrr+WQZ*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Dxy_S_G2xyz_S_vrr = QCZ*I_ERI_Dxy_S_F2xy_S_vrr+WQZ*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Dxz_S_G2xyz_S_vrr = QCZ*I_ERI_Dxz_S_F2xy_S_vrr+WQZ*I_ERI_Dxz_S_F2xy_S_M1_vrr+oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_G2xyz_S_vrr = QCZ*I_ERI_D2y_S_F2xy_S_vrr+WQZ*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Dyz_S_G2xyz_S_vrr = QCZ*I_ERI_Dyz_S_F2xy_S_vrr+WQZ*I_ERI_Dyz_S_F2xy_S_M1_vrr+oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_G2xyz_S_vrr = QCZ*I_ERI_D2z_S_F2xy_S_vrr+WQZ*I_ERI_D2z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2z_S_vrr = QCZ*I_ERI_D2x_S_F2xz_S_vrr+WQZ*I_ERI_D2x_S_F2xz_S_M1_vrr+oned2e*I_ERI_D2x_S_D2x_S_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_G2x2z_S_vrr = QCZ*I_ERI_Dxy_S_F2xz_S_vrr+WQZ*I_ERI_Dxy_S_F2xz_S_M1_vrr+oned2e*I_ERI_Dxy_S_D2x_S_vrr-rhod2esq*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_G2x2z_S_vrr = QCZ*I_ERI_Dxz_S_F2xz_S_vrr+WQZ*I_ERI_Dxz_S_F2xz_S_M1_vrr+oned2e*I_ERI_Dxz_S_D2x_S_vrr-rhod2esq*I_ERI_Dxz_S_D2x_S_M1_vrr+oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2z_S_vrr = QCZ*I_ERI_D2y_S_F2xz_S_vrr+WQZ*I_ERI_D2y_S_F2xz_S_M1_vrr+oned2e*I_ERI_D2y_S_D2x_S_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_G2x2z_S_vrr = QCZ*I_ERI_Dyz_S_F2xz_S_vrr+WQZ*I_ERI_Dyz_S_F2xz_S_M1_vrr+oned2e*I_ERI_Dyz_S_D2x_S_vrr-rhod2esq*I_ERI_Dyz_S_D2x_S_M1_vrr+oned2k*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2z_S_vrr = QCZ*I_ERI_D2z_S_F2xz_S_vrr+WQZ*I_ERI_D2z_S_F2xz_S_M1_vrr+oned2e*I_ERI_D2z_S_D2x_S_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3y_S_vrr = QCX*I_ERI_D2x_S_F3y_S_vrr+WQX*I_ERI_D2x_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx3y_S_vrr = QCX*I_ERI_Dxy_S_F3y_S_vrr+WQX*I_ERI_Dxy_S_F3y_S_M1_vrr+oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx3y_S_vrr = QCX*I_ERI_Dxz_S_F3y_S_vrr+WQX*I_ERI_Dxz_S_F3y_S_M1_vrr+oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3y_S_vrr = QCX*I_ERI_D2y_S_F3y_S_vrr+WQX*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx3y_S_vrr = QCX*I_ERI_Dyz_S_F3y_S_vrr+WQX*I_ERI_Dyz_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3y_S_vrr = QCX*I_ERI_D2z_S_F3y_S_vrr+WQX*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_vrr = QCZ*I_ERI_D2x_S_Fx2y_S_vrr+WQZ*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx2yz_S_vrr = QCZ*I_ERI_Dxy_S_Fx2y_S_vrr+WQZ*I_ERI_Dxy_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx2yz_S_vrr = QCZ*I_ERI_Dxz_S_Fx2y_S_vrr+WQZ*I_ERI_Dxz_S_Fx2y_S_M1_vrr+oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_vrr = QCZ*I_ERI_D2y_S_Fx2y_S_vrr+WQZ*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx2yz_S_vrr = QCZ*I_ERI_Dyz_S_Fx2y_S_vrr+WQZ*I_ERI_Dyz_S_Fx2y_S_M1_vrr+oned2k*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_vrr = QCZ*I_ERI_D2z_S_Fx2y_S_vrr+WQZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_vrr = QCY*I_ERI_D2x_S_Fx2z_S_vrr+WQY*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gxy2z_S_vrr = QCY*I_ERI_Dxy_S_Fx2z_S_vrr+WQY*I_ERI_Dxy_S_Fx2z_S_M1_vrr+oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gxy2z_S_vrr = QCY*I_ERI_Dxz_S_Fx2z_S_vrr+WQY*I_ERI_Dxz_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_vrr = QCY*I_ERI_D2y_S_Fx2z_S_vrr+WQY*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gxy2z_S_vrr = QCY*I_ERI_Dyz_S_Fx2z_S_vrr+WQY*I_ERI_Dyz_S_Fx2z_S_M1_vrr+oned2k*I_ERI_Pz_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_vrr = QCY*I_ERI_D2z_S_Fx2z_S_vrr+WQY*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3z_S_vrr = QCX*I_ERI_D2x_S_F3z_S_vrr+WQX*I_ERI_D2x_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx3z_S_vrr = QCX*I_ERI_Dxy_S_F3z_S_vrr+WQX*I_ERI_Dxy_S_F3z_S_M1_vrr+oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx3z_S_vrr = QCX*I_ERI_Dxz_S_F3z_S_vrr+WQX*I_ERI_Dxz_S_F3z_S_M1_vrr+oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3z_S_vrr = QCX*I_ERI_D2y_S_F3z_S_vrr+WQX*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx3z_S_vrr = QCX*I_ERI_Dyz_S_F3z_S_vrr+WQX*I_ERI_Dyz_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3z_S_vrr = QCX*I_ERI_D2z_S_F3z_S_vrr+WQX*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_D2x_S_G4y_S_vrr = QCY*I_ERI_D2x_S_F3y_S_vrr+WQY*I_ERI_D2x_S_F3y_S_M1_vrr+3*oned2e*I_ERI_D2x_S_D2y_S_vrr-3*rhod2esq*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_G4y_S_vrr = QCY*I_ERI_Dxy_S_F3y_S_vrr+WQY*I_ERI_Dxy_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Dxy_S_D2y_S_vrr-3*rhod2esq*I_ERI_Dxy_S_D2y_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_Dxz_S_G4y_S_vrr = QCY*I_ERI_Dxz_S_F3y_S_vrr+WQY*I_ERI_Dxz_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Dxz_S_D2y_S_vrr-3*rhod2esq*I_ERI_Dxz_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_G4y_S_vrr = QCY*I_ERI_D2y_S_F3y_S_vrr+WQY*I_ERI_D2y_S_F3y_S_M1_vrr+3*oned2e*I_ERI_D2y_S_D2y_S_vrr-3*rhod2esq*I_ERI_D2y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_Dyz_S_G4y_S_vrr = QCY*I_ERI_Dyz_S_F3y_S_vrr+WQY*I_ERI_Dyz_S_F3y_S_M1_vrr+3*oned2e*I_ERI_Dyz_S_D2y_S_vrr-3*rhod2esq*I_ERI_Dyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_G4y_S_vrr = QCY*I_ERI_D2z_S_F3y_S_vrr+WQY*I_ERI_D2z_S_F3y_S_M1_vrr+3*oned2e*I_ERI_D2z_S_D2y_S_vrr-3*rhod2esq*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_G3yz_S_vrr = QCZ*I_ERI_D2x_S_F3y_S_vrr+WQZ*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Dxy_S_G3yz_S_vrr = QCZ*I_ERI_Dxy_S_F3y_S_vrr+WQZ*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Dxz_S_G3yz_S_vrr = QCZ*I_ERI_Dxz_S_F3y_S_vrr+WQZ*I_ERI_Dxz_S_F3y_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_G3yz_S_vrr = QCZ*I_ERI_D2y_S_F3y_S_vrr+WQZ*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Dyz_S_G3yz_S_vrr = QCZ*I_ERI_Dyz_S_F3y_S_vrr+WQZ*I_ERI_Dyz_S_F3y_S_M1_vrr+oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_G3yz_S_vrr = QCZ*I_ERI_D2z_S_F3y_S_vrr+WQZ*I_ERI_D2z_S_F3y_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_G2y2z_S_vrr = QCZ*I_ERI_D2x_S_F2yz_S_vrr+WQZ*I_ERI_D2x_S_F2yz_S_M1_vrr+oned2e*I_ERI_D2x_S_D2y_S_vrr-rhod2esq*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_G2y2z_S_vrr = QCZ*I_ERI_Dxy_S_F2yz_S_vrr+WQZ*I_ERI_Dxy_S_F2yz_S_M1_vrr+oned2e*I_ERI_Dxy_S_D2y_S_vrr-rhod2esq*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_G2y2z_S_vrr = QCZ*I_ERI_Dxz_S_F2yz_S_vrr+WQZ*I_ERI_Dxz_S_F2yz_S_M1_vrr+oned2e*I_ERI_Dxz_S_D2y_S_vrr-rhod2esq*I_ERI_Dxz_S_D2y_S_M1_vrr+oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_D2y_S_G2y2z_S_vrr = QCZ*I_ERI_D2y_S_F2yz_S_vrr+WQZ*I_ERI_D2y_S_F2yz_S_M1_vrr+oned2e*I_ERI_D2y_S_D2y_S_vrr-rhod2esq*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_G2y2z_S_vrr = QCZ*I_ERI_Dyz_S_F2yz_S_vrr+WQZ*I_ERI_Dyz_S_F2yz_S_M1_vrr+oned2e*I_ERI_Dyz_S_D2y_S_vrr-rhod2esq*I_ERI_Dyz_S_D2y_S_M1_vrr+oned2k*I_ERI_Py_S_F2yz_S_M1_vrr;
      Double I_ERI_D2z_S_G2y2z_S_vrr = QCZ*I_ERI_D2z_S_F2yz_S_vrr+WQZ*I_ERI_D2z_S_F2yz_S_M1_vrr+oned2e*I_ERI_D2z_S_D2y_S_vrr-rhod2esq*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M1_vrr;
      Double I_ERI_D2x_S_Gy3z_S_vrr = QCY*I_ERI_D2x_S_F3z_S_vrr+WQY*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gy3z_S_vrr = QCY*I_ERI_Dxy_S_F3z_S_vrr+WQY*I_ERI_Dxy_S_F3z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gy3z_S_vrr = QCY*I_ERI_Dxz_S_F3z_S_vrr+WQY*I_ERI_Dxz_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gy3z_S_vrr = QCY*I_ERI_D2y_S_F3z_S_vrr+WQY*I_ERI_D2y_S_F3z_S_M1_vrr+2*oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gy3z_S_vrr = QCY*I_ERI_Dyz_S_F3z_S_vrr+WQY*I_ERI_Dyz_S_F3z_S_M1_vrr+oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gy3z_S_vrr = QCY*I_ERI_D2z_S_F3z_S_vrr+WQY*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_D2x_S_G4z_S_vrr = QCZ*I_ERI_D2x_S_F3z_S_vrr+WQZ*I_ERI_D2x_S_F3z_S_M1_vrr+3*oned2e*I_ERI_D2x_S_D2z_S_vrr-3*rhod2esq*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_G4z_S_vrr = QCZ*I_ERI_Dxy_S_F3z_S_vrr+WQZ*I_ERI_Dxy_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Dxy_S_D2z_S_vrr-3*rhod2esq*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_G4z_S_vrr = QCZ*I_ERI_Dxz_S_F3z_S_vrr+WQZ*I_ERI_Dxz_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Dxz_S_D2z_S_vrr-3*rhod2esq*I_ERI_Dxz_S_D2z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_G4z_S_vrr = QCZ*I_ERI_D2y_S_F3z_S_vrr+WQZ*I_ERI_D2y_S_F3z_S_M1_vrr+3*oned2e*I_ERI_D2y_S_D2z_S_vrr-3*rhod2esq*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_G4z_S_vrr = QCZ*I_ERI_Dyz_S_F3z_S_vrr+WQZ*I_ERI_Dyz_S_F3z_S_M1_vrr+3*oned2e*I_ERI_Dyz_S_D2z_S_vrr-3*rhod2esq*I_ERI_Dyz_S_D2z_S_M1_vrr+oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_G4z_S_vrr = QCZ*I_ERI_D2z_S_F3z_S_vrr+WQZ*I_ERI_D2z_S_F3z_S_M1_vrr+3*oned2e*I_ERI_D2z_S_D2z_S_vrr-3*rhod2esq*I_ERI_D2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;

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
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_vrr = PAX*I_ERI_F3x_S_F3x_S_vrr+WPX*I_ERI_F3x_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F3x_S_vrr = PAY*I_ERI_F3x_S_F3x_S_vrr+WPY*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xz_S_F3x_S_vrr = PAZ*I_ERI_F3x_S_F3x_S_vrr+WPZ*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3x_S_vrr = PAY*I_ERI_F2xy_S_F3x_S_vrr+WPY*I_ERI_F2xy_S_F3x_S_M1_vrr+oned2z*I_ERI_D2x_S_F3x_S_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3x_S_vrr = PAZ*I_ERI_F2xy_S_F3x_S_vrr+WPZ*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3x_S_vrr = PAZ*I_ERI_F2xz_S_F3x_S_vrr+WPZ*I_ERI_F2xz_S_F3x_S_M1_vrr+oned2z*I_ERI_D2x_S_F3x_S_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3x_S_vrr = PAX*I_ERI_F3y_S_F3x_S_vrr+WPX*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3x_S_vrr = PAZ*I_ERI_Fx2y_S_F3x_S_vrr+WPZ*I_ERI_Fx2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3x_S_vrr = PAY*I_ERI_Fx2z_S_F3x_S_vrr+WPY*I_ERI_Fx2z_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3x_S_vrr = PAX*I_ERI_F3z_S_F3x_S_vrr+WPX*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F3x_S_vrr = PAY*I_ERI_F3y_S_F3x_S_vrr+WPY*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_G3yz_S_F3x_S_vrr = PAZ*I_ERI_F3y_S_F3x_S_vrr+WPZ*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3x_S_vrr = PAZ*I_ERI_F2yz_S_F3x_S_vrr+WPZ*I_ERI_F2yz_S_F3x_S_M1_vrr+oned2z*I_ERI_D2y_S_F3x_S_vrr-rhod2zsq*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3x_S_vrr = PAY*I_ERI_F3z_S_F3x_S_vrr+WPY*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4z_S_F3x_S_vrr = PAZ*I_ERI_F3z_S_F3x_S_vrr+WPZ*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_G4x_S_F2xy_S_vrr = PAX*I_ERI_F3x_S_F2xy_S_vrr+WPX*I_ERI_F3x_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xy_S_vrr = PAY*I_ERI_F3x_S_F2xy_S_vrr+WPY*I_ERI_F3x_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xy_S_vrr = PAZ*I_ERI_F3x_S_F2xy_S_vrr+WPZ*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_vrr = PAY*I_ERI_F2xy_S_F2xy_S_vrr+WPY*I_ERI_F2xy_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xy_S_vrr = PAZ*I_ERI_F2xy_S_F2xy_S_vrr+WPZ*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xy_S_vrr = PAZ*I_ERI_F2xz_S_F2xy_S_vrr+WPZ*I_ERI_F2xz_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_vrr = PAX*I_ERI_F3y_S_F2xy_S_vrr+WPX*I_ERI_F3y_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xy_S_vrr = PAZ*I_ERI_Fx2y_S_F2xy_S_vrr+WPZ*I_ERI_Fx2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xy_S_vrr = PAY*I_ERI_Fx2z_S_F2xy_S_vrr+WPY*I_ERI_Fx2z_S_F2xy_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_vrr = PAX*I_ERI_F3z_S_F2xy_S_vrr+WPX*I_ERI_F3z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4y_S_F2xy_S_vrr = PAY*I_ERI_F3y_S_F2xy_S_vrr+WPY*I_ERI_F3y_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xy_S_vrr = PAZ*I_ERI_F3y_S_F2xy_S_vrr+WPZ*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xy_S_vrr = PAZ*I_ERI_F2yz_S_F2xy_S_vrr+WPZ*I_ERI_F2yz_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2y_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_vrr = PAY*I_ERI_F3z_S_F2xy_S_vrr+WPY*I_ERI_F3z_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xy_S_vrr = PAZ*I_ERI_F3z_S_F2xy_S_vrr+WPZ*I_ERI_F3z_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4x_S_F2xz_S_vrr = PAX*I_ERI_F3x_S_F2xz_S_vrr+WPX*I_ERI_F3x_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xz_S_vrr = PAY*I_ERI_F3x_S_F2xz_S_vrr+WPY*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xz_S_vrr = PAZ*I_ERI_F3x_S_F2xz_S_vrr+WPZ*I_ERI_F3x_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_vrr = PAY*I_ERI_F2xy_S_F2xz_S_vrr+WPY*I_ERI_F2xy_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xz_S_vrr = PAZ*I_ERI_F2xy_S_F2xz_S_vrr+WPZ*I_ERI_F2xy_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xz_S_vrr = PAZ*I_ERI_F2xz_S_F2xz_S_vrr+WPZ*I_ERI_F2xz_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_vrr = PAX*I_ERI_F3y_S_F2xz_S_vrr+WPX*I_ERI_F3y_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xz_S_vrr = PAZ*I_ERI_Fx2y_S_F2xz_S_vrr+WPZ*I_ERI_Fx2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xz_S_vrr = PAY*I_ERI_Fx2z_S_F2xz_S_vrr+WPY*I_ERI_Fx2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_vrr = PAX*I_ERI_F3z_S_F2xz_S_vrr+WPX*I_ERI_F3z_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4y_S_F2xz_S_vrr = PAY*I_ERI_F3y_S_F2xz_S_vrr+WPY*I_ERI_F3y_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xz_S_vrr = PAZ*I_ERI_F3y_S_F2xz_S_vrr+WPZ*I_ERI_F3y_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xz_S_vrr = PAZ*I_ERI_F2yz_S_F2xz_S_vrr+WPZ*I_ERI_F2yz_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2y_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_vrr = PAY*I_ERI_F3z_S_F2xz_S_vrr+WPY*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4z_S_F2xz_S_vrr = PAZ*I_ERI_F3z_S_F2xz_S_vrr+WPZ*I_ERI_F3z_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2y_S_vrr = PAX*I_ERI_F3x_S_Fx2y_S_vrr+WPX*I_ERI_F3x_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_vrr = PAY*I_ERI_F3x_S_Fx2y_S_vrr+WPY*I_ERI_F3x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_vrr = PAZ*I_ERI_F3x_S_Fx2y_S_vrr+WPZ*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_vrr = PAY*I_ERI_F2xy_S_Fx2y_S_vrr+WPY*I_ERI_F2xy_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2y_S_vrr = PAZ*I_ERI_F2xy_S_Fx2y_S_vrr+WPZ*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2y_S_vrr = PAZ*I_ERI_F2xz_S_Fx2y_S_vrr+WPZ*I_ERI_F2xz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_vrr = PAX*I_ERI_F3y_S_Fx2y_S_vrr+WPX*I_ERI_F3y_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2y_S_vrr = PAZ*I_ERI_Fx2y_S_Fx2y_S_vrr+WPZ*I_ERI_Fx2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2y_S_vrr = PAY*I_ERI_Fx2z_S_Fx2y_S_vrr+WPY*I_ERI_Fx2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_vrr = PAX*I_ERI_F3z_S_Fx2y_S_vrr+WPX*I_ERI_F3z_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2y_S_vrr = PAY*I_ERI_F3y_S_Fx2y_S_vrr+WPY*I_ERI_F3y_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_vrr = PAZ*I_ERI_F3y_S_Fx2y_S_vrr+WPZ*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2y_S_vrr = PAZ*I_ERI_F2yz_S_Fx2y_S_vrr+WPZ*I_ERI_F2yz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2y_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_vrr = PAY*I_ERI_F3z_S_Fx2y_S_vrr+WPY*I_ERI_F3z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2y_S_vrr = PAZ*I_ERI_F3z_S_Fx2y_S_vrr+WPZ*I_ERI_F3z_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fxyz_S_vrr = PAX*I_ERI_F3x_S_Fxyz_S_vrr+WPX*I_ERI_F3x_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_vrr = PAY*I_ERI_F3x_S_Fxyz_S_vrr+WPY*I_ERI_F3x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_vrr = PAZ*I_ERI_F3x_S_Fxyz_S_vrr+WPZ*I_ERI_F3x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_vrr = PAY*I_ERI_F2xy_S_Fxyz_S_vrr+WPY*I_ERI_F2xy_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fxyz_S_vrr = PAZ*I_ERI_F2xy_S_Fxyz_S_vrr+WPZ*I_ERI_F2xy_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fxyz_S_vrr = PAZ*I_ERI_F2xz_S_Fxyz_S_vrr+WPZ*I_ERI_F2xz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_vrr = PAX*I_ERI_F3y_S_Fxyz_S_vrr+WPX*I_ERI_F3y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fxyz_S_vrr = PAZ*I_ERI_Fx2y_S_Fxyz_S_vrr+WPZ*I_ERI_Fx2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fxyz_S_vrr = PAY*I_ERI_Fx2z_S_Fxyz_S_vrr+WPY*I_ERI_Fx2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_vrr = PAX*I_ERI_F3z_S_Fxyz_S_vrr+WPX*I_ERI_F3z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4y_S_Fxyz_S_vrr = PAY*I_ERI_F3y_S_Fxyz_S_vrr+WPY*I_ERI_F3y_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_vrr = PAZ*I_ERI_F3y_S_Fxyz_S_vrr+WPZ*I_ERI_F3y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fxyz_S_vrr = PAZ*I_ERI_F2yz_S_Fxyz_S_vrr+WPZ*I_ERI_F2yz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2y_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_vrr = PAY*I_ERI_F3z_S_Fxyz_S_vrr+WPY*I_ERI_F3z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4z_S_Fxyz_S_vrr = PAZ*I_ERI_F3z_S_Fxyz_S_vrr+WPZ*I_ERI_F3z_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2z_S_vrr = PAX*I_ERI_F3x_S_Fx2z_S_vrr+WPX*I_ERI_F3x_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_vrr = PAY*I_ERI_F3x_S_Fx2z_S_vrr+WPY*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_vrr = PAZ*I_ERI_F3x_S_Fx2z_S_vrr+WPZ*I_ERI_F3x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_vrr = PAY*I_ERI_F2xy_S_Fx2z_S_vrr+WPY*I_ERI_F2xy_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2z_S_vrr = PAZ*I_ERI_F2xy_S_Fx2z_S_vrr+WPZ*I_ERI_F2xy_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dxz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2z_S_vrr = PAZ*I_ERI_F2xz_S_Fx2z_S_vrr+WPZ*I_ERI_F2xz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_vrr = PAX*I_ERI_F3y_S_Fx2z_S_vrr+WPX*I_ERI_F3y_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2z_S_vrr = PAZ*I_ERI_Fx2y_S_Fx2z_S_vrr+WPZ*I_ERI_Fx2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2z_S_vrr = PAY*I_ERI_Fx2z_S_Fx2z_S_vrr+WPY*I_ERI_Fx2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_vrr = PAX*I_ERI_F3z_S_Fx2z_S_vrr+WPX*I_ERI_F3z_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2z_S_vrr = PAY*I_ERI_F3y_S_Fx2z_S_vrr+WPY*I_ERI_F3y_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_vrr = PAZ*I_ERI_F3y_S_Fx2z_S_vrr+WPZ*I_ERI_F3y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2z_S_vrr = PAZ*I_ERI_F2yz_S_Fx2z_S_vrr+WPZ*I_ERI_F2yz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dxz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_vrr = PAY*I_ERI_F3z_S_Fx2z_S_vrr+WPY*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2z_S_vrr = PAZ*I_ERI_F3z_S_Fx2z_S_vrr+WPZ*I_ERI_F3z_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4x_S_F3y_S_vrr = PAX*I_ERI_F3x_S_F3y_S_vrr+WPX*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xy_S_F3y_S_vrr = PAY*I_ERI_F3x_S_F3y_S_vrr+WPY*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F3y_S_vrr = PAZ*I_ERI_F3x_S_F3y_S_vrr+WPZ*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3y_S_vrr = PAY*I_ERI_F2xy_S_F3y_S_vrr+WPY*I_ERI_F2xy_S_F3y_S_M1_vrr+oned2z*I_ERI_D2x_S_F3y_S_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3y_S_vrr = PAZ*I_ERI_F2xy_S_F3y_S_vrr+WPZ*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3y_S_vrr = PAZ*I_ERI_F2xz_S_F3y_S_vrr+WPZ*I_ERI_F2xz_S_F3y_S_M1_vrr+oned2z*I_ERI_D2x_S_F3y_S_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3y_S_vrr = PAX*I_ERI_F3y_S_F3y_S_vrr+WPX*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3y_S_vrr = PAZ*I_ERI_Fx2y_S_F3y_S_vrr+WPZ*I_ERI_Fx2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3y_S_vrr = PAY*I_ERI_Fx2z_S_F3y_S_vrr+WPY*I_ERI_Fx2z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3y_S_vrr = PAX*I_ERI_F3z_S_F3y_S_vrr+WPX*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4y_S_F3y_S_vrr = PAY*I_ERI_F3y_S_F3y_S_vrr+WPY*I_ERI_F3y_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F3y_S_vrr = PAZ*I_ERI_F3y_S_F3y_S_vrr+WPZ*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3y_S_vrr = PAZ*I_ERI_F2yz_S_F3y_S_vrr+WPZ*I_ERI_F2yz_S_F3y_S_M1_vrr+oned2z*I_ERI_D2y_S_F3y_S_vrr-rhod2zsq*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3y_S_vrr = PAY*I_ERI_F3z_S_F3y_S_vrr+WPY*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F3y_S_vrr = PAZ*I_ERI_F3z_S_F3y_S_vrr+WPZ*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_G4x_S_F2yz_S_vrr = PAX*I_ERI_F3x_S_F2yz_S_vrr+WPX*I_ERI_F3x_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xy_S_F2yz_S_vrr = PAY*I_ERI_F3x_S_F2yz_S_vrr+WPY*I_ERI_F3x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G3xz_S_F2yz_S_vrr = PAZ*I_ERI_F3x_S_F2yz_S_vrr+WPZ*I_ERI_F3x_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_vrr = PAY*I_ERI_F2xy_S_F2yz_S_vrr+WPY*I_ERI_F2xy_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2yz_S_vrr = PAZ*I_ERI_F2xy_S_F2yz_S_vrr+WPZ*I_ERI_F2xy_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2yz_S_vrr = PAZ*I_ERI_F2xz_S_F2yz_S_vrr+WPZ*I_ERI_F2xz_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_vrr = PAX*I_ERI_F3y_S_F2yz_S_vrr+WPX*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2yz_S_vrr = PAZ*I_ERI_Fx2y_S_F2yz_S_vrr+WPZ*I_ERI_Fx2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2yz_S_vrr = PAY*I_ERI_Fx2z_S_F2yz_S_vrr+WPY*I_ERI_Fx2z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_vrr = PAX*I_ERI_F3z_S_F2yz_S_vrr+WPX*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4y_S_F2yz_S_vrr = PAY*I_ERI_F3y_S_F2yz_S_vrr+WPY*I_ERI_F3y_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_G3yz_S_F2yz_S_vrr = PAZ*I_ERI_F3y_S_F2yz_S_vrr+WPZ*I_ERI_F3y_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2yz_S_vrr = PAZ*I_ERI_F2yz_S_F2yz_S_vrr+WPZ*I_ERI_F2yz_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2y_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_vrr = PAY*I_ERI_F3z_S_F2yz_S_vrr+WPY*I_ERI_F3z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4z_S_F2yz_S_vrr = PAZ*I_ERI_F3z_S_F2yz_S_vrr+WPZ*I_ERI_F3z_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fy2z_S_vrr = PAX*I_ERI_F3x_S_Fy2z_S_vrr+WPX*I_ERI_F3x_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_vrr = PAY*I_ERI_F3x_S_Fy2z_S_vrr+WPY*I_ERI_F3x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_vrr = PAZ*I_ERI_F3x_S_Fy2z_S_vrr+WPZ*I_ERI_F3x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_vrr = PAY*I_ERI_F2xy_S_Fy2z_S_vrr+WPY*I_ERI_F2xy_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fy2z_S_vrr = PAZ*I_ERI_F2xy_S_Fy2z_S_vrr+WPZ*I_ERI_F2xy_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fy2z_S_vrr = PAZ*I_ERI_F2xz_S_Fy2z_S_vrr+WPZ*I_ERI_F2xz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_vrr = PAX*I_ERI_F3y_S_Fy2z_S_vrr+WPX*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fy2z_S_vrr = PAZ*I_ERI_Fx2y_S_Fy2z_S_vrr+WPZ*I_ERI_Fx2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fy2z_S_vrr = PAY*I_ERI_Fx2z_S_Fy2z_S_vrr+WPY*I_ERI_Fx2z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_vrr = PAX*I_ERI_F3z_S_Fy2z_S_vrr+WPX*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fy2z_S_vrr = PAY*I_ERI_F3y_S_Fy2z_S_vrr+WPY*I_ERI_F3y_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_vrr = PAZ*I_ERI_F3y_S_Fy2z_S_vrr+WPZ*I_ERI_F3y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fy2z_S_vrr = PAZ*I_ERI_F2yz_S_Fy2z_S_vrr+WPZ*I_ERI_F2yz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dyz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_vrr = PAY*I_ERI_F3z_S_Fy2z_S_vrr+WPY*I_ERI_F3z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fy2z_S_vrr = PAZ*I_ERI_F3z_S_Fy2z_S_vrr+WPZ*I_ERI_F3z_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4x_S_F3z_S_vrr = PAX*I_ERI_F3x_S_F3z_S_vrr+WPX*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xy_S_F3z_S_vrr = PAY*I_ERI_F3x_S_F3z_S_vrr+WPY*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xz_S_F3z_S_vrr = PAZ*I_ERI_F3x_S_F3z_S_vrr+WPZ*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3z_S_vrr = PAY*I_ERI_F2xy_S_F3z_S_vrr+WPY*I_ERI_F2xy_S_F3z_S_M1_vrr+oned2z*I_ERI_D2x_S_F3z_S_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3z_S_vrr = PAZ*I_ERI_F2xy_S_F3z_S_vrr+WPZ*I_ERI_F2xy_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3z_S_vrr = PAZ*I_ERI_F2xz_S_F3z_S_vrr+WPZ*I_ERI_F2xz_S_F3z_S_M1_vrr+oned2z*I_ERI_D2x_S_F3z_S_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3z_S_vrr = PAX*I_ERI_F3y_S_F3z_S_vrr+WPX*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3z_S_vrr = PAZ*I_ERI_Fx2y_S_F3z_S_vrr+WPZ*I_ERI_Fx2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3z_S_vrr = PAY*I_ERI_Fx2z_S_F3z_S_vrr+WPY*I_ERI_Fx2z_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3z_S_vrr = PAX*I_ERI_F3z_S_F3z_S_vrr+WPX*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4y_S_F3z_S_vrr = PAY*I_ERI_F3y_S_F3z_S_vrr+WPY*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_G3yz_S_F3z_S_vrr = PAZ*I_ERI_F3y_S_F3z_S_vrr+WPZ*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3z_S_vrr = PAZ*I_ERI_F2yz_S_F3z_S_vrr+WPZ*I_ERI_F2yz_S_F3z_S_M1_vrr+oned2z*I_ERI_D2y_S_F3z_S_vrr-rhod2zsq*I_ERI_D2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3z_S_vrr = PAY*I_ERI_F3z_S_F3z_S_vrr+WPY*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4z_S_F3z_S_vrr = PAZ*I_ERI_F3z_S_F3z_S_vrr+WPZ*I_ERI_F3z_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_F3x_S += I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S += I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S += I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S += I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S += I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S += I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S += I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S += I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S += I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S += I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S += I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S += I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S += I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S += I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S += I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S += I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S += I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S += I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S += I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S += I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S += I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S += I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S += I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S += I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S += I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S += I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S += I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S += I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S += I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S += I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S += I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S += I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S += I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S += I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S += I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S += I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S += I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S += I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S += I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S += I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S += I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S += I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S += I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S += I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S += I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S += I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S += I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S += I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S += I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S += I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S += I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S += I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S += I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S += I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S += I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S += I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S += I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S += I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S += I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S += I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_a_coefs = alpha;
      I_ERI_G4x_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_a += SQ_ERI_G_S_F_S_a_coefs*I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_a_coefs = alpha;
      I_ERI_F3x_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_a += SQ_ERI_F_S_F_S_a_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_F3x_S += I_ERI_Px_S_F3x_S_vrr;
      I_ERI_Py_S_F3x_S += I_ERI_Py_S_F3x_S_vrr;
      I_ERI_Pz_S_F3x_S += I_ERI_Pz_S_F3x_S_vrr;
      I_ERI_Px_S_F2xy_S += I_ERI_Px_S_F2xy_S_vrr;
      I_ERI_Py_S_F2xy_S += I_ERI_Py_S_F2xy_S_vrr;
      I_ERI_Pz_S_F2xy_S += I_ERI_Pz_S_F2xy_S_vrr;
      I_ERI_Px_S_F2xz_S += I_ERI_Px_S_F2xz_S_vrr;
      I_ERI_Py_S_F2xz_S += I_ERI_Py_S_F2xz_S_vrr;
      I_ERI_Pz_S_F2xz_S += I_ERI_Pz_S_F2xz_S_vrr;
      I_ERI_Px_S_Fx2y_S += I_ERI_Px_S_Fx2y_S_vrr;
      I_ERI_Py_S_Fx2y_S += I_ERI_Py_S_Fx2y_S_vrr;
      I_ERI_Pz_S_Fx2y_S += I_ERI_Pz_S_Fx2y_S_vrr;
      I_ERI_Px_S_Fxyz_S += I_ERI_Px_S_Fxyz_S_vrr;
      I_ERI_Py_S_Fxyz_S += I_ERI_Py_S_Fxyz_S_vrr;
      I_ERI_Pz_S_Fxyz_S += I_ERI_Pz_S_Fxyz_S_vrr;
      I_ERI_Px_S_Fx2z_S += I_ERI_Px_S_Fx2z_S_vrr;
      I_ERI_Py_S_Fx2z_S += I_ERI_Py_S_Fx2z_S_vrr;
      I_ERI_Pz_S_Fx2z_S += I_ERI_Pz_S_Fx2z_S_vrr;
      I_ERI_Px_S_F3y_S += I_ERI_Px_S_F3y_S_vrr;
      I_ERI_Py_S_F3y_S += I_ERI_Py_S_F3y_S_vrr;
      I_ERI_Pz_S_F3y_S += I_ERI_Pz_S_F3y_S_vrr;
      I_ERI_Px_S_F2yz_S += I_ERI_Px_S_F2yz_S_vrr;
      I_ERI_Py_S_F2yz_S += I_ERI_Py_S_F2yz_S_vrr;
      I_ERI_Pz_S_F2yz_S += I_ERI_Pz_S_F2yz_S_vrr;
      I_ERI_Px_S_Fy2z_S += I_ERI_Px_S_Fy2z_S_vrr;
      I_ERI_Py_S_Fy2z_S += I_ERI_Py_S_Fy2z_S_vrr;
      I_ERI_Pz_S_Fy2z_S += I_ERI_Pz_S_Fy2z_S_vrr;
      I_ERI_Px_S_F3z_S += I_ERI_Px_S_F3z_S_vrr;
      I_ERI_Py_S_F3z_S += I_ERI_Py_S_F3z_S_vrr;
      I_ERI_Pz_S_F3z_S += I_ERI_Pz_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_b_coefs = beta;
      I_ERI_G4x_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_b += SQ_ERI_G_S_F_S_b_coefs*I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_b_coefs = beta;
      I_ERI_F3x_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_b += SQ_ERI_F_S_F_S_b_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_b_coefs = beta;
      I_ERI_D2x_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_b += SQ_ERI_D_S_F_S_b_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_G_S_d_coefs = delta;
      I_ERI_F3x_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G4x_S_vrr;
      I_ERI_F2xy_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G4x_S_vrr;
      I_ERI_F2xz_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G4x_S_vrr;
      I_ERI_Fx2y_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G4x_S_vrr;
      I_ERI_Fxyz_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G4x_S_vrr;
      I_ERI_Fx2z_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G4x_S_vrr;
      I_ERI_F3y_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G4x_S_vrr;
      I_ERI_F2yz_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G4x_S_vrr;
      I_ERI_Fy2z_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G4x_S_vrr;
      I_ERI_F3z_S_G4x_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G4x_S_vrr;
      I_ERI_F3x_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G3xy_S_vrr;
      I_ERI_F2xy_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G3xy_S_vrr;
      I_ERI_F2xz_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G3xy_S_vrr;
      I_ERI_Fx2y_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G3xy_S_vrr;
      I_ERI_Fxyz_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G3xy_S_vrr;
      I_ERI_Fx2z_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G3xy_S_vrr;
      I_ERI_F3y_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G3xy_S_vrr;
      I_ERI_F2yz_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G3xy_S_vrr;
      I_ERI_Fy2z_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G3xy_S_vrr;
      I_ERI_F3z_S_G3xy_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G3xy_S_vrr;
      I_ERI_F3x_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G3xz_S_vrr;
      I_ERI_F2xy_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G3xz_S_vrr;
      I_ERI_F2xz_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G3xz_S_vrr;
      I_ERI_Fx2y_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G3xz_S_vrr;
      I_ERI_Fxyz_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G3xz_S_vrr;
      I_ERI_Fx2z_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G3xz_S_vrr;
      I_ERI_F3y_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G3xz_S_vrr;
      I_ERI_F2yz_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G3xz_S_vrr;
      I_ERI_Fy2z_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G3xz_S_vrr;
      I_ERI_F3z_S_G3xz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G3xz_S_vrr;
      I_ERI_F3x_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G2x2y_S_vrr;
      I_ERI_F2xy_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G2x2y_S_vrr;
      I_ERI_F2xz_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G2x2y_S_vrr;
      I_ERI_Fx2y_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G2x2y_S_vrr;
      I_ERI_Fxyz_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G2x2y_S_vrr;
      I_ERI_Fx2z_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G2x2y_S_vrr;
      I_ERI_F3y_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G2x2y_S_vrr;
      I_ERI_F2yz_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G2x2y_S_vrr;
      I_ERI_Fy2z_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G2x2y_S_vrr;
      I_ERI_F3z_S_G2x2y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G2x2y_S_vrr;
      I_ERI_F3x_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G2xyz_S_vrr;
      I_ERI_F2xy_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G2xyz_S_vrr;
      I_ERI_F2xz_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G2xyz_S_vrr;
      I_ERI_Fx2y_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G2xyz_S_vrr;
      I_ERI_Fxyz_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G2xyz_S_vrr;
      I_ERI_Fx2z_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G2xyz_S_vrr;
      I_ERI_F3y_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G2xyz_S_vrr;
      I_ERI_F2yz_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G2xyz_S_vrr;
      I_ERI_Fy2z_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G2xyz_S_vrr;
      I_ERI_F3z_S_G2xyz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G2xyz_S_vrr;
      I_ERI_F3x_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G2x2z_S_vrr;
      I_ERI_F2xy_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G2x2z_S_vrr;
      I_ERI_F2xz_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G2x2z_S_vrr;
      I_ERI_Fx2y_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G2x2z_S_vrr;
      I_ERI_Fxyz_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G2x2z_S_vrr;
      I_ERI_Fx2z_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G2x2z_S_vrr;
      I_ERI_F3y_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G2x2z_S_vrr;
      I_ERI_F2yz_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G2x2z_S_vrr;
      I_ERI_Fy2z_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G2x2z_S_vrr;
      I_ERI_F3z_S_G2x2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G2x2z_S_vrr;
      I_ERI_F3x_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_Gx3y_S_vrr;
      I_ERI_F2xy_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_Gx3y_S_vrr;
      I_ERI_F2xz_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_Gx3y_S_vrr;
      I_ERI_Fx2y_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_Gx3y_S_vrr;
      I_ERI_Fxyz_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_Gx3y_S_vrr;
      I_ERI_Fx2z_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_Gx3y_S_vrr;
      I_ERI_F3y_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_Gx3y_S_vrr;
      I_ERI_F2yz_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_Gx3y_S_vrr;
      I_ERI_Fy2z_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_Gx3y_S_vrr;
      I_ERI_F3z_S_Gx3y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_Gx3y_S_vrr;
      I_ERI_F3x_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_Gx2yz_S_vrr;
      I_ERI_F2xy_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_Gx2yz_S_vrr;
      I_ERI_F2xz_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_Gx2yz_S_vrr;
      I_ERI_Fx2y_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_Gx2yz_S_vrr;
      I_ERI_Fxyz_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_Gx2yz_S_vrr;
      I_ERI_Fx2z_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_Gx2yz_S_vrr;
      I_ERI_F3y_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_Gx2yz_S_vrr;
      I_ERI_F2yz_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_Gx2yz_S_vrr;
      I_ERI_Fy2z_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_Gx2yz_S_vrr;
      I_ERI_F3z_S_Gx2yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_Gx2yz_S_vrr;
      I_ERI_F3x_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_Gxy2z_S_vrr;
      I_ERI_F2xy_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_Gxy2z_S_vrr;
      I_ERI_F2xz_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_Gxy2z_S_vrr;
      I_ERI_Fx2y_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_Gxy2z_S_vrr;
      I_ERI_Fxyz_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_Gxy2z_S_vrr;
      I_ERI_Fx2z_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_Gxy2z_S_vrr;
      I_ERI_F3y_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_Gxy2z_S_vrr;
      I_ERI_F2yz_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_Gxy2z_S_vrr;
      I_ERI_Fy2z_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_Gxy2z_S_vrr;
      I_ERI_F3z_S_Gxy2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_Gxy2z_S_vrr;
      I_ERI_F3x_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_Gx3z_S_vrr;
      I_ERI_F2xy_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_Gx3z_S_vrr;
      I_ERI_F2xz_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_Gx3z_S_vrr;
      I_ERI_Fx2y_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_Gx3z_S_vrr;
      I_ERI_Fxyz_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_Gx3z_S_vrr;
      I_ERI_Fx2z_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_Gx3z_S_vrr;
      I_ERI_F3y_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_Gx3z_S_vrr;
      I_ERI_F2yz_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_Gx3z_S_vrr;
      I_ERI_Fy2z_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_Gx3z_S_vrr;
      I_ERI_F3z_S_Gx3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_Gx3z_S_vrr;
      I_ERI_F3x_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G4y_S_vrr;
      I_ERI_F2xy_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G4y_S_vrr;
      I_ERI_F2xz_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G4y_S_vrr;
      I_ERI_Fx2y_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G4y_S_vrr;
      I_ERI_Fxyz_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G4y_S_vrr;
      I_ERI_Fx2z_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G4y_S_vrr;
      I_ERI_F3y_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G4y_S_vrr;
      I_ERI_F2yz_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G4y_S_vrr;
      I_ERI_Fy2z_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G4y_S_vrr;
      I_ERI_F3z_S_G4y_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G4y_S_vrr;
      I_ERI_F3x_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G3yz_S_vrr;
      I_ERI_F2xy_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G3yz_S_vrr;
      I_ERI_F2xz_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G3yz_S_vrr;
      I_ERI_Fx2y_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G3yz_S_vrr;
      I_ERI_Fxyz_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G3yz_S_vrr;
      I_ERI_Fx2z_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G3yz_S_vrr;
      I_ERI_F3y_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G3yz_S_vrr;
      I_ERI_F2yz_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G3yz_S_vrr;
      I_ERI_Fy2z_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G3yz_S_vrr;
      I_ERI_F3z_S_G3yz_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G3yz_S_vrr;
      I_ERI_F3x_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G2y2z_S_vrr;
      I_ERI_F2xy_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G2y2z_S_vrr;
      I_ERI_F2xz_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G2y2z_S_vrr;
      I_ERI_Fx2y_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G2y2z_S_vrr;
      I_ERI_Fxyz_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G2y2z_S_vrr;
      I_ERI_Fx2z_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G2y2z_S_vrr;
      I_ERI_F3y_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G2y2z_S_vrr;
      I_ERI_F2yz_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G2y2z_S_vrr;
      I_ERI_Fy2z_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G2y2z_S_vrr;
      I_ERI_F3z_S_G2y2z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G2y2z_S_vrr;
      I_ERI_F3x_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_Gy3z_S_vrr;
      I_ERI_F2xy_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_Gy3z_S_vrr;
      I_ERI_F2xz_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_Gy3z_S_vrr;
      I_ERI_Fx2y_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_Gy3z_S_vrr;
      I_ERI_Fxyz_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_Gy3z_S_vrr;
      I_ERI_Fx2z_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_Gy3z_S_vrr;
      I_ERI_F3y_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_Gy3z_S_vrr;
      I_ERI_F2yz_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_Gy3z_S_vrr;
      I_ERI_Fy2z_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_Gy3z_S_vrr;
      I_ERI_F3z_S_Gy3z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_Gy3z_S_vrr;
      I_ERI_F3x_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3x_S_G4z_S_vrr;
      I_ERI_F2xy_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xy_S_G4z_S_vrr;
      I_ERI_F2xz_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2xz_S_G4z_S_vrr;
      I_ERI_Fx2y_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2y_S_G4z_S_vrr;
      I_ERI_Fxyz_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fxyz_S_G4z_S_vrr;
      I_ERI_Fx2z_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fx2z_S_G4z_S_vrr;
      I_ERI_F3y_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3y_S_G4z_S_vrr;
      I_ERI_F2yz_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F2yz_S_G4z_S_vrr;
      I_ERI_Fy2z_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_Fy2z_S_G4z_S_vrr;
      I_ERI_F3z_S_G4z_S_d += SQ_ERI_F_S_G_S_d_coefs*I_ERI_F3z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_d_coefs = delta;
      I_ERI_F3x_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_d += SQ_ERI_F_S_F_S_d_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_G_S_d_coefs = delta;
      I_ERI_D2x_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G4x_S_vrr;
      I_ERI_Dxy_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G4x_S_vrr;
      I_ERI_Dxz_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G4x_S_vrr;
      I_ERI_D2y_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G4x_S_vrr;
      I_ERI_Dyz_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G4x_S_vrr;
      I_ERI_D2z_S_G4x_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G4x_S_vrr;
      I_ERI_D2x_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G3xy_S_vrr;
      I_ERI_Dxy_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G3xy_S_vrr;
      I_ERI_Dxz_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G3xy_S_vrr;
      I_ERI_D2y_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G3xy_S_vrr;
      I_ERI_Dyz_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G3xy_S_vrr;
      I_ERI_D2z_S_G3xy_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G3xy_S_vrr;
      I_ERI_D2x_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G3xz_S_vrr;
      I_ERI_Dxy_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G3xz_S_vrr;
      I_ERI_Dxz_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G3xz_S_vrr;
      I_ERI_D2y_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G3xz_S_vrr;
      I_ERI_Dyz_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G3xz_S_vrr;
      I_ERI_D2z_S_G3xz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G3xz_S_vrr;
      I_ERI_D2x_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G2x2y_S_vrr;
      I_ERI_Dxy_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G2x2y_S_vrr;
      I_ERI_Dxz_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G2x2y_S_vrr;
      I_ERI_D2y_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G2x2y_S_vrr;
      I_ERI_Dyz_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G2x2y_S_vrr;
      I_ERI_D2z_S_G2x2y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G2x2y_S_vrr;
      I_ERI_D2x_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G2xyz_S_vrr;
      I_ERI_Dxy_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G2xyz_S_vrr;
      I_ERI_Dxz_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G2xyz_S_vrr;
      I_ERI_D2y_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G2xyz_S_vrr;
      I_ERI_Dyz_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G2xyz_S_vrr;
      I_ERI_D2z_S_G2xyz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G2xyz_S_vrr;
      I_ERI_D2x_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G2x2z_S_vrr;
      I_ERI_Dxy_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G2x2z_S_vrr;
      I_ERI_Dxz_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G2x2z_S_vrr;
      I_ERI_D2y_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G2x2z_S_vrr;
      I_ERI_Dyz_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G2x2z_S_vrr;
      I_ERI_D2z_S_G2x2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G2x2z_S_vrr;
      I_ERI_D2x_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_Gx3y_S_vrr;
      I_ERI_Dxy_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_Gx3y_S_vrr;
      I_ERI_Dxz_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_Gx3y_S_vrr;
      I_ERI_D2y_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_Gx3y_S_vrr;
      I_ERI_Dyz_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_Gx3y_S_vrr;
      I_ERI_D2z_S_Gx3y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_Gx3y_S_vrr;
      I_ERI_D2x_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_Gx2yz_S_vrr;
      I_ERI_Dxy_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_Gx2yz_S_vrr;
      I_ERI_Dxz_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_Gx2yz_S_vrr;
      I_ERI_D2y_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_Gx2yz_S_vrr;
      I_ERI_Dyz_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_Gx2yz_S_vrr;
      I_ERI_D2z_S_Gx2yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_Gx2yz_S_vrr;
      I_ERI_D2x_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_Gxy2z_S_vrr;
      I_ERI_Dxy_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_Gxy2z_S_vrr;
      I_ERI_Dxz_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_Gxy2z_S_vrr;
      I_ERI_D2y_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_Gxy2z_S_vrr;
      I_ERI_Dyz_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_Gxy2z_S_vrr;
      I_ERI_D2z_S_Gxy2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_Gxy2z_S_vrr;
      I_ERI_D2x_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_Gx3z_S_vrr;
      I_ERI_Dxy_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_Gx3z_S_vrr;
      I_ERI_Dxz_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_Gx3z_S_vrr;
      I_ERI_D2y_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_Gx3z_S_vrr;
      I_ERI_Dyz_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_Gx3z_S_vrr;
      I_ERI_D2z_S_Gx3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_Gx3z_S_vrr;
      I_ERI_D2x_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G4y_S_vrr;
      I_ERI_Dxy_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G4y_S_vrr;
      I_ERI_Dxz_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G4y_S_vrr;
      I_ERI_D2y_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G4y_S_vrr;
      I_ERI_Dyz_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G4y_S_vrr;
      I_ERI_D2z_S_G4y_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G4y_S_vrr;
      I_ERI_D2x_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G3yz_S_vrr;
      I_ERI_Dxy_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G3yz_S_vrr;
      I_ERI_Dxz_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G3yz_S_vrr;
      I_ERI_D2y_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G3yz_S_vrr;
      I_ERI_Dyz_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G3yz_S_vrr;
      I_ERI_D2z_S_G3yz_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G3yz_S_vrr;
      I_ERI_D2x_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G2y2z_S_vrr;
      I_ERI_Dxy_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G2y2z_S_vrr;
      I_ERI_Dxz_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G2y2z_S_vrr;
      I_ERI_D2y_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G2y2z_S_vrr;
      I_ERI_Dyz_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G2y2z_S_vrr;
      I_ERI_D2z_S_G2y2z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G2y2z_S_vrr;
      I_ERI_D2x_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_Gy3z_S_vrr;
      I_ERI_Dxy_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_Gy3z_S_vrr;
      I_ERI_Dxz_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_Gy3z_S_vrr;
      I_ERI_D2y_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_Gy3z_S_vrr;
      I_ERI_Dyz_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_Gy3z_S_vrr;
      I_ERI_D2z_S_Gy3z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_Gy3z_S_vrr;
      I_ERI_D2x_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2x_S_G4z_S_vrr;
      I_ERI_Dxy_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxy_S_G4z_S_vrr;
      I_ERI_Dxz_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dxz_S_G4z_S_vrr;
      I_ERI_D2y_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2y_S_G4z_S_vrr;
      I_ERI_Dyz_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_Dyz_S_G4z_S_vrr;
      I_ERI_D2z_S_G4z_S_d += SQ_ERI_D_S_G_S_d_coefs*I_ERI_D2z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_d
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_d_coefs = delta;
      I_ERI_D2x_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_d += SQ_ERI_D_S_F_S_d_coefs*I_ERI_D2z_S_F3z_S_vrr;
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
   * shell quartet name: SQ_ERI_D_S_F_P_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_G_S_d
   * RHS shell quartet name: SQ_ERI_D_S_F_S_d
   ************************************************************/
  Double I_ERI_D2x_S_F3x_Px_d = I_ERI_D2x_S_G4x_S_d+CDX*I_ERI_D2x_S_F3x_S_d;
  Double I_ERI_Dxy_S_F3x_Px_d = I_ERI_Dxy_S_G4x_S_d+CDX*I_ERI_Dxy_S_F3x_S_d;
  Double I_ERI_Dxz_S_F3x_Px_d = I_ERI_Dxz_S_G4x_S_d+CDX*I_ERI_Dxz_S_F3x_S_d;
  Double I_ERI_D2y_S_F3x_Px_d = I_ERI_D2y_S_G4x_S_d+CDX*I_ERI_D2y_S_F3x_S_d;
  Double I_ERI_Dyz_S_F3x_Px_d = I_ERI_Dyz_S_G4x_S_d+CDX*I_ERI_Dyz_S_F3x_S_d;
  Double I_ERI_D2z_S_F3x_Px_d = I_ERI_D2z_S_G4x_S_d+CDX*I_ERI_D2z_S_F3x_S_d;
  Double I_ERI_D2x_S_F2xy_Px_d = I_ERI_D2x_S_G3xy_S_d+CDX*I_ERI_D2x_S_F2xy_S_d;
  Double I_ERI_Dxy_S_F2xy_Px_d = I_ERI_Dxy_S_G3xy_S_d+CDX*I_ERI_Dxy_S_F2xy_S_d;
  Double I_ERI_Dxz_S_F2xy_Px_d = I_ERI_Dxz_S_G3xy_S_d+CDX*I_ERI_Dxz_S_F2xy_S_d;
  Double I_ERI_D2y_S_F2xy_Px_d = I_ERI_D2y_S_G3xy_S_d+CDX*I_ERI_D2y_S_F2xy_S_d;
  Double I_ERI_Dyz_S_F2xy_Px_d = I_ERI_Dyz_S_G3xy_S_d+CDX*I_ERI_Dyz_S_F2xy_S_d;
  Double I_ERI_D2z_S_F2xy_Px_d = I_ERI_D2z_S_G3xy_S_d+CDX*I_ERI_D2z_S_F2xy_S_d;
  Double I_ERI_D2x_S_F2xz_Px_d = I_ERI_D2x_S_G3xz_S_d+CDX*I_ERI_D2x_S_F2xz_S_d;
  Double I_ERI_Dxy_S_F2xz_Px_d = I_ERI_Dxy_S_G3xz_S_d+CDX*I_ERI_Dxy_S_F2xz_S_d;
  Double I_ERI_Dxz_S_F2xz_Px_d = I_ERI_Dxz_S_G3xz_S_d+CDX*I_ERI_Dxz_S_F2xz_S_d;
  Double I_ERI_D2y_S_F2xz_Px_d = I_ERI_D2y_S_G3xz_S_d+CDX*I_ERI_D2y_S_F2xz_S_d;
  Double I_ERI_Dyz_S_F2xz_Px_d = I_ERI_Dyz_S_G3xz_S_d+CDX*I_ERI_Dyz_S_F2xz_S_d;
  Double I_ERI_D2z_S_F2xz_Px_d = I_ERI_D2z_S_G3xz_S_d+CDX*I_ERI_D2z_S_F2xz_S_d;
  Double I_ERI_D2x_S_Fx2y_Px_d = I_ERI_D2x_S_G2x2y_S_d+CDX*I_ERI_D2x_S_Fx2y_S_d;
  Double I_ERI_Dxy_S_Fx2y_Px_d = I_ERI_Dxy_S_G2x2y_S_d+CDX*I_ERI_Dxy_S_Fx2y_S_d;
  Double I_ERI_Dxz_S_Fx2y_Px_d = I_ERI_Dxz_S_G2x2y_S_d+CDX*I_ERI_Dxz_S_Fx2y_S_d;
  Double I_ERI_D2y_S_Fx2y_Px_d = I_ERI_D2y_S_G2x2y_S_d+CDX*I_ERI_D2y_S_Fx2y_S_d;
  Double I_ERI_Dyz_S_Fx2y_Px_d = I_ERI_Dyz_S_G2x2y_S_d+CDX*I_ERI_Dyz_S_Fx2y_S_d;
  Double I_ERI_D2z_S_Fx2y_Px_d = I_ERI_D2z_S_G2x2y_S_d+CDX*I_ERI_D2z_S_Fx2y_S_d;
  Double I_ERI_D2x_S_Fxyz_Px_d = I_ERI_D2x_S_G2xyz_S_d+CDX*I_ERI_D2x_S_Fxyz_S_d;
  Double I_ERI_Dxy_S_Fxyz_Px_d = I_ERI_Dxy_S_G2xyz_S_d+CDX*I_ERI_Dxy_S_Fxyz_S_d;
  Double I_ERI_Dxz_S_Fxyz_Px_d = I_ERI_Dxz_S_G2xyz_S_d+CDX*I_ERI_Dxz_S_Fxyz_S_d;
  Double I_ERI_D2y_S_Fxyz_Px_d = I_ERI_D2y_S_G2xyz_S_d+CDX*I_ERI_D2y_S_Fxyz_S_d;
  Double I_ERI_Dyz_S_Fxyz_Px_d = I_ERI_Dyz_S_G2xyz_S_d+CDX*I_ERI_Dyz_S_Fxyz_S_d;
  Double I_ERI_D2z_S_Fxyz_Px_d = I_ERI_D2z_S_G2xyz_S_d+CDX*I_ERI_D2z_S_Fxyz_S_d;
  Double I_ERI_D2x_S_Fx2z_Px_d = I_ERI_D2x_S_G2x2z_S_d+CDX*I_ERI_D2x_S_Fx2z_S_d;
  Double I_ERI_Dxy_S_Fx2z_Px_d = I_ERI_Dxy_S_G2x2z_S_d+CDX*I_ERI_Dxy_S_Fx2z_S_d;
  Double I_ERI_Dxz_S_Fx2z_Px_d = I_ERI_Dxz_S_G2x2z_S_d+CDX*I_ERI_Dxz_S_Fx2z_S_d;
  Double I_ERI_D2y_S_Fx2z_Px_d = I_ERI_D2y_S_G2x2z_S_d+CDX*I_ERI_D2y_S_Fx2z_S_d;
  Double I_ERI_Dyz_S_Fx2z_Px_d = I_ERI_Dyz_S_G2x2z_S_d+CDX*I_ERI_Dyz_S_Fx2z_S_d;
  Double I_ERI_D2z_S_Fx2z_Px_d = I_ERI_D2z_S_G2x2z_S_d+CDX*I_ERI_D2z_S_Fx2z_S_d;
  Double I_ERI_D2x_S_F3y_Px_d = I_ERI_D2x_S_Gx3y_S_d+CDX*I_ERI_D2x_S_F3y_S_d;
  Double I_ERI_Dxy_S_F3y_Px_d = I_ERI_Dxy_S_Gx3y_S_d+CDX*I_ERI_Dxy_S_F3y_S_d;
  Double I_ERI_Dxz_S_F3y_Px_d = I_ERI_Dxz_S_Gx3y_S_d+CDX*I_ERI_Dxz_S_F3y_S_d;
  Double I_ERI_D2y_S_F3y_Px_d = I_ERI_D2y_S_Gx3y_S_d+CDX*I_ERI_D2y_S_F3y_S_d;
  Double I_ERI_Dyz_S_F3y_Px_d = I_ERI_Dyz_S_Gx3y_S_d+CDX*I_ERI_Dyz_S_F3y_S_d;
  Double I_ERI_D2z_S_F3y_Px_d = I_ERI_D2z_S_Gx3y_S_d+CDX*I_ERI_D2z_S_F3y_S_d;
  Double I_ERI_D2x_S_F2yz_Px_d = I_ERI_D2x_S_Gx2yz_S_d+CDX*I_ERI_D2x_S_F2yz_S_d;
  Double I_ERI_Dxy_S_F2yz_Px_d = I_ERI_Dxy_S_Gx2yz_S_d+CDX*I_ERI_Dxy_S_F2yz_S_d;
  Double I_ERI_Dxz_S_F2yz_Px_d = I_ERI_Dxz_S_Gx2yz_S_d+CDX*I_ERI_Dxz_S_F2yz_S_d;
  Double I_ERI_D2y_S_F2yz_Px_d = I_ERI_D2y_S_Gx2yz_S_d+CDX*I_ERI_D2y_S_F2yz_S_d;
  Double I_ERI_Dyz_S_F2yz_Px_d = I_ERI_Dyz_S_Gx2yz_S_d+CDX*I_ERI_Dyz_S_F2yz_S_d;
  Double I_ERI_D2z_S_F2yz_Px_d = I_ERI_D2z_S_Gx2yz_S_d+CDX*I_ERI_D2z_S_F2yz_S_d;
  Double I_ERI_D2x_S_Fy2z_Px_d = I_ERI_D2x_S_Gxy2z_S_d+CDX*I_ERI_D2x_S_Fy2z_S_d;
  Double I_ERI_Dxy_S_Fy2z_Px_d = I_ERI_Dxy_S_Gxy2z_S_d+CDX*I_ERI_Dxy_S_Fy2z_S_d;
  Double I_ERI_Dxz_S_Fy2z_Px_d = I_ERI_Dxz_S_Gxy2z_S_d+CDX*I_ERI_Dxz_S_Fy2z_S_d;
  Double I_ERI_D2y_S_Fy2z_Px_d = I_ERI_D2y_S_Gxy2z_S_d+CDX*I_ERI_D2y_S_Fy2z_S_d;
  Double I_ERI_Dyz_S_Fy2z_Px_d = I_ERI_Dyz_S_Gxy2z_S_d+CDX*I_ERI_Dyz_S_Fy2z_S_d;
  Double I_ERI_D2z_S_Fy2z_Px_d = I_ERI_D2z_S_Gxy2z_S_d+CDX*I_ERI_D2z_S_Fy2z_S_d;
  Double I_ERI_D2x_S_F3z_Px_d = I_ERI_D2x_S_Gx3z_S_d+CDX*I_ERI_D2x_S_F3z_S_d;
  Double I_ERI_Dxy_S_F3z_Px_d = I_ERI_Dxy_S_Gx3z_S_d+CDX*I_ERI_Dxy_S_F3z_S_d;
  Double I_ERI_Dxz_S_F3z_Px_d = I_ERI_Dxz_S_Gx3z_S_d+CDX*I_ERI_Dxz_S_F3z_S_d;
  Double I_ERI_D2y_S_F3z_Px_d = I_ERI_D2y_S_Gx3z_S_d+CDX*I_ERI_D2y_S_F3z_S_d;
  Double I_ERI_Dyz_S_F3z_Px_d = I_ERI_Dyz_S_Gx3z_S_d+CDX*I_ERI_Dyz_S_F3z_S_d;
  Double I_ERI_D2z_S_F3z_Px_d = I_ERI_D2z_S_Gx3z_S_d+CDX*I_ERI_D2z_S_F3z_S_d;
  Double I_ERI_D2x_S_F3x_Py_d = I_ERI_D2x_S_G3xy_S_d+CDY*I_ERI_D2x_S_F3x_S_d;
  Double I_ERI_Dxy_S_F3x_Py_d = I_ERI_Dxy_S_G3xy_S_d+CDY*I_ERI_Dxy_S_F3x_S_d;
  Double I_ERI_Dxz_S_F3x_Py_d = I_ERI_Dxz_S_G3xy_S_d+CDY*I_ERI_Dxz_S_F3x_S_d;
  Double I_ERI_D2y_S_F3x_Py_d = I_ERI_D2y_S_G3xy_S_d+CDY*I_ERI_D2y_S_F3x_S_d;
  Double I_ERI_Dyz_S_F3x_Py_d = I_ERI_Dyz_S_G3xy_S_d+CDY*I_ERI_Dyz_S_F3x_S_d;
  Double I_ERI_D2z_S_F3x_Py_d = I_ERI_D2z_S_G3xy_S_d+CDY*I_ERI_D2z_S_F3x_S_d;
  Double I_ERI_D2x_S_F2xy_Py_d = I_ERI_D2x_S_G2x2y_S_d+CDY*I_ERI_D2x_S_F2xy_S_d;
  Double I_ERI_Dxy_S_F2xy_Py_d = I_ERI_Dxy_S_G2x2y_S_d+CDY*I_ERI_Dxy_S_F2xy_S_d;
  Double I_ERI_Dxz_S_F2xy_Py_d = I_ERI_Dxz_S_G2x2y_S_d+CDY*I_ERI_Dxz_S_F2xy_S_d;
  Double I_ERI_D2y_S_F2xy_Py_d = I_ERI_D2y_S_G2x2y_S_d+CDY*I_ERI_D2y_S_F2xy_S_d;
  Double I_ERI_Dyz_S_F2xy_Py_d = I_ERI_Dyz_S_G2x2y_S_d+CDY*I_ERI_Dyz_S_F2xy_S_d;
  Double I_ERI_D2z_S_F2xy_Py_d = I_ERI_D2z_S_G2x2y_S_d+CDY*I_ERI_D2z_S_F2xy_S_d;
  Double I_ERI_D2x_S_F2xz_Py_d = I_ERI_D2x_S_G2xyz_S_d+CDY*I_ERI_D2x_S_F2xz_S_d;
  Double I_ERI_Dxy_S_F2xz_Py_d = I_ERI_Dxy_S_G2xyz_S_d+CDY*I_ERI_Dxy_S_F2xz_S_d;
  Double I_ERI_Dxz_S_F2xz_Py_d = I_ERI_Dxz_S_G2xyz_S_d+CDY*I_ERI_Dxz_S_F2xz_S_d;
  Double I_ERI_D2y_S_F2xz_Py_d = I_ERI_D2y_S_G2xyz_S_d+CDY*I_ERI_D2y_S_F2xz_S_d;
  Double I_ERI_Dyz_S_F2xz_Py_d = I_ERI_Dyz_S_G2xyz_S_d+CDY*I_ERI_Dyz_S_F2xz_S_d;
  Double I_ERI_D2z_S_F2xz_Py_d = I_ERI_D2z_S_G2xyz_S_d+CDY*I_ERI_D2z_S_F2xz_S_d;
  Double I_ERI_D2x_S_Fx2y_Py_d = I_ERI_D2x_S_Gx3y_S_d+CDY*I_ERI_D2x_S_Fx2y_S_d;
  Double I_ERI_Dxy_S_Fx2y_Py_d = I_ERI_Dxy_S_Gx3y_S_d+CDY*I_ERI_Dxy_S_Fx2y_S_d;
  Double I_ERI_Dxz_S_Fx2y_Py_d = I_ERI_Dxz_S_Gx3y_S_d+CDY*I_ERI_Dxz_S_Fx2y_S_d;
  Double I_ERI_D2y_S_Fx2y_Py_d = I_ERI_D2y_S_Gx3y_S_d+CDY*I_ERI_D2y_S_Fx2y_S_d;
  Double I_ERI_Dyz_S_Fx2y_Py_d = I_ERI_Dyz_S_Gx3y_S_d+CDY*I_ERI_Dyz_S_Fx2y_S_d;
  Double I_ERI_D2z_S_Fx2y_Py_d = I_ERI_D2z_S_Gx3y_S_d+CDY*I_ERI_D2z_S_Fx2y_S_d;
  Double I_ERI_D2x_S_Fxyz_Py_d = I_ERI_D2x_S_Gx2yz_S_d+CDY*I_ERI_D2x_S_Fxyz_S_d;
  Double I_ERI_Dxy_S_Fxyz_Py_d = I_ERI_Dxy_S_Gx2yz_S_d+CDY*I_ERI_Dxy_S_Fxyz_S_d;
  Double I_ERI_Dxz_S_Fxyz_Py_d = I_ERI_Dxz_S_Gx2yz_S_d+CDY*I_ERI_Dxz_S_Fxyz_S_d;
  Double I_ERI_D2y_S_Fxyz_Py_d = I_ERI_D2y_S_Gx2yz_S_d+CDY*I_ERI_D2y_S_Fxyz_S_d;
  Double I_ERI_Dyz_S_Fxyz_Py_d = I_ERI_Dyz_S_Gx2yz_S_d+CDY*I_ERI_Dyz_S_Fxyz_S_d;
  Double I_ERI_D2z_S_Fxyz_Py_d = I_ERI_D2z_S_Gx2yz_S_d+CDY*I_ERI_D2z_S_Fxyz_S_d;
  Double I_ERI_D2x_S_Fx2z_Py_d = I_ERI_D2x_S_Gxy2z_S_d+CDY*I_ERI_D2x_S_Fx2z_S_d;
  Double I_ERI_Dxy_S_Fx2z_Py_d = I_ERI_Dxy_S_Gxy2z_S_d+CDY*I_ERI_Dxy_S_Fx2z_S_d;
  Double I_ERI_Dxz_S_Fx2z_Py_d = I_ERI_Dxz_S_Gxy2z_S_d+CDY*I_ERI_Dxz_S_Fx2z_S_d;
  Double I_ERI_D2y_S_Fx2z_Py_d = I_ERI_D2y_S_Gxy2z_S_d+CDY*I_ERI_D2y_S_Fx2z_S_d;
  Double I_ERI_Dyz_S_Fx2z_Py_d = I_ERI_Dyz_S_Gxy2z_S_d+CDY*I_ERI_Dyz_S_Fx2z_S_d;
  Double I_ERI_D2z_S_Fx2z_Py_d = I_ERI_D2z_S_Gxy2z_S_d+CDY*I_ERI_D2z_S_Fx2z_S_d;
  Double I_ERI_D2x_S_F3y_Py_d = I_ERI_D2x_S_G4y_S_d+CDY*I_ERI_D2x_S_F3y_S_d;
  Double I_ERI_Dxy_S_F3y_Py_d = I_ERI_Dxy_S_G4y_S_d+CDY*I_ERI_Dxy_S_F3y_S_d;
  Double I_ERI_Dxz_S_F3y_Py_d = I_ERI_Dxz_S_G4y_S_d+CDY*I_ERI_Dxz_S_F3y_S_d;
  Double I_ERI_D2y_S_F3y_Py_d = I_ERI_D2y_S_G4y_S_d+CDY*I_ERI_D2y_S_F3y_S_d;
  Double I_ERI_Dyz_S_F3y_Py_d = I_ERI_Dyz_S_G4y_S_d+CDY*I_ERI_Dyz_S_F3y_S_d;
  Double I_ERI_D2z_S_F3y_Py_d = I_ERI_D2z_S_G4y_S_d+CDY*I_ERI_D2z_S_F3y_S_d;
  Double I_ERI_D2x_S_F2yz_Py_d = I_ERI_D2x_S_G3yz_S_d+CDY*I_ERI_D2x_S_F2yz_S_d;
  Double I_ERI_Dxy_S_F2yz_Py_d = I_ERI_Dxy_S_G3yz_S_d+CDY*I_ERI_Dxy_S_F2yz_S_d;
  Double I_ERI_Dxz_S_F2yz_Py_d = I_ERI_Dxz_S_G3yz_S_d+CDY*I_ERI_Dxz_S_F2yz_S_d;
  Double I_ERI_D2y_S_F2yz_Py_d = I_ERI_D2y_S_G3yz_S_d+CDY*I_ERI_D2y_S_F2yz_S_d;
  Double I_ERI_Dyz_S_F2yz_Py_d = I_ERI_Dyz_S_G3yz_S_d+CDY*I_ERI_Dyz_S_F2yz_S_d;
  Double I_ERI_D2z_S_F2yz_Py_d = I_ERI_D2z_S_G3yz_S_d+CDY*I_ERI_D2z_S_F2yz_S_d;
  Double I_ERI_D2x_S_Fy2z_Py_d = I_ERI_D2x_S_G2y2z_S_d+CDY*I_ERI_D2x_S_Fy2z_S_d;
  Double I_ERI_Dxy_S_Fy2z_Py_d = I_ERI_Dxy_S_G2y2z_S_d+CDY*I_ERI_Dxy_S_Fy2z_S_d;
  Double I_ERI_Dxz_S_Fy2z_Py_d = I_ERI_Dxz_S_G2y2z_S_d+CDY*I_ERI_Dxz_S_Fy2z_S_d;
  Double I_ERI_D2y_S_Fy2z_Py_d = I_ERI_D2y_S_G2y2z_S_d+CDY*I_ERI_D2y_S_Fy2z_S_d;
  Double I_ERI_Dyz_S_Fy2z_Py_d = I_ERI_Dyz_S_G2y2z_S_d+CDY*I_ERI_Dyz_S_Fy2z_S_d;
  Double I_ERI_D2z_S_Fy2z_Py_d = I_ERI_D2z_S_G2y2z_S_d+CDY*I_ERI_D2z_S_Fy2z_S_d;
  Double I_ERI_D2x_S_F3z_Py_d = I_ERI_D2x_S_Gy3z_S_d+CDY*I_ERI_D2x_S_F3z_S_d;
  Double I_ERI_Dxy_S_F3z_Py_d = I_ERI_Dxy_S_Gy3z_S_d+CDY*I_ERI_Dxy_S_F3z_S_d;
  Double I_ERI_Dxz_S_F3z_Py_d = I_ERI_Dxz_S_Gy3z_S_d+CDY*I_ERI_Dxz_S_F3z_S_d;
  Double I_ERI_D2y_S_F3z_Py_d = I_ERI_D2y_S_Gy3z_S_d+CDY*I_ERI_D2y_S_F3z_S_d;
  Double I_ERI_Dyz_S_F3z_Py_d = I_ERI_Dyz_S_Gy3z_S_d+CDY*I_ERI_Dyz_S_F3z_S_d;
  Double I_ERI_D2z_S_F3z_Py_d = I_ERI_D2z_S_Gy3z_S_d+CDY*I_ERI_D2z_S_F3z_S_d;
  Double I_ERI_D2x_S_F3x_Pz_d = I_ERI_D2x_S_G3xz_S_d+CDZ*I_ERI_D2x_S_F3x_S_d;
  Double I_ERI_Dxy_S_F3x_Pz_d = I_ERI_Dxy_S_G3xz_S_d+CDZ*I_ERI_Dxy_S_F3x_S_d;
  Double I_ERI_Dxz_S_F3x_Pz_d = I_ERI_Dxz_S_G3xz_S_d+CDZ*I_ERI_Dxz_S_F3x_S_d;
  Double I_ERI_D2y_S_F3x_Pz_d = I_ERI_D2y_S_G3xz_S_d+CDZ*I_ERI_D2y_S_F3x_S_d;
  Double I_ERI_Dyz_S_F3x_Pz_d = I_ERI_Dyz_S_G3xz_S_d+CDZ*I_ERI_Dyz_S_F3x_S_d;
  Double I_ERI_D2z_S_F3x_Pz_d = I_ERI_D2z_S_G3xz_S_d+CDZ*I_ERI_D2z_S_F3x_S_d;
  Double I_ERI_D2x_S_F2xy_Pz_d = I_ERI_D2x_S_G2xyz_S_d+CDZ*I_ERI_D2x_S_F2xy_S_d;
  Double I_ERI_Dxy_S_F2xy_Pz_d = I_ERI_Dxy_S_G2xyz_S_d+CDZ*I_ERI_Dxy_S_F2xy_S_d;
  Double I_ERI_Dxz_S_F2xy_Pz_d = I_ERI_Dxz_S_G2xyz_S_d+CDZ*I_ERI_Dxz_S_F2xy_S_d;
  Double I_ERI_D2y_S_F2xy_Pz_d = I_ERI_D2y_S_G2xyz_S_d+CDZ*I_ERI_D2y_S_F2xy_S_d;
  Double I_ERI_Dyz_S_F2xy_Pz_d = I_ERI_Dyz_S_G2xyz_S_d+CDZ*I_ERI_Dyz_S_F2xy_S_d;
  Double I_ERI_D2z_S_F2xy_Pz_d = I_ERI_D2z_S_G2xyz_S_d+CDZ*I_ERI_D2z_S_F2xy_S_d;
  Double I_ERI_D2x_S_F2xz_Pz_d = I_ERI_D2x_S_G2x2z_S_d+CDZ*I_ERI_D2x_S_F2xz_S_d;
  Double I_ERI_Dxy_S_F2xz_Pz_d = I_ERI_Dxy_S_G2x2z_S_d+CDZ*I_ERI_Dxy_S_F2xz_S_d;
  Double I_ERI_Dxz_S_F2xz_Pz_d = I_ERI_Dxz_S_G2x2z_S_d+CDZ*I_ERI_Dxz_S_F2xz_S_d;
  Double I_ERI_D2y_S_F2xz_Pz_d = I_ERI_D2y_S_G2x2z_S_d+CDZ*I_ERI_D2y_S_F2xz_S_d;
  Double I_ERI_Dyz_S_F2xz_Pz_d = I_ERI_Dyz_S_G2x2z_S_d+CDZ*I_ERI_Dyz_S_F2xz_S_d;
  Double I_ERI_D2z_S_F2xz_Pz_d = I_ERI_D2z_S_G2x2z_S_d+CDZ*I_ERI_D2z_S_F2xz_S_d;
  Double I_ERI_D2x_S_Fx2y_Pz_d = I_ERI_D2x_S_Gx2yz_S_d+CDZ*I_ERI_D2x_S_Fx2y_S_d;
  Double I_ERI_Dxy_S_Fx2y_Pz_d = I_ERI_Dxy_S_Gx2yz_S_d+CDZ*I_ERI_Dxy_S_Fx2y_S_d;
  Double I_ERI_Dxz_S_Fx2y_Pz_d = I_ERI_Dxz_S_Gx2yz_S_d+CDZ*I_ERI_Dxz_S_Fx2y_S_d;
  Double I_ERI_D2y_S_Fx2y_Pz_d = I_ERI_D2y_S_Gx2yz_S_d+CDZ*I_ERI_D2y_S_Fx2y_S_d;
  Double I_ERI_Dyz_S_Fx2y_Pz_d = I_ERI_Dyz_S_Gx2yz_S_d+CDZ*I_ERI_Dyz_S_Fx2y_S_d;
  Double I_ERI_D2z_S_Fx2y_Pz_d = I_ERI_D2z_S_Gx2yz_S_d+CDZ*I_ERI_D2z_S_Fx2y_S_d;
  Double I_ERI_D2x_S_Fxyz_Pz_d = I_ERI_D2x_S_Gxy2z_S_d+CDZ*I_ERI_D2x_S_Fxyz_S_d;
  Double I_ERI_Dxy_S_Fxyz_Pz_d = I_ERI_Dxy_S_Gxy2z_S_d+CDZ*I_ERI_Dxy_S_Fxyz_S_d;
  Double I_ERI_Dxz_S_Fxyz_Pz_d = I_ERI_Dxz_S_Gxy2z_S_d+CDZ*I_ERI_Dxz_S_Fxyz_S_d;
  Double I_ERI_D2y_S_Fxyz_Pz_d = I_ERI_D2y_S_Gxy2z_S_d+CDZ*I_ERI_D2y_S_Fxyz_S_d;
  Double I_ERI_Dyz_S_Fxyz_Pz_d = I_ERI_Dyz_S_Gxy2z_S_d+CDZ*I_ERI_Dyz_S_Fxyz_S_d;
  Double I_ERI_D2z_S_Fxyz_Pz_d = I_ERI_D2z_S_Gxy2z_S_d+CDZ*I_ERI_D2z_S_Fxyz_S_d;
  Double I_ERI_D2x_S_Fx2z_Pz_d = I_ERI_D2x_S_Gx3z_S_d+CDZ*I_ERI_D2x_S_Fx2z_S_d;
  Double I_ERI_Dxy_S_Fx2z_Pz_d = I_ERI_Dxy_S_Gx3z_S_d+CDZ*I_ERI_Dxy_S_Fx2z_S_d;
  Double I_ERI_Dxz_S_Fx2z_Pz_d = I_ERI_Dxz_S_Gx3z_S_d+CDZ*I_ERI_Dxz_S_Fx2z_S_d;
  Double I_ERI_D2y_S_Fx2z_Pz_d = I_ERI_D2y_S_Gx3z_S_d+CDZ*I_ERI_D2y_S_Fx2z_S_d;
  Double I_ERI_Dyz_S_Fx2z_Pz_d = I_ERI_Dyz_S_Gx3z_S_d+CDZ*I_ERI_Dyz_S_Fx2z_S_d;
  Double I_ERI_D2z_S_Fx2z_Pz_d = I_ERI_D2z_S_Gx3z_S_d+CDZ*I_ERI_D2z_S_Fx2z_S_d;
  Double I_ERI_D2x_S_F3y_Pz_d = I_ERI_D2x_S_G3yz_S_d+CDZ*I_ERI_D2x_S_F3y_S_d;
  Double I_ERI_Dxy_S_F3y_Pz_d = I_ERI_Dxy_S_G3yz_S_d+CDZ*I_ERI_Dxy_S_F3y_S_d;
  Double I_ERI_Dxz_S_F3y_Pz_d = I_ERI_Dxz_S_G3yz_S_d+CDZ*I_ERI_Dxz_S_F3y_S_d;
  Double I_ERI_D2y_S_F3y_Pz_d = I_ERI_D2y_S_G3yz_S_d+CDZ*I_ERI_D2y_S_F3y_S_d;
  Double I_ERI_Dyz_S_F3y_Pz_d = I_ERI_Dyz_S_G3yz_S_d+CDZ*I_ERI_Dyz_S_F3y_S_d;
  Double I_ERI_D2z_S_F3y_Pz_d = I_ERI_D2z_S_G3yz_S_d+CDZ*I_ERI_D2z_S_F3y_S_d;
  Double I_ERI_D2x_S_F2yz_Pz_d = I_ERI_D2x_S_G2y2z_S_d+CDZ*I_ERI_D2x_S_F2yz_S_d;
  Double I_ERI_Dxy_S_F2yz_Pz_d = I_ERI_Dxy_S_G2y2z_S_d+CDZ*I_ERI_Dxy_S_F2yz_S_d;
  Double I_ERI_Dxz_S_F2yz_Pz_d = I_ERI_Dxz_S_G2y2z_S_d+CDZ*I_ERI_Dxz_S_F2yz_S_d;
  Double I_ERI_D2y_S_F2yz_Pz_d = I_ERI_D2y_S_G2y2z_S_d+CDZ*I_ERI_D2y_S_F2yz_S_d;
  Double I_ERI_Dyz_S_F2yz_Pz_d = I_ERI_Dyz_S_G2y2z_S_d+CDZ*I_ERI_Dyz_S_F2yz_S_d;
  Double I_ERI_D2z_S_F2yz_Pz_d = I_ERI_D2z_S_G2y2z_S_d+CDZ*I_ERI_D2z_S_F2yz_S_d;
  Double I_ERI_D2x_S_Fy2z_Pz_d = I_ERI_D2x_S_Gy3z_S_d+CDZ*I_ERI_D2x_S_Fy2z_S_d;
  Double I_ERI_Dxy_S_Fy2z_Pz_d = I_ERI_Dxy_S_Gy3z_S_d+CDZ*I_ERI_Dxy_S_Fy2z_S_d;
  Double I_ERI_Dxz_S_Fy2z_Pz_d = I_ERI_Dxz_S_Gy3z_S_d+CDZ*I_ERI_Dxz_S_Fy2z_S_d;
  Double I_ERI_D2y_S_Fy2z_Pz_d = I_ERI_D2y_S_Gy3z_S_d+CDZ*I_ERI_D2y_S_Fy2z_S_d;
  Double I_ERI_Dyz_S_Fy2z_Pz_d = I_ERI_Dyz_S_Gy3z_S_d+CDZ*I_ERI_Dyz_S_Fy2z_S_d;
  Double I_ERI_D2z_S_Fy2z_Pz_d = I_ERI_D2z_S_Gy3z_S_d+CDZ*I_ERI_D2z_S_Fy2z_S_d;
  Double I_ERI_D2x_S_F3z_Pz_d = I_ERI_D2x_S_G4z_S_d+CDZ*I_ERI_D2x_S_F3z_S_d;
  Double I_ERI_Dxy_S_F3z_Pz_d = I_ERI_Dxy_S_G4z_S_d+CDZ*I_ERI_Dxy_S_F3z_S_d;
  Double I_ERI_Dxz_S_F3z_Pz_d = I_ERI_Dxz_S_G4z_S_d+CDZ*I_ERI_Dxz_S_F3z_S_d;
  Double I_ERI_D2y_S_F3z_Pz_d = I_ERI_D2y_S_G4z_S_d+CDZ*I_ERI_D2y_S_F3z_S_d;
  Double I_ERI_Dyz_S_F3z_Pz_d = I_ERI_Dyz_S_G4z_S_d+CDZ*I_ERI_Dyz_S_F3z_S_d;
  Double I_ERI_D2z_S_F3z_Pz_d = I_ERI_D2z_S_G4z_S_d+CDZ*I_ERI_D2z_S_F3z_S_d;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_F_P_d
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S_d
   * RHS shell quartet name: SQ_ERI_F_S_F_S_d
   ************************************************************/
  Double I_ERI_F3x_S_F3x_Px_d = I_ERI_F3x_S_G4x_S_d+CDX*I_ERI_F3x_S_F3x_S_d;
  Double I_ERI_F2xy_S_F3x_Px_d = I_ERI_F2xy_S_G4x_S_d+CDX*I_ERI_F2xy_S_F3x_S_d;
  Double I_ERI_F2xz_S_F3x_Px_d = I_ERI_F2xz_S_G4x_S_d+CDX*I_ERI_F2xz_S_F3x_S_d;
  Double I_ERI_Fx2y_S_F3x_Px_d = I_ERI_Fx2y_S_G4x_S_d+CDX*I_ERI_Fx2y_S_F3x_S_d;
  Double I_ERI_Fxyz_S_F3x_Px_d = I_ERI_Fxyz_S_G4x_S_d+CDX*I_ERI_Fxyz_S_F3x_S_d;
  Double I_ERI_Fx2z_S_F3x_Px_d = I_ERI_Fx2z_S_G4x_S_d+CDX*I_ERI_Fx2z_S_F3x_S_d;
  Double I_ERI_F3y_S_F3x_Px_d = I_ERI_F3y_S_G4x_S_d+CDX*I_ERI_F3y_S_F3x_S_d;
  Double I_ERI_F2yz_S_F3x_Px_d = I_ERI_F2yz_S_G4x_S_d+CDX*I_ERI_F2yz_S_F3x_S_d;
  Double I_ERI_Fy2z_S_F3x_Px_d = I_ERI_Fy2z_S_G4x_S_d+CDX*I_ERI_Fy2z_S_F3x_S_d;
  Double I_ERI_F3z_S_F3x_Px_d = I_ERI_F3z_S_G4x_S_d+CDX*I_ERI_F3z_S_F3x_S_d;
  Double I_ERI_F3x_S_F2xy_Px_d = I_ERI_F3x_S_G3xy_S_d+CDX*I_ERI_F3x_S_F2xy_S_d;
  Double I_ERI_F2xy_S_F2xy_Px_d = I_ERI_F2xy_S_G3xy_S_d+CDX*I_ERI_F2xy_S_F2xy_S_d;
  Double I_ERI_F2xz_S_F2xy_Px_d = I_ERI_F2xz_S_G3xy_S_d+CDX*I_ERI_F2xz_S_F2xy_S_d;
  Double I_ERI_Fx2y_S_F2xy_Px_d = I_ERI_Fx2y_S_G3xy_S_d+CDX*I_ERI_Fx2y_S_F2xy_S_d;
  Double I_ERI_Fxyz_S_F2xy_Px_d = I_ERI_Fxyz_S_G3xy_S_d+CDX*I_ERI_Fxyz_S_F2xy_S_d;
  Double I_ERI_Fx2z_S_F2xy_Px_d = I_ERI_Fx2z_S_G3xy_S_d+CDX*I_ERI_Fx2z_S_F2xy_S_d;
  Double I_ERI_F3y_S_F2xy_Px_d = I_ERI_F3y_S_G3xy_S_d+CDX*I_ERI_F3y_S_F2xy_S_d;
  Double I_ERI_F2yz_S_F2xy_Px_d = I_ERI_F2yz_S_G3xy_S_d+CDX*I_ERI_F2yz_S_F2xy_S_d;
  Double I_ERI_Fy2z_S_F2xy_Px_d = I_ERI_Fy2z_S_G3xy_S_d+CDX*I_ERI_Fy2z_S_F2xy_S_d;
  Double I_ERI_F3z_S_F2xy_Px_d = I_ERI_F3z_S_G3xy_S_d+CDX*I_ERI_F3z_S_F2xy_S_d;
  Double I_ERI_F3x_S_F2xz_Px_d = I_ERI_F3x_S_G3xz_S_d+CDX*I_ERI_F3x_S_F2xz_S_d;
  Double I_ERI_F2xy_S_F2xz_Px_d = I_ERI_F2xy_S_G3xz_S_d+CDX*I_ERI_F2xy_S_F2xz_S_d;
  Double I_ERI_F2xz_S_F2xz_Px_d = I_ERI_F2xz_S_G3xz_S_d+CDX*I_ERI_F2xz_S_F2xz_S_d;
  Double I_ERI_Fx2y_S_F2xz_Px_d = I_ERI_Fx2y_S_G3xz_S_d+CDX*I_ERI_Fx2y_S_F2xz_S_d;
  Double I_ERI_Fxyz_S_F2xz_Px_d = I_ERI_Fxyz_S_G3xz_S_d+CDX*I_ERI_Fxyz_S_F2xz_S_d;
  Double I_ERI_Fx2z_S_F2xz_Px_d = I_ERI_Fx2z_S_G3xz_S_d+CDX*I_ERI_Fx2z_S_F2xz_S_d;
  Double I_ERI_F3y_S_F2xz_Px_d = I_ERI_F3y_S_G3xz_S_d+CDX*I_ERI_F3y_S_F2xz_S_d;
  Double I_ERI_F2yz_S_F2xz_Px_d = I_ERI_F2yz_S_G3xz_S_d+CDX*I_ERI_F2yz_S_F2xz_S_d;
  Double I_ERI_Fy2z_S_F2xz_Px_d = I_ERI_Fy2z_S_G3xz_S_d+CDX*I_ERI_Fy2z_S_F2xz_S_d;
  Double I_ERI_F3z_S_F2xz_Px_d = I_ERI_F3z_S_G3xz_S_d+CDX*I_ERI_F3z_S_F2xz_S_d;
  Double I_ERI_F3x_S_Fx2y_Px_d = I_ERI_F3x_S_G2x2y_S_d+CDX*I_ERI_F3x_S_Fx2y_S_d;
  Double I_ERI_F2xy_S_Fx2y_Px_d = I_ERI_F2xy_S_G2x2y_S_d+CDX*I_ERI_F2xy_S_Fx2y_S_d;
  Double I_ERI_F2xz_S_Fx2y_Px_d = I_ERI_F2xz_S_G2x2y_S_d+CDX*I_ERI_F2xz_S_Fx2y_S_d;
  Double I_ERI_Fx2y_S_Fx2y_Px_d = I_ERI_Fx2y_S_G2x2y_S_d+CDX*I_ERI_Fx2y_S_Fx2y_S_d;
  Double I_ERI_Fxyz_S_Fx2y_Px_d = I_ERI_Fxyz_S_G2x2y_S_d+CDX*I_ERI_Fxyz_S_Fx2y_S_d;
  Double I_ERI_Fx2z_S_Fx2y_Px_d = I_ERI_Fx2z_S_G2x2y_S_d+CDX*I_ERI_Fx2z_S_Fx2y_S_d;
  Double I_ERI_F3y_S_Fx2y_Px_d = I_ERI_F3y_S_G2x2y_S_d+CDX*I_ERI_F3y_S_Fx2y_S_d;
  Double I_ERI_F2yz_S_Fx2y_Px_d = I_ERI_F2yz_S_G2x2y_S_d+CDX*I_ERI_F2yz_S_Fx2y_S_d;
  Double I_ERI_Fy2z_S_Fx2y_Px_d = I_ERI_Fy2z_S_G2x2y_S_d+CDX*I_ERI_Fy2z_S_Fx2y_S_d;
  Double I_ERI_F3z_S_Fx2y_Px_d = I_ERI_F3z_S_G2x2y_S_d+CDX*I_ERI_F3z_S_Fx2y_S_d;
  Double I_ERI_F3x_S_Fxyz_Px_d = I_ERI_F3x_S_G2xyz_S_d+CDX*I_ERI_F3x_S_Fxyz_S_d;
  Double I_ERI_F2xy_S_Fxyz_Px_d = I_ERI_F2xy_S_G2xyz_S_d+CDX*I_ERI_F2xy_S_Fxyz_S_d;
  Double I_ERI_F2xz_S_Fxyz_Px_d = I_ERI_F2xz_S_G2xyz_S_d+CDX*I_ERI_F2xz_S_Fxyz_S_d;
  Double I_ERI_Fx2y_S_Fxyz_Px_d = I_ERI_Fx2y_S_G2xyz_S_d+CDX*I_ERI_Fx2y_S_Fxyz_S_d;
  Double I_ERI_Fxyz_S_Fxyz_Px_d = I_ERI_Fxyz_S_G2xyz_S_d+CDX*I_ERI_Fxyz_S_Fxyz_S_d;
  Double I_ERI_Fx2z_S_Fxyz_Px_d = I_ERI_Fx2z_S_G2xyz_S_d+CDX*I_ERI_Fx2z_S_Fxyz_S_d;
  Double I_ERI_F3y_S_Fxyz_Px_d = I_ERI_F3y_S_G2xyz_S_d+CDX*I_ERI_F3y_S_Fxyz_S_d;
  Double I_ERI_F2yz_S_Fxyz_Px_d = I_ERI_F2yz_S_G2xyz_S_d+CDX*I_ERI_F2yz_S_Fxyz_S_d;
  Double I_ERI_Fy2z_S_Fxyz_Px_d = I_ERI_Fy2z_S_G2xyz_S_d+CDX*I_ERI_Fy2z_S_Fxyz_S_d;
  Double I_ERI_F3z_S_Fxyz_Px_d = I_ERI_F3z_S_G2xyz_S_d+CDX*I_ERI_F3z_S_Fxyz_S_d;
  Double I_ERI_F3x_S_Fx2z_Px_d = I_ERI_F3x_S_G2x2z_S_d+CDX*I_ERI_F3x_S_Fx2z_S_d;
  Double I_ERI_F2xy_S_Fx2z_Px_d = I_ERI_F2xy_S_G2x2z_S_d+CDX*I_ERI_F2xy_S_Fx2z_S_d;
  Double I_ERI_F2xz_S_Fx2z_Px_d = I_ERI_F2xz_S_G2x2z_S_d+CDX*I_ERI_F2xz_S_Fx2z_S_d;
  Double I_ERI_Fx2y_S_Fx2z_Px_d = I_ERI_Fx2y_S_G2x2z_S_d+CDX*I_ERI_Fx2y_S_Fx2z_S_d;
  Double I_ERI_Fxyz_S_Fx2z_Px_d = I_ERI_Fxyz_S_G2x2z_S_d+CDX*I_ERI_Fxyz_S_Fx2z_S_d;
  Double I_ERI_Fx2z_S_Fx2z_Px_d = I_ERI_Fx2z_S_G2x2z_S_d+CDX*I_ERI_Fx2z_S_Fx2z_S_d;
  Double I_ERI_F3y_S_Fx2z_Px_d = I_ERI_F3y_S_G2x2z_S_d+CDX*I_ERI_F3y_S_Fx2z_S_d;
  Double I_ERI_F2yz_S_Fx2z_Px_d = I_ERI_F2yz_S_G2x2z_S_d+CDX*I_ERI_F2yz_S_Fx2z_S_d;
  Double I_ERI_Fy2z_S_Fx2z_Px_d = I_ERI_Fy2z_S_G2x2z_S_d+CDX*I_ERI_Fy2z_S_Fx2z_S_d;
  Double I_ERI_F3z_S_Fx2z_Px_d = I_ERI_F3z_S_G2x2z_S_d+CDX*I_ERI_F3z_S_Fx2z_S_d;
  Double I_ERI_F3x_S_F3y_Px_d = I_ERI_F3x_S_Gx3y_S_d+CDX*I_ERI_F3x_S_F3y_S_d;
  Double I_ERI_F2xy_S_F3y_Px_d = I_ERI_F2xy_S_Gx3y_S_d+CDX*I_ERI_F2xy_S_F3y_S_d;
  Double I_ERI_F2xz_S_F3y_Px_d = I_ERI_F2xz_S_Gx3y_S_d+CDX*I_ERI_F2xz_S_F3y_S_d;
  Double I_ERI_Fx2y_S_F3y_Px_d = I_ERI_Fx2y_S_Gx3y_S_d+CDX*I_ERI_Fx2y_S_F3y_S_d;
  Double I_ERI_Fxyz_S_F3y_Px_d = I_ERI_Fxyz_S_Gx3y_S_d+CDX*I_ERI_Fxyz_S_F3y_S_d;
  Double I_ERI_Fx2z_S_F3y_Px_d = I_ERI_Fx2z_S_Gx3y_S_d+CDX*I_ERI_Fx2z_S_F3y_S_d;
  Double I_ERI_F3y_S_F3y_Px_d = I_ERI_F3y_S_Gx3y_S_d+CDX*I_ERI_F3y_S_F3y_S_d;
  Double I_ERI_F2yz_S_F3y_Px_d = I_ERI_F2yz_S_Gx3y_S_d+CDX*I_ERI_F2yz_S_F3y_S_d;
  Double I_ERI_Fy2z_S_F3y_Px_d = I_ERI_Fy2z_S_Gx3y_S_d+CDX*I_ERI_Fy2z_S_F3y_S_d;
  Double I_ERI_F3z_S_F3y_Px_d = I_ERI_F3z_S_Gx3y_S_d+CDX*I_ERI_F3z_S_F3y_S_d;
  Double I_ERI_F3x_S_F2yz_Px_d = I_ERI_F3x_S_Gx2yz_S_d+CDX*I_ERI_F3x_S_F2yz_S_d;
  Double I_ERI_F2xy_S_F2yz_Px_d = I_ERI_F2xy_S_Gx2yz_S_d+CDX*I_ERI_F2xy_S_F2yz_S_d;
  Double I_ERI_F2xz_S_F2yz_Px_d = I_ERI_F2xz_S_Gx2yz_S_d+CDX*I_ERI_F2xz_S_F2yz_S_d;
  Double I_ERI_Fx2y_S_F2yz_Px_d = I_ERI_Fx2y_S_Gx2yz_S_d+CDX*I_ERI_Fx2y_S_F2yz_S_d;
  Double I_ERI_Fxyz_S_F2yz_Px_d = I_ERI_Fxyz_S_Gx2yz_S_d+CDX*I_ERI_Fxyz_S_F2yz_S_d;
  Double I_ERI_Fx2z_S_F2yz_Px_d = I_ERI_Fx2z_S_Gx2yz_S_d+CDX*I_ERI_Fx2z_S_F2yz_S_d;
  Double I_ERI_F3y_S_F2yz_Px_d = I_ERI_F3y_S_Gx2yz_S_d+CDX*I_ERI_F3y_S_F2yz_S_d;
  Double I_ERI_F2yz_S_F2yz_Px_d = I_ERI_F2yz_S_Gx2yz_S_d+CDX*I_ERI_F2yz_S_F2yz_S_d;
  Double I_ERI_Fy2z_S_F2yz_Px_d = I_ERI_Fy2z_S_Gx2yz_S_d+CDX*I_ERI_Fy2z_S_F2yz_S_d;
  Double I_ERI_F3z_S_F2yz_Px_d = I_ERI_F3z_S_Gx2yz_S_d+CDX*I_ERI_F3z_S_F2yz_S_d;
  Double I_ERI_F3x_S_Fy2z_Px_d = I_ERI_F3x_S_Gxy2z_S_d+CDX*I_ERI_F3x_S_Fy2z_S_d;
  Double I_ERI_F2xy_S_Fy2z_Px_d = I_ERI_F2xy_S_Gxy2z_S_d+CDX*I_ERI_F2xy_S_Fy2z_S_d;
  Double I_ERI_F2xz_S_Fy2z_Px_d = I_ERI_F2xz_S_Gxy2z_S_d+CDX*I_ERI_F2xz_S_Fy2z_S_d;
  Double I_ERI_Fx2y_S_Fy2z_Px_d = I_ERI_Fx2y_S_Gxy2z_S_d+CDX*I_ERI_Fx2y_S_Fy2z_S_d;
  Double I_ERI_Fxyz_S_Fy2z_Px_d = I_ERI_Fxyz_S_Gxy2z_S_d+CDX*I_ERI_Fxyz_S_Fy2z_S_d;
  Double I_ERI_Fx2z_S_Fy2z_Px_d = I_ERI_Fx2z_S_Gxy2z_S_d+CDX*I_ERI_Fx2z_S_Fy2z_S_d;
  Double I_ERI_F3y_S_Fy2z_Px_d = I_ERI_F3y_S_Gxy2z_S_d+CDX*I_ERI_F3y_S_Fy2z_S_d;
  Double I_ERI_F2yz_S_Fy2z_Px_d = I_ERI_F2yz_S_Gxy2z_S_d+CDX*I_ERI_F2yz_S_Fy2z_S_d;
  Double I_ERI_Fy2z_S_Fy2z_Px_d = I_ERI_Fy2z_S_Gxy2z_S_d+CDX*I_ERI_Fy2z_S_Fy2z_S_d;
  Double I_ERI_F3z_S_Fy2z_Px_d = I_ERI_F3z_S_Gxy2z_S_d+CDX*I_ERI_F3z_S_Fy2z_S_d;
  Double I_ERI_F3x_S_F3z_Px_d = I_ERI_F3x_S_Gx3z_S_d+CDX*I_ERI_F3x_S_F3z_S_d;
  Double I_ERI_F2xy_S_F3z_Px_d = I_ERI_F2xy_S_Gx3z_S_d+CDX*I_ERI_F2xy_S_F3z_S_d;
  Double I_ERI_F2xz_S_F3z_Px_d = I_ERI_F2xz_S_Gx3z_S_d+CDX*I_ERI_F2xz_S_F3z_S_d;
  Double I_ERI_Fx2y_S_F3z_Px_d = I_ERI_Fx2y_S_Gx3z_S_d+CDX*I_ERI_Fx2y_S_F3z_S_d;
  Double I_ERI_Fxyz_S_F3z_Px_d = I_ERI_Fxyz_S_Gx3z_S_d+CDX*I_ERI_Fxyz_S_F3z_S_d;
  Double I_ERI_Fx2z_S_F3z_Px_d = I_ERI_Fx2z_S_Gx3z_S_d+CDX*I_ERI_Fx2z_S_F3z_S_d;
  Double I_ERI_F3y_S_F3z_Px_d = I_ERI_F3y_S_Gx3z_S_d+CDX*I_ERI_F3y_S_F3z_S_d;
  Double I_ERI_F2yz_S_F3z_Px_d = I_ERI_F2yz_S_Gx3z_S_d+CDX*I_ERI_F2yz_S_F3z_S_d;
  Double I_ERI_Fy2z_S_F3z_Px_d = I_ERI_Fy2z_S_Gx3z_S_d+CDX*I_ERI_Fy2z_S_F3z_S_d;
  Double I_ERI_F3z_S_F3z_Px_d = I_ERI_F3z_S_Gx3z_S_d+CDX*I_ERI_F3z_S_F3z_S_d;
  Double I_ERI_F3x_S_F3x_Py_d = I_ERI_F3x_S_G3xy_S_d+CDY*I_ERI_F3x_S_F3x_S_d;
  Double I_ERI_F2xy_S_F3x_Py_d = I_ERI_F2xy_S_G3xy_S_d+CDY*I_ERI_F2xy_S_F3x_S_d;
  Double I_ERI_F2xz_S_F3x_Py_d = I_ERI_F2xz_S_G3xy_S_d+CDY*I_ERI_F2xz_S_F3x_S_d;
  Double I_ERI_Fx2y_S_F3x_Py_d = I_ERI_Fx2y_S_G3xy_S_d+CDY*I_ERI_Fx2y_S_F3x_S_d;
  Double I_ERI_Fxyz_S_F3x_Py_d = I_ERI_Fxyz_S_G3xy_S_d+CDY*I_ERI_Fxyz_S_F3x_S_d;
  Double I_ERI_Fx2z_S_F3x_Py_d = I_ERI_Fx2z_S_G3xy_S_d+CDY*I_ERI_Fx2z_S_F3x_S_d;
  Double I_ERI_F3y_S_F3x_Py_d = I_ERI_F3y_S_G3xy_S_d+CDY*I_ERI_F3y_S_F3x_S_d;
  Double I_ERI_F2yz_S_F3x_Py_d = I_ERI_F2yz_S_G3xy_S_d+CDY*I_ERI_F2yz_S_F3x_S_d;
  Double I_ERI_Fy2z_S_F3x_Py_d = I_ERI_Fy2z_S_G3xy_S_d+CDY*I_ERI_Fy2z_S_F3x_S_d;
  Double I_ERI_F3z_S_F3x_Py_d = I_ERI_F3z_S_G3xy_S_d+CDY*I_ERI_F3z_S_F3x_S_d;
  Double I_ERI_F3x_S_F2xy_Py_d = I_ERI_F3x_S_G2x2y_S_d+CDY*I_ERI_F3x_S_F2xy_S_d;
  Double I_ERI_F2xy_S_F2xy_Py_d = I_ERI_F2xy_S_G2x2y_S_d+CDY*I_ERI_F2xy_S_F2xy_S_d;
  Double I_ERI_F2xz_S_F2xy_Py_d = I_ERI_F2xz_S_G2x2y_S_d+CDY*I_ERI_F2xz_S_F2xy_S_d;
  Double I_ERI_Fx2y_S_F2xy_Py_d = I_ERI_Fx2y_S_G2x2y_S_d+CDY*I_ERI_Fx2y_S_F2xy_S_d;
  Double I_ERI_Fxyz_S_F2xy_Py_d = I_ERI_Fxyz_S_G2x2y_S_d+CDY*I_ERI_Fxyz_S_F2xy_S_d;
  Double I_ERI_Fx2z_S_F2xy_Py_d = I_ERI_Fx2z_S_G2x2y_S_d+CDY*I_ERI_Fx2z_S_F2xy_S_d;
  Double I_ERI_F3y_S_F2xy_Py_d = I_ERI_F3y_S_G2x2y_S_d+CDY*I_ERI_F3y_S_F2xy_S_d;
  Double I_ERI_F2yz_S_F2xy_Py_d = I_ERI_F2yz_S_G2x2y_S_d+CDY*I_ERI_F2yz_S_F2xy_S_d;
  Double I_ERI_Fy2z_S_F2xy_Py_d = I_ERI_Fy2z_S_G2x2y_S_d+CDY*I_ERI_Fy2z_S_F2xy_S_d;
  Double I_ERI_F3z_S_F2xy_Py_d = I_ERI_F3z_S_G2x2y_S_d+CDY*I_ERI_F3z_S_F2xy_S_d;
  Double I_ERI_F3x_S_F2xz_Py_d = I_ERI_F3x_S_G2xyz_S_d+CDY*I_ERI_F3x_S_F2xz_S_d;
  Double I_ERI_F2xy_S_F2xz_Py_d = I_ERI_F2xy_S_G2xyz_S_d+CDY*I_ERI_F2xy_S_F2xz_S_d;
  Double I_ERI_F2xz_S_F2xz_Py_d = I_ERI_F2xz_S_G2xyz_S_d+CDY*I_ERI_F2xz_S_F2xz_S_d;
  Double I_ERI_Fx2y_S_F2xz_Py_d = I_ERI_Fx2y_S_G2xyz_S_d+CDY*I_ERI_Fx2y_S_F2xz_S_d;
  Double I_ERI_Fxyz_S_F2xz_Py_d = I_ERI_Fxyz_S_G2xyz_S_d+CDY*I_ERI_Fxyz_S_F2xz_S_d;
  Double I_ERI_Fx2z_S_F2xz_Py_d = I_ERI_Fx2z_S_G2xyz_S_d+CDY*I_ERI_Fx2z_S_F2xz_S_d;
  Double I_ERI_F3y_S_F2xz_Py_d = I_ERI_F3y_S_G2xyz_S_d+CDY*I_ERI_F3y_S_F2xz_S_d;
  Double I_ERI_F2yz_S_F2xz_Py_d = I_ERI_F2yz_S_G2xyz_S_d+CDY*I_ERI_F2yz_S_F2xz_S_d;
  Double I_ERI_Fy2z_S_F2xz_Py_d = I_ERI_Fy2z_S_G2xyz_S_d+CDY*I_ERI_Fy2z_S_F2xz_S_d;
  Double I_ERI_F3z_S_F2xz_Py_d = I_ERI_F3z_S_G2xyz_S_d+CDY*I_ERI_F3z_S_F2xz_S_d;
  Double I_ERI_F3x_S_Fx2y_Py_d = I_ERI_F3x_S_Gx3y_S_d+CDY*I_ERI_F3x_S_Fx2y_S_d;
  Double I_ERI_F2xy_S_Fx2y_Py_d = I_ERI_F2xy_S_Gx3y_S_d+CDY*I_ERI_F2xy_S_Fx2y_S_d;
  Double I_ERI_F2xz_S_Fx2y_Py_d = I_ERI_F2xz_S_Gx3y_S_d+CDY*I_ERI_F2xz_S_Fx2y_S_d;
  Double I_ERI_Fx2y_S_Fx2y_Py_d = I_ERI_Fx2y_S_Gx3y_S_d+CDY*I_ERI_Fx2y_S_Fx2y_S_d;
  Double I_ERI_Fxyz_S_Fx2y_Py_d = I_ERI_Fxyz_S_Gx3y_S_d+CDY*I_ERI_Fxyz_S_Fx2y_S_d;
  Double I_ERI_Fx2z_S_Fx2y_Py_d = I_ERI_Fx2z_S_Gx3y_S_d+CDY*I_ERI_Fx2z_S_Fx2y_S_d;
  Double I_ERI_F3y_S_Fx2y_Py_d = I_ERI_F3y_S_Gx3y_S_d+CDY*I_ERI_F3y_S_Fx2y_S_d;
  Double I_ERI_F2yz_S_Fx2y_Py_d = I_ERI_F2yz_S_Gx3y_S_d+CDY*I_ERI_F2yz_S_Fx2y_S_d;
  Double I_ERI_Fy2z_S_Fx2y_Py_d = I_ERI_Fy2z_S_Gx3y_S_d+CDY*I_ERI_Fy2z_S_Fx2y_S_d;
  Double I_ERI_F3z_S_Fx2y_Py_d = I_ERI_F3z_S_Gx3y_S_d+CDY*I_ERI_F3z_S_Fx2y_S_d;
  Double I_ERI_F3x_S_Fxyz_Py_d = I_ERI_F3x_S_Gx2yz_S_d+CDY*I_ERI_F3x_S_Fxyz_S_d;
  Double I_ERI_F2xy_S_Fxyz_Py_d = I_ERI_F2xy_S_Gx2yz_S_d+CDY*I_ERI_F2xy_S_Fxyz_S_d;
  Double I_ERI_F2xz_S_Fxyz_Py_d = I_ERI_F2xz_S_Gx2yz_S_d+CDY*I_ERI_F2xz_S_Fxyz_S_d;
  Double I_ERI_Fx2y_S_Fxyz_Py_d = I_ERI_Fx2y_S_Gx2yz_S_d+CDY*I_ERI_Fx2y_S_Fxyz_S_d;
  Double I_ERI_Fxyz_S_Fxyz_Py_d = I_ERI_Fxyz_S_Gx2yz_S_d+CDY*I_ERI_Fxyz_S_Fxyz_S_d;
  Double I_ERI_Fx2z_S_Fxyz_Py_d = I_ERI_Fx2z_S_Gx2yz_S_d+CDY*I_ERI_Fx2z_S_Fxyz_S_d;
  Double I_ERI_F3y_S_Fxyz_Py_d = I_ERI_F3y_S_Gx2yz_S_d+CDY*I_ERI_F3y_S_Fxyz_S_d;
  Double I_ERI_F2yz_S_Fxyz_Py_d = I_ERI_F2yz_S_Gx2yz_S_d+CDY*I_ERI_F2yz_S_Fxyz_S_d;
  Double I_ERI_Fy2z_S_Fxyz_Py_d = I_ERI_Fy2z_S_Gx2yz_S_d+CDY*I_ERI_Fy2z_S_Fxyz_S_d;
  Double I_ERI_F3z_S_Fxyz_Py_d = I_ERI_F3z_S_Gx2yz_S_d+CDY*I_ERI_F3z_S_Fxyz_S_d;
  Double I_ERI_F3x_S_Fx2z_Py_d = I_ERI_F3x_S_Gxy2z_S_d+CDY*I_ERI_F3x_S_Fx2z_S_d;
  Double I_ERI_F2xy_S_Fx2z_Py_d = I_ERI_F2xy_S_Gxy2z_S_d+CDY*I_ERI_F2xy_S_Fx2z_S_d;
  Double I_ERI_F2xz_S_Fx2z_Py_d = I_ERI_F2xz_S_Gxy2z_S_d+CDY*I_ERI_F2xz_S_Fx2z_S_d;
  Double I_ERI_Fx2y_S_Fx2z_Py_d = I_ERI_Fx2y_S_Gxy2z_S_d+CDY*I_ERI_Fx2y_S_Fx2z_S_d;
  Double I_ERI_Fxyz_S_Fx2z_Py_d = I_ERI_Fxyz_S_Gxy2z_S_d+CDY*I_ERI_Fxyz_S_Fx2z_S_d;
  Double I_ERI_Fx2z_S_Fx2z_Py_d = I_ERI_Fx2z_S_Gxy2z_S_d+CDY*I_ERI_Fx2z_S_Fx2z_S_d;
  Double I_ERI_F3y_S_Fx2z_Py_d = I_ERI_F3y_S_Gxy2z_S_d+CDY*I_ERI_F3y_S_Fx2z_S_d;
  Double I_ERI_F2yz_S_Fx2z_Py_d = I_ERI_F2yz_S_Gxy2z_S_d+CDY*I_ERI_F2yz_S_Fx2z_S_d;
  Double I_ERI_Fy2z_S_Fx2z_Py_d = I_ERI_Fy2z_S_Gxy2z_S_d+CDY*I_ERI_Fy2z_S_Fx2z_S_d;
  Double I_ERI_F3z_S_Fx2z_Py_d = I_ERI_F3z_S_Gxy2z_S_d+CDY*I_ERI_F3z_S_Fx2z_S_d;
  Double I_ERI_F3x_S_F3y_Py_d = I_ERI_F3x_S_G4y_S_d+CDY*I_ERI_F3x_S_F3y_S_d;
  Double I_ERI_F2xy_S_F3y_Py_d = I_ERI_F2xy_S_G4y_S_d+CDY*I_ERI_F2xy_S_F3y_S_d;
  Double I_ERI_F2xz_S_F3y_Py_d = I_ERI_F2xz_S_G4y_S_d+CDY*I_ERI_F2xz_S_F3y_S_d;
  Double I_ERI_Fx2y_S_F3y_Py_d = I_ERI_Fx2y_S_G4y_S_d+CDY*I_ERI_Fx2y_S_F3y_S_d;
  Double I_ERI_Fxyz_S_F3y_Py_d = I_ERI_Fxyz_S_G4y_S_d+CDY*I_ERI_Fxyz_S_F3y_S_d;
  Double I_ERI_Fx2z_S_F3y_Py_d = I_ERI_Fx2z_S_G4y_S_d+CDY*I_ERI_Fx2z_S_F3y_S_d;
  Double I_ERI_F3y_S_F3y_Py_d = I_ERI_F3y_S_G4y_S_d+CDY*I_ERI_F3y_S_F3y_S_d;
  Double I_ERI_F2yz_S_F3y_Py_d = I_ERI_F2yz_S_G4y_S_d+CDY*I_ERI_F2yz_S_F3y_S_d;
  Double I_ERI_Fy2z_S_F3y_Py_d = I_ERI_Fy2z_S_G4y_S_d+CDY*I_ERI_Fy2z_S_F3y_S_d;
  Double I_ERI_F3z_S_F3y_Py_d = I_ERI_F3z_S_G4y_S_d+CDY*I_ERI_F3z_S_F3y_S_d;
  Double I_ERI_F3x_S_F2yz_Py_d = I_ERI_F3x_S_G3yz_S_d+CDY*I_ERI_F3x_S_F2yz_S_d;
  Double I_ERI_F2xy_S_F2yz_Py_d = I_ERI_F2xy_S_G3yz_S_d+CDY*I_ERI_F2xy_S_F2yz_S_d;
  Double I_ERI_F2xz_S_F2yz_Py_d = I_ERI_F2xz_S_G3yz_S_d+CDY*I_ERI_F2xz_S_F2yz_S_d;
  Double I_ERI_Fx2y_S_F2yz_Py_d = I_ERI_Fx2y_S_G3yz_S_d+CDY*I_ERI_Fx2y_S_F2yz_S_d;
  Double I_ERI_Fxyz_S_F2yz_Py_d = I_ERI_Fxyz_S_G3yz_S_d+CDY*I_ERI_Fxyz_S_F2yz_S_d;
  Double I_ERI_Fx2z_S_F2yz_Py_d = I_ERI_Fx2z_S_G3yz_S_d+CDY*I_ERI_Fx2z_S_F2yz_S_d;
  Double I_ERI_F3y_S_F2yz_Py_d = I_ERI_F3y_S_G3yz_S_d+CDY*I_ERI_F3y_S_F2yz_S_d;
  Double I_ERI_F2yz_S_F2yz_Py_d = I_ERI_F2yz_S_G3yz_S_d+CDY*I_ERI_F2yz_S_F2yz_S_d;
  Double I_ERI_Fy2z_S_F2yz_Py_d = I_ERI_Fy2z_S_G3yz_S_d+CDY*I_ERI_Fy2z_S_F2yz_S_d;
  Double I_ERI_F3z_S_F2yz_Py_d = I_ERI_F3z_S_G3yz_S_d+CDY*I_ERI_F3z_S_F2yz_S_d;
  Double I_ERI_F3x_S_Fy2z_Py_d = I_ERI_F3x_S_G2y2z_S_d+CDY*I_ERI_F3x_S_Fy2z_S_d;
  Double I_ERI_F2xy_S_Fy2z_Py_d = I_ERI_F2xy_S_G2y2z_S_d+CDY*I_ERI_F2xy_S_Fy2z_S_d;
  Double I_ERI_F2xz_S_Fy2z_Py_d = I_ERI_F2xz_S_G2y2z_S_d+CDY*I_ERI_F2xz_S_Fy2z_S_d;
  Double I_ERI_Fx2y_S_Fy2z_Py_d = I_ERI_Fx2y_S_G2y2z_S_d+CDY*I_ERI_Fx2y_S_Fy2z_S_d;
  Double I_ERI_Fxyz_S_Fy2z_Py_d = I_ERI_Fxyz_S_G2y2z_S_d+CDY*I_ERI_Fxyz_S_Fy2z_S_d;
  Double I_ERI_Fx2z_S_Fy2z_Py_d = I_ERI_Fx2z_S_G2y2z_S_d+CDY*I_ERI_Fx2z_S_Fy2z_S_d;
  Double I_ERI_F3y_S_Fy2z_Py_d = I_ERI_F3y_S_G2y2z_S_d+CDY*I_ERI_F3y_S_Fy2z_S_d;
  Double I_ERI_F2yz_S_Fy2z_Py_d = I_ERI_F2yz_S_G2y2z_S_d+CDY*I_ERI_F2yz_S_Fy2z_S_d;
  Double I_ERI_Fy2z_S_Fy2z_Py_d = I_ERI_Fy2z_S_G2y2z_S_d+CDY*I_ERI_Fy2z_S_Fy2z_S_d;
  Double I_ERI_F3z_S_Fy2z_Py_d = I_ERI_F3z_S_G2y2z_S_d+CDY*I_ERI_F3z_S_Fy2z_S_d;
  Double I_ERI_F3x_S_F3z_Py_d = I_ERI_F3x_S_Gy3z_S_d+CDY*I_ERI_F3x_S_F3z_S_d;
  Double I_ERI_F2xy_S_F3z_Py_d = I_ERI_F2xy_S_Gy3z_S_d+CDY*I_ERI_F2xy_S_F3z_S_d;
  Double I_ERI_F2xz_S_F3z_Py_d = I_ERI_F2xz_S_Gy3z_S_d+CDY*I_ERI_F2xz_S_F3z_S_d;
  Double I_ERI_Fx2y_S_F3z_Py_d = I_ERI_Fx2y_S_Gy3z_S_d+CDY*I_ERI_Fx2y_S_F3z_S_d;
  Double I_ERI_Fxyz_S_F3z_Py_d = I_ERI_Fxyz_S_Gy3z_S_d+CDY*I_ERI_Fxyz_S_F3z_S_d;
  Double I_ERI_Fx2z_S_F3z_Py_d = I_ERI_Fx2z_S_Gy3z_S_d+CDY*I_ERI_Fx2z_S_F3z_S_d;
  Double I_ERI_F3y_S_F3z_Py_d = I_ERI_F3y_S_Gy3z_S_d+CDY*I_ERI_F3y_S_F3z_S_d;
  Double I_ERI_F2yz_S_F3z_Py_d = I_ERI_F2yz_S_Gy3z_S_d+CDY*I_ERI_F2yz_S_F3z_S_d;
  Double I_ERI_Fy2z_S_F3z_Py_d = I_ERI_Fy2z_S_Gy3z_S_d+CDY*I_ERI_Fy2z_S_F3z_S_d;
  Double I_ERI_F3z_S_F3z_Py_d = I_ERI_F3z_S_Gy3z_S_d+CDY*I_ERI_F3z_S_F3z_S_d;
  Double I_ERI_F3x_S_F3x_Pz_d = I_ERI_F3x_S_G3xz_S_d+CDZ*I_ERI_F3x_S_F3x_S_d;
  Double I_ERI_F2xy_S_F3x_Pz_d = I_ERI_F2xy_S_G3xz_S_d+CDZ*I_ERI_F2xy_S_F3x_S_d;
  Double I_ERI_F2xz_S_F3x_Pz_d = I_ERI_F2xz_S_G3xz_S_d+CDZ*I_ERI_F2xz_S_F3x_S_d;
  Double I_ERI_Fx2y_S_F3x_Pz_d = I_ERI_Fx2y_S_G3xz_S_d+CDZ*I_ERI_Fx2y_S_F3x_S_d;
  Double I_ERI_Fxyz_S_F3x_Pz_d = I_ERI_Fxyz_S_G3xz_S_d+CDZ*I_ERI_Fxyz_S_F3x_S_d;
  Double I_ERI_Fx2z_S_F3x_Pz_d = I_ERI_Fx2z_S_G3xz_S_d+CDZ*I_ERI_Fx2z_S_F3x_S_d;
  Double I_ERI_F3y_S_F3x_Pz_d = I_ERI_F3y_S_G3xz_S_d+CDZ*I_ERI_F3y_S_F3x_S_d;
  Double I_ERI_F2yz_S_F3x_Pz_d = I_ERI_F2yz_S_G3xz_S_d+CDZ*I_ERI_F2yz_S_F3x_S_d;
  Double I_ERI_Fy2z_S_F3x_Pz_d = I_ERI_Fy2z_S_G3xz_S_d+CDZ*I_ERI_Fy2z_S_F3x_S_d;
  Double I_ERI_F3z_S_F3x_Pz_d = I_ERI_F3z_S_G3xz_S_d+CDZ*I_ERI_F3z_S_F3x_S_d;
  Double I_ERI_F3x_S_F2xy_Pz_d = I_ERI_F3x_S_G2xyz_S_d+CDZ*I_ERI_F3x_S_F2xy_S_d;
  Double I_ERI_F2xy_S_F2xy_Pz_d = I_ERI_F2xy_S_G2xyz_S_d+CDZ*I_ERI_F2xy_S_F2xy_S_d;
  Double I_ERI_F2xz_S_F2xy_Pz_d = I_ERI_F2xz_S_G2xyz_S_d+CDZ*I_ERI_F2xz_S_F2xy_S_d;
  Double I_ERI_Fx2y_S_F2xy_Pz_d = I_ERI_Fx2y_S_G2xyz_S_d+CDZ*I_ERI_Fx2y_S_F2xy_S_d;
  Double I_ERI_Fxyz_S_F2xy_Pz_d = I_ERI_Fxyz_S_G2xyz_S_d+CDZ*I_ERI_Fxyz_S_F2xy_S_d;
  Double I_ERI_Fx2z_S_F2xy_Pz_d = I_ERI_Fx2z_S_G2xyz_S_d+CDZ*I_ERI_Fx2z_S_F2xy_S_d;
  Double I_ERI_F3y_S_F2xy_Pz_d = I_ERI_F3y_S_G2xyz_S_d+CDZ*I_ERI_F3y_S_F2xy_S_d;
  Double I_ERI_F2yz_S_F2xy_Pz_d = I_ERI_F2yz_S_G2xyz_S_d+CDZ*I_ERI_F2yz_S_F2xy_S_d;
  Double I_ERI_Fy2z_S_F2xy_Pz_d = I_ERI_Fy2z_S_G2xyz_S_d+CDZ*I_ERI_Fy2z_S_F2xy_S_d;
  Double I_ERI_F3z_S_F2xy_Pz_d = I_ERI_F3z_S_G2xyz_S_d+CDZ*I_ERI_F3z_S_F2xy_S_d;
  Double I_ERI_F3x_S_F2xz_Pz_d = I_ERI_F3x_S_G2x2z_S_d+CDZ*I_ERI_F3x_S_F2xz_S_d;
  Double I_ERI_F2xy_S_F2xz_Pz_d = I_ERI_F2xy_S_G2x2z_S_d+CDZ*I_ERI_F2xy_S_F2xz_S_d;
  Double I_ERI_F2xz_S_F2xz_Pz_d = I_ERI_F2xz_S_G2x2z_S_d+CDZ*I_ERI_F2xz_S_F2xz_S_d;
  Double I_ERI_Fx2y_S_F2xz_Pz_d = I_ERI_Fx2y_S_G2x2z_S_d+CDZ*I_ERI_Fx2y_S_F2xz_S_d;
  Double I_ERI_Fxyz_S_F2xz_Pz_d = I_ERI_Fxyz_S_G2x2z_S_d+CDZ*I_ERI_Fxyz_S_F2xz_S_d;
  Double I_ERI_Fx2z_S_F2xz_Pz_d = I_ERI_Fx2z_S_G2x2z_S_d+CDZ*I_ERI_Fx2z_S_F2xz_S_d;
  Double I_ERI_F3y_S_F2xz_Pz_d = I_ERI_F3y_S_G2x2z_S_d+CDZ*I_ERI_F3y_S_F2xz_S_d;
  Double I_ERI_F2yz_S_F2xz_Pz_d = I_ERI_F2yz_S_G2x2z_S_d+CDZ*I_ERI_F2yz_S_F2xz_S_d;
  Double I_ERI_Fy2z_S_F2xz_Pz_d = I_ERI_Fy2z_S_G2x2z_S_d+CDZ*I_ERI_Fy2z_S_F2xz_S_d;
  Double I_ERI_F3z_S_F2xz_Pz_d = I_ERI_F3z_S_G2x2z_S_d+CDZ*I_ERI_F3z_S_F2xz_S_d;
  Double I_ERI_F3x_S_Fx2y_Pz_d = I_ERI_F3x_S_Gx2yz_S_d+CDZ*I_ERI_F3x_S_Fx2y_S_d;
  Double I_ERI_F2xy_S_Fx2y_Pz_d = I_ERI_F2xy_S_Gx2yz_S_d+CDZ*I_ERI_F2xy_S_Fx2y_S_d;
  Double I_ERI_F2xz_S_Fx2y_Pz_d = I_ERI_F2xz_S_Gx2yz_S_d+CDZ*I_ERI_F2xz_S_Fx2y_S_d;
  Double I_ERI_Fx2y_S_Fx2y_Pz_d = I_ERI_Fx2y_S_Gx2yz_S_d+CDZ*I_ERI_Fx2y_S_Fx2y_S_d;
  Double I_ERI_Fxyz_S_Fx2y_Pz_d = I_ERI_Fxyz_S_Gx2yz_S_d+CDZ*I_ERI_Fxyz_S_Fx2y_S_d;
  Double I_ERI_Fx2z_S_Fx2y_Pz_d = I_ERI_Fx2z_S_Gx2yz_S_d+CDZ*I_ERI_Fx2z_S_Fx2y_S_d;
  Double I_ERI_F3y_S_Fx2y_Pz_d = I_ERI_F3y_S_Gx2yz_S_d+CDZ*I_ERI_F3y_S_Fx2y_S_d;
  Double I_ERI_F2yz_S_Fx2y_Pz_d = I_ERI_F2yz_S_Gx2yz_S_d+CDZ*I_ERI_F2yz_S_Fx2y_S_d;
  Double I_ERI_Fy2z_S_Fx2y_Pz_d = I_ERI_Fy2z_S_Gx2yz_S_d+CDZ*I_ERI_Fy2z_S_Fx2y_S_d;
  Double I_ERI_F3z_S_Fx2y_Pz_d = I_ERI_F3z_S_Gx2yz_S_d+CDZ*I_ERI_F3z_S_Fx2y_S_d;
  Double I_ERI_F3x_S_Fxyz_Pz_d = I_ERI_F3x_S_Gxy2z_S_d+CDZ*I_ERI_F3x_S_Fxyz_S_d;
  Double I_ERI_F2xy_S_Fxyz_Pz_d = I_ERI_F2xy_S_Gxy2z_S_d+CDZ*I_ERI_F2xy_S_Fxyz_S_d;
  Double I_ERI_F2xz_S_Fxyz_Pz_d = I_ERI_F2xz_S_Gxy2z_S_d+CDZ*I_ERI_F2xz_S_Fxyz_S_d;
  Double I_ERI_Fx2y_S_Fxyz_Pz_d = I_ERI_Fx2y_S_Gxy2z_S_d+CDZ*I_ERI_Fx2y_S_Fxyz_S_d;
  Double I_ERI_Fxyz_S_Fxyz_Pz_d = I_ERI_Fxyz_S_Gxy2z_S_d+CDZ*I_ERI_Fxyz_S_Fxyz_S_d;
  Double I_ERI_Fx2z_S_Fxyz_Pz_d = I_ERI_Fx2z_S_Gxy2z_S_d+CDZ*I_ERI_Fx2z_S_Fxyz_S_d;
  Double I_ERI_F3y_S_Fxyz_Pz_d = I_ERI_F3y_S_Gxy2z_S_d+CDZ*I_ERI_F3y_S_Fxyz_S_d;
  Double I_ERI_F2yz_S_Fxyz_Pz_d = I_ERI_F2yz_S_Gxy2z_S_d+CDZ*I_ERI_F2yz_S_Fxyz_S_d;
  Double I_ERI_Fy2z_S_Fxyz_Pz_d = I_ERI_Fy2z_S_Gxy2z_S_d+CDZ*I_ERI_Fy2z_S_Fxyz_S_d;
  Double I_ERI_F3z_S_Fxyz_Pz_d = I_ERI_F3z_S_Gxy2z_S_d+CDZ*I_ERI_F3z_S_Fxyz_S_d;
  Double I_ERI_F3x_S_Fx2z_Pz_d = I_ERI_F3x_S_Gx3z_S_d+CDZ*I_ERI_F3x_S_Fx2z_S_d;
  Double I_ERI_F2xy_S_Fx2z_Pz_d = I_ERI_F2xy_S_Gx3z_S_d+CDZ*I_ERI_F2xy_S_Fx2z_S_d;
  Double I_ERI_F2xz_S_Fx2z_Pz_d = I_ERI_F2xz_S_Gx3z_S_d+CDZ*I_ERI_F2xz_S_Fx2z_S_d;
  Double I_ERI_Fx2y_S_Fx2z_Pz_d = I_ERI_Fx2y_S_Gx3z_S_d+CDZ*I_ERI_Fx2y_S_Fx2z_S_d;
  Double I_ERI_Fxyz_S_Fx2z_Pz_d = I_ERI_Fxyz_S_Gx3z_S_d+CDZ*I_ERI_Fxyz_S_Fx2z_S_d;
  Double I_ERI_Fx2z_S_Fx2z_Pz_d = I_ERI_Fx2z_S_Gx3z_S_d+CDZ*I_ERI_Fx2z_S_Fx2z_S_d;
  Double I_ERI_F3y_S_Fx2z_Pz_d = I_ERI_F3y_S_Gx3z_S_d+CDZ*I_ERI_F3y_S_Fx2z_S_d;
  Double I_ERI_F2yz_S_Fx2z_Pz_d = I_ERI_F2yz_S_Gx3z_S_d+CDZ*I_ERI_F2yz_S_Fx2z_S_d;
  Double I_ERI_Fy2z_S_Fx2z_Pz_d = I_ERI_Fy2z_S_Gx3z_S_d+CDZ*I_ERI_Fy2z_S_Fx2z_S_d;
  Double I_ERI_F3z_S_Fx2z_Pz_d = I_ERI_F3z_S_Gx3z_S_d+CDZ*I_ERI_F3z_S_Fx2z_S_d;
  Double I_ERI_F3x_S_F3y_Pz_d = I_ERI_F3x_S_G3yz_S_d+CDZ*I_ERI_F3x_S_F3y_S_d;
  Double I_ERI_F2xy_S_F3y_Pz_d = I_ERI_F2xy_S_G3yz_S_d+CDZ*I_ERI_F2xy_S_F3y_S_d;
  Double I_ERI_F2xz_S_F3y_Pz_d = I_ERI_F2xz_S_G3yz_S_d+CDZ*I_ERI_F2xz_S_F3y_S_d;
  Double I_ERI_Fx2y_S_F3y_Pz_d = I_ERI_Fx2y_S_G3yz_S_d+CDZ*I_ERI_Fx2y_S_F3y_S_d;
  Double I_ERI_Fxyz_S_F3y_Pz_d = I_ERI_Fxyz_S_G3yz_S_d+CDZ*I_ERI_Fxyz_S_F3y_S_d;
  Double I_ERI_Fx2z_S_F3y_Pz_d = I_ERI_Fx2z_S_G3yz_S_d+CDZ*I_ERI_Fx2z_S_F3y_S_d;
  Double I_ERI_F3y_S_F3y_Pz_d = I_ERI_F3y_S_G3yz_S_d+CDZ*I_ERI_F3y_S_F3y_S_d;
  Double I_ERI_F2yz_S_F3y_Pz_d = I_ERI_F2yz_S_G3yz_S_d+CDZ*I_ERI_F2yz_S_F3y_S_d;
  Double I_ERI_Fy2z_S_F3y_Pz_d = I_ERI_Fy2z_S_G3yz_S_d+CDZ*I_ERI_Fy2z_S_F3y_S_d;
  Double I_ERI_F3z_S_F3y_Pz_d = I_ERI_F3z_S_G3yz_S_d+CDZ*I_ERI_F3z_S_F3y_S_d;
  Double I_ERI_F3x_S_F2yz_Pz_d = I_ERI_F3x_S_G2y2z_S_d+CDZ*I_ERI_F3x_S_F2yz_S_d;
  Double I_ERI_F2xy_S_F2yz_Pz_d = I_ERI_F2xy_S_G2y2z_S_d+CDZ*I_ERI_F2xy_S_F2yz_S_d;
  Double I_ERI_F2xz_S_F2yz_Pz_d = I_ERI_F2xz_S_G2y2z_S_d+CDZ*I_ERI_F2xz_S_F2yz_S_d;
  Double I_ERI_Fx2y_S_F2yz_Pz_d = I_ERI_Fx2y_S_G2y2z_S_d+CDZ*I_ERI_Fx2y_S_F2yz_S_d;
  Double I_ERI_Fxyz_S_F2yz_Pz_d = I_ERI_Fxyz_S_G2y2z_S_d+CDZ*I_ERI_Fxyz_S_F2yz_S_d;
  Double I_ERI_Fx2z_S_F2yz_Pz_d = I_ERI_Fx2z_S_G2y2z_S_d+CDZ*I_ERI_Fx2z_S_F2yz_S_d;
  Double I_ERI_F3y_S_F2yz_Pz_d = I_ERI_F3y_S_G2y2z_S_d+CDZ*I_ERI_F3y_S_F2yz_S_d;
  Double I_ERI_F2yz_S_F2yz_Pz_d = I_ERI_F2yz_S_G2y2z_S_d+CDZ*I_ERI_F2yz_S_F2yz_S_d;
  Double I_ERI_Fy2z_S_F2yz_Pz_d = I_ERI_Fy2z_S_G2y2z_S_d+CDZ*I_ERI_Fy2z_S_F2yz_S_d;
  Double I_ERI_F3z_S_F2yz_Pz_d = I_ERI_F3z_S_G2y2z_S_d+CDZ*I_ERI_F3z_S_F2yz_S_d;
  Double I_ERI_F3x_S_Fy2z_Pz_d = I_ERI_F3x_S_Gy3z_S_d+CDZ*I_ERI_F3x_S_Fy2z_S_d;
  Double I_ERI_F2xy_S_Fy2z_Pz_d = I_ERI_F2xy_S_Gy3z_S_d+CDZ*I_ERI_F2xy_S_Fy2z_S_d;
  Double I_ERI_F2xz_S_Fy2z_Pz_d = I_ERI_F2xz_S_Gy3z_S_d+CDZ*I_ERI_F2xz_S_Fy2z_S_d;
  Double I_ERI_Fx2y_S_Fy2z_Pz_d = I_ERI_Fx2y_S_Gy3z_S_d+CDZ*I_ERI_Fx2y_S_Fy2z_S_d;
  Double I_ERI_Fxyz_S_Fy2z_Pz_d = I_ERI_Fxyz_S_Gy3z_S_d+CDZ*I_ERI_Fxyz_S_Fy2z_S_d;
  Double I_ERI_Fx2z_S_Fy2z_Pz_d = I_ERI_Fx2z_S_Gy3z_S_d+CDZ*I_ERI_Fx2z_S_Fy2z_S_d;
  Double I_ERI_F3y_S_Fy2z_Pz_d = I_ERI_F3y_S_Gy3z_S_d+CDZ*I_ERI_F3y_S_Fy2z_S_d;
  Double I_ERI_F2yz_S_Fy2z_Pz_d = I_ERI_F2yz_S_Gy3z_S_d+CDZ*I_ERI_F2yz_S_Fy2z_S_d;
  Double I_ERI_Fy2z_S_Fy2z_Pz_d = I_ERI_Fy2z_S_Gy3z_S_d+CDZ*I_ERI_Fy2z_S_Fy2z_S_d;
  Double I_ERI_F3z_S_Fy2z_Pz_d = I_ERI_F3z_S_Gy3z_S_d+CDZ*I_ERI_F3z_S_Fy2z_S_d;
  Double I_ERI_F3x_S_F3z_Pz_d = I_ERI_F3x_S_G4z_S_d+CDZ*I_ERI_F3x_S_F3z_S_d;
  Double I_ERI_F2xy_S_F3z_Pz_d = I_ERI_F2xy_S_G4z_S_d+CDZ*I_ERI_F2xy_S_F3z_S_d;
  Double I_ERI_F2xz_S_F3z_Pz_d = I_ERI_F2xz_S_G4z_S_d+CDZ*I_ERI_F2xz_S_F3z_S_d;
  Double I_ERI_Fx2y_S_F3z_Pz_d = I_ERI_Fx2y_S_G4z_S_d+CDZ*I_ERI_Fx2y_S_F3z_S_d;
  Double I_ERI_Fxyz_S_F3z_Pz_d = I_ERI_Fxyz_S_G4z_S_d+CDZ*I_ERI_Fxyz_S_F3z_S_d;
  Double I_ERI_Fx2z_S_F3z_Pz_d = I_ERI_Fx2z_S_G4z_S_d+CDZ*I_ERI_Fx2z_S_F3z_S_d;
  Double I_ERI_F3y_S_F3z_Pz_d = I_ERI_F3y_S_G4z_S_d+CDZ*I_ERI_F3y_S_F3z_S_d;
  Double I_ERI_F2yz_S_F3z_Pz_d = I_ERI_F2yz_S_G4z_S_d+CDZ*I_ERI_F2yz_S_F3z_S_d;
  Double I_ERI_Fy2z_S_F3z_Pz_d = I_ERI_Fy2z_S_G4z_S_d+CDZ*I_ERI_Fy2z_S_F3z_S_d;
  Double I_ERI_F3z_S_F3z_Pz_d = I_ERI_F3z_S_G4z_S_d+CDZ*I_ERI_F3z_S_F3z_S_d;

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
   * shell quartet name: SQ_ERI_P_P_F_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S
   * RHS shell quartet name: SQ_ERI_P_S_F_S
   ************************************************************/
  Double I_ERI_Px_Px_F3x_S = I_ERI_D2x_S_F3x_S+ABX*I_ERI_Px_S_F3x_S;
  Double I_ERI_Py_Px_F3x_S = I_ERI_Dxy_S_F3x_S+ABX*I_ERI_Py_S_F3x_S;
  Double I_ERI_Pz_Px_F3x_S = I_ERI_Dxz_S_F3x_S+ABX*I_ERI_Pz_S_F3x_S;
  Double I_ERI_Px_Py_F3x_S = I_ERI_Dxy_S_F3x_S+ABY*I_ERI_Px_S_F3x_S;
  Double I_ERI_Py_Py_F3x_S = I_ERI_D2y_S_F3x_S+ABY*I_ERI_Py_S_F3x_S;
  Double I_ERI_Pz_Py_F3x_S = I_ERI_Dyz_S_F3x_S+ABY*I_ERI_Pz_S_F3x_S;
  Double I_ERI_Px_Pz_F3x_S = I_ERI_Dxz_S_F3x_S+ABZ*I_ERI_Px_S_F3x_S;
  Double I_ERI_Py_Pz_F3x_S = I_ERI_Dyz_S_F3x_S+ABZ*I_ERI_Py_S_F3x_S;
  Double I_ERI_Pz_Pz_F3x_S = I_ERI_D2z_S_F3x_S+ABZ*I_ERI_Pz_S_F3x_S;
  Double I_ERI_Px_Px_F2xy_S = I_ERI_D2x_S_F2xy_S+ABX*I_ERI_Px_S_F2xy_S;
  Double I_ERI_Py_Px_F2xy_S = I_ERI_Dxy_S_F2xy_S+ABX*I_ERI_Py_S_F2xy_S;
  Double I_ERI_Pz_Px_F2xy_S = I_ERI_Dxz_S_F2xy_S+ABX*I_ERI_Pz_S_F2xy_S;
  Double I_ERI_Px_Py_F2xy_S = I_ERI_Dxy_S_F2xy_S+ABY*I_ERI_Px_S_F2xy_S;
  Double I_ERI_Py_Py_F2xy_S = I_ERI_D2y_S_F2xy_S+ABY*I_ERI_Py_S_F2xy_S;
  Double I_ERI_Pz_Py_F2xy_S = I_ERI_Dyz_S_F2xy_S+ABY*I_ERI_Pz_S_F2xy_S;
  Double I_ERI_Px_Pz_F2xy_S = I_ERI_Dxz_S_F2xy_S+ABZ*I_ERI_Px_S_F2xy_S;
  Double I_ERI_Py_Pz_F2xy_S = I_ERI_Dyz_S_F2xy_S+ABZ*I_ERI_Py_S_F2xy_S;
  Double I_ERI_Pz_Pz_F2xy_S = I_ERI_D2z_S_F2xy_S+ABZ*I_ERI_Pz_S_F2xy_S;
  Double I_ERI_Px_Px_F2xz_S = I_ERI_D2x_S_F2xz_S+ABX*I_ERI_Px_S_F2xz_S;
  Double I_ERI_Py_Px_F2xz_S = I_ERI_Dxy_S_F2xz_S+ABX*I_ERI_Py_S_F2xz_S;
  Double I_ERI_Pz_Px_F2xz_S = I_ERI_Dxz_S_F2xz_S+ABX*I_ERI_Pz_S_F2xz_S;
  Double I_ERI_Px_Py_F2xz_S = I_ERI_Dxy_S_F2xz_S+ABY*I_ERI_Px_S_F2xz_S;
  Double I_ERI_Py_Py_F2xz_S = I_ERI_D2y_S_F2xz_S+ABY*I_ERI_Py_S_F2xz_S;
  Double I_ERI_Pz_Py_F2xz_S = I_ERI_Dyz_S_F2xz_S+ABY*I_ERI_Pz_S_F2xz_S;
  Double I_ERI_Px_Pz_F2xz_S = I_ERI_Dxz_S_F2xz_S+ABZ*I_ERI_Px_S_F2xz_S;
  Double I_ERI_Py_Pz_F2xz_S = I_ERI_Dyz_S_F2xz_S+ABZ*I_ERI_Py_S_F2xz_S;
  Double I_ERI_Pz_Pz_F2xz_S = I_ERI_D2z_S_F2xz_S+ABZ*I_ERI_Pz_S_F2xz_S;
  Double I_ERI_Px_Px_Fx2y_S = I_ERI_D2x_S_Fx2y_S+ABX*I_ERI_Px_S_Fx2y_S;
  Double I_ERI_Py_Px_Fx2y_S = I_ERI_Dxy_S_Fx2y_S+ABX*I_ERI_Py_S_Fx2y_S;
  Double I_ERI_Pz_Px_Fx2y_S = I_ERI_Dxz_S_Fx2y_S+ABX*I_ERI_Pz_S_Fx2y_S;
  Double I_ERI_Px_Py_Fx2y_S = I_ERI_Dxy_S_Fx2y_S+ABY*I_ERI_Px_S_Fx2y_S;
  Double I_ERI_Py_Py_Fx2y_S = I_ERI_D2y_S_Fx2y_S+ABY*I_ERI_Py_S_Fx2y_S;
  Double I_ERI_Pz_Py_Fx2y_S = I_ERI_Dyz_S_Fx2y_S+ABY*I_ERI_Pz_S_Fx2y_S;
  Double I_ERI_Px_Pz_Fx2y_S = I_ERI_Dxz_S_Fx2y_S+ABZ*I_ERI_Px_S_Fx2y_S;
  Double I_ERI_Py_Pz_Fx2y_S = I_ERI_Dyz_S_Fx2y_S+ABZ*I_ERI_Py_S_Fx2y_S;
  Double I_ERI_Pz_Pz_Fx2y_S = I_ERI_D2z_S_Fx2y_S+ABZ*I_ERI_Pz_S_Fx2y_S;
  Double I_ERI_Px_Px_Fxyz_S = I_ERI_D2x_S_Fxyz_S+ABX*I_ERI_Px_S_Fxyz_S;
  Double I_ERI_Py_Px_Fxyz_S = I_ERI_Dxy_S_Fxyz_S+ABX*I_ERI_Py_S_Fxyz_S;
  Double I_ERI_Pz_Px_Fxyz_S = I_ERI_Dxz_S_Fxyz_S+ABX*I_ERI_Pz_S_Fxyz_S;
  Double I_ERI_Px_Py_Fxyz_S = I_ERI_Dxy_S_Fxyz_S+ABY*I_ERI_Px_S_Fxyz_S;
  Double I_ERI_Py_Py_Fxyz_S = I_ERI_D2y_S_Fxyz_S+ABY*I_ERI_Py_S_Fxyz_S;
  Double I_ERI_Pz_Py_Fxyz_S = I_ERI_Dyz_S_Fxyz_S+ABY*I_ERI_Pz_S_Fxyz_S;
  Double I_ERI_Px_Pz_Fxyz_S = I_ERI_Dxz_S_Fxyz_S+ABZ*I_ERI_Px_S_Fxyz_S;
  Double I_ERI_Py_Pz_Fxyz_S = I_ERI_Dyz_S_Fxyz_S+ABZ*I_ERI_Py_S_Fxyz_S;
  Double I_ERI_Pz_Pz_Fxyz_S = I_ERI_D2z_S_Fxyz_S+ABZ*I_ERI_Pz_S_Fxyz_S;
  Double I_ERI_Px_Px_Fx2z_S = I_ERI_D2x_S_Fx2z_S+ABX*I_ERI_Px_S_Fx2z_S;
  Double I_ERI_Py_Px_Fx2z_S = I_ERI_Dxy_S_Fx2z_S+ABX*I_ERI_Py_S_Fx2z_S;
  Double I_ERI_Pz_Px_Fx2z_S = I_ERI_Dxz_S_Fx2z_S+ABX*I_ERI_Pz_S_Fx2z_S;
  Double I_ERI_Px_Py_Fx2z_S = I_ERI_Dxy_S_Fx2z_S+ABY*I_ERI_Px_S_Fx2z_S;
  Double I_ERI_Py_Py_Fx2z_S = I_ERI_D2y_S_Fx2z_S+ABY*I_ERI_Py_S_Fx2z_S;
  Double I_ERI_Pz_Py_Fx2z_S = I_ERI_Dyz_S_Fx2z_S+ABY*I_ERI_Pz_S_Fx2z_S;
  Double I_ERI_Px_Pz_Fx2z_S = I_ERI_Dxz_S_Fx2z_S+ABZ*I_ERI_Px_S_Fx2z_S;
  Double I_ERI_Py_Pz_Fx2z_S = I_ERI_Dyz_S_Fx2z_S+ABZ*I_ERI_Py_S_Fx2z_S;
  Double I_ERI_Pz_Pz_Fx2z_S = I_ERI_D2z_S_Fx2z_S+ABZ*I_ERI_Pz_S_Fx2z_S;
  Double I_ERI_Px_Px_F3y_S = I_ERI_D2x_S_F3y_S+ABX*I_ERI_Px_S_F3y_S;
  Double I_ERI_Py_Px_F3y_S = I_ERI_Dxy_S_F3y_S+ABX*I_ERI_Py_S_F3y_S;
  Double I_ERI_Pz_Px_F3y_S = I_ERI_Dxz_S_F3y_S+ABX*I_ERI_Pz_S_F3y_S;
  Double I_ERI_Px_Py_F3y_S = I_ERI_Dxy_S_F3y_S+ABY*I_ERI_Px_S_F3y_S;
  Double I_ERI_Py_Py_F3y_S = I_ERI_D2y_S_F3y_S+ABY*I_ERI_Py_S_F3y_S;
  Double I_ERI_Pz_Py_F3y_S = I_ERI_Dyz_S_F3y_S+ABY*I_ERI_Pz_S_F3y_S;
  Double I_ERI_Px_Pz_F3y_S = I_ERI_Dxz_S_F3y_S+ABZ*I_ERI_Px_S_F3y_S;
  Double I_ERI_Py_Pz_F3y_S = I_ERI_Dyz_S_F3y_S+ABZ*I_ERI_Py_S_F3y_S;
  Double I_ERI_Pz_Pz_F3y_S = I_ERI_D2z_S_F3y_S+ABZ*I_ERI_Pz_S_F3y_S;
  Double I_ERI_Px_Px_F2yz_S = I_ERI_D2x_S_F2yz_S+ABX*I_ERI_Px_S_F2yz_S;
  Double I_ERI_Py_Px_F2yz_S = I_ERI_Dxy_S_F2yz_S+ABX*I_ERI_Py_S_F2yz_S;
  Double I_ERI_Pz_Px_F2yz_S = I_ERI_Dxz_S_F2yz_S+ABX*I_ERI_Pz_S_F2yz_S;
  Double I_ERI_Px_Py_F2yz_S = I_ERI_Dxy_S_F2yz_S+ABY*I_ERI_Px_S_F2yz_S;
  Double I_ERI_Py_Py_F2yz_S = I_ERI_D2y_S_F2yz_S+ABY*I_ERI_Py_S_F2yz_S;
  Double I_ERI_Pz_Py_F2yz_S = I_ERI_Dyz_S_F2yz_S+ABY*I_ERI_Pz_S_F2yz_S;
  Double I_ERI_Px_Pz_F2yz_S = I_ERI_Dxz_S_F2yz_S+ABZ*I_ERI_Px_S_F2yz_S;
  Double I_ERI_Py_Pz_F2yz_S = I_ERI_Dyz_S_F2yz_S+ABZ*I_ERI_Py_S_F2yz_S;
  Double I_ERI_Pz_Pz_F2yz_S = I_ERI_D2z_S_F2yz_S+ABZ*I_ERI_Pz_S_F2yz_S;
  Double I_ERI_Px_Px_Fy2z_S = I_ERI_D2x_S_Fy2z_S+ABX*I_ERI_Px_S_Fy2z_S;
  Double I_ERI_Py_Px_Fy2z_S = I_ERI_Dxy_S_Fy2z_S+ABX*I_ERI_Py_S_Fy2z_S;
  Double I_ERI_Pz_Px_Fy2z_S = I_ERI_Dxz_S_Fy2z_S+ABX*I_ERI_Pz_S_Fy2z_S;
  Double I_ERI_Px_Py_Fy2z_S = I_ERI_Dxy_S_Fy2z_S+ABY*I_ERI_Px_S_Fy2z_S;
  Double I_ERI_Py_Py_Fy2z_S = I_ERI_D2y_S_Fy2z_S+ABY*I_ERI_Py_S_Fy2z_S;
  Double I_ERI_Pz_Py_Fy2z_S = I_ERI_Dyz_S_Fy2z_S+ABY*I_ERI_Pz_S_Fy2z_S;
  Double I_ERI_Px_Pz_Fy2z_S = I_ERI_Dxz_S_Fy2z_S+ABZ*I_ERI_Px_S_Fy2z_S;
  Double I_ERI_Py_Pz_Fy2z_S = I_ERI_Dyz_S_Fy2z_S+ABZ*I_ERI_Py_S_Fy2z_S;
  Double I_ERI_Pz_Pz_Fy2z_S = I_ERI_D2z_S_Fy2z_S+ABZ*I_ERI_Pz_S_Fy2z_S;
  Double I_ERI_Px_Px_F3z_S = I_ERI_D2x_S_F3z_S+ABX*I_ERI_Px_S_F3z_S;
  Double I_ERI_Py_Px_F3z_S = I_ERI_Dxy_S_F3z_S+ABX*I_ERI_Py_S_F3z_S;
  Double I_ERI_Pz_Px_F3z_S = I_ERI_Dxz_S_F3z_S+ABX*I_ERI_Pz_S_F3z_S;
  Double I_ERI_Px_Py_F3z_S = I_ERI_Dxy_S_F3z_S+ABY*I_ERI_Px_S_F3z_S;
  Double I_ERI_Py_Py_F3z_S = I_ERI_D2y_S_F3z_S+ABY*I_ERI_Py_S_F3z_S;
  Double I_ERI_Pz_Py_F3z_S = I_ERI_Dyz_S_F3z_S+ABY*I_ERI_Pz_S_F3z_S;
  Double I_ERI_Px_Pz_F3z_S = I_ERI_Dxz_S_F3z_S+ABZ*I_ERI_Px_S_F3z_S;
  Double I_ERI_Py_Pz_F3z_S = I_ERI_Dyz_S_F3z_S+ABZ*I_ERI_Py_S_F3z_S;
  Double I_ERI_Pz_Pz_F3z_S = I_ERI_D2z_S_F3z_S+ABZ*I_ERI_Pz_S_F3z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_F_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_a
   * RHS shell quartet name: SQ_ERI_F_S_F_S_a
   ************************************************************/
  Double I_ERI_F3x_Px_F3x_S_a = I_ERI_G4x_S_F3x_S_a+ABX*I_ERI_F3x_S_F3x_S_a;
  Double I_ERI_F2xy_Px_F3x_S_a = I_ERI_G3xy_S_F3x_S_a+ABX*I_ERI_F2xy_S_F3x_S_a;
  Double I_ERI_F2xz_Px_F3x_S_a = I_ERI_G3xz_S_F3x_S_a+ABX*I_ERI_F2xz_S_F3x_S_a;
  Double I_ERI_Fx2y_Px_F3x_S_a = I_ERI_G2x2y_S_F3x_S_a+ABX*I_ERI_Fx2y_S_F3x_S_a;
  Double I_ERI_Fxyz_Px_F3x_S_a = I_ERI_G2xyz_S_F3x_S_a+ABX*I_ERI_Fxyz_S_F3x_S_a;
  Double I_ERI_Fx2z_Px_F3x_S_a = I_ERI_G2x2z_S_F3x_S_a+ABX*I_ERI_Fx2z_S_F3x_S_a;
  Double I_ERI_F3y_Px_F3x_S_a = I_ERI_Gx3y_S_F3x_S_a+ABX*I_ERI_F3y_S_F3x_S_a;
  Double I_ERI_F2yz_Px_F3x_S_a = I_ERI_Gx2yz_S_F3x_S_a+ABX*I_ERI_F2yz_S_F3x_S_a;
  Double I_ERI_Fy2z_Px_F3x_S_a = I_ERI_Gxy2z_S_F3x_S_a+ABX*I_ERI_Fy2z_S_F3x_S_a;
  Double I_ERI_F3z_Px_F3x_S_a = I_ERI_Gx3z_S_F3x_S_a+ABX*I_ERI_F3z_S_F3x_S_a;
  Double I_ERI_F3x_Py_F3x_S_a = I_ERI_G3xy_S_F3x_S_a+ABY*I_ERI_F3x_S_F3x_S_a;
  Double I_ERI_F2xy_Py_F3x_S_a = I_ERI_G2x2y_S_F3x_S_a+ABY*I_ERI_F2xy_S_F3x_S_a;
  Double I_ERI_F2xz_Py_F3x_S_a = I_ERI_G2xyz_S_F3x_S_a+ABY*I_ERI_F2xz_S_F3x_S_a;
  Double I_ERI_Fx2y_Py_F3x_S_a = I_ERI_Gx3y_S_F3x_S_a+ABY*I_ERI_Fx2y_S_F3x_S_a;
  Double I_ERI_Fxyz_Py_F3x_S_a = I_ERI_Gx2yz_S_F3x_S_a+ABY*I_ERI_Fxyz_S_F3x_S_a;
  Double I_ERI_Fx2z_Py_F3x_S_a = I_ERI_Gxy2z_S_F3x_S_a+ABY*I_ERI_Fx2z_S_F3x_S_a;
  Double I_ERI_F3y_Py_F3x_S_a = I_ERI_G4y_S_F3x_S_a+ABY*I_ERI_F3y_S_F3x_S_a;
  Double I_ERI_F2yz_Py_F3x_S_a = I_ERI_G3yz_S_F3x_S_a+ABY*I_ERI_F2yz_S_F3x_S_a;
  Double I_ERI_Fy2z_Py_F3x_S_a = I_ERI_G2y2z_S_F3x_S_a+ABY*I_ERI_Fy2z_S_F3x_S_a;
  Double I_ERI_F3z_Py_F3x_S_a = I_ERI_Gy3z_S_F3x_S_a+ABY*I_ERI_F3z_S_F3x_S_a;
  Double I_ERI_F3x_Pz_F3x_S_a = I_ERI_G3xz_S_F3x_S_a+ABZ*I_ERI_F3x_S_F3x_S_a;
  Double I_ERI_F2xy_Pz_F3x_S_a = I_ERI_G2xyz_S_F3x_S_a+ABZ*I_ERI_F2xy_S_F3x_S_a;
  Double I_ERI_F2xz_Pz_F3x_S_a = I_ERI_G2x2z_S_F3x_S_a+ABZ*I_ERI_F2xz_S_F3x_S_a;
  Double I_ERI_Fx2y_Pz_F3x_S_a = I_ERI_Gx2yz_S_F3x_S_a+ABZ*I_ERI_Fx2y_S_F3x_S_a;
  Double I_ERI_Fxyz_Pz_F3x_S_a = I_ERI_Gxy2z_S_F3x_S_a+ABZ*I_ERI_Fxyz_S_F3x_S_a;
  Double I_ERI_Fx2z_Pz_F3x_S_a = I_ERI_Gx3z_S_F3x_S_a+ABZ*I_ERI_Fx2z_S_F3x_S_a;
  Double I_ERI_F3y_Pz_F3x_S_a = I_ERI_G3yz_S_F3x_S_a+ABZ*I_ERI_F3y_S_F3x_S_a;
  Double I_ERI_F2yz_Pz_F3x_S_a = I_ERI_G2y2z_S_F3x_S_a+ABZ*I_ERI_F2yz_S_F3x_S_a;
  Double I_ERI_Fy2z_Pz_F3x_S_a = I_ERI_Gy3z_S_F3x_S_a+ABZ*I_ERI_Fy2z_S_F3x_S_a;
  Double I_ERI_F3z_Pz_F3x_S_a = I_ERI_G4z_S_F3x_S_a+ABZ*I_ERI_F3z_S_F3x_S_a;
  Double I_ERI_F3x_Px_F2xy_S_a = I_ERI_G4x_S_F2xy_S_a+ABX*I_ERI_F3x_S_F2xy_S_a;
  Double I_ERI_F2xy_Px_F2xy_S_a = I_ERI_G3xy_S_F2xy_S_a+ABX*I_ERI_F2xy_S_F2xy_S_a;
  Double I_ERI_F2xz_Px_F2xy_S_a = I_ERI_G3xz_S_F2xy_S_a+ABX*I_ERI_F2xz_S_F2xy_S_a;
  Double I_ERI_Fx2y_Px_F2xy_S_a = I_ERI_G2x2y_S_F2xy_S_a+ABX*I_ERI_Fx2y_S_F2xy_S_a;
  Double I_ERI_Fxyz_Px_F2xy_S_a = I_ERI_G2xyz_S_F2xy_S_a+ABX*I_ERI_Fxyz_S_F2xy_S_a;
  Double I_ERI_Fx2z_Px_F2xy_S_a = I_ERI_G2x2z_S_F2xy_S_a+ABX*I_ERI_Fx2z_S_F2xy_S_a;
  Double I_ERI_F3y_Px_F2xy_S_a = I_ERI_Gx3y_S_F2xy_S_a+ABX*I_ERI_F3y_S_F2xy_S_a;
  Double I_ERI_F2yz_Px_F2xy_S_a = I_ERI_Gx2yz_S_F2xy_S_a+ABX*I_ERI_F2yz_S_F2xy_S_a;
  Double I_ERI_Fy2z_Px_F2xy_S_a = I_ERI_Gxy2z_S_F2xy_S_a+ABX*I_ERI_Fy2z_S_F2xy_S_a;
  Double I_ERI_F3z_Px_F2xy_S_a = I_ERI_Gx3z_S_F2xy_S_a+ABX*I_ERI_F3z_S_F2xy_S_a;
  Double I_ERI_F3x_Py_F2xy_S_a = I_ERI_G3xy_S_F2xy_S_a+ABY*I_ERI_F3x_S_F2xy_S_a;
  Double I_ERI_F2xy_Py_F2xy_S_a = I_ERI_G2x2y_S_F2xy_S_a+ABY*I_ERI_F2xy_S_F2xy_S_a;
  Double I_ERI_F2xz_Py_F2xy_S_a = I_ERI_G2xyz_S_F2xy_S_a+ABY*I_ERI_F2xz_S_F2xy_S_a;
  Double I_ERI_Fx2y_Py_F2xy_S_a = I_ERI_Gx3y_S_F2xy_S_a+ABY*I_ERI_Fx2y_S_F2xy_S_a;
  Double I_ERI_Fxyz_Py_F2xy_S_a = I_ERI_Gx2yz_S_F2xy_S_a+ABY*I_ERI_Fxyz_S_F2xy_S_a;
  Double I_ERI_Fx2z_Py_F2xy_S_a = I_ERI_Gxy2z_S_F2xy_S_a+ABY*I_ERI_Fx2z_S_F2xy_S_a;
  Double I_ERI_F3y_Py_F2xy_S_a = I_ERI_G4y_S_F2xy_S_a+ABY*I_ERI_F3y_S_F2xy_S_a;
  Double I_ERI_F2yz_Py_F2xy_S_a = I_ERI_G3yz_S_F2xy_S_a+ABY*I_ERI_F2yz_S_F2xy_S_a;
  Double I_ERI_Fy2z_Py_F2xy_S_a = I_ERI_G2y2z_S_F2xy_S_a+ABY*I_ERI_Fy2z_S_F2xy_S_a;
  Double I_ERI_F3z_Py_F2xy_S_a = I_ERI_Gy3z_S_F2xy_S_a+ABY*I_ERI_F3z_S_F2xy_S_a;
  Double I_ERI_F3x_Pz_F2xy_S_a = I_ERI_G3xz_S_F2xy_S_a+ABZ*I_ERI_F3x_S_F2xy_S_a;
  Double I_ERI_F2xy_Pz_F2xy_S_a = I_ERI_G2xyz_S_F2xy_S_a+ABZ*I_ERI_F2xy_S_F2xy_S_a;
  Double I_ERI_F2xz_Pz_F2xy_S_a = I_ERI_G2x2z_S_F2xy_S_a+ABZ*I_ERI_F2xz_S_F2xy_S_a;
  Double I_ERI_Fx2y_Pz_F2xy_S_a = I_ERI_Gx2yz_S_F2xy_S_a+ABZ*I_ERI_Fx2y_S_F2xy_S_a;
  Double I_ERI_Fxyz_Pz_F2xy_S_a = I_ERI_Gxy2z_S_F2xy_S_a+ABZ*I_ERI_Fxyz_S_F2xy_S_a;
  Double I_ERI_Fx2z_Pz_F2xy_S_a = I_ERI_Gx3z_S_F2xy_S_a+ABZ*I_ERI_Fx2z_S_F2xy_S_a;
  Double I_ERI_F3y_Pz_F2xy_S_a = I_ERI_G3yz_S_F2xy_S_a+ABZ*I_ERI_F3y_S_F2xy_S_a;
  Double I_ERI_F2yz_Pz_F2xy_S_a = I_ERI_G2y2z_S_F2xy_S_a+ABZ*I_ERI_F2yz_S_F2xy_S_a;
  Double I_ERI_Fy2z_Pz_F2xy_S_a = I_ERI_Gy3z_S_F2xy_S_a+ABZ*I_ERI_Fy2z_S_F2xy_S_a;
  Double I_ERI_F3z_Pz_F2xy_S_a = I_ERI_G4z_S_F2xy_S_a+ABZ*I_ERI_F3z_S_F2xy_S_a;
  Double I_ERI_F3x_Px_F2xz_S_a = I_ERI_G4x_S_F2xz_S_a+ABX*I_ERI_F3x_S_F2xz_S_a;
  Double I_ERI_F2xy_Px_F2xz_S_a = I_ERI_G3xy_S_F2xz_S_a+ABX*I_ERI_F2xy_S_F2xz_S_a;
  Double I_ERI_F2xz_Px_F2xz_S_a = I_ERI_G3xz_S_F2xz_S_a+ABX*I_ERI_F2xz_S_F2xz_S_a;
  Double I_ERI_Fx2y_Px_F2xz_S_a = I_ERI_G2x2y_S_F2xz_S_a+ABX*I_ERI_Fx2y_S_F2xz_S_a;
  Double I_ERI_Fxyz_Px_F2xz_S_a = I_ERI_G2xyz_S_F2xz_S_a+ABX*I_ERI_Fxyz_S_F2xz_S_a;
  Double I_ERI_Fx2z_Px_F2xz_S_a = I_ERI_G2x2z_S_F2xz_S_a+ABX*I_ERI_Fx2z_S_F2xz_S_a;
  Double I_ERI_F3y_Px_F2xz_S_a = I_ERI_Gx3y_S_F2xz_S_a+ABX*I_ERI_F3y_S_F2xz_S_a;
  Double I_ERI_F2yz_Px_F2xz_S_a = I_ERI_Gx2yz_S_F2xz_S_a+ABX*I_ERI_F2yz_S_F2xz_S_a;
  Double I_ERI_Fy2z_Px_F2xz_S_a = I_ERI_Gxy2z_S_F2xz_S_a+ABX*I_ERI_Fy2z_S_F2xz_S_a;
  Double I_ERI_F3z_Px_F2xz_S_a = I_ERI_Gx3z_S_F2xz_S_a+ABX*I_ERI_F3z_S_F2xz_S_a;
  Double I_ERI_F3x_Py_F2xz_S_a = I_ERI_G3xy_S_F2xz_S_a+ABY*I_ERI_F3x_S_F2xz_S_a;
  Double I_ERI_F2xy_Py_F2xz_S_a = I_ERI_G2x2y_S_F2xz_S_a+ABY*I_ERI_F2xy_S_F2xz_S_a;
  Double I_ERI_F2xz_Py_F2xz_S_a = I_ERI_G2xyz_S_F2xz_S_a+ABY*I_ERI_F2xz_S_F2xz_S_a;
  Double I_ERI_Fx2y_Py_F2xz_S_a = I_ERI_Gx3y_S_F2xz_S_a+ABY*I_ERI_Fx2y_S_F2xz_S_a;
  Double I_ERI_Fxyz_Py_F2xz_S_a = I_ERI_Gx2yz_S_F2xz_S_a+ABY*I_ERI_Fxyz_S_F2xz_S_a;
  Double I_ERI_Fx2z_Py_F2xz_S_a = I_ERI_Gxy2z_S_F2xz_S_a+ABY*I_ERI_Fx2z_S_F2xz_S_a;
  Double I_ERI_F3y_Py_F2xz_S_a = I_ERI_G4y_S_F2xz_S_a+ABY*I_ERI_F3y_S_F2xz_S_a;
  Double I_ERI_F2yz_Py_F2xz_S_a = I_ERI_G3yz_S_F2xz_S_a+ABY*I_ERI_F2yz_S_F2xz_S_a;
  Double I_ERI_Fy2z_Py_F2xz_S_a = I_ERI_G2y2z_S_F2xz_S_a+ABY*I_ERI_Fy2z_S_F2xz_S_a;
  Double I_ERI_F3z_Py_F2xz_S_a = I_ERI_Gy3z_S_F2xz_S_a+ABY*I_ERI_F3z_S_F2xz_S_a;
  Double I_ERI_F3x_Pz_F2xz_S_a = I_ERI_G3xz_S_F2xz_S_a+ABZ*I_ERI_F3x_S_F2xz_S_a;
  Double I_ERI_F2xy_Pz_F2xz_S_a = I_ERI_G2xyz_S_F2xz_S_a+ABZ*I_ERI_F2xy_S_F2xz_S_a;
  Double I_ERI_F2xz_Pz_F2xz_S_a = I_ERI_G2x2z_S_F2xz_S_a+ABZ*I_ERI_F2xz_S_F2xz_S_a;
  Double I_ERI_Fx2y_Pz_F2xz_S_a = I_ERI_Gx2yz_S_F2xz_S_a+ABZ*I_ERI_Fx2y_S_F2xz_S_a;
  Double I_ERI_Fxyz_Pz_F2xz_S_a = I_ERI_Gxy2z_S_F2xz_S_a+ABZ*I_ERI_Fxyz_S_F2xz_S_a;
  Double I_ERI_Fx2z_Pz_F2xz_S_a = I_ERI_Gx3z_S_F2xz_S_a+ABZ*I_ERI_Fx2z_S_F2xz_S_a;
  Double I_ERI_F3y_Pz_F2xz_S_a = I_ERI_G3yz_S_F2xz_S_a+ABZ*I_ERI_F3y_S_F2xz_S_a;
  Double I_ERI_F2yz_Pz_F2xz_S_a = I_ERI_G2y2z_S_F2xz_S_a+ABZ*I_ERI_F2yz_S_F2xz_S_a;
  Double I_ERI_Fy2z_Pz_F2xz_S_a = I_ERI_Gy3z_S_F2xz_S_a+ABZ*I_ERI_Fy2z_S_F2xz_S_a;
  Double I_ERI_F3z_Pz_F2xz_S_a = I_ERI_G4z_S_F2xz_S_a+ABZ*I_ERI_F3z_S_F2xz_S_a;
  Double I_ERI_F3x_Px_Fx2y_S_a = I_ERI_G4x_S_Fx2y_S_a+ABX*I_ERI_F3x_S_Fx2y_S_a;
  Double I_ERI_F2xy_Px_Fx2y_S_a = I_ERI_G3xy_S_Fx2y_S_a+ABX*I_ERI_F2xy_S_Fx2y_S_a;
  Double I_ERI_F2xz_Px_Fx2y_S_a = I_ERI_G3xz_S_Fx2y_S_a+ABX*I_ERI_F2xz_S_Fx2y_S_a;
  Double I_ERI_Fx2y_Px_Fx2y_S_a = I_ERI_G2x2y_S_Fx2y_S_a+ABX*I_ERI_Fx2y_S_Fx2y_S_a;
  Double I_ERI_Fxyz_Px_Fx2y_S_a = I_ERI_G2xyz_S_Fx2y_S_a+ABX*I_ERI_Fxyz_S_Fx2y_S_a;
  Double I_ERI_Fx2z_Px_Fx2y_S_a = I_ERI_G2x2z_S_Fx2y_S_a+ABX*I_ERI_Fx2z_S_Fx2y_S_a;
  Double I_ERI_F3y_Px_Fx2y_S_a = I_ERI_Gx3y_S_Fx2y_S_a+ABX*I_ERI_F3y_S_Fx2y_S_a;
  Double I_ERI_F2yz_Px_Fx2y_S_a = I_ERI_Gx2yz_S_Fx2y_S_a+ABX*I_ERI_F2yz_S_Fx2y_S_a;
  Double I_ERI_Fy2z_Px_Fx2y_S_a = I_ERI_Gxy2z_S_Fx2y_S_a+ABX*I_ERI_Fy2z_S_Fx2y_S_a;
  Double I_ERI_F3z_Px_Fx2y_S_a = I_ERI_Gx3z_S_Fx2y_S_a+ABX*I_ERI_F3z_S_Fx2y_S_a;
  Double I_ERI_F3x_Py_Fx2y_S_a = I_ERI_G3xy_S_Fx2y_S_a+ABY*I_ERI_F3x_S_Fx2y_S_a;
  Double I_ERI_F2xy_Py_Fx2y_S_a = I_ERI_G2x2y_S_Fx2y_S_a+ABY*I_ERI_F2xy_S_Fx2y_S_a;
  Double I_ERI_F2xz_Py_Fx2y_S_a = I_ERI_G2xyz_S_Fx2y_S_a+ABY*I_ERI_F2xz_S_Fx2y_S_a;
  Double I_ERI_Fx2y_Py_Fx2y_S_a = I_ERI_Gx3y_S_Fx2y_S_a+ABY*I_ERI_Fx2y_S_Fx2y_S_a;
  Double I_ERI_Fxyz_Py_Fx2y_S_a = I_ERI_Gx2yz_S_Fx2y_S_a+ABY*I_ERI_Fxyz_S_Fx2y_S_a;
  Double I_ERI_Fx2z_Py_Fx2y_S_a = I_ERI_Gxy2z_S_Fx2y_S_a+ABY*I_ERI_Fx2z_S_Fx2y_S_a;
  Double I_ERI_F3y_Py_Fx2y_S_a = I_ERI_G4y_S_Fx2y_S_a+ABY*I_ERI_F3y_S_Fx2y_S_a;
  Double I_ERI_F2yz_Py_Fx2y_S_a = I_ERI_G3yz_S_Fx2y_S_a+ABY*I_ERI_F2yz_S_Fx2y_S_a;
  Double I_ERI_Fy2z_Py_Fx2y_S_a = I_ERI_G2y2z_S_Fx2y_S_a+ABY*I_ERI_Fy2z_S_Fx2y_S_a;
  Double I_ERI_F3z_Py_Fx2y_S_a = I_ERI_Gy3z_S_Fx2y_S_a+ABY*I_ERI_F3z_S_Fx2y_S_a;
  Double I_ERI_F3x_Pz_Fx2y_S_a = I_ERI_G3xz_S_Fx2y_S_a+ABZ*I_ERI_F3x_S_Fx2y_S_a;
  Double I_ERI_F2xy_Pz_Fx2y_S_a = I_ERI_G2xyz_S_Fx2y_S_a+ABZ*I_ERI_F2xy_S_Fx2y_S_a;
  Double I_ERI_F2xz_Pz_Fx2y_S_a = I_ERI_G2x2z_S_Fx2y_S_a+ABZ*I_ERI_F2xz_S_Fx2y_S_a;
  Double I_ERI_Fx2y_Pz_Fx2y_S_a = I_ERI_Gx2yz_S_Fx2y_S_a+ABZ*I_ERI_Fx2y_S_Fx2y_S_a;
  Double I_ERI_Fxyz_Pz_Fx2y_S_a = I_ERI_Gxy2z_S_Fx2y_S_a+ABZ*I_ERI_Fxyz_S_Fx2y_S_a;
  Double I_ERI_Fx2z_Pz_Fx2y_S_a = I_ERI_Gx3z_S_Fx2y_S_a+ABZ*I_ERI_Fx2z_S_Fx2y_S_a;
  Double I_ERI_F3y_Pz_Fx2y_S_a = I_ERI_G3yz_S_Fx2y_S_a+ABZ*I_ERI_F3y_S_Fx2y_S_a;
  Double I_ERI_F2yz_Pz_Fx2y_S_a = I_ERI_G2y2z_S_Fx2y_S_a+ABZ*I_ERI_F2yz_S_Fx2y_S_a;
  Double I_ERI_Fy2z_Pz_Fx2y_S_a = I_ERI_Gy3z_S_Fx2y_S_a+ABZ*I_ERI_Fy2z_S_Fx2y_S_a;
  Double I_ERI_F3z_Pz_Fx2y_S_a = I_ERI_G4z_S_Fx2y_S_a+ABZ*I_ERI_F3z_S_Fx2y_S_a;
  Double I_ERI_F3x_Px_Fxyz_S_a = I_ERI_G4x_S_Fxyz_S_a+ABX*I_ERI_F3x_S_Fxyz_S_a;
  Double I_ERI_F2xy_Px_Fxyz_S_a = I_ERI_G3xy_S_Fxyz_S_a+ABX*I_ERI_F2xy_S_Fxyz_S_a;
  Double I_ERI_F2xz_Px_Fxyz_S_a = I_ERI_G3xz_S_Fxyz_S_a+ABX*I_ERI_F2xz_S_Fxyz_S_a;
  Double I_ERI_Fx2y_Px_Fxyz_S_a = I_ERI_G2x2y_S_Fxyz_S_a+ABX*I_ERI_Fx2y_S_Fxyz_S_a;
  Double I_ERI_Fxyz_Px_Fxyz_S_a = I_ERI_G2xyz_S_Fxyz_S_a+ABX*I_ERI_Fxyz_S_Fxyz_S_a;
  Double I_ERI_Fx2z_Px_Fxyz_S_a = I_ERI_G2x2z_S_Fxyz_S_a+ABX*I_ERI_Fx2z_S_Fxyz_S_a;
  Double I_ERI_F3y_Px_Fxyz_S_a = I_ERI_Gx3y_S_Fxyz_S_a+ABX*I_ERI_F3y_S_Fxyz_S_a;
  Double I_ERI_F2yz_Px_Fxyz_S_a = I_ERI_Gx2yz_S_Fxyz_S_a+ABX*I_ERI_F2yz_S_Fxyz_S_a;
  Double I_ERI_Fy2z_Px_Fxyz_S_a = I_ERI_Gxy2z_S_Fxyz_S_a+ABX*I_ERI_Fy2z_S_Fxyz_S_a;
  Double I_ERI_F3z_Px_Fxyz_S_a = I_ERI_Gx3z_S_Fxyz_S_a+ABX*I_ERI_F3z_S_Fxyz_S_a;
  Double I_ERI_F3x_Py_Fxyz_S_a = I_ERI_G3xy_S_Fxyz_S_a+ABY*I_ERI_F3x_S_Fxyz_S_a;
  Double I_ERI_F2xy_Py_Fxyz_S_a = I_ERI_G2x2y_S_Fxyz_S_a+ABY*I_ERI_F2xy_S_Fxyz_S_a;
  Double I_ERI_F2xz_Py_Fxyz_S_a = I_ERI_G2xyz_S_Fxyz_S_a+ABY*I_ERI_F2xz_S_Fxyz_S_a;
  Double I_ERI_Fx2y_Py_Fxyz_S_a = I_ERI_Gx3y_S_Fxyz_S_a+ABY*I_ERI_Fx2y_S_Fxyz_S_a;
  Double I_ERI_Fxyz_Py_Fxyz_S_a = I_ERI_Gx2yz_S_Fxyz_S_a+ABY*I_ERI_Fxyz_S_Fxyz_S_a;
  Double I_ERI_Fx2z_Py_Fxyz_S_a = I_ERI_Gxy2z_S_Fxyz_S_a+ABY*I_ERI_Fx2z_S_Fxyz_S_a;
  Double I_ERI_F3y_Py_Fxyz_S_a = I_ERI_G4y_S_Fxyz_S_a+ABY*I_ERI_F3y_S_Fxyz_S_a;
  Double I_ERI_F2yz_Py_Fxyz_S_a = I_ERI_G3yz_S_Fxyz_S_a+ABY*I_ERI_F2yz_S_Fxyz_S_a;
  Double I_ERI_Fy2z_Py_Fxyz_S_a = I_ERI_G2y2z_S_Fxyz_S_a+ABY*I_ERI_Fy2z_S_Fxyz_S_a;
  Double I_ERI_F3z_Py_Fxyz_S_a = I_ERI_Gy3z_S_Fxyz_S_a+ABY*I_ERI_F3z_S_Fxyz_S_a;
  Double I_ERI_F3x_Pz_Fxyz_S_a = I_ERI_G3xz_S_Fxyz_S_a+ABZ*I_ERI_F3x_S_Fxyz_S_a;
  Double I_ERI_F2xy_Pz_Fxyz_S_a = I_ERI_G2xyz_S_Fxyz_S_a+ABZ*I_ERI_F2xy_S_Fxyz_S_a;
  Double I_ERI_F2xz_Pz_Fxyz_S_a = I_ERI_G2x2z_S_Fxyz_S_a+ABZ*I_ERI_F2xz_S_Fxyz_S_a;
  Double I_ERI_Fx2y_Pz_Fxyz_S_a = I_ERI_Gx2yz_S_Fxyz_S_a+ABZ*I_ERI_Fx2y_S_Fxyz_S_a;
  Double I_ERI_Fxyz_Pz_Fxyz_S_a = I_ERI_Gxy2z_S_Fxyz_S_a+ABZ*I_ERI_Fxyz_S_Fxyz_S_a;
  Double I_ERI_Fx2z_Pz_Fxyz_S_a = I_ERI_Gx3z_S_Fxyz_S_a+ABZ*I_ERI_Fx2z_S_Fxyz_S_a;
  Double I_ERI_F3y_Pz_Fxyz_S_a = I_ERI_G3yz_S_Fxyz_S_a+ABZ*I_ERI_F3y_S_Fxyz_S_a;
  Double I_ERI_F2yz_Pz_Fxyz_S_a = I_ERI_G2y2z_S_Fxyz_S_a+ABZ*I_ERI_F2yz_S_Fxyz_S_a;
  Double I_ERI_Fy2z_Pz_Fxyz_S_a = I_ERI_Gy3z_S_Fxyz_S_a+ABZ*I_ERI_Fy2z_S_Fxyz_S_a;
  Double I_ERI_F3z_Pz_Fxyz_S_a = I_ERI_G4z_S_Fxyz_S_a+ABZ*I_ERI_F3z_S_Fxyz_S_a;
  Double I_ERI_F3x_Px_Fx2z_S_a = I_ERI_G4x_S_Fx2z_S_a+ABX*I_ERI_F3x_S_Fx2z_S_a;
  Double I_ERI_F2xy_Px_Fx2z_S_a = I_ERI_G3xy_S_Fx2z_S_a+ABX*I_ERI_F2xy_S_Fx2z_S_a;
  Double I_ERI_F2xz_Px_Fx2z_S_a = I_ERI_G3xz_S_Fx2z_S_a+ABX*I_ERI_F2xz_S_Fx2z_S_a;
  Double I_ERI_Fx2y_Px_Fx2z_S_a = I_ERI_G2x2y_S_Fx2z_S_a+ABX*I_ERI_Fx2y_S_Fx2z_S_a;
  Double I_ERI_Fxyz_Px_Fx2z_S_a = I_ERI_G2xyz_S_Fx2z_S_a+ABX*I_ERI_Fxyz_S_Fx2z_S_a;
  Double I_ERI_Fx2z_Px_Fx2z_S_a = I_ERI_G2x2z_S_Fx2z_S_a+ABX*I_ERI_Fx2z_S_Fx2z_S_a;
  Double I_ERI_F3y_Px_Fx2z_S_a = I_ERI_Gx3y_S_Fx2z_S_a+ABX*I_ERI_F3y_S_Fx2z_S_a;
  Double I_ERI_F2yz_Px_Fx2z_S_a = I_ERI_Gx2yz_S_Fx2z_S_a+ABX*I_ERI_F2yz_S_Fx2z_S_a;
  Double I_ERI_Fy2z_Px_Fx2z_S_a = I_ERI_Gxy2z_S_Fx2z_S_a+ABX*I_ERI_Fy2z_S_Fx2z_S_a;
  Double I_ERI_F3z_Px_Fx2z_S_a = I_ERI_Gx3z_S_Fx2z_S_a+ABX*I_ERI_F3z_S_Fx2z_S_a;
  Double I_ERI_F3x_Py_Fx2z_S_a = I_ERI_G3xy_S_Fx2z_S_a+ABY*I_ERI_F3x_S_Fx2z_S_a;
  Double I_ERI_F2xy_Py_Fx2z_S_a = I_ERI_G2x2y_S_Fx2z_S_a+ABY*I_ERI_F2xy_S_Fx2z_S_a;
  Double I_ERI_F2xz_Py_Fx2z_S_a = I_ERI_G2xyz_S_Fx2z_S_a+ABY*I_ERI_F2xz_S_Fx2z_S_a;
  Double I_ERI_Fx2y_Py_Fx2z_S_a = I_ERI_Gx3y_S_Fx2z_S_a+ABY*I_ERI_Fx2y_S_Fx2z_S_a;
  Double I_ERI_Fxyz_Py_Fx2z_S_a = I_ERI_Gx2yz_S_Fx2z_S_a+ABY*I_ERI_Fxyz_S_Fx2z_S_a;
  Double I_ERI_Fx2z_Py_Fx2z_S_a = I_ERI_Gxy2z_S_Fx2z_S_a+ABY*I_ERI_Fx2z_S_Fx2z_S_a;
  Double I_ERI_F3y_Py_Fx2z_S_a = I_ERI_G4y_S_Fx2z_S_a+ABY*I_ERI_F3y_S_Fx2z_S_a;
  Double I_ERI_F2yz_Py_Fx2z_S_a = I_ERI_G3yz_S_Fx2z_S_a+ABY*I_ERI_F2yz_S_Fx2z_S_a;
  Double I_ERI_Fy2z_Py_Fx2z_S_a = I_ERI_G2y2z_S_Fx2z_S_a+ABY*I_ERI_Fy2z_S_Fx2z_S_a;
  Double I_ERI_F3z_Py_Fx2z_S_a = I_ERI_Gy3z_S_Fx2z_S_a+ABY*I_ERI_F3z_S_Fx2z_S_a;
  Double I_ERI_F3x_Pz_Fx2z_S_a = I_ERI_G3xz_S_Fx2z_S_a+ABZ*I_ERI_F3x_S_Fx2z_S_a;
  Double I_ERI_F2xy_Pz_Fx2z_S_a = I_ERI_G2xyz_S_Fx2z_S_a+ABZ*I_ERI_F2xy_S_Fx2z_S_a;
  Double I_ERI_F2xz_Pz_Fx2z_S_a = I_ERI_G2x2z_S_Fx2z_S_a+ABZ*I_ERI_F2xz_S_Fx2z_S_a;
  Double I_ERI_Fx2y_Pz_Fx2z_S_a = I_ERI_Gx2yz_S_Fx2z_S_a+ABZ*I_ERI_Fx2y_S_Fx2z_S_a;
  Double I_ERI_Fxyz_Pz_Fx2z_S_a = I_ERI_Gxy2z_S_Fx2z_S_a+ABZ*I_ERI_Fxyz_S_Fx2z_S_a;
  Double I_ERI_Fx2z_Pz_Fx2z_S_a = I_ERI_Gx3z_S_Fx2z_S_a+ABZ*I_ERI_Fx2z_S_Fx2z_S_a;
  Double I_ERI_F3y_Pz_Fx2z_S_a = I_ERI_G3yz_S_Fx2z_S_a+ABZ*I_ERI_F3y_S_Fx2z_S_a;
  Double I_ERI_F2yz_Pz_Fx2z_S_a = I_ERI_G2y2z_S_Fx2z_S_a+ABZ*I_ERI_F2yz_S_Fx2z_S_a;
  Double I_ERI_Fy2z_Pz_Fx2z_S_a = I_ERI_Gy3z_S_Fx2z_S_a+ABZ*I_ERI_Fy2z_S_Fx2z_S_a;
  Double I_ERI_F3z_Pz_Fx2z_S_a = I_ERI_G4z_S_Fx2z_S_a+ABZ*I_ERI_F3z_S_Fx2z_S_a;
  Double I_ERI_F3x_Px_F3y_S_a = I_ERI_G4x_S_F3y_S_a+ABX*I_ERI_F3x_S_F3y_S_a;
  Double I_ERI_F2xy_Px_F3y_S_a = I_ERI_G3xy_S_F3y_S_a+ABX*I_ERI_F2xy_S_F3y_S_a;
  Double I_ERI_F2xz_Px_F3y_S_a = I_ERI_G3xz_S_F3y_S_a+ABX*I_ERI_F2xz_S_F3y_S_a;
  Double I_ERI_Fx2y_Px_F3y_S_a = I_ERI_G2x2y_S_F3y_S_a+ABX*I_ERI_Fx2y_S_F3y_S_a;
  Double I_ERI_Fxyz_Px_F3y_S_a = I_ERI_G2xyz_S_F3y_S_a+ABX*I_ERI_Fxyz_S_F3y_S_a;
  Double I_ERI_Fx2z_Px_F3y_S_a = I_ERI_G2x2z_S_F3y_S_a+ABX*I_ERI_Fx2z_S_F3y_S_a;
  Double I_ERI_F3y_Px_F3y_S_a = I_ERI_Gx3y_S_F3y_S_a+ABX*I_ERI_F3y_S_F3y_S_a;
  Double I_ERI_F2yz_Px_F3y_S_a = I_ERI_Gx2yz_S_F3y_S_a+ABX*I_ERI_F2yz_S_F3y_S_a;
  Double I_ERI_Fy2z_Px_F3y_S_a = I_ERI_Gxy2z_S_F3y_S_a+ABX*I_ERI_Fy2z_S_F3y_S_a;
  Double I_ERI_F3z_Px_F3y_S_a = I_ERI_Gx3z_S_F3y_S_a+ABX*I_ERI_F3z_S_F3y_S_a;
  Double I_ERI_F3x_Py_F3y_S_a = I_ERI_G3xy_S_F3y_S_a+ABY*I_ERI_F3x_S_F3y_S_a;
  Double I_ERI_F2xy_Py_F3y_S_a = I_ERI_G2x2y_S_F3y_S_a+ABY*I_ERI_F2xy_S_F3y_S_a;
  Double I_ERI_F2xz_Py_F3y_S_a = I_ERI_G2xyz_S_F3y_S_a+ABY*I_ERI_F2xz_S_F3y_S_a;
  Double I_ERI_Fx2y_Py_F3y_S_a = I_ERI_Gx3y_S_F3y_S_a+ABY*I_ERI_Fx2y_S_F3y_S_a;
  Double I_ERI_Fxyz_Py_F3y_S_a = I_ERI_Gx2yz_S_F3y_S_a+ABY*I_ERI_Fxyz_S_F3y_S_a;
  Double I_ERI_Fx2z_Py_F3y_S_a = I_ERI_Gxy2z_S_F3y_S_a+ABY*I_ERI_Fx2z_S_F3y_S_a;
  Double I_ERI_F3y_Py_F3y_S_a = I_ERI_G4y_S_F3y_S_a+ABY*I_ERI_F3y_S_F3y_S_a;
  Double I_ERI_F2yz_Py_F3y_S_a = I_ERI_G3yz_S_F3y_S_a+ABY*I_ERI_F2yz_S_F3y_S_a;
  Double I_ERI_Fy2z_Py_F3y_S_a = I_ERI_G2y2z_S_F3y_S_a+ABY*I_ERI_Fy2z_S_F3y_S_a;
  Double I_ERI_F3z_Py_F3y_S_a = I_ERI_Gy3z_S_F3y_S_a+ABY*I_ERI_F3z_S_F3y_S_a;
  Double I_ERI_F3x_Pz_F3y_S_a = I_ERI_G3xz_S_F3y_S_a+ABZ*I_ERI_F3x_S_F3y_S_a;
  Double I_ERI_F2xy_Pz_F3y_S_a = I_ERI_G2xyz_S_F3y_S_a+ABZ*I_ERI_F2xy_S_F3y_S_a;
  Double I_ERI_F2xz_Pz_F3y_S_a = I_ERI_G2x2z_S_F3y_S_a+ABZ*I_ERI_F2xz_S_F3y_S_a;
  Double I_ERI_Fx2y_Pz_F3y_S_a = I_ERI_Gx2yz_S_F3y_S_a+ABZ*I_ERI_Fx2y_S_F3y_S_a;
  Double I_ERI_Fxyz_Pz_F3y_S_a = I_ERI_Gxy2z_S_F3y_S_a+ABZ*I_ERI_Fxyz_S_F3y_S_a;
  Double I_ERI_Fx2z_Pz_F3y_S_a = I_ERI_Gx3z_S_F3y_S_a+ABZ*I_ERI_Fx2z_S_F3y_S_a;
  Double I_ERI_F3y_Pz_F3y_S_a = I_ERI_G3yz_S_F3y_S_a+ABZ*I_ERI_F3y_S_F3y_S_a;
  Double I_ERI_F2yz_Pz_F3y_S_a = I_ERI_G2y2z_S_F3y_S_a+ABZ*I_ERI_F2yz_S_F3y_S_a;
  Double I_ERI_Fy2z_Pz_F3y_S_a = I_ERI_Gy3z_S_F3y_S_a+ABZ*I_ERI_Fy2z_S_F3y_S_a;
  Double I_ERI_F3z_Pz_F3y_S_a = I_ERI_G4z_S_F3y_S_a+ABZ*I_ERI_F3z_S_F3y_S_a;
  Double I_ERI_F3x_Px_F2yz_S_a = I_ERI_G4x_S_F2yz_S_a+ABX*I_ERI_F3x_S_F2yz_S_a;
  Double I_ERI_F2xy_Px_F2yz_S_a = I_ERI_G3xy_S_F2yz_S_a+ABX*I_ERI_F2xy_S_F2yz_S_a;
  Double I_ERI_F2xz_Px_F2yz_S_a = I_ERI_G3xz_S_F2yz_S_a+ABX*I_ERI_F2xz_S_F2yz_S_a;
  Double I_ERI_Fx2y_Px_F2yz_S_a = I_ERI_G2x2y_S_F2yz_S_a+ABX*I_ERI_Fx2y_S_F2yz_S_a;
  Double I_ERI_Fxyz_Px_F2yz_S_a = I_ERI_G2xyz_S_F2yz_S_a+ABX*I_ERI_Fxyz_S_F2yz_S_a;
  Double I_ERI_Fx2z_Px_F2yz_S_a = I_ERI_G2x2z_S_F2yz_S_a+ABX*I_ERI_Fx2z_S_F2yz_S_a;
  Double I_ERI_F3y_Px_F2yz_S_a = I_ERI_Gx3y_S_F2yz_S_a+ABX*I_ERI_F3y_S_F2yz_S_a;
  Double I_ERI_F2yz_Px_F2yz_S_a = I_ERI_Gx2yz_S_F2yz_S_a+ABX*I_ERI_F2yz_S_F2yz_S_a;
  Double I_ERI_Fy2z_Px_F2yz_S_a = I_ERI_Gxy2z_S_F2yz_S_a+ABX*I_ERI_Fy2z_S_F2yz_S_a;
  Double I_ERI_F3z_Px_F2yz_S_a = I_ERI_Gx3z_S_F2yz_S_a+ABX*I_ERI_F3z_S_F2yz_S_a;
  Double I_ERI_F3x_Py_F2yz_S_a = I_ERI_G3xy_S_F2yz_S_a+ABY*I_ERI_F3x_S_F2yz_S_a;
  Double I_ERI_F2xy_Py_F2yz_S_a = I_ERI_G2x2y_S_F2yz_S_a+ABY*I_ERI_F2xy_S_F2yz_S_a;
  Double I_ERI_F2xz_Py_F2yz_S_a = I_ERI_G2xyz_S_F2yz_S_a+ABY*I_ERI_F2xz_S_F2yz_S_a;
  Double I_ERI_Fx2y_Py_F2yz_S_a = I_ERI_Gx3y_S_F2yz_S_a+ABY*I_ERI_Fx2y_S_F2yz_S_a;
  Double I_ERI_Fxyz_Py_F2yz_S_a = I_ERI_Gx2yz_S_F2yz_S_a+ABY*I_ERI_Fxyz_S_F2yz_S_a;
  Double I_ERI_Fx2z_Py_F2yz_S_a = I_ERI_Gxy2z_S_F2yz_S_a+ABY*I_ERI_Fx2z_S_F2yz_S_a;
  Double I_ERI_F3y_Py_F2yz_S_a = I_ERI_G4y_S_F2yz_S_a+ABY*I_ERI_F3y_S_F2yz_S_a;
  Double I_ERI_F2yz_Py_F2yz_S_a = I_ERI_G3yz_S_F2yz_S_a+ABY*I_ERI_F2yz_S_F2yz_S_a;
  Double I_ERI_Fy2z_Py_F2yz_S_a = I_ERI_G2y2z_S_F2yz_S_a+ABY*I_ERI_Fy2z_S_F2yz_S_a;
  Double I_ERI_F3z_Py_F2yz_S_a = I_ERI_Gy3z_S_F2yz_S_a+ABY*I_ERI_F3z_S_F2yz_S_a;
  Double I_ERI_F3x_Pz_F2yz_S_a = I_ERI_G3xz_S_F2yz_S_a+ABZ*I_ERI_F3x_S_F2yz_S_a;
  Double I_ERI_F2xy_Pz_F2yz_S_a = I_ERI_G2xyz_S_F2yz_S_a+ABZ*I_ERI_F2xy_S_F2yz_S_a;
  Double I_ERI_F2xz_Pz_F2yz_S_a = I_ERI_G2x2z_S_F2yz_S_a+ABZ*I_ERI_F2xz_S_F2yz_S_a;
  Double I_ERI_Fx2y_Pz_F2yz_S_a = I_ERI_Gx2yz_S_F2yz_S_a+ABZ*I_ERI_Fx2y_S_F2yz_S_a;
  Double I_ERI_Fxyz_Pz_F2yz_S_a = I_ERI_Gxy2z_S_F2yz_S_a+ABZ*I_ERI_Fxyz_S_F2yz_S_a;
  Double I_ERI_Fx2z_Pz_F2yz_S_a = I_ERI_Gx3z_S_F2yz_S_a+ABZ*I_ERI_Fx2z_S_F2yz_S_a;
  Double I_ERI_F3y_Pz_F2yz_S_a = I_ERI_G3yz_S_F2yz_S_a+ABZ*I_ERI_F3y_S_F2yz_S_a;
  Double I_ERI_F2yz_Pz_F2yz_S_a = I_ERI_G2y2z_S_F2yz_S_a+ABZ*I_ERI_F2yz_S_F2yz_S_a;
  Double I_ERI_Fy2z_Pz_F2yz_S_a = I_ERI_Gy3z_S_F2yz_S_a+ABZ*I_ERI_Fy2z_S_F2yz_S_a;
  Double I_ERI_F3z_Pz_F2yz_S_a = I_ERI_G4z_S_F2yz_S_a+ABZ*I_ERI_F3z_S_F2yz_S_a;
  Double I_ERI_F3x_Px_Fy2z_S_a = I_ERI_G4x_S_Fy2z_S_a+ABX*I_ERI_F3x_S_Fy2z_S_a;
  Double I_ERI_F2xy_Px_Fy2z_S_a = I_ERI_G3xy_S_Fy2z_S_a+ABX*I_ERI_F2xy_S_Fy2z_S_a;
  Double I_ERI_F2xz_Px_Fy2z_S_a = I_ERI_G3xz_S_Fy2z_S_a+ABX*I_ERI_F2xz_S_Fy2z_S_a;
  Double I_ERI_Fx2y_Px_Fy2z_S_a = I_ERI_G2x2y_S_Fy2z_S_a+ABX*I_ERI_Fx2y_S_Fy2z_S_a;
  Double I_ERI_Fxyz_Px_Fy2z_S_a = I_ERI_G2xyz_S_Fy2z_S_a+ABX*I_ERI_Fxyz_S_Fy2z_S_a;
  Double I_ERI_Fx2z_Px_Fy2z_S_a = I_ERI_G2x2z_S_Fy2z_S_a+ABX*I_ERI_Fx2z_S_Fy2z_S_a;
  Double I_ERI_F3y_Px_Fy2z_S_a = I_ERI_Gx3y_S_Fy2z_S_a+ABX*I_ERI_F3y_S_Fy2z_S_a;
  Double I_ERI_F2yz_Px_Fy2z_S_a = I_ERI_Gx2yz_S_Fy2z_S_a+ABX*I_ERI_F2yz_S_Fy2z_S_a;
  Double I_ERI_Fy2z_Px_Fy2z_S_a = I_ERI_Gxy2z_S_Fy2z_S_a+ABX*I_ERI_Fy2z_S_Fy2z_S_a;
  Double I_ERI_F3z_Px_Fy2z_S_a = I_ERI_Gx3z_S_Fy2z_S_a+ABX*I_ERI_F3z_S_Fy2z_S_a;
  Double I_ERI_F3x_Py_Fy2z_S_a = I_ERI_G3xy_S_Fy2z_S_a+ABY*I_ERI_F3x_S_Fy2z_S_a;
  Double I_ERI_F2xy_Py_Fy2z_S_a = I_ERI_G2x2y_S_Fy2z_S_a+ABY*I_ERI_F2xy_S_Fy2z_S_a;
  Double I_ERI_F2xz_Py_Fy2z_S_a = I_ERI_G2xyz_S_Fy2z_S_a+ABY*I_ERI_F2xz_S_Fy2z_S_a;
  Double I_ERI_Fx2y_Py_Fy2z_S_a = I_ERI_Gx3y_S_Fy2z_S_a+ABY*I_ERI_Fx2y_S_Fy2z_S_a;
  Double I_ERI_Fxyz_Py_Fy2z_S_a = I_ERI_Gx2yz_S_Fy2z_S_a+ABY*I_ERI_Fxyz_S_Fy2z_S_a;
  Double I_ERI_Fx2z_Py_Fy2z_S_a = I_ERI_Gxy2z_S_Fy2z_S_a+ABY*I_ERI_Fx2z_S_Fy2z_S_a;
  Double I_ERI_F3y_Py_Fy2z_S_a = I_ERI_G4y_S_Fy2z_S_a+ABY*I_ERI_F3y_S_Fy2z_S_a;
  Double I_ERI_F2yz_Py_Fy2z_S_a = I_ERI_G3yz_S_Fy2z_S_a+ABY*I_ERI_F2yz_S_Fy2z_S_a;
  Double I_ERI_Fy2z_Py_Fy2z_S_a = I_ERI_G2y2z_S_Fy2z_S_a+ABY*I_ERI_Fy2z_S_Fy2z_S_a;
  Double I_ERI_F3z_Py_Fy2z_S_a = I_ERI_Gy3z_S_Fy2z_S_a+ABY*I_ERI_F3z_S_Fy2z_S_a;
  Double I_ERI_F3x_Pz_Fy2z_S_a = I_ERI_G3xz_S_Fy2z_S_a+ABZ*I_ERI_F3x_S_Fy2z_S_a;
  Double I_ERI_F2xy_Pz_Fy2z_S_a = I_ERI_G2xyz_S_Fy2z_S_a+ABZ*I_ERI_F2xy_S_Fy2z_S_a;
  Double I_ERI_F2xz_Pz_Fy2z_S_a = I_ERI_G2x2z_S_Fy2z_S_a+ABZ*I_ERI_F2xz_S_Fy2z_S_a;
  Double I_ERI_Fx2y_Pz_Fy2z_S_a = I_ERI_Gx2yz_S_Fy2z_S_a+ABZ*I_ERI_Fx2y_S_Fy2z_S_a;
  Double I_ERI_Fxyz_Pz_Fy2z_S_a = I_ERI_Gxy2z_S_Fy2z_S_a+ABZ*I_ERI_Fxyz_S_Fy2z_S_a;
  Double I_ERI_Fx2z_Pz_Fy2z_S_a = I_ERI_Gx3z_S_Fy2z_S_a+ABZ*I_ERI_Fx2z_S_Fy2z_S_a;
  Double I_ERI_F3y_Pz_Fy2z_S_a = I_ERI_G3yz_S_Fy2z_S_a+ABZ*I_ERI_F3y_S_Fy2z_S_a;
  Double I_ERI_F2yz_Pz_Fy2z_S_a = I_ERI_G2y2z_S_Fy2z_S_a+ABZ*I_ERI_F2yz_S_Fy2z_S_a;
  Double I_ERI_Fy2z_Pz_Fy2z_S_a = I_ERI_Gy3z_S_Fy2z_S_a+ABZ*I_ERI_Fy2z_S_Fy2z_S_a;
  Double I_ERI_F3z_Pz_Fy2z_S_a = I_ERI_G4z_S_Fy2z_S_a+ABZ*I_ERI_F3z_S_Fy2z_S_a;
  Double I_ERI_F3x_Px_F3z_S_a = I_ERI_G4x_S_F3z_S_a+ABX*I_ERI_F3x_S_F3z_S_a;
  Double I_ERI_F2xy_Px_F3z_S_a = I_ERI_G3xy_S_F3z_S_a+ABX*I_ERI_F2xy_S_F3z_S_a;
  Double I_ERI_F2xz_Px_F3z_S_a = I_ERI_G3xz_S_F3z_S_a+ABX*I_ERI_F2xz_S_F3z_S_a;
  Double I_ERI_Fx2y_Px_F3z_S_a = I_ERI_G2x2y_S_F3z_S_a+ABX*I_ERI_Fx2y_S_F3z_S_a;
  Double I_ERI_Fxyz_Px_F3z_S_a = I_ERI_G2xyz_S_F3z_S_a+ABX*I_ERI_Fxyz_S_F3z_S_a;
  Double I_ERI_Fx2z_Px_F3z_S_a = I_ERI_G2x2z_S_F3z_S_a+ABX*I_ERI_Fx2z_S_F3z_S_a;
  Double I_ERI_F3y_Px_F3z_S_a = I_ERI_Gx3y_S_F3z_S_a+ABX*I_ERI_F3y_S_F3z_S_a;
  Double I_ERI_F2yz_Px_F3z_S_a = I_ERI_Gx2yz_S_F3z_S_a+ABX*I_ERI_F2yz_S_F3z_S_a;
  Double I_ERI_Fy2z_Px_F3z_S_a = I_ERI_Gxy2z_S_F3z_S_a+ABX*I_ERI_Fy2z_S_F3z_S_a;
  Double I_ERI_F3z_Px_F3z_S_a = I_ERI_Gx3z_S_F3z_S_a+ABX*I_ERI_F3z_S_F3z_S_a;
  Double I_ERI_F3x_Py_F3z_S_a = I_ERI_G3xy_S_F3z_S_a+ABY*I_ERI_F3x_S_F3z_S_a;
  Double I_ERI_F2xy_Py_F3z_S_a = I_ERI_G2x2y_S_F3z_S_a+ABY*I_ERI_F2xy_S_F3z_S_a;
  Double I_ERI_F2xz_Py_F3z_S_a = I_ERI_G2xyz_S_F3z_S_a+ABY*I_ERI_F2xz_S_F3z_S_a;
  Double I_ERI_Fx2y_Py_F3z_S_a = I_ERI_Gx3y_S_F3z_S_a+ABY*I_ERI_Fx2y_S_F3z_S_a;
  Double I_ERI_Fxyz_Py_F3z_S_a = I_ERI_Gx2yz_S_F3z_S_a+ABY*I_ERI_Fxyz_S_F3z_S_a;
  Double I_ERI_Fx2z_Py_F3z_S_a = I_ERI_Gxy2z_S_F3z_S_a+ABY*I_ERI_Fx2z_S_F3z_S_a;
  Double I_ERI_F3y_Py_F3z_S_a = I_ERI_G4y_S_F3z_S_a+ABY*I_ERI_F3y_S_F3z_S_a;
  Double I_ERI_F2yz_Py_F3z_S_a = I_ERI_G3yz_S_F3z_S_a+ABY*I_ERI_F2yz_S_F3z_S_a;
  Double I_ERI_Fy2z_Py_F3z_S_a = I_ERI_G2y2z_S_F3z_S_a+ABY*I_ERI_Fy2z_S_F3z_S_a;
  Double I_ERI_F3z_Py_F3z_S_a = I_ERI_Gy3z_S_F3z_S_a+ABY*I_ERI_F3z_S_F3z_S_a;
  Double I_ERI_F3x_Pz_F3z_S_a = I_ERI_G3xz_S_F3z_S_a+ABZ*I_ERI_F3x_S_F3z_S_a;
  Double I_ERI_F2xy_Pz_F3z_S_a = I_ERI_G2xyz_S_F3z_S_a+ABZ*I_ERI_F2xy_S_F3z_S_a;
  Double I_ERI_F2xz_Pz_F3z_S_a = I_ERI_G2x2z_S_F3z_S_a+ABZ*I_ERI_F2xz_S_F3z_S_a;
  Double I_ERI_Fx2y_Pz_F3z_S_a = I_ERI_Gx2yz_S_F3z_S_a+ABZ*I_ERI_Fx2y_S_F3z_S_a;
  Double I_ERI_Fxyz_Pz_F3z_S_a = I_ERI_Gxy2z_S_F3z_S_a+ABZ*I_ERI_Fxyz_S_F3z_S_a;
  Double I_ERI_Fx2z_Pz_F3z_S_a = I_ERI_Gx3z_S_F3z_S_a+ABZ*I_ERI_Fx2z_S_F3z_S_a;
  Double I_ERI_F3y_Pz_F3z_S_a = I_ERI_G3yz_S_F3z_S_a+ABZ*I_ERI_F3y_S_F3z_S_a;
  Double I_ERI_F2yz_Pz_F3z_S_a = I_ERI_G2y2z_S_F3z_S_a+ABZ*I_ERI_F2yz_S_F3z_S_a;
  Double I_ERI_Fy2z_Pz_F3z_S_a = I_ERI_Gy3z_S_F3z_S_a+ABZ*I_ERI_Fy2z_S_F3z_S_a;
  Double I_ERI_F3z_Pz_F3z_S_a = I_ERI_G4z_S_F3z_S_a+ABZ*I_ERI_F3z_S_F3z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_b
   * RHS shell quartet name: SQ_ERI_D_S_F_S_b
   ************************************************************/
  Double I_ERI_D2x_Px_F3x_S_b = I_ERI_F3x_S_F3x_S_b+ABX*I_ERI_D2x_S_F3x_S_b;
  Double I_ERI_Dxy_Px_F3x_S_b = I_ERI_F2xy_S_F3x_S_b+ABX*I_ERI_Dxy_S_F3x_S_b;
  Double I_ERI_Dxz_Px_F3x_S_b = I_ERI_F2xz_S_F3x_S_b+ABX*I_ERI_Dxz_S_F3x_S_b;
  Double I_ERI_D2y_Px_F3x_S_b = I_ERI_Fx2y_S_F3x_S_b+ABX*I_ERI_D2y_S_F3x_S_b;
  Double I_ERI_Dyz_Px_F3x_S_b = I_ERI_Fxyz_S_F3x_S_b+ABX*I_ERI_Dyz_S_F3x_S_b;
  Double I_ERI_D2z_Px_F3x_S_b = I_ERI_Fx2z_S_F3x_S_b+ABX*I_ERI_D2z_S_F3x_S_b;
  Double I_ERI_D2x_Py_F3x_S_b = I_ERI_F2xy_S_F3x_S_b+ABY*I_ERI_D2x_S_F3x_S_b;
  Double I_ERI_Dxy_Py_F3x_S_b = I_ERI_Fx2y_S_F3x_S_b+ABY*I_ERI_Dxy_S_F3x_S_b;
  Double I_ERI_Dxz_Py_F3x_S_b = I_ERI_Fxyz_S_F3x_S_b+ABY*I_ERI_Dxz_S_F3x_S_b;
  Double I_ERI_D2y_Py_F3x_S_b = I_ERI_F3y_S_F3x_S_b+ABY*I_ERI_D2y_S_F3x_S_b;
  Double I_ERI_Dyz_Py_F3x_S_b = I_ERI_F2yz_S_F3x_S_b+ABY*I_ERI_Dyz_S_F3x_S_b;
  Double I_ERI_D2z_Py_F3x_S_b = I_ERI_Fy2z_S_F3x_S_b+ABY*I_ERI_D2z_S_F3x_S_b;
  Double I_ERI_D2x_Pz_F3x_S_b = I_ERI_F2xz_S_F3x_S_b+ABZ*I_ERI_D2x_S_F3x_S_b;
  Double I_ERI_Dxy_Pz_F3x_S_b = I_ERI_Fxyz_S_F3x_S_b+ABZ*I_ERI_Dxy_S_F3x_S_b;
  Double I_ERI_Dxz_Pz_F3x_S_b = I_ERI_Fx2z_S_F3x_S_b+ABZ*I_ERI_Dxz_S_F3x_S_b;
  Double I_ERI_D2y_Pz_F3x_S_b = I_ERI_F2yz_S_F3x_S_b+ABZ*I_ERI_D2y_S_F3x_S_b;
  Double I_ERI_Dyz_Pz_F3x_S_b = I_ERI_Fy2z_S_F3x_S_b+ABZ*I_ERI_Dyz_S_F3x_S_b;
  Double I_ERI_D2z_Pz_F3x_S_b = I_ERI_F3z_S_F3x_S_b+ABZ*I_ERI_D2z_S_F3x_S_b;
  Double I_ERI_D2x_Px_F2xy_S_b = I_ERI_F3x_S_F2xy_S_b+ABX*I_ERI_D2x_S_F2xy_S_b;
  Double I_ERI_Dxy_Px_F2xy_S_b = I_ERI_F2xy_S_F2xy_S_b+ABX*I_ERI_Dxy_S_F2xy_S_b;
  Double I_ERI_Dxz_Px_F2xy_S_b = I_ERI_F2xz_S_F2xy_S_b+ABX*I_ERI_Dxz_S_F2xy_S_b;
  Double I_ERI_D2y_Px_F2xy_S_b = I_ERI_Fx2y_S_F2xy_S_b+ABX*I_ERI_D2y_S_F2xy_S_b;
  Double I_ERI_Dyz_Px_F2xy_S_b = I_ERI_Fxyz_S_F2xy_S_b+ABX*I_ERI_Dyz_S_F2xy_S_b;
  Double I_ERI_D2z_Px_F2xy_S_b = I_ERI_Fx2z_S_F2xy_S_b+ABX*I_ERI_D2z_S_F2xy_S_b;
  Double I_ERI_D2x_Py_F2xy_S_b = I_ERI_F2xy_S_F2xy_S_b+ABY*I_ERI_D2x_S_F2xy_S_b;
  Double I_ERI_Dxy_Py_F2xy_S_b = I_ERI_Fx2y_S_F2xy_S_b+ABY*I_ERI_Dxy_S_F2xy_S_b;
  Double I_ERI_Dxz_Py_F2xy_S_b = I_ERI_Fxyz_S_F2xy_S_b+ABY*I_ERI_Dxz_S_F2xy_S_b;
  Double I_ERI_D2y_Py_F2xy_S_b = I_ERI_F3y_S_F2xy_S_b+ABY*I_ERI_D2y_S_F2xy_S_b;
  Double I_ERI_Dyz_Py_F2xy_S_b = I_ERI_F2yz_S_F2xy_S_b+ABY*I_ERI_Dyz_S_F2xy_S_b;
  Double I_ERI_D2z_Py_F2xy_S_b = I_ERI_Fy2z_S_F2xy_S_b+ABY*I_ERI_D2z_S_F2xy_S_b;
  Double I_ERI_D2x_Pz_F2xy_S_b = I_ERI_F2xz_S_F2xy_S_b+ABZ*I_ERI_D2x_S_F2xy_S_b;
  Double I_ERI_Dxy_Pz_F2xy_S_b = I_ERI_Fxyz_S_F2xy_S_b+ABZ*I_ERI_Dxy_S_F2xy_S_b;
  Double I_ERI_Dxz_Pz_F2xy_S_b = I_ERI_Fx2z_S_F2xy_S_b+ABZ*I_ERI_Dxz_S_F2xy_S_b;
  Double I_ERI_D2y_Pz_F2xy_S_b = I_ERI_F2yz_S_F2xy_S_b+ABZ*I_ERI_D2y_S_F2xy_S_b;
  Double I_ERI_Dyz_Pz_F2xy_S_b = I_ERI_Fy2z_S_F2xy_S_b+ABZ*I_ERI_Dyz_S_F2xy_S_b;
  Double I_ERI_D2z_Pz_F2xy_S_b = I_ERI_F3z_S_F2xy_S_b+ABZ*I_ERI_D2z_S_F2xy_S_b;
  Double I_ERI_D2x_Px_F2xz_S_b = I_ERI_F3x_S_F2xz_S_b+ABX*I_ERI_D2x_S_F2xz_S_b;
  Double I_ERI_Dxy_Px_F2xz_S_b = I_ERI_F2xy_S_F2xz_S_b+ABX*I_ERI_Dxy_S_F2xz_S_b;
  Double I_ERI_Dxz_Px_F2xz_S_b = I_ERI_F2xz_S_F2xz_S_b+ABX*I_ERI_Dxz_S_F2xz_S_b;
  Double I_ERI_D2y_Px_F2xz_S_b = I_ERI_Fx2y_S_F2xz_S_b+ABX*I_ERI_D2y_S_F2xz_S_b;
  Double I_ERI_Dyz_Px_F2xz_S_b = I_ERI_Fxyz_S_F2xz_S_b+ABX*I_ERI_Dyz_S_F2xz_S_b;
  Double I_ERI_D2z_Px_F2xz_S_b = I_ERI_Fx2z_S_F2xz_S_b+ABX*I_ERI_D2z_S_F2xz_S_b;
  Double I_ERI_D2x_Py_F2xz_S_b = I_ERI_F2xy_S_F2xz_S_b+ABY*I_ERI_D2x_S_F2xz_S_b;
  Double I_ERI_Dxy_Py_F2xz_S_b = I_ERI_Fx2y_S_F2xz_S_b+ABY*I_ERI_Dxy_S_F2xz_S_b;
  Double I_ERI_Dxz_Py_F2xz_S_b = I_ERI_Fxyz_S_F2xz_S_b+ABY*I_ERI_Dxz_S_F2xz_S_b;
  Double I_ERI_D2y_Py_F2xz_S_b = I_ERI_F3y_S_F2xz_S_b+ABY*I_ERI_D2y_S_F2xz_S_b;
  Double I_ERI_Dyz_Py_F2xz_S_b = I_ERI_F2yz_S_F2xz_S_b+ABY*I_ERI_Dyz_S_F2xz_S_b;
  Double I_ERI_D2z_Py_F2xz_S_b = I_ERI_Fy2z_S_F2xz_S_b+ABY*I_ERI_D2z_S_F2xz_S_b;
  Double I_ERI_D2x_Pz_F2xz_S_b = I_ERI_F2xz_S_F2xz_S_b+ABZ*I_ERI_D2x_S_F2xz_S_b;
  Double I_ERI_Dxy_Pz_F2xz_S_b = I_ERI_Fxyz_S_F2xz_S_b+ABZ*I_ERI_Dxy_S_F2xz_S_b;
  Double I_ERI_Dxz_Pz_F2xz_S_b = I_ERI_Fx2z_S_F2xz_S_b+ABZ*I_ERI_Dxz_S_F2xz_S_b;
  Double I_ERI_D2y_Pz_F2xz_S_b = I_ERI_F2yz_S_F2xz_S_b+ABZ*I_ERI_D2y_S_F2xz_S_b;
  Double I_ERI_Dyz_Pz_F2xz_S_b = I_ERI_Fy2z_S_F2xz_S_b+ABZ*I_ERI_Dyz_S_F2xz_S_b;
  Double I_ERI_D2z_Pz_F2xz_S_b = I_ERI_F3z_S_F2xz_S_b+ABZ*I_ERI_D2z_S_F2xz_S_b;
  Double I_ERI_D2x_Px_Fx2y_S_b = I_ERI_F3x_S_Fx2y_S_b+ABX*I_ERI_D2x_S_Fx2y_S_b;
  Double I_ERI_Dxy_Px_Fx2y_S_b = I_ERI_F2xy_S_Fx2y_S_b+ABX*I_ERI_Dxy_S_Fx2y_S_b;
  Double I_ERI_Dxz_Px_Fx2y_S_b = I_ERI_F2xz_S_Fx2y_S_b+ABX*I_ERI_Dxz_S_Fx2y_S_b;
  Double I_ERI_D2y_Px_Fx2y_S_b = I_ERI_Fx2y_S_Fx2y_S_b+ABX*I_ERI_D2y_S_Fx2y_S_b;
  Double I_ERI_Dyz_Px_Fx2y_S_b = I_ERI_Fxyz_S_Fx2y_S_b+ABX*I_ERI_Dyz_S_Fx2y_S_b;
  Double I_ERI_D2z_Px_Fx2y_S_b = I_ERI_Fx2z_S_Fx2y_S_b+ABX*I_ERI_D2z_S_Fx2y_S_b;
  Double I_ERI_D2x_Py_Fx2y_S_b = I_ERI_F2xy_S_Fx2y_S_b+ABY*I_ERI_D2x_S_Fx2y_S_b;
  Double I_ERI_Dxy_Py_Fx2y_S_b = I_ERI_Fx2y_S_Fx2y_S_b+ABY*I_ERI_Dxy_S_Fx2y_S_b;
  Double I_ERI_Dxz_Py_Fx2y_S_b = I_ERI_Fxyz_S_Fx2y_S_b+ABY*I_ERI_Dxz_S_Fx2y_S_b;
  Double I_ERI_D2y_Py_Fx2y_S_b = I_ERI_F3y_S_Fx2y_S_b+ABY*I_ERI_D2y_S_Fx2y_S_b;
  Double I_ERI_Dyz_Py_Fx2y_S_b = I_ERI_F2yz_S_Fx2y_S_b+ABY*I_ERI_Dyz_S_Fx2y_S_b;
  Double I_ERI_D2z_Py_Fx2y_S_b = I_ERI_Fy2z_S_Fx2y_S_b+ABY*I_ERI_D2z_S_Fx2y_S_b;
  Double I_ERI_D2x_Pz_Fx2y_S_b = I_ERI_F2xz_S_Fx2y_S_b+ABZ*I_ERI_D2x_S_Fx2y_S_b;
  Double I_ERI_Dxy_Pz_Fx2y_S_b = I_ERI_Fxyz_S_Fx2y_S_b+ABZ*I_ERI_Dxy_S_Fx2y_S_b;
  Double I_ERI_Dxz_Pz_Fx2y_S_b = I_ERI_Fx2z_S_Fx2y_S_b+ABZ*I_ERI_Dxz_S_Fx2y_S_b;
  Double I_ERI_D2y_Pz_Fx2y_S_b = I_ERI_F2yz_S_Fx2y_S_b+ABZ*I_ERI_D2y_S_Fx2y_S_b;
  Double I_ERI_Dyz_Pz_Fx2y_S_b = I_ERI_Fy2z_S_Fx2y_S_b+ABZ*I_ERI_Dyz_S_Fx2y_S_b;
  Double I_ERI_D2z_Pz_Fx2y_S_b = I_ERI_F3z_S_Fx2y_S_b+ABZ*I_ERI_D2z_S_Fx2y_S_b;
  Double I_ERI_D2x_Px_Fxyz_S_b = I_ERI_F3x_S_Fxyz_S_b+ABX*I_ERI_D2x_S_Fxyz_S_b;
  Double I_ERI_Dxy_Px_Fxyz_S_b = I_ERI_F2xy_S_Fxyz_S_b+ABX*I_ERI_Dxy_S_Fxyz_S_b;
  Double I_ERI_Dxz_Px_Fxyz_S_b = I_ERI_F2xz_S_Fxyz_S_b+ABX*I_ERI_Dxz_S_Fxyz_S_b;
  Double I_ERI_D2y_Px_Fxyz_S_b = I_ERI_Fx2y_S_Fxyz_S_b+ABX*I_ERI_D2y_S_Fxyz_S_b;
  Double I_ERI_Dyz_Px_Fxyz_S_b = I_ERI_Fxyz_S_Fxyz_S_b+ABX*I_ERI_Dyz_S_Fxyz_S_b;
  Double I_ERI_D2z_Px_Fxyz_S_b = I_ERI_Fx2z_S_Fxyz_S_b+ABX*I_ERI_D2z_S_Fxyz_S_b;
  Double I_ERI_D2x_Py_Fxyz_S_b = I_ERI_F2xy_S_Fxyz_S_b+ABY*I_ERI_D2x_S_Fxyz_S_b;
  Double I_ERI_Dxy_Py_Fxyz_S_b = I_ERI_Fx2y_S_Fxyz_S_b+ABY*I_ERI_Dxy_S_Fxyz_S_b;
  Double I_ERI_Dxz_Py_Fxyz_S_b = I_ERI_Fxyz_S_Fxyz_S_b+ABY*I_ERI_Dxz_S_Fxyz_S_b;
  Double I_ERI_D2y_Py_Fxyz_S_b = I_ERI_F3y_S_Fxyz_S_b+ABY*I_ERI_D2y_S_Fxyz_S_b;
  Double I_ERI_Dyz_Py_Fxyz_S_b = I_ERI_F2yz_S_Fxyz_S_b+ABY*I_ERI_Dyz_S_Fxyz_S_b;
  Double I_ERI_D2z_Py_Fxyz_S_b = I_ERI_Fy2z_S_Fxyz_S_b+ABY*I_ERI_D2z_S_Fxyz_S_b;
  Double I_ERI_D2x_Pz_Fxyz_S_b = I_ERI_F2xz_S_Fxyz_S_b+ABZ*I_ERI_D2x_S_Fxyz_S_b;
  Double I_ERI_Dxy_Pz_Fxyz_S_b = I_ERI_Fxyz_S_Fxyz_S_b+ABZ*I_ERI_Dxy_S_Fxyz_S_b;
  Double I_ERI_Dxz_Pz_Fxyz_S_b = I_ERI_Fx2z_S_Fxyz_S_b+ABZ*I_ERI_Dxz_S_Fxyz_S_b;
  Double I_ERI_D2y_Pz_Fxyz_S_b = I_ERI_F2yz_S_Fxyz_S_b+ABZ*I_ERI_D2y_S_Fxyz_S_b;
  Double I_ERI_Dyz_Pz_Fxyz_S_b = I_ERI_Fy2z_S_Fxyz_S_b+ABZ*I_ERI_Dyz_S_Fxyz_S_b;
  Double I_ERI_D2z_Pz_Fxyz_S_b = I_ERI_F3z_S_Fxyz_S_b+ABZ*I_ERI_D2z_S_Fxyz_S_b;
  Double I_ERI_D2x_Px_Fx2z_S_b = I_ERI_F3x_S_Fx2z_S_b+ABX*I_ERI_D2x_S_Fx2z_S_b;
  Double I_ERI_Dxy_Px_Fx2z_S_b = I_ERI_F2xy_S_Fx2z_S_b+ABX*I_ERI_Dxy_S_Fx2z_S_b;
  Double I_ERI_Dxz_Px_Fx2z_S_b = I_ERI_F2xz_S_Fx2z_S_b+ABX*I_ERI_Dxz_S_Fx2z_S_b;
  Double I_ERI_D2y_Px_Fx2z_S_b = I_ERI_Fx2y_S_Fx2z_S_b+ABX*I_ERI_D2y_S_Fx2z_S_b;
  Double I_ERI_Dyz_Px_Fx2z_S_b = I_ERI_Fxyz_S_Fx2z_S_b+ABX*I_ERI_Dyz_S_Fx2z_S_b;
  Double I_ERI_D2z_Px_Fx2z_S_b = I_ERI_Fx2z_S_Fx2z_S_b+ABX*I_ERI_D2z_S_Fx2z_S_b;
  Double I_ERI_D2x_Py_Fx2z_S_b = I_ERI_F2xy_S_Fx2z_S_b+ABY*I_ERI_D2x_S_Fx2z_S_b;
  Double I_ERI_Dxy_Py_Fx2z_S_b = I_ERI_Fx2y_S_Fx2z_S_b+ABY*I_ERI_Dxy_S_Fx2z_S_b;
  Double I_ERI_Dxz_Py_Fx2z_S_b = I_ERI_Fxyz_S_Fx2z_S_b+ABY*I_ERI_Dxz_S_Fx2z_S_b;
  Double I_ERI_D2y_Py_Fx2z_S_b = I_ERI_F3y_S_Fx2z_S_b+ABY*I_ERI_D2y_S_Fx2z_S_b;
  Double I_ERI_Dyz_Py_Fx2z_S_b = I_ERI_F2yz_S_Fx2z_S_b+ABY*I_ERI_Dyz_S_Fx2z_S_b;
  Double I_ERI_D2z_Py_Fx2z_S_b = I_ERI_Fy2z_S_Fx2z_S_b+ABY*I_ERI_D2z_S_Fx2z_S_b;
  Double I_ERI_D2x_Pz_Fx2z_S_b = I_ERI_F2xz_S_Fx2z_S_b+ABZ*I_ERI_D2x_S_Fx2z_S_b;
  Double I_ERI_Dxy_Pz_Fx2z_S_b = I_ERI_Fxyz_S_Fx2z_S_b+ABZ*I_ERI_Dxy_S_Fx2z_S_b;
  Double I_ERI_Dxz_Pz_Fx2z_S_b = I_ERI_Fx2z_S_Fx2z_S_b+ABZ*I_ERI_Dxz_S_Fx2z_S_b;
  Double I_ERI_D2y_Pz_Fx2z_S_b = I_ERI_F2yz_S_Fx2z_S_b+ABZ*I_ERI_D2y_S_Fx2z_S_b;
  Double I_ERI_Dyz_Pz_Fx2z_S_b = I_ERI_Fy2z_S_Fx2z_S_b+ABZ*I_ERI_Dyz_S_Fx2z_S_b;
  Double I_ERI_D2z_Pz_Fx2z_S_b = I_ERI_F3z_S_Fx2z_S_b+ABZ*I_ERI_D2z_S_Fx2z_S_b;
  Double I_ERI_D2x_Px_F3y_S_b = I_ERI_F3x_S_F3y_S_b+ABX*I_ERI_D2x_S_F3y_S_b;
  Double I_ERI_Dxy_Px_F3y_S_b = I_ERI_F2xy_S_F3y_S_b+ABX*I_ERI_Dxy_S_F3y_S_b;
  Double I_ERI_Dxz_Px_F3y_S_b = I_ERI_F2xz_S_F3y_S_b+ABX*I_ERI_Dxz_S_F3y_S_b;
  Double I_ERI_D2y_Px_F3y_S_b = I_ERI_Fx2y_S_F3y_S_b+ABX*I_ERI_D2y_S_F3y_S_b;
  Double I_ERI_Dyz_Px_F3y_S_b = I_ERI_Fxyz_S_F3y_S_b+ABX*I_ERI_Dyz_S_F3y_S_b;
  Double I_ERI_D2z_Px_F3y_S_b = I_ERI_Fx2z_S_F3y_S_b+ABX*I_ERI_D2z_S_F3y_S_b;
  Double I_ERI_D2x_Py_F3y_S_b = I_ERI_F2xy_S_F3y_S_b+ABY*I_ERI_D2x_S_F3y_S_b;
  Double I_ERI_Dxy_Py_F3y_S_b = I_ERI_Fx2y_S_F3y_S_b+ABY*I_ERI_Dxy_S_F3y_S_b;
  Double I_ERI_Dxz_Py_F3y_S_b = I_ERI_Fxyz_S_F3y_S_b+ABY*I_ERI_Dxz_S_F3y_S_b;
  Double I_ERI_D2y_Py_F3y_S_b = I_ERI_F3y_S_F3y_S_b+ABY*I_ERI_D2y_S_F3y_S_b;
  Double I_ERI_Dyz_Py_F3y_S_b = I_ERI_F2yz_S_F3y_S_b+ABY*I_ERI_Dyz_S_F3y_S_b;
  Double I_ERI_D2z_Py_F3y_S_b = I_ERI_Fy2z_S_F3y_S_b+ABY*I_ERI_D2z_S_F3y_S_b;
  Double I_ERI_D2x_Pz_F3y_S_b = I_ERI_F2xz_S_F3y_S_b+ABZ*I_ERI_D2x_S_F3y_S_b;
  Double I_ERI_Dxy_Pz_F3y_S_b = I_ERI_Fxyz_S_F3y_S_b+ABZ*I_ERI_Dxy_S_F3y_S_b;
  Double I_ERI_Dxz_Pz_F3y_S_b = I_ERI_Fx2z_S_F3y_S_b+ABZ*I_ERI_Dxz_S_F3y_S_b;
  Double I_ERI_D2y_Pz_F3y_S_b = I_ERI_F2yz_S_F3y_S_b+ABZ*I_ERI_D2y_S_F3y_S_b;
  Double I_ERI_Dyz_Pz_F3y_S_b = I_ERI_Fy2z_S_F3y_S_b+ABZ*I_ERI_Dyz_S_F3y_S_b;
  Double I_ERI_D2z_Pz_F3y_S_b = I_ERI_F3z_S_F3y_S_b+ABZ*I_ERI_D2z_S_F3y_S_b;
  Double I_ERI_D2x_Px_F2yz_S_b = I_ERI_F3x_S_F2yz_S_b+ABX*I_ERI_D2x_S_F2yz_S_b;
  Double I_ERI_Dxy_Px_F2yz_S_b = I_ERI_F2xy_S_F2yz_S_b+ABX*I_ERI_Dxy_S_F2yz_S_b;
  Double I_ERI_Dxz_Px_F2yz_S_b = I_ERI_F2xz_S_F2yz_S_b+ABX*I_ERI_Dxz_S_F2yz_S_b;
  Double I_ERI_D2y_Px_F2yz_S_b = I_ERI_Fx2y_S_F2yz_S_b+ABX*I_ERI_D2y_S_F2yz_S_b;
  Double I_ERI_Dyz_Px_F2yz_S_b = I_ERI_Fxyz_S_F2yz_S_b+ABX*I_ERI_Dyz_S_F2yz_S_b;
  Double I_ERI_D2z_Px_F2yz_S_b = I_ERI_Fx2z_S_F2yz_S_b+ABX*I_ERI_D2z_S_F2yz_S_b;
  Double I_ERI_D2x_Py_F2yz_S_b = I_ERI_F2xy_S_F2yz_S_b+ABY*I_ERI_D2x_S_F2yz_S_b;
  Double I_ERI_Dxy_Py_F2yz_S_b = I_ERI_Fx2y_S_F2yz_S_b+ABY*I_ERI_Dxy_S_F2yz_S_b;
  Double I_ERI_Dxz_Py_F2yz_S_b = I_ERI_Fxyz_S_F2yz_S_b+ABY*I_ERI_Dxz_S_F2yz_S_b;
  Double I_ERI_D2y_Py_F2yz_S_b = I_ERI_F3y_S_F2yz_S_b+ABY*I_ERI_D2y_S_F2yz_S_b;
  Double I_ERI_Dyz_Py_F2yz_S_b = I_ERI_F2yz_S_F2yz_S_b+ABY*I_ERI_Dyz_S_F2yz_S_b;
  Double I_ERI_D2z_Py_F2yz_S_b = I_ERI_Fy2z_S_F2yz_S_b+ABY*I_ERI_D2z_S_F2yz_S_b;
  Double I_ERI_D2x_Pz_F2yz_S_b = I_ERI_F2xz_S_F2yz_S_b+ABZ*I_ERI_D2x_S_F2yz_S_b;
  Double I_ERI_Dxy_Pz_F2yz_S_b = I_ERI_Fxyz_S_F2yz_S_b+ABZ*I_ERI_Dxy_S_F2yz_S_b;
  Double I_ERI_Dxz_Pz_F2yz_S_b = I_ERI_Fx2z_S_F2yz_S_b+ABZ*I_ERI_Dxz_S_F2yz_S_b;
  Double I_ERI_D2y_Pz_F2yz_S_b = I_ERI_F2yz_S_F2yz_S_b+ABZ*I_ERI_D2y_S_F2yz_S_b;
  Double I_ERI_Dyz_Pz_F2yz_S_b = I_ERI_Fy2z_S_F2yz_S_b+ABZ*I_ERI_Dyz_S_F2yz_S_b;
  Double I_ERI_D2z_Pz_F2yz_S_b = I_ERI_F3z_S_F2yz_S_b+ABZ*I_ERI_D2z_S_F2yz_S_b;
  Double I_ERI_D2x_Px_Fy2z_S_b = I_ERI_F3x_S_Fy2z_S_b+ABX*I_ERI_D2x_S_Fy2z_S_b;
  Double I_ERI_Dxy_Px_Fy2z_S_b = I_ERI_F2xy_S_Fy2z_S_b+ABX*I_ERI_Dxy_S_Fy2z_S_b;
  Double I_ERI_Dxz_Px_Fy2z_S_b = I_ERI_F2xz_S_Fy2z_S_b+ABX*I_ERI_Dxz_S_Fy2z_S_b;
  Double I_ERI_D2y_Px_Fy2z_S_b = I_ERI_Fx2y_S_Fy2z_S_b+ABX*I_ERI_D2y_S_Fy2z_S_b;
  Double I_ERI_Dyz_Px_Fy2z_S_b = I_ERI_Fxyz_S_Fy2z_S_b+ABX*I_ERI_Dyz_S_Fy2z_S_b;
  Double I_ERI_D2z_Px_Fy2z_S_b = I_ERI_Fx2z_S_Fy2z_S_b+ABX*I_ERI_D2z_S_Fy2z_S_b;
  Double I_ERI_D2x_Py_Fy2z_S_b = I_ERI_F2xy_S_Fy2z_S_b+ABY*I_ERI_D2x_S_Fy2z_S_b;
  Double I_ERI_Dxy_Py_Fy2z_S_b = I_ERI_Fx2y_S_Fy2z_S_b+ABY*I_ERI_Dxy_S_Fy2z_S_b;
  Double I_ERI_Dxz_Py_Fy2z_S_b = I_ERI_Fxyz_S_Fy2z_S_b+ABY*I_ERI_Dxz_S_Fy2z_S_b;
  Double I_ERI_D2y_Py_Fy2z_S_b = I_ERI_F3y_S_Fy2z_S_b+ABY*I_ERI_D2y_S_Fy2z_S_b;
  Double I_ERI_Dyz_Py_Fy2z_S_b = I_ERI_F2yz_S_Fy2z_S_b+ABY*I_ERI_Dyz_S_Fy2z_S_b;
  Double I_ERI_D2z_Py_Fy2z_S_b = I_ERI_Fy2z_S_Fy2z_S_b+ABY*I_ERI_D2z_S_Fy2z_S_b;
  Double I_ERI_D2x_Pz_Fy2z_S_b = I_ERI_F2xz_S_Fy2z_S_b+ABZ*I_ERI_D2x_S_Fy2z_S_b;
  Double I_ERI_Dxy_Pz_Fy2z_S_b = I_ERI_Fxyz_S_Fy2z_S_b+ABZ*I_ERI_Dxy_S_Fy2z_S_b;
  Double I_ERI_Dxz_Pz_Fy2z_S_b = I_ERI_Fx2z_S_Fy2z_S_b+ABZ*I_ERI_Dxz_S_Fy2z_S_b;
  Double I_ERI_D2y_Pz_Fy2z_S_b = I_ERI_F2yz_S_Fy2z_S_b+ABZ*I_ERI_D2y_S_Fy2z_S_b;
  Double I_ERI_Dyz_Pz_Fy2z_S_b = I_ERI_Fy2z_S_Fy2z_S_b+ABZ*I_ERI_Dyz_S_Fy2z_S_b;
  Double I_ERI_D2z_Pz_Fy2z_S_b = I_ERI_F3z_S_Fy2z_S_b+ABZ*I_ERI_D2z_S_Fy2z_S_b;
  Double I_ERI_D2x_Px_F3z_S_b = I_ERI_F3x_S_F3z_S_b+ABX*I_ERI_D2x_S_F3z_S_b;
  Double I_ERI_Dxy_Px_F3z_S_b = I_ERI_F2xy_S_F3z_S_b+ABX*I_ERI_Dxy_S_F3z_S_b;
  Double I_ERI_Dxz_Px_F3z_S_b = I_ERI_F2xz_S_F3z_S_b+ABX*I_ERI_Dxz_S_F3z_S_b;
  Double I_ERI_D2y_Px_F3z_S_b = I_ERI_Fx2y_S_F3z_S_b+ABX*I_ERI_D2y_S_F3z_S_b;
  Double I_ERI_Dyz_Px_F3z_S_b = I_ERI_Fxyz_S_F3z_S_b+ABX*I_ERI_Dyz_S_F3z_S_b;
  Double I_ERI_D2z_Px_F3z_S_b = I_ERI_Fx2z_S_F3z_S_b+ABX*I_ERI_D2z_S_F3z_S_b;
  Double I_ERI_D2x_Py_F3z_S_b = I_ERI_F2xy_S_F3z_S_b+ABY*I_ERI_D2x_S_F3z_S_b;
  Double I_ERI_Dxy_Py_F3z_S_b = I_ERI_Fx2y_S_F3z_S_b+ABY*I_ERI_Dxy_S_F3z_S_b;
  Double I_ERI_Dxz_Py_F3z_S_b = I_ERI_Fxyz_S_F3z_S_b+ABY*I_ERI_Dxz_S_F3z_S_b;
  Double I_ERI_D2y_Py_F3z_S_b = I_ERI_F3y_S_F3z_S_b+ABY*I_ERI_D2y_S_F3z_S_b;
  Double I_ERI_Dyz_Py_F3z_S_b = I_ERI_F2yz_S_F3z_S_b+ABY*I_ERI_Dyz_S_F3z_S_b;
  Double I_ERI_D2z_Py_F3z_S_b = I_ERI_Fy2z_S_F3z_S_b+ABY*I_ERI_D2z_S_F3z_S_b;
  Double I_ERI_D2x_Pz_F3z_S_b = I_ERI_F2xz_S_F3z_S_b+ABZ*I_ERI_D2x_S_F3z_S_b;
  Double I_ERI_Dxy_Pz_F3z_S_b = I_ERI_Fxyz_S_F3z_S_b+ABZ*I_ERI_Dxy_S_F3z_S_b;
  Double I_ERI_Dxz_Pz_F3z_S_b = I_ERI_Fx2z_S_F3z_S_b+ABZ*I_ERI_Dxz_S_F3z_S_b;
  Double I_ERI_D2y_Pz_F3z_S_b = I_ERI_F2yz_S_F3z_S_b+ABZ*I_ERI_D2y_S_F3z_S_b;
  Double I_ERI_Dyz_Pz_F3z_S_b = I_ERI_Fy2z_S_F3z_S_b+ABZ*I_ERI_Dyz_S_F3z_S_b;
  Double I_ERI_D2z_Pz_F3z_S_b = I_ERI_F3z_S_F3z_S_b+ABZ*I_ERI_D2z_S_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_F_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 50 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_b
   * RHS shell quartet name: SQ_ERI_F_S_F_S_b
   ************************************************************/
  Double I_ERI_F3x_Px_F3x_S_b = I_ERI_G4x_S_F3x_S_b+ABX*I_ERI_F3x_S_F3x_S_b;
  Double I_ERI_F2xy_Px_F3x_S_b = I_ERI_G3xy_S_F3x_S_b+ABX*I_ERI_F2xy_S_F3x_S_b;
  Double I_ERI_F2xz_Px_F3x_S_b = I_ERI_G3xz_S_F3x_S_b+ABX*I_ERI_F2xz_S_F3x_S_b;
  Double I_ERI_Fx2y_Px_F3x_S_b = I_ERI_G2x2y_S_F3x_S_b+ABX*I_ERI_Fx2y_S_F3x_S_b;
  Double I_ERI_Fxyz_Px_F3x_S_b = I_ERI_G2xyz_S_F3x_S_b+ABX*I_ERI_Fxyz_S_F3x_S_b;
  Double I_ERI_Fx2z_Px_F3x_S_b = I_ERI_G2x2z_S_F3x_S_b+ABX*I_ERI_Fx2z_S_F3x_S_b;
  Double I_ERI_F3y_Px_F3x_S_b = I_ERI_Gx3y_S_F3x_S_b+ABX*I_ERI_F3y_S_F3x_S_b;
  Double I_ERI_F2yz_Px_F3x_S_b = I_ERI_Gx2yz_S_F3x_S_b+ABX*I_ERI_F2yz_S_F3x_S_b;
  Double I_ERI_Fy2z_Px_F3x_S_b = I_ERI_Gxy2z_S_F3x_S_b+ABX*I_ERI_Fy2z_S_F3x_S_b;
  Double I_ERI_F3z_Px_F3x_S_b = I_ERI_Gx3z_S_F3x_S_b+ABX*I_ERI_F3z_S_F3x_S_b;
  Double I_ERI_F2xy_Py_F3x_S_b = I_ERI_G2x2y_S_F3x_S_b+ABY*I_ERI_F2xy_S_F3x_S_b;
  Double I_ERI_F2xz_Py_F3x_S_b = I_ERI_G2xyz_S_F3x_S_b+ABY*I_ERI_F2xz_S_F3x_S_b;
  Double I_ERI_Fx2y_Py_F3x_S_b = I_ERI_Gx3y_S_F3x_S_b+ABY*I_ERI_Fx2y_S_F3x_S_b;
  Double I_ERI_Fxyz_Py_F3x_S_b = I_ERI_Gx2yz_S_F3x_S_b+ABY*I_ERI_Fxyz_S_F3x_S_b;
  Double I_ERI_Fx2z_Py_F3x_S_b = I_ERI_Gxy2z_S_F3x_S_b+ABY*I_ERI_Fx2z_S_F3x_S_b;
  Double I_ERI_F3y_Py_F3x_S_b = I_ERI_G4y_S_F3x_S_b+ABY*I_ERI_F3y_S_F3x_S_b;
  Double I_ERI_F2yz_Py_F3x_S_b = I_ERI_G3yz_S_F3x_S_b+ABY*I_ERI_F2yz_S_F3x_S_b;
  Double I_ERI_Fy2z_Py_F3x_S_b = I_ERI_G2y2z_S_F3x_S_b+ABY*I_ERI_Fy2z_S_F3x_S_b;
  Double I_ERI_F3z_Py_F3x_S_b = I_ERI_Gy3z_S_F3x_S_b+ABY*I_ERI_F3z_S_F3x_S_b;
  Double I_ERI_F2xz_Pz_F3x_S_b = I_ERI_G2x2z_S_F3x_S_b+ABZ*I_ERI_F2xz_S_F3x_S_b;
  Double I_ERI_Fxyz_Pz_F3x_S_b = I_ERI_Gxy2z_S_F3x_S_b+ABZ*I_ERI_Fxyz_S_F3x_S_b;
  Double I_ERI_Fx2z_Pz_F3x_S_b = I_ERI_Gx3z_S_F3x_S_b+ABZ*I_ERI_Fx2z_S_F3x_S_b;
  Double I_ERI_F2yz_Pz_F3x_S_b = I_ERI_G2y2z_S_F3x_S_b+ABZ*I_ERI_F2yz_S_F3x_S_b;
  Double I_ERI_Fy2z_Pz_F3x_S_b = I_ERI_Gy3z_S_F3x_S_b+ABZ*I_ERI_Fy2z_S_F3x_S_b;
  Double I_ERI_F3z_Pz_F3x_S_b = I_ERI_G4z_S_F3x_S_b+ABZ*I_ERI_F3z_S_F3x_S_b;
  Double I_ERI_F3x_Px_F2xy_S_b = I_ERI_G4x_S_F2xy_S_b+ABX*I_ERI_F3x_S_F2xy_S_b;
  Double I_ERI_F2xy_Px_F2xy_S_b = I_ERI_G3xy_S_F2xy_S_b+ABX*I_ERI_F2xy_S_F2xy_S_b;
  Double I_ERI_F2xz_Px_F2xy_S_b = I_ERI_G3xz_S_F2xy_S_b+ABX*I_ERI_F2xz_S_F2xy_S_b;
  Double I_ERI_Fx2y_Px_F2xy_S_b = I_ERI_G2x2y_S_F2xy_S_b+ABX*I_ERI_Fx2y_S_F2xy_S_b;
  Double I_ERI_Fxyz_Px_F2xy_S_b = I_ERI_G2xyz_S_F2xy_S_b+ABX*I_ERI_Fxyz_S_F2xy_S_b;
  Double I_ERI_Fx2z_Px_F2xy_S_b = I_ERI_G2x2z_S_F2xy_S_b+ABX*I_ERI_Fx2z_S_F2xy_S_b;
  Double I_ERI_F3y_Px_F2xy_S_b = I_ERI_Gx3y_S_F2xy_S_b+ABX*I_ERI_F3y_S_F2xy_S_b;
  Double I_ERI_F2yz_Px_F2xy_S_b = I_ERI_Gx2yz_S_F2xy_S_b+ABX*I_ERI_F2yz_S_F2xy_S_b;
  Double I_ERI_Fy2z_Px_F2xy_S_b = I_ERI_Gxy2z_S_F2xy_S_b+ABX*I_ERI_Fy2z_S_F2xy_S_b;
  Double I_ERI_F3z_Px_F2xy_S_b = I_ERI_Gx3z_S_F2xy_S_b+ABX*I_ERI_F3z_S_F2xy_S_b;
  Double I_ERI_F2xy_Py_F2xy_S_b = I_ERI_G2x2y_S_F2xy_S_b+ABY*I_ERI_F2xy_S_F2xy_S_b;
  Double I_ERI_F2xz_Py_F2xy_S_b = I_ERI_G2xyz_S_F2xy_S_b+ABY*I_ERI_F2xz_S_F2xy_S_b;
  Double I_ERI_Fx2y_Py_F2xy_S_b = I_ERI_Gx3y_S_F2xy_S_b+ABY*I_ERI_Fx2y_S_F2xy_S_b;
  Double I_ERI_Fxyz_Py_F2xy_S_b = I_ERI_Gx2yz_S_F2xy_S_b+ABY*I_ERI_Fxyz_S_F2xy_S_b;
  Double I_ERI_Fx2z_Py_F2xy_S_b = I_ERI_Gxy2z_S_F2xy_S_b+ABY*I_ERI_Fx2z_S_F2xy_S_b;
  Double I_ERI_F3y_Py_F2xy_S_b = I_ERI_G4y_S_F2xy_S_b+ABY*I_ERI_F3y_S_F2xy_S_b;
  Double I_ERI_F2yz_Py_F2xy_S_b = I_ERI_G3yz_S_F2xy_S_b+ABY*I_ERI_F2yz_S_F2xy_S_b;
  Double I_ERI_Fy2z_Py_F2xy_S_b = I_ERI_G2y2z_S_F2xy_S_b+ABY*I_ERI_Fy2z_S_F2xy_S_b;
  Double I_ERI_F3z_Py_F2xy_S_b = I_ERI_Gy3z_S_F2xy_S_b+ABY*I_ERI_F3z_S_F2xy_S_b;
  Double I_ERI_F2xz_Pz_F2xy_S_b = I_ERI_G2x2z_S_F2xy_S_b+ABZ*I_ERI_F2xz_S_F2xy_S_b;
  Double I_ERI_Fxyz_Pz_F2xy_S_b = I_ERI_Gxy2z_S_F2xy_S_b+ABZ*I_ERI_Fxyz_S_F2xy_S_b;
  Double I_ERI_Fx2z_Pz_F2xy_S_b = I_ERI_Gx3z_S_F2xy_S_b+ABZ*I_ERI_Fx2z_S_F2xy_S_b;
  Double I_ERI_F2yz_Pz_F2xy_S_b = I_ERI_G2y2z_S_F2xy_S_b+ABZ*I_ERI_F2yz_S_F2xy_S_b;
  Double I_ERI_Fy2z_Pz_F2xy_S_b = I_ERI_Gy3z_S_F2xy_S_b+ABZ*I_ERI_Fy2z_S_F2xy_S_b;
  Double I_ERI_F3z_Pz_F2xy_S_b = I_ERI_G4z_S_F2xy_S_b+ABZ*I_ERI_F3z_S_F2xy_S_b;
  Double I_ERI_F3x_Px_F2xz_S_b = I_ERI_G4x_S_F2xz_S_b+ABX*I_ERI_F3x_S_F2xz_S_b;
  Double I_ERI_F2xy_Px_F2xz_S_b = I_ERI_G3xy_S_F2xz_S_b+ABX*I_ERI_F2xy_S_F2xz_S_b;
  Double I_ERI_F2xz_Px_F2xz_S_b = I_ERI_G3xz_S_F2xz_S_b+ABX*I_ERI_F2xz_S_F2xz_S_b;
  Double I_ERI_Fx2y_Px_F2xz_S_b = I_ERI_G2x2y_S_F2xz_S_b+ABX*I_ERI_Fx2y_S_F2xz_S_b;
  Double I_ERI_Fxyz_Px_F2xz_S_b = I_ERI_G2xyz_S_F2xz_S_b+ABX*I_ERI_Fxyz_S_F2xz_S_b;
  Double I_ERI_Fx2z_Px_F2xz_S_b = I_ERI_G2x2z_S_F2xz_S_b+ABX*I_ERI_Fx2z_S_F2xz_S_b;
  Double I_ERI_F3y_Px_F2xz_S_b = I_ERI_Gx3y_S_F2xz_S_b+ABX*I_ERI_F3y_S_F2xz_S_b;
  Double I_ERI_F2yz_Px_F2xz_S_b = I_ERI_Gx2yz_S_F2xz_S_b+ABX*I_ERI_F2yz_S_F2xz_S_b;
  Double I_ERI_Fy2z_Px_F2xz_S_b = I_ERI_Gxy2z_S_F2xz_S_b+ABX*I_ERI_Fy2z_S_F2xz_S_b;
  Double I_ERI_F3z_Px_F2xz_S_b = I_ERI_Gx3z_S_F2xz_S_b+ABX*I_ERI_F3z_S_F2xz_S_b;
  Double I_ERI_F2xy_Py_F2xz_S_b = I_ERI_G2x2y_S_F2xz_S_b+ABY*I_ERI_F2xy_S_F2xz_S_b;
  Double I_ERI_F2xz_Py_F2xz_S_b = I_ERI_G2xyz_S_F2xz_S_b+ABY*I_ERI_F2xz_S_F2xz_S_b;
  Double I_ERI_Fx2y_Py_F2xz_S_b = I_ERI_Gx3y_S_F2xz_S_b+ABY*I_ERI_Fx2y_S_F2xz_S_b;
  Double I_ERI_Fxyz_Py_F2xz_S_b = I_ERI_Gx2yz_S_F2xz_S_b+ABY*I_ERI_Fxyz_S_F2xz_S_b;
  Double I_ERI_Fx2z_Py_F2xz_S_b = I_ERI_Gxy2z_S_F2xz_S_b+ABY*I_ERI_Fx2z_S_F2xz_S_b;
  Double I_ERI_F3y_Py_F2xz_S_b = I_ERI_G4y_S_F2xz_S_b+ABY*I_ERI_F3y_S_F2xz_S_b;
  Double I_ERI_F2yz_Py_F2xz_S_b = I_ERI_G3yz_S_F2xz_S_b+ABY*I_ERI_F2yz_S_F2xz_S_b;
  Double I_ERI_Fy2z_Py_F2xz_S_b = I_ERI_G2y2z_S_F2xz_S_b+ABY*I_ERI_Fy2z_S_F2xz_S_b;
  Double I_ERI_F3z_Py_F2xz_S_b = I_ERI_Gy3z_S_F2xz_S_b+ABY*I_ERI_F3z_S_F2xz_S_b;
  Double I_ERI_F2xz_Pz_F2xz_S_b = I_ERI_G2x2z_S_F2xz_S_b+ABZ*I_ERI_F2xz_S_F2xz_S_b;
  Double I_ERI_Fxyz_Pz_F2xz_S_b = I_ERI_Gxy2z_S_F2xz_S_b+ABZ*I_ERI_Fxyz_S_F2xz_S_b;
  Double I_ERI_Fx2z_Pz_F2xz_S_b = I_ERI_Gx3z_S_F2xz_S_b+ABZ*I_ERI_Fx2z_S_F2xz_S_b;
  Double I_ERI_F2yz_Pz_F2xz_S_b = I_ERI_G2y2z_S_F2xz_S_b+ABZ*I_ERI_F2yz_S_F2xz_S_b;
  Double I_ERI_Fy2z_Pz_F2xz_S_b = I_ERI_Gy3z_S_F2xz_S_b+ABZ*I_ERI_Fy2z_S_F2xz_S_b;
  Double I_ERI_F3z_Pz_F2xz_S_b = I_ERI_G4z_S_F2xz_S_b+ABZ*I_ERI_F3z_S_F2xz_S_b;
  Double I_ERI_F3x_Px_Fx2y_S_b = I_ERI_G4x_S_Fx2y_S_b+ABX*I_ERI_F3x_S_Fx2y_S_b;
  Double I_ERI_F2xy_Px_Fx2y_S_b = I_ERI_G3xy_S_Fx2y_S_b+ABX*I_ERI_F2xy_S_Fx2y_S_b;
  Double I_ERI_F2xz_Px_Fx2y_S_b = I_ERI_G3xz_S_Fx2y_S_b+ABX*I_ERI_F2xz_S_Fx2y_S_b;
  Double I_ERI_Fx2y_Px_Fx2y_S_b = I_ERI_G2x2y_S_Fx2y_S_b+ABX*I_ERI_Fx2y_S_Fx2y_S_b;
  Double I_ERI_Fxyz_Px_Fx2y_S_b = I_ERI_G2xyz_S_Fx2y_S_b+ABX*I_ERI_Fxyz_S_Fx2y_S_b;
  Double I_ERI_Fx2z_Px_Fx2y_S_b = I_ERI_G2x2z_S_Fx2y_S_b+ABX*I_ERI_Fx2z_S_Fx2y_S_b;
  Double I_ERI_F3y_Px_Fx2y_S_b = I_ERI_Gx3y_S_Fx2y_S_b+ABX*I_ERI_F3y_S_Fx2y_S_b;
  Double I_ERI_F2yz_Px_Fx2y_S_b = I_ERI_Gx2yz_S_Fx2y_S_b+ABX*I_ERI_F2yz_S_Fx2y_S_b;
  Double I_ERI_Fy2z_Px_Fx2y_S_b = I_ERI_Gxy2z_S_Fx2y_S_b+ABX*I_ERI_Fy2z_S_Fx2y_S_b;
  Double I_ERI_F3z_Px_Fx2y_S_b = I_ERI_Gx3z_S_Fx2y_S_b+ABX*I_ERI_F3z_S_Fx2y_S_b;
  Double I_ERI_F2xy_Py_Fx2y_S_b = I_ERI_G2x2y_S_Fx2y_S_b+ABY*I_ERI_F2xy_S_Fx2y_S_b;
  Double I_ERI_F2xz_Py_Fx2y_S_b = I_ERI_G2xyz_S_Fx2y_S_b+ABY*I_ERI_F2xz_S_Fx2y_S_b;
  Double I_ERI_Fx2y_Py_Fx2y_S_b = I_ERI_Gx3y_S_Fx2y_S_b+ABY*I_ERI_Fx2y_S_Fx2y_S_b;
  Double I_ERI_Fxyz_Py_Fx2y_S_b = I_ERI_Gx2yz_S_Fx2y_S_b+ABY*I_ERI_Fxyz_S_Fx2y_S_b;
  Double I_ERI_Fx2z_Py_Fx2y_S_b = I_ERI_Gxy2z_S_Fx2y_S_b+ABY*I_ERI_Fx2z_S_Fx2y_S_b;
  Double I_ERI_F3y_Py_Fx2y_S_b = I_ERI_G4y_S_Fx2y_S_b+ABY*I_ERI_F3y_S_Fx2y_S_b;
  Double I_ERI_F2yz_Py_Fx2y_S_b = I_ERI_G3yz_S_Fx2y_S_b+ABY*I_ERI_F2yz_S_Fx2y_S_b;
  Double I_ERI_Fy2z_Py_Fx2y_S_b = I_ERI_G2y2z_S_Fx2y_S_b+ABY*I_ERI_Fy2z_S_Fx2y_S_b;
  Double I_ERI_F3z_Py_Fx2y_S_b = I_ERI_Gy3z_S_Fx2y_S_b+ABY*I_ERI_F3z_S_Fx2y_S_b;
  Double I_ERI_F2xz_Pz_Fx2y_S_b = I_ERI_G2x2z_S_Fx2y_S_b+ABZ*I_ERI_F2xz_S_Fx2y_S_b;
  Double I_ERI_Fxyz_Pz_Fx2y_S_b = I_ERI_Gxy2z_S_Fx2y_S_b+ABZ*I_ERI_Fxyz_S_Fx2y_S_b;
  Double I_ERI_Fx2z_Pz_Fx2y_S_b = I_ERI_Gx3z_S_Fx2y_S_b+ABZ*I_ERI_Fx2z_S_Fx2y_S_b;
  Double I_ERI_F2yz_Pz_Fx2y_S_b = I_ERI_G2y2z_S_Fx2y_S_b+ABZ*I_ERI_F2yz_S_Fx2y_S_b;
  Double I_ERI_Fy2z_Pz_Fx2y_S_b = I_ERI_Gy3z_S_Fx2y_S_b+ABZ*I_ERI_Fy2z_S_Fx2y_S_b;
  Double I_ERI_F3z_Pz_Fx2y_S_b = I_ERI_G4z_S_Fx2y_S_b+ABZ*I_ERI_F3z_S_Fx2y_S_b;
  Double I_ERI_F3x_Px_Fxyz_S_b = I_ERI_G4x_S_Fxyz_S_b+ABX*I_ERI_F3x_S_Fxyz_S_b;
  Double I_ERI_F2xy_Px_Fxyz_S_b = I_ERI_G3xy_S_Fxyz_S_b+ABX*I_ERI_F2xy_S_Fxyz_S_b;
  Double I_ERI_F2xz_Px_Fxyz_S_b = I_ERI_G3xz_S_Fxyz_S_b+ABX*I_ERI_F2xz_S_Fxyz_S_b;
  Double I_ERI_Fx2y_Px_Fxyz_S_b = I_ERI_G2x2y_S_Fxyz_S_b+ABX*I_ERI_Fx2y_S_Fxyz_S_b;
  Double I_ERI_Fxyz_Px_Fxyz_S_b = I_ERI_G2xyz_S_Fxyz_S_b+ABX*I_ERI_Fxyz_S_Fxyz_S_b;
  Double I_ERI_Fx2z_Px_Fxyz_S_b = I_ERI_G2x2z_S_Fxyz_S_b+ABX*I_ERI_Fx2z_S_Fxyz_S_b;
  Double I_ERI_F3y_Px_Fxyz_S_b = I_ERI_Gx3y_S_Fxyz_S_b+ABX*I_ERI_F3y_S_Fxyz_S_b;
  Double I_ERI_F2yz_Px_Fxyz_S_b = I_ERI_Gx2yz_S_Fxyz_S_b+ABX*I_ERI_F2yz_S_Fxyz_S_b;
  Double I_ERI_Fy2z_Px_Fxyz_S_b = I_ERI_Gxy2z_S_Fxyz_S_b+ABX*I_ERI_Fy2z_S_Fxyz_S_b;
  Double I_ERI_F3z_Px_Fxyz_S_b = I_ERI_Gx3z_S_Fxyz_S_b+ABX*I_ERI_F3z_S_Fxyz_S_b;
  Double I_ERI_F2xy_Py_Fxyz_S_b = I_ERI_G2x2y_S_Fxyz_S_b+ABY*I_ERI_F2xy_S_Fxyz_S_b;
  Double I_ERI_F2xz_Py_Fxyz_S_b = I_ERI_G2xyz_S_Fxyz_S_b+ABY*I_ERI_F2xz_S_Fxyz_S_b;
  Double I_ERI_Fx2y_Py_Fxyz_S_b = I_ERI_Gx3y_S_Fxyz_S_b+ABY*I_ERI_Fx2y_S_Fxyz_S_b;
  Double I_ERI_Fxyz_Py_Fxyz_S_b = I_ERI_Gx2yz_S_Fxyz_S_b+ABY*I_ERI_Fxyz_S_Fxyz_S_b;
  Double I_ERI_Fx2z_Py_Fxyz_S_b = I_ERI_Gxy2z_S_Fxyz_S_b+ABY*I_ERI_Fx2z_S_Fxyz_S_b;
  Double I_ERI_F3y_Py_Fxyz_S_b = I_ERI_G4y_S_Fxyz_S_b+ABY*I_ERI_F3y_S_Fxyz_S_b;
  Double I_ERI_F2yz_Py_Fxyz_S_b = I_ERI_G3yz_S_Fxyz_S_b+ABY*I_ERI_F2yz_S_Fxyz_S_b;
  Double I_ERI_Fy2z_Py_Fxyz_S_b = I_ERI_G2y2z_S_Fxyz_S_b+ABY*I_ERI_Fy2z_S_Fxyz_S_b;
  Double I_ERI_F3z_Py_Fxyz_S_b = I_ERI_Gy3z_S_Fxyz_S_b+ABY*I_ERI_F3z_S_Fxyz_S_b;
  Double I_ERI_F2xz_Pz_Fxyz_S_b = I_ERI_G2x2z_S_Fxyz_S_b+ABZ*I_ERI_F2xz_S_Fxyz_S_b;
  Double I_ERI_Fxyz_Pz_Fxyz_S_b = I_ERI_Gxy2z_S_Fxyz_S_b+ABZ*I_ERI_Fxyz_S_Fxyz_S_b;
  Double I_ERI_Fx2z_Pz_Fxyz_S_b = I_ERI_Gx3z_S_Fxyz_S_b+ABZ*I_ERI_Fx2z_S_Fxyz_S_b;
  Double I_ERI_F2yz_Pz_Fxyz_S_b = I_ERI_G2y2z_S_Fxyz_S_b+ABZ*I_ERI_F2yz_S_Fxyz_S_b;
  Double I_ERI_Fy2z_Pz_Fxyz_S_b = I_ERI_Gy3z_S_Fxyz_S_b+ABZ*I_ERI_Fy2z_S_Fxyz_S_b;
  Double I_ERI_F3z_Pz_Fxyz_S_b = I_ERI_G4z_S_Fxyz_S_b+ABZ*I_ERI_F3z_S_Fxyz_S_b;
  Double I_ERI_F3x_Px_Fx2z_S_b = I_ERI_G4x_S_Fx2z_S_b+ABX*I_ERI_F3x_S_Fx2z_S_b;
  Double I_ERI_F2xy_Px_Fx2z_S_b = I_ERI_G3xy_S_Fx2z_S_b+ABX*I_ERI_F2xy_S_Fx2z_S_b;
  Double I_ERI_F2xz_Px_Fx2z_S_b = I_ERI_G3xz_S_Fx2z_S_b+ABX*I_ERI_F2xz_S_Fx2z_S_b;
  Double I_ERI_Fx2y_Px_Fx2z_S_b = I_ERI_G2x2y_S_Fx2z_S_b+ABX*I_ERI_Fx2y_S_Fx2z_S_b;
  Double I_ERI_Fxyz_Px_Fx2z_S_b = I_ERI_G2xyz_S_Fx2z_S_b+ABX*I_ERI_Fxyz_S_Fx2z_S_b;
  Double I_ERI_Fx2z_Px_Fx2z_S_b = I_ERI_G2x2z_S_Fx2z_S_b+ABX*I_ERI_Fx2z_S_Fx2z_S_b;
  Double I_ERI_F3y_Px_Fx2z_S_b = I_ERI_Gx3y_S_Fx2z_S_b+ABX*I_ERI_F3y_S_Fx2z_S_b;
  Double I_ERI_F2yz_Px_Fx2z_S_b = I_ERI_Gx2yz_S_Fx2z_S_b+ABX*I_ERI_F2yz_S_Fx2z_S_b;
  Double I_ERI_Fy2z_Px_Fx2z_S_b = I_ERI_Gxy2z_S_Fx2z_S_b+ABX*I_ERI_Fy2z_S_Fx2z_S_b;
  Double I_ERI_F3z_Px_Fx2z_S_b = I_ERI_Gx3z_S_Fx2z_S_b+ABX*I_ERI_F3z_S_Fx2z_S_b;
  Double I_ERI_F2xy_Py_Fx2z_S_b = I_ERI_G2x2y_S_Fx2z_S_b+ABY*I_ERI_F2xy_S_Fx2z_S_b;
  Double I_ERI_F2xz_Py_Fx2z_S_b = I_ERI_G2xyz_S_Fx2z_S_b+ABY*I_ERI_F2xz_S_Fx2z_S_b;
  Double I_ERI_Fx2y_Py_Fx2z_S_b = I_ERI_Gx3y_S_Fx2z_S_b+ABY*I_ERI_Fx2y_S_Fx2z_S_b;
  Double I_ERI_Fxyz_Py_Fx2z_S_b = I_ERI_Gx2yz_S_Fx2z_S_b+ABY*I_ERI_Fxyz_S_Fx2z_S_b;
  Double I_ERI_Fx2z_Py_Fx2z_S_b = I_ERI_Gxy2z_S_Fx2z_S_b+ABY*I_ERI_Fx2z_S_Fx2z_S_b;
  Double I_ERI_F3y_Py_Fx2z_S_b = I_ERI_G4y_S_Fx2z_S_b+ABY*I_ERI_F3y_S_Fx2z_S_b;
  Double I_ERI_F2yz_Py_Fx2z_S_b = I_ERI_G3yz_S_Fx2z_S_b+ABY*I_ERI_F2yz_S_Fx2z_S_b;
  Double I_ERI_Fy2z_Py_Fx2z_S_b = I_ERI_G2y2z_S_Fx2z_S_b+ABY*I_ERI_Fy2z_S_Fx2z_S_b;
  Double I_ERI_F3z_Py_Fx2z_S_b = I_ERI_Gy3z_S_Fx2z_S_b+ABY*I_ERI_F3z_S_Fx2z_S_b;
  Double I_ERI_F2xz_Pz_Fx2z_S_b = I_ERI_G2x2z_S_Fx2z_S_b+ABZ*I_ERI_F2xz_S_Fx2z_S_b;
  Double I_ERI_Fxyz_Pz_Fx2z_S_b = I_ERI_Gxy2z_S_Fx2z_S_b+ABZ*I_ERI_Fxyz_S_Fx2z_S_b;
  Double I_ERI_Fx2z_Pz_Fx2z_S_b = I_ERI_Gx3z_S_Fx2z_S_b+ABZ*I_ERI_Fx2z_S_Fx2z_S_b;
  Double I_ERI_F2yz_Pz_Fx2z_S_b = I_ERI_G2y2z_S_Fx2z_S_b+ABZ*I_ERI_F2yz_S_Fx2z_S_b;
  Double I_ERI_Fy2z_Pz_Fx2z_S_b = I_ERI_Gy3z_S_Fx2z_S_b+ABZ*I_ERI_Fy2z_S_Fx2z_S_b;
  Double I_ERI_F3z_Pz_Fx2z_S_b = I_ERI_G4z_S_Fx2z_S_b+ABZ*I_ERI_F3z_S_Fx2z_S_b;
  Double I_ERI_F3x_Px_F3y_S_b = I_ERI_G4x_S_F3y_S_b+ABX*I_ERI_F3x_S_F3y_S_b;
  Double I_ERI_F2xy_Px_F3y_S_b = I_ERI_G3xy_S_F3y_S_b+ABX*I_ERI_F2xy_S_F3y_S_b;
  Double I_ERI_F2xz_Px_F3y_S_b = I_ERI_G3xz_S_F3y_S_b+ABX*I_ERI_F2xz_S_F3y_S_b;
  Double I_ERI_Fx2y_Px_F3y_S_b = I_ERI_G2x2y_S_F3y_S_b+ABX*I_ERI_Fx2y_S_F3y_S_b;
  Double I_ERI_Fxyz_Px_F3y_S_b = I_ERI_G2xyz_S_F3y_S_b+ABX*I_ERI_Fxyz_S_F3y_S_b;
  Double I_ERI_Fx2z_Px_F3y_S_b = I_ERI_G2x2z_S_F3y_S_b+ABX*I_ERI_Fx2z_S_F3y_S_b;
  Double I_ERI_F3y_Px_F3y_S_b = I_ERI_Gx3y_S_F3y_S_b+ABX*I_ERI_F3y_S_F3y_S_b;
  Double I_ERI_F2yz_Px_F3y_S_b = I_ERI_Gx2yz_S_F3y_S_b+ABX*I_ERI_F2yz_S_F3y_S_b;
  Double I_ERI_Fy2z_Px_F3y_S_b = I_ERI_Gxy2z_S_F3y_S_b+ABX*I_ERI_Fy2z_S_F3y_S_b;
  Double I_ERI_F3z_Px_F3y_S_b = I_ERI_Gx3z_S_F3y_S_b+ABX*I_ERI_F3z_S_F3y_S_b;
  Double I_ERI_F2xy_Py_F3y_S_b = I_ERI_G2x2y_S_F3y_S_b+ABY*I_ERI_F2xy_S_F3y_S_b;
  Double I_ERI_F2xz_Py_F3y_S_b = I_ERI_G2xyz_S_F3y_S_b+ABY*I_ERI_F2xz_S_F3y_S_b;
  Double I_ERI_Fx2y_Py_F3y_S_b = I_ERI_Gx3y_S_F3y_S_b+ABY*I_ERI_Fx2y_S_F3y_S_b;
  Double I_ERI_Fxyz_Py_F3y_S_b = I_ERI_Gx2yz_S_F3y_S_b+ABY*I_ERI_Fxyz_S_F3y_S_b;
  Double I_ERI_Fx2z_Py_F3y_S_b = I_ERI_Gxy2z_S_F3y_S_b+ABY*I_ERI_Fx2z_S_F3y_S_b;
  Double I_ERI_F3y_Py_F3y_S_b = I_ERI_G4y_S_F3y_S_b+ABY*I_ERI_F3y_S_F3y_S_b;
  Double I_ERI_F2yz_Py_F3y_S_b = I_ERI_G3yz_S_F3y_S_b+ABY*I_ERI_F2yz_S_F3y_S_b;
  Double I_ERI_Fy2z_Py_F3y_S_b = I_ERI_G2y2z_S_F3y_S_b+ABY*I_ERI_Fy2z_S_F3y_S_b;
  Double I_ERI_F3z_Py_F3y_S_b = I_ERI_Gy3z_S_F3y_S_b+ABY*I_ERI_F3z_S_F3y_S_b;
  Double I_ERI_F2xz_Pz_F3y_S_b = I_ERI_G2x2z_S_F3y_S_b+ABZ*I_ERI_F2xz_S_F3y_S_b;
  Double I_ERI_Fxyz_Pz_F3y_S_b = I_ERI_Gxy2z_S_F3y_S_b+ABZ*I_ERI_Fxyz_S_F3y_S_b;
  Double I_ERI_Fx2z_Pz_F3y_S_b = I_ERI_Gx3z_S_F3y_S_b+ABZ*I_ERI_Fx2z_S_F3y_S_b;
  Double I_ERI_F2yz_Pz_F3y_S_b = I_ERI_G2y2z_S_F3y_S_b+ABZ*I_ERI_F2yz_S_F3y_S_b;
  Double I_ERI_Fy2z_Pz_F3y_S_b = I_ERI_Gy3z_S_F3y_S_b+ABZ*I_ERI_Fy2z_S_F3y_S_b;
  Double I_ERI_F3z_Pz_F3y_S_b = I_ERI_G4z_S_F3y_S_b+ABZ*I_ERI_F3z_S_F3y_S_b;
  Double I_ERI_F3x_Px_F2yz_S_b = I_ERI_G4x_S_F2yz_S_b+ABX*I_ERI_F3x_S_F2yz_S_b;
  Double I_ERI_F2xy_Px_F2yz_S_b = I_ERI_G3xy_S_F2yz_S_b+ABX*I_ERI_F2xy_S_F2yz_S_b;
  Double I_ERI_F2xz_Px_F2yz_S_b = I_ERI_G3xz_S_F2yz_S_b+ABX*I_ERI_F2xz_S_F2yz_S_b;
  Double I_ERI_Fx2y_Px_F2yz_S_b = I_ERI_G2x2y_S_F2yz_S_b+ABX*I_ERI_Fx2y_S_F2yz_S_b;
  Double I_ERI_Fxyz_Px_F2yz_S_b = I_ERI_G2xyz_S_F2yz_S_b+ABX*I_ERI_Fxyz_S_F2yz_S_b;
  Double I_ERI_Fx2z_Px_F2yz_S_b = I_ERI_G2x2z_S_F2yz_S_b+ABX*I_ERI_Fx2z_S_F2yz_S_b;
  Double I_ERI_F3y_Px_F2yz_S_b = I_ERI_Gx3y_S_F2yz_S_b+ABX*I_ERI_F3y_S_F2yz_S_b;
  Double I_ERI_F2yz_Px_F2yz_S_b = I_ERI_Gx2yz_S_F2yz_S_b+ABX*I_ERI_F2yz_S_F2yz_S_b;
  Double I_ERI_Fy2z_Px_F2yz_S_b = I_ERI_Gxy2z_S_F2yz_S_b+ABX*I_ERI_Fy2z_S_F2yz_S_b;
  Double I_ERI_F3z_Px_F2yz_S_b = I_ERI_Gx3z_S_F2yz_S_b+ABX*I_ERI_F3z_S_F2yz_S_b;
  Double I_ERI_F2xy_Py_F2yz_S_b = I_ERI_G2x2y_S_F2yz_S_b+ABY*I_ERI_F2xy_S_F2yz_S_b;
  Double I_ERI_F2xz_Py_F2yz_S_b = I_ERI_G2xyz_S_F2yz_S_b+ABY*I_ERI_F2xz_S_F2yz_S_b;
  Double I_ERI_Fx2y_Py_F2yz_S_b = I_ERI_Gx3y_S_F2yz_S_b+ABY*I_ERI_Fx2y_S_F2yz_S_b;
  Double I_ERI_Fxyz_Py_F2yz_S_b = I_ERI_Gx2yz_S_F2yz_S_b+ABY*I_ERI_Fxyz_S_F2yz_S_b;
  Double I_ERI_Fx2z_Py_F2yz_S_b = I_ERI_Gxy2z_S_F2yz_S_b+ABY*I_ERI_Fx2z_S_F2yz_S_b;
  Double I_ERI_F3y_Py_F2yz_S_b = I_ERI_G4y_S_F2yz_S_b+ABY*I_ERI_F3y_S_F2yz_S_b;
  Double I_ERI_F2yz_Py_F2yz_S_b = I_ERI_G3yz_S_F2yz_S_b+ABY*I_ERI_F2yz_S_F2yz_S_b;
  Double I_ERI_Fy2z_Py_F2yz_S_b = I_ERI_G2y2z_S_F2yz_S_b+ABY*I_ERI_Fy2z_S_F2yz_S_b;
  Double I_ERI_F3z_Py_F2yz_S_b = I_ERI_Gy3z_S_F2yz_S_b+ABY*I_ERI_F3z_S_F2yz_S_b;
  Double I_ERI_F2xz_Pz_F2yz_S_b = I_ERI_G2x2z_S_F2yz_S_b+ABZ*I_ERI_F2xz_S_F2yz_S_b;
  Double I_ERI_Fxyz_Pz_F2yz_S_b = I_ERI_Gxy2z_S_F2yz_S_b+ABZ*I_ERI_Fxyz_S_F2yz_S_b;
  Double I_ERI_Fx2z_Pz_F2yz_S_b = I_ERI_Gx3z_S_F2yz_S_b+ABZ*I_ERI_Fx2z_S_F2yz_S_b;
  Double I_ERI_F2yz_Pz_F2yz_S_b = I_ERI_G2y2z_S_F2yz_S_b+ABZ*I_ERI_F2yz_S_F2yz_S_b;
  Double I_ERI_Fy2z_Pz_F2yz_S_b = I_ERI_Gy3z_S_F2yz_S_b+ABZ*I_ERI_Fy2z_S_F2yz_S_b;
  Double I_ERI_F3z_Pz_F2yz_S_b = I_ERI_G4z_S_F2yz_S_b+ABZ*I_ERI_F3z_S_F2yz_S_b;
  Double I_ERI_F3x_Px_Fy2z_S_b = I_ERI_G4x_S_Fy2z_S_b+ABX*I_ERI_F3x_S_Fy2z_S_b;
  Double I_ERI_F2xy_Px_Fy2z_S_b = I_ERI_G3xy_S_Fy2z_S_b+ABX*I_ERI_F2xy_S_Fy2z_S_b;
  Double I_ERI_F2xz_Px_Fy2z_S_b = I_ERI_G3xz_S_Fy2z_S_b+ABX*I_ERI_F2xz_S_Fy2z_S_b;
  Double I_ERI_Fx2y_Px_Fy2z_S_b = I_ERI_G2x2y_S_Fy2z_S_b+ABX*I_ERI_Fx2y_S_Fy2z_S_b;
  Double I_ERI_Fxyz_Px_Fy2z_S_b = I_ERI_G2xyz_S_Fy2z_S_b+ABX*I_ERI_Fxyz_S_Fy2z_S_b;
  Double I_ERI_Fx2z_Px_Fy2z_S_b = I_ERI_G2x2z_S_Fy2z_S_b+ABX*I_ERI_Fx2z_S_Fy2z_S_b;
  Double I_ERI_F3y_Px_Fy2z_S_b = I_ERI_Gx3y_S_Fy2z_S_b+ABX*I_ERI_F3y_S_Fy2z_S_b;
  Double I_ERI_F2yz_Px_Fy2z_S_b = I_ERI_Gx2yz_S_Fy2z_S_b+ABX*I_ERI_F2yz_S_Fy2z_S_b;
  Double I_ERI_Fy2z_Px_Fy2z_S_b = I_ERI_Gxy2z_S_Fy2z_S_b+ABX*I_ERI_Fy2z_S_Fy2z_S_b;
  Double I_ERI_F3z_Px_Fy2z_S_b = I_ERI_Gx3z_S_Fy2z_S_b+ABX*I_ERI_F3z_S_Fy2z_S_b;
  Double I_ERI_F2xy_Py_Fy2z_S_b = I_ERI_G2x2y_S_Fy2z_S_b+ABY*I_ERI_F2xy_S_Fy2z_S_b;
  Double I_ERI_F2xz_Py_Fy2z_S_b = I_ERI_G2xyz_S_Fy2z_S_b+ABY*I_ERI_F2xz_S_Fy2z_S_b;
  Double I_ERI_Fx2y_Py_Fy2z_S_b = I_ERI_Gx3y_S_Fy2z_S_b+ABY*I_ERI_Fx2y_S_Fy2z_S_b;
  Double I_ERI_Fxyz_Py_Fy2z_S_b = I_ERI_Gx2yz_S_Fy2z_S_b+ABY*I_ERI_Fxyz_S_Fy2z_S_b;
  Double I_ERI_Fx2z_Py_Fy2z_S_b = I_ERI_Gxy2z_S_Fy2z_S_b+ABY*I_ERI_Fx2z_S_Fy2z_S_b;
  Double I_ERI_F3y_Py_Fy2z_S_b = I_ERI_G4y_S_Fy2z_S_b+ABY*I_ERI_F3y_S_Fy2z_S_b;
  Double I_ERI_F2yz_Py_Fy2z_S_b = I_ERI_G3yz_S_Fy2z_S_b+ABY*I_ERI_F2yz_S_Fy2z_S_b;
  Double I_ERI_Fy2z_Py_Fy2z_S_b = I_ERI_G2y2z_S_Fy2z_S_b+ABY*I_ERI_Fy2z_S_Fy2z_S_b;
  Double I_ERI_F3z_Py_Fy2z_S_b = I_ERI_Gy3z_S_Fy2z_S_b+ABY*I_ERI_F3z_S_Fy2z_S_b;
  Double I_ERI_F2xz_Pz_Fy2z_S_b = I_ERI_G2x2z_S_Fy2z_S_b+ABZ*I_ERI_F2xz_S_Fy2z_S_b;
  Double I_ERI_Fxyz_Pz_Fy2z_S_b = I_ERI_Gxy2z_S_Fy2z_S_b+ABZ*I_ERI_Fxyz_S_Fy2z_S_b;
  Double I_ERI_Fx2z_Pz_Fy2z_S_b = I_ERI_Gx3z_S_Fy2z_S_b+ABZ*I_ERI_Fx2z_S_Fy2z_S_b;
  Double I_ERI_F2yz_Pz_Fy2z_S_b = I_ERI_G2y2z_S_Fy2z_S_b+ABZ*I_ERI_F2yz_S_Fy2z_S_b;
  Double I_ERI_Fy2z_Pz_Fy2z_S_b = I_ERI_Gy3z_S_Fy2z_S_b+ABZ*I_ERI_Fy2z_S_Fy2z_S_b;
  Double I_ERI_F3z_Pz_Fy2z_S_b = I_ERI_G4z_S_Fy2z_S_b+ABZ*I_ERI_F3z_S_Fy2z_S_b;
  Double I_ERI_F3x_Px_F3z_S_b = I_ERI_G4x_S_F3z_S_b+ABX*I_ERI_F3x_S_F3z_S_b;
  Double I_ERI_F2xy_Px_F3z_S_b = I_ERI_G3xy_S_F3z_S_b+ABX*I_ERI_F2xy_S_F3z_S_b;
  Double I_ERI_F2xz_Px_F3z_S_b = I_ERI_G3xz_S_F3z_S_b+ABX*I_ERI_F2xz_S_F3z_S_b;
  Double I_ERI_Fx2y_Px_F3z_S_b = I_ERI_G2x2y_S_F3z_S_b+ABX*I_ERI_Fx2y_S_F3z_S_b;
  Double I_ERI_Fxyz_Px_F3z_S_b = I_ERI_G2xyz_S_F3z_S_b+ABX*I_ERI_Fxyz_S_F3z_S_b;
  Double I_ERI_Fx2z_Px_F3z_S_b = I_ERI_G2x2z_S_F3z_S_b+ABX*I_ERI_Fx2z_S_F3z_S_b;
  Double I_ERI_F3y_Px_F3z_S_b = I_ERI_Gx3y_S_F3z_S_b+ABX*I_ERI_F3y_S_F3z_S_b;
  Double I_ERI_F2yz_Px_F3z_S_b = I_ERI_Gx2yz_S_F3z_S_b+ABX*I_ERI_F2yz_S_F3z_S_b;
  Double I_ERI_Fy2z_Px_F3z_S_b = I_ERI_Gxy2z_S_F3z_S_b+ABX*I_ERI_Fy2z_S_F3z_S_b;
  Double I_ERI_F3z_Px_F3z_S_b = I_ERI_Gx3z_S_F3z_S_b+ABX*I_ERI_F3z_S_F3z_S_b;
  Double I_ERI_F2xy_Py_F3z_S_b = I_ERI_G2x2y_S_F3z_S_b+ABY*I_ERI_F2xy_S_F3z_S_b;
  Double I_ERI_F2xz_Py_F3z_S_b = I_ERI_G2xyz_S_F3z_S_b+ABY*I_ERI_F2xz_S_F3z_S_b;
  Double I_ERI_Fx2y_Py_F3z_S_b = I_ERI_Gx3y_S_F3z_S_b+ABY*I_ERI_Fx2y_S_F3z_S_b;
  Double I_ERI_Fxyz_Py_F3z_S_b = I_ERI_Gx2yz_S_F3z_S_b+ABY*I_ERI_Fxyz_S_F3z_S_b;
  Double I_ERI_Fx2z_Py_F3z_S_b = I_ERI_Gxy2z_S_F3z_S_b+ABY*I_ERI_Fx2z_S_F3z_S_b;
  Double I_ERI_F3y_Py_F3z_S_b = I_ERI_G4y_S_F3z_S_b+ABY*I_ERI_F3y_S_F3z_S_b;
  Double I_ERI_F2yz_Py_F3z_S_b = I_ERI_G3yz_S_F3z_S_b+ABY*I_ERI_F2yz_S_F3z_S_b;
  Double I_ERI_Fy2z_Py_F3z_S_b = I_ERI_G2y2z_S_F3z_S_b+ABY*I_ERI_Fy2z_S_F3z_S_b;
  Double I_ERI_F3z_Py_F3z_S_b = I_ERI_Gy3z_S_F3z_S_b+ABY*I_ERI_F3z_S_F3z_S_b;
  Double I_ERI_F2xz_Pz_F3z_S_b = I_ERI_G2x2z_S_F3z_S_b+ABZ*I_ERI_F2xz_S_F3z_S_b;
  Double I_ERI_Fxyz_Pz_F3z_S_b = I_ERI_Gxy2z_S_F3z_S_b+ABZ*I_ERI_Fxyz_S_F3z_S_b;
  Double I_ERI_Fx2z_Pz_F3z_S_b = I_ERI_Gx3z_S_F3z_S_b+ABZ*I_ERI_Fx2z_S_F3z_S_b;
  Double I_ERI_F2yz_Pz_F3z_S_b = I_ERI_G2y2z_S_F3z_S_b+ABZ*I_ERI_F2yz_S_F3z_S_b;
  Double I_ERI_Fy2z_Pz_F3z_S_b = I_ERI_Gy3z_S_F3z_S_b+ABZ*I_ERI_Fy2z_S_F3z_S_b;
  Double I_ERI_F3z_Pz_F3z_S_b = I_ERI_G4z_S_F3z_S_b+ABZ*I_ERI_F3z_S_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_F_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_b
   * RHS shell quartet name: SQ_ERI_D_P_F_S_b
   ************************************************************/
  Double I_ERI_D2x_D2x_F3x_S_b = I_ERI_F3x_Px_F3x_S_b+ABX*I_ERI_D2x_Px_F3x_S_b;
  Double I_ERI_Dxy_D2x_F3x_S_b = I_ERI_F2xy_Px_F3x_S_b+ABX*I_ERI_Dxy_Px_F3x_S_b;
  Double I_ERI_Dxz_D2x_F3x_S_b = I_ERI_F2xz_Px_F3x_S_b+ABX*I_ERI_Dxz_Px_F3x_S_b;
  Double I_ERI_D2y_D2x_F3x_S_b = I_ERI_Fx2y_Px_F3x_S_b+ABX*I_ERI_D2y_Px_F3x_S_b;
  Double I_ERI_Dyz_D2x_F3x_S_b = I_ERI_Fxyz_Px_F3x_S_b+ABX*I_ERI_Dyz_Px_F3x_S_b;
  Double I_ERI_D2z_D2x_F3x_S_b = I_ERI_Fx2z_Px_F3x_S_b+ABX*I_ERI_D2z_Px_F3x_S_b;
  Double I_ERI_D2x_Dxy_F3x_S_b = I_ERI_F2xy_Px_F3x_S_b+ABY*I_ERI_D2x_Px_F3x_S_b;
  Double I_ERI_Dxy_Dxy_F3x_S_b = I_ERI_Fx2y_Px_F3x_S_b+ABY*I_ERI_Dxy_Px_F3x_S_b;
  Double I_ERI_Dxz_Dxy_F3x_S_b = I_ERI_Fxyz_Px_F3x_S_b+ABY*I_ERI_Dxz_Px_F3x_S_b;
  Double I_ERI_D2y_Dxy_F3x_S_b = I_ERI_F3y_Px_F3x_S_b+ABY*I_ERI_D2y_Px_F3x_S_b;
  Double I_ERI_Dyz_Dxy_F3x_S_b = I_ERI_F2yz_Px_F3x_S_b+ABY*I_ERI_Dyz_Px_F3x_S_b;
  Double I_ERI_D2z_Dxy_F3x_S_b = I_ERI_Fy2z_Px_F3x_S_b+ABY*I_ERI_D2z_Px_F3x_S_b;
  Double I_ERI_D2x_Dxz_F3x_S_b = I_ERI_F2xz_Px_F3x_S_b+ABZ*I_ERI_D2x_Px_F3x_S_b;
  Double I_ERI_Dxy_Dxz_F3x_S_b = I_ERI_Fxyz_Px_F3x_S_b+ABZ*I_ERI_Dxy_Px_F3x_S_b;
  Double I_ERI_Dxz_Dxz_F3x_S_b = I_ERI_Fx2z_Px_F3x_S_b+ABZ*I_ERI_Dxz_Px_F3x_S_b;
  Double I_ERI_D2y_Dxz_F3x_S_b = I_ERI_F2yz_Px_F3x_S_b+ABZ*I_ERI_D2y_Px_F3x_S_b;
  Double I_ERI_Dyz_Dxz_F3x_S_b = I_ERI_Fy2z_Px_F3x_S_b+ABZ*I_ERI_Dyz_Px_F3x_S_b;
  Double I_ERI_D2z_Dxz_F3x_S_b = I_ERI_F3z_Px_F3x_S_b+ABZ*I_ERI_D2z_Px_F3x_S_b;
  Double I_ERI_D2x_D2y_F3x_S_b = I_ERI_F2xy_Py_F3x_S_b+ABY*I_ERI_D2x_Py_F3x_S_b;
  Double I_ERI_Dxy_D2y_F3x_S_b = I_ERI_Fx2y_Py_F3x_S_b+ABY*I_ERI_Dxy_Py_F3x_S_b;
  Double I_ERI_Dxz_D2y_F3x_S_b = I_ERI_Fxyz_Py_F3x_S_b+ABY*I_ERI_Dxz_Py_F3x_S_b;
  Double I_ERI_D2y_D2y_F3x_S_b = I_ERI_F3y_Py_F3x_S_b+ABY*I_ERI_D2y_Py_F3x_S_b;
  Double I_ERI_Dyz_D2y_F3x_S_b = I_ERI_F2yz_Py_F3x_S_b+ABY*I_ERI_Dyz_Py_F3x_S_b;
  Double I_ERI_D2z_D2y_F3x_S_b = I_ERI_Fy2z_Py_F3x_S_b+ABY*I_ERI_D2z_Py_F3x_S_b;
  Double I_ERI_D2x_Dyz_F3x_S_b = I_ERI_F2xz_Py_F3x_S_b+ABZ*I_ERI_D2x_Py_F3x_S_b;
  Double I_ERI_Dxy_Dyz_F3x_S_b = I_ERI_Fxyz_Py_F3x_S_b+ABZ*I_ERI_Dxy_Py_F3x_S_b;
  Double I_ERI_Dxz_Dyz_F3x_S_b = I_ERI_Fx2z_Py_F3x_S_b+ABZ*I_ERI_Dxz_Py_F3x_S_b;
  Double I_ERI_D2y_Dyz_F3x_S_b = I_ERI_F2yz_Py_F3x_S_b+ABZ*I_ERI_D2y_Py_F3x_S_b;
  Double I_ERI_Dyz_Dyz_F3x_S_b = I_ERI_Fy2z_Py_F3x_S_b+ABZ*I_ERI_Dyz_Py_F3x_S_b;
  Double I_ERI_D2z_Dyz_F3x_S_b = I_ERI_F3z_Py_F3x_S_b+ABZ*I_ERI_D2z_Py_F3x_S_b;
  Double I_ERI_D2x_D2z_F3x_S_b = I_ERI_F2xz_Pz_F3x_S_b+ABZ*I_ERI_D2x_Pz_F3x_S_b;
  Double I_ERI_Dxy_D2z_F3x_S_b = I_ERI_Fxyz_Pz_F3x_S_b+ABZ*I_ERI_Dxy_Pz_F3x_S_b;
  Double I_ERI_Dxz_D2z_F3x_S_b = I_ERI_Fx2z_Pz_F3x_S_b+ABZ*I_ERI_Dxz_Pz_F3x_S_b;
  Double I_ERI_D2y_D2z_F3x_S_b = I_ERI_F2yz_Pz_F3x_S_b+ABZ*I_ERI_D2y_Pz_F3x_S_b;
  Double I_ERI_Dyz_D2z_F3x_S_b = I_ERI_Fy2z_Pz_F3x_S_b+ABZ*I_ERI_Dyz_Pz_F3x_S_b;
  Double I_ERI_D2z_D2z_F3x_S_b = I_ERI_F3z_Pz_F3x_S_b+ABZ*I_ERI_D2z_Pz_F3x_S_b;
  Double I_ERI_D2x_D2x_F2xy_S_b = I_ERI_F3x_Px_F2xy_S_b+ABX*I_ERI_D2x_Px_F2xy_S_b;
  Double I_ERI_Dxy_D2x_F2xy_S_b = I_ERI_F2xy_Px_F2xy_S_b+ABX*I_ERI_Dxy_Px_F2xy_S_b;
  Double I_ERI_Dxz_D2x_F2xy_S_b = I_ERI_F2xz_Px_F2xy_S_b+ABX*I_ERI_Dxz_Px_F2xy_S_b;
  Double I_ERI_D2y_D2x_F2xy_S_b = I_ERI_Fx2y_Px_F2xy_S_b+ABX*I_ERI_D2y_Px_F2xy_S_b;
  Double I_ERI_Dyz_D2x_F2xy_S_b = I_ERI_Fxyz_Px_F2xy_S_b+ABX*I_ERI_Dyz_Px_F2xy_S_b;
  Double I_ERI_D2z_D2x_F2xy_S_b = I_ERI_Fx2z_Px_F2xy_S_b+ABX*I_ERI_D2z_Px_F2xy_S_b;
  Double I_ERI_D2x_Dxy_F2xy_S_b = I_ERI_F2xy_Px_F2xy_S_b+ABY*I_ERI_D2x_Px_F2xy_S_b;
  Double I_ERI_Dxy_Dxy_F2xy_S_b = I_ERI_Fx2y_Px_F2xy_S_b+ABY*I_ERI_Dxy_Px_F2xy_S_b;
  Double I_ERI_Dxz_Dxy_F2xy_S_b = I_ERI_Fxyz_Px_F2xy_S_b+ABY*I_ERI_Dxz_Px_F2xy_S_b;
  Double I_ERI_D2y_Dxy_F2xy_S_b = I_ERI_F3y_Px_F2xy_S_b+ABY*I_ERI_D2y_Px_F2xy_S_b;
  Double I_ERI_Dyz_Dxy_F2xy_S_b = I_ERI_F2yz_Px_F2xy_S_b+ABY*I_ERI_Dyz_Px_F2xy_S_b;
  Double I_ERI_D2z_Dxy_F2xy_S_b = I_ERI_Fy2z_Px_F2xy_S_b+ABY*I_ERI_D2z_Px_F2xy_S_b;
  Double I_ERI_D2x_Dxz_F2xy_S_b = I_ERI_F2xz_Px_F2xy_S_b+ABZ*I_ERI_D2x_Px_F2xy_S_b;
  Double I_ERI_Dxy_Dxz_F2xy_S_b = I_ERI_Fxyz_Px_F2xy_S_b+ABZ*I_ERI_Dxy_Px_F2xy_S_b;
  Double I_ERI_Dxz_Dxz_F2xy_S_b = I_ERI_Fx2z_Px_F2xy_S_b+ABZ*I_ERI_Dxz_Px_F2xy_S_b;
  Double I_ERI_D2y_Dxz_F2xy_S_b = I_ERI_F2yz_Px_F2xy_S_b+ABZ*I_ERI_D2y_Px_F2xy_S_b;
  Double I_ERI_Dyz_Dxz_F2xy_S_b = I_ERI_Fy2z_Px_F2xy_S_b+ABZ*I_ERI_Dyz_Px_F2xy_S_b;
  Double I_ERI_D2z_Dxz_F2xy_S_b = I_ERI_F3z_Px_F2xy_S_b+ABZ*I_ERI_D2z_Px_F2xy_S_b;
  Double I_ERI_D2x_D2y_F2xy_S_b = I_ERI_F2xy_Py_F2xy_S_b+ABY*I_ERI_D2x_Py_F2xy_S_b;
  Double I_ERI_Dxy_D2y_F2xy_S_b = I_ERI_Fx2y_Py_F2xy_S_b+ABY*I_ERI_Dxy_Py_F2xy_S_b;
  Double I_ERI_Dxz_D2y_F2xy_S_b = I_ERI_Fxyz_Py_F2xy_S_b+ABY*I_ERI_Dxz_Py_F2xy_S_b;
  Double I_ERI_D2y_D2y_F2xy_S_b = I_ERI_F3y_Py_F2xy_S_b+ABY*I_ERI_D2y_Py_F2xy_S_b;
  Double I_ERI_Dyz_D2y_F2xy_S_b = I_ERI_F2yz_Py_F2xy_S_b+ABY*I_ERI_Dyz_Py_F2xy_S_b;
  Double I_ERI_D2z_D2y_F2xy_S_b = I_ERI_Fy2z_Py_F2xy_S_b+ABY*I_ERI_D2z_Py_F2xy_S_b;
  Double I_ERI_D2x_Dyz_F2xy_S_b = I_ERI_F2xz_Py_F2xy_S_b+ABZ*I_ERI_D2x_Py_F2xy_S_b;
  Double I_ERI_Dxy_Dyz_F2xy_S_b = I_ERI_Fxyz_Py_F2xy_S_b+ABZ*I_ERI_Dxy_Py_F2xy_S_b;
  Double I_ERI_Dxz_Dyz_F2xy_S_b = I_ERI_Fx2z_Py_F2xy_S_b+ABZ*I_ERI_Dxz_Py_F2xy_S_b;
  Double I_ERI_D2y_Dyz_F2xy_S_b = I_ERI_F2yz_Py_F2xy_S_b+ABZ*I_ERI_D2y_Py_F2xy_S_b;
  Double I_ERI_Dyz_Dyz_F2xy_S_b = I_ERI_Fy2z_Py_F2xy_S_b+ABZ*I_ERI_Dyz_Py_F2xy_S_b;
  Double I_ERI_D2z_Dyz_F2xy_S_b = I_ERI_F3z_Py_F2xy_S_b+ABZ*I_ERI_D2z_Py_F2xy_S_b;
  Double I_ERI_D2x_D2z_F2xy_S_b = I_ERI_F2xz_Pz_F2xy_S_b+ABZ*I_ERI_D2x_Pz_F2xy_S_b;
  Double I_ERI_Dxy_D2z_F2xy_S_b = I_ERI_Fxyz_Pz_F2xy_S_b+ABZ*I_ERI_Dxy_Pz_F2xy_S_b;
  Double I_ERI_Dxz_D2z_F2xy_S_b = I_ERI_Fx2z_Pz_F2xy_S_b+ABZ*I_ERI_Dxz_Pz_F2xy_S_b;
  Double I_ERI_D2y_D2z_F2xy_S_b = I_ERI_F2yz_Pz_F2xy_S_b+ABZ*I_ERI_D2y_Pz_F2xy_S_b;
  Double I_ERI_Dyz_D2z_F2xy_S_b = I_ERI_Fy2z_Pz_F2xy_S_b+ABZ*I_ERI_Dyz_Pz_F2xy_S_b;
  Double I_ERI_D2z_D2z_F2xy_S_b = I_ERI_F3z_Pz_F2xy_S_b+ABZ*I_ERI_D2z_Pz_F2xy_S_b;
  Double I_ERI_D2x_D2x_F2xz_S_b = I_ERI_F3x_Px_F2xz_S_b+ABX*I_ERI_D2x_Px_F2xz_S_b;
  Double I_ERI_Dxy_D2x_F2xz_S_b = I_ERI_F2xy_Px_F2xz_S_b+ABX*I_ERI_Dxy_Px_F2xz_S_b;
  Double I_ERI_Dxz_D2x_F2xz_S_b = I_ERI_F2xz_Px_F2xz_S_b+ABX*I_ERI_Dxz_Px_F2xz_S_b;
  Double I_ERI_D2y_D2x_F2xz_S_b = I_ERI_Fx2y_Px_F2xz_S_b+ABX*I_ERI_D2y_Px_F2xz_S_b;
  Double I_ERI_Dyz_D2x_F2xz_S_b = I_ERI_Fxyz_Px_F2xz_S_b+ABX*I_ERI_Dyz_Px_F2xz_S_b;
  Double I_ERI_D2z_D2x_F2xz_S_b = I_ERI_Fx2z_Px_F2xz_S_b+ABX*I_ERI_D2z_Px_F2xz_S_b;
  Double I_ERI_D2x_Dxy_F2xz_S_b = I_ERI_F2xy_Px_F2xz_S_b+ABY*I_ERI_D2x_Px_F2xz_S_b;
  Double I_ERI_Dxy_Dxy_F2xz_S_b = I_ERI_Fx2y_Px_F2xz_S_b+ABY*I_ERI_Dxy_Px_F2xz_S_b;
  Double I_ERI_Dxz_Dxy_F2xz_S_b = I_ERI_Fxyz_Px_F2xz_S_b+ABY*I_ERI_Dxz_Px_F2xz_S_b;
  Double I_ERI_D2y_Dxy_F2xz_S_b = I_ERI_F3y_Px_F2xz_S_b+ABY*I_ERI_D2y_Px_F2xz_S_b;
  Double I_ERI_Dyz_Dxy_F2xz_S_b = I_ERI_F2yz_Px_F2xz_S_b+ABY*I_ERI_Dyz_Px_F2xz_S_b;
  Double I_ERI_D2z_Dxy_F2xz_S_b = I_ERI_Fy2z_Px_F2xz_S_b+ABY*I_ERI_D2z_Px_F2xz_S_b;
  Double I_ERI_D2x_Dxz_F2xz_S_b = I_ERI_F2xz_Px_F2xz_S_b+ABZ*I_ERI_D2x_Px_F2xz_S_b;
  Double I_ERI_Dxy_Dxz_F2xz_S_b = I_ERI_Fxyz_Px_F2xz_S_b+ABZ*I_ERI_Dxy_Px_F2xz_S_b;
  Double I_ERI_Dxz_Dxz_F2xz_S_b = I_ERI_Fx2z_Px_F2xz_S_b+ABZ*I_ERI_Dxz_Px_F2xz_S_b;
  Double I_ERI_D2y_Dxz_F2xz_S_b = I_ERI_F2yz_Px_F2xz_S_b+ABZ*I_ERI_D2y_Px_F2xz_S_b;
  Double I_ERI_Dyz_Dxz_F2xz_S_b = I_ERI_Fy2z_Px_F2xz_S_b+ABZ*I_ERI_Dyz_Px_F2xz_S_b;
  Double I_ERI_D2z_Dxz_F2xz_S_b = I_ERI_F3z_Px_F2xz_S_b+ABZ*I_ERI_D2z_Px_F2xz_S_b;
  Double I_ERI_D2x_D2y_F2xz_S_b = I_ERI_F2xy_Py_F2xz_S_b+ABY*I_ERI_D2x_Py_F2xz_S_b;
  Double I_ERI_Dxy_D2y_F2xz_S_b = I_ERI_Fx2y_Py_F2xz_S_b+ABY*I_ERI_Dxy_Py_F2xz_S_b;
  Double I_ERI_Dxz_D2y_F2xz_S_b = I_ERI_Fxyz_Py_F2xz_S_b+ABY*I_ERI_Dxz_Py_F2xz_S_b;
  Double I_ERI_D2y_D2y_F2xz_S_b = I_ERI_F3y_Py_F2xz_S_b+ABY*I_ERI_D2y_Py_F2xz_S_b;
  Double I_ERI_Dyz_D2y_F2xz_S_b = I_ERI_F2yz_Py_F2xz_S_b+ABY*I_ERI_Dyz_Py_F2xz_S_b;
  Double I_ERI_D2z_D2y_F2xz_S_b = I_ERI_Fy2z_Py_F2xz_S_b+ABY*I_ERI_D2z_Py_F2xz_S_b;
  Double I_ERI_D2x_Dyz_F2xz_S_b = I_ERI_F2xz_Py_F2xz_S_b+ABZ*I_ERI_D2x_Py_F2xz_S_b;
  Double I_ERI_Dxy_Dyz_F2xz_S_b = I_ERI_Fxyz_Py_F2xz_S_b+ABZ*I_ERI_Dxy_Py_F2xz_S_b;
  Double I_ERI_Dxz_Dyz_F2xz_S_b = I_ERI_Fx2z_Py_F2xz_S_b+ABZ*I_ERI_Dxz_Py_F2xz_S_b;
  Double I_ERI_D2y_Dyz_F2xz_S_b = I_ERI_F2yz_Py_F2xz_S_b+ABZ*I_ERI_D2y_Py_F2xz_S_b;
  Double I_ERI_Dyz_Dyz_F2xz_S_b = I_ERI_Fy2z_Py_F2xz_S_b+ABZ*I_ERI_Dyz_Py_F2xz_S_b;
  Double I_ERI_D2z_Dyz_F2xz_S_b = I_ERI_F3z_Py_F2xz_S_b+ABZ*I_ERI_D2z_Py_F2xz_S_b;
  Double I_ERI_D2x_D2z_F2xz_S_b = I_ERI_F2xz_Pz_F2xz_S_b+ABZ*I_ERI_D2x_Pz_F2xz_S_b;
  Double I_ERI_Dxy_D2z_F2xz_S_b = I_ERI_Fxyz_Pz_F2xz_S_b+ABZ*I_ERI_Dxy_Pz_F2xz_S_b;
  Double I_ERI_Dxz_D2z_F2xz_S_b = I_ERI_Fx2z_Pz_F2xz_S_b+ABZ*I_ERI_Dxz_Pz_F2xz_S_b;
  Double I_ERI_D2y_D2z_F2xz_S_b = I_ERI_F2yz_Pz_F2xz_S_b+ABZ*I_ERI_D2y_Pz_F2xz_S_b;
  Double I_ERI_Dyz_D2z_F2xz_S_b = I_ERI_Fy2z_Pz_F2xz_S_b+ABZ*I_ERI_Dyz_Pz_F2xz_S_b;
  Double I_ERI_D2z_D2z_F2xz_S_b = I_ERI_F3z_Pz_F2xz_S_b+ABZ*I_ERI_D2z_Pz_F2xz_S_b;
  Double I_ERI_D2x_D2x_Fx2y_S_b = I_ERI_F3x_Px_Fx2y_S_b+ABX*I_ERI_D2x_Px_Fx2y_S_b;
  Double I_ERI_Dxy_D2x_Fx2y_S_b = I_ERI_F2xy_Px_Fx2y_S_b+ABX*I_ERI_Dxy_Px_Fx2y_S_b;
  Double I_ERI_Dxz_D2x_Fx2y_S_b = I_ERI_F2xz_Px_Fx2y_S_b+ABX*I_ERI_Dxz_Px_Fx2y_S_b;
  Double I_ERI_D2y_D2x_Fx2y_S_b = I_ERI_Fx2y_Px_Fx2y_S_b+ABX*I_ERI_D2y_Px_Fx2y_S_b;
  Double I_ERI_Dyz_D2x_Fx2y_S_b = I_ERI_Fxyz_Px_Fx2y_S_b+ABX*I_ERI_Dyz_Px_Fx2y_S_b;
  Double I_ERI_D2z_D2x_Fx2y_S_b = I_ERI_Fx2z_Px_Fx2y_S_b+ABX*I_ERI_D2z_Px_Fx2y_S_b;
  Double I_ERI_D2x_Dxy_Fx2y_S_b = I_ERI_F2xy_Px_Fx2y_S_b+ABY*I_ERI_D2x_Px_Fx2y_S_b;
  Double I_ERI_Dxy_Dxy_Fx2y_S_b = I_ERI_Fx2y_Px_Fx2y_S_b+ABY*I_ERI_Dxy_Px_Fx2y_S_b;
  Double I_ERI_Dxz_Dxy_Fx2y_S_b = I_ERI_Fxyz_Px_Fx2y_S_b+ABY*I_ERI_Dxz_Px_Fx2y_S_b;
  Double I_ERI_D2y_Dxy_Fx2y_S_b = I_ERI_F3y_Px_Fx2y_S_b+ABY*I_ERI_D2y_Px_Fx2y_S_b;
  Double I_ERI_Dyz_Dxy_Fx2y_S_b = I_ERI_F2yz_Px_Fx2y_S_b+ABY*I_ERI_Dyz_Px_Fx2y_S_b;
  Double I_ERI_D2z_Dxy_Fx2y_S_b = I_ERI_Fy2z_Px_Fx2y_S_b+ABY*I_ERI_D2z_Px_Fx2y_S_b;
  Double I_ERI_D2x_Dxz_Fx2y_S_b = I_ERI_F2xz_Px_Fx2y_S_b+ABZ*I_ERI_D2x_Px_Fx2y_S_b;
  Double I_ERI_Dxy_Dxz_Fx2y_S_b = I_ERI_Fxyz_Px_Fx2y_S_b+ABZ*I_ERI_Dxy_Px_Fx2y_S_b;
  Double I_ERI_Dxz_Dxz_Fx2y_S_b = I_ERI_Fx2z_Px_Fx2y_S_b+ABZ*I_ERI_Dxz_Px_Fx2y_S_b;
  Double I_ERI_D2y_Dxz_Fx2y_S_b = I_ERI_F2yz_Px_Fx2y_S_b+ABZ*I_ERI_D2y_Px_Fx2y_S_b;
  Double I_ERI_Dyz_Dxz_Fx2y_S_b = I_ERI_Fy2z_Px_Fx2y_S_b+ABZ*I_ERI_Dyz_Px_Fx2y_S_b;
  Double I_ERI_D2z_Dxz_Fx2y_S_b = I_ERI_F3z_Px_Fx2y_S_b+ABZ*I_ERI_D2z_Px_Fx2y_S_b;
  Double I_ERI_D2x_D2y_Fx2y_S_b = I_ERI_F2xy_Py_Fx2y_S_b+ABY*I_ERI_D2x_Py_Fx2y_S_b;
  Double I_ERI_Dxy_D2y_Fx2y_S_b = I_ERI_Fx2y_Py_Fx2y_S_b+ABY*I_ERI_Dxy_Py_Fx2y_S_b;
  Double I_ERI_Dxz_D2y_Fx2y_S_b = I_ERI_Fxyz_Py_Fx2y_S_b+ABY*I_ERI_Dxz_Py_Fx2y_S_b;
  Double I_ERI_D2y_D2y_Fx2y_S_b = I_ERI_F3y_Py_Fx2y_S_b+ABY*I_ERI_D2y_Py_Fx2y_S_b;
  Double I_ERI_Dyz_D2y_Fx2y_S_b = I_ERI_F2yz_Py_Fx2y_S_b+ABY*I_ERI_Dyz_Py_Fx2y_S_b;
  Double I_ERI_D2z_D2y_Fx2y_S_b = I_ERI_Fy2z_Py_Fx2y_S_b+ABY*I_ERI_D2z_Py_Fx2y_S_b;
  Double I_ERI_D2x_Dyz_Fx2y_S_b = I_ERI_F2xz_Py_Fx2y_S_b+ABZ*I_ERI_D2x_Py_Fx2y_S_b;
  Double I_ERI_Dxy_Dyz_Fx2y_S_b = I_ERI_Fxyz_Py_Fx2y_S_b+ABZ*I_ERI_Dxy_Py_Fx2y_S_b;
  Double I_ERI_Dxz_Dyz_Fx2y_S_b = I_ERI_Fx2z_Py_Fx2y_S_b+ABZ*I_ERI_Dxz_Py_Fx2y_S_b;
  Double I_ERI_D2y_Dyz_Fx2y_S_b = I_ERI_F2yz_Py_Fx2y_S_b+ABZ*I_ERI_D2y_Py_Fx2y_S_b;
  Double I_ERI_Dyz_Dyz_Fx2y_S_b = I_ERI_Fy2z_Py_Fx2y_S_b+ABZ*I_ERI_Dyz_Py_Fx2y_S_b;
  Double I_ERI_D2z_Dyz_Fx2y_S_b = I_ERI_F3z_Py_Fx2y_S_b+ABZ*I_ERI_D2z_Py_Fx2y_S_b;
  Double I_ERI_D2x_D2z_Fx2y_S_b = I_ERI_F2xz_Pz_Fx2y_S_b+ABZ*I_ERI_D2x_Pz_Fx2y_S_b;
  Double I_ERI_Dxy_D2z_Fx2y_S_b = I_ERI_Fxyz_Pz_Fx2y_S_b+ABZ*I_ERI_Dxy_Pz_Fx2y_S_b;
  Double I_ERI_Dxz_D2z_Fx2y_S_b = I_ERI_Fx2z_Pz_Fx2y_S_b+ABZ*I_ERI_Dxz_Pz_Fx2y_S_b;
  Double I_ERI_D2y_D2z_Fx2y_S_b = I_ERI_F2yz_Pz_Fx2y_S_b+ABZ*I_ERI_D2y_Pz_Fx2y_S_b;
  Double I_ERI_Dyz_D2z_Fx2y_S_b = I_ERI_Fy2z_Pz_Fx2y_S_b+ABZ*I_ERI_Dyz_Pz_Fx2y_S_b;
  Double I_ERI_D2z_D2z_Fx2y_S_b = I_ERI_F3z_Pz_Fx2y_S_b+ABZ*I_ERI_D2z_Pz_Fx2y_S_b;
  Double I_ERI_D2x_D2x_Fxyz_S_b = I_ERI_F3x_Px_Fxyz_S_b+ABX*I_ERI_D2x_Px_Fxyz_S_b;
  Double I_ERI_Dxy_D2x_Fxyz_S_b = I_ERI_F2xy_Px_Fxyz_S_b+ABX*I_ERI_Dxy_Px_Fxyz_S_b;
  Double I_ERI_Dxz_D2x_Fxyz_S_b = I_ERI_F2xz_Px_Fxyz_S_b+ABX*I_ERI_Dxz_Px_Fxyz_S_b;
  Double I_ERI_D2y_D2x_Fxyz_S_b = I_ERI_Fx2y_Px_Fxyz_S_b+ABX*I_ERI_D2y_Px_Fxyz_S_b;
  Double I_ERI_Dyz_D2x_Fxyz_S_b = I_ERI_Fxyz_Px_Fxyz_S_b+ABX*I_ERI_Dyz_Px_Fxyz_S_b;
  Double I_ERI_D2z_D2x_Fxyz_S_b = I_ERI_Fx2z_Px_Fxyz_S_b+ABX*I_ERI_D2z_Px_Fxyz_S_b;
  Double I_ERI_D2x_Dxy_Fxyz_S_b = I_ERI_F2xy_Px_Fxyz_S_b+ABY*I_ERI_D2x_Px_Fxyz_S_b;
  Double I_ERI_Dxy_Dxy_Fxyz_S_b = I_ERI_Fx2y_Px_Fxyz_S_b+ABY*I_ERI_Dxy_Px_Fxyz_S_b;
  Double I_ERI_Dxz_Dxy_Fxyz_S_b = I_ERI_Fxyz_Px_Fxyz_S_b+ABY*I_ERI_Dxz_Px_Fxyz_S_b;
  Double I_ERI_D2y_Dxy_Fxyz_S_b = I_ERI_F3y_Px_Fxyz_S_b+ABY*I_ERI_D2y_Px_Fxyz_S_b;
  Double I_ERI_Dyz_Dxy_Fxyz_S_b = I_ERI_F2yz_Px_Fxyz_S_b+ABY*I_ERI_Dyz_Px_Fxyz_S_b;
  Double I_ERI_D2z_Dxy_Fxyz_S_b = I_ERI_Fy2z_Px_Fxyz_S_b+ABY*I_ERI_D2z_Px_Fxyz_S_b;
  Double I_ERI_D2x_Dxz_Fxyz_S_b = I_ERI_F2xz_Px_Fxyz_S_b+ABZ*I_ERI_D2x_Px_Fxyz_S_b;
  Double I_ERI_Dxy_Dxz_Fxyz_S_b = I_ERI_Fxyz_Px_Fxyz_S_b+ABZ*I_ERI_Dxy_Px_Fxyz_S_b;
  Double I_ERI_Dxz_Dxz_Fxyz_S_b = I_ERI_Fx2z_Px_Fxyz_S_b+ABZ*I_ERI_Dxz_Px_Fxyz_S_b;
  Double I_ERI_D2y_Dxz_Fxyz_S_b = I_ERI_F2yz_Px_Fxyz_S_b+ABZ*I_ERI_D2y_Px_Fxyz_S_b;
  Double I_ERI_Dyz_Dxz_Fxyz_S_b = I_ERI_Fy2z_Px_Fxyz_S_b+ABZ*I_ERI_Dyz_Px_Fxyz_S_b;
  Double I_ERI_D2z_Dxz_Fxyz_S_b = I_ERI_F3z_Px_Fxyz_S_b+ABZ*I_ERI_D2z_Px_Fxyz_S_b;
  Double I_ERI_D2x_D2y_Fxyz_S_b = I_ERI_F2xy_Py_Fxyz_S_b+ABY*I_ERI_D2x_Py_Fxyz_S_b;
  Double I_ERI_Dxy_D2y_Fxyz_S_b = I_ERI_Fx2y_Py_Fxyz_S_b+ABY*I_ERI_Dxy_Py_Fxyz_S_b;
  Double I_ERI_Dxz_D2y_Fxyz_S_b = I_ERI_Fxyz_Py_Fxyz_S_b+ABY*I_ERI_Dxz_Py_Fxyz_S_b;
  Double I_ERI_D2y_D2y_Fxyz_S_b = I_ERI_F3y_Py_Fxyz_S_b+ABY*I_ERI_D2y_Py_Fxyz_S_b;
  Double I_ERI_Dyz_D2y_Fxyz_S_b = I_ERI_F2yz_Py_Fxyz_S_b+ABY*I_ERI_Dyz_Py_Fxyz_S_b;
  Double I_ERI_D2z_D2y_Fxyz_S_b = I_ERI_Fy2z_Py_Fxyz_S_b+ABY*I_ERI_D2z_Py_Fxyz_S_b;
  Double I_ERI_D2x_Dyz_Fxyz_S_b = I_ERI_F2xz_Py_Fxyz_S_b+ABZ*I_ERI_D2x_Py_Fxyz_S_b;
  Double I_ERI_Dxy_Dyz_Fxyz_S_b = I_ERI_Fxyz_Py_Fxyz_S_b+ABZ*I_ERI_Dxy_Py_Fxyz_S_b;
  Double I_ERI_Dxz_Dyz_Fxyz_S_b = I_ERI_Fx2z_Py_Fxyz_S_b+ABZ*I_ERI_Dxz_Py_Fxyz_S_b;
  Double I_ERI_D2y_Dyz_Fxyz_S_b = I_ERI_F2yz_Py_Fxyz_S_b+ABZ*I_ERI_D2y_Py_Fxyz_S_b;
  Double I_ERI_Dyz_Dyz_Fxyz_S_b = I_ERI_Fy2z_Py_Fxyz_S_b+ABZ*I_ERI_Dyz_Py_Fxyz_S_b;
  Double I_ERI_D2z_Dyz_Fxyz_S_b = I_ERI_F3z_Py_Fxyz_S_b+ABZ*I_ERI_D2z_Py_Fxyz_S_b;
  Double I_ERI_D2x_D2z_Fxyz_S_b = I_ERI_F2xz_Pz_Fxyz_S_b+ABZ*I_ERI_D2x_Pz_Fxyz_S_b;
  Double I_ERI_Dxy_D2z_Fxyz_S_b = I_ERI_Fxyz_Pz_Fxyz_S_b+ABZ*I_ERI_Dxy_Pz_Fxyz_S_b;
  Double I_ERI_Dxz_D2z_Fxyz_S_b = I_ERI_Fx2z_Pz_Fxyz_S_b+ABZ*I_ERI_Dxz_Pz_Fxyz_S_b;
  Double I_ERI_D2y_D2z_Fxyz_S_b = I_ERI_F2yz_Pz_Fxyz_S_b+ABZ*I_ERI_D2y_Pz_Fxyz_S_b;
  Double I_ERI_Dyz_D2z_Fxyz_S_b = I_ERI_Fy2z_Pz_Fxyz_S_b+ABZ*I_ERI_Dyz_Pz_Fxyz_S_b;
  Double I_ERI_D2z_D2z_Fxyz_S_b = I_ERI_F3z_Pz_Fxyz_S_b+ABZ*I_ERI_D2z_Pz_Fxyz_S_b;
  Double I_ERI_D2x_D2x_Fx2z_S_b = I_ERI_F3x_Px_Fx2z_S_b+ABX*I_ERI_D2x_Px_Fx2z_S_b;
  Double I_ERI_Dxy_D2x_Fx2z_S_b = I_ERI_F2xy_Px_Fx2z_S_b+ABX*I_ERI_Dxy_Px_Fx2z_S_b;
  Double I_ERI_Dxz_D2x_Fx2z_S_b = I_ERI_F2xz_Px_Fx2z_S_b+ABX*I_ERI_Dxz_Px_Fx2z_S_b;
  Double I_ERI_D2y_D2x_Fx2z_S_b = I_ERI_Fx2y_Px_Fx2z_S_b+ABX*I_ERI_D2y_Px_Fx2z_S_b;
  Double I_ERI_Dyz_D2x_Fx2z_S_b = I_ERI_Fxyz_Px_Fx2z_S_b+ABX*I_ERI_Dyz_Px_Fx2z_S_b;
  Double I_ERI_D2z_D2x_Fx2z_S_b = I_ERI_Fx2z_Px_Fx2z_S_b+ABX*I_ERI_D2z_Px_Fx2z_S_b;
  Double I_ERI_D2x_Dxy_Fx2z_S_b = I_ERI_F2xy_Px_Fx2z_S_b+ABY*I_ERI_D2x_Px_Fx2z_S_b;
  Double I_ERI_Dxy_Dxy_Fx2z_S_b = I_ERI_Fx2y_Px_Fx2z_S_b+ABY*I_ERI_Dxy_Px_Fx2z_S_b;
  Double I_ERI_Dxz_Dxy_Fx2z_S_b = I_ERI_Fxyz_Px_Fx2z_S_b+ABY*I_ERI_Dxz_Px_Fx2z_S_b;
  Double I_ERI_D2y_Dxy_Fx2z_S_b = I_ERI_F3y_Px_Fx2z_S_b+ABY*I_ERI_D2y_Px_Fx2z_S_b;
  Double I_ERI_Dyz_Dxy_Fx2z_S_b = I_ERI_F2yz_Px_Fx2z_S_b+ABY*I_ERI_Dyz_Px_Fx2z_S_b;
  Double I_ERI_D2z_Dxy_Fx2z_S_b = I_ERI_Fy2z_Px_Fx2z_S_b+ABY*I_ERI_D2z_Px_Fx2z_S_b;
  Double I_ERI_D2x_Dxz_Fx2z_S_b = I_ERI_F2xz_Px_Fx2z_S_b+ABZ*I_ERI_D2x_Px_Fx2z_S_b;
  Double I_ERI_Dxy_Dxz_Fx2z_S_b = I_ERI_Fxyz_Px_Fx2z_S_b+ABZ*I_ERI_Dxy_Px_Fx2z_S_b;
  Double I_ERI_Dxz_Dxz_Fx2z_S_b = I_ERI_Fx2z_Px_Fx2z_S_b+ABZ*I_ERI_Dxz_Px_Fx2z_S_b;
  Double I_ERI_D2y_Dxz_Fx2z_S_b = I_ERI_F2yz_Px_Fx2z_S_b+ABZ*I_ERI_D2y_Px_Fx2z_S_b;
  Double I_ERI_Dyz_Dxz_Fx2z_S_b = I_ERI_Fy2z_Px_Fx2z_S_b+ABZ*I_ERI_Dyz_Px_Fx2z_S_b;
  Double I_ERI_D2z_Dxz_Fx2z_S_b = I_ERI_F3z_Px_Fx2z_S_b+ABZ*I_ERI_D2z_Px_Fx2z_S_b;
  Double I_ERI_D2x_D2y_Fx2z_S_b = I_ERI_F2xy_Py_Fx2z_S_b+ABY*I_ERI_D2x_Py_Fx2z_S_b;
  Double I_ERI_Dxy_D2y_Fx2z_S_b = I_ERI_Fx2y_Py_Fx2z_S_b+ABY*I_ERI_Dxy_Py_Fx2z_S_b;
  Double I_ERI_Dxz_D2y_Fx2z_S_b = I_ERI_Fxyz_Py_Fx2z_S_b+ABY*I_ERI_Dxz_Py_Fx2z_S_b;
  Double I_ERI_D2y_D2y_Fx2z_S_b = I_ERI_F3y_Py_Fx2z_S_b+ABY*I_ERI_D2y_Py_Fx2z_S_b;
  Double I_ERI_Dyz_D2y_Fx2z_S_b = I_ERI_F2yz_Py_Fx2z_S_b+ABY*I_ERI_Dyz_Py_Fx2z_S_b;
  Double I_ERI_D2z_D2y_Fx2z_S_b = I_ERI_Fy2z_Py_Fx2z_S_b+ABY*I_ERI_D2z_Py_Fx2z_S_b;
  Double I_ERI_D2x_Dyz_Fx2z_S_b = I_ERI_F2xz_Py_Fx2z_S_b+ABZ*I_ERI_D2x_Py_Fx2z_S_b;
  Double I_ERI_Dxy_Dyz_Fx2z_S_b = I_ERI_Fxyz_Py_Fx2z_S_b+ABZ*I_ERI_Dxy_Py_Fx2z_S_b;
  Double I_ERI_Dxz_Dyz_Fx2z_S_b = I_ERI_Fx2z_Py_Fx2z_S_b+ABZ*I_ERI_Dxz_Py_Fx2z_S_b;
  Double I_ERI_D2y_Dyz_Fx2z_S_b = I_ERI_F2yz_Py_Fx2z_S_b+ABZ*I_ERI_D2y_Py_Fx2z_S_b;
  Double I_ERI_Dyz_Dyz_Fx2z_S_b = I_ERI_Fy2z_Py_Fx2z_S_b+ABZ*I_ERI_Dyz_Py_Fx2z_S_b;
  Double I_ERI_D2z_Dyz_Fx2z_S_b = I_ERI_F3z_Py_Fx2z_S_b+ABZ*I_ERI_D2z_Py_Fx2z_S_b;
  Double I_ERI_D2x_D2z_Fx2z_S_b = I_ERI_F2xz_Pz_Fx2z_S_b+ABZ*I_ERI_D2x_Pz_Fx2z_S_b;
  Double I_ERI_Dxy_D2z_Fx2z_S_b = I_ERI_Fxyz_Pz_Fx2z_S_b+ABZ*I_ERI_Dxy_Pz_Fx2z_S_b;
  Double I_ERI_Dxz_D2z_Fx2z_S_b = I_ERI_Fx2z_Pz_Fx2z_S_b+ABZ*I_ERI_Dxz_Pz_Fx2z_S_b;
  Double I_ERI_D2y_D2z_Fx2z_S_b = I_ERI_F2yz_Pz_Fx2z_S_b+ABZ*I_ERI_D2y_Pz_Fx2z_S_b;
  Double I_ERI_Dyz_D2z_Fx2z_S_b = I_ERI_Fy2z_Pz_Fx2z_S_b+ABZ*I_ERI_Dyz_Pz_Fx2z_S_b;
  Double I_ERI_D2z_D2z_Fx2z_S_b = I_ERI_F3z_Pz_Fx2z_S_b+ABZ*I_ERI_D2z_Pz_Fx2z_S_b;
  Double I_ERI_D2x_D2x_F3y_S_b = I_ERI_F3x_Px_F3y_S_b+ABX*I_ERI_D2x_Px_F3y_S_b;
  Double I_ERI_Dxy_D2x_F3y_S_b = I_ERI_F2xy_Px_F3y_S_b+ABX*I_ERI_Dxy_Px_F3y_S_b;
  Double I_ERI_Dxz_D2x_F3y_S_b = I_ERI_F2xz_Px_F3y_S_b+ABX*I_ERI_Dxz_Px_F3y_S_b;
  Double I_ERI_D2y_D2x_F3y_S_b = I_ERI_Fx2y_Px_F3y_S_b+ABX*I_ERI_D2y_Px_F3y_S_b;
  Double I_ERI_Dyz_D2x_F3y_S_b = I_ERI_Fxyz_Px_F3y_S_b+ABX*I_ERI_Dyz_Px_F3y_S_b;
  Double I_ERI_D2z_D2x_F3y_S_b = I_ERI_Fx2z_Px_F3y_S_b+ABX*I_ERI_D2z_Px_F3y_S_b;
  Double I_ERI_D2x_Dxy_F3y_S_b = I_ERI_F2xy_Px_F3y_S_b+ABY*I_ERI_D2x_Px_F3y_S_b;
  Double I_ERI_Dxy_Dxy_F3y_S_b = I_ERI_Fx2y_Px_F3y_S_b+ABY*I_ERI_Dxy_Px_F3y_S_b;
  Double I_ERI_Dxz_Dxy_F3y_S_b = I_ERI_Fxyz_Px_F3y_S_b+ABY*I_ERI_Dxz_Px_F3y_S_b;
  Double I_ERI_D2y_Dxy_F3y_S_b = I_ERI_F3y_Px_F3y_S_b+ABY*I_ERI_D2y_Px_F3y_S_b;
  Double I_ERI_Dyz_Dxy_F3y_S_b = I_ERI_F2yz_Px_F3y_S_b+ABY*I_ERI_Dyz_Px_F3y_S_b;
  Double I_ERI_D2z_Dxy_F3y_S_b = I_ERI_Fy2z_Px_F3y_S_b+ABY*I_ERI_D2z_Px_F3y_S_b;
  Double I_ERI_D2x_Dxz_F3y_S_b = I_ERI_F2xz_Px_F3y_S_b+ABZ*I_ERI_D2x_Px_F3y_S_b;
  Double I_ERI_Dxy_Dxz_F3y_S_b = I_ERI_Fxyz_Px_F3y_S_b+ABZ*I_ERI_Dxy_Px_F3y_S_b;
  Double I_ERI_Dxz_Dxz_F3y_S_b = I_ERI_Fx2z_Px_F3y_S_b+ABZ*I_ERI_Dxz_Px_F3y_S_b;
  Double I_ERI_D2y_Dxz_F3y_S_b = I_ERI_F2yz_Px_F3y_S_b+ABZ*I_ERI_D2y_Px_F3y_S_b;
  Double I_ERI_Dyz_Dxz_F3y_S_b = I_ERI_Fy2z_Px_F3y_S_b+ABZ*I_ERI_Dyz_Px_F3y_S_b;
  Double I_ERI_D2z_Dxz_F3y_S_b = I_ERI_F3z_Px_F3y_S_b+ABZ*I_ERI_D2z_Px_F3y_S_b;
  Double I_ERI_D2x_D2y_F3y_S_b = I_ERI_F2xy_Py_F3y_S_b+ABY*I_ERI_D2x_Py_F3y_S_b;
  Double I_ERI_Dxy_D2y_F3y_S_b = I_ERI_Fx2y_Py_F3y_S_b+ABY*I_ERI_Dxy_Py_F3y_S_b;
  Double I_ERI_Dxz_D2y_F3y_S_b = I_ERI_Fxyz_Py_F3y_S_b+ABY*I_ERI_Dxz_Py_F3y_S_b;
  Double I_ERI_D2y_D2y_F3y_S_b = I_ERI_F3y_Py_F3y_S_b+ABY*I_ERI_D2y_Py_F3y_S_b;
  Double I_ERI_Dyz_D2y_F3y_S_b = I_ERI_F2yz_Py_F3y_S_b+ABY*I_ERI_Dyz_Py_F3y_S_b;
  Double I_ERI_D2z_D2y_F3y_S_b = I_ERI_Fy2z_Py_F3y_S_b+ABY*I_ERI_D2z_Py_F3y_S_b;
  Double I_ERI_D2x_Dyz_F3y_S_b = I_ERI_F2xz_Py_F3y_S_b+ABZ*I_ERI_D2x_Py_F3y_S_b;
  Double I_ERI_Dxy_Dyz_F3y_S_b = I_ERI_Fxyz_Py_F3y_S_b+ABZ*I_ERI_Dxy_Py_F3y_S_b;
  Double I_ERI_Dxz_Dyz_F3y_S_b = I_ERI_Fx2z_Py_F3y_S_b+ABZ*I_ERI_Dxz_Py_F3y_S_b;
  Double I_ERI_D2y_Dyz_F3y_S_b = I_ERI_F2yz_Py_F3y_S_b+ABZ*I_ERI_D2y_Py_F3y_S_b;
  Double I_ERI_Dyz_Dyz_F3y_S_b = I_ERI_Fy2z_Py_F3y_S_b+ABZ*I_ERI_Dyz_Py_F3y_S_b;
  Double I_ERI_D2z_Dyz_F3y_S_b = I_ERI_F3z_Py_F3y_S_b+ABZ*I_ERI_D2z_Py_F3y_S_b;
  Double I_ERI_D2x_D2z_F3y_S_b = I_ERI_F2xz_Pz_F3y_S_b+ABZ*I_ERI_D2x_Pz_F3y_S_b;
  Double I_ERI_Dxy_D2z_F3y_S_b = I_ERI_Fxyz_Pz_F3y_S_b+ABZ*I_ERI_Dxy_Pz_F3y_S_b;
  Double I_ERI_Dxz_D2z_F3y_S_b = I_ERI_Fx2z_Pz_F3y_S_b+ABZ*I_ERI_Dxz_Pz_F3y_S_b;
  Double I_ERI_D2y_D2z_F3y_S_b = I_ERI_F2yz_Pz_F3y_S_b+ABZ*I_ERI_D2y_Pz_F3y_S_b;
  Double I_ERI_Dyz_D2z_F3y_S_b = I_ERI_Fy2z_Pz_F3y_S_b+ABZ*I_ERI_Dyz_Pz_F3y_S_b;
  Double I_ERI_D2z_D2z_F3y_S_b = I_ERI_F3z_Pz_F3y_S_b+ABZ*I_ERI_D2z_Pz_F3y_S_b;
  Double I_ERI_D2x_D2x_F2yz_S_b = I_ERI_F3x_Px_F2yz_S_b+ABX*I_ERI_D2x_Px_F2yz_S_b;
  Double I_ERI_Dxy_D2x_F2yz_S_b = I_ERI_F2xy_Px_F2yz_S_b+ABX*I_ERI_Dxy_Px_F2yz_S_b;
  Double I_ERI_Dxz_D2x_F2yz_S_b = I_ERI_F2xz_Px_F2yz_S_b+ABX*I_ERI_Dxz_Px_F2yz_S_b;
  Double I_ERI_D2y_D2x_F2yz_S_b = I_ERI_Fx2y_Px_F2yz_S_b+ABX*I_ERI_D2y_Px_F2yz_S_b;
  Double I_ERI_Dyz_D2x_F2yz_S_b = I_ERI_Fxyz_Px_F2yz_S_b+ABX*I_ERI_Dyz_Px_F2yz_S_b;
  Double I_ERI_D2z_D2x_F2yz_S_b = I_ERI_Fx2z_Px_F2yz_S_b+ABX*I_ERI_D2z_Px_F2yz_S_b;
  Double I_ERI_D2x_Dxy_F2yz_S_b = I_ERI_F2xy_Px_F2yz_S_b+ABY*I_ERI_D2x_Px_F2yz_S_b;
  Double I_ERI_Dxy_Dxy_F2yz_S_b = I_ERI_Fx2y_Px_F2yz_S_b+ABY*I_ERI_Dxy_Px_F2yz_S_b;
  Double I_ERI_Dxz_Dxy_F2yz_S_b = I_ERI_Fxyz_Px_F2yz_S_b+ABY*I_ERI_Dxz_Px_F2yz_S_b;
  Double I_ERI_D2y_Dxy_F2yz_S_b = I_ERI_F3y_Px_F2yz_S_b+ABY*I_ERI_D2y_Px_F2yz_S_b;
  Double I_ERI_Dyz_Dxy_F2yz_S_b = I_ERI_F2yz_Px_F2yz_S_b+ABY*I_ERI_Dyz_Px_F2yz_S_b;
  Double I_ERI_D2z_Dxy_F2yz_S_b = I_ERI_Fy2z_Px_F2yz_S_b+ABY*I_ERI_D2z_Px_F2yz_S_b;
  Double I_ERI_D2x_Dxz_F2yz_S_b = I_ERI_F2xz_Px_F2yz_S_b+ABZ*I_ERI_D2x_Px_F2yz_S_b;
  Double I_ERI_Dxy_Dxz_F2yz_S_b = I_ERI_Fxyz_Px_F2yz_S_b+ABZ*I_ERI_Dxy_Px_F2yz_S_b;
  Double I_ERI_Dxz_Dxz_F2yz_S_b = I_ERI_Fx2z_Px_F2yz_S_b+ABZ*I_ERI_Dxz_Px_F2yz_S_b;
  Double I_ERI_D2y_Dxz_F2yz_S_b = I_ERI_F2yz_Px_F2yz_S_b+ABZ*I_ERI_D2y_Px_F2yz_S_b;
  Double I_ERI_Dyz_Dxz_F2yz_S_b = I_ERI_Fy2z_Px_F2yz_S_b+ABZ*I_ERI_Dyz_Px_F2yz_S_b;
  Double I_ERI_D2z_Dxz_F2yz_S_b = I_ERI_F3z_Px_F2yz_S_b+ABZ*I_ERI_D2z_Px_F2yz_S_b;
  Double I_ERI_D2x_D2y_F2yz_S_b = I_ERI_F2xy_Py_F2yz_S_b+ABY*I_ERI_D2x_Py_F2yz_S_b;
  Double I_ERI_Dxy_D2y_F2yz_S_b = I_ERI_Fx2y_Py_F2yz_S_b+ABY*I_ERI_Dxy_Py_F2yz_S_b;
  Double I_ERI_Dxz_D2y_F2yz_S_b = I_ERI_Fxyz_Py_F2yz_S_b+ABY*I_ERI_Dxz_Py_F2yz_S_b;
  Double I_ERI_D2y_D2y_F2yz_S_b = I_ERI_F3y_Py_F2yz_S_b+ABY*I_ERI_D2y_Py_F2yz_S_b;
  Double I_ERI_Dyz_D2y_F2yz_S_b = I_ERI_F2yz_Py_F2yz_S_b+ABY*I_ERI_Dyz_Py_F2yz_S_b;
  Double I_ERI_D2z_D2y_F2yz_S_b = I_ERI_Fy2z_Py_F2yz_S_b+ABY*I_ERI_D2z_Py_F2yz_S_b;
  Double I_ERI_D2x_Dyz_F2yz_S_b = I_ERI_F2xz_Py_F2yz_S_b+ABZ*I_ERI_D2x_Py_F2yz_S_b;
  Double I_ERI_Dxy_Dyz_F2yz_S_b = I_ERI_Fxyz_Py_F2yz_S_b+ABZ*I_ERI_Dxy_Py_F2yz_S_b;
  Double I_ERI_Dxz_Dyz_F2yz_S_b = I_ERI_Fx2z_Py_F2yz_S_b+ABZ*I_ERI_Dxz_Py_F2yz_S_b;
  Double I_ERI_D2y_Dyz_F2yz_S_b = I_ERI_F2yz_Py_F2yz_S_b+ABZ*I_ERI_D2y_Py_F2yz_S_b;
  Double I_ERI_Dyz_Dyz_F2yz_S_b = I_ERI_Fy2z_Py_F2yz_S_b+ABZ*I_ERI_Dyz_Py_F2yz_S_b;
  Double I_ERI_D2z_Dyz_F2yz_S_b = I_ERI_F3z_Py_F2yz_S_b+ABZ*I_ERI_D2z_Py_F2yz_S_b;
  Double I_ERI_D2x_D2z_F2yz_S_b = I_ERI_F2xz_Pz_F2yz_S_b+ABZ*I_ERI_D2x_Pz_F2yz_S_b;
  Double I_ERI_Dxy_D2z_F2yz_S_b = I_ERI_Fxyz_Pz_F2yz_S_b+ABZ*I_ERI_Dxy_Pz_F2yz_S_b;
  Double I_ERI_Dxz_D2z_F2yz_S_b = I_ERI_Fx2z_Pz_F2yz_S_b+ABZ*I_ERI_Dxz_Pz_F2yz_S_b;
  Double I_ERI_D2y_D2z_F2yz_S_b = I_ERI_F2yz_Pz_F2yz_S_b+ABZ*I_ERI_D2y_Pz_F2yz_S_b;
  Double I_ERI_Dyz_D2z_F2yz_S_b = I_ERI_Fy2z_Pz_F2yz_S_b+ABZ*I_ERI_Dyz_Pz_F2yz_S_b;
  Double I_ERI_D2z_D2z_F2yz_S_b = I_ERI_F3z_Pz_F2yz_S_b+ABZ*I_ERI_D2z_Pz_F2yz_S_b;
  Double I_ERI_D2x_D2x_Fy2z_S_b = I_ERI_F3x_Px_Fy2z_S_b+ABX*I_ERI_D2x_Px_Fy2z_S_b;
  Double I_ERI_Dxy_D2x_Fy2z_S_b = I_ERI_F2xy_Px_Fy2z_S_b+ABX*I_ERI_Dxy_Px_Fy2z_S_b;
  Double I_ERI_Dxz_D2x_Fy2z_S_b = I_ERI_F2xz_Px_Fy2z_S_b+ABX*I_ERI_Dxz_Px_Fy2z_S_b;
  Double I_ERI_D2y_D2x_Fy2z_S_b = I_ERI_Fx2y_Px_Fy2z_S_b+ABX*I_ERI_D2y_Px_Fy2z_S_b;
  Double I_ERI_Dyz_D2x_Fy2z_S_b = I_ERI_Fxyz_Px_Fy2z_S_b+ABX*I_ERI_Dyz_Px_Fy2z_S_b;
  Double I_ERI_D2z_D2x_Fy2z_S_b = I_ERI_Fx2z_Px_Fy2z_S_b+ABX*I_ERI_D2z_Px_Fy2z_S_b;
  Double I_ERI_D2x_Dxy_Fy2z_S_b = I_ERI_F2xy_Px_Fy2z_S_b+ABY*I_ERI_D2x_Px_Fy2z_S_b;
  Double I_ERI_Dxy_Dxy_Fy2z_S_b = I_ERI_Fx2y_Px_Fy2z_S_b+ABY*I_ERI_Dxy_Px_Fy2z_S_b;
  Double I_ERI_Dxz_Dxy_Fy2z_S_b = I_ERI_Fxyz_Px_Fy2z_S_b+ABY*I_ERI_Dxz_Px_Fy2z_S_b;
  Double I_ERI_D2y_Dxy_Fy2z_S_b = I_ERI_F3y_Px_Fy2z_S_b+ABY*I_ERI_D2y_Px_Fy2z_S_b;
  Double I_ERI_Dyz_Dxy_Fy2z_S_b = I_ERI_F2yz_Px_Fy2z_S_b+ABY*I_ERI_Dyz_Px_Fy2z_S_b;
  Double I_ERI_D2z_Dxy_Fy2z_S_b = I_ERI_Fy2z_Px_Fy2z_S_b+ABY*I_ERI_D2z_Px_Fy2z_S_b;
  Double I_ERI_D2x_Dxz_Fy2z_S_b = I_ERI_F2xz_Px_Fy2z_S_b+ABZ*I_ERI_D2x_Px_Fy2z_S_b;
  Double I_ERI_Dxy_Dxz_Fy2z_S_b = I_ERI_Fxyz_Px_Fy2z_S_b+ABZ*I_ERI_Dxy_Px_Fy2z_S_b;
  Double I_ERI_Dxz_Dxz_Fy2z_S_b = I_ERI_Fx2z_Px_Fy2z_S_b+ABZ*I_ERI_Dxz_Px_Fy2z_S_b;
  Double I_ERI_D2y_Dxz_Fy2z_S_b = I_ERI_F2yz_Px_Fy2z_S_b+ABZ*I_ERI_D2y_Px_Fy2z_S_b;
  Double I_ERI_Dyz_Dxz_Fy2z_S_b = I_ERI_Fy2z_Px_Fy2z_S_b+ABZ*I_ERI_Dyz_Px_Fy2z_S_b;
  Double I_ERI_D2z_Dxz_Fy2z_S_b = I_ERI_F3z_Px_Fy2z_S_b+ABZ*I_ERI_D2z_Px_Fy2z_S_b;
  Double I_ERI_D2x_D2y_Fy2z_S_b = I_ERI_F2xy_Py_Fy2z_S_b+ABY*I_ERI_D2x_Py_Fy2z_S_b;
  Double I_ERI_Dxy_D2y_Fy2z_S_b = I_ERI_Fx2y_Py_Fy2z_S_b+ABY*I_ERI_Dxy_Py_Fy2z_S_b;
  Double I_ERI_Dxz_D2y_Fy2z_S_b = I_ERI_Fxyz_Py_Fy2z_S_b+ABY*I_ERI_Dxz_Py_Fy2z_S_b;
  Double I_ERI_D2y_D2y_Fy2z_S_b = I_ERI_F3y_Py_Fy2z_S_b+ABY*I_ERI_D2y_Py_Fy2z_S_b;
  Double I_ERI_Dyz_D2y_Fy2z_S_b = I_ERI_F2yz_Py_Fy2z_S_b+ABY*I_ERI_Dyz_Py_Fy2z_S_b;
  Double I_ERI_D2z_D2y_Fy2z_S_b = I_ERI_Fy2z_Py_Fy2z_S_b+ABY*I_ERI_D2z_Py_Fy2z_S_b;
  Double I_ERI_D2x_Dyz_Fy2z_S_b = I_ERI_F2xz_Py_Fy2z_S_b+ABZ*I_ERI_D2x_Py_Fy2z_S_b;
  Double I_ERI_Dxy_Dyz_Fy2z_S_b = I_ERI_Fxyz_Py_Fy2z_S_b+ABZ*I_ERI_Dxy_Py_Fy2z_S_b;
  Double I_ERI_Dxz_Dyz_Fy2z_S_b = I_ERI_Fx2z_Py_Fy2z_S_b+ABZ*I_ERI_Dxz_Py_Fy2z_S_b;
  Double I_ERI_D2y_Dyz_Fy2z_S_b = I_ERI_F2yz_Py_Fy2z_S_b+ABZ*I_ERI_D2y_Py_Fy2z_S_b;
  Double I_ERI_Dyz_Dyz_Fy2z_S_b = I_ERI_Fy2z_Py_Fy2z_S_b+ABZ*I_ERI_Dyz_Py_Fy2z_S_b;
  Double I_ERI_D2z_Dyz_Fy2z_S_b = I_ERI_F3z_Py_Fy2z_S_b+ABZ*I_ERI_D2z_Py_Fy2z_S_b;
  Double I_ERI_D2x_D2z_Fy2z_S_b = I_ERI_F2xz_Pz_Fy2z_S_b+ABZ*I_ERI_D2x_Pz_Fy2z_S_b;
  Double I_ERI_Dxy_D2z_Fy2z_S_b = I_ERI_Fxyz_Pz_Fy2z_S_b+ABZ*I_ERI_Dxy_Pz_Fy2z_S_b;
  Double I_ERI_Dxz_D2z_Fy2z_S_b = I_ERI_Fx2z_Pz_Fy2z_S_b+ABZ*I_ERI_Dxz_Pz_Fy2z_S_b;
  Double I_ERI_D2y_D2z_Fy2z_S_b = I_ERI_F2yz_Pz_Fy2z_S_b+ABZ*I_ERI_D2y_Pz_Fy2z_S_b;
  Double I_ERI_Dyz_D2z_Fy2z_S_b = I_ERI_Fy2z_Pz_Fy2z_S_b+ABZ*I_ERI_Dyz_Pz_Fy2z_S_b;
  Double I_ERI_D2z_D2z_Fy2z_S_b = I_ERI_F3z_Pz_Fy2z_S_b+ABZ*I_ERI_D2z_Pz_Fy2z_S_b;
  Double I_ERI_D2x_D2x_F3z_S_b = I_ERI_F3x_Px_F3z_S_b+ABX*I_ERI_D2x_Px_F3z_S_b;
  Double I_ERI_Dxy_D2x_F3z_S_b = I_ERI_F2xy_Px_F3z_S_b+ABX*I_ERI_Dxy_Px_F3z_S_b;
  Double I_ERI_Dxz_D2x_F3z_S_b = I_ERI_F2xz_Px_F3z_S_b+ABX*I_ERI_Dxz_Px_F3z_S_b;
  Double I_ERI_D2y_D2x_F3z_S_b = I_ERI_Fx2y_Px_F3z_S_b+ABX*I_ERI_D2y_Px_F3z_S_b;
  Double I_ERI_Dyz_D2x_F3z_S_b = I_ERI_Fxyz_Px_F3z_S_b+ABX*I_ERI_Dyz_Px_F3z_S_b;
  Double I_ERI_D2z_D2x_F3z_S_b = I_ERI_Fx2z_Px_F3z_S_b+ABX*I_ERI_D2z_Px_F3z_S_b;
  Double I_ERI_D2x_Dxy_F3z_S_b = I_ERI_F2xy_Px_F3z_S_b+ABY*I_ERI_D2x_Px_F3z_S_b;
  Double I_ERI_Dxy_Dxy_F3z_S_b = I_ERI_Fx2y_Px_F3z_S_b+ABY*I_ERI_Dxy_Px_F3z_S_b;
  Double I_ERI_Dxz_Dxy_F3z_S_b = I_ERI_Fxyz_Px_F3z_S_b+ABY*I_ERI_Dxz_Px_F3z_S_b;
  Double I_ERI_D2y_Dxy_F3z_S_b = I_ERI_F3y_Px_F3z_S_b+ABY*I_ERI_D2y_Px_F3z_S_b;
  Double I_ERI_Dyz_Dxy_F3z_S_b = I_ERI_F2yz_Px_F3z_S_b+ABY*I_ERI_Dyz_Px_F3z_S_b;
  Double I_ERI_D2z_Dxy_F3z_S_b = I_ERI_Fy2z_Px_F3z_S_b+ABY*I_ERI_D2z_Px_F3z_S_b;
  Double I_ERI_D2x_Dxz_F3z_S_b = I_ERI_F2xz_Px_F3z_S_b+ABZ*I_ERI_D2x_Px_F3z_S_b;
  Double I_ERI_Dxy_Dxz_F3z_S_b = I_ERI_Fxyz_Px_F3z_S_b+ABZ*I_ERI_Dxy_Px_F3z_S_b;
  Double I_ERI_Dxz_Dxz_F3z_S_b = I_ERI_Fx2z_Px_F3z_S_b+ABZ*I_ERI_Dxz_Px_F3z_S_b;
  Double I_ERI_D2y_Dxz_F3z_S_b = I_ERI_F2yz_Px_F3z_S_b+ABZ*I_ERI_D2y_Px_F3z_S_b;
  Double I_ERI_Dyz_Dxz_F3z_S_b = I_ERI_Fy2z_Px_F3z_S_b+ABZ*I_ERI_Dyz_Px_F3z_S_b;
  Double I_ERI_D2z_Dxz_F3z_S_b = I_ERI_F3z_Px_F3z_S_b+ABZ*I_ERI_D2z_Px_F3z_S_b;
  Double I_ERI_D2x_D2y_F3z_S_b = I_ERI_F2xy_Py_F3z_S_b+ABY*I_ERI_D2x_Py_F3z_S_b;
  Double I_ERI_Dxy_D2y_F3z_S_b = I_ERI_Fx2y_Py_F3z_S_b+ABY*I_ERI_Dxy_Py_F3z_S_b;
  Double I_ERI_Dxz_D2y_F3z_S_b = I_ERI_Fxyz_Py_F3z_S_b+ABY*I_ERI_Dxz_Py_F3z_S_b;
  Double I_ERI_D2y_D2y_F3z_S_b = I_ERI_F3y_Py_F3z_S_b+ABY*I_ERI_D2y_Py_F3z_S_b;
  Double I_ERI_Dyz_D2y_F3z_S_b = I_ERI_F2yz_Py_F3z_S_b+ABY*I_ERI_Dyz_Py_F3z_S_b;
  Double I_ERI_D2z_D2y_F3z_S_b = I_ERI_Fy2z_Py_F3z_S_b+ABY*I_ERI_D2z_Py_F3z_S_b;
  Double I_ERI_D2x_Dyz_F3z_S_b = I_ERI_F2xz_Py_F3z_S_b+ABZ*I_ERI_D2x_Py_F3z_S_b;
  Double I_ERI_Dxy_Dyz_F3z_S_b = I_ERI_Fxyz_Py_F3z_S_b+ABZ*I_ERI_Dxy_Py_F3z_S_b;
  Double I_ERI_Dxz_Dyz_F3z_S_b = I_ERI_Fx2z_Py_F3z_S_b+ABZ*I_ERI_Dxz_Py_F3z_S_b;
  Double I_ERI_D2y_Dyz_F3z_S_b = I_ERI_F2yz_Py_F3z_S_b+ABZ*I_ERI_D2y_Py_F3z_S_b;
  Double I_ERI_Dyz_Dyz_F3z_S_b = I_ERI_Fy2z_Py_F3z_S_b+ABZ*I_ERI_Dyz_Py_F3z_S_b;
  Double I_ERI_D2z_Dyz_F3z_S_b = I_ERI_F3z_Py_F3z_S_b+ABZ*I_ERI_D2z_Py_F3z_S_b;
  Double I_ERI_D2x_D2z_F3z_S_b = I_ERI_F2xz_Pz_F3z_S_b+ABZ*I_ERI_D2x_Pz_F3z_S_b;
  Double I_ERI_Dxy_D2z_F3z_S_b = I_ERI_Fxyz_Pz_F3z_S_b+ABZ*I_ERI_Dxy_Pz_F3z_S_b;
  Double I_ERI_Dxz_D2z_F3z_S_b = I_ERI_Fx2z_Pz_F3z_S_b+ABZ*I_ERI_Dxz_Pz_F3z_S_b;
  Double I_ERI_D2y_D2z_F3z_S_b = I_ERI_F2yz_Pz_F3z_S_b+ABZ*I_ERI_D2y_Pz_F3z_S_b;
  Double I_ERI_Dyz_D2z_F3z_S_b = I_ERI_Fy2z_Pz_F3z_S_b+ABZ*I_ERI_Dyz_Pz_F3z_S_b;
  Double I_ERI_D2z_D2z_F3z_S_b = I_ERI_F3z_Pz_F3z_S_b+ABZ*I_ERI_D2z_Pz_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_P_d
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P_d
   * RHS shell quartet name: SQ_ERI_D_S_F_P_d
   ************************************************************/
  Double I_ERI_D2x_Px_F3x_Px_d = I_ERI_F3x_S_F3x_Px_d+ABX*I_ERI_D2x_S_F3x_Px_d;
  Double I_ERI_Dxy_Px_F3x_Px_d = I_ERI_F2xy_S_F3x_Px_d+ABX*I_ERI_Dxy_S_F3x_Px_d;
  Double I_ERI_Dxz_Px_F3x_Px_d = I_ERI_F2xz_S_F3x_Px_d+ABX*I_ERI_Dxz_S_F3x_Px_d;
  Double I_ERI_D2y_Px_F3x_Px_d = I_ERI_Fx2y_S_F3x_Px_d+ABX*I_ERI_D2y_S_F3x_Px_d;
  Double I_ERI_Dyz_Px_F3x_Px_d = I_ERI_Fxyz_S_F3x_Px_d+ABX*I_ERI_Dyz_S_F3x_Px_d;
  Double I_ERI_D2z_Px_F3x_Px_d = I_ERI_Fx2z_S_F3x_Px_d+ABX*I_ERI_D2z_S_F3x_Px_d;
  Double I_ERI_D2x_Py_F3x_Px_d = I_ERI_F2xy_S_F3x_Px_d+ABY*I_ERI_D2x_S_F3x_Px_d;
  Double I_ERI_Dxy_Py_F3x_Px_d = I_ERI_Fx2y_S_F3x_Px_d+ABY*I_ERI_Dxy_S_F3x_Px_d;
  Double I_ERI_Dxz_Py_F3x_Px_d = I_ERI_Fxyz_S_F3x_Px_d+ABY*I_ERI_Dxz_S_F3x_Px_d;
  Double I_ERI_D2y_Py_F3x_Px_d = I_ERI_F3y_S_F3x_Px_d+ABY*I_ERI_D2y_S_F3x_Px_d;
  Double I_ERI_Dyz_Py_F3x_Px_d = I_ERI_F2yz_S_F3x_Px_d+ABY*I_ERI_Dyz_S_F3x_Px_d;
  Double I_ERI_D2z_Py_F3x_Px_d = I_ERI_Fy2z_S_F3x_Px_d+ABY*I_ERI_D2z_S_F3x_Px_d;
  Double I_ERI_D2x_Pz_F3x_Px_d = I_ERI_F2xz_S_F3x_Px_d+ABZ*I_ERI_D2x_S_F3x_Px_d;
  Double I_ERI_Dxy_Pz_F3x_Px_d = I_ERI_Fxyz_S_F3x_Px_d+ABZ*I_ERI_Dxy_S_F3x_Px_d;
  Double I_ERI_Dxz_Pz_F3x_Px_d = I_ERI_Fx2z_S_F3x_Px_d+ABZ*I_ERI_Dxz_S_F3x_Px_d;
  Double I_ERI_D2y_Pz_F3x_Px_d = I_ERI_F2yz_S_F3x_Px_d+ABZ*I_ERI_D2y_S_F3x_Px_d;
  Double I_ERI_Dyz_Pz_F3x_Px_d = I_ERI_Fy2z_S_F3x_Px_d+ABZ*I_ERI_Dyz_S_F3x_Px_d;
  Double I_ERI_D2z_Pz_F3x_Px_d = I_ERI_F3z_S_F3x_Px_d+ABZ*I_ERI_D2z_S_F3x_Px_d;
  Double I_ERI_D2x_Px_F2xy_Px_d = I_ERI_F3x_S_F2xy_Px_d+ABX*I_ERI_D2x_S_F2xy_Px_d;
  Double I_ERI_Dxy_Px_F2xy_Px_d = I_ERI_F2xy_S_F2xy_Px_d+ABX*I_ERI_Dxy_S_F2xy_Px_d;
  Double I_ERI_Dxz_Px_F2xy_Px_d = I_ERI_F2xz_S_F2xy_Px_d+ABX*I_ERI_Dxz_S_F2xy_Px_d;
  Double I_ERI_D2y_Px_F2xy_Px_d = I_ERI_Fx2y_S_F2xy_Px_d+ABX*I_ERI_D2y_S_F2xy_Px_d;
  Double I_ERI_Dyz_Px_F2xy_Px_d = I_ERI_Fxyz_S_F2xy_Px_d+ABX*I_ERI_Dyz_S_F2xy_Px_d;
  Double I_ERI_D2z_Px_F2xy_Px_d = I_ERI_Fx2z_S_F2xy_Px_d+ABX*I_ERI_D2z_S_F2xy_Px_d;
  Double I_ERI_D2x_Py_F2xy_Px_d = I_ERI_F2xy_S_F2xy_Px_d+ABY*I_ERI_D2x_S_F2xy_Px_d;
  Double I_ERI_Dxy_Py_F2xy_Px_d = I_ERI_Fx2y_S_F2xy_Px_d+ABY*I_ERI_Dxy_S_F2xy_Px_d;
  Double I_ERI_Dxz_Py_F2xy_Px_d = I_ERI_Fxyz_S_F2xy_Px_d+ABY*I_ERI_Dxz_S_F2xy_Px_d;
  Double I_ERI_D2y_Py_F2xy_Px_d = I_ERI_F3y_S_F2xy_Px_d+ABY*I_ERI_D2y_S_F2xy_Px_d;
  Double I_ERI_Dyz_Py_F2xy_Px_d = I_ERI_F2yz_S_F2xy_Px_d+ABY*I_ERI_Dyz_S_F2xy_Px_d;
  Double I_ERI_D2z_Py_F2xy_Px_d = I_ERI_Fy2z_S_F2xy_Px_d+ABY*I_ERI_D2z_S_F2xy_Px_d;
  Double I_ERI_D2x_Pz_F2xy_Px_d = I_ERI_F2xz_S_F2xy_Px_d+ABZ*I_ERI_D2x_S_F2xy_Px_d;
  Double I_ERI_Dxy_Pz_F2xy_Px_d = I_ERI_Fxyz_S_F2xy_Px_d+ABZ*I_ERI_Dxy_S_F2xy_Px_d;
  Double I_ERI_Dxz_Pz_F2xy_Px_d = I_ERI_Fx2z_S_F2xy_Px_d+ABZ*I_ERI_Dxz_S_F2xy_Px_d;
  Double I_ERI_D2y_Pz_F2xy_Px_d = I_ERI_F2yz_S_F2xy_Px_d+ABZ*I_ERI_D2y_S_F2xy_Px_d;
  Double I_ERI_Dyz_Pz_F2xy_Px_d = I_ERI_Fy2z_S_F2xy_Px_d+ABZ*I_ERI_Dyz_S_F2xy_Px_d;
  Double I_ERI_D2z_Pz_F2xy_Px_d = I_ERI_F3z_S_F2xy_Px_d+ABZ*I_ERI_D2z_S_F2xy_Px_d;
  Double I_ERI_D2x_Px_F2xz_Px_d = I_ERI_F3x_S_F2xz_Px_d+ABX*I_ERI_D2x_S_F2xz_Px_d;
  Double I_ERI_Dxy_Px_F2xz_Px_d = I_ERI_F2xy_S_F2xz_Px_d+ABX*I_ERI_Dxy_S_F2xz_Px_d;
  Double I_ERI_Dxz_Px_F2xz_Px_d = I_ERI_F2xz_S_F2xz_Px_d+ABX*I_ERI_Dxz_S_F2xz_Px_d;
  Double I_ERI_D2y_Px_F2xz_Px_d = I_ERI_Fx2y_S_F2xz_Px_d+ABX*I_ERI_D2y_S_F2xz_Px_d;
  Double I_ERI_Dyz_Px_F2xz_Px_d = I_ERI_Fxyz_S_F2xz_Px_d+ABX*I_ERI_Dyz_S_F2xz_Px_d;
  Double I_ERI_D2z_Px_F2xz_Px_d = I_ERI_Fx2z_S_F2xz_Px_d+ABX*I_ERI_D2z_S_F2xz_Px_d;
  Double I_ERI_D2x_Py_F2xz_Px_d = I_ERI_F2xy_S_F2xz_Px_d+ABY*I_ERI_D2x_S_F2xz_Px_d;
  Double I_ERI_Dxy_Py_F2xz_Px_d = I_ERI_Fx2y_S_F2xz_Px_d+ABY*I_ERI_Dxy_S_F2xz_Px_d;
  Double I_ERI_Dxz_Py_F2xz_Px_d = I_ERI_Fxyz_S_F2xz_Px_d+ABY*I_ERI_Dxz_S_F2xz_Px_d;
  Double I_ERI_D2y_Py_F2xz_Px_d = I_ERI_F3y_S_F2xz_Px_d+ABY*I_ERI_D2y_S_F2xz_Px_d;
  Double I_ERI_Dyz_Py_F2xz_Px_d = I_ERI_F2yz_S_F2xz_Px_d+ABY*I_ERI_Dyz_S_F2xz_Px_d;
  Double I_ERI_D2z_Py_F2xz_Px_d = I_ERI_Fy2z_S_F2xz_Px_d+ABY*I_ERI_D2z_S_F2xz_Px_d;
  Double I_ERI_D2x_Pz_F2xz_Px_d = I_ERI_F2xz_S_F2xz_Px_d+ABZ*I_ERI_D2x_S_F2xz_Px_d;
  Double I_ERI_Dxy_Pz_F2xz_Px_d = I_ERI_Fxyz_S_F2xz_Px_d+ABZ*I_ERI_Dxy_S_F2xz_Px_d;
  Double I_ERI_Dxz_Pz_F2xz_Px_d = I_ERI_Fx2z_S_F2xz_Px_d+ABZ*I_ERI_Dxz_S_F2xz_Px_d;
  Double I_ERI_D2y_Pz_F2xz_Px_d = I_ERI_F2yz_S_F2xz_Px_d+ABZ*I_ERI_D2y_S_F2xz_Px_d;
  Double I_ERI_Dyz_Pz_F2xz_Px_d = I_ERI_Fy2z_S_F2xz_Px_d+ABZ*I_ERI_Dyz_S_F2xz_Px_d;
  Double I_ERI_D2z_Pz_F2xz_Px_d = I_ERI_F3z_S_F2xz_Px_d+ABZ*I_ERI_D2z_S_F2xz_Px_d;
  Double I_ERI_D2x_Px_Fx2y_Px_d = I_ERI_F3x_S_Fx2y_Px_d+ABX*I_ERI_D2x_S_Fx2y_Px_d;
  Double I_ERI_Dxy_Px_Fx2y_Px_d = I_ERI_F2xy_S_Fx2y_Px_d+ABX*I_ERI_Dxy_S_Fx2y_Px_d;
  Double I_ERI_Dxz_Px_Fx2y_Px_d = I_ERI_F2xz_S_Fx2y_Px_d+ABX*I_ERI_Dxz_S_Fx2y_Px_d;
  Double I_ERI_D2y_Px_Fx2y_Px_d = I_ERI_Fx2y_S_Fx2y_Px_d+ABX*I_ERI_D2y_S_Fx2y_Px_d;
  Double I_ERI_Dyz_Px_Fx2y_Px_d = I_ERI_Fxyz_S_Fx2y_Px_d+ABX*I_ERI_Dyz_S_Fx2y_Px_d;
  Double I_ERI_D2z_Px_Fx2y_Px_d = I_ERI_Fx2z_S_Fx2y_Px_d+ABX*I_ERI_D2z_S_Fx2y_Px_d;
  Double I_ERI_D2x_Py_Fx2y_Px_d = I_ERI_F2xy_S_Fx2y_Px_d+ABY*I_ERI_D2x_S_Fx2y_Px_d;
  Double I_ERI_Dxy_Py_Fx2y_Px_d = I_ERI_Fx2y_S_Fx2y_Px_d+ABY*I_ERI_Dxy_S_Fx2y_Px_d;
  Double I_ERI_Dxz_Py_Fx2y_Px_d = I_ERI_Fxyz_S_Fx2y_Px_d+ABY*I_ERI_Dxz_S_Fx2y_Px_d;
  Double I_ERI_D2y_Py_Fx2y_Px_d = I_ERI_F3y_S_Fx2y_Px_d+ABY*I_ERI_D2y_S_Fx2y_Px_d;
  Double I_ERI_Dyz_Py_Fx2y_Px_d = I_ERI_F2yz_S_Fx2y_Px_d+ABY*I_ERI_Dyz_S_Fx2y_Px_d;
  Double I_ERI_D2z_Py_Fx2y_Px_d = I_ERI_Fy2z_S_Fx2y_Px_d+ABY*I_ERI_D2z_S_Fx2y_Px_d;
  Double I_ERI_D2x_Pz_Fx2y_Px_d = I_ERI_F2xz_S_Fx2y_Px_d+ABZ*I_ERI_D2x_S_Fx2y_Px_d;
  Double I_ERI_Dxy_Pz_Fx2y_Px_d = I_ERI_Fxyz_S_Fx2y_Px_d+ABZ*I_ERI_Dxy_S_Fx2y_Px_d;
  Double I_ERI_Dxz_Pz_Fx2y_Px_d = I_ERI_Fx2z_S_Fx2y_Px_d+ABZ*I_ERI_Dxz_S_Fx2y_Px_d;
  Double I_ERI_D2y_Pz_Fx2y_Px_d = I_ERI_F2yz_S_Fx2y_Px_d+ABZ*I_ERI_D2y_S_Fx2y_Px_d;
  Double I_ERI_Dyz_Pz_Fx2y_Px_d = I_ERI_Fy2z_S_Fx2y_Px_d+ABZ*I_ERI_Dyz_S_Fx2y_Px_d;
  Double I_ERI_D2z_Pz_Fx2y_Px_d = I_ERI_F3z_S_Fx2y_Px_d+ABZ*I_ERI_D2z_S_Fx2y_Px_d;
  Double I_ERI_D2x_Px_Fxyz_Px_d = I_ERI_F3x_S_Fxyz_Px_d+ABX*I_ERI_D2x_S_Fxyz_Px_d;
  Double I_ERI_Dxy_Px_Fxyz_Px_d = I_ERI_F2xy_S_Fxyz_Px_d+ABX*I_ERI_Dxy_S_Fxyz_Px_d;
  Double I_ERI_Dxz_Px_Fxyz_Px_d = I_ERI_F2xz_S_Fxyz_Px_d+ABX*I_ERI_Dxz_S_Fxyz_Px_d;
  Double I_ERI_D2y_Px_Fxyz_Px_d = I_ERI_Fx2y_S_Fxyz_Px_d+ABX*I_ERI_D2y_S_Fxyz_Px_d;
  Double I_ERI_Dyz_Px_Fxyz_Px_d = I_ERI_Fxyz_S_Fxyz_Px_d+ABX*I_ERI_Dyz_S_Fxyz_Px_d;
  Double I_ERI_D2z_Px_Fxyz_Px_d = I_ERI_Fx2z_S_Fxyz_Px_d+ABX*I_ERI_D2z_S_Fxyz_Px_d;
  Double I_ERI_D2x_Py_Fxyz_Px_d = I_ERI_F2xy_S_Fxyz_Px_d+ABY*I_ERI_D2x_S_Fxyz_Px_d;
  Double I_ERI_Dxy_Py_Fxyz_Px_d = I_ERI_Fx2y_S_Fxyz_Px_d+ABY*I_ERI_Dxy_S_Fxyz_Px_d;
  Double I_ERI_Dxz_Py_Fxyz_Px_d = I_ERI_Fxyz_S_Fxyz_Px_d+ABY*I_ERI_Dxz_S_Fxyz_Px_d;
  Double I_ERI_D2y_Py_Fxyz_Px_d = I_ERI_F3y_S_Fxyz_Px_d+ABY*I_ERI_D2y_S_Fxyz_Px_d;
  Double I_ERI_Dyz_Py_Fxyz_Px_d = I_ERI_F2yz_S_Fxyz_Px_d+ABY*I_ERI_Dyz_S_Fxyz_Px_d;
  Double I_ERI_D2z_Py_Fxyz_Px_d = I_ERI_Fy2z_S_Fxyz_Px_d+ABY*I_ERI_D2z_S_Fxyz_Px_d;
  Double I_ERI_D2x_Pz_Fxyz_Px_d = I_ERI_F2xz_S_Fxyz_Px_d+ABZ*I_ERI_D2x_S_Fxyz_Px_d;
  Double I_ERI_Dxy_Pz_Fxyz_Px_d = I_ERI_Fxyz_S_Fxyz_Px_d+ABZ*I_ERI_Dxy_S_Fxyz_Px_d;
  Double I_ERI_Dxz_Pz_Fxyz_Px_d = I_ERI_Fx2z_S_Fxyz_Px_d+ABZ*I_ERI_Dxz_S_Fxyz_Px_d;
  Double I_ERI_D2y_Pz_Fxyz_Px_d = I_ERI_F2yz_S_Fxyz_Px_d+ABZ*I_ERI_D2y_S_Fxyz_Px_d;
  Double I_ERI_Dyz_Pz_Fxyz_Px_d = I_ERI_Fy2z_S_Fxyz_Px_d+ABZ*I_ERI_Dyz_S_Fxyz_Px_d;
  Double I_ERI_D2z_Pz_Fxyz_Px_d = I_ERI_F3z_S_Fxyz_Px_d+ABZ*I_ERI_D2z_S_Fxyz_Px_d;
  Double I_ERI_D2x_Px_Fx2z_Px_d = I_ERI_F3x_S_Fx2z_Px_d+ABX*I_ERI_D2x_S_Fx2z_Px_d;
  Double I_ERI_Dxy_Px_Fx2z_Px_d = I_ERI_F2xy_S_Fx2z_Px_d+ABX*I_ERI_Dxy_S_Fx2z_Px_d;
  Double I_ERI_Dxz_Px_Fx2z_Px_d = I_ERI_F2xz_S_Fx2z_Px_d+ABX*I_ERI_Dxz_S_Fx2z_Px_d;
  Double I_ERI_D2y_Px_Fx2z_Px_d = I_ERI_Fx2y_S_Fx2z_Px_d+ABX*I_ERI_D2y_S_Fx2z_Px_d;
  Double I_ERI_Dyz_Px_Fx2z_Px_d = I_ERI_Fxyz_S_Fx2z_Px_d+ABX*I_ERI_Dyz_S_Fx2z_Px_d;
  Double I_ERI_D2z_Px_Fx2z_Px_d = I_ERI_Fx2z_S_Fx2z_Px_d+ABX*I_ERI_D2z_S_Fx2z_Px_d;
  Double I_ERI_D2x_Py_Fx2z_Px_d = I_ERI_F2xy_S_Fx2z_Px_d+ABY*I_ERI_D2x_S_Fx2z_Px_d;
  Double I_ERI_Dxy_Py_Fx2z_Px_d = I_ERI_Fx2y_S_Fx2z_Px_d+ABY*I_ERI_Dxy_S_Fx2z_Px_d;
  Double I_ERI_Dxz_Py_Fx2z_Px_d = I_ERI_Fxyz_S_Fx2z_Px_d+ABY*I_ERI_Dxz_S_Fx2z_Px_d;
  Double I_ERI_D2y_Py_Fx2z_Px_d = I_ERI_F3y_S_Fx2z_Px_d+ABY*I_ERI_D2y_S_Fx2z_Px_d;
  Double I_ERI_Dyz_Py_Fx2z_Px_d = I_ERI_F2yz_S_Fx2z_Px_d+ABY*I_ERI_Dyz_S_Fx2z_Px_d;
  Double I_ERI_D2z_Py_Fx2z_Px_d = I_ERI_Fy2z_S_Fx2z_Px_d+ABY*I_ERI_D2z_S_Fx2z_Px_d;
  Double I_ERI_D2x_Pz_Fx2z_Px_d = I_ERI_F2xz_S_Fx2z_Px_d+ABZ*I_ERI_D2x_S_Fx2z_Px_d;
  Double I_ERI_Dxy_Pz_Fx2z_Px_d = I_ERI_Fxyz_S_Fx2z_Px_d+ABZ*I_ERI_Dxy_S_Fx2z_Px_d;
  Double I_ERI_Dxz_Pz_Fx2z_Px_d = I_ERI_Fx2z_S_Fx2z_Px_d+ABZ*I_ERI_Dxz_S_Fx2z_Px_d;
  Double I_ERI_D2y_Pz_Fx2z_Px_d = I_ERI_F2yz_S_Fx2z_Px_d+ABZ*I_ERI_D2y_S_Fx2z_Px_d;
  Double I_ERI_Dyz_Pz_Fx2z_Px_d = I_ERI_Fy2z_S_Fx2z_Px_d+ABZ*I_ERI_Dyz_S_Fx2z_Px_d;
  Double I_ERI_D2z_Pz_Fx2z_Px_d = I_ERI_F3z_S_Fx2z_Px_d+ABZ*I_ERI_D2z_S_Fx2z_Px_d;
  Double I_ERI_D2x_Px_F3y_Px_d = I_ERI_F3x_S_F3y_Px_d+ABX*I_ERI_D2x_S_F3y_Px_d;
  Double I_ERI_Dxy_Px_F3y_Px_d = I_ERI_F2xy_S_F3y_Px_d+ABX*I_ERI_Dxy_S_F3y_Px_d;
  Double I_ERI_Dxz_Px_F3y_Px_d = I_ERI_F2xz_S_F3y_Px_d+ABX*I_ERI_Dxz_S_F3y_Px_d;
  Double I_ERI_D2y_Px_F3y_Px_d = I_ERI_Fx2y_S_F3y_Px_d+ABX*I_ERI_D2y_S_F3y_Px_d;
  Double I_ERI_Dyz_Px_F3y_Px_d = I_ERI_Fxyz_S_F3y_Px_d+ABX*I_ERI_Dyz_S_F3y_Px_d;
  Double I_ERI_D2z_Px_F3y_Px_d = I_ERI_Fx2z_S_F3y_Px_d+ABX*I_ERI_D2z_S_F3y_Px_d;
  Double I_ERI_D2x_Py_F3y_Px_d = I_ERI_F2xy_S_F3y_Px_d+ABY*I_ERI_D2x_S_F3y_Px_d;
  Double I_ERI_Dxy_Py_F3y_Px_d = I_ERI_Fx2y_S_F3y_Px_d+ABY*I_ERI_Dxy_S_F3y_Px_d;
  Double I_ERI_Dxz_Py_F3y_Px_d = I_ERI_Fxyz_S_F3y_Px_d+ABY*I_ERI_Dxz_S_F3y_Px_d;
  Double I_ERI_D2y_Py_F3y_Px_d = I_ERI_F3y_S_F3y_Px_d+ABY*I_ERI_D2y_S_F3y_Px_d;
  Double I_ERI_Dyz_Py_F3y_Px_d = I_ERI_F2yz_S_F3y_Px_d+ABY*I_ERI_Dyz_S_F3y_Px_d;
  Double I_ERI_D2z_Py_F3y_Px_d = I_ERI_Fy2z_S_F3y_Px_d+ABY*I_ERI_D2z_S_F3y_Px_d;
  Double I_ERI_D2x_Pz_F3y_Px_d = I_ERI_F2xz_S_F3y_Px_d+ABZ*I_ERI_D2x_S_F3y_Px_d;
  Double I_ERI_Dxy_Pz_F3y_Px_d = I_ERI_Fxyz_S_F3y_Px_d+ABZ*I_ERI_Dxy_S_F3y_Px_d;
  Double I_ERI_Dxz_Pz_F3y_Px_d = I_ERI_Fx2z_S_F3y_Px_d+ABZ*I_ERI_Dxz_S_F3y_Px_d;
  Double I_ERI_D2y_Pz_F3y_Px_d = I_ERI_F2yz_S_F3y_Px_d+ABZ*I_ERI_D2y_S_F3y_Px_d;
  Double I_ERI_Dyz_Pz_F3y_Px_d = I_ERI_Fy2z_S_F3y_Px_d+ABZ*I_ERI_Dyz_S_F3y_Px_d;
  Double I_ERI_D2z_Pz_F3y_Px_d = I_ERI_F3z_S_F3y_Px_d+ABZ*I_ERI_D2z_S_F3y_Px_d;
  Double I_ERI_D2x_Px_F2yz_Px_d = I_ERI_F3x_S_F2yz_Px_d+ABX*I_ERI_D2x_S_F2yz_Px_d;
  Double I_ERI_Dxy_Px_F2yz_Px_d = I_ERI_F2xy_S_F2yz_Px_d+ABX*I_ERI_Dxy_S_F2yz_Px_d;
  Double I_ERI_Dxz_Px_F2yz_Px_d = I_ERI_F2xz_S_F2yz_Px_d+ABX*I_ERI_Dxz_S_F2yz_Px_d;
  Double I_ERI_D2y_Px_F2yz_Px_d = I_ERI_Fx2y_S_F2yz_Px_d+ABX*I_ERI_D2y_S_F2yz_Px_d;
  Double I_ERI_Dyz_Px_F2yz_Px_d = I_ERI_Fxyz_S_F2yz_Px_d+ABX*I_ERI_Dyz_S_F2yz_Px_d;
  Double I_ERI_D2z_Px_F2yz_Px_d = I_ERI_Fx2z_S_F2yz_Px_d+ABX*I_ERI_D2z_S_F2yz_Px_d;
  Double I_ERI_D2x_Py_F2yz_Px_d = I_ERI_F2xy_S_F2yz_Px_d+ABY*I_ERI_D2x_S_F2yz_Px_d;
  Double I_ERI_Dxy_Py_F2yz_Px_d = I_ERI_Fx2y_S_F2yz_Px_d+ABY*I_ERI_Dxy_S_F2yz_Px_d;
  Double I_ERI_Dxz_Py_F2yz_Px_d = I_ERI_Fxyz_S_F2yz_Px_d+ABY*I_ERI_Dxz_S_F2yz_Px_d;
  Double I_ERI_D2y_Py_F2yz_Px_d = I_ERI_F3y_S_F2yz_Px_d+ABY*I_ERI_D2y_S_F2yz_Px_d;
  Double I_ERI_Dyz_Py_F2yz_Px_d = I_ERI_F2yz_S_F2yz_Px_d+ABY*I_ERI_Dyz_S_F2yz_Px_d;
  Double I_ERI_D2z_Py_F2yz_Px_d = I_ERI_Fy2z_S_F2yz_Px_d+ABY*I_ERI_D2z_S_F2yz_Px_d;
  Double I_ERI_D2x_Pz_F2yz_Px_d = I_ERI_F2xz_S_F2yz_Px_d+ABZ*I_ERI_D2x_S_F2yz_Px_d;
  Double I_ERI_Dxy_Pz_F2yz_Px_d = I_ERI_Fxyz_S_F2yz_Px_d+ABZ*I_ERI_Dxy_S_F2yz_Px_d;
  Double I_ERI_Dxz_Pz_F2yz_Px_d = I_ERI_Fx2z_S_F2yz_Px_d+ABZ*I_ERI_Dxz_S_F2yz_Px_d;
  Double I_ERI_D2y_Pz_F2yz_Px_d = I_ERI_F2yz_S_F2yz_Px_d+ABZ*I_ERI_D2y_S_F2yz_Px_d;
  Double I_ERI_Dyz_Pz_F2yz_Px_d = I_ERI_Fy2z_S_F2yz_Px_d+ABZ*I_ERI_Dyz_S_F2yz_Px_d;
  Double I_ERI_D2z_Pz_F2yz_Px_d = I_ERI_F3z_S_F2yz_Px_d+ABZ*I_ERI_D2z_S_F2yz_Px_d;
  Double I_ERI_D2x_Px_Fy2z_Px_d = I_ERI_F3x_S_Fy2z_Px_d+ABX*I_ERI_D2x_S_Fy2z_Px_d;
  Double I_ERI_Dxy_Px_Fy2z_Px_d = I_ERI_F2xy_S_Fy2z_Px_d+ABX*I_ERI_Dxy_S_Fy2z_Px_d;
  Double I_ERI_Dxz_Px_Fy2z_Px_d = I_ERI_F2xz_S_Fy2z_Px_d+ABX*I_ERI_Dxz_S_Fy2z_Px_d;
  Double I_ERI_D2y_Px_Fy2z_Px_d = I_ERI_Fx2y_S_Fy2z_Px_d+ABX*I_ERI_D2y_S_Fy2z_Px_d;
  Double I_ERI_Dyz_Px_Fy2z_Px_d = I_ERI_Fxyz_S_Fy2z_Px_d+ABX*I_ERI_Dyz_S_Fy2z_Px_d;
  Double I_ERI_D2z_Px_Fy2z_Px_d = I_ERI_Fx2z_S_Fy2z_Px_d+ABX*I_ERI_D2z_S_Fy2z_Px_d;
  Double I_ERI_D2x_Py_Fy2z_Px_d = I_ERI_F2xy_S_Fy2z_Px_d+ABY*I_ERI_D2x_S_Fy2z_Px_d;
  Double I_ERI_Dxy_Py_Fy2z_Px_d = I_ERI_Fx2y_S_Fy2z_Px_d+ABY*I_ERI_Dxy_S_Fy2z_Px_d;
  Double I_ERI_Dxz_Py_Fy2z_Px_d = I_ERI_Fxyz_S_Fy2z_Px_d+ABY*I_ERI_Dxz_S_Fy2z_Px_d;
  Double I_ERI_D2y_Py_Fy2z_Px_d = I_ERI_F3y_S_Fy2z_Px_d+ABY*I_ERI_D2y_S_Fy2z_Px_d;
  Double I_ERI_Dyz_Py_Fy2z_Px_d = I_ERI_F2yz_S_Fy2z_Px_d+ABY*I_ERI_Dyz_S_Fy2z_Px_d;
  Double I_ERI_D2z_Py_Fy2z_Px_d = I_ERI_Fy2z_S_Fy2z_Px_d+ABY*I_ERI_D2z_S_Fy2z_Px_d;
  Double I_ERI_D2x_Pz_Fy2z_Px_d = I_ERI_F2xz_S_Fy2z_Px_d+ABZ*I_ERI_D2x_S_Fy2z_Px_d;
  Double I_ERI_Dxy_Pz_Fy2z_Px_d = I_ERI_Fxyz_S_Fy2z_Px_d+ABZ*I_ERI_Dxy_S_Fy2z_Px_d;
  Double I_ERI_Dxz_Pz_Fy2z_Px_d = I_ERI_Fx2z_S_Fy2z_Px_d+ABZ*I_ERI_Dxz_S_Fy2z_Px_d;
  Double I_ERI_D2y_Pz_Fy2z_Px_d = I_ERI_F2yz_S_Fy2z_Px_d+ABZ*I_ERI_D2y_S_Fy2z_Px_d;
  Double I_ERI_Dyz_Pz_Fy2z_Px_d = I_ERI_Fy2z_S_Fy2z_Px_d+ABZ*I_ERI_Dyz_S_Fy2z_Px_d;
  Double I_ERI_D2z_Pz_Fy2z_Px_d = I_ERI_F3z_S_Fy2z_Px_d+ABZ*I_ERI_D2z_S_Fy2z_Px_d;
  Double I_ERI_D2x_Px_F3z_Px_d = I_ERI_F3x_S_F3z_Px_d+ABX*I_ERI_D2x_S_F3z_Px_d;
  Double I_ERI_Dxy_Px_F3z_Px_d = I_ERI_F2xy_S_F3z_Px_d+ABX*I_ERI_Dxy_S_F3z_Px_d;
  Double I_ERI_Dxz_Px_F3z_Px_d = I_ERI_F2xz_S_F3z_Px_d+ABX*I_ERI_Dxz_S_F3z_Px_d;
  Double I_ERI_D2y_Px_F3z_Px_d = I_ERI_Fx2y_S_F3z_Px_d+ABX*I_ERI_D2y_S_F3z_Px_d;
  Double I_ERI_Dyz_Px_F3z_Px_d = I_ERI_Fxyz_S_F3z_Px_d+ABX*I_ERI_Dyz_S_F3z_Px_d;
  Double I_ERI_D2z_Px_F3z_Px_d = I_ERI_Fx2z_S_F3z_Px_d+ABX*I_ERI_D2z_S_F3z_Px_d;
  Double I_ERI_D2x_Py_F3z_Px_d = I_ERI_F2xy_S_F3z_Px_d+ABY*I_ERI_D2x_S_F3z_Px_d;
  Double I_ERI_Dxy_Py_F3z_Px_d = I_ERI_Fx2y_S_F3z_Px_d+ABY*I_ERI_Dxy_S_F3z_Px_d;
  Double I_ERI_Dxz_Py_F3z_Px_d = I_ERI_Fxyz_S_F3z_Px_d+ABY*I_ERI_Dxz_S_F3z_Px_d;
  Double I_ERI_D2y_Py_F3z_Px_d = I_ERI_F3y_S_F3z_Px_d+ABY*I_ERI_D2y_S_F3z_Px_d;
  Double I_ERI_Dyz_Py_F3z_Px_d = I_ERI_F2yz_S_F3z_Px_d+ABY*I_ERI_Dyz_S_F3z_Px_d;
  Double I_ERI_D2z_Py_F3z_Px_d = I_ERI_Fy2z_S_F3z_Px_d+ABY*I_ERI_D2z_S_F3z_Px_d;
  Double I_ERI_D2x_Pz_F3z_Px_d = I_ERI_F2xz_S_F3z_Px_d+ABZ*I_ERI_D2x_S_F3z_Px_d;
  Double I_ERI_Dxy_Pz_F3z_Px_d = I_ERI_Fxyz_S_F3z_Px_d+ABZ*I_ERI_Dxy_S_F3z_Px_d;
  Double I_ERI_Dxz_Pz_F3z_Px_d = I_ERI_Fx2z_S_F3z_Px_d+ABZ*I_ERI_Dxz_S_F3z_Px_d;
  Double I_ERI_D2y_Pz_F3z_Px_d = I_ERI_F2yz_S_F3z_Px_d+ABZ*I_ERI_D2y_S_F3z_Px_d;
  Double I_ERI_Dyz_Pz_F3z_Px_d = I_ERI_Fy2z_S_F3z_Px_d+ABZ*I_ERI_Dyz_S_F3z_Px_d;
  Double I_ERI_D2z_Pz_F3z_Px_d = I_ERI_F3z_S_F3z_Px_d+ABZ*I_ERI_D2z_S_F3z_Px_d;
  Double I_ERI_D2x_Px_F3x_Py_d = I_ERI_F3x_S_F3x_Py_d+ABX*I_ERI_D2x_S_F3x_Py_d;
  Double I_ERI_Dxy_Px_F3x_Py_d = I_ERI_F2xy_S_F3x_Py_d+ABX*I_ERI_Dxy_S_F3x_Py_d;
  Double I_ERI_Dxz_Px_F3x_Py_d = I_ERI_F2xz_S_F3x_Py_d+ABX*I_ERI_Dxz_S_F3x_Py_d;
  Double I_ERI_D2y_Px_F3x_Py_d = I_ERI_Fx2y_S_F3x_Py_d+ABX*I_ERI_D2y_S_F3x_Py_d;
  Double I_ERI_Dyz_Px_F3x_Py_d = I_ERI_Fxyz_S_F3x_Py_d+ABX*I_ERI_Dyz_S_F3x_Py_d;
  Double I_ERI_D2z_Px_F3x_Py_d = I_ERI_Fx2z_S_F3x_Py_d+ABX*I_ERI_D2z_S_F3x_Py_d;
  Double I_ERI_D2x_Py_F3x_Py_d = I_ERI_F2xy_S_F3x_Py_d+ABY*I_ERI_D2x_S_F3x_Py_d;
  Double I_ERI_Dxy_Py_F3x_Py_d = I_ERI_Fx2y_S_F3x_Py_d+ABY*I_ERI_Dxy_S_F3x_Py_d;
  Double I_ERI_Dxz_Py_F3x_Py_d = I_ERI_Fxyz_S_F3x_Py_d+ABY*I_ERI_Dxz_S_F3x_Py_d;
  Double I_ERI_D2y_Py_F3x_Py_d = I_ERI_F3y_S_F3x_Py_d+ABY*I_ERI_D2y_S_F3x_Py_d;
  Double I_ERI_Dyz_Py_F3x_Py_d = I_ERI_F2yz_S_F3x_Py_d+ABY*I_ERI_Dyz_S_F3x_Py_d;
  Double I_ERI_D2z_Py_F3x_Py_d = I_ERI_Fy2z_S_F3x_Py_d+ABY*I_ERI_D2z_S_F3x_Py_d;
  Double I_ERI_D2x_Pz_F3x_Py_d = I_ERI_F2xz_S_F3x_Py_d+ABZ*I_ERI_D2x_S_F3x_Py_d;
  Double I_ERI_Dxy_Pz_F3x_Py_d = I_ERI_Fxyz_S_F3x_Py_d+ABZ*I_ERI_Dxy_S_F3x_Py_d;
  Double I_ERI_Dxz_Pz_F3x_Py_d = I_ERI_Fx2z_S_F3x_Py_d+ABZ*I_ERI_Dxz_S_F3x_Py_d;
  Double I_ERI_D2y_Pz_F3x_Py_d = I_ERI_F2yz_S_F3x_Py_d+ABZ*I_ERI_D2y_S_F3x_Py_d;
  Double I_ERI_Dyz_Pz_F3x_Py_d = I_ERI_Fy2z_S_F3x_Py_d+ABZ*I_ERI_Dyz_S_F3x_Py_d;
  Double I_ERI_D2z_Pz_F3x_Py_d = I_ERI_F3z_S_F3x_Py_d+ABZ*I_ERI_D2z_S_F3x_Py_d;
  Double I_ERI_D2x_Px_F2xy_Py_d = I_ERI_F3x_S_F2xy_Py_d+ABX*I_ERI_D2x_S_F2xy_Py_d;
  Double I_ERI_Dxy_Px_F2xy_Py_d = I_ERI_F2xy_S_F2xy_Py_d+ABX*I_ERI_Dxy_S_F2xy_Py_d;
  Double I_ERI_Dxz_Px_F2xy_Py_d = I_ERI_F2xz_S_F2xy_Py_d+ABX*I_ERI_Dxz_S_F2xy_Py_d;
  Double I_ERI_D2y_Px_F2xy_Py_d = I_ERI_Fx2y_S_F2xy_Py_d+ABX*I_ERI_D2y_S_F2xy_Py_d;
  Double I_ERI_Dyz_Px_F2xy_Py_d = I_ERI_Fxyz_S_F2xy_Py_d+ABX*I_ERI_Dyz_S_F2xy_Py_d;
  Double I_ERI_D2z_Px_F2xy_Py_d = I_ERI_Fx2z_S_F2xy_Py_d+ABX*I_ERI_D2z_S_F2xy_Py_d;
  Double I_ERI_D2x_Py_F2xy_Py_d = I_ERI_F2xy_S_F2xy_Py_d+ABY*I_ERI_D2x_S_F2xy_Py_d;
  Double I_ERI_Dxy_Py_F2xy_Py_d = I_ERI_Fx2y_S_F2xy_Py_d+ABY*I_ERI_Dxy_S_F2xy_Py_d;
  Double I_ERI_Dxz_Py_F2xy_Py_d = I_ERI_Fxyz_S_F2xy_Py_d+ABY*I_ERI_Dxz_S_F2xy_Py_d;
  Double I_ERI_D2y_Py_F2xy_Py_d = I_ERI_F3y_S_F2xy_Py_d+ABY*I_ERI_D2y_S_F2xy_Py_d;
  Double I_ERI_Dyz_Py_F2xy_Py_d = I_ERI_F2yz_S_F2xy_Py_d+ABY*I_ERI_Dyz_S_F2xy_Py_d;
  Double I_ERI_D2z_Py_F2xy_Py_d = I_ERI_Fy2z_S_F2xy_Py_d+ABY*I_ERI_D2z_S_F2xy_Py_d;
  Double I_ERI_D2x_Pz_F2xy_Py_d = I_ERI_F2xz_S_F2xy_Py_d+ABZ*I_ERI_D2x_S_F2xy_Py_d;
  Double I_ERI_Dxy_Pz_F2xy_Py_d = I_ERI_Fxyz_S_F2xy_Py_d+ABZ*I_ERI_Dxy_S_F2xy_Py_d;
  Double I_ERI_Dxz_Pz_F2xy_Py_d = I_ERI_Fx2z_S_F2xy_Py_d+ABZ*I_ERI_Dxz_S_F2xy_Py_d;
  Double I_ERI_D2y_Pz_F2xy_Py_d = I_ERI_F2yz_S_F2xy_Py_d+ABZ*I_ERI_D2y_S_F2xy_Py_d;
  Double I_ERI_Dyz_Pz_F2xy_Py_d = I_ERI_Fy2z_S_F2xy_Py_d+ABZ*I_ERI_Dyz_S_F2xy_Py_d;
  Double I_ERI_D2z_Pz_F2xy_Py_d = I_ERI_F3z_S_F2xy_Py_d+ABZ*I_ERI_D2z_S_F2xy_Py_d;
  Double I_ERI_D2x_Px_F2xz_Py_d = I_ERI_F3x_S_F2xz_Py_d+ABX*I_ERI_D2x_S_F2xz_Py_d;
  Double I_ERI_Dxy_Px_F2xz_Py_d = I_ERI_F2xy_S_F2xz_Py_d+ABX*I_ERI_Dxy_S_F2xz_Py_d;
  Double I_ERI_Dxz_Px_F2xz_Py_d = I_ERI_F2xz_S_F2xz_Py_d+ABX*I_ERI_Dxz_S_F2xz_Py_d;
  Double I_ERI_D2y_Px_F2xz_Py_d = I_ERI_Fx2y_S_F2xz_Py_d+ABX*I_ERI_D2y_S_F2xz_Py_d;
  Double I_ERI_Dyz_Px_F2xz_Py_d = I_ERI_Fxyz_S_F2xz_Py_d+ABX*I_ERI_Dyz_S_F2xz_Py_d;
  Double I_ERI_D2z_Px_F2xz_Py_d = I_ERI_Fx2z_S_F2xz_Py_d+ABX*I_ERI_D2z_S_F2xz_Py_d;
  Double I_ERI_D2x_Py_F2xz_Py_d = I_ERI_F2xy_S_F2xz_Py_d+ABY*I_ERI_D2x_S_F2xz_Py_d;
  Double I_ERI_Dxy_Py_F2xz_Py_d = I_ERI_Fx2y_S_F2xz_Py_d+ABY*I_ERI_Dxy_S_F2xz_Py_d;
  Double I_ERI_Dxz_Py_F2xz_Py_d = I_ERI_Fxyz_S_F2xz_Py_d+ABY*I_ERI_Dxz_S_F2xz_Py_d;
  Double I_ERI_D2y_Py_F2xz_Py_d = I_ERI_F3y_S_F2xz_Py_d+ABY*I_ERI_D2y_S_F2xz_Py_d;
  Double I_ERI_Dyz_Py_F2xz_Py_d = I_ERI_F2yz_S_F2xz_Py_d+ABY*I_ERI_Dyz_S_F2xz_Py_d;
  Double I_ERI_D2z_Py_F2xz_Py_d = I_ERI_Fy2z_S_F2xz_Py_d+ABY*I_ERI_D2z_S_F2xz_Py_d;
  Double I_ERI_D2x_Pz_F2xz_Py_d = I_ERI_F2xz_S_F2xz_Py_d+ABZ*I_ERI_D2x_S_F2xz_Py_d;
  Double I_ERI_Dxy_Pz_F2xz_Py_d = I_ERI_Fxyz_S_F2xz_Py_d+ABZ*I_ERI_Dxy_S_F2xz_Py_d;
  Double I_ERI_Dxz_Pz_F2xz_Py_d = I_ERI_Fx2z_S_F2xz_Py_d+ABZ*I_ERI_Dxz_S_F2xz_Py_d;
  Double I_ERI_D2y_Pz_F2xz_Py_d = I_ERI_F2yz_S_F2xz_Py_d+ABZ*I_ERI_D2y_S_F2xz_Py_d;
  Double I_ERI_Dyz_Pz_F2xz_Py_d = I_ERI_Fy2z_S_F2xz_Py_d+ABZ*I_ERI_Dyz_S_F2xz_Py_d;
  Double I_ERI_D2z_Pz_F2xz_Py_d = I_ERI_F3z_S_F2xz_Py_d+ABZ*I_ERI_D2z_S_F2xz_Py_d;
  Double I_ERI_D2x_Px_Fx2y_Py_d = I_ERI_F3x_S_Fx2y_Py_d+ABX*I_ERI_D2x_S_Fx2y_Py_d;
  Double I_ERI_Dxy_Px_Fx2y_Py_d = I_ERI_F2xy_S_Fx2y_Py_d+ABX*I_ERI_Dxy_S_Fx2y_Py_d;
  Double I_ERI_Dxz_Px_Fx2y_Py_d = I_ERI_F2xz_S_Fx2y_Py_d+ABX*I_ERI_Dxz_S_Fx2y_Py_d;
  Double I_ERI_D2y_Px_Fx2y_Py_d = I_ERI_Fx2y_S_Fx2y_Py_d+ABX*I_ERI_D2y_S_Fx2y_Py_d;
  Double I_ERI_Dyz_Px_Fx2y_Py_d = I_ERI_Fxyz_S_Fx2y_Py_d+ABX*I_ERI_Dyz_S_Fx2y_Py_d;
  Double I_ERI_D2z_Px_Fx2y_Py_d = I_ERI_Fx2z_S_Fx2y_Py_d+ABX*I_ERI_D2z_S_Fx2y_Py_d;
  Double I_ERI_D2x_Py_Fx2y_Py_d = I_ERI_F2xy_S_Fx2y_Py_d+ABY*I_ERI_D2x_S_Fx2y_Py_d;
  Double I_ERI_Dxy_Py_Fx2y_Py_d = I_ERI_Fx2y_S_Fx2y_Py_d+ABY*I_ERI_Dxy_S_Fx2y_Py_d;
  Double I_ERI_Dxz_Py_Fx2y_Py_d = I_ERI_Fxyz_S_Fx2y_Py_d+ABY*I_ERI_Dxz_S_Fx2y_Py_d;
  Double I_ERI_D2y_Py_Fx2y_Py_d = I_ERI_F3y_S_Fx2y_Py_d+ABY*I_ERI_D2y_S_Fx2y_Py_d;
  Double I_ERI_Dyz_Py_Fx2y_Py_d = I_ERI_F2yz_S_Fx2y_Py_d+ABY*I_ERI_Dyz_S_Fx2y_Py_d;
  Double I_ERI_D2z_Py_Fx2y_Py_d = I_ERI_Fy2z_S_Fx2y_Py_d+ABY*I_ERI_D2z_S_Fx2y_Py_d;
  Double I_ERI_D2x_Pz_Fx2y_Py_d = I_ERI_F2xz_S_Fx2y_Py_d+ABZ*I_ERI_D2x_S_Fx2y_Py_d;
  Double I_ERI_Dxy_Pz_Fx2y_Py_d = I_ERI_Fxyz_S_Fx2y_Py_d+ABZ*I_ERI_Dxy_S_Fx2y_Py_d;
  Double I_ERI_Dxz_Pz_Fx2y_Py_d = I_ERI_Fx2z_S_Fx2y_Py_d+ABZ*I_ERI_Dxz_S_Fx2y_Py_d;
  Double I_ERI_D2y_Pz_Fx2y_Py_d = I_ERI_F2yz_S_Fx2y_Py_d+ABZ*I_ERI_D2y_S_Fx2y_Py_d;
  Double I_ERI_Dyz_Pz_Fx2y_Py_d = I_ERI_Fy2z_S_Fx2y_Py_d+ABZ*I_ERI_Dyz_S_Fx2y_Py_d;
  Double I_ERI_D2z_Pz_Fx2y_Py_d = I_ERI_F3z_S_Fx2y_Py_d+ABZ*I_ERI_D2z_S_Fx2y_Py_d;
  Double I_ERI_D2x_Px_Fxyz_Py_d = I_ERI_F3x_S_Fxyz_Py_d+ABX*I_ERI_D2x_S_Fxyz_Py_d;
  Double I_ERI_Dxy_Px_Fxyz_Py_d = I_ERI_F2xy_S_Fxyz_Py_d+ABX*I_ERI_Dxy_S_Fxyz_Py_d;
  Double I_ERI_Dxz_Px_Fxyz_Py_d = I_ERI_F2xz_S_Fxyz_Py_d+ABX*I_ERI_Dxz_S_Fxyz_Py_d;
  Double I_ERI_D2y_Px_Fxyz_Py_d = I_ERI_Fx2y_S_Fxyz_Py_d+ABX*I_ERI_D2y_S_Fxyz_Py_d;
  Double I_ERI_Dyz_Px_Fxyz_Py_d = I_ERI_Fxyz_S_Fxyz_Py_d+ABX*I_ERI_Dyz_S_Fxyz_Py_d;
  Double I_ERI_D2z_Px_Fxyz_Py_d = I_ERI_Fx2z_S_Fxyz_Py_d+ABX*I_ERI_D2z_S_Fxyz_Py_d;
  Double I_ERI_D2x_Py_Fxyz_Py_d = I_ERI_F2xy_S_Fxyz_Py_d+ABY*I_ERI_D2x_S_Fxyz_Py_d;
  Double I_ERI_Dxy_Py_Fxyz_Py_d = I_ERI_Fx2y_S_Fxyz_Py_d+ABY*I_ERI_Dxy_S_Fxyz_Py_d;
  Double I_ERI_Dxz_Py_Fxyz_Py_d = I_ERI_Fxyz_S_Fxyz_Py_d+ABY*I_ERI_Dxz_S_Fxyz_Py_d;
  Double I_ERI_D2y_Py_Fxyz_Py_d = I_ERI_F3y_S_Fxyz_Py_d+ABY*I_ERI_D2y_S_Fxyz_Py_d;
  Double I_ERI_Dyz_Py_Fxyz_Py_d = I_ERI_F2yz_S_Fxyz_Py_d+ABY*I_ERI_Dyz_S_Fxyz_Py_d;
  Double I_ERI_D2z_Py_Fxyz_Py_d = I_ERI_Fy2z_S_Fxyz_Py_d+ABY*I_ERI_D2z_S_Fxyz_Py_d;
  Double I_ERI_D2x_Pz_Fxyz_Py_d = I_ERI_F2xz_S_Fxyz_Py_d+ABZ*I_ERI_D2x_S_Fxyz_Py_d;
  Double I_ERI_Dxy_Pz_Fxyz_Py_d = I_ERI_Fxyz_S_Fxyz_Py_d+ABZ*I_ERI_Dxy_S_Fxyz_Py_d;
  Double I_ERI_Dxz_Pz_Fxyz_Py_d = I_ERI_Fx2z_S_Fxyz_Py_d+ABZ*I_ERI_Dxz_S_Fxyz_Py_d;
  Double I_ERI_D2y_Pz_Fxyz_Py_d = I_ERI_F2yz_S_Fxyz_Py_d+ABZ*I_ERI_D2y_S_Fxyz_Py_d;
  Double I_ERI_Dyz_Pz_Fxyz_Py_d = I_ERI_Fy2z_S_Fxyz_Py_d+ABZ*I_ERI_Dyz_S_Fxyz_Py_d;
  Double I_ERI_D2z_Pz_Fxyz_Py_d = I_ERI_F3z_S_Fxyz_Py_d+ABZ*I_ERI_D2z_S_Fxyz_Py_d;
  Double I_ERI_D2x_Px_Fx2z_Py_d = I_ERI_F3x_S_Fx2z_Py_d+ABX*I_ERI_D2x_S_Fx2z_Py_d;
  Double I_ERI_Dxy_Px_Fx2z_Py_d = I_ERI_F2xy_S_Fx2z_Py_d+ABX*I_ERI_Dxy_S_Fx2z_Py_d;
  Double I_ERI_Dxz_Px_Fx2z_Py_d = I_ERI_F2xz_S_Fx2z_Py_d+ABX*I_ERI_Dxz_S_Fx2z_Py_d;
  Double I_ERI_D2y_Px_Fx2z_Py_d = I_ERI_Fx2y_S_Fx2z_Py_d+ABX*I_ERI_D2y_S_Fx2z_Py_d;
  Double I_ERI_Dyz_Px_Fx2z_Py_d = I_ERI_Fxyz_S_Fx2z_Py_d+ABX*I_ERI_Dyz_S_Fx2z_Py_d;
  Double I_ERI_D2z_Px_Fx2z_Py_d = I_ERI_Fx2z_S_Fx2z_Py_d+ABX*I_ERI_D2z_S_Fx2z_Py_d;
  Double I_ERI_D2x_Py_Fx2z_Py_d = I_ERI_F2xy_S_Fx2z_Py_d+ABY*I_ERI_D2x_S_Fx2z_Py_d;
  Double I_ERI_Dxy_Py_Fx2z_Py_d = I_ERI_Fx2y_S_Fx2z_Py_d+ABY*I_ERI_Dxy_S_Fx2z_Py_d;
  Double I_ERI_Dxz_Py_Fx2z_Py_d = I_ERI_Fxyz_S_Fx2z_Py_d+ABY*I_ERI_Dxz_S_Fx2z_Py_d;
  Double I_ERI_D2y_Py_Fx2z_Py_d = I_ERI_F3y_S_Fx2z_Py_d+ABY*I_ERI_D2y_S_Fx2z_Py_d;
  Double I_ERI_Dyz_Py_Fx2z_Py_d = I_ERI_F2yz_S_Fx2z_Py_d+ABY*I_ERI_Dyz_S_Fx2z_Py_d;
  Double I_ERI_D2z_Py_Fx2z_Py_d = I_ERI_Fy2z_S_Fx2z_Py_d+ABY*I_ERI_D2z_S_Fx2z_Py_d;
  Double I_ERI_D2x_Pz_Fx2z_Py_d = I_ERI_F2xz_S_Fx2z_Py_d+ABZ*I_ERI_D2x_S_Fx2z_Py_d;
  Double I_ERI_Dxy_Pz_Fx2z_Py_d = I_ERI_Fxyz_S_Fx2z_Py_d+ABZ*I_ERI_Dxy_S_Fx2z_Py_d;
  Double I_ERI_Dxz_Pz_Fx2z_Py_d = I_ERI_Fx2z_S_Fx2z_Py_d+ABZ*I_ERI_Dxz_S_Fx2z_Py_d;
  Double I_ERI_D2y_Pz_Fx2z_Py_d = I_ERI_F2yz_S_Fx2z_Py_d+ABZ*I_ERI_D2y_S_Fx2z_Py_d;
  Double I_ERI_Dyz_Pz_Fx2z_Py_d = I_ERI_Fy2z_S_Fx2z_Py_d+ABZ*I_ERI_Dyz_S_Fx2z_Py_d;
  Double I_ERI_D2z_Pz_Fx2z_Py_d = I_ERI_F3z_S_Fx2z_Py_d+ABZ*I_ERI_D2z_S_Fx2z_Py_d;
  Double I_ERI_D2x_Px_F3y_Py_d = I_ERI_F3x_S_F3y_Py_d+ABX*I_ERI_D2x_S_F3y_Py_d;
  Double I_ERI_Dxy_Px_F3y_Py_d = I_ERI_F2xy_S_F3y_Py_d+ABX*I_ERI_Dxy_S_F3y_Py_d;
  Double I_ERI_Dxz_Px_F3y_Py_d = I_ERI_F2xz_S_F3y_Py_d+ABX*I_ERI_Dxz_S_F3y_Py_d;
  Double I_ERI_D2y_Px_F3y_Py_d = I_ERI_Fx2y_S_F3y_Py_d+ABX*I_ERI_D2y_S_F3y_Py_d;
  Double I_ERI_Dyz_Px_F3y_Py_d = I_ERI_Fxyz_S_F3y_Py_d+ABX*I_ERI_Dyz_S_F3y_Py_d;
  Double I_ERI_D2z_Px_F3y_Py_d = I_ERI_Fx2z_S_F3y_Py_d+ABX*I_ERI_D2z_S_F3y_Py_d;
  Double I_ERI_D2x_Py_F3y_Py_d = I_ERI_F2xy_S_F3y_Py_d+ABY*I_ERI_D2x_S_F3y_Py_d;
  Double I_ERI_Dxy_Py_F3y_Py_d = I_ERI_Fx2y_S_F3y_Py_d+ABY*I_ERI_Dxy_S_F3y_Py_d;
  Double I_ERI_Dxz_Py_F3y_Py_d = I_ERI_Fxyz_S_F3y_Py_d+ABY*I_ERI_Dxz_S_F3y_Py_d;
  Double I_ERI_D2y_Py_F3y_Py_d = I_ERI_F3y_S_F3y_Py_d+ABY*I_ERI_D2y_S_F3y_Py_d;
  Double I_ERI_Dyz_Py_F3y_Py_d = I_ERI_F2yz_S_F3y_Py_d+ABY*I_ERI_Dyz_S_F3y_Py_d;
  Double I_ERI_D2z_Py_F3y_Py_d = I_ERI_Fy2z_S_F3y_Py_d+ABY*I_ERI_D2z_S_F3y_Py_d;
  Double I_ERI_D2x_Pz_F3y_Py_d = I_ERI_F2xz_S_F3y_Py_d+ABZ*I_ERI_D2x_S_F3y_Py_d;
  Double I_ERI_Dxy_Pz_F3y_Py_d = I_ERI_Fxyz_S_F3y_Py_d+ABZ*I_ERI_Dxy_S_F3y_Py_d;
  Double I_ERI_Dxz_Pz_F3y_Py_d = I_ERI_Fx2z_S_F3y_Py_d+ABZ*I_ERI_Dxz_S_F3y_Py_d;
  Double I_ERI_D2y_Pz_F3y_Py_d = I_ERI_F2yz_S_F3y_Py_d+ABZ*I_ERI_D2y_S_F3y_Py_d;
  Double I_ERI_Dyz_Pz_F3y_Py_d = I_ERI_Fy2z_S_F3y_Py_d+ABZ*I_ERI_Dyz_S_F3y_Py_d;
  Double I_ERI_D2z_Pz_F3y_Py_d = I_ERI_F3z_S_F3y_Py_d+ABZ*I_ERI_D2z_S_F3y_Py_d;
  Double I_ERI_D2x_Px_F2yz_Py_d = I_ERI_F3x_S_F2yz_Py_d+ABX*I_ERI_D2x_S_F2yz_Py_d;
  Double I_ERI_Dxy_Px_F2yz_Py_d = I_ERI_F2xy_S_F2yz_Py_d+ABX*I_ERI_Dxy_S_F2yz_Py_d;
  Double I_ERI_Dxz_Px_F2yz_Py_d = I_ERI_F2xz_S_F2yz_Py_d+ABX*I_ERI_Dxz_S_F2yz_Py_d;
  Double I_ERI_D2y_Px_F2yz_Py_d = I_ERI_Fx2y_S_F2yz_Py_d+ABX*I_ERI_D2y_S_F2yz_Py_d;
  Double I_ERI_Dyz_Px_F2yz_Py_d = I_ERI_Fxyz_S_F2yz_Py_d+ABX*I_ERI_Dyz_S_F2yz_Py_d;
  Double I_ERI_D2z_Px_F2yz_Py_d = I_ERI_Fx2z_S_F2yz_Py_d+ABX*I_ERI_D2z_S_F2yz_Py_d;
  Double I_ERI_D2x_Py_F2yz_Py_d = I_ERI_F2xy_S_F2yz_Py_d+ABY*I_ERI_D2x_S_F2yz_Py_d;
  Double I_ERI_Dxy_Py_F2yz_Py_d = I_ERI_Fx2y_S_F2yz_Py_d+ABY*I_ERI_Dxy_S_F2yz_Py_d;
  Double I_ERI_Dxz_Py_F2yz_Py_d = I_ERI_Fxyz_S_F2yz_Py_d+ABY*I_ERI_Dxz_S_F2yz_Py_d;
  Double I_ERI_D2y_Py_F2yz_Py_d = I_ERI_F3y_S_F2yz_Py_d+ABY*I_ERI_D2y_S_F2yz_Py_d;
  Double I_ERI_Dyz_Py_F2yz_Py_d = I_ERI_F2yz_S_F2yz_Py_d+ABY*I_ERI_Dyz_S_F2yz_Py_d;
  Double I_ERI_D2z_Py_F2yz_Py_d = I_ERI_Fy2z_S_F2yz_Py_d+ABY*I_ERI_D2z_S_F2yz_Py_d;
  Double I_ERI_D2x_Pz_F2yz_Py_d = I_ERI_F2xz_S_F2yz_Py_d+ABZ*I_ERI_D2x_S_F2yz_Py_d;
  Double I_ERI_Dxy_Pz_F2yz_Py_d = I_ERI_Fxyz_S_F2yz_Py_d+ABZ*I_ERI_Dxy_S_F2yz_Py_d;
  Double I_ERI_Dxz_Pz_F2yz_Py_d = I_ERI_Fx2z_S_F2yz_Py_d+ABZ*I_ERI_Dxz_S_F2yz_Py_d;
  Double I_ERI_D2y_Pz_F2yz_Py_d = I_ERI_F2yz_S_F2yz_Py_d+ABZ*I_ERI_D2y_S_F2yz_Py_d;
  Double I_ERI_Dyz_Pz_F2yz_Py_d = I_ERI_Fy2z_S_F2yz_Py_d+ABZ*I_ERI_Dyz_S_F2yz_Py_d;
  Double I_ERI_D2z_Pz_F2yz_Py_d = I_ERI_F3z_S_F2yz_Py_d+ABZ*I_ERI_D2z_S_F2yz_Py_d;
  Double I_ERI_D2x_Px_Fy2z_Py_d = I_ERI_F3x_S_Fy2z_Py_d+ABX*I_ERI_D2x_S_Fy2z_Py_d;
  Double I_ERI_Dxy_Px_Fy2z_Py_d = I_ERI_F2xy_S_Fy2z_Py_d+ABX*I_ERI_Dxy_S_Fy2z_Py_d;
  Double I_ERI_Dxz_Px_Fy2z_Py_d = I_ERI_F2xz_S_Fy2z_Py_d+ABX*I_ERI_Dxz_S_Fy2z_Py_d;
  Double I_ERI_D2y_Px_Fy2z_Py_d = I_ERI_Fx2y_S_Fy2z_Py_d+ABX*I_ERI_D2y_S_Fy2z_Py_d;
  Double I_ERI_Dyz_Px_Fy2z_Py_d = I_ERI_Fxyz_S_Fy2z_Py_d+ABX*I_ERI_Dyz_S_Fy2z_Py_d;
  Double I_ERI_D2z_Px_Fy2z_Py_d = I_ERI_Fx2z_S_Fy2z_Py_d+ABX*I_ERI_D2z_S_Fy2z_Py_d;
  Double I_ERI_D2x_Py_Fy2z_Py_d = I_ERI_F2xy_S_Fy2z_Py_d+ABY*I_ERI_D2x_S_Fy2z_Py_d;
  Double I_ERI_Dxy_Py_Fy2z_Py_d = I_ERI_Fx2y_S_Fy2z_Py_d+ABY*I_ERI_Dxy_S_Fy2z_Py_d;
  Double I_ERI_Dxz_Py_Fy2z_Py_d = I_ERI_Fxyz_S_Fy2z_Py_d+ABY*I_ERI_Dxz_S_Fy2z_Py_d;
  Double I_ERI_D2y_Py_Fy2z_Py_d = I_ERI_F3y_S_Fy2z_Py_d+ABY*I_ERI_D2y_S_Fy2z_Py_d;
  Double I_ERI_Dyz_Py_Fy2z_Py_d = I_ERI_F2yz_S_Fy2z_Py_d+ABY*I_ERI_Dyz_S_Fy2z_Py_d;
  Double I_ERI_D2z_Py_Fy2z_Py_d = I_ERI_Fy2z_S_Fy2z_Py_d+ABY*I_ERI_D2z_S_Fy2z_Py_d;
  Double I_ERI_D2x_Pz_Fy2z_Py_d = I_ERI_F2xz_S_Fy2z_Py_d+ABZ*I_ERI_D2x_S_Fy2z_Py_d;
  Double I_ERI_Dxy_Pz_Fy2z_Py_d = I_ERI_Fxyz_S_Fy2z_Py_d+ABZ*I_ERI_Dxy_S_Fy2z_Py_d;
  Double I_ERI_Dxz_Pz_Fy2z_Py_d = I_ERI_Fx2z_S_Fy2z_Py_d+ABZ*I_ERI_Dxz_S_Fy2z_Py_d;
  Double I_ERI_D2y_Pz_Fy2z_Py_d = I_ERI_F2yz_S_Fy2z_Py_d+ABZ*I_ERI_D2y_S_Fy2z_Py_d;
  Double I_ERI_Dyz_Pz_Fy2z_Py_d = I_ERI_Fy2z_S_Fy2z_Py_d+ABZ*I_ERI_Dyz_S_Fy2z_Py_d;
  Double I_ERI_D2z_Pz_Fy2z_Py_d = I_ERI_F3z_S_Fy2z_Py_d+ABZ*I_ERI_D2z_S_Fy2z_Py_d;
  Double I_ERI_D2x_Px_F3z_Py_d = I_ERI_F3x_S_F3z_Py_d+ABX*I_ERI_D2x_S_F3z_Py_d;
  Double I_ERI_Dxy_Px_F3z_Py_d = I_ERI_F2xy_S_F3z_Py_d+ABX*I_ERI_Dxy_S_F3z_Py_d;
  Double I_ERI_Dxz_Px_F3z_Py_d = I_ERI_F2xz_S_F3z_Py_d+ABX*I_ERI_Dxz_S_F3z_Py_d;
  Double I_ERI_D2y_Px_F3z_Py_d = I_ERI_Fx2y_S_F3z_Py_d+ABX*I_ERI_D2y_S_F3z_Py_d;
  Double I_ERI_Dyz_Px_F3z_Py_d = I_ERI_Fxyz_S_F3z_Py_d+ABX*I_ERI_Dyz_S_F3z_Py_d;
  Double I_ERI_D2z_Px_F3z_Py_d = I_ERI_Fx2z_S_F3z_Py_d+ABX*I_ERI_D2z_S_F3z_Py_d;
  Double I_ERI_D2x_Py_F3z_Py_d = I_ERI_F2xy_S_F3z_Py_d+ABY*I_ERI_D2x_S_F3z_Py_d;
  Double I_ERI_Dxy_Py_F3z_Py_d = I_ERI_Fx2y_S_F3z_Py_d+ABY*I_ERI_Dxy_S_F3z_Py_d;
  Double I_ERI_Dxz_Py_F3z_Py_d = I_ERI_Fxyz_S_F3z_Py_d+ABY*I_ERI_Dxz_S_F3z_Py_d;
  Double I_ERI_D2y_Py_F3z_Py_d = I_ERI_F3y_S_F3z_Py_d+ABY*I_ERI_D2y_S_F3z_Py_d;
  Double I_ERI_Dyz_Py_F3z_Py_d = I_ERI_F2yz_S_F3z_Py_d+ABY*I_ERI_Dyz_S_F3z_Py_d;
  Double I_ERI_D2z_Py_F3z_Py_d = I_ERI_Fy2z_S_F3z_Py_d+ABY*I_ERI_D2z_S_F3z_Py_d;
  Double I_ERI_D2x_Pz_F3z_Py_d = I_ERI_F2xz_S_F3z_Py_d+ABZ*I_ERI_D2x_S_F3z_Py_d;
  Double I_ERI_Dxy_Pz_F3z_Py_d = I_ERI_Fxyz_S_F3z_Py_d+ABZ*I_ERI_Dxy_S_F3z_Py_d;
  Double I_ERI_Dxz_Pz_F3z_Py_d = I_ERI_Fx2z_S_F3z_Py_d+ABZ*I_ERI_Dxz_S_F3z_Py_d;
  Double I_ERI_D2y_Pz_F3z_Py_d = I_ERI_F2yz_S_F3z_Py_d+ABZ*I_ERI_D2y_S_F3z_Py_d;
  Double I_ERI_Dyz_Pz_F3z_Py_d = I_ERI_Fy2z_S_F3z_Py_d+ABZ*I_ERI_Dyz_S_F3z_Py_d;
  Double I_ERI_D2z_Pz_F3z_Py_d = I_ERI_F3z_S_F3z_Py_d+ABZ*I_ERI_D2z_S_F3z_Py_d;
  Double I_ERI_D2x_Px_F3x_Pz_d = I_ERI_F3x_S_F3x_Pz_d+ABX*I_ERI_D2x_S_F3x_Pz_d;
  Double I_ERI_Dxy_Px_F3x_Pz_d = I_ERI_F2xy_S_F3x_Pz_d+ABX*I_ERI_Dxy_S_F3x_Pz_d;
  Double I_ERI_Dxz_Px_F3x_Pz_d = I_ERI_F2xz_S_F3x_Pz_d+ABX*I_ERI_Dxz_S_F3x_Pz_d;
  Double I_ERI_D2y_Px_F3x_Pz_d = I_ERI_Fx2y_S_F3x_Pz_d+ABX*I_ERI_D2y_S_F3x_Pz_d;
  Double I_ERI_Dyz_Px_F3x_Pz_d = I_ERI_Fxyz_S_F3x_Pz_d+ABX*I_ERI_Dyz_S_F3x_Pz_d;
  Double I_ERI_D2z_Px_F3x_Pz_d = I_ERI_Fx2z_S_F3x_Pz_d+ABX*I_ERI_D2z_S_F3x_Pz_d;
  Double I_ERI_D2x_Py_F3x_Pz_d = I_ERI_F2xy_S_F3x_Pz_d+ABY*I_ERI_D2x_S_F3x_Pz_d;
  Double I_ERI_Dxy_Py_F3x_Pz_d = I_ERI_Fx2y_S_F3x_Pz_d+ABY*I_ERI_Dxy_S_F3x_Pz_d;
  Double I_ERI_Dxz_Py_F3x_Pz_d = I_ERI_Fxyz_S_F3x_Pz_d+ABY*I_ERI_Dxz_S_F3x_Pz_d;
  Double I_ERI_D2y_Py_F3x_Pz_d = I_ERI_F3y_S_F3x_Pz_d+ABY*I_ERI_D2y_S_F3x_Pz_d;
  Double I_ERI_Dyz_Py_F3x_Pz_d = I_ERI_F2yz_S_F3x_Pz_d+ABY*I_ERI_Dyz_S_F3x_Pz_d;
  Double I_ERI_D2z_Py_F3x_Pz_d = I_ERI_Fy2z_S_F3x_Pz_d+ABY*I_ERI_D2z_S_F3x_Pz_d;
  Double I_ERI_D2x_Pz_F3x_Pz_d = I_ERI_F2xz_S_F3x_Pz_d+ABZ*I_ERI_D2x_S_F3x_Pz_d;
  Double I_ERI_Dxy_Pz_F3x_Pz_d = I_ERI_Fxyz_S_F3x_Pz_d+ABZ*I_ERI_Dxy_S_F3x_Pz_d;
  Double I_ERI_Dxz_Pz_F3x_Pz_d = I_ERI_Fx2z_S_F3x_Pz_d+ABZ*I_ERI_Dxz_S_F3x_Pz_d;
  Double I_ERI_D2y_Pz_F3x_Pz_d = I_ERI_F2yz_S_F3x_Pz_d+ABZ*I_ERI_D2y_S_F3x_Pz_d;
  Double I_ERI_Dyz_Pz_F3x_Pz_d = I_ERI_Fy2z_S_F3x_Pz_d+ABZ*I_ERI_Dyz_S_F3x_Pz_d;
  Double I_ERI_D2z_Pz_F3x_Pz_d = I_ERI_F3z_S_F3x_Pz_d+ABZ*I_ERI_D2z_S_F3x_Pz_d;
  Double I_ERI_D2x_Px_F2xy_Pz_d = I_ERI_F3x_S_F2xy_Pz_d+ABX*I_ERI_D2x_S_F2xy_Pz_d;
  Double I_ERI_Dxy_Px_F2xy_Pz_d = I_ERI_F2xy_S_F2xy_Pz_d+ABX*I_ERI_Dxy_S_F2xy_Pz_d;
  Double I_ERI_Dxz_Px_F2xy_Pz_d = I_ERI_F2xz_S_F2xy_Pz_d+ABX*I_ERI_Dxz_S_F2xy_Pz_d;
  Double I_ERI_D2y_Px_F2xy_Pz_d = I_ERI_Fx2y_S_F2xy_Pz_d+ABX*I_ERI_D2y_S_F2xy_Pz_d;
  Double I_ERI_Dyz_Px_F2xy_Pz_d = I_ERI_Fxyz_S_F2xy_Pz_d+ABX*I_ERI_Dyz_S_F2xy_Pz_d;
  Double I_ERI_D2z_Px_F2xy_Pz_d = I_ERI_Fx2z_S_F2xy_Pz_d+ABX*I_ERI_D2z_S_F2xy_Pz_d;
  Double I_ERI_D2x_Py_F2xy_Pz_d = I_ERI_F2xy_S_F2xy_Pz_d+ABY*I_ERI_D2x_S_F2xy_Pz_d;
  Double I_ERI_Dxy_Py_F2xy_Pz_d = I_ERI_Fx2y_S_F2xy_Pz_d+ABY*I_ERI_Dxy_S_F2xy_Pz_d;
  Double I_ERI_Dxz_Py_F2xy_Pz_d = I_ERI_Fxyz_S_F2xy_Pz_d+ABY*I_ERI_Dxz_S_F2xy_Pz_d;
  Double I_ERI_D2y_Py_F2xy_Pz_d = I_ERI_F3y_S_F2xy_Pz_d+ABY*I_ERI_D2y_S_F2xy_Pz_d;
  Double I_ERI_Dyz_Py_F2xy_Pz_d = I_ERI_F2yz_S_F2xy_Pz_d+ABY*I_ERI_Dyz_S_F2xy_Pz_d;
  Double I_ERI_D2z_Py_F2xy_Pz_d = I_ERI_Fy2z_S_F2xy_Pz_d+ABY*I_ERI_D2z_S_F2xy_Pz_d;
  Double I_ERI_D2x_Pz_F2xy_Pz_d = I_ERI_F2xz_S_F2xy_Pz_d+ABZ*I_ERI_D2x_S_F2xy_Pz_d;
  Double I_ERI_Dxy_Pz_F2xy_Pz_d = I_ERI_Fxyz_S_F2xy_Pz_d+ABZ*I_ERI_Dxy_S_F2xy_Pz_d;
  Double I_ERI_Dxz_Pz_F2xy_Pz_d = I_ERI_Fx2z_S_F2xy_Pz_d+ABZ*I_ERI_Dxz_S_F2xy_Pz_d;
  Double I_ERI_D2y_Pz_F2xy_Pz_d = I_ERI_F2yz_S_F2xy_Pz_d+ABZ*I_ERI_D2y_S_F2xy_Pz_d;
  Double I_ERI_Dyz_Pz_F2xy_Pz_d = I_ERI_Fy2z_S_F2xy_Pz_d+ABZ*I_ERI_Dyz_S_F2xy_Pz_d;
  Double I_ERI_D2z_Pz_F2xy_Pz_d = I_ERI_F3z_S_F2xy_Pz_d+ABZ*I_ERI_D2z_S_F2xy_Pz_d;
  Double I_ERI_D2x_Px_F2xz_Pz_d = I_ERI_F3x_S_F2xz_Pz_d+ABX*I_ERI_D2x_S_F2xz_Pz_d;
  Double I_ERI_Dxy_Px_F2xz_Pz_d = I_ERI_F2xy_S_F2xz_Pz_d+ABX*I_ERI_Dxy_S_F2xz_Pz_d;
  Double I_ERI_Dxz_Px_F2xz_Pz_d = I_ERI_F2xz_S_F2xz_Pz_d+ABX*I_ERI_Dxz_S_F2xz_Pz_d;
  Double I_ERI_D2y_Px_F2xz_Pz_d = I_ERI_Fx2y_S_F2xz_Pz_d+ABX*I_ERI_D2y_S_F2xz_Pz_d;
  Double I_ERI_Dyz_Px_F2xz_Pz_d = I_ERI_Fxyz_S_F2xz_Pz_d+ABX*I_ERI_Dyz_S_F2xz_Pz_d;
  Double I_ERI_D2z_Px_F2xz_Pz_d = I_ERI_Fx2z_S_F2xz_Pz_d+ABX*I_ERI_D2z_S_F2xz_Pz_d;
  Double I_ERI_D2x_Py_F2xz_Pz_d = I_ERI_F2xy_S_F2xz_Pz_d+ABY*I_ERI_D2x_S_F2xz_Pz_d;
  Double I_ERI_Dxy_Py_F2xz_Pz_d = I_ERI_Fx2y_S_F2xz_Pz_d+ABY*I_ERI_Dxy_S_F2xz_Pz_d;
  Double I_ERI_Dxz_Py_F2xz_Pz_d = I_ERI_Fxyz_S_F2xz_Pz_d+ABY*I_ERI_Dxz_S_F2xz_Pz_d;
  Double I_ERI_D2y_Py_F2xz_Pz_d = I_ERI_F3y_S_F2xz_Pz_d+ABY*I_ERI_D2y_S_F2xz_Pz_d;
  Double I_ERI_Dyz_Py_F2xz_Pz_d = I_ERI_F2yz_S_F2xz_Pz_d+ABY*I_ERI_Dyz_S_F2xz_Pz_d;
  Double I_ERI_D2z_Py_F2xz_Pz_d = I_ERI_Fy2z_S_F2xz_Pz_d+ABY*I_ERI_D2z_S_F2xz_Pz_d;
  Double I_ERI_D2x_Pz_F2xz_Pz_d = I_ERI_F2xz_S_F2xz_Pz_d+ABZ*I_ERI_D2x_S_F2xz_Pz_d;
  Double I_ERI_Dxy_Pz_F2xz_Pz_d = I_ERI_Fxyz_S_F2xz_Pz_d+ABZ*I_ERI_Dxy_S_F2xz_Pz_d;
  Double I_ERI_Dxz_Pz_F2xz_Pz_d = I_ERI_Fx2z_S_F2xz_Pz_d+ABZ*I_ERI_Dxz_S_F2xz_Pz_d;
  Double I_ERI_D2y_Pz_F2xz_Pz_d = I_ERI_F2yz_S_F2xz_Pz_d+ABZ*I_ERI_D2y_S_F2xz_Pz_d;
  Double I_ERI_Dyz_Pz_F2xz_Pz_d = I_ERI_Fy2z_S_F2xz_Pz_d+ABZ*I_ERI_Dyz_S_F2xz_Pz_d;
  Double I_ERI_D2z_Pz_F2xz_Pz_d = I_ERI_F3z_S_F2xz_Pz_d+ABZ*I_ERI_D2z_S_F2xz_Pz_d;
  Double I_ERI_D2x_Px_Fx2y_Pz_d = I_ERI_F3x_S_Fx2y_Pz_d+ABX*I_ERI_D2x_S_Fx2y_Pz_d;
  Double I_ERI_Dxy_Px_Fx2y_Pz_d = I_ERI_F2xy_S_Fx2y_Pz_d+ABX*I_ERI_Dxy_S_Fx2y_Pz_d;
  Double I_ERI_Dxz_Px_Fx2y_Pz_d = I_ERI_F2xz_S_Fx2y_Pz_d+ABX*I_ERI_Dxz_S_Fx2y_Pz_d;
  Double I_ERI_D2y_Px_Fx2y_Pz_d = I_ERI_Fx2y_S_Fx2y_Pz_d+ABX*I_ERI_D2y_S_Fx2y_Pz_d;
  Double I_ERI_Dyz_Px_Fx2y_Pz_d = I_ERI_Fxyz_S_Fx2y_Pz_d+ABX*I_ERI_Dyz_S_Fx2y_Pz_d;
  Double I_ERI_D2z_Px_Fx2y_Pz_d = I_ERI_Fx2z_S_Fx2y_Pz_d+ABX*I_ERI_D2z_S_Fx2y_Pz_d;
  Double I_ERI_D2x_Py_Fx2y_Pz_d = I_ERI_F2xy_S_Fx2y_Pz_d+ABY*I_ERI_D2x_S_Fx2y_Pz_d;
  Double I_ERI_Dxy_Py_Fx2y_Pz_d = I_ERI_Fx2y_S_Fx2y_Pz_d+ABY*I_ERI_Dxy_S_Fx2y_Pz_d;
  Double I_ERI_Dxz_Py_Fx2y_Pz_d = I_ERI_Fxyz_S_Fx2y_Pz_d+ABY*I_ERI_Dxz_S_Fx2y_Pz_d;
  Double I_ERI_D2y_Py_Fx2y_Pz_d = I_ERI_F3y_S_Fx2y_Pz_d+ABY*I_ERI_D2y_S_Fx2y_Pz_d;
  Double I_ERI_Dyz_Py_Fx2y_Pz_d = I_ERI_F2yz_S_Fx2y_Pz_d+ABY*I_ERI_Dyz_S_Fx2y_Pz_d;
  Double I_ERI_D2z_Py_Fx2y_Pz_d = I_ERI_Fy2z_S_Fx2y_Pz_d+ABY*I_ERI_D2z_S_Fx2y_Pz_d;
  Double I_ERI_D2x_Pz_Fx2y_Pz_d = I_ERI_F2xz_S_Fx2y_Pz_d+ABZ*I_ERI_D2x_S_Fx2y_Pz_d;
  Double I_ERI_Dxy_Pz_Fx2y_Pz_d = I_ERI_Fxyz_S_Fx2y_Pz_d+ABZ*I_ERI_Dxy_S_Fx2y_Pz_d;
  Double I_ERI_Dxz_Pz_Fx2y_Pz_d = I_ERI_Fx2z_S_Fx2y_Pz_d+ABZ*I_ERI_Dxz_S_Fx2y_Pz_d;
  Double I_ERI_D2y_Pz_Fx2y_Pz_d = I_ERI_F2yz_S_Fx2y_Pz_d+ABZ*I_ERI_D2y_S_Fx2y_Pz_d;
  Double I_ERI_Dyz_Pz_Fx2y_Pz_d = I_ERI_Fy2z_S_Fx2y_Pz_d+ABZ*I_ERI_Dyz_S_Fx2y_Pz_d;
  Double I_ERI_D2z_Pz_Fx2y_Pz_d = I_ERI_F3z_S_Fx2y_Pz_d+ABZ*I_ERI_D2z_S_Fx2y_Pz_d;
  Double I_ERI_D2x_Px_Fxyz_Pz_d = I_ERI_F3x_S_Fxyz_Pz_d+ABX*I_ERI_D2x_S_Fxyz_Pz_d;
  Double I_ERI_Dxy_Px_Fxyz_Pz_d = I_ERI_F2xy_S_Fxyz_Pz_d+ABX*I_ERI_Dxy_S_Fxyz_Pz_d;
  Double I_ERI_Dxz_Px_Fxyz_Pz_d = I_ERI_F2xz_S_Fxyz_Pz_d+ABX*I_ERI_Dxz_S_Fxyz_Pz_d;
  Double I_ERI_D2y_Px_Fxyz_Pz_d = I_ERI_Fx2y_S_Fxyz_Pz_d+ABX*I_ERI_D2y_S_Fxyz_Pz_d;
  Double I_ERI_Dyz_Px_Fxyz_Pz_d = I_ERI_Fxyz_S_Fxyz_Pz_d+ABX*I_ERI_Dyz_S_Fxyz_Pz_d;
  Double I_ERI_D2z_Px_Fxyz_Pz_d = I_ERI_Fx2z_S_Fxyz_Pz_d+ABX*I_ERI_D2z_S_Fxyz_Pz_d;
  Double I_ERI_D2x_Py_Fxyz_Pz_d = I_ERI_F2xy_S_Fxyz_Pz_d+ABY*I_ERI_D2x_S_Fxyz_Pz_d;
  Double I_ERI_Dxy_Py_Fxyz_Pz_d = I_ERI_Fx2y_S_Fxyz_Pz_d+ABY*I_ERI_Dxy_S_Fxyz_Pz_d;
  Double I_ERI_Dxz_Py_Fxyz_Pz_d = I_ERI_Fxyz_S_Fxyz_Pz_d+ABY*I_ERI_Dxz_S_Fxyz_Pz_d;
  Double I_ERI_D2y_Py_Fxyz_Pz_d = I_ERI_F3y_S_Fxyz_Pz_d+ABY*I_ERI_D2y_S_Fxyz_Pz_d;
  Double I_ERI_Dyz_Py_Fxyz_Pz_d = I_ERI_F2yz_S_Fxyz_Pz_d+ABY*I_ERI_Dyz_S_Fxyz_Pz_d;
  Double I_ERI_D2z_Py_Fxyz_Pz_d = I_ERI_Fy2z_S_Fxyz_Pz_d+ABY*I_ERI_D2z_S_Fxyz_Pz_d;
  Double I_ERI_D2x_Pz_Fxyz_Pz_d = I_ERI_F2xz_S_Fxyz_Pz_d+ABZ*I_ERI_D2x_S_Fxyz_Pz_d;
  Double I_ERI_Dxy_Pz_Fxyz_Pz_d = I_ERI_Fxyz_S_Fxyz_Pz_d+ABZ*I_ERI_Dxy_S_Fxyz_Pz_d;
  Double I_ERI_Dxz_Pz_Fxyz_Pz_d = I_ERI_Fx2z_S_Fxyz_Pz_d+ABZ*I_ERI_Dxz_S_Fxyz_Pz_d;
  Double I_ERI_D2y_Pz_Fxyz_Pz_d = I_ERI_F2yz_S_Fxyz_Pz_d+ABZ*I_ERI_D2y_S_Fxyz_Pz_d;
  Double I_ERI_Dyz_Pz_Fxyz_Pz_d = I_ERI_Fy2z_S_Fxyz_Pz_d+ABZ*I_ERI_Dyz_S_Fxyz_Pz_d;
  Double I_ERI_D2z_Pz_Fxyz_Pz_d = I_ERI_F3z_S_Fxyz_Pz_d+ABZ*I_ERI_D2z_S_Fxyz_Pz_d;
  Double I_ERI_D2x_Px_Fx2z_Pz_d = I_ERI_F3x_S_Fx2z_Pz_d+ABX*I_ERI_D2x_S_Fx2z_Pz_d;
  Double I_ERI_Dxy_Px_Fx2z_Pz_d = I_ERI_F2xy_S_Fx2z_Pz_d+ABX*I_ERI_Dxy_S_Fx2z_Pz_d;
  Double I_ERI_Dxz_Px_Fx2z_Pz_d = I_ERI_F2xz_S_Fx2z_Pz_d+ABX*I_ERI_Dxz_S_Fx2z_Pz_d;
  Double I_ERI_D2y_Px_Fx2z_Pz_d = I_ERI_Fx2y_S_Fx2z_Pz_d+ABX*I_ERI_D2y_S_Fx2z_Pz_d;
  Double I_ERI_Dyz_Px_Fx2z_Pz_d = I_ERI_Fxyz_S_Fx2z_Pz_d+ABX*I_ERI_Dyz_S_Fx2z_Pz_d;
  Double I_ERI_D2z_Px_Fx2z_Pz_d = I_ERI_Fx2z_S_Fx2z_Pz_d+ABX*I_ERI_D2z_S_Fx2z_Pz_d;
  Double I_ERI_D2x_Py_Fx2z_Pz_d = I_ERI_F2xy_S_Fx2z_Pz_d+ABY*I_ERI_D2x_S_Fx2z_Pz_d;
  Double I_ERI_Dxy_Py_Fx2z_Pz_d = I_ERI_Fx2y_S_Fx2z_Pz_d+ABY*I_ERI_Dxy_S_Fx2z_Pz_d;
  Double I_ERI_Dxz_Py_Fx2z_Pz_d = I_ERI_Fxyz_S_Fx2z_Pz_d+ABY*I_ERI_Dxz_S_Fx2z_Pz_d;
  Double I_ERI_D2y_Py_Fx2z_Pz_d = I_ERI_F3y_S_Fx2z_Pz_d+ABY*I_ERI_D2y_S_Fx2z_Pz_d;
  Double I_ERI_Dyz_Py_Fx2z_Pz_d = I_ERI_F2yz_S_Fx2z_Pz_d+ABY*I_ERI_Dyz_S_Fx2z_Pz_d;
  Double I_ERI_D2z_Py_Fx2z_Pz_d = I_ERI_Fy2z_S_Fx2z_Pz_d+ABY*I_ERI_D2z_S_Fx2z_Pz_d;
  Double I_ERI_D2x_Pz_Fx2z_Pz_d = I_ERI_F2xz_S_Fx2z_Pz_d+ABZ*I_ERI_D2x_S_Fx2z_Pz_d;
  Double I_ERI_Dxy_Pz_Fx2z_Pz_d = I_ERI_Fxyz_S_Fx2z_Pz_d+ABZ*I_ERI_Dxy_S_Fx2z_Pz_d;
  Double I_ERI_Dxz_Pz_Fx2z_Pz_d = I_ERI_Fx2z_S_Fx2z_Pz_d+ABZ*I_ERI_Dxz_S_Fx2z_Pz_d;
  Double I_ERI_D2y_Pz_Fx2z_Pz_d = I_ERI_F2yz_S_Fx2z_Pz_d+ABZ*I_ERI_D2y_S_Fx2z_Pz_d;
  Double I_ERI_Dyz_Pz_Fx2z_Pz_d = I_ERI_Fy2z_S_Fx2z_Pz_d+ABZ*I_ERI_Dyz_S_Fx2z_Pz_d;
  Double I_ERI_D2z_Pz_Fx2z_Pz_d = I_ERI_F3z_S_Fx2z_Pz_d+ABZ*I_ERI_D2z_S_Fx2z_Pz_d;
  Double I_ERI_D2x_Px_F3y_Pz_d = I_ERI_F3x_S_F3y_Pz_d+ABX*I_ERI_D2x_S_F3y_Pz_d;
  Double I_ERI_Dxy_Px_F3y_Pz_d = I_ERI_F2xy_S_F3y_Pz_d+ABX*I_ERI_Dxy_S_F3y_Pz_d;
  Double I_ERI_Dxz_Px_F3y_Pz_d = I_ERI_F2xz_S_F3y_Pz_d+ABX*I_ERI_Dxz_S_F3y_Pz_d;
  Double I_ERI_D2y_Px_F3y_Pz_d = I_ERI_Fx2y_S_F3y_Pz_d+ABX*I_ERI_D2y_S_F3y_Pz_d;
  Double I_ERI_Dyz_Px_F3y_Pz_d = I_ERI_Fxyz_S_F3y_Pz_d+ABX*I_ERI_Dyz_S_F3y_Pz_d;
  Double I_ERI_D2z_Px_F3y_Pz_d = I_ERI_Fx2z_S_F3y_Pz_d+ABX*I_ERI_D2z_S_F3y_Pz_d;
  Double I_ERI_D2x_Py_F3y_Pz_d = I_ERI_F2xy_S_F3y_Pz_d+ABY*I_ERI_D2x_S_F3y_Pz_d;
  Double I_ERI_Dxy_Py_F3y_Pz_d = I_ERI_Fx2y_S_F3y_Pz_d+ABY*I_ERI_Dxy_S_F3y_Pz_d;
  Double I_ERI_Dxz_Py_F3y_Pz_d = I_ERI_Fxyz_S_F3y_Pz_d+ABY*I_ERI_Dxz_S_F3y_Pz_d;
  Double I_ERI_D2y_Py_F3y_Pz_d = I_ERI_F3y_S_F3y_Pz_d+ABY*I_ERI_D2y_S_F3y_Pz_d;
  Double I_ERI_Dyz_Py_F3y_Pz_d = I_ERI_F2yz_S_F3y_Pz_d+ABY*I_ERI_Dyz_S_F3y_Pz_d;
  Double I_ERI_D2z_Py_F3y_Pz_d = I_ERI_Fy2z_S_F3y_Pz_d+ABY*I_ERI_D2z_S_F3y_Pz_d;
  Double I_ERI_D2x_Pz_F3y_Pz_d = I_ERI_F2xz_S_F3y_Pz_d+ABZ*I_ERI_D2x_S_F3y_Pz_d;
  Double I_ERI_Dxy_Pz_F3y_Pz_d = I_ERI_Fxyz_S_F3y_Pz_d+ABZ*I_ERI_Dxy_S_F3y_Pz_d;
  Double I_ERI_Dxz_Pz_F3y_Pz_d = I_ERI_Fx2z_S_F3y_Pz_d+ABZ*I_ERI_Dxz_S_F3y_Pz_d;
  Double I_ERI_D2y_Pz_F3y_Pz_d = I_ERI_F2yz_S_F3y_Pz_d+ABZ*I_ERI_D2y_S_F3y_Pz_d;
  Double I_ERI_Dyz_Pz_F3y_Pz_d = I_ERI_Fy2z_S_F3y_Pz_d+ABZ*I_ERI_Dyz_S_F3y_Pz_d;
  Double I_ERI_D2z_Pz_F3y_Pz_d = I_ERI_F3z_S_F3y_Pz_d+ABZ*I_ERI_D2z_S_F3y_Pz_d;
  Double I_ERI_D2x_Px_F2yz_Pz_d = I_ERI_F3x_S_F2yz_Pz_d+ABX*I_ERI_D2x_S_F2yz_Pz_d;
  Double I_ERI_Dxy_Px_F2yz_Pz_d = I_ERI_F2xy_S_F2yz_Pz_d+ABX*I_ERI_Dxy_S_F2yz_Pz_d;
  Double I_ERI_Dxz_Px_F2yz_Pz_d = I_ERI_F2xz_S_F2yz_Pz_d+ABX*I_ERI_Dxz_S_F2yz_Pz_d;
  Double I_ERI_D2y_Px_F2yz_Pz_d = I_ERI_Fx2y_S_F2yz_Pz_d+ABX*I_ERI_D2y_S_F2yz_Pz_d;
  Double I_ERI_Dyz_Px_F2yz_Pz_d = I_ERI_Fxyz_S_F2yz_Pz_d+ABX*I_ERI_Dyz_S_F2yz_Pz_d;
  Double I_ERI_D2z_Px_F2yz_Pz_d = I_ERI_Fx2z_S_F2yz_Pz_d+ABX*I_ERI_D2z_S_F2yz_Pz_d;
  Double I_ERI_D2x_Py_F2yz_Pz_d = I_ERI_F2xy_S_F2yz_Pz_d+ABY*I_ERI_D2x_S_F2yz_Pz_d;
  Double I_ERI_Dxy_Py_F2yz_Pz_d = I_ERI_Fx2y_S_F2yz_Pz_d+ABY*I_ERI_Dxy_S_F2yz_Pz_d;
  Double I_ERI_Dxz_Py_F2yz_Pz_d = I_ERI_Fxyz_S_F2yz_Pz_d+ABY*I_ERI_Dxz_S_F2yz_Pz_d;
  Double I_ERI_D2y_Py_F2yz_Pz_d = I_ERI_F3y_S_F2yz_Pz_d+ABY*I_ERI_D2y_S_F2yz_Pz_d;
  Double I_ERI_Dyz_Py_F2yz_Pz_d = I_ERI_F2yz_S_F2yz_Pz_d+ABY*I_ERI_Dyz_S_F2yz_Pz_d;
  Double I_ERI_D2z_Py_F2yz_Pz_d = I_ERI_Fy2z_S_F2yz_Pz_d+ABY*I_ERI_D2z_S_F2yz_Pz_d;
  Double I_ERI_D2x_Pz_F2yz_Pz_d = I_ERI_F2xz_S_F2yz_Pz_d+ABZ*I_ERI_D2x_S_F2yz_Pz_d;
  Double I_ERI_Dxy_Pz_F2yz_Pz_d = I_ERI_Fxyz_S_F2yz_Pz_d+ABZ*I_ERI_Dxy_S_F2yz_Pz_d;
  Double I_ERI_Dxz_Pz_F2yz_Pz_d = I_ERI_Fx2z_S_F2yz_Pz_d+ABZ*I_ERI_Dxz_S_F2yz_Pz_d;
  Double I_ERI_D2y_Pz_F2yz_Pz_d = I_ERI_F2yz_S_F2yz_Pz_d+ABZ*I_ERI_D2y_S_F2yz_Pz_d;
  Double I_ERI_Dyz_Pz_F2yz_Pz_d = I_ERI_Fy2z_S_F2yz_Pz_d+ABZ*I_ERI_Dyz_S_F2yz_Pz_d;
  Double I_ERI_D2z_Pz_F2yz_Pz_d = I_ERI_F3z_S_F2yz_Pz_d+ABZ*I_ERI_D2z_S_F2yz_Pz_d;
  Double I_ERI_D2x_Px_Fy2z_Pz_d = I_ERI_F3x_S_Fy2z_Pz_d+ABX*I_ERI_D2x_S_Fy2z_Pz_d;
  Double I_ERI_Dxy_Px_Fy2z_Pz_d = I_ERI_F2xy_S_Fy2z_Pz_d+ABX*I_ERI_Dxy_S_Fy2z_Pz_d;
  Double I_ERI_Dxz_Px_Fy2z_Pz_d = I_ERI_F2xz_S_Fy2z_Pz_d+ABX*I_ERI_Dxz_S_Fy2z_Pz_d;
  Double I_ERI_D2y_Px_Fy2z_Pz_d = I_ERI_Fx2y_S_Fy2z_Pz_d+ABX*I_ERI_D2y_S_Fy2z_Pz_d;
  Double I_ERI_Dyz_Px_Fy2z_Pz_d = I_ERI_Fxyz_S_Fy2z_Pz_d+ABX*I_ERI_Dyz_S_Fy2z_Pz_d;
  Double I_ERI_D2z_Px_Fy2z_Pz_d = I_ERI_Fx2z_S_Fy2z_Pz_d+ABX*I_ERI_D2z_S_Fy2z_Pz_d;
  Double I_ERI_D2x_Py_Fy2z_Pz_d = I_ERI_F2xy_S_Fy2z_Pz_d+ABY*I_ERI_D2x_S_Fy2z_Pz_d;
  Double I_ERI_Dxy_Py_Fy2z_Pz_d = I_ERI_Fx2y_S_Fy2z_Pz_d+ABY*I_ERI_Dxy_S_Fy2z_Pz_d;
  Double I_ERI_Dxz_Py_Fy2z_Pz_d = I_ERI_Fxyz_S_Fy2z_Pz_d+ABY*I_ERI_Dxz_S_Fy2z_Pz_d;
  Double I_ERI_D2y_Py_Fy2z_Pz_d = I_ERI_F3y_S_Fy2z_Pz_d+ABY*I_ERI_D2y_S_Fy2z_Pz_d;
  Double I_ERI_Dyz_Py_Fy2z_Pz_d = I_ERI_F2yz_S_Fy2z_Pz_d+ABY*I_ERI_Dyz_S_Fy2z_Pz_d;
  Double I_ERI_D2z_Py_Fy2z_Pz_d = I_ERI_Fy2z_S_Fy2z_Pz_d+ABY*I_ERI_D2z_S_Fy2z_Pz_d;
  Double I_ERI_D2x_Pz_Fy2z_Pz_d = I_ERI_F2xz_S_Fy2z_Pz_d+ABZ*I_ERI_D2x_S_Fy2z_Pz_d;
  Double I_ERI_Dxy_Pz_Fy2z_Pz_d = I_ERI_Fxyz_S_Fy2z_Pz_d+ABZ*I_ERI_Dxy_S_Fy2z_Pz_d;
  Double I_ERI_Dxz_Pz_Fy2z_Pz_d = I_ERI_Fx2z_S_Fy2z_Pz_d+ABZ*I_ERI_Dxz_S_Fy2z_Pz_d;
  Double I_ERI_D2y_Pz_Fy2z_Pz_d = I_ERI_F2yz_S_Fy2z_Pz_d+ABZ*I_ERI_D2y_S_Fy2z_Pz_d;
  Double I_ERI_Dyz_Pz_Fy2z_Pz_d = I_ERI_Fy2z_S_Fy2z_Pz_d+ABZ*I_ERI_Dyz_S_Fy2z_Pz_d;
  Double I_ERI_D2z_Pz_Fy2z_Pz_d = I_ERI_F3z_S_Fy2z_Pz_d+ABZ*I_ERI_D2z_S_Fy2z_Pz_d;
  Double I_ERI_D2x_Px_F3z_Pz_d = I_ERI_F3x_S_F3z_Pz_d+ABX*I_ERI_D2x_S_F3z_Pz_d;
  Double I_ERI_Dxy_Px_F3z_Pz_d = I_ERI_F2xy_S_F3z_Pz_d+ABX*I_ERI_Dxy_S_F3z_Pz_d;
  Double I_ERI_Dxz_Px_F3z_Pz_d = I_ERI_F2xz_S_F3z_Pz_d+ABX*I_ERI_Dxz_S_F3z_Pz_d;
  Double I_ERI_D2y_Px_F3z_Pz_d = I_ERI_Fx2y_S_F3z_Pz_d+ABX*I_ERI_D2y_S_F3z_Pz_d;
  Double I_ERI_Dyz_Px_F3z_Pz_d = I_ERI_Fxyz_S_F3z_Pz_d+ABX*I_ERI_Dyz_S_F3z_Pz_d;
  Double I_ERI_D2z_Px_F3z_Pz_d = I_ERI_Fx2z_S_F3z_Pz_d+ABX*I_ERI_D2z_S_F3z_Pz_d;
  Double I_ERI_D2x_Py_F3z_Pz_d = I_ERI_F2xy_S_F3z_Pz_d+ABY*I_ERI_D2x_S_F3z_Pz_d;
  Double I_ERI_Dxy_Py_F3z_Pz_d = I_ERI_Fx2y_S_F3z_Pz_d+ABY*I_ERI_Dxy_S_F3z_Pz_d;
  Double I_ERI_Dxz_Py_F3z_Pz_d = I_ERI_Fxyz_S_F3z_Pz_d+ABY*I_ERI_Dxz_S_F3z_Pz_d;
  Double I_ERI_D2y_Py_F3z_Pz_d = I_ERI_F3y_S_F3z_Pz_d+ABY*I_ERI_D2y_S_F3z_Pz_d;
  Double I_ERI_Dyz_Py_F3z_Pz_d = I_ERI_F2yz_S_F3z_Pz_d+ABY*I_ERI_Dyz_S_F3z_Pz_d;
  Double I_ERI_D2z_Py_F3z_Pz_d = I_ERI_Fy2z_S_F3z_Pz_d+ABY*I_ERI_D2z_S_F3z_Pz_d;
  Double I_ERI_D2x_Pz_F3z_Pz_d = I_ERI_F2xz_S_F3z_Pz_d+ABZ*I_ERI_D2x_S_F3z_Pz_d;
  Double I_ERI_Dxy_Pz_F3z_Pz_d = I_ERI_Fxyz_S_F3z_Pz_d+ABZ*I_ERI_Dxy_S_F3z_Pz_d;
  Double I_ERI_Dxz_Pz_F3z_Pz_d = I_ERI_Fx2z_S_F3z_Pz_d+ABZ*I_ERI_Dxz_S_F3z_Pz_d;
  Double I_ERI_D2y_Pz_F3z_Pz_d = I_ERI_F2yz_S_F3z_Pz_d+ABZ*I_ERI_D2y_S_F3z_Pz_d;
  Double I_ERI_Dyz_Pz_F3z_Pz_d = I_ERI_Fy2z_S_F3z_Pz_d+ABZ*I_ERI_Dyz_S_F3z_Pz_d;
  Double I_ERI_D2z_Pz_F3z_Pz_d = I_ERI_F3z_S_F3z_Pz_d+ABZ*I_ERI_D2z_S_F3z_Pz_d;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_a
   * RHS shell quartet name: SQ_ERI_P_P_F_S
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_F3x_Px_F3x_S_a-2*I_ERI_Px_Px_F3x_S;
  abcd[1] = 2.0E0*I_ERI_F2xy_Px_F3x_S_a-1*I_ERI_Py_Px_F3x_S;
  abcd[2] = 2.0E0*I_ERI_F2xz_Px_F3x_S_a-1*I_ERI_Pz_Px_F3x_S;
  abcd[3] = 2.0E0*I_ERI_Fx2y_Px_F3x_S_a;
  abcd[4] = 2.0E0*I_ERI_Fxyz_Px_F3x_S_a;
  abcd[5] = 2.0E0*I_ERI_Fx2z_Px_F3x_S_a;
  abcd[6] = 2.0E0*I_ERI_F3x_Py_F3x_S_a-2*I_ERI_Px_Py_F3x_S;
  abcd[7] = 2.0E0*I_ERI_F2xy_Py_F3x_S_a-1*I_ERI_Py_Py_F3x_S;
  abcd[8] = 2.0E0*I_ERI_F2xz_Py_F3x_S_a-1*I_ERI_Pz_Py_F3x_S;
  abcd[9] = 2.0E0*I_ERI_Fx2y_Py_F3x_S_a;
  abcd[10] = 2.0E0*I_ERI_Fxyz_Py_F3x_S_a;
  abcd[11] = 2.0E0*I_ERI_Fx2z_Py_F3x_S_a;
  abcd[12] = 2.0E0*I_ERI_F3x_Pz_F3x_S_a-2*I_ERI_Px_Pz_F3x_S;
  abcd[13] = 2.0E0*I_ERI_F2xy_Pz_F3x_S_a-1*I_ERI_Py_Pz_F3x_S;
  abcd[14] = 2.0E0*I_ERI_F2xz_Pz_F3x_S_a-1*I_ERI_Pz_Pz_F3x_S;
  abcd[15] = 2.0E0*I_ERI_Fx2y_Pz_F3x_S_a;
  abcd[16] = 2.0E0*I_ERI_Fxyz_Pz_F3x_S_a;
  abcd[17] = 2.0E0*I_ERI_Fx2z_Pz_F3x_S_a;
  abcd[18] = 2.0E0*I_ERI_F3x_Px_F2xy_S_a-2*I_ERI_Px_Px_F2xy_S;
  abcd[19] = 2.0E0*I_ERI_F2xy_Px_F2xy_S_a-1*I_ERI_Py_Px_F2xy_S;
  abcd[20] = 2.0E0*I_ERI_F2xz_Px_F2xy_S_a-1*I_ERI_Pz_Px_F2xy_S;
  abcd[21] = 2.0E0*I_ERI_Fx2y_Px_F2xy_S_a;
  abcd[22] = 2.0E0*I_ERI_Fxyz_Px_F2xy_S_a;
  abcd[23] = 2.0E0*I_ERI_Fx2z_Px_F2xy_S_a;
  abcd[24] = 2.0E0*I_ERI_F3x_Py_F2xy_S_a-2*I_ERI_Px_Py_F2xy_S;
  abcd[25] = 2.0E0*I_ERI_F2xy_Py_F2xy_S_a-1*I_ERI_Py_Py_F2xy_S;
  abcd[26] = 2.0E0*I_ERI_F2xz_Py_F2xy_S_a-1*I_ERI_Pz_Py_F2xy_S;
  abcd[27] = 2.0E0*I_ERI_Fx2y_Py_F2xy_S_a;
  abcd[28] = 2.0E0*I_ERI_Fxyz_Py_F2xy_S_a;
  abcd[29] = 2.0E0*I_ERI_Fx2z_Py_F2xy_S_a;
  abcd[30] = 2.0E0*I_ERI_F3x_Pz_F2xy_S_a-2*I_ERI_Px_Pz_F2xy_S;
  abcd[31] = 2.0E0*I_ERI_F2xy_Pz_F2xy_S_a-1*I_ERI_Py_Pz_F2xy_S;
  abcd[32] = 2.0E0*I_ERI_F2xz_Pz_F2xy_S_a-1*I_ERI_Pz_Pz_F2xy_S;
  abcd[33] = 2.0E0*I_ERI_Fx2y_Pz_F2xy_S_a;
  abcd[34] = 2.0E0*I_ERI_Fxyz_Pz_F2xy_S_a;
  abcd[35] = 2.0E0*I_ERI_Fx2z_Pz_F2xy_S_a;
  abcd[36] = 2.0E0*I_ERI_F3x_Px_F2xz_S_a-2*I_ERI_Px_Px_F2xz_S;
  abcd[37] = 2.0E0*I_ERI_F2xy_Px_F2xz_S_a-1*I_ERI_Py_Px_F2xz_S;
  abcd[38] = 2.0E0*I_ERI_F2xz_Px_F2xz_S_a-1*I_ERI_Pz_Px_F2xz_S;
  abcd[39] = 2.0E0*I_ERI_Fx2y_Px_F2xz_S_a;
  abcd[40] = 2.0E0*I_ERI_Fxyz_Px_F2xz_S_a;
  abcd[41] = 2.0E0*I_ERI_Fx2z_Px_F2xz_S_a;
  abcd[42] = 2.0E0*I_ERI_F3x_Py_F2xz_S_a-2*I_ERI_Px_Py_F2xz_S;
  abcd[43] = 2.0E0*I_ERI_F2xy_Py_F2xz_S_a-1*I_ERI_Py_Py_F2xz_S;
  abcd[44] = 2.0E0*I_ERI_F2xz_Py_F2xz_S_a-1*I_ERI_Pz_Py_F2xz_S;
  abcd[45] = 2.0E0*I_ERI_Fx2y_Py_F2xz_S_a;
  abcd[46] = 2.0E0*I_ERI_Fxyz_Py_F2xz_S_a;
  abcd[47] = 2.0E0*I_ERI_Fx2z_Py_F2xz_S_a;
  abcd[48] = 2.0E0*I_ERI_F3x_Pz_F2xz_S_a-2*I_ERI_Px_Pz_F2xz_S;
  abcd[49] = 2.0E0*I_ERI_F2xy_Pz_F2xz_S_a-1*I_ERI_Py_Pz_F2xz_S;
  abcd[50] = 2.0E0*I_ERI_F2xz_Pz_F2xz_S_a-1*I_ERI_Pz_Pz_F2xz_S;
  abcd[51] = 2.0E0*I_ERI_Fx2y_Pz_F2xz_S_a;
  abcd[52] = 2.0E0*I_ERI_Fxyz_Pz_F2xz_S_a;
  abcd[53] = 2.0E0*I_ERI_Fx2z_Pz_F2xz_S_a;
  abcd[54] = 2.0E0*I_ERI_F3x_Px_Fx2y_S_a-2*I_ERI_Px_Px_Fx2y_S;
  abcd[55] = 2.0E0*I_ERI_F2xy_Px_Fx2y_S_a-1*I_ERI_Py_Px_Fx2y_S;
  abcd[56] = 2.0E0*I_ERI_F2xz_Px_Fx2y_S_a-1*I_ERI_Pz_Px_Fx2y_S;
  abcd[57] = 2.0E0*I_ERI_Fx2y_Px_Fx2y_S_a;
  abcd[58] = 2.0E0*I_ERI_Fxyz_Px_Fx2y_S_a;
  abcd[59] = 2.0E0*I_ERI_Fx2z_Px_Fx2y_S_a;
  abcd[60] = 2.0E0*I_ERI_F3x_Py_Fx2y_S_a-2*I_ERI_Px_Py_Fx2y_S;
  abcd[61] = 2.0E0*I_ERI_F2xy_Py_Fx2y_S_a-1*I_ERI_Py_Py_Fx2y_S;
  abcd[62] = 2.0E0*I_ERI_F2xz_Py_Fx2y_S_a-1*I_ERI_Pz_Py_Fx2y_S;
  abcd[63] = 2.0E0*I_ERI_Fx2y_Py_Fx2y_S_a;
  abcd[64] = 2.0E0*I_ERI_Fxyz_Py_Fx2y_S_a;
  abcd[65] = 2.0E0*I_ERI_Fx2z_Py_Fx2y_S_a;
  abcd[66] = 2.0E0*I_ERI_F3x_Pz_Fx2y_S_a-2*I_ERI_Px_Pz_Fx2y_S;
  abcd[67] = 2.0E0*I_ERI_F2xy_Pz_Fx2y_S_a-1*I_ERI_Py_Pz_Fx2y_S;
  abcd[68] = 2.0E0*I_ERI_F2xz_Pz_Fx2y_S_a-1*I_ERI_Pz_Pz_Fx2y_S;
  abcd[69] = 2.0E0*I_ERI_Fx2y_Pz_Fx2y_S_a;
  abcd[70] = 2.0E0*I_ERI_Fxyz_Pz_Fx2y_S_a;
  abcd[71] = 2.0E0*I_ERI_Fx2z_Pz_Fx2y_S_a;
  abcd[72] = 2.0E0*I_ERI_F3x_Px_Fxyz_S_a-2*I_ERI_Px_Px_Fxyz_S;
  abcd[73] = 2.0E0*I_ERI_F2xy_Px_Fxyz_S_a-1*I_ERI_Py_Px_Fxyz_S;
  abcd[74] = 2.0E0*I_ERI_F2xz_Px_Fxyz_S_a-1*I_ERI_Pz_Px_Fxyz_S;
  abcd[75] = 2.0E0*I_ERI_Fx2y_Px_Fxyz_S_a;
  abcd[76] = 2.0E0*I_ERI_Fxyz_Px_Fxyz_S_a;
  abcd[77] = 2.0E0*I_ERI_Fx2z_Px_Fxyz_S_a;
  abcd[78] = 2.0E0*I_ERI_F3x_Py_Fxyz_S_a-2*I_ERI_Px_Py_Fxyz_S;
  abcd[79] = 2.0E0*I_ERI_F2xy_Py_Fxyz_S_a-1*I_ERI_Py_Py_Fxyz_S;
  abcd[80] = 2.0E0*I_ERI_F2xz_Py_Fxyz_S_a-1*I_ERI_Pz_Py_Fxyz_S;
  abcd[81] = 2.0E0*I_ERI_Fx2y_Py_Fxyz_S_a;
  abcd[82] = 2.0E0*I_ERI_Fxyz_Py_Fxyz_S_a;
  abcd[83] = 2.0E0*I_ERI_Fx2z_Py_Fxyz_S_a;
  abcd[84] = 2.0E0*I_ERI_F3x_Pz_Fxyz_S_a-2*I_ERI_Px_Pz_Fxyz_S;
  abcd[85] = 2.0E0*I_ERI_F2xy_Pz_Fxyz_S_a-1*I_ERI_Py_Pz_Fxyz_S;
  abcd[86] = 2.0E0*I_ERI_F2xz_Pz_Fxyz_S_a-1*I_ERI_Pz_Pz_Fxyz_S;
  abcd[87] = 2.0E0*I_ERI_Fx2y_Pz_Fxyz_S_a;
  abcd[88] = 2.0E0*I_ERI_Fxyz_Pz_Fxyz_S_a;
  abcd[89] = 2.0E0*I_ERI_Fx2z_Pz_Fxyz_S_a;
  abcd[90] = 2.0E0*I_ERI_F3x_Px_Fx2z_S_a-2*I_ERI_Px_Px_Fx2z_S;
  abcd[91] = 2.0E0*I_ERI_F2xy_Px_Fx2z_S_a-1*I_ERI_Py_Px_Fx2z_S;
  abcd[92] = 2.0E0*I_ERI_F2xz_Px_Fx2z_S_a-1*I_ERI_Pz_Px_Fx2z_S;
  abcd[93] = 2.0E0*I_ERI_Fx2y_Px_Fx2z_S_a;
  abcd[94] = 2.0E0*I_ERI_Fxyz_Px_Fx2z_S_a;
  abcd[95] = 2.0E0*I_ERI_Fx2z_Px_Fx2z_S_a;
  abcd[96] = 2.0E0*I_ERI_F3x_Py_Fx2z_S_a-2*I_ERI_Px_Py_Fx2z_S;
  abcd[97] = 2.0E0*I_ERI_F2xy_Py_Fx2z_S_a-1*I_ERI_Py_Py_Fx2z_S;
  abcd[98] = 2.0E0*I_ERI_F2xz_Py_Fx2z_S_a-1*I_ERI_Pz_Py_Fx2z_S;
  abcd[99] = 2.0E0*I_ERI_Fx2y_Py_Fx2z_S_a;
  abcd[100] = 2.0E0*I_ERI_Fxyz_Py_Fx2z_S_a;
  abcd[101] = 2.0E0*I_ERI_Fx2z_Py_Fx2z_S_a;
  abcd[102] = 2.0E0*I_ERI_F3x_Pz_Fx2z_S_a-2*I_ERI_Px_Pz_Fx2z_S;
  abcd[103] = 2.0E0*I_ERI_F2xy_Pz_Fx2z_S_a-1*I_ERI_Py_Pz_Fx2z_S;
  abcd[104] = 2.0E0*I_ERI_F2xz_Pz_Fx2z_S_a-1*I_ERI_Pz_Pz_Fx2z_S;
  abcd[105] = 2.0E0*I_ERI_Fx2y_Pz_Fx2z_S_a;
  abcd[106] = 2.0E0*I_ERI_Fxyz_Pz_Fx2z_S_a;
  abcd[107] = 2.0E0*I_ERI_Fx2z_Pz_Fx2z_S_a;
  abcd[108] = 2.0E0*I_ERI_F3x_Px_F3y_S_a-2*I_ERI_Px_Px_F3y_S;
  abcd[109] = 2.0E0*I_ERI_F2xy_Px_F3y_S_a-1*I_ERI_Py_Px_F3y_S;
  abcd[110] = 2.0E0*I_ERI_F2xz_Px_F3y_S_a-1*I_ERI_Pz_Px_F3y_S;
  abcd[111] = 2.0E0*I_ERI_Fx2y_Px_F3y_S_a;
  abcd[112] = 2.0E0*I_ERI_Fxyz_Px_F3y_S_a;
  abcd[113] = 2.0E0*I_ERI_Fx2z_Px_F3y_S_a;
  abcd[114] = 2.0E0*I_ERI_F3x_Py_F3y_S_a-2*I_ERI_Px_Py_F3y_S;
  abcd[115] = 2.0E0*I_ERI_F2xy_Py_F3y_S_a-1*I_ERI_Py_Py_F3y_S;
  abcd[116] = 2.0E0*I_ERI_F2xz_Py_F3y_S_a-1*I_ERI_Pz_Py_F3y_S;
  abcd[117] = 2.0E0*I_ERI_Fx2y_Py_F3y_S_a;
  abcd[118] = 2.0E0*I_ERI_Fxyz_Py_F3y_S_a;
  abcd[119] = 2.0E0*I_ERI_Fx2z_Py_F3y_S_a;
  abcd[120] = 2.0E0*I_ERI_F3x_Pz_F3y_S_a-2*I_ERI_Px_Pz_F3y_S;
  abcd[121] = 2.0E0*I_ERI_F2xy_Pz_F3y_S_a-1*I_ERI_Py_Pz_F3y_S;
  abcd[122] = 2.0E0*I_ERI_F2xz_Pz_F3y_S_a-1*I_ERI_Pz_Pz_F3y_S;
  abcd[123] = 2.0E0*I_ERI_Fx2y_Pz_F3y_S_a;
  abcd[124] = 2.0E0*I_ERI_Fxyz_Pz_F3y_S_a;
  abcd[125] = 2.0E0*I_ERI_Fx2z_Pz_F3y_S_a;
  abcd[126] = 2.0E0*I_ERI_F3x_Px_F2yz_S_a-2*I_ERI_Px_Px_F2yz_S;
  abcd[127] = 2.0E0*I_ERI_F2xy_Px_F2yz_S_a-1*I_ERI_Py_Px_F2yz_S;
  abcd[128] = 2.0E0*I_ERI_F2xz_Px_F2yz_S_a-1*I_ERI_Pz_Px_F2yz_S;
  abcd[129] = 2.0E0*I_ERI_Fx2y_Px_F2yz_S_a;
  abcd[130] = 2.0E0*I_ERI_Fxyz_Px_F2yz_S_a;
  abcd[131] = 2.0E0*I_ERI_Fx2z_Px_F2yz_S_a;
  abcd[132] = 2.0E0*I_ERI_F3x_Py_F2yz_S_a-2*I_ERI_Px_Py_F2yz_S;
  abcd[133] = 2.0E0*I_ERI_F2xy_Py_F2yz_S_a-1*I_ERI_Py_Py_F2yz_S;
  abcd[134] = 2.0E0*I_ERI_F2xz_Py_F2yz_S_a-1*I_ERI_Pz_Py_F2yz_S;
  abcd[135] = 2.0E0*I_ERI_Fx2y_Py_F2yz_S_a;
  abcd[136] = 2.0E0*I_ERI_Fxyz_Py_F2yz_S_a;
  abcd[137] = 2.0E0*I_ERI_Fx2z_Py_F2yz_S_a;
  abcd[138] = 2.0E0*I_ERI_F3x_Pz_F2yz_S_a-2*I_ERI_Px_Pz_F2yz_S;
  abcd[139] = 2.0E0*I_ERI_F2xy_Pz_F2yz_S_a-1*I_ERI_Py_Pz_F2yz_S;
  abcd[140] = 2.0E0*I_ERI_F2xz_Pz_F2yz_S_a-1*I_ERI_Pz_Pz_F2yz_S;
  abcd[141] = 2.0E0*I_ERI_Fx2y_Pz_F2yz_S_a;
  abcd[142] = 2.0E0*I_ERI_Fxyz_Pz_F2yz_S_a;
  abcd[143] = 2.0E0*I_ERI_Fx2z_Pz_F2yz_S_a;
  abcd[144] = 2.0E0*I_ERI_F3x_Px_Fy2z_S_a-2*I_ERI_Px_Px_Fy2z_S;
  abcd[145] = 2.0E0*I_ERI_F2xy_Px_Fy2z_S_a-1*I_ERI_Py_Px_Fy2z_S;
  abcd[146] = 2.0E0*I_ERI_F2xz_Px_Fy2z_S_a-1*I_ERI_Pz_Px_Fy2z_S;
  abcd[147] = 2.0E0*I_ERI_Fx2y_Px_Fy2z_S_a;
  abcd[148] = 2.0E0*I_ERI_Fxyz_Px_Fy2z_S_a;
  abcd[149] = 2.0E0*I_ERI_Fx2z_Px_Fy2z_S_a;
  abcd[150] = 2.0E0*I_ERI_F3x_Py_Fy2z_S_a-2*I_ERI_Px_Py_Fy2z_S;
  abcd[151] = 2.0E0*I_ERI_F2xy_Py_Fy2z_S_a-1*I_ERI_Py_Py_Fy2z_S;
  abcd[152] = 2.0E0*I_ERI_F2xz_Py_Fy2z_S_a-1*I_ERI_Pz_Py_Fy2z_S;
  abcd[153] = 2.0E0*I_ERI_Fx2y_Py_Fy2z_S_a;
  abcd[154] = 2.0E0*I_ERI_Fxyz_Py_Fy2z_S_a;
  abcd[155] = 2.0E0*I_ERI_Fx2z_Py_Fy2z_S_a;
  abcd[156] = 2.0E0*I_ERI_F3x_Pz_Fy2z_S_a-2*I_ERI_Px_Pz_Fy2z_S;
  abcd[157] = 2.0E0*I_ERI_F2xy_Pz_Fy2z_S_a-1*I_ERI_Py_Pz_Fy2z_S;
  abcd[158] = 2.0E0*I_ERI_F2xz_Pz_Fy2z_S_a-1*I_ERI_Pz_Pz_Fy2z_S;
  abcd[159] = 2.0E0*I_ERI_Fx2y_Pz_Fy2z_S_a;
  abcd[160] = 2.0E0*I_ERI_Fxyz_Pz_Fy2z_S_a;
  abcd[161] = 2.0E0*I_ERI_Fx2z_Pz_Fy2z_S_a;
  abcd[162] = 2.0E0*I_ERI_F3x_Px_F3z_S_a-2*I_ERI_Px_Px_F3z_S;
  abcd[163] = 2.0E0*I_ERI_F2xy_Px_F3z_S_a-1*I_ERI_Py_Px_F3z_S;
  abcd[164] = 2.0E0*I_ERI_F2xz_Px_F3z_S_a-1*I_ERI_Pz_Px_F3z_S;
  abcd[165] = 2.0E0*I_ERI_Fx2y_Px_F3z_S_a;
  abcd[166] = 2.0E0*I_ERI_Fxyz_Px_F3z_S_a;
  abcd[167] = 2.0E0*I_ERI_Fx2z_Px_F3z_S_a;
  abcd[168] = 2.0E0*I_ERI_F3x_Py_F3z_S_a-2*I_ERI_Px_Py_F3z_S;
  abcd[169] = 2.0E0*I_ERI_F2xy_Py_F3z_S_a-1*I_ERI_Py_Py_F3z_S;
  abcd[170] = 2.0E0*I_ERI_F2xz_Py_F3z_S_a-1*I_ERI_Pz_Py_F3z_S;
  abcd[171] = 2.0E0*I_ERI_Fx2y_Py_F3z_S_a;
  abcd[172] = 2.0E0*I_ERI_Fxyz_Py_F3z_S_a;
  abcd[173] = 2.0E0*I_ERI_Fx2z_Py_F3z_S_a;
  abcd[174] = 2.0E0*I_ERI_F3x_Pz_F3z_S_a-2*I_ERI_Px_Pz_F3z_S;
  abcd[175] = 2.0E0*I_ERI_F2xy_Pz_F3z_S_a-1*I_ERI_Py_Pz_F3z_S;
  abcd[176] = 2.0E0*I_ERI_F2xz_Pz_F3z_S_a-1*I_ERI_Pz_Pz_F3z_S;
  abcd[177] = 2.0E0*I_ERI_Fx2y_Pz_F3z_S_a;
  abcd[178] = 2.0E0*I_ERI_Fxyz_Pz_F3z_S_a;
  abcd[179] = 2.0E0*I_ERI_Fx2z_Pz_F3z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_a
   * RHS shell quartet name: SQ_ERI_P_P_F_S
   ************************************************************/
  abcd[180] = 2.0E0*I_ERI_F2xy_Px_F3x_S_a;
  abcd[181] = 2.0E0*I_ERI_Fx2y_Px_F3x_S_a-1*I_ERI_Px_Px_F3x_S;
  abcd[182] = 2.0E0*I_ERI_Fxyz_Px_F3x_S_a;
  abcd[183] = 2.0E0*I_ERI_F3y_Px_F3x_S_a-2*I_ERI_Py_Px_F3x_S;
  abcd[184] = 2.0E0*I_ERI_F2yz_Px_F3x_S_a-1*I_ERI_Pz_Px_F3x_S;
  abcd[185] = 2.0E0*I_ERI_Fy2z_Px_F3x_S_a;
  abcd[186] = 2.0E0*I_ERI_F2xy_Py_F3x_S_a;
  abcd[187] = 2.0E0*I_ERI_Fx2y_Py_F3x_S_a-1*I_ERI_Px_Py_F3x_S;
  abcd[188] = 2.0E0*I_ERI_Fxyz_Py_F3x_S_a;
  abcd[189] = 2.0E0*I_ERI_F3y_Py_F3x_S_a-2*I_ERI_Py_Py_F3x_S;
  abcd[190] = 2.0E0*I_ERI_F2yz_Py_F3x_S_a-1*I_ERI_Pz_Py_F3x_S;
  abcd[191] = 2.0E0*I_ERI_Fy2z_Py_F3x_S_a;
  abcd[192] = 2.0E0*I_ERI_F2xy_Pz_F3x_S_a;
  abcd[193] = 2.0E0*I_ERI_Fx2y_Pz_F3x_S_a-1*I_ERI_Px_Pz_F3x_S;
  abcd[194] = 2.0E0*I_ERI_Fxyz_Pz_F3x_S_a;
  abcd[195] = 2.0E0*I_ERI_F3y_Pz_F3x_S_a-2*I_ERI_Py_Pz_F3x_S;
  abcd[196] = 2.0E0*I_ERI_F2yz_Pz_F3x_S_a-1*I_ERI_Pz_Pz_F3x_S;
  abcd[197] = 2.0E0*I_ERI_Fy2z_Pz_F3x_S_a;
  abcd[198] = 2.0E0*I_ERI_F2xy_Px_F2xy_S_a;
  abcd[199] = 2.0E0*I_ERI_Fx2y_Px_F2xy_S_a-1*I_ERI_Px_Px_F2xy_S;
  abcd[200] = 2.0E0*I_ERI_Fxyz_Px_F2xy_S_a;
  abcd[201] = 2.0E0*I_ERI_F3y_Px_F2xy_S_a-2*I_ERI_Py_Px_F2xy_S;
  abcd[202] = 2.0E0*I_ERI_F2yz_Px_F2xy_S_a-1*I_ERI_Pz_Px_F2xy_S;
  abcd[203] = 2.0E0*I_ERI_Fy2z_Px_F2xy_S_a;
  abcd[204] = 2.0E0*I_ERI_F2xy_Py_F2xy_S_a;
  abcd[205] = 2.0E0*I_ERI_Fx2y_Py_F2xy_S_a-1*I_ERI_Px_Py_F2xy_S;
  abcd[206] = 2.0E0*I_ERI_Fxyz_Py_F2xy_S_a;
  abcd[207] = 2.0E0*I_ERI_F3y_Py_F2xy_S_a-2*I_ERI_Py_Py_F2xy_S;
  abcd[208] = 2.0E0*I_ERI_F2yz_Py_F2xy_S_a-1*I_ERI_Pz_Py_F2xy_S;
  abcd[209] = 2.0E0*I_ERI_Fy2z_Py_F2xy_S_a;
  abcd[210] = 2.0E0*I_ERI_F2xy_Pz_F2xy_S_a;
  abcd[211] = 2.0E0*I_ERI_Fx2y_Pz_F2xy_S_a-1*I_ERI_Px_Pz_F2xy_S;
  abcd[212] = 2.0E0*I_ERI_Fxyz_Pz_F2xy_S_a;
  abcd[213] = 2.0E0*I_ERI_F3y_Pz_F2xy_S_a-2*I_ERI_Py_Pz_F2xy_S;
  abcd[214] = 2.0E0*I_ERI_F2yz_Pz_F2xy_S_a-1*I_ERI_Pz_Pz_F2xy_S;
  abcd[215] = 2.0E0*I_ERI_Fy2z_Pz_F2xy_S_a;
  abcd[216] = 2.0E0*I_ERI_F2xy_Px_F2xz_S_a;
  abcd[217] = 2.0E0*I_ERI_Fx2y_Px_F2xz_S_a-1*I_ERI_Px_Px_F2xz_S;
  abcd[218] = 2.0E0*I_ERI_Fxyz_Px_F2xz_S_a;
  abcd[219] = 2.0E0*I_ERI_F3y_Px_F2xz_S_a-2*I_ERI_Py_Px_F2xz_S;
  abcd[220] = 2.0E0*I_ERI_F2yz_Px_F2xz_S_a-1*I_ERI_Pz_Px_F2xz_S;
  abcd[221] = 2.0E0*I_ERI_Fy2z_Px_F2xz_S_a;
  abcd[222] = 2.0E0*I_ERI_F2xy_Py_F2xz_S_a;
  abcd[223] = 2.0E0*I_ERI_Fx2y_Py_F2xz_S_a-1*I_ERI_Px_Py_F2xz_S;
  abcd[224] = 2.0E0*I_ERI_Fxyz_Py_F2xz_S_a;
  abcd[225] = 2.0E0*I_ERI_F3y_Py_F2xz_S_a-2*I_ERI_Py_Py_F2xz_S;
  abcd[226] = 2.0E0*I_ERI_F2yz_Py_F2xz_S_a-1*I_ERI_Pz_Py_F2xz_S;
  abcd[227] = 2.0E0*I_ERI_Fy2z_Py_F2xz_S_a;
  abcd[228] = 2.0E0*I_ERI_F2xy_Pz_F2xz_S_a;
  abcd[229] = 2.0E0*I_ERI_Fx2y_Pz_F2xz_S_a-1*I_ERI_Px_Pz_F2xz_S;
  abcd[230] = 2.0E0*I_ERI_Fxyz_Pz_F2xz_S_a;
  abcd[231] = 2.0E0*I_ERI_F3y_Pz_F2xz_S_a-2*I_ERI_Py_Pz_F2xz_S;
  abcd[232] = 2.0E0*I_ERI_F2yz_Pz_F2xz_S_a-1*I_ERI_Pz_Pz_F2xz_S;
  abcd[233] = 2.0E0*I_ERI_Fy2z_Pz_F2xz_S_a;
  abcd[234] = 2.0E0*I_ERI_F2xy_Px_Fx2y_S_a;
  abcd[235] = 2.0E0*I_ERI_Fx2y_Px_Fx2y_S_a-1*I_ERI_Px_Px_Fx2y_S;
  abcd[236] = 2.0E0*I_ERI_Fxyz_Px_Fx2y_S_a;
  abcd[237] = 2.0E0*I_ERI_F3y_Px_Fx2y_S_a-2*I_ERI_Py_Px_Fx2y_S;
  abcd[238] = 2.0E0*I_ERI_F2yz_Px_Fx2y_S_a-1*I_ERI_Pz_Px_Fx2y_S;
  abcd[239] = 2.0E0*I_ERI_Fy2z_Px_Fx2y_S_a;
  abcd[240] = 2.0E0*I_ERI_F2xy_Py_Fx2y_S_a;
  abcd[241] = 2.0E0*I_ERI_Fx2y_Py_Fx2y_S_a-1*I_ERI_Px_Py_Fx2y_S;
  abcd[242] = 2.0E0*I_ERI_Fxyz_Py_Fx2y_S_a;
  abcd[243] = 2.0E0*I_ERI_F3y_Py_Fx2y_S_a-2*I_ERI_Py_Py_Fx2y_S;
  abcd[244] = 2.0E0*I_ERI_F2yz_Py_Fx2y_S_a-1*I_ERI_Pz_Py_Fx2y_S;
  abcd[245] = 2.0E0*I_ERI_Fy2z_Py_Fx2y_S_a;
  abcd[246] = 2.0E0*I_ERI_F2xy_Pz_Fx2y_S_a;
  abcd[247] = 2.0E0*I_ERI_Fx2y_Pz_Fx2y_S_a-1*I_ERI_Px_Pz_Fx2y_S;
  abcd[248] = 2.0E0*I_ERI_Fxyz_Pz_Fx2y_S_a;
  abcd[249] = 2.0E0*I_ERI_F3y_Pz_Fx2y_S_a-2*I_ERI_Py_Pz_Fx2y_S;
  abcd[250] = 2.0E0*I_ERI_F2yz_Pz_Fx2y_S_a-1*I_ERI_Pz_Pz_Fx2y_S;
  abcd[251] = 2.0E0*I_ERI_Fy2z_Pz_Fx2y_S_a;
  abcd[252] = 2.0E0*I_ERI_F2xy_Px_Fxyz_S_a;
  abcd[253] = 2.0E0*I_ERI_Fx2y_Px_Fxyz_S_a-1*I_ERI_Px_Px_Fxyz_S;
  abcd[254] = 2.0E0*I_ERI_Fxyz_Px_Fxyz_S_a;
  abcd[255] = 2.0E0*I_ERI_F3y_Px_Fxyz_S_a-2*I_ERI_Py_Px_Fxyz_S;
  abcd[256] = 2.0E0*I_ERI_F2yz_Px_Fxyz_S_a-1*I_ERI_Pz_Px_Fxyz_S;
  abcd[257] = 2.0E0*I_ERI_Fy2z_Px_Fxyz_S_a;
  abcd[258] = 2.0E0*I_ERI_F2xy_Py_Fxyz_S_a;
  abcd[259] = 2.0E0*I_ERI_Fx2y_Py_Fxyz_S_a-1*I_ERI_Px_Py_Fxyz_S;
  abcd[260] = 2.0E0*I_ERI_Fxyz_Py_Fxyz_S_a;
  abcd[261] = 2.0E0*I_ERI_F3y_Py_Fxyz_S_a-2*I_ERI_Py_Py_Fxyz_S;
  abcd[262] = 2.0E0*I_ERI_F2yz_Py_Fxyz_S_a-1*I_ERI_Pz_Py_Fxyz_S;
  abcd[263] = 2.0E0*I_ERI_Fy2z_Py_Fxyz_S_a;
  abcd[264] = 2.0E0*I_ERI_F2xy_Pz_Fxyz_S_a;
  abcd[265] = 2.0E0*I_ERI_Fx2y_Pz_Fxyz_S_a-1*I_ERI_Px_Pz_Fxyz_S;
  abcd[266] = 2.0E0*I_ERI_Fxyz_Pz_Fxyz_S_a;
  abcd[267] = 2.0E0*I_ERI_F3y_Pz_Fxyz_S_a-2*I_ERI_Py_Pz_Fxyz_S;
  abcd[268] = 2.0E0*I_ERI_F2yz_Pz_Fxyz_S_a-1*I_ERI_Pz_Pz_Fxyz_S;
  abcd[269] = 2.0E0*I_ERI_Fy2z_Pz_Fxyz_S_a;
  abcd[270] = 2.0E0*I_ERI_F2xy_Px_Fx2z_S_a;
  abcd[271] = 2.0E0*I_ERI_Fx2y_Px_Fx2z_S_a-1*I_ERI_Px_Px_Fx2z_S;
  abcd[272] = 2.0E0*I_ERI_Fxyz_Px_Fx2z_S_a;
  abcd[273] = 2.0E0*I_ERI_F3y_Px_Fx2z_S_a-2*I_ERI_Py_Px_Fx2z_S;
  abcd[274] = 2.0E0*I_ERI_F2yz_Px_Fx2z_S_a-1*I_ERI_Pz_Px_Fx2z_S;
  abcd[275] = 2.0E0*I_ERI_Fy2z_Px_Fx2z_S_a;
  abcd[276] = 2.0E0*I_ERI_F2xy_Py_Fx2z_S_a;
  abcd[277] = 2.0E0*I_ERI_Fx2y_Py_Fx2z_S_a-1*I_ERI_Px_Py_Fx2z_S;
  abcd[278] = 2.0E0*I_ERI_Fxyz_Py_Fx2z_S_a;
  abcd[279] = 2.0E0*I_ERI_F3y_Py_Fx2z_S_a-2*I_ERI_Py_Py_Fx2z_S;
  abcd[280] = 2.0E0*I_ERI_F2yz_Py_Fx2z_S_a-1*I_ERI_Pz_Py_Fx2z_S;
  abcd[281] = 2.0E0*I_ERI_Fy2z_Py_Fx2z_S_a;
  abcd[282] = 2.0E0*I_ERI_F2xy_Pz_Fx2z_S_a;
  abcd[283] = 2.0E0*I_ERI_Fx2y_Pz_Fx2z_S_a-1*I_ERI_Px_Pz_Fx2z_S;
  abcd[284] = 2.0E0*I_ERI_Fxyz_Pz_Fx2z_S_a;
  abcd[285] = 2.0E0*I_ERI_F3y_Pz_Fx2z_S_a-2*I_ERI_Py_Pz_Fx2z_S;
  abcd[286] = 2.0E0*I_ERI_F2yz_Pz_Fx2z_S_a-1*I_ERI_Pz_Pz_Fx2z_S;
  abcd[287] = 2.0E0*I_ERI_Fy2z_Pz_Fx2z_S_a;
  abcd[288] = 2.0E0*I_ERI_F2xy_Px_F3y_S_a;
  abcd[289] = 2.0E0*I_ERI_Fx2y_Px_F3y_S_a-1*I_ERI_Px_Px_F3y_S;
  abcd[290] = 2.0E0*I_ERI_Fxyz_Px_F3y_S_a;
  abcd[291] = 2.0E0*I_ERI_F3y_Px_F3y_S_a-2*I_ERI_Py_Px_F3y_S;
  abcd[292] = 2.0E0*I_ERI_F2yz_Px_F3y_S_a-1*I_ERI_Pz_Px_F3y_S;
  abcd[293] = 2.0E0*I_ERI_Fy2z_Px_F3y_S_a;
  abcd[294] = 2.0E0*I_ERI_F2xy_Py_F3y_S_a;
  abcd[295] = 2.0E0*I_ERI_Fx2y_Py_F3y_S_a-1*I_ERI_Px_Py_F3y_S;
  abcd[296] = 2.0E0*I_ERI_Fxyz_Py_F3y_S_a;
  abcd[297] = 2.0E0*I_ERI_F3y_Py_F3y_S_a-2*I_ERI_Py_Py_F3y_S;
  abcd[298] = 2.0E0*I_ERI_F2yz_Py_F3y_S_a-1*I_ERI_Pz_Py_F3y_S;
  abcd[299] = 2.0E0*I_ERI_Fy2z_Py_F3y_S_a;
  abcd[300] = 2.0E0*I_ERI_F2xy_Pz_F3y_S_a;
  abcd[301] = 2.0E0*I_ERI_Fx2y_Pz_F3y_S_a-1*I_ERI_Px_Pz_F3y_S;
  abcd[302] = 2.0E0*I_ERI_Fxyz_Pz_F3y_S_a;
  abcd[303] = 2.0E0*I_ERI_F3y_Pz_F3y_S_a-2*I_ERI_Py_Pz_F3y_S;
  abcd[304] = 2.0E0*I_ERI_F2yz_Pz_F3y_S_a-1*I_ERI_Pz_Pz_F3y_S;
  abcd[305] = 2.0E0*I_ERI_Fy2z_Pz_F3y_S_a;
  abcd[306] = 2.0E0*I_ERI_F2xy_Px_F2yz_S_a;
  abcd[307] = 2.0E0*I_ERI_Fx2y_Px_F2yz_S_a-1*I_ERI_Px_Px_F2yz_S;
  abcd[308] = 2.0E0*I_ERI_Fxyz_Px_F2yz_S_a;
  abcd[309] = 2.0E0*I_ERI_F3y_Px_F2yz_S_a-2*I_ERI_Py_Px_F2yz_S;
  abcd[310] = 2.0E0*I_ERI_F2yz_Px_F2yz_S_a-1*I_ERI_Pz_Px_F2yz_S;
  abcd[311] = 2.0E0*I_ERI_Fy2z_Px_F2yz_S_a;
  abcd[312] = 2.0E0*I_ERI_F2xy_Py_F2yz_S_a;
  abcd[313] = 2.0E0*I_ERI_Fx2y_Py_F2yz_S_a-1*I_ERI_Px_Py_F2yz_S;
  abcd[314] = 2.0E0*I_ERI_Fxyz_Py_F2yz_S_a;
  abcd[315] = 2.0E0*I_ERI_F3y_Py_F2yz_S_a-2*I_ERI_Py_Py_F2yz_S;
  abcd[316] = 2.0E0*I_ERI_F2yz_Py_F2yz_S_a-1*I_ERI_Pz_Py_F2yz_S;
  abcd[317] = 2.0E0*I_ERI_Fy2z_Py_F2yz_S_a;
  abcd[318] = 2.0E0*I_ERI_F2xy_Pz_F2yz_S_a;
  abcd[319] = 2.0E0*I_ERI_Fx2y_Pz_F2yz_S_a-1*I_ERI_Px_Pz_F2yz_S;
  abcd[320] = 2.0E0*I_ERI_Fxyz_Pz_F2yz_S_a;
  abcd[321] = 2.0E0*I_ERI_F3y_Pz_F2yz_S_a-2*I_ERI_Py_Pz_F2yz_S;
  abcd[322] = 2.0E0*I_ERI_F2yz_Pz_F2yz_S_a-1*I_ERI_Pz_Pz_F2yz_S;
  abcd[323] = 2.0E0*I_ERI_Fy2z_Pz_F2yz_S_a;
  abcd[324] = 2.0E0*I_ERI_F2xy_Px_Fy2z_S_a;
  abcd[325] = 2.0E0*I_ERI_Fx2y_Px_Fy2z_S_a-1*I_ERI_Px_Px_Fy2z_S;
  abcd[326] = 2.0E0*I_ERI_Fxyz_Px_Fy2z_S_a;
  abcd[327] = 2.0E0*I_ERI_F3y_Px_Fy2z_S_a-2*I_ERI_Py_Px_Fy2z_S;
  abcd[328] = 2.0E0*I_ERI_F2yz_Px_Fy2z_S_a-1*I_ERI_Pz_Px_Fy2z_S;
  abcd[329] = 2.0E0*I_ERI_Fy2z_Px_Fy2z_S_a;
  abcd[330] = 2.0E0*I_ERI_F2xy_Py_Fy2z_S_a;
  abcd[331] = 2.0E0*I_ERI_Fx2y_Py_Fy2z_S_a-1*I_ERI_Px_Py_Fy2z_S;
  abcd[332] = 2.0E0*I_ERI_Fxyz_Py_Fy2z_S_a;
  abcd[333] = 2.0E0*I_ERI_F3y_Py_Fy2z_S_a-2*I_ERI_Py_Py_Fy2z_S;
  abcd[334] = 2.0E0*I_ERI_F2yz_Py_Fy2z_S_a-1*I_ERI_Pz_Py_Fy2z_S;
  abcd[335] = 2.0E0*I_ERI_Fy2z_Py_Fy2z_S_a;
  abcd[336] = 2.0E0*I_ERI_F2xy_Pz_Fy2z_S_a;
  abcd[337] = 2.0E0*I_ERI_Fx2y_Pz_Fy2z_S_a-1*I_ERI_Px_Pz_Fy2z_S;
  abcd[338] = 2.0E0*I_ERI_Fxyz_Pz_Fy2z_S_a;
  abcd[339] = 2.0E0*I_ERI_F3y_Pz_Fy2z_S_a-2*I_ERI_Py_Pz_Fy2z_S;
  abcd[340] = 2.0E0*I_ERI_F2yz_Pz_Fy2z_S_a-1*I_ERI_Pz_Pz_Fy2z_S;
  abcd[341] = 2.0E0*I_ERI_Fy2z_Pz_Fy2z_S_a;
  abcd[342] = 2.0E0*I_ERI_F2xy_Px_F3z_S_a;
  abcd[343] = 2.0E0*I_ERI_Fx2y_Px_F3z_S_a-1*I_ERI_Px_Px_F3z_S;
  abcd[344] = 2.0E0*I_ERI_Fxyz_Px_F3z_S_a;
  abcd[345] = 2.0E0*I_ERI_F3y_Px_F3z_S_a-2*I_ERI_Py_Px_F3z_S;
  abcd[346] = 2.0E0*I_ERI_F2yz_Px_F3z_S_a-1*I_ERI_Pz_Px_F3z_S;
  abcd[347] = 2.0E0*I_ERI_Fy2z_Px_F3z_S_a;
  abcd[348] = 2.0E0*I_ERI_F2xy_Py_F3z_S_a;
  abcd[349] = 2.0E0*I_ERI_Fx2y_Py_F3z_S_a-1*I_ERI_Px_Py_F3z_S;
  abcd[350] = 2.0E0*I_ERI_Fxyz_Py_F3z_S_a;
  abcd[351] = 2.0E0*I_ERI_F3y_Py_F3z_S_a-2*I_ERI_Py_Py_F3z_S;
  abcd[352] = 2.0E0*I_ERI_F2yz_Py_F3z_S_a-1*I_ERI_Pz_Py_F3z_S;
  abcd[353] = 2.0E0*I_ERI_Fy2z_Py_F3z_S_a;
  abcd[354] = 2.0E0*I_ERI_F2xy_Pz_F3z_S_a;
  abcd[355] = 2.0E0*I_ERI_Fx2y_Pz_F3z_S_a-1*I_ERI_Px_Pz_F3z_S;
  abcd[356] = 2.0E0*I_ERI_Fxyz_Pz_F3z_S_a;
  abcd[357] = 2.0E0*I_ERI_F3y_Pz_F3z_S_a-2*I_ERI_Py_Pz_F3z_S;
  abcd[358] = 2.0E0*I_ERI_F2yz_Pz_F3z_S_a-1*I_ERI_Pz_Pz_F3z_S;
  abcd[359] = 2.0E0*I_ERI_Fy2z_Pz_F3z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_a
   * RHS shell quartet name: SQ_ERI_P_P_F_S
   ************************************************************/
  abcd[360] = 2.0E0*I_ERI_F2xz_Px_F3x_S_a;
  abcd[361] = 2.0E0*I_ERI_Fxyz_Px_F3x_S_a;
  abcd[362] = 2.0E0*I_ERI_Fx2z_Px_F3x_S_a-1*I_ERI_Px_Px_F3x_S;
  abcd[363] = 2.0E0*I_ERI_F2yz_Px_F3x_S_a;
  abcd[364] = 2.0E0*I_ERI_Fy2z_Px_F3x_S_a-1*I_ERI_Py_Px_F3x_S;
  abcd[365] = 2.0E0*I_ERI_F3z_Px_F3x_S_a-2*I_ERI_Pz_Px_F3x_S;
  abcd[366] = 2.0E0*I_ERI_F2xz_Py_F3x_S_a;
  abcd[367] = 2.0E0*I_ERI_Fxyz_Py_F3x_S_a;
  abcd[368] = 2.0E0*I_ERI_Fx2z_Py_F3x_S_a-1*I_ERI_Px_Py_F3x_S;
  abcd[369] = 2.0E0*I_ERI_F2yz_Py_F3x_S_a;
  abcd[370] = 2.0E0*I_ERI_Fy2z_Py_F3x_S_a-1*I_ERI_Py_Py_F3x_S;
  abcd[371] = 2.0E0*I_ERI_F3z_Py_F3x_S_a-2*I_ERI_Pz_Py_F3x_S;
  abcd[372] = 2.0E0*I_ERI_F2xz_Pz_F3x_S_a;
  abcd[373] = 2.0E0*I_ERI_Fxyz_Pz_F3x_S_a;
  abcd[374] = 2.0E0*I_ERI_Fx2z_Pz_F3x_S_a-1*I_ERI_Px_Pz_F3x_S;
  abcd[375] = 2.0E0*I_ERI_F2yz_Pz_F3x_S_a;
  abcd[376] = 2.0E0*I_ERI_Fy2z_Pz_F3x_S_a-1*I_ERI_Py_Pz_F3x_S;
  abcd[377] = 2.0E0*I_ERI_F3z_Pz_F3x_S_a-2*I_ERI_Pz_Pz_F3x_S;
  abcd[378] = 2.0E0*I_ERI_F2xz_Px_F2xy_S_a;
  abcd[379] = 2.0E0*I_ERI_Fxyz_Px_F2xy_S_a;
  abcd[380] = 2.0E0*I_ERI_Fx2z_Px_F2xy_S_a-1*I_ERI_Px_Px_F2xy_S;
  abcd[381] = 2.0E0*I_ERI_F2yz_Px_F2xy_S_a;
  abcd[382] = 2.0E0*I_ERI_Fy2z_Px_F2xy_S_a-1*I_ERI_Py_Px_F2xy_S;
  abcd[383] = 2.0E0*I_ERI_F3z_Px_F2xy_S_a-2*I_ERI_Pz_Px_F2xy_S;
  abcd[384] = 2.0E0*I_ERI_F2xz_Py_F2xy_S_a;
  abcd[385] = 2.0E0*I_ERI_Fxyz_Py_F2xy_S_a;
  abcd[386] = 2.0E0*I_ERI_Fx2z_Py_F2xy_S_a-1*I_ERI_Px_Py_F2xy_S;
  abcd[387] = 2.0E0*I_ERI_F2yz_Py_F2xy_S_a;
  abcd[388] = 2.0E0*I_ERI_Fy2z_Py_F2xy_S_a-1*I_ERI_Py_Py_F2xy_S;
  abcd[389] = 2.0E0*I_ERI_F3z_Py_F2xy_S_a-2*I_ERI_Pz_Py_F2xy_S;
  abcd[390] = 2.0E0*I_ERI_F2xz_Pz_F2xy_S_a;
  abcd[391] = 2.0E0*I_ERI_Fxyz_Pz_F2xy_S_a;
  abcd[392] = 2.0E0*I_ERI_Fx2z_Pz_F2xy_S_a-1*I_ERI_Px_Pz_F2xy_S;
  abcd[393] = 2.0E0*I_ERI_F2yz_Pz_F2xy_S_a;
  abcd[394] = 2.0E0*I_ERI_Fy2z_Pz_F2xy_S_a-1*I_ERI_Py_Pz_F2xy_S;
  abcd[395] = 2.0E0*I_ERI_F3z_Pz_F2xy_S_a-2*I_ERI_Pz_Pz_F2xy_S;
  abcd[396] = 2.0E0*I_ERI_F2xz_Px_F2xz_S_a;
  abcd[397] = 2.0E0*I_ERI_Fxyz_Px_F2xz_S_a;
  abcd[398] = 2.0E0*I_ERI_Fx2z_Px_F2xz_S_a-1*I_ERI_Px_Px_F2xz_S;
  abcd[399] = 2.0E0*I_ERI_F2yz_Px_F2xz_S_a;
  abcd[400] = 2.0E0*I_ERI_Fy2z_Px_F2xz_S_a-1*I_ERI_Py_Px_F2xz_S;
  abcd[401] = 2.0E0*I_ERI_F3z_Px_F2xz_S_a-2*I_ERI_Pz_Px_F2xz_S;
  abcd[402] = 2.0E0*I_ERI_F2xz_Py_F2xz_S_a;
  abcd[403] = 2.0E0*I_ERI_Fxyz_Py_F2xz_S_a;
  abcd[404] = 2.0E0*I_ERI_Fx2z_Py_F2xz_S_a-1*I_ERI_Px_Py_F2xz_S;
  abcd[405] = 2.0E0*I_ERI_F2yz_Py_F2xz_S_a;
  abcd[406] = 2.0E0*I_ERI_Fy2z_Py_F2xz_S_a-1*I_ERI_Py_Py_F2xz_S;
  abcd[407] = 2.0E0*I_ERI_F3z_Py_F2xz_S_a-2*I_ERI_Pz_Py_F2xz_S;
  abcd[408] = 2.0E0*I_ERI_F2xz_Pz_F2xz_S_a;
  abcd[409] = 2.0E0*I_ERI_Fxyz_Pz_F2xz_S_a;
  abcd[410] = 2.0E0*I_ERI_Fx2z_Pz_F2xz_S_a-1*I_ERI_Px_Pz_F2xz_S;
  abcd[411] = 2.0E0*I_ERI_F2yz_Pz_F2xz_S_a;
  abcd[412] = 2.0E0*I_ERI_Fy2z_Pz_F2xz_S_a-1*I_ERI_Py_Pz_F2xz_S;
  abcd[413] = 2.0E0*I_ERI_F3z_Pz_F2xz_S_a-2*I_ERI_Pz_Pz_F2xz_S;
  abcd[414] = 2.0E0*I_ERI_F2xz_Px_Fx2y_S_a;
  abcd[415] = 2.0E0*I_ERI_Fxyz_Px_Fx2y_S_a;
  abcd[416] = 2.0E0*I_ERI_Fx2z_Px_Fx2y_S_a-1*I_ERI_Px_Px_Fx2y_S;
  abcd[417] = 2.0E0*I_ERI_F2yz_Px_Fx2y_S_a;
  abcd[418] = 2.0E0*I_ERI_Fy2z_Px_Fx2y_S_a-1*I_ERI_Py_Px_Fx2y_S;
  abcd[419] = 2.0E0*I_ERI_F3z_Px_Fx2y_S_a-2*I_ERI_Pz_Px_Fx2y_S;
  abcd[420] = 2.0E0*I_ERI_F2xz_Py_Fx2y_S_a;
  abcd[421] = 2.0E0*I_ERI_Fxyz_Py_Fx2y_S_a;
  abcd[422] = 2.0E0*I_ERI_Fx2z_Py_Fx2y_S_a-1*I_ERI_Px_Py_Fx2y_S;
  abcd[423] = 2.0E0*I_ERI_F2yz_Py_Fx2y_S_a;
  abcd[424] = 2.0E0*I_ERI_Fy2z_Py_Fx2y_S_a-1*I_ERI_Py_Py_Fx2y_S;
  abcd[425] = 2.0E0*I_ERI_F3z_Py_Fx2y_S_a-2*I_ERI_Pz_Py_Fx2y_S;
  abcd[426] = 2.0E0*I_ERI_F2xz_Pz_Fx2y_S_a;
  abcd[427] = 2.0E0*I_ERI_Fxyz_Pz_Fx2y_S_a;
  abcd[428] = 2.0E0*I_ERI_Fx2z_Pz_Fx2y_S_a-1*I_ERI_Px_Pz_Fx2y_S;
  abcd[429] = 2.0E0*I_ERI_F2yz_Pz_Fx2y_S_a;
  abcd[430] = 2.0E0*I_ERI_Fy2z_Pz_Fx2y_S_a-1*I_ERI_Py_Pz_Fx2y_S;
  abcd[431] = 2.0E0*I_ERI_F3z_Pz_Fx2y_S_a-2*I_ERI_Pz_Pz_Fx2y_S;
  abcd[432] = 2.0E0*I_ERI_F2xz_Px_Fxyz_S_a;
  abcd[433] = 2.0E0*I_ERI_Fxyz_Px_Fxyz_S_a;
  abcd[434] = 2.0E0*I_ERI_Fx2z_Px_Fxyz_S_a-1*I_ERI_Px_Px_Fxyz_S;
  abcd[435] = 2.0E0*I_ERI_F2yz_Px_Fxyz_S_a;
  abcd[436] = 2.0E0*I_ERI_Fy2z_Px_Fxyz_S_a-1*I_ERI_Py_Px_Fxyz_S;
  abcd[437] = 2.0E0*I_ERI_F3z_Px_Fxyz_S_a-2*I_ERI_Pz_Px_Fxyz_S;
  abcd[438] = 2.0E0*I_ERI_F2xz_Py_Fxyz_S_a;
  abcd[439] = 2.0E0*I_ERI_Fxyz_Py_Fxyz_S_a;
  abcd[440] = 2.0E0*I_ERI_Fx2z_Py_Fxyz_S_a-1*I_ERI_Px_Py_Fxyz_S;
  abcd[441] = 2.0E0*I_ERI_F2yz_Py_Fxyz_S_a;
  abcd[442] = 2.0E0*I_ERI_Fy2z_Py_Fxyz_S_a-1*I_ERI_Py_Py_Fxyz_S;
  abcd[443] = 2.0E0*I_ERI_F3z_Py_Fxyz_S_a-2*I_ERI_Pz_Py_Fxyz_S;
  abcd[444] = 2.0E0*I_ERI_F2xz_Pz_Fxyz_S_a;
  abcd[445] = 2.0E0*I_ERI_Fxyz_Pz_Fxyz_S_a;
  abcd[446] = 2.0E0*I_ERI_Fx2z_Pz_Fxyz_S_a-1*I_ERI_Px_Pz_Fxyz_S;
  abcd[447] = 2.0E0*I_ERI_F2yz_Pz_Fxyz_S_a;
  abcd[448] = 2.0E0*I_ERI_Fy2z_Pz_Fxyz_S_a-1*I_ERI_Py_Pz_Fxyz_S;
  abcd[449] = 2.0E0*I_ERI_F3z_Pz_Fxyz_S_a-2*I_ERI_Pz_Pz_Fxyz_S;
  abcd[450] = 2.0E0*I_ERI_F2xz_Px_Fx2z_S_a;
  abcd[451] = 2.0E0*I_ERI_Fxyz_Px_Fx2z_S_a;
  abcd[452] = 2.0E0*I_ERI_Fx2z_Px_Fx2z_S_a-1*I_ERI_Px_Px_Fx2z_S;
  abcd[453] = 2.0E0*I_ERI_F2yz_Px_Fx2z_S_a;
  abcd[454] = 2.0E0*I_ERI_Fy2z_Px_Fx2z_S_a-1*I_ERI_Py_Px_Fx2z_S;
  abcd[455] = 2.0E0*I_ERI_F3z_Px_Fx2z_S_a-2*I_ERI_Pz_Px_Fx2z_S;
  abcd[456] = 2.0E0*I_ERI_F2xz_Py_Fx2z_S_a;
  abcd[457] = 2.0E0*I_ERI_Fxyz_Py_Fx2z_S_a;
  abcd[458] = 2.0E0*I_ERI_Fx2z_Py_Fx2z_S_a-1*I_ERI_Px_Py_Fx2z_S;
  abcd[459] = 2.0E0*I_ERI_F2yz_Py_Fx2z_S_a;
  abcd[460] = 2.0E0*I_ERI_Fy2z_Py_Fx2z_S_a-1*I_ERI_Py_Py_Fx2z_S;
  abcd[461] = 2.0E0*I_ERI_F3z_Py_Fx2z_S_a-2*I_ERI_Pz_Py_Fx2z_S;
  abcd[462] = 2.0E0*I_ERI_F2xz_Pz_Fx2z_S_a;
  abcd[463] = 2.0E0*I_ERI_Fxyz_Pz_Fx2z_S_a;
  abcd[464] = 2.0E0*I_ERI_Fx2z_Pz_Fx2z_S_a-1*I_ERI_Px_Pz_Fx2z_S;
  abcd[465] = 2.0E0*I_ERI_F2yz_Pz_Fx2z_S_a;
  abcd[466] = 2.0E0*I_ERI_Fy2z_Pz_Fx2z_S_a-1*I_ERI_Py_Pz_Fx2z_S;
  abcd[467] = 2.0E0*I_ERI_F3z_Pz_Fx2z_S_a-2*I_ERI_Pz_Pz_Fx2z_S;
  abcd[468] = 2.0E0*I_ERI_F2xz_Px_F3y_S_a;
  abcd[469] = 2.0E0*I_ERI_Fxyz_Px_F3y_S_a;
  abcd[470] = 2.0E0*I_ERI_Fx2z_Px_F3y_S_a-1*I_ERI_Px_Px_F3y_S;
  abcd[471] = 2.0E0*I_ERI_F2yz_Px_F3y_S_a;
  abcd[472] = 2.0E0*I_ERI_Fy2z_Px_F3y_S_a-1*I_ERI_Py_Px_F3y_S;
  abcd[473] = 2.0E0*I_ERI_F3z_Px_F3y_S_a-2*I_ERI_Pz_Px_F3y_S;
  abcd[474] = 2.0E0*I_ERI_F2xz_Py_F3y_S_a;
  abcd[475] = 2.0E0*I_ERI_Fxyz_Py_F3y_S_a;
  abcd[476] = 2.0E0*I_ERI_Fx2z_Py_F3y_S_a-1*I_ERI_Px_Py_F3y_S;
  abcd[477] = 2.0E0*I_ERI_F2yz_Py_F3y_S_a;
  abcd[478] = 2.0E0*I_ERI_Fy2z_Py_F3y_S_a-1*I_ERI_Py_Py_F3y_S;
  abcd[479] = 2.0E0*I_ERI_F3z_Py_F3y_S_a-2*I_ERI_Pz_Py_F3y_S;
  abcd[480] = 2.0E0*I_ERI_F2xz_Pz_F3y_S_a;
  abcd[481] = 2.0E0*I_ERI_Fxyz_Pz_F3y_S_a;
  abcd[482] = 2.0E0*I_ERI_Fx2z_Pz_F3y_S_a-1*I_ERI_Px_Pz_F3y_S;
  abcd[483] = 2.0E0*I_ERI_F2yz_Pz_F3y_S_a;
  abcd[484] = 2.0E0*I_ERI_Fy2z_Pz_F3y_S_a-1*I_ERI_Py_Pz_F3y_S;
  abcd[485] = 2.0E0*I_ERI_F3z_Pz_F3y_S_a-2*I_ERI_Pz_Pz_F3y_S;
  abcd[486] = 2.0E0*I_ERI_F2xz_Px_F2yz_S_a;
  abcd[487] = 2.0E0*I_ERI_Fxyz_Px_F2yz_S_a;
  abcd[488] = 2.0E0*I_ERI_Fx2z_Px_F2yz_S_a-1*I_ERI_Px_Px_F2yz_S;
  abcd[489] = 2.0E0*I_ERI_F2yz_Px_F2yz_S_a;
  abcd[490] = 2.0E0*I_ERI_Fy2z_Px_F2yz_S_a-1*I_ERI_Py_Px_F2yz_S;
  abcd[491] = 2.0E0*I_ERI_F3z_Px_F2yz_S_a-2*I_ERI_Pz_Px_F2yz_S;
  abcd[492] = 2.0E0*I_ERI_F2xz_Py_F2yz_S_a;
  abcd[493] = 2.0E0*I_ERI_Fxyz_Py_F2yz_S_a;
  abcd[494] = 2.0E0*I_ERI_Fx2z_Py_F2yz_S_a-1*I_ERI_Px_Py_F2yz_S;
  abcd[495] = 2.0E0*I_ERI_F2yz_Py_F2yz_S_a;
  abcd[496] = 2.0E0*I_ERI_Fy2z_Py_F2yz_S_a-1*I_ERI_Py_Py_F2yz_S;
  abcd[497] = 2.0E0*I_ERI_F3z_Py_F2yz_S_a-2*I_ERI_Pz_Py_F2yz_S;
  abcd[498] = 2.0E0*I_ERI_F2xz_Pz_F2yz_S_a;
  abcd[499] = 2.0E0*I_ERI_Fxyz_Pz_F2yz_S_a;
  abcd[500] = 2.0E0*I_ERI_Fx2z_Pz_F2yz_S_a-1*I_ERI_Px_Pz_F2yz_S;
  abcd[501] = 2.0E0*I_ERI_F2yz_Pz_F2yz_S_a;
  abcd[502] = 2.0E0*I_ERI_Fy2z_Pz_F2yz_S_a-1*I_ERI_Py_Pz_F2yz_S;
  abcd[503] = 2.0E0*I_ERI_F3z_Pz_F2yz_S_a-2*I_ERI_Pz_Pz_F2yz_S;
  abcd[504] = 2.0E0*I_ERI_F2xz_Px_Fy2z_S_a;
  abcd[505] = 2.0E0*I_ERI_Fxyz_Px_Fy2z_S_a;
  abcd[506] = 2.0E0*I_ERI_Fx2z_Px_Fy2z_S_a-1*I_ERI_Px_Px_Fy2z_S;
  abcd[507] = 2.0E0*I_ERI_F2yz_Px_Fy2z_S_a;
  abcd[508] = 2.0E0*I_ERI_Fy2z_Px_Fy2z_S_a-1*I_ERI_Py_Px_Fy2z_S;
  abcd[509] = 2.0E0*I_ERI_F3z_Px_Fy2z_S_a-2*I_ERI_Pz_Px_Fy2z_S;
  abcd[510] = 2.0E0*I_ERI_F2xz_Py_Fy2z_S_a;
  abcd[511] = 2.0E0*I_ERI_Fxyz_Py_Fy2z_S_a;
  abcd[512] = 2.0E0*I_ERI_Fx2z_Py_Fy2z_S_a-1*I_ERI_Px_Py_Fy2z_S;
  abcd[513] = 2.0E0*I_ERI_F2yz_Py_Fy2z_S_a;
  abcd[514] = 2.0E0*I_ERI_Fy2z_Py_Fy2z_S_a-1*I_ERI_Py_Py_Fy2z_S;
  abcd[515] = 2.0E0*I_ERI_F3z_Py_Fy2z_S_a-2*I_ERI_Pz_Py_Fy2z_S;
  abcd[516] = 2.0E0*I_ERI_F2xz_Pz_Fy2z_S_a;
  abcd[517] = 2.0E0*I_ERI_Fxyz_Pz_Fy2z_S_a;
  abcd[518] = 2.0E0*I_ERI_Fx2z_Pz_Fy2z_S_a-1*I_ERI_Px_Pz_Fy2z_S;
  abcd[519] = 2.0E0*I_ERI_F2yz_Pz_Fy2z_S_a;
  abcd[520] = 2.0E0*I_ERI_Fy2z_Pz_Fy2z_S_a-1*I_ERI_Py_Pz_Fy2z_S;
  abcd[521] = 2.0E0*I_ERI_F3z_Pz_Fy2z_S_a-2*I_ERI_Pz_Pz_Fy2z_S;
  abcd[522] = 2.0E0*I_ERI_F2xz_Px_F3z_S_a;
  abcd[523] = 2.0E0*I_ERI_Fxyz_Px_F3z_S_a;
  abcd[524] = 2.0E0*I_ERI_Fx2z_Px_F3z_S_a-1*I_ERI_Px_Px_F3z_S;
  abcd[525] = 2.0E0*I_ERI_F2yz_Px_F3z_S_a;
  abcd[526] = 2.0E0*I_ERI_Fy2z_Px_F3z_S_a-1*I_ERI_Py_Px_F3z_S;
  abcd[527] = 2.0E0*I_ERI_F3z_Px_F3z_S_a-2*I_ERI_Pz_Px_F3z_S;
  abcd[528] = 2.0E0*I_ERI_F2xz_Py_F3z_S_a;
  abcd[529] = 2.0E0*I_ERI_Fxyz_Py_F3z_S_a;
  abcd[530] = 2.0E0*I_ERI_Fx2z_Py_F3z_S_a-1*I_ERI_Px_Py_F3z_S;
  abcd[531] = 2.0E0*I_ERI_F2yz_Py_F3z_S_a;
  abcd[532] = 2.0E0*I_ERI_Fy2z_Py_F3z_S_a-1*I_ERI_Py_Py_F3z_S;
  abcd[533] = 2.0E0*I_ERI_F3z_Py_F3z_S_a-2*I_ERI_Pz_Py_F3z_S;
  abcd[534] = 2.0E0*I_ERI_F2xz_Pz_F3z_S_a;
  abcd[535] = 2.0E0*I_ERI_Fxyz_Pz_F3z_S_a;
  abcd[536] = 2.0E0*I_ERI_Fx2z_Pz_F3z_S_a-1*I_ERI_Px_Pz_F3z_S;
  abcd[537] = 2.0E0*I_ERI_F2yz_Pz_F3z_S_a;
  abcd[538] = 2.0E0*I_ERI_Fy2z_Pz_F3z_S_a-1*I_ERI_Py_Pz_F3z_S;
  abcd[539] = 2.0E0*I_ERI_F3z_Pz_F3z_S_a-2*I_ERI_Pz_Pz_F3z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_b
   * RHS shell quartet name: SQ_ERI_D_S_F_S
   ************************************************************/
  abcd[540] = 2.0E0*I_ERI_D2x_D2x_F3x_S_b-1*I_ERI_D2x_S_F3x_S;
  abcd[541] = 2.0E0*I_ERI_Dxy_D2x_F3x_S_b-1*I_ERI_Dxy_S_F3x_S;
  abcd[542] = 2.0E0*I_ERI_Dxz_D2x_F3x_S_b-1*I_ERI_Dxz_S_F3x_S;
  abcd[543] = 2.0E0*I_ERI_D2y_D2x_F3x_S_b-1*I_ERI_D2y_S_F3x_S;
  abcd[544] = 2.0E0*I_ERI_Dyz_D2x_F3x_S_b-1*I_ERI_Dyz_S_F3x_S;
  abcd[545] = 2.0E0*I_ERI_D2z_D2x_F3x_S_b-1*I_ERI_D2z_S_F3x_S;
  abcd[546] = 2.0E0*I_ERI_D2x_Dxy_F3x_S_b;
  abcd[547] = 2.0E0*I_ERI_Dxy_Dxy_F3x_S_b;
  abcd[548] = 2.0E0*I_ERI_Dxz_Dxy_F3x_S_b;
  abcd[549] = 2.0E0*I_ERI_D2y_Dxy_F3x_S_b;
  abcd[550] = 2.0E0*I_ERI_Dyz_Dxy_F3x_S_b;
  abcd[551] = 2.0E0*I_ERI_D2z_Dxy_F3x_S_b;
  abcd[552] = 2.0E0*I_ERI_D2x_Dxz_F3x_S_b;
  abcd[553] = 2.0E0*I_ERI_Dxy_Dxz_F3x_S_b;
  abcd[554] = 2.0E0*I_ERI_Dxz_Dxz_F3x_S_b;
  abcd[555] = 2.0E0*I_ERI_D2y_Dxz_F3x_S_b;
  abcd[556] = 2.0E0*I_ERI_Dyz_Dxz_F3x_S_b;
  abcd[557] = 2.0E0*I_ERI_D2z_Dxz_F3x_S_b;
  abcd[558] = 2.0E0*I_ERI_D2x_D2x_F2xy_S_b-1*I_ERI_D2x_S_F2xy_S;
  abcd[559] = 2.0E0*I_ERI_Dxy_D2x_F2xy_S_b-1*I_ERI_Dxy_S_F2xy_S;
  abcd[560] = 2.0E0*I_ERI_Dxz_D2x_F2xy_S_b-1*I_ERI_Dxz_S_F2xy_S;
  abcd[561] = 2.0E0*I_ERI_D2y_D2x_F2xy_S_b-1*I_ERI_D2y_S_F2xy_S;
  abcd[562] = 2.0E0*I_ERI_Dyz_D2x_F2xy_S_b-1*I_ERI_Dyz_S_F2xy_S;
  abcd[563] = 2.0E0*I_ERI_D2z_D2x_F2xy_S_b-1*I_ERI_D2z_S_F2xy_S;
  abcd[564] = 2.0E0*I_ERI_D2x_Dxy_F2xy_S_b;
  abcd[565] = 2.0E0*I_ERI_Dxy_Dxy_F2xy_S_b;
  abcd[566] = 2.0E0*I_ERI_Dxz_Dxy_F2xy_S_b;
  abcd[567] = 2.0E0*I_ERI_D2y_Dxy_F2xy_S_b;
  abcd[568] = 2.0E0*I_ERI_Dyz_Dxy_F2xy_S_b;
  abcd[569] = 2.0E0*I_ERI_D2z_Dxy_F2xy_S_b;
  abcd[570] = 2.0E0*I_ERI_D2x_Dxz_F2xy_S_b;
  abcd[571] = 2.0E0*I_ERI_Dxy_Dxz_F2xy_S_b;
  abcd[572] = 2.0E0*I_ERI_Dxz_Dxz_F2xy_S_b;
  abcd[573] = 2.0E0*I_ERI_D2y_Dxz_F2xy_S_b;
  abcd[574] = 2.0E0*I_ERI_Dyz_Dxz_F2xy_S_b;
  abcd[575] = 2.0E0*I_ERI_D2z_Dxz_F2xy_S_b;
  abcd[576] = 2.0E0*I_ERI_D2x_D2x_F2xz_S_b-1*I_ERI_D2x_S_F2xz_S;
  abcd[577] = 2.0E0*I_ERI_Dxy_D2x_F2xz_S_b-1*I_ERI_Dxy_S_F2xz_S;
  abcd[578] = 2.0E0*I_ERI_Dxz_D2x_F2xz_S_b-1*I_ERI_Dxz_S_F2xz_S;
  abcd[579] = 2.0E0*I_ERI_D2y_D2x_F2xz_S_b-1*I_ERI_D2y_S_F2xz_S;
  abcd[580] = 2.0E0*I_ERI_Dyz_D2x_F2xz_S_b-1*I_ERI_Dyz_S_F2xz_S;
  abcd[581] = 2.0E0*I_ERI_D2z_D2x_F2xz_S_b-1*I_ERI_D2z_S_F2xz_S;
  abcd[582] = 2.0E0*I_ERI_D2x_Dxy_F2xz_S_b;
  abcd[583] = 2.0E0*I_ERI_Dxy_Dxy_F2xz_S_b;
  abcd[584] = 2.0E0*I_ERI_Dxz_Dxy_F2xz_S_b;
  abcd[585] = 2.0E0*I_ERI_D2y_Dxy_F2xz_S_b;
  abcd[586] = 2.0E0*I_ERI_Dyz_Dxy_F2xz_S_b;
  abcd[587] = 2.0E0*I_ERI_D2z_Dxy_F2xz_S_b;
  abcd[588] = 2.0E0*I_ERI_D2x_Dxz_F2xz_S_b;
  abcd[589] = 2.0E0*I_ERI_Dxy_Dxz_F2xz_S_b;
  abcd[590] = 2.0E0*I_ERI_Dxz_Dxz_F2xz_S_b;
  abcd[591] = 2.0E0*I_ERI_D2y_Dxz_F2xz_S_b;
  abcd[592] = 2.0E0*I_ERI_Dyz_Dxz_F2xz_S_b;
  abcd[593] = 2.0E0*I_ERI_D2z_Dxz_F2xz_S_b;
  abcd[594] = 2.0E0*I_ERI_D2x_D2x_Fx2y_S_b-1*I_ERI_D2x_S_Fx2y_S;
  abcd[595] = 2.0E0*I_ERI_Dxy_D2x_Fx2y_S_b-1*I_ERI_Dxy_S_Fx2y_S;
  abcd[596] = 2.0E0*I_ERI_Dxz_D2x_Fx2y_S_b-1*I_ERI_Dxz_S_Fx2y_S;
  abcd[597] = 2.0E0*I_ERI_D2y_D2x_Fx2y_S_b-1*I_ERI_D2y_S_Fx2y_S;
  abcd[598] = 2.0E0*I_ERI_Dyz_D2x_Fx2y_S_b-1*I_ERI_Dyz_S_Fx2y_S;
  abcd[599] = 2.0E0*I_ERI_D2z_D2x_Fx2y_S_b-1*I_ERI_D2z_S_Fx2y_S;
  abcd[600] = 2.0E0*I_ERI_D2x_Dxy_Fx2y_S_b;
  abcd[601] = 2.0E0*I_ERI_Dxy_Dxy_Fx2y_S_b;
  abcd[602] = 2.0E0*I_ERI_Dxz_Dxy_Fx2y_S_b;
  abcd[603] = 2.0E0*I_ERI_D2y_Dxy_Fx2y_S_b;
  abcd[604] = 2.0E0*I_ERI_Dyz_Dxy_Fx2y_S_b;
  abcd[605] = 2.0E0*I_ERI_D2z_Dxy_Fx2y_S_b;
  abcd[606] = 2.0E0*I_ERI_D2x_Dxz_Fx2y_S_b;
  abcd[607] = 2.0E0*I_ERI_Dxy_Dxz_Fx2y_S_b;
  abcd[608] = 2.0E0*I_ERI_Dxz_Dxz_Fx2y_S_b;
  abcd[609] = 2.0E0*I_ERI_D2y_Dxz_Fx2y_S_b;
  abcd[610] = 2.0E0*I_ERI_Dyz_Dxz_Fx2y_S_b;
  abcd[611] = 2.0E0*I_ERI_D2z_Dxz_Fx2y_S_b;
  abcd[612] = 2.0E0*I_ERI_D2x_D2x_Fxyz_S_b-1*I_ERI_D2x_S_Fxyz_S;
  abcd[613] = 2.0E0*I_ERI_Dxy_D2x_Fxyz_S_b-1*I_ERI_Dxy_S_Fxyz_S;
  abcd[614] = 2.0E0*I_ERI_Dxz_D2x_Fxyz_S_b-1*I_ERI_Dxz_S_Fxyz_S;
  abcd[615] = 2.0E0*I_ERI_D2y_D2x_Fxyz_S_b-1*I_ERI_D2y_S_Fxyz_S;
  abcd[616] = 2.0E0*I_ERI_Dyz_D2x_Fxyz_S_b-1*I_ERI_Dyz_S_Fxyz_S;
  abcd[617] = 2.0E0*I_ERI_D2z_D2x_Fxyz_S_b-1*I_ERI_D2z_S_Fxyz_S;
  abcd[618] = 2.0E0*I_ERI_D2x_Dxy_Fxyz_S_b;
  abcd[619] = 2.0E0*I_ERI_Dxy_Dxy_Fxyz_S_b;
  abcd[620] = 2.0E0*I_ERI_Dxz_Dxy_Fxyz_S_b;
  abcd[621] = 2.0E0*I_ERI_D2y_Dxy_Fxyz_S_b;
  abcd[622] = 2.0E0*I_ERI_Dyz_Dxy_Fxyz_S_b;
  abcd[623] = 2.0E0*I_ERI_D2z_Dxy_Fxyz_S_b;
  abcd[624] = 2.0E0*I_ERI_D2x_Dxz_Fxyz_S_b;
  abcd[625] = 2.0E0*I_ERI_Dxy_Dxz_Fxyz_S_b;
  abcd[626] = 2.0E0*I_ERI_Dxz_Dxz_Fxyz_S_b;
  abcd[627] = 2.0E0*I_ERI_D2y_Dxz_Fxyz_S_b;
  abcd[628] = 2.0E0*I_ERI_Dyz_Dxz_Fxyz_S_b;
  abcd[629] = 2.0E0*I_ERI_D2z_Dxz_Fxyz_S_b;
  abcd[630] = 2.0E0*I_ERI_D2x_D2x_Fx2z_S_b-1*I_ERI_D2x_S_Fx2z_S;
  abcd[631] = 2.0E0*I_ERI_Dxy_D2x_Fx2z_S_b-1*I_ERI_Dxy_S_Fx2z_S;
  abcd[632] = 2.0E0*I_ERI_Dxz_D2x_Fx2z_S_b-1*I_ERI_Dxz_S_Fx2z_S;
  abcd[633] = 2.0E0*I_ERI_D2y_D2x_Fx2z_S_b-1*I_ERI_D2y_S_Fx2z_S;
  abcd[634] = 2.0E0*I_ERI_Dyz_D2x_Fx2z_S_b-1*I_ERI_Dyz_S_Fx2z_S;
  abcd[635] = 2.0E0*I_ERI_D2z_D2x_Fx2z_S_b-1*I_ERI_D2z_S_Fx2z_S;
  abcd[636] = 2.0E0*I_ERI_D2x_Dxy_Fx2z_S_b;
  abcd[637] = 2.0E0*I_ERI_Dxy_Dxy_Fx2z_S_b;
  abcd[638] = 2.0E0*I_ERI_Dxz_Dxy_Fx2z_S_b;
  abcd[639] = 2.0E0*I_ERI_D2y_Dxy_Fx2z_S_b;
  abcd[640] = 2.0E0*I_ERI_Dyz_Dxy_Fx2z_S_b;
  abcd[641] = 2.0E0*I_ERI_D2z_Dxy_Fx2z_S_b;
  abcd[642] = 2.0E0*I_ERI_D2x_Dxz_Fx2z_S_b;
  abcd[643] = 2.0E0*I_ERI_Dxy_Dxz_Fx2z_S_b;
  abcd[644] = 2.0E0*I_ERI_Dxz_Dxz_Fx2z_S_b;
  abcd[645] = 2.0E0*I_ERI_D2y_Dxz_Fx2z_S_b;
  abcd[646] = 2.0E0*I_ERI_Dyz_Dxz_Fx2z_S_b;
  abcd[647] = 2.0E0*I_ERI_D2z_Dxz_Fx2z_S_b;
  abcd[648] = 2.0E0*I_ERI_D2x_D2x_F3y_S_b-1*I_ERI_D2x_S_F3y_S;
  abcd[649] = 2.0E0*I_ERI_Dxy_D2x_F3y_S_b-1*I_ERI_Dxy_S_F3y_S;
  abcd[650] = 2.0E0*I_ERI_Dxz_D2x_F3y_S_b-1*I_ERI_Dxz_S_F3y_S;
  abcd[651] = 2.0E0*I_ERI_D2y_D2x_F3y_S_b-1*I_ERI_D2y_S_F3y_S;
  abcd[652] = 2.0E0*I_ERI_Dyz_D2x_F3y_S_b-1*I_ERI_Dyz_S_F3y_S;
  abcd[653] = 2.0E0*I_ERI_D2z_D2x_F3y_S_b-1*I_ERI_D2z_S_F3y_S;
  abcd[654] = 2.0E0*I_ERI_D2x_Dxy_F3y_S_b;
  abcd[655] = 2.0E0*I_ERI_Dxy_Dxy_F3y_S_b;
  abcd[656] = 2.0E0*I_ERI_Dxz_Dxy_F3y_S_b;
  abcd[657] = 2.0E0*I_ERI_D2y_Dxy_F3y_S_b;
  abcd[658] = 2.0E0*I_ERI_Dyz_Dxy_F3y_S_b;
  abcd[659] = 2.0E0*I_ERI_D2z_Dxy_F3y_S_b;
  abcd[660] = 2.0E0*I_ERI_D2x_Dxz_F3y_S_b;
  abcd[661] = 2.0E0*I_ERI_Dxy_Dxz_F3y_S_b;
  abcd[662] = 2.0E0*I_ERI_Dxz_Dxz_F3y_S_b;
  abcd[663] = 2.0E0*I_ERI_D2y_Dxz_F3y_S_b;
  abcd[664] = 2.0E0*I_ERI_Dyz_Dxz_F3y_S_b;
  abcd[665] = 2.0E0*I_ERI_D2z_Dxz_F3y_S_b;
  abcd[666] = 2.0E0*I_ERI_D2x_D2x_F2yz_S_b-1*I_ERI_D2x_S_F2yz_S;
  abcd[667] = 2.0E0*I_ERI_Dxy_D2x_F2yz_S_b-1*I_ERI_Dxy_S_F2yz_S;
  abcd[668] = 2.0E0*I_ERI_Dxz_D2x_F2yz_S_b-1*I_ERI_Dxz_S_F2yz_S;
  abcd[669] = 2.0E0*I_ERI_D2y_D2x_F2yz_S_b-1*I_ERI_D2y_S_F2yz_S;
  abcd[670] = 2.0E0*I_ERI_Dyz_D2x_F2yz_S_b-1*I_ERI_Dyz_S_F2yz_S;
  abcd[671] = 2.0E0*I_ERI_D2z_D2x_F2yz_S_b-1*I_ERI_D2z_S_F2yz_S;
  abcd[672] = 2.0E0*I_ERI_D2x_Dxy_F2yz_S_b;
  abcd[673] = 2.0E0*I_ERI_Dxy_Dxy_F2yz_S_b;
  abcd[674] = 2.0E0*I_ERI_Dxz_Dxy_F2yz_S_b;
  abcd[675] = 2.0E0*I_ERI_D2y_Dxy_F2yz_S_b;
  abcd[676] = 2.0E0*I_ERI_Dyz_Dxy_F2yz_S_b;
  abcd[677] = 2.0E0*I_ERI_D2z_Dxy_F2yz_S_b;
  abcd[678] = 2.0E0*I_ERI_D2x_Dxz_F2yz_S_b;
  abcd[679] = 2.0E0*I_ERI_Dxy_Dxz_F2yz_S_b;
  abcd[680] = 2.0E0*I_ERI_Dxz_Dxz_F2yz_S_b;
  abcd[681] = 2.0E0*I_ERI_D2y_Dxz_F2yz_S_b;
  abcd[682] = 2.0E0*I_ERI_Dyz_Dxz_F2yz_S_b;
  abcd[683] = 2.0E0*I_ERI_D2z_Dxz_F2yz_S_b;
  abcd[684] = 2.0E0*I_ERI_D2x_D2x_Fy2z_S_b-1*I_ERI_D2x_S_Fy2z_S;
  abcd[685] = 2.0E0*I_ERI_Dxy_D2x_Fy2z_S_b-1*I_ERI_Dxy_S_Fy2z_S;
  abcd[686] = 2.0E0*I_ERI_Dxz_D2x_Fy2z_S_b-1*I_ERI_Dxz_S_Fy2z_S;
  abcd[687] = 2.0E0*I_ERI_D2y_D2x_Fy2z_S_b-1*I_ERI_D2y_S_Fy2z_S;
  abcd[688] = 2.0E0*I_ERI_Dyz_D2x_Fy2z_S_b-1*I_ERI_Dyz_S_Fy2z_S;
  abcd[689] = 2.0E0*I_ERI_D2z_D2x_Fy2z_S_b-1*I_ERI_D2z_S_Fy2z_S;
  abcd[690] = 2.0E0*I_ERI_D2x_Dxy_Fy2z_S_b;
  abcd[691] = 2.0E0*I_ERI_Dxy_Dxy_Fy2z_S_b;
  abcd[692] = 2.0E0*I_ERI_Dxz_Dxy_Fy2z_S_b;
  abcd[693] = 2.0E0*I_ERI_D2y_Dxy_Fy2z_S_b;
  abcd[694] = 2.0E0*I_ERI_Dyz_Dxy_Fy2z_S_b;
  abcd[695] = 2.0E0*I_ERI_D2z_Dxy_Fy2z_S_b;
  abcd[696] = 2.0E0*I_ERI_D2x_Dxz_Fy2z_S_b;
  abcd[697] = 2.0E0*I_ERI_Dxy_Dxz_Fy2z_S_b;
  abcd[698] = 2.0E0*I_ERI_Dxz_Dxz_Fy2z_S_b;
  abcd[699] = 2.0E0*I_ERI_D2y_Dxz_Fy2z_S_b;
  abcd[700] = 2.0E0*I_ERI_Dyz_Dxz_Fy2z_S_b;
  abcd[701] = 2.0E0*I_ERI_D2z_Dxz_Fy2z_S_b;
  abcd[702] = 2.0E0*I_ERI_D2x_D2x_F3z_S_b-1*I_ERI_D2x_S_F3z_S;
  abcd[703] = 2.0E0*I_ERI_Dxy_D2x_F3z_S_b-1*I_ERI_Dxy_S_F3z_S;
  abcd[704] = 2.0E0*I_ERI_Dxz_D2x_F3z_S_b-1*I_ERI_Dxz_S_F3z_S;
  abcd[705] = 2.0E0*I_ERI_D2y_D2x_F3z_S_b-1*I_ERI_D2y_S_F3z_S;
  abcd[706] = 2.0E0*I_ERI_Dyz_D2x_F3z_S_b-1*I_ERI_Dyz_S_F3z_S;
  abcd[707] = 2.0E0*I_ERI_D2z_D2x_F3z_S_b-1*I_ERI_D2z_S_F3z_S;
  abcd[708] = 2.0E0*I_ERI_D2x_Dxy_F3z_S_b;
  abcd[709] = 2.0E0*I_ERI_Dxy_Dxy_F3z_S_b;
  abcd[710] = 2.0E0*I_ERI_Dxz_Dxy_F3z_S_b;
  abcd[711] = 2.0E0*I_ERI_D2y_Dxy_F3z_S_b;
  abcd[712] = 2.0E0*I_ERI_Dyz_Dxy_F3z_S_b;
  abcd[713] = 2.0E0*I_ERI_D2z_Dxy_F3z_S_b;
  abcd[714] = 2.0E0*I_ERI_D2x_Dxz_F3z_S_b;
  abcd[715] = 2.0E0*I_ERI_Dxy_Dxz_F3z_S_b;
  abcd[716] = 2.0E0*I_ERI_Dxz_Dxz_F3z_S_b;
  abcd[717] = 2.0E0*I_ERI_D2y_Dxz_F3z_S_b;
  abcd[718] = 2.0E0*I_ERI_Dyz_Dxz_F3z_S_b;
  abcd[719] = 2.0E0*I_ERI_D2z_Dxz_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_b
   * RHS shell quartet name: SQ_ERI_D_S_F_S
   ************************************************************/
  abcd[720] = 2.0E0*I_ERI_D2x_Dxy_F3x_S_b;
  abcd[721] = 2.0E0*I_ERI_Dxy_Dxy_F3x_S_b;
  abcd[722] = 2.0E0*I_ERI_Dxz_Dxy_F3x_S_b;
  abcd[723] = 2.0E0*I_ERI_D2y_Dxy_F3x_S_b;
  abcd[724] = 2.0E0*I_ERI_Dyz_Dxy_F3x_S_b;
  abcd[725] = 2.0E0*I_ERI_D2z_Dxy_F3x_S_b;
  abcd[726] = 2.0E0*I_ERI_D2x_D2y_F3x_S_b-1*I_ERI_D2x_S_F3x_S;
  abcd[727] = 2.0E0*I_ERI_Dxy_D2y_F3x_S_b-1*I_ERI_Dxy_S_F3x_S;
  abcd[728] = 2.0E0*I_ERI_Dxz_D2y_F3x_S_b-1*I_ERI_Dxz_S_F3x_S;
  abcd[729] = 2.0E0*I_ERI_D2y_D2y_F3x_S_b-1*I_ERI_D2y_S_F3x_S;
  abcd[730] = 2.0E0*I_ERI_Dyz_D2y_F3x_S_b-1*I_ERI_Dyz_S_F3x_S;
  abcd[731] = 2.0E0*I_ERI_D2z_D2y_F3x_S_b-1*I_ERI_D2z_S_F3x_S;
  abcd[732] = 2.0E0*I_ERI_D2x_Dyz_F3x_S_b;
  abcd[733] = 2.0E0*I_ERI_Dxy_Dyz_F3x_S_b;
  abcd[734] = 2.0E0*I_ERI_Dxz_Dyz_F3x_S_b;
  abcd[735] = 2.0E0*I_ERI_D2y_Dyz_F3x_S_b;
  abcd[736] = 2.0E0*I_ERI_Dyz_Dyz_F3x_S_b;
  abcd[737] = 2.0E0*I_ERI_D2z_Dyz_F3x_S_b;
  abcd[738] = 2.0E0*I_ERI_D2x_Dxy_F2xy_S_b;
  abcd[739] = 2.0E0*I_ERI_Dxy_Dxy_F2xy_S_b;
  abcd[740] = 2.0E0*I_ERI_Dxz_Dxy_F2xy_S_b;
  abcd[741] = 2.0E0*I_ERI_D2y_Dxy_F2xy_S_b;
  abcd[742] = 2.0E0*I_ERI_Dyz_Dxy_F2xy_S_b;
  abcd[743] = 2.0E0*I_ERI_D2z_Dxy_F2xy_S_b;
  abcd[744] = 2.0E0*I_ERI_D2x_D2y_F2xy_S_b-1*I_ERI_D2x_S_F2xy_S;
  abcd[745] = 2.0E0*I_ERI_Dxy_D2y_F2xy_S_b-1*I_ERI_Dxy_S_F2xy_S;
  abcd[746] = 2.0E0*I_ERI_Dxz_D2y_F2xy_S_b-1*I_ERI_Dxz_S_F2xy_S;
  abcd[747] = 2.0E0*I_ERI_D2y_D2y_F2xy_S_b-1*I_ERI_D2y_S_F2xy_S;
  abcd[748] = 2.0E0*I_ERI_Dyz_D2y_F2xy_S_b-1*I_ERI_Dyz_S_F2xy_S;
  abcd[749] = 2.0E0*I_ERI_D2z_D2y_F2xy_S_b-1*I_ERI_D2z_S_F2xy_S;
  abcd[750] = 2.0E0*I_ERI_D2x_Dyz_F2xy_S_b;
  abcd[751] = 2.0E0*I_ERI_Dxy_Dyz_F2xy_S_b;
  abcd[752] = 2.0E0*I_ERI_Dxz_Dyz_F2xy_S_b;
  abcd[753] = 2.0E0*I_ERI_D2y_Dyz_F2xy_S_b;
  abcd[754] = 2.0E0*I_ERI_Dyz_Dyz_F2xy_S_b;
  abcd[755] = 2.0E0*I_ERI_D2z_Dyz_F2xy_S_b;
  abcd[756] = 2.0E0*I_ERI_D2x_Dxy_F2xz_S_b;
  abcd[757] = 2.0E0*I_ERI_Dxy_Dxy_F2xz_S_b;
  abcd[758] = 2.0E0*I_ERI_Dxz_Dxy_F2xz_S_b;
  abcd[759] = 2.0E0*I_ERI_D2y_Dxy_F2xz_S_b;
  abcd[760] = 2.0E0*I_ERI_Dyz_Dxy_F2xz_S_b;
  abcd[761] = 2.0E0*I_ERI_D2z_Dxy_F2xz_S_b;
  abcd[762] = 2.0E0*I_ERI_D2x_D2y_F2xz_S_b-1*I_ERI_D2x_S_F2xz_S;
  abcd[763] = 2.0E0*I_ERI_Dxy_D2y_F2xz_S_b-1*I_ERI_Dxy_S_F2xz_S;
  abcd[764] = 2.0E0*I_ERI_Dxz_D2y_F2xz_S_b-1*I_ERI_Dxz_S_F2xz_S;
  abcd[765] = 2.0E0*I_ERI_D2y_D2y_F2xz_S_b-1*I_ERI_D2y_S_F2xz_S;
  abcd[766] = 2.0E0*I_ERI_Dyz_D2y_F2xz_S_b-1*I_ERI_Dyz_S_F2xz_S;
  abcd[767] = 2.0E0*I_ERI_D2z_D2y_F2xz_S_b-1*I_ERI_D2z_S_F2xz_S;
  abcd[768] = 2.0E0*I_ERI_D2x_Dyz_F2xz_S_b;
  abcd[769] = 2.0E0*I_ERI_Dxy_Dyz_F2xz_S_b;
  abcd[770] = 2.0E0*I_ERI_Dxz_Dyz_F2xz_S_b;
  abcd[771] = 2.0E0*I_ERI_D2y_Dyz_F2xz_S_b;
  abcd[772] = 2.0E0*I_ERI_Dyz_Dyz_F2xz_S_b;
  abcd[773] = 2.0E0*I_ERI_D2z_Dyz_F2xz_S_b;
  abcd[774] = 2.0E0*I_ERI_D2x_Dxy_Fx2y_S_b;
  abcd[775] = 2.0E0*I_ERI_Dxy_Dxy_Fx2y_S_b;
  abcd[776] = 2.0E0*I_ERI_Dxz_Dxy_Fx2y_S_b;
  abcd[777] = 2.0E0*I_ERI_D2y_Dxy_Fx2y_S_b;
  abcd[778] = 2.0E0*I_ERI_Dyz_Dxy_Fx2y_S_b;
  abcd[779] = 2.0E0*I_ERI_D2z_Dxy_Fx2y_S_b;
  abcd[780] = 2.0E0*I_ERI_D2x_D2y_Fx2y_S_b-1*I_ERI_D2x_S_Fx2y_S;
  abcd[781] = 2.0E0*I_ERI_Dxy_D2y_Fx2y_S_b-1*I_ERI_Dxy_S_Fx2y_S;
  abcd[782] = 2.0E0*I_ERI_Dxz_D2y_Fx2y_S_b-1*I_ERI_Dxz_S_Fx2y_S;
  abcd[783] = 2.0E0*I_ERI_D2y_D2y_Fx2y_S_b-1*I_ERI_D2y_S_Fx2y_S;
  abcd[784] = 2.0E0*I_ERI_Dyz_D2y_Fx2y_S_b-1*I_ERI_Dyz_S_Fx2y_S;
  abcd[785] = 2.0E0*I_ERI_D2z_D2y_Fx2y_S_b-1*I_ERI_D2z_S_Fx2y_S;
  abcd[786] = 2.0E0*I_ERI_D2x_Dyz_Fx2y_S_b;
  abcd[787] = 2.0E0*I_ERI_Dxy_Dyz_Fx2y_S_b;
  abcd[788] = 2.0E0*I_ERI_Dxz_Dyz_Fx2y_S_b;
  abcd[789] = 2.0E0*I_ERI_D2y_Dyz_Fx2y_S_b;
  abcd[790] = 2.0E0*I_ERI_Dyz_Dyz_Fx2y_S_b;
  abcd[791] = 2.0E0*I_ERI_D2z_Dyz_Fx2y_S_b;
  abcd[792] = 2.0E0*I_ERI_D2x_Dxy_Fxyz_S_b;
  abcd[793] = 2.0E0*I_ERI_Dxy_Dxy_Fxyz_S_b;
  abcd[794] = 2.0E0*I_ERI_Dxz_Dxy_Fxyz_S_b;
  abcd[795] = 2.0E0*I_ERI_D2y_Dxy_Fxyz_S_b;
  abcd[796] = 2.0E0*I_ERI_Dyz_Dxy_Fxyz_S_b;
  abcd[797] = 2.0E0*I_ERI_D2z_Dxy_Fxyz_S_b;
  abcd[798] = 2.0E0*I_ERI_D2x_D2y_Fxyz_S_b-1*I_ERI_D2x_S_Fxyz_S;
  abcd[799] = 2.0E0*I_ERI_Dxy_D2y_Fxyz_S_b-1*I_ERI_Dxy_S_Fxyz_S;
  abcd[800] = 2.0E0*I_ERI_Dxz_D2y_Fxyz_S_b-1*I_ERI_Dxz_S_Fxyz_S;
  abcd[801] = 2.0E0*I_ERI_D2y_D2y_Fxyz_S_b-1*I_ERI_D2y_S_Fxyz_S;
  abcd[802] = 2.0E0*I_ERI_Dyz_D2y_Fxyz_S_b-1*I_ERI_Dyz_S_Fxyz_S;
  abcd[803] = 2.0E0*I_ERI_D2z_D2y_Fxyz_S_b-1*I_ERI_D2z_S_Fxyz_S;
  abcd[804] = 2.0E0*I_ERI_D2x_Dyz_Fxyz_S_b;
  abcd[805] = 2.0E0*I_ERI_Dxy_Dyz_Fxyz_S_b;
  abcd[806] = 2.0E0*I_ERI_Dxz_Dyz_Fxyz_S_b;
  abcd[807] = 2.0E0*I_ERI_D2y_Dyz_Fxyz_S_b;
  abcd[808] = 2.0E0*I_ERI_Dyz_Dyz_Fxyz_S_b;
  abcd[809] = 2.0E0*I_ERI_D2z_Dyz_Fxyz_S_b;
  abcd[810] = 2.0E0*I_ERI_D2x_Dxy_Fx2z_S_b;
  abcd[811] = 2.0E0*I_ERI_Dxy_Dxy_Fx2z_S_b;
  abcd[812] = 2.0E0*I_ERI_Dxz_Dxy_Fx2z_S_b;
  abcd[813] = 2.0E0*I_ERI_D2y_Dxy_Fx2z_S_b;
  abcd[814] = 2.0E0*I_ERI_Dyz_Dxy_Fx2z_S_b;
  abcd[815] = 2.0E0*I_ERI_D2z_Dxy_Fx2z_S_b;
  abcd[816] = 2.0E0*I_ERI_D2x_D2y_Fx2z_S_b-1*I_ERI_D2x_S_Fx2z_S;
  abcd[817] = 2.0E0*I_ERI_Dxy_D2y_Fx2z_S_b-1*I_ERI_Dxy_S_Fx2z_S;
  abcd[818] = 2.0E0*I_ERI_Dxz_D2y_Fx2z_S_b-1*I_ERI_Dxz_S_Fx2z_S;
  abcd[819] = 2.0E0*I_ERI_D2y_D2y_Fx2z_S_b-1*I_ERI_D2y_S_Fx2z_S;
  abcd[820] = 2.0E0*I_ERI_Dyz_D2y_Fx2z_S_b-1*I_ERI_Dyz_S_Fx2z_S;
  abcd[821] = 2.0E0*I_ERI_D2z_D2y_Fx2z_S_b-1*I_ERI_D2z_S_Fx2z_S;
  abcd[822] = 2.0E0*I_ERI_D2x_Dyz_Fx2z_S_b;
  abcd[823] = 2.0E0*I_ERI_Dxy_Dyz_Fx2z_S_b;
  abcd[824] = 2.0E0*I_ERI_Dxz_Dyz_Fx2z_S_b;
  abcd[825] = 2.0E0*I_ERI_D2y_Dyz_Fx2z_S_b;
  abcd[826] = 2.0E0*I_ERI_Dyz_Dyz_Fx2z_S_b;
  abcd[827] = 2.0E0*I_ERI_D2z_Dyz_Fx2z_S_b;
  abcd[828] = 2.0E0*I_ERI_D2x_Dxy_F3y_S_b;
  abcd[829] = 2.0E0*I_ERI_Dxy_Dxy_F3y_S_b;
  abcd[830] = 2.0E0*I_ERI_Dxz_Dxy_F3y_S_b;
  abcd[831] = 2.0E0*I_ERI_D2y_Dxy_F3y_S_b;
  abcd[832] = 2.0E0*I_ERI_Dyz_Dxy_F3y_S_b;
  abcd[833] = 2.0E0*I_ERI_D2z_Dxy_F3y_S_b;
  abcd[834] = 2.0E0*I_ERI_D2x_D2y_F3y_S_b-1*I_ERI_D2x_S_F3y_S;
  abcd[835] = 2.0E0*I_ERI_Dxy_D2y_F3y_S_b-1*I_ERI_Dxy_S_F3y_S;
  abcd[836] = 2.0E0*I_ERI_Dxz_D2y_F3y_S_b-1*I_ERI_Dxz_S_F3y_S;
  abcd[837] = 2.0E0*I_ERI_D2y_D2y_F3y_S_b-1*I_ERI_D2y_S_F3y_S;
  abcd[838] = 2.0E0*I_ERI_Dyz_D2y_F3y_S_b-1*I_ERI_Dyz_S_F3y_S;
  abcd[839] = 2.0E0*I_ERI_D2z_D2y_F3y_S_b-1*I_ERI_D2z_S_F3y_S;
  abcd[840] = 2.0E0*I_ERI_D2x_Dyz_F3y_S_b;
  abcd[841] = 2.0E0*I_ERI_Dxy_Dyz_F3y_S_b;
  abcd[842] = 2.0E0*I_ERI_Dxz_Dyz_F3y_S_b;
  abcd[843] = 2.0E0*I_ERI_D2y_Dyz_F3y_S_b;
  abcd[844] = 2.0E0*I_ERI_Dyz_Dyz_F3y_S_b;
  abcd[845] = 2.0E0*I_ERI_D2z_Dyz_F3y_S_b;
  abcd[846] = 2.0E0*I_ERI_D2x_Dxy_F2yz_S_b;
  abcd[847] = 2.0E0*I_ERI_Dxy_Dxy_F2yz_S_b;
  abcd[848] = 2.0E0*I_ERI_Dxz_Dxy_F2yz_S_b;
  abcd[849] = 2.0E0*I_ERI_D2y_Dxy_F2yz_S_b;
  abcd[850] = 2.0E0*I_ERI_Dyz_Dxy_F2yz_S_b;
  abcd[851] = 2.0E0*I_ERI_D2z_Dxy_F2yz_S_b;
  abcd[852] = 2.0E0*I_ERI_D2x_D2y_F2yz_S_b-1*I_ERI_D2x_S_F2yz_S;
  abcd[853] = 2.0E0*I_ERI_Dxy_D2y_F2yz_S_b-1*I_ERI_Dxy_S_F2yz_S;
  abcd[854] = 2.0E0*I_ERI_Dxz_D2y_F2yz_S_b-1*I_ERI_Dxz_S_F2yz_S;
  abcd[855] = 2.0E0*I_ERI_D2y_D2y_F2yz_S_b-1*I_ERI_D2y_S_F2yz_S;
  abcd[856] = 2.0E0*I_ERI_Dyz_D2y_F2yz_S_b-1*I_ERI_Dyz_S_F2yz_S;
  abcd[857] = 2.0E0*I_ERI_D2z_D2y_F2yz_S_b-1*I_ERI_D2z_S_F2yz_S;
  abcd[858] = 2.0E0*I_ERI_D2x_Dyz_F2yz_S_b;
  abcd[859] = 2.0E0*I_ERI_Dxy_Dyz_F2yz_S_b;
  abcd[860] = 2.0E0*I_ERI_Dxz_Dyz_F2yz_S_b;
  abcd[861] = 2.0E0*I_ERI_D2y_Dyz_F2yz_S_b;
  abcd[862] = 2.0E0*I_ERI_Dyz_Dyz_F2yz_S_b;
  abcd[863] = 2.0E0*I_ERI_D2z_Dyz_F2yz_S_b;
  abcd[864] = 2.0E0*I_ERI_D2x_Dxy_Fy2z_S_b;
  abcd[865] = 2.0E0*I_ERI_Dxy_Dxy_Fy2z_S_b;
  abcd[866] = 2.0E0*I_ERI_Dxz_Dxy_Fy2z_S_b;
  abcd[867] = 2.0E0*I_ERI_D2y_Dxy_Fy2z_S_b;
  abcd[868] = 2.0E0*I_ERI_Dyz_Dxy_Fy2z_S_b;
  abcd[869] = 2.0E0*I_ERI_D2z_Dxy_Fy2z_S_b;
  abcd[870] = 2.0E0*I_ERI_D2x_D2y_Fy2z_S_b-1*I_ERI_D2x_S_Fy2z_S;
  abcd[871] = 2.0E0*I_ERI_Dxy_D2y_Fy2z_S_b-1*I_ERI_Dxy_S_Fy2z_S;
  abcd[872] = 2.0E0*I_ERI_Dxz_D2y_Fy2z_S_b-1*I_ERI_Dxz_S_Fy2z_S;
  abcd[873] = 2.0E0*I_ERI_D2y_D2y_Fy2z_S_b-1*I_ERI_D2y_S_Fy2z_S;
  abcd[874] = 2.0E0*I_ERI_Dyz_D2y_Fy2z_S_b-1*I_ERI_Dyz_S_Fy2z_S;
  abcd[875] = 2.0E0*I_ERI_D2z_D2y_Fy2z_S_b-1*I_ERI_D2z_S_Fy2z_S;
  abcd[876] = 2.0E0*I_ERI_D2x_Dyz_Fy2z_S_b;
  abcd[877] = 2.0E0*I_ERI_Dxy_Dyz_Fy2z_S_b;
  abcd[878] = 2.0E0*I_ERI_Dxz_Dyz_Fy2z_S_b;
  abcd[879] = 2.0E0*I_ERI_D2y_Dyz_Fy2z_S_b;
  abcd[880] = 2.0E0*I_ERI_Dyz_Dyz_Fy2z_S_b;
  abcd[881] = 2.0E0*I_ERI_D2z_Dyz_Fy2z_S_b;
  abcd[882] = 2.0E0*I_ERI_D2x_Dxy_F3z_S_b;
  abcd[883] = 2.0E0*I_ERI_Dxy_Dxy_F3z_S_b;
  abcd[884] = 2.0E0*I_ERI_Dxz_Dxy_F3z_S_b;
  abcd[885] = 2.0E0*I_ERI_D2y_Dxy_F3z_S_b;
  abcd[886] = 2.0E0*I_ERI_Dyz_Dxy_F3z_S_b;
  abcd[887] = 2.0E0*I_ERI_D2z_Dxy_F3z_S_b;
  abcd[888] = 2.0E0*I_ERI_D2x_D2y_F3z_S_b-1*I_ERI_D2x_S_F3z_S;
  abcd[889] = 2.0E0*I_ERI_Dxy_D2y_F3z_S_b-1*I_ERI_Dxy_S_F3z_S;
  abcd[890] = 2.0E0*I_ERI_Dxz_D2y_F3z_S_b-1*I_ERI_Dxz_S_F3z_S;
  abcd[891] = 2.0E0*I_ERI_D2y_D2y_F3z_S_b-1*I_ERI_D2y_S_F3z_S;
  abcd[892] = 2.0E0*I_ERI_Dyz_D2y_F3z_S_b-1*I_ERI_Dyz_S_F3z_S;
  abcd[893] = 2.0E0*I_ERI_D2z_D2y_F3z_S_b-1*I_ERI_D2z_S_F3z_S;
  abcd[894] = 2.0E0*I_ERI_D2x_Dyz_F3z_S_b;
  abcd[895] = 2.0E0*I_ERI_Dxy_Dyz_F3z_S_b;
  abcd[896] = 2.0E0*I_ERI_Dxz_Dyz_F3z_S_b;
  abcd[897] = 2.0E0*I_ERI_D2y_Dyz_F3z_S_b;
  abcd[898] = 2.0E0*I_ERI_Dyz_Dyz_F3z_S_b;
  abcd[899] = 2.0E0*I_ERI_D2z_Dyz_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_b
   * RHS shell quartet name: SQ_ERI_D_S_F_S
   ************************************************************/
  abcd[900] = 2.0E0*I_ERI_D2x_Dxz_F3x_S_b;
  abcd[901] = 2.0E0*I_ERI_Dxy_Dxz_F3x_S_b;
  abcd[902] = 2.0E0*I_ERI_Dxz_Dxz_F3x_S_b;
  abcd[903] = 2.0E0*I_ERI_D2y_Dxz_F3x_S_b;
  abcd[904] = 2.0E0*I_ERI_Dyz_Dxz_F3x_S_b;
  abcd[905] = 2.0E0*I_ERI_D2z_Dxz_F3x_S_b;
  abcd[906] = 2.0E0*I_ERI_D2x_Dyz_F3x_S_b;
  abcd[907] = 2.0E0*I_ERI_Dxy_Dyz_F3x_S_b;
  abcd[908] = 2.0E0*I_ERI_Dxz_Dyz_F3x_S_b;
  abcd[909] = 2.0E0*I_ERI_D2y_Dyz_F3x_S_b;
  abcd[910] = 2.0E0*I_ERI_Dyz_Dyz_F3x_S_b;
  abcd[911] = 2.0E0*I_ERI_D2z_Dyz_F3x_S_b;
  abcd[912] = 2.0E0*I_ERI_D2x_D2z_F3x_S_b-1*I_ERI_D2x_S_F3x_S;
  abcd[913] = 2.0E0*I_ERI_Dxy_D2z_F3x_S_b-1*I_ERI_Dxy_S_F3x_S;
  abcd[914] = 2.0E0*I_ERI_Dxz_D2z_F3x_S_b-1*I_ERI_Dxz_S_F3x_S;
  abcd[915] = 2.0E0*I_ERI_D2y_D2z_F3x_S_b-1*I_ERI_D2y_S_F3x_S;
  abcd[916] = 2.0E0*I_ERI_Dyz_D2z_F3x_S_b-1*I_ERI_Dyz_S_F3x_S;
  abcd[917] = 2.0E0*I_ERI_D2z_D2z_F3x_S_b-1*I_ERI_D2z_S_F3x_S;
  abcd[918] = 2.0E0*I_ERI_D2x_Dxz_F2xy_S_b;
  abcd[919] = 2.0E0*I_ERI_Dxy_Dxz_F2xy_S_b;
  abcd[920] = 2.0E0*I_ERI_Dxz_Dxz_F2xy_S_b;
  abcd[921] = 2.0E0*I_ERI_D2y_Dxz_F2xy_S_b;
  abcd[922] = 2.0E0*I_ERI_Dyz_Dxz_F2xy_S_b;
  abcd[923] = 2.0E0*I_ERI_D2z_Dxz_F2xy_S_b;
  abcd[924] = 2.0E0*I_ERI_D2x_Dyz_F2xy_S_b;
  abcd[925] = 2.0E0*I_ERI_Dxy_Dyz_F2xy_S_b;
  abcd[926] = 2.0E0*I_ERI_Dxz_Dyz_F2xy_S_b;
  abcd[927] = 2.0E0*I_ERI_D2y_Dyz_F2xy_S_b;
  abcd[928] = 2.0E0*I_ERI_Dyz_Dyz_F2xy_S_b;
  abcd[929] = 2.0E0*I_ERI_D2z_Dyz_F2xy_S_b;
  abcd[930] = 2.0E0*I_ERI_D2x_D2z_F2xy_S_b-1*I_ERI_D2x_S_F2xy_S;
  abcd[931] = 2.0E0*I_ERI_Dxy_D2z_F2xy_S_b-1*I_ERI_Dxy_S_F2xy_S;
  abcd[932] = 2.0E0*I_ERI_Dxz_D2z_F2xy_S_b-1*I_ERI_Dxz_S_F2xy_S;
  abcd[933] = 2.0E0*I_ERI_D2y_D2z_F2xy_S_b-1*I_ERI_D2y_S_F2xy_S;
  abcd[934] = 2.0E0*I_ERI_Dyz_D2z_F2xy_S_b-1*I_ERI_Dyz_S_F2xy_S;
  abcd[935] = 2.0E0*I_ERI_D2z_D2z_F2xy_S_b-1*I_ERI_D2z_S_F2xy_S;
  abcd[936] = 2.0E0*I_ERI_D2x_Dxz_F2xz_S_b;
  abcd[937] = 2.0E0*I_ERI_Dxy_Dxz_F2xz_S_b;
  abcd[938] = 2.0E0*I_ERI_Dxz_Dxz_F2xz_S_b;
  abcd[939] = 2.0E0*I_ERI_D2y_Dxz_F2xz_S_b;
  abcd[940] = 2.0E0*I_ERI_Dyz_Dxz_F2xz_S_b;
  abcd[941] = 2.0E0*I_ERI_D2z_Dxz_F2xz_S_b;
  abcd[942] = 2.0E0*I_ERI_D2x_Dyz_F2xz_S_b;
  abcd[943] = 2.0E0*I_ERI_Dxy_Dyz_F2xz_S_b;
  abcd[944] = 2.0E0*I_ERI_Dxz_Dyz_F2xz_S_b;
  abcd[945] = 2.0E0*I_ERI_D2y_Dyz_F2xz_S_b;
  abcd[946] = 2.0E0*I_ERI_Dyz_Dyz_F2xz_S_b;
  abcd[947] = 2.0E0*I_ERI_D2z_Dyz_F2xz_S_b;
  abcd[948] = 2.0E0*I_ERI_D2x_D2z_F2xz_S_b-1*I_ERI_D2x_S_F2xz_S;
  abcd[949] = 2.0E0*I_ERI_Dxy_D2z_F2xz_S_b-1*I_ERI_Dxy_S_F2xz_S;
  abcd[950] = 2.0E0*I_ERI_Dxz_D2z_F2xz_S_b-1*I_ERI_Dxz_S_F2xz_S;
  abcd[951] = 2.0E0*I_ERI_D2y_D2z_F2xz_S_b-1*I_ERI_D2y_S_F2xz_S;
  abcd[952] = 2.0E0*I_ERI_Dyz_D2z_F2xz_S_b-1*I_ERI_Dyz_S_F2xz_S;
  abcd[953] = 2.0E0*I_ERI_D2z_D2z_F2xz_S_b-1*I_ERI_D2z_S_F2xz_S;
  abcd[954] = 2.0E0*I_ERI_D2x_Dxz_Fx2y_S_b;
  abcd[955] = 2.0E0*I_ERI_Dxy_Dxz_Fx2y_S_b;
  abcd[956] = 2.0E0*I_ERI_Dxz_Dxz_Fx2y_S_b;
  abcd[957] = 2.0E0*I_ERI_D2y_Dxz_Fx2y_S_b;
  abcd[958] = 2.0E0*I_ERI_Dyz_Dxz_Fx2y_S_b;
  abcd[959] = 2.0E0*I_ERI_D2z_Dxz_Fx2y_S_b;
  abcd[960] = 2.0E0*I_ERI_D2x_Dyz_Fx2y_S_b;
  abcd[961] = 2.0E0*I_ERI_Dxy_Dyz_Fx2y_S_b;
  abcd[962] = 2.0E0*I_ERI_Dxz_Dyz_Fx2y_S_b;
  abcd[963] = 2.0E0*I_ERI_D2y_Dyz_Fx2y_S_b;
  abcd[964] = 2.0E0*I_ERI_Dyz_Dyz_Fx2y_S_b;
  abcd[965] = 2.0E0*I_ERI_D2z_Dyz_Fx2y_S_b;
  abcd[966] = 2.0E0*I_ERI_D2x_D2z_Fx2y_S_b-1*I_ERI_D2x_S_Fx2y_S;
  abcd[967] = 2.0E0*I_ERI_Dxy_D2z_Fx2y_S_b-1*I_ERI_Dxy_S_Fx2y_S;
  abcd[968] = 2.0E0*I_ERI_Dxz_D2z_Fx2y_S_b-1*I_ERI_Dxz_S_Fx2y_S;
  abcd[969] = 2.0E0*I_ERI_D2y_D2z_Fx2y_S_b-1*I_ERI_D2y_S_Fx2y_S;
  abcd[970] = 2.0E0*I_ERI_Dyz_D2z_Fx2y_S_b-1*I_ERI_Dyz_S_Fx2y_S;
  abcd[971] = 2.0E0*I_ERI_D2z_D2z_Fx2y_S_b-1*I_ERI_D2z_S_Fx2y_S;
  abcd[972] = 2.0E0*I_ERI_D2x_Dxz_Fxyz_S_b;
  abcd[973] = 2.0E0*I_ERI_Dxy_Dxz_Fxyz_S_b;
  abcd[974] = 2.0E0*I_ERI_Dxz_Dxz_Fxyz_S_b;
  abcd[975] = 2.0E0*I_ERI_D2y_Dxz_Fxyz_S_b;
  abcd[976] = 2.0E0*I_ERI_Dyz_Dxz_Fxyz_S_b;
  abcd[977] = 2.0E0*I_ERI_D2z_Dxz_Fxyz_S_b;
  abcd[978] = 2.0E0*I_ERI_D2x_Dyz_Fxyz_S_b;
  abcd[979] = 2.0E0*I_ERI_Dxy_Dyz_Fxyz_S_b;
  abcd[980] = 2.0E0*I_ERI_Dxz_Dyz_Fxyz_S_b;
  abcd[981] = 2.0E0*I_ERI_D2y_Dyz_Fxyz_S_b;
  abcd[982] = 2.0E0*I_ERI_Dyz_Dyz_Fxyz_S_b;
  abcd[983] = 2.0E0*I_ERI_D2z_Dyz_Fxyz_S_b;
  abcd[984] = 2.0E0*I_ERI_D2x_D2z_Fxyz_S_b-1*I_ERI_D2x_S_Fxyz_S;
  abcd[985] = 2.0E0*I_ERI_Dxy_D2z_Fxyz_S_b-1*I_ERI_Dxy_S_Fxyz_S;
  abcd[986] = 2.0E0*I_ERI_Dxz_D2z_Fxyz_S_b-1*I_ERI_Dxz_S_Fxyz_S;
  abcd[987] = 2.0E0*I_ERI_D2y_D2z_Fxyz_S_b-1*I_ERI_D2y_S_Fxyz_S;
  abcd[988] = 2.0E0*I_ERI_Dyz_D2z_Fxyz_S_b-1*I_ERI_Dyz_S_Fxyz_S;
  abcd[989] = 2.0E0*I_ERI_D2z_D2z_Fxyz_S_b-1*I_ERI_D2z_S_Fxyz_S;
  abcd[990] = 2.0E0*I_ERI_D2x_Dxz_Fx2z_S_b;
  abcd[991] = 2.0E0*I_ERI_Dxy_Dxz_Fx2z_S_b;
  abcd[992] = 2.0E0*I_ERI_Dxz_Dxz_Fx2z_S_b;
  abcd[993] = 2.0E0*I_ERI_D2y_Dxz_Fx2z_S_b;
  abcd[994] = 2.0E0*I_ERI_Dyz_Dxz_Fx2z_S_b;
  abcd[995] = 2.0E0*I_ERI_D2z_Dxz_Fx2z_S_b;
  abcd[996] = 2.0E0*I_ERI_D2x_Dyz_Fx2z_S_b;
  abcd[997] = 2.0E0*I_ERI_Dxy_Dyz_Fx2z_S_b;
  abcd[998] = 2.0E0*I_ERI_Dxz_Dyz_Fx2z_S_b;
  abcd[999] = 2.0E0*I_ERI_D2y_Dyz_Fx2z_S_b;
  abcd[1000] = 2.0E0*I_ERI_Dyz_Dyz_Fx2z_S_b;
  abcd[1001] = 2.0E0*I_ERI_D2z_Dyz_Fx2z_S_b;
  abcd[1002] = 2.0E0*I_ERI_D2x_D2z_Fx2z_S_b-1*I_ERI_D2x_S_Fx2z_S;
  abcd[1003] = 2.0E0*I_ERI_Dxy_D2z_Fx2z_S_b-1*I_ERI_Dxy_S_Fx2z_S;
  abcd[1004] = 2.0E0*I_ERI_Dxz_D2z_Fx2z_S_b-1*I_ERI_Dxz_S_Fx2z_S;
  abcd[1005] = 2.0E0*I_ERI_D2y_D2z_Fx2z_S_b-1*I_ERI_D2y_S_Fx2z_S;
  abcd[1006] = 2.0E0*I_ERI_Dyz_D2z_Fx2z_S_b-1*I_ERI_Dyz_S_Fx2z_S;
  abcd[1007] = 2.0E0*I_ERI_D2z_D2z_Fx2z_S_b-1*I_ERI_D2z_S_Fx2z_S;
  abcd[1008] = 2.0E0*I_ERI_D2x_Dxz_F3y_S_b;
  abcd[1009] = 2.0E0*I_ERI_Dxy_Dxz_F3y_S_b;
  abcd[1010] = 2.0E0*I_ERI_Dxz_Dxz_F3y_S_b;
  abcd[1011] = 2.0E0*I_ERI_D2y_Dxz_F3y_S_b;
  abcd[1012] = 2.0E0*I_ERI_Dyz_Dxz_F3y_S_b;
  abcd[1013] = 2.0E0*I_ERI_D2z_Dxz_F3y_S_b;
  abcd[1014] = 2.0E0*I_ERI_D2x_Dyz_F3y_S_b;
  abcd[1015] = 2.0E0*I_ERI_Dxy_Dyz_F3y_S_b;
  abcd[1016] = 2.0E0*I_ERI_Dxz_Dyz_F3y_S_b;
  abcd[1017] = 2.0E0*I_ERI_D2y_Dyz_F3y_S_b;
  abcd[1018] = 2.0E0*I_ERI_Dyz_Dyz_F3y_S_b;
  abcd[1019] = 2.0E0*I_ERI_D2z_Dyz_F3y_S_b;
  abcd[1020] = 2.0E0*I_ERI_D2x_D2z_F3y_S_b-1*I_ERI_D2x_S_F3y_S;
  abcd[1021] = 2.0E0*I_ERI_Dxy_D2z_F3y_S_b-1*I_ERI_Dxy_S_F3y_S;
  abcd[1022] = 2.0E0*I_ERI_Dxz_D2z_F3y_S_b-1*I_ERI_Dxz_S_F3y_S;
  abcd[1023] = 2.0E0*I_ERI_D2y_D2z_F3y_S_b-1*I_ERI_D2y_S_F3y_S;
  abcd[1024] = 2.0E0*I_ERI_Dyz_D2z_F3y_S_b-1*I_ERI_Dyz_S_F3y_S;
  abcd[1025] = 2.0E0*I_ERI_D2z_D2z_F3y_S_b-1*I_ERI_D2z_S_F3y_S;
  abcd[1026] = 2.0E0*I_ERI_D2x_Dxz_F2yz_S_b;
  abcd[1027] = 2.0E0*I_ERI_Dxy_Dxz_F2yz_S_b;
  abcd[1028] = 2.0E0*I_ERI_Dxz_Dxz_F2yz_S_b;
  abcd[1029] = 2.0E0*I_ERI_D2y_Dxz_F2yz_S_b;
  abcd[1030] = 2.0E0*I_ERI_Dyz_Dxz_F2yz_S_b;
  abcd[1031] = 2.0E0*I_ERI_D2z_Dxz_F2yz_S_b;
  abcd[1032] = 2.0E0*I_ERI_D2x_Dyz_F2yz_S_b;
  abcd[1033] = 2.0E0*I_ERI_Dxy_Dyz_F2yz_S_b;
  abcd[1034] = 2.0E0*I_ERI_Dxz_Dyz_F2yz_S_b;
  abcd[1035] = 2.0E0*I_ERI_D2y_Dyz_F2yz_S_b;
  abcd[1036] = 2.0E0*I_ERI_Dyz_Dyz_F2yz_S_b;
  abcd[1037] = 2.0E0*I_ERI_D2z_Dyz_F2yz_S_b;
  abcd[1038] = 2.0E0*I_ERI_D2x_D2z_F2yz_S_b-1*I_ERI_D2x_S_F2yz_S;
  abcd[1039] = 2.0E0*I_ERI_Dxy_D2z_F2yz_S_b-1*I_ERI_Dxy_S_F2yz_S;
  abcd[1040] = 2.0E0*I_ERI_Dxz_D2z_F2yz_S_b-1*I_ERI_Dxz_S_F2yz_S;
  abcd[1041] = 2.0E0*I_ERI_D2y_D2z_F2yz_S_b-1*I_ERI_D2y_S_F2yz_S;
  abcd[1042] = 2.0E0*I_ERI_Dyz_D2z_F2yz_S_b-1*I_ERI_Dyz_S_F2yz_S;
  abcd[1043] = 2.0E0*I_ERI_D2z_D2z_F2yz_S_b-1*I_ERI_D2z_S_F2yz_S;
  abcd[1044] = 2.0E0*I_ERI_D2x_Dxz_Fy2z_S_b;
  abcd[1045] = 2.0E0*I_ERI_Dxy_Dxz_Fy2z_S_b;
  abcd[1046] = 2.0E0*I_ERI_Dxz_Dxz_Fy2z_S_b;
  abcd[1047] = 2.0E0*I_ERI_D2y_Dxz_Fy2z_S_b;
  abcd[1048] = 2.0E0*I_ERI_Dyz_Dxz_Fy2z_S_b;
  abcd[1049] = 2.0E0*I_ERI_D2z_Dxz_Fy2z_S_b;
  abcd[1050] = 2.0E0*I_ERI_D2x_Dyz_Fy2z_S_b;
  abcd[1051] = 2.0E0*I_ERI_Dxy_Dyz_Fy2z_S_b;
  abcd[1052] = 2.0E0*I_ERI_Dxz_Dyz_Fy2z_S_b;
  abcd[1053] = 2.0E0*I_ERI_D2y_Dyz_Fy2z_S_b;
  abcd[1054] = 2.0E0*I_ERI_Dyz_Dyz_Fy2z_S_b;
  abcd[1055] = 2.0E0*I_ERI_D2z_Dyz_Fy2z_S_b;
  abcd[1056] = 2.0E0*I_ERI_D2x_D2z_Fy2z_S_b-1*I_ERI_D2x_S_Fy2z_S;
  abcd[1057] = 2.0E0*I_ERI_Dxy_D2z_Fy2z_S_b-1*I_ERI_Dxy_S_Fy2z_S;
  abcd[1058] = 2.0E0*I_ERI_Dxz_D2z_Fy2z_S_b-1*I_ERI_Dxz_S_Fy2z_S;
  abcd[1059] = 2.0E0*I_ERI_D2y_D2z_Fy2z_S_b-1*I_ERI_D2y_S_Fy2z_S;
  abcd[1060] = 2.0E0*I_ERI_Dyz_D2z_Fy2z_S_b-1*I_ERI_Dyz_S_Fy2z_S;
  abcd[1061] = 2.0E0*I_ERI_D2z_D2z_Fy2z_S_b-1*I_ERI_D2z_S_Fy2z_S;
  abcd[1062] = 2.0E0*I_ERI_D2x_Dxz_F3z_S_b;
  abcd[1063] = 2.0E0*I_ERI_Dxy_Dxz_F3z_S_b;
  abcd[1064] = 2.0E0*I_ERI_Dxz_Dxz_F3z_S_b;
  abcd[1065] = 2.0E0*I_ERI_D2y_Dxz_F3z_S_b;
  abcd[1066] = 2.0E0*I_ERI_Dyz_Dxz_F3z_S_b;
  abcd[1067] = 2.0E0*I_ERI_D2z_Dxz_F3z_S_b;
  abcd[1068] = 2.0E0*I_ERI_D2x_Dyz_F3z_S_b;
  abcd[1069] = 2.0E0*I_ERI_Dxy_Dyz_F3z_S_b;
  abcd[1070] = 2.0E0*I_ERI_Dxz_Dyz_F3z_S_b;
  abcd[1071] = 2.0E0*I_ERI_D2y_Dyz_F3z_S_b;
  abcd[1072] = 2.0E0*I_ERI_Dyz_Dyz_F3z_S_b;
  abcd[1073] = 2.0E0*I_ERI_D2z_Dyz_F3z_S_b;
  abcd[1074] = 2.0E0*I_ERI_D2x_D2z_F3z_S_b-1*I_ERI_D2x_S_F3z_S;
  abcd[1075] = 2.0E0*I_ERI_Dxy_D2z_F3z_S_b-1*I_ERI_Dxy_S_F3z_S;
  abcd[1076] = 2.0E0*I_ERI_Dxz_D2z_F3z_S_b-1*I_ERI_Dxz_S_F3z_S;
  abcd[1077] = 2.0E0*I_ERI_D2y_D2z_F3z_S_b-1*I_ERI_D2y_S_F3z_S;
  abcd[1078] = 2.0E0*I_ERI_Dyz_D2z_F3z_S_b-1*I_ERI_Dyz_S_F3z_S;
  abcd[1079] = 2.0E0*I_ERI_D2z_D2z_F3z_S_b-1*I_ERI_D2z_S_F3z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_ddx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_P_d
   ************************************************************/
  abcd[1080] = 2.0E0*I_ERI_D2x_Px_F3x_Px_d;
  abcd[1081] = 2.0E0*I_ERI_Dxy_Px_F3x_Px_d;
  abcd[1082] = 2.0E0*I_ERI_Dxz_Px_F3x_Px_d;
  abcd[1083] = 2.0E0*I_ERI_D2y_Px_F3x_Px_d;
  abcd[1084] = 2.0E0*I_ERI_Dyz_Px_F3x_Px_d;
  abcd[1085] = 2.0E0*I_ERI_D2z_Px_F3x_Px_d;
  abcd[1086] = 2.0E0*I_ERI_D2x_Py_F3x_Px_d;
  abcd[1087] = 2.0E0*I_ERI_Dxy_Py_F3x_Px_d;
  abcd[1088] = 2.0E0*I_ERI_Dxz_Py_F3x_Px_d;
  abcd[1089] = 2.0E0*I_ERI_D2y_Py_F3x_Px_d;
  abcd[1090] = 2.0E0*I_ERI_Dyz_Py_F3x_Px_d;
  abcd[1091] = 2.0E0*I_ERI_D2z_Py_F3x_Px_d;
  abcd[1092] = 2.0E0*I_ERI_D2x_Pz_F3x_Px_d;
  abcd[1093] = 2.0E0*I_ERI_Dxy_Pz_F3x_Px_d;
  abcd[1094] = 2.0E0*I_ERI_Dxz_Pz_F3x_Px_d;
  abcd[1095] = 2.0E0*I_ERI_D2y_Pz_F3x_Px_d;
  abcd[1096] = 2.0E0*I_ERI_Dyz_Pz_F3x_Px_d;
  abcd[1097] = 2.0E0*I_ERI_D2z_Pz_F3x_Px_d;
  abcd[1098] = 2.0E0*I_ERI_D2x_Px_F2xy_Px_d;
  abcd[1099] = 2.0E0*I_ERI_Dxy_Px_F2xy_Px_d;
  abcd[1100] = 2.0E0*I_ERI_Dxz_Px_F2xy_Px_d;
  abcd[1101] = 2.0E0*I_ERI_D2y_Px_F2xy_Px_d;
  abcd[1102] = 2.0E0*I_ERI_Dyz_Px_F2xy_Px_d;
  abcd[1103] = 2.0E0*I_ERI_D2z_Px_F2xy_Px_d;
  abcd[1104] = 2.0E0*I_ERI_D2x_Py_F2xy_Px_d;
  abcd[1105] = 2.0E0*I_ERI_Dxy_Py_F2xy_Px_d;
  abcd[1106] = 2.0E0*I_ERI_Dxz_Py_F2xy_Px_d;
  abcd[1107] = 2.0E0*I_ERI_D2y_Py_F2xy_Px_d;
  abcd[1108] = 2.0E0*I_ERI_Dyz_Py_F2xy_Px_d;
  abcd[1109] = 2.0E0*I_ERI_D2z_Py_F2xy_Px_d;
  abcd[1110] = 2.0E0*I_ERI_D2x_Pz_F2xy_Px_d;
  abcd[1111] = 2.0E0*I_ERI_Dxy_Pz_F2xy_Px_d;
  abcd[1112] = 2.0E0*I_ERI_Dxz_Pz_F2xy_Px_d;
  abcd[1113] = 2.0E0*I_ERI_D2y_Pz_F2xy_Px_d;
  abcd[1114] = 2.0E0*I_ERI_Dyz_Pz_F2xy_Px_d;
  abcd[1115] = 2.0E0*I_ERI_D2z_Pz_F2xy_Px_d;
  abcd[1116] = 2.0E0*I_ERI_D2x_Px_F2xz_Px_d;
  abcd[1117] = 2.0E0*I_ERI_Dxy_Px_F2xz_Px_d;
  abcd[1118] = 2.0E0*I_ERI_Dxz_Px_F2xz_Px_d;
  abcd[1119] = 2.0E0*I_ERI_D2y_Px_F2xz_Px_d;
  abcd[1120] = 2.0E0*I_ERI_Dyz_Px_F2xz_Px_d;
  abcd[1121] = 2.0E0*I_ERI_D2z_Px_F2xz_Px_d;
  abcd[1122] = 2.0E0*I_ERI_D2x_Py_F2xz_Px_d;
  abcd[1123] = 2.0E0*I_ERI_Dxy_Py_F2xz_Px_d;
  abcd[1124] = 2.0E0*I_ERI_Dxz_Py_F2xz_Px_d;
  abcd[1125] = 2.0E0*I_ERI_D2y_Py_F2xz_Px_d;
  abcd[1126] = 2.0E0*I_ERI_Dyz_Py_F2xz_Px_d;
  abcd[1127] = 2.0E0*I_ERI_D2z_Py_F2xz_Px_d;
  abcd[1128] = 2.0E0*I_ERI_D2x_Pz_F2xz_Px_d;
  abcd[1129] = 2.0E0*I_ERI_Dxy_Pz_F2xz_Px_d;
  abcd[1130] = 2.0E0*I_ERI_Dxz_Pz_F2xz_Px_d;
  abcd[1131] = 2.0E0*I_ERI_D2y_Pz_F2xz_Px_d;
  abcd[1132] = 2.0E0*I_ERI_Dyz_Pz_F2xz_Px_d;
  abcd[1133] = 2.0E0*I_ERI_D2z_Pz_F2xz_Px_d;
  abcd[1134] = 2.0E0*I_ERI_D2x_Px_Fx2y_Px_d;
  abcd[1135] = 2.0E0*I_ERI_Dxy_Px_Fx2y_Px_d;
  abcd[1136] = 2.0E0*I_ERI_Dxz_Px_Fx2y_Px_d;
  abcd[1137] = 2.0E0*I_ERI_D2y_Px_Fx2y_Px_d;
  abcd[1138] = 2.0E0*I_ERI_Dyz_Px_Fx2y_Px_d;
  abcd[1139] = 2.0E0*I_ERI_D2z_Px_Fx2y_Px_d;
  abcd[1140] = 2.0E0*I_ERI_D2x_Py_Fx2y_Px_d;
  abcd[1141] = 2.0E0*I_ERI_Dxy_Py_Fx2y_Px_d;
  abcd[1142] = 2.0E0*I_ERI_Dxz_Py_Fx2y_Px_d;
  abcd[1143] = 2.0E0*I_ERI_D2y_Py_Fx2y_Px_d;
  abcd[1144] = 2.0E0*I_ERI_Dyz_Py_Fx2y_Px_d;
  abcd[1145] = 2.0E0*I_ERI_D2z_Py_Fx2y_Px_d;
  abcd[1146] = 2.0E0*I_ERI_D2x_Pz_Fx2y_Px_d;
  abcd[1147] = 2.0E0*I_ERI_Dxy_Pz_Fx2y_Px_d;
  abcd[1148] = 2.0E0*I_ERI_Dxz_Pz_Fx2y_Px_d;
  abcd[1149] = 2.0E0*I_ERI_D2y_Pz_Fx2y_Px_d;
  abcd[1150] = 2.0E0*I_ERI_Dyz_Pz_Fx2y_Px_d;
  abcd[1151] = 2.0E0*I_ERI_D2z_Pz_Fx2y_Px_d;
  abcd[1152] = 2.0E0*I_ERI_D2x_Px_Fxyz_Px_d;
  abcd[1153] = 2.0E0*I_ERI_Dxy_Px_Fxyz_Px_d;
  abcd[1154] = 2.0E0*I_ERI_Dxz_Px_Fxyz_Px_d;
  abcd[1155] = 2.0E0*I_ERI_D2y_Px_Fxyz_Px_d;
  abcd[1156] = 2.0E0*I_ERI_Dyz_Px_Fxyz_Px_d;
  abcd[1157] = 2.0E0*I_ERI_D2z_Px_Fxyz_Px_d;
  abcd[1158] = 2.0E0*I_ERI_D2x_Py_Fxyz_Px_d;
  abcd[1159] = 2.0E0*I_ERI_Dxy_Py_Fxyz_Px_d;
  abcd[1160] = 2.0E0*I_ERI_Dxz_Py_Fxyz_Px_d;
  abcd[1161] = 2.0E0*I_ERI_D2y_Py_Fxyz_Px_d;
  abcd[1162] = 2.0E0*I_ERI_Dyz_Py_Fxyz_Px_d;
  abcd[1163] = 2.0E0*I_ERI_D2z_Py_Fxyz_Px_d;
  abcd[1164] = 2.0E0*I_ERI_D2x_Pz_Fxyz_Px_d;
  abcd[1165] = 2.0E0*I_ERI_Dxy_Pz_Fxyz_Px_d;
  abcd[1166] = 2.0E0*I_ERI_Dxz_Pz_Fxyz_Px_d;
  abcd[1167] = 2.0E0*I_ERI_D2y_Pz_Fxyz_Px_d;
  abcd[1168] = 2.0E0*I_ERI_Dyz_Pz_Fxyz_Px_d;
  abcd[1169] = 2.0E0*I_ERI_D2z_Pz_Fxyz_Px_d;
  abcd[1170] = 2.0E0*I_ERI_D2x_Px_Fx2z_Px_d;
  abcd[1171] = 2.0E0*I_ERI_Dxy_Px_Fx2z_Px_d;
  abcd[1172] = 2.0E0*I_ERI_Dxz_Px_Fx2z_Px_d;
  abcd[1173] = 2.0E0*I_ERI_D2y_Px_Fx2z_Px_d;
  abcd[1174] = 2.0E0*I_ERI_Dyz_Px_Fx2z_Px_d;
  abcd[1175] = 2.0E0*I_ERI_D2z_Px_Fx2z_Px_d;
  abcd[1176] = 2.0E0*I_ERI_D2x_Py_Fx2z_Px_d;
  abcd[1177] = 2.0E0*I_ERI_Dxy_Py_Fx2z_Px_d;
  abcd[1178] = 2.0E0*I_ERI_Dxz_Py_Fx2z_Px_d;
  abcd[1179] = 2.0E0*I_ERI_D2y_Py_Fx2z_Px_d;
  abcd[1180] = 2.0E0*I_ERI_Dyz_Py_Fx2z_Px_d;
  abcd[1181] = 2.0E0*I_ERI_D2z_Py_Fx2z_Px_d;
  abcd[1182] = 2.0E0*I_ERI_D2x_Pz_Fx2z_Px_d;
  abcd[1183] = 2.0E0*I_ERI_Dxy_Pz_Fx2z_Px_d;
  abcd[1184] = 2.0E0*I_ERI_Dxz_Pz_Fx2z_Px_d;
  abcd[1185] = 2.0E0*I_ERI_D2y_Pz_Fx2z_Px_d;
  abcd[1186] = 2.0E0*I_ERI_Dyz_Pz_Fx2z_Px_d;
  abcd[1187] = 2.0E0*I_ERI_D2z_Pz_Fx2z_Px_d;
  abcd[1188] = 2.0E0*I_ERI_D2x_Px_F3y_Px_d;
  abcd[1189] = 2.0E0*I_ERI_Dxy_Px_F3y_Px_d;
  abcd[1190] = 2.0E0*I_ERI_Dxz_Px_F3y_Px_d;
  abcd[1191] = 2.0E0*I_ERI_D2y_Px_F3y_Px_d;
  abcd[1192] = 2.0E0*I_ERI_Dyz_Px_F3y_Px_d;
  abcd[1193] = 2.0E0*I_ERI_D2z_Px_F3y_Px_d;
  abcd[1194] = 2.0E0*I_ERI_D2x_Py_F3y_Px_d;
  abcd[1195] = 2.0E0*I_ERI_Dxy_Py_F3y_Px_d;
  abcd[1196] = 2.0E0*I_ERI_Dxz_Py_F3y_Px_d;
  abcd[1197] = 2.0E0*I_ERI_D2y_Py_F3y_Px_d;
  abcd[1198] = 2.0E0*I_ERI_Dyz_Py_F3y_Px_d;
  abcd[1199] = 2.0E0*I_ERI_D2z_Py_F3y_Px_d;
  abcd[1200] = 2.0E0*I_ERI_D2x_Pz_F3y_Px_d;
  abcd[1201] = 2.0E0*I_ERI_Dxy_Pz_F3y_Px_d;
  abcd[1202] = 2.0E0*I_ERI_Dxz_Pz_F3y_Px_d;
  abcd[1203] = 2.0E0*I_ERI_D2y_Pz_F3y_Px_d;
  abcd[1204] = 2.0E0*I_ERI_Dyz_Pz_F3y_Px_d;
  abcd[1205] = 2.0E0*I_ERI_D2z_Pz_F3y_Px_d;
  abcd[1206] = 2.0E0*I_ERI_D2x_Px_F2yz_Px_d;
  abcd[1207] = 2.0E0*I_ERI_Dxy_Px_F2yz_Px_d;
  abcd[1208] = 2.0E0*I_ERI_Dxz_Px_F2yz_Px_d;
  abcd[1209] = 2.0E0*I_ERI_D2y_Px_F2yz_Px_d;
  abcd[1210] = 2.0E0*I_ERI_Dyz_Px_F2yz_Px_d;
  abcd[1211] = 2.0E0*I_ERI_D2z_Px_F2yz_Px_d;
  abcd[1212] = 2.0E0*I_ERI_D2x_Py_F2yz_Px_d;
  abcd[1213] = 2.0E0*I_ERI_Dxy_Py_F2yz_Px_d;
  abcd[1214] = 2.0E0*I_ERI_Dxz_Py_F2yz_Px_d;
  abcd[1215] = 2.0E0*I_ERI_D2y_Py_F2yz_Px_d;
  abcd[1216] = 2.0E0*I_ERI_Dyz_Py_F2yz_Px_d;
  abcd[1217] = 2.0E0*I_ERI_D2z_Py_F2yz_Px_d;
  abcd[1218] = 2.0E0*I_ERI_D2x_Pz_F2yz_Px_d;
  abcd[1219] = 2.0E0*I_ERI_Dxy_Pz_F2yz_Px_d;
  abcd[1220] = 2.0E0*I_ERI_Dxz_Pz_F2yz_Px_d;
  abcd[1221] = 2.0E0*I_ERI_D2y_Pz_F2yz_Px_d;
  abcd[1222] = 2.0E0*I_ERI_Dyz_Pz_F2yz_Px_d;
  abcd[1223] = 2.0E0*I_ERI_D2z_Pz_F2yz_Px_d;
  abcd[1224] = 2.0E0*I_ERI_D2x_Px_Fy2z_Px_d;
  abcd[1225] = 2.0E0*I_ERI_Dxy_Px_Fy2z_Px_d;
  abcd[1226] = 2.0E0*I_ERI_Dxz_Px_Fy2z_Px_d;
  abcd[1227] = 2.0E0*I_ERI_D2y_Px_Fy2z_Px_d;
  abcd[1228] = 2.0E0*I_ERI_Dyz_Px_Fy2z_Px_d;
  abcd[1229] = 2.0E0*I_ERI_D2z_Px_Fy2z_Px_d;
  abcd[1230] = 2.0E0*I_ERI_D2x_Py_Fy2z_Px_d;
  abcd[1231] = 2.0E0*I_ERI_Dxy_Py_Fy2z_Px_d;
  abcd[1232] = 2.0E0*I_ERI_Dxz_Py_Fy2z_Px_d;
  abcd[1233] = 2.0E0*I_ERI_D2y_Py_Fy2z_Px_d;
  abcd[1234] = 2.0E0*I_ERI_Dyz_Py_Fy2z_Px_d;
  abcd[1235] = 2.0E0*I_ERI_D2z_Py_Fy2z_Px_d;
  abcd[1236] = 2.0E0*I_ERI_D2x_Pz_Fy2z_Px_d;
  abcd[1237] = 2.0E0*I_ERI_Dxy_Pz_Fy2z_Px_d;
  abcd[1238] = 2.0E0*I_ERI_Dxz_Pz_Fy2z_Px_d;
  abcd[1239] = 2.0E0*I_ERI_D2y_Pz_Fy2z_Px_d;
  abcd[1240] = 2.0E0*I_ERI_Dyz_Pz_Fy2z_Px_d;
  abcd[1241] = 2.0E0*I_ERI_D2z_Pz_Fy2z_Px_d;
  abcd[1242] = 2.0E0*I_ERI_D2x_Px_F3z_Px_d;
  abcd[1243] = 2.0E0*I_ERI_Dxy_Px_F3z_Px_d;
  abcd[1244] = 2.0E0*I_ERI_Dxz_Px_F3z_Px_d;
  abcd[1245] = 2.0E0*I_ERI_D2y_Px_F3z_Px_d;
  abcd[1246] = 2.0E0*I_ERI_Dyz_Px_F3z_Px_d;
  abcd[1247] = 2.0E0*I_ERI_D2z_Px_F3z_Px_d;
  abcd[1248] = 2.0E0*I_ERI_D2x_Py_F3z_Px_d;
  abcd[1249] = 2.0E0*I_ERI_Dxy_Py_F3z_Px_d;
  abcd[1250] = 2.0E0*I_ERI_Dxz_Py_F3z_Px_d;
  abcd[1251] = 2.0E0*I_ERI_D2y_Py_F3z_Px_d;
  abcd[1252] = 2.0E0*I_ERI_Dyz_Py_F3z_Px_d;
  abcd[1253] = 2.0E0*I_ERI_D2z_Py_F3z_Px_d;
  abcd[1254] = 2.0E0*I_ERI_D2x_Pz_F3z_Px_d;
  abcd[1255] = 2.0E0*I_ERI_Dxy_Pz_F3z_Px_d;
  abcd[1256] = 2.0E0*I_ERI_Dxz_Pz_F3z_Px_d;
  abcd[1257] = 2.0E0*I_ERI_D2y_Pz_F3z_Px_d;
  abcd[1258] = 2.0E0*I_ERI_Dyz_Pz_F3z_Px_d;
  abcd[1259] = 2.0E0*I_ERI_D2z_Pz_F3z_Px_d;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_ddy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_P_d
   ************************************************************/
  abcd[1260] = 2.0E0*I_ERI_D2x_Px_F3x_Py_d;
  abcd[1261] = 2.0E0*I_ERI_Dxy_Px_F3x_Py_d;
  abcd[1262] = 2.0E0*I_ERI_Dxz_Px_F3x_Py_d;
  abcd[1263] = 2.0E0*I_ERI_D2y_Px_F3x_Py_d;
  abcd[1264] = 2.0E0*I_ERI_Dyz_Px_F3x_Py_d;
  abcd[1265] = 2.0E0*I_ERI_D2z_Px_F3x_Py_d;
  abcd[1266] = 2.0E0*I_ERI_D2x_Py_F3x_Py_d;
  abcd[1267] = 2.0E0*I_ERI_Dxy_Py_F3x_Py_d;
  abcd[1268] = 2.0E0*I_ERI_Dxz_Py_F3x_Py_d;
  abcd[1269] = 2.0E0*I_ERI_D2y_Py_F3x_Py_d;
  abcd[1270] = 2.0E0*I_ERI_Dyz_Py_F3x_Py_d;
  abcd[1271] = 2.0E0*I_ERI_D2z_Py_F3x_Py_d;
  abcd[1272] = 2.0E0*I_ERI_D2x_Pz_F3x_Py_d;
  abcd[1273] = 2.0E0*I_ERI_Dxy_Pz_F3x_Py_d;
  abcd[1274] = 2.0E0*I_ERI_Dxz_Pz_F3x_Py_d;
  abcd[1275] = 2.0E0*I_ERI_D2y_Pz_F3x_Py_d;
  abcd[1276] = 2.0E0*I_ERI_Dyz_Pz_F3x_Py_d;
  abcd[1277] = 2.0E0*I_ERI_D2z_Pz_F3x_Py_d;
  abcd[1278] = 2.0E0*I_ERI_D2x_Px_F2xy_Py_d;
  abcd[1279] = 2.0E0*I_ERI_Dxy_Px_F2xy_Py_d;
  abcd[1280] = 2.0E0*I_ERI_Dxz_Px_F2xy_Py_d;
  abcd[1281] = 2.0E0*I_ERI_D2y_Px_F2xy_Py_d;
  abcd[1282] = 2.0E0*I_ERI_Dyz_Px_F2xy_Py_d;
  abcd[1283] = 2.0E0*I_ERI_D2z_Px_F2xy_Py_d;
  abcd[1284] = 2.0E0*I_ERI_D2x_Py_F2xy_Py_d;
  abcd[1285] = 2.0E0*I_ERI_Dxy_Py_F2xy_Py_d;
  abcd[1286] = 2.0E0*I_ERI_Dxz_Py_F2xy_Py_d;
  abcd[1287] = 2.0E0*I_ERI_D2y_Py_F2xy_Py_d;
  abcd[1288] = 2.0E0*I_ERI_Dyz_Py_F2xy_Py_d;
  abcd[1289] = 2.0E0*I_ERI_D2z_Py_F2xy_Py_d;
  abcd[1290] = 2.0E0*I_ERI_D2x_Pz_F2xy_Py_d;
  abcd[1291] = 2.0E0*I_ERI_Dxy_Pz_F2xy_Py_d;
  abcd[1292] = 2.0E0*I_ERI_Dxz_Pz_F2xy_Py_d;
  abcd[1293] = 2.0E0*I_ERI_D2y_Pz_F2xy_Py_d;
  abcd[1294] = 2.0E0*I_ERI_Dyz_Pz_F2xy_Py_d;
  abcd[1295] = 2.0E0*I_ERI_D2z_Pz_F2xy_Py_d;
  abcd[1296] = 2.0E0*I_ERI_D2x_Px_F2xz_Py_d;
  abcd[1297] = 2.0E0*I_ERI_Dxy_Px_F2xz_Py_d;
  abcd[1298] = 2.0E0*I_ERI_Dxz_Px_F2xz_Py_d;
  abcd[1299] = 2.0E0*I_ERI_D2y_Px_F2xz_Py_d;
  abcd[1300] = 2.0E0*I_ERI_Dyz_Px_F2xz_Py_d;
  abcd[1301] = 2.0E0*I_ERI_D2z_Px_F2xz_Py_d;
  abcd[1302] = 2.0E0*I_ERI_D2x_Py_F2xz_Py_d;
  abcd[1303] = 2.0E0*I_ERI_Dxy_Py_F2xz_Py_d;
  abcd[1304] = 2.0E0*I_ERI_Dxz_Py_F2xz_Py_d;
  abcd[1305] = 2.0E0*I_ERI_D2y_Py_F2xz_Py_d;
  abcd[1306] = 2.0E0*I_ERI_Dyz_Py_F2xz_Py_d;
  abcd[1307] = 2.0E0*I_ERI_D2z_Py_F2xz_Py_d;
  abcd[1308] = 2.0E0*I_ERI_D2x_Pz_F2xz_Py_d;
  abcd[1309] = 2.0E0*I_ERI_Dxy_Pz_F2xz_Py_d;
  abcd[1310] = 2.0E0*I_ERI_Dxz_Pz_F2xz_Py_d;
  abcd[1311] = 2.0E0*I_ERI_D2y_Pz_F2xz_Py_d;
  abcd[1312] = 2.0E0*I_ERI_Dyz_Pz_F2xz_Py_d;
  abcd[1313] = 2.0E0*I_ERI_D2z_Pz_F2xz_Py_d;
  abcd[1314] = 2.0E0*I_ERI_D2x_Px_Fx2y_Py_d;
  abcd[1315] = 2.0E0*I_ERI_Dxy_Px_Fx2y_Py_d;
  abcd[1316] = 2.0E0*I_ERI_Dxz_Px_Fx2y_Py_d;
  abcd[1317] = 2.0E0*I_ERI_D2y_Px_Fx2y_Py_d;
  abcd[1318] = 2.0E0*I_ERI_Dyz_Px_Fx2y_Py_d;
  abcd[1319] = 2.0E0*I_ERI_D2z_Px_Fx2y_Py_d;
  abcd[1320] = 2.0E0*I_ERI_D2x_Py_Fx2y_Py_d;
  abcd[1321] = 2.0E0*I_ERI_Dxy_Py_Fx2y_Py_d;
  abcd[1322] = 2.0E0*I_ERI_Dxz_Py_Fx2y_Py_d;
  abcd[1323] = 2.0E0*I_ERI_D2y_Py_Fx2y_Py_d;
  abcd[1324] = 2.0E0*I_ERI_Dyz_Py_Fx2y_Py_d;
  abcd[1325] = 2.0E0*I_ERI_D2z_Py_Fx2y_Py_d;
  abcd[1326] = 2.0E0*I_ERI_D2x_Pz_Fx2y_Py_d;
  abcd[1327] = 2.0E0*I_ERI_Dxy_Pz_Fx2y_Py_d;
  abcd[1328] = 2.0E0*I_ERI_Dxz_Pz_Fx2y_Py_d;
  abcd[1329] = 2.0E0*I_ERI_D2y_Pz_Fx2y_Py_d;
  abcd[1330] = 2.0E0*I_ERI_Dyz_Pz_Fx2y_Py_d;
  abcd[1331] = 2.0E0*I_ERI_D2z_Pz_Fx2y_Py_d;
  abcd[1332] = 2.0E0*I_ERI_D2x_Px_Fxyz_Py_d;
  abcd[1333] = 2.0E0*I_ERI_Dxy_Px_Fxyz_Py_d;
  abcd[1334] = 2.0E0*I_ERI_Dxz_Px_Fxyz_Py_d;
  abcd[1335] = 2.0E0*I_ERI_D2y_Px_Fxyz_Py_d;
  abcd[1336] = 2.0E0*I_ERI_Dyz_Px_Fxyz_Py_d;
  abcd[1337] = 2.0E0*I_ERI_D2z_Px_Fxyz_Py_d;
  abcd[1338] = 2.0E0*I_ERI_D2x_Py_Fxyz_Py_d;
  abcd[1339] = 2.0E0*I_ERI_Dxy_Py_Fxyz_Py_d;
  abcd[1340] = 2.0E0*I_ERI_Dxz_Py_Fxyz_Py_d;
  abcd[1341] = 2.0E0*I_ERI_D2y_Py_Fxyz_Py_d;
  abcd[1342] = 2.0E0*I_ERI_Dyz_Py_Fxyz_Py_d;
  abcd[1343] = 2.0E0*I_ERI_D2z_Py_Fxyz_Py_d;
  abcd[1344] = 2.0E0*I_ERI_D2x_Pz_Fxyz_Py_d;
  abcd[1345] = 2.0E0*I_ERI_Dxy_Pz_Fxyz_Py_d;
  abcd[1346] = 2.0E0*I_ERI_Dxz_Pz_Fxyz_Py_d;
  abcd[1347] = 2.0E0*I_ERI_D2y_Pz_Fxyz_Py_d;
  abcd[1348] = 2.0E0*I_ERI_Dyz_Pz_Fxyz_Py_d;
  abcd[1349] = 2.0E0*I_ERI_D2z_Pz_Fxyz_Py_d;
  abcd[1350] = 2.0E0*I_ERI_D2x_Px_Fx2z_Py_d;
  abcd[1351] = 2.0E0*I_ERI_Dxy_Px_Fx2z_Py_d;
  abcd[1352] = 2.0E0*I_ERI_Dxz_Px_Fx2z_Py_d;
  abcd[1353] = 2.0E0*I_ERI_D2y_Px_Fx2z_Py_d;
  abcd[1354] = 2.0E0*I_ERI_Dyz_Px_Fx2z_Py_d;
  abcd[1355] = 2.0E0*I_ERI_D2z_Px_Fx2z_Py_d;
  abcd[1356] = 2.0E0*I_ERI_D2x_Py_Fx2z_Py_d;
  abcd[1357] = 2.0E0*I_ERI_Dxy_Py_Fx2z_Py_d;
  abcd[1358] = 2.0E0*I_ERI_Dxz_Py_Fx2z_Py_d;
  abcd[1359] = 2.0E0*I_ERI_D2y_Py_Fx2z_Py_d;
  abcd[1360] = 2.0E0*I_ERI_Dyz_Py_Fx2z_Py_d;
  abcd[1361] = 2.0E0*I_ERI_D2z_Py_Fx2z_Py_d;
  abcd[1362] = 2.0E0*I_ERI_D2x_Pz_Fx2z_Py_d;
  abcd[1363] = 2.0E0*I_ERI_Dxy_Pz_Fx2z_Py_d;
  abcd[1364] = 2.0E0*I_ERI_Dxz_Pz_Fx2z_Py_d;
  abcd[1365] = 2.0E0*I_ERI_D2y_Pz_Fx2z_Py_d;
  abcd[1366] = 2.0E0*I_ERI_Dyz_Pz_Fx2z_Py_d;
  abcd[1367] = 2.0E0*I_ERI_D2z_Pz_Fx2z_Py_d;
  abcd[1368] = 2.0E0*I_ERI_D2x_Px_F3y_Py_d;
  abcd[1369] = 2.0E0*I_ERI_Dxy_Px_F3y_Py_d;
  abcd[1370] = 2.0E0*I_ERI_Dxz_Px_F3y_Py_d;
  abcd[1371] = 2.0E0*I_ERI_D2y_Px_F3y_Py_d;
  abcd[1372] = 2.0E0*I_ERI_Dyz_Px_F3y_Py_d;
  abcd[1373] = 2.0E0*I_ERI_D2z_Px_F3y_Py_d;
  abcd[1374] = 2.0E0*I_ERI_D2x_Py_F3y_Py_d;
  abcd[1375] = 2.0E0*I_ERI_Dxy_Py_F3y_Py_d;
  abcd[1376] = 2.0E0*I_ERI_Dxz_Py_F3y_Py_d;
  abcd[1377] = 2.0E0*I_ERI_D2y_Py_F3y_Py_d;
  abcd[1378] = 2.0E0*I_ERI_Dyz_Py_F3y_Py_d;
  abcd[1379] = 2.0E0*I_ERI_D2z_Py_F3y_Py_d;
  abcd[1380] = 2.0E0*I_ERI_D2x_Pz_F3y_Py_d;
  abcd[1381] = 2.0E0*I_ERI_Dxy_Pz_F3y_Py_d;
  abcd[1382] = 2.0E0*I_ERI_Dxz_Pz_F3y_Py_d;
  abcd[1383] = 2.0E0*I_ERI_D2y_Pz_F3y_Py_d;
  abcd[1384] = 2.0E0*I_ERI_Dyz_Pz_F3y_Py_d;
  abcd[1385] = 2.0E0*I_ERI_D2z_Pz_F3y_Py_d;
  abcd[1386] = 2.0E0*I_ERI_D2x_Px_F2yz_Py_d;
  abcd[1387] = 2.0E0*I_ERI_Dxy_Px_F2yz_Py_d;
  abcd[1388] = 2.0E0*I_ERI_Dxz_Px_F2yz_Py_d;
  abcd[1389] = 2.0E0*I_ERI_D2y_Px_F2yz_Py_d;
  abcd[1390] = 2.0E0*I_ERI_Dyz_Px_F2yz_Py_d;
  abcd[1391] = 2.0E0*I_ERI_D2z_Px_F2yz_Py_d;
  abcd[1392] = 2.0E0*I_ERI_D2x_Py_F2yz_Py_d;
  abcd[1393] = 2.0E0*I_ERI_Dxy_Py_F2yz_Py_d;
  abcd[1394] = 2.0E0*I_ERI_Dxz_Py_F2yz_Py_d;
  abcd[1395] = 2.0E0*I_ERI_D2y_Py_F2yz_Py_d;
  abcd[1396] = 2.0E0*I_ERI_Dyz_Py_F2yz_Py_d;
  abcd[1397] = 2.0E0*I_ERI_D2z_Py_F2yz_Py_d;
  abcd[1398] = 2.0E0*I_ERI_D2x_Pz_F2yz_Py_d;
  abcd[1399] = 2.0E0*I_ERI_Dxy_Pz_F2yz_Py_d;
  abcd[1400] = 2.0E0*I_ERI_Dxz_Pz_F2yz_Py_d;
  abcd[1401] = 2.0E0*I_ERI_D2y_Pz_F2yz_Py_d;
  abcd[1402] = 2.0E0*I_ERI_Dyz_Pz_F2yz_Py_d;
  abcd[1403] = 2.0E0*I_ERI_D2z_Pz_F2yz_Py_d;
  abcd[1404] = 2.0E0*I_ERI_D2x_Px_Fy2z_Py_d;
  abcd[1405] = 2.0E0*I_ERI_Dxy_Px_Fy2z_Py_d;
  abcd[1406] = 2.0E0*I_ERI_Dxz_Px_Fy2z_Py_d;
  abcd[1407] = 2.0E0*I_ERI_D2y_Px_Fy2z_Py_d;
  abcd[1408] = 2.0E0*I_ERI_Dyz_Px_Fy2z_Py_d;
  abcd[1409] = 2.0E0*I_ERI_D2z_Px_Fy2z_Py_d;
  abcd[1410] = 2.0E0*I_ERI_D2x_Py_Fy2z_Py_d;
  abcd[1411] = 2.0E0*I_ERI_Dxy_Py_Fy2z_Py_d;
  abcd[1412] = 2.0E0*I_ERI_Dxz_Py_Fy2z_Py_d;
  abcd[1413] = 2.0E0*I_ERI_D2y_Py_Fy2z_Py_d;
  abcd[1414] = 2.0E0*I_ERI_Dyz_Py_Fy2z_Py_d;
  abcd[1415] = 2.0E0*I_ERI_D2z_Py_Fy2z_Py_d;
  abcd[1416] = 2.0E0*I_ERI_D2x_Pz_Fy2z_Py_d;
  abcd[1417] = 2.0E0*I_ERI_Dxy_Pz_Fy2z_Py_d;
  abcd[1418] = 2.0E0*I_ERI_Dxz_Pz_Fy2z_Py_d;
  abcd[1419] = 2.0E0*I_ERI_D2y_Pz_Fy2z_Py_d;
  abcd[1420] = 2.0E0*I_ERI_Dyz_Pz_Fy2z_Py_d;
  abcd[1421] = 2.0E0*I_ERI_D2z_Pz_Fy2z_Py_d;
  abcd[1422] = 2.0E0*I_ERI_D2x_Px_F3z_Py_d;
  abcd[1423] = 2.0E0*I_ERI_Dxy_Px_F3z_Py_d;
  abcd[1424] = 2.0E0*I_ERI_Dxz_Px_F3z_Py_d;
  abcd[1425] = 2.0E0*I_ERI_D2y_Px_F3z_Py_d;
  abcd[1426] = 2.0E0*I_ERI_Dyz_Px_F3z_Py_d;
  abcd[1427] = 2.0E0*I_ERI_D2z_Px_F3z_Py_d;
  abcd[1428] = 2.0E0*I_ERI_D2x_Py_F3z_Py_d;
  abcd[1429] = 2.0E0*I_ERI_Dxy_Py_F3z_Py_d;
  abcd[1430] = 2.0E0*I_ERI_Dxz_Py_F3z_Py_d;
  abcd[1431] = 2.0E0*I_ERI_D2y_Py_F3z_Py_d;
  abcd[1432] = 2.0E0*I_ERI_Dyz_Py_F3z_Py_d;
  abcd[1433] = 2.0E0*I_ERI_D2z_Py_F3z_Py_d;
  abcd[1434] = 2.0E0*I_ERI_D2x_Pz_F3z_Py_d;
  abcd[1435] = 2.0E0*I_ERI_Dxy_Pz_F3z_Py_d;
  abcd[1436] = 2.0E0*I_ERI_Dxz_Pz_F3z_Py_d;
  abcd[1437] = 2.0E0*I_ERI_D2y_Pz_F3z_Py_d;
  abcd[1438] = 2.0E0*I_ERI_Dyz_Pz_F3z_Py_d;
  abcd[1439] = 2.0E0*I_ERI_D2z_Pz_F3z_Py_d;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_ddz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_P_d
   ************************************************************/
  abcd[1440] = 2.0E0*I_ERI_D2x_Px_F3x_Pz_d;
  abcd[1441] = 2.0E0*I_ERI_Dxy_Px_F3x_Pz_d;
  abcd[1442] = 2.0E0*I_ERI_Dxz_Px_F3x_Pz_d;
  abcd[1443] = 2.0E0*I_ERI_D2y_Px_F3x_Pz_d;
  abcd[1444] = 2.0E0*I_ERI_Dyz_Px_F3x_Pz_d;
  abcd[1445] = 2.0E0*I_ERI_D2z_Px_F3x_Pz_d;
  abcd[1446] = 2.0E0*I_ERI_D2x_Py_F3x_Pz_d;
  abcd[1447] = 2.0E0*I_ERI_Dxy_Py_F3x_Pz_d;
  abcd[1448] = 2.0E0*I_ERI_Dxz_Py_F3x_Pz_d;
  abcd[1449] = 2.0E0*I_ERI_D2y_Py_F3x_Pz_d;
  abcd[1450] = 2.0E0*I_ERI_Dyz_Py_F3x_Pz_d;
  abcd[1451] = 2.0E0*I_ERI_D2z_Py_F3x_Pz_d;
  abcd[1452] = 2.0E0*I_ERI_D2x_Pz_F3x_Pz_d;
  abcd[1453] = 2.0E0*I_ERI_Dxy_Pz_F3x_Pz_d;
  abcd[1454] = 2.0E0*I_ERI_Dxz_Pz_F3x_Pz_d;
  abcd[1455] = 2.0E0*I_ERI_D2y_Pz_F3x_Pz_d;
  abcd[1456] = 2.0E0*I_ERI_Dyz_Pz_F3x_Pz_d;
  abcd[1457] = 2.0E0*I_ERI_D2z_Pz_F3x_Pz_d;
  abcd[1458] = 2.0E0*I_ERI_D2x_Px_F2xy_Pz_d;
  abcd[1459] = 2.0E0*I_ERI_Dxy_Px_F2xy_Pz_d;
  abcd[1460] = 2.0E0*I_ERI_Dxz_Px_F2xy_Pz_d;
  abcd[1461] = 2.0E0*I_ERI_D2y_Px_F2xy_Pz_d;
  abcd[1462] = 2.0E0*I_ERI_Dyz_Px_F2xy_Pz_d;
  abcd[1463] = 2.0E0*I_ERI_D2z_Px_F2xy_Pz_d;
  abcd[1464] = 2.0E0*I_ERI_D2x_Py_F2xy_Pz_d;
  abcd[1465] = 2.0E0*I_ERI_Dxy_Py_F2xy_Pz_d;
  abcd[1466] = 2.0E0*I_ERI_Dxz_Py_F2xy_Pz_d;
  abcd[1467] = 2.0E0*I_ERI_D2y_Py_F2xy_Pz_d;
  abcd[1468] = 2.0E0*I_ERI_Dyz_Py_F2xy_Pz_d;
  abcd[1469] = 2.0E0*I_ERI_D2z_Py_F2xy_Pz_d;
  abcd[1470] = 2.0E0*I_ERI_D2x_Pz_F2xy_Pz_d;
  abcd[1471] = 2.0E0*I_ERI_Dxy_Pz_F2xy_Pz_d;
  abcd[1472] = 2.0E0*I_ERI_Dxz_Pz_F2xy_Pz_d;
  abcd[1473] = 2.0E0*I_ERI_D2y_Pz_F2xy_Pz_d;
  abcd[1474] = 2.0E0*I_ERI_Dyz_Pz_F2xy_Pz_d;
  abcd[1475] = 2.0E0*I_ERI_D2z_Pz_F2xy_Pz_d;
  abcd[1476] = 2.0E0*I_ERI_D2x_Px_F2xz_Pz_d;
  abcd[1477] = 2.0E0*I_ERI_Dxy_Px_F2xz_Pz_d;
  abcd[1478] = 2.0E0*I_ERI_Dxz_Px_F2xz_Pz_d;
  abcd[1479] = 2.0E0*I_ERI_D2y_Px_F2xz_Pz_d;
  abcd[1480] = 2.0E0*I_ERI_Dyz_Px_F2xz_Pz_d;
  abcd[1481] = 2.0E0*I_ERI_D2z_Px_F2xz_Pz_d;
  abcd[1482] = 2.0E0*I_ERI_D2x_Py_F2xz_Pz_d;
  abcd[1483] = 2.0E0*I_ERI_Dxy_Py_F2xz_Pz_d;
  abcd[1484] = 2.0E0*I_ERI_Dxz_Py_F2xz_Pz_d;
  abcd[1485] = 2.0E0*I_ERI_D2y_Py_F2xz_Pz_d;
  abcd[1486] = 2.0E0*I_ERI_Dyz_Py_F2xz_Pz_d;
  abcd[1487] = 2.0E0*I_ERI_D2z_Py_F2xz_Pz_d;
  abcd[1488] = 2.0E0*I_ERI_D2x_Pz_F2xz_Pz_d;
  abcd[1489] = 2.0E0*I_ERI_Dxy_Pz_F2xz_Pz_d;
  abcd[1490] = 2.0E0*I_ERI_Dxz_Pz_F2xz_Pz_d;
  abcd[1491] = 2.0E0*I_ERI_D2y_Pz_F2xz_Pz_d;
  abcd[1492] = 2.0E0*I_ERI_Dyz_Pz_F2xz_Pz_d;
  abcd[1493] = 2.0E0*I_ERI_D2z_Pz_F2xz_Pz_d;
  abcd[1494] = 2.0E0*I_ERI_D2x_Px_Fx2y_Pz_d;
  abcd[1495] = 2.0E0*I_ERI_Dxy_Px_Fx2y_Pz_d;
  abcd[1496] = 2.0E0*I_ERI_Dxz_Px_Fx2y_Pz_d;
  abcd[1497] = 2.0E0*I_ERI_D2y_Px_Fx2y_Pz_d;
  abcd[1498] = 2.0E0*I_ERI_Dyz_Px_Fx2y_Pz_d;
  abcd[1499] = 2.0E0*I_ERI_D2z_Px_Fx2y_Pz_d;
  abcd[1500] = 2.0E0*I_ERI_D2x_Py_Fx2y_Pz_d;
  abcd[1501] = 2.0E0*I_ERI_Dxy_Py_Fx2y_Pz_d;
  abcd[1502] = 2.0E0*I_ERI_Dxz_Py_Fx2y_Pz_d;
  abcd[1503] = 2.0E0*I_ERI_D2y_Py_Fx2y_Pz_d;
  abcd[1504] = 2.0E0*I_ERI_Dyz_Py_Fx2y_Pz_d;
  abcd[1505] = 2.0E0*I_ERI_D2z_Py_Fx2y_Pz_d;
  abcd[1506] = 2.0E0*I_ERI_D2x_Pz_Fx2y_Pz_d;
  abcd[1507] = 2.0E0*I_ERI_Dxy_Pz_Fx2y_Pz_d;
  abcd[1508] = 2.0E0*I_ERI_Dxz_Pz_Fx2y_Pz_d;
  abcd[1509] = 2.0E0*I_ERI_D2y_Pz_Fx2y_Pz_d;
  abcd[1510] = 2.0E0*I_ERI_Dyz_Pz_Fx2y_Pz_d;
  abcd[1511] = 2.0E0*I_ERI_D2z_Pz_Fx2y_Pz_d;
  abcd[1512] = 2.0E0*I_ERI_D2x_Px_Fxyz_Pz_d;
  abcd[1513] = 2.0E0*I_ERI_Dxy_Px_Fxyz_Pz_d;
  abcd[1514] = 2.0E0*I_ERI_Dxz_Px_Fxyz_Pz_d;
  abcd[1515] = 2.0E0*I_ERI_D2y_Px_Fxyz_Pz_d;
  abcd[1516] = 2.0E0*I_ERI_Dyz_Px_Fxyz_Pz_d;
  abcd[1517] = 2.0E0*I_ERI_D2z_Px_Fxyz_Pz_d;
  abcd[1518] = 2.0E0*I_ERI_D2x_Py_Fxyz_Pz_d;
  abcd[1519] = 2.0E0*I_ERI_Dxy_Py_Fxyz_Pz_d;
  abcd[1520] = 2.0E0*I_ERI_Dxz_Py_Fxyz_Pz_d;
  abcd[1521] = 2.0E0*I_ERI_D2y_Py_Fxyz_Pz_d;
  abcd[1522] = 2.0E0*I_ERI_Dyz_Py_Fxyz_Pz_d;
  abcd[1523] = 2.0E0*I_ERI_D2z_Py_Fxyz_Pz_d;
  abcd[1524] = 2.0E0*I_ERI_D2x_Pz_Fxyz_Pz_d;
  abcd[1525] = 2.0E0*I_ERI_Dxy_Pz_Fxyz_Pz_d;
  abcd[1526] = 2.0E0*I_ERI_Dxz_Pz_Fxyz_Pz_d;
  abcd[1527] = 2.0E0*I_ERI_D2y_Pz_Fxyz_Pz_d;
  abcd[1528] = 2.0E0*I_ERI_Dyz_Pz_Fxyz_Pz_d;
  abcd[1529] = 2.0E0*I_ERI_D2z_Pz_Fxyz_Pz_d;
  abcd[1530] = 2.0E0*I_ERI_D2x_Px_Fx2z_Pz_d;
  abcd[1531] = 2.0E0*I_ERI_Dxy_Px_Fx2z_Pz_d;
  abcd[1532] = 2.0E0*I_ERI_Dxz_Px_Fx2z_Pz_d;
  abcd[1533] = 2.0E0*I_ERI_D2y_Px_Fx2z_Pz_d;
  abcd[1534] = 2.0E0*I_ERI_Dyz_Px_Fx2z_Pz_d;
  abcd[1535] = 2.0E0*I_ERI_D2z_Px_Fx2z_Pz_d;
  abcd[1536] = 2.0E0*I_ERI_D2x_Py_Fx2z_Pz_d;
  abcd[1537] = 2.0E0*I_ERI_Dxy_Py_Fx2z_Pz_d;
  abcd[1538] = 2.0E0*I_ERI_Dxz_Py_Fx2z_Pz_d;
  abcd[1539] = 2.0E0*I_ERI_D2y_Py_Fx2z_Pz_d;
  abcd[1540] = 2.0E0*I_ERI_Dyz_Py_Fx2z_Pz_d;
  abcd[1541] = 2.0E0*I_ERI_D2z_Py_Fx2z_Pz_d;
  abcd[1542] = 2.0E0*I_ERI_D2x_Pz_Fx2z_Pz_d;
  abcd[1543] = 2.0E0*I_ERI_Dxy_Pz_Fx2z_Pz_d;
  abcd[1544] = 2.0E0*I_ERI_Dxz_Pz_Fx2z_Pz_d;
  abcd[1545] = 2.0E0*I_ERI_D2y_Pz_Fx2z_Pz_d;
  abcd[1546] = 2.0E0*I_ERI_Dyz_Pz_Fx2z_Pz_d;
  abcd[1547] = 2.0E0*I_ERI_D2z_Pz_Fx2z_Pz_d;
  abcd[1548] = 2.0E0*I_ERI_D2x_Px_F3y_Pz_d;
  abcd[1549] = 2.0E0*I_ERI_Dxy_Px_F3y_Pz_d;
  abcd[1550] = 2.0E0*I_ERI_Dxz_Px_F3y_Pz_d;
  abcd[1551] = 2.0E0*I_ERI_D2y_Px_F3y_Pz_d;
  abcd[1552] = 2.0E0*I_ERI_Dyz_Px_F3y_Pz_d;
  abcd[1553] = 2.0E0*I_ERI_D2z_Px_F3y_Pz_d;
  abcd[1554] = 2.0E0*I_ERI_D2x_Py_F3y_Pz_d;
  abcd[1555] = 2.0E0*I_ERI_Dxy_Py_F3y_Pz_d;
  abcd[1556] = 2.0E0*I_ERI_Dxz_Py_F3y_Pz_d;
  abcd[1557] = 2.0E0*I_ERI_D2y_Py_F3y_Pz_d;
  abcd[1558] = 2.0E0*I_ERI_Dyz_Py_F3y_Pz_d;
  abcd[1559] = 2.0E0*I_ERI_D2z_Py_F3y_Pz_d;
  abcd[1560] = 2.0E0*I_ERI_D2x_Pz_F3y_Pz_d;
  abcd[1561] = 2.0E0*I_ERI_Dxy_Pz_F3y_Pz_d;
  abcd[1562] = 2.0E0*I_ERI_Dxz_Pz_F3y_Pz_d;
  abcd[1563] = 2.0E0*I_ERI_D2y_Pz_F3y_Pz_d;
  abcd[1564] = 2.0E0*I_ERI_Dyz_Pz_F3y_Pz_d;
  abcd[1565] = 2.0E0*I_ERI_D2z_Pz_F3y_Pz_d;
  abcd[1566] = 2.0E0*I_ERI_D2x_Px_F2yz_Pz_d;
  abcd[1567] = 2.0E0*I_ERI_Dxy_Px_F2yz_Pz_d;
  abcd[1568] = 2.0E0*I_ERI_Dxz_Px_F2yz_Pz_d;
  abcd[1569] = 2.0E0*I_ERI_D2y_Px_F2yz_Pz_d;
  abcd[1570] = 2.0E0*I_ERI_Dyz_Px_F2yz_Pz_d;
  abcd[1571] = 2.0E0*I_ERI_D2z_Px_F2yz_Pz_d;
  abcd[1572] = 2.0E0*I_ERI_D2x_Py_F2yz_Pz_d;
  abcd[1573] = 2.0E0*I_ERI_Dxy_Py_F2yz_Pz_d;
  abcd[1574] = 2.0E0*I_ERI_Dxz_Py_F2yz_Pz_d;
  abcd[1575] = 2.0E0*I_ERI_D2y_Py_F2yz_Pz_d;
  abcd[1576] = 2.0E0*I_ERI_Dyz_Py_F2yz_Pz_d;
  abcd[1577] = 2.0E0*I_ERI_D2z_Py_F2yz_Pz_d;
  abcd[1578] = 2.0E0*I_ERI_D2x_Pz_F2yz_Pz_d;
  abcd[1579] = 2.0E0*I_ERI_Dxy_Pz_F2yz_Pz_d;
  abcd[1580] = 2.0E0*I_ERI_Dxz_Pz_F2yz_Pz_d;
  abcd[1581] = 2.0E0*I_ERI_D2y_Pz_F2yz_Pz_d;
  abcd[1582] = 2.0E0*I_ERI_Dyz_Pz_F2yz_Pz_d;
  abcd[1583] = 2.0E0*I_ERI_D2z_Pz_F2yz_Pz_d;
  abcd[1584] = 2.0E0*I_ERI_D2x_Px_Fy2z_Pz_d;
  abcd[1585] = 2.0E0*I_ERI_Dxy_Px_Fy2z_Pz_d;
  abcd[1586] = 2.0E0*I_ERI_Dxz_Px_Fy2z_Pz_d;
  abcd[1587] = 2.0E0*I_ERI_D2y_Px_Fy2z_Pz_d;
  abcd[1588] = 2.0E0*I_ERI_Dyz_Px_Fy2z_Pz_d;
  abcd[1589] = 2.0E0*I_ERI_D2z_Px_Fy2z_Pz_d;
  abcd[1590] = 2.0E0*I_ERI_D2x_Py_Fy2z_Pz_d;
  abcd[1591] = 2.0E0*I_ERI_Dxy_Py_Fy2z_Pz_d;
  abcd[1592] = 2.0E0*I_ERI_Dxz_Py_Fy2z_Pz_d;
  abcd[1593] = 2.0E0*I_ERI_D2y_Py_Fy2z_Pz_d;
  abcd[1594] = 2.0E0*I_ERI_Dyz_Py_Fy2z_Pz_d;
  abcd[1595] = 2.0E0*I_ERI_D2z_Py_Fy2z_Pz_d;
  abcd[1596] = 2.0E0*I_ERI_D2x_Pz_Fy2z_Pz_d;
  abcd[1597] = 2.0E0*I_ERI_Dxy_Pz_Fy2z_Pz_d;
  abcd[1598] = 2.0E0*I_ERI_Dxz_Pz_Fy2z_Pz_d;
  abcd[1599] = 2.0E0*I_ERI_D2y_Pz_Fy2z_Pz_d;
  abcd[1600] = 2.0E0*I_ERI_Dyz_Pz_Fy2z_Pz_d;
  abcd[1601] = 2.0E0*I_ERI_D2z_Pz_Fy2z_Pz_d;
  abcd[1602] = 2.0E0*I_ERI_D2x_Px_F3z_Pz_d;
  abcd[1603] = 2.0E0*I_ERI_Dxy_Px_F3z_Pz_d;
  abcd[1604] = 2.0E0*I_ERI_Dxz_Px_F3z_Pz_d;
  abcd[1605] = 2.0E0*I_ERI_D2y_Px_F3z_Pz_d;
  abcd[1606] = 2.0E0*I_ERI_Dyz_Px_F3z_Pz_d;
  abcd[1607] = 2.0E0*I_ERI_D2z_Px_F3z_Pz_d;
  abcd[1608] = 2.0E0*I_ERI_D2x_Py_F3z_Pz_d;
  abcd[1609] = 2.0E0*I_ERI_Dxy_Py_F3z_Pz_d;
  abcd[1610] = 2.0E0*I_ERI_Dxz_Py_F3z_Pz_d;
  abcd[1611] = 2.0E0*I_ERI_D2y_Py_F3z_Pz_d;
  abcd[1612] = 2.0E0*I_ERI_Dyz_Py_F3z_Pz_d;
  abcd[1613] = 2.0E0*I_ERI_D2z_Py_F3z_Pz_d;
  abcd[1614] = 2.0E0*I_ERI_D2x_Pz_F3z_Pz_d;
  abcd[1615] = 2.0E0*I_ERI_Dxy_Pz_F3z_Pz_d;
  abcd[1616] = 2.0E0*I_ERI_Dxz_Pz_F3z_Pz_d;
  abcd[1617] = 2.0E0*I_ERI_D2y_Pz_F3z_Pz_d;
  abcd[1618] = 2.0E0*I_ERI_Dyz_Pz_F3z_Pz_d;
  abcd[1619] = 2.0E0*I_ERI_D2z_Pz_F3z_Pz_d;
}
