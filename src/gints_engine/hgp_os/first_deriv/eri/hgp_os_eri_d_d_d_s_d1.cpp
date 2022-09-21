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
// BRA1 as redundant position, total RHS integrals evaluated as: 22757
// BRA2 as redundant position, total RHS integrals evaluated as: 22757
// KET1 as redundant position, total RHS integrals evaluated as: 23163
// KET2 as redundant position, total RHS integrals evaluated as: 21581
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

void hgp_os_eri_d_d_d_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_D2x_S = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S = 0.0E0;
  Double I_ERI_F3y_S_D2x_S = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S = 0.0E0;
  Double I_ERI_F3z_S_D2x_S = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3x_S_D2y_S = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S = 0.0E0;
  Double I_ERI_F3y_S_D2y_S = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S = 0.0E0;
  Double I_ERI_F3z_S_D2y_S = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3x_S_D2z_S = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S = 0.0E0;
  Double I_ERI_F3y_S_D2z_S = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S = 0.0E0;
  Double I_ERI_F3z_S_D2z_S = 0.0E0;
  Double I_ERI_D2x_S_D2x_S = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S = 0.0E0;
  Double I_ERI_D2y_S_D2x_S = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S = 0.0E0;
  Double I_ERI_D2z_S_D2x_S = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S = 0.0E0;
  Double I_ERI_D2x_S_D2y_S = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S = 0.0E0;
  Double I_ERI_D2y_S_D2y_S = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S = 0.0E0;
  Double I_ERI_D2z_S_D2y_S = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S = 0.0E0;
  Double I_ERI_D2x_S_D2z_S = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S = 0.0E0;
  Double I_ERI_D2y_S_D2z_S = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S = 0.0E0;
  Double I_ERI_D2z_S_D2z_S = 0.0E0;
  Double I_ERI_H5x_S_D2x_S_a = 0.0E0;
  Double I_ERI_H4xy_S_D2x_S_a = 0.0E0;
  Double I_ERI_H4xz_S_D2x_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_D2x_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_D2x_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_D2x_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H5y_S_D2x_S_a = 0.0E0;
  Double I_ERI_H4yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H5z_S_D2x_S_a = 0.0E0;
  Double I_ERI_H5x_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H4xy_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H4xz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H5y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H4yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H5z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_H5x_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H4xy_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H4xz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H5y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H4yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H5z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_H5x_S_D2y_S_a = 0.0E0;
  Double I_ERI_H4xy_S_D2y_S_a = 0.0E0;
  Double I_ERI_H4xz_S_D2y_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_D2y_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_D2y_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_D2y_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H5y_S_D2y_S_a = 0.0E0;
  Double I_ERI_H4yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H5z_S_D2y_S_a = 0.0E0;
  Double I_ERI_H5x_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H4xy_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H4xz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H5y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H4yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H5z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_H5x_S_D2z_S_a = 0.0E0;
  Double I_ERI_H4xy_S_D2z_S_a = 0.0E0;
  Double I_ERI_H4xz_S_D2z_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_D2z_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_D2z_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_D2z_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_D2z_S_a = 0.0E0;
  Double I_ERI_H5y_S_D2z_S_a = 0.0E0;
  Double I_ERI_H4yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_D2z_S_a = 0.0E0;
  Double I_ERI_H5z_S_D2z_S_a = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_a = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_a = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_a = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_a = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_a = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_a = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_a = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_a = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_a = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_a = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_a = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_a = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_a = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_a = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_a = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_a = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_a = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_a = 0.0E0;
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
  Double I_ERI_G4x_S_F3x_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_c = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_c = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_c = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_c = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_c = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_c = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_c = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_c = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_c = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_c = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_c = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_c = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_c = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_c = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_c = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_c = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_c = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_c = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_c = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_c = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_c = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_c = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_c = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_c = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_c = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_c = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_c = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_c = 0.0E0;
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
  Double I_ERI_G4x_S_Px_S = 0.0E0;
  Double I_ERI_G3xy_S_Px_S = 0.0E0;
  Double I_ERI_G3xz_S_Px_S = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S = 0.0E0;
  Double I_ERI_G4y_S_Px_S = 0.0E0;
  Double I_ERI_G3yz_S_Px_S = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S = 0.0E0;
  Double I_ERI_G4z_S_Px_S = 0.0E0;
  Double I_ERI_G4x_S_Py_S = 0.0E0;
  Double I_ERI_G3xy_S_Py_S = 0.0E0;
  Double I_ERI_G3xz_S_Py_S = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S = 0.0E0;
  Double I_ERI_G4y_S_Py_S = 0.0E0;
  Double I_ERI_G3yz_S_Py_S = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S = 0.0E0;
  Double I_ERI_G4z_S_Py_S = 0.0E0;
  Double I_ERI_G4x_S_Pz_S = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S = 0.0E0;
  Double I_ERI_G4y_S_Pz_S = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S = 0.0E0;
  Double I_ERI_G4z_S_Pz_S = 0.0E0;
  Double I_ERI_F3x_S_Px_S = 0.0E0;
  Double I_ERI_F2xy_S_Px_S = 0.0E0;
  Double I_ERI_F2xz_S_Px_S = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S = 0.0E0;
  Double I_ERI_F3y_S_Px_S = 0.0E0;
  Double I_ERI_F2yz_S_Px_S = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S = 0.0E0;
  Double I_ERI_F3z_S_Px_S = 0.0E0;
  Double I_ERI_F3x_S_Py_S = 0.0E0;
  Double I_ERI_F2xy_S_Py_S = 0.0E0;
  Double I_ERI_F2xz_S_Py_S = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S = 0.0E0;
  Double I_ERI_F3y_S_Py_S = 0.0E0;
  Double I_ERI_F2yz_S_Py_S = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S = 0.0E0;
  Double I_ERI_F3z_S_Py_S = 0.0E0;
  Double I_ERI_F3x_S_Pz_S = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S = 0.0E0;
  Double I_ERI_F3y_S_Pz_S = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S = 0.0E0;
  Double I_ERI_F3z_S_Pz_S = 0.0E0;
  Double I_ERI_D2x_S_Px_S = 0.0E0;
  Double I_ERI_Dxy_S_Px_S = 0.0E0;
  Double I_ERI_Dxz_S_Px_S = 0.0E0;
  Double I_ERI_D2y_S_Px_S = 0.0E0;
  Double I_ERI_Dyz_S_Px_S = 0.0E0;
  Double I_ERI_D2z_S_Px_S = 0.0E0;
  Double I_ERI_D2x_S_Py_S = 0.0E0;
  Double I_ERI_Dxy_S_Py_S = 0.0E0;
  Double I_ERI_Dxz_S_Py_S = 0.0E0;
  Double I_ERI_D2y_S_Py_S = 0.0E0;
  Double I_ERI_Dyz_S_Py_S = 0.0E0;
  Double I_ERI_D2z_S_Py_S = 0.0E0;
  Double I_ERI_D2x_S_Pz_S = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S = 0.0E0;
  Double I_ERI_D2y_S_Pz_S = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S = 0.0E0;
  Double I_ERI_D2z_S_Pz_S = 0.0E0;
  Double I_ERI_H5x_S_D2x_S_b = 0.0E0;
  Double I_ERI_H4xy_S_D2x_S_b = 0.0E0;
  Double I_ERI_H4xz_S_D2x_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_D2x_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_D2x_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_D2x_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H5y_S_D2x_S_b = 0.0E0;
  Double I_ERI_H4yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H5z_S_D2x_S_b = 0.0E0;
  Double I_ERI_H5x_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H4xy_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H4xz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H5y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H4yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H5z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_H5x_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H4xy_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H4xz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H5y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H4yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H5z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_H5x_S_D2y_S_b = 0.0E0;
  Double I_ERI_H4xy_S_D2y_S_b = 0.0E0;
  Double I_ERI_H4xz_S_D2y_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_D2y_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_D2y_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_D2y_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H5y_S_D2y_S_b = 0.0E0;
  Double I_ERI_H4yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H5z_S_D2y_S_b = 0.0E0;
  Double I_ERI_H5x_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H4xy_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H4xz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H5y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H4yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H5z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_H5x_S_D2z_S_b = 0.0E0;
  Double I_ERI_H4xy_S_D2z_S_b = 0.0E0;
  Double I_ERI_H4xz_S_D2z_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_D2z_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_D2z_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_D2z_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_D2z_S_b = 0.0E0;
  Double I_ERI_H5y_S_D2z_S_b = 0.0E0;
  Double I_ERI_H4yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_D2z_S_b = 0.0E0;
  Double I_ERI_H5z_S_D2z_S_b = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_b = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_b = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_b = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_b = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_b = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_b = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_b = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_b = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_b = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_b = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_b = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_b = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_b = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_b = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_b = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_b = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_b = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_b = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_b = 0.0E0;
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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M5_vrr = PAX*I_ERI_Px_S_S_S_M5_vrr+WPX*I_ERI_Px_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_D2y_S_S_S_M5_vrr = PAY*I_ERI_Py_S_S_S_M5_vrr+WPY*I_ERI_Py_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_D2z_S_S_S_M5_vrr = PAZ*I_ERI_Pz_S_S_S_M5_vrr+WPZ*I_ERI_Pz_S_S_S_M6_vrr+oned2z*I_ERI_S_S_S_S_M5_vrr-rhod2zsq*I_ERI_S_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F2xz_S_S_S_M4_vrr = PAZ*I_ERI_D2x_S_S_S_M4_vrr+WPZ*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_Fx2y_S_S_S_M4_vrr = PAX*I_ERI_D2y_S_S_S_M4_vrr+WPX*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_Fx2z_S_S_S_M4_vrr = PAX*I_ERI_D2z_S_S_S_M4_vrr+WPX*I_ERI_D2z_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F2yz_S_S_S_M4_vrr = PAZ*I_ERI_D2y_S_S_S_M4_vrr+WPZ*I_ERI_D2y_S_S_S_M5_vrr;
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
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G2x2y_S_S_S_M3_vrr = PAY*I_ERI_F2xy_S_S_S_M3_vrr+WPY*I_ERI_F2xy_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G2xyz_S_S_S_M3_vrr = PAZ*I_ERI_F2xy_S_S_S_M3_vrr+WPZ*I_ERI_F2xy_S_S_S_M4_vrr;
      Double I_ERI_G2x2z_S_S_S_M3_vrr = PAZ*I_ERI_F2xz_S_S_S_M3_vrr+WPZ*I_ERI_F2xz_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx2yz_S_S_S_M3_vrr = PAZ*I_ERI_Fx2y_S_S_S_M3_vrr+WPZ*I_ERI_Fx2y_S_S_S_M4_vrr;
      Double I_ERI_Gxy2z_S_S_S_M3_vrr = PAY*I_ERI_Fx2z_S_S_S_M3_vrr+WPY*I_ERI_Fx2z_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_G2y2z_S_S_S_M3_vrr = PAZ*I_ERI_F2yz_S_S_S_M3_vrr+WPZ*I_ERI_F2yz_S_S_S_M4_vrr+oned2z*I_ERI_D2y_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Gy3z_S_S_S_M3_vrr = PAY*I_ERI_F3z_S_S_S_M3_vrr+WPY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M2_vrr = PAX*I_ERI_F3x_S_S_S_M2_vrr+WPX*I_ERI_F3x_S_S_S_M3_vrr+3*oned2z*I_ERI_D2x_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_S_S_M2_vrr = PAY*I_ERI_F3x_S_S_S_M2_vrr+WPY*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_S_S_M2_vrr = PAZ*I_ERI_F3x_S_S_S_M2_vrr+WPZ*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_S_S_M2_vrr = PAY*I_ERI_F2xy_S_S_S_M2_vrr+WPY*I_ERI_F2xy_S_S_S_M3_vrr+oned2z*I_ERI_D2x_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_S_S_M2_vrr = PAZ*I_ERI_F2xy_S_S_S_M2_vrr+WPZ*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_S_S_M2_vrr = PAZ*I_ERI_F2xz_S_S_S_M2_vrr+WPZ*I_ERI_F2xz_S_S_S_M3_vrr+oned2z*I_ERI_D2x_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_S_S_M2_vrr = PAX*I_ERI_F3y_S_S_S_M2_vrr+WPX*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_S_S_M2_vrr = PAZ*I_ERI_Fx2y_S_S_S_M2_vrr+WPZ*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_S_S_M2_vrr = PAY*I_ERI_Fx2z_S_S_S_M2_vrr+WPY*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_S_S_M2_vrr = PAX*I_ERI_F3z_S_S_S_M2_vrr+WPX*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_S_S_M2_vrr = PAY*I_ERI_F3y_S_S_S_M2_vrr+WPY*I_ERI_F3y_S_S_S_M3_vrr+3*oned2z*I_ERI_D2y_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_S_S_M2_vrr = PAZ*I_ERI_F3y_S_S_S_M2_vrr+WPZ*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_S_S_M2_vrr = PAZ*I_ERI_F2yz_S_S_S_M2_vrr+WPZ*I_ERI_F2yz_S_S_S_M3_vrr+oned2z*I_ERI_D2y_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_S_S_M2_vrr = PAY*I_ERI_F3z_S_S_S_M2_vrr+WPY*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_S_S_M2_vrr = PAZ*I_ERI_F3z_S_S_S_M2_vrr+WPZ*I_ERI_F3z_S_S_S_M3_vrr+3*oned2z*I_ERI_D2z_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M2_vrr = QCX*I_ERI_G4x_S_S_S_M2_vrr+WQX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Px_S_M2_vrr = QCX*I_ERI_G3xy_S_S_S_M2_vrr+WQX*I_ERI_G3xy_S_S_S_M3_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Px_S_M2_vrr = QCX*I_ERI_G3xz_S_S_S_M2_vrr+WQX*I_ERI_G3xz_S_S_S_M3_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Px_S_M2_vrr = QCX*I_ERI_G2x2y_S_S_S_M2_vrr+WQX*I_ERI_G2x2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Px_S_M2_vrr = QCX*I_ERI_G2xyz_S_S_S_M2_vrr+WQX*I_ERI_G2xyz_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Px_S_M2_vrr = QCX*I_ERI_G2x2z_S_S_S_M2_vrr+WQX*I_ERI_G2x2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Px_S_M2_vrr = QCX*I_ERI_Gx3y_S_S_S_M2_vrr+WQX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M2_vrr = QCX*I_ERI_Gx2yz_S_S_S_M2_vrr+WQX*I_ERI_Gx2yz_S_S_S_M3_vrr+oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M2_vrr = QCX*I_ERI_Gxy2z_S_S_S_M2_vrr+WQX*I_ERI_Gxy2z_S_S_S_M3_vrr+oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Px_S_M2_vrr = QCX*I_ERI_Gx3z_S_S_S_M2_vrr+WQX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2k*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Px_S_M2_vrr = QCX*I_ERI_G4y_S_S_S_M2_vrr+WQX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Px_S_M2_vrr = QCX*I_ERI_G3yz_S_S_S_M2_vrr+WQX*I_ERI_G3yz_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Px_S_M2_vrr = QCX*I_ERI_G2y2z_S_S_S_M2_vrr+WQX*I_ERI_G2y2z_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Px_S_M2_vrr = QCX*I_ERI_Gy3z_S_S_S_M2_vrr+WQX*I_ERI_Gy3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Px_S_M2_vrr = QCX*I_ERI_G4z_S_S_S_M2_vrr+WQX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_G4x_S_Py_S_M2_vrr = QCY*I_ERI_G4x_S_S_S_M2_vrr+WQY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Py_S_M2_vrr = QCY*I_ERI_G3xy_S_S_S_M2_vrr+WQY*I_ERI_G3xy_S_S_S_M3_vrr+oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Py_S_M2_vrr = QCY*I_ERI_G3xz_S_S_S_M2_vrr+WQY*I_ERI_G3xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Py_S_M2_vrr = QCY*I_ERI_G2x2y_S_S_S_M2_vrr+WQY*I_ERI_G2x2y_S_S_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Py_S_M2_vrr = QCY*I_ERI_G2xyz_S_S_S_M2_vrr+WQY*I_ERI_G2xyz_S_S_S_M3_vrr+oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Py_S_M2_vrr = QCY*I_ERI_G2x2z_S_S_S_M2_vrr+WQY*I_ERI_G2x2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Py_S_M2_vrr = QCY*I_ERI_Gx3y_S_S_S_M2_vrr+WQY*I_ERI_Gx3y_S_S_S_M3_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M2_vrr = QCY*I_ERI_Gx2yz_S_S_S_M2_vrr+WQY*I_ERI_Gx2yz_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M2_vrr = QCY*I_ERI_Gxy2z_S_S_S_M2_vrr+WQY*I_ERI_Gxy2z_S_S_S_M3_vrr+oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Py_S_M2_vrr = QCY*I_ERI_Gx3z_S_S_S_M2_vrr+WQY*I_ERI_Gx3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Py_S_M2_vrr = QCY*I_ERI_G4y_S_S_S_M2_vrr+WQY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Py_S_M2_vrr = QCY*I_ERI_G3yz_S_S_S_M2_vrr+WQY*I_ERI_G3yz_S_S_S_M3_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Py_S_M2_vrr = QCY*I_ERI_G2y2z_S_S_S_M2_vrr+WQY*I_ERI_G2y2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Py_S_M2_vrr = QCY*I_ERI_Gy3z_S_S_S_M2_vrr+WQY*I_ERI_Gy3z_S_S_S_M3_vrr+oned2k*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Py_S_M2_vrr = QCY*I_ERI_G4z_S_S_S_M2_vrr+WQY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_G4x_S_Pz_S_M2_vrr = QCZ*I_ERI_G4x_S_S_S_M2_vrr+WQZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Pz_S_M2_vrr = QCZ*I_ERI_G3xy_S_S_S_M2_vrr+WQZ*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Pz_S_M2_vrr = QCZ*I_ERI_G3xz_S_S_S_M2_vrr+WQZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M2_vrr = QCZ*I_ERI_G2x2y_S_S_S_M2_vrr+WQZ*I_ERI_G2x2y_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M2_vrr = QCZ*I_ERI_G2xyz_S_S_S_M2_vrr+WQZ*I_ERI_G2xyz_S_S_S_M3_vrr+oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M2_vrr = QCZ*I_ERI_G2x2z_S_S_S_M2_vrr+WQZ*I_ERI_G2x2z_S_S_S_M3_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx3y_S_S_S_M2_vrr+WQZ*I_ERI_Gx3y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx2yz_S_S_S_M2_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M3_vrr+oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gxy2z_S_S_S_M2_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx3z_S_S_S_M2_vrr+WQZ*I_ERI_Gx3z_S_S_S_M3_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Pz_S_M2_vrr = QCZ*I_ERI_G4y_S_S_S_M2_vrr+WQZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Pz_S_M2_vrr = QCZ*I_ERI_G3yz_S_S_S_M2_vrr+WQZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M2_vrr = QCZ*I_ERI_G2y2z_S_S_S_M2_vrr+WQZ*I_ERI_G2y2z_S_S_S_M3_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gy3z_S_S_S_M2_vrr+WQZ*I_ERI_Gy3z_S_S_S_M3_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Pz_S_M2_vrr = QCZ*I_ERI_G4z_S_S_S_M2_vrr+WQZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2k*I_ERI_F3z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_P_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M1_vrr = QCX*I_ERI_Px_S_Px_S_M1_vrr+WQX*I_ERI_Px_S_Px_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_D2x_S_M1_vrr = QCX*I_ERI_Py_S_Px_S_M1_vrr+WQX*I_ERI_Py_S_Px_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_D2x_S_M1_vrr = QCX*I_ERI_Pz_S_Px_S_M1_vrr+WQX*I_ERI_Pz_S_Px_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_Px_S_Dxy_S_M1_vrr = QCY*I_ERI_Px_S_Px_S_M1_vrr+WQY*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_Dxy_S_M1_vrr = QCY*I_ERI_Py_S_Px_S_M1_vrr+WQY*I_ERI_Py_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxy_S_M1_vrr = QCY*I_ERI_Pz_S_Px_S_M1_vrr+WQY*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_Dxz_S_M1_vrr = QCZ*I_ERI_Px_S_Px_S_M1_vrr+WQZ*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_Dxz_S_M1_vrr = QCZ*I_ERI_Py_S_Px_S_M1_vrr+WQZ*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxz_S_M1_vrr = QCZ*I_ERI_Pz_S_Px_S_M1_vrr+WQZ*I_ERI_Pz_S_Px_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_D2y_S_M1_vrr = QCY*I_ERI_Px_S_Py_S_M1_vrr+WQY*I_ERI_Px_S_Py_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Py_S_D2y_S_M1_vrr = QCY*I_ERI_Py_S_Py_S_M1_vrr+WQY*I_ERI_Py_S_Py_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_D2y_S_M1_vrr = QCY*I_ERI_Pz_S_Py_S_M1_vrr+WQY*I_ERI_Pz_S_Py_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr;
      Double I_ERI_Px_S_Dyz_S_M1_vrr = QCZ*I_ERI_Px_S_Py_S_M1_vrr+WQZ*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Dyz_S_M1_vrr = QCZ*I_ERI_Py_S_Py_S_M1_vrr+WQZ*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_Dyz_S_M1_vrr = QCZ*I_ERI_Pz_S_Py_S_M1_vrr+WQZ*I_ERI_Pz_S_Py_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_D2z_S_M1_vrr = QCZ*I_ERI_Px_S_Pz_S_M1_vrr+WQZ*I_ERI_Px_S_Pz_S_M2_vrr+oned2e*I_ERI_Px_S_S_S_M1_vrr-rhod2esq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Py_S_D2z_S_M1_vrr = QCZ*I_ERI_Py_S_Pz_S_M1_vrr+WQZ*I_ERI_Py_S_Pz_S_M2_vrr+oned2e*I_ERI_Py_S_S_S_M1_vrr-rhod2esq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_Pz_S_D2z_S_M1_vrr = QCZ*I_ERI_Pz_S_Pz_S_M1_vrr+WQZ*I_ERI_Pz_S_Pz_S_M2_vrr+oned2e*I_ERI_Pz_S_S_S_M1_vrr-rhod2esq*I_ERI_Pz_S_S_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;

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
       * expanding position: KET1
       * code section is: VRR
       * totally 4 integrals are omitted 
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
      Double I_ERI_Dxy_S_Dxz_S_M1_vrr = QCZ*I_ERI_Dxy_S_Px_S_M1_vrr+WQZ*I_ERI_Dxy_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2y_S_Px_S_M1_vrr+WQZ*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_D2z_S_Px_S_M1_vrr+WQZ*I_ERI_D2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_D2y_S_M1_vrr = QCY*I_ERI_D2x_S_Py_S_M1_vrr+WQY*I_ERI_D2x_S_Py_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = QCY*I_ERI_Dxy_S_Py_S_M1_vrr+WQY*I_ERI_Dxy_S_Py_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = QCY*I_ERI_Dxz_S_Py_S_M1_vrr+WQY*I_ERI_Dxz_S_Py_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = QCY*I_ERI_D2y_S_Py_S_M1_vrr+WQY*I_ERI_D2y_S_Py_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = QCY*I_ERI_Dyz_S_Py_S_M1_vrr+WQY*I_ERI_Dyz_S_Py_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = QCY*I_ERI_D2z_S_Py_S_M1_vrr+WQY*I_ERI_D2z_S_Py_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_Dyz_S_M1_vrr = QCZ*I_ERI_D2x_S_Py_S_M1_vrr+WQZ*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Dyz_S_M1_vrr = QCZ*I_ERI_Dxy_S_Py_S_M1_vrr+WQZ*I_ERI_Dxy_S_Py_S_M2_vrr;
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
       * shell quartet name: SQ_ERI_G_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_M1_vrr = QCX*I_ERI_G4x_S_Px_S_M1_vrr+WQX*I_ERI_G4x_S_Px_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr+4*oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_D2x_S_M1_vrr = QCX*I_ERI_G3xy_S_Px_S_M1_vrr+WQX*I_ERI_G3xy_S_Px_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_D2x_S_M1_vrr = QCX*I_ERI_G3xz_S_Px_S_M1_vrr+WQX*I_ERI_G3xz_S_Px_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2x_S_M1_vrr = QCX*I_ERI_G2x2y_S_Px_S_M1_vrr+WQX*I_ERI_G2x2y_S_Px_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2x_S_M1_vrr = QCX*I_ERI_G2xyz_S_Px_S_M1_vrr+WQX*I_ERI_G2xyz_S_Px_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2x_S_M1_vrr = QCX*I_ERI_G2x2z_S_Px_S_M1_vrr+WQX*I_ERI_G2x2z_S_Px_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2x_S_M1_vrr = QCX*I_ERI_Gx3y_S_Px_S_M1_vrr+WQX*I_ERI_Gx3y_S_Px_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2x_S_M1_vrr = QCX*I_ERI_Gx2yz_S_Px_S_M1_vrr+WQX*I_ERI_Gx2yz_S_Px_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2x_S_M1_vrr = QCX*I_ERI_Gxy2z_S_Px_S_M1_vrr+WQX*I_ERI_Gxy2z_S_Px_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2x_S_M1_vrr = QCX*I_ERI_Gx3z_S_Px_S_M1_vrr+WQX*I_ERI_Gx3z_S_Px_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_D2x_S_M1_vrr = QCX*I_ERI_G4y_S_Px_S_M1_vrr+WQX*I_ERI_G4y_S_Px_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_D2x_S_M1_vrr = QCX*I_ERI_G3yz_S_Px_S_M1_vrr+WQX*I_ERI_G3yz_S_Px_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2x_S_M1_vrr = QCX*I_ERI_G2y2z_S_Px_S_M1_vrr+WQX*I_ERI_G2y2z_S_Px_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2x_S_M1_vrr = QCX*I_ERI_Gy3z_S_Px_S_M1_vrr+WQX*I_ERI_Gy3z_S_Px_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_D2x_S_M1_vrr = QCX*I_ERI_G4z_S_Px_S_M1_vrr+WQX*I_ERI_G4z_S_Px_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Dxy_S_M1_vrr = QCY*I_ERI_G4x_S_Px_S_M1_vrr+WQY*I_ERI_G4x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_Dxy_S_M1_vrr = QCY*I_ERI_G3xy_S_Px_S_M1_vrr+WQY*I_ERI_G3xy_S_Px_S_M2_vrr+oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Dxy_S_M1_vrr = QCY*I_ERI_G3xz_S_Px_S_M1_vrr+WQY*I_ERI_G3xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_M1_vrr = QCY*I_ERI_G2x2y_S_Px_S_M1_vrr+WQY*I_ERI_G2x2y_S_Px_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_Dxy_S_M1_vrr = QCY*I_ERI_G2xyz_S_Px_S_M1_vrr+WQY*I_ERI_G2xyz_S_Px_S_M2_vrr+oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dxy_S_M1_vrr = QCY*I_ERI_G2x2z_S_Px_S_M1_vrr+WQY*I_ERI_G2x2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx3y_S_Px_S_M1_vrr+WQY*I_ERI_Gx3y_S_Px_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx2yz_S_Px_S_M1_vrr+WQY*I_ERI_Gx2yz_S_Px_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gxy2z_S_Px_S_M1_vrr+WQY*I_ERI_Gxy2z_S_Px_S_M2_vrr+oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx3z_S_Px_S_M1_vrr+WQY*I_ERI_Gx3z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_Dxy_S_M1_vrr = QCY*I_ERI_G4y_S_Px_S_M1_vrr+WQY*I_ERI_G4y_S_Px_S_M2_vrr+4*oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Dxy_S_M1_vrr = QCY*I_ERI_G3yz_S_Px_S_M1_vrr+WQY*I_ERI_G3yz_S_Px_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dxy_S_M1_vrr = QCY*I_ERI_G2y2z_S_Px_S_M1_vrr+WQY*I_ERI_G2y2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gy3z_S_Px_S_M1_vrr+WQY*I_ERI_Gy3z_S_Px_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Dxy_S_M1_vrr = QCY*I_ERI_G4z_S_Px_S_M1_vrr+WQY*I_ERI_G4z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_Dxz_S_M1_vrr = QCZ*I_ERI_G4x_S_Px_S_M1_vrr+WQZ*I_ERI_G4x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_Dxz_S_M1_vrr = QCZ*I_ERI_G3xy_S_Px_S_M1_vrr+WQZ*I_ERI_G3xy_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Dxz_S_M1_vrr = QCZ*I_ERI_G3xz_S_Px_S_M1_vrr+WQZ*I_ERI_G3xz_S_Px_S_M2_vrr+oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dxz_S_M1_vrr = QCZ*I_ERI_G2x2y_S_Px_S_M1_vrr+WQZ*I_ERI_G2x2y_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_G2x2z_S_Px_S_M1_vrr+WQZ*I_ERI_G2x2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dxz_S_M1_vrr = QCZ*I_ERI_Gx3y_S_Px_S_M1_vrr+WQZ*I_ERI_Gx3y_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dxz_S_M1_vrr = QCZ*I_ERI_Gx3z_S_Px_S_M1_vrr+WQZ*I_ERI_Gx3z_S_Px_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_Dxz_S_M1_vrr = QCZ*I_ERI_G4y_S_Px_S_M1_vrr+WQZ*I_ERI_G4y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Dxz_S_M1_vrr = QCZ*I_ERI_G3yz_S_Px_S_M1_vrr+WQZ*I_ERI_G3yz_S_Px_S_M2_vrr+oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dxz_S_M1_vrr = QCZ*I_ERI_G2y2z_S_Px_S_M1_vrr+WQZ*I_ERI_G2y2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dxz_S_M1_vrr = QCZ*I_ERI_Gy3z_S_Px_S_M1_vrr+WQZ*I_ERI_Gy3z_S_Px_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Dxz_S_M1_vrr = QCZ*I_ERI_G4z_S_Px_S_M1_vrr+WQZ*I_ERI_G4z_S_Px_S_M2_vrr+4*oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_D2y_S_M1_vrr = QCY*I_ERI_G4x_S_Py_S_M1_vrr+WQY*I_ERI_G4x_S_Py_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_D2y_S_M1_vrr = QCY*I_ERI_G3xy_S_Py_S_M1_vrr+WQY*I_ERI_G3xy_S_Py_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G3xz_S_D2y_S_M1_vrr = QCY*I_ERI_G3xz_S_Py_S_M1_vrr+WQY*I_ERI_G3xz_S_Py_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2y_S_M1_vrr = QCY*I_ERI_G2x2y_S_Py_S_M1_vrr+WQY*I_ERI_G2x2y_S_Py_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2y_S_M1_vrr = QCY*I_ERI_G2xyz_S_Py_S_M1_vrr+WQY*I_ERI_G2xyz_S_Py_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xz_S_Py_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2y_S_M1_vrr = QCY*I_ERI_G2x2z_S_Py_S_M1_vrr+WQY*I_ERI_G2x2z_S_Py_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2y_S_M1_vrr = QCY*I_ERI_Gx3y_S_Py_S_M1_vrr+WQY*I_ERI_Gx3y_S_Py_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_Py_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2y_S_M1_vrr = QCY*I_ERI_Gx2yz_S_Py_S_M1_vrr+WQY*I_ERI_Gx2yz_S_Py_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Py_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2y_S_M1_vrr = QCY*I_ERI_Gxy2z_S_Py_S_M1_vrr+WQY*I_ERI_Gxy2z_S_Py_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fx2z_S_Py_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2y_S_M1_vrr = QCY*I_ERI_Gx3z_S_Py_S_M1_vrr+WQY*I_ERI_Gx3z_S_Py_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_D2y_S_M1_vrr = QCY*I_ERI_G4y_S_Py_S_M1_vrr+WQY*I_ERI_G4y_S_Py_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr+4*oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G3yz_S_D2y_S_M1_vrr = QCY*I_ERI_G3yz_S_Py_S_M1_vrr+WQY*I_ERI_G3yz_S_Py_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2y_S_M1_vrr = QCY*I_ERI_G2y2z_S_Py_S_M1_vrr+WQY*I_ERI_G2y2z_S_Py_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2y_S_M1_vrr = QCY*I_ERI_Gy3z_S_Py_S_M1_vrr+WQY*I_ERI_Gy3z_S_Py_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4z_S_D2y_S_M1_vrr = QCY*I_ERI_G4z_S_Py_S_M1_vrr+WQY*I_ERI_G4z_S_Py_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Dyz_S_M1_vrr = QCZ*I_ERI_G4x_S_Py_S_M1_vrr+WQZ*I_ERI_G4x_S_Py_S_M2_vrr;
      Double I_ERI_G3xy_S_Dyz_S_M1_vrr = QCZ*I_ERI_G3xy_S_Py_S_M1_vrr+WQZ*I_ERI_G3xy_S_Py_S_M2_vrr;
      Double I_ERI_G3xz_S_Dyz_S_M1_vrr = QCZ*I_ERI_G3xz_S_Py_S_M1_vrr+WQZ*I_ERI_G3xz_S_Py_S_M2_vrr+oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dyz_S_M1_vrr = QCZ*I_ERI_G2x2y_S_Py_S_M1_vrr+WQZ*I_ERI_G2x2y_S_Py_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dyz_S_M1_vrr = QCZ*I_ERI_G2x2z_S_Py_S_M1_vrr+WQZ*I_ERI_G2x2z_S_Py_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Py_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dyz_S_M1_vrr = QCZ*I_ERI_Gx3y_S_Py_S_M1_vrr+WQZ*I_ERI_Gx3y_S_Py_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dyz_S_M1_vrr = QCZ*I_ERI_Gx3z_S_Py_S_M1_vrr+WQZ*I_ERI_Gx3z_S_Py_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_Py_S_M2_vrr;
      Double I_ERI_G4y_S_Dyz_S_M1_vrr = QCZ*I_ERI_G4y_S_Py_S_M1_vrr+WQZ*I_ERI_G4y_S_Py_S_M2_vrr;
      Double I_ERI_G3yz_S_Dyz_S_M1_vrr = QCZ*I_ERI_G3yz_S_Py_S_M1_vrr+WQZ*I_ERI_G3yz_S_Py_S_M2_vrr+oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dyz_S_M1_vrr = QCZ*I_ERI_G2y2z_S_Py_S_M1_vrr+WQZ*I_ERI_G2y2z_S_Py_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dyz_S_M1_vrr = QCZ*I_ERI_Gy3z_S_Py_S_M1_vrr+WQZ*I_ERI_Gy3z_S_Py_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_Py_S_M2_vrr;
      Double I_ERI_G4z_S_Dyz_S_M1_vrr = QCZ*I_ERI_G4z_S_Py_S_M1_vrr+WQZ*I_ERI_G4z_S_Py_S_M2_vrr+4*oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4x_S_D2z_S_M1_vrr = QCZ*I_ERI_G4x_S_Pz_S_M1_vrr+WQZ*I_ERI_G4x_S_Pz_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_D2z_S_M1_vrr = QCZ*I_ERI_G3xy_S_Pz_S_M1_vrr+WQZ*I_ERI_G3xy_S_Pz_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_D2z_S_M1_vrr = QCZ*I_ERI_G3xz_S_Pz_S_M1_vrr+WQZ*I_ERI_G3xz_S_Pz_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2z_S_M1_vrr = QCZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+WQZ*I_ERI_G2x2y_S_Pz_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2z_S_M1_vrr = QCZ*I_ERI_G2xyz_S_Pz_S_M1_vrr+WQZ*I_ERI_G2xyz_S_Pz_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2z_S_M1_vrr = QCZ*I_ERI_G2x2z_S_Pz_S_M1_vrr+WQZ*I_ERI_G2x2z_S_Pz_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx3y_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx2yz_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx2yz_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_Fx2y_S_Pz_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gxy2z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gxy2z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx3z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx3z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_Pz_S_M2_vrr;
      Double I_ERI_G4y_S_D2z_S_M1_vrr = QCZ*I_ERI_G4y_S_Pz_S_M1_vrr+WQZ*I_ERI_G4y_S_Pz_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_D2z_S_M1_vrr = QCZ*I_ERI_G3yz_S_Pz_S_M1_vrr+WQZ*I_ERI_G3yz_S_Pz_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2z_S_M1_vrr = QCZ*I_ERI_G2y2z_S_Pz_S_M1_vrr+WQZ*I_ERI_G2y2z_S_Pz_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Pz_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gy3z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gy3z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_Pz_S_M2_vrr;
      Double I_ERI_G4z_S_D2z_S_M1_vrr = QCZ*I_ERI_G4z_S_Pz_S_M1_vrr+WQZ*I_ERI_G4z_S_Pz_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr+4*oned2k*I_ERI_F3z_S_Pz_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_F_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_D2x_S += I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S += I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S += I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S += I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S += I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S += I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S += I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S += I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S += I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S += I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S += I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S += I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S += I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S += I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S += I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S += I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S += I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S += I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S += I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S += I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S += I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S += I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S += I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S += I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S += I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S += I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S += I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S += I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S += I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S += I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S += I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S += I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S += I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S += I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S += I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S += I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S += I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S += I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S += I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S += I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S += I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S += I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S += I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S += I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S += I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S += I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S += I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S += I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S += I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S += I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S += I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S += I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S += I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S += I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S += I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S += I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S += I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S += I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S += I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S += I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_D2x_S += I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S += I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S += I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S += I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S += I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S += I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S += I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S += I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S += I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S += I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S += I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S += I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S += I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S += I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S += I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S += I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S += I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S += I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S += I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S += I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S += I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S += I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S += I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S += I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S += I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S += I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S += I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S += I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S += I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S += I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S += I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S += I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S += I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S += I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S += I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S += I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_D_S_a_coefs = alpha;
      I_ERI_H5x_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_D2x_S_vrr;
      I_ERI_H4xy_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_D2x_S_vrr;
      I_ERI_H4xz_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_D2x_S_vrr;
      I_ERI_H3x2y_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_D2x_S_vrr;
      I_ERI_H3xyz_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_D2x_S_vrr;
      I_ERI_H3x2z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_D2x_S_vrr;
      I_ERI_H2x3y_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_D2x_S_vrr;
      I_ERI_H2x2yz_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_D2x_S_vrr;
      I_ERI_H2xy2z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_D2x_S_vrr;
      I_ERI_H2x3z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_D2x_S_vrr;
      I_ERI_Hx4y_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_D2x_S_vrr;
      I_ERI_Hx3yz_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_D2x_S_vrr;
      I_ERI_Hx2y2z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_D2x_S_vrr;
      I_ERI_Hxy3z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_D2x_S_vrr;
      I_ERI_Hx4z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_D2x_S_vrr;
      I_ERI_H5y_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_D2x_S_vrr;
      I_ERI_H4yz_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_D2x_S_vrr;
      I_ERI_H3y2z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_D2x_S_vrr;
      I_ERI_H2y3z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_D2x_S_vrr;
      I_ERI_Hy4z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_D2x_S_vrr;
      I_ERI_H5z_S_D2x_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_D2x_S_vrr;
      I_ERI_H5x_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_Dxy_S_vrr;
      I_ERI_H4xy_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_Dxy_S_vrr;
      I_ERI_H4xz_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_Dxy_S_vrr;
      I_ERI_H3x2y_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_Dxy_S_vrr;
      I_ERI_H3xyz_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_Dxy_S_vrr;
      I_ERI_H3x2z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_Dxy_S_vrr;
      I_ERI_H2x3y_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_Dxy_S_vrr;
      I_ERI_H2x2yz_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_Dxy_S_vrr;
      I_ERI_H2xy2z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_Dxy_S_vrr;
      I_ERI_H2x3z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_Dxy_S_vrr;
      I_ERI_Hx4y_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_Dxy_S_vrr;
      I_ERI_Hx3yz_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_Dxy_S_vrr;
      I_ERI_Hx2y2z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_Dxy_S_vrr;
      I_ERI_Hxy3z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_Dxy_S_vrr;
      I_ERI_Hx4z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_Dxy_S_vrr;
      I_ERI_H5y_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_Dxy_S_vrr;
      I_ERI_H4yz_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_Dxy_S_vrr;
      I_ERI_H3y2z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_Dxy_S_vrr;
      I_ERI_H2y3z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_Dxy_S_vrr;
      I_ERI_Hy4z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_Dxy_S_vrr;
      I_ERI_H5z_S_Dxy_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_Dxy_S_vrr;
      I_ERI_H5x_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_Dxz_S_vrr;
      I_ERI_H4xy_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_Dxz_S_vrr;
      I_ERI_H4xz_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_Dxz_S_vrr;
      I_ERI_H3x2y_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_Dxz_S_vrr;
      I_ERI_H3xyz_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_Dxz_S_vrr;
      I_ERI_H3x2z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_Dxz_S_vrr;
      I_ERI_H2x3y_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_Dxz_S_vrr;
      I_ERI_H2x2yz_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_Dxz_S_vrr;
      I_ERI_H2xy2z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_Dxz_S_vrr;
      I_ERI_H2x3z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_Dxz_S_vrr;
      I_ERI_Hx4y_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_Dxz_S_vrr;
      I_ERI_Hx3yz_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_Dxz_S_vrr;
      I_ERI_Hx2y2z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_Dxz_S_vrr;
      I_ERI_Hxy3z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_Dxz_S_vrr;
      I_ERI_Hx4z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_Dxz_S_vrr;
      I_ERI_H5y_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_Dxz_S_vrr;
      I_ERI_H4yz_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_Dxz_S_vrr;
      I_ERI_H3y2z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_Dxz_S_vrr;
      I_ERI_H2y3z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_Dxz_S_vrr;
      I_ERI_Hy4z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_Dxz_S_vrr;
      I_ERI_H5z_S_Dxz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_Dxz_S_vrr;
      I_ERI_H5x_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_D2y_S_vrr;
      I_ERI_H4xy_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_D2y_S_vrr;
      I_ERI_H4xz_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_D2y_S_vrr;
      I_ERI_H3x2y_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_D2y_S_vrr;
      I_ERI_H3xyz_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_D2y_S_vrr;
      I_ERI_H3x2z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_D2y_S_vrr;
      I_ERI_H2x3y_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_D2y_S_vrr;
      I_ERI_H2x2yz_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_D2y_S_vrr;
      I_ERI_H2xy2z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_D2y_S_vrr;
      I_ERI_H2x3z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_D2y_S_vrr;
      I_ERI_Hx4y_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_D2y_S_vrr;
      I_ERI_Hx3yz_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_D2y_S_vrr;
      I_ERI_Hx2y2z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_D2y_S_vrr;
      I_ERI_Hxy3z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_D2y_S_vrr;
      I_ERI_Hx4z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_D2y_S_vrr;
      I_ERI_H5y_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_D2y_S_vrr;
      I_ERI_H4yz_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_D2y_S_vrr;
      I_ERI_H3y2z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_D2y_S_vrr;
      I_ERI_H2y3z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_D2y_S_vrr;
      I_ERI_Hy4z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_D2y_S_vrr;
      I_ERI_H5z_S_D2y_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_D2y_S_vrr;
      I_ERI_H5x_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_Dyz_S_vrr;
      I_ERI_H4xy_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_Dyz_S_vrr;
      I_ERI_H4xz_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_Dyz_S_vrr;
      I_ERI_H3x2y_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_Dyz_S_vrr;
      I_ERI_H3xyz_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_Dyz_S_vrr;
      I_ERI_H3x2z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_Dyz_S_vrr;
      I_ERI_H2x3y_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_Dyz_S_vrr;
      I_ERI_H2x2yz_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_Dyz_S_vrr;
      I_ERI_H2xy2z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_Dyz_S_vrr;
      I_ERI_H2x3z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_Dyz_S_vrr;
      I_ERI_Hx4y_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_Dyz_S_vrr;
      I_ERI_Hx3yz_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_Dyz_S_vrr;
      I_ERI_Hx2y2z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_Dyz_S_vrr;
      I_ERI_Hxy3z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_Dyz_S_vrr;
      I_ERI_Hx4z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_Dyz_S_vrr;
      I_ERI_H5y_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_Dyz_S_vrr;
      I_ERI_H4yz_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_Dyz_S_vrr;
      I_ERI_H3y2z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_Dyz_S_vrr;
      I_ERI_H2y3z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_Dyz_S_vrr;
      I_ERI_Hy4z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_Dyz_S_vrr;
      I_ERI_H5z_S_Dyz_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_Dyz_S_vrr;
      I_ERI_H5x_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5x_S_D2z_S_vrr;
      I_ERI_H4xy_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xy_S_D2z_S_vrr;
      I_ERI_H4xz_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4xz_S_D2z_S_vrr;
      I_ERI_H3x2y_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2y_S_D2z_S_vrr;
      I_ERI_H3xyz_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3xyz_S_D2z_S_vrr;
      I_ERI_H3x2z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3x2z_S_D2z_S_vrr;
      I_ERI_H2x3y_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3y_S_D2z_S_vrr;
      I_ERI_H2x2yz_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x2yz_S_D2z_S_vrr;
      I_ERI_H2xy2z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2xy2z_S_D2z_S_vrr;
      I_ERI_H2x3z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2x3z_S_D2z_S_vrr;
      I_ERI_Hx4y_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4y_S_D2z_S_vrr;
      I_ERI_Hx3yz_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx3yz_S_D2z_S_vrr;
      I_ERI_Hx2y2z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx2y2z_S_D2z_S_vrr;
      I_ERI_Hxy3z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hxy3z_S_D2z_S_vrr;
      I_ERI_Hx4z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hx4z_S_D2z_S_vrr;
      I_ERI_H5y_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5y_S_D2z_S_vrr;
      I_ERI_H4yz_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H4yz_S_D2z_S_vrr;
      I_ERI_H3y2z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H3y2z_S_D2z_S_vrr;
      I_ERI_H2y3z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H2y3z_S_D2z_S_vrr;
      I_ERI_Hy4z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_Hy4z_S_D2z_S_vrr;
      I_ERI_H5z_S_D2z_S_a += SQ_ERI_H_S_D_S_a_coefs*I_ERI_H5z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_a_coefs = alpha;
      I_ERI_G4x_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_a += SQ_ERI_G_S_D_S_a_coefs*I_ERI_G4z_S_D2z_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_F_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_c_coefs = gamma;
      I_ERI_G4x_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_c += SQ_ERI_G_S_F_S_c_coefs*I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_c_coefs = gamma;
      I_ERI_F3x_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_c += SQ_ERI_F_S_F_S_c_coefs*I_ERI_F3z_S_F3z_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_Px_S += I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S += I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S += I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S += I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S += I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S += I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S += I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S += I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S += I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S += I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S += I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S += I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S += I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S += I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S += I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S += I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S += I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S += I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S += I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S += I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S += I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S += I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S += I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S += I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S += I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S += I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S += I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S += I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S += I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S += I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S += I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S += I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S += I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S += I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S += I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S += I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S += I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S += I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S += I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S += I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S += I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S += I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S += I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S += I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S += I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_Px_S += I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S += I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S += I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S += I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S += I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S += I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S += I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S += I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S += I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S += I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S += I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S += I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S += I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S += I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S += I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S += I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S += I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S += I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S += I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S += I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S += I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S += I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S += I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S += I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S += I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S += I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S += I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S += I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S += I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S += I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_Px_S += I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S += I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S += I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S += I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S += I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S += I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S += I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S += I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S += I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S += I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S += I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S += I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S += I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S += I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S += I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S += I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S += I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S += I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_D_S_b_coefs = beta;
      I_ERI_H5x_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_D2x_S_vrr;
      I_ERI_H4xy_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_D2x_S_vrr;
      I_ERI_H4xz_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_D2x_S_vrr;
      I_ERI_H3x2y_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_D2x_S_vrr;
      I_ERI_H3xyz_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_D2x_S_vrr;
      I_ERI_H3x2z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_D2x_S_vrr;
      I_ERI_H2x3y_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_D2x_S_vrr;
      I_ERI_H2x2yz_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_D2x_S_vrr;
      I_ERI_H2xy2z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_D2x_S_vrr;
      I_ERI_H2x3z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_D2x_S_vrr;
      I_ERI_Hx4y_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_D2x_S_vrr;
      I_ERI_Hx3yz_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_D2x_S_vrr;
      I_ERI_Hx2y2z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_D2x_S_vrr;
      I_ERI_Hxy3z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_D2x_S_vrr;
      I_ERI_Hx4z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_D2x_S_vrr;
      I_ERI_H5y_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_D2x_S_vrr;
      I_ERI_H4yz_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_D2x_S_vrr;
      I_ERI_H3y2z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_D2x_S_vrr;
      I_ERI_H2y3z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_D2x_S_vrr;
      I_ERI_Hy4z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_D2x_S_vrr;
      I_ERI_H5z_S_D2x_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_D2x_S_vrr;
      I_ERI_H5x_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_Dxy_S_vrr;
      I_ERI_H4xy_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_Dxy_S_vrr;
      I_ERI_H4xz_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_Dxy_S_vrr;
      I_ERI_H3x2y_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_Dxy_S_vrr;
      I_ERI_H3xyz_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_Dxy_S_vrr;
      I_ERI_H3x2z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_Dxy_S_vrr;
      I_ERI_H2x3y_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_Dxy_S_vrr;
      I_ERI_H2x2yz_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_Dxy_S_vrr;
      I_ERI_H2xy2z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_Dxy_S_vrr;
      I_ERI_H2x3z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_Dxy_S_vrr;
      I_ERI_Hx4y_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_Dxy_S_vrr;
      I_ERI_Hx3yz_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_Dxy_S_vrr;
      I_ERI_Hx2y2z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_Dxy_S_vrr;
      I_ERI_Hxy3z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_Dxy_S_vrr;
      I_ERI_Hx4z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_Dxy_S_vrr;
      I_ERI_H5y_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_Dxy_S_vrr;
      I_ERI_H4yz_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_Dxy_S_vrr;
      I_ERI_H3y2z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_Dxy_S_vrr;
      I_ERI_H2y3z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_Dxy_S_vrr;
      I_ERI_Hy4z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_Dxy_S_vrr;
      I_ERI_H5z_S_Dxy_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_Dxy_S_vrr;
      I_ERI_H5x_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_Dxz_S_vrr;
      I_ERI_H4xy_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_Dxz_S_vrr;
      I_ERI_H4xz_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_Dxz_S_vrr;
      I_ERI_H3x2y_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_Dxz_S_vrr;
      I_ERI_H3xyz_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_Dxz_S_vrr;
      I_ERI_H3x2z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_Dxz_S_vrr;
      I_ERI_H2x3y_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_Dxz_S_vrr;
      I_ERI_H2x2yz_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_Dxz_S_vrr;
      I_ERI_H2xy2z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_Dxz_S_vrr;
      I_ERI_H2x3z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_Dxz_S_vrr;
      I_ERI_Hx4y_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_Dxz_S_vrr;
      I_ERI_Hx3yz_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_Dxz_S_vrr;
      I_ERI_Hx2y2z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_Dxz_S_vrr;
      I_ERI_Hxy3z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_Dxz_S_vrr;
      I_ERI_Hx4z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_Dxz_S_vrr;
      I_ERI_H5y_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_Dxz_S_vrr;
      I_ERI_H4yz_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_Dxz_S_vrr;
      I_ERI_H3y2z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_Dxz_S_vrr;
      I_ERI_H2y3z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_Dxz_S_vrr;
      I_ERI_Hy4z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_Dxz_S_vrr;
      I_ERI_H5z_S_Dxz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_Dxz_S_vrr;
      I_ERI_H5x_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_D2y_S_vrr;
      I_ERI_H4xy_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_D2y_S_vrr;
      I_ERI_H4xz_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_D2y_S_vrr;
      I_ERI_H3x2y_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_D2y_S_vrr;
      I_ERI_H3xyz_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_D2y_S_vrr;
      I_ERI_H3x2z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_D2y_S_vrr;
      I_ERI_H2x3y_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_D2y_S_vrr;
      I_ERI_H2x2yz_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_D2y_S_vrr;
      I_ERI_H2xy2z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_D2y_S_vrr;
      I_ERI_H2x3z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_D2y_S_vrr;
      I_ERI_Hx4y_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_D2y_S_vrr;
      I_ERI_Hx3yz_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_D2y_S_vrr;
      I_ERI_Hx2y2z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_D2y_S_vrr;
      I_ERI_Hxy3z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_D2y_S_vrr;
      I_ERI_Hx4z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_D2y_S_vrr;
      I_ERI_H5y_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_D2y_S_vrr;
      I_ERI_H4yz_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_D2y_S_vrr;
      I_ERI_H3y2z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_D2y_S_vrr;
      I_ERI_H2y3z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_D2y_S_vrr;
      I_ERI_Hy4z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_D2y_S_vrr;
      I_ERI_H5z_S_D2y_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_D2y_S_vrr;
      I_ERI_H5x_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_Dyz_S_vrr;
      I_ERI_H4xy_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_Dyz_S_vrr;
      I_ERI_H4xz_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_Dyz_S_vrr;
      I_ERI_H3x2y_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_Dyz_S_vrr;
      I_ERI_H3xyz_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_Dyz_S_vrr;
      I_ERI_H3x2z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_Dyz_S_vrr;
      I_ERI_H2x3y_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_Dyz_S_vrr;
      I_ERI_H2x2yz_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_Dyz_S_vrr;
      I_ERI_H2xy2z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_Dyz_S_vrr;
      I_ERI_H2x3z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_Dyz_S_vrr;
      I_ERI_Hx4y_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_Dyz_S_vrr;
      I_ERI_Hx3yz_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_Dyz_S_vrr;
      I_ERI_Hx2y2z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_Dyz_S_vrr;
      I_ERI_Hxy3z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_Dyz_S_vrr;
      I_ERI_Hx4z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_Dyz_S_vrr;
      I_ERI_H5y_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_Dyz_S_vrr;
      I_ERI_H4yz_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_Dyz_S_vrr;
      I_ERI_H3y2z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_Dyz_S_vrr;
      I_ERI_H2y3z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_Dyz_S_vrr;
      I_ERI_Hy4z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_Dyz_S_vrr;
      I_ERI_H5z_S_Dyz_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_Dyz_S_vrr;
      I_ERI_H5x_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5x_S_D2z_S_vrr;
      I_ERI_H4xy_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xy_S_D2z_S_vrr;
      I_ERI_H4xz_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4xz_S_D2z_S_vrr;
      I_ERI_H3x2y_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2y_S_D2z_S_vrr;
      I_ERI_H3xyz_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3xyz_S_D2z_S_vrr;
      I_ERI_H3x2z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3x2z_S_D2z_S_vrr;
      I_ERI_H2x3y_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3y_S_D2z_S_vrr;
      I_ERI_H2x2yz_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x2yz_S_D2z_S_vrr;
      I_ERI_H2xy2z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2xy2z_S_D2z_S_vrr;
      I_ERI_H2x3z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2x3z_S_D2z_S_vrr;
      I_ERI_Hx4y_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4y_S_D2z_S_vrr;
      I_ERI_Hx3yz_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx3yz_S_D2z_S_vrr;
      I_ERI_Hx2y2z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx2y2z_S_D2z_S_vrr;
      I_ERI_Hxy3z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hxy3z_S_D2z_S_vrr;
      I_ERI_Hx4z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hx4z_S_D2z_S_vrr;
      I_ERI_H5y_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5y_S_D2z_S_vrr;
      I_ERI_H4yz_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H4yz_S_D2z_S_vrr;
      I_ERI_H3y2z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H3y2z_S_D2z_S_vrr;
      I_ERI_H2y3z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H2y3z_S_D2z_S_vrr;
      I_ERI_Hy4z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_Hy4z_S_D2z_S_vrr;
      I_ERI_H5z_S_D2z_S_b += SQ_ERI_H_S_D_S_b_coefs*I_ERI_H5z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_b_coefs = beta;
      I_ERI_G4x_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_b += SQ_ERI_G_S_D_S_b_coefs*I_ERI_G4z_S_D2z_S_vrr;

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
   * shell quartet name: SQ_ERI_P_P_D_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S
   * RHS shell quartet name: SQ_ERI_P_S_D_S
   ************************************************************/
  Double I_ERI_Px_Px_D2x_S = I_ERI_D2x_S_D2x_S+ABX*I_ERI_Px_S_D2x_S;
  Double I_ERI_Py_Px_D2x_S = I_ERI_Dxy_S_D2x_S+ABX*I_ERI_Py_S_D2x_S;
  Double I_ERI_Pz_Px_D2x_S = I_ERI_Dxz_S_D2x_S+ABX*I_ERI_Pz_S_D2x_S;
  Double I_ERI_Px_Py_D2x_S = I_ERI_Dxy_S_D2x_S+ABY*I_ERI_Px_S_D2x_S;
  Double I_ERI_Py_Py_D2x_S = I_ERI_D2y_S_D2x_S+ABY*I_ERI_Py_S_D2x_S;
  Double I_ERI_Pz_Py_D2x_S = I_ERI_Dyz_S_D2x_S+ABY*I_ERI_Pz_S_D2x_S;
  Double I_ERI_Px_Pz_D2x_S = I_ERI_Dxz_S_D2x_S+ABZ*I_ERI_Px_S_D2x_S;
  Double I_ERI_Py_Pz_D2x_S = I_ERI_Dyz_S_D2x_S+ABZ*I_ERI_Py_S_D2x_S;
  Double I_ERI_Pz_Pz_D2x_S = I_ERI_D2z_S_D2x_S+ABZ*I_ERI_Pz_S_D2x_S;
  Double I_ERI_Px_Px_Dxy_S = I_ERI_D2x_S_Dxy_S+ABX*I_ERI_Px_S_Dxy_S;
  Double I_ERI_Py_Px_Dxy_S = I_ERI_Dxy_S_Dxy_S+ABX*I_ERI_Py_S_Dxy_S;
  Double I_ERI_Pz_Px_Dxy_S = I_ERI_Dxz_S_Dxy_S+ABX*I_ERI_Pz_S_Dxy_S;
  Double I_ERI_Px_Py_Dxy_S = I_ERI_Dxy_S_Dxy_S+ABY*I_ERI_Px_S_Dxy_S;
  Double I_ERI_Py_Py_Dxy_S = I_ERI_D2y_S_Dxy_S+ABY*I_ERI_Py_S_Dxy_S;
  Double I_ERI_Pz_Py_Dxy_S = I_ERI_Dyz_S_Dxy_S+ABY*I_ERI_Pz_S_Dxy_S;
  Double I_ERI_Px_Pz_Dxy_S = I_ERI_Dxz_S_Dxy_S+ABZ*I_ERI_Px_S_Dxy_S;
  Double I_ERI_Py_Pz_Dxy_S = I_ERI_Dyz_S_Dxy_S+ABZ*I_ERI_Py_S_Dxy_S;
  Double I_ERI_Pz_Pz_Dxy_S = I_ERI_D2z_S_Dxy_S+ABZ*I_ERI_Pz_S_Dxy_S;
  Double I_ERI_Px_Px_Dxz_S = I_ERI_D2x_S_Dxz_S+ABX*I_ERI_Px_S_Dxz_S;
  Double I_ERI_Py_Px_Dxz_S = I_ERI_Dxy_S_Dxz_S+ABX*I_ERI_Py_S_Dxz_S;
  Double I_ERI_Pz_Px_Dxz_S = I_ERI_Dxz_S_Dxz_S+ABX*I_ERI_Pz_S_Dxz_S;
  Double I_ERI_Px_Py_Dxz_S = I_ERI_Dxy_S_Dxz_S+ABY*I_ERI_Px_S_Dxz_S;
  Double I_ERI_Py_Py_Dxz_S = I_ERI_D2y_S_Dxz_S+ABY*I_ERI_Py_S_Dxz_S;
  Double I_ERI_Pz_Py_Dxz_S = I_ERI_Dyz_S_Dxz_S+ABY*I_ERI_Pz_S_Dxz_S;
  Double I_ERI_Px_Pz_Dxz_S = I_ERI_Dxz_S_Dxz_S+ABZ*I_ERI_Px_S_Dxz_S;
  Double I_ERI_Py_Pz_Dxz_S = I_ERI_Dyz_S_Dxz_S+ABZ*I_ERI_Py_S_Dxz_S;
  Double I_ERI_Pz_Pz_Dxz_S = I_ERI_D2z_S_Dxz_S+ABZ*I_ERI_Pz_S_Dxz_S;
  Double I_ERI_Px_Px_D2y_S = I_ERI_D2x_S_D2y_S+ABX*I_ERI_Px_S_D2y_S;
  Double I_ERI_Py_Px_D2y_S = I_ERI_Dxy_S_D2y_S+ABX*I_ERI_Py_S_D2y_S;
  Double I_ERI_Pz_Px_D2y_S = I_ERI_Dxz_S_D2y_S+ABX*I_ERI_Pz_S_D2y_S;
  Double I_ERI_Px_Py_D2y_S = I_ERI_Dxy_S_D2y_S+ABY*I_ERI_Px_S_D2y_S;
  Double I_ERI_Py_Py_D2y_S = I_ERI_D2y_S_D2y_S+ABY*I_ERI_Py_S_D2y_S;
  Double I_ERI_Pz_Py_D2y_S = I_ERI_Dyz_S_D2y_S+ABY*I_ERI_Pz_S_D2y_S;
  Double I_ERI_Px_Pz_D2y_S = I_ERI_Dxz_S_D2y_S+ABZ*I_ERI_Px_S_D2y_S;
  Double I_ERI_Py_Pz_D2y_S = I_ERI_Dyz_S_D2y_S+ABZ*I_ERI_Py_S_D2y_S;
  Double I_ERI_Pz_Pz_D2y_S = I_ERI_D2z_S_D2y_S+ABZ*I_ERI_Pz_S_D2y_S;
  Double I_ERI_Px_Px_Dyz_S = I_ERI_D2x_S_Dyz_S+ABX*I_ERI_Px_S_Dyz_S;
  Double I_ERI_Py_Px_Dyz_S = I_ERI_Dxy_S_Dyz_S+ABX*I_ERI_Py_S_Dyz_S;
  Double I_ERI_Pz_Px_Dyz_S = I_ERI_Dxz_S_Dyz_S+ABX*I_ERI_Pz_S_Dyz_S;
  Double I_ERI_Px_Py_Dyz_S = I_ERI_Dxy_S_Dyz_S+ABY*I_ERI_Px_S_Dyz_S;
  Double I_ERI_Py_Py_Dyz_S = I_ERI_D2y_S_Dyz_S+ABY*I_ERI_Py_S_Dyz_S;
  Double I_ERI_Pz_Py_Dyz_S = I_ERI_Dyz_S_Dyz_S+ABY*I_ERI_Pz_S_Dyz_S;
  Double I_ERI_Px_Pz_Dyz_S = I_ERI_Dxz_S_Dyz_S+ABZ*I_ERI_Px_S_Dyz_S;
  Double I_ERI_Py_Pz_Dyz_S = I_ERI_Dyz_S_Dyz_S+ABZ*I_ERI_Py_S_Dyz_S;
  Double I_ERI_Pz_Pz_Dyz_S = I_ERI_D2z_S_Dyz_S+ABZ*I_ERI_Pz_S_Dyz_S;
  Double I_ERI_Px_Px_D2z_S = I_ERI_D2x_S_D2z_S+ABX*I_ERI_Px_S_D2z_S;
  Double I_ERI_Py_Px_D2z_S = I_ERI_Dxy_S_D2z_S+ABX*I_ERI_Py_S_D2z_S;
  Double I_ERI_Pz_Px_D2z_S = I_ERI_Dxz_S_D2z_S+ABX*I_ERI_Pz_S_D2z_S;
  Double I_ERI_Px_Py_D2z_S = I_ERI_Dxy_S_D2z_S+ABY*I_ERI_Px_S_D2z_S;
  Double I_ERI_Py_Py_D2z_S = I_ERI_D2y_S_D2z_S+ABY*I_ERI_Py_S_D2z_S;
  Double I_ERI_Pz_Py_D2z_S = I_ERI_Dyz_S_D2z_S+ABY*I_ERI_Pz_S_D2z_S;
  Double I_ERI_Px_Pz_D2z_S = I_ERI_Dxz_S_D2z_S+ABZ*I_ERI_Px_S_D2z_S;
  Double I_ERI_Py_Pz_D2z_S = I_ERI_Dyz_S_D2z_S+ABZ*I_ERI_Py_S_D2z_S;
  Double I_ERI_Pz_Pz_D2z_S = I_ERI_D2z_S_D2z_S+ABZ*I_ERI_Pz_S_D2z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S
   * RHS shell quartet name: SQ_ERI_D_S_P_S
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S = I_ERI_F3x_S_Px_S+ABX*I_ERI_D2x_S_Px_S;
  Double I_ERI_Dxy_Px_Px_S = I_ERI_F2xy_S_Px_S+ABX*I_ERI_Dxy_S_Px_S;
  Double I_ERI_Dxz_Px_Px_S = I_ERI_F2xz_S_Px_S+ABX*I_ERI_Dxz_S_Px_S;
  Double I_ERI_D2y_Px_Px_S = I_ERI_Fx2y_S_Px_S+ABX*I_ERI_D2y_S_Px_S;
  Double I_ERI_Dyz_Px_Px_S = I_ERI_Fxyz_S_Px_S+ABX*I_ERI_Dyz_S_Px_S;
  Double I_ERI_D2z_Px_Px_S = I_ERI_Fx2z_S_Px_S+ABX*I_ERI_D2z_S_Px_S;
  Double I_ERI_D2x_Py_Px_S = I_ERI_F2xy_S_Px_S+ABY*I_ERI_D2x_S_Px_S;
  Double I_ERI_Dxy_Py_Px_S = I_ERI_Fx2y_S_Px_S+ABY*I_ERI_Dxy_S_Px_S;
  Double I_ERI_Dxz_Py_Px_S = I_ERI_Fxyz_S_Px_S+ABY*I_ERI_Dxz_S_Px_S;
  Double I_ERI_D2y_Py_Px_S = I_ERI_F3y_S_Px_S+ABY*I_ERI_D2y_S_Px_S;
  Double I_ERI_Dyz_Py_Px_S = I_ERI_F2yz_S_Px_S+ABY*I_ERI_Dyz_S_Px_S;
  Double I_ERI_D2z_Py_Px_S = I_ERI_Fy2z_S_Px_S+ABY*I_ERI_D2z_S_Px_S;
  Double I_ERI_D2x_Pz_Px_S = I_ERI_F2xz_S_Px_S+ABZ*I_ERI_D2x_S_Px_S;
  Double I_ERI_Dxy_Pz_Px_S = I_ERI_Fxyz_S_Px_S+ABZ*I_ERI_Dxy_S_Px_S;
  Double I_ERI_Dxz_Pz_Px_S = I_ERI_Fx2z_S_Px_S+ABZ*I_ERI_Dxz_S_Px_S;
  Double I_ERI_D2y_Pz_Px_S = I_ERI_F2yz_S_Px_S+ABZ*I_ERI_D2y_S_Px_S;
  Double I_ERI_Dyz_Pz_Px_S = I_ERI_Fy2z_S_Px_S+ABZ*I_ERI_Dyz_S_Px_S;
  Double I_ERI_D2z_Pz_Px_S = I_ERI_F3z_S_Px_S+ABZ*I_ERI_D2z_S_Px_S;
  Double I_ERI_D2x_Px_Py_S = I_ERI_F3x_S_Py_S+ABX*I_ERI_D2x_S_Py_S;
  Double I_ERI_Dxy_Px_Py_S = I_ERI_F2xy_S_Py_S+ABX*I_ERI_Dxy_S_Py_S;
  Double I_ERI_Dxz_Px_Py_S = I_ERI_F2xz_S_Py_S+ABX*I_ERI_Dxz_S_Py_S;
  Double I_ERI_D2y_Px_Py_S = I_ERI_Fx2y_S_Py_S+ABX*I_ERI_D2y_S_Py_S;
  Double I_ERI_Dyz_Px_Py_S = I_ERI_Fxyz_S_Py_S+ABX*I_ERI_Dyz_S_Py_S;
  Double I_ERI_D2z_Px_Py_S = I_ERI_Fx2z_S_Py_S+ABX*I_ERI_D2z_S_Py_S;
  Double I_ERI_D2x_Py_Py_S = I_ERI_F2xy_S_Py_S+ABY*I_ERI_D2x_S_Py_S;
  Double I_ERI_Dxy_Py_Py_S = I_ERI_Fx2y_S_Py_S+ABY*I_ERI_Dxy_S_Py_S;
  Double I_ERI_Dxz_Py_Py_S = I_ERI_Fxyz_S_Py_S+ABY*I_ERI_Dxz_S_Py_S;
  Double I_ERI_D2y_Py_Py_S = I_ERI_F3y_S_Py_S+ABY*I_ERI_D2y_S_Py_S;
  Double I_ERI_Dyz_Py_Py_S = I_ERI_F2yz_S_Py_S+ABY*I_ERI_Dyz_S_Py_S;
  Double I_ERI_D2z_Py_Py_S = I_ERI_Fy2z_S_Py_S+ABY*I_ERI_D2z_S_Py_S;
  Double I_ERI_D2x_Pz_Py_S = I_ERI_F2xz_S_Py_S+ABZ*I_ERI_D2x_S_Py_S;
  Double I_ERI_Dxy_Pz_Py_S = I_ERI_Fxyz_S_Py_S+ABZ*I_ERI_Dxy_S_Py_S;
  Double I_ERI_Dxz_Pz_Py_S = I_ERI_Fx2z_S_Py_S+ABZ*I_ERI_Dxz_S_Py_S;
  Double I_ERI_D2y_Pz_Py_S = I_ERI_F2yz_S_Py_S+ABZ*I_ERI_D2y_S_Py_S;
  Double I_ERI_Dyz_Pz_Py_S = I_ERI_Fy2z_S_Py_S+ABZ*I_ERI_Dyz_S_Py_S;
  Double I_ERI_D2z_Pz_Py_S = I_ERI_F3z_S_Py_S+ABZ*I_ERI_D2z_S_Py_S;
  Double I_ERI_D2x_Px_Pz_S = I_ERI_F3x_S_Pz_S+ABX*I_ERI_D2x_S_Pz_S;
  Double I_ERI_Dxy_Px_Pz_S = I_ERI_F2xy_S_Pz_S+ABX*I_ERI_Dxy_S_Pz_S;
  Double I_ERI_Dxz_Px_Pz_S = I_ERI_F2xz_S_Pz_S+ABX*I_ERI_Dxz_S_Pz_S;
  Double I_ERI_D2y_Px_Pz_S = I_ERI_Fx2y_S_Pz_S+ABX*I_ERI_D2y_S_Pz_S;
  Double I_ERI_Dyz_Px_Pz_S = I_ERI_Fxyz_S_Pz_S+ABX*I_ERI_Dyz_S_Pz_S;
  Double I_ERI_D2z_Px_Pz_S = I_ERI_Fx2z_S_Pz_S+ABX*I_ERI_D2z_S_Pz_S;
  Double I_ERI_D2x_Py_Pz_S = I_ERI_F2xy_S_Pz_S+ABY*I_ERI_D2x_S_Pz_S;
  Double I_ERI_Dxy_Py_Pz_S = I_ERI_Fx2y_S_Pz_S+ABY*I_ERI_Dxy_S_Pz_S;
  Double I_ERI_Dxz_Py_Pz_S = I_ERI_Fxyz_S_Pz_S+ABY*I_ERI_Dxz_S_Pz_S;
  Double I_ERI_D2y_Py_Pz_S = I_ERI_F3y_S_Pz_S+ABY*I_ERI_D2y_S_Pz_S;
  Double I_ERI_Dyz_Py_Pz_S = I_ERI_F2yz_S_Pz_S+ABY*I_ERI_Dyz_S_Pz_S;
  Double I_ERI_D2z_Py_Pz_S = I_ERI_Fy2z_S_Pz_S+ABY*I_ERI_D2z_S_Pz_S;
  Double I_ERI_D2x_Pz_Pz_S = I_ERI_F2xz_S_Pz_S+ABZ*I_ERI_D2x_S_Pz_S;
  Double I_ERI_Dxy_Pz_Pz_S = I_ERI_Fxyz_S_Pz_S+ABZ*I_ERI_Dxy_S_Pz_S;
  Double I_ERI_Dxz_Pz_Pz_S = I_ERI_Fx2z_S_Pz_S+ABZ*I_ERI_Dxz_S_Pz_S;
  Double I_ERI_D2y_Pz_Pz_S = I_ERI_F2yz_S_Pz_S+ABZ*I_ERI_D2y_S_Pz_S;
  Double I_ERI_Dyz_Pz_Pz_S = I_ERI_Fy2z_S_Pz_S+ABZ*I_ERI_Dyz_S_Pz_S;
  Double I_ERI_D2z_Pz_Pz_S = I_ERI_F3z_S_Pz_S+ABZ*I_ERI_D2z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S
   * RHS shell quartet name: SQ_ERI_D_S_D_S
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S = I_ERI_F3x_S_D2x_S+ABX*I_ERI_D2x_S_D2x_S;
  Double I_ERI_Dxy_Px_D2x_S = I_ERI_F2xy_S_D2x_S+ABX*I_ERI_Dxy_S_D2x_S;
  Double I_ERI_Dxz_Px_D2x_S = I_ERI_F2xz_S_D2x_S+ABX*I_ERI_Dxz_S_D2x_S;
  Double I_ERI_D2y_Px_D2x_S = I_ERI_Fx2y_S_D2x_S+ABX*I_ERI_D2y_S_D2x_S;
  Double I_ERI_Dyz_Px_D2x_S = I_ERI_Fxyz_S_D2x_S+ABX*I_ERI_Dyz_S_D2x_S;
  Double I_ERI_D2z_Px_D2x_S = I_ERI_Fx2z_S_D2x_S+ABX*I_ERI_D2z_S_D2x_S;
  Double I_ERI_D2x_Py_D2x_S = I_ERI_F2xy_S_D2x_S+ABY*I_ERI_D2x_S_D2x_S;
  Double I_ERI_Dxy_Py_D2x_S = I_ERI_Fx2y_S_D2x_S+ABY*I_ERI_Dxy_S_D2x_S;
  Double I_ERI_Dxz_Py_D2x_S = I_ERI_Fxyz_S_D2x_S+ABY*I_ERI_Dxz_S_D2x_S;
  Double I_ERI_D2y_Py_D2x_S = I_ERI_F3y_S_D2x_S+ABY*I_ERI_D2y_S_D2x_S;
  Double I_ERI_Dyz_Py_D2x_S = I_ERI_F2yz_S_D2x_S+ABY*I_ERI_Dyz_S_D2x_S;
  Double I_ERI_D2z_Py_D2x_S = I_ERI_Fy2z_S_D2x_S+ABY*I_ERI_D2z_S_D2x_S;
  Double I_ERI_D2x_Pz_D2x_S = I_ERI_F2xz_S_D2x_S+ABZ*I_ERI_D2x_S_D2x_S;
  Double I_ERI_Dxy_Pz_D2x_S = I_ERI_Fxyz_S_D2x_S+ABZ*I_ERI_Dxy_S_D2x_S;
  Double I_ERI_Dxz_Pz_D2x_S = I_ERI_Fx2z_S_D2x_S+ABZ*I_ERI_Dxz_S_D2x_S;
  Double I_ERI_D2y_Pz_D2x_S = I_ERI_F2yz_S_D2x_S+ABZ*I_ERI_D2y_S_D2x_S;
  Double I_ERI_Dyz_Pz_D2x_S = I_ERI_Fy2z_S_D2x_S+ABZ*I_ERI_Dyz_S_D2x_S;
  Double I_ERI_D2z_Pz_D2x_S = I_ERI_F3z_S_D2x_S+ABZ*I_ERI_D2z_S_D2x_S;
  Double I_ERI_D2x_Px_Dxy_S = I_ERI_F3x_S_Dxy_S+ABX*I_ERI_D2x_S_Dxy_S;
  Double I_ERI_Dxy_Px_Dxy_S = I_ERI_F2xy_S_Dxy_S+ABX*I_ERI_Dxy_S_Dxy_S;
  Double I_ERI_Dxz_Px_Dxy_S = I_ERI_F2xz_S_Dxy_S+ABX*I_ERI_Dxz_S_Dxy_S;
  Double I_ERI_D2y_Px_Dxy_S = I_ERI_Fx2y_S_Dxy_S+ABX*I_ERI_D2y_S_Dxy_S;
  Double I_ERI_Dyz_Px_Dxy_S = I_ERI_Fxyz_S_Dxy_S+ABX*I_ERI_Dyz_S_Dxy_S;
  Double I_ERI_D2z_Px_Dxy_S = I_ERI_Fx2z_S_Dxy_S+ABX*I_ERI_D2z_S_Dxy_S;
  Double I_ERI_D2x_Py_Dxy_S = I_ERI_F2xy_S_Dxy_S+ABY*I_ERI_D2x_S_Dxy_S;
  Double I_ERI_Dxy_Py_Dxy_S = I_ERI_Fx2y_S_Dxy_S+ABY*I_ERI_Dxy_S_Dxy_S;
  Double I_ERI_Dxz_Py_Dxy_S = I_ERI_Fxyz_S_Dxy_S+ABY*I_ERI_Dxz_S_Dxy_S;
  Double I_ERI_D2y_Py_Dxy_S = I_ERI_F3y_S_Dxy_S+ABY*I_ERI_D2y_S_Dxy_S;
  Double I_ERI_Dyz_Py_Dxy_S = I_ERI_F2yz_S_Dxy_S+ABY*I_ERI_Dyz_S_Dxy_S;
  Double I_ERI_D2z_Py_Dxy_S = I_ERI_Fy2z_S_Dxy_S+ABY*I_ERI_D2z_S_Dxy_S;
  Double I_ERI_D2x_Pz_Dxy_S = I_ERI_F2xz_S_Dxy_S+ABZ*I_ERI_D2x_S_Dxy_S;
  Double I_ERI_Dxy_Pz_Dxy_S = I_ERI_Fxyz_S_Dxy_S+ABZ*I_ERI_Dxy_S_Dxy_S;
  Double I_ERI_Dxz_Pz_Dxy_S = I_ERI_Fx2z_S_Dxy_S+ABZ*I_ERI_Dxz_S_Dxy_S;
  Double I_ERI_D2y_Pz_Dxy_S = I_ERI_F2yz_S_Dxy_S+ABZ*I_ERI_D2y_S_Dxy_S;
  Double I_ERI_Dyz_Pz_Dxy_S = I_ERI_Fy2z_S_Dxy_S+ABZ*I_ERI_Dyz_S_Dxy_S;
  Double I_ERI_D2z_Pz_Dxy_S = I_ERI_F3z_S_Dxy_S+ABZ*I_ERI_D2z_S_Dxy_S;
  Double I_ERI_D2x_Px_Dxz_S = I_ERI_F3x_S_Dxz_S+ABX*I_ERI_D2x_S_Dxz_S;
  Double I_ERI_Dxy_Px_Dxz_S = I_ERI_F2xy_S_Dxz_S+ABX*I_ERI_Dxy_S_Dxz_S;
  Double I_ERI_Dxz_Px_Dxz_S = I_ERI_F2xz_S_Dxz_S+ABX*I_ERI_Dxz_S_Dxz_S;
  Double I_ERI_D2y_Px_Dxz_S = I_ERI_Fx2y_S_Dxz_S+ABX*I_ERI_D2y_S_Dxz_S;
  Double I_ERI_Dyz_Px_Dxz_S = I_ERI_Fxyz_S_Dxz_S+ABX*I_ERI_Dyz_S_Dxz_S;
  Double I_ERI_D2z_Px_Dxz_S = I_ERI_Fx2z_S_Dxz_S+ABX*I_ERI_D2z_S_Dxz_S;
  Double I_ERI_D2x_Py_Dxz_S = I_ERI_F2xy_S_Dxz_S+ABY*I_ERI_D2x_S_Dxz_S;
  Double I_ERI_Dxy_Py_Dxz_S = I_ERI_Fx2y_S_Dxz_S+ABY*I_ERI_Dxy_S_Dxz_S;
  Double I_ERI_Dxz_Py_Dxz_S = I_ERI_Fxyz_S_Dxz_S+ABY*I_ERI_Dxz_S_Dxz_S;
  Double I_ERI_D2y_Py_Dxz_S = I_ERI_F3y_S_Dxz_S+ABY*I_ERI_D2y_S_Dxz_S;
  Double I_ERI_Dyz_Py_Dxz_S = I_ERI_F2yz_S_Dxz_S+ABY*I_ERI_Dyz_S_Dxz_S;
  Double I_ERI_D2z_Py_Dxz_S = I_ERI_Fy2z_S_Dxz_S+ABY*I_ERI_D2z_S_Dxz_S;
  Double I_ERI_D2x_Pz_Dxz_S = I_ERI_F2xz_S_Dxz_S+ABZ*I_ERI_D2x_S_Dxz_S;
  Double I_ERI_Dxy_Pz_Dxz_S = I_ERI_Fxyz_S_Dxz_S+ABZ*I_ERI_Dxy_S_Dxz_S;
  Double I_ERI_Dxz_Pz_Dxz_S = I_ERI_Fx2z_S_Dxz_S+ABZ*I_ERI_Dxz_S_Dxz_S;
  Double I_ERI_D2y_Pz_Dxz_S = I_ERI_F2yz_S_Dxz_S+ABZ*I_ERI_D2y_S_Dxz_S;
  Double I_ERI_Dyz_Pz_Dxz_S = I_ERI_Fy2z_S_Dxz_S+ABZ*I_ERI_Dyz_S_Dxz_S;
  Double I_ERI_D2z_Pz_Dxz_S = I_ERI_F3z_S_Dxz_S+ABZ*I_ERI_D2z_S_Dxz_S;
  Double I_ERI_D2x_Px_D2y_S = I_ERI_F3x_S_D2y_S+ABX*I_ERI_D2x_S_D2y_S;
  Double I_ERI_Dxy_Px_D2y_S = I_ERI_F2xy_S_D2y_S+ABX*I_ERI_Dxy_S_D2y_S;
  Double I_ERI_Dxz_Px_D2y_S = I_ERI_F2xz_S_D2y_S+ABX*I_ERI_Dxz_S_D2y_S;
  Double I_ERI_D2y_Px_D2y_S = I_ERI_Fx2y_S_D2y_S+ABX*I_ERI_D2y_S_D2y_S;
  Double I_ERI_Dyz_Px_D2y_S = I_ERI_Fxyz_S_D2y_S+ABX*I_ERI_Dyz_S_D2y_S;
  Double I_ERI_D2z_Px_D2y_S = I_ERI_Fx2z_S_D2y_S+ABX*I_ERI_D2z_S_D2y_S;
  Double I_ERI_D2x_Py_D2y_S = I_ERI_F2xy_S_D2y_S+ABY*I_ERI_D2x_S_D2y_S;
  Double I_ERI_Dxy_Py_D2y_S = I_ERI_Fx2y_S_D2y_S+ABY*I_ERI_Dxy_S_D2y_S;
  Double I_ERI_Dxz_Py_D2y_S = I_ERI_Fxyz_S_D2y_S+ABY*I_ERI_Dxz_S_D2y_S;
  Double I_ERI_D2y_Py_D2y_S = I_ERI_F3y_S_D2y_S+ABY*I_ERI_D2y_S_D2y_S;
  Double I_ERI_Dyz_Py_D2y_S = I_ERI_F2yz_S_D2y_S+ABY*I_ERI_Dyz_S_D2y_S;
  Double I_ERI_D2z_Py_D2y_S = I_ERI_Fy2z_S_D2y_S+ABY*I_ERI_D2z_S_D2y_S;
  Double I_ERI_D2x_Pz_D2y_S = I_ERI_F2xz_S_D2y_S+ABZ*I_ERI_D2x_S_D2y_S;
  Double I_ERI_Dxy_Pz_D2y_S = I_ERI_Fxyz_S_D2y_S+ABZ*I_ERI_Dxy_S_D2y_S;
  Double I_ERI_Dxz_Pz_D2y_S = I_ERI_Fx2z_S_D2y_S+ABZ*I_ERI_Dxz_S_D2y_S;
  Double I_ERI_D2y_Pz_D2y_S = I_ERI_F2yz_S_D2y_S+ABZ*I_ERI_D2y_S_D2y_S;
  Double I_ERI_Dyz_Pz_D2y_S = I_ERI_Fy2z_S_D2y_S+ABZ*I_ERI_Dyz_S_D2y_S;
  Double I_ERI_D2z_Pz_D2y_S = I_ERI_F3z_S_D2y_S+ABZ*I_ERI_D2z_S_D2y_S;
  Double I_ERI_D2x_Px_Dyz_S = I_ERI_F3x_S_Dyz_S+ABX*I_ERI_D2x_S_Dyz_S;
  Double I_ERI_Dxy_Px_Dyz_S = I_ERI_F2xy_S_Dyz_S+ABX*I_ERI_Dxy_S_Dyz_S;
  Double I_ERI_Dxz_Px_Dyz_S = I_ERI_F2xz_S_Dyz_S+ABX*I_ERI_Dxz_S_Dyz_S;
  Double I_ERI_D2y_Px_Dyz_S = I_ERI_Fx2y_S_Dyz_S+ABX*I_ERI_D2y_S_Dyz_S;
  Double I_ERI_Dyz_Px_Dyz_S = I_ERI_Fxyz_S_Dyz_S+ABX*I_ERI_Dyz_S_Dyz_S;
  Double I_ERI_D2z_Px_Dyz_S = I_ERI_Fx2z_S_Dyz_S+ABX*I_ERI_D2z_S_Dyz_S;
  Double I_ERI_D2x_Py_Dyz_S = I_ERI_F2xy_S_Dyz_S+ABY*I_ERI_D2x_S_Dyz_S;
  Double I_ERI_Dxy_Py_Dyz_S = I_ERI_Fx2y_S_Dyz_S+ABY*I_ERI_Dxy_S_Dyz_S;
  Double I_ERI_Dxz_Py_Dyz_S = I_ERI_Fxyz_S_Dyz_S+ABY*I_ERI_Dxz_S_Dyz_S;
  Double I_ERI_D2y_Py_Dyz_S = I_ERI_F3y_S_Dyz_S+ABY*I_ERI_D2y_S_Dyz_S;
  Double I_ERI_Dyz_Py_Dyz_S = I_ERI_F2yz_S_Dyz_S+ABY*I_ERI_Dyz_S_Dyz_S;
  Double I_ERI_D2z_Py_Dyz_S = I_ERI_Fy2z_S_Dyz_S+ABY*I_ERI_D2z_S_Dyz_S;
  Double I_ERI_D2x_Pz_Dyz_S = I_ERI_F2xz_S_Dyz_S+ABZ*I_ERI_D2x_S_Dyz_S;
  Double I_ERI_Dxy_Pz_Dyz_S = I_ERI_Fxyz_S_Dyz_S+ABZ*I_ERI_Dxy_S_Dyz_S;
  Double I_ERI_Dxz_Pz_Dyz_S = I_ERI_Fx2z_S_Dyz_S+ABZ*I_ERI_Dxz_S_Dyz_S;
  Double I_ERI_D2y_Pz_Dyz_S = I_ERI_F2yz_S_Dyz_S+ABZ*I_ERI_D2y_S_Dyz_S;
  Double I_ERI_Dyz_Pz_Dyz_S = I_ERI_Fy2z_S_Dyz_S+ABZ*I_ERI_Dyz_S_Dyz_S;
  Double I_ERI_D2z_Pz_Dyz_S = I_ERI_F3z_S_Dyz_S+ABZ*I_ERI_D2z_S_Dyz_S;
  Double I_ERI_D2x_Px_D2z_S = I_ERI_F3x_S_D2z_S+ABX*I_ERI_D2x_S_D2z_S;
  Double I_ERI_Dxy_Px_D2z_S = I_ERI_F2xy_S_D2z_S+ABX*I_ERI_Dxy_S_D2z_S;
  Double I_ERI_Dxz_Px_D2z_S = I_ERI_F2xz_S_D2z_S+ABX*I_ERI_Dxz_S_D2z_S;
  Double I_ERI_D2y_Px_D2z_S = I_ERI_Fx2y_S_D2z_S+ABX*I_ERI_D2y_S_D2z_S;
  Double I_ERI_Dyz_Px_D2z_S = I_ERI_Fxyz_S_D2z_S+ABX*I_ERI_Dyz_S_D2z_S;
  Double I_ERI_D2z_Px_D2z_S = I_ERI_Fx2z_S_D2z_S+ABX*I_ERI_D2z_S_D2z_S;
  Double I_ERI_D2x_Py_D2z_S = I_ERI_F2xy_S_D2z_S+ABY*I_ERI_D2x_S_D2z_S;
  Double I_ERI_Dxy_Py_D2z_S = I_ERI_Fx2y_S_D2z_S+ABY*I_ERI_Dxy_S_D2z_S;
  Double I_ERI_Dxz_Py_D2z_S = I_ERI_Fxyz_S_D2z_S+ABY*I_ERI_Dxz_S_D2z_S;
  Double I_ERI_D2y_Py_D2z_S = I_ERI_F3y_S_D2z_S+ABY*I_ERI_D2y_S_D2z_S;
  Double I_ERI_Dyz_Py_D2z_S = I_ERI_F2yz_S_D2z_S+ABY*I_ERI_Dyz_S_D2z_S;
  Double I_ERI_D2z_Py_D2z_S = I_ERI_Fy2z_S_D2z_S+ABY*I_ERI_D2z_S_D2z_S;
  Double I_ERI_D2x_Pz_D2z_S = I_ERI_F2xz_S_D2z_S+ABZ*I_ERI_D2x_S_D2z_S;
  Double I_ERI_Dxy_Pz_D2z_S = I_ERI_Fxyz_S_D2z_S+ABZ*I_ERI_Dxy_S_D2z_S;
  Double I_ERI_Dxz_Pz_D2z_S = I_ERI_Fx2z_S_D2z_S+ABZ*I_ERI_Dxz_S_D2z_S;
  Double I_ERI_D2y_Pz_D2z_S = I_ERI_F2yz_S_D2z_S+ABZ*I_ERI_D2y_S_D2z_S;
  Double I_ERI_Dyz_Pz_D2z_S = I_ERI_Fy2z_S_D2z_S+ABZ*I_ERI_Dyz_S_D2z_S;
  Double I_ERI_D2z_Pz_D2z_S = I_ERI_F3z_S_D2z_S+ABZ*I_ERI_D2z_S_D2z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_D_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_D_S
   * RHS shell quartet name: SQ_ERI_P_P_D_S
   ************************************************************/
  Double I_ERI_Px_D2x_D2x_S = I_ERI_D2x_Px_D2x_S+ABX*I_ERI_Px_Px_D2x_S;
  Double I_ERI_Py_D2x_D2x_S = I_ERI_Dxy_Px_D2x_S+ABX*I_ERI_Py_Px_D2x_S;
  Double I_ERI_Pz_D2x_D2x_S = I_ERI_Dxz_Px_D2x_S+ABX*I_ERI_Pz_Px_D2x_S;
  Double I_ERI_Px_Dxy_D2x_S = I_ERI_Dxy_Px_D2x_S+ABY*I_ERI_Px_Px_D2x_S;
  Double I_ERI_Py_Dxy_D2x_S = I_ERI_D2y_Px_D2x_S+ABY*I_ERI_Py_Px_D2x_S;
  Double I_ERI_Pz_Dxy_D2x_S = I_ERI_Dyz_Px_D2x_S+ABY*I_ERI_Pz_Px_D2x_S;
  Double I_ERI_Px_Dxz_D2x_S = I_ERI_Dxz_Px_D2x_S+ABZ*I_ERI_Px_Px_D2x_S;
  Double I_ERI_Py_Dxz_D2x_S = I_ERI_Dyz_Px_D2x_S+ABZ*I_ERI_Py_Px_D2x_S;
  Double I_ERI_Pz_Dxz_D2x_S = I_ERI_D2z_Px_D2x_S+ABZ*I_ERI_Pz_Px_D2x_S;
  Double I_ERI_Px_D2y_D2x_S = I_ERI_Dxy_Py_D2x_S+ABY*I_ERI_Px_Py_D2x_S;
  Double I_ERI_Py_D2y_D2x_S = I_ERI_D2y_Py_D2x_S+ABY*I_ERI_Py_Py_D2x_S;
  Double I_ERI_Pz_D2y_D2x_S = I_ERI_Dyz_Py_D2x_S+ABY*I_ERI_Pz_Py_D2x_S;
  Double I_ERI_Px_Dyz_D2x_S = I_ERI_Dxz_Py_D2x_S+ABZ*I_ERI_Px_Py_D2x_S;
  Double I_ERI_Py_Dyz_D2x_S = I_ERI_Dyz_Py_D2x_S+ABZ*I_ERI_Py_Py_D2x_S;
  Double I_ERI_Pz_Dyz_D2x_S = I_ERI_D2z_Py_D2x_S+ABZ*I_ERI_Pz_Py_D2x_S;
  Double I_ERI_Px_D2z_D2x_S = I_ERI_Dxz_Pz_D2x_S+ABZ*I_ERI_Px_Pz_D2x_S;
  Double I_ERI_Py_D2z_D2x_S = I_ERI_Dyz_Pz_D2x_S+ABZ*I_ERI_Py_Pz_D2x_S;
  Double I_ERI_Pz_D2z_D2x_S = I_ERI_D2z_Pz_D2x_S+ABZ*I_ERI_Pz_Pz_D2x_S;
  Double I_ERI_Px_D2x_Dxy_S = I_ERI_D2x_Px_Dxy_S+ABX*I_ERI_Px_Px_Dxy_S;
  Double I_ERI_Py_D2x_Dxy_S = I_ERI_Dxy_Px_Dxy_S+ABX*I_ERI_Py_Px_Dxy_S;
  Double I_ERI_Pz_D2x_Dxy_S = I_ERI_Dxz_Px_Dxy_S+ABX*I_ERI_Pz_Px_Dxy_S;
  Double I_ERI_Px_Dxy_Dxy_S = I_ERI_Dxy_Px_Dxy_S+ABY*I_ERI_Px_Px_Dxy_S;
  Double I_ERI_Py_Dxy_Dxy_S = I_ERI_D2y_Px_Dxy_S+ABY*I_ERI_Py_Px_Dxy_S;
  Double I_ERI_Pz_Dxy_Dxy_S = I_ERI_Dyz_Px_Dxy_S+ABY*I_ERI_Pz_Px_Dxy_S;
  Double I_ERI_Px_Dxz_Dxy_S = I_ERI_Dxz_Px_Dxy_S+ABZ*I_ERI_Px_Px_Dxy_S;
  Double I_ERI_Py_Dxz_Dxy_S = I_ERI_Dyz_Px_Dxy_S+ABZ*I_ERI_Py_Px_Dxy_S;
  Double I_ERI_Pz_Dxz_Dxy_S = I_ERI_D2z_Px_Dxy_S+ABZ*I_ERI_Pz_Px_Dxy_S;
  Double I_ERI_Px_D2y_Dxy_S = I_ERI_Dxy_Py_Dxy_S+ABY*I_ERI_Px_Py_Dxy_S;
  Double I_ERI_Py_D2y_Dxy_S = I_ERI_D2y_Py_Dxy_S+ABY*I_ERI_Py_Py_Dxy_S;
  Double I_ERI_Pz_D2y_Dxy_S = I_ERI_Dyz_Py_Dxy_S+ABY*I_ERI_Pz_Py_Dxy_S;
  Double I_ERI_Px_Dyz_Dxy_S = I_ERI_Dxz_Py_Dxy_S+ABZ*I_ERI_Px_Py_Dxy_S;
  Double I_ERI_Py_Dyz_Dxy_S = I_ERI_Dyz_Py_Dxy_S+ABZ*I_ERI_Py_Py_Dxy_S;
  Double I_ERI_Pz_Dyz_Dxy_S = I_ERI_D2z_Py_Dxy_S+ABZ*I_ERI_Pz_Py_Dxy_S;
  Double I_ERI_Px_D2z_Dxy_S = I_ERI_Dxz_Pz_Dxy_S+ABZ*I_ERI_Px_Pz_Dxy_S;
  Double I_ERI_Py_D2z_Dxy_S = I_ERI_Dyz_Pz_Dxy_S+ABZ*I_ERI_Py_Pz_Dxy_S;
  Double I_ERI_Pz_D2z_Dxy_S = I_ERI_D2z_Pz_Dxy_S+ABZ*I_ERI_Pz_Pz_Dxy_S;
  Double I_ERI_Px_D2x_Dxz_S = I_ERI_D2x_Px_Dxz_S+ABX*I_ERI_Px_Px_Dxz_S;
  Double I_ERI_Py_D2x_Dxz_S = I_ERI_Dxy_Px_Dxz_S+ABX*I_ERI_Py_Px_Dxz_S;
  Double I_ERI_Pz_D2x_Dxz_S = I_ERI_Dxz_Px_Dxz_S+ABX*I_ERI_Pz_Px_Dxz_S;
  Double I_ERI_Px_Dxy_Dxz_S = I_ERI_Dxy_Px_Dxz_S+ABY*I_ERI_Px_Px_Dxz_S;
  Double I_ERI_Py_Dxy_Dxz_S = I_ERI_D2y_Px_Dxz_S+ABY*I_ERI_Py_Px_Dxz_S;
  Double I_ERI_Pz_Dxy_Dxz_S = I_ERI_Dyz_Px_Dxz_S+ABY*I_ERI_Pz_Px_Dxz_S;
  Double I_ERI_Px_Dxz_Dxz_S = I_ERI_Dxz_Px_Dxz_S+ABZ*I_ERI_Px_Px_Dxz_S;
  Double I_ERI_Py_Dxz_Dxz_S = I_ERI_Dyz_Px_Dxz_S+ABZ*I_ERI_Py_Px_Dxz_S;
  Double I_ERI_Pz_Dxz_Dxz_S = I_ERI_D2z_Px_Dxz_S+ABZ*I_ERI_Pz_Px_Dxz_S;
  Double I_ERI_Px_D2y_Dxz_S = I_ERI_Dxy_Py_Dxz_S+ABY*I_ERI_Px_Py_Dxz_S;
  Double I_ERI_Py_D2y_Dxz_S = I_ERI_D2y_Py_Dxz_S+ABY*I_ERI_Py_Py_Dxz_S;
  Double I_ERI_Pz_D2y_Dxz_S = I_ERI_Dyz_Py_Dxz_S+ABY*I_ERI_Pz_Py_Dxz_S;
  Double I_ERI_Px_Dyz_Dxz_S = I_ERI_Dxz_Py_Dxz_S+ABZ*I_ERI_Px_Py_Dxz_S;
  Double I_ERI_Py_Dyz_Dxz_S = I_ERI_Dyz_Py_Dxz_S+ABZ*I_ERI_Py_Py_Dxz_S;
  Double I_ERI_Pz_Dyz_Dxz_S = I_ERI_D2z_Py_Dxz_S+ABZ*I_ERI_Pz_Py_Dxz_S;
  Double I_ERI_Px_D2z_Dxz_S = I_ERI_Dxz_Pz_Dxz_S+ABZ*I_ERI_Px_Pz_Dxz_S;
  Double I_ERI_Py_D2z_Dxz_S = I_ERI_Dyz_Pz_Dxz_S+ABZ*I_ERI_Py_Pz_Dxz_S;
  Double I_ERI_Pz_D2z_Dxz_S = I_ERI_D2z_Pz_Dxz_S+ABZ*I_ERI_Pz_Pz_Dxz_S;
  Double I_ERI_Px_D2x_D2y_S = I_ERI_D2x_Px_D2y_S+ABX*I_ERI_Px_Px_D2y_S;
  Double I_ERI_Py_D2x_D2y_S = I_ERI_Dxy_Px_D2y_S+ABX*I_ERI_Py_Px_D2y_S;
  Double I_ERI_Pz_D2x_D2y_S = I_ERI_Dxz_Px_D2y_S+ABX*I_ERI_Pz_Px_D2y_S;
  Double I_ERI_Px_Dxy_D2y_S = I_ERI_Dxy_Px_D2y_S+ABY*I_ERI_Px_Px_D2y_S;
  Double I_ERI_Py_Dxy_D2y_S = I_ERI_D2y_Px_D2y_S+ABY*I_ERI_Py_Px_D2y_S;
  Double I_ERI_Pz_Dxy_D2y_S = I_ERI_Dyz_Px_D2y_S+ABY*I_ERI_Pz_Px_D2y_S;
  Double I_ERI_Px_Dxz_D2y_S = I_ERI_Dxz_Px_D2y_S+ABZ*I_ERI_Px_Px_D2y_S;
  Double I_ERI_Py_Dxz_D2y_S = I_ERI_Dyz_Px_D2y_S+ABZ*I_ERI_Py_Px_D2y_S;
  Double I_ERI_Pz_Dxz_D2y_S = I_ERI_D2z_Px_D2y_S+ABZ*I_ERI_Pz_Px_D2y_S;
  Double I_ERI_Px_D2y_D2y_S = I_ERI_Dxy_Py_D2y_S+ABY*I_ERI_Px_Py_D2y_S;
  Double I_ERI_Py_D2y_D2y_S = I_ERI_D2y_Py_D2y_S+ABY*I_ERI_Py_Py_D2y_S;
  Double I_ERI_Pz_D2y_D2y_S = I_ERI_Dyz_Py_D2y_S+ABY*I_ERI_Pz_Py_D2y_S;
  Double I_ERI_Px_Dyz_D2y_S = I_ERI_Dxz_Py_D2y_S+ABZ*I_ERI_Px_Py_D2y_S;
  Double I_ERI_Py_Dyz_D2y_S = I_ERI_Dyz_Py_D2y_S+ABZ*I_ERI_Py_Py_D2y_S;
  Double I_ERI_Pz_Dyz_D2y_S = I_ERI_D2z_Py_D2y_S+ABZ*I_ERI_Pz_Py_D2y_S;
  Double I_ERI_Px_D2z_D2y_S = I_ERI_Dxz_Pz_D2y_S+ABZ*I_ERI_Px_Pz_D2y_S;
  Double I_ERI_Py_D2z_D2y_S = I_ERI_Dyz_Pz_D2y_S+ABZ*I_ERI_Py_Pz_D2y_S;
  Double I_ERI_Pz_D2z_D2y_S = I_ERI_D2z_Pz_D2y_S+ABZ*I_ERI_Pz_Pz_D2y_S;
  Double I_ERI_Px_D2x_Dyz_S = I_ERI_D2x_Px_Dyz_S+ABX*I_ERI_Px_Px_Dyz_S;
  Double I_ERI_Py_D2x_Dyz_S = I_ERI_Dxy_Px_Dyz_S+ABX*I_ERI_Py_Px_Dyz_S;
  Double I_ERI_Pz_D2x_Dyz_S = I_ERI_Dxz_Px_Dyz_S+ABX*I_ERI_Pz_Px_Dyz_S;
  Double I_ERI_Px_Dxy_Dyz_S = I_ERI_Dxy_Px_Dyz_S+ABY*I_ERI_Px_Px_Dyz_S;
  Double I_ERI_Py_Dxy_Dyz_S = I_ERI_D2y_Px_Dyz_S+ABY*I_ERI_Py_Px_Dyz_S;
  Double I_ERI_Pz_Dxy_Dyz_S = I_ERI_Dyz_Px_Dyz_S+ABY*I_ERI_Pz_Px_Dyz_S;
  Double I_ERI_Px_Dxz_Dyz_S = I_ERI_Dxz_Px_Dyz_S+ABZ*I_ERI_Px_Px_Dyz_S;
  Double I_ERI_Py_Dxz_Dyz_S = I_ERI_Dyz_Px_Dyz_S+ABZ*I_ERI_Py_Px_Dyz_S;
  Double I_ERI_Pz_Dxz_Dyz_S = I_ERI_D2z_Px_Dyz_S+ABZ*I_ERI_Pz_Px_Dyz_S;
  Double I_ERI_Px_D2y_Dyz_S = I_ERI_Dxy_Py_Dyz_S+ABY*I_ERI_Px_Py_Dyz_S;
  Double I_ERI_Py_D2y_Dyz_S = I_ERI_D2y_Py_Dyz_S+ABY*I_ERI_Py_Py_Dyz_S;
  Double I_ERI_Pz_D2y_Dyz_S = I_ERI_Dyz_Py_Dyz_S+ABY*I_ERI_Pz_Py_Dyz_S;
  Double I_ERI_Px_Dyz_Dyz_S = I_ERI_Dxz_Py_Dyz_S+ABZ*I_ERI_Px_Py_Dyz_S;
  Double I_ERI_Py_Dyz_Dyz_S = I_ERI_Dyz_Py_Dyz_S+ABZ*I_ERI_Py_Py_Dyz_S;
  Double I_ERI_Pz_Dyz_Dyz_S = I_ERI_D2z_Py_Dyz_S+ABZ*I_ERI_Pz_Py_Dyz_S;
  Double I_ERI_Px_D2z_Dyz_S = I_ERI_Dxz_Pz_Dyz_S+ABZ*I_ERI_Px_Pz_Dyz_S;
  Double I_ERI_Py_D2z_Dyz_S = I_ERI_Dyz_Pz_Dyz_S+ABZ*I_ERI_Py_Pz_Dyz_S;
  Double I_ERI_Pz_D2z_Dyz_S = I_ERI_D2z_Pz_Dyz_S+ABZ*I_ERI_Pz_Pz_Dyz_S;
  Double I_ERI_Px_D2x_D2z_S = I_ERI_D2x_Px_D2z_S+ABX*I_ERI_Px_Px_D2z_S;
  Double I_ERI_Py_D2x_D2z_S = I_ERI_Dxy_Px_D2z_S+ABX*I_ERI_Py_Px_D2z_S;
  Double I_ERI_Pz_D2x_D2z_S = I_ERI_Dxz_Px_D2z_S+ABX*I_ERI_Pz_Px_D2z_S;
  Double I_ERI_Px_Dxy_D2z_S = I_ERI_Dxy_Px_D2z_S+ABY*I_ERI_Px_Px_D2z_S;
  Double I_ERI_Py_Dxy_D2z_S = I_ERI_D2y_Px_D2z_S+ABY*I_ERI_Py_Px_D2z_S;
  Double I_ERI_Pz_Dxy_D2z_S = I_ERI_Dyz_Px_D2z_S+ABY*I_ERI_Pz_Px_D2z_S;
  Double I_ERI_Px_Dxz_D2z_S = I_ERI_Dxz_Px_D2z_S+ABZ*I_ERI_Px_Px_D2z_S;
  Double I_ERI_Py_Dxz_D2z_S = I_ERI_Dyz_Px_D2z_S+ABZ*I_ERI_Py_Px_D2z_S;
  Double I_ERI_Pz_Dxz_D2z_S = I_ERI_D2z_Px_D2z_S+ABZ*I_ERI_Pz_Px_D2z_S;
  Double I_ERI_Px_D2y_D2z_S = I_ERI_Dxy_Py_D2z_S+ABY*I_ERI_Px_Py_D2z_S;
  Double I_ERI_Py_D2y_D2z_S = I_ERI_D2y_Py_D2z_S+ABY*I_ERI_Py_Py_D2z_S;
  Double I_ERI_Pz_D2y_D2z_S = I_ERI_Dyz_Py_D2z_S+ABY*I_ERI_Pz_Py_D2z_S;
  Double I_ERI_Px_Dyz_D2z_S = I_ERI_Dxz_Py_D2z_S+ABZ*I_ERI_Px_Py_D2z_S;
  Double I_ERI_Py_Dyz_D2z_S = I_ERI_Dyz_Py_D2z_S+ABZ*I_ERI_Py_Py_D2z_S;
  Double I_ERI_Pz_Dyz_D2z_S = I_ERI_D2z_Py_D2z_S+ABZ*I_ERI_Pz_Py_D2z_S;
  Double I_ERI_Px_D2z_D2z_S = I_ERI_Dxz_Pz_D2z_S+ABZ*I_ERI_Px_Pz_D2z_S;
  Double I_ERI_Py_D2z_D2z_S = I_ERI_Dyz_Pz_D2z_S+ABZ*I_ERI_Py_Pz_D2z_S;
  Double I_ERI_Pz_D2z_D2z_S = I_ERI_D2z_Pz_D2z_S+ABZ*I_ERI_Pz_Pz_D2z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S
   * RHS shell quartet name: SQ_ERI_F_S_P_S
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S = I_ERI_G4x_S_Px_S+ABX*I_ERI_F3x_S_Px_S;
  Double I_ERI_F2xy_Px_Px_S = I_ERI_G3xy_S_Px_S+ABX*I_ERI_F2xy_S_Px_S;
  Double I_ERI_F2xz_Px_Px_S = I_ERI_G3xz_S_Px_S+ABX*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fx2y_Px_Px_S = I_ERI_G2x2y_S_Px_S+ABX*I_ERI_Fx2y_S_Px_S;
  Double I_ERI_Fxyz_Px_Px_S = I_ERI_G2xyz_S_Px_S+ABX*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_Px_Px_S = I_ERI_G2x2z_S_Px_S+ABX*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F3y_Px_Px_S = I_ERI_Gx3y_S_Px_S+ABX*I_ERI_F3y_S_Px_S;
  Double I_ERI_F2yz_Px_Px_S = I_ERI_Gx2yz_S_Px_S+ABX*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_Px_Px_S = I_ERI_Gxy2z_S_Px_S+ABX*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_Px_Px_S = I_ERI_Gx3z_S_Px_S+ABX*I_ERI_F3z_S_Px_S;
  Double I_ERI_F2xy_Py_Px_S = I_ERI_G2x2y_S_Px_S+ABY*I_ERI_F2xy_S_Px_S;
  Double I_ERI_F2xz_Py_Px_S = I_ERI_G2xyz_S_Px_S+ABY*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fx2y_Py_Px_S = I_ERI_Gx3y_S_Px_S+ABY*I_ERI_Fx2y_S_Px_S;
  Double I_ERI_Fxyz_Py_Px_S = I_ERI_Gx2yz_S_Px_S+ABY*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_Py_Px_S = I_ERI_Gxy2z_S_Px_S+ABY*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F3y_Py_Px_S = I_ERI_G4y_S_Px_S+ABY*I_ERI_F3y_S_Px_S;
  Double I_ERI_F2yz_Py_Px_S = I_ERI_G3yz_S_Px_S+ABY*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_Py_Px_S = I_ERI_G2y2z_S_Px_S+ABY*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_Py_Px_S = I_ERI_Gy3z_S_Px_S+ABY*I_ERI_F3z_S_Px_S;
  Double I_ERI_F2xz_Pz_Px_S = I_ERI_G2x2z_S_Px_S+ABZ*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fxyz_Pz_Px_S = I_ERI_Gxy2z_S_Px_S+ABZ*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_Pz_Px_S = I_ERI_Gx3z_S_Px_S+ABZ*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F2yz_Pz_Px_S = I_ERI_G2y2z_S_Px_S+ABZ*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_Pz_Px_S = I_ERI_Gy3z_S_Px_S+ABZ*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_Pz_Px_S = I_ERI_G4z_S_Px_S+ABZ*I_ERI_F3z_S_Px_S;
  Double I_ERI_F3x_Px_Py_S = I_ERI_G4x_S_Py_S+ABX*I_ERI_F3x_S_Py_S;
  Double I_ERI_F2xy_Px_Py_S = I_ERI_G3xy_S_Py_S+ABX*I_ERI_F2xy_S_Py_S;
  Double I_ERI_F2xz_Px_Py_S = I_ERI_G3xz_S_Py_S+ABX*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fx2y_Px_Py_S = I_ERI_G2x2y_S_Py_S+ABX*I_ERI_Fx2y_S_Py_S;
  Double I_ERI_Fxyz_Px_Py_S = I_ERI_G2xyz_S_Py_S+ABX*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_Px_Py_S = I_ERI_G2x2z_S_Py_S+ABX*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F3y_Px_Py_S = I_ERI_Gx3y_S_Py_S+ABX*I_ERI_F3y_S_Py_S;
  Double I_ERI_F2yz_Px_Py_S = I_ERI_Gx2yz_S_Py_S+ABX*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_Px_Py_S = I_ERI_Gxy2z_S_Py_S+ABX*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_Px_Py_S = I_ERI_Gx3z_S_Py_S+ABX*I_ERI_F3z_S_Py_S;
  Double I_ERI_F2xy_Py_Py_S = I_ERI_G2x2y_S_Py_S+ABY*I_ERI_F2xy_S_Py_S;
  Double I_ERI_F2xz_Py_Py_S = I_ERI_G2xyz_S_Py_S+ABY*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fx2y_Py_Py_S = I_ERI_Gx3y_S_Py_S+ABY*I_ERI_Fx2y_S_Py_S;
  Double I_ERI_Fxyz_Py_Py_S = I_ERI_Gx2yz_S_Py_S+ABY*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_Py_Py_S = I_ERI_Gxy2z_S_Py_S+ABY*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F3y_Py_Py_S = I_ERI_G4y_S_Py_S+ABY*I_ERI_F3y_S_Py_S;
  Double I_ERI_F2yz_Py_Py_S = I_ERI_G3yz_S_Py_S+ABY*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_Py_Py_S = I_ERI_G2y2z_S_Py_S+ABY*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_Py_Py_S = I_ERI_Gy3z_S_Py_S+ABY*I_ERI_F3z_S_Py_S;
  Double I_ERI_F2xz_Pz_Py_S = I_ERI_G2x2z_S_Py_S+ABZ*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fxyz_Pz_Py_S = I_ERI_Gxy2z_S_Py_S+ABZ*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_Pz_Py_S = I_ERI_Gx3z_S_Py_S+ABZ*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F2yz_Pz_Py_S = I_ERI_G2y2z_S_Py_S+ABZ*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_Pz_Py_S = I_ERI_Gy3z_S_Py_S+ABZ*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_Pz_Py_S = I_ERI_G4z_S_Py_S+ABZ*I_ERI_F3z_S_Py_S;
  Double I_ERI_F3x_Px_Pz_S = I_ERI_G4x_S_Pz_S+ABX*I_ERI_F3x_S_Pz_S;
  Double I_ERI_F2xy_Px_Pz_S = I_ERI_G3xy_S_Pz_S+ABX*I_ERI_F2xy_S_Pz_S;
  Double I_ERI_F2xz_Px_Pz_S = I_ERI_G3xz_S_Pz_S+ABX*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fx2y_Px_Pz_S = I_ERI_G2x2y_S_Pz_S+ABX*I_ERI_Fx2y_S_Pz_S;
  Double I_ERI_Fxyz_Px_Pz_S = I_ERI_G2xyz_S_Pz_S+ABX*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_Px_Pz_S = I_ERI_G2x2z_S_Pz_S+ABX*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F3y_Px_Pz_S = I_ERI_Gx3y_S_Pz_S+ABX*I_ERI_F3y_S_Pz_S;
  Double I_ERI_F2yz_Px_Pz_S = I_ERI_Gx2yz_S_Pz_S+ABX*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_Px_Pz_S = I_ERI_Gxy2z_S_Pz_S+ABX*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_Px_Pz_S = I_ERI_Gx3z_S_Pz_S+ABX*I_ERI_F3z_S_Pz_S;
  Double I_ERI_F2xy_Py_Pz_S = I_ERI_G2x2y_S_Pz_S+ABY*I_ERI_F2xy_S_Pz_S;
  Double I_ERI_F2xz_Py_Pz_S = I_ERI_G2xyz_S_Pz_S+ABY*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fx2y_Py_Pz_S = I_ERI_Gx3y_S_Pz_S+ABY*I_ERI_Fx2y_S_Pz_S;
  Double I_ERI_Fxyz_Py_Pz_S = I_ERI_Gx2yz_S_Pz_S+ABY*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_Py_Pz_S = I_ERI_Gxy2z_S_Pz_S+ABY*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F3y_Py_Pz_S = I_ERI_G4y_S_Pz_S+ABY*I_ERI_F3y_S_Pz_S;
  Double I_ERI_F2yz_Py_Pz_S = I_ERI_G3yz_S_Pz_S+ABY*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_Py_Pz_S = I_ERI_G2y2z_S_Pz_S+ABY*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_Py_Pz_S = I_ERI_Gy3z_S_Pz_S+ABY*I_ERI_F3z_S_Pz_S;
  Double I_ERI_F2xz_Pz_Pz_S = I_ERI_G2x2z_S_Pz_S+ABZ*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fxyz_Pz_Pz_S = I_ERI_Gxy2z_S_Pz_S+ABZ*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_Pz_Pz_S = I_ERI_Gx3z_S_Pz_S+ABZ*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F2yz_Pz_Pz_S = I_ERI_G2y2z_S_Pz_S+ABZ*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_Pz_Pz_S = I_ERI_Gy3z_S_Pz_S+ABZ*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_Pz_Pz_S = I_ERI_G4z_S_Pz_S+ABZ*I_ERI_F3z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S
   * RHS shell quartet name: SQ_ERI_D_P_P_S
   ************************************************************/
  Double I_ERI_D2x_D2x_Px_S = I_ERI_F3x_Px_Px_S+ABX*I_ERI_D2x_Px_Px_S;
  Double I_ERI_Dxy_D2x_Px_S = I_ERI_F2xy_Px_Px_S+ABX*I_ERI_Dxy_Px_Px_S;
  Double I_ERI_Dxz_D2x_Px_S = I_ERI_F2xz_Px_Px_S+ABX*I_ERI_Dxz_Px_Px_S;
  Double I_ERI_D2y_D2x_Px_S = I_ERI_Fx2y_Px_Px_S+ABX*I_ERI_D2y_Px_Px_S;
  Double I_ERI_Dyz_D2x_Px_S = I_ERI_Fxyz_Px_Px_S+ABX*I_ERI_Dyz_Px_Px_S;
  Double I_ERI_D2z_D2x_Px_S = I_ERI_Fx2z_Px_Px_S+ABX*I_ERI_D2z_Px_Px_S;
  Double I_ERI_D2x_Dxy_Px_S = I_ERI_F2xy_Px_Px_S+ABY*I_ERI_D2x_Px_Px_S;
  Double I_ERI_Dxy_Dxy_Px_S = I_ERI_Fx2y_Px_Px_S+ABY*I_ERI_Dxy_Px_Px_S;
  Double I_ERI_Dxz_Dxy_Px_S = I_ERI_Fxyz_Px_Px_S+ABY*I_ERI_Dxz_Px_Px_S;
  Double I_ERI_D2y_Dxy_Px_S = I_ERI_F3y_Px_Px_S+ABY*I_ERI_D2y_Px_Px_S;
  Double I_ERI_Dyz_Dxy_Px_S = I_ERI_F2yz_Px_Px_S+ABY*I_ERI_Dyz_Px_Px_S;
  Double I_ERI_D2z_Dxy_Px_S = I_ERI_Fy2z_Px_Px_S+ABY*I_ERI_D2z_Px_Px_S;
  Double I_ERI_D2x_Dxz_Px_S = I_ERI_F2xz_Px_Px_S+ABZ*I_ERI_D2x_Px_Px_S;
  Double I_ERI_Dxy_Dxz_Px_S = I_ERI_Fxyz_Px_Px_S+ABZ*I_ERI_Dxy_Px_Px_S;
  Double I_ERI_Dxz_Dxz_Px_S = I_ERI_Fx2z_Px_Px_S+ABZ*I_ERI_Dxz_Px_Px_S;
  Double I_ERI_D2y_Dxz_Px_S = I_ERI_F2yz_Px_Px_S+ABZ*I_ERI_D2y_Px_Px_S;
  Double I_ERI_Dyz_Dxz_Px_S = I_ERI_Fy2z_Px_Px_S+ABZ*I_ERI_Dyz_Px_Px_S;
  Double I_ERI_D2z_Dxz_Px_S = I_ERI_F3z_Px_Px_S+ABZ*I_ERI_D2z_Px_Px_S;
  Double I_ERI_D2x_D2y_Px_S = I_ERI_F2xy_Py_Px_S+ABY*I_ERI_D2x_Py_Px_S;
  Double I_ERI_Dxy_D2y_Px_S = I_ERI_Fx2y_Py_Px_S+ABY*I_ERI_Dxy_Py_Px_S;
  Double I_ERI_Dxz_D2y_Px_S = I_ERI_Fxyz_Py_Px_S+ABY*I_ERI_Dxz_Py_Px_S;
  Double I_ERI_D2y_D2y_Px_S = I_ERI_F3y_Py_Px_S+ABY*I_ERI_D2y_Py_Px_S;
  Double I_ERI_Dyz_D2y_Px_S = I_ERI_F2yz_Py_Px_S+ABY*I_ERI_Dyz_Py_Px_S;
  Double I_ERI_D2z_D2y_Px_S = I_ERI_Fy2z_Py_Px_S+ABY*I_ERI_D2z_Py_Px_S;
  Double I_ERI_D2x_Dyz_Px_S = I_ERI_F2xz_Py_Px_S+ABZ*I_ERI_D2x_Py_Px_S;
  Double I_ERI_Dxy_Dyz_Px_S = I_ERI_Fxyz_Py_Px_S+ABZ*I_ERI_Dxy_Py_Px_S;
  Double I_ERI_Dxz_Dyz_Px_S = I_ERI_Fx2z_Py_Px_S+ABZ*I_ERI_Dxz_Py_Px_S;
  Double I_ERI_D2y_Dyz_Px_S = I_ERI_F2yz_Py_Px_S+ABZ*I_ERI_D2y_Py_Px_S;
  Double I_ERI_Dyz_Dyz_Px_S = I_ERI_Fy2z_Py_Px_S+ABZ*I_ERI_Dyz_Py_Px_S;
  Double I_ERI_D2z_Dyz_Px_S = I_ERI_F3z_Py_Px_S+ABZ*I_ERI_D2z_Py_Px_S;
  Double I_ERI_D2x_D2z_Px_S = I_ERI_F2xz_Pz_Px_S+ABZ*I_ERI_D2x_Pz_Px_S;
  Double I_ERI_Dxy_D2z_Px_S = I_ERI_Fxyz_Pz_Px_S+ABZ*I_ERI_Dxy_Pz_Px_S;
  Double I_ERI_Dxz_D2z_Px_S = I_ERI_Fx2z_Pz_Px_S+ABZ*I_ERI_Dxz_Pz_Px_S;
  Double I_ERI_D2y_D2z_Px_S = I_ERI_F2yz_Pz_Px_S+ABZ*I_ERI_D2y_Pz_Px_S;
  Double I_ERI_Dyz_D2z_Px_S = I_ERI_Fy2z_Pz_Px_S+ABZ*I_ERI_Dyz_Pz_Px_S;
  Double I_ERI_D2z_D2z_Px_S = I_ERI_F3z_Pz_Px_S+ABZ*I_ERI_D2z_Pz_Px_S;
  Double I_ERI_D2x_D2x_Py_S = I_ERI_F3x_Px_Py_S+ABX*I_ERI_D2x_Px_Py_S;
  Double I_ERI_Dxy_D2x_Py_S = I_ERI_F2xy_Px_Py_S+ABX*I_ERI_Dxy_Px_Py_S;
  Double I_ERI_Dxz_D2x_Py_S = I_ERI_F2xz_Px_Py_S+ABX*I_ERI_Dxz_Px_Py_S;
  Double I_ERI_D2y_D2x_Py_S = I_ERI_Fx2y_Px_Py_S+ABX*I_ERI_D2y_Px_Py_S;
  Double I_ERI_Dyz_D2x_Py_S = I_ERI_Fxyz_Px_Py_S+ABX*I_ERI_Dyz_Px_Py_S;
  Double I_ERI_D2z_D2x_Py_S = I_ERI_Fx2z_Px_Py_S+ABX*I_ERI_D2z_Px_Py_S;
  Double I_ERI_D2x_Dxy_Py_S = I_ERI_F2xy_Px_Py_S+ABY*I_ERI_D2x_Px_Py_S;
  Double I_ERI_Dxy_Dxy_Py_S = I_ERI_Fx2y_Px_Py_S+ABY*I_ERI_Dxy_Px_Py_S;
  Double I_ERI_Dxz_Dxy_Py_S = I_ERI_Fxyz_Px_Py_S+ABY*I_ERI_Dxz_Px_Py_S;
  Double I_ERI_D2y_Dxy_Py_S = I_ERI_F3y_Px_Py_S+ABY*I_ERI_D2y_Px_Py_S;
  Double I_ERI_Dyz_Dxy_Py_S = I_ERI_F2yz_Px_Py_S+ABY*I_ERI_Dyz_Px_Py_S;
  Double I_ERI_D2z_Dxy_Py_S = I_ERI_Fy2z_Px_Py_S+ABY*I_ERI_D2z_Px_Py_S;
  Double I_ERI_D2x_Dxz_Py_S = I_ERI_F2xz_Px_Py_S+ABZ*I_ERI_D2x_Px_Py_S;
  Double I_ERI_Dxy_Dxz_Py_S = I_ERI_Fxyz_Px_Py_S+ABZ*I_ERI_Dxy_Px_Py_S;
  Double I_ERI_Dxz_Dxz_Py_S = I_ERI_Fx2z_Px_Py_S+ABZ*I_ERI_Dxz_Px_Py_S;
  Double I_ERI_D2y_Dxz_Py_S = I_ERI_F2yz_Px_Py_S+ABZ*I_ERI_D2y_Px_Py_S;
  Double I_ERI_Dyz_Dxz_Py_S = I_ERI_Fy2z_Px_Py_S+ABZ*I_ERI_Dyz_Px_Py_S;
  Double I_ERI_D2z_Dxz_Py_S = I_ERI_F3z_Px_Py_S+ABZ*I_ERI_D2z_Px_Py_S;
  Double I_ERI_D2x_D2y_Py_S = I_ERI_F2xy_Py_Py_S+ABY*I_ERI_D2x_Py_Py_S;
  Double I_ERI_Dxy_D2y_Py_S = I_ERI_Fx2y_Py_Py_S+ABY*I_ERI_Dxy_Py_Py_S;
  Double I_ERI_Dxz_D2y_Py_S = I_ERI_Fxyz_Py_Py_S+ABY*I_ERI_Dxz_Py_Py_S;
  Double I_ERI_D2y_D2y_Py_S = I_ERI_F3y_Py_Py_S+ABY*I_ERI_D2y_Py_Py_S;
  Double I_ERI_Dyz_D2y_Py_S = I_ERI_F2yz_Py_Py_S+ABY*I_ERI_Dyz_Py_Py_S;
  Double I_ERI_D2z_D2y_Py_S = I_ERI_Fy2z_Py_Py_S+ABY*I_ERI_D2z_Py_Py_S;
  Double I_ERI_D2x_Dyz_Py_S = I_ERI_F2xz_Py_Py_S+ABZ*I_ERI_D2x_Py_Py_S;
  Double I_ERI_Dxy_Dyz_Py_S = I_ERI_Fxyz_Py_Py_S+ABZ*I_ERI_Dxy_Py_Py_S;
  Double I_ERI_Dxz_Dyz_Py_S = I_ERI_Fx2z_Py_Py_S+ABZ*I_ERI_Dxz_Py_Py_S;
  Double I_ERI_D2y_Dyz_Py_S = I_ERI_F2yz_Py_Py_S+ABZ*I_ERI_D2y_Py_Py_S;
  Double I_ERI_Dyz_Dyz_Py_S = I_ERI_Fy2z_Py_Py_S+ABZ*I_ERI_Dyz_Py_Py_S;
  Double I_ERI_D2z_Dyz_Py_S = I_ERI_F3z_Py_Py_S+ABZ*I_ERI_D2z_Py_Py_S;
  Double I_ERI_D2x_D2z_Py_S = I_ERI_F2xz_Pz_Py_S+ABZ*I_ERI_D2x_Pz_Py_S;
  Double I_ERI_Dxy_D2z_Py_S = I_ERI_Fxyz_Pz_Py_S+ABZ*I_ERI_Dxy_Pz_Py_S;
  Double I_ERI_Dxz_D2z_Py_S = I_ERI_Fx2z_Pz_Py_S+ABZ*I_ERI_Dxz_Pz_Py_S;
  Double I_ERI_D2y_D2z_Py_S = I_ERI_F2yz_Pz_Py_S+ABZ*I_ERI_D2y_Pz_Py_S;
  Double I_ERI_Dyz_D2z_Py_S = I_ERI_Fy2z_Pz_Py_S+ABZ*I_ERI_Dyz_Pz_Py_S;
  Double I_ERI_D2z_D2z_Py_S = I_ERI_F3z_Pz_Py_S+ABZ*I_ERI_D2z_Pz_Py_S;
  Double I_ERI_D2x_D2x_Pz_S = I_ERI_F3x_Px_Pz_S+ABX*I_ERI_D2x_Px_Pz_S;
  Double I_ERI_Dxy_D2x_Pz_S = I_ERI_F2xy_Px_Pz_S+ABX*I_ERI_Dxy_Px_Pz_S;
  Double I_ERI_Dxz_D2x_Pz_S = I_ERI_F2xz_Px_Pz_S+ABX*I_ERI_Dxz_Px_Pz_S;
  Double I_ERI_D2y_D2x_Pz_S = I_ERI_Fx2y_Px_Pz_S+ABX*I_ERI_D2y_Px_Pz_S;
  Double I_ERI_Dyz_D2x_Pz_S = I_ERI_Fxyz_Px_Pz_S+ABX*I_ERI_Dyz_Px_Pz_S;
  Double I_ERI_D2z_D2x_Pz_S = I_ERI_Fx2z_Px_Pz_S+ABX*I_ERI_D2z_Px_Pz_S;
  Double I_ERI_D2x_Dxy_Pz_S = I_ERI_F2xy_Px_Pz_S+ABY*I_ERI_D2x_Px_Pz_S;
  Double I_ERI_Dxy_Dxy_Pz_S = I_ERI_Fx2y_Px_Pz_S+ABY*I_ERI_Dxy_Px_Pz_S;
  Double I_ERI_Dxz_Dxy_Pz_S = I_ERI_Fxyz_Px_Pz_S+ABY*I_ERI_Dxz_Px_Pz_S;
  Double I_ERI_D2y_Dxy_Pz_S = I_ERI_F3y_Px_Pz_S+ABY*I_ERI_D2y_Px_Pz_S;
  Double I_ERI_Dyz_Dxy_Pz_S = I_ERI_F2yz_Px_Pz_S+ABY*I_ERI_Dyz_Px_Pz_S;
  Double I_ERI_D2z_Dxy_Pz_S = I_ERI_Fy2z_Px_Pz_S+ABY*I_ERI_D2z_Px_Pz_S;
  Double I_ERI_D2x_Dxz_Pz_S = I_ERI_F2xz_Px_Pz_S+ABZ*I_ERI_D2x_Px_Pz_S;
  Double I_ERI_Dxy_Dxz_Pz_S = I_ERI_Fxyz_Px_Pz_S+ABZ*I_ERI_Dxy_Px_Pz_S;
  Double I_ERI_Dxz_Dxz_Pz_S = I_ERI_Fx2z_Px_Pz_S+ABZ*I_ERI_Dxz_Px_Pz_S;
  Double I_ERI_D2y_Dxz_Pz_S = I_ERI_F2yz_Px_Pz_S+ABZ*I_ERI_D2y_Px_Pz_S;
  Double I_ERI_Dyz_Dxz_Pz_S = I_ERI_Fy2z_Px_Pz_S+ABZ*I_ERI_Dyz_Px_Pz_S;
  Double I_ERI_D2z_Dxz_Pz_S = I_ERI_F3z_Px_Pz_S+ABZ*I_ERI_D2z_Px_Pz_S;
  Double I_ERI_D2x_D2y_Pz_S = I_ERI_F2xy_Py_Pz_S+ABY*I_ERI_D2x_Py_Pz_S;
  Double I_ERI_Dxy_D2y_Pz_S = I_ERI_Fx2y_Py_Pz_S+ABY*I_ERI_Dxy_Py_Pz_S;
  Double I_ERI_Dxz_D2y_Pz_S = I_ERI_Fxyz_Py_Pz_S+ABY*I_ERI_Dxz_Py_Pz_S;
  Double I_ERI_D2y_D2y_Pz_S = I_ERI_F3y_Py_Pz_S+ABY*I_ERI_D2y_Py_Pz_S;
  Double I_ERI_Dyz_D2y_Pz_S = I_ERI_F2yz_Py_Pz_S+ABY*I_ERI_Dyz_Py_Pz_S;
  Double I_ERI_D2z_D2y_Pz_S = I_ERI_Fy2z_Py_Pz_S+ABY*I_ERI_D2z_Py_Pz_S;
  Double I_ERI_D2x_Dyz_Pz_S = I_ERI_F2xz_Py_Pz_S+ABZ*I_ERI_D2x_Py_Pz_S;
  Double I_ERI_Dxy_Dyz_Pz_S = I_ERI_Fxyz_Py_Pz_S+ABZ*I_ERI_Dxy_Py_Pz_S;
  Double I_ERI_Dxz_Dyz_Pz_S = I_ERI_Fx2z_Py_Pz_S+ABZ*I_ERI_Dxz_Py_Pz_S;
  Double I_ERI_D2y_Dyz_Pz_S = I_ERI_F2yz_Py_Pz_S+ABZ*I_ERI_D2y_Py_Pz_S;
  Double I_ERI_Dyz_Dyz_Pz_S = I_ERI_Fy2z_Py_Pz_S+ABZ*I_ERI_Dyz_Py_Pz_S;
  Double I_ERI_D2z_Dyz_Pz_S = I_ERI_F3z_Py_Pz_S+ABZ*I_ERI_D2z_Py_Pz_S;
  Double I_ERI_D2x_D2z_Pz_S = I_ERI_F2xz_Pz_Pz_S+ABZ*I_ERI_D2x_Pz_Pz_S;
  Double I_ERI_Dxy_D2z_Pz_S = I_ERI_Fxyz_Pz_Pz_S+ABZ*I_ERI_Dxy_Pz_Pz_S;
  Double I_ERI_Dxz_D2z_Pz_S = I_ERI_Fx2z_Pz_Pz_S+ABZ*I_ERI_Dxz_Pz_Pz_S;
  Double I_ERI_D2y_D2z_Pz_S = I_ERI_F2yz_Pz_Pz_S+ABZ*I_ERI_D2y_Pz_Pz_S;
  Double I_ERI_Dyz_D2z_Pz_S = I_ERI_Fy2z_Pz_Pz_S+ABZ*I_ERI_Dyz_Pz_Pz_S;
  Double I_ERI_D2z_D2z_Pz_S = I_ERI_F3z_Pz_Pz_S+ABZ*I_ERI_D2z_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_a
   * RHS shell quartet name: SQ_ERI_F_S_D_S_a
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_a = I_ERI_G4x_S_D2x_S_a+ABX*I_ERI_F3x_S_D2x_S_a;
  Double I_ERI_F2xy_Px_D2x_S_a = I_ERI_G3xy_S_D2x_S_a+ABX*I_ERI_F2xy_S_D2x_S_a;
  Double I_ERI_F2xz_Px_D2x_S_a = I_ERI_G3xz_S_D2x_S_a+ABX*I_ERI_F2xz_S_D2x_S_a;
  Double I_ERI_Fx2y_Px_D2x_S_a = I_ERI_G2x2y_S_D2x_S_a+ABX*I_ERI_Fx2y_S_D2x_S_a;
  Double I_ERI_Fxyz_Px_D2x_S_a = I_ERI_G2xyz_S_D2x_S_a+ABX*I_ERI_Fxyz_S_D2x_S_a;
  Double I_ERI_Fx2z_Px_D2x_S_a = I_ERI_G2x2z_S_D2x_S_a+ABX*I_ERI_Fx2z_S_D2x_S_a;
  Double I_ERI_F3y_Px_D2x_S_a = I_ERI_Gx3y_S_D2x_S_a+ABX*I_ERI_F3y_S_D2x_S_a;
  Double I_ERI_F2yz_Px_D2x_S_a = I_ERI_Gx2yz_S_D2x_S_a+ABX*I_ERI_F2yz_S_D2x_S_a;
  Double I_ERI_Fy2z_Px_D2x_S_a = I_ERI_Gxy2z_S_D2x_S_a+ABX*I_ERI_Fy2z_S_D2x_S_a;
  Double I_ERI_F3z_Px_D2x_S_a = I_ERI_Gx3z_S_D2x_S_a+ABX*I_ERI_F3z_S_D2x_S_a;
  Double I_ERI_F3x_Py_D2x_S_a = I_ERI_G3xy_S_D2x_S_a+ABY*I_ERI_F3x_S_D2x_S_a;
  Double I_ERI_F2xy_Py_D2x_S_a = I_ERI_G2x2y_S_D2x_S_a+ABY*I_ERI_F2xy_S_D2x_S_a;
  Double I_ERI_F2xz_Py_D2x_S_a = I_ERI_G2xyz_S_D2x_S_a+ABY*I_ERI_F2xz_S_D2x_S_a;
  Double I_ERI_Fx2y_Py_D2x_S_a = I_ERI_Gx3y_S_D2x_S_a+ABY*I_ERI_Fx2y_S_D2x_S_a;
  Double I_ERI_Fxyz_Py_D2x_S_a = I_ERI_Gx2yz_S_D2x_S_a+ABY*I_ERI_Fxyz_S_D2x_S_a;
  Double I_ERI_Fx2z_Py_D2x_S_a = I_ERI_Gxy2z_S_D2x_S_a+ABY*I_ERI_Fx2z_S_D2x_S_a;
  Double I_ERI_F3y_Py_D2x_S_a = I_ERI_G4y_S_D2x_S_a+ABY*I_ERI_F3y_S_D2x_S_a;
  Double I_ERI_F2yz_Py_D2x_S_a = I_ERI_G3yz_S_D2x_S_a+ABY*I_ERI_F2yz_S_D2x_S_a;
  Double I_ERI_Fy2z_Py_D2x_S_a = I_ERI_G2y2z_S_D2x_S_a+ABY*I_ERI_Fy2z_S_D2x_S_a;
  Double I_ERI_F3z_Py_D2x_S_a = I_ERI_Gy3z_S_D2x_S_a+ABY*I_ERI_F3z_S_D2x_S_a;
  Double I_ERI_F3x_Pz_D2x_S_a = I_ERI_G3xz_S_D2x_S_a+ABZ*I_ERI_F3x_S_D2x_S_a;
  Double I_ERI_F2xy_Pz_D2x_S_a = I_ERI_G2xyz_S_D2x_S_a+ABZ*I_ERI_F2xy_S_D2x_S_a;
  Double I_ERI_F2xz_Pz_D2x_S_a = I_ERI_G2x2z_S_D2x_S_a+ABZ*I_ERI_F2xz_S_D2x_S_a;
  Double I_ERI_Fx2y_Pz_D2x_S_a = I_ERI_Gx2yz_S_D2x_S_a+ABZ*I_ERI_Fx2y_S_D2x_S_a;
  Double I_ERI_Fxyz_Pz_D2x_S_a = I_ERI_Gxy2z_S_D2x_S_a+ABZ*I_ERI_Fxyz_S_D2x_S_a;
  Double I_ERI_Fx2z_Pz_D2x_S_a = I_ERI_Gx3z_S_D2x_S_a+ABZ*I_ERI_Fx2z_S_D2x_S_a;
  Double I_ERI_F3y_Pz_D2x_S_a = I_ERI_G3yz_S_D2x_S_a+ABZ*I_ERI_F3y_S_D2x_S_a;
  Double I_ERI_F2yz_Pz_D2x_S_a = I_ERI_G2y2z_S_D2x_S_a+ABZ*I_ERI_F2yz_S_D2x_S_a;
  Double I_ERI_Fy2z_Pz_D2x_S_a = I_ERI_Gy3z_S_D2x_S_a+ABZ*I_ERI_Fy2z_S_D2x_S_a;
  Double I_ERI_F3z_Pz_D2x_S_a = I_ERI_G4z_S_D2x_S_a+ABZ*I_ERI_F3z_S_D2x_S_a;
  Double I_ERI_F3x_Px_Dxy_S_a = I_ERI_G4x_S_Dxy_S_a+ABX*I_ERI_F3x_S_Dxy_S_a;
  Double I_ERI_F2xy_Px_Dxy_S_a = I_ERI_G3xy_S_Dxy_S_a+ABX*I_ERI_F2xy_S_Dxy_S_a;
  Double I_ERI_F2xz_Px_Dxy_S_a = I_ERI_G3xz_S_Dxy_S_a+ABX*I_ERI_F2xz_S_Dxy_S_a;
  Double I_ERI_Fx2y_Px_Dxy_S_a = I_ERI_G2x2y_S_Dxy_S_a+ABX*I_ERI_Fx2y_S_Dxy_S_a;
  Double I_ERI_Fxyz_Px_Dxy_S_a = I_ERI_G2xyz_S_Dxy_S_a+ABX*I_ERI_Fxyz_S_Dxy_S_a;
  Double I_ERI_Fx2z_Px_Dxy_S_a = I_ERI_G2x2z_S_Dxy_S_a+ABX*I_ERI_Fx2z_S_Dxy_S_a;
  Double I_ERI_F3y_Px_Dxy_S_a = I_ERI_Gx3y_S_Dxy_S_a+ABX*I_ERI_F3y_S_Dxy_S_a;
  Double I_ERI_F2yz_Px_Dxy_S_a = I_ERI_Gx2yz_S_Dxy_S_a+ABX*I_ERI_F2yz_S_Dxy_S_a;
  Double I_ERI_Fy2z_Px_Dxy_S_a = I_ERI_Gxy2z_S_Dxy_S_a+ABX*I_ERI_Fy2z_S_Dxy_S_a;
  Double I_ERI_F3z_Px_Dxy_S_a = I_ERI_Gx3z_S_Dxy_S_a+ABX*I_ERI_F3z_S_Dxy_S_a;
  Double I_ERI_F3x_Py_Dxy_S_a = I_ERI_G3xy_S_Dxy_S_a+ABY*I_ERI_F3x_S_Dxy_S_a;
  Double I_ERI_F2xy_Py_Dxy_S_a = I_ERI_G2x2y_S_Dxy_S_a+ABY*I_ERI_F2xy_S_Dxy_S_a;
  Double I_ERI_F2xz_Py_Dxy_S_a = I_ERI_G2xyz_S_Dxy_S_a+ABY*I_ERI_F2xz_S_Dxy_S_a;
  Double I_ERI_Fx2y_Py_Dxy_S_a = I_ERI_Gx3y_S_Dxy_S_a+ABY*I_ERI_Fx2y_S_Dxy_S_a;
  Double I_ERI_Fxyz_Py_Dxy_S_a = I_ERI_Gx2yz_S_Dxy_S_a+ABY*I_ERI_Fxyz_S_Dxy_S_a;
  Double I_ERI_Fx2z_Py_Dxy_S_a = I_ERI_Gxy2z_S_Dxy_S_a+ABY*I_ERI_Fx2z_S_Dxy_S_a;
  Double I_ERI_F3y_Py_Dxy_S_a = I_ERI_G4y_S_Dxy_S_a+ABY*I_ERI_F3y_S_Dxy_S_a;
  Double I_ERI_F2yz_Py_Dxy_S_a = I_ERI_G3yz_S_Dxy_S_a+ABY*I_ERI_F2yz_S_Dxy_S_a;
  Double I_ERI_Fy2z_Py_Dxy_S_a = I_ERI_G2y2z_S_Dxy_S_a+ABY*I_ERI_Fy2z_S_Dxy_S_a;
  Double I_ERI_F3z_Py_Dxy_S_a = I_ERI_Gy3z_S_Dxy_S_a+ABY*I_ERI_F3z_S_Dxy_S_a;
  Double I_ERI_F3x_Pz_Dxy_S_a = I_ERI_G3xz_S_Dxy_S_a+ABZ*I_ERI_F3x_S_Dxy_S_a;
  Double I_ERI_F2xy_Pz_Dxy_S_a = I_ERI_G2xyz_S_Dxy_S_a+ABZ*I_ERI_F2xy_S_Dxy_S_a;
  Double I_ERI_F2xz_Pz_Dxy_S_a = I_ERI_G2x2z_S_Dxy_S_a+ABZ*I_ERI_F2xz_S_Dxy_S_a;
  Double I_ERI_Fx2y_Pz_Dxy_S_a = I_ERI_Gx2yz_S_Dxy_S_a+ABZ*I_ERI_Fx2y_S_Dxy_S_a;
  Double I_ERI_Fxyz_Pz_Dxy_S_a = I_ERI_Gxy2z_S_Dxy_S_a+ABZ*I_ERI_Fxyz_S_Dxy_S_a;
  Double I_ERI_Fx2z_Pz_Dxy_S_a = I_ERI_Gx3z_S_Dxy_S_a+ABZ*I_ERI_Fx2z_S_Dxy_S_a;
  Double I_ERI_F3y_Pz_Dxy_S_a = I_ERI_G3yz_S_Dxy_S_a+ABZ*I_ERI_F3y_S_Dxy_S_a;
  Double I_ERI_F2yz_Pz_Dxy_S_a = I_ERI_G2y2z_S_Dxy_S_a+ABZ*I_ERI_F2yz_S_Dxy_S_a;
  Double I_ERI_Fy2z_Pz_Dxy_S_a = I_ERI_Gy3z_S_Dxy_S_a+ABZ*I_ERI_Fy2z_S_Dxy_S_a;
  Double I_ERI_F3z_Pz_Dxy_S_a = I_ERI_G4z_S_Dxy_S_a+ABZ*I_ERI_F3z_S_Dxy_S_a;
  Double I_ERI_F3x_Px_Dxz_S_a = I_ERI_G4x_S_Dxz_S_a+ABX*I_ERI_F3x_S_Dxz_S_a;
  Double I_ERI_F2xy_Px_Dxz_S_a = I_ERI_G3xy_S_Dxz_S_a+ABX*I_ERI_F2xy_S_Dxz_S_a;
  Double I_ERI_F2xz_Px_Dxz_S_a = I_ERI_G3xz_S_Dxz_S_a+ABX*I_ERI_F2xz_S_Dxz_S_a;
  Double I_ERI_Fx2y_Px_Dxz_S_a = I_ERI_G2x2y_S_Dxz_S_a+ABX*I_ERI_Fx2y_S_Dxz_S_a;
  Double I_ERI_Fxyz_Px_Dxz_S_a = I_ERI_G2xyz_S_Dxz_S_a+ABX*I_ERI_Fxyz_S_Dxz_S_a;
  Double I_ERI_Fx2z_Px_Dxz_S_a = I_ERI_G2x2z_S_Dxz_S_a+ABX*I_ERI_Fx2z_S_Dxz_S_a;
  Double I_ERI_F3y_Px_Dxz_S_a = I_ERI_Gx3y_S_Dxz_S_a+ABX*I_ERI_F3y_S_Dxz_S_a;
  Double I_ERI_F2yz_Px_Dxz_S_a = I_ERI_Gx2yz_S_Dxz_S_a+ABX*I_ERI_F2yz_S_Dxz_S_a;
  Double I_ERI_Fy2z_Px_Dxz_S_a = I_ERI_Gxy2z_S_Dxz_S_a+ABX*I_ERI_Fy2z_S_Dxz_S_a;
  Double I_ERI_F3z_Px_Dxz_S_a = I_ERI_Gx3z_S_Dxz_S_a+ABX*I_ERI_F3z_S_Dxz_S_a;
  Double I_ERI_F3x_Py_Dxz_S_a = I_ERI_G3xy_S_Dxz_S_a+ABY*I_ERI_F3x_S_Dxz_S_a;
  Double I_ERI_F2xy_Py_Dxz_S_a = I_ERI_G2x2y_S_Dxz_S_a+ABY*I_ERI_F2xy_S_Dxz_S_a;
  Double I_ERI_F2xz_Py_Dxz_S_a = I_ERI_G2xyz_S_Dxz_S_a+ABY*I_ERI_F2xz_S_Dxz_S_a;
  Double I_ERI_Fx2y_Py_Dxz_S_a = I_ERI_Gx3y_S_Dxz_S_a+ABY*I_ERI_Fx2y_S_Dxz_S_a;
  Double I_ERI_Fxyz_Py_Dxz_S_a = I_ERI_Gx2yz_S_Dxz_S_a+ABY*I_ERI_Fxyz_S_Dxz_S_a;
  Double I_ERI_Fx2z_Py_Dxz_S_a = I_ERI_Gxy2z_S_Dxz_S_a+ABY*I_ERI_Fx2z_S_Dxz_S_a;
  Double I_ERI_F3y_Py_Dxz_S_a = I_ERI_G4y_S_Dxz_S_a+ABY*I_ERI_F3y_S_Dxz_S_a;
  Double I_ERI_F2yz_Py_Dxz_S_a = I_ERI_G3yz_S_Dxz_S_a+ABY*I_ERI_F2yz_S_Dxz_S_a;
  Double I_ERI_Fy2z_Py_Dxz_S_a = I_ERI_G2y2z_S_Dxz_S_a+ABY*I_ERI_Fy2z_S_Dxz_S_a;
  Double I_ERI_F3z_Py_Dxz_S_a = I_ERI_Gy3z_S_Dxz_S_a+ABY*I_ERI_F3z_S_Dxz_S_a;
  Double I_ERI_F3x_Pz_Dxz_S_a = I_ERI_G3xz_S_Dxz_S_a+ABZ*I_ERI_F3x_S_Dxz_S_a;
  Double I_ERI_F2xy_Pz_Dxz_S_a = I_ERI_G2xyz_S_Dxz_S_a+ABZ*I_ERI_F2xy_S_Dxz_S_a;
  Double I_ERI_F2xz_Pz_Dxz_S_a = I_ERI_G2x2z_S_Dxz_S_a+ABZ*I_ERI_F2xz_S_Dxz_S_a;
  Double I_ERI_Fx2y_Pz_Dxz_S_a = I_ERI_Gx2yz_S_Dxz_S_a+ABZ*I_ERI_Fx2y_S_Dxz_S_a;
  Double I_ERI_Fxyz_Pz_Dxz_S_a = I_ERI_Gxy2z_S_Dxz_S_a+ABZ*I_ERI_Fxyz_S_Dxz_S_a;
  Double I_ERI_Fx2z_Pz_Dxz_S_a = I_ERI_Gx3z_S_Dxz_S_a+ABZ*I_ERI_Fx2z_S_Dxz_S_a;
  Double I_ERI_F3y_Pz_Dxz_S_a = I_ERI_G3yz_S_Dxz_S_a+ABZ*I_ERI_F3y_S_Dxz_S_a;
  Double I_ERI_F2yz_Pz_Dxz_S_a = I_ERI_G2y2z_S_Dxz_S_a+ABZ*I_ERI_F2yz_S_Dxz_S_a;
  Double I_ERI_Fy2z_Pz_Dxz_S_a = I_ERI_Gy3z_S_Dxz_S_a+ABZ*I_ERI_Fy2z_S_Dxz_S_a;
  Double I_ERI_F3z_Pz_Dxz_S_a = I_ERI_G4z_S_Dxz_S_a+ABZ*I_ERI_F3z_S_Dxz_S_a;
  Double I_ERI_F3x_Px_D2y_S_a = I_ERI_G4x_S_D2y_S_a+ABX*I_ERI_F3x_S_D2y_S_a;
  Double I_ERI_F2xy_Px_D2y_S_a = I_ERI_G3xy_S_D2y_S_a+ABX*I_ERI_F2xy_S_D2y_S_a;
  Double I_ERI_F2xz_Px_D2y_S_a = I_ERI_G3xz_S_D2y_S_a+ABX*I_ERI_F2xz_S_D2y_S_a;
  Double I_ERI_Fx2y_Px_D2y_S_a = I_ERI_G2x2y_S_D2y_S_a+ABX*I_ERI_Fx2y_S_D2y_S_a;
  Double I_ERI_Fxyz_Px_D2y_S_a = I_ERI_G2xyz_S_D2y_S_a+ABX*I_ERI_Fxyz_S_D2y_S_a;
  Double I_ERI_Fx2z_Px_D2y_S_a = I_ERI_G2x2z_S_D2y_S_a+ABX*I_ERI_Fx2z_S_D2y_S_a;
  Double I_ERI_F3y_Px_D2y_S_a = I_ERI_Gx3y_S_D2y_S_a+ABX*I_ERI_F3y_S_D2y_S_a;
  Double I_ERI_F2yz_Px_D2y_S_a = I_ERI_Gx2yz_S_D2y_S_a+ABX*I_ERI_F2yz_S_D2y_S_a;
  Double I_ERI_Fy2z_Px_D2y_S_a = I_ERI_Gxy2z_S_D2y_S_a+ABX*I_ERI_Fy2z_S_D2y_S_a;
  Double I_ERI_F3z_Px_D2y_S_a = I_ERI_Gx3z_S_D2y_S_a+ABX*I_ERI_F3z_S_D2y_S_a;
  Double I_ERI_F3x_Py_D2y_S_a = I_ERI_G3xy_S_D2y_S_a+ABY*I_ERI_F3x_S_D2y_S_a;
  Double I_ERI_F2xy_Py_D2y_S_a = I_ERI_G2x2y_S_D2y_S_a+ABY*I_ERI_F2xy_S_D2y_S_a;
  Double I_ERI_F2xz_Py_D2y_S_a = I_ERI_G2xyz_S_D2y_S_a+ABY*I_ERI_F2xz_S_D2y_S_a;
  Double I_ERI_Fx2y_Py_D2y_S_a = I_ERI_Gx3y_S_D2y_S_a+ABY*I_ERI_Fx2y_S_D2y_S_a;
  Double I_ERI_Fxyz_Py_D2y_S_a = I_ERI_Gx2yz_S_D2y_S_a+ABY*I_ERI_Fxyz_S_D2y_S_a;
  Double I_ERI_Fx2z_Py_D2y_S_a = I_ERI_Gxy2z_S_D2y_S_a+ABY*I_ERI_Fx2z_S_D2y_S_a;
  Double I_ERI_F3y_Py_D2y_S_a = I_ERI_G4y_S_D2y_S_a+ABY*I_ERI_F3y_S_D2y_S_a;
  Double I_ERI_F2yz_Py_D2y_S_a = I_ERI_G3yz_S_D2y_S_a+ABY*I_ERI_F2yz_S_D2y_S_a;
  Double I_ERI_Fy2z_Py_D2y_S_a = I_ERI_G2y2z_S_D2y_S_a+ABY*I_ERI_Fy2z_S_D2y_S_a;
  Double I_ERI_F3z_Py_D2y_S_a = I_ERI_Gy3z_S_D2y_S_a+ABY*I_ERI_F3z_S_D2y_S_a;
  Double I_ERI_F3x_Pz_D2y_S_a = I_ERI_G3xz_S_D2y_S_a+ABZ*I_ERI_F3x_S_D2y_S_a;
  Double I_ERI_F2xy_Pz_D2y_S_a = I_ERI_G2xyz_S_D2y_S_a+ABZ*I_ERI_F2xy_S_D2y_S_a;
  Double I_ERI_F2xz_Pz_D2y_S_a = I_ERI_G2x2z_S_D2y_S_a+ABZ*I_ERI_F2xz_S_D2y_S_a;
  Double I_ERI_Fx2y_Pz_D2y_S_a = I_ERI_Gx2yz_S_D2y_S_a+ABZ*I_ERI_Fx2y_S_D2y_S_a;
  Double I_ERI_Fxyz_Pz_D2y_S_a = I_ERI_Gxy2z_S_D2y_S_a+ABZ*I_ERI_Fxyz_S_D2y_S_a;
  Double I_ERI_Fx2z_Pz_D2y_S_a = I_ERI_Gx3z_S_D2y_S_a+ABZ*I_ERI_Fx2z_S_D2y_S_a;
  Double I_ERI_F3y_Pz_D2y_S_a = I_ERI_G3yz_S_D2y_S_a+ABZ*I_ERI_F3y_S_D2y_S_a;
  Double I_ERI_F2yz_Pz_D2y_S_a = I_ERI_G2y2z_S_D2y_S_a+ABZ*I_ERI_F2yz_S_D2y_S_a;
  Double I_ERI_Fy2z_Pz_D2y_S_a = I_ERI_Gy3z_S_D2y_S_a+ABZ*I_ERI_Fy2z_S_D2y_S_a;
  Double I_ERI_F3z_Pz_D2y_S_a = I_ERI_G4z_S_D2y_S_a+ABZ*I_ERI_F3z_S_D2y_S_a;
  Double I_ERI_F3x_Px_Dyz_S_a = I_ERI_G4x_S_Dyz_S_a+ABX*I_ERI_F3x_S_Dyz_S_a;
  Double I_ERI_F2xy_Px_Dyz_S_a = I_ERI_G3xy_S_Dyz_S_a+ABX*I_ERI_F2xy_S_Dyz_S_a;
  Double I_ERI_F2xz_Px_Dyz_S_a = I_ERI_G3xz_S_Dyz_S_a+ABX*I_ERI_F2xz_S_Dyz_S_a;
  Double I_ERI_Fx2y_Px_Dyz_S_a = I_ERI_G2x2y_S_Dyz_S_a+ABX*I_ERI_Fx2y_S_Dyz_S_a;
  Double I_ERI_Fxyz_Px_Dyz_S_a = I_ERI_G2xyz_S_Dyz_S_a+ABX*I_ERI_Fxyz_S_Dyz_S_a;
  Double I_ERI_Fx2z_Px_Dyz_S_a = I_ERI_G2x2z_S_Dyz_S_a+ABX*I_ERI_Fx2z_S_Dyz_S_a;
  Double I_ERI_F3y_Px_Dyz_S_a = I_ERI_Gx3y_S_Dyz_S_a+ABX*I_ERI_F3y_S_Dyz_S_a;
  Double I_ERI_F2yz_Px_Dyz_S_a = I_ERI_Gx2yz_S_Dyz_S_a+ABX*I_ERI_F2yz_S_Dyz_S_a;
  Double I_ERI_Fy2z_Px_Dyz_S_a = I_ERI_Gxy2z_S_Dyz_S_a+ABX*I_ERI_Fy2z_S_Dyz_S_a;
  Double I_ERI_F3z_Px_Dyz_S_a = I_ERI_Gx3z_S_Dyz_S_a+ABX*I_ERI_F3z_S_Dyz_S_a;
  Double I_ERI_F3x_Py_Dyz_S_a = I_ERI_G3xy_S_Dyz_S_a+ABY*I_ERI_F3x_S_Dyz_S_a;
  Double I_ERI_F2xy_Py_Dyz_S_a = I_ERI_G2x2y_S_Dyz_S_a+ABY*I_ERI_F2xy_S_Dyz_S_a;
  Double I_ERI_F2xz_Py_Dyz_S_a = I_ERI_G2xyz_S_Dyz_S_a+ABY*I_ERI_F2xz_S_Dyz_S_a;
  Double I_ERI_Fx2y_Py_Dyz_S_a = I_ERI_Gx3y_S_Dyz_S_a+ABY*I_ERI_Fx2y_S_Dyz_S_a;
  Double I_ERI_Fxyz_Py_Dyz_S_a = I_ERI_Gx2yz_S_Dyz_S_a+ABY*I_ERI_Fxyz_S_Dyz_S_a;
  Double I_ERI_Fx2z_Py_Dyz_S_a = I_ERI_Gxy2z_S_Dyz_S_a+ABY*I_ERI_Fx2z_S_Dyz_S_a;
  Double I_ERI_F3y_Py_Dyz_S_a = I_ERI_G4y_S_Dyz_S_a+ABY*I_ERI_F3y_S_Dyz_S_a;
  Double I_ERI_F2yz_Py_Dyz_S_a = I_ERI_G3yz_S_Dyz_S_a+ABY*I_ERI_F2yz_S_Dyz_S_a;
  Double I_ERI_Fy2z_Py_Dyz_S_a = I_ERI_G2y2z_S_Dyz_S_a+ABY*I_ERI_Fy2z_S_Dyz_S_a;
  Double I_ERI_F3z_Py_Dyz_S_a = I_ERI_Gy3z_S_Dyz_S_a+ABY*I_ERI_F3z_S_Dyz_S_a;
  Double I_ERI_F3x_Pz_Dyz_S_a = I_ERI_G3xz_S_Dyz_S_a+ABZ*I_ERI_F3x_S_Dyz_S_a;
  Double I_ERI_F2xy_Pz_Dyz_S_a = I_ERI_G2xyz_S_Dyz_S_a+ABZ*I_ERI_F2xy_S_Dyz_S_a;
  Double I_ERI_F2xz_Pz_Dyz_S_a = I_ERI_G2x2z_S_Dyz_S_a+ABZ*I_ERI_F2xz_S_Dyz_S_a;
  Double I_ERI_Fx2y_Pz_Dyz_S_a = I_ERI_Gx2yz_S_Dyz_S_a+ABZ*I_ERI_Fx2y_S_Dyz_S_a;
  Double I_ERI_Fxyz_Pz_Dyz_S_a = I_ERI_Gxy2z_S_Dyz_S_a+ABZ*I_ERI_Fxyz_S_Dyz_S_a;
  Double I_ERI_Fx2z_Pz_Dyz_S_a = I_ERI_Gx3z_S_Dyz_S_a+ABZ*I_ERI_Fx2z_S_Dyz_S_a;
  Double I_ERI_F3y_Pz_Dyz_S_a = I_ERI_G3yz_S_Dyz_S_a+ABZ*I_ERI_F3y_S_Dyz_S_a;
  Double I_ERI_F2yz_Pz_Dyz_S_a = I_ERI_G2y2z_S_Dyz_S_a+ABZ*I_ERI_F2yz_S_Dyz_S_a;
  Double I_ERI_Fy2z_Pz_Dyz_S_a = I_ERI_Gy3z_S_Dyz_S_a+ABZ*I_ERI_Fy2z_S_Dyz_S_a;
  Double I_ERI_F3z_Pz_Dyz_S_a = I_ERI_G4z_S_Dyz_S_a+ABZ*I_ERI_F3z_S_Dyz_S_a;
  Double I_ERI_F3x_Px_D2z_S_a = I_ERI_G4x_S_D2z_S_a+ABX*I_ERI_F3x_S_D2z_S_a;
  Double I_ERI_F2xy_Px_D2z_S_a = I_ERI_G3xy_S_D2z_S_a+ABX*I_ERI_F2xy_S_D2z_S_a;
  Double I_ERI_F2xz_Px_D2z_S_a = I_ERI_G3xz_S_D2z_S_a+ABX*I_ERI_F2xz_S_D2z_S_a;
  Double I_ERI_Fx2y_Px_D2z_S_a = I_ERI_G2x2y_S_D2z_S_a+ABX*I_ERI_Fx2y_S_D2z_S_a;
  Double I_ERI_Fxyz_Px_D2z_S_a = I_ERI_G2xyz_S_D2z_S_a+ABX*I_ERI_Fxyz_S_D2z_S_a;
  Double I_ERI_Fx2z_Px_D2z_S_a = I_ERI_G2x2z_S_D2z_S_a+ABX*I_ERI_Fx2z_S_D2z_S_a;
  Double I_ERI_F3y_Px_D2z_S_a = I_ERI_Gx3y_S_D2z_S_a+ABX*I_ERI_F3y_S_D2z_S_a;
  Double I_ERI_F2yz_Px_D2z_S_a = I_ERI_Gx2yz_S_D2z_S_a+ABX*I_ERI_F2yz_S_D2z_S_a;
  Double I_ERI_Fy2z_Px_D2z_S_a = I_ERI_Gxy2z_S_D2z_S_a+ABX*I_ERI_Fy2z_S_D2z_S_a;
  Double I_ERI_F3z_Px_D2z_S_a = I_ERI_Gx3z_S_D2z_S_a+ABX*I_ERI_F3z_S_D2z_S_a;
  Double I_ERI_F3x_Py_D2z_S_a = I_ERI_G3xy_S_D2z_S_a+ABY*I_ERI_F3x_S_D2z_S_a;
  Double I_ERI_F2xy_Py_D2z_S_a = I_ERI_G2x2y_S_D2z_S_a+ABY*I_ERI_F2xy_S_D2z_S_a;
  Double I_ERI_F2xz_Py_D2z_S_a = I_ERI_G2xyz_S_D2z_S_a+ABY*I_ERI_F2xz_S_D2z_S_a;
  Double I_ERI_Fx2y_Py_D2z_S_a = I_ERI_Gx3y_S_D2z_S_a+ABY*I_ERI_Fx2y_S_D2z_S_a;
  Double I_ERI_Fxyz_Py_D2z_S_a = I_ERI_Gx2yz_S_D2z_S_a+ABY*I_ERI_Fxyz_S_D2z_S_a;
  Double I_ERI_Fx2z_Py_D2z_S_a = I_ERI_Gxy2z_S_D2z_S_a+ABY*I_ERI_Fx2z_S_D2z_S_a;
  Double I_ERI_F3y_Py_D2z_S_a = I_ERI_G4y_S_D2z_S_a+ABY*I_ERI_F3y_S_D2z_S_a;
  Double I_ERI_F2yz_Py_D2z_S_a = I_ERI_G3yz_S_D2z_S_a+ABY*I_ERI_F2yz_S_D2z_S_a;
  Double I_ERI_Fy2z_Py_D2z_S_a = I_ERI_G2y2z_S_D2z_S_a+ABY*I_ERI_Fy2z_S_D2z_S_a;
  Double I_ERI_F3z_Py_D2z_S_a = I_ERI_Gy3z_S_D2z_S_a+ABY*I_ERI_F3z_S_D2z_S_a;
  Double I_ERI_F3x_Pz_D2z_S_a = I_ERI_G3xz_S_D2z_S_a+ABZ*I_ERI_F3x_S_D2z_S_a;
  Double I_ERI_F2xy_Pz_D2z_S_a = I_ERI_G2xyz_S_D2z_S_a+ABZ*I_ERI_F2xy_S_D2z_S_a;
  Double I_ERI_F2xz_Pz_D2z_S_a = I_ERI_G2x2z_S_D2z_S_a+ABZ*I_ERI_F2xz_S_D2z_S_a;
  Double I_ERI_Fx2y_Pz_D2z_S_a = I_ERI_Gx2yz_S_D2z_S_a+ABZ*I_ERI_Fx2y_S_D2z_S_a;
  Double I_ERI_Fxyz_Pz_D2z_S_a = I_ERI_Gxy2z_S_D2z_S_a+ABZ*I_ERI_Fxyz_S_D2z_S_a;
  Double I_ERI_Fx2z_Pz_D2z_S_a = I_ERI_Gx3z_S_D2z_S_a+ABZ*I_ERI_Fx2z_S_D2z_S_a;
  Double I_ERI_F3y_Pz_D2z_S_a = I_ERI_G3yz_S_D2z_S_a+ABZ*I_ERI_F3y_S_D2z_S_a;
  Double I_ERI_F2yz_Pz_D2z_S_a = I_ERI_G2y2z_S_D2z_S_a+ABZ*I_ERI_F2yz_S_D2z_S_a;
  Double I_ERI_Fy2z_Pz_D2z_S_a = I_ERI_Gy3z_S_D2z_S_a+ABZ*I_ERI_Fy2z_S_D2z_S_a;
  Double I_ERI_F3z_Pz_D2z_S_a = I_ERI_G4z_S_D2z_S_a+ABZ*I_ERI_F3z_S_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_D_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_D_S_a
   * RHS shell quartet name: SQ_ERI_G_S_D_S_a
   ************************************************************/
  Double I_ERI_G4x_Px_D2x_S_a = I_ERI_H5x_S_D2x_S_a+ABX*I_ERI_G4x_S_D2x_S_a;
  Double I_ERI_G3xy_Px_D2x_S_a = I_ERI_H4xy_S_D2x_S_a+ABX*I_ERI_G3xy_S_D2x_S_a;
  Double I_ERI_G3xz_Px_D2x_S_a = I_ERI_H4xz_S_D2x_S_a+ABX*I_ERI_G3xz_S_D2x_S_a;
  Double I_ERI_G2x2y_Px_D2x_S_a = I_ERI_H3x2y_S_D2x_S_a+ABX*I_ERI_G2x2y_S_D2x_S_a;
  Double I_ERI_G2xyz_Px_D2x_S_a = I_ERI_H3xyz_S_D2x_S_a+ABX*I_ERI_G2xyz_S_D2x_S_a;
  Double I_ERI_G2x2z_Px_D2x_S_a = I_ERI_H3x2z_S_D2x_S_a+ABX*I_ERI_G2x2z_S_D2x_S_a;
  Double I_ERI_Gx3y_Px_D2x_S_a = I_ERI_H2x3y_S_D2x_S_a+ABX*I_ERI_Gx3y_S_D2x_S_a;
  Double I_ERI_Gx2yz_Px_D2x_S_a = I_ERI_H2x2yz_S_D2x_S_a+ABX*I_ERI_Gx2yz_S_D2x_S_a;
  Double I_ERI_Gxy2z_Px_D2x_S_a = I_ERI_H2xy2z_S_D2x_S_a+ABX*I_ERI_Gxy2z_S_D2x_S_a;
  Double I_ERI_Gx3z_Px_D2x_S_a = I_ERI_H2x3z_S_D2x_S_a+ABX*I_ERI_Gx3z_S_D2x_S_a;
  Double I_ERI_G4y_Px_D2x_S_a = I_ERI_Hx4y_S_D2x_S_a+ABX*I_ERI_G4y_S_D2x_S_a;
  Double I_ERI_G3yz_Px_D2x_S_a = I_ERI_Hx3yz_S_D2x_S_a+ABX*I_ERI_G3yz_S_D2x_S_a;
  Double I_ERI_G2y2z_Px_D2x_S_a = I_ERI_Hx2y2z_S_D2x_S_a+ABX*I_ERI_G2y2z_S_D2x_S_a;
  Double I_ERI_Gy3z_Px_D2x_S_a = I_ERI_Hxy3z_S_D2x_S_a+ABX*I_ERI_Gy3z_S_D2x_S_a;
  Double I_ERI_G4z_Px_D2x_S_a = I_ERI_Hx4z_S_D2x_S_a+ABX*I_ERI_G4z_S_D2x_S_a;
  Double I_ERI_G3xy_Py_D2x_S_a = I_ERI_H3x2y_S_D2x_S_a+ABY*I_ERI_G3xy_S_D2x_S_a;
  Double I_ERI_G3xz_Py_D2x_S_a = I_ERI_H3xyz_S_D2x_S_a+ABY*I_ERI_G3xz_S_D2x_S_a;
  Double I_ERI_G2x2y_Py_D2x_S_a = I_ERI_H2x3y_S_D2x_S_a+ABY*I_ERI_G2x2y_S_D2x_S_a;
  Double I_ERI_G2xyz_Py_D2x_S_a = I_ERI_H2x2yz_S_D2x_S_a+ABY*I_ERI_G2xyz_S_D2x_S_a;
  Double I_ERI_G2x2z_Py_D2x_S_a = I_ERI_H2xy2z_S_D2x_S_a+ABY*I_ERI_G2x2z_S_D2x_S_a;
  Double I_ERI_Gx3y_Py_D2x_S_a = I_ERI_Hx4y_S_D2x_S_a+ABY*I_ERI_Gx3y_S_D2x_S_a;
  Double I_ERI_Gx2yz_Py_D2x_S_a = I_ERI_Hx3yz_S_D2x_S_a+ABY*I_ERI_Gx2yz_S_D2x_S_a;
  Double I_ERI_Gxy2z_Py_D2x_S_a = I_ERI_Hx2y2z_S_D2x_S_a+ABY*I_ERI_Gxy2z_S_D2x_S_a;
  Double I_ERI_Gx3z_Py_D2x_S_a = I_ERI_Hxy3z_S_D2x_S_a+ABY*I_ERI_Gx3z_S_D2x_S_a;
  Double I_ERI_G4y_Py_D2x_S_a = I_ERI_H5y_S_D2x_S_a+ABY*I_ERI_G4y_S_D2x_S_a;
  Double I_ERI_G3yz_Py_D2x_S_a = I_ERI_H4yz_S_D2x_S_a+ABY*I_ERI_G3yz_S_D2x_S_a;
  Double I_ERI_G2y2z_Py_D2x_S_a = I_ERI_H3y2z_S_D2x_S_a+ABY*I_ERI_G2y2z_S_D2x_S_a;
  Double I_ERI_Gy3z_Py_D2x_S_a = I_ERI_H2y3z_S_D2x_S_a+ABY*I_ERI_Gy3z_S_D2x_S_a;
  Double I_ERI_G4z_Py_D2x_S_a = I_ERI_Hy4z_S_D2x_S_a+ABY*I_ERI_G4z_S_D2x_S_a;
  Double I_ERI_G3xz_Pz_D2x_S_a = I_ERI_H3x2z_S_D2x_S_a+ABZ*I_ERI_G3xz_S_D2x_S_a;
  Double I_ERI_G2xyz_Pz_D2x_S_a = I_ERI_H2xy2z_S_D2x_S_a+ABZ*I_ERI_G2xyz_S_D2x_S_a;
  Double I_ERI_G2x2z_Pz_D2x_S_a = I_ERI_H2x3z_S_D2x_S_a+ABZ*I_ERI_G2x2z_S_D2x_S_a;
  Double I_ERI_Gx2yz_Pz_D2x_S_a = I_ERI_Hx2y2z_S_D2x_S_a+ABZ*I_ERI_Gx2yz_S_D2x_S_a;
  Double I_ERI_Gxy2z_Pz_D2x_S_a = I_ERI_Hxy3z_S_D2x_S_a+ABZ*I_ERI_Gxy2z_S_D2x_S_a;
  Double I_ERI_Gx3z_Pz_D2x_S_a = I_ERI_Hx4z_S_D2x_S_a+ABZ*I_ERI_Gx3z_S_D2x_S_a;
  Double I_ERI_G3yz_Pz_D2x_S_a = I_ERI_H3y2z_S_D2x_S_a+ABZ*I_ERI_G3yz_S_D2x_S_a;
  Double I_ERI_G2y2z_Pz_D2x_S_a = I_ERI_H2y3z_S_D2x_S_a+ABZ*I_ERI_G2y2z_S_D2x_S_a;
  Double I_ERI_Gy3z_Pz_D2x_S_a = I_ERI_Hy4z_S_D2x_S_a+ABZ*I_ERI_Gy3z_S_D2x_S_a;
  Double I_ERI_G4z_Pz_D2x_S_a = I_ERI_H5z_S_D2x_S_a+ABZ*I_ERI_G4z_S_D2x_S_a;
  Double I_ERI_G4x_Px_Dxy_S_a = I_ERI_H5x_S_Dxy_S_a+ABX*I_ERI_G4x_S_Dxy_S_a;
  Double I_ERI_G3xy_Px_Dxy_S_a = I_ERI_H4xy_S_Dxy_S_a+ABX*I_ERI_G3xy_S_Dxy_S_a;
  Double I_ERI_G3xz_Px_Dxy_S_a = I_ERI_H4xz_S_Dxy_S_a+ABX*I_ERI_G3xz_S_Dxy_S_a;
  Double I_ERI_G2x2y_Px_Dxy_S_a = I_ERI_H3x2y_S_Dxy_S_a+ABX*I_ERI_G2x2y_S_Dxy_S_a;
  Double I_ERI_G2xyz_Px_Dxy_S_a = I_ERI_H3xyz_S_Dxy_S_a+ABX*I_ERI_G2xyz_S_Dxy_S_a;
  Double I_ERI_G2x2z_Px_Dxy_S_a = I_ERI_H3x2z_S_Dxy_S_a+ABX*I_ERI_G2x2z_S_Dxy_S_a;
  Double I_ERI_Gx3y_Px_Dxy_S_a = I_ERI_H2x3y_S_Dxy_S_a+ABX*I_ERI_Gx3y_S_Dxy_S_a;
  Double I_ERI_Gx2yz_Px_Dxy_S_a = I_ERI_H2x2yz_S_Dxy_S_a+ABX*I_ERI_Gx2yz_S_Dxy_S_a;
  Double I_ERI_Gxy2z_Px_Dxy_S_a = I_ERI_H2xy2z_S_Dxy_S_a+ABX*I_ERI_Gxy2z_S_Dxy_S_a;
  Double I_ERI_Gx3z_Px_Dxy_S_a = I_ERI_H2x3z_S_Dxy_S_a+ABX*I_ERI_Gx3z_S_Dxy_S_a;
  Double I_ERI_G4y_Px_Dxy_S_a = I_ERI_Hx4y_S_Dxy_S_a+ABX*I_ERI_G4y_S_Dxy_S_a;
  Double I_ERI_G3yz_Px_Dxy_S_a = I_ERI_Hx3yz_S_Dxy_S_a+ABX*I_ERI_G3yz_S_Dxy_S_a;
  Double I_ERI_G2y2z_Px_Dxy_S_a = I_ERI_Hx2y2z_S_Dxy_S_a+ABX*I_ERI_G2y2z_S_Dxy_S_a;
  Double I_ERI_Gy3z_Px_Dxy_S_a = I_ERI_Hxy3z_S_Dxy_S_a+ABX*I_ERI_Gy3z_S_Dxy_S_a;
  Double I_ERI_G4z_Px_Dxy_S_a = I_ERI_Hx4z_S_Dxy_S_a+ABX*I_ERI_G4z_S_Dxy_S_a;
  Double I_ERI_G3xy_Py_Dxy_S_a = I_ERI_H3x2y_S_Dxy_S_a+ABY*I_ERI_G3xy_S_Dxy_S_a;
  Double I_ERI_G3xz_Py_Dxy_S_a = I_ERI_H3xyz_S_Dxy_S_a+ABY*I_ERI_G3xz_S_Dxy_S_a;
  Double I_ERI_G2x2y_Py_Dxy_S_a = I_ERI_H2x3y_S_Dxy_S_a+ABY*I_ERI_G2x2y_S_Dxy_S_a;
  Double I_ERI_G2xyz_Py_Dxy_S_a = I_ERI_H2x2yz_S_Dxy_S_a+ABY*I_ERI_G2xyz_S_Dxy_S_a;
  Double I_ERI_G2x2z_Py_Dxy_S_a = I_ERI_H2xy2z_S_Dxy_S_a+ABY*I_ERI_G2x2z_S_Dxy_S_a;
  Double I_ERI_Gx3y_Py_Dxy_S_a = I_ERI_Hx4y_S_Dxy_S_a+ABY*I_ERI_Gx3y_S_Dxy_S_a;
  Double I_ERI_Gx2yz_Py_Dxy_S_a = I_ERI_Hx3yz_S_Dxy_S_a+ABY*I_ERI_Gx2yz_S_Dxy_S_a;
  Double I_ERI_Gxy2z_Py_Dxy_S_a = I_ERI_Hx2y2z_S_Dxy_S_a+ABY*I_ERI_Gxy2z_S_Dxy_S_a;
  Double I_ERI_Gx3z_Py_Dxy_S_a = I_ERI_Hxy3z_S_Dxy_S_a+ABY*I_ERI_Gx3z_S_Dxy_S_a;
  Double I_ERI_G4y_Py_Dxy_S_a = I_ERI_H5y_S_Dxy_S_a+ABY*I_ERI_G4y_S_Dxy_S_a;
  Double I_ERI_G3yz_Py_Dxy_S_a = I_ERI_H4yz_S_Dxy_S_a+ABY*I_ERI_G3yz_S_Dxy_S_a;
  Double I_ERI_G2y2z_Py_Dxy_S_a = I_ERI_H3y2z_S_Dxy_S_a+ABY*I_ERI_G2y2z_S_Dxy_S_a;
  Double I_ERI_Gy3z_Py_Dxy_S_a = I_ERI_H2y3z_S_Dxy_S_a+ABY*I_ERI_Gy3z_S_Dxy_S_a;
  Double I_ERI_G4z_Py_Dxy_S_a = I_ERI_Hy4z_S_Dxy_S_a+ABY*I_ERI_G4z_S_Dxy_S_a;
  Double I_ERI_G3xz_Pz_Dxy_S_a = I_ERI_H3x2z_S_Dxy_S_a+ABZ*I_ERI_G3xz_S_Dxy_S_a;
  Double I_ERI_G2xyz_Pz_Dxy_S_a = I_ERI_H2xy2z_S_Dxy_S_a+ABZ*I_ERI_G2xyz_S_Dxy_S_a;
  Double I_ERI_G2x2z_Pz_Dxy_S_a = I_ERI_H2x3z_S_Dxy_S_a+ABZ*I_ERI_G2x2z_S_Dxy_S_a;
  Double I_ERI_Gx2yz_Pz_Dxy_S_a = I_ERI_Hx2y2z_S_Dxy_S_a+ABZ*I_ERI_Gx2yz_S_Dxy_S_a;
  Double I_ERI_Gxy2z_Pz_Dxy_S_a = I_ERI_Hxy3z_S_Dxy_S_a+ABZ*I_ERI_Gxy2z_S_Dxy_S_a;
  Double I_ERI_Gx3z_Pz_Dxy_S_a = I_ERI_Hx4z_S_Dxy_S_a+ABZ*I_ERI_Gx3z_S_Dxy_S_a;
  Double I_ERI_G3yz_Pz_Dxy_S_a = I_ERI_H3y2z_S_Dxy_S_a+ABZ*I_ERI_G3yz_S_Dxy_S_a;
  Double I_ERI_G2y2z_Pz_Dxy_S_a = I_ERI_H2y3z_S_Dxy_S_a+ABZ*I_ERI_G2y2z_S_Dxy_S_a;
  Double I_ERI_Gy3z_Pz_Dxy_S_a = I_ERI_Hy4z_S_Dxy_S_a+ABZ*I_ERI_Gy3z_S_Dxy_S_a;
  Double I_ERI_G4z_Pz_Dxy_S_a = I_ERI_H5z_S_Dxy_S_a+ABZ*I_ERI_G4z_S_Dxy_S_a;
  Double I_ERI_G4x_Px_Dxz_S_a = I_ERI_H5x_S_Dxz_S_a+ABX*I_ERI_G4x_S_Dxz_S_a;
  Double I_ERI_G3xy_Px_Dxz_S_a = I_ERI_H4xy_S_Dxz_S_a+ABX*I_ERI_G3xy_S_Dxz_S_a;
  Double I_ERI_G3xz_Px_Dxz_S_a = I_ERI_H4xz_S_Dxz_S_a+ABX*I_ERI_G3xz_S_Dxz_S_a;
  Double I_ERI_G2x2y_Px_Dxz_S_a = I_ERI_H3x2y_S_Dxz_S_a+ABX*I_ERI_G2x2y_S_Dxz_S_a;
  Double I_ERI_G2xyz_Px_Dxz_S_a = I_ERI_H3xyz_S_Dxz_S_a+ABX*I_ERI_G2xyz_S_Dxz_S_a;
  Double I_ERI_G2x2z_Px_Dxz_S_a = I_ERI_H3x2z_S_Dxz_S_a+ABX*I_ERI_G2x2z_S_Dxz_S_a;
  Double I_ERI_Gx3y_Px_Dxz_S_a = I_ERI_H2x3y_S_Dxz_S_a+ABX*I_ERI_Gx3y_S_Dxz_S_a;
  Double I_ERI_Gx2yz_Px_Dxz_S_a = I_ERI_H2x2yz_S_Dxz_S_a+ABX*I_ERI_Gx2yz_S_Dxz_S_a;
  Double I_ERI_Gxy2z_Px_Dxz_S_a = I_ERI_H2xy2z_S_Dxz_S_a+ABX*I_ERI_Gxy2z_S_Dxz_S_a;
  Double I_ERI_Gx3z_Px_Dxz_S_a = I_ERI_H2x3z_S_Dxz_S_a+ABX*I_ERI_Gx3z_S_Dxz_S_a;
  Double I_ERI_G4y_Px_Dxz_S_a = I_ERI_Hx4y_S_Dxz_S_a+ABX*I_ERI_G4y_S_Dxz_S_a;
  Double I_ERI_G3yz_Px_Dxz_S_a = I_ERI_Hx3yz_S_Dxz_S_a+ABX*I_ERI_G3yz_S_Dxz_S_a;
  Double I_ERI_G2y2z_Px_Dxz_S_a = I_ERI_Hx2y2z_S_Dxz_S_a+ABX*I_ERI_G2y2z_S_Dxz_S_a;
  Double I_ERI_Gy3z_Px_Dxz_S_a = I_ERI_Hxy3z_S_Dxz_S_a+ABX*I_ERI_Gy3z_S_Dxz_S_a;
  Double I_ERI_G4z_Px_Dxz_S_a = I_ERI_Hx4z_S_Dxz_S_a+ABX*I_ERI_G4z_S_Dxz_S_a;
  Double I_ERI_G3xy_Py_Dxz_S_a = I_ERI_H3x2y_S_Dxz_S_a+ABY*I_ERI_G3xy_S_Dxz_S_a;
  Double I_ERI_G3xz_Py_Dxz_S_a = I_ERI_H3xyz_S_Dxz_S_a+ABY*I_ERI_G3xz_S_Dxz_S_a;
  Double I_ERI_G2x2y_Py_Dxz_S_a = I_ERI_H2x3y_S_Dxz_S_a+ABY*I_ERI_G2x2y_S_Dxz_S_a;
  Double I_ERI_G2xyz_Py_Dxz_S_a = I_ERI_H2x2yz_S_Dxz_S_a+ABY*I_ERI_G2xyz_S_Dxz_S_a;
  Double I_ERI_G2x2z_Py_Dxz_S_a = I_ERI_H2xy2z_S_Dxz_S_a+ABY*I_ERI_G2x2z_S_Dxz_S_a;
  Double I_ERI_Gx3y_Py_Dxz_S_a = I_ERI_Hx4y_S_Dxz_S_a+ABY*I_ERI_Gx3y_S_Dxz_S_a;
  Double I_ERI_Gx2yz_Py_Dxz_S_a = I_ERI_Hx3yz_S_Dxz_S_a+ABY*I_ERI_Gx2yz_S_Dxz_S_a;
  Double I_ERI_Gxy2z_Py_Dxz_S_a = I_ERI_Hx2y2z_S_Dxz_S_a+ABY*I_ERI_Gxy2z_S_Dxz_S_a;
  Double I_ERI_Gx3z_Py_Dxz_S_a = I_ERI_Hxy3z_S_Dxz_S_a+ABY*I_ERI_Gx3z_S_Dxz_S_a;
  Double I_ERI_G4y_Py_Dxz_S_a = I_ERI_H5y_S_Dxz_S_a+ABY*I_ERI_G4y_S_Dxz_S_a;
  Double I_ERI_G3yz_Py_Dxz_S_a = I_ERI_H4yz_S_Dxz_S_a+ABY*I_ERI_G3yz_S_Dxz_S_a;
  Double I_ERI_G2y2z_Py_Dxz_S_a = I_ERI_H3y2z_S_Dxz_S_a+ABY*I_ERI_G2y2z_S_Dxz_S_a;
  Double I_ERI_Gy3z_Py_Dxz_S_a = I_ERI_H2y3z_S_Dxz_S_a+ABY*I_ERI_Gy3z_S_Dxz_S_a;
  Double I_ERI_G4z_Py_Dxz_S_a = I_ERI_Hy4z_S_Dxz_S_a+ABY*I_ERI_G4z_S_Dxz_S_a;
  Double I_ERI_G3xz_Pz_Dxz_S_a = I_ERI_H3x2z_S_Dxz_S_a+ABZ*I_ERI_G3xz_S_Dxz_S_a;
  Double I_ERI_G2xyz_Pz_Dxz_S_a = I_ERI_H2xy2z_S_Dxz_S_a+ABZ*I_ERI_G2xyz_S_Dxz_S_a;
  Double I_ERI_G2x2z_Pz_Dxz_S_a = I_ERI_H2x3z_S_Dxz_S_a+ABZ*I_ERI_G2x2z_S_Dxz_S_a;
  Double I_ERI_Gx2yz_Pz_Dxz_S_a = I_ERI_Hx2y2z_S_Dxz_S_a+ABZ*I_ERI_Gx2yz_S_Dxz_S_a;
  Double I_ERI_Gxy2z_Pz_Dxz_S_a = I_ERI_Hxy3z_S_Dxz_S_a+ABZ*I_ERI_Gxy2z_S_Dxz_S_a;
  Double I_ERI_Gx3z_Pz_Dxz_S_a = I_ERI_Hx4z_S_Dxz_S_a+ABZ*I_ERI_Gx3z_S_Dxz_S_a;
  Double I_ERI_G3yz_Pz_Dxz_S_a = I_ERI_H3y2z_S_Dxz_S_a+ABZ*I_ERI_G3yz_S_Dxz_S_a;
  Double I_ERI_G2y2z_Pz_Dxz_S_a = I_ERI_H2y3z_S_Dxz_S_a+ABZ*I_ERI_G2y2z_S_Dxz_S_a;
  Double I_ERI_Gy3z_Pz_Dxz_S_a = I_ERI_Hy4z_S_Dxz_S_a+ABZ*I_ERI_Gy3z_S_Dxz_S_a;
  Double I_ERI_G4z_Pz_Dxz_S_a = I_ERI_H5z_S_Dxz_S_a+ABZ*I_ERI_G4z_S_Dxz_S_a;
  Double I_ERI_G4x_Px_D2y_S_a = I_ERI_H5x_S_D2y_S_a+ABX*I_ERI_G4x_S_D2y_S_a;
  Double I_ERI_G3xy_Px_D2y_S_a = I_ERI_H4xy_S_D2y_S_a+ABX*I_ERI_G3xy_S_D2y_S_a;
  Double I_ERI_G3xz_Px_D2y_S_a = I_ERI_H4xz_S_D2y_S_a+ABX*I_ERI_G3xz_S_D2y_S_a;
  Double I_ERI_G2x2y_Px_D2y_S_a = I_ERI_H3x2y_S_D2y_S_a+ABX*I_ERI_G2x2y_S_D2y_S_a;
  Double I_ERI_G2xyz_Px_D2y_S_a = I_ERI_H3xyz_S_D2y_S_a+ABX*I_ERI_G2xyz_S_D2y_S_a;
  Double I_ERI_G2x2z_Px_D2y_S_a = I_ERI_H3x2z_S_D2y_S_a+ABX*I_ERI_G2x2z_S_D2y_S_a;
  Double I_ERI_Gx3y_Px_D2y_S_a = I_ERI_H2x3y_S_D2y_S_a+ABX*I_ERI_Gx3y_S_D2y_S_a;
  Double I_ERI_Gx2yz_Px_D2y_S_a = I_ERI_H2x2yz_S_D2y_S_a+ABX*I_ERI_Gx2yz_S_D2y_S_a;
  Double I_ERI_Gxy2z_Px_D2y_S_a = I_ERI_H2xy2z_S_D2y_S_a+ABX*I_ERI_Gxy2z_S_D2y_S_a;
  Double I_ERI_Gx3z_Px_D2y_S_a = I_ERI_H2x3z_S_D2y_S_a+ABX*I_ERI_Gx3z_S_D2y_S_a;
  Double I_ERI_G4y_Px_D2y_S_a = I_ERI_Hx4y_S_D2y_S_a+ABX*I_ERI_G4y_S_D2y_S_a;
  Double I_ERI_G3yz_Px_D2y_S_a = I_ERI_Hx3yz_S_D2y_S_a+ABX*I_ERI_G3yz_S_D2y_S_a;
  Double I_ERI_G2y2z_Px_D2y_S_a = I_ERI_Hx2y2z_S_D2y_S_a+ABX*I_ERI_G2y2z_S_D2y_S_a;
  Double I_ERI_Gy3z_Px_D2y_S_a = I_ERI_Hxy3z_S_D2y_S_a+ABX*I_ERI_Gy3z_S_D2y_S_a;
  Double I_ERI_G4z_Px_D2y_S_a = I_ERI_Hx4z_S_D2y_S_a+ABX*I_ERI_G4z_S_D2y_S_a;
  Double I_ERI_G3xy_Py_D2y_S_a = I_ERI_H3x2y_S_D2y_S_a+ABY*I_ERI_G3xy_S_D2y_S_a;
  Double I_ERI_G3xz_Py_D2y_S_a = I_ERI_H3xyz_S_D2y_S_a+ABY*I_ERI_G3xz_S_D2y_S_a;
  Double I_ERI_G2x2y_Py_D2y_S_a = I_ERI_H2x3y_S_D2y_S_a+ABY*I_ERI_G2x2y_S_D2y_S_a;
  Double I_ERI_G2xyz_Py_D2y_S_a = I_ERI_H2x2yz_S_D2y_S_a+ABY*I_ERI_G2xyz_S_D2y_S_a;
  Double I_ERI_G2x2z_Py_D2y_S_a = I_ERI_H2xy2z_S_D2y_S_a+ABY*I_ERI_G2x2z_S_D2y_S_a;
  Double I_ERI_Gx3y_Py_D2y_S_a = I_ERI_Hx4y_S_D2y_S_a+ABY*I_ERI_Gx3y_S_D2y_S_a;
  Double I_ERI_Gx2yz_Py_D2y_S_a = I_ERI_Hx3yz_S_D2y_S_a+ABY*I_ERI_Gx2yz_S_D2y_S_a;
  Double I_ERI_Gxy2z_Py_D2y_S_a = I_ERI_Hx2y2z_S_D2y_S_a+ABY*I_ERI_Gxy2z_S_D2y_S_a;
  Double I_ERI_Gx3z_Py_D2y_S_a = I_ERI_Hxy3z_S_D2y_S_a+ABY*I_ERI_Gx3z_S_D2y_S_a;
  Double I_ERI_G4y_Py_D2y_S_a = I_ERI_H5y_S_D2y_S_a+ABY*I_ERI_G4y_S_D2y_S_a;
  Double I_ERI_G3yz_Py_D2y_S_a = I_ERI_H4yz_S_D2y_S_a+ABY*I_ERI_G3yz_S_D2y_S_a;
  Double I_ERI_G2y2z_Py_D2y_S_a = I_ERI_H3y2z_S_D2y_S_a+ABY*I_ERI_G2y2z_S_D2y_S_a;
  Double I_ERI_Gy3z_Py_D2y_S_a = I_ERI_H2y3z_S_D2y_S_a+ABY*I_ERI_Gy3z_S_D2y_S_a;
  Double I_ERI_G4z_Py_D2y_S_a = I_ERI_Hy4z_S_D2y_S_a+ABY*I_ERI_G4z_S_D2y_S_a;
  Double I_ERI_G3xz_Pz_D2y_S_a = I_ERI_H3x2z_S_D2y_S_a+ABZ*I_ERI_G3xz_S_D2y_S_a;
  Double I_ERI_G2xyz_Pz_D2y_S_a = I_ERI_H2xy2z_S_D2y_S_a+ABZ*I_ERI_G2xyz_S_D2y_S_a;
  Double I_ERI_G2x2z_Pz_D2y_S_a = I_ERI_H2x3z_S_D2y_S_a+ABZ*I_ERI_G2x2z_S_D2y_S_a;
  Double I_ERI_Gx2yz_Pz_D2y_S_a = I_ERI_Hx2y2z_S_D2y_S_a+ABZ*I_ERI_Gx2yz_S_D2y_S_a;
  Double I_ERI_Gxy2z_Pz_D2y_S_a = I_ERI_Hxy3z_S_D2y_S_a+ABZ*I_ERI_Gxy2z_S_D2y_S_a;
  Double I_ERI_Gx3z_Pz_D2y_S_a = I_ERI_Hx4z_S_D2y_S_a+ABZ*I_ERI_Gx3z_S_D2y_S_a;
  Double I_ERI_G3yz_Pz_D2y_S_a = I_ERI_H3y2z_S_D2y_S_a+ABZ*I_ERI_G3yz_S_D2y_S_a;
  Double I_ERI_G2y2z_Pz_D2y_S_a = I_ERI_H2y3z_S_D2y_S_a+ABZ*I_ERI_G2y2z_S_D2y_S_a;
  Double I_ERI_Gy3z_Pz_D2y_S_a = I_ERI_Hy4z_S_D2y_S_a+ABZ*I_ERI_Gy3z_S_D2y_S_a;
  Double I_ERI_G4z_Pz_D2y_S_a = I_ERI_H5z_S_D2y_S_a+ABZ*I_ERI_G4z_S_D2y_S_a;
  Double I_ERI_G4x_Px_Dyz_S_a = I_ERI_H5x_S_Dyz_S_a+ABX*I_ERI_G4x_S_Dyz_S_a;
  Double I_ERI_G3xy_Px_Dyz_S_a = I_ERI_H4xy_S_Dyz_S_a+ABX*I_ERI_G3xy_S_Dyz_S_a;
  Double I_ERI_G3xz_Px_Dyz_S_a = I_ERI_H4xz_S_Dyz_S_a+ABX*I_ERI_G3xz_S_Dyz_S_a;
  Double I_ERI_G2x2y_Px_Dyz_S_a = I_ERI_H3x2y_S_Dyz_S_a+ABX*I_ERI_G2x2y_S_Dyz_S_a;
  Double I_ERI_G2xyz_Px_Dyz_S_a = I_ERI_H3xyz_S_Dyz_S_a+ABX*I_ERI_G2xyz_S_Dyz_S_a;
  Double I_ERI_G2x2z_Px_Dyz_S_a = I_ERI_H3x2z_S_Dyz_S_a+ABX*I_ERI_G2x2z_S_Dyz_S_a;
  Double I_ERI_Gx3y_Px_Dyz_S_a = I_ERI_H2x3y_S_Dyz_S_a+ABX*I_ERI_Gx3y_S_Dyz_S_a;
  Double I_ERI_Gx2yz_Px_Dyz_S_a = I_ERI_H2x2yz_S_Dyz_S_a+ABX*I_ERI_Gx2yz_S_Dyz_S_a;
  Double I_ERI_Gxy2z_Px_Dyz_S_a = I_ERI_H2xy2z_S_Dyz_S_a+ABX*I_ERI_Gxy2z_S_Dyz_S_a;
  Double I_ERI_Gx3z_Px_Dyz_S_a = I_ERI_H2x3z_S_Dyz_S_a+ABX*I_ERI_Gx3z_S_Dyz_S_a;
  Double I_ERI_G4y_Px_Dyz_S_a = I_ERI_Hx4y_S_Dyz_S_a+ABX*I_ERI_G4y_S_Dyz_S_a;
  Double I_ERI_G3yz_Px_Dyz_S_a = I_ERI_Hx3yz_S_Dyz_S_a+ABX*I_ERI_G3yz_S_Dyz_S_a;
  Double I_ERI_G2y2z_Px_Dyz_S_a = I_ERI_Hx2y2z_S_Dyz_S_a+ABX*I_ERI_G2y2z_S_Dyz_S_a;
  Double I_ERI_Gy3z_Px_Dyz_S_a = I_ERI_Hxy3z_S_Dyz_S_a+ABX*I_ERI_Gy3z_S_Dyz_S_a;
  Double I_ERI_G4z_Px_Dyz_S_a = I_ERI_Hx4z_S_Dyz_S_a+ABX*I_ERI_G4z_S_Dyz_S_a;
  Double I_ERI_G3xy_Py_Dyz_S_a = I_ERI_H3x2y_S_Dyz_S_a+ABY*I_ERI_G3xy_S_Dyz_S_a;
  Double I_ERI_G3xz_Py_Dyz_S_a = I_ERI_H3xyz_S_Dyz_S_a+ABY*I_ERI_G3xz_S_Dyz_S_a;
  Double I_ERI_G2x2y_Py_Dyz_S_a = I_ERI_H2x3y_S_Dyz_S_a+ABY*I_ERI_G2x2y_S_Dyz_S_a;
  Double I_ERI_G2xyz_Py_Dyz_S_a = I_ERI_H2x2yz_S_Dyz_S_a+ABY*I_ERI_G2xyz_S_Dyz_S_a;
  Double I_ERI_G2x2z_Py_Dyz_S_a = I_ERI_H2xy2z_S_Dyz_S_a+ABY*I_ERI_G2x2z_S_Dyz_S_a;
  Double I_ERI_Gx3y_Py_Dyz_S_a = I_ERI_Hx4y_S_Dyz_S_a+ABY*I_ERI_Gx3y_S_Dyz_S_a;
  Double I_ERI_Gx2yz_Py_Dyz_S_a = I_ERI_Hx3yz_S_Dyz_S_a+ABY*I_ERI_Gx2yz_S_Dyz_S_a;
  Double I_ERI_Gxy2z_Py_Dyz_S_a = I_ERI_Hx2y2z_S_Dyz_S_a+ABY*I_ERI_Gxy2z_S_Dyz_S_a;
  Double I_ERI_Gx3z_Py_Dyz_S_a = I_ERI_Hxy3z_S_Dyz_S_a+ABY*I_ERI_Gx3z_S_Dyz_S_a;
  Double I_ERI_G4y_Py_Dyz_S_a = I_ERI_H5y_S_Dyz_S_a+ABY*I_ERI_G4y_S_Dyz_S_a;
  Double I_ERI_G3yz_Py_Dyz_S_a = I_ERI_H4yz_S_Dyz_S_a+ABY*I_ERI_G3yz_S_Dyz_S_a;
  Double I_ERI_G2y2z_Py_Dyz_S_a = I_ERI_H3y2z_S_Dyz_S_a+ABY*I_ERI_G2y2z_S_Dyz_S_a;
  Double I_ERI_Gy3z_Py_Dyz_S_a = I_ERI_H2y3z_S_Dyz_S_a+ABY*I_ERI_Gy3z_S_Dyz_S_a;
  Double I_ERI_G4z_Py_Dyz_S_a = I_ERI_Hy4z_S_Dyz_S_a+ABY*I_ERI_G4z_S_Dyz_S_a;
  Double I_ERI_G3xz_Pz_Dyz_S_a = I_ERI_H3x2z_S_Dyz_S_a+ABZ*I_ERI_G3xz_S_Dyz_S_a;
  Double I_ERI_G2xyz_Pz_Dyz_S_a = I_ERI_H2xy2z_S_Dyz_S_a+ABZ*I_ERI_G2xyz_S_Dyz_S_a;
  Double I_ERI_G2x2z_Pz_Dyz_S_a = I_ERI_H2x3z_S_Dyz_S_a+ABZ*I_ERI_G2x2z_S_Dyz_S_a;
  Double I_ERI_Gx2yz_Pz_Dyz_S_a = I_ERI_Hx2y2z_S_Dyz_S_a+ABZ*I_ERI_Gx2yz_S_Dyz_S_a;
  Double I_ERI_Gxy2z_Pz_Dyz_S_a = I_ERI_Hxy3z_S_Dyz_S_a+ABZ*I_ERI_Gxy2z_S_Dyz_S_a;
  Double I_ERI_Gx3z_Pz_Dyz_S_a = I_ERI_Hx4z_S_Dyz_S_a+ABZ*I_ERI_Gx3z_S_Dyz_S_a;
  Double I_ERI_G3yz_Pz_Dyz_S_a = I_ERI_H3y2z_S_Dyz_S_a+ABZ*I_ERI_G3yz_S_Dyz_S_a;
  Double I_ERI_G2y2z_Pz_Dyz_S_a = I_ERI_H2y3z_S_Dyz_S_a+ABZ*I_ERI_G2y2z_S_Dyz_S_a;
  Double I_ERI_Gy3z_Pz_Dyz_S_a = I_ERI_Hy4z_S_Dyz_S_a+ABZ*I_ERI_Gy3z_S_Dyz_S_a;
  Double I_ERI_G4z_Pz_Dyz_S_a = I_ERI_H5z_S_Dyz_S_a+ABZ*I_ERI_G4z_S_Dyz_S_a;
  Double I_ERI_G4x_Px_D2z_S_a = I_ERI_H5x_S_D2z_S_a+ABX*I_ERI_G4x_S_D2z_S_a;
  Double I_ERI_G3xy_Px_D2z_S_a = I_ERI_H4xy_S_D2z_S_a+ABX*I_ERI_G3xy_S_D2z_S_a;
  Double I_ERI_G3xz_Px_D2z_S_a = I_ERI_H4xz_S_D2z_S_a+ABX*I_ERI_G3xz_S_D2z_S_a;
  Double I_ERI_G2x2y_Px_D2z_S_a = I_ERI_H3x2y_S_D2z_S_a+ABX*I_ERI_G2x2y_S_D2z_S_a;
  Double I_ERI_G2xyz_Px_D2z_S_a = I_ERI_H3xyz_S_D2z_S_a+ABX*I_ERI_G2xyz_S_D2z_S_a;
  Double I_ERI_G2x2z_Px_D2z_S_a = I_ERI_H3x2z_S_D2z_S_a+ABX*I_ERI_G2x2z_S_D2z_S_a;
  Double I_ERI_Gx3y_Px_D2z_S_a = I_ERI_H2x3y_S_D2z_S_a+ABX*I_ERI_Gx3y_S_D2z_S_a;
  Double I_ERI_Gx2yz_Px_D2z_S_a = I_ERI_H2x2yz_S_D2z_S_a+ABX*I_ERI_Gx2yz_S_D2z_S_a;
  Double I_ERI_Gxy2z_Px_D2z_S_a = I_ERI_H2xy2z_S_D2z_S_a+ABX*I_ERI_Gxy2z_S_D2z_S_a;
  Double I_ERI_Gx3z_Px_D2z_S_a = I_ERI_H2x3z_S_D2z_S_a+ABX*I_ERI_Gx3z_S_D2z_S_a;
  Double I_ERI_G4y_Px_D2z_S_a = I_ERI_Hx4y_S_D2z_S_a+ABX*I_ERI_G4y_S_D2z_S_a;
  Double I_ERI_G3yz_Px_D2z_S_a = I_ERI_Hx3yz_S_D2z_S_a+ABX*I_ERI_G3yz_S_D2z_S_a;
  Double I_ERI_G2y2z_Px_D2z_S_a = I_ERI_Hx2y2z_S_D2z_S_a+ABX*I_ERI_G2y2z_S_D2z_S_a;
  Double I_ERI_Gy3z_Px_D2z_S_a = I_ERI_Hxy3z_S_D2z_S_a+ABX*I_ERI_Gy3z_S_D2z_S_a;
  Double I_ERI_G4z_Px_D2z_S_a = I_ERI_Hx4z_S_D2z_S_a+ABX*I_ERI_G4z_S_D2z_S_a;
  Double I_ERI_G3xy_Py_D2z_S_a = I_ERI_H3x2y_S_D2z_S_a+ABY*I_ERI_G3xy_S_D2z_S_a;
  Double I_ERI_G3xz_Py_D2z_S_a = I_ERI_H3xyz_S_D2z_S_a+ABY*I_ERI_G3xz_S_D2z_S_a;
  Double I_ERI_G2x2y_Py_D2z_S_a = I_ERI_H2x3y_S_D2z_S_a+ABY*I_ERI_G2x2y_S_D2z_S_a;
  Double I_ERI_G2xyz_Py_D2z_S_a = I_ERI_H2x2yz_S_D2z_S_a+ABY*I_ERI_G2xyz_S_D2z_S_a;
  Double I_ERI_G2x2z_Py_D2z_S_a = I_ERI_H2xy2z_S_D2z_S_a+ABY*I_ERI_G2x2z_S_D2z_S_a;
  Double I_ERI_Gx3y_Py_D2z_S_a = I_ERI_Hx4y_S_D2z_S_a+ABY*I_ERI_Gx3y_S_D2z_S_a;
  Double I_ERI_Gx2yz_Py_D2z_S_a = I_ERI_Hx3yz_S_D2z_S_a+ABY*I_ERI_Gx2yz_S_D2z_S_a;
  Double I_ERI_Gxy2z_Py_D2z_S_a = I_ERI_Hx2y2z_S_D2z_S_a+ABY*I_ERI_Gxy2z_S_D2z_S_a;
  Double I_ERI_Gx3z_Py_D2z_S_a = I_ERI_Hxy3z_S_D2z_S_a+ABY*I_ERI_Gx3z_S_D2z_S_a;
  Double I_ERI_G4y_Py_D2z_S_a = I_ERI_H5y_S_D2z_S_a+ABY*I_ERI_G4y_S_D2z_S_a;
  Double I_ERI_G3yz_Py_D2z_S_a = I_ERI_H4yz_S_D2z_S_a+ABY*I_ERI_G3yz_S_D2z_S_a;
  Double I_ERI_G2y2z_Py_D2z_S_a = I_ERI_H3y2z_S_D2z_S_a+ABY*I_ERI_G2y2z_S_D2z_S_a;
  Double I_ERI_Gy3z_Py_D2z_S_a = I_ERI_H2y3z_S_D2z_S_a+ABY*I_ERI_Gy3z_S_D2z_S_a;
  Double I_ERI_G4z_Py_D2z_S_a = I_ERI_Hy4z_S_D2z_S_a+ABY*I_ERI_G4z_S_D2z_S_a;
  Double I_ERI_G3xz_Pz_D2z_S_a = I_ERI_H3x2z_S_D2z_S_a+ABZ*I_ERI_G3xz_S_D2z_S_a;
  Double I_ERI_G2xyz_Pz_D2z_S_a = I_ERI_H2xy2z_S_D2z_S_a+ABZ*I_ERI_G2xyz_S_D2z_S_a;
  Double I_ERI_G2x2z_Pz_D2z_S_a = I_ERI_H2x3z_S_D2z_S_a+ABZ*I_ERI_G2x2z_S_D2z_S_a;
  Double I_ERI_Gx2yz_Pz_D2z_S_a = I_ERI_Hx2y2z_S_D2z_S_a+ABZ*I_ERI_Gx2yz_S_D2z_S_a;
  Double I_ERI_Gxy2z_Pz_D2z_S_a = I_ERI_Hxy3z_S_D2z_S_a+ABZ*I_ERI_Gxy2z_S_D2z_S_a;
  Double I_ERI_Gx3z_Pz_D2z_S_a = I_ERI_Hx4z_S_D2z_S_a+ABZ*I_ERI_Gx3z_S_D2z_S_a;
  Double I_ERI_G3yz_Pz_D2z_S_a = I_ERI_H3y2z_S_D2z_S_a+ABZ*I_ERI_G3yz_S_D2z_S_a;
  Double I_ERI_G2y2z_Pz_D2z_S_a = I_ERI_H2y3z_S_D2z_S_a+ABZ*I_ERI_G2y2z_S_D2z_S_a;
  Double I_ERI_Gy3z_Pz_D2z_S_a = I_ERI_Hy4z_S_D2z_S_a+ABZ*I_ERI_Gy3z_S_D2z_S_a;
  Double I_ERI_G4z_Pz_D2z_S_a = I_ERI_H5z_S_D2z_S_a+ABZ*I_ERI_G4z_S_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_D_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_D_S_a
   * RHS shell quartet name: SQ_ERI_F_P_D_S_a
   ************************************************************/
  Double I_ERI_F3x_D2x_D2x_S_a = I_ERI_G4x_Px_D2x_S_a+ABX*I_ERI_F3x_Px_D2x_S_a;
  Double I_ERI_F2xy_D2x_D2x_S_a = I_ERI_G3xy_Px_D2x_S_a+ABX*I_ERI_F2xy_Px_D2x_S_a;
  Double I_ERI_F2xz_D2x_D2x_S_a = I_ERI_G3xz_Px_D2x_S_a+ABX*I_ERI_F2xz_Px_D2x_S_a;
  Double I_ERI_Fx2y_D2x_D2x_S_a = I_ERI_G2x2y_Px_D2x_S_a+ABX*I_ERI_Fx2y_Px_D2x_S_a;
  Double I_ERI_Fxyz_D2x_D2x_S_a = I_ERI_G2xyz_Px_D2x_S_a+ABX*I_ERI_Fxyz_Px_D2x_S_a;
  Double I_ERI_Fx2z_D2x_D2x_S_a = I_ERI_G2x2z_Px_D2x_S_a+ABX*I_ERI_Fx2z_Px_D2x_S_a;
  Double I_ERI_F3y_D2x_D2x_S_a = I_ERI_Gx3y_Px_D2x_S_a+ABX*I_ERI_F3y_Px_D2x_S_a;
  Double I_ERI_F2yz_D2x_D2x_S_a = I_ERI_Gx2yz_Px_D2x_S_a+ABX*I_ERI_F2yz_Px_D2x_S_a;
  Double I_ERI_Fy2z_D2x_D2x_S_a = I_ERI_Gxy2z_Px_D2x_S_a+ABX*I_ERI_Fy2z_Px_D2x_S_a;
  Double I_ERI_F3z_D2x_D2x_S_a = I_ERI_Gx3z_Px_D2x_S_a+ABX*I_ERI_F3z_Px_D2x_S_a;
  Double I_ERI_F3x_Dxy_D2x_S_a = I_ERI_G3xy_Px_D2x_S_a+ABY*I_ERI_F3x_Px_D2x_S_a;
  Double I_ERI_F2xy_Dxy_D2x_S_a = I_ERI_G2x2y_Px_D2x_S_a+ABY*I_ERI_F2xy_Px_D2x_S_a;
  Double I_ERI_F2xz_Dxy_D2x_S_a = I_ERI_G2xyz_Px_D2x_S_a+ABY*I_ERI_F2xz_Px_D2x_S_a;
  Double I_ERI_Fx2y_Dxy_D2x_S_a = I_ERI_Gx3y_Px_D2x_S_a+ABY*I_ERI_Fx2y_Px_D2x_S_a;
  Double I_ERI_Fxyz_Dxy_D2x_S_a = I_ERI_Gx2yz_Px_D2x_S_a+ABY*I_ERI_Fxyz_Px_D2x_S_a;
  Double I_ERI_Fx2z_Dxy_D2x_S_a = I_ERI_Gxy2z_Px_D2x_S_a+ABY*I_ERI_Fx2z_Px_D2x_S_a;
  Double I_ERI_F3y_Dxy_D2x_S_a = I_ERI_G4y_Px_D2x_S_a+ABY*I_ERI_F3y_Px_D2x_S_a;
  Double I_ERI_F2yz_Dxy_D2x_S_a = I_ERI_G3yz_Px_D2x_S_a+ABY*I_ERI_F2yz_Px_D2x_S_a;
  Double I_ERI_Fy2z_Dxy_D2x_S_a = I_ERI_G2y2z_Px_D2x_S_a+ABY*I_ERI_Fy2z_Px_D2x_S_a;
  Double I_ERI_F3z_Dxy_D2x_S_a = I_ERI_Gy3z_Px_D2x_S_a+ABY*I_ERI_F3z_Px_D2x_S_a;
  Double I_ERI_F3x_Dxz_D2x_S_a = I_ERI_G3xz_Px_D2x_S_a+ABZ*I_ERI_F3x_Px_D2x_S_a;
  Double I_ERI_F2xy_Dxz_D2x_S_a = I_ERI_G2xyz_Px_D2x_S_a+ABZ*I_ERI_F2xy_Px_D2x_S_a;
  Double I_ERI_F2xz_Dxz_D2x_S_a = I_ERI_G2x2z_Px_D2x_S_a+ABZ*I_ERI_F2xz_Px_D2x_S_a;
  Double I_ERI_Fx2y_Dxz_D2x_S_a = I_ERI_Gx2yz_Px_D2x_S_a+ABZ*I_ERI_Fx2y_Px_D2x_S_a;
  Double I_ERI_Fxyz_Dxz_D2x_S_a = I_ERI_Gxy2z_Px_D2x_S_a+ABZ*I_ERI_Fxyz_Px_D2x_S_a;
  Double I_ERI_Fx2z_Dxz_D2x_S_a = I_ERI_Gx3z_Px_D2x_S_a+ABZ*I_ERI_Fx2z_Px_D2x_S_a;
  Double I_ERI_F3y_Dxz_D2x_S_a = I_ERI_G3yz_Px_D2x_S_a+ABZ*I_ERI_F3y_Px_D2x_S_a;
  Double I_ERI_F2yz_Dxz_D2x_S_a = I_ERI_G2y2z_Px_D2x_S_a+ABZ*I_ERI_F2yz_Px_D2x_S_a;
  Double I_ERI_Fy2z_Dxz_D2x_S_a = I_ERI_Gy3z_Px_D2x_S_a+ABZ*I_ERI_Fy2z_Px_D2x_S_a;
  Double I_ERI_F3z_Dxz_D2x_S_a = I_ERI_G4z_Px_D2x_S_a+ABZ*I_ERI_F3z_Px_D2x_S_a;
  Double I_ERI_F3x_D2y_D2x_S_a = I_ERI_G3xy_Py_D2x_S_a+ABY*I_ERI_F3x_Py_D2x_S_a;
  Double I_ERI_F2xy_D2y_D2x_S_a = I_ERI_G2x2y_Py_D2x_S_a+ABY*I_ERI_F2xy_Py_D2x_S_a;
  Double I_ERI_F2xz_D2y_D2x_S_a = I_ERI_G2xyz_Py_D2x_S_a+ABY*I_ERI_F2xz_Py_D2x_S_a;
  Double I_ERI_Fx2y_D2y_D2x_S_a = I_ERI_Gx3y_Py_D2x_S_a+ABY*I_ERI_Fx2y_Py_D2x_S_a;
  Double I_ERI_Fxyz_D2y_D2x_S_a = I_ERI_Gx2yz_Py_D2x_S_a+ABY*I_ERI_Fxyz_Py_D2x_S_a;
  Double I_ERI_Fx2z_D2y_D2x_S_a = I_ERI_Gxy2z_Py_D2x_S_a+ABY*I_ERI_Fx2z_Py_D2x_S_a;
  Double I_ERI_F3y_D2y_D2x_S_a = I_ERI_G4y_Py_D2x_S_a+ABY*I_ERI_F3y_Py_D2x_S_a;
  Double I_ERI_F2yz_D2y_D2x_S_a = I_ERI_G3yz_Py_D2x_S_a+ABY*I_ERI_F2yz_Py_D2x_S_a;
  Double I_ERI_Fy2z_D2y_D2x_S_a = I_ERI_G2y2z_Py_D2x_S_a+ABY*I_ERI_Fy2z_Py_D2x_S_a;
  Double I_ERI_F3z_D2y_D2x_S_a = I_ERI_Gy3z_Py_D2x_S_a+ABY*I_ERI_F3z_Py_D2x_S_a;
  Double I_ERI_F3x_Dyz_D2x_S_a = I_ERI_G3xz_Py_D2x_S_a+ABZ*I_ERI_F3x_Py_D2x_S_a;
  Double I_ERI_F2xy_Dyz_D2x_S_a = I_ERI_G2xyz_Py_D2x_S_a+ABZ*I_ERI_F2xy_Py_D2x_S_a;
  Double I_ERI_F2xz_Dyz_D2x_S_a = I_ERI_G2x2z_Py_D2x_S_a+ABZ*I_ERI_F2xz_Py_D2x_S_a;
  Double I_ERI_Fx2y_Dyz_D2x_S_a = I_ERI_Gx2yz_Py_D2x_S_a+ABZ*I_ERI_Fx2y_Py_D2x_S_a;
  Double I_ERI_Fxyz_Dyz_D2x_S_a = I_ERI_Gxy2z_Py_D2x_S_a+ABZ*I_ERI_Fxyz_Py_D2x_S_a;
  Double I_ERI_Fx2z_Dyz_D2x_S_a = I_ERI_Gx3z_Py_D2x_S_a+ABZ*I_ERI_Fx2z_Py_D2x_S_a;
  Double I_ERI_F3y_Dyz_D2x_S_a = I_ERI_G3yz_Py_D2x_S_a+ABZ*I_ERI_F3y_Py_D2x_S_a;
  Double I_ERI_F2yz_Dyz_D2x_S_a = I_ERI_G2y2z_Py_D2x_S_a+ABZ*I_ERI_F2yz_Py_D2x_S_a;
  Double I_ERI_Fy2z_Dyz_D2x_S_a = I_ERI_Gy3z_Py_D2x_S_a+ABZ*I_ERI_Fy2z_Py_D2x_S_a;
  Double I_ERI_F3z_Dyz_D2x_S_a = I_ERI_G4z_Py_D2x_S_a+ABZ*I_ERI_F3z_Py_D2x_S_a;
  Double I_ERI_F3x_D2z_D2x_S_a = I_ERI_G3xz_Pz_D2x_S_a+ABZ*I_ERI_F3x_Pz_D2x_S_a;
  Double I_ERI_F2xy_D2z_D2x_S_a = I_ERI_G2xyz_Pz_D2x_S_a+ABZ*I_ERI_F2xy_Pz_D2x_S_a;
  Double I_ERI_F2xz_D2z_D2x_S_a = I_ERI_G2x2z_Pz_D2x_S_a+ABZ*I_ERI_F2xz_Pz_D2x_S_a;
  Double I_ERI_Fx2y_D2z_D2x_S_a = I_ERI_Gx2yz_Pz_D2x_S_a+ABZ*I_ERI_Fx2y_Pz_D2x_S_a;
  Double I_ERI_Fxyz_D2z_D2x_S_a = I_ERI_Gxy2z_Pz_D2x_S_a+ABZ*I_ERI_Fxyz_Pz_D2x_S_a;
  Double I_ERI_Fx2z_D2z_D2x_S_a = I_ERI_Gx3z_Pz_D2x_S_a+ABZ*I_ERI_Fx2z_Pz_D2x_S_a;
  Double I_ERI_F3y_D2z_D2x_S_a = I_ERI_G3yz_Pz_D2x_S_a+ABZ*I_ERI_F3y_Pz_D2x_S_a;
  Double I_ERI_F2yz_D2z_D2x_S_a = I_ERI_G2y2z_Pz_D2x_S_a+ABZ*I_ERI_F2yz_Pz_D2x_S_a;
  Double I_ERI_Fy2z_D2z_D2x_S_a = I_ERI_Gy3z_Pz_D2x_S_a+ABZ*I_ERI_Fy2z_Pz_D2x_S_a;
  Double I_ERI_F3z_D2z_D2x_S_a = I_ERI_G4z_Pz_D2x_S_a+ABZ*I_ERI_F3z_Pz_D2x_S_a;
  Double I_ERI_F3x_D2x_Dxy_S_a = I_ERI_G4x_Px_Dxy_S_a+ABX*I_ERI_F3x_Px_Dxy_S_a;
  Double I_ERI_F2xy_D2x_Dxy_S_a = I_ERI_G3xy_Px_Dxy_S_a+ABX*I_ERI_F2xy_Px_Dxy_S_a;
  Double I_ERI_F2xz_D2x_Dxy_S_a = I_ERI_G3xz_Px_Dxy_S_a+ABX*I_ERI_F2xz_Px_Dxy_S_a;
  Double I_ERI_Fx2y_D2x_Dxy_S_a = I_ERI_G2x2y_Px_Dxy_S_a+ABX*I_ERI_Fx2y_Px_Dxy_S_a;
  Double I_ERI_Fxyz_D2x_Dxy_S_a = I_ERI_G2xyz_Px_Dxy_S_a+ABX*I_ERI_Fxyz_Px_Dxy_S_a;
  Double I_ERI_Fx2z_D2x_Dxy_S_a = I_ERI_G2x2z_Px_Dxy_S_a+ABX*I_ERI_Fx2z_Px_Dxy_S_a;
  Double I_ERI_F3y_D2x_Dxy_S_a = I_ERI_Gx3y_Px_Dxy_S_a+ABX*I_ERI_F3y_Px_Dxy_S_a;
  Double I_ERI_F2yz_D2x_Dxy_S_a = I_ERI_Gx2yz_Px_Dxy_S_a+ABX*I_ERI_F2yz_Px_Dxy_S_a;
  Double I_ERI_Fy2z_D2x_Dxy_S_a = I_ERI_Gxy2z_Px_Dxy_S_a+ABX*I_ERI_Fy2z_Px_Dxy_S_a;
  Double I_ERI_F3z_D2x_Dxy_S_a = I_ERI_Gx3z_Px_Dxy_S_a+ABX*I_ERI_F3z_Px_Dxy_S_a;
  Double I_ERI_F3x_Dxy_Dxy_S_a = I_ERI_G3xy_Px_Dxy_S_a+ABY*I_ERI_F3x_Px_Dxy_S_a;
  Double I_ERI_F2xy_Dxy_Dxy_S_a = I_ERI_G2x2y_Px_Dxy_S_a+ABY*I_ERI_F2xy_Px_Dxy_S_a;
  Double I_ERI_F2xz_Dxy_Dxy_S_a = I_ERI_G2xyz_Px_Dxy_S_a+ABY*I_ERI_F2xz_Px_Dxy_S_a;
  Double I_ERI_Fx2y_Dxy_Dxy_S_a = I_ERI_Gx3y_Px_Dxy_S_a+ABY*I_ERI_Fx2y_Px_Dxy_S_a;
  Double I_ERI_Fxyz_Dxy_Dxy_S_a = I_ERI_Gx2yz_Px_Dxy_S_a+ABY*I_ERI_Fxyz_Px_Dxy_S_a;
  Double I_ERI_Fx2z_Dxy_Dxy_S_a = I_ERI_Gxy2z_Px_Dxy_S_a+ABY*I_ERI_Fx2z_Px_Dxy_S_a;
  Double I_ERI_F3y_Dxy_Dxy_S_a = I_ERI_G4y_Px_Dxy_S_a+ABY*I_ERI_F3y_Px_Dxy_S_a;
  Double I_ERI_F2yz_Dxy_Dxy_S_a = I_ERI_G3yz_Px_Dxy_S_a+ABY*I_ERI_F2yz_Px_Dxy_S_a;
  Double I_ERI_Fy2z_Dxy_Dxy_S_a = I_ERI_G2y2z_Px_Dxy_S_a+ABY*I_ERI_Fy2z_Px_Dxy_S_a;
  Double I_ERI_F3z_Dxy_Dxy_S_a = I_ERI_Gy3z_Px_Dxy_S_a+ABY*I_ERI_F3z_Px_Dxy_S_a;
  Double I_ERI_F3x_Dxz_Dxy_S_a = I_ERI_G3xz_Px_Dxy_S_a+ABZ*I_ERI_F3x_Px_Dxy_S_a;
  Double I_ERI_F2xy_Dxz_Dxy_S_a = I_ERI_G2xyz_Px_Dxy_S_a+ABZ*I_ERI_F2xy_Px_Dxy_S_a;
  Double I_ERI_F2xz_Dxz_Dxy_S_a = I_ERI_G2x2z_Px_Dxy_S_a+ABZ*I_ERI_F2xz_Px_Dxy_S_a;
  Double I_ERI_Fx2y_Dxz_Dxy_S_a = I_ERI_Gx2yz_Px_Dxy_S_a+ABZ*I_ERI_Fx2y_Px_Dxy_S_a;
  Double I_ERI_Fxyz_Dxz_Dxy_S_a = I_ERI_Gxy2z_Px_Dxy_S_a+ABZ*I_ERI_Fxyz_Px_Dxy_S_a;
  Double I_ERI_Fx2z_Dxz_Dxy_S_a = I_ERI_Gx3z_Px_Dxy_S_a+ABZ*I_ERI_Fx2z_Px_Dxy_S_a;
  Double I_ERI_F3y_Dxz_Dxy_S_a = I_ERI_G3yz_Px_Dxy_S_a+ABZ*I_ERI_F3y_Px_Dxy_S_a;
  Double I_ERI_F2yz_Dxz_Dxy_S_a = I_ERI_G2y2z_Px_Dxy_S_a+ABZ*I_ERI_F2yz_Px_Dxy_S_a;
  Double I_ERI_Fy2z_Dxz_Dxy_S_a = I_ERI_Gy3z_Px_Dxy_S_a+ABZ*I_ERI_Fy2z_Px_Dxy_S_a;
  Double I_ERI_F3z_Dxz_Dxy_S_a = I_ERI_G4z_Px_Dxy_S_a+ABZ*I_ERI_F3z_Px_Dxy_S_a;
  Double I_ERI_F3x_D2y_Dxy_S_a = I_ERI_G3xy_Py_Dxy_S_a+ABY*I_ERI_F3x_Py_Dxy_S_a;
  Double I_ERI_F2xy_D2y_Dxy_S_a = I_ERI_G2x2y_Py_Dxy_S_a+ABY*I_ERI_F2xy_Py_Dxy_S_a;
  Double I_ERI_F2xz_D2y_Dxy_S_a = I_ERI_G2xyz_Py_Dxy_S_a+ABY*I_ERI_F2xz_Py_Dxy_S_a;
  Double I_ERI_Fx2y_D2y_Dxy_S_a = I_ERI_Gx3y_Py_Dxy_S_a+ABY*I_ERI_Fx2y_Py_Dxy_S_a;
  Double I_ERI_Fxyz_D2y_Dxy_S_a = I_ERI_Gx2yz_Py_Dxy_S_a+ABY*I_ERI_Fxyz_Py_Dxy_S_a;
  Double I_ERI_Fx2z_D2y_Dxy_S_a = I_ERI_Gxy2z_Py_Dxy_S_a+ABY*I_ERI_Fx2z_Py_Dxy_S_a;
  Double I_ERI_F3y_D2y_Dxy_S_a = I_ERI_G4y_Py_Dxy_S_a+ABY*I_ERI_F3y_Py_Dxy_S_a;
  Double I_ERI_F2yz_D2y_Dxy_S_a = I_ERI_G3yz_Py_Dxy_S_a+ABY*I_ERI_F2yz_Py_Dxy_S_a;
  Double I_ERI_Fy2z_D2y_Dxy_S_a = I_ERI_G2y2z_Py_Dxy_S_a+ABY*I_ERI_Fy2z_Py_Dxy_S_a;
  Double I_ERI_F3z_D2y_Dxy_S_a = I_ERI_Gy3z_Py_Dxy_S_a+ABY*I_ERI_F3z_Py_Dxy_S_a;
  Double I_ERI_F3x_Dyz_Dxy_S_a = I_ERI_G3xz_Py_Dxy_S_a+ABZ*I_ERI_F3x_Py_Dxy_S_a;
  Double I_ERI_F2xy_Dyz_Dxy_S_a = I_ERI_G2xyz_Py_Dxy_S_a+ABZ*I_ERI_F2xy_Py_Dxy_S_a;
  Double I_ERI_F2xz_Dyz_Dxy_S_a = I_ERI_G2x2z_Py_Dxy_S_a+ABZ*I_ERI_F2xz_Py_Dxy_S_a;
  Double I_ERI_Fx2y_Dyz_Dxy_S_a = I_ERI_Gx2yz_Py_Dxy_S_a+ABZ*I_ERI_Fx2y_Py_Dxy_S_a;
  Double I_ERI_Fxyz_Dyz_Dxy_S_a = I_ERI_Gxy2z_Py_Dxy_S_a+ABZ*I_ERI_Fxyz_Py_Dxy_S_a;
  Double I_ERI_Fx2z_Dyz_Dxy_S_a = I_ERI_Gx3z_Py_Dxy_S_a+ABZ*I_ERI_Fx2z_Py_Dxy_S_a;
  Double I_ERI_F3y_Dyz_Dxy_S_a = I_ERI_G3yz_Py_Dxy_S_a+ABZ*I_ERI_F3y_Py_Dxy_S_a;
  Double I_ERI_F2yz_Dyz_Dxy_S_a = I_ERI_G2y2z_Py_Dxy_S_a+ABZ*I_ERI_F2yz_Py_Dxy_S_a;
  Double I_ERI_Fy2z_Dyz_Dxy_S_a = I_ERI_Gy3z_Py_Dxy_S_a+ABZ*I_ERI_Fy2z_Py_Dxy_S_a;
  Double I_ERI_F3z_Dyz_Dxy_S_a = I_ERI_G4z_Py_Dxy_S_a+ABZ*I_ERI_F3z_Py_Dxy_S_a;
  Double I_ERI_F3x_D2z_Dxy_S_a = I_ERI_G3xz_Pz_Dxy_S_a+ABZ*I_ERI_F3x_Pz_Dxy_S_a;
  Double I_ERI_F2xy_D2z_Dxy_S_a = I_ERI_G2xyz_Pz_Dxy_S_a+ABZ*I_ERI_F2xy_Pz_Dxy_S_a;
  Double I_ERI_F2xz_D2z_Dxy_S_a = I_ERI_G2x2z_Pz_Dxy_S_a+ABZ*I_ERI_F2xz_Pz_Dxy_S_a;
  Double I_ERI_Fx2y_D2z_Dxy_S_a = I_ERI_Gx2yz_Pz_Dxy_S_a+ABZ*I_ERI_Fx2y_Pz_Dxy_S_a;
  Double I_ERI_Fxyz_D2z_Dxy_S_a = I_ERI_Gxy2z_Pz_Dxy_S_a+ABZ*I_ERI_Fxyz_Pz_Dxy_S_a;
  Double I_ERI_Fx2z_D2z_Dxy_S_a = I_ERI_Gx3z_Pz_Dxy_S_a+ABZ*I_ERI_Fx2z_Pz_Dxy_S_a;
  Double I_ERI_F3y_D2z_Dxy_S_a = I_ERI_G3yz_Pz_Dxy_S_a+ABZ*I_ERI_F3y_Pz_Dxy_S_a;
  Double I_ERI_F2yz_D2z_Dxy_S_a = I_ERI_G2y2z_Pz_Dxy_S_a+ABZ*I_ERI_F2yz_Pz_Dxy_S_a;
  Double I_ERI_Fy2z_D2z_Dxy_S_a = I_ERI_Gy3z_Pz_Dxy_S_a+ABZ*I_ERI_Fy2z_Pz_Dxy_S_a;
  Double I_ERI_F3z_D2z_Dxy_S_a = I_ERI_G4z_Pz_Dxy_S_a+ABZ*I_ERI_F3z_Pz_Dxy_S_a;
  Double I_ERI_F3x_D2x_Dxz_S_a = I_ERI_G4x_Px_Dxz_S_a+ABX*I_ERI_F3x_Px_Dxz_S_a;
  Double I_ERI_F2xy_D2x_Dxz_S_a = I_ERI_G3xy_Px_Dxz_S_a+ABX*I_ERI_F2xy_Px_Dxz_S_a;
  Double I_ERI_F2xz_D2x_Dxz_S_a = I_ERI_G3xz_Px_Dxz_S_a+ABX*I_ERI_F2xz_Px_Dxz_S_a;
  Double I_ERI_Fx2y_D2x_Dxz_S_a = I_ERI_G2x2y_Px_Dxz_S_a+ABX*I_ERI_Fx2y_Px_Dxz_S_a;
  Double I_ERI_Fxyz_D2x_Dxz_S_a = I_ERI_G2xyz_Px_Dxz_S_a+ABX*I_ERI_Fxyz_Px_Dxz_S_a;
  Double I_ERI_Fx2z_D2x_Dxz_S_a = I_ERI_G2x2z_Px_Dxz_S_a+ABX*I_ERI_Fx2z_Px_Dxz_S_a;
  Double I_ERI_F3y_D2x_Dxz_S_a = I_ERI_Gx3y_Px_Dxz_S_a+ABX*I_ERI_F3y_Px_Dxz_S_a;
  Double I_ERI_F2yz_D2x_Dxz_S_a = I_ERI_Gx2yz_Px_Dxz_S_a+ABX*I_ERI_F2yz_Px_Dxz_S_a;
  Double I_ERI_Fy2z_D2x_Dxz_S_a = I_ERI_Gxy2z_Px_Dxz_S_a+ABX*I_ERI_Fy2z_Px_Dxz_S_a;
  Double I_ERI_F3z_D2x_Dxz_S_a = I_ERI_Gx3z_Px_Dxz_S_a+ABX*I_ERI_F3z_Px_Dxz_S_a;
  Double I_ERI_F3x_Dxy_Dxz_S_a = I_ERI_G3xy_Px_Dxz_S_a+ABY*I_ERI_F3x_Px_Dxz_S_a;
  Double I_ERI_F2xy_Dxy_Dxz_S_a = I_ERI_G2x2y_Px_Dxz_S_a+ABY*I_ERI_F2xy_Px_Dxz_S_a;
  Double I_ERI_F2xz_Dxy_Dxz_S_a = I_ERI_G2xyz_Px_Dxz_S_a+ABY*I_ERI_F2xz_Px_Dxz_S_a;
  Double I_ERI_Fx2y_Dxy_Dxz_S_a = I_ERI_Gx3y_Px_Dxz_S_a+ABY*I_ERI_Fx2y_Px_Dxz_S_a;
  Double I_ERI_Fxyz_Dxy_Dxz_S_a = I_ERI_Gx2yz_Px_Dxz_S_a+ABY*I_ERI_Fxyz_Px_Dxz_S_a;
  Double I_ERI_Fx2z_Dxy_Dxz_S_a = I_ERI_Gxy2z_Px_Dxz_S_a+ABY*I_ERI_Fx2z_Px_Dxz_S_a;
  Double I_ERI_F3y_Dxy_Dxz_S_a = I_ERI_G4y_Px_Dxz_S_a+ABY*I_ERI_F3y_Px_Dxz_S_a;
  Double I_ERI_F2yz_Dxy_Dxz_S_a = I_ERI_G3yz_Px_Dxz_S_a+ABY*I_ERI_F2yz_Px_Dxz_S_a;
  Double I_ERI_Fy2z_Dxy_Dxz_S_a = I_ERI_G2y2z_Px_Dxz_S_a+ABY*I_ERI_Fy2z_Px_Dxz_S_a;
  Double I_ERI_F3z_Dxy_Dxz_S_a = I_ERI_Gy3z_Px_Dxz_S_a+ABY*I_ERI_F3z_Px_Dxz_S_a;
  Double I_ERI_F3x_Dxz_Dxz_S_a = I_ERI_G3xz_Px_Dxz_S_a+ABZ*I_ERI_F3x_Px_Dxz_S_a;
  Double I_ERI_F2xy_Dxz_Dxz_S_a = I_ERI_G2xyz_Px_Dxz_S_a+ABZ*I_ERI_F2xy_Px_Dxz_S_a;
  Double I_ERI_F2xz_Dxz_Dxz_S_a = I_ERI_G2x2z_Px_Dxz_S_a+ABZ*I_ERI_F2xz_Px_Dxz_S_a;
  Double I_ERI_Fx2y_Dxz_Dxz_S_a = I_ERI_Gx2yz_Px_Dxz_S_a+ABZ*I_ERI_Fx2y_Px_Dxz_S_a;
  Double I_ERI_Fxyz_Dxz_Dxz_S_a = I_ERI_Gxy2z_Px_Dxz_S_a+ABZ*I_ERI_Fxyz_Px_Dxz_S_a;
  Double I_ERI_Fx2z_Dxz_Dxz_S_a = I_ERI_Gx3z_Px_Dxz_S_a+ABZ*I_ERI_Fx2z_Px_Dxz_S_a;
  Double I_ERI_F3y_Dxz_Dxz_S_a = I_ERI_G3yz_Px_Dxz_S_a+ABZ*I_ERI_F3y_Px_Dxz_S_a;
  Double I_ERI_F2yz_Dxz_Dxz_S_a = I_ERI_G2y2z_Px_Dxz_S_a+ABZ*I_ERI_F2yz_Px_Dxz_S_a;
  Double I_ERI_Fy2z_Dxz_Dxz_S_a = I_ERI_Gy3z_Px_Dxz_S_a+ABZ*I_ERI_Fy2z_Px_Dxz_S_a;
  Double I_ERI_F3z_Dxz_Dxz_S_a = I_ERI_G4z_Px_Dxz_S_a+ABZ*I_ERI_F3z_Px_Dxz_S_a;
  Double I_ERI_F3x_D2y_Dxz_S_a = I_ERI_G3xy_Py_Dxz_S_a+ABY*I_ERI_F3x_Py_Dxz_S_a;
  Double I_ERI_F2xy_D2y_Dxz_S_a = I_ERI_G2x2y_Py_Dxz_S_a+ABY*I_ERI_F2xy_Py_Dxz_S_a;
  Double I_ERI_F2xz_D2y_Dxz_S_a = I_ERI_G2xyz_Py_Dxz_S_a+ABY*I_ERI_F2xz_Py_Dxz_S_a;
  Double I_ERI_Fx2y_D2y_Dxz_S_a = I_ERI_Gx3y_Py_Dxz_S_a+ABY*I_ERI_Fx2y_Py_Dxz_S_a;
  Double I_ERI_Fxyz_D2y_Dxz_S_a = I_ERI_Gx2yz_Py_Dxz_S_a+ABY*I_ERI_Fxyz_Py_Dxz_S_a;
  Double I_ERI_Fx2z_D2y_Dxz_S_a = I_ERI_Gxy2z_Py_Dxz_S_a+ABY*I_ERI_Fx2z_Py_Dxz_S_a;
  Double I_ERI_F3y_D2y_Dxz_S_a = I_ERI_G4y_Py_Dxz_S_a+ABY*I_ERI_F3y_Py_Dxz_S_a;
  Double I_ERI_F2yz_D2y_Dxz_S_a = I_ERI_G3yz_Py_Dxz_S_a+ABY*I_ERI_F2yz_Py_Dxz_S_a;
  Double I_ERI_Fy2z_D2y_Dxz_S_a = I_ERI_G2y2z_Py_Dxz_S_a+ABY*I_ERI_Fy2z_Py_Dxz_S_a;
  Double I_ERI_F3z_D2y_Dxz_S_a = I_ERI_Gy3z_Py_Dxz_S_a+ABY*I_ERI_F3z_Py_Dxz_S_a;
  Double I_ERI_F3x_Dyz_Dxz_S_a = I_ERI_G3xz_Py_Dxz_S_a+ABZ*I_ERI_F3x_Py_Dxz_S_a;
  Double I_ERI_F2xy_Dyz_Dxz_S_a = I_ERI_G2xyz_Py_Dxz_S_a+ABZ*I_ERI_F2xy_Py_Dxz_S_a;
  Double I_ERI_F2xz_Dyz_Dxz_S_a = I_ERI_G2x2z_Py_Dxz_S_a+ABZ*I_ERI_F2xz_Py_Dxz_S_a;
  Double I_ERI_Fx2y_Dyz_Dxz_S_a = I_ERI_Gx2yz_Py_Dxz_S_a+ABZ*I_ERI_Fx2y_Py_Dxz_S_a;
  Double I_ERI_Fxyz_Dyz_Dxz_S_a = I_ERI_Gxy2z_Py_Dxz_S_a+ABZ*I_ERI_Fxyz_Py_Dxz_S_a;
  Double I_ERI_Fx2z_Dyz_Dxz_S_a = I_ERI_Gx3z_Py_Dxz_S_a+ABZ*I_ERI_Fx2z_Py_Dxz_S_a;
  Double I_ERI_F3y_Dyz_Dxz_S_a = I_ERI_G3yz_Py_Dxz_S_a+ABZ*I_ERI_F3y_Py_Dxz_S_a;
  Double I_ERI_F2yz_Dyz_Dxz_S_a = I_ERI_G2y2z_Py_Dxz_S_a+ABZ*I_ERI_F2yz_Py_Dxz_S_a;
  Double I_ERI_Fy2z_Dyz_Dxz_S_a = I_ERI_Gy3z_Py_Dxz_S_a+ABZ*I_ERI_Fy2z_Py_Dxz_S_a;
  Double I_ERI_F3z_Dyz_Dxz_S_a = I_ERI_G4z_Py_Dxz_S_a+ABZ*I_ERI_F3z_Py_Dxz_S_a;
  Double I_ERI_F3x_D2z_Dxz_S_a = I_ERI_G3xz_Pz_Dxz_S_a+ABZ*I_ERI_F3x_Pz_Dxz_S_a;
  Double I_ERI_F2xy_D2z_Dxz_S_a = I_ERI_G2xyz_Pz_Dxz_S_a+ABZ*I_ERI_F2xy_Pz_Dxz_S_a;
  Double I_ERI_F2xz_D2z_Dxz_S_a = I_ERI_G2x2z_Pz_Dxz_S_a+ABZ*I_ERI_F2xz_Pz_Dxz_S_a;
  Double I_ERI_Fx2y_D2z_Dxz_S_a = I_ERI_Gx2yz_Pz_Dxz_S_a+ABZ*I_ERI_Fx2y_Pz_Dxz_S_a;
  Double I_ERI_Fxyz_D2z_Dxz_S_a = I_ERI_Gxy2z_Pz_Dxz_S_a+ABZ*I_ERI_Fxyz_Pz_Dxz_S_a;
  Double I_ERI_Fx2z_D2z_Dxz_S_a = I_ERI_Gx3z_Pz_Dxz_S_a+ABZ*I_ERI_Fx2z_Pz_Dxz_S_a;
  Double I_ERI_F3y_D2z_Dxz_S_a = I_ERI_G3yz_Pz_Dxz_S_a+ABZ*I_ERI_F3y_Pz_Dxz_S_a;
  Double I_ERI_F2yz_D2z_Dxz_S_a = I_ERI_G2y2z_Pz_Dxz_S_a+ABZ*I_ERI_F2yz_Pz_Dxz_S_a;
  Double I_ERI_Fy2z_D2z_Dxz_S_a = I_ERI_Gy3z_Pz_Dxz_S_a+ABZ*I_ERI_Fy2z_Pz_Dxz_S_a;
  Double I_ERI_F3z_D2z_Dxz_S_a = I_ERI_G4z_Pz_Dxz_S_a+ABZ*I_ERI_F3z_Pz_Dxz_S_a;
  Double I_ERI_F3x_D2x_D2y_S_a = I_ERI_G4x_Px_D2y_S_a+ABX*I_ERI_F3x_Px_D2y_S_a;
  Double I_ERI_F2xy_D2x_D2y_S_a = I_ERI_G3xy_Px_D2y_S_a+ABX*I_ERI_F2xy_Px_D2y_S_a;
  Double I_ERI_F2xz_D2x_D2y_S_a = I_ERI_G3xz_Px_D2y_S_a+ABX*I_ERI_F2xz_Px_D2y_S_a;
  Double I_ERI_Fx2y_D2x_D2y_S_a = I_ERI_G2x2y_Px_D2y_S_a+ABX*I_ERI_Fx2y_Px_D2y_S_a;
  Double I_ERI_Fxyz_D2x_D2y_S_a = I_ERI_G2xyz_Px_D2y_S_a+ABX*I_ERI_Fxyz_Px_D2y_S_a;
  Double I_ERI_Fx2z_D2x_D2y_S_a = I_ERI_G2x2z_Px_D2y_S_a+ABX*I_ERI_Fx2z_Px_D2y_S_a;
  Double I_ERI_F3y_D2x_D2y_S_a = I_ERI_Gx3y_Px_D2y_S_a+ABX*I_ERI_F3y_Px_D2y_S_a;
  Double I_ERI_F2yz_D2x_D2y_S_a = I_ERI_Gx2yz_Px_D2y_S_a+ABX*I_ERI_F2yz_Px_D2y_S_a;
  Double I_ERI_Fy2z_D2x_D2y_S_a = I_ERI_Gxy2z_Px_D2y_S_a+ABX*I_ERI_Fy2z_Px_D2y_S_a;
  Double I_ERI_F3z_D2x_D2y_S_a = I_ERI_Gx3z_Px_D2y_S_a+ABX*I_ERI_F3z_Px_D2y_S_a;
  Double I_ERI_F3x_Dxy_D2y_S_a = I_ERI_G3xy_Px_D2y_S_a+ABY*I_ERI_F3x_Px_D2y_S_a;
  Double I_ERI_F2xy_Dxy_D2y_S_a = I_ERI_G2x2y_Px_D2y_S_a+ABY*I_ERI_F2xy_Px_D2y_S_a;
  Double I_ERI_F2xz_Dxy_D2y_S_a = I_ERI_G2xyz_Px_D2y_S_a+ABY*I_ERI_F2xz_Px_D2y_S_a;
  Double I_ERI_Fx2y_Dxy_D2y_S_a = I_ERI_Gx3y_Px_D2y_S_a+ABY*I_ERI_Fx2y_Px_D2y_S_a;
  Double I_ERI_Fxyz_Dxy_D2y_S_a = I_ERI_Gx2yz_Px_D2y_S_a+ABY*I_ERI_Fxyz_Px_D2y_S_a;
  Double I_ERI_Fx2z_Dxy_D2y_S_a = I_ERI_Gxy2z_Px_D2y_S_a+ABY*I_ERI_Fx2z_Px_D2y_S_a;
  Double I_ERI_F3y_Dxy_D2y_S_a = I_ERI_G4y_Px_D2y_S_a+ABY*I_ERI_F3y_Px_D2y_S_a;
  Double I_ERI_F2yz_Dxy_D2y_S_a = I_ERI_G3yz_Px_D2y_S_a+ABY*I_ERI_F2yz_Px_D2y_S_a;
  Double I_ERI_Fy2z_Dxy_D2y_S_a = I_ERI_G2y2z_Px_D2y_S_a+ABY*I_ERI_Fy2z_Px_D2y_S_a;
  Double I_ERI_F3z_Dxy_D2y_S_a = I_ERI_Gy3z_Px_D2y_S_a+ABY*I_ERI_F3z_Px_D2y_S_a;
  Double I_ERI_F3x_Dxz_D2y_S_a = I_ERI_G3xz_Px_D2y_S_a+ABZ*I_ERI_F3x_Px_D2y_S_a;
  Double I_ERI_F2xy_Dxz_D2y_S_a = I_ERI_G2xyz_Px_D2y_S_a+ABZ*I_ERI_F2xy_Px_D2y_S_a;
  Double I_ERI_F2xz_Dxz_D2y_S_a = I_ERI_G2x2z_Px_D2y_S_a+ABZ*I_ERI_F2xz_Px_D2y_S_a;
  Double I_ERI_Fx2y_Dxz_D2y_S_a = I_ERI_Gx2yz_Px_D2y_S_a+ABZ*I_ERI_Fx2y_Px_D2y_S_a;
  Double I_ERI_Fxyz_Dxz_D2y_S_a = I_ERI_Gxy2z_Px_D2y_S_a+ABZ*I_ERI_Fxyz_Px_D2y_S_a;
  Double I_ERI_Fx2z_Dxz_D2y_S_a = I_ERI_Gx3z_Px_D2y_S_a+ABZ*I_ERI_Fx2z_Px_D2y_S_a;
  Double I_ERI_F3y_Dxz_D2y_S_a = I_ERI_G3yz_Px_D2y_S_a+ABZ*I_ERI_F3y_Px_D2y_S_a;
  Double I_ERI_F2yz_Dxz_D2y_S_a = I_ERI_G2y2z_Px_D2y_S_a+ABZ*I_ERI_F2yz_Px_D2y_S_a;
  Double I_ERI_Fy2z_Dxz_D2y_S_a = I_ERI_Gy3z_Px_D2y_S_a+ABZ*I_ERI_Fy2z_Px_D2y_S_a;
  Double I_ERI_F3z_Dxz_D2y_S_a = I_ERI_G4z_Px_D2y_S_a+ABZ*I_ERI_F3z_Px_D2y_S_a;
  Double I_ERI_F3x_D2y_D2y_S_a = I_ERI_G3xy_Py_D2y_S_a+ABY*I_ERI_F3x_Py_D2y_S_a;
  Double I_ERI_F2xy_D2y_D2y_S_a = I_ERI_G2x2y_Py_D2y_S_a+ABY*I_ERI_F2xy_Py_D2y_S_a;
  Double I_ERI_F2xz_D2y_D2y_S_a = I_ERI_G2xyz_Py_D2y_S_a+ABY*I_ERI_F2xz_Py_D2y_S_a;
  Double I_ERI_Fx2y_D2y_D2y_S_a = I_ERI_Gx3y_Py_D2y_S_a+ABY*I_ERI_Fx2y_Py_D2y_S_a;
  Double I_ERI_Fxyz_D2y_D2y_S_a = I_ERI_Gx2yz_Py_D2y_S_a+ABY*I_ERI_Fxyz_Py_D2y_S_a;
  Double I_ERI_Fx2z_D2y_D2y_S_a = I_ERI_Gxy2z_Py_D2y_S_a+ABY*I_ERI_Fx2z_Py_D2y_S_a;
  Double I_ERI_F3y_D2y_D2y_S_a = I_ERI_G4y_Py_D2y_S_a+ABY*I_ERI_F3y_Py_D2y_S_a;
  Double I_ERI_F2yz_D2y_D2y_S_a = I_ERI_G3yz_Py_D2y_S_a+ABY*I_ERI_F2yz_Py_D2y_S_a;
  Double I_ERI_Fy2z_D2y_D2y_S_a = I_ERI_G2y2z_Py_D2y_S_a+ABY*I_ERI_Fy2z_Py_D2y_S_a;
  Double I_ERI_F3z_D2y_D2y_S_a = I_ERI_Gy3z_Py_D2y_S_a+ABY*I_ERI_F3z_Py_D2y_S_a;
  Double I_ERI_F3x_Dyz_D2y_S_a = I_ERI_G3xz_Py_D2y_S_a+ABZ*I_ERI_F3x_Py_D2y_S_a;
  Double I_ERI_F2xy_Dyz_D2y_S_a = I_ERI_G2xyz_Py_D2y_S_a+ABZ*I_ERI_F2xy_Py_D2y_S_a;
  Double I_ERI_F2xz_Dyz_D2y_S_a = I_ERI_G2x2z_Py_D2y_S_a+ABZ*I_ERI_F2xz_Py_D2y_S_a;
  Double I_ERI_Fx2y_Dyz_D2y_S_a = I_ERI_Gx2yz_Py_D2y_S_a+ABZ*I_ERI_Fx2y_Py_D2y_S_a;
  Double I_ERI_Fxyz_Dyz_D2y_S_a = I_ERI_Gxy2z_Py_D2y_S_a+ABZ*I_ERI_Fxyz_Py_D2y_S_a;
  Double I_ERI_Fx2z_Dyz_D2y_S_a = I_ERI_Gx3z_Py_D2y_S_a+ABZ*I_ERI_Fx2z_Py_D2y_S_a;
  Double I_ERI_F3y_Dyz_D2y_S_a = I_ERI_G3yz_Py_D2y_S_a+ABZ*I_ERI_F3y_Py_D2y_S_a;
  Double I_ERI_F2yz_Dyz_D2y_S_a = I_ERI_G2y2z_Py_D2y_S_a+ABZ*I_ERI_F2yz_Py_D2y_S_a;
  Double I_ERI_Fy2z_Dyz_D2y_S_a = I_ERI_Gy3z_Py_D2y_S_a+ABZ*I_ERI_Fy2z_Py_D2y_S_a;
  Double I_ERI_F3z_Dyz_D2y_S_a = I_ERI_G4z_Py_D2y_S_a+ABZ*I_ERI_F3z_Py_D2y_S_a;
  Double I_ERI_F3x_D2z_D2y_S_a = I_ERI_G3xz_Pz_D2y_S_a+ABZ*I_ERI_F3x_Pz_D2y_S_a;
  Double I_ERI_F2xy_D2z_D2y_S_a = I_ERI_G2xyz_Pz_D2y_S_a+ABZ*I_ERI_F2xy_Pz_D2y_S_a;
  Double I_ERI_F2xz_D2z_D2y_S_a = I_ERI_G2x2z_Pz_D2y_S_a+ABZ*I_ERI_F2xz_Pz_D2y_S_a;
  Double I_ERI_Fx2y_D2z_D2y_S_a = I_ERI_Gx2yz_Pz_D2y_S_a+ABZ*I_ERI_Fx2y_Pz_D2y_S_a;
  Double I_ERI_Fxyz_D2z_D2y_S_a = I_ERI_Gxy2z_Pz_D2y_S_a+ABZ*I_ERI_Fxyz_Pz_D2y_S_a;
  Double I_ERI_Fx2z_D2z_D2y_S_a = I_ERI_Gx3z_Pz_D2y_S_a+ABZ*I_ERI_Fx2z_Pz_D2y_S_a;
  Double I_ERI_F3y_D2z_D2y_S_a = I_ERI_G3yz_Pz_D2y_S_a+ABZ*I_ERI_F3y_Pz_D2y_S_a;
  Double I_ERI_F2yz_D2z_D2y_S_a = I_ERI_G2y2z_Pz_D2y_S_a+ABZ*I_ERI_F2yz_Pz_D2y_S_a;
  Double I_ERI_Fy2z_D2z_D2y_S_a = I_ERI_Gy3z_Pz_D2y_S_a+ABZ*I_ERI_Fy2z_Pz_D2y_S_a;
  Double I_ERI_F3z_D2z_D2y_S_a = I_ERI_G4z_Pz_D2y_S_a+ABZ*I_ERI_F3z_Pz_D2y_S_a;
  Double I_ERI_F3x_D2x_Dyz_S_a = I_ERI_G4x_Px_Dyz_S_a+ABX*I_ERI_F3x_Px_Dyz_S_a;
  Double I_ERI_F2xy_D2x_Dyz_S_a = I_ERI_G3xy_Px_Dyz_S_a+ABX*I_ERI_F2xy_Px_Dyz_S_a;
  Double I_ERI_F2xz_D2x_Dyz_S_a = I_ERI_G3xz_Px_Dyz_S_a+ABX*I_ERI_F2xz_Px_Dyz_S_a;
  Double I_ERI_Fx2y_D2x_Dyz_S_a = I_ERI_G2x2y_Px_Dyz_S_a+ABX*I_ERI_Fx2y_Px_Dyz_S_a;
  Double I_ERI_Fxyz_D2x_Dyz_S_a = I_ERI_G2xyz_Px_Dyz_S_a+ABX*I_ERI_Fxyz_Px_Dyz_S_a;
  Double I_ERI_Fx2z_D2x_Dyz_S_a = I_ERI_G2x2z_Px_Dyz_S_a+ABX*I_ERI_Fx2z_Px_Dyz_S_a;
  Double I_ERI_F3y_D2x_Dyz_S_a = I_ERI_Gx3y_Px_Dyz_S_a+ABX*I_ERI_F3y_Px_Dyz_S_a;
  Double I_ERI_F2yz_D2x_Dyz_S_a = I_ERI_Gx2yz_Px_Dyz_S_a+ABX*I_ERI_F2yz_Px_Dyz_S_a;
  Double I_ERI_Fy2z_D2x_Dyz_S_a = I_ERI_Gxy2z_Px_Dyz_S_a+ABX*I_ERI_Fy2z_Px_Dyz_S_a;
  Double I_ERI_F3z_D2x_Dyz_S_a = I_ERI_Gx3z_Px_Dyz_S_a+ABX*I_ERI_F3z_Px_Dyz_S_a;
  Double I_ERI_F3x_Dxy_Dyz_S_a = I_ERI_G3xy_Px_Dyz_S_a+ABY*I_ERI_F3x_Px_Dyz_S_a;
  Double I_ERI_F2xy_Dxy_Dyz_S_a = I_ERI_G2x2y_Px_Dyz_S_a+ABY*I_ERI_F2xy_Px_Dyz_S_a;
  Double I_ERI_F2xz_Dxy_Dyz_S_a = I_ERI_G2xyz_Px_Dyz_S_a+ABY*I_ERI_F2xz_Px_Dyz_S_a;
  Double I_ERI_Fx2y_Dxy_Dyz_S_a = I_ERI_Gx3y_Px_Dyz_S_a+ABY*I_ERI_Fx2y_Px_Dyz_S_a;
  Double I_ERI_Fxyz_Dxy_Dyz_S_a = I_ERI_Gx2yz_Px_Dyz_S_a+ABY*I_ERI_Fxyz_Px_Dyz_S_a;
  Double I_ERI_Fx2z_Dxy_Dyz_S_a = I_ERI_Gxy2z_Px_Dyz_S_a+ABY*I_ERI_Fx2z_Px_Dyz_S_a;
  Double I_ERI_F3y_Dxy_Dyz_S_a = I_ERI_G4y_Px_Dyz_S_a+ABY*I_ERI_F3y_Px_Dyz_S_a;
  Double I_ERI_F2yz_Dxy_Dyz_S_a = I_ERI_G3yz_Px_Dyz_S_a+ABY*I_ERI_F2yz_Px_Dyz_S_a;
  Double I_ERI_Fy2z_Dxy_Dyz_S_a = I_ERI_G2y2z_Px_Dyz_S_a+ABY*I_ERI_Fy2z_Px_Dyz_S_a;
  Double I_ERI_F3z_Dxy_Dyz_S_a = I_ERI_Gy3z_Px_Dyz_S_a+ABY*I_ERI_F3z_Px_Dyz_S_a;
  Double I_ERI_F3x_Dxz_Dyz_S_a = I_ERI_G3xz_Px_Dyz_S_a+ABZ*I_ERI_F3x_Px_Dyz_S_a;
  Double I_ERI_F2xy_Dxz_Dyz_S_a = I_ERI_G2xyz_Px_Dyz_S_a+ABZ*I_ERI_F2xy_Px_Dyz_S_a;
  Double I_ERI_F2xz_Dxz_Dyz_S_a = I_ERI_G2x2z_Px_Dyz_S_a+ABZ*I_ERI_F2xz_Px_Dyz_S_a;
  Double I_ERI_Fx2y_Dxz_Dyz_S_a = I_ERI_Gx2yz_Px_Dyz_S_a+ABZ*I_ERI_Fx2y_Px_Dyz_S_a;
  Double I_ERI_Fxyz_Dxz_Dyz_S_a = I_ERI_Gxy2z_Px_Dyz_S_a+ABZ*I_ERI_Fxyz_Px_Dyz_S_a;
  Double I_ERI_Fx2z_Dxz_Dyz_S_a = I_ERI_Gx3z_Px_Dyz_S_a+ABZ*I_ERI_Fx2z_Px_Dyz_S_a;
  Double I_ERI_F3y_Dxz_Dyz_S_a = I_ERI_G3yz_Px_Dyz_S_a+ABZ*I_ERI_F3y_Px_Dyz_S_a;
  Double I_ERI_F2yz_Dxz_Dyz_S_a = I_ERI_G2y2z_Px_Dyz_S_a+ABZ*I_ERI_F2yz_Px_Dyz_S_a;
  Double I_ERI_Fy2z_Dxz_Dyz_S_a = I_ERI_Gy3z_Px_Dyz_S_a+ABZ*I_ERI_Fy2z_Px_Dyz_S_a;
  Double I_ERI_F3z_Dxz_Dyz_S_a = I_ERI_G4z_Px_Dyz_S_a+ABZ*I_ERI_F3z_Px_Dyz_S_a;
  Double I_ERI_F3x_D2y_Dyz_S_a = I_ERI_G3xy_Py_Dyz_S_a+ABY*I_ERI_F3x_Py_Dyz_S_a;
  Double I_ERI_F2xy_D2y_Dyz_S_a = I_ERI_G2x2y_Py_Dyz_S_a+ABY*I_ERI_F2xy_Py_Dyz_S_a;
  Double I_ERI_F2xz_D2y_Dyz_S_a = I_ERI_G2xyz_Py_Dyz_S_a+ABY*I_ERI_F2xz_Py_Dyz_S_a;
  Double I_ERI_Fx2y_D2y_Dyz_S_a = I_ERI_Gx3y_Py_Dyz_S_a+ABY*I_ERI_Fx2y_Py_Dyz_S_a;
  Double I_ERI_Fxyz_D2y_Dyz_S_a = I_ERI_Gx2yz_Py_Dyz_S_a+ABY*I_ERI_Fxyz_Py_Dyz_S_a;
  Double I_ERI_Fx2z_D2y_Dyz_S_a = I_ERI_Gxy2z_Py_Dyz_S_a+ABY*I_ERI_Fx2z_Py_Dyz_S_a;
  Double I_ERI_F3y_D2y_Dyz_S_a = I_ERI_G4y_Py_Dyz_S_a+ABY*I_ERI_F3y_Py_Dyz_S_a;
  Double I_ERI_F2yz_D2y_Dyz_S_a = I_ERI_G3yz_Py_Dyz_S_a+ABY*I_ERI_F2yz_Py_Dyz_S_a;
  Double I_ERI_Fy2z_D2y_Dyz_S_a = I_ERI_G2y2z_Py_Dyz_S_a+ABY*I_ERI_Fy2z_Py_Dyz_S_a;
  Double I_ERI_F3z_D2y_Dyz_S_a = I_ERI_Gy3z_Py_Dyz_S_a+ABY*I_ERI_F3z_Py_Dyz_S_a;
  Double I_ERI_F3x_Dyz_Dyz_S_a = I_ERI_G3xz_Py_Dyz_S_a+ABZ*I_ERI_F3x_Py_Dyz_S_a;
  Double I_ERI_F2xy_Dyz_Dyz_S_a = I_ERI_G2xyz_Py_Dyz_S_a+ABZ*I_ERI_F2xy_Py_Dyz_S_a;
  Double I_ERI_F2xz_Dyz_Dyz_S_a = I_ERI_G2x2z_Py_Dyz_S_a+ABZ*I_ERI_F2xz_Py_Dyz_S_a;
  Double I_ERI_Fx2y_Dyz_Dyz_S_a = I_ERI_Gx2yz_Py_Dyz_S_a+ABZ*I_ERI_Fx2y_Py_Dyz_S_a;
  Double I_ERI_Fxyz_Dyz_Dyz_S_a = I_ERI_Gxy2z_Py_Dyz_S_a+ABZ*I_ERI_Fxyz_Py_Dyz_S_a;
  Double I_ERI_Fx2z_Dyz_Dyz_S_a = I_ERI_Gx3z_Py_Dyz_S_a+ABZ*I_ERI_Fx2z_Py_Dyz_S_a;
  Double I_ERI_F3y_Dyz_Dyz_S_a = I_ERI_G3yz_Py_Dyz_S_a+ABZ*I_ERI_F3y_Py_Dyz_S_a;
  Double I_ERI_F2yz_Dyz_Dyz_S_a = I_ERI_G2y2z_Py_Dyz_S_a+ABZ*I_ERI_F2yz_Py_Dyz_S_a;
  Double I_ERI_Fy2z_Dyz_Dyz_S_a = I_ERI_Gy3z_Py_Dyz_S_a+ABZ*I_ERI_Fy2z_Py_Dyz_S_a;
  Double I_ERI_F3z_Dyz_Dyz_S_a = I_ERI_G4z_Py_Dyz_S_a+ABZ*I_ERI_F3z_Py_Dyz_S_a;
  Double I_ERI_F3x_D2z_Dyz_S_a = I_ERI_G3xz_Pz_Dyz_S_a+ABZ*I_ERI_F3x_Pz_Dyz_S_a;
  Double I_ERI_F2xy_D2z_Dyz_S_a = I_ERI_G2xyz_Pz_Dyz_S_a+ABZ*I_ERI_F2xy_Pz_Dyz_S_a;
  Double I_ERI_F2xz_D2z_Dyz_S_a = I_ERI_G2x2z_Pz_Dyz_S_a+ABZ*I_ERI_F2xz_Pz_Dyz_S_a;
  Double I_ERI_Fx2y_D2z_Dyz_S_a = I_ERI_Gx2yz_Pz_Dyz_S_a+ABZ*I_ERI_Fx2y_Pz_Dyz_S_a;
  Double I_ERI_Fxyz_D2z_Dyz_S_a = I_ERI_Gxy2z_Pz_Dyz_S_a+ABZ*I_ERI_Fxyz_Pz_Dyz_S_a;
  Double I_ERI_Fx2z_D2z_Dyz_S_a = I_ERI_Gx3z_Pz_Dyz_S_a+ABZ*I_ERI_Fx2z_Pz_Dyz_S_a;
  Double I_ERI_F3y_D2z_Dyz_S_a = I_ERI_G3yz_Pz_Dyz_S_a+ABZ*I_ERI_F3y_Pz_Dyz_S_a;
  Double I_ERI_F2yz_D2z_Dyz_S_a = I_ERI_G2y2z_Pz_Dyz_S_a+ABZ*I_ERI_F2yz_Pz_Dyz_S_a;
  Double I_ERI_Fy2z_D2z_Dyz_S_a = I_ERI_Gy3z_Pz_Dyz_S_a+ABZ*I_ERI_Fy2z_Pz_Dyz_S_a;
  Double I_ERI_F3z_D2z_Dyz_S_a = I_ERI_G4z_Pz_Dyz_S_a+ABZ*I_ERI_F3z_Pz_Dyz_S_a;
  Double I_ERI_F3x_D2x_D2z_S_a = I_ERI_G4x_Px_D2z_S_a+ABX*I_ERI_F3x_Px_D2z_S_a;
  Double I_ERI_F2xy_D2x_D2z_S_a = I_ERI_G3xy_Px_D2z_S_a+ABX*I_ERI_F2xy_Px_D2z_S_a;
  Double I_ERI_F2xz_D2x_D2z_S_a = I_ERI_G3xz_Px_D2z_S_a+ABX*I_ERI_F2xz_Px_D2z_S_a;
  Double I_ERI_Fx2y_D2x_D2z_S_a = I_ERI_G2x2y_Px_D2z_S_a+ABX*I_ERI_Fx2y_Px_D2z_S_a;
  Double I_ERI_Fxyz_D2x_D2z_S_a = I_ERI_G2xyz_Px_D2z_S_a+ABX*I_ERI_Fxyz_Px_D2z_S_a;
  Double I_ERI_Fx2z_D2x_D2z_S_a = I_ERI_G2x2z_Px_D2z_S_a+ABX*I_ERI_Fx2z_Px_D2z_S_a;
  Double I_ERI_F3y_D2x_D2z_S_a = I_ERI_Gx3y_Px_D2z_S_a+ABX*I_ERI_F3y_Px_D2z_S_a;
  Double I_ERI_F2yz_D2x_D2z_S_a = I_ERI_Gx2yz_Px_D2z_S_a+ABX*I_ERI_F2yz_Px_D2z_S_a;
  Double I_ERI_Fy2z_D2x_D2z_S_a = I_ERI_Gxy2z_Px_D2z_S_a+ABX*I_ERI_Fy2z_Px_D2z_S_a;
  Double I_ERI_F3z_D2x_D2z_S_a = I_ERI_Gx3z_Px_D2z_S_a+ABX*I_ERI_F3z_Px_D2z_S_a;
  Double I_ERI_F3x_Dxy_D2z_S_a = I_ERI_G3xy_Px_D2z_S_a+ABY*I_ERI_F3x_Px_D2z_S_a;
  Double I_ERI_F2xy_Dxy_D2z_S_a = I_ERI_G2x2y_Px_D2z_S_a+ABY*I_ERI_F2xy_Px_D2z_S_a;
  Double I_ERI_F2xz_Dxy_D2z_S_a = I_ERI_G2xyz_Px_D2z_S_a+ABY*I_ERI_F2xz_Px_D2z_S_a;
  Double I_ERI_Fx2y_Dxy_D2z_S_a = I_ERI_Gx3y_Px_D2z_S_a+ABY*I_ERI_Fx2y_Px_D2z_S_a;
  Double I_ERI_Fxyz_Dxy_D2z_S_a = I_ERI_Gx2yz_Px_D2z_S_a+ABY*I_ERI_Fxyz_Px_D2z_S_a;
  Double I_ERI_Fx2z_Dxy_D2z_S_a = I_ERI_Gxy2z_Px_D2z_S_a+ABY*I_ERI_Fx2z_Px_D2z_S_a;
  Double I_ERI_F3y_Dxy_D2z_S_a = I_ERI_G4y_Px_D2z_S_a+ABY*I_ERI_F3y_Px_D2z_S_a;
  Double I_ERI_F2yz_Dxy_D2z_S_a = I_ERI_G3yz_Px_D2z_S_a+ABY*I_ERI_F2yz_Px_D2z_S_a;
  Double I_ERI_Fy2z_Dxy_D2z_S_a = I_ERI_G2y2z_Px_D2z_S_a+ABY*I_ERI_Fy2z_Px_D2z_S_a;
  Double I_ERI_F3z_Dxy_D2z_S_a = I_ERI_Gy3z_Px_D2z_S_a+ABY*I_ERI_F3z_Px_D2z_S_a;
  Double I_ERI_F3x_Dxz_D2z_S_a = I_ERI_G3xz_Px_D2z_S_a+ABZ*I_ERI_F3x_Px_D2z_S_a;
  Double I_ERI_F2xy_Dxz_D2z_S_a = I_ERI_G2xyz_Px_D2z_S_a+ABZ*I_ERI_F2xy_Px_D2z_S_a;
  Double I_ERI_F2xz_Dxz_D2z_S_a = I_ERI_G2x2z_Px_D2z_S_a+ABZ*I_ERI_F2xz_Px_D2z_S_a;
  Double I_ERI_Fx2y_Dxz_D2z_S_a = I_ERI_Gx2yz_Px_D2z_S_a+ABZ*I_ERI_Fx2y_Px_D2z_S_a;
  Double I_ERI_Fxyz_Dxz_D2z_S_a = I_ERI_Gxy2z_Px_D2z_S_a+ABZ*I_ERI_Fxyz_Px_D2z_S_a;
  Double I_ERI_Fx2z_Dxz_D2z_S_a = I_ERI_Gx3z_Px_D2z_S_a+ABZ*I_ERI_Fx2z_Px_D2z_S_a;
  Double I_ERI_F3y_Dxz_D2z_S_a = I_ERI_G3yz_Px_D2z_S_a+ABZ*I_ERI_F3y_Px_D2z_S_a;
  Double I_ERI_F2yz_Dxz_D2z_S_a = I_ERI_G2y2z_Px_D2z_S_a+ABZ*I_ERI_F2yz_Px_D2z_S_a;
  Double I_ERI_Fy2z_Dxz_D2z_S_a = I_ERI_Gy3z_Px_D2z_S_a+ABZ*I_ERI_Fy2z_Px_D2z_S_a;
  Double I_ERI_F3z_Dxz_D2z_S_a = I_ERI_G4z_Px_D2z_S_a+ABZ*I_ERI_F3z_Px_D2z_S_a;
  Double I_ERI_F3x_D2y_D2z_S_a = I_ERI_G3xy_Py_D2z_S_a+ABY*I_ERI_F3x_Py_D2z_S_a;
  Double I_ERI_F2xy_D2y_D2z_S_a = I_ERI_G2x2y_Py_D2z_S_a+ABY*I_ERI_F2xy_Py_D2z_S_a;
  Double I_ERI_F2xz_D2y_D2z_S_a = I_ERI_G2xyz_Py_D2z_S_a+ABY*I_ERI_F2xz_Py_D2z_S_a;
  Double I_ERI_Fx2y_D2y_D2z_S_a = I_ERI_Gx3y_Py_D2z_S_a+ABY*I_ERI_Fx2y_Py_D2z_S_a;
  Double I_ERI_Fxyz_D2y_D2z_S_a = I_ERI_Gx2yz_Py_D2z_S_a+ABY*I_ERI_Fxyz_Py_D2z_S_a;
  Double I_ERI_Fx2z_D2y_D2z_S_a = I_ERI_Gxy2z_Py_D2z_S_a+ABY*I_ERI_Fx2z_Py_D2z_S_a;
  Double I_ERI_F3y_D2y_D2z_S_a = I_ERI_G4y_Py_D2z_S_a+ABY*I_ERI_F3y_Py_D2z_S_a;
  Double I_ERI_F2yz_D2y_D2z_S_a = I_ERI_G3yz_Py_D2z_S_a+ABY*I_ERI_F2yz_Py_D2z_S_a;
  Double I_ERI_Fy2z_D2y_D2z_S_a = I_ERI_G2y2z_Py_D2z_S_a+ABY*I_ERI_Fy2z_Py_D2z_S_a;
  Double I_ERI_F3z_D2y_D2z_S_a = I_ERI_Gy3z_Py_D2z_S_a+ABY*I_ERI_F3z_Py_D2z_S_a;
  Double I_ERI_F3x_Dyz_D2z_S_a = I_ERI_G3xz_Py_D2z_S_a+ABZ*I_ERI_F3x_Py_D2z_S_a;
  Double I_ERI_F2xy_Dyz_D2z_S_a = I_ERI_G2xyz_Py_D2z_S_a+ABZ*I_ERI_F2xy_Py_D2z_S_a;
  Double I_ERI_F2xz_Dyz_D2z_S_a = I_ERI_G2x2z_Py_D2z_S_a+ABZ*I_ERI_F2xz_Py_D2z_S_a;
  Double I_ERI_Fx2y_Dyz_D2z_S_a = I_ERI_Gx2yz_Py_D2z_S_a+ABZ*I_ERI_Fx2y_Py_D2z_S_a;
  Double I_ERI_Fxyz_Dyz_D2z_S_a = I_ERI_Gxy2z_Py_D2z_S_a+ABZ*I_ERI_Fxyz_Py_D2z_S_a;
  Double I_ERI_Fx2z_Dyz_D2z_S_a = I_ERI_Gx3z_Py_D2z_S_a+ABZ*I_ERI_Fx2z_Py_D2z_S_a;
  Double I_ERI_F3y_Dyz_D2z_S_a = I_ERI_G3yz_Py_D2z_S_a+ABZ*I_ERI_F3y_Py_D2z_S_a;
  Double I_ERI_F2yz_Dyz_D2z_S_a = I_ERI_G2y2z_Py_D2z_S_a+ABZ*I_ERI_F2yz_Py_D2z_S_a;
  Double I_ERI_Fy2z_Dyz_D2z_S_a = I_ERI_Gy3z_Py_D2z_S_a+ABZ*I_ERI_Fy2z_Py_D2z_S_a;
  Double I_ERI_F3z_Dyz_D2z_S_a = I_ERI_G4z_Py_D2z_S_a+ABZ*I_ERI_F3z_Py_D2z_S_a;
  Double I_ERI_F3x_D2z_D2z_S_a = I_ERI_G3xz_Pz_D2z_S_a+ABZ*I_ERI_F3x_Pz_D2z_S_a;
  Double I_ERI_F2xy_D2z_D2z_S_a = I_ERI_G2xyz_Pz_D2z_S_a+ABZ*I_ERI_F2xy_Pz_D2z_S_a;
  Double I_ERI_F2xz_D2z_D2z_S_a = I_ERI_G2x2z_Pz_D2z_S_a+ABZ*I_ERI_F2xz_Pz_D2z_S_a;
  Double I_ERI_Fx2y_D2z_D2z_S_a = I_ERI_Gx2yz_Pz_D2z_S_a+ABZ*I_ERI_Fx2y_Pz_D2z_S_a;
  Double I_ERI_Fxyz_D2z_D2z_S_a = I_ERI_Gxy2z_Pz_D2z_S_a+ABZ*I_ERI_Fxyz_Pz_D2z_S_a;
  Double I_ERI_Fx2z_D2z_D2z_S_a = I_ERI_Gx3z_Pz_D2z_S_a+ABZ*I_ERI_Fx2z_Pz_D2z_S_a;
  Double I_ERI_F3y_D2z_D2z_S_a = I_ERI_G3yz_Pz_D2z_S_a+ABZ*I_ERI_F3y_Pz_D2z_S_a;
  Double I_ERI_F2yz_D2z_D2z_S_a = I_ERI_G2y2z_Pz_D2z_S_a+ABZ*I_ERI_F2yz_Pz_D2z_S_a;
  Double I_ERI_Fy2z_D2z_D2z_S_a = I_ERI_Gy3z_Pz_D2z_S_a+ABZ*I_ERI_Fy2z_Pz_D2z_S_a;
  Double I_ERI_F3z_D2z_D2z_S_a = I_ERI_G4z_Pz_D2z_S_a+ABZ*I_ERI_F3z_Pz_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_b
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S_b = I_ERI_F3x_S_D2x_S_b+ABX*I_ERI_D2x_S_D2x_S_b;
  Double I_ERI_Dxy_Px_D2x_S_b = I_ERI_F2xy_S_D2x_S_b+ABX*I_ERI_Dxy_S_D2x_S_b;
  Double I_ERI_Dxz_Px_D2x_S_b = I_ERI_F2xz_S_D2x_S_b+ABX*I_ERI_Dxz_S_D2x_S_b;
  Double I_ERI_D2y_Px_D2x_S_b = I_ERI_Fx2y_S_D2x_S_b+ABX*I_ERI_D2y_S_D2x_S_b;
  Double I_ERI_Dyz_Px_D2x_S_b = I_ERI_Fxyz_S_D2x_S_b+ABX*I_ERI_Dyz_S_D2x_S_b;
  Double I_ERI_D2z_Px_D2x_S_b = I_ERI_Fx2z_S_D2x_S_b+ABX*I_ERI_D2z_S_D2x_S_b;
  Double I_ERI_D2x_Py_D2x_S_b = I_ERI_F2xy_S_D2x_S_b+ABY*I_ERI_D2x_S_D2x_S_b;
  Double I_ERI_Dxy_Py_D2x_S_b = I_ERI_Fx2y_S_D2x_S_b+ABY*I_ERI_Dxy_S_D2x_S_b;
  Double I_ERI_Dxz_Py_D2x_S_b = I_ERI_Fxyz_S_D2x_S_b+ABY*I_ERI_Dxz_S_D2x_S_b;
  Double I_ERI_D2y_Py_D2x_S_b = I_ERI_F3y_S_D2x_S_b+ABY*I_ERI_D2y_S_D2x_S_b;
  Double I_ERI_Dyz_Py_D2x_S_b = I_ERI_F2yz_S_D2x_S_b+ABY*I_ERI_Dyz_S_D2x_S_b;
  Double I_ERI_D2z_Py_D2x_S_b = I_ERI_Fy2z_S_D2x_S_b+ABY*I_ERI_D2z_S_D2x_S_b;
  Double I_ERI_D2x_Pz_D2x_S_b = I_ERI_F2xz_S_D2x_S_b+ABZ*I_ERI_D2x_S_D2x_S_b;
  Double I_ERI_Dxy_Pz_D2x_S_b = I_ERI_Fxyz_S_D2x_S_b+ABZ*I_ERI_Dxy_S_D2x_S_b;
  Double I_ERI_Dxz_Pz_D2x_S_b = I_ERI_Fx2z_S_D2x_S_b+ABZ*I_ERI_Dxz_S_D2x_S_b;
  Double I_ERI_D2y_Pz_D2x_S_b = I_ERI_F2yz_S_D2x_S_b+ABZ*I_ERI_D2y_S_D2x_S_b;
  Double I_ERI_Dyz_Pz_D2x_S_b = I_ERI_Fy2z_S_D2x_S_b+ABZ*I_ERI_Dyz_S_D2x_S_b;
  Double I_ERI_D2z_Pz_D2x_S_b = I_ERI_F3z_S_D2x_S_b+ABZ*I_ERI_D2z_S_D2x_S_b;
  Double I_ERI_D2x_Px_Dxy_S_b = I_ERI_F3x_S_Dxy_S_b+ABX*I_ERI_D2x_S_Dxy_S_b;
  Double I_ERI_Dxy_Px_Dxy_S_b = I_ERI_F2xy_S_Dxy_S_b+ABX*I_ERI_Dxy_S_Dxy_S_b;
  Double I_ERI_Dxz_Px_Dxy_S_b = I_ERI_F2xz_S_Dxy_S_b+ABX*I_ERI_Dxz_S_Dxy_S_b;
  Double I_ERI_D2y_Px_Dxy_S_b = I_ERI_Fx2y_S_Dxy_S_b+ABX*I_ERI_D2y_S_Dxy_S_b;
  Double I_ERI_Dyz_Px_Dxy_S_b = I_ERI_Fxyz_S_Dxy_S_b+ABX*I_ERI_Dyz_S_Dxy_S_b;
  Double I_ERI_D2z_Px_Dxy_S_b = I_ERI_Fx2z_S_Dxy_S_b+ABX*I_ERI_D2z_S_Dxy_S_b;
  Double I_ERI_D2x_Py_Dxy_S_b = I_ERI_F2xy_S_Dxy_S_b+ABY*I_ERI_D2x_S_Dxy_S_b;
  Double I_ERI_Dxy_Py_Dxy_S_b = I_ERI_Fx2y_S_Dxy_S_b+ABY*I_ERI_Dxy_S_Dxy_S_b;
  Double I_ERI_Dxz_Py_Dxy_S_b = I_ERI_Fxyz_S_Dxy_S_b+ABY*I_ERI_Dxz_S_Dxy_S_b;
  Double I_ERI_D2y_Py_Dxy_S_b = I_ERI_F3y_S_Dxy_S_b+ABY*I_ERI_D2y_S_Dxy_S_b;
  Double I_ERI_Dyz_Py_Dxy_S_b = I_ERI_F2yz_S_Dxy_S_b+ABY*I_ERI_Dyz_S_Dxy_S_b;
  Double I_ERI_D2z_Py_Dxy_S_b = I_ERI_Fy2z_S_Dxy_S_b+ABY*I_ERI_D2z_S_Dxy_S_b;
  Double I_ERI_D2x_Pz_Dxy_S_b = I_ERI_F2xz_S_Dxy_S_b+ABZ*I_ERI_D2x_S_Dxy_S_b;
  Double I_ERI_Dxy_Pz_Dxy_S_b = I_ERI_Fxyz_S_Dxy_S_b+ABZ*I_ERI_Dxy_S_Dxy_S_b;
  Double I_ERI_Dxz_Pz_Dxy_S_b = I_ERI_Fx2z_S_Dxy_S_b+ABZ*I_ERI_Dxz_S_Dxy_S_b;
  Double I_ERI_D2y_Pz_Dxy_S_b = I_ERI_F2yz_S_Dxy_S_b+ABZ*I_ERI_D2y_S_Dxy_S_b;
  Double I_ERI_Dyz_Pz_Dxy_S_b = I_ERI_Fy2z_S_Dxy_S_b+ABZ*I_ERI_Dyz_S_Dxy_S_b;
  Double I_ERI_D2z_Pz_Dxy_S_b = I_ERI_F3z_S_Dxy_S_b+ABZ*I_ERI_D2z_S_Dxy_S_b;
  Double I_ERI_D2x_Px_Dxz_S_b = I_ERI_F3x_S_Dxz_S_b+ABX*I_ERI_D2x_S_Dxz_S_b;
  Double I_ERI_Dxy_Px_Dxz_S_b = I_ERI_F2xy_S_Dxz_S_b+ABX*I_ERI_Dxy_S_Dxz_S_b;
  Double I_ERI_Dxz_Px_Dxz_S_b = I_ERI_F2xz_S_Dxz_S_b+ABX*I_ERI_Dxz_S_Dxz_S_b;
  Double I_ERI_D2y_Px_Dxz_S_b = I_ERI_Fx2y_S_Dxz_S_b+ABX*I_ERI_D2y_S_Dxz_S_b;
  Double I_ERI_Dyz_Px_Dxz_S_b = I_ERI_Fxyz_S_Dxz_S_b+ABX*I_ERI_Dyz_S_Dxz_S_b;
  Double I_ERI_D2z_Px_Dxz_S_b = I_ERI_Fx2z_S_Dxz_S_b+ABX*I_ERI_D2z_S_Dxz_S_b;
  Double I_ERI_D2x_Py_Dxz_S_b = I_ERI_F2xy_S_Dxz_S_b+ABY*I_ERI_D2x_S_Dxz_S_b;
  Double I_ERI_Dxy_Py_Dxz_S_b = I_ERI_Fx2y_S_Dxz_S_b+ABY*I_ERI_Dxy_S_Dxz_S_b;
  Double I_ERI_Dxz_Py_Dxz_S_b = I_ERI_Fxyz_S_Dxz_S_b+ABY*I_ERI_Dxz_S_Dxz_S_b;
  Double I_ERI_D2y_Py_Dxz_S_b = I_ERI_F3y_S_Dxz_S_b+ABY*I_ERI_D2y_S_Dxz_S_b;
  Double I_ERI_Dyz_Py_Dxz_S_b = I_ERI_F2yz_S_Dxz_S_b+ABY*I_ERI_Dyz_S_Dxz_S_b;
  Double I_ERI_D2z_Py_Dxz_S_b = I_ERI_Fy2z_S_Dxz_S_b+ABY*I_ERI_D2z_S_Dxz_S_b;
  Double I_ERI_D2x_Pz_Dxz_S_b = I_ERI_F2xz_S_Dxz_S_b+ABZ*I_ERI_D2x_S_Dxz_S_b;
  Double I_ERI_Dxy_Pz_Dxz_S_b = I_ERI_Fxyz_S_Dxz_S_b+ABZ*I_ERI_Dxy_S_Dxz_S_b;
  Double I_ERI_Dxz_Pz_Dxz_S_b = I_ERI_Fx2z_S_Dxz_S_b+ABZ*I_ERI_Dxz_S_Dxz_S_b;
  Double I_ERI_D2y_Pz_Dxz_S_b = I_ERI_F2yz_S_Dxz_S_b+ABZ*I_ERI_D2y_S_Dxz_S_b;
  Double I_ERI_Dyz_Pz_Dxz_S_b = I_ERI_Fy2z_S_Dxz_S_b+ABZ*I_ERI_Dyz_S_Dxz_S_b;
  Double I_ERI_D2z_Pz_Dxz_S_b = I_ERI_F3z_S_Dxz_S_b+ABZ*I_ERI_D2z_S_Dxz_S_b;
  Double I_ERI_D2x_Px_D2y_S_b = I_ERI_F3x_S_D2y_S_b+ABX*I_ERI_D2x_S_D2y_S_b;
  Double I_ERI_Dxy_Px_D2y_S_b = I_ERI_F2xy_S_D2y_S_b+ABX*I_ERI_Dxy_S_D2y_S_b;
  Double I_ERI_Dxz_Px_D2y_S_b = I_ERI_F2xz_S_D2y_S_b+ABX*I_ERI_Dxz_S_D2y_S_b;
  Double I_ERI_D2y_Px_D2y_S_b = I_ERI_Fx2y_S_D2y_S_b+ABX*I_ERI_D2y_S_D2y_S_b;
  Double I_ERI_Dyz_Px_D2y_S_b = I_ERI_Fxyz_S_D2y_S_b+ABX*I_ERI_Dyz_S_D2y_S_b;
  Double I_ERI_D2z_Px_D2y_S_b = I_ERI_Fx2z_S_D2y_S_b+ABX*I_ERI_D2z_S_D2y_S_b;
  Double I_ERI_D2x_Py_D2y_S_b = I_ERI_F2xy_S_D2y_S_b+ABY*I_ERI_D2x_S_D2y_S_b;
  Double I_ERI_Dxy_Py_D2y_S_b = I_ERI_Fx2y_S_D2y_S_b+ABY*I_ERI_Dxy_S_D2y_S_b;
  Double I_ERI_Dxz_Py_D2y_S_b = I_ERI_Fxyz_S_D2y_S_b+ABY*I_ERI_Dxz_S_D2y_S_b;
  Double I_ERI_D2y_Py_D2y_S_b = I_ERI_F3y_S_D2y_S_b+ABY*I_ERI_D2y_S_D2y_S_b;
  Double I_ERI_Dyz_Py_D2y_S_b = I_ERI_F2yz_S_D2y_S_b+ABY*I_ERI_Dyz_S_D2y_S_b;
  Double I_ERI_D2z_Py_D2y_S_b = I_ERI_Fy2z_S_D2y_S_b+ABY*I_ERI_D2z_S_D2y_S_b;
  Double I_ERI_D2x_Pz_D2y_S_b = I_ERI_F2xz_S_D2y_S_b+ABZ*I_ERI_D2x_S_D2y_S_b;
  Double I_ERI_Dxy_Pz_D2y_S_b = I_ERI_Fxyz_S_D2y_S_b+ABZ*I_ERI_Dxy_S_D2y_S_b;
  Double I_ERI_Dxz_Pz_D2y_S_b = I_ERI_Fx2z_S_D2y_S_b+ABZ*I_ERI_Dxz_S_D2y_S_b;
  Double I_ERI_D2y_Pz_D2y_S_b = I_ERI_F2yz_S_D2y_S_b+ABZ*I_ERI_D2y_S_D2y_S_b;
  Double I_ERI_Dyz_Pz_D2y_S_b = I_ERI_Fy2z_S_D2y_S_b+ABZ*I_ERI_Dyz_S_D2y_S_b;
  Double I_ERI_D2z_Pz_D2y_S_b = I_ERI_F3z_S_D2y_S_b+ABZ*I_ERI_D2z_S_D2y_S_b;
  Double I_ERI_D2x_Px_Dyz_S_b = I_ERI_F3x_S_Dyz_S_b+ABX*I_ERI_D2x_S_Dyz_S_b;
  Double I_ERI_Dxy_Px_Dyz_S_b = I_ERI_F2xy_S_Dyz_S_b+ABX*I_ERI_Dxy_S_Dyz_S_b;
  Double I_ERI_Dxz_Px_Dyz_S_b = I_ERI_F2xz_S_Dyz_S_b+ABX*I_ERI_Dxz_S_Dyz_S_b;
  Double I_ERI_D2y_Px_Dyz_S_b = I_ERI_Fx2y_S_Dyz_S_b+ABX*I_ERI_D2y_S_Dyz_S_b;
  Double I_ERI_Dyz_Px_Dyz_S_b = I_ERI_Fxyz_S_Dyz_S_b+ABX*I_ERI_Dyz_S_Dyz_S_b;
  Double I_ERI_D2z_Px_Dyz_S_b = I_ERI_Fx2z_S_Dyz_S_b+ABX*I_ERI_D2z_S_Dyz_S_b;
  Double I_ERI_D2x_Py_Dyz_S_b = I_ERI_F2xy_S_Dyz_S_b+ABY*I_ERI_D2x_S_Dyz_S_b;
  Double I_ERI_Dxy_Py_Dyz_S_b = I_ERI_Fx2y_S_Dyz_S_b+ABY*I_ERI_Dxy_S_Dyz_S_b;
  Double I_ERI_Dxz_Py_Dyz_S_b = I_ERI_Fxyz_S_Dyz_S_b+ABY*I_ERI_Dxz_S_Dyz_S_b;
  Double I_ERI_D2y_Py_Dyz_S_b = I_ERI_F3y_S_Dyz_S_b+ABY*I_ERI_D2y_S_Dyz_S_b;
  Double I_ERI_Dyz_Py_Dyz_S_b = I_ERI_F2yz_S_Dyz_S_b+ABY*I_ERI_Dyz_S_Dyz_S_b;
  Double I_ERI_D2z_Py_Dyz_S_b = I_ERI_Fy2z_S_Dyz_S_b+ABY*I_ERI_D2z_S_Dyz_S_b;
  Double I_ERI_D2x_Pz_Dyz_S_b = I_ERI_F2xz_S_Dyz_S_b+ABZ*I_ERI_D2x_S_Dyz_S_b;
  Double I_ERI_Dxy_Pz_Dyz_S_b = I_ERI_Fxyz_S_Dyz_S_b+ABZ*I_ERI_Dxy_S_Dyz_S_b;
  Double I_ERI_Dxz_Pz_Dyz_S_b = I_ERI_Fx2z_S_Dyz_S_b+ABZ*I_ERI_Dxz_S_Dyz_S_b;
  Double I_ERI_D2y_Pz_Dyz_S_b = I_ERI_F2yz_S_Dyz_S_b+ABZ*I_ERI_D2y_S_Dyz_S_b;
  Double I_ERI_Dyz_Pz_Dyz_S_b = I_ERI_Fy2z_S_Dyz_S_b+ABZ*I_ERI_Dyz_S_Dyz_S_b;
  Double I_ERI_D2z_Pz_Dyz_S_b = I_ERI_F3z_S_Dyz_S_b+ABZ*I_ERI_D2z_S_Dyz_S_b;
  Double I_ERI_D2x_Px_D2z_S_b = I_ERI_F3x_S_D2z_S_b+ABX*I_ERI_D2x_S_D2z_S_b;
  Double I_ERI_Dxy_Px_D2z_S_b = I_ERI_F2xy_S_D2z_S_b+ABX*I_ERI_Dxy_S_D2z_S_b;
  Double I_ERI_Dxz_Px_D2z_S_b = I_ERI_F2xz_S_D2z_S_b+ABX*I_ERI_Dxz_S_D2z_S_b;
  Double I_ERI_D2y_Px_D2z_S_b = I_ERI_Fx2y_S_D2z_S_b+ABX*I_ERI_D2y_S_D2z_S_b;
  Double I_ERI_Dyz_Px_D2z_S_b = I_ERI_Fxyz_S_D2z_S_b+ABX*I_ERI_Dyz_S_D2z_S_b;
  Double I_ERI_D2z_Px_D2z_S_b = I_ERI_Fx2z_S_D2z_S_b+ABX*I_ERI_D2z_S_D2z_S_b;
  Double I_ERI_D2x_Py_D2z_S_b = I_ERI_F2xy_S_D2z_S_b+ABY*I_ERI_D2x_S_D2z_S_b;
  Double I_ERI_Dxy_Py_D2z_S_b = I_ERI_Fx2y_S_D2z_S_b+ABY*I_ERI_Dxy_S_D2z_S_b;
  Double I_ERI_Dxz_Py_D2z_S_b = I_ERI_Fxyz_S_D2z_S_b+ABY*I_ERI_Dxz_S_D2z_S_b;
  Double I_ERI_D2y_Py_D2z_S_b = I_ERI_F3y_S_D2z_S_b+ABY*I_ERI_D2y_S_D2z_S_b;
  Double I_ERI_Dyz_Py_D2z_S_b = I_ERI_F2yz_S_D2z_S_b+ABY*I_ERI_Dyz_S_D2z_S_b;
  Double I_ERI_D2z_Py_D2z_S_b = I_ERI_Fy2z_S_D2z_S_b+ABY*I_ERI_D2z_S_D2z_S_b;
  Double I_ERI_D2x_Pz_D2z_S_b = I_ERI_F2xz_S_D2z_S_b+ABZ*I_ERI_D2x_S_D2z_S_b;
  Double I_ERI_Dxy_Pz_D2z_S_b = I_ERI_Fxyz_S_D2z_S_b+ABZ*I_ERI_Dxy_S_D2z_S_b;
  Double I_ERI_Dxz_Pz_D2z_S_b = I_ERI_Fx2z_S_D2z_S_b+ABZ*I_ERI_Dxz_S_D2z_S_b;
  Double I_ERI_D2y_Pz_D2z_S_b = I_ERI_F2yz_S_D2z_S_b+ABZ*I_ERI_D2y_S_D2z_S_b;
  Double I_ERI_Dyz_Pz_D2z_S_b = I_ERI_Fy2z_S_D2z_S_b+ABZ*I_ERI_Dyz_S_D2z_S_b;
  Double I_ERI_D2z_Pz_D2z_S_b = I_ERI_F3z_S_D2z_S_b+ABZ*I_ERI_D2z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_b
   * RHS shell quartet name: SQ_ERI_F_S_D_S_b
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_b = I_ERI_G4x_S_D2x_S_b+ABX*I_ERI_F3x_S_D2x_S_b;
  Double I_ERI_F2xy_Px_D2x_S_b = I_ERI_G3xy_S_D2x_S_b+ABX*I_ERI_F2xy_S_D2x_S_b;
  Double I_ERI_F2xz_Px_D2x_S_b = I_ERI_G3xz_S_D2x_S_b+ABX*I_ERI_F2xz_S_D2x_S_b;
  Double I_ERI_Fx2y_Px_D2x_S_b = I_ERI_G2x2y_S_D2x_S_b+ABX*I_ERI_Fx2y_S_D2x_S_b;
  Double I_ERI_Fxyz_Px_D2x_S_b = I_ERI_G2xyz_S_D2x_S_b+ABX*I_ERI_Fxyz_S_D2x_S_b;
  Double I_ERI_Fx2z_Px_D2x_S_b = I_ERI_G2x2z_S_D2x_S_b+ABX*I_ERI_Fx2z_S_D2x_S_b;
  Double I_ERI_F3y_Px_D2x_S_b = I_ERI_Gx3y_S_D2x_S_b+ABX*I_ERI_F3y_S_D2x_S_b;
  Double I_ERI_F2yz_Px_D2x_S_b = I_ERI_Gx2yz_S_D2x_S_b+ABX*I_ERI_F2yz_S_D2x_S_b;
  Double I_ERI_Fy2z_Px_D2x_S_b = I_ERI_Gxy2z_S_D2x_S_b+ABX*I_ERI_Fy2z_S_D2x_S_b;
  Double I_ERI_F3z_Px_D2x_S_b = I_ERI_Gx3z_S_D2x_S_b+ABX*I_ERI_F3z_S_D2x_S_b;
  Double I_ERI_F3x_Py_D2x_S_b = I_ERI_G3xy_S_D2x_S_b+ABY*I_ERI_F3x_S_D2x_S_b;
  Double I_ERI_F2xy_Py_D2x_S_b = I_ERI_G2x2y_S_D2x_S_b+ABY*I_ERI_F2xy_S_D2x_S_b;
  Double I_ERI_F2xz_Py_D2x_S_b = I_ERI_G2xyz_S_D2x_S_b+ABY*I_ERI_F2xz_S_D2x_S_b;
  Double I_ERI_Fx2y_Py_D2x_S_b = I_ERI_Gx3y_S_D2x_S_b+ABY*I_ERI_Fx2y_S_D2x_S_b;
  Double I_ERI_Fxyz_Py_D2x_S_b = I_ERI_Gx2yz_S_D2x_S_b+ABY*I_ERI_Fxyz_S_D2x_S_b;
  Double I_ERI_Fx2z_Py_D2x_S_b = I_ERI_Gxy2z_S_D2x_S_b+ABY*I_ERI_Fx2z_S_D2x_S_b;
  Double I_ERI_F3y_Py_D2x_S_b = I_ERI_G4y_S_D2x_S_b+ABY*I_ERI_F3y_S_D2x_S_b;
  Double I_ERI_F2yz_Py_D2x_S_b = I_ERI_G3yz_S_D2x_S_b+ABY*I_ERI_F2yz_S_D2x_S_b;
  Double I_ERI_Fy2z_Py_D2x_S_b = I_ERI_G2y2z_S_D2x_S_b+ABY*I_ERI_Fy2z_S_D2x_S_b;
  Double I_ERI_F3z_Py_D2x_S_b = I_ERI_Gy3z_S_D2x_S_b+ABY*I_ERI_F3z_S_D2x_S_b;
  Double I_ERI_F3x_Pz_D2x_S_b = I_ERI_G3xz_S_D2x_S_b+ABZ*I_ERI_F3x_S_D2x_S_b;
  Double I_ERI_F2xy_Pz_D2x_S_b = I_ERI_G2xyz_S_D2x_S_b+ABZ*I_ERI_F2xy_S_D2x_S_b;
  Double I_ERI_F2xz_Pz_D2x_S_b = I_ERI_G2x2z_S_D2x_S_b+ABZ*I_ERI_F2xz_S_D2x_S_b;
  Double I_ERI_Fx2y_Pz_D2x_S_b = I_ERI_Gx2yz_S_D2x_S_b+ABZ*I_ERI_Fx2y_S_D2x_S_b;
  Double I_ERI_Fxyz_Pz_D2x_S_b = I_ERI_Gxy2z_S_D2x_S_b+ABZ*I_ERI_Fxyz_S_D2x_S_b;
  Double I_ERI_Fx2z_Pz_D2x_S_b = I_ERI_Gx3z_S_D2x_S_b+ABZ*I_ERI_Fx2z_S_D2x_S_b;
  Double I_ERI_F3y_Pz_D2x_S_b = I_ERI_G3yz_S_D2x_S_b+ABZ*I_ERI_F3y_S_D2x_S_b;
  Double I_ERI_F2yz_Pz_D2x_S_b = I_ERI_G2y2z_S_D2x_S_b+ABZ*I_ERI_F2yz_S_D2x_S_b;
  Double I_ERI_Fy2z_Pz_D2x_S_b = I_ERI_Gy3z_S_D2x_S_b+ABZ*I_ERI_Fy2z_S_D2x_S_b;
  Double I_ERI_F3z_Pz_D2x_S_b = I_ERI_G4z_S_D2x_S_b+ABZ*I_ERI_F3z_S_D2x_S_b;
  Double I_ERI_F3x_Px_Dxy_S_b = I_ERI_G4x_S_Dxy_S_b+ABX*I_ERI_F3x_S_Dxy_S_b;
  Double I_ERI_F2xy_Px_Dxy_S_b = I_ERI_G3xy_S_Dxy_S_b+ABX*I_ERI_F2xy_S_Dxy_S_b;
  Double I_ERI_F2xz_Px_Dxy_S_b = I_ERI_G3xz_S_Dxy_S_b+ABX*I_ERI_F2xz_S_Dxy_S_b;
  Double I_ERI_Fx2y_Px_Dxy_S_b = I_ERI_G2x2y_S_Dxy_S_b+ABX*I_ERI_Fx2y_S_Dxy_S_b;
  Double I_ERI_Fxyz_Px_Dxy_S_b = I_ERI_G2xyz_S_Dxy_S_b+ABX*I_ERI_Fxyz_S_Dxy_S_b;
  Double I_ERI_Fx2z_Px_Dxy_S_b = I_ERI_G2x2z_S_Dxy_S_b+ABX*I_ERI_Fx2z_S_Dxy_S_b;
  Double I_ERI_F3y_Px_Dxy_S_b = I_ERI_Gx3y_S_Dxy_S_b+ABX*I_ERI_F3y_S_Dxy_S_b;
  Double I_ERI_F2yz_Px_Dxy_S_b = I_ERI_Gx2yz_S_Dxy_S_b+ABX*I_ERI_F2yz_S_Dxy_S_b;
  Double I_ERI_Fy2z_Px_Dxy_S_b = I_ERI_Gxy2z_S_Dxy_S_b+ABX*I_ERI_Fy2z_S_Dxy_S_b;
  Double I_ERI_F3z_Px_Dxy_S_b = I_ERI_Gx3z_S_Dxy_S_b+ABX*I_ERI_F3z_S_Dxy_S_b;
  Double I_ERI_F3x_Py_Dxy_S_b = I_ERI_G3xy_S_Dxy_S_b+ABY*I_ERI_F3x_S_Dxy_S_b;
  Double I_ERI_F2xy_Py_Dxy_S_b = I_ERI_G2x2y_S_Dxy_S_b+ABY*I_ERI_F2xy_S_Dxy_S_b;
  Double I_ERI_F2xz_Py_Dxy_S_b = I_ERI_G2xyz_S_Dxy_S_b+ABY*I_ERI_F2xz_S_Dxy_S_b;
  Double I_ERI_Fx2y_Py_Dxy_S_b = I_ERI_Gx3y_S_Dxy_S_b+ABY*I_ERI_Fx2y_S_Dxy_S_b;
  Double I_ERI_Fxyz_Py_Dxy_S_b = I_ERI_Gx2yz_S_Dxy_S_b+ABY*I_ERI_Fxyz_S_Dxy_S_b;
  Double I_ERI_Fx2z_Py_Dxy_S_b = I_ERI_Gxy2z_S_Dxy_S_b+ABY*I_ERI_Fx2z_S_Dxy_S_b;
  Double I_ERI_F3y_Py_Dxy_S_b = I_ERI_G4y_S_Dxy_S_b+ABY*I_ERI_F3y_S_Dxy_S_b;
  Double I_ERI_F2yz_Py_Dxy_S_b = I_ERI_G3yz_S_Dxy_S_b+ABY*I_ERI_F2yz_S_Dxy_S_b;
  Double I_ERI_Fy2z_Py_Dxy_S_b = I_ERI_G2y2z_S_Dxy_S_b+ABY*I_ERI_Fy2z_S_Dxy_S_b;
  Double I_ERI_F3z_Py_Dxy_S_b = I_ERI_Gy3z_S_Dxy_S_b+ABY*I_ERI_F3z_S_Dxy_S_b;
  Double I_ERI_F3x_Pz_Dxy_S_b = I_ERI_G3xz_S_Dxy_S_b+ABZ*I_ERI_F3x_S_Dxy_S_b;
  Double I_ERI_F2xy_Pz_Dxy_S_b = I_ERI_G2xyz_S_Dxy_S_b+ABZ*I_ERI_F2xy_S_Dxy_S_b;
  Double I_ERI_F2xz_Pz_Dxy_S_b = I_ERI_G2x2z_S_Dxy_S_b+ABZ*I_ERI_F2xz_S_Dxy_S_b;
  Double I_ERI_Fx2y_Pz_Dxy_S_b = I_ERI_Gx2yz_S_Dxy_S_b+ABZ*I_ERI_Fx2y_S_Dxy_S_b;
  Double I_ERI_Fxyz_Pz_Dxy_S_b = I_ERI_Gxy2z_S_Dxy_S_b+ABZ*I_ERI_Fxyz_S_Dxy_S_b;
  Double I_ERI_Fx2z_Pz_Dxy_S_b = I_ERI_Gx3z_S_Dxy_S_b+ABZ*I_ERI_Fx2z_S_Dxy_S_b;
  Double I_ERI_F3y_Pz_Dxy_S_b = I_ERI_G3yz_S_Dxy_S_b+ABZ*I_ERI_F3y_S_Dxy_S_b;
  Double I_ERI_F2yz_Pz_Dxy_S_b = I_ERI_G2y2z_S_Dxy_S_b+ABZ*I_ERI_F2yz_S_Dxy_S_b;
  Double I_ERI_Fy2z_Pz_Dxy_S_b = I_ERI_Gy3z_S_Dxy_S_b+ABZ*I_ERI_Fy2z_S_Dxy_S_b;
  Double I_ERI_F3z_Pz_Dxy_S_b = I_ERI_G4z_S_Dxy_S_b+ABZ*I_ERI_F3z_S_Dxy_S_b;
  Double I_ERI_F3x_Px_Dxz_S_b = I_ERI_G4x_S_Dxz_S_b+ABX*I_ERI_F3x_S_Dxz_S_b;
  Double I_ERI_F2xy_Px_Dxz_S_b = I_ERI_G3xy_S_Dxz_S_b+ABX*I_ERI_F2xy_S_Dxz_S_b;
  Double I_ERI_F2xz_Px_Dxz_S_b = I_ERI_G3xz_S_Dxz_S_b+ABX*I_ERI_F2xz_S_Dxz_S_b;
  Double I_ERI_Fx2y_Px_Dxz_S_b = I_ERI_G2x2y_S_Dxz_S_b+ABX*I_ERI_Fx2y_S_Dxz_S_b;
  Double I_ERI_Fxyz_Px_Dxz_S_b = I_ERI_G2xyz_S_Dxz_S_b+ABX*I_ERI_Fxyz_S_Dxz_S_b;
  Double I_ERI_Fx2z_Px_Dxz_S_b = I_ERI_G2x2z_S_Dxz_S_b+ABX*I_ERI_Fx2z_S_Dxz_S_b;
  Double I_ERI_F3y_Px_Dxz_S_b = I_ERI_Gx3y_S_Dxz_S_b+ABX*I_ERI_F3y_S_Dxz_S_b;
  Double I_ERI_F2yz_Px_Dxz_S_b = I_ERI_Gx2yz_S_Dxz_S_b+ABX*I_ERI_F2yz_S_Dxz_S_b;
  Double I_ERI_Fy2z_Px_Dxz_S_b = I_ERI_Gxy2z_S_Dxz_S_b+ABX*I_ERI_Fy2z_S_Dxz_S_b;
  Double I_ERI_F3z_Px_Dxz_S_b = I_ERI_Gx3z_S_Dxz_S_b+ABX*I_ERI_F3z_S_Dxz_S_b;
  Double I_ERI_F3x_Py_Dxz_S_b = I_ERI_G3xy_S_Dxz_S_b+ABY*I_ERI_F3x_S_Dxz_S_b;
  Double I_ERI_F2xy_Py_Dxz_S_b = I_ERI_G2x2y_S_Dxz_S_b+ABY*I_ERI_F2xy_S_Dxz_S_b;
  Double I_ERI_F2xz_Py_Dxz_S_b = I_ERI_G2xyz_S_Dxz_S_b+ABY*I_ERI_F2xz_S_Dxz_S_b;
  Double I_ERI_Fx2y_Py_Dxz_S_b = I_ERI_Gx3y_S_Dxz_S_b+ABY*I_ERI_Fx2y_S_Dxz_S_b;
  Double I_ERI_Fxyz_Py_Dxz_S_b = I_ERI_Gx2yz_S_Dxz_S_b+ABY*I_ERI_Fxyz_S_Dxz_S_b;
  Double I_ERI_Fx2z_Py_Dxz_S_b = I_ERI_Gxy2z_S_Dxz_S_b+ABY*I_ERI_Fx2z_S_Dxz_S_b;
  Double I_ERI_F3y_Py_Dxz_S_b = I_ERI_G4y_S_Dxz_S_b+ABY*I_ERI_F3y_S_Dxz_S_b;
  Double I_ERI_F2yz_Py_Dxz_S_b = I_ERI_G3yz_S_Dxz_S_b+ABY*I_ERI_F2yz_S_Dxz_S_b;
  Double I_ERI_Fy2z_Py_Dxz_S_b = I_ERI_G2y2z_S_Dxz_S_b+ABY*I_ERI_Fy2z_S_Dxz_S_b;
  Double I_ERI_F3z_Py_Dxz_S_b = I_ERI_Gy3z_S_Dxz_S_b+ABY*I_ERI_F3z_S_Dxz_S_b;
  Double I_ERI_F3x_Pz_Dxz_S_b = I_ERI_G3xz_S_Dxz_S_b+ABZ*I_ERI_F3x_S_Dxz_S_b;
  Double I_ERI_F2xy_Pz_Dxz_S_b = I_ERI_G2xyz_S_Dxz_S_b+ABZ*I_ERI_F2xy_S_Dxz_S_b;
  Double I_ERI_F2xz_Pz_Dxz_S_b = I_ERI_G2x2z_S_Dxz_S_b+ABZ*I_ERI_F2xz_S_Dxz_S_b;
  Double I_ERI_Fx2y_Pz_Dxz_S_b = I_ERI_Gx2yz_S_Dxz_S_b+ABZ*I_ERI_Fx2y_S_Dxz_S_b;
  Double I_ERI_Fxyz_Pz_Dxz_S_b = I_ERI_Gxy2z_S_Dxz_S_b+ABZ*I_ERI_Fxyz_S_Dxz_S_b;
  Double I_ERI_Fx2z_Pz_Dxz_S_b = I_ERI_Gx3z_S_Dxz_S_b+ABZ*I_ERI_Fx2z_S_Dxz_S_b;
  Double I_ERI_F3y_Pz_Dxz_S_b = I_ERI_G3yz_S_Dxz_S_b+ABZ*I_ERI_F3y_S_Dxz_S_b;
  Double I_ERI_F2yz_Pz_Dxz_S_b = I_ERI_G2y2z_S_Dxz_S_b+ABZ*I_ERI_F2yz_S_Dxz_S_b;
  Double I_ERI_Fy2z_Pz_Dxz_S_b = I_ERI_Gy3z_S_Dxz_S_b+ABZ*I_ERI_Fy2z_S_Dxz_S_b;
  Double I_ERI_F3z_Pz_Dxz_S_b = I_ERI_G4z_S_Dxz_S_b+ABZ*I_ERI_F3z_S_Dxz_S_b;
  Double I_ERI_F3x_Px_D2y_S_b = I_ERI_G4x_S_D2y_S_b+ABX*I_ERI_F3x_S_D2y_S_b;
  Double I_ERI_F2xy_Px_D2y_S_b = I_ERI_G3xy_S_D2y_S_b+ABX*I_ERI_F2xy_S_D2y_S_b;
  Double I_ERI_F2xz_Px_D2y_S_b = I_ERI_G3xz_S_D2y_S_b+ABX*I_ERI_F2xz_S_D2y_S_b;
  Double I_ERI_Fx2y_Px_D2y_S_b = I_ERI_G2x2y_S_D2y_S_b+ABX*I_ERI_Fx2y_S_D2y_S_b;
  Double I_ERI_Fxyz_Px_D2y_S_b = I_ERI_G2xyz_S_D2y_S_b+ABX*I_ERI_Fxyz_S_D2y_S_b;
  Double I_ERI_Fx2z_Px_D2y_S_b = I_ERI_G2x2z_S_D2y_S_b+ABX*I_ERI_Fx2z_S_D2y_S_b;
  Double I_ERI_F3y_Px_D2y_S_b = I_ERI_Gx3y_S_D2y_S_b+ABX*I_ERI_F3y_S_D2y_S_b;
  Double I_ERI_F2yz_Px_D2y_S_b = I_ERI_Gx2yz_S_D2y_S_b+ABX*I_ERI_F2yz_S_D2y_S_b;
  Double I_ERI_Fy2z_Px_D2y_S_b = I_ERI_Gxy2z_S_D2y_S_b+ABX*I_ERI_Fy2z_S_D2y_S_b;
  Double I_ERI_F3z_Px_D2y_S_b = I_ERI_Gx3z_S_D2y_S_b+ABX*I_ERI_F3z_S_D2y_S_b;
  Double I_ERI_F3x_Py_D2y_S_b = I_ERI_G3xy_S_D2y_S_b+ABY*I_ERI_F3x_S_D2y_S_b;
  Double I_ERI_F2xy_Py_D2y_S_b = I_ERI_G2x2y_S_D2y_S_b+ABY*I_ERI_F2xy_S_D2y_S_b;
  Double I_ERI_F2xz_Py_D2y_S_b = I_ERI_G2xyz_S_D2y_S_b+ABY*I_ERI_F2xz_S_D2y_S_b;
  Double I_ERI_Fx2y_Py_D2y_S_b = I_ERI_Gx3y_S_D2y_S_b+ABY*I_ERI_Fx2y_S_D2y_S_b;
  Double I_ERI_Fxyz_Py_D2y_S_b = I_ERI_Gx2yz_S_D2y_S_b+ABY*I_ERI_Fxyz_S_D2y_S_b;
  Double I_ERI_Fx2z_Py_D2y_S_b = I_ERI_Gxy2z_S_D2y_S_b+ABY*I_ERI_Fx2z_S_D2y_S_b;
  Double I_ERI_F3y_Py_D2y_S_b = I_ERI_G4y_S_D2y_S_b+ABY*I_ERI_F3y_S_D2y_S_b;
  Double I_ERI_F2yz_Py_D2y_S_b = I_ERI_G3yz_S_D2y_S_b+ABY*I_ERI_F2yz_S_D2y_S_b;
  Double I_ERI_Fy2z_Py_D2y_S_b = I_ERI_G2y2z_S_D2y_S_b+ABY*I_ERI_Fy2z_S_D2y_S_b;
  Double I_ERI_F3z_Py_D2y_S_b = I_ERI_Gy3z_S_D2y_S_b+ABY*I_ERI_F3z_S_D2y_S_b;
  Double I_ERI_F3x_Pz_D2y_S_b = I_ERI_G3xz_S_D2y_S_b+ABZ*I_ERI_F3x_S_D2y_S_b;
  Double I_ERI_F2xy_Pz_D2y_S_b = I_ERI_G2xyz_S_D2y_S_b+ABZ*I_ERI_F2xy_S_D2y_S_b;
  Double I_ERI_F2xz_Pz_D2y_S_b = I_ERI_G2x2z_S_D2y_S_b+ABZ*I_ERI_F2xz_S_D2y_S_b;
  Double I_ERI_Fx2y_Pz_D2y_S_b = I_ERI_Gx2yz_S_D2y_S_b+ABZ*I_ERI_Fx2y_S_D2y_S_b;
  Double I_ERI_Fxyz_Pz_D2y_S_b = I_ERI_Gxy2z_S_D2y_S_b+ABZ*I_ERI_Fxyz_S_D2y_S_b;
  Double I_ERI_Fx2z_Pz_D2y_S_b = I_ERI_Gx3z_S_D2y_S_b+ABZ*I_ERI_Fx2z_S_D2y_S_b;
  Double I_ERI_F3y_Pz_D2y_S_b = I_ERI_G3yz_S_D2y_S_b+ABZ*I_ERI_F3y_S_D2y_S_b;
  Double I_ERI_F2yz_Pz_D2y_S_b = I_ERI_G2y2z_S_D2y_S_b+ABZ*I_ERI_F2yz_S_D2y_S_b;
  Double I_ERI_Fy2z_Pz_D2y_S_b = I_ERI_Gy3z_S_D2y_S_b+ABZ*I_ERI_Fy2z_S_D2y_S_b;
  Double I_ERI_F3z_Pz_D2y_S_b = I_ERI_G4z_S_D2y_S_b+ABZ*I_ERI_F3z_S_D2y_S_b;
  Double I_ERI_F3x_Px_Dyz_S_b = I_ERI_G4x_S_Dyz_S_b+ABX*I_ERI_F3x_S_Dyz_S_b;
  Double I_ERI_F2xy_Px_Dyz_S_b = I_ERI_G3xy_S_Dyz_S_b+ABX*I_ERI_F2xy_S_Dyz_S_b;
  Double I_ERI_F2xz_Px_Dyz_S_b = I_ERI_G3xz_S_Dyz_S_b+ABX*I_ERI_F2xz_S_Dyz_S_b;
  Double I_ERI_Fx2y_Px_Dyz_S_b = I_ERI_G2x2y_S_Dyz_S_b+ABX*I_ERI_Fx2y_S_Dyz_S_b;
  Double I_ERI_Fxyz_Px_Dyz_S_b = I_ERI_G2xyz_S_Dyz_S_b+ABX*I_ERI_Fxyz_S_Dyz_S_b;
  Double I_ERI_Fx2z_Px_Dyz_S_b = I_ERI_G2x2z_S_Dyz_S_b+ABX*I_ERI_Fx2z_S_Dyz_S_b;
  Double I_ERI_F3y_Px_Dyz_S_b = I_ERI_Gx3y_S_Dyz_S_b+ABX*I_ERI_F3y_S_Dyz_S_b;
  Double I_ERI_F2yz_Px_Dyz_S_b = I_ERI_Gx2yz_S_Dyz_S_b+ABX*I_ERI_F2yz_S_Dyz_S_b;
  Double I_ERI_Fy2z_Px_Dyz_S_b = I_ERI_Gxy2z_S_Dyz_S_b+ABX*I_ERI_Fy2z_S_Dyz_S_b;
  Double I_ERI_F3z_Px_Dyz_S_b = I_ERI_Gx3z_S_Dyz_S_b+ABX*I_ERI_F3z_S_Dyz_S_b;
  Double I_ERI_F3x_Py_Dyz_S_b = I_ERI_G3xy_S_Dyz_S_b+ABY*I_ERI_F3x_S_Dyz_S_b;
  Double I_ERI_F2xy_Py_Dyz_S_b = I_ERI_G2x2y_S_Dyz_S_b+ABY*I_ERI_F2xy_S_Dyz_S_b;
  Double I_ERI_F2xz_Py_Dyz_S_b = I_ERI_G2xyz_S_Dyz_S_b+ABY*I_ERI_F2xz_S_Dyz_S_b;
  Double I_ERI_Fx2y_Py_Dyz_S_b = I_ERI_Gx3y_S_Dyz_S_b+ABY*I_ERI_Fx2y_S_Dyz_S_b;
  Double I_ERI_Fxyz_Py_Dyz_S_b = I_ERI_Gx2yz_S_Dyz_S_b+ABY*I_ERI_Fxyz_S_Dyz_S_b;
  Double I_ERI_Fx2z_Py_Dyz_S_b = I_ERI_Gxy2z_S_Dyz_S_b+ABY*I_ERI_Fx2z_S_Dyz_S_b;
  Double I_ERI_F3y_Py_Dyz_S_b = I_ERI_G4y_S_Dyz_S_b+ABY*I_ERI_F3y_S_Dyz_S_b;
  Double I_ERI_F2yz_Py_Dyz_S_b = I_ERI_G3yz_S_Dyz_S_b+ABY*I_ERI_F2yz_S_Dyz_S_b;
  Double I_ERI_Fy2z_Py_Dyz_S_b = I_ERI_G2y2z_S_Dyz_S_b+ABY*I_ERI_Fy2z_S_Dyz_S_b;
  Double I_ERI_F3z_Py_Dyz_S_b = I_ERI_Gy3z_S_Dyz_S_b+ABY*I_ERI_F3z_S_Dyz_S_b;
  Double I_ERI_F3x_Pz_Dyz_S_b = I_ERI_G3xz_S_Dyz_S_b+ABZ*I_ERI_F3x_S_Dyz_S_b;
  Double I_ERI_F2xy_Pz_Dyz_S_b = I_ERI_G2xyz_S_Dyz_S_b+ABZ*I_ERI_F2xy_S_Dyz_S_b;
  Double I_ERI_F2xz_Pz_Dyz_S_b = I_ERI_G2x2z_S_Dyz_S_b+ABZ*I_ERI_F2xz_S_Dyz_S_b;
  Double I_ERI_Fx2y_Pz_Dyz_S_b = I_ERI_Gx2yz_S_Dyz_S_b+ABZ*I_ERI_Fx2y_S_Dyz_S_b;
  Double I_ERI_Fxyz_Pz_Dyz_S_b = I_ERI_Gxy2z_S_Dyz_S_b+ABZ*I_ERI_Fxyz_S_Dyz_S_b;
  Double I_ERI_Fx2z_Pz_Dyz_S_b = I_ERI_Gx3z_S_Dyz_S_b+ABZ*I_ERI_Fx2z_S_Dyz_S_b;
  Double I_ERI_F3y_Pz_Dyz_S_b = I_ERI_G3yz_S_Dyz_S_b+ABZ*I_ERI_F3y_S_Dyz_S_b;
  Double I_ERI_F2yz_Pz_Dyz_S_b = I_ERI_G2y2z_S_Dyz_S_b+ABZ*I_ERI_F2yz_S_Dyz_S_b;
  Double I_ERI_Fy2z_Pz_Dyz_S_b = I_ERI_Gy3z_S_Dyz_S_b+ABZ*I_ERI_Fy2z_S_Dyz_S_b;
  Double I_ERI_F3z_Pz_Dyz_S_b = I_ERI_G4z_S_Dyz_S_b+ABZ*I_ERI_F3z_S_Dyz_S_b;
  Double I_ERI_F3x_Px_D2z_S_b = I_ERI_G4x_S_D2z_S_b+ABX*I_ERI_F3x_S_D2z_S_b;
  Double I_ERI_F2xy_Px_D2z_S_b = I_ERI_G3xy_S_D2z_S_b+ABX*I_ERI_F2xy_S_D2z_S_b;
  Double I_ERI_F2xz_Px_D2z_S_b = I_ERI_G3xz_S_D2z_S_b+ABX*I_ERI_F2xz_S_D2z_S_b;
  Double I_ERI_Fx2y_Px_D2z_S_b = I_ERI_G2x2y_S_D2z_S_b+ABX*I_ERI_Fx2y_S_D2z_S_b;
  Double I_ERI_Fxyz_Px_D2z_S_b = I_ERI_G2xyz_S_D2z_S_b+ABX*I_ERI_Fxyz_S_D2z_S_b;
  Double I_ERI_Fx2z_Px_D2z_S_b = I_ERI_G2x2z_S_D2z_S_b+ABX*I_ERI_Fx2z_S_D2z_S_b;
  Double I_ERI_F3y_Px_D2z_S_b = I_ERI_Gx3y_S_D2z_S_b+ABX*I_ERI_F3y_S_D2z_S_b;
  Double I_ERI_F2yz_Px_D2z_S_b = I_ERI_Gx2yz_S_D2z_S_b+ABX*I_ERI_F2yz_S_D2z_S_b;
  Double I_ERI_Fy2z_Px_D2z_S_b = I_ERI_Gxy2z_S_D2z_S_b+ABX*I_ERI_Fy2z_S_D2z_S_b;
  Double I_ERI_F3z_Px_D2z_S_b = I_ERI_Gx3z_S_D2z_S_b+ABX*I_ERI_F3z_S_D2z_S_b;
  Double I_ERI_F3x_Py_D2z_S_b = I_ERI_G3xy_S_D2z_S_b+ABY*I_ERI_F3x_S_D2z_S_b;
  Double I_ERI_F2xy_Py_D2z_S_b = I_ERI_G2x2y_S_D2z_S_b+ABY*I_ERI_F2xy_S_D2z_S_b;
  Double I_ERI_F2xz_Py_D2z_S_b = I_ERI_G2xyz_S_D2z_S_b+ABY*I_ERI_F2xz_S_D2z_S_b;
  Double I_ERI_Fx2y_Py_D2z_S_b = I_ERI_Gx3y_S_D2z_S_b+ABY*I_ERI_Fx2y_S_D2z_S_b;
  Double I_ERI_Fxyz_Py_D2z_S_b = I_ERI_Gx2yz_S_D2z_S_b+ABY*I_ERI_Fxyz_S_D2z_S_b;
  Double I_ERI_Fx2z_Py_D2z_S_b = I_ERI_Gxy2z_S_D2z_S_b+ABY*I_ERI_Fx2z_S_D2z_S_b;
  Double I_ERI_F3y_Py_D2z_S_b = I_ERI_G4y_S_D2z_S_b+ABY*I_ERI_F3y_S_D2z_S_b;
  Double I_ERI_F2yz_Py_D2z_S_b = I_ERI_G3yz_S_D2z_S_b+ABY*I_ERI_F2yz_S_D2z_S_b;
  Double I_ERI_Fy2z_Py_D2z_S_b = I_ERI_G2y2z_S_D2z_S_b+ABY*I_ERI_Fy2z_S_D2z_S_b;
  Double I_ERI_F3z_Py_D2z_S_b = I_ERI_Gy3z_S_D2z_S_b+ABY*I_ERI_F3z_S_D2z_S_b;
  Double I_ERI_F3x_Pz_D2z_S_b = I_ERI_G3xz_S_D2z_S_b+ABZ*I_ERI_F3x_S_D2z_S_b;
  Double I_ERI_F2xy_Pz_D2z_S_b = I_ERI_G2xyz_S_D2z_S_b+ABZ*I_ERI_F2xy_S_D2z_S_b;
  Double I_ERI_F2xz_Pz_D2z_S_b = I_ERI_G2x2z_S_D2z_S_b+ABZ*I_ERI_F2xz_S_D2z_S_b;
  Double I_ERI_Fx2y_Pz_D2z_S_b = I_ERI_Gx2yz_S_D2z_S_b+ABZ*I_ERI_Fx2y_S_D2z_S_b;
  Double I_ERI_Fxyz_Pz_D2z_S_b = I_ERI_Gxy2z_S_D2z_S_b+ABZ*I_ERI_Fxyz_S_D2z_S_b;
  Double I_ERI_Fx2z_Pz_D2z_S_b = I_ERI_Gx3z_S_D2z_S_b+ABZ*I_ERI_Fx2z_S_D2z_S_b;
  Double I_ERI_F3y_Pz_D2z_S_b = I_ERI_G3yz_S_D2z_S_b+ABZ*I_ERI_F3y_S_D2z_S_b;
  Double I_ERI_F2yz_Pz_D2z_S_b = I_ERI_G2y2z_S_D2z_S_b+ABZ*I_ERI_F2yz_S_D2z_S_b;
  Double I_ERI_Fy2z_Pz_D2z_S_b = I_ERI_Gy3z_S_D2z_S_b+ABZ*I_ERI_Fy2z_S_D2z_S_b;
  Double I_ERI_F3z_Pz_D2z_S_b = I_ERI_G4z_S_D2z_S_b+ABZ*I_ERI_F3z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 72 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_b
   * RHS shell quartet name: SQ_ERI_D_P_D_S_b
   ************************************************************/
  Double I_ERI_D2x_D2x_D2x_S_b = I_ERI_F3x_Px_D2x_S_b+ABX*I_ERI_D2x_Px_D2x_S_b;
  Double I_ERI_Dxy_D2x_D2x_S_b = I_ERI_F2xy_Px_D2x_S_b+ABX*I_ERI_Dxy_Px_D2x_S_b;
  Double I_ERI_Dxz_D2x_D2x_S_b = I_ERI_F2xz_Px_D2x_S_b+ABX*I_ERI_Dxz_Px_D2x_S_b;
  Double I_ERI_D2y_D2x_D2x_S_b = I_ERI_Fx2y_Px_D2x_S_b+ABX*I_ERI_D2y_Px_D2x_S_b;
  Double I_ERI_Dyz_D2x_D2x_S_b = I_ERI_Fxyz_Px_D2x_S_b+ABX*I_ERI_Dyz_Px_D2x_S_b;
  Double I_ERI_D2z_D2x_D2x_S_b = I_ERI_Fx2z_Px_D2x_S_b+ABX*I_ERI_D2z_Px_D2x_S_b;
  Double I_ERI_D2x_Dxy_D2x_S_b = I_ERI_F2xy_Px_D2x_S_b+ABY*I_ERI_D2x_Px_D2x_S_b;
  Double I_ERI_Dxy_Dxy_D2x_S_b = I_ERI_Fx2y_Px_D2x_S_b+ABY*I_ERI_Dxy_Px_D2x_S_b;
  Double I_ERI_Dxz_Dxy_D2x_S_b = I_ERI_Fxyz_Px_D2x_S_b+ABY*I_ERI_Dxz_Px_D2x_S_b;
  Double I_ERI_D2y_Dxy_D2x_S_b = I_ERI_F3y_Px_D2x_S_b+ABY*I_ERI_D2y_Px_D2x_S_b;
  Double I_ERI_Dyz_Dxy_D2x_S_b = I_ERI_F2yz_Px_D2x_S_b+ABY*I_ERI_Dyz_Px_D2x_S_b;
  Double I_ERI_D2z_Dxy_D2x_S_b = I_ERI_Fy2z_Px_D2x_S_b+ABY*I_ERI_D2z_Px_D2x_S_b;
  Double I_ERI_D2x_D2y_D2x_S_b = I_ERI_F2xy_Py_D2x_S_b+ABY*I_ERI_D2x_Py_D2x_S_b;
  Double I_ERI_Dxy_D2y_D2x_S_b = I_ERI_Fx2y_Py_D2x_S_b+ABY*I_ERI_Dxy_Py_D2x_S_b;
  Double I_ERI_Dxz_D2y_D2x_S_b = I_ERI_Fxyz_Py_D2x_S_b+ABY*I_ERI_Dxz_Py_D2x_S_b;
  Double I_ERI_D2y_D2y_D2x_S_b = I_ERI_F3y_Py_D2x_S_b+ABY*I_ERI_D2y_Py_D2x_S_b;
  Double I_ERI_Dyz_D2y_D2x_S_b = I_ERI_F2yz_Py_D2x_S_b+ABY*I_ERI_Dyz_Py_D2x_S_b;
  Double I_ERI_D2z_D2y_D2x_S_b = I_ERI_Fy2z_Py_D2x_S_b+ABY*I_ERI_D2z_Py_D2x_S_b;
  Double I_ERI_D2x_D2z_D2x_S_b = I_ERI_F2xz_Pz_D2x_S_b+ABZ*I_ERI_D2x_Pz_D2x_S_b;
  Double I_ERI_Dxy_D2z_D2x_S_b = I_ERI_Fxyz_Pz_D2x_S_b+ABZ*I_ERI_Dxy_Pz_D2x_S_b;
  Double I_ERI_Dxz_D2z_D2x_S_b = I_ERI_Fx2z_Pz_D2x_S_b+ABZ*I_ERI_Dxz_Pz_D2x_S_b;
  Double I_ERI_D2y_D2z_D2x_S_b = I_ERI_F2yz_Pz_D2x_S_b+ABZ*I_ERI_D2y_Pz_D2x_S_b;
  Double I_ERI_Dyz_D2z_D2x_S_b = I_ERI_Fy2z_Pz_D2x_S_b+ABZ*I_ERI_Dyz_Pz_D2x_S_b;
  Double I_ERI_D2z_D2z_D2x_S_b = I_ERI_F3z_Pz_D2x_S_b+ABZ*I_ERI_D2z_Pz_D2x_S_b;
  Double I_ERI_D2x_D2x_Dxy_S_b = I_ERI_F3x_Px_Dxy_S_b+ABX*I_ERI_D2x_Px_Dxy_S_b;
  Double I_ERI_Dxy_D2x_Dxy_S_b = I_ERI_F2xy_Px_Dxy_S_b+ABX*I_ERI_Dxy_Px_Dxy_S_b;
  Double I_ERI_Dxz_D2x_Dxy_S_b = I_ERI_F2xz_Px_Dxy_S_b+ABX*I_ERI_Dxz_Px_Dxy_S_b;
  Double I_ERI_D2y_D2x_Dxy_S_b = I_ERI_Fx2y_Px_Dxy_S_b+ABX*I_ERI_D2y_Px_Dxy_S_b;
  Double I_ERI_Dyz_D2x_Dxy_S_b = I_ERI_Fxyz_Px_Dxy_S_b+ABX*I_ERI_Dyz_Px_Dxy_S_b;
  Double I_ERI_D2z_D2x_Dxy_S_b = I_ERI_Fx2z_Px_Dxy_S_b+ABX*I_ERI_D2z_Px_Dxy_S_b;
  Double I_ERI_D2x_Dxy_Dxy_S_b = I_ERI_F2xy_Px_Dxy_S_b+ABY*I_ERI_D2x_Px_Dxy_S_b;
  Double I_ERI_Dxy_Dxy_Dxy_S_b = I_ERI_Fx2y_Px_Dxy_S_b+ABY*I_ERI_Dxy_Px_Dxy_S_b;
  Double I_ERI_Dxz_Dxy_Dxy_S_b = I_ERI_Fxyz_Px_Dxy_S_b+ABY*I_ERI_Dxz_Px_Dxy_S_b;
  Double I_ERI_D2y_Dxy_Dxy_S_b = I_ERI_F3y_Px_Dxy_S_b+ABY*I_ERI_D2y_Px_Dxy_S_b;
  Double I_ERI_Dyz_Dxy_Dxy_S_b = I_ERI_F2yz_Px_Dxy_S_b+ABY*I_ERI_Dyz_Px_Dxy_S_b;
  Double I_ERI_D2z_Dxy_Dxy_S_b = I_ERI_Fy2z_Px_Dxy_S_b+ABY*I_ERI_D2z_Px_Dxy_S_b;
  Double I_ERI_D2x_D2y_Dxy_S_b = I_ERI_F2xy_Py_Dxy_S_b+ABY*I_ERI_D2x_Py_Dxy_S_b;
  Double I_ERI_Dxy_D2y_Dxy_S_b = I_ERI_Fx2y_Py_Dxy_S_b+ABY*I_ERI_Dxy_Py_Dxy_S_b;
  Double I_ERI_Dxz_D2y_Dxy_S_b = I_ERI_Fxyz_Py_Dxy_S_b+ABY*I_ERI_Dxz_Py_Dxy_S_b;
  Double I_ERI_D2y_D2y_Dxy_S_b = I_ERI_F3y_Py_Dxy_S_b+ABY*I_ERI_D2y_Py_Dxy_S_b;
  Double I_ERI_Dyz_D2y_Dxy_S_b = I_ERI_F2yz_Py_Dxy_S_b+ABY*I_ERI_Dyz_Py_Dxy_S_b;
  Double I_ERI_D2z_D2y_Dxy_S_b = I_ERI_Fy2z_Py_Dxy_S_b+ABY*I_ERI_D2z_Py_Dxy_S_b;
  Double I_ERI_D2x_D2z_Dxy_S_b = I_ERI_F2xz_Pz_Dxy_S_b+ABZ*I_ERI_D2x_Pz_Dxy_S_b;
  Double I_ERI_Dxy_D2z_Dxy_S_b = I_ERI_Fxyz_Pz_Dxy_S_b+ABZ*I_ERI_Dxy_Pz_Dxy_S_b;
  Double I_ERI_Dxz_D2z_Dxy_S_b = I_ERI_Fx2z_Pz_Dxy_S_b+ABZ*I_ERI_Dxz_Pz_Dxy_S_b;
  Double I_ERI_D2y_D2z_Dxy_S_b = I_ERI_F2yz_Pz_Dxy_S_b+ABZ*I_ERI_D2y_Pz_Dxy_S_b;
  Double I_ERI_Dyz_D2z_Dxy_S_b = I_ERI_Fy2z_Pz_Dxy_S_b+ABZ*I_ERI_Dyz_Pz_Dxy_S_b;
  Double I_ERI_D2z_D2z_Dxy_S_b = I_ERI_F3z_Pz_Dxy_S_b+ABZ*I_ERI_D2z_Pz_Dxy_S_b;
  Double I_ERI_D2x_D2x_Dxz_S_b = I_ERI_F3x_Px_Dxz_S_b+ABX*I_ERI_D2x_Px_Dxz_S_b;
  Double I_ERI_Dxy_D2x_Dxz_S_b = I_ERI_F2xy_Px_Dxz_S_b+ABX*I_ERI_Dxy_Px_Dxz_S_b;
  Double I_ERI_Dxz_D2x_Dxz_S_b = I_ERI_F2xz_Px_Dxz_S_b+ABX*I_ERI_Dxz_Px_Dxz_S_b;
  Double I_ERI_D2y_D2x_Dxz_S_b = I_ERI_Fx2y_Px_Dxz_S_b+ABX*I_ERI_D2y_Px_Dxz_S_b;
  Double I_ERI_Dyz_D2x_Dxz_S_b = I_ERI_Fxyz_Px_Dxz_S_b+ABX*I_ERI_Dyz_Px_Dxz_S_b;
  Double I_ERI_D2z_D2x_Dxz_S_b = I_ERI_Fx2z_Px_Dxz_S_b+ABX*I_ERI_D2z_Px_Dxz_S_b;
  Double I_ERI_D2x_Dxy_Dxz_S_b = I_ERI_F2xy_Px_Dxz_S_b+ABY*I_ERI_D2x_Px_Dxz_S_b;
  Double I_ERI_Dxy_Dxy_Dxz_S_b = I_ERI_Fx2y_Px_Dxz_S_b+ABY*I_ERI_Dxy_Px_Dxz_S_b;
  Double I_ERI_Dxz_Dxy_Dxz_S_b = I_ERI_Fxyz_Px_Dxz_S_b+ABY*I_ERI_Dxz_Px_Dxz_S_b;
  Double I_ERI_D2y_Dxy_Dxz_S_b = I_ERI_F3y_Px_Dxz_S_b+ABY*I_ERI_D2y_Px_Dxz_S_b;
  Double I_ERI_Dyz_Dxy_Dxz_S_b = I_ERI_F2yz_Px_Dxz_S_b+ABY*I_ERI_Dyz_Px_Dxz_S_b;
  Double I_ERI_D2z_Dxy_Dxz_S_b = I_ERI_Fy2z_Px_Dxz_S_b+ABY*I_ERI_D2z_Px_Dxz_S_b;
  Double I_ERI_D2x_D2y_Dxz_S_b = I_ERI_F2xy_Py_Dxz_S_b+ABY*I_ERI_D2x_Py_Dxz_S_b;
  Double I_ERI_Dxy_D2y_Dxz_S_b = I_ERI_Fx2y_Py_Dxz_S_b+ABY*I_ERI_Dxy_Py_Dxz_S_b;
  Double I_ERI_Dxz_D2y_Dxz_S_b = I_ERI_Fxyz_Py_Dxz_S_b+ABY*I_ERI_Dxz_Py_Dxz_S_b;
  Double I_ERI_D2y_D2y_Dxz_S_b = I_ERI_F3y_Py_Dxz_S_b+ABY*I_ERI_D2y_Py_Dxz_S_b;
  Double I_ERI_Dyz_D2y_Dxz_S_b = I_ERI_F2yz_Py_Dxz_S_b+ABY*I_ERI_Dyz_Py_Dxz_S_b;
  Double I_ERI_D2z_D2y_Dxz_S_b = I_ERI_Fy2z_Py_Dxz_S_b+ABY*I_ERI_D2z_Py_Dxz_S_b;
  Double I_ERI_D2x_D2z_Dxz_S_b = I_ERI_F2xz_Pz_Dxz_S_b+ABZ*I_ERI_D2x_Pz_Dxz_S_b;
  Double I_ERI_Dxy_D2z_Dxz_S_b = I_ERI_Fxyz_Pz_Dxz_S_b+ABZ*I_ERI_Dxy_Pz_Dxz_S_b;
  Double I_ERI_Dxz_D2z_Dxz_S_b = I_ERI_Fx2z_Pz_Dxz_S_b+ABZ*I_ERI_Dxz_Pz_Dxz_S_b;
  Double I_ERI_D2y_D2z_Dxz_S_b = I_ERI_F2yz_Pz_Dxz_S_b+ABZ*I_ERI_D2y_Pz_Dxz_S_b;
  Double I_ERI_Dyz_D2z_Dxz_S_b = I_ERI_Fy2z_Pz_Dxz_S_b+ABZ*I_ERI_Dyz_Pz_Dxz_S_b;
  Double I_ERI_D2z_D2z_Dxz_S_b = I_ERI_F3z_Pz_Dxz_S_b+ABZ*I_ERI_D2z_Pz_Dxz_S_b;
  Double I_ERI_D2x_D2x_D2y_S_b = I_ERI_F3x_Px_D2y_S_b+ABX*I_ERI_D2x_Px_D2y_S_b;
  Double I_ERI_Dxy_D2x_D2y_S_b = I_ERI_F2xy_Px_D2y_S_b+ABX*I_ERI_Dxy_Px_D2y_S_b;
  Double I_ERI_Dxz_D2x_D2y_S_b = I_ERI_F2xz_Px_D2y_S_b+ABX*I_ERI_Dxz_Px_D2y_S_b;
  Double I_ERI_D2y_D2x_D2y_S_b = I_ERI_Fx2y_Px_D2y_S_b+ABX*I_ERI_D2y_Px_D2y_S_b;
  Double I_ERI_Dyz_D2x_D2y_S_b = I_ERI_Fxyz_Px_D2y_S_b+ABX*I_ERI_Dyz_Px_D2y_S_b;
  Double I_ERI_D2z_D2x_D2y_S_b = I_ERI_Fx2z_Px_D2y_S_b+ABX*I_ERI_D2z_Px_D2y_S_b;
  Double I_ERI_D2x_Dxy_D2y_S_b = I_ERI_F2xy_Px_D2y_S_b+ABY*I_ERI_D2x_Px_D2y_S_b;
  Double I_ERI_Dxy_Dxy_D2y_S_b = I_ERI_Fx2y_Px_D2y_S_b+ABY*I_ERI_Dxy_Px_D2y_S_b;
  Double I_ERI_Dxz_Dxy_D2y_S_b = I_ERI_Fxyz_Px_D2y_S_b+ABY*I_ERI_Dxz_Px_D2y_S_b;
  Double I_ERI_D2y_Dxy_D2y_S_b = I_ERI_F3y_Px_D2y_S_b+ABY*I_ERI_D2y_Px_D2y_S_b;
  Double I_ERI_Dyz_Dxy_D2y_S_b = I_ERI_F2yz_Px_D2y_S_b+ABY*I_ERI_Dyz_Px_D2y_S_b;
  Double I_ERI_D2z_Dxy_D2y_S_b = I_ERI_Fy2z_Px_D2y_S_b+ABY*I_ERI_D2z_Px_D2y_S_b;
  Double I_ERI_D2x_D2y_D2y_S_b = I_ERI_F2xy_Py_D2y_S_b+ABY*I_ERI_D2x_Py_D2y_S_b;
  Double I_ERI_Dxy_D2y_D2y_S_b = I_ERI_Fx2y_Py_D2y_S_b+ABY*I_ERI_Dxy_Py_D2y_S_b;
  Double I_ERI_Dxz_D2y_D2y_S_b = I_ERI_Fxyz_Py_D2y_S_b+ABY*I_ERI_Dxz_Py_D2y_S_b;
  Double I_ERI_D2y_D2y_D2y_S_b = I_ERI_F3y_Py_D2y_S_b+ABY*I_ERI_D2y_Py_D2y_S_b;
  Double I_ERI_Dyz_D2y_D2y_S_b = I_ERI_F2yz_Py_D2y_S_b+ABY*I_ERI_Dyz_Py_D2y_S_b;
  Double I_ERI_D2z_D2y_D2y_S_b = I_ERI_Fy2z_Py_D2y_S_b+ABY*I_ERI_D2z_Py_D2y_S_b;
  Double I_ERI_D2x_D2z_D2y_S_b = I_ERI_F2xz_Pz_D2y_S_b+ABZ*I_ERI_D2x_Pz_D2y_S_b;
  Double I_ERI_Dxy_D2z_D2y_S_b = I_ERI_Fxyz_Pz_D2y_S_b+ABZ*I_ERI_Dxy_Pz_D2y_S_b;
  Double I_ERI_Dxz_D2z_D2y_S_b = I_ERI_Fx2z_Pz_D2y_S_b+ABZ*I_ERI_Dxz_Pz_D2y_S_b;
  Double I_ERI_D2y_D2z_D2y_S_b = I_ERI_F2yz_Pz_D2y_S_b+ABZ*I_ERI_D2y_Pz_D2y_S_b;
  Double I_ERI_Dyz_D2z_D2y_S_b = I_ERI_Fy2z_Pz_D2y_S_b+ABZ*I_ERI_Dyz_Pz_D2y_S_b;
  Double I_ERI_D2z_D2z_D2y_S_b = I_ERI_F3z_Pz_D2y_S_b+ABZ*I_ERI_D2z_Pz_D2y_S_b;
  Double I_ERI_D2x_D2x_Dyz_S_b = I_ERI_F3x_Px_Dyz_S_b+ABX*I_ERI_D2x_Px_Dyz_S_b;
  Double I_ERI_Dxy_D2x_Dyz_S_b = I_ERI_F2xy_Px_Dyz_S_b+ABX*I_ERI_Dxy_Px_Dyz_S_b;
  Double I_ERI_Dxz_D2x_Dyz_S_b = I_ERI_F2xz_Px_Dyz_S_b+ABX*I_ERI_Dxz_Px_Dyz_S_b;
  Double I_ERI_D2y_D2x_Dyz_S_b = I_ERI_Fx2y_Px_Dyz_S_b+ABX*I_ERI_D2y_Px_Dyz_S_b;
  Double I_ERI_Dyz_D2x_Dyz_S_b = I_ERI_Fxyz_Px_Dyz_S_b+ABX*I_ERI_Dyz_Px_Dyz_S_b;
  Double I_ERI_D2z_D2x_Dyz_S_b = I_ERI_Fx2z_Px_Dyz_S_b+ABX*I_ERI_D2z_Px_Dyz_S_b;
  Double I_ERI_D2x_Dxy_Dyz_S_b = I_ERI_F2xy_Px_Dyz_S_b+ABY*I_ERI_D2x_Px_Dyz_S_b;
  Double I_ERI_Dxy_Dxy_Dyz_S_b = I_ERI_Fx2y_Px_Dyz_S_b+ABY*I_ERI_Dxy_Px_Dyz_S_b;
  Double I_ERI_Dxz_Dxy_Dyz_S_b = I_ERI_Fxyz_Px_Dyz_S_b+ABY*I_ERI_Dxz_Px_Dyz_S_b;
  Double I_ERI_D2y_Dxy_Dyz_S_b = I_ERI_F3y_Px_Dyz_S_b+ABY*I_ERI_D2y_Px_Dyz_S_b;
  Double I_ERI_Dyz_Dxy_Dyz_S_b = I_ERI_F2yz_Px_Dyz_S_b+ABY*I_ERI_Dyz_Px_Dyz_S_b;
  Double I_ERI_D2z_Dxy_Dyz_S_b = I_ERI_Fy2z_Px_Dyz_S_b+ABY*I_ERI_D2z_Px_Dyz_S_b;
  Double I_ERI_D2x_D2y_Dyz_S_b = I_ERI_F2xy_Py_Dyz_S_b+ABY*I_ERI_D2x_Py_Dyz_S_b;
  Double I_ERI_Dxy_D2y_Dyz_S_b = I_ERI_Fx2y_Py_Dyz_S_b+ABY*I_ERI_Dxy_Py_Dyz_S_b;
  Double I_ERI_Dxz_D2y_Dyz_S_b = I_ERI_Fxyz_Py_Dyz_S_b+ABY*I_ERI_Dxz_Py_Dyz_S_b;
  Double I_ERI_D2y_D2y_Dyz_S_b = I_ERI_F3y_Py_Dyz_S_b+ABY*I_ERI_D2y_Py_Dyz_S_b;
  Double I_ERI_Dyz_D2y_Dyz_S_b = I_ERI_F2yz_Py_Dyz_S_b+ABY*I_ERI_Dyz_Py_Dyz_S_b;
  Double I_ERI_D2z_D2y_Dyz_S_b = I_ERI_Fy2z_Py_Dyz_S_b+ABY*I_ERI_D2z_Py_Dyz_S_b;
  Double I_ERI_D2x_D2z_Dyz_S_b = I_ERI_F2xz_Pz_Dyz_S_b+ABZ*I_ERI_D2x_Pz_Dyz_S_b;
  Double I_ERI_Dxy_D2z_Dyz_S_b = I_ERI_Fxyz_Pz_Dyz_S_b+ABZ*I_ERI_Dxy_Pz_Dyz_S_b;
  Double I_ERI_Dxz_D2z_Dyz_S_b = I_ERI_Fx2z_Pz_Dyz_S_b+ABZ*I_ERI_Dxz_Pz_Dyz_S_b;
  Double I_ERI_D2y_D2z_Dyz_S_b = I_ERI_F2yz_Pz_Dyz_S_b+ABZ*I_ERI_D2y_Pz_Dyz_S_b;
  Double I_ERI_Dyz_D2z_Dyz_S_b = I_ERI_Fy2z_Pz_Dyz_S_b+ABZ*I_ERI_Dyz_Pz_Dyz_S_b;
  Double I_ERI_D2z_D2z_Dyz_S_b = I_ERI_F3z_Pz_Dyz_S_b+ABZ*I_ERI_D2z_Pz_Dyz_S_b;
  Double I_ERI_D2x_D2x_D2z_S_b = I_ERI_F3x_Px_D2z_S_b+ABX*I_ERI_D2x_Px_D2z_S_b;
  Double I_ERI_Dxy_D2x_D2z_S_b = I_ERI_F2xy_Px_D2z_S_b+ABX*I_ERI_Dxy_Px_D2z_S_b;
  Double I_ERI_Dxz_D2x_D2z_S_b = I_ERI_F2xz_Px_D2z_S_b+ABX*I_ERI_Dxz_Px_D2z_S_b;
  Double I_ERI_D2y_D2x_D2z_S_b = I_ERI_Fx2y_Px_D2z_S_b+ABX*I_ERI_D2y_Px_D2z_S_b;
  Double I_ERI_Dyz_D2x_D2z_S_b = I_ERI_Fxyz_Px_D2z_S_b+ABX*I_ERI_Dyz_Px_D2z_S_b;
  Double I_ERI_D2z_D2x_D2z_S_b = I_ERI_Fx2z_Px_D2z_S_b+ABX*I_ERI_D2z_Px_D2z_S_b;
  Double I_ERI_D2x_Dxy_D2z_S_b = I_ERI_F2xy_Px_D2z_S_b+ABY*I_ERI_D2x_Px_D2z_S_b;
  Double I_ERI_Dxy_Dxy_D2z_S_b = I_ERI_Fx2y_Px_D2z_S_b+ABY*I_ERI_Dxy_Px_D2z_S_b;
  Double I_ERI_Dxz_Dxy_D2z_S_b = I_ERI_Fxyz_Px_D2z_S_b+ABY*I_ERI_Dxz_Px_D2z_S_b;
  Double I_ERI_D2y_Dxy_D2z_S_b = I_ERI_F3y_Px_D2z_S_b+ABY*I_ERI_D2y_Px_D2z_S_b;
  Double I_ERI_Dyz_Dxy_D2z_S_b = I_ERI_F2yz_Px_D2z_S_b+ABY*I_ERI_Dyz_Px_D2z_S_b;
  Double I_ERI_D2z_Dxy_D2z_S_b = I_ERI_Fy2z_Px_D2z_S_b+ABY*I_ERI_D2z_Px_D2z_S_b;
  Double I_ERI_D2x_D2y_D2z_S_b = I_ERI_F2xy_Py_D2z_S_b+ABY*I_ERI_D2x_Py_D2z_S_b;
  Double I_ERI_Dxy_D2y_D2z_S_b = I_ERI_Fx2y_Py_D2z_S_b+ABY*I_ERI_Dxy_Py_D2z_S_b;
  Double I_ERI_Dxz_D2y_D2z_S_b = I_ERI_Fxyz_Py_D2z_S_b+ABY*I_ERI_Dxz_Py_D2z_S_b;
  Double I_ERI_D2y_D2y_D2z_S_b = I_ERI_F3y_Py_D2z_S_b+ABY*I_ERI_D2y_Py_D2z_S_b;
  Double I_ERI_Dyz_D2y_D2z_S_b = I_ERI_F2yz_Py_D2z_S_b+ABY*I_ERI_Dyz_Py_D2z_S_b;
  Double I_ERI_D2z_D2y_D2z_S_b = I_ERI_Fy2z_Py_D2z_S_b+ABY*I_ERI_D2z_Py_D2z_S_b;
  Double I_ERI_D2x_D2z_D2z_S_b = I_ERI_F2xz_Pz_D2z_S_b+ABZ*I_ERI_D2x_Pz_D2z_S_b;
  Double I_ERI_Dxy_D2z_D2z_S_b = I_ERI_Fxyz_Pz_D2z_S_b+ABZ*I_ERI_Dxy_Pz_D2z_S_b;
  Double I_ERI_Dxz_D2z_D2z_S_b = I_ERI_Fx2z_Pz_D2z_S_b+ABZ*I_ERI_Dxz_Pz_D2z_S_b;
  Double I_ERI_D2y_D2z_D2z_S_b = I_ERI_F2yz_Pz_D2z_S_b+ABZ*I_ERI_D2y_Pz_D2z_S_b;
  Double I_ERI_Dyz_D2z_D2z_S_b = I_ERI_Fy2z_Pz_D2z_S_b+ABZ*I_ERI_Dyz_Pz_D2z_S_b;
  Double I_ERI_D2z_D2z_D2z_S_b = I_ERI_F3z_Pz_D2z_S_b+ABZ*I_ERI_D2z_Pz_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 72 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_D_S_b
   * RHS shell quartet name: SQ_ERI_G_S_D_S_b
   ************************************************************/
  Double I_ERI_G4x_Px_D2x_S_b = I_ERI_H5x_S_D2x_S_b+ABX*I_ERI_G4x_S_D2x_S_b;
  Double I_ERI_G3xy_Px_D2x_S_b = I_ERI_H4xy_S_D2x_S_b+ABX*I_ERI_G3xy_S_D2x_S_b;
  Double I_ERI_G3xz_Px_D2x_S_b = I_ERI_H4xz_S_D2x_S_b+ABX*I_ERI_G3xz_S_D2x_S_b;
  Double I_ERI_G2x2y_Px_D2x_S_b = I_ERI_H3x2y_S_D2x_S_b+ABX*I_ERI_G2x2y_S_D2x_S_b;
  Double I_ERI_G2xyz_Px_D2x_S_b = I_ERI_H3xyz_S_D2x_S_b+ABX*I_ERI_G2xyz_S_D2x_S_b;
  Double I_ERI_G2x2z_Px_D2x_S_b = I_ERI_H3x2z_S_D2x_S_b+ABX*I_ERI_G2x2z_S_D2x_S_b;
  Double I_ERI_Gx3y_Px_D2x_S_b = I_ERI_H2x3y_S_D2x_S_b+ABX*I_ERI_Gx3y_S_D2x_S_b;
  Double I_ERI_Gx2yz_Px_D2x_S_b = I_ERI_H2x2yz_S_D2x_S_b+ABX*I_ERI_Gx2yz_S_D2x_S_b;
  Double I_ERI_Gxy2z_Px_D2x_S_b = I_ERI_H2xy2z_S_D2x_S_b+ABX*I_ERI_Gxy2z_S_D2x_S_b;
  Double I_ERI_Gx3z_Px_D2x_S_b = I_ERI_H2x3z_S_D2x_S_b+ABX*I_ERI_Gx3z_S_D2x_S_b;
  Double I_ERI_G3yz_Px_D2x_S_b = I_ERI_Hx3yz_S_D2x_S_b+ABX*I_ERI_G3yz_S_D2x_S_b;
  Double I_ERI_G2y2z_Px_D2x_S_b = I_ERI_Hx2y2z_S_D2x_S_b+ABX*I_ERI_G2y2z_S_D2x_S_b;
  Double I_ERI_Gy3z_Px_D2x_S_b = I_ERI_Hxy3z_S_D2x_S_b+ABX*I_ERI_Gy3z_S_D2x_S_b;
  Double I_ERI_G3xy_Py_D2x_S_b = I_ERI_H3x2y_S_D2x_S_b+ABY*I_ERI_G3xy_S_D2x_S_b;
  Double I_ERI_G2x2y_Py_D2x_S_b = I_ERI_H2x3y_S_D2x_S_b+ABY*I_ERI_G2x2y_S_D2x_S_b;
  Double I_ERI_G2xyz_Py_D2x_S_b = I_ERI_H2x2yz_S_D2x_S_b+ABY*I_ERI_G2xyz_S_D2x_S_b;
  Double I_ERI_Gx3y_Py_D2x_S_b = I_ERI_Hx4y_S_D2x_S_b+ABY*I_ERI_Gx3y_S_D2x_S_b;
  Double I_ERI_Gx2yz_Py_D2x_S_b = I_ERI_Hx3yz_S_D2x_S_b+ABY*I_ERI_Gx2yz_S_D2x_S_b;
  Double I_ERI_Gxy2z_Py_D2x_S_b = I_ERI_Hx2y2z_S_D2x_S_b+ABY*I_ERI_Gxy2z_S_D2x_S_b;
  Double I_ERI_G4y_Py_D2x_S_b = I_ERI_H5y_S_D2x_S_b+ABY*I_ERI_G4y_S_D2x_S_b;
  Double I_ERI_G3yz_Py_D2x_S_b = I_ERI_H4yz_S_D2x_S_b+ABY*I_ERI_G3yz_S_D2x_S_b;
  Double I_ERI_G2y2z_Py_D2x_S_b = I_ERI_H3y2z_S_D2x_S_b+ABY*I_ERI_G2y2z_S_D2x_S_b;
  Double I_ERI_Gy3z_Py_D2x_S_b = I_ERI_H2y3z_S_D2x_S_b+ABY*I_ERI_Gy3z_S_D2x_S_b;
  Double I_ERI_G3xz_Pz_D2x_S_b = I_ERI_H3x2z_S_D2x_S_b+ABZ*I_ERI_G3xz_S_D2x_S_b;
  Double I_ERI_G2xyz_Pz_D2x_S_b = I_ERI_H2xy2z_S_D2x_S_b+ABZ*I_ERI_G2xyz_S_D2x_S_b;
  Double I_ERI_G2x2z_Pz_D2x_S_b = I_ERI_H2x3z_S_D2x_S_b+ABZ*I_ERI_G2x2z_S_D2x_S_b;
  Double I_ERI_Gx2yz_Pz_D2x_S_b = I_ERI_Hx2y2z_S_D2x_S_b+ABZ*I_ERI_Gx2yz_S_D2x_S_b;
  Double I_ERI_Gxy2z_Pz_D2x_S_b = I_ERI_Hxy3z_S_D2x_S_b+ABZ*I_ERI_Gxy2z_S_D2x_S_b;
  Double I_ERI_Gx3z_Pz_D2x_S_b = I_ERI_Hx4z_S_D2x_S_b+ABZ*I_ERI_Gx3z_S_D2x_S_b;
  Double I_ERI_G3yz_Pz_D2x_S_b = I_ERI_H3y2z_S_D2x_S_b+ABZ*I_ERI_G3yz_S_D2x_S_b;
  Double I_ERI_G2y2z_Pz_D2x_S_b = I_ERI_H2y3z_S_D2x_S_b+ABZ*I_ERI_G2y2z_S_D2x_S_b;
  Double I_ERI_Gy3z_Pz_D2x_S_b = I_ERI_Hy4z_S_D2x_S_b+ABZ*I_ERI_Gy3z_S_D2x_S_b;
  Double I_ERI_G4z_Pz_D2x_S_b = I_ERI_H5z_S_D2x_S_b+ABZ*I_ERI_G4z_S_D2x_S_b;
  Double I_ERI_G4x_Px_Dxy_S_b = I_ERI_H5x_S_Dxy_S_b+ABX*I_ERI_G4x_S_Dxy_S_b;
  Double I_ERI_G3xy_Px_Dxy_S_b = I_ERI_H4xy_S_Dxy_S_b+ABX*I_ERI_G3xy_S_Dxy_S_b;
  Double I_ERI_G3xz_Px_Dxy_S_b = I_ERI_H4xz_S_Dxy_S_b+ABX*I_ERI_G3xz_S_Dxy_S_b;
  Double I_ERI_G2x2y_Px_Dxy_S_b = I_ERI_H3x2y_S_Dxy_S_b+ABX*I_ERI_G2x2y_S_Dxy_S_b;
  Double I_ERI_G2xyz_Px_Dxy_S_b = I_ERI_H3xyz_S_Dxy_S_b+ABX*I_ERI_G2xyz_S_Dxy_S_b;
  Double I_ERI_G2x2z_Px_Dxy_S_b = I_ERI_H3x2z_S_Dxy_S_b+ABX*I_ERI_G2x2z_S_Dxy_S_b;
  Double I_ERI_Gx3y_Px_Dxy_S_b = I_ERI_H2x3y_S_Dxy_S_b+ABX*I_ERI_Gx3y_S_Dxy_S_b;
  Double I_ERI_Gx2yz_Px_Dxy_S_b = I_ERI_H2x2yz_S_Dxy_S_b+ABX*I_ERI_Gx2yz_S_Dxy_S_b;
  Double I_ERI_Gxy2z_Px_Dxy_S_b = I_ERI_H2xy2z_S_Dxy_S_b+ABX*I_ERI_Gxy2z_S_Dxy_S_b;
  Double I_ERI_Gx3z_Px_Dxy_S_b = I_ERI_H2x3z_S_Dxy_S_b+ABX*I_ERI_Gx3z_S_Dxy_S_b;
  Double I_ERI_G3yz_Px_Dxy_S_b = I_ERI_Hx3yz_S_Dxy_S_b+ABX*I_ERI_G3yz_S_Dxy_S_b;
  Double I_ERI_G2y2z_Px_Dxy_S_b = I_ERI_Hx2y2z_S_Dxy_S_b+ABX*I_ERI_G2y2z_S_Dxy_S_b;
  Double I_ERI_Gy3z_Px_Dxy_S_b = I_ERI_Hxy3z_S_Dxy_S_b+ABX*I_ERI_Gy3z_S_Dxy_S_b;
  Double I_ERI_G3xy_Py_Dxy_S_b = I_ERI_H3x2y_S_Dxy_S_b+ABY*I_ERI_G3xy_S_Dxy_S_b;
  Double I_ERI_G2x2y_Py_Dxy_S_b = I_ERI_H2x3y_S_Dxy_S_b+ABY*I_ERI_G2x2y_S_Dxy_S_b;
  Double I_ERI_G2xyz_Py_Dxy_S_b = I_ERI_H2x2yz_S_Dxy_S_b+ABY*I_ERI_G2xyz_S_Dxy_S_b;
  Double I_ERI_Gx3y_Py_Dxy_S_b = I_ERI_Hx4y_S_Dxy_S_b+ABY*I_ERI_Gx3y_S_Dxy_S_b;
  Double I_ERI_Gx2yz_Py_Dxy_S_b = I_ERI_Hx3yz_S_Dxy_S_b+ABY*I_ERI_Gx2yz_S_Dxy_S_b;
  Double I_ERI_Gxy2z_Py_Dxy_S_b = I_ERI_Hx2y2z_S_Dxy_S_b+ABY*I_ERI_Gxy2z_S_Dxy_S_b;
  Double I_ERI_G4y_Py_Dxy_S_b = I_ERI_H5y_S_Dxy_S_b+ABY*I_ERI_G4y_S_Dxy_S_b;
  Double I_ERI_G3yz_Py_Dxy_S_b = I_ERI_H4yz_S_Dxy_S_b+ABY*I_ERI_G3yz_S_Dxy_S_b;
  Double I_ERI_G2y2z_Py_Dxy_S_b = I_ERI_H3y2z_S_Dxy_S_b+ABY*I_ERI_G2y2z_S_Dxy_S_b;
  Double I_ERI_Gy3z_Py_Dxy_S_b = I_ERI_H2y3z_S_Dxy_S_b+ABY*I_ERI_Gy3z_S_Dxy_S_b;
  Double I_ERI_G3xz_Pz_Dxy_S_b = I_ERI_H3x2z_S_Dxy_S_b+ABZ*I_ERI_G3xz_S_Dxy_S_b;
  Double I_ERI_G2xyz_Pz_Dxy_S_b = I_ERI_H2xy2z_S_Dxy_S_b+ABZ*I_ERI_G2xyz_S_Dxy_S_b;
  Double I_ERI_G2x2z_Pz_Dxy_S_b = I_ERI_H2x3z_S_Dxy_S_b+ABZ*I_ERI_G2x2z_S_Dxy_S_b;
  Double I_ERI_Gx2yz_Pz_Dxy_S_b = I_ERI_Hx2y2z_S_Dxy_S_b+ABZ*I_ERI_Gx2yz_S_Dxy_S_b;
  Double I_ERI_Gxy2z_Pz_Dxy_S_b = I_ERI_Hxy3z_S_Dxy_S_b+ABZ*I_ERI_Gxy2z_S_Dxy_S_b;
  Double I_ERI_Gx3z_Pz_Dxy_S_b = I_ERI_Hx4z_S_Dxy_S_b+ABZ*I_ERI_Gx3z_S_Dxy_S_b;
  Double I_ERI_G3yz_Pz_Dxy_S_b = I_ERI_H3y2z_S_Dxy_S_b+ABZ*I_ERI_G3yz_S_Dxy_S_b;
  Double I_ERI_G2y2z_Pz_Dxy_S_b = I_ERI_H2y3z_S_Dxy_S_b+ABZ*I_ERI_G2y2z_S_Dxy_S_b;
  Double I_ERI_Gy3z_Pz_Dxy_S_b = I_ERI_Hy4z_S_Dxy_S_b+ABZ*I_ERI_Gy3z_S_Dxy_S_b;
  Double I_ERI_G4z_Pz_Dxy_S_b = I_ERI_H5z_S_Dxy_S_b+ABZ*I_ERI_G4z_S_Dxy_S_b;
  Double I_ERI_G4x_Px_Dxz_S_b = I_ERI_H5x_S_Dxz_S_b+ABX*I_ERI_G4x_S_Dxz_S_b;
  Double I_ERI_G3xy_Px_Dxz_S_b = I_ERI_H4xy_S_Dxz_S_b+ABX*I_ERI_G3xy_S_Dxz_S_b;
  Double I_ERI_G3xz_Px_Dxz_S_b = I_ERI_H4xz_S_Dxz_S_b+ABX*I_ERI_G3xz_S_Dxz_S_b;
  Double I_ERI_G2x2y_Px_Dxz_S_b = I_ERI_H3x2y_S_Dxz_S_b+ABX*I_ERI_G2x2y_S_Dxz_S_b;
  Double I_ERI_G2xyz_Px_Dxz_S_b = I_ERI_H3xyz_S_Dxz_S_b+ABX*I_ERI_G2xyz_S_Dxz_S_b;
  Double I_ERI_G2x2z_Px_Dxz_S_b = I_ERI_H3x2z_S_Dxz_S_b+ABX*I_ERI_G2x2z_S_Dxz_S_b;
  Double I_ERI_Gx3y_Px_Dxz_S_b = I_ERI_H2x3y_S_Dxz_S_b+ABX*I_ERI_Gx3y_S_Dxz_S_b;
  Double I_ERI_Gx2yz_Px_Dxz_S_b = I_ERI_H2x2yz_S_Dxz_S_b+ABX*I_ERI_Gx2yz_S_Dxz_S_b;
  Double I_ERI_Gxy2z_Px_Dxz_S_b = I_ERI_H2xy2z_S_Dxz_S_b+ABX*I_ERI_Gxy2z_S_Dxz_S_b;
  Double I_ERI_Gx3z_Px_Dxz_S_b = I_ERI_H2x3z_S_Dxz_S_b+ABX*I_ERI_Gx3z_S_Dxz_S_b;
  Double I_ERI_G3yz_Px_Dxz_S_b = I_ERI_Hx3yz_S_Dxz_S_b+ABX*I_ERI_G3yz_S_Dxz_S_b;
  Double I_ERI_G2y2z_Px_Dxz_S_b = I_ERI_Hx2y2z_S_Dxz_S_b+ABX*I_ERI_G2y2z_S_Dxz_S_b;
  Double I_ERI_Gy3z_Px_Dxz_S_b = I_ERI_Hxy3z_S_Dxz_S_b+ABX*I_ERI_Gy3z_S_Dxz_S_b;
  Double I_ERI_G3xy_Py_Dxz_S_b = I_ERI_H3x2y_S_Dxz_S_b+ABY*I_ERI_G3xy_S_Dxz_S_b;
  Double I_ERI_G2x2y_Py_Dxz_S_b = I_ERI_H2x3y_S_Dxz_S_b+ABY*I_ERI_G2x2y_S_Dxz_S_b;
  Double I_ERI_G2xyz_Py_Dxz_S_b = I_ERI_H2x2yz_S_Dxz_S_b+ABY*I_ERI_G2xyz_S_Dxz_S_b;
  Double I_ERI_Gx3y_Py_Dxz_S_b = I_ERI_Hx4y_S_Dxz_S_b+ABY*I_ERI_Gx3y_S_Dxz_S_b;
  Double I_ERI_Gx2yz_Py_Dxz_S_b = I_ERI_Hx3yz_S_Dxz_S_b+ABY*I_ERI_Gx2yz_S_Dxz_S_b;
  Double I_ERI_Gxy2z_Py_Dxz_S_b = I_ERI_Hx2y2z_S_Dxz_S_b+ABY*I_ERI_Gxy2z_S_Dxz_S_b;
  Double I_ERI_G4y_Py_Dxz_S_b = I_ERI_H5y_S_Dxz_S_b+ABY*I_ERI_G4y_S_Dxz_S_b;
  Double I_ERI_G3yz_Py_Dxz_S_b = I_ERI_H4yz_S_Dxz_S_b+ABY*I_ERI_G3yz_S_Dxz_S_b;
  Double I_ERI_G2y2z_Py_Dxz_S_b = I_ERI_H3y2z_S_Dxz_S_b+ABY*I_ERI_G2y2z_S_Dxz_S_b;
  Double I_ERI_Gy3z_Py_Dxz_S_b = I_ERI_H2y3z_S_Dxz_S_b+ABY*I_ERI_Gy3z_S_Dxz_S_b;
  Double I_ERI_G3xz_Pz_Dxz_S_b = I_ERI_H3x2z_S_Dxz_S_b+ABZ*I_ERI_G3xz_S_Dxz_S_b;
  Double I_ERI_G2xyz_Pz_Dxz_S_b = I_ERI_H2xy2z_S_Dxz_S_b+ABZ*I_ERI_G2xyz_S_Dxz_S_b;
  Double I_ERI_G2x2z_Pz_Dxz_S_b = I_ERI_H2x3z_S_Dxz_S_b+ABZ*I_ERI_G2x2z_S_Dxz_S_b;
  Double I_ERI_Gx2yz_Pz_Dxz_S_b = I_ERI_Hx2y2z_S_Dxz_S_b+ABZ*I_ERI_Gx2yz_S_Dxz_S_b;
  Double I_ERI_Gxy2z_Pz_Dxz_S_b = I_ERI_Hxy3z_S_Dxz_S_b+ABZ*I_ERI_Gxy2z_S_Dxz_S_b;
  Double I_ERI_Gx3z_Pz_Dxz_S_b = I_ERI_Hx4z_S_Dxz_S_b+ABZ*I_ERI_Gx3z_S_Dxz_S_b;
  Double I_ERI_G3yz_Pz_Dxz_S_b = I_ERI_H3y2z_S_Dxz_S_b+ABZ*I_ERI_G3yz_S_Dxz_S_b;
  Double I_ERI_G2y2z_Pz_Dxz_S_b = I_ERI_H2y3z_S_Dxz_S_b+ABZ*I_ERI_G2y2z_S_Dxz_S_b;
  Double I_ERI_Gy3z_Pz_Dxz_S_b = I_ERI_Hy4z_S_Dxz_S_b+ABZ*I_ERI_Gy3z_S_Dxz_S_b;
  Double I_ERI_G4z_Pz_Dxz_S_b = I_ERI_H5z_S_Dxz_S_b+ABZ*I_ERI_G4z_S_Dxz_S_b;
  Double I_ERI_G4x_Px_D2y_S_b = I_ERI_H5x_S_D2y_S_b+ABX*I_ERI_G4x_S_D2y_S_b;
  Double I_ERI_G3xy_Px_D2y_S_b = I_ERI_H4xy_S_D2y_S_b+ABX*I_ERI_G3xy_S_D2y_S_b;
  Double I_ERI_G3xz_Px_D2y_S_b = I_ERI_H4xz_S_D2y_S_b+ABX*I_ERI_G3xz_S_D2y_S_b;
  Double I_ERI_G2x2y_Px_D2y_S_b = I_ERI_H3x2y_S_D2y_S_b+ABX*I_ERI_G2x2y_S_D2y_S_b;
  Double I_ERI_G2xyz_Px_D2y_S_b = I_ERI_H3xyz_S_D2y_S_b+ABX*I_ERI_G2xyz_S_D2y_S_b;
  Double I_ERI_G2x2z_Px_D2y_S_b = I_ERI_H3x2z_S_D2y_S_b+ABX*I_ERI_G2x2z_S_D2y_S_b;
  Double I_ERI_Gx3y_Px_D2y_S_b = I_ERI_H2x3y_S_D2y_S_b+ABX*I_ERI_Gx3y_S_D2y_S_b;
  Double I_ERI_Gx2yz_Px_D2y_S_b = I_ERI_H2x2yz_S_D2y_S_b+ABX*I_ERI_Gx2yz_S_D2y_S_b;
  Double I_ERI_Gxy2z_Px_D2y_S_b = I_ERI_H2xy2z_S_D2y_S_b+ABX*I_ERI_Gxy2z_S_D2y_S_b;
  Double I_ERI_Gx3z_Px_D2y_S_b = I_ERI_H2x3z_S_D2y_S_b+ABX*I_ERI_Gx3z_S_D2y_S_b;
  Double I_ERI_G3yz_Px_D2y_S_b = I_ERI_Hx3yz_S_D2y_S_b+ABX*I_ERI_G3yz_S_D2y_S_b;
  Double I_ERI_G2y2z_Px_D2y_S_b = I_ERI_Hx2y2z_S_D2y_S_b+ABX*I_ERI_G2y2z_S_D2y_S_b;
  Double I_ERI_Gy3z_Px_D2y_S_b = I_ERI_Hxy3z_S_D2y_S_b+ABX*I_ERI_Gy3z_S_D2y_S_b;
  Double I_ERI_G3xy_Py_D2y_S_b = I_ERI_H3x2y_S_D2y_S_b+ABY*I_ERI_G3xy_S_D2y_S_b;
  Double I_ERI_G2x2y_Py_D2y_S_b = I_ERI_H2x3y_S_D2y_S_b+ABY*I_ERI_G2x2y_S_D2y_S_b;
  Double I_ERI_G2xyz_Py_D2y_S_b = I_ERI_H2x2yz_S_D2y_S_b+ABY*I_ERI_G2xyz_S_D2y_S_b;
  Double I_ERI_Gx3y_Py_D2y_S_b = I_ERI_Hx4y_S_D2y_S_b+ABY*I_ERI_Gx3y_S_D2y_S_b;
  Double I_ERI_Gx2yz_Py_D2y_S_b = I_ERI_Hx3yz_S_D2y_S_b+ABY*I_ERI_Gx2yz_S_D2y_S_b;
  Double I_ERI_Gxy2z_Py_D2y_S_b = I_ERI_Hx2y2z_S_D2y_S_b+ABY*I_ERI_Gxy2z_S_D2y_S_b;
  Double I_ERI_G4y_Py_D2y_S_b = I_ERI_H5y_S_D2y_S_b+ABY*I_ERI_G4y_S_D2y_S_b;
  Double I_ERI_G3yz_Py_D2y_S_b = I_ERI_H4yz_S_D2y_S_b+ABY*I_ERI_G3yz_S_D2y_S_b;
  Double I_ERI_G2y2z_Py_D2y_S_b = I_ERI_H3y2z_S_D2y_S_b+ABY*I_ERI_G2y2z_S_D2y_S_b;
  Double I_ERI_Gy3z_Py_D2y_S_b = I_ERI_H2y3z_S_D2y_S_b+ABY*I_ERI_Gy3z_S_D2y_S_b;
  Double I_ERI_G3xz_Pz_D2y_S_b = I_ERI_H3x2z_S_D2y_S_b+ABZ*I_ERI_G3xz_S_D2y_S_b;
  Double I_ERI_G2xyz_Pz_D2y_S_b = I_ERI_H2xy2z_S_D2y_S_b+ABZ*I_ERI_G2xyz_S_D2y_S_b;
  Double I_ERI_G2x2z_Pz_D2y_S_b = I_ERI_H2x3z_S_D2y_S_b+ABZ*I_ERI_G2x2z_S_D2y_S_b;
  Double I_ERI_Gx2yz_Pz_D2y_S_b = I_ERI_Hx2y2z_S_D2y_S_b+ABZ*I_ERI_Gx2yz_S_D2y_S_b;
  Double I_ERI_Gxy2z_Pz_D2y_S_b = I_ERI_Hxy3z_S_D2y_S_b+ABZ*I_ERI_Gxy2z_S_D2y_S_b;
  Double I_ERI_Gx3z_Pz_D2y_S_b = I_ERI_Hx4z_S_D2y_S_b+ABZ*I_ERI_Gx3z_S_D2y_S_b;
  Double I_ERI_G3yz_Pz_D2y_S_b = I_ERI_H3y2z_S_D2y_S_b+ABZ*I_ERI_G3yz_S_D2y_S_b;
  Double I_ERI_G2y2z_Pz_D2y_S_b = I_ERI_H2y3z_S_D2y_S_b+ABZ*I_ERI_G2y2z_S_D2y_S_b;
  Double I_ERI_Gy3z_Pz_D2y_S_b = I_ERI_Hy4z_S_D2y_S_b+ABZ*I_ERI_Gy3z_S_D2y_S_b;
  Double I_ERI_G4z_Pz_D2y_S_b = I_ERI_H5z_S_D2y_S_b+ABZ*I_ERI_G4z_S_D2y_S_b;
  Double I_ERI_G4x_Px_Dyz_S_b = I_ERI_H5x_S_Dyz_S_b+ABX*I_ERI_G4x_S_Dyz_S_b;
  Double I_ERI_G3xy_Px_Dyz_S_b = I_ERI_H4xy_S_Dyz_S_b+ABX*I_ERI_G3xy_S_Dyz_S_b;
  Double I_ERI_G3xz_Px_Dyz_S_b = I_ERI_H4xz_S_Dyz_S_b+ABX*I_ERI_G3xz_S_Dyz_S_b;
  Double I_ERI_G2x2y_Px_Dyz_S_b = I_ERI_H3x2y_S_Dyz_S_b+ABX*I_ERI_G2x2y_S_Dyz_S_b;
  Double I_ERI_G2xyz_Px_Dyz_S_b = I_ERI_H3xyz_S_Dyz_S_b+ABX*I_ERI_G2xyz_S_Dyz_S_b;
  Double I_ERI_G2x2z_Px_Dyz_S_b = I_ERI_H3x2z_S_Dyz_S_b+ABX*I_ERI_G2x2z_S_Dyz_S_b;
  Double I_ERI_Gx3y_Px_Dyz_S_b = I_ERI_H2x3y_S_Dyz_S_b+ABX*I_ERI_Gx3y_S_Dyz_S_b;
  Double I_ERI_Gx2yz_Px_Dyz_S_b = I_ERI_H2x2yz_S_Dyz_S_b+ABX*I_ERI_Gx2yz_S_Dyz_S_b;
  Double I_ERI_Gxy2z_Px_Dyz_S_b = I_ERI_H2xy2z_S_Dyz_S_b+ABX*I_ERI_Gxy2z_S_Dyz_S_b;
  Double I_ERI_Gx3z_Px_Dyz_S_b = I_ERI_H2x3z_S_Dyz_S_b+ABX*I_ERI_Gx3z_S_Dyz_S_b;
  Double I_ERI_G3yz_Px_Dyz_S_b = I_ERI_Hx3yz_S_Dyz_S_b+ABX*I_ERI_G3yz_S_Dyz_S_b;
  Double I_ERI_G2y2z_Px_Dyz_S_b = I_ERI_Hx2y2z_S_Dyz_S_b+ABX*I_ERI_G2y2z_S_Dyz_S_b;
  Double I_ERI_Gy3z_Px_Dyz_S_b = I_ERI_Hxy3z_S_Dyz_S_b+ABX*I_ERI_Gy3z_S_Dyz_S_b;
  Double I_ERI_G3xy_Py_Dyz_S_b = I_ERI_H3x2y_S_Dyz_S_b+ABY*I_ERI_G3xy_S_Dyz_S_b;
  Double I_ERI_G2x2y_Py_Dyz_S_b = I_ERI_H2x3y_S_Dyz_S_b+ABY*I_ERI_G2x2y_S_Dyz_S_b;
  Double I_ERI_G2xyz_Py_Dyz_S_b = I_ERI_H2x2yz_S_Dyz_S_b+ABY*I_ERI_G2xyz_S_Dyz_S_b;
  Double I_ERI_Gx3y_Py_Dyz_S_b = I_ERI_Hx4y_S_Dyz_S_b+ABY*I_ERI_Gx3y_S_Dyz_S_b;
  Double I_ERI_Gx2yz_Py_Dyz_S_b = I_ERI_Hx3yz_S_Dyz_S_b+ABY*I_ERI_Gx2yz_S_Dyz_S_b;
  Double I_ERI_Gxy2z_Py_Dyz_S_b = I_ERI_Hx2y2z_S_Dyz_S_b+ABY*I_ERI_Gxy2z_S_Dyz_S_b;
  Double I_ERI_G4y_Py_Dyz_S_b = I_ERI_H5y_S_Dyz_S_b+ABY*I_ERI_G4y_S_Dyz_S_b;
  Double I_ERI_G3yz_Py_Dyz_S_b = I_ERI_H4yz_S_Dyz_S_b+ABY*I_ERI_G3yz_S_Dyz_S_b;
  Double I_ERI_G2y2z_Py_Dyz_S_b = I_ERI_H3y2z_S_Dyz_S_b+ABY*I_ERI_G2y2z_S_Dyz_S_b;
  Double I_ERI_Gy3z_Py_Dyz_S_b = I_ERI_H2y3z_S_Dyz_S_b+ABY*I_ERI_Gy3z_S_Dyz_S_b;
  Double I_ERI_G3xz_Pz_Dyz_S_b = I_ERI_H3x2z_S_Dyz_S_b+ABZ*I_ERI_G3xz_S_Dyz_S_b;
  Double I_ERI_G2xyz_Pz_Dyz_S_b = I_ERI_H2xy2z_S_Dyz_S_b+ABZ*I_ERI_G2xyz_S_Dyz_S_b;
  Double I_ERI_G2x2z_Pz_Dyz_S_b = I_ERI_H2x3z_S_Dyz_S_b+ABZ*I_ERI_G2x2z_S_Dyz_S_b;
  Double I_ERI_Gx2yz_Pz_Dyz_S_b = I_ERI_Hx2y2z_S_Dyz_S_b+ABZ*I_ERI_Gx2yz_S_Dyz_S_b;
  Double I_ERI_Gxy2z_Pz_Dyz_S_b = I_ERI_Hxy3z_S_Dyz_S_b+ABZ*I_ERI_Gxy2z_S_Dyz_S_b;
  Double I_ERI_Gx3z_Pz_Dyz_S_b = I_ERI_Hx4z_S_Dyz_S_b+ABZ*I_ERI_Gx3z_S_Dyz_S_b;
  Double I_ERI_G3yz_Pz_Dyz_S_b = I_ERI_H3y2z_S_Dyz_S_b+ABZ*I_ERI_G3yz_S_Dyz_S_b;
  Double I_ERI_G2y2z_Pz_Dyz_S_b = I_ERI_H2y3z_S_Dyz_S_b+ABZ*I_ERI_G2y2z_S_Dyz_S_b;
  Double I_ERI_Gy3z_Pz_Dyz_S_b = I_ERI_Hy4z_S_Dyz_S_b+ABZ*I_ERI_Gy3z_S_Dyz_S_b;
  Double I_ERI_G4z_Pz_Dyz_S_b = I_ERI_H5z_S_Dyz_S_b+ABZ*I_ERI_G4z_S_Dyz_S_b;
  Double I_ERI_G4x_Px_D2z_S_b = I_ERI_H5x_S_D2z_S_b+ABX*I_ERI_G4x_S_D2z_S_b;
  Double I_ERI_G3xy_Px_D2z_S_b = I_ERI_H4xy_S_D2z_S_b+ABX*I_ERI_G3xy_S_D2z_S_b;
  Double I_ERI_G3xz_Px_D2z_S_b = I_ERI_H4xz_S_D2z_S_b+ABX*I_ERI_G3xz_S_D2z_S_b;
  Double I_ERI_G2x2y_Px_D2z_S_b = I_ERI_H3x2y_S_D2z_S_b+ABX*I_ERI_G2x2y_S_D2z_S_b;
  Double I_ERI_G2xyz_Px_D2z_S_b = I_ERI_H3xyz_S_D2z_S_b+ABX*I_ERI_G2xyz_S_D2z_S_b;
  Double I_ERI_G2x2z_Px_D2z_S_b = I_ERI_H3x2z_S_D2z_S_b+ABX*I_ERI_G2x2z_S_D2z_S_b;
  Double I_ERI_Gx3y_Px_D2z_S_b = I_ERI_H2x3y_S_D2z_S_b+ABX*I_ERI_Gx3y_S_D2z_S_b;
  Double I_ERI_Gx2yz_Px_D2z_S_b = I_ERI_H2x2yz_S_D2z_S_b+ABX*I_ERI_Gx2yz_S_D2z_S_b;
  Double I_ERI_Gxy2z_Px_D2z_S_b = I_ERI_H2xy2z_S_D2z_S_b+ABX*I_ERI_Gxy2z_S_D2z_S_b;
  Double I_ERI_Gx3z_Px_D2z_S_b = I_ERI_H2x3z_S_D2z_S_b+ABX*I_ERI_Gx3z_S_D2z_S_b;
  Double I_ERI_G3yz_Px_D2z_S_b = I_ERI_Hx3yz_S_D2z_S_b+ABX*I_ERI_G3yz_S_D2z_S_b;
  Double I_ERI_G2y2z_Px_D2z_S_b = I_ERI_Hx2y2z_S_D2z_S_b+ABX*I_ERI_G2y2z_S_D2z_S_b;
  Double I_ERI_Gy3z_Px_D2z_S_b = I_ERI_Hxy3z_S_D2z_S_b+ABX*I_ERI_Gy3z_S_D2z_S_b;
  Double I_ERI_G3xy_Py_D2z_S_b = I_ERI_H3x2y_S_D2z_S_b+ABY*I_ERI_G3xy_S_D2z_S_b;
  Double I_ERI_G2x2y_Py_D2z_S_b = I_ERI_H2x3y_S_D2z_S_b+ABY*I_ERI_G2x2y_S_D2z_S_b;
  Double I_ERI_G2xyz_Py_D2z_S_b = I_ERI_H2x2yz_S_D2z_S_b+ABY*I_ERI_G2xyz_S_D2z_S_b;
  Double I_ERI_Gx3y_Py_D2z_S_b = I_ERI_Hx4y_S_D2z_S_b+ABY*I_ERI_Gx3y_S_D2z_S_b;
  Double I_ERI_Gx2yz_Py_D2z_S_b = I_ERI_Hx3yz_S_D2z_S_b+ABY*I_ERI_Gx2yz_S_D2z_S_b;
  Double I_ERI_Gxy2z_Py_D2z_S_b = I_ERI_Hx2y2z_S_D2z_S_b+ABY*I_ERI_Gxy2z_S_D2z_S_b;
  Double I_ERI_G4y_Py_D2z_S_b = I_ERI_H5y_S_D2z_S_b+ABY*I_ERI_G4y_S_D2z_S_b;
  Double I_ERI_G3yz_Py_D2z_S_b = I_ERI_H4yz_S_D2z_S_b+ABY*I_ERI_G3yz_S_D2z_S_b;
  Double I_ERI_G2y2z_Py_D2z_S_b = I_ERI_H3y2z_S_D2z_S_b+ABY*I_ERI_G2y2z_S_D2z_S_b;
  Double I_ERI_Gy3z_Py_D2z_S_b = I_ERI_H2y3z_S_D2z_S_b+ABY*I_ERI_Gy3z_S_D2z_S_b;
  Double I_ERI_G3xz_Pz_D2z_S_b = I_ERI_H3x2z_S_D2z_S_b+ABZ*I_ERI_G3xz_S_D2z_S_b;
  Double I_ERI_G2xyz_Pz_D2z_S_b = I_ERI_H2xy2z_S_D2z_S_b+ABZ*I_ERI_G2xyz_S_D2z_S_b;
  Double I_ERI_G2x2z_Pz_D2z_S_b = I_ERI_H2x3z_S_D2z_S_b+ABZ*I_ERI_G2x2z_S_D2z_S_b;
  Double I_ERI_Gx2yz_Pz_D2z_S_b = I_ERI_Hx2y2z_S_D2z_S_b+ABZ*I_ERI_Gx2yz_S_D2z_S_b;
  Double I_ERI_Gxy2z_Pz_D2z_S_b = I_ERI_Hxy3z_S_D2z_S_b+ABZ*I_ERI_Gxy2z_S_D2z_S_b;
  Double I_ERI_Gx3z_Pz_D2z_S_b = I_ERI_Hx4z_S_D2z_S_b+ABZ*I_ERI_Gx3z_S_D2z_S_b;
  Double I_ERI_G3yz_Pz_D2z_S_b = I_ERI_H3y2z_S_D2z_S_b+ABZ*I_ERI_G3yz_S_D2z_S_b;
  Double I_ERI_G2y2z_Pz_D2z_S_b = I_ERI_H2y3z_S_D2z_S_b+ABZ*I_ERI_G2y2z_S_D2z_S_b;
  Double I_ERI_Gy3z_Pz_D2z_S_b = I_ERI_Hy4z_S_D2z_S_b+ABZ*I_ERI_Gy3z_S_D2z_S_b;
  Double I_ERI_G4z_Pz_D2z_S_b = I_ERI_H5z_S_D2z_S_b+ABZ*I_ERI_G4z_S_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 144 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_D_S_b
   * RHS shell quartet name: SQ_ERI_F_P_D_S_b
   ************************************************************/
  Double I_ERI_F3x_D2x_D2x_S_b = I_ERI_G4x_Px_D2x_S_b+ABX*I_ERI_F3x_Px_D2x_S_b;
  Double I_ERI_F2xy_D2x_D2x_S_b = I_ERI_G3xy_Px_D2x_S_b+ABX*I_ERI_F2xy_Px_D2x_S_b;
  Double I_ERI_F2xz_D2x_D2x_S_b = I_ERI_G3xz_Px_D2x_S_b+ABX*I_ERI_F2xz_Px_D2x_S_b;
  Double I_ERI_Fx2y_D2x_D2x_S_b = I_ERI_G2x2y_Px_D2x_S_b+ABX*I_ERI_Fx2y_Px_D2x_S_b;
  Double I_ERI_Fxyz_D2x_D2x_S_b = I_ERI_G2xyz_Px_D2x_S_b+ABX*I_ERI_Fxyz_Px_D2x_S_b;
  Double I_ERI_Fx2z_D2x_D2x_S_b = I_ERI_G2x2z_Px_D2x_S_b+ABX*I_ERI_Fx2z_Px_D2x_S_b;
  Double I_ERI_F3y_D2x_D2x_S_b = I_ERI_Gx3y_Px_D2x_S_b+ABX*I_ERI_F3y_Px_D2x_S_b;
  Double I_ERI_F2yz_D2x_D2x_S_b = I_ERI_Gx2yz_Px_D2x_S_b+ABX*I_ERI_F2yz_Px_D2x_S_b;
  Double I_ERI_Fy2z_D2x_D2x_S_b = I_ERI_Gxy2z_Px_D2x_S_b+ABX*I_ERI_Fy2z_Px_D2x_S_b;
  Double I_ERI_F3z_D2x_D2x_S_b = I_ERI_Gx3z_Px_D2x_S_b+ABX*I_ERI_F3z_Px_D2x_S_b;
  Double I_ERI_F2xz_Dxy_D2x_S_b = I_ERI_G2xyz_Px_D2x_S_b+ABY*I_ERI_F2xz_Px_D2x_S_b;
  Double I_ERI_Fxyz_Dxy_D2x_S_b = I_ERI_Gx2yz_Px_D2x_S_b+ABY*I_ERI_Fxyz_Px_D2x_S_b;
  Double I_ERI_Fx2z_Dxy_D2x_S_b = I_ERI_Gxy2z_Px_D2x_S_b+ABY*I_ERI_Fx2z_Px_D2x_S_b;
  Double I_ERI_F2yz_Dxy_D2x_S_b = I_ERI_G3yz_Px_D2x_S_b+ABY*I_ERI_F2yz_Px_D2x_S_b;
  Double I_ERI_Fy2z_Dxy_D2x_S_b = I_ERI_G2y2z_Px_D2x_S_b+ABY*I_ERI_Fy2z_Px_D2x_S_b;
  Double I_ERI_F3z_Dxy_D2x_S_b = I_ERI_Gy3z_Px_D2x_S_b+ABY*I_ERI_F3z_Px_D2x_S_b;
  Double I_ERI_F3x_D2y_D2x_S_b = I_ERI_G3xy_Py_D2x_S_b+ABY*I_ERI_F3x_Py_D2x_S_b;
  Double I_ERI_F2xy_D2y_D2x_S_b = I_ERI_G2x2y_Py_D2x_S_b+ABY*I_ERI_F2xy_Py_D2x_S_b;
  Double I_ERI_F2xz_D2y_D2x_S_b = I_ERI_G2xyz_Py_D2x_S_b+ABY*I_ERI_F2xz_Py_D2x_S_b;
  Double I_ERI_Fx2y_D2y_D2x_S_b = I_ERI_Gx3y_Py_D2x_S_b+ABY*I_ERI_Fx2y_Py_D2x_S_b;
  Double I_ERI_Fxyz_D2y_D2x_S_b = I_ERI_Gx2yz_Py_D2x_S_b+ABY*I_ERI_Fxyz_Py_D2x_S_b;
  Double I_ERI_Fx2z_D2y_D2x_S_b = I_ERI_Gxy2z_Py_D2x_S_b+ABY*I_ERI_Fx2z_Py_D2x_S_b;
  Double I_ERI_F3y_D2y_D2x_S_b = I_ERI_G4y_Py_D2x_S_b+ABY*I_ERI_F3y_Py_D2x_S_b;
  Double I_ERI_F2yz_D2y_D2x_S_b = I_ERI_G3yz_Py_D2x_S_b+ABY*I_ERI_F2yz_Py_D2x_S_b;
  Double I_ERI_Fy2z_D2y_D2x_S_b = I_ERI_G2y2z_Py_D2x_S_b+ABY*I_ERI_Fy2z_Py_D2x_S_b;
  Double I_ERI_F3z_D2y_D2x_S_b = I_ERI_Gy3z_Py_D2x_S_b+ABY*I_ERI_F3z_Py_D2x_S_b;
  Double I_ERI_F3x_D2z_D2x_S_b = I_ERI_G3xz_Pz_D2x_S_b+ABZ*I_ERI_F3x_Pz_D2x_S_b;
  Double I_ERI_F2xy_D2z_D2x_S_b = I_ERI_G2xyz_Pz_D2x_S_b+ABZ*I_ERI_F2xy_Pz_D2x_S_b;
  Double I_ERI_F2xz_D2z_D2x_S_b = I_ERI_G2x2z_Pz_D2x_S_b+ABZ*I_ERI_F2xz_Pz_D2x_S_b;
  Double I_ERI_Fx2y_D2z_D2x_S_b = I_ERI_Gx2yz_Pz_D2x_S_b+ABZ*I_ERI_Fx2y_Pz_D2x_S_b;
  Double I_ERI_Fxyz_D2z_D2x_S_b = I_ERI_Gxy2z_Pz_D2x_S_b+ABZ*I_ERI_Fxyz_Pz_D2x_S_b;
  Double I_ERI_Fx2z_D2z_D2x_S_b = I_ERI_Gx3z_Pz_D2x_S_b+ABZ*I_ERI_Fx2z_Pz_D2x_S_b;
  Double I_ERI_F3y_D2z_D2x_S_b = I_ERI_G3yz_Pz_D2x_S_b+ABZ*I_ERI_F3y_Pz_D2x_S_b;
  Double I_ERI_F2yz_D2z_D2x_S_b = I_ERI_G2y2z_Pz_D2x_S_b+ABZ*I_ERI_F2yz_Pz_D2x_S_b;
  Double I_ERI_Fy2z_D2z_D2x_S_b = I_ERI_Gy3z_Pz_D2x_S_b+ABZ*I_ERI_Fy2z_Pz_D2x_S_b;
  Double I_ERI_F3z_D2z_D2x_S_b = I_ERI_G4z_Pz_D2x_S_b+ABZ*I_ERI_F3z_Pz_D2x_S_b;
  Double I_ERI_F3x_D2x_Dxy_S_b = I_ERI_G4x_Px_Dxy_S_b+ABX*I_ERI_F3x_Px_Dxy_S_b;
  Double I_ERI_F2xy_D2x_Dxy_S_b = I_ERI_G3xy_Px_Dxy_S_b+ABX*I_ERI_F2xy_Px_Dxy_S_b;
  Double I_ERI_F2xz_D2x_Dxy_S_b = I_ERI_G3xz_Px_Dxy_S_b+ABX*I_ERI_F2xz_Px_Dxy_S_b;
  Double I_ERI_Fx2y_D2x_Dxy_S_b = I_ERI_G2x2y_Px_Dxy_S_b+ABX*I_ERI_Fx2y_Px_Dxy_S_b;
  Double I_ERI_Fxyz_D2x_Dxy_S_b = I_ERI_G2xyz_Px_Dxy_S_b+ABX*I_ERI_Fxyz_Px_Dxy_S_b;
  Double I_ERI_Fx2z_D2x_Dxy_S_b = I_ERI_G2x2z_Px_Dxy_S_b+ABX*I_ERI_Fx2z_Px_Dxy_S_b;
  Double I_ERI_F3y_D2x_Dxy_S_b = I_ERI_Gx3y_Px_Dxy_S_b+ABX*I_ERI_F3y_Px_Dxy_S_b;
  Double I_ERI_F2yz_D2x_Dxy_S_b = I_ERI_Gx2yz_Px_Dxy_S_b+ABX*I_ERI_F2yz_Px_Dxy_S_b;
  Double I_ERI_Fy2z_D2x_Dxy_S_b = I_ERI_Gxy2z_Px_Dxy_S_b+ABX*I_ERI_Fy2z_Px_Dxy_S_b;
  Double I_ERI_F3z_D2x_Dxy_S_b = I_ERI_Gx3z_Px_Dxy_S_b+ABX*I_ERI_F3z_Px_Dxy_S_b;
  Double I_ERI_F2xz_Dxy_Dxy_S_b = I_ERI_G2xyz_Px_Dxy_S_b+ABY*I_ERI_F2xz_Px_Dxy_S_b;
  Double I_ERI_Fxyz_Dxy_Dxy_S_b = I_ERI_Gx2yz_Px_Dxy_S_b+ABY*I_ERI_Fxyz_Px_Dxy_S_b;
  Double I_ERI_Fx2z_Dxy_Dxy_S_b = I_ERI_Gxy2z_Px_Dxy_S_b+ABY*I_ERI_Fx2z_Px_Dxy_S_b;
  Double I_ERI_F2yz_Dxy_Dxy_S_b = I_ERI_G3yz_Px_Dxy_S_b+ABY*I_ERI_F2yz_Px_Dxy_S_b;
  Double I_ERI_Fy2z_Dxy_Dxy_S_b = I_ERI_G2y2z_Px_Dxy_S_b+ABY*I_ERI_Fy2z_Px_Dxy_S_b;
  Double I_ERI_F3z_Dxy_Dxy_S_b = I_ERI_Gy3z_Px_Dxy_S_b+ABY*I_ERI_F3z_Px_Dxy_S_b;
  Double I_ERI_F3x_D2y_Dxy_S_b = I_ERI_G3xy_Py_Dxy_S_b+ABY*I_ERI_F3x_Py_Dxy_S_b;
  Double I_ERI_F2xy_D2y_Dxy_S_b = I_ERI_G2x2y_Py_Dxy_S_b+ABY*I_ERI_F2xy_Py_Dxy_S_b;
  Double I_ERI_F2xz_D2y_Dxy_S_b = I_ERI_G2xyz_Py_Dxy_S_b+ABY*I_ERI_F2xz_Py_Dxy_S_b;
  Double I_ERI_Fx2y_D2y_Dxy_S_b = I_ERI_Gx3y_Py_Dxy_S_b+ABY*I_ERI_Fx2y_Py_Dxy_S_b;
  Double I_ERI_Fxyz_D2y_Dxy_S_b = I_ERI_Gx2yz_Py_Dxy_S_b+ABY*I_ERI_Fxyz_Py_Dxy_S_b;
  Double I_ERI_Fx2z_D2y_Dxy_S_b = I_ERI_Gxy2z_Py_Dxy_S_b+ABY*I_ERI_Fx2z_Py_Dxy_S_b;
  Double I_ERI_F3y_D2y_Dxy_S_b = I_ERI_G4y_Py_Dxy_S_b+ABY*I_ERI_F3y_Py_Dxy_S_b;
  Double I_ERI_F2yz_D2y_Dxy_S_b = I_ERI_G3yz_Py_Dxy_S_b+ABY*I_ERI_F2yz_Py_Dxy_S_b;
  Double I_ERI_Fy2z_D2y_Dxy_S_b = I_ERI_G2y2z_Py_Dxy_S_b+ABY*I_ERI_Fy2z_Py_Dxy_S_b;
  Double I_ERI_F3z_D2y_Dxy_S_b = I_ERI_Gy3z_Py_Dxy_S_b+ABY*I_ERI_F3z_Py_Dxy_S_b;
  Double I_ERI_F3x_D2z_Dxy_S_b = I_ERI_G3xz_Pz_Dxy_S_b+ABZ*I_ERI_F3x_Pz_Dxy_S_b;
  Double I_ERI_F2xy_D2z_Dxy_S_b = I_ERI_G2xyz_Pz_Dxy_S_b+ABZ*I_ERI_F2xy_Pz_Dxy_S_b;
  Double I_ERI_F2xz_D2z_Dxy_S_b = I_ERI_G2x2z_Pz_Dxy_S_b+ABZ*I_ERI_F2xz_Pz_Dxy_S_b;
  Double I_ERI_Fx2y_D2z_Dxy_S_b = I_ERI_Gx2yz_Pz_Dxy_S_b+ABZ*I_ERI_Fx2y_Pz_Dxy_S_b;
  Double I_ERI_Fxyz_D2z_Dxy_S_b = I_ERI_Gxy2z_Pz_Dxy_S_b+ABZ*I_ERI_Fxyz_Pz_Dxy_S_b;
  Double I_ERI_Fx2z_D2z_Dxy_S_b = I_ERI_Gx3z_Pz_Dxy_S_b+ABZ*I_ERI_Fx2z_Pz_Dxy_S_b;
  Double I_ERI_F3y_D2z_Dxy_S_b = I_ERI_G3yz_Pz_Dxy_S_b+ABZ*I_ERI_F3y_Pz_Dxy_S_b;
  Double I_ERI_F2yz_D2z_Dxy_S_b = I_ERI_G2y2z_Pz_Dxy_S_b+ABZ*I_ERI_F2yz_Pz_Dxy_S_b;
  Double I_ERI_Fy2z_D2z_Dxy_S_b = I_ERI_Gy3z_Pz_Dxy_S_b+ABZ*I_ERI_Fy2z_Pz_Dxy_S_b;
  Double I_ERI_F3z_D2z_Dxy_S_b = I_ERI_G4z_Pz_Dxy_S_b+ABZ*I_ERI_F3z_Pz_Dxy_S_b;
  Double I_ERI_F3x_D2x_Dxz_S_b = I_ERI_G4x_Px_Dxz_S_b+ABX*I_ERI_F3x_Px_Dxz_S_b;
  Double I_ERI_F2xy_D2x_Dxz_S_b = I_ERI_G3xy_Px_Dxz_S_b+ABX*I_ERI_F2xy_Px_Dxz_S_b;
  Double I_ERI_F2xz_D2x_Dxz_S_b = I_ERI_G3xz_Px_Dxz_S_b+ABX*I_ERI_F2xz_Px_Dxz_S_b;
  Double I_ERI_Fx2y_D2x_Dxz_S_b = I_ERI_G2x2y_Px_Dxz_S_b+ABX*I_ERI_Fx2y_Px_Dxz_S_b;
  Double I_ERI_Fxyz_D2x_Dxz_S_b = I_ERI_G2xyz_Px_Dxz_S_b+ABX*I_ERI_Fxyz_Px_Dxz_S_b;
  Double I_ERI_Fx2z_D2x_Dxz_S_b = I_ERI_G2x2z_Px_Dxz_S_b+ABX*I_ERI_Fx2z_Px_Dxz_S_b;
  Double I_ERI_F3y_D2x_Dxz_S_b = I_ERI_Gx3y_Px_Dxz_S_b+ABX*I_ERI_F3y_Px_Dxz_S_b;
  Double I_ERI_F2yz_D2x_Dxz_S_b = I_ERI_Gx2yz_Px_Dxz_S_b+ABX*I_ERI_F2yz_Px_Dxz_S_b;
  Double I_ERI_Fy2z_D2x_Dxz_S_b = I_ERI_Gxy2z_Px_Dxz_S_b+ABX*I_ERI_Fy2z_Px_Dxz_S_b;
  Double I_ERI_F3z_D2x_Dxz_S_b = I_ERI_Gx3z_Px_Dxz_S_b+ABX*I_ERI_F3z_Px_Dxz_S_b;
  Double I_ERI_F2xz_Dxy_Dxz_S_b = I_ERI_G2xyz_Px_Dxz_S_b+ABY*I_ERI_F2xz_Px_Dxz_S_b;
  Double I_ERI_Fxyz_Dxy_Dxz_S_b = I_ERI_Gx2yz_Px_Dxz_S_b+ABY*I_ERI_Fxyz_Px_Dxz_S_b;
  Double I_ERI_Fx2z_Dxy_Dxz_S_b = I_ERI_Gxy2z_Px_Dxz_S_b+ABY*I_ERI_Fx2z_Px_Dxz_S_b;
  Double I_ERI_F2yz_Dxy_Dxz_S_b = I_ERI_G3yz_Px_Dxz_S_b+ABY*I_ERI_F2yz_Px_Dxz_S_b;
  Double I_ERI_Fy2z_Dxy_Dxz_S_b = I_ERI_G2y2z_Px_Dxz_S_b+ABY*I_ERI_Fy2z_Px_Dxz_S_b;
  Double I_ERI_F3z_Dxy_Dxz_S_b = I_ERI_Gy3z_Px_Dxz_S_b+ABY*I_ERI_F3z_Px_Dxz_S_b;
  Double I_ERI_F3x_D2y_Dxz_S_b = I_ERI_G3xy_Py_Dxz_S_b+ABY*I_ERI_F3x_Py_Dxz_S_b;
  Double I_ERI_F2xy_D2y_Dxz_S_b = I_ERI_G2x2y_Py_Dxz_S_b+ABY*I_ERI_F2xy_Py_Dxz_S_b;
  Double I_ERI_F2xz_D2y_Dxz_S_b = I_ERI_G2xyz_Py_Dxz_S_b+ABY*I_ERI_F2xz_Py_Dxz_S_b;
  Double I_ERI_Fx2y_D2y_Dxz_S_b = I_ERI_Gx3y_Py_Dxz_S_b+ABY*I_ERI_Fx2y_Py_Dxz_S_b;
  Double I_ERI_Fxyz_D2y_Dxz_S_b = I_ERI_Gx2yz_Py_Dxz_S_b+ABY*I_ERI_Fxyz_Py_Dxz_S_b;
  Double I_ERI_Fx2z_D2y_Dxz_S_b = I_ERI_Gxy2z_Py_Dxz_S_b+ABY*I_ERI_Fx2z_Py_Dxz_S_b;
  Double I_ERI_F3y_D2y_Dxz_S_b = I_ERI_G4y_Py_Dxz_S_b+ABY*I_ERI_F3y_Py_Dxz_S_b;
  Double I_ERI_F2yz_D2y_Dxz_S_b = I_ERI_G3yz_Py_Dxz_S_b+ABY*I_ERI_F2yz_Py_Dxz_S_b;
  Double I_ERI_Fy2z_D2y_Dxz_S_b = I_ERI_G2y2z_Py_Dxz_S_b+ABY*I_ERI_Fy2z_Py_Dxz_S_b;
  Double I_ERI_F3z_D2y_Dxz_S_b = I_ERI_Gy3z_Py_Dxz_S_b+ABY*I_ERI_F3z_Py_Dxz_S_b;
  Double I_ERI_F3x_D2z_Dxz_S_b = I_ERI_G3xz_Pz_Dxz_S_b+ABZ*I_ERI_F3x_Pz_Dxz_S_b;
  Double I_ERI_F2xy_D2z_Dxz_S_b = I_ERI_G2xyz_Pz_Dxz_S_b+ABZ*I_ERI_F2xy_Pz_Dxz_S_b;
  Double I_ERI_F2xz_D2z_Dxz_S_b = I_ERI_G2x2z_Pz_Dxz_S_b+ABZ*I_ERI_F2xz_Pz_Dxz_S_b;
  Double I_ERI_Fx2y_D2z_Dxz_S_b = I_ERI_Gx2yz_Pz_Dxz_S_b+ABZ*I_ERI_Fx2y_Pz_Dxz_S_b;
  Double I_ERI_Fxyz_D2z_Dxz_S_b = I_ERI_Gxy2z_Pz_Dxz_S_b+ABZ*I_ERI_Fxyz_Pz_Dxz_S_b;
  Double I_ERI_Fx2z_D2z_Dxz_S_b = I_ERI_Gx3z_Pz_Dxz_S_b+ABZ*I_ERI_Fx2z_Pz_Dxz_S_b;
  Double I_ERI_F3y_D2z_Dxz_S_b = I_ERI_G3yz_Pz_Dxz_S_b+ABZ*I_ERI_F3y_Pz_Dxz_S_b;
  Double I_ERI_F2yz_D2z_Dxz_S_b = I_ERI_G2y2z_Pz_Dxz_S_b+ABZ*I_ERI_F2yz_Pz_Dxz_S_b;
  Double I_ERI_Fy2z_D2z_Dxz_S_b = I_ERI_Gy3z_Pz_Dxz_S_b+ABZ*I_ERI_Fy2z_Pz_Dxz_S_b;
  Double I_ERI_F3z_D2z_Dxz_S_b = I_ERI_G4z_Pz_Dxz_S_b+ABZ*I_ERI_F3z_Pz_Dxz_S_b;
  Double I_ERI_F3x_D2x_D2y_S_b = I_ERI_G4x_Px_D2y_S_b+ABX*I_ERI_F3x_Px_D2y_S_b;
  Double I_ERI_F2xy_D2x_D2y_S_b = I_ERI_G3xy_Px_D2y_S_b+ABX*I_ERI_F2xy_Px_D2y_S_b;
  Double I_ERI_F2xz_D2x_D2y_S_b = I_ERI_G3xz_Px_D2y_S_b+ABX*I_ERI_F2xz_Px_D2y_S_b;
  Double I_ERI_Fx2y_D2x_D2y_S_b = I_ERI_G2x2y_Px_D2y_S_b+ABX*I_ERI_Fx2y_Px_D2y_S_b;
  Double I_ERI_Fxyz_D2x_D2y_S_b = I_ERI_G2xyz_Px_D2y_S_b+ABX*I_ERI_Fxyz_Px_D2y_S_b;
  Double I_ERI_Fx2z_D2x_D2y_S_b = I_ERI_G2x2z_Px_D2y_S_b+ABX*I_ERI_Fx2z_Px_D2y_S_b;
  Double I_ERI_F3y_D2x_D2y_S_b = I_ERI_Gx3y_Px_D2y_S_b+ABX*I_ERI_F3y_Px_D2y_S_b;
  Double I_ERI_F2yz_D2x_D2y_S_b = I_ERI_Gx2yz_Px_D2y_S_b+ABX*I_ERI_F2yz_Px_D2y_S_b;
  Double I_ERI_Fy2z_D2x_D2y_S_b = I_ERI_Gxy2z_Px_D2y_S_b+ABX*I_ERI_Fy2z_Px_D2y_S_b;
  Double I_ERI_F3z_D2x_D2y_S_b = I_ERI_Gx3z_Px_D2y_S_b+ABX*I_ERI_F3z_Px_D2y_S_b;
  Double I_ERI_F2xz_Dxy_D2y_S_b = I_ERI_G2xyz_Px_D2y_S_b+ABY*I_ERI_F2xz_Px_D2y_S_b;
  Double I_ERI_Fxyz_Dxy_D2y_S_b = I_ERI_Gx2yz_Px_D2y_S_b+ABY*I_ERI_Fxyz_Px_D2y_S_b;
  Double I_ERI_Fx2z_Dxy_D2y_S_b = I_ERI_Gxy2z_Px_D2y_S_b+ABY*I_ERI_Fx2z_Px_D2y_S_b;
  Double I_ERI_F2yz_Dxy_D2y_S_b = I_ERI_G3yz_Px_D2y_S_b+ABY*I_ERI_F2yz_Px_D2y_S_b;
  Double I_ERI_Fy2z_Dxy_D2y_S_b = I_ERI_G2y2z_Px_D2y_S_b+ABY*I_ERI_Fy2z_Px_D2y_S_b;
  Double I_ERI_F3z_Dxy_D2y_S_b = I_ERI_Gy3z_Px_D2y_S_b+ABY*I_ERI_F3z_Px_D2y_S_b;
  Double I_ERI_F3x_D2y_D2y_S_b = I_ERI_G3xy_Py_D2y_S_b+ABY*I_ERI_F3x_Py_D2y_S_b;
  Double I_ERI_F2xy_D2y_D2y_S_b = I_ERI_G2x2y_Py_D2y_S_b+ABY*I_ERI_F2xy_Py_D2y_S_b;
  Double I_ERI_F2xz_D2y_D2y_S_b = I_ERI_G2xyz_Py_D2y_S_b+ABY*I_ERI_F2xz_Py_D2y_S_b;
  Double I_ERI_Fx2y_D2y_D2y_S_b = I_ERI_Gx3y_Py_D2y_S_b+ABY*I_ERI_Fx2y_Py_D2y_S_b;
  Double I_ERI_Fxyz_D2y_D2y_S_b = I_ERI_Gx2yz_Py_D2y_S_b+ABY*I_ERI_Fxyz_Py_D2y_S_b;
  Double I_ERI_Fx2z_D2y_D2y_S_b = I_ERI_Gxy2z_Py_D2y_S_b+ABY*I_ERI_Fx2z_Py_D2y_S_b;
  Double I_ERI_F3y_D2y_D2y_S_b = I_ERI_G4y_Py_D2y_S_b+ABY*I_ERI_F3y_Py_D2y_S_b;
  Double I_ERI_F2yz_D2y_D2y_S_b = I_ERI_G3yz_Py_D2y_S_b+ABY*I_ERI_F2yz_Py_D2y_S_b;
  Double I_ERI_Fy2z_D2y_D2y_S_b = I_ERI_G2y2z_Py_D2y_S_b+ABY*I_ERI_Fy2z_Py_D2y_S_b;
  Double I_ERI_F3z_D2y_D2y_S_b = I_ERI_Gy3z_Py_D2y_S_b+ABY*I_ERI_F3z_Py_D2y_S_b;
  Double I_ERI_F3x_D2z_D2y_S_b = I_ERI_G3xz_Pz_D2y_S_b+ABZ*I_ERI_F3x_Pz_D2y_S_b;
  Double I_ERI_F2xy_D2z_D2y_S_b = I_ERI_G2xyz_Pz_D2y_S_b+ABZ*I_ERI_F2xy_Pz_D2y_S_b;
  Double I_ERI_F2xz_D2z_D2y_S_b = I_ERI_G2x2z_Pz_D2y_S_b+ABZ*I_ERI_F2xz_Pz_D2y_S_b;
  Double I_ERI_Fx2y_D2z_D2y_S_b = I_ERI_Gx2yz_Pz_D2y_S_b+ABZ*I_ERI_Fx2y_Pz_D2y_S_b;
  Double I_ERI_Fxyz_D2z_D2y_S_b = I_ERI_Gxy2z_Pz_D2y_S_b+ABZ*I_ERI_Fxyz_Pz_D2y_S_b;
  Double I_ERI_Fx2z_D2z_D2y_S_b = I_ERI_Gx3z_Pz_D2y_S_b+ABZ*I_ERI_Fx2z_Pz_D2y_S_b;
  Double I_ERI_F3y_D2z_D2y_S_b = I_ERI_G3yz_Pz_D2y_S_b+ABZ*I_ERI_F3y_Pz_D2y_S_b;
  Double I_ERI_F2yz_D2z_D2y_S_b = I_ERI_G2y2z_Pz_D2y_S_b+ABZ*I_ERI_F2yz_Pz_D2y_S_b;
  Double I_ERI_Fy2z_D2z_D2y_S_b = I_ERI_Gy3z_Pz_D2y_S_b+ABZ*I_ERI_Fy2z_Pz_D2y_S_b;
  Double I_ERI_F3z_D2z_D2y_S_b = I_ERI_G4z_Pz_D2y_S_b+ABZ*I_ERI_F3z_Pz_D2y_S_b;
  Double I_ERI_F3x_D2x_Dyz_S_b = I_ERI_G4x_Px_Dyz_S_b+ABX*I_ERI_F3x_Px_Dyz_S_b;
  Double I_ERI_F2xy_D2x_Dyz_S_b = I_ERI_G3xy_Px_Dyz_S_b+ABX*I_ERI_F2xy_Px_Dyz_S_b;
  Double I_ERI_F2xz_D2x_Dyz_S_b = I_ERI_G3xz_Px_Dyz_S_b+ABX*I_ERI_F2xz_Px_Dyz_S_b;
  Double I_ERI_Fx2y_D2x_Dyz_S_b = I_ERI_G2x2y_Px_Dyz_S_b+ABX*I_ERI_Fx2y_Px_Dyz_S_b;
  Double I_ERI_Fxyz_D2x_Dyz_S_b = I_ERI_G2xyz_Px_Dyz_S_b+ABX*I_ERI_Fxyz_Px_Dyz_S_b;
  Double I_ERI_Fx2z_D2x_Dyz_S_b = I_ERI_G2x2z_Px_Dyz_S_b+ABX*I_ERI_Fx2z_Px_Dyz_S_b;
  Double I_ERI_F3y_D2x_Dyz_S_b = I_ERI_Gx3y_Px_Dyz_S_b+ABX*I_ERI_F3y_Px_Dyz_S_b;
  Double I_ERI_F2yz_D2x_Dyz_S_b = I_ERI_Gx2yz_Px_Dyz_S_b+ABX*I_ERI_F2yz_Px_Dyz_S_b;
  Double I_ERI_Fy2z_D2x_Dyz_S_b = I_ERI_Gxy2z_Px_Dyz_S_b+ABX*I_ERI_Fy2z_Px_Dyz_S_b;
  Double I_ERI_F3z_D2x_Dyz_S_b = I_ERI_Gx3z_Px_Dyz_S_b+ABX*I_ERI_F3z_Px_Dyz_S_b;
  Double I_ERI_F2xz_Dxy_Dyz_S_b = I_ERI_G2xyz_Px_Dyz_S_b+ABY*I_ERI_F2xz_Px_Dyz_S_b;
  Double I_ERI_Fxyz_Dxy_Dyz_S_b = I_ERI_Gx2yz_Px_Dyz_S_b+ABY*I_ERI_Fxyz_Px_Dyz_S_b;
  Double I_ERI_Fx2z_Dxy_Dyz_S_b = I_ERI_Gxy2z_Px_Dyz_S_b+ABY*I_ERI_Fx2z_Px_Dyz_S_b;
  Double I_ERI_F2yz_Dxy_Dyz_S_b = I_ERI_G3yz_Px_Dyz_S_b+ABY*I_ERI_F2yz_Px_Dyz_S_b;
  Double I_ERI_Fy2z_Dxy_Dyz_S_b = I_ERI_G2y2z_Px_Dyz_S_b+ABY*I_ERI_Fy2z_Px_Dyz_S_b;
  Double I_ERI_F3z_Dxy_Dyz_S_b = I_ERI_Gy3z_Px_Dyz_S_b+ABY*I_ERI_F3z_Px_Dyz_S_b;
  Double I_ERI_F3x_D2y_Dyz_S_b = I_ERI_G3xy_Py_Dyz_S_b+ABY*I_ERI_F3x_Py_Dyz_S_b;
  Double I_ERI_F2xy_D2y_Dyz_S_b = I_ERI_G2x2y_Py_Dyz_S_b+ABY*I_ERI_F2xy_Py_Dyz_S_b;
  Double I_ERI_F2xz_D2y_Dyz_S_b = I_ERI_G2xyz_Py_Dyz_S_b+ABY*I_ERI_F2xz_Py_Dyz_S_b;
  Double I_ERI_Fx2y_D2y_Dyz_S_b = I_ERI_Gx3y_Py_Dyz_S_b+ABY*I_ERI_Fx2y_Py_Dyz_S_b;
  Double I_ERI_Fxyz_D2y_Dyz_S_b = I_ERI_Gx2yz_Py_Dyz_S_b+ABY*I_ERI_Fxyz_Py_Dyz_S_b;
  Double I_ERI_Fx2z_D2y_Dyz_S_b = I_ERI_Gxy2z_Py_Dyz_S_b+ABY*I_ERI_Fx2z_Py_Dyz_S_b;
  Double I_ERI_F3y_D2y_Dyz_S_b = I_ERI_G4y_Py_Dyz_S_b+ABY*I_ERI_F3y_Py_Dyz_S_b;
  Double I_ERI_F2yz_D2y_Dyz_S_b = I_ERI_G3yz_Py_Dyz_S_b+ABY*I_ERI_F2yz_Py_Dyz_S_b;
  Double I_ERI_Fy2z_D2y_Dyz_S_b = I_ERI_G2y2z_Py_Dyz_S_b+ABY*I_ERI_Fy2z_Py_Dyz_S_b;
  Double I_ERI_F3z_D2y_Dyz_S_b = I_ERI_Gy3z_Py_Dyz_S_b+ABY*I_ERI_F3z_Py_Dyz_S_b;
  Double I_ERI_F3x_D2z_Dyz_S_b = I_ERI_G3xz_Pz_Dyz_S_b+ABZ*I_ERI_F3x_Pz_Dyz_S_b;
  Double I_ERI_F2xy_D2z_Dyz_S_b = I_ERI_G2xyz_Pz_Dyz_S_b+ABZ*I_ERI_F2xy_Pz_Dyz_S_b;
  Double I_ERI_F2xz_D2z_Dyz_S_b = I_ERI_G2x2z_Pz_Dyz_S_b+ABZ*I_ERI_F2xz_Pz_Dyz_S_b;
  Double I_ERI_Fx2y_D2z_Dyz_S_b = I_ERI_Gx2yz_Pz_Dyz_S_b+ABZ*I_ERI_Fx2y_Pz_Dyz_S_b;
  Double I_ERI_Fxyz_D2z_Dyz_S_b = I_ERI_Gxy2z_Pz_Dyz_S_b+ABZ*I_ERI_Fxyz_Pz_Dyz_S_b;
  Double I_ERI_Fx2z_D2z_Dyz_S_b = I_ERI_Gx3z_Pz_Dyz_S_b+ABZ*I_ERI_Fx2z_Pz_Dyz_S_b;
  Double I_ERI_F3y_D2z_Dyz_S_b = I_ERI_G3yz_Pz_Dyz_S_b+ABZ*I_ERI_F3y_Pz_Dyz_S_b;
  Double I_ERI_F2yz_D2z_Dyz_S_b = I_ERI_G2y2z_Pz_Dyz_S_b+ABZ*I_ERI_F2yz_Pz_Dyz_S_b;
  Double I_ERI_Fy2z_D2z_Dyz_S_b = I_ERI_Gy3z_Pz_Dyz_S_b+ABZ*I_ERI_Fy2z_Pz_Dyz_S_b;
  Double I_ERI_F3z_D2z_Dyz_S_b = I_ERI_G4z_Pz_Dyz_S_b+ABZ*I_ERI_F3z_Pz_Dyz_S_b;
  Double I_ERI_F3x_D2x_D2z_S_b = I_ERI_G4x_Px_D2z_S_b+ABX*I_ERI_F3x_Px_D2z_S_b;
  Double I_ERI_F2xy_D2x_D2z_S_b = I_ERI_G3xy_Px_D2z_S_b+ABX*I_ERI_F2xy_Px_D2z_S_b;
  Double I_ERI_F2xz_D2x_D2z_S_b = I_ERI_G3xz_Px_D2z_S_b+ABX*I_ERI_F2xz_Px_D2z_S_b;
  Double I_ERI_Fx2y_D2x_D2z_S_b = I_ERI_G2x2y_Px_D2z_S_b+ABX*I_ERI_Fx2y_Px_D2z_S_b;
  Double I_ERI_Fxyz_D2x_D2z_S_b = I_ERI_G2xyz_Px_D2z_S_b+ABX*I_ERI_Fxyz_Px_D2z_S_b;
  Double I_ERI_Fx2z_D2x_D2z_S_b = I_ERI_G2x2z_Px_D2z_S_b+ABX*I_ERI_Fx2z_Px_D2z_S_b;
  Double I_ERI_F3y_D2x_D2z_S_b = I_ERI_Gx3y_Px_D2z_S_b+ABX*I_ERI_F3y_Px_D2z_S_b;
  Double I_ERI_F2yz_D2x_D2z_S_b = I_ERI_Gx2yz_Px_D2z_S_b+ABX*I_ERI_F2yz_Px_D2z_S_b;
  Double I_ERI_Fy2z_D2x_D2z_S_b = I_ERI_Gxy2z_Px_D2z_S_b+ABX*I_ERI_Fy2z_Px_D2z_S_b;
  Double I_ERI_F3z_D2x_D2z_S_b = I_ERI_Gx3z_Px_D2z_S_b+ABX*I_ERI_F3z_Px_D2z_S_b;
  Double I_ERI_F2xz_Dxy_D2z_S_b = I_ERI_G2xyz_Px_D2z_S_b+ABY*I_ERI_F2xz_Px_D2z_S_b;
  Double I_ERI_Fxyz_Dxy_D2z_S_b = I_ERI_Gx2yz_Px_D2z_S_b+ABY*I_ERI_Fxyz_Px_D2z_S_b;
  Double I_ERI_Fx2z_Dxy_D2z_S_b = I_ERI_Gxy2z_Px_D2z_S_b+ABY*I_ERI_Fx2z_Px_D2z_S_b;
  Double I_ERI_F2yz_Dxy_D2z_S_b = I_ERI_G3yz_Px_D2z_S_b+ABY*I_ERI_F2yz_Px_D2z_S_b;
  Double I_ERI_Fy2z_Dxy_D2z_S_b = I_ERI_G2y2z_Px_D2z_S_b+ABY*I_ERI_Fy2z_Px_D2z_S_b;
  Double I_ERI_F3z_Dxy_D2z_S_b = I_ERI_Gy3z_Px_D2z_S_b+ABY*I_ERI_F3z_Px_D2z_S_b;
  Double I_ERI_F3x_D2y_D2z_S_b = I_ERI_G3xy_Py_D2z_S_b+ABY*I_ERI_F3x_Py_D2z_S_b;
  Double I_ERI_F2xy_D2y_D2z_S_b = I_ERI_G2x2y_Py_D2z_S_b+ABY*I_ERI_F2xy_Py_D2z_S_b;
  Double I_ERI_F2xz_D2y_D2z_S_b = I_ERI_G2xyz_Py_D2z_S_b+ABY*I_ERI_F2xz_Py_D2z_S_b;
  Double I_ERI_Fx2y_D2y_D2z_S_b = I_ERI_Gx3y_Py_D2z_S_b+ABY*I_ERI_Fx2y_Py_D2z_S_b;
  Double I_ERI_Fxyz_D2y_D2z_S_b = I_ERI_Gx2yz_Py_D2z_S_b+ABY*I_ERI_Fxyz_Py_D2z_S_b;
  Double I_ERI_Fx2z_D2y_D2z_S_b = I_ERI_Gxy2z_Py_D2z_S_b+ABY*I_ERI_Fx2z_Py_D2z_S_b;
  Double I_ERI_F3y_D2y_D2z_S_b = I_ERI_G4y_Py_D2z_S_b+ABY*I_ERI_F3y_Py_D2z_S_b;
  Double I_ERI_F2yz_D2y_D2z_S_b = I_ERI_G3yz_Py_D2z_S_b+ABY*I_ERI_F2yz_Py_D2z_S_b;
  Double I_ERI_Fy2z_D2y_D2z_S_b = I_ERI_G2y2z_Py_D2z_S_b+ABY*I_ERI_Fy2z_Py_D2z_S_b;
  Double I_ERI_F3z_D2y_D2z_S_b = I_ERI_Gy3z_Py_D2z_S_b+ABY*I_ERI_F3z_Py_D2z_S_b;
  Double I_ERI_F3x_D2z_D2z_S_b = I_ERI_G3xz_Pz_D2z_S_b+ABZ*I_ERI_F3x_Pz_D2z_S_b;
  Double I_ERI_F2xy_D2z_D2z_S_b = I_ERI_G2xyz_Pz_D2z_S_b+ABZ*I_ERI_F2xy_Pz_D2z_S_b;
  Double I_ERI_F2xz_D2z_D2z_S_b = I_ERI_G2x2z_Pz_D2z_S_b+ABZ*I_ERI_F2xz_Pz_D2z_S_b;
  Double I_ERI_Fx2y_D2z_D2z_S_b = I_ERI_Gx2yz_Pz_D2z_S_b+ABZ*I_ERI_Fx2y_Pz_D2z_S_b;
  Double I_ERI_Fxyz_D2z_D2z_S_b = I_ERI_Gxy2z_Pz_D2z_S_b+ABZ*I_ERI_Fxyz_Pz_D2z_S_b;
  Double I_ERI_Fx2z_D2z_D2z_S_b = I_ERI_Gx3z_Pz_D2z_S_b+ABZ*I_ERI_Fx2z_Pz_D2z_S_b;
  Double I_ERI_F3y_D2z_D2z_S_b = I_ERI_G3yz_Pz_D2z_S_b+ABZ*I_ERI_F3y_Pz_D2z_S_b;
  Double I_ERI_F2yz_D2z_D2z_S_b = I_ERI_G2y2z_Pz_D2z_S_b+ABZ*I_ERI_F2yz_Pz_D2z_S_b;
  Double I_ERI_Fy2z_D2z_D2z_S_b = I_ERI_Gy3z_Pz_D2z_S_b+ABZ*I_ERI_Fy2z_Pz_D2z_S_b;
  Double I_ERI_F3z_D2z_D2z_S_b = I_ERI_G4z_Pz_D2z_S_b+ABZ*I_ERI_F3z_Pz_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_F_D_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_b
   * RHS shell quartet name: SQ_ERI_D_D_D_S_b
   ************************************************************/
  Double I_ERI_D2x_F3x_D2x_S_b = I_ERI_F3x_D2x_D2x_S_b+ABX*I_ERI_D2x_D2x_D2x_S_b;
  Double I_ERI_Dxy_F3x_D2x_S_b = I_ERI_F2xy_D2x_D2x_S_b+ABX*I_ERI_Dxy_D2x_D2x_S_b;
  Double I_ERI_Dxz_F3x_D2x_S_b = I_ERI_F2xz_D2x_D2x_S_b+ABX*I_ERI_Dxz_D2x_D2x_S_b;
  Double I_ERI_D2y_F3x_D2x_S_b = I_ERI_Fx2y_D2x_D2x_S_b+ABX*I_ERI_D2y_D2x_D2x_S_b;
  Double I_ERI_Dyz_F3x_D2x_S_b = I_ERI_Fxyz_D2x_D2x_S_b+ABX*I_ERI_Dyz_D2x_D2x_S_b;
  Double I_ERI_D2z_F3x_D2x_S_b = I_ERI_Fx2z_D2x_D2x_S_b+ABX*I_ERI_D2z_D2x_D2x_S_b;
  Double I_ERI_D2x_F2xy_D2x_S_b = I_ERI_F2xy_D2x_D2x_S_b+ABY*I_ERI_D2x_D2x_D2x_S_b;
  Double I_ERI_Dxy_F2xy_D2x_S_b = I_ERI_Fx2y_D2x_D2x_S_b+ABY*I_ERI_Dxy_D2x_D2x_S_b;
  Double I_ERI_Dxz_F2xy_D2x_S_b = I_ERI_Fxyz_D2x_D2x_S_b+ABY*I_ERI_Dxz_D2x_D2x_S_b;
  Double I_ERI_D2y_F2xy_D2x_S_b = I_ERI_F3y_D2x_D2x_S_b+ABY*I_ERI_D2y_D2x_D2x_S_b;
  Double I_ERI_Dyz_F2xy_D2x_S_b = I_ERI_F2yz_D2x_D2x_S_b+ABY*I_ERI_Dyz_D2x_D2x_S_b;
  Double I_ERI_D2z_F2xy_D2x_S_b = I_ERI_Fy2z_D2x_D2x_S_b+ABY*I_ERI_D2z_D2x_D2x_S_b;
  Double I_ERI_D2x_F2xz_D2x_S_b = I_ERI_F2xz_D2x_D2x_S_b+ABZ*I_ERI_D2x_D2x_D2x_S_b;
  Double I_ERI_Dxy_F2xz_D2x_S_b = I_ERI_Fxyz_D2x_D2x_S_b+ABZ*I_ERI_Dxy_D2x_D2x_S_b;
  Double I_ERI_Dxz_F2xz_D2x_S_b = I_ERI_Fx2z_D2x_D2x_S_b+ABZ*I_ERI_Dxz_D2x_D2x_S_b;
  Double I_ERI_D2y_F2xz_D2x_S_b = I_ERI_F2yz_D2x_D2x_S_b+ABZ*I_ERI_D2y_D2x_D2x_S_b;
  Double I_ERI_Dyz_F2xz_D2x_S_b = I_ERI_Fy2z_D2x_D2x_S_b+ABZ*I_ERI_Dyz_D2x_D2x_S_b;
  Double I_ERI_D2z_F2xz_D2x_S_b = I_ERI_F3z_D2x_D2x_S_b+ABZ*I_ERI_D2z_D2x_D2x_S_b;
  Double I_ERI_D2x_Fx2y_D2x_S_b = I_ERI_F3x_D2y_D2x_S_b+ABX*I_ERI_D2x_D2y_D2x_S_b;
  Double I_ERI_Dxy_Fx2y_D2x_S_b = I_ERI_F2xy_D2y_D2x_S_b+ABX*I_ERI_Dxy_D2y_D2x_S_b;
  Double I_ERI_Dxz_Fx2y_D2x_S_b = I_ERI_F2xz_D2y_D2x_S_b+ABX*I_ERI_Dxz_D2y_D2x_S_b;
  Double I_ERI_D2y_Fx2y_D2x_S_b = I_ERI_Fx2y_D2y_D2x_S_b+ABX*I_ERI_D2y_D2y_D2x_S_b;
  Double I_ERI_Dyz_Fx2y_D2x_S_b = I_ERI_Fxyz_D2y_D2x_S_b+ABX*I_ERI_Dyz_D2y_D2x_S_b;
  Double I_ERI_D2z_Fx2y_D2x_S_b = I_ERI_Fx2z_D2y_D2x_S_b+ABX*I_ERI_D2z_D2y_D2x_S_b;
  Double I_ERI_D2x_Fxyz_D2x_S_b = I_ERI_F2xz_Dxy_D2x_S_b+ABZ*I_ERI_D2x_Dxy_D2x_S_b;
  Double I_ERI_Dxy_Fxyz_D2x_S_b = I_ERI_Fxyz_Dxy_D2x_S_b+ABZ*I_ERI_Dxy_Dxy_D2x_S_b;
  Double I_ERI_Dxz_Fxyz_D2x_S_b = I_ERI_Fx2z_Dxy_D2x_S_b+ABZ*I_ERI_Dxz_Dxy_D2x_S_b;
  Double I_ERI_D2y_Fxyz_D2x_S_b = I_ERI_F2yz_Dxy_D2x_S_b+ABZ*I_ERI_D2y_Dxy_D2x_S_b;
  Double I_ERI_Dyz_Fxyz_D2x_S_b = I_ERI_Fy2z_Dxy_D2x_S_b+ABZ*I_ERI_Dyz_Dxy_D2x_S_b;
  Double I_ERI_D2z_Fxyz_D2x_S_b = I_ERI_F3z_Dxy_D2x_S_b+ABZ*I_ERI_D2z_Dxy_D2x_S_b;
  Double I_ERI_D2x_Fx2z_D2x_S_b = I_ERI_F3x_D2z_D2x_S_b+ABX*I_ERI_D2x_D2z_D2x_S_b;
  Double I_ERI_Dxy_Fx2z_D2x_S_b = I_ERI_F2xy_D2z_D2x_S_b+ABX*I_ERI_Dxy_D2z_D2x_S_b;
  Double I_ERI_Dxz_Fx2z_D2x_S_b = I_ERI_F2xz_D2z_D2x_S_b+ABX*I_ERI_Dxz_D2z_D2x_S_b;
  Double I_ERI_D2y_Fx2z_D2x_S_b = I_ERI_Fx2y_D2z_D2x_S_b+ABX*I_ERI_D2y_D2z_D2x_S_b;
  Double I_ERI_Dyz_Fx2z_D2x_S_b = I_ERI_Fxyz_D2z_D2x_S_b+ABX*I_ERI_Dyz_D2z_D2x_S_b;
  Double I_ERI_D2z_Fx2z_D2x_S_b = I_ERI_Fx2z_D2z_D2x_S_b+ABX*I_ERI_D2z_D2z_D2x_S_b;
  Double I_ERI_D2x_F3y_D2x_S_b = I_ERI_F2xy_D2y_D2x_S_b+ABY*I_ERI_D2x_D2y_D2x_S_b;
  Double I_ERI_Dxy_F3y_D2x_S_b = I_ERI_Fx2y_D2y_D2x_S_b+ABY*I_ERI_Dxy_D2y_D2x_S_b;
  Double I_ERI_Dxz_F3y_D2x_S_b = I_ERI_Fxyz_D2y_D2x_S_b+ABY*I_ERI_Dxz_D2y_D2x_S_b;
  Double I_ERI_D2y_F3y_D2x_S_b = I_ERI_F3y_D2y_D2x_S_b+ABY*I_ERI_D2y_D2y_D2x_S_b;
  Double I_ERI_Dyz_F3y_D2x_S_b = I_ERI_F2yz_D2y_D2x_S_b+ABY*I_ERI_Dyz_D2y_D2x_S_b;
  Double I_ERI_D2z_F3y_D2x_S_b = I_ERI_Fy2z_D2y_D2x_S_b+ABY*I_ERI_D2z_D2y_D2x_S_b;
  Double I_ERI_D2x_F2yz_D2x_S_b = I_ERI_F2xz_D2y_D2x_S_b+ABZ*I_ERI_D2x_D2y_D2x_S_b;
  Double I_ERI_Dxy_F2yz_D2x_S_b = I_ERI_Fxyz_D2y_D2x_S_b+ABZ*I_ERI_Dxy_D2y_D2x_S_b;
  Double I_ERI_Dxz_F2yz_D2x_S_b = I_ERI_Fx2z_D2y_D2x_S_b+ABZ*I_ERI_Dxz_D2y_D2x_S_b;
  Double I_ERI_D2y_F2yz_D2x_S_b = I_ERI_F2yz_D2y_D2x_S_b+ABZ*I_ERI_D2y_D2y_D2x_S_b;
  Double I_ERI_Dyz_F2yz_D2x_S_b = I_ERI_Fy2z_D2y_D2x_S_b+ABZ*I_ERI_Dyz_D2y_D2x_S_b;
  Double I_ERI_D2z_F2yz_D2x_S_b = I_ERI_F3z_D2y_D2x_S_b+ABZ*I_ERI_D2z_D2y_D2x_S_b;
  Double I_ERI_D2x_Fy2z_D2x_S_b = I_ERI_F2xy_D2z_D2x_S_b+ABY*I_ERI_D2x_D2z_D2x_S_b;
  Double I_ERI_Dxy_Fy2z_D2x_S_b = I_ERI_Fx2y_D2z_D2x_S_b+ABY*I_ERI_Dxy_D2z_D2x_S_b;
  Double I_ERI_Dxz_Fy2z_D2x_S_b = I_ERI_Fxyz_D2z_D2x_S_b+ABY*I_ERI_Dxz_D2z_D2x_S_b;
  Double I_ERI_D2y_Fy2z_D2x_S_b = I_ERI_F3y_D2z_D2x_S_b+ABY*I_ERI_D2y_D2z_D2x_S_b;
  Double I_ERI_Dyz_Fy2z_D2x_S_b = I_ERI_F2yz_D2z_D2x_S_b+ABY*I_ERI_Dyz_D2z_D2x_S_b;
  Double I_ERI_D2z_Fy2z_D2x_S_b = I_ERI_Fy2z_D2z_D2x_S_b+ABY*I_ERI_D2z_D2z_D2x_S_b;
  Double I_ERI_D2x_F3z_D2x_S_b = I_ERI_F2xz_D2z_D2x_S_b+ABZ*I_ERI_D2x_D2z_D2x_S_b;
  Double I_ERI_Dxy_F3z_D2x_S_b = I_ERI_Fxyz_D2z_D2x_S_b+ABZ*I_ERI_Dxy_D2z_D2x_S_b;
  Double I_ERI_Dxz_F3z_D2x_S_b = I_ERI_Fx2z_D2z_D2x_S_b+ABZ*I_ERI_Dxz_D2z_D2x_S_b;
  Double I_ERI_D2y_F3z_D2x_S_b = I_ERI_F2yz_D2z_D2x_S_b+ABZ*I_ERI_D2y_D2z_D2x_S_b;
  Double I_ERI_Dyz_F3z_D2x_S_b = I_ERI_Fy2z_D2z_D2x_S_b+ABZ*I_ERI_Dyz_D2z_D2x_S_b;
  Double I_ERI_D2z_F3z_D2x_S_b = I_ERI_F3z_D2z_D2x_S_b+ABZ*I_ERI_D2z_D2z_D2x_S_b;
  Double I_ERI_D2x_F3x_Dxy_S_b = I_ERI_F3x_D2x_Dxy_S_b+ABX*I_ERI_D2x_D2x_Dxy_S_b;
  Double I_ERI_Dxy_F3x_Dxy_S_b = I_ERI_F2xy_D2x_Dxy_S_b+ABX*I_ERI_Dxy_D2x_Dxy_S_b;
  Double I_ERI_Dxz_F3x_Dxy_S_b = I_ERI_F2xz_D2x_Dxy_S_b+ABX*I_ERI_Dxz_D2x_Dxy_S_b;
  Double I_ERI_D2y_F3x_Dxy_S_b = I_ERI_Fx2y_D2x_Dxy_S_b+ABX*I_ERI_D2y_D2x_Dxy_S_b;
  Double I_ERI_Dyz_F3x_Dxy_S_b = I_ERI_Fxyz_D2x_Dxy_S_b+ABX*I_ERI_Dyz_D2x_Dxy_S_b;
  Double I_ERI_D2z_F3x_Dxy_S_b = I_ERI_Fx2z_D2x_Dxy_S_b+ABX*I_ERI_D2z_D2x_Dxy_S_b;
  Double I_ERI_D2x_F2xy_Dxy_S_b = I_ERI_F2xy_D2x_Dxy_S_b+ABY*I_ERI_D2x_D2x_Dxy_S_b;
  Double I_ERI_Dxy_F2xy_Dxy_S_b = I_ERI_Fx2y_D2x_Dxy_S_b+ABY*I_ERI_Dxy_D2x_Dxy_S_b;
  Double I_ERI_Dxz_F2xy_Dxy_S_b = I_ERI_Fxyz_D2x_Dxy_S_b+ABY*I_ERI_Dxz_D2x_Dxy_S_b;
  Double I_ERI_D2y_F2xy_Dxy_S_b = I_ERI_F3y_D2x_Dxy_S_b+ABY*I_ERI_D2y_D2x_Dxy_S_b;
  Double I_ERI_Dyz_F2xy_Dxy_S_b = I_ERI_F2yz_D2x_Dxy_S_b+ABY*I_ERI_Dyz_D2x_Dxy_S_b;
  Double I_ERI_D2z_F2xy_Dxy_S_b = I_ERI_Fy2z_D2x_Dxy_S_b+ABY*I_ERI_D2z_D2x_Dxy_S_b;
  Double I_ERI_D2x_F2xz_Dxy_S_b = I_ERI_F2xz_D2x_Dxy_S_b+ABZ*I_ERI_D2x_D2x_Dxy_S_b;
  Double I_ERI_Dxy_F2xz_Dxy_S_b = I_ERI_Fxyz_D2x_Dxy_S_b+ABZ*I_ERI_Dxy_D2x_Dxy_S_b;
  Double I_ERI_Dxz_F2xz_Dxy_S_b = I_ERI_Fx2z_D2x_Dxy_S_b+ABZ*I_ERI_Dxz_D2x_Dxy_S_b;
  Double I_ERI_D2y_F2xz_Dxy_S_b = I_ERI_F2yz_D2x_Dxy_S_b+ABZ*I_ERI_D2y_D2x_Dxy_S_b;
  Double I_ERI_Dyz_F2xz_Dxy_S_b = I_ERI_Fy2z_D2x_Dxy_S_b+ABZ*I_ERI_Dyz_D2x_Dxy_S_b;
  Double I_ERI_D2z_F2xz_Dxy_S_b = I_ERI_F3z_D2x_Dxy_S_b+ABZ*I_ERI_D2z_D2x_Dxy_S_b;
  Double I_ERI_D2x_Fx2y_Dxy_S_b = I_ERI_F3x_D2y_Dxy_S_b+ABX*I_ERI_D2x_D2y_Dxy_S_b;
  Double I_ERI_Dxy_Fx2y_Dxy_S_b = I_ERI_F2xy_D2y_Dxy_S_b+ABX*I_ERI_Dxy_D2y_Dxy_S_b;
  Double I_ERI_Dxz_Fx2y_Dxy_S_b = I_ERI_F2xz_D2y_Dxy_S_b+ABX*I_ERI_Dxz_D2y_Dxy_S_b;
  Double I_ERI_D2y_Fx2y_Dxy_S_b = I_ERI_Fx2y_D2y_Dxy_S_b+ABX*I_ERI_D2y_D2y_Dxy_S_b;
  Double I_ERI_Dyz_Fx2y_Dxy_S_b = I_ERI_Fxyz_D2y_Dxy_S_b+ABX*I_ERI_Dyz_D2y_Dxy_S_b;
  Double I_ERI_D2z_Fx2y_Dxy_S_b = I_ERI_Fx2z_D2y_Dxy_S_b+ABX*I_ERI_D2z_D2y_Dxy_S_b;
  Double I_ERI_D2x_Fxyz_Dxy_S_b = I_ERI_F2xz_Dxy_Dxy_S_b+ABZ*I_ERI_D2x_Dxy_Dxy_S_b;
  Double I_ERI_Dxy_Fxyz_Dxy_S_b = I_ERI_Fxyz_Dxy_Dxy_S_b+ABZ*I_ERI_Dxy_Dxy_Dxy_S_b;
  Double I_ERI_Dxz_Fxyz_Dxy_S_b = I_ERI_Fx2z_Dxy_Dxy_S_b+ABZ*I_ERI_Dxz_Dxy_Dxy_S_b;
  Double I_ERI_D2y_Fxyz_Dxy_S_b = I_ERI_F2yz_Dxy_Dxy_S_b+ABZ*I_ERI_D2y_Dxy_Dxy_S_b;
  Double I_ERI_Dyz_Fxyz_Dxy_S_b = I_ERI_Fy2z_Dxy_Dxy_S_b+ABZ*I_ERI_Dyz_Dxy_Dxy_S_b;
  Double I_ERI_D2z_Fxyz_Dxy_S_b = I_ERI_F3z_Dxy_Dxy_S_b+ABZ*I_ERI_D2z_Dxy_Dxy_S_b;
  Double I_ERI_D2x_Fx2z_Dxy_S_b = I_ERI_F3x_D2z_Dxy_S_b+ABX*I_ERI_D2x_D2z_Dxy_S_b;
  Double I_ERI_Dxy_Fx2z_Dxy_S_b = I_ERI_F2xy_D2z_Dxy_S_b+ABX*I_ERI_Dxy_D2z_Dxy_S_b;
  Double I_ERI_Dxz_Fx2z_Dxy_S_b = I_ERI_F2xz_D2z_Dxy_S_b+ABX*I_ERI_Dxz_D2z_Dxy_S_b;
  Double I_ERI_D2y_Fx2z_Dxy_S_b = I_ERI_Fx2y_D2z_Dxy_S_b+ABX*I_ERI_D2y_D2z_Dxy_S_b;
  Double I_ERI_Dyz_Fx2z_Dxy_S_b = I_ERI_Fxyz_D2z_Dxy_S_b+ABX*I_ERI_Dyz_D2z_Dxy_S_b;
  Double I_ERI_D2z_Fx2z_Dxy_S_b = I_ERI_Fx2z_D2z_Dxy_S_b+ABX*I_ERI_D2z_D2z_Dxy_S_b;
  Double I_ERI_D2x_F3y_Dxy_S_b = I_ERI_F2xy_D2y_Dxy_S_b+ABY*I_ERI_D2x_D2y_Dxy_S_b;
  Double I_ERI_Dxy_F3y_Dxy_S_b = I_ERI_Fx2y_D2y_Dxy_S_b+ABY*I_ERI_Dxy_D2y_Dxy_S_b;
  Double I_ERI_Dxz_F3y_Dxy_S_b = I_ERI_Fxyz_D2y_Dxy_S_b+ABY*I_ERI_Dxz_D2y_Dxy_S_b;
  Double I_ERI_D2y_F3y_Dxy_S_b = I_ERI_F3y_D2y_Dxy_S_b+ABY*I_ERI_D2y_D2y_Dxy_S_b;
  Double I_ERI_Dyz_F3y_Dxy_S_b = I_ERI_F2yz_D2y_Dxy_S_b+ABY*I_ERI_Dyz_D2y_Dxy_S_b;
  Double I_ERI_D2z_F3y_Dxy_S_b = I_ERI_Fy2z_D2y_Dxy_S_b+ABY*I_ERI_D2z_D2y_Dxy_S_b;
  Double I_ERI_D2x_F2yz_Dxy_S_b = I_ERI_F2xz_D2y_Dxy_S_b+ABZ*I_ERI_D2x_D2y_Dxy_S_b;
  Double I_ERI_Dxy_F2yz_Dxy_S_b = I_ERI_Fxyz_D2y_Dxy_S_b+ABZ*I_ERI_Dxy_D2y_Dxy_S_b;
  Double I_ERI_Dxz_F2yz_Dxy_S_b = I_ERI_Fx2z_D2y_Dxy_S_b+ABZ*I_ERI_Dxz_D2y_Dxy_S_b;
  Double I_ERI_D2y_F2yz_Dxy_S_b = I_ERI_F2yz_D2y_Dxy_S_b+ABZ*I_ERI_D2y_D2y_Dxy_S_b;
  Double I_ERI_Dyz_F2yz_Dxy_S_b = I_ERI_Fy2z_D2y_Dxy_S_b+ABZ*I_ERI_Dyz_D2y_Dxy_S_b;
  Double I_ERI_D2z_F2yz_Dxy_S_b = I_ERI_F3z_D2y_Dxy_S_b+ABZ*I_ERI_D2z_D2y_Dxy_S_b;
  Double I_ERI_D2x_Fy2z_Dxy_S_b = I_ERI_F2xy_D2z_Dxy_S_b+ABY*I_ERI_D2x_D2z_Dxy_S_b;
  Double I_ERI_Dxy_Fy2z_Dxy_S_b = I_ERI_Fx2y_D2z_Dxy_S_b+ABY*I_ERI_Dxy_D2z_Dxy_S_b;
  Double I_ERI_Dxz_Fy2z_Dxy_S_b = I_ERI_Fxyz_D2z_Dxy_S_b+ABY*I_ERI_Dxz_D2z_Dxy_S_b;
  Double I_ERI_D2y_Fy2z_Dxy_S_b = I_ERI_F3y_D2z_Dxy_S_b+ABY*I_ERI_D2y_D2z_Dxy_S_b;
  Double I_ERI_Dyz_Fy2z_Dxy_S_b = I_ERI_F2yz_D2z_Dxy_S_b+ABY*I_ERI_Dyz_D2z_Dxy_S_b;
  Double I_ERI_D2z_Fy2z_Dxy_S_b = I_ERI_Fy2z_D2z_Dxy_S_b+ABY*I_ERI_D2z_D2z_Dxy_S_b;
  Double I_ERI_D2x_F3z_Dxy_S_b = I_ERI_F2xz_D2z_Dxy_S_b+ABZ*I_ERI_D2x_D2z_Dxy_S_b;
  Double I_ERI_Dxy_F3z_Dxy_S_b = I_ERI_Fxyz_D2z_Dxy_S_b+ABZ*I_ERI_Dxy_D2z_Dxy_S_b;
  Double I_ERI_Dxz_F3z_Dxy_S_b = I_ERI_Fx2z_D2z_Dxy_S_b+ABZ*I_ERI_Dxz_D2z_Dxy_S_b;
  Double I_ERI_D2y_F3z_Dxy_S_b = I_ERI_F2yz_D2z_Dxy_S_b+ABZ*I_ERI_D2y_D2z_Dxy_S_b;
  Double I_ERI_Dyz_F3z_Dxy_S_b = I_ERI_Fy2z_D2z_Dxy_S_b+ABZ*I_ERI_Dyz_D2z_Dxy_S_b;
  Double I_ERI_D2z_F3z_Dxy_S_b = I_ERI_F3z_D2z_Dxy_S_b+ABZ*I_ERI_D2z_D2z_Dxy_S_b;
  Double I_ERI_D2x_F3x_Dxz_S_b = I_ERI_F3x_D2x_Dxz_S_b+ABX*I_ERI_D2x_D2x_Dxz_S_b;
  Double I_ERI_Dxy_F3x_Dxz_S_b = I_ERI_F2xy_D2x_Dxz_S_b+ABX*I_ERI_Dxy_D2x_Dxz_S_b;
  Double I_ERI_Dxz_F3x_Dxz_S_b = I_ERI_F2xz_D2x_Dxz_S_b+ABX*I_ERI_Dxz_D2x_Dxz_S_b;
  Double I_ERI_D2y_F3x_Dxz_S_b = I_ERI_Fx2y_D2x_Dxz_S_b+ABX*I_ERI_D2y_D2x_Dxz_S_b;
  Double I_ERI_Dyz_F3x_Dxz_S_b = I_ERI_Fxyz_D2x_Dxz_S_b+ABX*I_ERI_Dyz_D2x_Dxz_S_b;
  Double I_ERI_D2z_F3x_Dxz_S_b = I_ERI_Fx2z_D2x_Dxz_S_b+ABX*I_ERI_D2z_D2x_Dxz_S_b;
  Double I_ERI_D2x_F2xy_Dxz_S_b = I_ERI_F2xy_D2x_Dxz_S_b+ABY*I_ERI_D2x_D2x_Dxz_S_b;
  Double I_ERI_Dxy_F2xy_Dxz_S_b = I_ERI_Fx2y_D2x_Dxz_S_b+ABY*I_ERI_Dxy_D2x_Dxz_S_b;
  Double I_ERI_Dxz_F2xy_Dxz_S_b = I_ERI_Fxyz_D2x_Dxz_S_b+ABY*I_ERI_Dxz_D2x_Dxz_S_b;
  Double I_ERI_D2y_F2xy_Dxz_S_b = I_ERI_F3y_D2x_Dxz_S_b+ABY*I_ERI_D2y_D2x_Dxz_S_b;
  Double I_ERI_Dyz_F2xy_Dxz_S_b = I_ERI_F2yz_D2x_Dxz_S_b+ABY*I_ERI_Dyz_D2x_Dxz_S_b;
  Double I_ERI_D2z_F2xy_Dxz_S_b = I_ERI_Fy2z_D2x_Dxz_S_b+ABY*I_ERI_D2z_D2x_Dxz_S_b;
  Double I_ERI_D2x_F2xz_Dxz_S_b = I_ERI_F2xz_D2x_Dxz_S_b+ABZ*I_ERI_D2x_D2x_Dxz_S_b;
  Double I_ERI_Dxy_F2xz_Dxz_S_b = I_ERI_Fxyz_D2x_Dxz_S_b+ABZ*I_ERI_Dxy_D2x_Dxz_S_b;
  Double I_ERI_Dxz_F2xz_Dxz_S_b = I_ERI_Fx2z_D2x_Dxz_S_b+ABZ*I_ERI_Dxz_D2x_Dxz_S_b;
  Double I_ERI_D2y_F2xz_Dxz_S_b = I_ERI_F2yz_D2x_Dxz_S_b+ABZ*I_ERI_D2y_D2x_Dxz_S_b;
  Double I_ERI_Dyz_F2xz_Dxz_S_b = I_ERI_Fy2z_D2x_Dxz_S_b+ABZ*I_ERI_Dyz_D2x_Dxz_S_b;
  Double I_ERI_D2z_F2xz_Dxz_S_b = I_ERI_F3z_D2x_Dxz_S_b+ABZ*I_ERI_D2z_D2x_Dxz_S_b;
  Double I_ERI_D2x_Fx2y_Dxz_S_b = I_ERI_F3x_D2y_Dxz_S_b+ABX*I_ERI_D2x_D2y_Dxz_S_b;
  Double I_ERI_Dxy_Fx2y_Dxz_S_b = I_ERI_F2xy_D2y_Dxz_S_b+ABX*I_ERI_Dxy_D2y_Dxz_S_b;
  Double I_ERI_Dxz_Fx2y_Dxz_S_b = I_ERI_F2xz_D2y_Dxz_S_b+ABX*I_ERI_Dxz_D2y_Dxz_S_b;
  Double I_ERI_D2y_Fx2y_Dxz_S_b = I_ERI_Fx2y_D2y_Dxz_S_b+ABX*I_ERI_D2y_D2y_Dxz_S_b;
  Double I_ERI_Dyz_Fx2y_Dxz_S_b = I_ERI_Fxyz_D2y_Dxz_S_b+ABX*I_ERI_Dyz_D2y_Dxz_S_b;
  Double I_ERI_D2z_Fx2y_Dxz_S_b = I_ERI_Fx2z_D2y_Dxz_S_b+ABX*I_ERI_D2z_D2y_Dxz_S_b;
  Double I_ERI_D2x_Fxyz_Dxz_S_b = I_ERI_F2xz_Dxy_Dxz_S_b+ABZ*I_ERI_D2x_Dxy_Dxz_S_b;
  Double I_ERI_Dxy_Fxyz_Dxz_S_b = I_ERI_Fxyz_Dxy_Dxz_S_b+ABZ*I_ERI_Dxy_Dxy_Dxz_S_b;
  Double I_ERI_Dxz_Fxyz_Dxz_S_b = I_ERI_Fx2z_Dxy_Dxz_S_b+ABZ*I_ERI_Dxz_Dxy_Dxz_S_b;
  Double I_ERI_D2y_Fxyz_Dxz_S_b = I_ERI_F2yz_Dxy_Dxz_S_b+ABZ*I_ERI_D2y_Dxy_Dxz_S_b;
  Double I_ERI_Dyz_Fxyz_Dxz_S_b = I_ERI_Fy2z_Dxy_Dxz_S_b+ABZ*I_ERI_Dyz_Dxy_Dxz_S_b;
  Double I_ERI_D2z_Fxyz_Dxz_S_b = I_ERI_F3z_Dxy_Dxz_S_b+ABZ*I_ERI_D2z_Dxy_Dxz_S_b;
  Double I_ERI_D2x_Fx2z_Dxz_S_b = I_ERI_F3x_D2z_Dxz_S_b+ABX*I_ERI_D2x_D2z_Dxz_S_b;
  Double I_ERI_Dxy_Fx2z_Dxz_S_b = I_ERI_F2xy_D2z_Dxz_S_b+ABX*I_ERI_Dxy_D2z_Dxz_S_b;
  Double I_ERI_Dxz_Fx2z_Dxz_S_b = I_ERI_F2xz_D2z_Dxz_S_b+ABX*I_ERI_Dxz_D2z_Dxz_S_b;
  Double I_ERI_D2y_Fx2z_Dxz_S_b = I_ERI_Fx2y_D2z_Dxz_S_b+ABX*I_ERI_D2y_D2z_Dxz_S_b;
  Double I_ERI_Dyz_Fx2z_Dxz_S_b = I_ERI_Fxyz_D2z_Dxz_S_b+ABX*I_ERI_Dyz_D2z_Dxz_S_b;
  Double I_ERI_D2z_Fx2z_Dxz_S_b = I_ERI_Fx2z_D2z_Dxz_S_b+ABX*I_ERI_D2z_D2z_Dxz_S_b;
  Double I_ERI_D2x_F3y_Dxz_S_b = I_ERI_F2xy_D2y_Dxz_S_b+ABY*I_ERI_D2x_D2y_Dxz_S_b;
  Double I_ERI_Dxy_F3y_Dxz_S_b = I_ERI_Fx2y_D2y_Dxz_S_b+ABY*I_ERI_Dxy_D2y_Dxz_S_b;
  Double I_ERI_Dxz_F3y_Dxz_S_b = I_ERI_Fxyz_D2y_Dxz_S_b+ABY*I_ERI_Dxz_D2y_Dxz_S_b;
  Double I_ERI_D2y_F3y_Dxz_S_b = I_ERI_F3y_D2y_Dxz_S_b+ABY*I_ERI_D2y_D2y_Dxz_S_b;
  Double I_ERI_Dyz_F3y_Dxz_S_b = I_ERI_F2yz_D2y_Dxz_S_b+ABY*I_ERI_Dyz_D2y_Dxz_S_b;
  Double I_ERI_D2z_F3y_Dxz_S_b = I_ERI_Fy2z_D2y_Dxz_S_b+ABY*I_ERI_D2z_D2y_Dxz_S_b;
  Double I_ERI_D2x_F2yz_Dxz_S_b = I_ERI_F2xz_D2y_Dxz_S_b+ABZ*I_ERI_D2x_D2y_Dxz_S_b;
  Double I_ERI_Dxy_F2yz_Dxz_S_b = I_ERI_Fxyz_D2y_Dxz_S_b+ABZ*I_ERI_Dxy_D2y_Dxz_S_b;
  Double I_ERI_Dxz_F2yz_Dxz_S_b = I_ERI_Fx2z_D2y_Dxz_S_b+ABZ*I_ERI_Dxz_D2y_Dxz_S_b;
  Double I_ERI_D2y_F2yz_Dxz_S_b = I_ERI_F2yz_D2y_Dxz_S_b+ABZ*I_ERI_D2y_D2y_Dxz_S_b;
  Double I_ERI_Dyz_F2yz_Dxz_S_b = I_ERI_Fy2z_D2y_Dxz_S_b+ABZ*I_ERI_Dyz_D2y_Dxz_S_b;
  Double I_ERI_D2z_F2yz_Dxz_S_b = I_ERI_F3z_D2y_Dxz_S_b+ABZ*I_ERI_D2z_D2y_Dxz_S_b;
  Double I_ERI_D2x_Fy2z_Dxz_S_b = I_ERI_F2xy_D2z_Dxz_S_b+ABY*I_ERI_D2x_D2z_Dxz_S_b;
  Double I_ERI_Dxy_Fy2z_Dxz_S_b = I_ERI_Fx2y_D2z_Dxz_S_b+ABY*I_ERI_Dxy_D2z_Dxz_S_b;
  Double I_ERI_Dxz_Fy2z_Dxz_S_b = I_ERI_Fxyz_D2z_Dxz_S_b+ABY*I_ERI_Dxz_D2z_Dxz_S_b;
  Double I_ERI_D2y_Fy2z_Dxz_S_b = I_ERI_F3y_D2z_Dxz_S_b+ABY*I_ERI_D2y_D2z_Dxz_S_b;
  Double I_ERI_Dyz_Fy2z_Dxz_S_b = I_ERI_F2yz_D2z_Dxz_S_b+ABY*I_ERI_Dyz_D2z_Dxz_S_b;
  Double I_ERI_D2z_Fy2z_Dxz_S_b = I_ERI_Fy2z_D2z_Dxz_S_b+ABY*I_ERI_D2z_D2z_Dxz_S_b;
  Double I_ERI_D2x_F3z_Dxz_S_b = I_ERI_F2xz_D2z_Dxz_S_b+ABZ*I_ERI_D2x_D2z_Dxz_S_b;
  Double I_ERI_Dxy_F3z_Dxz_S_b = I_ERI_Fxyz_D2z_Dxz_S_b+ABZ*I_ERI_Dxy_D2z_Dxz_S_b;
  Double I_ERI_Dxz_F3z_Dxz_S_b = I_ERI_Fx2z_D2z_Dxz_S_b+ABZ*I_ERI_Dxz_D2z_Dxz_S_b;
  Double I_ERI_D2y_F3z_Dxz_S_b = I_ERI_F2yz_D2z_Dxz_S_b+ABZ*I_ERI_D2y_D2z_Dxz_S_b;
  Double I_ERI_Dyz_F3z_Dxz_S_b = I_ERI_Fy2z_D2z_Dxz_S_b+ABZ*I_ERI_Dyz_D2z_Dxz_S_b;
  Double I_ERI_D2z_F3z_Dxz_S_b = I_ERI_F3z_D2z_Dxz_S_b+ABZ*I_ERI_D2z_D2z_Dxz_S_b;
  Double I_ERI_D2x_F3x_D2y_S_b = I_ERI_F3x_D2x_D2y_S_b+ABX*I_ERI_D2x_D2x_D2y_S_b;
  Double I_ERI_Dxy_F3x_D2y_S_b = I_ERI_F2xy_D2x_D2y_S_b+ABX*I_ERI_Dxy_D2x_D2y_S_b;
  Double I_ERI_Dxz_F3x_D2y_S_b = I_ERI_F2xz_D2x_D2y_S_b+ABX*I_ERI_Dxz_D2x_D2y_S_b;
  Double I_ERI_D2y_F3x_D2y_S_b = I_ERI_Fx2y_D2x_D2y_S_b+ABX*I_ERI_D2y_D2x_D2y_S_b;
  Double I_ERI_Dyz_F3x_D2y_S_b = I_ERI_Fxyz_D2x_D2y_S_b+ABX*I_ERI_Dyz_D2x_D2y_S_b;
  Double I_ERI_D2z_F3x_D2y_S_b = I_ERI_Fx2z_D2x_D2y_S_b+ABX*I_ERI_D2z_D2x_D2y_S_b;
  Double I_ERI_D2x_F2xy_D2y_S_b = I_ERI_F2xy_D2x_D2y_S_b+ABY*I_ERI_D2x_D2x_D2y_S_b;
  Double I_ERI_Dxy_F2xy_D2y_S_b = I_ERI_Fx2y_D2x_D2y_S_b+ABY*I_ERI_Dxy_D2x_D2y_S_b;
  Double I_ERI_Dxz_F2xy_D2y_S_b = I_ERI_Fxyz_D2x_D2y_S_b+ABY*I_ERI_Dxz_D2x_D2y_S_b;
  Double I_ERI_D2y_F2xy_D2y_S_b = I_ERI_F3y_D2x_D2y_S_b+ABY*I_ERI_D2y_D2x_D2y_S_b;
  Double I_ERI_Dyz_F2xy_D2y_S_b = I_ERI_F2yz_D2x_D2y_S_b+ABY*I_ERI_Dyz_D2x_D2y_S_b;
  Double I_ERI_D2z_F2xy_D2y_S_b = I_ERI_Fy2z_D2x_D2y_S_b+ABY*I_ERI_D2z_D2x_D2y_S_b;
  Double I_ERI_D2x_F2xz_D2y_S_b = I_ERI_F2xz_D2x_D2y_S_b+ABZ*I_ERI_D2x_D2x_D2y_S_b;
  Double I_ERI_Dxy_F2xz_D2y_S_b = I_ERI_Fxyz_D2x_D2y_S_b+ABZ*I_ERI_Dxy_D2x_D2y_S_b;
  Double I_ERI_Dxz_F2xz_D2y_S_b = I_ERI_Fx2z_D2x_D2y_S_b+ABZ*I_ERI_Dxz_D2x_D2y_S_b;
  Double I_ERI_D2y_F2xz_D2y_S_b = I_ERI_F2yz_D2x_D2y_S_b+ABZ*I_ERI_D2y_D2x_D2y_S_b;
  Double I_ERI_Dyz_F2xz_D2y_S_b = I_ERI_Fy2z_D2x_D2y_S_b+ABZ*I_ERI_Dyz_D2x_D2y_S_b;
  Double I_ERI_D2z_F2xz_D2y_S_b = I_ERI_F3z_D2x_D2y_S_b+ABZ*I_ERI_D2z_D2x_D2y_S_b;
  Double I_ERI_D2x_Fx2y_D2y_S_b = I_ERI_F3x_D2y_D2y_S_b+ABX*I_ERI_D2x_D2y_D2y_S_b;
  Double I_ERI_Dxy_Fx2y_D2y_S_b = I_ERI_F2xy_D2y_D2y_S_b+ABX*I_ERI_Dxy_D2y_D2y_S_b;
  Double I_ERI_Dxz_Fx2y_D2y_S_b = I_ERI_F2xz_D2y_D2y_S_b+ABX*I_ERI_Dxz_D2y_D2y_S_b;
  Double I_ERI_D2y_Fx2y_D2y_S_b = I_ERI_Fx2y_D2y_D2y_S_b+ABX*I_ERI_D2y_D2y_D2y_S_b;
  Double I_ERI_Dyz_Fx2y_D2y_S_b = I_ERI_Fxyz_D2y_D2y_S_b+ABX*I_ERI_Dyz_D2y_D2y_S_b;
  Double I_ERI_D2z_Fx2y_D2y_S_b = I_ERI_Fx2z_D2y_D2y_S_b+ABX*I_ERI_D2z_D2y_D2y_S_b;
  Double I_ERI_D2x_Fxyz_D2y_S_b = I_ERI_F2xz_Dxy_D2y_S_b+ABZ*I_ERI_D2x_Dxy_D2y_S_b;
  Double I_ERI_Dxy_Fxyz_D2y_S_b = I_ERI_Fxyz_Dxy_D2y_S_b+ABZ*I_ERI_Dxy_Dxy_D2y_S_b;
  Double I_ERI_Dxz_Fxyz_D2y_S_b = I_ERI_Fx2z_Dxy_D2y_S_b+ABZ*I_ERI_Dxz_Dxy_D2y_S_b;
  Double I_ERI_D2y_Fxyz_D2y_S_b = I_ERI_F2yz_Dxy_D2y_S_b+ABZ*I_ERI_D2y_Dxy_D2y_S_b;
  Double I_ERI_Dyz_Fxyz_D2y_S_b = I_ERI_Fy2z_Dxy_D2y_S_b+ABZ*I_ERI_Dyz_Dxy_D2y_S_b;
  Double I_ERI_D2z_Fxyz_D2y_S_b = I_ERI_F3z_Dxy_D2y_S_b+ABZ*I_ERI_D2z_Dxy_D2y_S_b;
  Double I_ERI_D2x_Fx2z_D2y_S_b = I_ERI_F3x_D2z_D2y_S_b+ABX*I_ERI_D2x_D2z_D2y_S_b;
  Double I_ERI_Dxy_Fx2z_D2y_S_b = I_ERI_F2xy_D2z_D2y_S_b+ABX*I_ERI_Dxy_D2z_D2y_S_b;
  Double I_ERI_Dxz_Fx2z_D2y_S_b = I_ERI_F2xz_D2z_D2y_S_b+ABX*I_ERI_Dxz_D2z_D2y_S_b;
  Double I_ERI_D2y_Fx2z_D2y_S_b = I_ERI_Fx2y_D2z_D2y_S_b+ABX*I_ERI_D2y_D2z_D2y_S_b;
  Double I_ERI_Dyz_Fx2z_D2y_S_b = I_ERI_Fxyz_D2z_D2y_S_b+ABX*I_ERI_Dyz_D2z_D2y_S_b;
  Double I_ERI_D2z_Fx2z_D2y_S_b = I_ERI_Fx2z_D2z_D2y_S_b+ABX*I_ERI_D2z_D2z_D2y_S_b;
  Double I_ERI_D2x_F3y_D2y_S_b = I_ERI_F2xy_D2y_D2y_S_b+ABY*I_ERI_D2x_D2y_D2y_S_b;
  Double I_ERI_Dxy_F3y_D2y_S_b = I_ERI_Fx2y_D2y_D2y_S_b+ABY*I_ERI_Dxy_D2y_D2y_S_b;
  Double I_ERI_Dxz_F3y_D2y_S_b = I_ERI_Fxyz_D2y_D2y_S_b+ABY*I_ERI_Dxz_D2y_D2y_S_b;
  Double I_ERI_D2y_F3y_D2y_S_b = I_ERI_F3y_D2y_D2y_S_b+ABY*I_ERI_D2y_D2y_D2y_S_b;
  Double I_ERI_Dyz_F3y_D2y_S_b = I_ERI_F2yz_D2y_D2y_S_b+ABY*I_ERI_Dyz_D2y_D2y_S_b;
  Double I_ERI_D2z_F3y_D2y_S_b = I_ERI_Fy2z_D2y_D2y_S_b+ABY*I_ERI_D2z_D2y_D2y_S_b;
  Double I_ERI_D2x_F2yz_D2y_S_b = I_ERI_F2xz_D2y_D2y_S_b+ABZ*I_ERI_D2x_D2y_D2y_S_b;
  Double I_ERI_Dxy_F2yz_D2y_S_b = I_ERI_Fxyz_D2y_D2y_S_b+ABZ*I_ERI_Dxy_D2y_D2y_S_b;
  Double I_ERI_Dxz_F2yz_D2y_S_b = I_ERI_Fx2z_D2y_D2y_S_b+ABZ*I_ERI_Dxz_D2y_D2y_S_b;
  Double I_ERI_D2y_F2yz_D2y_S_b = I_ERI_F2yz_D2y_D2y_S_b+ABZ*I_ERI_D2y_D2y_D2y_S_b;
  Double I_ERI_Dyz_F2yz_D2y_S_b = I_ERI_Fy2z_D2y_D2y_S_b+ABZ*I_ERI_Dyz_D2y_D2y_S_b;
  Double I_ERI_D2z_F2yz_D2y_S_b = I_ERI_F3z_D2y_D2y_S_b+ABZ*I_ERI_D2z_D2y_D2y_S_b;
  Double I_ERI_D2x_Fy2z_D2y_S_b = I_ERI_F2xy_D2z_D2y_S_b+ABY*I_ERI_D2x_D2z_D2y_S_b;
  Double I_ERI_Dxy_Fy2z_D2y_S_b = I_ERI_Fx2y_D2z_D2y_S_b+ABY*I_ERI_Dxy_D2z_D2y_S_b;
  Double I_ERI_Dxz_Fy2z_D2y_S_b = I_ERI_Fxyz_D2z_D2y_S_b+ABY*I_ERI_Dxz_D2z_D2y_S_b;
  Double I_ERI_D2y_Fy2z_D2y_S_b = I_ERI_F3y_D2z_D2y_S_b+ABY*I_ERI_D2y_D2z_D2y_S_b;
  Double I_ERI_Dyz_Fy2z_D2y_S_b = I_ERI_F2yz_D2z_D2y_S_b+ABY*I_ERI_Dyz_D2z_D2y_S_b;
  Double I_ERI_D2z_Fy2z_D2y_S_b = I_ERI_Fy2z_D2z_D2y_S_b+ABY*I_ERI_D2z_D2z_D2y_S_b;
  Double I_ERI_D2x_F3z_D2y_S_b = I_ERI_F2xz_D2z_D2y_S_b+ABZ*I_ERI_D2x_D2z_D2y_S_b;
  Double I_ERI_Dxy_F3z_D2y_S_b = I_ERI_Fxyz_D2z_D2y_S_b+ABZ*I_ERI_Dxy_D2z_D2y_S_b;
  Double I_ERI_Dxz_F3z_D2y_S_b = I_ERI_Fx2z_D2z_D2y_S_b+ABZ*I_ERI_Dxz_D2z_D2y_S_b;
  Double I_ERI_D2y_F3z_D2y_S_b = I_ERI_F2yz_D2z_D2y_S_b+ABZ*I_ERI_D2y_D2z_D2y_S_b;
  Double I_ERI_Dyz_F3z_D2y_S_b = I_ERI_Fy2z_D2z_D2y_S_b+ABZ*I_ERI_Dyz_D2z_D2y_S_b;
  Double I_ERI_D2z_F3z_D2y_S_b = I_ERI_F3z_D2z_D2y_S_b+ABZ*I_ERI_D2z_D2z_D2y_S_b;
  Double I_ERI_D2x_F3x_Dyz_S_b = I_ERI_F3x_D2x_Dyz_S_b+ABX*I_ERI_D2x_D2x_Dyz_S_b;
  Double I_ERI_Dxy_F3x_Dyz_S_b = I_ERI_F2xy_D2x_Dyz_S_b+ABX*I_ERI_Dxy_D2x_Dyz_S_b;
  Double I_ERI_Dxz_F3x_Dyz_S_b = I_ERI_F2xz_D2x_Dyz_S_b+ABX*I_ERI_Dxz_D2x_Dyz_S_b;
  Double I_ERI_D2y_F3x_Dyz_S_b = I_ERI_Fx2y_D2x_Dyz_S_b+ABX*I_ERI_D2y_D2x_Dyz_S_b;
  Double I_ERI_Dyz_F3x_Dyz_S_b = I_ERI_Fxyz_D2x_Dyz_S_b+ABX*I_ERI_Dyz_D2x_Dyz_S_b;
  Double I_ERI_D2z_F3x_Dyz_S_b = I_ERI_Fx2z_D2x_Dyz_S_b+ABX*I_ERI_D2z_D2x_Dyz_S_b;
  Double I_ERI_D2x_F2xy_Dyz_S_b = I_ERI_F2xy_D2x_Dyz_S_b+ABY*I_ERI_D2x_D2x_Dyz_S_b;
  Double I_ERI_Dxy_F2xy_Dyz_S_b = I_ERI_Fx2y_D2x_Dyz_S_b+ABY*I_ERI_Dxy_D2x_Dyz_S_b;
  Double I_ERI_Dxz_F2xy_Dyz_S_b = I_ERI_Fxyz_D2x_Dyz_S_b+ABY*I_ERI_Dxz_D2x_Dyz_S_b;
  Double I_ERI_D2y_F2xy_Dyz_S_b = I_ERI_F3y_D2x_Dyz_S_b+ABY*I_ERI_D2y_D2x_Dyz_S_b;
  Double I_ERI_Dyz_F2xy_Dyz_S_b = I_ERI_F2yz_D2x_Dyz_S_b+ABY*I_ERI_Dyz_D2x_Dyz_S_b;
  Double I_ERI_D2z_F2xy_Dyz_S_b = I_ERI_Fy2z_D2x_Dyz_S_b+ABY*I_ERI_D2z_D2x_Dyz_S_b;
  Double I_ERI_D2x_F2xz_Dyz_S_b = I_ERI_F2xz_D2x_Dyz_S_b+ABZ*I_ERI_D2x_D2x_Dyz_S_b;
  Double I_ERI_Dxy_F2xz_Dyz_S_b = I_ERI_Fxyz_D2x_Dyz_S_b+ABZ*I_ERI_Dxy_D2x_Dyz_S_b;
  Double I_ERI_Dxz_F2xz_Dyz_S_b = I_ERI_Fx2z_D2x_Dyz_S_b+ABZ*I_ERI_Dxz_D2x_Dyz_S_b;
  Double I_ERI_D2y_F2xz_Dyz_S_b = I_ERI_F2yz_D2x_Dyz_S_b+ABZ*I_ERI_D2y_D2x_Dyz_S_b;
  Double I_ERI_Dyz_F2xz_Dyz_S_b = I_ERI_Fy2z_D2x_Dyz_S_b+ABZ*I_ERI_Dyz_D2x_Dyz_S_b;
  Double I_ERI_D2z_F2xz_Dyz_S_b = I_ERI_F3z_D2x_Dyz_S_b+ABZ*I_ERI_D2z_D2x_Dyz_S_b;
  Double I_ERI_D2x_Fx2y_Dyz_S_b = I_ERI_F3x_D2y_Dyz_S_b+ABX*I_ERI_D2x_D2y_Dyz_S_b;
  Double I_ERI_Dxy_Fx2y_Dyz_S_b = I_ERI_F2xy_D2y_Dyz_S_b+ABX*I_ERI_Dxy_D2y_Dyz_S_b;
  Double I_ERI_Dxz_Fx2y_Dyz_S_b = I_ERI_F2xz_D2y_Dyz_S_b+ABX*I_ERI_Dxz_D2y_Dyz_S_b;
  Double I_ERI_D2y_Fx2y_Dyz_S_b = I_ERI_Fx2y_D2y_Dyz_S_b+ABX*I_ERI_D2y_D2y_Dyz_S_b;
  Double I_ERI_Dyz_Fx2y_Dyz_S_b = I_ERI_Fxyz_D2y_Dyz_S_b+ABX*I_ERI_Dyz_D2y_Dyz_S_b;
  Double I_ERI_D2z_Fx2y_Dyz_S_b = I_ERI_Fx2z_D2y_Dyz_S_b+ABX*I_ERI_D2z_D2y_Dyz_S_b;
  Double I_ERI_D2x_Fxyz_Dyz_S_b = I_ERI_F2xz_Dxy_Dyz_S_b+ABZ*I_ERI_D2x_Dxy_Dyz_S_b;
  Double I_ERI_Dxy_Fxyz_Dyz_S_b = I_ERI_Fxyz_Dxy_Dyz_S_b+ABZ*I_ERI_Dxy_Dxy_Dyz_S_b;
  Double I_ERI_Dxz_Fxyz_Dyz_S_b = I_ERI_Fx2z_Dxy_Dyz_S_b+ABZ*I_ERI_Dxz_Dxy_Dyz_S_b;
  Double I_ERI_D2y_Fxyz_Dyz_S_b = I_ERI_F2yz_Dxy_Dyz_S_b+ABZ*I_ERI_D2y_Dxy_Dyz_S_b;
  Double I_ERI_Dyz_Fxyz_Dyz_S_b = I_ERI_Fy2z_Dxy_Dyz_S_b+ABZ*I_ERI_Dyz_Dxy_Dyz_S_b;
  Double I_ERI_D2z_Fxyz_Dyz_S_b = I_ERI_F3z_Dxy_Dyz_S_b+ABZ*I_ERI_D2z_Dxy_Dyz_S_b;
  Double I_ERI_D2x_Fx2z_Dyz_S_b = I_ERI_F3x_D2z_Dyz_S_b+ABX*I_ERI_D2x_D2z_Dyz_S_b;
  Double I_ERI_Dxy_Fx2z_Dyz_S_b = I_ERI_F2xy_D2z_Dyz_S_b+ABX*I_ERI_Dxy_D2z_Dyz_S_b;
  Double I_ERI_Dxz_Fx2z_Dyz_S_b = I_ERI_F2xz_D2z_Dyz_S_b+ABX*I_ERI_Dxz_D2z_Dyz_S_b;
  Double I_ERI_D2y_Fx2z_Dyz_S_b = I_ERI_Fx2y_D2z_Dyz_S_b+ABX*I_ERI_D2y_D2z_Dyz_S_b;
  Double I_ERI_Dyz_Fx2z_Dyz_S_b = I_ERI_Fxyz_D2z_Dyz_S_b+ABX*I_ERI_Dyz_D2z_Dyz_S_b;
  Double I_ERI_D2z_Fx2z_Dyz_S_b = I_ERI_Fx2z_D2z_Dyz_S_b+ABX*I_ERI_D2z_D2z_Dyz_S_b;
  Double I_ERI_D2x_F3y_Dyz_S_b = I_ERI_F2xy_D2y_Dyz_S_b+ABY*I_ERI_D2x_D2y_Dyz_S_b;
  Double I_ERI_Dxy_F3y_Dyz_S_b = I_ERI_Fx2y_D2y_Dyz_S_b+ABY*I_ERI_Dxy_D2y_Dyz_S_b;
  Double I_ERI_Dxz_F3y_Dyz_S_b = I_ERI_Fxyz_D2y_Dyz_S_b+ABY*I_ERI_Dxz_D2y_Dyz_S_b;
  Double I_ERI_D2y_F3y_Dyz_S_b = I_ERI_F3y_D2y_Dyz_S_b+ABY*I_ERI_D2y_D2y_Dyz_S_b;
  Double I_ERI_Dyz_F3y_Dyz_S_b = I_ERI_F2yz_D2y_Dyz_S_b+ABY*I_ERI_Dyz_D2y_Dyz_S_b;
  Double I_ERI_D2z_F3y_Dyz_S_b = I_ERI_Fy2z_D2y_Dyz_S_b+ABY*I_ERI_D2z_D2y_Dyz_S_b;
  Double I_ERI_D2x_F2yz_Dyz_S_b = I_ERI_F2xz_D2y_Dyz_S_b+ABZ*I_ERI_D2x_D2y_Dyz_S_b;
  Double I_ERI_Dxy_F2yz_Dyz_S_b = I_ERI_Fxyz_D2y_Dyz_S_b+ABZ*I_ERI_Dxy_D2y_Dyz_S_b;
  Double I_ERI_Dxz_F2yz_Dyz_S_b = I_ERI_Fx2z_D2y_Dyz_S_b+ABZ*I_ERI_Dxz_D2y_Dyz_S_b;
  Double I_ERI_D2y_F2yz_Dyz_S_b = I_ERI_F2yz_D2y_Dyz_S_b+ABZ*I_ERI_D2y_D2y_Dyz_S_b;
  Double I_ERI_Dyz_F2yz_Dyz_S_b = I_ERI_Fy2z_D2y_Dyz_S_b+ABZ*I_ERI_Dyz_D2y_Dyz_S_b;
  Double I_ERI_D2z_F2yz_Dyz_S_b = I_ERI_F3z_D2y_Dyz_S_b+ABZ*I_ERI_D2z_D2y_Dyz_S_b;
  Double I_ERI_D2x_Fy2z_Dyz_S_b = I_ERI_F2xy_D2z_Dyz_S_b+ABY*I_ERI_D2x_D2z_Dyz_S_b;
  Double I_ERI_Dxy_Fy2z_Dyz_S_b = I_ERI_Fx2y_D2z_Dyz_S_b+ABY*I_ERI_Dxy_D2z_Dyz_S_b;
  Double I_ERI_Dxz_Fy2z_Dyz_S_b = I_ERI_Fxyz_D2z_Dyz_S_b+ABY*I_ERI_Dxz_D2z_Dyz_S_b;
  Double I_ERI_D2y_Fy2z_Dyz_S_b = I_ERI_F3y_D2z_Dyz_S_b+ABY*I_ERI_D2y_D2z_Dyz_S_b;
  Double I_ERI_Dyz_Fy2z_Dyz_S_b = I_ERI_F2yz_D2z_Dyz_S_b+ABY*I_ERI_Dyz_D2z_Dyz_S_b;
  Double I_ERI_D2z_Fy2z_Dyz_S_b = I_ERI_Fy2z_D2z_Dyz_S_b+ABY*I_ERI_D2z_D2z_Dyz_S_b;
  Double I_ERI_D2x_F3z_Dyz_S_b = I_ERI_F2xz_D2z_Dyz_S_b+ABZ*I_ERI_D2x_D2z_Dyz_S_b;
  Double I_ERI_Dxy_F3z_Dyz_S_b = I_ERI_Fxyz_D2z_Dyz_S_b+ABZ*I_ERI_Dxy_D2z_Dyz_S_b;
  Double I_ERI_Dxz_F3z_Dyz_S_b = I_ERI_Fx2z_D2z_Dyz_S_b+ABZ*I_ERI_Dxz_D2z_Dyz_S_b;
  Double I_ERI_D2y_F3z_Dyz_S_b = I_ERI_F2yz_D2z_Dyz_S_b+ABZ*I_ERI_D2y_D2z_Dyz_S_b;
  Double I_ERI_Dyz_F3z_Dyz_S_b = I_ERI_Fy2z_D2z_Dyz_S_b+ABZ*I_ERI_Dyz_D2z_Dyz_S_b;
  Double I_ERI_D2z_F3z_Dyz_S_b = I_ERI_F3z_D2z_Dyz_S_b+ABZ*I_ERI_D2z_D2z_Dyz_S_b;
  Double I_ERI_D2x_F3x_D2z_S_b = I_ERI_F3x_D2x_D2z_S_b+ABX*I_ERI_D2x_D2x_D2z_S_b;
  Double I_ERI_Dxy_F3x_D2z_S_b = I_ERI_F2xy_D2x_D2z_S_b+ABX*I_ERI_Dxy_D2x_D2z_S_b;
  Double I_ERI_Dxz_F3x_D2z_S_b = I_ERI_F2xz_D2x_D2z_S_b+ABX*I_ERI_Dxz_D2x_D2z_S_b;
  Double I_ERI_D2y_F3x_D2z_S_b = I_ERI_Fx2y_D2x_D2z_S_b+ABX*I_ERI_D2y_D2x_D2z_S_b;
  Double I_ERI_Dyz_F3x_D2z_S_b = I_ERI_Fxyz_D2x_D2z_S_b+ABX*I_ERI_Dyz_D2x_D2z_S_b;
  Double I_ERI_D2z_F3x_D2z_S_b = I_ERI_Fx2z_D2x_D2z_S_b+ABX*I_ERI_D2z_D2x_D2z_S_b;
  Double I_ERI_D2x_F2xy_D2z_S_b = I_ERI_F2xy_D2x_D2z_S_b+ABY*I_ERI_D2x_D2x_D2z_S_b;
  Double I_ERI_Dxy_F2xy_D2z_S_b = I_ERI_Fx2y_D2x_D2z_S_b+ABY*I_ERI_Dxy_D2x_D2z_S_b;
  Double I_ERI_Dxz_F2xy_D2z_S_b = I_ERI_Fxyz_D2x_D2z_S_b+ABY*I_ERI_Dxz_D2x_D2z_S_b;
  Double I_ERI_D2y_F2xy_D2z_S_b = I_ERI_F3y_D2x_D2z_S_b+ABY*I_ERI_D2y_D2x_D2z_S_b;
  Double I_ERI_Dyz_F2xy_D2z_S_b = I_ERI_F2yz_D2x_D2z_S_b+ABY*I_ERI_Dyz_D2x_D2z_S_b;
  Double I_ERI_D2z_F2xy_D2z_S_b = I_ERI_Fy2z_D2x_D2z_S_b+ABY*I_ERI_D2z_D2x_D2z_S_b;
  Double I_ERI_D2x_F2xz_D2z_S_b = I_ERI_F2xz_D2x_D2z_S_b+ABZ*I_ERI_D2x_D2x_D2z_S_b;
  Double I_ERI_Dxy_F2xz_D2z_S_b = I_ERI_Fxyz_D2x_D2z_S_b+ABZ*I_ERI_Dxy_D2x_D2z_S_b;
  Double I_ERI_Dxz_F2xz_D2z_S_b = I_ERI_Fx2z_D2x_D2z_S_b+ABZ*I_ERI_Dxz_D2x_D2z_S_b;
  Double I_ERI_D2y_F2xz_D2z_S_b = I_ERI_F2yz_D2x_D2z_S_b+ABZ*I_ERI_D2y_D2x_D2z_S_b;
  Double I_ERI_Dyz_F2xz_D2z_S_b = I_ERI_Fy2z_D2x_D2z_S_b+ABZ*I_ERI_Dyz_D2x_D2z_S_b;
  Double I_ERI_D2z_F2xz_D2z_S_b = I_ERI_F3z_D2x_D2z_S_b+ABZ*I_ERI_D2z_D2x_D2z_S_b;
  Double I_ERI_D2x_Fx2y_D2z_S_b = I_ERI_F3x_D2y_D2z_S_b+ABX*I_ERI_D2x_D2y_D2z_S_b;
  Double I_ERI_Dxy_Fx2y_D2z_S_b = I_ERI_F2xy_D2y_D2z_S_b+ABX*I_ERI_Dxy_D2y_D2z_S_b;
  Double I_ERI_Dxz_Fx2y_D2z_S_b = I_ERI_F2xz_D2y_D2z_S_b+ABX*I_ERI_Dxz_D2y_D2z_S_b;
  Double I_ERI_D2y_Fx2y_D2z_S_b = I_ERI_Fx2y_D2y_D2z_S_b+ABX*I_ERI_D2y_D2y_D2z_S_b;
  Double I_ERI_Dyz_Fx2y_D2z_S_b = I_ERI_Fxyz_D2y_D2z_S_b+ABX*I_ERI_Dyz_D2y_D2z_S_b;
  Double I_ERI_D2z_Fx2y_D2z_S_b = I_ERI_Fx2z_D2y_D2z_S_b+ABX*I_ERI_D2z_D2y_D2z_S_b;
  Double I_ERI_D2x_Fxyz_D2z_S_b = I_ERI_F2xz_Dxy_D2z_S_b+ABZ*I_ERI_D2x_Dxy_D2z_S_b;
  Double I_ERI_Dxy_Fxyz_D2z_S_b = I_ERI_Fxyz_Dxy_D2z_S_b+ABZ*I_ERI_Dxy_Dxy_D2z_S_b;
  Double I_ERI_Dxz_Fxyz_D2z_S_b = I_ERI_Fx2z_Dxy_D2z_S_b+ABZ*I_ERI_Dxz_Dxy_D2z_S_b;
  Double I_ERI_D2y_Fxyz_D2z_S_b = I_ERI_F2yz_Dxy_D2z_S_b+ABZ*I_ERI_D2y_Dxy_D2z_S_b;
  Double I_ERI_Dyz_Fxyz_D2z_S_b = I_ERI_Fy2z_Dxy_D2z_S_b+ABZ*I_ERI_Dyz_Dxy_D2z_S_b;
  Double I_ERI_D2z_Fxyz_D2z_S_b = I_ERI_F3z_Dxy_D2z_S_b+ABZ*I_ERI_D2z_Dxy_D2z_S_b;
  Double I_ERI_D2x_Fx2z_D2z_S_b = I_ERI_F3x_D2z_D2z_S_b+ABX*I_ERI_D2x_D2z_D2z_S_b;
  Double I_ERI_Dxy_Fx2z_D2z_S_b = I_ERI_F2xy_D2z_D2z_S_b+ABX*I_ERI_Dxy_D2z_D2z_S_b;
  Double I_ERI_Dxz_Fx2z_D2z_S_b = I_ERI_F2xz_D2z_D2z_S_b+ABX*I_ERI_Dxz_D2z_D2z_S_b;
  Double I_ERI_D2y_Fx2z_D2z_S_b = I_ERI_Fx2y_D2z_D2z_S_b+ABX*I_ERI_D2y_D2z_D2z_S_b;
  Double I_ERI_Dyz_Fx2z_D2z_S_b = I_ERI_Fxyz_D2z_D2z_S_b+ABX*I_ERI_Dyz_D2z_D2z_S_b;
  Double I_ERI_D2z_Fx2z_D2z_S_b = I_ERI_Fx2z_D2z_D2z_S_b+ABX*I_ERI_D2z_D2z_D2z_S_b;
  Double I_ERI_D2x_F3y_D2z_S_b = I_ERI_F2xy_D2y_D2z_S_b+ABY*I_ERI_D2x_D2y_D2z_S_b;
  Double I_ERI_Dxy_F3y_D2z_S_b = I_ERI_Fx2y_D2y_D2z_S_b+ABY*I_ERI_Dxy_D2y_D2z_S_b;
  Double I_ERI_Dxz_F3y_D2z_S_b = I_ERI_Fxyz_D2y_D2z_S_b+ABY*I_ERI_Dxz_D2y_D2z_S_b;
  Double I_ERI_D2y_F3y_D2z_S_b = I_ERI_F3y_D2y_D2z_S_b+ABY*I_ERI_D2y_D2y_D2z_S_b;
  Double I_ERI_Dyz_F3y_D2z_S_b = I_ERI_F2yz_D2y_D2z_S_b+ABY*I_ERI_Dyz_D2y_D2z_S_b;
  Double I_ERI_D2z_F3y_D2z_S_b = I_ERI_Fy2z_D2y_D2z_S_b+ABY*I_ERI_D2z_D2y_D2z_S_b;
  Double I_ERI_D2x_F2yz_D2z_S_b = I_ERI_F2xz_D2y_D2z_S_b+ABZ*I_ERI_D2x_D2y_D2z_S_b;
  Double I_ERI_Dxy_F2yz_D2z_S_b = I_ERI_Fxyz_D2y_D2z_S_b+ABZ*I_ERI_Dxy_D2y_D2z_S_b;
  Double I_ERI_Dxz_F2yz_D2z_S_b = I_ERI_Fx2z_D2y_D2z_S_b+ABZ*I_ERI_Dxz_D2y_D2z_S_b;
  Double I_ERI_D2y_F2yz_D2z_S_b = I_ERI_F2yz_D2y_D2z_S_b+ABZ*I_ERI_D2y_D2y_D2z_S_b;
  Double I_ERI_Dyz_F2yz_D2z_S_b = I_ERI_Fy2z_D2y_D2z_S_b+ABZ*I_ERI_Dyz_D2y_D2z_S_b;
  Double I_ERI_D2z_F2yz_D2z_S_b = I_ERI_F3z_D2y_D2z_S_b+ABZ*I_ERI_D2z_D2y_D2z_S_b;
  Double I_ERI_D2x_Fy2z_D2z_S_b = I_ERI_F2xy_D2z_D2z_S_b+ABY*I_ERI_D2x_D2z_D2z_S_b;
  Double I_ERI_Dxy_Fy2z_D2z_S_b = I_ERI_Fx2y_D2z_D2z_S_b+ABY*I_ERI_Dxy_D2z_D2z_S_b;
  Double I_ERI_Dxz_Fy2z_D2z_S_b = I_ERI_Fxyz_D2z_D2z_S_b+ABY*I_ERI_Dxz_D2z_D2z_S_b;
  Double I_ERI_D2y_Fy2z_D2z_S_b = I_ERI_F3y_D2z_D2z_S_b+ABY*I_ERI_D2y_D2z_D2z_S_b;
  Double I_ERI_Dyz_Fy2z_D2z_S_b = I_ERI_F2yz_D2z_D2z_S_b+ABY*I_ERI_Dyz_D2z_D2z_S_b;
  Double I_ERI_D2z_Fy2z_D2z_S_b = I_ERI_Fy2z_D2z_D2z_S_b+ABY*I_ERI_D2z_D2z_D2z_S_b;
  Double I_ERI_D2x_F3z_D2z_S_b = I_ERI_F2xz_D2z_D2z_S_b+ABZ*I_ERI_D2x_D2z_D2z_S_b;
  Double I_ERI_Dxy_F3z_D2z_S_b = I_ERI_Fxyz_D2z_D2z_S_b+ABZ*I_ERI_Dxy_D2z_D2z_S_b;
  Double I_ERI_Dxz_F3z_D2z_S_b = I_ERI_Fx2z_D2z_D2z_S_b+ABZ*I_ERI_Dxz_D2z_D2z_S_b;
  Double I_ERI_D2y_F3z_D2z_S_b = I_ERI_F2yz_D2z_D2z_S_b+ABZ*I_ERI_D2y_D2z_D2z_S_b;
  Double I_ERI_Dyz_F3z_D2z_S_b = I_ERI_Fy2z_D2z_D2z_S_b+ABZ*I_ERI_Dyz_D2z_D2z_S_b;
  Double I_ERI_D2z_F3z_D2z_S_b = I_ERI_F3z_D2z_D2z_S_b+ABZ*I_ERI_D2z_D2z_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_c
   * RHS shell quartet name: SQ_ERI_D_S_F_S_c
   ************************************************************/
  Double I_ERI_D2x_Px_F3x_S_c = I_ERI_F3x_S_F3x_S_c+ABX*I_ERI_D2x_S_F3x_S_c;
  Double I_ERI_Dxy_Px_F3x_S_c = I_ERI_F2xy_S_F3x_S_c+ABX*I_ERI_Dxy_S_F3x_S_c;
  Double I_ERI_Dxz_Px_F3x_S_c = I_ERI_F2xz_S_F3x_S_c+ABX*I_ERI_Dxz_S_F3x_S_c;
  Double I_ERI_D2y_Px_F3x_S_c = I_ERI_Fx2y_S_F3x_S_c+ABX*I_ERI_D2y_S_F3x_S_c;
  Double I_ERI_Dyz_Px_F3x_S_c = I_ERI_Fxyz_S_F3x_S_c+ABX*I_ERI_Dyz_S_F3x_S_c;
  Double I_ERI_D2z_Px_F3x_S_c = I_ERI_Fx2z_S_F3x_S_c+ABX*I_ERI_D2z_S_F3x_S_c;
  Double I_ERI_D2x_Py_F3x_S_c = I_ERI_F2xy_S_F3x_S_c+ABY*I_ERI_D2x_S_F3x_S_c;
  Double I_ERI_Dxy_Py_F3x_S_c = I_ERI_Fx2y_S_F3x_S_c+ABY*I_ERI_Dxy_S_F3x_S_c;
  Double I_ERI_Dxz_Py_F3x_S_c = I_ERI_Fxyz_S_F3x_S_c+ABY*I_ERI_Dxz_S_F3x_S_c;
  Double I_ERI_D2y_Py_F3x_S_c = I_ERI_F3y_S_F3x_S_c+ABY*I_ERI_D2y_S_F3x_S_c;
  Double I_ERI_Dyz_Py_F3x_S_c = I_ERI_F2yz_S_F3x_S_c+ABY*I_ERI_Dyz_S_F3x_S_c;
  Double I_ERI_D2z_Py_F3x_S_c = I_ERI_Fy2z_S_F3x_S_c+ABY*I_ERI_D2z_S_F3x_S_c;
  Double I_ERI_D2x_Pz_F3x_S_c = I_ERI_F2xz_S_F3x_S_c+ABZ*I_ERI_D2x_S_F3x_S_c;
  Double I_ERI_Dxy_Pz_F3x_S_c = I_ERI_Fxyz_S_F3x_S_c+ABZ*I_ERI_Dxy_S_F3x_S_c;
  Double I_ERI_Dxz_Pz_F3x_S_c = I_ERI_Fx2z_S_F3x_S_c+ABZ*I_ERI_Dxz_S_F3x_S_c;
  Double I_ERI_D2y_Pz_F3x_S_c = I_ERI_F2yz_S_F3x_S_c+ABZ*I_ERI_D2y_S_F3x_S_c;
  Double I_ERI_Dyz_Pz_F3x_S_c = I_ERI_Fy2z_S_F3x_S_c+ABZ*I_ERI_Dyz_S_F3x_S_c;
  Double I_ERI_D2z_Pz_F3x_S_c = I_ERI_F3z_S_F3x_S_c+ABZ*I_ERI_D2z_S_F3x_S_c;
  Double I_ERI_D2x_Px_F2xy_S_c = I_ERI_F3x_S_F2xy_S_c+ABX*I_ERI_D2x_S_F2xy_S_c;
  Double I_ERI_Dxy_Px_F2xy_S_c = I_ERI_F2xy_S_F2xy_S_c+ABX*I_ERI_Dxy_S_F2xy_S_c;
  Double I_ERI_Dxz_Px_F2xy_S_c = I_ERI_F2xz_S_F2xy_S_c+ABX*I_ERI_Dxz_S_F2xy_S_c;
  Double I_ERI_D2y_Px_F2xy_S_c = I_ERI_Fx2y_S_F2xy_S_c+ABX*I_ERI_D2y_S_F2xy_S_c;
  Double I_ERI_Dyz_Px_F2xy_S_c = I_ERI_Fxyz_S_F2xy_S_c+ABX*I_ERI_Dyz_S_F2xy_S_c;
  Double I_ERI_D2z_Px_F2xy_S_c = I_ERI_Fx2z_S_F2xy_S_c+ABX*I_ERI_D2z_S_F2xy_S_c;
  Double I_ERI_D2x_Py_F2xy_S_c = I_ERI_F2xy_S_F2xy_S_c+ABY*I_ERI_D2x_S_F2xy_S_c;
  Double I_ERI_Dxy_Py_F2xy_S_c = I_ERI_Fx2y_S_F2xy_S_c+ABY*I_ERI_Dxy_S_F2xy_S_c;
  Double I_ERI_Dxz_Py_F2xy_S_c = I_ERI_Fxyz_S_F2xy_S_c+ABY*I_ERI_Dxz_S_F2xy_S_c;
  Double I_ERI_D2y_Py_F2xy_S_c = I_ERI_F3y_S_F2xy_S_c+ABY*I_ERI_D2y_S_F2xy_S_c;
  Double I_ERI_Dyz_Py_F2xy_S_c = I_ERI_F2yz_S_F2xy_S_c+ABY*I_ERI_Dyz_S_F2xy_S_c;
  Double I_ERI_D2z_Py_F2xy_S_c = I_ERI_Fy2z_S_F2xy_S_c+ABY*I_ERI_D2z_S_F2xy_S_c;
  Double I_ERI_D2x_Pz_F2xy_S_c = I_ERI_F2xz_S_F2xy_S_c+ABZ*I_ERI_D2x_S_F2xy_S_c;
  Double I_ERI_Dxy_Pz_F2xy_S_c = I_ERI_Fxyz_S_F2xy_S_c+ABZ*I_ERI_Dxy_S_F2xy_S_c;
  Double I_ERI_Dxz_Pz_F2xy_S_c = I_ERI_Fx2z_S_F2xy_S_c+ABZ*I_ERI_Dxz_S_F2xy_S_c;
  Double I_ERI_D2y_Pz_F2xy_S_c = I_ERI_F2yz_S_F2xy_S_c+ABZ*I_ERI_D2y_S_F2xy_S_c;
  Double I_ERI_Dyz_Pz_F2xy_S_c = I_ERI_Fy2z_S_F2xy_S_c+ABZ*I_ERI_Dyz_S_F2xy_S_c;
  Double I_ERI_D2z_Pz_F2xy_S_c = I_ERI_F3z_S_F2xy_S_c+ABZ*I_ERI_D2z_S_F2xy_S_c;
  Double I_ERI_D2x_Px_F2xz_S_c = I_ERI_F3x_S_F2xz_S_c+ABX*I_ERI_D2x_S_F2xz_S_c;
  Double I_ERI_Dxy_Px_F2xz_S_c = I_ERI_F2xy_S_F2xz_S_c+ABX*I_ERI_Dxy_S_F2xz_S_c;
  Double I_ERI_Dxz_Px_F2xz_S_c = I_ERI_F2xz_S_F2xz_S_c+ABX*I_ERI_Dxz_S_F2xz_S_c;
  Double I_ERI_D2y_Px_F2xz_S_c = I_ERI_Fx2y_S_F2xz_S_c+ABX*I_ERI_D2y_S_F2xz_S_c;
  Double I_ERI_Dyz_Px_F2xz_S_c = I_ERI_Fxyz_S_F2xz_S_c+ABX*I_ERI_Dyz_S_F2xz_S_c;
  Double I_ERI_D2z_Px_F2xz_S_c = I_ERI_Fx2z_S_F2xz_S_c+ABX*I_ERI_D2z_S_F2xz_S_c;
  Double I_ERI_D2x_Py_F2xz_S_c = I_ERI_F2xy_S_F2xz_S_c+ABY*I_ERI_D2x_S_F2xz_S_c;
  Double I_ERI_Dxy_Py_F2xz_S_c = I_ERI_Fx2y_S_F2xz_S_c+ABY*I_ERI_Dxy_S_F2xz_S_c;
  Double I_ERI_Dxz_Py_F2xz_S_c = I_ERI_Fxyz_S_F2xz_S_c+ABY*I_ERI_Dxz_S_F2xz_S_c;
  Double I_ERI_D2y_Py_F2xz_S_c = I_ERI_F3y_S_F2xz_S_c+ABY*I_ERI_D2y_S_F2xz_S_c;
  Double I_ERI_Dyz_Py_F2xz_S_c = I_ERI_F2yz_S_F2xz_S_c+ABY*I_ERI_Dyz_S_F2xz_S_c;
  Double I_ERI_D2z_Py_F2xz_S_c = I_ERI_Fy2z_S_F2xz_S_c+ABY*I_ERI_D2z_S_F2xz_S_c;
  Double I_ERI_D2x_Pz_F2xz_S_c = I_ERI_F2xz_S_F2xz_S_c+ABZ*I_ERI_D2x_S_F2xz_S_c;
  Double I_ERI_Dxy_Pz_F2xz_S_c = I_ERI_Fxyz_S_F2xz_S_c+ABZ*I_ERI_Dxy_S_F2xz_S_c;
  Double I_ERI_Dxz_Pz_F2xz_S_c = I_ERI_Fx2z_S_F2xz_S_c+ABZ*I_ERI_Dxz_S_F2xz_S_c;
  Double I_ERI_D2y_Pz_F2xz_S_c = I_ERI_F2yz_S_F2xz_S_c+ABZ*I_ERI_D2y_S_F2xz_S_c;
  Double I_ERI_Dyz_Pz_F2xz_S_c = I_ERI_Fy2z_S_F2xz_S_c+ABZ*I_ERI_Dyz_S_F2xz_S_c;
  Double I_ERI_D2z_Pz_F2xz_S_c = I_ERI_F3z_S_F2xz_S_c+ABZ*I_ERI_D2z_S_F2xz_S_c;
  Double I_ERI_D2x_Px_Fx2y_S_c = I_ERI_F3x_S_Fx2y_S_c+ABX*I_ERI_D2x_S_Fx2y_S_c;
  Double I_ERI_Dxy_Px_Fx2y_S_c = I_ERI_F2xy_S_Fx2y_S_c+ABX*I_ERI_Dxy_S_Fx2y_S_c;
  Double I_ERI_Dxz_Px_Fx2y_S_c = I_ERI_F2xz_S_Fx2y_S_c+ABX*I_ERI_Dxz_S_Fx2y_S_c;
  Double I_ERI_D2y_Px_Fx2y_S_c = I_ERI_Fx2y_S_Fx2y_S_c+ABX*I_ERI_D2y_S_Fx2y_S_c;
  Double I_ERI_Dyz_Px_Fx2y_S_c = I_ERI_Fxyz_S_Fx2y_S_c+ABX*I_ERI_Dyz_S_Fx2y_S_c;
  Double I_ERI_D2z_Px_Fx2y_S_c = I_ERI_Fx2z_S_Fx2y_S_c+ABX*I_ERI_D2z_S_Fx2y_S_c;
  Double I_ERI_D2x_Py_Fx2y_S_c = I_ERI_F2xy_S_Fx2y_S_c+ABY*I_ERI_D2x_S_Fx2y_S_c;
  Double I_ERI_Dxy_Py_Fx2y_S_c = I_ERI_Fx2y_S_Fx2y_S_c+ABY*I_ERI_Dxy_S_Fx2y_S_c;
  Double I_ERI_Dxz_Py_Fx2y_S_c = I_ERI_Fxyz_S_Fx2y_S_c+ABY*I_ERI_Dxz_S_Fx2y_S_c;
  Double I_ERI_D2y_Py_Fx2y_S_c = I_ERI_F3y_S_Fx2y_S_c+ABY*I_ERI_D2y_S_Fx2y_S_c;
  Double I_ERI_Dyz_Py_Fx2y_S_c = I_ERI_F2yz_S_Fx2y_S_c+ABY*I_ERI_Dyz_S_Fx2y_S_c;
  Double I_ERI_D2z_Py_Fx2y_S_c = I_ERI_Fy2z_S_Fx2y_S_c+ABY*I_ERI_D2z_S_Fx2y_S_c;
  Double I_ERI_D2x_Pz_Fx2y_S_c = I_ERI_F2xz_S_Fx2y_S_c+ABZ*I_ERI_D2x_S_Fx2y_S_c;
  Double I_ERI_Dxy_Pz_Fx2y_S_c = I_ERI_Fxyz_S_Fx2y_S_c+ABZ*I_ERI_Dxy_S_Fx2y_S_c;
  Double I_ERI_Dxz_Pz_Fx2y_S_c = I_ERI_Fx2z_S_Fx2y_S_c+ABZ*I_ERI_Dxz_S_Fx2y_S_c;
  Double I_ERI_D2y_Pz_Fx2y_S_c = I_ERI_F2yz_S_Fx2y_S_c+ABZ*I_ERI_D2y_S_Fx2y_S_c;
  Double I_ERI_Dyz_Pz_Fx2y_S_c = I_ERI_Fy2z_S_Fx2y_S_c+ABZ*I_ERI_Dyz_S_Fx2y_S_c;
  Double I_ERI_D2z_Pz_Fx2y_S_c = I_ERI_F3z_S_Fx2y_S_c+ABZ*I_ERI_D2z_S_Fx2y_S_c;
  Double I_ERI_D2x_Px_Fxyz_S_c = I_ERI_F3x_S_Fxyz_S_c+ABX*I_ERI_D2x_S_Fxyz_S_c;
  Double I_ERI_Dxy_Px_Fxyz_S_c = I_ERI_F2xy_S_Fxyz_S_c+ABX*I_ERI_Dxy_S_Fxyz_S_c;
  Double I_ERI_Dxz_Px_Fxyz_S_c = I_ERI_F2xz_S_Fxyz_S_c+ABX*I_ERI_Dxz_S_Fxyz_S_c;
  Double I_ERI_D2y_Px_Fxyz_S_c = I_ERI_Fx2y_S_Fxyz_S_c+ABX*I_ERI_D2y_S_Fxyz_S_c;
  Double I_ERI_Dyz_Px_Fxyz_S_c = I_ERI_Fxyz_S_Fxyz_S_c+ABX*I_ERI_Dyz_S_Fxyz_S_c;
  Double I_ERI_D2z_Px_Fxyz_S_c = I_ERI_Fx2z_S_Fxyz_S_c+ABX*I_ERI_D2z_S_Fxyz_S_c;
  Double I_ERI_D2x_Py_Fxyz_S_c = I_ERI_F2xy_S_Fxyz_S_c+ABY*I_ERI_D2x_S_Fxyz_S_c;
  Double I_ERI_Dxy_Py_Fxyz_S_c = I_ERI_Fx2y_S_Fxyz_S_c+ABY*I_ERI_Dxy_S_Fxyz_S_c;
  Double I_ERI_Dxz_Py_Fxyz_S_c = I_ERI_Fxyz_S_Fxyz_S_c+ABY*I_ERI_Dxz_S_Fxyz_S_c;
  Double I_ERI_D2y_Py_Fxyz_S_c = I_ERI_F3y_S_Fxyz_S_c+ABY*I_ERI_D2y_S_Fxyz_S_c;
  Double I_ERI_Dyz_Py_Fxyz_S_c = I_ERI_F2yz_S_Fxyz_S_c+ABY*I_ERI_Dyz_S_Fxyz_S_c;
  Double I_ERI_D2z_Py_Fxyz_S_c = I_ERI_Fy2z_S_Fxyz_S_c+ABY*I_ERI_D2z_S_Fxyz_S_c;
  Double I_ERI_D2x_Pz_Fxyz_S_c = I_ERI_F2xz_S_Fxyz_S_c+ABZ*I_ERI_D2x_S_Fxyz_S_c;
  Double I_ERI_Dxy_Pz_Fxyz_S_c = I_ERI_Fxyz_S_Fxyz_S_c+ABZ*I_ERI_Dxy_S_Fxyz_S_c;
  Double I_ERI_Dxz_Pz_Fxyz_S_c = I_ERI_Fx2z_S_Fxyz_S_c+ABZ*I_ERI_Dxz_S_Fxyz_S_c;
  Double I_ERI_D2y_Pz_Fxyz_S_c = I_ERI_F2yz_S_Fxyz_S_c+ABZ*I_ERI_D2y_S_Fxyz_S_c;
  Double I_ERI_Dyz_Pz_Fxyz_S_c = I_ERI_Fy2z_S_Fxyz_S_c+ABZ*I_ERI_Dyz_S_Fxyz_S_c;
  Double I_ERI_D2z_Pz_Fxyz_S_c = I_ERI_F3z_S_Fxyz_S_c+ABZ*I_ERI_D2z_S_Fxyz_S_c;
  Double I_ERI_D2x_Px_Fx2z_S_c = I_ERI_F3x_S_Fx2z_S_c+ABX*I_ERI_D2x_S_Fx2z_S_c;
  Double I_ERI_Dxy_Px_Fx2z_S_c = I_ERI_F2xy_S_Fx2z_S_c+ABX*I_ERI_Dxy_S_Fx2z_S_c;
  Double I_ERI_Dxz_Px_Fx2z_S_c = I_ERI_F2xz_S_Fx2z_S_c+ABX*I_ERI_Dxz_S_Fx2z_S_c;
  Double I_ERI_D2y_Px_Fx2z_S_c = I_ERI_Fx2y_S_Fx2z_S_c+ABX*I_ERI_D2y_S_Fx2z_S_c;
  Double I_ERI_Dyz_Px_Fx2z_S_c = I_ERI_Fxyz_S_Fx2z_S_c+ABX*I_ERI_Dyz_S_Fx2z_S_c;
  Double I_ERI_D2z_Px_Fx2z_S_c = I_ERI_Fx2z_S_Fx2z_S_c+ABX*I_ERI_D2z_S_Fx2z_S_c;
  Double I_ERI_D2x_Py_Fx2z_S_c = I_ERI_F2xy_S_Fx2z_S_c+ABY*I_ERI_D2x_S_Fx2z_S_c;
  Double I_ERI_Dxy_Py_Fx2z_S_c = I_ERI_Fx2y_S_Fx2z_S_c+ABY*I_ERI_Dxy_S_Fx2z_S_c;
  Double I_ERI_Dxz_Py_Fx2z_S_c = I_ERI_Fxyz_S_Fx2z_S_c+ABY*I_ERI_Dxz_S_Fx2z_S_c;
  Double I_ERI_D2y_Py_Fx2z_S_c = I_ERI_F3y_S_Fx2z_S_c+ABY*I_ERI_D2y_S_Fx2z_S_c;
  Double I_ERI_Dyz_Py_Fx2z_S_c = I_ERI_F2yz_S_Fx2z_S_c+ABY*I_ERI_Dyz_S_Fx2z_S_c;
  Double I_ERI_D2z_Py_Fx2z_S_c = I_ERI_Fy2z_S_Fx2z_S_c+ABY*I_ERI_D2z_S_Fx2z_S_c;
  Double I_ERI_D2x_Pz_Fx2z_S_c = I_ERI_F2xz_S_Fx2z_S_c+ABZ*I_ERI_D2x_S_Fx2z_S_c;
  Double I_ERI_Dxy_Pz_Fx2z_S_c = I_ERI_Fxyz_S_Fx2z_S_c+ABZ*I_ERI_Dxy_S_Fx2z_S_c;
  Double I_ERI_Dxz_Pz_Fx2z_S_c = I_ERI_Fx2z_S_Fx2z_S_c+ABZ*I_ERI_Dxz_S_Fx2z_S_c;
  Double I_ERI_D2y_Pz_Fx2z_S_c = I_ERI_F2yz_S_Fx2z_S_c+ABZ*I_ERI_D2y_S_Fx2z_S_c;
  Double I_ERI_Dyz_Pz_Fx2z_S_c = I_ERI_Fy2z_S_Fx2z_S_c+ABZ*I_ERI_Dyz_S_Fx2z_S_c;
  Double I_ERI_D2z_Pz_Fx2z_S_c = I_ERI_F3z_S_Fx2z_S_c+ABZ*I_ERI_D2z_S_Fx2z_S_c;
  Double I_ERI_D2x_Px_F3y_S_c = I_ERI_F3x_S_F3y_S_c+ABX*I_ERI_D2x_S_F3y_S_c;
  Double I_ERI_Dxy_Px_F3y_S_c = I_ERI_F2xy_S_F3y_S_c+ABX*I_ERI_Dxy_S_F3y_S_c;
  Double I_ERI_Dxz_Px_F3y_S_c = I_ERI_F2xz_S_F3y_S_c+ABX*I_ERI_Dxz_S_F3y_S_c;
  Double I_ERI_D2y_Px_F3y_S_c = I_ERI_Fx2y_S_F3y_S_c+ABX*I_ERI_D2y_S_F3y_S_c;
  Double I_ERI_Dyz_Px_F3y_S_c = I_ERI_Fxyz_S_F3y_S_c+ABX*I_ERI_Dyz_S_F3y_S_c;
  Double I_ERI_D2z_Px_F3y_S_c = I_ERI_Fx2z_S_F3y_S_c+ABX*I_ERI_D2z_S_F3y_S_c;
  Double I_ERI_D2x_Py_F3y_S_c = I_ERI_F2xy_S_F3y_S_c+ABY*I_ERI_D2x_S_F3y_S_c;
  Double I_ERI_Dxy_Py_F3y_S_c = I_ERI_Fx2y_S_F3y_S_c+ABY*I_ERI_Dxy_S_F3y_S_c;
  Double I_ERI_Dxz_Py_F3y_S_c = I_ERI_Fxyz_S_F3y_S_c+ABY*I_ERI_Dxz_S_F3y_S_c;
  Double I_ERI_D2y_Py_F3y_S_c = I_ERI_F3y_S_F3y_S_c+ABY*I_ERI_D2y_S_F3y_S_c;
  Double I_ERI_Dyz_Py_F3y_S_c = I_ERI_F2yz_S_F3y_S_c+ABY*I_ERI_Dyz_S_F3y_S_c;
  Double I_ERI_D2z_Py_F3y_S_c = I_ERI_Fy2z_S_F3y_S_c+ABY*I_ERI_D2z_S_F3y_S_c;
  Double I_ERI_D2x_Pz_F3y_S_c = I_ERI_F2xz_S_F3y_S_c+ABZ*I_ERI_D2x_S_F3y_S_c;
  Double I_ERI_Dxy_Pz_F3y_S_c = I_ERI_Fxyz_S_F3y_S_c+ABZ*I_ERI_Dxy_S_F3y_S_c;
  Double I_ERI_Dxz_Pz_F3y_S_c = I_ERI_Fx2z_S_F3y_S_c+ABZ*I_ERI_Dxz_S_F3y_S_c;
  Double I_ERI_D2y_Pz_F3y_S_c = I_ERI_F2yz_S_F3y_S_c+ABZ*I_ERI_D2y_S_F3y_S_c;
  Double I_ERI_Dyz_Pz_F3y_S_c = I_ERI_Fy2z_S_F3y_S_c+ABZ*I_ERI_Dyz_S_F3y_S_c;
  Double I_ERI_D2z_Pz_F3y_S_c = I_ERI_F3z_S_F3y_S_c+ABZ*I_ERI_D2z_S_F3y_S_c;
  Double I_ERI_D2x_Px_F2yz_S_c = I_ERI_F3x_S_F2yz_S_c+ABX*I_ERI_D2x_S_F2yz_S_c;
  Double I_ERI_Dxy_Px_F2yz_S_c = I_ERI_F2xy_S_F2yz_S_c+ABX*I_ERI_Dxy_S_F2yz_S_c;
  Double I_ERI_Dxz_Px_F2yz_S_c = I_ERI_F2xz_S_F2yz_S_c+ABX*I_ERI_Dxz_S_F2yz_S_c;
  Double I_ERI_D2y_Px_F2yz_S_c = I_ERI_Fx2y_S_F2yz_S_c+ABX*I_ERI_D2y_S_F2yz_S_c;
  Double I_ERI_Dyz_Px_F2yz_S_c = I_ERI_Fxyz_S_F2yz_S_c+ABX*I_ERI_Dyz_S_F2yz_S_c;
  Double I_ERI_D2z_Px_F2yz_S_c = I_ERI_Fx2z_S_F2yz_S_c+ABX*I_ERI_D2z_S_F2yz_S_c;
  Double I_ERI_D2x_Py_F2yz_S_c = I_ERI_F2xy_S_F2yz_S_c+ABY*I_ERI_D2x_S_F2yz_S_c;
  Double I_ERI_Dxy_Py_F2yz_S_c = I_ERI_Fx2y_S_F2yz_S_c+ABY*I_ERI_Dxy_S_F2yz_S_c;
  Double I_ERI_Dxz_Py_F2yz_S_c = I_ERI_Fxyz_S_F2yz_S_c+ABY*I_ERI_Dxz_S_F2yz_S_c;
  Double I_ERI_D2y_Py_F2yz_S_c = I_ERI_F3y_S_F2yz_S_c+ABY*I_ERI_D2y_S_F2yz_S_c;
  Double I_ERI_Dyz_Py_F2yz_S_c = I_ERI_F2yz_S_F2yz_S_c+ABY*I_ERI_Dyz_S_F2yz_S_c;
  Double I_ERI_D2z_Py_F2yz_S_c = I_ERI_Fy2z_S_F2yz_S_c+ABY*I_ERI_D2z_S_F2yz_S_c;
  Double I_ERI_D2x_Pz_F2yz_S_c = I_ERI_F2xz_S_F2yz_S_c+ABZ*I_ERI_D2x_S_F2yz_S_c;
  Double I_ERI_Dxy_Pz_F2yz_S_c = I_ERI_Fxyz_S_F2yz_S_c+ABZ*I_ERI_Dxy_S_F2yz_S_c;
  Double I_ERI_Dxz_Pz_F2yz_S_c = I_ERI_Fx2z_S_F2yz_S_c+ABZ*I_ERI_Dxz_S_F2yz_S_c;
  Double I_ERI_D2y_Pz_F2yz_S_c = I_ERI_F2yz_S_F2yz_S_c+ABZ*I_ERI_D2y_S_F2yz_S_c;
  Double I_ERI_Dyz_Pz_F2yz_S_c = I_ERI_Fy2z_S_F2yz_S_c+ABZ*I_ERI_Dyz_S_F2yz_S_c;
  Double I_ERI_D2z_Pz_F2yz_S_c = I_ERI_F3z_S_F2yz_S_c+ABZ*I_ERI_D2z_S_F2yz_S_c;
  Double I_ERI_D2x_Px_Fy2z_S_c = I_ERI_F3x_S_Fy2z_S_c+ABX*I_ERI_D2x_S_Fy2z_S_c;
  Double I_ERI_Dxy_Px_Fy2z_S_c = I_ERI_F2xy_S_Fy2z_S_c+ABX*I_ERI_Dxy_S_Fy2z_S_c;
  Double I_ERI_Dxz_Px_Fy2z_S_c = I_ERI_F2xz_S_Fy2z_S_c+ABX*I_ERI_Dxz_S_Fy2z_S_c;
  Double I_ERI_D2y_Px_Fy2z_S_c = I_ERI_Fx2y_S_Fy2z_S_c+ABX*I_ERI_D2y_S_Fy2z_S_c;
  Double I_ERI_Dyz_Px_Fy2z_S_c = I_ERI_Fxyz_S_Fy2z_S_c+ABX*I_ERI_Dyz_S_Fy2z_S_c;
  Double I_ERI_D2z_Px_Fy2z_S_c = I_ERI_Fx2z_S_Fy2z_S_c+ABX*I_ERI_D2z_S_Fy2z_S_c;
  Double I_ERI_D2x_Py_Fy2z_S_c = I_ERI_F2xy_S_Fy2z_S_c+ABY*I_ERI_D2x_S_Fy2z_S_c;
  Double I_ERI_Dxy_Py_Fy2z_S_c = I_ERI_Fx2y_S_Fy2z_S_c+ABY*I_ERI_Dxy_S_Fy2z_S_c;
  Double I_ERI_Dxz_Py_Fy2z_S_c = I_ERI_Fxyz_S_Fy2z_S_c+ABY*I_ERI_Dxz_S_Fy2z_S_c;
  Double I_ERI_D2y_Py_Fy2z_S_c = I_ERI_F3y_S_Fy2z_S_c+ABY*I_ERI_D2y_S_Fy2z_S_c;
  Double I_ERI_Dyz_Py_Fy2z_S_c = I_ERI_F2yz_S_Fy2z_S_c+ABY*I_ERI_Dyz_S_Fy2z_S_c;
  Double I_ERI_D2z_Py_Fy2z_S_c = I_ERI_Fy2z_S_Fy2z_S_c+ABY*I_ERI_D2z_S_Fy2z_S_c;
  Double I_ERI_D2x_Pz_Fy2z_S_c = I_ERI_F2xz_S_Fy2z_S_c+ABZ*I_ERI_D2x_S_Fy2z_S_c;
  Double I_ERI_Dxy_Pz_Fy2z_S_c = I_ERI_Fxyz_S_Fy2z_S_c+ABZ*I_ERI_Dxy_S_Fy2z_S_c;
  Double I_ERI_Dxz_Pz_Fy2z_S_c = I_ERI_Fx2z_S_Fy2z_S_c+ABZ*I_ERI_Dxz_S_Fy2z_S_c;
  Double I_ERI_D2y_Pz_Fy2z_S_c = I_ERI_F2yz_S_Fy2z_S_c+ABZ*I_ERI_D2y_S_Fy2z_S_c;
  Double I_ERI_Dyz_Pz_Fy2z_S_c = I_ERI_Fy2z_S_Fy2z_S_c+ABZ*I_ERI_Dyz_S_Fy2z_S_c;
  Double I_ERI_D2z_Pz_Fy2z_S_c = I_ERI_F3z_S_Fy2z_S_c+ABZ*I_ERI_D2z_S_Fy2z_S_c;
  Double I_ERI_D2x_Px_F3z_S_c = I_ERI_F3x_S_F3z_S_c+ABX*I_ERI_D2x_S_F3z_S_c;
  Double I_ERI_Dxy_Px_F3z_S_c = I_ERI_F2xy_S_F3z_S_c+ABX*I_ERI_Dxy_S_F3z_S_c;
  Double I_ERI_Dxz_Px_F3z_S_c = I_ERI_F2xz_S_F3z_S_c+ABX*I_ERI_Dxz_S_F3z_S_c;
  Double I_ERI_D2y_Px_F3z_S_c = I_ERI_Fx2y_S_F3z_S_c+ABX*I_ERI_D2y_S_F3z_S_c;
  Double I_ERI_Dyz_Px_F3z_S_c = I_ERI_Fxyz_S_F3z_S_c+ABX*I_ERI_Dyz_S_F3z_S_c;
  Double I_ERI_D2z_Px_F3z_S_c = I_ERI_Fx2z_S_F3z_S_c+ABX*I_ERI_D2z_S_F3z_S_c;
  Double I_ERI_D2x_Py_F3z_S_c = I_ERI_F2xy_S_F3z_S_c+ABY*I_ERI_D2x_S_F3z_S_c;
  Double I_ERI_Dxy_Py_F3z_S_c = I_ERI_Fx2y_S_F3z_S_c+ABY*I_ERI_Dxy_S_F3z_S_c;
  Double I_ERI_Dxz_Py_F3z_S_c = I_ERI_Fxyz_S_F3z_S_c+ABY*I_ERI_Dxz_S_F3z_S_c;
  Double I_ERI_D2y_Py_F3z_S_c = I_ERI_F3y_S_F3z_S_c+ABY*I_ERI_D2y_S_F3z_S_c;
  Double I_ERI_Dyz_Py_F3z_S_c = I_ERI_F2yz_S_F3z_S_c+ABY*I_ERI_Dyz_S_F3z_S_c;
  Double I_ERI_D2z_Py_F3z_S_c = I_ERI_Fy2z_S_F3z_S_c+ABY*I_ERI_D2z_S_F3z_S_c;
  Double I_ERI_D2x_Pz_F3z_S_c = I_ERI_F2xz_S_F3z_S_c+ABZ*I_ERI_D2x_S_F3z_S_c;
  Double I_ERI_Dxy_Pz_F3z_S_c = I_ERI_Fxyz_S_F3z_S_c+ABZ*I_ERI_Dxy_S_F3z_S_c;
  Double I_ERI_Dxz_Pz_F3z_S_c = I_ERI_Fx2z_S_F3z_S_c+ABZ*I_ERI_Dxz_S_F3z_S_c;
  Double I_ERI_D2y_Pz_F3z_S_c = I_ERI_F2yz_S_F3z_S_c+ABZ*I_ERI_D2y_S_F3z_S_c;
  Double I_ERI_Dyz_Pz_F3z_S_c = I_ERI_Fy2z_S_F3z_S_c+ABZ*I_ERI_Dyz_S_F3z_S_c;
  Double I_ERI_D2z_Pz_F3z_S_c = I_ERI_F3z_S_F3z_S_c+ABZ*I_ERI_D2z_S_F3z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_F_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 50 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_c
   * RHS shell quartet name: SQ_ERI_F_S_F_S_c
   ************************************************************/
  Double I_ERI_F3x_Px_F3x_S_c = I_ERI_G4x_S_F3x_S_c+ABX*I_ERI_F3x_S_F3x_S_c;
  Double I_ERI_F2xy_Px_F3x_S_c = I_ERI_G3xy_S_F3x_S_c+ABX*I_ERI_F2xy_S_F3x_S_c;
  Double I_ERI_F2xz_Px_F3x_S_c = I_ERI_G3xz_S_F3x_S_c+ABX*I_ERI_F2xz_S_F3x_S_c;
  Double I_ERI_Fx2y_Px_F3x_S_c = I_ERI_G2x2y_S_F3x_S_c+ABX*I_ERI_Fx2y_S_F3x_S_c;
  Double I_ERI_Fxyz_Px_F3x_S_c = I_ERI_G2xyz_S_F3x_S_c+ABX*I_ERI_Fxyz_S_F3x_S_c;
  Double I_ERI_Fx2z_Px_F3x_S_c = I_ERI_G2x2z_S_F3x_S_c+ABX*I_ERI_Fx2z_S_F3x_S_c;
  Double I_ERI_F3y_Px_F3x_S_c = I_ERI_Gx3y_S_F3x_S_c+ABX*I_ERI_F3y_S_F3x_S_c;
  Double I_ERI_F2yz_Px_F3x_S_c = I_ERI_Gx2yz_S_F3x_S_c+ABX*I_ERI_F2yz_S_F3x_S_c;
  Double I_ERI_Fy2z_Px_F3x_S_c = I_ERI_Gxy2z_S_F3x_S_c+ABX*I_ERI_Fy2z_S_F3x_S_c;
  Double I_ERI_F3z_Px_F3x_S_c = I_ERI_Gx3z_S_F3x_S_c+ABX*I_ERI_F3z_S_F3x_S_c;
  Double I_ERI_F2xy_Py_F3x_S_c = I_ERI_G2x2y_S_F3x_S_c+ABY*I_ERI_F2xy_S_F3x_S_c;
  Double I_ERI_F2xz_Py_F3x_S_c = I_ERI_G2xyz_S_F3x_S_c+ABY*I_ERI_F2xz_S_F3x_S_c;
  Double I_ERI_Fx2y_Py_F3x_S_c = I_ERI_Gx3y_S_F3x_S_c+ABY*I_ERI_Fx2y_S_F3x_S_c;
  Double I_ERI_Fxyz_Py_F3x_S_c = I_ERI_Gx2yz_S_F3x_S_c+ABY*I_ERI_Fxyz_S_F3x_S_c;
  Double I_ERI_Fx2z_Py_F3x_S_c = I_ERI_Gxy2z_S_F3x_S_c+ABY*I_ERI_Fx2z_S_F3x_S_c;
  Double I_ERI_F3y_Py_F3x_S_c = I_ERI_G4y_S_F3x_S_c+ABY*I_ERI_F3y_S_F3x_S_c;
  Double I_ERI_F2yz_Py_F3x_S_c = I_ERI_G3yz_S_F3x_S_c+ABY*I_ERI_F2yz_S_F3x_S_c;
  Double I_ERI_Fy2z_Py_F3x_S_c = I_ERI_G2y2z_S_F3x_S_c+ABY*I_ERI_Fy2z_S_F3x_S_c;
  Double I_ERI_F3z_Py_F3x_S_c = I_ERI_Gy3z_S_F3x_S_c+ABY*I_ERI_F3z_S_F3x_S_c;
  Double I_ERI_F2xz_Pz_F3x_S_c = I_ERI_G2x2z_S_F3x_S_c+ABZ*I_ERI_F2xz_S_F3x_S_c;
  Double I_ERI_Fxyz_Pz_F3x_S_c = I_ERI_Gxy2z_S_F3x_S_c+ABZ*I_ERI_Fxyz_S_F3x_S_c;
  Double I_ERI_Fx2z_Pz_F3x_S_c = I_ERI_Gx3z_S_F3x_S_c+ABZ*I_ERI_Fx2z_S_F3x_S_c;
  Double I_ERI_F2yz_Pz_F3x_S_c = I_ERI_G2y2z_S_F3x_S_c+ABZ*I_ERI_F2yz_S_F3x_S_c;
  Double I_ERI_Fy2z_Pz_F3x_S_c = I_ERI_Gy3z_S_F3x_S_c+ABZ*I_ERI_Fy2z_S_F3x_S_c;
  Double I_ERI_F3z_Pz_F3x_S_c = I_ERI_G4z_S_F3x_S_c+ABZ*I_ERI_F3z_S_F3x_S_c;
  Double I_ERI_F3x_Px_F2xy_S_c = I_ERI_G4x_S_F2xy_S_c+ABX*I_ERI_F3x_S_F2xy_S_c;
  Double I_ERI_F2xy_Px_F2xy_S_c = I_ERI_G3xy_S_F2xy_S_c+ABX*I_ERI_F2xy_S_F2xy_S_c;
  Double I_ERI_F2xz_Px_F2xy_S_c = I_ERI_G3xz_S_F2xy_S_c+ABX*I_ERI_F2xz_S_F2xy_S_c;
  Double I_ERI_Fx2y_Px_F2xy_S_c = I_ERI_G2x2y_S_F2xy_S_c+ABX*I_ERI_Fx2y_S_F2xy_S_c;
  Double I_ERI_Fxyz_Px_F2xy_S_c = I_ERI_G2xyz_S_F2xy_S_c+ABX*I_ERI_Fxyz_S_F2xy_S_c;
  Double I_ERI_Fx2z_Px_F2xy_S_c = I_ERI_G2x2z_S_F2xy_S_c+ABX*I_ERI_Fx2z_S_F2xy_S_c;
  Double I_ERI_F3y_Px_F2xy_S_c = I_ERI_Gx3y_S_F2xy_S_c+ABX*I_ERI_F3y_S_F2xy_S_c;
  Double I_ERI_F2yz_Px_F2xy_S_c = I_ERI_Gx2yz_S_F2xy_S_c+ABX*I_ERI_F2yz_S_F2xy_S_c;
  Double I_ERI_Fy2z_Px_F2xy_S_c = I_ERI_Gxy2z_S_F2xy_S_c+ABX*I_ERI_Fy2z_S_F2xy_S_c;
  Double I_ERI_F3z_Px_F2xy_S_c = I_ERI_Gx3z_S_F2xy_S_c+ABX*I_ERI_F3z_S_F2xy_S_c;
  Double I_ERI_F2xy_Py_F2xy_S_c = I_ERI_G2x2y_S_F2xy_S_c+ABY*I_ERI_F2xy_S_F2xy_S_c;
  Double I_ERI_F2xz_Py_F2xy_S_c = I_ERI_G2xyz_S_F2xy_S_c+ABY*I_ERI_F2xz_S_F2xy_S_c;
  Double I_ERI_Fx2y_Py_F2xy_S_c = I_ERI_Gx3y_S_F2xy_S_c+ABY*I_ERI_Fx2y_S_F2xy_S_c;
  Double I_ERI_Fxyz_Py_F2xy_S_c = I_ERI_Gx2yz_S_F2xy_S_c+ABY*I_ERI_Fxyz_S_F2xy_S_c;
  Double I_ERI_Fx2z_Py_F2xy_S_c = I_ERI_Gxy2z_S_F2xy_S_c+ABY*I_ERI_Fx2z_S_F2xy_S_c;
  Double I_ERI_F3y_Py_F2xy_S_c = I_ERI_G4y_S_F2xy_S_c+ABY*I_ERI_F3y_S_F2xy_S_c;
  Double I_ERI_F2yz_Py_F2xy_S_c = I_ERI_G3yz_S_F2xy_S_c+ABY*I_ERI_F2yz_S_F2xy_S_c;
  Double I_ERI_Fy2z_Py_F2xy_S_c = I_ERI_G2y2z_S_F2xy_S_c+ABY*I_ERI_Fy2z_S_F2xy_S_c;
  Double I_ERI_F3z_Py_F2xy_S_c = I_ERI_Gy3z_S_F2xy_S_c+ABY*I_ERI_F3z_S_F2xy_S_c;
  Double I_ERI_F2xz_Pz_F2xy_S_c = I_ERI_G2x2z_S_F2xy_S_c+ABZ*I_ERI_F2xz_S_F2xy_S_c;
  Double I_ERI_Fxyz_Pz_F2xy_S_c = I_ERI_Gxy2z_S_F2xy_S_c+ABZ*I_ERI_Fxyz_S_F2xy_S_c;
  Double I_ERI_Fx2z_Pz_F2xy_S_c = I_ERI_Gx3z_S_F2xy_S_c+ABZ*I_ERI_Fx2z_S_F2xy_S_c;
  Double I_ERI_F2yz_Pz_F2xy_S_c = I_ERI_G2y2z_S_F2xy_S_c+ABZ*I_ERI_F2yz_S_F2xy_S_c;
  Double I_ERI_Fy2z_Pz_F2xy_S_c = I_ERI_Gy3z_S_F2xy_S_c+ABZ*I_ERI_Fy2z_S_F2xy_S_c;
  Double I_ERI_F3z_Pz_F2xy_S_c = I_ERI_G4z_S_F2xy_S_c+ABZ*I_ERI_F3z_S_F2xy_S_c;
  Double I_ERI_F3x_Px_F2xz_S_c = I_ERI_G4x_S_F2xz_S_c+ABX*I_ERI_F3x_S_F2xz_S_c;
  Double I_ERI_F2xy_Px_F2xz_S_c = I_ERI_G3xy_S_F2xz_S_c+ABX*I_ERI_F2xy_S_F2xz_S_c;
  Double I_ERI_F2xz_Px_F2xz_S_c = I_ERI_G3xz_S_F2xz_S_c+ABX*I_ERI_F2xz_S_F2xz_S_c;
  Double I_ERI_Fx2y_Px_F2xz_S_c = I_ERI_G2x2y_S_F2xz_S_c+ABX*I_ERI_Fx2y_S_F2xz_S_c;
  Double I_ERI_Fxyz_Px_F2xz_S_c = I_ERI_G2xyz_S_F2xz_S_c+ABX*I_ERI_Fxyz_S_F2xz_S_c;
  Double I_ERI_Fx2z_Px_F2xz_S_c = I_ERI_G2x2z_S_F2xz_S_c+ABX*I_ERI_Fx2z_S_F2xz_S_c;
  Double I_ERI_F3y_Px_F2xz_S_c = I_ERI_Gx3y_S_F2xz_S_c+ABX*I_ERI_F3y_S_F2xz_S_c;
  Double I_ERI_F2yz_Px_F2xz_S_c = I_ERI_Gx2yz_S_F2xz_S_c+ABX*I_ERI_F2yz_S_F2xz_S_c;
  Double I_ERI_Fy2z_Px_F2xz_S_c = I_ERI_Gxy2z_S_F2xz_S_c+ABX*I_ERI_Fy2z_S_F2xz_S_c;
  Double I_ERI_F3z_Px_F2xz_S_c = I_ERI_Gx3z_S_F2xz_S_c+ABX*I_ERI_F3z_S_F2xz_S_c;
  Double I_ERI_F2xy_Py_F2xz_S_c = I_ERI_G2x2y_S_F2xz_S_c+ABY*I_ERI_F2xy_S_F2xz_S_c;
  Double I_ERI_F2xz_Py_F2xz_S_c = I_ERI_G2xyz_S_F2xz_S_c+ABY*I_ERI_F2xz_S_F2xz_S_c;
  Double I_ERI_Fx2y_Py_F2xz_S_c = I_ERI_Gx3y_S_F2xz_S_c+ABY*I_ERI_Fx2y_S_F2xz_S_c;
  Double I_ERI_Fxyz_Py_F2xz_S_c = I_ERI_Gx2yz_S_F2xz_S_c+ABY*I_ERI_Fxyz_S_F2xz_S_c;
  Double I_ERI_Fx2z_Py_F2xz_S_c = I_ERI_Gxy2z_S_F2xz_S_c+ABY*I_ERI_Fx2z_S_F2xz_S_c;
  Double I_ERI_F3y_Py_F2xz_S_c = I_ERI_G4y_S_F2xz_S_c+ABY*I_ERI_F3y_S_F2xz_S_c;
  Double I_ERI_F2yz_Py_F2xz_S_c = I_ERI_G3yz_S_F2xz_S_c+ABY*I_ERI_F2yz_S_F2xz_S_c;
  Double I_ERI_Fy2z_Py_F2xz_S_c = I_ERI_G2y2z_S_F2xz_S_c+ABY*I_ERI_Fy2z_S_F2xz_S_c;
  Double I_ERI_F3z_Py_F2xz_S_c = I_ERI_Gy3z_S_F2xz_S_c+ABY*I_ERI_F3z_S_F2xz_S_c;
  Double I_ERI_F2xz_Pz_F2xz_S_c = I_ERI_G2x2z_S_F2xz_S_c+ABZ*I_ERI_F2xz_S_F2xz_S_c;
  Double I_ERI_Fxyz_Pz_F2xz_S_c = I_ERI_Gxy2z_S_F2xz_S_c+ABZ*I_ERI_Fxyz_S_F2xz_S_c;
  Double I_ERI_Fx2z_Pz_F2xz_S_c = I_ERI_Gx3z_S_F2xz_S_c+ABZ*I_ERI_Fx2z_S_F2xz_S_c;
  Double I_ERI_F2yz_Pz_F2xz_S_c = I_ERI_G2y2z_S_F2xz_S_c+ABZ*I_ERI_F2yz_S_F2xz_S_c;
  Double I_ERI_Fy2z_Pz_F2xz_S_c = I_ERI_Gy3z_S_F2xz_S_c+ABZ*I_ERI_Fy2z_S_F2xz_S_c;
  Double I_ERI_F3z_Pz_F2xz_S_c = I_ERI_G4z_S_F2xz_S_c+ABZ*I_ERI_F3z_S_F2xz_S_c;
  Double I_ERI_F3x_Px_Fx2y_S_c = I_ERI_G4x_S_Fx2y_S_c+ABX*I_ERI_F3x_S_Fx2y_S_c;
  Double I_ERI_F2xy_Px_Fx2y_S_c = I_ERI_G3xy_S_Fx2y_S_c+ABX*I_ERI_F2xy_S_Fx2y_S_c;
  Double I_ERI_F2xz_Px_Fx2y_S_c = I_ERI_G3xz_S_Fx2y_S_c+ABX*I_ERI_F2xz_S_Fx2y_S_c;
  Double I_ERI_Fx2y_Px_Fx2y_S_c = I_ERI_G2x2y_S_Fx2y_S_c+ABX*I_ERI_Fx2y_S_Fx2y_S_c;
  Double I_ERI_Fxyz_Px_Fx2y_S_c = I_ERI_G2xyz_S_Fx2y_S_c+ABX*I_ERI_Fxyz_S_Fx2y_S_c;
  Double I_ERI_Fx2z_Px_Fx2y_S_c = I_ERI_G2x2z_S_Fx2y_S_c+ABX*I_ERI_Fx2z_S_Fx2y_S_c;
  Double I_ERI_F3y_Px_Fx2y_S_c = I_ERI_Gx3y_S_Fx2y_S_c+ABX*I_ERI_F3y_S_Fx2y_S_c;
  Double I_ERI_F2yz_Px_Fx2y_S_c = I_ERI_Gx2yz_S_Fx2y_S_c+ABX*I_ERI_F2yz_S_Fx2y_S_c;
  Double I_ERI_Fy2z_Px_Fx2y_S_c = I_ERI_Gxy2z_S_Fx2y_S_c+ABX*I_ERI_Fy2z_S_Fx2y_S_c;
  Double I_ERI_F3z_Px_Fx2y_S_c = I_ERI_Gx3z_S_Fx2y_S_c+ABX*I_ERI_F3z_S_Fx2y_S_c;
  Double I_ERI_F2xy_Py_Fx2y_S_c = I_ERI_G2x2y_S_Fx2y_S_c+ABY*I_ERI_F2xy_S_Fx2y_S_c;
  Double I_ERI_F2xz_Py_Fx2y_S_c = I_ERI_G2xyz_S_Fx2y_S_c+ABY*I_ERI_F2xz_S_Fx2y_S_c;
  Double I_ERI_Fx2y_Py_Fx2y_S_c = I_ERI_Gx3y_S_Fx2y_S_c+ABY*I_ERI_Fx2y_S_Fx2y_S_c;
  Double I_ERI_Fxyz_Py_Fx2y_S_c = I_ERI_Gx2yz_S_Fx2y_S_c+ABY*I_ERI_Fxyz_S_Fx2y_S_c;
  Double I_ERI_Fx2z_Py_Fx2y_S_c = I_ERI_Gxy2z_S_Fx2y_S_c+ABY*I_ERI_Fx2z_S_Fx2y_S_c;
  Double I_ERI_F3y_Py_Fx2y_S_c = I_ERI_G4y_S_Fx2y_S_c+ABY*I_ERI_F3y_S_Fx2y_S_c;
  Double I_ERI_F2yz_Py_Fx2y_S_c = I_ERI_G3yz_S_Fx2y_S_c+ABY*I_ERI_F2yz_S_Fx2y_S_c;
  Double I_ERI_Fy2z_Py_Fx2y_S_c = I_ERI_G2y2z_S_Fx2y_S_c+ABY*I_ERI_Fy2z_S_Fx2y_S_c;
  Double I_ERI_F3z_Py_Fx2y_S_c = I_ERI_Gy3z_S_Fx2y_S_c+ABY*I_ERI_F3z_S_Fx2y_S_c;
  Double I_ERI_F2xz_Pz_Fx2y_S_c = I_ERI_G2x2z_S_Fx2y_S_c+ABZ*I_ERI_F2xz_S_Fx2y_S_c;
  Double I_ERI_Fxyz_Pz_Fx2y_S_c = I_ERI_Gxy2z_S_Fx2y_S_c+ABZ*I_ERI_Fxyz_S_Fx2y_S_c;
  Double I_ERI_Fx2z_Pz_Fx2y_S_c = I_ERI_Gx3z_S_Fx2y_S_c+ABZ*I_ERI_Fx2z_S_Fx2y_S_c;
  Double I_ERI_F2yz_Pz_Fx2y_S_c = I_ERI_G2y2z_S_Fx2y_S_c+ABZ*I_ERI_F2yz_S_Fx2y_S_c;
  Double I_ERI_Fy2z_Pz_Fx2y_S_c = I_ERI_Gy3z_S_Fx2y_S_c+ABZ*I_ERI_Fy2z_S_Fx2y_S_c;
  Double I_ERI_F3z_Pz_Fx2y_S_c = I_ERI_G4z_S_Fx2y_S_c+ABZ*I_ERI_F3z_S_Fx2y_S_c;
  Double I_ERI_F3x_Px_Fxyz_S_c = I_ERI_G4x_S_Fxyz_S_c+ABX*I_ERI_F3x_S_Fxyz_S_c;
  Double I_ERI_F2xy_Px_Fxyz_S_c = I_ERI_G3xy_S_Fxyz_S_c+ABX*I_ERI_F2xy_S_Fxyz_S_c;
  Double I_ERI_F2xz_Px_Fxyz_S_c = I_ERI_G3xz_S_Fxyz_S_c+ABX*I_ERI_F2xz_S_Fxyz_S_c;
  Double I_ERI_Fx2y_Px_Fxyz_S_c = I_ERI_G2x2y_S_Fxyz_S_c+ABX*I_ERI_Fx2y_S_Fxyz_S_c;
  Double I_ERI_Fxyz_Px_Fxyz_S_c = I_ERI_G2xyz_S_Fxyz_S_c+ABX*I_ERI_Fxyz_S_Fxyz_S_c;
  Double I_ERI_Fx2z_Px_Fxyz_S_c = I_ERI_G2x2z_S_Fxyz_S_c+ABX*I_ERI_Fx2z_S_Fxyz_S_c;
  Double I_ERI_F3y_Px_Fxyz_S_c = I_ERI_Gx3y_S_Fxyz_S_c+ABX*I_ERI_F3y_S_Fxyz_S_c;
  Double I_ERI_F2yz_Px_Fxyz_S_c = I_ERI_Gx2yz_S_Fxyz_S_c+ABX*I_ERI_F2yz_S_Fxyz_S_c;
  Double I_ERI_Fy2z_Px_Fxyz_S_c = I_ERI_Gxy2z_S_Fxyz_S_c+ABX*I_ERI_Fy2z_S_Fxyz_S_c;
  Double I_ERI_F3z_Px_Fxyz_S_c = I_ERI_Gx3z_S_Fxyz_S_c+ABX*I_ERI_F3z_S_Fxyz_S_c;
  Double I_ERI_F2xy_Py_Fxyz_S_c = I_ERI_G2x2y_S_Fxyz_S_c+ABY*I_ERI_F2xy_S_Fxyz_S_c;
  Double I_ERI_F2xz_Py_Fxyz_S_c = I_ERI_G2xyz_S_Fxyz_S_c+ABY*I_ERI_F2xz_S_Fxyz_S_c;
  Double I_ERI_Fx2y_Py_Fxyz_S_c = I_ERI_Gx3y_S_Fxyz_S_c+ABY*I_ERI_Fx2y_S_Fxyz_S_c;
  Double I_ERI_Fxyz_Py_Fxyz_S_c = I_ERI_Gx2yz_S_Fxyz_S_c+ABY*I_ERI_Fxyz_S_Fxyz_S_c;
  Double I_ERI_Fx2z_Py_Fxyz_S_c = I_ERI_Gxy2z_S_Fxyz_S_c+ABY*I_ERI_Fx2z_S_Fxyz_S_c;
  Double I_ERI_F3y_Py_Fxyz_S_c = I_ERI_G4y_S_Fxyz_S_c+ABY*I_ERI_F3y_S_Fxyz_S_c;
  Double I_ERI_F2yz_Py_Fxyz_S_c = I_ERI_G3yz_S_Fxyz_S_c+ABY*I_ERI_F2yz_S_Fxyz_S_c;
  Double I_ERI_Fy2z_Py_Fxyz_S_c = I_ERI_G2y2z_S_Fxyz_S_c+ABY*I_ERI_Fy2z_S_Fxyz_S_c;
  Double I_ERI_F3z_Py_Fxyz_S_c = I_ERI_Gy3z_S_Fxyz_S_c+ABY*I_ERI_F3z_S_Fxyz_S_c;
  Double I_ERI_F2xz_Pz_Fxyz_S_c = I_ERI_G2x2z_S_Fxyz_S_c+ABZ*I_ERI_F2xz_S_Fxyz_S_c;
  Double I_ERI_Fxyz_Pz_Fxyz_S_c = I_ERI_Gxy2z_S_Fxyz_S_c+ABZ*I_ERI_Fxyz_S_Fxyz_S_c;
  Double I_ERI_Fx2z_Pz_Fxyz_S_c = I_ERI_Gx3z_S_Fxyz_S_c+ABZ*I_ERI_Fx2z_S_Fxyz_S_c;
  Double I_ERI_F2yz_Pz_Fxyz_S_c = I_ERI_G2y2z_S_Fxyz_S_c+ABZ*I_ERI_F2yz_S_Fxyz_S_c;
  Double I_ERI_Fy2z_Pz_Fxyz_S_c = I_ERI_Gy3z_S_Fxyz_S_c+ABZ*I_ERI_Fy2z_S_Fxyz_S_c;
  Double I_ERI_F3z_Pz_Fxyz_S_c = I_ERI_G4z_S_Fxyz_S_c+ABZ*I_ERI_F3z_S_Fxyz_S_c;
  Double I_ERI_F3x_Px_Fx2z_S_c = I_ERI_G4x_S_Fx2z_S_c+ABX*I_ERI_F3x_S_Fx2z_S_c;
  Double I_ERI_F2xy_Px_Fx2z_S_c = I_ERI_G3xy_S_Fx2z_S_c+ABX*I_ERI_F2xy_S_Fx2z_S_c;
  Double I_ERI_F2xz_Px_Fx2z_S_c = I_ERI_G3xz_S_Fx2z_S_c+ABX*I_ERI_F2xz_S_Fx2z_S_c;
  Double I_ERI_Fx2y_Px_Fx2z_S_c = I_ERI_G2x2y_S_Fx2z_S_c+ABX*I_ERI_Fx2y_S_Fx2z_S_c;
  Double I_ERI_Fxyz_Px_Fx2z_S_c = I_ERI_G2xyz_S_Fx2z_S_c+ABX*I_ERI_Fxyz_S_Fx2z_S_c;
  Double I_ERI_Fx2z_Px_Fx2z_S_c = I_ERI_G2x2z_S_Fx2z_S_c+ABX*I_ERI_Fx2z_S_Fx2z_S_c;
  Double I_ERI_F3y_Px_Fx2z_S_c = I_ERI_Gx3y_S_Fx2z_S_c+ABX*I_ERI_F3y_S_Fx2z_S_c;
  Double I_ERI_F2yz_Px_Fx2z_S_c = I_ERI_Gx2yz_S_Fx2z_S_c+ABX*I_ERI_F2yz_S_Fx2z_S_c;
  Double I_ERI_Fy2z_Px_Fx2z_S_c = I_ERI_Gxy2z_S_Fx2z_S_c+ABX*I_ERI_Fy2z_S_Fx2z_S_c;
  Double I_ERI_F3z_Px_Fx2z_S_c = I_ERI_Gx3z_S_Fx2z_S_c+ABX*I_ERI_F3z_S_Fx2z_S_c;
  Double I_ERI_F2xy_Py_Fx2z_S_c = I_ERI_G2x2y_S_Fx2z_S_c+ABY*I_ERI_F2xy_S_Fx2z_S_c;
  Double I_ERI_F2xz_Py_Fx2z_S_c = I_ERI_G2xyz_S_Fx2z_S_c+ABY*I_ERI_F2xz_S_Fx2z_S_c;
  Double I_ERI_Fx2y_Py_Fx2z_S_c = I_ERI_Gx3y_S_Fx2z_S_c+ABY*I_ERI_Fx2y_S_Fx2z_S_c;
  Double I_ERI_Fxyz_Py_Fx2z_S_c = I_ERI_Gx2yz_S_Fx2z_S_c+ABY*I_ERI_Fxyz_S_Fx2z_S_c;
  Double I_ERI_Fx2z_Py_Fx2z_S_c = I_ERI_Gxy2z_S_Fx2z_S_c+ABY*I_ERI_Fx2z_S_Fx2z_S_c;
  Double I_ERI_F3y_Py_Fx2z_S_c = I_ERI_G4y_S_Fx2z_S_c+ABY*I_ERI_F3y_S_Fx2z_S_c;
  Double I_ERI_F2yz_Py_Fx2z_S_c = I_ERI_G3yz_S_Fx2z_S_c+ABY*I_ERI_F2yz_S_Fx2z_S_c;
  Double I_ERI_Fy2z_Py_Fx2z_S_c = I_ERI_G2y2z_S_Fx2z_S_c+ABY*I_ERI_Fy2z_S_Fx2z_S_c;
  Double I_ERI_F3z_Py_Fx2z_S_c = I_ERI_Gy3z_S_Fx2z_S_c+ABY*I_ERI_F3z_S_Fx2z_S_c;
  Double I_ERI_F2xz_Pz_Fx2z_S_c = I_ERI_G2x2z_S_Fx2z_S_c+ABZ*I_ERI_F2xz_S_Fx2z_S_c;
  Double I_ERI_Fxyz_Pz_Fx2z_S_c = I_ERI_Gxy2z_S_Fx2z_S_c+ABZ*I_ERI_Fxyz_S_Fx2z_S_c;
  Double I_ERI_Fx2z_Pz_Fx2z_S_c = I_ERI_Gx3z_S_Fx2z_S_c+ABZ*I_ERI_Fx2z_S_Fx2z_S_c;
  Double I_ERI_F2yz_Pz_Fx2z_S_c = I_ERI_G2y2z_S_Fx2z_S_c+ABZ*I_ERI_F2yz_S_Fx2z_S_c;
  Double I_ERI_Fy2z_Pz_Fx2z_S_c = I_ERI_Gy3z_S_Fx2z_S_c+ABZ*I_ERI_Fy2z_S_Fx2z_S_c;
  Double I_ERI_F3z_Pz_Fx2z_S_c = I_ERI_G4z_S_Fx2z_S_c+ABZ*I_ERI_F3z_S_Fx2z_S_c;
  Double I_ERI_F3x_Px_F3y_S_c = I_ERI_G4x_S_F3y_S_c+ABX*I_ERI_F3x_S_F3y_S_c;
  Double I_ERI_F2xy_Px_F3y_S_c = I_ERI_G3xy_S_F3y_S_c+ABX*I_ERI_F2xy_S_F3y_S_c;
  Double I_ERI_F2xz_Px_F3y_S_c = I_ERI_G3xz_S_F3y_S_c+ABX*I_ERI_F2xz_S_F3y_S_c;
  Double I_ERI_Fx2y_Px_F3y_S_c = I_ERI_G2x2y_S_F3y_S_c+ABX*I_ERI_Fx2y_S_F3y_S_c;
  Double I_ERI_Fxyz_Px_F3y_S_c = I_ERI_G2xyz_S_F3y_S_c+ABX*I_ERI_Fxyz_S_F3y_S_c;
  Double I_ERI_Fx2z_Px_F3y_S_c = I_ERI_G2x2z_S_F3y_S_c+ABX*I_ERI_Fx2z_S_F3y_S_c;
  Double I_ERI_F3y_Px_F3y_S_c = I_ERI_Gx3y_S_F3y_S_c+ABX*I_ERI_F3y_S_F3y_S_c;
  Double I_ERI_F2yz_Px_F3y_S_c = I_ERI_Gx2yz_S_F3y_S_c+ABX*I_ERI_F2yz_S_F3y_S_c;
  Double I_ERI_Fy2z_Px_F3y_S_c = I_ERI_Gxy2z_S_F3y_S_c+ABX*I_ERI_Fy2z_S_F3y_S_c;
  Double I_ERI_F3z_Px_F3y_S_c = I_ERI_Gx3z_S_F3y_S_c+ABX*I_ERI_F3z_S_F3y_S_c;
  Double I_ERI_F2xy_Py_F3y_S_c = I_ERI_G2x2y_S_F3y_S_c+ABY*I_ERI_F2xy_S_F3y_S_c;
  Double I_ERI_F2xz_Py_F3y_S_c = I_ERI_G2xyz_S_F3y_S_c+ABY*I_ERI_F2xz_S_F3y_S_c;
  Double I_ERI_Fx2y_Py_F3y_S_c = I_ERI_Gx3y_S_F3y_S_c+ABY*I_ERI_Fx2y_S_F3y_S_c;
  Double I_ERI_Fxyz_Py_F3y_S_c = I_ERI_Gx2yz_S_F3y_S_c+ABY*I_ERI_Fxyz_S_F3y_S_c;
  Double I_ERI_Fx2z_Py_F3y_S_c = I_ERI_Gxy2z_S_F3y_S_c+ABY*I_ERI_Fx2z_S_F3y_S_c;
  Double I_ERI_F3y_Py_F3y_S_c = I_ERI_G4y_S_F3y_S_c+ABY*I_ERI_F3y_S_F3y_S_c;
  Double I_ERI_F2yz_Py_F3y_S_c = I_ERI_G3yz_S_F3y_S_c+ABY*I_ERI_F2yz_S_F3y_S_c;
  Double I_ERI_Fy2z_Py_F3y_S_c = I_ERI_G2y2z_S_F3y_S_c+ABY*I_ERI_Fy2z_S_F3y_S_c;
  Double I_ERI_F3z_Py_F3y_S_c = I_ERI_Gy3z_S_F3y_S_c+ABY*I_ERI_F3z_S_F3y_S_c;
  Double I_ERI_F2xz_Pz_F3y_S_c = I_ERI_G2x2z_S_F3y_S_c+ABZ*I_ERI_F2xz_S_F3y_S_c;
  Double I_ERI_Fxyz_Pz_F3y_S_c = I_ERI_Gxy2z_S_F3y_S_c+ABZ*I_ERI_Fxyz_S_F3y_S_c;
  Double I_ERI_Fx2z_Pz_F3y_S_c = I_ERI_Gx3z_S_F3y_S_c+ABZ*I_ERI_Fx2z_S_F3y_S_c;
  Double I_ERI_F2yz_Pz_F3y_S_c = I_ERI_G2y2z_S_F3y_S_c+ABZ*I_ERI_F2yz_S_F3y_S_c;
  Double I_ERI_Fy2z_Pz_F3y_S_c = I_ERI_Gy3z_S_F3y_S_c+ABZ*I_ERI_Fy2z_S_F3y_S_c;
  Double I_ERI_F3z_Pz_F3y_S_c = I_ERI_G4z_S_F3y_S_c+ABZ*I_ERI_F3z_S_F3y_S_c;
  Double I_ERI_F3x_Px_F2yz_S_c = I_ERI_G4x_S_F2yz_S_c+ABX*I_ERI_F3x_S_F2yz_S_c;
  Double I_ERI_F2xy_Px_F2yz_S_c = I_ERI_G3xy_S_F2yz_S_c+ABX*I_ERI_F2xy_S_F2yz_S_c;
  Double I_ERI_F2xz_Px_F2yz_S_c = I_ERI_G3xz_S_F2yz_S_c+ABX*I_ERI_F2xz_S_F2yz_S_c;
  Double I_ERI_Fx2y_Px_F2yz_S_c = I_ERI_G2x2y_S_F2yz_S_c+ABX*I_ERI_Fx2y_S_F2yz_S_c;
  Double I_ERI_Fxyz_Px_F2yz_S_c = I_ERI_G2xyz_S_F2yz_S_c+ABX*I_ERI_Fxyz_S_F2yz_S_c;
  Double I_ERI_Fx2z_Px_F2yz_S_c = I_ERI_G2x2z_S_F2yz_S_c+ABX*I_ERI_Fx2z_S_F2yz_S_c;
  Double I_ERI_F3y_Px_F2yz_S_c = I_ERI_Gx3y_S_F2yz_S_c+ABX*I_ERI_F3y_S_F2yz_S_c;
  Double I_ERI_F2yz_Px_F2yz_S_c = I_ERI_Gx2yz_S_F2yz_S_c+ABX*I_ERI_F2yz_S_F2yz_S_c;
  Double I_ERI_Fy2z_Px_F2yz_S_c = I_ERI_Gxy2z_S_F2yz_S_c+ABX*I_ERI_Fy2z_S_F2yz_S_c;
  Double I_ERI_F3z_Px_F2yz_S_c = I_ERI_Gx3z_S_F2yz_S_c+ABX*I_ERI_F3z_S_F2yz_S_c;
  Double I_ERI_F2xy_Py_F2yz_S_c = I_ERI_G2x2y_S_F2yz_S_c+ABY*I_ERI_F2xy_S_F2yz_S_c;
  Double I_ERI_F2xz_Py_F2yz_S_c = I_ERI_G2xyz_S_F2yz_S_c+ABY*I_ERI_F2xz_S_F2yz_S_c;
  Double I_ERI_Fx2y_Py_F2yz_S_c = I_ERI_Gx3y_S_F2yz_S_c+ABY*I_ERI_Fx2y_S_F2yz_S_c;
  Double I_ERI_Fxyz_Py_F2yz_S_c = I_ERI_Gx2yz_S_F2yz_S_c+ABY*I_ERI_Fxyz_S_F2yz_S_c;
  Double I_ERI_Fx2z_Py_F2yz_S_c = I_ERI_Gxy2z_S_F2yz_S_c+ABY*I_ERI_Fx2z_S_F2yz_S_c;
  Double I_ERI_F3y_Py_F2yz_S_c = I_ERI_G4y_S_F2yz_S_c+ABY*I_ERI_F3y_S_F2yz_S_c;
  Double I_ERI_F2yz_Py_F2yz_S_c = I_ERI_G3yz_S_F2yz_S_c+ABY*I_ERI_F2yz_S_F2yz_S_c;
  Double I_ERI_Fy2z_Py_F2yz_S_c = I_ERI_G2y2z_S_F2yz_S_c+ABY*I_ERI_Fy2z_S_F2yz_S_c;
  Double I_ERI_F3z_Py_F2yz_S_c = I_ERI_Gy3z_S_F2yz_S_c+ABY*I_ERI_F3z_S_F2yz_S_c;
  Double I_ERI_F2xz_Pz_F2yz_S_c = I_ERI_G2x2z_S_F2yz_S_c+ABZ*I_ERI_F2xz_S_F2yz_S_c;
  Double I_ERI_Fxyz_Pz_F2yz_S_c = I_ERI_Gxy2z_S_F2yz_S_c+ABZ*I_ERI_Fxyz_S_F2yz_S_c;
  Double I_ERI_Fx2z_Pz_F2yz_S_c = I_ERI_Gx3z_S_F2yz_S_c+ABZ*I_ERI_Fx2z_S_F2yz_S_c;
  Double I_ERI_F2yz_Pz_F2yz_S_c = I_ERI_G2y2z_S_F2yz_S_c+ABZ*I_ERI_F2yz_S_F2yz_S_c;
  Double I_ERI_Fy2z_Pz_F2yz_S_c = I_ERI_Gy3z_S_F2yz_S_c+ABZ*I_ERI_Fy2z_S_F2yz_S_c;
  Double I_ERI_F3z_Pz_F2yz_S_c = I_ERI_G4z_S_F2yz_S_c+ABZ*I_ERI_F3z_S_F2yz_S_c;
  Double I_ERI_F3x_Px_Fy2z_S_c = I_ERI_G4x_S_Fy2z_S_c+ABX*I_ERI_F3x_S_Fy2z_S_c;
  Double I_ERI_F2xy_Px_Fy2z_S_c = I_ERI_G3xy_S_Fy2z_S_c+ABX*I_ERI_F2xy_S_Fy2z_S_c;
  Double I_ERI_F2xz_Px_Fy2z_S_c = I_ERI_G3xz_S_Fy2z_S_c+ABX*I_ERI_F2xz_S_Fy2z_S_c;
  Double I_ERI_Fx2y_Px_Fy2z_S_c = I_ERI_G2x2y_S_Fy2z_S_c+ABX*I_ERI_Fx2y_S_Fy2z_S_c;
  Double I_ERI_Fxyz_Px_Fy2z_S_c = I_ERI_G2xyz_S_Fy2z_S_c+ABX*I_ERI_Fxyz_S_Fy2z_S_c;
  Double I_ERI_Fx2z_Px_Fy2z_S_c = I_ERI_G2x2z_S_Fy2z_S_c+ABX*I_ERI_Fx2z_S_Fy2z_S_c;
  Double I_ERI_F3y_Px_Fy2z_S_c = I_ERI_Gx3y_S_Fy2z_S_c+ABX*I_ERI_F3y_S_Fy2z_S_c;
  Double I_ERI_F2yz_Px_Fy2z_S_c = I_ERI_Gx2yz_S_Fy2z_S_c+ABX*I_ERI_F2yz_S_Fy2z_S_c;
  Double I_ERI_Fy2z_Px_Fy2z_S_c = I_ERI_Gxy2z_S_Fy2z_S_c+ABX*I_ERI_Fy2z_S_Fy2z_S_c;
  Double I_ERI_F3z_Px_Fy2z_S_c = I_ERI_Gx3z_S_Fy2z_S_c+ABX*I_ERI_F3z_S_Fy2z_S_c;
  Double I_ERI_F2xy_Py_Fy2z_S_c = I_ERI_G2x2y_S_Fy2z_S_c+ABY*I_ERI_F2xy_S_Fy2z_S_c;
  Double I_ERI_F2xz_Py_Fy2z_S_c = I_ERI_G2xyz_S_Fy2z_S_c+ABY*I_ERI_F2xz_S_Fy2z_S_c;
  Double I_ERI_Fx2y_Py_Fy2z_S_c = I_ERI_Gx3y_S_Fy2z_S_c+ABY*I_ERI_Fx2y_S_Fy2z_S_c;
  Double I_ERI_Fxyz_Py_Fy2z_S_c = I_ERI_Gx2yz_S_Fy2z_S_c+ABY*I_ERI_Fxyz_S_Fy2z_S_c;
  Double I_ERI_Fx2z_Py_Fy2z_S_c = I_ERI_Gxy2z_S_Fy2z_S_c+ABY*I_ERI_Fx2z_S_Fy2z_S_c;
  Double I_ERI_F3y_Py_Fy2z_S_c = I_ERI_G4y_S_Fy2z_S_c+ABY*I_ERI_F3y_S_Fy2z_S_c;
  Double I_ERI_F2yz_Py_Fy2z_S_c = I_ERI_G3yz_S_Fy2z_S_c+ABY*I_ERI_F2yz_S_Fy2z_S_c;
  Double I_ERI_Fy2z_Py_Fy2z_S_c = I_ERI_G2y2z_S_Fy2z_S_c+ABY*I_ERI_Fy2z_S_Fy2z_S_c;
  Double I_ERI_F3z_Py_Fy2z_S_c = I_ERI_Gy3z_S_Fy2z_S_c+ABY*I_ERI_F3z_S_Fy2z_S_c;
  Double I_ERI_F2xz_Pz_Fy2z_S_c = I_ERI_G2x2z_S_Fy2z_S_c+ABZ*I_ERI_F2xz_S_Fy2z_S_c;
  Double I_ERI_Fxyz_Pz_Fy2z_S_c = I_ERI_Gxy2z_S_Fy2z_S_c+ABZ*I_ERI_Fxyz_S_Fy2z_S_c;
  Double I_ERI_Fx2z_Pz_Fy2z_S_c = I_ERI_Gx3z_S_Fy2z_S_c+ABZ*I_ERI_Fx2z_S_Fy2z_S_c;
  Double I_ERI_F2yz_Pz_Fy2z_S_c = I_ERI_G2y2z_S_Fy2z_S_c+ABZ*I_ERI_F2yz_S_Fy2z_S_c;
  Double I_ERI_Fy2z_Pz_Fy2z_S_c = I_ERI_Gy3z_S_Fy2z_S_c+ABZ*I_ERI_Fy2z_S_Fy2z_S_c;
  Double I_ERI_F3z_Pz_Fy2z_S_c = I_ERI_G4z_S_Fy2z_S_c+ABZ*I_ERI_F3z_S_Fy2z_S_c;
  Double I_ERI_F3x_Px_F3z_S_c = I_ERI_G4x_S_F3z_S_c+ABX*I_ERI_F3x_S_F3z_S_c;
  Double I_ERI_F2xy_Px_F3z_S_c = I_ERI_G3xy_S_F3z_S_c+ABX*I_ERI_F2xy_S_F3z_S_c;
  Double I_ERI_F2xz_Px_F3z_S_c = I_ERI_G3xz_S_F3z_S_c+ABX*I_ERI_F2xz_S_F3z_S_c;
  Double I_ERI_Fx2y_Px_F3z_S_c = I_ERI_G2x2y_S_F3z_S_c+ABX*I_ERI_Fx2y_S_F3z_S_c;
  Double I_ERI_Fxyz_Px_F3z_S_c = I_ERI_G2xyz_S_F3z_S_c+ABX*I_ERI_Fxyz_S_F3z_S_c;
  Double I_ERI_Fx2z_Px_F3z_S_c = I_ERI_G2x2z_S_F3z_S_c+ABX*I_ERI_Fx2z_S_F3z_S_c;
  Double I_ERI_F3y_Px_F3z_S_c = I_ERI_Gx3y_S_F3z_S_c+ABX*I_ERI_F3y_S_F3z_S_c;
  Double I_ERI_F2yz_Px_F3z_S_c = I_ERI_Gx2yz_S_F3z_S_c+ABX*I_ERI_F2yz_S_F3z_S_c;
  Double I_ERI_Fy2z_Px_F3z_S_c = I_ERI_Gxy2z_S_F3z_S_c+ABX*I_ERI_Fy2z_S_F3z_S_c;
  Double I_ERI_F3z_Px_F3z_S_c = I_ERI_Gx3z_S_F3z_S_c+ABX*I_ERI_F3z_S_F3z_S_c;
  Double I_ERI_F2xy_Py_F3z_S_c = I_ERI_G2x2y_S_F3z_S_c+ABY*I_ERI_F2xy_S_F3z_S_c;
  Double I_ERI_F2xz_Py_F3z_S_c = I_ERI_G2xyz_S_F3z_S_c+ABY*I_ERI_F2xz_S_F3z_S_c;
  Double I_ERI_Fx2y_Py_F3z_S_c = I_ERI_Gx3y_S_F3z_S_c+ABY*I_ERI_Fx2y_S_F3z_S_c;
  Double I_ERI_Fxyz_Py_F3z_S_c = I_ERI_Gx2yz_S_F3z_S_c+ABY*I_ERI_Fxyz_S_F3z_S_c;
  Double I_ERI_Fx2z_Py_F3z_S_c = I_ERI_Gxy2z_S_F3z_S_c+ABY*I_ERI_Fx2z_S_F3z_S_c;
  Double I_ERI_F3y_Py_F3z_S_c = I_ERI_G4y_S_F3z_S_c+ABY*I_ERI_F3y_S_F3z_S_c;
  Double I_ERI_F2yz_Py_F3z_S_c = I_ERI_G3yz_S_F3z_S_c+ABY*I_ERI_F2yz_S_F3z_S_c;
  Double I_ERI_Fy2z_Py_F3z_S_c = I_ERI_G2y2z_S_F3z_S_c+ABY*I_ERI_Fy2z_S_F3z_S_c;
  Double I_ERI_F3z_Py_F3z_S_c = I_ERI_Gy3z_S_F3z_S_c+ABY*I_ERI_F3z_S_F3z_S_c;
  Double I_ERI_F2xz_Pz_F3z_S_c = I_ERI_G2x2z_S_F3z_S_c+ABZ*I_ERI_F2xz_S_F3z_S_c;
  Double I_ERI_Fxyz_Pz_F3z_S_c = I_ERI_Gxy2z_S_F3z_S_c+ABZ*I_ERI_Fxyz_S_F3z_S_c;
  Double I_ERI_Fx2z_Pz_F3z_S_c = I_ERI_Gx3z_S_F3z_S_c+ABZ*I_ERI_Fx2z_S_F3z_S_c;
  Double I_ERI_F2yz_Pz_F3z_S_c = I_ERI_G2y2z_S_F3z_S_c+ABZ*I_ERI_F2yz_S_F3z_S_c;
  Double I_ERI_Fy2z_Pz_F3z_S_c = I_ERI_Gy3z_S_F3z_S_c+ABZ*I_ERI_Fy2z_S_F3z_S_c;
  Double I_ERI_F3z_Pz_F3z_S_c = I_ERI_G4z_S_F3z_S_c+ABZ*I_ERI_F3z_S_F3z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_F_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_S_c
   * RHS shell quartet name: SQ_ERI_D_P_F_S_c
   ************************************************************/
  Double I_ERI_D2x_D2x_F3x_S_c = I_ERI_F3x_Px_F3x_S_c+ABX*I_ERI_D2x_Px_F3x_S_c;
  Double I_ERI_Dxy_D2x_F3x_S_c = I_ERI_F2xy_Px_F3x_S_c+ABX*I_ERI_Dxy_Px_F3x_S_c;
  Double I_ERI_Dxz_D2x_F3x_S_c = I_ERI_F2xz_Px_F3x_S_c+ABX*I_ERI_Dxz_Px_F3x_S_c;
  Double I_ERI_D2y_D2x_F3x_S_c = I_ERI_Fx2y_Px_F3x_S_c+ABX*I_ERI_D2y_Px_F3x_S_c;
  Double I_ERI_Dyz_D2x_F3x_S_c = I_ERI_Fxyz_Px_F3x_S_c+ABX*I_ERI_Dyz_Px_F3x_S_c;
  Double I_ERI_D2z_D2x_F3x_S_c = I_ERI_Fx2z_Px_F3x_S_c+ABX*I_ERI_D2z_Px_F3x_S_c;
  Double I_ERI_D2x_Dxy_F3x_S_c = I_ERI_F2xy_Px_F3x_S_c+ABY*I_ERI_D2x_Px_F3x_S_c;
  Double I_ERI_Dxy_Dxy_F3x_S_c = I_ERI_Fx2y_Px_F3x_S_c+ABY*I_ERI_Dxy_Px_F3x_S_c;
  Double I_ERI_Dxz_Dxy_F3x_S_c = I_ERI_Fxyz_Px_F3x_S_c+ABY*I_ERI_Dxz_Px_F3x_S_c;
  Double I_ERI_D2y_Dxy_F3x_S_c = I_ERI_F3y_Px_F3x_S_c+ABY*I_ERI_D2y_Px_F3x_S_c;
  Double I_ERI_Dyz_Dxy_F3x_S_c = I_ERI_F2yz_Px_F3x_S_c+ABY*I_ERI_Dyz_Px_F3x_S_c;
  Double I_ERI_D2z_Dxy_F3x_S_c = I_ERI_Fy2z_Px_F3x_S_c+ABY*I_ERI_D2z_Px_F3x_S_c;
  Double I_ERI_D2x_Dxz_F3x_S_c = I_ERI_F2xz_Px_F3x_S_c+ABZ*I_ERI_D2x_Px_F3x_S_c;
  Double I_ERI_Dxy_Dxz_F3x_S_c = I_ERI_Fxyz_Px_F3x_S_c+ABZ*I_ERI_Dxy_Px_F3x_S_c;
  Double I_ERI_Dxz_Dxz_F3x_S_c = I_ERI_Fx2z_Px_F3x_S_c+ABZ*I_ERI_Dxz_Px_F3x_S_c;
  Double I_ERI_D2y_Dxz_F3x_S_c = I_ERI_F2yz_Px_F3x_S_c+ABZ*I_ERI_D2y_Px_F3x_S_c;
  Double I_ERI_Dyz_Dxz_F3x_S_c = I_ERI_Fy2z_Px_F3x_S_c+ABZ*I_ERI_Dyz_Px_F3x_S_c;
  Double I_ERI_D2z_Dxz_F3x_S_c = I_ERI_F3z_Px_F3x_S_c+ABZ*I_ERI_D2z_Px_F3x_S_c;
  Double I_ERI_D2x_D2y_F3x_S_c = I_ERI_F2xy_Py_F3x_S_c+ABY*I_ERI_D2x_Py_F3x_S_c;
  Double I_ERI_Dxy_D2y_F3x_S_c = I_ERI_Fx2y_Py_F3x_S_c+ABY*I_ERI_Dxy_Py_F3x_S_c;
  Double I_ERI_Dxz_D2y_F3x_S_c = I_ERI_Fxyz_Py_F3x_S_c+ABY*I_ERI_Dxz_Py_F3x_S_c;
  Double I_ERI_D2y_D2y_F3x_S_c = I_ERI_F3y_Py_F3x_S_c+ABY*I_ERI_D2y_Py_F3x_S_c;
  Double I_ERI_Dyz_D2y_F3x_S_c = I_ERI_F2yz_Py_F3x_S_c+ABY*I_ERI_Dyz_Py_F3x_S_c;
  Double I_ERI_D2z_D2y_F3x_S_c = I_ERI_Fy2z_Py_F3x_S_c+ABY*I_ERI_D2z_Py_F3x_S_c;
  Double I_ERI_D2x_Dyz_F3x_S_c = I_ERI_F2xz_Py_F3x_S_c+ABZ*I_ERI_D2x_Py_F3x_S_c;
  Double I_ERI_Dxy_Dyz_F3x_S_c = I_ERI_Fxyz_Py_F3x_S_c+ABZ*I_ERI_Dxy_Py_F3x_S_c;
  Double I_ERI_Dxz_Dyz_F3x_S_c = I_ERI_Fx2z_Py_F3x_S_c+ABZ*I_ERI_Dxz_Py_F3x_S_c;
  Double I_ERI_D2y_Dyz_F3x_S_c = I_ERI_F2yz_Py_F3x_S_c+ABZ*I_ERI_D2y_Py_F3x_S_c;
  Double I_ERI_Dyz_Dyz_F3x_S_c = I_ERI_Fy2z_Py_F3x_S_c+ABZ*I_ERI_Dyz_Py_F3x_S_c;
  Double I_ERI_D2z_Dyz_F3x_S_c = I_ERI_F3z_Py_F3x_S_c+ABZ*I_ERI_D2z_Py_F3x_S_c;
  Double I_ERI_D2x_D2z_F3x_S_c = I_ERI_F2xz_Pz_F3x_S_c+ABZ*I_ERI_D2x_Pz_F3x_S_c;
  Double I_ERI_Dxy_D2z_F3x_S_c = I_ERI_Fxyz_Pz_F3x_S_c+ABZ*I_ERI_Dxy_Pz_F3x_S_c;
  Double I_ERI_Dxz_D2z_F3x_S_c = I_ERI_Fx2z_Pz_F3x_S_c+ABZ*I_ERI_Dxz_Pz_F3x_S_c;
  Double I_ERI_D2y_D2z_F3x_S_c = I_ERI_F2yz_Pz_F3x_S_c+ABZ*I_ERI_D2y_Pz_F3x_S_c;
  Double I_ERI_Dyz_D2z_F3x_S_c = I_ERI_Fy2z_Pz_F3x_S_c+ABZ*I_ERI_Dyz_Pz_F3x_S_c;
  Double I_ERI_D2z_D2z_F3x_S_c = I_ERI_F3z_Pz_F3x_S_c+ABZ*I_ERI_D2z_Pz_F3x_S_c;
  Double I_ERI_D2x_D2x_F2xy_S_c = I_ERI_F3x_Px_F2xy_S_c+ABX*I_ERI_D2x_Px_F2xy_S_c;
  Double I_ERI_Dxy_D2x_F2xy_S_c = I_ERI_F2xy_Px_F2xy_S_c+ABX*I_ERI_Dxy_Px_F2xy_S_c;
  Double I_ERI_Dxz_D2x_F2xy_S_c = I_ERI_F2xz_Px_F2xy_S_c+ABX*I_ERI_Dxz_Px_F2xy_S_c;
  Double I_ERI_D2y_D2x_F2xy_S_c = I_ERI_Fx2y_Px_F2xy_S_c+ABX*I_ERI_D2y_Px_F2xy_S_c;
  Double I_ERI_Dyz_D2x_F2xy_S_c = I_ERI_Fxyz_Px_F2xy_S_c+ABX*I_ERI_Dyz_Px_F2xy_S_c;
  Double I_ERI_D2z_D2x_F2xy_S_c = I_ERI_Fx2z_Px_F2xy_S_c+ABX*I_ERI_D2z_Px_F2xy_S_c;
  Double I_ERI_D2x_Dxy_F2xy_S_c = I_ERI_F2xy_Px_F2xy_S_c+ABY*I_ERI_D2x_Px_F2xy_S_c;
  Double I_ERI_Dxy_Dxy_F2xy_S_c = I_ERI_Fx2y_Px_F2xy_S_c+ABY*I_ERI_Dxy_Px_F2xy_S_c;
  Double I_ERI_Dxz_Dxy_F2xy_S_c = I_ERI_Fxyz_Px_F2xy_S_c+ABY*I_ERI_Dxz_Px_F2xy_S_c;
  Double I_ERI_D2y_Dxy_F2xy_S_c = I_ERI_F3y_Px_F2xy_S_c+ABY*I_ERI_D2y_Px_F2xy_S_c;
  Double I_ERI_Dyz_Dxy_F2xy_S_c = I_ERI_F2yz_Px_F2xy_S_c+ABY*I_ERI_Dyz_Px_F2xy_S_c;
  Double I_ERI_D2z_Dxy_F2xy_S_c = I_ERI_Fy2z_Px_F2xy_S_c+ABY*I_ERI_D2z_Px_F2xy_S_c;
  Double I_ERI_D2x_Dxz_F2xy_S_c = I_ERI_F2xz_Px_F2xy_S_c+ABZ*I_ERI_D2x_Px_F2xy_S_c;
  Double I_ERI_Dxy_Dxz_F2xy_S_c = I_ERI_Fxyz_Px_F2xy_S_c+ABZ*I_ERI_Dxy_Px_F2xy_S_c;
  Double I_ERI_Dxz_Dxz_F2xy_S_c = I_ERI_Fx2z_Px_F2xy_S_c+ABZ*I_ERI_Dxz_Px_F2xy_S_c;
  Double I_ERI_D2y_Dxz_F2xy_S_c = I_ERI_F2yz_Px_F2xy_S_c+ABZ*I_ERI_D2y_Px_F2xy_S_c;
  Double I_ERI_Dyz_Dxz_F2xy_S_c = I_ERI_Fy2z_Px_F2xy_S_c+ABZ*I_ERI_Dyz_Px_F2xy_S_c;
  Double I_ERI_D2z_Dxz_F2xy_S_c = I_ERI_F3z_Px_F2xy_S_c+ABZ*I_ERI_D2z_Px_F2xy_S_c;
  Double I_ERI_D2x_D2y_F2xy_S_c = I_ERI_F2xy_Py_F2xy_S_c+ABY*I_ERI_D2x_Py_F2xy_S_c;
  Double I_ERI_Dxy_D2y_F2xy_S_c = I_ERI_Fx2y_Py_F2xy_S_c+ABY*I_ERI_Dxy_Py_F2xy_S_c;
  Double I_ERI_Dxz_D2y_F2xy_S_c = I_ERI_Fxyz_Py_F2xy_S_c+ABY*I_ERI_Dxz_Py_F2xy_S_c;
  Double I_ERI_D2y_D2y_F2xy_S_c = I_ERI_F3y_Py_F2xy_S_c+ABY*I_ERI_D2y_Py_F2xy_S_c;
  Double I_ERI_Dyz_D2y_F2xy_S_c = I_ERI_F2yz_Py_F2xy_S_c+ABY*I_ERI_Dyz_Py_F2xy_S_c;
  Double I_ERI_D2z_D2y_F2xy_S_c = I_ERI_Fy2z_Py_F2xy_S_c+ABY*I_ERI_D2z_Py_F2xy_S_c;
  Double I_ERI_D2x_Dyz_F2xy_S_c = I_ERI_F2xz_Py_F2xy_S_c+ABZ*I_ERI_D2x_Py_F2xy_S_c;
  Double I_ERI_Dxy_Dyz_F2xy_S_c = I_ERI_Fxyz_Py_F2xy_S_c+ABZ*I_ERI_Dxy_Py_F2xy_S_c;
  Double I_ERI_Dxz_Dyz_F2xy_S_c = I_ERI_Fx2z_Py_F2xy_S_c+ABZ*I_ERI_Dxz_Py_F2xy_S_c;
  Double I_ERI_D2y_Dyz_F2xy_S_c = I_ERI_F2yz_Py_F2xy_S_c+ABZ*I_ERI_D2y_Py_F2xy_S_c;
  Double I_ERI_Dyz_Dyz_F2xy_S_c = I_ERI_Fy2z_Py_F2xy_S_c+ABZ*I_ERI_Dyz_Py_F2xy_S_c;
  Double I_ERI_D2z_Dyz_F2xy_S_c = I_ERI_F3z_Py_F2xy_S_c+ABZ*I_ERI_D2z_Py_F2xy_S_c;
  Double I_ERI_D2x_D2z_F2xy_S_c = I_ERI_F2xz_Pz_F2xy_S_c+ABZ*I_ERI_D2x_Pz_F2xy_S_c;
  Double I_ERI_Dxy_D2z_F2xy_S_c = I_ERI_Fxyz_Pz_F2xy_S_c+ABZ*I_ERI_Dxy_Pz_F2xy_S_c;
  Double I_ERI_Dxz_D2z_F2xy_S_c = I_ERI_Fx2z_Pz_F2xy_S_c+ABZ*I_ERI_Dxz_Pz_F2xy_S_c;
  Double I_ERI_D2y_D2z_F2xy_S_c = I_ERI_F2yz_Pz_F2xy_S_c+ABZ*I_ERI_D2y_Pz_F2xy_S_c;
  Double I_ERI_Dyz_D2z_F2xy_S_c = I_ERI_Fy2z_Pz_F2xy_S_c+ABZ*I_ERI_Dyz_Pz_F2xy_S_c;
  Double I_ERI_D2z_D2z_F2xy_S_c = I_ERI_F3z_Pz_F2xy_S_c+ABZ*I_ERI_D2z_Pz_F2xy_S_c;
  Double I_ERI_D2x_D2x_F2xz_S_c = I_ERI_F3x_Px_F2xz_S_c+ABX*I_ERI_D2x_Px_F2xz_S_c;
  Double I_ERI_Dxy_D2x_F2xz_S_c = I_ERI_F2xy_Px_F2xz_S_c+ABX*I_ERI_Dxy_Px_F2xz_S_c;
  Double I_ERI_Dxz_D2x_F2xz_S_c = I_ERI_F2xz_Px_F2xz_S_c+ABX*I_ERI_Dxz_Px_F2xz_S_c;
  Double I_ERI_D2y_D2x_F2xz_S_c = I_ERI_Fx2y_Px_F2xz_S_c+ABX*I_ERI_D2y_Px_F2xz_S_c;
  Double I_ERI_Dyz_D2x_F2xz_S_c = I_ERI_Fxyz_Px_F2xz_S_c+ABX*I_ERI_Dyz_Px_F2xz_S_c;
  Double I_ERI_D2z_D2x_F2xz_S_c = I_ERI_Fx2z_Px_F2xz_S_c+ABX*I_ERI_D2z_Px_F2xz_S_c;
  Double I_ERI_D2x_Dxy_F2xz_S_c = I_ERI_F2xy_Px_F2xz_S_c+ABY*I_ERI_D2x_Px_F2xz_S_c;
  Double I_ERI_Dxy_Dxy_F2xz_S_c = I_ERI_Fx2y_Px_F2xz_S_c+ABY*I_ERI_Dxy_Px_F2xz_S_c;
  Double I_ERI_Dxz_Dxy_F2xz_S_c = I_ERI_Fxyz_Px_F2xz_S_c+ABY*I_ERI_Dxz_Px_F2xz_S_c;
  Double I_ERI_D2y_Dxy_F2xz_S_c = I_ERI_F3y_Px_F2xz_S_c+ABY*I_ERI_D2y_Px_F2xz_S_c;
  Double I_ERI_Dyz_Dxy_F2xz_S_c = I_ERI_F2yz_Px_F2xz_S_c+ABY*I_ERI_Dyz_Px_F2xz_S_c;
  Double I_ERI_D2z_Dxy_F2xz_S_c = I_ERI_Fy2z_Px_F2xz_S_c+ABY*I_ERI_D2z_Px_F2xz_S_c;
  Double I_ERI_D2x_Dxz_F2xz_S_c = I_ERI_F2xz_Px_F2xz_S_c+ABZ*I_ERI_D2x_Px_F2xz_S_c;
  Double I_ERI_Dxy_Dxz_F2xz_S_c = I_ERI_Fxyz_Px_F2xz_S_c+ABZ*I_ERI_Dxy_Px_F2xz_S_c;
  Double I_ERI_Dxz_Dxz_F2xz_S_c = I_ERI_Fx2z_Px_F2xz_S_c+ABZ*I_ERI_Dxz_Px_F2xz_S_c;
  Double I_ERI_D2y_Dxz_F2xz_S_c = I_ERI_F2yz_Px_F2xz_S_c+ABZ*I_ERI_D2y_Px_F2xz_S_c;
  Double I_ERI_Dyz_Dxz_F2xz_S_c = I_ERI_Fy2z_Px_F2xz_S_c+ABZ*I_ERI_Dyz_Px_F2xz_S_c;
  Double I_ERI_D2z_Dxz_F2xz_S_c = I_ERI_F3z_Px_F2xz_S_c+ABZ*I_ERI_D2z_Px_F2xz_S_c;
  Double I_ERI_D2x_D2y_F2xz_S_c = I_ERI_F2xy_Py_F2xz_S_c+ABY*I_ERI_D2x_Py_F2xz_S_c;
  Double I_ERI_Dxy_D2y_F2xz_S_c = I_ERI_Fx2y_Py_F2xz_S_c+ABY*I_ERI_Dxy_Py_F2xz_S_c;
  Double I_ERI_Dxz_D2y_F2xz_S_c = I_ERI_Fxyz_Py_F2xz_S_c+ABY*I_ERI_Dxz_Py_F2xz_S_c;
  Double I_ERI_D2y_D2y_F2xz_S_c = I_ERI_F3y_Py_F2xz_S_c+ABY*I_ERI_D2y_Py_F2xz_S_c;
  Double I_ERI_Dyz_D2y_F2xz_S_c = I_ERI_F2yz_Py_F2xz_S_c+ABY*I_ERI_Dyz_Py_F2xz_S_c;
  Double I_ERI_D2z_D2y_F2xz_S_c = I_ERI_Fy2z_Py_F2xz_S_c+ABY*I_ERI_D2z_Py_F2xz_S_c;
  Double I_ERI_D2x_Dyz_F2xz_S_c = I_ERI_F2xz_Py_F2xz_S_c+ABZ*I_ERI_D2x_Py_F2xz_S_c;
  Double I_ERI_Dxy_Dyz_F2xz_S_c = I_ERI_Fxyz_Py_F2xz_S_c+ABZ*I_ERI_Dxy_Py_F2xz_S_c;
  Double I_ERI_Dxz_Dyz_F2xz_S_c = I_ERI_Fx2z_Py_F2xz_S_c+ABZ*I_ERI_Dxz_Py_F2xz_S_c;
  Double I_ERI_D2y_Dyz_F2xz_S_c = I_ERI_F2yz_Py_F2xz_S_c+ABZ*I_ERI_D2y_Py_F2xz_S_c;
  Double I_ERI_Dyz_Dyz_F2xz_S_c = I_ERI_Fy2z_Py_F2xz_S_c+ABZ*I_ERI_Dyz_Py_F2xz_S_c;
  Double I_ERI_D2z_Dyz_F2xz_S_c = I_ERI_F3z_Py_F2xz_S_c+ABZ*I_ERI_D2z_Py_F2xz_S_c;
  Double I_ERI_D2x_D2z_F2xz_S_c = I_ERI_F2xz_Pz_F2xz_S_c+ABZ*I_ERI_D2x_Pz_F2xz_S_c;
  Double I_ERI_Dxy_D2z_F2xz_S_c = I_ERI_Fxyz_Pz_F2xz_S_c+ABZ*I_ERI_Dxy_Pz_F2xz_S_c;
  Double I_ERI_Dxz_D2z_F2xz_S_c = I_ERI_Fx2z_Pz_F2xz_S_c+ABZ*I_ERI_Dxz_Pz_F2xz_S_c;
  Double I_ERI_D2y_D2z_F2xz_S_c = I_ERI_F2yz_Pz_F2xz_S_c+ABZ*I_ERI_D2y_Pz_F2xz_S_c;
  Double I_ERI_Dyz_D2z_F2xz_S_c = I_ERI_Fy2z_Pz_F2xz_S_c+ABZ*I_ERI_Dyz_Pz_F2xz_S_c;
  Double I_ERI_D2z_D2z_F2xz_S_c = I_ERI_F3z_Pz_F2xz_S_c+ABZ*I_ERI_D2z_Pz_F2xz_S_c;
  Double I_ERI_D2x_D2x_Fx2y_S_c = I_ERI_F3x_Px_Fx2y_S_c+ABX*I_ERI_D2x_Px_Fx2y_S_c;
  Double I_ERI_Dxy_D2x_Fx2y_S_c = I_ERI_F2xy_Px_Fx2y_S_c+ABX*I_ERI_Dxy_Px_Fx2y_S_c;
  Double I_ERI_Dxz_D2x_Fx2y_S_c = I_ERI_F2xz_Px_Fx2y_S_c+ABX*I_ERI_Dxz_Px_Fx2y_S_c;
  Double I_ERI_D2y_D2x_Fx2y_S_c = I_ERI_Fx2y_Px_Fx2y_S_c+ABX*I_ERI_D2y_Px_Fx2y_S_c;
  Double I_ERI_Dyz_D2x_Fx2y_S_c = I_ERI_Fxyz_Px_Fx2y_S_c+ABX*I_ERI_Dyz_Px_Fx2y_S_c;
  Double I_ERI_D2z_D2x_Fx2y_S_c = I_ERI_Fx2z_Px_Fx2y_S_c+ABX*I_ERI_D2z_Px_Fx2y_S_c;
  Double I_ERI_D2x_Dxy_Fx2y_S_c = I_ERI_F2xy_Px_Fx2y_S_c+ABY*I_ERI_D2x_Px_Fx2y_S_c;
  Double I_ERI_Dxy_Dxy_Fx2y_S_c = I_ERI_Fx2y_Px_Fx2y_S_c+ABY*I_ERI_Dxy_Px_Fx2y_S_c;
  Double I_ERI_Dxz_Dxy_Fx2y_S_c = I_ERI_Fxyz_Px_Fx2y_S_c+ABY*I_ERI_Dxz_Px_Fx2y_S_c;
  Double I_ERI_D2y_Dxy_Fx2y_S_c = I_ERI_F3y_Px_Fx2y_S_c+ABY*I_ERI_D2y_Px_Fx2y_S_c;
  Double I_ERI_Dyz_Dxy_Fx2y_S_c = I_ERI_F2yz_Px_Fx2y_S_c+ABY*I_ERI_Dyz_Px_Fx2y_S_c;
  Double I_ERI_D2z_Dxy_Fx2y_S_c = I_ERI_Fy2z_Px_Fx2y_S_c+ABY*I_ERI_D2z_Px_Fx2y_S_c;
  Double I_ERI_D2x_Dxz_Fx2y_S_c = I_ERI_F2xz_Px_Fx2y_S_c+ABZ*I_ERI_D2x_Px_Fx2y_S_c;
  Double I_ERI_Dxy_Dxz_Fx2y_S_c = I_ERI_Fxyz_Px_Fx2y_S_c+ABZ*I_ERI_Dxy_Px_Fx2y_S_c;
  Double I_ERI_Dxz_Dxz_Fx2y_S_c = I_ERI_Fx2z_Px_Fx2y_S_c+ABZ*I_ERI_Dxz_Px_Fx2y_S_c;
  Double I_ERI_D2y_Dxz_Fx2y_S_c = I_ERI_F2yz_Px_Fx2y_S_c+ABZ*I_ERI_D2y_Px_Fx2y_S_c;
  Double I_ERI_Dyz_Dxz_Fx2y_S_c = I_ERI_Fy2z_Px_Fx2y_S_c+ABZ*I_ERI_Dyz_Px_Fx2y_S_c;
  Double I_ERI_D2z_Dxz_Fx2y_S_c = I_ERI_F3z_Px_Fx2y_S_c+ABZ*I_ERI_D2z_Px_Fx2y_S_c;
  Double I_ERI_D2x_D2y_Fx2y_S_c = I_ERI_F2xy_Py_Fx2y_S_c+ABY*I_ERI_D2x_Py_Fx2y_S_c;
  Double I_ERI_Dxy_D2y_Fx2y_S_c = I_ERI_Fx2y_Py_Fx2y_S_c+ABY*I_ERI_Dxy_Py_Fx2y_S_c;
  Double I_ERI_Dxz_D2y_Fx2y_S_c = I_ERI_Fxyz_Py_Fx2y_S_c+ABY*I_ERI_Dxz_Py_Fx2y_S_c;
  Double I_ERI_D2y_D2y_Fx2y_S_c = I_ERI_F3y_Py_Fx2y_S_c+ABY*I_ERI_D2y_Py_Fx2y_S_c;
  Double I_ERI_Dyz_D2y_Fx2y_S_c = I_ERI_F2yz_Py_Fx2y_S_c+ABY*I_ERI_Dyz_Py_Fx2y_S_c;
  Double I_ERI_D2z_D2y_Fx2y_S_c = I_ERI_Fy2z_Py_Fx2y_S_c+ABY*I_ERI_D2z_Py_Fx2y_S_c;
  Double I_ERI_D2x_Dyz_Fx2y_S_c = I_ERI_F2xz_Py_Fx2y_S_c+ABZ*I_ERI_D2x_Py_Fx2y_S_c;
  Double I_ERI_Dxy_Dyz_Fx2y_S_c = I_ERI_Fxyz_Py_Fx2y_S_c+ABZ*I_ERI_Dxy_Py_Fx2y_S_c;
  Double I_ERI_Dxz_Dyz_Fx2y_S_c = I_ERI_Fx2z_Py_Fx2y_S_c+ABZ*I_ERI_Dxz_Py_Fx2y_S_c;
  Double I_ERI_D2y_Dyz_Fx2y_S_c = I_ERI_F2yz_Py_Fx2y_S_c+ABZ*I_ERI_D2y_Py_Fx2y_S_c;
  Double I_ERI_Dyz_Dyz_Fx2y_S_c = I_ERI_Fy2z_Py_Fx2y_S_c+ABZ*I_ERI_Dyz_Py_Fx2y_S_c;
  Double I_ERI_D2z_Dyz_Fx2y_S_c = I_ERI_F3z_Py_Fx2y_S_c+ABZ*I_ERI_D2z_Py_Fx2y_S_c;
  Double I_ERI_D2x_D2z_Fx2y_S_c = I_ERI_F2xz_Pz_Fx2y_S_c+ABZ*I_ERI_D2x_Pz_Fx2y_S_c;
  Double I_ERI_Dxy_D2z_Fx2y_S_c = I_ERI_Fxyz_Pz_Fx2y_S_c+ABZ*I_ERI_Dxy_Pz_Fx2y_S_c;
  Double I_ERI_Dxz_D2z_Fx2y_S_c = I_ERI_Fx2z_Pz_Fx2y_S_c+ABZ*I_ERI_Dxz_Pz_Fx2y_S_c;
  Double I_ERI_D2y_D2z_Fx2y_S_c = I_ERI_F2yz_Pz_Fx2y_S_c+ABZ*I_ERI_D2y_Pz_Fx2y_S_c;
  Double I_ERI_Dyz_D2z_Fx2y_S_c = I_ERI_Fy2z_Pz_Fx2y_S_c+ABZ*I_ERI_Dyz_Pz_Fx2y_S_c;
  Double I_ERI_D2z_D2z_Fx2y_S_c = I_ERI_F3z_Pz_Fx2y_S_c+ABZ*I_ERI_D2z_Pz_Fx2y_S_c;
  Double I_ERI_D2x_D2x_Fxyz_S_c = I_ERI_F3x_Px_Fxyz_S_c+ABX*I_ERI_D2x_Px_Fxyz_S_c;
  Double I_ERI_Dxy_D2x_Fxyz_S_c = I_ERI_F2xy_Px_Fxyz_S_c+ABX*I_ERI_Dxy_Px_Fxyz_S_c;
  Double I_ERI_Dxz_D2x_Fxyz_S_c = I_ERI_F2xz_Px_Fxyz_S_c+ABX*I_ERI_Dxz_Px_Fxyz_S_c;
  Double I_ERI_D2y_D2x_Fxyz_S_c = I_ERI_Fx2y_Px_Fxyz_S_c+ABX*I_ERI_D2y_Px_Fxyz_S_c;
  Double I_ERI_Dyz_D2x_Fxyz_S_c = I_ERI_Fxyz_Px_Fxyz_S_c+ABX*I_ERI_Dyz_Px_Fxyz_S_c;
  Double I_ERI_D2z_D2x_Fxyz_S_c = I_ERI_Fx2z_Px_Fxyz_S_c+ABX*I_ERI_D2z_Px_Fxyz_S_c;
  Double I_ERI_D2x_Dxy_Fxyz_S_c = I_ERI_F2xy_Px_Fxyz_S_c+ABY*I_ERI_D2x_Px_Fxyz_S_c;
  Double I_ERI_Dxy_Dxy_Fxyz_S_c = I_ERI_Fx2y_Px_Fxyz_S_c+ABY*I_ERI_Dxy_Px_Fxyz_S_c;
  Double I_ERI_Dxz_Dxy_Fxyz_S_c = I_ERI_Fxyz_Px_Fxyz_S_c+ABY*I_ERI_Dxz_Px_Fxyz_S_c;
  Double I_ERI_D2y_Dxy_Fxyz_S_c = I_ERI_F3y_Px_Fxyz_S_c+ABY*I_ERI_D2y_Px_Fxyz_S_c;
  Double I_ERI_Dyz_Dxy_Fxyz_S_c = I_ERI_F2yz_Px_Fxyz_S_c+ABY*I_ERI_Dyz_Px_Fxyz_S_c;
  Double I_ERI_D2z_Dxy_Fxyz_S_c = I_ERI_Fy2z_Px_Fxyz_S_c+ABY*I_ERI_D2z_Px_Fxyz_S_c;
  Double I_ERI_D2x_Dxz_Fxyz_S_c = I_ERI_F2xz_Px_Fxyz_S_c+ABZ*I_ERI_D2x_Px_Fxyz_S_c;
  Double I_ERI_Dxy_Dxz_Fxyz_S_c = I_ERI_Fxyz_Px_Fxyz_S_c+ABZ*I_ERI_Dxy_Px_Fxyz_S_c;
  Double I_ERI_Dxz_Dxz_Fxyz_S_c = I_ERI_Fx2z_Px_Fxyz_S_c+ABZ*I_ERI_Dxz_Px_Fxyz_S_c;
  Double I_ERI_D2y_Dxz_Fxyz_S_c = I_ERI_F2yz_Px_Fxyz_S_c+ABZ*I_ERI_D2y_Px_Fxyz_S_c;
  Double I_ERI_Dyz_Dxz_Fxyz_S_c = I_ERI_Fy2z_Px_Fxyz_S_c+ABZ*I_ERI_Dyz_Px_Fxyz_S_c;
  Double I_ERI_D2z_Dxz_Fxyz_S_c = I_ERI_F3z_Px_Fxyz_S_c+ABZ*I_ERI_D2z_Px_Fxyz_S_c;
  Double I_ERI_D2x_D2y_Fxyz_S_c = I_ERI_F2xy_Py_Fxyz_S_c+ABY*I_ERI_D2x_Py_Fxyz_S_c;
  Double I_ERI_Dxy_D2y_Fxyz_S_c = I_ERI_Fx2y_Py_Fxyz_S_c+ABY*I_ERI_Dxy_Py_Fxyz_S_c;
  Double I_ERI_Dxz_D2y_Fxyz_S_c = I_ERI_Fxyz_Py_Fxyz_S_c+ABY*I_ERI_Dxz_Py_Fxyz_S_c;
  Double I_ERI_D2y_D2y_Fxyz_S_c = I_ERI_F3y_Py_Fxyz_S_c+ABY*I_ERI_D2y_Py_Fxyz_S_c;
  Double I_ERI_Dyz_D2y_Fxyz_S_c = I_ERI_F2yz_Py_Fxyz_S_c+ABY*I_ERI_Dyz_Py_Fxyz_S_c;
  Double I_ERI_D2z_D2y_Fxyz_S_c = I_ERI_Fy2z_Py_Fxyz_S_c+ABY*I_ERI_D2z_Py_Fxyz_S_c;
  Double I_ERI_D2x_Dyz_Fxyz_S_c = I_ERI_F2xz_Py_Fxyz_S_c+ABZ*I_ERI_D2x_Py_Fxyz_S_c;
  Double I_ERI_Dxy_Dyz_Fxyz_S_c = I_ERI_Fxyz_Py_Fxyz_S_c+ABZ*I_ERI_Dxy_Py_Fxyz_S_c;
  Double I_ERI_Dxz_Dyz_Fxyz_S_c = I_ERI_Fx2z_Py_Fxyz_S_c+ABZ*I_ERI_Dxz_Py_Fxyz_S_c;
  Double I_ERI_D2y_Dyz_Fxyz_S_c = I_ERI_F2yz_Py_Fxyz_S_c+ABZ*I_ERI_D2y_Py_Fxyz_S_c;
  Double I_ERI_Dyz_Dyz_Fxyz_S_c = I_ERI_Fy2z_Py_Fxyz_S_c+ABZ*I_ERI_Dyz_Py_Fxyz_S_c;
  Double I_ERI_D2z_Dyz_Fxyz_S_c = I_ERI_F3z_Py_Fxyz_S_c+ABZ*I_ERI_D2z_Py_Fxyz_S_c;
  Double I_ERI_D2x_D2z_Fxyz_S_c = I_ERI_F2xz_Pz_Fxyz_S_c+ABZ*I_ERI_D2x_Pz_Fxyz_S_c;
  Double I_ERI_Dxy_D2z_Fxyz_S_c = I_ERI_Fxyz_Pz_Fxyz_S_c+ABZ*I_ERI_Dxy_Pz_Fxyz_S_c;
  Double I_ERI_Dxz_D2z_Fxyz_S_c = I_ERI_Fx2z_Pz_Fxyz_S_c+ABZ*I_ERI_Dxz_Pz_Fxyz_S_c;
  Double I_ERI_D2y_D2z_Fxyz_S_c = I_ERI_F2yz_Pz_Fxyz_S_c+ABZ*I_ERI_D2y_Pz_Fxyz_S_c;
  Double I_ERI_Dyz_D2z_Fxyz_S_c = I_ERI_Fy2z_Pz_Fxyz_S_c+ABZ*I_ERI_Dyz_Pz_Fxyz_S_c;
  Double I_ERI_D2z_D2z_Fxyz_S_c = I_ERI_F3z_Pz_Fxyz_S_c+ABZ*I_ERI_D2z_Pz_Fxyz_S_c;
  Double I_ERI_D2x_D2x_Fx2z_S_c = I_ERI_F3x_Px_Fx2z_S_c+ABX*I_ERI_D2x_Px_Fx2z_S_c;
  Double I_ERI_Dxy_D2x_Fx2z_S_c = I_ERI_F2xy_Px_Fx2z_S_c+ABX*I_ERI_Dxy_Px_Fx2z_S_c;
  Double I_ERI_Dxz_D2x_Fx2z_S_c = I_ERI_F2xz_Px_Fx2z_S_c+ABX*I_ERI_Dxz_Px_Fx2z_S_c;
  Double I_ERI_D2y_D2x_Fx2z_S_c = I_ERI_Fx2y_Px_Fx2z_S_c+ABX*I_ERI_D2y_Px_Fx2z_S_c;
  Double I_ERI_Dyz_D2x_Fx2z_S_c = I_ERI_Fxyz_Px_Fx2z_S_c+ABX*I_ERI_Dyz_Px_Fx2z_S_c;
  Double I_ERI_D2z_D2x_Fx2z_S_c = I_ERI_Fx2z_Px_Fx2z_S_c+ABX*I_ERI_D2z_Px_Fx2z_S_c;
  Double I_ERI_D2x_Dxy_Fx2z_S_c = I_ERI_F2xy_Px_Fx2z_S_c+ABY*I_ERI_D2x_Px_Fx2z_S_c;
  Double I_ERI_Dxy_Dxy_Fx2z_S_c = I_ERI_Fx2y_Px_Fx2z_S_c+ABY*I_ERI_Dxy_Px_Fx2z_S_c;
  Double I_ERI_Dxz_Dxy_Fx2z_S_c = I_ERI_Fxyz_Px_Fx2z_S_c+ABY*I_ERI_Dxz_Px_Fx2z_S_c;
  Double I_ERI_D2y_Dxy_Fx2z_S_c = I_ERI_F3y_Px_Fx2z_S_c+ABY*I_ERI_D2y_Px_Fx2z_S_c;
  Double I_ERI_Dyz_Dxy_Fx2z_S_c = I_ERI_F2yz_Px_Fx2z_S_c+ABY*I_ERI_Dyz_Px_Fx2z_S_c;
  Double I_ERI_D2z_Dxy_Fx2z_S_c = I_ERI_Fy2z_Px_Fx2z_S_c+ABY*I_ERI_D2z_Px_Fx2z_S_c;
  Double I_ERI_D2x_Dxz_Fx2z_S_c = I_ERI_F2xz_Px_Fx2z_S_c+ABZ*I_ERI_D2x_Px_Fx2z_S_c;
  Double I_ERI_Dxy_Dxz_Fx2z_S_c = I_ERI_Fxyz_Px_Fx2z_S_c+ABZ*I_ERI_Dxy_Px_Fx2z_S_c;
  Double I_ERI_Dxz_Dxz_Fx2z_S_c = I_ERI_Fx2z_Px_Fx2z_S_c+ABZ*I_ERI_Dxz_Px_Fx2z_S_c;
  Double I_ERI_D2y_Dxz_Fx2z_S_c = I_ERI_F2yz_Px_Fx2z_S_c+ABZ*I_ERI_D2y_Px_Fx2z_S_c;
  Double I_ERI_Dyz_Dxz_Fx2z_S_c = I_ERI_Fy2z_Px_Fx2z_S_c+ABZ*I_ERI_Dyz_Px_Fx2z_S_c;
  Double I_ERI_D2z_Dxz_Fx2z_S_c = I_ERI_F3z_Px_Fx2z_S_c+ABZ*I_ERI_D2z_Px_Fx2z_S_c;
  Double I_ERI_D2x_D2y_Fx2z_S_c = I_ERI_F2xy_Py_Fx2z_S_c+ABY*I_ERI_D2x_Py_Fx2z_S_c;
  Double I_ERI_Dxy_D2y_Fx2z_S_c = I_ERI_Fx2y_Py_Fx2z_S_c+ABY*I_ERI_Dxy_Py_Fx2z_S_c;
  Double I_ERI_Dxz_D2y_Fx2z_S_c = I_ERI_Fxyz_Py_Fx2z_S_c+ABY*I_ERI_Dxz_Py_Fx2z_S_c;
  Double I_ERI_D2y_D2y_Fx2z_S_c = I_ERI_F3y_Py_Fx2z_S_c+ABY*I_ERI_D2y_Py_Fx2z_S_c;
  Double I_ERI_Dyz_D2y_Fx2z_S_c = I_ERI_F2yz_Py_Fx2z_S_c+ABY*I_ERI_Dyz_Py_Fx2z_S_c;
  Double I_ERI_D2z_D2y_Fx2z_S_c = I_ERI_Fy2z_Py_Fx2z_S_c+ABY*I_ERI_D2z_Py_Fx2z_S_c;
  Double I_ERI_D2x_Dyz_Fx2z_S_c = I_ERI_F2xz_Py_Fx2z_S_c+ABZ*I_ERI_D2x_Py_Fx2z_S_c;
  Double I_ERI_Dxy_Dyz_Fx2z_S_c = I_ERI_Fxyz_Py_Fx2z_S_c+ABZ*I_ERI_Dxy_Py_Fx2z_S_c;
  Double I_ERI_Dxz_Dyz_Fx2z_S_c = I_ERI_Fx2z_Py_Fx2z_S_c+ABZ*I_ERI_Dxz_Py_Fx2z_S_c;
  Double I_ERI_D2y_Dyz_Fx2z_S_c = I_ERI_F2yz_Py_Fx2z_S_c+ABZ*I_ERI_D2y_Py_Fx2z_S_c;
  Double I_ERI_Dyz_Dyz_Fx2z_S_c = I_ERI_Fy2z_Py_Fx2z_S_c+ABZ*I_ERI_Dyz_Py_Fx2z_S_c;
  Double I_ERI_D2z_Dyz_Fx2z_S_c = I_ERI_F3z_Py_Fx2z_S_c+ABZ*I_ERI_D2z_Py_Fx2z_S_c;
  Double I_ERI_D2x_D2z_Fx2z_S_c = I_ERI_F2xz_Pz_Fx2z_S_c+ABZ*I_ERI_D2x_Pz_Fx2z_S_c;
  Double I_ERI_Dxy_D2z_Fx2z_S_c = I_ERI_Fxyz_Pz_Fx2z_S_c+ABZ*I_ERI_Dxy_Pz_Fx2z_S_c;
  Double I_ERI_Dxz_D2z_Fx2z_S_c = I_ERI_Fx2z_Pz_Fx2z_S_c+ABZ*I_ERI_Dxz_Pz_Fx2z_S_c;
  Double I_ERI_D2y_D2z_Fx2z_S_c = I_ERI_F2yz_Pz_Fx2z_S_c+ABZ*I_ERI_D2y_Pz_Fx2z_S_c;
  Double I_ERI_Dyz_D2z_Fx2z_S_c = I_ERI_Fy2z_Pz_Fx2z_S_c+ABZ*I_ERI_Dyz_Pz_Fx2z_S_c;
  Double I_ERI_D2z_D2z_Fx2z_S_c = I_ERI_F3z_Pz_Fx2z_S_c+ABZ*I_ERI_D2z_Pz_Fx2z_S_c;
  Double I_ERI_D2x_D2x_F3y_S_c = I_ERI_F3x_Px_F3y_S_c+ABX*I_ERI_D2x_Px_F3y_S_c;
  Double I_ERI_Dxy_D2x_F3y_S_c = I_ERI_F2xy_Px_F3y_S_c+ABX*I_ERI_Dxy_Px_F3y_S_c;
  Double I_ERI_Dxz_D2x_F3y_S_c = I_ERI_F2xz_Px_F3y_S_c+ABX*I_ERI_Dxz_Px_F3y_S_c;
  Double I_ERI_D2y_D2x_F3y_S_c = I_ERI_Fx2y_Px_F3y_S_c+ABX*I_ERI_D2y_Px_F3y_S_c;
  Double I_ERI_Dyz_D2x_F3y_S_c = I_ERI_Fxyz_Px_F3y_S_c+ABX*I_ERI_Dyz_Px_F3y_S_c;
  Double I_ERI_D2z_D2x_F3y_S_c = I_ERI_Fx2z_Px_F3y_S_c+ABX*I_ERI_D2z_Px_F3y_S_c;
  Double I_ERI_D2x_Dxy_F3y_S_c = I_ERI_F2xy_Px_F3y_S_c+ABY*I_ERI_D2x_Px_F3y_S_c;
  Double I_ERI_Dxy_Dxy_F3y_S_c = I_ERI_Fx2y_Px_F3y_S_c+ABY*I_ERI_Dxy_Px_F3y_S_c;
  Double I_ERI_Dxz_Dxy_F3y_S_c = I_ERI_Fxyz_Px_F3y_S_c+ABY*I_ERI_Dxz_Px_F3y_S_c;
  Double I_ERI_D2y_Dxy_F3y_S_c = I_ERI_F3y_Px_F3y_S_c+ABY*I_ERI_D2y_Px_F3y_S_c;
  Double I_ERI_Dyz_Dxy_F3y_S_c = I_ERI_F2yz_Px_F3y_S_c+ABY*I_ERI_Dyz_Px_F3y_S_c;
  Double I_ERI_D2z_Dxy_F3y_S_c = I_ERI_Fy2z_Px_F3y_S_c+ABY*I_ERI_D2z_Px_F3y_S_c;
  Double I_ERI_D2x_Dxz_F3y_S_c = I_ERI_F2xz_Px_F3y_S_c+ABZ*I_ERI_D2x_Px_F3y_S_c;
  Double I_ERI_Dxy_Dxz_F3y_S_c = I_ERI_Fxyz_Px_F3y_S_c+ABZ*I_ERI_Dxy_Px_F3y_S_c;
  Double I_ERI_Dxz_Dxz_F3y_S_c = I_ERI_Fx2z_Px_F3y_S_c+ABZ*I_ERI_Dxz_Px_F3y_S_c;
  Double I_ERI_D2y_Dxz_F3y_S_c = I_ERI_F2yz_Px_F3y_S_c+ABZ*I_ERI_D2y_Px_F3y_S_c;
  Double I_ERI_Dyz_Dxz_F3y_S_c = I_ERI_Fy2z_Px_F3y_S_c+ABZ*I_ERI_Dyz_Px_F3y_S_c;
  Double I_ERI_D2z_Dxz_F3y_S_c = I_ERI_F3z_Px_F3y_S_c+ABZ*I_ERI_D2z_Px_F3y_S_c;
  Double I_ERI_D2x_D2y_F3y_S_c = I_ERI_F2xy_Py_F3y_S_c+ABY*I_ERI_D2x_Py_F3y_S_c;
  Double I_ERI_Dxy_D2y_F3y_S_c = I_ERI_Fx2y_Py_F3y_S_c+ABY*I_ERI_Dxy_Py_F3y_S_c;
  Double I_ERI_Dxz_D2y_F3y_S_c = I_ERI_Fxyz_Py_F3y_S_c+ABY*I_ERI_Dxz_Py_F3y_S_c;
  Double I_ERI_D2y_D2y_F3y_S_c = I_ERI_F3y_Py_F3y_S_c+ABY*I_ERI_D2y_Py_F3y_S_c;
  Double I_ERI_Dyz_D2y_F3y_S_c = I_ERI_F2yz_Py_F3y_S_c+ABY*I_ERI_Dyz_Py_F3y_S_c;
  Double I_ERI_D2z_D2y_F3y_S_c = I_ERI_Fy2z_Py_F3y_S_c+ABY*I_ERI_D2z_Py_F3y_S_c;
  Double I_ERI_D2x_Dyz_F3y_S_c = I_ERI_F2xz_Py_F3y_S_c+ABZ*I_ERI_D2x_Py_F3y_S_c;
  Double I_ERI_Dxy_Dyz_F3y_S_c = I_ERI_Fxyz_Py_F3y_S_c+ABZ*I_ERI_Dxy_Py_F3y_S_c;
  Double I_ERI_Dxz_Dyz_F3y_S_c = I_ERI_Fx2z_Py_F3y_S_c+ABZ*I_ERI_Dxz_Py_F3y_S_c;
  Double I_ERI_D2y_Dyz_F3y_S_c = I_ERI_F2yz_Py_F3y_S_c+ABZ*I_ERI_D2y_Py_F3y_S_c;
  Double I_ERI_Dyz_Dyz_F3y_S_c = I_ERI_Fy2z_Py_F3y_S_c+ABZ*I_ERI_Dyz_Py_F3y_S_c;
  Double I_ERI_D2z_Dyz_F3y_S_c = I_ERI_F3z_Py_F3y_S_c+ABZ*I_ERI_D2z_Py_F3y_S_c;
  Double I_ERI_D2x_D2z_F3y_S_c = I_ERI_F2xz_Pz_F3y_S_c+ABZ*I_ERI_D2x_Pz_F3y_S_c;
  Double I_ERI_Dxy_D2z_F3y_S_c = I_ERI_Fxyz_Pz_F3y_S_c+ABZ*I_ERI_Dxy_Pz_F3y_S_c;
  Double I_ERI_Dxz_D2z_F3y_S_c = I_ERI_Fx2z_Pz_F3y_S_c+ABZ*I_ERI_Dxz_Pz_F3y_S_c;
  Double I_ERI_D2y_D2z_F3y_S_c = I_ERI_F2yz_Pz_F3y_S_c+ABZ*I_ERI_D2y_Pz_F3y_S_c;
  Double I_ERI_Dyz_D2z_F3y_S_c = I_ERI_Fy2z_Pz_F3y_S_c+ABZ*I_ERI_Dyz_Pz_F3y_S_c;
  Double I_ERI_D2z_D2z_F3y_S_c = I_ERI_F3z_Pz_F3y_S_c+ABZ*I_ERI_D2z_Pz_F3y_S_c;
  Double I_ERI_D2x_D2x_F2yz_S_c = I_ERI_F3x_Px_F2yz_S_c+ABX*I_ERI_D2x_Px_F2yz_S_c;
  Double I_ERI_Dxy_D2x_F2yz_S_c = I_ERI_F2xy_Px_F2yz_S_c+ABX*I_ERI_Dxy_Px_F2yz_S_c;
  Double I_ERI_Dxz_D2x_F2yz_S_c = I_ERI_F2xz_Px_F2yz_S_c+ABX*I_ERI_Dxz_Px_F2yz_S_c;
  Double I_ERI_D2y_D2x_F2yz_S_c = I_ERI_Fx2y_Px_F2yz_S_c+ABX*I_ERI_D2y_Px_F2yz_S_c;
  Double I_ERI_Dyz_D2x_F2yz_S_c = I_ERI_Fxyz_Px_F2yz_S_c+ABX*I_ERI_Dyz_Px_F2yz_S_c;
  Double I_ERI_D2z_D2x_F2yz_S_c = I_ERI_Fx2z_Px_F2yz_S_c+ABX*I_ERI_D2z_Px_F2yz_S_c;
  Double I_ERI_D2x_Dxy_F2yz_S_c = I_ERI_F2xy_Px_F2yz_S_c+ABY*I_ERI_D2x_Px_F2yz_S_c;
  Double I_ERI_Dxy_Dxy_F2yz_S_c = I_ERI_Fx2y_Px_F2yz_S_c+ABY*I_ERI_Dxy_Px_F2yz_S_c;
  Double I_ERI_Dxz_Dxy_F2yz_S_c = I_ERI_Fxyz_Px_F2yz_S_c+ABY*I_ERI_Dxz_Px_F2yz_S_c;
  Double I_ERI_D2y_Dxy_F2yz_S_c = I_ERI_F3y_Px_F2yz_S_c+ABY*I_ERI_D2y_Px_F2yz_S_c;
  Double I_ERI_Dyz_Dxy_F2yz_S_c = I_ERI_F2yz_Px_F2yz_S_c+ABY*I_ERI_Dyz_Px_F2yz_S_c;
  Double I_ERI_D2z_Dxy_F2yz_S_c = I_ERI_Fy2z_Px_F2yz_S_c+ABY*I_ERI_D2z_Px_F2yz_S_c;
  Double I_ERI_D2x_Dxz_F2yz_S_c = I_ERI_F2xz_Px_F2yz_S_c+ABZ*I_ERI_D2x_Px_F2yz_S_c;
  Double I_ERI_Dxy_Dxz_F2yz_S_c = I_ERI_Fxyz_Px_F2yz_S_c+ABZ*I_ERI_Dxy_Px_F2yz_S_c;
  Double I_ERI_Dxz_Dxz_F2yz_S_c = I_ERI_Fx2z_Px_F2yz_S_c+ABZ*I_ERI_Dxz_Px_F2yz_S_c;
  Double I_ERI_D2y_Dxz_F2yz_S_c = I_ERI_F2yz_Px_F2yz_S_c+ABZ*I_ERI_D2y_Px_F2yz_S_c;
  Double I_ERI_Dyz_Dxz_F2yz_S_c = I_ERI_Fy2z_Px_F2yz_S_c+ABZ*I_ERI_Dyz_Px_F2yz_S_c;
  Double I_ERI_D2z_Dxz_F2yz_S_c = I_ERI_F3z_Px_F2yz_S_c+ABZ*I_ERI_D2z_Px_F2yz_S_c;
  Double I_ERI_D2x_D2y_F2yz_S_c = I_ERI_F2xy_Py_F2yz_S_c+ABY*I_ERI_D2x_Py_F2yz_S_c;
  Double I_ERI_Dxy_D2y_F2yz_S_c = I_ERI_Fx2y_Py_F2yz_S_c+ABY*I_ERI_Dxy_Py_F2yz_S_c;
  Double I_ERI_Dxz_D2y_F2yz_S_c = I_ERI_Fxyz_Py_F2yz_S_c+ABY*I_ERI_Dxz_Py_F2yz_S_c;
  Double I_ERI_D2y_D2y_F2yz_S_c = I_ERI_F3y_Py_F2yz_S_c+ABY*I_ERI_D2y_Py_F2yz_S_c;
  Double I_ERI_Dyz_D2y_F2yz_S_c = I_ERI_F2yz_Py_F2yz_S_c+ABY*I_ERI_Dyz_Py_F2yz_S_c;
  Double I_ERI_D2z_D2y_F2yz_S_c = I_ERI_Fy2z_Py_F2yz_S_c+ABY*I_ERI_D2z_Py_F2yz_S_c;
  Double I_ERI_D2x_Dyz_F2yz_S_c = I_ERI_F2xz_Py_F2yz_S_c+ABZ*I_ERI_D2x_Py_F2yz_S_c;
  Double I_ERI_Dxy_Dyz_F2yz_S_c = I_ERI_Fxyz_Py_F2yz_S_c+ABZ*I_ERI_Dxy_Py_F2yz_S_c;
  Double I_ERI_Dxz_Dyz_F2yz_S_c = I_ERI_Fx2z_Py_F2yz_S_c+ABZ*I_ERI_Dxz_Py_F2yz_S_c;
  Double I_ERI_D2y_Dyz_F2yz_S_c = I_ERI_F2yz_Py_F2yz_S_c+ABZ*I_ERI_D2y_Py_F2yz_S_c;
  Double I_ERI_Dyz_Dyz_F2yz_S_c = I_ERI_Fy2z_Py_F2yz_S_c+ABZ*I_ERI_Dyz_Py_F2yz_S_c;
  Double I_ERI_D2z_Dyz_F2yz_S_c = I_ERI_F3z_Py_F2yz_S_c+ABZ*I_ERI_D2z_Py_F2yz_S_c;
  Double I_ERI_D2x_D2z_F2yz_S_c = I_ERI_F2xz_Pz_F2yz_S_c+ABZ*I_ERI_D2x_Pz_F2yz_S_c;
  Double I_ERI_Dxy_D2z_F2yz_S_c = I_ERI_Fxyz_Pz_F2yz_S_c+ABZ*I_ERI_Dxy_Pz_F2yz_S_c;
  Double I_ERI_Dxz_D2z_F2yz_S_c = I_ERI_Fx2z_Pz_F2yz_S_c+ABZ*I_ERI_Dxz_Pz_F2yz_S_c;
  Double I_ERI_D2y_D2z_F2yz_S_c = I_ERI_F2yz_Pz_F2yz_S_c+ABZ*I_ERI_D2y_Pz_F2yz_S_c;
  Double I_ERI_Dyz_D2z_F2yz_S_c = I_ERI_Fy2z_Pz_F2yz_S_c+ABZ*I_ERI_Dyz_Pz_F2yz_S_c;
  Double I_ERI_D2z_D2z_F2yz_S_c = I_ERI_F3z_Pz_F2yz_S_c+ABZ*I_ERI_D2z_Pz_F2yz_S_c;
  Double I_ERI_D2x_D2x_Fy2z_S_c = I_ERI_F3x_Px_Fy2z_S_c+ABX*I_ERI_D2x_Px_Fy2z_S_c;
  Double I_ERI_Dxy_D2x_Fy2z_S_c = I_ERI_F2xy_Px_Fy2z_S_c+ABX*I_ERI_Dxy_Px_Fy2z_S_c;
  Double I_ERI_Dxz_D2x_Fy2z_S_c = I_ERI_F2xz_Px_Fy2z_S_c+ABX*I_ERI_Dxz_Px_Fy2z_S_c;
  Double I_ERI_D2y_D2x_Fy2z_S_c = I_ERI_Fx2y_Px_Fy2z_S_c+ABX*I_ERI_D2y_Px_Fy2z_S_c;
  Double I_ERI_Dyz_D2x_Fy2z_S_c = I_ERI_Fxyz_Px_Fy2z_S_c+ABX*I_ERI_Dyz_Px_Fy2z_S_c;
  Double I_ERI_D2z_D2x_Fy2z_S_c = I_ERI_Fx2z_Px_Fy2z_S_c+ABX*I_ERI_D2z_Px_Fy2z_S_c;
  Double I_ERI_D2x_Dxy_Fy2z_S_c = I_ERI_F2xy_Px_Fy2z_S_c+ABY*I_ERI_D2x_Px_Fy2z_S_c;
  Double I_ERI_Dxy_Dxy_Fy2z_S_c = I_ERI_Fx2y_Px_Fy2z_S_c+ABY*I_ERI_Dxy_Px_Fy2z_S_c;
  Double I_ERI_Dxz_Dxy_Fy2z_S_c = I_ERI_Fxyz_Px_Fy2z_S_c+ABY*I_ERI_Dxz_Px_Fy2z_S_c;
  Double I_ERI_D2y_Dxy_Fy2z_S_c = I_ERI_F3y_Px_Fy2z_S_c+ABY*I_ERI_D2y_Px_Fy2z_S_c;
  Double I_ERI_Dyz_Dxy_Fy2z_S_c = I_ERI_F2yz_Px_Fy2z_S_c+ABY*I_ERI_Dyz_Px_Fy2z_S_c;
  Double I_ERI_D2z_Dxy_Fy2z_S_c = I_ERI_Fy2z_Px_Fy2z_S_c+ABY*I_ERI_D2z_Px_Fy2z_S_c;
  Double I_ERI_D2x_Dxz_Fy2z_S_c = I_ERI_F2xz_Px_Fy2z_S_c+ABZ*I_ERI_D2x_Px_Fy2z_S_c;
  Double I_ERI_Dxy_Dxz_Fy2z_S_c = I_ERI_Fxyz_Px_Fy2z_S_c+ABZ*I_ERI_Dxy_Px_Fy2z_S_c;
  Double I_ERI_Dxz_Dxz_Fy2z_S_c = I_ERI_Fx2z_Px_Fy2z_S_c+ABZ*I_ERI_Dxz_Px_Fy2z_S_c;
  Double I_ERI_D2y_Dxz_Fy2z_S_c = I_ERI_F2yz_Px_Fy2z_S_c+ABZ*I_ERI_D2y_Px_Fy2z_S_c;
  Double I_ERI_Dyz_Dxz_Fy2z_S_c = I_ERI_Fy2z_Px_Fy2z_S_c+ABZ*I_ERI_Dyz_Px_Fy2z_S_c;
  Double I_ERI_D2z_Dxz_Fy2z_S_c = I_ERI_F3z_Px_Fy2z_S_c+ABZ*I_ERI_D2z_Px_Fy2z_S_c;
  Double I_ERI_D2x_D2y_Fy2z_S_c = I_ERI_F2xy_Py_Fy2z_S_c+ABY*I_ERI_D2x_Py_Fy2z_S_c;
  Double I_ERI_Dxy_D2y_Fy2z_S_c = I_ERI_Fx2y_Py_Fy2z_S_c+ABY*I_ERI_Dxy_Py_Fy2z_S_c;
  Double I_ERI_Dxz_D2y_Fy2z_S_c = I_ERI_Fxyz_Py_Fy2z_S_c+ABY*I_ERI_Dxz_Py_Fy2z_S_c;
  Double I_ERI_D2y_D2y_Fy2z_S_c = I_ERI_F3y_Py_Fy2z_S_c+ABY*I_ERI_D2y_Py_Fy2z_S_c;
  Double I_ERI_Dyz_D2y_Fy2z_S_c = I_ERI_F2yz_Py_Fy2z_S_c+ABY*I_ERI_Dyz_Py_Fy2z_S_c;
  Double I_ERI_D2z_D2y_Fy2z_S_c = I_ERI_Fy2z_Py_Fy2z_S_c+ABY*I_ERI_D2z_Py_Fy2z_S_c;
  Double I_ERI_D2x_Dyz_Fy2z_S_c = I_ERI_F2xz_Py_Fy2z_S_c+ABZ*I_ERI_D2x_Py_Fy2z_S_c;
  Double I_ERI_Dxy_Dyz_Fy2z_S_c = I_ERI_Fxyz_Py_Fy2z_S_c+ABZ*I_ERI_Dxy_Py_Fy2z_S_c;
  Double I_ERI_Dxz_Dyz_Fy2z_S_c = I_ERI_Fx2z_Py_Fy2z_S_c+ABZ*I_ERI_Dxz_Py_Fy2z_S_c;
  Double I_ERI_D2y_Dyz_Fy2z_S_c = I_ERI_F2yz_Py_Fy2z_S_c+ABZ*I_ERI_D2y_Py_Fy2z_S_c;
  Double I_ERI_Dyz_Dyz_Fy2z_S_c = I_ERI_Fy2z_Py_Fy2z_S_c+ABZ*I_ERI_Dyz_Py_Fy2z_S_c;
  Double I_ERI_D2z_Dyz_Fy2z_S_c = I_ERI_F3z_Py_Fy2z_S_c+ABZ*I_ERI_D2z_Py_Fy2z_S_c;
  Double I_ERI_D2x_D2z_Fy2z_S_c = I_ERI_F2xz_Pz_Fy2z_S_c+ABZ*I_ERI_D2x_Pz_Fy2z_S_c;
  Double I_ERI_Dxy_D2z_Fy2z_S_c = I_ERI_Fxyz_Pz_Fy2z_S_c+ABZ*I_ERI_Dxy_Pz_Fy2z_S_c;
  Double I_ERI_Dxz_D2z_Fy2z_S_c = I_ERI_Fx2z_Pz_Fy2z_S_c+ABZ*I_ERI_Dxz_Pz_Fy2z_S_c;
  Double I_ERI_D2y_D2z_Fy2z_S_c = I_ERI_F2yz_Pz_Fy2z_S_c+ABZ*I_ERI_D2y_Pz_Fy2z_S_c;
  Double I_ERI_Dyz_D2z_Fy2z_S_c = I_ERI_Fy2z_Pz_Fy2z_S_c+ABZ*I_ERI_Dyz_Pz_Fy2z_S_c;
  Double I_ERI_D2z_D2z_Fy2z_S_c = I_ERI_F3z_Pz_Fy2z_S_c+ABZ*I_ERI_D2z_Pz_Fy2z_S_c;
  Double I_ERI_D2x_D2x_F3z_S_c = I_ERI_F3x_Px_F3z_S_c+ABX*I_ERI_D2x_Px_F3z_S_c;
  Double I_ERI_Dxy_D2x_F3z_S_c = I_ERI_F2xy_Px_F3z_S_c+ABX*I_ERI_Dxy_Px_F3z_S_c;
  Double I_ERI_Dxz_D2x_F3z_S_c = I_ERI_F2xz_Px_F3z_S_c+ABX*I_ERI_Dxz_Px_F3z_S_c;
  Double I_ERI_D2y_D2x_F3z_S_c = I_ERI_Fx2y_Px_F3z_S_c+ABX*I_ERI_D2y_Px_F3z_S_c;
  Double I_ERI_Dyz_D2x_F3z_S_c = I_ERI_Fxyz_Px_F3z_S_c+ABX*I_ERI_Dyz_Px_F3z_S_c;
  Double I_ERI_D2z_D2x_F3z_S_c = I_ERI_Fx2z_Px_F3z_S_c+ABX*I_ERI_D2z_Px_F3z_S_c;
  Double I_ERI_D2x_Dxy_F3z_S_c = I_ERI_F2xy_Px_F3z_S_c+ABY*I_ERI_D2x_Px_F3z_S_c;
  Double I_ERI_Dxy_Dxy_F3z_S_c = I_ERI_Fx2y_Px_F3z_S_c+ABY*I_ERI_Dxy_Px_F3z_S_c;
  Double I_ERI_Dxz_Dxy_F3z_S_c = I_ERI_Fxyz_Px_F3z_S_c+ABY*I_ERI_Dxz_Px_F3z_S_c;
  Double I_ERI_D2y_Dxy_F3z_S_c = I_ERI_F3y_Px_F3z_S_c+ABY*I_ERI_D2y_Px_F3z_S_c;
  Double I_ERI_Dyz_Dxy_F3z_S_c = I_ERI_F2yz_Px_F3z_S_c+ABY*I_ERI_Dyz_Px_F3z_S_c;
  Double I_ERI_D2z_Dxy_F3z_S_c = I_ERI_Fy2z_Px_F3z_S_c+ABY*I_ERI_D2z_Px_F3z_S_c;
  Double I_ERI_D2x_Dxz_F3z_S_c = I_ERI_F2xz_Px_F3z_S_c+ABZ*I_ERI_D2x_Px_F3z_S_c;
  Double I_ERI_Dxy_Dxz_F3z_S_c = I_ERI_Fxyz_Px_F3z_S_c+ABZ*I_ERI_Dxy_Px_F3z_S_c;
  Double I_ERI_Dxz_Dxz_F3z_S_c = I_ERI_Fx2z_Px_F3z_S_c+ABZ*I_ERI_Dxz_Px_F3z_S_c;
  Double I_ERI_D2y_Dxz_F3z_S_c = I_ERI_F2yz_Px_F3z_S_c+ABZ*I_ERI_D2y_Px_F3z_S_c;
  Double I_ERI_Dyz_Dxz_F3z_S_c = I_ERI_Fy2z_Px_F3z_S_c+ABZ*I_ERI_Dyz_Px_F3z_S_c;
  Double I_ERI_D2z_Dxz_F3z_S_c = I_ERI_F3z_Px_F3z_S_c+ABZ*I_ERI_D2z_Px_F3z_S_c;
  Double I_ERI_D2x_D2y_F3z_S_c = I_ERI_F2xy_Py_F3z_S_c+ABY*I_ERI_D2x_Py_F3z_S_c;
  Double I_ERI_Dxy_D2y_F3z_S_c = I_ERI_Fx2y_Py_F3z_S_c+ABY*I_ERI_Dxy_Py_F3z_S_c;
  Double I_ERI_Dxz_D2y_F3z_S_c = I_ERI_Fxyz_Py_F3z_S_c+ABY*I_ERI_Dxz_Py_F3z_S_c;
  Double I_ERI_D2y_D2y_F3z_S_c = I_ERI_F3y_Py_F3z_S_c+ABY*I_ERI_D2y_Py_F3z_S_c;
  Double I_ERI_Dyz_D2y_F3z_S_c = I_ERI_F2yz_Py_F3z_S_c+ABY*I_ERI_Dyz_Py_F3z_S_c;
  Double I_ERI_D2z_D2y_F3z_S_c = I_ERI_Fy2z_Py_F3z_S_c+ABY*I_ERI_D2z_Py_F3z_S_c;
  Double I_ERI_D2x_Dyz_F3z_S_c = I_ERI_F2xz_Py_F3z_S_c+ABZ*I_ERI_D2x_Py_F3z_S_c;
  Double I_ERI_Dxy_Dyz_F3z_S_c = I_ERI_Fxyz_Py_F3z_S_c+ABZ*I_ERI_Dxy_Py_F3z_S_c;
  Double I_ERI_Dxz_Dyz_F3z_S_c = I_ERI_Fx2z_Py_F3z_S_c+ABZ*I_ERI_Dxz_Py_F3z_S_c;
  Double I_ERI_D2y_Dyz_F3z_S_c = I_ERI_F2yz_Py_F3z_S_c+ABZ*I_ERI_D2y_Py_F3z_S_c;
  Double I_ERI_Dyz_Dyz_F3z_S_c = I_ERI_Fy2z_Py_F3z_S_c+ABZ*I_ERI_Dyz_Py_F3z_S_c;
  Double I_ERI_D2z_Dyz_F3z_S_c = I_ERI_F3z_Py_F3z_S_c+ABZ*I_ERI_D2z_Py_F3z_S_c;
  Double I_ERI_D2x_D2z_F3z_S_c = I_ERI_F2xz_Pz_F3z_S_c+ABZ*I_ERI_D2x_Pz_F3z_S_c;
  Double I_ERI_Dxy_D2z_F3z_S_c = I_ERI_Fxyz_Pz_F3z_S_c+ABZ*I_ERI_Dxy_Pz_F3z_S_c;
  Double I_ERI_Dxz_D2z_F3z_S_c = I_ERI_Fx2z_Pz_F3z_S_c+ABZ*I_ERI_Dxz_Pz_F3z_S_c;
  Double I_ERI_D2y_D2z_F3z_S_c = I_ERI_F2yz_Pz_F3z_S_c+ABZ*I_ERI_D2y_Pz_F3z_S_c;
  Double I_ERI_Dyz_D2z_F3z_S_c = I_ERI_Fy2z_Pz_F3z_S_c+ABZ*I_ERI_Dyz_Pz_F3z_S_c;
  Double I_ERI_D2z_D2z_F3z_S_c = I_ERI_F3z_Pz_F3z_S_c+ABZ*I_ERI_D2z_Pz_F3z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_a
   * RHS shell quartet name: SQ_ERI_P_D_D_S
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_F3x_D2x_D2x_S_a-2*I_ERI_Px_D2x_D2x_S;
  abcd[1] = 2.0E0*I_ERI_F2xy_D2x_D2x_S_a-1*I_ERI_Py_D2x_D2x_S;
  abcd[2] = 2.0E0*I_ERI_F2xz_D2x_D2x_S_a-1*I_ERI_Pz_D2x_D2x_S;
  abcd[3] = 2.0E0*I_ERI_Fx2y_D2x_D2x_S_a;
  abcd[4] = 2.0E0*I_ERI_Fxyz_D2x_D2x_S_a;
  abcd[5] = 2.0E0*I_ERI_Fx2z_D2x_D2x_S_a;
  abcd[6] = 2.0E0*I_ERI_F3x_Dxy_D2x_S_a-2*I_ERI_Px_Dxy_D2x_S;
  abcd[7] = 2.0E0*I_ERI_F2xy_Dxy_D2x_S_a-1*I_ERI_Py_Dxy_D2x_S;
  abcd[8] = 2.0E0*I_ERI_F2xz_Dxy_D2x_S_a-1*I_ERI_Pz_Dxy_D2x_S;
  abcd[9] = 2.0E0*I_ERI_Fx2y_Dxy_D2x_S_a;
  abcd[10] = 2.0E0*I_ERI_Fxyz_Dxy_D2x_S_a;
  abcd[11] = 2.0E0*I_ERI_Fx2z_Dxy_D2x_S_a;
  abcd[12] = 2.0E0*I_ERI_F3x_Dxz_D2x_S_a-2*I_ERI_Px_Dxz_D2x_S;
  abcd[13] = 2.0E0*I_ERI_F2xy_Dxz_D2x_S_a-1*I_ERI_Py_Dxz_D2x_S;
  abcd[14] = 2.0E0*I_ERI_F2xz_Dxz_D2x_S_a-1*I_ERI_Pz_Dxz_D2x_S;
  abcd[15] = 2.0E0*I_ERI_Fx2y_Dxz_D2x_S_a;
  abcd[16] = 2.0E0*I_ERI_Fxyz_Dxz_D2x_S_a;
  abcd[17] = 2.0E0*I_ERI_Fx2z_Dxz_D2x_S_a;
  abcd[18] = 2.0E0*I_ERI_F3x_D2y_D2x_S_a-2*I_ERI_Px_D2y_D2x_S;
  abcd[19] = 2.0E0*I_ERI_F2xy_D2y_D2x_S_a-1*I_ERI_Py_D2y_D2x_S;
  abcd[20] = 2.0E0*I_ERI_F2xz_D2y_D2x_S_a-1*I_ERI_Pz_D2y_D2x_S;
  abcd[21] = 2.0E0*I_ERI_Fx2y_D2y_D2x_S_a;
  abcd[22] = 2.0E0*I_ERI_Fxyz_D2y_D2x_S_a;
  abcd[23] = 2.0E0*I_ERI_Fx2z_D2y_D2x_S_a;
  abcd[24] = 2.0E0*I_ERI_F3x_Dyz_D2x_S_a-2*I_ERI_Px_Dyz_D2x_S;
  abcd[25] = 2.0E0*I_ERI_F2xy_Dyz_D2x_S_a-1*I_ERI_Py_Dyz_D2x_S;
  abcd[26] = 2.0E0*I_ERI_F2xz_Dyz_D2x_S_a-1*I_ERI_Pz_Dyz_D2x_S;
  abcd[27] = 2.0E0*I_ERI_Fx2y_Dyz_D2x_S_a;
  abcd[28] = 2.0E0*I_ERI_Fxyz_Dyz_D2x_S_a;
  abcd[29] = 2.0E0*I_ERI_Fx2z_Dyz_D2x_S_a;
  abcd[30] = 2.0E0*I_ERI_F3x_D2z_D2x_S_a-2*I_ERI_Px_D2z_D2x_S;
  abcd[31] = 2.0E0*I_ERI_F2xy_D2z_D2x_S_a-1*I_ERI_Py_D2z_D2x_S;
  abcd[32] = 2.0E0*I_ERI_F2xz_D2z_D2x_S_a-1*I_ERI_Pz_D2z_D2x_S;
  abcd[33] = 2.0E0*I_ERI_Fx2y_D2z_D2x_S_a;
  abcd[34] = 2.0E0*I_ERI_Fxyz_D2z_D2x_S_a;
  abcd[35] = 2.0E0*I_ERI_Fx2z_D2z_D2x_S_a;
  abcd[36] = 2.0E0*I_ERI_F3x_D2x_Dxy_S_a-2*I_ERI_Px_D2x_Dxy_S;
  abcd[37] = 2.0E0*I_ERI_F2xy_D2x_Dxy_S_a-1*I_ERI_Py_D2x_Dxy_S;
  abcd[38] = 2.0E0*I_ERI_F2xz_D2x_Dxy_S_a-1*I_ERI_Pz_D2x_Dxy_S;
  abcd[39] = 2.0E0*I_ERI_Fx2y_D2x_Dxy_S_a;
  abcd[40] = 2.0E0*I_ERI_Fxyz_D2x_Dxy_S_a;
  abcd[41] = 2.0E0*I_ERI_Fx2z_D2x_Dxy_S_a;
  abcd[42] = 2.0E0*I_ERI_F3x_Dxy_Dxy_S_a-2*I_ERI_Px_Dxy_Dxy_S;
  abcd[43] = 2.0E0*I_ERI_F2xy_Dxy_Dxy_S_a-1*I_ERI_Py_Dxy_Dxy_S;
  abcd[44] = 2.0E0*I_ERI_F2xz_Dxy_Dxy_S_a-1*I_ERI_Pz_Dxy_Dxy_S;
  abcd[45] = 2.0E0*I_ERI_Fx2y_Dxy_Dxy_S_a;
  abcd[46] = 2.0E0*I_ERI_Fxyz_Dxy_Dxy_S_a;
  abcd[47] = 2.0E0*I_ERI_Fx2z_Dxy_Dxy_S_a;
  abcd[48] = 2.0E0*I_ERI_F3x_Dxz_Dxy_S_a-2*I_ERI_Px_Dxz_Dxy_S;
  abcd[49] = 2.0E0*I_ERI_F2xy_Dxz_Dxy_S_a-1*I_ERI_Py_Dxz_Dxy_S;
  abcd[50] = 2.0E0*I_ERI_F2xz_Dxz_Dxy_S_a-1*I_ERI_Pz_Dxz_Dxy_S;
  abcd[51] = 2.0E0*I_ERI_Fx2y_Dxz_Dxy_S_a;
  abcd[52] = 2.0E0*I_ERI_Fxyz_Dxz_Dxy_S_a;
  abcd[53] = 2.0E0*I_ERI_Fx2z_Dxz_Dxy_S_a;
  abcd[54] = 2.0E0*I_ERI_F3x_D2y_Dxy_S_a-2*I_ERI_Px_D2y_Dxy_S;
  abcd[55] = 2.0E0*I_ERI_F2xy_D2y_Dxy_S_a-1*I_ERI_Py_D2y_Dxy_S;
  abcd[56] = 2.0E0*I_ERI_F2xz_D2y_Dxy_S_a-1*I_ERI_Pz_D2y_Dxy_S;
  abcd[57] = 2.0E0*I_ERI_Fx2y_D2y_Dxy_S_a;
  abcd[58] = 2.0E0*I_ERI_Fxyz_D2y_Dxy_S_a;
  abcd[59] = 2.0E0*I_ERI_Fx2z_D2y_Dxy_S_a;
  abcd[60] = 2.0E0*I_ERI_F3x_Dyz_Dxy_S_a-2*I_ERI_Px_Dyz_Dxy_S;
  abcd[61] = 2.0E0*I_ERI_F2xy_Dyz_Dxy_S_a-1*I_ERI_Py_Dyz_Dxy_S;
  abcd[62] = 2.0E0*I_ERI_F2xz_Dyz_Dxy_S_a-1*I_ERI_Pz_Dyz_Dxy_S;
  abcd[63] = 2.0E0*I_ERI_Fx2y_Dyz_Dxy_S_a;
  abcd[64] = 2.0E0*I_ERI_Fxyz_Dyz_Dxy_S_a;
  abcd[65] = 2.0E0*I_ERI_Fx2z_Dyz_Dxy_S_a;
  abcd[66] = 2.0E0*I_ERI_F3x_D2z_Dxy_S_a-2*I_ERI_Px_D2z_Dxy_S;
  abcd[67] = 2.0E0*I_ERI_F2xy_D2z_Dxy_S_a-1*I_ERI_Py_D2z_Dxy_S;
  abcd[68] = 2.0E0*I_ERI_F2xz_D2z_Dxy_S_a-1*I_ERI_Pz_D2z_Dxy_S;
  abcd[69] = 2.0E0*I_ERI_Fx2y_D2z_Dxy_S_a;
  abcd[70] = 2.0E0*I_ERI_Fxyz_D2z_Dxy_S_a;
  abcd[71] = 2.0E0*I_ERI_Fx2z_D2z_Dxy_S_a;
  abcd[72] = 2.0E0*I_ERI_F3x_D2x_Dxz_S_a-2*I_ERI_Px_D2x_Dxz_S;
  abcd[73] = 2.0E0*I_ERI_F2xy_D2x_Dxz_S_a-1*I_ERI_Py_D2x_Dxz_S;
  abcd[74] = 2.0E0*I_ERI_F2xz_D2x_Dxz_S_a-1*I_ERI_Pz_D2x_Dxz_S;
  abcd[75] = 2.0E0*I_ERI_Fx2y_D2x_Dxz_S_a;
  abcd[76] = 2.0E0*I_ERI_Fxyz_D2x_Dxz_S_a;
  abcd[77] = 2.0E0*I_ERI_Fx2z_D2x_Dxz_S_a;
  abcd[78] = 2.0E0*I_ERI_F3x_Dxy_Dxz_S_a-2*I_ERI_Px_Dxy_Dxz_S;
  abcd[79] = 2.0E0*I_ERI_F2xy_Dxy_Dxz_S_a-1*I_ERI_Py_Dxy_Dxz_S;
  abcd[80] = 2.0E0*I_ERI_F2xz_Dxy_Dxz_S_a-1*I_ERI_Pz_Dxy_Dxz_S;
  abcd[81] = 2.0E0*I_ERI_Fx2y_Dxy_Dxz_S_a;
  abcd[82] = 2.0E0*I_ERI_Fxyz_Dxy_Dxz_S_a;
  abcd[83] = 2.0E0*I_ERI_Fx2z_Dxy_Dxz_S_a;
  abcd[84] = 2.0E0*I_ERI_F3x_Dxz_Dxz_S_a-2*I_ERI_Px_Dxz_Dxz_S;
  abcd[85] = 2.0E0*I_ERI_F2xy_Dxz_Dxz_S_a-1*I_ERI_Py_Dxz_Dxz_S;
  abcd[86] = 2.0E0*I_ERI_F2xz_Dxz_Dxz_S_a-1*I_ERI_Pz_Dxz_Dxz_S;
  abcd[87] = 2.0E0*I_ERI_Fx2y_Dxz_Dxz_S_a;
  abcd[88] = 2.0E0*I_ERI_Fxyz_Dxz_Dxz_S_a;
  abcd[89] = 2.0E0*I_ERI_Fx2z_Dxz_Dxz_S_a;
  abcd[90] = 2.0E0*I_ERI_F3x_D2y_Dxz_S_a-2*I_ERI_Px_D2y_Dxz_S;
  abcd[91] = 2.0E0*I_ERI_F2xy_D2y_Dxz_S_a-1*I_ERI_Py_D2y_Dxz_S;
  abcd[92] = 2.0E0*I_ERI_F2xz_D2y_Dxz_S_a-1*I_ERI_Pz_D2y_Dxz_S;
  abcd[93] = 2.0E0*I_ERI_Fx2y_D2y_Dxz_S_a;
  abcd[94] = 2.0E0*I_ERI_Fxyz_D2y_Dxz_S_a;
  abcd[95] = 2.0E0*I_ERI_Fx2z_D2y_Dxz_S_a;
  abcd[96] = 2.0E0*I_ERI_F3x_Dyz_Dxz_S_a-2*I_ERI_Px_Dyz_Dxz_S;
  abcd[97] = 2.0E0*I_ERI_F2xy_Dyz_Dxz_S_a-1*I_ERI_Py_Dyz_Dxz_S;
  abcd[98] = 2.0E0*I_ERI_F2xz_Dyz_Dxz_S_a-1*I_ERI_Pz_Dyz_Dxz_S;
  abcd[99] = 2.0E0*I_ERI_Fx2y_Dyz_Dxz_S_a;
  abcd[100] = 2.0E0*I_ERI_Fxyz_Dyz_Dxz_S_a;
  abcd[101] = 2.0E0*I_ERI_Fx2z_Dyz_Dxz_S_a;
  abcd[102] = 2.0E0*I_ERI_F3x_D2z_Dxz_S_a-2*I_ERI_Px_D2z_Dxz_S;
  abcd[103] = 2.0E0*I_ERI_F2xy_D2z_Dxz_S_a-1*I_ERI_Py_D2z_Dxz_S;
  abcd[104] = 2.0E0*I_ERI_F2xz_D2z_Dxz_S_a-1*I_ERI_Pz_D2z_Dxz_S;
  abcd[105] = 2.0E0*I_ERI_Fx2y_D2z_Dxz_S_a;
  abcd[106] = 2.0E0*I_ERI_Fxyz_D2z_Dxz_S_a;
  abcd[107] = 2.0E0*I_ERI_Fx2z_D2z_Dxz_S_a;
  abcd[108] = 2.0E0*I_ERI_F3x_D2x_D2y_S_a-2*I_ERI_Px_D2x_D2y_S;
  abcd[109] = 2.0E0*I_ERI_F2xy_D2x_D2y_S_a-1*I_ERI_Py_D2x_D2y_S;
  abcd[110] = 2.0E0*I_ERI_F2xz_D2x_D2y_S_a-1*I_ERI_Pz_D2x_D2y_S;
  abcd[111] = 2.0E0*I_ERI_Fx2y_D2x_D2y_S_a;
  abcd[112] = 2.0E0*I_ERI_Fxyz_D2x_D2y_S_a;
  abcd[113] = 2.0E0*I_ERI_Fx2z_D2x_D2y_S_a;
  abcd[114] = 2.0E0*I_ERI_F3x_Dxy_D2y_S_a-2*I_ERI_Px_Dxy_D2y_S;
  abcd[115] = 2.0E0*I_ERI_F2xy_Dxy_D2y_S_a-1*I_ERI_Py_Dxy_D2y_S;
  abcd[116] = 2.0E0*I_ERI_F2xz_Dxy_D2y_S_a-1*I_ERI_Pz_Dxy_D2y_S;
  abcd[117] = 2.0E0*I_ERI_Fx2y_Dxy_D2y_S_a;
  abcd[118] = 2.0E0*I_ERI_Fxyz_Dxy_D2y_S_a;
  abcd[119] = 2.0E0*I_ERI_Fx2z_Dxy_D2y_S_a;
  abcd[120] = 2.0E0*I_ERI_F3x_Dxz_D2y_S_a-2*I_ERI_Px_Dxz_D2y_S;
  abcd[121] = 2.0E0*I_ERI_F2xy_Dxz_D2y_S_a-1*I_ERI_Py_Dxz_D2y_S;
  abcd[122] = 2.0E0*I_ERI_F2xz_Dxz_D2y_S_a-1*I_ERI_Pz_Dxz_D2y_S;
  abcd[123] = 2.0E0*I_ERI_Fx2y_Dxz_D2y_S_a;
  abcd[124] = 2.0E0*I_ERI_Fxyz_Dxz_D2y_S_a;
  abcd[125] = 2.0E0*I_ERI_Fx2z_Dxz_D2y_S_a;
  abcd[126] = 2.0E0*I_ERI_F3x_D2y_D2y_S_a-2*I_ERI_Px_D2y_D2y_S;
  abcd[127] = 2.0E0*I_ERI_F2xy_D2y_D2y_S_a-1*I_ERI_Py_D2y_D2y_S;
  abcd[128] = 2.0E0*I_ERI_F2xz_D2y_D2y_S_a-1*I_ERI_Pz_D2y_D2y_S;
  abcd[129] = 2.0E0*I_ERI_Fx2y_D2y_D2y_S_a;
  abcd[130] = 2.0E0*I_ERI_Fxyz_D2y_D2y_S_a;
  abcd[131] = 2.0E0*I_ERI_Fx2z_D2y_D2y_S_a;
  abcd[132] = 2.0E0*I_ERI_F3x_Dyz_D2y_S_a-2*I_ERI_Px_Dyz_D2y_S;
  abcd[133] = 2.0E0*I_ERI_F2xy_Dyz_D2y_S_a-1*I_ERI_Py_Dyz_D2y_S;
  abcd[134] = 2.0E0*I_ERI_F2xz_Dyz_D2y_S_a-1*I_ERI_Pz_Dyz_D2y_S;
  abcd[135] = 2.0E0*I_ERI_Fx2y_Dyz_D2y_S_a;
  abcd[136] = 2.0E0*I_ERI_Fxyz_Dyz_D2y_S_a;
  abcd[137] = 2.0E0*I_ERI_Fx2z_Dyz_D2y_S_a;
  abcd[138] = 2.0E0*I_ERI_F3x_D2z_D2y_S_a-2*I_ERI_Px_D2z_D2y_S;
  abcd[139] = 2.0E0*I_ERI_F2xy_D2z_D2y_S_a-1*I_ERI_Py_D2z_D2y_S;
  abcd[140] = 2.0E0*I_ERI_F2xz_D2z_D2y_S_a-1*I_ERI_Pz_D2z_D2y_S;
  abcd[141] = 2.0E0*I_ERI_Fx2y_D2z_D2y_S_a;
  abcd[142] = 2.0E0*I_ERI_Fxyz_D2z_D2y_S_a;
  abcd[143] = 2.0E0*I_ERI_Fx2z_D2z_D2y_S_a;
  abcd[144] = 2.0E0*I_ERI_F3x_D2x_Dyz_S_a-2*I_ERI_Px_D2x_Dyz_S;
  abcd[145] = 2.0E0*I_ERI_F2xy_D2x_Dyz_S_a-1*I_ERI_Py_D2x_Dyz_S;
  abcd[146] = 2.0E0*I_ERI_F2xz_D2x_Dyz_S_a-1*I_ERI_Pz_D2x_Dyz_S;
  abcd[147] = 2.0E0*I_ERI_Fx2y_D2x_Dyz_S_a;
  abcd[148] = 2.0E0*I_ERI_Fxyz_D2x_Dyz_S_a;
  abcd[149] = 2.0E0*I_ERI_Fx2z_D2x_Dyz_S_a;
  abcd[150] = 2.0E0*I_ERI_F3x_Dxy_Dyz_S_a-2*I_ERI_Px_Dxy_Dyz_S;
  abcd[151] = 2.0E0*I_ERI_F2xy_Dxy_Dyz_S_a-1*I_ERI_Py_Dxy_Dyz_S;
  abcd[152] = 2.0E0*I_ERI_F2xz_Dxy_Dyz_S_a-1*I_ERI_Pz_Dxy_Dyz_S;
  abcd[153] = 2.0E0*I_ERI_Fx2y_Dxy_Dyz_S_a;
  abcd[154] = 2.0E0*I_ERI_Fxyz_Dxy_Dyz_S_a;
  abcd[155] = 2.0E0*I_ERI_Fx2z_Dxy_Dyz_S_a;
  abcd[156] = 2.0E0*I_ERI_F3x_Dxz_Dyz_S_a-2*I_ERI_Px_Dxz_Dyz_S;
  abcd[157] = 2.0E0*I_ERI_F2xy_Dxz_Dyz_S_a-1*I_ERI_Py_Dxz_Dyz_S;
  abcd[158] = 2.0E0*I_ERI_F2xz_Dxz_Dyz_S_a-1*I_ERI_Pz_Dxz_Dyz_S;
  abcd[159] = 2.0E0*I_ERI_Fx2y_Dxz_Dyz_S_a;
  abcd[160] = 2.0E0*I_ERI_Fxyz_Dxz_Dyz_S_a;
  abcd[161] = 2.0E0*I_ERI_Fx2z_Dxz_Dyz_S_a;
  abcd[162] = 2.0E0*I_ERI_F3x_D2y_Dyz_S_a-2*I_ERI_Px_D2y_Dyz_S;
  abcd[163] = 2.0E0*I_ERI_F2xy_D2y_Dyz_S_a-1*I_ERI_Py_D2y_Dyz_S;
  abcd[164] = 2.0E0*I_ERI_F2xz_D2y_Dyz_S_a-1*I_ERI_Pz_D2y_Dyz_S;
  abcd[165] = 2.0E0*I_ERI_Fx2y_D2y_Dyz_S_a;
  abcd[166] = 2.0E0*I_ERI_Fxyz_D2y_Dyz_S_a;
  abcd[167] = 2.0E0*I_ERI_Fx2z_D2y_Dyz_S_a;
  abcd[168] = 2.0E0*I_ERI_F3x_Dyz_Dyz_S_a-2*I_ERI_Px_Dyz_Dyz_S;
  abcd[169] = 2.0E0*I_ERI_F2xy_Dyz_Dyz_S_a-1*I_ERI_Py_Dyz_Dyz_S;
  abcd[170] = 2.0E0*I_ERI_F2xz_Dyz_Dyz_S_a-1*I_ERI_Pz_Dyz_Dyz_S;
  abcd[171] = 2.0E0*I_ERI_Fx2y_Dyz_Dyz_S_a;
  abcd[172] = 2.0E0*I_ERI_Fxyz_Dyz_Dyz_S_a;
  abcd[173] = 2.0E0*I_ERI_Fx2z_Dyz_Dyz_S_a;
  abcd[174] = 2.0E0*I_ERI_F3x_D2z_Dyz_S_a-2*I_ERI_Px_D2z_Dyz_S;
  abcd[175] = 2.0E0*I_ERI_F2xy_D2z_Dyz_S_a-1*I_ERI_Py_D2z_Dyz_S;
  abcd[176] = 2.0E0*I_ERI_F2xz_D2z_Dyz_S_a-1*I_ERI_Pz_D2z_Dyz_S;
  abcd[177] = 2.0E0*I_ERI_Fx2y_D2z_Dyz_S_a;
  abcd[178] = 2.0E0*I_ERI_Fxyz_D2z_Dyz_S_a;
  abcd[179] = 2.0E0*I_ERI_Fx2z_D2z_Dyz_S_a;
  abcd[180] = 2.0E0*I_ERI_F3x_D2x_D2z_S_a-2*I_ERI_Px_D2x_D2z_S;
  abcd[181] = 2.0E0*I_ERI_F2xy_D2x_D2z_S_a-1*I_ERI_Py_D2x_D2z_S;
  abcd[182] = 2.0E0*I_ERI_F2xz_D2x_D2z_S_a-1*I_ERI_Pz_D2x_D2z_S;
  abcd[183] = 2.0E0*I_ERI_Fx2y_D2x_D2z_S_a;
  abcd[184] = 2.0E0*I_ERI_Fxyz_D2x_D2z_S_a;
  abcd[185] = 2.0E0*I_ERI_Fx2z_D2x_D2z_S_a;
  abcd[186] = 2.0E0*I_ERI_F3x_Dxy_D2z_S_a-2*I_ERI_Px_Dxy_D2z_S;
  abcd[187] = 2.0E0*I_ERI_F2xy_Dxy_D2z_S_a-1*I_ERI_Py_Dxy_D2z_S;
  abcd[188] = 2.0E0*I_ERI_F2xz_Dxy_D2z_S_a-1*I_ERI_Pz_Dxy_D2z_S;
  abcd[189] = 2.0E0*I_ERI_Fx2y_Dxy_D2z_S_a;
  abcd[190] = 2.0E0*I_ERI_Fxyz_Dxy_D2z_S_a;
  abcd[191] = 2.0E0*I_ERI_Fx2z_Dxy_D2z_S_a;
  abcd[192] = 2.0E0*I_ERI_F3x_Dxz_D2z_S_a-2*I_ERI_Px_Dxz_D2z_S;
  abcd[193] = 2.0E0*I_ERI_F2xy_Dxz_D2z_S_a-1*I_ERI_Py_Dxz_D2z_S;
  abcd[194] = 2.0E0*I_ERI_F2xz_Dxz_D2z_S_a-1*I_ERI_Pz_Dxz_D2z_S;
  abcd[195] = 2.0E0*I_ERI_Fx2y_Dxz_D2z_S_a;
  abcd[196] = 2.0E0*I_ERI_Fxyz_Dxz_D2z_S_a;
  abcd[197] = 2.0E0*I_ERI_Fx2z_Dxz_D2z_S_a;
  abcd[198] = 2.0E0*I_ERI_F3x_D2y_D2z_S_a-2*I_ERI_Px_D2y_D2z_S;
  abcd[199] = 2.0E0*I_ERI_F2xy_D2y_D2z_S_a-1*I_ERI_Py_D2y_D2z_S;
  abcd[200] = 2.0E0*I_ERI_F2xz_D2y_D2z_S_a-1*I_ERI_Pz_D2y_D2z_S;
  abcd[201] = 2.0E0*I_ERI_Fx2y_D2y_D2z_S_a;
  abcd[202] = 2.0E0*I_ERI_Fxyz_D2y_D2z_S_a;
  abcd[203] = 2.0E0*I_ERI_Fx2z_D2y_D2z_S_a;
  abcd[204] = 2.0E0*I_ERI_F3x_Dyz_D2z_S_a-2*I_ERI_Px_Dyz_D2z_S;
  abcd[205] = 2.0E0*I_ERI_F2xy_Dyz_D2z_S_a-1*I_ERI_Py_Dyz_D2z_S;
  abcd[206] = 2.0E0*I_ERI_F2xz_Dyz_D2z_S_a-1*I_ERI_Pz_Dyz_D2z_S;
  abcd[207] = 2.0E0*I_ERI_Fx2y_Dyz_D2z_S_a;
  abcd[208] = 2.0E0*I_ERI_Fxyz_Dyz_D2z_S_a;
  abcd[209] = 2.0E0*I_ERI_Fx2z_Dyz_D2z_S_a;
  abcd[210] = 2.0E0*I_ERI_F3x_D2z_D2z_S_a-2*I_ERI_Px_D2z_D2z_S;
  abcd[211] = 2.0E0*I_ERI_F2xy_D2z_D2z_S_a-1*I_ERI_Py_D2z_D2z_S;
  abcd[212] = 2.0E0*I_ERI_F2xz_D2z_D2z_S_a-1*I_ERI_Pz_D2z_D2z_S;
  abcd[213] = 2.0E0*I_ERI_Fx2y_D2z_D2z_S_a;
  abcd[214] = 2.0E0*I_ERI_Fxyz_D2z_D2z_S_a;
  abcd[215] = 2.0E0*I_ERI_Fx2z_D2z_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_a
   * RHS shell quartet name: SQ_ERI_P_D_D_S
   ************************************************************/
  abcd[216] = 2.0E0*I_ERI_F2xy_D2x_D2x_S_a;
  abcd[217] = 2.0E0*I_ERI_Fx2y_D2x_D2x_S_a-1*I_ERI_Px_D2x_D2x_S;
  abcd[218] = 2.0E0*I_ERI_Fxyz_D2x_D2x_S_a;
  abcd[219] = 2.0E0*I_ERI_F3y_D2x_D2x_S_a-2*I_ERI_Py_D2x_D2x_S;
  abcd[220] = 2.0E0*I_ERI_F2yz_D2x_D2x_S_a-1*I_ERI_Pz_D2x_D2x_S;
  abcd[221] = 2.0E0*I_ERI_Fy2z_D2x_D2x_S_a;
  abcd[222] = 2.0E0*I_ERI_F2xy_Dxy_D2x_S_a;
  abcd[223] = 2.0E0*I_ERI_Fx2y_Dxy_D2x_S_a-1*I_ERI_Px_Dxy_D2x_S;
  abcd[224] = 2.0E0*I_ERI_Fxyz_Dxy_D2x_S_a;
  abcd[225] = 2.0E0*I_ERI_F3y_Dxy_D2x_S_a-2*I_ERI_Py_Dxy_D2x_S;
  abcd[226] = 2.0E0*I_ERI_F2yz_Dxy_D2x_S_a-1*I_ERI_Pz_Dxy_D2x_S;
  abcd[227] = 2.0E0*I_ERI_Fy2z_Dxy_D2x_S_a;
  abcd[228] = 2.0E0*I_ERI_F2xy_Dxz_D2x_S_a;
  abcd[229] = 2.0E0*I_ERI_Fx2y_Dxz_D2x_S_a-1*I_ERI_Px_Dxz_D2x_S;
  abcd[230] = 2.0E0*I_ERI_Fxyz_Dxz_D2x_S_a;
  abcd[231] = 2.0E0*I_ERI_F3y_Dxz_D2x_S_a-2*I_ERI_Py_Dxz_D2x_S;
  abcd[232] = 2.0E0*I_ERI_F2yz_Dxz_D2x_S_a-1*I_ERI_Pz_Dxz_D2x_S;
  abcd[233] = 2.0E0*I_ERI_Fy2z_Dxz_D2x_S_a;
  abcd[234] = 2.0E0*I_ERI_F2xy_D2y_D2x_S_a;
  abcd[235] = 2.0E0*I_ERI_Fx2y_D2y_D2x_S_a-1*I_ERI_Px_D2y_D2x_S;
  abcd[236] = 2.0E0*I_ERI_Fxyz_D2y_D2x_S_a;
  abcd[237] = 2.0E0*I_ERI_F3y_D2y_D2x_S_a-2*I_ERI_Py_D2y_D2x_S;
  abcd[238] = 2.0E0*I_ERI_F2yz_D2y_D2x_S_a-1*I_ERI_Pz_D2y_D2x_S;
  abcd[239] = 2.0E0*I_ERI_Fy2z_D2y_D2x_S_a;
  abcd[240] = 2.0E0*I_ERI_F2xy_Dyz_D2x_S_a;
  abcd[241] = 2.0E0*I_ERI_Fx2y_Dyz_D2x_S_a-1*I_ERI_Px_Dyz_D2x_S;
  abcd[242] = 2.0E0*I_ERI_Fxyz_Dyz_D2x_S_a;
  abcd[243] = 2.0E0*I_ERI_F3y_Dyz_D2x_S_a-2*I_ERI_Py_Dyz_D2x_S;
  abcd[244] = 2.0E0*I_ERI_F2yz_Dyz_D2x_S_a-1*I_ERI_Pz_Dyz_D2x_S;
  abcd[245] = 2.0E0*I_ERI_Fy2z_Dyz_D2x_S_a;
  abcd[246] = 2.0E0*I_ERI_F2xy_D2z_D2x_S_a;
  abcd[247] = 2.0E0*I_ERI_Fx2y_D2z_D2x_S_a-1*I_ERI_Px_D2z_D2x_S;
  abcd[248] = 2.0E0*I_ERI_Fxyz_D2z_D2x_S_a;
  abcd[249] = 2.0E0*I_ERI_F3y_D2z_D2x_S_a-2*I_ERI_Py_D2z_D2x_S;
  abcd[250] = 2.0E0*I_ERI_F2yz_D2z_D2x_S_a-1*I_ERI_Pz_D2z_D2x_S;
  abcd[251] = 2.0E0*I_ERI_Fy2z_D2z_D2x_S_a;
  abcd[252] = 2.0E0*I_ERI_F2xy_D2x_Dxy_S_a;
  abcd[253] = 2.0E0*I_ERI_Fx2y_D2x_Dxy_S_a-1*I_ERI_Px_D2x_Dxy_S;
  abcd[254] = 2.0E0*I_ERI_Fxyz_D2x_Dxy_S_a;
  abcd[255] = 2.0E0*I_ERI_F3y_D2x_Dxy_S_a-2*I_ERI_Py_D2x_Dxy_S;
  abcd[256] = 2.0E0*I_ERI_F2yz_D2x_Dxy_S_a-1*I_ERI_Pz_D2x_Dxy_S;
  abcd[257] = 2.0E0*I_ERI_Fy2z_D2x_Dxy_S_a;
  abcd[258] = 2.0E0*I_ERI_F2xy_Dxy_Dxy_S_a;
  abcd[259] = 2.0E0*I_ERI_Fx2y_Dxy_Dxy_S_a-1*I_ERI_Px_Dxy_Dxy_S;
  abcd[260] = 2.0E0*I_ERI_Fxyz_Dxy_Dxy_S_a;
  abcd[261] = 2.0E0*I_ERI_F3y_Dxy_Dxy_S_a-2*I_ERI_Py_Dxy_Dxy_S;
  abcd[262] = 2.0E0*I_ERI_F2yz_Dxy_Dxy_S_a-1*I_ERI_Pz_Dxy_Dxy_S;
  abcd[263] = 2.0E0*I_ERI_Fy2z_Dxy_Dxy_S_a;
  abcd[264] = 2.0E0*I_ERI_F2xy_Dxz_Dxy_S_a;
  abcd[265] = 2.0E0*I_ERI_Fx2y_Dxz_Dxy_S_a-1*I_ERI_Px_Dxz_Dxy_S;
  abcd[266] = 2.0E0*I_ERI_Fxyz_Dxz_Dxy_S_a;
  abcd[267] = 2.0E0*I_ERI_F3y_Dxz_Dxy_S_a-2*I_ERI_Py_Dxz_Dxy_S;
  abcd[268] = 2.0E0*I_ERI_F2yz_Dxz_Dxy_S_a-1*I_ERI_Pz_Dxz_Dxy_S;
  abcd[269] = 2.0E0*I_ERI_Fy2z_Dxz_Dxy_S_a;
  abcd[270] = 2.0E0*I_ERI_F2xy_D2y_Dxy_S_a;
  abcd[271] = 2.0E0*I_ERI_Fx2y_D2y_Dxy_S_a-1*I_ERI_Px_D2y_Dxy_S;
  abcd[272] = 2.0E0*I_ERI_Fxyz_D2y_Dxy_S_a;
  abcd[273] = 2.0E0*I_ERI_F3y_D2y_Dxy_S_a-2*I_ERI_Py_D2y_Dxy_S;
  abcd[274] = 2.0E0*I_ERI_F2yz_D2y_Dxy_S_a-1*I_ERI_Pz_D2y_Dxy_S;
  abcd[275] = 2.0E0*I_ERI_Fy2z_D2y_Dxy_S_a;
  abcd[276] = 2.0E0*I_ERI_F2xy_Dyz_Dxy_S_a;
  abcd[277] = 2.0E0*I_ERI_Fx2y_Dyz_Dxy_S_a-1*I_ERI_Px_Dyz_Dxy_S;
  abcd[278] = 2.0E0*I_ERI_Fxyz_Dyz_Dxy_S_a;
  abcd[279] = 2.0E0*I_ERI_F3y_Dyz_Dxy_S_a-2*I_ERI_Py_Dyz_Dxy_S;
  abcd[280] = 2.0E0*I_ERI_F2yz_Dyz_Dxy_S_a-1*I_ERI_Pz_Dyz_Dxy_S;
  abcd[281] = 2.0E0*I_ERI_Fy2z_Dyz_Dxy_S_a;
  abcd[282] = 2.0E0*I_ERI_F2xy_D2z_Dxy_S_a;
  abcd[283] = 2.0E0*I_ERI_Fx2y_D2z_Dxy_S_a-1*I_ERI_Px_D2z_Dxy_S;
  abcd[284] = 2.0E0*I_ERI_Fxyz_D2z_Dxy_S_a;
  abcd[285] = 2.0E0*I_ERI_F3y_D2z_Dxy_S_a-2*I_ERI_Py_D2z_Dxy_S;
  abcd[286] = 2.0E0*I_ERI_F2yz_D2z_Dxy_S_a-1*I_ERI_Pz_D2z_Dxy_S;
  abcd[287] = 2.0E0*I_ERI_Fy2z_D2z_Dxy_S_a;
  abcd[288] = 2.0E0*I_ERI_F2xy_D2x_Dxz_S_a;
  abcd[289] = 2.0E0*I_ERI_Fx2y_D2x_Dxz_S_a-1*I_ERI_Px_D2x_Dxz_S;
  abcd[290] = 2.0E0*I_ERI_Fxyz_D2x_Dxz_S_a;
  abcd[291] = 2.0E0*I_ERI_F3y_D2x_Dxz_S_a-2*I_ERI_Py_D2x_Dxz_S;
  abcd[292] = 2.0E0*I_ERI_F2yz_D2x_Dxz_S_a-1*I_ERI_Pz_D2x_Dxz_S;
  abcd[293] = 2.0E0*I_ERI_Fy2z_D2x_Dxz_S_a;
  abcd[294] = 2.0E0*I_ERI_F2xy_Dxy_Dxz_S_a;
  abcd[295] = 2.0E0*I_ERI_Fx2y_Dxy_Dxz_S_a-1*I_ERI_Px_Dxy_Dxz_S;
  abcd[296] = 2.0E0*I_ERI_Fxyz_Dxy_Dxz_S_a;
  abcd[297] = 2.0E0*I_ERI_F3y_Dxy_Dxz_S_a-2*I_ERI_Py_Dxy_Dxz_S;
  abcd[298] = 2.0E0*I_ERI_F2yz_Dxy_Dxz_S_a-1*I_ERI_Pz_Dxy_Dxz_S;
  abcd[299] = 2.0E0*I_ERI_Fy2z_Dxy_Dxz_S_a;
  abcd[300] = 2.0E0*I_ERI_F2xy_Dxz_Dxz_S_a;
  abcd[301] = 2.0E0*I_ERI_Fx2y_Dxz_Dxz_S_a-1*I_ERI_Px_Dxz_Dxz_S;
  abcd[302] = 2.0E0*I_ERI_Fxyz_Dxz_Dxz_S_a;
  abcd[303] = 2.0E0*I_ERI_F3y_Dxz_Dxz_S_a-2*I_ERI_Py_Dxz_Dxz_S;
  abcd[304] = 2.0E0*I_ERI_F2yz_Dxz_Dxz_S_a-1*I_ERI_Pz_Dxz_Dxz_S;
  abcd[305] = 2.0E0*I_ERI_Fy2z_Dxz_Dxz_S_a;
  abcd[306] = 2.0E0*I_ERI_F2xy_D2y_Dxz_S_a;
  abcd[307] = 2.0E0*I_ERI_Fx2y_D2y_Dxz_S_a-1*I_ERI_Px_D2y_Dxz_S;
  abcd[308] = 2.0E0*I_ERI_Fxyz_D2y_Dxz_S_a;
  abcd[309] = 2.0E0*I_ERI_F3y_D2y_Dxz_S_a-2*I_ERI_Py_D2y_Dxz_S;
  abcd[310] = 2.0E0*I_ERI_F2yz_D2y_Dxz_S_a-1*I_ERI_Pz_D2y_Dxz_S;
  abcd[311] = 2.0E0*I_ERI_Fy2z_D2y_Dxz_S_a;
  abcd[312] = 2.0E0*I_ERI_F2xy_Dyz_Dxz_S_a;
  abcd[313] = 2.0E0*I_ERI_Fx2y_Dyz_Dxz_S_a-1*I_ERI_Px_Dyz_Dxz_S;
  abcd[314] = 2.0E0*I_ERI_Fxyz_Dyz_Dxz_S_a;
  abcd[315] = 2.0E0*I_ERI_F3y_Dyz_Dxz_S_a-2*I_ERI_Py_Dyz_Dxz_S;
  abcd[316] = 2.0E0*I_ERI_F2yz_Dyz_Dxz_S_a-1*I_ERI_Pz_Dyz_Dxz_S;
  abcd[317] = 2.0E0*I_ERI_Fy2z_Dyz_Dxz_S_a;
  abcd[318] = 2.0E0*I_ERI_F2xy_D2z_Dxz_S_a;
  abcd[319] = 2.0E0*I_ERI_Fx2y_D2z_Dxz_S_a-1*I_ERI_Px_D2z_Dxz_S;
  abcd[320] = 2.0E0*I_ERI_Fxyz_D2z_Dxz_S_a;
  abcd[321] = 2.0E0*I_ERI_F3y_D2z_Dxz_S_a-2*I_ERI_Py_D2z_Dxz_S;
  abcd[322] = 2.0E0*I_ERI_F2yz_D2z_Dxz_S_a-1*I_ERI_Pz_D2z_Dxz_S;
  abcd[323] = 2.0E0*I_ERI_Fy2z_D2z_Dxz_S_a;
  abcd[324] = 2.0E0*I_ERI_F2xy_D2x_D2y_S_a;
  abcd[325] = 2.0E0*I_ERI_Fx2y_D2x_D2y_S_a-1*I_ERI_Px_D2x_D2y_S;
  abcd[326] = 2.0E0*I_ERI_Fxyz_D2x_D2y_S_a;
  abcd[327] = 2.0E0*I_ERI_F3y_D2x_D2y_S_a-2*I_ERI_Py_D2x_D2y_S;
  abcd[328] = 2.0E0*I_ERI_F2yz_D2x_D2y_S_a-1*I_ERI_Pz_D2x_D2y_S;
  abcd[329] = 2.0E0*I_ERI_Fy2z_D2x_D2y_S_a;
  abcd[330] = 2.0E0*I_ERI_F2xy_Dxy_D2y_S_a;
  abcd[331] = 2.0E0*I_ERI_Fx2y_Dxy_D2y_S_a-1*I_ERI_Px_Dxy_D2y_S;
  abcd[332] = 2.0E0*I_ERI_Fxyz_Dxy_D2y_S_a;
  abcd[333] = 2.0E0*I_ERI_F3y_Dxy_D2y_S_a-2*I_ERI_Py_Dxy_D2y_S;
  abcd[334] = 2.0E0*I_ERI_F2yz_Dxy_D2y_S_a-1*I_ERI_Pz_Dxy_D2y_S;
  abcd[335] = 2.0E0*I_ERI_Fy2z_Dxy_D2y_S_a;
  abcd[336] = 2.0E0*I_ERI_F2xy_Dxz_D2y_S_a;
  abcd[337] = 2.0E0*I_ERI_Fx2y_Dxz_D2y_S_a-1*I_ERI_Px_Dxz_D2y_S;
  abcd[338] = 2.0E0*I_ERI_Fxyz_Dxz_D2y_S_a;
  abcd[339] = 2.0E0*I_ERI_F3y_Dxz_D2y_S_a-2*I_ERI_Py_Dxz_D2y_S;
  abcd[340] = 2.0E0*I_ERI_F2yz_Dxz_D2y_S_a-1*I_ERI_Pz_Dxz_D2y_S;
  abcd[341] = 2.0E0*I_ERI_Fy2z_Dxz_D2y_S_a;
  abcd[342] = 2.0E0*I_ERI_F2xy_D2y_D2y_S_a;
  abcd[343] = 2.0E0*I_ERI_Fx2y_D2y_D2y_S_a-1*I_ERI_Px_D2y_D2y_S;
  abcd[344] = 2.0E0*I_ERI_Fxyz_D2y_D2y_S_a;
  abcd[345] = 2.0E0*I_ERI_F3y_D2y_D2y_S_a-2*I_ERI_Py_D2y_D2y_S;
  abcd[346] = 2.0E0*I_ERI_F2yz_D2y_D2y_S_a-1*I_ERI_Pz_D2y_D2y_S;
  abcd[347] = 2.0E0*I_ERI_Fy2z_D2y_D2y_S_a;
  abcd[348] = 2.0E0*I_ERI_F2xy_Dyz_D2y_S_a;
  abcd[349] = 2.0E0*I_ERI_Fx2y_Dyz_D2y_S_a-1*I_ERI_Px_Dyz_D2y_S;
  abcd[350] = 2.0E0*I_ERI_Fxyz_Dyz_D2y_S_a;
  abcd[351] = 2.0E0*I_ERI_F3y_Dyz_D2y_S_a-2*I_ERI_Py_Dyz_D2y_S;
  abcd[352] = 2.0E0*I_ERI_F2yz_Dyz_D2y_S_a-1*I_ERI_Pz_Dyz_D2y_S;
  abcd[353] = 2.0E0*I_ERI_Fy2z_Dyz_D2y_S_a;
  abcd[354] = 2.0E0*I_ERI_F2xy_D2z_D2y_S_a;
  abcd[355] = 2.0E0*I_ERI_Fx2y_D2z_D2y_S_a-1*I_ERI_Px_D2z_D2y_S;
  abcd[356] = 2.0E0*I_ERI_Fxyz_D2z_D2y_S_a;
  abcd[357] = 2.0E0*I_ERI_F3y_D2z_D2y_S_a-2*I_ERI_Py_D2z_D2y_S;
  abcd[358] = 2.0E0*I_ERI_F2yz_D2z_D2y_S_a-1*I_ERI_Pz_D2z_D2y_S;
  abcd[359] = 2.0E0*I_ERI_Fy2z_D2z_D2y_S_a;
  abcd[360] = 2.0E0*I_ERI_F2xy_D2x_Dyz_S_a;
  abcd[361] = 2.0E0*I_ERI_Fx2y_D2x_Dyz_S_a-1*I_ERI_Px_D2x_Dyz_S;
  abcd[362] = 2.0E0*I_ERI_Fxyz_D2x_Dyz_S_a;
  abcd[363] = 2.0E0*I_ERI_F3y_D2x_Dyz_S_a-2*I_ERI_Py_D2x_Dyz_S;
  abcd[364] = 2.0E0*I_ERI_F2yz_D2x_Dyz_S_a-1*I_ERI_Pz_D2x_Dyz_S;
  abcd[365] = 2.0E0*I_ERI_Fy2z_D2x_Dyz_S_a;
  abcd[366] = 2.0E0*I_ERI_F2xy_Dxy_Dyz_S_a;
  abcd[367] = 2.0E0*I_ERI_Fx2y_Dxy_Dyz_S_a-1*I_ERI_Px_Dxy_Dyz_S;
  abcd[368] = 2.0E0*I_ERI_Fxyz_Dxy_Dyz_S_a;
  abcd[369] = 2.0E0*I_ERI_F3y_Dxy_Dyz_S_a-2*I_ERI_Py_Dxy_Dyz_S;
  abcd[370] = 2.0E0*I_ERI_F2yz_Dxy_Dyz_S_a-1*I_ERI_Pz_Dxy_Dyz_S;
  abcd[371] = 2.0E0*I_ERI_Fy2z_Dxy_Dyz_S_a;
  abcd[372] = 2.0E0*I_ERI_F2xy_Dxz_Dyz_S_a;
  abcd[373] = 2.0E0*I_ERI_Fx2y_Dxz_Dyz_S_a-1*I_ERI_Px_Dxz_Dyz_S;
  abcd[374] = 2.0E0*I_ERI_Fxyz_Dxz_Dyz_S_a;
  abcd[375] = 2.0E0*I_ERI_F3y_Dxz_Dyz_S_a-2*I_ERI_Py_Dxz_Dyz_S;
  abcd[376] = 2.0E0*I_ERI_F2yz_Dxz_Dyz_S_a-1*I_ERI_Pz_Dxz_Dyz_S;
  abcd[377] = 2.0E0*I_ERI_Fy2z_Dxz_Dyz_S_a;
  abcd[378] = 2.0E0*I_ERI_F2xy_D2y_Dyz_S_a;
  abcd[379] = 2.0E0*I_ERI_Fx2y_D2y_Dyz_S_a-1*I_ERI_Px_D2y_Dyz_S;
  abcd[380] = 2.0E0*I_ERI_Fxyz_D2y_Dyz_S_a;
  abcd[381] = 2.0E0*I_ERI_F3y_D2y_Dyz_S_a-2*I_ERI_Py_D2y_Dyz_S;
  abcd[382] = 2.0E0*I_ERI_F2yz_D2y_Dyz_S_a-1*I_ERI_Pz_D2y_Dyz_S;
  abcd[383] = 2.0E0*I_ERI_Fy2z_D2y_Dyz_S_a;
  abcd[384] = 2.0E0*I_ERI_F2xy_Dyz_Dyz_S_a;
  abcd[385] = 2.0E0*I_ERI_Fx2y_Dyz_Dyz_S_a-1*I_ERI_Px_Dyz_Dyz_S;
  abcd[386] = 2.0E0*I_ERI_Fxyz_Dyz_Dyz_S_a;
  abcd[387] = 2.0E0*I_ERI_F3y_Dyz_Dyz_S_a-2*I_ERI_Py_Dyz_Dyz_S;
  abcd[388] = 2.0E0*I_ERI_F2yz_Dyz_Dyz_S_a-1*I_ERI_Pz_Dyz_Dyz_S;
  abcd[389] = 2.0E0*I_ERI_Fy2z_Dyz_Dyz_S_a;
  abcd[390] = 2.0E0*I_ERI_F2xy_D2z_Dyz_S_a;
  abcd[391] = 2.0E0*I_ERI_Fx2y_D2z_Dyz_S_a-1*I_ERI_Px_D2z_Dyz_S;
  abcd[392] = 2.0E0*I_ERI_Fxyz_D2z_Dyz_S_a;
  abcd[393] = 2.0E0*I_ERI_F3y_D2z_Dyz_S_a-2*I_ERI_Py_D2z_Dyz_S;
  abcd[394] = 2.0E0*I_ERI_F2yz_D2z_Dyz_S_a-1*I_ERI_Pz_D2z_Dyz_S;
  abcd[395] = 2.0E0*I_ERI_Fy2z_D2z_Dyz_S_a;
  abcd[396] = 2.0E0*I_ERI_F2xy_D2x_D2z_S_a;
  abcd[397] = 2.0E0*I_ERI_Fx2y_D2x_D2z_S_a-1*I_ERI_Px_D2x_D2z_S;
  abcd[398] = 2.0E0*I_ERI_Fxyz_D2x_D2z_S_a;
  abcd[399] = 2.0E0*I_ERI_F3y_D2x_D2z_S_a-2*I_ERI_Py_D2x_D2z_S;
  abcd[400] = 2.0E0*I_ERI_F2yz_D2x_D2z_S_a-1*I_ERI_Pz_D2x_D2z_S;
  abcd[401] = 2.0E0*I_ERI_Fy2z_D2x_D2z_S_a;
  abcd[402] = 2.0E0*I_ERI_F2xy_Dxy_D2z_S_a;
  abcd[403] = 2.0E0*I_ERI_Fx2y_Dxy_D2z_S_a-1*I_ERI_Px_Dxy_D2z_S;
  abcd[404] = 2.0E0*I_ERI_Fxyz_Dxy_D2z_S_a;
  abcd[405] = 2.0E0*I_ERI_F3y_Dxy_D2z_S_a-2*I_ERI_Py_Dxy_D2z_S;
  abcd[406] = 2.0E0*I_ERI_F2yz_Dxy_D2z_S_a-1*I_ERI_Pz_Dxy_D2z_S;
  abcd[407] = 2.0E0*I_ERI_Fy2z_Dxy_D2z_S_a;
  abcd[408] = 2.0E0*I_ERI_F2xy_Dxz_D2z_S_a;
  abcd[409] = 2.0E0*I_ERI_Fx2y_Dxz_D2z_S_a-1*I_ERI_Px_Dxz_D2z_S;
  abcd[410] = 2.0E0*I_ERI_Fxyz_Dxz_D2z_S_a;
  abcd[411] = 2.0E0*I_ERI_F3y_Dxz_D2z_S_a-2*I_ERI_Py_Dxz_D2z_S;
  abcd[412] = 2.0E0*I_ERI_F2yz_Dxz_D2z_S_a-1*I_ERI_Pz_Dxz_D2z_S;
  abcd[413] = 2.0E0*I_ERI_Fy2z_Dxz_D2z_S_a;
  abcd[414] = 2.0E0*I_ERI_F2xy_D2y_D2z_S_a;
  abcd[415] = 2.0E0*I_ERI_Fx2y_D2y_D2z_S_a-1*I_ERI_Px_D2y_D2z_S;
  abcd[416] = 2.0E0*I_ERI_Fxyz_D2y_D2z_S_a;
  abcd[417] = 2.0E0*I_ERI_F3y_D2y_D2z_S_a-2*I_ERI_Py_D2y_D2z_S;
  abcd[418] = 2.0E0*I_ERI_F2yz_D2y_D2z_S_a-1*I_ERI_Pz_D2y_D2z_S;
  abcd[419] = 2.0E0*I_ERI_Fy2z_D2y_D2z_S_a;
  abcd[420] = 2.0E0*I_ERI_F2xy_Dyz_D2z_S_a;
  abcd[421] = 2.0E0*I_ERI_Fx2y_Dyz_D2z_S_a-1*I_ERI_Px_Dyz_D2z_S;
  abcd[422] = 2.0E0*I_ERI_Fxyz_Dyz_D2z_S_a;
  abcd[423] = 2.0E0*I_ERI_F3y_Dyz_D2z_S_a-2*I_ERI_Py_Dyz_D2z_S;
  abcd[424] = 2.0E0*I_ERI_F2yz_Dyz_D2z_S_a-1*I_ERI_Pz_Dyz_D2z_S;
  abcd[425] = 2.0E0*I_ERI_Fy2z_Dyz_D2z_S_a;
  abcd[426] = 2.0E0*I_ERI_F2xy_D2z_D2z_S_a;
  abcd[427] = 2.0E0*I_ERI_Fx2y_D2z_D2z_S_a-1*I_ERI_Px_D2z_D2z_S;
  abcd[428] = 2.0E0*I_ERI_Fxyz_D2z_D2z_S_a;
  abcd[429] = 2.0E0*I_ERI_F3y_D2z_D2z_S_a-2*I_ERI_Py_D2z_D2z_S;
  abcd[430] = 2.0E0*I_ERI_F2yz_D2z_D2z_S_a-1*I_ERI_Pz_D2z_D2z_S;
  abcd[431] = 2.0E0*I_ERI_Fy2z_D2z_D2z_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_a
   * RHS shell quartet name: SQ_ERI_P_D_D_S
   ************************************************************/
  abcd[432] = 2.0E0*I_ERI_F2xz_D2x_D2x_S_a;
  abcd[433] = 2.0E0*I_ERI_Fxyz_D2x_D2x_S_a;
  abcd[434] = 2.0E0*I_ERI_Fx2z_D2x_D2x_S_a-1*I_ERI_Px_D2x_D2x_S;
  abcd[435] = 2.0E0*I_ERI_F2yz_D2x_D2x_S_a;
  abcd[436] = 2.0E0*I_ERI_Fy2z_D2x_D2x_S_a-1*I_ERI_Py_D2x_D2x_S;
  abcd[437] = 2.0E0*I_ERI_F3z_D2x_D2x_S_a-2*I_ERI_Pz_D2x_D2x_S;
  abcd[438] = 2.0E0*I_ERI_F2xz_Dxy_D2x_S_a;
  abcd[439] = 2.0E0*I_ERI_Fxyz_Dxy_D2x_S_a;
  abcd[440] = 2.0E0*I_ERI_Fx2z_Dxy_D2x_S_a-1*I_ERI_Px_Dxy_D2x_S;
  abcd[441] = 2.0E0*I_ERI_F2yz_Dxy_D2x_S_a;
  abcd[442] = 2.0E0*I_ERI_Fy2z_Dxy_D2x_S_a-1*I_ERI_Py_Dxy_D2x_S;
  abcd[443] = 2.0E0*I_ERI_F3z_Dxy_D2x_S_a-2*I_ERI_Pz_Dxy_D2x_S;
  abcd[444] = 2.0E0*I_ERI_F2xz_Dxz_D2x_S_a;
  abcd[445] = 2.0E0*I_ERI_Fxyz_Dxz_D2x_S_a;
  abcd[446] = 2.0E0*I_ERI_Fx2z_Dxz_D2x_S_a-1*I_ERI_Px_Dxz_D2x_S;
  abcd[447] = 2.0E0*I_ERI_F2yz_Dxz_D2x_S_a;
  abcd[448] = 2.0E0*I_ERI_Fy2z_Dxz_D2x_S_a-1*I_ERI_Py_Dxz_D2x_S;
  abcd[449] = 2.0E0*I_ERI_F3z_Dxz_D2x_S_a-2*I_ERI_Pz_Dxz_D2x_S;
  abcd[450] = 2.0E0*I_ERI_F2xz_D2y_D2x_S_a;
  abcd[451] = 2.0E0*I_ERI_Fxyz_D2y_D2x_S_a;
  abcd[452] = 2.0E0*I_ERI_Fx2z_D2y_D2x_S_a-1*I_ERI_Px_D2y_D2x_S;
  abcd[453] = 2.0E0*I_ERI_F2yz_D2y_D2x_S_a;
  abcd[454] = 2.0E0*I_ERI_Fy2z_D2y_D2x_S_a-1*I_ERI_Py_D2y_D2x_S;
  abcd[455] = 2.0E0*I_ERI_F3z_D2y_D2x_S_a-2*I_ERI_Pz_D2y_D2x_S;
  abcd[456] = 2.0E0*I_ERI_F2xz_Dyz_D2x_S_a;
  abcd[457] = 2.0E0*I_ERI_Fxyz_Dyz_D2x_S_a;
  abcd[458] = 2.0E0*I_ERI_Fx2z_Dyz_D2x_S_a-1*I_ERI_Px_Dyz_D2x_S;
  abcd[459] = 2.0E0*I_ERI_F2yz_Dyz_D2x_S_a;
  abcd[460] = 2.0E0*I_ERI_Fy2z_Dyz_D2x_S_a-1*I_ERI_Py_Dyz_D2x_S;
  abcd[461] = 2.0E0*I_ERI_F3z_Dyz_D2x_S_a-2*I_ERI_Pz_Dyz_D2x_S;
  abcd[462] = 2.0E0*I_ERI_F2xz_D2z_D2x_S_a;
  abcd[463] = 2.0E0*I_ERI_Fxyz_D2z_D2x_S_a;
  abcd[464] = 2.0E0*I_ERI_Fx2z_D2z_D2x_S_a-1*I_ERI_Px_D2z_D2x_S;
  abcd[465] = 2.0E0*I_ERI_F2yz_D2z_D2x_S_a;
  abcd[466] = 2.0E0*I_ERI_Fy2z_D2z_D2x_S_a-1*I_ERI_Py_D2z_D2x_S;
  abcd[467] = 2.0E0*I_ERI_F3z_D2z_D2x_S_a-2*I_ERI_Pz_D2z_D2x_S;
  abcd[468] = 2.0E0*I_ERI_F2xz_D2x_Dxy_S_a;
  abcd[469] = 2.0E0*I_ERI_Fxyz_D2x_Dxy_S_a;
  abcd[470] = 2.0E0*I_ERI_Fx2z_D2x_Dxy_S_a-1*I_ERI_Px_D2x_Dxy_S;
  abcd[471] = 2.0E0*I_ERI_F2yz_D2x_Dxy_S_a;
  abcd[472] = 2.0E0*I_ERI_Fy2z_D2x_Dxy_S_a-1*I_ERI_Py_D2x_Dxy_S;
  abcd[473] = 2.0E0*I_ERI_F3z_D2x_Dxy_S_a-2*I_ERI_Pz_D2x_Dxy_S;
  abcd[474] = 2.0E0*I_ERI_F2xz_Dxy_Dxy_S_a;
  abcd[475] = 2.0E0*I_ERI_Fxyz_Dxy_Dxy_S_a;
  abcd[476] = 2.0E0*I_ERI_Fx2z_Dxy_Dxy_S_a-1*I_ERI_Px_Dxy_Dxy_S;
  abcd[477] = 2.0E0*I_ERI_F2yz_Dxy_Dxy_S_a;
  abcd[478] = 2.0E0*I_ERI_Fy2z_Dxy_Dxy_S_a-1*I_ERI_Py_Dxy_Dxy_S;
  abcd[479] = 2.0E0*I_ERI_F3z_Dxy_Dxy_S_a-2*I_ERI_Pz_Dxy_Dxy_S;
  abcd[480] = 2.0E0*I_ERI_F2xz_Dxz_Dxy_S_a;
  abcd[481] = 2.0E0*I_ERI_Fxyz_Dxz_Dxy_S_a;
  abcd[482] = 2.0E0*I_ERI_Fx2z_Dxz_Dxy_S_a-1*I_ERI_Px_Dxz_Dxy_S;
  abcd[483] = 2.0E0*I_ERI_F2yz_Dxz_Dxy_S_a;
  abcd[484] = 2.0E0*I_ERI_Fy2z_Dxz_Dxy_S_a-1*I_ERI_Py_Dxz_Dxy_S;
  abcd[485] = 2.0E0*I_ERI_F3z_Dxz_Dxy_S_a-2*I_ERI_Pz_Dxz_Dxy_S;
  abcd[486] = 2.0E0*I_ERI_F2xz_D2y_Dxy_S_a;
  abcd[487] = 2.0E0*I_ERI_Fxyz_D2y_Dxy_S_a;
  abcd[488] = 2.0E0*I_ERI_Fx2z_D2y_Dxy_S_a-1*I_ERI_Px_D2y_Dxy_S;
  abcd[489] = 2.0E0*I_ERI_F2yz_D2y_Dxy_S_a;
  abcd[490] = 2.0E0*I_ERI_Fy2z_D2y_Dxy_S_a-1*I_ERI_Py_D2y_Dxy_S;
  abcd[491] = 2.0E0*I_ERI_F3z_D2y_Dxy_S_a-2*I_ERI_Pz_D2y_Dxy_S;
  abcd[492] = 2.0E0*I_ERI_F2xz_Dyz_Dxy_S_a;
  abcd[493] = 2.0E0*I_ERI_Fxyz_Dyz_Dxy_S_a;
  abcd[494] = 2.0E0*I_ERI_Fx2z_Dyz_Dxy_S_a-1*I_ERI_Px_Dyz_Dxy_S;
  abcd[495] = 2.0E0*I_ERI_F2yz_Dyz_Dxy_S_a;
  abcd[496] = 2.0E0*I_ERI_Fy2z_Dyz_Dxy_S_a-1*I_ERI_Py_Dyz_Dxy_S;
  abcd[497] = 2.0E0*I_ERI_F3z_Dyz_Dxy_S_a-2*I_ERI_Pz_Dyz_Dxy_S;
  abcd[498] = 2.0E0*I_ERI_F2xz_D2z_Dxy_S_a;
  abcd[499] = 2.0E0*I_ERI_Fxyz_D2z_Dxy_S_a;
  abcd[500] = 2.0E0*I_ERI_Fx2z_D2z_Dxy_S_a-1*I_ERI_Px_D2z_Dxy_S;
  abcd[501] = 2.0E0*I_ERI_F2yz_D2z_Dxy_S_a;
  abcd[502] = 2.0E0*I_ERI_Fy2z_D2z_Dxy_S_a-1*I_ERI_Py_D2z_Dxy_S;
  abcd[503] = 2.0E0*I_ERI_F3z_D2z_Dxy_S_a-2*I_ERI_Pz_D2z_Dxy_S;
  abcd[504] = 2.0E0*I_ERI_F2xz_D2x_Dxz_S_a;
  abcd[505] = 2.0E0*I_ERI_Fxyz_D2x_Dxz_S_a;
  abcd[506] = 2.0E0*I_ERI_Fx2z_D2x_Dxz_S_a-1*I_ERI_Px_D2x_Dxz_S;
  abcd[507] = 2.0E0*I_ERI_F2yz_D2x_Dxz_S_a;
  abcd[508] = 2.0E0*I_ERI_Fy2z_D2x_Dxz_S_a-1*I_ERI_Py_D2x_Dxz_S;
  abcd[509] = 2.0E0*I_ERI_F3z_D2x_Dxz_S_a-2*I_ERI_Pz_D2x_Dxz_S;
  abcd[510] = 2.0E0*I_ERI_F2xz_Dxy_Dxz_S_a;
  abcd[511] = 2.0E0*I_ERI_Fxyz_Dxy_Dxz_S_a;
  abcd[512] = 2.0E0*I_ERI_Fx2z_Dxy_Dxz_S_a-1*I_ERI_Px_Dxy_Dxz_S;
  abcd[513] = 2.0E0*I_ERI_F2yz_Dxy_Dxz_S_a;
  abcd[514] = 2.0E0*I_ERI_Fy2z_Dxy_Dxz_S_a-1*I_ERI_Py_Dxy_Dxz_S;
  abcd[515] = 2.0E0*I_ERI_F3z_Dxy_Dxz_S_a-2*I_ERI_Pz_Dxy_Dxz_S;
  abcd[516] = 2.0E0*I_ERI_F2xz_Dxz_Dxz_S_a;
  abcd[517] = 2.0E0*I_ERI_Fxyz_Dxz_Dxz_S_a;
  abcd[518] = 2.0E0*I_ERI_Fx2z_Dxz_Dxz_S_a-1*I_ERI_Px_Dxz_Dxz_S;
  abcd[519] = 2.0E0*I_ERI_F2yz_Dxz_Dxz_S_a;
  abcd[520] = 2.0E0*I_ERI_Fy2z_Dxz_Dxz_S_a-1*I_ERI_Py_Dxz_Dxz_S;
  abcd[521] = 2.0E0*I_ERI_F3z_Dxz_Dxz_S_a-2*I_ERI_Pz_Dxz_Dxz_S;
  abcd[522] = 2.0E0*I_ERI_F2xz_D2y_Dxz_S_a;
  abcd[523] = 2.0E0*I_ERI_Fxyz_D2y_Dxz_S_a;
  abcd[524] = 2.0E0*I_ERI_Fx2z_D2y_Dxz_S_a-1*I_ERI_Px_D2y_Dxz_S;
  abcd[525] = 2.0E0*I_ERI_F2yz_D2y_Dxz_S_a;
  abcd[526] = 2.0E0*I_ERI_Fy2z_D2y_Dxz_S_a-1*I_ERI_Py_D2y_Dxz_S;
  abcd[527] = 2.0E0*I_ERI_F3z_D2y_Dxz_S_a-2*I_ERI_Pz_D2y_Dxz_S;
  abcd[528] = 2.0E0*I_ERI_F2xz_Dyz_Dxz_S_a;
  abcd[529] = 2.0E0*I_ERI_Fxyz_Dyz_Dxz_S_a;
  abcd[530] = 2.0E0*I_ERI_Fx2z_Dyz_Dxz_S_a-1*I_ERI_Px_Dyz_Dxz_S;
  abcd[531] = 2.0E0*I_ERI_F2yz_Dyz_Dxz_S_a;
  abcd[532] = 2.0E0*I_ERI_Fy2z_Dyz_Dxz_S_a-1*I_ERI_Py_Dyz_Dxz_S;
  abcd[533] = 2.0E0*I_ERI_F3z_Dyz_Dxz_S_a-2*I_ERI_Pz_Dyz_Dxz_S;
  abcd[534] = 2.0E0*I_ERI_F2xz_D2z_Dxz_S_a;
  abcd[535] = 2.0E0*I_ERI_Fxyz_D2z_Dxz_S_a;
  abcd[536] = 2.0E0*I_ERI_Fx2z_D2z_Dxz_S_a-1*I_ERI_Px_D2z_Dxz_S;
  abcd[537] = 2.0E0*I_ERI_F2yz_D2z_Dxz_S_a;
  abcd[538] = 2.0E0*I_ERI_Fy2z_D2z_Dxz_S_a-1*I_ERI_Py_D2z_Dxz_S;
  abcd[539] = 2.0E0*I_ERI_F3z_D2z_Dxz_S_a-2*I_ERI_Pz_D2z_Dxz_S;
  abcd[540] = 2.0E0*I_ERI_F2xz_D2x_D2y_S_a;
  abcd[541] = 2.0E0*I_ERI_Fxyz_D2x_D2y_S_a;
  abcd[542] = 2.0E0*I_ERI_Fx2z_D2x_D2y_S_a-1*I_ERI_Px_D2x_D2y_S;
  abcd[543] = 2.0E0*I_ERI_F2yz_D2x_D2y_S_a;
  abcd[544] = 2.0E0*I_ERI_Fy2z_D2x_D2y_S_a-1*I_ERI_Py_D2x_D2y_S;
  abcd[545] = 2.0E0*I_ERI_F3z_D2x_D2y_S_a-2*I_ERI_Pz_D2x_D2y_S;
  abcd[546] = 2.0E0*I_ERI_F2xz_Dxy_D2y_S_a;
  abcd[547] = 2.0E0*I_ERI_Fxyz_Dxy_D2y_S_a;
  abcd[548] = 2.0E0*I_ERI_Fx2z_Dxy_D2y_S_a-1*I_ERI_Px_Dxy_D2y_S;
  abcd[549] = 2.0E0*I_ERI_F2yz_Dxy_D2y_S_a;
  abcd[550] = 2.0E0*I_ERI_Fy2z_Dxy_D2y_S_a-1*I_ERI_Py_Dxy_D2y_S;
  abcd[551] = 2.0E0*I_ERI_F3z_Dxy_D2y_S_a-2*I_ERI_Pz_Dxy_D2y_S;
  abcd[552] = 2.0E0*I_ERI_F2xz_Dxz_D2y_S_a;
  abcd[553] = 2.0E0*I_ERI_Fxyz_Dxz_D2y_S_a;
  abcd[554] = 2.0E0*I_ERI_Fx2z_Dxz_D2y_S_a-1*I_ERI_Px_Dxz_D2y_S;
  abcd[555] = 2.0E0*I_ERI_F2yz_Dxz_D2y_S_a;
  abcd[556] = 2.0E0*I_ERI_Fy2z_Dxz_D2y_S_a-1*I_ERI_Py_Dxz_D2y_S;
  abcd[557] = 2.0E0*I_ERI_F3z_Dxz_D2y_S_a-2*I_ERI_Pz_Dxz_D2y_S;
  abcd[558] = 2.0E0*I_ERI_F2xz_D2y_D2y_S_a;
  abcd[559] = 2.0E0*I_ERI_Fxyz_D2y_D2y_S_a;
  abcd[560] = 2.0E0*I_ERI_Fx2z_D2y_D2y_S_a-1*I_ERI_Px_D2y_D2y_S;
  abcd[561] = 2.0E0*I_ERI_F2yz_D2y_D2y_S_a;
  abcd[562] = 2.0E0*I_ERI_Fy2z_D2y_D2y_S_a-1*I_ERI_Py_D2y_D2y_S;
  abcd[563] = 2.0E0*I_ERI_F3z_D2y_D2y_S_a-2*I_ERI_Pz_D2y_D2y_S;
  abcd[564] = 2.0E0*I_ERI_F2xz_Dyz_D2y_S_a;
  abcd[565] = 2.0E0*I_ERI_Fxyz_Dyz_D2y_S_a;
  abcd[566] = 2.0E0*I_ERI_Fx2z_Dyz_D2y_S_a-1*I_ERI_Px_Dyz_D2y_S;
  abcd[567] = 2.0E0*I_ERI_F2yz_Dyz_D2y_S_a;
  abcd[568] = 2.0E0*I_ERI_Fy2z_Dyz_D2y_S_a-1*I_ERI_Py_Dyz_D2y_S;
  abcd[569] = 2.0E0*I_ERI_F3z_Dyz_D2y_S_a-2*I_ERI_Pz_Dyz_D2y_S;
  abcd[570] = 2.0E0*I_ERI_F2xz_D2z_D2y_S_a;
  abcd[571] = 2.0E0*I_ERI_Fxyz_D2z_D2y_S_a;
  abcd[572] = 2.0E0*I_ERI_Fx2z_D2z_D2y_S_a-1*I_ERI_Px_D2z_D2y_S;
  abcd[573] = 2.0E0*I_ERI_F2yz_D2z_D2y_S_a;
  abcd[574] = 2.0E0*I_ERI_Fy2z_D2z_D2y_S_a-1*I_ERI_Py_D2z_D2y_S;
  abcd[575] = 2.0E0*I_ERI_F3z_D2z_D2y_S_a-2*I_ERI_Pz_D2z_D2y_S;
  abcd[576] = 2.0E0*I_ERI_F2xz_D2x_Dyz_S_a;
  abcd[577] = 2.0E0*I_ERI_Fxyz_D2x_Dyz_S_a;
  abcd[578] = 2.0E0*I_ERI_Fx2z_D2x_Dyz_S_a-1*I_ERI_Px_D2x_Dyz_S;
  abcd[579] = 2.0E0*I_ERI_F2yz_D2x_Dyz_S_a;
  abcd[580] = 2.0E0*I_ERI_Fy2z_D2x_Dyz_S_a-1*I_ERI_Py_D2x_Dyz_S;
  abcd[581] = 2.0E0*I_ERI_F3z_D2x_Dyz_S_a-2*I_ERI_Pz_D2x_Dyz_S;
  abcd[582] = 2.0E0*I_ERI_F2xz_Dxy_Dyz_S_a;
  abcd[583] = 2.0E0*I_ERI_Fxyz_Dxy_Dyz_S_a;
  abcd[584] = 2.0E0*I_ERI_Fx2z_Dxy_Dyz_S_a-1*I_ERI_Px_Dxy_Dyz_S;
  abcd[585] = 2.0E0*I_ERI_F2yz_Dxy_Dyz_S_a;
  abcd[586] = 2.0E0*I_ERI_Fy2z_Dxy_Dyz_S_a-1*I_ERI_Py_Dxy_Dyz_S;
  abcd[587] = 2.0E0*I_ERI_F3z_Dxy_Dyz_S_a-2*I_ERI_Pz_Dxy_Dyz_S;
  abcd[588] = 2.0E0*I_ERI_F2xz_Dxz_Dyz_S_a;
  abcd[589] = 2.0E0*I_ERI_Fxyz_Dxz_Dyz_S_a;
  abcd[590] = 2.0E0*I_ERI_Fx2z_Dxz_Dyz_S_a-1*I_ERI_Px_Dxz_Dyz_S;
  abcd[591] = 2.0E0*I_ERI_F2yz_Dxz_Dyz_S_a;
  abcd[592] = 2.0E0*I_ERI_Fy2z_Dxz_Dyz_S_a-1*I_ERI_Py_Dxz_Dyz_S;
  abcd[593] = 2.0E0*I_ERI_F3z_Dxz_Dyz_S_a-2*I_ERI_Pz_Dxz_Dyz_S;
  abcd[594] = 2.0E0*I_ERI_F2xz_D2y_Dyz_S_a;
  abcd[595] = 2.0E0*I_ERI_Fxyz_D2y_Dyz_S_a;
  abcd[596] = 2.0E0*I_ERI_Fx2z_D2y_Dyz_S_a-1*I_ERI_Px_D2y_Dyz_S;
  abcd[597] = 2.0E0*I_ERI_F2yz_D2y_Dyz_S_a;
  abcd[598] = 2.0E0*I_ERI_Fy2z_D2y_Dyz_S_a-1*I_ERI_Py_D2y_Dyz_S;
  abcd[599] = 2.0E0*I_ERI_F3z_D2y_Dyz_S_a-2*I_ERI_Pz_D2y_Dyz_S;
  abcd[600] = 2.0E0*I_ERI_F2xz_Dyz_Dyz_S_a;
  abcd[601] = 2.0E0*I_ERI_Fxyz_Dyz_Dyz_S_a;
  abcd[602] = 2.0E0*I_ERI_Fx2z_Dyz_Dyz_S_a-1*I_ERI_Px_Dyz_Dyz_S;
  abcd[603] = 2.0E0*I_ERI_F2yz_Dyz_Dyz_S_a;
  abcd[604] = 2.0E0*I_ERI_Fy2z_Dyz_Dyz_S_a-1*I_ERI_Py_Dyz_Dyz_S;
  abcd[605] = 2.0E0*I_ERI_F3z_Dyz_Dyz_S_a-2*I_ERI_Pz_Dyz_Dyz_S;
  abcd[606] = 2.0E0*I_ERI_F2xz_D2z_Dyz_S_a;
  abcd[607] = 2.0E0*I_ERI_Fxyz_D2z_Dyz_S_a;
  abcd[608] = 2.0E0*I_ERI_Fx2z_D2z_Dyz_S_a-1*I_ERI_Px_D2z_Dyz_S;
  abcd[609] = 2.0E0*I_ERI_F2yz_D2z_Dyz_S_a;
  abcd[610] = 2.0E0*I_ERI_Fy2z_D2z_Dyz_S_a-1*I_ERI_Py_D2z_Dyz_S;
  abcd[611] = 2.0E0*I_ERI_F3z_D2z_Dyz_S_a-2*I_ERI_Pz_D2z_Dyz_S;
  abcd[612] = 2.0E0*I_ERI_F2xz_D2x_D2z_S_a;
  abcd[613] = 2.0E0*I_ERI_Fxyz_D2x_D2z_S_a;
  abcd[614] = 2.0E0*I_ERI_Fx2z_D2x_D2z_S_a-1*I_ERI_Px_D2x_D2z_S;
  abcd[615] = 2.0E0*I_ERI_F2yz_D2x_D2z_S_a;
  abcd[616] = 2.0E0*I_ERI_Fy2z_D2x_D2z_S_a-1*I_ERI_Py_D2x_D2z_S;
  abcd[617] = 2.0E0*I_ERI_F3z_D2x_D2z_S_a-2*I_ERI_Pz_D2x_D2z_S;
  abcd[618] = 2.0E0*I_ERI_F2xz_Dxy_D2z_S_a;
  abcd[619] = 2.0E0*I_ERI_Fxyz_Dxy_D2z_S_a;
  abcd[620] = 2.0E0*I_ERI_Fx2z_Dxy_D2z_S_a-1*I_ERI_Px_Dxy_D2z_S;
  abcd[621] = 2.0E0*I_ERI_F2yz_Dxy_D2z_S_a;
  abcd[622] = 2.0E0*I_ERI_Fy2z_Dxy_D2z_S_a-1*I_ERI_Py_Dxy_D2z_S;
  abcd[623] = 2.0E0*I_ERI_F3z_Dxy_D2z_S_a-2*I_ERI_Pz_Dxy_D2z_S;
  abcd[624] = 2.0E0*I_ERI_F2xz_Dxz_D2z_S_a;
  abcd[625] = 2.0E0*I_ERI_Fxyz_Dxz_D2z_S_a;
  abcd[626] = 2.0E0*I_ERI_Fx2z_Dxz_D2z_S_a-1*I_ERI_Px_Dxz_D2z_S;
  abcd[627] = 2.0E0*I_ERI_F2yz_Dxz_D2z_S_a;
  abcd[628] = 2.0E0*I_ERI_Fy2z_Dxz_D2z_S_a-1*I_ERI_Py_Dxz_D2z_S;
  abcd[629] = 2.0E0*I_ERI_F3z_Dxz_D2z_S_a-2*I_ERI_Pz_Dxz_D2z_S;
  abcd[630] = 2.0E0*I_ERI_F2xz_D2y_D2z_S_a;
  abcd[631] = 2.0E0*I_ERI_Fxyz_D2y_D2z_S_a;
  abcd[632] = 2.0E0*I_ERI_Fx2z_D2y_D2z_S_a-1*I_ERI_Px_D2y_D2z_S;
  abcd[633] = 2.0E0*I_ERI_F2yz_D2y_D2z_S_a;
  abcd[634] = 2.0E0*I_ERI_Fy2z_D2y_D2z_S_a-1*I_ERI_Py_D2y_D2z_S;
  abcd[635] = 2.0E0*I_ERI_F3z_D2y_D2z_S_a-2*I_ERI_Pz_D2y_D2z_S;
  abcd[636] = 2.0E0*I_ERI_F2xz_Dyz_D2z_S_a;
  abcd[637] = 2.0E0*I_ERI_Fxyz_Dyz_D2z_S_a;
  abcd[638] = 2.0E0*I_ERI_Fx2z_Dyz_D2z_S_a-1*I_ERI_Px_Dyz_D2z_S;
  abcd[639] = 2.0E0*I_ERI_F2yz_Dyz_D2z_S_a;
  abcd[640] = 2.0E0*I_ERI_Fy2z_Dyz_D2z_S_a-1*I_ERI_Py_Dyz_D2z_S;
  abcd[641] = 2.0E0*I_ERI_F3z_Dyz_D2z_S_a-2*I_ERI_Pz_Dyz_D2z_S;
  abcd[642] = 2.0E0*I_ERI_F2xz_D2z_D2z_S_a;
  abcd[643] = 2.0E0*I_ERI_Fxyz_D2z_D2z_S_a;
  abcd[644] = 2.0E0*I_ERI_Fx2z_D2z_D2z_S_a-1*I_ERI_Px_D2z_D2z_S;
  abcd[645] = 2.0E0*I_ERI_F2yz_D2z_D2z_S_a;
  abcd[646] = 2.0E0*I_ERI_Fy2z_D2z_D2z_S_a-1*I_ERI_Py_D2z_D2z_S;
  abcd[647] = 2.0E0*I_ERI_F3z_D2z_D2z_S_a-2*I_ERI_Pz_D2z_D2z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_D_S_b
   * RHS shell quartet name: SQ_ERI_D_P_D_S
   ************************************************************/
  abcd[648] = 2.0E0*I_ERI_D2x_F3x_D2x_S_b-2*I_ERI_D2x_Px_D2x_S;
  abcd[649] = 2.0E0*I_ERI_Dxy_F3x_D2x_S_b-2*I_ERI_Dxy_Px_D2x_S;
  abcd[650] = 2.0E0*I_ERI_Dxz_F3x_D2x_S_b-2*I_ERI_Dxz_Px_D2x_S;
  abcd[651] = 2.0E0*I_ERI_D2y_F3x_D2x_S_b-2*I_ERI_D2y_Px_D2x_S;
  abcd[652] = 2.0E0*I_ERI_Dyz_F3x_D2x_S_b-2*I_ERI_Dyz_Px_D2x_S;
  abcd[653] = 2.0E0*I_ERI_D2z_F3x_D2x_S_b-2*I_ERI_D2z_Px_D2x_S;
  abcd[654] = 2.0E0*I_ERI_D2x_F2xy_D2x_S_b-1*I_ERI_D2x_Py_D2x_S;
  abcd[655] = 2.0E0*I_ERI_Dxy_F2xy_D2x_S_b-1*I_ERI_Dxy_Py_D2x_S;
  abcd[656] = 2.0E0*I_ERI_Dxz_F2xy_D2x_S_b-1*I_ERI_Dxz_Py_D2x_S;
  abcd[657] = 2.0E0*I_ERI_D2y_F2xy_D2x_S_b-1*I_ERI_D2y_Py_D2x_S;
  abcd[658] = 2.0E0*I_ERI_Dyz_F2xy_D2x_S_b-1*I_ERI_Dyz_Py_D2x_S;
  abcd[659] = 2.0E0*I_ERI_D2z_F2xy_D2x_S_b-1*I_ERI_D2z_Py_D2x_S;
  abcd[660] = 2.0E0*I_ERI_D2x_F2xz_D2x_S_b-1*I_ERI_D2x_Pz_D2x_S;
  abcd[661] = 2.0E0*I_ERI_Dxy_F2xz_D2x_S_b-1*I_ERI_Dxy_Pz_D2x_S;
  abcd[662] = 2.0E0*I_ERI_Dxz_F2xz_D2x_S_b-1*I_ERI_Dxz_Pz_D2x_S;
  abcd[663] = 2.0E0*I_ERI_D2y_F2xz_D2x_S_b-1*I_ERI_D2y_Pz_D2x_S;
  abcd[664] = 2.0E0*I_ERI_Dyz_F2xz_D2x_S_b-1*I_ERI_Dyz_Pz_D2x_S;
  abcd[665] = 2.0E0*I_ERI_D2z_F2xz_D2x_S_b-1*I_ERI_D2z_Pz_D2x_S;
  abcd[666] = 2.0E0*I_ERI_D2x_Fx2y_D2x_S_b;
  abcd[667] = 2.0E0*I_ERI_Dxy_Fx2y_D2x_S_b;
  abcd[668] = 2.0E0*I_ERI_Dxz_Fx2y_D2x_S_b;
  abcd[669] = 2.0E0*I_ERI_D2y_Fx2y_D2x_S_b;
  abcd[670] = 2.0E0*I_ERI_Dyz_Fx2y_D2x_S_b;
  abcd[671] = 2.0E0*I_ERI_D2z_Fx2y_D2x_S_b;
  abcd[672] = 2.0E0*I_ERI_D2x_Fxyz_D2x_S_b;
  abcd[673] = 2.0E0*I_ERI_Dxy_Fxyz_D2x_S_b;
  abcd[674] = 2.0E0*I_ERI_Dxz_Fxyz_D2x_S_b;
  abcd[675] = 2.0E0*I_ERI_D2y_Fxyz_D2x_S_b;
  abcd[676] = 2.0E0*I_ERI_Dyz_Fxyz_D2x_S_b;
  abcd[677] = 2.0E0*I_ERI_D2z_Fxyz_D2x_S_b;
  abcd[678] = 2.0E0*I_ERI_D2x_Fx2z_D2x_S_b;
  abcd[679] = 2.0E0*I_ERI_Dxy_Fx2z_D2x_S_b;
  abcd[680] = 2.0E0*I_ERI_Dxz_Fx2z_D2x_S_b;
  abcd[681] = 2.0E0*I_ERI_D2y_Fx2z_D2x_S_b;
  abcd[682] = 2.0E0*I_ERI_Dyz_Fx2z_D2x_S_b;
  abcd[683] = 2.0E0*I_ERI_D2z_Fx2z_D2x_S_b;
  abcd[684] = 2.0E0*I_ERI_D2x_F3x_Dxy_S_b-2*I_ERI_D2x_Px_Dxy_S;
  abcd[685] = 2.0E0*I_ERI_Dxy_F3x_Dxy_S_b-2*I_ERI_Dxy_Px_Dxy_S;
  abcd[686] = 2.0E0*I_ERI_Dxz_F3x_Dxy_S_b-2*I_ERI_Dxz_Px_Dxy_S;
  abcd[687] = 2.0E0*I_ERI_D2y_F3x_Dxy_S_b-2*I_ERI_D2y_Px_Dxy_S;
  abcd[688] = 2.0E0*I_ERI_Dyz_F3x_Dxy_S_b-2*I_ERI_Dyz_Px_Dxy_S;
  abcd[689] = 2.0E0*I_ERI_D2z_F3x_Dxy_S_b-2*I_ERI_D2z_Px_Dxy_S;
  abcd[690] = 2.0E0*I_ERI_D2x_F2xy_Dxy_S_b-1*I_ERI_D2x_Py_Dxy_S;
  abcd[691] = 2.0E0*I_ERI_Dxy_F2xy_Dxy_S_b-1*I_ERI_Dxy_Py_Dxy_S;
  abcd[692] = 2.0E0*I_ERI_Dxz_F2xy_Dxy_S_b-1*I_ERI_Dxz_Py_Dxy_S;
  abcd[693] = 2.0E0*I_ERI_D2y_F2xy_Dxy_S_b-1*I_ERI_D2y_Py_Dxy_S;
  abcd[694] = 2.0E0*I_ERI_Dyz_F2xy_Dxy_S_b-1*I_ERI_Dyz_Py_Dxy_S;
  abcd[695] = 2.0E0*I_ERI_D2z_F2xy_Dxy_S_b-1*I_ERI_D2z_Py_Dxy_S;
  abcd[696] = 2.0E0*I_ERI_D2x_F2xz_Dxy_S_b-1*I_ERI_D2x_Pz_Dxy_S;
  abcd[697] = 2.0E0*I_ERI_Dxy_F2xz_Dxy_S_b-1*I_ERI_Dxy_Pz_Dxy_S;
  abcd[698] = 2.0E0*I_ERI_Dxz_F2xz_Dxy_S_b-1*I_ERI_Dxz_Pz_Dxy_S;
  abcd[699] = 2.0E0*I_ERI_D2y_F2xz_Dxy_S_b-1*I_ERI_D2y_Pz_Dxy_S;
  abcd[700] = 2.0E0*I_ERI_Dyz_F2xz_Dxy_S_b-1*I_ERI_Dyz_Pz_Dxy_S;
  abcd[701] = 2.0E0*I_ERI_D2z_F2xz_Dxy_S_b-1*I_ERI_D2z_Pz_Dxy_S;
  abcd[702] = 2.0E0*I_ERI_D2x_Fx2y_Dxy_S_b;
  abcd[703] = 2.0E0*I_ERI_Dxy_Fx2y_Dxy_S_b;
  abcd[704] = 2.0E0*I_ERI_Dxz_Fx2y_Dxy_S_b;
  abcd[705] = 2.0E0*I_ERI_D2y_Fx2y_Dxy_S_b;
  abcd[706] = 2.0E0*I_ERI_Dyz_Fx2y_Dxy_S_b;
  abcd[707] = 2.0E0*I_ERI_D2z_Fx2y_Dxy_S_b;
  abcd[708] = 2.0E0*I_ERI_D2x_Fxyz_Dxy_S_b;
  abcd[709] = 2.0E0*I_ERI_Dxy_Fxyz_Dxy_S_b;
  abcd[710] = 2.0E0*I_ERI_Dxz_Fxyz_Dxy_S_b;
  abcd[711] = 2.0E0*I_ERI_D2y_Fxyz_Dxy_S_b;
  abcd[712] = 2.0E0*I_ERI_Dyz_Fxyz_Dxy_S_b;
  abcd[713] = 2.0E0*I_ERI_D2z_Fxyz_Dxy_S_b;
  abcd[714] = 2.0E0*I_ERI_D2x_Fx2z_Dxy_S_b;
  abcd[715] = 2.0E0*I_ERI_Dxy_Fx2z_Dxy_S_b;
  abcd[716] = 2.0E0*I_ERI_Dxz_Fx2z_Dxy_S_b;
  abcd[717] = 2.0E0*I_ERI_D2y_Fx2z_Dxy_S_b;
  abcd[718] = 2.0E0*I_ERI_Dyz_Fx2z_Dxy_S_b;
  abcd[719] = 2.0E0*I_ERI_D2z_Fx2z_Dxy_S_b;
  abcd[720] = 2.0E0*I_ERI_D2x_F3x_Dxz_S_b-2*I_ERI_D2x_Px_Dxz_S;
  abcd[721] = 2.0E0*I_ERI_Dxy_F3x_Dxz_S_b-2*I_ERI_Dxy_Px_Dxz_S;
  abcd[722] = 2.0E0*I_ERI_Dxz_F3x_Dxz_S_b-2*I_ERI_Dxz_Px_Dxz_S;
  abcd[723] = 2.0E0*I_ERI_D2y_F3x_Dxz_S_b-2*I_ERI_D2y_Px_Dxz_S;
  abcd[724] = 2.0E0*I_ERI_Dyz_F3x_Dxz_S_b-2*I_ERI_Dyz_Px_Dxz_S;
  abcd[725] = 2.0E0*I_ERI_D2z_F3x_Dxz_S_b-2*I_ERI_D2z_Px_Dxz_S;
  abcd[726] = 2.0E0*I_ERI_D2x_F2xy_Dxz_S_b-1*I_ERI_D2x_Py_Dxz_S;
  abcd[727] = 2.0E0*I_ERI_Dxy_F2xy_Dxz_S_b-1*I_ERI_Dxy_Py_Dxz_S;
  abcd[728] = 2.0E0*I_ERI_Dxz_F2xy_Dxz_S_b-1*I_ERI_Dxz_Py_Dxz_S;
  abcd[729] = 2.0E0*I_ERI_D2y_F2xy_Dxz_S_b-1*I_ERI_D2y_Py_Dxz_S;
  abcd[730] = 2.0E0*I_ERI_Dyz_F2xy_Dxz_S_b-1*I_ERI_Dyz_Py_Dxz_S;
  abcd[731] = 2.0E0*I_ERI_D2z_F2xy_Dxz_S_b-1*I_ERI_D2z_Py_Dxz_S;
  abcd[732] = 2.0E0*I_ERI_D2x_F2xz_Dxz_S_b-1*I_ERI_D2x_Pz_Dxz_S;
  abcd[733] = 2.0E0*I_ERI_Dxy_F2xz_Dxz_S_b-1*I_ERI_Dxy_Pz_Dxz_S;
  abcd[734] = 2.0E0*I_ERI_Dxz_F2xz_Dxz_S_b-1*I_ERI_Dxz_Pz_Dxz_S;
  abcd[735] = 2.0E0*I_ERI_D2y_F2xz_Dxz_S_b-1*I_ERI_D2y_Pz_Dxz_S;
  abcd[736] = 2.0E0*I_ERI_Dyz_F2xz_Dxz_S_b-1*I_ERI_Dyz_Pz_Dxz_S;
  abcd[737] = 2.0E0*I_ERI_D2z_F2xz_Dxz_S_b-1*I_ERI_D2z_Pz_Dxz_S;
  abcd[738] = 2.0E0*I_ERI_D2x_Fx2y_Dxz_S_b;
  abcd[739] = 2.0E0*I_ERI_Dxy_Fx2y_Dxz_S_b;
  abcd[740] = 2.0E0*I_ERI_Dxz_Fx2y_Dxz_S_b;
  abcd[741] = 2.0E0*I_ERI_D2y_Fx2y_Dxz_S_b;
  abcd[742] = 2.0E0*I_ERI_Dyz_Fx2y_Dxz_S_b;
  abcd[743] = 2.0E0*I_ERI_D2z_Fx2y_Dxz_S_b;
  abcd[744] = 2.0E0*I_ERI_D2x_Fxyz_Dxz_S_b;
  abcd[745] = 2.0E0*I_ERI_Dxy_Fxyz_Dxz_S_b;
  abcd[746] = 2.0E0*I_ERI_Dxz_Fxyz_Dxz_S_b;
  abcd[747] = 2.0E0*I_ERI_D2y_Fxyz_Dxz_S_b;
  abcd[748] = 2.0E0*I_ERI_Dyz_Fxyz_Dxz_S_b;
  abcd[749] = 2.0E0*I_ERI_D2z_Fxyz_Dxz_S_b;
  abcd[750] = 2.0E0*I_ERI_D2x_Fx2z_Dxz_S_b;
  abcd[751] = 2.0E0*I_ERI_Dxy_Fx2z_Dxz_S_b;
  abcd[752] = 2.0E0*I_ERI_Dxz_Fx2z_Dxz_S_b;
  abcd[753] = 2.0E0*I_ERI_D2y_Fx2z_Dxz_S_b;
  abcd[754] = 2.0E0*I_ERI_Dyz_Fx2z_Dxz_S_b;
  abcd[755] = 2.0E0*I_ERI_D2z_Fx2z_Dxz_S_b;
  abcd[756] = 2.0E0*I_ERI_D2x_F3x_D2y_S_b-2*I_ERI_D2x_Px_D2y_S;
  abcd[757] = 2.0E0*I_ERI_Dxy_F3x_D2y_S_b-2*I_ERI_Dxy_Px_D2y_S;
  abcd[758] = 2.0E0*I_ERI_Dxz_F3x_D2y_S_b-2*I_ERI_Dxz_Px_D2y_S;
  abcd[759] = 2.0E0*I_ERI_D2y_F3x_D2y_S_b-2*I_ERI_D2y_Px_D2y_S;
  abcd[760] = 2.0E0*I_ERI_Dyz_F3x_D2y_S_b-2*I_ERI_Dyz_Px_D2y_S;
  abcd[761] = 2.0E0*I_ERI_D2z_F3x_D2y_S_b-2*I_ERI_D2z_Px_D2y_S;
  abcd[762] = 2.0E0*I_ERI_D2x_F2xy_D2y_S_b-1*I_ERI_D2x_Py_D2y_S;
  abcd[763] = 2.0E0*I_ERI_Dxy_F2xy_D2y_S_b-1*I_ERI_Dxy_Py_D2y_S;
  abcd[764] = 2.0E0*I_ERI_Dxz_F2xy_D2y_S_b-1*I_ERI_Dxz_Py_D2y_S;
  abcd[765] = 2.0E0*I_ERI_D2y_F2xy_D2y_S_b-1*I_ERI_D2y_Py_D2y_S;
  abcd[766] = 2.0E0*I_ERI_Dyz_F2xy_D2y_S_b-1*I_ERI_Dyz_Py_D2y_S;
  abcd[767] = 2.0E0*I_ERI_D2z_F2xy_D2y_S_b-1*I_ERI_D2z_Py_D2y_S;
  abcd[768] = 2.0E0*I_ERI_D2x_F2xz_D2y_S_b-1*I_ERI_D2x_Pz_D2y_S;
  abcd[769] = 2.0E0*I_ERI_Dxy_F2xz_D2y_S_b-1*I_ERI_Dxy_Pz_D2y_S;
  abcd[770] = 2.0E0*I_ERI_Dxz_F2xz_D2y_S_b-1*I_ERI_Dxz_Pz_D2y_S;
  abcd[771] = 2.0E0*I_ERI_D2y_F2xz_D2y_S_b-1*I_ERI_D2y_Pz_D2y_S;
  abcd[772] = 2.0E0*I_ERI_Dyz_F2xz_D2y_S_b-1*I_ERI_Dyz_Pz_D2y_S;
  abcd[773] = 2.0E0*I_ERI_D2z_F2xz_D2y_S_b-1*I_ERI_D2z_Pz_D2y_S;
  abcd[774] = 2.0E0*I_ERI_D2x_Fx2y_D2y_S_b;
  abcd[775] = 2.0E0*I_ERI_Dxy_Fx2y_D2y_S_b;
  abcd[776] = 2.0E0*I_ERI_Dxz_Fx2y_D2y_S_b;
  abcd[777] = 2.0E0*I_ERI_D2y_Fx2y_D2y_S_b;
  abcd[778] = 2.0E0*I_ERI_Dyz_Fx2y_D2y_S_b;
  abcd[779] = 2.0E0*I_ERI_D2z_Fx2y_D2y_S_b;
  abcd[780] = 2.0E0*I_ERI_D2x_Fxyz_D2y_S_b;
  abcd[781] = 2.0E0*I_ERI_Dxy_Fxyz_D2y_S_b;
  abcd[782] = 2.0E0*I_ERI_Dxz_Fxyz_D2y_S_b;
  abcd[783] = 2.0E0*I_ERI_D2y_Fxyz_D2y_S_b;
  abcd[784] = 2.0E0*I_ERI_Dyz_Fxyz_D2y_S_b;
  abcd[785] = 2.0E0*I_ERI_D2z_Fxyz_D2y_S_b;
  abcd[786] = 2.0E0*I_ERI_D2x_Fx2z_D2y_S_b;
  abcd[787] = 2.0E0*I_ERI_Dxy_Fx2z_D2y_S_b;
  abcd[788] = 2.0E0*I_ERI_Dxz_Fx2z_D2y_S_b;
  abcd[789] = 2.0E0*I_ERI_D2y_Fx2z_D2y_S_b;
  abcd[790] = 2.0E0*I_ERI_Dyz_Fx2z_D2y_S_b;
  abcd[791] = 2.0E0*I_ERI_D2z_Fx2z_D2y_S_b;
  abcd[792] = 2.0E0*I_ERI_D2x_F3x_Dyz_S_b-2*I_ERI_D2x_Px_Dyz_S;
  abcd[793] = 2.0E0*I_ERI_Dxy_F3x_Dyz_S_b-2*I_ERI_Dxy_Px_Dyz_S;
  abcd[794] = 2.0E0*I_ERI_Dxz_F3x_Dyz_S_b-2*I_ERI_Dxz_Px_Dyz_S;
  abcd[795] = 2.0E0*I_ERI_D2y_F3x_Dyz_S_b-2*I_ERI_D2y_Px_Dyz_S;
  abcd[796] = 2.0E0*I_ERI_Dyz_F3x_Dyz_S_b-2*I_ERI_Dyz_Px_Dyz_S;
  abcd[797] = 2.0E0*I_ERI_D2z_F3x_Dyz_S_b-2*I_ERI_D2z_Px_Dyz_S;
  abcd[798] = 2.0E0*I_ERI_D2x_F2xy_Dyz_S_b-1*I_ERI_D2x_Py_Dyz_S;
  abcd[799] = 2.0E0*I_ERI_Dxy_F2xy_Dyz_S_b-1*I_ERI_Dxy_Py_Dyz_S;
  abcd[800] = 2.0E0*I_ERI_Dxz_F2xy_Dyz_S_b-1*I_ERI_Dxz_Py_Dyz_S;
  abcd[801] = 2.0E0*I_ERI_D2y_F2xy_Dyz_S_b-1*I_ERI_D2y_Py_Dyz_S;
  abcd[802] = 2.0E0*I_ERI_Dyz_F2xy_Dyz_S_b-1*I_ERI_Dyz_Py_Dyz_S;
  abcd[803] = 2.0E0*I_ERI_D2z_F2xy_Dyz_S_b-1*I_ERI_D2z_Py_Dyz_S;
  abcd[804] = 2.0E0*I_ERI_D2x_F2xz_Dyz_S_b-1*I_ERI_D2x_Pz_Dyz_S;
  abcd[805] = 2.0E0*I_ERI_Dxy_F2xz_Dyz_S_b-1*I_ERI_Dxy_Pz_Dyz_S;
  abcd[806] = 2.0E0*I_ERI_Dxz_F2xz_Dyz_S_b-1*I_ERI_Dxz_Pz_Dyz_S;
  abcd[807] = 2.0E0*I_ERI_D2y_F2xz_Dyz_S_b-1*I_ERI_D2y_Pz_Dyz_S;
  abcd[808] = 2.0E0*I_ERI_Dyz_F2xz_Dyz_S_b-1*I_ERI_Dyz_Pz_Dyz_S;
  abcd[809] = 2.0E0*I_ERI_D2z_F2xz_Dyz_S_b-1*I_ERI_D2z_Pz_Dyz_S;
  abcd[810] = 2.0E0*I_ERI_D2x_Fx2y_Dyz_S_b;
  abcd[811] = 2.0E0*I_ERI_Dxy_Fx2y_Dyz_S_b;
  abcd[812] = 2.0E0*I_ERI_Dxz_Fx2y_Dyz_S_b;
  abcd[813] = 2.0E0*I_ERI_D2y_Fx2y_Dyz_S_b;
  abcd[814] = 2.0E0*I_ERI_Dyz_Fx2y_Dyz_S_b;
  abcd[815] = 2.0E0*I_ERI_D2z_Fx2y_Dyz_S_b;
  abcd[816] = 2.0E0*I_ERI_D2x_Fxyz_Dyz_S_b;
  abcd[817] = 2.0E0*I_ERI_Dxy_Fxyz_Dyz_S_b;
  abcd[818] = 2.0E0*I_ERI_Dxz_Fxyz_Dyz_S_b;
  abcd[819] = 2.0E0*I_ERI_D2y_Fxyz_Dyz_S_b;
  abcd[820] = 2.0E0*I_ERI_Dyz_Fxyz_Dyz_S_b;
  abcd[821] = 2.0E0*I_ERI_D2z_Fxyz_Dyz_S_b;
  abcd[822] = 2.0E0*I_ERI_D2x_Fx2z_Dyz_S_b;
  abcd[823] = 2.0E0*I_ERI_Dxy_Fx2z_Dyz_S_b;
  abcd[824] = 2.0E0*I_ERI_Dxz_Fx2z_Dyz_S_b;
  abcd[825] = 2.0E0*I_ERI_D2y_Fx2z_Dyz_S_b;
  abcd[826] = 2.0E0*I_ERI_Dyz_Fx2z_Dyz_S_b;
  abcd[827] = 2.0E0*I_ERI_D2z_Fx2z_Dyz_S_b;
  abcd[828] = 2.0E0*I_ERI_D2x_F3x_D2z_S_b-2*I_ERI_D2x_Px_D2z_S;
  abcd[829] = 2.0E0*I_ERI_Dxy_F3x_D2z_S_b-2*I_ERI_Dxy_Px_D2z_S;
  abcd[830] = 2.0E0*I_ERI_Dxz_F3x_D2z_S_b-2*I_ERI_Dxz_Px_D2z_S;
  abcd[831] = 2.0E0*I_ERI_D2y_F3x_D2z_S_b-2*I_ERI_D2y_Px_D2z_S;
  abcd[832] = 2.0E0*I_ERI_Dyz_F3x_D2z_S_b-2*I_ERI_Dyz_Px_D2z_S;
  abcd[833] = 2.0E0*I_ERI_D2z_F3x_D2z_S_b-2*I_ERI_D2z_Px_D2z_S;
  abcd[834] = 2.0E0*I_ERI_D2x_F2xy_D2z_S_b-1*I_ERI_D2x_Py_D2z_S;
  abcd[835] = 2.0E0*I_ERI_Dxy_F2xy_D2z_S_b-1*I_ERI_Dxy_Py_D2z_S;
  abcd[836] = 2.0E0*I_ERI_Dxz_F2xy_D2z_S_b-1*I_ERI_Dxz_Py_D2z_S;
  abcd[837] = 2.0E0*I_ERI_D2y_F2xy_D2z_S_b-1*I_ERI_D2y_Py_D2z_S;
  abcd[838] = 2.0E0*I_ERI_Dyz_F2xy_D2z_S_b-1*I_ERI_Dyz_Py_D2z_S;
  abcd[839] = 2.0E0*I_ERI_D2z_F2xy_D2z_S_b-1*I_ERI_D2z_Py_D2z_S;
  abcd[840] = 2.0E0*I_ERI_D2x_F2xz_D2z_S_b-1*I_ERI_D2x_Pz_D2z_S;
  abcd[841] = 2.0E0*I_ERI_Dxy_F2xz_D2z_S_b-1*I_ERI_Dxy_Pz_D2z_S;
  abcd[842] = 2.0E0*I_ERI_Dxz_F2xz_D2z_S_b-1*I_ERI_Dxz_Pz_D2z_S;
  abcd[843] = 2.0E0*I_ERI_D2y_F2xz_D2z_S_b-1*I_ERI_D2y_Pz_D2z_S;
  abcd[844] = 2.0E0*I_ERI_Dyz_F2xz_D2z_S_b-1*I_ERI_Dyz_Pz_D2z_S;
  abcd[845] = 2.0E0*I_ERI_D2z_F2xz_D2z_S_b-1*I_ERI_D2z_Pz_D2z_S;
  abcd[846] = 2.0E0*I_ERI_D2x_Fx2y_D2z_S_b;
  abcd[847] = 2.0E0*I_ERI_Dxy_Fx2y_D2z_S_b;
  abcd[848] = 2.0E0*I_ERI_Dxz_Fx2y_D2z_S_b;
  abcd[849] = 2.0E0*I_ERI_D2y_Fx2y_D2z_S_b;
  abcd[850] = 2.0E0*I_ERI_Dyz_Fx2y_D2z_S_b;
  abcd[851] = 2.0E0*I_ERI_D2z_Fx2y_D2z_S_b;
  abcd[852] = 2.0E0*I_ERI_D2x_Fxyz_D2z_S_b;
  abcd[853] = 2.0E0*I_ERI_Dxy_Fxyz_D2z_S_b;
  abcd[854] = 2.0E0*I_ERI_Dxz_Fxyz_D2z_S_b;
  abcd[855] = 2.0E0*I_ERI_D2y_Fxyz_D2z_S_b;
  abcd[856] = 2.0E0*I_ERI_Dyz_Fxyz_D2z_S_b;
  abcd[857] = 2.0E0*I_ERI_D2z_Fxyz_D2z_S_b;
  abcd[858] = 2.0E0*I_ERI_D2x_Fx2z_D2z_S_b;
  abcd[859] = 2.0E0*I_ERI_Dxy_Fx2z_D2z_S_b;
  abcd[860] = 2.0E0*I_ERI_Dxz_Fx2z_D2z_S_b;
  abcd[861] = 2.0E0*I_ERI_D2y_Fx2z_D2z_S_b;
  abcd[862] = 2.0E0*I_ERI_Dyz_Fx2z_D2z_S_b;
  abcd[863] = 2.0E0*I_ERI_D2z_Fx2z_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_D_S_b
   * RHS shell quartet name: SQ_ERI_D_P_D_S
   ************************************************************/
  abcd[864] = 2.0E0*I_ERI_D2x_F2xy_D2x_S_b;
  abcd[865] = 2.0E0*I_ERI_Dxy_F2xy_D2x_S_b;
  abcd[866] = 2.0E0*I_ERI_Dxz_F2xy_D2x_S_b;
  abcd[867] = 2.0E0*I_ERI_D2y_F2xy_D2x_S_b;
  abcd[868] = 2.0E0*I_ERI_Dyz_F2xy_D2x_S_b;
  abcd[869] = 2.0E0*I_ERI_D2z_F2xy_D2x_S_b;
  abcd[870] = 2.0E0*I_ERI_D2x_Fx2y_D2x_S_b-1*I_ERI_D2x_Px_D2x_S;
  abcd[871] = 2.0E0*I_ERI_Dxy_Fx2y_D2x_S_b-1*I_ERI_Dxy_Px_D2x_S;
  abcd[872] = 2.0E0*I_ERI_Dxz_Fx2y_D2x_S_b-1*I_ERI_Dxz_Px_D2x_S;
  abcd[873] = 2.0E0*I_ERI_D2y_Fx2y_D2x_S_b-1*I_ERI_D2y_Px_D2x_S;
  abcd[874] = 2.0E0*I_ERI_Dyz_Fx2y_D2x_S_b-1*I_ERI_Dyz_Px_D2x_S;
  abcd[875] = 2.0E0*I_ERI_D2z_Fx2y_D2x_S_b-1*I_ERI_D2z_Px_D2x_S;
  abcd[876] = 2.0E0*I_ERI_D2x_Fxyz_D2x_S_b;
  abcd[877] = 2.0E0*I_ERI_Dxy_Fxyz_D2x_S_b;
  abcd[878] = 2.0E0*I_ERI_Dxz_Fxyz_D2x_S_b;
  abcd[879] = 2.0E0*I_ERI_D2y_Fxyz_D2x_S_b;
  abcd[880] = 2.0E0*I_ERI_Dyz_Fxyz_D2x_S_b;
  abcd[881] = 2.0E0*I_ERI_D2z_Fxyz_D2x_S_b;
  abcd[882] = 2.0E0*I_ERI_D2x_F3y_D2x_S_b-2*I_ERI_D2x_Py_D2x_S;
  abcd[883] = 2.0E0*I_ERI_Dxy_F3y_D2x_S_b-2*I_ERI_Dxy_Py_D2x_S;
  abcd[884] = 2.0E0*I_ERI_Dxz_F3y_D2x_S_b-2*I_ERI_Dxz_Py_D2x_S;
  abcd[885] = 2.0E0*I_ERI_D2y_F3y_D2x_S_b-2*I_ERI_D2y_Py_D2x_S;
  abcd[886] = 2.0E0*I_ERI_Dyz_F3y_D2x_S_b-2*I_ERI_Dyz_Py_D2x_S;
  abcd[887] = 2.0E0*I_ERI_D2z_F3y_D2x_S_b-2*I_ERI_D2z_Py_D2x_S;
  abcd[888] = 2.0E0*I_ERI_D2x_F2yz_D2x_S_b-1*I_ERI_D2x_Pz_D2x_S;
  abcd[889] = 2.0E0*I_ERI_Dxy_F2yz_D2x_S_b-1*I_ERI_Dxy_Pz_D2x_S;
  abcd[890] = 2.0E0*I_ERI_Dxz_F2yz_D2x_S_b-1*I_ERI_Dxz_Pz_D2x_S;
  abcd[891] = 2.0E0*I_ERI_D2y_F2yz_D2x_S_b-1*I_ERI_D2y_Pz_D2x_S;
  abcd[892] = 2.0E0*I_ERI_Dyz_F2yz_D2x_S_b-1*I_ERI_Dyz_Pz_D2x_S;
  abcd[893] = 2.0E0*I_ERI_D2z_F2yz_D2x_S_b-1*I_ERI_D2z_Pz_D2x_S;
  abcd[894] = 2.0E0*I_ERI_D2x_Fy2z_D2x_S_b;
  abcd[895] = 2.0E0*I_ERI_Dxy_Fy2z_D2x_S_b;
  abcd[896] = 2.0E0*I_ERI_Dxz_Fy2z_D2x_S_b;
  abcd[897] = 2.0E0*I_ERI_D2y_Fy2z_D2x_S_b;
  abcd[898] = 2.0E0*I_ERI_Dyz_Fy2z_D2x_S_b;
  abcd[899] = 2.0E0*I_ERI_D2z_Fy2z_D2x_S_b;
  abcd[900] = 2.0E0*I_ERI_D2x_F2xy_Dxy_S_b;
  abcd[901] = 2.0E0*I_ERI_Dxy_F2xy_Dxy_S_b;
  abcd[902] = 2.0E0*I_ERI_Dxz_F2xy_Dxy_S_b;
  abcd[903] = 2.0E0*I_ERI_D2y_F2xy_Dxy_S_b;
  abcd[904] = 2.0E0*I_ERI_Dyz_F2xy_Dxy_S_b;
  abcd[905] = 2.0E0*I_ERI_D2z_F2xy_Dxy_S_b;
  abcd[906] = 2.0E0*I_ERI_D2x_Fx2y_Dxy_S_b-1*I_ERI_D2x_Px_Dxy_S;
  abcd[907] = 2.0E0*I_ERI_Dxy_Fx2y_Dxy_S_b-1*I_ERI_Dxy_Px_Dxy_S;
  abcd[908] = 2.0E0*I_ERI_Dxz_Fx2y_Dxy_S_b-1*I_ERI_Dxz_Px_Dxy_S;
  abcd[909] = 2.0E0*I_ERI_D2y_Fx2y_Dxy_S_b-1*I_ERI_D2y_Px_Dxy_S;
  abcd[910] = 2.0E0*I_ERI_Dyz_Fx2y_Dxy_S_b-1*I_ERI_Dyz_Px_Dxy_S;
  abcd[911] = 2.0E0*I_ERI_D2z_Fx2y_Dxy_S_b-1*I_ERI_D2z_Px_Dxy_S;
  abcd[912] = 2.0E0*I_ERI_D2x_Fxyz_Dxy_S_b;
  abcd[913] = 2.0E0*I_ERI_Dxy_Fxyz_Dxy_S_b;
  abcd[914] = 2.0E0*I_ERI_Dxz_Fxyz_Dxy_S_b;
  abcd[915] = 2.0E0*I_ERI_D2y_Fxyz_Dxy_S_b;
  abcd[916] = 2.0E0*I_ERI_Dyz_Fxyz_Dxy_S_b;
  abcd[917] = 2.0E0*I_ERI_D2z_Fxyz_Dxy_S_b;
  abcd[918] = 2.0E0*I_ERI_D2x_F3y_Dxy_S_b-2*I_ERI_D2x_Py_Dxy_S;
  abcd[919] = 2.0E0*I_ERI_Dxy_F3y_Dxy_S_b-2*I_ERI_Dxy_Py_Dxy_S;
  abcd[920] = 2.0E0*I_ERI_Dxz_F3y_Dxy_S_b-2*I_ERI_Dxz_Py_Dxy_S;
  abcd[921] = 2.0E0*I_ERI_D2y_F3y_Dxy_S_b-2*I_ERI_D2y_Py_Dxy_S;
  abcd[922] = 2.0E0*I_ERI_Dyz_F3y_Dxy_S_b-2*I_ERI_Dyz_Py_Dxy_S;
  abcd[923] = 2.0E0*I_ERI_D2z_F3y_Dxy_S_b-2*I_ERI_D2z_Py_Dxy_S;
  abcd[924] = 2.0E0*I_ERI_D2x_F2yz_Dxy_S_b-1*I_ERI_D2x_Pz_Dxy_S;
  abcd[925] = 2.0E0*I_ERI_Dxy_F2yz_Dxy_S_b-1*I_ERI_Dxy_Pz_Dxy_S;
  abcd[926] = 2.0E0*I_ERI_Dxz_F2yz_Dxy_S_b-1*I_ERI_Dxz_Pz_Dxy_S;
  abcd[927] = 2.0E0*I_ERI_D2y_F2yz_Dxy_S_b-1*I_ERI_D2y_Pz_Dxy_S;
  abcd[928] = 2.0E0*I_ERI_Dyz_F2yz_Dxy_S_b-1*I_ERI_Dyz_Pz_Dxy_S;
  abcd[929] = 2.0E0*I_ERI_D2z_F2yz_Dxy_S_b-1*I_ERI_D2z_Pz_Dxy_S;
  abcd[930] = 2.0E0*I_ERI_D2x_Fy2z_Dxy_S_b;
  abcd[931] = 2.0E0*I_ERI_Dxy_Fy2z_Dxy_S_b;
  abcd[932] = 2.0E0*I_ERI_Dxz_Fy2z_Dxy_S_b;
  abcd[933] = 2.0E0*I_ERI_D2y_Fy2z_Dxy_S_b;
  abcd[934] = 2.0E0*I_ERI_Dyz_Fy2z_Dxy_S_b;
  abcd[935] = 2.0E0*I_ERI_D2z_Fy2z_Dxy_S_b;
  abcd[936] = 2.0E0*I_ERI_D2x_F2xy_Dxz_S_b;
  abcd[937] = 2.0E0*I_ERI_Dxy_F2xy_Dxz_S_b;
  abcd[938] = 2.0E0*I_ERI_Dxz_F2xy_Dxz_S_b;
  abcd[939] = 2.0E0*I_ERI_D2y_F2xy_Dxz_S_b;
  abcd[940] = 2.0E0*I_ERI_Dyz_F2xy_Dxz_S_b;
  abcd[941] = 2.0E0*I_ERI_D2z_F2xy_Dxz_S_b;
  abcd[942] = 2.0E0*I_ERI_D2x_Fx2y_Dxz_S_b-1*I_ERI_D2x_Px_Dxz_S;
  abcd[943] = 2.0E0*I_ERI_Dxy_Fx2y_Dxz_S_b-1*I_ERI_Dxy_Px_Dxz_S;
  abcd[944] = 2.0E0*I_ERI_Dxz_Fx2y_Dxz_S_b-1*I_ERI_Dxz_Px_Dxz_S;
  abcd[945] = 2.0E0*I_ERI_D2y_Fx2y_Dxz_S_b-1*I_ERI_D2y_Px_Dxz_S;
  abcd[946] = 2.0E0*I_ERI_Dyz_Fx2y_Dxz_S_b-1*I_ERI_Dyz_Px_Dxz_S;
  abcd[947] = 2.0E0*I_ERI_D2z_Fx2y_Dxz_S_b-1*I_ERI_D2z_Px_Dxz_S;
  abcd[948] = 2.0E0*I_ERI_D2x_Fxyz_Dxz_S_b;
  abcd[949] = 2.0E0*I_ERI_Dxy_Fxyz_Dxz_S_b;
  abcd[950] = 2.0E0*I_ERI_Dxz_Fxyz_Dxz_S_b;
  abcd[951] = 2.0E0*I_ERI_D2y_Fxyz_Dxz_S_b;
  abcd[952] = 2.0E0*I_ERI_Dyz_Fxyz_Dxz_S_b;
  abcd[953] = 2.0E0*I_ERI_D2z_Fxyz_Dxz_S_b;
  abcd[954] = 2.0E0*I_ERI_D2x_F3y_Dxz_S_b-2*I_ERI_D2x_Py_Dxz_S;
  abcd[955] = 2.0E0*I_ERI_Dxy_F3y_Dxz_S_b-2*I_ERI_Dxy_Py_Dxz_S;
  abcd[956] = 2.0E0*I_ERI_Dxz_F3y_Dxz_S_b-2*I_ERI_Dxz_Py_Dxz_S;
  abcd[957] = 2.0E0*I_ERI_D2y_F3y_Dxz_S_b-2*I_ERI_D2y_Py_Dxz_S;
  abcd[958] = 2.0E0*I_ERI_Dyz_F3y_Dxz_S_b-2*I_ERI_Dyz_Py_Dxz_S;
  abcd[959] = 2.0E0*I_ERI_D2z_F3y_Dxz_S_b-2*I_ERI_D2z_Py_Dxz_S;
  abcd[960] = 2.0E0*I_ERI_D2x_F2yz_Dxz_S_b-1*I_ERI_D2x_Pz_Dxz_S;
  abcd[961] = 2.0E0*I_ERI_Dxy_F2yz_Dxz_S_b-1*I_ERI_Dxy_Pz_Dxz_S;
  abcd[962] = 2.0E0*I_ERI_Dxz_F2yz_Dxz_S_b-1*I_ERI_Dxz_Pz_Dxz_S;
  abcd[963] = 2.0E0*I_ERI_D2y_F2yz_Dxz_S_b-1*I_ERI_D2y_Pz_Dxz_S;
  abcd[964] = 2.0E0*I_ERI_Dyz_F2yz_Dxz_S_b-1*I_ERI_Dyz_Pz_Dxz_S;
  abcd[965] = 2.0E0*I_ERI_D2z_F2yz_Dxz_S_b-1*I_ERI_D2z_Pz_Dxz_S;
  abcd[966] = 2.0E0*I_ERI_D2x_Fy2z_Dxz_S_b;
  abcd[967] = 2.0E0*I_ERI_Dxy_Fy2z_Dxz_S_b;
  abcd[968] = 2.0E0*I_ERI_Dxz_Fy2z_Dxz_S_b;
  abcd[969] = 2.0E0*I_ERI_D2y_Fy2z_Dxz_S_b;
  abcd[970] = 2.0E0*I_ERI_Dyz_Fy2z_Dxz_S_b;
  abcd[971] = 2.0E0*I_ERI_D2z_Fy2z_Dxz_S_b;
  abcd[972] = 2.0E0*I_ERI_D2x_F2xy_D2y_S_b;
  abcd[973] = 2.0E0*I_ERI_Dxy_F2xy_D2y_S_b;
  abcd[974] = 2.0E0*I_ERI_Dxz_F2xy_D2y_S_b;
  abcd[975] = 2.0E0*I_ERI_D2y_F2xy_D2y_S_b;
  abcd[976] = 2.0E0*I_ERI_Dyz_F2xy_D2y_S_b;
  abcd[977] = 2.0E0*I_ERI_D2z_F2xy_D2y_S_b;
  abcd[978] = 2.0E0*I_ERI_D2x_Fx2y_D2y_S_b-1*I_ERI_D2x_Px_D2y_S;
  abcd[979] = 2.0E0*I_ERI_Dxy_Fx2y_D2y_S_b-1*I_ERI_Dxy_Px_D2y_S;
  abcd[980] = 2.0E0*I_ERI_Dxz_Fx2y_D2y_S_b-1*I_ERI_Dxz_Px_D2y_S;
  abcd[981] = 2.0E0*I_ERI_D2y_Fx2y_D2y_S_b-1*I_ERI_D2y_Px_D2y_S;
  abcd[982] = 2.0E0*I_ERI_Dyz_Fx2y_D2y_S_b-1*I_ERI_Dyz_Px_D2y_S;
  abcd[983] = 2.0E0*I_ERI_D2z_Fx2y_D2y_S_b-1*I_ERI_D2z_Px_D2y_S;
  abcd[984] = 2.0E0*I_ERI_D2x_Fxyz_D2y_S_b;
  abcd[985] = 2.0E0*I_ERI_Dxy_Fxyz_D2y_S_b;
  abcd[986] = 2.0E0*I_ERI_Dxz_Fxyz_D2y_S_b;
  abcd[987] = 2.0E0*I_ERI_D2y_Fxyz_D2y_S_b;
  abcd[988] = 2.0E0*I_ERI_Dyz_Fxyz_D2y_S_b;
  abcd[989] = 2.0E0*I_ERI_D2z_Fxyz_D2y_S_b;
  abcd[990] = 2.0E0*I_ERI_D2x_F3y_D2y_S_b-2*I_ERI_D2x_Py_D2y_S;
  abcd[991] = 2.0E0*I_ERI_Dxy_F3y_D2y_S_b-2*I_ERI_Dxy_Py_D2y_S;
  abcd[992] = 2.0E0*I_ERI_Dxz_F3y_D2y_S_b-2*I_ERI_Dxz_Py_D2y_S;
  abcd[993] = 2.0E0*I_ERI_D2y_F3y_D2y_S_b-2*I_ERI_D2y_Py_D2y_S;
  abcd[994] = 2.0E0*I_ERI_Dyz_F3y_D2y_S_b-2*I_ERI_Dyz_Py_D2y_S;
  abcd[995] = 2.0E0*I_ERI_D2z_F3y_D2y_S_b-2*I_ERI_D2z_Py_D2y_S;
  abcd[996] = 2.0E0*I_ERI_D2x_F2yz_D2y_S_b-1*I_ERI_D2x_Pz_D2y_S;
  abcd[997] = 2.0E0*I_ERI_Dxy_F2yz_D2y_S_b-1*I_ERI_Dxy_Pz_D2y_S;
  abcd[998] = 2.0E0*I_ERI_Dxz_F2yz_D2y_S_b-1*I_ERI_Dxz_Pz_D2y_S;
  abcd[999] = 2.0E0*I_ERI_D2y_F2yz_D2y_S_b-1*I_ERI_D2y_Pz_D2y_S;
  abcd[1000] = 2.0E0*I_ERI_Dyz_F2yz_D2y_S_b-1*I_ERI_Dyz_Pz_D2y_S;
  abcd[1001] = 2.0E0*I_ERI_D2z_F2yz_D2y_S_b-1*I_ERI_D2z_Pz_D2y_S;
  abcd[1002] = 2.0E0*I_ERI_D2x_Fy2z_D2y_S_b;
  abcd[1003] = 2.0E0*I_ERI_Dxy_Fy2z_D2y_S_b;
  abcd[1004] = 2.0E0*I_ERI_Dxz_Fy2z_D2y_S_b;
  abcd[1005] = 2.0E0*I_ERI_D2y_Fy2z_D2y_S_b;
  abcd[1006] = 2.0E0*I_ERI_Dyz_Fy2z_D2y_S_b;
  abcd[1007] = 2.0E0*I_ERI_D2z_Fy2z_D2y_S_b;
  abcd[1008] = 2.0E0*I_ERI_D2x_F2xy_Dyz_S_b;
  abcd[1009] = 2.0E0*I_ERI_Dxy_F2xy_Dyz_S_b;
  abcd[1010] = 2.0E0*I_ERI_Dxz_F2xy_Dyz_S_b;
  abcd[1011] = 2.0E0*I_ERI_D2y_F2xy_Dyz_S_b;
  abcd[1012] = 2.0E0*I_ERI_Dyz_F2xy_Dyz_S_b;
  abcd[1013] = 2.0E0*I_ERI_D2z_F2xy_Dyz_S_b;
  abcd[1014] = 2.0E0*I_ERI_D2x_Fx2y_Dyz_S_b-1*I_ERI_D2x_Px_Dyz_S;
  abcd[1015] = 2.0E0*I_ERI_Dxy_Fx2y_Dyz_S_b-1*I_ERI_Dxy_Px_Dyz_S;
  abcd[1016] = 2.0E0*I_ERI_Dxz_Fx2y_Dyz_S_b-1*I_ERI_Dxz_Px_Dyz_S;
  abcd[1017] = 2.0E0*I_ERI_D2y_Fx2y_Dyz_S_b-1*I_ERI_D2y_Px_Dyz_S;
  abcd[1018] = 2.0E0*I_ERI_Dyz_Fx2y_Dyz_S_b-1*I_ERI_Dyz_Px_Dyz_S;
  abcd[1019] = 2.0E0*I_ERI_D2z_Fx2y_Dyz_S_b-1*I_ERI_D2z_Px_Dyz_S;
  abcd[1020] = 2.0E0*I_ERI_D2x_Fxyz_Dyz_S_b;
  abcd[1021] = 2.0E0*I_ERI_Dxy_Fxyz_Dyz_S_b;
  abcd[1022] = 2.0E0*I_ERI_Dxz_Fxyz_Dyz_S_b;
  abcd[1023] = 2.0E0*I_ERI_D2y_Fxyz_Dyz_S_b;
  abcd[1024] = 2.0E0*I_ERI_Dyz_Fxyz_Dyz_S_b;
  abcd[1025] = 2.0E0*I_ERI_D2z_Fxyz_Dyz_S_b;
  abcd[1026] = 2.0E0*I_ERI_D2x_F3y_Dyz_S_b-2*I_ERI_D2x_Py_Dyz_S;
  abcd[1027] = 2.0E0*I_ERI_Dxy_F3y_Dyz_S_b-2*I_ERI_Dxy_Py_Dyz_S;
  abcd[1028] = 2.0E0*I_ERI_Dxz_F3y_Dyz_S_b-2*I_ERI_Dxz_Py_Dyz_S;
  abcd[1029] = 2.0E0*I_ERI_D2y_F3y_Dyz_S_b-2*I_ERI_D2y_Py_Dyz_S;
  abcd[1030] = 2.0E0*I_ERI_Dyz_F3y_Dyz_S_b-2*I_ERI_Dyz_Py_Dyz_S;
  abcd[1031] = 2.0E0*I_ERI_D2z_F3y_Dyz_S_b-2*I_ERI_D2z_Py_Dyz_S;
  abcd[1032] = 2.0E0*I_ERI_D2x_F2yz_Dyz_S_b-1*I_ERI_D2x_Pz_Dyz_S;
  abcd[1033] = 2.0E0*I_ERI_Dxy_F2yz_Dyz_S_b-1*I_ERI_Dxy_Pz_Dyz_S;
  abcd[1034] = 2.0E0*I_ERI_Dxz_F2yz_Dyz_S_b-1*I_ERI_Dxz_Pz_Dyz_S;
  abcd[1035] = 2.0E0*I_ERI_D2y_F2yz_Dyz_S_b-1*I_ERI_D2y_Pz_Dyz_S;
  abcd[1036] = 2.0E0*I_ERI_Dyz_F2yz_Dyz_S_b-1*I_ERI_Dyz_Pz_Dyz_S;
  abcd[1037] = 2.0E0*I_ERI_D2z_F2yz_Dyz_S_b-1*I_ERI_D2z_Pz_Dyz_S;
  abcd[1038] = 2.0E0*I_ERI_D2x_Fy2z_Dyz_S_b;
  abcd[1039] = 2.0E0*I_ERI_Dxy_Fy2z_Dyz_S_b;
  abcd[1040] = 2.0E0*I_ERI_Dxz_Fy2z_Dyz_S_b;
  abcd[1041] = 2.0E0*I_ERI_D2y_Fy2z_Dyz_S_b;
  abcd[1042] = 2.0E0*I_ERI_Dyz_Fy2z_Dyz_S_b;
  abcd[1043] = 2.0E0*I_ERI_D2z_Fy2z_Dyz_S_b;
  abcd[1044] = 2.0E0*I_ERI_D2x_F2xy_D2z_S_b;
  abcd[1045] = 2.0E0*I_ERI_Dxy_F2xy_D2z_S_b;
  abcd[1046] = 2.0E0*I_ERI_Dxz_F2xy_D2z_S_b;
  abcd[1047] = 2.0E0*I_ERI_D2y_F2xy_D2z_S_b;
  abcd[1048] = 2.0E0*I_ERI_Dyz_F2xy_D2z_S_b;
  abcd[1049] = 2.0E0*I_ERI_D2z_F2xy_D2z_S_b;
  abcd[1050] = 2.0E0*I_ERI_D2x_Fx2y_D2z_S_b-1*I_ERI_D2x_Px_D2z_S;
  abcd[1051] = 2.0E0*I_ERI_Dxy_Fx2y_D2z_S_b-1*I_ERI_Dxy_Px_D2z_S;
  abcd[1052] = 2.0E0*I_ERI_Dxz_Fx2y_D2z_S_b-1*I_ERI_Dxz_Px_D2z_S;
  abcd[1053] = 2.0E0*I_ERI_D2y_Fx2y_D2z_S_b-1*I_ERI_D2y_Px_D2z_S;
  abcd[1054] = 2.0E0*I_ERI_Dyz_Fx2y_D2z_S_b-1*I_ERI_Dyz_Px_D2z_S;
  abcd[1055] = 2.0E0*I_ERI_D2z_Fx2y_D2z_S_b-1*I_ERI_D2z_Px_D2z_S;
  abcd[1056] = 2.0E0*I_ERI_D2x_Fxyz_D2z_S_b;
  abcd[1057] = 2.0E0*I_ERI_Dxy_Fxyz_D2z_S_b;
  abcd[1058] = 2.0E0*I_ERI_Dxz_Fxyz_D2z_S_b;
  abcd[1059] = 2.0E0*I_ERI_D2y_Fxyz_D2z_S_b;
  abcd[1060] = 2.0E0*I_ERI_Dyz_Fxyz_D2z_S_b;
  abcd[1061] = 2.0E0*I_ERI_D2z_Fxyz_D2z_S_b;
  abcd[1062] = 2.0E0*I_ERI_D2x_F3y_D2z_S_b-2*I_ERI_D2x_Py_D2z_S;
  abcd[1063] = 2.0E0*I_ERI_Dxy_F3y_D2z_S_b-2*I_ERI_Dxy_Py_D2z_S;
  abcd[1064] = 2.0E0*I_ERI_Dxz_F3y_D2z_S_b-2*I_ERI_Dxz_Py_D2z_S;
  abcd[1065] = 2.0E0*I_ERI_D2y_F3y_D2z_S_b-2*I_ERI_D2y_Py_D2z_S;
  abcd[1066] = 2.0E0*I_ERI_Dyz_F3y_D2z_S_b-2*I_ERI_Dyz_Py_D2z_S;
  abcd[1067] = 2.0E0*I_ERI_D2z_F3y_D2z_S_b-2*I_ERI_D2z_Py_D2z_S;
  abcd[1068] = 2.0E0*I_ERI_D2x_F2yz_D2z_S_b-1*I_ERI_D2x_Pz_D2z_S;
  abcd[1069] = 2.0E0*I_ERI_Dxy_F2yz_D2z_S_b-1*I_ERI_Dxy_Pz_D2z_S;
  abcd[1070] = 2.0E0*I_ERI_Dxz_F2yz_D2z_S_b-1*I_ERI_Dxz_Pz_D2z_S;
  abcd[1071] = 2.0E0*I_ERI_D2y_F2yz_D2z_S_b-1*I_ERI_D2y_Pz_D2z_S;
  abcd[1072] = 2.0E0*I_ERI_Dyz_F2yz_D2z_S_b-1*I_ERI_Dyz_Pz_D2z_S;
  abcd[1073] = 2.0E0*I_ERI_D2z_F2yz_D2z_S_b-1*I_ERI_D2z_Pz_D2z_S;
  abcd[1074] = 2.0E0*I_ERI_D2x_Fy2z_D2z_S_b;
  abcd[1075] = 2.0E0*I_ERI_Dxy_Fy2z_D2z_S_b;
  abcd[1076] = 2.0E0*I_ERI_Dxz_Fy2z_D2z_S_b;
  abcd[1077] = 2.0E0*I_ERI_D2y_Fy2z_D2z_S_b;
  abcd[1078] = 2.0E0*I_ERI_Dyz_Fy2z_D2z_S_b;
  abcd[1079] = 2.0E0*I_ERI_D2z_Fy2z_D2z_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_D_S_b
   * RHS shell quartet name: SQ_ERI_D_P_D_S
   ************************************************************/
  abcd[1080] = 2.0E0*I_ERI_D2x_F2xz_D2x_S_b;
  abcd[1081] = 2.0E0*I_ERI_Dxy_F2xz_D2x_S_b;
  abcd[1082] = 2.0E0*I_ERI_Dxz_F2xz_D2x_S_b;
  abcd[1083] = 2.0E0*I_ERI_D2y_F2xz_D2x_S_b;
  abcd[1084] = 2.0E0*I_ERI_Dyz_F2xz_D2x_S_b;
  abcd[1085] = 2.0E0*I_ERI_D2z_F2xz_D2x_S_b;
  abcd[1086] = 2.0E0*I_ERI_D2x_Fxyz_D2x_S_b;
  abcd[1087] = 2.0E0*I_ERI_Dxy_Fxyz_D2x_S_b;
  abcd[1088] = 2.0E0*I_ERI_Dxz_Fxyz_D2x_S_b;
  abcd[1089] = 2.0E0*I_ERI_D2y_Fxyz_D2x_S_b;
  abcd[1090] = 2.0E0*I_ERI_Dyz_Fxyz_D2x_S_b;
  abcd[1091] = 2.0E0*I_ERI_D2z_Fxyz_D2x_S_b;
  abcd[1092] = 2.0E0*I_ERI_D2x_Fx2z_D2x_S_b-1*I_ERI_D2x_Px_D2x_S;
  abcd[1093] = 2.0E0*I_ERI_Dxy_Fx2z_D2x_S_b-1*I_ERI_Dxy_Px_D2x_S;
  abcd[1094] = 2.0E0*I_ERI_Dxz_Fx2z_D2x_S_b-1*I_ERI_Dxz_Px_D2x_S;
  abcd[1095] = 2.0E0*I_ERI_D2y_Fx2z_D2x_S_b-1*I_ERI_D2y_Px_D2x_S;
  abcd[1096] = 2.0E0*I_ERI_Dyz_Fx2z_D2x_S_b-1*I_ERI_Dyz_Px_D2x_S;
  abcd[1097] = 2.0E0*I_ERI_D2z_Fx2z_D2x_S_b-1*I_ERI_D2z_Px_D2x_S;
  abcd[1098] = 2.0E0*I_ERI_D2x_F2yz_D2x_S_b;
  abcd[1099] = 2.0E0*I_ERI_Dxy_F2yz_D2x_S_b;
  abcd[1100] = 2.0E0*I_ERI_Dxz_F2yz_D2x_S_b;
  abcd[1101] = 2.0E0*I_ERI_D2y_F2yz_D2x_S_b;
  abcd[1102] = 2.0E0*I_ERI_Dyz_F2yz_D2x_S_b;
  abcd[1103] = 2.0E0*I_ERI_D2z_F2yz_D2x_S_b;
  abcd[1104] = 2.0E0*I_ERI_D2x_Fy2z_D2x_S_b-1*I_ERI_D2x_Py_D2x_S;
  abcd[1105] = 2.0E0*I_ERI_Dxy_Fy2z_D2x_S_b-1*I_ERI_Dxy_Py_D2x_S;
  abcd[1106] = 2.0E0*I_ERI_Dxz_Fy2z_D2x_S_b-1*I_ERI_Dxz_Py_D2x_S;
  abcd[1107] = 2.0E0*I_ERI_D2y_Fy2z_D2x_S_b-1*I_ERI_D2y_Py_D2x_S;
  abcd[1108] = 2.0E0*I_ERI_Dyz_Fy2z_D2x_S_b-1*I_ERI_Dyz_Py_D2x_S;
  abcd[1109] = 2.0E0*I_ERI_D2z_Fy2z_D2x_S_b-1*I_ERI_D2z_Py_D2x_S;
  abcd[1110] = 2.0E0*I_ERI_D2x_F3z_D2x_S_b-2*I_ERI_D2x_Pz_D2x_S;
  abcd[1111] = 2.0E0*I_ERI_Dxy_F3z_D2x_S_b-2*I_ERI_Dxy_Pz_D2x_S;
  abcd[1112] = 2.0E0*I_ERI_Dxz_F3z_D2x_S_b-2*I_ERI_Dxz_Pz_D2x_S;
  abcd[1113] = 2.0E0*I_ERI_D2y_F3z_D2x_S_b-2*I_ERI_D2y_Pz_D2x_S;
  abcd[1114] = 2.0E0*I_ERI_Dyz_F3z_D2x_S_b-2*I_ERI_Dyz_Pz_D2x_S;
  abcd[1115] = 2.0E0*I_ERI_D2z_F3z_D2x_S_b-2*I_ERI_D2z_Pz_D2x_S;
  abcd[1116] = 2.0E0*I_ERI_D2x_F2xz_Dxy_S_b;
  abcd[1117] = 2.0E0*I_ERI_Dxy_F2xz_Dxy_S_b;
  abcd[1118] = 2.0E0*I_ERI_Dxz_F2xz_Dxy_S_b;
  abcd[1119] = 2.0E0*I_ERI_D2y_F2xz_Dxy_S_b;
  abcd[1120] = 2.0E0*I_ERI_Dyz_F2xz_Dxy_S_b;
  abcd[1121] = 2.0E0*I_ERI_D2z_F2xz_Dxy_S_b;
  abcd[1122] = 2.0E0*I_ERI_D2x_Fxyz_Dxy_S_b;
  abcd[1123] = 2.0E0*I_ERI_Dxy_Fxyz_Dxy_S_b;
  abcd[1124] = 2.0E0*I_ERI_Dxz_Fxyz_Dxy_S_b;
  abcd[1125] = 2.0E0*I_ERI_D2y_Fxyz_Dxy_S_b;
  abcd[1126] = 2.0E0*I_ERI_Dyz_Fxyz_Dxy_S_b;
  abcd[1127] = 2.0E0*I_ERI_D2z_Fxyz_Dxy_S_b;
  abcd[1128] = 2.0E0*I_ERI_D2x_Fx2z_Dxy_S_b-1*I_ERI_D2x_Px_Dxy_S;
  abcd[1129] = 2.0E0*I_ERI_Dxy_Fx2z_Dxy_S_b-1*I_ERI_Dxy_Px_Dxy_S;
  abcd[1130] = 2.0E0*I_ERI_Dxz_Fx2z_Dxy_S_b-1*I_ERI_Dxz_Px_Dxy_S;
  abcd[1131] = 2.0E0*I_ERI_D2y_Fx2z_Dxy_S_b-1*I_ERI_D2y_Px_Dxy_S;
  abcd[1132] = 2.0E0*I_ERI_Dyz_Fx2z_Dxy_S_b-1*I_ERI_Dyz_Px_Dxy_S;
  abcd[1133] = 2.0E0*I_ERI_D2z_Fx2z_Dxy_S_b-1*I_ERI_D2z_Px_Dxy_S;
  abcd[1134] = 2.0E0*I_ERI_D2x_F2yz_Dxy_S_b;
  abcd[1135] = 2.0E0*I_ERI_Dxy_F2yz_Dxy_S_b;
  abcd[1136] = 2.0E0*I_ERI_Dxz_F2yz_Dxy_S_b;
  abcd[1137] = 2.0E0*I_ERI_D2y_F2yz_Dxy_S_b;
  abcd[1138] = 2.0E0*I_ERI_Dyz_F2yz_Dxy_S_b;
  abcd[1139] = 2.0E0*I_ERI_D2z_F2yz_Dxy_S_b;
  abcd[1140] = 2.0E0*I_ERI_D2x_Fy2z_Dxy_S_b-1*I_ERI_D2x_Py_Dxy_S;
  abcd[1141] = 2.0E0*I_ERI_Dxy_Fy2z_Dxy_S_b-1*I_ERI_Dxy_Py_Dxy_S;
  abcd[1142] = 2.0E0*I_ERI_Dxz_Fy2z_Dxy_S_b-1*I_ERI_Dxz_Py_Dxy_S;
  abcd[1143] = 2.0E0*I_ERI_D2y_Fy2z_Dxy_S_b-1*I_ERI_D2y_Py_Dxy_S;
  abcd[1144] = 2.0E0*I_ERI_Dyz_Fy2z_Dxy_S_b-1*I_ERI_Dyz_Py_Dxy_S;
  abcd[1145] = 2.0E0*I_ERI_D2z_Fy2z_Dxy_S_b-1*I_ERI_D2z_Py_Dxy_S;
  abcd[1146] = 2.0E0*I_ERI_D2x_F3z_Dxy_S_b-2*I_ERI_D2x_Pz_Dxy_S;
  abcd[1147] = 2.0E0*I_ERI_Dxy_F3z_Dxy_S_b-2*I_ERI_Dxy_Pz_Dxy_S;
  abcd[1148] = 2.0E0*I_ERI_Dxz_F3z_Dxy_S_b-2*I_ERI_Dxz_Pz_Dxy_S;
  abcd[1149] = 2.0E0*I_ERI_D2y_F3z_Dxy_S_b-2*I_ERI_D2y_Pz_Dxy_S;
  abcd[1150] = 2.0E0*I_ERI_Dyz_F3z_Dxy_S_b-2*I_ERI_Dyz_Pz_Dxy_S;
  abcd[1151] = 2.0E0*I_ERI_D2z_F3z_Dxy_S_b-2*I_ERI_D2z_Pz_Dxy_S;
  abcd[1152] = 2.0E0*I_ERI_D2x_F2xz_Dxz_S_b;
  abcd[1153] = 2.0E0*I_ERI_Dxy_F2xz_Dxz_S_b;
  abcd[1154] = 2.0E0*I_ERI_Dxz_F2xz_Dxz_S_b;
  abcd[1155] = 2.0E0*I_ERI_D2y_F2xz_Dxz_S_b;
  abcd[1156] = 2.0E0*I_ERI_Dyz_F2xz_Dxz_S_b;
  abcd[1157] = 2.0E0*I_ERI_D2z_F2xz_Dxz_S_b;
  abcd[1158] = 2.0E0*I_ERI_D2x_Fxyz_Dxz_S_b;
  abcd[1159] = 2.0E0*I_ERI_Dxy_Fxyz_Dxz_S_b;
  abcd[1160] = 2.0E0*I_ERI_Dxz_Fxyz_Dxz_S_b;
  abcd[1161] = 2.0E0*I_ERI_D2y_Fxyz_Dxz_S_b;
  abcd[1162] = 2.0E0*I_ERI_Dyz_Fxyz_Dxz_S_b;
  abcd[1163] = 2.0E0*I_ERI_D2z_Fxyz_Dxz_S_b;
  abcd[1164] = 2.0E0*I_ERI_D2x_Fx2z_Dxz_S_b-1*I_ERI_D2x_Px_Dxz_S;
  abcd[1165] = 2.0E0*I_ERI_Dxy_Fx2z_Dxz_S_b-1*I_ERI_Dxy_Px_Dxz_S;
  abcd[1166] = 2.0E0*I_ERI_Dxz_Fx2z_Dxz_S_b-1*I_ERI_Dxz_Px_Dxz_S;
  abcd[1167] = 2.0E0*I_ERI_D2y_Fx2z_Dxz_S_b-1*I_ERI_D2y_Px_Dxz_S;
  abcd[1168] = 2.0E0*I_ERI_Dyz_Fx2z_Dxz_S_b-1*I_ERI_Dyz_Px_Dxz_S;
  abcd[1169] = 2.0E0*I_ERI_D2z_Fx2z_Dxz_S_b-1*I_ERI_D2z_Px_Dxz_S;
  abcd[1170] = 2.0E0*I_ERI_D2x_F2yz_Dxz_S_b;
  abcd[1171] = 2.0E0*I_ERI_Dxy_F2yz_Dxz_S_b;
  abcd[1172] = 2.0E0*I_ERI_Dxz_F2yz_Dxz_S_b;
  abcd[1173] = 2.0E0*I_ERI_D2y_F2yz_Dxz_S_b;
  abcd[1174] = 2.0E0*I_ERI_Dyz_F2yz_Dxz_S_b;
  abcd[1175] = 2.0E0*I_ERI_D2z_F2yz_Dxz_S_b;
  abcd[1176] = 2.0E0*I_ERI_D2x_Fy2z_Dxz_S_b-1*I_ERI_D2x_Py_Dxz_S;
  abcd[1177] = 2.0E0*I_ERI_Dxy_Fy2z_Dxz_S_b-1*I_ERI_Dxy_Py_Dxz_S;
  abcd[1178] = 2.0E0*I_ERI_Dxz_Fy2z_Dxz_S_b-1*I_ERI_Dxz_Py_Dxz_S;
  abcd[1179] = 2.0E0*I_ERI_D2y_Fy2z_Dxz_S_b-1*I_ERI_D2y_Py_Dxz_S;
  abcd[1180] = 2.0E0*I_ERI_Dyz_Fy2z_Dxz_S_b-1*I_ERI_Dyz_Py_Dxz_S;
  abcd[1181] = 2.0E0*I_ERI_D2z_Fy2z_Dxz_S_b-1*I_ERI_D2z_Py_Dxz_S;
  abcd[1182] = 2.0E0*I_ERI_D2x_F3z_Dxz_S_b-2*I_ERI_D2x_Pz_Dxz_S;
  abcd[1183] = 2.0E0*I_ERI_Dxy_F3z_Dxz_S_b-2*I_ERI_Dxy_Pz_Dxz_S;
  abcd[1184] = 2.0E0*I_ERI_Dxz_F3z_Dxz_S_b-2*I_ERI_Dxz_Pz_Dxz_S;
  abcd[1185] = 2.0E0*I_ERI_D2y_F3z_Dxz_S_b-2*I_ERI_D2y_Pz_Dxz_S;
  abcd[1186] = 2.0E0*I_ERI_Dyz_F3z_Dxz_S_b-2*I_ERI_Dyz_Pz_Dxz_S;
  abcd[1187] = 2.0E0*I_ERI_D2z_F3z_Dxz_S_b-2*I_ERI_D2z_Pz_Dxz_S;
  abcd[1188] = 2.0E0*I_ERI_D2x_F2xz_D2y_S_b;
  abcd[1189] = 2.0E0*I_ERI_Dxy_F2xz_D2y_S_b;
  abcd[1190] = 2.0E0*I_ERI_Dxz_F2xz_D2y_S_b;
  abcd[1191] = 2.0E0*I_ERI_D2y_F2xz_D2y_S_b;
  abcd[1192] = 2.0E0*I_ERI_Dyz_F2xz_D2y_S_b;
  abcd[1193] = 2.0E0*I_ERI_D2z_F2xz_D2y_S_b;
  abcd[1194] = 2.0E0*I_ERI_D2x_Fxyz_D2y_S_b;
  abcd[1195] = 2.0E0*I_ERI_Dxy_Fxyz_D2y_S_b;
  abcd[1196] = 2.0E0*I_ERI_Dxz_Fxyz_D2y_S_b;
  abcd[1197] = 2.0E0*I_ERI_D2y_Fxyz_D2y_S_b;
  abcd[1198] = 2.0E0*I_ERI_Dyz_Fxyz_D2y_S_b;
  abcd[1199] = 2.0E0*I_ERI_D2z_Fxyz_D2y_S_b;
  abcd[1200] = 2.0E0*I_ERI_D2x_Fx2z_D2y_S_b-1*I_ERI_D2x_Px_D2y_S;
  abcd[1201] = 2.0E0*I_ERI_Dxy_Fx2z_D2y_S_b-1*I_ERI_Dxy_Px_D2y_S;
  abcd[1202] = 2.0E0*I_ERI_Dxz_Fx2z_D2y_S_b-1*I_ERI_Dxz_Px_D2y_S;
  abcd[1203] = 2.0E0*I_ERI_D2y_Fx2z_D2y_S_b-1*I_ERI_D2y_Px_D2y_S;
  abcd[1204] = 2.0E0*I_ERI_Dyz_Fx2z_D2y_S_b-1*I_ERI_Dyz_Px_D2y_S;
  abcd[1205] = 2.0E0*I_ERI_D2z_Fx2z_D2y_S_b-1*I_ERI_D2z_Px_D2y_S;
  abcd[1206] = 2.0E0*I_ERI_D2x_F2yz_D2y_S_b;
  abcd[1207] = 2.0E0*I_ERI_Dxy_F2yz_D2y_S_b;
  abcd[1208] = 2.0E0*I_ERI_Dxz_F2yz_D2y_S_b;
  abcd[1209] = 2.0E0*I_ERI_D2y_F2yz_D2y_S_b;
  abcd[1210] = 2.0E0*I_ERI_Dyz_F2yz_D2y_S_b;
  abcd[1211] = 2.0E0*I_ERI_D2z_F2yz_D2y_S_b;
  abcd[1212] = 2.0E0*I_ERI_D2x_Fy2z_D2y_S_b-1*I_ERI_D2x_Py_D2y_S;
  abcd[1213] = 2.0E0*I_ERI_Dxy_Fy2z_D2y_S_b-1*I_ERI_Dxy_Py_D2y_S;
  abcd[1214] = 2.0E0*I_ERI_Dxz_Fy2z_D2y_S_b-1*I_ERI_Dxz_Py_D2y_S;
  abcd[1215] = 2.0E0*I_ERI_D2y_Fy2z_D2y_S_b-1*I_ERI_D2y_Py_D2y_S;
  abcd[1216] = 2.0E0*I_ERI_Dyz_Fy2z_D2y_S_b-1*I_ERI_Dyz_Py_D2y_S;
  abcd[1217] = 2.0E0*I_ERI_D2z_Fy2z_D2y_S_b-1*I_ERI_D2z_Py_D2y_S;
  abcd[1218] = 2.0E0*I_ERI_D2x_F3z_D2y_S_b-2*I_ERI_D2x_Pz_D2y_S;
  abcd[1219] = 2.0E0*I_ERI_Dxy_F3z_D2y_S_b-2*I_ERI_Dxy_Pz_D2y_S;
  abcd[1220] = 2.0E0*I_ERI_Dxz_F3z_D2y_S_b-2*I_ERI_Dxz_Pz_D2y_S;
  abcd[1221] = 2.0E0*I_ERI_D2y_F3z_D2y_S_b-2*I_ERI_D2y_Pz_D2y_S;
  abcd[1222] = 2.0E0*I_ERI_Dyz_F3z_D2y_S_b-2*I_ERI_Dyz_Pz_D2y_S;
  abcd[1223] = 2.0E0*I_ERI_D2z_F3z_D2y_S_b-2*I_ERI_D2z_Pz_D2y_S;
  abcd[1224] = 2.0E0*I_ERI_D2x_F2xz_Dyz_S_b;
  abcd[1225] = 2.0E0*I_ERI_Dxy_F2xz_Dyz_S_b;
  abcd[1226] = 2.0E0*I_ERI_Dxz_F2xz_Dyz_S_b;
  abcd[1227] = 2.0E0*I_ERI_D2y_F2xz_Dyz_S_b;
  abcd[1228] = 2.0E0*I_ERI_Dyz_F2xz_Dyz_S_b;
  abcd[1229] = 2.0E0*I_ERI_D2z_F2xz_Dyz_S_b;
  abcd[1230] = 2.0E0*I_ERI_D2x_Fxyz_Dyz_S_b;
  abcd[1231] = 2.0E0*I_ERI_Dxy_Fxyz_Dyz_S_b;
  abcd[1232] = 2.0E0*I_ERI_Dxz_Fxyz_Dyz_S_b;
  abcd[1233] = 2.0E0*I_ERI_D2y_Fxyz_Dyz_S_b;
  abcd[1234] = 2.0E0*I_ERI_Dyz_Fxyz_Dyz_S_b;
  abcd[1235] = 2.0E0*I_ERI_D2z_Fxyz_Dyz_S_b;
  abcd[1236] = 2.0E0*I_ERI_D2x_Fx2z_Dyz_S_b-1*I_ERI_D2x_Px_Dyz_S;
  abcd[1237] = 2.0E0*I_ERI_Dxy_Fx2z_Dyz_S_b-1*I_ERI_Dxy_Px_Dyz_S;
  abcd[1238] = 2.0E0*I_ERI_Dxz_Fx2z_Dyz_S_b-1*I_ERI_Dxz_Px_Dyz_S;
  abcd[1239] = 2.0E0*I_ERI_D2y_Fx2z_Dyz_S_b-1*I_ERI_D2y_Px_Dyz_S;
  abcd[1240] = 2.0E0*I_ERI_Dyz_Fx2z_Dyz_S_b-1*I_ERI_Dyz_Px_Dyz_S;
  abcd[1241] = 2.0E0*I_ERI_D2z_Fx2z_Dyz_S_b-1*I_ERI_D2z_Px_Dyz_S;
  abcd[1242] = 2.0E0*I_ERI_D2x_F2yz_Dyz_S_b;
  abcd[1243] = 2.0E0*I_ERI_Dxy_F2yz_Dyz_S_b;
  abcd[1244] = 2.0E0*I_ERI_Dxz_F2yz_Dyz_S_b;
  abcd[1245] = 2.0E0*I_ERI_D2y_F2yz_Dyz_S_b;
  abcd[1246] = 2.0E0*I_ERI_Dyz_F2yz_Dyz_S_b;
  abcd[1247] = 2.0E0*I_ERI_D2z_F2yz_Dyz_S_b;
  abcd[1248] = 2.0E0*I_ERI_D2x_Fy2z_Dyz_S_b-1*I_ERI_D2x_Py_Dyz_S;
  abcd[1249] = 2.0E0*I_ERI_Dxy_Fy2z_Dyz_S_b-1*I_ERI_Dxy_Py_Dyz_S;
  abcd[1250] = 2.0E0*I_ERI_Dxz_Fy2z_Dyz_S_b-1*I_ERI_Dxz_Py_Dyz_S;
  abcd[1251] = 2.0E0*I_ERI_D2y_Fy2z_Dyz_S_b-1*I_ERI_D2y_Py_Dyz_S;
  abcd[1252] = 2.0E0*I_ERI_Dyz_Fy2z_Dyz_S_b-1*I_ERI_Dyz_Py_Dyz_S;
  abcd[1253] = 2.0E0*I_ERI_D2z_Fy2z_Dyz_S_b-1*I_ERI_D2z_Py_Dyz_S;
  abcd[1254] = 2.0E0*I_ERI_D2x_F3z_Dyz_S_b-2*I_ERI_D2x_Pz_Dyz_S;
  abcd[1255] = 2.0E0*I_ERI_Dxy_F3z_Dyz_S_b-2*I_ERI_Dxy_Pz_Dyz_S;
  abcd[1256] = 2.0E0*I_ERI_Dxz_F3z_Dyz_S_b-2*I_ERI_Dxz_Pz_Dyz_S;
  abcd[1257] = 2.0E0*I_ERI_D2y_F3z_Dyz_S_b-2*I_ERI_D2y_Pz_Dyz_S;
  abcd[1258] = 2.0E0*I_ERI_Dyz_F3z_Dyz_S_b-2*I_ERI_Dyz_Pz_Dyz_S;
  abcd[1259] = 2.0E0*I_ERI_D2z_F3z_Dyz_S_b-2*I_ERI_D2z_Pz_Dyz_S;
  abcd[1260] = 2.0E0*I_ERI_D2x_F2xz_D2z_S_b;
  abcd[1261] = 2.0E0*I_ERI_Dxy_F2xz_D2z_S_b;
  abcd[1262] = 2.0E0*I_ERI_Dxz_F2xz_D2z_S_b;
  abcd[1263] = 2.0E0*I_ERI_D2y_F2xz_D2z_S_b;
  abcd[1264] = 2.0E0*I_ERI_Dyz_F2xz_D2z_S_b;
  abcd[1265] = 2.0E0*I_ERI_D2z_F2xz_D2z_S_b;
  abcd[1266] = 2.0E0*I_ERI_D2x_Fxyz_D2z_S_b;
  abcd[1267] = 2.0E0*I_ERI_Dxy_Fxyz_D2z_S_b;
  abcd[1268] = 2.0E0*I_ERI_Dxz_Fxyz_D2z_S_b;
  abcd[1269] = 2.0E0*I_ERI_D2y_Fxyz_D2z_S_b;
  abcd[1270] = 2.0E0*I_ERI_Dyz_Fxyz_D2z_S_b;
  abcd[1271] = 2.0E0*I_ERI_D2z_Fxyz_D2z_S_b;
  abcd[1272] = 2.0E0*I_ERI_D2x_Fx2z_D2z_S_b-1*I_ERI_D2x_Px_D2z_S;
  abcd[1273] = 2.0E0*I_ERI_Dxy_Fx2z_D2z_S_b-1*I_ERI_Dxy_Px_D2z_S;
  abcd[1274] = 2.0E0*I_ERI_Dxz_Fx2z_D2z_S_b-1*I_ERI_Dxz_Px_D2z_S;
  abcd[1275] = 2.0E0*I_ERI_D2y_Fx2z_D2z_S_b-1*I_ERI_D2y_Px_D2z_S;
  abcd[1276] = 2.0E0*I_ERI_Dyz_Fx2z_D2z_S_b-1*I_ERI_Dyz_Px_D2z_S;
  abcd[1277] = 2.0E0*I_ERI_D2z_Fx2z_D2z_S_b-1*I_ERI_D2z_Px_D2z_S;
  abcd[1278] = 2.0E0*I_ERI_D2x_F2yz_D2z_S_b;
  abcd[1279] = 2.0E0*I_ERI_Dxy_F2yz_D2z_S_b;
  abcd[1280] = 2.0E0*I_ERI_Dxz_F2yz_D2z_S_b;
  abcd[1281] = 2.0E0*I_ERI_D2y_F2yz_D2z_S_b;
  abcd[1282] = 2.0E0*I_ERI_Dyz_F2yz_D2z_S_b;
  abcd[1283] = 2.0E0*I_ERI_D2z_F2yz_D2z_S_b;
  abcd[1284] = 2.0E0*I_ERI_D2x_Fy2z_D2z_S_b-1*I_ERI_D2x_Py_D2z_S;
  abcd[1285] = 2.0E0*I_ERI_Dxy_Fy2z_D2z_S_b-1*I_ERI_Dxy_Py_D2z_S;
  abcd[1286] = 2.0E0*I_ERI_Dxz_Fy2z_D2z_S_b-1*I_ERI_Dxz_Py_D2z_S;
  abcd[1287] = 2.0E0*I_ERI_D2y_Fy2z_D2z_S_b-1*I_ERI_D2y_Py_D2z_S;
  abcd[1288] = 2.0E0*I_ERI_Dyz_Fy2z_D2z_S_b-1*I_ERI_Dyz_Py_D2z_S;
  abcd[1289] = 2.0E0*I_ERI_D2z_Fy2z_D2z_S_b-1*I_ERI_D2z_Py_D2z_S;
  abcd[1290] = 2.0E0*I_ERI_D2x_F3z_D2z_S_b-2*I_ERI_D2x_Pz_D2z_S;
  abcd[1291] = 2.0E0*I_ERI_Dxy_F3z_D2z_S_b-2*I_ERI_Dxy_Pz_D2z_S;
  abcd[1292] = 2.0E0*I_ERI_Dxz_F3z_D2z_S_b-2*I_ERI_Dxz_Pz_D2z_S;
  abcd[1293] = 2.0E0*I_ERI_D2y_F3z_D2z_S_b-2*I_ERI_D2y_Pz_D2z_S;
  abcd[1294] = 2.0E0*I_ERI_Dyz_F3z_D2z_S_b-2*I_ERI_Dyz_Pz_D2z_S;
  abcd[1295] = 2.0E0*I_ERI_D2z_F3z_D2z_S_b-2*I_ERI_D2z_Pz_D2z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_c
   * RHS shell quartet name: SQ_ERI_D_D_P_S
   ************************************************************/
  abcd[1296] = 2.0E0*I_ERI_D2x_D2x_F3x_S_c-2*I_ERI_D2x_D2x_Px_S;
  abcd[1297] = 2.0E0*I_ERI_Dxy_D2x_F3x_S_c-2*I_ERI_Dxy_D2x_Px_S;
  abcd[1298] = 2.0E0*I_ERI_Dxz_D2x_F3x_S_c-2*I_ERI_Dxz_D2x_Px_S;
  abcd[1299] = 2.0E0*I_ERI_D2y_D2x_F3x_S_c-2*I_ERI_D2y_D2x_Px_S;
  abcd[1300] = 2.0E0*I_ERI_Dyz_D2x_F3x_S_c-2*I_ERI_Dyz_D2x_Px_S;
  abcd[1301] = 2.0E0*I_ERI_D2z_D2x_F3x_S_c-2*I_ERI_D2z_D2x_Px_S;
  abcd[1302] = 2.0E0*I_ERI_D2x_Dxy_F3x_S_c-2*I_ERI_D2x_Dxy_Px_S;
  abcd[1303] = 2.0E0*I_ERI_Dxy_Dxy_F3x_S_c-2*I_ERI_Dxy_Dxy_Px_S;
  abcd[1304] = 2.0E0*I_ERI_Dxz_Dxy_F3x_S_c-2*I_ERI_Dxz_Dxy_Px_S;
  abcd[1305] = 2.0E0*I_ERI_D2y_Dxy_F3x_S_c-2*I_ERI_D2y_Dxy_Px_S;
  abcd[1306] = 2.0E0*I_ERI_Dyz_Dxy_F3x_S_c-2*I_ERI_Dyz_Dxy_Px_S;
  abcd[1307] = 2.0E0*I_ERI_D2z_Dxy_F3x_S_c-2*I_ERI_D2z_Dxy_Px_S;
  abcd[1308] = 2.0E0*I_ERI_D2x_Dxz_F3x_S_c-2*I_ERI_D2x_Dxz_Px_S;
  abcd[1309] = 2.0E0*I_ERI_Dxy_Dxz_F3x_S_c-2*I_ERI_Dxy_Dxz_Px_S;
  abcd[1310] = 2.0E0*I_ERI_Dxz_Dxz_F3x_S_c-2*I_ERI_Dxz_Dxz_Px_S;
  abcd[1311] = 2.0E0*I_ERI_D2y_Dxz_F3x_S_c-2*I_ERI_D2y_Dxz_Px_S;
  abcd[1312] = 2.0E0*I_ERI_Dyz_Dxz_F3x_S_c-2*I_ERI_Dyz_Dxz_Px_S;
  abcd[1313] = 2.0E0*I_ERI_D2z_Dxz_F3x_S_c-2*I_ERI_D2z_Dxz_Px_S;
  abcd[1314] = 2.0E0*I_ERI_D2x_D2y_F3x_S_c-2*I_ERI_D2x_D2y_Px_S;
  abcd[1315] = 2.0E0*I_ERI_Dxy_D2y_F3x_S_c-2*I_ERI_Dxy_D2y_Px_S;
  abcd[1316] = 2.0E0*I_ERI_Dxz_D2y_F3x_S_c-2*I_ERI_Dxz_D2y_Px_S;
  abcd[1317] = 2.0E0*I_ERI_D2y_D2y_F3x_S_c-2*I_ERI_D2y_D2y_Px_S;
  abcd[1318] = 2.0E0*I_ERI_Dyz_D2y_F3x_S_c-2*I_ERI_Dyz_D2y_Px_S;
  abcd[1319] = 2.0E0*I_ERI_D2z_D2y_F3x_S_c-2*I_ERI_D2z_D2y_Px_S;
  abcd[1320] = 2.0E0*I_ERI_D2x_Dyz_F3x_S_c-2*I_ERI_D2x_Dyz_Px_S;
  abcd[1321] = 2.0E0*I_ERI_Dxy_Dyz_F3x_S_c-2*I_ERI_Dxy_Dyz_Px_S;
  abcd[1322] = 2.0E0*I_ERI_Dxz_Dyz_F3x_S_c-2*I_ERI_Dxz_Dyz_Px_S;
  abcd[1323] = 2.0E0*I_ERI_D2y_Dyz_F3x_S_c-2*I_ERI_D2y_Dyz_Px_S;
  abcd[1324] = 2.0E0*I_ERI_Dyz_Dyz_F3x_S_c-2*I_ERI_Dyz_Dyz_Px_S;
  abcd[1325] = 2.0E0*I_ERI_D2z_Dyz_F3x_S_c-2*I_ERI_D2z_Dyz_Px_S;
  abcd[1326] = 2.0E0*I_ERI_D2x_D2z_F3x_S_c-2*I_ERI_D2x_D2z_Px_S;
  abcd[1327] = 2.0E0*I_ERI_Dxy_D2z_F3x_S_c-2*I_ERI_Dxy_D2z_Px_S;
  abcd[1328] = 2.0E0*I_ERI_Dxz_D2z_F3x_S_c-2*I_ERI_Dxz_D2z_Px_S;
  abcd[1329] = 2.0E0*I_ERI_D2y_D2z_F3x_S_c-2*I_ERI_D2y_D2z_Px_S;
  abcd[1330] = 2.0E0*I_ERI_Dyz_D2z_F3x_S_c-2*I_ERI_Dyz_D2z_Px_S;
  abcd[1331] = 2.0E0*I_ERI_D2z_D2z_F3x_S_c-2*I_ERI_D2z_D2z_Px_S;
  abcd[1332] = 2.0E0*I_ERI_D2x_D2x_F2xy_S_c-1*I_ERI_D2x_D2x_Py_S;
  abcd[1333] = 2.0E0*I_ERI_Dxy_D2x_F2xy_S_c-1*I_ERI_Dxy_D2x_Py_S;
  abcd[1334] = 2.0E0*I_ERI_Dxz_D2x_F2xy_S_c-1*I_ERI_Dxz_D2x_Py_S;
  abcd[1335] = 2.0E0*I_ERI_D2y_D2x_F2xy_S_c-1*I_ERI_D2y_D2x_Py_S;
  abcd[1336] = 2.0E0*I_ERI_Dyz_D2x_F2xy_S_c-1*I_ERI_Dyz_D2x_Py_S;
  abcd[1337] = 2.0E0*I_ERI_D2z_D2x_F2xy_S_c-1*I_ERI_D2z_D2x_Py_S;
  abcd[1338] = 2.0E0*I_ERI_D2x_Dxy_F2xy_S_c-1*I_ERI_D2x_Dxy_Py_S;
  abcd[1339] = 2.0E0*I_ERI_Dxy_Dxy_F2xy_S_c-1*I_ERI_Dxy_Dxy_Py_S;
  abcd[1340] = 2.0E0*I_ERI_Dxz_Dxy_F2xy_S_c-1*I_ERI_Dxz_Dxy_Py_S;
  abcd[1341] = 2.0E0*I_ERI_D2y_Dxy_F2xy_S_c-1*I_ERI_D2y_Dxy_Py_S;
  abcd[1342] = 2.0E0*I_ERI_Dyz_Dxy_F2xy_S_c-1*I_ERI_Dyz_Dxy_Py_S;
  abcd[1343] = 2.0E0*I_ERI_D2z_Dxy_F2xy_S_c-1*I_ERI_D2z_Dxy_Py_S;
  abcd[1344] = 2.0E0*I_ERI_D2x_Dxz_F2xy_S_c-1*I_ERI_D2x_Dxz_Py_S;
  abcd[1345] = 2.0E0*I_ERI_Dxy_Dxz_F2xy_S_c-1*I_ERI_Dxy_Dxz_Py_S;
  abcd[1346] = 2.0E0*I_ERI_Dxz_Dxz_F2xy_S_c-1*I_ERI_Dxz_Dxz_Py_S;
  abcd[1347] = 2.0E0*I_ERI_D2y_Dxz_F2xy_S_c-1*I_ERI_D2y_Dxz_Py_S;
  abcd[1348] = 2.0E0*I_ERI_Dyz_Dxz_F2xy_S_c-1*I_ERI_Dyz_Dxz_Py_S;
  abcd[1349] = 2.0E0*I_ERI_D2z_Dxz_F2xy_S_c-1*I_ERI_D2z_Dxz_Py_S;
  abcd[1350] = 2.0E0*I_ERI_D2x_D2y_F2xy_S_c-1*I_ERI_D2x_D2y_Py_S;
  abcd[1351] = 2.0E0*I_ERI_Dxy_D2y_F2xy_S_c-1*I_ERI_Dxy_D2y_Py_S;
  abcd[1352] = 2.0E0*I_ERI_Dxz_D2y_F2xy_S_c-1*I_ERI_Dxz_D2y_Py_S;
  abcd[1353] = 2.0E0*I_ERI_D2y_D2y_F2xy_S_c-1*I_ERI_D2y_D2y_Py_S;
  abcd[1354] = 2.0E0*I_ERI_Dyz_D2y_F2xy_S_c-1*I_ERI_Dyz_D2y_Py_S;
  abcd[1355] = 2.0E0*I_ERI_D2z_D2y_F2xy_S_c-1*I_ERI_D2z_D2y_Py_S;
  abcd[1356] = 2.0E0*I_ERI_D2x_Dyz_F2xy_S_c-1*I_ERI_D2x_Dyz_Py_S;
  abcd[1357] = 2.0E0*I_ERI_Dxy_Dyz_F2xy_S_c-1*I_ERI_Dxy_Dyz_Py_S;
  abcd[1358] = 2.0E0*I_ERI_Dxz_Dyz_F2xy_S_c-1*I_ERI_Dxz_Dyz_Py_S;
  abcd[1359] = 2.0E0*I_ERI_D2y_Dyz_F2xy_S_c-1*I_ERI_D2y_Dyz_Py_S;
  abcd[1360] = 2.0E0*I_ERI_Dyz_Dyz_F2xy_S_c-1*I_ERI_Dyz_Dyz_Py_S;
  abcd[1361] = 2.0E0*I_ERI_D2z_Dyz_F2xy_S_c-1*I_ERI_D2z_Dyz_Py_S;
  abcd[1362] = 2.0E0*I_ERI_D2x_D2z_F2xy_S_c-1*I_ERI_D2x_D2z_Py_S;
  abcd[1363] = 2.0E0*I_ERI_Dxy_D2z_F2xy_S_c-1*I_ERI_Dxy_D2z_Py_S;
  abcd[1364] = 2.0E0*I_ERI_Dxz_D2z_F2xy_S_c-1*I_ERI_Dxz_D2z_Py_S;
  abcd[1365] = 2.0E0*I_ERI_D2y_D2z_F2xy_S_c-1*I_ERI_D2y_D2z_Py_S;
  abcd[1366] = 2.0E0*I_ERI_Dyz_D2z_F2xy_S_c-1*I_ERI_Dyz_D2z_Py_S;
  abcd[1367] = 2.0E0*I_ERI_D2z_D2z_F2xy_S_c-1*I_ERI_D2z_D2z_Py_S;
  abcd[1368] = 2.0E0*I_ERI_D2x_D2x_F2xz_S_c-1*I_ERI_D2x_D2x_Pz_S;
  abcd[1369] = 2.0E0*I_ERI_Dxy_D2x_F2xz_S_c-1*I_ERI_Dxy_D2x_Pz_S;
  abcd[1370] = 2.0E0*I_ERI_Dxz_D2x_F2xz_S_c-1*I_ERI_Dxz_D2x_Pz_S;
  abcd[1371] = 2.0E0*I_ERI_D2y_D2x_F2xz_S_c-1*I_ERI_D2y_D2x_Pz_S;
  abcd[1372] = 2.0E0*I_ERI_Dyz_D2x_F2xz_S_c-1*I_ERI_Dyz_D2x_Pz_S;
  abcd[1373] = 2.0E0*I_ERI_D2z_D2x_F2xz_S_c-1*I_ERI_D2z_D2x_Pz_S;
  abcd[1374] = 2.0E0*I_ERI_D2x_Dxy_F2xz_S_c-1*I_ERI_D2x_Dxy_Pz_S;
  abcd[1375] = 2.0E0*I_ERI_Dxy_Dxy_F2xz_S_c-1*I_ERI_Dxy_Dxy_Pz_S;
  abcd[1376] = 2.0E0*I_ERI_Dxz_Dxy_F2xz_S_c-1*I_ERI_Dxz_Dxy_Pz_S;
  abcd[1377] = 2.0E0*I_ERI_D2y_Dxy_F2xz_S_c-1*I_ERI_D2y_Dxy_Pz_S;
  abcd[1378] = 2.0E0*I_ERI_Dyz_Dxy_F2xz_S_c-1*I_ERI_Dyz_Dxy_Pz_S;
  abcd[1379] = 2.0E0*I_ERI_D2z_Dxy_F2xz_S_c-1*I_ERI_D2z_Dxy_Pz_S;
  abcd[1380] = 2.0E0*I_ERI_D2x_Dxz_F2xz_S_c-1*I_ERI_D2x_Dxz_Pz_S;
  abcd[1381] = 2.0E0*I_ERI_Dxy_Dxz_F2xz_S_c-1*I_ERI_Dxy_Dxz_Pz_S;
  abcd[1382] = 2.0E0*I_ERI_Dxz_Dxz_F2xz_S_c-1*I_ERI_Dxz_Dxz_Pz_S;
  abcd[1383] = 2.0E0*I_ERI_D2y_Dxz_F2xz_S_c-1*I_ERI_D2y_Dxz_Pz_S;
  abcd[1384] = 2.0E0*I_ERI_Dyz_Dxz_F2xz_S_c-1*I_ERI_Dyz_Dxz_Pz_S;
  abcd[1385] = 2.0E0*I_ERI_D2z_Dxz_F2xz_S_c-1*I_ERI_D2z_Dxz_Pz_S;
  abcd[1386] = 2.0E0*I_ERI_D2x_D2y_F2xz_S_c-1*I_ERI_D2x_D2y_Pz_S;
  abcd[1387] = 2.0E0*I_ERI_Dxy_D2y_F2xz_S_c-1*I_ERI_Dxy_D2y_Pz_S;
  abcd[1388] = 2.0E0*I_ERI_Dxz_D2y_F2xz_S_c-1*I_ERI_Dxz_D2y_Pz_S;
  abcd[1389] = 2.0E0*I_ERI_D2y_D2y_F2xz_S_c-1*I_ERI_D2y_D2y_Pz_S;
  abcd[1390] = 2.0E0*I_ERI_Dyz_D2y_F2xz_S_c-1*I_ERI_Dyz_D2y_Pz_S;
  abcd[1391] = 2.0E0*I_ERI_D2z_D2y_F2xz_S_c-1*I_ERI_D2z_D2y_Pz_S;
  abcd[1392] = 2.0E0*I_ERI_D2x_Dyz_F2xz_S_c-1*I_ERI_D2x_Dyz_Pz_S;
  abcd[1393] = 2.0E0*I_ERI_Dxy_Dyz_F2xz_S_c-1*I_ERI_Dxy_Dyz_Pz_S;
  abcd[1394] = 2.0E0*I_ERI_Dxz_Dyz_F2xz_S_c-1*I_ERI_Dxz_Dyz_Pz_S;
  abcd[1395] = 2.0E0*I_ERI_D2y_Dyz_F2xz_S_c-1*I_ERI_D2y_Dyz_Pz_S;
  abcd[1396] = 2.0E0*I_ERI_Dyz_Dyz_F2xz_S_c-1*I_ERI_Dyz_Dyz_Pz_S;
  abcd[1397] = 2.0E0*I_ERI_D2z_Dyz_F2xz_S_c-1*I_ERI_D2z_Dyz_Pz_S;
  abcd[1398] = 2.0E0*I_ERI_D2x_D2z_F2xz_S_c-1*I_ERI_D2x_D2z_Pz_S;
  abcd[1399] = 2.0E0*I_ERI_Dxy_D2z_F2xz_S_c-1*I_ERI_Dxy_D2z_Pz_S;
  abcd[1400] = 2.0E0*I_ERI_Dxz_D2z_F2xz_S_c-1*I_ERI_Dxz_D2z_Pz_S;
  abcd[1401] = 2.0E0*I_ERI_D2y_D2z_F2xz_S_c-1*I_ERI_D2y_D2z_Pz_S;
  abcd[1402] = 2.0E0*I_ERI_Dyz_D2z_F2xz_S_c-1*I_ERI_Dyz_D2z_Pz_S;
  abcd[1403] = 2.0E0*I_ERI_D2z_D2z_F2xz_S_c-1*I_ERI_D2z_D2z_Pz_S;
  abcd[1404] = 2.0E0*I_ERI_D2x_D2x_Fx2y_S_c;
  abcd[1405] = 2.0E0*I_ERI_Dxy_D2x_Fx2y_S_c;
  abcd[1406] = 2.0E0*I_ERI_Dxz_D2x_Fx2y_S_c;
  abcd[1407] = 2.0E0*I_ERI_D2y_D2x_Fx2y_S_c;
  abcd[1408] = 2.0E0*I_ERI_Dyz_D2x_Fx2y_S_c;
  abcd[1409] = 2.0E0*I_ERI_D2z_D2x_Fx2y_S_c;
  abcd[1410] = 2.0E0*I_ERI_D2x_Dxy_Fx2y_S_c;
  abcd[1411] = 2.0E0*I_ERI_Dxy_Dxy_Fx2y_S_c;
  abcd[1412] = 2.0E0*I_ERI_Dxz_Dxy_Fx2y_S_c;
  abcd[1413] = 2.0E0*I_ERI_D2y_Dxy_Fx2y_S_c;
  abcd[1414] = 2.0E0*I_ERI_Dyz_Dxy_Fx2y_S_c;
  abcd[1415] = 2.0E0*I_ERI_D2z_Dxy_Fx2y_S_c;
  abcd[1416] = 2.0E0*I_ERI_D2x_Dxz_Fx2y_S_c;
  abcd[1417] = 2.0E0*I_ERI_Dxy_Dxz_Fx2y_S_c;
  abcd[1418] = 2.0E0*I_ERI_Dxz_Dxz_Fx2y_S_c;
  abcd[1419] = 2.0E0*I_ERI_D2y_Dxz_Fx2y_S_c;
  abcd[1420] = 2.0E0*I_ERI_Dyz_Dxz_Fx2y_S_c;
  abcd[1421] = 2.0E0*I_ERI_D2z_Dxz_Fx2y_S_c;
  abcd[1422] = 2.0E0*I_ERI_D2x_D2y_Fx2y_S_c;
  abcd[1423] = 2.0E0*I_ERI_Dxy_D2y_Fx2y_S_c;
  abcd[1424] = 2.0E0*I_ERI_Dxz_D2y_Fx2y_S_c;
  abcd[1425] = 2.0E0*I_ERI_D2y_D2y_Fx2y_S_c;
  abcd[1426] = 2.0E0*I_ERI_Dyz_D2y_Fx2y_S_c;
  abcd[1427] = 2.0E0*I_ERI_D2z_D2y_Fx2y_S_c;
  abcd[1428] = 2.0E0*I_ERI_D2x_Dyz_Fx2y_S_c;
  abcd[1429] = 2.0E0*I_ERI_Dxy_Dyz_Fx2y_S_c;
  abcd[1430] = 2.0E0*I_ERI_Dxz_Dyz_Fx2y_S_c;
  abcd[1431] = 2.0E0*I_ERI_D2y_Dyz_Fx2y_S_c;
  abcd[1432] = 2.0E0*I_ERI_Dyz_Dyz_Fx2y_S_c;
  abcd[1433] = 2.0E0*I_ERI_D2z_Dyz_Fx2y_S_c;
  abcd[1434] = 2.0E0*I_ERI_D2x_D2z_Fx2y_S_c;
  abcd[1435] = 2.0E0*I_ERI_Dxy_D2z_Fx2y_S_c;
  abcd[1436] = 2.0E0*I_ERI_Dxz_D2z_Fx2y_S_c;
  abcd[1437] = 2.0E0*I_ERI_D2y_D2z_Fx2y_S_c;
  abcd[1438] = 2.0E0*I_ERI_Dyz_D2z_Fx2y_S_c;
  abcd[1439] = 2.0E0*I_ERI_D2z_D2z_Fx2y_S_c;
  abcd[1440] = 2.0E0*I_ERI_D2x_D2x_Fxyz_S_c;
  abcd[1441] = 2.0E0*I_ERI_Dxy_D2x_Fxyz_S_c;
  abcd[1442] = 2.0E0*I_ERI_Dxz_D2x_Fxyz_S_c;
  abcd[1443] = 2.0E0*I_ERI_D2y_D2x_Fxyz_S_c;
  abcd[1444] = 2.0E0*I_ERI_Dyz_D2x_Fxyz_S_c;
  abcd[1445] = 2.0E0*I_ERI_D2z_D2x_Fxyz_S_c;
  abcd[1446] = 2.0E0*I_ERI_D2x_Dxy_Fxyz_S_c;
  abcd[1447] = 2.0E0*I_ERI_Dxy_Dxy_Fxyz_S_c;
  abcd[1448] = 2.0E0*I_ERI_Dxz_Dxy_Fxyz_S_c;
  abcd[1449] = 2.0E0*I_ERI_D2y_Dxy_Fxyz_S_c;
  abcd[1450] = 2.0E0*I_ERI_Dyz_Dxy_Fxyz_S_c;
  abcd[1451] = 2.0E0*I_ERI_D2z_Dxy_Fxyz_S_c;
  abcd[1452] = 2.0E0*I_ERI_D2x_Dxz_Fxyz_S_c;
  abcd[1453] = 2.0E0*I_ERI_Dxy_Dxz_Fxyz_S_c;
  abcd[1454] = 2.0E0*I_ERI_Dxz_Dxz_Fxyz_S_c;
  abcd[1455] = 2.0E0*I_ERI_D2y_Dxz_Fxyz_S_c;
  abcd[1456] = 2.0E0*I_ERI_Dyz_Dxz_Fxyz_S_c;
  abcd[1457] = 2.0E0*I_ERI_D2z_Dxz_Fxyz_S_c;
  abcd[1458] = 2.0E0*I_ERI_D2x_D2y_Fxyz_S_c;
  abcd[1459] = 2.0E0*I_ERI_Dxy_D2y_Fxyz_S_c;
  abcd[1460] = 2.0E0*I_ERI_Dxz_D2y_Fxyz_S_c;
  abcd[1461] = 2.0E0*I_ERI_D2y_D2y_Fxyz_S_c;
  abcd[1462] = 2.0E0*I_ERI_Dyz_D2y_Fxyz_S_c;
  abcd[1463] = 2.0E0*I_ERI_D2z_D2y_Fxyz_S_c;
  abcd[1464] = 2.0E0*I_ERI_D2x_Dyz_Fxyz_S_c;
  abcd[1465] = 2.0E0*I_ERI_Dxy_Dyz_Fxyz_S_c;
  abcd[1466] = 2.0E0*I_ERI_Dxz_Dyz_Fxyz_S_c;
  abcd[1467] = 2.0E0*I_ERI_D2y_Dyz_Fxyz_S_c;
  abcd[1468] = 2.0E0*I_ERI_Dyz_Dyz_Fxyz_S_c;
  abcd[1469] = 2.0E0*I_ERI_D2z_Dyz_Fxyz_S_c;
  abcd[1470] = 2.0E0*I_ERI_D2x_D2z_Fxyz_S_c;
  abcd[1471] = 2.0E0*I_ERI_Dxy_D2z_Fxyz_S_c;
  abcd[1472] = 2.0E0*I_ERI_Dxz_D2z_Fxyz_S_c;
  abcd[1473] = 2.0E0*I_ERI_D2y_D2z_Fxyz_S_c;
  abcd[1474] = 2.0E0*I_ERI_Dyz_D2z_Fxyz_S_c;
  abcd[1475] = 2.0E0*I_ERI_D2z_D2z_Fxyz_S_c;
  abcd[1476] = 2.0E0*I_ERI_D2x_D2x_Fx2z_S_c;
  abcd[1477] = 2.0E0*I_ERI_Dxy_D2x_Fx2z_S_c;
  abcd[1478] = 2.0E0*I_ERI_Dxz_D2x_Fx2z_S_c;
  abcd[1479] = 2.0E0*I_ERI_D2y_D2x_Fx2z_S_c;
  abcd[1480] = 2.0E0*I_ERI_Dyz_D2x_Fx2z_S_c;
  abcd[1481] = 2.0E0*I_ERI_D2z_D2x_Fx2z_S_c;
  abcd[1482] = 2.0E0*I_ERI_D2x_Dxy_Fx2z_S_c;
  abcd[1483] = 2.0E0*I_ERI_Dxy_Dxy_Fx2z_S_c;
  abcd[1484] = 2.0E0*I_ERI_Dxz_Dxy_Fx2z_S_c;
  abcd[1485] = 2.0E0*I_ERI_D2y_Dxy_Fx2z_S_c;
  abcd[1486] = 2.0E0*I_ERI_Dyz_Dxy_Fx2z_S_c;
  abcd[1487] = 2.0E0*I_ERI_D2z_Dxy_Fx2z_S_c;
  abcd[1488] = 2.0E0*I_ERI_D2x_Dxz_Fx2z_S_c;
  abcd[1489] = 2.0E0*I_ERI_Dxy_Dxz_Fx2z_S_c;
  abcd[1490] = 2.0E0*I_ERI_Dxz_Dxz_Fx2z_S_c;
  abcd[1491] = 2.0E0*I_ERI_D2y_Dxz_Fx2z_S_c;
  abcd[1492] = 2.0E0*I_ERI_Dyz_Dxz_Fx2z_S_c;
  abcd[1493] = 2.0E0*I_ERI_D2z_Dxz_Fx2z_S_c;
  abcd[1494] = 2.0E0*I_ERI_D2x_D2y_Fx2z_S_c;
  abcd[1495] = 2.0E0*I_ERI_Dxy_D2y_Fx2z_S_c;
  abcd[1496] = 2.0E0*I_ERI_Dxz_D2y_Fx2z_S_c;
  abcd[1497] = 2.0E0*I_ERI_D2y_D2y_Fx2z_S_c;
  abcd[1498] = 2.0E0*I_ERI_Dyz_D2y_Fx2z_S_c;
  abcd[1499] = 2.0E0*I_ERI_D2z_D2y_Fx2z_S_c;
  abcd[1500] = 2.0E0*I_ERI_D2x_Dyz_Fx2z_S_c;
  abcd[1501] = 2.0E0*I_ERI_Dxy_Dyz_Fx2z_S_c;
  abcd[1502] = 2.0E0*I_ERI_Dxz_Dyz_Fx2z_S_c;
  abcd[1503] = 2.0E0*I_ERI_D2y_Dyz_Fx2z_S_c;
  abcd[1504] = 2.0E0*I_ERI_Dyz_Dyz_Fx2z_S_c;
  abcd[1505] = 2.0E0*I_ERI_D2z_Dyz_Fx2z_S_c;
  abcd[1506] = 2.0E0*I_ERI_D2x_D2z_Fx2z_S_c;
  abcd[1507] = 2.0E0*I_ERI_Dxy_D2z_Fx2z_S_c;
  abcd[1508] = 2.0E0*I_ERI_Dxz_D2z_Fx2z_S_c;
  abcd[1509] = 2.0E0*I_ERI_D2y_D2z_Fx2z_S_c;
  abcd[1510] = 2.0E0*I_ERI_Dyz_D2z_Fx2z_S_c;
  abcd[1511] = 2.0E0*I_ERI_D2z_D2z_Fx2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_c
   * RHS shell quartet name: SQ_ERI_D_D_P_S
   ************************************************************/
  abcd[1512] = 2.0E0*I_ERI_D2x_D2x_F2xy_S_c;
  abcd[1513] = 2.0E0*I_ERI_Dxy_D2x_F2xy_S_c;
  abcd[1514] = 2.0E0*I_ERI_Dxz_D2x_F2xy_S_c;
  abcd[1515] = 2.0E0*I_ERI_D2y_D2x_F2xy_S_c;
  abcd[1516] = 2.0E0*I_ERI_Dyz_D2x_F2xy_S_c;
  abcd[1517] = 2.0E0*I_ERI_D2z_D2x_F2xy_S_c;
  abcd[1518] = 2.0E0*I_ERI_D2x_Dxy_F2xy_S_c;
  abcd[1519] = 2.0E0*I_ERI_Dxy_Dxy_F2xy_S_c;
  abcd[1520] = 2.0E0*I_ERI_Dxz_Dxy_F2xy_S_c;
  abcd[1521] = 2.0E0*I_ERI_D2y_Dxy_F2xy_S_c;
  abcd[1522] = 2.0E0*I_ERI_Dyz_Dxy_F2xy_S_c;
  abcd[1523] = 2.0E0*I_ERI_D2z_Dxy_F2xy_S_c;
  abcd[1524] = 2.0E0*I_ERI_D2x_Dxz_F2xy_S_c;
  abcd[1525] = 2.0E0*I_ERI_Dxy_Dxz_F2xy_S_c;
  abcd[1526] = 2.0E0*I_ERI_Dxz_Dxz_F2xy_S_c;
  abcd[1527] = 2.0E0*I_ERI_D2y_Dxz_F2xy_S_c;
  abcd[1528] = 2.0E0*I_ERI_Dyz_Dxz_F2xy_S_c;
  abcd[1529] = 2.0E0*I_ERI_D2z_Dxz_F2xy_S_c;
  abcd[1530] = 2.0E0*I_ERI_D2x_D2y_F2xy_S_c;
  abcd[1531] = 2.0E0*I_ERI_Dxy_D2y_F2xy_S_c;
  abcd[1532] = 2.0E0*I_ERI_Dxz_D2y_F2xy_S_c;
  abcd[1533] = 2.0E0*I_ERI_D2y_D2y_F2xy_S_c;
  abcd[1534] = 2.0E0*I_ERI_Dyz_D2y_F2xy_S_c;
  abcd[1535] = 2.0E0*I_ERI_D2z_D2y_F2xy_S_c;
  abcd[1536] = 2.0E0*I_ERI_D2x_Dyz_F2xy_S_c;
  abcd[1537] = 2.0E0*I_ERI_Dxy_Dyz_F2xy_S_c;
  abcd[1538] = 2.0E0*I_ERI_Dxz_Dyz_F2xy_S_c;
  abcd[1539] = 2.0E0*I_ERI_D2y_Dyz_F2xy_S_c;
  abcd[1540] = 2.0E0*I_ERI_Dyz_Dyz_F2xy_S_c;
  abcd[1541] = 2.0E0*I_ERI_D2z_Dyz_F2xy_S_c;
  abcd[1542] = 2.0E0*I_ERI_D2x_D2z_F2xy_S_c;
  abcd[1543] = 2.0E0*I_ERI_Dxy_D2z_F2xy_S_c;
  abcd[1544] = 2.0E0*I_ERI_Dxz_D2z_F2xy_S_c;
  abcd[1545] = 2.0E0*I_ERI_D2y_D2z_F2xy_S_c;
  abcd[1546] = 2.0E0*I_ERI_Dyz_D2z_F2xy_S_c;
  abcd[1547] = 2.0E0*I_ERI_D2z_D2z_F2xy_S_c;
  abcd[1548] = 2.0E0*I_ERI_D2x_D2x_Fx2y_S_c-1*I_ERI_D2x_D2x_Px_S;
  abcd[1549] = 2.0E0*I_ERI_Dxy_D2x_Fx2y_S_c-1*I_ERI_Dxy_D2x_Px_S;
  abcd[1550] = 2.0E0*I_ERI_Dxz_D2x_Fx2y_S_c-1*I_ERI_Dxz_D2x_Px_S;
  abcd[1551] = 2.0E0*I_ERI_D2y_D2x_Fx2y_S_c-1*I_ERI_D2y_D2x_Px_S;
  abcd[1552] = 2.0E0*I_ERI_Dyz_D2x_Fx2y_S_c-1*I_ERI_Dyz_D2x_Px_S;
  abcd[1553] = 2.0E0*I_ERI_D2z_D2x_Fx2y_S_c-1*I_ERI_D2z_D2x_Px_S;
  abcd[1554] = 2.0E0*I_ERI_D2x_Dxy_Fx2y_S_c-1*I_ERI_D2x_Dxy_Px_S;
  abcd[1555] = 2.0E0*I_ERI_Dxy_Dxy_Fx2y_S_c-1*I_ERI_Dxy_Dxy_Px_S;
  abcd[1556] = 2.0E0*I_ERI_Dxz_Dxy_Fx2y_S_c-1*I_ERI_Dxz_Dxy_Px_S;
  abcd[1557] = 2.0E0*I_ERI_D2y_Dxy_Fx2y_S_c-1*I_ERI_D2y_Dxy_Px_S;
  abcd[1558] = 2.0E0*I_ERI_Dyz_Dxy_Fx2y_S_c-1*I_ERI_Dyz_Dxy_Px_S;
  abcd[1559] = 2.0E0*I_ERI_D2z_Dxy_Fx2y_S_c-1*I_ERI_D2z_Dxy_Px_S;
  abcd[1560] = 2.0E0*I_ERI_D2x_Dxz_Fx2y_S_c-1*I_ERI_D2x_Dxz_Px_S;
  abcd[1561] = 2.0E0*I_ERI_Dxy_Dxz_Fx2y_S_c-1*I_ERI_Dxy_Dxz_Px_S;
  abcd[1562] = 2.0E0*I_ERI_Dxz_Dxz_Fx2y_S_c-1*I_ERI_Dxz_Dxz_Px_S;
  abcd[1563] = 2.0E0*I_ERI_D2y_Dxz_Fx2y_S_c-1*I_ERI_D2y_Dxz_Px_S;
  abcd[1564] = 2.0E0*I_ERI_Dyz_Dxz_Fx2y_S_c-1*I_ERI_Dyz_Dxz_Px_S;
  abcd[1565] = 2.0E0*I_ERI_D2z_Dxz_Fx2y_S_c-1*I_ERI_D2z_Dxz_Px_S;
  abcd[1566] = 2.0E0*I_ERI_D2x_D2y_Fx2y_S_c-1*I_ERI_D2x_D2y_Px_S;
  abcd[1567] = 2.0E0*I_ERI_Dxy_D2y_Fx2y_S_c-1*I_ERI_Dxy_D2y_Px_S;
  abcd[1568] = 2.0E0*I_ERI_Dxz_D2y_Fx2y_S_c-1*I_ERI_Dxz_D2y_Px_S;
  abcd[1569] = 2.0E0*I_ERI_D2y_D2y_Fx2y_S_c-1*I_ERI_D2y_D2y_Px_S;
  abcd[1570] = 2.0E0*I_ERI_Dyz_D2y_Fx2y_S_c-1*I_ERI_Dyz_D2y_Px_S;
  abcd[1571] = 2.0E0*I_ERI_D2z_D2y_Fx2y_S_c-1*I_ERI_D2z_D2y_Px_S;
  abcd[1572] = 2.0E0*I_ERI_D2x_Dyz_Fx2y_S_c-1*I_ERI_D2x_Dyz_Px_S;
  abcd[1573] = 2.0E0*I_ERI_Dxy_Dyz_Fx2y_S_c-1*I_ERI_Dxy_Dyz_Px_S;
  abcd[1574] = 2.0E0*I_ERI_Dxz_Dyz_Fx2y_S_c-1*I_ERI_Dxz_Dyz_Px_S;
  abcd[1575] = 2.0E0*I_ERI_D2y_Dyz_Fx2y_S_c-1*I_ERI_D2y_Dyz_Px_S;
  abcd[1576] = 2.0E0*I_ERI_Dyz_Dyz_Fx2y_S_c-1*I_ERI_Dyz_Dyz_Px_S;
  abcd[1577] = 2.0E0*I_ERI_D2z_Dyz_Fx2y_S_c-1*I_ERI_D2z_Dyz_Px_S;
  abcd[1578] = 2.0E0*I_ERI_D2x_D2z_Fx2y_S_c-1*I_ERI_D2x_D2z_Px_S;
  abcd[1579] = 2.0E0*I_ERI_Dxy_D2z_Fx2y_S_c-1*I_ERI_Dxy_D2z_Px_S;
  abcd[1580] = 2.0E0*I_ERI_Dxz_D2z_Fx2y_S_c-1*I_ERI_Dxz_D2z_Px_S;
  abcd[1581] = 2.0E0*I_ERI_D2y_D2z_Fx2y_S_c-1*I_ERI_D2y_D2z_Px_S;
  abcd[1582] = 2.0E0*I_ERI_Dyz_D2z_Fx2y_S_c-1*I_ERI_Dyz_D2z_Px_S;
  abcd[1583] = 2.0E0*I_ERI_D2z_D2z_Fx2y_S_c-1*I_ERI_D2z_D2z_Px_S;
  abcd[1584] = 2.0E0*I_ERI_D2x_D2x_Fxyz_S_c;
  abcd[1585] = 2.0E0*I_ERI_Dxy_D2x_Fxyz_S_c;
  abcd[1586] = 2.0E0*I_ERI_Dxz_D2x_Fxyz_S_c;
  abcd[1587] = 2.0E0*I_ERI_D2y_D2x_Fxyz_S_c;
  abcd[1588] = 2.0E0*I_ERI_Dyz_D2x_Fxyz_S_c;
  abcd[1589] = 2.0E0*I_ERI_D2z_D2x_Fxyz_S_c;
  abcd[1590] = 2.0E0*I_ERI_D2x_Dxy_Fxyz_S_c;
  abcd[1591] = 2.0E0*I_ERI_Dxy_Dxy_Fxyz_S_c;
  abcd[1592] = 2.0E0*I_ERI_Dxz_Dxy_Fxyz_S_c;
  abcd[1593] = 2.0E0*I_ERI_D2y_Dxy_Fxyz_S_c;
  abcd[1594] = 2.0E0*I_ERI_Dyz_Dxy_Fxyz_S_c;
  abcd[1595] = 2.0E0*I_ERI_D2z_Dxy_Fxyz_S_c;
  abcd[1596] = 2.0E0*I_ERI_D2x_Dxz_Fxyz_S_c;
  abcd[1597] = 2.0E0*I_ERI_Dxy_Dxz_Fxyz_S_c;
  abcd[1598] = 2.0E0*I_ERI_Dxz_Dxz_Fxyz_S_c;
  abcd[1599] = 2.0E0*I_ERI_D2y_Dxz_Fxyz_S_c;
  abcd[1600] = 2.0E0*I_ERI_Dyz_Dxz_Fxyz_S_c;
  abcd[1601] = 2.0E0*I_ERI_D2z_Dxz_Fxyz_S_c;
  abcd[1602] = 2.0E0*I_ERI_D2x_D2y_Fxyz_S_c;
  abcd[1603] = 2.0E0*I_ERI_Dxy_D2y_Fxyz_S_c;
  abcd[1604] = 2.0E0*I_ERI_Dxz_D2y_Fxyz_S_c;
  abcd[1605] = 2.0E0*I_ERI_D2y_D2y_Fxyz_S_c;
  abcd[1606] = 2.0E0*I_ERI_Dyz_D2y_Fxyz_S_c;
  abcd[1607] = 2.0E0*I_ERI_D2z_D2y_Fxyz_S_c;
  abcd[1608] = 2.0E0*I_ERI_D2x_Dyz_Fxyz_S_c;
  abcd[1609] = 2.0E0*I_ERI_Dxy_Dyz_Fxyz_S_c;
  abcd[1610] = 2.0E0*I_ERI_Dxz_Dyz_Fxyz_S_c;
  abcd[1611] = 2.0E0*I_ERI_D2y_Dyz_Fxyz_S_c;
  abcd[1612] = 2.0E0*I_ERI_Dyz_Dyz_Fxyz_S_c;
  abcd[1613] = 2.0E0*I_ERI_D2z_Dyz_Fxyz_S_c;
  abcd[1614] = 2.0E0*I_ERI_D2x_D2z_Fxyz_S_c;
  abcd[1615] = 2.0E0*I_ERI_Dxy_D2z_Fxyz_S_c;
  abcd[1616] = 2.0E0*I_ERI_Dxz_D2z_Fxyz_S_c;
  abcd[1617] = 2.0E0*I_ERI_D2y_D2z_Fxyz_S_c;
  abcd[1618] = 2.0E0*I_ERI_Dyz_D2z_Fxyz_S_c;
  abcd[1619] = 2.0E0*I_ERI_D2z_D2z_Fxyz_S_c;
  abcd[1620] = 2.0E0*I_ERI_D2x_D2x_F3y_S_c-2*I_ERI_D2x_D2x_Py_S;
  abcd[1621] = 2.0E0*I_ERI_Dxy_D2x_F3y_S_c-2*I_ERI_Dxy_D2x_Py_S;
  abcd[1622] = 2.0E0*I_ERI_Dxz_D2x_F3y_S_c-2*I_ERI_Dxz_D2x_Py_S;
  abcd[1623] = 2.0E0*I_ERI_D2y_D2x_F3y_S_c-2*I_ERI_D2y_D2x_Py_S;
  abcd[1624] = 2.0E0*I_ERI_Dyz_D2x_F3y_S_c-2*I_ERI_Dyz_D2x_Py_S;
  abcd[1625] = 2.0E0*I_ERI_D2z_D2x_F3y_S_c-2*I_ERI_D2z_D2x_Py_S;
  abcd[1626] = 2.0E0*I_ERI_D2x_Dxy_F3y_S_c-2*I_ERI_D2x_Dxy_Py_S;
  abcd[1627] = 2.0E0*I_ERI_Dxy_Dxy_F3y_S_c-2*I_ERI_Dxy_Dxy_Py_S;
  abcd[1628] = 2.0E0*I_ERI_Dxz_Dxy_F3y_S_c-2*I_ERI_Dxz_Dxy_Py_S;
  abcd[1629] = 2.0E0*I_ERI_D2y_Dxy_F3y_S_c-2*I_ERI_D2y_Dxy_Py_S;
  abcd[1630] = 2.0E0*I_ERI_Dyz_Dxy_F3y_S_c-2*I_ERI_Dyz_Dxy_Py_S;
  abcd[1631] = 2.0E0*I_ERI_D2z_Dxy_F3y_S_c-2*I_ERI_D2z_Dxy_Py_S;
  abcd[1632] = 2.0E0*I_ERI_D2x_Dxz_F3y_S_c-2*I_ERI_D2x_Dxz_Py_S;
  abcd[1633] = 2.0E0*I_ERI_Dxy_Dxz_F3y_S_c-2*I_ERI_Dxy_Dxz_Py_S;
  abcd[1634] = 2.0E0*I_ERI_Dxz_Dxz_F3y_S_c-2*I_ERI_Dxz_Dxz_Py_S;
  abcd[1635] = 2.0E0*I_ERI_D2y_Dxz_F3y_S_c-2*I_ERI_D2y_Dxz_Py_S;
  abcd[1636] = 2.0E0*I_ERI_Dyz_Dxz_F3y_S_c-2*I_ERI_Dyz_Dxz_Py_S;
  abcd[1637] = 2.0E0*I_ERI_D2z_Dxz_F3y_S_c-2*I_ERI_D2z_Dxz_Py_S;
  abcd[1638] = 2.0E0*I_ERI_D2x_D2y_F3y_S_c-2*I_ERI_D2x_D2y_Py_S;
  abcd[1639] = 2.0E0*I_ERI_Dxy_D2y_F3y_S_c-2*I_ERI_Dxy_D2y_Py_S;
  abcd[1640] = 2.0E0*I_ERI_Dxz_D2y_F3y_S_c-2*I_ERI_Dxz_D2y_Py_S;
  abcd[1641] = 2.0E0*I_ERI_D2y_D2y_F3y_S_c-2*I_ERI_D2y_D2y_Py_S;
  abcd[1642] = 2.0E0*I_ERI_Dyz_D2y_F3y_S_c-2*I_ERI_Dyz_D2y_Py_S;
  abcd[1643] = 2.0E0*I_ERI_D2z_D2y_F3y_S_c-2*I_ERI_D2z_D2y_Py_S;
  abcd[1644] = 2.0E0*I_ERI_D2x_Dyz_F3y_S_c-2*I_ERI_D2x_Dyz_Py_S;
  abcd[1645] = 2.0E0*I_ERI_Dxy_Dyz_F3y_S_c-2*I_ERI_Dxy_Dyz_Py_S;
  abcd[1646] = 2.0E0*I_ERI_Dxz_Dyz_F3y_S_c-2*I_ERI_Dxz_Dyz_Py_S;
  abcd[1647] = 2.0E0*I_ERI_D2y_Dyz_F3y_S_c-2*I_ERI_D2y_Dyz_Py_S;
  abcd[1648] = 2.0E0*I_ERI_Dyz_Dyz_F3y_S_c-2*I_ERI_Dyz_Dyz_Py_S;
  abcd[1649] = 2.0E0*I_ERI_D2z_Dyz_F3y_S_c-2*I_ERI_D2z_Dyz_Py_S;
  abcd[1650] = 2.0E0*I_ERI_D2x_D2z_F3y_S_c-2*I_ERI_D2x_D2z_Py_S;
  abcd[1651] = 2.0E0*I_ERI_Dxy_D2z_F3y_S_c-2*I_ERI_Dxy_D2z_Py_S;
  abcd[1652] = 2.0E0*I_ERI_Dxz_D2z_F3y_S_c-2*I_ERI_Dxz_D2z_Py_S;
  abcd[1653] = 2.0E0*I_ERI_D2y_D2z_F3y_S_c-2*I_ERI_D2y_D2z_Py_S;
  abcd[1654] = 2.0E0*I_ERI_Dyz_D2z_F3y_S_c-2*I_ERI_Dyz_D2z_Py_S;
  abcd[1655] = 2.0E0*I_ERI_D2z_D2z_F3y_S_c-2*I_ERI_D2z_D2z_Py_S;
  abcd[1656] = 2.0E0*I_ERI_D2x_D2x_F2yz_S_c-1*I_ERI_D2x_D2x_Pz_S;
  abcd[1657] = 2.0E0*I_ERI_Dxy_D2x_F2yz_S_c-1*I_ERI_Dxy_D2x_Pz_S;
  abcd[1658] = 2.0E0*I_ERI_Dxz_D2x_F2yz_S_c-1*I_ERI_Dxz_D2x_Pz_S;
  abcd[1659] = 2.0E0*I_ERI_D2y_D2x_F2yz_S_c-1*I_ERI_D2y_D2x_Pz_S;
  abcd[1660] = 2.0E0*I_ERI_Dyz_D2x_F2yz_S_c-1*I_ERI_Dyz_D2x_Pz_S;
  abcd[1661] = 2.0E0*I_ERI_D2z_D2x_F2yz_S_c-1*I_ERI_D2z_D2x_Pz_S;
  abcd[1662] = 2.0E0*I_ERI_D2x_Dxy_F2yz_S_c-1*I_ERI_D2x_Dxy_Pz_S;
  abcd[1663] = 2.0E0*I_ERI_Dxy_Dxy_F2yz_S_c-1*I_ERI_Dxy_Dxy_Pz_S;
  abcd[1664] = 2.0E0*I_ERI_Dxz_Dxy_F2yz_S_c-1*I_ERI_Dxz_Dxy_Pz_S;
  abcd[1665] = 2.0E0*I_ERI_D2y_Dxy_F2yz_S_c-1*I_ERI_D2y_Dxy_Pz_S;
  abcd[1666] = 2.0E0*I_ERI_Dyz_Dxy_F2yz_S_c-1*I_ERI_Dyz_Dxy_Pz_S;
  abcd[1667] = 2.0E0*I_ERI_D2z_Dxy_F2yz_S_c-1*I_ERI_D2z_Dxy_Pz_S;
  abcd[1668] = 2.0E0*I_ERI_D2x_Dxz_F2yz_S_c-1*I_ERI_D2x_Dxz_Pz_S;
  abcd[1669] = 2.0E0*I_ERI_Dxy_Dxz_F2yz_S_c-1*I_ERI_Dxy_Dxz_Pz_S;
  abcd[1670] = 2.0E0*I_ERI_Dxz_Dxz_F2yz_S_c-1*I_ERI_Dxz_Dxz_Pz_S;
  abcd[1671] = 2.0E0*I_ERI_D2y_Dxz_F2yz_S_c-1*I_ERI_D2y_Dxz_Pz_S;
  abcd[1672] = 2.0E0*I_ERI_Dyz_Dxz_F2yz_S_c-1*I_ERI_Dyz_Dxz_Pz_S;
  abcd[1673] = 2.0E0*I_ERI_D2z_Dxz_F2yz_S_c-1*I_ERI_D2z_Dxz_Pz_S;
  abcd[1674] = 2.0E0*I_ERI_D2x_D2y_F2yz_S_c-1*I_ERI_D2x_D2y_Pz_S;
  abcd[1675] = 2.0E0*I_ERI_Dxy_D2y_F2yz_S_c-1*I_ERI_Dxy_D2y_Pz_S;
  abcd[1676] = 2.0E0*I_ERI_Dxz_D2y_F2yz_S_c-1*I_ERI_Dxz_D2y_Pz_S;
  abcd[1677] = 2.0E0*I_ERI_D2y_D2y_F2yz_S_c-1*I_ERI_D2y_D2y_Pz_S;
  abcd[1678] = 2.0E0*I_ERI_Dyz_D2y_F2yz_S_c-1*I_ERI_Dyz_D2y_Pz_S;
  abcd[1679] = 2.0E0*I_ERI_D2z_D2y_F2yz_S_c-1*I_ERI_D2z_D2y_Pz_S;
  abcd[1680] = 2.0E0*I_ERI_D2x_Dyz_F2yz_S_c-1*I_ERI_D2x_Dyz_Pz_S;
  abcd[1681] = 2.0E0*I_ERI_Dxy_Dyz_F2yz_S_c-1*I_ERI_Dxy_Dyz_Pz_S;
  abcd[1682] = 2.0E0*I_ERI_Dxz_Dyz_F2yz_S_c-1*I_ERI_Dxz_Dyz_Pz_S;
  abcd[1683] = 2.0E0*I_ERI_D2y_Dyz_F2yz_S_c-1*I_ERI_D2y_Dyz_Pz_S;
  abcd[1684] = 2.0E0*I_ERI_Dyz_Dyz_F2yz_S_c-1*I_ERI_Dyz_Dyz_Pz_S;
  abcd[1685] = 2.0E0*I_ERI_D2z_Dyz_F2yz_S_c-1*I_ERI_D2z_Dyz_Pz_S;
  abcd[1686] = 2.0E0*I_ERI_D2x_D2z_F2yz_S_c-1*I_ERI_D2x_D2z_Pz_S;
  abcd[1687] = 2.0E0*I_ERI_Dxy_D2z_F2yz_S_c-1*I_ERI_Dxy_D2z_Pz_S;
  abcd[1688] = 2.0E0*I_ERI_Dxz_D2z_F2yz_S_c-1*I_ERI_Dxz_D2z_Pz_S;
  abcd[1689] = 2.0E0*I_ERI_D2y_D2z_F2yz_S_c-1*I_ERI_D2y_D2z_Pz_S;
  abcd[1690] = 2.0E0*I_ERI_Dyz_D2z_F2yz_S_c-1*I_ERI_Dyz_D2z_Pz_S;
  abcd[1691] = 2.0E0*I_ERI_D2z_D2z_F2yz_S_c-1*I_ERI_D2z_D2z_Pz_S;
  abcd[1692] = 2.0E0*I_ERI_D2x_D2x_Fy2z_S_c;
  abcd[1693] = 2.0E0*I_ERI_Dxy_D2x_Fy2z_S_c;
  abcd[1694] = 2.0E0*I_ERI_Dxz_D2x_Fy2z_S_c;
  abcd[1695] = 2.0E0*I_ERI_D2y_D2x_Fy2z_S_c;
  abcd[1696] = 2.0E0*I_ERI_Dyz_D2x_Fy2z_S_c;
  abcd[1697] = 2.0E0*I_ERI_D2z_D2x_Fy2z_S_c;
  abcd[1698] = 2.0E0*I_ERI_D2x_Dxy_Fy2z_S_c;
  abcd[1699] = 2.0E0*I_ERI_Dxy_Dxy_Fy2z_S_c;
  abcd[1700] = 2.0E0*I_ERI_Dxz_Dxy_Fy2z_S_c;
  abcd[1701] = 2.0E0*I_ERI_D2y_Dxy_Fy2z_S_c;
  abcd[1702] = 2.0E0*I_ERI_Dyz_Dxy_Fy2z_S_c;
  abcd[1703] = 2.0E0*I_ERI_D2z_Dxy_Fy2z_S_c;
  abcd[1704] = 2.0E0*I_ERI_D2x_Dxz_Fy2z_S_c;
  abcd[1705] = 2.0E0*I_ERI_Dxy_Dxz_Fy2z_S_c;
  abcd[1706] = 2.0E0*I_ERI_Dxz_Dxz_Fy2z_S_c;
  abcd[1707] = 2.0E0*I_ERI_D2y_Dxz_Fy2z_S_c;
  abcd[1708] = 2.0E0*I_ERI_Dyz_Dxz_Fy2z_S_c;
  abcd[1709] = 2.0E0*I_ERI_D2z_Dxz_Fy2z_S_c;
  abcd[1710] = 2.0E0*I_ERI_D2x_D2y_Fy2z_S_c;
  abcd[1711] = 2.0E0*I_ERI_Dxy_D2y_Fy2z_S_c;
  abcd[1712] = 2.0E0*I_ERI_Dxz_D2y_Fy2z_S_c;
  abcd[1713] = 2.0E0*I_ERI_D2y_D2y_Fy2z_S_c;
  abcd[1714] = 2.0E0*I_ERI_Dyz_D2y_Fy2z_S_c;
  abcd[1715] = 2.0E0*I_ERI_D2z_D2y_Fy2z_S_c;
  abcd[1716] = 2.0E0*I_ERI_D2x_Dyz_Fy2z_S_c;
  abcd[1717] = 2.0E0*I_ERI_Dxy_Dyz_Fy2z_S_c;
  abcd[1718] = 2.0E0*I_ERI_Dxz_Dyz_Fy2z_S_c;
  abcd[1719] = 2.0E0*I_ERI_D2y_Dyz_Fy2z_S_c;
  abcd[1720] = 2.0E0*I_ERI_Dyz_Dyz_Fy2z_S_c;
  abcd[1721] = 2.0E0*I_ERI_D2z_Dyz_Fy2z_S_c;
  abcd[1722] = 2.0E0*I_ERI_D2x_D2z_Fy2z_S_c;
  abcd[1723] = 2.0E0*I_ERI_Dxy_D2z_Fy2z_S_c;
  abcd[1724] = 2.0E0*I_ERI_Dxz_D2z_Fy2z_S_c;
  abcd[1725] = 2.0E0*I_ERI_D2y_D2z_Fy2z_S_c;
  abcd[1726] = 2.0E0*I_ERI_Dyz_D2z_Fy2z_S_c;
  abcd[1727] = 2.0E0*I_ERI_D2z_D2z_Fy2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_F_S_c
   * RHS shell quartet name: SQ_ERI_D_D_P_S
   ************************************************************/
  abcd[1728] = 2.0E0*I_ERI_D2x_D2x_F2xz_S_c;
  abcd[1729] = 2.0E0*I_ERI_Dxy_D2x_F2xz_S_c;
  abcd[1730] = 2.0E0*I_ERI_Dxz_D2x_F2xz_S_c;
  abcd[1731] = 2.0E0*I_ERI_D2y_D2x_F2xz_S_c;
  abcd[1732] = 2.0E0*I_ERI_Dyz_D2x_F2xz_S_c;
  abcd[1733] = 2.0E0*I_ERI_D2z_D2x_F2xz_S_c;
  abcd[1734] = 2.0E0*I_ERI_D2x_Dxy_F2xz_S_c;
  abcd[1735] = 2.0E0*I_ERI_Dxy_Dxy_F2xz_S_c;
  abcd[1736] = 2.0E0*I_ERI_Dxz_Dxy_F2xz_S_c;
  abcd[1737] = 2.0E0*I_ERI_D2y_Dxy_F2xz_S_c;
  abcd[1738] = 2.0E0*I_ERI_Dyz_Dxy_F2xz_S_c;
  abcd[1739] = 2.0E0*I_ERI_D2z_Dxy_F2xz_S_c;
  abcd[1740] = 2.0E0*I_ERI_D2x_Dxz_F2xz_S_c;
  abcd[1741] = 2.0E0*I_ERI_Dxy_Dxz_F2xz_S_c;
  abcd[1742] = 2.0E0*I_ERI_Dxz_Dxz_F2xz_S_c;
  abcd[1743] = 2.0E0*I_ERI_D2y_Dxz_F2xz_S_c;
  abcd[1744] = 2.0E0*I_ERI_Dyz_Dxz_F2xz_S_c;
  abcd[1745] = 2.0E0*I_ERI_D2z_Dxz_F2xz_S_c;
  abcd[1746] = 2.0E0*I_ERI_D2x_D2y_F2xz_S_c;
  abcd[1747] = 2.0E0*I_ERI_Dxy_D2y_F2xz_S_c;
  abcd[1748] = 2.0E0*I_ERI_Dxz_D2y_F2xz_S_c;
  abcd[1749] = 2.0E0*I_ERI_D2y_D2y_F2xz_S_c;
  abcd[1750] = 2.0E0*I_ERI_Dyz_D2y_F2xz_S_c;
  abcd[1751] = 2.0E0*I_ERI_D2z_D2y_F2xz_S_c;
  abcd[1752] = 2.0E0*I_ERI_D2x_Dyz_F2xz_S_c;
  abcd[1753] = 2.0E0*I_ERI_Dxy_Dyz_F2xz_S_c;
  abcd[1754] = 2.0E0*I_ERI_Dxz_Dyz_F2xz_S_c;
  abcd[1755] = 2.0E0*I_ERI_D2y_Dyz_F2xz_S_c;
  abcd[1756] = 2.0E0*I_ERI_Dyz_Dyz_F2xz_S_c;
  abcd[1757] = 2.0E0*I_ERI_D2z_Dyz_F2xz_S_c;
  abcd[1758] = 2.0E0*I_ERI_D2x_D2z_F2xz_S_c;
  abcd[1759] = 2.0E0*I_ERI_Dxy_D2z_F2xz_S_c;
  abcd[1760] = 2.0E0*I_ERI_Dxz_D2z_F2xz_S_c;
  abcd[1761] = 2.0E0*I_ERI_D2y_D2z_F2xz_S_c;
  abcd[1762] = 2.0E0*I_ERI_Dyz_D2z_F2xz_S_c;
  abcd[1763] = 2.0E0*I_ERI_D2z_D2z_F2xz_S_c;
  abcd[1764] = 2.0E0*I_ERI_D2x_D2x_Fxyz_S_c;
  abcd[1765] = 2.0E0*I_ERI_Dxy_D2x_Fxyz_S_c;
  abcd[1766] = 2.0E0*I_ERI_Dxz_D2x_Fxyz_S_c;
  abcd[1767] = 2.0E0*I_ERI_D2y_D2x_Fxyz_S_c;
  abcd[1768] = 2.0E0*I_ERI_Dyz_D2x_Fxyz_S_c;
  abcd[1769] = 2.0E0*I_ERI_D2z_D2x_Fxyz_S_c;
  abcd[1770] = 2.0E0*I_ERI_D2x_Dxy_Fxyz_S_c;
  abcd[1771] = 2.0E0*I_ERI_Dxy_Dxy_Fxyz_S_c;
  abcd[1772] = 2.0E0*I_ERI_Dxz_Dxy_Fxyz_S_c;
  abcd[1773] = 2.0E0*I_ERI_D2y_Dxy_Fxyz_S_c;
  abcd[1774] = 2.0E0*I_ERI_Dyz_Dxy_Fxyz_S_c;
  abcd[1775] = 2.0E0*I_ERI_D2z_Dxy_Fxyz_S_c;
  abcd[1776] = 2.0E0*I_ERI_D2x_Dxz_Fxyz_S_c;
  abcd[1777] = 2.0E0*I_ERI_Dxy_Dxz_Fxyz_S_c;
  abcd[1778] = 2.0E0*I_ERI_Dxz_Dxz_Fxyz_S_c;
  abcd[1779] = 2.0E0*I_ERI_D2y_Dxz_Fxyz_S_c;
  abcd[1780] = 2.0E0*I_ERI_Dyz_Dxz_Fxyz_S_c;
  abcd[1781] = 2.0E0*I_ERI_D2z_Dxz_Fxyz_S_c;
  abcd[1782] = 2.0E0*I_ERI_D2x_D2y_Fxyz_S_c;
  abcd[1783] = 2.0E0*I_ERI_Dxy_D2y_Fxyz_S_c;
  abcd[1784] = 2.0E0*I_ERI_Dxz_D2y_Fxyz_S_c;
  abcd[1785] = 2.0E0*I_ERI_D2y_D2y_Fxyz_S_c;
  abcd[1786] = 2.0E0*I_ERI_Dyz_D2y_Fxyz_S_c;
  abcd[1787] = 2.0E0*I_ERI_D2z_D2y_Fxyz_S_c;
  abcd[1788] = 2.0E0*I_ERI_D2x_Dyz_Fxyz_S_c;
  abcd[1789] = 2.0E0*I_ERI_Dxy_Dyz_Fxyz_S_c;
  abcd[1790] = 2.0E0*I_ERI_Dxz_Dyz_Fxyz_S_c;
  abcd[1791] = 2.0E0*I_ERI_D2y_Dyz_Fxyz_S_c;
  abcd[1792] = 2.0E0*I_ERI_Dyz_Dyz_Fxyz_S_c;
  abcd[1793] = 2.0E0*I_ERI_D2z_Dyz_Fxyz_S_c;
  abcd[1794] = 2.0E0*I_ERI_D2x_D2z_Fxyz_S_c;
  abcd[1795] = 2.0E0*I_ERI_Dxy_D2z_Fxyz_S_c;
  abcd[1796] = 2.0E0*I_ERI_Dxz_D2z_Fxyz_S_c;
  abcd[1797] = 2.0E0*I_ERI_D2y_D2z_Fxyz_S_c;
  abcd[1798] = 2.0E0*I_ERI_Dyz_D2z_Fxyz_S_c;
  abcd[1799] = 2.0E0*I_ERI_D2z_D2z_Fxyz_S_c;
  abcd[1800] = 2.0E0*I_ERI_D2x_D2x_Fx2z_S_c-1*I_ERI_D2x_D2x_Px_S;
  abcd[1801] = 2.0E0*I_ERI_Dxy_D2x_Fx2z_S_c-1*I_ERI_Dxy_D2x_Px_S;
  abcd[1802] = 2.0E0*I_ERI_Dxz_D2x_Fx2z_S_c-1*I_ERI_Dxz_D2x_Px_S;
  abcd[1803] = 2.0E0*I_ERI_D2y_D2x_Fx2z_S_c-1*I_ERI_D2y_D2x_Px_S;
  abcd[1804] = 2.0E0*I_ERI_Dyz_D2x_Fx2z_S_c-1*I_ERI_Dyz_D2x_Px_S;
  abcd[1805] = 2.0E0*I_ERI_D2z_D2x_Fx2z_S_c-1*I_ERI_D2z_D2x_Px_S;
  abcd[1806] = 2.0E0*I_ERI_D2x_Dxy_Fx2z_S_c-1*I_ERI_D2x_Dxy_Px_S;
  abcd[1807] = 2.0E0*I_ERI_Dxy_Dxy_Fx2z_S_c-1*I_ERI_Dxy_Dxy_Px_S;
  abcd[1808] = 2.0E0*I_ERI_Dxz_Dxy_Fx2z_S_c-1*I_ERI_Dxz_Dxy_Px_S;
  abcd[1809] = 2.0E0*I_ERI_D2y_Dxy_Fx2z_S_c-1*I_ERI_D2y_Dxy_Px_S;
  abcd[1810] = 2.0E0*I_ERI_Dyz_Dxy_Fx2z_S_c-1*I_ERI_Dyz_Dxy_Px_S;
  abcd[1811] = 2.0E0*I_ERI_D2z_Dxy_Fx2z_S_c-1*I_ERI_D2z_Dxy_Px_S;
  abcd[1812] = 2.0E0*I_ERI_D2x_Dxz_Fx2z_S_c-1*I_ERI_D2x_Dxz_Px_S;
  abcd[1813] = 2.0E0*I_ERI_Dxy_Dxz_Fx2z_S_c-1*I_ERI_Dxy_Dxz_Px_S;
  abcd[1814] = 2.0E0*I_ERI_Dxz_Dxz_Fx2z_S_c-1*I_ERI_Dxz_Dxz_Px_S;
  abcd[1815] = 2.0E0*I_ERI_D2y_Dxz_Fx2z_S_c-1*I_ERI_D2y_Dxz_Px_S;
  abcd[1816] = 2.0E0*I_ERI_Dyz_Dxz_Fx2z_S_c-1*I_ERI_Dyz_Dxz_Px_S;
  abcd[1817] = 2.0E0*I_ERI_D2z_Dxz_Fx2z_S_c-1*I_ERI_D2z_Dxz_Px_S;
  abcd[1818] = 2.0E0*I_ERI_D2x_D2y_Fx2z_S_c-1*I_ERI_D2x_D2y_Px_S;
  abcd[1819] = 2.0E0*I_ERI_Dxy_D2y_Fx2z_S_c-1*I_ERI_Dxy_D2y_Px_S;
  abcd[1820] = 2.0E0*I_ERI_Dxz_D2y_Fx2z_S_c-1*I_ERI_Dxz_D2y_Px_S;
  abcd[1821] = 2.0E0*I_ERI_D2y_D2y_Fx2z_S_c-1*I_ERI_D2y_D2y_Px_S;
  abcd[1822] = 2.0E0*I_ERI_Dyz_D2y_Fx2z_S_c-1*I_ERI_Dyz_D2y_Px_S;
  abcd[1823] = 2.0E0*I_ERI_D2z_D2y_Fx2z_S_c-1*I_ERI_D2z_D2y_Px_S;
  abcd[1824] = 2.0E0*I_ERI_D2x_Dyz_Fx2z_S_c-1*I_ERI_D2x_Dyz_Px_S;
  abcd[1825] = 2.0E0*I_ERI_Dxy_Dyz_Fx2z_S_c-1*I_ERI_Dxy_Dyz_Px_S;
  abcd[1826] = 2.0E0*I_ERI_Dxz_Dyz_Fx2z_S_c-1*I_ERI_Dxz_Dyz_Px_S;
  abcd[1827] = 2.0E0*I_ERI_D2y_Dyz_Fx2z_S_c-1*I_ERI_D2y_Dyz_Px_S;
  abcd[1828] = 2.0E0*I_ERI_Dyz_Dyz_Fx2z_S_c-1*I_ERI_Dyz_Dyz_Px_S;
  abcd[1829] = 2.0E0*I_ERI_D2z_Dyz_Fx2z_S_c-1*I_ERI_D2z_Dyz_Px_S;
  abcd[1830] = 2.0E0*I_ERI_D2x_D2z_Fx2z_S_c-1*I_ERI_D2x_D2z_Px_S;
  abcd[1831] = 2.0E0*I_ERI_Dxy_D2z_Fx2z_S_c-1*I_ERI_Dxy_D2z_Px_S;
  abcd[1832] = 2.0E0*I_ERI_Dxz_D2z_Fx2z_S_c-1*I_ERI_Dxz_D2z_Px_S;
  abcd[1833] = 2.0E0*I_ERI_D2y_D2z_Fx2z_S_c-1*I_ERI_D2y_D2z_Px_S;
  abcd[1834] = 2.0E0*I_ERI_Dyz_D2z_Fx2z_S_c-1*I_ERI_Dyz_D2z_Px_S;
  abcd[1835] = 2.0E0*I_ERI_D2z_D2z_Fx2z_S_c-1*I_ERI_D2z_D2z_Px_S;
  abcd[1836] = 2.0E0*I_ERI_D2x_D2x_F2yz_S_c;
  abcd[1837] = 2.0E0*I_ERI_Dxy_D2x_F2yz_S_c;
  abcd[1838] = 2.0E0*I_ERI_Dxz_D2x_F2yz_S_c;
  abcd[1839] = 2.0E0*I_ERI_D2y_D2x_F2yz_S_c;
  abcd[1840] = 2.0E0*I_ERI_Dyz_D2x_F2yz_S_c;
  abcd[1841] = 2.0E0*I_ERI_D2z_D2x_F2yz_S_c;
  abcd[1842] = 2.0E0*I_ERI_D2x_Dxy_F2yz_S_c;
  abcd[1843] = 2.0E0*I_ERI_Dxy_Dxy_F2yz_S_c;
  abcd[1844] = 2.0E0*I_ERI_Dxz_Dxy_F2yz_S_c;
  abcd[1845] = 2.0E0*I_ERI_D2y_Dxy_F2yz_S_c;
  abcd[1846] = 2.0E0*I_ERI_Dyz_Dxy_F2yz_S_c;
  abcd[1847] = 2.0E0*I_ERI_D2z_Dxy_F2yz_S_c;
  abcd[1848] = 2.0E0*I_ERI_D2x_Dxz_F2yz_S_c;
  abcd[1849] = 2.0E0*I_ERI_Dxy_Dxz_F2yz_S_c;
  abcd[1850] = 2.0E0*I_ERI_Dxz_Dxz_F2yz_S_c;
  abcd[1851] = 2.0E0*I_ERI_D2y_Dxz_F2yz_S_c;
  abcd[1852] = 2.0E0*I_ERI_Dyz_Dxz_F2yz_S_c;
  abcd[1853] = 2.0E0*I_ERI_D2z_Dxz_F2yz_S_c;
  abcd[1854] = 2.0E0*I_ERI_D2x_D2y_F2yz_S_c;
  abcd[1855] = 2.0E0*I_ERI_Dxy_D2y_F2yz_S_c;
  abcd[1856] = 2.0E0*I_ERI_Dxz_D2y_F2yz_S_c;
  abcd[1857] = 2.0E0*I_ERI_D2y_D2y_F2yz_S_c;
  abcd[1858] = 2.0E0*I_ERI_Dyz_D2y_F2yz_S_c;
  abcd[1859] = 2.0E0*I_ERI_D2z_D2y_F2yz_S_c;
  abcd[1860] = 2.0E0*I_ERI_D2x_Dyz_F2yz_S_c;
  abcd[1861] = 2.0E0*I_ERI_Dxy_Dyz_F2yz_S_c;
  abcd[1862] = 2.0E0*I_ERI_Dxz_Dyz_F2yz_S_c;
  abcd[1863] = 2.0E0*I_ERI_D2y_Dyz_F2yz_S_c;
  abcd[1864] = 2.0E0*I_ERI_Dyz_Dyz_F2yz_S_c;
  abcd[1865] = 2.0E0*I_ERI_D2z_Dyz_F2yz_S_c;
  abcd[1866] = 2.0E0*I_ERI_D2x_D2z_F2yz_S_c;
  abcd[1867] = 2.0E0*I_ERI_Dxy_D2z_F2yz_S_c;
  abcd[1868] = 2.0E0*I_ERI_Dxz_D2z_F2yz_S_c;
  abcd[1869] = 2.0E0*I_ERI_D2y_D2z_F2yz_S_c;
  abcd[1870] = 2.0E0*I_ERI_Dyz_D2z_F2yz_S_c;
  abcd[1871] = 2.0E0*I_ERI_D2z_D2z_F2yz_S_c;
  abcd[1872] = 2.0E0*I_ERI_D2x_D2x_Fy2z_S_c-1*I_ERI_D2x_D2x_Py_S;
  abcd[1873] = 2.0E0*I_ERI_Dxy_D2x_Fy2z_S_c-1*I_ERI_Dxy_D2x_Py_S;
  abcd[1874] = 2.0E0*I_ERI_Dxz_D2x_Fy2z_S_c-1*I_ERI_Dxz_D2x_Py_S;
  abcd[1875] = 2.0E0*I_ERI_D2y_D2x_Fy2z_S_c-1*I_ERI_D2y_D2x_Py_S;
  abcd[1876] = 2.0E0*I_ERI_Dyz_D2x_Fy2z_S_c-1*I_ERI_Dyz_D2x_Py_S;
  abcd[1877] = 2.0E0*I_ERI_D2z_D2x_Fy2z_S_c-1*I_ERI_D2z_D2x_Py_S;
  abcd[1878] = 2.0E0*I_ERI_D2x_Dxy_Fy2z_S_c-1*I_ERI_D2x_Dxy_Py_S;
  abcd[1879] = 2.0E0*I_ERI_Dxy_Dxy_Fy2z_S_c-1*I_ERI_Dxy_Dxy_Py_S;
  abcd[1880] = 2.0E0*I_ERI_Dxz_Dxy_Fy2z_S_c-1*I_ERI_Dxz_Dxy_Py_S;
  abcd[1881] = 2.0E0*I_ERI_D2y_Dxy_Fy2z_S_c-1*I_ERI_D2y_Dxy_Py_S;
  abcd[1882] = 2.0E0*I_ERI_Dyz_Dxy_Fy2z_S_c-1*I_ERI_Dyz_Dxy_Py_S;
  abcd[1883] = 2.0E0*I_ERI_D2z_Dxy_Fy2z_S_c-1*I_ERI_D2z_Dxy_Py_S;
  abcd[1884] = 2.0E0*I_ERI_D2x_Dxz_Fy2z_S_c-1*I_ERI_D2x_Dxz_Py_S;
  abcd[1885] = 2.0E0*I_ERI_Dxy_Dxz_Fy2z_S_c-1*I_ERI_Dxy_Dxz_Py_S;
  abcd[1886] = 2.0E0*I_ERI_Dxz_Dxz_Fy2z_S_c-1*I_ERI_Dxz_Dxz_Py_S;
  abcd[1887] = 2.0E0*I_ERI_D2y_Dxz_Fy2z_S_c-1*I_ERI_D2y_Dxz_Py_S;
  abcd[1888] = 2.0E0*I_ERI_Dyz_Dxz_Fy2z_S_c-1*I_ERI_Dyz_Dxz_Py_S;
  abcd[1889] = 2.0E0*I_ERI_D2z_Dxz_Fy2z_S_c-1*I_ERI_D2z_Dxz_Py_S;
  abcd[1890] = 2.0E0*I_ERI_D2x_D2y_Fy2z_S_c-1*I_ERI_D2x_D2y_Py_S;
  abcd[1891] = 2.0E0*I_ERI_Dxy_D2y_Fy2z_S_c-1*I_ERI_Dxy_D2y_Py_S;
  abcd[1892] = 2.0E0*I_ERI_Dxz_D2y_Fy2z_S_c-1*I_ERI_Dxz_D2y_Py_S;
  abcd[1893] = 2.0E0*I_ERI_D2y_D2y_Fy2z_S_c-1*I_ERI_D2y_D2y_Py_S;
  abcd[1894] = 2.0E0*I_ERI_Dyz_D2y_Fy2z_S_c-1*I_ERI_Dyz_D2y_Py_S;
  abcd[1895] = 2.0E0*I_ERI_D2z_D2y_Fy2z_S_c-1*I_ERI_D2z_D2y_Py_S;
  abcd[1896] = 2.0E0*I_ERI_D2x_Dyz_Fy2z_S_c-1*I_ERI_D2x_Dyz_Py_S;
  abcd[1897] = 2.0E0*I_ERI_Dxy_Dyz_Fy2z_S_c-1*I_ERI_Dxy_Dyz_Py_S;
  abcd[1898] = 2.0E0*I_ERI_Dxz_Dyz_Fy2z_S_c-1*I_ERI_Dxz_Dyz_Py_S;
  abcd[1899] = 2.0E0*I_ERI_D2y_Dyz_Fy2z_S_c-1*I_ERI_D2y_Dyz_Py_S;
  abcd[1900] = 2.0E0*I_ERI_Dyz_Dyz_Fy2z_S_c-1*I_ERI_Dyz_Dyz_Py_S;
  abcd[1901] = 2.0E0*I_ERI_D2z_Dyz_Fy2z_S_c-1*I_ERI_D2z_Dyz_Py_S;
  abcd[1902] = 2.0E0*I_ERI_D2x_D2z_Fy2z_S_c-1*I_ERI_D2x_D2z_Py_S;
  abcd[1903] = 2.0E0*I_ERI_Dxy_D2z_Fy2z_S_c-1*I_ERI_Dxy_D2z_Py_S;
  abcd[1904] = 2.0E0*I_ERI_Dxz_D2z_Fy2z_S_c-1*I_ERI_Dxz_D2z_Py_S;
  abcd[1905] = 2.0E0*I_ERI_D2y_D2z_Fy2z_S_c-1*I_ERI_D2y_D2z_Py_S;
  abcd[1906] = 2.0E0*I_ERI_Dyz_D2z_Fy2z_S_c-1*I_ERI_Dyz_D2z_Py_S;
  abcd[1907] = 2.0E0*I_ERI_D2z_D2z_Fy2z_S_c-1*I_ERI_D2z_D2z_Py_S;
  abcd[1908] = 2.0E0*I_ERI_D2x_D2x_F3z_S_c-2*I_ERI_D2x_D2x_Pz_S;
  abcd[1909] = 2.0E0*I_ERI_Dxy_D2x_F3z_S_c-2*I_ERI_Dxy_D2x_Pz_S;
  abcd[1910] = 2.0E0*I_ERI_Dxz_D2x_F3z_S_c-2*I_ERI_Dxz_D2x_Pz_S;
  abcd[1911] = 2.0E0*I_ERI_D2y_D2x_F3z_S_c-2*I_ERI_D2y_D2x_Pz_S;
  abcd[1912] = 2.0E0*I_ERI_Dyz_D2x_F3z_S_c-2*I_ERI_Dyz_D2x_Pz_S;
  abcd[1913] = 2.0E0*I_ERI_D2z_D2x_F3z_S_c-2*I_ERI_D2z_D2x_Pz_S;
  abcd[1914] = 2.0E0*I_ERI_D2x_Dxy_F3z_S_c-2*I_ERI_D2x_Dxy_Pz_S;
  abcd[1915] = 2.0E0*I_ERI_Dxy_Dxy_F3z_S_c-2*I_ERI_Dxy_Dxy_Pz_S;
  abcd[1916] = 2.0E0*I_ERI_Dxz_Dxy_F3z_S_c-2*I_ERI_Dxz_Dxy_Pz_S;
  abcd[1917] = 2.0E0*I_ERI_D2y_Dxy_F3z_S_c-2*I_ERI_D2y_Dxy_Pz_S;
  abcd[1918] = 2.0E0*I_ERI_Dyz_Dxy_F3z_S_c-2*I_ERI_Dyz_Dxy_Pz_S;
  abcd[1919] = 2.0E0*I_ERI_D2z_Dxy_F3z_S_c-2*I_ERI_D2z_Dxy_Pz_S;
  abcd[1920] = 2.0E0*I_ERI_D2x_Dxz_F3z_S_c-2*I_ERI_D2x_Dxz_Pz_S;
  abcd[1921] = 2.0E0*I_ERI_Dxy_Dxz_F3z_S_c-2*I_ERI_Dxy_Dxz_Pz_S;
  abcd[1922] = 2.0E0*I_ERI_Dxz_Dxz_F3z_S_c-2*I_ERI_Dxz_Dxz_Pz_S;
  abcd[1923] = 2.0E0*I_ERI_D2y_Dxz_F3z_S_c-2*I_ERI_D2y_Dxz_Pz_S;
  abcd[1924] = 2.0E0*I_ERI_Dyz_Dxz_F3z_S_c-2*I_ERI_Dyz_Dxz_Pz_S;
  abcd[1925] = 2.0E0*I_ERI_D2z_Dxz_F3z_S_c-2*I_ERI_D2z_Dxz_Pz_S;
  abcd[1926] = 2.0E0*I_ERI_D2x_D2y_F3z_S_c-2*I_ERI_D2x_D2y_Pz_S;
  abcd[1927] = 2.0E0*I_ERI_Dxy_D2y_F3z_S_c-2*I_ERI_Dxy_D2y_Pz_S;
  abcd[1928] = 2.0E0*I_ERI_Dxz_D2y_F3z_S_c-2*I_ERI_Dxz_D2y_Pz_S;
  abcd[1929] = 2.0E0*I_ERI_D2y_D2y_F3z_S_c-2*I_ERI_D2y_D2y_Pz_S;
  abcd[1930] = 2.0E0*I_ERI_Dyz_D2y_F3z_S_c-2*I_ERI_Dyz_D2y_Pz_S;
  abcd[1931] = 2.0E0*I_ERI_D2z_D2y_F3z_S_c-2*I_ERI_D2z_D2y_Pz_S;
  abcd[1932] = 2.0E0*I_ERI_D2x_Dyz_F3z_S_c-2*I_ERI_D2x_Dyz_Pz_S;
  abcd[1933] = 2.0E0*I_ERI_Dxy_Dyz_F3z_S_c-2*I_ERI_Dxy_Dyz_Pz_S;
  abcd[1934] = 2.0E0*I_ERI_Dxz_Dyz_F3z_S_c-2*I_ERI_Dxz_Dyz_Pz_S;
  abcd[1935] = 2.0E0*I_ERI_D2y_Dyz_F3z_S_c-2*I_ERI_D2y_Dyz_Pz_S;
  abcd[1936] = 2.0E0*I_ERI_Dyz_Dyz_F3z_S_c-2*I_ERI_Dyz_Dyz_Pz_S;
  abcd[1937] = 2.0E0*I_ERI_D2z_Dyz_F3z_S_c-2*I_ERI_D2z_Dyz_Pz_S;
  abcd[1938] = 2.0E0*I_ERI_D2x_D2z_F3z_S_c-2*I_ERI_D2x_D2z_Pz_S;
  abcd[1939] = 2.0E0*I_ERI_Dxy_D2z_F3z_S_c-2*I_ERI_Dxy_D2z_Pz_S;
  abcd[1940] = 2.0E0*I_ERI_Dxz_D2z_F3z_S_c-2*I_ERI_Dxz_D2z_Pz_S;
  abcd[1941] = 2.0E0*I_ERI_D2y_D2z_F3z_S_c-2*I_ERI_D2y_D2z_Pz_S;
  abcd[1942] = 2.0E0*I_ERI_Dyz_D2z_F3z_S_c-2*I_ERI_Dyz_D2z_Pz_S;
  abcd[1943] = 2.0E0*I_ERI_D2z_D2z_F3z_S_c-2*I_ERI_D2z_D2z_Pz_S;
}
